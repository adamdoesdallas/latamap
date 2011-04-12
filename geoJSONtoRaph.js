var allPolygons = [];

function setPolys(geom) {
    /* geom: a GeoJSON object.
       Parse geom, and add it as a child to the global Raphael object, R.
    */
    //geom.type either Polygon or MultiPolygon
    var translationFunction = function(coords) {
	    var projected = Projections.contiguous.forward(coords);
	    var x = 217+(projected[0]*190/2400000);
	    var y = 260+(projected[1]*-190/2400000);
	    return x + " " + y;
	};

    parseGeom(geom, translationFunction);
};
function parseGeom(geom, translationFunction) {
    /* geom: a GeoJSON object.
       translationFunction: a function that positions itself properly onto the Raphael canvas
       Parse geom, and add it as a child to the global Raphael object, R.
    */
    function parsePolygon(polygon) {
	//a polygon is an array of linestrings, each of which will become a subpath
	var svgstring = "";
	var nls = polygon.length;
	for (var i=0; i < nls; i++) {
	    //and each linestring is an array of [x, y] coords.
	    svgstring += "M "+translationFunction([polygon[i][0][0], polygon[i][0][1]]);
	    svgstring += "L "+$.map(polygon[i].slice(1),
				    translationFunction //join each [x, y] into "x' y'"
				   ).join(" "); //join each "x' y'" into "x1 y1 x2 y2 ..."
	    svgstring += " Z ";
	}
	var svgpath = R.path(svgstring);
	return svgpath;
    }
    if (geom.type == "Polygon") {
	var polygon = parsePolygon(geom.coordinates);
	allPolygons.push(polygon); // = this[this.mapType].polygons.concat([polygon]);
    } else if (geom.type == "MultiPolygon") {
	var polygons = $.map(geom.coordinates, function(el, i) {
	    return parsePolygon(el);
	});
	allPolygons = allPolygons.concat(polygons);
    }
};


$(document).ready(function() {
    R = Raphael(document.getElementById("map"), 640, 480);
    $.getJSON('/path/to/geojson/file.json', function(data, status) {
	setPolys(data)
    });
});


/**
 * // snipped from http://google-maps-utility-library-v3.googlecode.com/svn/trunk/arcgislink/src/arcgislink.js
 * // thanks to Nianwei Liu!
 *
 * Create a Albers Equal-Area Conic Projection based Spatial Reference. The <code>params</code> passed in construction should
 * include the following properties:<code>
 * <br/>-wkid: well-known id
 * <br/>-semi_major:  ellipsoidal semi-major axis in meter
 * <br/>-unit: meters per unit
 * <br/>-inverse_flattening: inverse of flattening of the ellipsoid where 1/f  =  a/(a - b)
 * <br/>-standard_parallel_1: phi1, latitude of the first standard parallel
 * <br/>-standard_parallel_2: phi2, latitude of the second standard parallel
 * <br/>-latitude_of_origin: phi0, latitude of the false origin
 * <br/>-central_meridian: lamda0, longitude of the false origin  (with respect to the prime meridian)
 * <br/>-false_easting: FE, false easting, the Eastings value assigned to the natural origin
 * <br/>-false_northing: FN, false northing, the Northings value assigned to the natural origin
 * </code>
 * <br/> e.g. 
 * <code> var albers  = new Albers({wkid:9999, semi_major: 6378206.4,inverse_flattening: 294.9786982,
 *   standard_parallel_1: 29.5, standard_parallel_2: 45.5,
 *   central_meridian: -96.0, latitude_of_origin: 23,false_easting: 0,
 *   'false_northing': 0, unit: 1 }); </code>
 * @name Albers
 * @class This class (<code>Albers</code>) represents a Spatial Reference System based on <a target=wiki href  = 'http://en.wikipedia.org/wiki/Albers_projection'>Albers Projection</a>. 
 * @extends SpatialReference
 * @constructor
 * @param {Object} params
 */
function Albers(params) {

    //http://pubs.er.usgs.gov/djvu/PP/PP_1395.pdf, page 101 &  292
    //for NAD_1983_Alaska_Albers: LatLng()<  === > Point();
    params = params || {};
    //SpatialReference.call(this, params);
    var f_i = params.inverse_flattening;
    var phi1 = params.standard_parallel_1 * this.RAD_DEG;
    var phi2 = params.standard_parallel_2 * this.RAD_DEG;
    var phi0 = params.latitude_of_origin * this.RAD_DEG;
    this.a_ = params.semi_major / params.unit;
    this.lamda0_ = params.central_meridian * this.RAD_DEG;
    this.FE_ = params.false_easting;
    this.FN_ = params.false_northing;
    
    var f = 1.0 / f_i; //e: eccentricity of the ellipsoid where e^2  =  2f - f^2 
    var es = 2 * f - f * f;
    this.e_ = Math.sqrt(es);
    var m1 = this.calc_m_(phi1, es);
    var m2 = this.calc_m_(phi2, es);
    
    var q1 = this.calc_q_(phi1, this.e_);
    var q2 = this.calc_q_(phi2, this.e_);
    var q0 = this.calc_q_(phi0, this.e_);
    
    this.n_ = (m1 * m1 - m2 * m2) / (q2 - q1);
    this.C_ = m1 * m1 + this.n_ * q1;
    this.rho0_ = this.calc_rho_(this.a_, this.C_, this.n_, q0);
}

Albers.prototype.RAD_DEG = Math.PI / 180;
/**
 * calc_m_
 * @param {number} phi
 * @param {number} es e square
 */
Albers.prototype.calc_m_ = function(phi, es) {
    var sinphi = Math.sin(phi);
    return Math.cos(phi) / Math.sqrt(1 - es * sinphi * sinphi);
};


/**
 * formular (3-12) page 101
 * @param {Object} phi
 * @param {Object} e
 */
Albers.prototype.calc_q_ = function(phi, e) {
    var esp = e * Math.sin(phi);
    return (1 - e * e) * (Math.sin(phi) / (1 - esp * esp) - (1 / (2 * e)) * Math.log((1 - esp) / (1 + esp)));
};

Albers.prototype.calc_rho_ = function(a, C, n, q) {
    return a * Math.sqrt(C - n * q) / n;
};

Albers.prototype.calc_phi_ = function(q, e, phi) {
    var esp = e * Math.sin(phi);
    return phi + (1 - esp * esp) * (1 - esp * esp) / (2 * Math.cos(phi)) * (q / (1 - e * e) - Math.sin(phi) / (1 - esp * esp) + Math.log((1 - esp) / (1 + esp)) / (2 * e));
};

Albers.prototype.solve_phi_ = function(q, e, init) {
    // iteration
    var i = 0;
    var phi = init;
    var newphi = this.calc_phi_(q, e, phi);
    while (Math.abs(newphi - phi) > 0.00000001 && i < 10) {
	i++;
	phi = newphi;
	newphi = this.calc_phi_(q, e, phi);
    }
    return newphi;
};

/** 
 * see {@link SpatialReference}
 * @param {Array.number} lnglat
 * @return {Array.number}
 */
Albers.prototype.forward = function(lnglat) {
    var phi = lnglat[1] * this.RAD_DEG;
    var lamda = lnglat[0] * this.RAD_DEG;
    var q = this.calc_q_(phi, this.e_);
    var rho = this.calc_rho_(this.a_, this.C_, this.n_, q);
    var theta = this.n_ * (lamda - this.lamda0_);
    var E = this.FE_ + rho * Math.sin(theta);
    var N = this.FN_ + this.rho0_ - rho * Math.cos(theta);
    return [E, N];
};
/**
 * see {@link SpatialReference}
 * @param {Array.number}  coords
 * @return {Array.number}
 */
Albers.prototype.inverse = function(coords) {
    var E = coords[0] - this.FE_;
    var N = coords[1] - this.FN_;
    var rho = Math.sqrt(E * E + (this.rho0_ - N) * (this.rho0_ - N)); 
    var adj = this.n_ > 0 ? 1 : -1;
    var theta = Math.atan(adj * E / (adj * this.rho0_  - adj * N));
    var q = (this.C_ - rho * rho * this.n_ * this.n_ / (this.a_ * this.a_)) / this.n_;
    var init = Math.asin(q / 2);
    var phi = this.solve_phi_(q, this.e_, init);
    var lamda = theta / this.n_ + this.lamda0_;
    return [lamda / this.RAD_DEG, phi / this.RAD_DEG];
};
/**
 *  see {@link SpatialReference}
 *  @return {number}
 */
Albers.prototype.getCircum = function() {
    return Math.PI * 2 * this.a_;
};

var Projections = {
    "contiguous": new Albers({
	unit: 1,
	semi_major: 6378137,
	inverse_flattening: 298.257222101,
	// in the wkt the below is called "latitude_of_center" and is 37.5 degrees
	//arcgislink.js's Util.registerSR function would fail to recognize it.
	//we want all our y values positive, so picking 23 as below all of continental US.
	latitude_of_origin: 23.0, 
	central_meridian: -96.0,
	false_easting: 0,
	false_northing: 0,
	standard_parallel_1: 29.5,
	standard_parallel_2: 45.5
    })
};
