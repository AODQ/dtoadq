////////////////////////////////////////////////////////////////
//
//                           HG_SDF
//
//     GLSL LIBRARY FOR BUILDING SIGNED DISTANCE BOUNDS
//     (I ported to opencl)
//
//     version 2016-01-10
//
//     Check http://mercury.sexy/hg_sdf for updates
//     and usage examples. Send feedback to spheretracing@mercury.sexy.
//
//     Brought to you by MERCURY http://mercury.sexy
//
//
//
// Released as Creative Commons Attribution-NonCommercial (CC BY-NC)
//
////////////////////////////////////////////////////////////////


// Rotate around axis by angle (ei: xz rotates x towards z)
float2 opRotate ( float2 p, float a  ) {
  return cos(a)*p + sin(a)*(float2)(p.y, -p.x);
}

float2 opRotate45 ( float2 p ) {
  return (p + (float2)(p.y, -p.x))*sqrt(0.5f);
}

// Repeat space along one axis
float opMod1 ( float p, float size, float* id ) {
  float hs = size*0.5f;
  if ( id ) *id = floor((p + hs)/size);
  return fmod(p + hs, size) - hs;
}
float2 opMod2 ( float2 p, float2 size, float2* id ) {
  float2 hs = size*0.5f;
  if ( id ) *id = floor((p + hs)/size);
  return fmod(p + hs, size) - hs;
}
float3 opMod3 ( float3 p, float3 size, float3* id ) {
  float3 hs = size*0.5f;
  if ( id ) *id = floor((p + hs)/size);
  return fmod(p + hs, size) - hs;
}

// Mirror every second cell to match boundaries
float opModMirror1 ( float p, float size, float* id ) {
  float hs = size*0.5f;
  float c = floor((p + hs)/size);
  if ( id ) {*id = c;}
  p = fmod(p + hs, size) - hs;
  return p*(fmod(c, 2.0f)*2.0f - 1.0f);
}
float2 opModMirror2 ( float2 p, float2 size, float2* id ) {
  float2 hs = size*0.5f;
  float2 c = floor((p + hs)/size);
  if ( id ) *id = c;
  p = fmod(p + hs, size) - hs;
  return p*(fmod(c, 2.0f)*2.0f - 1.0f);
}

// Repeat domain in positive direction
float opModSingle1 ( float p, float size, float* id ) {
  float hs = size*0.5f;
  if ( id ) *id = floor((p+hs)/size);
  if ( p >= 0.0f ) p = fmod(p+hs, size) - hs;
  return p;
}

// Repeat only a few times, from start to stop
float opModInterval1(float p, float size, float start, float stop, float* id ) {
  float hs = size*0.5;
  float c = floor((p+hs)/size);
  p = fmod(p+hs, size) - hs;
  if ( c > stop ) {
    p += size*(c - stop);
    c = stop;
  }
  if ( c < start ) {
    p += size*(c - start);
    c = start;
  }
  if ( id ) *id = c;
  return p;
}

// Repeat around origin by a fixed angle
float2 opModPolar ( float2 p, float repetitions, float* id ) {
  float angle = 2*PI/repetitions;
  float a = atan2(p.y, p.x) + angle/2.0f;
  float c = floor(a/angle);
  a = fmod(a, angle) - angle/2.0f;
  p = (float2)(cos(a), sin(a))*length(p);
  // For odd number of repetitions, fix cell index of the cell in -x direction
  if ( fabs(c) >= (repetitions/2.0f) ) c = fabs(c);
  if ( id ) *id = c;
  return p;
}

// similar to modmirror2 except mirrors at diagonal
float2 opModGrid2 ( float2 p, float2 size, float2* id ) {
  float2 hs = size*0.5f;
  float2 c = floor((p + hs)/hs);
  p  = fmod(p + hs, size) - hs;
  p *= fmod(c, (float2)(2.0f))*2.0f - (float2)(1.0f);
  p -= hs;
  if ( p.x > p.y ) p.xy = p.yx;
  if ( id ) *id = floor(c/2.0f);
  return p;
}

// Mirror at an axis-aligned plane at a specific distance from origin
float opMirror ( float p, float dist, float* id ) {
  if ( id ) *id = sign(p);
  return fabs(p) - dist;
}

// Mirror in both dimensions and at diagonal, yielding 1/8th of the space
float2 opMirrorOctant ( float2 p, float2 dist, float2* id ) {
  if ( id ) *id = sign(p);
  p.x = opMirror(p.x, dist.x, 0);
  p.y = opMirror(p.y, dist.y, 0);
  if ( p.y > p.x ) p.xy = p.yx;
  return p;
}


// Reflect space at a plane
float3 opReflect ( float3 p, float3 plane_normal, float offset, float* id ) {
  float t = dot(p, plane_normal) + offset;
  if ( t < 0.0f ) p = p - (2.0f*t)*plane_normal;
  if ( id ) *id = sign(t);
  return p;
}

// The "Chamfer" flavour makes a 45-degree chamfered edge
// (the diagonal of a square of size <r>):
float opUnionChamfer(float a, float b, float r) {
	return fmin(fmin(a, b), (a - r + b)*sqrt(0.5f));
}

// Intersection has to deal with what is normally the inside of the resulting object
// when using union, which we normally don't care about too much. Thus, intersection
// implementations sometimes differ from union implementations.
float opIntersectionChamfer(float a, float b, float r) {
	return fmax(fmax(a, b), (a + r + b)*sqrt(0.5f));
}

// Difference can be built from Intersection or Union:
float opSubtractChamfer (float a, float b, float r) {
	return opIntersectionChamfer(a, -b, r);
}

// The "Round" variant uses a quarter-circle to join the two objects smoothly:
float opUnionRound(float a, float b, float r) {
	float2 u = fmax((float2)(r - a,r - b), (float2)(0.0f));
	return fmax(r, fmin (a, b)) - length(u);
}

float opIntersectionRound(float a, float b, float r) {
	float2 u = fmax((float2)(r + a,r + b), (float2)(0.0f));
	return fmin(-r, fmax (a, b)) + length(u);
}

float opSubtractRound (float a, float b, float r) {
	return opIntersectionRound(a, -b, r);
}


// The "Columns" flavour makes n-1 circular columns at a 45 degree angle:
float opUnionColumns(float a, float b, float r, float n) {
	if ((a < r) && (b < r)) {
		float2 p = (float2)(a, b);
		float columnradius = r*sqrt(2.0f)/((n-1.0f)*2.0f+sqrt(2.0f));
		p = opRotate45(p);
		p.x -= sqrt(2.0f)/2.0f*r;
		p.x += columnradius*sqrt(2.0f);
		if ( fmod(n,2.0f) == 1.0f ) {
			p.y += columnradius;
		}
		// At this point, we have turned 45 degrees and moved at a point on the
		// diagonal that we want to place the columns on.
		// Now, repeat the domain along this direction and place a circle.
		p.y = opMod1(p.y, columnradius*2.0f, 0);
		float result = length(p) - columnradius;
		result = fmin(result, p.x);
		result = fmin(result, a);
		return fmin(result, b);
	} else {
		return fmin(a, b);
	}
}

float opSubtractColumns(float a, float b, float r, float n) {
	a = -a;
	float m = fmin(a, b);
	//avoid the expensive computation where not needed (produces discontinuity though)
	if ((a < r) && (b < r)) {
		float2 p = (float2)(a, b);
		float columnradius = r*sqrt(2.0f)/n/2.0f;
		columnradius = r*sqrt(2.0f)/((n-1.0f)*2.0f+sqrt(2.0f));

		opRotate45(p);
		p.y += columnradius;
		p.x -= sqrt(2.0f)/2.0f*r;
		p.x += -columnradius*sqrt(2.0f)/2.0f;

		if ( fmod(n, 2.0f) == 1.0f ) {
			p.y += columnradius;
		}
		p.y = opMod1(p.y,columnradius*2.0f, 0);

		float result = -length(p) + columnradius;
		result = fmax(result, p.x);
		result = fmin(result, a);
		return -fmin(result, b);
	} else {
		return -m;
	}
}

float opIntersectionColumns(float a, float b, float r, float n) {
	return opSubtractColumns(a,-b,r, n);
}

// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float opUnionStairs(float a, float b, float r, float n) {
	float s = r/n;
	float u = b-r;
	return fmin(fmin(a,b), 0.5f * (u + a + fabs ((fmod (u - a + s, 2.0f * s)) - s)));
}

// We can just call Union since stairs are symmetric.
float opIntersectionStairs(float a, float b, float r, float n) {
	return -opUnionStairs(-a, -b, r, n);
}

float opSubtractStairs(float a, float b, float r, float n) {
	return -opUnionStairs(-a, b, r, n);
}


// Similar to opUnionRound, but more lipschitz-y at acute angles
// (and less so at 90 degrees). Useful when fudging around too much
// by MediaMolecule, from Alex Evans' siggraph slides
float opUnionSoft(float a, float b, float r) {
	float e = fmax(r - fabs(a - b), 0.0f);
	return fmin(a, b) - e*e*0.25f/r;
}


// produces a cylindical pipe that runs along the intersection.
// No objects remain, only the pipe. This is not a boolean operator.
float opPipe(float a, float b, float r) {
	return length((float2)(a, b)) - r;
}

// first object gets a v-shaped engraving where it intersect the second
float opEngrave(float a, float b, float r) {
	return fmax(a, (a + r - fabs(b))*sqrt(0.5f));
}

// first object gets a capenter-style groove cut out
float opGroove(float a, float b, float ra, float rb) {
	return fmax(a, fmin(a + ra, rb - fabs(b)));
}

// first object gets a capenter-style tongue attached
float opTongue(float a, float b, float ra, float rb) {
	return fmin(a, fmax(a - ra, fabs(b) - rb));
}