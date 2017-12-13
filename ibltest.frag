#define iGlobalTime iTime
#define PI 3.14159265359

// control with time or mouse
#define REVOLVE false
// Buf B gallery needs to be set to this value
#define GALLERY true
// samples to collect per frame (also higher overall quality image)
#define SAMPLES 5
#define SPEED 0.4

float sqr ( float s ) { return s * s; }
vec3 sqr ( vec3 s ) { return s * s; }

// change above to change the BSDF/material
// The 'template' was generated from the elegant mathematical modelling tool
// CURV https://github.com/doug-moen/curv/
// and I added bsdf IBL, etc to it

float hash(float f) { // Forgot where i got it from but any hash will do
  vec2 p = vec2(f, fract(sin(f*23423.234)));
  #define HASHSCALE1 .1031
  vec3 p3  = fract(vec3(p.xyx) * HASHSCALE1);
  p3 += dot(p3, p3.yzx + 19.19);
  return fract((p3.x + p3.y) * p3.z);
}
float hash( vec2 co ) {
  return fract(sin(dot(co.xy, vec2(132.9898, 78.233))) * 43758.5453);
}
float hash(float f, float t) {
   return hash(vec2(f, t));
}

// ------------------------ random material properties ------------------------

float Rm_diffuse  ( float t ) { return hash(t); }
float Rm_specular ( float t ) {
  if ( hash(hash(hash(t))) < 0.4 )
    return (1.0 - Rm_diffuse(t))*hash(hash(t));
  return 0.0;
}
float Rm_glossy ( float t ) {
  return 1.0 - (Rm_diffuse(t) + Rm_specular(t));
}

float Rm_glossy_lobe ( float t ) { return 0.9 + 0.1*hash(t, 71341.2342); }
float Rm_fresnel     ( float t ) { return 0.7 + hash(t, 96123.46) *1.3; }
float Rm_roughness   ( float t ) { return hash(t, 523494.23423); }
float Rm_metallic    ( float t ) { return hash(t, 592623.2342); }
float Rm_subsurface  ( float t ) { return hash(t, 2340234.234); }
float Rm_anisotropic ( float t ) { return hash(t, 962823.2342934); }

// --------------------------- map functions ----------------------------------

// from mercury hg
void opRotate ( inout vec2 p, float a ) {
  p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}
vec3 opRep( vec3 p, vec3 c ) { return mod(p,c)-0.5*c; }
float opMod1 ( float p, float size, inout float id ) {
  float hs = size*0.5;
  id = floor((p + hs)/size);
  return mod(p + hs, size) - hs;
}
// hg_sdf
float OpUnionRound(float a, float b, float r) {
	vec2 u = max(vec2(r - a,r - b), vec2(0));
	return max(r, min (a, b)) - length(u);
}

float sdMetaballs ( vec3 p, int amt, float t, float ID ) {
    float dist = 993423499.0;
    for ( int i = 0; i != amt; ++ i ) {
        float fi = hash(sqrt(hash(float(i+int(ID))))*28523.2348);
        vec3 O = p.xyz;
        O.x += sin(t*hash(fi*0.85)+fi*23.0)*hash(sqr(fi)*0.2)*0.3;
        O.y += sin(t*hash(fi*0.25)+fi*53.0)*hash(sqr(fi)*0.2)*0.3;
        O.z += sin(t*hash(fi*0.55)+fi*-73.0)*hash(sqr(fi)*0.2)*0.3;
        float tdist = length(O) - 0.1;
        dist = OpUnionRound(dist, tdist, 0.07);
    }
    return dist*0.5; // lipschittz
}
float sdBox( vec3 p, vec3 b ){
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) +
         length(max(d,0.0));
}
float sbox(in vec3 p, in vec3 b)
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

vec2 foldPent(in vec2 p)
{
    p.x = abs(p.x);
    const vec2 pl1 = vec2(0.809, 0.5878);
    const vec2 pl2 = vec2(-0.309, 0.951);
   	const vec2 pl3 = vec2(-0.809, 0.5878);
    p -= pl1*2.*min(0., dot(p, pl1));
    p -= pl2*2.*min(0., dot(p, pl2));
    p -= pl3*2.*min(0., dot(p, pl3));
    return p;
}

float cyl(in vec3 p, in vec2 h)
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float smin(in float a, in float b)
{
    float k = .15;
	float h = clamp(0.5 + 0.5*(b-a)/k, 0., 1.);
	return mix( b, a, h ) - k*h*(1.0-h);
}

//---------------------------------------------------------------
// MODEL FROM::::
//https://www.shadertoy.com/view/Xtj3Dm

float ORANGEMAP(in vec3 p)
{
    float r =length(p);
    vec2 sph = vec2(acos(p.y/r), atan(p.x, p.z));
    
    float d = r-1.; 
    d += sin(sph.y*7.)*0.02;
    d += sin(sph.y*20.)*0.002;
    float gbh = sin((sph.x+sph.y)*7.+0.5)*0.5+0.5;
    d += sin(sph.y*40.)*0.001*gbh;
    d += sin(sph.x*1.85+2.7)*0.3;
    
    //Leaves
    vec3 p2 = p;
    float rxz2 = dot(p.xz,p.xz);
    float rxz = sqrt(rxz2);
    rxz = exp2(rxz*6.-5.);
    p2.xz = foldPent(p2.xz);
    p2.y -= sqrt(rxz)*0.17 + sin(rxz*2.+p.z*p.x*10.)*0.05;
    float leaves = sbox(p2+vec3(0,-.92-smoothstep(-0.01,.05,rxz2)*0.05,0),vec3(.07- rxz*0.1,0.002+p2.x*0.15,0.8));
    leaves = smin(leaves, cyl(p+vec3(sin(p.y*3.5 + 0.8)*0.3 + 0.3,-1.1,0),vec2(.05,.25))); //Tail

    d = min(d, leaves);
    
    float flor = p.y+.65;
    //d = min(d, flor);
    return d;
}

vec4 model(vec4 r0, float t, float ID) {
  float r1 = r0.x; vec3 r2 = vec3(0.0,0.0,-0.15); float r3 = r2.x; float r4 = r1-r3; float r5 = r0.y; vec3 r6 = vec3(0.0,0.0,-0.15); float r7 = r6.y; float r8 = r5-r7; float r9 = r0.z; vec3 r10 = vec3(0.0,0.0,-0.15); float r11 = r10.z; float r12 = r9-r11; float r13 = r0.w; vec4 r14 = vec4(r4,r8,r12,r13); float r15 = r14.x; vec3 r16 = vec3(0.0,0.0,0.02); float r17 = r16.x; float r18 = r15-r17; float r19 = r14.y; vec3 r20 = vec3(0.0,0.0,0.02); float r21 = r20.y; float r22 = r19-r21; float r23 = r14.z; vec3 r24 = vec3(0.0,0.0,0.02); float r25 = r24.z; float r26 = r23-r25; float r27 = r14.w; float r28 = abs(r26); float r29 = 0.1; float r30 = r28-r29; float r31 = 0.0; vec3 r32 = vec3(r18,r22,r31); vec3 r33 = abs(r32); float r34 = 0.2; vec3 r35 = r33-vec3(r34,r34,r34); float r36 = max(max(r35.x,r35.y),r35.z); float r37 = 0.0; vec3 r38 = vec3(r18,r22,r37); float r39 = length(r38); float r40 = 0.05; float r41 = r39-r40; float r42 = 0.1; float r43 = -r42; float r44 = 0.1; float r45 = max(r26,r43); float r46 = min(r45,r44); float r47 = 0.1; float r48 = r46+r47; float r49 = 0.2; float r50 = r48/r49; float r51 = 1.0; float r52 = r51-r50; float r53 = r36*r52; float r54 = r41*r50; float r55 = r53+r54; float r56 = max(r30,r55); float r57 = r14.x; vec3 r58 = vec3(0.0,0.0,-0.08); float r59 = r58.x; float r60 = r57-r59; float r61 = r14.y; vec3 r62 = vec3(0.0,0.0,-0.08); float r63 = r62.y; float r64 = r61-r63; float r65 = r14.z; vec3 r66 = vec3(0.0,0.0,-0.08); float r67 = r66.z; float r68 = r65-r67; float r69 = r14.w; float r70 = 0.0; vec2 r71 = vec2(r60,r64); float r72 = length(r71); float r73 = 0.225; float r74 = r72-r73; float r75 = 0.0; vec2 r76 = vec2(r74,r68); float r77 = length(r76); float r78 = 0.01; float r79 = r77-r78; float r80 = -r79; float r81 = 0.02; float r82 = -r56; float r83 = -r80; float r84 = 0.0; bool r85 =(r82 == r84); float r86 = 0.5; float r87 = 0.5; float r88 = r83-r82; float r89 = r87*r88; float r90 = r89/r81; float r91 = r86+r90; float r92 = 0.0; float r93 = 1.0; float r94 = max(r91,r92); float r95 = min(r94,r93); float r96=r95; float r97 = 1.0; float r98 = r97-r96; float r99 = r83*r98; float r100 = r82*r96; float r101 = r99+r100; float r102 = r81*r96; float r103 = 1.0; float r104 = r103-r96; float r105 = r102*r104; float r106 = r101-r105; float r107 =(r85 ? r83 : r106); float r108 = -r107; float r109 = r0.x; vec3 r110 = vec3(0.0,0.0,0.2); float r111 = r110.x; float r112 = r109-r111; float r113 = r0.y; vec3 r114 = vec3(0.0,0.0,0.2); float r115 = r114.y; float r116 = r113-r115; float r117 = r0.z; vec3 r118 = vec3(0.0,0.0,0.2); float r119 = r118.z; float r120 = r117-r119; float r121 = r0.w; float r122 = 3.141592653589793; vec3 r123 = vec3(0.0,1.0,0.0); vec3 r124 = vec3(r112,r116,r120); float r125 = cos(r122); vec3 r126 = r124*vec3(r125,r125,r125); float r127 = r123.y; float r128 = r124.z; float r129 = r127*r128; float r130 = r123.z; float r131 = r124.y; float r132 = r130*r131; float r133 = r129-r132; float r134 = r123.z; float r135 = r124.x; float r136 = r134*r135; float r137 = r123.x; float r138 = r124.z; float r139 = r137*r138; float r140 = r136-r139; float r141 = r123.x; float r142 = r124.y; float r143 = r141*r142; float r144 = r123.y; float r145 = r124.x; float r146 = r144*r145; float r147 = r143-r146; vec3 r148 = vec3(r133,r140,r147); float r149 = sin(r122); vec3 r150 = r148*vec3(r149,r149,r149); vec3 r151 = r126-r150; float r152 = 1.0; float r153 = cos(r122); float r154 = r152-r153; vec3 r155 = r124*vec3(r154,r154,r154); float r156 = dot(r123,r155); vec3 r157 = r123*vec3(r156,r156,r156); vec3 r158 = r151+r157; vec3 r159=r158; float r160 = r159.x; float r161 = r159.y; float r162 = r159.z; vec2 r163 = vec2(r160,r161); float r164 = length(r163); vec2 r165 = vec2(r164,r162); vec2 r166=r165; float r167 = 0.0; float r168 = 0.3; vec2 r169 = vec2(r167,r168); vec2 r170 = r166-r169; vec2 r171=r170; float r172 = 0.3; float r173 = 0.35; vec2 r174 = vec2(r172,r173); float r175 = length(r174); vec2 r176 = r174/vec2(r175,r175); vec2 r177=r176; float r178 = dot(r171,r177); float r179=r178; float r180 = r166.y; float r181 = -r180; float r182 = max(r179,r181); float r183=r182; float r184 = r177.y; float r185 = r177.x; float r186 = -r185; vec2 r187 = vec2(r184,r186); float r188 = dot(r171,r187); float r189=r188; float r190 = r166.y; float r191 = 0.3; bool r192 =(r190 > r191); float r193 = 0.0; bool r194 =(r189 < r193); bool r195 =(r192 && r194); if (r195) {float r196 = length(r171); float r197 = max(r183,r196); r183=r197;} float r198 = r166.x; float r199 = 0.35; bool r200 =(r198 > r199); float r201 = 0.3; float r202 = 0.35; vec2 r203 = vec2(r201,r202); float r204 = length(r203); bool r205 =(r189 > r204); bool r206 =(r200 && r205); if (r206) {float r207 = 0.35; float r208 = 0.0; vec2 r209 = vec2(r207,r208); vec2 r210 = r166-r209; float r211 = length(r210); float r212 = max(r183,r211); r183=r212;} float r213 = abs(r183); float r214 = 0.0; float r215 = r213-r214; float r216 = r0[0]; float r217 = r0[1]; float r218 = r0[2]; float r219 = r0[3]; vec3 r220 = vec3(r216,r217,r218); vec3 r221 = vec3(0.0,0.0,-1.0); float r222 = dot(r220,r221); float r223 = -0.1; float r224 = r222-r223; float r225 = -r224; float r226 = 0.01; float r227 = -r215; float r228 = -r225; float r229 = 0.0; bool r230 =(r227 == r229); float r231 = 0.5; float r232 = 0.5; float r233 = r228-r227; float r234 = r232*r233; float r235 = r234/r226; float r236 = r231+r235; float r237 = 0.0; float r238 = 1.0; float r239 = max(r236,r237); float r240 = min(r239,r238); float r241=r240; float r242 = 1.0; float r243 = r242-r241; float r244 = r228*r243; float r245 = r227*r241; float r246 = r244+r245; float r247 = r226*r241; float r248 = 1.0; float r249 = r248-r241; float r250 = r247*r249; float r251 = r246-r250; float r252 =(r230 ? r228 : r251); float r253 = -r252; float r254 = 0.01; float r255 = 0.0; bool r256 =(r108 == r255); float r257 = 0.5; float r258 = 0.5; float r259 = r253-r108; float r260 = r258*r259; float r261 = r260/r254; float r262 = r257+r261; float r263 = 0.0; float r264 = 1.0; float r265 = max(r262,r263); float r266 = min(r265,r264); float r267=r266; float r268 = 1.0; float r269 = r268-r267; float r270 = r253*r269; float r271 = r108*r267; float r272 = r270+r271; float r273 = r254*r267; float r274 = 1.0; float r275 = r274-r267; float r276 = r273*r275; float r277 = r272-r276; float r278 =(r256 ? r253 : r277); float r279 = r0.x; vec3 r280 = vec3(0.0,0.0,-0.15); float r281 = r280.x; float r282 = r279-r281; float r283 = r0.y; vec3 r284 = vec3(0.0,0.0,-0.15); float r285 = r284.y; float r286 = r283-r285; float r287 = r0.z; vec3 r288 = vec3(0.0,0.0,-0.15); float r289 = r288.z; float r290 = r287-r289; float r291 = r0.w; float r292 = 0.0; vec2 r293 = vec2(r282,r286); float r294 = length(r293); float r295 = 0.15; float r296 = r294-r295; float r297 = 0.0; vec2 r298 = vec2(r296,r290); float r299 = length(r298); float r300 = 0.075; float r301 = r299-r300; float r302 = r0.x; vec3 r303 = vec3(0.0,0.0,-0.04); float r304 = r303.x; float r305 = r302-r304; float r306 = r0.y; vec3 r307 = vec3(0.0,0.0,-0.04); float r308 = r307.y; float r309 = r306-r308; float r310 = r0.z; vec3 r311 = vec3(0.0,0.0,-0.04); float r312 = r311.z; float r313 = r310-r312; float r314 = r0.w; float r315 = 0.0; vec2 r316 = vec2(r305,r309); float r317 = length(r316); float r318 = 0.15; float r319 = r317-r318; float r320 = 0.0; vec2 r321 = vec2(r319,r313); float r322 = length(r321); float r323 = 0.075; float r324 = r322-r323; float r325 = 0.09; float r326 = 0.0; bool r327 =(r301 == r326); float r328 = 0.5; float r329 = 0.5; float r330 = r324-r301; float r331 = r329*r330; float r332 = r331/r325; float r333 = r328+r332; float r334 = 0.0; float r335 = 1.0; float r336 = max(r333,r334); float r337 = min(r336,r335); float r338=r337; float r339 = 1.0; float r340 = r339-r338; float r341 = r324*r340; float r342 = r301*r338; float r343 = r341+r342; float r344 = r325*r338; float r345 = 1.0; float r346 = r345-r338; float r347 = r344*r346; float r348 = r343-r347; float r349 =(r327 ? r324 : r348); float r350 = r0.x; vec3 r351 = vec3(0.0,0.2,-0.08); float r352 = r351.x; float r353 = r350-r352; float r354 = r0.y; vec3 r355 = vec3(0.0,0.2,-0.08); float r356 = r355.y; float r357 = r354-r356; float r358 = r0.z; vec3 r359 = vec3(0.0,0.2,-0.08); float r360 = r359.z; float r361 = r358-r360; float r362 = r0.w; float r363 = 0.7853981633974483; vec3 r364 = vec3(0.0,0.0,1.0); vec3 r365 = vec3(r353,r357,r361); float r366 = cos(r363); vec3 r367 = r365*vec3(r366,r366,r366); float r368 = r364.y; float r369 = r365.z; float r370 = r368*r369; float r371 = r364.z; float r372 = r365.y; float r373 = r371*r372; float r374 = r370-r373; float r375 = r364.z; float r376 = r365.x; float r377 = r375*r376; float r378 = r364.x; float r379 = r365.z; float r380 = r378*r379; float r381 = r377-r380; float r382 = r364.x; float r383 = r365.y; float r384 = r382*r383; float r385 = r364.y; float r386 = r365.x; float r387 = r385*r386; float r388 = r384-r387; vec3 r389 = vec3(r374,r381,r388); float r390 = sin(r363); vec3 r391 = r389*vec3(r390,r390,r390); vec3 r392 = r367-r391; float r393 = 1.0; float r394 = cos(r363); float r395 = r393-r394; vec3 r396 = r365*vec3(r395,r395,r395); float r397 = dot(r364,r396); vec3 r398 = r364*vec3(r397,r397,r397); vec3 r399 = r392+r398; vec3 r400=r399; float r401 = r400.x; float r402 = r400.y; float r403 = r400.z; vec4 r404 = vec4(r401,r402,r403,r362); vec3 r405 = r404.xyz; vec3 r406 = vec3(0.5773502691896258,-0.5773502691896258,0.5773502691896258); float r407 = dot(r405,r406); vec3 r408 = r404.xyz; vec3 r409 = vec3(-0.5773502691896258,0.5773502691896258,0.5773502691896258); float r410 = dot(r408,r409); vec3 r411 = r404.xyz; vec3 r412 = vec3(-0.5773502691896258,-0.5773502691896258,-0.5773502691896258); float r413 = dot(r411,r412); vec3 r414 = r404.xyz; vec3 r415 = vec3(0.5773502691896258,0.5773502691896258,-0.5773502691896258); float r416 = dot(r414,r415); float r417 = max(r407,max(r410,max(r413,r416))); float r418 = 0.11; float r419 = r417-r418; float r420 = -r419; float r421 = 0.01; float r422 = -r349; float r423 = -r420; float r424 = 0.0; bool r425 =(r422 == r424); float r426 = 0.5; float r427 = 0.5; float r428 = r423-r422; float r429 = r427*r428; float r430 = r429/r421; float r431 = r426+r430; float r432 = 0.0; float r433 = 1.0; float r434 = max(r431,r432); float r435 = min(r434,r433); float r436=r435; float r437 = 1.0; float r438 = r437-r436; float r439 = r423*r438; float r440 = r422*r436; float r441 = r439+r440; float r442 = r421*r436; float r443 = 1.0; float r444 = r443-r436; float r445 = r442*r444; float r446 = r441-r445; float r447 =(r425 ? r423 : r446); float r448 = -r447; float r449 = 0.01; float r450 = 0.0; bool r451 =(r278 == r450); float r452 = 0.5; float r453 = 0.5; float r454 = r448-r278; float r455 = r453*r454; float r456 = r455/r449; float r457 = r452+r456; float r458 = 0.0; float r459 = 1.0; float r460 = max(r457,r458); float r461 = min(r460,r459); float r462=r461; float r463 = 1.0; float r464 = r463-r462; float r465 = r448*r464; float r466 = r278*r462; float r467 = r465+r466; float r468 = r449*r462; float r469 = 1.0; float r470 = r469-r462; float r471 = r468*r470; float r472 = r467-r471; float r473 =(r451 ? r448 : r472); float r474 = r0.x; vec3 r475 = vec3(0.0,0.144,-0.15); float r476 = r475.x; float r477 = r474-r476; float r478 = r0.y; vec3 r479 = vec3(0.0,0.144,-0.15); float r480 = r479.y; float r481 = r478-r480; float r482 = r0.z; vec3 r483 = vec3(0.0,0.144,-0.15); float r484 = r483.z; float r485 = r482-r484; float r486 = r0.w; float r487 = 0.7853981633974483; vec3 r488 = vec3(1.0,0.0,0.0); vec3 r489 = vec3(r477,r481,r485); float r490 = cos(r487); vec3 r491 = r489*vec3(r490,r490,r490); float r492 = r488.y; float r493 = r489.z; float r494 = r492*r493; float r495 = r488.z; float r496 = r489.y; float r497 = r495*r496; float r498 = r494-r497; float r499 = r488.z; float r500 = r489.x; float r501 = r499*r500; float r502 = r488.x; float r503 = r489.z; float r504 = r502*r503; float r505 = r501-r504; float r506 = r488.x; float r507 = r489.y; float r508 = r506*r507; float r509 = r488.y; float r510 = r489.x; float r511 = r509*r510; float r512 = r508-r511; vec3 r513 = vec3(r498,r505,r512); float r514 = sin(r487); vec3 r515 = r513*vec3(r514,r514,r514); vec3 r516 = r491-r515; float r517 = 1.0; float r518 = cos(r487); float r519 = r517-r518; vec3 r520 = r489*vec3(r519,r519,r519); float r521 = dot(r488,r520); vec3 r522 = r488*vec3(r521,r521,r521); vec3 r523 = r516+r522; vec3 r524=r523; float r525 = r524.x; float r526 = r524.y; float r527 = r524.z; vec3 r528 = vec3(0.1,0.1,0.1); float r529 = r528.x; float r530 = r525/r529; vec3 r531 = vec3(0.1,0.1,0.1); float r532 = r531.y; float r533 = r526/r532; vec3 r534 = vec3(0.1,0.1,0.1); float r535 = r534.z; float r536 = r527/r535; float r537 = -1.5707963267948966; vec3 r538 = vec3(1.0,0.0,0.0); vec3 r539 = vec3(r530,r533,r536); float r540 = cos(r537); vec3 r541 = r539*vec3(r540,r540,r540); float r542 = r538.y; float r543 = r539.z; float r544 = r542*r543; float r545 = r538.z; float r546 = r539.y; float r547 = r545*r546; float r548 = r544-r547; float r549 = r538.z; float r550 = r539.x; float r551 = r549*r550; float r552 = r538.x; float r553 = r539.z; float r554 = r552*r553; float r555 = r551-r554; float r556 = r538.x; float r557 = r539.y; float r558 = r556*r557; float r559 = r538.y; float r560 = r539.x; float r561 = r559*r560; float r562 = r558-r561; vec3 r563 = vec3(r548,r555,r562); float r564 = sin(r537); vec3 r565 = r563*vec3(r564,r564,r564); vec3 r566 = r541-r565; float r567 = 1.0; float r568 = cos(r537); float r569 = r567-r568; vec3 r570 = r539*vec3(r569,r569,r569); float r571 = dot(r538,r570); vec3 r572 = r538*vec3(r571,r571,r571); vec3 r573 = r566+r572; vec3 r574=r573; float r575 = r574.x; float r576 = r574.y; float r577 = r574.z; vec4 r578 = vec4(r575,r576,r577,r486); float r579 = r578.x; vec3 r580 = vec3(-1.3,-1.8,0.0); float r581 = r580.x; float r582 = r579-r581; float r583 = r578.y; vec3 r584 = vec3(-1.3,-1.8,0.0); float r585 = r584.y; float r586 = r583-r585; float r587 = r578.z; vec3 r588 = vec3(-1.3,-1.8,0.0); float r589 = r588.z; float r590 = r587-r589; float r591 = r578.w; vec4 r592 = vec4(r582,r586,r590,r591); float r593 = r592.z; float r594 = abs(r593); float r595 = 0.1; float r596 = r594-r595; float r597=r596; float r598 = r592.x; float r599 = r592.y; float r600 = 0.0; float r601 = r592.w; vec4 r602 = vec4(r598,r599,r600,r601); vec2 r603 = r602.xy; vec2 r604 = vec2(1.0,1.0); vec2 r605 = r603-r604; vec2 r606=r605; vec2 r607 = vec2(-0.30000000000000004,0.10000000000000009); float r608 = dot(r606,r607); vec2 r609 = vec2(-0.30000000000000004,0.10000000000000009); vec2 r610 = vec2(-0.30000000000000004,0.10000000000000009); float r611 = dot(r609,r610); float r612 = r608/r611; float r613 = 0.0; float r614 = 1.0; float r615 = max(r612,r613); float r616 = min(r615,r614); float r617=r616; vec2 r618 = vec2(-0.30000000000000004,0.10000000000000009); vec2 r619 = r618*vec2(r617,r617); vec2 r620 = r606-r619; float r621 = length(r620); float r622 = 0.05; float r623 = r621-r622; vec2 r624 = r602.xy; vec2 r625 = vec2(2.0,2.0); vec2 r626 = r624-r625; vec2 r627=r626; vec2 r628 = vec2(-0.30000000000000004,0.5); float r629 = dot(r627,r628); vec2 r630 = vec2(-0.30000000000000004,0.5); vec2 r631 = vec2(-0.30000000000000004,0.5); float r632 = dot(r630,r631); float r633 = r629/r632; float r634 = 0.0; float r635 = 1.0; float r636 = max(r633,r634); float r637 = min(r636,r635); float r638=r637; vec2 r639 = vec2(-0.30000000000000004,0.5); vec2 r640 = r639*vec2(r638,r638); vec2 r641 = r627-r640; float r642 = length(r641); float r643 = 0.05; float r644 = r642-r643; float r645 = 0.73; float r646 = 0.0; bool r647 =(r623 == r646); float r648 = 0.5; float r649 = 0.5; float r650 = r644-r623; float r651 = r649*r650; float r652 = r651/r645; float r653 = r648+r652; float r654 = 0.0; float r655 = 1.0; float r656 = max(r653,r654); float r657 = min(r656,r655); float r658=r657; float r659 = 1.0; float r660 = r659-r658; float r661 = r644*r660; float r662 = r623*r658; float r663 = r661+r662; float r664 = r645*r658; float r665 = 1.0; float r666 = r665-r658; float r667 = r664*r666; float r668 = r663-r667; float r669 =(r647 ? r644 : r668); vec2 r670 = r602.xy; vec2 r671 = vec2(1.2,1.0); vec2 r672 = r670-r671; vec2 r673=r672; vec2 r674 = vec2(0.19999999999999996,1.7000000000000002); float r675 = dot(r673,r674); vec2 r676 = vec2(0.19999999999999996,1.7000000000000002); vec2 r677 = vec2(0.19999999999999996,1.7000000000000002); float r678 = dot(r676,r677); float r679 = r675/r678; float r680 = 0.0; float r681 = 1.0; float r682 = max(r679,r680); float r683 = min(r682,r681); float r684=r683; vec2 r685 = vec2(0.19999999999999996,1.7000000000000002); vec2 r686 = r685*vec2(r684,r684); vec2 r687 = r673-r686; float r688 = length(r687); float r689 = 0.05; float r690 = r688-r689; float r691 = 0.73; float r692 = 0.0; bool r693 =(r669 == r692); float r694 = 0.5; float r695 = 0.5; float r696 = r690-r669; float r697 = r695*r696; float r698 = r697/r691; float r699 = r694+r698; float r700 = 0.0; float r701 = 1.0; float r702 = max(r699,r700); float r703 = min(r702,r701); float r704=r703; float r705 = 1.0; float r706 = r705-r704; float r707 = r690*r706; float r708 = r669*r704; float r709 = r707+r708; float r710 = r691*r704; float r711 = 1.0; float r712 = r711-r704; float r713 = r710*r712; float r714 = r709-r713; float r715 =(r693 ? r690 : r714); float r716=r715; vec2 r717 = vec2(r597,r716); float r718 = 0.0; vec2 r719 = max(r717,r718); float r720 = length(r719); float r721 = max(r597,r716); float r722 = 0.0; float r723 = min(r721,r722); float r724 = r720+r723; vec3 r725 = vec3(0.1,0.1,0.1); float r726 = min(min(r725.x,r725.y),r725.z); float r727 = r724*r726; float r728 = 0.02; float r729 = 0.0; bool r730 =(r473 == r729); float r731 = 0.5; float r732 = 0.5; float r733 = r727-r473; float r734 = r732*r733; float r735 = r734/r728; float r736 = r731+r735; float r737 = 0.0; float r738 = 1.0; float r739 = max(r736,r737); float r740 = min(r739,r738); float r741=r740; float r742 = 1.0; float r743 = r742-r741; float r744 = r727*r743; float r745 = r473*r741; float r746 = r744+r745; float r747 = r728*r741; float r748 = 1.0; float r749 = r748-r741; float r750 = r747*r749; float r751 = r746-r750; float r752 =(r730 ? r727 : r751); float r753 = abs(r752); float r754 = 0.0; float r755 = r753-r754; float r756 = 4.0; float r757 = r755/r756; float r758 = r0.x; vec3 r759 = vec3(0.0,0.15,0.04); float r760 = r759.x; float r761 = r758-r760; float r762 = r0.y; vec3 r763 = vec3(0.0,0.15,0.04); float r764 = r763.y; float r765 = r762-r764; float r766 = r0.z; vec3 r767 = vec3(0.0,0.15,0.04); float r768 = r767.z; float r769 = r766-r768; float r770 = r0.w; float r771 = 0.39269908169872414; vec3 r772 = vec3(1.0,0.0,0.0); vec3 r773 = vec3(r761,r765,r769); float r774 = cos(r771); vec3 r775 = r773*vec3(r774,r774,r774); float r776 = r772.y; float r777 = r773.z; float r778 = r776*r777; float r779 = r772.z; float r780 = r773.y; float r781 = r779*r780; float r782 = r778-r781; float r783 = r772.z; float r784 = r773.x; float r785 = r783*r784; float r786 = r772.x; float r787 = r773.z; float r788 = r786*r787; float r789 = r785-r788; float r790 = r772.x; float r791 = r773.y; float r792 = r790*r791; float r793 = r772.y; float r794 = r773.x; float r795 = r793*r794; float r796 = r792-r795; vec3 r797 = vec3(r782,r789,r796); float r798 = sin(r771); vec3 r799 = r797*vec3(r798,r798,r798); vec3 r800 = r775-r799; float r801 = 1.0; float r802 = cos(r771); float r803 = r801-r802; vec3 r804 = r773*vec3(r803,r803,r803); float r805 = dot(r772,r804); vec3 r806 = r772*vec3(r805,r805,r805); vec3 r807 = r800+r806; vec3 r808=r807; float r809 = r808.x; float r810 = r808.y; float r811 = r808.z; vec4 r812 = vec4(r809,r810,r811,r770); float r813 = r812.x; vec3 r814 = vec3(0.0,0.0,0.3); float r815 = r814.x; float r816 = r813-r815; float r817 = r812.y; vec3 r818 = vec3(0.0,0.0,0.3); float r819 = r818.y; float r820 = r817-r819; float r821 = r812.z; vec3 r822 = vec3(0.0,0.0,0.3); float r823 = r822.z; float r824 = r821-r823; float r825 = r812.w; vec3 r826 = vec3(r816,r820,r824); float r827 = length(r826); float r828 = 0.25; float r829 = r827-r828; float r830 = r812.x; vec3 r831 = vec3(0.0,0.2,0.3); float r832 = r831.x; float r833 = r830-r832; float r834 = r812.y; vec3 r835 = vec3(0.0,0.2,0.3); float r836 = r835.y; float r837 = r834-r836; float r838 = r812.z; vec3 r839 = vec3(0.0,0.2,0.3); float r840 = r839.z; float r841 = r838-r840; float r842 = r812.w; float r843 = 0.7853981633974483; vec3 r844 = vec3(0.0,1.0,0.0); vec3 r845 = vec3(r833,r837,r841); float r846 = cos(r843); vec3 r847 = r845*vec3(r846,r846,r846); float r848 = r844.y; float r849 = r845.z; float r850 = r848*r849; float r851 = r844.z; float r852 = r845.y; float r853 = r851*r852; float r854 = r850-r853; float r855 = r844.z; float r856 = r845.x; float r857 = r855*r856; float r858 = r844.x; float r859 = r845.z; float r860 = r858*r859; float r861 = r857-r860; float r862 = r844.x; float r863 = r845.y; float r864 = r862*r863; float r865 = r844.y; float r866 = r845.x; float r867 = r865*r866; float r868 = r864-r867; vec3 r869 = vec3(r854,r861,r868); float r870 = sin(r843); vec3 r871 = r869*vec3(r870,r870,r870); vec3 r872 = r847-r871; float r873 = 1.0; float r874 = cos(r843); float r875 = r873-r874; vec3 r876 = r845*vec3(r875,r875,r875); float r877 = dot(r844,r876); vec3 r878 = r844*vec3(r877,r877,r877); vec3 r879 = r872+r878; vec3 r880=r879; float r881 = r880.x; float r882 = r880.y; float r883 = r880.z; vec3 r884 = vec3(r881,r882,r883); float r885 = length(r884); float r886 = 0.135; float r887 = r885-r886; float r888 = -r887; float r889 = 0.004; float r890 = -r829; float r891 = -r888; float r892 = 0.0; bool r893 =(r890 == r892); float r894 = 0.5; float r895 = 0.5; float r896 = r891-r890; float r897 = r895*r896; float r898 = r897/r889; float r899 = r894+r898; float r900 = 0.0; float r901 = 1.0; float r902 = max(r899,r900); float r903 = min(r902,r901); float r904=r903; float r905 = 1.0; float r906 = r905-r904; float r907 = r891*r906; float r908 = r890*r904; float r909 = r907+r908; float r910 = r889*r904; float r911 = 1.0; float r912 = r911-r904; float r913 = r910*r912; float r914 = r909-r913; float r915 =(r893 ? r891 : r914); float r916 = -r915; float r917 = r812.x; vec3 r918 = vec3(0.0,0.2,0.3); float r919 = r918.x; float r920 = r917-r919; float r921 = r812.y; vec3 r922 = vec3(0.0,0.2,0.3); float r923 = r922.y; float r924 = r921-r923; float r925 = r812.z; vec3 r926 = vec3(0.0,0.2,0.3); float r927 = r926.z; float r928 = r925-r927; float r929 = r812.w; float r930 = 0.7853981633974483; vec3 r931 = vec3(0.0,1.0,0.0); vec3 r932 = vec3(r920,r924,r928); float r933 = cos(r930); vec3 r934 = r932*vec3(r933,r933,r933); float r935 = r931.y; float r936 = r932.z; float r937 = r935*r936; float r938 = r931.z; float r939 = r932.y; float r940 = r938*r939; float r941 = r937-r940; float r942 = r931.z; float r943 = r932.x; float r944 = r942*r943; float r945 = r931.x; float r946 = r932.z; float r947 = r945*r946; float r948 = r944-r947; float r949 = r931.x; float r950 = r932.y; float r951 = r949*r950; float r952 = r931.y; float r953 = r932.x; float r954 = r952*r953; float r955 = r951-r954; vec3 r956 = vec3(r941,r948,r955); float r957 = sin(r930); vec3 r958 = r956*vec3(r957,r957,r957); vec3 r959 = r934-r958; float r960 = 1.0; float r961 = cos(r930); float r962 = r960-r961; vec3 r963 = r932*vec3(r962,r962,r962); float r964 = dot(r931,r963); vec3 r965 = r931*vec3(r964,r964,r964); vec3 r966 = r959+r965; vec3 r967=r966; float r968 = r967.x; float r969 = r967.y; float r970 = r967.z; vec3 r971 = vec3(r968,r969,r970); float r972 = length(r971); float r973 = 0.135; float r974 = r972-r973; float r975 = abs(r974); float r976 = 0.0; float r977 = r975-r976; float r978 = r812[0]; float r979 = r812[1]; float r980 = r812[2]; float r981 = r812[3]; vec3 r982 = vec3(r978,r979,r980); vec3 r983 = vec3(0.0,-1.0,0.0); float r984 = dot(r982,r983); float r985 = -0.24; float r986 = r984-r985; float r987 = -r986; float r988 = 0.001; float r989 = -r977; float r990 = -r987; float r991 = 0.0; bool r992 =(r989 == r991); float r993 = 0.5; float r994 = 0.5; float r995 = r990-r989; float r996 = r994*r995; float r997 = r996/r988; float r998 = r993+r997; float r999 = 0.0; float r1000 = 1.0; float r1001 = max(r998,r999); float r1002 = min(r1001,r1000); float r1003=r1002; float r1004 = 1.0; float r1005 = r1004-r1003; float r1006 = r990*r1005; float r1007 = r989*r1003; float r1008 = r1006+r1007; float r1009 = r988*r1003; float r1010 = 1.0; float r1011 = r1010-r1003; float r1012 = r1009*r1011; float r1013 = r1008-r1012; float r1014 =(r992 ? r990 : r1013); float r1015 = -r1014; float r1016 = 0.01; float r1017 = 0.0; bool r1018 =(r916 == r1017); float r1019 = 0.5; float r1020 = 0.5; float r1021 = r1015-r916; float r1022 = r1020*r1021; float r1023 = r1022/r1016; float r1024 = r1019+r1023; float r1025 = 0.0; float r1026 = 1.0; float r1027 = max(r1024,r1025); float r1028 = min(r1027,r1026); float r1029=r1028; float r1030 = 1.0; float r1031 = r1030-r1029; float r1032 = r1015*r1031; float r1033 = r916*r1029; float r1034 = r1032+r1033; float r1035 = r1016*r1029; float r1036 = 1.0; float r1037 = r1036-r1029; float r1038 = r1035*r1037; float r1039 = r1034-r1038; float r1040 =(r1018 ? r1015 : r1039); float r1041 = r812.x; vec3 r1042 = vec3(0.0,0.0,0.3); float r1043 = r1042.x; float r1044 = r1041-r1043; float r1045 = r812.y; vec3 r1046 = vec3(0.0,0.0,0.3); float r1047 = r1046.y; float r1048 = r1045-r1047; float r1049 = r812.z; vec3 r1050 = vec3(0.0,0.0,0.3); float r1051 = r1050.z; float r1052 = r1049-r1051; float r1053 = r812.w; float r1054 = 0.0; vec2 r1055 = vec2(r1044,r1048); float r1056 = length(r1055); float r1057 = 0.25; float r1058 = r1056-r1057; float r1059 = 0.0; vec2 r1060 = vec2(r1058,r1052); float r1061 = length(r1060); float r1062 = 0.006; float r1063 = r1061-r1062; float r1064 = -r1063; float r1065 = 0.001; float r1066 = -r1040; float r1067 = -r1064; float r1068 = 0.0; bool r1069 =(r1066 == r1068); float r1070 = 0.5; float r1071 = 0.5; float r1072 = r1067-r1066; float r1073 = r1071*r1072; float r1074 = r1073/r1065; float r1075 = r1070+r1074; float r1076 = 0.0; float r1077 = 1.0; float r1078 = max(r1075,r1076); float r1079 = min(r1078,r1077); float r1080=r1079; float r1081 = 1.0; float r1082 = r1081-r1080; float r1083 = r1067*r1082; float r1084 = r1066*r1080; float r1085 = r1083+r1084; float r1086 = r1065*r1080; float r1087 = 1.0; float r1088 = r1087-r1080; float r1089 = r1086*r1088; float r1090 = r1085-r1089; float r1091 =(r1069 ? r1067 : r1090); float r1092 = -r1091; float r1093 = r812.x; vec3 r1094 = vec3(-0.1,-0.2,0.3); float r1095 = r1094.x; float r1096 = r1093-r1095; float r1097 = r812.y; vec3 r1098 = vec3(-0.1,-0.2,0.3); float r1099 = r1098.y; float r1100 = r1097-r1099; float r1101 = r812.z; vec3 r1102 = vec3(-0.1,-0.2,0.3); float r1103 = r1102.z; float r1104 = r1101-r1103; float r1105 = r812.w; vec4 r1106 = vec4(r1096,r1100,r1104,r1105); vec3 r1107 = r1106.xyz; vec3 r1108 = vec3(0.0,0.85065080835204,0.5257311121191336); float r1109 = dot(r1107,r1108); float r1110 = abs(r1109); vec3 r1111 = r1106.xyz; vec3 r1112 = vec3(0.0,-0.85065080835204,0.5257311121191336); float r1113 = dot(r1111,r1112); float r1114 = abs(r1113); vec3 r1115 = r1106.xyz; vec3 r1116 = vec3(0.5257311121191336,0.0,0.85065080835204); float r1117 = dot(r1115,r1116); float r1118 = abs(r1117); vec3 r1119 = r1106.xyz; vec3 r1120 = vec3(-0.5257311121191336,0.0,0.85065080835204); float r1121 = dot(r1119,r1120); float r1122 = abs(r1121); vec3 r1123 = r1106.xyz; vec3 r1124 = vec3(0.85065080835204,0.5257311121191336,0.0); float r1125 = dot(r1123,r1124); float r1126 = abs(r1125); vec3 r1127 = r1106.xyz; vec3 r1128 = vec3(-0.85065080835204,0.5257311121191336,0.0); float r1129 = dot(r1127,r1128); float r1130 = abs(r1129); float r1131 = max(r1110,max(r1114,max(r1118,max(r1122,max(r1126,r1130))))); float r1132 = 0.095; float r1133 = r1131-r1132; float r1134 = 0.1; float r1135 = 0.0; bool r1136 =(r1092 == r1135); float r1137 = 0.5; float r1138 = 0.5; float r1139 = r1133-r1092; float r1140 = r1138*r1139; float r1141 = r1140/r1134; float r1142 = r1137+r1141; float r1143 = 0.0; float r1144 = 1.0; float r1145 = max(r1142,r1143); float r1146 = min(r1145,r1144); float r1147=r1146; float r1148 = 1.0; float r1149 = r1148-r1147; float r1150 = r1133*r1149; float r1151 = r1092*r1147; float r1152 = r1150+r1151; float r1153 = r1134*r1147; float r1154 = 1.0; float r1155 = r1154-r1147; float r1156 = r1153*r1155; float r1157 = r1152-r1156; float r1158 =(r1136 ? r1133 : r1157); float r1159 = abs(r1158); float r1160 = 0.0; float r1161 = r1159-r1160; float r1162 = 2.0; float r1163 = r1161/r1162; float r1164 = min(r757,r1163); ;

  vec3 r1165 = vec3(0.0);
    
  float tim = mod(t*-0.25f+0.8, 1.0);
  if ( tim < 0.75 )
    r1164 = ORANGEMAP(r0.xyz/0.2)*.2;
  if ( tim < 0.5 )
    r1164 = sdMetaballs(r0.xyz, 4, t, ID);
  if ( tim < 0.25 )
    r1164 = (sin(r0.x*2.0)*0.1 +
             cos(r0.z*128.0)*0.001 + sdBox(r0.xyz, vec3(0.2)))*0.5;
  
  return vec4(r1164,r1165);
}

vec4 map ( vec4 r0 ) {
  float t = r0.w;
  // rotate and shift objects around a bit
  r0.x += r0.w*1.72;
  float ID = 0.0, ZID;
  r0.x = opMod1(r0.x, 1.0, ID);
  r0.y = opMod1(r0.y, 1.0, ZID);
  r0.w += sin(r0.w)*2.3*ZID;
  ID = ID*4096.0+ZID;
  opRotate(r0.xy, r0.w*sin(4.8*ID)*PI*(0.5+0.5*abs(cos(hash(ID)))));
  opRotate(r0.yz, sin(r0.w*sin(ID))*PI*(0.5 + 0.5*(cos(ID))));
  vec4 res = model(r0, t, ID);
  return vec4(res.x, res.y, ID, res.w);
}
// IBL image functions from https://www.shadertoy.com/view/ld3SRr 
vec3 Remap(vec3 c) {
  float HDR = smoothstep(2.8, 3.0, c.x+c.y+c.z);
  vec3 red_lit = mix(c,4.3*vec3(4.5, 2.5, 3.0), HDR);
  vec3 blu_lit = mix(c,1.8*vec3(2.0, 2.5, 3.0), HDR);
  return mix(red_lit, blu_lit, 0.4);
  return c;
}
vec3 Diffuse_Map(vec3 ray) {
  return pow(Remap(texture(iChannel0, ray.yzx).rgb), vec3(2.2));
}

vec3 Specular_Map(vec3 ray) {
  ray = ray.yzx;
  vec3 sharp = Remap(pow(texture(iChannel0, ray).rgb, vec3(2.2)));
  vec3 blur  = Remap(pow(texture(iChannel1, ray).rgb, vec3(1.7)));
  return pow(mix(sharp, blur, 0.3), vec3(1.2));
}

vec4 castRay( in vec3 ro, in vec3 rd, in float time ) {
  float tmin = 1.0;
  float tmax = 200.0;
   
  float t = tmin;
  vec3 c = vec3(-1.0,-1.0,-1.0);
  for (int i=0; i<200; i++) {
    float precis = 0.0005*t;
    vec4 res = map( vec4(ro+rd*t, time) );
    if (res.x < precis) {
      c = res.yzw;
      break;
    }
    t += res.x;
    if (t > tmax) break;
  }
  return vec4( t, c );
}
vec3 calcNormal( in vec3 pos, float time ) {
  vec2 e = vec2(1.0,-1.0)*0.5773*0.0005;
  return normalize( e.xyy*map( vec4(pos + e.xyy,time) ).x + 
            e.yyx*map( vec4(pos + e.yyx,time) ).x + 
            e.yxy*map( vec4(pos + e.yxy,time) ).x + 
            e.xxx*map( vec4(pos + e.xxx,time) ).x );
}

vec2 hash2( vec2 n ) {
  return fract(sin(vec2( n.x*n.y, n.x+n.y))*vec2(225.1459123,312.3490423));
}

// --------------------------- sampling funcs ---------------------------------
// from PBRT 3rd ed
vec2 Concentric_Sample_Disk(vec2 rand) {
  vec2 u = hash2(rand);
  // maps to [-1, 1]^2
  vec2 offset = 2.0 * u - vec2(1.0, 1.0);
  if ( offset.x == 0.0 && offset.y == 0.0 )
  return vec2(0.0);

  float theta, r;
  if ( abs(offset.x) > abs(offset.y) ) {
  r = offset.x;
  theta = (PI/4.0) * (offset.y/offset.x);
  } else {
  r = offset.y;
  theta = (PI/2.0) * (offset.x/offset.y);
  }
  return r * vec2(cos(theta), sin(theta));
}

// from PBRT 3rd ed
vec3 Sample_Cosine_Hemisphere ( vec2 rand, inout float pdf ) {
  vec2 d = Concentric_Sample_Disk(rand);
  float phi = sqrt(max(0.0, 1.0 - d.x*d.x - d.y*d.y));
  pdf = 2.0 * phi * (1.0/PI);
  return vec3(d.x, d.y, phi);
}
vec3 Sample_Uniform_Cone ( vec2 rand, float theta ) {
  vec2 r = hash2(rand);
  float cos_theta = 1.0 - r.x + r.x*theta;
  float sin_theta = sqrt(1.0 - cos_theta*cos_theta);
  float phi = r.y * 2.0 * PI;
  return vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
}
vec3 Reorient_Hemisphere ( vec3 wo, vec3 N ) {
  vec3 binormal = (abs(N.x) < 1.0 ? vec3(1.0, 0.0, 0.0) :
                   vec3(0.0, 1.0, 0.0));
  binormal = normalize(cross(N, binormal));
  vec3 bitangent = cross(binormal, N);
  return bitangent*wo.x + binormal*wo.y + wo.z*N;
}
vec3 Sample_Cosine_Hemisphere_N ( vec3 N , vec2 rand, inout float pdf ) {
  vec3 wo = normalize(Sample_Cosine_Hemisphere(rand, pdf));
  return Reorient_Hemisphere(wo, N);
}
vec3 Sample_Uniform_Cone_N ( vec3 wi, vec3 N , float theta, vec2 rand ) {
  vec3 NN = reflect(wi, N);
  vec3 wo = normalize(Sample_Uniform_Cone(rand, theta));
  return Reorient_Hemisphere(wo, N);
}

float Uniform_Cone_PDF ( float theta ) {
  return (2.0 * PI * (1.0 - theta));
}

vec3 BRDF_Diffuse_F ( vec3 wi, vec3 N, vec2 rand, inout float pdf ) {
  return Sample_Cosine_Hemisphere_N(N, rand, pdf);
}

vec3 BRDF_Glossy_F ( vec3 wi, vec3 N, vec2 rand, float m_glossy_lobe,
                     inout float pdf ) {
  vec3 res = Sample_Uniform_Cone_N(wi, N, m_glossy_lobe, rand);
  pdf = Uniform_Cone_PDF(m_glossy_lobe);
  return res;
}

vec3 BRDF_Specular_F ( vec3 wi, vec3 N, vec2 rand, inout float pdf ) {
  pdf = 1.0;
  return reflect(wi, N);
}

// --------------------------- BRDF -------------------------------------------

float SchlickFresnel(float u) {
  float f = clamp(1.0 - u, 0.0, 1.0);
  float f2 = f*f;
  return f2*f2*f; // f^5
}

float smithG_GGX_Correlated ( float L, float R, float a2 ) {
  return L * sqrt((R - a2*sqr(R) + a2));
}
vec3 BRDF ( vec3 wi, vec3 wo, vec3 N, int ID, float t ) {
  // materials
  float mhash = hash(t);
  vec3 colour = vec3(mhash, hash(mhash), hash(hash(mhash)));
  float m_roughness     = Rm_roughness    (t);
  float m_metallic    = Rm_metallic     (t);
  float m_fresnel     = Rm_fresnel    (t);
  float m_diffuse     = Rm_diffuse    (t);
  float m_glossy_lobe   = Rm_glossy_lobe  (t);
  float m_subsurface    = Rm_subsurface   (t);
  float m_anisotropic   = Rm_anisotropic  (t);
  
  // get binormal/bitangent, half vec and L/V
  vec3 X = (abs(N.x) < 1.0 ? vec3(1.0, 0.0, 0.0) :
                             vec3(0.0, 1.0, 0.0));
  X = normalize(cross(N, X));
  vec3 Y = cross(X, N),
       L = wo,
       V = -wi,
       H = normalize(L+V);
  
  float cosNV = dot(N, V),
        cosNL = dot(N, L),
        cosNH = dot(N, H),
        cosHL = dot(H, L),
        cosHV = dot(H, V);
  
  // --------- diffuse, I just made it up 
  vec3 diffuse = colour * pow(dot(H, L)*dot(N, V), 0.5) * (1.0/(PI));
  // --------- microfacet
  vec3 F, G, D;
  float FL = SchlickFresnel(cosNL),
        FV = SchlickFresnel(cosNV);

  {// ---------- Fresnel
    // fresnel diffuse from disney, modified to use lambertian and f0
    float f0 = m_fresnel * (m_metallic);
    float Fd90 = f0 * dot(H, L)*dot(H, L);
    
    F = (1.0 - f0) * diffuse +
      (vec3(mix(1.0, Fd90, dot(H, L)) * mix(1.0, Fd90, FV)));
  }
  {// ---------- Geometry
    // Heits 2014, SmithGGXCorrelated with half vec combined with  anisotropic
    // from disney's GTR2_aniso
    float a2 = (0.5 + (m_roughness));
    float aspect = sqrt(1.0 - m_anisotropic*0.9);
    float ax = a2/aspect,
          ay = a2*aspect;
    float GGXV = smithG_GGX_Correlated(cosHL, cosNV, ax),
          GGXL = smithG_GGX_Correlated(cosNV, cosHL, ay);
    G = vec3(0.5 /  (GGXV*ax + GGXL*ay));
  }
  {// ---------- Distribution
    // Hyper-Cauchy distribution
    float param = 1.2+m_anisotropic,
          shp = (1.1 - (m_roughness)*0.55);
    float tanHL = length(cross(H, L))/dot(H, L);
    D = ((param - 1.0)*pow(sqrt(2.0), (2.0*param - 2.0)))/
        (PI*sqr(shp) * pow(cosHL, 4.0) *
         pow(2.0 + (sqr(tanHL)/sqr(shp)), param)) * vec3(1.0);
  }
  // If the microfacet is described with H, then the following energy
  // conservation model may be used [Edwards et al. 2006] with reciprocity
  // [Neumann et al. 1999]
  vec3 microfacet = F*G*D/(1.0 * dot(H, V) * max(dot(N, L), dot(N, V)));
  
  // --------- subsurface
  // based off the Henyey-Greenstein
  float r_term = (0.0 + 0.7*(1.0-m_roughness));
  float m_term = 0.2+(m_subsurface)*1.0;
  float rr_term = m_term * (1.0 - sqr(r_term))/(4.0*PI);
  rr_term *= 1.0/(pow(1.0 + sqr(r_term) - 2.0*r_term*dot(H, L), 3.0/2.0));
  vec3 ss_diff =  diffuse*vec3(1.0)*rr_term *
                 (FL + FV + FL*FV*(rr_term));
  vec3 rd = mix(diffuse, pow(ss_diff, vec3(0.5)), m_subsurface*0.5);
    
  return (microfacet + rd);
}


vec4 BRDF_F ( vec3 wi, vec3 N, int scount, vec2 rand, float t,
              inout float pdf ) {
  float fsample = float(hash(vec2(scount, scount)));
  int ID = -1;
  float remainder = 1.0;
  remainder -= Rm_diffuse(t);
  float rem_t = hash(rand);
  vec3 brdf_f;
  if ( remainder < rem_t ) {
    brdf_f = BRDF_Diffuse_F(wi, N, rand, pdf);
    ID = 0;
  }
  remainder -= Rm_specular(t);
  if ( ID == -1 && remainder < rem_t ) {
    brdf_f = BRDF_Specular_F(wi, N, rand, pdf);
    ID = 1;
  }
  remainder -= Rm_glossy(t);
  if ( ID == -1 && remainder < rem_t ) {
    brdf_f = BRDF_Glossy_F(wi, N, rand, Rm_glossy_lobe(t), pdf);
    ID = 2;
  }
  return vec4(brdf_f, float(ID));
}

// --------------------------- render -----------------------------------------


vec3 render( in vec3 ro, in vec3 rd, vec2 rand, float time) {
  vec4 res = castRay(ro,rd, time);
  float dist = res.x;
  vec3 mat_col = res.yzw;
  vec3 O = ro + dist*rd;
  vec3 N = calcNormal( O , time );
  vec3 col = vec3(0.0);
  rand *= hash(iGlobalTime);
  float t = res.z ;

  for ( int i = 0; i != SAMPLES; ++ i ) {
    vec3 wi = rd;
    float pdf;
    vec4 brdf_f = BRDF_F(wi, N, int(sin(float(iFrame))*15.0), rand, t, pdf);
    vec3 L = brdf_f.xyz;
    int ID = int(brdf_f.w);

    rand.x = 2.0 * hash(rand);
    rand.y = 2.0 * hash(rand);

    float mult;
    if ( ID == 0 ) mult = Rm_diffuse(t);
    if ( ID == 1 ) mult = Rm_specular(t);
    if ( ID == 2 ) mult = Rm_glossy(t);
    
    if( mat_col.x>=0.0 ) {
      vec3 brdf = BRDF(wi, L, N, ID, t);
    
      // F_s(wi, wo, N) * L_o(P, wo) * cosTheta_o do
      col += brdf /(pdf*(1.0/mult))  * dot(L, N) * Specular_Map(L) ;
      col += sqr(length(ro-O))*0.0001; // light fog
    } else {
      col += Diffuse_Map(rd)*1.5;
    }
  }


  col /= vec3(SAMPLES);

  return vec3( clamp(col,0.0,1.0) );
}
// Create a matrix to transform coordinates to look towards a given point.
// * `eye` is the position of the camera.
// * `centre` is the position to look towards.
// * `up` is the 'up' direction.
mat3 look_at(vec3 eye, vec3 centre, vec3 up) {
  vec3 ww = normalize(centre - eye);
  vec3 uu = normalize(cross(ww, up));
  vec3 vv = normalize(cross(uu, ww));
  return mat3(uu, vv, ww);
}
// Generate a ray direction for ray-casting.
// * `camera` is the camera look-at matrix.
// * `pos` is the screen position, normally in the range -1..1
// * `lens` is the lens length of the camera (encodes field-of-view).
//   0 is very wide, and 2 is a good default.
vec3 ray_direction(mat3 camera, vec2 pos, float lens) {
  return normalize(camera * vec3(pos, lens));
}
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  vec2 p = -1.0 + 2.0 * fragCoord.xy / iResolution.xy;
  p.x *= iResolution.x/iResolution.y;

  float time = 3.5 + iMouse.x/iResolution.x*64.0*PI*SPEED;
  bool refresh = false;
  if ( REVOLVE ) {
    time = iGlobalTime*0.25;
    refresh = true;
  }
  if ( GALLERY ) {
    float timr = 22.0+iTime*0.1;
    time = floor(mod(timr, 150.0));
  }
  vec3 eye = vec3(sin(time*0.2)*5.0,
          cos(time*0.2)*5.0,
          0.6 + (pow(abs(sin(time*0.4)), 2.0))*0.6);
  vec3 centre = vec3(0.0, 0.0, 0.1);
  vec3 up = vec3(0.0, 0.0, 1.0);
  mat3 camera = look_at(eye, centre, up);
  vec3 dir = ray_direction(camera, p, 2.5);

  vec2 rand = vec2(gl_FragCoord.x + iTime,
             gl_FragCoord.y + sin(float(iFrame)));

  vec3 col = render( eye, dir, rand, time );

  vec2 xy = gl_FragCoord.xy/iResolution.xy;
  
  float fDelta = float(iFrame) - texture(iChannel2, vec2(0.5)).z;
  
  if ( fDelta < 1.0 || refresh )
    fragColor = vec4(col,1.0);
  else
    fragColor = (texture(iChannel3, xy) * fDelta + vec4(col,1.0))
           / (fDelta+1.0);
}


