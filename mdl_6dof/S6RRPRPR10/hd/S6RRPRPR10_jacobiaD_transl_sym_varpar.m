% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR10
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR10_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR10_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:53
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (103->31), mult. (307->49), div. (0->0), fcn. (282->8), ass. (0->26)
	t199 = sin(pkin(11));
	t200 = sin(pkin(6));
	t201 = cos(pkin(11));
	t220 = t200 * (r_i_i_C(1) * t199 + r_i_i_C(2) * t201 + pkin(8));
	t219 = r_i_i_C(3) + qJ(3);
	t203 = sin(qJ(2));
	t204 = sin(qJ(1));
	t218 = t204 * t203;
	t205 = cos(qJ(2));
	t217 = t204 * t205;
	t206 = cos(qJ(1));
	t216 = t206 * t203;
	t215 = t206 * t205;
	t214 = qJD(2) * t203;
	t202 = cos(pkin(6));
	t213 = t202 * t218;
	t212 = t202 * t215;
	t211 = qJD(2) * t202 + qJD(1);
	t210 = t201 * r_i_i_C(1) - t199 * r_i_i_C(2) + pkin(2);
	t208 = t202 * t217 + t216;
	t207 = t202 * t216 + t217;
	t194 = -qJD(1) * t213 - t204 * t214 + t211 * t215;
	t193 = t208 * qJD(1) + t207 * qJD(2);
	t192 = t207 * qJD(1) + t208 * qJD(2);
	t191 = -qJD(1) * t212 - qJD(2) * t215 + t211 * t218;
	t1 = [-(-t212 + t218) * qJD(3) - t219 * t193 - t210 * t194 + (-t206 * pkin(1) - t204 * t220) * qJD(1), -(t213 - t215) * qJD(3) - t219 * t192 + t210 * t191, -t191, 0, 0, 0; t208 * qJD(3) - t219 * t191 - t210 * t192 + (-t204 * pkin(1) + t206 * t220) * qJD(1), t207 * qJD(3) - t210 * t193 + t219 * t194, t193, 0, 0, 0; 0, (t203 * qJD(3) + (-t210 * t203 + t219 * t205) * qJD(2)) * t200, t200 * t214, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:54
	% EndTime: 2019-10-10 10:20:54
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (222->61), mult. (494->100), div. (0->0), fcn. (467->10), ass. (0->45)
	t259 = sin(qJ(1));
	t256 = cos(pkin(6));
	t267 = qJD(2) * t256 + qJD(1);
	t258 = sin(qJ(2));
	t280 = t259 * t258;
	t272 = t256 * t280;
	t275 = qJD(2) * t258;
	t260 = cos(qJ(2));
	t261 = cos(qJ(1));
	t277 = t261 * t260;
	t239 = -qJD(1) * t272 - t259 * t275 + t267 * t277;
	t255 = sin(pkin(6));
	t281 = t255 * t261;
	t268 = qJD(4) * t281;
	t285 = t239 - t268;
	t253 = pkin(11) + qJ(4);
	t251 = sin(t253);
	t252 = cos(t253);
	t266 = r_i_i_C(1) * t251 + r_i_i_C(2) * t252;
	t263 = qJD(4) * t266;
	t284 = r_i_i_C(3) + pkin(9) + qJ(3);
	t283 = t255 * t258;
	t282 = t255 * t259;
	t279 = t259 * t260;
	t278 = t261 * t258;
	t276 = qJD(1) * t255;
	t274 = qJD(2) * t260;
	t241 = t256 * t278 + t279;
	t273 = qJD(4) * t241;
	t271 = t256 * t277;
	t270 = pkin(3) * sin(pkin(11)) + pkin(8);
	t269 = t261 * t276;
	t250 = cos(pkin(11)) * pkin(3) + pkin(2);
	t265 = t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + t250;
	t264 = t256 * t279 + t278;
	t262 = t259 * t276 - t273;
	t244 = t252 * t268;
	t243 = -t272 + t277;
	t240 = -t271 + t280;
	t238 = t264 * qJD(1) + t241 * qJD(2);
	t237 = t241 * qJD(1) + t264 * qJD(2);
	t236 = -qJD(1) * t271 - t261 * t274 + t267 * t280;
	t235 = t251 * t269 - t237 * t252 + (-t243 * t251 + t252 * t282) * qJD(4);
	t234 = t252 * t269 + t237 * t251 + (-t243 * t252 - t251 * t282) * qJD(4);
	t1 = [(-t239 * t252 + t251 * t273 + t244) * r_i_i_C(1) + (t285 * t251 + t252 * t273) * r_i_i_C(2) - t239 * t250 - t240 * qJD(3) - t284 * t238 + (-t261 * pkin(1) + (-t266 - t270) * t282) * qJD(1), t243 * qJD(3) + t265 * t236 - t284 * t237 + t263 * t264, -t236, t234 * r_i_i_C(1) - t235 * r_i_i_C(2), 0, 0; t235 * r_i_i_C(1) + t234 * r_i_i_C(2) + t264 * qJD(3) - t237 * t250 - t284 * t236 + (-pkin(1) * t259 + t270 * t281) * qJD(1), t241 * qJD(3) - t265 * t238 + t284 * t239 + t240 * t263, t238, t244 * r_i_i_C(2) + (t262 * r_i_i_C(1) - t239 * r_i_i_C(2)) * t252 + (-t285 * r_i_i_C(1) - t262 * r_i_i_C(2)) * t251, 0, 0; 0, (qJD(3) * t258 - t260 * t263 + (-t265 * t258 + t284 * t260) * qJD(2)) * t255, t255 * t275, -t266 * t255 * t274 + ((-t251 * t256 - t252 * t283) * r_i_i_C(1) + (t251 * t283 - t252 * t256) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:54
	% EndTime: 2019-10-10 10:20:55
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (429->69), mult. (881->107), div. (0->0), fcn. (857->10), ass. (0->51)
	t298 = pkin(11) + qJ(4);
	t296 = sin(t298);
	t297 = cos(t298);
	t335 = r_i_i_C(3) + qJ(5);
	t337 = pkin(4) - r_i_i_C(2);
	t342 = (t337 * t296 - t335 * t297) * qJD(4) - qJD(5) * t296;
	t304 = sin(qJ(1));
	t301 = cos(pkin(6));
	t316 = qJD(2) * t301 + qJD(1);
	t303 = sin(qJ(2));
	t330 = t304 * t303;
	t322 = t301 * t330;
	t325 = qJD(2) * t303;
	t305 = cos(qJ(2));
	t306 = cos(qJ(1));
	t327 = t306 * t305;
	t281 = -qJD(1) * t322 - t304 * t325 + t316 * t327;
	t328 = t306 * t303;
	t329 = t304 * t305;
	t286 = t301 * t328 + t329;
	t300 = sin(pkin(6));
	t331 = t300 * t306;
	t315 = t286 * t296 + t297 * t331;
	t326 = qJD(1) * t300;
	t319 = t304 * t326;
	t341 = t315 * qJD(4) - t281 * t297 - t296 * t319;
	t314 = -t286 * t297 + t296 * t331;
	t340 = t314 * qJD(4) - t281 * t296 + t297 * t319;
	t295 = cos(pkin(11)) * pkin(3) + pkin(2);
	t338 = t335 * t296 + t337 * t297 + t295;
	t336 = r_i_i_C(1) + pkin(9) + qJ(3);
	t288 = -t322 + t327;
	t334 = t288 * t296;
	t333 = t300 * t303;
	t332 = t300 * t304;
	t324 = qJD(2) * t305;
	t321 = t301 * t327;
	t320 = pkin(3) * sin(pkin(11)) + pkin(8);
	t318 = t306 * t326;
	t317 = t300 * t324;
	t313 = t288 * t297 + t296 * t332;
	t312 = t301 * t296 + t297 * t333;
	t311 = t301 * t329 + t328;
	t285 = -t321 + t330;
	t282 = t312 * qJD(4) + t296 * t317;
	t280 = t311 * qJD(1) + t286 * qJD(2);
	t279 = t286 * qJD(1) + t311 * qJD(2);
	t278 = -qJD(1) * t321 - t306 * t324 + t316 * t330;
	t273 = t296 * t318 - qJD(4) * t334 + (qJD(4) * t332 - t279) * t297;
	t272 = t313 * qJD(4) - t279 * t296 - t297 * t318;
	t1 = [-t315 * qJD(5) - t281 * t295 - t285 * qJD(3) - t336 * t280 + t337 * t341 + t335 * t340 + (-t306 * pkin(1) - t320 * t332) * qJD(1), t288 * qJD(3) + t338 * t278 - t336 * t279 + t311 * t342, -t278, t313 * qJD(5) - t337 * t272 + t335 * t273, t272, 0; -(t297 * t332 - t334) * qJD(5) - t279 * t295 + t311 * qJD(3) - t336 * t278 + t337 * t273 + t335 * t272 + (-t304 * pkin(1) + t320 * t331) * qJD(1), t286 * qJD(3) - t280 * t338 + t336 * t281 + t342 * t285, t280, -t314 * qJD(5) - t335 * t341 + t337 * t340, -t340, 0; 0, (qJD(3) * t303 - t342 * t305 + (-t303 * t338 + t336 * t305) * qJD(2)) * t300, t300 * t325, t312 * qJD(5) + t335 * (t297 * t317 + (-t296 * t333 + t297 * t301) * qJD(4)) - t337 * t282, t282, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:55
	% EndTime: 2019-10-10 10:20:56
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (806->98), mult. (1673->158), div. (0->0), fcn. (1681->12), ass. (0->66)
	t388 = sin(qJ(2));
	t389 = sin(qJ(1));
	t391 = cos(qJ(2));
	t392 = cos(qJ(1));
	t431 = cos(pkin(6));
	t408 = t392 * t431;
	t368 = t388 * t408 + t389 * t391;
	t383 = pkin(11) + qJ(4);
	t381 = sin(t383);
	t382 = cos(t383);
	t385 = sin(pkin(6));
	t422 = t385 * t392;
	t440 = t368 * t382 - t381 * t422;
	t406 = t391 * t408;
	t421 = t388 * t389;
	t367 = -t406 + t421;
	t387 = sin(qJ(6));
	t390 = cos(qJ(6));
	t436 = t368 * t381 + t382 * t422;
	t439 = t367 * t390 + t387 * t436;
	t438 = t367 * t387 - t390 * t436;
	t405 = t390 * r_i_i_C(1) - t387 * r_i_i_C(2);
	t395 = t405 * qJD(6) + qJD(5);
	t404 = -t387 * r_i_i_C(1) - t390 * r_i_i_C(2);
	t402 = qJ(5) - t404;
	t416 = pkin(4) + pkin(10) + r_i_i_C(3);
	t393 = (t416 * t381 - t402 * t382) * qJD(4) - t395 * t381;
	t403 = qJD(2) * t431 + qJD(1);
	t409 = t389 * t431;
	t407 = t388 * t409;
	t418 = qJD(2) * t388;
	t420 = t392 * t391;
	t356 = -qJD(1) * t407 - t389 * t418 + t403 * t420;
	t419 = qJD(1) * t385;
	t413 = t389 * t419;
	t435 = qJD(4) * t436 - t356 * t382 - t381 * t413;
	t380 = cos(pkin(11)) * pkin(3) + pkin(2);
	t433 = t402 * t381 + t416 * t382 + t380;
	t432 = -pkin(5) - pkin(9) - qJ(3);
	t370 = -t407 + t420;
	t426 = t370 * t381;
	t425 = t385 * t388;
	t424 = t385 * t389;
	t423 = t385 * t391;
	t417 = qJD(2) * t391;
	t414 = pkin(3) * sin(pkin(11)) + pkin(8);
	t412 = t392 * t419;
	t411 = t385 * t417;
	t410 = t385 * t418;
	t400 = t370 * t382 + t381 * t424;
	t399 = t405 - t432;
	t369 = t392 * t388 + t391 * t409;
	t365 = t381 * t425 - t431 * t382;
	t397 = t431 * t381 + t382 * t425;
	t396 = t404 * qJD(6) + qJD(3);
	t349 = t440 * qJD(4) + t356 * t381 - t382 * t413;
	t362 = -t382 * t424 + t426;
	t357 = t397 * qJD(4) + t381 * t411;
	t355 = t369 * qJD(1) + t368 * qJD(2);
	t354 = t368 * qJD(1) + t369 * qJD(2);
	t353 = -qJD(1) * t406 - t392 * t417 + t403 * t421;
	t348 = t381 * t412 - qJD(4) * t426 + (qJD(4) * t424 - t354) * t382;
	t347 = t400 * qJD(4) - t354 * t381 - t382 * t412;
	t346 = t347 * t387 - t353 * t390 + (t362 * t390 - t369 * t387) * qJD(6);
	t345 = t347 * t390 + t353 * t387 + (-t362 * t387 - t369 * t390) * qJD(6);
	t1 = [-t367 * qJD(3) - t436 * qJD(5) - t356 * t380 - t402 * t349 - t399 * t355 + (t438 * r_i_i_C(1) + t439 * r_i_i_C(2)) * qJD(6) + t416 * t435 + (-t392 * pkin(1) - t414 * t424) * qJD(1), t433 * t353 - t399 * t354 + t393 * t369 + t396 * t370, -t353, -t416 * t347 + t402 * t348 + t395 * t400, t347, r_i_i_C(1) * t345 - r_i_i_C(2) * t346; t346 * r_i_i_C(1) + t345 * r_i_i_C(2) + t347 * qJ(5) + t369 * qJD(3) + t362 * qJD(5) - t354 * t380 + t432 * t353 + t416 * t348 + (-pkin(1) * t389 + t414 * t422) * qJD(1), -t355 * t433 + t399 * t356 + t393 * t367 + t396 * t368, t355, -t416 * t349 + t395 * t440 - t402 * t435, t349, (t349 * t390 - t355 * t387) * r_i_i_C(1) + (-t349 * t387 - t355 * t390) * r_i_i_C(2) + (-t439 * r_i_i_C(1) + t438 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t433 + t396) * t388 + (t399 * qJD(2) - t393) * t391) * t385, t410, t395 * t397 + t402 * (-t365 * qJD(4) + t382 * t411) - t416 * t357, t357, (t357 * t390 - t387 * t410) * r_i_i_C(1) + (-t357 * t387 - t390 * t410) * r_i_i_C(2) + ((-t365 * t387 + t390 * t423) * r_i_i_C(1) + (-t365 * t390 - t387 * t423) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end