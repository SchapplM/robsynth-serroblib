% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR9
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR9_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR9_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR9_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
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
	% StartTime: 2019-10-10 12:44:10
	% EndTime: 2019-10-10 12:44:10
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
	% StartTime: 2019-10-10 12:44:11
	% EndTime: 2019-10-10 12:44:11
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (145->51), mult. (441->91), div. (0->0), fcn. (412->8), ass. (0->38)
	t242 = cos(pkin(6));
	t244 = sin(qJ(2));
	t245 = sin(qJ(1));
	t264 = t245 * t244;
	t258 = t242 * t264;
	t247 = cos(qJ(2));
	t248 = cos(qJ(1));
	t261 = t248 * t247;
	t233 = -qJD(1) * t258 - qJD(2) * t264 + (qJD(2) * t242 + qJD(1)) * t261;
	t241 = sin(pkin(6));
	t265 = t241 * t248;
	t256 = qJD(3) * t265;
	t269 = t233 - t256;
	t268 = -r_i_i_C(3) - pkin(9);
	t267 = t241 * t245;
	t246 = cos(qJ(3));
	t266 = t241 * t246;
	t263 = t245 * t247;
	t262 = t248 * t244;
	t260 = qJD(1) * t241;
	t235 = t242 * t262 + t263;
	t259 = qJD(3) * t235;
	t257 = t248 * t260;
	t243 = sin(qJ(3));
	t255 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
	t254 = t246 * r_i_i_C(1) - t243 * r_i_i_C(2) + pkin(2);
	t253 = t242 * t261 - t264;
	t252 = t242 * t263 + t262;
	t251 = t258 - t261;
	t250 = qJD(3) * t255;
	t249 = t245 * t260 - t259;
	t238 = t246 * t256;
	t232 = t252 * qJD(1) + t235 * qJD(2);
	t231 = t235 * qJD(1) + t252 * qJD(2);
	t230 = -t253 * qJD(1) + t251 * qJD(2);
	t229 = t243 * t257 - t231 * t246 + (t243 * t251 + t245 * t266) * qJD(3);
	t228 = t246 * t257 + t231 * t243 + (-t243 * t267 + t246 * t251) * qJD(3);
	t1 = [(-t233 * t246 + t243 * t259 + t238) * r_i_i_C(1) + (t269 * t243 + t246 * t259) * r_i_i_C(2) - t233 * pkin(2) + t268 * t232 + (-t248 * pkin(1) + (-pkin(8) + t255) * t267) * qJD(1), t254 * t230 + t268 * t231 - t252 * t250, t228 * r_i_i_C(1) - t229 * r_i_i_C(2), 0, 0, 0; -t231 * pkin(2) + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) + t268 * t230 + (-pkin(1) * t245 + pkin(8) * t265) * qJD(1), -t254 * t232 - t268 * t233 + t253 * t250, t238 * r_i_i_C(2) + (t249 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t246 + (-t269 * r_i_i_C(1) - t249 * r_i_i_C(2)) * t243, 0, 0, 0; 0, (t247 * t250 + (-t254 * t244 - t268 * t247) * qJD(2)) * t241, t255 * t247 * t241 * qJD(2) + ((-t242 * t243 - t244 * t266) * r_i_i_C(1) + (t241 * t243 * t244 - t242 * t246) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:11
	% EndTime: 2019-10-10 12:44:11
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (332->69), mult. (663->114), div. (0->0), fcn. (618->10), ass. (0->55)
	t280 = cos(pkin(6));
	t282 = sin(qJ(2));
	t283 = sin(qJ(1));
	t314 = t283 * t282;
	t302 = t280 * t314;
	t285 = cos(qJ(2));
	t286 = cos(qJ(1));
	t311 = t286 * t285;
	t266 = -qJD(1) * t302 - qJD(2) * t314 + (qJD(2) * t280 + qJD(1)) * t311;
	t277 = qJD(3) + qJD(4);
	t279 = sin(pkin(6));
	t315 = t279 * t286;
	t322 = t277 * t315 - t266;
	t321 = r_i_i_C(3) + pkin(10) + pkin(9);
	t320 = pkin(3) * qJD(3);
	t312 = t286 * t282;
	t313 = t283 * t285;
	t268 = t280 * t312 + t313;
	t319 = t268 * t277;
	t278 = qJ(3) + qJ(4);
	t275 = sin(t278);
	t318 = t275 * t277;
	t317 = t279 * t283;
	t284 = cos(qJ(3));
	t316 = t279 * t284;
	t276 = cos(t278);
	t292 = t302 - t311;
	t306 = qJD(1) * t286;
	t300 = t279 * t306;
	t290 = t277 * t292 + t300;
	t293 = t280 * t313 + t312;
	t264 = t268 * qJD(1) + t293 * qJD(2);
	t296 = t277 * t317 - t264;
	t259 = -t296 * t275 + t290 * t276;
	t260 = t290 * t275 + t296 * t276;
	t310 = t259 * r_i_i_C(1) - t260 * r_i_i_C(2);
	t307 = qJD(1) * t283;
	t301 = t279 * t307;
	t291 = t301 - t319;
	t298 = t322 * t276;
	t309 = (t322 * t275 + t291 * t276) * r_i_i_C(1) + (-t291 * t275 + t298) * r_i_i_C(2);
	t299 = qJD(2) * t279 * t285;
	t289 = -t277 * t280 - t299;
	t304 = t277 * t279 * t282;
	t308 = (t289 * t275 - t276 * t304) * r_i_i_C(1) + (t275 * t304 + t289 * t276) * r_i_i_C(2);
	t281 = sin(qJ(3));
	t305 = t281 * t320;
	t297 = -r_i_i_C(1) * t275 - r_i_i_C(2) * t276;
	t274 = t284 * pkin(3) + pkin(2);
	t295 = t276 * r_i_i_C(1) - t275 * r_i_i_C(2) + t274;
	t294 = t280 * t311 - t314;
	t288 = t297 * t277 - t305;
	t265 = t293 * qJD(1) + t268 * qJD(2);
	t263 = -t294 * qJD(1) + t292 * qJD(2);
	t1 = [(t268 * t318 + t298) * r_i_i_C(1) + (t266 * t275 + t276 * t319) * r_i_i_C(2) - t266 * t274 + t268 * t305 - pkin(1) * t306 - t321 * t265 + ((-r_i_i_C(2) * t318 + t284 * t320) * t286 + (-pkin(3) * t281 - pkin(8) + t297) * t307) * t279, t295 * t263 - t321 * t264 - t288 * t293, (t284 * t300 + t264 * t281 + (-t281 * t317 + t284 * t292) * qJD(3)) * pkin(3) + t310, t310, 0, 0; t260 * r_i_i_C(1) + t259 * r_i_i_C(2) - t264 * t274 - t321 * t263 + (-pkin(1) * t283 + pkin(8) * t315) * qJD(1) + (t281 * t300 + (t281 * t292 + t283 * t316) * qJD(3)) * pkin(3), -t295 * t265 + t321 * t266 + t288 * t294, (t284 * t301 - t266 * t281 + (-t268 * t284 + t281 * t315) * qJD(3)) * pkin(3) + t309, t309, 0, 0; 0, (t288 * t285 + (-t295 * t282 + t321 * t285) * qJD(2)) * t279, (-t281 * t299 + (-t280 * t281 - t282 * t316) * qJD(3)) * pkin(3) + t308, t308, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:13
	% EndTime: 2019-10-10 12:44:13
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (797->95), mult. (1447->150), div. (0->0), fcn. (1413->12), ass. (0->62)
	t431 = r_i_i_C(3) + qJ(5);
	t387 = sin(pkin(12));
	t389 = cos(pkin(12));
	t406 = r_i_i_C(1) * t389 - r_i_i_C(2) * t387 + pkin(4);
	t394 = cos(qJ(3));
	t382 = t394 * pkin(3) + pkin(2);
	t386 = qJ(3) + qJ(4);
	t383 = sin(t386);
	t384 = cos(t386);
	t434 = t431 * t383 + t406 * t384 + t382;
	t385 = qJD(3) + qJD(4);
	t430 = t383 * t385;
	t429 = t384 * t385;
	t388 = sin(pkin(6));
	t392 = sin(qJ(2));
	t428 = t388 * t392;
	t393 = sin(qJ(1));
	t427 = t388 * t393;
	t426 = t388 * t394;
	t396 = cos(qJ(1));
	t425 = t388 * t396;
	t424 = t393 * t392;
	t395 = cos(qJ(2));
	t423 = t393 * t395;
	t422 = t396 * t392;
	t421 = t396 * t395;
	t420 = qJD(1) * t388;
	t419 = qJD(2) * t392;
	t418 = qJD(2) * t395;
	t417 = t384 * t428;
	t390 = cos(pkin(6));
	t416 = t390 * t424;
	t415 = t383 * t425;
	t414 = t384 * t425;
	t413 = t390 * t421;
	t412 = t393 * t420;
	t411 = t396 * t420;
	t410 = t388 * t418;
	t408 = qJD(2) * t390 + qJD(1);
	t363 = -qJD(1) * t416 - t393 * t419 + t408 * t421;
	t409 = -t363 * t384 + t385 * t414;
	t370 = t390 * t422 + t423;
	t404 = t390 * t423 + t422;
	t361 = t370 * qJD(1) + t404 * qJD(2);
	t407 = t385 * t427 - t361;
	t397 = -pkin(10) - pkin(9);
	t405 = t387 * r_i_i_C(1) + t389 * r_i_i_C(2) - t397;
	t403 = t385 * t390 + t410;
	t351 = t363 * t383 + t370 * t429 - t384 * t412 - t385 * t415;
	t372 = -t416 + t421;
	t349 = t372 * t429 + t407 * t383 - t384 * t411;
	t350 = -t372 * t430 + t383 * t411 + t407 * t384;
	t402 = -(-t372 * t384 - t383 * t427) * qJD(5) + t431 * t350 - t406 * t349;
	t401 = -(-t370 * t384 + t415) * qJD(5) + t431 * (-t370 * t430 + t383 * t412 - t409) - t406 * t351;
	t358 = t403 * t383 + t385 * t417;
	t400 = -(-t390 * t383 - t417) * qJD(5) + t431 * (t403 * t384 - t428 * t430) - t406 * t358;
	t391 = sin(qJ(3));
	t398 = -qJD(3) * t391 * pkin(3) + t383 * qJD(5) + (-t406 * t383 + t431 * t384) * t385;
	t362 = t404 * qJD(1) + t370 * qJD(2);
	t360 = -qJD(1) * t413 - t396 * t418 + t408 * t424;
	t354 = (t370 * t385 - t412) * t383 + t409;
	t1 = [(t354 * t389 - t362 * t387) * r_i_i_C(1) + (-t354 * t387 - t362 * t389) * r_i_i_C(2) + t354 * pkin(4) - (t370 * t383 + t414) * qJD(5) - t363 * t382 + t362 * t397 - t431 * t351 + (-t396 * pkin(1) - pkin(8) * t427) * qJD(1) + (-t391 * t412 + (t370 * t391 + t394 * t425) * qJD(3)) * pkin(3), t434 * t360 - t405 * t361 - t398 * t404, (t394 * t411 + t361 * t391 + (-t372 * t394 - t391 * t427) * qJD(3)) * pkin(3) + t402, t402, t349, 0; (t350 * t389 - t360 * t387) * r_i_i_C(1) + (-t350 * t387 - t360 * t389) * r_i_i_C(2) + t350 * pkin(4) - (-t372 * t383 + t384 * t427) * qJD(5) - t361 * t382 + t360 * t397 + t431 * t349 + (-t393 * pkin(1) + pkin(8) * t425) * qJD(1) + (t391 * t411 + (-t372 * t391 + t393 * t426) * qJD(3)) * pkin(3), t405 * t363 - t434 * t362 + t398 * (t413 - t424), (t394 * t412 - t363 * t391 + (-t370 * t394 + t391 * t425) * qJD(3)) * pkin(3) + t401, t401, t351, 0; 0, (-t434 * t419 + (t405 * qJD(2) + t398) * t395) * t388, (-t391 * t410 + (-t390 * t391 - t392 * t426) * qJD(3)) * pkin(3) + t400, t400, t358, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:44:13
	% EndTime: 2019-10-10 12:44:14
	% DurationCPUTime: 1.14s
	% Computational Cost: add. (1171->121), mult. (1953->194), div. (0->0), fcn. (1954->14), ass. (0->79)
	t443 = pkin(12) + qJ(6);
	t439 = sin(t443);
	t440 = cos(t443);
	t470 = r_i_i_C(1) * t439 + r_i_i_C(2) * t440;
	t510 = t470 * qJD(6);
	t501 = r_i_i_C(3) + pkin(11) + qJ(5);
	t437 = cos(pkin(12)) * pkin(5) + pkin(4);
	t471 = r_i_i_C(1) * t440 - r_i_i_C(2) * t439;
	t468 = t437 + t471;
	t448 = cos(pkin(6));
	t452 = sin(qJ(1));
	t454 = cos(qJ(2));
	t490 = t452 * t454;
	t451 = sin(qJ(2));
	t455 = cos(qJ(1));
	t491 = t451 * t455;
	t424 = t448 * t491 + t490;
	t445 = qJ(3) + qJ(4);
	t441 = sin(t445);
	t447 = sin(pkin(6));
	t493 = t447 * t455;
	t431 = t441 * t493;
	t442 = cos(t445);
	t417 = t424 * t442 - t431;
	t489 = t454 * t455;
	t479 = t448 * t489;
	t492 = t451 * t452;
	t423 = -t479 + t492;
	t509 = t417 * t439 - t423 * t440;
	t508 = t417 * t440 + t423 * t439;
	t444 = qJD(3) + qJD(4);
	t450 = sin(qJ(3));
	t457 = (t468 * t441 - t501 * t442) * t444 + qJD(3) * t450 * pkin(3) + t442 * t510 - t441 * qJD(5);
	t472 = qJD(2) * t448 + qJD(1);
	t480 = t448 * t492;
	t487 = qJD(2) * t451;
	t415 = -qJD(1) * t480 - t452 * t487 + t472 * t489;
	t488 = qJD(1) * t447;
	t478 = t452 * t488;
	t482 = t442 * t493;
	t403 = (-t424 * t444 + t478) * t441 + t415 * t442 - t444 * t482;
	t453 = cos(qJ(3));
	t438 = pkin(3) * t453 + pkin(2);
	t504 = t501 * t441 + t468 * t442 + t438;
	t498 = t442 * t444;
	t497 = t447 * t451;
	t496 = t447 * t452;
	t495 = t447 * t453;
	t494 = t447 * t454;
	t486 = qJD(2) * t454;
	t481 = t444 * t497;
	t477 = t455 * t488;
	t476 = t447 * t487;
	t475 = t447 * t486;
	t474 = -pkin(5) * sin(pkin(12)) - pkin(10) - pkin(9);
	t425 = t448 * t490 + t491;
	t413 = t424 * qJD(1) + t425 * qJD(2);
	t469 = t444 * t496 - t413;
	t466 = t424 * t441 + t482;
	t465 = t471 * qJD(6);
	t463 = t444 * t448 + t475;
	t462 = t470 - t474;
	t402 = t415 * t441 + t424 * t498 - t444 * t431 - t442 * t478;
	t460 = qJD(5) * t417 - t402 * t468 + t501 * t403 + t510 * t466;
	t426 = -t480 + t489;
	t400 = t426 * t498 + t469 * t441 - t442 * t477;
	t401 = t469 * t442 + (-t426 * t444 + t477) * t441;
	t419 = -t426 * t441 + t442 * t496;
	t420 = t426 * t442 + t441 * t496;
	t459 = qJD(5) * t420 - t400 * t468 + t501 * t401 - t510 * t419;
	t410 = t463 * t441 + t442 * t481;
	t411 = -t441 * t481 + t463 * t442;
	t422 = t441 * t448 + t442 * t497;
	t458 = qJD(5) * t422 - t510 * (-t441 * t497 + t442 * t448) + t501 * t411 - t468 * t410;
	t414 = t425 * qJD(1) + t424 * qJD(2);
	t412 = -qJD(1) * t479 - t455 * t486 + t472 * t492;
	t391 = t401 * t440 - t412 * t439 + (-t420 * t439 + t425 * t440) * qJD(6);
	t390 = -t401 * t439 - t412 * t440 + (-t420 * t440 - t425 * t439) * qJD(6);
	t1 = [-t466 * qJD(5) - t415 * t438 - t468 * t403 - t462 * t414 - t501 * t402 + (t509 * r_i_i_C(1) + t508 * r_i_i_C(2)) * qJD(6) + (-pkin(1) * t455 - pkin(8) * t496) * qJD(1) + (-t450 * t478 + (t424 * t450 + t453 * t493) * qJD(3)) * pkin(3), t504 * t412 - t462 * t413 + t457 * t425 + t426 * t465, (t453 * t477 + t413 * t450 + (-t426 * t453 - t450 * t496) * qJD(3)) * pkin(3) + t459, t459, t400, r_i_i_C(1) * t390 - r_i_i_C(2) * t391; t391 * r_i_i_C(1) + t390 * r_i_i_C(2) - t419 * qJD(5) + t401 * t437 - t413 * t438 + t474 * t412 + t501 * t400 + (-pkin(1) * t452 + pkin(8) * t493) * qJD(1) + (t450 * t477 + (-t426 * t450 + t452 * t495) * qJD(3)) * pkin(3), -t414 * t504 + t462 * t415 + t457 * t423 + t424 * t465, (t453 * t478 - t415 * t450 + (-t424 * t453 + t450 * t493) * qJD(3)) * pkin(3) + t460, t460, t402, (-t403 * t439 + t414 * t440) * r_i_i_C(1) + (-t403 * t440 - t414 * t439) * r_i_i_C(2) + (-t508 * r_i_i_C(1) + t509 * r_i_i_C(2)) * qJD(6); 0, ((-qJD(2) * t504 + t465) * t451 + (t462 * qJD(2) - t457) * t454) * t447, (-t450 * t475 + (-t448 * t450 - t451 * t495) * qJD(3)) * pkin(3) + t458, t458, t410, (-t411 * t439 + t440 * t476) * r_i_i_C(1) + (-t411 * t440 - t439 * t476) * r_i_i_C(2) + ((-t422 * t440 + t439 * t494) * r_i_i_C(1) + (t422 * t439 + t440 * t494) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end