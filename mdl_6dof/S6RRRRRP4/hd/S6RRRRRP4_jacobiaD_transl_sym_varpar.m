% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:02
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:02:07
	% EndTime: 2019-10-10 13:02:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:02:07
	% EndTime: 2019-10-10 13:02:07
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
	% StartTime: 2019-10-10 13:02:07
	% EndTime: 2019-10-10 13:02:07
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (19->15), mult. (64->29), div. (0->0), fcn. (40->4), ass. (0->13)
	t28 = pkin(7) + r_i_i_C(3);
	t18 = sin(qJ(1));
	t27 = qJD(1) * t18;
	t20 = cos(qJ(1));
	t26 = qJD(1) * t20;
	t25 = qJD(2) * t18;
	t24 = qJD(2) * t20;
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t23 = r_i_i_C(1) * t17 + r_i_i_C(2) * t19;
	t22 = -r_i_i_C(1) * t19 + r_i_i_C(2) * t17 - pkin(1);
	t21 = t23 * qJD(2);
	t1 = [t23 * t25 + (-t28 * t18 + t22 * t20) * qJD(1), (t17 * t24 + t19 * t27) * r_i_i_C(2) + (t17 * t27 - t19 * t24) * r_i_i_C(1), 0, 0, 0, 0; -t20 * t21 + (t22 * t18 + t28 * t20) * qJD(1), (t17 * t25 - t19 * t26) * r_i_i_C(2) + (-t17 * t26 - t19 * t25) * r_i_i_C(1), 0, 0, 0, 0; 0, -t21, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:02:07
	% EndTime: 2019-10-10 13:02:07
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (81->26), mult. (114->37), div. (0->0), fcn. (73->6), ass. (0->27)
	t38 = qJD(2) + qJD(3);
	t39 = qJ(2) + qJ(3);
	t37 = cos(t39);
	t59 = r_i_i_C(2) * t37;
	t36 = sin(t39);
	t61 = r_i_i_C(1) * t36;
	t49 = t59 + t61;
	t47 = t49 * t38;
	t40 = sin(qJ(2));
	t62 = pkin(2) * t40;
	t63 = qJD(2) * t62 + t47;
	t60 = r_i_i_C(2) * t36;
	t58 = r_i_i_C(3) + pkin(8) + pkin(7);
	t57 = t37 * t38;
	t41 = sin(qJ(1));
	t56 = qJD(1) * t41;
	t43 = cos(qJ(1));
	t55 = qJD(1) * t43;
	t42 = cos(qJ(2));
	t54 = qJD(2) * t42;
	t53 = r_i_i_C(1) * t57;
	t52 = t38 * t60;
	t51 = qJD(1) * t59;
	t48 = -t42 * pkin(2) - r_i_i_C(1) * t37 - pkin(1) + t60;
	t46 = t41 * t51 + t56 * t61 + (t52 - t53) * t43;
	t31 = t41 * t52;
	t1 = [t63 * t41 + (-t58 * t41 + t48 * t43) * qJD(1), (t40 * t56 - t43 * t54) * pkin(2) + t46, t46, 0, 0, 0; -t63 * t43 + (t48 * t41 + t58 * t43) * qJD(1), t31 + (-pkin(2) * t54 - t53) * t41 + (-t49 - t62) * t55, -t43 * t51 + t31 + (-t36 * t55 - t41 * t57) * r_i_i_C(1), 0, 0, 0; 0, -t63, -t47, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:02:08
	% EndTime: 2019-10-10 13:02:09
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (270->58), mult. (394->95), div. (0->0), fcn. (297->8), ass. (0->54)
	t258 = qJ(2) + qJ(3);
	t255 = sin(t258);
	t262 = cos(qJ(4));
	t310 = r_i_i_C(1) * t262 + pkin(3);
	t313 = t255 * t310;
	t259 = sin(qJ(4));
	t293 = qJD(4) * t262;
	t256 = cos(t258);
	t257 = qJD(2) + qJD(3);
	t301 = t256 * t257;
	t312 = t255 * t293 + t259 * t301;
	t307 = pkin(9) + r_i_i_C(3);
	t288 = t307 * t256;
	t260 = sin(qJ(2));
	t302 = pkin(2) * qJD(2);
	t290 = t260 * t302;
	t305 = pkin(3) * t255;
	t311 = -t290 + (t288 - t305) * t257;
	t294 = qJD(4) * t259;
	t281 = t255 * t294;
	t308 = r_i_i_C(1) * t281 + t312 * r_i_i_C(2);
	t306 = pkin(2) * t260;
	t303 = r_i_i_C(2) * t259;
	t261 = sin(qJ(1));
	t300 = t257 * t261;
	t299 = t257 * t262;
	t264 = cos(qJ(1));
	t298 = t257 * t264;
	t297 = t262 * t264;
	t296 = qJD(1) * t261;
	t295 = qJD(1) * t264;
	t292 = t255 * t303;
	t291 = qJD(1) * t303;
	t289 = t307 * t255;
	t287 = t307 * t261;
	t286 = t255 * t299;
	t276 = qJD(4) * t256 - qJD(1);
	t275 = qJD(1) * t256 - qJD(4);
	t274 = t310 * t256;
	t273 = t310 * t264;
	t272 = t308 * t264 + t296 * t313;
	t271 = t276 * t259;
	t270 = t264 * t255 * t291 + t308 * t261 + t295 * t288;
	t263 = cos(qJ(2));
	t269 = -t263 * pkin(2) - pkin(3) * t256 - pkin(1) - t289;
	t268 = t255 * t298 + t275 * t261;
	t267 = -t263 * t302 + (-t274 - t289) * t257;
	t266 = -t256 * r_i_i_C(2) * t293 + (-t256 * t294 - t286) * r_i_i_C(1) + t307 * t301 + (-t305 + t292) * t257;
	t265 = -pkin(8) - pkin(7);
	t239 = -t275 * t297 + (t271 + t286) * t261;
	t238 = t276 * t262 * t261 + (-t255 * t300 + t275 * t264) * t259;
	t237 = t268 * t262 + t264 * t271;
	t236 = t268 * t259 - t276 * t297;
	t1 = [t239 * r_i_i_C(1) + t238 * r_i_i_C(2) - t311 * t261 + (t261 * t265 + t269 * t264) * qJD(1), (-t288 - t292 + t306) * t296 + t267 * t264 + t272, (-t261 * t291 - t307 * t298) * t255 + (-qJD(1) * t287 - t257 * t273) * t256 + t272, t236 * r_i_i_C(1) + t237 * r_i_i_C(2), 0, 0; -t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + t311 * t264 + (t269 * t261 - t264 * t265) * qJD(1), (-t306 - t313) * t295 + t267 * t261 + t270, -t274 * t300 + (-qJD(1) * t273 - t257 * t287) * t255 + t270, -t238 * r_i_i_C(1) + t239 * r_i_i_C(2), 0, 0; 0, t266 - t290, t266, (-t256 * t299 + t281) * r_i_i_C(2) - t312 * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:02:08
	% EndTime: 2019-10-10 13:02:09
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (524->77), mult. (560->110), div. (0->0), fcn. (431->10), ass. (0->69)
	t291 = qJ(4) + qJ(5);
	t285 = sin(t291);
	t292 = qJ(2) + qJ(3);
	t286 = sin(t292);
	t288 = cos(t292);
	t290 = qJD(2) + qJD(3);
	t343 = t288 * t290;
	t287 = cos(t291);
	t289 = qJD(4) + qJD(5);
	t344 = t287 * t289;
	t364 = t285 * t343 + t286 * t344;
	t299 = -pkin(10) - pkin(9);
	t293 = sin(qJ(4));
	t347 = pkin(4) * qJD(4);
	t333 = t293 * t347;
	t363 = t290 * t299 + t333;
	t296 = cos(qJ(4));
	t283 = t296 * pkin(4) + pkin(3);
	t361 = r_i_i_C(1) * t287 + t283;
	t362 = t286 * t361 + t288 * t299;
	t294 = sin(qJ(2));
	t348 = pkin(2) * qJD(2);
	t330 = t294 * t348;
	t345 = t286 * t290;
	t349 = r_i_i_C(3) - t299;
	t353 = pkin(4) * t293;
	t360 = (t349 * t290 - t333) * t288 + (pkin(8) + pkin(7) + t353) * qJD(1) - t283 * t345 - t330;
	t346 = t285 * t286;
	t329 = t289 * t346;
	t359 = r_i_i_C(1) * t329 + t364 * r_i_i_C(2) + t363 * t286;
	t298 = cos(qJ(1));
	t317 = t288 * t289 - qJD(1);
	t358 = t298 * t317;
	t337 = qJD(1) * t288;
	t316 = -t289 + t337;
	t295 = sin(qJ(1));
	t326 = t295 * t345;
	t356 = t316 * t298 - t326;
	t354 = pkin(2) * t294;
	t351 = r_i_i_C(2) * t287;
	t350 = r_i_i_C(3) * t288;
	t341 = t290 * t298;
	t325 = t286 * t341;
	t304 = t316 * t295 + t325;
	t261 = t304 * t285 - t287 * t358;
	t262 = t285 * t358 + t304 * t287;
	t339 = t261 * r_i_i_C(1) + t262 * r_i_i_C(2);
	t311 = t317 * t295;
	t263 = t356 * t285 + t287 * t311;
	t264 = t285 * t311 - t356 * t287;
	t338 = -t263 * r_i_i_C(1) + t264 * r_i_i_C(2);
	t336 = qJD(1) * t295;
	t335 = qJD(1) * t298;
	t334 = r_i_i_C(2) * t346;
	t332 = t296 * t347;
	t331 = r_i_i_C(2) * qJD(1) * t285;
	t314 = -qJD(4) + t337;
	t313 = -r_i_i_C(1) * t285 - t351;
	t312 = t361 * t290;
	t310 = (-qJD(4) * t288 + qJD(1)) * t296;
	t309 = t298 * t286 * t331 + t359 * t295 + t335 * t350;
	t306 = t359 * t298 + t362 * t336;
	t305 = (-r_i_i_C(3) * t286 - t288 * t361) * t290;
	t297 = cos(qJ(2));
	t303 = t332 + (-t297 * pkin(2) - t283 * t288 - t349 * t286 - pkin(1)) * qJD(1);
	t302 = -t297 * t348 + t305;
	t301 = t290 * t334 + r_i_i_C(3) * t343 - t286 * t312 + (t313 * t289 - t363) * t288;
	t274 = r_i_i_C(2) * t329;
	t1 = [t264 * r_i_i_C(1) + t263 * r_i_i_C(2) - t360 * t295 + t303 * t298, (-t334 - t350 + t354) * t336 + t302 * t298 + t306, (-r_i_i_C(3) * t341 - t295 * t331) * t286 + (-r_i_i_C(3) * t336 - t298 * t312) * t288 + t306, (t298 * t310 + (t314 * t295 + t325) * t293) * pkin(4) + t339, t339, 0; -t262 * r_i_i_C(1) + t261 * r_i_i_C(2) + t303 * t295 + t360 * t298, t302 * t295 + (-t362 - t354) * t335 + t309, t295 * t305 - t335 * t362 + t309, (t295 * t310 + (-t314 * t298 + t326) * t293) * pkin(4) + t338, t338, 0; 0, t301 - t330, t301, t274 + (-r_i_i_C(1) * t344 - t332) * t286 + (t313 - t353) * t343, -t364 * r_i_i_C(1) - t343 * t351 + t274, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:02:09
	% EndTime: 2019-10-10 13:02:10
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (959->97), mult. (946->129), div. (0->0), fcn. (771->10), ass. (0->76)
	t356 = qJ(4) + qJ(5);
	t350 = sin(t356);
	t433 = r_i_i_C(3) + qJ(6);
	t436 = t350 * t433;
	t352 = cos(t356);
	t354 = qJD(4) + qJD(5);
	t363 = cos(qJ(1));
	t412 = t363 * t350;
	t397 = t354 * t412;
	t360 = sin(qJ(1));
	t410 = qJD(1) * t360;
	t435 = t352 * t410 + t397;
	t434 = pkin(5) + r_i_i_C(1);
	t357 = qJ(2) + qJ(3);
	t353 = cos(t357);
	t391 = t353 * t410;
	t351 = sin(t357);
	t355 = qJD(2) + qJD(3);
	t418 = t351 * t355;
	t395 = t363 * t418;
	t431 = t391 + t395;
	t358 = sin(qJ(4));
	t419 = pkin(4) * qJD(4);
	t406 = t358 * t419;
	t408 = qJD(6) * t350;
	t430 = -t406 + t408;
	t361 = cos(qJ(4));
	t348 = t361 * pkin(4) + pkin(3);
	t359 = sin(qJ(2));
	t420 = pkin(2) * qJD(2);
	t405 = t359 * t420;
	t364 = -pkin(10) - pkin(9);
	t422 = r_i_i_C(2) - t364;
	t424 = pkin(4) * t358;
	t429 = (t422 * t355 + t430) * t353 + (pkin(8) + pkin(7) + t424) * qJD(1) - t348 * t418 - t405;
	t425 = pkin(2) * t359;
	t423 = r_i_i_C(2) * t353;
	t417 = t352 * t354;
	t416 = t353 * t355;
	t415 = t355 * t360;
	t414 = t360 * t350;
	t413 = t360 * t352;
	t411 = t363 * t352;
	t409 = qJD(1) * t363;
	t407 = t352 * qJD(6);
	t404 = t361 * t419;
	t403 = t434 * t350;
	t402 = t351 * t417;
	t399 = t354 * t414;
	t398 = t351 * t415;
	t396 = t354 * t411;
	t394 = t433 * t352 * t416 + t351 * t407;
	t390 = t433 * t354;
	t388 = t351 * t406;
	t387 = qJD(1) * t353 - qJD(4);
	t382 = t434 * t351 * t399 + t360 * t388 + t364 * t398 + t409 * t423;
	t381 = (-qJD(4) * t353 + qJD(1)) * t361;
	t379 = t353 * t411 + t414;
	t377 = t350 * t409 + t354 * t413;
	t376 = -t352 * t434 - t436;
	t375 = -t352 * t390 - t408;
	t310 = t350 * t395 - t353 * t396 - t399 + (t353 * t414 + t411) * qJD(1);
	t311 = t431 * t352 + t353 * t397 - t377;
	t374 = qJD(6) * t379 + t434 * t310 - t433 * t311;
	t312 = -t350 * t398 + t377 * t353 - t435;
	t313 = t379 * qJD(1) - t352 * t398 - t353 * t399 - t396;
	t373 = -(-t353 * t413 + t412) * qJD(6) + t433 * t313 - t434 * t312;
	t372 = -t348 + t376;
	t371 = t363 * t388 + t431 * t364 + (t434 * t435 + (t348 + t436) * t410) * t351;
	t370 = t372 * t351 - t353 * t364;
	t362 = cos(qJ(2));
	t369 = t404 - t407 + (-t362 * pkin(2) - t348 * t353 - t422 * t351 - pkin(1)) * qJD(1);
	t368 = t375 * t351 + (-r_i_i_C(2) * t351 + t372 * t353) * t355;
	t367 = -t362 * t420 + t368;
	t366 = r_i_i_C(2) * t416 + t370 * t355 + (-t354 * t403 + t433 * t417 + t430) * t353;
	t1 = [-t312 * t433 - t313 * t434 - t429 * t360 + t369 * t363, (-t423 + t425) * t410 + t367 * t363 + t371, -r_i_i_C(2) * t391 + t368 * t363 + t371, (t363 * t381 + (t387 * t360 + t395) * t358) * pkin(4) + t374, t374, -t310; -t310 * t433 - t311 * t434 + t369 * t360 + t429 * t363, (t370 - t425) * t409 + t367 * t360 + t382, (-t364 * t409 + t372 * t415) * t353 + ((-r_i_i_C(2) * t355 + t375) * t360 + t372 * t409) * t351 + t382, (t360 * t381 + (-t387 * t363 + t398) * t358) * pkin(4) + t373, t373, t312; 0, t366 - t405, t366, (-t403 - t424) * t416 + (t376 * t354 - t404) * t351 + t394, -t434 * t402 + (-t351 * t390 - t416 * t434) * t350 + t394, t350 * t416 + t402;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end