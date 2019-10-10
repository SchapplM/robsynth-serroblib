% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 12:38:37
	% EndTime: 2019-10-10 12:38:37
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
	% StartTime: 2019-10-10 12:38:39
	% EndTime: 2019-10-10 12:38:39
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
	% StartTime: 2019-10-10 12:38:39
	% EndTime: 2019-10-10 12:38:40
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (474->71), mult. (702->99), div. (0->0), fcn. (565->8), ass. (0->56)
	t316 = sin(qJ(4));
	t319 = cos(qJ(4));
	t321 = cos(qJ(1));
	t357 = qJD(4) * t321;
	t318 = sin(qJ(1));
	t360 = qJD(1) * t318;
	t330 = t316 * t357 + t319 * t360;
	t369 = -r_i_i_C(1) - pkin(4);
	t367 = r_i_i_C(3) + qJ(5);
	t374 = t367 * t316;
	t315 = qJ(2) + qJ(3);
	t313 = cos(t315);
	t370 = pkin(9) + r_i_i_C(2);
	t352 = t370 * t313;
	t314 = qJD(2) + qJD(3);
	t351 = t370 * t314;
	t317 = sin(qJ(2));
	t366 = pkin(2) * qJD(2);
	t354 = t317 * t366;
	t356 = qJD(5) * t316;
	t312 = sin(t315);
	t365 = t312 * t314;
	t373 = -pkin(3) * t365 + (t351 + t356) * t313 - t354;
	t358 = qJD(4) * t319;
	t328 = -t367 * t358 - t356;
	t368 = pkin(2) * t317;
	t364 = t313 * t314;
	t363 = t313 * t318;
	t362 = t316 * t318;
	t361 = t319 * t321;
	t359 = qJD(1) * t321;
	t355 = t319 * qJD(5);
	t353 = t370 * t312;
	t350 = t369 * t316;
	t349 = t318 * t365;
	t348 = t321 * t365;
	t343 = qJD(4) * t362;
	t341 = t319 * t357;
	t340 = -t369 * t312 * t343 + t359 * t352;
	t335 = t313 * t361 + t362;
	t320 = cos(qJ(2));
	t333 = -pkin(2) * t320 - pkin(3) * t313 - pkin(1) - t353;
	t332 = (-t369 * t330 + (pkin(3) + t374) * t360) * t312;
	t331 = t369 * t319 - t374;
	t329 = t316 * t359 + t318 * t358;
	t327 = -pkin(3) + t331;
	t326 = t327 * t314;
	t325 = t312 * t326 + t370 * t364 + (qJD(4) * t350 - t328) * t313;
	t324 = t328 * t312 + (t327 * t313 - t353) * t314;
	t323 = -t320 * t366 + t324;
	t322 = -pkin(8) - pkin(7);
	t284 = t335 * qJD(1) - t313 * t343 - t319 * t349 - t341;
	t283 = t329 * t313 - t316 * t349 - t330;
	t282 = t330 * t313 + t319 * t348 - t329;
	t281 = t316 * t348 - t313 * t341 - t343 + (t313 * t362 + t361) * qJD(1);
	t1 = [-t321 * t355 + t369 * t284 - t367 * t283 - t373 * t318 + (t318 * t322 + t333 * t321) * qJD(1), (-t352 + t368) * t360 + t323 * t321 + t332, t324 * t321 - t352 * t360 + t332, t335 * qJD(5) - t369 * t281 - t367 * t282, -t281, 0; -t318 * t355 + t369 * t282 - t367 * t281 + t373 * t321 + (t333 * t318 - t321 * t322) * qJD(1), (t327 * t312 - t368) * t359 + t323 * t318 + t340, t326 * t363 + ((-t351 + t328) * t318 + t327 * t359) * t312 + t340, -(t316 * t321 - t319 * t363) * qJD(5) + t367 * t284 + t369 * t283, t283, 0; 0, t325 - t354, t325, (t367 * t319 + t350) * t364 + (t331 * qJD(4) + t355) * t312, t312 * t358 + t316 * t364, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:38:39
	% EndTime: 2019-10-10 12:38:40
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (972->126), mult. (1504->185), div. (0->0), fcn. (1371->10), ass. (0->81)
	t368 = sin(qJ(4));
	t372 = cos(qJ(4));
	t374 = cos(qJ(1));
	t422 = qJD(4) * t374;
	t370 = sin(qJ(1));
	t425 = qJD(1) * t370;
	t388 = t368 * t422 + t372 * t425;
	t448 = pkin(4) + pkin(5);
	t447 = pkin(10) + r_i_i_C(3);
	t367 = sin(qJ(6));
	t371 = cos(qJ(6));
	t397 = t367 * t368 + t371 * t372;
	t441 = qJD(4) - qJD(6);
	t379 = t441 * t397;
	t398 = t367 * t372 - t368 * t371;
	t446 = t441 * t398;
	t423 = qJD(4) * t372;
	t393 = -qJ(5) * t423 - qJD(5) * t368;
	t366 = qJ(2) + qJ(3);
	t364 = cos(t366);
	t426 = t374 * t372;
	t429 = t370 * t368;
	t328 = t364 * t429 + t426;
	t427 = t374 * t368;
	t428 = t370 * t372;
	t329 = t364 * t428 - t427;
	t401 = t328 * t367 + t329 * t371;
	t402 = t328 * t371 - t329 * t367;
	t443 = (t401 * r_i_i_C(1) + t402 * r_i_i_C(2)) * qJD(6);
	t363 = sin(t366);
	t365 = qJD(2) + qJD(3);
	t369 = sin(qJ(2));
	t433 = pkin(2) * qJD(2);
	t419 = t369 * t433;
	t420 = pkin(9) - t447;
	t442 = -t419 + (-pkin(3) * t363 + t420 * t364) * t365;
	t436 = pkin(2) * t369;
	t434 = pkin(9) * t364;
	t432 = qJ(5) * t368;
	t424 = qJD(1) * t374;
	t387 = t368 * t424 + t370 * t423;
	t430 = t365 * t370;
	t417 = t363 * t430;
	t326 = t387 * t364 - t368 * t417 - t388;
	t323 = t326 * t371;
	t431 = t364 * t365;
	t418 = t448 * t368;
	t416 = t374 * t365 * t363;
	t413 = t364 * t425;
	t411 = qJD(4) * t429;
	t409 = t372 * t422;
	t408 = -t367 * r_i_i_C(1) - qJ(5);
	t407 = t367 * r_i_i_C(2) - t448;
	t392 = t398 * t365;
	t384 = t364 * t392;
	t403 = -(-t379 * t363 + t384) * r_i_i_C(1) - (t363 * t446 + t397 * t431) * r_i_i_C(2);
	t330 = t364 * t427 - t428;
	t331 = t364 * t426 + t429;
	t400 = t330 * t371 - t331 * t367;
	t399 = t330 * t367 + t331 * t371;
	t396 = t371 * r_i_i_C(2) - t408;
	t395 = t371 * r_i_i_C(1) - t407;
	t394 = -t372 * t448 - t432;
	t391 = t397 * t365;
	t390 = -pkin(3) + t394;
	t385 = t364 * t391;
	t389 = (t370 * t384 + (-t379 * t370 + t398 * t424) * t363) * r_i_i_C(2) + (-t370 * t385 + (-t370 * t446 - t397 * t424) * t363) * r_i_i_C(1) + t424 * t434 + t447 * t417 + t448 * t363 * t411;
	t373 = cos(qJ(2));
	t386 = -t373 * pkin(2) - pkin(3) * t364 - t420 * t363 - pkin(1);
	t383 = t447 * (t416 + t413) + (-t385 * r_i_i_C(1) + t384 * r_i_i_C(2)) * t374 + ((-t379 * t374 - t398 * t425) * r_i_i_C(2) + (-t374 * t446 + t397 * t425) * r_i_i_C(1) + (t432 + pkin(3)) * t425 + t448 * t388) * t363;
	t382 = t390 * t363 - t364 * t447;
	t378 = t393 * t363 + (-pkin(9) * t363 + t390 * t364) * t365;
	t377 = -t373 * t433 + t378;
	t376 = pkin(9) * t431 + t382 * t365 + (-t391 * r_i_i_C(1) + t392 * r_i_i_C(2)) * t363 + (r_i_i_C(1) * t446 + r_i_i_C(2) * t379 - qJD(4) * t418 - t393) * t364;
	t375 = -pkin(8) - pkin(7);
	t327 = t331 * qJD(1) - t364 * t411 - t372 * t417 - t409;
	t325 = t388 * t364 + t372 * t416 - t387;
	t324 = t328 * qJD(1) - t364 * t409 + t368 * t416 - t411;
	t314 = t400 * qJD(6) - t324 * t367 - t325 * t371;
	t313 = -t399 * qJD(6) - t324 * t371 + t325 * t367;
	t1 = [-t323 * r_i_i_C(2) - t328 * qJD(5) + t408 * t326 - t395 * t327 + (-t402 * r_i_i_C(1) + t401 * r_i_i_C(2)) * qJD(6) - t442 * t370 + (t370 * t375 + t386 * t374) * qJD(1), t383 + (-t434 + t436) * t425 + t377 * t374, -pkin(9) * t413 + t378 * t374 + t383, t331 * qJD(5) - t396 * t325 + t395 * t324 + (t399 * r_i_i_C(1) + t400 * r_i_i_C(2)) * qJD(6), -t324, t313 * r_i_i_C(1) - t314 * r_i_i_C(2); t314 * r_i_i_C(1) + t313 * r_i_i_C(2) - t324 * qJ(5) + t330 * qJD(5) - t448 * t325 + t442 * t374 + (t386 * t370 - t374 * t375) * qJD(1), (t382 - t436) * t424 + t377 * t370 + t389, (t390 * t430 - t424 * t447) * t364 + ((-pkin(9) * t365 + t393) * t370 + t390 * t424) * t363 + t389, -t323 * r_i_i_C(1) + t329 * qJD(5) + t407 * t326 + t396 * t327 + t443, t326, (-t327 * t367 + t323) * r_i_i_C(1) + (-t326 * t367 - t327 * t371) * r_i_i_C(2) - t443; 0, t376 - t419, t376, (qJ(5) * t372 - t418) * t431 + (t394 * qJD(4) + qJD(5) * t372) * t363 - t403, t363 * t423 + t368 * t431, t403;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end