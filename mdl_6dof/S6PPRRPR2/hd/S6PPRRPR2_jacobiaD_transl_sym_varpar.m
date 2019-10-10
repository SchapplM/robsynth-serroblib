% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (23->19), mult. (86->45), div. (0->0), fcn. (88->10), ass. (0->22)
	t111 = sin(pkin(11));
	t117 = cos(pkin(6));
	t126 = t111 * t117;
	t112 = sin(pkin(7));
	t113 = sin(pkin(6));
	t125 = t112 * t113;
	t115 = cos(pkin(11));
	t124 = t115 * t117;
	t116 = cos(pkin(7));
	t118 = sin(qJ(3));
	t123 = t116 * t118;
	t119 = cos(qJ(3));
	t122 = t116 * t119;
	t121 = t118 * t125;
	t120 = t119 * t125;
	t114 = cos(pkin(12));
	t110 = sin(pkin(12));
	t109 = -t110 * t126 + t115 * t114;
	t108 = -t115 * t110 - t114 * t126;
	t107 = t110 * t124 + t111 * t114;
	t106 = -t111 * t110 + t114 * t124;
	t1 = [0, 0, ((-t108 * t123 - t109 * t119 - t111 * t121) * r_i_i_C(1) + (-t108 * t122 + t109 * t118 - t111 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-t106 * t123 - t107 * t119 + t115 * t121) * r_i_i_C(1) + (-t106 * t122 + t107 * t118 + t115 * t120) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, 0, ((-r_i_i_C(1) * t118 - r_i_i_C(2) * t119) * t117 * t112 + ((-t110 * t119 - t114 * t123) * r_i_i_C(1) + (t110 * t118 - t114 * t122) * r_i_i_C(2)) * t113) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:35
	% EndTime: 2019-10-09 21:12:35
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (153->45), mult. (509->95), div. (0->0), fcn. (560->12), ass. (0->45)
	t281 = sin(pkin(12));
	t284 = sin(pkin(6));
	t290 = sin(qJ(3));
	t292 = cos(qJ(3));
	t285 = cos(pkin(12));
	t287 = cos(pkin(7));
	t306 = t285 * t287;
	t283 = sin(pkin(7));
	t288 = cos(pkin(6));
	t309 = t283 * t288;
	t268 = (t281 * t292 + t290 * t306) * t284 + t290 * t309;
	t286 = cos(pkin(11));
	t282 = sin(pkin(11));
	t311 = t282 * t288;
	t277 = -t281 * t311 + t286 * t285;
	t276 = -t286 * t281 - t285 * t311;
	t310 = t283 * t284;
	t296 = t276 * t287 + t282 * t310;
	t264 = t277 * t292 + t296 * t290;
	t305 = t286 * t288;
	t275 = t281 * t305 + t282 * t285;
	t274 = -t282 * t281 + t285 * t305;
	t302 = t286 * t310;
	t297 = -t274 * t287 + t302;
	t316 = -t275 * t292 + t297 * t290;
	t315 = -pkin(9) - r_i_i_C(3);
	t314 = t275 * t290;
	t308 = t284 * t285;
	t307 = t284 * t287;
	t304 = qJD(3) * t290;
	t303 = qJD(3) * t292;
	t300 = t283 * t303;
	t299 = t287 * t303;
	t289 = sin(qJ(4));
	t291 = cos(qJ(4));
	t298 = t289 * r_i_i_C(1) + t291 * r_i_i_C(2);
	t294 = qJD(4) * t298;
	t293 = qJD(3) * (t291 * r_i_i_C(1) - t289 * r_i_i_C(2) + pkin(3));
	t273 = -t283 * t308 + t288 * t287;
	t270 = -t276 * t283 + t282 * t307;
	t269 = -t274 * t283 - t286 * t307;
	t265 = t284 * t281 * t304 - t288 * t300 - t299 * t308;
	t259 = -t282 * t284 * t300 - t276 * t299 + t277 * t304;
	t257 = -t274 * t299 + (t292 * t302 + t314) * qJD(3);
	t1 = [0, 0, t315 * t259 - (-t277 * t290 + t296 * t292) * t294 - t264 * t293, t298 * t259 + ((-t264 * t291 - t270 * t289) * r_i_i_C(1) + (t264 * t289 - t270 * t291) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t315 * t257 - (-t297 * t292 - t314) * t294 + t316 * t293, t298 * t257 + ((-t269 * t289 + t291 * t316) * r_i_i_C(1) + (-t269 * t291 - t289 * t316) * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, t315 * t265 - (t292 * t309 + (-t281 * t290 + t292 * t306) * t284) * t294 - t268 * t293, t298 * t265 + ((-t268 * t291 - t273 * t289) * r_i_i_C(1) + (t268 * t289 - t273 * t291) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:36
	% EndTime: 2019-10-09 21:12:36
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (338->54), mult. (1082->105), div. (0->0), fcn. (1228->12), ass. (0->52)
	t325 = sin(pkin(12));
	t328 = sin(pkin(6));
	t334 = sin(qJ(3));
	t336 = cos(qJ(3));
	t329 = cos(pkin(12));
	t331 = cos(pkin(7));
	t352 = t329 * t331;
	t327 = sin(pkin(7));
	t332 = cos(pkin(6));
	t355 = t327 * t332;
	t312 = (t325 * t336 + t334 * t352) * t328 + t334 * t355;
	t330 = cos(pkin(11));
	t326 = sin(pkin(11));
	t357 = t326 * t332;
	t321 = -t325 * t357 + t330 * t329;
	t320 = -t330 * t325 - t329 * t357;
	t356 = t327 * t328;
	t340 = t320 * t331 + t326 * t356;
	t308 = t321 * t336 + t340 * t334;
	t351 = t330 * t332;
	t319 = t325 * t351 + t326 * t329;
	t318 = -t326 * t325 + t329 * t351;
	t348 = t330 * t356;
	t341 = -t318 * t331 + t348;
	t364 = -t319 * t336 + t341 * t334;
	t363 = pkin(4) - r_i_i_C(2);
	t362 = -pkin(9) - r_i_i_C(1);
	t361 = r_i_i_C(3) + qJ(5);
	t360 = t319 * t334;
	t354 = t328 * t329;
	t353 = t328 * t331;
	t350 = qJD(3) * t334;
	t349 = qJD(3) * t336;
	t346 = t327 * t349;
	t345 = t331 * t349;
	t313 = -t318 * t327 - t330 * t353;
	t333 = sin(qJ(4));
	t335 = cos(qJ(4));
	t344 = t313 * t333 - t335 * t364;
	t314 = -t320 * t327 + t326 * t353;
	t343 = t308 * t335 + t314 * t333;
	t317 = -t327 * t354 + t332 * t331;
	t342 = t312 * t335 + t317 * t333;
	t338 = qJD(3) * (t361 * t333 + t363 * t335 + pkin(3));
	t337 = qJD(5) * t333 + (-t363 * t333 + t361 * t335) * qJD(4);
	t309 = t328 * t325 * t350 - t332 * t346 - t345 * t354;
	t303 = -t326 * t328 * t346 - t320 * t345 + t321 * t350;
	t301 = -t318 * t345 + (t336 * t348 + t360) * qJD(3);
	t299 = t342 * qJD(4) - t309 * t333;
	t297 = t343 * qJD(4) - t303 * t333;
	t295 = t344 * qJD(4) - t301 * t333;
	t1 = [0, 0, t362 * t303 + t337 * (-t321 * t334 + t340 * t336) - t308 * t338, t343 * qJD(5) + t361 * (-t303 * t335 + (-t308 * t333 + t314 * t335) * qJD(4)) - t363 * t297, t297, 0; 0, 0, t362 * t301 + t337 * (-t341 * t336 - t360) + t364 * t338, t344 * qJD(5) + t361 * (-t301 * t335 + (t313 * t335 + t333 * t364) * qJD(4)) - t363 * t295, t295, 0; 0, 0, t362 * t309 + t337 * (t336 * t355 + (-t325 * t334 + t336 * t352) * t328) - t312 * t338, t342 * qJD(5) + t361 * (-t309 * t335 + (-t312 * t333 + t317 * t335) * qJD(4)) - t363 * t299, t299, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:37
	% EndTime: 2019-10-09 21:12:38
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (737->77), mult. (2326->142), div. (0->0), fcn. (2700->14), ass. (0->68)
	t470 = sin(pkin(12));
	t471 = sin(pkin(11));
	t456 = t471 * t470;
	t474 = cos(pkin(12));
	t475 = cos(pkin(11));
	t463 = t475 * t474;
	t477 = cos(pkin(6));
	t438 = -t477 * t463 + t456;
	t476 = cos(pkin(7));
	t436 = t438 * t476;
	t472 = sin(pkin(7));
	t473 = sin(pkin(6));
	t461 = t473 * t472;
	t447 = t475 * t461;
	t482 = t436 + t447;
	t458 = t471 * t474;
	t462 = t475 * t470;
	t439 = t477 * t458 + t462;
	t457 = t471 * t473;
	t481 = t439 * t476 - t472 * t457;
	t464 = t476 * t473;
	t480 = t474 * t464 + t472 * t477;
	t419 = t477 * t462 + t458;
	t430 = sin(qJ(3));
	t478 = cos(qJ(3));
	t402 = t419 * t430 + t482 * t478;
	t479 = t481 * t478;
	t460 = t473 * t470;
	t410 = t430 * t460 - t480 * t478;
	t428 = sin(qJ(6));
	t431 = cos(qJ(6));
	t466 = t431 * r_i_i_C(1) - t428 * r_i_i_C(2);
	t444 = t466 * qJD(6) + qJD(5);
	t469 = qJD(3) * t430;
	t468 = pkin(4) + pkin(10) + r_i_i_C(3);
	t467 = t419 * t478;
	t465 = -t428 * r_i_i_C(1) - t431 * r_i_i_C(2);
	t403 = -t482 * t430 + t467;
	t412 = t438 * t472 - t475 * t464;
	t429 = sin(qJ(4));
	t432 = cos(qJ(4));
	t455 = t403 * t432 + t412 * t429;
	t394 = t403 * t429 - t412 * t432;
	t420 = -t477 * t456 + t463;
	t405 = t420 * t478 - t481 * t430;
	t413 = t439 * t472 + t476 * t457;
	t454 = t405 * t432 + t413 * t429;
	t396 = t405 * t429 - t413 * t432;
	t411 = t480 * t430 + t478 * t460;
	t418 = -t474 * t461 + t477 * t476;
	t453 = t411 * t432 + t418 * t429;
	t406 = t411 * t429 - t418 * t432;
	t452 = qJ(5) - t465;
	t450 = -pkin(5) - pkin(9) - t466;
	t449 = qJD(6) * t465;
	t440 = -t452 * t429 - t468 * t432 - pkin(3);
	t433 = -t444 * t429 + (t468 * t429 - t452 * t432) * qJD(4);
	t409 = t411 * qJD(3);
	t408 = t410 * qJD(3);
	t404 = t420 * t430 + t479;
	t401 = t405 * qJD(3);
	t400 = qJD(3) * t479 + t420 * t469;
	t399 = -t447 * t469 + (-t430 * t436 + t467) * qJD(3);
	t398 = t402 * qJD(3);
	t392 = t453 * qJD(4) - t408 * t429;
	t390 = t454 * qJD(4) - t400 * t429;
	t388 = t455 * qJD(4) - t398 * t429;
	t1 = [0, 0, t450 * t400 + t440 * t401 + t433 * t404 + t405 * t449, t444 * t454 + t452 * (-t396 * qJD(4) - t400 * t432) - t468 * t390, t390, (t390 * t431 - t401 * t428) * r_i_i_C(1) + (-t390 * t428 - t401 * t431) * r_i_i_C(2) + ((-t396 * t428 - t404 * t431) * r_i_i_C(1) + (-t396 * t431 + t404 * t428) * r_i_i_C(2)) * qJD(6); 0, 0, t450 * t398 + t440 * t399 + t433 * t402 + t403 * t449, t444 * t455 + t452 * (-t394 * qJD(4) - t398 * t432) - t468 * t388, t388, (t388 * t431 - t399 * t428) * r_i_i_C(1) + (-t388 * t428 - t399 * t431) * r_i_i_C(2) + ((-t394 * t428 - t402 * t431) * r_i_i_C(1) + (-t394 * t431 + t402 * t428) * r_i_i_C(2)) * qJD(6); 0, 0, t450 * t408 + t440 * t409 + t433 * t410 + t411 * t449, t444 * t453 + t452 * (-t406 * qJD(4) - t408 * t432) - t468 * t392, t392, (t392 * t431 - t409 * t428) * r_i_i_C(1) + (-t392 * t428 - t409 * t431) * r_i_i_C(2) + ((-t406 * t428 - t410 * t431) * r_i_i_C(1) + (-t406 * t431 + t410 * t428) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end