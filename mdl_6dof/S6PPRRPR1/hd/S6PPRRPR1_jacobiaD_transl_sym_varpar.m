% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiaD_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
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
	% StartTime: 2019-10-09 21:10:41
	% EndTime: 2019-10-09 21:10:41
	% DurationCPUTime: 0.26s
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
	% StartTime: 2019-10-09 21:10:43
	% EndTime: 2019-10-09 21:10:43
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (411->56), mult. (1327->109), div. (0->0), fcn. (1510->14), ass. (0->54)
	t361 = sin(pkin(12));
	t364 = sin(pkin(6));
	t371 = sin(qJ(3));
	t373 = cos(qJ(3));
	t366 = cos(pkin(12));
	t368 = cos(pkin(7));
	t391 = t366 * t368;
	t363 = sin(pkin(7));
	t369 = cos(pkin(6));
	t394 = t363 * t369;
	t347 = (t361 * t373 + t371 * t391) * t364 + t371 * t394;
	t367 = cos(pkin(11));
	t362 = sin(pkin(11));
	t396 = t362 * t369;
	t356 = -t361 * t396 + t366 * t367;
	t355 = -t361 * t367 - t366 * t396;
	t395 = t363 * t364;
	t377 = t355 * t368 + t362 * t395;
	t343 = t356 * t373 + t371 * t377;
	t390 = t367 * t369;
	t354 = t361 * t390 + t362 * t366;
	t353 = -t361 * t362 + t366 * t390;
	t387 = t367 * t395;
	t378 = -t353 * t368 + t387;
	t401 = -t354 * t373 + t371 * t378;
	t400 = r_i_i_C(3) + qJ(5);
	t399 = t354 * t371;
	t393 = t364 * t366;
	t392 = t364 * t368;
	t389 = qJD(3) * t371;
	t388 = qJD(3) * t373;
	t385 = t363 * t388;
	t384 = t368 * t388;
	t348 = -t353 * t363 - t367 * t392;
	t370 = sin(qJ(4));
	t372 = cos(qJ(4));
	t383 = t348 * t370 - t372 * t401;
	t349 = -t355 * t363 + t362 * t392;
	t382 = t343 * t372 + t349 * t370;
	t352 = -t363 * t393 + t368 * t369;
	t381 = t347 * t372 + t352 * t370;
	t360 = sin(pkin(13));
	t365 = cos(pkin(13));
	t380 = r_i_i_C(1) * t365 - r_i_i_C(2) * t360 + pkin(4);
	t379 = -r_i_i_C(1) * t360 - r_i_i_C(2) * t365 - pkin(9);
	t375 = qJD(3) * (t370 * t400 + t372 * t380 + pkin(3));
	t374 = t370 * qJD(5) + (-t370 * t380 + t372 * t400) * qJD(4);
	t344 = t361 * t364 * t389 - t369 * t385 - t384 * t393;
	t338 = -t362 * t364 * t385 - t355 * t384 + t356 * t389;
	t336 = -t353 * t384 + (t373 * t387 + t399) * qJD(3);
	t334 = qJD(4) * t381 - t344 * t370;
	t332 = qJD(4) * t382 - t338 * t370;
	t330 = qJD(4) * t383 - t336 * t370;
	t1 = [0, 0, t379 * t338 + t374 * (-t356 * t371 + t373 * t377) - t343 * t375, t382 * qJD(5) + t400 * (-t338 * t372 + (-t343 * t370 + t349 * t372) * qJD(4)) - t380 * t332, t332, 0; 0, 0, t379 * t336 + t374 * (-t373 * t378 - t399) + t401 * t375, t383 * qJD(5) + t400 * (-t336 * t372 + (t348 * t372 + t370 * t401) * qJD(4)) - t380 * t330, t330, 0; 0, 0, t379 * t344 + t374 * (t373 * t394 + (-t361 * t371 + t373 * t391) * t364) - t347 * t375, t381 * qJD(5) + t400 * (-t344 * t372 + (-t347 * t370 + t352 * t372) * qJD(4)) - t380 * t334, t334, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:43
	% EndTime: 2019-10-09 21:10:44
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (761->83), mult. (2148->148), div. (0->0), fcn. (2498->16), ass. (0->72)
	t478 = sin(pkin(12));
	t479 = sin(pkin(11));
	t465 = t479 * t478;
	t482 = cos(pkin(12));
	t483 = cos(pkin(11));
	t471 = t483 * t482;
	t485 = cos(pkin(6));
	t447 = -t485 * t471 + t465;
	t484 = cos(pkin(7));
	t445 = t447 * t484;
	t480 = sin(pkin(7));
	t481 = sin(pkin(6));
	t469 = t481 * t480;
	t456 = t483 * t469;
	t491 = t445 + t456;
	t467 = t479 * t482;
	t470 = t483 * t478;
	t448 = t485 * t467 + t470;
	t466 = t479 * t481;
	t490 = t448 * t484 - t480 * t466;
	t472 = t484 * t481;
	t489 = t482 * t472 + t485 * t480;
	t424 = t485 * t470 + t467;
	t440 = sin(qJ(3));
	t487 = cos(qJ(3));
	t407 = t424 * t440 + t491 * t487;
	t488 = t490 * t487;
	t436 = pkin(13) + qJ(6);
	t434 = sin(t436);
	t435 = cos(t436);
	t474 = t434 * r_i_i_C(1) + t435 * r_i_i_C(2);
	t458 = qJD(6) * t474;
	t468 = t481 * t478;
	t415 = t440 * t468 - t489 * t487;
	t486 = r_i_i_C(3) + pkin(10) + qJ(5);
	t477 = qJD(3) * t440;
	t476 = t424 * t487;
	t475 = t435 * r_i_i_C(1) - t434 * r_i_i_C(2);
	t408 = -t491 * t440 + t476;
	t417 = t447 * t480 - t483 * t472;
	t439 = sin(qJ(4));
	t441 = cos(qJ(4));
	t400 = t408 * t441 + t417 * t439;
	t464 = -t408 * t439 + t417 * t441;
	t425 = -t485 * t465 + t471;
	t410 = t425 * t487 - t490 * t440;
	t418 = t448 * t480 + t484 * t466;
	t402 = t410 * t441 + t418 * t439;
	t463 = -t410 * t439 + t418 * t441;
	t416 = t489 * t440 + t487 * t468;
	t423 = -t482 * t469 + t485 * t484;
	t412 = t416 * t441 + t423 * t439;
	t462 = -t416 * t439 + t423 * t441;
	t461 = cos(pkin(13)) * pkin(5) + pkin(4) + t475;
	t459 = qJD(6) * t475;
	t453 = -sin(pkin(13)) * pkin(5) - pkin(9) - t474;
	t450 = -t486 * t439 - t461 * t441 - pkin(3);
	t442 = -t439 * qJD(5) + t441 * t458 + (t461 * t439 - t486 * t441) * qJD(4);
	t414 = t416 * qJD(3);
	t413 = t415 * qJD(3);
	t409 = t425 * t440 + t488;
	t406 = t410 * qJD(3);
	t405 = t488 * qJD(3) + t425 * t477;
	t404 = -t456 * t477 + (-t440 * t445 + t476) * qJD(3);
	t403 = t407 * qJD(3);
	t398 = t462 * qJD(4) - t413 * t441;
	t397 = t412 * qJD(4) - t413 * t439;
	t396 = t463 * qJD(4) - t405 * t441;
	t395 = t402 * qJD(4) - t405 * t439;
	t394 = t464 * qJD(4) - t403 * t441;
	t393 = t400 * qJD(4) - t403 * t439;
	t1 = [0, 0, t453 * t405 + t450 * t406 + t442 * t409 + t410 * t459, t402 * qJD(5) - t461 * t395 + t486 * t396 - t463 * t458, t395, (-t396 * t434 + t406 * t435) * r_i_i_C(1) + (-t396 * t435 - t406 * t434) * r_i_i_C(2) + ((-t402 * t435 - t409 * t434) * r_i_i_C(1) + (t402 * t434 - t409 * t435) * r_i_i_C(2)) * qJD(6); 0, 0, t453 * t403 + t450 * t404 + t442 * t407 + t408 * t459, t400 * qJD(5) - t461 * t393 + t486 * t394 - t464 * t458, t393, (-t394 * t434 + t404 * t435) * r_i_i_C(1) + (-t394 * t435 - t404 * t434) * r_i_i_C(2) + ((-t400 * t435 - t407 * t434) * r_i_i_C(1) + (t400 * t434 - t407 * t435) * r_i_i_C(2)) * qJD(6); 0, 0, t453 * t413 + t450 * t414 + t442 * t415 + t416 * t459, t412 * qJD(5) - t461 * t397 + t486 * t398 - t462 * t458, t397, (-t398 * t434 + t414 * t435) * r_i_i_C(1) + (-t398 * t435 - t414 * t434) * r_i_i_C(2) + ((-t412 * t435 - t415 * t434) * r_i_i_C(1) + (t412 * t434 - t415 * t435) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end