% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:28
	% EndTime: 2019-10-09 21:14:28
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
	% StartTime: 2019-10-09 21:14:29
	% EndTime: 2019-10-09 21:14:30
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
	% StartTime: 2019-10-09 21:14:31
	% EndTime: 2019-10-09 21:14:32
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (599->84), mult. (1900->158), div. (0->0), fcn. (2202->14), ass. (0->66)
	t461 = sin(pkin(12));
	t462 = sin(pkin(11));
	t447 = t462 * t461;
	t465 = cos(pkin(12));
	t466 = cos(pkin(11));
	t453 = t466 * t465;
	t468 = cos(pkin(6));
	t431 = -t468 * t453 + t447;
	t467 = cos(pkin(7));
	t429 = t431 * t467;
	t463 = sin(pkin(7));
	t464 = sin(pkin(6));
	t451 = t464 * t463;
	t439 = t466 * t451;
	t474 = t429 + t439;
	t449 = t462 * t465;
	t452 = t466 * t461;
	t432 = t468 * t449 + t452;
	t448 = t462 * t464;
	t473 = t432 * t467 - t463 * t448;
	t454 = t467 * t464;
	t472 = t465 * t454 + t468 * t463;
	t412 = t468 * t452 + t449;
	t423 = sin(qJ(3));
	t469 = cos(qJ(3));
	t395 = t412 * t423 + t474 * t469;
	t471 = t473 * t469;
	t421 = sin(qJ(5));
	t424 = cos(qJ(5));
	t441 = qJD(5) * (t421 * r_i_i_C(1) + t424 * r_i_i_C(2));
	t450 = t464 * t461;
	t403 = t423 * t450 - t472 * t469;
	t470 = pkin(10) + r_i_i_C(3);
	t460 = qJD(3) * t423;
	t459 = qJD(5) * t421;
	t458 = qJD(5) * t424;
	t457 = t412 * t469;
	t396 = -t474 * t423 + t457;
	t405 = t431 * t463 - t466 * t454;
	t422 = sin(qJ(4));
	t425 = cos(qJ(4));
	t388 = t396 * t425 + t405 * t422;
	t446 = -t396 * t422 + t405 * t425;
	t413 = -t468 * t447 + t453;
	t398 = t413 * t469 - t473 * t423;
	t406 = t432 * t463 + t467 * t448;
	t390 = t398 * t425 + t406 * t422;
	t445 = -t398 * t422 + t406 * t425;
	t404 = t472 * t423 + t469 * t450;
	t411 = -t465 * t451 + t468 * t467;
	t400 = t404 * t425 + t411 * t422;
	t444 = -t404 * t422 + t411 * t425;
	t443 = t424 * r_i_i_C(1) - t421 * r_i_i_C(2) + pkin(4);
	t434 = -t470 * t422 - t443 * t425 - pkin(3);
	t426 = t425 * t441 + (t443 * t422 - t470 * t425) * qJD(4);
	t402 = t404 * qJD(3);
	t401 = t403 * qJD(3);
	t397 = t413 * t423 + t471;
	t394 = t398 * qJD(3);
	t393 = t471 * qJD(3) + t413 * t460;
	t392 = -t439 * t460 + (-t423 * t429 + t457) * qJD(3);
	t391 = t395 * qJD(3);
	t386 = t444 * qJD(4) - t401 * t425;
	t384 = t445 * qJD(4) - t393 * t425;
	t382 = t446 * qJD(4) - t391 * t425;
	t1 = [0, 0, (-t393 * t421 + t398 * t458) * r_i_i_C(1) + (-t393 * t424 - t398 * t459) * r_i_i_C(2) - t393 * pkin(9) + t434 * t394 + t426 * t397, t470 * t384 - t445 * t441 + t443 * (-t390 * qJD(4) + t393 * t422), (-t384 * t421 + t394 * t424) * r_i_i_C(1) + (-t384 * t424 - t394 * t421) * r_i_i_C(2) + ((-t390 * t424 - t397 * t421) * r_i_i_C(1) + (t390 * t421 - t397 * t424) * r_i_i_C(2)) * qJD(5), 0; 0, 0, (-t391 * t421 + t396 * t458) * r_i_i_C(1) + (-t391 * t424 - t396 * t459) * r_i_i_C(2) - t391 * pkin(9) + t434 * t392 + t426 * t395, t470 * t382 - t446 * t441 + t443 * (-t388 * qJD(4) + t391 * t422), (-t382 * t421 + t392 * t424) * r_i_i_C(1) + (-t382 * t424 - t392 * t421) * r_i_i_C(2) + ((-t388 * t424 - t395 * t421) * r_i_i_C(1) + (t388 * t421 - t395 * t424) * r_i_i_C(2)) * qJD(5), 0; 0, 0, (-t401 * t421 + t404 * t458) * r_i_i_C(1) + (-t401 * t424 - t404 * t459) * r_i_i_C(2) - t401 * pkin(9) + t434 * t402 + t426 * t403, t470 * t386 - t444 * t441 + t443 * (-t400 * qJD(4) + t401 * t422), (-t386 * t421 + t402 * t424) * r_i_i_C(1) + (-t386 * t424 - t402 * t421) * r_i_i_C(2) + ((-t400 * t424 - t403 * t421) * r_i_i_C(1) + (t400 * t421 - t403 * t424) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:14:32
	% EndTime: 2019-10-09 21:14:32
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (823->95), mult. (2546->168), div. (0->0), fcn. (2964->14), ass. (0->76)
	t482 = sin(pkin(12));
	t483 = sin(pkin(11));
	t469 = t483 * t482;
	t486 = cos(pkin(12));
	t487 = cos(pkin(11));
	t475 = t487 * t486;
	t489 = cos(pkin(6));
	t444 = -t489 * t475 + t469;
	t488 = cos(pkin(7));
	t442 = t444 * t488;
	t484 = sin(pkin(7));
	t485 = sin(pkin(6));
	t473 = t485 * t484;
	t453 = t487 * t473;
	t498 = t442 + t453;
	t471 = t483 * t486;
	t474 = t487 * t482;
	t445 = t489 * t471 + t474;
	t470 = t483 * t485;
	t497 = t445 * t488 - t484 * t470;
	t476 = t488 * t485;
	t496 = t486 * t476 + t489 * t484;
	t493 = pkin(5) + r_i_i_C(1);
	t434 = sin(qJ(5));
	t437 = cos(qJ(5));
	t450 = qJD(5) * (t437 * r_i_i_C(2) + t493 * t434);
	t423 = t489 * t474 + t471;
	t436 = sin(qJ(3));
	t492 = cos(qJ(3));
	t406 = t423 * t436 + t498 * t492;
	t495 = t497 * t492;
	t472 = t485 * t482;
	t414 = t436 * t472 - t496 * t492;
	t490 = r_i_i_C(3) + qJ(6) + pkin(10);
	t481 = qJD(3) * t436;
	t480 = qJD(5) * t434;
	t479 = qJD(5) * t437;
	t478 = t423 * t492;
	t402 = t406 * qJD(3);
	t438 = cos(qJ(4));
	t407 = -t498 * t436 + t478;
	t416 = t444 * t484 - t487 * t476;
	t435 = sin(qJ(4));
	t463 = -t407 * t435 + t416 * t438;
	t393 = t463 * qJD(4) - t402 * t438;
	t403 = -t453 * t481 + (-t436 * t442 + t478) * qJD(3);
	t468 = -t393 * t434 + t403 * t437;
	t424 = -t489 * t469 + t475;
	t404 = t495 * qJD(3) + t424 * t481;
	t409 = t424 * t492 - t497 * t436;
	t417 = t445 * t484 + t488 * t470;
	t462 = -t409 * t435 + t417 * t438;
	t395 = t462 * qJD(4) - t404 * t438;
	t405 = t409 * qJD(3);
	t467 = -t395 * t434 + t405 * t437;
	t412 = t414 * qJD(3);
	t415 = t496 * t436 + t492 * t472;
	t422 = -t486 * t473 + t489 * t488;
	t460 = -t415 * t435 + t422 * t438;
	t397 = t460 * qJD(4) - t412 * t438;
	t413 = t415 * qJD(3);
	t466 = -t397 * t434 + t413 * t437;
	t399 = t407 * t438 + t416 * t435;
	t465 = -t399 * t437 - t406 * t434;
	t401 = t409 * t438 + t417 * t435;
	t408 = t424 * t436 + t495;
	t464 = -t401 * t437 - t408 * t434;
	t411 = t415 * t438 + t422 * t435;
	t461 = -t411 * t437 - t414 * t434;
	t459 = -t434 * r_i_i_C(2) + t493 * t437 + pkin(4);
	t447 = -t490 * t435 - t459 * t438 - pkin(3);
	t439 = -t435 * qJD(6) + t438 * t450 + (t459 * t435 - t490 * t438) * qJD(4);
	t396 = t411 * qJD(4) - t412 * t435;
	t394 = t401 * qJD(4) - t404 * t435;
	t392 = t399 * qJD(4) - t402 * t435;
	t1 = [0, 0, (-t404 * t437 - t409 * t480) * r_i_i_C(2) - t404 * pkin(9) + t447 * t405 + t439 * t408 + t493 * (-t404 * t434 + t409 * t479), t401 * qJD(6) - t459 * t394 + t490 * t395 - t462 * t450, t467 * r_i_i_C(1) + (-t395 * t437 - t405 * t434) * r_i_i_C(2) + (t464 * r_i_i_C(1) + (t401 * t434 - t408 * t437) * r_i_i_C(2)) * qJD(5) + (t464 * qJD(5) + t467) * pkin(5), t394; 0, 0, (-t402 * t437 - t407 * t480) * r_i_i_C(2) - t402 * pkin(9) + t447 * t403 + t439 * t406 + t493 * (-t402 * t434 + t407 * t479), t399 * qJD(6) - t459 * t392 + t490 * t393 - t463 * t450, t468 * r_i_i_C(1) + (-t393 * t437 - t403 * t434) * r_i_i_C(2) + (t465 * r_i_i_C(1) + (t399 * t434 - t406 * t437) * r_i_i_C(2)) * qJD(5) + (t465 * qJD(5) + t468) * pkin(5), t392; 0, 0, (-t412 * t437 - t415 * t480) * r_i_i_C(2) - t412 * pkin(9) + t447 * t413 + t439 * t414 + t493 * (-t412 * t434 + t415 * t479), t411 * qJD(6) - t459 * t396 + t490 * t397 - t460 * t450, t466 * r_i_i_C(1) + (-t397 * t437 - t413 * t434) * r_i_i_C(2) + (t461 * r_i_i_C(1) + (t411 * t434 - t414 * t437) * r_i_i_C(2)) * qJD(5) + (t461 * qJD(5) + t466) * pkin(5), t396;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end