% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRRP2
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
% Datum: 2019-10-09 21:16
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
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
	% StartTime: 2019-10-09 21:16:22
	% EndTime: 2019-10-09 21:16:23
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
	% StartTime: 2019-10-09 21:16:24
	% EndTime: 2019-10-09 21:16:25
	% DurationCPUTime: 0.71s
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
	% StartTime: 2019-10-09 21:16:27
	% EndTime: 2019-10-09 21:16:27
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (1160->104), mult. (3610->180), div. (0->0), fcn. (4263->14), ass. (0->81)
	t538 = sin(pkin(12));
	t539 = sin(pkin(11));
	t521 = t539 * t538;
	t542 = cos(pkin(12));
	t543 = cos(pkin(11));
	t527 = t543 * t542;
	t545 = cos(pkin(6));
	t498 = -t545 * t527 + t521;
	t544 = cos(pkin(7));
	t496 = t498 * t544;
	t540 = sin(pkin(7));
	t541 = sin(pkin(6));
	t525 = t541 * t540;
	t511 = t543 * t525;
	t553 = t496 + t511;
	t523 = t539 * t542;
	t526 = t543 * t538;
	t499 = t545 * t523 + t526;
	t522 = t539 * t541;
	t552 = t499 * t544 - t540 * t522;
	t528 = t544 * t541;
	t551 = t542 * t528 + t545 * t540;
	t480 = t545 * t526 + t523;
	t491 = sin(qJ(3));
	t547 = cos(qJ(3));
	t462 = t480 * t491 + t553 * t547;
	t550 = t552 * t547;
	t524 = t541 * t538;
	t471 = t491 * t524 - t551 * t547;
	t549 = pkin(10) + r_i_i_C(2);
	t548 = r_i_i_C(1) + pkin(5);
	t546 = r_i_i_C(3) + qJ(6);
	t489 = sin(qJ(5));
	t493 = cos(qJ(4));
	t537 = t489 * t493;
	t536 = qJD(3) * t491;
	t490 = sin(qJ(4));
	t535 = qJD(4) * t490;
	t534 = qJD(5) * t493;
	t533 = t480 * t547;
	t458 = t462 * qJD(3);
	t532 = t462 * t534 - t458;
	t481 = -t545 * t521 + t527;
	t460 = t550 * qJD(3) + t481 * t536;
	t464 = t481 * t491 + t550;
	t531 = t464 * t534 - t460;
	t469 = t471 * qJD(3);
	t530 = t471 * t534 - t469;
	t463 = -t553 * t491 + t533;
	t473 = t498 * t540 - t543 * t528;
	t453 = t463 * t493 + t473 * t490;
	t492 = cos(qJ(5));
	t520 = t453 * t492 + t462 * t489;
	t465 = t481 * t547 - t552 * t491;
	t474 = t499 * t540 + t544 * t522;
	t455 = t465 * t493 + t474 * t490;
	t519 = t455 * t492 + t464 * t489;
	t518 = -t463 * t490 + t473 * t493;
	t517 = -t465 * t490 + t474 * t493;
	t472 = t551 * t491 + t547 * t524;
	t479 = -t542 * t525 + t545 * t544;
	t467 = t472 * t493 + t479 * t490;
	t516 = t467 * t492 + t471 * t489;
	t515 = -t472 * t490 + t479 * t493;
	t513 = -t493 * pkin(4) - t549 * t490 - pkin(3);
	t508 = qJD(4) * (pkin(4) * t490 - t549 * t493);
	t507 = t546 * t489 + t548 * t492 + pkin(4);
	t459 = -t511 * t536 + (-t491 * t496 + t533) * qJD(3);
	t506 = qJD(5) * t463 - t459 * t493 + t462 * t535;
	t461 = t465 * qJD(3);
	t505 = qJD(5) * t465 - t461 * t493 + t464 * t535;
	t470 = t472 * qJD(3);
	t504 = qJD(5) * t472 - t470 * t493 + t471 * t535;
	t501 = qJD(6) * t489 + (-t548 * t489 + t546 * t492) * qJD(5);
	t451 = t515 * qJD(4) - t469 * t493;
	t449 = t517 * qJD(4) - t460 * t493;
	t447 = t518 * qJD(4) - t458 * t493;
	t442 = t516 * qJD(5) + t451 * t489 - t470 * t492;
	t436 = t519 * qJD(5) + t449 * t489 - t461 * t492;
	t434 = t520 * qJD(5) + t447 * t489 - t459 * t492;
	t1 = [0, 0, -(t464 * t537 + t465 * t492) * qJD(6) - t460 * pkin(9) + t548 * (t531 * t489 + t505 * t492) + t546 * (t505 * t489 - t531 * t492) + t464 * t508 + t513 * t461, t549 * t449 + t501 * t517 + t507 * (-t455 * qJD(4) + t460 * t490), t519 * qJD(6) + t546 * (t449 * t492 + t461 * t489 + (-t455 * t489 + t464 * t492) * qJD(5)) - t548 * t436, t436; 0, 0, -(t462 * t537 + t463 * t492) * qJD(6) - t458 * pkin(9) + t548 * (t532 * t489 + t506 * t492) + t546 * (t506 * t489 - t532 * t492) + t462 * t508 + t513 * t459, t549 * t447 + t501 * t518 + t507 * (-t453 * qJD(4) + t458 * t490), t520 * qJD(6) + t546 * (t447 * t492 + t459 * t489 + (-t453 * t489 + t462 * t492) * qJD(5)) - t548 * t434, t434; 0, 0, -(t471 * t537 + t472 * t492) * qJD(6) - t469 * pkin(9) + t548 * (t530 * t489 + t504 * t492) + t546 * (t504 * t489 - t530 * t492) + t471 * t508 + t513 * t470, t549 * t451 + t501 * t515 + t507 * (-t467 * qJD(4) + t469 * t490), t516 * qJD(6) + t546 * (t451 * t492 + t470 * t489 + (-t467 * t489 + t471 * t492) * qJD(5)) - t548 * t442, t442;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end