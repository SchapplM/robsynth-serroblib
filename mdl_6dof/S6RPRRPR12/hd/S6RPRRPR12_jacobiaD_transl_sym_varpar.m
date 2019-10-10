% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR12_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR12_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
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
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
	t114 = sin(pkin(6));
	t124 = t114 * (r_i_i_C(3) + qJ(2));
	t116 = cos(pkin(6));
	t117 = sin(qJ(1));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t116 * t118;
	t120 = qJD(1) * t114;
	t119 = t114 * qJD(2);
	t115 = cos(pkin(12));
	t113 = sin(pkin(12));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (103->49), mult. (358->92), div. (0->0), fcn. (352->10), ass. (0->47)
	t248 = cos(pkin(7));
	t278 = r_i_i_C(3) + pkin(9);
	t279 = t278 * t248 + qJ(2);
	t250 = sin(qJ(3));
	t277 = t250 * r_i_i_C(2);
	t246 = sin(pkin(6));
	t251 = sin(qJ(1));
	t276 = t246 * t251;
	t253 = cos(qJ(1));
	t275 = t246 * t253;
	t274 = t248 * t250;
	t252 = cos(qJ(3));
	t273 = t248 * t252;
	t244 = sin(pkin(12));
	t272 = t251 * t244;
	t247 = cos(pkin(12));
	t271 = t251 * t247;
	t270 = t253 * t244;
	t269 = t253 * t247;
	t268 = qJD(1) * t251;
	t267 = qJD(1) * t253;
	t245 = sin(pkin(7));
	t266 = qJD(3) * t245;
	t265 = t278 * t245;
	t249 = cos(pkin(6));
	t264 = t249 * t272;
	t263 = t246 * t268;
	t262 = t246 * t267;
	t261 = t266 * t275;
	t260 = r_i_i_C(1) * t250 + r_i_i_C(2) * t252;
	t236 = -qJD(1) * t264 + t247 * t267;
	t237 = -t249 * t269 + t272;
	t259 = qJD(3) * t237 * t248 - t236;
	t258 = t260 * t245;
	t256 = t249 * t271 + t270;
	t257 = t245 * t276 - t248 * t256;
	t238 = t249 * t270 + t271;
	t233 = t237 * qJD(1);
	t255 = t233 * t248 + t245 * t262;
	t235 = t256 * qJD(1);
	t254 = -qJD(3) * t238 - t235 * t248 + t245 * t263;
	t241 = t252 * t261;
	t240 = -t264 + t269;
	t234 = t238 * qJD(1);
	t232 = -t234 * t252 + t255 * t250 + (-t240 * t250 + t257 * t252) * qJD(3);
	t231 = t234 * t250 + t255 * t252 + (-t240 * t252 - t257 * t250) * qJD(3);
	t1 = [-pkin(1) * t267 + t241 * r_i_i_C(1) + (-t252 * r_i_i_C(1) - pkin(2) + t277) * t236 + (t260 * t248 - t265) * t235 + ((t237 * t273 + t238 * t250) * r_i_i_C(1) + (-t237 * t274 + t238 * t252) * r_i_i_C(2)) * qJD(3) + ((-t266 * t277 + qJD(2)) * t253 + (-t258 - t279) * t268) * t246, t262, t231 * r_i_i_C(1) - t232 * r_i_i_C(2), 0, 0, 0; qJD(2) * t276 - t234 * pkin(2) + t232 * r_i_i_C(1) + t231 * r_i_i_C(2) - t233 * t265 + (-t251 * pkin(1) + t279 * t275) * qJD(1), t263, t241 * r_i_i_C(2) + (t254 * r_i_i_C(1) + t259 * r_i_i_C(2)) * t252 + ((t259 + t261) * r_i_i_C(1) - t254 * r_i_i_C(2)) * t250, 0, 0, 0; 0, 0, (-t249 * t258 + ((-t244 * t252 - t247 * t274) * r_i_i_C(1) + (t244 * t250 - t247 * t273) * r_i_i_C(2)) * t246) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:28
	% EndTime: 2019-10-10 01:43:28
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (393->78), mult. (1301->142), div. (0->0), fcn. (1376->12), ass. (0->64)
	t408 = sin(qJ(3));
	t402 = sin(pkin(7));
	t403 = sin(pkin(6));
	t412 = cos(qJ(1));
	t442 = t403 * t412;
	t431 = t402 * t442;
	t405 = cos(pkin(7));
	t404 = cos(pkin(12));
	t406 = cos(pkin(6));
	t440 = t406 * t412;
	t401 = sin(pkin(12));
	t409 = sin(qJ(1));
	t445 = t401 * t409;
	t421 = t404 * t440 - t445;
	t454 = t421 * t405;
	t424 = -t454 + t431;
	t463 = t424 * t408;
	t432 = t406 * t445;
	t438 = qJD(1) * t412;
	t389 = -qJD(1) * t432 + t404 * t438;
	t439 = t409 * t404;
	t392 = t401 * t440 + t439;
	t411 = cos(qJ(3));
	t420 = t401 * t412 + t406 * t439;
	t388 = t420 * qJD(1);
	t443 = t403 * t409;
	t429 = qJD(1) * t443;
	t418 = -t388 * t405 + t402 * t429;
	t428 = t411 * t431;
	t368 = -qJD(3) * t428 + (-qJD(3) * t392 + t418) * t408 - (-qJD(3) * t454 - t389) * t411;
	t449 = t388 * t402;
	t380 = t405 * t429 + t449;
	t407 = sin(qJ(4));
	t410 = cos(qJ(4));
	t462 = t368 * t407 - t380 * t410;
	t461 = -t368 * t410 - t380 * t407;
	t373 = -t392 * t411 + t463;
	t453 = t421 * t402 + t405 * t442;
	t460 = -t373 * t407 + t453 * t410;
	t459 = t373 * t410 + t453 * t407;
	t441 = t404 * t405;
	t444 = t402 * t406;
	t452 = (-t401 * t408 + t411 * t441) * t403 + t411 * t444;
	t382 = (t401 * t411 + t408 * t441) * t403 + t408 * t444;
	t450 = r_i_i_C(3) + pkin(10);
	t436 = t403 * qJD(2);
	t427 = t407 * r_i_i_C(1) + t410 * r_i_i_C(2);
	t425 = r_i_i_C(1) * t410 - r_i_i_C(2) * t407 + pkin(3);
	t423 = t402 * t443 - t405 * t420;
	t419 = qJD(4) * t427;
	t394 = t404 * t412 - t432;
	t415 = -t394 * t408 + t423 * t411;
	t375 = t394 * t411 + t423 * t408;
	t413 = t373 * qJD(3) - t389 * t408 + t418 * t411;
	t390 = -t402 * t403 * t404 + t405 * t406;
	t387 = t392 * qJD(1);
	t385 = t402 * t420 + t405 * t443;
	t378 = t453 * qJD(1);
	t376 = t452 * qJD(3);
	t366 = qJD(1) * t463 + t415 * qJD(3) - t387 * t411;
	t365 = t375 * qJD(3) - t387 * t408 + (t411 * t454 - t428) * qJD(1);
	t364 = t366 * t410 + t378 * t407 + (-t375 * t407 + t385 * t410) * qJD(4);
	t363 = -t366 * t407 + t378 * t410 + (-t375 * t410 - t385 * t407) * qJD(4);
	t1 = [t461 * r_i_i_C(1) + t462 * r_i_i_C(2) - t368 * pkin(3) - t389 * pkin(2) - pkin(9) * t449 + t412 * t436 + t450 * t413 + (t460 * r_i_i_C(1) - t459 * r_i_i_C(2)) * qJD(4) + (-t412 * pkin(1) + (-pkin(9) * t405 - qJ(2)) * t443) * qJD(1), t403 * t438, -t425 * t365 + t450 * t366 - t415 * t419, r_i_i_C(1) * t363 - r_i_i_C(2) * t364, 0, 0; t409 * t436 - t387 * pkin(2) + t366 * pkin(3) + t364 * r_i_i_C(1) + t363 * r_i_i_C(2) + t450 * t365 + (-t409 * pkin(1) + pkin(9) * t453 + qJ(2) * t442) * qJD(1), t429, t450 * t368 - (-t392 * t408 - t424 * t411) * t419 + t425 * t413, -t462 * r_i_i_C(1) + t461 * r_i_i_C(2) + (t459 * r_i_i_C(1) + t460 * r_i_i_C(2)) * qJD(4), 0, 0; 0, 0, -t425 * t382 * qJD(3) + t450 * t376 - t452 * t419, -t427 * t376 + ((-t382 * t410 - t390 * t407) * r_i_i_C(1) + (t382 * t407 - t390 * t410) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:30
	% EndTime: 2019-10-10 01:43:31
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (746->86), mult. (2410->145), div. (0->0), fcn. (2622->12), ass. (0->70)
	t468 = cos(pkin(12));
	t470 = cos(pkin(6));
	t465 = sin(pkin(12));
	t473 = sin(qJ(1));
	t510 = t473 * t465;
	t500 = t470 * t510;
	t476 = cos(qJ(1));
	t506 = qJD(1) * t476;
	t453 = -qJD(1) * t500 + t468 * t506;
	t508 = t476 * t465;
	t509 = t473 * t468;
	t456 = t470 * t508 + t509;
	t472 = sin(qJ(3));
	t475 = cos(qJ(3));
	t486 = t470 * t509 + t508;
	t452 = t486 * qJD(1);
	t466 = sin(pkin(7));
	t469 = cos(pkin(7));
	t467 = sin(pkin(6));
	t513 = t467 * t473;
	t497 = qJD(1) * t513;
	t485 = -t452 * t469 + t466 * t497;
	t512 = t467 * t476;
	t499 = t466 * t512;
	t496 = t475 * t499;
	t507 = t476 * t468;
	t487 = t470 * t507 - t510;
	t528 = t487 * t469;
	t429 = -qJD(3) * t496 + (-qJD(3) * t456 + t485) * t472 - (-qJD(3) * t528 - t453) * t475;
	t518 = t452 * t466;
	t444 = t469 * t497 + t518;
	t471 = sin(qJ(4));
	t474 = cos(qJ(4));
	t490 = -t528 + t499;
	t537 = t490 * t472;
	t437 = -t456 * t475 + t537;
	t527 = t487 * t466 + t469 * t512;
	t533 = t437 * t474 + t527 * t471;
	t539 = t533 * qJD(4) - t429 * t471 + t444 * t474;
	t534 = t437 * t471 - t527 * t474;
	t538 = t534 * qJD(4) + t429 * t474 + t444 * t471;
	t511 = t468 * t469;
	t514 = t466 * t470;
	t526 = (-t465 * t472 + t475 * t511) * t467 + t475 * t514;
	t458 = -t500 + t507;
	t489 = t466 * t513 - t469 * t486;
	t439 = t458 * t475 + t489 * t472;
	t449 = t466 * t486 + t469 * t513;
	t524 = -t439 * t471 + t449 * t474;
	t446 = (t465 * t475 + t472 * t511) * t467 + t472 * t514;
	t521 = r_i_i_C(3) + qJ(5);
	t522 = r_i_i_C(2) - pkin(4);
	t483 = t521 * t471 - t522 * t474 + pkin(3);
	t523 = r_i_i_C(1) + pkin(10);
	t504 = t467 * qJD(2);
	t492 = t439 * t474 + t449 * t471;
	t454 = -t467 * t468 * t466 + t470 * t469;
	t491 = t446 * t474 + t454 * t471;
	t481 = -t458 * t472 + t489 * t475;
	t479 = qJD(5) * t471 + (t522 * t471 + t521 * t474) * qJD(4);
	t478 = t527 * qJD(1);
	t477 = t437 * qJD(3) - t453 * t472 + t485 * t475;
	t451 = t456 * qJD(1);
	t441 = t526 * qJD(3);
	t432 = t491 * qJD(4) + t441 * t471;
	t427 = qJD(1) * t537 + t481 * qJD(3) - t451 * t475;
	t426 = t439 * qJD(3) - t451 * t472 + (t475 * t528 - t496) * qJD(1);
	t421 = t524 * qJD(4) + t427 * t474 + t471 * t478;
	t420 = t492 * qJD(4) + t427 * t471 - t474 * t478;
	t1 = [t534 * qJD(5) - t429 * pkin(3) - t453 * pkin(2) - pkin(9) * t518 + t476 * t504 + t523 * t477 + t522 * t538 + t521 * t539 + (-t476 * pkin(1) + (-pkin(9) * t469 - qJ(2)) * t513) * qJD(1), t467 * t506, -t483 * t426 + t523 * t427 + t479 * t481, t492 * qJD(5) + t522 * t420 + t521 * t421, t420, 0; -t524 * qJD(5) + t427 * pkin(3) - t451 * pkin(2) + t473 * t504 + t523 * t426 - t522 * t421 + t521 * t420 + (-t473 * pkin(1) + pkin(9) * t527 + qJ(2) * t512) * qJD(1), t497, t523 * t429 + t479 * (-t456 * t472 - t490 * t475) + t483 * t477, -t533 * qJD(5) + t521 * t538 - t522 * t539, -t539, 0; 0, 0, -t483 * t446 * qJD(3) + t523 * t441 + t479 * t526, t491 * qJD(5) + t521 * (t441 * t474 + (-t446 * t471 + t454 * t474) * qJD(4)) + t522 * t432, t432, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:31
	% EndTime: 2019-10-10 01:43:33
	% DurationCPUTime: 1.56s
	% Computational Cost: add. (1445->121), mult. (4618->199), div. (0->0), fcn. (5138->14), ass. (0->86)
	t572 = sin(qJ(3));
	t635 = sin(pkin(7));
	t636 = sin(pkin(6));
	t615 = t636 * t635;
	t643 = cos(qJ(1));
	t599 = t643 * t615;
	t637 = cos(pkin(12));
	t639 = cos(pkin(6));
	t619 = t639 * t637;
	t634 = sin(pkin(12));
	t641 = sin(qJ(1));
	t554 = -t643 * t619 + t641 * t634;
	t638 = cos(pkin(7));
	t627 = t638 * t554;
	t617 = t639 * t634;
	t555 = t643 * t617 + t641 * t637;
	t642 = cos(qJ(3));
	t628 = t555 * t642;
	t534 = (t599 + t627) * t572 - t628;
	t616 = t638 * t636;
	t545 = t554 * t635 - t643 * t616;
	t571 = sin(qJ(4));
	t574 = cos(qJ(4));
	t525 = t534 * t571 + t545 * t574;
	t589 = t642 * t599;
	t650 = -t642 * t627 - t589;
	t531 = t555 * t572 - t650;
	t570 = sin(qJ(6));
	t573 = cos(qJ(6));
	t662 = t525 * t570 - t531 * t573;
	t661 = t525 * t573 + t531 * t570;
	t556 = -t641 * t617 + t643 * t637;
	t552 = t556 * qJD(1);
	t595 = t641 * t615;
	t583 = t641 * t619 + t643 * t634;
	t630 = t583 * qJD(1);
	t647 = t630 * t638;
	t520 = -t650 * qJD(3) - t552 * t642 + (-qJD(1) * t595 + t555 * qJD(3) + t647) * t572;
	t596 = t641 * t616;
	t611 = t630 * t635;
	t581 = qJD(1) * t596 + t611;
	t660 = t525 * qJD(4) - t520 * t574 + t581 * t571;
	t655 = t534 * t574 - t545 * t571;
	t659 = t655 * qJD(4) + t520 * t571 + t581 * t574;
	t656 = t637 * t616 + t639 * t635;
	t588 = t642 * t595;
	t593 = t572 * t599;
	t517 = t642 * t647 - qJD(1) * t588 + t552 * t572 - (t572 * t627 + t593 - t628) * qJD(3);
	t648 = pkin(5) + pkin(10);
	t580 = t583 * t638;
	t536 = t556 * t642 + (-t580 + t595) * t572;
	t546 = t583 * t635 + t596;
	t646 = -t536 * t571 + t546 * t574;
	t629 = pkin(4) + pkin(11) + r_i_i_C(3);
	t645 = -t556 * t572 - t642 * t580 + t588;
	t614 = t636 * t634;
	t542 = t572 * t614 - t656 * t642;
	t575 = qJD(1) * t545;
	t622 = t573 * r_i_i_C(1) - t570 * r_i_i_C(2);
	t592 = t622 * qJD(6) + qJD(5);
	t625 = t643 * t636;
	t623 = t641 * t636;
	t621 = -t570 * r_i_i_C(1) - t573 * r_i_i_C(2);
	t609 = t536 * t574 + t546 * t571;
	t543 = t656 * t572 + t642 * t614;
	t553 = -t637 * t615 + t639 * t638;
	t608 = t543 * t574 + t553 * t571;
	t529 = t543 * t571 - t553 * t574;
	t606 = qJ(5) - t621;
	t605 = qJD(1) * t625;
	t603 = t622 + t648;
	t602 = qJD(6) * t621;
	t582 = -t606 * t571 - t629 * t574 - pkin(3);
	t578 = qJD(1) * t627;
	t576 = -t592 * t571 + (t629 * t571 - t606 * t574) * qJD(4);
	t551 = t555 * qJD(1);
	t540 = t543 * qJD(3);
	t539 = t542 * qJD(3);
	t521 = t608 * qJD(4) - t539 * t571;
	t516 = qJD(1) * t593 + t645 * qJD(3) - t551 * t642 + t572 * t578;
	t515 = -qJD(1) * t589 + t536 * qJD(3) - t551 * t572 - t642 * t578;
	t510 = t646 * qJD(4) + t516 * t574 - t571 * t575;
	t509 = t609 * qJD(4) + t516 * t571 + t574 * t575;
	t508 = t509 * t570 + t515 * t573 + (t570 * t645 - t573 * t646) * qJD(6);
	t507 = t509 * t573 - t515 * t570 + (t570 * t646 + t573 * t645) * qJD(6);
	t1 = [t525 * qJD(5) + t520 * pkin(3) - t552 * pkin(2) - pkin(9) * t611 + qJD(2) * t625 + t606 * t659 - t603 * t517 + (t661 * r_i_i_C(1) - t662 * r_i_i_C(2)) * qJD(6) - t629 * t660 + (-t643 * pkin(1) - pkin(9) * t596 - qJ(2) * t623) * qJD(1), t605, t582 * t515 + t603 * t516 + t536 * t602 - t576 * t645, -t629 * t509 + t606 * t510 + t592 * t609, t509, t507 * r_i_i_C(1) - t508 * r_i_i_C(2); -qJD(1) * t641 * pkin(1) - t551 * pkin(2) + t516 * pkin(3) - pkin(9) * t575 + t508 * r_i_i_C(1) + t507 * r_i_i_C(2) + qJ(2) * t605 + t509 * qJ(5) + qJD(2) * t623 - qJD(5) * t646 + t629 * t510 + t648 * t515, qJD(1) * t623, t582 * t517 - t520 * t603 + t576 * t531 - t534 * t602, -t592 * t655 + t606 * t660 + t629 * t659, -t659, (-t517 * t570 - t573 * t659) * r_i_i_C(1) + (-t517 * t573 + t570 * t659) * r_i_i_C(2) + (t662 * r_i_i_C(1) + t661 * r_i_i_C(2)) * qJD(6); 0, 0, -t603 * t539 + t582 * t540 + t576 * t542 + t543 * t602, t592 * t608 + t606 * (-t529 * qJD(4) - t539 * t574) - t629 * t521, t521, (t521 * t573 - t540 * t570) * r_i_i_C(1) + (-t521 * t570 - t540 * t573) * r_i_i_C(2) + ((-t529 * t570 - t542 * t573) * r_i_i_C(1) + (-t529 * t573 + t542 * t570) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end