% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PRPPR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:44
% EndTime: 2019-12-31 17:38:48
% DurationCPUTime: 4.54s
% Computational Cost: add. (12645->383), mult. (20883->461), div. (0->0), fcn. (11122->8), ass. (0->233)
t551 = sin(pkin(7));
t556 = sin(qJ(2));
t558 = cos(qJ(2));
t560 = qJD(2) ^ 2;
t520 = qJDD(2) * t556 + t558 * t560;
t553 = cos(pkin(7));
t524 = g(1) * t551 - t553 * g(2);
t594 = -pkin(5) * t520 + t558 * t524;
t648 = t551 * t594;
t647 = t553 * t594;
t550 = sin(pkin(8));
t552 = cos(pkin(8));
t517 = -t550 * qJDD(2) + t552 * t560;
t518 = qJDD(2) * t552 + t550 * t560;
t457 = t558 * t517 + t518 * t556;
t515 = qJDD(4) + t524;
t468 = qJ(4) * t517 + t515 * t552;
t574 = qJ(4) * t518 + t515 * t550;
t646 = -pkin(5) * t457 + t468 * t558 + t556 * t574;
t525 = g(1) * t553 + g(2) * t551;
t548 = g(3) - qJDD(1);
t495 = -t558 * t525 - t556 * t548;
t602 = (qJD(3) * qJD(2));
t543 = 2 * t602;
t601 = qJDD(2) * qJ(3);
t569 = t543 + t601 + t495;
t631 = pkin(2) + pkin(3);
t451 = -t631 * t560 + t569;
t549 = qJDD(2) * pkin(2);
t494 = -t525 * t556 + t558 * t548;
t575 = -qJDD(3) - t494;
t471 = -t560 * qJ(3) - t549 - t575;
t561 = -qJDD(2) * pkin(3) + t471;
t416 = t451 * t550 - t552 * t561;
t417 = t552 * t451 + t550 * t561;
t376 = t416 * t552 - t417 * t550;
t573 = t416 * t550 + t417 * t552;
t357 = t376 * t556 - t558 * t573;
t645 = t376 * t558 + t556 * t573;
t585 = -t517 * t556 + t558 * t518;
t400 = pkin(5) * t585 + t468 * t556 - t558 * t574;
t415 = -pkin(4) * t560 - qJDD(2) * pkin(6) + t417;
t555 = sin(qJ(5));
t557 = cos(qJ(5));
t392 = t415 * t555 - t557 * t515;
t393 = t415 * t557 + t515 * t555;
t371 = t555 * t392 + t557 * t393;
t502 = t553 * t524;
t640 = -t525 * t551 + t502;
t414 = qJDD(2) * pkin(4) - t560 * pkin(6) + t416;
t363 = t371 * t550 - t414 * t552;
t364 = t371 * t552 + t414 * t550;
t606 = pkin(4) * t414 - pkin(6) * t371;
t639 = qJ(3) * t364 - t631 * t363 + t606;
t638 = qJ(3) * t573 + t376 * t631;
t547 = t557 ^ 2;
t541 = t547 * t560;
t559 = qJD(5) ^ 2;
t531 = -t541 - t559;
t532 = t555 * t560 * t557;
t526 = qJDD(5) + t532;
t615 = t526 * t555;
t486 = t531 * t557 - t615;
t603 = qJD(2) * qJD(5);
t536 = t555 * t603;
t598 = t557 * qJDD(2);
t512 = 0.2e1 * t536 - t598;
t433 = t486 * t550 + t512 * t552;
t435 = t486 * t552 - t512 * t550;
t404 = t557 * t414;
t596 = -pkin(4) * t512 - pkin(6) * t486 + t404;
t637 = qJ(3) * t435 - t631 * t433 + t596;
t546 = t555 ^ 2;
t611 = t546 * t560;
t529 = -t559 - t611;
t527 = qJDD(5) - t532;
t612 = t527 * t557;
t488 = -t529 * t555 - t612;
t595 = t557 * t603;
t599 = t555 * qJDD(2);
t510 = 0.2e1 * t595 + t599;
t434 = t488 * t550 + t510 * t552;
t436 = t488 * t552 - t510 * t550;
t403 = t555 * t414;
t597 = pkin(4) * t510 + pkin(6) * t488 + t403;
t636 = qJ(3) * t436 - t631 * t434 - t597;
t604 = -t546 - t547;
t519 = t604 * qJDD(2);
t522 = t541 + t611;
t467 = t519 * t550 + t522 * t552;
t470 = t519 * t552 - t522 * t550;
t580 = pkin(4) * t522 + pkin(6) * t519 + t371;
t635 = qJ(3) * t470 - t631 * t467 - t580;
t634 = qJ(3) * t518 + t631 * t517 + t417;
t633 = -qJ(3) * t517 + t631 * t518 + t416;
t630 = pkin(1) * t520;
t521 = qJDD(2) * t558 - t556 * t560;
t629 = pkin(1) * t521;
t370 = t392 * t557 - t393 * t555;
t628 = pkin(4) * t370;
t626 = pkin(6) * t552;
t625 = qJ(1) * t520;
t624 = qJ(1) * t521;
t622 = qJ(4) * t364;
t621 = qJ(4) * t376;
t620 = qJ(4) * t573;
t617 = t524 * t551;
t614 = t526 * t557;
t613 = t527 * t555;
t610 = t551 * t520;
t609 = t551 * t548;
t608 = t553 * t520;
t607 = t553 * t548;
t474 = pkin(5) * t521 + t524 * t556;
t605 = -qJ(1) * t608 - t551 * t474;
t600 = t553 * qJDD(2);
t422 = -t467 * t558 + t470 * t556;
t352 = -pkin(1) * t422 - t635;
t423 = t467 * t556 + t470 * t558;
t593 = qJ(1) * t423 + t352;
t378 = -pkin(1) * t457 - t634;
t592 = qJ(1) * t585 + t378;
t379 = -pkin(1) * t585 - t633;
t591 = -qJ(1) * t457 + t379;
t576 = 0.2e1 * t601 + t495;
t437 = -t576 - (2 * t602) - t630;
t590 = t437 + t624;
t452 = t495 + t630;
t589 = -t452 + t624;
t588 = -qJ(4) * t363 - t550 * t628;
t456 = -pkin(2) * t560 + t569;
t587 = t558 * t456 + t471 * t556;
t586 = t494 * t556 + t558 * t495;
t583 = -t553 * t525 - t617;
t582 = t550 * t532;
t581 = t552 * t532;
t577 = -pkin(2) * t471 + qJ(3) * t456;
t419 = t456 * t556 - t471 * t558;
t431 = t494 * t558 - t495 * t556;
t572 = -pkin(4) * t552 - pkin(6) * t550 - pkin(3);
t568 = 0.2e1 * t549 + t575;
t482 = t531 * t555 + t614;
t381 = -pkin(4) * t482 + t392;
t385 = -pkin(6) * t482 + t403;
t566 = -qJ(4) * t433 - t381 * t550 + t552 * t385;
t484 = t529 * t557 - t613;
t382 = -pkin(4) * t484 + t393;
t386 = -pkin(6) * t484 + t404;
t565 = -qJ(4) * t434 - t382 * t550 + t552 * t386;
t563 = -qJ(4) * t435 - t381 * t552 - t385 * t550;
t562 = -qJ(4) * t436 - t382 * t552 - t386 * t550;
t538 = t551 * qJDD(2);
t530 = t541 - t559;
t528 = t559 - t611;
t523 = -t541 + t611;
t516 = pkin(1) * t524;
t513 = t536 - t598;
t511 = -t595 - t599;
t508 = t604 * t603;
t501 = t553 * t521;
t500 = t551 * t521;
t493 = t511 * t557 + t546 * t603;
t492 = -t513 * t555 + t547 * t603;
t491 = qJDD(5) * t550 + t508 * t552;
t490 = -t552 * qJDD(5) + t508 * t550;
t487 = -t528 * t555 + t614;
t485 = t530 * t557 - t613;
t483 = t528 * t557 + t615;
t481 = t530 * t555 + t612;
t480 = (t511 - t595) * t555;
t479 = (t513 + t536) * t557;
t464 = t553 * t474;
t455 = t510 * t555 + t512 * t557;
t454 = -t510 * t557 + t512 * t555;
t453 = t494 - t629;
t447 = t493 * t552 - t582;
t446 = t492 * t552 + t582;
t445 = t493 * t550 + t581;
t444 = t492 * t550 - t581;
t443 = t487 * t552 - t550 * t599;
t442 = t485 * t552 - t550 * t598;
t441 = t487 * t550 + t552 * t599;
t440 = t485 * t550 + t552 * t598;
t438 = -t568 - t629;
t428 = t455 * t552 + t523 * t550;
t427 = t455 * t550 - t523 * t552;
t426 = t490 * t556 + t491 * t558;
t425 = -t490 * t558 + t491 * t556;
t424 = pkin(5) * t586 + t516;
t412 = t445 * t556 + t447 * t558;
t411 = t444 * t556 + t446 * t558;
t410 = -t445 * t558 + t447 * t556;
t409 = -t444 * t558 + t446 * t556;
t408 = t441 * t556 + t443 * t558;
t407 = t440 * t556 + t442 * t558;
t406 = -t441 * t558 + t443 * t556;
t405 = -t440 * t558 + t442 * t556;
t398 = t434 * t556 + t436 * t558;
t397 = t433 * t556 + t435 * t558;
t396 = -t434 * t558 + t436 * t556;
t395 = -t433 * t558 + t435 * t556;
t394 = -pkin(5) * t419 + (-pkin(2) * t556 + qJ(3) * t558) * t524;
t389 = t427 * t556 + t428 * t558;
t388 = -t427 * t558 + t428 * t556;
t387 = pkin(5) * t587 + t516 + (pkin(2) * t558 + qJ(3) * t556) * t524;
t380 = -pkin(1) * t419 - t577;
t373 = qJ(3) * t515 + t621;
t372 = t631 * t515 - t620;
t366 = -qJ(4) * t467 + t370 * t552;
t365 = qJ(4) * t470 + t370 * t550;
t362 = qJ(3) * t484 + t565;
t361 = qJ(3) * t482 + t566;
t356 = t631 * t484 + t562;
t355 = t631 * t482 + t563;
t354 = -pkin(1) * t396 - t636;
t353 = -pkin(1) * t395 - t637;
t351 = -pkin(5) * t422 + t365 * t556 + t366 * t558;
t350 = pkin(5) * t423 - t365 * t558 + t366 * t556;
t349 = t363 * t556 + t364 * t558;
t348 = -t363 * t558 + t364 * t556;
t347 = -pkin(5) * t645 - t372 * t556 + t373 * t558;
t346 = pkin(1) * t515 - pkin(5) * t357 + t372 * t558 + t373 * t556;
t345 = -pkin(5) * t396 - t356 * t556 + t362 * t558;
t344 = -pkin(5) * t395 - t355 * t556 + t361 * t558;
t343 = -pkin(1) * t645 - t638;
t342 = pkin(1) * t484 + pkin(5) * t398 + t356 * t558 + t362 * t556;
t341 = pkin(1) * t482 + pkin(5) * t397 + t355 * t558 + t361 * t556;
t340 = -(qJ(3) - t626) * t370 + t588;
t339 = -t622 - (pkin(2) - t572) * t370;
t338 = -pkin(1) * t348 - t639;
t337 = -pkin(5) * t348 - t339 * t556 + t340 * t558;
t336 = -pkin(1) * t370 + pkin(5) * t349 + t339 * t558 + t340 * t556;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t609, -t607, -t640, -qJ(1) * t640, 0, 0, t501, 0, -t608, t538, -t464 + (-t453 + t625) * t551, t551 * t589 - t647, t553 * t431, -qJ(1) * (t551 * t586 + t502) - (pkin(1) * t551 - pkin(5) * t553) * t431, 0, t501, 0, t538, t608, 0, -t464 + (-t438 + t625) * t551, -t553 * t419, -t551 * t590 + t647, t553 * t394 - t551 * t380 - qJ(1) * (t551 * t587 + t502), 0, 0, -t553 * t585, 0, -t553 * t457, t538, -t400 * t553 - t551 * t591, -t551 * t592 + t553 * t646, t553 * t645, t553 * t347 - t551 * t343 - qJ(1) * (-t357 * t551 + t515 * t553), t412 * t553 - t480 * t551, t389 * t553 - t454 * t551, t408 * t553 - t483 * t551, t411 * t553 - t479 * t551, t407 * t553 - t481 * t551, t553 * t426, t553 * t344 - t551 * t353 - qJ(1) * (t397 * t551 + t482 * t553), t553 * t345 - t551 * t354 - qJ(1) * (t398 * t551 + t484 * t553), t351 * t553 - t551 * t593, t553 * t337 - t551 * t338 - qJ(1) * (t349 * t551 - t370 * t553); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t607, -t609, t583, qJ(1) * t583, 0, 0, t500, 0, -t610, -t600, t453 * t553 + t605, -t553 * t589 - t648, t551 * t431, qJ(1) * (t553 * t586 - t617) - (-pkin(1) * t553 - pkin(5) * t551) * t431, 0, t500, 0, -t600, t610, 0, t438 * t553 + t605, -t551 * t419, t553 * t590 + t648, t551 * t394 + t553 * t380 + qJ(1) * (t553 * t587 - t617), 0, 0, -t551 * t585, 0, -t551 * t457, -t600, -t400 * t551 + t553 * t591, t551 * t646 + t553 * t592, t551 * t645, t551 * t347 + t553 * t343 + qJ(1) * (-t357 * t553 - t515 * t551), t412 * t551 + t480 * t553, t389 * t551 + t454 * t553, t408 * t551 + t483 * t553, t411 * t551 + t479 * t553, t407 * t551 + t481 * t553, t551 * t426, t551 * t344 + t553 * t353 + qJ(1) * (t397 * t553 - t482 * t551), t551 * t345 + t553 * t354 + qJ(1) * (t398 * t553 - t484 * t551), t351 * t551 + t553 * t593, t551 * t337 + t553 * t338 + qJ(1) * (t349 * t553 + t370 * t551); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t524, t525, 0, 0, 0, 0, t520, 0, t521, 0, t594, -t474, t586, t424, 0, t520, 0, 0, -t521, 0, t594, t587, t474, t387, 0, 0, -t457, 0, t585, 0, t646, t400, t357, t346, t410, t388, t406, t409, t405, t425, t341, t342, t350, t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t548, -t524, 0, 0, 0, t521, 0, -t520, 0, -t474, -t594, t431, pkin(5) * t431, 0, t521, 0, 0, t520, 0, -t474, -t419, t594, t394, 0, 0, -t585, 0, -t457, 0, -t400, t646, t645, t347, t412, t389, t408, t411, t407, t426, t344, t345, t351, t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t548, 0, -t525, 0, 0, 0, 0, 0, 0, -qJDD(2), t453, t452, 0, pkin(1) * t431, 0, 0, 0, -qJDD(2), 0, 0, t438, 0, t437, t380, 0, 0, 0, 0, 0, -qJDD(2), t379, t378, 0, t343, t480, t454, t483, t479, t481, 0, t353, t354, t352, t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t524, t525, 0, 0, 0, 0, t520, 0, t521, 0, t594, -t474, t586, t424, 0, t520, 0, 0, -t521, 0, t594, t587, t474, t387, 0, 0, -t457, 0, t585, 0, t646, t400, t357, t346, t410, t388, t406, t409, t405, t425, t341, t342, t350, t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t560, 0, 0, -t524, t494, 0, 0, qJDD(2), 0, 0, t560, 0, 0, t471, t524, qJ(3) * t524, 0, 0, -t518, 0, -t517, 0, t574, t468, t376, t373, t447, t428, t443, t446, t442, t491, t361, t362, t366, t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t560, 0, qJDD(2), 0, t524, 0, t495, 0, 0, t560, 0, 0, -qJDD(2), 0, t524, t456, 0, pkin(2) * t524, 0, 0, -t517, 0, t518, 0, t468, -t574, -t573, t372, -t445, -t427, -t441, -t444, -t440, -t490, t355, t356, -t365, t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t494, -t495, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t568, 0, t543 + t576, t577, 0, 0, 0, 0, 0, qJDD(2), t633, t634, 0, t638, -t480, -t454, -t483, -t479, -t481, 0, t637, t636, t635, t639; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t560, 0, 0, t471, t524, 0, 0, 0, -t518, 0, -t517, 0, t574, t468, t376, t621, t447, t428, t443, t446, t442, t491, t566, t565, t366, t370 * t626 + t588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, -t471, 0, t456, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(3) * t518 + t416, pkin(3) * t517 + t417, 0, pkin(3) * t376, -t480, -t454, -t483, -t479, -t481, 0, -pkin(3) * t433 + t596, -pkin(3) * t434 - t597, -pkin(3) * t467 - t580, -pkin(3) * t363 + t606; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t560, 0, 0, qJDD(2), 0, -t524, -t456, 0, 0, 0, 0, t517, 0, -t518, 0, -t468, t574, t573, -pkin(3) * t515 + t620, t445, t427, t441, t444, t440, t490, -pkin(3) * t482 - t563, -pkin(3) * t484 - t562, t365, -t370 * t572 + t622; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), 0, -t560, 0, 0, t515, t416, 0, t493, t455, t487, t492, t485, t508, t385, t386, t370, pkin(6) * t370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t560, 0, -qJDD(2), 0, -t515, 0, t417, 0, t532, -t523, t599, -t532, t598, -qJDD(5), t381, t382, 0, t628; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2), -t416, -t417, 0, 0, t480, t454, t483, t479, t481, 0, -t596, t597, t580, -t606; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t511, t512, t526, t595, t530, -t595, 0, t414, t392, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t536, -t510, t528, t513, t527, t536, -t414, 0, t393, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t532, t523, -t599, t532, -t598, qJDD(5), -t392, -t393, 0, 0;];
m_new_reg = t1;
