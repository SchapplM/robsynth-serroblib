% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPR10_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:13
% EndTime: 2019-12-31 18:26:19
% DurationCPUTime: 6.15s
% Computational Cost: add. (25026->390), mult. (35502->454), div. (0->0), fcn. (14858->8), ass. (0->229)
t589 = -qJD(1) + qJD(3);
t587 = t589 ^ 2;
t588 = qJDD(1) - qJDD(3);
t594 = sin(pkin(8));
t595 = cos(pkin(8));
t553 = t595 * t587 - t594 * t588;
t555 = t594 * t587 + t595 * t588;
t598 = sin(qJ(3));
t601 = cos(qJ(3));
t499 = t601 * t553 - t598 * t555;
t593 = g(3) + qJDD(4);
t528 = qJ(4) * t553 + t595 * t593;
t676 = qJ(4) * t555 + t594 * t593;
t465 = pkin(6) * t499 + t601 * t528 - t598 * t676;
t599 = sin(qJ(1));
t602 = cos(qJ(1));
t503 = t598 * t553 + t601 * t555;
t673 = t599 * t499 - t602 * t503;
t686 = pkin(6) * t503 + t598 * t528 + t601 * t676;
t697 = -pkin(5) * t673 + t599 * t465 - t602 * t686;
t687 = t602 * t499 + t599 * t503;
t696 = -pkin(5) * t687 + t602 * t465 + t599 * t686;
t605 = qJD(1) ^ 2;
t574 = t599 * g(1) - t602 * g(2);
t630 = -qJDD(2) + t574;
t615 = -t605 * qJ(2) - t630;
t662 = pkin(1) + pkin(2);
t608 = -t662 * qJDD(1) + t615;
t643 = (qJD(2) * qJD(1));
t585 = 2 * t643;
t575 = t602 * g(1) + t599 * g(2);
t590 = qJDD(1) * qJ(2);
t621 = -t575 + t590;
t609 = -t662 * t605 + t585 + t621;
t487 = t598 * t608 + t601 * t609;
t481 = -t587 * pkin(3) + t487;
t607 = -t598 * t609 + t601 * t608;
t606 = -t588 * pkin(3) + t607;
t444 = t594 * t481 - t595 * t606;
t445 = t595 * t481 + t594 * t606;
t639 = t594 * t444 + t595 * t445;
t426 = t595 * t444 - t594 * t445;
t647 = t601 * t426;
t404 = -t598 * t639 + t647;
t650 = t598 * t426;
t683 = t601 * t639 + t650;
t694 = t599 * t404 - t602 * t683;
t693 = t602 * t404 + t599 * t683;
t624 = t598 * t587 + t601 * t588;
t663 = pkin(6) * t624 + t598 * g(3);
t558 = t601 * t587 - t598 * t588;
t672 = t602 * t558 + t599 * t624;
t677 = pkin(6) * t558 + t601 * g(3);
t688 = -pkin(5) * t672 + t599 * t663 + t602 * t677;
t626 = t599 * t558 - t602 * t624;
t681 = -pkin(5) * t626 + t599 * t677 - t602 * t663;
t454 = -t598 * t487 - t601 * t607;
t629 = t601 * t487 - t598 * t607;
t675 = t599 * t454 - t602 * t629;
t674 = t602 * t454 + t599 * t629;
t442 = -t587 * pkin(4) - t588 * pkin(7) + t445;
t597 = sin(qJ(5));
t600 = cos(qJ(5));
t434 = t597 * t442 - t600 * t593;
t435 = t600 * t442 + t597 * t593;
t420 = t597 * t434 + t600 * t435;
t661 = pkin(3) * t426;
t660 = pkin(6) * t454;
t657 = qJDD(1) * pkin(1);
t591 = t597 ^ 2;
t656 = t591 * t587;
t441 = t588 * pkin(4) - t587 * pkin(7) + t444;
t436 = t597 * t441;
t573 = t600 * t587 * t597;
t565 = qJDD(5) + t573;
t653 = t597 * t565;
t566 = qJDD(5) - t573;
t652 = t597 * t566;
t651 = t597 * t588;
t437 = t600 * t441;
t649 = t600 * t565;
t648 = t600 * t566;
t578 = t600 * t588;
t646 = -pkin(4) * t441 + pkin(7) * t420;
t592 = t600 ^ 2;
t645 = t591 + t592;
t644 = qJD(5) * t589;
t604 = qJD(5) ^ 2;
t570 = -t604 - t656;
t520 = -t597 * t570 - t648;
t577 = t600 * t644;
t547 = 0.2e1 * t577 - t651;
t642 = -pkin(4) * t547 + pkin(7) * t520 + t436;
t580 = t592 * t587;
t572 = -t580 - t604;
t518 = t600 * t572 - t653;
t640 = t597 * t644;
t550 = -t578 - 0.2e1 * t640;
t641 = pkin(4) * t550 + pkin(7) * t518 - t437;
t537 = t605 * pkin(1) - t621 - (2 * t643);
t541 = -t615 + t657;
t638 = -t602 * t537 - t599 * t541;
t637 = -t599 * t574 - t602 * t575;
t636 = t594 * t573;
t635 = t595 * t573;
t556 = t645 * t588;
t561 = t580 + t656;
t634 = pkin(4) * t561 - pkin(7) * t556 + t420;
t633 = -pkin(2) * g(3) + pkin(6) * t629;
t408 = t594 * t420 - t595 * t441;
t632 = pkin(3) * t408 + t646;
t631 = -pkin(3) * t553 - t445;
t567 = t599 * qJDD(1) + t602 * t605;
t544 = -pkin(5) * t567 + t602 * g(3);
t568 = t602 * qJDD(1) - t599 * t605;
t543 = pkin(5) * t568 + t599 * g(3);
t419 = t600 * t434 - t597 * t435;
t628 = t599 * t537 - t602 * t541;
t625 = t602 * t574 - t599 * t575;
t483 = t594 * t520 - t595 * t547;
t623 = pkin(3) * t483 + t642;
t482 = t594 * t518 + t595 * t550;
t622 = pkin(3) * t482 + t641;
t620 = -pkin(3) * t555 - t444;
t409 = t595 * t420 + t594 * t441;
t386 = qJ(4) * t409 - (-pkin(4) * t595 - pkin(7) * t594 - pkin(3)) * t419;
t390 = -qJ(4) * t408 - (pkin(4) * t594 - pkin(7) * t595) * t419;
t395 = t601 * t408 + t598 * t409;
t619 = -pkin(6) * t395 - t598 * t386 + t601 * t390;
t514 = t597 * t572 + t649;
t428 = -pkin(4) * t514 + t434;
t430 = -pkin(7) * t514 + t436;
t484 = t595 * t518 - t594 * t550;
t406 = -pkin(3) * t514 + qJ(4) * t484 + t595 * t428 + t594 * t430;
t412 = -qJ(4) * t482 - t594 * t428 + t595 * t430;
t447 = t601 * t482 + t598 * t484;
t618 = -pkin(6) * t447 - t598 * t406 + t601 * t412;
t516 = t600 * t570 - t652;
t429 = -pkin(4) * t516 + t435;
t431 = -pkin(7) * t516 + t437;
t485 = t595 * t520 + t594 * t547;
t407 = -pkin(3) * t516 + qJ(4) * t485 + t595 * t429 + t594 * t431;
t413 = -qJ(4) * t483 - t594 * t429 + t595 * t431;
t448 = t601 * t483 + t598 * t485;
t617 = -pkin(6) * t448 - t598 * t407 + t601 * t413;
t505 = -t594 * t556 + t595 * t561;
t616 = pkin(3) * t505 + t634;
t396 = -t598 * t408 + t601 * t409;
t614 = -pkin(6) * t396 - t601 * t386 - t598 * t390;
t449 = -t598 * t482 + t601 * t484;
t613 = -pkin(6) * t449 - t601 * t406 - t598 * t412;
t450 = -t598 * t483 + t601 * t485;
t612 = -pkin(6) * t450 - t601 * t407 - t598 * t413;
t423 = -pkin(3) * t593 + qJ(4) * t639;
t611 = pkin(6) * t404 + qJ(4) * t647 - t598 * t423;
t610 = -pkin(6) * t683 - qJ(4) * t650 - t601 * t423;
t603 = pkin(1) * g(3);
t596 = qJ(2) * g(3);
t571 = t580 - t604;
t569 = t604 - t656;
t563 = t630 + 0.2e1 * t657;
t562 = -t580 + t656;
t549 = -t578 - t640;
t548 = t577 - t651;
t546 = -t575 + t585 + 0.2e1 * t590;
t542 = t645 * t644;
t524 = t594 * qJDD(5) + t595 * t542;
t523 = -t595 * qJDD(5) + t594 * t542;
t522 = t600 * t548 - t591 * t644;
t521 = -t597 * t549 - t592 * t644;
t519 = -t597 * t569 + t649;
t517 = t600 * t571 - t652;
t515 = t600 * t569 + t653;
t513 = t597 * t571 + t648;
t510 = (t548 + t577) * t597;
t509 = -t600 * t549 + t597 * t577;
t506 = -t595 * t556 - t594 * t561;
t498 = -t597 * t547 + t600 * t550;
t497 = t600 * t547 + t597 * t550;
t496 = pkin(1) * t541 - qJ(2) * t537;
t495 = t595 * t519 - t594 * t651;
t494 = t595 * t517 - t594 * t578;
t493 = t594 * t519 + t595 * t651;
t492 = t594 * t517 + t595 * t578;
t491 = t595 * t522 - t636;
t490 = t595 * t521 + t636;
t489 = t594 * t522 + t635;
t488 = t594 * t521 - t635;
t477 = -t598 * t523 + t601 * t524;
t476 = t601 * t523 + t598 * t524;
t475 = t595 * t498 + t594 * t562;
t474 = t594 * t498 - t595 * t562;
t473 = -t598 * t505 + t601 * t506;
t472 = t601 * t505 + t598 * t506;
t469 = -t598 * t493 + t601 * t495;
t468 = -t598 * t492 + t601 * t494;
t467 = t601 * t493 + t598 * t495;
t466 = t601 * t492 + t598 * t494;
t461 = -t598 * t489 + t601 * t491;
t460 = -t598 * t488 + t601 * t490;
t459 = t601 * t489 + t598 * t491;
t458 = t601 * t488 + t598 * t490;
t457 = -qJ(2) * t558 + t624 * t662 - t607;
t456 = qJ(2) * t624 + t558 * t662 + t487;
t451 = t596 + t660;
t446 = t603 - t633;
t439 = -t598 * t474 + t601 * t475;
t438 = t601 * t474 + t598 * t475;
t422 = -qJ(2) * t499 + t503 * t662 - t620;
t421 = qJ(2) * t503 + t499 * t662 - t631;
t416 = qJ(2) * t629 + t454 * t662;
t415 = -qJ(4) * t505 + t595 * t419;
t414 = qJ(4) * t506 + t594 * t419;
t401 = qJ(2) * t450 - t662 * t448 - t623;
t400 = qJ(2) * t449 - t662 * t447 - t622;
t399 = qJ(2) * t473 - t662 * t472 - t616;
t398 = -pkin(6) * t472 - t598 * t414 + t601 * t415;
t397 = pkin(6) * t473 + t601 * t414 + t598 * t415;
t394 = qJ(2) * t516 + t617;
t393 = qJ(2) * t514 + t618;
t392 = t662 * t516 + t612;
t391 = t662 * t514 + t613;
t389 = qJ(2) * t593 + t611;
t387 = t662 * t593 + t610;
t385 = qJ(2) * t683 + t404 * t662 + t661;
t384 = qJ(2) * t396 - t662 * t395 - t632;
t383 = -qJ(2) * t419 + t619;
t382 = -t419 * t662 + t614;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t568, 0, -t567, 0, -t543, -t544, -t625, -pkin(5) * t625, 0, t568, 0, 0, t567, 0, -t543, t628, t544, pkin(5) * t628 + (-t599 * pkin(1) + t602 * qJ(2)) * g(3), 0, 0, t626, 0, -t672, 0, -t681, t688, t674, -pkin(5) * t674 - t599 * t446 + t602 * t451, 0, 0, t673, 0, -t687, 0, -t697, t696, t693, -pkin(5) * t693 - t599 * t387 + t602 * t389, t599 * t459 + t602 * t461, t599 * t438 + t602 * t439, t599 * t467 + t602 * t469, t599 * t458 + t602 * t460, t599 * t466 + t602 * t468, t599 * t476 + t602 * t477, t602 * t393 - t599 * t391 - pkin(5) * (-t602 * t447 + t599 * t449), t602 * t394 - t599 * t392 - pkin(5) * (-t602 * t448 + t599 * t450), t602 * t398 + t599 * t397 - pkin(5) * (-t602 * t472 + t599 * t473), t602 * t383 - t599 * t382 - pkin(5) * (-t602 * t395 + t599 * t396); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t567, 0, t568, 0, t544, -t543, t637, pkin(5) * t637, 0, t567, 0, 0, -t568, 0, t544, t638, t543, pkin(5) * t638 + (t602 * pkin(1) + t599 * qJ(2)) * g(3), 0, 0, -t672, 0, -t626, 0, t688, t681, t675, -pkin(5) * t675 + t602 * t446 + t599 * t451, 0, 0, -t687, 0, -t673, 0, t696, t697, t694, -pkin(5) * t694 + t602 * t387 + t599 * t389, -t602 * t459 + t599 * t461, -t602 * t438 + t599 * t439, -t602 * t467 + t599 * t469, -t602 * t458 + t599 * t460, -t602 * t466 + t599 * t468, -t602 * t476 + t599 * t477, t599 * t393 + t602 * t391 + pkin(5) * (t599 * t447 + t602 * t449), t599 * t394 + t602 * t392 + pkin(5) * (t599 * t448 + t602 * t450), t599 * t398 - t602 * t397 + pkin(5) * (t599 * t472 + t602 * t473), t599 * t383 + t602 * t382 + pkin(5) * (t599 * t395 + t602 * t396); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t574, t575, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t563, 0, t546, t496, 0, 0, 0, 0, 0, t588, t457, t456, 0, t416, 0, 0, 0, 0, 0, t588, t422, t421, 0, t385, -t510, -t497, -t515, t509, -t513, 0, t400, t401, t399, t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t605, 0, 0, -g(3), -t574, 0, 0, qJDD(1), 0, 0, t605, 0, 0, -t541, g(3), t596, 0, 0, -t624, 0, -t558, 0, t663, t677, t454, t451, 0, 0, -t503, 0, -t499, 0, t686, t465, t404, t389, t461, t439, t469, t460, t468, t477, t393, t394, t398, t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, 0, qJDD(1), 0, g(3), 0, -t575, 0, 0, t605, 0, 0, -qJDD(1), 0, g(3), -t537, 0, t603, 0, 0, -t558, 0, t624, 0, t677, -t663, -t629, t446, 0, 0, -t499, 0, t503, 0, t465, -t686, -t683, t387, -t459, -t438, -t467, -t458, -t466, -t476, t391, t392, -t397, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t574, t575, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t563, 0, t546, t496, 0, 0, 0, 0, 0, t588, t457, t456, 0, t416, 0, 0, 0, 0, 0, t588, t422, t421, 0, t385, -t510, -t497, -t515, t509, -t513, 0, t400, t401, t399, t384; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t605, 0, 0, -t541, g(3), 0, 0, 0, -t624, 0, -t558, 0, t663, t677, t454, t660, 0, 0, -t503, 0, -t499, 0, t686, t465, t404, t611, t461, t439, t469, t460, t468, t477, t618, t617, t398, t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, t541, 0, -t537, 0, 0, 0, 0, 0, 0, t588, pkin(2) * t624 - t607, pkin(2) * t558 + t487, 0, pkin(2) * t454, 0, 0, 0, 0, 0, t588, pkin(2) * t503 - t620, pkin(2) * t499 - t631, 0, pkin(2) * t404 + t661, -t510, -t497, -t515, t509, -t513, 0, -pkin(2) * t447 - t622, -pkin(2) * t448 - t623, -pkin(2) * t472 - t616, -pkin(2) * t395 - t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t605, 0, 0, qJDD(1), 0, -g(3), t537, 0, 0, 0, 0, t558, 0, -t624, 0, -t677, t663, t629, t633, 0, 0, t499, 0, -t503, 0, -t465, t686, t683, -pkin(2) * t593 - t610, t459, t438, t467, t458, t466, t476, -pkin(2) * t514 - t613, -pkin(2) * t516 - t612, t397, pkin(2) * t419 - t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, 0, -t587, 0, 0, g(3), -t607, 0, 0, 0, -t555, 0, -t553, 0, t676, t528, t426, qJ(4) * t426, t491, t475, t495, t490, t494, t524, t412, t413, t415, t390; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t587, 0, -t588, 0, -g(3), 0, t487, 0, 0, 0, t553, 0, -t555, 0, -t528, t676, t639, t423, t489, t474, t493, t488, t492, t523, t406, t407, t414, t386; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, t607, -t487, 0, 0, 0, 0, 0, 0, 0, -t588, t620, t631, 0, -t661, t510, t497, t515, -t509, t513, 0, t622, t623, t616, t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, 0, -t587, 0, 0, t593, t444, 0, t522, t498, t519, t521, t517, t542, t430, t431, t419, pkin(7) * t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t587, 0, -t588, 0, -t593, 0, t445, 0, t573, -t562, t651, -t573, t578, -qJDD(5), t428, t429, 0, pkin(4) * t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t588, -t444, -t445, 0, 0, t510, t497, t515, -t509, t513, 0, t641, t642, t634, t646; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t548, t550, t565, -t577, t571, t577, 0, t441, t434, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t640, t547, t569, t549, t566, -t640, -t441, 0, t435, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t573, t562, -t651, t573, -t578, qJDD(5), -t434, -t435, 0, 0;];
m_new_reg = t1;