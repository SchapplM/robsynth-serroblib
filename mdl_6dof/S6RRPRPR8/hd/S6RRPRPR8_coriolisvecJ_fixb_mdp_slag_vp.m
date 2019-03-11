% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:32
% EndTime: 2019-03-09 10:53:43
% DurationCPUTime: 7.30s
% Computational Cost: add. (5610->523), mult. (13915->685), div. (0->0), fcn. (10115->8), ass. (0->228)
t559 = sin(pkin(10));
t562 = sin(qJ(4));
t560 = cos(pkin(10));
t681 = cos(qJ(4));
t629 = t681 * t560;
t585 = -t562 * t559 + t629;
t565 = cos(qJ(2));
t652 = qJD(1) * t565;
t624 = qJD(4) * t681;
t645 = qJD(4) * t562;
t688 = -t559 * t645 + t560 * t624;
t657 = -t585 * t652 + t688;
t563 = sin(qJ(2));
t653 = qJD(1) * t563;
t625 = t559 * t653;
t648 = qJD(2) * t560;
t512 = -t625 + t648;
t628 = t560 * t653;
t649 = qJD(2) * t559;
t513 = t628 + t649;
t458 = -t681 * t512 + t513 * t562;
t561 = sin(qJ(6));
t564 = cos(qJ(6));
t586 = -t562 * t512 - t513 * t681;
t394 = t458 * t561 - t564 * t586;
t593 = -t564 * t458 - t561 * t586;
t646 = qJD(2) * t565;
t590 = t629 * t646;
t635 = qJD(1) * qJD(2);
t622 = t565 * t635;
t608 = t559 * t622;
t410 = t562 * (qJD(4) * t513 + t608) - qJD(1) * t590 - t512 * t624;
t518 = t559 * t681 + t562 * t560;
t574 = t565 * t518;
t573 = qJD(2) * t574;
t687 = qJD(4) * t586;
t411 = qJD(1) * t573 - t687;
t643 = qJD(6) * t564;
t644 = qJD(6) * t561;
t611 = t564 * t410 - t561 * t411 - t458 * t643 - t586 * t644;
t641 = t563 * MDP(30);
t620 = qJD(2) * t641;
t546 = -qJD(4) + t652;
t638 = -qJD(6) - t546;
t699 = t593 * t638;
t703 = (t394 ^ 2 - t593 ^ 2) * MDP(27) + MDP(26) * t394 * t593 - (t611 + t699) * MDP(28) - qJD(1) * t620;
t702 = t458 ^ 2;
t683 = t586 ^ 2;
t701 = t394 * t638;
t674 = t458 * t546;
t700 = t546 * t586;
t598 = pkin(2) * t563 - qJ(3) * t565;
t520 = t598 * qJD(1);
t478 = pkin(7) * t625 + t560 * t520;
t670 = t560 * t565;
t588 = pkin(3) * t563 - pkin(8) * t670;
t446 = qJD(1) * t588 + t478;
t503 = t559 * t520;
t671 = t560 * t563;
t672 = t559 * t565;
t584 = -pkin(7) * t671 - pkin(8) * t672;
t465 = qJD(1) * t584 + t503;
t680 = pkin(8) + qJ(3);
t528 = t680 * t559;
t529 = t680 * t560;
t698 = qJD(3) * t629 - t681 * t465 - t528 * t624 + (-qJD(3) * t559 - qJD(4) * t529 - t446) * t562;
t502 = t518 * qJD(4);
t656 = -qJD(1) * t574 + t502;
t552 = pkin(7) * t653;
t679 = qJD(2) * pkin(2);
t618 = -qJD(3) + t679;
t605 = -t552 + t618;
t579 = pkin(3) * t512 + t605;
t571 = -qJ(5) * t586 + t579;
t682 = pkin(4) + pkin(5);
t374 = -t458 * t682 + t571;
t532 = t546 * qJD(5);
t623 = t563 * t635;
t542 = qJ(5) * t623;
t498 = qJD(2) * t598 - qJD(3) * t563;
t488 = t498 * qJD(1);
t523 = (qJD(3) - t552) * qJD(2);
t444 = t560 * t488 - t523 * t559;
t578 = t588 * qJD(2);
t414 = qJD(1) * t578 + t444;
t525 = -pkin(2) * t565 - qJ(3) * t563 - pkin(1);
t505 = t525 * qJD(1);
t553 = pkin(7) * t652;
t531 = qJD(2) * qJ(3) + t553;
t468 = t560 * t505 - t531 * t559;
t633 = pkin(3) * t652;
t419 = -pkin(8) * t513 + t468 - t633;
t445 = t559 * t488 + t560 * t523;
t426 = -pkin(8) * t608 + t445;
t469 = t559 * t505 + t560 * t531;
t429 = pkin(8) * t512 + t469;
t581 = -t562 * t414 - t419 * t624 - t681 * t426 + t429 * t645;
t356 = -t532 + t542 - t581;
t354 = pkin(9) * t411 + t356;
t631 = t682 * t563;
t610 = qJD(2) * t631;
t612 = -t681 * t414 + t419 * t645 + t562 * t426 + t429 * t624;
t355 = pkin(9) * t410 - qJD(1) * t610 + t612;
t617 = t561 * t354 - t564 * t355;
t696 = t374 * t394 + t617;
t693 = -0.2e1 * t635;
t692 = MDP(4) * t563;
t691 = MDP(5) * (t563 ^ 2 - t565 ^ 2);
t666 = qJ(5) * t653 - t698;
t476 = -t562 * t528 + t681 * t529;
t690 = t518 * qJD(3) + qJD(4) * t476 + t446 * t681 - t562 * t465;
t506 = t559 * t633 + t553;
t689 = qJ(5) * t657 + qJD(5) * t518 + t506;
t380 = t681 * t419 - t562 * t429;
t639 = qJD(5) - t380;
t640 = -pkin(9) * t586 - t639;
t368 = t546 * t682 - t640;
t630 = t564 * t354 + t561 * t355 + t368 * t643;
t686 = -t374 * t593 + t630;
t616 = t561 * t410 + t564 * t411;
t362 = t394 * qJD(6) - t616;
t678 = qJ(5) * t458;
t381 = t562 * t419 + t681 * t429;
t371 = pkin(9) * t458 + t381;
t533 = t546 * qJ(5);
t369 = t371 - t533;
t677 = t369 * t561;
t385 = pkin(4) * t458 - t571;
t676 = t385 * t586;
t675 = t458 * t586;
t673 = t559 * t563;
t567 = qJD(2) ^ 2;
t669 = t563 * t567;
t668 = t565 * t567;
t568 = qJD(1) ^ 2;
t667 = t565 * t568;
t665 = pkin(4) * t653 + t690;
t591 = -t518 * t561 - t564 * t585;
t664 = qJD(6) * t591 + t561 * t656 + t564 * t657;
t464 = t518 * t564 - t561 * t585;
t663 = qJD(6) * t464 + t561 * t657 - t564 * t656;
t662 = -t656 * t682 + t689;
t661 = -pkin(4) * t656 + t689;
t511 = t560 * t525;
t467 = -pkin(8) * t671 + t511 + (-pkin(7) * t559 - pkin(3)) * t565;
t544 = pkin(7) * t670;
t484 = t559 * t525 + t544;
t474 = -pkin(8) * t673 + t484;
t659 = t562 * t467 + t681 * t474;
t647 = qJD(2) * t563;
t632 = pkin(7) * t647;
t472 = t560 * t498 + t559 * t632;
t545 = pkin(7) * t622;
t497 = pkin(3) * t608 + t545;
t554 = pkin(7) * t646;
t627 = t559 * t646;
t507 = pkin(3) * t627 + t554;
t521 = pkin(3) * t673 + t563 * pkin(7);
t475 = t528 * t681 + t562 * t529;
t651 = qJD(2) * t475;
t650 = qJD(2) * t476;
t642 = t563 * MDP(19);
t637 = MDP(11) * qJD(1);
t636 = MDP(12) * qJD(1);
t634 = pkin(7) * t672;
t549 = -pkin(3) * t560 - pkin(2);
t621 = qJD(2) * t642;
t619 = pkin(1) * t693;
t615 = -t512 + t648;
t614 = -t513 + t649;
t613 = pkin(4) * t623;
t606 = t605 - t618;
t397 = -qJ(5) * t565 + t659;
t495 = t585 * t563;
t603 = qJ(5) * t495 - t521;
t601 = t467 * t681 - t562 * t474;
t439 = -pkin(9) * t585 + t476;
t600 = pkin(9) * t657 - qJD(1) * t631 + qJD(6) * t439 - t690;
t438 = -t518 * pkin(9) + t475;
t599 = -pkin(9) * t656 - qJD(6) * t438 + t666;
t597 = qJ(5) * t564 - t561 * t682;
t596 = qJ(5) * t561 + t564 * t682;
t350 = t368 * t561 + t369 * t564;
t398 = t565 * pkin(4) - t601;
t384 = t565 * pkin(5) - t495 * pkin(9) + t398;
t494 = t518 * t563;
t386 = pkin(9) * t494 + t397;
t595 = t384 * t564 - t386 * t561;
t594 = t384 * t561 + t386 * t564;
t592 = t564 * t494 - t495 * t561;
t435 = t494 * t561 + t495 * t564;
t587 = qJ(5) * t518 - t549;
t583 = -t381 * t546 - t612;
t436 = t578 + t472;
t489 = t559 * t498;
t449 = qJD(2) * t584 + t489;
t582 = -t436 * t681 + t562 * t449 + t467 * t645 + t474 * t624;
t580 = t562 * t436 + t681 * t449 + t467 * t624 - t474 * t645;
t577 = -qJ(5) * t410 - qJD(5) * t586 - t497;
t440 = t502 * t563 + t562 * t627 - t590;
t576 = -qJ(5) * t440 + qJD(5) * t495 - t507;
t358 = t612 - t613;
t572 = -t380 * t546 + t581;
t363 = pkin(4) * t411 - t577;
t364 = qJ(5) * t647 - qJD(5) * t565 + t580;
t483 = t511 - t634;
t479 = -pkin(7) * t628 + t503;
t473 = -t560 * t632 + t489;
t456 = -pkin(4) * t585 - t587;
t441 = t563 * t688 + t573;
t425 = t585 * t682 + t587;
t420 = pkin(4) * t494 - t603;
t400 = -pkin(4) * t586 + t678;
t399 = -t494 * t682 + t603;
t383 = -t410 - t674;
t382 = t586 * t682 - t678;
t377 = -t533 + t381;
t376 = pkin(4) * t546 + t639;
t375 = pkin(4) * t441 - t576;
t373 = qJD(6) * t435 - t440 * t561 - t564 * t441;
t372 = qJD(6) * t592 - t440 * t564 + t441 * t561;
t367 = -t441 * t682 + t576;
t366 = -pkin(4) * t647 + t582;
t360 = pkin(9) * t441 + t364;
t359 = t440 * pkin(9) + t582 - t610;
t357 = -t411 * t682 + t577;
t349 = t368 * t564 - t677;
t1 = [(-t546 - t652) * t621 + (-t472 * t513 + t473 * t512 + (-t444 * t560 - t445 * t559) * t563 + (-t468 * t560 - t469 * t559 + (-t483 * t560 - t484 * t559) * qJD(1)) * t646) * MDP(13) + (t411 * t565 + t441 * t546 + (-qJD(1) * t494 - t458) * t647) * MDP(18) + (t358 * t565 + t363 * t494 + t366 * t546 + t375 * t458 + t385 * t441 + t411 * t420 + (-qJD(1) * t398 - t376) * t647) * MDP(22) + 0.2e1 * t622 * t692 + t691 * t693 - MDP(7) * t669 + (pkin(7) * t669 + t565 * t619) * MDP(10) + MDP(6) * t668 + (t356 * t397 + t358 * t398 + t363 * t420 + t364 * t377 + t366 * t376 + t375 * t385) * MDP(25) + (-pkin(7) * t668 + t563 * t619) * MDP(9) + (t580 * t546 - t581 * t565 - t507 * t586 - t521 * t410 + t497 * t495 + t579 * t440 + (-qJD(1) * t659 - t381) * t647) * MDP(21) + (t582 * t546 + t612 * t565 + t507 * t458 + t521 * t411 + t497 * t494 - t579 * t441 + (qJD(1) * t601 + t380) * t647) * MDP(20) + ((qJD(1) * t473 + t445) * t565 + ((pkin(7) * t513 - t560 * t605) * t565 + (-t469 + (-t484 + 0.2e1 * t544) * qJD(1)) * t563) * qJD(2)) * MDP(12) + ((-qJD(1) * t472 - t444) * t565 + ((-pkin(7) * t512 - t559 * t605) * t565 + (t468 + (t483 + 0.2e1 * t634) * qJD(1)) * t563) * qJD(2)) * MDP(11) + (t444 * t483 + t445 * t484 + t468 * t472 + t469 * t473 + (-t605 + t552) * t554) * MDP(14) + (-t611 * t565 - t372 * t638 + (-qJD(1) * t435 - t394) * t647) * MDP(28) + ((qJD(6) * t595 + t359 * t561 + t360 * t564) * t638 - (-t369 * t644 + t630) * t565 + t367 * t394 - t399 * t611 + t357 * t435 + t374 * t372 + (qJD(1) * t594 + t350) * t647) * MDP(32) + (-t362 * t565 + t373 * t638 + (-qJD(1) * t592 + t593) * t647) * MDP(29) + (-(t359 * t564 - t360 * t561) * t638 - t617 * t565 + t367 * t593 + t399 * t362 - t357 * t592 + t374 * t373 + (-t350 * t565 + t594 * t638) * qJD(6) + (-qJD(1) * t595 - t349) * t647) * MDP(31) + (t638 - t652) * t620 + (t372 * t394 - t435 * t611) * MDP(26) + (-t362 * t435 - t372 * t593 - t373 * t394 - t592 * t611) * MDP(27) + (t410 * t565 + t440 * t546 + (qJD(1) * t495 - t586) * t647) * MDP(17) + (-t356 * t565 - t363 * t495 - t364 * t546 + t375 * t586 + t385 * t440 + t410 * t420 + (qJD(1) * t397 + t377) * t647) * MDP(24) + (t410 * t494 - t411 * t495 + t440 * t458 + t441 * t586) * MDP(16) + (-t356 * t494 + t358 * t495 - t364 * t458 - t366 * t586 - t376 * t440 - t377 * t441 - t397 * t411 - t398 * t410) * MDP(23) + (-t410 * t495 + t440 * t586) * MDP(15); (t478 * t513 - t479 * t512 + (qJD(3) * t512 + t468 * t652 + t445) * t560 + (qJD(3) * t513 + t469 * t652 - t444) * t559) * MDP(13) + ((-qJ(3) * t648 + t469) * t563 + (pkin(7) * t614 + t560 * t606 - t479) * t565) * t636 + ((-qJ(3) * t649 - t468) * t563 + (-pkin(7) * t615 + t559 * t606 + t478) * t565) * t637 + (-t468 * t478 - t469 * t479 + (-t468 * t559 + t469 * t560) * qJD(3) + (-t444 * t559 + t445 * t560) * qJ(3) + (t605 - t679) * t553) * MDP(14) - t667 * t692 + t568 * t691 + (t663 * t638 + (-qJD(2) * t591 - t593) * t653) * MDP(29) + (t394 * t664 - t464 * t611) * MDP(26) + (-t664 * t638 + (-qJD(2) * t464 + t394) * t653) * MDP(28) + (t357 * t464 - t425 * t611 - (t561 * t600 + t564 * t599) * t638 + t662 * t394 + t664 * t374 + ((t438 * t561 + t439 * t564) * qJD(2) - t350) * t653) * MDP(32) + (-t362 * t464 - t394 * t663 - t591 * t611 - t593 * t664) * MDP(27) + (-t363 * t585 + t411 * t456 + t665 * t546 - t661 * t458 + t656 * t385 + (t376 - t651) * t653) * MDP(22) + (t356 * t585 + t358 * t518 + t376 * t657 - t377 * t656 - t410 * t475 - t411 * t476 + t458 * t666 - t586 * t665) * MDP(23) + (-t363 * t518 + t410 * t456 + t666 * t546 - t661 * t586 - t657 * t385 + (-t377 + t650) * t653) * MDP(24) + (t356 * t476 + t358 * t475 + t363 * t456 + t376 * t665 - t377 * t666 - t385 * t661) * MDP(25) + (t656 * t546 + (qJD(2) * t585 + t458) * t653) * MDP(18) + (t549 * t411 - t506 * t458 - t497 * t585 + t690 * t546 - t656 * t579 + (-t380 - t651) * t653) * MDP(20) + (-t410 * t518 - t586 * t657) * MDP(15) + (-t657 * t546 + (qJD(2) * t518 + t586) * t653) * MDP(17) + (-t410 * t585 - t411 * t518 - t458 * t657 + t586 * t656) * MDP(16) + (-t549 * t410 + t506 * t586 + t497 * t518 + t698 * t546 - t657 * t579 + (t381 - t650) * t653) * MDP(21) + (-t357 * t591 + t425 * t362 - (t561 * t599 - t564 * t600) * t638 + t662 * t593 + t663 * t374 + (-(t438 * t564 - t439 * t561) * qJD(2) + t349) * t653) * MDP(31) + (t546 * t642 - t638 * t641) * qJD(1) + (MDP(9) * t563 * t568 + MDP(10) * t667) * pkin(1); (-t512 ^ 2 - t513 ^ 2) * MDP(13) + (t468 * t513 - t469 * t512 + t545) * MDP(14) + (-t683 - t702) * MDP(23) + (t376 * t586 + t377 * t458 + t363) * MDP(25) + (-t362 + t701) * MDP(31) + (t611 - t699) * MDP(32) + (MDP(20) + MDP(22)) * (t411 + t700) + (t614 * t637 + t615 * t636) * t565 + (-MDP(21) + MDP(24)) * (t410 - t674); -MDP(15) * t675 + (t683 - t702) * MDP(16) + t383 * MDP(17) + (-t518 * t622 + t687 + t700) * MDP(18) + qJD(1) * t621 + (-t579 * t586 + t583) * MDP(20) + (-t458 * t579 + t572) * MDP(21) + (-t400 * t458 + t583 + 0.2e1 * t613 + t676) * MDP(22) + (pkin(4) * t410 - qJ(5) * t411 - (t377 - t381) * t586 + (t376 - t639) * t458) * MDP(23) + (-t385 * t458 - t400 * t586 - 0.2e1 * t532 + 0.2e1 * t542 - t572) * MDP(24) + (-pkin(4) * t358 + qJ(5) * t356 - t376 * t381 + t377 * t639 - t385 * t400) * MDP(25) + (t362 + t701) * MDP(29) + (t596 * t623 - t382 * t593 - (-t371 * t564 + t561 * t640) * t638 + (t597 * t638 + t350) * qJD(6) + t696) * MDP(31) + (t597 * t623 - t382 * t394 - (t371 * t561 + t564 * t640) * t638 + (-t596 * t638 - t677) * qJD(6) + t686) * MDP(32) - t703; (-t623 - t675) * MDP(22) + t383 * MDP(23) + (-t546 ^ 2 - t683) * MDP(24) + (t377 * t546 + t358 - t676) * MDP(25) + (-t564 * t623 + t586 * t593) * MDP(31) + (t394 * t586 + t561 * t623) * MDP(32) - (MDP(31) * t561 + MDP(32) * t564) * t638 ^ 2; (t616 - t701) * MDP(29) + (-t350 * t638 - t696) * MDP(31) + (-t349 * t638 - t686) * MDP(32) + (-MDP(29) * t394 - MDP(31) * t350 + MDP(32) * t677) * qJD(6) + t703;];
tauc  = t1;
