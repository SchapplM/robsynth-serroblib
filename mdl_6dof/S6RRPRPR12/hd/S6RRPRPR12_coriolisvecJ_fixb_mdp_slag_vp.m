% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:22:01
% EndTime: 2019-03-09 11:22:14
% DurationCPUTime: 7.50s
% Computational Cost: add. (7217->522), mult. (18561->718), div. (0->0), fcn. (14017->10), ass. (0->235)
t565 = cos(qJ(4));
t566 = cos(qJ(2));
t559 = sin(pkin(6));
t563 = sin(qJ(2));
t650 = qJD(1) * t563;
t631 = t559 * t650;
t536 = pkin(2) * t631;
t597 = pkin(9) * t563 - qJ(3) * t566;
t652 = qJD(1) * t559;
t476 = t597 * t652 + t536;
t630 = t566 * t652;
t560 = cos(pkin(6));
t651 = qJD(1) * t560;
t638 = pkin(1) * t651;
t497 = pkin(8) * t630 + t563 * t638;
t478 = pkin(3) * t630 + t497;
t562 = sin(qJ(4));
t612 = -t476 * t562 + t565 * t478;
t646 = qJD(4) * t562;
t567 = -pkin(2) - pkin(9);
t662 = qJ(5) - t567;
t664 = t562 * t563;
t696 = -(pkin(4) * t566 - qJ(5) * t664) * t652 - t612 - qJD(5) * t565 + t646 * t662;
t603 = t565 * t631;
t618 = t662 * t565;
t657 = t565 * t476 + t562 * t478;
t695 = qJ(5) * t603 + qJD(4) * t618 + qJD(5) * t562 + t657;
t524 = qJD(4) + t631;
t564 = cos(qJ(6));
t561 = sin(qJ(6));
t540 = qJD(2) + t651;
t486 = t540 * t562 + t565 * t630;
t488 = t540 * t565 - t562 * t630;
t558 = sin(pkin(11));
t680 = cos(pkin(11));
t582 = -t558 * t486 + t488 * t680;
t675 = t582 * t561;
t414 = -t564 * t524 + t675;
t611 = -t680 * t486 - t488 * t558;
t686 = qJD(6) - t611;
t694 = t414 * t686;
t416 = t524 * t561 + t564 * t582;
t693 = t416 * t686;
t619 = qJD(4) * t680;
t668 = t558 * t562;
t656 = -t558 * t646 + t565 * t619 + t603 * t680 - t631 * t668;
t581 = t558 * t565 + t562 * t680;
t645 = qJD(4) * t565;
t655 = t558 * t645 + t562 * t619 + t581 * t631;
t610 = t686 * t564;
t639 = qJD(1) * qJD(2);
t621 = t559 * t639;
t602 = t563 * t621;
t442 = -qJD(4) * t486 + t562 * t602;
t649 = qJD(2) * t563;
t576 = t565 * t649 + t566 * t646;
t443 = t540 * t645 - t576 * t652;
t409 = t442 * t558 + t443 * t680;
t665 = t561 * t409;
t692 = -t686 * t610 - t665;
t496 = pkin(8) * t631 - t566 * t638;
t687 = qJD(3) + t496;
t691 = t687 + pkin(4) * t645 + t631 * (pkin(4) * t565 + pkin(3));
t556 = t563 ^ 2;
t690 = MDP(5) * (-t566 ^ 2 + t556);
t661 = t695 * t558 + t696 * t680;
t659 = t696 * t558 - t695 * t680;
t667 = t559 * t563;
t542 = pkin(8) * t667;
t633 = -pkin(1) * t566 - pkin(2);
t456 = pkin(3) * t667 + t542 + (-pkin(9) + t633) * t560;
t620 = -qJ(3) * t563 - pkin(1);
t474 = (t566 * t567 + t620) * t559;
t658 = t562 * t456 + t565 * t474;
t681 = pkin(1) * t563;
t545 = t560 * t681;
t666 = t559 * t566;
t688 = pkin(8) * t666 + t545;
t642 = pkin(3) * t631 + t687;
t526 = t566 * t621;
t429 = t540 * t567 + t642;
t452 = qJD(1) * t474;
t403 = t429 * t562 + t452 * t565;
t519 = pkin(2) * t602;
t647 = qJD(3) * t563;
t569 = (qJD(2) * t597 - t647) * t559;
t436 = qJD(1) * t569 + t519;
t637 = pkin(1) * qJD(2) * t560;
t605 = qJD(1) * t637;
t489 = pkin(8) * t526 + t563 * t605;
t457 = pkin(3) * t526 + t489;
t570 = -qJD(4) * t403 - t562 * t436 + t565 * t457;
t364 = pkin(4) * t526 - qJ(5) * t442 - qJD(5) * t488 + t570;
t580 = t429 * t645 + t565 * t436 - t452 * t646 + t562 * t457;
t368 = -qJ(5) * t443 - qJD(5) * t486 + t580;
t351 = t364 * t680 - t558 * t368;
t349 = -pkin(5) * t526 - t351;
t549 = pkin(4) * t558 + pkin(10);
t685 = t686 * (pkin(4) * t488 + pkin(5) * t582 - pkin(10) * t611 + qJD(6) * t549) + t349;
t352 = t558 * t364 + t680 * t368;
t402 = t565 * t429 - t452 * t562;
t395 = -qJ(5) * t488 + t402;
t392 = pkin(4) * t524 + t395;
t396 = -qJ(5) * t486 + t403;
t669 = t558 * t396;
t365 = t392 * t680 - t669;
t393 = t680 * t396;
t366 = t558 * t392 + t393;
t505 = t565 * t680 - t668;
t683 = t351 * t505 + t352 * t581 - t365 * t655 + t366 * t656;
t682 = pkin(3) + pkin(8);
t410 = t442 * t680 - t558 * t443;
t643 = qJD(6) * t564;
t634 = t564 * t410 + t524 * t643 + t561 * t526;
t644 = qJD(6) * t561;
t382 = -t582 * t644 + t634;
t679 = t382 * t561;
t678 = t409 * t581;
t677 = t414 * t582;
t676 = t416 * t582;
t674 = t486 * t524;
t673 = t488 * t524;
t672 = t505 * t564;
t671 = t524 * t567;
t555 = t559 ^ 2;
t670 = t555 * qJD(1) ^ 2;
t406 = t564 * t409;
t663 = t562 * pkin(4) + qJ(3);
t501 = t560 * t562 + t565 * t666;
t629 = t559 * t649;
t468 = -qJD(4) * t501 + t562 * t629;
t502 = t560 * t565 - t562 * t666;
t531 = pkin(2) * t629;
t449 = t531 + t569;
t479 = (t666 * t682 + t545) * qJD(2);
t571 = -qJD(4) * t658 - t449 * t562 + t565 * t479;
t648 = qJD(2) * t566;
t628 = t559 * t648;
t373 = pkin(4) * t628 - qJ(5) * t468 - qJD(5) * t502 + t571;
t469 = -t559 * t576 + t560 * t645;
t579 = t565 * t449 + t456 * t645 - t474 * t646 + t562 * t479;
t377 = -qJ(5) * t469 - qJD(5) * t501 + t579;
t358 = t558 * t373 + t680 * t377;
t613 = t565 * t456 - t474 * t562;
t400 = pkin(4) * t667 - qJ(5) * t502 + t613;
t405 = -qJ(5) * t501 + t658;
t379 = t558 * t400 + t680 * t405;
t660 = pkin(5) * t630 - t661;
t654 = pkin(8) * t602 - t566 * t605;
t636 = t524 * t563 * t565;
t635 = t566 * t670;
t463 = -t540 * qJD(3) + t654;
t490 = -t560 * qJ(3) - t688;
t624 = t524 * t645;
t622 = t555 * t639;
t350 = pkin(10) * t526 + t352;
t437 = -pkin(3) * t602 - t463;
t411 = pkin(4) * t443 + t437;
t362 = pkin(5) * t409 - pkin(10) * t410 + t411;
t617 = -t350 * t561 + t564 * t362;
t616 = t410 * t561 - t564 * t526;
t615 = t561 * t655 - t564 * t630;
t614 = t561 * t630 + t564 * t655;
t607 = -qJD(6) * t581 - t540;
t604 = t563 * t635;
t527 = t540 * qJ(3);
t446 = t527 + t478;
t473 = pkin(3) * t666 - t490;
t601 = -0.2e1 * pkin(1) * t622;
t600 = t497 * t540 - t489;
t445 = pkin(5) * t581 - pkin(10) * t505 + t663;
t599 = pkin(10) * t630 - qJD(6) * t445 - t659;
t513 = t662 * t562;
t461 = -t513 * t680 - t558 * t618;
t598 = -pkin(5) * t656 - pkin(10) * t655 + qJD(6) * t461 - t691;
t596 = t350 * t564 + t362 * t561;
t360 = pkin(10) * t524 + t366;
t419 = pkin(4) * t486 + qJD(5) + t446;
t384 = -pkin(5) * t611 - pkin(10) * t582 + t419;
t354 = t360 * t564 + t384 * t561;
t594 = t360 * t561 - t384 * t564;
t375 = pkin(10) * t667 + t379;
t438 = t501 * t680 + t502 * t558;
t439 = -t558 * t501 + t502 * t680;
t585 = pkin(4) * t501 + t473;
t390 = pkin(5) * t438 - pkin(10) * t439 + t585;
t593 = t375 * t564 + t390 * t561;
t592 = -t375 * t561 + t390 * t564;
t498 = t688 * qJD(2);
t591 = t489 * t560 + t498 * t540;
t535 = t566 * t637;
t590 = -pkin(8) * t629 + t535;
t589 = -t540 * t630 + t526;
t588 = t406 + (t561 * t611 - t644) * t686;
t587 = -t439 * t561 + t564 * t667;
t421 = t439 * t564 + t561 * t667;
t586 = t446 * t563 + t567 * t648;
t584 = t524 * t562;
t357 = t373 * t680 - t558 * t377;
t378 = t400 * t680 - t558 * t405;
t491 = (-pkin(2) * t566 + t620) * t559;
t578 = t505 * t643 - t615;
t577 = -t505 * t644 - t614;
t552 = t560 * qJD(3);
t455 = -t629 * t682 + t535 + t552;
t575 = (-qJ(3) * t648 - t647) * t559;
t359 = -t524 * pkin(5) - t365;
t370 = t395 * t680 - t669;
t573 = -t549 * t409 + (t359 + t370) * t686;
t572 = pkin(4) * t469 + t455;
t550 = -pkin(4) * t680 - pkin(5);
t512 = t565 * t526;
t495 = -qJ(3) * t630 + t536;
t492 = t560 * t633 + t542;
t484 = -t552 - t590;
t483 = qJD(1) * t491;
t480 = t531 + t575;
t475 = -t527 - t497;
t470 = -pkin(2) * t540 + t687;
t460 = -t513 * t558 + t618 * t680;
t459 = qJD(1) * t575 + t519;
t454 = t483 * t631;
t418 = t468 * t680 - t558 * t469;
t417 = t468 * t558 + t469 * t680;
t388 = qJD(6) * t421 + t418 * t561 - t564 * t628;
t387 = qJD(6) * t587 + t418 * t564 + t561 * t628;
t383 = qJD(6) * t416 + t616;
t374 = -pkin(5) * t667 - t378;
t371 = pkin(5) * t417 - pkin(10) * t418 + t572;
t369 = t395 * t558 + t393;
t356 = pkin(10) * t628 + t358;
t355 = -pkin(5) * t628 - t357;
t348 = -qJD(6) * t354 + t617;
t347 = -qJD(6) * t594 + t596;
t1 = [(t468 * t524 + (t442 * t563 + (qJD(1) * t502 + t488) * t648) * t559) * MDP(17) + (-t469 * t524 + (-t443 * t563 + (-qJD(1) * t501 - t486) * t648) * t559) * MDP(18) + (t571 * t524 + t455 * t486 + t473 * t443 + t437 * t501 + t446 * t469 + (t570 * t563 + (qJD(1) * t613 + t402) * t648) * t559) * MDP(20) + (-t463 * t560 - t484 * t540 + (-t483 * t648 - t459 * t563 + (-t480 * t563 - t491 * t648) * qJD(1)) * t559) * MDP(13) + (-t579 * t524 + t455 * t488 + t473 * t442 + t437 * t502 + t446 * t468 + (-t580 * t563 + (-qJD(1) * t658 - t403) * t648) * t559) * MDP(21) + (t524 * t559 + t555 * t650) * MDP(19) * t648 + (t563 * t601 - t591) * MDP(9) + (-(qJD(6) * t592 + t356 * t564 + t371 * t561) * t686 - t593 * t409 - t347 * t438 - t354 * t417 + t355 * t416 + t374 * t382 + t349 * t421 + t359 * t387) * MDP(30) + ((-qJD(6) * t593 - t356 * t561 + t371 * t564) * t686 + t592 * t409 + t348 * t438 - t594 * t417 + t355 * t414 + t374 * t383 - t349 * t587 + t359 * t388) * MDP(29) + (t351 * t378 + t352 * t379 + t365 * t357 + t366 * t358 + t411 * t585 + t419 * t572) * MDP(23) + (-t540 * t590 + t560 * t654 + t566 * t601) * MDP(10) + (-t442 * t501 - t443 * t502 - t468 * t486 - t469 * t488) * MDP(16) + (t442 * t502 + t468 * t488) * MDP(15) + (t459 * t491 + t463 * t490 + t470 * t498 + t475 * t484 + t480 * t483 + t489 * t492) * MDP(14) + (-t463 * t566 + t489 * t563 + (t470 * t566 + t475 * t563) * qJD(2) + (-t484 * t566 + t498 * t563 + (t490 * t563 + t492 * t566) * qJD(2)) * qJD(1)) * t559 * MDP(11) + (t382 * t421 + t387 * t416) * MDP(24) + (t382 * t587 - t383 * t421 - t387 * t414 - t388 * t416) * MDP(25) + ((-t483 * t649 + t459 * t566 + (t480 * t566 - t491 * t649) * qJD(1)) * t559 + t591) * MDP(12) + (t409 * t438 + t417 * t686) * MDP(28) + (-t383 * t438 - t388 * t686 + t409 * t587 - t414 * t417) * MDP(27) + (t382 * t438 + t387 * t686 + t409 * t421 + t416 * t417) * MDP(26) + (-t351 * t439 - t352 * t438 - t357 * t582 + t358 * t611 - t365 * t418 - t366 * t417 - t378 * t410 - t379 * t409) * MDP(22) + (MDP(6) * t628 - MDP(7) * t629) * (t540 + t651) + 0.2e1 * (t563 * t566 * MDP(4) - t690) * t622; (-t351 * t460 + t352 * t461 + t661 * t365 + t659 * t366 + t411 * t663 + t419 * t691) * MDP(23) + (t656 * t686 + t678) * MDP(28) + ((t445 * t564 - t461 * t561) * t409 + t348 * t581 + t460 * t383 + t349 * t561 * t505 + (t561 * t599 - t564 * t598) * t686 + t660 * t414 - t656 * t594 + t578 * t359) * MDP(29) + (-(t445 * t561 + t461 * t564) * t409 - t347 * t581 + t460 * t382 + t349 * t672 + (t561 * t598 + t564 * t599) * t686 + t660 * t416 - t656 * t354 + t577 * t359) * MDP(30) + (-t383 * t581 - t414 * t656 - t505 * t665 - t578 * t686) * MDP(27) + (t382 * t581 + t406 * t505 + t416 * t656 + t577 * t686) * MDP(26) + (-pkin(2) * t489 - qJ(3) * t463 - t470 * t497 - t475 * t687 - t483 * t495) * MDP(14) + (t687 * t540 + (t483 * t566 + t495 * t563) * t652 - t463) * MDP(13) + (-t461 * t409 + t410 * t460 - t582 * t661 + t611 * t659 - t683) * MDP(22) + (-qJD(2) + t540) * MDP(7) * t631 + (t442 * t565 - t488 * t584) * MDP(15) + ((-qJ(3) * qJD(2) - t475 - t497) * t563 + (-pkin(2) * qJD(2) - t470 + t687) * t566) * MDP(11) * t652 + (t670 * t681 + t600) * MDP(9) + (t615 * t416 + t614 * t414 + (-t679 - t383 * t564 + (t414 * t561 - t416 * t564) * qJD(6)) * t505) * MDP(25) + ((-t443 - t673) * t565 + (-t442 + t674) * t562) * MDP(16) + (qJ(3) * t442 + t437 * t565 + t657 * t524 + t642 * t488 + (-t446 * t562 - t565 * t671) * qJD(4) + (t403 * t566 - t562 * t586) * t652) * MDP(21) + (qJ(3) * t443 + t437 * t562 - t612 * t524 + t642 * t486 + (t446 * t565 - t562 * t671) * qJD(4) + (-t402 * t566 + t565 * t586) * t652) * MDP(20) + (t382 * t672 + t416 * t577) * MDP(24) + (-t524 * t646 + t512 + (-t488 * t566 - t524 * t664) * t652) * MDP(17) + (-t495 * t630 + t454 - t600) * MDP(12) - MDP(4) * t604 - t524 * MDP(19) * t630 + t589 * MDP(6) + (pkin(1) * t635 - t496 * t540 + t654) * MDP(10) + t670 * t690 + (-t624 + (-t636 + (-qJD(2) * t562 + t486) * t566) * t652) * MDP(18); t589 * MDP(11) + MDP(12) * t604 + (-t540 ^ 2 - t556 * t670) * MDP(13) + (t475 * t540 + t454 + t489) * MDP(14) + (-t486 * t540 - t524 * t584 + t512) * MDP(20) + (-t624 - t488 * t540 + (-t562 * t648 - t636) * t652) * MDP(21) + (-t410 * t505 + t582 * t655 + t611 * t656 - t678) * MDP(22) + (-t419 * t540 + t683) * MDP(23) + (-t581 * t665 - t505 * t383 + t655 * t414 + (-t561 * t656 + t564 * t607) * t686) * MDP(29) + (-t581 * t406 - t505 * t382 + t655 * t416 + (-t561 * t607 - t564 * t656) * t686) * MDP(30); t488 * t486 * MDP(15) + (-t486 ^ 2 + t488 ^ 2) * MDP(16) + (t442 + t674) * MDP(17) + (-t443 + t673) * MDP(18) + MDP(19) * t526 + (t403 * t524 - t446 * t488 + t570) * MDP(20) + (t402 * t524 + t446 * t486 - t580) * MDP(21) + ((-t409 * t558 - t410 * t680) * pkin(4) + (t365 - t370) * t611 + (t366 - t369) * t582) * MDP(22) + (t365 * t369 - t366 * t370 + (t351 * t680 + t352 * t558 - t419 * t488) * pkin(4)) * MDP(23) + (t416 * t610 + t679) * MDP(24) + ((t382 - t694) * t564 + (-t383 - t693) * t561) * MDP(25) + (-t676 - t692) * MDP(26) + (t588 + t677) * MDP(27) - t686 * t582 * MDP(28) + (-t369 * t414 + t550 * t383 + t573 * t561 - t564 * t685 + t582 * t594) * MDP(29) + (t354 * t582 - t369 * t416 + t550 * t382 + t561 * t685 + t573 * t564) * MDP(30); (-t582 ^ 2 - t611 ^ 2) * MDP(22) + (t365 * t582 - t366 * t611 + t411) * MDP(23) + (t588 - t677) * MDP(29) + (-t676 + t692) * MDP(30); t416 * t414 * MDP(24) + (-t414 ^ 2 + t416 ^ 2) * MDP(25) + (t634 + t694) * MDP(26) + (-t616 + t693) * MDP(27) + t409 * MDP(28) + (t354 * t686 - t359 * t416 + t617) * MDP(29) + (t359 * t414 - t594 * t686 - t596) * MDP(30) + (-MDP(26) * t675 - MDP(27) * t416 - MDP(29) * t354 + MDP(30) * t594) * qJD(6);];
tauc  = t1;
