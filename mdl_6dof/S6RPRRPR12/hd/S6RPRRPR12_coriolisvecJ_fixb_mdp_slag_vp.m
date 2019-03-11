% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR12_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR12_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR12_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:53:13
% EndTime: 2019-03-09 05:53:29
% DurationCPUTime: 9.31s
% Computational Cost: add. (11346->564), mult. (37575->752), div. (0->0), fcn. (32135->12), ass. (0->228)
t537 = sin(pkin(12));
t539 = cos(pkin(12));
t542 = sin(qJ(3));
t670 = cos(pkin(7));
t675 = cos(qJ(3));
t605 = t670 * t675;
t695 = t537 * t605 + t539 * t542;
t538 = sin(pkin(6));
t622 = t542 * t670;
t562 = t537 * t622 - t539 * t675;
t559 = t538 * t562;
t669 = sin(pkin(7));
t604 = t669 * t675;
t694 = qJD(1) * t559 + qJD(3) * t604;
t541 = sin(qJ(4));
t671 = cos(pkin(6));
t595 = t671 * t669;
t489 = t538 * (t537 * t675 + t539 * t622) + t542 * t595;
t484 = qJD(1) * t489;
t544 = cos(qJ(4));
t647 = t544 * t484;
t596 = t671 * t670;
t618 = qJD(1) * t669;
t603 = t538 * t618;
t517 = t539 * t603;
t632 = qJD(3) - t517;
t682 = qJD(1) * t596 + t632;
t451 = t541 * t682 + t647;
t448 = qJD(6) + t451;
t637 = qJD(2) * t538;
t693 = t695 * t637;
t569 = t675 * t595;
t586 = t539 * t605;
t651 = t537 * t542;
t482 = (t569 + (t586 - t651) * t538) * qJD(3);
t548 = qJD(1) * t482;
t692 = -qJD(4) * t682 - t548;
t691 = -t538 * t586 - t569;
t690 = MDP(4) * t537 + MDP(5) * t539;
t638 = qJD(1) * t538;
t626 = t537 * t638;
t480 = qJD(1) * t691 + t542 * t626;
t475 = qJD(4) + t480;
t619 = qJD(1) * t671;
t608 = pkin(1) * t619;
t527 = t539 * t608;
t652 = t537 * t538;
t554 = t671 * pkin(2) + (-pkin(9) * t670 - qJ(2)) * t652;
t477 = qJD(1) * t554 + t527;
t500 = (-pkin(9) * t537 * t669 - pkin(2) * t539 - pkin(1)) * t538;
t494 = qJD(1) * t500 + qJD(2);
t445 = -t477 * t669 + t670 * t494;
t400 = t480 * pkin(3) - t484 * pkin(10) + t445;
t650 = t538 * t539;
t528 = qJ(2) * t650;
t504 = qJD(1) * t528 + t537 * t608;
t553 = (t650 * t670 + t595) * pkin(9);
t470 = qJD(1) * t553 + t504;
t621 = t542 * t669;
t421 = t675 * t470 + t477 * t622 + t494 * t621;
t403 = pkin(10) * t682 + t421;
t372 = -t544 * t400 + t403 * t541;
t630 = qJD(5) + t372;
t687 = (t537 ^ 2 + t539 ^ 2) * MDP(6) * t538 ^ 2;
t686 = -t542 * t470 + t477 * t605 + t494 * t604;
t508 = t541 * t670 + t544 * t621;
t585 = t537 * t603;
t643 = -qJD(4) * t508 - t541 * t694 - t544 * t585;
t507 = t541 * t621 - t544 * t670;
t642 = -qJD(4) * t507 - t541 * t585 + t544 * t694;
t627 = pkin(1) * t671;
t640 = t537 * t627 + t528;
t486 = t553 + t640;
t531 = t539 * t627;
t490 = t531 + t554;
t685 = -t542 * t486 + t490 * t605 + t500 * t604;
t601 = qJD(3) * t621;
t579 = -t638 * t695 + t601;
t636 = qJD(4) * t541;
t654 = t480 * t541;
t684 = qJD(5) * t541 + t421 + (-t636 - t654) * pkin(4);
t483 = t489 * qJD(3);
t468 = qJD(1) * t483;
t631 = pkin(5) * t451 + t630;
t373 = t541 * t400 + t544 * t403;
t370 = -qJ(5) * t475 - t373;
t449 = t484 * t541 - t544 * t682;
t673 = pkin(5) * t449;
t362 = -t370 - t673;
t426 = t484 * t636 + t544 * t692;
t677 = pkin(4) + pkin(11);
t683 = t677 * t426 + (t362 - t373 + t673) * t448;
t679 = t451 ^ 2;
t678 = t475 ^ 2;
t676 = pkin(5) + pkin(10);
t674 = pkin(4) * t468;
t672 = pkin(10) * t468;
t667 = qJ(5) * t449;
t666 = qJ(5) * t544;
t665 = t370 * t475;
t664 = t373 * t475;
t653 = t482 * t541;
t427 = (t541 * t632 + t647) * qJD(4) + (t596 * t636 + t653) * qJD(1);
t540 = sin(qJ(6));
t543 = cos(qJ(6));
t633 = qJD(6) * t543;
t628 = t540 * t427 + t449 * t633 + t543 * t468;
t634 = qJD(6) * t540;
t379 = -t475 * t634 + t628;
t663 = t379 * t543;
t655 = t475 * t540;
t428 = -t543 * t449 + t655;
t662 = t428 * t448;
t661 = t428 * t475;
t430 = t449 * t540 + t475 * t543;
t660 = t430 * t448;
t659 = t430 * t475;
t658 = t449 * t475;
t657 = t451 * t449;
t656 = t451 * t475;
t466 = t468 * qJ(5);
t648 = t540 * t426;
t423 = t543 * t426;
t635 = qJD(4) * t544;
t646 = qJ(5) * t635 + t480 * t666 + t684;
t453 = -t490 * t669 + t670 * t500;
t488 = t538 * t651 + t691;
t412 = t488 * pkin(3) - t489 * pkin(10) + t453;
t623 = t538 * t669;
t506 = t539 * t623 - t596;
t419 = t675 * t486 - t506 * pkin(10) + (t490 * t670 + t500 * t669) * t542;
t645 = t541 * t412 + t544 * t419;
t444 = pkin(3) * t484 + pkin(10) * t480;
t644 = t541 * t444 + t544 * t686;
t381 = -qJ(5) * t484 - t644;
t641 = -pkin(5) * t654 - t676 * t636 + t381;
t629 = pkin(10) * t636;
t525 = t676 * t544;
t625 = qJD(3) * t675;
t624 = -qJ(5) * t541 - pkin(3);
t555 = qJD(2) * t559;
t392 = -qJD(1) * t555 + qJD(3) * t686;
t584 = t537 * qJD(2) * t623;
t572 = qJD(1) * t584;
t436 = t468 * pkin(3) - pkin(10) * t548 + t572;
t610 = t541 * t392 + t400 * t636 + t403 * t635 - t544 * t436;
t350 = -pkin(5) * t426 - t468 * t677 + t610;
t602 = qJD(3) * t622;
t393 = qJD(1) * t693 + t470 * t625 + t477 * t602 + t494 * t601;
t558 = qJ(5) * t426 - qJD(5) * t451 + t393;
t357 = t427 * t677 + t558;
t617 = t543 * t350 - t357 * t540;
t616 = t412 * t544 - t541 * t419;
t615 = -t543 * t427 + t468 * t540;
t614 = t475 * t544;
t613 = t448 * t540;
t611 = t544 * t392 + t400 * t635 - t403 * t636 + t541 * t436;
t375 = -qJ(5) * t488 - t645;
t441 = t484 * t540 + t543 * t654;
t600 = t543 * t636 + t441;
t442 = t484 * t543 - t540 * t654;
t599 = t540 * t636 - t442;
t416 = t541 * t686;
t510 = -t544 * t677 + t624;
t598 = qJD(6) * t510 + t416 + (-pkin(5) * t480 - t444) * t544 - t677 * t484 - qJD(4) * t525;
t524 = t676 * t541;
t597 = -qJD(6) * t524 + t684 - t475 * (pkin(11) * t541 - t666);
t407 = t486 * t625 + t490 * t602 + t500 * t601 + t693;
t593 = t350 * t540 + t357 * t543;
t360 = -t475 * t677 + t631;
t402 = -pkin(3) * t682 - t686;
t547 = -t451 * qJ(5) + t402;
t368 = t449 * t677 + t547;
t347 = t360 * t543 - t368 * t540;
t348 = t360 * t540 + t368 * t543;
t456 = t489 * t544 - t506 * t541;
t364 = pkin(5) * t456 - t488 * t677 - t616;
t455 = t489 * t541 + t506 * t544;
t418 = t506 * pkin(3) - t685;
t549 = -t456 * qJ(5) + t418;
t374 = t455 * t677 + t549;
t592 = t364 * t543 - t374 * t540;
t591 = t364 * t540 + t374 * t543;
t589 = t455 * t543 - t488 * t540;
t438 = t455 * t540 + t488 * t543;
t588 = (-qJ(2) * t626 + t527) * t537 - t504 * t539;
t406 = qJD(3) * t685 - t555;
t440 = t483 * pkin(3) - t482 * pkin(10) + t584;
t583 = -t541 * t406 - t412 * t636 - t419 * t635 + t440 * t544;
t354 = -t475 * qJD(5) - t466 - t611;
t577 = -t448 * t613 - t423;
t575 = t402 * t475 - t672;
t378 = t449 * pkin(4) + t547;
t574 = -t378 * t475 + t672;
t573 = t378 * t451 + t610;
t571 = t544 * t406 + t412 * t635 - t419 * t636 + t541 * t440;
t351 = -pkin(5) * t427 - t354;
t570 = t351 + (t448 * t677 + t667) * t448;
t568 = -t448 ^ 2 * t543 + t648;
t567 = -t426 + t658;
t564 = t543 * t507 + t540 * t604;
t563 = -t540 * t507 + t543 * t604;
t435 = -t489 * t636 + (-qJD(4) * t506 + t482) * t544;
t557 = -qJ(5) * t435 - qJD(5) * t456 + t407;
t356 = -qJ(5) * t483 - qJD(5) * t488 - t571;
t520 = -pkin(4) * t544 + t624;
t434 = qJD(4) * t456 + t653;
t425 = t426 * t541;
t410 = pkin(4) * t451 + t667;
t389 = t426 * t456;
t386 = qJD(6) * t589 + t434 * t540 + t483 * t543;
t385 = qJD(6) * t438 - t434 * t543 + t483 * t540;
t384 = t455 * pkin(4) + t549;
t383 = -pkin(4) * t484 - t444 * t544 + t416;
t380 = qJD(6) * t430 + t615;
t376 = -pkin(4) * t488 - t616;
t369 = -pkin(4) * t475 + t630;
t367 = -pkin(5) * t455 - t375;
t363 = pkin(4) * t434 + t557;
t361 = pkin(4) * t427 + t558;
t359 = t434 * t677 + t557;
t358 = -pkin(4) * t483 - t583;
t355 = t610 - t674;
t353 = -pkin(5) * t434 - t356;
t352 = pkin(5) * t435 - t483 * t677 - t583;
t346 = -qJD(6) * t348 + t617;
t345 = qJD(6) * t347 + t593;
t1 = [(-(qJD(6) * t592 + t352 * t540 + t359 * t543) * t448 + t591 * t426 - t345 * t456 - t348 * t435 + t353 * t430 + t367 * t379 + t351 * t438 + t362 * t386) * MDP(32) + (-t489 * t468 - t482 * t480 - t484 * t483 - t488 * t548) * MDP(9) + (t484 * t482 + t489 * t548) * MDP(8) + 0.2e1 * qJD(2) * qJD(1) * t687 + (-t426 * t488 + t435 * t475 + t451 * t483 + t456 * t468) * MDP(17) + (-t354 * t488 - t356 * t475 - t361 * t456 - t363 * t451 - t370 * t483 - t375 * t468 - t378 * t435 + t384 * t426) * MDP(24) + (-0.2e1 * t690 * t619 + (t480 * t669 + t488 * t618) * t537 * MDP(13) + ((t539 * t640 + (qJ(2) * t652 - t531) * t537) * qJD(1) - t588) * MDP(7)) * t637 + (t354 * t375 + t355 * t376 + t356 * t370 + t358 * t369 + t361 * t384 + t363 * t378) * MDP(25) + (t435 * t451 - t389) * MDP(15) + (t435 * t448 - t389) * MDP(30) + (t468 * t488 + t475 * t483) * MDP(19) + (t355 * t488 + t358 * t475 - t361 * t455 - t363 * t449 + t369 * t483 + t376 * t468 - t378 * t434 - t384 * t427) * MDP(23) + (-t427 * t488 - t434 * t475 - t449 * t483 - t455 * t468) * MDP(18) + ((-qJD(6) * t591 + t352 * t543 - t359 * t540) * t448 - t592 * t426 + t346 * t456 + t347 * t435 + t353 * t428 + t367 * t380 - t351 * t589 + t362 * t385) * MDP(31) + (t379 * t589 - t380 * t438 - t385 * t430 - t386 * t428) * MDP(27) + (-t380 * t456 - t385 * t448 - t426 * t589 - t428 * t435) * MDP(29) + (t379 * t438 + t386 * t430) * MDP(26) + (t379 * t456 + t386 * t448 - t426 * t438 + t430 * t435) * MDP(28) + (t354 * t455 + t355 * t456 + t356 * t449 + t358 * t451 + t369 * t435 + t370 * t434 + t375 * t427 - t376 * t426) * MDP(22) + (t426 * t455 - t427 * t456 - t434 * t451 - t435 * t449) * MDP(16) + (t632 + (-t506 + t596) * qJD(1)) * t482 * MDP(10) + (-t372 * t483 + t393 * t455 + t402 * t434 + t407 * t449 + t418 * t427 + t468 * t616 + t475 * t583 - t488 * t610) * MDP(20) + (t393 * t506 - t407 * t682 + t445 * t483 + t453 * t468) * MDP(13) + (t392 * t506 - t406 * t682 + t445 * t482 + t453 * t548 + t484 * t584 + t489 * t572) * MDP(14) + (t468 * t506 - t483 * t682) * MDP(11) + (-t373 * t483 + t393 * t456 + t402 * t435 + t407 * t451 - t418 * t426 - t468 * t645 - t475 * t571 - t488 * t611) * MDP(21); t588 * MDP(7) * t638 + (t468 * t670 - t480 * t585 - t579 * t682) * MDP(13) + (-t484 * t585 + t548 * t670 - t682 * t694) * MDP(14) + (-t507 * t426 - t427 * t508 - t449 * t642 - t451 * t643) * MDP(22) + (-t354 * t508 + t355 * t507 - t361 * t604 - t369 * t643 - t370 * t642 + t378 * t579) * MDP(25) + (t508 * t380 - t564 * t426 + (t563 * qJD(6) - t540 * t579 - t543 * t643) * t448 + t642 * t428) * MDP(31) + (t508 * t379 - t563 * t426 + (-t564 * qJD(6) + t540 * t643 - t543 * t579) * t448 + t642 * t430) * MDP(32) + (MDP(20) - MDP(23)) * (-t427 * t604 + t449 * t579 - t507 * t468 + t475 * t643) + (MDP(21) - MDP(24)) * (t426 * t604 + t451 * t579 - t508 * t468 - t475 * t642) + (t538 * t671 * t690 - t687) * qJD(1) ^ 2; -t480 ^ 2 * MDP(9) + (t480 * t632 + (t480 * t596 + t482) * qJD(1)) * MDP(10) - MDP(11) * t468 + (t421 * t682 - t393) * MDP(13) + (-t686 * t517 + t445 * t480 + (t562 * t637 + t596 * t686) * qJD(1)) * MDP(14) + (t451 * t614 - t425) * MDP(15) + ((-t426 - t658) * t544 + (-t427 - t656) * t541) * MDP(16) + (t468 * t541 + t475 * t614) * MDP(17) + (t468 * t544 - t541 * t678) * MDP(18) + (-pkin(3) * t427 - t393 * t544 - t421 * t449 + (t416 + (-pkin(10) * qJD(4) - t444) * t544) * t475 + t575 * t541) * MDP(20) + (pkin(3) * t426 + t393 * t541 - t421 * t451 + (t629 + t644) * t475 + t575 * t544) * MDP(21) + (-t381 * t449 - t383 * t451 + (-t354 + t475 * t369 + (qJD(4) * t451 - t427) * pkin(10)) * t544 + (t355 + t665 + (qJD(4) * t449 - t426) * pkin(10)) * t541) * MDP(22) + (t361 * t544 - t427 * t520 + (pkin(10) * t635 - t383) * t475 + t646 * t449 + t574 * t541) * MDP(23) + (-t361 * t541 + t426 * t520 + (t381 - t629) * t475 + t646 * t451 + t574 * t544) * MDP(24) + (t361 * t520 - t369 * t383 - t370 * t381 - t646 * t378 + (-t354 * t544 + t355 * t541 + (t369 * t544 + t370 * t541) * qJD(4)) * pkin(10)) * MDP(25) + (-t379 * t540 * t544 + (-t544 * t633 + t599) * t430) * MDP(26) + (t428 * t442 + t430 * t441 + (-t428 * t540 + t430 * t543) * t636 + (-t663 + t380 * t540 + (t428 * t543 + t430 * t540) * qJD(6)) * t544) * MDP(27) + (t379 * t541 + t599 * t448 + (-t448 * t633 + t648 + t659) * t544) * MDP(28) + (-t380 * t541 + t600 * t448 + (t448 * t634 + t423 - t661) * t544) * MDP(29) + (t448 * t614 - t425) * MDP(30) + (-(-t510 * t540 + t524 * t543) * t426 + t346 * t541 + t525 * t380 + (t540 * t597 - t543 * t598) * t448 + t641 * t428 - t600 * t362 + (t347 * t475 + t351 * t543 - t362 * t634) * t544) * MDP(31) + ((t510 * t543 + t524 * t540) * t426 - t345 * t541 + t525 * t379 + (t540 * t598 + t543 * t597) * t448 + t641 * t430 + t599 * t362 + (-t348 * t475 - t351 * t540 - t362 * t633) * t544) * MDP(32) + (MDP(11) * t682 - t445 * MDP(13) - t451 * MDP(17) + t449 * MDP(18) - t475 * MDP(19) + t372 * MDP(20) + t373 * MDP(21) - t369 * MDP(23) + t370 * MDP(24) + t480 * MDP(8) + MDP(9) * t484) * t484; MDP(15) * t657 + (-t449 ^ 2 + t679) * MDP(16) + t567 * MDP(17) + (-t484 * t635 + t541 * t692 + t656) * MDP(18) + t468 * MDP(19) + (-t402 * t451 - t610 + t664) * MDP(20) + (-t372 * t475 + t402 * t449 - t611) * MDP(21) + (pkin(4) * t426 - qJ(5) * t427 + (-t370 - t373) * t451 + (t369 - t630) * t449) * MDP(22) + (t410 * t449 + t573 - t664 - 0.2e1 * t674) * MDP(23) + (-t378 * t449 + t410 * t451 + t475 * t630 - t354 + t466) * MDP(24) + (-pkin(4) * t355 - qJ(5) * t354 - t369 * t373 - t370 * t630 - t378 * t410) * MDP(25) + (-t430 * t613 + t663) * MDP(26) + ((-t380 - t660) * t543 + (-t379 + t662) * t540) * MDP(27) + (t430 * t449 + t577) * MDP(28) + (-t428 * t449 + t568) * MDP(29) + t448 * t449 * MDP(30) + (qJ(5) * t380 + t347 * t449 + t631 * t428 + t570 * t540 + t543 * t683) * MDP(31) + (qJ(5) * t379 - t348 * t449 + t631 * t430 - t540 * t683 + t570 * t543) * MDP(32); t567 * MDP(22) + (t468 - t657) * MDP(23) + (-t678 - t679) * MDP(24) + (t573 + t665 - t674) * MDP(25) + (t577 - t661) * MDP(31) + (t568 - t659) * MDP(32); t430 * t428 * MDP(26) + (-t428 ^ 2 + t430 ^ 2) * MDP(27) + (t628 + t662) * MDP(28) + (-t615 + t660) * MDP(29) - t426 * MDP(30) + (t348 * t448 - t362 * t430 + t617) * MDP(31) + (t347 * t448 + t362 * t428 - t593) * MDP(32) + (-MDP(28) * t655 - MDP(29) * t430 - MDP(31) * t348 - MDP(32) * t347) * qJD(6);];
tauc  = t1;
