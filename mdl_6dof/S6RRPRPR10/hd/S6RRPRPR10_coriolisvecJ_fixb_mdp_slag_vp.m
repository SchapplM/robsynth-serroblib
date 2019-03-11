% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:47
% EndTime: 2019-03-09 11:10:01
% DurationCPUTime: 9.25s
% Computational Cost: add. (8295->573), mult. (22274->758), div. (0->0), fcn. (17912->10), ass. (0->242)
t575 = cos(pkin(6));
t673 = qJD(1) * t575;
t564 = qJD(2) + t673;
t572 = sin(pkin(11));
t574 = cos(pkin(11));
t578 = sin(qJ(2));
t573 = sin(pkin(6));
t674 = qJD(1) * t573;
t652 = t578 * t674;
t518 = t564 * t574 - t572 * t652;
t519 = t564 * t572 + t574 * t652;
t577 = sin(qJ(4));
t707 = cos(qJ(4));
t469 = -t707 * t518 + t519 * t577;
t579 = cos(qJ(6));
t580 = cos(qJ(2));
t672 = qJD(1) * t580;
t651 = t573 * t672;
t555 = -qJD(4) + t651;
t576 = sin(qJ(6));
t693 = t555 * t576;
t443 = -t579 * t469 - t693;
t606 = -t577 * t518 - t519 * t707;
t717 = qJD(6) - t606;
t724 = t443 * t717;
t445 = t469 * t576 - t555 * t579;
t723 = t445 * t717;
t653 = t707 * t574;
t689 = t573 * t580;
t622 = t653 * t689;
t607 = qJD(1) * t622;
t634 = t572 * t651;
t648 = qJD(4) * t707;
t668 = qJD(4) * t577;
t678 = t572 * t668 - t574 * t648 - t577 * t634 + t607;
t721 = t469 ^ 2;
t709 = t606 ^ 2;
t720 = t443 * t469;
t719 = t445 * t469;
t696 = t469 * t555;
t698 = t606 * t555;
t625 = pkin(2) * t578 - qJ(3) * t580;
t530 = t625 * t674;
t658 = pkin(1) * t673;
t531 = -pkin(8) * t652 + t580 * t658;
t479 = t574 * t530 - t531 * t572;
t687 = t574 * t580;
t598 = (pkin(3) * t578 - pkin(9) * t687) * t573;
t451 = qJD(1) * t598 + t479;
t480 = t572 * t530 + t574 * t531;
t461 = -pkin(9) * t634 + t480;
t704 = pkin(9) + qJ(3);
t550 = t704 * t572;
t551 = t704 * t574;
t718 = qJD(3) * t653 - t707 * t461 - t550 * t648 + (-qJD(3) * t572 - qJD(4) * t551 - t451) * t577;
t540 = t572 * t707 + t577 * t574;
t590 = t540 * t689;
t677 = -qJD(1) * t590 + t540 * qJD(4);
t532 = pkin(8) * t651 + t578 * t658;
t508 = qJ(3) * t564 + t532;
t527 = (-pkin(2) * t580 - qJ(3) * t578 - pkin(1)) * t573;
t513 = qJD(1) * t527;
t454 = -t508 * t572 + t574 * t513;
t413 = -pkin(3) * t651 - pkin(9) * t519 + t454;
t455 = t574 * t508 + t572 * t513;
t426 = pkin(9) * t518 + t455;
t387 = -t707 * t413 + t426 * t577;
t664 = -qJD(5) - t387;
t569 = t573 ^ 2;
t660 = qJD(1) * qJD(2);
t716 = -0.2e1 * t569 * t660;
t715 = MDP(5) * (t578 ^ 2 - t580 ^ 2);
t714 = t717 ^ 2;
t685 = qJ(5) * t652 - t718;
t499 = pkin(3) * t634 + t532;
t713 = -qJ(5) * t678 + qJD(5) * t540 + t499;
t505 = -t577 * t550 + t551 * t707;
t712 = -qJD(3) * t540 - qJD(4) * t505 - t451 * t707 + t577 * t461;
t711 = qJD(4) * t606;
t663 = -pkin(5) * t606 - t664;
t388 = t577 * t413 + t707 * t426;
t386 = qJ(5) * t555 - t388;
t705 = pkin(5) * t469;
t378 = -t386 - t705;
t646 = t573 * t660;
t631 = t580 * t646;
t615 = t572 * t631;
t424 = t577 * (qJD(4) * t519 + t615) - qJD(2) * t607 - t518 * t648;
t708 = pkin(4) + pkin(10);
t710 = t708 * t424 + (t378 - t388 + t705) * t717;
t511 = (qJD(2) * t625 - qJD(3) * t578) * t573;
t669 = qJD(2) * t578;
t650 = t573 * t669;
t657 = pkin(1) * qJD(2) * t575;
t612 = -pkin(8) * t650 + t580 * t657;
t517 = qJD(3) * t575 + t612;
t458 = t574 * t511 - t517 * t572;
t595 = qJD(2) * t598;
t431 = t458 + t595;
t526 = pkin(8) * t689 + (pkin(1) * t578 + qJ(3)) * t575;
t475 = -t526 * t572 + t574 * t527;
t690 = t573 * t578;
t534 = t572 * t575 + t574 * t690;
t433 = -pkin(3) * t689 - pkin(9) * t534 + t475;
t459 = t572 * t511 + t574 * t517;
t649 = qJD(2) * t689;
t633 = t572 * t649;
t446 = -pkin(9) * t633 + t459;
t476 = t574 * t526 + t572 * t527;
t656 = t572 * t690;
t688 = t574 * t575;
t608 = -t656 + t688;
t450 = pkin(9) * t608 + t476;
t601 = -t577 * t431 - t433 * t648 - t707 * t446 + t450 * t668;
t373 = -t573 * (qJ(5) * t669 - qJD(5) * t580) + t601;
t706 = pkin(1) * t580;
t703 = qJ(5) * t469;
t589 = qJD(2) * t590;
t425 = qJD(1) * t589 - t711;
t632 = t578 * t646;
t666 = qJD(6) * t579;
t654 = t576 * t425 + t469 * t666 + t579 * t632;
t667 = qJD(6) * t576;
t390 = t555 * t667 + t654;
t702 = t390 * t579;
t501 = -pkin(2) * t564 + qJD(3) - t531;
t466 = -pkin(3) * t518 + t501;
t586 = qJ(5) * t606 + t466;
t395 = pkin(4) * t469 + t586;
t701 = t395 * t606;
t700 = t454 * t578;
t699 = t455 * t578;
t697 = t469 * t606;
t638 = qJD(1) * t657;
t525 = pkin(8) * t631 + t578 * t638;
t695 = t525 * t572;
t539 = t572 * t577 - t653;
t694 = t539 * t576;
t582 = qJD(1) ^ 2;
t692 = t569 * t582;
t691 = t572 * t580;
t686 = t576 * t424;
t419 = t579 * t424;
t684 = -pkin(4) * t652 + t712;
t683 = -pkin(4) * t677 + t713;
t682 = t577 * t433 + t707 * t450;
t680 = -pkin(5) * t677 - t685;
t495 = qJD(1) * t511;
t604 = -pkin(8) * t632 + t580 * t638;
t496 = qJD(3) * t564 + t604;
t448 = t572 * t495 + t574 * t496;
t533 = pkin(8) * t649 + t578 * t657;
t504 = t550 * t707 + t577 * t551;
t671 = qJD(2) * t504;
t670 = qJD(2) * t505;
t665 = qJD(2) - t564;
t659 = pkin(1) * t692;
t655 = t579 * t689;
t491 = pkin(3) * t615 + t525;
t500 = pkin(3) * t633 + t533;
t568 = -pkin(3) * t574 - pkin(2);
t639 = t708 * t690;
t621 = qJD(2) * t639;
t447 = t574 * t495 - t496 * t572;
t414 = qJD(1) * t595 + t447;
t427 = -pkin(9) * t615 + t448;
t635 = t413 * t668 - t707 * t414 + t426 * t648 + t577 * t427;
t367 = -pkin(5) * t424 - qJD(1) * t621 + t635;
t594 = qJ(5) * t424 + qJD(5) * t606 + t491;
t372 = t425 * t708 + t594;
t645 = t579 * t367 - t372 * t576;
t644 = t576 * t652 + t579 * t677;
t643 = -t576 * t677 + t579 * t652;
t642 = t717 * t576;
t641 = t717 * t579;
t637 = t569 * t580 * t578 * MDP(4);
t636 = -t413 * t648 - t577 * t414 + t426 * t668 - t707 * t427;
t630 = pkin(1) * t716;
t629 = t433 * t707 - t577 * t450;
t611 = -qJ(5) * t540 + t568;
t474 = t539 * t708 + t611;
t627 = pkin(5) * t678 - qJD(1) * t639 + qJD(6) * t474 + t712;
t481 = t540 * pkin(5) + t504;
t626 = -qJD(6) * t481 - t677 * t708 + t713;
t624 = pkin(4) * t632;
t623 = -qJD(5) * t555 - t636;
t620 = t367 * t576 + t372 * t579;
t376 = t555 * t708 + t663;
t382 = t469 * t708 + t586;
t363 = t376 * t579 - t382 * t576;
t364 = t376 * t576 + t382 * t579;
t397 = pkin(4) * t689 - t629;
t486 = t534 * t707 + t577 * t608;
t383 = t486 * pkin(5) + pkin(10) * t689 + t397;
t591 = t707 * t608;
t485 = t534 * t577 - t591;
t565 = pkin(8) * t690;
t487 = pkin(3) * t656 + t565 + (t568 - t706) * t575;
t585 = -t486 * qJ(5) + t487;
t393 = t485 * t708 + t585;
t618 = t383 * t579 - t393 * t576;
t617 = t383 * t576 + t393 * t579;
t616 = qJ(5) * t632;
t396 = qJ(5) * t689 - t682;
t463 = t485 * t579 + t576 * t689;
t603 = -t388 * t555 - t635;
t602 = -t431 * t707 + t433 * t668 + t577 * t446 + t450 * t648;
t370 = -t616 - t623;
t366 = -pkin(5) * t425 - t370;
t600 = t366 + (t708 * t717 + t703) * t717;
t599 = -t579 * t425 + t576 * t632;
t597 = t539 * t667 - t644;
t596 = t539 * t666 - t643;
t438 = -qJD(4) * t591 - qJD(2) * t622 + (qJD(4) * t534 + t633) * t577;
t593 = qJ(5) * t438 - qJD(5) * t486 + t500;
t371 = -t624 + t635;
t588 = -qJ(3) * t669 + (-pkin(2) * qJD(2) + qJD(3) - t501) * t580;
t377 = pkin(4) * t425 + t594;
t584 = -t424 - t696;
t583 = -t540 * t631 + t711;
t529 = t565 + (-pkin(2) - t706) * t575;
t488 = pkin(4) * t539 + t611;
t482 = -t539 * pkin(5) + t505;
t464 = t485 * t576 - t655;
t439 = qJD(4) * t486 + t589;
t406 = t424 * t540;
t405 = -pkin(4) * t606 + t703;
t403 = t485 * pkin(4) + t585;
t402 = t424 * t486;
t399 = qJD(6) * t463 + t439 * t576 + t579 * t650;
t398 = -t579 * t439 - qJD(6) * t655 + (qJD(6) * t485 + t650) * t576;
t391 = qJD(6) * t445 + t599;
t385 = pkin(4) * t555 - t664;
t384 = -pkin(5) * t485 - t396;
t381 = pkin(4) * t439 + t593;
t375 = t439 * t708 + t593;
t374 = -pkin(4) * t650 + t602;
t369 = -pkin(5) * t439 - t373;
t368 = -t438 * pkin(5) + t602 - t621;
t362 = -qJD(6) * t364 + t645;
t361 = qJD(6) * t363 + t620;
t1 = [(t438 * t606 - t402) * MDP(15) + (t370 * t485 + t371 * t486 + t373 * t469 - t374 * t606 - t385 * t438 + t386 * t439 + t396 * t425 - t397 * t424) * MDP(22) + (t424 * t485 - t425 * t486 + t438 * t469 + t439 * t606) * MDP(16) + (-t601 * t555 - t500 * t606 - t487 * t424 + t491 * t486 - t466 * t438 + (-t636 * t580 + (-qJD(1) * t682 - t388) * t669) * t573) * MDP(21) + (t438 * t555 + (t424 * t580 + (qJD(1) * t486 - t606) * t669) * t573) * MDP(17) + (t373 * t555 - t377 * t486 + t381 * t606 + t395 * t438 + t403 * t424 + (t370 * t580 + (-qJD(1) * t396 - t386) * t669) * t573) * MDP(24) + (t390 * t463 - t391 * t464 - t398 * t445 - t399 * t443) * MDP(27) + (t390 * t464 + t399 * t445) * MDP(26) + (-t438 * t717 - t402) * MDP(30) + (-t391 * t486 - t398 * t717 - t424 * t463 + t438 * t443) * MDP(29) + (t390 * t486 + t399 * t717 - t424 * t464 - t438 * t445) * MDP(28) + (-(qJD(6) * t618 + t368 * t576 + t375 * t579) * t717 + t617 * t424 - t361 * t486 + t364 * t438 + t369 * t445 + t384 * t390 + t366 * t464 + t378 * t399) * MDP(32) + ((-qJD(6) * t617 + t368 * t579 - t375 * t576) * t717 - t618 * t424 + t362 * t486 - t363 * t438 + t369 * t443 + t384 * t391 - t366 * t463 + t378 * t398) * MDP(31) + (MDP(6) * t649 - MDP(7) * t650) * (t564 + t673) + (t370 * t396 + t371 * t397 + t373 * t386 + t374 * t385 + t377 * t403 + t381 * t395) * MDP(25) + (t447 * t475 + t448 * t476 + t454 * t458 + t455 * t459 + t501 * t533 + t525 * t529) * MDP(14) + (t519 * t533 + t525 * t534 + ((qJD(1) * t459 + t448) * t580 + (t501 * t687 - t699 + (-t476 * t578 + t529 * t687) * qJD(1)) * qJD(2)) * t573) * MDP(12) + (-t525 * t688 - t533 * t518 + (t578 * t695 + (-qJD(1) * t458 - t447) * t580 + (t501 * t691 + t700 + (t475 * t578 + t529 * t691) * qJD(1)) * qJD(2)) * t573) * MDP(11) + (-t374 * t555 - t377 * t485 - t381 * t469 - t395 * t439 - t403 * t425 + (-t371 * t580 + (qJD(1) * t397 + t385) * t669) * t573) * MDP(23) + (t439 * t555 + (t425 * t580 + (-qJD(1) * t485 - t469) * t669) * t573) * MDP(18) + (t602 * t555 + t500 * t469 + t487 * t425 + t491 * t485 + t466 * t439 + (t635 * t580 + (qJD(1) * t629 - t387) * t669) * t573) * MDP(20) + 0.2e1 * t637 * t660 + (t459 * t518 + t448 * t608 - t458 * t519 - t447 * t534 + (-t454 * t574 - t455 * t572 + (-t475 * t574 - t476 * t572) * qJD(1)) * t649) * MDP(13) + (-t555 * t573 - t569 * t672) * MDP(19) * t669 + (-t564 * t612 - t575 * t604 + t580 * t630) * MDP(10) + (-t525 * t575 - t533 * t564 + t578 * t630) * MDP(9) + t715 * t716; (-t377 * t539 - t677 * t395 - t425 * t488 + t683 * t469) * MDP(23) + (t568 * t425 + t677 * t466 - t499 * t469 + t491 * t539) * MDP(20) + (-pkin(2) * t525 - t454 * t479 - t455 * t480 - t501 * t532 + (-t454 * t572 + t455 * t574) * qJD(3) + (-t447 * t572 + t448 * t574) * qJ(3)) * MDP(14) + (-t377 * t540 + t678 * t395 + t424 * t488 - t606 * t683) * MDP(24) + (-t568 * t424 - t678 * t466 + t491 * t540 + t499 * t606) * MDP(21) + (t555 * MDP(19) + (t386 + t670) * MDP(24) + (t388 - t670) * MDP(21) + (-t385 + t671) * MDP(23) + (qJD(2) * t540 + t606) * MDP(17) + (t387 - t671) * MDP(20) + (-qJD(2) * t539 + t469) * MDP(18) - t665 * MDP(7)) * t652 + (t370 * t539 + t371 * t540 - t385 * t678 + t386 * t677 - t424 * t504 - t425 * t505 + t469 * t685 + t606 * t684) * MDP(22) + (t424 * t539 - t425 * t540 + t469 * t678 + t606 * t677) * MDP(16) + (t606 * t678 - t406) * MDP(15) + ((t474 * t579 + t481 * t576) * t424 - t361 * t540 + t482 * t390 + t366 * t694 + (t576 * t627 + t579 * t626) * t717 + t680 * t445 + t678 * t364 + t596 * t378) * MDP(32) + (t390 * t540 - t445 * t678 - t539 * t686 + t596 * t717) * MDP(28) + (-t391 * t540 - t419 * t539 + t443 * t678 - t597 * t717) * MDP(29) + (-(-t474 * t576 + t481 * t579) * t424 + t362 * t540 + t482 * t391 - t366 * t579 * t539 + (t576 * t626 - t579 * t627) * t717 + t680 * t443 - t678 * t363 + t597 * t378) * MDP(31) + (-t678 * t717 - t406) * MDP(30) + (t644 * t445 + t643 * t443 + (t702 - t391 * t576 + (-t443 * t579 - t445 * t576) * qJD(6)) * t539) * MDP(27) + (-t519 * t532 + t695 + (-t480 * t580 + t574 * t588 + t699) * t674) * MDP(12) + (t518 * t532 - t525 * t574 + (t479 * t580 + t572 * t588 - t700) * t674) * MDP(11) + (t390 * t694 + t445 * t596) * MDP(26) + (-t370 * t505 + t371 * t504 + t377 * t488 - t385 * t684 + t386 * t685 - t395 * t683) * MDP(25) + (t532 * t564 + t578 * t659 - t525) * MDP(9) + (t531 * t564 + t580 * t659 - t604) * MDP(10) + (t479 * t519 - t480 * t518 + (qJD(3) * t518 + t454 * t651 + t448) * t574 + (qJD(3) * t519 + t455 * t651 - t447) * t572) * MDP(13) - t582 * t637 + t665 * MDP(6) * t651 + (t678 * MDP(17) + t677 * MDP(18) - t712 * MDP(20) + MDP(21) * t718 + t684 * MDP(23) + t685 * MDP(24)) * t555 + t692 * t715; (-t518 ^ 2 - t519 ^ 2) * MDP(13) + (t454 * t519 - t455 * t518 + t525) * MDP(14) + (t425 + t698) * MDP(20) + (-t709 - t721) * MDP(22) + (t583 - t698) * MDP(23) + (t385 * t606 - t386 * t469 + t377) * MDP(25) + (-t579 * t714 + t686 + t720) * MDP(31) + (t576 * t714 + t419 + t719) * MDP(32) + ((qJD(2) * t572 - t519) * MDP(11) + (qJD(2) * t574 - t518) * MDP(12)) * t651 + (-MDP(21) + MDP(24)) * (t424 - t696); -MDP(15) * t697 + (t709 - t721) * MDP(16) + t584 * MDP(17) + (t583 + t698) * MDP(18) + MDP(19) * t632 + (t466 * t606 + t603) * MDP(20) + (t387 * t555 + t466 * t469 + t636) * MDP(21) + (pkin(4) * t424 - qJ(5) * t425 - (-t386 - t388) * t606 + (t385 + t664) * t469) * MDP(22) + (t405 * t469 - t603 - 0.2e1 * t624 - t701) * MDP(23) + (-t395 * t469 - t405 * t606 + t555 * t664 + 0.2e1 * t616 + t623) * MDP(24) + (-pkin(4) * t371 - qJ(5) * t370 - t385 * t388 + t386 * t664 - t395 * t405) * MDP(25) + (-t445 * t642 + t702) * MDP(26) + ((-t391 - t723) * t579 + (-t390 + t724) * t576) * MDP(27) + (-t642 * t717 - t419 + t719) * MDP(28) + (-t641 * t717 + t686 - t720) * MDP(29) + t717 * t469 * MDP(30) + (qJ(5) * t391 + t363 * t469 + t663 * t443 + t600 * t576 + t579 * t710) * MDP(31) + (qJ(5) * t390 - t364 * t469 + t663 * t445 - t576 * t710 + t600 * t579) * MDP(32); t584 * MDP(22) + (t632 + t697) * MDP(23) + (-t555 ^ 2 - t709) * MDP(24) + (-t386 * t555 + t371 - t701) * MDP(25) + (t443 * t555 - t419) * MDP(31) + (t445 * t555 + t686) * MDP(32) + (-MDP(31) * t642 - MDP(32) * t641) * t717; t445 * t443 * MDP(26) + (-t443 ^ 2 + t445 ^ 2) * MDP(27) + (t654 + t724) * MDP(28) + (-t599 + t723) * MDP(29) - t424 * MDP(30) + (t364 * t717 - t378 * t445 + t645) * MDP(31) + (t363 * t717 + t378 * t443 - t620) * MDP(32) + (MDP(28) * t693 - MDP(29) * t445 - MDP(31) * t364 - MDP(32) * t363) * qJD(6);];
tauc  = t1;
