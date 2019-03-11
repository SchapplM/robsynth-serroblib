% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:31
% EndTime: 2019-03-09 12:25:46
% DurationCPUTime: 9.51s
% Computational Cost: add. (10173->545), mult. (25022->699), div. (0->0), fcn. (18497->8), ass. (0->232)
t582 = sin(pkin(10));
t585 = sin(qJ(4));
t583 = cos(pkin(10));
t587 = cos(qJ(4));
t690 = t583 * t587;
t542 = t582 * t585 - t690;
t588 = cos(qJ(2));
t605 = t542 * t588;
t672 = qJD(1) * t605 - t542 * qJD(4);
t586 = sin(qJ(2));
t667 = qJD(1) * t586;
t645 = t582 * t667;
t656 = t583 * qJD(2);
t535 = t645 - t656;
t644 = t583 * t667;
t663 = qJD(2) * t582;
t537 = t644 + t663;
t479 = t535 * t585 - t537 * t587;
t584 = sin(qJ(5));
t619 = -t535 * t587 - t585 * t537;
t699 = cos(qJ(5));
t429 = t584 * t479 + t619 * t699;
t735 = t429 ^ 2;
t666 = qJD(1) * t588;
t571 = -qJD(4) + t666;
t562 = -qJD(5) + t571;
t734 = t429 * t562;
t692 = t583 * t585;
t543 = t582 * t587 + t692;
t713 = t543 * t588;
t504 = qJD(1) * t713;
t600 = t543 * qJD(4);
t705 = t504 - t600;
t621 = pkin(2) * t586 - qJ(3) * t588;
t545 = t621 * qJD(1);
t495 = pkin(7) * t645 + t583 * t545;
t689 = t583 * t588;
t617 = pkin(3) * t586 - pkin(8) * t689;
t469 = qJD(1) * t617 + t495;
t526 = t582 * t545;
t691 = t583 * t586;
t693 = t582 * t588;
t610 = -pkin(7) * t691 - pkin(8) * t693;
t485 = qJD(1) * t610 + t526;
t697 = pkin(8) + qJ(3);
t554 = t697 * t582;
t555 = t697 * t583;
t670 = -t585 * t554 + t587 * t555;
t733 = t543 * qJD(3) + qJD(4) * t670 + t587 * t469 - t485 * t585;
t658 = qJD(4) * t587;
t732 = qJD(3) * t690 - t585 * t469 - t587 * t485 - t554 * t658;
t721 = -t479 * t699 + t584 * t619;
t700 = t721 ^ 2;
t731 = t429 * t721;
t730 = t562 * t721;
t729 = pkin(4) * t667 + pkin(9) * t672 + t733;
t660 = qJD(3) * t582;
t686 = t585 * t555;
t698 = pkin(9) * t543;
t728 = pkin(9) * t504 - t585 * t660 + (-t686 - t698) * qJD(4) + t732;
t651 = qJD(1) * qJD(2);
t641 = t588 * t651;
t624 = t587 * t641;
t625 = t582 * t641;
t442 = -t535 * t658 + t583 * t624 + (-qJD(4) * t537 - t625) * t585;
t659 = qJD(4) * t585;
t626 = -t535 * t659 + t537 * t658 + t582 * t624 + t641 * t692;
t643 = qJD(5) * t699;
t657 = qJD(5) * t584;
t391 = -t699 * t442 - t479 * t657 + t584 * t626 - t619 * t643;
t380 = -t391 + t734;
t392 = qJD(5) * t721 + t584 * t442 + t699 * t626;
t642 = t586 * t651;
t727 = MDP(26) * t642 + (-t392 - t730) * MDP(25) - MDP(22) * t731 + (t700 - t735) * MDP(23) + t380 * MDP(24);
t714 = t542 * t586;
t726 = pkin(9) * t714;
t575 = pkin(7) * t667;
t696 = qJD(2) * pkin(2);
t638 = qJD(3) - t696;
t551 = t575 + t638;
t494 = pkin(3) * t535 + t551;
t439 = -pkin(4) * t619 + t494;
t386 = -pkin(5) * t429 - qJ(6) * t721 + t439;
t725 = t386 * t429;
t724 = t429 * t439;
t723 = t479 * t571;
t676 = t542 * t643 + t543 * t657 - t584 * t705 - t672 * t699;
t484 = -t584 * t542 + t543 * t699;
t675 = qJD(5) * t484 + t584 * t672 - t699 * t705;
t399 = pkin(5) * t721 - qJ(6) * t429;
t719 = -0.2e1 * t651;
t717 = MDP(4) * t586;
t716 = MDP(5) * (t586 ^ 2 - t588 ^ 2);
t695 = t386 * t721;
t715 = t439 * t721;
t712 = t571 * t619;
t632 = -t587 * t554 - t686;
t463 = t632 - t698;
t464 = -pkin(9) * t542 + t670;
t612 = t463 * t699 - t584 * t464;
t709 = qJD(5) * t612 - t584 * t729 + t699 * t728;
t421 = t584 * t463 + t464 * t699;
t708 = qJD(5) * t421 + t584 * t728 + t699 * t729;
t552 = -pkin(2) * t588 - qJ(3) * t586 - pkin(1);
t534 = t583 * t552;
t486 = -pkin(8) * t691 + t534 + (-pkin(7) * t582 - pkin(3)) * t588;
t569 = pkin(7) * t689;
t502 = t582 * t552 + t569;
t694 = t582 * t586;
t493 = -pkin(8) * t694 + t502;
t687 = t585 * t493;
t633 = t587 * t486 - t687;
t417 = -pkin(4) * t588 + t633 + t726;
t515 = t543 * t586;
t673 = t585 * t486 + t587 * t493;
t422 = -pkin(9) * t515 + t673;
t707 = t584 * t417 + t699 * t422;
t576 = pkin(7) * t666;
t649 = pkin(3) * t666;
t529 = t582 * t649 + t576;
t706 = -pkin(4) * t705 - t529;
t465 = -qJD(2) * t605 - t586 * t600;
t522 = qJD(2) * t621 - qJD(3) * t586;
t662 = qJD(2) * t586;
t648 = pkin(7) * t662;
t491 = t583 * t522 + t582 * t648;
t599 = t617 * qJD(2);
t460 = t599 + t491;
t508 = t582 * t522;
t471 = qJD(2) * t610 + t508;
t635 = t587 * t460 - t471 * t585;
t388 = pkin(4) * t662 - pkin(9) * t465 - qJD(4) * t673 + t635;
t597 = qJD(2) * t713;
t647 = t585 * t460 + t587 * t471 + t486 * t658;
t393 = -pkin(9) * t597 + (-t687 + t726) * qJD(4) + t647;
t702 = -qJD(5) * t707 + t388 * t699 - t584 * t393;
t528 = t552 * qJD(1);
t557 = qJD(2) * qJ(3) + t576;
t487 = t583 * t528 - t557 * t582;
t448 = -pkin(8) * t537 + t487 - t649;
t488 = t582 * t528 + t583 * t557;
t452 = -pkin(8) * t535 + t488;
t412 = t448 * t585 + t452 * t587;
t404 = pkin(9) * t619 + t412;
t688 = t584 * t404;
t589 = qJD(2) ^ 2;
t685 = t586 * t589;
t684 = t588 * t589;
t590 = qJD(1) ^ 2;
t683 = t588 * t590;
t411 = t587 * t448 - t452 * t585;
t403 = pkin(9) * t479 + t411;
t376 = t403 * t699 - t688;
t682 = -pkin(4) * t643 - qJD(6) + t376;
t681 = -qJ(6) * t667 + t709;
t680 = pkin(5) * t667 + t708;
t679 = pkin(5) * t675 + qJ(6) * t676 - t484 * qJD(6) + t706;
t506 = t522 * qJD(1);
t550 = (qJD(3) - t575) * qJD(2);
t468 = t582 * t506 + t583 * t550;
t570 = pkin(7) * t641;
t521 = pkin(3) * t625 + t570;
t661 = qJD(2) * t588;
t577 = pkin(7) * t661;
t530 = t582 * pkin(3) * t661 + t577;
t546 = pkin(3) * t694 + t586 * pkin(7);
t665 = qJD(2) * t612;
t664 = qJD(2) * t421;
t655 = t586 * MDP(19);
t400 = -pkin(4) * t571 + t403;
t371 = t400 * t699 - t688;
t654 = qJD(6) - t371;
t653 = MDP(11) * qJD(1);
t652 = MDP(12) * qJD(1);
t650 = pkin(7) * t693;
t573 = -pkin(3) * t583 - pkin(2);
t646 = t699 * t404;
t640 = qJD(2) * t655;
t639 = pkin(1) * t719;
t467 = t583 * t506 - t550 * t582;
t445 = qJD(1) * t599 + t467;
t449 = -pkin(8) * t625 + t468;
t636 = t587 * t445 - t449 * t585;
t631 = t535 + t656;
t630 = -t537 + t663;
t629 = pkin(5) * t642;
t594 = -qJD(4) * t412 + t636;
t377 = pkin(4) * t642 - pkin(9) * t442 + t594;
t603 = t585 * t445 + t448 * t658 + t587 * t449 - t452 * t659;
t381 = -pkin(9) * t626 + t603;
t628 = -t584 * t377 - t699 * t381 - t400 * t643 + t404 * t657;
t627 = -t699 * t377 + t584 * t381 + t400 * t657 + t404 * t643;
t490 = pkin(4) * t515 + t546;
t375 = t584 * t403 + t646;
t623 = pkin(4) * t657 - t375;
t622 = -t551 + t638;
t507 = pkin(4) * t542 + t573;
t553 = t562 * qJD(6);
t567 = qJ(6) * t642;
t363 = t567 - t553 - t628;
t615 = t417 * t699 - t584 * t422;
t372 = t584 * t400 + t646;
t459 = -t584 * t515 - t699 * t714;
t609 = -t371 * t562 + t628;
t608 = -t372 * t562 - t627;
t602 = t584 * t388 + t699 * t393 + t417 * t643 - t422 * t657;
t596 = qJD(4) * t714;
t423 = pkin(4) * t626 + t521;
t364 = t627 - t629;
t367 = t392 * pkin(5) + t391 * qJ(6) - qJD(6) * t721 + t423;
t591 = t596 - t597;
t443 = -pkin(4) * t591 + t530;
t574 = -pkin(4) * t699 - pkin(5);
t572 = pkin(4) * t584 + qJ(6);
t501 = t534 - t650;
t496 = -pkin(7) * t644 + t526;
t492 = -t583 * t648 + t508;
t483 = t542 * t699 + t543 * t584;
t458 = t515 * t699 - t584 * t714;
t419 = pkin(5) * t483 - qJ(6) * t484 + t507;
t408 = pkin(5) * t458 - qJ(6) * t459 + t490;
t407 = qJD(5) * t459 + t584 * t465 - t591 * t699;
t406 = -t465 * t699 + t515 * t643 - t584 * t591 - t657 * t714;
t396 = -pkin(4) * t479 + t399;
t395 = t588 * pkin(5) - t615;
t394 = -qJ(6) * t588 + t707;
t370 = -t562 * qJ(6) + t372;
t369 = t562 * pkin(5) + t654;
t368 = t407 * pkin(5) + t406 * qJ(6) - t459 * qJD(6) + t443;
t366 = -pkin(5) * t662 - t702;
t365 = qJ(6) * t662 - qJD(6) * t588 + t602;
t1 = [(t364 * t588 + t366 * t562 + t367 * t458 - t368 * t429 + t386 * t407 + t392 * t408 + (-qJD(1) * t395 - t369) * t662) * MDP(29) + (t371 * t662 + t490 * t392 + t439 * t407 + t423 * t458 - t429 * t443 - t562 * t702 + t627 * t588 + t615 * t642) * MDP(27) + (-t363 * t458 + t364 * t459 + t365 * t429 + t366 * t721 - t369 * t406 - t370 * t407 - t391 * t395 - t392 * t394) * MDP(30) + (t391 * t458 - t392 * t459 - t406 * t429 - t407 * t721) * MDP(23) + (t392 * t588 + t407 * t562 + (-qJD(1) * t458 + t429) * t662) * MDP(25) + 0.2e1 * t641 * t717 + t716 * t719 + (t391 * t588 + t406 * t562 + (qJD(1) * t459 + t721) * t662) * MDP(24) + (-t363 * t588 - t365 * t562 - t367 * t459 - t368 * t721 + t386 * t406 + t391 * t408 + (qJD(1) * t394 + t370) * t662) * MDP(31) + (-t391 * t459 - t406 * t721) * MDP(22) + (t602 * t562 - t628 * t588 + t443 * t721 - t490 * t391 + t423 * t459 - t439 * t406 + (-qJD(1) * t707 - t372) * t662) * MDP(28) + (-t635 * t571 - t636 * t588 - t530 * t619 + t546 * t626 + t521 * t515 + (t412 * t588 - t494 * t714 + t571 * t673) * qJD(4) + (t494 * t713 + (qJD(1) * t633 + t411) * t586) * qJD(2)) * MDP(20) - MDP(7) * t685 + (pkin(7) * t685 + t588 * t639) * MDP(10) + (-pkin(7) * t684 + t586 * t639) * MDP(9) + ((qJD(1) * t492 + t468) * t588 + ((pkin(7) * t537 + t551 * t583) * t588 + (-t488 + (-t502 + 0.2e1 * t569) * qJD(1)) * t586) * qJD(2)) * MDP(12) + ((-qJD(1) * t491 - t467) * t588 + ((pkin(7) * t535 + t551 * t582) * t588 + (t487 + (t501 + 0.2e1 * t650) * qJD(1)) * t586) * qJD(2)) * MDP(11) + (t467 * t501 + t468 * t502 + t487 * t491 + t488 * t492 + (t551 + t575) * t577) * MDP(14) + (-t571 - t666) * t640 + (-t491 * t537 - t492 * t535 + (-t467 * t583 - t468 * t582) * t586 + (-t487 * t583 - t488 * t582 + (-t501 * t583 - t502 * t582) * qJD(1)) * t661) * MDP(13) + (-t442 * t515 + t465 * t619 - t479 * t591 + t626 * t714) * MDP(16) + ((-t493 * t659 + t647) * t571 + t603 * t588 - t530 * t479 + t546 * t442 - t521 * t714 + t494 * t465 + (-qJD(1) * t673 - t412) * t662) * MDP(21) + (-t442 * t588 - t465 * t571 + (-qJD(1) * t714 - t479) * t662) * MDP(17) + (-t442 * t714 - t465 * t479) * MDP(15) + (t626 * t588 - t571 * t596 + (t571 * t713 + (-t515 * qJD(1) + t619) * t586) * qJD(2)) * MDP(18) + (t363 * t394 + t364 * t395 + t365 * t370 + t366 * t369 + t367 * t408 + t368 * t386) * MDP(32) + MDP(6) * t684 + (-t562 - t666) * MDP(26) * t662; (-t363 * t483 + t364 * t484 - t369 * t676 - t370 * t675 + t391 * t612 - t392 * t421 + t429 * t681 + t680 * t721) * MDP(30) + (t391 * t483 - t392 * t484 - t429 * t676 - t675 * t721) * MDP(23) + (t367 * t483 + t392 * t419 + t680 * t562 - t679 * t429 + t675 * t386 + (t369 + t665) * t667) * MDP(29) + (t675 * t562 + (-qJD(2) * t483 - t429) * t667) * MDP(25) + (t507 * t392 + t423 * t483 + t708 * t562 + t675 * t439 + (-t371 + t665) * t667 - t706 * t429) * MDP(27) - t683 * t717 + (MDP(9) * t586 * t590 + MDP(10) * t683) * pkin(1) + t590 * t716 + (-t391 * t484 - t676 * t721) * MDP(22) + (-t507 * t391 + t423 * t484 + t709 * t562 - t676 * t439 + (t372 - t664) * t667 + t706 * t721) * MDP(28) + (-t705 * t571 + (-qJD(2) * t542 - t619) * t667) * MDP(18) + (-t367 * t484 + t391 * t419 - t681 * t562 - t679 * t721 + t676 * t386 + (-t370 + t664) * t667) * MDP(31) + (t676 * t562 + (qJD(2) * t484 - t721) * t667) * MDP(24) + t571 * qJD(1) * t655 + (-t487 * t495 - t488 * t496 + (-t487 * t582 + t488 * t583) * qJD(3) + (-t467 * t582 + t468 * t583) * qJ(3) + (-t551 - t696) * t576) * MDP(14) + (t573 * t442 + t529 * t479 + t521 * t543 + ((-qJD(4) * t555 - t660) * t585 + t732) * t571 + t672 * t494 + (-qJD(2) * t670 + t412) * t667) * MDP(21) + (t573 * t626 + t521 * t542 + t529 * t619 + (qJD(2) * t632 - t411) * t667 - t705 * t494 + t733 * t571) * MDP(20) + ((-qJ(3) * t663 - t487) * t586 + (-pkin(7) * t631 + t582 * t622 + t495) * t588) * t653 + (t495 * t537 + t496 * t535 + (-qJD(3) * t535 + t487 * t666 + t468) * t583 + (qJD(3) * t537 + t488 * t666 - t467) * t582) * MDP(13) + ((-qJ(3) * t656 + t488) * t586 + (pkin(7) * t630 + t583 * t622 - t496) * t588) * t652 + (t363 * t421 - t364 * t612 + t367 * t419 + t369 * t680 + t370 * t681 + t386 * t679) * MDP(32) + (t442 * t543 - t479 * t672) * MDP(15) + (-t672 * t571 + (qJD(2) * t543 + t479) * t667) * MDP(17) + (-t442 * t542 - t479 * t705 - t543 * t626 + t619 * t672) * MDP(16) + t562 * MDP(26) * t667; (-t535 ^ 2 - t537 ^ 2) * MDP(13) + (t487 * t537 + t488 * t535 + t570) * MDP(14) + (t626 + t723) * MDP(20) + (t442 - t712) * MDP(21) + (-t700 - t735) * MDP(30) + (-t369 * t721 - t370 * t429 + t367) * MDP(32) + (-MDP(28) + MDP(31)) * (t391 + t734) + (MDP(27) + MDP(29)) * (t392 - t730) + (t630 * t653 + t631 * t652) * t588; t479 * t619 * MDP(15) + (t479 ^ 2 - t619 ^ 2) * MDP(16) + (t442 + t712) * MDP(17) + (-t626 + t723) * MDP(18) + qJD(1) * t640 + (-t412 * t571 + t479 * t494 + t594) * MDP(20) + (-t411 * t571 - t494 * t619 - t603) * MDP(21) + (-t375 * t562 - t715 + (-t429 * t479 + t562 * t657 + t642 * t699) * pkin(4) - t627) * MDP(27) + (-t376 * t562 - t724 + (t479 * t721 + t562 * t643 - t584 * t642) * pkin(4) + t628) * MDP(28) + (-t695 + t396 * t429 + t623 * t562 + (pkin(5) - t574) * t642 - t627) * MDP(29) + (-t391 * t574 - t392 * t572 + (t370 + t623) * t721 - (t369 + t682) * t429) * MDP(30) + (t396 * t721 + t562 * t682 + t572 * t642 + t363 + t725) * MDP(31) + (t363 * t572 + t364 * t574 + t369 * t623 - t370 * t682 - t386 * t396) * MDP(32) + t727; (t608 - t715) * MDP(27) + (t609 - t724) * MDP(28) + (t399 * t429 + t608 + 0.2e1 * t629 - t695) * MDP(29) + (pkin(5) * t391 - qJ(6) * t392 + (t370 - t372) * t721 - (t369 - t654) * t429) * MDP(30) + (t399 * t721 - 0.2e1 * t553 + 0.2e1 * t567 - t609 + t725) * MDP(31) + (-pkin(5) * t364 + qJ(6) * t363 - t369 * t372 + t370 * t654 - t386 * t399) * MDP(32) + t727; (-t642 - t731) * MDP(29) + t380 * MDP(30) + (-t562 ^ 2 - t700) * MDP(31) + (t370 * t562 + t364 + t695) * MDP(32);];
tauc  = t1;
