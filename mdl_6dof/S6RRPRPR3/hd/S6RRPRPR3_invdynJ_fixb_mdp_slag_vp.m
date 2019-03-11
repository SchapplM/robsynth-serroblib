% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:27
% EndTime: 2019-03-09 10:19:42
% DurationCPUTime: 11.59s
% Computational Cost: add. (9931->570), mult. (23454->745), div. (0->0), fcn. (18124->16), ass. (0->259)
t644 = sin(pkin(10));
t651 = sin(qJ(2));
t716 = qJD(1) * t651;
t646 = cos(pkin(10));
t655 = cos(qJ(2));
t735 = t646 * t655;
t595 = qJD(1) * t735 - t644 * t716;
t585 = qJD(4) - t595;
t653 = cos(qJ(6));
t610 = t644 * t655 + t646 * t651;
t598 = t610 * qJD(1);
t650 = sin(qJ(4));
t654 = cos(qJ(4));
t711 = t654 * qJD(2);
t564 = t598 * t650 - t711;
t566 = qJD(2) * t650 + t598 * t654;
t643 = sin(pkin(11));
t645 = cos(pkin(11));
t680 = -t564 * t645 - t566 * t643;
t730 = t653 * t680;
t503 = t564 * t643 - t566 * t645;
t649 = sin(qJ(6));
t752 = t503 * t649;
t455 = t730 + t752;
t579 = qJD(6) + t585;
t755 = t455 * t579;
t708 = qJDD(1) * t655;
t709 = qJDD(1) * t651;
t682 = -t644 * t709 + t646 * t708;
t550 = -qJD(2) * t598 + t682;
t710 = qJD(1) * qJD(2);
t700 = t655 * t710;
t701 = t651 * t710;
t551 = qJDD(1) * t610 - t644 * t701 + t646 * t700;
t768 = pkin(2) * t655;
t632 = pkin(1) + t768;
t665 = pkin(2) * t701 - qJDD(1) * t632 + qJDD(3);
t479 = -pkin(3) * t550 - pkin(8) * t551 + t665;
t474 = t654 * t479;
t616 = -qJD(1) * t632 + qJD(3);
t514 = -pkin(3) * t595 - pkin(8) * t598 + t616;
t648 = -qJ(3) - pkin(7);
t617 = t648 * t651;
t614 = qJD(1) * t617;
t758 = qJD(2) * pkin(2);
t604 = t614 + t758;
t618 = t648 * t655;
t615 = qJD(1) * t618;
t736 = t646 * t615;
t549 = t644 * t604 - t736;
t536 = qJD(2) * pkin(8) + t549;
t476 = t514 * t650 + t536 * t654;
t697 = qJD(2) * t648;
t592 = -qJD(3) * t651 + t655 * t697;
t544 = qJDD(2) * pkin(2) + qJD(1) * t592 + qJDD(1) * t617;
t591 = qJD(3) * t655 + t651 * t697;
t554 = qJD(1) * t591 - qJDD(1) * t618;
t491 = t644 * t544 + t646 * t554;
t489 = qJDD(2) * pkin(8) + t491;
t714 = qJD(4) * t650;
t495 = qJD(4) * t711 + t650 * qJDD(2) + t654 * t551 - t598 * t714;
t596 = t610 * qJD(2);
t545 = qJD(1) * t596 + qJDD(4) - t682;
t411 = pkin(4) * t545 - qJ(5) * t495 - qJD(4) * t476 - qJD(5) * t566 - t489 * t650 + t474;
t496 = qJD(4) * t566 - t654 * qJDD(2) + t551 * t650;
t713 = qJD(4) * t654;
t668 = t650 * t479 + t654 * t489 + t514 * t713 - t536 * t714;
t413 = -qJ(5) * t496 - qJD(5) * t564 + t668;
t399 = t645 * t411 - t413 * t643;
t445 = t495 * t645 - t496 * t643;
t397 = pkin(5) * t545 - pkin(9) * t445 + t399;
t400 = t643 * t411 + t645 * t413;
t444 = -t495 * t643 - t496 * t645;
t398 = pkin(9) * t444 + t400;
t475 = t654 * t514 - t536 * t650;
t462 = -qJ(5) * t566 + t475;
t448 = pkin(4) * t585 + t462;
t463 = -qJ(5) * t564 + t476;
t738 = t645 * t463;
t425 = t643 * t448 + t738;
t776 = pkin(9) * t680;
t416 = t425 + t776;
t712 = qJD(6) * t649;
t415 = t416 * t712;
t601 = t644 * t615;
t548 = t604 * t646 + t601;
t535 = -qJD(2) * pkin(3) - t548;
t499 = pkin(4) * t564 + qJD(5) + t535;
t449 = -pkin(5) * t680 + t499;
t637 = qJ(4) + pkin(11) + qJ(6);
t625 = sin(t637);
t626 = cos(t637);
t656 = cos(qJ(1));
t640 = qJ(2) + pkin(10);
t634 = cos(t640);
t652 = sin(qJ(1));
t741 = t634 * t652;
t568 = t625 * t656 - t626 * t741;
t740 = t634 * t656;
t570 = t625 * t652 + t626 * t740;
t633 = sin(t640);
t761 = g(3) * t633;
t789 = g(1) * t570 - g(2) * t568 - t649 * t397 - t653 * t398 - t449 * t455 + t626 * t761 + t415;
t538 = qJDD(6) + t545;
t777 = -t653 * t503 + t649 * t680;
t788 = t538 * MDP(26) + (-t455 ^ 2 + t777 ^ 2) * MDP(23) - t455 * MDP(22) * t777;
t609 = t643 * t654 + t645 * t650;
t780 = t585 * t609;
t677 = t643 * t650 - t645 * t654;
t779 = t585 * t677;
t757 = t777 * t579;
t525 = pkin(2) * t716 + pkin(3) * t598 - pkin(8) * t595;
t520 = t654 * t525;
t556 = t614 * t646 + t601;
t627 = pkin(2) * t644 + pkin(8);
t728 = qJ(5) + t627;
t693 = qJD(4) * t728;
t786 = -pkin(4) * t598 - t520 + (qJ(5) * t595 - t693) * t654 + (-qJD(5) + t556) * t650;
t721 = t650 * t525 + t654 * t556;
t745 = t595 * t650;
t785 = -qJ(5) * t745 - qJD(5) * t654 + t650 * t693 + t721;
t490 = t544 * t646 - t644 * t554;
t488 = -qJDD(2) * pkin(3) - t490;
t690 = g(1) * t656 + g(2) * t652;
t664 = -g(3) * t634 + t633 * t690;
t784 = -qJD(4) * t627 * t585 - t488 + t664;
t783 = t714 - t745;
t567 = t625 * t741 + t626 * t656;
t569 = -t625 * t740 + t626 * t652;
t696 = t653 * t397 - t649 * t398;
t782 = -g(1) * t569 + g(2) * t567 - t449 * t777 + t625 * t761 + t696;
t781 = pkin(9) * t503;
t553 = t609 * t653 - t649 * t677;
t726 = qJD(6) * t553 - t779 * t649 + t653 * t780;
t608 = t644 * t651 - t735;
t600 = t608 * qJD(2);
t703 = t610 * t713;
t778 = -t600 * t650 + t703;
t725 = t643 * t785 + t645 * t786;
t724 = t643 * t786 - t645 * t785;
t773 = g(1) * t652 - g(2) * t656;
t555 = t614 * t644 - t736;
t772 = pkin(4) * t783 - t555;
t729 = t654 * t656;
t734 = t650 * t652;
t586 = t634 * t734 + t729;
t732 = t652 * t654;
t733 = t650 * t656;
t588 = -t634 * t733 + t732;
t771 = -g(1) * t588 + g(2) * t586;
t679 = -t609 * t649 - t653 * t677;
t727 = qJD(6) * t679 - t649 * t780 - t653 * t779;
t770 = -t538 * t553 - t579 * t727;
t695 = -t653 * t444 + t445 * t649;
t408 = qJD(6) * t777 + t695;
t769 = pkin(2) * t646;
t767 = pkin(4) * t643;
t759 = g(3) * t655;
t756 = t455 * t598;
t754 = t777 * t598;
t753 = t495 * t650;
t750 = t545 * t650;
t749 = t564 * t585;
t748 = t564 * t598;
t747 = t566 * t585;
t746 = t566 * t598;
t743 = t610 * t650;
t742 = t610 * t654;
t457 = t643 * t463;
t424 = t645 * t448 - t457;
t414 = pkin(5) * t585 + t424 + t781;
t731 = t653 * t414;
t532 = t654 * t545;
t559 = t617 * t644 - t618 * t646;
t557 = t654 * t559;
t707 = t651 * t758;
t526 = pkin(3) * t596 + pkin(8) * t600 + t707;
t521 = t654 * t526;
t524 = t591 * t646 + t592 * t644;
t547 = pkin(3) * t608 - pkin(8) * t610 - t632;
t676 = qJ(5) * t600 - qJD(5) * t610;
t432 = pkin(4) * t596 - t524 * t650 + t521 + t676 * t654 + (-t557 + (qJ(5) * t610 - t547) * t650) * qJD(4);
t704 = t654 * t524 + t650 * t526 + t547 * t713;
t436 = -qJ(5) * t703 + (-qJD(4) * t559 + t676) * t650 + t704;
t406 = t643 * t432 + t645 * t436;
t428 = t645 * t462 - t457;
t534 = t654 * t547;
t467 = pkin(4) * t608 - qJ(5) * t742 - t559 * t650 + t534;
t720 = t650 * t547 + t557;
t480 = -qJ(5) * t743 + t720;
t438 = t643 * t467 + t645 * t480;
t719 = pkin(5) * t780 + t772;
t605 = t728 * t650;
t606 = t728 * t654;
t540 = -t643 * t605 + t645 * t606;
t641 = t651 ^ 2;
t718 = -t655 ^ 2 + t641;
t715 = qJD(4) * t610;
t705 = qJD(6) * t730 + t649 * t444 + t653 * t445;
t631 = pkin(4) * t654 + pkin(3);
t699 = pkin(4) * t650 - t648;
t405 = t645 * t432 - t436 * t643;
t427 = -t462 * t643 - t738;
t437 = t645 * t467 - t480 * t643;
t523 = t591 * t644 - t646 * t592;
t539 = -t645 * t605 - t606 * t643;
t558 = -t646 * t617 - t618 * t644;
t692 = t585 * t654;
t691 = -qJD(4) * t514 - t489;
t688 = t679 * t538 - t579 * t726;
t686 = pkin(4) * t743 + t558;
t685 = -t536 * t713 + t474;
t510 = -pkin(9) * t677 + t540;
t684 = pkin(5) * t598 - pkin(9) * t779 + qJD(6) * t510 - t725;
t509 = -pkin(9) * t609 + t539;
t683 = pkin(9) * t780 - qJD(6) * t509 - t724;
t402 = t649 * t414 + t653 * t416;
t528 = t609 * t610;
t529 = t677 * t610;
t681 = -t653 * t528 + t529 * t649;
t484 = -t528 * t649 - t529 * t653;
t647 = -qJ(5) - pkin(8);
t678 = t631 * t634 - t633 * t647;
t675 = -t631 - t769;
t674 = pkin(4) * t778 + t523;
t673 = -t585 * t783 + t532;
t628 = pkin(4) * t645 + pkin(5);
t672 = t628 * t649 + t653 * t767;
t671 = t628 * t653 - t649 * t767;
t670 = -0.2e1 * pkin(1) * t710 - pkin(7) * qJDD(2);
t669 = -t600 * t654 - t610 * t714;
t407 = t503 * t712 + t705;
t667 = t535 * t585 - t627 * t545;
t657 = qJD(2) ^ 2;
t661 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t657 + t773;
t658 = qJD(1) ^ 2;
t660 = pkin(1) * t658 - pkin(7) * qJDD(1) + t690;
t440 = pkin(4) * t496 + qJDD(5) + t488;
t629 = -pkin(3) - t769;
t620 = t656 * t632;
t589 = t634 * t729 + t734;
t587 = -t634 * t732 + t733;
t563 = pkin(5) * t677 + t675;
t486 = pkin(5) * t528 + t686;
t485 = pkin(4) * t566 - pkin(5) * t503;
t482 = -t600 * t677 + t609 * t715;
t481 = t600 * t609 + t677 * t715;
t439 = -pkin(5) * t481 + t674;
t431 = -pkin(9) * t528 + t438;
t426 = pkin(5) * t608 + pkin(9) * t529 + t437;
t423 = qJD(6) * t484 - t653 * t481 - t482 * t649;
t422 = qJD(6) * t681 + t481 * t649 - t482 * t653;
t419 = t428 + t781;
t418 = t427 - t776;
t417 = -pkin(5) * t444 + t440;
t404 = pkin(9) * t481 + t406;
t403 = pkin(5) * t596 + pkin(9) * t482 + t405;
t401 = -t416 * t649 + t731;
t1 = [(t399 * t529 - t400 * t528 + t405 * t503 + t406 * t680 + t424 * t482 + t425 * t481 - t437 * t445 + t438 * t444 + t633 * t773) * MDP(20) + t773 * MDP(2) + (-t490 * t610 - t491 * t608 + t523 * t598 + t524 * t595 + t548 * t600 - t549 * t596 + t550 * t559 + t551 * t558 - t690) * MDP(11) + ((-t559 * t713 + t521) * t585 + t534 * t545 + t685 * t608 + t475 * t596 + t523 * t564 + t558 * t496 + t535 * t703 - g(1) * t587 - g(2) * t589 + ((-qJD(4) * t547 - t524) * t585 - t559 * t545 + t691 * t608 + t488 * t610 - t535 * t600) * t650) * MDP(18) + (-(-t564 * t654 - t566 * t650) * t600 + (-t753 - t496 * t654 + (t564 * t650 - t566 * t654) * qJD(4)) * t610) * MDP(14) + (t407 * t608 + t422 * t579 + t484 * t538 + t596 * t777) * MDP(24) + (-g(1) * t567 - g(2) * t569 - t402 * t596 + t486 * t407 + t415 * t608 + t417 * t484 + t449 * t422 + t439 * t777 + (-(-qJD(6) * t431 + t403) * t579 - t426 * t538 - t397 * t608) * t649 + (-(qJD(6) * t426 + t404) * t579 - t431 * t538 - (qJD(6) * t414 + t398) * t608) * t653) * MDP(28) + (t407 * t484 + t422 * t777) * MDP(22) + t690 * MDP(3) + (-t408 * t608 - t423 * t579 + t455 * t596 + t538 * t681) * MDP(25) + ((t403 * t653 - t404 * t649) * t579 + (t426 * t653 - t431 * t649) * t538 + t696 * t608 + t401 * t596 - t439 * t455 + t486 * t408 - t417 * t681 + t449 * t423 - g(1) * t568 - g(2) * t570 + ((-t426 * t649 - t431 * t653) * t579 - t402 * t608) * qJD(6)) * MDP(27) + (t407 * t681 - t408 * t484 + t422 * t455 - t423 * t777) * MDP(23) + (t538 * t608 + t579 * t596) * MDP(26) + (t545 * t608 + t585 * t596) * MDP(17) + (qJDD(2) * t651 + t655 * t657) * MDP(6) + (qJDD(2) * t655 - t651 * t657) * MDP(7) + (t651 * t670 + t655 * t661) * MDP(9) + (-t651 * t661 + t655 * t670) * MDP(10) + qJDD(1) * MDP(1) + (t400 * t438 + t425 * t406 + t399 * t437 + t424 * t405 + t440 * t686 + t499 * t674 - g(2) * t620 + (-g(1) * t699 - g(2) * t678) * t656 + (-g(1) * (-t632 - t678) - g(2) * t699) * t652) * MDP(21) + (qJDD(1) * t641 + 0.2e1 * t651 * t700) * MDP(4) + (t491 * t559 + t549 * t524 - t490 * t558 - t548 * t523 - t665 * t632 + t616 * t707 - g(1) * (-t632 * t652 - t648 * t656) - g(2) * (-t648 * t652 + t620)) * MDP(12) + (-t496 * t608 - t545 * t743 - t564 * t596 - t585 * t778) * MDP(16) + 0.2e1 * (t651 * t708 - t710 * t718) * MDP(5) + (t495 * t608 + t532 * t610 + t566 * t596 + t585 * t669) * MDP(15) + (t495 * t742 + t566 * t669) * MDP(13) + (-(-t559 * t714 + t704) * t585 - t720 * t545 - t668 * t608 - t476 * t596 + t523 * t566 + t558 * t495 + t488 * t742 - g(1) * t586 - g(2) * t588 + t669 * t535) * MDP(19); (t688 - t756) * MDP(25) + (-t399 * t609 - t400 * t677 + t424 * t779 - t425 * t780 + t444 * t540 - t445 * t539 + t503 * t725 - t634 * t690 + t680 * t724 - t761) * MDP(20) + (t585 * t692 - t746 + t750) * MDP(15) + (t673 + t748) * MDP(16) + ((t509 * t653 - t510 * t649) * t538 + t563 * t408 - t417 * t679 + (t649 * t683 - t653 * t684) * t579 - t719 * t455 + t726 * t449 + t664 * t626) * MDP(27) + ((t548 - t556) * t595 + (t550 * t644 - t551 * t646) * pkin(2)) * MDP(11) + (-(t509 * t649 + t510 * t653) * t538 + t563 * t407 + t417 * t553 + (t649 * t684 + t653 * t683) * t579 + t719 * t777 + t727 * t449 - t664 * t625) * MDP(28) + (t407 * t679 - t408 * t553 + t455 * t727 - t726 * t777) * MDP(23) + (-t754 - t770) * MDP(24) + (g(3) * t651 + t655 * t660) * MDP(10) + MDP(7) * t708 + MDP(6) * t709 + (t629 * t495 - t555 * t566 + t721 * t585 - t650 * t784 + t667 * t654) * MDP(19) + (t629 * t496 - t520 * t585 - t555 * t564 + (t556 * t585 + t667) * t650 + t784 * t654) * MDP(18) + qJDD(2) * MDP(8) + (t407 * t553 + t727 * t777) * MDP(22) + ((t495 - t749) * t654 + (-t496 - t747) * t650) * MDP(14) + (t566 * t692 + t753) * MDP(13) + (t548 * t555 - t549 * t556 + (-t759 + t490 * t646 + t491 * t644 + (-qJD(1) * t616 + t690) * t651) * pkin(2)) * MDP(12) + (t651 * t660 - t759) * MDP(9) + (t400 * t540 + t399 * t539 + t440 * t675 - g(3) * (t678 + t768) + t772 * t499 + t724 * t425 + t725 * t424 + t690 * (pkin(2) * t651 + t631 * t633 + t634 * t647)) * MDP(21) + (-MDP(4) * t651 * t655 + MDP(5) * t718) * t658 + (-t401 * MDP(27) + (t549 - t555) * MDP(11) + t402 * MDP(28) + t476 * MDP(19) - t475 * MDP(18) - t585 * MDP(17) - t579 * MDP(26)) * t598; (-t595 ^ 2 - t598 ^ 2) * MDP(11) + (t548 * t598 - t549 * t595 + t665 - t773) * MDP(12) + (t673 - t748) * MDP(18) + (-t585 ^ 2 * t654 - t746 - t750) * MDP(19) + (t444 * t609 + t445 * t677 - t503 * t780 - t680 * t779) * MDP(20) + (-t399 * t677 + t400 * t609 - t424 * t780 - t425 * t779 - t499 * t598 - t773) * MDP(21) + (t688 + t756) * MDP(27) + (-t754 + t770) * MDP(28); t566 * t564 * MDP(13) + (-t564 ^ 2 + t566 ^ 2) * MDP(14) + (t495 + t749) * MDP(15) + (-t496 + t747) * MDP(16) + t545 * MDP(17) + (t476 * t585 - t535 * t566 + (t691 + t761) * t650 + t685 + t771) * MDP(18) + (g(1) * t589 - g(2) * t587 + t475 * t585 + t535 * t564 + t654 * t761 - t668) * MDP(19) + ((t444 * t643 - t445 * t645) * pkin(4) + (t424 - t428) * t680 + (-t425 - t427) * t503) * MDP(20) + (-t424 * t427 - t425 * t428 + (t399 * t645 + t400 * t643 - t499 * t566 + t650 * t761 + t771) * pkin(4)) * MDP(21) + (t407 - t755) * MDP(24) + (-t408 + t757) * MDP(25) + (t671 * t538 - (t418 * t653 - t419 * t649) * t579 + t485 * t455 + (-t579 * t672 - t402) * qJD(6) + t782) * MDP(27) + (-t672 * t538 + (t418 * t649 + t419 * t653) * t579 - t485 * t777 + (-t579 * t671 - t731) * qJD(6) + t789) * MDP(28) + t788; (-t503 ^ 2 - t680 ^ 2) * MDP(20) + (-t424 * t503 - t425 * t680 + t440 - t664) * MDP(21) + (t408 + t757) * MDP(27) + (t407 + t755) * MDP(28); (t705 - t755) * MDP(24) + (-t695 + t757) * MDP(25) + (t402 * t579 + t782) * MDP(27) + (t401 * t579 + t789) * MDP(28) + (MDP(24) * t752 - MDP(25) * t777 - MDP(27) * t402 - MDP(28) * t731) * qJD(6) + t788;];
tau  = t1;
