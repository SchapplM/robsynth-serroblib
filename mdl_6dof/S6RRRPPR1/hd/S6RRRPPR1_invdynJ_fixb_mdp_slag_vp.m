% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:25
% EndTime: 2019-03-09 15:23:41
% DurationCPUTime: 11.85s
% Computational Cost: add. (11874->586), mult. (28536->760), div. (0->0), fcn. (21729->18), ass. (0->285)
t680 = qJD(2) + qJD(3);
t688 = sin(qJ(3));
t692 = cos(qJ(3));
t689 = sin(qJ(2));
t779 = qJD(1) * t689;
t760 = t688 * t779;
t693 = cos(qJ(2));
t769 = qJDD(1) * t693;
t772 = qJD(1) * qJD(2);
t758 = t693 * t772;
t770 = qJDD(1) * t689;
t720 = -t758 - t770;
t771 = qJD(1) * qJD(3);
t829 = -t693 * t771 + t720;
t547 = -t680 * t760 + t688 * t769 - t692 * t829;
t778 = qJD(1) * t693;
t599 = -t688 * t778 - t692 * t779;
t677 = qJDD(2) + qJDD(3);
t823 = pkin(7) + pkin(8);
t572 = qJDD(2) * pkin(2) + t720 * t823;
t759 = t689 * t772;
t719 = -t759 + t769;
t576 = t823 * t719;
t636 = t823 * t693;
t626 = qJD(1) * t636;
t608 = t692 * t626;
t635 = t823 * t689;
t624 = qJD(1) * t635;
t810 = qJD(2) * pkin(2);
t611 = -t624 + t810;
t728 = -t611 * t688 - t608;
t710 = qJD(3) * t728 + t692 * t572 - t688 * t576;
t478 = pkin(3) * t677 - qJ(4) * t547 + qJD(4) * t599 + t710;
t777 = qJD(3) * t688;
t594 = t626 * t777;
t598 = t692 * t778 - t760;
t744 = -qJD(3) * t611 - t576;
t482 = t598 * qJD(4) - t594 + (qJ(4) * t829 + t572) * t688 + ((-t689 * t771 + t719) * qJ(4) - t744) * t692;
t685 = sin(pkin(10));
t809 = cos(pkin(10));
t443 = t478 * t809 - t685 * t482;
t442 = -t677 * pkin(4) + qJDD(5) - t443;
t683 = qJ(2) + qJ(3);
t671 = pkin(10) + t683;
t658 = sin(t671);
t690 = sin(qJ(1));
t694 = cos(qJ(1));
t740 = g(1) * t694 + g(2) * t690;
t838 = t740 * t658;
t841 = t442 - t838;
t564 = t809 * t598 + t599 * t685;
t559 = qJD(6) - t564;
t684 = sin(pkin(11));
t686 = cos(pkin(11));
t721 = t685 * t598 - t599 * t809;
t548 = -t686 * t680 + t684 * t721;
t691 = cos(qJ(6));
t550 = t680 * t684 + t686 * t721;
t687 = sin(qJ(6));
t801 = t550 * t687;
t837 = -t691 * t548 - t801;
t840 = t559 * t837;
t659 = cos(t671);
t839 = -g(3) * t659 - t841;
t617 = t688 * t689 - t692 * t693;
t675 = t693 * pkin(2);
t811 = pkin(1) + t675;
t743 = pkin(3) * t617 - t811;
t672 = sin(t683);
t673 = cos(t683);
t836 = -g(3) * t673 + t672 * t740;
t634 = qJD(1) * t811;
t835 = qJDD(1) * t811;
t729 = t548 * t687 - t550 * t691;
t834 = t559 * t729;
t618 = t688 * t693 + t689 * t692;
t571 = t680 * t618;
t614 = t684 * t691 + t686 * t687;
t597 = t614 * qJD(6);
t831 = t614 * t564 - t597;
t613 = t684 * t687 - t691 * t686;
t830 = t559 * t613;
t808 = qJ(4) * t598;
t546 = -t728 + t808;
t539 = t685 * t546;
t592 = t599 * qJ(4);
t604 = t688 * t626;
t747 = t692 * t611 - t604;
t545 = t592 + t747;
t502 = t545 * t809 - t539;
t822 = pkin(3) * t599;
t520 = pkin(4) * t721 - qJ(5) * t564 - t822;
t467 = t686 * t502 + t684 * t520;
t827 = -qJD(5) * t686 + t467;
t746 = t624 * t688 - t608;
t553 = t746 - t808;
t783 = -t692 * t624 - t604;
t554 = t592 + t783;
t516 = t685 * t553 + t554 * t809;
t667 = pkin(2) * t779;
t517 = t520 + t667;
t468 = -t516 * t684 + t686 * t517;
t652 = t685 * t688 * pkin(2);
t776 = qJD(3) * t692;
t589 = t809 * pkin(2) * t776 - qJD(3) * t652;
t582 = qJD(5) + t589;
t826 = -t582 * t684 - t468;
t469 = t686 * t516 + t684 * t517;
t749 = t582 * t686 - t469;
t751 = t809 * t688;
t785 = t809 * t553 - t554 * t685 + (t685 * t692 + t751) * qJD(3) * pkin(2);
t782 = -t688 * t635 + t692 * t636;
t737 = t659 * pkin(4) + t658 * qJ(5);
t825 = g(1) * t690 - g(2) * t694;
t701 = qJD(1) * t571;
t698 = -t617 * qJDD(1) - t701;
t655 = pkin(2) * t759;
t697 = pkin(3) * t701 + qJDD(1) * t743 + qJDD(4) + t655;
t503 = t547 * t685 - t809 * t698;
t500 = qJDD(6) + t503;
t824 = t500 * t614 - t559 * t830;
t820 = pkin(3) * t672;
t819 = pkin(4) * t658;
t679 = pkin(11) + qJ(6);
t669 = sin(t679);
t815 = g(3) * t669;
t670 = cos(t679);
t814 = g(3) * t670;
t812 = t686 * pkin(5);
t674 = t686 * pkin(9);
t807 = qJ(5) * t659;
t444 = t685 * t478 + t809 * t482;
t438 = qJ(5) * t677 + qJD(5) * t680 + t444;
t504 = t809 * t547 + t685 * t698;
t447 = t503 * pkin(4) - t504 * qJ(5) - qJD(5) * t721 + t697;
t424 = t686 * t438 + t684 * t447;
t422 = t424 * t686;
t537 = pkin(3) * t680 + t545;
t495 = t537 * t809 - t539;
t491 = -t680 * pkin(4) + qJD(5) - t495;
t806 = t491 * t564;
t804 = t837 * t721;
t803 = t729 * t721;
t570 = t680 * t617;
t534 = -t570 * t809 - t685 * t571;
t802 = t534 * t684;
t800 = t564 * t684;
t799 = t564 * t686;
t568 = -t685 * t617 + t618 * t809;
t798 = t568 * t684;
t797 = t568 * t686;
t665 = pkin(2) * t692 + pkin(3);
t591 = pkin(2) * t751 + t685 * t665;
t585 = qJ(5) + t591;
t794 = t585 * t686;
t656 = pkin(3) * t685 + qJ(5);
t793 = t656 * t686;
t792 = t659 * t684;
t791 = t669 * t690;
t790 = t669 * t694;
t789 = t670 * t690;
t788 = t670 * t694;
t761 = qJD(2) * t823;
t625 = t689 * t761;
t627 = t693 * t761;
t718 = -t692 * t625 - t688 * t627 - t635 * t776 - t636 * t777;
t513 = -qJ(4) * t571 - qJD(4) * t617 + t718;
t709 = -qJD(3) * t782 + t625 * t688 - t692 * t627;
t514 = qJ(4) * t570 - qJD(4) * t618 + t709;
t465 = t513 * t809 + t685 * t514;
t533 = -t570 * t685 + t571 * t809;
t668 = t689 * t810;
t756 = pkin(3) * t571 + t668;
t472 = pkin(4) * t533 - qJ(5) * t534 - qJD(5) * t568 + t756;
t436 = t686 * t465 + t684 * t472;
t752 = t809 * t546;
t496 = t685 * t537 + t752;
t492 = qJ(5) * t680 + t496;
t573 = -pkin(3) * t598 + qJD(4) - t634;
t506 = -pkin(4) * t564 - qJ(5) * t721 + t573;
t459 = t686 * t492 + t684 * t506;
t567 = t617 * t809 + t618 * t685;
t526 = pkin(4) * t567 - qJ(5) * t568 + t743;
t745 = -t692 * t635 - t636 * t688;
t560 = -qJ(4) * t618 + t745;
t561 = -qJ(4) * t617 + t782;
t528 = t685 * t560 + t561 * t809;
t476 = t684 * t526 + t686 * t528;
t556 = pkin(5) * t800;
t786 = -t556 + t785;
t784 = t589 - t516;
t662 = pkin(3) * t673;
t781 = t662 + t675;
t681 = t689 ^ 2;
t780 = -t693 ^ 2 + t681;
t774 = qJD(6) * t687;
t773 = qJD(6) * t691;
t767 = pkin(9) * t800;
t493 = t504 * t684 - t686 * t677;
t494 = t504 * t686 + t677 * t684;
t764 = -t687 * t493 + t691 * t494 - t548 * t773;
t763 = t662 + t737;
t628 = -pkin(2) * t689 - t820;
t755 = t628 - t819;
t423 = -t438 * t684 + t686 * t447;
t420 = pkin(5) * t503 - pkin(9) * t494 + t423;
t421 = -pkin(9) * t493 + t424;
t750 = t691 * t420 - t421 * t687;
t435 = -t465 * t684 + t686 * t472;
t748 = t691 * t493 + t687 * t494;
t458 = -t492 * t684 + t686 * t506;
t466 = -t502 * t684 + t686 * t520;
t464 = t513 * t685 - t809 * t514;
t475 = t686 * t526 - t528 * t684;
t501 = t545 * t685 + t752;
t527 = -t809 * t560 + t561 * t685;
t742 = pkin(5) * t721 - pkin(9) * t799;
t660 = -pkin(3) * t809 - pkin(4);
t741 = -t819 - t820;
t738 = -t613 * t500 + t559 * t831;
t590 = t665 * t809 - t652;
t736 = t420 * t687 + t421 * t691;
t735 = -t423 * t686 - t424 * t684;
t734 = -t423 * t684 + t422;
t448 = -pkin(5) * t564 - pkin(9) * t550 + t458;
t452 = -pkin(9) * t548 + t459;
t426 = t448 * t691 - t452 * t687;
t427 = t448 * t687 + t452 * t691;
t457 = pkin(5) * t567 - pkin(9) * t797 + t475;
t462 = -pkin(9) * t798 + t476;
t733 = t457 * t691 - t462 * t687;
t732 = t457 * t687 + t462 * t691;
t731 = t458 * t684 - t459 * t686;
t730 = t495 * t564 + t496 * t721;
t586 = -pkin(4) - t590;
t726 = -0.2e1 * pkin(1) * t772 - pkin(7) * qJDD(2);
t575 = t674 + t794;
t725 = qJD(6) * t575 + t742 - t826;
t574 = (-pkin(9) - t585) * t684;
t724 = -qJD(6) * t574 - t749 - t767;
t601 = t674 + t793;
t723 = qJD(5) * t684 + qJD(6) * t601 + t466 + t742;
t600 = (-pkin(9) - t656) * t684;
t722 = -qJD(6) * t600 - t767 + t827;
t440 = -t550 * t774 + t764;
t716 = -t503 * t585 + t564 * t582 - t806;
t715 = qJD(5) * t564 - t503 * t656 - t806;
t695 = qJD(2) ^ 2;
t714 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t695 + t825;
t696 = qJD(1) ^ 2;
t713 = pkin(1) * t696 - pkin(7) * qJDD(1) + t740;
t712 = t442 * t568 + t491 * t534 - t740;
t711 = -t458 * t721 + t686 * t839;
t441 = -qJD(6) * t729 + t748;
t707 = -g(3) * t658 + t458 * t799 + t459 * t800 - t659 * t740 + t422;
t706 = g(3) * t792 + t459 * t721 + t684 * t841;
t432 = t493 * pkin(5) + t442;
t479 = t548 * pkin(5) + t491;
t705 = -t426 * t721 + t432 * t613 - t479 * t831 - t659 * t814 + t670 * t838;
t704 = g(3) * t672 - t688 * t572 + t634 * t598 + t673 * t740 + t744 * t692 + t594;
t703 = t599 * t598 * MDP(11) - t559 * t721 * MDP(28) + (-t440 * t613 - t441 * t614 - t729 * t831 - t830 * t837) * MDP(25) + (t440 * t614 + t729 * t830) * MDP(24) + (t803 + t824) * MDP(26) + (t738 - t804) * MDP(27) + (-t598 * t680 + t547) * MDP(13) + (-t599 * t680 + t698) * MDP(14) + (-t598 ^ 2 + t599 ^ 2) * MDP(12) + t677 * MDP(15);
t702 = t427 * t721 + t432 * t614 - t479 * t830 + t659 * t815 - t669 * t838;
t699 = -t634 * t599 + t710 + t836;
t678 = -qJ(4) - t823;
t631 = t694 * t807;
t630 = t690 * t807;
t629 = t660 - t812;
t623 = pkin(1) + t781;
t610 = t694 * t623;
t593 = t655 - t835;
t581 = t586 - t812;
t580 = t659 * t788 + t791;
t579 = -t659 * t790 + t789;
t578 = -t659 * t789 + t790;
t577 = t659 * t791 + t788;
t530 = t613 * t568;
t529 = t614 * t568;
t485 = pkin(5) * t798 + t527;
t483 = t501 + t556;
t461 = t534 * t614 + t773 * t797 - t774 * t798;
t460 = -t534 * t613 - t568 * t597;
t451 = pkin(5) * t802 + t464;
t429 = -pkin(9) * t802 + t436;
t428 = pkin(5) * t533 - t534 * t674 + t435;
t1 = [(-t598 * t668 + t825 * t673 + t745 * t677 + t709 * t680 + (t593 - t835) * t617 - 0.2e1 * t634 * t571) * MDP(16) + qJDD(1) * MDP(1) + (-t440 * t529 + t441 * t530 + t460 * t837 + t461 * t729) * MDP(25) + (-t441 * t567 - t461 * t559 - t500 * t529 + t533 * t837) * MDP(27) + ((t428 * t691 - t429 * t687) * t559 + t733 * t500 + t750 * t567 + t426 * t533 - t451 * t837 + t485 * t441 + t432 * t529 + t479 * t461 - g(1) * t578 - g(2) * t580 + (-t427 * t567 - t559 * t732) * qJD(6)) * MDP(29) + (-t443 * t568 - t444 * t567 + t464 * t721 + t465 * t564 - t495 * t534 - t496 * t533 - t503 * t528 + t504 * t527 - t740) * MDP(18) + (-t424 * t567 + t436 * t564 - t459 * t533 + t464 * t550 - t476 * t503 + t527 * t494 + t686 * t712 - t792 * t825) * MDP(21) + (t659 * t686 * t825 + t423 * t567 - t435 * t564 + t458 * t533 + t464 * t548 + t475 * t503 + t527 * t493 + t684 * t712) * MDP(20) + (-t435 * t550 - t436 * t548 - t475 * t494 - t476 * t493 + t825 * t658 + t735 * t568 + (-t458 * t686 - t459 * t684) * t534) * MDP(22) + t825 * MDP(2) + (-t547 * t617 - t570 * t598 + t599 * t571 + t618 * t698) * MDP(12) + (-t570 * t680 + t618 * t677) * MDP(13) + (t547 * t618 + t570 * t599) * MDP(11) + (t440 * t567 + t460 * t559 - t500 * t530 - t533 * t729) * MDP(26) + (-(t428 * t687 + t429 * t691) * t559 - t732 * t500 - t736 * t567 - t427 * t533 - t451 * t729 + t485 * t440 - t432 * t530 + t479 * t460 - g(1) * t577 - g(2) * t579 + (-t426 * t567 - t559 * t733) * qJD(6)) * MDP(30) + (-t440 * t530 - t460 * t729) * MDP(24) + (-t547 * t811 + t634 * t570 + t593 * t618 - t599 * t668 - t672 * t825 - t677 * t782 - t680 * t718) * MDP(17) + (qJDD(2) * t689 + t693 * t695) * MDP(6) + (qJDD(2) * t693 - t689 * t695) * MDP(7) + (-t571 * t680 - t617 * t677) * MDP(14) + 0.2e1 * (t689 * t769 - t772 * t780) * MDP(5) + (qJDD(1) * t681 + 0.2e1 * t689 * t758) * MDP(4) + (t500 * t567 + t533 * t559) * MDP(28) + (t444 * t528 + t496 * t465 - t443 * t527 - t495 * t464 + t697 * t743 + t573 * t756 - g(1) * (-t623 * t690 - t678 * t694) - g(2) * (-t678 * t690 + t610)) * MDP(19) + (t689 * t726 + t693 * t714) * MDP(9) + (-t689 * t714 + t693 * t726) * MDP(10) + (-g(2) * t610 + t423 * t475 + t424 * t476 + t458 * t435 + t459 * t436 + t442 * t527 + t491 * t464 + (g(1) * t678 - g(2) * t737) * t694 + (-g(1) * (-t623 - t737) + g(2) * t678) * t690) * MDP(23) + t740 * MDP(3); (t444 * t591 + t443 * t590 - t573 * (t667 - t822) - g(3) * t781 - t740 * t628 + t784 * t496 - t785 * t495) * MDP(19) + (-t503 * t591 - t504 * t590 + t564 * t784 + t721 * t785 + t730) * MDP(18) + (t442 * t586 - g(1) * (t694 * t755 + t631) - g(2) * (t690 * t755 + t630) - g(3) * (t675 + t763) + t734 * t585 + t785 * t491 + t749 * t459 + t826 * t458) * MDP(23) + (-t469 * t564 + t494 * t586 + t550 * t785 + t686 * t716 + t706) * MDP(21) + (-t746 * t680 + (t598 * t779 + t677 * t692 - t680 * t777) * pkin(2) + t699) * MDP(16) + (-(t574 * t687 + t575 * t691) * t500 + t581 * t440 + (t687 * t725 + t691 * t724) * t559 - t786 * t729 + t702) * MDP(30) + (-g(3) * t693 + t689 * t713) * MDP(9) + ((t574 * t691 - t575 * t687) * t500 + t581 * t441 + (t687 * t724 - t691 * t725) * t559 - t786 * t837 + t705) * MDP(29) + qJDD(2) * MDP(8) + (t468 * t564 + t493 * t586 + t548 * t785 + t684 * t716 + t711) * MDP(20) + (g(3) * t689 + t693 * t713) * MDP(10) + MDP(7) * t769 + (-t493 * t794 + t468 * t550 - t749 * t548 + (t494 * t585 + t550 * t582 - t423) * t684 + t707) * MDP(22) + (t783 * t680 + (t599 * t779 - t677 * t688 - t680 * t776) * pkin(2) + t704) * MDP(17) + MDP(6) * t770 + t703 + (-MDP(4) * t689 * t693 + MDP(5) * t780) * t696; (-t493 * t793 + t466 * t550 + t827 * t548 + (qJD(5) * t550 + t494 * t656 - t423) * t684 + t707) * MDP(22) + ((t600 * t691 - t601 * t687) * t500 + t629 * t441 + t483 * t837 + (t687 * t722 - t691 * t723) * t559 + t705) * MDP(29) + (t495 * t501 - t496 * t502 + (t443 * t809 + t444 * t685 + t573 * t599 + t836) * pkin(3)) * MDP(19) + (-t467 * t564 + t494 * t660 - t501 * t550 + t686 * t715 + t706) * MDP(21) + (-t680 * t728 + t699) * MDP(16) + (-(t600 * t687 + t601 * t691) * t500 + t629 * t440 + t483 * t729 + (t687 * t723 + t691 * t722) * t559 + t702) * MDP(30) + (t466 * t564 + t493 * t660 - t501 * t548 + t684 * t715 + t711) * MDP(20) + (t680 * t747 + t704) * MDP(17) + (t442 * t660 - t459 * t467 - t458 * t466 - t491 * t501 - g(1) * (t694 * t741 + t631) - g(2) * (t690 * t741 + t630) - g(3) * t763 + t734 * t656 - t731 * qJD(5)) * MDP(23) + (-t501 * t721 - t502 * t564 + (-t503 * t685 - t504 * t809) * pkin(3) + t730) * MDP(18) + t703; -t721 ^ 2 * MDP(18) + (t503 * t686 - t548 * t721) * MDP(20) + (-t503 * t684 - t550 * t721) * MDP(21) + (-t493 * t684 - t494 * t686) * MDP(22) + (-t491 * t721 - t735 - t825) * MDP(23) + (t738 + t804) * MDP(29) + (t803 - t824) * MDP(30) + ((t548 * t686 - t550 * t684) * MDP(22) + t731 * MDP(23) - (MDP(20) * t684 + MDP(21) * t686 + MDP(18)) * t564) * t564 + (t495 * t721 - t496 * t564 + t697 - t825) * MDP(19); (-t550 * t564 + t493) * MDP(20) + (t548 * t564 + t494) * MDP(21) + (-t548 ^ 2 - t550 ^ 2) * MDP(22) + (t458 * t550 + t459 * t548 - t839) * MDP(23) + (t441 - t834) * MDP(29) + (t440 + t840) * MDP(30); t729 * t837 * MDP(24) + (t729 ^ 2 - t837 ^ 2) * MDP(25) + (t764 - t840) * MDP(26) + (-t748 - t834) * MDP(27) + t500 * MDP(28) + (-g(1) * t579 + g(2) * t577 + t427 * t559 + t479 * t729 + t658 * t815 + t750) * MDP(29) + (g(1) * t580 - g(2) * t578 + t426 * t559 - t479 * t837 + t658 * t814 - t736) * MDP(30) + (-MDP(26) * t801 + MDP(27) * t729 - MDP(29) * t427 - MDP(30) * t426) * qJD(6);];
tau  = t1;
