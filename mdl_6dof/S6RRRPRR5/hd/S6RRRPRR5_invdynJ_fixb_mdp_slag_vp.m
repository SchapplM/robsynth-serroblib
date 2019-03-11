% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:35
% EndTime: 2019-03-09 18:23:49
% DurationCPUTime: 11.68s
% Computational Cost: add. (7490->620), mult. (16148->764), div. (0->0), fcn. (11939->14), ass. (0->283)
t668 = sin(qJ(3));
t669 = sin(qJ(2));
t673 = cos(qJ(2));
t813 = cos(qJ(3));
t604 = t668 * t673 + t669 * t813;
t588 = t604 * qJD(1);
t666 = sin(qJ(6));
t667 = sin(qJ(5));
t671 = cos(qJ(6));
t672 = cos(qJ(5));
t601 = t666 * t672 + t667 * t671;
t819 = qJD(5) + qJD(6);
t774 = (t588 + t819) * t601;
t660 = qJD(2) + qJD(3);
t784 = t668 * t669;
t717 = t660 * t784;
t747 = t813 * t673;
t727 = qJD(1) * t747;
t741 = qJDD(1) * t813;
t753 = qJDD(1) * t673;
t728 = t660 * t727 + t668 * t753 + t669 * t741;
t515 = qJD(1) * t717 - t728;
t508 = -qJDD(5) + t515;
t507 = -qJDD(6) + t508;
t828 = qJD(5) + t588;
t573 = qJD(6) + t828;
t760 = qJD(6) * t666;
t763 = qJD(5) * t667;
t782 = t671 * t672;
t820 = -t666 * t667 + t782;
t829 = t820 * t588 - t666 * t763 - t667 * t760 + t782 * t819;
t835 = t601 * t507 - t573 * t829;
t766 = qJD(1) * t669;
t746 = t668 * t766;
t586 = -t727 + t746;
t555 = -t672 * t586 + t660 * t667;
t834 = t555 * t828;
t602 = -t747 + t784;
t533 = t820 * t602;
t833 = -t507 * t820 - t573 * t774;
t557 = t586 * t667 + t660 * t672;
t799 = t557 * t666;
t491 = t671 * t555 + t799;
t832 = t491 * t573;
t713 = t555 * t666 - t671 * t557;
t831 = t573 * t713;
t536 = pkin(3) * t588 + qJ(4) * t586;
t581 = t588 * pkin(9);
t509 = t536 + t581;
t815 = pkin(3) + pkin(9);
t761 = qJD(5) * t815;
t830 = t509 + t761;
t814 = pkin(8) + pkin(7);
t615 = t814 * t673;
t607 = qJD(1) * t615;
t590 = t668 * t607;
t614 = t814 * t669;
t605 = qJD(1) * t614;
t804 = qJD(2) * pkin(2);
t594 = -t605 + t804;
t538 = -t813 * t594 + t590;
t757 = qJD(4) + t538;
t807 = t588 * pkin(4);
t758 = t807 + t757;
t486 = -t660 * t815 + t758;
t812 = pkin(2) * t673;
t646 = pkin(1) + t812;
t613 = t646 * qJD(1);
t695 = -qJ(4) * t588 - t613;
t489 = t586 * t815 + t695;
t457 = t486 * t667 + t489 * t672;
t446 = -pkin(10) * t555 + t457;
t443 = t446 * t760;
t593 = t813 * t607;
t539 = t668 * t594 + t593;
t808 = t586 * pkin(4);
t514 = t539 - t808;
t653 = t660 * qJ(4);
t498 = t514 + t653;
t477 = pkin(5) * t555 + t498;
t664 = qJ(5) + qJ(6);
t654 = sin(t664);
t656 = cos(t664);
t670 = sin(qJ(1));
t792 = t656 * t670;
t665 = qJ(2) + qJ(3);
t655 = sin(t665);
t674 = cos(qJ(1));
t793 = t655 * t674;
t564 = t654 * t793 + t792;
t791 = t656 * t674;
t794 = t655 * t670;
t566 = -t654 * t794 + t791;
t657 = cos(t665);
t642 = g(3) * t657;
t827 = g(1) * t564 - g(2) * t566 + t477 * t491 - t642 * t654 + t443;
t563 = -t654 * t670 + t655 * t791;
t565 = t654 * t674 + t655 * t792;
t548 = t660 * t604;
t754 = qJDD(1) * t669;
t716 = t668 * t754 - t673 * t741;
t516 = qJD(1) * t548 + t716;
t659 = qJDD(2) + qJDD(3);
t762 = qJD(5) * t672;
t474 = t667 * t516 + t586 * t762 + t672 * t659 - t660 * t763;
t755 = qJD(1) * qJD(2);
t744 = t669 * t755;
t583 = pkin(2) * t744 - qJDD(1) * t646;
t683 = qJ(4) * t515 - qJD(4) * t588 + t583;
t435 = t516 * t815 + t683;
t743 = t673 * t755;
t553 = qJDD(2) * pkin(2) + t814 * (-t743 - t754);
t554 = t814 * (-t744 + t753);
t745 = qJD(3) * t813;
t765 = qJD(3) * t668;
t730 = -t813 * t553 + t668 * t554 + t594 * t765 + t607 * t745;
t711 = qJDD(4) + t730;
t442 = -pkin(4) * t515 - t659 * t815 + t711;
t739 = -t435 * t667 + t672 * t442;
t685 = -qJD(5) * t457 + t739;
t419 = -pkin(5) * t508 - pkin(10) * t474 + t685;
t736 = -t672 * t516 + t659 * t667;
t475 = qJD(5) * t557 + t736;
t751 = -t672 * t435 - t667 * t442 - t486 * t762;
t699 = -t489 * t763 - t751;
t420 = -pkin(10) * t475 + t699;
t740 = t671 * t419 - t666 * t420;
t826 = -g(1) * t563 - g(2) * t565 + t477 * t713 + t656 * t642 + t740;
t825 = (-t491 ^ 2 + t713 ^ 2) * MDP(30) - t507 * MDP(33) - t491 * MDP(29) * t713;
t534 = t601 * t602;
t710 = -qJ(4) * t604 - t646;
t519 = t602 * t815 + t710;
t560 = t614 * t813 + t668 * t615;
t530 = t604 * pkin(4) + t560;
t526 = t667 * t530;
t777 = t672 * t519 + t526;
t542 = -t668 * t605 + t593;
t520 = t542 - t808;
t752 = pkin(2) * t765;
t824 = (-t520 + t752) * t672;
t823 = -t586 * pkin(5) - pkin(10) * t763;
t543 = -t605 * t813 - t590;
t772 = -pkin(2) * t745 - qJD(4) + t543;
t769 = t657 * pkin(3) + t655 * qJ(4);
t647 = t659 * qJ(4);
t822 = -t660 * qJD(4) - t647;
t749 = -pkin(5) * t672 - pkin(4);
t821 = pkin(5) * t762 - t749 * t588;
t748 = qJD(2) * t814;
t606 = t669 * t748;
t608 = t673 * t748;
t495 = t813 * t606 + t668 * t608 + t614 * t745 + t615 * t765;
t561 = -t668 * t614 + t615 * t813;
t722 = g(1) * t670 - g(2) * t674;
t818 = -t495 * t660 + t561 * t659 + t655 * t722;
t496 = qJD(3) * t561 - t668 * t606 + t608 * t813;
t817 = t496 * t660 + t560 * t659 - t657 * t722;
t737 = t474 * t666 + t671 * t475;
t433 = -qJD(6) * t713 + t737;
t816 = t588 ^ 2;
t811 = pkin(3) * t659;
t810 = pkin(10) * t588;
t809 = pkin(10) * t602;
t641 = g(3) * t655;
t645 = -pkin(2) * t813 - pkin(3);
t636 = -pkin(9) + t645;
t806 = -pkin(10) + t636;
t805 = -pkin(10) - t815;
t456 = t672 * t486 - t489 * t667;
t445 = -pkin(10) * t557 + t456;
t436 = pkin(5) * t828 + t445;
t803 = t436 * t671;
t802 = t446 * t671;
t801 = t474 * t672;
t800 = t539 * t660;
t798 = t586 * t588;
t797 = t588 * t672;
t795 = t602 * t667;
t790 = t657 * t670;
t789 = t657 * t674;
t787 = t667 * t508;
t786 = t667 * t670;
t785 = t667 * t674;
t783 = t670 * t672;
t501 = t672 * t508;
t781 = t672 * t674;
t780 = t815 * t508;
t529 = pkin(2) * t766 + t536;
t497 = t529 + t581;
t779 = t672 * t497 + t667 * t520;
t778 = t672 * t509 + t667 * t514;
t773 = t821 - t772;
t771 = t807 - t772;
t770 = t757 + t821;
t662 = t669 ^ 2;
t768 = -t673 ^ 2 + t662;
t764 = qJD(5) * t636;
t759 = qJD(6) * t671;
t651 = t669 * t804;
t750 = t671 * t474 - t666 * t475 - t555 * t759;
t596 = t806 * t672;
t610 = t805 * t672;
t638 = pkin(2) * t668 + qJ(4);
t742 = -t519 - t809;
t729 = -t668 * t553 - t813 * t554 - t594 * t745 + t607 * t765;
t465 = t729 + t822;
t444 = -pkin(4) * t516 - t465;
t738 = t444 * t672 - t457 * t586;
t735 = t828 ^ 2;
t734 = t828 * t498;
t733 = t828 * t557;
t732 = qJD(6) * t436 + t420;
t726 = -t542 + t752;
t724 = -pkin(2) * t669 - pkin(3) * t655;
t723 = g(1) * t674 + g(2) * t670;
t595 = t806 * t667;
t721 = qJD(6) * t595 + (-t497 - t810) * t667 + t636 * t763 + t823 - t824;
t504 = t672 * t514;
t609 = t805 * t667;
t720 = qJD(6) * t609 + t504 + t823 + (-t810 - t830) * t667;
t568 = pkin(10) * t797;
t719 = -t596 * t819 - t667 * t752 + t568 + t779;
t718 = -t610 * t819 + t568 + t778;
t425 = t436 * t666 + t802;
t709 = t444 * t667 + t456 * t586 + (t762 + t797) * t498;
t708 = t646 + t769;
t707 = -0.2e1 * pkin(1) * t755 - pkin(7) * qJDD(2);
t706 = -t667 * t735 - t501;
t705 = t548 * t667 + t602 * t762;
t704 = -t548 * t672 + t602 * t763;
t547 = -qJD(2) * t747 - t673 * t745 + t717;
t703 = qJ(4) * t547 - qJD(4) * t604 + t651;
t461 = t548 * t815 + t703;
t473 = -t547 * pkin(4) + t496;
t698 = t672 * t461 + t667 * t473 - t519 * t763 + t530 * t762;
t432 = -t557 * t760 + t750;
t697 = -g(1) * t789 - g(2) * t790 - t641 - t729;
t696 = -g(1) * t793 - g(2) * t794 + t642 + t730;
t694 = -t672 * t735 + t787;
t692 = -t657 * t723 - t641;
t472 = -pkin(4) * t548 - t495;
t677 = qJD(2) ^ 2;
t690 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t677 + t722;
t678 = qJD(1) ^ 2;
t689 = pkin(1) * t678 - pkin(7) * qJDD(1) + t723;
t688 = -t613 * t586 - t697;
t687 = t613 * t588 - t696;
t523 = pkin(3) * t586 + t695;
t686 = t523 * t588 + qJDD(4) + t696;
t684 = -t523 * t586 + t697 - t822;
t682 = t728 + (t586 - t746) * t660;
t424 = -t446 * t666 + t803;
t431 = pkin(5) * t475 + t444;
t681 = t424 * t586 + t431 * t601 + t477 * t829 + t654 * t692;
t680 = -t425 * t586 + t431 * t820 - t477 * t774 + t656 * t692;
t679 = MDP(11) * t798 + (-t432 * t601 - t433 * t820 + t491 * t774 + t713 * t829) * MDP(30) + (t432 * t820 + t713 * t774) * MDP(29) + ((-t475 - t733) * t672 + (-t474 + t834) * t667) * MDP(23) + (-t586 * t713 + t833) * MDP(31) + (-t491 * t586 + t835) * MDP(32) + (-t667 * t733 + t801) * MDP(22) + (t557 * t586 + t706) * MDP(24) + (-t555 * t586 + t694) * MDP(25) + t682 * MDP(13) - t716 * MDP(14) + (-t586 ^ 2 + t816) * MDP(12) + t659 * MDP(15) + (MDP(26) * t828 + t573 * MDP(33)) * t586;
t658 = t667 * pkin(5);
t639 = qJ(4) + t658;
t617 = qJ(4) * t789;
t616 = qJ(4) * t790;
t612 = t638 + t658;
t579 = -t655 * t786 + t781;
t578 = t655 * t783 + t785;
t577 = t655 * t785 + t783;
t576 = t655 * t781 - t786;
t537 = pkin(3) * t602 + t710;
t532 = -t653 - t539;
t531 = -t602 * pkin(4) + t561;
t528 = -pkin(3) * t660 + t757;
t527 = t672 * t530;
t499 = t602 * t749 + t561;
t476 = pkin(3) * t548 + t703;
t468 = t672 * t473;
t466 = t711 - t811;
t464 = t672 * t809 + t777;
t459 = pkin(5) * t604 + t667 * t742 + t527;
t454 = pkin(3) * t516 + t683;
t453 = pkin(5) * t704 + t472;
t450 = t534 * t819 - t548 * t820;
t449 = t533 * t819 + t601 * t548;
t423 = -pkin(10) * t704 + t698;
t422 = -pkin(5) * t547 + t468 + (-pkin(10) * t548 - t461) * t667 + (t672 * t742 - t526) * qJD(5);
t1 = [(-t516 * t646 - t548 * t613 + t583 * t602 + t586 * t651 - t817) * MDP(16) + (-t454 * t602 - t476 * t586 - t516 * t537 - t523 * t548 + t817) * MDP(19) + (t515 * t646 + t547 * t613 + t583 * t604 + t588 * t651 - t818) * MDP(17) + (-t454 * t604 - t476 * t588 + t515 * t537 + t523 * t547 + t818) * MDP(20) + (t465 * t602 + t466 * t604 + t495 * t586 + t496 * t588 - t515 * t560 - t516 * t561 - t528 * t547 + t532 * t548 - t723) * MDP(18) + t722 * MDP(2) + t723 * MDP(3) + ((-t461 * t667 + t468) * t828 - (-t519 * t667 + t527) * t508 + t739 * t604 - t456 * t547 + t472 * t555 + t531 * t475 - g(1) * t579 - g(2) * t577 + (-t444 * t602 - t498 * t548) * t672 + (-t457 * t604 + t498 * t795 - t777 * t828) * qJD(5)) * MDP(27) + (-t508 * t604 - t547 * t828) * MDP(26) + (-t475 * t604 - t501 * t602 + t547 * t555 - t704 * t828) * MDP(25) + (t474 * t604 - t547 * t557 - t602 * t787 + t705 * t828) * MDP(24) + (g(1) * t578 - g(2) * t576 + t444 * t795 + t457 * t547 + t472 * t557 + t531 * t474 + t498 * t705 + t508 * t777 - t604 * t699 - t698 * t828) * MDP(28) + (t669 * t707 + t673 * t690) * MDP(9) + (-t669 * t690 + t673 * t707) * MDP(10) + qJDD(1) * MDP(1) + (t432 * t533 - t433 * t534 - t449 * t491 + t450 * t713) * MDP(30) + (-t433 * t604 - t450 * t573 + t491 * t547 - t507 * t533) * MDP(32) + ((t422 * t671 - t423 * t666) * t573 - (t459 * t671 - t464 * t666) * t507 + t740 * t604 - t424 * t547 + t453 * t491 + t499 * t433 - t431 * t533 + t477 * t450 - g(1) * t566 - g(2) * t564 + ((-t459 * t666 - t464 * t671) * t573 - t425 * t604) * qJD(6)) * MDP(34) + (t432 * t604 + t449 * t573 - t507 * t534 + t547 * t713) * MDP(31) + (g(1) * t565 - g(2) * t563 + t425 * t547 + t431 * t534 + t499 * t432 + t443 * t604 + t477 * t449 - t453 * t713 + (-(-qJD(6) * t464 + t422) * t573 + t459 * t507 - t419 * t604) * t666 + (-(qJD(6) * t459 + t423) * t573 + t464 * t507 - t732 * t604) * t671) * MDP(35) + (t432 * t534 - t449 * t713) * MDP(29) + (t454 * t537 - t465 * t561 + t466 * t560 + t523 * t476 + t532 * t495 + t528 * t496 + (-g(1) * t814 - g(2) * t708) * t674 + (g(1) * t708 - g(2) * t814) * t670) * MDP(21) + (qJDD(2) * t669 + t673 * t677) * MDP(6) + (qJDD(2) * t673 - t669 * t677) * MDP(7) + (-t548 * t660 - t602 * t659) * MDP(14) + (-t547 * t660 + t604 * t659) * MDP(13) + (t515 * t602 - t516 * t604 + t547 * t586 - t548 * t588) * MDP(12) + (-t515 * t604 - t547 * t588) * MDP(11) + (-t507 * t604 - t547 * t573) * MDP(33) + (qJDD(1) * t662 + 0.2e1 * t669 * t743) * MDP(4) + 0.2e1 * (t669 * t753 - t755 * t768) * MDP(5) + (t474 * t795 + t557 * t705) * MDP(22) + ((-t555 * t667 + t557 * t672) * t548 + (t801 - t475 * t667 + (-t555 * t672 - t557 * t667) * qJD(5)) * t602) * MDP(23); (g(3) * t669 + t673 * t689) * MDP(10) + (-t636 * t501 + t638 * t475 + t824 * t828 + t771 * t555 + ((t497 - t764) * t828 + t692) * t667 + t709) * MDP(27) + (-g(3) * t673 + t669 * t689) * MDP(9) + MDP(7) * t753 + (-t465 * t638 + t466 * t645 - t523 * t529 - g(1) * (t674 * t724 + t617) - g(2) * (t670 * t724 + t616) - g(3) * (t769 + t812) + t772 * t532 + t726 * t528) * MDP(21) + t679 + (-t515 * t645 - t516 * t638 + (-t532 + t726) * t588 + (t528 + t772) * t586) * MDP(18) + (t542 * t660 + (-t586 * t766 + t659 * t813 - t660 * t765) * pkin(2) + t687) * MDP(16) + (t638 * t474 + t779 * t828 + t771 * t557 + (-t764 * t828 + t692) * t672 + (t636 * t508 - t752 * t828 - t734) * t667 + t738) * MDP(28) + (t529 * t588 + t638 * t659 - t660 * t772 + t684) * MDP(20) + (t529 * t586 + t726 * t660 + (-pkin(3) + t645) * t659 + t686) * MDP(19) + (t543 * t660 + (-t588 * t766 - t659 * t668 - t660 * t745) * pkin(2) + t688) * MDP(17) + ((t595 * t671 + t596 * t666) * t507 + t612 * t432 + (t666 * t721 + t671 * t719) * t573 - t773 * t713 + t680) * MDP(35) + qJDD(2) * MDP(8) + (-(-t595 * t666 + t596 * t671) * t507 + t612 * t433 + (t666 * t719 - t671 * t721) * t573 + t773 * t491 + t681) * MDP(34) + MDP(6) * t754 + (-MDP(4) * t669 * t673 + MDP(5) * t768) * t678; ((t609 * t671 + t610 * t666) * t507 + t639 * t432 + (t666 * t720 + t671 * t718) * t573 - t770 * t713 + t680) * MDP(35) + (-t465 * qJ(4) - t466 * pkin(3) - t523 * t536 - t528 * t539 - g(1) * (-pkin(3) * t793 + t617) - g(2) * (-pkin(3) * t794 + t616) - g(3) * t769 - t757 * t532) * MDP(21) + t679 + (t672 * t780 + qJ(4) * t475 - t504 * t828 + t758 * t555 + (t828 * t830 + t692) * t667 + t709) * MDP(27) + (pkin(3) * t515 - qJ(4) * t516 + (-t532 - t539) * t588 + (t528 - t757) * t586) * MDP(18) + (qJ(4) * t474 + t778 * t828 + t758 * t557 + (-t734 - t780) * t667 + (t761 * t828 + t692) * t672 + t738) * MDP(28) + (t536 * t588 + t660 * t757 + t647 + t684) * MDP(20) + (t536 * t586 + t686 - t800 - 0.2e1 * t811) * MDP(19) + (t687 + t800) * MDP(16) + (-(-t609 * t666 + t610 * t671) * t507 + t639 * t433 + (t666 * t718 - t671 * t720) * t573 + t770 * t491 + t681) * MDP(34) + (-t538 * t660 + t688) * MDP(17); t682 * MDP(18) + (t659 - t798) * MDP(19) + (-t660 ^ 2 - t816) * MDP(20) + (t532 * t660 + t686 - t811) * MDP(21) + (-t555 * t660 + t706) * MDP(27) + (-t557 * t660 + t694) * MDP(28) + (-t660 * t491 + t833) * MDP(34) + (t660 * t713 + t835) * MDP(35); t557 * t555 * MDP(22) + (-t555 ^ 2 + t557 ^ 2) * MDP(23) + (t474 + t834) * MDP(24) + (-t736 + (-qJD(5) + t828) * t557) * MDP(25) - t508 * MDP(26) + (-g(1) * t576 - g(2) * t578 + t457 * t828 - t498 * t557 + t642 * t672 + t685) * MDP(27) + (g(1) * t577 - g(2) * t579 + t456 * t828 + t498 * t555 + (qJD(5) * t489 - t642) * t667 + t751) * MDP(28) + (t432 + t832) * MDP(31) + (-t433 - t831) * MDP(32) + (-(-t445 * t666 - t802) * t573 - t425 * qJD(6) + (-t491 * t557 - t507 * t671 - t573 * t760) * pkin(5) + t826) * MDP(34) + ((-t446 * t573 - t419) * t666 + (t445 * t573 - t732) * t671 + (t507 * t666 + t557 * t713 - t573 * t759) * pkin(5) + t827) * MDP(35) + t825; (t750 + t832) * MDP(31) + (-t737 - t831) * MDP(32) + (t425 * t573 + t826) * MDP(34) + (-t666 * t419 - t671 * t420 + t424 * t573 + t827) * MDP(35) + (-MDP(31) * t799 + MDP(32) * t713 - MDP(34) * t425 - MDP(35) * t803) * qJD(6) + t825;];
tau  = t1;
