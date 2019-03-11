% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:20:01
% EndTime: 2019-03-09 17:20:18
% DurationCPUTime: 13.01s
% Computational Cost: add. (5976->673), mult. (12755->823), div. (0->0), fcn. (8578->8), ass. (0->265)
t683 = sin(qJ(1));
t685 = cos(qJ(3));
t686 = cos(qJ(2));
t681 = sin(qJ(3));
t687 = cos(qJ(1));
t808 = t687 * t681;
t566 = -t683 * t685 + t686 * t808;
t809 = t686 * t687;
t567 = t681 * t683 + t685 * t809;
t680 = sin(qJ(5));
t684 = cos(qJ(5));
t511 = t566 * t684 - t567 * t680;
t811 = t683 * t686;
t564 = t681 * t811 + t685 * t687;
t810 = t685 * t686;
t565 = t683 * t810 - t808;
t858 = t564 * t684 - t565 * t680;
t883 = -g(1) * t511 - g(2) * t858;
t771 = t685 * qJD(2);
t682 = sin(qJ(2));
t784 = qJD(1) * t682;
t588 = t681 * t784 - t771;
t760 = t685 * t784;
t781 = qJD(2) * t681;
t590 = t760 + t781;
t613 = -qJD(2) * pkin(2) + pkin(7) * t784;
t509 = t588 * pkin(3) - t590 * qJ(4) + t613;
t485 = -pkin(4) * t588 - t509;
t522 = t588 * t680 + t590 * t684;
t814 = t682 * t685;
t817 = t681 * t684;
t553 = t680 * t814 - t682 * t817;
t783 = qJD(1) * t686;
t634 = -qJD(3) + t783;
t845 = pkin(3) + pkin(4);
t787 = -t686 * pkin(2) - t682 * pkin(8);
t866 = -pkin(1) + t787;
t575 = t866 * qJD(1);
t659 = pkin(7) * t783;
t614 = qJD(2) * pkin(8) + t659;
t527 = t685 * t575 - t681 * t614;
t770 = qJD(4) - t527;
t876 = pkin(9) * t590 - t770;
t477 = t634 * t845 - t876;
t528 = t681 * t575 + t685 * t614;
t493 = pkin(9) * t588 + t528;
t619 = t634 * qJ(4);
t483 = t493 - t619;
t457 = t477 * t680 + t483 * t684;
t768 = qJD(1) * qJD(2);
t749 = t686 * t768;
t767 = qJDD(1) * t682;
t777 = qJD(3) * t682;
t870 = qJD(1) * t777 - qJDD(2);
t513 = -qJD(3) * t771 + (-t749 - t767) * t685 + t870 * t681;
t662 = t686 * qJDD(1);
t860 = -t682 * t768 + t662;
t586 = qJDD(3) - t860;
t731 = pkin(2) * t682 - pkin(8) * t686;
t603 = t731 * qJD(2);
t535 = qJD(1) * t603 + qJDD(1) * t866;
t561 = pkin(7) * t860 + qJDD(2) * pkin(8);
t776 = qJD(3) * t685;
t778 = qJD(3) * t681;
t735 = -t685 * t535 + t681 * t561 + t575 * t778 + t614 * t776;
t719 = qJDD(4) + t735;
t450 = pkin(9) * t513 - t586 * t845 + t719;
t573 = t586 * qJ(4);
t617 = t634 * qJD(4);
t706 = t681 * t535 + t685 * t561 + t575 * t776 - t614 * t778;
t462 = t573 - t617 + t706;
t514 = t681 * (qJD(2) * (qJD(3) + t783) + t767) + t870 * t685;
t452 = pkin(9) * t514 + t462;
t745 = t684 * t450 - t680 * t452;
t693 = -t457 * qJD(5) + t745;
t882 = g(3) * t553 - t485 * t522 + t693 + t883;
t881 = -t685 * t783 + t776;
t761 = t681 * t783;
t880 = t761 - t778;
t775 = qJD(5) * t522;
t460 = -t513 * t680 - t684 * t514 + t775;
t721 = -t684 * t588 + t590 * t680;
t517 = t721 ^ 2;
t574 = -qJDD(5) + t586;
t846 = t522 ^ 2;
t856 = qJD(5) + t634;
t773 = qJD(5) * t684;
t774 = qJD(5) * t680;
t459 = t684 * t513 - t680 * t514 - t588 * t773 + t590 * t774;
t875 = -t721 * t856 + t459;
t879 = -(t517 - t846) * MDP(23) + MDP(22) * t522 * t721 - MDP(24) * t875 + (t522 * t856 - t460) * MDP(25) - t574 * MDP(26);
t844 = pkin(8) - pkin(9);
t616 = t844 * t685;
t763 = -pkin(7) * t681 - pkin(3);
t696 = -pkin(9) * t810 + (-pkin(4) + t763) * t682;
t600 = t731 * qJD(1);
t820 = t600 * t685;
t878 = qJD(1) * t696 - qJD(3) * t616 - t820;
t569 = t681 * t600;
t795 = qJ(4) * t784 + t569;
t816 = t681 * t686;
t872 = (-pkin(7) * t814 + pkin(9) * t816) * qJD(1) + t795 + t844 * t778;
t873 = qJ(6) * t522;
t592 = t680 * t681 + t684 * t685;
t530 = t592 * qJD(5) - t680 * t778 - t684 * t776;
t708 = t592 * t686;
t801 = qJD(1) * t708 + t530;
t800 = t680 * t881 + t681 * t773 + t684 * t880 - t685 * t774;
t813 = t682 * t687;
t815 = t682 * t683;
t871 = g(1) * t813 + g(2) * t815;
t865 = qJ(6) * t721;
t728 = g(1) * t687 + g(2) * t683;
t853 = t682 * t728;
t636 = pkin(7) * t816;
t670 = t686 * pkin(3);
t515 = pkin(4) * t686 + t636 + t670 + (-pkin(9) * t682 - t866) * t685;
t638 = pkin(7) * t810;
t790 = t681 * t866 + t638;
t537 = -qJ(4) * t686 + t790;
t818 = t681 * t682;
t526 = pkin(9) * t818 + t537;
t802 = t680 * t515 + t684 * t526;
t746 = -qJ(4) * t680 - t684 * t845;
t864 = -qJD(5) * t746 + t680 * t493 + t684 * t876;
t604 = qJ(4) * t684 - t680 * t845;
t863 = -qJD(5) * t604 - t684 * t493 + t680 * t876;
t862 = t878 * t684;
t615 = t844 * t681;
t793 = t680 * t615 + t684 * t616;
t861 = -t615 * t773 + t616 * t774 + t680 * t878 + t684 * t872;
t859 = qJ(4) * t881 + t681 * qJD(4) + t659;
t665 = t681 * qJ(4);
t855 = t685 * pkin(3) + pkin(2) + t665;
t832 = g(3) * t686;
t854 = -t832 + t871;
t840 = pkin(8) * t586;
t851 = t509 * t634 + t840;
t766 = t845 * t681;
t733 = -pkin(7) - t766;
t779 = qJD(2) * t686;
t755 = t686 * t771;
t791 = qJ(4) * t755 + qJD(4) * t814;
t478 = (-t685 * t845 - t665) * t777 + t733 * t779 + t791;
t824 = t566 * t680;
t512 = t567 * t684 + t824;
t554 = t592 * t682;
t705 = -t680 * t450 - t684 * t452 - t477 * t773 + t483 * t774;
t826 = t564 * t680;
t722 = t565 * t684 + t826;
t848 = -g(1) * t512 - g(2) * t722 - g(3) * t554 - t485 * t721 - t705;
t847 = -0.2e1 * pkin(1);
t843 = pkin(3) * t586;
t842 = pkin(3) * t634;
t841 = pkin(5) * t680;
t837 = g(1) * t683;
t834 = g(2) * t687;
t674 = g(3) * t682;
t655 = pkin(5) * t684 + pkin(4);
t831 = -pkin(3) - t655;
t830 = pkin(8) * qJD(3);
t499 = -t619 + t528;
t828 = t499 * t634;
t827 = t513 * t681;
t823 = t588 * t634;
t822 = t590 * t634;
t821 = t590 * t685;
t819 = t680 * t685;
t690 = qJD(1) ^ 2;
t812 = t682 * t690;
t456 = t684 * t477 - t483 * t680;
t446 = t456 - t873;
t445 = pkin(5) * t856 + t446;
t807 = -t446 + t445;
t806 = -qJ(6) * t800 - qJD(6) * t592 - t861;
t720 = -t817 + t819;
t805 = pkin(5) * t784 + qJ(6) * t801 - qJD(5) * t793 + qJD(6) * t720 + t680 * t872 - t862;
t799 = -qJD(3) * t766 + t761 * t845 + t859;
t798 = -t864 - t873;
t797 = t863 + t865;
t796 = -pkin(3) * t880 - t859;
t794 = t681 * t603 + t776 * t866;
t534 = t590 * pkin(3) + t588 * qJ(4);
t792 = t871 * t685;
t675 = t682 ^ 2;
t786 = -t686 ^ 2 + t675;
t782 = qJD(2) * t590;
t780 = qJD(2) * t682;
t765 = qJ(4) * t780 + t794;
t656 = pkin(7) * t767;
t562 = -qJDD(2) * pkin(2) + pkin(7) * t749 + t656;
t764 = g(1) * t809 + g(2) * t811 + t674;
t762 = pkin(3) * t681 + pkin(7);
t758 = t634 * t781;
t757 = t634 * t771;
t756 = t681 * t779;
t754 = t634 * t778;
t753 = t681 * t777;
t747 = qJ(4) + t841;
t726 = -qJD(3) * t638 + t603 * t685 - t778 * t866;
t473 = pkin(9) * t753 + qJD(2) * t696 - t726;
t474 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t814 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t681) * t686 + t765;
t744 = t684 * t473 - t474 * t680;
t741 = t684 * t515 - t526 * t680;
t739 = t684 * t615 - t616 * t680;
t738 = t685 * t866 - t636;
t737 = t856 ^ 2;
t581 = t685 * pkin(4) + t855;
t734 = pkin(3) * t810 + qJ(4) * t816 - t787;
t494 = -pkin(4) * t590 - t534;
t732 = t763 * t682;
t730 = -g(1) * t564 + g(2) * t566;
t729 = g(1) * t565 - g(2) * t567;
t727 = -t565 * pkin(3) + t687 * pkin(7) - qJ(4) * t564;
t725 = qJD(3) * t613 - t840;
t497 = t842 + t770;
t724 = t497 * t685 - t499 * t681;
t716 = -g(1) * t566 - g(2) * t564 - g(3) * t818;
t715 = -t634 * t830 + t832;
t463 = t514 * pkin(3) + t513 * qJ(4) - t590 * qJD(4) + t562;
t713 = -pkin(7) * qJDD(2) + t768 * t847;
t712 = t586 * t681 - t634 * t776;
t711 = t586 * t685 + t754;
t709 = g(3) * t720;
t707 = -t463 - t715;
t704 = t680 * t473 + t684 * t474 + t515 * t773 - t526 * t774;
t689 = qJD(2) ^ 2;
t700 = pkin(7) * t689 + qJDD(1) * t847 - t837;
t699 = t687 * pkin(1) + pkin(2) * t809 + t567 * pkin(3) + t683 * pkin(7) + pkin(8) * t813 + qJ(4) * t566;
t453 = -pkin(4) * t514 - t463;
t695 = -t716 - t735;
t444 = pkin(5) * t460 + qJDD(6) + t453;
t692 = -t528 * t634 + t695;
t691 = g(1) * t567 + g(2) * t565 + g(3) * t814 - t527 * t634 - t706;
t679 = -qJ(6) - pkin(9);
t646 = g(2) * t813;
t641 = pkin(8) * t809;
t637 = pkin(8) * t811;
t631 = qJ(4) * t814;
t599 = -pkin(5) + t746;
t559 = t566 * pkin(3);
t557 = t564 * pkin(3);
t547 = t682 * t762 - t631;
t538 = t670 - t738;
t536 = t682 * t733 + t631;
t533 = qJD(1) * t732 - t820;
t532 = -pkin(7) * t760 + t795;
t501 = -qJ(6) * t592 + t793;
t500 = qJ(6) * t720 + t739;
t489 = qJ(4) * t753 + pkin(7) * t779 + (t682 * t776 + t756) * pkin(3) - t791;
t484 = qJD(2) * t732 - t726;
t482 = -t513 - t823;
t481 = -qJD(4) * t686 + (-t682 * t771 - t686 * t778) * pkin(7) + t765;
t480 = qJD(2) * t708 + (qJD(3) - qJD(5)) * t682 * t720;
t479 = t530 * t682 + t680 * t755 - t684 * t756;
t469 = pkin(5) * t721 + qJD(6) + t485;
t468 = -qJ(6) * t553 + t802;
t467 = pkin(5) * t686 - qJ(6) * t554 + t741;
t466 = t719 - t843;
t447 = t457 - t865;
t443 = -qJ(6) * t479 - qJD(6) * t553 + t704;
t442 = -pkin(5) * t780 - qJ(6) * t480 - qJD(5) * t802 - qJD(6) * t554 + t744;
t441 = -qJ(6) * t460 - qJD(6) * t721 - t705;
t440 = -pkin(5) * t574 + qJ(6) * t459 - qJD(6) * t522 + t693;
t1 = [t728 * MDP(3) + (qJDD(1) * t675 + 0.2e1 * t682 * t749) * MDP(4) + (t744 * t856 - t741 * t574 + t745 * t686 - t456 * t780 + t478 * t721 + t536 * t460 + t453 * t553 + t485 * t479 + g(1) * t722 - g(2) * t512 + (-t457 * t686 - t802 * t856) * qJD(5)) * MDP(27) + (g(1) * t858 - g(2) * t511 + t453 * t554 + t457 * t780 - t536 * t459 + t478 * t522 + t485 * t480 + t574 * t802 + t686 * t705 - t704 * t856) * MDP(28) + (-t460 * t686 - t479 * t856 + t553 * t574 + t721 * t780) * MDP(25) + (-t574 * t686 - t780 * t856) * MDP(26) + (-t459 * t686 + t480 * t856 - t522 * t780 - t554 * t574) * MDP(24) + 0.2e1 * (t662 * t682 - t768 * t786) * MDP(5) + (t682 * t700 + t686 * t713 + t646) * MDP(10) + qJDD(1) * MDP(1) + (-t513 * t814 + (-t753 + t755) * t590) * MDP(11) + (-t726 * t634 + t738 * t586 + ((pkin(7) * t588 + t613 * t681) * qJD(2) + t735) * t686 + (t613 * t776 + t527 * qJD(2) + t562 * t681 + (t514 - t758) * pkin(7)) * t682 + t729) * MDP(16) + (-t481 * t634 - t489 * t590 + t513 * t547 + t537 * t586 + (-t509 * t771 - t462) * t686 + (qJD(2) * t499 - t463 * t685 + t509 * t778) * t682 - t730) * MDP(20) + ((t514 + t758) * t686 + (-qJD(2) * t588 - t712) * t682) * MDP(14) + (-g(1) * t727 - g(2) * t699 + t462 * t537 + t463 * t547 + t466 * t538 + t499 * t481 + t497 * t484 + t509 * t489 - t837 * t866) * MDP(21) + (-g(1) * t815 - t440 * t554 - t441 * t553 - t442 * t522 - t443 * t721 - t445 * t480 - t447 * t479 + t459 * t467 - t460 * t468 + t646) * MDP(29) + (t459 * t553 - t460 * t554 - t479 * t522 - t480 * t721) * MDP(23) + (-t586 * t686 - t634 * t780) * MDP(15) + (t484 * t634 + t489 * t588 + t514 * t547 - t538 * t586 + (t509 * t781 + t466) * t686 + (-qJD(2) * t497 + t463 * t681 + t509 * t776) * t682 + t729) * MDP(18) + ((t513 - t757) * t686 + (t711 + t782) * t682) * MDP(13) + (-t481 * t588 + t484 * t590 - t513 * t538 - t514 * t537 - t646 + t724 * t779 + (t837 - t462 * t681 + t466 * t685 + (-t497 * t681 - t499 * t685) * qJD(3)) * t682) * MDP(19) + (-t834 + t837) * MDP(2) + ((-t588 * t685 - t590 * t681) * t779 + (t827 - t514 * t685 + (t588 * t681 - t821) * qJD(3)) * t682) * MDP(12) + (t713 * t682 + (-t700 - t834) * t686) * MDP(9) + (t441 * t468 + t447 * t443 + t440 * t467 + t445 * t442 + t444 * (pkin(5) * t553 + t631) - g(1) * (-pkin(1) * t683 - pkin(2) * t811 - pkin(5) * t826 - t565 * t655 + t727) - g(2) * (pkin(5) * t824 + t567 * t655 + t699) + (t444 * (-pkin(4) * t681 - t762) - t679 * t834 - (-pkin(8) - t679) * t837) * t682 + (pkin(5) * t479 + t478) * t469) * MDP(30) + (t794 * t634 - t790 * t586 + (t613 * t771 + (-t754 + t782) * pkin(7) + t706) * t686 + (-t613 * t778 - t528 * qJD(2) + t562 * t685 + (-t513 - t757) * pkin(7)) * t682 + t730) * MDP(17) + (qJDD(2) * t682 + t686 * t689) * MDP(6) + (qJDD(2) * t686 - t682 * t689) * MDP(7) + (-t459 * t554 + t480 * t522) * MDP(22); (-t739 * t574 + t581 * t460 - g(3) * t708 + (-t616 * t773 + (-qJD(5) * t615 + t872) * t680 - t862) * t856 + t799 * t721 + t800 * t485 + (t453 + t853) * t592) * MDP(27) + (-t514 * t855 - t533 * t634 + t796 * t588 - t681 * t851 + t707 * t685 + t792) * MDP(18) + (-t513 * t855 + t532 * t634 - t796 * t590 + t851 * t685 + (t707 + t853) * t681) * MDP(20) + (t574 * t592 - t800 * t856) * MDP(25) + (t634 * MDP(15) + t497 * MDP(18) - t499 * MDP(20) + t522 * MDP(24) - MDP(25) * t721 + MDP(26) * t856 + t456 * MDP(27)) * t784 + (t574 * t720 - t801 * t856) * MDP(24) + (t793 * t574 - t581 * t459 - t453 * t720 + t686 * t709 + t861 * t856 + t799 * t522 - t801 * t485 + (-t457 * qJD(1) - t720 * t728) * t682) * MDP(28) + ((pkin(1) * t690 - pkin(7) * qJDD(1)) * t686 + t764) * MDP(10) + qJDD(2) * MDP(8) + ((t588 * t682 - t634 * t816) * qJD(1) + t711) * MDP(14) + ((-t590 * t682 + t634 * t810) * qJD(1) + t712) * MDP(13) + (pkin(2) * t513 - t569 * t634 + t725 * t685 + (-t613 * t810 + t528 * t682 + (-t590 * t686 + t634 * t814) * pkin(7)) * qJD(1) + (t562 + t715 - t853) * t681) * MDP(17) + (t440 * t720 - t441 * t592 + t445 * t801 - t447 * t800 + t459 * t500 - t460 * t501 - t522 * t805 - t721 * t806 + t764) * MDP(29) + (t459 * t592 + t460 * t720 - t522 * t800 + t721 * t801) * MDP(23) + (t459 * t720 - t522 * t801) * MDP(22) + (t441 * t501 + t440 * t500 + t444 * (pkin(5) * t592 + t581) - g(1) * (t679 * t809 + t641) - g(2) * (t679 * t811 + t637) - g(3) * (t655 * t810 + t816 * t841 + t734) + (t800 * pkin(5) + (pkin(4) * t634 + t842) * t681 + t859) * t469 + t806 * t447 + t805 * t445 + (-g(3) * t679 + t728 * (t681 * t747 - t685 * t831 + pkin(2))) * t682) * MDP(30) + (pkin(1) * t812 - t656 + t854) * MDP(9) - t686 * MDP(4) * t812 + t786 * t690 * MDP(5) + (-t634 * t821 - t827) * MDP(11) + (t532 * t588 - t533 * t590 + (t462 - t634 * t497 + (qJD(3) * t590 - t514) * pkin(8)) * t685 + (t466 + t828 + (qJD(3) * t588 - t513) * pkin(8)) * t681 - t764) * MDP(19) + (-pkin(2) * t514 + t725 * t681 + (-t832 - t562 + (t600 + t830) * t634) * t685 + (-t613 * t816 - t527 * t682 + (-t588 * t686 + t634 * t818) * pkin(7)) * qJD(1) + t792) * MDP(16) + ((-t513 + t823) * t685 + (-t514 + t822) * t681) * MDP(12) + (-t499 * t532 - t497 * t533 - g(1) * t641 - g(2) * t637 - g(3) * t734 + t796 * t509 + (qJD(3) * t724 + t462 * t685 + t466 * t681) * pkin(8) + (-t463 + t853) * t855) * MDP(21) + MDP(7) * t662 + MDP(6) * t767; t590 * t588 * MDP(11) + (-t588 ^ 2 + t590 ^ 2) * MDP(12) + t482 * MDP(13) + (-t514 - t822) * MDP(14) + t586 * MDP(15) + (-t590 * t613 + t692) * MDP(16) + (t588 * t613 + t691) * MDP(17) + (-t509 * t590 - t534 * t588 - qJDD(4) + t692 + 0.2e1 * t843) * MDP(18) + (pkin(3) * t513 - qJ(4) * t514 + (t499 - t528) * t590 + (t497 - t770) * t588) * MDP(19) + (-t509 * t588 + t534 * t590 + 0.2e1 * t573 - 0.2e1 * t617 - t691) * MDP(20) + (t462 * qJ(4) - t466 * pkin(3) - t509 * t534 - t497 * t528 - g(1) * (qJ(4) * t567 - t559) - g(2) * (qJ(4) * t565 - t557) - g(3) * (-pkin(3) * t818 + t631) + t770 * t499) * MDP(21) + (-t494 * t721 - t746 * t574 + t856 * t863 - t882) * MDP(27) + (-t494 * t522 + t604 * t574 + t856 * t864 + t848) * MDP(28) + (t459 * t599 - t460 * t604 + (t445 - t798) * t721 + (-t447 - t797) * t522) * MDP(29) + (t441 * t604 + t440 * t599 - t469 * (-pkin(5) * t522 + t494) - g(1) * (-t566 * t655 + t567 * t747 - t559) - g(2) * (-t564 * t655 + t565 * t747 - t557) - g(3) * t631 - (pkin(5) * t819 + t681 * t831) * t674 + t798 * t447 + t797 * t445) * MDP(30) - t879; -t586 * MDP(18) + t482 * MDP(19) - t634 ^ 2 * MDP(20) + (qJDD(4) - t695 + t828 - t843) * MDP(21) + t716 * MDP(30) + (MDP(18) * t588 - MDP(20) * t590 + MDP(21) * t509 - MDP(27) * t721 - MDP(28) * t522 - t469 * MDP(30)) * t590 + (-t574 * MDP(27) + t875 * MDP(29) + (t447 * t856 + t440) * MDP(30) - MDP(28) * t737) * t684 + (t574 * MDP(28) + (t522 * t634 - t460 + t775) * MDP(29) + (-t445 * t856 + t441) * MDP(30) - MDP(27) * t737) * t680; (t457 * t856 + t882) * MDP(27) + (t456 * t856 - t848) * MDP(28) + (pkin(5) * t459 - t721 * t807) * MDP(29) + (t807 * t447 + (-t469 * t522 + t682 * t709 + t440 + t883) * pkin(5)) * MDP(30) + t879; (-t517 - t846) * MDP(29) + (t445 * t522 + t447 * t721 + t444 + t854) * MDP(30);];
tau  = t1;
