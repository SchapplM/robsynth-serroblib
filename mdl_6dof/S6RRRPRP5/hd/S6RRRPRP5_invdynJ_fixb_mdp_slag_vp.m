% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:58
% EndTime: 2019-03-09 16:52:21
% DurationCPUTime: 17.40s
% Computational Cost: add. (13002->702), mult. (29028->863), div. (0->0), fcn. (20974->14), ass. (0->311)
t735 = sin(qJ(3));
t736 = sin(qJ(2));
t851 = qJD(1) * t736;
t822 = t735 * t851;
t738 = cos(qJ(3));
t837 = t738 * qJD(2);
t675 = t822 - t837;
t848 = qJD(2) * t735;
t677 = t738 * t851 + t848;
t732 = sin(pkin(10));
t893 = cos(pkin(10));
t605 = -t732 * t675 + t677 * t893;
t734 = sin(qJ(5));
t739 = cos(qJ(2));
t850 = qJD(1) * t739;
t810 = qJD(3) + t850;
t833 = qJDD(1) * t736;
t765 = qJD(2) * t810 + t833;
t755 = t765 * t735;
t834 = qJD(1) * qJD(3);
t800 = t736 * t834 - qJDD(2);
t597 = t738 * t800 + t755;
t782 = t800 * t735;
t750 = t765 * t738 - t782;
t746 = -t732 * t597 + t750 * t893;
t747 = -t597 * t893 - t732 * t750;
t784 = t675 * t893 + t732 * t677;
t906 = cos(qJ(5));
t762 = t906 * t784;
t840 = qJD(5) * t734;
t492 = qJD(5) * t762 + t605 * t840 - t734 * t747 - t906 * t746;
t556 = t605 * t734 + t762;
t705 = -qJD(3) + t850;
t696 = -qJD(5) + t705;
t885 = t556 * t696;
t483 = -t492 - t885;
t924 = t605 * t906 - t734 * t784;
t493 = qJD(5) * t924 + t734 * t746 - t906 * t747;
t719 = t739 * qJDD(1);
t835 = qJD(1) * qJD(2);
t913 = -t736 * t835 + t719;
t666 = qJDD(3) - t913;
t663 = qJDD(5) + t666;
t886 = t556 ^ 2;
t887 = t924 * t696;
t936 = t924 ^ 2;
t946 = t556 * t924;
t947 = t483 * MDP(22) + MDP(20) * t946 + t663 * MDP(24) + (-t493 - t887) * MDP(23) + (-t886 + t936) * MDP(21);
t874 = t738 * t739;
t904 = pkin(3) * t736;
t795 = -qJ(4) * t874 + t904;
t733 = -qJ(4) - pkin(8);
t815 = qJD(3) * t733;
t806 = pkin(2) * t736 - pkin(8) * t739;
t679 = t806 * qJD(1);
t859 = pkin(7) * t822 + t738 * t679;
t945 = -qJD(1) * t795 - qJD(4) * t735 + t738 * t815 - t859;
t659 = t735 * t679;
t841 = qJD(4) * t738;
t877 = t736 * t738;
t880 = t735 * t739;
t944 = t659 + (-pkin(7) * t877 - qJ(4) * t880) * qJD(1) - t735 * t815 - t841;
t894 = qJD(2) * pkin(2);
t694 = pkin(7) * t851 - t894;
t618 = t675 * pkin(3) + qJD(4) + t694;
t565 = pkin(4) * t784 + t618;
t499 = t556 * pkin(5) - qJ(6) * t924 + t565;
t942 = t499 * t556;
t941 = t556 * t565;
t667 = t732 * t738 + t735 * t893;
t771 = t739 * t667;
t630 = qJD(1) * t771;
t769 = qJD(3) * t667;
t931 = t630 - t769;
t783 = t732 * t735 - t738 * t893;
t772 = t739 * t783;
t631 = qJD(1) * t772;
t911 = qJD(3) * t783;
t939 = t631 - t911;
t514 = pkin(5) * t924 + qJ(6) * t556;
t918 = t732 * t944 + t893 * t945;
t917 = t732 * t945 - t893 * t944;
t689 = t733 * t735;
t690 = t733 * t738;
t616 = t893 * t689 + t690 * t732;
t590 = -pkin(9) * t667 + t616;
t617 = t732 * t689 - t893 * t690;
t591 = -pkin(9) * t783 + t617;
t542 = t734 * t590 + t591 * t906;
t729 = qJ(3) + pkin(10);
t720 = qJ(5) + t729;
t710 = sin(t720);
t727 = g(3) * t739;
t737 = sin(qJ(1));
t740 = cos(qJ(1));
t803 = g(1) * t740 + g(2) * t737;
t767 = t803 * t736 - t727;
t938 = t542 * t663 + t710 * t767;
t934 = t565 * t924;
t933 = pkin(4) * t851 + pkin(9) * t939 - t918;
t932 = pkin(9) * t931 + t917;
t654 = t663 * pkin(5);
t910 = qJDD(6) - t654;
t930 = t499 * t924 + t910;
t929 = pkin(9) * t605;
t928 = pkin(9) * t784;
t761 = t906 * t783;
t866 = t667 * t840 - (-qJD(3) - qJD(5)) * t761 - t631 * t906 - t931 * t734;
t774 = t734 * t783;
t818 = qJD(5) * t906;
t824 = t906 * t667;
t865 = qJD(3) * t824 - qJD(5) * t774 - t630 * t906 + t667 * t818 + t734 * t939;
t716 = pkin(7) * t850;
t830 = pkin(3) * t880;
t844 = qJD(3) * t735;
t926 = pkin(3) * t844 - qJD(1) * t830 - t716;
t876 = t736 * t740;
t878 = t736 * t737;
t925 = g(1) * t876 + g(2) * t878 - t727;
t688 = -pkin(2) * t739 - pkin(8) * t736 - pkin(1);
t664 = t688 * qJD(1);
t695 = qJD(2) * pkin(8) + t716;
t612 = t738 * t664 - t695 * t735;
t579 = -qJ(4) * t677 + t612;
t613 = t735 * t664 + t738 * t695;
t580 = -qJ(4) * t675 + t613;
t814 = t893 * t580;
t534 = -t579 * t732 - t814;
t525 = t534 + t928;
t575 = t732 * t580;
t535 = t893 * t579 - t575;
t526 = t535 - t929;
t823 = t893 * pkin(3);
t712 = t823 + pkin(4);
t905 = pkin(3) * t732;
t827 = t734 * t905;
t922 = -qJD(5) * t827 - t734 * t525 - t526 * t906 + t712 * t818;
t787 = t590 * t906 - t734 * t591;
t921 = qJD(5) * t787 - t734 * t933 + t906 * t932;
t920 = qJD(5) * t542 + t734 * t932 + t906 * t933;
t669 = t738 * t688;
t903 = pkin(7) * t735;
t609 = -qJ(4) * t877 + t669 + (-pkin(3) - t903) * t739;
t707 = pkin(7) * t874;
t857 = t735 * t688 + t707;
t881 = t735 * t736;
t615 = -qJ(4) * t881 + t857;
t562 = t893 * t609 - t615 * t732;
t643 = t783 * t736;
t538 = -pkin(4) * t739 + pkin(9) * t643 + t562;
t563 = t732 * t609 + t893 * t615;
t773 = t736 * t667;
t543 = -pkin(9) * t773 + t563;
t919 = t734 * t538 + t906 * t543;
t864 = -pkin(4) * t931 + t926;
t650 = t663 * qJ(6);
t687 = t696 * qJD(6);
t916 = t650 - t687;
t714 = pkin(7) * t833;
t816 = t739 * t835;
t891 = qJDD(2) * pkin(2);
t652 = pkin(7) * t816 + t714 - t891;
t915 = qJD(3) * pkin(8) * t705 - t652;
t858 = t734 * t712 + t906 * t905;
t875 = t737 * t739;
t655 = t735 * t875 + t738 * t740;
t873 = t739 * t740;
t657 = -t735 * t873 + t737 * t738;
t912 = -g(1) * t657 + g(2) * t655;
t680 = t806 * qJD(2);
t847 = qJD(2) * t736;
t860 = t738 * t680 + t847 * t903;
t551 = -t736 * t841 + t795 * qJD(2) + (-t707 + (qJ(4) * t736 - t688) * t735) * qJD(3) + t860;
t842 = qJD(3) * t738;
t861 = t735 * t680 + t688 * t842;
t561 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t877 + (-qJD(4) * t736 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t739) * t735 + t861;
t511 = t893 * t551 - t732 * t561;
t760 = t736 * t769;
t501 = pkin(9) * t760 + (t736 * pkin(4) + pkin(9) * t772) * qJD(2) + t511;
t512 = t732 * t551 + t893 * t561;
t748 = -qJD(2) * t771 + t736 * t911;
t504 = pkin(9) * t748 + t512;
t909 = -qJD(5) * t919 + t501 * t906 - t734 * t504;
t908 = -0.2e1 * pkin(1);
t901 = g(1) * t737;
t897 = g(2) * t740;
t896 = g(3) * t736;
t895 = t735 * pkin(3);
t722 = t736 * pkin(7);
t724 = t738 * pkin(3);
t713 = t724 + pkin(2);
t571 = -pkin(3) * t705 + t579;
t529 = t893 * t571 - t575;
t515 = -pkin(4) * t705 + t529 - t929;
t530 = t732 * t571 + t814;
t519 = t530 - t928;
t489 = t734 * t515 + t519 * t906;
t890 = t489 * t696;
t884 = t677 * t705;
t883 = t710 * t736;
t711 = cos(t720);
t882 = t711 * t736;
t879 = t735 * t740;
t872 = t740 * t710;
t871 = -qJ(6) * t851 + t921;
t870 = pkin(5) * t851 + t920;
t602 = t824 - t774;
t869 = pkin(5) * t865 + qJ(6) * t866 - t602 * qJD(6) + t864;
t614 = qJD(1) * t680 + qJDD(1) * t688;
t608 = t738 * t614;
t651 = pkin(7) * t913 + qJDD(2) * pkin(8);
t831 = t735 * qJDD(2);
t832 = qJDD(1) * t738;
t513 = -t735 * t651 + t608 - (t736 * t832 + t738 * t816 + t831) * qJ(4) - t677 * qJD(4) + t666 * pkin(3) - t580 * qJD(3);
t780 = t735 * t614 + t738 * t651 + t664 * t842 - t695 * t844;
t521 = -qJ(4) * t597 - qJD(4) * t675 + t780;
t487 = t732 * t513 + t893 * t521;
t863 = qJD(6) + t922;
t862 = qJD(5) * t858 + t525 * t906 - t734 * t526;
t682 = pkin(4) * cos(t729) + t724;
t855 = pkin(3) * t881 + t722;
t854 = t740 * pkin(1) + t737 * pkin(7);
t730 = t736 ^ 2;
t853 = -t739 ^ 2 + t730;
t849 = qJD(2) * t675;
t846 = qJD(2) * t739;
t845 = qJD(3) * t675;
t843 = qJD(3) * t736;
t839 = t677 * qJD(2);
t838 = t694 * qJD(3);
t488 = t515 * t906 - t734 * t519;
t836 = qJD(6) - t488;
t825 = pkin(7) * t846 + qJD(2) * t830 + t842 * t904;
t821 = t705 * t842;
t820 = t705 * t844;
t819 = t735 * t843;
t624 = t710 * t875 + t711 * t740;
t625 = t711 * t875 - t872;
t813 = -t624 * pkin(5) + qJ(6) * t625;
t626 = -t737 * t711 + t739 * t872;
t627 = t710 * t737 + t711 * t873;
t812 = -t626 * pkin(5) + qJ(6) * t627;
t486 = t893 * t513 - t732 * t521;
t811 = -qJD(3) * t664 - t651;
t478 = t666 * pkin(4) - pkin(9) * t746 + t486;
t482 = pkin(9) * t747 + t487;
t809 = t734 * t478 + t906 * t482 + t515 * t818 - t519 * t840;
t808 = -t906 * t478 + t734 * t482 + t515 * t840 + t519 * t818;
t708 = g(1) * t878;
t807 = -g(2) * t876 + t708;
t572 = pkin(3) * t677 + pkin(4) * t605;
t805 = -g(1) * t624 + g(2) * t626;
t804 = g(1) * t625 - g(2) * t627;
t802 = t695 * t842 - t608;
t801 = -pkin(8) * t666 + t838;
t799 = t713 * t739 - t733 * t736;
t796 = qJD(2) * qJD(3) + t833;
t678 = pkin(2) + t682;
t794 = pkin(5) * t711 + qJ(6) * t710 + t678;
t792 = -pkin(7) * qJDD(2) + t835 * t908;
t790 = t538 * t906 - t734 * t543;
t786 = t666 * t735 - t821;
t785 = t738 * t666 + t820;
t781 = -t735 * t834 + t832;
t779 = t734 * t501 + t906 * t504 + t538 * t818 - t543 * t840;
t742 = qJD(1) ^ 2;
t778 = pkin(1) * t742 + t803;
t777 = t712 * t906 - t827;
t741 = qJD(2) ^ 2;
t776 = pkin(7) * t741 + qJDD(1) * t908 + t897;
t775 = t787 * t663 + t711 * t925;
t570 = t597 * pkin(3) + qJDD(4) + t652;
t764 = t734 * t773;
t758 = g(1) * t627 + g(2) * t625 + g(3) * t882 - t809;
t757 = g(1) * t626 + g(2) * t624 + g(3) * t883 - t808;
t632 = pkin(4) * t783 - t713;
t756 = t736 * t824;
t610 = pkin(4) * t773 + t855;
t753 = -t488 * t696 + t758;
t752 = t696 * t862 + t757;
t751 = -t757 + t930;
t749 = -qJD(2) * t772 - t760;
t564 = -pkin(4) * t748 + t825;
t522 = -pkin(4) * t747 + t570;
t473 = t493 * pkin(5) + t492 * qJ(6) - qJD(6) * t924 + t522;
t728 = -pkin(9) + t733;
t725 = t740 * pkin(7);
t684 = qJ(6) * t882;
t681 = pkin(4) * sin(t729) + t895;
t658 = t735 * t737 + t738 * t873;
t656 = -t737 * t874 + t879;
t644 = -pkin(5) - t777;
t642 = qJ(6) + t858;
t601 = t667 * t734 + t761;
t585 = -t643 * t906 - t764;
t584 = -t643 * t734 + t756;
t539 = t601 * pkin(5) - t602 * qJ(6) + t632;
t527 = t584 * pkin(5) - t585 * qJ(6) + t610;
t524 = -qJD(5) * t764 - t643 * t818 + t734 * t749 - t748 * t906;
t523 = qJD(5) * t756 - t643 * t840 - t734 * t748 - t749 * t906;
t506 = t739 * pkin(5) - t790;
t505 = -qJ(6) * t739 + t919;
t503 = t514 + t572;
t485 = -t696 * qJ(6) + t489;
t484 = t696 * pkin(5) + t836;
t479 = t524 * pkin(5) + t523 * qJ(6) - t585 * qJD(6) + t564;
t475 = -pkin(5) * t847 - t909;
t474 = qJ(6) * t847 - qJD(6) * t739 + t779;
t472 = t808 + t910;
t471 = t809 + t916;
t1 = [(-t489 * t847 - t610 * t492 + t522 * t585 - t565 * t523 + t564 * t924 - t663 * t919 + t696 * t779 + t739 * t809 + t805) * MDP(26) + (-t492 * t585 - t523 * t924) * MDP(20) + (t492 * t584 - t493 * t585 + t523 * t556 - t524 * t924) * MDP(21) + (t492 * t739 + t523 * t696 + t585 * t663 + t847 * t924) * MDP(22) + (-t471 * t739 - t473 * t585 - t474 * t696 - t479 * t924 + t485 * t847 + t492 * t527 + t499 * t523 + t505 * t663 - t805) * MDP(29) + (-t471 * t584 + t472 * t585 - t474 * t556 + t475 * t924 - t484 * t523 - t485 * t524 - t492 * t506 - t493 * t505 + t807) * MDP(28) + ((-t675 * t738 - t677 * t735) * t846 + ((-t677 * qJD(3) - t597) * t738 + (-t750 + t845) * t735) * t736) * MDP(12) + qJDD(1) * MDP(1) + (-t512 * t784 + t563 * t747 - t487 * t773 + t530 * t748 - t511 * t605 - t562 * t746 + t486 * t643 + t529 * (t667 * t843 + t783 * t846) + t807) * MDP(18) + (-t677 * t819 + (t739 * t839 + (t796 + t816) * t877 - t736 * t782) * t738) * MDP(11) + (qJDD(2) * t736 + t739 * t741) * MDP(6) + (qJDD(2) * t739 - t736 * t741) * MDP(7) + (-t897 + t901) * MDP(2) + (t792 * t736 + (-t776 + t901) * t739) * MDP(9) + (t487 * t563 + t530 * t512 + t486 * t562 + t529 * t511 + t570 * t855 + t618 * t825 - g(1) * (pkin(3) * t879 + t725) - g(2) * (t713 * t873 - t733 * t876 + t854) + (-g(1) * (-pkin(1) - t799) - g(2) * t895) * t737) * MDP(19) + (t471 * t505 + t485 * t474 + t473 * t527 + t499 * t479 + t472 * t506 + t484 * t475 - g(1) * (-pkin(5) * t625 - qJ(6) * t624 + t681 * t740 + t725) - g(2) * (pkin(5) * t627 + qJ(6) * t626 + t678 * t873 - t728 * t876 + t854) + (-g(1) * (-t678 * t739 + t728 * t736 - pkin(1)) - g(2) * t681) * t737) * MDP(30) + 0.2e1 * (t719 * t736 - t835 * t853) * MDP(5) + (-(-t688 * t844 + t860) * t705 + t669 * t666 - g(1) * t656 - g(2) * t658 + ((t821 + t849) * pkin(7) + (-pkin(7) * t666 + qJD(2) * t694 - t811) * t735 + t802) * t739 + (pkin(7) * t597 + qJD(2) * t612 + t652 * t735 + t738 * t838) * t736) * MDP(16) + (t861 * t705 - t857 * t666 - g(1) * t655 - g(2) * t657 + (t694 * t837 + (-t820 + t839) * pkin(7) + t780) * t739 + (-t735 * t838 - t613 * qJD(2) + t652 * t738 + (t831 + t781 * t736 + (-t705 + t810) * t837) * pkin(7)) * t736) * MDP(17) + ((t705 * t848 + t597) * t739 + (-t786 - t849) * t736) * MDP(14) + (-t666 * t739 - t705 * t847) * MDP(15) + (-t663 * t739 - t696 * t847) * MDP(24) + (t493 * t739 + t524 * t696 - t556 * t847 - t584 * t663) * MDP(23) + (t472 * t739 + t473 * t584 + t475 * t696 + t479 * t556 - t484 * t847 + t493 * t527 + t499 * t524 - t506 * t663 + t804) * MDP(27) + ((-t831 + (-t705 - t810) * t837) * t739 + (-t739 * t781 + t785 + t839) * t736) * MDP(13) + (qJDD(1) * t730 + 0.2e1 * t736 * t816) * MDP(4) + t803 * MDP(3) + (t736 * t776 + t739 * t792 - t708) * MDP(10) + (t488 * t847 + t610 * t493 + t522 * t584 + t565 * t524 + t564 * t556 + t790 * t663 - t696 * t909 + t808 * t739 + t804) * MDP(25); (-t471 * t601 + t472 * t602 - t484 * t866 - t485 * t865 + t492 * t787 - t493 * t542 - t556 * t871 - t739 * t803 + t870 * t924 - t896) * MDP(28) + (t705 * MDP(15) - MDP(22) * t924 + t556 * MDP(23) + t696 * MDP(24) - t488 * MDP(25) + t489 * MDP(26) + t484 * MDP(27) - t485 * MDP(29)) * t851 + (t492 * t601 - t493 * t602 + t556 * t866 - t865 * t924) * MDP(21) + (-t492 * t602 - t866 * t924) * MDP(20) + (t471 * t542 - t472 * t787 + t473 * t539 + t869 * t499 + t871 * t485 + t870 * t484 + (-g(3) * t794 + t728 * t803) * t739 + (g(3) * t728 + t794 * t803) * t736) * MDP(30) + (t473 * t601 + t493 * t539 + t499 * t865 + t556 * t869 + t696 * t870 + t775) * MDP(27) + (t632 * t493 + t522 * t601 + t864 * t556 + t865 * t565 + t696 * t920 + t775) * MDP(25) + (t602 * t663 + t696 * t866) * MDP(22) + (-t632 * t492 + t522 * t602 - t866 * t565 + t696 * t921 + t864 * t924 - t938) * MDP(26) + (-t473 * t602 + t492 * t539 + t499 * t866 - t696 * t871 - t869 * t924 + t938) * MDP(29) + (-t601 * t663 + t696 * t865) * MDP(23) + qJDD(2) * MDP(8) + MDP(7) * t719 + MDP(6) * t833 + (t896 + (-pkin(7) * qJDD(1) + t778) * t739) * MDP(10) + (-t800 * t735 ^ 2 + (t755 - t884) * t738) * MDP(11) + ((-t597 + t884) * t735 + (-t845 + t831 + t796 * t738 + (-t819 + (t675 + t837) * t739) * qJD(1)) * t738) * MDP(12) + ((t675 * t736 - t705 * t880) * qJD(1) + t785) * MDP(14) + ((-t677 * t736 + t705 * t874) * qJD(1) + t786) * MDP(13) + (-g(1) * t873 - g(2) * t875 - t486 * t667 - t487 * t783 - t529 * t939 + t931 * t530 - t918 * t605 - t616 * t746 + t617 * t747 - t917 * t784 - t896) * MDP(18) + (t736 * t778 - t714 - t727) * MDP(9) + (-MDP(4) * t736 * t739 + MDP(5) * t853) * t742 + (t487 * t617 + t486 * t616 - t570 * t713 - g(3) * t799 + t926 * t618 + t917 * t530 + t918 * t529 + t803 * (t713 * t736 + t733 * t739)) * MDP(19) + (-t659 * t705 + (-pkin(7) * t677 * t739 + t613 * t736) * qJD(1) + (-pkin(2) * t796 + (t705 * t722 + (-t694 - t894) * t739) * qJD(1) + t801) * t738 + (-t891 + t727 + (pkin(2) * t834 - t803) * t736 - t915) * t735) * MDP(17) + (-pkin(2) * t597 + t859 * t705 + t801 * t735 + (-t612 * t736 + (-pkin(7) * t675 - t694 * t735) * t739) * qJD(1) + (t767 + t915) * t738) * MDP(16); t677 * t675 * MDP(11) + (-t675 ^ 2 + t677 ^ 2) * MDP(12) + (-t675 * t705 + t750) * MDP(13) + (-t884 - t597) * MDP(14) + t666 * MDP(15) + (-t613 * t705 - t677 * t694 + (t811 + t896) * t735 - t802 + t912) * MDP(16) + (g(1) * t658 - g(2) * t656 + g(3) * t877 - t612 * t705 + t675 * t694 - t780) * MDP(17) + (-t746 * t823 + t747 * t905 + (-t529 + t535) * t784 + (t530 + t534) * t605) * MDP(18) + (-t529 * t534 - t530 * t535 + (g(3) * t881 + t486 * t893 + t487 * t732 - t618 * t677 + t912) * pkin(3)) * MDP(19) + (-t572 * t556 + t663 * t777 + t752 - t934) * MDP(25) + (-t572 * t924 - t858 * t663 + t696 * t922 + t758 + t941) * MDP(26) + (-t503 * t556 - t644 * t663 + t752 - t930) * MDP(27) + (-t492 * t644 - t493 * t642 + (t485 + t862) * t924 + (t484 - t863) * t556) * MDP(28) + (t503 * t924 + t642 * t663 - t696 * t863 - t758 + t916 - t942) * MDP(29) + (t471 * t642 + t472 * t644 - t499 * t503 - g(1) * (-t681 * t873 + t682 * t737 + t812) - g(2) * (-t681 * t875 - t682 * t740 + t813) - g(3) * (t684 + (-pkin(5) * t710 - t681) * t736) + t863 * t485 + t862 * t484) * MDP(30) + t947; (-t605 ^ 2 - t784 ^ 2) * MDP(18) + (t529 * t605 + t530 * t784 + t570 - t767) * MDP(19) + (-t886 - t936) * MDP(28) + (-t484 * t924 + t485 * t556 + t473 - t925) * MDP(30) + (-MDP(26) + MDP(29)) * (t492 - t885) + (MDP(25) + MDP(27)) * (t493 - t887); (t757 - t890 - t934) * MDP(25) + (t753 + t941) * MDP(26) + (-t514 * t556 + t654 - t751 - t890) * MDP(27) + (pkin(5) * t492 - qJ(6) * t493 + (t485 - t489) * t924 + (t484 - t836) * t556) * MDP(28) + (t514 * t924 + 0.2e1 * t650 - 0.2e1 * t687 - t753 - t942) * MDP(29) + (t471 * qJ(6) - t472 * pkin(5) - t499 * t514 - t484 * t489 - g(1) * t812 - g(2) * t813 - g(3) * (-pkin(5) * t883 + t684) + t836 * t485) * MDP(30) + t947; (-t663 + t946) * MDP(27) + t483 * MDP(28) + (-t696 ^ 2 - t936) * MDP(29) + (t485 * t696 + t751) * MDP(30);];
tau  = t1;
