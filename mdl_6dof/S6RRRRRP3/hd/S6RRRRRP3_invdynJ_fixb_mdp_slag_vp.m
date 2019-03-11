% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:54
% EndTime: 2019-03-10 01:10:11
% DurationCPUTime: 12.84s
% Computational Cost: add. (12440->629), mult. (27220->784), div. (0->0), fcn. (20230->14), ass. (0->288)
t789 = cos(qJ(2));
t941 = pkin(8) + pkin(7);
t741 = t941 * t789;
t724 = qJD(1) * t741;
t784 = sin(qJ(3));
t699 = t784 * t724;
t785 = sin(qJ(2));
t739 = t941 * t785;
t722 = qJD(1) * t739;
t939 = cos(qJ(3));
t642 = -t722 * t939 - t699;
t860 = qJD(3) * t939;
t971 = -pkin(2) * t860 + t642;
t861 = qJD(1) * t939;
t882 = qJD(1) * t785;
t693 = t784 * t882 - t789 * t861;
t907 = t784 * t789;
t695 = -qJD(1) * t907 - t785 * t861;
t633 = -pkin(3) * t695 + pkin(9) * t693;
t619 = pkin(2) * t882 + t633;
t783 = sin(qJ(4));
t788 = cos(qJ(4));
t972 = -t788 * t619 + t783 * t971;
t970 = t783 * t619 + t788 * t971;
t773 = t788 * pkin(10);
t841 = -t695 * pkin(4) + t693 * t773;
t762 = pkin(2) * t784 + pkin(9);
t928 = -pkin(10) - t762;
t856 = qJD(4) * t928;
t969 = -t788 * t856 + t841 - t972;
t626 = t788 * t633;
t927 = qJD(2) * pkin(2);
t703 = -t722 + t927;
t637 = t703 * t939 - t699;
t940 = -pkin(9) - pkin(10);
t866 = qJD(4) * t940;
t968 = -t637 * t783 - t788 * t866 + t626 + t841;
t919 = t693 * t783;
t872 = pkin(10) * t919;
t964 = -t783 * t856 + t872 + t970;
t892 = t783 * t633 + t788 * t637;
t963 = -t783 * t866 + t872 + t892;
t787 = cos(qJ(5));
t782 = sin(qJ(5));
t910 = t782 * t788;
t715 = t783 * t787 + t910;
t873 = qJD(4) + qJD(5);
t645 = t873 * t715;
t962 = t715 * t693 + t645;
t911 = t782 * t783;
t713 = -t787 * t788 + t911;
t877 = qJD(5) * t787;
t879 = qJD(4) * t788;
t893 = t713 * t693 - t787 * t879 - t788 * t877 + t873 * t911;
t775 = qJDD(2) + qJDD(3);
t777 = qJD(2) + qJD(3);
t716 = t785 * t939 + t907;
t864 = t939 * t789;
t820 = -t784 * t785 + t864;
t646 = t777 * t820;
t796 = t646 * qJD(1);
t794 = t716 * qJDD(1) + t796;
t880 = qJD(4) * t783;
t576 = t695 * t880 + t783 * t775 + t777 * t879 + t788 * t794;
t657 = t783 * t695 + t777 * t788;
t658 = -t695 * t788 + t777 * t783;
t869 = -t695 * t879 + t777 * t880 + t783 * t794;
t825 = t788 * t775 - t869;
t878 = qJD(5) * t782;
t523 = -t787 * t576 - t657 * t877 + t658 * t878 - t782 * t825;
t830 = t657 * t782 + t787 * t658;
t524 = qJD(5) * t830 + t576 * t782 - t787 * t825;
t590 = -t787 * t657 + t658 * t782;
t588 = t590 ^ 2;
t647 = t777 * t716;
t875 = qJDD(1) * t785;
t834 = -qJDD(1) * t864 + t784 * t875;
t603 = qJD(1) * t647 + t834;
t602 = qJDD(4) + t603;
t600 = qJDD(5) + t602;
t688 = qJD(4) + t693;
t675 = qJD(5) + t688;
t942 = t830 ^ 2;
t967 = t600 * MDP(29) + (t675 * t830 - t524) * MDP(28) + t590 * MDP(25) * t830 + (t590 * t675 - t523) * MDP(27) + (-t588 + t942) * MDP(26);
t966 = qJ(6) * t590;
t876 = qJD(1) * qJD(2);
t858 = t789 * t876;
t651 = qJDD(2) * pkin(2) + t941 * (-t858 - t875);
t859 = t785 * t876;
t874 = qJDD(1) * t789;
t656 = t941 * (-t859 + t874);
t881 = qJD(3) * t784;
t842 = -t939 * t651 + t784 * t656 + t703 * t881 + t724 * t860;
t931 = t775 * pkin(3);
t562 = t842 - t931;
t781 = qJ(2) + qJ(3);
t772 = cos(t781);
t933 = g(3) * t772;
t965 = t562 + t933;
t770 = sin(t781);
t790 = cos(qJ(1));
t914 = t770 * t790;
t786 = sin(qJ(1));
t915 = t770 * t786;
t961 = g(1) * t914 + g(2) * t915;
t923 = t646 * t783;
t960 = t716 * t879 + t923;
t945 = -t933 + t961;
t622 = -t777 * pkin(3) - t637;
t586 = -t657 * pkin(4) + t622;
t780 = qJ(4) + qJ(5);
t769 = sin(t780);
t771 = cos(t780);
t913 = t772 * t786;
t665 = t769 * t790 - t771 * t913;
t912 = t772 * t790;
t667 = t769 * t786 + t771 * t912;
t757 = g(3) * t770;
t756 = pkin(2) * t859;
t930 = t789 * pkin(2);
t765 = pkin(1) + t930;
t951 = -pkin(9) * t716 - t765;
t544 = t603 * pkin(3) - pkin(9) * t796 + qJDD(1) * t951 + t756;
t542 = t788 * t544;
t802 = t784 * t651 + t656 * t939 + t703 * t860 - t724 * t881;
t561 = t775 * pkin(9) + t802;
t737 = t765 * qJD(1);
t617 = pkin(3) * t693 + pkin(9) * t695 - t737;
t702 = t939 * t724;
t638 = t784 * t703 + t702;
t623 = pkin(9) * t777 + t638;
t578 = t617 * t783 + t623 * t788;
t497 = pkin(4) * t602 - pkin(10) * t576 - qJD(4) * t578 - t561 * t783 + t542;
t815 = t783 * t544 + t788 * t561 + t617 * t879 - t623 * t880;
t501 = pkin(10) * t825 + t815;
t577 = t788 * t617 - t623 * t783;
t550 = -pkin(10) * t658 + t577;
t539 = pkin(4) * t688 + t550;
t551 = pkin(10) * t657 + t578;
t843 = -t782 * t497 - t787 * t501 - t539 * t877 + t551 * t878;
t959 = g(1) * t667 - g(2) * t665 + t586 * t590 + t771 * t757 + t843;
t957 = qJ(6) * t830;
t540 = t590 * pkin(5) + qJD(6) + t586;
t956 = t540 * t830;
t636 = -pkin(3) * t820 + t951;
t630 = t788 * t636;
t660 = -t784 * t739 + t741 * t939;
t917 = t716 * t788;
t569 = -pkin(4) * t820 - pkin(10) * t917 - t660 * t783 + t630;
t652 = t788 * t660;
t889 = t783 * t636 + t652;
t918 = t716 * t783;
t580 = -pkin(10) * t918 + t889;
t895 = t782 * t569 + t787 * t580;
t955 = -qJ(6) * t962 - t713 * qJD(6);
t641 = -t784 * t722 + t702;
t840 = pkin(2) * t881 - t641;
t954 = t969 * t787;
t953 = (t880 + t919) * pkin(4);
t706 = t928 * t783;
t707 = t762 * t788 + t773;
t888 = t782 * t706 + t787 * t707;
t774 = t788 * pkin(4);
t732 = pkin(5) * t771 + t774;
t721 = pkin(3) + t732;
t776 = -qJ(6) + t940;
t848 = t772 * t721 - t770 * t776;
t952 = t968 * t787;
t738 = t940 * t783;
t740 = pkin(9) * t788 + t773;
t887 = t782 * t738 + t787 * t740;
t950 = -t738 * t877 + t740 * t878 + t782 * t968 + t787 * t963;
t949 = -t706 * t877 + t707 * t878 + t782 * t969 + t964 * t787;
t948 = t695 * pkin(5) + qJ(6) * t893 - qJD(6) * t715;
t947 = -t939 * t739 - t784 * t741;
t838 = g(1) * t790 + g(2) * t786;
t664 = t769 * t913 + t771 * t790;
t666 = -t769 * t912 + t771 * t786;
t946 = -g(1) * t666 + g(2) * t664 + t769 * t757;
t547 = t787 * t551;
t518 = t539 * t782 + t547;
t854 = t787 * t497 - t782 * t501;
t803 = -qJD(5) * t518 + t854;
t944 = -t586 * t830 + t803 + t946;
t932 = t713 * pkin(5);
t929 = pkin(3) + t774;
t926 = qJ(6) * t715;
t925 = t576 * t783;
t924 = t622 * t693;
t922 = t646 * t788;
t921 = t657 * t688;
t920 = t658 * t688;
t545 = t782 * t551;
t909 = t783 * t786;
t908 = t783 * t790;
t906 = t786 * t788;
t904 = t788 * t790;
t517 = t787 * t539 - t545;
t506 = t517 - t957;
t504 = pkin(5) * t675 + t506;
t903 = -t506 + t504;
t902 = -t949 + t955;
t901 = -qJD(5) * t888 + t782 * t964 + t948 - t954;
t900 = t787 * t550 - t545;
t899 = -t950 + t955;
t898 = -qJD(5) * t887 + t782 * t963 + t948 - t952;
t886 = t953 + t840;
t731 = pkin(4) * t783 + pkin(5) * t769;
t885 = t731 + t941;
t778 = t785 ^ 2;
t884 = -t789 ^ 2 + t778;
t871 = t785 * t927;
t870 = qJD(4) * pkin(9) * t688;
t868 = g(1) * t912 + g(2) * t913 + t757;
t867 = qJD(2) * t941;
t863 = t716 * t880;
t612 = t622 * t879;
t585 = pkin(3) * t647 - pkin(9) * t646 + t871;
t582 = t788 * t585;
t723 = t785 * t867;
t725 = t789 * t867;
t595 = qJD(3) * t947 - t939 * t723 - t784 * t725;
t512 = -pkin(10) * t922 + pkin(4) * t647 - t595 * t783 + t582 + (-t652 + (pkin(10) * t716 - t636) * t783) * qJD(4);
t814 = t783 * t585 + t788 * t595 + t636 * t879 - t660 * t880;
t520 = -pkin(10) * t960 + t814;
t853 = t787 * t512 - t520 * t782;
t852 = -t550 * t782 - t547;
t851 = t787 * t569 - t580 * t782;
t849 = t787 * t706 - t707 * t782;
t847 = t787 * t738 - t740 * t782;
t846 = t688 * t788;
t845 = -qJD(4) * t617 - t561;
t764 = -pkin(2) * t939 - pkin(3);
t839 = -t638 + t953;
t837 = g(1) * t786 - g(2) * t790;
t836 = -t623 * t879 + t542;
t835 = -pkin(9) * t602 + t924;
t736 = -t774 + t764;
t831 = -t602 * t762 + t924;
t829 = t721 * t770 + t772 * t776;
t827 = -t578 * t695 + t783 * t965 + t612;
t826 = t577 * t695 + t622 * t880 + t788 * t961;
t618 = pkin(4) * t918 - t947;
t824 = pkin(5) * t962 + t953;
t823 = t838 * t770;
t822 = -0.2e1 * pkin(1) * t876 - pkin(7) * qJDD(2);
t821 = t765 + t848;
t818 = -t863 + t922;
t596 = -t784 * t723 + t725 * t939 - t739 * t881 + t741 * t860;
t813 = t782 * t512 + t787 * t520 + t569 * t877 - t580 * t878;
t556 = pkin(4) * t960 + t596;
t792 = qJD(2) ^ 2;
t808 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t792 + t837;
t793 = qJD(1) ^ 2;
t807 = pkin(1) * t793 - pkin(7) * qJDD(1) + t838;
t805 = -t737 * t695 - t842 + t945;
t489 = pkin(5) * t600 + qJ(6) * t523 - qJD(6) * t830 + t803;
t491 = -qJ(6) * t524 - qJD(6) * t590 - t843;
t507 = t518 - t966;
t801 = -t489 * t715 - t491 * t713 + t504 * t893 - t507 * t962 - t868;
t529 = -pkin(4) * t825 + t562;
t800 = t517 * t695 + t529 * t713 + t586 * t962 + t771 * t945;
t799 = -t518 * t695 + t529 * t715 + (-t823 + t933) * t769 - t893 * t586;
t499 = t524 * pkin(5) + qJDD(6) + t529;
t797 = -t737 * t693 - t802 + t868;
t795 = (t523 * t713 - t524 * t715 + t590 * t893 - t830 * t962) * MDP(26) + (-t523 * t715 - t830 * t893) * MDP(25) + ((t576 + t921) * t788 + (t825 - t920) * t783) * MDP(19) + (t600 * t715 - t675 * t893 + t695 * t830) * MDP(27) + (-t590 * t695 - t600 * t713 - t675 * t962) * MDP(28) + (t658 * t846 + t925) * MDP(18) + (-t688 ^ 2 * t783 + t602 * t788 + t657 * t695) * MDP(21) + (t602 * t783 + t658 * t695 + t688 * t846) * MDP(20) + (t693 * t777 + t794) * MDP(13) + (-t834 + (-qJD(1) * t716 - t695) * t777) * MDP(14) + (-t693 ^ 2 + t695 ^ 2) * MDP(12) + t775 * MDP(15) + (-MDP(11) * t693 + MDP(22) * t688 + MDP(29) * t675) * t695;
t763 = pkin(4) * t787 + pkin(5);
t705 = t713 * qJ(6);
t690 = -qJDD(1) * t765 + t756;
t684 = t772 * t904 + t909;
t683 = -t772 * t908 + t906;
t682 = -t772 * t906 + t908;
t681 = t772 * t909 + t904;
t628 = t713 * t716;
t627 = t715 * t716;
t621 = -t705 + t887;
t620 = t847 - t926;
t607 = -t705 + t888;
t606 = t849 - t926;
t534 = t646 * t910 - t782 * t863 - t878 * t918 + (t873 * t917 + t923) * t787;
t533 = t645 * t716 + t646 * t713;
t526 = -qJ(6) * t627 + t895;
t525 = -pkin(5) * t820 + qJ(6) * t628 + t851;
t509 = t900 - t957;
t508 = t852 + t966;
t493 = -qJ(6) * t534 - qJD(6) * t627 + t813;
t492 = pkin(5) * t647 + qJ(6) * t533 - qJD(5) * t895 + qJD(6) * t628 + t853;
t1 = [(-g(1) * t681 - g(2) * t683 + t562 * t917 - t576 * t947 - t578 * t647 + t596 * t658 - t602 * t889 + t622 * t818 - t688 * t814 + t815 * t820) * MDP(24) + ((-t660 * t879 + t582) * t688 + t630 * t602 - t836 * t820 + t577 * t647 - t596 * t657 + t947 * t825 + t716 * t612 - g(1) * t682 - g(2) * t684 + ((-qJD(4) * t636 - t595) * t688 - t660 * t602 - t845 * t820 + t562 * t716 + t622 * t646) * t783) * MDP(23) + (-t596 * t777 - t603 * t765 - t647 * t737 - t690 * t820 + t693 * t871 + t772 * t837 + t775 * t947) * MDP(16) + (t853 * t675 + t851 * t600 - t854 * t820 + t517 * t647 + t556 * t590 + t618 * t524 + t529 * t627 + t586 * t534 - g(1) * t665 - g(2) * t667 + (t518 * t820 - t675 * t895) * qJD(5)) * MDP(30) + (-t602 * t918 + t657 * t647 - t688 * t960 - t820 * t825) * MDP(21) + (-t600 * t820 + t647 * t675) * MDP(29) + (-t647 * t777 + t775 * t820) * MDP(14) + (t523 * t820 - t533 * t675 - t600 * t628 + t647 * t830) * MDP(27) + (-g(1) * t664 - g(2) * t666 - t518 * t647 - t618 * t523 - t529 * t628 - t586 * t533 + t556 * t830 - t600 * t895 - t675 * t813 - t820 * t843) * MDP(31) + (-t576 * t820 + t602 * t917 + t647 * t658 + t688 * t818) * MDP(20) + (-t716 * t603 - t646 * t693 + t695 * t647 + t794 * t820) * MDP(12) + (qJDD(2) * t785 + t789 * t792) * MDP(6) + (qJDD(2) * t789 - t785 * t792) * MDP(7) + (-g(1) * t915 + g(2) * t914 - t595 * t777 - t737 * t646 - t660 * t775 + t690 * t716 - t695 * t871 - t765 * t794) * MDP(17) + (-t602 * t820 + t647 * t688) * MDP(22) + (t524 * t820 - t534 * t675 - t590 * t647 - t600 * t627) * MDP(28) + (t646 * t777 + t716 * t775) * MDP(13) + (t523 * t627 + t524 * t628 + t533 * t590 - t534 * t830) * MDP(26) + (t523 * t628 - t533 * t830) * MDP(25) + (t489 * t628 - t491 * t627 - t492 * t830 - t493 * t590 + t504 * t533 - t507 * t534 + t523 * t525 - t524 * t526 + t770 * t837) * MDP(32) + (-t695 * t646 + t716 * t794) * MDP(11) + ((t657 * t788 - t658 * t783) * t646 + (t788 * t825 - t925 + (-t657 * t783 - t658 * t788) * qJD(4)) * t716) * MDP(19) + (t576 * t917 + t658 * t818) * MDP(18) + (t491 * t526 + t507 * t493 + t489 * t525 + t504 * t492 + t499 * (t627 * pkin(5) + t618) + t540 * (t534 * pkin(5) + t556) + (-g(1) * t885 - g(2) * t821) * t790 + (g(1) * t821 - g(2) * t885) * t786) * MDP(33) + 0.2e1 * (t785 * t874 - t876 * t884) * MDP(5) + qJDD(1) * MDP(1) + (qJDD(1) * t778 + 0.2e1 * t785 * t858) * MDP(4) + (t785 * t822 + t789 * t808) * MDP(9) + (-t785 * t808 + t789 * t822) * MDP(10) + t837 * MDP(2) + t838 * MDP(3); (t764 * t576 + t831 * t788 - t783 * t823 + t840 * t658 + (t762 * t880 + t970) * t688 + t827) * MDP(24) + t795 + (t523 * t606 - t524 * t607 - t590 * t902 - t830 * t901 + t801) * MDP(32) + (-t736 * t523 - t888 * t600 + t675 * t949 + t886 * t830 + t799) * MDP(31) + (t764 * t869 + (-t764 * t775 - t965) * t788 + t831 * t783 - t840 * t657 + (-t762 * t879 + t972) * t688 + t826) * MDP(23) + (t491 * t607 + t489 * t606 + t499 * (t736 + t932) - g(3) * (t848 + t930) + (t824 + t840) * t540 + t902 * t507 + t901 * t504 + t838 * (pkin(2) * t785 + t829)) * MDP(33) + (t642 * t777 + (t695 * t882 - t775 * t784 - t777 * t860) * pkin(2) + t797) * MDP(17) + (t849 * t600 + t736 * t524 + (-t707 * t877 + (-qJD(5) * t706 + t964) * t782 - t954) * t675 + t886 * t590 + t800) * MDP(30) + (-g(3) * t789 + t785 * t807) * MDP(9) + (g(3) * t785 + t789 * t807) * MDP(10) + (t641 * t777 + (-t693 * t882 + t775 * t939 - t777 * t881) * pkin(2) + t805) * MDP(16) + qJDD(2) * MDP(8) + MDP(7) * t874 + MDP(6) * t875 + (-MDP(4) * t785 * t789 + MDP(5) * t884) * t793; t795 + (t637 * t777 + t797) * MDP(17) + (t523 * t620 - t524 * t621 - t590 * t899 - t830 * t898 + t801) * MDP(32) + (t523 * t929 - t887 * t600 + t675 * t950 + t839 * t830 + t799) * MDP(31) + (-pkin(3) * t869 - t626 * t688 + t638 * t657 + (t637 * t688 + t835) * t783 + (-t965 - t870 + t931) * t788 + t826) * MDP(23) + (t491 * t621 + t489 * t620 + t499 * (-t929 + t932) - g(3) * t848 + (t824 - t638) * t540 + t899 * t507 + t898 * t504 + t838 * t829) * MDP(33) + (-pkin(3) * t576 + t892 * t688 - t638 * t658 + t835 * t788 + (-t823 + t870) * t783 + t827) * MDP(24) + (t638 * t777 + t805) * MDP(16) + (t847 * t600 - t929 * t524 + (-t740 * t877 + (-qJD(5) * t738 + t963) * t782 - t952) * t675 + t839 * t590 + t800) * MDP(30); -t658 * t657 * MDP(18) + (-t657 ^ 2 + t658 ^ 2) * MDP(19) + (t576 - t921) * MDP(20) + (t825 + t920) * MDP(21) + t602 * MDP(22) + (-g(1) * t683 + g(2) * t681 + t578 * t688 - t622 * t658 + (t845 + t757) * t783 + t836) * MDP(23) + (g(1) * t684 - g(2) * t682 + t577 * t688 - t622 * t657 + t757 * t788 - t815) * MDP(24) + (-t852 * t675 + (-t590 * t658 + t600 * t787 - t675 * t878) * pkin(4) + t944) * MDP(30) + (t900 * t675 + (-t600 * t782 - t658 * t830 - t675 * t877) * pkin(4) + t959) * MDP(31) + (-t504 * t590 + t507 * t830 + t508 * t830 + t509 * t590 + t523 * t763 + (-t524 * t782 + (-t590 * t787 + t782 * t830) * qJD(5)) * pkin(4)) * MDP(32) + (t489 * t763 - t507 * t509 - t504 * t508 - pkin(5) * t956 - g(1) * (-t731 * t912 + t732 * t786) - g(2) * (-t731 * t913 - t732 * t790) + t731 * t757 + (t491 * t782 - t540 * t658 + (-t504 * t782 + t507 * t787) * qJD(5)) * pkin(4)) * MDP(33) + t967; (t518 * t675 + t944) * MDP(30) + (t517 * t675 + t959) * MDP(31) + (pkin(5) * t523 - t590 * t903) * MDP(32) + (t903 * t507 + (t489 + t946 - t956) * pkin(5)) * MDP(33) + t967; (-t588 - t942) * MDP(32) + (t504 * t830 + t507 * t590 + t499 - t945) * MDP(33);];
tau  = t1;
