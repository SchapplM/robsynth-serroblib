% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PRPRR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:57
% EndTime: 2019-12-05 15:48:09
% DurationCPUTime: 12.47s
% Computational Cost: add. (42082->538), mult. (78250->770), div. (0->0), fcn. (53006->10), ass. (0->364)
t909 = sin(qJ(2));
t912 = cos(qJ(2));
t903 = sin(pkin(8));
t905 = cos(pkin(8));
t962 = t905 * g(1) + t903 * g(2);
t964 = g(3) - qJDD(1);
t842 = -t909 * t964 - t912 * t962;
t986 = qJD(2) ^ 2;
t834 = -t986 * pkin(2) + t842;
t902 = sin(pkin(9));
t904 = cos(pkin(9));
t841 = -t909 * t962 + t912 * t964;
t914 = qJDD(2) * pkin(2) - t841;
t763 = t902 * t834 - t904 * t914;
t764 = t904 * t834 + t902 * t914;
t942 = t902 * t763 + t904 * t764;
t694 = t904 * t763 - t902 * t764;
t965 = t912 * t694;
t640 = -t909 * t942 + t965;
t973 = t909 * t694;
t1001 = t912 * t942 + t973;
t876 = t903 * g(1) - t905 * g(2);
t868 = -qJDD(3) + t876;
t869 = t902 * qJDD(2) + t904 * t986;
t811 = qJ(3) * t869 - t904 * t868;
t870 = t904 * qJDD(2) - t902 * t986;
t927 = -qJ(3) * t870 - t902 * t868;
t939 = -t909 * t869 + t912 * t870;
t1000 = -pkin(5) * t939 + t909 * t811 + t912 * t927;
t907 = sin(qJ(5));
t910 = cos(qJ(5));
t911 = cos(qJ(4));
t908 = sin(qJ(4));
t961 = qJD(2) * t908;
t847 = -t910 * t911 * qJD(2) + t907 * t961;
t849 = (t907 * t911 + t908 * t910) * qJD(2);
t798 = t849 * t847;
t898 = qJDD(4) + qJDD(5);
t991 = -t798 + t898;
t999 = t907 * t991;
t998 = t910 * t991;
t801 = t912 * t869 + t909 * t870;
t716 = pkin(5) * t801 + t912 * t811 - t909 * t927;
t996 = t903 * t964;
t995 = t905 * t964;
t761 = -t986 * pkin(3) + qJDD(2) * pkin(6) + t764;
t731 = t908 * t761 + t911 * t868;
t732 = t911 * t761 - t908 * t868;
t673 = t908 * t731 + t911 * t732;
t958 = qJD(2) * qJD(4);
t888 = t911 * t958;
t956 = t908 * qJDD(2);
t865 = t888 + t956;
t949 = t908 * t958;
t955 = t911 * qJDD(2);
t917 = t949 - t955;
t771 = -t847 * qJD(5) + t910 * t865 - t907 * t917;
t899 = qJD(4) + qJD(5);
t844 = t899 * t847;
t990 = -t844 + t771;
t856 = t905 * t876;
t989 = t903 * t962 - t856;
t940 = t907 * t865 + t910 * t917;
t742 = (qJD(5) - t899) * t849 + t940;
t845 = t847 ^ 2;
t846 = t849 ^ 2;
t897 = t899 ^ 2;
t987 = t911 ^ 2;
t985 = pkin(2) * t694;
t884 = t908 * t986 * t911;
t877 = qJDD(4) + t884;
t697 = (-t865 + t888) * pkin(7) + t877 * pkin(4) - t731;
t879 = qJD(4) * pkin(4) - pkin(7) * t961;
t893 = t987 * t986;
t698 = -pkin(4) * t893 - pkin(7) * t917 - qJD(4) * t879 + t732;
t649 = -t910 * t697 + t907 * t698;
t650 = t907 * t697 + t910 * t698;
t605 = -t910 * t649 + t907 * t650;
t984 = pkin(4) * t605;
t746 = t844 + t771;
t682 = -t742 * t907 - t910 * t746;
t983 = pkin(4) * t682;
t982 = t849 * t899;
t981 = t899 * t907;
t980 = t899 * t910;
t979 = t903 * t876;
t760 = -qJDD(2) * pkin(3) - t986 * pkin(6) + t763;
t718 = pkin(4) * t917 - pkin(7) * t893 + t879 * t961 + t760;
t978 = t907 * t718;
t793 = t798 + t898;
t977 = t907 * t793;
t976 = t908 * t605;
t751 = t908 * t760;
t975 = t908 * t877;
t878 = qJDD(4) - t884;
t974 = t908 * t878;
t970 = t910 * t718;
t969 = t910 * t793;
t968 = t911 * t605;
t752 = t911 * t760;
t967 = t911 * t877;
t966 = t911 * t878;
t963 = -pkin(3) * t760 + pkin(6) * t673;
t900 = t908 ^ 2;
t960 = t900 * t986;
t957 = t905 * qJDD(2);
t954 = t900 + t987;
t953 = t902 * t798;
t952 = t904 * t798;
t913 = qJD(4) ^ 2;
t881 = -t913 - t960;
t831 = -t908 * t881 - t966;
t864 = 0.2e1 * t888 + t956;
t951 = -pkin(3) * t864 + pkin(6) * t831 + t751;
t883 = -t893 - t913;
t829 = t911 * t883 - t975;
t866 = -0.2e1 * t949 + t955;
t950 = pkin(3) * t866 + pkin(6) * t829 - t752;
t871 = t954 * qJDD(2);
t874 = t893 + t960;
t809 = t902 * t871 + t904 * t874;
t812 = t904 * t871 - t902 * t874;
t749 = t912 * t809 + t909 * t812;
t933 = pkin(3) * t874 + pkin(6) * t871 + t673;
t918 = pkin(2) * t809 + t933;
t629 = -pkin(1) * t749 - t918;
t750 = -t909 * t809 + t912 * t812;
t947 = qJ(1) * t750 + t629;
t929 = -pkin(2) * t869 - t764;
t705 = pkin(1) * t801 - t929;
t946 = -qJ(1) * t939 + t705;
t924 = pkin(2) * t870 - t763;
t706 = -pkin(1) * t939 - t924;
t945 = -qJ(1) * t801 + t706;
t872 = t909 * qJDD(2) + t912 * t986;
t796 = pkin(1) * t872 + t842;
t873 = t912 * qJDD(2) - t909 * t986;
t944 = qJ(1) * t873 - t796;
t797 = -pkin(1) * t873 + t841;
t943 = qJ(1) * t872 - t797;
t606 = t907 * t649 + t910 * t650;
t941 = t909 * t841 + t912 * t842;
t937 = -t905 * t962 - t979;
t684 = -t742 * t910 + t907 * t746;
t772 = -t845 - t846;
t584 = -pkin(4) * t772 + pkin(7) * t684 + t606;
t590 = -pkin(7) * t682 - t605;
t628 = -t908 * t682 + t911 * t684;
t936 = -pkin(3) * t772 + pkin(6) * t628 + t911 * t584 + t908 * t590;
t791 = -t897 - t845;
t734 = t910 * t791 - t999;
t741 = (qJD(5) + t899) * t849 + t940;
t635 = -pkin(4) * t741 + pkin(7) * t734 - t970;
t733 = t907 * t791 + t998;
t666 = -pkin(7) * t733 + t978;
t675 = -t908 * t733 + t911 * t734;
t935 = -pkin(3) * t741 + pkin(6) * t675 + t911 * t635 + t908 * t666;
t832 = -t846 - t897;
t748 = -t907 * t832 - t969;
t645 = -pkin(4) * t990 + pkin(7) * t748 + t978;
t747 = t910 * t832 - t977;
t671 = -pkin(7) * t747 + t970;
t686 = -t908 * t747 + t911 * t748;
t934 = -pkin(3) * t990 + pkin(6) * t686 + t911 * t645 + t908 * t671;
t932 = t902 * t884;
t931 = t904 * t884;
t646 = t902 * t673 - t904 * t760;
t930 = pkin(2) * t646 + t963;
t821 = pkin(5) * t872 - t912 * t876;
t928 = -pkin(5) * t873 - t909 * t876;
t672 = t911 * t731 - t908 * t732;
t775 = t912 * t841 - t909 * t842;
t780 = t902 * t831 - t904 * t864;
t926 = pkin(2) * t780 + t951;
t779 = t902 * t829 + t904 * t866;
t925 = pkin(2) * t779 + t950;
t923 = pkin(4) * t733 - t649;
t577 = t911 * t606 - t976;
t594 = -pkin(4) * t718 + pkin(7) * t606;
t922 = -pkin(3) * t718 + pkin(6) * t577 - pkin(7) * t976 + t911 * t594;
t620 = t902 * t628 - t904 * t772;
t921 = pkin(2) * t620 + t936;
t642 = t902 * t675 - t904 * t741;
t920 = pkin(2) * t642 + t935;
t653 = t902 * t686 - t904 * t990;
t919 = pkin(2) * t653 + t934;
t916 = pkin(4) * t747 - t650;
t572 = t902 * t577 - t904 * t718;
t915 = pkin(2) * t572 + t922;
t890 = t903 * qJDD(2);
t882 = t893 - t913;
t880 = t913 - t960;
t875 = -t893 + t960;
t862 = t954 * t958;
t840 = -t846 + t897;
t839 = t845 - t897;
t838 = t911 * t865 - t900 * t958;
t837 = t908 * t917 - t987 * t958;
t836 = t902 * qJDD(4) + t904 * t862;
t835 = -t904 * qJDD(4) + t902 * t862;
t830 = -t908 * t880 + t967;
t828 = t911 * t882 - t974;
t827 = t911 * t881 - t974;
t826 = t911 * t880 + t975;
t825 = t908 * t883 + t967;
t824 = t908 * t882 + t966;
t823 = (t865 + t888) * t908;
t822 = t908 * t888 + t911 * t917;
t800 = -t908 * t864 + t911 * t866;
t799 = t911 * t864 + t908 * t866;
t795 = t846 - t845;
t790 = t904 * t838 - t932;
t789 = t904 * t837 + t932;
t788 = t902 * t838 + t931;
t787 = t902 * t837 - t931;
t786 = t904 * t830 + t902 * t956;
t785 = t904 * t828 + t902 * t955;
t784 = t902 * t830 - t904 * t956;
t783 = t902 * t828 - t904 * t955;
t782 = t904 * t831 + t902 * t864;
t781 = t904 * t829 - t902 * t866;
t778 = (-t847 * t910 + t849 * t907) * t899;
t777 = (-t847 * t907 - t849 * t910) * t899;
t770 = -t849 * qJD(5) - t940;
t769 = t904 * t800 + t902 * t875;
t768 = t902 * t800 - t904 * t875;
t766 = -t909 * t835 + t912 * t836;
t765 = t912 * t835 + t909 * t836;
t759 = pkin(1) * t876 + pkin(5) * t941;
t757 = t910 * t839 - t977;
t756 = -t907 * t840 + t998;
t755 = t907 * t839 + t969;
t754 = t910 * t840 + t999;
t738 = t910 * t771 - t849 * t981;
t737 = t907 * t771 + t849 * t980;
t736 = -t907 * t770 + t847 * t980;
t735 = t910 * t770 + t847 * t981;
t730 = -t909 * t788 + t912 * t790;
t729 = -t909 * t787 + t912 * t789;
t728 = t912 * t788 + t909 * t790;
t727 = t912 * t787 + t909 * t789;
t726 = -t909 * t784 + t912 * t786;
t725 = -t909 * t783 + t912 * t785;
t724 = t912 * t784 + t909 * t786;
t723 = t912 * t783 + t909 * t785;
t720 = -pkin(6) * t827 + t752;
t719 = -pkin(6) * t825 + t751;
t712 = -t909 * t780 + t912 * t782;
t711 = -t909 * t779 + t912 * t781;
t710 = t912 * t780 + t909 * t782;
t709 = t912 * t779 + t909 * t781;
t708 = -t908 * t777 + t911 * t778;
t707 = t911 * t777 + t908 * t778;
t704 = -t909 * t768 + t912 * t769;
t703 = t912 * t768 + t909 * t769;
t702 = -pkin(3) * t827 + t732;
t701 = -pkin(3) * t825 + t731;
t700 = t904 * t708 + t902 * t898;
t699 = t902 * t708 - t904 * t898;
t691 = pkin(2) * t868 + qJ(3) * t942;
t690 = -t908 * t755 + t911 * t757;
t689 = -t908 * t754 + t911 * t756;
t688 = t911 * t755 + t908 * t757;
t687 = t911 * t754 + t908 * t756;
t685 = t911 * t747 + t908 * t748;
t683 = -t910 * t741 - t907 * t990;
t681 = -t907 * t741 + t910 * t990;
t679 = -t908 * t737 + t911 * t738;
t678 = -t908 * t735 + t911 * t736;
t677 = t911 * t737 + t908 * t738;
t676 = t911 * t735 + t908 * t736;
t674 = t911 * t733 + t908 * t734;
t664 = t904 * t679 + t953;
t663 = t904 * t678 - t953;
t662 = t902 * t679 - t952;
t661 = t902 * t678 + t952;
t660 = -qJ(3) * t809 + t904 * t672;
t659 = qJ(3) * t812 + t902 * t672;
t658 = t904 * t690 - t902 * t742;
t657 = t904 * t689 + t902 * t746;
t656 = t902 * t690 + t904 * t742;
t655 = t902 * t689 - t904 * t746;
t654 = t904 * t686 + t902 * t990;
t652 = -t909 * t699 + t912 * t700;
t651 = t912 * t699 + t909 * t700;
t647 = t904 * t673 + t902 * t760;
t643 = t904 * t675 + t902 * t741;
t637 = -pkin(1) * t710 - t926;
t636 = -pkin(1) * t709 - t925;
t634 = -qJ(3) * t780 - t902 * t702 + t904 * t720;
t633 = -qJ(3) * t779 - t902 * t701 + t904 * t719;
t631 = -pkin(2) * t827 + qJ(3) * t782 + t904 * t702 + t902 * t720;
t630 = -pkin(2) * t825 + qJ(3) * t781 + t904 * t701 + t902 * t719;
t627 = -t908 * t681 + t911 * t683;
t626 = t911 * t682 + t908 * t684;
t625 = t911 * t681 + t908 * t683;
t623 = t904 * t627 + t902 * t795;
t622 = t902 * t627 - t904 * t795;
t621 = t904 * t628 + t902 * t772;
t619 = pkin(1) * t640 + t985;
t618 = -t909 * t662 + t912 * t664;
t617 = -t909 * t661 + t912 * t663;
t616 = t912 * t662 + t909 * t664;
t615 = t912 * t661 + t909 * t663;
t614 = -t909 * t656 + t912 * t658;
t613 = -t909 * t655 + t912 * t657;
t612 = t912 * t656 + t909 * t658;
t611 = t912 * t655 + t909 * t657;
t610 = -pkin(3) * t626 - t983;
t609 = -pkin(3) * t685 - t916;
t608 = -t909 * t653 + t912 * t654;
t607 = t912 * t653 + t909 * t654;
t604 = -pkin(3) * t674 - t923;
t603 = -pkin(5) * t749 - t909 * t659 + t912 * t660;
t602 = pkin(5) * t750 + t912 * t659 + t909 * t660;
t601 = -t909 * t646 + t912 * t647;
t600 = t912 * t646 + t909 * t647;
t599 = -t909 * t642 + t912 * t643;
t598 = t912 * t642 + t909 * t643;
t597 = pkin(5) * t640 + qJ(3) * t965 - t909 * t691;
t596 = pkin(1) * t868 + pkin(5) * t1001 + qJ(3) * t973 + t912 * t691;
t595 = -pkin(6) * t685 - t908 * t645 + t911 * t671;
t592 = -qJ(3) * t646 - (pkin(3) * t902 - pkin(6) * t904) * t672;
t591 = -pkin(6) * t674 - t908 * t635 + t911 * t666;
t588 = -pkin(5) * t710 - t909 * t631 + t912 * t634;
t587 = -pkin(5) * t709 - t909 * t630 + t912 * t633;
t586 = -pkin(1) * t827 + pkin(5) * t712 + t912 * t631 + t909 * t634;
t585 = -pkin(1) * t825 + pkin(5) * t711 + t912 * t630 + t909 * t633;
t582 = -t909 * t622 + t912 * t623;
t581 = t912 * t622 + t909 * t623;
t580 = -t909 * t620 + t912 * t621;
t579 = t912 * t620 + t909 * t621;
t578 = qJ(3) * t647 - (-pkin(3) * t904 - pkin(6) * t902 - pkin(2)) * t672;
t576 = t908 * t606 + t968;
t574 = -pkin(1) * t600 - t930;
t573 = t904 * t577 + t902 * t718;
t571 = -qJ(3) * t653 + t904 * t595 - t902 * t609;
t570 = -pkin(1) * t607 - t919;
t569 = -qJ(3) * t642 + t904 * t591 - t902 * t604;
t568 = -pkin(2) * t685 + qJ(3) * t654 + t902 * t595 + t904 * t609;
t567 = -pkin(3) * t576 - t984;
t566 = -pkin(1) * t598 - t920;
t565 = -pkin(2) * t674 + qJ(3) * t643 + t902 * t591 + t904 * t604;
t564 = -pkin(6) * t626 - t908 * t584 + t911 * t590;
t563 = -pkin(5) * t600 - t909 * t578 + t912 * t592;
t562 = -pkin(6) * t576 - pkin(7) * t968 - t908 * t594;
t561 = pkin(1) * t672 + pkin(5) * t601 + t912 * t578 + t909 * t592;
t560 = -t909 * t572 + t912 * t573;
t559 = t912 * t572 + t909 * t573;
t558 = -qJ(3) * t620 + t904 * t564 - t902 * t610;
t557 = -pkin(2) * t626 + qJ(3) * t621 + t902 * t564 + t904 * t610;
t556 = -pkin(1) * t579 - t921;
t555 = -pkin(5) * t607 - t909 * t568 + t912 * t571;
t554 = -pkin(1) * t685 + pkin(5) * t608 + t912 * t568 + t909 * t571;
t553 = -pkin(5) * t598 - t909 * t565 + t912 * t569;
t552 = -pkin(1) * t674 + pkin(5) * t599 + t912 * t565 + t909 * t569;
t551 = -qJ(3) * t572 + t904 * t562 - t902 * t567;
t550 = -pkin(1) * t559 - t915;
t549 = -pkin(2) * t576 + qJ(3) * t573 + t902 * t562 + t904 * t567;
t548 = -pkin(5) * t579 - t909 * t557 + t912 * t558;
t547 = -pkin(1) * t626 + pkin(5) * t580 + t912 * t557 + t909 * t558;
t546 = -pkin(5) * t559 - t909 * t549 + t912 * t551;
t545 = -pkin(1) * t576 + pkin(5) * t560 + t912 * t549 + t909 * t551;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t996, -t995, t989, qJ(1) * t989, 0, 0, t905 * t873, 0, -t905 * t872, t890, t943 * t903 + t905 * t928, t905 * t821 + t944 * t903, t905 * t775, -qJ(1) * (t903 * t941 + t856) - (t903 * pkin(1) - t905 * pkin(5)) * t775, 0, 0, t905 * t939, 0, -t905 * t801, t890, t1000 * t905 - t903 * t945, t905 * t716 - t903 * t946, t905 * t640, t905 * t597 - t903 * t619 - qJ(1) * (t1001 * t903 + t905 * t868), t905 * t730 + t903 * t823, t905 * t704 + t903 * t799, t905 * t726 + t903 * t826, t905 * t729 - t903 * t822, t905 * t725 + t903 * t824, t905 * t766, t905 * t587 - t903 * t636 - qJ(1) * (t903 * t711 - t905 * t825), t905 * t588 - t903 * t637 - qJ(1) * (t903 * t712 - t905 * t827), t905 * t603 - t903 * t947, t905 * t563 - t903 * t574 - qJ(1) * (t903 * t601 + t672 * t905), t905 * t618 + t903 * t677, t905 * t582 + t903 * t625, t905 * t613 + t903 * t687, t905 * t617 + t903 * t676, t905 * t614 + t903 * t688, t905 * t652 + t903 * t707, t905 * t553 - t903 * t566 - qJ(1) * (t903 * t599 - t905 * t674), t905 * t555 - t903 * t570 - qJ(1) * (t903 * t608 - t905 * t685), t905 * t548 - t903 * t556 - qJ(1) * (t903 * t580 - t905 * t626), t905 * t546 - t903 * t550 - qJ(1) * (t903 * t560 - t905 * t576); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t995, -t996, t937, qJ(1) * t937, 0, 0, t903 * t873, 0, -t903 * t872, -t957, t903 * t928 - t943 * t905, t903 * t821 - t944 * t905, t903 * t775, qJ(1) * (t905 * t941 - t979) - (-t905 * pkin(1) - t903 * pkin(5)) * t775, 0, 0, t903 * t939, 0, -t903 * t801, -t957, t1000 * t903 + t905 * t945, t903 * t716 + t905 * t946, t903 * t640, t903 * t597 + t905 * t619 + qJ(1) * (t1001 * t905 - t903 * t868), t903 * t730 - t905 * t823, t903 * t704 - t905 * t799, t903 * t726 - t905 * t826, t903 * t729 + t905 * t822, t903 * t725 - t905 * t824, t903 * t766, t903 * t587 + t905 * t636 + qJ(1) * (t905 * t711 + t903 * t825), t903 * t588 + t905 * t637 + qJ(1) * (t905 * t712 + t903 * t827), t903 * t603 + t905 * t947, t903 * t563 + t905 * t574 + qJ(1) * (t905 * t601 - t672 * t903), t903 * t618 - t905 * t677, t903 * t582 - t905 * t625, t903 * t613 - t905 * t687, t903 * t617 - t905 * t676, t903 * t614 - t905 * t688, t903 * t652 - t905 * t707, t903 * t553 + t905 * t566 + qJ(1) * (t905 * t599 + t903 * t674), t903 * t555 + t905 * t570 + qJ(1) * (t905 * t608 + t903 * t685), t903 * t548 + t905 * t556 + qJ(1) * (t905 * t580 + t903 * t626), t903 * t546 + t905 * t550 + qJ(1) * (t905 * t560 + t903 * t576); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t876, t962, 0, 0, 0, 0, t872, 0, t873, 0, -t821, t928, t941, t759, 0, 0, t801, 0, t939, 0, -t716, t1000, t1001, t596, t728, t703, t724, t727, t723, t765, t585, t586, t602, t561, t616, t581, t611, t615, t612, t651, t552, t554, t547, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t964, -t876, 0, 0, 0, t873, 0, -t872, 0, t928, t821, t775, pkin(5) * t775, 0, 0, t939, 0, -t801, 0, t1000, t716, t640, t597, t730, t704, t726, t729, t725, t766, t587, t588, t603, t563, t618, t582, t613, t617, t614, t652, t553, t555, t548, t546; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t964, 0, -t962, 0, 0, 0, 0, 0, 0, -qJDD(2), t797, t796, 0, pkin(1) * t775, 0, 0, 0, 0, 0, -qJDD(2), t706, t705, 0, t619, -t823, -t799, -t826, t822, -t824, 0, t636, t637, t629, t574, -t677, -t625, -t687, -t676, -t688, -t707, t566, t570, t556, t550; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t876, t962, 0, 0, 0, 0, t872, 0, t873, 0, -t821, t928, t941, t759, 0, 0, t801, 0, t939, 0, -t716, t1000, t1001, t596, t728, t703, t724, t727, t723, t765, t585, t586, t602, t561, t616, t581, t611, t615, t612, t651, t552, t554, t547, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t986, 0, 0, -t876, t841, 0, 0, 0, t870, 0, -t869, 0, t927, t811, t694, qJ(3) * t694, t790, t769, t786, t789, t785, t836, t633, t634, t660, t592, t664, t623, t657, t663, t658, t700, t569, t571, t558, t551; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t986, 0, qJDD(2), 0, t876, 0, t842, 0, 0, 0, t869, 0, t870, 0, -t811, t927, t942, t691, t788, t768, t784, t787, t783, t835, t630, t631, t659, t578, t662, t622, t655, t661, t656, t699, t565, t568, t557, t549; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t841, -t842, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t924, t929, 0, -t985, t823, t799, t826, -t822, t824, 0, t925, t926, t918, t930, t677, t625, t687, t676, t688, t707, t920, t919, t921, t915; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t986, 0, 0, -t868, t763, 0, t838, t800, t830, t837, t828, t862, t719, t720, t672, pkin(6) * t672, t679, t627, t689, t678, t690, t708, t591, t595, t564, t562; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t986, 0, qJDD(2), 0, t868, 0, t764, 0, t884, -t875, -t956, -t884, -t955, -qJDD(4), t701, t702, 0, pkin(3) * t672, -t798, -t795, -t746, t798, t742, -t898, t604, t609, t610, t567; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t763, -t764, 0, 0, t823, t799, t826, -t822, t824, 0, t950, t951, t933, t963, t677, t625, t687, t676, t688, t707, t935, t934, t936, t922; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t865, t866, t877, -t888, t882, t888, 0, t760, t731, 0, t738, t683, t756, t736, t757, t778, t666, t671, t590, -pkin(7) * t605; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t949, t864, t880, -t917, t878, -t949, -t760, 0, t732, 0, t737, t681, t754, t735, t755, t777, t635, t645, t584, t594; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t884, t875, t956, t884, t955, qJDD(4), -t731, -t732, 0, 0, t798, t795, t746, -t798, -t742, t898, t923, t916, t983, t984; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t771, -t741, t991, t844, t839, -t844, 0, t718, t649, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t982, t990, t840, t770, t793, -t982, -t718, 0, t650, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t798, t795, t746, -t798, -t742, t898, -t649, -t650, 0, 0;];
m_new_reg = t1;