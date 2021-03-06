% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 16:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RPRPPR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:25:40
% EndTime: 2019-05-05 16:26:08
% DurationCPUTime: 27.62s
% Computational Cost: add. (181854->876), mult. (411248->1356), div. (0->0), fcn. (285793->12), ass. (0->599)
t927 = sin(qJ(1));
t1004 = qJD(1) ^ 2;
t930 = cos(qJ(1));
t891 = g(1) * t930 + g(2) * t927;
t876 = -pkin(1) * t1004 - t891;
t920 = sin(pkin(9));
t923 = cos(pkin(9));
t890 = g(1) * t927 - t930 * g(2);
t936 = qJDD(1) * pkin(1) + t890;
t819 = t920 * t876 - t923 * t936;
t820 = t923 * t876 + t920 * t936;
t945 = t819 * t920 + t923 * t820;
t747 = t819 * t923 - t820 * t920;
t986 = t747 * t930;
t667 = t927 * t945 - t986;
t987 = t747 * t927;
t668 = t930 * t945 + t987;
t880 = qJDD(1) * t920 + t1004 * t923;
t881 = qJDD(1) * t923 - t1004 * t920;
t824 = -t880 * t927 + t930 * t881;
t915 = g(3) - qJDD(2);
t853 = qJ(2) * t880 - t915 * t923;
t937 = -qJ(2) * t881 - t915 * t920;
t1024 = -pkin(6) * t824 + t853 * t927 + t930 * t937;
t919 = sin(pkin(10));
t922 = cos(pkin(10));
t929 = cos(qJ(3));
t966 = qJD(1) * t929;
t926 = sin(qJ(3));
t967 = qJD(1) * t926;
t865 = t919 * t967 - t922 * t966;
t867 = t919 * t966 + t922 * t967;
t818 = t867 * t865;
t1008 = qJDD(3) - t818;
t1023 = t1008 * t919;
t1022 = t1008 * t922;
t918 = sin(pkin(11));
t921 = cos(pkin(11));
t835 = -t921 * qJD(3) + t867 * t918;
t838 = qJD(3) * t918 + t867 * t921;
t772 = t838 * t835;
t959 = qJD(1) * qJD(3);
t948 = t929 * t959;
t958 = qJDD(1) * t926;
t878 = t948 + t958;
t908 = t929 * qJDD(1);
t949 = t926 * t959;
t939 = -t908 + t949;
t821 = t878 * t919 + t922 * t939;
t1010 = -t772 + t821;
t1021 = t1010 * t918;
t1020 = t1010 * t921;
t925 = sin(qJ(6));
t928 = cos(qJ(6));
t763 = t928 * t835 + t838 * t925;
t765 = -t835 * t925 + t838 * t928;
t699 = t765 * t763;
t938 = qJDD(6) + t821;
t1012 = -t699 + t938;
t1019 = t1012 * t925;
t1018 = t1012 * t928;
t1009 = t930 * t880 + t881 * t927;
t1017 = pkin(6) * t1009 + t853 * t930 - t927 * t937;
t860 = qJD(6) + t865;
t740 = t860 * t763;
t822 = t922 * t878 - t919 * t939;
t796 = t921 * qJDD(3) - t822 * t918;
t797 = t918 * qJDD(3) + t921 * t822;
t935 = qJD(6) * t763 - t796 * t925 - t797 * t928;
t1011 = -t740 - t935;
t794 = t865 * t835;
t728 = -t794 - t797;
t729 = -t797 + t794;
t963 = qJD(3) * t867;
t774 = t821 + t963;
t946 = -t928 * t796 + t925 * t797;
t619 = (qJD(6) - t860) * t765 + t946;
t761 = t763 ^ 2;
t762 = t765 ^ 2;
t1006 = t835 ^ 2;
t834 = t838 ^ 2;
t859 = t860 ^ 2;
t1005 = t865 ^ 2;
t864 = t867 ^ 2;
t1003 = 2 * qJD(4);
t1002 = pkin(4) * t919;
t934 = qJDD(1) * pkin(7) + t820;
t933 = -pkin(2) * t1004 + t934;
t782 = -t926 * t915 + t929 * t933;
t887 = qJD(3) * pkin(3) - qJ(4) * t967;
t914 = t929 ^ 2;
t910 = t914 * t1004;
t734 = -pkin(3) * t910 - qJ(4) * t939 - qJD(3) * t887 + t782;
t969 = t929 * t915;
t932 = -t926 * t934 - t969 - t878 * qJ(4) + qJDD(3) * pkin(3) + (qJD(3) * t929 * qJ(4) + (pkin(3) * t929 + pkin(2)) * t967) * qJD(1);
t647 = -0.2e1 * qJD(4) * t865 + t922 * t734 + t919 * t932;
t807 = pkin(4) * t865 - qJ(5) * t867;
t931 = qJD(3) ^ 2;
t605 = -pkin(4) * t931 + qJDD(3) * qJ(5) - t807 * t865 + t647;
t798 = -qJDD(1) * pkin(2) - t1004 * pkin(7) + t819;
t744 = t939 * pkin(3) - qJ(4) * t910 + t887 * t967 + qJDD(4) + t798;
t964 = qJD(3) * t865;
t943 = -t822 + t964;
t650 = pkin(4) * t774 + t943 * qJ(5) + t744;
t554 = 0.2e1 * qJD(5) * t838 + t918 * t605 - t921 * t650;
t515 = pkin(5) * t1010 + pkin(8) * t728 - t554;
t555 = -0.2e1 * qJD(5) * t835 + t921 * t605 + t918 * t650;
t783 = pkin(5) * t865 - pkin(8) * t838;
t529 = -pkin(5) * t1006 + pkin(8) * t796 - t783 * t865 + t555;
t464 = -t928 * t515 + t529 * t925;
t465 = t925 * t515 + t928 * t529;
t412 = -t464 * t928 + t465 * t925;
t1001 = t412 * t918;
t1000 = t412 * t921;
t947 = t919 * t734 - t922 * t932;
t604 = qJDD(5) - t931 * qJ(5) - qJDD(3) * pkin(4) + (t1003 + t807) * t867 + t947;
t572 = -t796 * pkin(5) - pkin(8) * t1006 + t838 * t783 + t604;
t999 = t572 * t925;
t998 = t572 * t928;
t646 = t1003 * t867 + t947;
t576 = -t646 * t922 + t647 * t919;
t997 = t576 * t926;
t996 = t576 * t929;
t995 = t604 * t918;
t994 = t604 * t921;
t674 = t699 + t938;
t993 = t674 * t925;
t992 = t674 * t928;
t731 = t772 + t821;
t991 = t731 * t918;
t990 = t731 * t921;
t989 = t744 * t919;
t988 = t744 * t922;
t985 = t798 * t926;
t984 = t798 * t929;
t811 = qJDD(3) + t818;
t983 = t811 * t919;
t982 = t811 * t922;
t981 = t821 * t919;
t980 = t838 * t865;
t979 = t860 * t925;
t978 = t860 * t928;
t977 = t865 * t918;
t976 = t865 * t921;
t879 = t908 - 0.2e1 * t949;
t975 = t879 * t929;
t899 = t929 * t1004 * t926;
t888 = qJDD(3) + t899;
t972 = t888 * t926;
t889 = qJDD(3) - t899;
t971 = t889 * t926;
t970 = t889 * t929;
t913 = t926 ^ 2;
t968 = t913 + t914;
t965 = t1004 * t913;
t962 = qJD(3) * t919;
t961 = qJD(3) * t922;
t957 = qJDD(3) * t923;
t956 = t919 * t699;
t955 = t922 * t699;
t954 = t919 * t772;
t953 = t922 * t772;
t952 = t920 * t818;
t951 = t923 * t818;
t950 = -pkin(4) * t922 - pkin(3);
t413 = t464 * t925 + t928 * t465;
t577 = t646 * t919 + t922 * t647;
t780 = t926 * t933 + t969;
t706 = t780 * t926 + t929 * t782;
t839 = -t890 * t927 - t930 * t891;
t942 = t920 * t899;
t941 = t923 * t899;
t884 = qJDD(1) * t930 - t1004 * t927;
t940 = -pkin(6) * t884 - g(3) * t927;
t482 = -t554 * t921 + t555 * t918;
t483 = t554 * t918 + t555 * t921;
t705 = t780 * t929 - t782 * t926;
t837 = t890 * t930 - t891 * t927;
t724 = -t796 - t980;
t776 = -t821 + t963;
t904 = t920 * qJDD(3);
t898 = -t910 - t931;
t897 = t910 - t931;
t896 = -t931 - t965;
t895 = t931 - t965;
t886 = t910 - t965;
t885 = t910 + t965;
t883 = qJDD(1) * t927 + t1004 * t930;
t882 = t968 * qJDD(1);
t877 = 0.2e1 * t948 + t958;
t874 = t929 * t888;
t873 = t968 * t959;
t861 = -pkin(6) * t883 + g(3) * t930;
t856 = -t864 - t931;
t855 = -t864 + t931;
t854 = t1005 - t931;
t849 = t878 * t929 - t913 * t959;
t848 = -t914 * t959 + t926 * t939;
t847 = t873 * t923 + t904;
t846 = t873 * t920 - t957;
t845 = -t896 * t926 - t970;
t844 = -t895 * t926 + t874;
t843 = t898 * t929 - t972;
t842 = t897 * t929 - t971;
t841 = t896 * t929 - t971;
t840 = t898 * t926 + t874;
t830 = t882 * t923 - t885 * t920;
t829 = t882 * t920 + t885 * t923;
t823 = -t877 * t926 + t975;
t815 = -t864 + t1005;
t814 = t922 * t821;
t809 = -t1005 - t931;
t806 = t849 * t923 - t942;
t805 = t848 * t923 + t942;
t804 = t849 * t920 + t941;
t803 = t848 * t920 - t941;
t802 = t844 * t923 + t920 * t958;
t801 = t842 * t923 + t908 * t920;
t800 = t844 * t920 - t923 * t958;
t799 = t842 * t920 - t908 * t923;
t792 = (-t865 * t922 + t867 * t919) * qJD(3);
t791 = (-t865 * t919 - t867 * t922) * qJD(3);
t790 = -t834 + t1005;
t789 = -t1005 + t1006;
t788 = t845 * t923 + t877 * t920;
t787 = t843 * t923 - t879 * t920;
t786 = t845 * t920 - t877 * t923;
t785 = t843 * t920 + t879 * t923;
t781 = t823 * t923 - t886 * t920;
t779 = t823 * t920 + t886 * t923;
t778 = -t822 - t964;
t773 = -t1005 - t864;
t770 = -t834 + t1006;
t769 = t822 * t922 - t867 * t962;
t768 = t822 * t919 + t867 * t961;
t767 = t865 * t961 + t981;
t766 = t865 * t962 - t814;
t760 = -t834 - t1005;
t759 = -t856 * t919 - t982;
t758 = -t855 * t919 + t1022;
t757 = t854 * t922 - t983;
t756 = t856 * t922 - t983;
t755 = t855 * t922 + t1023;
t754 = t854 * t919 + t982;
t753 = -t829 * t927 + t830 * t930;
t752 = t829 * t930 + t830 * t927;
t751 = -pkin(7) * t841 + t984;
t750 = -pkin(7) * t840 + t985;
t749 = -t1005 - t1006;
t743 = -pkin(2) * t841 + t782;
t742 = -pkin(2) * t840 + t780;
t739 = -t762 + t859;
t738 = t761 - t859;
t737 = t834 + t1006;
t736 = t809 * t922 - t1023;
t735 = t809 * t919 + t1022;
t733 = pkin(1) * t915 + qJ(2) * t945;
t725 = t796 - t980;
t720 = t797 * t921 - t838 * t977;
t719 = -t797 * t918 - t838 * t976;
t718 = -t796 * t918 + t835 * t976;
t717 = -t796 * t921 - t835 * t977;
t716 = (-t835 * t921 + t838 * t918) * t865;
t715 = (t835 * t918 + t838 * t921) * t865;
t714 = -t791 * t926 + t792 * t929;
t713 = -t762 - t859;
t712 = -t786 * t927 + t788 * t930;
t711 = -t785 * t927 + t787 * t930;
t710 = t786 * t930 + t788 * t927;
t709 = t785 * t930 + t787 * t927;
t708 = t714 * t923 + t904;
t707 = t714 * t920 - t957;
t703 = t776 * t922 - t778 * t919;
t702 = -t774 * t922 + t919 * t943;
t701 = t776 * t919 + t778 * t922;
t700 = -t774 * t919 - t922 * t943;
t698 = -t762 + t761;
t697 = -t768 * t926 + t769 * t929;
t696 = -t766 * t926 + t767 * t929;
t695 = t716 * t922 + t981;
t694 = t716 * t919 - t814;
t693 = -t756 * t926 + t759 * t929;
t692 = -t755 * t926 + t758 * t929;
t691 = -t754 * t926 + t757 * t929;
t690 = t756 * t929 + t759 * t926;
t689 = t789 * t921 - t991;
t688 = -t790 * t918 + t1020;
t687 = -t789 * t918 - t990;
t686 = -t790 * t921 - t1021;
t685 = -t859 - t761;
t684 = -qJ(2) * t829 + t705 * t923;
t683 = qJ(2) * t830 + t705 * t920;
t682 = -qJ(4) * t756 + t988;
t681 = t720 * t922 + t954;
t680 = t718 * t922 - t954;
t679 = t720 * t919 - t953;
t678 = t718 * t919 + t953;
t676 = -qJD(6) * t765 - t946;
t672 = -t760 * t918 - t990;
t671 = t760 * t921 - t991;
t670 = t706 * t923 + t798 * t920;
t669 = t706 * t920 - t798 * t923;
t666 = -qJ(4) * t735 + t989;
t665 = (-t763 * t928 + t765 * t925) * t860;
t664 = (-t763 * t925 - t765 * t928) * t860;
t663 = t749 * t921 - t1021;
t662 = t749 * t918 + t1020;
t661 = t697 * t923 + t952;
t660 = t696 * t923 - t952;
t659 = t697 * t920 - t951;
t658 = t696 * t920 + t951;
t657 = -t735 * t926 + t736 * t929;
t656 = t735 * t929 + t736 * t926;
t655 = -t724 * t921 - t728 * t918;
t654 = t725 * t921 + t729 * t918;
t653 = -t724 * t918 + t728 * t921;
t652 = -t725 * t918 + t729 * t921;
t651 = -t761 - t762;
t645 = t693 * t923 - t920 * t943;
t644 = t692 * t923 - t778 * t920;
t643 = t691 * t923 + t776 * t920;
t642 = t693 * t920 + t923 * t943;
t641 = t692 * t920 + t778 * t923;
t640 = t691 * t920 - t776 * t923;
t638 = pkin(3) * t943 + qJ(4) * t759 + t989;
t637 = -qJ(2) * t786 - t743 * t920 + t751 * t923;
t636 = -qJ(2) * t785 - t742 * t920 + t750 * t923;
t635 = -pkin(3) * t774 + qJ(4) * t736 - t988;
t634 = -t701 * t926 + t703 * t929;
t633 = -t700 * t926 + t702 * t929;
t632 = t701 * t929 + t703 * t926;
t631 = t657 * t923 + t774 * t920;
t630 = t657 * t920 - t774 * t923;
t629 = -pkin(1) * t841 + qJ(2) * t788 + t743 * t923 + t751 * t920;
t628 = -pkin(1) * t840 + qJ(2) * t787 + t742 * t923 + t750 * t920;
t627 = t689 * t922 - t724 * t919;
t626 = t688 * t922 - t728 * t919;
t625 = t689 * t919 + t724 * t922;
t624 = t688 * t919 + t728 * t922;
t623 = -t740 + t935;
t618 = (qJD(6) + t860) * t765 + t946;
t617 = -t765 * t979 - t928 * t935;
t616 = t765 * t978 - t925 * t935;
t615 = -t676 * t925 + t763 * t978;
t614 = t676 * t928 + t763 * t979;
t613 = t738 * t928 - t993;
t612 = -t739 * t925 + t1018;
t611 = t738 * t925 + t992;
t610 = t739 * t928 + t1019;
t609 = t654 * t922 - t770 * t919;
t608 = t654 * t919 + t770 * t922;
t607 = t672 * t922 - t729 * t919;
t606 = t672 * t919 + t729 * t922;
t602 = t663 * t922 - t725 * t919;
t601 = t663 * t919 + t725 * t922;
t600 = -t713 * t925 - t992;
t599 = t713 * t928 - t993;
t598 = t655 * t922 - t737 * t919;
t597 = t655 * t919 + t737 * t922;
t596 = t633 * t923 - t815 * t920;
t595 = t633 * t920 + t815 * t923;
t594 = -t694 * t926 + t695 * t929;
t593 = t634 * t923 + t773 * t920;
t592 = t634 * t920 - t773 * t923;
t591 = t685 * t928 - t1019;
t590 = t685 * t925 + t1018;
t589 = -t679 * t926 + t681 * t929;
t588 = -t678 * t926 + t680 * t929;
t587 = -t669 * t927 + t670 * t930;
t586 = t669 * t930 + t670 * t927;
t585 = -t664 * t918 + t665 * t921;
t584 = -t664 * t921 - t665 * t918;
t583 = -pkin(2) * t632 - pkin(3) * t701;
t582 = t585 * t922 + t919 * t938;
t581 = t585 * t919 - t922 * t938;
t580 = t594 * t923 - t715 * t920;
t579 = t594 * t920 + t715 * t923;
t578 = -qJ(2) * t669 - (pkin(2) * t920 - pkin(7) * t923) * t705;
t575 = -pkin(2) * t690 - pkin(3) * t756 + t647;
t574 = -t642 * t927 + t645 * t930;
t573 = t642 * t930 + t645 * t927;
t571 = t589 * t923 - t719 * t920;
t570 = t588 * t923 - t717 * t920;
t569 = t589 * t920 + t719 * t923;
t568 = t588 * t920 + t717 * t923;
t567 = -qJ(5) * t671 + t994;
t566 = -qJ(5) * t662 + t995;
t565 = -pkin(2) * t656 - pkin(3) * t735 + t646;
t564 = -t630 * t927 + t631 * t930;
t563 = t630 * t930 + t631 * t927;
t562 = -t625 * t926 + t627 * t929;
t561 = -t624 * t926 + t626 * t929;
t560 = -t619 * t928 - t623 * t925;
t559 = -t1011 * t925 - t618 * t928;
t558 = -t619 * t925 + t623 * t928;
t557 = t1011 * t928 - t618 * t925;
t556 = -pkin(3) * t744 + qJ(4) * t577;
t553 = -t616 * t918 + t617 * t921;
t552 = -t614 * t918 + t615 * t921;
t551 = -t616 * t921 - t617 * t918;
t550 = -t614 * t921 - t615 * t918;
t549 = -pkin(7) * t690 - t638 * t926 + t682 * t929;
t548 = -t608 * t926 + t609 * t929;
t547 = -t611 * t918 + t613 * t921;
t546 = -t610 * t918 + t612 * t921;
t545 = -t611 * t921 - t613 * t918;
t544 = -t610 * t921 - t612 * t918;
t543 = qJ(2) * t670 - (-pkin(2) * t923 - pkin(7) * t920 - pkin(1)) * t705;
t542 = -qJ(4) * t701 - t576;
t541 = -t606 * t926 + t607 * t929;
t540 = t606 * t929 + t607 * t926;
t539 = -t601 * t926 + t602 * t929;
t538 = t601 * t929 + t602 * t926;
t537 = -t599 * t918 + t600 * t921;
t536 = t599 * t921 + t600 * t918;
t535 = -t597 * t926 + t598 * t929;
t534 = t597 * t929 + t598 * t926;
t533 = -pkin(7) * t656 - t635 * t926 + t666 * t929;
t532 = -pkin(3) * t773 + qJ(4) * t703 + t577;
t531 = -t592 * t927 + t593 * t930;
t530 = t592 * t930 + t593 * t927;
t527 = -t590 * t918 + t591 * t921;
t526 = t590 * t921 + t591 * t918;
t525 = t553 * t922 + t956;
t524 = t552 * t922 - t956;
t523 = t553 * t919 - t955;
t522 = t552 * t919 + t955;
t521 = t562 * t923 - t687 * t920;
t520 = t561 * t923 - t686 * t920;
t519 = t562 * t920 + t687 * t923;
t518 = t561 * t920 + t686 * t923;
t517 = -pkin(4) * t671 + t555;
t516 = -pkin(4) * t662 + t554;
t512 = -pkin(8) * t599 + t998;
t511 = t541 * t923 + t671 * t920;
t510 = t541 * t920 - t671 * t923;
t509 = t548 * t923 - t652 * t920;
t508 = t548 * t920 + t652 * t923;
t507 = t539 * t923 + t662 * t920;
t506 = t539 * t920 - t662 * t923;
t505 = -pkin(8) * t590 + t999;
t504 = t535 * t923 + t653 * t920;
t503 = t535 * t920 - t653 * t923;
t502 = t547 * t922 - t619 * t919;
t501 = t546 * t922 - t623 * t919;
t500 = t547 * t919 + t619 * t922;
t499 = t546 * t919 + t623 * t922;
t498 = -t581 * t926 + t582 * t929;
t497 = t1011 * t919 + t537 * t922;
t496 = -t1011 * t922 + t537 * t919;
t495 = t577 * t929 - t997;
t494 = t577 * t926 + t996;
t493 = t527 * t922 + t618 * t919;
t492 = t527 * t919 - t618 * t922;
t491 = t495 * t923 + t744 * t920;
t490 = t495 * t920 - t744 * t923;
t489 = -pkin(5) * t1011 + pkin(8) * t600 + t999;
t488 = -pkin(5) * t618 + pkin(8) * t591 - t998;
t487 = -t558 * t918 + t560 * t921;
t486 = -t557 * t918 + t559 * t921;
t485 = t558 * t921 + t560 * t918;
t484 = -t557 * t921 - t559 * t918;
t481 = t486 * t922 - t698 * t919;
t480 = t486 * t919 + t698 * t922;
t479 = -qJ(2) * t642 + t549 * t923 - t575 * t920;
t478 = t498 * t923 - t584 * t920;
t477 = t498 * t920 + t584 * t923;
t476 = -qJ(5) * t653 - t482;
t475 = t487 * t922 + t651 * t919;
t474 = t487 * t919 - t651 * t922;
t473 = -pkin(2) * t494 - pkin(3) * t576;
t472 = t483 * t922 + t604 * t919;
t471 = t483 * t919 - t604 * t922;
t470 = -pkin(1) * t690 + qJ(2) * t645 + t549 * t920 + t575 * t923;
t469 = -t523 * t926 + t525 * t929;
t468 = -t522 * t926 + t524 * t929;
t467 = -qJ(2) * t630 + t533 * t923 - t565 * t920;
t466 = -pkin(2) * t540 - pkin(3) * t606 - pkin(4) * t729 - qJ(5) * t672 - t995;
t462 = -pkin(2) * t538 - pkin(3) * t601 - pkin(4) * t725 - qJ(5) * t663 + t994;
t461 = -pkin(7) * t632 - t532 * t926 + t542 * t929;
t460 = -qJ(4) * t606 - t517 * t919 + t567 * t922;
t459 = -qJ(4) * t601 - t516 * t919 + t566 * t922;
t458 = -t510 * t927 + t511 * t930;
t457 = t510 * t930 + t511 * t927;
t456 = -t506 * t927 + t507 * t930;
t455 = t506 * t930 + t507 * t927;
t454 = -pkin(1) * t656 + qJ(2) * t631 + t533 * t920 + t565 * t923;
t453 = -pkin(4) * t485 - pkin(5) * t558;
t452 = -pkin(3) * t671 + qJ(4) * t607 + t517 * t922 + t567 * t919;
t451 = -t503 * t927 + t504 * t930;
t450 = t503 * t930 + t504 * t927;
t449 = -t500 * t926 + t502 * t929;
t448 = -t499 * t926 + t501 * t929;
t447 = -pkin(3) * t662 + qJ(4) * t602 + t516 * t922 + t566 * t919;
t446 = -t496 * t926 + t497 * t929;
t445 = t496 * t929 + t497 * t926;
t444 = -qJ(4) * t597 + t1002 * t653 + t476 * t922;
t443 = -pkin(7) * t494 - qJ(4) * t996 - t556 * t926;
t442 = -t492 * t926 + t493 * t929;
t441 = t492 * t929 + t493 * t926;
t440 = t469 * t923 - t551 * t920;
t439 = t468 * t923 - t550 * t920;
t438 = t469 * t920 + t551 * t923;
t437 = t468 * t920 + t550 * t923;
t436 = -t490 * t927 + t491 * t930;
t435 = t490 * t930 + t491 * t927;
t434 = qJ(4) * t598 + t919 * t476 + t653 * t950;
t433 = -qJ(2) * t592 + t461 * t923 - t583 * t920;
t432 = t449 * t923 - t545 * t920;
t431 = t448 * t923 - t544 * t920;
t430 = t449 * t920 + t545 * t923;
t429 = t448 * t920 + t544 * t923;
t428 = -pkin(2) * t534 - pkin(3) * t597 - pkin(4) * t737 - qJ(5) * t655 - t483;
t427 = -pkin(4) * t536 - pkin(5) * t599 + t465;
t426 = t446 * t923 + t536 * t920;
t425 = t446 * t920 - t536 * t923;
t424 = -qJ(5) * t536 - t489 * t918 + t512 * t921;
t423 = -pkin(1) * t632 + qJ(2) * t593 + t461 * t920 + t583 * t923;
t422 = -t480 * t926 + t481 * t929;
t421 = -pkin(4) * t526 - pkin(5) * t590 + t464;
t420 = -qJ(5) * t526 - t488 * t918 + t505 * t921;
t419 = -t474 * t926 + t475 * t929;
t418 = t474 * t929 + t475 * t926;
t417 = t442 * t923 + t526 * t920;
t416 = t442 * t920 - t526 * t923;
t415 = -t471 * t926 + t472 * t929;
t414 = t471 * t929 + t472 * t926;
t411 = -pkin(5) * t572 + pkin(8) * t413;
t410 = -pkin(8) * t558 - t412;
t409 = t422 * t923 - t484 * t920;
t408 = t422 * t920 + t484 * t923;
t407 = -qJ(4) * t471 + (-qJ(5) * t922 + t1002) * t482;
t406 = -pkin(5) * t651 + pkin(8) * t560 + t413;
t405 = t419 * t923 + t485 * t920;
t404 = t419 * t920 - t485 * t923;
t403 = -pkin(7) * t540 - t452 * t926 + t460 * t929;
t402 = -pkin(7) * t538 - t447 * t926 + t459 * t929;
t401 = t415 * t923 + t482 * t920;
t400 = t415 * t920 - t482 * t923;
t399 = -qJ(2) * t490 + t443 * t923 - t473 * t920;
t398 = -pkin(7) * t534 - t434 * t926 + t444 * t929;
t397 = qJ(4) * t472 + (-qJ(5) * t919 + t950) * t482;
t396 = -t425 * t927 + t426 * t930;
t395 = t425 * t930 + t426 * t927;
t394 = -pkin(1) * t494 + qJ(2) * t491 + t443 * t920 + t473 * t923;
t393 = -pkin(2) * t445 - pkin(3) * t496 + pkin(4) * t1011 - qJ(5) * t537 - t489 * t921 - t512 * t918;
t392 = -t416 * t927 + t417 * t930;
t391 = t416 * t930 + t417 * t927;
t390 = -pkin(2) * t414 - pkin(3) * t471 + pkin(4) * t604 - qJ(5) * t483;
t389 = -qJ(4) * t496 + t424 * t922 - t427 * t919;
t388 = -pkin(2) * t441 - pkin(3) * t492 + pkin(4) * t618 - qJ(5) * t527 - t488 * t921 - t505 * t918;
t387 = t413 * t921 - t1001;
t386 = t413 * t918 + t1000;
t385 = -qJ(2) * t510 + t403 * t923 - t466 * t920;
t384 = -qJ(4) * t492 + t420 * t922 - t421 * t919;
t383 = -qJ(2) * t506 + t402 * t923 - t462 * t920;
t382 = -pkin(3) * t536 + qJ(4) * t497 + t424 * t919 + t427 * t922;
t381 = t387 * t922 + t572 * t919;
t380 = t387 * t919 - t572 * t922;
t379 = -pkin(1) * t540 + qJ(2) * t511 + t403 * t920 + t466 * t923;
t378 = -pkin(1) * t538 + qJ(2) * t507 + t402 * t920 + t462 * t923;
t377 = -pkin(3) * t526 + qJ(4) * t493 + t420 * t919 + t421 * t922;
t376 = -t404 * t927 + t405 * t930;
t375 = t404 * t930 + t405 * t927;
t374 = -qJ(2) * t503 + t398 * t923 - t428 * t920;
t373 = -t400 * t927 + t401 * t930;
t372 = t400 * t930 + t401 * t927;
t371 = -qJ(5) * t485 - t406 * t918 + t410 * t921;
t370 = -pkin(1) * t534 + qJ(2) * t504 + t398 * t920 + t428 * t923;
t369 = -pkin(4) * t386 - pkin(5) * t412;
t368 = -pkin(7) * t414 - t397 * t926 + t407 * t929;
t367 = -qJ(4) * t474 + t371 * t922 - t453 * t919;
t366 = -pkin(8) * t1000 - qJ(5) * t386 - t411 * t918;
t365 = -pkin(3) * t485 + qJ(4) * t475 + t371 * t919 + t453 * t922;
t364 = -pkin(2) * t418 - pkin(3) * t474 + pkin(4) * t651 - qJ(5) * t487 - t406 * t921 - t410 * t918;
t363 = -t380 * t926 + t381 * t929;
t362 = t380 * t929 + t381 * t926;
t361 = -pkin(7) * t445 - t382 * t926 + t389 * t929;
t360 = -pkin(7) * t441 - t377 * t926 + t384 * t929;
t359 = t363 * t923 + t386 * t920;
t358 = t363 * t920 - t386 * t923;
t357 = -qJ(2) * t400 + t368 * t923 - t390 * t920;
t356 = -qJ(2) * t425 + t361 * t923 - t393 * t920;
t355 = -pkin(1) * t445 + qJ(2) * t426 + t361 * t920 + t393 * t923;
t354 = -qJ(2) * t416 + t360 * t923 - t388 * t920;
t353 = -pkin(1) * t414 + qJ(2) * t401 + t368 * t920 + t390 * t923;
t352 = -pkin(1) * t441 + qJ(2) * t417 + t360 * t920 + t388 * t923;
t351 = -pkin(7) * t418 - t365 * t926 + t367 * t929;
t350 = -qJ(4) * t380 + t366 * t922 - t369 * t919;
t349 = -pkin(2) * t362 - pkin(3) * t380 + pkin(4) * t572 + pkin(8) * t1001 - qJ(5) * t387 - t411 * t921;
t348 = -pkin(3) * t386 + qJ(4) * t381 + t366 * t919 + t369 * t922;
t347 = -t358 * t927 + t359 * t930;
t346 = t358 * t930 + t359 * t927;
t345 = -qJ(2) * t404 + t351 * t923 - t364 * t920;
t344 = -pkin(1) * t418 + qJ(2) * t405 + t351 * t920 + t364 * t923;
t343 = -pkin(7) * t362 - t348 * t926 + t350 * t929;
t342 = -qJ(2) * t358 + t343 * t923 - t349 * t920;
t341 = -pkin(1) * t362 + qJ(2) * t359 + t343 * t920 + t349 * t923;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t883, -t884, 0, t839, 0, 0, 0, 0, 0, 0, -t1009, -t824, 0, t668, 0, 0, 0, 0, 0, 0, t711, t712, t753, t587, 0, 0, 0, 0, 0, 0, t564, t574, t531, t436, 0, 0, 0, 0, 0, 0, t456, t458, t451, t373, 0, 0, 0, 0, 0, 0, t392, t396, t376, t347; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t884, -t883, 0, t837, 0, 0, 0, 0, 0, 0, t824, -t1009, 0, t667, 0, 0, 0, 0, 0, 0, t709, t710, t752, t586, 0, 0, 0, 0, 0, 0, t563, t573, t530, t435, 0, 0, 0, 0, 0, 0, t455, t457, t450, t372, 0, 0, 0, 0, 0, 0, t391, t395, t375, t346; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t915, 0, 0, 0, 0, 0, 0, t840, t841, 0, -t705, 0, 0, 0, 0, 0, 0, t656, t690, t632, t494, 0, 0, 0, 0, 0, 0, t538, t540, t534, t414, 0, 0, 0, 0, 0, 0, t441, t445, t418, t362; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t884, 0, -t883, 0, t940, -t861, -t837, -pkin(6) * t837, 0, 0, t824, 0, -t1009, 0, t1024, t1017, -t667, -pkin(6) * t667 + qJ(2) * t986 - t733 * t927, -t804 * t927 + t806 * t930, -t779 * t927 + t781 * t930, -t800 * t927 + t802 * t930, -t803 * t927 + t805 * t930, -t799 * t927 + t801 * t930, -t846 * t927 + t847 * t930, -pkin(6) * t709 - t628 * t927 + t636 * t930, -pkin(6) * t710 - t629 * t927 + t637 * t930, -pkin(6) * t752 - t683 * t927 + t684 * t930, -pkin(6) * t586 - t543 * t927 + t578 * t930, -t659 * t927 + t661 * t930, -t595 * t927 + t596 * t930, -t641 * t927 + t644 * t930, -t658 * t927 + t660 * t930, -t640 * t927 + t643 * t930, -t707 * t927 + t708 * t930, -pkin(6) * t563 - t454 * t927 + t467 * t930, -pkin(6) * t573 - t470 * t927 + t479 * t930, -pkin(6) * t530 - t423 * t927 + t433 * t930, -pkin(6) * t435 - t394 * t927 + t399 * t930, -t569 * t927 + t571 * t930, -t508 * t927 + t509 * t930, -t518 * t927 + t520 * t930, -t568 * t927 + t570 * t930, -t519 * t927 + t521 * t930, -t579 * t927 + t580 * t930, -pkin(6) * t455 - t378 * t927 + t383 * t930, -pkin(6) * t457 - t379 * t927 + t385 * t930, -pkin(6) * t450 - t370 * t927 + t374 * t930, -pkin(6) * t372 - t353 * t927 + t357 * t930, -t438 * t927 + t440 * t930, -t408 * t927 + t409 * t930, -t429 * t927 + t431 * t930, -t437 * t927 + t439 * t930, -t430 * t927 + t432 * t930, -t477 * t927 + t478 * t930, -pkin(6) * t391 - t352 * t927 + t354 * t930, -pkin(6) * t395 - t355 * t927 + t356 * t930, -pkin(6) * t375 - t344 * t927 + t345 * t930, -pkin(6) * t346 - t341 * t927 + t342 * t930; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t883, 0, t884, 0, t861, t940, t839, pkin(6) * t839, 0, 0, t1009, 0, t824, 0, -t1017, t1024, t668, pkin(6) * t668 + qJ(2) * t987 + t733 * t930, t804 * t930 + t806 * t927, t779 * t930 + t781 * t927, t800 * t930 + t802 * t927, t803 * t930 + t805 * t927, t799 * t930 + t801 * t927, t846 * t930 + t847 * t927, pkin(6) * t711 + t628 * t930 + t636 * t927, pkin(6) * t712 + t629 * t930 + t637 * t927, pkin(6) * t753 + t683 * t930 + t684 * t927, pkin(6) * t587 + t543 * t930 + t578 * t927, t659 * t930 + t661 * t927, t595 * t930 + t596 * t927, t641 * t930 + t644 * t927, t658 * t930 + t660 * t927, t640 * t930 + t643 * t927, t707 * t930 + t708 * t927, pkin(6) * t564 + t454 * t930 + t467 * t927, pkin(6) * t574 + t470 * t930 + t479 * t927, pkin(6) * t531 + t423 * t930 + t433 * t927, pkin(6) * t436 + t394 * t930 + t399 * t927, t569 * t930 + t571 * t927, t508 * t930 + t509 * t927, t518 * t930 + t520 * t927, t568 * t930 + t570 * t927, t519 * t930 + t521 * t927, t579 * t930 + t580 * t927, pkin(6) * t456 + t378 * t930 + t383 * t927, pkin(6) * t458 + t379 * t930 + t385 * t927, pkin(6) * t451 + t370 * t930 + t374 * t927, pkin(6) * t373 + t353 * t930 + t357 * t927, t438 * t930 + t440 * t927, t408 * t930 + t409 * t927, t429 * t930 + t431 * t927, t437 * t930 + t439 * t927, t430 * t930 + t432 * t927, t477 * t930 + t478 * t927, pkin(6) * t392 + t352 * t930 + t354 * t927, pkin(6) * t396 + t355 * t930 + t356 * t927, pkin(6) * t376 + t344 * t930 + t345 * t927, pkin(6) * t347 + t341 * t930 + t342 * t927; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t890, t891, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t881 - t819, -pkin(1) * t880 - t820, 0, -pkin(1) * t747, (t878 + t948) * t926, t877 * t929 + t879 * t926, t895 * t929 + t972, t975, t897 * t926 + t970, 0, pkin(1) * t785 + pkin(2) * t879 + pkin(7) * t843 - t984, pkin(1) * t786 - pkin(2) * t877 + pkin(7) * t845 + t985, pkin(1) * t829 + pkin(2) * t885 + pkin(7) * t882 + t706, pkin(1) * t669 - pkin(2) * t798 + pkin(7) * t706, t768 * t929 + t769 * t926, t700 * t929 + t702 * t926, t755 * t929 + t758 * t926, t766 * t929 + t767 * t926, t754 * t929 + t757 * t926, t791 * t929 + t792 * t926, pkin(1) * t630 - pkin(2) * t774 + pkin(7) * t657 + t635 * t929 + t666 * t926, pkin(1) * t642 + pkin(2) * t943 + pkin(7) * t693 + t638 * t929 + t682 * t926, pkin(1) * t592 - pkin(2) * t773 + pkin(7) * t634 + t532 * t929 + t542 * t926, pkin(1) * t490 - pkin(2) * t744 + pkin(7) * t495 - qJ(4) * t997 + t556 * t929, t679 * t929 + t681 * t926, t608 * t929 + t609 * t926, t624 * t929 + t626 * t926, t678 * t929 + t680 * t926, t625 * t929 + t627 * t926, t694 * t929 + t695 * t926, pkin(1) * t506 - pkin(2) * t662 + pkin(7) * t539 + t447 * t929 + t459 * t926, pkin(1) * t510 - pkin(2) * t671 + pkin(7) * t541 + t452 * t929 + t460 * t926, pkin(1) * t503 - pkin(2) * t653 + pkin(7) * t535 + t434 * t929 + t444 * t926, pkin(1) * t400 - pkin(2) * t482 + pkin(7) * t415 + t397 * t929 + t407 * t926, t523 * t929 + t525 * t926, t480 * t929 + t481 * t926, t499 * t929 + t501 * t926, t522 * t929 + t524 * t926, t500 * t929 + t502 * t926, t581 * t929 + t582 * t926, pkin(1) * t416 - pkin(2) * t526 + pkin(7) * t442 + t377 * t929 + t384 * t926, pkin(1) * t425 - pkin(2) * t536 + pkin(7) * t446 + t382 * t929 + t389 * t926, pkin(1) * t404 - pkin(2) * t485 + pkin(7) * t419 + t365 * t929 + t367 * t926, pkin(1) * t358 - pkin(2) * t386 + pkin(7) * t363 + t348 * t929 + t350 * t926;];
tauB_reg  = t1;
