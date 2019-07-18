% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 16:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRPRPR11_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:03:47
% EndTime: 2019-05-06 16:04:22
% DurationCPUTime: 28.02s
% Computational Cost: add. (156732->777), mult. (334711->1126), div. (0->0), fcn. (219067->10), ass. (0->537)
t908 = sin(qJ(2));
t900 = t908 ^ 2;
t915 = qJD(1) ^ 2;
t893 = t900 * t915;
t914 = qJD(2) ^ 2;
t879 = -t893 - t914;
t912 = cos(qJ(2));
t945 = t908 * t912 * t915;
t872 = -qJDD(2) + t945;
t967 = t912 * t872;
t817 = -t879 * t908 + t967;
t957 = qJD(1) * qJD(2);
t889 = t912 * t957;
t955 = qJDD(1) * t908;
t861 = 0.2e1 * t889 + t955;
t909 = sin(qJ(1));
t913 = cos(qJ(1));
t766 = t817 * t909 - t861 * t913;
t1043 = pkin(6) * t766;
t770 = t817 * t913 + t861 * t909;
t1042 = pkin(6) * t770;
t901 = t912 ^ 2;
t974 = t901 * t915;
t880 = t914 + t974;
t871 = qJDD(2) + t945;
t986 = t871 * t908;
t819 = t880 * t912 + t986;
t892 = t912 * qJDD(1);
t939 = t908 * t957;
t864 = t892 - 0.2e1 * t939;
t767 = t819 * t909 - t864 * t913;
t1041 = pkin(6) * t767;
t771 = t819 * t913 + t864 * t909;
t1040 = pkin(6) * t771;
t968 = t912 * t871;
t814 = -t880 * t908 + t968;
t1039 = pkin(1) * t814;
t1038 = pkin(7) * t814;
t1037 = pkin(7) * t817;
t878 = -t893 + t914;
t816 = -t878 * t908 + t968;
t953 = qJDD(1) * t913;
t1036 = t816 * t909 - t908 * t953;
t954 = qJDD(1) * t909;
t1035 = t816 * t913 + t908 * t954;
t1034 = 2 * qJD(3);
t1033 = -2 * qJD(5);
t985 = t872 * t908;
t811 = t879 * t912 + t985;
t1032 = pkin(1) * t811;
t1031 = pkin(7) * t811;
t1030 = pkin(7) * t819;
t907 = sin(qJ(4));
t911 = cos(qJ(4));
t961 = qJD(1) * t912;
t854 = qJD(2) * t907 + t911 * t961;
t856 = qJD(2) * t911 - t907 * t961;
t807 = t856 * t854;
t862 = t889 + t955;
t847 = qJDD(4) + t862;
t1017 = -t807 + t847;
t1029 = t1017 * t907;
t1028 = t1017 * t911;
t903 = sin(pkin(10));
t904 = cos(pkin(10));
t796 = t904 * t854 + t856 * t903;
t798 = -t854 * t903 + t856 * t904;
t735 = t798 * t796;
t1018 = -t735 + t847;
t1027 = t1018 * t903;
t1026 = t1018 * t904;
t906 = sin(qJ(6));
t910 = cos(qJ(6));
t727 = t910 * t796 + t798 * t906;
t729 = -t796 * t906 + t798 * t910;
t659 = t729 * t727;
t843 = qJDD(6) + t847;
t1020 = -t659 + t843;
t1025 = t1020 * t906;
t1024 = t1020 * t910;
t881 = -t914 + t974;
t820 = t881 * t912 + t985;
t1023 = t820 * t909 - t913 * t892;
t1022 = t820 * t913 + t909 * t892;
t928 = t862 + t889;
t1021 = qJ(3) * t928;
t863 = t892 - t939;
t935 = t907 * qJDD(2) + t911 * t863;
t790 = -qJD(4) * t856 - t935;
t791 = -t854 * qJD(4) + t911 * qJDD(2) - t907 * t863;
t717 = t790 * t903 + t791 * t904;
t962 = qJD(1) * t908;
t885 = qJD(4) + t962;
t777 = t885 * t796;
t1019 = t717 - t777;
t685 = t717 + t777;
t993 = t854 * t885;
t750 = t791 + t993;
t1016 = -pkin(2) * t939 + t1034 * t962;
t1012 = pkin(2) + pkin(8);
t873 = pkin(3) * t962 - qJD(2) * pkin(8);
t874 = t909 * g(1) - t913 * g(2);
t924 = -qJDD(1) * pkin(1) - t874;
t917 = -t1016 + t924 - t1021;
t695 = -t873 * t962 + (-pkin(3) * t901 - pkin(7)) * t915 - t1012 * t863 + t917;
t875 = g(1) * t913 + g(2) * t909;
t840 = -pkin(1) * t915 + qJDD(1) * pkin(7) - t875;
t825 = t912 * g(3) + t908 * t840;
t1008 = qJ(3) * t908;
t1011 = pkin(2) * t912;
t930 = -t1008 - t1011;
t859 = t930 * qJD(1);
t922 = -qJDD(2) * pkin(2) - t914 * qJ(3) + t859 * t962 + qJDD(3) + t825;
t715 = -t871 * pkin(8) + (t862 - t889) * pkin(3) + t922;
t641 = t907 * t695 - t911 * t715;
t605 = pkin(4) * t1017 - qJ(5) * t750 - t641;
t1013 = t854 ^ 2;
t642 = t911 * t695 + t907 * t715;
t824 = pkin(4) * t885 - qJ(5) * t856;
t611 = -pkin(4) * t1013 + qJ(5) * t790 - t824 * t885 + t642;
t931 = t1033 * t798 + t904 * t605 - t903 * t611;
t547 = t1033 * t796 + t903 * t605 + t904 * t611;
t1015 = -t881 * t908 + t967;
t896 = t908 * g(3);
t1014 = -(qJD(1) * t859 + t840) * t912 + t914 * pkin(2) + t896;
t725 = t727 ^ 2;
t726 = t729 ^ 2;
t794 = t796 ^ 2;
t795 = t798 ^ 2;
t846 = t856 ^ 2;
t877 = qJD(6) + t885;
t876 = t877 ^ 2;
t882 = t885 ^ 2;
t963 = t900 + t901;
t866 = t963 * qJDD(1);
t869 = t893 + t974;
t803 = t866 * t909 + t869 * t913;
t1010 = pkin(6) * t803;
t1009 = t915 * pkin(7);
t520 = pkin(5) * t1018 - pkin(9) * t685 + t931;
t716 = t904 * t790 - t791 * t903;
t765 = pkin(5) * t885 - pkin(9) * t798;
t522 = -pkin(5) * t794 + pkin(9) * t716 - t765 * t885 + t547;
t471 = -t910 * t520 + t522 * t906;
t472 = t520 * t906 + t522 * t910;
t443 = -t471 * t910 + t472 * t906;
t1007 = t443 * t903;
t1006 = t443 * t904;
t499 = t547 * t903 + t904 * t931;
t1005 = t499 * t911;
t918 = qJDD(2) * qJ(3) - t1014;
t713 = t863 * pkin(3) - pkin(8) * t974 + (t1034 + t873) * qJD(2) + t918;
t647 = -t790 * pkin(4) - t1013 * qJ(5) + t856 * t824 + qJDD(5) + t713;
t574 = -t716 * pkin(5) - t794 * pkin(9) + t798 * t765 + t647;
t1004 = t574 * t906;
t1003 = t574 * t910;
t1002 = t647 * t903;
t1001 = t647 * t904;
t654 = t659 + t843;
t1000 = t654 * t906;
t999 = t654 * t910;
t711 = t735 + t847;
t998 = t711 * t903;
t997 = t711 * t904;
t996 = t798 * t885;
t839 = -t924 + t1009;
t995 = t839 * t908;
t994 = t839 * t912;
t992 = t861 * t908;
t988 = t864 * t912;
t984 = t877 * t906;
t983 = t877 * t910;
t978 = t885 * t903;
t977 = t885 * t904;
t976 = t885 * t907;
t975 = t885 * t911;
t973 = t907 * t499;
t972 = t907 * t713;
t780 = t807 + t847;
t971 = t907 * t780;
t970 = t911 * t713;
t969 = t911 * t780;
t965 = pkin(1) * t869 + pkin(7) * t866;
t958 = qJD(6) + t877;
t956 = qJD(3) * qJD(2);
t952 = -t846 - t882;
t951 = t908 * t659;
t950 = t912 * t659;
t949 = t908 * t735;
t948 = t912 * t735;
t947 = t908 * t807;
t946 = t912 * t807;
t444 = t471 * t906 + t910 * t472;
t500 = t904 * t547 - t903 * t931;
t937 = -t910 * t716 + t906 * t717;
t826 = t912 * t840 - t896;
t753 = t825 * t908 + t912 * t826;
t809 = -t874 * t909 - t913 * t875;
t934 = t909 * t945;
t933 = t913 * t945;
t868 = -t909 * t915 + t953;
t932 = -pkin(6) * t868 - g(3) * t909;
t929 = pkin(2) * t908 - qJ(3) * t912;
t572 = -t911 * t641 + t907 * t642;
t573 = t907 * t641 + t911 * t642;
t927 = t906 * t716 + t910 * t717;
t752 = t825 * t912 - t826 * t908;
t926 = t878 * t912 + t986;
t808 = t874 * t913 - t875 * t909;
t925 = t791 - t993;
t923 = t716 + t996;
t920 = (-qJD(6) + t877) * t729 - t937;
t629 = -qJD(6) * t727 + t927;
t919 = (-qJD(4) + t885) * t856 - t935;
t754 = t918 + 0.2e1 * t956;
t916 = t863 * pkin(2) + t1016 + t839;
t870 = -t893 + t974;
t867 = t913 * t915 + t954;
t857 = t929 * qJDD(1);
t849 = t963 * t957;
t837 = t912 * t847;
t836 = t908 * t847;
t835 = -pkin(6) * t867 + g(3) * t913;
t832 = -t846 + t882;
t831 = -t882 + t1013;
t830 = qJDD(2) * t909 + t849 * t913;
t829 = t862 * t912 - t900 * t957;
t828 = -qJDD(2) * t913 + t849 * t909;
t827 = -t863 * t908 - t901 * t957;
t812 = t928 * t908;
t810 = (t863 - t939) * t912;
t805 = t846 - t1013;
t804 = t866 * t913 - t869 * t909;
t802 = pkin(6) * t804;
t801 = t988 - t992;
t800 = t861 * t912 + t864 * t908;
t793 = -t882 - t1013;
t789 = t829 * t913 - t934;
t788 = t827 * t913 + t934;
t787 = t829 * t909 + t933;
t786 = t827 * t909 - t933;
t778 = -t846 - t1013;
t776 = -t795 + t882;
t775 = t794 - t882;
t764 = -t994 - t1031;
t763 = -t995 - t1038;
t762 = (t854 * t911 - t856 * t907) * t885;
t761 = (t854 * t907 + t856 * t911) * t885;
t760 = t801 * t913 - t870 * t909;
t759 = t801 * t909 + t870 * t913;
t757 = -t795 - t882;
t756 = t825 - t1039;
t755 = t826 - t1032;
t745 = (qJD(4) + t885) * t856 + t935;
t743 = qJ(3) * t869 + t922;
t742 = pkin(2) * t869 + t754;
t741 = -t791 * t911 + t856 * t976;
t740 = -t791 * t907 - t856 * t975;
t739 = t790 * t907 - t854 * t975;
t738 = -t790 * t911 - t854 * t976;
t737 = t916 + t1021;
t736 = -t761 * t908 + t837;
t734 = -t831 * t911 + t971;
t733 = t832 * t907 - t1028;
t732 = -t831 * t907 - t969;
t731 = -t832 * t911 - t1029;
t730 = t795 - t794;
t724 = -t907 * t952 - t969;
t723 = t911 * t952 - t971;
t722 = -t1009 + (-t863 - t864) * pkin(2) + t917;
t721 = -t882 - t794;
t720 = (t861 + t928) * qJ(3) + t916;
t719 = t753 * t913 - t839 * t909;
t718 = t753 * t909 + t839 * t913;
t709 = t911 * t793 - t1029;
t708 = t907 * t793 + t1028;
t706 = t877 * t727;
t703 = -t726 + t876;
t702 = t725 - t876;
t701 = (-t796 * t904 + t798 * t903) * t885;
t700 = (-t796 * t903 - t798 * t904) * t885;
t699 = -t726 - t876;
t698 = pkin(2) * t871 - qJ(3) * t880 + t1039 - t922;
t697 = -t740 * t908 + t946;
t696 = -t738 * t908 - t946;
t694 = t1032 + pkin(2) * t879 - 0.2e1 * t956 + (-qJDD(2) + t872) * qJ(3) + t1014;
t692 = -t794 - t795;
t691 = t754 * t912 + t908 * t922;
t690 = t754 * t908 - t912 * t922;
t689 = t907 * t750 + t911 * t919;
t688 = t745 * t911 + t907 * t925;
t687 = -t750 * t911 + t907 * t919;
t686 = t745 * t907 - t911 * t925;
t680 = -t716 + t996;
t679 = -pkin(2) * t992 + t720 * t912 + t1031;
t678 = -qJ(3) * t988 - t722 * t908 + t1038;
t677 = t717 * t904 - t798 * t978;
t676 = t717 * t903 + t798 * t977;
t675 = -t716 * t903 + t796 * t977;
t674 = t716 * t904 + t796 * t978;
t673 = -t742 * t908 + t743 * t912;
t672 = -t731 * t908 + t750 * t912;
t671 = -t732 * t908 + t912 * t919;
t670 = t775 * t904 - t998;
t669 = -t776 * t903 + t1026;
t668 = t775 * t903 + t997;
t667 = t776 * t904 + t1027;
t666 = t723 * t908 + t912 * t925;
t665 = -t757 * t903 - t997;
t664 = -t723 * t912 + t908 * t925;
t663 = t757 * t904 - t998;
t661 = t708 * t908 + t745 * t912;
t660 = -t708 * t912 + t745 * t908;
t658 = t726 - t725;
t657 = -t686 * t908 + t805 * t912;
t656 = -t876 - t725;
t652 = t721 * t904 - t1027;
t651 = t721 * t903 + t1026;
t650 = t687 * t908 + t778 * t912;
t649 = -t687 * t912 + t778 * t908;
t646 = (-t727 * t910 + t729 * t906) * t877;
t645 = (-t727 * t906 - t729 * t910) * t877;
t644 = t691 * t913 - t737 * t909;
t643 = t691 * t909 + t737 * t913;
t640 = t700 * t907 - t701 * t911;
t639 = -t700 * t911 - t701 * t907;
t637 = -t639 * t908 + t837;
t636 = -pkin(1) * t690 + pkin(2) * t922 - qJ(3) * t754;
t635 = t666 * t913 + t724 * t909;
t634 = t666 * t909 - t724 * t913;
t633 = -t725 - t726;
t632 = t661 * t913 + t709 * t909;
t631 = t661 * t909 - t709 * t913;
t630 = pkin(3) * t687 - qJ(3) * t689;
t628 = -qJD(6) * t729 - t937;
t627 = -pkin(7) * t690 - t737 * t929;
t626 = t685 * t903 + t904 * t923;
t625 = -t1019 * t903 - t680 * t904;
t624 = -t685 * t904 + t903 * t923;
t623 = t1019 * t904 - t680 * t903;
t621 = t702 * t910 - t1000;
t620 = -t703 * t906 + t1024;
t619 = t702 * t906 + t999;
t618 = t703 * t910 + t1025;
t617 = t676 * t907 - t677 * t911;
t616 = t674 * t907 - t675 * t911;
t615 = -t676 * t911 - t677 * t907;
t614 = -t674 * t911 - t675 * t907;
t613 = -t699 * t906 - t999;
t612 = t699 * t910 - t1000;
t609 = t668 * t907 - t670 * t911;
t608 = t667 * t907 - t669 * t911;
t607 = -t668 * t911 - t670 * t907;
t606 = -t667 * t911 - t669 * t907;
t604 = -t907 * t663 + t665 * t911;
t603 = t663 * t911 + t665 * t907;
t600 = t650 * t913 + t689 * t909;
t599 = t650 * t909 - t689 * t913;
t598 = pkin(3) * t925 - t1012 * t724 - t972;
t597 = pkin(3) * t745 - t1012 * t709 + t970;
t596 = -qJ(5) * t663 + t1001;
t595 = -t727 * t958 + t927;
t594 = t629 + t706;
t593 = t629 - t706;
t590 = t729 * t958 + t937;
t589 = t656 * t910 - t1025;
t588 = t656 * t906 + t1024;
t587 = t629 * t910 - t729 * t984;
t586 = t629 * t906 + t729 * t983;
t585 = -t628 * t906 + t727 * t983;
t584 = t628 * t910 + t727 * t984;
t583 = -t615 * t908 + t948;
t582 = -t614 * t908 - t948;
t581 = pkin(3) * t723 - qJ(3) * t724 - t642;
t580 = -t907 * t651 + t652 * t911;
t579 = t651 * t911 + t652 * t907;
t578 = -qJ(5) * t651 + t1002;
t577 = pkin(3) * t708 - qJ(3) * t709 - t641;
t576 = -t645 * t903 + t646 * t904;
t575 = t645 * t904 + t646 * t903;
t571 = -t606 * t908 + t685 * t912;
t570 = -t607 * t908 + t912 * t923;
t569 = t1019 * t912 + t603 * t908;
t568 = t1019 * t908 - t603 * t912;
t567 = -pkin(1) * t664 - qJ(3) * t925 + t1012 * t723 - t970;
t566 = -pkin(4) * t1019 + qJ(5) * t665 + t1002;
t565 = -pkin(1) * t660 - qJ(3) * t745 + t1012 * t708 - t972;
t564 = t572 * t908 + t713 * t912;
t563 = -t572 * t912 + t713 * t908;
t562 = -pkin(4) * t680 + qJ(5) * t652 - t1001;
t561 = t579 * t908 + t680 * t912;
t560 = -t579 * t912 + t680 * t908;
t559 = -t907 * t624 + t626 * t911;
t558 = t623 * t907 - t625 * t911;
t557 = t624 * t911 + t626 * t907;
t556 = -t623 * t911 - t625 * t907;
t555 = -t619 * t903 + t621 * t904;
t554 = -t618 * t903 + t620 * t904;
t553 = t619 * t904 + t621 * t903;
t552 = t618 * t904 + t620 * t903;
t551 = -t612 * t903 + t613 * t904;
t550 = t612 * t904 + t613 * t903;
t548 = -t556 * t908 + t730 * t912;
t544 = t557 * t908 + t692 * t912;
t543 = -t557 * t912 + t692 * t908;
t542 = -pkin(9) * t612 + t1003;
t541 = pkin(3) * t778 - t1012 * t689 - t573;
t540 = t594 * t906 + t910 * t920;
t539 = -t590 * t910 - t593 * t906;
t538 = -t594 * t910 + t906 * t920;
t537 = -t590 * t906 + t593 * t910;
t536 = -t588 * t903 + t589 * t904;
t535 = t588 * t904 + t589 * t903;
t534 = -t586 * t903 + t587 * t904;
t533 = -t584 * t903 + t585 * t904;
t532 = t586 * t904 + t587 * t903;
t531 = t584 * t904 + t585 * t903;
t529 = -pkin(9) * t588 + t1004;
t528 = t569 * t913 + t604 * t909;
t527 = t569 * t909 - t604 * t913;
t526 = t575 * t907 - t576 * t911;
t525 = -t575 * t911 - t576 * t907;
t524 = -t525 * t908 + t843 * t912;
t523 = pkin(3) * t572 - qJ(3) * t573;
t521 = -pkin(7) * t664 + t581 * t912 - t598 * t908;
t518 = -pkin(7) * t660 + t577 * t912 - t597 * t908;
t517 = t561 * t913 + t580 * t909;
t516 = t561 * t909 - t580 * t913;
t515 = -pkin(1) * t649 - qJ(3) * t778 + t1012 * t687 + t572;
t514 = pkin(3) * t713 - t1012 * t573;
t513 = t564 * t913 + t573 * t909;
t512 = t564 * t909 - t573 * t913;
t511 = -pkin(5) * t595 + pkin(9) * t613 + t1004;
t510 = -pkin(5) * t590 + pkin(9) * t589 - t1003;
t509 = t553 * t907 - t555 * t911;
t508 = t552 * t907 - t554 * t911;
t507 = -t553 * t911 - t555 * t907;
t506 = -t552 * t911 - t554 * t907;
t505 = -pkin(7) * t649 - t541 * t908 + t630 * t912;
t504 = -t907 * t550 + t551 * t911;
t503 = t550 * t911 + t551 * t907;
t502 = t544 * t913 + t559 * t909;
t501 = t544 * t909 - t559 * t913;
t497 = pkin(3) * t557 + pkin(4) * t624 - qJ(3) * t559;
t496 = -pkin(1) * t563 - qJ(3) * t713 + t1012 * t572;
t495 = pkin(3) * t603 + pkin(4) * t663 - qJ(3) * t604 - t547;
t494 = -t538 * t903 + t540 * t904;
t493 = -t537 * t903 + t539 * t904;
t492 = t538 * t904 + t540 * t903;
t491 = t537 * t904 + t539 * t903;
t489 = -t907 * t535 + t536 * t911;
t488 = t535 * t911 + t536 * t907;
t487 = t532 * t907 - t534 * t911;
t486 = t531 * t907 - t533 * t911;
t485 = -t532 * t911 - t534 * t907;
t484 = -t531 * t911 - t533 * t907;
t483 = -pkin(4) * t647 + qJ(5) * t500;
t482 = -t506 * t908 + t594 * t912;
t481 = -t507 * t908 + t912 * t920;
t480 = t503 * t908 + t595 * t912;
t479 = -t503 * t912 + t595 * t908;
t478 = pkin(3) * t579 + pkin(4) * t651 - qJ(3) * t580 + t931;
t477 = -qJ(5) * t624 - t499;
t476 = pkin(3) * t1019 - t1012 * t604 - t911 * t566 - t907 * t596;
t475 = -t485 * t908 + t950;
t474 = -t484 * t908 - t950;
t473 = -pkin(4) * t692 + qJ(5) * t626 + t500;
t469 = pkin(3) * t680 - t1012 * t580 - t911 * t562 - t907 * t578;
t468 = t488 * t908 + t590 * t912;
t467 = -t488 * t912 + t590 * t908;
t466 = -pkin(1) * t568 - qJ(3) * t1019 + t1012 * t603 + t907 * t566 - t911 * t596;
t465 = -pkin(7) * t563 - t514 * t908 + t523 * t912;
t464 = -pkin(1) * t560 - qJ(3) * t680 + t1012 * t579 + t907 * t562 - t911 * t578;
t463 = -qJ(5) * t550 - t511 * t903 + t542 * t904;
t462 = t500 * t911 - t973;
t461 = t500 * t907 + t1005;
t460 = -qJ(5) * t535 - t510 * t903 + t529 * t904;
t459 = -pkin(4) * t595 + qJ(5) * t551 + t511 * t904 + t542 * t903;
t458 = t480 * t913 + t504 * t909;
t457 = t480 * t909 - t504 * t913;
t456 = t461 * t908 + t647 * t912;
t455 = -t461 * t912 + t647 * t908;
t454 = -t907 * t492 + t494 * t911;
t453 = t491 * t907 - t493 * t911;
t452 = t492 * t911 + t494 * t907;
t451 = -t491 * t911 - t493 * t907;
t450 = -pkin(4) * t590 + qJ(5) * t536 + t510 * t904 + t529 * t903;
t449 = -t451 * t908 + t658 * t912;
t448 = t452 * t908 + t633 * t912;
t447 = -t452 * t912 + t633 * t908;
t446 = t468 * t913 + t489 * t909;
t445 = t468 * t909 - t489 * t913;
t442 = -pkin(7) * t568 - t476 * t908 + t495 * t912;
t441 = -pkin(5) * t574 + pkin(9) * t444;
t440 = -pkin(7) * t560 - t469 * t908 + t478 * t912;
t439 = -pkin(9) * t538 - t443;
t438 = -pkin(5) * t633 + pkin(9) * t540 + t444;
t437 = pkin(3) * t692 - t1012 * t559 - t911 * t473 - t907 * t477;
t436 = pkin(3) * t503 + pkin(4) * t550 + pkin(5) * t612 - qJ(3) * t504 - t472;
t435 = -pkin(1) * t543 - qJ(3) * t692 + t1012 * t557 + t907 * t473 - t911 * t477;
t434 = t456 * t913 + t462 * t909;
t433 = t456 * t909 - t462 * t913;
t432 = pkin(3) * t488 + pkin(4) * t535 + pkin(5) * t588 - qJ(3) * t489 - t471;
t431 = t448 * t913 + t454 * t909;
t430 = t448 * t909 - t454 * t913;
t429 = pkin(3) * t461 + pkin(4) * t499 - qJ(3) * t462;
t428 = -pkin(7) * t543 - t437 * t908 + t497 * t912;
t427 = t444 * t904 - t1007;
t426 = t444 * t903 + t1006;
t424 = pkin(3) * t595 - t1012 * t504 - t911 * t459 - t907 * t463;
t423 = pkin(3) * t452 + pkin(4) * t492 + pkin(5) * t538 - qJ(3) * t454;
t422 = pkin(3) * t647 + qJ(5) * t973 - t1012 * t462 - t911 * t483;
t421 = pkin(3) * t590 - t1012 * t489 - t911 * t450 - t907 * t460;
t420 = -pkin(1) * t479 - qJ(3) * t595 + t1012 * t503 + t907 * t459 - t911 * t463;
t419 = -qJ(5) * t492 - t438 * t903 + t439 * t904;
t418 = -pkin(4) * t633 + qJ(5) * t494 + t438 * t904 + t439 * t903;
t417 = -pkin(1) * t467 - qJ(3) * t590 + t1012 * t488 + t907 * t450 - t911 * t460;
t416 = -pkin(1) * t455 - qJ(3) * t647 + qJ(5) * t1005 + t1012 * t461 + t907 * t483;
t415 = -t907 * t426 + t427 * t911;
t414 = t426 * t911 + t427 * t907;
t413 = -pkin(9) * t1006 - qJ(5) * t426 - t441 * t903;
t412 = t414 * t908 + t574 * t912;
t411 = -t414 * t912 + t574 * t908;
t410 = -pkin(7) * t479 - t424 * t908 + t436 * t912;
t409 = -pkin(4) * t574 - pkin(9) * t1007 + qJ(5) * t427 + t441 * t904;
t408 = -pkin(7) * t467 - t421 * t908 + t432 * t912;
t407 = -pkin(7) * t455 - t422 * t908 + t429 * t912;
t406 = pkin(3) * t633 - t1012 * t454 - t911 * t418 - t907 * t419;
t405 = t412 * t913 + t415 * t909;
t404 = t412 * t909 - t415 * t913;
t403 = -pkin(1) * t447 - qJ(3) * t633 + t1012 * t452 + t907 * t418 - t911 * t419;
t402 = -pkin(7) * t447 - t406 * t908 + t423 * t912;
t401 = pkin(3) * t414 + pkin(4) * t426 + pkin(5) * t443 - qJ(3) * t415;
t400 = pkin(3) * t574 - t1012 * t415 - t911 * t409 - t907 * t413;
t399 = -pkin(1) * t411 - qJ(3) * t574 + t1012 * t414 + t907 * t409 - t911 * t413;
t398 = -pkin(7) * t411 - t400 * t908 + t401 * t912;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t867, -t868, 0, t809, 0, 0, 0, 0, 0, 0, -t771, t770, t804, t719, 0, 0, 0, 0, 0, 0, t804, t771, -t770, t644, 0, 0, 0, 0, 0, 0, t632, t635, t600, t513, 0, 0, 0, 0, 0, 0, t517, t528, t502, t434, 0, 0, 0, 0, 0, 0, t446, t458, t431, t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t868, -t867, 0, t808, 0, 0, 0, 0, 0, 0, -t767, t766, t803, t718, 0, 0, 0, 0, 0, 0, t803, t767, -t766, t643, 0, 0, 0, 0, 0, 0, t631, t634, t599, t512, 0, 0, 0, 0, 0, 0, t516, t527, t501, t433, 0, 0, 0, 0, 0, 0, t445, t457, t430, t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t814, t811, 0, -t752, 0, 0, 0, 0, 0, 0, 0, -t814, -t811, t690, 0, 0, 0, 0, 0, 0, t660, t664, t649, t563, 0, 0, 0, 0, 0, 0, t560, t568, t543, t455, 0, 0, 0, 0, 0, 0, t467, t479, t447, t411; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t868, 0, -t867, 0, t932, -t835, -t808, -pkin(6) * t808, t789, t760, t1035, t788, t1022, t830, -t756 * t909 + t763 * t913 + t1041, -t755 * t909 + t764 * t913 - t1043, t752 * t913 - t1010, -pkin(6) * t718 - (pkin(1) * t909 - pkin(7) * t913) * t752, t830, -t1035, -t1022, t789, t760, t788, t673 * t913 - t857 * t909 - t1010, t678 * t913 - t698 * t909 - t1041, t679 * t913 - t694 * t909 + t1043, -pkin(6) * t643 + t627 * t913 - t636 * t909, t697 * t913 - t741 * t909, t657 * t913 - t688 * t909, t672 * t913 - t733 * t909, t696 * t913 - t739 * t909, t671 * t913 - t734 * t909, t736 * t913 - t762 * t909, -pkin(6) * t631 + t518 * t913 - t565 * t909, -pkin(6) * t634 + t521 * t913 - t567 * t909, -pkin(6) * t599 + t505 * t913 - t515 * t909, -pkin(6) * t512 + t465 * t913 - t496 * t909, t583 * t913 - t617 * t909, t548 * t913 - t558 * t909, t571 * t913 - t608 * t909, t582 * t913 - t616 * t909, t570 * t913 - t609 * t909, t637 * t913 - t640 * t909, -pkin(6) * t516 + t440 * t913 - t464 * t909, -pkin(6) * t527 + t442 * t913 - t466 * t909, -pkin(6) * t501 + t428 * t913 - t435 * t909, -pkin(6) * t433 + t407 * t913 - t416 * t909, t475 * t913 - t487 * t909, t449 * t913 - t453 * t909, t482 * t913 - t508 * t909, t474 * t913 - t486 * t909, t481 * t913 - t509 * t909, t524 * t913 - t526 * t909, -pkin(6) * t445 + t408 * t913 - t417 * t909, -pkin(6) * t457 + t410 * t913 - t420 * t909, -pkin(6) * t430 + t402 * t913 - t403 * t909, -pkin(6) * t404 + t398 * t913 - t399 * t909; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t867, 0, t868, 0, t835, t932, t809, pkin(6) * t809, t787, t759, t1036, t786, t1023, t828, t756 * t913 + t763 * t909 - t1040, t755 * t913 + t764 * t909 + t1042, t752 * t909 + t802, pkin(6) * t719 - (-pkin(1) * t913 - pkin(7) * t909) * t752, t828, -t1036, -t1023, t787, t759, t786, t673 * t909 + t857 * t913 + t802, t678 * t909 + t698 * t913 + t1040, t679 * t909 + t694 * t913 - t1042, pkin(6) * t644 + t627 * t909 + t636 * t913, t697 * t909 + t741 * t913, t657 * t909 + t688 * t913, t672 * t909 + t733 * t913, t696 * t909 + t739 * t913, t671 * t909 + t734 * t913, t736 * t909 + t762 * t913, pkin(6) * t632 + t518 * t909 + t565 * t913, pkin(6) * t635 + t521 * t909 + t567 * t913, pkin(6) * t600 + t505 * t909 + t515 * t913, pkin(6) * t513 + t465 * t909 + t496 * t913, t583 * t909 + t617 * t913, t548 * t909 + t558 * t913, t571 * t909 + t608 * t913, t582 * t909 + t616 * t913, t570 * t909 + t609 * t913, t637 * t909 + t640 * t913, pkin(6) * t517 + t440 * t909 + t464 * t913, pkin(6) * t528 + t442 * t909 + t466 * t913, pkin(6) * t502 + t428 * t909 + t435 * t913, pkin(6) * t434 + t407 * t909 + t416 * t913, t475 * t909 + t487 * t913, t449 * t909 + t453 * t913, t482 * t909 + t508 * t913, t474 * t909 + t486 * t913, t481 * t909 + t509 * t913, t524 * t909 + t526 * t913, pkin(6) * t446 + t408 * t909 + t417 * t913, pkin(6) * t458 + t410 * t909 + t420 * t913, pkin(6) * t431 + t402 * t909 + t403 * t913, pkin(6) * t405 + t398 * t909 + t399 * t913; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t874, t875, 0, 0, t812, t800, t926, t810, -t1015, 0, pkin(1) * t864 - t1030 + t994, -pkin(1) * t861 + t1037 - t995, t753 + t965, pkin(1) * t839 + pkin(7) * t753, 0, -t926, t1015, t812, t800, t810, t742 * t912 + t743 * t908 + t965, t1030 + t912 * t722 + (-pkin(1) - t1008) * t864, -t1037 + t908 * t720 + (pkin(1) + t1011) * t861, pkin(7) * t691 + (pkin(1) - t930) * t737, t740 * t912 + t947, t686 * t912 + t805 * t908, t731 * t912 + t750 * t908, t738 * t912 - t947, t732 * t912 + t908 * t919, t761 * t912 + t836, -pkin(1) * t709 + pkin(7) * t661 + t577 * t908 + t597 * t912, -pkin(1) * t724 + pkin(7) * t666 + t581 * t908 + t598 * t912, -pkin(1) * t689 + pkin(7) * t650 + t541 * t912 + t630 * t908, -pkin(1) * t573 + pkin(7) * t564 + t514 * t912 + t523 * t908, t615 * t912 + t949, t556 * t912 + t730 * t908, t606 * t912 + t685 * t908, t614 * t912 - t949, t607 * t912 + t908 * t923, t639 * t912 + t836, -pkin(1) * t580 + pkin(7) * t561 + t469 * t912 + t478 * t908, -pkin(1) * t604 + pkin(7) * t569 + t476 * t912 + t495 * t908, -pkin(1) * t559 + pkin(7) * t544 + t437 * t912 + t497 * t908, -pkin(1) * t462 + pkin(7) * t456 + t422 * t912 + t429 * t908, t485 * t912 + t951, t451 * t912 + t658 * t908, t506 * t912 + t594 * t908, t484 * t912 - t951, t507 * t912 + t908 * t920, t525 * t912 + t843 * t908, -pkin(1) * t489 + pkin(7) * t468 + t421 * t912 + t432 * t908, -pkin(1) * t504 + pkin(7) * t480 + t424 * t912 + t436 * t908, -pkin(1) * t454 + pkin(7) * t448 + t406 * t912 + t423 * t908, -pkin(1) * t415 + pkin(7) * t412 + t400 * t912 + t401 * t908;];
tauB_reg  = t1;