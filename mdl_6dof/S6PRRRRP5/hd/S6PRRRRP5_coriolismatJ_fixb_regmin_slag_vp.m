% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x27]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:57
% EndTime: 2019-03-09 00:26:39
% DurationCPUTime: 23.37s
% Computational Cost: add. (18189->1009), mult. (48305->1481), div. (0->0), fcn. (53961->12), ass. (0->701)
t1032 = cos(pkin(7));
t677 = sin(qJ(3));
t909 = pkin(2) * t1032;
t848 = t677 * t909;
t674 = sin(pkin(7));
t681 = cos(qJ(3));
t985 = t674 * t681;
t601 = pkin(9) * t985 + t848;
t581 = pkin(10) * t1032 + t601;
t680 = cos(qJ(4));
t970 = t680 * t581;
t826 = -pkin(3) * t681 - pkin(10) * t677;
t582 = (-pkin(2) + t826) * t674;
t676 = sin(qJ(4));
t977 = t676 * t582;
t391 = t970 + t977;
t370 = -pkin(11) * t985 + t391;
t675 = sin(qJ(5));
t679 = cos(qJ(5));
t986 = t674 * t677;
t599 = pkin(9) * t986 - t681 * t909;
t580 = -pkin(3) * t1032 + t599;
t596 = -t1032 * t680 + t676 * t986;
t859 = t676 * t1032;
t598 = t680 * t986 + t859;
t825 = pkin(4) * t596 - pkin(11) * t598;
t691 = t580 + t825;
t193 = t370 * t679 + t675 * t691;
t971 = t679 * t681;
t912 = t674 * t971;
t981 = t675 * t598;
t509 = t912 + t981;
t151 = -qJ(6) * t509 + t193;
t824 = -pkin(4) * t680 - pkin(11) * t676;
t636 = -pkin(3) + t824;
t613 = t679 * t636;
t976 = t676 * t679;
t823 = -qJ(6) * t976 + t613;
t493 = (-pkin(10) * t675 - pkin(5)) * t680 + t823;
t1084 = -t493 / 0.2e1;
t978 = t675 * t680;
t924 = pkin(10) * t978;
t514 = t823 - t924;
t863 = t514 / 0.2e1 + t1084;
t1104 = t863 * t151;
t1033 = cos(pkin(6));
t1048 = cos(qJ(2));
t1031 = sin(pkin(6));
t814 = t1032 * t1031;
t1099 = t1033 * t674 + t1048 * t814;
t678 = sin(qJ(2));
t858 = t678 * t1031;
t832 = t677 * t858;
t504 = -t1099 * t681 + t832;
t1008 = t504 * t675;
t830 = t681 * t858;
t505 = t1099 * t677 + t830;
t1003 = t505 * t680;
t829 = t1031 * t1048;
t693 = t1032 * t1033 - t674 * t829;
t689 = t693 * t676;
t366 = t689 + t1003;
t1018 = t366 * t679;
t243 = t1008 + t1018;
t1103 = t863 * t243;
t390 = t581 * t676 - t582 * t680;
t369 = pkin(4) * t985 + t390;
t1102 = (t390 / 0.2e1 - t369 / 0.2e1) * t675;
t1101 = -t493 + t514;
t690 = t693 * t674;
t777 = -t829 / 0.2e1;
t1100 = t690 / 0.2e1 + t777;
t416 = t679 * t596;
t871 = -t416 / 0.2e1;
t914 = t675 * t985;
t974 = t679 * t598;
t511 = -t914 + t974;
t999 = t511 * t680;
t1098 = t676 * t871 + t999 / 0.2e1;
t1002 = t509 * t680;
t980 = t675 * t676;
t876 = -t980 / 0.2e1;
t736 = t596 * t876 + t1002 / 0.2e1;
t669 = t675 ^ 2;
t671 = t679 ^ 2;
t653 = t671 - t669;
t1097 = qJD(4) * t653;
t1049 = t680 / 0.2e1;
t1054 = -t676 / 0.2e1;
t881 = t596 * t1054;
t751 = t1049 * t598 + t881;
t1096 = t751 * qJD(4);
t926 = t680 * qJD(3);
t657 = t676 * t926;
t853 = qJD(2) * t751 + t657;
t932 = t596 * qJD(2);
t747 = -qJD(3) * t751 + t598 * t932;
t1095 = t511 ^ 2;
t1094 = t596 ^ 2;
t463 = pkin(4) * t598 + pkin(11) * t596;
t432 = t679 * t463;
t984 = t675 * t390;
t855 = t432 + t984;
t161 = pkin(5) * t598 + qJ(6) * t416 + t855;
t1093 = t161 / 0.2e1;
t1038 = t680 * pkin(11);
t1041 = t676 * pkin(4);
t640 = -t1038 + t1041;
t472 = t848 + (pkin(9) + t640) * t985;
t445 = t679 * t472;
t600 = (pkin(3) * t677 - pkin(10) * t681) * t674;
t587 = t676 * t600;
t588 = t680 * t599;
t964 = t587 - t588;
t398 = pkin(11) * t986 + t964;
t983 = t675 * t398;
t259 = t445 - t983;
t967 = t680 * t681;
t979 = t675 * t677;
t563 = (t679 * t967 + t979) * t674;
t915 = t676 * t985;
t184 = pkin(5) * t915 - qJ(6) * t563 + t259;
t1092 = t184 / 0.2e1;
t1006 = t504 * t679;
t1019 = t366 * t675;
t242 = -t1006 + t1019;
t1091 = -t242 / 0.2e1;
t1090 = t243 / 0.2e1;
t1004 = t505 * t679;
t302 = t504 * t978 + t1004;
t1089 = t302 / 0.2e1;
t643 = t679 * t986;
t562 = t680 * t914 - t643;
t1044 = t562 * pkin(5);
t586 = t676 * t599;
t968 = t680 * t600;
t854 = t586 + t968;
t397 = -pkin(4) * t986 - t854;
t317 = t397 + t1044;
t1088 = t317 / 0.2e1;
t833 = t676 * t858;
t792 = t674 * t833;
t775 = t678 * t814;
t744 = t677 * t775;
t556 = t681 * t829 - t744;
t993 = t556 * t680;
t448 = t792 + t993;
t1012 = t448 * t675;
t743 = t681 * t775;
t555 = t677 * t829 + t743;
t996 = t555 * t679;
t321 = t996 - t1012;
t1087 = t321 / 0.2e1;
t1086 = t366 / 0.2e1;
t1085 = t432 / 0.2e1;
t494 = t511 * t675;
t885 = -t494 / 0.2e1;
t623 = t679 * t640;
t662 = pkin(10) * t980;
t963 = t662 + t623;
t972 = t679 * t680;
t502 = pkin(5) * t676 - qJ(6) * t972 + t963;
t1083 = t502 / 0.2e1;
t1082 = -t509 / 0.2e1;
t1081 = t509 / 0.2e1;
t1080 = -t511 / 0.2e1;
t1079 = t511 / 0.2e1;
t923 = pkin(10) * t972;
t565 = t636 * t675 + t923;
t515 = -qJ(6) * t980 + t565;
t1077 = -t515 / 0.2e1;
t1076 = t515 / 0.2e1;
t622 = t675 * t640;
t922 = pkin(10) * t976;
t834 = t622 - t922;
t516 = -qJ(6) * t978 + t834;
t1075 = t516 / 0.2e1;
t1074 = -t562 / 0.2e1;
t1073 = t562 / 0.2e1;
t1072 = -t563 / 0.2e1;
t1071 = t563 / 0.2e1;
t1070 = -t596 / 0.2e1;
t1069 = t596 / 0.2e1;
t1068 = -t598 / 0.2e1;
t1067 = t598 / 0.2e1;
t1066 = t622 / 0.2e1;
t1042 = t675 * pkin(5);
t910 = pkin(10) + t1042;
t625 = t910 * t676;
t1065 = -t625 / 0.2e1;
t1064 = t625 / 0.2e1;
t626 = t910 * t680;
t1063 = t626 / 0.2e1;
t1036 = -qJ(6) - pkin(11);
t637 = t1036 * t675;
t1062 = -t637 / 0.2e1;
t1061 = t637 / 0.2e1;
t638 = t1036 * t679;
t1060 = t638 / 0.2e1;
t1059 = -t638 / 0.2e1;
t1058 = t643 / 0.2e1;
t666 = -pkin(5) * t679 - pkin(4);
t1057 = -t666 / 0.2e1;
t1056 = -t675 / 0.2e1;
t1055 = t675 / 0.2e1;
t1053 = t676 / 0.2e1;
t1052 = -t679 / 0.2e1;
t1051 = t679 / 0.2e1;
t1050 = -t680 / 0.2e1;
t1047 = pkin(5) * t511;
t1046 = pkin(5) * t563;
t1045 = pkin(10) * t509;
t1043 = t596 * pkin(5);
t1040 = t680 * pkin(5);
t1039 = t680 * pkin(10);
t1037 = -qJD(5) / 0.2e1;
t1035 = pkin(5) * qJD(5);
t1034 = pkin(5) * qJD(6);
t192 = t370 * t675 - t679 * t691;
t150 = -qJ(6) * t511 - t192;
t128 = t150 + t1043;
t966 = t128 - t150;
t21 = t966 * t509;
t1030 = qJD(2) * t21;
t1029 = t128 * t675;
t1028 = t128 * t679;
t1027 = t151 * t675;
t1026 = t242 * t675;
t1025 = t242 * t680;
t1024 = t243 * t679;
t1023 = t243 * t680;
t1022 = t321 * t680;
t1011 = t448 * t679;
t998 = t555 * t675;
t322 = t998 + t1011;
t1021 = t322 * t680;
t365 = t505 * t676 - t680 * t693;
t1020 = t365 * t511;
t186 = t365 * t675;
t1017 = t369 * t679;
t1016 = t397 * t675;
t1015 = t397 * t679;
t831 = t680 * t858;
t994 = t556 * t676;
t447 = -t674 * t831 + t994;
t1014 = t447 * t675;
t1013 = t447 * t679;
t1007 = t504 * t676;
t1005 = t505 * t675;
t303 = -t504 * t972 + t1005;
t49 = -t1007 * t365 - t242 * t302 + t243 * t303;
t1010 = t49 * qJD(1);
t50 = -t242 * t321 + t243 * t322 + t365 * t447;
t1009 = t50 * qJD(1);
t51 = (-t1024 + t366 - t1026) * t365;
t1001 = t51 * qJD(1);
t1000 = t511 * t679;
t997 = t555 * t676;
t995 = t555 * t680;
t992 = t563 * t675;
t991 = t565 * t598;
t990 = t580 * t680;
t668 = t674 ^ 2;
t673 = t681 ^ 2;
t989 = t668 * t673;
t988 = t668 * t677;
t670 = t676 ^ 2;
t987 = t670 * t675;
t982 = t675 * t509;
t412 = t675 * t596;
t975 = t679 * t509;
t973 = t679 * t670;
t969 = t680 * t596;
t379 = t679 * t390;
t431 = t675 * t463;
t965 = t379 - t431;
t389 = t679 * t398;
t444 = t675 * t472;
t260 = t389 + t444;
t652 = t671 + t669;
t672 = t680 ^ 2;
t654 = t672 - t670;
t869 = t972 / 0.2e1;
t717 = t511 * t869 + t671 * t881;
t268 = -t992 / 0.2e1 + t717;
t962 = qJD(2) * t268;
t275 = -t1094 * t675 + t509 * t598;
t961 = qJD(2) * t275;
t276 = -t1094 * t679 + t511 * t598;
t960 = qJD(2) * t276;
t878 = -t985 / 0.2e1;
t839 = t675 * t878;
t870 = t416 / 0.2e1;
t278 = t1058 + t676 * t870 + (t839 + t1080) * t680;
t959 = qJD(2) * t278;
t868 = t971 / 0.2e1;
t716 = (t680 * t868 + t979 / 0.2e1) * t674;
t279 = t716 - t736;
t958 = qJD(2) * t279;
t957 = qJD(2) * t412;
t956 = qJD(2) * t416;
t955 = qJD(2) * t511;
t954 = qJD(2) * t681;
t953 = qJD(3) * t674;
t952 = qJD(3) * t681;
t951 = qJD(4) * t365;
t950 = qJD(4) * t675;
t949 = qJD(4) * t676;
t948 = qJD(4) * t679;
t947 = qJD(4) * t680;
t946 = qJD(5) * t243;
t945 = qJD(5) * t509;
t944 = qJD(5) * t596;
t943 = qJD(5) * t675;
t942 = qJD(5) * t679;
t941 = qJD(5) * t680;
t530 = t416 * t980;
t149 = t530 + (-t1002 / 0.2e1 + t1072) * t679 + (-t999 / 0.2e1 + t1073) * t675;
t940 = t149 * qJD(2);
t219 = -t509 * t563 - t511 * t562;
t939 = t219 * qJD(2);
t238 = (t975 + t494) * t596;
t938 = t238 * qJD(2);
t286 = -t509 * t915 - t562 * t596;
t937 = t286 * qJD(2);
t287 = t511 * t915 + t563 * t596;
t936 = t287 * qJD(2);
t349 = -t598 * t676 - t969;
t376 = t349 * t985;
t935 = t376 * qJD(2);
t495 = -t596 * t986 + t676 * t989;
t934 = t495 * qJD(2);
t496 = -t598 * t986 + t680 * t989;
t933 = t496 * qJD(2);
t931 = t596 * qJD(4);
t930 = t598 * qJD(4);
t608 = (-t677 ^ 2 + t673) * t668;
t929 = t608 * qJD(2);
t928 = t674 * qJD(4);
t927 = t676 * qJD(3);
t925 = pkin(5) * t976;
t921 = pkin(5) * t943;
t920 = pkin(5) * t942;
t919 = t1047 / 0.2e1;
t918 = -t1043 / 0.2e1;
t917 = t1042 / 0.2e1;
t916 = -t1040 / 0.2e1;
t913 = t675 * t969;
t911 = t679 * t969;
t908 = t511 * t932;
t906 = t674 * t954;
t905 = t674 * t952;
t904 = t679 * t927;
t903 = t681 * t928;
t902 = t675 * t948;
t901 = t676 * t948;
t900 = t675 * t941;
t899 = t679 * t941;
t898 = t668 * t954;
t897 = t677 * t953;
t896 = t675 * t942;
t895 = t676 * t947;
t894 = t151 * t1052;
t893 = t1020 / 0.2e1;
t892 = t369 * t1055;
t891 = -t1014 / 0.2e1;
t890 = -t1013 / 0.2e1;
t889 = t504 * t1070;
t888 = t504 * t1067;
t466 = -t1007 / 0.2e1;
t467 = t1007 / 0.2e1;
t887 = t509 * t1054;
t884 = t511 * t1053;
t882 = t994 / 0.2e1;
t880 = -t986 / 0.2e1;
t879 = t986 / 0.2e1;
t877 = t985 / 0.2e1;
t875 = t980 / 0.2e1;
t874 = t978 / 0.2e1;
t873 = -t976 / 0.2e1;
t872 = t976 / 0.2e1;
t867 = t150 / 0.2e1 - t128 / 0.2e1;
t866 = t379 / 0.2e1 - t431 / 0.2e1;
t865 = -t389 / 0.2e1 - t444 / 0.2e1;
t862 = t587 / 0.2e1 - t588 / 0.2e1;
t861 = t504 * t1032;
t860 = t505 * t1032;
t856 = t1032 * qJD(2);
t852 = pkin(5) * t872;
t851 = pkin(10) * t877;
t850 = -qJD(5) - t932;
t849 = -qJD(5) + t926;
t847 = t675 * t901;
t846 = t952 * t988;
t845 = t677 * t898;
t844 = t675 * t904;
t843 = t680 * t906;
t842 = t365 * t872;
t841 = t676 * t878;
t840 = t676 * t877;
t838 = t675 * t877;
t836 = t674 * t868;
t828 = t674 * t856;
t827 = t1032 * t953;
t822 = t867 * t243;
t821 = 0.2e1 * t844;
t820 = t445 / 0.2e1 - t983 / 0.2e1;
t819 = t505 / 0.2e1 + t365 * t1054;
t818 = t586 / 0.2e1 + t968 / 0.2e1;
t817 = t1008 / 0.2e1 - t243 / 0.2e1;
t816 = t1006 / 0.2e1 + t242 / 0.2e1;
t815 = -t1054 * t638 + t1084;
t813 = pkin(11) * t841;
t812 = -qJD(4) + t906;
t208 = -qJ(6) * t562 + t260;
t266 = pkin(5) * t509 + t369;
t683 = t184 * t1091 + t208 * t1090 + t128 * t1089 + t303 * t151 / 0.2e1 + t365 * t1088 + t266 * t466;
t707 = t1065 * t447 + t1077 * t322 + t1084 * t321;
t2 = t683 + t707;
t24 = t128 * t184 + t151 * t208 + t266 * t317;
t811 = qJD(1) * t2 + qJD(2) * t24;
t22 = t1047 * t266 - t151 * t966;
t8 = t822 + (t893 - t321 / 0.2e1) * pkin(5);
t810 = qJD(1) * t8 + qJD(2) * t22;
t207 = qJ(6) * t412 - t965;
t319 = -pkin(5) * t412 + t391;
t23 = t128 * t161 + t151 * t207 + t266 * t319;
t685 = (t894 + t1029 / 0.2e1 + t319 / 0.2e1) * t365 + t161 * t1091 + t207 * t1090 + t266 * t1086;
t706 = t1057 * t447 + t1060 * t322 + t1062 * t321;
t4 = t685 + t706;
t809 = qJD(1) * t4 + qJD(2) * t23;
t25 = -t161 * t511 - t207 * t509 + (t1027 + t1028) * t596;
t756 = t1070 * t242 + t1081 * t365;
t761 = t243 * t1069 - t1020 / 0.2e1;
t36 = (-t322 / 0.2e1 + t756) * t679 + (t1087 + t761) * t675;
t808 = t36 * qJD(1) + t25 * qJD(2);
t26 = -t128 * t563 - t151 * t562 - t184 * t511 - t208 * t509;
t692 = t1071 * t242 + t1074 * t243 + t1080 * t302 + t1082 * t303;
t731 = (t1052 * t321 + t1056 * t322) * t676;
t29 = t731 - t692;
t807 = -t29 * qJD(1) + t26 * qJD(2);
t735 = t1082 * t504 + t242 * t878;
t758 = t1069 * t302 + t1073 * t365;
t45 = t1022 / 0.2e1 + (t891 + t735) * t676 + t758;
t58 = -t192 * t915 + t259 * t596 + t369 * t562 + t397 * t509;
t806 = t45 * qJD(1) + t58 * qJD(2);
t734 = t1080 * t504 + t243 * t878;
t757 = t1070 * t303 + t1071 * t365;
t48 = -t1021 / 0.2e1 + (t890 + t734) * t676 + t757;
t59 = -t193 * t915 - t260 * t596 + t369 * t563 + t397 * t511;
t805 = t48 * qJD(1) + t59 * qJD(2);
t52 = -t192 * t598 + t391 * t509 + (t432 + (-t369 + t390) * t675) * t596;
t763 = t1068 * t242 + t1081 * t366;
t54 = t1013 / 0.2e1 + t763;
t804 = t54 * qJD(1) + t52 * qJD(2);
t53 = -t193 * t598 + t391 * t511 + (t965 - t1017) * t596;
t760 = t1068 * t243 + t1079 * t366;
t57 = t891 + t760;
t803 = t57 * qJD(1) + t53 * qJD(2);
t61 = -t128 * t511 - t151 * t509;
t778 = -t831 / 0.2e1;
t703 = t674 * t778 + t882;
t764 = t1079 * t242 + t1082 * t243;
t72 = t703 - t764;
t802 = -qJD(1) * t72 + qJD(2) * t61;
t801 = t637 * t675 + t638 * t679;
t800 = t856 + qJD(3);
t100 = t192 * t596 - t369 * t509;
t754 = -t1011 / 0.2e1 - t998 / 0.2e1;
t75 = t754 + t756;
t799 = qJD(1) * t75 - qJD(2) * t100;
t101 = -t193 * t596 + t369 * t511;
t755 = -t1012 / 0.2e1 + t996 / 0.2e1;
t74 = t755 + t761;
t798 = qJD(1) * t74 - qJD(2) * t101;
t133 = t390 * t986 - t580 * t915 - t596 * t601 + t854 * t985;
t738 = t1069 * t505 + t365 * t880;
t96 = t995 / 0.2e1 + t738;
t797 = t96 * qJD(1) - t133 * qJD(2);
t134 = t601 * t598 + (-t391 * t677 + (t964 + t990) * t681) * t674;
t737 = t1067 * t505 + t366 * t880;
t99 = -t997 / 0.2e1 + t737;
t796 = t99 * qJD(1) + t134 * qJD(2);
t121 = (t1040 / 0.2e1 + t863) * t679;
t17 = (t1043 / 0.2e1 - t867) * t679;
t795 = qJD(2) * t17 - qJD(3) * t121;
t155 = t1101 * t980;
t694 = -t509 * t863 - t867 * t980;
t16 = t1046 / 0.2e1 + t694;
t794 = -qJD(2) * t16 + qJD(3) * t155;
t793 = t849 * t676;
t791 = pkin(11) * t878 - t391 / 0.2e1;
t144 = t889 + t993 / 0.2e1 + (-t365 * t681 / 0.2e1 + t833 / 0.2e1) * t674;
t250 = -t390 * t985 - t580 * t596;
t790 = -qJD(1) * t144 - qJD(2) * t250;
t146 = t888 + t882 + (t1086 * t681 + t778) * t674;
t251 = -t391 * t985 - t580 * t598;
t789 = -qJD(1) * t146 + qJD(2) * t251;
t682 = t777 - t690 / 0.2e1;
t720 = -t743 / 0.2e1;
t271 = t720 + t860 / 0.2e1 + t682 * t677;
t518 = -pkin(2) * t988 - t1032 * t601;
t788 = qJD(1) * t271 - qJD(2) * t518;
t721 = t744 / 0.2e1;
t272 = t721 - t861 / 0.2e1 + t682 * t681;
t517 = pkin(2) * t668 * t681 - t1032 * t599;
t787 = qJD(1) * t272 + qJD(2) * t517;
t197 = -t911 + (t838 - t974 / 0.2e1 + t1080) * t676;
t621 = t672 * t679 - t973;
t786 = -qJD(2) * t197 - qJD(3) * t621;
t198 = t913 + (t836 + t981 / 0.2e1 + t1081) * t676;
t620 = t654 * t675;
t785 = -qJD(2) * t198 + qJD(3) * t620;
t308 = (-t982 - t1000) * t676;
t618 = t652 * t670;
t784 = qJD(2) * t308 - qJD(3) * t618;
t378 = -t598 ^ 2 + t1094;
t783 = qJD(2) * t378 + qJD(3) * t349;
t782 = qJD(2) * t349 + qJD(3) * t654;
t781 = t926 - t932;
t780 = qJD(2) * t598 + t927;
t323 = t494 - t975;
t779 = qJD(2) * t323 + qJD(4) * t652;
t776 = t1038 / 0.2e1 - t1041 / 0.2e1;
t774 = pkin(3) * t1068 + t1053 * t580;
t773 = pkin(4) * t1074 - t1015 / 0.2e1;
t772 = pkin(4) * t1072 + t1016 / 0.2e1;
t771 = pkin(10) * t1080 - t1017 / 0.2e1;
t589 = t680 * t879 + t859 / 0.2e1;
t770 = t589 * qJD(2) + t927 / 0.2e1;
t769 = t674 * t826;
t696 = pkin(3) * t1069 + t990 / 0.2e1 + pkin(10) * t841;
t246 = t696 + t862;
t768 = pkin(3) * t926 - qJD(2) * t246;
t767 = (t1054 * t637 + t1076) * t679;
t564 = -t613 + t924;
t766 = t1049 * t192 + t1070 * t564;
t765 = t1050 * t193 + t1069 * t565;
t762 = t1026 / 0.2e1 + t1024 / 0.2e1;
t759 = t1052 * t242 + t1055 * t243;
t753 = t1076 * t509 + t1079 * t493;
t752 = t1080 * t502 + t1082 * t516;
t750 = t1059 * t562 + t1061 * t563;
t749 = t977 / 0.2e1 + t970 / 0.2e1;
t748 = -t975 / 0.2e1 + t885;
t746 = -t904 - t950;
t745 = pkin(4) * t1080 + pkin(11) * t871;
t741 = t1037 * t676 + t853;
t740 = t776 * t675;
t739 = t776 * t679;
t152 = t493 * t502 + t515 * t516 + t625 * t626;
t684 = (t1052 * t515 + t1055 * t493 + t1063) * t365 + t502 * t1091 + t243 * t1075 + t366 * t1064;
t698 = t1060 * t303 + t1062 * t302 + t467 * t666;
t20 = t684 + t698;
t686 = t1063 * t266 + t1064 * t319 + t1075 * t151 + t1076 * t207 + t1083 * t128 + t1093 * t493;
t708 = t1059 * t208 + t1061 * t184 + t1088 * t666;
t5 = -t686 + t708;
t733 = qJD(1) * t20 - qJD(2) * t5 + qJD(3) * t152;
t732 = t759 * t676;
t730 = t1067 * t564 + t1070 * t963;
t704 = t1050 * t151 + t1054 * t207 + t1069 * t515;
t705 = t1050 * t128 + t1054 * t161 + t1069 * t493;
t10 = (-t208 / 0.2e1 + t705) * t679 + (t1092 + t704) * t675 + t750 + t752;
t153 = (t493 * t680 + t502 * t676) * t679 + (t515 * t680 + t516 * t676) * t675;
t63 = (t1025 / 0.2e1 - t303 / 0.2e1) * t679 + (-t1023 / 0.2e1 + t1089) * t675;
t729 = t63 * qJD(1) + t10 * qJD(2) - t153 * qJD(3);
t11 = -t867 * t515 - t1104 + (t1065 * t511 + t266 * t873 + t1092) * pkin(5);
t157 = t1101 * t515 + t625 * t925;
t27 = -t1103 + (t365 * t873 + t1089) * pkin(5);
t728 = -qJD(1) * t27 - qJD(2) * t11 + qJD(3) * t157;
t299 = (t493 * t679 + t515 * t675) * t676;
t695 = t1044 / 0.2e1 + pkin(4) * t880 - t818;
t39 = (t1027 / 0.2e1 + t1028 / 0.2e1) * t676 + t695 + t753;
t86 = t467 - t732;
t727 = qJD(1) * t86 - qJD(2) * t39 - qJD(3) * t299;
t31 = (t1085 - t1045 / 0.2e1 + t1102) * t680 + (t192 / 0.2e1 + (pkin(10) * t1069 + t791) * t675) * t676 + t730 + t773;
t363 = t564 * t676 + (-t662 + t623) * t680;
t69 = (-t1019 / 0.2e1 + t816) * t676;
t726 = -t69 * qJD(1) - t31 * qJD(2) - t363 * qJD(3);
t32 = t596 * t1066 + t991 / 0.2e1 + (t771 + t866) * t680 + (t193 / 0.2e1 + t791 * t679) * t676 + t772;
t364 = t622 * t680 + (-t565 + t923) * t676;
t68 = (-t1018 / 0.2e1 - t817) * t676;
t725 = -t68 * qJD(1) - t32 * qJD(2) + t364 * qJD(3);
t499 = -pkin(10) * t987 - t564 * t680;
t65 = (t1045 / 0.2e1 + t892) * t676 + t766 + t865;
t88 = -t675 * t819 + t680 * t816;
t724 = qJD(1) * t88 + qJD(2) * t65 - qJD(3) * t499;
t500 = -pkin(10) * t973 - t565 * t680;
t64 = t676 * t771 + t765 + t820;
t87 = t679 * t819 + t680 * t817;
t723 = qJD(1) * t87 + qJD(2) * t64 + qJD(3) * t500;
t13 = t867 * t638 + (t1056 * t266 + t1057 * t511 + t1093) * pkin(5);
t377 = t666 * t1042;
t91 = t863 * t638 + (t1056 * t625 + t666 * t873 + t1083) * pkin(5);
t722 = -qJD(2) * t13 - qJD(3) * t91 + qJD(4) * t377;
t210 = (-t982 + t1000) * t676;
t228 = -t975 + 0.2e1 * t885;
t506 = t509 ^ 2;
t255 = t506 - t1095;
t719 = qJD(2) * t255 - qJD(3) * t210 + qJD(4) * t228;
t350 = t506 + t1095;
t718 = qJD(2) * t350 - qJD(3) * t308 + qJD(4) * t323;
t248 = -t586 / 0.2e1 + (t851 - t600 / 0.2e1) * t680 + t774;
t309 = t467 + t466;
t715 = pkin(3) * t927 - qJD(1) * t309 - qJD(2) * t248;
t702 = pkin(4) * t1081 + t1017 / 0.2e1 + pkin(11) * t412 / 0.2e1;
t107 = t702 - t866;
t507 = t1066 - t740;
t714 = pkin(4) * t948 - qJD(2) * t107 - qJD(3) * t507;
t109 = -t432 / 0.2e1 - t1102 + t745;
t508 = -t623 / 0.2e1 + t739;
t713 = pkin(4) * t950 - qJD(2) * t109 - qJD(3) * t508;
t301 = t748 * t676;
t320 = -t982 / 0.2e1 + t1000 / 0.2e1;
t712 = -qJD(3) * t301 - qJD(4) * t320 + t509 * t955;
t602 = (t669 / 0.2e1 - t671 / 0.2e1) * t676;
t711 = qJD(2) * t320 - qJD(3) * t602 + t902;
t710 = qJD(5) * t589 + t747;
t709 = t1059 * t509 + t1061 * t511 + t894;
t619 = t653 * t670;
t701 = qJD(2) * t210 + qJD(3) * t619 + 0.2e1 * t847;
t700 = -qJD(2) * t228 - t1097 + t821;
t699 = qJD(3) * t675 * t973 - qJD(2) * t301 + qJD(4) * t602;
t196 = -t1039 / 0.2e1 + t767 + (t916 + t815) * t675;
t43 = (t918 + t128 / 0.2e1) * t675 + t709 + t749;
t687 = t1003 / 0.2e1 + t689 / 0.2e1;
t79 = t687 - t762;
t697 = -qJD(1) * t79 - qJD(2) * t43 + qJD(3) * t196 - qJD(4) * t801;
t664 = t949 / 0.2e1;
t641 = qJD(3) * t879;
t629 = -0.2e1 * t676 * t896;
t595 = t602 * qJD(5);
t583 = (t898 - t928 / 0.2e1) * t677;
t503 = qJD(3) * t840 + qJD(4) * t589;
t462 = t662 + t623 / 0.2e1 + t739;
t461 = t922 - t622 / 0.2e1 - t740;
t392 = (t746 - t955) * pkin(5);
t342 = t349 * qJD(4);
t318 = t323 * qJD(6);
t316 = t320 * qJD(5);
t312 = t504 * t680;
t310 = 0.2e1 * t467;
t300 = t308 * qJD(6);
t296 = t301 * qJD(5);
t281 = t680 * t839 + t1058 + t1098;
t280 = t716 + t736;
t274 = -t860 / 0.2e1 + t720 + t1100 * t677;
t273 = t861 / 0.2e1 + t721 + t1100 * t681;
t267 = t992 / 0.2e1 + t717;
t249 = t680 * t851 + t774 + t818;
t247 = t696 - t862;
t220 = t228 * qJD(5);
t209 = t210 * qJD(5);
t200 = t598 * t872 + t676 * t838 + t884 + t911;
t199 = t598 * t876 + t676 * t836 + t887 - t913;
t195 = t1039 / 0.2e1 + pkin(5) * t874 + t767 + t815 * t675;
t188 = t365 * t679;
t148 = t1051 * t563 + t1056 * t562 + t680 * t748 + t530;
t147 = t366 * t877 - t703 + t888;
t145 = t365 * t878 + t889 - t993 / 0.2e1 - t792 / 0.2e1;
t122 = (t863 + t916) * t679;
t110 = t892 + t984 / 0.2e1 + t1085 + t745;
t108 = t702 + t866;
t98 = t997 / 0.2e1 + t737;
t97 = -t995 / 0.2e1 + t738;
t92 = pkin(5) * t1083 + t1059 * t514 + t1060 * t493 + t625 * t917 + t666 * t852;
t90 = t1023 / 0.2e1 + t842 + t504 * t874 + t1004 / 0.2e1;
t89 = -t1025 / 0.2e1 + t365 * t876 + t504 * t869 - t1005 / 0.2e1;
t85 = t466 - t732;
t80 = t687 + t762;
t77 = t755 - t761;
t76 = t754 - t756;
t73 = t703 + t764;
t71 = t1054 * t243 + t366 * t872 + t504 * t876;
t70 = t1054 * t242 + t366 * t875 + t504 * t872;
t67 = pkin(10) * t884 + t369 * t872 - t765 + t820;
t66 = pkin(10) * t887 + t369 * t876 - t766 + t865;
t62 = t1051 * t303 + t1056 * t302 - t680 * t759;
t56 = t1014 / 0.2e1 + t760;
t55 = t890 + t763;
t47 = t1021 / 0.2e1 + t447 * t872 + t734 * t676 + t757;
t46 = -t1022 / 0.2e1 + t447 * t875 + t735 * t676 + t758;
t44 = -t1029 / 0.2e1 + t675 * t918 - t709 + t749;
t40 = t128 * t873 + t151 * t876 + t695 - t753;
t38 = pkin(5) * t186;
t35 = t322 * t1051 + t321 * t1056 + t759 * t596 + (t975 / 0.2e1 + t885) * t365;
t34 = t834 * t1070 - t991 / 0.2e1 + t965 * t1050 + t193 * t1054 + t391 * t872 + t369 * t869 + t679 * t813 + t772 + t1098 * pkin(10);
t33 = pkin(10) * t736 + t1050 * t855 + t1054 * t192 + t369 * t874 + t391 * t875 + t675 * t813 - t730 + t773;
t30 = t731 + t692;
t28 = t1103 + (t1089 + t842) * pkin(5);
t19 = t684 - t698;
t18 = t150 * t1051 - t1028 / 0.2e1 + pkin(5) * t870;
t15 = -t1046 / 0.2e1 + t694;
t14 = pkin(5) * t1093 + t1059 * t150 + t1060 * t128 + t266 * t917 + t666 * t919;
t12 = pkin(5) * t1092 + t150 * t1076 + t128 * t1077 + t266 * t852 + t625 * t919 + t1104;
t9 = t1051 * t208 + t1056 * t184 + t675 * t704 + t679 * t705 - t750 + t752;
t7 = t822 + (t1087 + t893) * pkin(5);
t6 = t686 + t708;
t3 = t685 - t706;
t1 = t683 - t707;
t37 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t50 + qJD(3) * t49 + qJD(4) * t51; 0, 0, -qJD(2) * t858, -qJD(2) * t829, 0, 0, 0, 0, 0 (-t1032 * t555 - t668 * t830) * qJD(2) + t274 * qJD(3) (-t1032 * t556 + t668 * t832) * qJD(2) + t273 * qJD(3), 0, 0, 0, 0, 0 (t447 * t985 + t555 * t596) * qJD(2) + t97 * qJD(3) + t147 * qJD(4) (t448 * t985 + t555 * t598) * qJD(2) + t98 * qJD(3) + t145 * qJD(4), 0, 0, 0, 0, 0 (t321 * t596 + t447 * t509) * qJD(2) + t46 * qJD(3) + t55 * qJD(4) + t77 * qJD(5) (-t322 * t596 + t447 * t511) * qJD(2) + t47 * qJD(3) + t56 * qJD(4) + t76 * qJD(5) (-t321 * t511 - t322 * t509) * qJD(2) + t30 * qJD(3) + t35 * qJD(4), t1009 + (t128 * t321 + t151 * t322 + t266 * t447) * qJD(2) + t1 * qJD(3) + t3 * qJD(4) + t7 * qJD(5) + t73 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t274 - qJD(3) * t505, qJD(2) * t273 + qJD(3) * t504, 0, 0, 0, 0, 0, qJD(2) * t97 + qJD(4) * t310 - t505 * t926, qJD(2) * t98 + qJD(4) * t312 + t505 * t927, 0, 0, 0, 0, 0, t46 * qJD(2) + (-t302 * t680 - t504 * t987) * qJD(3) + t70 * qJD(4) + t90 * qJD(5), t47 * qJD(2) + (t303 * t680 - t504 * t973) * qJD(3) + t71 * qJD(4) + t89 * qJD(5), t30 * qJD(2) + t62 * qJD(4) + (-t302 * t679 - t303 * t675) * t927, t1010 + t1 * qJD(2) + (-t1007 * t625 + t302 * t493 + t303 * t515) * qJD(3) + t19 * qJD(4) + t28 * qJD(5) + t85 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t147 + qJD(3) * t310 - qJD(4) * t366, qJD(2) * t145 + qJD(3) * t312 + t951, 0, 0, 0, 0, 0, qJD(2) * t55 + qJD(3) * t70 + qJD(5) * t186 - t366 * t948, qJD(2) * t56 + qJD(3) * t71 + qJD(5) * t188 + t366 * t950, t35 * qJD(2) + t62 * qJD(3) - t652 * t951, t1001 + t3 * qJD(2) + t19 * qJD(3) + (t365 * t801 + t366 * t666) * qJD(4) + t38 * qJD(5) + t80 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t77 + qJD(3) * t90 + qJD(4) * t186 - t946, qJD(2) * t76 + qJD(3) * t89 + qJD(4) * t188 + qJD(5) * t242, 0, -pkin(5) * t946 + qJD(2) * t7 + qJD(3) * t28 + qJD(4) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t73 + qJD(3) * t85 + qJD(4) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271 * qJD(3), -t272 * qJD(3), 0, 0, 0, 0, 0, qJD(3) * t96 + qJD(4) * t146, qJD(3) * t99 + qJD(4) * t144, 0, 0, 0, 0, 0, qJD(3) * t45 + qJD(4) * t54 - qJD(5) * t74, qJD(3) * t48 + qJD(4) * t57 - qJD(5) * t75, -qJD(3) * t29 + qJD(4) * t36, qJD(3) * t2 + qJD(4) * t4 + qJD(5) * t8 - qJD(6) * t72 - t1009; 0, 0, 0, 0, t846, t608 * qJD(3), t681 * t827, -t677 * t827, 0, t518 * qJD(3), -t517 * qJD(3) (t680 * t905 - t931) * t598, qJD(3) * t376 + qJD(4) * t378, -qJD(3) * t496 + t596 * t903, qJD(3) * t495 + t598 * t903, -t846, -qJD(3) * t133 - qJD(4) * t251, qJD(3) * t134 + qJD(4) * t250 (qJD(3) * t563 - t679 * t931 - t945) * t511, qJD(3) * t219 + qJD(4) * t238 + qJD(5) * t255, qJD(3) * t287 + qJD(4) * t276 - t509 * t944, qJD(3) * t286 - qJD(4) * t275 - t511 * t944 (t676 * t905 + t930) * t596, qJD(3) * t58 + qJD(4) * t52 + qJD(5) * t101, qJD(3) * t59 + qJD(4) * t53 + qJD(5) * t100, qJD(3) * t26 + qJD(4) * t25 + qJD(5) * t21 + qJD(6) * t350, qJD(3) * t24 + qJD(4) * t23 + qJD(5) * t22 + qJD(6) * t61; 0, 0, 0, 0, t845, t929, t800 * t985, -t800 * t986, 0, -qJD(3) * t601 - t788, qJD(3) * t599 - t787, t674 * t780 * t967 + t1096, t654 * t905 + t342 + t935, t676 * t897 - t933, t680 * t897 + t934, -t583 (-t601 * t680 + t676 * t769) * qJD(3) + t249 * qJD(4) + t797 (t601 * t676 + t680 * t769) * qJD(3) + t247 * qJD(4) + t796, qJD(4) * t267 + t296 + (t904 + t955) * t563, t939 + t148 * qJD(4) - t209 + (-t562 * t679 - t992) * t927, t936 + (-t563 * t680 + t670 * t912) * qJD(3) + t200 * qJD(4) + t280 * qJD(5), t937 + (t562 * t680 - t670 * t914) * qJD(3) + t199 * qJD(4) + t281 * qJD(5), -t1096 + (qJD(5) / 0.2e1 - t781) * t915, -t259 * t926 + t33 * qJD(4) + t67 * qJD(5) + (pkin(10) * t562 - t564 * t985 + t1016) * t927 + t806, t260 * t926 + t34 * qJD(4) + t66 * qJD(5) + (pkin(10) * t563 - t565 * t985 + t1015) * t927 + t805 (-t493 * t563 - t515 * t562 + (-t184 * t679 - t208 * t675) * t676) * qJD(3) + t9 * qJD(4) + t15 * qJD(5) - t300 + t807 (t184 * t493 + t208 * t515 + t317 * t625) * qJD(3) + t6 * qJD(4) + t12 * qJD(5) + t40 * qJD(6) + t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t747, t783, t812 * t596, t812 * t598, t641, qJD(3) * t249 - qJD(4) * t391 - t789, qJD(3) * t247 + qJD(4) * t390 - t790, qJD(3) * t267 + t316 + (-t950 - t955) * t416, t148 * qJD(3) - t653 * t931 + t220 + t938, qJD(3) * t200 + t675 * t930 + t960, qJD(3) * t199 + t679 * t930 - t961, t710, t33 * qJD(3) + (-t391 * t679 + t675 * t825) * qJD(4) + t110 * qJD(5) + t804, t34 * qJD(3) + (t391 * t675 + t679 * t825) * qJD(4) + t108 * qJD(5) + t803, t9 * qJD(3) + (-t161 * t675 + t207 * t679 + (t637 * t679 - t638 * t675) * t596) * qJD(4) + t18 * qJD(5) + t318 + t808, t6 * qJD(3) + (t161 * t637 - t207 * t638 + t319 * t666) * qJD(4) + t14 * qJD(5) + t44 * qJD(6) + t809; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t712, t719, qJD(3) * t280 + t509 * t850, qJD(3) * t281 + t511 * t850, t503, qJD(3) * t67 + qJD(4) * t110 - qJD(5) * t193 - t798, qJD(3) * t66 + qJD(4) * t108 + qJD(5) * t192 - t799, pkin(5) * t945 + qJD(3) * t15 + qJD(4) * t18 + t1030, qJD(3) * t12 + qJD(4) * t14 - t1035 * t151 + t810; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t718, qJD(3) * t40 + qJD(4) * t44 + t802; 0, 0, 0, 0, 0, 0, 0, 0, 0, t271 * qJD(2), t272 * qJD(2), 0, 0, 0, 0, 0, -qJD(2) * t96 + qJD(4) * t309, -qJD(2) * t99, 0, 0, 0, 0, 0, -qJD(2) * t45 - qJD(4) * t69 - qJD(5) * t87, -qJD(2) * t48 - qJD(4) * t68 - qJD(5) * t88, qJD(2) * t29 + qJD(4) * t63, -qJD(2) * t2 + qJD(4) * t20 - qJD(5) * t27 + qJD(6) * t86 - t1010; 0, 0, 0, 0, -t845, -t929, -t681 * t828, t677 * t828, 0, t788, t787, -t598 * t843 + t1096, t342 - t935, -t680 * t903 + t933, t676 * t903 - t934, t583, qJD(4) * t248 - t797, qJD(4) * t246 - t796, qJD(4) * t268 - t563 * t955 + t296, qJD(4) * t149 - t209 - t939, -qJD(4) * t197 - qJD(5) * t279 - t936, -qJD(4) * t198 - qJD(5) * t278 - t937, -t1096 + (-t932 + t1037) * t915, -qJD(4) * t31 - qJD(5) * t64 - t806, -qJD(4) * t32 - qJD(5) * t65 - t805, qJD(4) * t10 + qJD(5) * t16 - t300 - t807, -qJD(4) * t5 - qJD(5) * t11 - qJD(6) * t39 - t811; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t895, t654 * qJD(4), 0, 0, 0, -pkin(3) * t949, -pkin(3) * t947, -t670 * t896 + t671 * t895, -qJD(5) * t619 - 0.2e1 * t680 * t847, -qJD(4) * t621 + t676 * t900, qJD(4) * t620 + t676 * t899, -t895, -qJD(4) * t363 - qJD(5) * t500, qJD(4) * t364 + qJD(5) * t499, -qJD(4) * t153 - qJD(5) * t155 + qJD(6) * t618, qJD(4) * t152 + qJD(5) * t157 - qJD(6) * t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t853, t782, -t812 * t680, t812 * t676, qJD(2) * t880, -pkin(10) * t947 - t715, pkin(10) * t949 - t768, t962 - t595 + (t671 * t927 + t902) * t680, t940 + t629 + (-0.2e1 * t844 + t1097) * t680, t675 * t949 + t786, t785 + t901, -t741 (t675 * t824 - t923) * qJD(4) + t462 * qJD(5) + t726 (t679 * t824 + t924) * qJD(4) + t461 * qJD(5) + t725 ((-t637 * t680 + t516) * t679 + (t638 * t680 - t502) * t675) * qJD(4) + t122 * qJD(5) + t729 (t502 * t637 - t516 * t638 + t626 * t666) * qJD(4) + t92 * qJD(5) + t195 * qJD(6) + t733; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t699, -t701, t675 * t793 - t958, t679 * t793 - t959, qJD(2) * t841 + t664, qJD(4) * t462 - qJD(5) * t565 - t723, qJD(4) * t461 + qJD(5) * t564 - t724, qJD(4) * t122 + t676 * t921 - t794, qJD(4) * t92 - t1035 * t515 + t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t784, qJD(4) * t195 + t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t146 - qJD(3) * t309, -qJD(2) * t144, 0, 0, 0, 0, 0, -qJD(2) * t54 + qJD(3) * t69, -qJD(2) * t57 + qJD(3) * t68, -qJD(2) * t36 - qJD(3) * t63, -qJD(2) * t4 - qJD(3) * t20 - qJD(6) * t79 - t1001; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t747, -t783, t781 * t985, -t780 * t985, t641, -qJD(3) * t248 + t789, -qJD(3) * t246 + t790, -qJD(3) * t268 + t679 * t908 + t316, -qJD(3) * t149 + t220 - t938, qJD(3) * t197 + qJD(5) * t416 - t960, qJD(3) * t198 - qJD(5) * t412 + t961, -t710, qJD(3) * t31 + qJD(5) * t109 - t804, qJD(3) * t32 + qJD(5) * t107 - t803, -qJD(3) * t10 - qJD(5) * t17 + t318 - t808, qJD(3) * t5 - qJD(5) * t13 - qJD(6) * t43 - t809; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t853, -t782, t843, -t676 * t906, qJD(2) * t879, t715, t768, -t657 * t671 - t595 - t962, t680 * t821 + t629 - t940, -t786 - t899, -t785 + t900, t741, qJD(5) * t508 - t726, qJD(5) * t507 - t725, qJD(5) * t121 - t729, -qJD(5) * t91 + qJD(6) * t196 - t733; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t896, t653 * qJD(5), 0, 0, 0, -pkin(4) * t943, -pkin(4) * t942, qJD(6) * t652, qJD(5) * t377 - qJD(6) * t801; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t711, -t700, -t679 * t849 + t956, t675 * t849 - t957, -t770, -pkin(11) * t942 - t713, pkin(11) * t943 - t714, -t795 - t920, t1035 * t638 + t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t779, t697; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t74 + qJD(3) * t87, qJD(2) * t75 + qJD(3) * t88, 0, -qJD(2) * t8 + qJD(3) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t712, -t719, qJD(3) * t279 - qJD(4) * t416 + t509 * t932, qJD(3) * t278 + qJD(4) * t412 + t908, t503, qJD(3) * t64 - qJD(4) * t109 + t798, qJD(3) * t65 - qJD(4) * t107 + t799, -qJD(3) * t16 + qJD(4) * t17 - t1030, qJD(3) * t11 + qJD(4) * t13 - t1034 * t511 - t810; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t699, t701, t958 + (-t675 * t927 + t948) * t680, t680 * t746 + t959, qJD(2) * t840 + t664, -qJD(4) * t508 + t723, -qJD(4) * t507 + t724, -qJD(4) * t121 + t794, qJD(4) * t91 - qJD(6) * t925 - t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t711, t700, t679 * t926 - t956, -t675 * t926 + t957, t770, t713, t714, t795, -t1034 * t675 - t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t72 - qJD(3) * t86 + qJD(4) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t718, qJD(3) * t39 + qJD(4) * t43 + t1035 * t511 - t802; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t784, -qJD(4) * t196 + t676 * t920 - t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t779, -t697 + t921; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t37;