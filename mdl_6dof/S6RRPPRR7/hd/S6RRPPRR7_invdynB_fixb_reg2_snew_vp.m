% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRPPRR7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:07:31
% EndTime: 2019-05-06 11:08:17
% DurationCPUTime: 30.08s
% Computational Cost: add. (86302->805), mult. (201340->1204), div. (0->0), fcn. (141027->10), ass. (0->580)
t900 = cos(pkin(6));
t889 = t900 * qJD(1) + qJD(2);
t1065 = t889 ^ 2;
t899 = sin(pkin(6));
t896 = t899 ^ 2;
t909 = qJD(1) ^ 2;
t1040 = t896 * t909;
t903 = sin(qJ(2));
t897 = t903 ^ 2;
t886 = t897 * t1040;
t1069 = t886 + t1065;
t907 = cos(qJ(2));
t1086 = t1040 * t907;
t864 = t903 * t1086;
t888 = qJDD(1) * t900 + qJDD(2);
t832 = -t864 + t888;
t816 = t907 * t832;
t743 = t1069 * t903 - t816;
t1103 = pkin(8) * t743;
t898 = t907 ^ 2;
t887 = t898 * t1040;
t805 = t887 + t1065;
t831 = t864 + t888;
t744 = t805 * t907 + t831 * t903;
t1094 = pkin(8) * t744;
t904 = sin(qJ(1));
t1100 = t743 * t904;
t908 = cos(qJ(1));
t1099 = t743 * t908;
t1090 = t744 * t904;
t1089 = t744 * t908;
t1044 = t832 * t903;
t1068 = t887 - t1065;
t757 = t1068 * t907 - t1044;
t1088 = t757 * t904;
t1087 = t757 * t908;
t1070 = -t886 + t1065;
t815 = t907 * t831;
t1079 = -t1070 * t903 + t815;
t1102 = t1079 * t904;
t1101 = t1079 * t908;
t1027 = t900 * t903;
t1098 = -t1027 * t805 + t900 * t815;
t1033 = t899 * t903;
t1097 = -t1033 * t805 + t899 * t815;
t843 = t886 - t887;
t1035 = t899 * t843;
t1008 = qJDD(1) * t903;
t1012 = qJD(1) * t907;
t777 = (t1008 + (qJD(2) + t889) * t1012) * t899;
t1049 = t777 * t907;
t1007 = qJDD(1) * t907;
t1013 = qJD(1) * t903;
t996 = t899 * t1013;
t875 = qJD(2) * t996;
t834 = t1007 * t899 - t875;
t849 = t889 * t996;
t776 = t834 - t849;
t954 = t776 * t903 + t1049;
t675 = t900 * t954 - t1035;
t1050 = t777 * t903;
t716 = t776 * t907 - t1050;
t1096 = t675 * t904 - t716 * t908;
t1095 = t675 * t908 + t716 * t904;
t833 = (qJD(2) * t1012 + t1008) * t899;
t995 = t899 * t1012;
t848 = t889 * t995;
t779 = t848 + t833;
t1093 = qJ(3) * t779;
t902 = sin(qJ(5));
t906 = cos(qJ(5));
t811 = -t889 * t906 + t902 * t995;
t812 = t902 * t889 + t906 * t995;
t751 = t811 * t812;
t824 = qJDD(5) + t833;
t1071 = -t751 + t824;
t1085 = t1071 * t902;
t1084 = t1071 * t906;
t867 = qJD(5) + t996;
t901 = sin(qJ(6));
t905 = cos(qJ(6));
t763 = -t812 * t901 - t905 * t867;
t765 = -t812 * t905 + t867 * t901;
t711 = t765 * t763;
t988 = t902 * t834 - t906 * t888;
t731 = qJD(5) * t812 + t988;
t730 = -qJDD(6) + t731;
t1072 = -t711 - t730;
t1083 = t1072 * t901;
t1082 = t1072 * t905;
t1046 = t811 * t867;
t999 = -t811 * qJD(5) + t906 * t834 + t902 * t888;
t1081 = t999 - t1046;
t1026 = t900 * t907;
t1078 = t1026 * t1069 + t832 * t1027;
t1077 = t1026 * t1070 + t831 * t1027;
t1076 = t1027 * t1068 + t900 * t816;
t1032 = t899 * t907;
t1075 = t1032 * t1069 + t832 * t1033;
t1074 = t1032 * t1070 + t831 * t1033;
t1073 = t1033 * t1068 + t899 * t816;
t651 = -t763 * qJD(6) + t901 * t824 - t905 * t999;
t806 = -qJD(6) + t811;
t728 = t763 * t806;
t605 = t728 + t651;
t820 = t900 * t843;
t1067 = t899 * t954 + t820;
t1041 = t888 * t899;
t1014 = qJD(1) * t899;
t998 = t889 * t1014;
t948 = (t897 + t898) * t998;
t747 = t904 * t1041 + t908 * t948;
t1058 = t907 * g(3);
t982 = -pkin(2) * t907 - qJ(3) * t903;
t830 = t982 * t1014;
t1061 = pkin(8) * t899;
t873 = g(1) * t904 - t908 * g(2);
t828 = qJDD(1) * pkin(1) + t1061 * t909 + t873;
t874 = g(1) * t908 + g(2) * t904;
t829 = -pkin(1) * t909 + qJDD(1) * t1061 - t874;
t989 = -t828 * t1026 + t903 * t829;
t938 = -t888 * pkin(2) - qJ(3) * t1065 + qJDD(3) + t989;
t666 = -t899 * (t1013 * t830 + t1058) - t938;
t746 = -t908 * t1041 + t904 * t948;
t990 = -t905 * t824 - t901 * t999;
t602 = (qJD(6) + t806) * t765 + t990;
t985 = qJD(3) * t996;
t870 = 0.2e1 * t985;
t827 = -pkin(3) * t889 - qJ(4) * t996;
t1059 = t900 * g(3);
t941 = pkin(2) * t776 + t1059 + t1093;
t914 = t834 * pkin(3) - qJ(4) * t887 + t827 * t996 + qJDD(4) + t941;
t571 = t833 * pkin(4) + t834 * pkin(9) + t870 + (t828 + (pkin(4) * t907 - pkin(9) * t903) * t889 * qJD(1)) * t899 + t914;
t1064 = -2 * qJD(4);
t915 = -pkin(3) * t831 - t833 * qJ(4) + t938;
t1000 = t889 * t907 * qJ(4);
t944 = -t903 * t830 - t1000;
t984 = pkin(4) * t903 + pkin(9) * t907;
t586 = -t1065 * pkin(4) - t888 * pkin(9) + (t1058 + (t1064 * t903 - t984 * t996 - t944) * qJD(1)) * t899 + t915;
t520 = t902 * t571 + t906 * t586;
t991 = -t906 * t571 + t902 * t586;
t463 = t906 * t520 + t902 * t991;
t462 = t902 * t520 - t906 * t991;
t924 = t1027 * t848 + t899 * t864;
t919 = t834 * t1026 - t924;
t937 = -t834 * t903 - t898 * t998;
t667 = t904 * t937 + t908 * t919;
t669 = -t904 * t919 + t908 * t937;
t920 = t833 * t1027 + t924;
t947 = t907 * t833 - t897 * t998;
t668 = t904 * t947 + t908 * t920;
t670 = -t904 * t920 + t908 * t947;
t761 = t763 ^ 2;
t762 = t765 ^ 2;
t803 = t806 ^ 2;
t810 = t811 ^ 2;
t1066 = t812 ^ 2;
t860 = t867 ^ 2;
t1063 = pkin(2) + pkin(3);
t844 = -t886 - t887;
t1034 = t899 * t844;
t775 = t834 + t849;
t778 = (t1008 + (qJD(2) - t889) * t1012) * t899;
t955 = t775 * t903 - t778 * t907;
t674 = t900 * t955 - t1034;
t715 = t775 * t907 + t778 * t903;
t612 = t674 * t908 + t715 * t904;
t1062 = pkin(7) * t612;
t1060 = pkin(8) * t900;
t1057 = pkin(4) + qJ(3);
t1056 = qJ(3) * t805;
t1055 = qJ(3) * t844;
t635 = t711 - t730;
t1054 = t635 * t901;
t1053 = t635 * t905;
t997 = t889 * t1013;
t772 = -t875 + (-t997 + t1007) * t899;
t1052 = t772 * t903;
t1051 = t772 * t907;
t1048 = t806 * t901;
t1047 = t806 * t905;
t1043 = t867 * t902;
t1042 = t867 * t906;
t1039 = t899 * t775;
t1038 = t899 * t776;
t1037 = t899 * t779;
t1036 = t899 * t828;
t1031 = t900 * t775;
t1030 = t900 * t776;
t1029 = t900 * t779;
t1028 = t900 * t844;
t748 = -pkin(5) * t811 + pkin(10) * t812;
t496 = -t824 * pkin(5) - t860 * pkin(10) - t748 * t812 + t991;
t1025 = t901 * t496;
t1011 = qJD(3) * t889;
t865 = 0.2e1 * t1011;
t742 = -g(3) * t1033 + t828 * t1027 + t907 * t829;
t946 = pkin(2) * t1065 - t888 * qJ(3) - t830 * t995 - t742;
t987 = 0.2e1 * qJD(4) * t1014;
t911 = pkin(3) * t887 - t889 * t827 + t907 * t987 + t946;
t617 = -t834 * qJ(4) + t865 - t911;
t581 = t888 * pkin(4) - pkin(9) * t1065 - t1086 * t984 + t617;
t1023 = t902 * t581;
t727 = t751 + t824;
t1022 = t902 * t727;
t782 = t1036 + t1059;
t1021 = t903 * t782;
t1020 = t905 * t496;
t1018 = t906 * t581;
t1017 = t906 * t727;
t1016 = t907 * t782;
t497 = -pkin(5) * t860 + pkin(10) * t824 + t748 * t811 + t520;
t532 = t581 + t1081 * pkin(10) + (-t812 * t867 - t731) * pkin(5);
t461 = t905 * t497 + t901 * t532;
t1015 = pkin(1) * t674 + t715 * t1061;
t1009 = pkin(9) + t1063;
t1006 = t761 + t762;
t1005 = t902 * t711;
t1004 = t906 * t711;
t1003 = t903 * t751;
t1002 = t907 * t751;
t994 = -t860 - t1066;
t671 = t899 * t955 + t1028;
t993 = -pkin(1) * t671 + t715 * t1060;
t992 = -pkin(5) * t902 - qJ(4);
t460 = t497 * t901 - t905 * t532;
t438 = t901 * t460 + t905 * t461;
t786 = -t873 * t904 - t908 * t874;
t986 = pkin(5) * t906 + t1057;
t862 = qJDD(1) * t908 - t904 * t909;
t983 = -pkin(7) * t862 - g(3) * t904;
t981 = -t1021 + t1103;
t980 = t1016 - t1094;
t426 = t438 * t906 + t902 * t496;
t437 = -t460 * t905 + t461 * t901;
t979 = -t426 * t907 + t437 * t903;
t978 = -t463 * t907 + t581 * t903;
t606 = t728 - t651;
t536 = -t602 * t905 - t606 * t901;
t507 = -t1006 * t902 + t906 * t536;
t534 = -t602 * t901 + t606 * t905;
t977 = -t507 * t907 + t534 * t903;
t603 = (-qJD(6) + t806) * t765 - t990;
t535 = t603 * t905 - t605 * t901;
t710 = -t762 + t761;
t515 = -t535 * t906 + t710 * t902;
t533 = t603 * t901 + t605 * t905;
t976 = t515 * t907 + t533 * t903;
t664 = -t803 - t761;
t565 = t664 * t905 - t1083;
t522 = t565 * t906 - t902 * t603;
t564 = t664 * t901 + t1082;
t975 = -t522 * t907 + t564 * t903;
t697 = -t762 - t803;
t573 = -t697 * t901 - t1053;
t525 = t906 * t573 + t605 * t902;
t572 = t697 * t905 - t1054;
t974 = -t525 * t907 + t572 * t903;
t724 = t761 - t803;
t591 = t724 * t905 - t1054;
t526 = -t591 * t906 + t602 * t902;
t589 = t724 * t901 + t1053;
t973 = t526 * t907 + t589 * t903;
t725 = -t762 + t803;
t590 = -t725 * t901 + t1082;
t527 = -t590 * t906 + t606 * t902;
t588 = t725 * t905 + t1083;
t972 = t527 * t907 + t588 * t903;
t599 = t1048 * t765 + t651 * t905;
t553 = -t599 * t906 - t1005;
t598 = -t1047 * t765 + t651 * t901;
t971 = t553 * t907 + t598 * t903;
t650 = -qJD(6) * t765 - t990;
t597 = -t1047 * t763 - t650 * t901;
t554 = -t597 * t906 + t1005;
t596 = -t1048 * t763 + t650 * t905;
t970 = t554 * t907 + t596 * t903;
t646 = (t763 * t905 - t765 * t901) * t806;
t594 = -t646 * t906 + t730 * t902;
t645 = (t763 * t901 + t765 * t905) * t806;
t969 = t594 * t907 + t645 * t903;
t691 = (-qJD(5) - t867) * t812 - t988;
t609 = -t1081 * t902 + t691 * t906;
t750 = -t810 + t1066;
t968 = t609 * t907 + t750 * t903;
t695 = t999 + t1046;
t923 = (qJD(5) - t867) * t812 + t988;
t610 = -t902 * t695 + t906 * t923;
t722 = -t810 - t1066;
t967 = -t610 * t907 + t722 * t903;
t618 = (t1058 + (t1000 + (t1064 + t830) * t903) * qJD(1)) * t899 + t915;
t966 = t617 * t903 - t618 * t907;
t739 = -t860 - t810;
t648 = t906 * t739 - t1085;
t965 = -t648 * t907 + t691 * t903;
t657 = t865 - t946;
t649 = -pkin(2) * t844 + t657;
t653 = -t1055 - t666;
t964 = t649 * t907 + t653 * t903;
t656 = -t902 * t994 - t1017;
t963 = -t1081 * t903 - t656 * t907;
t962 = t657 * t903 + t666 * t907;
t770 = t810 - t860;
t660 = -t770 * t906 + t1022;
t961 = t660 * t907 + t903 * t923;
t771 = t860 - t1066;
t662 = t771 * t902 - t1084;
t960 = t662 * t907 - t695 * t903;
t720 = (-t811 * t906 + t812 * t902) * t867;
t959 = t720 * t907 + t824 * t903;
t741 = g(3) * t1032 + t989;
t958 = -t907 * t741 + t903 * t742;
t654 = t741 * t903 + t742 * t907;
t957 = t779 * t907 + t1052;
t773 = -t875 + (t997 + t1007) * t899;
t774 = -t848 + t833;
t956 = -t773 * t903 + t774 * t907;
t953 = -t1069 * t907 - t1044;
t952 = t805 * t903 - t815;
t951 = -t1068 * t903 - t816;
t785 = t873 * t908 - t874 * t904;
t680 = -t1043 * t812 + t906 * t999;
t943 = t680 * t907 + t1003;
t681 = t1042 * t811 + t731 * t902;
t942 = t681 * t907 - t1003;
t940 = (-t671 * t899 - t674 * t900) * pkin(8);
t939 = t896 * t907 * t997 - t900 * t864;
t425 = t438 * t902 - t496 * t906;
t395 = (pkin(10) * t906 + t992) * t437 + t1009 * t425;
t402 = -pkin(5) * t496 + pkin(10) * t438 - qJ(4) * t426 + t1057 * t425;
t411 = t426 * t903 + t437 * t907;
t936 = pkin(8) * t411 + t395 * t907 + t402 * t903;
t430 = -pkin(10) * t534 - t437;
t506 = t1006 * t906 + t902 * t536;
t409 = t1009 * t506 - t906 * t430 + t534 * t992;
t412 = pkin(5) * t1006 + pkin(10) * t536 - qJ(4) * t507 + t1057 * t506 + t438;
t465 = t507 * t903 + t534 * t907;
t935 = pkin(8) * t465 + t409 * t907 + t412 * t903;
t449 = -pkin(5) * t564 + t460;
t470 = -pkin(10) * t564 + t1025;
t521 = t565 * t902 + t603 * t906;
t415 = -qJ(4) * t564 + t1009 * t521 + t902 * t449 - t906 * t470;
t427 = pkin(5) * t603 + pkin(10) * t565 - qJ(4) * t522 + t1057 * t521 - t1020;
t481 = t522 * t903 + t564 * t907;
t934 = pkin(8) * t481 + t415 * t907 + t427 * t903;
t450 = -pkin(5) * t572 + t461;
t472 = -pkin(10) * t572 + t1020;
t524 = t902 * t573 - t605 * t906;
t416 = -qJ(4) * t572 + t1009 * t524 + t902 * t450 - t906 * t472;
t431 = -pkin(5) * t605 + pkin(10) * t573 - qJ(4) * t525 + t1057 * t524 + t1025;
t485 = t525 * t903 + t572 * t907;
t933 = pkin(8) * t485 + t416 * t907 + t431 * t903;
t421 = -qJ(4) * t581 + t1009 * t462;
t424 = -qJ(4) * t463 + t1057 * t462;
t453 = t463 * t903 + t581 * t907;
t932 = pkin(8) * t453 + t421 * t907 + t424 * t903;
t607 = t695 * t906 + t902 * t923;
t445 = -qJ(4) * t722 + t1009 * t607 + t462;
t491 = -qJ(4) * t610 + t1057 * t607;
t560 = t610 * t903 + t722 * t907;
t931 = pkin(8) * t560 + t445 * t907 + t491 * t903;
t647 = t902 * t739 + t1084;
t469 = -qJ(4) * t648 + t1057 * t647 - t991;
t482 = -qJ(4) * t691 + t1009 * t647 - t1023;
t578 = t648 * t903 + t691 * t907;
t930 = pkin(8) * t578 + t469 * t903 + t482 * t907;
t655 = t906 * t994 - t1022;
t471 = -qJ(4) * t656 + t1057 * t655 - t520;
t486 = qJ(4) * t1081 + t1009 * t655 - t1018;
t585 = -t1081 * t907 + t656 * t903;
t929 = pkin(8) * t585 + t471 * t903 + t486 * t907;
t871 = -0.2e1 * t985;
t912 = -t914 - t1036;
t619 = t871 + t912;
t508 = -qJ(4) * t617 - t1063 * t619;
t541 = t617 * t907 + t618 * t903;
t545 = -qJ(3) * t619 - qJ(4) * t618;
t928 = pkin(8) * t541 + t508 * t907 + t545 * t903;
t561 = -qJ(4) * t805 - t1063 * t776 + t619;
t740 = -qJ(3) * t776 - qJ(4) * t831;
t927 = t561 * t907 + t740 * t903 + t1094;
t566 = -0.2e1 * t1011 + t1063 * t844 + (t773 + t834) * qJ(4) + t911;
t579 = t1055 + qJ(4) * t774 + t903 * t987 + (qJD(1) * t944 - t1058) * t899 - t915;
t717 = -t773 * t907 - t774 * t903;
t926 = pkin(8) * t717 + t566 * t907 + t579 * t903;
t580 = qJ(4) * t1069 + t1093 + t870 - t912;
t690 = -qJ(4) * t832 + t1063 * t779;
t925 = t580 * t903 + t690 * t907 - t1103;
t736 = t833 * t1033 + t939;
t735 = t834 * t1032 - t939;
t918 = t941 + t1036;
t916 = t870 + t918;
t632 = qJ(3) * t777 + t916;
t922 = pkin(2) * t1049 + t632 * t903 - t1103;
t633 = pkin(2) * t772 + t916;
t921 = qJ(3) * t1052 + t633 * t907 - t1094;
t582 = t657 * t907 - t666 * t903;
t658 = t871 - t918;
t917 = pkin(8) * t582 + t658 * t982;
t861 = qJDD(1) * t904 + t908 * t909;
t856 = t900 * t888;
t835 = -pkin(7) * t861 + g(3) * t908;
t821 = qJ(3) * t832;
t721 = -pkin(2) * t778 + qJ(3) * t775;
t719 = (-t811 * t902 - t812 * t906) * t867;
t714 = -t779 * t903 + t1051;
t709 = t1038 + t1098;
t708 = t900 * t951 + t1039;
t707 = -t773 * t899 + t1076;
t706 = -t1039 + t1076;
t705 = t1037 + t1078;
t704 = t777 * t899 + t1078;
t703 = -t774 * t899 + t1077;
t702 = -t778 * t899 + t1077;
t701 = -t1030 + t1097;
t700 = -t1029 + t1075;
t699 = -t777 * t900 + t1075;
t698 = t778 * t900 + t1074;
t689 = t772 * t899 + t1098;
t688 = t900 * t952 - t1038;
t687 = t900 * t953 - t1037;
t686 = -t772 * t900 + t1097;
t685 = t899 * t952 + t1030;
t684 = t899 * t953 + t1029;
t682 = t1042 * t812 + t902 * t999;
t679 = t1043 * t811 - t731 * t906;
t678 = -t720 * t903 + t824 * t907;
t676 = t900 * t956 + t1034;
t673 = t900 * t957 - t1035;
t672 = t899 * t956 - t1028;
t663 = -qJ(3) * t773 + t1063 * t774;
t661 = -t770 * t902 - t1017;
t659 = -t771 * t906 - t1085;
t644 = -t709 * t904 - t1089;
t643 = -t705 * t904 - t1099;
t642 = -t704 * t904 - t1099;
t641 = -t702 * t904 + t1101;
t640 = t709 * t908 - t1090;
t639 = t705 * t908 - t1100;
t638 = t704 * t908 - t1100;
t637 = t702 * t908 + t1102;
t631 = pkin(2) * t831 - t1056 + t666;
t630 = pkin(2) * t1069 + t657 + t821;
t629 = -t680 * t903 + t1002;
t628 = -t681 * t903 - t1002;
t627 = -t689 * t904 - t1089;
t626 = -t688 * t904 + t1089;
t625 = -t687 * t904 + t1099;
t624 = t689 * t908 - t1090;
t623 = t688 * t908 + t1090;
t622 = t687 * t908 + t1100;
t621 = t899 * t782 + t900 * t958;
t620 = -t900 * t782 + t899 * t958;
t616 = -t899 * t719 + t900 * t959;
t615 = -t676 * t904 + t717 * t908;
t614 = -t674 * t904 + t715 * t908;
t613 = t676 * t908 + t717 * t904;
t611 = pkin(7) * t614;
t608 = t1081 * t906 + t691 * t902;
t595 = -t646 * t902 - t730 * t906;
t593 = -t660 * t903 + t907 * t923;
t592 = -t662 * t903 - t695 * t907;
t587 = -t1021 + (-t701 * t899 - t709 * t900) * pkin(8);
t584 = pkin(2) * t666 + qJ(3) * t657;
t577 = -t1016 + (-t684 * t899 - t687 * t900) * pkin(8);
t576 = -pkin(1) * t701 + t899 * t741 + t900 * t980;
t575 = t1063 * t1069 + t617 + t821;
t574 = -t1063 * t831 + t1056 + t618;
t568 = -t609 * t903 + t750 * t907;
t567 = -pkin(1) * t684 + t899 * t742 + t900 * t981;
t563 = -t899 * t682 + t900 * t943;
t562 = -t899 * t679 + t900 * t942;
t559 = -pkin(1) * t620 + t1060 * t654;
t558 = -t621 * t904 + t654 * t908;
t557 = t621 * t908 + t654 * t904;
t556 = -t599 * t902 + t1004;
t555 = -t597 * t902 - t1004;
t552 = t654 * t900 + t993;
t551 = -t899 * t661 + t900 * t961;
t550 = -t899 * t659 + t900 * t960;
t549 = (-t620 * t899 - t621 * t900) * pkin(8);
t548 = t940 - t958;
t547 = -t899 * t658 + t900 * t962;
t546 = t900 * t658 + t899 * t962;
t544 = t899 * t655 + t900 * t963;
t543 = -t900 * t655 + t899 * t963;
t542 = -t594 * t903 + t645 * t907;
t540 = t899 * t647 + t900 * t965;
t539 = -t900 * t647 + t899 * t965;
t538 = -pkin(2) * t1050 + t907 * t632 + (-t699 * t899 - t704 * t900) * pkin(8);
t537 = qJ(3) * t1051 - t903 * t633 + (-t686 * t899 - t689 * t900) * pkin(8);
t529 = -t591 * t902 - t602 * t906;
t528 = -t590 * t902 - t606 * t906;
t523 = -t899 * t608 + t900 * t968;
t519 = t899 * t607 + t900 * t967;
t518 = -t900 * t607 + t899 * t967;
t516 = -t903 * t649 + t907 * t653 + t940;
t514 = -t535 * t902 - t710 * t906;
t513 = -pkin(1) * t699 - t899 * t630 + t900 * t922;
t512 = -t553 * t903 + t598 * t907;
t511 = -t554 * t903 + t596 * t907;
t510 = -pkin(1) * t686 - t899 * t631 + t900 * t921;
t509 = t907 * t580 - t903 * t690 + (-t700 * t899 - t705 * t900) * pkin(8);
t505 = -t899 * t721 + t900 * t964 + t993;
t504 = qJ(3) * t617 - t1063 * t618;
t503 = -t903 * t561 + t907 * t740 + (-t685 * t899 - t688 * t900) * pkin(8);
t502 = -t899 * t619 + t900 * t966;
t501 = t900 * t619 + t899 * t966;
t500 = -t899 * t595 + t900 * t969;
t499 = -t544 * t904 + t585 * t908;
t498 = t544 * t908 + t585 * t904;
t495 = -t547 * t904 + t582 * t908;
t494 = t547 * t908 + t582 * t904;
t490 = -t540 * t904 + t578 * t908;
t489 = t540 * t908 + t578 * t904;
t488 = -t526 * t903 + t589 * t907;
t487 = -t527 * t903 + t588 * t907;
t484 = -t903 * t566 + t907 * t579 + (-t672 * t899 - t676 * t900) * pkin(8);
t483 = -pkin(1) * t700 - t899 * t575 + t900 * t925;
t480 = -pkin(1) * t685 - t899 * t574 + t900 * t927;
t479 = -t519 * t904 + t560 * t908;
t478 = t519 * t908 + t560 * t904;
t477 = -t1009 * t656 - t1057 * t1081 - t1023;
t476 = -pkin(1) * t672 - t899 * t663 + t900 * t926;
t475 = -t899 * t556 + t900 * t971;
t474 = -t899 * t555 + t900 * t970;
t473 = -t1009 * t648 + t1057 * t691 + t1018;
t468 = -t515 * t903 + t533 * t907;
t467 = -t502 * t904 + t541 * t908;
t466 = t502 * t908 + t541 * t904;
t464 = (pkin(2) * t903 - qJ(3) * t907) * t658 + (-t546 * t899 - t547 * t900) * pkin(8);
t458 = -pkin(1) * t546 - t899 * t584 + t900 * t917;
t457 = -t899 * t529 + t900 * t973;
t456 = -t899 * t528 + t900 * t972;
t455 = t899 * t524 + t900 * t974;
t454 = -t900 * t524 + t899 * t974;
t452 = t899 * t521 + t900 * t975;
t451 = -t900 * t521 + t899 * t975;
t448 = -t899 * t514 + t900 * t976;
t447 = t899 * t506 + t900 * t977;
t446 = -t900 * t506 + t899 * t977;
t444 = -t1009 * t610 + t1057 * t722 - t463;
t443 = -t455 * t904 + t485 * t908;
t442 = t455 * t908 + t485 * t904;
t441 = -t452 * t904 + t481 * t908;
t440 = t452 * t908 + t481 * t904;
t439 = -t903 * t508 + t907 * t545 + (-t501 * t899 - t502 * t900) * pkin(8);
t436 = t899 * t462 + t900 * t978;
t435 = -t900 * t462 + t899 * t978;
t434 = t907 * t471 - t903 * t486 + (-t543 * t899 - t544 * t900) * pkin(8);
t433 = -pkin(1) * t501 - t899 * t504 + t900 * t928;
t432 = t907 * t469 - t903 * t482 + (-t539 * t899 - t540 * t900) * pkin(8);
t429 = -t447 * t904 + t465 * t908;
t428 = t447 * t908 + t465 * t904;
t423 = -pkin(1) * t543 - t899 * t477 + t900 * t929;
t422 = -pkin(1) * t539 - t899 * t473 + t900 * t930;
t420 = -t903 * t445 + t907 * t491 + (-t518 * t899 - t519 * t900) * pkin(8);
t419 = -t1009 * t463 + t1057 * t581;
t418 = -t436 * t904 + t453 * t908;
t417 = t436 * t908 + t453 * t904;
t414 = -t1009 * t525 + t1057 * t572 - t906 * t450 - t902 * t472;
t413 = -t1009 * t522 + t1057 * t564 - t906 * t449 - t902 * t470;
t410 = -pkin(1) * t518 - t899 * t444 + t900 * t931;
t408 = -t1009 * t507 - t902 * t430 + t534 * t986;
t407 = t899 * t425 + t900 * t979;
t406 = -t900 * t425 + t899 * t979;
t405 = -t903 * t416 + t907 * t431 + (-t454 * t899 - t455 * t900) * pkin(8);
t404 = -t903 * t415 + t907 * t427 + (-t451 * t899 - t452 * t900) * pkin(8);
t403 = -t903 * t421 + t907 * t424 + (-t435 * t899 - t436 * t900) * pkin(8);
t401 = -pkin(1) * t454 - t899 * t414 + t900 * t933;
t400 = -t407 * t904 + t411 * t908;
t399 = t407 * t908 + t411 * t904;
t398 = -t903 * t409 + t907 * t412 + (-t446 * t899 - t447 * t900) * pkin(8);
t397 = -pkin(1) * t451 - t899 * t413 + t900 * t934;
t396 = -pkin(1) * t435 - t899 * t419 + t900 * t932;
t394 = -t1009 * t426 + (pkin(10) * t902 + t986) * t437;
t393 = -pkin(1) * t446 - t899 * t408 + t900 * t935;
t392 = -t903 * t395 + t907 * t402 + (-t406 * t899 - t407 * t900) * pkin(8);
t391 = -pkin(1) * t406 - t899 * t394 + t900 * t936;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t861, -t862, 0, t786, 0, 0, 0, 0, 0, 0, t644, t625, t614, t558, 0, 0, 0, 0, 0, 0, t627, t614, t642, t495, 0, 0, 0, 0, 0, 0, t643, t626, t615, t467, 0, 0, 0, 0, 0, 0, t490, t499, t479, t418, 0, 0, 0, 0, 0, 0, t441, t443, t429, t400; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t862, -t861, 0, t785, 0, 0, 0, 0, 0, 0, t640, t622, t612, t557, 0, 0, 0, 0, 0, 0, t624, t612, t638, t494, 0, 0, 0, 0, 0, 0, t639, t623, t613, t466, 0, 0, 0, 0, 0, 0, t489, t498, t478, t417, 0, 0, 0, 0, 0, 0, t440, t442, t428, t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t701, t684, t671, t620, 0, 0, 0, 0, 0, 0, t686, t671, t699, t546, 0, 0, 0, 0, 0, 0, t700, t685, t672, t501, 0, 0, 0, 0, 0, 0, t539, t543, t518, t435, 0, 0, 0, 0, 0, 0, t451, t454, t446, t406; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t862, 0, -t861, 0, t983, -t835, -t785, -pkin(7) * t785, t670, -t1096, t641, t669, -t706 * t904 + t1087, t747, -pkin(7) * t640 - t576 * t904 + t587 * t908, -pkin(7) * t622 - t567 * t904 + t577 * t908, t548 * t908 - t552 * t904 - t1062, -pkin(7) * t557 + t549 * t908 - t559 * t904, t670, t641, t1096, t747, -t708 * t904 - t1087, t669, -pkin(7) * t624 - t510 * t904 + t537 * t908, -t505 * t904 + t516 * t908 - t1062, -pkin(7) * t638 - t513 * t904 + t538 * t908, -pkin(7) * t494 - t458 * t904 + t464 * t908, t669, -t673 * t904 + t714 * t908, -t707 * t904 + t1087, t670, -t703 * t904 + t1101, t747, -pkin(7) * t639 - t483 * t904 + t509 * t908, -pkin(7) * t623 - t480 * t904 + t503 * t908, -pkin(7) * t613 - t476 * t904 + t484 * t908, -pkin(7) * t466 - t433 * t904 + t439 * t908, -t563 * t904 + t629 * t908, -t523 * t904 + t568 * t908, -t550 * t904 + t592 * t908, -t562 * t904 + t628 * t908, -t551 * t904 + t593 * t908, -t616 * t904 + t678 * t908, -pkin(7) * t489 - t422 * t904 + t432 * t908, -pkin(7) * t498 - t423 * t904 + t434 * t908, -pkin(7) * t478 - t410 * t904 + t420 * t908, -pkin(7) * t417 - t396 * t904 + t403 * t908, -t475 * t904 + t512 * t908, -t448 * t904 + t468 * t908, -t456 * t904 + t487 * t908, -t474 * t904 + t511 * t908, -t457 * t904 + t488 * t908, -t500 * t904 + t542 * t908, -pkin(7) * t440 - t397 * t904 + t404 * t908, -pkin(7) * t442 - t401 * t904 + t405 * t908, -pkin(7) * t428 - t393 * t904 + t398 * t908, -pkin(7) * t399 - t391 * t904 + t392 * t908; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t861, 0, t862, 0, t835, t983, t786, pkin(7) * t786, t668, t1095, t637, t667, t706 * t908 + t1088, t746, pkin(7) * t644 + t576 * t908 + t587 * t904, pkin(7) * t625 + t567 * t908 + t577 * t904, t548 * t904 + t552 * t908 + t611, pkin(7) * t558 + t549 * t904 + t559 * t908, t668, t637, -t1095, t746, t708 * t908 - t1088, t667, pkin(7) * t627 + t510 * t908 + t537 * t904, t505 * t908 + t516 * t904 + t611, pkin(7) * t642 + t513 * t908 + t538 * t904, pkin(7) * t495 + t458 * t908 + t464 * t904, t667, t673 * t908 + t714 * t904, t707 * t908 + t1088, t668, t703 * t908 + t1102, t746, pkin(7) * t643 + t483 * t908 + t509 * t904, pkin(7) * t626 + t480 * t908 + t503 * t904, pkin(7) * t615 + t476 * t908 + t484 * t904, pkin(7) * t467 + t433 * t908 + t439 * t904, t563 * t908 + t629 * t904, t523 * t908 + t568 * t904, t550 * t908 + t592 * t904, t562 * t908 + t628 * t904, t551 * t908 + t593 * t904, t616 * t908 + t678 * t904, pkin(7) * t490 + t422 * t908 + t432 * t904, pkin(7) * t499 + t423 * t908 + t434 * t904, pkin(7) * t479 + t410 * t908 + t420 * t904, pkin(7) * t418 + t396 * t908 + t403 * t904, t475 * t908 + t512 * t904, t448 * t908 + t468 * t904, t456 * t908 + t487 * t904, t474 * t908 + t511 * t904, t457 * t908 + t488 * t904, t500 * t908 + t542 * t904, pkin(7) * t441 + t397 * t908 + t404 * t904, pkin(7) * t443 + t401 * t908 + t405 * t904, pkin(7) * t429 + t393 * t908 + t398 * t904, pkin(7) * t400 + t391 * t908 + t392 * t904; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t873, t874, 0, 0, t736, t1067, t698, t735, t1031 + t1073, t856, pkin(1) * t709 - t900 * t741 + t899 * t980, pkin(1) * t687 - t900 * t742 + t899 * t981, t654 * t899 + t1015, pkin(1) * t621 + t1061 * t654, t736, t698, -t1067, t856, t899 * t951 - t1031, t735, pkin(1) * t689 + t900 * t631 + t899 * t921, t900 * t721 + t899 * t964 + t1015, pkin(1) * t704 + t900 * t630 + t899 * t922, pkin(1) * t547 + t900 * t584 + t899 * t917, t735, t899 * t957 + t820, t773 * t900 + t1073, t736, t774 * t900 + t1074, t856, pkin(1) * t705 + t900 * t575 + t899 * t925, pkin(1) * t688 + t900 * t574 + t899 * t927, pkin(1) * t676 + t900 * t663 + t899 * t926, pkin(1) * t502 + t900 * t504 + t899 * t928, t900 * t682 + t899 * t943, t900 * t608 + t899 * t968, t900 * t659 + t899 * t960, t900 * t679 + t899 * t942, t900 * t661 + t899 * t961, t900 * t719 + t899 * t959, pkin(1) * t540 + t900 * t473 + t899 * t930, pkin(1) * t544 + t900 * t477 + t899 * t929, pkin(1) * t519 + t900 * t444 + t899 * t931, pkin(1) * t436 + t900 * t419 + t899 * t932, t900 * t556 + t899 * t971, t900 * t514 + t899 * t976, t900 * t528 + t899 * t972, t900 * t555 + t899 * t970, t900 * t529 + t899 * t973, t900 * t595 + t899 * t969, pkin(1) * t452 + t900 * t413 + t899 * t934, pkin(1) * t455 + t900 * t414 + t899 * t933, pkin(1) * t447 + t900 * t408 + t899 * t935, pkin(1) * t407 + t900 * t394 + t899 * t936;];
tauB_reg  = t1;