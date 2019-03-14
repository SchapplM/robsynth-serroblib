% Calculate minimal parameter regressor of coriolis matrix for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRRRPP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:39
% EndTime: 2019-03-09 21:03:32
% DurationCPUTime: 37.18s
% Computational Cost: add. (32104->990), mult. (66161->1265), div. (0->0), fcn. (71453->8), ass. (0->722)
t1065 = cos(qJ(4));
t662 = sin(qJ(3));
t882 = t1065 * t662;
t661 = sin(qJ(4));
t664 = cos(qJ(3));
t960 = t661 * t664;
t746 = t882 + t960;
t580 = t746 * qJ(5);
t1098 = pkin(8) + pkin(9);
t825 = t1098 * t1065;
t803 = t662 * t825;
t827 = t1098 * t960;
t693 = -t803 - t827;
t1125 = -t580 + t693;
t660 = sin(pkin(10));
t1133 = t660 * t1125;
t1050 = cos(pkin(10));
t881 = t1065 * t664;
t801 = t1098 * t881;
t888 = t661 * t1098;
t829 = t662 * t888;
t518 = t801 - t829;
t961 = t661 * t662;
t745 = t881 - t961;
t983 = t745 * qJ(5);
t437 = t518 + t983;
t427 = t1050 * t437;
t1138 = t1133 / 0.2e1 + t427 / 0.2e1;
t802 = t664 * t825;
t692 = -t802 + t829;
t680 = t692 - t983;
t425 = t1050 * t680;
t828 = t664 * t888;
t694 = -t828 - t803;
t436 = -t580 + t694;
t972 = t660 * t436;
t812 = t425 / 0.2e1 - t972 / 0.2e1;
t171 = t812 - t1138;
t1148 = qJD(3) * t171;
t836 = t1050 * t745 - t660 * t746;
t1143 = t836 / 0.2e1;
t665 = cos(qJ(2));
t1141 = t665 * t1143;
t562 = t745 * t665;
t844 = t1050 * t562;
t561 = t746 * t665;
t969 = t660 * t561;
t943 = t969 / 0.2e1 - t844 / 0.2e1;
t313 = t1141 + t943;
t1147 = t313 * qJD(2);
t1142 = -t1138 - t812;
t1146 = qJD(2) * t1142;
t312 = t1141 - t943;
t663 = sin(qJ(2));
t959 = t662 * t663;
t558 = t661 * t959 - t663 * t881;
t559 = t746 * t663;
t837 = -t1050 * t559 + t660 * t558;
t904 = t665 * qJD(1);
t879 = t837 * t904;
t1145 = t312 * qJD(2) + t879;
t1126 = t1050 * t746 + t660 * t745;
t465 = t1126 * qJD(5);
t1144 = -qJD(3) * t1142 - t465;
t1140 = t427 + t1133;
t1092 = t1140 / 0.2e1;
t1086 = -t836 / 0.2e1;
t1135 = t1126 / 0.2e1;
t472 = t1126 ^ 2;
t1127 = -t1050 * t558 - t660 * t559;
t1137 = t1127 ^ 2;
t1064 = pkin(7) * t662;
t816 = -pkin(2) * t665 - t663 * pkin(8);
t604 = -pkin(1) + t816;
t589 = t664 * t604;
t957 = t663 * t664;
t820 = -pkin(9) * t957 + t589;
t489 = (-pkin(3) - t1064) * t665 + t820;
t955 = t664 * t665;
t898 = pkin(7) * t955;
t543 = t662 * t604 + t898;
t508 = -pkin(9) * t959 + t543;
t884 = t1065 * t508;
t338 = t661 * t489 + t884;
t553 = t559 * qJ(5);
t730 = -t553 + t338;
t1116 = t1050 * t730;
t467 = t1065 * t489;
t963 = t661 * t508;
t337 = -t467 + t963;
t552 = t558 * qJ(5);
t305 = -t337 + t552;
t284 = -pkin(4) * t665 + t305;
t266 = t660 * t284;
t160 = t1116 + t266;
t954 = t665 * qJ(6);
t153 = t160 - t954;
t1096 = t153 / 0.2e1;
t1090 = -t837 / 0.2e1;
t1136 = -t1116 / 0.2e1;
t1089 = t1127 / 0.2e1;
t847 = t1050 * t1125;
t971 = t660 * t437;
t296 = t847 - t971;
t942 = pkin(3) * t959 + t663 * pkin(7);
t494 = pkin(4) * t559 + t942;
t236 = -pkin(5) * t837 - qJ(6) * t1127 + t494;
t1028 = t236 * t1127;
t902 = qJD(3) + qJD(4);
t1128 = t902 * t665;
t659 = t665 ^ 2;
t395 = t659 + t1137;
t929 = qJD(2) * t1126;
t1132 = -t395 * qJD(1) - t1127 * t929 + t1128;
t935 = qJD(1) * t1127;
t918 = qJD(6) * t1126;
t1131 = qJD(6) * t1127;
t402 = t1127 * qJD(5);
t1053 = t663 * pkin(5);
t1054 = t663 * pkin(3);
t1052 = t665 * pkin(8);
t1055 = t663 * pkin(2);
t609 = -t1052 + t1055;
t597 = t664 * t609;
t632 = pkin(7) * t959;
t496 = -pkin(9) * t955 + t1054 + t597 + t632;
t471 = t1065 * t496;
t596 = t662 * t609;
t897 = pkin(7) * t957;
t958 = t662 * t665;
t512 = -pkin(9) * t958 + t596 - t897;
t962 = t661 * t512;
t838 = t471 - t962;
t290 = pkin(4) * t663 - qJ(5) * t562 + t838;
t851 = t1050 * t290;
t470 = t661 * t496;
t504 = t1065 * t512;
t945 = t504 + t470;
t307 = -qJ(5) * t561 + t945;
t975 = t660 * t307;
t722 = t975 / 0.2e1 - t851 / 0.2e1;
t688 = -t1053 / 0.2e1 + t722;
t713 = t660 * t730;
t159 = t1050 * t284 - t713;
t1123 = t160 * t1135 + t1143 * t159;
t154 = t665 * pkin(5) - t159;
t1122 = t1086 * t154 + t153 * t1135;
t1130 = t929 + t935;
t1076 = -t746 / 0.2e1;
t1129 = t836 * t902;
t424 = t844 - t969;
t1004 = t424 * qJ(6);
t539 = t1050 * t561;
t968 = t660 * t562;
t420 = t539 + t968;
t1062 = t420 * pkin(5);
t1124 = -t1004 / 0.2e1 + t1062 / 0.2e1;
t294 = -t425 + t972;
t678 = t660 * t680;
t846 = t1050 * t436;
t297 = t846 + t678;
t690 = (-t296 + t297) * t836 + (-t1140 + t294) * t1126;
t1121 = qJD(2) * t690;
t786 = -t1126 * t296 + t1140 * t836;
t1119 = qJD(2) * t786;
t785 = t837 ^ 2 + t1137;
t1118 = qJD(5) * t785;
t1117 = qJD(5) * t786;
t681 = 0.2e1 * t1127 * t1135 + t836 * t837;
t1115 = t681 * qJD(5);
t1017 = t296 * t837;
t864 = -t1017 / 0.2e1;
t1114 = t864 - t1122;
t1113 = t864 - t1123;
t860 = -t266 / 0.2e1 + t1136;
t648 = -pkin(3) * t664 - pkin(2);
t1111 = t942 * t746 / 0.2e1 - t648 * t558 / 0.2e1;
t1110 = -t942 * t745 / 0.2e1 + t648 * t559 / 0.2e1;
t899 = pkin(7) * t958;
t507 = t820 - t899;
t964 = t661 * t507;
t357 = -t884 - t964;
t731 = t553 + t357;
t306 = t1050 * t731;
t885 = t1065 * t507;
t358 = t885 - t963;
t314 = t552 + t358;
t974 = t660 * t314;
t951 = t974 / 0.2e1 - t306 / 0.2e1;
t1109 = -t296 * t1089 + t1092 * t837;
t549 = -pkin(4) * t745 + t648;
t287 = -pkin(5) * t836 - qJ(6) * t1126 + t549;
t1108 = t1086 * t236 + t1090 * t287;
t1107 = t1086 * t160 + t1135 * t159;
t1106 = t1096 * t836 + t1135 * t154;
t656 = t662 ^ 2;
t658 = t664 ^ 2;
t623 = t658 - t656;
t905 = t663 * qJD(1);
t878 = t664 * t905;
t1105 = qJD(2) * t623 - 0.2e1 * t662 * t878;
t179 = -t306 + t974;
t714 = t660 * t731;
t849 = t1050 * t314;
t180 = t849 + t714;
t1104 = t1086 * t180 - t1089 * t294 + t1090 * t297 - t1135 * t179;
t976 = t660 * t305;
t164 = t1116 + t976;
t850 = t1050 * t305;
t165 = t850 - t713;
t1103 = t1086 * t165 - t1089 * t1140 + t1090 * t296 - t1135 * t164;
t784 = t836 ^ 2 + t472;
t1102 = qJD(1) * t681 + qJD(2) * t784;
t1101 = qJD(1) * t785 + qJD(2) * t681;
t1100 = qJD(3) * t690 + qJD(5) * t784;
t1099 = -pkin(5) / 0.2e1;
t1097 = qJ(6) / 0.2e1;
t1095 = -t160 / 0.2e1;
t1094 = t236 / 0.2e1;
t1093 = t287 / 0.2e1;
t1091 = t297 / 0.2e1;
t843 = t1050 * t661;
t896 = t1065 * pkin(3);
t647 = t896 + pkin(4);
t966 = t660 * t647;
t568 = pkin(3) * t843 + t966;
t556 = qJ(6) + t568;
t1084 = -t556 / 0.2e1;
t1083 = t556 / 0.2e1;
t965 = t660 * t661;
t627 = pkin(3) * t965;
t567 = t1050 * t647 - t627;
t557 = -pkin(5) - t567;
t1082 = t557 / 0.2e1;
t1081 = t558 / 0.2e1;
t1080 = -t567 / 0.2e1;
t1079 = t567 / 0.2e1;
t1078 = -t568 / 0.2e1;
t883 = t1065 * t660;
t578 = (t843 + t883) * pkin(3);
t1077 = t578 / 0.2e1;
t633 = pkin(3) * t957;
t1075 = t633 / 0.2e1;
t1057 = t660 * pkin(4);
t634 = qJ(6) + t1057;
t1074 = -t634 / 0.2e1;
t1073 = t634 / 0.2e1;
t635 = -pkin(4) * t1050 - pkin(5);
t1072 = -t635 / 0.2e1;
t1071 = t635 / 0.2e1;
t652 = t662 * pkin(3);
t1070 = t652 / 0.2e1;
t1069 = t660 / 0.2e1;
t1068 = -t663 / 0.2e1;
t1067 = -t665 / 0.2e1;
t1066 = t665 / 0.2e1;
t1063 = t1127 * pkin(5);
t1061 = t1126 * pkin(5);
t1060 = t558 * pkin(4);
t1059 = t561 * pkin(4);
t1058 = t746 * pkin(4);
t1056 = t661 * pkin(3);
t1051 = pkin(4) * qJD(4);
t1006 = t837 * qJ(6);
t244 = -t1006 - t1060 + t1063;
t55 = t164 * t665 - t244 * t837 + t1028;
t1049 = qJD(1) * t55;
t242 = t244 + t633;
t57 = t179 * t665 - t242 * t837 + t1028;
t1048 = qJD(1) * t57;
t63 = t1127 * t154 + t153 * t837;
t1047 = qJD(1) * t63;
t64 = -t1127 * t159 + t160 * t837;
t1046 = qJD(1) * t64;
t75 = t153 * t665 + t1028;
t1045 = qJD(1) * t75;
t1036 = t164 * t296;
t1035 = t165 * t1140;
t1033 = t180 * t1140;
t276 = t660 * t290;
t292 = t1050 * t307;
t162 = t292 + t276;
t157 = t663 * qJ(6) + t162;
t161 = t851 - t975;
t158 = -t161 - t1053;
t630 = pkin(3) * t958;
t654 = t665 * pkin(7);
t598 = t654 + t630;
t495 = t598 + t1059;
t237 = -t1004 + t495 + t1062;
t21 = t153 * t157 + t154 * t158 + t236 * t237;
t1031 = t21 * qJD(1);
t22 = t153 * t165 + t154 * t164 + t236 * t244;
t1030 = t22 * qJD(1);
t23 = t153 * t180 + t154 * t179 + t236 * t242;
t1029 = t23 * qJD(1);
t24 = t1127 * t158 - t153 * t420 + t154 * t424 + t157 * t837;
t1026 = t24 * qJD(1);
t791 = t1127 * t164 + t165 * t837;
t793 = -t1127 * t153 + t154 * t837;
t25 = t791 + t793;
t1025 = t25 * qJD(1);
t26 = -t1127 * t161 - t159 * t424 - t160 * t420 + t162 * t837;
t1024 = t26 * qJD(1);
t792 = -t1127 * t160 - t159 * t837;
t27 = t791 + t792;
t1023 = t27 * qJD(1);
t790 = t1127 * t179 + t180 * t837;
t28 = t790 + t793;
t1022 = t28 * qJD(1);
t1020 = t287 * t1126;
t1019 = t287 * t836;
t29 = t790 + t792;
t1018 = t29 * qJD(1);
t1016 = t296 * t578;
t1013 = t1140 * t1127;
t579 = t1050 * t896 - t627;
t1012 = t1140 * t579;
t40 = t159 * t161 + t160 * t162 + t494 * t495;
t1011 = t40 * qJD(1);
t41 = -t1060 * t494 - t159 * t164 + t160 * t165;
t1010 = t41 * qJD(1);
t1008 = t1127 * t660;
t509 = t633 - t1060;
t42 = -t159 * t179 + t160 * t180 + t494 * t509;
t1007 = t42 * qJD(1);
t45 = t1127 * t237 - t153 * t663 + t157 * t665 + t236 * t424;
t1003 = t45 * qJD(1);
t46 = -t154 * t663 + t158 * t665 + t236 * t420 - t237 * t837;
t1002 = t46 * qJD(1);
t1001 = t1126 * t660;
t1000 = t1126 * t665;
t999 = t836 * qJ(6);
t997 = t693 * t665;
t996 = t518 * t665;
t995 = t556 * t1127;
t994 = t556 * t1126;
t993 = t557 * t837;
t992 = t557 * t836;
t183 = t236 * t837;
t56 = t1127 * t244 + t165 * t665 + t183;
t991 = t56 * qJD(1);
t990 = t567 * t837;
t989 = t567 * t836;
t988 = t568 * t1127;
t987 = t568 * t1126;
t986 = t579 * t837;
t985 = t579 * t836;
t58 = t1127 * t242 + t180 * t665 + t183;
t984 = t58 * qJD(1);
t982 = t634 * t1127;
t981 = t634 * t1126;
t980 = t635 * t837;
t979 = t635 * t836;
t657 = t663 ^ 2;
t956 = t664 * t657;
t263 = t558 * t746 - t559 * t745;
t952 = t902 * t263;
t344 = t1076 * t559 - t1081 * t745;
t949 = t902 * t344;
t624 = t659 - t657;
t819 = pkin(3) * t1066 - t489 / 0.2e1;
t208 = (t507 / 0.2e1 + t819) * t661;
t941 = qJD(1) * t208;
t502 = t942 * t558;
t238 = t357 * t665 - t559 * t633 + t502;
t940 = qJD(1) * t238;
t842 = t942 * t559;
t239 = t358 * t665 - t558 * t633 - t842;
t939 = qJD(1) * t239;
t256 = -t337 * t665 - t842;
t938 = qJD(1) * t256;
t257 = -t338 * t665 + t502;
t937 = qJD(1) * t257;
t936 = qJD(1) * t312;
t428 = t559 * t663 - t561 * t665;
t934 = qJD(1) * t428;
t429 = t558 * t663 + t562 * t665;
t933 = qJD(1) * t429;
t542 = -t589 + t899;
t490 = -t1064 * t657 - t542 * t665;
t932 = qJD(1) * t490;
t491 = -pkin(7) * t956 - t543 * t665;
t931 = qJD(1) * t491;
t930 = qJD(1) * t558;
t928 = qJD(2) * t746;
t927 = qJD(2) * t648;
t926 = qJD(2) * t662;
t925 = qJD(2) * t663;
t924 = qJD(2) * t664;
t923 = qJD(2) * t665;
t922 = qJD(3) * t662;
t921 = qJD(3) * t664;
t920 = qJD(3) * t665;
t919 = qJD(4) * t648;
t917 = qJD(6) * t665;
t142 = t337 * t663 - t598 * t559 - t561 * t942 + t665 * t838;
t916 = t142 * qJD(1);
t143 = -t338 * t663 - t598 * t558 + t562 * t942 + t665 * t945;
t915 = t143 * qJD(1);
t826 = t896 / 0.2e1;
t737 = -t467 / 0.2e1 + t665 * t826;
t210 = t885 / 0.2e1 + t737;
t914 = t210 * qJD(1);
t810 = t539 / 0.2e1 + t968 / 0.2e1;
t311 = -t1000 / 0.2e1 + t810;
t913 = t311 * qJD(1);
t912 = t313 * qJD(1);
t331 = t558 * t561 - t559 * t562;
t911 = t331 * qJD(1);
t364 = t542 * t663 + (-t632 + t597) * t665;
t910 = t364 * qJD(1);
t365 = t596 * t665 + (-t543 + t898) * t663;
t909 = t365 * qJD(1);
t403 = t837 * qJD(5);
t729 = -t960 / 0.2e1 - t882 / 0.2e1;
t457 = (t1076 + t729) * t665;
t448 = t457 * qJD(1);
t728 = t881 / 0.2e1 - t961 / 0.2e1;
t458 = (t745 / 0.2e1 + t728) * t665;
t449 = t458 * qJD(1);
t594 = t624 * t662;
t908 = t594 * qJD(1);
t595 = t659 * t664 - t956;
t907 = t595 * qJD(1);
t906 = t624 * qJD(1);
t903 = t579 * qJD(4) + qJD(6);
t547 = -t1060 / 0.2e1;
t901 = t547 + t1075;
t577 = t1058 / 0.2e1;
t900 = t577 + t1070;
t895 = pkin(1) * t905;
t894 = pkin(1) * t904;
t892 = -t652 / 0.2e1;
t891 = -t1054 / 0.2e1;
t890 = t1082 + t1099;
t889 = t1071 + t1099;
t887 = -t1050 / 0.2e1;
t886 = -t954 / 0.2e1 - t860;
t880 = t837 * t935;
t877 = t662 * t924;
t876 = t663 * t924;
t875 = t662 * t920;
t874 = t664 * t920;
t873 = t662 * t921;
t872 = t663 * t923;
t871 = t663 * t904;
t866 = -t1035 / 0.2e1;
t865 = -t179 * t296 / 0.2e1;
t862 = t1127 * t1077;
t861 = t1126 * t1077;
t859 = t276 / 0.2e1 + t292 / 0.2e1;
t858 = -t470 / 0.2e1 - t504 / 0.2e1;
t857 = t1084 + t1073;
t856 = t1083 + t1097;
t855 = t1082 + t1072;
t854 = t1073 + t1097;
t853 = t1065 * qJD(3);
t852 = t1065 * qJD(4);
t848 = t1050 * t837;
t845 = t1050 * t836;
t839 = -t1140 * t420 - t296 * t424;
t835 = t403 - t917;
t833 = t902 * t746;
t831 = t661 * t891;
t830 = -qJD(3) + t904;
t824 = t662 * t876;
t363 = t986 / 0.2e1;
t822 = t363 + t862;
t387 = t985 / 0.2e1;
t821 = t387 + t861;
t818 = t904 - qJD(3) / 0.2e1;
t817 = t849 / 0.2e1;
t815 = t856 * t663;
t814 = t854 * t663;
t813 = t1136 - t976 / 0.2e1;
t811 = t471 / 0.2e1 - t962 / 0.2e1;
t809 = t663 * t826;
t308 = t1058 - t999 + t1061;
t291 = t308 + t652;
t675 = t153 * t1091 + t154 * t294 / 0.2e1 + t865 + t291 * t1094 + t242 * t1093;
t761 = t1082 * t158 + t1083 * t157;
t1 = -t1033 / 0.2e1 - t675 + t761;
t788 = t1140 * t297 - t294 * t296;
t53 = t287 * t291 + t788;
t808 = -t1 * qJD(1) + t53 * qJD(2);
t676 = t296 * t1096 + t154 * t1092 - t1036 / 0.2e1 + t308 * t1094 + t244 * t1093;
t760 = t1071 * t158 + t1073 * t157;
t3 = t866 - t676 + t760;
t54 = t287 * t308;
t807 = -t3 * qJD(1) + t54 * qJD(2);
t101 = t1033 / 0.2e1;
t554 = t652 + t1058;
t672 = t101 - t159 * t294 / 0.2e1 + t160 * t1091 + t865 + t494 * t554 / 0.2e1 + t509 * t549 / 0.2e1;
t759 = t1078 * t162 + t1080 * t161;
t6 = t672 + t759;
t65 = t549 * t554 + t788;
t806 = t6 * qJD(1) + t65 * qJD(2);
t753 = t1017 / 0.2e1 + t1013 / 0.2e1;
t670 = t753 + t1104;
t749 = t1082 * t424 + t1084 * t420;
t9 = t670 + t749 + t1122;
t805 = -t9 * qJD(1) + t1121;
t66 = t1058 * t549;
t703 = t159 * t1092 + t296 * t1095 + t1036 / 0.2e1;
t724 = t162 * t1069 + t161 * t1050 / 0.2e1;
t7 = t866 + (t1076 * t494 + t1081 * t549 + t724) * pkin(4) + t703;
t804 = -t7 * qJD(1) + t66 * qJD(2);
t612 = -qJD(4) + t830;
t800 = t859 + t1108;
t799 = t858 - t1110;
t671 = t753 + t1103;
t747 = t1071 * t424 + t1074 * t420;
t11 = t671 + t747 + t1122;
t798 = t11 * qJD(1);
t748 = t1078 * t420 + t1080 * t424;
t13 = t670 + t748 + t1123;
t797 = -t13 * qJD(1) + t1121;
t706 = (-t660 * t420 / 0.2e1 + t424 * t887) * pkin(4);
t15 = t706 + t671 + t1123;
t796 = t15 * qJD(1);
t77 = t1127 * t857 + t837 * t855 + t822;
t85 = t1126 * t857 + t836 * t855 + t821;
t795 = qJD(1) * t77 + qJD(2) * t85;
t698 = -t990 / 0.2e1 - t988 / 0.2e1 + t862;
t710 = (-t1008 / 0.2e1 - t848 / 0.2e1) * pkin(4);
t78 = -t986 / 0.2e1 + t710 - t698;
t697 = -t989 / 0.2e1 - t987 / 0.2e1 + t861;
t708 = (-t1001 / 0.2e1 - t845 / 0.2e1) * pkin(4);
t86 = -t985 / 0.2e1 + t708 - t697;
t794 = qJD(1) * t78 + qJD(2) * t86;
t104 = -t291 * t836 + t1020;
t755 = t1093 * t1127 + t1135 * t236;
t674 = t1066 * t294 + t1086 * t242 + t1090 * t291 + t755;
t32 = t663 * t890 + t674 + t722;
t783 = -qJD(1) * t32 - qJD(2) * t104;
t105 = -t1126 * t291 - t1019;
t696 = t859 - t1108;
t701 = t1066 * t297 + t1089 * t291 + t1135 * t242;
t30 = t815 + t696 + t701;
t782 = qJD(1) * t30 - qJD(2) * t105;
t764 = t630 / 0.2e1 + t654 / 0.2e1 + t1059 / 0.2e1;
t685 = t764 - t1109;
t38 = t685 - t1106 + t1124;
t781 = qJD(1) * t38 - t1119;
t43 = t685 + t1107;
t780 = qJD(1) * t43 - t1119;
t108 = -t308 * t836 + t1020;
t673 = t1066 * t1140 + t1086 * t244 + t1090 * t308 + t755;
t36 = t663 * t889 + t673 + t722;
t779 = -qJD(1) * t36 - qJD(2) * t108;
t109 = -t1126 * t308 - t1019;
t700 = t1066 * t296 + t1089 * t308 + t1135 * t244;
t34 = t814 + t696 + t700;
t778 = qJD(1) * t34 - qJD(2) * t109;
t773 = t830 * t663;
t102 = -t1127 * t890 - t837 * t856 + t901;
t144 = -t1126 * t890 - t836 * t856 + t900;
t772 = qJD(1) * t102 + qJD(2) * t144;
t121 = -t1127 * t889 - t837 * t854 + t547;
t181 = -t1126 * t889 - t836 * t854 + t577;
t771 = qJD(1) * t121 + qJD(2) * t181;
t667 = -t664 * t746 * t891 + t1066 * t694 + t558 * t892;
t686 = t858 + t1110;
t146 = t831 - t667 + t686;
t452 = t648 * t745 + t652 * t746;
t770 = t146 * qJD(1) - t452 * qJD(2);
t666 = t1067 * t692 - t1075 * t745 - t559 * t892;
t684 = t811 - t1111;
t147 = t809 - t666 + t684;
t451 = t648 * t746 - t652 * t745;
t769 = t147 * qJD(1) - t451 * qJD(2);
t752 = t1078 * t837 + t1079 * t1127;
t196 = t752 + t901;
t750 = t1078 * t836 + t1079 * t1126;
t234 = t750 + t900;
t768 = qJD(1) * t196 + qJD(2) * t234;
t709 = (t1069 * t837 + t1127 * t887) * pkin(4);
t240 = t1060 / 0.2e1 + t709;
t707 = (t1069 * t836 + t1126 * t887) * pkin(4);
t277 = -t1058 / 0.2e1 + t707;
t767 = qJD(1) * t240 + qJD(2) * t277;
t366 = -t558 ^ 2 + t559 ^ 2;
t193 = qJD(1) * t366 + qJD(2) * t263;
t438 = t745 ^ 2 - t746 ^ 2;
t206 = qJD(1) * t263 + qJD(2) * t438;
t320 = qJD(1) * t837 + qJD(2) * t836;
t765 = t1052 / 0.2e1 - t1055 / 0.2e1;
t733 = t765 * t662;
t505 = t596 / 0.2e1 - t733;
t763 = pkin(2) * t924 - qJD(1) * t505;
t732 = t765 * t664;
t506 = -t597 / 0.2e1 + t732;
t762 = pkin(2) * t926 - qJD(1) * t506;
t756 = t1072 * t179 + t1074 * t180;
t754 = t1072 * t294 + t1074 * t297;
t751 = t1086 * t1127 - t1135 * t837;
t702 = -t236 * t1126 / 0.2e1 - t287 * t1127 / 0.2e1 + t1140 * t1067;
t51 = t688 - t702;
t744 = qJD(1) * t51 + t287 * t929;
t743 = t664 * t773;
t188 = t1068 + t751;
t742 = qJD(1) * t188 - t836 * t929;
t200 = -t996 / 0.2e1 + t684;
t741 = qJD(1) * t200 - t746 * t927;
t201 = -t997 / 0.2e1 + t686;
t740 = qJD(1) * t201 - t745 * t927;
t253 = -qJD(2) * t344 - t559 * t930;
t264 = qJD(1) * t344 + t745 * t928;
t581 = (t656 / 0.2e1 - t658 / 0.2e1) * t663;
t738 = -qJD(1) * t581 + t877;
t736 = t811 + t1111;
t735 = t547 + t1063 / 0.2e1 - t1006 / 0.2e1;
t734 = t577 + t1061 / 0.2e1 - t999 / 0.2e1;
t727 = qJD(1) * t662 * t956 + qJD(2) * t581;
t593 = t623 * t657;
t726 = qJD(1) * t593 + 0.2e1 * t824;
t723 = t579 * t1066 - t850 / 0.2e1;
t721 = -t971 / 0.2e1 + t847 / 0.2e1;
t720 = t1067 * t578 - t813;
t683 = t1077 * t154 + t1082 * t164 + t1083 * t165 + t1096 * t579;
t18 = t683 + t756;
t333 = t556 * t579 + t557 * t578;
t251 = t1012 / 0.2e1;
t689 = t251 + t1140 * t1082 - t1016 / 0.2e1 + t296 * t1083;
t48 = t689 + t754;
t719 = t18 * qJD(1) + t48 * qJD(2) + t333 * qJD(3);
t682 = t1077 * t159 + t1078 * t165 + t1079 * t164 + t1095 * t579;
t712 = (t1069 * t180 + t179 * t887) * pkin(4);
t19 = t712 + t682;
t339 = -t567 * t578 + t568 * t579;
t699 = t1140 * t1079 + t1016 / 0.2e1 + t296 * t1078;
t711 = (t1069 * t297 + t294 * t887) * pkin(4);
t49 = -t1012 / 0.2e1 + t711 + t699;
t718 = t19 * qJD(1) + t49 * qJD(2) - t339 * qJD(3);
t71 = t720 - t951;
t717 = qJD(1) * t71 + qJD(3) * t578 - t1146;
t67 = t665 * t856 + t860 + t951;
t716 = -qJD(1) * t67 + qJD(3) * t556 - t1146;
t669 = -t678 / 0.2e1 - t846 / 0.2e1;
t174 = t669 + t721;
t72 = t817 + (-t507 / 0.2e1 + t489 / 0.2e1) * t965 + t723;
t715 = qJD(1) * t72 - qJD(2) * t174 - qJD(3) * t579;
t229 = -t1013 / 0.2e1;
t705 = t229 - t1103;
t704 = t229 - t1104;
t695 = t764 + t1109;
t519 = -qJ(6) + (t826 - pkin(4) / 0.2e1 - t647 / 0.2e1) * t660;
t69 = t665 * t854 - t813 + t860;
t687 = qJD(1) * t69 + qJD(3) * t519 - qJD(4) * t634;
t639 = t925 / 0.2e1;
t638 = -t905 / 0.2e1;
t637 = t905 / 0.2e1;
t588 = t818 * t663;
t575 = t581 * qJD(3);
t571 = t578 * qJD(4);
t555 = (-qJD(4) / 0.2e1 + t818) * t663;
t482 = t1057 / 0.2e1 + qJ(6) + t966 / 0.2e1 + (t843 + t883 / 0.2e1) * pkin(3);
t466 = t836 * qJD(5);
t462 = t632 + t597 / 0.2e1 + t732;
t461 = t897 - t596 / 0.2e1 - t733;
t460 = -t1067 * t746 + t665 * t729;
t459 = -t1066 * t745 + t665 * t728;
t440 = -t802 / 0.2e1 - t801 / 0.2e1 + t829;
t439 = t828 / 0.2e1 + t827 / 0.2e1 + t803;
t371 = t836 * qJD(6);
t360 = t458 * qJD(2) - t559 * t904;
t359 = t457 * qJD(2) + t558 * t904;
t336 = -t833 - t448;
t335 = t745 * t902 - t449;
t334 = t837 * qJD(6);
t330 = t1127 * t918;
t316 = t1000 / 0.2e1 + t810;
t278 = t577 + t707;
t255 = t459 * qJD(2) + t559 * t612;
t254 = t460 * qJD(2) - t558 * t612;
t243 = qJD(2) * t472 + t1126 * t935;
t241 = t547 + t709;
t235 = -t750 + t900;
t211 = t963 - t885 / 0.2e1 + t737;
t209 = -t884 - t964 / 0.2e1 + t819 * t661;
t203 = t996 / 0.2e1 + t736;
t202 = t997 / 0.2e1 + t799;
t197 = -t752 + t901;
t189 = t1068 - t751;
t182 = t1071 * t1126 + t1073 * t836 + t734;
t175 = 0.2e1 * t1138;
t173 = -t669 + t721;
t149 = t831 + t667 + t799;
t148 = t809 + t666 + t736;
t145 = t1082 * t1126 + t1083 * t836 + t1070 + t734;
t122 = t1071 * t1127 + t1073 * t837 + t735;
t103 = t1082 * t1127 + t1083 * t837 + t1075 + t735;
t98 = t1035 / 0.2e1;
t87 = t387 + t708 + t697;
t84 = -t994 / 0.2e1 + t992 / 0.2e1 - t981 / 0.2e1 + t979 / 0.2e1 + t821;
t79 = t363 + t710 + t698;
t76 = -t995 / 0.2e1 + t993 / 0.2e1 - t982 / 0.2e1 + t980 / 0.2e1 + t822;
t74 = -t713 / 0.2e1 + t817 + t714 / 0.2e1 - t723;
t73 = -t720 - t951;
t70 = t1067 * t634 - t813 + t886;
t68 = t1067 * t556 + t886 + t951;
t52 = t688 + t702;
t50 = t251 + t711 - t699;
t47 = t689 - t754;
t44 = t695 - t1107;
t39 = t695 + t1106 + t1124;
t37 = t1068 * t635 + t673 - t688;
t35 = t814 - t700 + t800;
t33 = t1068 * t557 + t674 - t688;
t31 = t815 - t701 + t800;
t20 = t712 - t682;
t17 = t683 - t756;
t16 = t706 + t705 + t1113;
t14 = t704 + t748 + t1113;
t12 = t705 + t747 + t1114;
t10 = t704 + t749 + t1114;
t8 = pkin(4) * t724 + t494 * t577 + t547 * t549 - t703 + t98;
t5 = t672 - t759;
t4 = t98 + t676 + t760;
t2 = t101 + t675 + t761;
t59 = [0, 0, 0, t872, t624 * qJD(2), 0, 0, 0, -pkin(1) * t925, -pkin(1) * t923, -t657 * t873 + t658 * t872, -t593 * qJD(3) - 0.2e1 * t665 * t824, -t595 * qJD(2) + t663 * t875, t594 * qJD(2) + t663 * t874, -t872, -qJD(2) * t364 - qJD(3) * t491, qJD(2) * t365 + qJD(3) * t490 (-qJD(2) * t562 + t559 * t902) * t558, qJD(2) * t331 + t366 * t902, -t429 * qJD(2) + t1128 * t559, -t428 * qJD(2) - t1128 * t558, -t872, -qJD(2) * t142 - qJD(3) * t238 - qJD(4) * t257, qJD(2) * t143 + qJD(3) * t239 + qJD(4) * t256, qJD(2) * t26 + qJD(3) * t29 + qJD(4) * t27 + t1118, qJD(2) * t40 + qJD(3) * t42 + qJD(4) * t41 + qJD(5) * t64, t46 * qJD(2) + t57 * qJD(3) + t55 * qJD(4) + t1131 * t837 + t402 * t665, t24 * qJD(2) + t28 * qJD(3) + t25 * qJD(4) - t837 * t917 + t1118, -t45 * qJD(2) - t58 * qJD(3) - t56 * qJD(4) + t395 * qJD(6) - t403 * t665, qJD(2) * t21 + qJD(3) * t23 + qJD(4) * t22 + qJD(5) * t63 - qJD(6) * t75; 0, 0, 0, t871, t906, t923, -t925, 0, -pkin(7) * t923 - t895, pkin(7) * t925 - t894, -t575 + (t658 * t905 + t877) * t665, t1105 * t665 - 0.2e1 * t663 * t873, t662 * t925 - t907, t876 + t908, -t588, -t910 + (t662 * t816 - t898) * qJD(2) + t462 * qJD(3), t909 + (t664 * t816 + t899) * qJD(2) + t461 * qJD(3) (t928 - t930) * t562 + t949, t911 + (-t561 * t746 + t562 * t745) * qJD(2) + t952, t459 * t902 + t746 * t925 - t933, t460 * t902 + t745 * t925 - t934, -t555, -t916 + (t561 * t648 - t598 * t745 + t663 * t693) * qJD(2) + t148 * qJD(3) + t203 * qJD(4), t915 + (-t518 * t663 + t562 * t648 + t598 * t746) * qJD(2) + t149 * qJD(3) + t202 * qJD(4), t1024 + (-t1126 * t161 + t162 * t836 + t839) * qJD(2) + t14 * qJD(3) + t16 * qJD(4) + t1115, t1011 + (t1140 * t162 + t161 * t296 + t495 * t549) * qJD(2) + t5 * qJD(3) + t8 * qJD(4) + t44 * qJD(5), t1002 + (-t237 * t836 + t287 * t420 + t296 * t663) * qJD(2) + t33 * qJD(3) + t37 * qJD(4) + t316 * qJD(5) + t189 * qJD(6), t1026 + (t1126 * t158 + t157 * t836 + t839) * qJD(2) + t10 * qJD(3) + t12 * qJD(4) + t1115 - t313 * qJD(6), -t1003 + (-t1126 * t237 + t1140 * t663 - t287 * t424) * qJD(2) + t31 * qJD(3) + t35 * qJD(4) - t312 * qJD(5) + t330, t1031 + (t1140 * t157 - t158 * t296 + t237 * t287) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t39 * qJD(5) + t52 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t727, -t726, t662 * t773, t743, t639, qJD(2) * t462 - qJD(3) * t543 - t931, qJD(2) * t461 + qJD(3) * t542 + t932, -t253, t193, t255, t254, t639, qJD(2) * t148 + qJD(3) * t357 + qJD(4) * t209 - t940, qJD(2) * t149 - qJD(3) * t358 + qJD(4) * t211 + t939, t1018 + t14 * qJD(2) + (-t988 - t990) * qJD(3) + t79 * qJD(4), t1007 + t5 * qJD(2) + (-t179 * t567 + t180 * t568) * qJD(3) + t20 * qJD(4) + t197 * qJD(5), qJD(2) * t33 - qJD(3) * t179 + qJD(4) * t73 + t1048, t1022 + t10 * qJD(2) + (t993 - t995) * qJD(3) + t76 * qJD(4) + t334, t31 * qJD(2) + t180 * qJD(3) + t74 * qJD(4) - t917 - t984, t1029 + t2 * qJD(2) + (t179 * t557 + t180 * t556) * qJD(3) + t17 * qJD(4) + t103 * qJD(5) + t68 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t193, t255, t254, t639, qJD(2) * t203 + qJD(3) * t209 - qJD(4) * t338 - t937, qJD(2) * t202 + qJD(3) * t211 + qJD(4) * t337 + t938, t1023 + t16 * qJD(2) + t79 * qJD(3) + (-t848 - t1008) * t1051, t1010 + t8 * qJD(2) + t20 * qJD(3) + (-t1050 * t164 + t165 * t660) * t1051 + t241 * qJD(5), qJD(2) * t37 + qJD(3) * t73 - qJD(4) * t164 + t1049, t1025 + t12 * qJD(2) + t76 * qJD(3) + (t980 - t982) * qJD(4) + t334, t35 * qJD(2) + t74 * qJD(3) + t165 * qJD(4) - t917 - t991, t1030 + t4 * qJD(2) + t17 * qJD(3) + (t164 * t635 + t165 * t634) * qJD(4) + t122 * qJD(5) + t70 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1101, qJD(2) * t44 + qJD(3) * t197 + qJD(4) * t241 + t1046, t316 * qJD(2) + t1127 * t904, t1101, -t1145, qJD(2) * t39 + qJD(3) * t103 + qJD(4) * t122 + t1047; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t189 + t880, t837 * t902 - t1147 - t879, -t1132, qJD(2) * t52 + qJD(3) * t68 + qJD(4) * t70 - t1045; 0, 0, 0, -t871, -t906, 0, 0, 0, t895, t894, -t658 * t871 - t575, 0.2e1 * t662 * t743, -t874 + t907, t875 - t908, t588, qJD(3) * t506 + t910, qJD(3) * t505 - t909, t562 * t930 + t949, -t911 + t952, -t458 * t902 + t933, -t457 * t902 + t934, t555, -qJD(3) * t147 - qJD(4) * t200 + t916, -qJD(3) * t146 - qJD(4) * t201 - t915, -qJD(3) * t13 - qJD(4) * t15 - t1024 + t1115, qJD(3) * t6 - qJD(4) * t7 - qJD(5) * t43 - t1011, qJD(3) * t32 + qJD(4) * t36 - qJD(5) * t311 - qJD(6) * t188 - t1002, -qJD(3) * t9 - qJD(4) * t11 - qJD(6) * t312 - t1026 + t1115, -qJD(3) * t30 - qJD(4) * t34 - qJD(5) * t313 + t1003 + t330, -qJD(3) * t1 - qJD(4) * t3 - qJD(5) * t38 - qJD(6) * t51 - t1031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t873, t623 * qJD(3), 0, 0, 0, -pkin(2) * t922, -pkin(2) * t921, t745 * t833, t902 * t438, 0, 0, 0, qJD(3) * t451 + t746 * t919, qJD(3) * t452 + t745 * t919, t1100, qJD(3) * t65 + qJD(4) * t66 + t1117, qJD(3) * t104 + qJD(4) * t108 + t836 * t918, t1100, qJD(3) * t105 + qJD(4) * t109 + qJD(6) * t472, qJD(3) * t53 + qJD(4) * t54 - t287 * t918 + t1117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t738, t1105, -t830 * t664, t830 * t662, t638, -pkin(8) * t921 - t762, pkin(8) * t922 - t763, t264, t206, t335, t336, t638, qJD(3) * t692 + t440 * qJD(4) - t769, -qJD(3) * t694 + t439 * qJD(4) - t770 (-t987 - t989) * qJD(3) + t87 * qJD(4) + t797 (-t294 * t567 + t297 * t568) * qJD(3) + t50 * qJD(4) + t235 * qJD(5) + t806, -qJD(3) * t294 + qJD(4) * t171 - t783 (t992 - t994) * qJD(3) + t84 * qJD(4) + t371 + t805, qJD(3) * t297 + qJD(4) * t173 - t782 (t294 * t557 + t297 * t556) * qJD(3) + t47 * qJD(4) + t145 * qJD(5) - t171 * qJD(6) + t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t206, t335, t336, t638, qJD(3) * t440 - qJD(4) * t518 - t741, qJD(3) * t439 - qJD(4) * t693 - t740, t87 * qJD(3) + (-t845 - t1001) * t1051 - t796, t50 * qJD(3) + (-t1050 * t1140 + t296 * t660) * t1051 + t278 * qJD(5) + t804, -qJD(4) * t1140 + t1148 - t779, t84 * qJD(3) + (t979 - t981) * qJD(4) + t371 - t798, qJD(3) * t173 + qJD(4) * t296 - t778, t47 * qJD(3) + (t1140 * t635 + t296 * t634) * qJD(4) + t182 * qJD(5) + t175 * qJD(6) + t807; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1102, qJD(3) * t235 + qJD(4) * t278 - t780, -t913, t1102, -t912, qJD(3) * t145 + qJD(4) * t182 - t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t742, -t936 + t1129, t243, qJD(4) * t175 - t1148 - t744; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t727, t726 (-t662 * t905 + t924) * t665 (-t878 - t926) * t665, t639, -qJD(2) * t506 + t931, -qJD(2) * t505 - t932, t253, -t193, t360, t359, t639, qJD(2) * t147 + qJD(4) * t208 + t940, qJD(2) * t146 + qJD(4) * t210 - t939, qJD(2) * t13 - qJD(4) * t78 - t1018, -qJD(2) * t6 - qJD(4) * t19 - qJD(5) * t196 - t1007, -qJD(2) * t32 - qJD(4) * t71 - t1048 - t402, qJD(2) * t9 + qJD(4) * t77 - t1022, t30 * qJD(2) - t72 * qJD(4) + t835 + t984, qJD(2) * t1 + qJD(4) * t18 - qJD(5) * t102 - qJD(6) * t67 - t1029; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t738, -t1105, t664 * t904, -t662 * t904, t637, t762, t763, -t264, -t206, t449, t448, t637, t769, t770, -qJD(4) * t86 - t797, -qJD(4) * t49 - qJD(5) * t234 - t806, qJD(4) * t1142 - t465 + t783, qJD(4) * t85 - t805, qJD(4) * t174 + t466 + t782, qJD(4) * t48 - qJD(5) * t144 - qJD(6) * t1142 - t808; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t1056, -pkin(3) * t852, 0, t339 * qJD(4), -t571, 0, t903, qJD(4) * t333 + qJD(6) * t556; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1056 * t902 + t941, t914 + (-t853 - t852) * pkin(3), -t794 (-t1050 * t578 + t579 * t660) * t1051 - t718, -t571 - t717, t795, -t715 + t903 (t578 * t635 + t579 * t634) * qJD(4) + t482 * qJD(6) + t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t768, -t1130, 0, t320, -t772; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t612, qJD(4) * t482 + t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t193, t360, t359, t639, qJD(2) * t200 - qJD(3) * t208 + t937, qJD(2) * t201 - qJD(3) * t210 - t938, qJD(2) * t15 + qJD(3) * t78 - t1023, qJD(2) * t7 + qJD(3) * t19 + qJD(5) * t240 - t1010, -qJD(2) * t36 + qJD(3) * t71 - t1049 - t402, qJD(2) * t11 - qJD(3) * t77 - t1025, t34 * qJD(2) + t72 * qJD(3) + t835 + t991, qJD(2) * t3 - qJD(3) * t18 - qJD(5) * t121 - qJD(6) * t69 - t1030; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t264, -t206, t449, t448, t637, t741, t740, qJD(3) * t86 + t796, qJD(3) * t49 + qJD(5) * t277 - t804, t1144 + t779, -qJD(3) * t85 + t798, -qJD(3) * t174 + t466 + t778, -qJD(3) * t48 - qJD(5) * t181 - t807; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t1056 - t941, pkin(3) * t853 - t914, t794, t718, t717, -t795, qJD(6) + t715, -qJD(6) * t519 - t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), t634 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t767, -t1130, 0, t320, -t771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t612, -t687; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1101, qJD(2) * t43 + qJD(3) * t196 - qJD(4) * t240 - t1046, t311 * qJD(2) - t1127 * t612, -t1101, t612 * t837 + t1147, qJD(2) * t38 + qJD(3) * t102 + qJD(4) * t121 - t1047 - t1131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1102, qJD(3) * t234 - qJD(4) * t277 + t780, t1126 * t902 + t913, -t1102, t912 - t1129, qJD(3) * t144 + qJD(4) * t181 + t781 - t918; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t768, t1130, 0, -t320, t772; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t767, t1130, 0, -t320, t771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t188 - t880, t1145, t1132, qJD(2) * t51 + qJD(3) * t67 + qJD(4) * t69 + t1045 + t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t742, t936, -t243, -t1144 + t744; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t612, qJD(4) * t519 - t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t612, t687; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t59;