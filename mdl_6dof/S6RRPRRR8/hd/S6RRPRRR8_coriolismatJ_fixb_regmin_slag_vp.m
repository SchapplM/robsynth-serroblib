% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x35]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRRR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:35
% EndTime: 2019-03-09 14:06:36
% DurationCPUTime: 41.03s
% Computational Cost: add. (26577->870), mult. (59146->1106), div. (0->0), fcn. (69061->10), ass. (0->654)
t601 = sin(pkin(11));
t976 = cos(qJ(4));
t799 = t976 * t601;
t602 = cos(pkin(11));
t605 = sin(qJ(4));
t897 = t605 * t602;
t552 = t799 + t897;
t607 = cos(qJ(5));
t539 = t607 * t552;
t604 = sin(qJ(5));
t798 = t976 * t602;
t898 = t605 * t601;
t674 = t798 - t898;
t900 = t604 * t674;
t1030 = t539 + t900;
t966 = pkin(8) + qJ(3);
t564 = t966 * t601;
t565 = t966 * t602;
t486 = t976 * t564 + t565 * t605;
t681 = -pkin(9) * t552 - t486;
t1041 = t607 * t681;
t487 = -t605 * t564 + t565 * t976;
t428 = pkin(9) * t674 + t487;
t904 = t604 * t428;
t719 = t1041 - t904;
t1122 = -pkin(10) * t1030 + t719;
t603 = sin(qJ(6));
t1141 = t603 * t1122;
t1031 = -t604 * t552 + t607 * t674;
t300 = t607 * t428 + t604 * t681;
t248 = pkin(10) * t1031 + t300;
t975 = cos(qJ(6));
t1153 = t975 * t248;
t240 = -t1153 / 0.2e1;
t1169 = -t1141 / 0.2e1 + t240;
t1170 = 0.2e1 * t1169;
t1118 = qJD(4) + qJD(5);
t1060 = qJD(5) + qJD(6);
t1159 = qJD(4) + t1060;
t606 = sin(qJ(2));
t526 = t552 * t606;
t529 = t674 * t606;
t717 = -t526 * t607 - t604 * t529;
t1043 = t603 * t717;
t511 = t607 * t529;
t903 = t604 * t526;
t1032 = t511 - t903;
t1069 = t975 * t1032;
t1101 = t1069 + t1043;
t1071 = t975 * t1030;
t1073 = t603 * t1031;
t1100 = t1073 + t1071;
t1128 = -t1100 / 0.2e1;
t934 = t601 * t606;
t555 = pkin(3) * t934 + t606 * pkin(7);
t476 = pkin(4) * t526 + t555;
t331 = -pkin(5) * t717 + t476;
t587 = -pkin(3) * t602 - pkin(2);
t518 = -pkin(4) * t674 + t587;
t372 = -pkin(5) * t1031 + t518;
t993 = -t372 / 0.2e1;
t1094 = t1101 * t993 + t1128 * t331;
t608 = cos(qJ(2));
t977 = t608 / 0.2e1;
t99 = t1153 + t1141;
t612 = -t977 * t99 + t1094;
t1140 = t975 * t1122;
t746 = t1140 / 0.2e1;
t747 = -t1140 / 0.2e1;
t1162 = t747 + t746;
t1168 = qJD(2) * t1162;
t1150 = -t603 * t1032 + t975 * t717;
t1154 = t603 * t248;
t1158 = t1140 - t1154;
t978 = -t608 / 0.2e1;
t1117 = t1150 * t993 + t1158 * t978;
t1149 = -t603 * t1030 + t975 * t1031;
t1145 = -t1149 / 0.2e1;
t611 = t1145 * t331 + t1117;
t1133 = -t1100 ^ 2 + t1149 ^ 2;
t1102 = qJD(2) * t1133;
t998 = t1149 / 0.2e1;
t1157 = 0.2e1 * t1128 * t1101 - t1150 * (t1145 - t998);
t1167 = qJD(1) * t1157 + t1102;
t1134 = -t1101 ^ 2 + t1150 ^ 2;
t1105 = qJD(1) * t1134;
t1166 = qJD(2) * t1157 + t1105;
t1023 = qJD(3) * t1149;
t1165 = qJD(6) * t1162 - t1023;
t1108 = t1154 / 0.2e1 + t747;
t1164 = 0.2e1 * t1108;
t1163 = t1153 / 0.2e1;
t1161 = t1101 * t998;
t801 = t1159 * t1149;
t1021 = qJD(3) * t1150;
t1160 = qJD(4) * t1150;
t1156 = qJD(1) * t1150 + qJD(2) * t1149;
t850 = qJD(2) * t1100;
t1083 = t1149 * t850;
t1131 = t1161 + t1150 * t1100 / 0.2e1;
t1148 = qJD(1) * t1131 + t1083;
t864 = qJD(1) * t1101;
t1087 = t1150 * t864;
t1147 = qJD(2) * t1131 + t1087;
t812 = t608 * qJD(1);
t1144 = t1150 * t812;
t1085 = t331 * t1150;
t1142 = t1149 * t372;
t395 = -t1069 / 0.2e1;
t745 = t1069 / 0.2e1;
t1126 = t395 + t745;
t1139 = qJD(1) * t1126;
t440 = -t1071 / 0.2e1;
t1011 = -t1073 / 0.2e1 + t440;
t744 = t1071 / 0.2e1;
t1082 = t744 + t1073 / 0.2e1;
t1029 = -t1011 + t1082;
t1078 = t1100 * t608;
t528 = t552 * t608;
t510 = t607 * t528;
t530 = t674 * t608;
t901 = t604 * t530;
t420 = t510 + t901;
t402 = t975 * t420;
t894 = t607 * t530;
t902 = t604 * t528;
t423 = t894 - t902;
t914 = t603 * t423;
t736 = -t402 / 0.2e1 - t914 / 0.2e1;
t1099 = -t1078 / 0.2e1 + t736;
t1120 = qJD(1) * t1099;
t1138 = -qJD(4) * t1029 - t1120;
t1028 = t1011 + t1082;
t1137 = qJD(4) * t1028 + t1120;
t1084 = t1149 * t608;
t785 = t975 * t423;
t917 = t603 * t420;
t879 = -t917 / 0.2e1 + t785 / 0.2e1;
t1097 = t1084 / 0.2e1 + t879;
t1121 = qJD(1) * t1097;
t1136 = -t1121 + t801;
t1132 = qJD(4) * t1164;
t1127 = t300 * t977;
t696 = t850 + t864;
t1125 = t1100 * t372;
t1077 = t1101 * t331;
t260 = t1101 * t812;
t1092 = t1031 / 0.2e1;
t986 = t1030 / 0.2e1;
t1093 = t1032 * t1092 + t717 * t986;
t847 = qJD(2) * t1030;
t1116 = qJD(1) * t1093 + t1031 * t847;
t1096 = -t1084 / 0.2e1 + t879;
t1115 = qJD(2) * t1096 + t1160;
t854 = qJD(1) * t1032;
t1114 = qJD(2) * t1093 + t717 * t854;
t1113 = qJD(2) * t1099 + qJD(4) * t1126 - t260;
t1112 = qJD(2) * t1097 + t1144;
t1027 = 0.2e1 * t745 + t1043;
t1098 = t1078 / 0.2e1 + t736;
t1111 = qJD(2) * t1098 - qJD(4) * t1027;
t1110 = qJD(3) * t1028 + qJD(4) * t1170;
t1109 = t1118 * t300;
t1058 = -t1032 ^ 2 + t717 ^ 2;
t1104 = qJD(1) * t1058;
t1055 = -t1030 ^ 2 + t1031 ^ 2;
t1103 = qJD(2) * t1055;
t991 = -t1041 / 0.2e1;
t970 = pkin(5) * t1030;
t971 = pkin(5) * t1032;
t1088 = pkin(10) * t1032;
t880 = t1118 * t1031;
t690 = t847 + t854;
t1081 = qJD(5) * t1150;
t1079 = t1032 * t476;
t399 = t1032 * t812;
t1075 = t1150 * qJD(6);
t1068 = qJD(3) * t1126;
t1067 = qJD(3) * t1027;
t1065 = qJD(3) * t1029;
t1064 = qJD(3) * t1031;
t984 = -t518 / 0.2e1;
t985 = -t476 / 0.2e1;
t1053 = t1031 * t985 + t717 * t984;
t1052 = qJD(1) * t717 + qJD(2) * t1031;
t1050 = qJD(2) * t1028 + t1139;
t1049 = qJD(1) * t1027 + qJD(2) * t1029;
t982 = -t552 / 0.2e1;
t1047 = t587 / 0.2e1;
t1046 = pkin(10) * t717;
t1045 = t476 * t717;
t1040 = t717 * t812;
t1020 = qJD(3) * t717;
t1015 = qJD(4) * t717;
t506 = t603 * pkin(5);
t1014 = qJD(5) * t506;
t1013 = qJD(5) * t717;
t1012 = qJD(6) * t506;
t877 = -t902 / 0.2e1 + t894 / 0.2e1;
t732 = -pkin(2) * t608 - t606 * qJ(3);
t561 = -pkin(1) + t732;
t546 = t602 * t561;
t932 = t602 * t606;
t470 = -pkin(8) * t932 + t546 + (-pkin(7) * t601 - pkin(3)) * t608;
t931 = t602 * t608;
t584 = pkin(7) * t931;
t509 = t601 * t561 + t584;
t484 = -pkin(8) * t934 + t509;
t349 = -t976 * t470 + t605 * t484;
t303 = -t529 * pkin(9) - t349;
t968 = t608 * pkin(4);
t295 = t303 - t968;
t268 = t607 * t295;
t1003 = -t268 / 0.2e1;
t1001 = -t295 / 0.2e1;
t350 = t605 * t470 + t484 * t976;
t304 = -t526 * pkin(9) + t350;
t296 = t607 * t304;
t771 = -t296 / 0.2e1;
t973 = pkin(4) * t552;
t380 = t973 + t970;
t992 = -t380 / 0.2e1;
t990 = -t1032 / 0.2e1;
t987 = -t1030 / 0.2e1;
t770 = t511 / 0.2e1;
t769 = t539 / 0.2e1;
t983 = t674 / 0.2e1;
t981 = t555 / 0.2e1;
t597 = t608 * pkin(7);
t980 = t597 / 0.2e1;
t979 = -t606 / 0.2e1;
t974 = pkin(4) * t529;
t972 = pkin(4) * t607;
t969 = t606 * pkin(5);
t967 = t608 * pkin(5);
t965 = pkin(4) * qJD(4);
t964 = pkin(4) * qJD(5);
t279 = t402 + t914;
t933 = t601 * t608;
t556 = pkin(3) * t933 + t597;
t477 = pkin(4) * t528 + t556;
t332 = pkin(5) * t420 + t477;
t908 = t604 * t304;
t143 = -t268 + t908;
t680 = -t143 - t1088;
t116 = t680 - t967;
t797 = t975 * t116;
t144 = t604 * t295 + t296;
t125 = t144 + t1046;
t928 = t603 * t125;
t60 = -t797 + t928;
t566 = t606 * pkin(2) - qJ(3) * t608;
t583 = pkin(7) * t934;
t519 = t602 * t566 + t583;
t473 = t606 * pkin(3) - pkin(8) * t931 + t519;
t459 = t976 * t473;
t520 = -pkin(7) * t932 + t601 * t566;
t485 = -pkin(8) * t933 + t520;
t899 = t605 * t485;
t760 = t459 - t899;
t301 = t606 * pkin(4) - pkin(9) * t530 + t760;
t287 = t607 * t301;
t458 = t605 * t473;
t478 = t976 * t485;
t878 = t478 + t458;
t305 = -pkin(9) * t528 + t878;
t907 = t604 * t305;
t762 = t287 - t907;
t120 = -pkin(10) * t423 + t762 + t969;
t796 = t975 * t120;
t286 = t604 * t301;
t302 = t607 * t305;
t881 = t302 + t286;
t126 = -pkin(10) * t420 + t881;
t927 = t603 * t126;
t5 = (t796 - t927) * t608 + t60 * t606 + t332 * t1150 - t331 * t279;
t963 = t5 * qJD(1);
t283 = t785 - t917;
t795 = t975 * t125;
t930 = t603 * t116;
t61 = t795 + t930;
t794 = t975 * t126;
t929 = t603 * t120;
t6 = (t794 + t929) * t608 - t61 * t606 + t332 * t1101 + t331 * t283;
t962 = t6 * qJD(1);
t117 = -t1030 * t1032 + t1031 * t717;
t86 = t1030 * t990 + t1032 * t987 + 0.2e1 * t1092 * t717;
t961 = t86 * qJD(4) + t117 * qJD(5);
t909 = t604 * t303;
t149 = -t296 - t909;
t127 = t149 - t1046;
t793 = t975 * t127;
t896 = t607 * t303;
t150 = t896 - t908;
t128 = t150 - t1088;
t925 = t603 * t128;
t64 = t793 - t925;
t739 = t974 + t971;
t27 = t1150 * t739 + t64 * t608 - t1077;
t960 = qJD(1) * t27;
t792 = t975 * t128;
t926 = t603 * t127;
t65 = t792 + t926;
t28 = t1101 * t739 + t65 * t608 + t1085;
t959 = qJD(1) * t28;
t647 = t603 * t680;
t62 = -t647 - t795;
t29 = t1150 * t971 + t608 * t62 - t1077;
t958 = qJD(1) * t29;
t631 = t975 * t680;
t63 = t631 - t928;
t30 = t1101 * t971 + t608 * t63 + t1085;
t957 = qJD(1) * t30;
t37 = -t60 * t608 + t1085;
t956 = qJD(1) * t37;
t38 = -t608 * t61 - t1077;
t955 = qJD(1) * t38;
t83 = t149 * t608 + t717 * t974 - t1079;
t954 = qJD(1) * t83;
t84 = t1032 * t974 + t150 * t608 + t1045;
t953 = qJD(1) * t84;
t90 = -t143 * t608 + t1045;
t952 = qJD(1) * t90;
t91 = -t144 * t608 - t1079;
t951 = qJD(1) * t91;
t43 = t143 * t606 - t476 * t420 + t477 * t717 + t608 * t762;
t946 = t43 * qJD(1);
t44 = t1032 * t477 - t144 * t606 + t476 * t423 + t608 * t881;
t945 = t44 * qJD(1);
t944 = t1031 * t608;
t943 = t1030 * t608;
t942 = t529 * t1031;
t941 = t529 * t1030;
t800 = pkin(5) + t972;
t567 = t975 * t800;
t585 = t603 * t604 * pkin(4);
t536 = t585 - t567;
t940 = t536 * t608;
t740 = t603 * t800;
t782 = t975 * t604;
t537 = pkin(4) * t782 + t740;
t939 = t537 * t608;
t910 = t603 * t607;
t543 = (t782 + t910) * pkin(4);
t938 = t543 * t608;
t758 = t975 * t972;
t544 = t758 - t585;
t937 = t544 * t608;
t936 = t552 * t717;
t935 = t552 * t1032;
t76 = -t1101 * t279 + t1150 * t283;
t892 = t76 * qJD(1);
t752 = -t795 / 0.2e1;
t891 = -t930 / 0.2e1 + t752;
t774 = t928 / 0.2e1;
t890 = -t797 / 0.2e1 + t774;
t889 = -t647 / 0.2e1 + t752;
t888 = t774 - t631 / 0.2e1;
t887 = t1118 * t1093;
t579 = t601 ^ 2 + t602 ^ 2;
t132 = -t1150 * t606 - t279 * t608;
t876 = qJD(1) * t132;
t133 = -t1101 * t606 + t283 * t608;
t875 = qJD(1) * t133;
t778 = t1100 * t978;
t154 = t778 - t736;
t874 = qJD(1) * t154;
t156 = t778 + t736;
t872 = qJD(1) * t156;
t777 = t1149 * t977;
t163 = t777 + t879;
t868 = qJD(1) * t163;
t164 = t777 - t879;
t867 = qJD(1) * t164;
t256 = -t349 * t608 - t555 * t526;
t866 = qJD(1) * t256;
t257 = -t350 * t608 - t555 * t529;
t865 = qJD(1) * t257;
t284 = -t420 * t608 - t606 * t717;
t863 = qJD(1) * t284;
t285 = -t1032 * t606 + t423 * t608;
t862 = qJD(1) * t285;
t733 = -t510 / 0.2e1 - t901 / 0.2e1;
t776 = t1030 * t978;
t307 = t776 - t733;
t861 = qJD(1) * t307;
t308 = t776 + t733;
t860 = qJD(1) * t308;
t309 = -t943 / 0.2e1 + t733;
t859 = qJD(1) * t309;
t310 = t944 / 0.2e1 + t877;
t858 = qJD(1) * t310;
t775 = t1031 * t977;
t312 = t775 + t877;
t857 = qJD(1) * t312;
t313 = t775 - t877;
t856 = qJD(1) * t313;
t508 = -pkin(7) * t933 + t546;
t371 = (t508 * t602 + t509 * t601) * t606;
t855 = qJD(1) * t371;
t424 = t526 * t606 - t528 * t608;
t853 = qJD(1) * t424;
t425 = -t529 * t606 + t530 * t608;
t852 = qJD(1) * t425;
t851 = qJD(1) * t529;
t849 = qJD(2) * t372;
t846 = qJD(2) * t518;
t845 = qJD(2) * t552;
t844 = qJD(2) * t587;
t843 = qJD(3) * t608;
t842 = qJD(4) * t1101;
t841 = qJD(4) * t1032;
t840 = qJD(4) * t529;
t839 = qJD(4) * t552;
t461 = t769 - t539 / 0.2e1;
t837 = qJD(5) * t461;
t835 = qJD(5) * t518;
t833 = qJD(6) * t372;
t130 = t349 * t606 - t556 * t526 - t555 * t528 + t608 * t760;
t832 = t130 * qJD(1);
t131 = -t350 * t606 + t556 * t529 + t555 * t530 + t608 * t878;
t831 = t131 * qJD(1);
t153 = -t1032 * t420 + t423 * t717;
t830 = t153 * qJD(1);
t622 = t1047 * t529 + t487 * t977 + t552 * t981;
t734 = t459 / 0.2e1 - t899 / 0.2e1;
t208 = -t622 + t734;
t829 = t208 * qJD(1);
t623 = t1047 * t526 + t486 * t977 - t674 * t981;
t766 = -t458 / 0.2e1 - t478 / 0.2e1;
t209 = t623 + t766;
t828 = t209 * qJD(1);
t253 = (t508 * t608 + t519 * t606) * t602 + (t509 * t608 + t520 * t606) * t601;
t827 = t253 * qJD(1);
t288 = pkin(7) ^ 2 * t606 * t608 + t508 * t519 + t509 * t520;
t826 = t288 * qJD(1);
t340 = -t526 * t530 - t528 * t529;
t825 = t340 * qJD(1);
t366 = -t508 * t606 + (t519 - 0.2e1 * t583) * t608;
t824 = t366 * qJD(1);
t367 = t520 * t608 + (-t509 + 0.2e1 * t584) * t606;
t823 = t367 * qJD(1);
t640 = -t897 / 0.2e1 - t799 / 0.2e1;
t447 = (t982 - t640) * t608;
t822 = t447 * qJD(1);
t448 = (t982 + t640) * t608;
t821 = t448 * qJD(1);
t449 = t798 * t977 + (-t674 + t898) * t978;
t820 = t449 * qJD(1);
t450 = (t983 - t798 / 0.2e1 + t898 / 0.2e1) * t608;
t819 = t450 * qJD(1);
t818 = t526 * qJD(4);
t600 = t606 ^ 2;
t549 = t579 * t600;
t817 = t549 * qJD(1);
t542 = t674 * qJD(4);
t816 = t579 * qJD(2);
t580 = t608 ^ 2 - t600;
t815 = t580 * qJD(1);
t814 = t606 * qJD(1);
t813 = t606 * qJD(2);
t811 = t608 * qJD(2);
t58 = -t1100 * t1101 + t1149 * t1150;
t810 = t58 * qJD(6) + t1118 * t1157;
t82 = -t1128 * t1150 + t1161;
t809 = t82 * qJD(6) + t1118 * t1131;
t808 = t975 / 0.2e1;
t807 = pkin(1) * t814;
t806 = pkin(1) * t812;
t805 = pkin(7) * t811;
t804 = t971 / 0.2e1;
t803 = t970 / 0.2e1;
t802 = t968 / 0.2e1;
t781 = t606 * t843;
t780 = t606 * t811;
t779 = t606 * t812;
t772 = t604 * t979;
t768 = t607 * t606 / 0.2e1;
t767 = -t286 / 0.2e1 - t302 / 0.2e1;
t765 = t975 * qJD(5);
t764 = t975 * qJD(6);
t763 = pkin(4) * t1118;
t505 = t529 * t812;
t759 = qJD(2) * t448 - t505;
t757 = -qJD(4) + t812;
t756 = -qJD(5) + t812;
t755 = -qJD(6) + t812;
t754 = qJD(3) + t844;
t753 = pkin(5) * t808;
t751 = t795 / 0.2e1;
t750 = -t794 / 0.2e1;
t749 = t793 / 0.2e1;
t748 = -t792 / 0.2e1;
t743 = t606 * t808;
t742 = t802 + t303 / 0.2e1;
t741 = t812 - qJD(4) / 0.2e1;
t738 = t287 / 0.2e1 - t907 / 0.2e1;
t731 = -t758 / 0.2e1;
t610 = -t1150 * t992 + t739 * t998 + t612;
t645 = -t927 / 0.2e1 + t796 / 0.2e1;
t618 = t536 * t979 + t645;
t1 = t610 + t618;
t100 = -t1149 * t380 + t1125;
t728 = -t1 * qJD(1) + t100 * qJD(2);
t101 = t1100 * t380 + t1142;
t609 = t1101 * t992 + t1128 * t739 + t611;
t646 = -t929 / 0.2e1 + t750;
t617 = t537 * t979 + t646;
t2 = t609 + t617;
t727 = -t2 * qJD(1) + t101 * qJD(2);
t105 = t1149 * t970 - t1125;
t8 = (-t1149 * t990 - t1150 * t987 + t743) * pkin(5) + t612 + t645;
t726 = -t8 * qJD(1) - t105 * qJD(2);
t106 = -t1100 * t970 - t1142;
t7 = (t1100 * t990 + t1101 * t987 + t603 * t979) * pkin(5) + t611 + t646;
t725 = -t7 * qJD(1) - t106 * qJD(2);
t720 = qJD(2) * t58 + t1105;
t718 = -t519 * t601 + t520 * t602;
t716 = qJD(1) * t58 + t1102;
t70 = t771 + t296 / 0.2e1 + (t1001 + t742) * t604;
t715 = t70 * qJD(1);
t146 = t991 + t1041 / 0.2e1;
t72 = t607 * t742 + t1003;
t714 = t72 * qJD(1) + t146 * qJD(2);
t713 = qJD(2) * t86 + t1104;
t712 = qJD(1) * t86 + t1103;
t254 = -t1030 * t518 + t1031 * t973;
t615 = t1030 * t985 + t1032 * t984 - t1127;
t40 = (t768 + t936 / 0.2e1 + t942 / 0.2e1) * pkin(4) + t615 + t738;
t711 = -t40 * qJD(1) - t254 * qJD(2);
t255 = -t1030 * t973 - t1031 * t518;
t613 = t719 * t978 + t1053;
t39 = (t772 - t935 / 0.2e1 - t941 / 0.2e1) * pkin(4) + t613 + t767;
t710 = -t39 * qJD(1) - t255 * qJD(2);
t709 = qJD(2) * t117 + t1104;
t708 = qJD(1) * t117 + t1103;
t261 = -t526 * t674 - t552 * t529;
t370 = t526 ^ 2 - t529 ^ 2;
t699 = qJD(1) * t370 + qJD(2) * t261;
t429 = -t552 ^ 2 + t674 ^ 2;
t698 = qJD(1) * t261 + qJD(2) * t429;
t339 = t440 + t744;
t697 = qJD(2) * t339 + t1139;
t357 = 0.2e1 * t770 - t903;
t379 = 0.2e1 * t769 + t900;
t693 = qJD(1) * t357 + qJD(2) * t379;
t677 = t508 * t601 / 0.2e1 - t509 * t602 / 0.2e1;
t368 = t980 + t677;
t557 = t579 * qJ(3);
t692 = qJD(1) * t368 - qJD(2) * t557;
t415 = t770 - t511 / 0.2e1;
t691 = qJD(1) * t415 + qJD(2) * t461;
t688 = -qJD(1) * t526 + qJD(2) * t674;
t687 = t845 + t851;
t686 = qJD(4) * t1030 + qJD(5) * t379;
t685 = -t1013 - t1015;
t684 = t1060 * t1028;
t683 = t567 / 0.2e1 + t753;
t682 = t603 * t967 / 0.2e1 + t891;
t678 = -qJD(5) / 0.2e1 + t741;
t626 = t978 * t99 + t1094;
t11 = t626 + t645;
t671 = qJD(1) * t11 - t1100 * t849;
t627 = t331 * t998 - t1117;
t12 = -t627 + t646;
t670 = qJD(1) * t12 - t1149 * t849;
t669 = -qJD(2) * t82 - t1087;
t666 = -qJD(1) * t82 - t1083;
t625 = -t719 * t977 + t1053;
t67 = t625 + t767;
t663 = qJD(1) * t67 - t1031 * t846;
t624 = t1127 + t476 * t986 + t518 * t1032 / 0.2e1;
t66 = -t624 + t738;
t662 = qJD(1) * t66 - t1030 * t846;
t362 = t526 * t982 + t529 * t983;
t657 = qJD(2) * t362 - t526 * t851;
t656 = -qJD(1) * t362 - t674 * t845;
t655 = qJD(2) * t447 - t505 + t840;
t654 = qJD(2) * t308 - qJD(5) * t415 - t399;
t648 = t608 * t753 + t890;
t644 = -t926 / 0.2e1 + t748;
t643 = -t925 / 0.2e1 + t749;
t15 = t749 - t939 / 0.2e1 + t751 + (-t128 / 0.2e1 + t116 / 0.2e1) * t603;
t48 = t1163 + t1141 / 0.2e1 + t1169;
t639 = qJD(1) * t15 + qJD(2) * t48 + qJD(4) * t537;
t621 = (-t127 / 0.2e1 - t125 / 0.2e1) * t603 + t748;
t16 = t940 / 0.2e1 + t797 / 0.2e1 + t621;
t638 = qJD(1) * t16 - qJD(4) * t536 - t1168;
t614 = t647 / 0.2e1 + t751;
t19 = t614 + t682;
t637 = t19 * qJD(1) - t506 * qJD(4);
t620 = t631 / 0.2e1;
t21 = t620 - t928 / 0.2e1 + t648;
t507 = t731 + t683;
t636 = -t21 * qJD(1) + t507 * qJD(4) - t1168;
t23 = -t938 / 0.2e1 + t614 + t643;
t45 = t240 + t1163;
t635 = qJD(1) * t23 - qJD(2) * t45 + qJD(4) * t543;
t24 = -t937 / 0.2e1 + t620 + t621;
t53 = t746 - t1154 / 0.2e1 + t1108;
t634 = qJD(1) * t24 - qJD(2) * t53 + qJD(4) * t544;
t633 = qJD(4) * t1100 + t1029 * t1060;
t632 = -t1075 - t1081 - t1160;
t629 = qJD(2) * t307 + qJD(5) * t357 - t399 + t841;
t628 = qJD(2) * t156 - t1060 * t1126 - t260;
t619 = qJD(2) * t732 + t843;
t616 = qJD(2) * t154 + t1027 * t1060 - t260 + t842;
t592 = -t814 / 0.2e1;
t591 = t814 / 0.2e1;
t590 = t813 / 0.2e1;
t547 = t741 * t606;
t541 = t544 * qJD(5);
t540 = t543 * qJD(5);
t525 = t537 * qJD(6);
t524 = t536 * qJD(6);
t523 = t678 * t606;
t495 = (-qJD(6) / 0.2e1 + t678) * t606;
t483 = t585 + t731 - t683;
t482 = -t506 / 0.2e1 - t740 / 0.2e1 + (-t782 - t910 / 0.2e1) * pkin(4);
t369 = t980 - t677;
t363 = t449 * qJD(2) - t526 * t812;
t351 = t362 * qJD(4);
t314 = t943 / 0.2e1 + t733;
t311 = -t944 / 0.2e1 + t877;
t306 = -t450 * qJD(2) + t526 * t757;
t259 = t261 * qJD(4);
t252 = -t1073 + 0.2e1 * t440;
t223 = -t1043 + 0.2e1 * t395;
t211 = t622 + t734;
t210 = -t623 + t766;
t151 = t312 * qJD(2) + t1040;
t147 = t904 + 0.2e1 * t991;
t102 = -t313 * qJD(2) - t717 * t757 + t1013;
t73 = t607 * t802 + t908 + t1003 - t896 / 0.2e1;
t71 = 0.2e1 * t771 - t909 / 0.2e1 + (t802 + t1001) * t604;
t69 = t624 + t738;
t68 = -t625 + t767;
t42 = pkin(4) * t772 - t613 + t767 + (t935 + t941) * pkin(4) / 0.2e1;
t41 = pkin(4) * t768 - t615 + t738 - (t936 + t942) * pkin(4) / 0.2e1;
t36 = t163 * qJD(2) + t1144;
t31 = -t164 * qJD(2) + (t1060 - t757) * t1150;
t26 = t937 / 0.2e1 + t644 + t888;
t25 = t938 / 0.2e1 + t643 + t889;
t22 = t648 + t888;
t20 = t682 + t889;
t18 = t939 / 0.2e1 + t643 + t891;
t17 = -t940 / 0.2e1 + t644 + t890;
t14 = -t626 + t645;
t13 = t627 + t646;
t10 = t1101 * t803 + t1100 * t804 + t750 + (-t969 / 0.2e1 - t120 / 0.2e1) * t603 - t611;
t9 = pkin(5) * t743 - t1149 * t804 - t1150 * t803 - t612 + t645;
t4 = -t609 + t617;
t3 = -t610 + t618;
t32 = [0, 0, 0, t780, t580 * qJD(2), 0, 0, 0, -pkin(1) * t813, -pkin(1) * t811, -t366 * qJD(2) + t602 * t781, t367 * qJD(2) - t601 * t781, -qJD(2) * t253 + qJD(3) * t549, qJD(2) * t288 - qJD(3) * t371 (qJD(2) * t530 - t818) * t529, qJD(2) * t340 + qJD(4) * t370, -t425 * qJD(2) + t608 * t818, -t424 * qJD(2) + t608 * t840, -t780, -t130 * qJD(2) - t257 * qJD(4) + t529 * t843, t131 * qJD(2) + t256 * qJD(4) - t526 * t843 (qJD(2) * t423 - t685) * t1032, qJD(2) * t153 + t1058 * t1118, -t285 * qJD(2) + t608 * t685, -t284 * qJD(2) + (qJD(5) * t1032 + t841) * t608, -t780, -t43 * qJD(2) - t83 * qJD(4) - t91 * qJD(5) + t1032 * t843, t44 * qJD(2) + t84 * qJD(4) + t90 * qJD(5) + t717 * t843 (qJD(2) * t283 - t632) * t1101, qJD(2) * t76 + t1134 * t1159, -t133 * qJD(2) + t608 * t632, -t132 * qJD(2) + (t1060 * t1101 + t842) * t608, -t780, -t5 * qJD(2) - t27 * qJD(4) - t29 * qJD(5) - t38 * qJD(6) + t1101 * t843, t6 * qJD(2) + t28 * qJD(4) + t30 * qJD(5) + t37 * qJD(6) + t1150 * t843; 0, 0, 0, t779, t815, t811, -t813, 0, -t805 - t807, pkin(7) * t813 - t806, -t584 * qJD(2) + t601 * t619 - t824, t601 * t805 + t602 * t619 + t823, qJD(2) * t718 - t827, t826 + (-pkin(2) * t597 + qJ(3) * t718) * qJD(2) + t369 * qJD(3), t530 * t687 + t351, t825 + (-t528 * t552 + t530 * t674) * qJD(2) + t259, -qJD(4) * t450 + t552 * t813 - t852, -qJD(4) * t447 + t674 * t813 - t853, -t547, -t832 + (-t486 * t606 + t528 * t587 - t556 * t674) * qJD(2) - t448 * qJD(3) + t211 * qJD(4), t831 + (-t487 * t606 + t530 * t587 + t552 * t556) * qJD(2) + t449 * qJD(3) + t210 * qJD(4), t423 * t690 + t887, t830 + (-t1030 * t420 + t1031 * t423) * qJD(2) + t961, -qJD(4) * t313 + qJD(5) * t311 + t1030 * t813 - t862, -qJD(4) * t307 + qJD(5) * t314 + t1031 * t813 - t863, -t523, -t946 + (-t1031 * t477 + t420 * t518 + t606 * t719) * qJD(2) - t308 * qJD(3) + t41 * qJD(4) + t69 * qJD(5), t945 + (t1030 * t477 - t300 * t606 + t423 * t518) * qJD(2) + t312 * qJD(3) + t42 * qJD(4) + t68 * qJD(5), t283 * t696 + t809, t892 + (-t1100 * t279 + t1149 * t283) * qJD(2) + t810, -qJD(4) * t164 + t1060 * t1096 + t1100 * t813 - t875, -qJD(4) * t154 + t1060 * t1098 + t1149 * t813 - t876, -t495, -t963 + (-t1149 * t332 + t1158 * t606 + t279 * t372) * qJD(2) - t156 * qJD(3) + t3 * qJD(4) + t9 * qJD(5) + t14 * qJD(6), t962 + (t1100 * t332 + t283 * t372 - t606 * t99) * qJD(2) + t163 * qJD(3) + t4 * qJD(4) + t10 * qJD(5) + t13 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (qJD(2) * t601 + t602 * t814) * t608 (qJD(2) * t602 - t601 * t814) * t608, t817, qJD(2) * t369 - t855, 0, 0, 0, 0, 0, -t759, t363, 0, 0, 0, 0, 0, -t654, t151, 0, 0, 0, 0, 0, -t628, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t657, t699, t306, -t655, t590, qJD(2) * t211 - qJD(4) * t350 - t865, qJD(2) * t210 + qJD(4) * t349 + t866, t1114, t713, t102, -t629, t590, qJD(2) * t41 + qJD(4) * t149 + qJD(5) * t71 - t954, qJD(2) * t42 - qJD(4) * t150 + qJD(5) * t73 + t953, t1147, t1166, t31, -t616, t590, qJD(2) * t3 + qJD(4) * t64 + qJD(5) * t25 + qJD(6) * t18 - t960, qJD(2) * t4 - qJD(4) * t65 + qJD(5) * t26 + qJD(6) * t17 + t959; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1114, t709, t311 * qJD(2) - t717 * t756 + t1015, t314 * qJD(2) - t357 * qJD(4) + t1032 * t756, t590, qJD(2) * t69 + qJD(3) * t415 + qJD(4) * t71 - qJD(5) * t144 - t951, qJD(2) * t68 + qJD(4) * t73 + qJD(5) * t143 + t952, t1147, t1166, -t1150 * t756 + t1075 + t1115, t223 * qJD(6) + t1101 * t756 + t1111, t590, qJD(2) * t9 + qJD(4) * t25 + qJD(5) * t62 + qJD(6) * t20 + t1068 - t958, qJD(2) * t10 + qJD(4) * t26 - qJD(5) * t63 + qJD(6) * t22 + t957; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t669, t720, -t1150 * t755 + t1081 + t1115, t223 * qJD(5) + t1101 * t755 + t1111, t590, qJD(2) * t14 + qJD(4) * t18 + qJD(5) * t20 - qJD(6) * t61 + t1068 - t955, qJD(2) * t13 + qJD(4) * t17 + qJD(5) * t22 + qJD(6) * t60 + t956; 0, 0, 0, -t779, -t815, 0, 0, 0, t807, t806, t824, -t823, t827, -qJD(3) * t368 - t826, -t530 * t851 + t351, t259 - t825, -qJD(4) * t449 + t852, -qJD(4) * t448 + t853, t547, -qJD(3) * t447 - qJD(4) * t208 + t832, qJD(3) * t450 - qJD(4) * t209 - t831, -t423 * t854 + t887, -t830 + t961, -qJD(4) * t312 - qJD(5) * t310 + t862, -qJD(4) * t308 - qJD(5) * t309 + t863, t523, -qJD(3) * t307 - qJD(4) * t40 - qJD(5) * t66 + t946, qJD(3) * t313 - qJD(4) * t39 - qJD(5) * t67 - t945, -t283 * t864 + t809, t810 - t892, -qJD(4) * t163 - t1060 * t1097 + t875, -qJD(4) * t156 - t1060 * t1099 + t876, t495, -qJD(3) * t154 - qJD(4) * t1 - qJD(5) * t8 - qJD(6) * t11 + t963, qJD(3) * t164 - qJD(4) * t2 - qJD(5) * t7 - qJD(6) * t12 - t962; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t579 * qJD(3), t557 * qJD(3), t552 * t542, t429 * qJD(4), 0, 0, 0, t587 * t839, t587 * t542, t880 * t1030, t1118 * t1055, 0, 0, 0, -qJD(4) * t254 + t1030 * t835, -qJD(4) * t255 + t1031 * t835, t801 * t1100, t1159 * t1133, 0, 0, 0, qJD(4) * t100 - qJD(5) * t105 + t1100 * t833, qJD(4) * t101 - qJD(5) * t106 + t1149 * t833; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t816, -t692, 0, 0, 0, 0, 0, -t822, t819, 0, 0, 0, 0, 0, t837 - t861, t856, 0, 0, 0, 0, 0, t684 - t874, t867; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t656, t698, t542 - t820, -t821 - t839, t592, -qJD(4) * t487 + t552 * t844 - t829, qJD(4) * t486 + t674 * t844 - t828, t1116, t712, -t857 + t880, -t686 - t860, t592, -t1109 + t711, -qJD(4) * t719 + t147 * qJD(5) + t710, t1148, t1167, t801 - t868, -t633 - t872, t592, -qJD(4) * t99 + t1060 * t1170 + t728, -qJD(4) * t1158 + t1060 * t1164 + t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1116, t708, -t858 + t880, -qJD(4) * t379 - qJD(5) * t1030 - t859, t592, qJD(3) * t461 - t1109 - t662, qJD(4) * t147 - qJD(5) * t719 - t663, t1148, t1167, t1136, -qJD(5) * t1100 + qJD(6) * t252 + t1138, t592, -qJD(5) * t99 + qJD(6) * t1170 + t1110 + t726, -qJD(5) * t1158 + qJD(6) * t1164 + t1132 + t725; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t666, t716, t1136, qJD(5) * t252 - qJD(6) * t1100 + t1138, t592, qJD(5) * t1170 - qJD(6) * t99 + t1110 - t671, qJD(5) * t1164 - qJD(6) * t1158 + t1132 - t670; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t602 * t779, t601 * t779, -t817, qJD(2) * t368 + t855, 0, 0, 0, 0, 0, t655, t306, 0, 0, 0, 0, 0, t629, t102, 0, 0, 0, 0, 0, t616, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t816, t692, 0, 0, 0, 0, 0, t822 + t839, t542 - t819, 0, 0, 0, 0, 0, t686 + t861, -t856 + t880, 0, 0, 0, 0, 0, t633 + t874, t801 - t867; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t687, t688, 0, 0, 0, 0, 0, t690, t1052, 0, 0, 0, 0, 0, t696, t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t693, t1052, 0, 0, 0, 0, 0, t1049, t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1049, t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t657, -t699, t363, t759, t590, qJD(2) * t208 - qJD(3) * t529 + t865, qJD(2) * t209 + qJD(3) * t526 - t866, -t1114, -t713, t151, t654, t590, qJD(2) * t40 - qJD(3) * t1032 + qJD(5) * t70 + t954, qJD(2) * t39 + qJD(5) * t72 - t1020 - t953, -t1147, -t1166, t36, t628, t590, qJD(2) * t1 - qJD(3) * t1101 - qJD(5) * t23 - qJD(6) * t15 + t960, qJD(2) * t2 - qJD(5) * t24 - qJD(6) * t16 - t1021 - t959; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t656, -t698, t820, t821, t591, -t552 * t754 + t829, -t674 * t754 + t828, -t1116, -t712, t857, -t837 + t860, t591, -qJD(3) * t1030 - t711, qJD(5) * t146 - t1064 - t710, -t1148, -t1167, t868, -t684 + t872, t591, -qJD(3) * t1100 + qJD(5) * t45 - qJD(6) * t48 - t728, qJD(5) * t53 + t1165 - t727; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t687, -t688, 0, 0, 0, 0, 0, -t690, -t1052, 0, 0, 0, 0, 0, -t696, -t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t604 * t964, -t607 * t964, 0, 0, 0, 0, 0, -t540 - t525, -t541 + t524; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t691, 0, -t604 * t763 + t715, -t607 * t763 + t714, 0, 0, 0, -t1050, 0, qJD(6) * t482 - t540 - t635, qJD(6) * t483 - t541 - t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1050, 0, qJD(5) * t482 - t525 - t639, qJD(5) * t483 + t524 - t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1114, -t709, t310 * qJD(2) + t1040, t309 * qJD(2) + t415 * qJD(4) - t399, t590, qJD(2) * t66 - qJD(3) * t357 - qJD(4) * t70 + t951, qJD(2) * t67 - qJD(4) * t72 - t1020 - t952, -t1147, -t1166, t1112, qJD(6) * t1126 + t1113, t590, qJD(2) * t8 + qJD(4) * t23 + qJD(6) * t19 - t1067 + t958, qJD(2) * t7 + qJD(4) * t24 + qJD(6) * t21 - t1021 - t957; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1116, -t708, t858, qJD(4) * t461 + t859, t591, -qJD(3) * t379 + t662, -qJD(4) * t146 - t1064 + t663, -t1148, -t1167, t1121, qJD(6) * t339 + t1137, t591, -qJD(4) * t45 - t1065 - t726, -qJD(4) * t53 + t1165 - t725; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t693, -t1052, 0, 0, 0, 0, 0, -t1049, -t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t691, 0, t604 * t965 - t715, t607 * t965 - t714, 0, 0, 0, t1050, 0, t635 - t1012, -qJD(6) * t507 + t634; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1012, -pkin(5) * t764; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t697, 0, -t1060 * t506 + t637 (-t765 - t764) * pkin(5) - t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t669, -t720, t1112, -qJD(5) * t1126 + t1113, t590, qJD(2) * t11 + qJD(4) * t15 - qJD(5) * t19 - t1067 + t955, qJD(2) * t12 + qJD(4) * t16 - qJD(5) * t21 - t1021 - t956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t666, -t716, t1121, -qJD(5) * t339 + t1137, t591, qJD(4) * t48 - t1065 + t671, -t1118 * t1162 - t1023 + t670; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1049, -t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1050, 0, t639 + t1014, qJD(5) * t507 + t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t697, 0, -t637 + t1014, pkin(5) * t765 + t636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t32;
