% Calculate minimal parameter regressor of coriolis matrix for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x33]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRRRPR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:24
% EndTime: 2019-03-09 22:10:17
% DurationCPUTime: 35.86s
% Computational Cost: add. (29346->880), mult. (59424->1085), div. (0->0), fcn. (69505->10), ass. (0->692)
t1050 = sin(qJ(3));
t1051 = cos(qJ(3));
t681 = sin(qJ(2));
t684 = cos(qJ(2));
t639 = t1050 * t681 - t1051 * t684;
t1061 = t639 / 0.2e1;
t678 = cos(pkin(11));
t683 = cos(qJ(4));
t677 = sin(pkin(11));
t680 = sin(qJ(4));
t970 = t677 * t680;
t630 = t678 * t683 - t970;
t1044 = t630 * pkin(5);
t1039 = t683 * pkin(4);
t670 = -pkin(3) - t1039;
t570 = t670 - t1044;
t1066 = -t570 / 0.2e1;
t665 = t1050 * t684;
t666 = t1051 * t681;
t641 = -t666 - t665;
t926 = t683 * t641;
t966 = t678 * t680;
t468 = t641 * t966 + t677 * t926;
t470 = t641 * t970 - t678 * t926;
t679 = sin(qJ(6));
t682 = cos(qJ(6));
t315 = -t682 * t468 + t679 * t470;
t1038 = -qJ(5) - pkin(9);
t650 = t1038 * t680;
t651 = t1038 * t683;
t1114 = t678 * t650 + t677 * t651;
t631 = t677 * t683 + t966;
t618 = t631 * pkin(10);
t1137 = t1114 - t618;
t1169 = t682 * t1137;
t540 = t677 * t650 - t678 * t651;
t617 = t630 * pkin(10);
t453 = t617 + t540;
t1187 = t679 * t453;
t779 = t1169 - t1187;
t702 = t1061 * t779 - t1066 * t315;
t1085 = -pkin(8) - pkin(7);
t654 = t1085 * t681;
t655 = t1085 * t684;
t793 = -t1050 * t655 - t1051 * t654;
t1151 = t793 * t680;
t503 = t683 * t639;
t1168 = -t641 * pkin(4) + qJ(5) * t503 + t1151;
t1041 = t641 * pkin(3);
t1042 = t639 * pkin(9);
t529 = -t1041 + t1042;
t511 = t683 * t529;
t299 = t511 + t1168;
t498 = t680 * t639;
t856 = qJ(5) * t498;
t1150 = t793 * t683;
t510 = t680 * t529;
t916 = t1150 - t510;
t330 = t856 - t916;
t177 = t677 * t299 + t678 * t330;
t467 = t631 * t639;
t459 = t467 * pkin(10);
t130 = t459 + t177;
t940 = t682 * t130;
t176 = t678 * t299 - t330 * t677;
t469 = t498 * t677 - t503 * t678;
t814 = -t641 * pkin(5) - pkin(10) * t469;
t107 = t176 + t814;
t965 = t679 * t107;
t749 = -t965 / 0.2e1 - t940 / 0.2e1;
t1197 = t702 + t749;
t862 = t1051 * pkin(2);
t669 = -t862 - pkin(3);
t649 = t669 - t1039;
t552 = t649 - t1044;
t1068 = -t552 / 0.2e1;
t861 = t1050 * pkin(2);
t668 = t861 + pkin(9);
t920 = qJ(5) + t668;
t619 = t920 * t680;
t620 = t920 * t683;
t1116 = -t678 * t619 - t677 * t620;
t1136 = t1116 - t618;
t482 = -t677 * t619 + t678 * t620;
t398 = t617 + t482;
t781 = t682 * t1136 - t679 * t398;
t704 = t1061 * t781 - t1068 * t315;
t1049 = pkin(2) * t681;
t506 = t529 + t1049;
t486 = t683 * t506;
t294 = t486 + t1168;
t485 = t680 * t506;
t917 = -t485 + t1150;
t324 = t856 - t917;
t168 = t678 * t294 - t324 * t677;
t95 = t168 + t814;
t1032 = t679 * t95;
t169 = t677 * t294 + t678 * t324;
t113 = t459 + t169;
t941 = t682 * t113;
t765 = -t1032 / 0.2e1 - t941 / 0.2e1;
t1196 = t704 + t765;
t1195 = t749 - t702;
t1194 = t765 - t704;
t1178 = qJD(4) + qJD(6);
t1062 = -t639 / 0.2e1;
t265 = t679 * t1136 + t682 * t398;
t1193 = t1062 * t265;
t301 = t679 * t1137 + t682 * t453;
t1192 = t1062 * t301;
t1191 = t1178 * t265;
t1190 = t1178 * t301;
t1189 = t1178 * t781;
t495 = -t682 * t630 + t679 * t631;
t1153 = t495 * t639;
t930 = t682 * t469;
t951 = t679 * t467;
t747 = t951 / 0.2e1 + t930 / 0.2e1;
t1163 = t1153 / 0.2e1 + t747;
t1179 = qJD(1) * t1163;
t918 = t1178 * t495;
t1184 = -t1179 - t918;
t1083 = -t1169 / 0.2e1;
t1077 = -t482 / 0.2e1;
t1176 = t482 / 0.2e1;
t827 = t1176 + t1077;
t1183 = t470 * t827;
t643 = t1050 * t654;
t645 = t1051 * t655;
t915 = -t645 + t643;
t1119 = t915 * t680;
t671 = -pkin(2) * t684 - pkin(1);
t1043 = t639 * pkin(3);
t796 = t641 * pkin(9) + t1043;
t729 = t671 + t796;
t350 = -t683 * t729 + t1119;
t320 = qJ(5) * t926 - t350;
t293 = pkin(4) * t639 + t320;
t1118 = t915 * t683;
t351 = t680 * t729 + t1118;
t971 = t641 * t680;
t321 = qJ(5) * t971 + t351;
t969 = t678 * t321;
t167 = t677 * t293 + t969;
t182 = -t320 * t677 - t969;
t1149 = t167 + t182;
t778 = t630 * t679 + t682 * t631;
t1074 = -t778 / 0.2e1;
t1075 = -t495 / 0.2e1;
t456 = t682 * t470;
t950 = t679 * t468;
t1117 = t456 + t950;
t1160 = t1074 * t315 + t1075 * t1117;
t866 = qJD(2) + qJD(3);
t1182 = t1160 * t866;
t1181 = t1163 * t866;
t1164 = -t1153 / 0.2e1 + t747;
t1180 = t1164 * t866;
t1120 = t866 * t778;
t1177 = -qJD(1) * t1160 + t1120 * t495;
t840 = t1150 / 0.2e1;
t1132 = t495 ^ 2 - t778 ^ 2;
t1175 = t1132 * t866;
t609 = pkin(4) * t498;
t1115 = -t609 + t915;
t1133 = -t467 * pkin(5) + t1115;
t1174 = t1133 * t495;
t1173 = t1133 * t778;
t877 = t315 * qJD(1);
t1167 = t495 * t866 + t877;
t1131 = -t1117 ^ 2 + t315 ^ 2;
t1166 = qJD(1) * t1131;
t1069 = -t540 / 0.2e1;
t1073 = t540 / 0.2e1;
t1165 = t1069 + t1073;
t1134 = t1117 * qJD(1);
t1135 = t1134 + t1120;
t454 = t682 * t467;
t949 = t679 * t469;
t314 = -t454 + t949;
t440 = -pkin(4) * t971 + t793;
t328 = -t468 * pkin(5) + t440;
t1045 = t470 * pkin(10);
t304 = t677 * t321;
t166 = t678 * t293 - t304;
t708 = t639 * pkin(5) - t1045 + t166;
t89 = t682 * t708;
t1047 = pkin(10) * t468;
t110 = t167 + t1047;
t964 = t679 * t110;
t65 = -t89 + t964;
t1162 = t1133 * t315 + t328 * t314 + t65 * t641;
t318 = t930 + t951;
t700 = t679 * t708;
t942 = t682 * t110;
t66 = t700 + t942;
t1161 = t1133 * t1117 + t328 * t318 + t66 * t641;
t1159 = t167 / 0.2e1;
t1158 = -t631 / 0.2e1;
t1004 = t328 * t495;
t251 = -t1004 / 0.2e1;
t1157 = -t1114 / 0.2e1;
t1156 = -t1116 / 0.2e1;
t1123 = t328 * t778;
t842 = -t1123 / 0.2e1;
t1155 = t315 * t328;
t1154 = t315 * t495;
t623 = t666 / 0.2e1 + t665 / 0.2e1;
t902 = qJD(1) * t641;
t845 = t639 * t902;
t429 = -t623 * qJD(4) + t845;
t1148 = qJD(5) * t315;
t1144 = t1115 * t440;
t1143 = t1115 * t649;
t1142 = t1117 * t328;
t310 = t315 * qJD(4);
t1140 = t495 * qJD(5);
t1121 = t866 * t641;
t1139 = t639 * t1121;
t1059 = -t645 / 0.2e1;
t1138 = t1059 - t609 / 0.2e1 + t643 / 0.2e1;
t1091 = t1178 * t778;
t853 = t683 * t1051;
t854 = t680 * t1051;
t572 = (-t677 * t854 + t678 * t853) * pkin(2);
t928 = t682 * t572;
t571 = (-t677 * t853 - t678 * t854) * pkin(2);
t948 = t679 * t571;
t746 = -t948 / 0.2e1 - t928 / 0.2e1;
t1065 = t570 / 0.2e1;
t1067 = t552 / 0.2e1;
t820 = t1065 + t1067;
t193 = t495 * t820 + t746;
t1130 = (t350 - t1119) * t641;
t1129 = (t351 - t1118) * t641;
t771 = -qJD(6) * t315 - t310;
t1128 = t866 * t793;
t1126 = 0.2e1 * t641;
t1125 = 0.2e1 * t680;
t1124 = t1117 * t778;
t1122 = t778 * t820;
t527 = t866 * t639;
t972 = t641 * t668;
t973 = t639 * t669;
t1113 = t1043 / 0.2e1 + t972 / 0.2e1 - t973 / 0.2e1 + (t1051 * t1062 - t1050 * t641 / 0.2e1) * pkin(2);
t365 = qJD(6) * t623 - t429;
t1110 = t778 * qJD(5);
t1057 = t669 / 0.2e1;
t803 = -t862 / 0.2e1;
t1109 = t1057 + t803;
t1106 = t166 * t1158 + t1159 * t630;
t1108 = -t1106 + t1138;
t1107 = t1106 + t1138;
t633 = t641 ^ 2;
t865 = -t639 ^ 2 + t633;
t1105 = -t1117 * t877 + t1182;
t1104 = -t1134 * t315 + t1182;
t81 = -t1124 + t1154;
t1103 = t81 * t866 + t1166;
t68 = t1154 / 0.2e1 + t1117 * t1074 - t315 * t1075 - t1124 / 0.2e1;
t1102 = t68 * t866 + t1166;
t727 = t778 * t1061;
t794 = t454 / 0.2e1 - t949 / 0.2e1;
t208 = t727 - t794;
t835 = t456 / 0.2e1;
t225 = 0.2e1 * t835 + t950;
t292 = t639 * t1134;
t893 = qJD(4) * t1117;
t1101 = qJD(6) * t225 + t208 * t866 + t292 + t893;
t209 = t727 + t794;
t312 = t835 - t456 / 0.2e1;
t1100 = -t312 * qJD(6) + t209 * t866 + t292;
t280 = t468 * t630 + t470 * t631;
t331 = t468 ^ 2 + t470 ^ 2;
t1099 = t331 * qJD(1) + t280 * t866;
t61 = t81 * qJD(1) + t1175;
t1054 = -t678 / 0.2e1;
t1055 = t677 / 0.2e1;
t751 = t1054 * t470 + t1055 * t468;
t834 = t926 / 0.2e1;
t279 = (t834 + t751) * pkin(4);
t1052 = -t680 / 0.2e1;
t750 = t1054 * t631 + t1055 * t630;
t438 = (t1052 + t750) * pkin(4);
t1097 = t279 * qJD(1) + t438 * t866;
t1087 = t680 ^ 2;
t676 = t683 ^ 2;
t497 = (-t1087 / 0.2e1 + t676 / 0.2e1) * t641;
t944 = t680 * t683;
t850 = qJD(1) * t944;
t1096 = t497 * t866 + t633 * t850;
t532 = t630 ^ 2 + t631 ^ 2;
t190 = t280 * qJD(1) + t532 * t866;
t864 = -t676 + t1087;
t461 = -0.2e1 * t641 * t850 + t864 * t866;
t848 = t639 * t877;
t1094 = t771 - t848 + t1180;
t1093 = t848 + t1181;
t1092 = t225 * qJD(1) + t1120;
t52 = qJD(1) * t68 + t1175;
t1040 = t680 * pkin(4);
t573 = pkin(5) * t631 + t1040;
t863 = pkin(4) * t926;
t769 = pkin(5) * t470 - t863;
t1090 = t769 * t1075 + t842 - t573 * t315 / 0.2e1;
t1089 = t251 + t769 * t778 / 0.2e1 + t573 * t1117 / 0.2e1;
t1086 = -pkin(3) / 0.2e1;
t1084 = t166 / 0.2e1;
t1082 = t468 / 0.2e1;
t1079 = t1116 / 0.2e1;
t1071 = t1114 / 0.2e1;
t1064 = t572 / 0.2e1;
t1060 = t641 / 0.2e1;
t628 = t645 / 0.2e1;
t1058 = t668 / 0.2e1;
t1056 = -t670 / 0.2e1;
t1053 = t678 / 0.2e1;
t1048 = pkin(4) * t677;
t1037 = t68 * qJD(4) + t81 * qJD(6);
t1036 = t1178 * t1160;
t1035 = pkin(3) * qJD(3);
t1034 = pkin(4) * qJD(4);
t1031 = t682 * t95;
t963 = t679 * t113;
t6 = (-t963 + t1031) * t639 + t1162;
t1033 = t6 * qJD(1);
t7 = -(t941 + t1032) * t639 + t1161;
t1030 = t7 * qJD(1);
t943 = t682 * t107;
t962 = t679 * t130;
t8 = (t943 - t962) * t639 + t1162;
t1029 = t8 * qJD(1);
t9 = -(t940 + t965) * t639 + t1161;
t1028 = t9 * qJD(1);
t141 = t182 - t1047;
t939 = t682 * t141;
t183 = t678 * t320 - t304;
t142 = t183 - t1045;
t960 = t679 * t142;
t72 = t939 - t960;
t50 = t315 * t769 + t72 * t639 + t1142;
t1027 = qJD(1) * t50;
t938 = t682 * t142;
t961 = t679 * t141;
t73 = t938 + t961;
t51 = t1117 * t769 - t73 * t639 - t1155;
t1026 = qJD(1) * t51;
t54 = t639 * t65 - t1155;
t1025 = qJD(1) * t54;
t55 = -t639 * t66 + t1142;
t1024 = qJD(1) * t55;
t71 = -t166 * t470 + t167 * t468;
t1022 = qJD(1) * t71;
t83 = -t1117 * t314 - t315 * t318;
t1020 = qJD(1) * t83;
t1016 = t168 * t631;
t1015 = t169 * t630;
t1014 = t176 * t631;
t1013 = t177 * t630;
t1012 = t781 * t641;
t1011 = t265 * t641;
t1010 = t779 * t641;
t1009 = t301 * t641;
t783 = -t166 * t469 + t167 * t467;
t31 = -t168 * t470 + t169 * t468 + t783;
t1008 = t31 * qJD(1);
t32 = -t176 * t470 + t177 * t468 + t783;
t1007 = t32 * qJD(1);
t34 = -t1149 * t470 + (-t166 + t183) * t468;
t999 = t34 * qJD(1);
t35 = t166 * t168 + t167 * t169 + t1144;
t998 = t35 * qJD(1);
t44 = t166 * t176 + t167 * t177 + t1144;
t997 = t44 * qJD(1);
t45 = t166 * t182 + t167 * t183 - t440 * t863;
t996 = t45 * qJD(1);
t995 = t1116 * t469;
t994 = t1116 * t631;
t993 = t482 * t467;
t992 = t482 * t630;
t989 = t778 * t639;
t988 = t1114 * t469;
t987 = t1114 * t631;
t986 = t540 * t467;
t985 = t540 * t630;
t981 = t552 * t314;
t980 = t552 * t318;
t979 = t552 * t495;
t978 = t570 * t314;
t977 = t570 * t318;
t976 = t570 * t495;
t947 = t679 * t572;
t929 = t682 * t571;
t93 = t486 * t639 + t1130;
t924 = t93 * qJD(1);
t94 = t1129 + (t917 - t1150) * t639;
t923 = t94 * qJD(1);
t341 = -t571 * t631 + t572 * t630;
t516 = t532 * qJD(5);
t919 = t341 * qJD(3) + t516;
t156 = -t314 * t639 + t315 * t641;
t914 = qJD(1) * t156;
t157 = -t1117 * t641 + t318 * t639;
t913 = qJD(1) * t157;
t210 = t989 / 0.2e1 + t794;
t185 = qJD(1) * t210;
t911 = qJD(1) * t1164;
t233 = t350 * t639 + t793 * t971;
t910 = qJD(1) * t233;
t234 = -t351 * t639 - t793 * t926;
t909 = qJD(1) * t234;
t417 = t865 * t680;
t907 = qJD(1) * t417;
t418 = t865 * t683;
t906 = qJD(1) * t418;
t483 = t1049 * t639 - t641 * t671;
t905 = qJD(1) * t483;
t484 = -t1049 * t641 - t639 * t671;
t904 = qJD(1) * t484;
t903 = qJD(1) * t639;
t901 = qJD(1) * t671;
t900 = qJD(1) * t684;
t899 = qJD(2) * t552;
t898 = qJD(2) * t669;
t897 = qJD(2) * t681;
t896 = qJD(2) * t684;
t895 = qJD(3) * t570;
t894 = qJD(3) * t671;
t892 = qJD(4) * t680;
t891 = qJD(4) * t683;
t890 = qJD(5) * t639;
t888 = qJD(6) * t552;
t887 = qJD(6) * t570;
t102 = t511 * t639 + t1130;
t886 = t102 * qJD(1);
t103 = t1129 + (t916 - t1150) * t639;
t885 = t103 * qJD(1);
t770 = -t669 / 0.2e1 + t803;
t802 = -t861 / 0.2e1;
t688 = (t1058 + t802 - pkin(9) / 0.2e1) * t641 + (t1086 + t770) * t639;
t178 = t680 * t688;
t884 = t178 * qJD(1);
t188 = t208 * qJD(1);
t184 = t209 * qJD(1);
t879 = t312 * qJD(1);
t874 = t865 * qJD(1);
t873 = t497 * qJD(1);
t872 = t498 * qJD(1);
t492 = t503 * qJD(1);
t509 = t864 * t633;
t871 = t509 * qJD(1);
t541 = t628 + t1059;
t870 = t541 * qJD(1);
t869 = t623 * qJD(1);
t660 = -t681 ^ 2 + t684 ^ 2;
t867 = t660 * qJD(1);
t860 = pkin(1) * t681 * qJD(1);
t859 = pkin(1) * t900;
t857 = t1040 / 0.2e1;
t855 = pkin(4) * t678 + pkin(5);
t852 = t639 * t901;
t851 = t641 * t901;
t849 = qJD(4) * t639 * t641;
t663 = t680 * t891;
t844 = t681 * t900;
t843 = t1004 / 0.2e1;
t841 = t440 * t1052;
t600 = t1048 * t679 - t682 * t855;
t839 = t600 * t1060;
t601 = t1048 * t682 + t679 * t855;
t838 = t601 * t1060;
t837 = -t960 / 0.2e1;
t836 = -t938 / 0.2e1;
t833 = t141 / 0.2e1 + t110 / 0.2e1;
t828 = t1156 + t1079;
t826 = -t485 / 0.2e1 + t840;
t825 = t510 / 0.2e1 - t1150 / 0.2e1;
t824 = t1157 + t1071;
t821 = t1065 + t1068;
t818 = t1051 * qJD(2);
t817 = t1051 * qJD(3);
t816 = t1050 * qJD(2);
t815 = t1050 * qJD(3);
t809 = t866 * t683;
t807 = -t863 / 0.2e1;
t806 = -qJD(6) - t903;
t805 = pkin(2) * t815;
t804 = pkin(2) * t816;
t801 = t861 / 0.2e1;
t797 = t89 / 0.2e1 + t600 * t1062;
t791 = t680 * t809;
t789 = t866 * t944;
t197 = t1040 * t649;
t697 = t1077 * t183 + t1149 * t1156 + t1176 * t166;
t762 = t1053 * t168 + t1055 * t169;
t4 = (t649 * t834 + t762 + t841) * pkin(4) + t697;
t785 = -t4 * qJD(1) + t197 * qJD(2);
t218 = t1116 * t571 + t482 * t572 + t649 * t861;
t686 = t571 * t1084 + t167 * t1064 + t176 * t1079 + t177 * t1176 + t1143 / 0.2e1 + t440 * t801;
t712 = t1056 * t1115 + t1069 * t169 + t1157 * t168;
t3 = t686 + t712;
t784 = t3 * qJD(1) + t218 * qJD(2);
t777 = t972 - t973;
t754 = -t571 * t470 / 0.2e1 + t468 * t1064;
t13 = (-t176 / 0.2e1 + t168 / 0.2e1) * t631 + (t177 / 0.2e1 - t169 / 0.2e1) * t630 + (t1156 + t1071) * t469 + (t1176 + t1069) * t467 + t754;
t776 = -qJD(1) * t13 - qJD(2) * t341;
t685 = t839 + t1090;
t705 = t1068 * t1117 - t1193;
t766 = -t963 / 0.2e1 + t1031 / 0.2e1;
t19 = t685 + t705 + t766;
t382 = t573 * t495;
t737 = t552 * t778;
t256 = t382 + t737;
t775 = -t19 * qJD(1) + t256 * qJD(2);
t689 = t838 - t1089;
t20 = t689 + t1196;
t383 = t573 * t778;
t257 = t383 - t979;
t774 = -t20 * qJD(1) + t257 * qJD(2);
t288 = t992 - t994;
t756 = t1077 * t468 + t1079 * t470;
t57 = t756 + t1108;
t773 = qJD(1) * t57 - qJD(2) * t288;
t772 = t641 * (qJD(4) + t903);
t768 = t1042 / 0.2e1 - t1041 / 0.2e1;
t764 = t839 - t1090;
t763 = t838 + t1089;
t761 = t1053 * t176 + t1055 * t177;
t759 = t1067 * t1117 + t1193;
t757 = t1065 * t1117 + t1192;
t755 = t1069 * t468 + t1071 * t470;
t753 = t1053 * t571 + t1055 * t572;
t752 = t1057 * t641 + t1058 * t639;
t748 = -t962 / 0.2e1 + t943 / 0.2e1;
t745 = -t947 / 0.2e1 + t929 / 0.2e1;
t40 = t842 - t759 + t766;
t744 = qJD(1) * t40 - t778 * t899;
t41 = t843 + t1196;
t743 = qJD(1) * t41 + t495 * t899;
t742 = t683 * t772;
t699 = t680 * t752 + t840;
t198 = t699 - t826;
t741 = -qJD(1) * t198 - t683 * t898;
t728 = t752 * t683;
t200 = -t486 / 0.2e1 - t728;
t740 = -qJD(1) * t200 - t680 * t898;
t736 = t570 * t778;
t735 = pkin(3) / 0.2e1 + t770;
t732 = t768 * t683;
t731 = t791 * t1126;
t730 = (t1054 * t469 + t1055 * t467) * pkin(4);
t696 = t1069 * t183 + t1073 * t166 + t1149 * t1157;
t10 = (t670 * t834 + t761 + t841) * pkin(4) + t696;
t245 = t1040 * t670;
t695 = t1165 * t1116;
t84 = ((t1056 - t649 / 0.2e1) * t680 + t753) * pkin(4) + t695;
t726 = -t10 * qJD(1) - t84 * qJD(2) + t245 * qJD(3);
t687 = (t1159 + t182 / 0.2e1) * t631 + (-t183 / 0.2e1 + t1084) * t630 + t730;
t14 = t468 * t828 + t1183 + t687;
t79 = -t827 * t631 + (-t824 - t828) * t630;
t725 = qJD(1) * t14 - qJD(3) * t79;
t17 = t468 * t824 + t687;
t724 = qJD(1) * t17 - qJD(2) * t79;
t691 = t382 + t737 / 0.2e1 + t736 / 0.2e1;
t131 = -t691 + t745;
t703 = t1066 * t1117 - t1192;
t23 = t685 + t703 + t748;
t266 = t382 + t736;
t723 = -t23 * qJD(1) - t131 * qJD(2) + t266 * qJD(3);
t132 = t193 - t383;
t24 = t689 + t1197;
t267 = t383 - t976;
t722 = -t24 * qJD(1) - t132 * qJD(2) + t267 * qJD(3);
t162 = t801 + (t1071 + t1079) * t631 + (t1069 + t1077) * t630;
t326 = t985 - t987;
t59 = t755 + t1108;
t721 = -qJD(1) * t59 - qJD(2) * t162 + qJD(3) * t326;
t694 = (t929 - t947) * t1062 + t315 * t802;
t36 = (-t779 / 0.2e1 + t781 / 0.2e1) * t641 + t821 * t314 + t694;
t720 = t36 * qJD(1) - t495 * t804;
t693 = (t928 + t948) * t1061 + t1117 * t802;
t38 = (t301 / 0.2e1 - t265 / 0.2e1) * t641 + t821 * t318 + t693;
t719 = t38 * qJD(1) - t778 * t804;
t181 = t683 * t688;
t718 = -t181 * qJD(1) - t680 * t804;
t701 = t680 * t768 + t840;
t220 = t701 + t825;
t563 = t735 * t683;
t717 = -qJD(1) * t220 + qJD(2) * t563 + t1035 * t683;
t222 = -t511 / 0.2e1 - t732;
t562 = t735 * t680;
t716 = -qJD(1) * t222 + qJD(2) * t562 + t1035 * t680;
t715 = (-t816 - t815) * pkin(2);
t192 = t745 - t1122;
t46 = t842 + t748 - t757;
t714 = qJD(1) * t46 + qJD(2) * t192 - t778 * t895;
t47 = t843 + t1197;
t713 = qJD(1) * t47 + qJD(2) * t193 + t495 * t895;
t692 = t601 * t1061 + t700 / 0.2e1;
t27 = t682 * t833 + t692 + t837;
t707 = qJD(1) * t27 + qJD(4) * t601;
t173 = t1083 + t1169 / 0.2e1;
t28 = -t679 * t833 + t797 + t836;
t706 = qJD(1) * t28 - qJD(3) * t173 - qJD(4) * t600;
t698 = t730 + (-t166 / 0.2e1 + t183 / 0.2e1) * t630 + t1149 * t1158;
t658 = t680 * t805;
t653 = t864 * qJD(4);
t575 = t601 * qJD(6);
t574 = t600 * qJD(6);
t565 = (t1086 + t1109) * t683;
t564 = pkin(3) * t1052 + t1109 * t680;
t507 = t866 * t623;
t488 = t497 * qJD(4);
t472 = (-t630 * t678 - t631 * t677) * t1034;
t466 = 0.2e1 * t628 - t643;
t464 = t778 * t805;
t463 = t495 * t805;
t458 = t492 + t891;
t457 = -t872 - t892;
t437 = pkin(4) * t750 + t857;
t423 = t438 * qJD(4);
t422 = t438 * qJD(5);
t421 = t437 * qJD(4);
t420 = t437 * qJD(5);
t407 = t440 * t857;
t400 = t789 - t873;
t399 = -t791 + t873;
t379 = t742 * t1125;
t352 = -t676 * t845 - t488;
t323 = qJD(4) * t503 - t906;
t322 = -qJD(4) * t498 + t907;
t278 = pkin(4) * t751 + t807;
t275 = -t488 + (t676 * t902 - t791) * t639;
t268 = t280 * qJD(5);
t253 = t1123 / 0.2e1;
t247 = -t1121 * t680 + t906;
t246 = -t641 * t809 - t907;
t232 = (qJD(4) - t903) * t926 * t1125 + t864 * t527;
t223 = t1151 + t511 / 0.2e1 - t732;
t221 = t701 - t825;
t217 = -t989 / 0.2e1 + t794;
t207 = t1164 * qJD(5);
t206 = t1163 * qJD(5);
t205 = t208 * qJD(5);
t204 = t209 * qJD(5);
t202 = t918 * t778;
t201 = t1151 + t486 / 0.2e1 - t728;
t199 = t699 + t826;
t195 = t745 + t1122;
t194 = t746 + (t1066 + t1068) * t495;
t180 = pkin(9) * t834 + t1113 * t683 + t1119;
t179 = -t1118 + pkin(9) * t971 / 0.2e1 + t1113 * t680;
t174 = t1187 + 0.2e1 * t1083;
t163 = t985 / 0.2e1 + t992 / 0.2e1 - t987 / 0.2e1 - t994 / 0.2e1 + t801;
t134 = t383 - t976 / 0.2e1 - t979 / 0.2e1 + t746;
t133 = t691 + t745;
t119 = -t185 - t1091;
t116 = -t1091 - t184;
t115 = t911 - t918;
t114 = t1091 + t188;
t92 = t1178 * t1132;
t85 = pkin(4) * t753 - t695 + (t649 + t670) * t857;
t78 = t79 * qJD(4);
t70 = -t1163 * t1178 - t913;
t69 = -qJD(4) * t209 - qJD(6) * t210 - t914;
t64 = -t1121 * t778 + t1164 * t1178 + t913;
t63 = -qJD(4) * t208 + qJD(6) * t217 + t1121 * t495 + t914;
t60 = -t755 + t1107;
t58 = -t756 + t1107;
t56 = -t1134 * t318 + t1036;
t49 = t253 + t748 + t757;
t48 = t251 + t1195;
t43 = t253 + t759 + t766;
t42 = t251 + t1194;
t39 = t1011 / 0.2e1 + t980 / 0.2e1 + t1173 + t1009 / 0.2e1 + t977 / 0.2e1 - t693;
t37 = -t1012 / 0.2e1 + t981 / 0.2e1 + t1174 - t1010 / 0.2e1 + t978 / 0.2e1 - t694;
t33 = t1135 * t318 + t1036;
t30 = -t942 / 0.2e1 + t837 + t939 / 0.2e1 - t692;
t29 = t964 / 0.2e1 + t836 - t961 / 0.2e1 - t797;
t26 = t763 + t1195;
t25 = -t703 + t748 + t764;
t22 = t763 + t1194;
t21 = -t705 + t764 + t766;
t18 = t1114 * t1082 + t468 * t1157 + t1165 * t470 + t698;
t16 = -t1020 + t1037;
t15 = t1116 * t1082 + t468 * t1156 + t1183 + t698;
t12 = t993 / 0.2e1 + t1013 / 0.2e1 - t995 / 0.2e1 - t1014 / 0.2e1 + t986 / 0.2e1 + t1015 / 0.2e1 - t988 / 0.2e1 - t1016 / 0.2e1 + t754;
t11 = pkin(4) * t761 + t670 * t807 + t407 - t696;
t5 = pkin(4) * t762 + t649 * t807 + t407 - t697;
t2 = t686 - t712;
t1 = t1020 + t866 * (-t314 * t778 - t318 * t495) + t1037;
t53 = [0, 0, 0, t681 * t896, t660 * qJD(2), 0, 0, 0, -pkin(1) * t897, -pkin(1) * t896, t1139, -t866 * t865, 0, 0, 0, qJD(2) * t483 - t641 * t894, qJD(2) * t484 - t639 * t894, t1139 * t676 - t633 * t663, t509 * qJD(4) - t639 * t731, t418 * t866 + t680 * t849, -t417 * t866 + t683 * t849, -t1139, qJD(2) * t93 + qJD(3) * t102 + qJD(4) * t234, qJD(2) * t94 + qJD(3) * t103 + qJD(4) * t233, qJD(2) * t31 + qJD(3) * t32 + qJD(4) * t34 + qJD(5) * t331, qJD(2) * t35 + qJD(3) * t44 + qJD(4) * t45 + qJD(5) * t71 (t318 * t866 + t771) * t1117, t1131 * t1178 + t83 * t866, t157 * t866 + t639 * t771 (-qJD(6) * t1117 - t893) * t639 + t866 * t156, -t1139, qJD(2) * t6 + qJD(3) * t8 + qJD(4) * t50 + qJD(6) * t55 - t1117 * t890, qJD(2) * t7 + qJD(3) * t9 + qJD(4) * t51 + qJD(6) * t54 + t315 * t890; 0, 0, 0, t844, t867, t896, -t897, 0, -pkin(7) * t896 - t860, pkin(7) * t897 - t859, t845, -t874, -t527, t1121, 0, -qJD(2) * t915 + qJD(3) * t466 + t905, t1128 + t904, t275, t232, t247, t246, -t429, t924 + (t680 * t777 - t1118) * qJD(2) + t179 * qJD(3) + t201 * qJD(4), t923 + (t683 * t777 + t1119) * qJD(2) + t180 * qJD(3) + t199 * qJD(4), t1008 + (t1015 + t993 - t995 - t1016) * qJD(2) + t12 * qJD(3) + t15 * qJD(4) + t268, t998 + (t1116 * t168 + t169 * t482 + t1143) * qJD(2) + t2 * qJD(3) + t5 * qJD(4) + t58 * qJD(5), t33, t1, t64, t63, t365, t1033 + (t1174 + t981 - t1012) * qJD(2) + t37 * qJD(3) + t21 * qJD(4) - t204 + t43 * qJD(6), t1030 + (t1173 + t980 + t1011) * qJD(2) + t39 * qJD(3) + t22 * qJD(4) + t206 + t42 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t845, -t874, -t527, t1121, 0, qJD(2) * t466 - qJD(3) * t915 - t851, t1128 - t852, t275, t232, t247, t246, -t429, t886 + t179 * qJD(2) + (t680 * t796 - t1118) * qJD(3) + t223 * qJD(4), t885 + t180 * qJD(2) + (t683 * t796 + t1119) * qJD(3) + t221 * qJD(4), t1007 + t12 * qJD(2) + (t1013 + t986 - t988 - t1014) * qJD(3) + t18 * qJD(4) + t268, t997 + t2 * qJD(2) + (t1114 * t176 + t1115 * t670 + t177 * t540) * qJD(3) + t11 * qJD(4) + t60 * qJD(5), t33, t1, t64, t63, t365, t1029 + t37 * qJD(2) + (t1174 + t978 - t1010) * qJD(3) + t25 * qJD(4) - t204 + t49 * qJD(6), t1028 + t39 * qJD(2) + (t1173 + t977 + t1009) * qJD(3) + t26 * qJD(4) + t206 + t48 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1096, t1126 * t789 + t871, t680 * t772, t742, t507, qJD(2) * t201 + qJD(3) * t223 - qJD(4) * t351 + t909, qJD(2) * t199 + qJD(3) * t221 + qJD(4) * t350 + t910, t999 + t15 * qJD(2) + t18 * qJD(3) + (-t468 * t678 - t470 * t677) * t1034, t996 + t5 * qJD(2) + t11 * qJD(3) + t278 * qJD(5) + (t182 * t678 + t183 * t677) * t1034, t1105, t1102, t1094, -t1101, t507, qJD(2) * t21 + qJD(3) * t25 + qJD(4) * t72 + qJD(6) * t30 + t1027, qJD(2) * t22 + qJD(3) * t26 - qJD(4) * t73 + qJD(6) * t29 + t1026; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1099, qJD(2) * t58 + qJD(3) * t60 + qJD(4) * t278 + t1022, 0, 0, 0, 0, 0, -t1100, t1093; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1104, t1103, t315 * t806 + t1180 - t310, -qJD(4) * t225 + t1117 * t806 + t217 * t866, t507, qJD(2) * t43 + qJD(3) * t49 + qJD(4) * t30 + qJD(5) * t312 - qJD(6) * t66 + t1024, qJD(2) * t42 + qJD(3) * t48 + qJD(4) * t29 + qJD(6) * t65 + t1025; 0, 0, 0, -t844, -t867, 0, 0, 0, t860, t859, -t845, t874, 0, 0, 0, qJD(3) * t541 - t905, -t904, t352, t379, t323, t322, t429, qJD(3) * t178 + qJD(4) * t200 - t924, qJD(3) * t181 + qJD(4) * t198 - t923, qJD(3) * t13 - qJD(4) * t14 - t1008 + t268, qJD(3) * t3 - qJD(4) * t4 - qJD(5) * t57 - t998, t56, t16, t70, t69, -t365, -qJD(3) * t36 - qJD(4) * t19 - qJD(6) * t40 - t1033 - t205, -qJD(3) * t38 - qJD(4) * t20 - qJD(6) * t41 - t1030 - t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t805, -pkin(2) * t817, t663, -t653, 0, 0, 0, t669 * t892 - t683 * t805, t669 * t891 + t658, t919, qJD(3) * t218 + qJD(4) * t197 + qJD(5) * t288, -t202, t92, 0, 0, 0, qJD(4) * t256 + t778 * t888 + t463, qJD(4) * t257 - t495 * t888 + t464; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t715 + t870 (-t818 - t817) * pkin(2), t663, -t653, 0, 0, 0, t564 * qJD(4) + t683 * t715 + t884, t565 * qJD(4) + t658 - t718, -t776 + t78 + t919 (t1114 * t571 + t572 * t540 + t670 * t861) * qJD(3) + t85 * qJD(4) + t163 * qJD(5) + t784, -t202, t92, 0, 0, 0, t133 * qJD(4) + t195 * qJD(6) + t463 - t720, t134 * qJD(4) + t194 * qJD(6) + t464 - t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t400, -t461, t458, t457, -t869, qJD(3) * t564 - t668 * t891 - t740, qJD(3) * t565 + t668 * t892 - t741, t472 - t725, t85 * qJD(3) + t420 + (t1116 * t677 - t482 * t678) * t1034 + t785, -t1177, t52, t1184, t116, -t869, t133 * qJD(3) - t1191 + t775, t134 * qJD(3) - t1189 + t774; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, qJD(3) * t163 + t421 - t773, 0, 0, 0, 0, 0, -t188, -t911; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1177, t61, t1184, t119, -t869, qJD(3) * t195 - t1191 - t744, qJD(3) * t194 - t1189 - t743; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t845, t874, 0, 0, 0, -qJD(2) * t541 + t851, t852, t352, t379, t323, t322, t429, -qJD(2) * t178 + qJD(4) * t222 - t886, -qJD(2) * t181 + qJD(4) * t220 - t885, -qJD(2) * t13 - qJD(4) * t17 - t1007 + t268, -qJD(2) * t3 - qJD(4) * t10 - qJD(5) * t59 - t997, t56, t16, t70, t69, -t365, qJD(2) * t36 - qJD(4) * t23 - qJD(6) * t46 - t1029 - t205, qJD(2) * t38 - qJD(4) * t24 - qJD(6) * t47 - t1028 - t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t804 - t870, pkin(2) * t818, t663, -t653, 0, 0, 0, -t562 * qJD(4) + t683 * t804 - t884, -t563 * qJD(4) + t718, t516 + t776 + t78, -qJD(4) * t84 - qJD(5) * t162 - t784, -t202, t92, 0, 0, 0, -t131 * qJD(4) - t192 * qJD(6) + t720, -t132 * qJD(4) - t193 * qJD(6) + t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t663, -t653, 0, 0, 0, -pkin(3) * t892, -pkin(3) * t891, t516, qJD(4) * t245 + qJD(5) * t326, -t202, t92, 0, 0, 0, qJD(4) * t266 + t778 * t887, qJD(4) * t267 - t495 * t887; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t400, -t461, t458, t457, -t869, -pkin(9) * t891 - t716, pkin(9) * t892 - t717, t472 - t724, t420 + (t1114 * t677 - t540 * t678) * t1034 + t726, -t1177, t52, t1184, t116, -t869, -t1190 + t723, -qJD(4) * t779 + t174 * qJD(6) + t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, t421 + t721, 0, 0, 0, 0, 0, -t188, -t911; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1177, t61, t1184, t119, -t869, -t1190 - t714, qJD(4) * t174 - qJD(6) * t779 - t713; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1096, -t731 - t871, -t503 * t866 - t680 * t845, t498 * t866 - t683 * t845, t507, -qJD(2) * t200 - qJD(3) * t222 - t909, -qJD(2) * t198 - qJD(3) * t220 - t910, qJD(2) * t14 + qJD(3) * t17 - t999, qJD(2) * t4 + qJD(3) * t10 + qJD(5) * t279 - t996, -t1105, -t1102, t1093, t1100, t507, qJD(2) * t19 + qJD(3) * t23 - qJD(5) * t1117 - qJD(6) * t27 - t1027, qJD(2) * t20 + qJD(3) * t24 - qJD(6) * t28 - t1026 + t1148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t399, t461, -t492, t872, t869, qJD(3) * t562 + t740, qJD(3) * t563 + t741, t725, qJD(3) * t84 + t422 - t785, t1177, -t52, t1179, t184, t869, qJD(3) * t131 - t1110 - t775, qJD(3) * t132 + t1140 - t774; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t399, t461, -t492, t872, t869, t716, t717, t724, t422 - t726, t1177, -t52, t1179, t184, t869, -t1110 - t723, qJD(6) * t173 + t1140 - t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t575, t574; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1097, 0, 0, 0, 0, 0, -t1135, t1167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t879, 0, -t575 - t707, t574 - t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1099, qJD(2) * t57 + qJD(3) * t59 - qJD(4) * t279 - t1022, 0, 0, 0, 0, 0, t1101, t1094; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, qJD(3) * t162 - t423 + t773, 0, 0, 0, 0, 0, t114, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190, -t423 - t721, 0, 0, 0, 0, 0, t114, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1097, 0, 0, 0, 0, 0, t1135, -t1167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1092, -t1167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1104, -t1103, t315 * t903 + t1181, t312 * qJD(4) + t1117 * t903 + t210 * t866, t507, qJD(2) * t40 + qJD(3) * t46 + qJD(4) * t27 - qJD(5) * t225 - t1024, qJD(2) * t41 + qJD(3) * t47 + qJD(4) * t28 - t1025 + t1148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1177, -t61, t1179, t185, t869, qJD(3) * t192 - t1110 + t744, qJD(3) * t193 + t1140 + t743; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1177, -t61, t1179, t185, t869, -t1110 + t714, -qJD(4) * t173 + t1140 + t713; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t879, 0, t707, t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1092, t1167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t53;
