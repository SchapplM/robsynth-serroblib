% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:23
% EndTime: 2019-03-09 05:11:12
% DurationCPUTime: 37.79s
% Computational Cost: add. (44381->859), mult. (86031->1104), div. (0->0), fcn. (105968->10), ass. (0->657)
t814 = qJD(3) + qJD(4);
t1042 = cos(qJ(4));
t1040 = sin(qJ(3));
t1043 = cos(qJ(3));
t644 = sin(pkin(10));
t646 = cos(pkin(10));
t611 = -t1040 * t644 + t1043 * t646;
t613 = t1040 * t646 + t1043 * t644;
t648 = sin(qJ(4));
t558 = -t1042 * t611 + t613 * t648;
t1041 = cos(qJ(6));
t643 = sin(pkin(11));
t802 = t1041 * t643;
t645 = cos(pkin(11));
t647 = sin(qJ(6));
t943 = t647 * t645;
t612 = t802 + t943;
t1125 = t612 * t558;
t1152 = t612 * t1125;
t944 = t647 * t643;
t1126 = t558 * t944;
t801 = t1041 * t645;
t1142 = -t558 * t801 + t1126;
t609 = -t801 + t944;
t1172 = t1142 * t609;
t1186 = -t1172 / 0.2e1;
t1188 = t1186 + t1152 / 0.2e1;
t1193 = t1172 / 0.2e1 - t1152 / 0.2e1 + t1188;
t594 = t1042 * t613;
t941 = t648 * t611;
t1089 = t594 + t941;
t1115 = t1089 * t609;
t1159 = -t1115 / 0.2e1;
t1114 = t1089 * t612;
t1160 = t1114 / 0.2e1;
t301 = 0.2e1 * t1160;
t513 = t1089 * t944;
t745 = t1089 * t801;
t378 = t745 - t513;
t118 = (t378 / 0.2e1 + t1159) * t612 + t609 * t301;
t116 = t118 * qJD(5);
t1192 = qJD(2) * t1193 + t116;
t1191 = 0.2e1 * t1186;
t978 = t1114 * t1142;
t981 = t1125 * t378;
t1176 = -t978 + t981;
t1181 = t1176 * qJD(1);
t136 = t1114 * t609 - t612 * t378;
t135 = t136 * qJD(6);
t1189 = t1181 + t135;
t1187 = t135 - t1181;
t1134 = t1089 ^ 2;
t966 = t558 ^ 2;
t1149 = t1134 - t966;
t1170 = t1149 * t645;
t1185 = qJD(1) * t1170;
t1171 = t1149 * t643;
t1184 = qJD(1) * t1171;
t972 = t378 * t1089;
t975 = t1142 * t558;
t1177 = t972 + t975;
t1183 = qJD(1) * t1177;
t1169 = t1159 + t1115 / 0.2e1;
t1182 = t1169 * qJD(2);
t1113 = t1089 * t643;
t539 = t1113 / 0.2e1;
t779 = -t1113 / 0.2e1;
t1178 = t779 + t539;
t1180 = t1178 * qJD(2);
t1179 = t1152 - t1172;
t1033 = pkin(7) + qJ(2);
t757 = t1033 * t646;
t758 = t1033 * t644;
t574 = -t1040 * t758 + t1043 * t757;
t482 = t611 * pkin(8) + t574;
t465 = t1042 * t482;
t572 = t1040 * t757 + t1043 * t758;
t481 = -t613 * pkin(8) - t572;
t661 = t648 * t481;
t1090 = t465 + t661;
t1123 = t643 * t558;
t1141 = -pkin(5) * t1123 + t1090;
t348 = t1042 * t481 - t648 * t482;
t256 = pkin(5) * t1113 - t348;
t1175 = t1141 * t256;
t1174 = t1141 * t609;
t1173 = t1141 * t612;
t977 = t1114 * t1089;
t980 = t1125 * t558;
t1155 = -t977 + t980;
t1168 = qJD(1) * t1155;
t1167 = t1149 * qJD(1);
t738 = t594 / 0.2e1;
t1088 = t738 + t941 / 0.2e1;
t1099 = t1089 * qJD(1);
t1153 = t558 * t1099;
t1166 = qJD(6) * t1088 + t1153;
t1083 = t661 / 0.2e1 + t465 / 0.2e1;
t1158 = -t1123 / 0.2e1;
t1165 = pkin(5) * t1158 + t1083;
t1111 = t1090 * t645;
t1022 = qJ(5) * t1089;
t1037 = t558 * pkin(4);
t635 = -t646 * pkin(2) - pkin(1);
t579 = -t611 * pkin(3) + t635;
t675 = t579 + t1037;
t662 = t675 - t1022;
t183 = t643 * t662 + t1111;
t158 = -pkin(9) * t1113 + t183;
t1032 = -pkin(9) - qJ(5);
t1112 = t1090 * t643;
t653 = t558 * pkin(5) - t1112 + (t1032 * t1089 + t675) * t645;
t90 = -t1041 * t653 + t647 * t158;
t91 = t1041 * t158 + t647 * t653;
t1164 = t1125 * t91 + t1142 * t90;
t1053 = -t609 / 0.2e1;
t1051 = t612 / 0.2e1;
t1085 = t1051 * t1114;
t189 = t1053 * t378 - t1085;
t1109 = t189 * qJD(6);
t903 = qJD(1) * t378;
t1163 = -t1142 * t903 + t1109;
t1162 = -t1089 * t90 + t1114 * t1141 - t1125 * t256;
t1161 = -t1089 * t91 + t1141 * t378 + t1142 * t256;
t1072 = -t1111 / 0.2e1;
t641 = t645 ^ 2;
t1124 = t641 * t558;
t532 = -t1124 / 0.2e1;
t533 = t1124 / 0.2e1;
t1071 = -t1125 / 0.2e1;
t1157 = pkin(9) * t1123;
t983 = t348 * t643;
t336 = t645 * t348;
t987 = t1090 * t348;
t1034 = t648 * pkin(3);
t633 = qJ(5) + t1034;
t639 = t643 ^ 2;
t1046 = t639 / 0.2e1;
t760 = t641 / 0.2e1 + t1046;
t1154 = (qJ(5) + t633) * t760;
t1029 = t1051 * t90 + t1053 * t91;
t1021 = qJ(5) * t645;
t638 = t645 * pkin(9);
t624 = t638 + t1021;
t759 = t643 * t1032;
t573 = t1041 * t624 + t647 * t759;
t1056 = t573 / 0.2e1;
t571 = -t1041 * t759 + t624 * t647;
t1059 = -t571 / 0.2e1;
t1151 = t1056 * t1114 - t1059 * t1115 - t1029;
t605 = t645 * t633 + t638;
t756 = (-pkin(9) - t633) * t643;
t546 = t1041 * t605 + t647 * t756;
t1066 = -t546 / 0.2e1;
t545 = -t1041 * t756 + t605 * t647;
t1067 = t545 / 0.2e1;
t1150 = t1066 * t1114 - t1067 * t1115 + t1029;
t736 = -t801 / 0.2e1;
t1128 = t558 * t736;
t687 = -t1126 / 0.2e1 - t1128;
t960 = t609 * t558;
t778 = t960 / 0.2e1;
t249 = t778 - t687;
t1108 = t249 * qJD(1);
t599 = t609 * qJD(6);
t1148 = t599 + t1108;
t1007 = t183 * t645;
t182 = t645 * t662 - t1112;
t1008 = t182 * t643;
t706 = t1008 / 0.2e1 - t1007 / 0.2e1;
t1106 = t1083 + t706;
t1147 = qJD(1) * t1106;
t1146 = qJD(5) * t1106;
t1105 = t1083 - t706;
t1145 = t1105 * qJD(5);
t1144 = t1114 * qJD(1);
t1093 = t814 * t612;
t1143 = -t189 * qJD(1) + t609 * t1093;
t1038 = t1089 * pkin(4);
t406 = qJ(5) * t558 + t1038;
t1122 = t645 * t558;
t1140 = pkin(5) * t1089 + pkin(9) * t1122;
t723 = t1089 * t1090 + t348 * t558;
t1138 = t1125 * t1144 - t1109;
t1039 = t1090 * pkin(4);
t1036 = t613 * pkin(3);
t360 = t1036 + t406;
t195 = t643 * t360 + t336;
t1003 = t195 * t645;
t194 = t645 * t360 - t983;
t1006 = t194 * t643;
t1082 = t1003 / 0.2e1 - t1006 / 0.2e1;
t1137 = t1082 * qJ(5) - t1039 / 0.2e1;
t1136 = t814 * t348;
t793 = t643 * t1099;
t743 = t645 * t793;
t728 = t558 * t743;
t1135 = -0.2e1 * t728;
t1061 = t558 / 0.2e1;
t1133 = -t1089 / 0.2e1;
t726 = t182 * t645 + t183 * t643;
t1121 = t726 * t558;
t1119 = t814 * t558;
t1118 = qJD(1) * t558;
t1117 = qJD(2) * t558;
t1116 = qJD(6) * t249;
t768 = t944 / 0.2e1;
t212 = t513 / 0.2e1 + (t768 - t801) * t1089;
t1104 = -t212 * qJD(1) + t1093;
t444 = -t960 / 0.2e1;
t1091 = t444 - t687;
t1103 = qJD(3) * t1091;
t1101 = qJD(6) * t1091;
t1100 = t1088 * qJD(1);
t1097 = -t1114 * t903 + t189 * t814;
t798 = t1114 * t1118;
t1096 = qJD(3) * t249 + t798;
t1095 = t1114 ^ 2;
t1035 = t645 * pkin(5);
t634 = -pkin(4) - t1035;
t1048 = t634 / 0.2e1;
t636 = -pkin(3) * t1042 - pkin(4);
t623 = t636 - t1035;
t1050 = t623 / 0.2e1;
t761 = t1048 + t1050;
t1094 = t612 * t761;
t627 = t639 + t641;
t1092 = t814 * t627;
t677 = t943 / 0.2e1 + t802 / 0.2e1;
t741 = 0.2e1 * t1159;
t1087 = qJD(4) * t741;
t1086 = qJD(6) * t1114;
t131 = -t1115 * t378 + t1095;
t1081 = qJD(1) * t131 + t118 * t814;
t151 = -t378 ^ 2 + t1095;
t1080 = qJD(1) * t151 + t136 * t814;
t1074 = t612 ^ 2;
t1075 = t609 ^ 2;
t488 = -t1074 + t1075;
t128 = t136 * qJD(1) + t488 * t814;
t563 = t1074 + t1075;
t109 = t118 * qJD(1) + t563 * t814;
t1077 = t609 * t814 + t1144;
t1073 = t613 ^ 2;
t1069 = t1126 / 0.2e1;
t1068 = -t545 / 0.2e1;
t1065 = t546 / 0.2e1;
t1063 = t1089 / 0.2e1;
t1062 = -t558 / 0.2e1;
t1058 = t571 / 0.2e1;
t1057 = -t573 / 0.2e1;
t739 = t1042 * t1041;
t807 = t647 * t1042;
t746 = t645 * t807;
t582 = (-t643 * t739 - t746) * pkin(3);
t1055 = -t582 / 0.2e1;
t747 = t643 * t807;
t583 = (t645 * t739 - t747) * pkin(3);
t1054 = t583 / 0.2e1;
t1052 = t609 / 0.2e1;
t1049 = -t634 / 0.2e1;
t1047 = t636 / 0.2e1;
t1044 = t648 / 0.2e1;
t1028 = pkin(3) * qJD(4);
t1027 = qJD(3) * pkin(3);
t142 = t1140 + t194;
t806 = t1041 * t142;
t163 = t195 + t1157;
t946 = t647 * t163;
t94 = t806 - t946;
t1024 = t94 * t612;
t804 = t1041 * t163;
t948 = t647 * t142;
t95 = t804 + t948;
t1023 = t95 * t609;
t664 = t1050 * t1089 + t1065 * t1142 - t1067 * t1125;
t708 = t94 * t1052 - t95 * t612 / 0.2e1;
t32 = t664 + t708;
t1020 = qJD(1) * t32;
t36 = -t1114 * t91 - t1115 * t90;
t1019 = qJD(1) * t36;
t49 = -t1114 * t256 + t558 * t90;
t1018 = qJD(1) * t49;
t50 = t256 * t378 - t558 * t91;
t1017 = qJD(1) * t50;
t984 = t348 * t1089;
t61 = -t984 - (t1007 - t1008) * t558;
t1016 = qJD(1) * t61;
t10 = -t1114 * t95 - t378 * t94 + t1164;
t1014 = t10 * qJD(1);
t206 = t645 * t406 - t983;
t150 = t1140 + t206;
t805 = t1041 * t150;
t207 = t643 * t406 + t336;
t167 = t207 + t1157;
t945 = t647 * t167;
t106 = t805 - t945;
t1013 = t106 * t612;
t803 = t1041 * t167;
t947 = t647 * t150;
t107 = t803 + t947;
t1012 = t107 * t609;
t11 = -t106 * t378 - t107 * t1114 + t1164;
t1011 = t11 * qJD(1);
t12 = -t90 * t94 + t91 * t95 + t1175;
t1010 = t12 * qJD(1);
t13 = -t106 * t90 + t107 * t91 + t1175;
t1009 = t13 * qJD(1);
t1005 = t194 * t645;
t1004 = t195 * t643;
t1002 = t206 * t643;
t1001 = t206 * t645;
t1000 = t207 * t643;
t999 = t207 * t645;
t21 = t558 * t94 + t1162;
t998 = t21 * qJD(1);
t22 = -t558 * t95 + t1161;
t997 = t22 * qJD(1);
t23 = t106 * t558 + t1162;
t996 = t23 * qJD(1);
t24 = -t107 * t558 + t1161;
t995 = t24 * qJD(1);
t992 = t256 * t609;
t991 = t256 * t612;
t31 = t1089 * t256 - t1125 * t90 + t1142 * t91;
t988 = t31 * qJD(1);
t37 = t182 * t194 + t183 * t195 - t987;
t982 = t37 * qJD(1);
t979 = t1125 * t609;
t973 = t1142 * t612;
t38 = t182 * t206 + t183 * t207 - t987;
t971 = t38 * qJD(1);
t39 = (t1004 + t1005) * t1089 - t1121;
t970 = t39 * qJD(1);
t40 = (t1000 + t1001) * t1089 - t1121;
t969 = t40 * qJD(1);
t968 = t546 * t609;
t57 = t1089 * t182 + t194 * t558 + t643 * t723;
t965 = t57 * qJD(1);
t964 = t573 * t609;
t58 = -t1089 * t183 - t195 * t558 + t645 * t723;
t963 = t58 * qJD(1);
t59 = (t182 + t1112) * t1089 + (t206 + t983) * t558;
t962 = t59 * qJD(1);
t60 = (-t183 + t1111) * t1089 + (-t207 + t336) * t558;
t961 = t60 * qJD(1);
t957 = t612 * t545;
t955 = t612 * t571;
t733 = t760 * t558;
t665 = t1047 * t1089 - t633 * t733;
t705 = t1005 / 0.2e1 + t1004 / 0.2e1;
t82 = t665 - t705;
t940 = t82 * qJD(1);
t411 = -t582 * t612 - t583 * t609;
t562 = t563 * qJD(5);
t935 = t411 * qJD(4) + t562;
t762 = 0.2e1 * t1133;
t394 = t762 * t645;
t890 = qJD(3) * t645;
t930 = -t394 * qJD(4) + t1089 * t890;
t676 = t627 * t1042;
t604 = t676 * pkin(3);
t625 = t627 * qJD(5);
t929 = t604 * qJD(4) + t625;
t628 = t644 ^ 2 + t646 ^ 2;
t108 = t726 * t1089;
t928 = qJD(1) * t108;
t120 = 0.2e1 * t1071 * t612 + t1191;
t927 = qJD(1) * t120;
t121 = t1191 - t1152;
t926 = qJD(1) * t121;
t925 = qJD(1) * t1193;
t133 = -t978 - t981;
t923 = qJD(1) * t133;
t143 = -t1090 * t558 - t984;
t922 = qJD(1) * t143;
t154 = t977 + t980;
t918 = qJD(1) * t154;
t157 = t972 - t975;
t915 = qJD(1) * t157;
t185 = (t1125 / 0.2e1 + t1071) * t609;
t914 = qJD(1) * t185;
t236 = t966 + t1134;
t230 = t236 * t643;
t911 = qJD(1) * t230;
t232 = t236 * t645;
t909 = qJD(1) * t232;
t284 = -t1036 * t558 - t1089 * t579;
t906 = qJD(1) * t284;
t285 = -t1036 * t1089 + t558 * t579;
t905 = qJD(1) * t285;
t900 = qJD(1) * t579;
t897 = qJD(2) * t1089;
t763 = t1062 + t1061;
t396 = t763 * t645;
t369 = t643 * t396;
t896 = qJD(3) * t369;
t895 = qJD(3) * t1089;
t893 = qJD(3) * t609;
t892 = qJD(3) * t612;
t891 = qJD(3) * t623;
t889 = qJD(4) * t369;
t887 = qJD(4) * t1089;
t886 = qJD(4) * t579;
t885 = qJD(4) * t609;
t884 = qJD(4) * t612;
t883 = qJD(4) * t634;
t882 = qJD(4) * t643;
t881 = qJD(4) * t645;
t880 = qJD(5) * t558;
t879 = qJD(6) * t558;
t667 = -qJ(5) * t733 - t1038 / 0.2e1;
t703 = t1001 / 0.2e1 + t1000 / 0.2e1;
t104 = t667 - t703;
t878 = t104 * qJD(1);
t812 = -t1042 / 0.2e1;
t655 = t633 * t1133 + t636 * t1062 + (t1044 * t1089 + t558 * t812) * pkin(3);
t650 = qJ(5) * t1063 - t1037 / 0.2e1 + t655;
t111 = t1072 + t1111 / 0.2e1 + t650 * t643;
t877 = t111 * qJD(1);
t115 = t1036 * t579;
t876 = t115 * qJD(1);
t699 = t1053 * t1115 - t1085;
t161 = -t699 + t1088;
t871 = t161 * qJD(1);
t543 = t639 * t558;
t523 = -t543 / 0.2e1;
t775 = t558 * t1046;
t732 = t533 + t775;
t216 = t532 + t523 + t732;
t864 = t216 * qJD(1);
t524 = t543 / 0.2e1;
t219 = t533 + t524 + t732;
t863 = t219 * qJD(1);
t228 = t627 * t1134;
t862 = t228 * qJD(1);
t861 = t236 * qJD(1);
t243 = (t1051 - t677) * t558;
t859 = t243 * qJD(1);
t244 = (t1051 + t677) * t558;
t858 = t244 * qJD(1);
t776 = t558 * t1051;
t245 = -t558 * t677 + t776;
t857 = t245 * qJD(1);
t770 = t1122 / 0.2e1;
t773 = t1123 / 0.2e1;
t246 = t1041 * t773 + t647 * t770 + t776;
t239 = t246 * qJD(1);
t247 = t1069 + (t736 + t1052) * t558;
t856 = t247 * qJD(1);
t248 = t1069 + (t736 + t1053) * t558;
t855 = t248 * qJD(1);
t853 = t1091 * qJD(1);
t591 = -t594 / 0.2e1;
t734 = t760 * t1089;
t271 = -t941 / 0.2e1 + t591 - t734;
t851 = t271 * qJD(1);
t764 = t1063 + t1133;
t291 = t764 * t609;
t849 = t291 * qJD(1);
t848 = t1115 * qJD(1);
t847 = t741 * qJD(1);
t298 = t764 * t612;
t845 = t298 * qJD(1);
t299 = t762 * t612;
t844 = t299 * qJD(1);
t682 = -t1044 * t558 + t1089 * t812;
t340 = (-t613 / 0.2e1 + t682) * pkin(3);
t843 = t340 * qJD(1);
t367 = t572 * t613 + t574 * t611;
t842 = t367 * qJD(1);
t368 = t1122 * t643;
t841 = t368 * qJD(3);
t840 = t368 * qJD(4);
t839 = t369 * qJD(1);
t838 = t1113 * qJD(1);
t388 = t763 * t643;
t837 = t388 * qJD(1);
t390 = t764 * t643;
t836 = t390 * qJD(1);
t392 = t762 * t643;
t835 = t392 * qJD(1);
t393 = t764 * t645;
t834 = t393 * qJD(1);
t833 = t394 * qJD(1);
t832 = t396 * qJD(1);
t397 = t543 + t1124;
t831 = t397 * qJD(1);
t442 = 0.2e1 * t738 + t941;
t829 = t442 * qJD(1);
t606 = t611 ^ 2;
t489 = t606 - t1073;
t828 = t489 * qJD(1);
t555 = t738 + t591;
t825 = t555 * qJD(1);
t824 = t555 * qJD(4);
t564 = t606 + t1073;
t821 = t564 * qJD(1);
t820 = t611 * qJD(1);
t600 = t611 * qJD(3);
t819 = t612 * qJD(6);
t818 = t613 * qJD(1);
t817 = t613 * qJD(3);
t621 = t628 * qJ(2);
t816 = t621 * qJD(1);
t815 = t628 * qJD(1);
t811 = t648 * t1027;
t810 = t648 * t1028;
t809 = t1034 / 0.2e1;
t408 = t1051 * t583 + t1053 * t582;
t808 = t408 * qJD(4);
t800 = t1115 * t1118;
t797 = t558 * t900;
t796 = t1089 * t900;
t795 = t639 * t1099;
t794 = t641 * t1099;
t792 = t643 * t890;
t791 = t643 * t881;
t790 = t1089 * t880;
t789 = t1089 * t1118;
t788 = t609 * t819;
t787 = t611 * t818;
t786 = t611 * t817;
t785 = t645 * t1099;
t784 = t992 / 0.2e1;
t783 = -t991 / 0.2e1;
t755 = t1042 * qJD(3);
t754 = t1042 * qJD(4);
t749 = -qJD(1) * t635 - qJD(2);
t748 = -qJD(6) - t1118;
t744 = t645 * t1153;
t740 = t532 + t775;
t735 = t814 * t1034;
t729 = -t1022 + t1037;
t660 = -t1054 * t1114 + t1055 * t378 - t1066 * t1125 + t1067 * t1142;
t695 = -t1056 * t1125 + t1059 * t1142;
t4 = (-t106 / 0.2e1 + t94 / 0.2e1) * t612 + (-t107 / 0.2e1 + t95 / 0.2e1) * t609 + t660 + t695;
t727 = -qJD(1) * t4 - qJD(3) * t411;
t725 = t1003 - t1006;
t724 = t999 - t1002;
t721 = -t1089 * t633 - t636 * t558;
t702 = t999 / 0.2e1 - t1002 / 0.2e1;
t649 = t702 * t633 + (-t1042 * t706 - t1044 * t348) * pkin(3) + t1090 * t1047;
t26 = t649 - t1137;
t548 = (t633 * t676 + t636 * t648) * pkin(3);
t720 = t26 * qJD(1) + t548 * qJD(3);
t27 = t1165 - t1150;
t287 = t957 - t968;
t719 = qJD(1) * t27 - qJD(3) * t287;
t663 = t1048 * t1089 + t1056 * t1142 - t1058 * t1125;
t707 = t1051 * t107 + t1053 * t106;
t34 = t663 - t707;
t718 = -qJD(1) * t34 + qJD(3) * t408;
t56 = (t207 / 0.2e1 - t195 / 0.2e1) * t645 + (-t206 / 0.2e1 + t194 / 0.2e1) * t643;
t717 = -qJD(1) * t56 - qJD(3) * t604;
t584 = t627 * t633;
t716 = -qJD(3) * t584 + t1147;
t715 = -t739 / 0.2e1;
t714 = -qJD(3) * t245 - qJD(4) * t243;
t713 = qJD(3) * t246 + qJD(4) * t244;
t712 = qJD(4) * t442 + t895;
t697 = -t1050 * t1114 + t1061 * t545;
t696 = t1050 * t378 + t1062 * t546;
t694 = -t1049 * t1125 + t1058 * t1089;
t693 = -t1048 * t1114 + t1058 * t558;
t692 = t1049 * t1142 + t1056 * t1089;
t691 = t1048 * t378 + t1057 * t558;
t679 = -t946 / 0.2e1 + t806 / 0.2e1;
t41 = t783 + t679 - t696;
t690 = qJD(1) * t41 - t612 * t891;
t681 = -t948 / 0.2e1 - t804 / 0.2e1;
t42 = t784 + t681 - t697;
t689 = qJD(1) * t42 + t609 * t891;
t688 = t1069 + t1128;
t234 = (t887 + t895) * t558;
t686 = t1119 * t1089;
t657 = t1034 * t1160 - t1050 * t1125 + t1061 * t582 + t1133 * t545;
t52 = t657 + t694;
t685 = -t52 * qJD(1) - t609 * t811;
t656 = t1050 * t1142 + t1062 * t583 + t1133 * t546 + t378 * t809;
t54 = t656 + t692;
t684 = -t54 * qJD(1) - t612 * t811;
t114 = t650 * t645;
t683 = -t114 * qJD(1) - t643 * t811;
t680 = -t947 / 0.2e1 - t803 / 0.2e1;
t678 = -t945 / 0.2e1 + t805 / 0.2e1;
t652 = t1050 * t1141 + t1054 * t91 + t1055 * t90 + t106 * t1068 + t1065 * t107 + t256 * t809;
t666 = t1049 * t1141 + t1057 * t95 + t1058 * t94;
t2 = t652 + t666;
t242 = t1034 * t623 - t545 * t582 + t546 * t583;
t673 = t2 * qJD(1) + t408 * qJD(2) + t242 * qJD(3);
t176 = t809 + (t1059 + t1068) * t612 + (t1056 + t1065) * t609;
t29 = t1165 + t1151;
t366 = t955 - t964;
t671 = -qJD(1) * t29 - qJD(3) * t176 + qJD(4) * t366;
t541 = t809 - t1154;
t622 = t627 * qJ(5);
t670 = qJD(3) * t541 - qJD(4) * t622 + t1147;
t659 = (-t746 / 0.2e1 + t643 * t715) * pkin(3);
t412 = t659 - t1094;
t45 = t783 + t678 - t691;
t669 = qJD(1) * t45 + qJD(3) * t412 - t612 * t883;
t658 = (t645 * t715 + t747 / 0.2e1) * pkin(3);
t413 = t609 * t761 + t658;
t46 = t784 + t680 - t693;
t668 = qJD(1) * t46 + qJD(3) * t413 + t609 * t883;
t626 = t643 * t810;
t603 = t612 * qJD(5);
t598 = t609 * qJD(5);
t586 = t612 * t810;
t585 = t609 * t810;
t542 = t809 + t1154;
t472 = t488 * qJD(6);
t415 = t659 + t1094;
t414 = t1049 * t609 + t1053 * t623 + t658;
t400 = t814 * t1088;
t389 = 0.2e1 * t539;
t387 = 0.2e1 * t1158;
t381 = t393 * qJD(4);
t341 = t1036 / 0.2e1 + t682 * pkin(3);
t302 = t1160 - t1114 / 0.2e1;
t270 = -t734 + t1088;
t254 = t778 + t688;
t252 = t444 + t688;
t241 = t245 * qJD(6);
t240 = t246 * qJD(6);
t225 = t991 / 0.2e1;
t224 = -t992 / 0.2e1;
t221 = -t819 - t239;
t218 = t533 + t523 + t740;
t217 = t532 + t524 + t740;
t215 = t745 / 0.2e1 - t513 / 0.2e1 + (t768 + t736) * t1089;
t199 = -0.2e1 * t1083;
t186 = t1142 * t1051 + t973 / 0.2e1;
t184 = -t1125 * t1052 - t979 / 0.2e1;
t177 = t809 + t955 / 0.2e1 - t964 / 0.2e1 + t957 / 0.2e1 - t968 / 0.2e1;
t162 = t699 + t1088;
t122 = 0.2e1 * t1188;
t113 = pkin(4) * t770 + t1021 * t1133 + t655 * t645 + t1112;
t112 = pkin(4) * t773 + qJ(5) * t779 + t655 * t643 + 0.2e1 * t1072;
t105 = t667 + t703;
t83 = t665 + t705;
t55 = t702 + t1082;
t53 = t1173 + t656 - t692;
t51 = t1174 + t657 - t694;
t48 = t225 + t678 + t691;
t47 = t224 + t680 + t693;
t44 = t225 + t679 + t696;
t43 = t224 + t681 + t697;
t35 = t663 + t707;
t33 = t664 - t708;
t30 = t1165 - t1151;
t28 = t1165 + t1150;
t25 = t649 + t1137;
t3 = -t1012 / 0.2e1 - t1013 / 0.2e1 - t1023 / 0.2e1 - t1024 / 0.2e1 + t660 - t695;
t1 = t652 - t666;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t628 * qJD(2), t621 * qJD(2), t786, t489 * qJD(3), 0, -t786, 0, 0, t635 * t817, t635 * t600, qJD(2) * t564, qJD(2) * t367, -t686, -t814 * t1149, 0, t234, 0, 0, -qJD(3) * t284 + t1089 * t886, -qJD(3) * t285 - t558 * t886, qJD(2) * t236, qJD(2) * t143 + qJD(3) * t115, -t641 * t686, 0.2e1 * t1119 * t645 * t1113, t814 * t1170, -t639 * t686, -t814 * t1171, t234, qJD(2) * t230 + qJD(3) * t57 + qJD(4) * t59 - t645 * t790, qJD(2) * t232 + qJD(3) * t58 + qJD(4) * t60 + t643 * t790, -qJD(3) * t39 - qJD(4) * t40 + qJD(5) * t228, qJD(2) * t61 + qJD(3) * t37 + qJD(4) * t38 - qJD(5) * t108 (t1142 * t814 - t1086) * t378, qJD(6) * t151 + t1176 * t814, -t1114 * t879 + t1177 * t814 (qJD(6) * t378 - t1125 * t814) * t1114, t1155 * t814 - t378 * t879, t234, qJD(2) * t154 + qJD(3) * t21 + qJD(4) * t23 + qJD(6) * t50 + t1115 * t880, qJD(2) * t157 + qJD(3) * t22 + qJD(4) * t24 + qJD(6) * t49 + t1114 * t880, qJD(2) * t133 + qJD(3) * t10 + qJD(4) * t11 + qJD(5) * t131, qJD(2) * t31 + qJD(3) * t12 + qJD(4) * t13 + qJD(5) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t815, t816, 0, 0, 0, 0, 0, 0, 0, 0, t821, t842, 0, 0, 0, 0, 0, 0, t824, 0, t861, qJD(3) * t341 + t922, 0, 0, 0, 0, 0, 0, -t381 + t911, t1178 * t814 + t909, t218 * qJD(4), qJD(3) * t83 + qJD(4) * t105 + qJD(5) * t270 + t1016, 0, 0, 0, 0, 0, 0, t1169 * t814 - t241 + t918, qJD(4) * t302 - t1101 + t915, t1193 * t814 + t923, t988 + (t973 - t979) * qJD(2) + t33 * qJD(3) + t35 * qJD(4) + t162 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t787, t828, t600, -t787, -t817, 0, -qJD(3) * t574 + t635 * t818, qJD(3) * t572 + t635 * t820, 0, 0, -t789, -t1167, -t1119, t1153, -t712, 0, -qJD(3) * t1090 + qJD(4) * t199 - t906, -t1136 - t905 (t1042 * t558 - t1089 * t648) * t1027, t876 + t341 * qJD(2) + (-t1042 * t1090 + t348 * t648) * t1027, -t840 - (t792 + t794) * t558, 0.2e1 * t728 + (t543 - t1124) * qJD(3) + t217 * qJD(4), qJD(4) * t389 + t643 * t895 + t1185, t840 - (-t792 + t795) * t558, -t1184 + t930, t1153, t965 + (t643 * t721 - t1111) * qJD(3) + t112 * qJD(4) + t387 * qJD(5), t963 + t1180 + (t645 * t721 + t1112) * qJD(3) + t113 * qJD(4) - t1122 * qJD(5), qJD(3) * t725 + t55 * qJD(4) - t970, t982 + t83 * qJD(2) + (t1090 * t636 + t633 * t725) * qJD(3) + t25 * qJD(4) + t1145, t186 * qJD(4) + t1109 + (t892 + t903) * t1142, qJD(3) * t1179 + t122 * qJD(4) + t1189, qJD(4) * t301 + t1089 * t892 + t1101 + t1183, t184 * qJD(4) - t1109 - (t893 + t1144) * t1125, -t1089 * t893 + t1087 + t1168 - t241, t1166, t998 + t1182 + (-t1089 * t545 - t1125 * t623 + t1174) * qJD(3) + t51 * qJD(4) - t246 * qJD(5) + t44 * qJD(6), t997 + (-t1089 * t546 + t1142 * t623 + t1173) * qJD(3) + t53 * qJD(4) + t249 * qJD(5) + t43 * qJD(6), t1014 + (t1125 * t546 + t1142 * t545 - t1023 - t1024) * qJD(3) + t3 * qJD(4) + t1192, t1010 + t33 * qJD(2) + (t1141 * t623 - t545 * t94 + t546 * t95) * qJD(3) + t1 * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1153, -t1167, -t1119, t1153, -qJD(3) * t442 - t887, 0, qJD(2) * t555 + qJD(3) * t199 - qJD(4) * t1090 + t796, -t1136 - t797, 0, 0, -t841 + (-t791 - t794) * t558, t217 * qJD(3) + (0.2e1 * t743 + (t639 - t641) * qJD(4)) * t558, qJD(3) * t389 + t1089 * t882 + t1185, t841 + (t791 - t795) * t558, -qJD(3) * t394 + t1089 * t881 - t1184, t1153, t962 - t393 * qJD(2) + t112 * qJD(3) + (t643 * t729 - t1111) * qJD(4) - t1123 * qJD(5), t1090 * t882 + t961 + t1180 + t113 * qJD(3) + (qJD(4) * t729 - t880) * t645, t218 * qJD(2) + t55 * qJD(3) + qJD(4) * t724 - t969, t971 + t105 * qJD(2) + t25 * qJD(3) + (qJ(5) * t724 - t1039) * qJD(4) + t1145, t186 * qJD(3) + t1109 + (t884 + t903) * t1142, t122 * qJD(3) + qJD(4) * t1179 + t1189, qJD(3) * t301 + qJD(6) * t252 + t1089 * t884 + t1183, t184 * qJD(3) - t1109 - (t885 + t1144) * t1125, qJD(3) * t741 - qJD(6) * t243 - t1089 * t885 + t1168, t1166, t996 + t1182 + t51 * qJD(3) + (-t1089 * t571 - t1125 * t634 + t1174) * qJD(4) - t244 * qJD(5) + t48 * qJD(6), t995 + t302 * qJD(2) + t53 * qJD(3) + (-t1089 * t573 + t1142 * t634 + t1173) * qJD(4) + t254 * qJD(5) + t47 * qJD(6), t1011 + t3 * qJD(3) + (t1125 * t573 + t1142 * t571 - t1012 - t1013) * qJD(4) + t1192, t1009 + t35 * qJD(2) + t1 * qJD(3) + (-t106 * t571 + t107 * t573 + t1141 * t634) * qJD(4) + t30 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t387 - qJD(4) * t1123 - t744, -t1122 * qJD(3) + (t793 - t881) * t558, t862, qJD(2) * t270 + t1105 * t814 - t928, 0, 0, 0, 0, 0, 0, qJD(6) * t215 - t713 + t800, qJD(4) * t254 + t1096, t1081, qJD(2) * t162 + qJD(3) * t28 + qJD(4) * t30 + t1019; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1097, t1080, qJD(4) * t252 + t1114 * t748 + t1103, -t1097, t378 * t748 + t714, t400, -qJD(2) * t245 + qJD(3) * t44 + qJD(4) * t48 + qJD(5) * t215 - qJD(6) * t91 + t1017, -qJD(2) * t1091 + qJD(3) * t43 + qJD(4) * t47 + qJD(6) * t90 + t1018, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t815, -t816, 0, 0, 0, 0, 0, 0, t817, t600, -t821, -t842, 0, 0, 0, 0, 0, 0, t712, -t1119, -t861, -qJD(3) * t340 - t922, 0, 0, 0, 0, 0, 0, -t911 + t930, -qJD(3) * t1113 + qJD(4) * t392 - t909, qJD(3) * t397 + qJD(4) * t219, -qJD(3) * t82 - qJD(4) * t104 + qJD(5) * t271 - t1016, 0, 0, 0, 0, 0, 0, -qJD(3) * t1115 + t1087 - t240 - t918, -qJD(3) * t1114 + qJD(4) * t299 + t1116 - t915, -qJD(3) * t121 - qJD(4) * t120 - t923, -qJD(3) * t32 - qJD(4) * t34 - qJD(5) * t161 - t988; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t818, t820, 0, 0, 0, 0, 0, 0, 0, 0, t1099, -t1118, 0, -t843, 0, 0, 0, 0, 0, 0, t785, -t838, t831, -t940, 0, 0, 0, 0, 0, 0, -t848, -t1144, -t926, t808 - t1020; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t829, -t1118, 0, 0, 0, 0, 0, 0, 0, 0, -t833, t835, t863, -t878, 0, 0, 0, 0, 0, 0, t847, t844, -t927, t718; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t851, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t871; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t1148, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t787, -t828, 0, t787, 0, 0, t749 * t613, t749 * t611, 0, 0, t789, t1167, 0, -t1153, -t824, 0, -t897 + t906, t905 + t1117, 0, qJD(2) * t340 - t876, t641 * t789 + t889, qJD(4) * t216 + t1135, -qJD(4) * t390 - t1185, t639 * t789 - t889, -t381 + t1184, -t1153, qJD(4) * t111 + qJD(5) * t388 - t645 * t897 - t965, qJD(2) * t1113 + qJD(4) * t114 + qJD(5) * t396 - t963, -qJD(2) * t397 + qJD(4) * t56 + t970, qJD(2) * t82 + qJD(4) * t26 - t1146 - t982, t1163, qJD(4) * t1193 + t1187, -qJD(4) * t298 - t1116 - t1183, qJD(4) * t185 + t1138, qJD(4) * t291 - t1168 - t240, -t1166, qJD(2) * t1115 + qJD(4) * t52 - qJD(5) * t245 - qJD(6) * t41 - t998, qJD(2) * t1114 + qJD(4) * t54 - qJD(5) * t1091 - qJD(6) * t42 - t997, qJD(2) * t121 + qJD(4) * t4 - t1014 + t116, qJD(2) * t32 + qJD(4) * t2 - qJD(5) * t27 - t1010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t818, -t820, 0, 0, 0, 0, 0, 0, 0, 0, -t1099, t1118, 0, t843, 0, 0, 0, 0, 0, 0, -t785, t838, -t831, t940, 0, 0, 0, 0, 0, 0, t848, t1144, t926, t808 + t1020; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t810, -pkin(3) * t754, 0, 0, 0, 0, 0, 0, 0, 0, -t645 * t810, t626, t929, qJD(4) * t548 + qJD(5) * t584, -t788, t472, 0, t788, 0, 0, t623 * t819 + t585, -t599 * t623 + t586, t935, qJD(4) * t242 + qJD(5) * t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t825, 0, -t735 (-t755 - t754) * pkin(3), 0, 0, t839, t864, -t836, -t839, -t834, 0, -t645 * t735 + t877, t626 - t683, -t717 + t929 (-pkin(4) * t648 + qJ(5) * t676) * t1028 + t542 * qJD(5) + t720, -t788, t472 + t925, -t845, t788 + t914, t849, 0, t415 * qJD(6) + t585 - t685, t414 * qJD(6) + t586 - t684, -t727 + t935 (t1034 * t634 - t582 * t571 + t583 * t573) * qJD(4) + t177 * qJD(5) + t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t837, t832, t1092, qJD(4) * t542 - t716, 0, 0, 0, 0, 0, 0, -t857, -t853, t109, qJD(4) * t177 - t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1143, t128, -t1148, t1143, t221, -t1100, qJD(4) * t415 - qJD(6) * t546 - t690, qJD(4) * t414 + qJD(6) * t545 - t689, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1153, t1167, 0, -t1153, t555 * qJD(3), 0, -qJD(2) * t442 - t796, t797 + t1117, 0, 0, t1153 * t641 - t896, -qJD(3) * t216 + t1135, qJD(3) * t390 - t1185, t1153 * t639 + t896, qJD(3) * t393 + t1184, -t1153, qJD(2) * t394 - qJD(3) * t111 - t962, -qJD(2) * t392 - qJD(3) * t114 - t961, -qJD(2) * t219 - qJD(3) * t56 + t969, qJD(2) * t104 - qJD(3) * t26 - t1146 - t971, t1163, -qJD(3) * t1193 + t1187, qJD(3) * t298 - qJD(6) * t247 - t1183, -qJD(3) * t185 + t1138, -qJD(3) * t291 - qJD(6) * t244 - t1168, -t1166, -qJD(2) * t741 - qJD(3) * t52 - qJD(5) * t243 - qJD(6) * t45 - t996, -qJD(2) * t299 - qJD(3) * t54 - qJD(5) * t248 - qJD(6) * t46 - t995, qJD(2) * t120 - qJD(3) * t4 - t1011 + t116, qJD(2) * t34 - qJD(3) * t2 - qJD(5) * t29 - t1009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t829, t1118, 0, 0, 0, 0, 0, 0, 0, 0, t833, -t835, -t863, t878, 0, 0, 0, 0, 0, 0, -t847, -t844, t927, -t718; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t825, 0, t811, pkin(3) * t755, 0, 0, -t839, -t864, t836, t839, t834, 0, t645 * t811 - t877, t683, t625 + t717, -qJD(5) * t541 - t720, -t788, t472 - t925, t845, t788 - t914, -t849, 0, -t412 * qJD(6) + t685, -t413 * qJD(6) + t684, t562 + t727, -qJD(5) * t176 - t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t625, t622 * qJD(5), -t788, t472, 0, t788, 0, 0, t634 * t819, -t634 * t599, t562, qJD(5) * t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1092, -t670, 0, 0, 0, 0, 0, 0, -t859, -t855, t109, t671; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1143, t128, -t599 - t856, t1143, -t819 - t858, -t1100, -qJD(6) * t573 - t669, qJD(6) * t571 - t668, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t388 + t744, -qJD(3) * t396 - t1153 * t643, -t862, -qJD(2) * t271 + t1106 * t814 + t928, 0, 0, 0, 0, 0, 0, -qJD(6) * t212 - t714 - t800, qJD(4) * t248 - t1086 + t1103 - t798, -t1081, qJD(2) * t161 + qJD(3) * t27 + qJD(4) * t29 - t1019; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t851, 0, 0, 0, 0, 0, 0, 0, 0, 0, t871; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t837, -t832, -t1092, qJD(4) * t541 + t716, 0, 0, 0, 0, 0, 0, t819 + t857, -t599 + t853, -t109, qJD(4) * t176 + t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1092, t670, 0, 0, 0, 0, 0, 0, t819 + t859, -t599 + t855, -t109, -t671; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1104, -t1077, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1097, -t1080, qJD(4) * t247 + t1096, t1097, t1118 * t378 + t713, t400, qJD(2) * t246 + qJD(3) * t41 + qJD(4) * t45 + qJD(5) * t212 - t1017, -qJD(2) * t249 + qJD(3) * t42 + qJD(4) * t46 + qJD(5) * t1114 - t1018, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239, -t1108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1143, -t128, t1108, -t1143, t239, t1100, qJD(4) * t412 - t603 + t690, qJD(4) * t413 + t598 + t689, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1143, -t128, t856, -t1143, t858, t1100, -t603 + t669, t598 + t668, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1104, t1077, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;
