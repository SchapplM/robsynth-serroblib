% Calculate inertial parameters regressor of coriolis matrix for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRRRPR3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:07
% EndTime: 2019-03-09 22:05:02
% DurationCPUTime: 46.64s
% Computational Cost: add. (34673->1060), mult. (65812->1224), div. (0->0), fcn. (76331->8), ass. (0->761)
t966 = qJD(2) + qJD(3);
t714 = qJD(4) + t966;
t1182 = sin(qJ(3));
t726 = sin(qJ(2));
t1184 = cos(qJ(3));
t728 = cos(qJ(2));
t950 = t1184 * t728;
t667 = -t1182 * t726 + t950;
t725 = sin(qJ(4));
t1101 = t725 * t667;
t1183 = cos(qJ(4));
t790 = t1182 * t728 + t1184 * t726;
t1240 = t1183 * t790;
t1258 = t1240 + t1101;
t1219 = t1258 ^ 2;
t1249 = t725 * t790;
t656 = t1183 * t667;
t553 = t656 - t1249;
t547 = t553 ^ 2;
t1231 = t547 - t1219;
t727 = cos(qJ(6));
t1301 = t1231 * t727;
t1306 = qJD(1) * t1301;
t724 = sin(qJ(6));
t1302 = t1231 * t724;
t1305 = qJD(1) * t1302;
t1216 = -pkin(8) - pkin(7);
t1228 = t1216 * t1182;
t881 = t1184 * t1216;
t573 = t728 * t1228 + t726 * t881;
t1291 = -t790 * pkin(9) + t573;
t1293 = t1183 * t1291;
t1304 = t1293 / 0.2e1;
t1295 = t725 * t1291;
t1303 = -t1295 / 0.2e1;
t1300 = t1231 * qJD(1);
t1285 = t727 * t553;
t1292 = t1285 / 0.2e1;
t1298 = 0.2e1 * t1292;
t340 = t1298 * qJD(5);
t1276 = t724 * t1258;
t1286 = t1276 / 0.2e1;
t513 = -t1276 / 0.2e1;
t1290 = t513 + t1286;
t1299 = qJD(4) * t1298 + qJD(6) * t1290;
t668 = t728 * t881;
t854 = t726 * t1228;
t571 = t668 - t854;
t857 = t1216 * t950;
t574 = -t857 + t854;
t1297 = -t574 - t571;
t876 = pkin(2) * t1184 + pkin(3);
t821 = t725 * t876;
t874 = t1182 * t1183;
t858 = pkin(2) * t874;
t646 = t858 + t821;
t723 = t727 ^ 2;
t1191 = t723 / 0.2e1;
t722 = t724 ^ 2;
t900 = t722 / 0.2e1 + t1191;
t1070 = t900 * t646;
t963 = t1183 * pkin(3);
t704 = -t963 - pkin(4);
t699 = -pkin(10) + t704;
t1296 = t1070 * t699;
t1130 = t553 * qJ(5);
t366 = pkin(4) * t1258 - t1130;
t1289 = t966 * t573;
t1288 = pkin(4) * t553;
t1287 = pkin(5) * t553;
t1278 = qJ(5) * t1258;
t705 = -pkin(2) * t728 - pkin(1);
t611 = -pkin(3) * t667 + t705;
t785 = t611 - t1278;
t316 = t785 - t1288;
t1135 = t316 * t553;
t677 = t1183 * t876;
t951 = t725 * t1182;
t886 = pkin(2) * t951;
t645 = t886 - t677;
t1203 = -t645 / 0.2e1;
t638 = -pkin(4) + t645;
t632 = -pkin(10) + t638;
t1201 = -t646 / 0.2e1;
t635 = qJ(5) + t646;
t1205 = t635 / 0.2e1;
t903 = t1201 + t1205;
t868 = t903 * t1258;
t757 = t868 - t553 * (t632 / 0.2e1 + t1203);
t1217 = pkin(4) + pkin(10);
t1186 = t1217 / 0.2e1;
t910 = t553 * t1186;
t928 = t1278 / 0.2e1;
t812 = t928 + t910;
t1284 = t757 - t812;
t1039 = qJD(1) * t1258;
t1283 = t553 * t1039;
t1189 = t724 / 0.2e1;
t1209 = -t1258 / 0.2e1;
t1208 = t1258 / 0.2e1;
t1178 = t1258 * pkin(5);
t1175 = t725 * pkin(3);
t711 = qJD(4) * t1175;
t1176 = t667 * pkin(9);
t743 = t574 + t1176;
t1241 = t1183 * t743;
t1265 = t1241 / 0.2e1;
t734 = (-pkin(9) * t1182 + t1228) * t728 + (-pkin(9) * t1184 + t881) * t726;
t732 = t725 * t734;
t731 = t732 / 0.2e1 + t1265;
t863 = -t1241 / 0.2e1 + t1303;
t194 = t731 + t863;
t1005 = t194 * qJD(1);
t1174 = pkin(3) * qJD(3);
t672 = t821 / 0.2e1;
t706 = t1175 / 0.2e1;
t1069 = t672 + t706;
t952 = t725 * t1184;
t887 = pkin(2) * t952;
t590 = -t887 / 0.2e1 + t1069;
t956 = t590 * qJD(2) + t725 * t1174 + t1005;
t1280 = t711 + t956;
t1062 = t722 + t723;
t1188 = t725 / 0.2e1;
t700 = qJ(5) + t1175;
t1195 = t700 / 0.2e1;
t1197 = t699 / 0.2e1;
t964 = -t1183 / 0.2e1;
t872 = t553 * t964;
t740 = pkin(3) * (t1188 * t1258 - t872) - t1195 * t1258 + t1197 * t553;
t1279 = t740 + t812;
t1277 = t316 * t1258;
t484 = t571 - t1176;
t1104 = t725 * t484;
t1075 = t1104 / 0.2e1 + t1304;
t1250 = t725 * t743;
t733 = -t1250 / 0.2e1 + t1304;
t192 = t733 - t1075;
t875 = t1184 * t1183;
t661 = (t875 - t951) * pkin(2);
t974 = t661 * qJD(2);
t1080 = t192 * qJD(1) + t974;
t884 = t963 / 0.2e1;
t823 = t677 / 0.2e1 + t884;
t859 = pkin(2) * t875;
t591 = -t859 / 0.2e1 + t823;
t1274 = -qJD(4) * t591 + t1080;
t465 = t1183 * t484;
t1210 = t465 / 0.2e1;
t186 = t1265 + t1210;
t660 = (t952 + t874) * pkin(2);
t1084 = t186 * qJD(1) + t660 * qJD(2);
t1273 = qJD(4) * t590 - t1084;
t562 = -t886 + t859 / 0.2e1 + t823;
t1272 = qJD(4) * t562 + t1080;
t466 = t1183 * t734;
t736 = -t466 / 0.2e1 + t1250 / 0.2e1;
t189 = t736 + t1075;
t980 = t645 * qJD(2);
t1082 = t189 * qJD(1) + t980;
t1271 = qJD(3) * t591 - t1082;
t196 = t1303 + t1210 + t731;
t1004 = t196 * qJD(1);
t1078 = t646 * qJD(2) + t1004;
t1270 = qJD(3) * t590 + t1078;
t681 = t887 / 0.2e1;
t561 = t858 + t681 + t1069;
t977 = t646 * qJD(4);
t1269 = qJD(3) * t561 + t1078 + t977;
t979 = t645 * qJD(4);
t1268 = -qJD(3) * t562 + t1082 + t979;
t975 = t660 * qJD(3);
t1267 = qJD(4) * t561 + t1084 + t975;
t1251 = t714 * t727;
t1263 = t724 * t1251;
t1192 = -t722 / 0.2e1;
t899 = t1191 + t1192;
t353 = t899 * t553;
t291 = qJD(1) * t353 + t1263;
t695 = t722 - t723;
t364 = t695 * t547;
t1266 = qJD(1) * t364 - 0.2e1 * t1263 * t553;
t751 = t1240 / 0.2e1;
t763 = t1249 / 0.2e1;
t320 = t1250 - t466;
t735 = t320 + t1178;
t1144 = t724 * t735;
t1143 = t727 * t735;
t1264 = t1217 * t553;
t1207 = -t632 / 0.2e1;
t1196 = -t700 / 0.2e1;
t1198 = -t661 / 0.2e1;
t1199 = t660 / 0.2e1;
t803 = -t1198 * t553 + t1199 * t1258;
t749 = -t1258 * (t1196 + t1205) + t803;
t1260 = t749 - t553 * (t1197 + t1207);
t423 = 0.2e1 * t763 - t656;
t889 = t966 * t553;
t246 = -qJD(4) * t423 + t889;
t424 = 0.2e1 * t751 + t1101;
t890 = t966 * t1258;
t244 = qJD(4) * t424 + t890;
t241 = t785 - t1264;
t134 = t241 * t724 - t1143;
t135 = t727 * t241 + t1144;
t1245 = t135 * t1189 - t134 * t727 / 0.2e1;
t1257 = t714 * t1231;
t1252 = t695 * t714;
t39 = (t134 * t724 + t135 * t727) * t1258;
t546 = t751 - t1240 / 0.2e1;
t1248 = t966 * t546;
t1068 = t1062 * t706;
t773 = t790 * pkin(3);
t324 = t773 + t366;
t255 = pkin(10) * t1258 + t324;
t1097 = t727 * t255;
t318 = t1241 + t1295;
t249 = t318 + t1287;
t1111 = t724 * t249;
t143 = t1097 + t1111;
t1154 = t143 * t724;
t1099 = t727 * t249;
t1109 = t724 * t255;
t142 = t1099 - t1109;
t1155 = t142 * t727;
t1086 = t1155 / 0.2e1 + t1154 / 0.2e1;
t730 = t1287 / 0.2e1 + t731;
t29 = t730 - t1086;
t707 = -t1175 / 0.2e1;
t1227 = qJ(5) + t858 / 0.2e1;
t673 = -t821 / 0.2e1;
t782 = t673 - t1227;
t867 = t900 * t660;
t380 = t707 + t867 + t782;
t971 = t700 * qJD(3);
t1246 = -qJD(1) * t29 + qJD(2) * t380 - t971;
t1244 = t975 + t977;
t1015 = qJD(6) * t727;
t698 = t724 * t1015;
t330 = (t1208 + t1209) * t727 * t724;
t996 = t330 * qJD(1);
t312 = t996 - t698;
t311 = t996 + t698;
t1242 = t1298 * t966;
t319 = -t465 + t1295;
t322 = t1293 + t1104;
t323 = t1241 + t732;
t840 = -t1258 * t323 + t320 * t553;
t761 = t1258 * t319 + t322 * t553 + t840;
t1239 = t761 * qJD(1);
t321 = t1293 - t1250;
t762 = t1258 * t318 + t321 * t553 + t840;
t1238 = t762 * qJD(1);
t981 = t635 * qJD(2);
t1237 = -t1004 - t981;
t993 = t353 * qJD(6);
t1235 = qJD(4) * t330 - t993;
t1095 = t727 * t1258;
t913 = t1095 / 0.2e1;
t329 = t1286 * t727 + t724 * t913;
t1234 = -qJD(4) * t329 + t993;
t720 = t726 * pkin(2);
t248 = t720 + t255;
t1100 = t727 * t248;
t250 = t319 + t1287;
t1110 = t724 * t250;
t141 = t1100 + t1110;
t1156 = t141 * t724;
t1098 = t727 * t250;
t1112 = t724 * t248;
t140 = t1098 - t1112;
t1157 = t140 * t727;
t811 = t1156 / 0.2e1 + t1157 / 0.2e1;
t27 = t730 - t811;
t1232 = -qJD(1) * t27 - t981;
t1023 = qJD(4) * t1258;
t1230 = t890 + t1023;
t281 = t1217 * t1258 - t1130;
t1096 = t727 * t281;
t254 = t323 + t1287;
t1137 = t254 * t724;
t153 = t1096 + t1137;
t1152 = t153 * t724;
t1108 = t724 * t281;
t240 = t254 * t727;
t152 = t240 - t1108;
t1153 = t152 * t727;
t843 = t1152 + t1153;
t193 = t731 - t863;
t1073 = t1295 / 0.2e1 - t465 / 0.2e1;
t195 = t731 + t1073;
t957 = t195 * qJD(2) + t193 * qJD(3) + t323 * qJD(4);
t187 = t733 + t736;
t710 = qJD(3) * t963;
t959 = t187 * qJD(1) - t591 * qJD(2) - t710;
t1065 = t707 - qJ(5);
t582 = t681 + t673 + t1065;
t1226 = qJD(2) * t582 - t1005 - t971;
t1225 = -t966 * t329 + t993;
t1224 = -t330 * t966 - t993;
t1035 = qJD(1) * t727;
t938 = t724 * t1035;
t113 = -t353 * t714 + t547 * t938;
t715 = t724 * qJD(4);
t1190 = -t724 / 0.2e1;
t918 = t553 * t1190;
t1106 = t724 * t553;
t919 = -t1106 / 0.2e1;
t800 = t919 + t918;
t1223 = -t553 * t715 + t800 * t966;
t1222 = t424 * t966 + t1023;
t1024 = qJD(4) * t553;
t1221 = t423 * t966 - t1024;
t1220 = qJD(2) * t761 + qJD(3) * t762;
t1218 = pkin(4) / 0.2e1;
t1215 = -qJ(5) / 0.2e1;
t1214 = qJ(5) / 0.2e1;
t1213 = t240 / 0.2e1;
t1212 = -t735 / 0.2e1;
t1211 = -t254 / 0.2e1;
t1206 = -t635 / 0.2e1;
t1204 = t638 / 0.2e1;
t1202 = t645 / 0.2e1;
t1200 = -t660 / 0.2e1;
t1194 = -t704 / 0.2e1;
t1193 = t704 / 0.2e1;
t1187 = t727 / 0.2e1;
t1185 = -t1217 / 0.2e1;
t1173 = qJD(2) * pkin(2);
t253 = t322 - t1178;
t7 = -t134 * t140 + t135 * t141 + t253 * t254;
t1172 = t7 * qJD(1);
t252 = t321 - t1178;
t8 = -t134 * t142 + t135 * t143 + t252 * t254;
t1171 = t8 * qJD(1);
t9 = -t134 * t152 + t135 * t153 - t254 * t735;
t1170 = t9 * qJD(1);
t1166 = qJD(1) * t39;
t82 = t1258 * t134 - t1285 * t254;
t1165 = qJD(1) * t82;
t83 = -t1106 * t254 - t1258 * t135;
t1164 = qJD(1) * t83;
t870 = t1095 * t135 + t1276 * t134;
t10 = -(-t140 * t724 + t141 * t727) * t553 + t870;
t1163 = t10 * qJD(1);
t11 = -(-t142 * t724 + t143 * t727) * t553 + t870;
t1162 = t11 * qJD(1);
t12 = t39 - (-t152 * t724 + t153 * t727) * t553;
t1161 = t12 * qJD(1);
t1159 = t135 * t553;
t893 = -t254 * t1095 - t134 * t553;
t20 = t1258 * t140 + t1285 * t253 + t893;
t1151 = t20 * qJD(1);
t1138 = t254 * t1258;
t21 = -t1159 - t141 * t1258 + (-t253 * t553 + t1138) * t724;
t1150 = t21 * qJD(1);
t22 = t1258 * t142 + t1285 * t252 + t893;
t1149 = t22 * qJD(1);
t23 = -t1159 - t143 * t1258 + (-t252 * t553 + t1138) * t724;
t1148 = t23 * qJD(1);
t24 = (t152 - t240) * t1258 - (t134 + t1143) * t553;
t1147 = t24 * qJD(1);
t25 = (-t153 + t1137) * t1258 - (t135 - t1144) * t553;
t1146 = t25 * qJD(1);
t1145 = t735 * t700;
t1142 = t252 * t724;
t1141 = t252 * t727;
t1140 = t253 * t724;
t1139 = t253 * t727;
t317 = t720 + t324;
t841 = t319 * t320 + t322 * t323;
t52 = t316 * t317 + t841;
t1133 = t52 * qJD(1);
t842 = t320 * t318 + t323 * t321;
t55 = t316 * t324 + t842;
t1132 = t55 * qJD(1);
t1131 = t1258 * t725;
t1129 = t553 * t699;
t63 = t316 * t366;
t1128 = t63 * qJD(1);
t1127 = t635 * t1258;
t1126 = t635 * t645;
t1125 = t635 * t661;
t1124 = t638 * t553;
t1123 = t645 * qJ(5);
t1122 = t645 * t553;
t1121 = t645 * t700;
t1120 = t646 * t1258;
t1119 = t661 * qJ(5);
t1118 = t661 * t700;
t1117 = t700 * t1258;
t1116 = t704 * t553;
t1115 = t705 * t667;
t615 = t720 + t773;
t87 = t611 * t615 + t841;
t1088 = t87 * qJD(1);
t88 = t611 * t773 + t842;
t1087 = t88 * qJD(1);
t810 = -t1153 / 0.2e1 - t1152 / 0.2e1;
t1077 = -t196 * qJD(2) - t194 * qJD(3);
t1071 = -t1278 / 0.2e1 - t1288 / 0.2e1;
t712 = qJD(4) * t963;
t718 = qJD(5) * t724;
t1067 = t724 * t712 + t718;
t719 = qJD(5) * t727;
t1066 = t727 * t712 + t719;
t1064 = -t645 * t715 + t718;
t1063 = -t727 * t979 + t719;
t1061 = qJ(5) * qJD(4);
t144 = t317 * t553 - t1277;
t1060 = qJD(1) * t144;
t145 = -t1258 * t317 - t1135;
t1059 = qJD(1) * t145;
t148 = t324 * t553 - t1277;
t1058 = qJD(1) * t148;
t149 = -t1258 * t324 - t1135;
t1057 = qJD(1) * t149;
t158 = t366 * t553 - t1277;
t1056 = qJD(1) * t158;
t159 = -t1258 * t366 - t1135;
t1055 = qJD(1) * t159;
t400 = t611 * t1258;
t265 = -t553 * t615 + t400;
t1051 = qJD(1) * t265;
t401 = t611 * t553;
t266 = t1258 * t615 + t401;
t1050 = qJD(1) * t266;
t286 = t553 * t773 - t400;
t1047 = qJD(1) * t286;
t287 = -t1258 * t773 - t401;
t1046 = qJD(1) * t287;
t524 = -t667 * t720 + t705 * t790;
t1042 = qJD(1) * t524;
t525 = t720 * t790 + t1115;
t1041 = qJD(1) * t525;
t1040 = qJD(1) * t553;
t1038 = qJD(1) * t611;
t1037 = qJD(1) * t667;
t1036 = qJD(1) * t724;
t1034 = qJD(1) * t728;
t1018 = qJD(4) * t611;
t1017 = qJD(4) * t727;
t1016 = qJD(6) * t724;
t1014 = qJD(6) * t1217;
t100 = (t1204 + t1194) * t553 + t749;
t1013 = t100 * qJD(1);
t904 = t1202 - t638 / 0.2e1;
t758 = -t553 * t904 - t868;
t818 = t1288 / 0.2e1 + t928;
t101 = t758 + t818;
t1012 = t101 * qJD(1);
t744 = t1122 / 0.2e1 - t1120 / 0.2e1 + t803;
t947 = t1183 * t553;
t774 = (-t1131 / 0.2e1 - t947 / 0.2e1) * pkin(3);
t105 = t774 - t744;
t1011 = t105 * qJD(1);
t871 = t1196 + t706;
t755 = -pkin(3) * t872 - t1194 * t553 + t1258 * t871;
t156 = t755 + t818;
t1010 = t156 * qJD(1);
t201 = t1297 * t790;
t1003 = t201 * qJD(1);
t209 = t1062 * t1258 * t553;
t1002 = t209 * qJD(1);
t256 = -t1297 * t573 + t705 * t720;
t999 = t256 * qJD(1);
t907 = 0.2e1 * t1208;
t345 = t907 * t724;
t332 = t345 * qJD(1);
t351 = 0.2e1 * t513;
t994 = t351 * qJD(1);
t354 = t907 * t727;
t339 = t354 * qJD(1);
t360 = t1292 - t1285 / 0.2e1;
t992 = t360 * qJD(1);
t481 = t667 ^ 2 - t790 ^ 2;
t991 = t481 * qJD(1);
t544 = t763 - t656 / 0.2e1;
t990 = t544 * qJD(1);
t545 = t656 / 0.2e1 - t1249 / 0.2e1;
t989 = t545 * qJD(1);
t988 = t546 * qJD(1);
t987 = t546 * qJD(4);
t986 = t1219 * qJD(1);
t985 = t553 * qJD(5);
t813 = t857 / 0.2e1;
t568 = -t668 / 0.2e1 + t813;
t984 = t568 * qJD(1);
t973 = t661 * qJD(3);
t696 = -t726 ^ 2 + t728 ^ 2;
t972 = t696 * qJD(1);
t717 = t724 * qJD(2);
t716 = t724 * qJD(3);
t970 = t726 * qJD(2);
t969 = t728 * qJD(2);
t968 = qJD(5) + t973;
t967 = t712 + qJD(5);
t962 = pkin(1) * t726 * qJD(1);
t961 = pkin(1) * t1034;
t960 = t672 + t1227;
t954 = t715 + t716 + t717;
t953 = t635 * t1183;
t946 = t1183 * t700;
t945 = t316 * t1039;
t944 = t553 * t1038;
t943 = t722 * t1040;
t942 = t723 * t1040;
t941 = t1258 * t1038;
t940 = t1219 * t1036;
t939 = t705 * t1037;
t936 = t727 * t715;
t935 = t1258 * t1016;
t934 = t1258 * t1015;
t933 = t1258 * t1040;
t930 = t726 * t969;
t929 = t1258 * t1035;
t925 = t253 * t1214;
t924 = t253 * t1195;
t923 = -t1137 / 0.2e1;
t922 = -t1121 / 0.2e1;
t921 = -t1117 / 0.2e1;
t916 = t632 * t1189;
t915 = t660 * t1190;
t914 = -t1095 / 0.2e1;
t912 = t632 * t1187;
t911 = t660 * t1187;
t909 = t1212 - t252 / 0.2e1;
t908 = t1212 - t253 / 0.2e1;
t898 = t1184 * qJD(2);
t897 = t1184 * qJD(3);
t896 = t1182 * qJD(2);
t895 = t1182 * qJD(3);
t477 = t1062 * t646;
t526 = t1062 * t660;
t894 = t1062 * t699;
t892 = qJD(6) * t544 - t1283;
t279 = qJD(6) * t545 + t1283;
t891 = t966 * t724;
t888 = t966 * t727;
t885 = qJD(6) + t1039;
t883 = t1062 * t707 - t1070;
t878 = t547 * t698;
t877 = t553 * t938;
t659 = t1062 * t1175;
t873 = t953 / 0.2e1;
t865 = 0.2e1 * t877;
t864 = -0.2e1 * t877;
t862 = t1215 - t903;
t861 = t1199 + t1196 + t1206;
t852 = t1258 * t877;
t851 = t724 * t910 + qJ(5) * t1286 - t1143 / 0.2e1;
t850 = -t1217 * t1292 + qJ(5) * t914 - t1144 / 0.2e1;
t849 = t1278 + t1264;
t808 = t1198 * t254 + t1206 * t252;
t1 = t924 + (t1197 * t140 + t1199 * t134 + t1207 * t142) * t727 + (t1197 * t141 + t1200 * t135 + t1207 * t143) * t724 + t808;
t267 = t526 * t632 + t1125;
t848 = -t1 * qJD(1) + t267 * qJD(2);
t257 = t477 * t632 - t1126;
t809 = t1202 * t254 - t1206 * t735;
t3 = t925 + (t140 * t1185 + t152 * t1207 + t134 * t646 / 0.2e1) * t727 + (t1185 * t141 + t1201 * t135 + t1207 * t153) * t724 + t809;
t847 = -t3 * qJD(1) + t257 * qJD(2);
t845 = t1156 + t1157;
t844 = t1154 + t1155;
t838 = -t553 * t632 + t1127;
t837 = t1117 - t1129;
t13 = -t811 + t1086;
t836 = qJD(1) * t13 + qJD(2) * t526;
t16 = (-t140 / 0.2e1 + t152 / 0.2e1) * t727 + (-t141 / 0.2e1 + t153 / 0.2e1) * t724;
t835 = qJD(1) * t16 + qJD(2) * t477;
t806 = t320 * t1199 + t323 * t661 / 0.2e1;
t746 = t1204 * t318 + t1205 * t321 + t806;
t807 = t1194 * t319 + t1196 * t322;
t31 = t746 + t807;
t376 = t638 * t660 + t1125;
t833 = t31 * qJD(1) + t376 * qJD(2);
t760 = -t320 * t903 - t323 * t904;
t816 = t1215 * t322 + t1218 * t319;
t32 = t760 + t816;
t365 = t638 * t646 - t1126;
t832 = t32 * qJD(1) + t365 * qJD(2);
t745 = t1201 * t321 + t1203 * t318 - t806;
t775 = (t1188 * t322 + t319 * t964) * pkin(3);
t34 = t775 + t745;
t377 = t645 * t660 + t646 * t661;
t831 = -t34 * qJD(1) + t377 * qJD(2);
t830 = t885 * t724;
t829 = t885 * t727;
t828 = qJD(2) * t189 + qJD(3) * t187;
t827 = -qJD(2) * t186 + qJD(4) * t194;
t826 = -qJD(2) * t192 - qJD(4) * t187;
t825 = qJD(3) * t186 + qJD(4) * t196;
t824 = qJD(3) * t192 - qJD(4) * t189;
t822 = t884 + t1203;
t817 = t1215 * t321 + t1218 * t318;
t815 = pkin(4) * t1199 - t1119 / 0.2e1;
t814 = t1215 + t871;
t805 = t1129 / 0.2e1 + t921;
t804 = -t1205 * t553 + t1207 * t1258;
t801 = -t1195 * t553 + t1209 * t699;
t799 = t1258 * t1186 - t1130 / 0.2e1;
t234 = t1139 / 0.2e1;
t43 = t234 - t1141 / 0.2e1 + t1260 * t724;
t798 = qJD(1) * t43 - t727 * t974;
t233 = t1140 / 0.2e1;
t45 = t233 - t1142 / 0.2e1 - t1260 * t727;
t797 = qJD(1) * t45 - t661 * t717;
t49 = t1284 * t724 + t908 * t727;
t796 = -qJD(1) * t49 + t727 * t980;
t51 = -t1284 * t727 + t908 * t724;
t795 = -qJD(1) * t51 + t645 * t717;
t789 = t248 / 0.2e1 + t804;
t68 = t1213 - t1098 / 0.2e1 + t789 * t724;
t794 = -qJD(1) * t68 - t727 * t981;
t70 = (t1211 + t250 / 0.2e1) * t724 + t789 * t727;
t793 = -qJD(1) * t70 + t635 * t717;
t792 = t1198 + t822;
t788 = t255 / 0.2e1 + t801;
t787 = t281 / 0.2e1 + t799;
t766 = -t1217 * t867 + t1119 / 0.2e1;
t160 = t1121 / 0.2e1 - t1296 + (-t953 / 0.2e1 - t900 * t725 * t632) * pkin(3) + t766;
t750 = -t1086 * t1217 + t1214 * t252;
t5 = t1145 / 0.2e1 + t810 * t699 + (-t1245 * t725 + t254 * t964) * pkin(3) + t750;
t503 = (t725 * t894 + t946) * pkin(3);
t786 = -t5 * qJD(1) - t160 * qJD(2) + t503 * qJD(3);
t18 = -t810 - t1086;
t296 = t867 + t883;
t784 = qJD(1) * t18 - qJD(2) * t296 + qJD(3) * t659;
t738 = (t1188 * t638 + t873) * pkin(3) + t922 + t646 * t1193;
t220 = t738 + t815;
t756 = t320 * t871 + (t1193 + t884) * t323;
t40 = t756 + t817;
t610 = -pkin(3) * t946 - t1175 * t704;
t783 = t40 * qJD(1) + t220 * qJD(2) - t610 * qJD(3);
t185 = -t863 + t1073;
t781 = qJD(2) * t185 + qJD(3) * t318 + qJD(4) * t193;
t780 = qJD(2) * t319 + qJD(3) * t185 + qJD(4) * t195;
t188 = -t733 + t736;
t190 = t736 - t1075;
t779 = -qJD(2) * t190 - qJD(3) * t188 - qJD(4) * t320;
t191 = t733 + t1075;
t778 = qJD(2) * t191 + qJD(3) * t321 - qJD(4) * t188;
t777 = qJD(2) * t322 + qJD(3) * t191 - qJD(4) * t190;
t164 = (t889 + t1024) * t1258;
t776 = t1230 * t553;
t419 = t861 * t724;
t80 = (t1211 + t249 / 0.2e1) * t724 + t788 * t727;
t770 = -qJD(1) * t80 - qJD(2) * t419 + t700 * t716;
t420 = t861 * t727;
t78 = t1213 - t1099 / 0.2e1 + t788 * t724;
t769 = -qJD(1) * t78 + qJD(2) * t420 - t727 * t971;
t768 = qJD(1) * t790;
t767 = t790 * qJD(3);
t764 = t705 * t768;
t37 = t730 + t810;
t382 = t782 + t1070;
t596 = t1065 + t1068;
t754 = qJD(1) * t37 - qJD(2) * t382 - qJD(3) * t596 + t1061;
t434 = t792 * t724;
t67 = t1279 * t727 + t909 * t724;
t753 = -t67 * qJD(1) - t434 * qJD(2) - t710 * t724;
t436 = t792 * t727;
t65 = -t1279 * t724 + t909 * t727;
t752 = -t65 * qJD(1) - t436 * qJD(2) - t710 * t727;
t166 = t787 * t724;
t430 = t862 * t727;
t584 = t814 * t727;
t748 = -qJ(5) * t1017 - t166 * qJD(1) + t430 * qJD(2) + t584 * qJD(3);
t167 = t787 * t727;
t429 = t862 * t724;
t583 = t814 * t724;
t747 = qJ(5) * t715 - t167 * qJD(1) - t429 * qJD(2) - t583 * qJD(3);
t560 = t966 * t790;
t721 = qJ(5) * qJD(5);
t713 = qJ(5) * t963;
t702 = qJ(5) * t1187;
t701 = qJ(5) * t1190;
t697 = t726 * t1034;
t692 = t700 * qJD(5);
t676 = t695 * qJD(6);
t671 = t714 * qJ(5);
t670 = t700 * t1187;
t669 = t700 * t1190;
t647 = t659 * qJD(4);
t623 = t727 * t973;
t622 = t661 * t716;
t612 = t635 * qJD(5);
t597 = qJ(5) + t706 + t1068;
t593 = t635 * t1187;
t592 = t635 * t1190;
t586 = t706 * t727 + t670 + t702;
t585 = t707 * t724 + t669 + t701;
t566 = t667 * t768;
t565 = t790 * t1037;
t559 = t966 * t667;
t558 = qJ(5) + t561;
t502 = t526 * qJD(3);
t492 = t668 / 0.2e1 + t813 - t854;
t478 = 0.2e1 * t553 * t698;
t468 = t477 * qJD(4);
t435 = t1187 * t661 + t727 * t822;
t433 = t1189 * t661 + t724 * t822;
t432 = t1187 * t646 + t593 + t702;
t431 = t1190 * t646 + t592 + t701;
t422 = t670 + t593 + t911;
t421 = t669 + t592 + t915;
t406 = t423 * qJD(5);
t383 = t960 + t1070;
t381 = t706 + t867 + t960;
t368 = t865 + t1252;
t367 = t864 - t1252;
t361 = t1209 * t727 + t913;
t333 = t800 * qJD(5);
t326 = -t339 - t1015;
t325 = -t332 - t1016;
t297 = t1192 * t660 + t1200 * t723 + t883;
t258 = -qJD(4) * t544 + t545 * t966;
t232 = t1141 / 0.2e1;
t231 = t1142 / 0.2e1;
t221 = t738 - t815;
t203 = 0.2e1 * t899 * t1258;
t161 = pkin(3) * t873 + t1068 * t632 + t1296 + t766 + t922;
t157 = t755 + t1071;
t155 = t723 * t933 - t1235;
t154 = t722 * t933 + t1235;
t117 = t478 + 0.2e1 * t852;
t106 = t774 + t744;
t102 = t758 + t1071;
t99 = -t1127 / 0.2e1 + t1124 / 0.2e1 + t921 + t1116 / 0.2e1 + t803;
t97 = -qJD(6) * t354 + t1306;
t96 = qJD(4) * t360 - qJD(6) * t345 + t1305;
t95 = (-t727 * t891 - t942) * t1258 + t1234;
t94 = (t724 * t888 - t943) * t1258 - t1234;
t93 = -t1137 - t1096 / 0.2e1 + t799 * t727;
t92 = t240 - t1108 / 0.2e1 + t799 * t724;
t90 = qJD(4) * t800 + qJD(6) * t361 - t553 * t891 - t1306;
t89 = t553 * t888 + t1299 - t1305;
t86 = qJD(4) * t203 - t695 * t890 + t478 - 0.2e1 * t852;
t81 = t923 - t1097 / 0.2e1 - t1111 / 0.2e1 + t801 * t727;
t79 = t1213 - t1109 / 0.2e1 + t1099 / 0.2e1 + t801 * t724;
t71 = t923 - t1100 / 0.2e1 - t1110 / 0.2e1 + t804 * t727;
t69 = t1213 - t1112 / 0.2e1 + t1098 / 0.2e1 + t804 * t724;
t66 = t727 * t740 + t231 + t850;
t64 = -t724 * t740 + t232 + t851;
t50 = -t727 * t757 + t233 + t850;
t48 = t724 * t757 + t234 + t851;
t46 = t1258 * t911 + t1292 * t661 + t553 * t912 + t635 * t914 + t727 * t805 + t231 + t233;
t44 = t1258 * t915 + t1286 * t635 + t632 * t918 + t661 * t919 - t724 * t805 + t232 + t234;
t41 = t756 - t817;
t36 = t730 - t810;
t35 = t775 - t745;
t33 = t760 - t816;
t30 = t746 - t807;
t28 = t730 + t1086;
t26 = t730 + t811;
t19 = t810 - t1086;
t17 = -t811 + t810;
t14 = -t811 - t1086;
t6 = -t1145 / 0.2e1 + t254 * t884 + t750 + t843 * t1197 + t1245 * t1175;
t4 = -t1217 * t811 + t1245 * t646 + t152 * t912 + t153 * t916 - t809 + t925;
t2 = t1245 * t660 + t142 * t912 + t143 * t916 + t699 * t811 - t808 + t924;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t930, t696 * qJD(2), 0, -t930, 0, 0, -pkin(1) * t970, -pkin(1) * t969, 0, 0 (qJD(2) * t790 + t767) * t667, t966 * t481, 0, -t560 * t667, 0, 0, t524 * qJD(2) + t705 * t767, qJD(2) * t525 + qJD(3) * t1115, qJD(2) * t201, qJD(2) * t256, t164, t1257, 0, -t776, 0, 0, qJD(2) * t265 - qJD(3) * t286 + t1018 * t1258, qJD(2) * t266 - qJD(3) * t287 + t1018 * t553, t1220, qJD(2) * t87 + qJD(3) * t88, 0, 0, 0, t164, t1257, -t776, t1220, qJD(2) * t144 + qJD(3) * t148 + qJD(4) * t158 - t1258 * t985, qJD(2) * t145 + qJD(3) * t149 + qJD(4) * t159 + qJD(5) * t1219, qJD(2) * t52 + qJD(3) * t55 + qJD(4) * t63 - qJD(5) * t1277, -t722 * t776 + t878, -0.2e1 * t1230 * t1285 * t724 - qJD(6) * t364, -t1302 * t714 - t553 * t934, -t723 * t776 - t878, -t1301 * t714 + t553 * t935, t164, qJD(2) * t20 + qJD(3) * t22 + qJD(4) * t24 + qJD(6) * t83 + t1219 * t718, qJD(2) * t21 + qJD(3) * t23 + qJD(4) * t25 + qJD(6) * t82 + t1219 * t719, qJD(2) * t10 + qJD(3) * t11 + qJD(4) * t12 + qJD(5) * t209, qJD(2) * t7 + qJD(3) * t8 + qJD(4) * t9 - qJD(5) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t697, t972, t969, -t697, -t970, 0, -pkin(7) * t969 - t962, pkin(7) * t970 - t961, 0, 0, t566, t991, t559, -t565, -t560, 0, qJD(2) * t571 + qJD(3) * t492 + t1042, t1041 - t1289, t1003 + (-t1182 * t790 - t1184 * t667) * t1173, t999 + (t1182 * t573 + t1184 * t571) * t1173, t1283, t1300, t246, -t933, -t244, 0, -t780 + t1051, -t777 + t1050, t1239 + (-t1120 + t1122) * qJD(2) + t106 * qJD(3), t1088 + (t319 * t645 + t322 * t646) * qJD(2) + t35 * qJD(3), 0, -t246, t244, t1283, t1300, -t933, t1239 + (t1124 - t1127) * qJD(2) + t99 * qJD(3) + t102 * qJD(4) - t406, t780 + t1060, t777 + t1059, t1133 + (t319 * t638 + t322 * t635) * qJD(2) + t30 * qJD(3) + t33 * qJD(4) + t195 * qJD(5), t94, t86, t89, t95, t90, t279, t1151 + (-t727 * t838 + t1140) * qJD(2) + t46 * qJD(3) + t50 * qJD(4) + t340 + t69 * qJD(6), t1150 + (t724 * t838 + t1139) * qJD(2) + t44 * qJD(3) + t48 * qJD(4) + t333 + t71 * qJD(6), -qJD(2) * t845 + t14 * qJD(3) + t17 * qJD(4) + t1163, t1172 + (t253 * t635 + t632 * t845) * qJD(2) + t2 * qJD(3) + t4 * qJD(4) + t26 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t566, t991, t559, -t565, -t560, 0, t492 * qJD(2) - t574 * qJD(3) + t764, -t1289 + t939, 0, 0, t1283, t1300, t246, -t933, -t244, 0, -t781 - t1047, -t778 - t1046, t1238 + t106 * qJD(2) + (-t947 - t1131) * t1174, t1087 + t35 * qJD(2) + (-t1183 * t318 + t321 * t725) * t1174, 0, -t246, t244, t1283, t1300, -t933, t1238 + t99 * qJD(2) + (t1116 - t1117) * qJD(3) + t157 * qJD(4) - t406, t781 + t1058, t778 + t1057, t1132 + t30 * qJD(2) + (t318 * t704 + t321 * t700) * qJD(3) + t41 * qJD(4) + t193 * qJD(5), t94, t86, t89, t95, t90, t279, t1149 + t46 * qJD(2) + (-t727 * t837 + t1142) * qJD(3) + t66 * qJD(4) + t340 + t79 * qJD(6), t1148 + t44 * qJD(2) + (t724 * t837 + t1141) * qJD(3) + t64 * qJD(4) + t333 + t81 * qJD(6), t14 * qJD(2) - qJD(3) * t844 + t19 * qJD(4) + t1162, t1171 + t2 * qJD(2) + (t252 * t700 + t699 * t844) * qJD(3) + t6 * qJD(4) + t28 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1283, t1300, -t1221, -t1283, -t1222, 0, t941 - t957, -t779 + t944, 0, 0, 0, t1221, t1222, t1283, t1300, -t1283, t102 * qJD(2) + t157 * qJD(3) + (-t1278 - t1288) * qJD(4) + t985, t957 + t1056, t779 + t1055, t1128 + t33 * qJD(2) + t41 * qJD(3) + (-pkin(4) * t323 - qJ(5) * t320) * qJD(4) + t323 * qJD(5) (t936 - t943) * t1258 - t1225, t478 + t966 * t203 + (-qJD(4) * t695 + t864) * t1258, t1017 * t553 + t1242 - t1305 (-t936 - t942) * t1258 + t1225, -t1306 + t1223, -t892, t1147 + t50 * qJD(2) + t66 * qJD(3) + (-t727 * t849 - t1144) * qJD(4) + t340 + t92 * qJD(6), -t735 * t1017 + t1146 + t48 * qJD(2) + t64 * qJD(3) + t93 * qJD(6) + (qJD(4) * t849 - t985) * t724, t17 * qJD(2) + t19 * qJD(3) - qJD(4) * t843 + t1161, t1170 + t4 * qJD(2) + t6 * qJD(3) + (-qJ(5) * t735 - t1217 * t843) * qJD(4) + t36 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1221, -t1283, t986, -t945 + t957, 0, 0, 0, 0, 0, 0, t1242 + t1299 + t940, t727 * t986 + t1223, t1002, qJD(2) * t26 + qJD(3) * t28 + qJD(4) * t36 - t1166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, -t1266, -t1285 * t885 + t1290 * t966, -t113, t361 * t966 + t553 * t830, t258, qJD(2) * t69 + qJD(3) * t79 + qJD(4) * t92 + qJD(5) * t1290 - qJD(6) * t135 + t1164, qJD(2) * t71 + qJD(3) * t81 + qJD(4) * t93 + qJD(6) * t134 + t1165, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t697, -t972, 0, t697, 0, 0, t962, t961, 0, 0, -t566, -t991, 0, t565, 0, 0, qJD(3) * t568 - t1042, -t1041, -t1003, -t999, -t1283, -t1300, 0, t933, -t987, 0, -t825 - t1051, -t824 - t1050, -qJD(3) * t105 - t1239, -qJD(3) * t34 - t1088, 0, 0, t987, -t1283, -t1300, t933, qJD(3) * t100 + qJD(4) * t101 - t1239, t825 - t1060, t824 - t1059, qJD(3) * t31 + qJD(4) * t32 + qJD(5) * t196 - t1133, t154, t117, t96, t155, t97, -t279, -qJD(3) * t45 + qJD(4) * t51 + qJD(6) * t68 - t1151, -qJD(3) * t43 + qJD(4) * t49 + qJD(6) * t70 - t1150, -qJD(3) * t13 - qJD(4) * t16 - t1163, -qJD(3) * t1 - qJD(4) * t3 + qJD(5) * t27 - t1172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -pkin(2) * t895, -pkin(2) * t897, 0, 0, 0, 0, 0, 0, 0, 0, -t1244, -t973 + t979, 0, qJD(3) * t377, 0, 0, 0, 0, 0, 0, 0, t1244, -t979 + t968, qJD(3) * t376 + qJD(4) * t365 + t612, -t698, t676, 0, t698, 0, 0, t1015 * t635 + t1064 + t622, -t1016 * t635 + t1063 + t623, -t502 - t468, qJD(3) * t267 + qJD(4) * t257 + t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t984 + (-t896 - t895) * pkin(2) (-t898 - t897) * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t1267, -t973 - t1272, -t1011 (-t1183 * t660 + t661 * t725) * t1174 + t831, 0, 0, 0, 0, 0, 0, t1013, t1267, t968 + t1272 (t660 * t704 + t1118) * qJD(3) + t221 * qJD(4) + t558 * qJD(5) + t833, -t698, t676, 0, t698, 0, 0, qJD(4) * t433 + qJD(6) * t422 + t622 + t718 - t797, qJD(4) * t435 + qJD(6) * t421 + t623 + t719 - t798, qJD(4) * t297 - t502 - t836 (t660 * t894 + t1118) * qJD(3) + t161 * qJD(4) + t381 * qJD(5) + t848; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t988, 0, -t1269, t1268, 0, 0, 0, 0, t988, 0, 0, 0, t1012, t1269, qJD(5) - t1268, t221 * qJD(3) + (-pkin(4) * t646 - t1123) * qJD(4) + t612 + t832, t312, t676, t992, -t312, 0, 0, qJD(3) * t433 + qJD(6) * t432 + t1064 - t795, qJD(3) * t435 + qJD(6) * t431 + t1063 - t796, qJD(3) * t297 - t468 - t835, t161 * qJD(3) + (-t1217 * t477 - t1123) * qJD(4) + t383 * qJD(5) + t847; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, qJD(3) * t558 + qJD(4) * t635 - t1237, 0, 0, 0, 0, 0, 0, t954, t1251, 0, qJD(3) * t381 + qJD(4) * t383 - t1232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t368, t325, t291, t326, -t989, qJD(3) * t422 + qJD(4) * t432 - t1016 * t632 - t794, qJD(3) * t421 + qJD(4) * t431 - t1015 * t632 - t793, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t566, -t991, 0, t565, 0, 0, -t568 * qJD(2) - t764, -t939, 0, 0, -t1283, -t1300, 0, t933, -t987, 0, -t827 + t1047, -t826 + t1046, qJD(2) * t105 - t1238, qJD(2) * t34 - t1087, 0, 0, t987, -t1283, -t1300, t933, -qJD(2) * t100 + qJD(4) * t156 - t1238, t827 - t1058, t826 - t1057, -qJD(2) * t31 + qJD(4) * t40 + qJD(5) * t194 - t1132, t154, t117, t96, t155, t97, -t279, qJD(2) * t45 + qJD(4) * t67 + qJD(6) * t78 - t1149, qJD(2) * t43 + qJD(4) * t65 + qJD(6) * t80 - t1148, qJD(2) * t13 - qJD(4) * t18 - t1162, qJD(2) * t1 - qJD(4) * t5 + qJD(5) * t29 - t1171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pkin(2) * t896 - t984, pkin(2) * t898, 0, 0, 0, 0, 0, 0, 0, 0, -t1273, t1274, t1011, -t831, 0, 0, 0, 0, 0, 0, -t1013, t1273, qJD(5) - t1274, qJD(4) * t220 - qJD(5) * t582 - t833, -t698, t676, 0, t698, 0, 0, qJD(4) * t434 - qJD(6) * t420 + t718 + t797, qJD(4) * t436 + qJD(6) * t419 + t719 + t798, qJD(4) * t296 + t836, -qJD(4) * t160 - qJD(5) * t380 - t848; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t711, -t712, 0, 0, 0, 0, 0, 0, 0, 0, 0, t711, t967, -qJD(4) * t610 + t692, -t698, t676, 0, t698, 0, 0, t1015 * t700 + t1067, -t1016 * t700 + t1066, -t647, qJD(4) * t503 + t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t988, 0, -t1280, -t712 + t959, 0, 0, 0, 0, t988, 0, 0, 0, t1010, t1280, -t959 + t967 (-pkin(4) * t1175 + t713) * qJD(4) + t692 + t783, t312, t676, t992, -t312, 0, 0, t586 * qJD(6) + t1067 - t753, t585 * qJD(6) + t1066 - t752, -t647 - t784 (-t1217 * t659 + t713) * qJD(4) + t597 * qJD(5) + t786; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, qJD(4) * t700 - t1226, 0, 0, 0, 0, 0, 0, t954, t1251, 0, qJD(4) * t597 - t1246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t368, t325, t291, t326, -t989, qJD(4) * t586 - t1016 * t699 - t769, qJD(4) * t585 - t1015 * t699 - t770, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1283, -t1300, 0, t1283, t1248, 0, -t941 - t1077, -t828 - t944, 0, 0, 0, 0, -t1248, -t1283, -t1300, t1283, -qJD(2) * t101 - qJD(3) * t156, -t1056 + t1077, t828 - t1055, -qJD(2) * t32 - qJD(3) * t40 - t1128, t1283 * t722 + t1224, t1258 * t865 + t478, -t360 * t966 + t1305 - t935, t1283 * t723 - t1224, t1306 - t934, t892, -qJD(2) * t51 - qJD(3) * t67 + qJD(6) * t166 - t1147, -qJD(2) * t49 - qJD(3) * t65 + qJD(6) * t167 - t1146, qJD(2) * t16 + qJD(3) * t18 - t1161, qJD(2) * t3 + qJD(3) * t5 + qJD(5) * t37 - t1170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t988, 0, t1270, t1271, 0, 0, 0, 0, -t988, 0, 0, 0, -t1012, -t1270, qJD(5) - t1271, -qJD(3) * t220 + t721 - t832, -t311, t676, -t992, t311, 0, 0, -qJD(3) * t434 - qJD(6) * t430 + t718 + t795, -qJD(3) * t436 + qJD(6) * t429 + t719 + t796, -qJD(3) * t296 + t835, qJD(3) * t160 - qJD(5) * t382 - t847; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t988, 0, t956, -t959, 0, 0, 0, 0, -t988, 0, 0, 0, -t1010, -t956, qJD(5) + t959, t721 - t783, -t311, t676, -t992, t311, 0, 0, -t584 * qJD(6) + t718 + t753, t583 * qJD(6) + t719 + t752, t784, -qJD(5) * t596 - t786; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t721, -t698, t676, 0, t698, 0, 0, qJ(5) * t1015 + t718, -qJ(5) * t1016 + t719, 0, t721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, t671, 0, 0, 0, 0, 0, 0, t954, t1251, 0, t754; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t368, -t830, t291, -t829, t990, t1014 * t724 - t748, t1014 * t727 - t747, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1283, -t986, t945 + t1077, 0, 0, 0, 0, 0, 0, qJD(6) * t351 - t940 (-qJD(6) * t1258 - t986) * t727, -t1002, -qJD(2) * t27 - qJD(3) * t29 - qJD(4) * t37 + t1166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t714, qJD(3) * t582 - t1061 + t1237, 0, 0, 0, 0, 0, 0, -t954, -t1251, 0, qJD(3) * t380 + qJD(4) * t382 + t1232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t714, -t1061 + t1226, 0, 0, 0, 0, 0, 0, -t954, -t1251, 0, qJD(4) * t596 + t1246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t714, -t671, 0, 0, 0, 0, 0, 0, -t954, -t1251, 0, -t754; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t994 - t1016, -t829, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t1266 (t1035 * t553 + t715) * t1258 + t966 * t345, t113 (-t1036 * t553 + t1017) * t1258 + t966 * t354, t258, -qJD(2) * t68 - qJD(3) * t78 - qJD(4) * t166 - qJD(5) * t351 - t1164, -qJD(2) * t70 - qJD(3) * t80 - qJD(4) * t167 + t1258 * t719 - t1165, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, t367, t332, -t291, t339, t989, qJD(3) * t420 + qJD(4) * t430 + t794, -qJD(3) * t419 - qJD(4) * t429 + t793, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, t367, t332, -t291, t339, t989, qJD(4) * t584 + t769, -qJD(4) * t583 + t770, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, t367, t1258 * t1036, -t291, t929, -t990, t748, t747, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t994, t929, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t15;
