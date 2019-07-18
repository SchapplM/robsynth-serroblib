% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 10:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRPPRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:24:32
% EndTime: 2019-05-06 10:25:46
% DurationCPUTime: 77.36s
% Computational Cost: add. (233165->946), mult. (631453->1457), div. (0->0), fcn. (492388->12), ass. (0->691)
t1121 = sin(qJ(1));
t1125 = cos(qJ(1));
t1115 = sin(pkin(6));
t1117 = cos(pkin(6));
t1120 = sin(qJ(2));
t1124 = cos(qJ(2));
t1106 = qJD(1) * t1117 + qJD(2);
t1295 = t1106 ^ 2;
t1114 = sin(pkin(11));
t1116 = cos(pkin(11));
t1265 = qJD(1) * t1115;
t1070 = (t1114 * t1124 + t1116 * t1120) * t1265;
t1296 = t1070 ^ 2;
t1027 = -t1296 - t1295;
t1105 = qJDD(1) * t1117 + qJDD(2);
t1242 = t1115 * t1124;
t1228 = qJD(1) * t1242;
t1243 = t1115 * t1120;
t1229 = qJD(1) * t1243;
t1068 = t1114 * t1229 - t1116 * t1228;
t1254 = t1070 * t1068;
t1351 = t1254 + t1105;
t1364 = t1114 * t1351;
t926 = -t1027 * t1116 + t1364;
t1362 = t1116 * t1351;
t928 = t1027 * t1114 + t1362;
t1177 = t1120 * t928 + t1124 * t926;
t1234 = qJDD(1) * t1120;
t1264 = qJD(1) * t1124;
t1079 = (qJD(2) * t1264 + t1234) * t1115;
t1235 = qJDD(1) * t1115;
t1080 = -qJD(2) * t1229 + t1124 * t1235;
t1020 = t1116 * t1079 + t1114 * t1080;
t1255 = t1068 * t1106;
t1309 = -t1020 + t1255;
t813 = -t1115 * t1309 + t1117 * t1177;
t854 = t1120 * t926 - t1124 * t928;
t733 = -t1121 * t854 + t1125 * t813;
t1430 = pkin(7) * t733;
t736 = t1121 * t813 + t1125 * t854;
t1429 = pkin(7) * t736;
t811 = t1115 * t1177 + t1117 * t1309;
t1428 = pkin(8) * (t1115 * t811 + t1117 * t813);
t1427 = pkin(1) * t811;
t1426 = pkin(1) * t813;
t1297 = t1068 ^ 2;
t1042 = t1297 - t1295;
t939 = t1042 * t1114 + t1362;
t944 = -t1042 * t1116 + t1364;
t1172 = t1120 * t944 - t1124 * t939;
t1019 = t1079 * t1114 - t1116 * t1080;
t1248 = t1106 * t1070;
t1352 = -t1019 + t1248;
t822 = t1115 * t1352 + t1117 * t1172;
t860 = t1120 * t939 + t1124 * t944;
t1425 = t1121 * t822 - t1125 * t860;
t1043 = -t1296 + t1295;
t1000 = t1254 - t1105;
t1373 = t1000 * t1114;
t937 = -t1043 * t1116 + t1373;
t1372 = t1000 * t1116;
t941 = t1043 * t1114 + t1372;
t1175 = t1120 * t941 + t1124 * t937;
t1308 = t1020 + t1255;
t819 = t1115 * t1308 + t1117 * t1175;
t858 = t1120 * t937 - t1124 * t941;
t1424 = t1121 * t819 + t1125 * t858;
t1423 = t1121 * t860 + t1125 * t822;
t1422 = -t1121 * t858 + t1125 * t819;
t995 = -t1295 - t1297;
t912 = t1114 * t995 - t1372;
t915 = -t1116 * t995 - t1373;
t1178 = t1120 * t915 - t1124 * t912;
t964 = t1019 + t1248;
t803 = t1115 * t964 + t1117 * t1178;
t834 = t1120 * t912 + t1124 * t915;
t722 = t1121 * t834 + t1125 * t803;
t1420 = pkin(7) * t722;
t723 = t1121 * t803 - t1125 * t834;
t1419 = pkin(7) * t723;
t1418 = pkin(8) * t854;
t801 = t1115 * t1178 - t1117 * t964;
t1417 = pkin(8) * (t1115 * t801 + t1117 * t803);
t1410 = t1115 * t1172 - t1117 * t1352;
t1409 = t1115 * t1175 - t1117 * t1308;
t1370 = -t1114 * t964 - t1116 * t1309;
t1280 = t1116 * t964;
t1379 = t1114 * t1309 - t1280;
t1390 = -t1120 * t1370 + t1124 * t1379;
t1305 = t1296 - t1297;
t1392 = t1120 * t1379 + t1124 * t1370;
t1399 = -t1115 * t1305 + t1117 * t1392;
t1408 = -t1121 * t1399 + t1125 * t1390;
t1407 = t1121 * t1390 + t1125 * t1399;
t1406 = pkin(1) * t801;
t1405 = pkin(1) * t803;
t1359 = t1114 * t1308 + t1116 * t1352;
t1360 = t1114 * t1352 - t1116 * t1308;
t1369 = -t1120 * t1360 + t1124 * t1359;
t1306 = -t1296 - t1297;
t1368 = t1120 * t1359 + t1124 * t1360;
t1377 = -t1115 * t1306 + t1117 * t1368;
t1389 = t1121 * t1369 + t1125 * t1377;
t1404 = pkin(7) * t1389;
t1391 = -t1121 * t1377 + t1125 * t1369;
t1403 = pkin(7) * t1391;
t1400 = t1115 * t1392 + t1117 * t1305;
t1397 = pkin(2) * t926;
t1396 = pkin(1) * t1377;
t1378 = t1115 * t1368 + t1117 * t1306;
t1395 = pkin(1) * t1378;
t1394 = qJ(3) * t926;
t1393 = qJ(3) * t928;
t1388 = (-t1115 * t1378 - t1117 * t1377) * pkin(8);
t1387 = pkin(8) * t834;
t1386 = pkin(8) * t1369;
t1376 = qJ(3) * t1360;
t1371 = -pkin(2) * t1306 + qJ(3) * t1359;
t1367 = pkin(2) * t912;
t1366 = qJ(3) * t912;
t1365 = qJ(3) * t915;
t1153 = (-t1068 * t1114 - t1070 * t1116) * t1106;
t1154 = (-t1068 * t1116 + t1070 * t1114) * t1106;
t1303 = -t1120 * t1153 + t1124 * t1154;
t1244 = t1115 * t1105;
t1304 = t1120 * t1154 + t1124 * t1153;
t1318 = t1117 * t1304 - t1244;
t1350 = -t1121 * t1318 + t1125 * t1303;
t1246 = t1106 * t1116;
t1155 = t1019 * t1114 + t1068 * t1246;
t1247 = t1106 * t1114;
t1210 = -t1116 * t1019 + t1068 * t1247;
t1302 = -t1120 * t1210 + t1124 * t1155;
t1222 = t1115 * t1254;
t1301 = t1120 * t1155 + t1124 * t1210;
t1319 = t1117 * t1301 + t1222;
t1349 = -t1121 * t1319 + t1125 * t1302;
t1209 = t1116 * t1020 - t1070 * t1247;
t1211 = t1114 * t1020 + t1070 * t1246;
t1300 = -t1120 * t1211 + t1124 * t1209;
t1299 = t1120 * t1209 + t1124 * t1211;
t1320 = t1117 * t1299 - t1222;
t1348 = -t1121 * t1320 + t1125 * t1300;
t1347 = t1121 * t1303 + t1125 * t1318;
t1346 = t1121 * t1302 + t1125 * t1319;
t1345 = t1121 * t1300 + t1125 * t1320;
t1344 = 2 * qJD(4);
t1342 = qJ(4) * t1309;
t1118 = sin(qJ(6));
t1119 = sin(qJ(5));
t1123 = cos(qJ(5));
t1034 = t1068 * t1119 + t1106 * t1123;
t1216 = -t1123 * t1019 + t1119 * t1105;
t935 = -qJD(5) * t1034 - t1216;
t1134 = qJDD(6) - t935;
t1062 = qJD(5) + t1070;
t1122 = cos(qJ(6));
t982 = t1034 * t1118 - t1122 * t1062;
t984 = t1034 * t1122 + t1062 * t1118;
t903 = t984 * t982;
t1317 = t1134 - t903;
t1333 = t1118 * t1317;
t1014 = qJDD(5) + t1020;
t1032 = -t1123 * t1068 + t1106 * t1119;
t963 = t1034 * t1032;
t1316 = -t963 + t1014;
t1332 = t1119 * t1316;
t1328 = t1122 * t1317;
t1327 = t1123 * t1316;
t1091 = t1117 * t1105;
t1323 = t1115 * t1304 + t1091;
t1221 = t1117 * t1254;
t1322 = t1115 * t1299 + t1221;
t1321 = t1115 * t1301 - t1221;
t936 = -t1032 * qJD(5) + t1119 * t1019 + t1123 * t1105;
t847 = -t982 * qJD(6) + t1118 * t1014 + t1122 * t936;
t1028 = qJD(6) + t1032;
t924 = t1028 * t982;
t790 = -t924 + t847;
t1258 = t1032 * t1062;
t1156 = t936 - t1258;
t1111 = t1115 ^ 2;
t1126 = qJD(1) ^ 2;
t1266 = qJD(1) * t1106;
t1315 = t1111 * (-t1117 * t1126 + t1266);
t1088 = t1106 * t1229;
t1041 = t1080 - t1088;
t1311 = t1117 * t1041;
t1089 = t1106 * t1228;
t1040 = t1089 + t1079;
t1100 = g(1) * t1121 - t1125 * g(2);
t1289 = pkin(8) * t1115;
t1075 = qJDD(1) * pkin(1) + t1126 * t1289 + t1100;
t1049 = t1117 * g(3) + t1115 * t1075;
t1152 = pkin(2) * t1106 - qJ(3) * t1229;
t1113 = t1124 ^ 2;
t1245 = t1111 * t1126;
t1224 = t1113 * t1245;
t952 = t1080 * pkin(2) + qJ(3) * t1224 - t1152 * t1229 - qJDD(3) + t1049;
t1307 = -pkin(3) * t1248 + t1070 * t1344 + t952;
t1101 = g(1) * t1125 + g(2) * t1121;
t1076 = -pkin(1) * t1126 + pkin(8) * t1235 - t1101;
t1253 = t1075 * t1117;
t1214 = -t1120 * t1076 + t1124 * t1253;
t1240 = t1120 * t1126;
t1223 = t1111 * t1240;
t920 = t1105 * pkin(2) - t1079 * qJ(3) + (pkin(2) * t1223 + (qJ(3) * t1266 - g(3)) * t1115) * t1124 + t1214;
t994 = -g(3) * t1243 + t1124 * t1076 + t1120 * t1253;
t921 = -pkin(2) * t1224 + t1080 * qJ(3) - t1106 * t1152 + t994;
t827 = -0.2e1 * qJD(3) * t1068 + t1114 * t920 + t1116 * t921;
t1217 = -t1122 * t1014 + t1118 * t936;
t787 = (qJD(6) - t1028) * t984 + t1217;
t980 = t982 ^ 2;
t981 = t984 ^ 2;
t1025 = t1028 ^ 2;
t1029 = t1032 ^ 2;
t1030 = t1034 ^ 2;
t1298 = t1062 ^ 2;
t1294 = pkin(3) + pkin(9);
t1293 = pkin(2) * t1115;
t1292 = pkin(2) * t1117;
t1291 = pkin(3) * t1114;
t1290 = pkin(3) * t1116;
t1288 = t1019 * pkin(3);
t1287 = t1105 * pkin(3);
t1218 = t1114 * t921 - t1116 * t920;
t1212 = qJDD(4) + t1218;
t1135 = -qJ(4) * t1295 + t1212;
t1003 = pkin(3) * t1068 - qJ(4) * t1070;
t1238 = 0.2e1 * qJD(3) + t1003;
t1128 = -t1294 * t1105 + t1308 * pkin(4) + (pkin(9) * t1068 + t1238) * t1070 + t1135;
t1036 = pkin(4) * t1070 - pkin(9) * t1106;
t1129 = -t1307 + t1342;
t763 = -pkin(4) * t1297 + t1019 * t1294 - t1036 * t1070 + t1129;
t679 = t1119 * t1128 + t1123 * t763;
t961 = pkin(5) * t1032 - pkin(10) * t1034;
t660 = -pkin(5) * t1298 + pkin(10) * t1014 - t1032 * t961 + t679;
t775 = -pkin(3) * t1295 + t1105 * qJ(4) - t1068 * t1003 + t1106 * t1344 + t827;
t745 = -t1019 * pkin(4) - pkin(9) * t1297 + t1106 * t1036 + t775;
t690 = -t1156 * pkin(10) + (t1034 * t1062 - t935) * pkin(5) + t745;
t615 = t1118 * t690 + t1122 * t660;
t1285 = t1114 * t952;
t1281 = t1116 * t952;
t678 = t1119 * t763 - t1123 * t1128;
t659 = -t1014 * pkin(5) - pkin(10) * t1298 + t1034 * t961 + t678;
t1276 = t1118 * t659;
t836 = t1134 + t903;
t1275 = t1118 * t836;
t1274 = t1119 * t745;
t909 = t963 + t1014;
t1273 = t1119 * t909;
t1262 = qJD(3) * t1070;
t826 = t1218 + 0.2e1 * t1262;
t730 = t1114 * t827 - t1116 * t826;
t1272 = t1120 * t730;
t1271 = t1122 * t659;
t1270 = t1122 * t836;
t1269 = t1123 * t745;
t1268 = t1123 * t909;
t1267 = t1124 * t730;
t1260 = t1028 * t1118;
t1259 = t1028 * t1122;
t1257 = t1062 * t1119;
t1256 = t1062 * t1123;
t1096 = t1124 * t1223;
t1077 = t1096 + t1105;
t1252 = t1077 * t1120;
t1251 = t1077 * t1124;
t1078 = -t1096 + t1105;
t1250 = t1078 * t1120;
t1249 = t1078 * t1124;
t1241 = t1120 * t1049;
t1239 = t1124 * t1049;
t1112 = t1120 ^ 2;
t1236 = t1112 + t1113;
t1233 = t1119 * t903;
t1232 = t1123 * t903;
t1231 = pkin(5) * t1123 + pkin(4);
t1230 = t1106 * t1265;
t1227 = t1114 * t963;
t1226 = t1116 * t963;
t1225 = t1112 * t1245;
t1220 = -qJ(4) * t1114 - pkin(2);
t1219 = pkin(5) * t1119 + qJ(4);
t731 = t1114 * t826 + t1116 * t827;
t614 = t1118 * t660 - t1122 * t690;
t571 = t1118 * t614 + t1122 * t615;
t1053 = -t1100 * t1121 - t1125 * t1101;
t1095 = qJDD(1) * t1125 - t1121 * t1126;
t1213 = -pkin(7) * t1095 - g(3) * t1121;
t1063 = -t1295 - t1225;
t1002 = -t1063 * t1120 - t1249;
t1207 = pkin(8) * t1002 - t1241;
t1084 = -t1295 - t1224;
t1023 = t1084 * t1124 - t1252;
t1206 = pkin(8) * t1023 + t1239;
t570 = t1118 * t615 - t1122 * t614;
t621 = t1119 * t679 - t1123 * t678;
t622 = t1119 * t678 + t1123 * t679;
t560 = t1119 * t571 - t1123 * t659;
t539 = t1114 * t570 - t1116 * t560;
t540 = t1114 * t560 + t1116 * t570;
t1205 = t1120 * t540 + t1124 * t539;
t602 = t1114 * t745 - t1116 * t621;
t603 = t1114 * t621 + t1116 * t745;
t1204 = t1120 * t603 + t1124 * t602;
t791 = -t924 - t847;
t715 = -t1118 * t791 - t1122 * t787;
t855 = t980 + t981;
t676 = t1119 * t715 + t1123 * t855;
t713 = -t1118 * t787 + t1122 * t791;
t628 = t1114 * t713 - t1116 * t676;
t629 = t1114 * t676 + t1116 * t713;
t1203 = t1120 * t629 + t1124 * t628;
t788 = (-qJD(6) - t1028) * t984 - t1217;
t714 = -t1118 * t790 + t1122 * t788;
t902 = -t981 + t980;
t685 = -t1119 * t714 - t1123 * t902;
t712 = t1118 * t788 + t1122 * t790;
t633 = t1114 * t712 + t1116 * t685;
t634 = -t1114 * t685 + t1116 * t712;
t1202 = t1120 * t634 + t1124 * t633;
t870 = -t1025 - t980;
t750 = t1122 * t870 - t1333;
t692 = t1119 * t750 + t1123 * t788;
t749 = t1118 * t870 + t1328;
t647 = t1114 * t749 - t1116 * t692;
t648 = t1114 * t692 + t1116 * t749;
t1201 = t1120 * t648 + t1124 * t647;
t899 = -t981 - t1025;
t762 = -t1118 * t899 - t1270;
t695 = t1119 * t762 - t1123 * t790;
t761 = t1122 * t899 - t1275;
t653 = t1114 * t761 - t1116 * t695;
t654 = t1114 * t695 + t1116 * t761;
t1200 = t1120 * t654 + t1124 * t653;
t923 = -t981 + t1025;
t778 = -t1118 * t923 + t1328;
t705 = -t1119 * t778 - t1123 * t791;
t776 = t1122 * t923 + t1333;
t663 = t1114 * t776 + t1116 * t705;
t665 = -t1114 * t705 + t1116 * t776;
t1199 = t1120 * t665 + t1124 * t663;
t922 = t980 - t1025;
t779 = t1122 * t922 - t1275;
t706 = -t1119 * t779 - t1123 * t787;
t777 = t1118 * t922 + t1270;
t664 = t1114 * t777 + t1116 * t706;
t666 = -t1114 * t706 + t1116 * t777;
t1198 = t1120 * t666 + t1124 * t664;
t846 = -qJD(6) * t984 - t1217;
t784 = -t1118 * t846 + t1259 * t982;
t741 = -t1119 * t784 - t1232;
t783 = t1122 * t846 + t1260 * t982;
t681 = t1114 * t783 + t1116 * t741;
t683 = -t1114 * t741 + t1116 * t783;
t1197 = t1120 * t683 + t1124 * t681;
t786 = t1122 * t847 - t1260 * t984;
t742 = -t1119 * t786 + t1232;
t785 = t1118 * t847 + t1259 * t984;
t682 = t1114 * t785 + t1116 * t742;
t684 = -t1114 * t742 + t1116 * t785;
t1196 = t1120 * t684 + t1124 * t682;
t1059 = -0.2e1 * t1262;
t780 = -t1070 * t1003 + t1059 - t1135 + t1287;
t703 = t1114 * t775 + t1116 * t780;
t704 = -t1114 * t780 + t1116 * t775;
t1195 = t1120 * t704 + t1124 * t703;
t849 = (t1118 * t984 - t1122 * t982) * t1028;
t793 = -t1119 * t849 + t1123 * t1134;
t848 = (-t1118 * t982 - t1122 * t984) * t1028;
t725 = t1114 * t848 + t1116 * t793;
t726 = -t1114 * t793 + t1116 * t848;
t1194 = t1120 * t726 + t1124 * t725;
t1193 = t1120 * t731 + t1267;
t1133 = (-qJD(5) + t1062) * t1034 - t1216;
t883 = -t936 - t1258;
t796 = t1119 * t1133 + t1123 * t883;
t919 = -t1029 - t1030;
t746 = t1114 * t919 - t1116 * t796;
t747 = t1114 * t796 + t1116 * t919;
t1192 = t1120 * t747 + t1124 * t746;
t879 = (qJD(5) + t1062) * t1034 + t1216;
t795 = t1119 * t879 - t1123 * t1156;
t962 = t1030 - t1029;
t752 = t1114 * t962 + t1116 * t795;
t753 = -t1114 * t795 + t1116 * t962;
t1191 = t1120 * t753 + t1124 * t752;
t930 = -t1298 - t1029;
t838 = t1119 * t930 + t1327;
t758 = t1114 * t879 - t1116 * t838;
t759 = t1114 * t838 + t1116 * t879;
t1190 = t1120 * t759 + t1124 * t758;
t947 = -t1030 - t1298;
t850 = t1123 * t947 - t1273;
t764 = t1114 * t1156 - t1116 * t850;
t765 = t1114 * t850 + t1116 * t1156;
t1189 = t1120 * t765 + t1124 * t764;
t990 = -t1030 + t1298;
t861 = -t1123 * t990 - t1332;
t770 = -t1114 * t883 + t1116 * t861;
t772 = -t1114 * t861 - t1116 * t883;
t1188 = t1120 * t772 + t1124 * t770;
t989 = t1029 - t1298;
t862 = -t1119 * t989 - t1268;
t771 = t1114 * t1133 + t1116 * t862;
t773 = -t1114 * t862 + t1116 * t1133;
t1187 = t1120 * t773 + t1124 * t771;
t875 = -t1032 * t1257 - t1123 * t935;
t828 = t1116 * t875 - t1227;
t830 = -t1114 * t875 - t1226;
t1186 = t1120 * t830 + t1124 * t828;
t877 = -t1034 * t1256 - t1119 * t936;
t829 = t1116 * t877 + t1227;
t831 = -t1114 * t877 + t1226;
t1185 = t1120 * t831 + t1124 * t829;
t906 = (t1032 * t1119 + t1034 * t1123) * t1062;
t866 = t1014 * t1114 + t1116 * t906;
t867 = t1014 * t1116 - t1114 * t906;
t1184 = t1120 * t867 + t1124 * t866;
t993 = g(3) * t1242 - t1214;
t1165 = t1120 * t994 - t1124 * t993;
t911 = t1120 * t993 + t1124 * t994;
t1037 = t1080 + t1088;
t1038 = -t1089 + t1079;
t1163 = t1037 * t1120 - t1038 * t1124;
t1162 = t1040 * t1124 + t1041 * t1120;
t1161 = t1063 * t1124 - t1250;
t1082 = t1295 - t1225;
t1160 = t1082 * t1124 + t1252;
t1159 = t1084 * t1120 + t1251;
t1083 = -t1295 + t1224;
t1158 = t1083 * t1120 + t1249;
t1052 = t1100 * t1125 - t1101 * t1121;
t561 = t1119 * t659 + t1123 * t571;
t521 = -t1294 * t561 + (pkin(10) * t1119 + t1231) * t570;
t525 = pkin(4) * t560 - pkin(5) * t659 + pkin(10) * t571 - qJ(4) * t561;
t506 = -pkin(2) * t561 + qJ(3) * t540 + t1114 * t525 + t1116 * t521;
t508 = -qJ(3) * t539 - t1114 * t521 + t1116 * t525;
t522 = -t1120 * t539 + t1124 * t540;
t1151 = pkin(8) * t522 + t1120 * t508 + t1124 * t506;
t563 = -pkin(10) * t713 - t570;
t677 = -t1119 * t855 + t1123 * t715;
t543 = -t1119 * t563 + t1231 * t713 - t1294 * t677;
t549 = pkin(4) * t676 + pkin(5) * t855 + pkin(10) * t715 - qJ(4) * t677 + t571;
t523 = -pkin(2) * t677 + qJ(3) * t629 + t1114 * t549 + t1116 * t543;
t524 = -qJ(3) * t628 - t1114 * t543 + t1116 * t549;
t582 = -t1120 * t628 + t1124 * t629;
t1150 = pkin(8) * t582 + t1120 * t524 + t1124 * t523;
t566 = pkin(4) * t745 - t1294 * t622;
t573 = pkin(4) * t621 - qJ(4) * t622;
t527 = -pkin(2) * t622 + qJ(3) * t603 + t1114 * t573 + t1116 * t566;
t535 = -qJ(3) * t602 - t1114 * t566 + t1116 * t573;
t562 = -t1120 * t602 + t1124 * t603;
t1149 = pkin(8) * t562 + t1120 * t535 + t1124 * t527;
t598 = -pkin(5) * t749 + t614;
t631 = -pkin(10) * t749 + t1276;
t693 = -t1119 * t788 + t1123 * t750;
t554 = pkin(4) * t749 - t1119 * t631 - t1123 * t598 - t1294 * t693;
t575 = pkin(4) * t692 + pkin(5) * t788 + pkin(10) * t750 - qJ(4) * t693 - t1271;
t528 = -pkin(2) * t693 + qJ(3) * t648 + t1114 * t575 + t1116 * t554;
t532 = -qJ(3) * t647 - t1114 * t554 + t1116 * t575;
t595 = -t1120 * t647 + t1124 * t648;
t1148 = pkin(8) * t595 + t1120 * t532 + t1124 * t528;
t600 = -pkin(5) * t761 + t615;
t635 = -pkin(10) * t761 + t1271;
t696 = t1119 * t790 + t1123 * t762;
t555 = pkin(4) * t761 - t1119 * t635 - t1123 * t600 - t1294 * t696;
t576 = pkin(4) * t695 - pkin(5) * t790 + pkin(10) * t762 - qJ(4) * t696 + t1276;
t531 = -pkin(2) * t696 + qJ(3) * t654 + t1114 * t576 + t1116 * t555;
t534 = -qJ(3) * t653 - t1114 * t555 + t1116 * t576;
t599 = -t1120 * t653 + t1124 * t654;
t1147 = pkin(8) * t599 + t1120 * t534 + t1124 * t531;
t798 = -t1119 * t883 + t1123 * t1133;
t594 = pkin(4) * t919 - t1294 * t798 - t622;
t716 = pkin(4) * t796 - qJ(4) * t798;
t567 = -pkin(2) * t798 + qJ(3) * t747 + t1114 * t716 + t1116 * t594;
t572 = -qJ(3) * t746 - t1114 * t594 + t1116 * t716;
t675 = -t1120 * t746 + t1124 * t747;
t1146 = pkin(8) * t675 + t1120 * t572 + t1124 * t567;
t839 = t1123 * t930 - t1332;
t642 = pkin(4) * t838 - qJ(4) * t839 - t678;
t656 = pkin(4) * t879 - t1294 * t839 + t1269;
t583 = -pkin(2) * t839 + qJ(3) * t759 + t1114 * t642 + t1116 * t656;
t590 = -qJ(3) * t758 - t1114 * t656 + t1116 * t642;
t687 = -t1120 * t758 + t1124 * t759;
t1145 = pkin(8) * t687 + t1120 * t590 + t1124 * t583;
t851 = -t1119 * t947 - t1268;
t643 = pkin(4) * t850 - qJ(4) * t851 - t679;
t661 = pkin(4) * t1156 - t1294 * t851 - t1274;
t588 = -pkin(2) * t851 + qJ(3) * t765 + t1114 * t643 + t1116 * t661;
t593 = -qJ(3) * t764 - t1114 * t661 + t1116 * t643;
t691 = -t1120 * t764 + t1124 * t765;
t1144 = pkin(8) * t691 + t1120 * t593 + t1124 * t588;
t824 = t1129 + t1288;
t630 = qJ(3) * t704 + (t1220 - t1290) * t824;
t639 = -t1120 * t703 + t1124 * t704;
t644 = -qJ(3) * t703 + (-qJ(4) * t1116 + t1291) * t824;
t1143 = pkin(8) * t639 + t1120 * t644 + t1124 * t630;
t748 = -pkin(3) * t1306 + t775;
t1131 = t1070 * t1238 + t1212;
t751 = -t1287 + (-t1295 - t1306) * qJ(4) + t1131;
t657 = t1114 * t751 + t1116 * t748 + t1371;
t670 = -t1114 * t748 + t1116 * t751 - t1376;
t1142 = t1120 * t670 + t1124 * t657 + t1386;
t694 = t1371 + t731;
t711 = -t730 - t1376;
t1141 = t1120 * t711 + t1124 * t694 + t1386;
t782 = (t1019 + t964) * pkin(3) + t1129;
t717 = t1116 * t782 - t1220 * t964 + t1365;
t728 = qJ(4) * t1280 - t1114 * t782 + t1366;
t1140 = t1120 * t728 + t1124 * t717 + t1387;
t781 = -t1288 + t1307 - 0.2e1 * t1342;
t719 = t1393 + t1114 * t781 + (-pkin(2) - t1290) * t1309;
t729 = t1116 * t781 + t1291 * t1309 - t1394;
t1139 = t1120 * t729 + t1124 * t719 - t1418;
t815 = -pkin(2) * t964 + t1281 - t1365;
t856 = -t1285 - t1366;
t1138 = t1120 * t856 + t1124 * t815 - t1387;
t823 = pkin(2) * t1309 - t1285 - t1393;
t865 = -t1281 + t1394;
t1137 = t1120 * t865 + t1124 * t823 + t1418;
t974 = t1037 * t1124 + t1038 * t1120;
t1136 = pkin(8) * t974 + t911;
t667 = t1124 * t731 - t1272;
t718 = pkin(2) * t952 + qJ(3) * t731;
t1132 = pkin(8) * t667 - qJ(3) * t1272 + t1124 * t718;
t1110 = t1115 * t1111;
t1094 = qJDD(1) * t1121 + t1125 * t1126;
t1086 = t1236 * t1245;
t1085 = (t1112 - t1113) * t1245;
t1081 = -pkin(7) * t1094 + g(3) * t1125;
t1047 = t1236 * t1230;
t1039 = (t1234 + (qJD(2) + t1106) * t1264) * t1115;
t1035 = t1079 * t1124 - t1112 * t1230;
t1031 = -t1080 * t1120 - t1113 * t1230;
t1022 = t1083 * t1124 - t1250;
t1021 = -t1082 * t1120 + t1251;
t992 = (t1110 * t1124 * t1126 + t1040 * t1117) * t1120;
t991 = (-t1110 * t1240 + t1311) * t1124;
t975 = -t1040 * t1120 + t1041 * t1124;
t951 = t1115 * t1041 + t1117 * t1159;
t950 = -t1115 * t1037 + t1117 * t1158;
t949 = -t1115 * t1038 + t1117 * t1160;
t948 = t1115 * t1159 - t1311;
t946 = -t1115 * t1039 + t1117 * t1161;
t945 = t1117 * t1039 + t1115 * t1161;
t934 = -t1115 * t1085 + t1117 * t1162;
t933 = t1115 * t1086 + t1117 * t1163;
t932 = -t1117 * t1086 + t1115 * t1163;
t907 = (-t1032 * t1123 + t1034 * t1119) * t1062;
t905 = t1023 * t1125 - t1121 * t951;
t904 = t1023 * t1121 + t1125 * t951;
t898 = t1002 * t1125 - t1121 * t946;
t897 = t1002 * t1121 + t1125 * t946;
t896 = t1115 * t1049 + t1117 * t1165;
t895 = -t1117 * t1049 + t1115 * t1165;
t878 = -t1034 * t1257 + t1123 * t936;
t876 = t1032 * t1256 - t1119 * t935;
t869 = -t1121 * t933 + t1125 * t974;
t868 = t1121 * t974 + t1125 * t933;
t864 = t1123 * t989 - t1273;
t863 = -t1119 * t990 + t1327;
t852 = -t1241 + (-t1115 * t948 - t1117 * t951) * pkin(8);
t845 = -t1239 + (-t1115 * t945 - t1117 * t946) * pkin(8);
t840 = -pkin(1) * t948 + t1115 * t993 + t1117 * t1206;
t832 = -pkin(1) * t945 + t1115 * t994 + t1117 * t1207;
t818 = pkin(8) * t1117 * t911 - pkin(1) * t895;
t817 = -t1121 * t896 + t1125 * t911;
t816 = t1121 * t911 + t1125 * t896;
t810 = -pkin(1) * t932 + t1117 * t1136;
t809 = pkin(2) * t1360 - pkin(3) * t1308 + qJ(4) * t1352;
t808 = (-t1115 * t895 - t1117 * t896) * pkin(8);
t799 = (-t1115 * t932 - t1117 * t933) * pkin(8) - t1165;
t797 = -t1119 * t1156 - t1123 * t879;
t794 = t1119 * t1134 + t1123 * t849;
t774 = -t827 - t1397;
t769 = t1059 - t1218 + t1367;
t768 = -t1120 * t866 + t1124 * t867;
t744 = t1123 * t786 + t1233;
t743 = t1123 * t784 - t1233;
t738 = -t1120 * t829 + t1124 * t831;
t737 = -t1120 * t828 + t1124 * t830;
t732 = -t1115 * t907 + t1117 * t1184;
t727 = -pkin(3) * t1027 + qJ(4) * t1351 + t1397 + t775;
t720 = -t1367 + (-t1295 - t995) * qJ(4) + (-t1105 + t1000) * pkin(3) + t1131;
t710 = -t1115 * t878 + t1117 * t1185;
t709 = -t1115 * t876 + t1117 * t1186;
t708 = -t1119 * t787 + t1123 * t779;
t707 = -t1119 * t791 + t1123 * t778;
t702 = -t1120 * t771 + t1124 * t773;
t701 = -t1120 * t770 + t1124 * t772;
t686 = -t1119 * t902 + t1123 * t714;
t680 = -t1120 * t752 + t1124 * t753;
t674 = -t1115 * t864 + t1117 * t1187;
t673 = -t1115 * t863 + t1117 * t1188;
t672 = -t1115 * t851 + t1117 * t1189;
t671 = t1115 * t1189 + t1117 * t851;
t669 = -t1115 * t839 + t1117 * t1190;
t668 = t1115 * t1190 + t1117 * t839;
t662 = -t1120 * t823 + t1124 * t865 + t1428;
t655 = -t1120 * t725 + t1124 * t726;
t652 = t1115 * t952 + t1117 * t1193;
t651 = t1115 * t1193 - t1117 * t952;
t650 = -t1120 * t815 + t1124 * t856 + t1417;
t649 = -t1115 * t797 + t1117 * t1191;
t646 = -t1115 * t798 + t1117 * t1192;
t645 = t1115 * t1192 + t1117 * t798;
t641 = -t1115 * t774 + t1117 * t1137 + t1427;
t640 = pkin(2) * t703 + pkin(3) * t780 + qJ(4) * t775;
t638 = -t1115 * t769 + t1117 * t1138 + t1406;
t637 = pkin(2) * t764 + qJ(4) * t1156 - t1294 * t850 + t1269;
t636 = -t1115 * t794 + t1117 * t1194;
t632 = pkin(2) * t758 + qJ(4) * t879 - t1294 * t838 + t1274;
t627 = -t1120 * t682 + t1124 * t684;
t626 = -t1120 * t681 + t1124 * t683;
t625 = -t1115 * t824 + t1117 * t1195;
t624 = t1115 * t1195 + t1117 * t824;
t623 = -t1120 * t719 + t1124 * t729 - t1428;
t620 = -t1121 * t672 + t1125 * t691;
t619 = t1121 * t691 + t1125 * t672;
t618 = -t1120 * t717 + t1124 * t728 - t1417;
t617 = -t1121 * t669 + t1125 * t687;
t616 = t1121 * t687 + t1125 * t669;
t612 = -t1120 * t664 + t1124 * t666;
t611 = -t1120 * t663 + t1124 * t665;
t610 = -t1115 * t744 + t1117 * t1196;
t609 = -t1115 * t743 + t1117 * t1197;
t608 = -t1121 * t646 + t1125 * t675;
t607 = t1121 * t675 + t1125 * t646;
t606 = -t1121 * t652 + t1125 * t667;
t605 = t1121 * t667 + t1125 * t652;
t604 = -t1115 * t727 + t1117 * t1139 - t1427;
t601 = -t1120 * t694 + t1124 * t711 + t1388;
t597 = -t1115 * t720 + t1117 * t1140 - t1406;
t596 = t1117 * t1141 - t1293 * t1360 - t1395;
t592 = -t1115 * t708 + t1117 * t1198;
t591 = -t1115 * t707 + t1117 * t1199;
t589 = -t1120 * t657 + t1124 * t670 + t1388;
t587 = -t1120 * t633 + t1124 * t634;
t586 = pkin(2) * t746 + qJ(4) * t919 - t1294 * t796 - t621;
t585 = -t1115 * t696 + t1117 * t1200;
t584 = t1115 * t1200 + t1117 * t696;
t581 = -t1121 * t625 + t1125 * t639;
t580 = t1121 * t639 + t1125 * t625;
t579 = -t1115 * t809 + t1117 * t1142 - t1395;
t578 = -t1115 * t693 + t1117 * t1201;
t577 = t1115 * t1201 + t1117 * t693;
t574 = -qJ(3) * t1267 - t1120 * t718 + (-t1115 * t651 - t1117 * t652) * pkin(8);
t569 = -t1115 * t686 + t1117 * t1202;
t568 = -pkin(1) * t651 + t1117 * t1132 - t1293 * t730;
t565 = -t1115 * t677 + t1117 * t1203;
t564 = t1115 * t1203 + t1117 * t677;
t559 = -t1121 * t585 + t1125 * t599;
t558 = t1121 * t599 + t1125 * t585;
t557 = -t1121 * t578 + t1125 * t595;
t556 = t1121 * t595 + t1125 * t578;
t553 = -t1120 * t630 + t1124 * t644 + (-t1115 * t624 - t1117 * t625) * pkin(8);
t552 = pkin(2) * t602 + qJ(4) * t745 - t1294 * t621;
t551 = -t1115 * t622 + t1117 * t1204;
t550 = t1115 * t1204 + t1117 * t622;
t548 = pkin(2) * t653 + qJ(4) * t761 - t1119 * t600 + t1123 * t635 - t1294 * t695;
t547 = pkin(2) * t647 + qJ(4) * t749 - t1119 * t598 + t1123 * t631 - t1294 * t692;
t546 = -t1121 * t565 + t1125 * t582;
t545 = t1121 * t582 + t1125 * t565;
t544 = -pkin(1) * t624 - t1115 * t640 + t1117 * t1143;
t542 = -t1120 * t588 + t1124 * t593 + (-t1115 * t671 - t1117 * t672) * pkin(8);
t541 = -t1120 * t583 + t1124 * t590 + (-t1115 * t668 - t1117 * t669) * pkin(8);
t538 = -pkin(1) * t671 - t1115 * t637 + t1117 * t1144;
t537 = pkin(2) * t628 + t1123 * t563 + t1219 * t713 - t1294 * t676;
t536 = -pkin(1) * t668 - t1115 * t632 + t1117 * t1145;
t533 = -t1120 * t567 + t1124 * t572 + (-t1115 * t645 - t1117 * t646) * pkin(8);
t530 = -t1121 * t551 + t1125 * t562;
t529 = t1121 * t562 + t1125 * t551;
t526 = -pkin(1) * t645 - t1115 * t586 + t1117 * t1146;
t520 = -t1115 * t561 + t1117 * t1205;
t519 = t1115 * t1205 + t1117 * t561;
t518 = -t1120 * t531 + t1124 * t534 + (-t1115 * t584 - t1117 * t585) * pkin(8);
t517 = -t1120 * t528 + t1124 * t532 + (-t1115 * t577 - t1117 * t578) * pkin(8);
t516 = pkin(2) * t539 - t1294 * t560 + (-pkin(10) * t1123 + t1219) * t570;
t515 = -pkin(1) * t584 - t1115 * t548 + t1117 * t1147;
t514 = -t1120 * t527 + t1124 * t535 + (-t1115 * t550 - t1117 * t551) * pkin(8);
t513 = -pkin(1) * t577 - t1115 * t547 + t1117 * t1148;
t512 = -pkin(1) * t550 - t1115 * t552 + t1117 * t1149;
t511 = -t1121 * t520 + t1125 * t522;
t510 = t1121 * t522 + t1125 * t520;
t509 = -t1120 * t523 + t1124 * t524 + (-t1115 * t564 - t1117 * t565) * pkin(8);
t507 = -pkin(1) * t564 - t1115 * t537 + t1117 * t1150;
t505 = -t1120 * t506 + t1124 * t508 + (-t1115 * t519 - t1117 * t520) * pkin(8);
t504 = -pkin(1) * t519 - t1115 * t516 + t1117 * t1151;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1094, -t1095, 0, t1053, 0, 0, 0, 0, 0, 0, t905, t898, t869, t817, 0, 0, 0, 0, 0, 0, t723, t736, t1391, t606, 0, 0, 0, 0, 0, 0, t1391, -t723, -t736, t581, 0, 0, 0, 0, 0, 0, t617, t620, t608, t530, 0, 0, 0, 0, 0, 0, t557, t559, t546, t511; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1095, -t1094, 0, t1052, 0, 0, 0, 0, 0, 0, t904, t897, t868, t816, 0, 0, 0, 0, 0, 0, -t722, -t733, t1389, t605, 0, 0, 0, 0, 0, 0, t1389, t722, t733, t580, 0, 0, 0, 0, 0, 0, t616, t619, t607, t529, 0, 0, 0, 0, 0, 0, t556, t558, t545, t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t948, t945, t932, t895, 0, 0, 0, 0, 0, 0, -t801, -t811, t1378, t651, 0, 0, 0, 0, 0, 0, t1378, t801, t811, t624, 0, 0, 0, 0, 0, 0, t668, t671, t645, t550, 0, 0, 0, 0, 0, 0, t577, t584, t564, t519; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1095, 0, -t1094, 0, t1213, -t1081, -t1052, -pkin(7) * t1052, t1035 * t1125 - t1121 * t992, -t1121 * t934 + t1125 * t975, t1021 * t1125 - t1121 * t949, t1031 * t1125 - t1121 * t991, t1022 * t1125 - t1121 * t950, t1047 * t1125 + t1121 * t1244, -pkin(7) * t904 - t1121 * t840 + t1125 * t852, -pkin(7) * t897 - t1121 * t832 + t1125 * t845, -pkin(7) * t868 - t1121 * t810 + t1125 * t799, -pkin(7) * t816 - t1121 * t818 + t1125 * t808, t1348, t1408, t1424, t1349, t1425, t1350, -t1121 * t638 + t1125 * t650 + t1420, -t1121 * t641 + t1125 * t662 + t1430, -t1121 * t596 + t1125 * t601 - t1404, -pkin(7) * t605 - t1121 * t568 + t1125 * t574, t1350, -t1424, -t1425, t1348, t1408, t1349, -t1121 * t579 + t1125 * t589 - t1404, -t1121 * t597 + t1125 * t618 - t1420, -t1121 * t604 + t1125 * t623 - t1430, -pkin(7) * t580 - t1121 * t544 + t1125 * t553, -t1121 * t710 + t1125 * t738, -t1121 * t649 + t1125 * t680, -t1121 * t673 + t1125 * t701, -t1121 * t709 + t1125 * t737, -t1121 * t674 + t1125 * t702, -t1121 * t732 + t1125 * t768, -pkin(7) * t616 - t1121 * t536 + t1125 * t541, -pkin(7) * t619 - t1121 * t538 + t1125 * t542, -pkin(7) * t607 - t1121 * t526 + t1125 * t533, -pkin(7) * t529 - t1121 * t512 + t1125 * t514, -t1121 * t610 + t1125 * t627, -t1121 * t569 + t1125 * t587, -t1121 * t591 + t1125 * t611, -t1121 * t609 + t1125 * t626, -t1121 * t592 + t1125 * t612, -t1121 * t636 + t1125 * t655, -pkin(7) * t556 - t1121 * t513 + t1125 * t517, -pkin(7) * t558 - t1121 * t515 + t1125 * t518, -pkin(7) * t545 - t1121 * t507 + t1125 * t509, -pkin(7) * t510 - t1121 * t504 + t1125 * t505; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1094, 0, t1095, 0, t1081, t1213, t1053, pkin(7) * t1053, t1035 * t1121 + t1125 * t992, t1121 * t975 + t1125 * t934, t1021 * t1121 + t1125 * t949, t1031 * t1121 + t1125 * t991, t1022 * t1121 + t1125 * t950, t1047 * t1121 - t1125 * t1244, pkin(7) * t905 + t1121 * t852 + t1125 * t840, pkin(7) * t898 + t1121 * t845 + t1125 * t832, pkin(7) * t869 + t1121 * t799 + t1125 * t810, pkin(7) * t817 + t1121 * t808 + t1125 * t818, t1345, t1407, -t1422, t1346, -t1423, t1347, t1121 * t650 + t1125 * t638 + t1419, t1121 * t662 + t1125 * t641 + t1429, t1121 * t601 + t1125 * t596 + t1403, pkin(7) * t606 + t1121 * t574 + t1125 * t568, t1347, t1422, t1423, t1345, t1407, t1346, t1121 * t589 + t1125 * t579 + t1403, t1121 * t618 + t1125 * t597 - t1419, t1121 * t623 + t1125 * t604 - t1429, pkin(7) * t581 + t1121 * t553 + t1125 * t544, t1121 * t738 + t1125 * t710, t1121 * t680 + t1125 * t649, t1121 * t701 + t1125 * t673, t1121 * t737 + t1125 * t709, t1121 * t702 + t1125 * t674, t1121 * t768 + t1125 * t732, pkin(7) * t617 + t1121 * t541 + t1125 * t536, pkin(7) * t620 + t1121 * t542 + t1125 * t538, pkin(7) * t608 + t1121 * t533 + t1125 * t526, pkin(7) * t530 + t1121 * t514 + t1125 * t512, t1121 * t627 + t1125 * t610, t1121 * t587 + t1125 * t569, t1121 * t611 + t1125 * t591, t1121 * t626 + t1125 * t609, t1121 * t612 + t1125 * t592, t1121 * t655 + t1125 * t636, pkin(7) * t557 + t1121 * t517 + t1125 * t513, pkin(7) * t559 + t1121 * t518 + t1125 * t515, pkin(7) * t546 + t1121 * t509 + t1125 * t507, pkin(7) * t511 + t1121 * t505 + t1125 * t504; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1100, t1101, 0, 0, (t1079 * t1115 + t1124 * t1315) * t1120, t1117 * t1085 + t1115 * t1162, t1117 * t1038 + t1115 * t1160, (t1080 * t1115 - t1120 * t1315) * t1124, t1117 * t1037 + t1115 * t1158, t1091, pkin(1) * t951 + t1115 * t1206 - t1117 * t993, pkin(1) * t946 + t1115 * t1207 - t1117 * t994, pkin(1) * t933 + t1115 * t1136, pkin(1) * t896 + t1289 * t911, t1322, t1400, -t1409, t1321, -t1410, t1323, t1115 * t1138 + t1117 * t769 - t1405, t1115 * t1137 + t1117 * t774 - t1426, t1115 * t1141 + t1292 * t1360 + t1396, pkin(1) * t652 + t1115 * t1132 + t1292 * t730, t1323, t1409, t1410, t1322, t1400, t1321, t1115 * t1142 + t1117 * t809 + t1396, t1115 * t1140 + t1117 * t720 + t1405, t1115 * t1139 + t1117 * t727 + t1426, pkin(1) * t625 + t1115 * t1143 + t1117 * t640, t1115 * t1185 + t1117 * t878, t1115 * t1191 + t1117 * t797, t1115 * t1188 + t1117 * t863, t1115 * t1186 + t1117 * t876, t1115 * t1187 + t1117 * t864, t1115 * t1184 + t1117 * t907, pkin(1) * t669 + t1115 * t1145 + t1117 * t632, pkin(1) * t672 + t1115 * t1144 + t1117 * t637, pkin(1) * t646 + t1115 * t1146 + t1117 * t586, pkin(1) * t551 + t1115 * t1149 + t1117 * t552, t1115 * t1196 + t1117 * t744, t1115 * t1202 + t1117 * t686, t1115 * t1199 + t1117 * t707, t1115 * t1197 + t1117 * t743, t1115 * t1198 + t1117 * t708, t1115 * t1194 + t1117 * t794, pkin(1) * t578 + t1115 * t1148 + t1117 * t547, pkin(1) * t585 + t1115 * t1147 + t1117 * t548, pkin(1) * t565 + t1115 * t1150 + t1117 * t537, pkin(1) * t520 + t1115 * t1151 + t1117 * t516;];
tauB_reg  = t1;