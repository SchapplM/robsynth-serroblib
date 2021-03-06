% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 11:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRPRR7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynB_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:39:55
% EndTime: 2019-05-07 11:43:01
% DurationCPUTime: 98.97s
% Computational Cost: add. (1028478->1122), mult. (2284777->1822), div. (0->0), fcn. (1888392->14), ass. (0->810)
t1251 = sin(pkin(6));
t1416 = qJD(1) * t1251;
t1254 = sin(qJ(6));
t1256 = sin(qJ(3));
t1261 = cos(qJ(3));
t1253 = cos(pkin(6));
t1415 = qJD(1) * t1253;
t1357 = qJD(2) + t1415;
t1257 = sin(qJ(2));
t1384 = t1251 * t1257;
t1371 = qJD(1) * t1384;
t1205 = t1256 * t1371 - t1261 * t1357;
t1206 = t1256 * t1357 + t1261 * t1371;
t1250 = sin(pkin(12));
t1252 = cos(pkin(12));
t1171 = -t1205 * t1252 - t1250 * t1206;
t1172 = -t1205 * t1250 + t1206 * t1252;
t1255 = sin(qJ(5));
t1260 = cos(qJ(5));
t1113 = t1171 * t1255 + t1172 * t1260;
t1262 = cos(qJ(2));
t1383 = t1251 * t1262;
t1236 = qJD(1) * t1383 - qJD(3);
t1228 = -qJD(5) + t1236;
t1259 = cos(qJ(6));
t1087 = t1113 * t1254 + t1259 * t1228;
t1089 = t1113 * t1259 - t1228 * t1254;
t1018 = t1089 * t1087;
t1413 = qJD(1) * t1262;
t1359 = qJD(2) * t1413;
t1374 = t1257 * qJDD(1);
t1215 = (t1359 + t1374) * t1251;
t1352 = qJDD(1) * t1253 + qJDD(2);
t1161 = -t1205 * qJD(3) + t1261 * t1215 + t1256 * t1352;
t1276 = -t1256 * t1215 + t1261 * t1352;
t1270 = t1206 * qJD(3) - t1276;
t1102 = t1252 * t1161 - t1250 * t1270;
t1355 = t1161 * t1250 + t1252 * t1270;
t1356 = t1102 * t1255 + t1260 * t1355;
t992 = -qJD(5) * t1113 - t1356;
t1280 = qJDD(6) - t992;
t1442 = -t1018 + t1280;
t1451 = t1254 * t1442;
t1450 = t1259 * t1442;
t1111 = -t1260 * t1171 + t1172 * t1255;
t1099 = t1111 * t1228;
t1279 = qJD(5) * t1111 - t1102 * t1260 + t1255 * t1355;
t1449 = t1279 - t1099;
t1351 = t1357 ^ 2;
t1108 = qJD(6) + t1111;
t1030 = t1108 * t1087;
t1375 = qJDD(1) * t1251;
t1216 = -qJD(2) * t1371 + t1262 * t1375;
t1210 = -qJDD(3) + t1216;
t1208 = -qJDD(5) + t1210;
t937 = -t1087 * qJD(6) - t1254 * t1208 - t1259 * t1279;
t905 = -t937 + t1030;
t1399 = t1171 * t1172;
t1275 = -t1210 + t1399;
t1448 = t1250 * t1275;
t1447 = t1252 * t1275;
t1404 = t1111 * t1113;
t1272 = -t1208 - t1404;
t1446 = t1255 * t1272;
t1397 = t1205 * t1206;
t1274 = -t1210 - t1397;
t1445 = t1256 * t1274;
t1444 = t1260 * t1272;
t1443 = t1261 * t1274;
t1157 = t1171 * t1236;
t1068 = -t1157 - t1102;
t1440 = -t1157 + t1102;
t1189 = t1205 * t1236;
t1140 = t1189 - t1161;
t1138 = t1189 + t1161;
t1247 = t1251 ^ 2;
t1439 = (t1253 ^ 2 + t1247) * t1416 + qJD(2) * t1253 * t1251;
t1358 = t1259 * t1208 - t1254 * t1279;
t900 = (qJD(6) - t1108) * t1089 + t1358;
t954 = (qJD(5) + t1228) * t1113 + t1356;
t1135 = (qJD(3) + t1236) * t1206 - t1276;
t1258 = sin(qJ(1));
t1263 = cos(qJ(1));
t1239 = g(1) * t1263 + t1258 * g(2);
t1437 = qJD(1) ^ 2;
t1211 = -pkin(1) * t1437 + pkin(8) * t1375 - t1239;
t1430 = pkin(2) * t1262;
t1350 = -pkin(9) * t1257 - t1430;
t1214 = t1350 * t1416;
t1238 = t1258 * g(1) - t1263 * g(2);
t1428 = pkin(8) * t1251;
t1273 = qJDD(1) * pkin(1) + t1428 * t1437 + t1238;
t1271 = t1253 * t1273;
t1269 = -g(3) * t1384 + t1257 * t1271;
t1126 = t1352 * pkin(9) - t1351 * pkin(2) + (t1214 * t1416 + t1211) * t1262 + t1269;
t1320 = qJD(1) * t1357;
t1298 = t1262 * t1320;
t1299 = t1257 * t1320;
t1427 = t1253 * g(3);
t1127 = -t1216 * pkin(2) - t1215 * pkin(9) - t1427 + (pkin(2) * t1299 - pkin(9) * t1298 - t1273) * t1251;
t1047 = t1126 * t1256 - t1261 * t1127;
t1438 = qJ(4) * t1140 - t1047;
t1085 = t1087 ^ 2;
t1086 = t1089 ^ 2;
t1107 = t1108 ^ 2;
t1109 = t1111 ^ 2;
t1110 = t1113 ^ 2;
t1436 = t1171 ^ 2;
t1170 = t1172 ^ 2;
t1435 = t1205 ^ 2;
t1204 = t1206 ^ 2;
t1434 = t1228 ^ 2;
t1433 = t1236 ^ 2;
t1432 = 2 * qJD(4);
t1431 = pkin(2) * t1257;
t1429 = pkin(5) * t1255;
t1037 = pkin(5) * t1111 - pkin(11) * t1113;
t1048 = t1261 * t1126 + t1256 * t1127;
t1181 = -pkin(3) * t1236 - qJ(4) * t1206;
t1267 = -qJ(4) * t1270 + t1236 * t1181 + t1048;
t1266 = -pkin(3) * t1435 + t1267;
t1268 = t1274 * pkin(3) + t1438;
t1373 = t1172 * t1432;
t915 = t1250 * t1266 - t1252 * t1268 + t1373;
t1265 = t1275 * pkin(4) + pkin(10) * t1068 - t915;
t1151 = -pkin(4) * t1236 - pkin(10) * t1172;
t916 = t1171 * t1432 + t1250 * t1268 + t1252 * t1266;
t864 = -pkin(4) * t1436 - pkin(10) * t1355 + t1151 * t1236 + t916;
t780 = t1255 * t1265 + t1260 * t864;
t769 = -pkin(5) * t1434 - pkin(11) * t1208 - t1037 * t1111 + t780;
t1354 = t1211 * t1257 - t1262 * t1271;
t1414 = qJD(1) * t1257;
t1125 = -t1352 * pkin(2) - t1351 * pkin(9) + (g(3) * t1262 + t1214 * t1414) * t1251 + t1354;
t1031 = t1270 * pkin(3) - t1435 * qJ(4) + t1181 * t1206 + qJDD(4) + t1125;
t935 = pkin(4) * t1355 - t1436 * pkin(10) + t1151 * t1172 + t1031;
t838 = t935 + (-t1113 * t1228 - t992) * pkin(5) + t1449 * pkin(11);
t730 = t1254 * t838 + t1259 * t769;
t779 = t1255 * t864 - t1260 * t1265;
t725 = t1255 * t780 - t1260 * t779;
t1426 = t1250 * t725;
t1425 = t1252 * t725;
t768 = t1208 * pkin(5) - pkin(11) * t1434 + t1037 * t1113 + t779;
t1424 = t1254 * t768;
t924 = t1018 + t1280;
t1423 = t1254 * t924;
t1422 = t1255 * t935;
t832 = t1250 * t916 - t1252 * t915;
t1421 = t1256 * t832;
t1420 = t1259 * t768;
t1419 = t1259 * t924;
t1418 = t1260 * t935;
t1417 = t1261 * t832;
t1027 = t1208 - t1404;
t1412 = t1027 * t1255;
t1411 = t1027 * t1260;
t1410 = t1031 * t1250;
t1409 = t1031 * t1252;
t1095 = t1210 + t1399;
t1408 = t1095 * t1250;
t1407 = t1095 * t1252;
t1406 = t1108 * t1254;
t1405 = t1108 * t1259;
t1403 = t1125 * t1256;
t1402 = t1125 * t1261;
t1152 = t1210 - t1397;
t1401 = t1152 * t1256;
t1400 = t1152 * t1261;
t1398 = t1172 * t1236;
t1396 = t1210 * t1257;
t1385 = t1247 * t1437;
t1235 = t1262 * t1257 * t1385;
t1212 = -t1235 + t1352;
t1395 = t1212 * t1257;
t1394 = t1212 * t1262;
t1213 = t1235 + t1352;
t1393 = t1213 * t1257;
t1392 = t1213 * t1262;
t1391 = t1228 * t1255;
t1390 = t1228 * t1260;
t1389 = t1236 * t1250;
t1388 = t1236 * t1252;
t1387 = t1236 * t1256;
t1386 = t1236 * t1261;
t1382 = t1253 * t1257;
t1191 = t1251 * t1273 + t1427;
t1381 = t1257 * t1191;
t1380 = t1262 * t1191;
t1248 = t1257 ^ 2;
t1249 = t1262 ^ 2;
t1376 = t1248 + t1249;
t1372 = -pkin(5) * t1260 - pkin(4);
t1369 = t1255 * t1018;
t1368 = t1260 * t1018;
t1367 = t1257 * t1404;
t1366 = t1262 * t1404;
t1365 = t1257 * t1399;
t1364 = t1262 * t1399;
t1363 = t1257 * t1397;
t1362 = t1262 * t1397;
t1361 = t1248 * t1385;
t1360 = t1249 * t1385;
t833 = t1250 * t915 + t1252 * t916;
t729 = t1254 * t769 - t1259 * t838;
t726 = t1255 * t779 + t1260 * t780;
t968 = t1047 * t1256 + t1261 * t1048;
t1193 = -t1238 * t1258 - t1263 * t1239;
t1234 = qJDD(1) * t1263 - t1258 * t1437;
t1349 = -pkin(7) * t1234 - g(3) * t1258;
t1201 = -t1361 - t1351;
t1173 = -t1201 * t1257 - t1394;
t1348 = pkin(8) * t1173 - t1381;
t1220 = -t1351 - t1360;
t1178 = t1220 * t1262 - t1393;
t1347 = pkin(8) * t1178 + t1380;
t674 = t1254 * t730 - t1259 * t729;
t675 = t1254 * t729 + t1259 * t730;
t651 = t1255 * t675 - t1260 * t768;
t652 = t1255 * t768 + t1260 * t675;
t623 = t1250 * t652 + t1252 * t651;
t624 = -t1250 * t651 + t1252 * t652;
t597 = -t1256 * t623 + t1261 * t624;
t1346 = t1257 * t597 - t1262 * t674;
t666 = t1250 * t726 + t1425;
t667 = t1252 * t726 - t1426;
t632 = -t1256 * t666 + t1261 * t667;
t1345 = t1257 * t632 - t1262 * t935;
t904 = -t1030 - t937;
t827 = -t1254 * t904 - t1259 * t900;
t973 = t1085 + t1086;
t795 = t1255 * t827 + t1260 * t973;
t796 = -t1255 * t973 + t1260 * t827;
t736 = t1250 * t796 + t1252 * t795;
t737 = -t1250 * t795 + t1252 * t796;
t677 = -t1256 * t736 + t1261 * t737;
t825 = -t1254 * t900 + t1259 * t904;
t1344 = t1257 * t677 - t1262 * t825;
t1016 = -t1086 + t1085;
t901 = (-qJD(6) - t1108) * t1089 - t1358;
t826 = t1254 * t905 + t1259 * t901;
t801 = t1016 * t1260 + t1255 * t826;
t802 = -t1016 * t1255 + t1260 * t826;
t740 = t1250 * t802 + t1252 * t801;
t741 = -t1250 * t801 + t1252 * t802;
t682 = -t1256 * t740 + t1261 * t741;
t824 = -t1254 * t901 + t1259 * t905;
t1343 = t1257 * t682 + t1262 * t824;
t988 = -t1107 - t1085;
t866 = t1259 * t988 - t1451;
t803 = t1255 * t866 + t1260 * t901;
t804 = -t1255 * t901 + t1260 * t866;
t743 = t1250 * t804 + t1252 * t803;
t744 = -t1250 * t803 + t1252 * t804;
t684 = -t1256 * t743 + t1261 * t744;
t865 = t1254 * t988 + t1450;
t1342 = t1257 * t684 - t1262 * t865;
t996 = -t1086 - t1107;
t874 = -t1254 * t996 - t1419;
t806 = t1255 * t874 + t1260 * t905;
t807 = -t1255 * t905 + t1260 * t874;
t747 = t1250 * t807 + t1252 * t806;
t748 = -t1250 * t806 + t1252 * t807;
t686 = -t1256 * t747 + t1261 * t748;
t873 = t1259 * t996 - t1423;
t1341 = t1257 * t686 - t1262 * t873;
t1026 = -t1086 + t1107;
t889 = -t1026 * t1254 + t1450;
t817 = t1255 * t889 + t1260 * t904;
t819 = -t1255 * t904 + t1260 * t889;
t754 = t1250 * t819 + t1252 * t817;
t756 = -t1250 * t817 + t1252 * t819;
t693 = -t1256 * t754 + t1261 * t756;
t887 = -t1026 * t1259 - t1451;
t1340 = t1257 * t693 + t1262 * t887;
t1025 = t1085 - t1107;
t890 = t1025 * t1259 - t1423;
t818 = t1255 * t890 + t1260 * t900;
t820 = -t1255 * t900 + t1260 * t890;
t755 = t1250 * t820 + t1252 * t818;
t757 = -t1250 * t818 + t1252 * t820;
t694 = -t1256 * t755 + t1261 * t757;
t888 = -t1025 * t1254 - t1419;
t1339 = t1257 * t694 + t1262 * t888;
t936 = -qJD(6) * t1089 - t1358;
t897 = t1087 * t1405 - t1254 * t936;
t854 = t1255 * t897 + t1368;
t856 = t1260 * t897 - t1369;
t775 = t1250 * t856 + t1252 * t854;
t777 = -t1250 * t854 + t1252 * t856;
t723 = -t1256 * t775 + t1261 * t777;
t896 = -t1087 * t1406 - t1259 * t936;
t1338 = t1257 * t723 + t1262 * t896;
t899 = -t1089 * t1406 + t1259 * t937;
t855 = t1255 * t899 - t1368;
t857 = t1260 * t899 + t1369;
t776 = t1250 * t857 + t1252 * t855;
t778 = -t1250 * t855 + t1252 * t857;
t724 = -t1256 * t776 + t1261 * t778;
t898 = -t1089 * t1405 - t1254 * t937;
t1337 = t1257 * t724 + t1262 * t898;
t952 = (-t1087 * t1259 + t1089 * t1254) * t1108;
t892 = t1255 * t952 - t1260 * t1280;
t893 = t1255 * t1280 + t1260 * t952;
t808 = t1250 * t893 + t1252 * t892;
t809 = -t1250 * t892 + t1252 * t893;
t751 = -t1256 * t808 + t1261 * t809;
t951 = (t1087 * t1254 + t1089 * t1259) * t1108;
t1336 = t1257 * t751 + t1262 * t951;
t1034 = -t1434 - t1109;
t959 = t1034 * t1255 + t1444;
t960 = t1034 * t1260 - t1446;
t879 = t1250 * t960 + t1252 * t959;
t880 = -t1250 * t959 + t1252 * t960;
t794 = -t1256 * t879 + t1261 * t880;
t953 = (qJD(5) - t1228) * t1113 + t1356;
t1335 = t1257 * t794 - t1262 * t953;
t1079 = -t1110 - t1434;
t978 = t1079 * t1260 + t1412;
t979 = -t1079 * t1255 + t1411;
t894 = t1250 * t979 + t1252 * t978;
t895 = -t1250 * t978 + t1252 * t979;
t823 = -t1256 * t894 + t1261 * t895;
t1334 = t1257 * t823 + t1262 * t1449;
t1094 = -t1110 + t1434;
t984 = t1094 * t1260 + t1446;
t986 = -t1094 * t1255 + t1444;
t911 = t1250 * t986 + t1252 * t984;
t913 = -t1250 * t984 + t1252 * t986;
t830 = -t1256 * t911 + t1261 * t913;
t958 = t1099 + t1279;
t1333 = t1257 * t830 + t1262 * t958;
t1093 = t1109 - t1434;
t985 = t1093 * t1255 - t1411;
t987 = t1093 * t1260 + t1412;
t912 = t1250 * t987 + t1252 * t985;
t914 = -t1250 * t985 + t1252 * t987;
t831 = -t1256 * t912 + t1261 * t914;
t1332 = t1257 * t831 + t1262 * t954;
t1000 = -t1109 - t1110;
t876 = -t1255 * t954 + t1260 * t958;
t878 = -t1255 * t958 - t1260 * t954;
t790 = t1250 * t878 + t1252 * t876;
t792 = -t1250 * t876 + t1252 * t878;
t735 = -t1256 * t790 + t1261 * t792;
t1331 = -t1000 * t1262 + t1257 * t735;
t761 = t1261 * t833 - t1421;
t1330 = -t1031 * t1262 + t1257 * t761;
t1039 = -t1110 + t1109;
t875 = -t1255 * t953 - t1260 * t1449;
t877 = t1255 * t1449 - t1260 * t953;
t789 = t1250 * t877 + t1252 * t875;
t791 = -t1250 * t875 + t1252 * t877;
t734 = -t1256 * t789 + t1261 * t791;
t1329 = t1039 * t1262 + t1257 * t734;
t1063 = t1355 - t1398;
t1104 = -t1433 - t1436;
t1023 = t1104 * t1250 + t1447;
t1024 = t1104 * t1252 - t1448;
t943 = -t1023 * t1256 + t1024 * t1261;
t1328 = -t1063 * t1262 + t1257 * t943;
t1064 = t1355 + t1398;
t1155 = -t1433 + t1436;
t1052 = t1155 * t1250 - t1407;
t1054 = t1155 * t1252 + t1408;
t972 = -t1052 * t1256 + t1054 * t1261;
t1327 = t1064 * t1262 + t1257 * t972;
t1145 = -t1170 - t1433;
t1045 = t1145 * t1252 + t1408;
t1046 = -t1145 * t1250 + t1407;
t966 = -t1045 * t1256 + t1046 * t1261;
t1326 = t1257 * t966 - t1262 * t1440;
t1156 = -t1170 + t1433;
t1051 = t1156 * t1252 + t1448;
t1053 = -t1156 * t1250 + t1447;
t971 = -t1051 * t1256 + t1053 * t1261;
t1325 = t1068 * t1262 + t1257 * t971;
t1072 = -t1170 - t1436;
t981 = -t1064 * t1250 + t1068 * t1252;
t983 = -t1064 * t1252 - t1068 * t1250;
t910 = -t1256 * t981 + t1261 * t983;
t1324 = -t1072 * t1262 + t1257 * t910;
t1114 = -t1170 + t1436;
t980 = -t1063 * t1250 + t1252 * t1440;
t982 = -t1063 * t1252 - t1250 * t1440;
t909 = -t1256 * t980 + t1261 * t982;
t1323 = t1114 * t1262 + t1257 * t909;
t1322 = -t1125 * t1262 + t1257 * t968;
t1021 = (t1111 * t1255 + t1113 * t1260) * t1228;
t1022 = (t1111 * t1260 - t1113 * t1255) * t1228;
t940 = t1021 * t1252 + t1022 * t1250;
t941 = -t1021 * t1250 + t1022 * t1252;
t861 = -t1256 * t940 + t1261 * t941;
t1321 = t1208 * t1262 + t1257 * t861;
t1319 = t1247 * t1257 * t1359;
t967 = -t1047 * t1261 + t1048 * t1256;
t1136 = (-qJD(3) + t1236) * t1206 + t1276;
t1061 = t1136 * t1261 - t1138 * t1256;
t1174 = -t1204 + t1435;
t1318 = t1061 * t1257 + t1174 * t1262;
t1062 = -t1135 * t1261 - t1140 * t1256;
t1150 = t1204 + t1435;
t1317 = t1062 * t1257 + t1150 * t1262;
t1164 = -t1433 - t1435;
t1091 = t1164 * t1261 - t1445;
t1316 = t1091 * t1257 + t1136 * t1262;
t1175 = -t1204 - t1433;
t1106 = -t1175 * t1256 + t1400;
t1315 = t1106 * t1257 - t1138 * t1262;
t1183 = -t1204 + t1433;
t1117 = -t1183 * t1256 + t1443;
t1314 = t1117 * t1257 + t1140 * t1262;
t1182 = -t1433 + t1435;
t1118 = t1182 * t1261 + t1401;
t1313 = t1118 * t1257 + t1135 * t1262;
t1166 = g(3) * t1383 + t1354;
t1167 = t1262 * t1211 + t1269;
t1312 = -t1166 * t1262 + t1167 * t1257;
t1103 = t1166 * t1257 + t1167 * t1262;
t1224 = t1251 * t1298;
t1185 = t1224 + t1215;
t1223 = t1251 * t1299;
t1188 = t1216 - t1223;
t1311 = t1185 * t1262 + t1188 * t1257;
t1186 = -t1224 + t1215;
t1187 = t1216 + t1223;
t1310 = -t1186 * t1262 + t1187 * t1257;
t1309 = t1201 * t1262 - t1395;
t1219 = -t1351 + t1360;
t1308 = t1219 * t1257 + t1394;
t1218 = t1351 - t1361;
t1307 = t1218 * t1262 + t1393;
t1306 = t1220 * t1257 + t1392;
t1192 = t1238 * t1263 - t1258 * t1239;
t1305 = t1251 * t1352;
t947 = -t1111 * t1391 + t1260 * t992;
t948 = -t1111 * t1390 - t1255 * t992;
t869 = t1250 * t948 + t1252 * t947;
t871 = -t1250 * t947 + t1252 * t948;
t787 = -t1256 * t869 + t1261 * t871;
t1304 = t1257 * t787 + t1366;
t949 = -t1113 * t1390 - t1255 * t1279;
t950 = t1113 * t1391 - t1260 * t1279;
t870 = t1250 * t950 + t1252 * t949;
t872 = -t1250 * t949 + t1252 * t950;
t788 = -t1256 * t870 + t1261 * t872;
t1303 = t1257 * t788 - t1366;
t1055 = t1171 * t1389 - t1252 * t1355;
t1056 = t1171 * t1388 + t1250 * t1355;
t976 = -t1055 * t1256 + t1056 * t1261;
t1302 = t1257 * t976 - t1364;
t1057 = t1102 * t1250 - t1172 * t1388;
t1058 = t1102 * t1252 + t1172 * t1389;
t977 = -t1057 * t1256 + t1058 * t1261;
t1301 = t1257 * t977 + t1364;
t1300 = t1251 * t1320;
t1130 = -t1205 * t1386 + t1256 * t1270;
t1297 = t1130 * t1257 + t1362;
t1132 = t1161 * t1261 + t1206 * t1387;
t1296 = t1132 * t1257 - t1362;
t598 = pkin(10) * t652 + (-pkin(11) * t1255 + t1372) * t674;
t611 = -pkin(10) * t651 + (-pkin(11) * t1260 + t1429) * t674;
t580 = -pkin(3) * t674 + qJ(4) * t624 + t1250 * t611 + t1252 * t598;
t582 = -qJ(4) * t623 - t1250 * t598 + t1252 * t611;
t596 = t1256 * t624 + t1261 * t623;
t567 = -pkin(9) * t596 - t1256 * t580 + t1261 * t582;
t581 = -pkin(2) * t596 - pkin(3) * t623 - pkin(4) * t651 + pkin(5) * t768 - pkin(11) * t675;
t595 = t1257 * t674 + t1262 * t597;
t1295 = pkin(8) * t595 + t1257 * t567 + t1262 * t581;
t659 = -pkin(11) * t825 - t674;
t640 = pkin(10) * t796 + t1255 * t659 + t1372 * t825;
t644 = -pkin(10) * t795 + t1260 * t659 + t1429 * t825;
t601 = -pkin(3) * t825 + qJ(4) * t737 + t1250 * t644 + t1252 * t640;
t603 = -qJ(4) * t736 - t1250 * t640 + t1252 * t644;
t676 = t1256 * t737 + t1261 * t736;
t585 = -pkin(9) * t676 - t1256 * t601 + t1261 * t603;
t618 = -pkin(2) * t676 - pkin(3) * t736 - pkin(4) * t795 - pkin(5) * t973 - pkin(11) * t827 - t675;
t660 = t1257 * t825 + t1262 * t677;
t1294 = pkin(8) * t660 + t1257 * t585 + t1262 * t618;
t707 = -pkin(5) * t865 + t729;
t752 = -pkin(11) * t865 + t1424;
t649 = -pkin(4) * t865 + pkin(10) * t804 + t1255 * t752 + t1260 * t707;
t653 = -pkin(10) * t803 - t1255 * t707 + t1260 * t752;
t612 = -pkin(3) * t865 + qJ(4) * t744 + t1250 * t653 + t1252 * t649;
t614 = -qJ(4) * t743 - t1250 * t649 + t1252 * t653;
t683 = t1256 * t744 + t1261 * t743;
t588 = -pkin(9) * t683 - t1256 * t612 + t1261 * t614;
t633 = -pkin(2) * t683 - pkin(3) * t743 - pkin(4) * t803 - pkin(5) * t901 - pkin(11) * t866 + t1420;
t671 = t1257 * t865 + t1262 * t684;
t1293 = pkin(8) * t671 + t1257 * t588 + t1262 * t633;
t708 = -pkin(5) * t873 + t730;
t753 = -pkin(11) * t873 + t1420;
t650 = -pkin(4) * t873 + pkin(10) * t807 + t1255 * t753 + t1260 * t708;
t654 = -pkin(10) * t806 - t1255 * t708 + t1260 * t753;
t613 = -pkin(3) * t873 + qJ(4) * t748 + t1250 * t654 + t1252 * t650;
t615 = -qJ(4) * t747 - t1250 * t650 + t1252 * t654;
t685 = t1256 * t748 + t1261 * t747;
t589 = -pkin(9) * t685 - t1256 * t613 + t1261 * t615;
t634 = -pkin(2) * t685 - pkin(3) * t747 - pkin(4) * t806 - pkin(5) * t905 - pkin(11) * t874 - t1424;
t673 = t1257 * t873 + t1262 * t686;
t1292 = pkin(8) * t673 + t1257 * t589 + t1262 * t634;
t714 = -pkin(4) * t935 + pkin(10) * t726;
t626 = -pkin(3) * t935 - pkin(10) * t1426 + qJ(4) * t667 + t1252 * t714;
t628 = -pkin(10) * t1425 - qJ(4) * t666 - t1250 * t714;
t631 = t1256 * t667 + t1261 * t666;
t592 = -pkin(9) * t631 - t1256 * t626 + t1261 * t628;
t606 = -pkin(2) * t631 - pkin(3) * t666 - pkin(4) * t725;
t627 = t1257 * t935 + t1262 * t632;
t1291 = pkin(8) * t627 + t1257 * t592 + t1262 * t606;
t704 = -pkin(4) * t1000 + pkin(10) * t878 + t726;
t705 = -pkin(10) * t876 - t725;
t645 = -pkin(3) * t1000 + qJ(4) * t792 + t1250 * t705 + t1252 * t704;
t646 = -qJ(4) * t790 - t1250 * t704 + t1252 * t705;
t733 = t1256 * t792 + t1261 * t790;
t605 = -pkin(9) * t733 - t1256 * t645 + t1261 * t646;
t690 = -pkin(2) * t733 - pkin(3) * t790 - pkin(4) * t876;
t718 = t1000 * t1257 + t1262 * t735;
t1290 = pkin(8) * t718 + t1257 * t605 + t1262 * t690;
t837 = -pkin(4) * t953 + pkin(10) * t960 - t1418;
t862 = -pkin(10) * t959 + t1422;
t728 = -pkin(3) * t953 + qJ(4) * t880 + t1250 * t862 + t1252 * t837;
t749 = -qJ(4) * t879 - t1250 * t837 + t1252 * t862;
t793 = t1256 * t880 + t1261 * t879;
t658 = -pkin(9) * t793 - t1256 * t728 + t1261 * t749;
t703 = -pkin(2) * t793 - pkin(3) * t879 - pkin(4) * t959 + t779;
t767 = t1257 * t953 + t1262 * t794;
t1289 = pkin(8) * t767 + t1257 * t658 + t1262 * t703;
t844 = pkin(4) * t1449 + pkin(10) * t979 + t1422;
t883 = -pkin(10) * t978 + t1418;
t742 = pkin(3) * t1449 + qJ(4) * t895 + t1250 * t883 + t1252 * t844;
t759 = -qJ(4) * t894 - t1250 * t844 + t1252 * t883;
t822 = t1256 * t895 + t1261 * t894;
t668 = -pkin(9) * t822 - t1256 * t742 + t1261 * t759;
t712 = -pkin(2) * t822 - pkin(3) * t894 - pkin(4) * t978 + t780;
t781 = -t1257 * t1449 + t1262 * t823;
t1288 = pkin(8) * t781 + t1257 * t668 + t1262 * t712;
t760 = t1256 * t833 + t1417;
t805 = -pkin(3) * t1031 + qJ(4) * t833;
t698 = -pkin(9) * t760 - qJ(4) * t1417 - t1256 * t805;
t715 = -pkin(2) * t760 - pkin(3) * t832;
t758 = t1031 * t1257 + t1262 * t761;
t1287 = pkin(8) * t758 + t1257 * t698 + t1262 * t715;
t782 = -pkin(3) * t1072 + qJ(4) * t983 + t833;
t799 = -qJ(4) * t981 - t832;
t908 = t1256 * t983 + t1261 * t981;
t713 = -pkin(9) * t908 - t1256 * t782 + t1261 * t799;
t850 = -pkin(2) * t908 - pkin(3) * t981;
t884 = t1072 * t1257 + t1262 * t910;
t1286 = pkin(8) * t884 + t1257 * t713 + t1262 * t850;
t920 = -pkin(3) * t1063 + qJ(4) * t1024 - t1409;
t942 = t1023 * t1261 + t1024 * t1256;
t946 = -qJ(4) * t1023 + t1410;
t814 = -pkin(9) * t942 - t1256 * t920 + t1261 * t946;
t839 = -pkin(2) * t942 + t1250 * t1267 - t1252 * t1438 + t1373 + (-t1250 * t1435 - t1252 * t1274 - t1023) * pkin(3);
t919 = t1063 * t1257 + t1262 * t943;
t1285 = pkin(8) * t919 + t1257 * t814 + t1262 * t839;
t922 = -pkin(3) * t1440 + qJ(4) * t1046 + t1410;
t963 = -qJ(4) * t1045 + t1409;
t965 = t1045 * t1261 + t1046 * t1256;
t834 = -pkin(9) * t965 - t1256 * t922 + t1261 * t963;
t845 = -pkin(2) * t965 - pkin(3) * t1045 + t916;
t926 = t1257 * t1440 + t1262 * t966;
t1284 = pkin(8) * t926 + t1257 * t834 + t1262 * t845;
t1090 = t1164 * t1256 + t1443;
t1001 = -pkin(2) * t1090 + t1047;
t1035 = -pkin(9) * t1090 + t1403;
t1040 = t1091 * t1262 - t1136 * t1257;
t1283 = pkin(8) * t1040 + t1001 * t1262 + t1035 * t1257;
t1105 = t1175 * t1261 + t1401;
t1004 = -pkin(2) * t1105 + t1048;
t1041 = -pkin(9) * t1105 + t1402;
t1042 = t1106 * t1262 + t1138 * t1257;
t1282 = pkin(8) * t1042 + t1004 * t1262 + t1041 * t1257;
t1146 = t1186 * t1257 + t1187 * t1262;
t1281 = pkin(8) * t1146 + t1103;
t1017 = t1062 * t1262 - t1150 * t1257;
t1060 = -t1135 * t1256 + t1140 * t1261;
t927 = -pkin(9) * t1060 - t967;
t1278 = pkin(8) * t1017 - t1060 * t1430 + t1257 * t927;
t938 = t1125 * t1257 + t1262 * t968;
t1277 = pkin(8) * t938 + t1350 * t967;
t1233 = t1258 * qJDD(1) + t1263 * t1437;
t1222 = t1376 * t1385;
t1221 = (t1248 - t1249) * t1385;
t1217 = -pkin(7) * t1233 + g(3) * t1263;
t1195 = t1253 * t1262 * t1210;
t1194 = t1210 * t1383;
t1190 = t1357 * t1376 * t1416;
t1184 = (t1374 + (0.2e1 * qJD(2) + t1415) * t1413) * t1251;
t1180 = t1262 * t1215 - t1248 * t1300;
t1179 = -t1257 * t1216 - t1249 * t1300;
t1177 = t1219 * t1262 - t1395;
t1176 = -t1218 * t1257 + t1392;
t1163 = (t1253 * t1215 + t1413 * t1439) * t1257;
t1162 = (t1253 * t1216 - t1414 * t1439) * t1262;
t1149 = (t1205 * t1261 - t1206 * t1256) * t1236;
t1148 = (t1205 * t1256 + t1206 * t1261) * t1236;
t1147 = -t1185 * t1257 + t1188 * t1262;
t1144 = t1188 * t1251 + t1253 * t1306;
t1143 = -t1187 * t1251 + t1253 * t1308;
t1142 = -t1186 * t1251 + t1253 * t1307;
t1141 = -t1188 * t1253 + t1251 * t1306;
t1134 = -t1184 * t1251 + t1253 * t1309;
t1133 = t1184 * t1253 + t1251 * t1309;
t1131 = t1161 * t1256 - t1206 * t1386;
t1129 = -t1205 * t1387 - t1261 * t1270;
t1128 = t1149 * t1262 - t1396;
t1124 = -t1221 * t1251 + t1253 * t1311;
t1123 = t1222 * t1251 + t1253 * t1310;
t1122 = -t1222 * t1253 + t1251 * t1310;
t1116 = t1182 * t1256 - t1400;
t1115 = t1183 * t1261 + t1445;
t1084 = (-t1171 * t1252 - t1172 * t1250) * t1236;
t1083 = (-t1171 * t1250 + t1172 * t1252) * t1236;
t1081 = -t1258 * t1144 + t1178 * t1263;
t1080 = t1144 * t1263 + t1258 * t1178;
t1078 = t1132 * t1262 + t1363;
t1077 = t1130 * t1262 - t1363;
t1076 = -t1258 * t1134 + t1173 * t1263;
t1075 = t1134 * t1263 + t1258 * t1173;
t1074 = t1191 * t1251 + t1253 * t1312;
t1073 = -t1191 * t1253 + t1251 * t1312;
t1071 = -t1148 * t1251 + t1149 * t1382 + t1195;
t1070 = -t1258 * t1123 + t1146 * t1263;
t1069 = t1123 * t1263 + t1258 * t1146;
t1059 = t1136 * t1256 + t1138 * t1261;
t1050 = t1118 * t1262 - t1135 * t1257;
t1049 = t1117 * t1262 - t1140 * t1257;
t1043 = -t1381 + (-t1141 * t1251 - t1144 * t1253) * pkin(8);
t1038 = -t1380 + (-t1133 * t1251 - t1134 * t1253) * pkin(8);
t1036 = -pkin(1) * t1141 + t1166 * t1251 + t1253 * t1347;
t1033 = t1061 * t1262 - t1174 * t1257;
t1032 = -pkin(1) * t1133 + t1167 * t1251 + t1253 * t1348;
t1020 = -t1131 * t1251 + t1253 * t1296;
t1019 = -t1129 * t1251 + t1253 * t1297;
t1015 = -t1083 * t1256 + t1084 * t1261;
t1014 = t1083 * t1261 + t1084 * t1256;
t1013 = pkin(8) * t1103 * t1253 - pkin(1) * t1073;
t1012 = -t1258 * t1074 + t1103 * t1263;
t1011 = t1074 * t1263 + t1258 * t1103;
t1010 = t1015 * t1262 - t1396;
t1009 = -pkin(1) * t1122 + t1253 * t1281;
t1008 = -pkin(2) * t1138 + pkin(9) * t1106 + t1403;
t1007 = -t1251 * t1116 + t1253 * t1313;
t1006 = -t1251 * t1115 + t1253 * t1314;
t1005 = (-t1073 * t1251 - t1074 * t1253) * pkin(8);
t1003 = pkin(2) * t1136 + pkin(9) * t1091 - t1402;
t1002 = (-t1122 * t1251 - t1123 * t1253) * pkin(8) - t1312;
t999 = -t1251 * t1105 + t1253 * t1315;
t998 = t1253 * t1105 + t1251 * t1315;
t995 = -t1251 * t1090 + t1253 * t1316;
t994 = t1253 * t1090 + t1251 * t1316;
t975 = t1057 * t1261 + t1058 * t1256;
t974 = t1055 * t1261 + t1056 * t1256;
t970 = t1052 * t1261 + t1054 * t1256;
t969 = t1051 * t1261 + t1053 * t1256;
t964 = -t1251 * t1059 + t1253 * t1318;
t962 = -t1251 * t1060 + t1253 * t1317;
t961 = t1253 * t1060 + t1251 * t1317;
t945 = t1262 * t977 - t1365;
t944 = t1262 * t976 + t1365;
t939 = -pkin(2) * t1125 + pkin(9) * t968;
t934 = t1042 * t1263 - t1258 * t999;
t933 = t1258 * t1042 + t1263 * t999;
t932 = -t1014 * t1251 + t1015 * t1382 + t1195;
t931 = -t1064 * t1257 + t1262 * t972;
t930 = -t1068 * t1257 + t1262 * t971;
t929 = t1040 * t1263 - t1258 * t995;
t928 = t1258 * t1040 + t1263 * t995;
t921 = pkin(2) * t1150 + pkin(9) * t1062 + t968;
t918 = t1017 * t1263 - t1258 * t962;
t917 = t1258 * t1017 + t1263 * t962;
t907 = t1256 * t982 + t1261 * t980;
t891 = -t1114 * t1257 + t1262 * t909;
t886 = -t1251 * t975 + t1253 * t1301;
t885 = -t1251 * t974 + t1253 * t1302;
t882 = -t1251 * t967 + t1253 * t1322;
t881 = t1251 * t1322 + t1253 * t967;
t868 = -t1251 * t970 + t1253 * t1327;
t867 = -t1251 * t969 + t1253 * t1325;
t860 = t1256 * t941 + t1261 * t940;
t859 = -t1251 * t965 + t1253 * t1326;
t858 = t1251 * t1326 + t1253 * t965;
t853 = -t1208 * t1257 + t1262 * t861;
t849 = -t1251 * t942 + t1253 * t1328;
t848 = t1251 * t1328 + t1253 * t942;
t847 = -t1004 * t1257 + t1041 * t1262 + (-t1251 * t998 - t1253 * t999) * pkin(8);
t846 = -t1001 * t1257 + t1035 * t1262 + (-t1251 * t994 - t1253 * t995) * pkin(8);
t843 = -pkin(1) * t998 - t1008 * t1251 + t1253 * t1282;
t842 = -t1258 * t882 + t1263 * t938;
t841 = t1258 * t938 + t1263 * t882;
t840 = -pkin(1) * t994 - t1003 * t1251 + t1253 * t1283;
t829 = t1256 * t914 + t1261 * t912;
t828 = t1256 * t913 + t1261 * t911;
t821 = -t1251 * t907 + t1253 * t1323;
t816 = -t1258 * t859 + t1263 * t926;
t815 = t1258 * t926 + t1263 * t859;
t813 = -t1251 * t908 + t1253 * t1324;
t812 = t1251 * t1324 + t1253 * t908;
t811 = -pkin(2) * t1440 + pkin(9) * t966 + t1256 * t963 + t1261 * t922;
t810 = t1060 * t1431 + t1262 * t927 + (-t1251 * t961 - t1253 * t962) * pkin(8);
t800 = -pkin(2) * t1063 + pkin(9) * t943 + t1256 * t946 + t1261 * t920;
t798 = -t1258 * t849 + t1263 * t919;
t797 = t1258 * t919 + t1263 * t849;
t786 = t1256 * t872 + t1261 * t870;
t785 = t1256 * t871 + t1261 * t869;
t784 = -t1257 * t954 + t1262 * t831;
t783 = -t1257 * t958 + t1262 * t830;
t773 = -t1251 * t860 + t1253 * t1321;
t772 = -pkin(1) * t961 - t1251 * t921 + t1253 * t1278;
t771 = t1262 * t788 + t1367;
t770 = t1262 * t787 - t1367;
t765 = -t1258 * t813 + t1263 * t884;
t764 = t1258 * t884 + t1263 * t813;
t763 = (-pkin(9) * t1262 + t1431) * t967 + (-t1251 * t881 - t1253 * t882) * pkin(8);
t762 = -pkin(1) * t881 - t1251 * t939 + t1253 * t1277;
t750 = t1256 * t809 + t1261 * t808;
t746 = -t1251 * t829 + t1253 * t1332;
t745 = -t1251 * t828 + t1253 * t1333;
t739 = -t1251 * t822 + t1253 * t1334;
t738 = t1251 * t1334 + t1253 * t822;
t732 = t1256 * t791 + t1261 * t789;
t731 = -t1257 * t951 + t1262 * t751;
t727 = -t1039 * t1257 + t1262 * t734;
t722 = t1256 * t778 + t1261 * t776;
t721 = t1256 * t777 + t1261 * t775;
t720 = -t1251 * t786 + t1253 * t1303;
t719 = -t1251 * t785 + t1253 * t1304;
t717 = -t1251 * t793 + t1253 * t1335;
t716 = t1251 * t1335 + t1253 * t793;
t711 = -pkin(2) * t1072 + pkin(9) * t910 + t1256 * t799 + t1261 * t782;
t710 = -t1257 * t898 + t1262 * t724;
t709 = -t1257 * t896 + t1262 * t723;
t706 = -t1257 * t845 + t1262 * t834 + (-t1251 * t858 - t1253 * t859) * pkin(8);
t702 = -t1257 * t839 + t1262 * t814 + (-t1251 * t848 - t1253 * t849) * pkin(8);
t701 = -t1258 * t739 + t1263 * t781;
t700 = t1258 * t781 + t1263 * t739;
t699 = -pkin(1) * t858 - t1251 * t811 + t1253 * t1284;
t697 = -t1251 * t760 + t1253 * t1330;
t696 = t1251 * t1330 + t1253 * t760;
t695 = -pkin(2) * t1031 + pkin(9) * t761 - qJ(4) * t1421 + t1261 * t805;
t692 = t1256 * t757 + t1261 * t755;
t691 = t1256 * t756 + t1261 * t754;
t689 = -t1258 * t717 + t1263 * t767;
t688 = t1258 * t767 + t1263 * t717;
t687 = -pkin(1) * t848 - t1251 * t800 + t1253 * t1285;
t681 = t1256 * t741 + t1261 * t740;
t680 = -t1251 * t750 + t1253 * t1336;
t679 = -t1257 * t888 + t1262 * t694;
t678 = -t1257 * t887 + t1262 * t693;
t672 = -t1251 * t732 + t1253 * t1329;
t670 = -t1251 * t733 + t1253 * t1331;
t669 = t1251 * t1331 + t1253 * t733;
t665 = -t1257 * t824 + t1262 * t682;
t664 = pkin(2) * t1449 + pkin(9) * t823 + t1256 * t759 + t1261 * t742;
t663 = -t1251 * t722 + t1253 * t1337;
t662 = -t1251 * t721 + t1253 * t1338;
t661 = -t1257 * t850 + t1262 * t713 + (-t1251 * t812 - t1253 * t813) * pkin(8);
t657 = -t1258 * t697 + t1263 * t758;
t656 = t1258 * t758 + t1263 * t697;
t655 = -pkin(2) * t953 + pkin(9) * t794 + t1256 * t749 + t1261 * t728;
t648 = -t1258 * t670 + t1263 * t718;
t647 = t1258 * t718 + t1263 * t670;
t643 = -pkin(1) * t812 - t1251 * t711 + t1253 * t1286;
t642 = -t1251 * t692 + t1253 * t1339;
t641 = -t1251 * t691 + t1253 * t1340;
t639 = -t1251 * t685 + t1253 * t1341;
t638 = t1251 * t1341 + t1253 * t685;
t637 = -t1251 * t683 + t1253 * t1342;
t636 = t1251 * t1342 + t1253 * t683;
t635 = -t1251 * t681 + t1253 * t1343;
t630 = -t1251 * t676 + t1253 * t1344;
t629 = t1251 * t1344 + t1253 * t676;
t625 = -t1257 * t712 + t1262 * t668 + (-t1251 * t738 - t1253 * t739) * pkin(8);
t622 = -t1258 * t639 + t1263 * t673;
t621 = t1258 * t673 + t1263 * t639;
t620 = -t1258 * t637 + t1263 * t671;
t619 = t1258 * t671 + t1263 * t637;
t617 = -t1257 * t715 + t1262 * t698 + (-t1251 * t696 - t1253 * t697) * pkin(8);
t616 = -t1257 * t703 + t1262 * t658 + (-t1251 * t716 - t1253 * t717) * pkin(8);
t610 = -t1258 * t630 + t1263 * t660;
t609 = t1258 * t660 + t1263 * t630;
t608 = -pkin(1) * t696 - t1251 * t695 + t1253 * t1287;
t607 = -pkin(1) * t738 - t1251 * t664 + t1253 * t1288;
t604 = -pkin(2) * t1000 + pkin(9) * t735 + t1256 * t646 + t1261 * t645;
t602 = -pkin(1) * t716 - t1251 * t655 + t1253 * t1289;
t600 = -t1251 * t631 + t1253 * t1345;
t599 = t1251 * t1345 + t1253 * t631;
t594 = -t1258 * t600 + t1263 * t627;
t593 = t1258 * t627 + t1263 * t600;
t591 = -pkin(2) * t935 + pkin(9) * t632 + t1256 * t628 + t1261 * t626;
t590 = -t1257 * t690 + t1262 * t605 + (-t1251 * t669 - t1253 * t670) * pkin(8);
t587 = -pkin(2) * t873 + pkin(9) * t686 + t1256 * t615 + t1261 * t613;
t586 = -pkin(2) * t865 + pkin(9) * t684 + t1256 * t614 + t1261 * t612;
t584 = -pkin(2) * t825 + pkin(9) * t677 + t1256 * t603 + t1261 * t601;
t583 = -pkin(1) * t669 - t1251 * t604 + t1253 * t1290;
t579 = -t1251 * t596 + t1253 * t1346;
t578 = t1251 * t1346 + t1253 * t596;
t577 = -t1257 * t634 + t1262 * t589 + (-t1251 * t638 - t1253 * t639) * pkin(8);
t576 = -t1257 * t633 + t1262 * t588 + (-t1251 * t636 - t1253 * t637) * pkin(8);
t575 = -t1257 * t618 + t1262 * t585 + (-t1251 * t629 - t1253 * t630) * pkin(8);
t574 = -t1258 * t579 + t1263 * t595;
t573 = t1258 * t595 + t1263 * t579;
t572 = -t1257 * t606 + t1262 * t592 + (-t1251 * t599 - t1253 * t600) * pkin(8);
t571 = -pkin(1) * t638 - t1251 * t587 + t1253 * t1292;
t570 = -pkin(1) * t636 - t1251 * t586 + t1253 * t1293;
t569 = -pkin(1) * t599 - t1251 * t591 + t1253 * t1291;
t568 = -pkin(1) * t629 - t1251 * t584 + t1253 * t1294;
t566 = -pkin(2) * t674 + pkin(9) * t597 + t1256 * t582 + t1261 * t580;
t565 = -t1257 * t581 + t1262 * t567 + (-t1251 * t578 - t1253 * t579) * pkin(8);
t564 = -pkin(1) * t578 - t1251 * t566 + t1253 * t1295;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1233, -t1234, 0, t1193, 0, 0, 0, 0, 0, 0, t1081, t1076, t1070, t1012, 0, 0, 0, 0, 0, 0, t929, t934, t918, t842, 0, 0, 0, 0, 0, 0, t798, t816, t765, t657, 0, 0, 0, 0, 0, 0, t689, t701, t648, t594, 0, 0, 0, 0, 0, 0, t620, t622, t610, t574; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1234, -t1233, 0, t1192, 0, 0, 0, 0, 0, 0, t1080, t1075, t1069, t1011, 0, 0, 0, 0, 0, 0, t928, t933, t917, t841, 0, 0, 0, 0, 0, 0, t797, t815, t764, t656, 0, 0, 0, 0, 0, 0, t688, t700, t647, t593, 0, 0, 0, 0, 0, 0, t619, t621, t609, t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1141, t1133, t1122, t1073, 0, 0, 0, 0, 0, 0, t994, t998, t961, t881, 0, 0, 0, 0, 0, 0, t848, t858, t812, t696, 0, 0, 0, 0, 0, 0, t716, t738, t669, t599, 0, 0, 0, 0, 0, 0, t636, t638, t629, t578; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1234, 0, -t1233, 0, t1349, -t1217, -t1192, -pkin(7) * t1192, -t1258 * t1163 + t1180 * t1263, -t1258 * t1124 + t1147 * t1263, -t1258 * t1142 + t1176 * t1263, -t1258 * t1162 + t1179 * t1263, -t1258 * t1143 + t1177 * t1263, t1263 * t1190 + t1258 * t1305, -pkin(7) * t1080 - t1258 * t1036 + t1043 * t1263, -pkin(7) * t1075 - t1258 * t1032 + t1038 * t1263, -pkin(7) * t1069 + t1002 * t1263 - t1258 * t1009, -pkin(7) * t1011 + t1005 * t1263 - t1258 * t1013, -t1258 * t1020 + t1078 * t1263, t1033 * t1263 - t1258 * t964, -t1258 * t1006 + t1049 * t1263, -t1258 * t1019 + t1077 * t1263, -t1258 * t1007 + t1050 * t1263, -t1258 * t1071 + t1128 * t1263, -pkin(7) * t928 - t1258 * t840 + t1263 * t846, -pkin(7) * t933 - t1258 * t843 + t1263 * t847, -pkin(7) * t917 - t1258 * t772 + t1263 * t810, -pkin(7) * t841 - t1258 * t762 + t1263 * t763, -t1258 * t886 + t1263 * t945, -t1258 * t821 + t1263 * t891, -t1258 * t867 + t1263 * t930, -t1258 * t885 + t1263 * t944, -t1258 * t868 + t1263 * t931, t1010 * t1263 - t1258 * t932, -pkin(7) * t797 - t1258 * t687 + t1263 * t702, -pkin(7) * t815 - t1258 * t699 + t1263 * t706, -pkin(7) * t764 - t1258 * t643 + t1263 * t661, -pkin(7) * t656 - t1258 * t608 + t1263 * t617, -t1258 * t720 + t1263 * t771, -t1258 * t672 + t1263 * t727, -t1258 * t745 + t1263 * t783, -t1258 * t719 + t1263 * t770, -t1258 * t746 + t1263 * t784, -t1258 * t773 + t1263 * t853, -pkin(7) * t688 - t1258 * t602 + t1263 * t616, -pkin(7) * t700 - t1258 * t607 + t1263 * t625, -pkin(7) * t647 - t1258 * t583 + t1263 * t590, -pkin(7) * t593 - t1258 * t569 + t1263 * t572, -t1258 * t663 + t1263 * t710, -t1258 * t635 + t1263 * t665, -t1258 * t641 + t1263 * t678, -t1258 * t662 + t1263 * t709, -t1258 * t642 + t1263 * t679, -t1258 * t680 + t1263 * t731, -pkin(7) * t619 - t1258 * t570 + t1263 * t576, -pkin(7) * t621 - t1258 * t571 + t1263 * t577, -pkin(7) * t609 - t1258 * t568 + t1263 * t575, -pkin(7) * t573 - t1258 * t564 + t1263 * t565; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1233, 0, t1234, 0, t1217, t1349, t1193, pkin(7) * t1193, t1163 * t1263 + t1258 * t1180, t1124 * t1263 + t1258 * t1147, t1142 * t1263 + t1258 * t1176, t1162 * t1263 + t1258 * t1179, t1143 * t1263 + t1258 * t1177, t1258 * t1190 - t1263 * t1305, pkin(7) * t1081 + t1036 * t1263 + t1258 * t1043, pkin(7) * t1076 + t1032 * t1263 + t1258 * t1038, pkin(7) * t1070 + t1258 * t1002 + t1009 * t1263, pkin(7) * t1012 + t1258 * t1005 + t1013 * t1263, t1020 * t1263 + t1258 * t1078, t1258 * t1033 + t1263 * t964, t1006 * t1263 + t1258 * t1049, t1019 * t1263 + t1258 * t1077, t1007 * t1263 + t1258 * t1050, t1071 * t1263 + t1258 * t1128, pkin(7) * t929 + t1258 * t846 + t1263 * t840, pkin(7) * t934 + t1258 * t847 + t1263 * t843, pkin(7) * t918 + t1258 * t810 + t1263 * t772, pkin(7) * t842 + t1258 * t763 + t1263 * t762, t1258 * t945 + t1263 * t886, t1258 * t891 + t1263 * t821, t1258 * t930 + t1263 * t867, t1258 * t944 + t1263 * t885, t1258 * t931 + t1263 * t868, t1258 * t1010 + t1263 * t932, pkin(7) * t798 + t1258 * t702 + t1263 * t687, pkin(7) * t816 + t1258 * t706 + t1263 * t699, pkin(7) * t765 + t1258 * t661 + t1263 * t643, pkin(7) * t657 + t1258 * t617 + t1263 * t608, t1258 * t771 + t1263 * t720, t1258 * t727 + t1263 * t672, t1258 * t783 + t1263 * t745, t1258 * t770 + t1263 * t719, t1258 * t784 + t1263 * t746, t1258 * t853 + t1263 * t773, pkin(7) * t689 + t1258 * t616 + t1263 * t602, pkin(7) * t701 + t1258 * t625 + t1263 * t607, pkin(7) * t648 + t1258 * t590 + t1263 * t583, pkin(7) * t594 + t1258 * t572 + t1263 * t569, t1258 * t710 + t1263 * t663, t1258 * t665 + t1263 * t635, t1258 * t678 + t1263 * t641, t1258 * t709 + t1263 * t662, t1258 * t679 + t1263 * t642, t1258 * t731 + t1263 * t680, pkin(7) * t620 + t1258 * t576 + t1263 * t570, pkin(7) * t622 + t1258 * t577 + t1263 * t571, pkin(7) * t610 + t1258 * t575 + t1263 * t568, pkin(7) * t574 + t1258 * t565 + t1263 * t564; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1238, t1239, 0, 0, t1215 * t1384 + t1319, t1221 * t1253 + t1251 * t1311, t1186 * t1253 + t1251 * t1307, t1216 * t1383 - t1319, t1187 * t1253 + t1251 * t1308, t1253 * t1352, pkin(1) * t1144 - t1166 * t1253 + t1251 * t1347, pkin(1) * t1134 - t1167 * t1253 + t1251 * t1348, pkin(1) * t1123 + t1251 * t1281, pkin(1) * t1074 + t1103 * t1428, t1131 * t1253 + t1251 * t1296, t1253 * t1059 + t1251 * t1318, t1253 * t1115 + t1251 * t1314, t1129 * t1253 + t1251 * t1297, t1253 * t1116 + t1251 * t1313, t1148 * t1253 + t1149 * t1384 + t1194, pkin(1) * t995 + t1003 * t1253 + t1251 * t1283, pkin(1) * t999 + t1008 * t1253 + t1251 * t1282, pkin(1) * t962 + t1251 * t1278 + t1253 * t921, pkin(1) * t882 + t1251 * t1277 + t1253 * t939, t1251 * t1301 + t1253 * t975, t1251 * t1323 + t1253 * t907, t1251 * t1325 + t1253 * t969, t1251 * t1302 + t1253 * t974, t1251 * t1327 + t1253 * t970, t1014 * t1253 + t1015 * t1384 + t1194, pkin(1) * t849 + t1251 * t1285 + t1253 * t800, pkin(1) * t859 + t1251 * t1284 + t1253 * t811, pkin(1) * t813 + t1251 * t1286 + t1253 * t711, pkin(1) * t697 + t1251 * t1287 + t1253 * t695, t1251 * t1303 + t1253 * t786, t1251 * t1329 + t1253 * t732, t1251 * t1333 + t1253 * t828, t1251 * t1304 + t1253 * t785, t1251 * t1332 + t1253 * t829, t1251 * t1321 + t1253 * t860, pkin(1) * t717 + t1251 * t1289 + t1253 * t655, pkin(1) * t739 + t1251 * t1288 + t1253 * t664, pkin(1) * t670 + t1251 * t1290 + t1253 * t604, pkin(1) * t600 + t1251 * t1291 + t1253 * t591, t1251 * t1337 + t1253 * t722, t1251 * t1343 + t1253 * t681, t1251 * t1340 + t1253 * t691, t1251 * t1338 + t1253 * t721, t1251 * t1339 + t1253 * t692, t1251 * t1336 + t1253 * t750, pkin(1) * t637 + t1251 * t1293 + t1253 * t586, pkin(1) * t639 + t1251 * t1292 + t1253 * t587, pkin(1) * t630 + t1251 * t1294 + t1253 * t584, pkin(1) * t579 + t1251 * t1295 + t1253 * t566;];
tauB_reg  = t1;
