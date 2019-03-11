% Calculate inertial parameters regressor of coriolis matrix for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRRPRR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:51
% EndTime: 2019-03-09 18:09:42
% DurationCPUTime: 44.13s
% Computational Cost: add. (67928->1064), mult. (128691->1363), div. (0->0), fcn. (155727->10), ass. (0->830)
t1360 = cos(qJ(3));
t879 = cos(qJ(2));
t1088 = t1360 * t879;
t876 = sin(qJ(3));
t877 = sin(qJ(2));
t846 = t876 * t877 - t1088;
t873 = sin(pkin(11));
t1220 = t873 * t846;
t1347 = cos(pkin(11));
t1205 = t876 * t879;
t849 = t1360 * t877 + t1205;
t832 = t1347 * t849;
t1414 = t832 - t1220;
t875 = sin(qJ(5));
t1434 = t875 * t1414;
t839 = t849 * qJ(4);
t1408 = -pkin(8) - pkin(7);
t1016 = t1408 * t1205;
t1013 = t1408 * t1360;
t850 = t877 * t1013;
t969 = t850 + t1016;
t1420 = -t839 + t969;
t1222 = t873 * t1420;
t1234 = t846 * qJ(4);
t1094 = t876 * t1408;
t1015 = t877 * t1094;
t993 = t1408 * t1088;
t796 = -t993 + t1015;
t898 = -t796 + t1234;
t688 = t1347 * t898;
t1459 = -t688 + t1222;
t878 = cos(qJ(5));
t1292 = t1459 * t878;
t771 = t1347 * t846 + t849 * t873;
t868 = -t879 * pkin(2) - pkin(1);
t806 = t846 * pkin(3) + t868;
t936 = t771 * pkin(4) + t806;
t907 = -pkin(9) * t1414 + t936;
t312 = t875 * t907 + t1292;
t271 = -pkin(10) * t1434 + t312;
t874 = sin(qJ(6));
t1216 = t874 * t271;
t1359 = cos(qJ(6));
t1244 = t1414 * t878;
t1356 = t771 * pkin(5);
t1471 = t1459 * t875;
t311 = -t878 * t907 + t1471;
t882 = -pkin(10) * t1244 + t1356 - t311;
t243 = t1359 * t882;
t170 = -t243 + t1216;
t270 = -t1471 + ((-pkin(10) - pkin(9)) * t1414 + t936) * t878;
t1085 = t1359 * t270;
t183 = t1085 - t1216;
t1465 = t170 + t183;
t1084 = t1359 * t271;
t881 = t874 * t882;
t171 = t1084 + t881;
t1217 = t874 * t270;
t182 = -t1084 - t1217;
t1464 = t171 + t182;
t1033 = t1347 * t1420;
t895 = t873 * t898;
t511 = t1033 + t895;
t1284 = t511 * t771;
t1014 = t879 * t1094;
t795 = t850 + t1014;
t705 = -t839 + t795;
t1221 = t873 * t705;
t992 = t879 * t1013;
t793 = t992 - t1015;
t897 = t793 + t1234;
t686 = t1347 * t897;
t509 = -t686 + t1221;
t983 = t509 * t1414 + t1284;
t1445 = t771 / 0.2e1;
t861 = t873 * pkin(3) + pkin(9);
t1351 = pkin(10) + t861;
t1038 = t1351 * t875;
t840 = t1351 * t878;
t1028 = t1359 * t1038 + t874 * t840;
t1444 = t1028 / 0.2e1;
t1030 = t1347 * t876;
t1100 = t1360 * pkin(2);
t867 = t1100 + pkin(3);
t827 = pkin(2) * t1030 + t873 * t867;
t816 = pkin(9) + t827;
t1352 = pkin(10) + t816;
t1039 = t1352 * t875;
t792 = t1352 * t878;
t1029 = t1359 * t1039 + t874 * t792;
t1401 = t1029 / 0.2e1;
t1080 = t1359 * t878;
t1213 = t874 * t875;
t955 = t1080 - t1213;
t1433 = t955 * t1414;
t1472 = -t1433 / 0.2e1;
t1081 = t1359 * t875;
t1212 = t874 * t878;
t848 = t1081 + t1212;
t1435 = t848 * t1414;
t1466 = -t1435 / 0.2e1;
t1441 = t848 * t771;
t1232 = t848 * t1441;
t1208 = t875 * t771;
t715 = t874 * t1208;
t537 = -t1080 * t771 + t715;
t1277 = t537 * t955;
t213 = -t1232 + t1277;
t1470 = t213 * qJD(4);
t1294 = t1459 * t1414;
t1108 = qJD(2) + qJD(3);
t1107 = qJD(5) + qJD(6);
t1374 = -t848 / 0.2e1;
t1376 = t955 / 0.2e1;
t1469 = t1464 * t1374 + t1376 * t1465;
t757 = -t1038 * t874 + t1359 * t840;
t1353 = t878 * pkin(5);
t862 = -pkin(3) * t1347 - pkin(4);
t853 = t862 - t1353;
t1468 = t1445 * t757 + t853 * t1472;
t1453 = t1444 * t771 + t1466 * t853;
t603 = -t1039 * t874 + t1359 * t792;
t1358 = pkin(2) * t876;
t860 = t873 * t1358;
t826 = t1347 * t867 - t860;
t815 = -pkin(4) - t826;
t801 = t815 - t1353;
t1467 = t1445 * t603 + t801 * t1472;
t1454 = t1401 * t771 + t1466 * t801;
t1440 = t878 * t771;
t729 = -t1440 / 0.2e1;
t1373 = t848 / 0.2e1;
t309 = -t1373 * t1435 + t1376 * t1433;
t1462 = t1107 * t309;
t1432 = t1414 * qJD(1);
t1461 = t771 * t1432;
t1023 = t1108 * t848;
t1436 = t309 * qJD(1);
t1458 = -t955 * t1023 - t1436;
t1457 = t1433 * qJD(1);
t1456 = t1435 * qJD(1);
t1439 = qJD(1) * t771;
t1017 = -qJD(5) - t1439;
t1026 = t1108 * t955;
t1455 = t848 * t1026 + t1436;
t388 = pkin(5) * t1434 - t511;
t1296 = t388 * t848;
t359 = t1296 / 0.2e1;
t767 = t1414 * pkin(5);
t1004 = pkin(10) * t1440 + t767;
t1211 = t875 * t511;
t1355 = t849 * pkin(3);
t525 = pkin(4) * t1414 + pkin(9) * t771 + t1355;
t514 = t878 * t525;
t319 = t514 - t1211;
t247 = t1004 + t319;
t1086 = t1359 * t247;
t1104 = pkin(10) * t1208;
t501 = t878 * t511;
t513 = t875 * t525;
t320 = t501 + t513;
t277 = t1104 + t320;
t1214 = t874 * t277;
t938 = -t1214 / 0.2e1 + t1086 / 0.2e1;
t79 = t359 + t938 - t1468;
t1032 = t1347 * t705;
t894 = t873 * t897;
t512 = t1032 + t894;
t1210 = t875 * t512;
t870 = t877 * pkin(2);
t507 = t525 + t870;
t496 = t878 * t507;
t317 = t496 - t1210;
t246 = t1004 + t317;
t1087 = t1359 * t246;
t495 = t875 * t507;
t502 = t878 * t512;
t318 = t502 + t495;
t274 = t1104 + t318;
t1215 = t874 * t274;
t939 = -t1215 / 0.2e1 + t1087 / 0.2e1;
t67 = t359 + t939 - t1467;
t1083 = t1359 * t274;
t1009 = -t1083 / 0.2e1;
t1219 = t874 * t246;
t941 = -t1219 / 0.2e1 + t1009;
t1452 = t941 - t1454;
t1451 = t939 + t1467;
t1082 = t1359 * t277;
t1008 = -t1082 / 0.2e1;
t1218 = t874 * t247;
t940 = -t1218 / 0.2e1 + t1008;
t1450 = t940 - t1453;
t1449 = t938 + t1468;
t1448 = t1414 ^ 2;
t1426 = 0.2e1 * t1414;
t1447 = -t603 / 0.2e1;
t1446 = -t757 / 0.2e1;
t1391 = -t1028 / 0.2e1;
t768 = -t1220 / 0.2e1 + t832 / 0.2e1;
t419 = t768 * qJD(5) + t1461;
t1431 = t1108 * t309 - t1435 * t1457;
t1430 = t1107 * t1028;
t1429 = t1107 * t1029;
t1428 = t1107 * t603;
t1427 = t1107 * t757;
t1206 = t875 * t878;
t1425 = t1206 * t1426;
t1424 = t171 / 0.2e1;
t835 = (t1360 * t873 + t1030) * pkin(2);
t1379 = t835 / 0.2e1;
t1019 = t1107 * t955;
t1419 = t1108 * t1414;
t337 = qJD(6) * t768 + t419;
t1101 = t1359 / 0.2e1;
t1007 = t1414 * t1101;
t1105 = pkin(5) * t1244;
t1018 = t1105 / 0.2e1;
t1209 = t875 * t1435;
t1411 = -t955 * t1018 + (t1209 / 0.2e1 + t1007) * pkin(5);
t1092 = -t1347 / 0.2e1;
t1367 = t873 / 0.2e1;
t1409 = (t509 * t1092 + t512 * t1367) * pkin(3);
t1075 = qJD(1) * t1206;
t871 = t875 ^ 2;
t872 = t878 ^ 2;
t856 = t872 - t871;
t599 = -0.2e1 * t1075 * t1414 + t1108 * t856;
t1368 = t871 / 0.2e1;
t554 = (t1368 - t872 / 0.2e1) * t1414;
t290 = t1075 * t1448 + t1108 * t554;
t236 = -t1433 * t848 - t1435 * t955;
t255 = -t1433 ^ 2 + t1435 ^ 2;
t74 = qJD(1) * t255 + t1108 * t236;
t706 = -t848 ^ 2 + t955 ^ 2;
t214 = qJD(1) * t236 + t1108 * t706;
t176 = t1086 - t1214;
t1407 = -t176 / 0.2e1;
t1406 = -t318 / 0.2e1;
t1405 = -t511 / 0.2e1;
t1404 = t511 / 0.2e1;
t1400 = -t1029 / 0.2e1;
t1398 = t603 / 0.2e1;
t1397 = t688 / 0.2e1;
t836 = t1100 * t1347 - t860;
t689 = t848 * t836;
t1396 = t689 / 0.2e1;
t690 = t955 * t836;
t1395 = t690 / 0.2e1;
t1389 = t757 / 0.2e1;
t1388 = t1414 / 0.2e1;
t1387 = -t771 / 0.2e1;
t1384 = t801 / 0.2e1;
t1383 = -t815 / 0.2e1;
t1382 = t815 / 0.2e1;
t1381 = -t816 / 0.2e1;
t1380 = t816 / 0.2e1;
t1378 = -t836 / 0.2e1;
t1377 = t836 / 0.2e1;
t1375 = -t955 / 0.2e1;
t1372 = -t853 / 0.2e1;
t1371 = t853 / 0.2e1;
t1370 = t861 / 0.2e1;
t1369 = -t862 / 0.2e1;
t1366 = -t874 / 0.2e1;
t1365 = t874 / 0.2e1;
t1364 = -t875 / 0.2e1;
t1363 = t875 / 0.2e1;
t1362 = -t878 / 0.2e1;
t1361 = t878 / 0.2e1;
t1357 = pkin(5) * t875;
t1354 = t874 * pkin(5);
t1350 = pkin(3) * qJD(3);
t1349 = pkin(5) * qJD(5);
t1348 = qJD(2) * pkin(2);
t911 = t1384 * t1414 + t1398 * t537 - t1401 * t1441;
t172 = t1087 - t1215;
t173 = t1083 + t1219;
t968 = t1374 * t173 + t1375 * t172;
t40 = t911 + t968;
t1343 = qJD(1) * t40;
t1095 = -t1356 / 0.2e1;
t880 = -t881 / 0.2e1 + t874 * t1095;
t51 = t1217 / 0.2e1 + t880;
t1342 = qJD(1) * t51;
t1298 = t388 * t1433;
t72 = t1105 * t1435 + t182 * t771 + t1298;
t1341 = qJD(1) * t72;
t1299 = t388 * t1435;
t73 = t1105 * t1433 - t183 * t771 - t1299;
t1340 = qJD(1) * t73;
t83 = t170 * t771 - t1299;
t1339 = qJD(1) * t83;
t84 = -t171 * t771 + t1298;
t1338 = qJD(1) * t84;
t1285 = t511 * t1414;
t1316 = t312 * t878;
t1318 = t311 * t875;
t95 = -t1285 - (t1316 + t1318) * t771;
t1337 = qJD(1) * t95;
t1331 = t172 * t848;
t1330 = t173 * t955;
t1329 = t176 * t848;
t177 = t1082 + t1218;
t1328 = t177 * t955;
t989 = t1441 * t171 + t170 * t537;
t18 = -t1433 * t172 - t1435 * t173 + t989;
t1327 = t18 * qJD(1);
t19 = -t1433 * t176 - t1435 * t177 + t989;
t1324 = t19 * qJD(1);
t741 = pkin(5) * t1208;
t387 = t509 - t741;
t20 = -t170 * t172 + t171 * t173 + t387 * t388;
t1323 = t20 * qJD(1);
t21 = -t1433 * t1464 - t1435 * t1465;
t1322 = t21 * qJD(1);
t386 = t1459 - t741;
t24 = -t170 * t176 + t171 * t177 + t386 * t388;
t1321 = t24 * qJD(1);
t29 = t1105 * t388 - t170 * t182 + t171 * t183;
t1320 = t29 * qJD(1);
t1319 = t311 * t1414;
t1317 = t312 * t1414;
t1315 = t317 * t875;
t1314 = t317 * t878;
t1313 = t318 * t875;
t1312 = t318 * t878;
t1311 = t319 * t875;
t1310 = t319 * t878;
t1309 = t320 * t875;
t1308 = t320 * t878;
t988 = -t1414 * t170 - t1441 * t388;
t35 = t1435 * t387 + t172 * t771 + t988;
t1307 = t35 * qJD(1);
t987 = -t1414 * t171 + t388 * t537;
t36 = t1433 * t387 - t173 * t771 + t987;
t1306 = t36 * qJD(1);
t37 = t1435 * t386 + t176 * t771 + t988;
t1305 = t37 * qJD(1);
t38 = t1433 * t386 - t177 * t771 + t987;
t1304 = t38 * qJD(1);
t1303 = t386 * t955;
t1302 = t386 * t848;
t1301 = t387 * t955;
t1300 = t387 * t848;
t1297 = t388 * t955;
t41 = t1414 * t388 - t1441 * t170 + t171 * t537;
t1295 = t41 * qJD(1);
t1290 = t509 * t862;
t1289 = t509 * t875;
t1288 = t509 * t878;
t1286 = t511 * t509;
t947 = -t243 / 0.2e1 + t1359 * t1095;
t53 = t1085 / 0.2e1 + t947;
t1283 = t53 * qJD(1);
t1282 = t1441 * t1433;
t1281 = t1441 * t771;
t1280 = t1435 * t537;
t1279 = t1435 * t1414;
t1278 = t537 * t771;
t1276 = t1433 * t1414;
t59 = -t311 * t317 + t312 * t318 - t1286;
t1275 = t59 * qJD(1);
t60 = -t511 * t1459 - t311 * t319 + t312 * t320;
t1274 = t60 * qJD(1);
t1271 = t1029 * t537;
t1270 = t1029 * t1414;
t1267 = t603 * t1441;
t1265 = t603 * t1414;
t946 = (-t311 * t878 + t312 * t875) * t771;
t61 = (t1313 + t1314) * t1414 - t946;
t1263 = t61 * qJD(1);
t63 = (t1309 + t1310) * t1414 - t946;
t1262 = t63 * qJD(1);
t1259 = t1028 * t537;
t1258 = t1028 * t1414;
t1255 = t757 * t1441;
t1253 = t757 * t1414;
t1250 = t1414 * t816;
t1249 = t1414 * t861;
t1248 = t1414 * t873;
t1247 = t771 ^ 2;
t1246 = t771 * t815;
t1243 = t801 * t1441;
t1241 = t801 * t537;
t1239 = t801 * t955;
t1238 = t801 * t848;
t1237 = t826 * t771;
t1236 = t827 * t1414;
t1235 = t955 * t771;
t1233 = t846 * t849;
t1230 = t853 * t1441;
t1228 = t853 * t537;
t1226 = t853 * t955;
t1225 = t853 * t848;
t1224 = t862 * t875;
t87 = t317 * t771 + t875 * t983 - t1319;
t1223 = t87 * qJD(1);
t1207 = t875 * t836;
t88 = -t318 * t771 + t878 * t983 - t1317;
t1202 = t88 * qJD(1);
t984 = t1284 + t1294;
t89 = t319 * t771 + t875 * t984 - t1319;
t1201 = t89 * qJD(1);
t90 = -t320 * t771 + t878 * t984 - t1317;
t1200 = t90 * qJD(1);
t1043 = t1389 + t1446;
t1044 = t1391 + t1444;
t1045 = t1398 + t1447;
t1047 = t1400 + t1401;
t149 = (-t1043 - t1045) * t848 - (-t1044 - t1047) * t955;
t1197 = t149 * qJD(5);
t751 = t871 * t771;
t752 = t872 * t771;
t577 = t751 + t752;
t1196 = t577 * qJD(4);
t1194 = t1107 * t236;
t287 = t1043 * t955 + t1044 * t848;
t1193 = t287 * qJD(5);
t937 = -t1212 / 0.2e1 - t1081 / 0.2e1;
t920 = t937 * t771;
t380 = t1441 / 0.2e1 - t920;
t1190 = t1107 * t380;
t385 = -t1441 / 0.2e1 - t920;
t1189 = t1107 * t385;
t1006 = -t1080 / 0.2e1;
t1184 = -t715 / 0.2e1 - t771 * t1006;
t1005 = t1080 / 0.2e1;
t1183 = t715 / 0.2e1 - t771 * t1005;
t1148 = qJD(3) * t878;
t1151 = qJD(2) * t878;
t1182 = (t1148 + t1151) * t1414;
t1181 = t1239 / 0.2e1 + t1226 / 0.2e1;
t1180 = t836 * t1006 + t1207 * t1365;
t1179 = qJD(1) * t213;
t220 = t1434 * t511 + t311 * t771;
t1178 = qJD(1) * t220;
t221 = -t1244 * t511 - t312 * t771;
t1177 = qJD(1) * t221;
t222 = -t1280 - t1282;
t1176 = qJD(1) * t222;
t223 = -t1280 + t1282;
t1175 = qJD(1) * t223;
t249 = -t1459 * t771 - t1285;
t1173 = qJD(1) * t249;
t265 = t1279 + t1281;
t1171 = qJD(1) * t265;
t267 = t1276 - t1278;
t1170 = qJD(1) * t267;
t982 = -t1247 + t1448;
t350 = t982 * t875;
t1169 = qJD(1) * t350;
t364 = t1247 + t1448;
t351 = t364 * t875;
t1168 = qJD(1) * t351;
t353 = t364 * t878;
t1167 = qJD(1) * t353;
t623 = t806 * t1414;
t810 = t870 + t1355;
t415 = t771 * t810 + t623;
t1166 = qJD(1) * t415;
t624 = t806 * t771;
t416 = t1414 * t810 - t624;
t1165 = qJD(1) * t416;
t433 = -t1355 * t771 - t623;
t1164 = qJD(1) * t433;
t434 = -t1355 * t1414 + t624;
t1163 = qJD(1) * t434;
t749 = t846 * t870 + t849 * t868;
t1160 = qJD(1) * t749;
t750 = -t846 * t868 + t849 * t870;
t1159 = qJD(1) * t750;
t1156 = qJD(1) * t868;
t1155 = qJD(1) * t879;
t1154 = qJD(2) * t955;
t1153 = qJD(2) * t848;
t1152 = qJD(2) * t875;
t1150 = qJD(3) * t853;
t1149 = qJD(3) * t868;
t1147 = qJD(5) * t875;
t1146 = qJD(5) * t878;
t1145 = qJD(6) * t801;
t1144 = qJD(6) * t853;
t1002 = t771 * (t872 / 0.2e1 + t1368);
t914 = -t1002 * t816 + t1382 * t1414;
t966 = -t1314 / 0.2e1 - t1313 / 0.2e1;
t107 = t914 + t966;
t1143 = t107 * qJD(1);
t913 = -t1002 * t861 + t1388 * t862;
t965 = -t1310 / 0.2e1 - t1309 / 0.2e1;
t151 = t913 + t965;
t1142 = t151 * qJD(1);
t163 = -t512 * t771 - t1294 + t983;
t1140 = t163 * qJD(1);
t1049 = t1459 / 0.2e1 - t509 / 0.2e1;
t957 = t1378 * t771 + t1379 * t1414;
t891 = -(t1382 + t1369) * t771 + (t1381 + t1370) * t1414 + t957;
t166 = -t1049 * t878 + t875 * t891;
t1139 = t166 * qJD(1);
t187 = t1459 * t512 + t806 * t810 - t1286;
t1138 = t187 * qJD(1);
t192 = t1355 * t806;
t1137 = t192 * qJD(1);
t900 = t1237 / 0.2e1 - t1236 / 0.2e1 + t957;
t1031 = t1347 * t771;
t921 = (-t1248 / 0.2e1 + t1031 / 0.2e1) * pkin(3);
t234 = t921 - t900;
t1136 = t234 * qJD(1);
t264 = -t1279 + t1281;
t1135 = t264 * qJD(1);
t266 = t1276 + t1278;
t1134 = t266 * qJD(1);
t332 = (-t793 - t796) * t849 + (t969 - t795) * t846;
t1131 = t332 * qJD(1);
t352 = t982 * t878;
t1130 = t352 * qJD(1);
t1129 = t364 * qJD(1);
t1128 = t982 * qJD(1);
t370 = t380 * qJD(1);
t1057 = -t1235 / 0.2e1;
t381 = t1057 + t1183;
t371 = t381 * qJD(1);
t1058 = t1235 / 0.2e1;
t382 = t1058 + t1184;
t373 = t382 * qJD(1);
t834 = t1355 / 0.2e1;
t1106 = t834 + t870 / 0.2e1;
t959 = t1388 * t826 + t1445 * t827;
t389 = t959 + t1106;
t1127 = t389 * qJD(1);
t410 = t793 * t969 + t795 * t796 + t868 * t870;
t1126 = t410 * qJD(1);
t935 = t1092 * t1414 - t1367 * t771;
t493 = (-t849 / 0.2e1 + t935) * pkin(3);
t1123 = t493 * qJD(1);
t1122 = t554 * qJD(1);
t1121 = t554 * qJD(5);
t1120 = t1434 * qJD(1);
t1042 = 0.2e1 * t1445;
t557 = t1042 * t875;
t549 = t557 * qJD(1);
t560 = t1042 * t878;
t550 = t560 * qJD(1);
t561 = 0.2e1 * t729;
t1119 = t561 * qJD(1);
t1118 = t577 * qJD(1);
t578 = t856 * t1448;
t1117 = t578 * qJD(1);
t707 = t846 ^ 2 - t849 ^ 2;
t1116 = t707 * qJD(1);
t1115 = t768 * qJD(1);
t762 = t1414 * qJD(4);
t828 = t835 * qJD(3);
t857 = -t877 ^ 2 + t879 ^ 2;
t1111 = t857 * qJD(1);
t1110 = t877 * qJD(2);
t1109 = t879 * qJD(2);
t1103 = t848 * t1357;
t1102 = -t1359 / 0.2e1;
t1099 = pkin(1) * t877 * qJD(1);
t1098 = pkin(1) * t1155;
t1097 = t1357 / 0.2e1;
t1096 = -t767 / 0.2e1;
t227 = t1045 * t955 + t1047 * t848;
t411 = t1373 * t690 - t1376 * t689;
t1091 = t411 * qJD(3) + t227 * qJD(5);
t358 = t1297 / 0.2e1;
t1090 = t848 * t1018 + t1097 * t1433 + t358;
t1089 = (-t1359 * t955 - t848 * t874) * t1349;
t1077 = t846 * t1156;
t1076 = t849 * t1156;
t1074 = t878 * t762;
t1073 = qJD(5) * t771 * t1414;
t589 = t1414 * t1439;
t1071 = qJD(1) * t1233;
t858 = t875 * t1146;
t1070 = t877 * t1109;
t1069 = t1312 / 0.2e1;
t1068 = -t1297 / 0.2e1;
t1067 = -t1296 / 0.2e1;
t1066 = t388 * t1364;
t1065 = t511 * t1379;
t1064 = t511 * t1363;
t1063 = t511 * t1361;
t1061 = -t1244 / 0.2e1;
t1054 = t1224 / 0.2e1;
t1053 = t861 * t1364;
t1052 = t386 / 0.2e1 - t387 / 0.2e1;
t1051 = -t495 / 0.2e1 - t502 / 0.2e1;
t1050 = -t501 / 0.2e1 - t513 / 0.2e1;
t1048 = t1401 + t1391;
t1046 = t1447 + t1389;
t1041 = t1384 + t1372;
t1040 = t1371 + t1384;
t1037 = t1360 * qJD(2);
t1036 = t1360 * qJD(3);
t1035 = t1359 * qJD(5);
t1034 = t1359 * qJD(6);
t740 = (t871 + t872) * t836;
t781 = t1108 * t849;
t1025 = t1108 * t878;
t580 = t1108 * t771;
t1022 = t1108 * t875;
t1021 = t1107 * t771;
t1020 = t1107 * t848;
t1012 = t1448 * t858;
t1001 = t1378 + t1369 + t1383;
t999 = t875 * t1025;
t361 = t771 * t1419;
t998 = t1414 * t580;
t996 = t878 * t1022;
t994 = -qJD(6) + t1017;
t412 = t689 * t848 + t690 * t955;
t962 = -t1395 * t1435 + t1396 * t1433;
t8 = (t1407 + t172 / 0.2e1) * t848 - (-t177 / 0.2e1 + t173 / 0.2e1) * t955 + t1048 * t537 - t1046 * t1441 + t962;
t991 = -qJD(1) * t8 - qJD(2) * t412;
t986 = t1312 - t1315;
t985 = t1308 - t1311;
t981 = -t1246 - t1250;
t980 = -t771 * t862 - t1249;
t963 = t1382 * t1459 - t1065;
t34 = -t1290 / 0.2e1 + (t1377 * t312 + t1380 * t320 + t1406 * t861) * t878 + (t1370 * t317 + t1377 * t311 + t1381 * t319) * t875 + t963;
t414 = t740 * t816 + t815 * t835;
t979 = -t34 * qJD(1) - t414 * qJD(2);
t910 = t1371 * t1414 + t1389 * t537 - t1441 * t1444;
t967 = t1374 * t177 + t1375 * t176;
t43 = t910 + t967;
t978 = -qJD(1) * t43 + qJD(2) * t411;
t893 = (t1061 * t848 + t1364 * t1433 + t1366 * t1414) * pkin(5) + t1068;
t45 = t893 + t1452;
t627 = t1103 + t1239;
t977 = qJD(1) * t45 - qJD(2) * t627;
t892 = (t1007 - t1209 / 0.2e1 - t955 * t1061) * pkin(5) + t1067;
t46 = t892 + t1451;
t830 = t955 * t1357;
t626 = -t830 + t1238;
t976 = qJD(1) * t46 - qJD(2) * t626;
t590 = -t826 * t835 + t827 * t836;
t883 = t827 * t1404 - t1065 + (-t826 / 0.2e1 + t1377) * t1459;
t91 = -t883 + t1409;
t975 = t91 * qJD(1) - t590 * qJD(2);
t81 = (t320 / 0.2e1 + t1406) * t878 + (-t319 / 0.2e1 + t317 / 0.2e1) * t875;
t974 = -qJD(1) * t81 - qJD(2) * t740;
t973 = t1103 + t1181;
t972 = t1414 * t1017;
t321 = t1397 - t686 / 0.2e1 + (-t1420 / 0.2e1 + t705 / 0.2e1) * t873;
t971 = -qJD(1) * t321 + qJD(2) * t835;
t887 = -t894 / 0.2e1 - t1032 / 0.2e1;
t888 = -t895 / 0.2e1 - t1033 / 0.2e1;
t323 = -t887 + t888;
t970 = -qJD(1) * t323 + qJD(2) * t836;
t964 = t1308 / 0.2e1 - t1311 / 0.2e1;
t961 = t1379 * t1435 - t1445 * t689;
t960 = t1379 * t1433 + t1387 * t690;
t958 = t1380 * t771 + t1383 * t1414;
t956 = t1369 * t1414 + t1370 * t771;
t64 = t1067 + t1451;
t954 = qJD(1) * t64 - t1153 * t801;
t65 = t1068 + t1452;
t953 = qJD(1) * t65 - t1154 * t801;
t69 = -t1041 * t1441 - t1048 * t1414 - t1052 * t955 + t961;
t952 = -qJD(1) * t69 + t1154 * t835;
t71 = t1041 * t537 + t1046 * t1414 + t1052 * t848 + t960;
t951 = -qJD(1) * t71 - t1153 * t835;
t169 = t1049 * t875 + t878 * t891;
t950 = -qJD(1) * t169 - t1152 * t835;
t905 = t875 * t958 - t1063;
t193 = t905 - t1051;
t949 = -qJD(1) * t193 - t1151 * t815;
t934 = t958 * t878;
t195 = -t496 / 0.2e1 - t934 + (t1405 + t512 / 0.2e1) * t875;
t948 = -qJD(1) * t195 - t1152 * t815;
t945 = t999 * t1426;
t944 = t1101 * t172 + t1365 * t173;
t943 = t1101 * t176 + t1365 * t177;
t942 = -t1101 * t689 + t1365 * t690;
t933 = t956 * t878;
t922 = (t1102 * t537 - t1366 * t1441) * pkin(5);
t890 = t922 - t1469;
t9 = t1045 * t1433 + t1047 * t1435 + t890;
t932 = qJD(1) * t9 - qJD(3) * t149;
t889 = t1029 * t1407 + t1379 * t388 + t1384 * t386 + t1395 * t171 + t1396 * t170 + t1398 * t177;
t912 = t1372 * t387 + t1444 * t172 + t1446 * t173;
t2 = t889 + t912;
t248 = t1029 * t689 + t603 * t690 + t801 * t835;
t931 = t2 * qJD(1) + t248 * qJD(2) + t411 * qJD(4);
t237 = t1357 * t801;
t903 = t1029 * t1424 + t182 * t1401 + t1447 * t1465;
t3 = (t1061 * t801 + t1066 + t944) * pkin(5) + t903;
t930 = -t3 * qJD(1) + t237 * qJD(2) + t227 * qJD(4);
t12 = t1043 * t1433 + t1044 * t1435 + t890;
t929 = qJD(1) * t12 - qJD(2) * t149;
t923 = (-t1102 * t1441 + t1365 * t537) * pkin(5);
t25 = (-t183 / 0.2e1 - t170 / 0.2e1) * t848 - (t1424 + t182 / 0.2e1) * t955 + t923;
t928 = -t25 * qJD(1) + t227 * qJD(2) + t287 * qJD(3);
t919 = t937 * t836;
t424 = -t1040 * t848 + t919;
t393 = t830 + t424;
t56 = t892 + t1449;
t696 = -t830 + t1225;
t927 = qJD(1) * t56 + qJD(2) * t393 - qJD(3) * t696;
t395 = (t1005 - t1213 / 0.2e1) * t836 + t973;
t55 = t893 + t1450;
t697 = t1103 + t1226;
t926 = qJD(1) * t55 - qJD(2) * t395 - qJD(3) * t697;
t76 = t1067 + t1449;
t918 = qJD(1) * t76 + qJD(2) * t424 - t1150 * t848;
t425 = -t1040 * t955 + t1180;
t77 = t1068 + t1450;
t917 = qJD(1) * t77 + qJD(2) * t425 - t1150 * t955;
t904 = t875 * t956 - t1063;
t216 = t904 - t1050;
t638 = t1001 * t878;
t916 = -qJD(1) * t216 + qJD(2) * t638 - t1148 * t862;
t218 = -t514 / 0.2e1 - t933 + (t1405 + t1404) * t875;
t637 = t1001 * t875;
t915 = -qJD(1) * t218 + qJD(2) * t637 - qJD(3) * t1224;
t120 = (-t1040 * t875 + t942) * pkin(5);
t327 = t1357 * t853;
t902 = t1444 * t1464 + t1446 * t1465;
t5 = (t1061 * t853 + t1066 + t943) * pkin(5) + t902;
t906 = -t5 * qJD(1) - t120 * qJD(2) + t327 * qJD(3) + t287 * qJD(4);
t427 = t1225 / 0.2e1 + t1238 / 0.2e1 + t919;
t899 = -t1250 / 0.2e1 - t1246 / 0.2e1 + t957;
t896 = t922 + t1469;
t859 = t877 * t1155;
t855 = t856 * qJD(5);
t829 = t836 * qJD(3);
t809 = t875 * t828;
t783 = t848 * t828;
t782 = t955 * t828;
t780 = t1108 * t846;
t765 = t771 * qJD(4);
t732 = t878 * t1432;
t710 = t992 / 0.2e1 + t993 / 0.2e1 - t1015;
t709 = -t850 - t1014 / 0.2e1 - t1016 / 0.2e1;
t708 = t740 * qJD(3);
t640 = t836 * t1362 + (t815 + t862) * t1361;
t639 = t1054 + t815 * t1363 - t1207 / 0.2e1;
t618 = t955 * t1020;
t617 = t848 * t1019;
t575 = t1108 * t768;
t562 = (t1387 + t1445) * t878;
t559 = t1440 / 0.2e1 + t729;
t558 = t771 * t1364 + t1208 / 0.2e1;
t548 = t558 * qJD(5);
t547 = t557 * qJD(5);
t544 = t1434 * qJD(4);
t531 = t550 + t1146;
t530 = -t549 - t1147;
t515 = t1107 * t706;
t494 = pkin(3) * t935 + t834;
t492 = t996 - t1122;
t491 = -t999 + t1122;
t448 = t1435 * qJD(4);
t444 = t1433 * qJD(4);
t426 = t1180 + t1181;
t409 = t412 * qJD(3);
t397 = t1017 * t1425;
t396 = t973 + t1180;
t394 = -t830 + t427;
t390 = -t959 + t1106;
t384 = t1057 + t1184;
t383 = t1058 + t1183;
t379 = t385 * qJD(4);
t374 = t384 * qJD(4);
t372 = t382 * qJD(4);
t369 = t380 * qJD(4);
t367 = t388 * t1097;
t335 = t589 * t872 - t1121;
t334 = t589 * t871 + t1121;
t330 = -t1020 - t370;
t329 = -t1019 - t373;
t328 = -t371 + t1019;
t324 = t887 + t888;
t322 = -t1222 / 0.2e1 + t1397 - t1221 / 0.2e1 + t686 / 0.2e1;
t276 = qJD(5) * t560 - t1130;
t275 = -t547 + t1169;
t257 = -t1121 - (t1432 * t872 + t996) * t771;
t256 = t1121 - (t1432 * t871 - t999) * t771;
t238 = (-qJD(5) + t1439) * t1425 + t1108 * (t751 - t752);
t235 = t921 + t900;
t229 = t559 * qJD(5) + t1022 * t1414 + t1130;
t228 = t548 - t1169 + t1182;
t219 = -t1064 - t1211 / 0.2e1 + t514 / 0.2e1 - t933;
t217 = t904 + t1050;
t210 = t1108 * t380 + t1433 * t1439;
t209 = t1108 * t381 + t1435 * t1439;
t196 = -t1064 - t1210 / 0.2e1 + t496 / 0.2e1 - t934;
t194 = t905 + t1051;
t168 = t1471 / 0.2e1 + t1249 * t1362 + t862 * t729 + t1289 / 0.2e1 + t899 * t878;
t167 = -t1292 / 0.2e1 + t1414 * t1053 - t771 * t1054 - t1288 / 0.2e1 + t899 * t875;
t165 = -t1107 * t381 - t1134;
t164 = -t1135 - t1190;
t150 = t913 - t965;
t121 = pkin(5) * t942 + (t801 + t853) * t1097;
t117 = t1108 * t385 + t1433 * t994;
t116 = t1108 * t383 + t1435 * t994;
t115 = -t1457 * t537 + t1462;
t110 = t1441 * t1456 - t1462;
t106 = t914 - t966;
t94 = t1023 * t1414 + t1107 * t383 + t1134;
t93 = t1026 * t1414 + t1135 + t1189;
t92 = t1409 + t883;
t86 = (t1023 + t1457) * t537 + t1462;
t85 = -(t1456 - t1026) * t1441 - t1462;
t80 = t1069 - t1315 / 0.2e1 + t964;
t78 = t358 + t940 + t1453;
t70 = -t1265 / 0.2e1 + t1241 / 0.2e1 + t1302 / 0.2e1 - t1253 / 0.2e1 + t1228 / 0.2e1 + t1300 / 0.2e1 + t960;
t68 = -t1270 / 0.2e1 - t1243 / 0.2e1 - t1303 / 0.2e1 - t1258 / 0.2e1 - t1230 / 0.2e1 - t1301 / 0.2e1 + t961;
t66 = t358 + t941 + t1454;
t62 = -t1175 + t1194;
t58 = t1008 + (t1096 - t247 / 0.2e1) * t874 + t1090 + t1453;
t57 = t1411 + t79;
t54 = t1216 - t1085 / 0.2e1 + t947;
t52 = -t1084 - t1217 / 0.2e1 + t880;
t48 = t1009 + (t1096 - t246 / 0.2e1) * t874 + t1090 + t1454;
t47 = t1411 + t67;
t44 = t1175 + t1108 * (t1232 + t1277) + t1194;
t42 = t910 - t967;
t39 = t911 - t968;
t33 = t861 * t1069 + t317 * t1053 + t1290 / 0.2e1 + (t1316 / 0.2e1 + t1318 / 0.2e1) * t836 + t964 * t816 + t963;
t26 = t1373 * t1465 + t1376 * t1464 + t923;
t13 = t1028 * t1466 - t1391 * t1435 + t896;
t10 = t1029 * t1466 - t1400 * t1435 + t896;
t7 = t1267 / 0.2e1 + t1328 / 0.2e1 + t1271 / 0.2e1 - t1329 / 0.2e1 + t1255 / 0.2e1 + t1330 / 0.2e1 + t1259 / 0.2e1 - t1331 / 0.2e1 + t962;
t6 = pkin(5) * t943 + t1018 * t853 + t367 - t902;
t4 = pkin(5) * t944 + t1018 * t801 + t367 - t903;
t1 = t889 - t912;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1070, t857 * qJD(2), 0, -t1070, 0, 0, -pkin(1) * t1110, -pkin(1) * t1109, 0, 0, -t846 * t781, t1108 * t707, 0, t1108 * t1233, 0, 0, qJD(2) * t749 + t1149 * t849, qJD(2) * t750 - t1149 * t846, qJD(2) * t332, qJD(2) * t410, -t998, -t1108 * t982, 0, t361, 0, 0, qJD(2) * t415 - qJD(3) * t433, qJD(2) * t416 - qJD(3) * t434, qJD(2) * t163 + qJD(4) * t364, qJD(2) * t187 + qJD(3) * t192 + qJD(4) * t249, -t872 * t998 - t1012, -t578 * qJD(5) + t771 * t945, -t1073 * t875 + t1108 * t352, -t871 * t998 + t1012, -t1073 * t878 - t1108 * t350, t361, qJD(2) * t87 + qJD(3) * t89 + qJD(4) * t351 + qJD(5) * t221, qJD(2) * t88 + qJD(3) * t90 + qJD(4) * t353 + qJD(5) * t220, -qJD(2) * t61 - qJD(3) * t63, qJD(2) * t59 + qJD(3) * t60 + qJD(4) * t95 (-t1107 * t1435 + t1108 * t537) * t1433, t1107 * t255 + t1108 * t223, -t1021 * t1435 + t1108 * t266 (t1107 * t1433 - t1108 * t1441) * t1435, -t1021 * t1433 + t1108 * t264, t361, qJD(2) * t35 + qJD(3) * t37 + qJD(4) * t265 + qJD(5) * t72 + qJD(6) * t84, qJD(2) * t36 + qJD(3) * t38 + qJD(4) * t267 + qJD(5) * t73 + qJD(6) * t83, qJD(2) * t18 + qJD(3) * t19 + qJD(4) * t222 + qJD(5) * t21, qJD(2) * t20 + qJD(3) * t24 + qJD(4) * t41 + qJD(5) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t859, t1111, t1109, -t859, -t1110, 0, -pkin(7) * t1109 - t1099, pkin(7) * t1110 - t1098, 0, 0, -t1071, t1116, -t780, t1071, -t781, 0, qJD(2) * t793 + qJD(3) * t710 + t1160, -qJD(2) * t795 + qJD(3) * t709 + t1159, t1131 + (t1360 * t846 - t849 * t876) * t1348, t1126 + (t1360 * t793 + t795 * t876) * t1348, -t589, -t1128, -t580, t1461, -t1419, 0, -qJD(2) * t509 + qJD(3) * t322 + t1166, -qJD(2) * t512 + qJD(3) * t324 + t1165, t1140 + (-t1236 + t1237) * qJD(2) + t235 * qJD(3), t1138 + (-t509 * t826 + t512 * t827) * qJD(2) + t92 * qJD(3) + t390 * qJD(4), t257, t238, t229, t256, t228, t419, t1223 + (t875 * t981 - t1288) * qJD(2) + t167 * qJD(3) + t196 * qJD(5), t1202 + (t878 * t981 + t1289) * qJD(2) + t168 * qJD(3) + t194 * qJD(5), qJD(2) * t986 + t80 * qJD(3) - t1263, t1275 + (t509 * t815 + t816 * t986) * qJD(2) + t33 * qJD(3) + t106 * qJD(4), t86, t44, t94, t85, t93, t337, t1307 + (-t1243 - t1270 - t1301) * qJD(2) + t68 * qJD(3) + t47 * qJD(5) + t67 * qJD(6), t1306 + (t1241 - t1265 + t1300) * qJD(2) + t70 * qJD(3) + t48 * qJD(5) + t66 * qJD(6), t1327 + (t1267 + t1271 + t1330 - t1331) * qJD(2) + t7 * qJD(3) + t10 * qJD(5), t1323 + (-t1029 * t172 + t173 * t603 + t387 * t801) * qJD(2) + t1 * qJD(3) + t39 * qJD(4) + t4 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1071, t1116, -t780, t1071, -t781, 0, qJD(2) * t710 - qJD(3) * t796 + t1076, qJD(2) * t709 - qJD(3) * t969 - t1077, 0, 0, -t589, -t1128, -t580, t1461, -t1419, 0, qJD(2) * t322 - qJD(3) * t1459 - t1164, qJD(2) * t324 - qJD(3) * t511 - t1163, t235 * qJD(2) + (t1031 - t1248) * t1350, t1137 + t92 * qJD(2) + (-t1347 * t1459 + t511 * t873) * t1350 + t494 * qJD(4), t257, t238, t229, t256, t228, t419, t1201 + t167 * qJD(2) + (t875 * t980 - t1292) * qJD(3) + t219 * qJD(5), t1200 + t168 * qJD(2) + (t878 * t980 + t1471) * qJD(3) + t217 * qJD(5), t80 * qJD(2) + qJD(3) * t985 - t1262, t1274 + t33 * qJD(2) + (t1459 * t862 + t861 * t985) * qJD(3) + t150 * qJD(4), t86, t44, t94, t85, t93, t337, t1305 + t68 * qJD(2) + (-t1230 - t1258 - t1303) * qJD(3) + t57 * qJD(5) + t79 * qJD(6), t1304 + t70 * qJD(2) + (t1228 - t1253 + t1302) * qJD(3) + t58 * qJD(5) + t78 * qJD(6), t1324 + t7 * qJD(2) + (t1255 + t1259 + t1328 - t1329) * qJD(3) + t13 * qJD(5), t1321 + t1 * qJD(2) + (-t1028 * t176 + t177 * t757 + t386 * t853) * qJD(3) + t42 * qJD(4) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1129, qJD(2) * t390 + qJD(3) * t494 + t1173, 0, 0, 0, 0, 0, 0, t548 + t1168, qJD(5) * t562 + t1167, 0, qJD(2) * t106 + qJD(3) * t150 + t1337, 0, 0, 0, 0, 0, 0, t1171 + t1189, t1107 * t384 + t1170, t1176, t1295 + t39 * qJD(2) + t42 * qJD(3) + (t1441 * t955 + t537 * t848) * qJD(4) + t26 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290, -t945 - t1117, t1108 * t559 + t875 * t972, t290, t1108 * t558 + t878 * t972, t575, qJD(2) * t196 + qJD(3) * t219 + qJD(4) * t558 - qJD(5) * t312 + t1177, qJD(2) * t194 + qJD(3) * t217 + qJD(4) * t562 + qJD(5) * t311 + t1178, 0, 0, t1431, t74, t116, -t1431, t117, t575, qJD(2) * t47 + qJD(3) * t57 + qJD(5) * t182 + qJD(6) * t52 + t1341 + t379, qJD(2) * t48 + qJD(3) * t58 - qJD(5) * t183 + qJD(6) * t54 + t1340 + t374, t1322 + t10 * qJD(2) + t13 * qJD(3) + (t1359 * t1435 - t1433 * t874) * t1349, t1320 + t4 * qJD(2) + t6 * qJD(3) + t26 * qJD(4) + (t1359 * t182 + t183 * t874) * t1349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1431, t74, t116, -t1431, t117, t575, qJD(2) * t67 + qJD(3) * t79 + qJD(5) * t52 - qJD(6) * t171 + t1338 + t379, qJD(2) * t66 + qJD(3) * t78 + qJD(5) * t54 + qJD(6) * t170 + t1339 + t374, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t859, -t1111, 0, t859, 0, 0, t1099, t1098, 0, 0, t1071, -t1116, 0, -t1071, 0, 0, -t1160, -t1159, -t1131, -t1126, t589, t1128, 0, -t1461, 0, 0, qJD(3) * t321 - t1166 - t762, qJD(3) * t323 - t1165 + t765, -qJD(3) * t234 - t1140, -qJD(3) * t91 - qJD(4) * t389 - t1138, t335, t397, t276, t334, t275, -t419, qJD(3) * t166 + qJD(5) * t195 - t1074 - t1223, qJD(3) * t169 + qJD(5) * t193 - t1202 + t544, qJD(3) * t81 - t1196 + t1263, qJD(3) * t34 + qJD(4) * t107 - t1275, t115, t62, t165, t110, t164, -t337, qJD(3) * t69 - qJD(5) * t46 - qJD(6) * t64 - t1307 - t444, qJD(3) * t71 - qJD(5) * t45 - qJD(6) * t65 - t1306 + t448, qJD(3) * t8 - qJD(5) * t9 - t1327 + t1470, qJD(3) * t2 + qJD(4) * t40 - qJD(5) * t3 - t1323; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t1358, -pkin(2) * t1036, 0, 0, 0, 0, 0, 0, 0, 0, -t828, -t829, 0, t590 * qJD(3), t858, t855, 0, -t858, 0, 0, t1147 * t815 - t828 * t878, t1146 * t815 + t809, t708, t414 * qJD(3), t618, t515, 0, -t617, 0, 0, qJD(5) * t626 + t1145 * t848 - t782, qJD(5) * t627 + t1145 * t955 + t783, t409, qJD(3) * t248 + qJD(5) * t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1108 * t1358 (-t1037 - t1036) * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t828 - t971, -t829 - t970, -t1136 (-t1347 * t835 + t836 * t873) * t1350 - t975, t858, t855, 0, -t858, 0, 0, t639 * qJD(5) - t1025 * t835 + t1139, qJD(5) * t640 + t809 - t950, t708 - t974 (t740 * t861 + t835 * t862) * qJD(3) - t979, t618, t515, 0, -t617, 0, 0, qJD(5) * t394 + qJD(6) * t427 - t782 - t952, qJD(5) * t396 + qJD(6) * t426 + t783 - t951, t409 - t991 + t1197 (t1028 * t689 + t690 * t757 + t835 * t853) * qJD(3) + t121 * qJD(5) + t931; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1432, t1439, 0, -t1127, 0, 0, 0, 0, 0, 0, -t732, t1120, -t1118, t1143, 0, 0, 0, 0, 0, 0, -t1457, t1456, t1179, t1091 + t1343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t492, t599, t531, t491, t530, -t1115, qJD(3) * t639 - t1146 * t816 - t948, qJD(3) * t640 + t1147 * t816 - t949, 0, 0, t1455, t214, t328, t1458, t330, -t1115, qJD(3) * t394 - t1428 - t976, qJD(3) * t396 + t1429 - t977, -t932 + t1089, t121 * qJD(3) + (-t1029 * t874 - t1359 * t603) * t1349 + t930; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1455, t214, t328, t1458, t330, -t1115, qJD(3) * t427 - t1428 - t954, qJD(3) * t426 + t1429 - t953, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1071, -t1116, 0, -t1071, 0, 0, -t1076, t1077, 0, 0, t589, t1128, 0, -t1461, 0, 0, -qJD(2) * t321 + t1164 - t762, -qJD(2) * t323 + t1163 + t765, qJD(2) * t234, qJD(2) * t91 + qJD(4) * t493 - t1137, t335, t397, t276, t334, t275, -t419, -qJD(2) * t166 + qJD(5) * t218 - t1074 - t1201, -qJD(2) * t169 + qJD(5) * t216 - t1200 + t544, -qJD(2) * t81 - t1196 + t1262, -qJD(2) * t34 + qJD(4) * t151 - t1274, t115, t62, t165, t110, t164, -t337, -qJD(2) * t69 - qJD(5) * t56 - qJD(6) * t76 - t1305 - t444, -qJD(2) * t71 - qJD(5) * t55 - qJD(6) * t77 - t1304 + t448, -qJD(2) * t8 - qJD(5) * t12 - t1324 + t1470, -qJD(2) * t2 + qJD(4) * t43 - qJD(5) * t5 - t1321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t876 * t1348, pkin(2) * t1037, 0, 0, 0, 0, 0, 0, 0, 0, t971, t970, t1136, t975, t858, t855, 0, -t858, 0, 0, -qJD(5) * t637 + t1151 * t835 - t1139, -qJD(5) * t638 + t950, t974, t979, t618, t515, 0, -t617, 0, 0, -qJD(5) * t393 - qJD(6) * t424 + t952, qJD(5) * t395 - qJD(6) * t425 + t951, t991 + t1197, -qJD(5) * t120 - t931; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t858, t855, 0, -t858, 0, 0, t862 * t1147, t862 * t1146, 0, 0, t618, t515, 0, -t617, 0, 0, qJD(5) * t696 + t1144 * t848, qJD(5) * t697 + t1144 * t955, 0, qJD(5) * t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1432, t1439, 0, t1123, 0, 0, 0, 0, 0, 0, -t732, t1120, -t1118, t1142, 0, 0, 0, 0, 0, 0, -t1457, t1456, t1179, -t978 + t1193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t492, t599, t531, t491, t530, -t1115, -t1146 * t861 - t915, t1147 * t861 - t916, 0, 0, t1455, t214, t328, t1458, t330, -t1115, -t1427 - t927, t1430 - t926, -t929 + t1089 (-t1028 * t874 - t1359 * t757) * t1349 + t906; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1455, t214, t328, t1458, t330, -t1115, -t1427 - t918, t1430 - t917, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1419, -t580, -t1129, qJD(2) * t389 - qJD(3) * t493 - t1173, 0, 0, 0, 0, 0, 0, -t547 - t1168 + t1182, qJD(5) * t561 - t1108 * t1434 - t1167, t1108 * t577, -qJD(2) * t107 - qJD(3) * t151 - t1337, 0, 0, 0, 0, 0, 0, t1108 * t1433 - t1171 - t1190, -t1107 * t382 - t1108 * t1435 - t1170, -t1108 * t213 - t1176, -qJD(2) * t40 - qJD(3) * t43 - qJD(5) * t25 - t1295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1432, -t1439, 0, t1127, 0, 0, 0, 0, 0, 0, t732, -t1120, t1118, -t1143, 0, 0, 0, 0, 0, 0, t1457, -t1456, -t1179, t1091 - t1343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1432, -t1439, 0, -t1123, 0, 0, 0, 0, 0, 0, t732, -t1120, t1118, -t1142, 0, 0, 0, 0, 0, 0, t1457, -t1456, -t1179, t978 + t1193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t530, t1119 - t1146, 0, 0, 0, 0, 0, 0, 0, 0, t330, t329, 0 (-t1359 * t848 + t874 * t955) * t1349 + t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t330, t329, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t1426 * t996 + t1117, -t1108 * t560 + t1461 * t875, -t290, t1108 * t557 + t1461 * t878, t575, -qJD(2) * t195 - qJD(3) * t218 + qJD(4) * t557 - t1177, -qJD(2) * t193 - qJD(3) * t216 - qJD(4) * t561 - t1178, 0, 0, -t1431, -t74, t209, t1431, t210, t575, qJD(2) * t46 + qJD(3) * t56 + qJD(6) * t51 - t1341 + t369, qJD(2) * t45 + qJD(3) * t55 + qJD(6) * t53 - t1340 + t372, qJD(2) * t9 + qJD(3) * t12 - t1322, qJD(2) * t3 + qJD(3) * t5 + qJD(4) * t25 - t1320; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t491, -t599, -t550, t492, t549, t1115, qJD(3) * t637 + t948, qJD(3) * t638 + t949, 0, 0, t1458, -t214, t371, t1455, t370, t1115, qJD(3) * t393 + t976, -qJD(3) * t395 + t977, t932, qJD(3) * t120 - t930; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t491, -t599, -t550, t492, t549, t1115, t915, t916, 0, 0, t1458, -t214, t371, t1455, t370, t1115, t927, t926, t929, -t906; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t549, -t1119, 0, 0, 0, 0, 0, 0, 0, 0, t370, t373, 0, -t928; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t1354, -pkin(5) * t1034, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1107 * t1354 + t1342, t1283 + (-t1035 - t1034) * pkin(5), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1431, -t74, t209, t1431, t210, t575, qJD(2) * t64 + qJD(3) * t76 - qJD(5) * t51 - t1338 + t369, qJD(2) * t65 + qJD(3) * t77 - qJD(5) * t53 - t1339 + t372, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1458, -t214, t371, t1455, t370, t1115, qJD(3) * t424 + t954, qJD(3) * t425 + t953, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1458, -t214, t371, t1455, t370, t1115, t918, t917, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t370, t373, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1349 * t874 - t1342, pkin(5) * t1035 - t1283, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t11;
