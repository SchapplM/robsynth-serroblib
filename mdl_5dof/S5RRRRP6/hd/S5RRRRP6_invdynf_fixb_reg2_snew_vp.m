% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRRP6_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:55:11
% EndTime: 2019-12-31 21:55:14
% DurationCPUTime: 2.68s
% Computational Cost: add. (12651->209), mult. (25966->277), div. (0->0), fcn. (18414->8), ass. (0->160)
t1451 = qJD(2) + qJD(3);
t1446 = t1451 ^ 2;
t1422 = sin(qJ(3));
t1426 = cos(qJ(3));
t1427 = cos(qJ(2));
t1456 = qJD(1) * t1427;
t1423 = sin(qJ(2));
t1457 = qJD(1) * t1423;
t1390 = t1422 * t1457 - t1426 * t1456;
t1389 = qJD(4) + t1390;
t1465 = qJD(4) + t1389;
t1443 = qJD(2) * t1456;
t1448 = t1423 * qJDD(1);
t1397 = t1443 + t1448;
t1417 = t1427 * qJDD(1);
t1444 = qJD(2) * t1457;
t1398 = t1417 - t1444;
t1440 = t1426 * t1397 + t1422 * t1398;
t1367 = -t1390 * qJD(3) + t1440;
t1464 = t1451 * t1390 - t1367;
t1420 = t1427 ^ 2;
t1430 = qJD(1) ^ 2;
t1424 = sin(qJ(1));
t1428 = cos(qJ(1));
t1406 = t1424 * g(1) - t1428 * g(2);
t1438 = qJDD(1) * pkin(1) + t1406;
t1439 = qJD(2) * pkin(2) - pkin(7) * t1457;
t1369 = t1398 * pkin(2) + (t1420 * pkin(7) + pkin(6)) * t1430 - t1439 * t1457 + t1438;
t1392 = (t1422 * t1427 + t1423 * t1426) * qJD(1);
t1421 = sin(qJ(4));
t1425 = cos(qJ(4));
t1377 = t1421 * t1392 - t1425 * t1451;
t1463 = t1377 ^ 2;
t1379 = t1425 * t1392 + t1421 * t1451;
t1462 = t1379 ^ 2;
t1461 = t1389 ^ 2;
t1460 = t1390 ^ 2;
t1459 = t1392 ^ 2;
t1458 = -2 * qJD(5);
t1455 = t1379 * t1377;
t1454 = t1392 * t1390;
t1453 = t1420 * t1430;
t1452 = t1423 * t1430;
t1450 = qJD(4) - t1389;
t1441 = t1422 * t1397 - t1426 * t1398;
t1354 = (0.2e1 * qJD(3) + qJD(2)) * t1392 + t1441;
t1327 = t1354 * pkin(3) + t1464 * pkin(8) - t1369;
t1407 = -t1428 * g(1) - t1424 * g(2);
t1433 = -t1430 * pkin(1) + qJDD(1) * pkin(6) + t1407;
t1386 = -t1423 * g(3) + t1427 * t1433;
t1366 = -pkin(2) * t1453 + t1398 * pkin(7) - qJD(2) * t1439 + t1386;
t1432 = t1423 * t1433;
t1431 = -t1432 - t1397 * pkin(7) + qJDD(2) * pkin(2) + (qJD(2) * pkin(7) * qJD(1) + pkin(2) * t1452 - g(3)) * t1427;
t1346 = t1426 * t1366 + t1422 * t1431;
t1375 = t1390 * pkin(3) - t1392 * pkin(8);
t1447 = qJDD(2) + qJDD(3);
t1332 = -t1446 * pkin(3) + t1447 * pkin(8) - t1390 * t1375 + t1346;
t1311 = t1421 * t1327 + t1425 * t1332;
t1419 = t1423 ^ 2;
t1449 = t1419 + t1420;
t1310 = t1425 * t1327 - t1421 * t1332;
t1345 = -t1422 * t1366 + t1426 * t1431;
t1442 = t1421 * t1367 - t1425 * t1447;
t1437 = -t1379 * qJD(4) - t1442;
t1436 = -t1425 * t1367 - t1421 * t1447;
t1435 = -t1392 * qJD(3) - qJDD(4) - t1441;
t1331 = -t1447 * pkin(3) - t1446 * pkin(8) + t1392 * t1375 - t1345;
t1347 = -t1435 - t1455;
t1344 = t1450 * t1377 + t1436;
t1429 = qJD(2) ^ 2;
t1411 = t1427 * t1452;
t1409 = -t1429 - t1453;
t1408 = -t1419 * t1430 - t1429;
t1405 = -qJDD(2) + t1411;
t1404 = qJDD(2) + t1411;
t1403 = t1449 * t1430;
t1402 = -t1424 * qJDD(1) - t1428 * t1430;
t1401 = t1428 * qJDD(1) - t1424 * t1430;
t1400 = t1449 * qJDD(1);
t1399 = t1417 - 0.2e1 * t1444;
t1396 = 0.2e1 * t1443 + t1448;
t1393 = t1430 * pkin(6) + t1438;
t1385 = -t1427 * g(3) - t1432;
t1384 = -t1459 - t1446;
t1383 = t1427 * t1405 - t1423 * t1408;
t1382 = -t1423 * t1404 + t1427 * t1409;
t1381 = t1423 * t1405 + t1427 * t1408;
t1380 = t1427 * t1404 + t1423 * t1409;
t1374 = -t1447 - t1454;
t1373 = t1447 - t1454;
t1372 = -t1446 - t1460;
t1370 = t1389 * pkin(4) - t1379 * qJ(5);
t1368 = -t1459 - t1460;
t1365 = -t1423 * t1385 + t1427 * t1386;
t1364 = t1427 * t1385 + t1423 * t1386;
t1360 = -t1461 - t1462;
t1359 = t1426 * t1374 - t1422 * t1384;
t1358 = t1422 * t1374 + t1426 * t1384;
t1357 = -qJD(2) * t1390 - t1440;
t1355 = qJD(2) * t1392 - t1441;
t1353 = -t1461 - t1463;
t1352 = t1426 * t1372 - t1422 * t1373;
t1351 = t1422 * t1372 + t1426 * t1373;
t1350 = -t1462 - t1463;
t1348 = t1435 - t1455;
t1343 = -t1465 * t1377 - t1436;
t1342 = -t1450 * t1379 - t1442;
t1341 = t1465 * t1379 + t1442;
t1340 = -t1423 * t1358 + t1427 * t1359;
t1339 = t1427 * t1358 + t1423 * t1359;
t1338 = t1426 * t1355 - t1422 * t1357;
t1337 = t1422 * t1355 + t1426 * t1357;
t1336 = -t1423 * t1351 + t1427 * t1352;
t1335 = t1427 * t1351 + t1423 * t1352;
t1334 = t1425 * t1348 - t1421 * t1360;
t1333 = t1421 * t1348 + t1425 * t1360;
t1329 = -t1421 * t1347 + t1425 * t1353;
t1328 = t1425 * t1347 + t1421 * t1353;
t1324 = -t1422 * t1345 + t1426 * t1346;
t1323 = t1426 * t1345 + t1422 * t1346;
t1322 = t1425 * t1342 - t1421 * t1344;
t1321 = t1421 * t1342 + t1425 * t1344;
t1320 = -t1423 * t1337 + t1427 * t1338;
t1319 = t1427 * t1337 + t1423 * t1338;
t1318 = t1426 * t1334 + t1422 * t1343;
t1317 = t1422 * t1334 - t1426 * t1343;
t1316 = t1426 * t1329 + t1422 * t1341;
t1315 = t1422 * t1329 - t1426 * t1341;
t1314 = t1426 * t1322 + t1422 * t1350;
t1313 = t1422 * t1322 - t1426 * t1350;
t1312 = -t1437 * pkin(4) - t1463 * qJ(5) + t1379 * t1370 + qJDD(5) + t1331;
t1309 = -t1423 * t1323 + t1427 * t1324;
t1308 = t1427 * t1323 + t1423 * t1324;
t1307 = -t1389 * t1370 + t1437 * qJ(5) + (-pkin(4) * t1377 + t1458) * t1377 + t1311;
t1306 = -t1423 * t1317 + t1427 * t1318;
t1305 = t1427 * t1317 + t1423 * t1318;
t1304 = t1347 * pkin(4) + t1344 * qJ(5) + t1379 * t1458 + t1310;
t1303 = -t1423 * t1315 + t1427 * t1316;
t1302 = t1427 * t1315 + t1423 * t1316;
t1301 = -t1423 * t1313 + t1427 * t1314;
t1300 = t1427 * t1313 + t1423 * t1314;
t1299 = -t1421 * t1310 + t1425 * t1311;
t1298 = t1425 * t1310 + t1421 * t1311;
t1297 = t1428 * t1306 + t1424 * t1333;
t1296 = t1424 * t1306 - t1428 * t1333;
t1295 = t1428 * t1303 + t1424 * t1328;
t1294 = t1424 * t1303 - t1428 * t1328;
t1293 = t1426 * t1299 + t1422 * t1331;
t1292 = t1422 * t1299 - t1426 * t1331;
t1291 = t1428 * t1301 + t1424 * t1321;
t1290 = t1424 * t1301 - t1428 * t1321;
t1289 = -t1421 * t1304 + t1425 * t1307;
t1288 = t1425 * t1304 + t1421 * t1307;
t1287 = t1426 * t1289 + t1422 * t1312;
t1286 = t1422 * t1289 - t1426 * t1312;
t1285 = -t1423 * t1292 + t1427 * t1293;
t1284 = t1427 * t1292 + t1423 * t1293;
t1283 = -t1423 * t1286 + t1427 * t1287;
t1282 = t1427 * t1286 + t1423 * t1287;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1402, -t1401, 0, -t1424 * t1406 + t1428 * t1407, 0, 0, 0, 0, 0, 0, t1428 * t1382 - t1424 * t1399, t1428 * t1383 + t1424 * t1396, t1428 * t1400 - t1424 * t1403, t1428 * t1365 - t1424 * t1393, 0, 0, 0, 0, 0, 0, t1428 * t1336 + t1424 * t1354, t1428 * t1340 - t1424 * t1464, t1428 * t1320 + t1424 * t1368, t1428 * t1309 - t1424 * t1369, 0, 0, 0, 0, 0, 0, t1295, t1297, t1291, t1428 * t1285 + t1424 * t1298, 0, 0, 0, 0, 0, 0, t1295, t1297, t1291, t1428 * t1283 + t1424 * t1288; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1401, t1402, 0, t1428 * t1406 + t1424 * t1407, 0, 0, 0, 0, 0, 0, t1424 * t1382 + t1428 * t1399, t1424 * t1383 - t1428 * t1396, t1424 * t1400 + t1428 * t1403, t1424 * t1365 + t1428 * t1393, 0, 0, 0, 0, 0, 0, t1424 * t1336 - t1428 * t1354, t1424 * t1340 + t1428 * t1464, t1424 * t1320 - t1428 * t1368, t1424 * t1309 + t1428 * t1369, 0, 0, 0, 0, 0, 0, t1294, t1296, t1290, t1424 * t1285 - t1428 * t1298, 0, 0, 0, 0, 0, 0, t1294, t1296, t1290, t1424 * t1283 - t1428 * t1288; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1380, t1381, 0, t1364, 0, 0, 0, 0, 0, 0, t1335, t1339, t1319, t1308, 0, 0, 0, 0, 0, 0, t1302, t1305, t1300, t1284, 0, 0, 0, 0, 0, 0, t1302, t1305, t1300, t1282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1430, -qJDD(1), 0, t1407, 0, 0, 0, 0, 0, 0, t1382, t1383, t1400, t1365, 0, 0, 0, 0, 0, 0, t1336, t1340, t1320, t1309, 0, 0, 0, 0, 0, 0, t1303, t1306, t1301, t1285, 0, 0, 0, 0, 0, 0, t1303, t1306, t1301, t1283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1430, 0, t1406, 0, 0, 0, 0, 0, 0, t1399, -t1396, t1403, t1393, 0, 0, 0, 0, 0, 0, -t1354, t1464, -t1368, t1369, 0, 0, 0, 0, 0, 0, -t1328, -t1333, -t1321, -t1298, 0, 0, 0, 0, 0, 0, -t1328, -t1333, -t1321, -t1288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1380, t1381, 0, t1364, 0, 0, 0, 0, 0, 0, t1335, t1339, t1319, t1308, 0, 0, 0, 0, 0, 0, t1302, t1305, t1300, t1284, 0, 0, 0, 0, 0, 0, t1302, t1305, t1300, t1282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1409, t1405, t1417, t1386, 0, 0, 0, 0, 0, 0, t1352, t1359, t1338, t1324, 0, 0, 0, 0, 0, 0, t1316, t1318, t1314, t1293, 0, 0, 0, 0, 0, 0, t1316, t1318, t1314, t1287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1404, t1408, -t1448, t1385, 0, 0, 0, 0, 0, 0, t1351, t1358, t1337, t1323, 0, 0, 0, 0, 0, 0, t1315, t1317, t1313, t1292, 0, 0, 0, 0, 0, 0, t1315, t1317, t1313, t1286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1399, t1396, -t1403, -t1393, 0, 0, 0, 0, 0, 0, t1354, -t1464, t1368, -t1369, 0, 0, 0, 0, 0, 0, t1328, t1333, t1321, t1298, 0, 0, 0, 0, 0, 0, t1328, t1333, t1321, t1288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1372, t1374, t1355, t1346, 0, 0, 0, 0, 0, 0, t1329, t1334, t1322, t1299, 0, 0, 0, 0, 0, 0, t1329, t1334, t1322, t1289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1373, t1384, t1357, t1345, 0, 0, 0, 0, 0, 0, -t1341, -t1343, -t1350, -t1331, 0, 0, 0, 0, 0, 0, -t1341, -t1343, -t1350, -t1312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1354, -t1464, t1368, -t1369, 0, 0, 0, 0, 0, 0, t1328, t1333, t1321, t1298, 0, 0, 0, 0, 0, 0, t1328, t1333, t1321, t1288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1353, t1348, t1342, t1311, 0, 0, 0, 0, 0, 0, t1353, t1348, t1342, t1307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1347, t1360, t1344, t1310, 0, 0, 0, 0, 0, 0, t1347, t1360, t1344, t1304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1341, t1343, t1350, t1331, 0, 0, 0, 0, 0, 0, t1341, t1343, t1350, t1312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1353, t1348, t1342, t1307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1347, t1360, t1344, t1304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1341, t1343, t1350, t1312;];
f_new_reg = t1;
