% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRPRP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:45
% EndTime: 2019-12-31 19:51:48
% DurationCPUTime: 2.74s
% Computational Cost: add. (7462->175), mult. (10742->215), div. (0->0), fcn. (7047->8), ass. (0->121)
t1494 = qJD(1) + qJD(2);
t1497 = sin(pkin(8));
t1498 = cos(pkin(8));
t1499 = sin(qJ(4));
t1502 = cos(qJ(4));
t1548 = -t1497 * t1499 + t1498 * t1502;
t1463 = t1548 * t1494;
t1513 = t1497 * t1502 + t1498 * t1499;
t1465 = t1513 * t1494;
t1529 = t1465 * t1463;
t1439 = -qJDD(4) + t1529;
t1462 = t1465 ^ 2;
t1505 = qJD(4) ^ 2;
t1536 = -t1462 - t1505;
t1414 = t1499 * t1439 + t1502 * t1536;
t1416 = t1502 * t1439 - t1499 * t1536;
t1386 = t1497 * t1414 - t1498 * t1416;
t1500 = sin(qJ(2));
t1503 = cos(qJ(2));
t1492 = qJDD(1) + qJDD(2);
t1461 = t1513 * t1492;
t1550 = 2 * qJD(4);
t1511 = t1463 * t1550 + t1461;
t1370 = t1500 * t1386 + t1503 * t1511;
t1371 = t1503 * t1386 - t1500 * t1511;
t1501 = sin(qJ(1));
t1504 = cos(qJ(1));
t1562 = t1504 * t1370 + t1501 * t1371;
t1561 = t1501 * t1370 - t1504 * t1371;
t1515 = t1548 * t1492;
t1530 = t1465 * qJD(4);
t1440 = -t1515 + 0.2e1 * t1530;
t1438 = qJDD(4) + t1529;
t1443 = t1463 ^ 2;
t1537 = -t1443 - t1505;
t1544 = -t1499 * t1438 + t1502 * t1537;
t1547 = t1502 * t1438 + t1499 * t1537;
t1552 = -t1497 * t1547 + t1498 * t1544;
t1557 = t1500 * t1440 + t1503 * t1552;
t1558 = -t1503 * t1440 + t1500 * t1552;
t1560 = -t1501 * t1558 + t1504 * t1557;
t1559 = t1501 * t1557 + t1504 * t1558;
t1396 = t1498 * t1414 + t1497 * t1416;
t1423 = t1462 + t1443;
t1534 = t1499 * t1461 + t1502 * t1515;
t1535 = -t1502 * t1461 + t1499 * t1515;
t1543 = -t1497 * t1535 + t1498 * t1534;
t1553 = -t1500 * t1423 + t1503 * t1543;
t1554 = t1503 * t1423 + t1500 * t1543;
t1556 = -t1501 * t1554 + t1504 * t1553;
t1555 = t1501 * t1553 + t1504 * t1554;
t1551 = t1497 * t1544 + t1498 * t1547;
t1484 = -t1504 * g(1) - t1501 * g(2);
t1506 = qJD(1) ^ 2;
t1477 = -t1506 * pkin(1) + t1484;
t1483 = t1501 * g(1) - t1504 * g(2);
t1512 = qJDD(1) * pkin(1) + t1483;
t1451 = t1503 * t1477 + t1500 * t1512;
t1491 = t1494 ^ 2;
t1549 = -t1491 * pkin(2) + t1492 * qJ(3) + 0.2e1 * qJD(3) * t1494 + t1451;
t1519 = t1503 * t1492;
t1475 = t1500 * t1491 - t1519;
t1520 = t1500 * t1492;
t1514 = -t1503 * t1491 - t1520;
t1546 = t1501 * t1475 + t1504 * t1514;
t1545 = t1504 * t1475 - t1501 * t1514;
t1542 = t1497 * t1534 + t1498 * t1535;
t1493 = t1498 ^ 2;
t1507 = t1497 ^ 2;
t1517 = t1493 + t1507;
t1472 = t1517 * t1491;
t1533 = t1498 * g(3);
t1528 = t1491 * t1498;
t1527 = t1493 * t1491;
t1526 = t1497 * t1492;
t1488 = t1498 * t1492;
t1518 = t1507 * t1491;
t1425 = -t1497 * g(3) + t1549 * t1498;
t1412 = -pkin(3) * t1527 + pkin(7) * t1488 + t1425;
t1510 = -t1533 + (pkin(3) * t1528 - pkin(7) * t1492 - t1549) * t1497;
t1394 = t1502 * t1412 + t1499 * t1510;
t1393 = -t1499 * t1412 + t1502 * t1510;
t1450 = -t1500 * t1477 + t1503 * t1512;
t1434 = -t1492 * pkin(2) - t1491 * qJ(3) + qJDD(3) - t1450;
t1419 = -pkin(3) * t1488 + t1434 + (-t1518 - t1527) * pkin(7);
t1480 = t1497 * t1528;
t1479 = -t1501 * qJDD(1) - t1504 * t1506;
t1478 = t1504 * qJDD(1) - t1501 * t1506;
t1468 = t1517 * t1492;
t1467 = t1498 * t1472;
t1466 = t1497 * t1472;
t1449 = -t1503 * t1467 - t1498 * t1520;
t1448 = t1503 * t1466 + t1497 * t1520;
t1447 = -t1500 * t1467 + t1498 * t1519;
t1446 = t1500 * t1466 - t1497 * t1519;
t1445 = t1503 * t1468 - t1500 * t1472;
t1444 = t1500 * t1468 + t1503 * t1472;
t1433 = -t1463 * pkin(4) - t1465 * qJ(5);
t1424 = -t1497 * t1549 - t1533;
t1418 = -t1500 * t1450 + t1503 * t1451;
t1417 = t1503 * t1450 + t1500 * t1451;
t1400 = -t1497 * t1424 + t1498 * t1425;
t1399 = t1498 * t1424 + t1497 * t1425;
t1390 = t1503 * t1400 + t1500 * t1434;
t1389 = t1500 * t1400 - t1503 * t1434;
t1378 = -(t1515 - t1530) * pkin(4) + (pkin(4) * qJD(4) - (2 * qJD(5))) * t1465 + t1419 - t1511 * qJ(5);
t1377 = -qJDD(4) * pkin(4) - t1505 * qJ(5) + t1465 * t1433 + qJDD(5) - t1393;
t1376 = -t1505 * pkin(4) + qJDD(4) * qJ(5) + (qJD(5) * t1550) + t1463 * t1433 + t1394;
t1369 = -t1499 * t1393 + t1502 * t1394;
t1368 = t1502 * t1393 + t1499 * t1394;
t1367 = t1502 * t1376 + t1499 * t1377;
t1366 = t1499 * t1376 - t1502 * t1377;
t1365 = -t1497 * t1368 + t1498 * t1369;
t1364 = t1498 * t1368 + t1497 * t1369;
t1363 = t1503 * t1365 + t1500 * t1419;
t1362 = t1500 * t1365 - t1503 * t1419;
t1361 = -t1497 * t1366 + t1498 * t1367;
t1360 = t1498 * t1366 + t1497 * t1367;
t1359 = t1503 * t1361 + t1500 * t1378;
t1358 = t1500 * t1361 - t1503 * t1378;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1479, -t1478, 0, -t1501 * t1483 + t1504 * t1484, 0, 0, 0, 0, 0, 0, t1546, t1545, 0, -t1501 * t1417 + t1504 * t1418, 0, 0, 0, 0, 0, 0, -t1501 * t1447 + t1504 * t1449, -t1501 * t1446 + t1504 * t1448, -t1501 * t1444 + t1504 * t1445, -t1501 * t1389 + t1504 * t1390, 0, 0, 0, 0, 0, 0, t1560, t1561, t1556, -t1501 * t1362 + t1504 * t1363, 0, 0, 0, 0, 0, 0, t1560, t1556, -t1561, -t1501 * t1358 + t1504 * t1359; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1478, t1479, 0, t1504 * t1483 + t1501 * t1484, 0, 0, 0, 0, 0, 0, -t1545, t1546, 0, t1504 * t1417 + t1501 * t1418, 0, 0, 0, 0, 0, 0, t1504 * t1447 + t1501 * t1449, t1504 * t1446 + t1501 * t1448, t1504 * t1444 + t1501 * t1445, t1504 * t1389 + t1501 * t1390, 0, 0, 0, 0, 0, 0, t1559, -t1562, t1555, t1504 * t1362 + t1501 * t1363, 0, 0, 0, 0, 0, 0, t1559, t1555, t1562, t1504 * t1358 + t1501 * t1359; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1399, 0, 0, 0, 0, 0, 0, t1551, t1396, t1542, t1364, 0, 0, 0, 0, 0, 0, t1551, t1542, -t1396, t1360; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1506, -qJDD(1), 0, t1484, 0, 0, 0, 0, 0, 0, t1514, t1475, 0, t1418, 0, 0, 0, 0, 0, 0, t1449, t1448, t1445, t1390, 0, 0, 0, 0, 0, 0, t1557, -t1371, t1553, t1363, 0, 0, 0, 0, 0, 0, t1557, t1553, t1371, t1359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1506, 0, t1483, 0, 0, 0, 0, 0, 0, -t1475, t1514, 0, t1417, 0, 0, 0, 0, 0, 0, t1447, t1446, t1444, t1389, 0, 0, 0, 0, 0, 0, t1558, -t1370, t1554, t1362, 0, 0, 0, 0, 0, 0, t1558, t1554, t1370, t1358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1399, 0, 0, 0, 0, 0, 0, t1551, t1396, t1542, t1364, 0, 0, 0, 0, 0, 0, t1551, t1542, -t1396, t1360; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1491, -t1492, 0, t1451, 0, 0, 0, 0, 0, 0, -t1467, t1466, t1468, t1400, 0, 0, 0, 0, 0, 0, t1552, -t1386, t1543, t1365, 0, 0, 0, 0, 0, 0, t1552, t1543, t1386, t1361; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1492, -t1491, 0, t1450, 0, 0, 0, 0, 0, 0, t1488, -t1526, t1472, -t1434, 0, 0, 0, 0, 0, 0, -t1440, -t1511, t1423, -t1419, 0, 0, 0, 0, 0, 0, -t1440, t1423, t1511, -t1378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1399, 0, 0, 0, 0, 0, 0, t1551, t1396, t1542, t1364, 0, 0, 0, 0, 0, 0, t1551, t1542, -t1396, t1360; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1527, t1480, t1488, t1425, 0, 0, 0, 0, 0, 0, t1544, t1416, t1534, t1369, 0, 0, 0, 0, 0, 0, t1544, t1534, -t1416, t1367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1480, -t1518, -t1526, t1424, 0, 0, 0, 0, 0, 0, t1547, t1414, t1535, t1368, 0, 0, 0, 0, 0, 0, t1547, t1535, -t1414, t1366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1488, t1526, -t1472, t1434, 0, 0, 0, 0, 0, 0, t1440, t1511, -t1423, t1419, 0, 0, 0, 0, 0, 0, t1440, -t1423, -t1511, t1378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1537, t1439, t1515, t1394, 0, 0, 0, 0, 0, 0, t1537, t1515, -t1439, t1376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1438, t1536, -t1461, t1393, 0, 0, 0, 0, 0, 0, t1438, -t1461, -t1536, -t1377; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1440, t1511, -t1423, t1419, 0, 0, 0, 0, 0, 0, t1440, -t1423, -t1511, t1378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1537, t1515, -t1439, t1376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1440, -t1423, -t1511, t1378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1438, t1461, t1536, t1377;];
f_new_reg = t1;
