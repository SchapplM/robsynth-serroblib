% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:27
% EndTime: 2019-03-09 15:29:28
% DurationCPUTime: 51.04s
% Computational Cost: add. (111321->1131), mult. (141966->1452), div. (0->0), fcn. (147188->8), ass. (0->677)
t979 = qJ(2) + qJ(3);
t960 = cos(t979);
t985 = cos(qJ(1));
t1281 = t960 * t985;
t982 = sin(qJ(1));
t1283 = t960 * t982;
t1464 = -pkin(8) - pkin(7);
t984 = cos(qJ(2));
t1357 = pkin(2) * t984;
t956 = pkin(1) + t1357;
t1166 = -t985 * t1464 - t982 * t956;
t1342 = rSges(5,3) + qJ(4);
t1443 = rSges(5,1) + pkin(3);
t959 = sin(t979);
t1484 = -t1342 * t959 - t1443 * t960;
t975 = t985 * rSges(5,2);
t609 = t1484 * t982 + t1166 + t975;
t1348 = t982 * rSges(5,2);
t957 = t982 * t1464;
t610 = t1348 - t957 + (-t1484 + t956) * t985;
t1209 = t609 * t1281 + t610 * t1283;
t1343 = rSges(6,1) + qJ(4);
t1465 = pkin(3) + pkin(4);
t1483 = -t1343 * t959 - t1465 * t960;
t938 = rSges(6,2) * t1283;
t574 = t938 + (-rSges(6,3) - qJ(5)) * t985 + t1483 * t982 + t1166;
t961 = t982 * qJ(5);
t1161 = -t957 - t961;
t1164 = -rSges(6,2) * t1281 - t982 * rSges(6,3);
t575 = (-t1483 + t956) * t985 + t1161 + t1164;
t1210 = t574 * t1281 + t575 * t1283;
t1322 = qJ(5) * t985;
t1124 = rSges(7,3) + pkin(9) + t1465;
t1355 = pkin(5) + qJ(4);
t1485 = -t1124 * t960 - t1355 * t959;
t983 = cos(qJ(6));
t1256 = t983 * t985;
t980 = sin(qJ(6));
t1259 = t982 * t980;
t872 = t1259 * t959 - t1256;
t1258 = t982 * t983;
t1280 = t980 * t985;
t873 = t1258 * t959 + t1280;
t1498 = -t873 * rSges(7,1) + t872 * rSges(7,2);
t483 = t1485 * t982 + t1166 - t1322 + t1498;
t874 = -t1280 * t959 - t1258;
t875 = t1256 * t959 - t1259;
t1066 = t875 * rSges(7,1) + t874 * rSges(7,2);
t484 = (-t1485 + t956) * t985 + t1066 + t1161;
t1218 = t483 * t1281 + t484 * t1283;
t1292 = t959 * t985;
t1293 = t959 * t982;
t1466 = m(7) / 0.2e1;
t1467 = m(6) / 0.2e1;
t1468 = m(5) / 0.2e1;
t1356 = pkin(4) * t959;
t901 = pkin(3) * t959 - qJ(4) * t960;
t1106 = t901 + t1356;
t981 = sin(qJ(2));
t1358 = pkin(2) * t981;
t1032 = t1106 + t1358;
t1351 = rSges(7,1) * t983;
t1065 = -rSges(7,2) * t980 + t1351;
t1350 = rSges(7,3) * t959;
t755 = -t1065 * t960 + t1350;
t1501 = -pkin(5) * t960 + pkin(9) * t959 + t755;
t1016 = t1032 + t1501;
t553 = t1016 * t982;
t555 = t1016 * t985;
t1068 = rSges(6,1) * t960 + rSges(6,2) * t959;
t1024 = t1032 - t1068;
t625 = t1024 * t982;
t627 = t1024 * t985;
t1171 = rSges(5,1) * t959 - rSges(5,3) * t960 + t901;
t1077 = t1171 + t1358;
t662 = t1077 * t982;
t664 = t1077 * t985;
t1134 = (-t553 * t1292 + t1293 * t555 + t1218) * t1466 + (-t625 * t1292 + t1293 * t627 + t1210) * t1467 + (-t662 * t1292 + t1293 * t664 + t1209) * t1468;
t1007 = t1356 + (-t1351 - t1355) * t960;
t1284 = t960 * t980;
t1152 = rSges(7,2) * t1284;
t1168 = rSges(7,3) * t1293 + t982 * t1152;
t944 = pkin(9) * t1293;
t945 = pkin(3) * t1293;
t1080 = t944 + t945 + t1168;
t544 = (t1007 + t1358) * t982 + t1080;
t1006 = -t1124 * t959 - t1152;
t1282 = t960 * t983;
t1153 = rSges(7,1) * t1282;
t914 = t985 * t1153;
t924 = qJ(4) * t1281;
t946 = pkin(5) * t1281;
t1126 = t914 + t924 + t946;
t545 = (t1006 - t1358) * t985 + t1126;
t1015 = -t1343 * t960 + (-rSges(6,2) + pkin(4)) * t959;
t614 = t945 + (t1015 + t1358) * t982;
t1163 = rSges(6,1) * t1281 + rSges(6,2) * t1292;
t1125 = t924 + t1163;
t1145 = t1465 * t959;
t615 = (-t1145 - t1358) * t985 + t1125;
t1104 = t1342 * t960;
t937 = rSges(5,1) * t1293;
t1165 = t937 + t945;
t650 = (-t1104 + t1358) * t982 + t1165;
t1123 = t1443 * t959;
t941 = rSges(5,3) * t1281;
t1167 = t924 + t941;
t651 = (-t1123 - t1358) * t985 + t1167;
t1135 = ((t544 * t985 + t545 * t982) * t959 + t1218) * t1466 + ((t614 * t985 + t615 * t982) * t959 + t1210) * t1467 + ((t650 * t985 + t651 * t982) * t959 + t1209) * t1468;
t24 = t1135 - t1134;
t1575 = t24 * qJD(1);
t1449 = t959 / 0.2e1;
t1537 = t960 / 0.2e1;
t1051 = Icges(7,5) * t983 - Icges(7,6) * t980;
t1001 = -Icges(7,3) * t959 + t1051 * t960;
t1331 = Icges(7,4) * t875;
t637 = Icges(7,2) * t874 + Icges(7,6) * t1281 + t1331;
t839 = Icges(7,4) * t874;
t640 = Icges(7,1) * t875 + Icges(7,5) * t1281 + t839;
t1046 = t637 * t980 - t640 * t983;
t1026 = t1001 * t985 - t1046;
t634 = Icges(7,5) * t875 + Icges(7,6) * t874 + Icges(7,3) * t1281;
t1330 = Icges(7,4) * t983;
t1057 = -Icges(7,2) * t980 + t1330;
t1002 = -Icges(7,6) * t959 + t1057 * t960;
t696 = t1002 * t985;
t1060 = Icges(7,1) * t983 - Icges(7,4) * t980;
t1328 = Icges(7,5) * t959;
t1003 = t1060 * t960 - t1328;
t698 = t1003 * t985;
t291 = (t696 * t980 - t698 * t983 + t634) * t960 + t1026 * t959;
t750 = Icges(7,6) * t960 + t1057 * t959;
t752 = Icges(7,5) * t960 + t1060 * t959;
t930 = Icges(7,4) * t1284;
t753 = -Icges(7,1) * t1282 + t1328 + t930;
t1044 = -t1002 * t980 - t753 * t983;
t1025 = Icges(7,3) * t960 + t1051 * t959 - t1044;
t1312 = t1001 * t959;
t997 = t1025 * t960 + t1312;
t340 = t874 * t750 + t875 * t752 + t985 * t997;
t357 = (t750 * t980 - t752 * t983 - t1001) * t960 + t1025 * t959;
t632 = Icges(7,5) * t873 - Icges(7,6) * t872 + Icges(7,3) * t1283;
t1314 = t632 * t959;
t838 = Icges(7,4) * t873;
t636 = Icges(7,2) * t872 - Icges(7,6) * t1283 - t838;
t837 = Icges(7,4) * t872;
t638 = Icges(7,1) * t873 + Icges(7,5) * t1283 - t837;
t1540 = t636 * t980 + t638 * t983;
t416 = t1540 * t960 - t1314;
t1313 = t634 * t959;
t417 = t1046 * t960 + t1313;
t458 = -t1001 * t1283 + t1002 * t872 + t753 * t873;
t460 = -t1001 * t1281 - t1002 * t874 + t875 * t753;
t502 = t1044 * t960 - t1312;
t1574 = t1449 * t357 + t1537 * t502 - (-t416 + t458) * t1293 / 0.4e1 - (t417 + t460) * t1292 / 0.4e1 + t1281 * (t291 + t340) / 0.4e1;
t840 = (Icges(7,5) * t980 + Icges(7,6) * t983) * t960;
t1301 = t959 * t840;
t672 = -rSges(7,1) * t872 - rSges(7,2) * t873;
t673 = rSges(7,1) * t874 - rSges(7,2) * t875;
t847 = Icges(7,2) * t1282 + t930;
t854 = (Icges(7,1) * t980 + t1330) * t960;
t191 = (-(t854 / 0.2e1 + t1002 / 0.2e1) * t983 + (t753 / 0.2e1 + t847 / 0.2e1) * t980) * t960 + m(7) * (-t483 * t672 + t484 * t673) + t1301 / 0.2e1;
t1573 = t191 * qJD(1);
t952 = Icges(5,5) * t959;
t1492 = Icges(5,1) * t960 + t952;
t1334 = Icges(4,4) * t959;
t900 = Icges(4,1) * t960 - t1334;
t778 = Icges(4,5) * t982 + t900 * t985;
t1572 = Icges(5,4) * t982 + t1492 * t985 + t778;
t1332 = Icges(6,4) * t960;
t953 = Icges(4,4) * t960;
t899 = Icges(4,1) * t959 + t953;
t1571 = -Icges(6,2) * t959 - t1332 - t899;
t1035 = t1106 + t1501;
t570 = t1035 * t982;
t572 = t1035 * t985;
t1075 = t1106 - t1068;
t652 = t1075 * t982;
t654 = t1075 * t985;
t689 = t1171 * t982;
t691 = t1171 * t985;
t1132 = (-t570 * t1292 + t1293 * t572 + t1218) * t1466 + (-t652 * t1292 + t1293 * t654 + t1210) * t1467 + (-t689 * t1292 + t1293 * t691 + t1209) * t1468;
t557 = t1007 * t982 + t1080;
t558 = t1006 * t985 + t1126;
t647 = t1015 * t982 + t945;
t648 = -t1145 * t985 + t1125;
t676 = -t1104 * t982 + t1165;
t677 = -t1123 * t985 + t1167;
t1133 = ((t557 * t985 + t558 * t982) * t959 + t1218) * t1466 + ((t647 * t985 + t648 * t982) * t959 + t1210) * t1467 + ((t676 * t985 + t677 * t982) * t959 + t1209) * t1468;
t1518 = t1132 - t1133;
t1570 = t1518 * qJD(1);
t381 = t1283 * t632 + t636 * t872 + t638 * t873;
t382 = t634 * t1283 - t872 * t637 + t873 * t640;
t1050 = t381 * t982 + t382 * t985;
t168 = t1050 * t960 + t458 * t959;
t1329 = Icges(5,5) * t960;
t888 = Icges(5,3) * t959 + t1329;
t763 = -Icges(5,6) * t985 + t888 * t982;
t931 = Icges(6,4) * t1283;
t773 = Icges(6,1) * t1293 + Icges(6,5) * t985 - t931;
t1569 = t763 + t773;
t767 = Icges(6,4) * t1293 - Icges(6,2) * t1283 + Icges(6,6) * t985;
t775 = -Icges(5,4) * t985 + t1492 * t982;
t1568 = t767 - t775;
t229 = -t381 * t985 + t382 * t982;
t1116 = t1283 / 0.4e1;
t1118 = -t1283 / 0.4e1;
t383 = t632 * t1281 - t874 * t636 + t875 * t638;
t384 = t634 * t1281 + t874 * t637 + t875 * t640;
t1551 = -t383 * t985 + t384 * t982;
t1566 = (t1116 + t1118) * t1551;
t1448 = -t982 / 0.2e1;
t1445 = -t985 / 0.2e1;
t1262 = t982 * t625;
t1267 = t982 * t553;
t419 = t555 * t985 + t1267;
t1223 = t419 * t1466 + (t627 * t985 + t1262) * t1467;
t1225 = (-t982 * t544 + t545 * t985) * t1466 + (-t982 * t614 + t615 * t985) * t1467;
t88 = t1225 - t1223;
t1565 = t88 * qJD(1);
t895 = Icges(6,1) * t959 - t1332;
t1564 = t888 + t895;
t1261 = t982 * t652;
t1266 = t982 * t570;
t440 = t572 * t985 + t1266;
t1219 = t440 * t1466 + (t654 * t985 + t1261) * t1467;
t1222 = (-t982 * t557 + t558 * t985) * t1466 + (-t982 * t647 + t648 * t985) * t1467;
t106 = t1222 - t1219;
t1563 = qJD(1) * t106;
t892 = Icges(5,4) * t960 + Icges(5,6) * t959;
t1305 = t892 * t985;
t929 = Icges(5,5) * t1281;
t764 = Icges(5,6) * t982 + Icges(5,3) * t1292 + t929;
t889 = Icges(4,5) * t960 - Icges(4,6) * t959;
t1306 = t889 * t985;
t766 = Icges(4,3) * t982 + t1306;
t1562 = t764 * t1292 + (Icges(5,2) * t982 + t1305 + t766) * t982 + t1572 * t1281;
t887 = -Icges(5,3) * t960 + t952;
t1561 = t887 + t1492 + t900;
t1333 = Icges(6,4) * t959;
t1061 = Icges(6,1) * t960 + t1333;
t893 = Icges(4,2) * t960 + t1334;
t1560 = -Icges(6,2) * t960 + t1061 + t1333 + t893;
t765 = Icges(4,5) * t1283 - Icges(4,6) * t1293 - Icges(4,3) * t985;
t933 = Icges(4,4) * t1293;
t777 = Icges(4,1) * t1283 - Icges(4,5) * t985 - t933;
t1559 = -t773 * t1292 - t982 * t765 + (t767 - t777) * t1281;
t1041 = -t767 * t960 + t773 * t959;
t774 = -Icges(6,5) * t982 + t895 * t985;
t713 = t774 * t1292;
t886 = Icges(6,5) * t959 - Icges(6,6) * t960;
t1307 = t886 * t985;
t762 = -Icges(6,3) * t982 + t1307;
t1091 = t982 * t762 - t713;
t709 = t778 * t1283;
t1092 = t985 * t766 - t709;
t771 = Icges(4,4) * t1283 - Icges(4,2) * t1293 - Icges(4,6) * t985;
t1309 = t771 * t959;
t1444 = t985 / 0.2e1;
t1447 = t982 / 0.2e1;
t1495 = (t763 * t959 + t775 * t960) * t982;
t1525 = t1445 * t1551 + t1448 * t229;
t1260 = t982 * t892;
t728 = t982 * (-Icges(5,2) * t985 + t1260);
t527 = t775 * t1281 + t763 * t1292 + t728;
t1526 = t527 * t985;
t761 = Icges(6,5) * t1293 - Icges(6,6) * t1283 + Icges(6,3) * t985;
t1541 = t1041 * t982 + t985 * t761 + t1495 + t1562;
t932 = Icges(6,4) * t1292;
t768 = -Icges(6,2) * t1281 - Icges(6,6) * t982 + t932;
t894 = -Icges(4,2) * t959 + t953;
t772 = Icges(4,6) * t982 + t894 * t985;
t1542 = -t1283 * t768 - t1293 * t772 - t1092;
t1543 = -t1281 * t768 - t1292 * t772 - t1091 + t1562;
t1544 = t768 * t960 + t772 * t959 + t761 - t765;
t1549 = t1292 * t771 + t982 * t761 + t1559;
t993 = (t1543 * t982 + t1549 * t985 - t1526) * t1444 + (-t1526 + t1551 + (-t709 + (t766 + t1309) * t985 + t1542 + t1559) * t985 + (-t1495 + t713 + (-t1041 - t762) * t982 + t1541) * t982) * t1445 + (-t229 + (t1091 + (t765 + t1544) * t985 - t1541 + t1543) * t985 + (-(t777 * t960 - t1309) * t985 + t1092 + t527 - t728 + t1544 * t982 + t1542 - t1549) * t982) * t1447 - t1525;
t897 = Icges(5,1) * t959 - t1329;
t1558 = t894 + t897 - t1571;
t1557 = (-Icges(6,5) - Icges(4,6) + Icges(5,6)) * t960 + (-Icges(5,4) - Icges(4,5) - Icges(6,6)) * t959;
t1556 = -Icges(6,1) * t1281 - t768 - t932 + (t887 - t893) * t985 + t1572;
t1555 = -Icges(4,2) * t1283 - t1568 + t777 - t933 + (-t1061 + t887) * t982;
t1554 = -Icges(5,1) * t1292 + t1571 * t985 + t764 - t772 + t774 + t929;
t1553 = -Icges(6,2) * t1293 + t1569 - t771 - t931 + (-t897 - t899) * t982;
t641 = rSges(7,3) * t1283 - t1498;
t534 = t755 * t1283 - t641 * t959;
t1049 = t982 * t383 + t384 * t985;
t1552 = t1049 * t960 + t460 * t959;
t977 = t982 ^ 2;
t978 = t985 ^ 2;
t1160 = t977 + t978;
t338 = -t483 * t985 - t484 * t982;
t228 = m(6) * (-t574 * t985 - t575 * t982) + m(7) * t338;
t1547 = qJD(1) * t228;
t1469 = m(4) / 0.2e1;
t1539 = -t959 / 0.2e1;
t1538 = -t960 / 0.2e1;
t1536 = t982 / 0.4e1;
t1535 = -t985 / 0.4e1;
t1027 = t1001 * t982 + t1540;
t695 = t1002 * t982;
t697 = t1003 * t982;
t290 = (t695 * t980 - t697 * t983 + t632) * t960 + t1027 * t959;
t339 = -t750 * t872 + t752 * t873 + t982 * t997;
t1524 = t290 + t339;
t1335 = Icges(3,4) * t981;
t916 = Icges(3,2) * t984 + t1335;
t919 = Icges(3,1) * t984 - t1335;
t1521 = (t916 / 0.2e1 - t919 / 0.2e1) * t981;
t1158 = qJD(2) + qJD(3);
t781 = rSges(4,1) * t1283 - rSges(4,2) * t1293 - t985 * rSges(4,3);
t674 = -t781 + t1166;
t1103 = -rSges(4,2) * t1292 + t982 * rSges(4,3);
t1352 = rSges(4,1) * t960;
t675 = -t957 + (t956 + t1352) * t985 + t1103;
t908 = -rSges(4,2) * t959 + t1352;
t1021 = (-t674 * t985 - t675 * t982) * t908;
t905 = pkin(3) * t960 + qJ(4) * t959;
t907 = rSges(5,1) * t960 + rSges(5,3) * t959;
t1170 = -t905 - t907;
t690 = t1170 * t982;
t692 = t1170 * t985;
t1214 = t692 * t609 + t690 * t610;
t1105 = -pkin(4) * t960 - t905;
t902 = rSges(6,1) * t959 - rSges(6,2) * t960;
t1074 = t1105 - t902;
t653 = t1074 * t982;
t655 = t1074 * t985;
t1217 = t655 * t574 + t653 * t575;
t1349 = rSges(7,3) * t960;
t754 = t1065 * t959 + t1349;
t909 = pkin(5) * t959 + pkin(9) * t960;
t1502 = -t754 - t909;
t1036 = t1105 + t1502;
t571 = t1036 * t982;
t573 = t1036 * t985;
t1226 = t573 * t483 + t571 * t484;
t904 = rSges(4,1) * t959 + rSges(4,2) * t960;
t1018 = t904 + t1358;
t1496 = t1018 * t985;
t1497 = t1018 * t982;
t1084 = (-t544 * t572 - t545 * t570 + t1226) * t1466 + (-t614 * t654 - t615 * t652 + t1217) * t1467 + (-t650 * t691 - t651 * t689 + t1214) * t1468 + (t1021 + (t1496 * t982 - t1497 * t985) * t904) * t1469;
t863 = t904 * t982;
t866 = t904 * t985;
t1085 = (-t553 * t558 - t555 * t557 + t1226) * t1466 + (-t625 * t648 - t627 * t647 + t1217) * t1467 + (-t662 * t677 - t664 * t676 + t1214) * t1468 + (-t1496 * t863 + t1497 * t866 + t1021) * t1469;
t1519 = t1084 - t1085;
t1212 = -t960 * t1267 - t555 * t1281;
t1493 = t1160 * t959;
t976 = t985 * pkin(7);
t1196 = -t982 * (pkin(1) * t982 + t1166 - t976) + t985 * (-t982 * pkin(7) - t957 + (-pkin(1) + t956) * t985);
t1009 = rSges(7,3) * t1281 + t1066;
t1178 = t1160 * t905;
t1082 = t982 * (pkin(4) * t1283 + t1322) + t1178 + t985 * (pkin(4) * t1281 - t961);
t342 = t985 * t1009 + t978 * t909 + t1082 + (t909 * t982 + t641) * t982;
t294 = t342 + t1196;
t182 = t1493 * t294 + t1212;
t1211 = -t960 * t1266 - t572 * t1281;
t212 = t1493 * t342 + t1211;
t1208 = -t960 * t1262 - t627 * t1281;
t446 = t982 * (rSges(6,1) * t1293 + rSges(6,3) * t985 - t938) + t1082 + t985 * (rSges(6,1) * t1292 + t1164);
t394 = t446 + t1196;
t260 = t1493 * t394 + t1208;
t1206 = -t664 * t1281 - t662 * t1283;
t506 = t982 * (t907 * t982 - t975) + t1178 + t985 * (t907 * t985 + t1348);
t433 = t506 + t1196;
t302 = t1493 * t433 + t1206;
t1207 = -t960 * t1261 - t654 * t1281;
t314 = t1493 * t446 + t1207;
t1201 = -t691 * t1281 - t689 * t1283;
t353 = t1493 * t506 + t1201;
t1147 = (t314 + t260) * t1467 + (t353 + t302) * t1468 + (t212 + t182) * t1466;
t1179 = t982 * (qJ(4) * t1283 - t945) + t985 * (-pkin(3) * t1292 + t924);
t1500 = t982 * (rSges(5,3) * t1283 - t937) + t985 * (-rSges(5,1) * t1292 + t941);
t532 = t1179 + t1500;
t1070 = t692 * t1292 + t690 * t1293 - t532 * t960;
t1079 = pkin(4) * t1493;
t1020 = -t1079 + t1179;
t1499 = t977 * t1068 + t985 * t1163;
t495 = t1020 + t1499;
t1071 = t655 * t1292 + t653 * t1293 - t495 * t960;
t699 = t1153 * t982 - t1168;
t700 = t914 + (-t1152 - t1350) * t985;
t1487 = (-pkin(9) * t1292 + t700 + t946) * t985 + (pkin(5) * t1283 + t699 - t944) * t982;
t403 = t1020 + t1487;
t1072 = t573 * t1292 + t571 * t1293 - t403 * t960;
t1491 = (t342 * t959 + t1072 + t1211) * t1466 + (t446 * t959 + t1071 + t1207) * t1467 + (t506 * t959 + t1070 + t1201) * t1468;
t1516 = t1147 - t1491;
t1515 = t1557 * t982;
t1514 = t1557 * t985;
t1513 = (-t1560 + t1561) * t960 + (-t1558 + t1564) * t959;
t1512 = (t771 + t1569) * t960 + (t777 + t1568) * t959;
t1511 = (-t1553 * t985 + t1554 * t982) * t960 + (t1555 * t985 - t1556 * t982) * t959;
t1156 = m(7) / 0.4e1 + m(6) / 0.4e1;
t1294 = t959 * t960;
t1169 = t1160 * t1294;
t1494 = (m(5) / 0.4e1 + t1156) * (t1169 - t1294);
t971 = Icges(3,4) * t984;
t917 = -Icges(3,2) * t981 + t971;
t918 = Icges(3,1) * t981 + t971;
t1203 = -Icges(7,2) * t873 + t638 - t837;
t1205 = Icges(7,1) * t872 - t636 + t838;
t666 = -Icges(7,5) * t872 - Icges(7,6) * t873;
t277 = -t1203 * t872 - t1205 * t873 + t1283 * t666;
t1202 = -Icges(7,2) * t875 + t640 + t839;
t1204 = -Icges(7,1) * t874 + t1331 + t637;
t667 = Icges(7,5) * t874 - Icges(7,6) * t875;
t278 = -t1202 * t872 - t1204 * t873 + t1283 * t667;
t133 = -t277 * t985 + t278 * t982;
t279 = t1203 * t874 - t1205 * t875 + t1281 * t666;
t280 = t1202 * t874 - t1204 * t875 + t1281 * t667;
t134 = -t279 * t985 + t280 * t982;
t1234 = t133 * t1445 + t134 * t1447;
t1488 = t1009 * t982 - t985 * t641;
t1127 = t959 * t433 + t1206;
t1128 = t959 * t394 + t1208;
t1129 = t959 * t294 + t1212;
t1146 = (t1071 + t1128) * t1467 + (t1070 + t1127) * t1468 + (t1072 + t1129) * t1466;
t317 = t667 * t959 + (t1202 * t980 + t1204 * t983) * t960;
t1318 = t317 * t982;
t316 = t666 * t959 + (t1203 * t980 + t1205 * t983) * t960;
t1319 = t316 * t985;
t1192 = t753 + t847;
t1193 = -t1002 - t854;
t391 = -t1192 * t872 - t1193 * t873 + t1283 * t840;
t392 = t1192 * t874 - t1193 * t875 + t1281 * t840;
t1083 = -t1319 / 0.4e1 + t1318 / 0.4e1 + t392 * t1536 + t391 * t1535;
t1250 = t985 * t1552;
t1273 = t982 * t168;
t1481 = t1250 / 0.4e1 + t1273 / 0.4e1 + t1552 * t1535 - t168 * t1536;
t835 = Icges(3,5) * t982 + t919 * t985;
t1173 = -t916 * t985 + t835;
t1257 = t982 * t984;
t1279 = t981 * t982;
t949 = Icges(3,4) * t1279;
t834 = Icges(3,1) * t1257 - Icges(3,5) * t985 - t949;
t1174 = -Icges(3,2) * t1257 + t834 - t949;
t833 = Icges(3,6) * t982 + t917 * t985;
t1175 = -t918 * t985 - t833;
t832 = Icges(3,4) * t1257 - Icges(3,2) * t1279 - Icges(3,6) * t985;
t1176 = t918 * t982 + t832;
t1480 = (-t1173 * t982 + t1174 * t985) * t981 + (t1175 * t982 + t1176 * t985) * t984;
t992 = t750 * t1284 / 0.2e1 - t752 * t1282 / 0.2e1 + t1560 * t1539 + (t1001 + t1564) * t1538 + t1558 * t1537 + (t1025 + t1561) * t1449;
t1475 = 0.2e1 * t1493;
t1474 = 0.4e1 * qJD(1);
t1473 = 2 * qJD(2);
t1472 = 4 * qJD(2);
t1471 = 2 * qJD(3);
t1470 = 4 * qJD(3);
t472 = t1488 * t960;
t536 = -t1009 * t959 + t1281 * t755;
t1130 = -t472 * t403 + t534 * t573 - t536 * t571;
t360 = (t699 * t985 - t700 * t982) * t960 + t1488 * t959;
t424 = (t754 * t982 - t641) * t960 + (-t755 * t982 - t699) * t959;
t425 = ((-t754 + t1349) * t985 + t1066) * t960 + (t755 * t985 + t700) * t959;
t1137 = t360 * t294 - t424 * t555 - t425 * t553;
t1460 = m(7) * (t1130 + t1137);
t38 = t342 * t360 - t424 * t572 - t425 * t570 + t1130;
t1459 = m(7) * t38;
t1093 = t672 * t982 + t985 * t673;
t239 = t294 * t1093;
t268 = t342 * t1093;
t865 = (rSges(7,1) * t980 + rSges(7,2) * t983) * t960;
t1456 = m(7) * (t239 + t268 + ((t555 + t572) * t985 + (t553 + t570) * t982) * t865);
t1227 = t424 * t483 + t425 * t484;
t1455 = m(7) * (t534 * t544 - t536 * t545 + t1227);
t1453 = m(7) * (t534 * t557 - t536 * t558 + t1227);
t1215 = t534 * t1281 - t536 * t1283;
t1451 = m(7) * (-t360 * t960 + (t424 * t985 + t425 * t982 - t472) * t959 + t1215);
t1353 = rSges(3,1) * t984;
t1110 = pkin(1) + t1353;
t1162 = rSges(3,2) * t1279 + t985 * rSges(3,3);
t720 = -t1110 * t982 + t1162 + t976;
t1278 = t981 * t985;
t951 = rSges(3,2) * t1278;
t721 = -t951 + t1110 * t985 + (rSges(3,3) + pkin(7)) * t982;
t920 = rSges(3,1) * t981 + rSges(3,2) * t984;
t883 = t920 * t982;
t884 = t920 * t985;
t1442 = m(3) * (t720 * t883 - t721 * t884);
t611 = t982 * t781 + t985 * (rSges(4,1) * t1281 + t1103);
t491 = t611 + t1196;
t649 = -t982 * t863 - t985 * t866;
t325 = t491 * t649 + (t1496 * t985 + t1497 * t982) * t908;
t1439 = m(4) * t325;
t421 = t1160 * t904 * t908 + t611 * t649;
t1436 = m(4) * t421;
t1435 = m(4) * (-t1496 * t675 + t1497 * t674);
t1434 = m(4) * (t674 * t863 - t675 * t866);
t220 = t433 * t532 - t662 * t690 - t664 * t692;
t1430 = m(5) * t220;
t257 = t506 * t532 - t689 * t690 - t691 * t692;
t1426 = m(5) * t257;
t1425 = m(5) * t302;
t1419 = m(5) * (t609 * t650 + t610 * t651);
t1418 = m(5) * t353;
t1416 = m(5) * (t609 * t676 + t610 * t677);
t1415 = m(5) * (t610 * t1292 - t1293 * t609);
t178 = t394 * t495 - t625 * t653 - t627 * t655;
t1411 = m(6) * t178;
t205 = t446 * t495 - t652 * t653 - t654 * t655;
t1407 = m(6) * t205;
t1406 = m(6) * t260;
t1400 = m(6) * t314;
t1399 = m(6) * (t574 * t614 + t575 * t615);
t1397 = m(6) * (t574 * t647 + t575 * t648);
t1396 = m(6) * (t575 * t1292 - t1293 * t574);
t1390 = m(7) * (-t360 * t472 + t424 * t534 - t425 * t536);
t109 = t294 * t403 - t553 * t571 - t555 * t573;
t1389 = m(7) * t109;
t113 = t342 * t403 - t570 * t571 - t572 * t573;
t1387 = m(7) * t113;
t1022 = t338 * t865;
t1384 = m(7) * (-t553 * t673 + t555 * t672 + t1022);
t1382 = m(7) * (-t570 * t673 + t572 * t672 + t1022);
t1381 = m(7) * t182;
t1376 = m(7) * t212;
t1374 = m(7) * (t483 * t544 + t484 * t545);
t1372 = m(7) * (t483 * t557 + t484 * t558);
t1370 = m(7) * (-t424 * t982 + t425 * t985);
t1369 = m(7) * (-t1493 * t472 + t1215);
t1368 = m(7) * (t484 * t1292 - t1293 * t483);
t1367 = m(7) * (-t536 * t1292 - t1293 * t534);
t401 = -t534 * t985 + t536 * t982;
t1366 = m(7) * t401;
t1361 = m(7) * (-t1093 * t960 - t1493 * t865);
t1045 = t672 * t985 - t673 * t982;
t1360 = m(7) * t1045 * t959;
t1359 = m(7) * t1093;
t1354 = m(7) * qJD(6);
t1150 = -t1370 / 0.2e1;
t258 = t1369 / 0.2e1;
t441 = t1361 / 0.2e1;
t54 = t441 + t258 - t1451 / 0.2e1;
t1341 = t54 * qJD(4) + qJD(5) * t1150;
t1157 = t1466 + t1467;
t1340 = t54 * qJD(6) + (-0.4e1 * t1494 + 0.2e1 * (t1468 + t1157) * (-t1493 * t960 + t1169)) * qJD(4);
t462 = 0.4e1 * t1494;
t91 = t1451 / 0.2e1;
t53 = t441 + t91 - t1369 / 0.2e1;
t1339 = t462 * qJD(4) + t53 * qJD(6);
t1321 = t290 * t985;
t1320 = t291 * t982;
t1316 = t416 * t982;
t1308 = t832 * t981;
t1255 = t984 * t985;
t1149 = t1370 / 0.2e1;
t315 = m(7) * (t571 * t985 - t982 * t573) + m(6) * (t653 * t985 - t982 * t655);
t1229 = t315 * qJD(3) + qJD(6) * t1149;
t1228 = t1158 * t1149;
t830 = Icges(3,5) * t1257 - Icges(3,6) * t1279 - Icges(3,3) * t985;
t1195 = -t834 * t1255 - t982 * t830;
t1056 = Icges(3,5) * t984 - Icges(3,6) * t981;
t831 = Icges(3,3) * t982 + t1056 * t985;
t1194 = t835 * t1255 + t982 * t831;
t546 = t1157 * t1475;
t1159 = t546 * qJD(1);
t1048 = t417 * t985 - t1316;
t193 = t1048 * t960 + t502 * t959;
t999 = t1027 * t960 - t1314;
t262 = -t695 * t872 + t697 * t873 + t982 * t999;
t998 = t1026 * t960 - t1313;
t263 = -t696 * t872 + t698 * t873 + t982 * t998;
t45 = (t262 * t982 + t263 * t985 + t458) * t960 + (-t1050 + t339) * t959;
t264 = t874 * t695 + t875 * t697 + t985 * t999;
t265 = t874 * t696 + t875 * t698 + t985 * t998;
t46 = (t264 * t982 + t265 * t985 + t460) * t960 + (-t1049 + t340) * t959;
t64 = (t290 * t982 + t291 * t985 + t502) * t960 + (-t1048 + t357) * t959;
t14 = t1390 + (t46 * t1444 + t45 * t1447 + t193 / 0.2e1) * t960 + (-t1250 / 0.2e1 - t1273 / 0.2e1 + t64 / 0.2e1) * t959;
t55 = t258 + t91 - t1361 / 0.2e1;
t1154 = t55 * qJD(4) + qJD(5) * t1149 + t14 * qJD(6);
t1148 = t1456 / 0.2e1 + t1234;
t1117 = t1283 / 0.2e1;
t1114 = t1281 / 0.2e1;
t1107 = -t908 - t1357;
t743 = t835 * t1257;
t1090 = t985 * t831 - t743;
t1087 = t833 * t981 - t830;
t1078 = t1160 * t1358;
t1076 = t1170 - t1357;
t1055 = -Icges(3,5) * t981 - Icges(3,6) * t984;
t1031 = t1105 - t1357;
t1023 = t1031 - t902;
t1019 = -t1078 + t1179;
t1017 = t1031 + t1502;
t122 = -t262 * t985 + t263 * t982;
t123 = -t264 * t985 + t265 * t982;
t1014 = (t123 + (-t1515 * t982 + t1511) * t985 + t1514 * t977) * t1447 + (t122 + (-t1514 * t985 + t1511) * t982 + t1515 * t978) * t1445;
t1013 = -t1093 * t472 + t401 * t865;
t1008 = t1481 + t1566;
t1004 = t123 * t1114 + t122 * t1117 + t45 * t1445 + t46 * t1447 + (t1320 - t1321) * t1449 + (t416 * t985 + t417 * t982) * t1537 - t1234 + t1525 * t959;
t1000 = t1116 * t1524 + t1574;
t100 = t392 * t959 + (t279 * t982 + t280 * t985) * t960;
t99 = t391 * t959 + (t277 * t982 + t278 * t985) * t960;
t996 = t133 * t1117 + t134 * t1114 + t193 * t1538 - t45 * t1283 / 0.2e1 - t46 * t1281 / 0.2e1 + t64 * t1539 + t100 * t1447 + t99 * t1445 - t1390 + (t1273 + t1250 + t1318 - t1319) * t1449;
t995 = -t1079 + t1019;
t990 = t1000 + t1008 - t1083;
t989 = t1118 * t1524 + t1008 + t1083 - t1574;
t988 = t1000 + t1083 - t1481 + t1566;
t987 = t1320 / 0.2e1 - t1321 / 0.2e1 + (t1260 + t340 + (-t886 + t889) * t982 + t1513 * t985 + t1556 * t960 + t1554 * t959) * t1447 + (t1513 * t982 + t1553 * t959 + t1555 * t960 - t1305 - t1306 + t1307 + t339) * t1445 - t993;
t986 = -t992 - t1316 / 0.2e1 + t1512 * t1448 + (t416 + t1512) * t1447;
t922 = -rSges(3,2) * t981 + t1353;
t878 = t1055 * t985;
t877 = t1055 * t982;
t759 = t1107 * t985;
t757 = t1107 * t982;
t665 = t1076 * t985;
t663 = t1076 * t982;
t628 = t1023 * t985;
t626 = t1023 * t982;
t602 = -t1078 + t649;
t569 = -t1281 * t865 + t959 * t673;
t568 = t1283 * t865 - t672 * t959;
t567 = -t1278 * t833 + t1194;
t566 = -t1278 * t832 - t1195;
t565 = -t1279 * t833 - t1090;
t556 = t1017 * t985;
t554 = t1017 * t982;
t547 = t1156 * t1475 - (m(6) + m(7)) * t1493 / 0.2e1;
t531 = t1359 / 0.2e1;
t512 = t1045 * t960;
t500 = -t1360 / 0.2e1;
t499 = t1019 + t1500;
t461 = t995 + t1499;
t435 = (t1301 + (t1192 * t980 + t1193 * t983) * t960) * t959;
t427 = -t566 * t985 + t567 * t982;
t426 = -(-(-t834 * t984 + t1308) * t982 - t985 * t830) * t985 + t565 * t982;
t398 = t1366 / 0.2e1;
t379 = t995 + t1487;
t359 = t1367 / 0.2e1;
t252 = qJD(6) * t1150;
t185 = (t565 - t743 + (t831 + t1308) * t985 + t1195) * t985 + t1194 * t982;
t184 = (t1087 * t985 - t1194 + t567) * t985 + (t1087 * t982 + t1090 + t566) * t982;
t176 = t1382 / 0.2e1;
t174 = t440 * t865 + t268;
t171 = t1384 / 0.2e1;
t159 = t1368 + t1396 + t1415;
t140 = t419 * t865 + t239;
t126 = t398 - t1359 / 0.2e1;
t125 = t531 + t398;
t124 = t531 - t1366 / 0.2e1;
t107 = t1219 + t1222;
t105 = t359 + t1360 / 0.2e1;
t104 = t500 + t359;
t103 = t500 - t1367 / 0.2e1;
t89 = t1223 + t1225;
t87 = t1376 + t1400 + t1418;
t82 = t1453 / 0.2e1;
t80 = t1455 / 0.2e1;
t75 = t1381 + t1406 + t1425;
t47 = t1372 + t1397 + t1416 + t992 + t1434;
t37 = t1459 / 0.2e1;
t35 = t1460 / 0.2e1;
t34 = t992 + t1442 + t1435 + t1419 + t1399 + t1374 - t1521 + (t917 / 0.2e1 + t918 / 0.2e1) * t984;
t31 = m(7) * t174 + t1234;
t30 = m(7) * t140 + t1234;
t28 = t1132 + t1133;
t25 = t1134 + t1135;
t22 = t1014 + t1387 + t1407 + t1426 + t1436;
t21 = t22 * qJD(3);
t20 = t1014 + t1389 + t1411 + t1430 + t1439;
t19 = t1491 + t1147 - t1146;
t18 = t1146 + t1516;
t17 = t1146 - t1516;
t16 = t37 - t1460 / 0.2e1 + t1148;
t15 = t35 - t1459 / 0.2e1 + t1148;
t11 = t35 + t37 - t1456 / 0.2e1 + t1004;
t10 = t993 + (t427 / 0.2e1 - t185 / 0.2e1) * t985 + (t426 / 0.2e1 + t184 / 0.2e1) * t982;
t9 = t176 + t82 + t988;
t8 = -t1453 / 0.2e1 + t176 + t989;
t7 = -t1382 / 0.2e1 + t82 + t990;
t6 = t988 + t80 + t171;
t5 = -t1455 / 0.2e1 + t989 + t171;
t4 = -t1384 / 0.2e1 + t990 + t80;
t3 = t993 - t1519;
t2 = t993 + t1519;
t1 = t987 + t1084 + t1085;
t12 = [t34 * qJD(2) + t47 * qJD(3) + t159 * qJD(4) + t228 * qJD(5) + t191 * qJD(6), t34 * qJD(1) + t1 * qJD(3) + t25 * qJD(4) + t89 * qJD(5) + t6 * qJD(6) + (m(3) * ((-t720 * t985 - t721 * t982) * t922 + (-t883 * t985 + t884 * t982) * t920) / 0.2e1 + (t674 * t759 + t675 * t757) * t1469 + (t609 * t665 + t610 * t663 - t650 * t664 - t651 * t662) * t1468 + (t574 * t628 + t575 * t626 - t614 * t627 - t615 * t625) * t1467 + (t483 * t556 + t484 * t554 - t544 * t555 - t545 * t553) * t1466) * t1473 + (t987 + t185 * t1444 + (t1173 * t984 + t1175 * t981) * t1447 + (t426 + t184) * t1448 + (t1174 * t984 - t1176 * t981 + t427) * t1445 + (t978 / 0.2e1 + t977 / 0.2e1) * t1056) * qJD(2), t47 * qJD(1) + t1 * qJD(2) + t987 * qJD(3) + t28 * qJD(4) + t107 * qJD(5) + t9 * qJD(6) + ((t1021 + (-t863 * t985 + t866 * t982) * t904) * t1469 + (-t676 * t691 - t677 * t689 + t1214) * t1468 + (-t647 * t654 - t648 * t652 + t1217) * t1467 + (-t557 * t572 - t558 * t570 + t1226) * t1466) * t1471, qJD(1) * t159 + qJD(2) * t25 + qJD(3) * t28 + qJD(5) * t547 + qJD(6) * t104, qJD(2) * t89 + qJD(3) * t107 + qJD(4) * t547 + qJD(6) * t125 + t1547, t1573 + t6 * qJD(2) + t9 * qJD(3) + t104 * qJD(4) + t125 * qJD(5) + (t435 + m(7) * (t483 * t568 + t484 * t569 - t534 * t672 - t536 * t673) + ((t317 / 0.2e1 + t392 / 0.2e1) * t985 + (t316 / 0.2e1 + t391 / 0.2e1) * t982) * t960) * qJD(6); (t986 + t1521 - (t917 + t918) * t984 / 0.2e1) * qJD(1) + t10 * qJD(2) + t3 * qJD(3) - t24 * qJD(4) - t88 * qJD(5) + t5 * qJD(6) + (-t1442 / 0.4e1 - t1435 / 0.4e1 - t1419 / 0.4e1 - t1399 / 0.4e1 - t1374 / 0.4e1) * t1474, t10 * qJD(1) + (m(4) * (-t1496 * t759 - t1497 * t757 + t491 * t602) + m(3) * ((t982 * (rSges(3,1) * t1257 - t1162) + t985 * (rSges(3,1) * t1255 + t982 * rSges(3,3) - t951)) * (-t982 * t883 - t884 * t985) + t1160 * t922 * t920) + (t977 * t878 + (-t982 * t877 + t1480) * t985) * t1447 + (t978 * t877 + (-t985 * t878 + t1480) * t982) * t1445 + m(7) * (t294 * t379 - t553 * t554 - t555 * t556) + m(6) * (t394 * t461 - t625 * t626 - t627 * t628) + m(5) * (t433 * t499 - t662 * t663 - t664 * t665) + t1014) * qJD(2) + t20 * qJD(3) + t75 * qJD(4) + t30 * qJD(6), t3 * qJD(1) + t20 * qJD(2) + t1014 * qJD(3) + t18 * qJD(4) + t15 * qJD(6) + (-t1387 / 0.4e1 - t1407 / 0.4e1 - t1426 / 0.4e1 - t1436 / 0.4e1) * t1470 + ((t325 + t421) * t1469 + (t113 + t109) * t1466 + (t205 + t178) * t1467 + (t257 + t220) * t1468) * t1471, qJD(2) * t75 + qJD(3) * t18 + t1340 - t1575, t252 - t1565, t5 * qJD(1) + t30 * qJD(2) + t15 * qJD(3) + (m(7) * (t512 * t294 - t569 * t553 - t568 * t555 + t1013) + t996) * qJD(6) + t1341; t986 * qJD(1) + t2 * qJD(2) + t993 * qJD(3) + t1518 * qJD(4) - t106 * qJD(5) + t8 * qJD(6) + (-t1434 / 0.4e1 - t1416 / 0.4e1 - t1397 / 0.4e1 - t1372 / 0.4e1) * t1474, t2 * qJD(1) + t1014 * qJD(2) + t21 + t19 * qJD(4) + t16 * qJD(6) + (-t1389 / 0.4e1 - t1411 / 0.4e1 - t1430 / 0.4e1 - t1439 / 0.4e1) * t1472 + ((t342 * t379 - t554 * t570 - t556 * t572 + t109) * t1466 + (t446 * t461 - t626 * t652 - t628 * t654 + t178) * t1467 + (t499 * t506 - t663 * t689 - t665 * t691 + t220) * t1468 + (t611 * t602 + (-t757 * t982 - t759 * t985) * t904 + t325) * t1469) * t1473, qJD(1) * t993 + qJD(2) * t22 + qJD(4) * t87 + qJD(6) * t31 + t21, qJD(2) * t19 + qJD(3) * t87 + t1340 + t1570, t252 - t1563, t8 * qJD(1) + t16 * qJD(2) + t31 * qJD(3) + (m(7) * (t512 * t342 - t568 * t572 - t569 * t570 + t1013) + t996) * qJD(6) + t1341; t24 * qJD(2) - t1518 * qJD(3) - t546 * qJD(5) + t103 * qJD(6) + (-t1368 / 0.4e1 - t1415 / 0.4e1 - t1396 / 0.4e1) * t1474, t1575 + t17 * qJD(3) + (-t1381 / 0.4e1 - t1406 / 0.4e1 - t1425 / 0.4e1) * t1472 + ((-t960 * t379 + (t554 * t982 + t556 * t985) * t959 + t1129) * t1466 + (-t960 * t461 + (t626 * t982 + t628 * t985) * t959 + t1128) * t1467 + (-t960 * t499 + (t663 * t982 + t665 * t985) * t959 + t1127) * t1468) * t1473 + t1339, -t1570 + t17 * qJD(2) + (-t1376 / 0.4e1 - t1400 / 0.4e1 - t1418 / 0.4e1) * t1470 + t1491 * t1471 + t1339, t1158 * t462, -t1159, t103 * qJD(1) + (-t512 * t960 + (t568 * t985 + t569 * t982) * t959) * t1354 + t1158 * t53; t88 * qJD(2) + t106 * qJD(3) + t546 * qJD(4) + t124 * qJD(6) - t1547, t1565 + ((t554 * t985 - t982 * t556) * t1466 + (t626 * t985 - t982 * t628) * t1467) * t1473 + t1229, qJD(2) * t315 + t1229 + t1563, t1159, 0, t124 * qJD(1) + (-t568 * t982 + t569 * t985) * t1354 + t1228; t4 * qJD(2) + t7 * qJD(3) + t105 * qJD(4) + t126 * qJD(5) - t1573, t4 * qJD(1) + ((-t379 * t472 + t534 * t556 - t536 * t554 + t1137 - t140) * m(7) + t1004) * qJD(2) + t11 * qJD(3) + t1154, t7 * qJD(1) + t11 * qJD(2) + ((t38 - t174) * m(7) + t1004) * qJD(3) + t1154, qJD(1) * t105 + t1158 * t55, qJD(1) * t126 + t1228 (t435 * t1449 + m(7) * (-t472 * t512 + t534 * t568 - t536 * t569) + (t100 * t1444 + t99 * t1447 + (t316 * t982 + t317 * t985) * t1449) * t960) * qJD(6) + t1158 * t14;];
Cq  = t12;