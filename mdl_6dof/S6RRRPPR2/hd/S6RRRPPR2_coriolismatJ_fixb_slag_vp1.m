% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:59
% EndTime: 2019-03-09 15:26:10
% DurationCPUTime: 60.03s
% Computational Cost: add. (152603->1119), mult. (138210->1433), div. (0->0), fcn. (143694->10), ass. (0->666)
t1589 = -Icges(6,4) + Icges(5,5);
t1588 = Icges(6,5) - Icges(5,6);
t979 = qJ(2) + qJ(3);
t959 = pkin(10) + t979;
t955 = sin(t959);
t956 = cos(t959);
t960 = sin(t979);
t961 = cos(t979);
t1591 = Icges(4,5) * t961 - Icges(4,6) * t960 + t1588 * t955 + t1589 * t956;
t985 = cos(qJ(1));
t1289 = t956 * t985;
t982 = sin(qJ(1));
t1290 = t956 * t982;
t1297 = t955 * t982;
t1104 = rSges(6,1) * t985 - rSges(6,3) * t1297;
t986 = -pkin(8) - pkin(7);
t1152 = -qJ(4) + t986;
t1373 = pkin(3) * t961;
t984 = cos(qJ(2));
t975 = t984 * pkin(2);
t957 = t975 + pkin(1);
t911 = t957 + t1373;
t1163 = -t985 * t1152 - t982 * t911;
t1335 = qJ(5) * t955;
t588 = (-t1335 + (rSges(6,2) - pkin(4)) * t956) * t982 + t1104 + t1163;
t1103 = rSges(6,1) * t982 - rSges(6,2) * t1289;
t1357 = rSges(6,3) + qJ(5);
t1371 = pkin(4) * t956;
t952 = t982 * t1152;
t589 = -t952 + (t1357 * t955 + t1371 + t911) * t985 + t1103;
t1209 = t588 * t1289 + t589 * t1290;
t1456 = -rSges(7,3) - pkin(4);
t983 = cos(qJ(6));
t1263 = t982 * t983;
t980 = sin(qJ(6));
t1279 = t980 * t985;
t864 = t1263 * t955 + t1279;
t1261 = t983 * t985;
t1280 = t980 * t982;
t865 = -t1280 * t955 + t1261;
t1507 = t865 * rSges(7,1) - t864 * rSges(7,2);
t931 = pkin(9) * t1290;
t890 = pkin(5) * t985 - t931;
t1527 = t1507 + (t1456 * t956 - t1335) * t982 + t890 + t1163;
t862 = t1261 * t955 - t1280;
t863 = t1279 * t955 + t1263;
t1067 = rSges(7,1) * t863 + rSges(7,2) * t862;
t1153 = pkin(9) - t1456;
t1370 = pkin(5) * t982;
t483 = t1370 - t952 + (t1153 * t956 + t1335 + t911) * t985 + t1067;
t1219 = t1289 * t1527 + t483 * t1290;
t1296 = t955 * t985;
t1479 = m(7) / 0.2e1;
t1481 = m(6) / 0.2e1;
t1374 = pkin(3) * t960;
t1372 = pkin(4) * t955;
t882 = -qJ(5) * t956 + t1372;
t1109 = t882 + t1374;
t1369 = pkin(9) * t955;
t1065 = rSges(7,1) * t980 + rSges(7,2) * t983;
t726 = rSges(7,3) * t955 - t1065 * t956;
t1020 = t1109 + t726 + t1369;
t981 = sin(qJ(2));
t1375 = pkin(2) * t981;
t1007 = t1020 + t1375;
t548 = t1007 * t982;
t550 = t1007 * t985;
t1363 = rSges(6,2) * t955;
t1064 = rSges(6,3) * t956 + t1363;
t1070 = t1374 + t1375;
t1023 = t1070 + t882 - t1064;
t620 = t1023 * t982;
t622 = t1023 * t985;
t1234 = (-t548 * t1296 + t1297 * t550 + t1219) * t1479 + (-t620 * t1296 + t1297 * t622 + t1209) * t1481;
t1083 = t1153 * t955;
t1162 = (rSges(7,1) * t1279 + rSges(7,2) * t1261) * t956;
t913 = qJ(5) * t1289;
t1125 = t913 + t1162;
t900 = t985 * t1070;
t535 = -t1083 * t985 + t1125 - t900;
t1000 = t1369 + (-qJ(5) - t1065) * t956;
t921 = rSges(7,3) * t1297;
t930 = pkin(4) * t1297;
t1161 = t921 + t930;
t536 = (t1070 + t1000) * t982 + t1161;
t1074 = -pkin(4) * t1296 + t913;
t1159 = rSges(6,2) * t1296 + rSges(6,3) * t1289;
t618 = t1074 - t900 + t1159;
t1024 = -t1357 * t956 - t1363;
t619 = t930 + (t1070 + t1024) * t982;
t1235 = ((t535 * t982 + t536 * t985) * t955 + t1219) * t1479 + ((t618 * t982 + t619 * t985) * t955 + t1209) * t1481;
t23 = t1235 - t1234;
t1590 = t23 * qJD(1);
t1587 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t1462 = t955 / 0.2e1;
t1545 = t956 / 0.2e1;
t1050 = Icges(7,5) * t980 + Icges(7,6) * t983;
t1002 = -Icges(7,3) * t955 + t1050 * t956;
t1346 = Icges(7,4) * t863;
t633 = Icges(7,2) * t862 + Icges(7,6) * t1289 + t1346;
t835 = Icges(7,4) * t862;
t636 = Icges(7,1) * t863 + Icges(7,5) * t1289 + t835;
t1045 = -t633 * t983 - t636 * t980;
t1027 = t1002 * t985 - t1045;
t630 = Icges(7,5) * t863 + Icges(7,6) * t862 + Icges(7,3) * t1289;
t1345 = Icges(7,4) * t980;
t1055 = Icges(7,2) * t983 + t1345;
t1003 = -Icges(7,6) * t955 + t1055 * t956;
t687 = t1003 * t985;
t1344 = Icges(7,4) * t983;
t1057 = Icges(7,1) * t980 + t1344;
t1004 = -Icges(7,5) * t955 + t1057 * t956;
t689 = t1004 * t985;
t277 = (-t687 * t983 - t689 * t980 + t630) * t956 + t1027 * t955;
t715 = Icges(7,6) * t956 + t1055 * t955;
t717 = Icges(7,5) * t956 + t1057 * t955;
t1042 = t1003 * t983 + t1004 * t980;
t713 = Icges(7,3) * t956 + t1050 * t955;
t1025 = -t1042 + t713;
t1321 = t1002 * t955;
t997 = t1025 * t956 + t1321;
t324 = t715 * t862 + t717 * t863 + t985 * t997;
t1549 = -t715 * t983 - t717 * t980 - t1002;
t338 = t1025 * t955 + t1549 * t956;
t1323 = t630 * t955;
t394 = t1045 * t956 + t1323;
t837 = Icges(7,4) * t865;
t634 = -Icges(7,2) * t864 - Icges(7,6) * t1290 + t837;
t836 = Icges(7,4) * t864;
t638 = -Icges(7,1) * t865 + Icges(7,5) * t1290 + t836;
t1044 = t634 * t983 - t638 * t980;
t632 = -Icges(7,5) * t865 + Icges(7,6) * t864 + Icges(7,3) * t1290;
t1322 = t632 * t955;
t396 = t1044 * t956 + t1322;
t443 = -t1002 * t1289 - t1003 * t862 - t1004 * t863;
t445 = -t1002 * t1290 - t1003 * t864 + t1004 * t865;
t471 = t1042 * t956 - t1321;
t1586 = t1462 * t338 + t1545 * t471 - (t394 + t443) * t1296 / 0.4e1 - (t396 + t445) * t1297 / 0.4e1 + t1289 * (t277 + t324) / 0.4e1;
t885 = t1335 + t1371;
t1107 = -t885 - t1373;
t1362 = rSges(7,3) * t956;
t725 = t1065 * t955 + t1362;
t1073 = t1107 - t725;
t573 = t1073 * t982 - t931;
t1021 = -pkin(9) * t956 + t1073;
t575 = t1021 * t985;
t1221 = -t548 * t573 - t550 * t575;
t1241 = t985 * t986;
t1368 = pkin(1) - t957;
t958 = t982 * t986;
t976 = t985 * pkin(7);
t1187 = -t982 * (t1368 * t982 - t1241 - t976) + t985 * (-pkin(7) * t982 - t1368 * t985 - t958);
t1010 = rSges(7,3) * t1289 + t1067;
t977 = t982 ^ 2;
t978 = t985 ^ 2;
t1156 = t977 + t978;
t1197 = -t982 * (t982 * t957 + t1163 + t1241) + t985 * (-t952 + t958 + (t911 - t957) * t985);
t1078 = t1156 * t885 + t1197;
t641 = rSges(7,3) * t1290 - t1507;
t308 = t1078 + (pkin(9) * t1289 + t1010 + t1370) * t985 + (t641 - t890) * t982;
t264 = t308 + t1187;
t1175 = t982 * (qJ(5) * t1290 - t930) + t985 * t1074;
t1502 = t1156 * t955;
t693 = t1065 * t1290 - t921;
t694 = -rSges(7,3) * t1296 + t1162;
t1008 = -pkin(9) * t1502 + t982 * t693 + t985 * t694 + t1175;
t1076 = t1156 * t1374;
t409 = -t1076 + t1008;
t110 = t264 * t409 + t1221;
t886 = -rSges(6,2) * t956 + rSges(6,3) * t955;
t1071 = t1107 - t886;
t646 = t1071 * t982;
t648 = t1071 * t985;
t1215 = -t620 * t646 - t622 * t648;
t388 = -t982 * (rSges(6,2) * t1290 + t1104) + t1078 + t985 * (rSges(6,3) * t1296 + t1103);
t331 = t388 + t1187;
t1077 = t1064 * t977 + t985 * t1159 + t1175;
t470 = -t1076 + t1077;
t165 = t331 * t470 + t1215;
t884 = rSges(5,1) * t955 + rSges(5,2) * t956;
t1032 = t1070 + t884;
t695 = t1032 * t982;
t697 = t1032 * t985;
t1364 = rSges(5,1) * t956;
t887 = -rSges(5,2) * t955 + t1364;
t1106 = -t887 - t1373;
t736 = t1106 * t982;
t738 = t1106 * t985;
t1207 = -t695 * t736 - t697 * t738;
t1101 = -rSges(5,2) * t1296 + rSges(5,3) * t982;
t763 = rSges(5,1) * t1290 - rSges(5,2) * t1297 - t985 * rSges(5,3);
t441 = t982 * t763 + t1197 + t985 * (rSges(5,1) * t1289 + t1101);
t376 = t441 + t1187;
t1160 = rSges(5,1) * t1297 + rSges(5,2) * t1290;
t1174 = -t982 * t1160 - t978 * t884;
t590 = -t1076 + t1174;
t234 = t376 * t590 + t1207;
t906 = rSges(4,1) * t960 + rSges(4,2) * t961;
t1016 = t906 + t1375;
t1505 = t1016 * t985;
t1506 = t1016 * t982;
t1285 = t960 * t985;
t1102 = -rSges(4,2) * t1285 + rSges(4,3) * t982;
t1286 = t960 * t982;
t1158 = rSges(4,2) * t1286 + t985 * rSges(4,3);
t1281 = t961 * t985;
t1282 = t961 * t982;
t628 = t982 * (rSges(4,1) * t1282 - t1158) + t985 * (rSges(4,1) * t1281 + t1102);
t490 = t628 + t1187;
t866 = t906 * t982;
t867 = t906 * t985;
t655 = -t982 * t866 - t985 * t867;
t1365 = rSges(4,1) * t961;
t907 = -rSges(4,2) * t960 + t1365;
t337 = t490 * t655 + (t1505 * t985 + t1506 * t982) * t907;
t1585 = m(4) * t337 + m(5) * t234 + m(6) * t165 + m(7) * t110;
t1584 = t1591 * t985;
t1338 = Icges(6,6) * t956;
t946 = Icges(5,4) * t956;
t1583 = t1338 + t946 + (Icges(5,1) + Icges(6,2)) * t955;
t572 = t1020 * t982;
t574 = t1020 * t985;
t1072 = t1109 - t1064;
t645 = t1072 * t982;
t647 = t1072 * t985;
t1228 = (-t572 * t1296 + t1297 * t574 + t1219) * t1479 + (-t645 * t1296 + t1297 * t647 + t1209) * t1481;
t561 = (t1000 + t1374) * t982 + t1161;
t562 = (-t1083 - t1374) * t985 + t1125;
t643 = t930 + (t1024 + t1374) * t982;
t644 = t913 + (-t1372 - t1374) * t985 + t1159;
t1233 = ((t561 * t985 + t562 * t982) * t955 + t1219) * t1479 + ((t643 * t985 + t644 * t982) * t955 + t1209) * t1481;
t1525 = t1228 - t1233;
t1582 = t1525 * qJD(1);
t369 = t630 * t1290 + t864 * t633 - t865 * t636;
t370 = t1290 * t632 - t634 * t864 - t638 * t865;
t1047 = t369 * t985 + t370 * t982;
t164 = t1047 * t956 + t445 * t955;
t1564 = Icges(4,5) * t1282 - Icges(4,6) * t1286 + t1290 * t1589 + t1297 * t1588 - t1587 * t985;
t1581 = t1587 * t982 + t1584;
t749 = Icges(5,4) * t1290 - Icges(5,2) * t1297 - Icges(5,6) * t985;
t754 = Icges(6,5) * t985 + Icges(6,6) * t1290 - Icges(6,3) * t1297;
t1580 = t749 + t754;
t920 = Icges(5,4) * t1297;
t751 = Icges(5,1) * t1290 - Icges(5,5) * t985 - t920;
t915 = Icges(6,6) * t1297;
t756 = Icges(6,4) * t985 + Icges(6,2) * t1290 - t915;
t1579 = t751 + t756;
t820 = (-Icges(7,5) * t983 + Icges(7,6) * t980) * t956;
t1303 = t955 * t820;
t662 = rSges(7,1) * t862 - rSges(7,2) * t863;
t663 = rSges(7,1) * t864 + rSges(7,2) * t865;
t823 = (Icges(7,2) * t980 - t1344) * t956;
t828 = (-Icges(7,1) * t983 + t1345) * t956;
t191 = (-(t828 / 0.2e1 + t1003 / 0.2e1) * t980 - (-t1004 / 0.2e1 + t823 / 0.2e1) * t983) * t956 + m(7) * (-t1527 * t663 + t483 * t662) + t1303 / 0.2e1;
t1577 = t191 * qJD(1);
t1347 = Icges(5,4) * t955;
t881 = Icges(5,1) * t956 - t1347;
t752 = Icges(5,5) * t982 + t881 * t985;
t873 = Icges(6,3) * t955 - t1338;
t753 = Icges(6,5) * t982 + t873 * t985;
t1348 = Icges(4,4) * t960;
t905 = Icges(4,1) * t961 - t1348;
t812 = Icges(4,5) * t982 + t905 * t985;
t1576 = -t812 * t1282 - t752 * t1290 - t753 * t1297;
t217 = t369 * t982 - t370 * t985;
t1560 = t726 * t1290 - t641 * t955;
t514 = -t1010 * t955 + t1289 * t726;
t1508 = -t1560 * t985 + t514 * t982;
t1118 = t1290 / 0.4e1;
t1120 = -t1290 / 0.4e1;
t367 = t630 * t1289 + t862 * t633 + t863 * t636;
t368 = t632 * t1289 - t634 * t862 + t863 * t638;
t1558 = t367 * t982 - t368 * t985;
t1575 = (t1118 + t1120) * t1558;
t1461 = -t982 / 0.2e1;
t1458 = -t985 / 0.2e1;
t1482 = m(5) / 0.2e1;
t563 = t982 * t572;
t1500 = t574 * t985 + t563;
t1108 = t884 + t1374;
t1068 = t1108 * t985;
t1504 = t1068 * t985;
t1548 = -m(7) / 0.2e1;
t639 = t982 * t645;
t735 = t1108 * t982;
t1129 = t1500 * t1548 + (-t647 * t985 - t639) * t1481 + (-t982 * t735 - t1504) * t1482;
t711 = pkin(3) * t1286 + t1160;
t1131 = (t561 * t982 - t562 * t985) * t1479 + (t643 * t982 - t644 * t985) * t1481 + (t711 * t982 + t1504) * t1482;
t81 = t1131 - t1129;
t1574 = qJD(1) * t81;
t537 = t982 * t548;
t1501 = t550 * t985 + t537;
t607 = t982 * t620;
t1132 = t1501 * t1548 + (-t622 * t985 - t607) * t1481 + (-t982 * t695 - t697 * t985) * t1482;
t699 = t1070 * t982 + t1160;
t700 = -t884 * t985 - t900;
t1134 = (-t535 * t985 + t536 * t982) * t1479 + (-t618 * t985 + t619 * t982) * t1481 + (t699 * t982 - t700 * t985) * t1482;
t78 = t1134 - t1132;
t1573 = t78 * qJD(1);
t1339 = Icges(6,6) * t955;
t1572 = -Icges(6,3) * t956 - t1339 + t881;
t878 = Icges(5,2) * t956 + t1347;
t1571 = -Icges(6,2) * t956 + t1339 + t878;
t229 = m(6) * (t589 * t1296 - t1297 * t588) + m(7) * (t483 * t1296 - t1297 * t1527);
t1570 = qJD(1) * t229;
t879 = -Icges(5,2) * t955 + t946;
t1569 = t879 + t1583;
t1568 = t1579 - t915 - t920 + (-Icges(5,2) - Icges(6,3)) * t1290;
t916 = Icges(6,6) * t1296;
t755 = Icges(6,4) * t982 - Icges(6,2) * t1289 + t916;
t1567 = Icges(6,3) * t1289 + t878 * t985 - t752 + t755 + t916;
t1566 = t1583 * t982 + t1580;
t750 = Icges(5,6) * t982 + t879 * t985;
t1565 = -t1583 * t985 - t750 + t753;
t1457 = t985 / 0.2e1;
t1460 = t982 / 0.2e1;
t1533 = t1458 * t1558 + t1461 * t217;
t809 = Icges(4,4) * t1282 - Icges(4,2) * t1286 - Icges(4,6) * t985;
t1550 = t749 * t955 - t756 * t956 + t809 * t960;
t1551 = t812 * t1281 + t752 * t1289 + t753 * t1296 + t1581 * t982;
t954 = Icges(4,4) * t961;
t903 = -Icges(4,2) * t960 + t954;
t810 = Icges(4,6) * t982 + t903 * t985;
t1552 = t750 * t955 + t755 * t956 + t810 * t960 - t1564;
t1553 = -t1285 * t810 - t1289 * t755 - t1296 * t750 + t1551;
t943 = Icges(4,4) * t1286;
t811 = Icges(4,1) * t1282 - Icges(4,5) * t985 - t943;
t1561 = -t811 * t1281 - t751 * t1289 + t754 * t1296 - t1564 * t982;
t1554 = t1285 * t809 - t1289 * t756 + t1296 * t749 + t1561;
t1562 = t1581 * t985 + t1576;
t1555 = -t1286 * t810 - t1290 * t755 - t1297 * t750 - t1562;
t994 = (t1553 * t982 + t1554 * t985) * t1457 + (t1558 + t1551 * t982 + ((t1550 + t1581) * t985 + t1555 + t1561 + t1576) * t985) * t1458 + (-t217 + (t1552 * t982 - t1554 + t1555 + t1562) * t982 + ((t1552 + t1564) * t985 + (-t751 * t956 + t754 * t955 - t811 * t961 + t1550) * t982 - t1551 + t1553) * t985) * t1460 - t1533;
t1563 = Icges(4,5) * t960 + Icges(4,6) * t961 - t1588 * t956 + t1589 * t955;
t1048 = t367 * t985 + t368 * t982;
t1559 = t1048 * t956 + t443 * t955;
t904 = Icges(4,1) * t960 + t954;
t1556 = t903 + t904;
t1509 = -t1527 * t985 - t483 * t982;
t1483 = m(4) / 0.2e1;
t1547 = -t955 / 0.2e1;
t1546 = -t956 / 0.2e1;
t1543 = t982 / 0.4e1;
t1542 = -t985 / 0.4e1;
t1026 = t1002 * t982 - t1044;
t686 = t1003 * t982;
t688 = t1004 * t982;
t278 = (-t686 * t983 - t688 * t980 + t632) * t956 + t1026 * t955;
t323 = t715 * t864 - t717 * t865 + t982 * t997;
t1531 = t278 + t323;
t1349 = Icges(3,4) * t981;
t933 = Icges(3,2) * t984 + t1349;
t936 = Icges(3,1) * t984 - t1349;
t1529 = (t936 / 0.2e1 - t933 / 0.2e1) * t981;
t1154 = qJD(2) + qJD(3);
t1105 = t957 + t1365;
t690 = -t1105 * t982 + t1158 - t1241;
t691 = t1105 * t985 + t1102 - t958;
t1018 = (-t690 * t985 - t691 * t982) * t907;
t649 = -t763 + t1163;
t650 = -t952 + (t911 + t1364) * t985 + t1101;
t1208 = t738 * t649 + t736 * t650;
t1218 = t648 * t588 + t646 * t589;
t1223 = t1527 * t575 + t573 * t483;
t1081 = (-t535 * t572 - t536 * t574 + t1223) * t1479 + (-t618 * t645 - t619 * t647 + t1218) * t1481 + (-t1068 * t699 - t700 * t735 + t1208) * t1482 + (t1018 + (t1505 * t982 - t1506 * t985) * t906) * t1483;
t1082 = (-t548 * t562 - t550 * t561 + t1223) * t1479 + (-t620 * t644 - t622 * t643 + t1218) * t1481 + (t1068 * t695 - t697 * t711 + t1208) * t1482 + (-t1505 * t866 + t1506 * t867 + t1018) * t1483;
t1526 = t1081 - t1082;
t1206 = -t622 * t1289 - t956 * t607;
t1135 = t955 * t331 + t1206;
t1212 = -t550 * t1289 - t956 * t537;
t1136 = t955 * t264 + t1212;
t1204 = t648 * t1296 + t646 * t1297;
t1211 = t575 * t1296 + t573 * t1297;
t1327 = t470 * t956;
t1328 = t409 * t956;
t1240 = (t1136 + t1211 - t1328) * t1479 + (t1135 + t1204 - t1327) * t1481;
t179 = t1502 * t264 + t1212;
t1210 = -t574 * t1289 - t956 * t563;
t192 = t1502 * t308 + t1210;
t233 = t1502 * t331 + t1206;
t1203 = -t647 * t1289 - t956 * t639;
t274 = t1502 * t388 + t1203;
t1356 = (t274 + t233) * t1481 + (t192 + t179) * t1479;
t1523 = t1240 - t1356;
t1522 = t1563 * t982;
t1521 = t1563 * t985;
t902 = Icges(4,2) * t961 + t1348;
t1520 = (-t902 + t905) * t961 - t1556 * t960 + (-t1571 + t1572) * t956 + (t873 - t1569) * t955;
t1519 = t1579 * t955 + t1580 * t956 + t809 * t961 + t811 * t960;
t1170 = -t902 * t985 + t812;
t1171 = -Icges(4,2) * t1282 + t811 - t943;
t1172 = -t904 * t985 - t810;
t1173 = t904 * t982 + t809;
t1518 = (-t1170 * t982 + t1171 * t985) * t960 + (t1172 * t982 + t1173 * t985) * t961 + (t1565 * t982 + t1566 * t985) * t956 + (t1567 * t982 + t1568 * t985) * t955;
t1516 = -0.2e1 * t1502;
t1478 = m(7) / 0.4e1;
t1480 = m(6) / 0.4e1;
t1151 = t1480 + t1478;
t1298 = t955 * t956;
t1164 = t1156 * t1298;
t1503 = t1151 * (t1164 - t1298);
t971 = Icges(3,4) * t984;
t934 = -Icges(3,2) * t981 + t971;
t935 = Icges(3,1) * t981 + t971;
t1200 = -Icges(7,2) * t863 + t636 + t835;
t1202 = -Icges(7,1) * t862 + t1346 + t633;
t656 = Icges(7,5) * t862 - Icges(7,6) * t863;
t280 = t1200 * t862 - t1202 * t863 + t1289 * t656;
t1199 = Icges(7,2) * t865 + t638 + t836;
t1201 = -Icges(7,1) * t864 - t634 - t837;
t657 = Icges(7,5) * t864 + Icges(7,6) * t865;
t281 = t1199 * t862 - t1201 * t863 + t1289 * t657;
t136 = t280 * t982 - t281 * t985;
t282 = t1200 * t864 + t1202 * t865 + t1290 * t656;
t283 = t1199 * t864 + t1201 * t865 + t1290 * t657;
t137 = t282 * t982 - t283 * t985;
t1238 = t136 * t1460 + t137 * t1458;
t1061 = t955 * t388 + t1203 + t1204;
t1062 = t955 * t308 + t1210 + t1211;
t1277 = t981 * t985;
t1176 = t977 * (-t1070 + t1375) + t985 * (pkin(2) * t1277 - t900);
t372 = t1008 + t1176;
t430 = t1077 + t1176;
t1239 = (-t372 * t956 + t1062) * t1479 + (-t430 * t956 + t1061) * t1481;
t1497 = t1010 * t982 - t985 * t641;
t313 = t657 * t955 + (-t1199 * t983 + t1201 * t980) * t956;
t1331 = t313 * t985;
t312 = t656 * t955 + (-t1200 * t983 + t1202 * t980) * t956;
t1332 = t312 * t982;
t1189 = -t1004 + t823;
t1190 = -t1003 - t828;
t358 = t1189 * t862 - t1190 * t863 + t1289 * t820;
t359 = t1189 * t864 + t1190 * t865 + t1290 * t820;
t1080 = t1332 / 0.4e1 - t1331 / 0.4e1 + t358 * t1543 + t359 * t1542;
t1256 = t985 * t1559;
t1273 = t982 * t164;
t1494 = t1256 / 0.4e1 + t1273 / 0.4e1 + t1559 * t1542 - t164 * t1543;
t853 = Icges(3,5) * t982 + t936 * t985;
t1165 = -t933 * t985 + t853;
t1262 = t982 * t984;
t1278 = t981 * t982;
t949 = Icges(3,4) * t1278;
t852 = Icges(3,1) * t1262 - Icges(3,5) * t985 - t949;
t1166 = -Icges(3,2) * t1262 + t852 - t949;
t851 = Icges(3,6) * t982 + t934 * t985;
t1167 = -t935 * t985 - t851;
t850 = Icges(3,4) * t1262 - Icges(3,2) * t1278 - Icges(3,6) * t985;
t1168 = t935 * t982 + t850;
t1493 = (-t1165 * t982 + t1166 * t985) * t981 + (t1167 * t982 + t1168 * t985) * t984;
t993 = (-t902 / 0.2e1 + t905 / 0.2e1) * t960 + t1556 * t961 / 0.2e1 + t1569 * t1545 + (t713 + t1572) * t1462 + (t1042 + t1571) * t1547 + (t873 - t1549) * t1546;
t1488 = 0.4e1 * qJD(1);
t1487 = 2 * qJD(2);
t1485 = 2 * qJD(3);
t1484 = 4 * qJD(3);
t346 = (t693 * t985 - t694 * t982) * t956 + t1497 * t955;
t399 = (t725 * t982 - t641) * t956 + (-t726 * t982 - t693) * t955;
t400 = ((-t725 + t1362) * t985 + t1067) * t956 + (t726 * t985 + t694) * t955;
t1137 = t346 * t264 - t399 * t550 - t400 * t548;
t1222 = t1560 * t575 - t514 * t573;
t459 = t1497 * t956;
t291 = t459 * t409;
t1474 = m(7) * (-t291 + t1137 + t1222);
t1063 = t346 * t308 - t399 * t574 - t400 * t572 + t1222;
t1473 = m(7) * (-t372 * t459 + t1063);
t522 = t662 * t985 + t663 * t982;
t209 = t264 * t522;
t250 = t308 * t522;
t834 = (-rSges(7,1) * t983 + rSges(7,2) * t980) * t956;
t1470 = m(7) * (t209 + t250 + ((t550 + t574) * t985 + (t548 + t572) * t982) * t834);
t1216 = t1289 * t1560 - t514 * t1290;
t1467 = m(7) * (-t346 * t956 + (t399 * t985 + t400 * t982 - t459) * t955 + t1216);
t1225 = t1527 * t399 + t400 * t483;
t1466 = m(7) * (t1560 * t536 - t514 * t535 + t1225);
t1465 = m(7) * (t1560 * t561 - t514 * t562 + t1225);
t1464 = m(7) * (t1560 * t399 - t346 * t459 - t400 * t514);
t1366 = rSges(3,1) * t984;
t1112 = pkin(1) + t1366;
t1157 = rSges(3,2) * t1278 + t985 * rSges(3,3);
t739 = -t1112 * t982 + t1157 + t976;
t951 = rSges(3,2) * t1277;
t740 = -t951 + t1112 * t985 + (rSges(3,3) + pkin(7)) * t982;
t937 = rSges(3,1) * t981 + rSges(3,2) * t984;
t898 = t937 * t982;
t899 = t937 * t985;
t1455 = m(3) * (t739 * t898 - t740 * t899);
t419 = t1156 * t906 * t907 + t628 * t655;
t415 = m(4) * t419;
t1449 = m(4) * (-t1505 * t691 + t1506 * t690);
t1448 = m(4) * (t690 * t866 - t691 * t867);
t1205 = -t1068 * t738 - t735 * t736;
t510 = t1174 + t1176;
t1444 = m(5) * (t441 * t510 + t1205);
t1441 = m(5) * (t649 * t699 + t650 * t700);
t1440 = m(5) * (-t1068 * t650 + t649 * t711);
t1439 = m(5) * (t649 * t985 + t650 * t982);
t1213 = -t645 * t646 - t647 * t648;
t1426 = m(6) * (t388 * t430 + t1213);
t1417 = m(6) * (t588 * t619 + t589 * t618);
t1416 = m(6) * (t588 * t643 + t589 * t644);
t1415 = m(6) * (t588 * t985 + t589 * t982);
t1220 = -t572 * t573 - t574 * t575;
t1404 = m(7) * (t308 * t372 + t1220);
t1019 = t1509 * t834;
t1399 = m(7) * (-t548 * t662 + t550 * t663 + t1019);
t1398 = m(7) * (-t572 * t662 + t574 * t663 + t1019);
t1391 = m(7) * (t1527 * t536 + t483 * t535);
t1390 = m(7) * (t1527 * t561 + t483 * t562);
t1389 = m(7) * (t399 * t982 - t400 * t985);
t1388 = m(7) * (-t1502 * t459 + t1216);
t1387 = m(7) * t1509;
t1386 = m(7) * (-t514 * t1296 - t1297 * t1560);
t1385 = m(7) * t1508;
t1380 = m(7) * (-t1502 * t834 - t522 * t956);
t1043 = t662 * t982 - t663 * t985;
t1377 = m(7) * t1043 * t955;
t1376 = m(7) * t522;
t1367 = m(7) * qJD(6);
t1146 = -t1389 / 0.2e1;
t252 = t1388 / 0.2e1;
t413 = t1380 / 0.2e1;
t53 = t413 + t252 - t1467 / 0.2e1;
t1355 = qJD(4) * t1146 + t53 * qJD(5);
t1150 = t1479 + t1481;
t1354 = t53 * qJD(6) + (-0.4e1 * t1503 + 0.2e1 * t1150 * (-t1502 * t956 + t1164)) * qJD(5);
t515 = 0.4e1 * t1503;
t87 = t1467 / 0.2e1;
t52 = t413 + t87 - t1388 / 0.2e1;
t1353 = t515 * qJD(5) + t52 * qJD(6);
t1334 = t277 * t982;
t1333 = t278 * t985;
t1329 = t396 * t982;
t1310 = t850 * t981;
t1260 = t984 * t985;
t1145 = t1389 / 0.2e1;
t263 = m(7) * (-t573 * t985 + t575 * t982) + m(6) * (-t646 * t985 + t648 * t982) + m(5) * (-t736 * t985 + t738 * t982);
t1227 = t263 * qJD(3) + qJD(6) * t1145;
t1226 = t1154 * t1145;
t848 = Icges(3,5) * t1262 - Icges(3,6) * t1278 - Icges(3,3) * t985;
t1178 = -t852 * t1260 - t982 * t848;
t1054 = Icges(3,5) * t984 - Icges(3,6) * t981;
t849 = Icges(3,3) * t982 + t1054 * t985;
t1177 = t853 * t1260 + t982 * t849;
t520 = t1150 * t1516;
t1155 = t520 * qJD(1);
t1046 = t394 * t985 + t1329;
t185 = t1046 * t956 + t471 * t955;
t999 = t1027 * t956 - t1323;
t257 = t687 * t864 - t689 * t865 + t982 * t999;
t998 = t1026 * t956 - t1322;
t258 = t686 * t864 - t688 * t865 + t982 * t998;
t44 = (t257 * t985 + t258 * t982 + t445) * t956 + (-t1047 + t323) * t955;
t259 = t687 * t862 + t689 * t863 + t985 * t999;
t260 = t686 * t862 + t688 * t863 + t985 * t998;
t45 = (t259 * t985 + t260 * t982 + t443) * t956 + (-t1048 + t324) * t955;
t63 = (t277 * t985 + t278 * t982 + t471) * t956 + (-t1046 + t338) * t955;
t14 = t1464 + (t45 * t1457 + t44 * t1460 + t185 / 0.2e1) * t956 + (-t1256 / 0.2e1 - t1273 / 0.2e1 + t63 / 0.2e1) * t955;
t54 = t252 + t87 - t1380 / 0.2e1;
t1148 = qJD(4) * t1145 + t54 * qJD(5) + t14 * qJD(6);
t1144 = t1470 / 0.2e1 + t1238;
t1128 = t308 * t409 + t1220;
t1127 = t388 * t470 + t1213;
t1126 = t441 * t590 + t1205;
t1119 = t1290 / 0.2e1;
t1116 = t1289 / 0.2e1;
t1110 = -t907 - t975;
t765 = t853 * t1262;
t1086 = t849 * t985 - t765;
t1084 = t851 * t981 - t848;
t1075 = t1156 * t1375;
t94 = m(6) * t233 + m(7) * t179;
t111 = m(6) * t274 + m(7) * t192;
t1069 = -t975 - t1373;
t1053 = -Icges(3,5) * t981 - Icges(3,6) * t984;
t1031 = t1069 - t885;
t1030 = t1069 - t887;
t1022 = t1031 - t886;
t1017 = -t1075 + t1176;
t123 = t257 * t982 - t258 * t985;
t124 = t259 * t982 - t260 * t985;
t1015 = (t124 + (t1522 * t982 + t1518) * t985 - t1521 * t977) * t1460 + (t123 + (t1521 * t985 + t1518) * t982 - t1522 * t978) * t1458;
t1014 = t1508 * t834 - t459 * t522;
t1009 = t1494 + t1575;
t1006 = t415 + t1015;
t1005 = t124 * t1116 + t123 * t1119 + t44 * t1458 + t45 * t1460 + (-t1333 + t1334) * t1462 + (t394 * t982 - t396 * t985) * t1545 - t1238 + t1533 * t955;
t1001 = t1118 * t1531 + t1586;
t103 = t358 * t955 + (t280 * t985 + t281 * t982) * t956;
t104 = t359 * t955 + (t282 * t985 + t283 * t982) * t956;
t996 = t103 * t1460 + t104 * t1458 + t137 * t1119 + t136 * t1116 + t185 * t1546 - t44 * t1290 / 0.2e1 - t45 * t1289 / 0.2e1 + t63 * t1547 - t1464 + (t1273 + t1256 - t1331 + t1332) * t1462;
t991 = t1001 + t1009 - t1080;
t990 = t1120 * t1531 + t1009 + t1080 - t1586;
t989 = t1001 + t1080 - t1494 + t1575;
t988 = -t1333 / 0.2e1 + t1334 / 0.2e1 + (t1170 * t961 + t1172 * t960 + t1520 * t985 + t1565 * t955 - t1567 * t956 + t1591 * t982 + t324) * t1460 + (t1171 * t961 - t1173 * t960 + t1520 * t982 - t1566 * t955 + t1568 * t956 - t1584 + t323) * t1458 - t994;
t987 = t1329 / 0.2e1 - t993 + t1519 * t1461 + (-t396 + t1519) * t1460;
t939 = -rSges(3,2) * t981 + t1366;
t893 = t1053 * t985;
t892 = t1053 * t982;
t805 = t1110 * t985;
t803 = t1110 * t982;
t698 = t1030 * t985;
t696 = t1030 * t982;
t623 = t1022 * t985;
t621 = t1022 * t982;
t617 = -t1075 + t655;
t567 = -t1277 * t851 + t1177;
t566 = -t1277 * t850 - t1178;
t565 = -t1278 * t851 - t1086;
t551 = (t1021 - t975) * t985;
t549 = -t931 + (t1031 - t725) * t982;
t544 = -t1289 * t834 + t662 * t955;
t543 = t1290 * t834 - t663 * t955;
t521 = -t1376 / 0.2e1;
t519 = t1151 * t1516 + (m(6) + m(7)) * t1502 / 0.2e1;
t499 = t1043 * t956;
t492 = t1377 / 0.2e1;
t473 = t1017 + t1174;
t424 = -t566 * t985 + t567 * t982;
t423 = -(-(-t852 * t984 + t1310) * t982 - t848 * t985) * t985 + t565 * t982;
t414 = t1017 + t1077;
t402 = (t1303 + (-t1189 * t983 + t1190 * t980) * t956) * t955;
t366 = -t1385 / 0.2e1;
t355 = -t1075 + t372;
t341 = t1386 / 0.2e1;
t243 = qJD(6) * t1146;
t194 = (t565 - t765 + (t849 + t1310) * t985 + t1178) * t985 + t1177 * t982;
t193 = (t1084 * t985 - t1177 + t567) * t985 + (t1084 * t982 + t1086 + t566) * t982;
t190 = -t1387 + t1415 + t1439;
t187 = t1398 / 0.2e1;
t183 = t1399 / 0.2e1;
t168 = t1500 * t834 + t250;
t139 = t1501 * t834 + t209;
t122 = t366 + t1376 / 0.2e1;
t121 = t521 + t366;
t120 = t521 + t1385 / 0.2e1;
t100 = t341 - t1377 / 0.2e1;
t99 = t492 + t341;
t98 = t492 - t1386 / 0.2e1;
t92 = t1465 / 0.2e1;
t88 = t1466 / 0.2e1;
t82 = t1129 + t1131;
t79 = t1132 + t1134;
t62 = t1390 + t1416 + t1440 + t993 + t1448;
t39 = t1529 + (t935 / 0.2e1 + t934 / 0.2e1) * t984 + t1455 + t1449 + t1441 + t1417 + t1391 + t993;
t35 = t1473 / 0.2e1;
t33 = t1474 / 0.2e1;
t30 = m(7) * t168 + t1238;
t29 = m(7) * t139 + t1238;
t27 = t1228 + t1233;
t24 = t1234 + t1235;
t22 = t1239 - t1523;
t21 = t1240 + t1356 - t1239;
t20 = t1239 + t1523;
t18 = t1006 + t1404 + t1426 + t1444;
t17 = t1015 + t1585;
t16 = t35 - t1474 / 0.2e1 + t1144;
t15 = t33 - t1473 / 0.2e1 + t1144;
t11 = (t423 / 0.2e1 + t193 / 0.2e1) * t982 + (t424 / 0.2e1 - t194 / 0.2e1) * t985 + t994;
t10 = t33 + t35 - t1470 / 0.2e1 + t1005;
t9 = t92 + t187 + t989;
t8 = t187 + t990 - t1465 / 0.2e1;
t7 = t991 + t92 - t1398 / 0.2e1;
t6 = t183 + t88 + t989;
t5 = t183 + t990 - t1466 / 0.2e1;
t4 = t88 + t991 - t1399 / 0.2e1;
t3 = t994 + t1526;
t2 = t994 - t1526;
t1 = t988 + t1081 + t1082;
t12 = [t39 * qJD(2) + t62 * qJD(3) + t190 * qJD(4) + t229 * qJD(5) + t191 * qJD(6), t39 * qJD(1) + t1 * qJD(3) + t79 * qJD(4) + t24 * qJD(5) + t6 * qJD(6) + (m(3) * ((-t739 * t985 - t740 * t982) * t939 + (-t898 * t985 + t899 * t982) * t937) / 0.2e1 + (t690 * t805 + t691 * t803) * t1483 + (t649 * t698 + t650 * t696 - t695 * t700 - t697 * t699) * t1482 + (t588 * t623 + t589 * t621 - t618 * t620 - t619 * t622) * t1481 + (t1527 * t551 + t483 * t549 - t535 * t548 - t536 * t550) * t1479) * t1487 + (t988 + (t1165 * t984 + t1167 * t981) * t1460 + t194 * t1457 + (t423 + t193) * t1461 + (t1166 * t984 - t1168 * t981 + t424) * t1458 + (t977 / 0.2e1 + t978 / 0.2e1) * t1054) * qJD(2), t62 * qJD(1) + t1 * qJD(2) + t988 * qJD(3) + t82 * qJD(4) + t27 * qJD(5) + t9 * qJD(6) + ((t1018 + (-t866 * t985 + t867 * t982) * t906) * t1483 + (t1208 + (-t711 + t735) * t1068) * t1482 + (-t643 * t647 - t644 * t645 + t1218) * t1481 + (-t561 * t574 - t562 * t572 + t1223) * t1479) * t1485, qJD(1) * t190 + qJD(2) * t79 + qJD(3) * t82 + qJD(5) * t519 + qJD(6) * t121, qJD(2) * t24 + qJD(3) * t27 + qJD(4) * t519 + qJD(6) * t99 + t1570, t1577 + t6 * qJD(2) + t9 * qJD(3) + t121 * qJD(4) + t99 * qJD(5) + (m(7) * (t1527 * t543 - t1560 * t663 + t483 * t544 - t514 * t662) + t402 + ((t312 / 0.2e1 + t358 / 0.2e1) * t985 + (t313 / 0.2e1 + t359 / 0.2e1) * t982) * t956) * qJD(6); (t987 - t1529 - (t935 + t934) * t984 / 0.2e1) * qJD(1) + t11 * qJD(2) + t2 * qJD(3) - t78 * qJD(4) - t23 * qJD(5) + t5 * qJD(6) + (-t1455 / 0.4e1 - t1449 / 0.4e1 - t1441 / 0.4e1 - t1417 / 0.4e1 - t1391 / 0.4e1) * t1488, t11 * qJD(1) + (m(7) * (t264 * t355 - t548 * t549 - t550 * t551) + m(6) * (t331 * t414 - t620 * t621 - t622 * t623) + m(5) * (t376 * t473 - t695 * t696 - t697 * t698) + m(4) * (-t1505 * t805 - t1506 * t803 + t490 * t617) + (t977 * t893 + (-t982 * t892 + t1493) * t985) * t1460 + (t978 * t892 + (-t985 * t893 + t1493) * t982) * t1458 + m(3) * ((t982 * (rSges(3,1) * t1262 - t1157) + t985 * (rSges(3,1) * t1260 + rSges(3,3) * t982 - t951)) * (-t898 * t982 - t899 * t985) + t1156 * t939 * t937) + t1015) * qJD(2) + t17 * qJD(3) + t94 * qJD(5) + t29 * qJD(6), t2 * qJD(1) + t17 * qJD(2) + t1015 * qJD(3) + t21 * qJD(5) + t15 * qJD(6) + (-t1404 / 0.4e1 - t1426 / 0.4e1 - t1444 / 0.4e1 - t415 / 0.4e1) * t1484 + ((t1128 + t110) * t1479 + (t1127 + t165) * t1481 + (t1126 + t234) * t1482 + (t337 + t419) * t1483) * t1485, t243 - t1573, qJD(2) * t94 + qJD(3) * t21 + t1354 - t1590, t5 * qJD(1) + t29 * qJD(2) + t15 * qJD(3) + (m(7) * (-t264 * t499 - t543 * t550 - t544 * t548 + t1014) + t996) * qJD(6) + t1355; t987 * qJD(1) + t3 * qJD(2) + t994 * qJD(3) - t81 * qJD(4) + t1525 * qJD(5) + t8 * qJD(6) + (-t1448 / 0.4e1 - t1440 / 0.4e1 - t1416 / 0.4e1 - t1390 / 0.4e1) * t1488, t3 * qJD(1) + t18 * qJD(3) + t22 * qJD(5) + t16 * qJD(6) + ((t264 * t372 + t308 * t355 - t549 * t572 - t551 * t574 + t1221) * t1479 + (t331 * t430 + t388 * t414 - t621 * t645 - t623 * t647 + t1215) * t1481 + (-t1068 * t698 + t376 * t510 + t441 * t473 - t696 * t735 + t1207) * t1482 + (t617 * t628 + (-t803 * t982 - t805 * t985) * t906 + t337) * t1483) * t1487 + (t1015 - t1585) * qJD(2), t994 * qJD(1) + t18 * qJD(2) + t1006 * qJD(3) + t111 * qJD(5) + t30 * qJD(6) + (m(5) * t1126 / 0.4e1 + t1128 * t1478 + t1127 * t1480) * t1484, t243 - t1574, qJD(2) * t22 + qJD(3) * t111 + t1354 + t1582, t8 * qJD(1) + t16 * qJD(2) + t30 * qJD(3) + (m(7) * (-t308 * t499 - t543 * t574 - t544 * t572 + t1014) + t996) * qJD(6) + t1355; t78 * qJD(2) + t81 * qJD(3) + t520 * qJD(5) + t120 * qJD(6) + (t1387 / 0.4e1 - t1415 / 0.4e1 - t1439 / 0.4e1) * t1488, t1573 + ((-t549 * t985 + t551 * t982) * t1479 + (-t621 * t985 + t623 * t982) * t1481 + (-t696 * t985 + t698 * t982) * t1482) * t1487 + t1227, qJD(2) * t263 + t1227 + t1574, 0, t1155, t120 * qJD(1) + (t543 * t982 - t544 * t985) * t1367 + t1226; t23 * qJD(2) - qJD(3) * t1525 - t520 * qJD(4) + t98 * qJD(6) - t1570, t1590 + (m(7) * (-t355 * t956 + (t549 * t982 + t551 * t985) * t955 + t1136) + m(6) * (-t414 * t956 + (t621 * t982 + t623 * t985) * t955 + t1135) - t94) * qJD(2) + t20 * qJD(3) + t1353, -t1582 + t20 * qJD(2) + (m(7) * (t1062 - t1328) + m(6) * (t1061 - t1327) - t111) * qJD(3) + t1353, -t1155, t1154 * t515, t98 * qJD(1) + (t499 * t956 + (t543 * t985 + t544 * t982) * t955) * t1367 + t1154 * t52; t4 * qJD(2) + t7 * qJD(3) + t122 * qJD(4) + t100 * qJD(5) - t1577, t4 * qJD(1) + ((t1560 * t551 - t355 * t459 - t514 * t549 + t1137 - t139) * m(7) + t1005) * qJD(2) + t10 * qJD(3) + t1148, t7 * qJD(1) + t10 * qJD(2) + ((-t291 - t168 + t1063) * m(7) + t1005) * qJD(3) + t1148, qJD(1) * t122 + t1226, qJD(1) * t100 + t1154 * t54 (t402 * t1462 + m(7) * (t1560 * t543 + t459 * t499 - t514 * t544) + (t103 * t1457 + t104 * t1460 + (t312 * t985 + t313 * t982) * t1462) * t956) * qJD(6) + t1154 * t14;];
Cq  = t12;
