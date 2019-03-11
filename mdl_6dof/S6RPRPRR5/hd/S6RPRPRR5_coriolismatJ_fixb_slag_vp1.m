% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:48:02
% EndTime: 2019-03-09 03:49:20
% DurationCPUTime: 71.70s
% Computational Cost: add. (190344->1163), mult. (254485->1609), div. (0->0), fcn. (306954->9), ass. (0->729)
t1281 = cos(qJ(5));
t768 = cos(qJ(1));
t1000 = t768 * t1281;
t763 = pkin(10) + qJ(3);
t751 = sin(t763);
t752 = cos(t763);
t1279 = sin(qJ(5));
t998 = t768 * t1279;
t678 = -t1000 * t751 + t752 * t998;
t1278 = sin(qJ(6));
t1280 = cos(qJ(6));
t677 = -t1000 * t752 - t751 * t998;
t767 = sin(qJ(1));
t879 = t1278 * t677 - t1280 * t767;
t930 = t767 * t1278 + t1280 * t677;
t490 = -Icges(7,5) * t930 + Icges(7,6) * t879 + Icges(7,3) * t678;
t1008 = Icges(7,4) * t1280;
t881 = Icges(7,2) * t1278 - t1008;
t516 = -Icges(7,6) * t677 + t678 * t881;
t1007 = Icges(7,4) * t1278;
t882 = -Icges(7,1) * t1280 + t1007;
t518 = -Icges(7,5) * t677 + t678 * t882;
t1378 = Icges(7,4) * t930;
t492 = Icges(7,2) * t879 + Icges(7,6) * t678 - t1378;
t619 = Icges(7,4) * t879;
t494 = -Icges(7,1) * t930 + Icges(7,5) * t678 + t619;
t1353 = t1278 * t492 - t1280 * t494;
t880 = -Icges(7,5) * t1280 + Icges(7,6) * t1278;
t868 = -Icges(7,3) * t677 + t678 * t880 + t1353;
t210 = -t677 * t490 + t516 * t879 - t518 * t930 + t678 * t868;
t707 = -t752 * t1279 + t751 * t1281;
t676 = t707 * t767;
t1005 = t676 * t1278;
t1006 = t676 * t1280;
t1170 = Icges(7,6) * t676;
t878 = -t1279 * t751 - t1281 * t752;
t849 = t878 * t1280;
t997 = t768 * t1278;
t819 = -t767 * t849 + t997;
t848 = t878 * t1278;
t999 = t768 * t1280;
t820 = t767 * t848 + t999;
t791 = Icges(7,4) * t819 + Icges(7,2) * t820 - t1170;
t1176 = Icges(7,5) * t676;
t792 = Icges(7,1) * t819 + Icges(7,4) * t820 - t1176;
t1360 = t791 * t1278 - t792 * t1280;
t1377 = t767 * t878;
t1526 = -Icges(7,5) * t1006 + Icges(7,6) * t1005 + Icges(7,3) * t1377 - t1360;
t1166 = Icges(7,3) * t676;
t790 = Icges(7,5) * t819 + Icges(7,6) * t820 - t1166;
t803 = -Icges(7,4) * t1006 + Icges(7,2) * t1005 + Icges(7,6) * t1377;
t805 = -Icges(7,1) * t1006 + Icges(7,4) * t1005 + Icges(7,5) * t1377;
t773 = t1526 * t678 + t677 * t790 + t803 * t879 - t805 * t930;
t101 = t210 * t767 + t768 * t773;
t643 = Icges(6,4) * t678;
t563 = -Icges(6,1) * t677 - Icges(6,5) * t767 - t643;
t1065 = Icges(6,2) * t677 + t563 - t643;
t644 = Icges(6,4) * t677;
t561 = -Icges(6,2) * t678 - Icges(6,6) * t767 - t644;
t1414 = Icges(6,1) * t678 + t561 - t644;
t641 = Icges(6,4) * t676;
t1397 = -Icges(6,2) * t1377 - t641;
t922 = Icges(6,5) * t768 + t641;
t562 = Icges(6,1) * t1377 - t922;
t1494 = t1397 + t562;
t1413 = Icges(6,4) * t1377;
t919 = Icges(6,2) * t676 + Icges(6,6) * t768;
t560 = t1413 - t919;
t1545 = Icges(6,1) * t676 + t1413 + t560;
t568 = Icges(6,5) * t678 - Icges(6,6) * t677;
t826 = Icges(6,5) * t676 + Icges(6,6) * t1377;
t192 = (t1377 * t1545 + t1494 * t676 - t768 * t826) * t768 + (t1065 * t676 + t1377 * t1414 - t768 * t568) * t767;
t194 = (-t1494 * t678 + t1545 * t677 + t767 * t826) * t768 - (t1065 * t678 - t1414 * t677 - t568 * t767) * t767;
t207 = t1377 * t490 - t516 * t820 - t518 * t819 + t676 * t868;
t1420 = t790 * t878;
t771 = -t1420 * t767 + t1526 * t676 - t803 * t820 - t805 * t819;
t98 = t207 * t767 + t768 * t771;
t841 = (t101 / 0.2e1 + t194 / 0.2e1) * t767 + (t98 / 0.2e1 - t192 / 0.2e1) * t768;
t1284 = t768 / 0.2e1;
t1290 = -t767 / 0.2e1;
t1340 = (t101 + t194) * t1290 + (t192 - t98) * t1284;
t330 = -t1360 * t707 - t1420;
t332 = -t1353 * t707 - t490 * t878;
t1109 = t707 * (-t330 * t768 + t332 * t767);
t780 = t678 * t790 + t791 * t879 - t792 * t930;
t1470 = t768 * t780;
t315 = t678 * t490 + t492 * t879 - t494 * t930;
t180 = t315 * t767 - t1470;
t1125 = t677 * t180;
t545 = -Icges(7,3) * t878 - t707 * t880;
t546 = -Icges(7,6) * t707 + t878 * t881;
t549 = -Icges(7,5) * t707 + t878 * t882;
t548 = -Icges(7,6) * t878 - t707 * t881;
t987 = t1278 * t548;
t551 = -Icges(7,5) * t878 - t707 * t882;
t992 = t1280 * t551;
t1354 = t987 - t992;
t543 = -Icges(7,3) * t707 + t878 * t880;
t867 = t543 - t1354;
t251 = -t1377 * t545 - t546 * t820 - t549 * t819 + t676 * t867;
t1128 = t676 * t490;
t314 = t492 * t820 + t494 * t819 - t1128;
t1127 = t676 * t545;
t356 = t548 * t820 + t551 * t819 - t1127;
t789 = t676 * t790;
t777 = t791 * t820 + t792 * t819 - t789;
t776 = t777 * t767;
t28 = t676 * t771 + t207 * t678 + t314 * t677 - t356 * t707 + (t251 + t776) * t878;
t1560 = t768 * t28;
t1561 = t676 * t98;
t1003 = t707 * t1278;
t1004 = t707 * t1280;
t229 = t1003 * t803 - t1004 * t805 + t1526 * t878 + t707 * t790;
t1532 = t229 * t768;
t228 = t868 * t878 + (t1278 * t516 - t1280 * t518 - t490) * t707;
t1546 = t228 * t767;
t1571 = t878 * (-t1546 - t1532);
t1574 = t678 * t101;
t252 = t677 * t545 + t546 * t879 - t549 * t930 + t678 * t867;
t358 = t678 * t545 + t548 * t879 - t551 * t930;
t779 = t767 * t780;
t31 = t210 * t678 - t315 * t677 + t358 * t707 + t676 * t773 + (t252 - t779) * t878;
t1575 = t767 * t31;
t1288 = t767 / 0.2e1;
t1367 = Icges(7,2) * t930 + t494 + t619;
t1368 = -Icges(7,1) * t879 - t1378 + t492;
t529 = Icges(7,5) * t879 + Icges(7,6) * t930;
t243 = t1367 * t820 - t1368 * t819 - t676 * t529;
t1358 = t792 - t881 * t768 + (Icges(7,4) * t848 + Icges(7,2) * t849) * t767;
t1359 = t791 + t882 * t768 - (Icges(7,1) * t848 + Icges(7,4) * t849) * t767;
t796 = -t880 * t768 + (Icges(7,5) * t848 + Icges(7,6) * t849) * t767;
t769 = t1358 * t820 - t1359 * t819 - t676 * t796;
t131 = t243 * t767 - t768 * t769;
t244 = t1367 * t879 + t1368 * t930 + t678 * t529;
t772 = t1358 * t879 + t1359 * t930 + t678 * t796;
t132 = t244 * t767 - t768 * t772;
t1380 = -t768 / 0.2e1;
t891 = t1288 * t132 + t131 * t1380;
t1577 = t1125 / 0.2e1 - t1109 / 0.2e1 - t891 - t1560 / 0.2e1 - t1561 / 0.2e1 + t1571 / 0.2e1 - t1574 / 0.2e1 - t1575 / 0.2e1;
t1576 = t1560 / 0.4e1 + t1561 / 0.4e1 - t1571 / 0.4e1 + t1574 / 0.4e1 + t1575 / 0.4e1;
t658 = t676 * rSges(7,3);
t495 = rSges(7,1) * t819 + rSges(7,2) * t820 - t658;
t519 = -rSges(7,1) * t1006 + rSges(7,2) * t1005 + rSges(7,3) * t1377;
t884 = -t1280 * rSges(7,1) + t1278 * rSges(7,2);
t552 = -t707 * rSges(7,3) + t878 * t884;
t554 = -rSges(7,3) * t878 - t707 * t884;
t304 = t1377 * t554 + t707 * t495 + t519 * t878 - t676 * t552;
t521 = t677 * rSges(7,3) - t678 * t884;
t809 = -rSges(7,1) * t930 + rSges(7,2) * t879 + t678 * rSges(7,3);
t306 = t521 * t878 + t678 * t552 + t677 * t554 + t707 * t809;
t1541 = 0.2e1 * t304 * t767 + 0.2e1 * t306 * t768;
t1295 = -t1541 / 0.4e1;
t378 = -t1354 * t707 - t545 * t878;
t1149 = t378 * t707;
t1150 = t332 * t677;
t1292 = -t878 / 0.2e1;
t1293 = t678 / 0.2e1;
t1294 = -t676 / 0.2e1;
t1307 = m(7) / 0.4e1;
t1465 = t1280 * t549;
t1466 = t1278 * t546;
t1528 = ((t1466 - t1465 + t545) * t707 + t867 * t878) * t878;
t1533 = t229 * t676;
t1547 = t228 * t678;
t1563 = t707 / 0.2e1;
t1454 = -t358 * t878 - t676 * t780;
t138 = t315 * t678 + t1454;
t1565 = -t138 / 0.2e1;
t1419 = t809 * t878;
t276 = -t1419 * t767 + t677 * t495 + t678 * t519 + t676 * t521;
t382 = t678 * t495 + t676 * t809;
t1138 = t554 * t676;
t410 = t495 * t878 - t1138;
t412 = t554 * t678 + t1419;
t1566 = -t1292 * (t1377 * t330 - t1149 + t1150 + t1528 + t1533 + t1547) + t31 * t1293 - t28 * t1294 - 0.4e1 * (t276 * t382 + t304 * t410 + t306 * t412) * t1307 + t677 * t1565 + (-t330 * t676 + t332 * t678 - t378 * t878) * t1563;
t1564 = t678 / 0.4e1;
t1562 = -t1546 / 0.4e1;
t1399 = -pkin(5) * t1377 - t676 * pkin(9);
t1075 = t1399 + t495;
t807 = -t677 * pkin(5) + t678 * pkin(9) + t809;
t341 = -t1075 * t767 - t768 * t807;
t1069 = pkin(5) * t707 - pkin(9) * t878 + t554;
t477 = t1069 * t767;
t480 = t1069 * t768;
t1074 = -t676 * pkin(5) + pkin(9) * t1377 + t519;
t1409 = t678 * pkin(5) + t677 * pkin(9) + t521;
t1492 = t1074 * t767 + t1409 * t768;
t1071 = -pkin(5) * t878 - pkin(9) * t707 + t552;
t1497 = t1071 * t768;
t1498 = t1071 * t767;
t873 = -t1492 * t382 + t1497 * t410 - t1498 * t412;
t1558 = -t276 * t341 - t304 * t480 + t306 * t477 - t873;
t1209 = m(7) * qJD(6);
t1557 = t1209 * t1295;
t1314 = -0.2e1 * t752;
t1556 = t1314 * t1492;
t1212 = m(7) * qJD(2);
t1555 = t1541 * t1212 / 0.4e1;
t1434 = t676 / 0.4e1;
t1552 = -t251 * t1434 - t1150 / 0.4e1 + t1149 / 0.2e1 - t677 * t358 / 0.4e1 - t1528 / 0.2e1 - t1533 / 0.4e1 - t252 * t1564 - t1547 / 0.4e1;
t1549 = t1492 * t341 - t1497 * t480 - t1498 * t477;
t764 = t767 ^ 2;
t765 = t768 ^ 2;
t1048 = t764 + t765;
t1107 = t751 * qJ(4);
t724 = t752 * pkin(3) + t1107;
t1052 = t1048 * t724;
t1103 = t752 * t768;
t1104 = t752 * t767;
t720 = pkin(4) * t1104 + t768 * pkin(8);
t762 = t767 * pkin(8);
t939 = t767 * t720 + t1052 + t768 * (pkin(4) * t1103 - t762);
t307 = -t341 + t939;
t1222 = pkin(4) * t751;
t721 = pkin(3) * t751 - qJ(4) * t752;
t969 = t721 + t1222;
t898 = t969 + t1069;
t450 = t898 * t767;
t452 = t898 * t768;
t1548 = -t1492 * t307 - t452 * t1497 - t450 * t1498;
t701 = Icges(6,4) * t707;
t600 = Icges(6,2) * t878 + t701;
t601 = -Icges(6,1) * t878 + t701;
t1363 = t600 + t601;
t700 = Icges(6,4) * t878;
t599 = -Icges(6,2) * t707 + t700;
t603 = Icges(6,1) * t707 + t700;
t1500 = t599 + t603;
t595 = -Icges(6,5) * t878 + Icges(6,6) * t707;
t351 = t1363 * t677 - t1500 * t678 + t595 * t767;
t364 = t1065 * t878 - t1414 * t707;
t1544 = t252 - t351 - t364;
t349 = t1363 * t1377 + t1500 * t676 - t768 * t595;
t1418 = t562 * t878;
t363 = -t1397 * t878 + t1545 * t707 - t1418;
t1543 = t363 + t251 + t349;
t1151 = t306 * t767;
t1152 = t304 * t768;
t1542 = t1151 - t1152;
t1540 = t251 / 0.2e1 + t349 / 0.2e1 + t363 / 0.2e1;
t1539 = t363 / 0.4e1 + t251 / 0.4e1 + t349 / 0.4e1;
t1536 = -t364 / 0.2e1 - t351 / 0.2e1 + t228 / 0.2e1 + t252 / 0.2e1;
t1375 = t768 * (t767 * (Icges(6,5) * t1377 - Icges(6,6) * t676 - Icges(6,3) * t768) + t678 * t560 + t677 * t562);
t399 = -t678 * t561 - t677 * t563 + t767 * (Icges(6,5) * t677 + Icges(6,6) * t678 + Icges(6,3) * t767);
t116 = -t1375 + (-t1377 * t562 + t1418 * t767 + t399) * t767;
t1283 = t768 / 0.4e1;
t1285 = -t768 / 0.4e1;
t1289 = -t767 / 0.4e1;
t1089 = -t307 + t341;
t622 = -t1278 * t1377 - t999;
t623 = t1280 * t1377 - t997;
t496 = t623 * rSges(7,1) + t622 * rSges(7,2) + t658;
t1366 = -t1399 + t496;
t1374 = (t1075 + t1366) * t768;
t1040 = 0.2e1 * m(7);
t1429 = -t1040 / 0.4e1;
t1346 = t1089 * t1429 * t1374;
t179 = t314 * t767 - t768 * t777;
t1370 = -t622 * t791 - t623 * t792 - t789;
t489 = Icges(7,5) * t623 + Icges(7,6) * t622 + t1166;
t491 = Icges(7,4) * t623 + Icges(7,2) * t622 + t1170;
t493 = Icges(7,1) * t623 + Icges(7,4) * t622 + t1176;
t183 = t678 * t489 + t491 * t879 + t622 * t492 - t493 * t930 + t623 * t494 + t1128;
t66 = t1370 * t768 + t183 * t767 + t779;
t1472 = t179 + t66;
t283 = t399 * t767 - t1375;
t184 = -t676 * t489 + t491 * t820 + t493 * t819 + t315;
t67 = t184 * t767 - t1470 + t776;
t1534 = -t1346 + 0.2e1 * (t116 + t67) * t1283 + 0.2e1 * (t180 + t283) * t1285 + 0.2e1 * t1472 * t1289 + (t252 / 0.4e1 + t228 / 0.4e1 - t364 / 0.4e1 - t351 / 0.4e1) * t767;
t1527 = t276 * t1314;
t1444 = t599 / 0.2e1 + t543 / 0.2e1 + t603 / 0.2e1 + t992 / 0.2e1 - t987 / 0.2e1;
t824 = t1465 / 0.2e1 - t1466 / 0.2e1 - t545 / 0.2e1 + t600 / 0.2e1 + t601 / 0.2e1;
t1520 = t1444 * t878 - t707 * t824;
t1106 = t751 * t767;
t748 = cos(pkin(10)) * pkin(2) + pkin(1);
t1217 = -pkin(7) - qJ(2);
t971 = t768 * t1217;
t808 = -pkin(3) * t1104 - qJ(4) * t1106 - t767 * t748 - t720 - t971;
t424 = t808 - t1075;
t1306 = pkin(3) + pkin(4);
t1338 = t1306 * t752 + t1107 + t748;
t749 = t767 * t1217;
t837 = t1338 * t768 - t749 - t762;
t425 = t837 + t807;
t1518 = -t304 * t424 + t306 * t425;
t1046 = qJD(1) * t707;
t1047 = qJD(1) * t878;
t1517 = t824 * t1046 - t1444 * t1047;
t913 = -t1497 * t424 - t1498 * t425;
t909 = -t1497 * t768 - t1498 * t767;
t910 = t1497 * t767 - t1498 * t768;
t1512 = -t276 * t307 + t304 * t452 - t306 * t450;
t1310 = m(6) / 0.2e1;
t1185 = Icges(5,1) * t752;
t746 = Icges(5,5) * t751;
t924 = t746 + t1185;
t654 = Icges(5,4) * t767 + t768 * t924;
t1180 = Icges(4,4) * t751;
t719 = Icges(4,1) * t752 - t1180;
t656 = Icges(4,5) * t767 + t719 * t768;
t1499 = t654 + t656;
t575 = t676 * rSges(6,1) + rSges(6,2) * t1377;
t576 = t678 * rSges(6,1) - t677 * rSges(6,2);
t1474 = t767 * t575 - t576 * t768;
t1054 = rSges(6,1) * t1377 - t676 * rSges(6,2);
t564 = t768 * rSges(6,3) - t1054;
t566 = -t677 * rSges(6,1) - t678 * rSges(6,2) - t767 * rSges(6,3);
t468 = -t767 * t564 - t566 * t768;
t416 = -t468 + t939;
t1495 = t1474 * t416;
t1493 = t1474 * t1314;
t411 = -t496 * t878 - t1138;
t1083 = t410 - t411;
t381 = (-t495 - t496) * t676;
t1324 = -0.2e1 * t381;
t331 = -t878 * t489 + (-t1278 * t491 + t1280 * t493) * t707;
t1350 = (t330 + t331) * t878;
t137 = t314 * t678 - t356 * t878 - t676 * t777;
t1431 = t767 / 0.4e1;
t357 = t548 * t622 + t551 * t623 + t1127;
t49 = -t357 * t878 + (t183 + t780) * t678 + t1370 * t676;
t50 = (t184 + t777) * t678 + t1454;
t1372 = -t138 * t1283 - t50 * t1285 + (t1350 - t137 - t49) * t1431;
t1405 = 0.2e1 * t382 * t1374;
t1479 = (0.2e1 * t1083 * t450 + t1324 * t307 + t1405) * t1307 - t1372;
t1478 = t1372 + (-0.2e1 * t1083 * t477 + t1324 * t341 - t1405) * t1307;
t1428 = t1040 / 0.4e1;
t1467 = (-Icges(4,6) + Icges(5,6)) * t752 + (-Icges(5,4) - Icges(4,5)) * t751;
t1388 = 0.2e1 * t1048;
t606 = rSges(6,1) * t707 + rSges(6,2) * t878;
t1462 = t1388 * t606;
t1105 = t751 * t768;
t736 = Icges(5,5) * t1103;
t646 = Icges(5,6) * t767 + Icges(5,3) * t1105 + t736;
t712 = Icges(4,5) * t752 - Icges(4,6) * t751;
t648 = Icges(4,3) * t767 + t712 * t768;
t713 = Icges(5,4) * t752 + Icges(5,6) * t751;
t650 = Icges(5,2) * t767 + t713 * t768;
t1460 = t646 * t1105 + (t648 + t650) * t767 + t1499 * t1103;
t1459 = -(t1074 * t410 + t1409 * t412 + t1518) * t1429 + t1552;
t1369 = t356 + t330;
t1426 = t1377 / 0.4e1;
t1458 = t1369 * t1426;
t1427 = -t1377 / 0.4e1;
t1457 = t1369 * t1427;
t714 = Icges(4,2) * t752 + t1180;
t1167 = Icges(5,3) * t752;
t920 = t1167 - t746;
t1456 = (-t714 - t920) * t768 + t1499;
t536 = -t884 * t768 + (rSges(7,1) * t848 + rSges(7,2) * t849) * t767;
t1317 = -0.2e1 * t536;
t1318 = 0.2e1 * t477;
t591 = (-rSges(7,1) * t1278 - rSges(7,2) * t1280) * t707;
t1038 = 0.2e1 * t591;
t912 = t424 * t768 + t425 * t767;
t1352 = t1038 * t912;
t537 = rSges(7,1) * t879 + rSges(7,2) * t930;
t273 = -t878 * t529 + (-t1278 * t1367 - t1280 * t1368) * t707;
t1153 = t273 * t767;
t272 = -t1003 * t1358 - t1004 * t1359 - t796 * t878;
t1154 = t272 * t768;
t589 = (-Icges(7,2) * t1280 - t1007) * t707;
t1364 = t551 + t589;
t590 = (-Icges(7,1) * t1278 - t1008) * t707;
t1365 = t548 - t590;
t588 = (-Icges(7,5) * t1278 - Icges(7,6) * t1280) * t707;
t293 = t1364 * t820 - t1365 * t819 - t676 * t588;
t294 = t1364 * t879 + t1365 * t930 + t678 * t588;
t944 = t1154 / 0.4e1 - t1153 / 0.4e1 + t294 * t1289 + t293 * t1283;
t928 = (t1317 * t480 + t1318 * t537 + t1352) * t1307 + t944;
t1319 = -0.2e1 * t450;
t929 = (-t1317 * t452 + t1319 * t537 - t1352) * t1307 - t944;
t1355 = t1109 / 0.4e1 - t1125 / 0.4e1;
t1079 = -t452 - t480;
t1080 = -t450 - t477;
t442 = t767 * t536 + t537 * t768;
t1400 = (t1089 * t442 + (t1079 * t768 + t1080 * t767) * t591) * t1428 - t891;
t1448 = t1426 * t179 - t1355 + t1400;
t1411 = t751 * t1388;
t1447 = t1314 * t442 - t591 * t1411;
t1445 = t1427 * t179 + t1355;
t1313 = -m(4) / 0.4e1;
t1437 = 0.4e1 * t1313;
t1436 = m(5) / 0.2e1;
t1311 = m(5) / 0.4e1;
t1309 = m(6) / 0.4e1;
t1435 = -t67 / 0.4e1;
t1430 = t878 / 0.2e1;
t604 = -rSges(6,1) * t878 + rSges(6,2) * t707;
t485 = t808 - t564;
t486 = t566 + t837;
t908 = t485 * t768 + t486 * t767;
t1423 = t604 * t908;
t935 = 0.4e1 * t1048;
t1422 = t604 * t606 * t935;
t933 = t606 + t969;
t534 = t933 * t768;
t1139 = t534 * t768;
t532 = t933 * t767;
t1140 = t532 * t767;
t907 = t1139 + t1140;
t1417 = t907 * t604;
t647 = Icges(4,5) * t1104 - Icges(4,6) * t1106 - Icges(4,3) * t768;
t737 = Icges(4,4) * t1106;
t655 = Icges(4,1) * t1104 - Icges(4,5) * t768 - t737;
t1062 = -t655 * t1103 - t767 * t647;
t651 = Icges(4,4) * t1104 - Icges(4,2) * t1106 - Icges(4,6) * t768;
t1173 = Icges(4,2) * t751;
t747 = Icges(4,4) * t752;
t652 = Icges(4,6) * t767 + (t747 - t1173) * t768;
t626 = t656 * t1104;
t952 = t648 * t768 - t626;
t1416 = -t1105 * t651 - t1106 * t652 - t1062 - t952;
t1415 = -t1105 * t652 + t1460;
t1412 = t1388 * t604;
t1134 = (-Icges(5,2) * t768 + t767 * t713) * t768;
t1410 = t1134 + t1460;
t741 = pkin(3) * t1106;
t932 = qJ(4) * t1104 - t741;
t872 = pkin(4) * t1106 - t932;
t437 = t872 - t1074;
t733 = qJ(4) * t1103;
t897 = -t1105 * t1306 + t733;
t438 = t897 + t1409;
t1408 = -(t410 * t437 - t412 * t438 - t1518) * t1429 - t1552;
t1207 = rSges(5,1) * t751;
t722 = -rSges(5,3) * t752 + t1207;
t1051 = t721 + t722;
t614 = t1051 * t767;
t1137 = t614 * t767;
t1144 = t450 * t767;
t449 = -0.2e1 * t1144;
t528 = -0.2e1 * t1140;
t592 = -0.2e1 * t1137;
t1009 = (t449 + 0.2e1 * t1144) * t1307 + (t528 + 0.2e1 * t1140) * t1309 + (t592 + 0.2e1 * t1137) * t1311;
t617 = t1051 * t768;
t1136 = t617 * t768;
t1143 = t452 * t768;
t723 = rSges(4,1) * t751 + rSges(4,2) * t752;
t941 = (t449 - 0.2e1 * t1143) * t1307 + (t528 - 0.2e1 * t1139) * t1309 + (t592 - 0.2e1 * t1136) * t1311 + t723 * t1388 * t1313;
t1407 = t941 - t1009;
t1406 = t1314 * t1374;
t1032 = -0.2e1 * t1106;
t1030 = 0.2e1 * t1103;
t1031 = 0.2e1 * t1104;
t1187 = rSges(5,3) + qJ(4);
t1282 = rSges(5,1) + pkin(3);
t1333 = t1187 * t751 + t1282 * t752 + t748;
t761 = t768 * rSges(5,2);
t555 = -t1333 * t767 + t761 - t971;
t1198 = t767 * rSges(5,2);
t556 = t1333 * t768 + t1198 - t749;
t1072 = t555 * t1030 + t556 * t1031;
t1077 = t485 * t1030 + t486 * t1031;
t1082 = t424 * t1030 + t425 * t1031;
t1023 = t450 * t1105;
t448 = -0.2e1 * t1023;
t1022 = t532 * t1105;
t525 = -0.2e1 * t1022;
t1021 = t614 * t1105;
t587 = -0.2e1 * t1021;
t1011 = (-t1032 * t452 + t1082 + t448) * t1307 + (-t1032 * t534 + t1077 + t525) * t1309 + (-t1032 * t617 + t1072 + t587) * t1311;
t1036 = 0.2e1 * t751;
t611 = t741 + (-t1187 * t752 + t1207) * t767;
t740 = rSges(5,3) * t1103;
t612 = -t1105 * t1282 + t733 + t740;
t526 = t575 + t872;
t527 = t897 + t576;
t883 = 0.2e1 * t526 * t768 + 0.2e1 * t527 * t767;
t1012 = ((t437 * t768 + t438 * t767) * t1036 + t1082) * t1307 + (t751 * t883 + t1077) * t1309 + ((t611 * t768 + t612 * t767) * t1036 + t1072) * t1311;
t1404 = t1011 - t1012;
t1141 = t480 * t768;
t1142 = t477 * t767;
t470 = 0.2e1 * t1142;
t1085 = (t470 + 0.2e1 * t1141) * t1307 + t1309 * t1462;
t1254 = m(6) * t1474;
t1086 = t1492 * t1428 - t1254 / 0.2e1;
t1403 = t1085 - t1086;
t1379 = m(7) * t1036;
t1087 = (t1074 * t768 - t1409 * t767) * t1379 / 0.4e1 + (-t575 * t768 - t576 * t767) * t1036 * t1309;
t1024 = t477 * t1105;
t465 = 0.2e1 * t1024;
t1091 = (-0.2e1 * t1024 + t465 - t1406) * t1307;
t1402 = t1087 - t1091;
t1063 = t752 * t1462;
t1078 = t480 * t1030 + t477 * t1031;
t709 = t1048 * t751;
t1316 = 0.2e1 * t709;
t1092 = (t1316 * t341 + t1078) * t1307 + (t1316 * t468 + t1063) * t1309;
t1094 = (-t1556 + (t341 - t909) * t1036 + t1078) * t1307 + (t1493 + (0.2e1 * t468 + t1412) * t751 + t1063) * t1309;
t1401 = t1092 - t1094;
t1389 = -0.4e1 * t341;
t1383 = t179 / 0.2e1;
t844 = t1377 / 0.2e1;
t1179 = Icges(5,5) * t752;
t711 = Icges(5,3) * t751 + t1179;
t645 = -Icges(5,6) * t768 + t711 * t767;
t653 = -Icges(5,4) * t768 + t767 * t924;
t1376 = (t645 * t751 + t653 * t752) * t767;
t1034 = 0.4e1 * t752;
t948 = t1311 + t1309 + t1307;
t459 = t948 * (-0.1e1 + t1048) * t751 * t1034;
t1362 = t1467 * t768;
t1361 = t1467 * t767;
t1356 = t1456 * t767;
t1090 = (-t1036 * t909 - t1556) * t1307 + (t1411 * t604 + t1493) * t1309;
t1347 = (-t1074 * t452 + t1409 * t450 - t913) * t1429 - m(6) * (t532 * t576 + t534 * t575 + t1423) / 0.2e1;
t1146 = t412 * t767;
t1147 = t410 * t768;
t1337 = t1038 * (-t1146 + t1147) - 0.2e1 * t382 * t442;
t1336 = (t437 * t480 + t438 * t477 - t913) * t1429 - m(6) * (t606 * t883 + 0.2e1 * t1423) / 0.4e1 + t1532 / 0.4e1;
t1010 = (0.2e1 * t1023 + t448 + t1406) * t1307 + (0.2e1 * t1022 + t525) * t1309 + (0.2e1 * t1021 + t587) * t1311;
t1335 = -(t180 / 0.4e1 + t1435) * t676 - (t179 / 0.4e1 + t66 / 0.4e1) * t678;
t1186 = Icges(4,1) * t751;
t1049 = rSges(4,2) * t1106 + t768 * rSges(4,3);
t1208 = rSges(4,1) * t752;
t962 = t748 + t1208;
t593 = -t767 * t962 + t1049 - t971;
t961 = -rSges(4,2) * t1105 + t767 * rSges(4,3);
t594 = t768 * t962 - t749 + t961;
t697 = t723 * t767;
t699 = t723 * t768;
t716 = Icges(5,1) * t751 - t1179;
t1334 = -t751 * (t719 / 0.2e1 - t714 / 0.2e1 + t746 + t1185 / 0.2e1 - t1167 / 0.2e1) - t752 * (t747 + t1186 / 0.2e1 - t1173 / 0.2e1 + t716 / 0.2e1 - t711 / 0.2e1) - (t555 * t611 + t556 * t612) * m(5) + (t593 * t697 - t594 * t699) * t1437;
t1332 = t180 * t1434 + t676 * t1435 + t1472 * t1564;
t900 = t767 * t697 + t699 * t768;
t942 = (t767 * t437 - t438 * t768) * t1428 + (t767 * t526 - t527 * t768) * t1310 + (t767 * t611 - t612 * t768) * t1436 + m(4) * t900 / 0.2e1;
t693 = t716 * t767;
t694 = -Icges(5,1) * t1105 + t736;
t925 = -t747 - t1186;
t695 = t925 * t767;
t696 = t925 * t768;
t1331 = ((t651 - t695 - t645 + t693) * t768 + (-t652 + t696 + t646 + t694) * t767) * t752;
t85 = t243 * t678 - t293 * t878 - t676 * t769;
t86 = t244 * t678 - t294 * t878 - t676 * t772;
t1330 = (t1153 - t1154) * t1430 - t678 * t132 / 0.2e1 + t676 * t131 / 0.2e1 + t85 * t1284 + t86 * t1290;
t1329 = (t1512 + t873) * t1428 + t1576;
t1328 = -t1558 * t1428 + t1576;
t1325 = 0.4e1 * t307;
t427 = t536 * t678 + t537 * t676;
t1321 = 0.2e1 * t427;
t435 = t536 * t878 - t591 * t676;
t1320 = 0.2e1 * t435;
t1315 = 0.4e1 * t709;
t1308 = m(7) / 0.2e1;
t409 = -0.2e1 * t1146;
t1084 = t410 * t1030 + t752 * t409;
t1300 = m(7) * (t1527 + (t382 - t1542) * t1036 + t1084);
t1296 = t137 / 0.2e1;
t1035 = 0.4e1 * t751;
t1272 = m(5) * (-t555 * t767 + t556 * t768) * t1035;
t295 = t410 * t1032 - 0.2e1 * t412 * t1105;
t1249 = m(7) * (-t381 * t1314 + (t411 * t767 + t412 * t768) * t1036 + t295);
t1037 = 0.4e1 * t591;
t911 = t1143 + t1144;
t1247 = m(7) * (t1037 * t911 + t1325 * t442);
t297 = 0.2e1 * t1147 + t409;
t1244 = m(7) * (-0.2e1 * t411 * t768 + 0.2e1 * t1146 + t297);
t1243 = m(7) * (t442 * t1389 + (t1141 + t1142) * t1037);
t1239 = m(7) * (t1316 * t382 + t1084);
t1235 = m(7) * t295;
t1234 = m(7) * t297;
t1228 = m(7) * t1447;
t1227 = m(7) * (t1032 * t480 + t465);
t1226 = m(7) * (t470 - 0.2e1 * t1142);
t1224 = (-t536 * t768 + t537 * t767) * t1379;
t1223 = t442 * t1040;
t1216 = m(5) * qJD(3);
t1215 = m(6) * qJD(1);
t1214 = m(6) * qJD(3);
t1213 = m(7) * qJD(1);
t1211 = m(7) * qJD(3);
t1210 = m(7) * qJD(5);
t1200 = t678 * t1350;
t1133 = t651 * t751;
t1110 = t878 * t588;
t1027 = t1541 * t1307;
t1093 = qJD(3) * t1027 + t1210 * t1295;
t1039 = -pkin(8) - t1217;
t864 = t1338 * t767;
t1081 = -t1039 * t768 - t1366 + t424 + t864;
t1076 = (-rSges(6,3) + t1039) * t768 - t864 + t1054 - t485;
t1058 = -t920 * t767 + t653;
t1056 = -Icges(4,2) * t1104 + t655 - t737;
t1053 = t767 * t932 + t768 * (-pkin(3) * t1105 + t733);
t725 = t752 * rSges(5,1) + t751 * rSges(5,3);
t1050 = -t724 - t725;
t1045 = qJD(1) * t767;
t947 = t1308 + t1310 + t1436;
t463 = t947 * t1411;
t1042 = t463 * qJD(1);
t1016 = t1565 + t50 / 0.2e1;
t1013 = t49 / 0.2e1 + t1296;
t995 = t1280 * t548;
t991 = t1280 * t590;
t986 = t1278 * t551;
t985 = t1278 * t589;
t973 = t713 / 0.2e1 + t712 / 0.2e1;
t968 = -t752 * pkin(4) - t724;
t967 = 0.2e1 * t1310;
t965 = t1040 / 0.2e1;
t958 = t1058 * t768;
t956 = t1056 * t768;
t951 = t652 * t751 - t647;
t934 = -t604 + t968;
t931 = -t654 * t1104 - t646 * t1106 + t650 * t768;
t533 = t934 * t767;
t535 = t934 * t768;
t906 = t533 * t767 + t535 * t768;
t905 = t593 * t768 + t594 * t767;
t903 = -t1136 - t1137;
t899 = t968 - t1071;
t894 = t1383 + t66 / 0.2e1;
t893 = t67 / 0.2e1 - t180 / 0.2e1 - t283 / 0.2e1 + t116 / 0.2e1;
t890 = t330 / 0.2e1 + t357 / 0.2e1 + t331 / 0.2e1 + t356 / 0.2e1;
t88 = ((t442 / 0.2e1 - t276 / 0.2e1) * t752 + (t1152 / 0.2e1 - t1151 / 0.2e1 + (t764 / 0.2e1 + t765 / 0.2e1) * t591) * t751) * m(7);
t887 = -t88 * qJD(4) - t1555;
t885 = -t1048 * t1222 + t1053;
t871 = t1328 + t1445;
t870 = t1329 + t1445;
t865 = -(t555 * t768 + t556 * t767) * m(5) + t905 * t1437 - m(3) * t1048 * (rSges(3,3) + qJ(2));
t860 = t1332 + t1479;
t859 = -t1332 + t1478;
t851 = -t1543 * t1283 + t1544 * t1289 - t1336 + t1562;
t850 = -t1532 / 0.4e1 + t1562 - t1347 - t1544 * t1431 + t1543 * t1285;
t840 = t1408 + t1458;
t839 = t1457 + t1459;
t497 = -t1134 + t1376;
t833 = t931 * t1290 + t894 + (-t767 * (-t655 * t752 + t1133) - t647 * t768 + t497) * t1380 + (t768 * t951 - t1410 + t1415) * t1284 + (t653 * t1103 + t645 * t1105 + t767 * t951 + t1416 + t931 + t952) * t1288;
t832 = -t893 + (t497 - t1376 + t1410) * t1290 + t1415 * t1288 + (-t626 + (t648 + t1133) * t768 + t1062 + t1416) * t1380;
t829 = t890 + (t919 + t560) * t1430 + (t922 + t562) * t1563 + (-t1377 + t844) * t603;
t815 = t767 * t894 - t768 * t893;
t814 = t986 / 0.2e1 - t991 / 0.2e1 + t995 / 0.2e1 + t985 / 0.2e1;
t15 = -t1377 * t1296 + t1566;
t14 = t137 * t844 - t1566;
t726 = -rSges(4,2) * t751 + t1208;
t618 = t1050 * t768;
t615 = t1050 * t767;
t505 = t768 * (-rSges(5,1) * t1105 + t740) - t722 * t764 + t1053;
t487 = t767 * (t725 * t767 - t761) + (t725 * t768 + t1198) * t768 + t1052;
t462 = -t948 * t1411 + (m(5) + m(6) + m(7)) * t1411 / 0.4e1;
t453 = t899 * t768;
t451 = t899 * t767;
t439 = -t1223 / 0.4e1;
t436 = -t537 * t878 - t591 * t678;
t433 = t1224 / 0.4e1;
t428 = -t1474 + t885;
t406 = t1034 * t903 + t1315 * t487;
t400 = 0.4e1 * t908;
t395 = (-t485 * t767 + t486 * t768) * t1035;
t385 = t1226 / 0.4e1;
t379 = t1227 / 0.4e1;
t376 = t1228 / 0.4e1;
t340 = -0.4e1 * t485 * t575 - 0.4e1 * t486 * t576;
t333 = t1492 + t885;
t322 = 0.4e1 * t485 * t526 + 0.4e1 * t486 * t527;
t316 = -t1034 * t907 + t1315 * t416;
t309 = (-t1110 + (-t986 - t985 - t995 + t991) * t707) * t878;
t302 = 0.4e1 * t1474 * t468 + t1422;
t301 = 0.4e1 * t912;
t298 = (-t424 * t767 + t425 * t768) * t1035;
t296 = t1234 / 0.4e1;
t292 = t1235 / 0.4e1;
t275 = 0.2e1 * t910 * t1428;
t261 = -0.4e1 * t424 * t536 + 0.4e1 * t425 * t537;
t260 = -0.4e1 * t1417 + 0.4e1 * t1495;
t231 = 0.4e1 * t1074 * t424 - 0.4e1 * t1409 * t425;
t216 = t1239 / 0.4e1;
t211 = -t1034 * t911 + t1315 * t307;
t201 = 0.4e1 * t424 * t437 + 0.4e1 * t425 * t438;
t178 = t298 * t1307 + t395 * t1309 + t1272 / 0.4e1;
t159 = t1244 / 0.4e1;
t154 = t1307 * t301 + t1309 * t400 - t865;
t152 = t261 * t1307 - t1110 / 0.2e1 - t814 * t707;
t149 = t385 + t1403;
t148 = -t1226 / 0.4e1 + t1085 + t1086;
t147 = t385 - t1403;
t144 = 0.4e1 * t1549;
t143 = t1374 * t1389;
t139 = t1249 / 0.4e1;
t122 = 0.4e1 * t1548;
t120 = 0.4e1 * t1083 * t412 - 0.4e1 * t381 * t382;
t119 = t1307 * t211 + t1309 * t316 + t1311 * t406;
t114 = t379 - t1402;
t113 = t379 + t1402;
t112 = -t1227 / 0.4e1 + t1087 + t1091;
t107 = t1374 * t1325;
t104 = t296 + t159 + t1223 / 0.4e1;
t103 = t439 + t296 - t1244 / 0.4e1;
t102 = t439 + t159 - t1234 / 0.4e1;
t87 = (t1036 * t1542 - t1447 - t1527) * t1307;
t84 = t292 + t139 - t1224 / 0.4e1;
t83 = t433 + t292 - t1249 / 0.4e1;
t82 = t433 + t139 - t1235 / 0.4e1;
t79 = t941 + t1009 - t942;
t78 = t942 - t1407;
t77 = t942 + t1407;
t76 = t1307 * t231 + t1309 * t340 + t1520;
t68 = t1300 / 0.4e1;
t61 = t1307 * t201 + t1309 * t322 - t1334 - t1520;
t44 = t216 + t68 - t1228 / 0.4e1;
t43 = t376 + t216 - t1300 / 0.4e1;
t42 = t376 + t68 - t1239 / 0.4e1;
t35 = t1092 + t1094 - t1090;
t34 = t1090 - t1401;
t33 = t1090 + t1401;
t32 = t1243 / 0.4e1 + t891;
t23 = t1247 / 0.4e1 + t891;
t22 = t1010 + t1404;
t21 = t1011 + t1012 - t1010;
t20 = t1010 - t1404;
t19 = -t144 * t1307 + t302 * t1309 - t841;
t18 = t122 * t1307 + t1309 * t260 + t841;
t17 = -t1200 / 0.2e1 + t120 * t1307 + t1013 * t678 - t1016 * t676;
t16 = t1307 * t143 + t815;
t13 = t107 * t1307 + t767 * t833 + t768 * t832;
t12 = t839 + t928 - t1335 - t1478;
t11 = t859 + t928 + t1458 - t1459;
t10 = t839 + t859 - t928;
t9 = t840 + t929 + t1335 - t1479;
t8 = t840 + t860 - t929;
t7 = -t1408 + t860 + t929 + t1457;
t6 = t850 + t851 + t815 + t1346;
t5 = t1539 * t768 + t1336 + t1534 + t850;
t4 = (t229 / 0.4e1 + t1539) * t768 + t851 + t1347 + t1534;
t3 = t870 + t871 - t1400;
t2 = -t1328 + t1448 + t870;
t1 = -t1329 + t1448 + t871;
t24 = [(m(6) * t1076 * t486 - m(7) * t1081 * t425) * qJD(1) + t154 * qJD(2) + t61 * qJD(3) + t178 * qJD(4) + t76 * qJD(5) + t152 * qJD(6), qJD(1) * t154 + qJD(3) * t77 + qJD(4) * t462 + qJD(5) * t148 + qJD(6) * t103, t61 * qJD(1) + t77 * qJD(2) + t21 * qJD(4) + t6 * qJD(5) + t9 * qJD(6) + (t424 * t453 + t425 * t451 - t437 * t452 - t438 * t450 - t107 / 0.4e1) * t1211 + (t485 * t535 + t486 * t533 - t526 * t534 - t527 * t532) * t1214 + (t555 * t618 + t556 * t615 - t611 * t617 - t612 * t614) * t1216 + (((-t645 / 0.2e1 + t693 / 0.2e1 + t651 / 0.2e1 - t695 / 0.2e1) * t768 + (t646 / 0.2e1 + t694 / 0.2e1 - t652 / 0.2e1 + t696 / 0.2e1) * t767) * t751 + t1532 / 0.2e1 + (t973 * t768 + t1540 - t832) * t768 + (t973 * t767 + t1536 - t833) * t767 + (-t958 / 0.2e1 - t956 / 0.2e1 + t1456 * t1288) * t752 + (-t905 * t726 + (-t697 * t768 + t699 * t767) * t723) * m(4)) * qJD(3), qJD(1) * t178 + qJD(2) * t462 + qJD(3) * t21 + qJD(5) * t113 + qJD(6) * t83, t76 * qJD(1) + t148 * qJD(2) + t6 * qJD(3) + t113 * qJD(4) + t12 * qJD(6) + (-t143 / 0.4e1 + t1074 * t480 - t1409 * t477 + t913) * t1210 + (((-t485 * t604 - t575 * t606) * t967 + t229 / 0.2e1 + t893 + t1540) * t768 + ((-t486 * t604 - t576 * t606) * t967 - t894 + t1536) * t767) * qJD(5), t152 * qJD(1) + t103 * qJD(2) + t9 * qJD(3) + t83 * qJD(4) + t12 * qJD(5) + (-t410 * t536 - t412 * t537 + t424 * t435 + t425 * t436 - t120 / 0.4e1) * t1209 + (-t309 + t1200 / 0.2e1 + (t294 / 0.2e1 + t273 / 0.2e1 - t1013) * t678 - (t272 / 0.2e1 + t293 / 0.2e1 - t1016) * t676) * qJD(6); t865 * qJD(1) + t78 * qJD(3) - t463 * qJD(4) + t147 * qJD(5) + t102 * qJD(6) + (-t301 / 0.4e1 + t1081 * t768) * t1213 + (-t400 / 0.4e1 - t1076 * t768) * t1215, 0, t78 * qJD(1) + ((-m(5) * t615 - m(6) * t533 - m(7) * t451) * t768 + (m(5) * t618 + m(6) * t535 + m(7) * t453) * t767) * qJD(3) + t275 * qJD(5) + qJD(6) * t1027, -t1042, t147 * qJD(1) + t275 * qJD(3) + (-t910 * qJD(5) + qJD(6) * t1295) * m(7), t102 * qJD(1) + (t435 * t767 - t436 * t768) * qJD(6) * t965 + t1093; t79 * qJD(2) + t13 * qJD(3) + t22 * qJD(4) + t5 * qJD(5) + t7 * qJD(6) + (-t201 / 0.4e1 + t1081 * t450) * t1213 + (-t322 / 0.4e1 - t1076 * t532) * t1215 + t829 * t1045 + t1334 * qJD(1) - t1517, t79 * qJD(1) + t1557, t13 * qJD(1) + ((t307 * t333 - t450 * t451 - t452 * t453) * m(7) + (t416 * t428 - t532 * t533 - t534 * t535) * m(6) + m(5) * (t487 * t505 - t614 * t615 - t617 * t618) + m(4) * (t723 * t726 * t935 - 0.4e1 * (t767 * (rSges(4,1) * t1104 - t1049) + t768 * (rSges(4,1) * t1103 + t961)) * t900) / 0.4e1 + t1340 + ((-t1361 * t767 + (t956 + t958 - t1356) * t751 + t1331) * t768 + t1362 * t764) * t1288 + ((-t1362 * t768 + t1331 + ((t1056 + t1058) * t768 - t1356) * t751) * t767 + t1361 * t765) * t1380) * qJD(3) + t119 * qJD(4) + t18 * qJD(5) + t23 * qJD(6), t22 * qJD(1) + t119 * qJD(3) + t33 * qJD(5) + t43 * qJD(6) + (-t459 + t947 * (-0.2e1 * t709 + t1411) * t752) * qJD(4), t5 * qJD(1) + t18 * qJD(3) + t33 * qJD(4) + t2 * qJD(6) + (t144 / 0.4e1 - t1079 * t1497 - t1080 * t1498 - t1089 * t1492) * t1210 + ((-t302 / 0.4e1 + (-t416 + t468) * t1474 - (-t1048 * t606 - t907) * t604) * m(6) + t1340) * qJD(5), t7 * qJD(1) + t1212 * t1295 + t23 * qJD(3) + t43 * qJD(4) + t2 * qJD(5) + ((t1319 * t436 - t1320 * t452 + t1321 * t307 - t1337) * t1308 - t14 - t1330) * qJD(6); (-t1272 / 0.4e1 + (-t298 / 0.4e1 - t1081 * t1106) * m(7) + (-t395 / 0.4e1 + t1076 * t1106) * m(6)) * qJD(1) + t463 * qJD(2) + t20 * qJD(3) + t112 * qJD(5) + t82 * qJD(6), t1042, t20 * qJD(1) + t459 * qJD(4) + t34 * qJD(5) + t42 * qJD(6) + (-t211 / 0.4e1 + (-t333 - t911) * t752 + (t451 * t767 + t453 * t768 + t307) * t751) * t1211 + (-t316 / 0.4e1 + (-t428 - t907) * t752 + (t416 + t906) * t751) * t1214 + (-t406 / 0.4e1 + (-t505 + t903) * t752 + (t615 * t767 + t618 * t768 + t487) * t751) * t1216, t459 * qJD(3), t112 * qJD(1) + t34 * qJD(3) + ((-m(7) * t1492 + t1254) * t752 + (-t1310 * t1412 + t909 * t965) * t751) * qJD(5) + t87 * qJD(6), t82 * qJD(1) + t42 * qJD(3) + t87 * qJD(5) + (t427 * t1314 + (t435 * t768 + t436 * t767) * t1036) * t1209 / 0.2e1; -t340 * t1215 / 0.4e1 + t149 * qJD(2) + t4 * qJD(3) + t114 * qJD(4) + t16 * qJD(5) + t11 * qJD(6) + (-t231 / 0.4e1 - t1081 * t477) * t1213 + (t1076 * t606 * t967 - t829) * t1045 + t1517, t149 * qJD(1) - t1557, t4 * qJD(1) - t1340 * qJD(3) + t35 * qJD(4) + t19 * qJD(5) + t1 * qJD(6) + (-t122 / 0.4e1 + t333 * t341 + t451 * t477 + t453 * t480 + t1548) * t1211 + (t1495 + t468 * t428 - t260 / 0.4e1 + t906 * t606 - t1417) * t1214, qJD(1) * t114 + qJD(3) * t35 + qJD(6) * t88, t16 * qJD(1) + t19 * qJD(3) + (m(7) * t1549 - t468 * t1254 - t1309 * t1422 + t841) * qJD(5) + t32 * qJD(6), t11 * qJD(1) + t1 * qJD(3) + t32 * qJD(5) + ((t1318 * t436 + t1320 * t480 + t1321 * t341 + t1337) * t1308 - t15 + t1330) * qJD(6) - t887; t588 * t1047 / 0.2e1 + t104 * qJD(2) + t8 * qJD(3) + t84 * qJD(4) + t10 * qJD(5) + t17 * qJD(6) + t814 * t1046 + t890 * qJD(1) * t678 + (-t261 / 0.4e1 - t1083 * t425 + t1081 * t412) * t1213, qJD(1) * t104 + t1093, t8 * qJD(1) + t1555 + ((t333 * t382 + t410 * t453 - t412 * t451 - t1512) * t965 + t179 * t844 - t1247 / 0.4e1 + t1577) * qJD(3) + t44 * qJD(4) + t3 * qJD(5) + t14 * qJD(6), qJD(1) * t84 + qJD(3) * t44 - qJD(5) * t88, t10 * qJD(1) + t3 * qJD(3) + (t1377 * t1383 + t1558 * t965 - t1243 / 0.4e1 + t1577) * qJD(5) + t15 * qJD(6) + t887, t17 * qJD(1) + t14 * qJD(3) + t15 * qJD(5) + (t86 * t1293 + t85 * t1294 + (-t272 * t676 + t273 * t678 - t309) * t1292 + (t382 * t427 + t410 * t435 - t412 * t436) * m(7)) * qJD(6);];
Cq  = t24;
