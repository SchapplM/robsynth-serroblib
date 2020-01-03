% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP7_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:11
% EndTime: 2019-12-31 21:57:31
% DurationCPUTime: 66.21s
% Computational Cost: add. (149163->1161), mult. (199923->1637), div. (0->0), fcn. (219379->8), ass. (0->749)
t942 = qJ(2) + qJ(3);
t929 = cos(t942);
t948 = cos(qJ(1));
t1237 = t929 * t948;
t945 = sin(qJ(1));
t1239 = t929 * t945;
t928 = sin(t942);
t1247 = t928 * t948;
t1382 = -t929 / 0.2e1;
t946 = cos(qJ(4));
t1218 = t946 * t948;
t943 = sin(qJ(4));
t1221 = t945 * t943;
t872 = t1218 * t929 + t1221;
t837 = Icges(6,5) * t872;
t1199 = t948 * t943;
t1220 = t945 * t946;
t871 = t1199 * t929 - t1220;
t646 = Icges(6,6) * t1247 + Icges(6,3) * t871 + t837;
t1304 = Icges(6,5) * t871;
t658 = Icges(6,1) * t872 + Icges(6,4) * t1247 + t1304;
t1011 = t646 * t943 + t658 * t946;
t652 = Icges(6,4) * t872 + Icges(6,2) * t1247 + Icges(6,6) * t871;
t1267 = t652 * t929;
t476 = t1011 * t928 - t1267;
t1309 = Icges(5,4) * t872;
t655 = -Icges(5,2) * t871 + Icges(5,6) * t1247 + t1309;
t840 = Icges(5,4) * t871;
t661 = Icges(5,1) * t872 + Icges(5,5) * t1247 - t840;
t1009 = -t655 * t943 + t661 * t946;
t649 = Icges(5,5) * t872 - Icges(5,6) * t871 + Icges(5,3) * t1247;
t1269 = t649 * t929;
t479 = t1009 * t928 - t1269;
t1493 = t476 + t479;
t1019 = Icges(6,5) * t946 + Icges(6,3) * t943;
t780 = Icges(6,6) * t928 + t1019 * t929;
t1303 = Icges(6,5) * t943;
t1026 = Icges(6,1) * t946 + t1303;
t788 = Icges(6,4) * t928 + t1026 * t929;
t1024 = Icges(6,4) * t946 + Icges(6,6) * t943;
t783 = -Icges(6,2) * t929 + t1024 * t928;
t1265 = t783 * t929;
t1250 = t928 * t943;
t1298 = Icges(6,6) * t929;
t1248 = t928 * t946;
t907 = Icges(6,5) * t1248;
t779 = Icges(6,3) * t1250 - t1298 + t907;
t787 = -Icges(6,4) * t929 + t1026 * t928;
t1007 = t779 * t943 + t787 * t946;
t784 = Icges(6,2) * t928 + t1024 * t929;
t987 = t1007 - t784;
t961 = -t928 * t987 + t1265;
t399 = t871 * t780 + t872 * t788 + t948 * t961;
t1307 = Icges(5,4) * t946;
t1025 = -Icges(5,2) * t943 + t1307;
t786 = Icges(5,6) * t928 + t1025 * t929;
t1308 = Icges(5,4) * t943;
t1027 = Icges(5,1) * t946 - t1308;
t790 = Icges(5,5) * t928 + t1027 * t929;
t1020 = Icges(5,5) * t946 - Icges(5,6) * t943;
t781 = -Icges(5,3) * t929 + t1020 * t928;
t1266 = t781 * t929;
t785 = -Icges(5,6) * t929 + t1025 * t928;
t789 = -Icges(5,5) * t929 + t1027 * t928;
t1006 = -t785 * t943 + t789 * t946;
t782 = Icges(5,3) * t928 + t1020 * t929;
t986 = -t1006 + t782;
t962 = t928 * t986 + t1266;
t400 = -t871 * t786 + t872 * t790 + t948 * t962;
t1495 = t399 + t400;
t1522 = t928 / 0.2e1;
t1249 = t928 * t945;
t869 = t1221 * t929 + t1218;
t870 = t1220 * t929 - t1199;
t650 = Icges(6,4) * t870 + Icges(6,2) * t1249 + Icges(6,6) * t869;
t1268 = t650 * t929;
t836 = Icges(6,5) * t870;
t645 = -Icges(6,6) * t1249 - Icges(6,3) * t869 - t836;
t835 = Icges(6,5) * t869;
t656 = Icges(6,1) * t870 + Icges(6,4) * t1249 + t835;
t1525 = t645 * t943 - t656 * t946;
t475 = t1525 * t928 + t1268;
t647 = Icges(5,5) * t870 - Icges(5,6) * t869 + Icges(5,3) * t1249;
t1270 = t647 * t929;
t839 = Icges(5,4) * t870;
t653 = -Icges(5,2) * t869 + Icges(5,6) * t1249 + t839;
t838 = Icges(5,4) * t869;
t660 = -Icges(5,1) * t870 - Icges(5,5) * t1249 + t838;
t1524 = t653 * t943 + t660 * t946;
t478 = t1524 * t928 + t1270;
t1573 = t478 + t475;
t971 = -t1019 * t928 + t1298;
t731 = t971 * t948;
t739 = t787 * t948;
t990 = t783 * t948 + t1011;
t355 = t990 * t929 + (t731 * t943 - t739 * t946 + t652) * t928;
t737 = t785 * t948;
t741 = t789 * t948;
t988 = -t781 * t948 - t1009;
t357 = -t988 * t929 + (t737 * t943 - t741 * t946 + t649) * t928;
t419 = t987 * t929 + (t780 * t943 + t788 * t946 + t783) * t928;
t420 = -t986 * t929 + (-t786 * t943 + t790 * t946 + t781) * t928;
t519 = t1249 * t783 + t779 * t869 + t787 * t870;
t520 = t1249 * t781 - t785 * t869 + t789 * t870;
t523 = t1247 * t783 + t871 * t779 + t872 * t787;
t524 = t1247 * t781 - t871 * t785 + t872 * t789;
t548 = t1007 * t928 - t1265;
t549 = t1006 * t928 - t1266;
t1579 = (t548 + t549) * t1522 + (t419 + t420) * t1382 + (t519 + t520 - t1573) * t1239 / 0.4e1 + (t523 + t524 + t1493) * t1237 / 0.4e1 + t1247 * (t355 + t357 + t1495) / 0.4e1;
t1328 = m(6) * qJD(1);
t1578 = -t1328 / 0.4e1;
t1399 = -pkin(7) - pkin(6);
t947 = cos(qJ(2));
t1337 = pkin(2) * t947;
t924 = pkin(1) + t1337;
t1124 = -t948 * t1399 - t945 * t924;
t1336 = pkin(3) * t929;
t1374 = rSges(6,2) + pkin(8);
t1458 = t1374 * t928;
t1375 = rSges(6,1) + pkin(4);
t1513 = rSges(6,3) + qJ(5);
t1544 = t1375 * t870 + t1513 * t869;
t525 = -t1544 + (-t1336 - t1458) * t945 + t1124;
t1062 = t924 + t1336;
t1553 = t1375 * t872 + t1513 * t871;
t925 = t945 * t1399;
t526 = -t925 + (t1062 + t1458) * t948 + t1553;
t363 = -0.4e1 * t525 * t869 + 0.4e1 * t526 * t871;
t1577 = t363 * t1328;
t1159 = -t1375 * t869 + t1513 * t870;
t1406 = -0.2e1 * t948;
t863 = (-rSges(5,1) * t943 - rSges(5,2) * t946) * t928;
t1117 = t863 * t1406;
t1262 = t863 * t945;
t1373 = rSges(5,3) + pkin(8);
t1459 = t1373 * t928;
t1461 = -t870 * rSges(5,1) + t869 * rSges(5,2);
t572 = t1461 + (-t1336 - t1459) * t945 + t1124;
t1034 = t872 * rSges(5,1) - t871 * rSges(5,2);
t573 = -t925 + (t1062 + t1459) * t948 + t1034;
t1178 = t572 * t1117 - 0.2e1 * t573 * t1262;
t1128 = (-t1375 * t943 + t1513 * t946) * t928;
t677 = t1128 * t948;
t1416 = -0.2e1 * t677;
t1438 = 0.2e1 * t526;
t675 = t1128 * t945;
t1183 = t525 * t1416 - t675 * t1438;
t1400 = m(6) / 0.4e1;
t1402 = m(5) / 0.4e1;
t704 = -rSges(5,1) * t871 - rSges(5,2) * t872;
t1413 = -0.2e1 * t704;
t699 = -rSges(5,1) * t869 - rSges(5,2) * t870;
t1414 = -0.2e1 * t699;
t1543 = t1375 * t946 + t1513 * t943;
t1145 = -rSges(6,2) * t929 + t1543 * t928;
t888 = pkin(3) * t928 - pkin(8) * t929;
t1096 = t888 + t1145;
t595 = t1096 * t948;
t1426 = -0.2e1 * t595;
t593 = t1096 * t945;
t1427 = -0.2e1 * t593;
t577 = -t1375 * t871 + t1513 * t872;
t1321 = rSges(5,1) * t946;
t1033 = -rSges(5,2) * t943 + t1321;
t1320 = rSges(5,3) * t929;
t792 = t1033 * t928 - t1320;
t1144 = t792 + t888;
t671 = t1144 * t945;
t673 = t1144 * t948;
t1463 = (-t1159 * t1426 + t1427 * t577 + t1183) * t1400 + (t1413 * t671 - t1414 * t673 + t1178) * t1402;
t1432 = 0.2e1 * t573;
t663 = rSges(5,3) * t1249 - t1461;
t1112 = rSges(5,2) * t1250;
t1107 = t928 * t1220;
t893 = rSges(5,1) * t1107;
t743 = -t893 + (t1112 + t1320) * t945;
t794 = rSges(5,3) * t928 + t1033 * t929;
t483 = (t792 * t945 + t743) * t929 + (t794 * t945 - t663) * t928;
t667 = rSges(5,3) * t1247 + t1034;
t1108 = t928 * t1199;
t1125 = rSges(5,2) * t1108 + rSges(5,3) * t1237;
t745 = -rSges(5,1) * t1218 * t928 + t1125;
t484 = (-t792 * t948 - t745) * t929 + (-t794 * t948 + t667) * t928;
t1189 = t484 * t1432 + 0.2e1 * t483 * t572;
t1143 = rSges(6,2) * t928 + t1543 * t929;
t1109 = t928 * t1221;
t1545 = t1375 * t1107 + t1513 * t1109;
t1158 = rSges(6,2) * t1239 - t1545;
t1534 = rSges(6,2) * t1249 + t1544;
t360 = (t1145 * t945 + t1158) * t929 + (t1143 * t945 - t1534) * t928;
t1160 = rSges(6,2) * t1247 + t1553;
t913 = rSges(6,2) * t1237;
t1533 = -t1543 * t1247 + t913;
t361 = (-t1145 * t948 - t1533) * t929 + (-t1143 * t948 + t1160) * t928;
t1198 = t361 * t1438 + 0.2e1 * t360 * t525;
t1039 = (-pkin(3) - t1321) * t928;
t915 = pkin(8) * t1237;
t1094 = t915 + t1125;
t643 = t1039 * t948 + t1094;
t1420 = -0.2e1 * t643;
t914 = pkin(3) * t1249;
t1126 = t893 + t914;
t982 = -t1373 * t929 - t1112;
t642 = t945 * t982 + t1126;
t1421 = 0.2e1 * t642;
t1030 = t914 + t1545;
t1093 = t1374 * t929;
t587 = -t1093 * t945 + t1030;
t1428 = 0.2e1 * t587;
t1123 = t913 + t915;
t969 = (-pkin(3) - t1543) * t928;
t586 = t948 * t969 + t1123;
t1429 = -0.2e1 * t586;
t1542 = t1145 * t1249 + t1534 * t929;
t489 = t1145 * t1247 + t1160 * t929;
t558 = t792 * t1249 + t663 * t929;
t560 = t1247 * t792 + t929 * t667;
t1471 = (t1420 * t560 + t1421 * t558 + t1189) * t1402 + (t1428 * t1542 + t1429 * t489 + t1198) * t1400;
t1576 = t1471 - t1463;
t944 = sin(qJ(2));
t1338 = pkin(2) * t944;
t1064 = t888 + t1338;
t1002 = t1064 + t1145;
t584 = t1002 * t948;
t1430 = -0.2e1 * t584;
t582 = t1002 * t945;
t1431 = -0.2e1 * t582;
t1042 = t1064 + t792;
t629 = t1042 * t945;
t631 = t1042 * t948;
t1464 = (-t1159 * t1430 + t1431 * t577 + t1183) * t1400 + (t1413 * t629 - t1414 * t631 + t1178) * t1402;
t1434 = -0.2e1 * t560;
t1439 = -0.2e1 * t489;
t569 = (t969 - t1338) * t948 + t1123;
t570 = (-t1093 + t1338) * t945 + t1030;
t617 = (t982 + t1338) * t945 + t1126;
t618 = (t1039 - t1338) * t948 + t1094;
t1472 = (t1434 * t618 + 0.2e1 * t558 * t617 + t1189) * t1402 + (t1439 * t569 + 0.2e1 * t1542 * t570 + t1198) * t1400;
t1575 = t1472 - t1464;
t435 = t1249 * t650 - t645 * t869 + t656 * t870;
t436 = t652 * t1249 + t869 * t646 + t870 * t658;
t437 = t1249 * t647 - t653 * t869 - t660 * t870;
t438 = t649 * t1249 - t869 * t655 + t870 * t661;
t1539 = (-t438 - t436) * t945 + (t435 + t437) * t948;
t687 = -Icges(5,5) * t869 - Icges(5,6) * t870;
t689 = -Icges(6,4) * t869 + Icges(6,6) * t870;
t1572 = t687 + t689;
t688 = -Icges(5,5) * t871 - Icges(5,6) * t872;
t690 = -Icges(6,4) * t871 + Icges(6,6) * t872;
t1571 = t688 + t690;
t1163 = Icges(5,2) * t872 - t661 + t840;
t1165 = Icges(6,3) * t872 - t1304 - t658;
t1570 = t1163 + t1165;
t1164 = Icges(5,2) * t870 + t660 + t838;
t1166 = Icges(6,3) * t870 - t656 - t835;
t1569 = t1164 + t1166;
t1167 = -Icges(5,1) * t871 - t1309 - t655;
t1169 = -Icges(6,1) * t871 + t646 + t837;
t1568 = t1167 + t1169;
t1168 = -Icges(5,1) * t869 - t653 - t839;
t1170 = -Icges(6,1) * t869 - t645 + t836;
t1567 = t1168 + t1170;
t841 = (Icges(6,3) * t946 - t1303) * t928;
t842 = (-Icges(5,5) * t943 - Icges(5,6) * t946) * t928;
t845 = (-Icges(6,4) * t943 + Icges(6,6) * t946) * t928;
t846 = (-Icges(5,2) * t946 - t1308) * t928;
t849 = -Icges(6,1) * t1250 + t907;
t850 = (-Icges(5,1) * t943 - t1307) * t928;
t1566 = -(-(t789 / 0.2e1 + t846 / 0.2e1 + t787 / 0.2e1 - t841 / 0.2e1) * t943 + (t850 / 0.2e1 - t785 / 0.2e1 + t849 / 0.2e1 + t779 / 0.2e1) * t946) * t928 + (t842 / 0.2e1 + t845 / 0.2e1) * t929;
t1017 = t437 * t945 + t438 * t948;
t82 = -t1017 * t928 + t520 * t929;
t1018 = t435 * t945 + t436 * t948;
t81 = -t1018 * t928 + t519 * t929;
t439 = t650 * t1247 - t871 * t645 + t872 * t656;
t440 = t652 * t1247 + t871 * t646 + t872 * t658;
t1016 = t945 * t439 + t440 * t948;
t1554 = t1016 * t928 - t523 * t929;
t441 = t647 * t1247 - t871 * t653 - t872 * t660;
t442 = t649 * t1247 - t871 * t655 + t872 * t661;
t1015 = t945 * t441 + t442 * t948;
t1555 = t1015 * t928 - t524 * t929;
t1535 = t1554 + t1555;
t1540 = -m(5) * qJD(1) / 0.4e1;
t1551 = (-t439 - t441) * t948 + (t440 + t442) * t945;
t1563 = t1572 * t1249 + t1567 * t870 + t1569 * t869;
t1562 = t1571 * t1249 + t1568 * t870 + t1570 * t869;
t1561 = t1572 * t1247 + t1567 * t872 + t1569 * t871;
t1560 = t1571 * t1247 + t1568 * t872 + t1570 * t871;
t1087 = t1249 / 0.4e1;
t1089 = -t1249 / 0.4e1;
t1552 = (t1087 + t1089) * t1551;
t1379 = t945 / 0.2e1;
t1376 = t948 / 0.2e1;
t1147 = -t787 + t841;
t1149 = t779 + t849;
t445 = t1147 * t869 + t1149 * t870 + t1249 * t845;
t1146 = -t789 - t846;
t1148 = -t785 + t850;
t446 = t1146 * t869 + t1148 * t870 + t1249 * t842;
t1550 = t445 + t446;
t447 = t1147 * t871 + t1149 * t872 + t1247 * t845;
t448 = t1146 * t871 + t1148 * t872 + t1247 * t842;
t1549 = t447 + t448;
t1541 = m(5) / 0.2e1;
t1521 = t929 / 0.2e1;
t1538 = t1562 * t945 - t1563 * t948;
t1537 = t1560 * t945 - t1561 * t948;
t921 = Icges(4,4) * t929;
t883 = -Icges(4,2) * t928 + t921;
t884 = Icges(4,1) * t928 + t921;
t1532 = t883 + t884;
t940 = t945 ^ 2;
t941 = t948 ^ 2;
t1121 = t940 + t941;
t802 = Icges(4,5) * t1239 - Icges(4,6) * t1249 - Icges(4,3) * t948;
t805 = Icges(4,6) * t945 + t883 * t948;
t1049 = t805 * t928 - t802;
t1310 = Icges(4,4) * t928;
t885 = Icges(4,1) * t929 - t1310;
t807 = Icges(4,5) * t945 + t885 * t948;
t756 = t807 * t1239;
t881 = Icges(4,5) * t929 - Icges(4,6) * t928;
t1256 = t881 * t948;
t803 = Icges(4,3) * t945 + t1256;
t1051 = t948 * t803 - t756;
t1154 = t807 * t1237 + t945 * t803;
t908 = Icges(4,4) * t1249;
t806 = Icges(4,1) * t1239 - Icges(4,5) * t948 - t908;
t1155 = -t806 * t1237 - t945 * t802;
t804 = Icges(4,4) * t1239 - Icges(4,2) * t1249 - Icges(4,6) * t948;
t1264 = t804 * t928;
t1377 = -t948 / 0.2e1;
t1481 = t1551 * t1376 - t1379 * t1539;
t554 = -t1249 * t805 - t1051;
t555 = -t1247 * t804 - t1155;
t556 = -t1247 * t805 + t1154;
t960 = (-t555 * t948 + t556 * t945) * t1376 + ((t554 - t756 + (t803 + t1264) * t948 + t1155) * t948 + t1154 * t945 + t1551) * t1377 + ((t1049 * t945 + t1051 + t554 + t555) * t945 + ((-t806 * t929 + t1264) * t945 - t1154 + t556 + (t802 + t1049) * t948) * t948 + t1539) * t1379 + t1481;
t1404 = m(4) / 0.4e1;
t1523 = -t928 / 0.2e1;
t1519 = -t948 / 0.4e1;
t730 = t971 * t945;
t738 = t787 * t945;
t991 = t783 * t945 - t1525;
t964 = -t928 * t991 + t1268;
t322 = t730 * t869 - t738 * t870 + t945 * t964;
t963 = -t928 * t990 + t1267;
t323 = t731 * t869 - t739 * t870 + t945 * t963;
t397 = t780 * t869 + t788 * t870 + t945 * t961;
t65 = (t1018 - t397) * t929 + (t322 * t945 + t323 * t948 + t519) * t928;
t736 = t785 * t945;
t740 = t789 * t945;
t989 = -t781 * t945 + t1524;
t966 = t928 * t989 + t1270;
t324 = t736 * t869 - t740 * t870 + t945 * t966;
t965 = t928 * t988 + t1269;
t325 = t737 * t869 - t741 * t870 + t945 * t965;
t398 = -t786 * t869 + t790 * t870 + t945 * t962;
t66 = (t1017 - t398) * t929 + (t324 * t945 + t325 * t948 + t520) * t928;
t1515 = t65 + t66;
t326 = t871 * t730 - t872 * t738 + t948 * t964;
t327 = t871 * t731 - t872 * t739 + t948 * t963;
t67 = (t1016 - t399) * t929 + (t326 * t945 + t327 * t948 + t523) * t928;
t328 = t871 * t736 - t872 * t740 + t948 * t966;
t329 = t871 * t737 - t872 * t741 + t948 * t965;
t68 = (t1015 - t400) * t929 + (t328 * t945 + t329 * t948 + t524) * t928;
t1514 = t68 + t67;
t1324 = m(6) * qJD(5);
t1501 = -t1550 * t929 + (t1562 * t948 + t1563 * t945) * t928;
t1500 = -t1549 * t929 + (t1560 * t948 + t1561 * t945) * t928;
t1499 = (-t322 - t324) * t948 + (t323 + t325) * t945;
t1498 = (-t326 - t328) * t948 + (t327 + t329) * t945;
t1496 = t397 + t398;
t1322 = rSges(4,1) * t929;
t887 = -rSges(4,2) * t928 + t1322;
t1255 = t887 * t945;
t1059 = -rSges(4,2) * t1247 + t945 * rSges(4,3);
t707 = -t925 + (t924 + t1322) * t948 + t1059;
t625 = -0.2e1 * t707 * t1255;
t1116 = t887 * t1406;
t808 = rSges(4,1) * t1239 - rSges(4,2) * t1249 - t948 * rSges(4,3);
t706 = -t808 + t1124;
t626 = t706 * t1116;
t1171 = t625 + t626;
t889 = pkin(8) * t928 + t1336;
t1142 = -t794 - t889;
t672 = t1142 * t945;
t513 = t672 * t1432;
t674 = t1142 * t948;
t1417 = 0.2e1 * t674;
t514 = t572 * t1417;
t1180 = t513 + t514;
t1095 = -t889 - t1143;
t594 = t1095 * t945;
t450 = t594 * t1438;
t596 = t1095 * t948;
t1425 = 0.2e1 * t596;
t451 = t525 * t1425;
t1188 = t450 + t451;
t1418 = -0.2e1 * t673;
t1419 = -0.2e1 * t671;
t886 = rSges(4,1) * t928 + rSges(4,2) * t929;
t981 = t886 + t1338;
t1465 = t981 * t948;
t1466 = t981 * t945;
t1104 = (t1426 * t570 + t1427 * t569 + t1188) * t1400 + (t1418 * t617 + t1419 * t618 + t1180) * t1402 + (0.2e1 * (t1465 * t945 - t1466 * t948) * t886 + t1171) * t1404;
t865 = t886 * t945;
t867 = t886 * t948;
t1105 = (-t1428 * t584 + t1429 * t582 + t1188) * t1400 + (t1420 * t629 - t1421 * t631 + t1180) * t1402 + (-0.2e1 * t1465 * t865 + 0.2e1 * t1466 * t867 + t1171) * t1404;
t1485 = t1104 - t1105;
t977 = -t1160 * t945 + t1534 * t948;
t406 = t977 * t928;
t1442 = 0.2e1 * t406;
t1137 = t945 * (pkin(8) * t1239 - t914) + t948 * (-pkin(3) * t1247 + t915);
t1453 = t1158 * t945 + t1533 * t948;
t470 = t1137 + t1453;
t1103 = t1425 * t1542 + t594 * t1439 + t470 * t1442;
t310 = t977 * t929 + (t1158 * t948 - t1533 * t945) * t928;
t1444 = 0.2e1 * t310;
t1136 = t1121 * t889;
t394 = t1160 * t948 + t1534 * t945 + t1136;
t35 = t1426 * t360 + t1427 * t361 + t1444 * t394 + t1103;
t1008 = t948 * t663 - t667 * t945;
t531 = t1008 * t928;
t1436 = 0.2e1 * t531;
t1462 = t945 * t743 + t948 * t745;
t528 = t1137 + t1462;
t1100 = t558 * t1417 + t672 * t1434 + t528 * t1436;
t422 = t1008 * t929 + (t743 * t948 - t745 * t945) * t928;
t1441 = 0.2e1 * t422;
t491 = t945 * t663 + t948 * t667 + t1136;
t95 = t1418 * t483 + t1419 * t484 + t1441 * t491 + t1100;
t1334 = t35 * t1400 + t95 * t1402;
t939 = t948 * pkin(6);
t1152 = -t945 * (pkin(1) * t945 + t1124 - t939) + t948 * (-t945 * pkin(6) - t925 + (-pkin(1) + t924) * t948);
t358 = t394 + t1152;
t182 = t358 * t1444;
t413 = t491 + t1152;
t266 = t413 * t1441;
t308 = t360 * t1430;
t309 = t361 * t1431;
t1423 = -0.2e1 * t631;
t403 = t483 * t1423;
t1424 = -0.2e1 * t629;
t404 = t484 * t1424;
t1335 = (t308 + t309 + t182 + t1103) * t1400 + (t403 + t404 + t266 + t1100) * t1402;
t1484 = t1334 - t1335;
t1405 = m(3) / 0.4e1;
t1323 = rSges(3,1) * t947;
t1066 = pkin(1) + t1323;
t1236 = t944 * t945;
t1122 = rSges(3,2) * t1236 + t948 * rSges(3,3);
t759 = -t1066 * t945 + t1122 + t939;
t1235 = t944 * t948;
t920 = rSges(3,2) * t1235;
t760 = -t920 + t1066 * t948 + (rSges(3,3) + pkin(6)) * t945;
t901 = rSges(3,1) * t944 + rSges(3,2) * t947;
t879 = t901 * t945;
t880 = t901 * t948;
t1311 = Icges(3,4) * t944;
t897 = Icges(3,2) * t947 + t1311;
t900 = Icges(3,1) * t947 - t1311;
t1483 = -(t900 / 0.2e1 - t897 / 0.2e1) * t944 - 0.4e1 * (-t1465 * t707 + t1466 * t706) * t1404 - 0.4e1 * (t759 * t879 - t760 * t880) * t1405;
t354 = t991 * t929 + (t730 * t943 - t738 * t946 + t650) * t928;
t356 = -t989 * t929 + (t736 * t943 - t740 * t946 + t647) * t928;
t1480 = t354 + t356 + t1496;
t882 = Icges(4,2) * t929 + t1310;
t1478 = t882 * t1523 + (t782 + t784) * t1382 + (t781 + t783 + t885) * t1522 + (t790 + t788) * t1248 / 0.2e1 + (t1532 + (t789 + t787) * t946) * t1521;
t936 = Icges(3,4) * t947;
t898 = -Icges(3,2) * t944 + t936;
t899 = Icges(3,1) * t944 + t936;
t627 = t945 * t808 + t948 * (rSges(4,1) * t1237 + t1059);
t676 = -t945 * t865 - t948 * t867;
t1457 = t1121 * t887 * t886 + t627 * t676;
t1455 = m(5) * (t422 * t531 + t483 * t558 - t484 * t560) + m(6) * (t1542 * t360 + t310 * t406 - t361 * t489);
t1046 = 0.2e1 * t358 + 0.2e1 * t394;
t471 = t1159 * t945 + t577 * t948;
t567 = t945 * t699 + t704 * t948;
t1454 = ((t413 + t491) * t567 + ((t631 + t673) * t948 + (t629 + t671) * t945) * t863) * t1541 + (t471 * t1046 + (-t584 - t595) * t1416 - 0.2e1 * (-t582 - t593) * t675) * t1400;
t1037 = t413 * t528 - t629 * t672 - t631 * t674;
t1038 = t358 * t470 - t582 * t594 - t584 * t596;
t1047 = t1538 * t1377 + t1537 * t1379;
t832 = Icges(3,5) * t945 + t900 * t948;
t1131 = -t897 * t948 + t832;
t1219 = t945 * t947;
t918 = Icges(3,4) * t1236;
t831 = Icges(3,1) * t1219 - Icges(3,5) * t948 - t918;
t1132 = -Icges(3,2) * t1219 + t831 - t918;
t830 = Icges(3,6) * t945 + t898 * t948;
t1133 = -t899 * t948 - t830;
t829 = Icges(3,4) * t1219 - Icges(3,2) * t1236 - Icges(3,6) * t948;
t1134 = t899 * t945 + t829;
t1449 = (-t1131 * t945 + t1132 * t948) * t944 + (t1133 * t945 + t1134 * t948) * t947;
t1138 = -t882 * t948 + t807;
t1139 = -Icges(4,2) * t1239 + t806 - t908;
t1140 = -t884 * t948 - t805;
t1141 = t884 * t945 + t804;
t1448 = (-t1138 * t945 + t1139 * t948) * t928 + (t1140 * t945 + t1141 * t948) * t929;
t376 = -t688 * t929 + (t1163 * t943 + t1167 * t946) * t928;
t1287 = t376 * t945;
t375 = -t687 * t929 + (t1164 * t943 + t1168 * t946) * t928;
t1288 = t375 * t948;
t374 = -t690 * t929 + (t1165 * t943 + t1169 * t946) * t928;
t1289 = t374 * t945;
t373 = -t689 * t929 + (t1166 * t943 + t1170 * t946) * t928;
t1290 = t373 * t948;
t978 = -t1290 / 0.4e1 + t1289 / 0.4e1 - t1288 / 0.4e1 + t1287 / 0.4e1 + t1549 * t945 / 0.4e1 + t1550 * t1519;
t1445 = (t1519 + t948 / 0.4e1) * t1535;
t449 = (t1159 * t948 - t577 * t945) * t928;
t1440 = 0.2e1 * t449;
t550 = (t699 * t948 - t704 * t945) * t928;
t1435 = 0.2e1 * t550;
t1433 = 0.4e1 * t567;
t633 = (t869 * t948 - t871 * t945) * t928;
t1422 = 0.2e1 * t633;
t682 = t945 * t869 + t871 * t948;
t1415 = 0.4e1 * t682;
t927 = t928 ^ 2;
t728 = t1221 * t927 + t869 * t929;
t1412 = 0.2e1 * t728;
t729 = -t1199 * t927 - t929 * t871;
t1411 = -0.2e1 * t729;
t1410 = 0.2e1 * t869;
t1409 = -0.2e1 * t870;
t1408 = 0.2e1 * t871;
t1407 = 0.2e1 * t872;
t1401 = m(6) / 0.2e1;
t1052 = -0.2e1 * t1108;
t1053 = 0.2e1 * t1109;
t1186 = t1052 * t1542 + t489 * t1053;
t1390 = m(6) * (t360 * t1408 + t361 * t1410 + 0.2e1 * (t310 * t928 + t406 * t929) * t943 + t1186);
t181 = 0.4e1 * t1542 * t728 + 0.4e1 * t406 * t633 - 0.4e1 * t489 * t729;
t1385 = t181 / 0.4e1;
t1380 = -t945 / 0.2e1;
t1369 = 0.4e1 * m(4) * (t706 * t865 - t707 * t867);
t1101 = t682 * t1442 + t1186;
t1358 = m(6) * (t1411 * t582 - t1412 * t584 + t1422 * t358 + t1101);
t1357 = m(6) * (t1411 * t593 - t1412 * t595 + t1422 * t394 + t1101);
t1172 = -t595 * t1052 + t593 * t1053;
t565 = t582 * t1053;
t566 = t584 * t1052;
t1174 = t565 - t566;
t1356 = m(6) * (t1046 * t682 + t1172 + t1174);
t1259 = t871 * t489;
t1187 = -0.2e1 * t1542 * t869 - 0.2e1 * t1259;
t1355 = m(6) * (t1410 * t1542 + t1187 + 0.2e1 * t1259);
t1115 = 0.2e1 * t1250;
t1098 = t470 * t1115 + t596 * t1408 + t594 * t1410;
t1240 = t929 * t943;
t1113 = 0.2e1 * t1240;
t347 = t358 * t1113;
t1354 = m(6) * (t347 + t1098 + t1174);
t1097 = t471 * t1115 - t677 * t1408 - t675 * t1410;
t1114 = 0.2e1 * t1248;
t1353 = m(6) * (t1114 * t358 - t1407 * t584 + t1409 * t582 + t1097);
t154 = t1113 * t394 + t1098 + t1172;
t1352 = m(6) * t154;
t1351 = m(6) * (t1114 * t394 - t1407 * t595 + t1409 * t593 + t1097);
t1258 = t871 * t582;
t1261 = t869 * t584;
t1177 = -0.2e1 * t1258 + 0.2e1 * t1261;
t1348 = m(6) * (t1177 + 0.2e1 * t1258 - 0.2e1 * t1261);
t1347 = m(6) * (t1412 * t525 + t1438 * t729 + t1187);
t1257 = t871 * t593;
t1260 = t869 * t595;
t1175 = -0.2e1 * t1257 + 0.2e1 * t1260;
t1344 = m(6) * (t1175 + 0.2e1 * t1257 - 0.2e1 * t1260);
t1343 = 0.2e1 * (-t1159 * t871 + t525 * t872 + t526 * t870 + t577 * t869) * m(6);
t1181 = t525 * t1052 - 0.2e1 * t526 * t1109;
t1342 = m(6) * (t1408 * t570 + t1410 * t569 + t1181);
t1341 = m(6) * (t1177 + t1181);
t1340 = m(6) * (t1408 * t587 + t1410 * t586 + t1181);
t1339 = m(6) * (t1175 + t1181);
t1333 = m(4) * qJD(2);
t1332 = m(4) * qJD(3);
t1330 = m(5) * qJD(2);
t1329 = m(5) * qJD(3);
t1327 = m(6) * qJD(2);
t1326 = m(6) * qJD(3);
t1325 = m(6) * qJD(4);
t1295 = t354 * t948;
t1294 = t355 * t945;
t1293 = t356 * t948;
t1292 = t357 * t945;
t1283 = t475 * t945;
t1282 = t478 * t945;
t538 = t627 + t1152;
t1279 = t538 * t676;
t1263 = t829 * t944;
t1217 = t947 * t948;
t827 = Icges(3,5) * t1219 - Icges(3,6) * t1236 - Icges(3,3) * t948;
t1151 = -t831 * t1217 - t945 * t827;
t1023 = Icges(3,5) * t947 - Icges(3,6) * t944;
t828 = Icges(3,3) * t945 + t1023 * t948;
t1150 = t832 * t1217 + t945 * t828;
t1120 = qJD(4) * t928;
t1118 = 0.4e1 * t863;
t1106 = t1324 / 0.4e1;
t1102 = t1416 * t1542 - t675 * t1439 + t471 * t1442;
t1099 = t558 * t1117 + 0.2e1 * t560 * t1262 + t567 * t1436;
t1092 = -t1250 / 0.2e1;
t1091 = t1250 / 0.2e1;
t1088 = t1249 / 0.2e1;
t1082 = t1247 / 0.2e1;
t1080 = -t1240 / 0.2e1;
t1079 = t1240 / 0.2e1;
t1065 = -t887 - t1337;
t1063 = -t889 - t1337;
t1054 = 0.4e1 * t1250;
t774 = t832 * t1219;
t1050 = t948 * t828 - t774;
t1048 = t830 * t944 - t827;
t1043 = t1121 * t1338;
t1041 = t1063 - t794;
t1036 = t1279 + t1466 * t1255 - t1465 * t1116 / 0.2e1;
t1022 = -Icges(3,5) * t944 - Icges(3,6) * t947;
t1021 = Icges(4,5) * t928 + Icges(4,6) * t929;
t1014 = t476 * t948 - t1283;
t1013 = t479 * t948 - t1282;
t1005 = t804 * t929 + t806 * t928;
t1001 = t1063 - t1143;
t999 = t1047 + t1454;
t843 = t1021 * t945;
t844 = t948 * t1021;
t996 = (-t940 * t844 + (t945 * t843 + t1448) * t948 + t1498) * t1379 + (-t941 * t843 + (t948 * t844 + t1448) * t945 + t1499) * t1377;
t983 = -t1043 + t1137;
t980 = t394 * t470 - t593 * t594 - t595 * t596;
t979 = t491 * t528 - t671 * t672 - t673 * t674;
t976 = (-t882 + t885) * t929 - t1532 * t928;
t958 = t779 * t1079 + t785 * t1080 + t780 * t1091 + t786 * t1092 + t1478;
t957 = t1445 + t1552;
t956 = -t1047 + (t1493 * t945 + t1573 * t948) * t1522 + (t1294 - t1295 + t1292 - t1293) * t1382 + t1514 * t1379 + t1515 * t1377 + t1499 * t1088 + t1498 * t1082 + t1481 * t929;
t955 = t1087 * t1480 + t1579;
t250 = t1014 * t928 - t548 * t929;
t251 = t1013 * t928 - t549 * t929;
t92 = (t1014 - t419) * t929 + (t354 * t945 + t355 * t948 + t548) * t928;
t93 = (t1013 - t420) * t929 + (t356 * t945 + t357 * t948 + t549) * t928;
t954 = -t1455 + (t251 + t250) * t1523 + (t1289 - t1290 + t1287 - t1288) * t1382 + (t93 + t92) * t1521 + t1500 * t1379 + t1501 * t1377 - t1515 * t1249 / 0.2e1 + t1538 * t1088 - t1514 * t1247 / 0.2e1 + t1537 * t1082 - (-t82 - t81) * t1239 / 0.2e1 - t1535 * t1237 / 0.2e1;
t953 = t1292 / 0.2e1 - t1293 / 0.2e1 + t1294 / 0.2e1 - t1295 / 0.2e1 + (t1138 * t929 + t1140 * t928 + t945 * t881 + t948 * t976 + t1495) * t1379 + (t1139 * t929 - t1141 * t928 + t945 * t976 - t1256 + t1496) * t1377 - t960;
t952 = t786 * t1091 + t780 * t1092 + t785 * t1079 + t779 * t1080 - t1283 / 0.2e1 + t1005 * t1380 - t1282 / 0.2e1 + (t1005 + t1573) * t1379 - t1478;
t951 = t1089 * t1480 - t1579 + t957 + t978;
t950 = -t1445 + t1552 + t955 + t978;
t949 = t955 + t957 - t978;
t903 = -rSges(3,2) * t944 + t1323;
t874 = t1022 * t948;
t873 = t1022 * t945;
t798 = t1065 * t948;
t796 = t1065 * t945;
t632 = t1041 * t948;
t630 = t1041 * t945;
t616 = -t1043 + t676;
t592 = (-t682 + t1240) * t1054;
t590 = -t1247 * t863 - t929 * t704;
t589 = t1249 * t863 + t699 * t929;
t588 = t592 * t1106;
t585 = t1001 * t948;
t583 = t1001 * t945;
t581 = -t1235 * t830 + t1150;
t580 = -t1235 * t829 - t1151;
t579 = -t1236 * t830 - t1050;
t568 = 0.4e1 * t927 * t943 * t946 + 0.4e1 * t869 * t870 + 0.4e1 * t871 * t872;
t509 = t983 + t1462;
t498 = -t577 * t929 - t677 * t928;
t497 = t1128 * t1249 + t1159 * t929;
t494 = -t929 * t842 + (t1146 * t943 + t1148 * t946) * t928;
t493 = -t929 * t845 + (t1147 * t943 + t1149 * t946) * t928;
t486 = -t580 * t948 + t581 * t945;
t485 = -(-(-t831 * t947 + t1263) * t945 - t948 * t827) * t948 + t579 * t945;
t482 = 0.4e1 * t1457;
t453 = -t592 * t1324 / 0.4e1;
t452 = t983 + t1453;
t381 = 0.4e1 * t1279 + 0.4e1 * (t1465 * t948 + t1466 * t945) * t887;
t378 = -0.4e1 * t572 * t699 + 0.4e1 * t573 * t704;
t364 = 0.4e1 * t572 * t642 + 0.4e1 * t573 * t643;
t362 = 0.4e1 * t572 * t617 + 0.4e1 * t573 * t618;
t307 = t491 * t1433 + (t671 * t945 + t673 * t948) * t1118;
t290 = 0.4e1 * t525 * t587 + 0.4e1 * t526 * t586;
t279 = -0.4e1 * t1159 * t525 + 0.4e1 * t526 * t577;
t265 = 0.4e1 * t525 * t570 + 0.4e1 * t526 * t569;
t259 = t394 * t1415 + (t593 * t945 + t595 * t948) * t1054;
t258 = t1339 / 0.4e1;
t257 = t413 * t1433 + (t629 * t945 + t631 * t948) * t1118;
t254 = t1340 / 0.4e1;
t253 = t1341 / 0.4e1;
t252 = 0.4e1 * t979;
t247 = t1342 / 0.4e1;
t239 = (t579 - t774 + (t828 + t1263) * t948 + t1151) * t948 + t1150 * t945;
t238 = (t1048 * t948 - t1150 + t581) * t948 + (t1048 * t945 + t1050 + t580) * t945;
t237 = t1343 / 0.4e1;
t235 = t358 * t1415 + (t582 * t945 + t584 * t948) * t1054;
t234 = 0.4e1 * t1037;
t207 = t1344 / 0.4e1;
t185 = 0.4e1 * t394 * t471 + 0.4e1 * t593 * t675 + 0.4e1 * t595 * t677;
t167 = t1347 / 0.4e1;
t165 = t1348 / 0.4e1;
t162 = 0.4e1 * t980;
t159 = 0.4e1 * t358 * t471 + 0.4e1 * t582 * t675 + 0.4e1 * t584 * t677;
t157 = t1351 / 0.4e1;
t155 = 0.4e1 * t1038;
t153 = t1352 / 0.4e1;
t147 = t1353 / 0.4e1;
t133 = t1354 / 0.4e1;
t130 = t1355 / 0.4e1;
t128 = t1356 / 0.4e1;
t127 = t279 * t1400 + t378 * t1402 - t1566;
t122 = t1357 / 0.4e1;
t118 = t1358 / 0.4e1;
t91 = t958 + t1369 / 0.4e1 + t364 * t1402 + t290 * t1400;
t56 = t958 + (t899 / 0.2e1 + t898 / 0.2e1) * t947 + t265 * t1400 + t362 * t1402 - t1483;
t53 = t254 + t207 - t1339 / 0.4e1;
t52 = t258 + t254 - t1344 / 0.4e1;
t51 = t258 + t207 - t1340 / 0.4e1;
t49 = -t1390 / 0.4e1;
t48 = t1390 / 0.4e1;
t47 = t247 + t165 - t1341 / 0.4e1;
t46 = t253 + t247 - t1348 / 0.4e1;
t45 = t253 + t165 - t1342 / 0.4e1;
t38 = t237 + t130 - t1347 / 0.4e1;
t37 = t167 + t237 - t1355 / 0.4e1;
t36 = t167 + t130 - t1343 / 0.4e1;
t31 = t133 + t153 - t1356 / 0.4e1;
t30 = t128 + t153 - t1354 / 0.4e1;
t29 = t128 + t133 - t1352 / 0.4e1;
t28 = t157 + t48 - t1357 / 0.4e1;
t27 = t122 + t157 + t49;
t26 = t122 + t48 - t1351 / 0.4e1;
t25 = t1400 * t185 + t1402 * t307 + t1047;
t24 = t147 + t48 - t1358 / 0.4e1;
t23 = t118 + t147 + t49;
t22 = t118 + t48 - t1353 / 0.4e1;
t21 = t1400 * t159 + t1402 * t257 + t1047;
t20 = t1400 * t162 + t1402 * t252 + t1404 * t482 + t996;
t19 = t20 * qJD(3);
t18 = t1400 * t155 + t1402 * t234 + t1404 * t381 + t996;
t15 = (t486 / 0.2e1 - t239 / 0.2e1) * t948 + (t238 / 0.2e1 + t485 / 0.2e1) * t945 + t960;
t14 = t999 + t1484;
t13 = t999 - t1484;
t12 = (-t93 / 0.2e1 - t92 / 0.2e1 + (t1554 / 0.2e1 + t1555 / 0.2e1) * t948 + (-t81 / 0.2e1 - t82 / 0.2e1) * t945) * t929 + (t251 / 0.2e1 + t250 / 0.2e1 + (t67 / 0.2e1 + t68 / 0.2e1) * t948 + (t65 / 0.2e1 + t66 / 0.2e1) * t945) * t928 + t1455;
t11 = t12 * qJD(4);
t10 = t960 - t1485;
t9 = t960 + t1485;
t8 = t953 + t1104 + t1105;
t7 = t956 + t1334 + t1335 - t1454;
t6 = t950 + t1471 + t1463;
t5 = t951 - t1576;
t4 = t949 + t1576;
t3 = t950 + t1472 + t1464;
t2 = t951 - t1575;
t1 = t949 + t1575;
t16 = [t56 * qJD(2) + t91 * qJD(3) + t127 * qJD(4) + t363 * t1106, t56 * qJD(1) + t8 * qJD(3) + t3 * qJD(4) + t46 * qJD(5) + (t525 * t585 + t526 * t583 - t569 * t582 - t570 * t584) * t1327 + (t572 * t632 + t573 * t630 - t617 * t631 - t618 * t629) * t1330 + (t706 * t798 + t707 * t796) * t1333 + (t953 + (t1131 * t947 + t1133 * t944) * t1379 + t239 * t1376 + (t238 + t485) * t1380 + (t1132 * t947 - t1134 * t944 + t486) * t1377 + (t941 / 0.2e1 + t940 / 0.2e1) * t1023 + ((-t759 * t948 - t760 * t945) * t903 + (-t879 * t948 + t880 * t945) * t901) * m(3)) * qJD(2), t91 * qJD(1) + t8 * qJD(2) + t953 * qJD(3) + t6 * qJD(4) + t52 * qJD(5) + (-t586 * t593 - t587 * t595 + t450 / 0.2e1 + t451 / 0.2e1) * t1326 + (-t642 * t673 - t643 * t671 + t513 / 0.2e1 + t514 / 0.2e1) * t1329 + (t625 / 0.2e1 + t626 / 0.2e1 + (-t865 * t948 + t867 * t945) * t886) * t1332, t127 * qJD(1) + t3 * qJD(2) + t6 * qJD(3) + t37 * qJD(5) + (-t1159 * t1542 - t489 * t577 + t497 * t525 + t498 * t526) * t1325 + ((t376 / 0.2e1 + t374 / 0.2e1 + t448 / 0.2e1 + t447 / 0.2e1) * t948 + (t445 / 0.2e1 + t375 / 0.2e1 + t373 / 0.2e1 + t446 / 0.2e1) * t945) * t1120 + ((-t494 - t493) * t929 + (-t558 * t699 - t560 * t704 + t589 * t572 + t590 * t573) * m(5)) * qJD(4), t46 * qJD(2) + t52 * qJD(3) + t37 * qJD(4) + t1577 / 0.4e1; (t952 - (t899 + t898) * t947 / 0.2e1 + t1483) * qJD(1) + t15 * qJD(2) + t10 * qJD(3) + t2 * qJD(4) + t45 * qJD(5) + t265 * t1578 + t362 * t1540, t15 * qJD(1) + ((t358 * t452 - t582 * t583 - t584 * t585) * m(6) + (t413 * t509 - t629 * t630 - t631 * t632) * m(5) + m(4) * (-t1465 * t798 - t1466 * t796 + t538 * t616) + 0.4e1 * ((t945 * (rSges(3,1) * t1219 - t1122) + t948 * (rSges(3,1) * t1217 + t945 * rSges(3,3) - t920)) * (-t945 * t879 - t880 * t948) + t1121 * t903 * t901) * t1405 + (t940 * t874 + (-t945 * t873 + t1449) * t948) * t1379 + (t941 * t873 + (-t948 * t874 + t1449) * t945) * t1377 + t996) * qJD(2) + t18 * qJD(3) + t21 * qJD(4) + t235 * t1106, t10 * qJD(1) + t18 * qJD(2) + t996 * qJD(3) + t13 * qJD(4) + t29 * qJD(5) + (-t162 / 0.4e1 + t980 + t1038) * t1326 + (-t252 / 0.4e1 + t979 + t1037) * t1329 + (-t482 / 0.4e1 + t1036 + t1457) * t1332, t2 * qJD(1) + t21 * qJD(2) + t13 * qJD(3) + (t954 + (t1423 * t589 + t1424 * t590 + t1435 * t413 + t1099) * t1541 + (t1430 * t497 + t1431 * t498 + t1440 * t358 + t1102) * t1401) * qJD(4) + t23 * qJD(5), t45 * qJD(1) + t235 * t1327 / 0.4e1 + t29 * qJD(3) + t23 * qJD(4) + t453; (-t1369 / 0.4e1 + t952) * qJD(1) + t9 * qJD(2) + t960 * qJD(3) + t5 * qJD(4) + t51 * qJD(5) + t290 * t1578 + t364 * t1540, t9 * qJD(1) + t996 * qJD(2) + t19 + t14 * qJD(4) + t30 * qJD(5) + (t394 * t452 - t583 * t593 - t585 * t595 - t155 / 0.4e1 + t1038) * t1327 + (t491 * t509 - t630 * t671 - t632 * t673 - t234 / 0.4e1 + t1037) * t1330 + (-t381 / 0.4e1 + t627 * t616 + (-t796 * t945 - t798 * t948) * t886 + t1036) * t1333, qJD(1) * t960 + t20 * qJD(2) + t25 * qJD(4) + t1106 * t259 + t19, t5 * qJD(1) + t14 * qJD(2) + t25 * qJD(3) + (t954 + (t1418 * t589 + t1419 * t590 + t1435 * t491 + t1099) * t1541 + (t1426 * t497 + t1427 * t498 + t1440 * t394 + t1102) * t1401) * qJD(4) + t27 * qJD(5), t51 * qJD(1) + t30 * qJD(2) + t259 * t1326 / 0.4e1 + t27 * qJD(4) + t453; t1566 * qJD(1) + t1 * qJD(2) + t4 * qJD(3) + t36 * qJD(5) + t378 * t1540 + t1578 * t279, t1 * qJD(1) + t956 * qJD(2) + t7 * qJD(3) + t11 + t22 * qJD(5) + (-t159 / 0.4e1 + t406 * t452 + t1542 * t585 - t489 * t583 + t182 / 0.2e1 + t308 / 0.2e1 + t309 / 0.2e1) * t1327 + (-t257 / 0.4e1 + t531 * t509 + t558 * t632 - t560 * t630 + t266 / 0.2e1 + t403 / 0.2e1 + t404 / 0.2e1) * t1330, t4 * qJD(1) + t7 * qJD(2) + (t956 + (-t307 / 0.4e1 + t95 / 0.2e1) * m(5) + (-t185 / 0.4e1 + t35 / 0.2e1) * m(6)) * qJD(3) + t11 + t26 * qJD(5), (qJD(2) + qJD(3)) * t12 + (((t374 + t376) * t948 + (t373 + t375) * t945) * t1382 + t1501 * t1379 + t1500 * t1376) * t1120 + t1385 * t1324 + ((t531 * t550 + t558 * t589 - t560 * t590) * m(5) + (t494 / 0.2e1 + t493 / 0.2e1) * t929 ^ 2 + (t1542 * t497 + t406 * t449 - t489 * t498) * m(6)) * qJD(4), t36 * qJD(1) + t22 * qJD(2) + t26 * qJD(3) + t1325 * t1385 + (t633 * t1250 + t871 * t728 + t869 * t729 - t568 / 0.4e1) * t1324; -t1577 / 0.4e1 + t47 * qJD(2) + t53 * qJD(3) + t38 * qJD(4), t47 * qJD(1) + (t452 * t1250 + t869 * t583 + t871 * t585 + t347 / 0.2e1 + t565 / 0.2e1 - t566 / 0.2e1 - t235 / 0.4e1) * t1327 + t31 * qJD(3) + t24 * qJD(4) + t588, t53 * qJD(1) + t31 * qJD(2) + (t154 / 0.2e1 - t259 / 0.4e1) * t1326 + t28 * qJD(4) + t588, t38 * qJD(1) + t24 * qJD(2) + t28 * qJD(3) + (t872 * t1542 - t870 * t489 + t871 * t497 + t869 * t498 - t181 / 0.4e1 + (t406 * t946 + t449 * t943) * t928) * t1325 + t568 * t1106, (t568 * qJD(4) / 0.4e1 + (qJD(2) / 0.4e1 + qJD(3) / 0.4e1) * t592) * m(6);];
Cq = t16;
