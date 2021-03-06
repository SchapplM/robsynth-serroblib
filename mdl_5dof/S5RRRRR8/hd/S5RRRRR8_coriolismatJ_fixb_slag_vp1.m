% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:26
% EndTime: 2019-12-31 22:25:45
% DurationCPUTime: 65.00s
% Computational Cost: add. (285882->1281), mult. (288055->1763), div. (0->0), fcn. (315805->10), ass. (0->855)
t1021 = qJ(2) + qJ(3);
t1005 = sin(t1021);
t1027 = cos(qJ(1));
t1243 = t1005 * t1027;
t1024 = sin(qJ(1));
t1245 = t1005 * t1024;
t1174 = -t1245 / 0.2e1;
t1007 = cos(t1021);
t1020 = qJ(4) + qJ(5);
t1004 = sin(t1020);
t1006 = cos(t1020);
t1242 = t1007 * t1024;
t903 = t1004 * t1242 + t1006 * t1027;
t1232 = t1027 * t1004;
t904 = t1006 * t1242 - t1232;
t688 = Icges(6,5) * t904 - Icges(6,6) * t903 + Icges(6,3) * t1245;
t885 = Icges(6,4) * t904;
t692 = Icges(6,2) * t903 - Icges(6,6) * t1245 - t885;
t884 = Icges(6,4) * t903;
t694 = Icges(6,1) * t904 + Icges(6,5) * t1245 - t884;
t905 = t1006 * t1024 - t1007 * t1232;
t1241 = t1007 * t1027;
t906 = t1004 * t1024 + t1006 * t1241;
t467 = t1243 * t688 - t692 * t905 + t694 * t906;
t690 = Icges(6,5) * t906 + Icges(6,6) * t905 + Icges(6,3) * t1243;
t1436 = Icges(6,4) * t906;
t693 = Icges(6,2) * t905 + Icges(6,6) * t1243 + t1436;
t886 = Icges(6,4) * t905;
t696 = Icges(6,1) * t906 + Icges(6,5) * t1243 + t886;
t468 = t1243 * t690 + t693 * t905 + t696 * t906;
t1097 = t1024 * t467 + t1027 * t468;
t1108 = Icges(6,5) * t1006 - Icges(6,6) * t1004;
t837 = -Icges(6,3) * t1007 + t1005 * t1108;
t1343 = Icges(6,4) * t1006;
t1113 = -Icges(6,2) * t1004 + t1343;
t839 = -Icges(6,6) * t1007 + t1005 * t1113;
t1344 = Icges(6,4) * t1004;
t1115 = Icges(6,1) * t1006 - t1344;
t841 = -Icges(6,5) * t1007 + t1005 * t1115;
t563 = t1243 * t837 + t839 * t905 + t841 * t906;
t1667 = t1005 * t1097 - t1007 * t563;
t1675 = t1667 * t1174;
t1171 = t1245 / 0.4e1;
t1173 = -t1245 / 0.4e1;
t1674 = t1171 + t1173;
t1025 = cos(qJ(4));
t1234 = t1025 * t1027;
t1022 = sin(qJ(4));
t1239 = t1022 * t1024;
t942 = t1007 * t1239 + t1234;
t1231 = t1027 * t1022;
t1236 = t1024 * t1025;
t943 = t1007 * t1236 - t1231;
t724 = Icges(5,5) * t943 - Icges(5,6) * t942 + Icges(5,3) * t1245;
t923 = Icges(5,4) * t943;
t728 = Icges(5,2) * t942 - Icges(5,6) * t1245 - t923;
t922 = Icges(5,4) * t942;
t730 = Icges(5,1) * t943 + Icges(5,5) * t1245 - t922;
t508 = t1245 * t724 + t728 * t942 + t730 * t943;
t944 = -t1007 * t1231 + t1236;
t945 = t1007 * t1234 + t1239;
t726 = Icges(5,5) * t945 + Icges(5,6) * t944 + Icges(5,3) * t1243;
t1437 = Icges(5,4) * t945;
t729 = Icges(5,2) * t944 + Icges(5,6) * t1243 + t1437;
t924 = Icges(5,4) * t944;
t732 = Icges(5,1) * t945 + Icges(5,5) * t1243 + t924;
t509 = t1245 * t726 - t729 * t942 + t732 * t943;
t1095 = t1024 * t508 + t1027 * t509;
t1109 = Icges(5,5) * t1025 - Icges(5,6) * t1022;
t862 = -Icges(5,3) * t1007 + t1005 * t1109;
t1345 = Icges(5,4) * t1025;
t1114 = -Icges(5,2) * t1022 + t1345;
t864 = -Icges(5,6) * t1007 + t1005 * t1114;
t1346 = Icges(5,4) * t1022;
t1116 = Icges(5,1) * t1025 - t1346;
t866 = -Icges(5,5) * t1007 + t1005 * t1116;
t584 = t1245 * t862 - t864 * t942 + t866 * t943;
t130 = -t1005 * t1095 + t1007 * t584;
t465 = t1245 * t688 + t692 * t903 + t694 * t904;
t466 = t1245 * t690 - t693 * t903 + t696 * t904;
t1098 = t1024 * t465 + t1027 * t466;
t561 = t1245 * t837 - t839 * t903 + t841 * t904;
t239 = t1005 * t1098 - t1007 * t561;
t1454 = m(6) * qJD(1);
t1658 = -t1454 / 0.4e1;
t1657 = -m(5) * qJD(1) / 0.4e1;
t1601 = -t904 * rSges(6,1) + t903 * rSges(6,2);
t697 = rSges(6,3) * t1245 - t1601;
t999 = pkin(4) * t1025 + pkin(3);
t1136 = pkin(4) * t1231 - t1242 * t999;
t1448 = t1007 * pkin(3);
t1028 = -pkin(9) - pkin(8);
t1447 = pkin(8) + t1028;
t1595 = t1005 * t1447;
t722 = (t1448 + t1595) * t1024 + t1136;
t1631 = t697 - t722;
t1317 = t1007 * t722;
t1441 = rSges(6,1) * t1006;
t1119 = -rSges(6,2) * t1004 + t1441;
t1438 = rSges(6,3) * t1007;
t843 = t1005 * t1119 - t1438;
t795 = t843 * t1245;
t616 = t1007 * t697 + t795;
t1460 = pkin(3) - t999;
t1594 = t1005 * t1460;
t829 = t1007 * t1447 - t1594;
t774 = t829 * t1245;
t506 = t774 - t1317 + t616;
t1408 = t1007 * t1631 - t506 + t774 + t795;
t1143 = 0.2e1 * t1408;
t1518 = m(6) / 0.4e1;
t1669 = t1143 * t1518;
t510 = t1243 * t724 - t728 * t944 + t730 * t945;
t511 = t1243 * t726 + t729 * t944 + t732 * t945;
t1094 = t1024 * t510 + t1027 * t511;
t586 = t1243 * t862 + t864 * t944 + t866 * t945;
t1668 = t1005 * t1094 - t1007 * t586;
t1659 = t1024 * t468 - t1027 * t467;
t1660 = t1024 * t511 - t1027 * t510;
t1666 = t1659 + t1660;
t927 = (-Icges(5,5) * t1022 - Icges(5,6) * t1025) * t1005;
t1309 = t1007 * t927;
t930 = (-Icges(5,2) * t1025 - t1346) * t1005;
t933 = (-Icges(5,1) * t1022 - t1345) * t1005;
t1665 = t1309 / 0.2e1 - (-(t866 / 0.2e1 + t930 / 0.2e1) * t1022 + (t933 / 0.2e1 - t864 / 0.2e1) * t1025) * t1005;
t1664 = t1674 * t1659;
t1663 = t1674 * t1660;
t1600 = -t943 * rSges(5,1) + t942 * rSges(5,2);
t733 = rSges(5,3) * t1245 - t1600;
t1442 = rSges(5,1) * t1025;
t1120 = -rSges(5,2) * t1022 + t1442;
t1440 = rSges(5,3) * t1007;
t868 = t1005 * t1120 - t1440;
t629 = t1007 * t733 + t1245 * t868;
t1520 = m(5) / 0.4e1;
t1647 = -t1005 / 0.2e1;
t1646 = t1007 / 0.2e1;
t1645 = t1024 / 0.4e1;
t1644 = -t1027 / 0.4e1;
t1165 = -t1242 / 0.4e1;
t1102 = t1022 * t728 + t1025 * t730;
t1075 = -t1024 * t862 - t1102;
t818 = t864 * t1024;
t820 = t866 * t1024;
t406 = -t1075 * t1007 + (t1022 * t818 - t1025 * t820 + t724) * t1005;
t1287 = t1025 * t866;
t1306 = t1022 * t864;
t1100 = t1287 - t1306;
t863 = Icges(5,3) * t1005 + t1007 * t1109;
t1073 = -t1100 + t863;
t1312 = t1007 * t862;
t1046 = t1005 * t1073 + t1312;
t865 = Icges(5,6) * t1005 + t1007 * t1114;
t867 = Icges(5,5) * t1005 + t1007 * t1116;
t461 = t1024 * t1046 - t865 * t942 + t867 * t943;
t1655 = t406 + t461;
t1101 = -t1022 * t729 + t1025 * t732;
t1074 = -t1027 * t862 - t1101;
t819 = t864 * t1027;
t821 = t866 * t1027;
t407 = -t1074 * t1007 + (t1022 * t819 - t1025 * t821 + t726) * t1005;
t462 = t1027 * t1046 + t865 * t944 + t867 * t945;
t1654 = t407 + t462;
t1316 = t1007 * t724;
t537 = t1005 * t1102 - t1316;
t1653 = t537 + t584;
t1315 = t1007 * t726;
t539 = t1005 * t1101 - t1315;
t1652 = t539 + t586;
t1240 = t1007 * t1028;
t1247 = t1005 * t1006;
t1199 = rSges(6,1) * t1247;
t1248 = t1004 * t1005;
t1197 = rSges(6,2) * t1248;
t1353 = rSges(6,3) * t1241 + t1027 * t1197;
t794 = -t1027 * t1199 + t1353;
t989 = pkin(8) * t1241;
t1651 = -t989 + (-t1240 + t1594) * t1027 + t794;
t997 = Icges(4,4) * t1007;
t956 = -Icges(4,2) * t1005 + t997;
t957 = Icges(4,1) * t1005 + t997;
t1650 = -t956 - t957;
t1018 = t1024 ^ 2;
t1019 = t1027 ^ 2;
t1229 = t1018 + t1019;
t1026 = cos(qJ(2));
t1450 = pkin(2) * t1026;
t1000 = pkin(1) + t1450;
t1517 = -pkin(7) - pkin(6);
t1001 = t1024 * t1517;
t1125 = rSges(5,1) * t945 + rSges(5,2) * t944;
t1501 = rSges(5,3) + pkin(8);
t1580 = -t1005 * t1501 - t1448;
t641 = -t1001 + (t1000 - t1580) * t1027 + t1125;
t1544 = 0.2e1 * t641;
t1246 = t1005 * t1022;
t1198 = rSges(5,2) * t1246;
t1244 = t1005 * t1025;
t1200 = rSges(5,1) * t1244;
t967 = t1024 * t1200;
t822 = -t967 + (t1198 + t1440) * t1024;
t869 = rSges(5,3) * t1005 + t1007 * t1120;
t542 = (t1024 * t868 + t822) * t1007 + (t1024 * t869 - t733) * t1005;
t735 = rSges(5,3) * t1243 + t1125;
t1351 = rSges(5,3) * t1241 + t1027 * t1198;
t823 = -t1027 * t1200 + t1351;
t543 = (-t1027 * t868 - t823) * t1007 + (-t1027 * t869 + t735) * t1005;
t1250 = -t1000 * t1024 - t1027 * t1517;
t640 = t1024 * t1580 + t1250 + t1600;
t1410 = t1544 * t543 + 0.2e1 * t542 * t640;
t1123 = rSges(6,1) * t906 + rSges(6,2) * t905;
t1220 = pkin(4) * t1239;
t1137 = -t1028 * t1243 + t1220;
t1439 = rSges(6,3) * t1005;
t611 = -t1001 + (t1007 * t999 + t1000 + t1439) * t1027 + t1123 + t1137;
t1558 = 0.2e1 * t611;
t610 = (-rSges(6,3) + t1028) * t1245 + t1136 + t1250 + t1601;
t1559 = 0.2e1 * t610;
t1084 = t1197 + t1438;
t1591 = t1024 * t1084;
t961 = t1024 * t1199;
t793 = -t961 + t1591;
t844 = t1007 * t1119 + t1439;
t1205 = t1007 * t793 + t1242 * t843 + t1245 * t844;
t988 = pkin(3) * t1245;
t1135 = pkin(8) * t1242 - t988;
t1354 = t1024 * t1240 + t1245 * t999;
t755 = -t1135 - t1354;
t1182 = t1460 * t1007;
t830 = -t1182 - t1595;
t352 = (t1024 * t829 + t755) * t1007 + (t1024 * t830 - t1631) * t1005 + t1205;
t1380 = -t830 - t844;
t1381 = t829 + t843;
t699 = rSges(6,3) * t1243 + t1123;
t676 = t1005 * t699;
t1449 = t1005 * pkin(8);
t723 = (-t1182 - t1449) * t1027 + t1137;
t353 = t676 + (t1027 * t1380 + t723) * t1005 + (-t1027 * t1381 - t1651) * t1007;
t1420 = t1558 * t353 + t1559 * t352;
t1122 = (-pkin(3) - t1442) * t1005;
t1201 = t989 + t1351;
t719 = t1027 * t1122 + t1201;
t1534 = -0.2e1 * t719;
t1071 = -t1007 * t1501 - t1198;
t1352 = t967 + t988;
t718 = t1024 * t1071 + t1352;
t1535 = 0.2e1 * t718;
t1065 = -t1240 + (-t999 - t1441) * t1005;
t675 = t1027 * t1065 + t1353;
t1538 = -0.2e1 * t675;
t1202 = t961 + t1354;
t674 = t1202 - t1591;
t1539 = 0.2e1 * t674;
t1389 = t699 + t723;
t1593 = t1007 * t1389;
t507 = t1243 * t1381 + t1593;
t631 = t1007 * t735 + t1243 * t868;
t1605 = (t1538 * t507 + t1539 * t506 + t1420) * t1518 + (t1534 * t631 + t1535 * t629 + t1410) * t1520;
t965 = pkin(3) * t1005 - pkin(8) * t1007;
t1204 = t965 + t1381;
t633 = t1204 * t1024;
t1607 = t633 * t1669;
t1649 = t1605 + t1607;
t1023 = sin(qJ(2));
t1451 = pkin(2) * t1023;
t655 = (t1065 - t1451) * t1027 + t1353;
t1542 = -0.2e1 * t655;
t654 = (-t1084 + t1451) * t1024 + t1202;
t1543 = 0.2e1 * t654;
t1552 = -0.2e1 * t631;
t683 = (t1071 + t1451) * t1024 + t1352;
t684 = (t1122 - t1451) * t1027 + t1201;
t1606 = (t1552 * t684 + 0.2e1 * t629 * t683 + t1410) * t1520 + (t1542 * t507 + t1543 * t506 + t1420) * t1518;
t1158 = t965 + t1451;
t1091 = t1158 + t1381;
t620 = t1091 * t1024;
t1608 = t620 * t1669;
t1648 = t1606 + t1608;
t1107 = t1004 * t692 + t1006 * t694;
t1320 = t1007 * t688;
t494 = t1005 * t1107 - t1320;
t1522 = m(4) / 0.4e1;
t1521 = m(5) / 0.2e1;
t1467 = t1005 / 0.2e1;
t1466 = -t1007 / 0.2e1;
t1464 = t1024 / 0.2e1;
t1462 = -t1027 / 0.2e1;
t1461 = t1027 / 0.2e1;
t1078 = -t1024 * t837 - t1107;
t789 = t839 * t1024;
t791 = t841 * t1024;
t374 = -t1078 * t1007 + (t1004 * t789 - t1006 * t791 + t688) * t1005;
t1327 = t1006 * t841;
t1334 = t1004 * t839;
t1105 = t1327 - t1334;
t838 = Icges(6,3) * t1005 + t1007 * t1108;
t1076 = -t1105 + t838;
t1314 = t1007 * t837;
t1047 = t1005 * t1076 + t1314;
t840 = Icges(6,6) * t1005 + t1007 * t1113;
t842 = Icges(6,5) * t1005 + t1007 * t1115;
t427 = t1024 * t1047 - t840 * t903 + t842 * t904;
t1637 = t374 + t427;
t1106 = -t1004 * t693 + t1006 * t696;
t1077 = -t1027 * t837 - t1106;
t790 = t839 * t1027;
t792 = t841 * t1027;
t375 = -t1077 * t1007 + (t1004 * t790 - t1006 * t792 + t690) * t1005;
t428 = t1027 * t1047 + t840 * t905 + t842 * t906;
t1636 = t375 + t428;
t1634 = t494 + t561;
t1319 = t1007 * t690;
t496 = t1005 * t1106 - t1319;
t1633 = t496 + t563;
t1376 = t841 + (-Icges(6,2) * t1006 - t1344) * t1005;
t1443 = rSges(4,1) * t1007;
t964 = -rSges(4,2) * t1005 + t1443;
t1288 = t1024 * t964;
t1155 = -rSges(4,2) * t1243 + rSges(4,3) * t1024;
t784 = -t1001 + (t1000 + t1443) * t1027 + t1155;
t700 = -0.2e1 * t784 * t1288;
t1524 = -0.2e1 * t1027;
t1221 = t964 * t1524;
t881 = rSges(4,1) * t1242 - rSges(4,2) * t1245 - rSges(4,3) * t1027;
t783 = -t881 + t1250;
t701 = t783 * t1221;
t1388 = t700 + t701;
t966 = t1448 + t1449;
t1368 = -t869 - t966;
t750 = t1368 * t1024;
t578 = t750 * t1544;
t752 = t1368 * t1027;
t1529 = 0.2e1 * t752;
t579 = t640 * t1529;
t1402 = t578 + t579;
t1203 = -t966 + t1380;
t634 = t1203 * t1024;
t516 = t634 * t1558;
t636 = t1203 * t1027;
t1547 = 0.2e1 * t636;
t517 = t610 * t1547;
t1407 = t516 + t517;
t1369 = t868 + t965;
t751 = t1369 * t1027;
t1530 = -0.2e1 * t751;
t749 = t1369 * t1024;
t1531 = -0.2e1 * t749;
t963 = rSges(4,1) * t1005 + rSges(4,2) * t1007;
t1070 = t963 + t1451;
t1589 = t1070 * t1027;
t1590 = t1070 * t1024;
t635 = t1204 * t1027;
t1216 = (t1542 * t633 - t1543 * t635 + t1407) * t1518 + (t1530 * t683 + t1531 * t684 + t1402) * t1520 + (0.2e1 * (t1024 * t1589 - t1027 * t1590) * t963 + t1388) * t1522;
t622 = t1091 * t1027;
t1140 = t1158 + t868;
t709 = t1140 * t1024;
t711 = t1140 * t1027;
t939 = t963 * t1024;
t940 = t963 * t1027;
t1217 = (t1538 * t620 - t1539 * t622 + t1407) * t1518 + (t1534 * t709 - t1535 * t711 + t1402) * t1520 + (-0.2e1 * t1589 * t939 + 0.2e1 * t1590 * t940 + t1388) * t1522;
t1615 = t1216 - t1217;
t1092 = -t1024 * t735 + t1027 * t733;
t595 = t1092 * t1005;
t1561 = 0.2e1 * t595;
t1362 = t1024 * t1135 + t1027 * (-pkin(3) * t1243 + t989);
t1602 = t1024 * t822 + t1027 * t823;
t591 = t1362 + t1602;
t1208 = t1529 * t629 + t1552 * t750 + t1561 * t591;
t493 = t1092 * t1007 + (-t1024 * t823 + t1027 * t822) * t1005;
t1566 = 0.2e1 * t493;
t1361 = t1229 * t966;
t553 = t1024 * t733 + t1027 * t735 + t1361;
t140 = t1530 * t542 + t1531 * t543 + t1566 * t553 + t1208;
t1549 = -0.2e1 * t634;
t669 = t697 * t1243;
t436 = t669 + (-t1024 * t1389 - t1027 * t722) * t1005;
t1568 = 0.2e1 * t436;
t1581 = t1651 * t1027 + (t755 + t793) * t1024;
t472 = t1362 + t1581;
t1215 = t1547 * t506 + t1549 * t507 + t1568 * t472;
t1548 = -0.2e1 * t635;
t1550 = -0.2e1 * t633;
t1394 = t1241 * t697 + t1243 * t793;
t314 = (t1005 * t755 - t1317) * t1027 + (-t1005 * t1651 - t1593) * t1024 + t1394;
t1570 = 0.2e1 * t314;
t422 = t1024 * t1631 + t1027 * t1389 + t1361;
t53 = t1548 * t352 + t1550 * t353 + t1570 * t422 + t1215;
t1445 = t140 * t1520 + t1518 * t53;
t1017 = t1027 * pkin(6);
t1374 = -t1024 * (pkin(1) * t1024 - t1017 + t1250) + t1027 * (-pkin(6) * t1024 - t1001 + (-pkin(1) + t1000) * t1027);
t377 = t422 + t1374;
t190 = t377 * t1570;
t1553 = -0.2e1 * t622;
t303 = t352 * t1553;
t1554 = -0.2e1 * t620;
t304 = t353 * t1554;
t483 = t553 + t1374;
t323 = t483 * t1566;
t1536 = -0.2e1 * t711;
t473 = t542 * t1536;
t1537 = -0.2e1 * t709;
t474 = t543 * t1537;
t1446 = (t473 + t474 + t323 + t1208) * t1520 + (t303 + t304 + t190 + t1215) * t1518;
t1614 = t1445 - t1446;
t491 = -t1073 * t1007 + (-t1022 * t865 + t1025 * t867 + t862) * t1005;
t615 = t1005 * t1100 - t1312;
t1613 = t1646 * t491 + t1647 * t615;
t938 = (-rSges(5,1) * t1022 - rSges(5,2) * t1025) * t1005;
t1222 = t938 * t1524;
t1289 = t1024 * t938;
t1401 = t1222 * t640 - 0.2e1 * t1289 * t641;
t894 = (-rSges(6,1) * t1004 - rSges(6,2) * t1006) * t1005;
t1138 = pkin(4) * t1246 - t894;
t786 = t1138 * t1027;
t1525 = 0.2e1 * t786;
t785 = t1138 * t1024;
t1404 = t1525 * t610 + t1558 * t785;
t782 = rSges(5,1) * t944 - rSges(5,2) * t945;
t1527 = -0.2e1 * t782;
t781 = -rSges(5,1) * t942 - rSges(5,2) * t943;
t1528 = -0.2e1 * t781;
t743 = rSges(6,1) * t905 - rSges(6,2) * t906;
t926 = t944 * pkin(4);
t1383 = -t743 - t926;
t1540 = 0.2e1 * t1383;
t742 = -rSges(6,1) * t903 - rSges(6,2) * t904;
t925 = t942 * pkin(4);
t664 = -t742 + t925;
t1541 = 0.2e1 * t664;
t1603 = (t1540 * t633 - t1541 * t635 + t1404) * t1518 + (t1527 * t749 - t1528 * t751 + t1401) * t1520;
t1604 = (t1540 * t620 - t1541 * t622 + t1404) * t1518 + (t1527 * t709 - t1528 * t711 + t1401) * t1520;
t1384 = -Icges(5,2) * t945 + t732 + t924;
t1386 = Icges(5,1) * t944 - t1437 - t729;
t776 = Icges(5,5) * t944 - Icges(5,6) * t945;
t434 = -t1007 * t776 + (-t1022 * t1384 + t1025 * t1386) * t1005;
t1259 = t434 * t1024;
t1385 = -Icges(5,2) * t943 + t730 - t922;
t1387 = -Icges(5,1) * t942 + t728 - t923;
t775 = -Icges(5,5) * t942 - Icges(5,6) * t943;
t433 = -t1007 * t775 + (-t1022 * t1385 + t1025 * t1387) * t1005;
t1260 = t433 * t1027;
t1370 = t866 + t930;
t1371 = -t864 + t933;
t514 = t1245 * t927 - t1370 * t942 + t1371 * t943;
t515 = t1243 * t927 + t1370 * t944 + t1371 * t945;
t1582 = -t1260 / 0.4e1 + t1259 / 0.4e1 + t515 * t1645 + t514 * t1644;
t1279 = t1027 * t1668;
t1297 = t1024 * t130;
t1612 = t1279 / 0.4e1 + t1668 * t1644 - t1297 / 0.4e1 + t130 * t1645;
t1610 = 0.4e1 * m(6);
t1526 = -0.2e1 * t785;
t1014 = Icges(3,4) * t1026;
t973 = -Icges(3,2) * t1023 + t1014;
t974 = Icges(3,1) * t1023 + t1014;
t703 = t1024 * t881 + t1027 * (rSges(4,1) * t1241 + t1155);
t753 = -t1024 * t939 - t1027 * t940;
t1588 = t1229 * t963 * t964 + t703 * t753;
t460 = -t1076 * t1007 + (-t1004 * t840 + t1006 * t842 + t837) * t1005;
t593 = t1005 * t1105 - t1314;
t1587 = t1466 * t460 + t593 * t1467;
t398 = t1245 * t775 - t1385 * t942 + t1387 * t943;
t399 = t1245 * t776 - t1384 * t942 + t1386 * t943;
t247 = t1024 * t399 - t1027 * t398;
t400 = t1243 * t775 + t1385 * t944 + t1387 * t945;
t401 = t1243 * t776 + t1384 * t944 + t1386 * t945;
t248 = t1024 * t401 - t1027 * t400;
t1586 = t1462 * t247 + t1464 * t248;
t1391 = -Icges(6,2) * t904 + t694 - t884;
t1393 = -Icges(6,1) * t903 + t692 - t885;
t736 = -Icges(6,5) * t903 - Icges(6,6) * t904;
t356 = t1245 * t736 - t1391 * t903 + t1393 * t904;
t1390 = -Icges(6,2) * t906 + t696 + t886;
t1392 = Icges(6,1) * t905 - t1436 - t693;
t737 = Icges(6,5) * t905 - Icges(6,6) * t906;
t357 = t1245 * t737 - t1390 * t903 + t1392 * t904;
t205 = t1024 * t357 - t1027 * t356;
t358 = t1243 * t736 + t1391 * t905 + t1393 * t906;
t359 = t1243 * t737 + t1390 * t905 + t1392 * t906;
t206 = t1024 * t359 - t1027 * t358;
t1421 = t1462 * t205 + t1464 * t206;
t1256 = t494 * t1024;
t1096 = t1027 * t496 + t1256;
t1585 = ((t1096 - t460) * t1007 + (t1024 * t374 + t1027 * t375 + t593) * t1005) * t1466 + (t1005 * t1096 - t1007 * t593) * t1467;
t1144 = 0.2e1 * t377 + 0.2e1 * t422;
t1225 = 0.2e1 * t1027;
t1226 = 0.2e1 * t1024;
t1397 = -t622 - t635;
t1398 = t620 + t633;
t612 = t1024 * t742 + t1027 * t743;
t558 = -t1024 * t925 + t1027 * t926 + t612;
t632 = t1024 * t781 + t1027 * t782;
t1584 = (t558 * t1144 + t1397 * t1525 + t1398 * t1526) * t1518 + (0.2e1 * (t483 + t553) * t632 + ((t711 + t751) * t1225 + (t709 + t749) * t1226) * t938) * t1520;
t1583 = m(5) * (t493 * t595 + t542 * t629 - t543 * t631) + (t314 * t436 + t352 * t506 - t353 * t507) * t1610 / 0.4e1;
t1131 = t483 * t591 - t709 * t750 - t711 * t752;
t1133 = t377 * t472 - t620 * t634 - t622 * t636;
t702 = t742 * t1243;
t598 = -t1245 * t743 + t702;
t637 = t1007 * t742 + t1245 * t894;
t638 = -t1007 * t743 - t1243 * t894;
t1132 = t436 * t598 + t506 * t637 - t507 * t638;
t1523 = m(3) / 0.4e1;
t1444 = rSges(3,1) * t1026;
t1183 = pkin(1) + t1444;
t1238 = t1023 * t1024;
t1249 = rSges(3,2) * t1238 + rSges(3,3) * t1027;
t845 = -t1024 * t1183 + t1017 + t1249;
t1237 = t1023 * t1027;
t995 = rSges(3,2) * t1237;
t846 = -t995 + t1183 * t1027 + (rSges(3,3) + pkin(6)) * t1024;
t976 = rSges(3,1) * t1023 + rSges(3,2) * t1026;
t952 = t976 * t1024;
t953 = t976 * t1027;
t1348 = Icges(3,4) * t1023;
t972 = Icges(3,2) * t1026 + t1348;
t975 = Icges(3,1) * t1026 - t1348;
t1579 = -(t975 / 0.2e1 - t972 / 0.2e1) * t1023 - 0.4e1 * (-t1589 * t784 + t1590 * t783) * t1522 - 0.4e1 * (t845 * t952 - t846 * t953) * t1523;
t392 = -t1007 * t737 + (-t1004 * t1390 + t1006 * t1392) * t1005;
t1263 = t392 * t1024;
t391 = -t1007 * t736 + (-t1004 * t1391 + t1006 * t1393) * t1005;
t1264 = t391 * t1027;
t889 = (-Icges(6,1) * t1004 - t1343) * t1005;
t1377 = -t839 + t889;
t887 = (-Icges(6,5) * t1004 - Icges(6,6) * t1006) * t1005;
t463 = t1245 * t887 - t1376 * t903 + t1377 * t904;
t464 = t1243 * t887 + t1376 * t905 + t1377 * t906;
t1145 = -t1264 / 0.4e1 + t1263 / 0.4e1 + t464 * t1645 + t463 * t1644;
t1576 = (t1027 / 0.4e1 + t1644) * t1667;
t1575 = t1612 + t1663;
t1159 = t1241 / 0.4e1;
t1163 = t1242 / 0.4e1;
t1167 = t1243 / 0.4e1;
t1574 = t1159 * t1652 + t1163 * t1653 + t1167 * t1654 + t1171 * t1655 - t1613;
t920 = Icges(3,5) * t1024 + t1027 * t975;
t1356 = -t1027 * t972 + t920;
t1235 = t1024 * t1026;
t993 = Icges(3,4) * t1238;
t919 = Icges(3,1) * t1235 - Icges(3,5) * t1027 - t993;
t1357 = -Icges(3,2) * t1235 + t919 - t993;
t918 = Icges(3,6) * t1024 + t1027 * t973;
t1358 = -t1027 * t974 - t918;
t917 = Icges(3,4) * t1235 - Icges(3,2) * t1238 - Icges(3,6) * t1027;
t1359 = t1024 * t974 + t917;
t1573 = (-t1024 * t1356 + t1027 * t1357) * t1023 + (t1024 * t1358 + t1027 * t1359) * t1026;
t1347 = Icges(4,4) * t1005;
t958 = Icges(4,1) * t1007 - t1347;
t880 = Icges(4,5) * t1024 + t1027 * t958;
t955 = Icges(4,2) * t1007 + t1347;
t1364 = -t1027 * t955 + t880;
t982 = Icges(4,4) * t1245;
t879 = Icges(4,1) * t1242 - Icges(4,5) * t1027 - t982;
t1365 = -Icges(4,2) * t1242 + t879 - t982;
t878 = Icges(4,6) * t1024 + t1027 * t956;
t1366 = -t1027 * t957 - t878;
t877 = Icges(4,4) * t1242 - Icges(4,2) * t1245 - Icges(4,6) * t1027;
t1367 = t1024 * t957 + t877;
t1572 = (-t1024 * t1364 + t1027 * t1365) * t1005 + (t1024 * t1366 + t1027 * t1367) * t1007;
t831 = t880 * t1242;
t954 = Icges(4,5) * t1007 - Icges(4,6) * t1005;
t1267 = t1027 * t954;
t876 = Icges(4,3) * t1024 + t1267;
t1148 = t1027 * t876 - t831;
t875 = Icges(4,5) * t1242 - Icges(4,6) * t1245 - Icges(4,3) * t1027;
t1150 = t1005 * t878 - t875;
t1330 = t1005 * t877;
t1378 = t1024 * t876 + t1241 * t880;
t1379 = -t1024 * t875 - t1241 * t879;
t625 = -t1245 * t878 - t1148;
t626 = -t1243 * t877 - t1379;
t627 = -t1243 * t878 + t1378;
t1042 = ((t625 - t831 + (t876 + t1330) * t1027 + t1379) * t1027 + t1378 * t1024 + t1666) * t1462 + (t1024 * t627 - t1027 * t626 + t1666) * t1461 + ((t1024 * t1150 + t1148 + t625 + t626) * t1024 + (-t1378 + t627 + t1024 * (-t1007 * t879 + t1330) + (t1150 + t875) * t1027) * t1027) * t1464;
t1571 = t867 * t1244 / 0.2e1 - t865 * t1246 / 0.2e1 + t955 * t1647 + (t958 + t862 + t837) * t1467 + (t1287 + t1327 - t1650) * t1646 + (t1306 + t1334 + t863 + t838) * t1466;
t1318 = t1007 * t699;
t471 = (-t1005 * t794 - t1318) * t1024 + t1394;
t1567 = 0.2e1 * t471;
t540 = t702 + (t1024 * t1383 - t1027 * t925) * t1005;
t1565 = 0.2e1 * t540;
t577 = -t1245 * t699 + t669;
t1563 = 0.2e1 * t577;
t1560 = 0.2e1 * t598;
t1557 = 0.2e1 * t612;
t1556 = 0.4e1 * t612;
t619 = (-t1024 * t782 + t1027 * t781) * t1005;
t1555 = 0.2e1 * t619;
t1551 = 0.4e1 * t632;
t1546 = 0.2e1 * t637;
t1545 = -0.2e1 * t638;
t1533 = -0.2e1 * t742;
t1532 = -0.2e1 * t743;
t1519 = m(6) / 0.2e1;
t524 = -t1005 * t697 + t1205;
t525 = -t844 * t1243 + t676 + (-t1027 * t843 - t794) * t1007;
t618 = t1243 * t843 + t1318;
t1516 = 0.2e1 * (t314 * t577 + t352 * t616 - t353 * t618 + t436 * t471 + t506 * t524 - t507 * t525) * m(6);
t1511 = m(6) * t1143 * t618;
t1211 = t1547 * t616 + t1549 * t618 + t1563 * t472;
t261 = t377 * t1567;
t420 = t524 * t1553;
t421 = t525 * t1554;
t1509 = m(6) * (t420 + t421 + t261 + t1211);
t69 = t1548 * t524 + t1550 * t525 + t1567 * t422 + t1211;
t1507 = m(6) * t69;
t1223 = t894 * t1524;
t1290 = t1024 * t894;
t1224 = 0.2e1 * t1290;
t1210 = t1223 * t506 + t1224 * t507 + t1557 * t436;
t1214 = t1545 * t620 - t1546 * t622 + t1560 * t377;
t1505 = m(6) * (t1210 + t1214);
t1212 = t1545 * t633 - t1546 * t635 + t1560 * t422;
t1503 = m(6) * (t1210 + t1212);
t1497 = 0.4e1 * m(4) * (t783 * t939 - t784 * t940);
t1209 = t1525 * t616 + t1526 * t618 + t1563 * t558;
t1481 = m(6) * (t1209 + t1214);
t1480 = m(6) * (t1209 + t1212);
t1479 = m(6) * (t612 * t1144 + (-t1225 * t1397 + t1226 * t1398) * t894);
t1478 = (t471 * t577 + t524 * t616 - t525 * t618) * t1610;
t1411 = t1558 * t525 + t1559 * t524;
t1477 = m(6) * (t1542 * t618 + t1543 * t616 + t1411);
t1476 = m(6) * (t1538 * t618 + t1539 * t616 + t1411);
t518 = t610 * t1546;
t519 = t638 * t1558;
t1406 = t518 + t519;
t1475 = m(6) * (t1532 * t507 + t1533 * t506 + t1406);
t1474 = m(6) * (t1540 * t618 + t1541 * t616 + t1406);
t1403 = t1223 * t610 - 0.2e1 * t1290 * t611;
t1469 = m(6) * (t1532 * t620 - t1533 * t622 + t1403);
t1468 = m(6) * (t1532 * t633 - t1533 * t635 + t1403);
t1465 = -t1024 / 0.2e1;
t1459 = m(4) * qJD(2);
t1458 = m(4) * qJD(3);
t1456 = m(5) * qJD(2);
t1455 = m(5) * qJD(3);
t1453 = m(6) * qJD(2);
t1452 = m(6) * qJD(3);
t603 = t703 + t1374;
t1428 = t603 * t753;
t1233 = t1026 * t1027;
t915 = Icges(3,5) * t1235 - Icges(3,6) * t1238 - Icges(3,3) * t1027;
t1373 = -t1024 * t915 - t1233 * t919;
t1112 = Icges(3,5) * t1026 - Icges(3,6) * t1023;
t916 = Icges(3,3) * t1024 + t1027 * t1112;
t1372 = t1024 * t916 + t1233 * t920;
t1310 = t1007 * t887;
t1305 = t1023 * t917;
t1266 = t374 * t1027;
t1265 = t375 * t1024;
t1262 = t406 * t1027;
t1261 = t407 * t1024;
t532 = -t1310 + (-t1004 * t1376 + t1006 * t1377) * t1005;
t1255 = t532 * t1007;
t1254 = t537 * t1024;
t1230 = qJD(2) + qJD(3);
t1228 = 0.4e1 * t894;
t1227 = 0.4e1 * t938;
t1219 = t1479 / 0.4e1 + t1421;
t1168 = t1243 / 0.2e1;
t1172 = t1245 / 0.2e1;
t163 = -t1007 * t463 + (t1024 * t356 + t1027 * t357) * t1005;
t164 = -t1007 * t464 + (t1024 * t358 + t1027 * t359) * t1005;
t1218 = t164 * t1168 + t163 * t1172 + (-t1255 + (t1024 * t391 + t1027 * t392) * t1005) * t1466;
t1213 = t1525 * t506 + t1526 * t507 + t1568 * t558;
t1207 = t1223 * t616 + t1224 * t618 + t1557 * t577;
t1206 = t1222 * t629 + 0.2e1 * t1289 * t631 + t1561 * t632;
t1178 = -t1248 / 0.2e1;
t1177 = t1248 / 0.2e1;
t1176 = -t1247 / 0.2e1;
t1175 = t1247 / 0.2e1;
t1170 = -t1243 / 0.2e1;
t1169 = -t1243 / 0.4e1;
t1166 = -t1242 / 0.2e1;
t1164 = t1242 / 0.2e1;
t1162 = -t1241 / 0.2e1;
t1161 = -t1241 / 0.4e1;
t1160 = t1241 / 0.2e1;
t1157 = -t964 - t1450;
t1156 = -t966 - t1450;
t1149 = t1023 * t918 - t915;
t857 = t920 * t1235;
t1147 = t1027 * t916 - t857;
t1146 = t1421 + t1586;
t1139 = t1156 - t869;
t1134 = t1229 * t1451;
t1130 = t1428 + t1590 * t1288 - t1589 * t1221 / 0.2e1;
t1129 = t1172 * t1667 + t1675;
t1121 = t889 * t1175 + t839 * t1176 - t1310 / 0.2e1 + t1376 * t1178;
t1111 = -Icges(3,5) * t1023 - Icges(3,6) * t1026;
t1110 = Icges(4,5) * t1005 + Icges(4,6) * t1007;
t1104 = t1005 * t879 + t1007 * t877;
t1093 = t1027 * t539 + t1254;
t1090 = t1156 + t1380;
t1089 = -t1511 / 0.4e1 + t1129;
t1051 = t1005 * t1078 + t1320;
t348 = t1024 * t1051 + t789 * t903 - t791 * t904;
t1050 = t1005 * t1077 + t1319;
t349 = t1024 * t1050 + t790 * t903 - t792 * t904;
t78 = (t1098 - t427) * t1007 + (t1024 * t348 + t1027 * t349 + t561) * t1005;
t350 = t1027 * t1051 - t789 * t905 - t791 * t906;
t351 = t1027 * t1050 - t790 * t905 - t792 * t906;
t79 = (t1097 - t428) * t1007 + (t1024 * t350 + t1027 * t351 + t563) * t1005;
t1088 = t1160 * t1667 + t1164 * t239 + t1168 * t79 + t1172 * t78 + t1585;
t1086 = t1146 + t1584;
t196 = t1024 * t349 - t1027 * t348;
t197 = t1024 * t351 - t1027 * t350;
t1049 = t1005 * t1075 + t1316;
t380 = t1024 * t1049 + t818 * t942 - t820 * t943;
t1048 = t1005 * t1074 + t1315;
t381 = t1024 * t1048 + t819 * t942 - t821 * t943;
t221 = t1024 * t381 - t1027 * t380;
t382 = t1027 * t1049 - t818 * t944 - t820 * t945;
t383 = t1027 * t1048 - t819 * t944 - t821 * t945;
t222 = t1024 * t383 - t1027 * t382;
t928 = t1110 * t1024;
t929 = t1027 * t1110;
t1085 = (t197 + t222 - t1018 * t929 + (t1024 * t928 + t1572) * t1027) * t1464 + (t196 + t221 - t1019 * t928 + (t1027 * t929 + t1572) * t1024) * t1462;
t1083 = t1516 / 0.4e1 + t1088;
t1072 = -t1134 + t1362;
t1069 = t422 * t472 - t633 * t634 - t635 * t636;
t1068 = t553 * t591 - t749 * t750 - t751 * t752;
t1067 = t577 * t598 + t616 * t637 - t618 * t638;
t1066 = (-t955 + t958) * t1007 + t1650 * t1005;
t1064 = t1576 + t1664;
t1063 = t1675 + (t392 + t464) * t1168 + (t1667 + t391 + t463) * t1172;
t1062 = t839 * t1175 + t889 * t1176 + t1310 / 0.2e1 + t1376 * t1177;
t1061 = t197 * t1168 + t196 * t1172 + t78 * t1462 + t79 * t1464 + (t1265 - t1266) * t1466 + (t1024 * t466 - t1027 * t465) * t1164 + t1659 * t1160 + (t1024 * t496 - t1027 * t494) * t1467 - t1421;
t1060 = t1159 * t1633 + t1163 * t1634 + t1167 * t1636 + t1171 * t1637 + t1587;
t1053 = t79 * t1170 + t78 * t1174 + t164 * t1464 + t163 * t1462 + t205 * t1172 + t206 * t1168 + t239 * t1166 + t1667 * t1162 + (t1263 - t1264) * t1466 - t1585;
t1045 = -t1516 / 0.4e1 + t1053;
t1044 = t1063 - t1255;
t1043 = -t1478 / 0.4e1 + t1053;
t1039 = t1175 * t842 + t1178 * t840 + t1571;
t108 = (t1095 - t461) * t1007 + (t1024 * t380 + t1027 * t381 + t584) * t1005;
t109 = (t1094 - t462) * t1007 + (t1024 * t382 + t1027 * t383 + t586) * t1005;
t1038 = t108 * t1462 + t109 * t1464 + t222 * t1168 + t221 * t1172 + t1061 + (t1261 - t1262) * t1466 + (t1024 * t509 - t1027 * t508) * t1164 + t1660 * t1160 + (t1024 * t539 - t1027 * t537) * t1467 - t1586;
t1037 = t1060 + t1064 - t1145;
t1036 = t1161 * t1633 + t1165 * t1634 + t1169 * t1636 + t1173 * t1637 + t1064 + t1145 - t1587;
t1035 = t1060 + t1145 - t1576 + t1664;
t138 = (t1093 - t491) * t1007 + (t1024 * t406 + t1027 * t407 + t615) * t1005;
t183 = -t1007 * t514 + (t1024 * t398 + t1027 * t399) * t1005;
t184 = -t1007 * t515 + (t1024 * t400 + t1027 * t401) * t1005;
t298 = t1005 * t1093 - t1007 * t615;
t1034 = t109 * t1170 + (t1259 - t1260) * t1466 + t298 * t1647 + t248 * t1168 + t183 * t1462 + t1053 + t1668 * t1162 + t108 * t1174 + t247 * t1172 + t138 * t1646 - t130 * t1166 + t184 * t1464 - t1583;
t1033 = t1261 / 0.2e1 - t1042 - t1262 / 0.2e1 + t1265 / 0.2e1 - t1266 / 0.2e1 + (t1005 * t1366 + t1007 * t1364 + t1024 * t954 + t1027 * t1066 + t428 + t462) * t1464 + (-t1005 * t1367 + t1007 * t1365 + t1024 * t1066 - t1267 + t427 + t461) * t1462;
t1032 = t1254 / 0.2e1 + t1256 / 0.2e1 + t1104 * t1465 + t842 * t1176 + t840 * t1177 - t1571 + (t1104 - t537 - t494) * t1464;
t1031 = t1575 + t1574 + t1037 - t1582;
t1030 = t1161 * t1652 + t1165 * t1653 + t1169 * t1654 + t1173 * t1655 + t1036 + t1575 + t1582 + t1613;
t1029 = t1035 + t1574 + t1582 - t1612 + t1663;
t1003 = t1005 ^ 2;
t978 = -rSges(3,2) * t1023 + t1444;
t947 = t1111 * t1027;
t946 = t1111 * t1024;
t873 = t1157 * t1027;
t871 = t1157 * t1024;
t712 = t1139 * t1027;
t710 = t1139 * t1024;
t682 = -t1134 + t753;
t653 = -t1007 * t782 - t1243 * t938;
t652 = t1007 * t781 + t1245 * t938;
t651 = -t1237 * t918 + t1372;
t650 = -t1237 * t917 - t1373;
t649 = -t1238 * t918 - t1147;
t623 = t1090 * t1027;
t621 = t1090 * t1024;
t590 = (pkin(4) * t1003 * t1022 - t1005 * t894) * t1027 + t1383 * t1007;
t589 = -t1003 * t1220 - t1007 * t925 + t637;
t572 = t1072 + t1602;
t555 = -t1309 + (-t1022 * t1370 + t1025 * t1371) * t1005;
t545 = t1024 * t651 - t1027 * t650;
t544 = t1024 * t649 - t1027 * (-t1024 * (-t1026 * t919 + t1305) - t1027 * t915);
t541 = 0.4e1 * t1588;
t454 = t1072 + t1581;
t440 = 0.4e1 * t1428 + 0.4e1 * (t1024 * t1590 + t1027 * t1589) * t964;
t437 = -0.4e1 * t640 * t781 + 0.4e1 * t641 * t782;
t426 = 0.4e1 * t640 * t718 + 0.4e1 * t641 * t719;
t412 = 0.4e1 * t640 * t683 + 0.4e1 * t641 * t684;
t393 = -0.4e1 * t610 * t742 + 0.4e1 * t611 * t743;
t376 = 0.4e1 * t610 * t674 + 0.4e1 * t611 * t675;
t373 = -0.4e1 * t1383 * t611 + 0.4e1 * t610 * t664;
t365 = 0.4e1 * t610 * t654 + 0.4e1 * t611 * t655;
t347 = t553 * t1551 + (t1024 * t749 + t1027 * t751) * t1227;
t312 = t483 * t1551 + (t1024 * t709 + t1027 * t711) * t1227;
t307 = 0.4e1 * t1068;
t283 = (t649 - t857 + (t916 + t1305) * t1027 + t1373) * t1027 + t1372 * t1024;
t282 = (t1027 * t1149 - t1372 + t651) * t1027 + (t1024 * t1149 + t1147 + t650) * t1024;
t281 = t1518 * t393 + t1121;
t279 = 0.4e1 * t1131;
t274 = t1468 / 0.4e1;
t263 = t1469 / 0.4e1;
t251 = t422 * t1556 + (t1024 * t633 + t1027 * t635) * t1228;
t241 = 0.4e1 * t1067;
t215 = t377 * t1556 + (t1024 * t620 + t1027 * t622) * t1228;
t214 = 0.4e1 * t422 * t558 - 0.4e1 * t633 * t785 - 0.4e1 * t635 * t786;
t207 = t1474 / 0.4e1;
t198 = 0.4e1 * t377 * t558 - 0.4e1 * t620 * t785 - 0.4e1 * t622 * t786;
t177 = 0.4e1 * t1069;
t175 = t1475 / 0.4e1;
t172 = t1476 / 0.4e1;
t171 = 0.4e1 * t1133;
t169 = t1477 / 0.4e1;
t166 = t1518 * t373 + t1520 * t437 + t1121 - t1665;
t165 = 0.4e1 * t1132;
t141 = t1039 + t1497 / 0.4e1 + t426 * t1520 + t376 * t1518;
t134 = 0.4e1 * t1408 * t507;
t132 = t1480 / 0.4e1;
t122 = t1481 / 0.4e1;
t107 = (t974 / 0.2e1 + t973 / 0.2e1) * t1026 + t1039 + t412 * t1520 + t365 * t1518 - t1579;
t94 = t1503 / 0.4e1;
t80 = t1505 / 0.4e1;
t68 = t1507 / 0.4e1;
t64 = t1509 / 0.4e1;
t59 = t1518 * t251 + t1421;
t58 = t1518 * t215 + t1421;
t46 = t1518 * t241 + t1218;
t45 = t46 * qJD(5);
t44 = t1518 * t165 + t1218;
t43 = t1518 * t214 + t1520 * t347 + t1146;
t42 = t1518 * t198 + t1520 * t312 + t1146;
t41 = t1518 * t177 + t1520 * t307 + t1522 * t541 + t1085;
t40 = t41 * qJD(3);
t39 = t1518 * t171 + t1520 * t279 + t1522 * t440 + t1085;
t37 = t68 - t1509 / 0.4e1 + t1219;
t36 = t64 - t1507 / 0.4e1 + t1219;
t35 = t1478 / 0.4e1 + t1088;
t34 = t35 * qJD(5);
t33 = t207 - t1475 / 0.4e1 + t1089;
t32 = t175 - t1474 / 0.4e1 + t1089;
t31 = t175 + t207 + t1511 / 0.4e1 + t1044;
t29 = -t134 * t1518 + t1129;
t28 = t1042 + (t282 / 0.2e1 + t544 / 0.2e1) * t1024 + (t545 / 0.2e1 - t283 / 0.2e1) * t1027;
t27 = t132 - t1503 / 0.4e1 + t1083;
t26 = t94 - t1480 / 0.4e1 + t1083;
t25 = t122 - t1505 / 0.4e1 + t1083;
t24 = t80 - t1481 / 0.4e1 + t1083;
t23 = t1086 - t1614;
t22 = t1086 + t1614;
t21 = t64 + t68 - t1479 / 0.4e1 + t1061;
t20 = t1042 + t1615;
t19 = t1042 - t1615;
t18 = t1035 + t172 + t274;
t17 = t172 + t1037 - t1468 / 0.4e1;
t16 = t1036 - t1476 / 0.4e1 + t274;
t15 = (-t138 / 0.2e1 + t1279 / 0.2e1 - t1297 / 0.2e1) * t1007 + (t298 / 0.2e1 + t109 * t1461 + t108 * t1464) * t1005 + t1088 + t1583;
t14 = t15 * qJD(4);
t13 = t1035 + t169 + t263;
t12 = t169 + t1037 - t1469 / 0.4e1;
t11 = t1036 - t1477 / 0.4e1 + t263;
t10 = t1045 + t132 + t94;
t9 = t122 + t1045 + t80;
t8 = t1033 + t1216 + t1217;
t7 = t1038 + t1445 + t1446 - t1584;
t6 = t1029 + t1603 + t1649;
t5 = t1031 + t1605 - t1607 - t1603;
t4 = t1030 + t1603 - t1649;
t3 = t1029 + t1604 + t1648;
t2 = t1031 + t1606 - t1608 - t1604;
t1 = t1030 + t1604 - t1648;
t30 = [t107 * qJD(2) + t141 * qJD(3) + t166 * qJD(4) + t281 * qJD(5), t107 * qJD(1) + t8 * qJD(3) + t3 * qJD(4) + t13 * qJD(5) + (t610 * t623 + t611 * t621 - t620 * t655 - t622 * t654) * t1453 + (t640 * t712 + t641 * t710 - t683 * t711 - t684 * t709) * t1456 + (t873 * t783 + t871 * t784) * t1459 + (t1033 + (t1023 * t1358 + t1026 * t1356) * t1464 + t283 * t1461 + (t282 + t544) * t1465 + (-t1023 * t1359 + t1026 * t1357 + t545) * t1462 + (t1018 / 0.2e1 + t1019 / 0.2e1) * t1112 + ((-t845 * t978 - t952 * t976) * t1027 + (-t846 * t978 + t953 * t976) * t1024) * m(3)) * qJD(2), t141 * qJD(1) + t8 * qJD(2) + t1033 * qJD(3) + t6 * qJD(4) + t18 * qJD(5) + (-t633 * t675 - t635 * t674 + t516 / 0.2e1 + t517 / 0.2e1) * t1452 + (-t751 * t718 - t749 * t719 + t578 / 0.2e1 + t579 / 0.2e1) * t1455 + (t700 / 0.2e1 + t701 / 0.2e1 + (t1024 * t940 - t1027 * t939) * t963) * t1458, t166 * qJD(1) + t3 * qJD(2) + t6 * qJD(3) + t31 * qJD(5) + (t1063 + (-t532 - t555) * t1007 + (t134 / 0.4e1 + t506 * t664 + t507 * t1383 + t589 * t610 + t590 * t611) * m(6) + (-t629 * t781 - t631 * t782 + t652 * t640 + t653 * t641) * m(5) + ((t515 / 0.2e1 + t434 / 0.2e1) * t1027 + (t433 / 0.2e1 + t514 / 0.2e1) * t1024) * t1005) * qJD(4), t281 * qJD(1) + t13 * qJD(2) + t18 * qJD(3) + t31 * qJD(4) + (t1044 + (-t616 * t742 - t618 * t743 + t518 / 0.2e1 + t519 / 0.2e1) * m(6)) * qJD(5); t28 * qJD(2) + t19 * qJD(3) + t1 * qJD(4) + t11 * qJD(5) + t365 * t1658 + t412 * t1657 + (t1032 - (t973 + t974) * t1026 / 0.2e1 + t1579) * qJD(1), t28 * qJD(1) + ((t377 * t454 - t620 * t621 - t622 * t623) * m(6) + (t483 * t572 - t709 * t710 - t711 * t712) * m(5) + m(4) * (-t1589 * t873 - t1590 * t871 + t603 * t682) + (t1018 * t947 + (-t1024 * t946 + t1573) * t1027) * t1464 + 0.4e1 * ((t1024 * (rSges(3,1) * t1235 - t1249) + t1027 * (rSges(3,1) * t1233 + rSges(3,3) * t1024 - t995)) * (-t1024 * t952 - t1027 * t953) + t1229 * t978 * t976) * t1523 + (t1019 * t946 + (-t1027 * t947 + t1573) * t1024) * t1462 + t1085) * qJD(2) + t39 * qJD(3) + t42 * qJD(4) + t58 * qJD(5), t19 * qJD(1) + t39 * qJD(2) + t1085 * qJD(3) + t23 * qJD(4) + t36 * qJD(5) + (-t177 / 0.4e1 + t1069 + t1133) * t1452 + (-t307 / 0.4e1 + t1068 + t1131) * t1455 + (-t541 / 0.4e1 + t1130 + t1588) * t1458, t1 * qJD(1) + t42 * qJD(2) + t23 * qJD(3) + (t1034 + (t1536 * t652 + t1537 * t653 + t1555 * t483 + t1206) * t1521 + (t1553 * t589 + t1554 * t590 + t1565 * t377 + t1213) * t1519) * qJD(4) + t9 * qJD(5), t11 * qJD(1) + t58 * qJD(2) + t36 * qJD(3) + t9 * qJD(4) + (t1043 + (t1207 + t1214) * t1519) * qJD(5); (-t1497 / 0.4e1 + t1032) * qJD(1) + t20 * qJD(2) + t1042 * qJD(3) + t4 * qJD(4) + t16 * qJD(5) + t376 * t1658 + t426 * t1657, t20 * qJD(1) + t1085 * qJD(2) + t40 + t22 * qJD(4) + t37 * qJD(5) + (-t171 / 0.4e1 + t422 * t454 - t621 * t633 - t623 * t635 + t1133) * t1453 + (-t279 / 0.4e1 + t553 * t572 - t749 * t710 - t751 * t712 + t1131) * t1456 + (-t440 / 0.4e1 + t703 * t682 + (-t1024 * t871 - t1027 * t873) * t963 + t1130) * t1459, qJD(1) * t1042 + qJD(2) * t41 + qJD(4) * t43 + qJD(5) * t59 + t40, t4 * qJD(1) + t22 * qJD(2) + t43 * qJD(3) + (t1034 + (t1530 * t652 + t1531 * t653 + t1555 * t553 + t1206) * t1521 + (t1548 * t589 + t1550 * t590 + t1565 * t422 + t1213) * t1519) * qJD(4) + t10 * qJD(5), t16 * qJD(1) + t37 * qJD(2) + t59 * qJD(3) + t10 * qJD(4) + (t1043 + (t1207 + t1212) * t1519) * qJD(5); t2 * qJD(2) + t5 * qJD(3) + t29 * qJD(4) + t32 * qJD(5) + (-t373 / 0.4e1 + t1408 * t611) * t1454 + t437 * t1657 + (t1062 + t1665) * qJD(1), t2 * qJD(1) + t1038 * qJD(2) + t7 * qJD(3) + t14 + t24 * qJD(5) + (t436 * t454 + t506 * t623 - t507 * t621 + t190 / 0.2e1 + t303 / 0.2e1 + t304 / 0.2e1 - t198 / 0.4e1) * t1453 + (t572 * t595 + t629 * t712 - t631 * t710 + t323 / 0.2e1 + t473 / 0.2e1 + t474 / 0.2e1 - t312 / 0.4e1) * t1456, t5 * qJD(1) + t7 * qJD(2) + (t1038 + (-t347 / 0.4e1 + t140 / 0.2e1) * m(5) + (-t214 / 0.4e1 + t53 / 0.2e1) * m(6)) * qJD(3) + t14 + t26 * qJD(5), t29 * qJD(1) + ((t436 * t540 + t506 * t589 - t507 * t590) * m(6) + t1007 ^ 2 * t555 / 0.2e1 + (t595 * t619 + t629 * t652 - t631 * t653) * m(5) + (t184 * t1461 + t183 * t1464 + (t433 * t1024 + t434 * t1027) * t1466) * t1005 + t1218) * qJD(4) + t44 * qJD(5) + t1230 * t15, t32 * qJD(1) + t24 * qJD(2) + t26 * qJD(3) + t44 * qJD(4) + ((-t241 / 0.4e1 + t1067 + t1132) * m(6) + t1218) * qJD(5); qJD(1) * t1062 + qJD(2) * t12 + qJD(3) * t17 + qJD(4) * t33 + qJD(5) * t1129 + t1658 * t393, t12 * qJD(1) + t1061 * qJD(2) + t21 * qJD(3) + t25 * qJD(4) + t34 + (t454 * t577 + t616 * t623 - t618 * t621 + t261 / 0.2e1 + t420 / 0.2e1 + t421 / 0.2e1 - t215 / 0.4e1) * t1453, t17 * qJD(1) + t21 * qJD(2) + ((t69 / 0.2e1 - t251 / 0.4e1) * m(6) + t1061) * qJD(3) + t27 * qJD(4) + t34, t33 * qJD(1) + t25 * qJD(2) + t27 * qJD(3) + ((t540 * t577 + t589 * t616 - t590 * t618 - t165 / 0.4e1 + t1132) * m(6) + t1218) * qJD(4) + t45, qJD(1) * t1129 + qJD(4) * t46 + t1230 * t35 + t45;];
Cq = t30;
