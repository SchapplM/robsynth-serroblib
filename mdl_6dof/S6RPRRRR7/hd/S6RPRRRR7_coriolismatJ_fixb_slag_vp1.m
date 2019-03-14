% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRRR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR7_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:59
% EndTime: 2019-03-09 07:17:53
% DurationCPUTime: 45.28s
% Computational Cost: add. (208911->1185), mult. (194053->1604), div. (0->0), fcn. (205951->10), ass. (0->771)
t878 = qJ(3) + qJ(4);
t871 = qJ(5) + t878;
t857 = cos(t871);
t879 = sin(qJ(6));
t882 = cos(qJ(6));
t740 = (-Icges(7,5) * t879 - Icges(7,6) * t882) * t857;
t856 = sin(t871);
t1157 = t856 * t740;
t1215 = Icges(7,4) * t882;
t953 = -Icges(7,2) * t879 + t1215;
t682 = Icges(7,6) * t856 + t857 * t953;
t1216 = Icges(7,4) * t879;
t957 = Icges(7,1) * t882 - t1216;
t684 = Icges(7,5) * t856 + t857 * t957;
t743 = (-Icges(7,2) * t882 - t1216) * t857;
t746 = (-Icges(7,1) * t879 - t1215) * t857;
t1446 = -t1157 / 0.2e1 - (-(t684 / 0.2e1 + t743 / 0.2e1) * t879 + (t746 / 0.2e1 - t682 / 0.2e1) * t882) * t857;
t1440 = -m(7) * qJD(1) / 0.4e1;
t884 = cos(qJ(1));
t1146 = t857 * t884;
t881 = sin(qJ(1));
t1148 = t857 * t881;
t1126 = t882 * t884;
t1129 = t881 * t879;
t788 = -t1129 * t856 + t1126;
t1128 = t881 * t882;
t1139 = t879 * t884;
t789 = t1128 * t856 + t1139;
t584 = Icges(7,5) * t789 + Icges(7,6) * t788 - Icges(7,3) * t1148;
t1217 = Icges(7,4) * t789;
t587 = Icges(7,2) * t788 - Icges(7,6) * t1148 + t1217;
t764 = Icges(7,4) * t788;
t590 = Icges(7,1) * t789 - Icges(7,5) * t1148 + t764;
t790 = t1139 * t856 + t1128;
t791 = t1126 * t856 - t1129;
t355 = t584 * t1146 + t790 * t587 - t791 * t590;
t1444 = t355 * t881;
t1443 = t355 * t884;
t586 = -Icges(7,5) * t791 + Icges(7,6) * t790 + Icges(7,3) * t1146;
t1178 = t586 * t856;
t766 = Icges(7,4) * t791;
t589 = Icges(7,2) * t790 + Icges(7,6) * t1146 - t766;
t765 = Icges(7,4) * t790;
t591 = Icges(7,1) * t791 - Icges(7,5) * t1146 - t765;
t1434 = t589 * t879 + t591 * t882;
t376 = t1434 * t857 - t1178;
t1306 = t881 / 0.2e1;
t1305 = -t884 / 0.2e1;
t1088 = t788 * t587 + t789 * t590;
t941 = -t589 * t790 - t591 * t791;
t1442 = t1088 + (-t584 * t881 - t586 * t884) * t857 + t941;
t1304 = t884 / 0.2e1;
t1307 = -t881 / 0.2e1;
t1219 = Icges(6,4) * t856;
t954 = Icges(6,2) * t857 + t1219;
t700 = -Icges(6,6) * t881 + t884 * t954;
t1218 = Icges(6,4) * t857;
t958 = Icges(6,1) * t856 + t1218;
t702 = -Icges(6,5) * t881 + t884 * t958;
t1393 = (t700 * t857 + t702 * t856) * t884;
t354 = -t1148 * t586 + t788 * t589 - t591 * t789;
t1398 = t354 * t881;
t353 = -t1148 * t584 + t1088;
t232 = t353 * t884 + t1398;
t356 = t1146 * t586 - t941;
t233 = t356 * t881 + t1443;
t1420 = t233 * t1305 + t1306 * t232;
t699 = Icges(6,6) * t884 + t881 * t954;
t1154 = t856 * t881;
t826 = Icges(6,4) * t1148;
t701 = Icges(6,1) * t1154 + Icges(6,5) * t884 + t826;
t936 = -t699 * t857 - t701 * t856;
t1439 = t884 * t936;
t947 = Icges(6,5) * t856 + Icges(6,6) * t857;
t1424 = t881 * t947;
t475 = t699 * t1148 + t701 * t1154 + t884 * (Icges(6,3) * t884 + t1424);
t1422 = t884 * t947;
t698 = -Icges(6,3) * t881 + t1422;
t476 = -t700 * t1148 - t702 * t1154 - t884 * t698;
t478 = -t881 * t698 + t1393;
t876 = t881 ^ 2;
t1107 = t356 + t1442;
t91 = t1107 * t884 + t1398;
t1104 = -t353 + t1442;
t92 = t1104 * t881 - t1443;
t909 = (t475 * t884 + t476 * t881) * t1307 + ((t476 - t1439) * t881 + (t478 - t1393 + (t698 + t936) * t881 + t475) * t884 + t91) * t1306 + (t478 * t881 + t876 * t698 + t92 + (t476 + (t698 - t936) * t884 + t1439) * t884) * t1304 - t1420;
t1391 = t791 * rSges(7,1) - t790 * rSges(7,2);
t599 = rSges(7,3) * t1146 - t1391;
t964 = rSges(7,1) * t882 - rSges(7,2) * t879;
t687 = rSges(7,3) * t856 + t857 * t964;
t1441 = t687 * t1146 - t599 * t856;
t1328 = m(6) / 0.4e1;
t1326 = m(7) / 0.4e1;
t859 = sin(t878);
t1221 = Icges(5,4) * t859;
t860 = cos(t878);
t955 = Icges(5,2) * t860 + t1221;
t729 = Icges(5,6) * t884 + t881 * t955;
t1143 = t859 * t881;
t1140 = t860 * t881;
t847 = Icges(5,4) * t1140;
t731 = Icges(5,1) * t1143 + Icges(5,5) * t884 + t847;
t933 = -t729 * t860 - t731 * t859;
t1438 = t884 * t933;
t798 = -Icges(6,2) * t856 + t1218;
t1437 = t798 + t958;
t1220 = Icges(5,4) * t860;
t818 = -Icges(5,2) * t859 + t1220;
t959 = Icges(5,1) * t859 + t1220;
t1436 = t818 + t959;
t880 = sin(qJ(3));
t1244 = pkin(3) * t880;
t967 = rSges(5,1) * t859 + rSges(5,2) * t860;
t1416 = t967 + t1244;
t715 = t1416 * t884;
t1034 = pkin(4) * t1143;
t885 = -pkin(8) - pkin(7);
t1046 = -pkin(9) + t885;
t854 = t884 * t1046;
t858 = t884 * t885;
t664 = -t854 + t858 + t1034;
t966 = rSges(6,1) * t856 + rSges(6,2) * t857;
t776 = t881 * t966;
t707 = t884 * rSges(6,3) + t776;
t1074 = -t664 - t707;
t1242 = pkin(3) * t884;
t1033 = t880 * t1242;
t1241 = pkin(4) * t859;
t971 = t1241 + t1244;
t1054 = t881 * t1046 + t884 * t971;
t637 = t884 * (t881 * t885 + t1033 - t1054);
t1162 = t884 * t966;
t908 = -t881 * rSges(6,3) + t1162;
t685 = t884 * t908;
t1076 = t637 - t685;
t752 = t884 * (t1033 + (pkin(7) + t885) * t881);
t1138 = t880 * t881;
t795 = pkin(3) * t1138 - pkin(7) * t884 - t858;
t383 = -t752 + (-t795 + t1074) * t881 + t1076;
t1153 = t856 * t884;
t831 = rSges(6,1) * t1146;
t751 = -rSges(6,2) * t1153 + t831;
t712 = t884 * t751;
t749 = rSges(6,1) * t1148 - rSges(6,2) * t1154;
t601 = -t749 * t881 - t712;
t1199 = t383 * t601;
t331 = 0.2e1 * t1199;
t1035 = -0.2e1 * t776;
t1230 = rSges(6,2) * t856;
t802 = rSges(6,1) * t857 - t1230;
t777 = t881 * t802;
t849 = pkin(4) * t1140;
t689 = t849 + t777;
t883 = cos(qJ(3));
t1127 = t881 * t883;
t852 = pkin(3) * t1127;
t654 = t852 + t689;
t577 = t654 * t1035;
t1037 = -0.2e1 * t1162;
t1240 = pkin(4) * t860;
t1243 = pkin(3) * t883;
t825 = t1240 + t1243;
t656 = (t802 + t825) * t884;
t578 = t656 * t1037;
t1008 = t331 + t577 + t578;
t597 = t789 * rSges(7,1) + t788 * rSges(7,2) - rSges(7,3) * t1148;
t1392 = -pkin(10) * t1148 + t597;
t1080 = -pkin(5) * t1154 - t1392;
t1000 = -t664 + t1080;
t836 = pkin(10) * t1146;
t1085 = (-pkin(5) * t1153 + t599 + t836) * t884;
t1001 = t637 + t1085;
t299 = -t752 + (-t795 + t1000) * t881 + t1001;
t1149 = t857 * t879;
t1031 = rSges(7,2) * t1149;
t1147 = t857 * t882;
t1032 = rSges(7,1) * t1147;
t649 = rSges(7,3) * t1154 + (-t1031 + t1032) * t881;
t1075 = -pkin(5) * t1148 - pkin(10) * t1154 - t649;
t1050 = -pkin(5) * t1146 - pkin(10) * t1153;
t1051 = rSges(7,3) * t1153 + t884 * t1032;
t650 = t1031 * t884 - t1051;
t571 = -t650 - t1050;
t1077 = t571 * t884;
t430 = t1075 * t881 - t1077;
t1209 = t299 * t430;
t238 = 0.2e1 * t1209;
t669 = t881 * t687;
t804 = pkin(5) * t857 + pkin(10) * t856;
t594 = t881 * t804 + t669;
t559 = t849 + t594;
t537 = t852 + t559;
t1229 = rSges(7,3) * t857;
t686 = -t856 * t964 + t1229;
t668 = t881 * t686;
t1239 = pkin(5) * t856;
t803 = pkin(10) * t857 - t1239;
t593 = t881 * t803 + t668;
t1191 = t537 * t593;
t451 = 0.2e1 * t1191;
t1070 = t687 + t804;
t540 = (t825 + t1070) * t884;
t1071 = -t686 - t803;
t595 = t1071 * t884;
t1189 = t540 * t595;
t452 = -0.2e1 * t1189;
t1015 = t238 + t451 + t452;
t1333 = -0.2e1 * t884;
t688 = -t776 - t1034;
t1337 = 0.2e1 * t688;
t850 = t884 * t1241;
t690 = t850 + t1162;
t1078 = t802 * t690 * t1333 + t777 * t1337;
t596 = t1070 * t884;
t1340 = -0.2e1 * t596;
t1342 = 0.2e1 * t594;
t558 = -t1034 + t593;
t560 = t850 + t595;
t1093 = t560 * t1340 + t558 * t1342;
t391 = t1080 * t881 + t1085;
t1359 = 0.2e1 * t391;
t877 = t884 ^ 2;
t1049 = t876 + t877;
t974 = t1049 * t1240;
t405 = -t974 + t430;
t281 = t405 * t1359;
t569 = -t707 * t881 - t685;
t1346 = 0.2e1 * t569;
t538 = t601 - t974;
t434 = t538 * t1346;
t1401 = (t434 + t1008 + t1078) * t1328 + (t281 + t1015 + t1093) * t1326;
t425 = t1074 * t881 + t1076;
t1196 = t425 * t601;
t691 = (t802 + t1240) * t884;
t1006 = t689 * t1035 + t691 * t1037 + 0.2e1 * t1196;
t328 = t1000 * t881 + t1001;
t561 = (t1070 + t1240) * t884;
t1378 = t328 * t430 + t559 * t593 - t561 * t595;
t1013 = 0.2e1 * t1378;
t1402 = (t1006 + t1008) * t1328 + (t1013 + t1015) * t1326;
t1435 = t1401 - t1402;
t1330 = m(5) / 0.4e1;
t1329 = m(6) / 0.2e1;
t1433 = -t856 / 0.2e1;
t1309 = t856 / 0.2e1;
t1432 = -t857 / 0.2e1;
t1431 = t857 / 0.2e1;
t1429 = t881 / 0.4e1;
t1428 = t884 / 0.4e1;
t1045 = 0.2e1 * m(7);
t1427 = t1045 / 0.4e1;
t1426 = t233 + t92;
t946 = Icges(7,5) * t882 - Icges(7,6) * t879;
t680 = Icges(7,3) * t856 + t857 * t946;
t423 = t1146 * t680 + t682 * t790 - t684 * t791;
t1425 = t423 * t856;
t949 = Icges(5,5) * t859 + Icges(5,6) * t860;
t1423 = t881 * t949;
t1421 = t884 * t949;
t1161 = t884 * t967;
t640 = t682 * t881;
t642 = t684 * t881;
t942 = -t587 * t879 + t590 * t882;
t914 = t680 * t881 - t942;
t275 = (-t640 * t879 + t642 * t882 + t584) * t857 + t914 * t856;
t1171 = t680 * t856;
t679 = Icges(7,3) * t857 - t856 * t946;
t1169 = t684 * t882;
t1170 = t682 * t879;
t937 = t1169 - t1170;
t912 = t679 - t937;
t1377 = t857 * t912 - t1171;
t681 = Icges(7,6) * t857 - t856 * t953;
t683 = Icges(7,5) * t857 - t856 * t957;
t313 = -t1377 * t881 + t681 * t788 + t683 * t789;
t1419 = t275 + t313;
t641 = t682 * t884;
t643 = t684 * t884;
t913 = -t680 * t884 + t1434;
t276 = (t641 * t879 - t643 * t882 + t586) * t857 + t913 * t856;
t314 = t1377 * t884 + t681 * t790 - t683 * t791;
t1418 = t276 + t314;
t861 = t884 * qJ(2);
t1383 = -t881 * pkin(1) + t861;
t924 = t1054 + t1383;
t472 = (-t1229 + t1239) * t884 - t836 + t924 + t1391;
t1044 = qJD(3) + qJD(4);
t982 = qJD(5) + t1044;
t822 = rSges(5,1) * t860 - rSges(5,2) * t859;
t969 = (t822 + t1243) * t884;
t1165 = t969 * t881;
t1172 = t656 * t881;
t1188 = t540 * t881;
t521 = -0.2e1 * t1188;
t635 = -0.2e1 * t1172;
t692 = -0.2e1 * t1165;
t1004 = (t521 + 0.2e1 * t1188) * t1326 + (t635 + 0.2e1 * t1172) * t1328 + (t692 + 0.2e1 * t1165) * t1330;
t1332 = 0.2e1 * t884;
t714 = t822 * t881 + t852;
t1005 = (t1332 * t537 + t521) * t1326 + (t1332 * t654 + t635) * t1328 + (t1332 * t714 + t692) * t1330;
t1415 = t1004 - t1005;
t1334 = 0.2e1 * t881;
t1036 = t967 * t1334;
t873 = t881 * rSges(5,3);
t612 = t861 - t873 + (-pkin(1) + t885) * t881 + t715;
t553 = t612 * t1036;
t613 = -t858 + (rSges(5,3) + pkin(1)) * t884 + (qJ(2) + t1416) * t881;
t554 = 0.2e1 * t613 * t1161;
t1087 = -t553 + t554;
t575 = t908 + t924;
t500 = t575 * t1337;
t923 = qJ(2) + t971;
t576 = -t854 + (pkin(1) + rSges(6,3)) * t884 + (t923 + t966) * t881;
t1403 = 0.2e1 * t576;
t501 = t690 * t1403;
t1090 = t500 + t501;
t1350 = 0.2e1 * t558;
t386 = t472 * t1350;
t1348 = 0.2e1 * t560;
t1238 = t884 * pkin(1);
t473 = t1238 - t854 + (t923 + t1239) * t881 + t1392;
t387 = t473 * t1348;
t1101 = t386 + t387;
t1347 = -0.2e1 * t561;
t1349 = 0.2e1 * t559;
t813 = t881 * t825;
t522 = -t1075 + t813;
t815 = t884 * t825;
t523 = t571 + t815;
t786 = rSges(5,1) * t1140 - rSges(5,2) * t1143;
t708 = t786 + t852;
t568 = -0.2e1 * t708 * t884 + 0.2e1 * t1165;
t657 = t749 + t813;
t658 = t751 + t815;
t1019 = (t1347 * t522 + t1349 * t523 + t1101) * t1326 + (-0.2e1 * t657 * t691 + 0.2e1 * t658 * t689 + t1090) * t1328 + (t568 * t822 + t1087) * t1330;
t547 = (-t1031 + t1240) * t884 - t1050 + t1051;
t1351 = 0.2e1 * t547;
t1352 = -0.2e1 * t540;
t546 = -t1075 + t849;
t674 = t749 + t849;
t675 = t831 + (-t1230 + t1240) * t884;
t787 = t822 * t884;
t1020 = (t1351 * t537 + t1352 * t546 + t1101) * t1326 + (0.2e1 * t654 * t675 - 0.2e1 * t656 * t674 + t1090) * t1328 + (0.2e1 * t714 * t787 - 0.2e1 * t786 * t969 + t1087) * t1330;
t1414 = t1019 - t1020;
t1166 = t691 * t881;
t1183 = t561 * t881;
t550 = -0.2e1 * t1183;
t665 = -0.2e1 * t1166;
t1097 = (t550 + 0.2e1 * t1183) * t1326 + (t665 + 0.2e1 * t1166) * t1328;
t1098 = (t1332 * t559 + t550) * t1326 + (t1332 * t689 + t665) * t1328;
t1413 = t1097 - t1098;
t514 = t575 * t1035;
t515 = t1162 * t1403;
t1089 = t514 + t515;
t1343 = 0.2e1 * t593;
t399 = t472 * t1343;
t1341 = 0.2e1 * t595;
t400 = t473 * t1341;
t1099 = t399 + t400;
t541 = -0.2e1 * t674 * t884 + 0.2e1 * t675 * t881;
t1109 = (t1340 * t546 + t1342 * t547 + t1099) * t1326 + (t541 * t802 + t1089) * t1328;
t1335 = 0.2e1 * t751;
t1336 = -0.2e1 * t749;
t1344 = 0.2e1 * t571;
t1345 = 0.2e1 * t1075;
t1110 = (t1344 * t559 + t1345 * t561 + t1099) * t1326 + (t1335 * t689 + t1336 * t691 + t1089) * t1328;
t1412 = t1109 - t1110;
t507 = -0.2e1 * t657 * t884 + 0.2e1 * t658 * t881;
t1111 = (t1340 * t522 + t1342 * t523 + t1099) * t1326 + (t507 * t802 + t1089) * t1328;
t1112 = (t1344 * t537 + t1345 * t540 + t1099) * t1326 + (t1335 * t654 + t1336 * t656 + t1089) * t1328;
t1411 = t1111 - t1112;
t1331 = m(4) / 0.4e1;
t1325 = -pkin(1) - pkin(7);
t968 = rSges(4,1) * t880 + rSges(4,2) * t883;
t907 = -t881 * rSges(4,3) + t884 * t968;
t666 = t1325 * t881 + t861 + t907;
t667 = (rSges(4,3) - t1325) * t884 + (qJ(2) + t968) * t881;
t845 = rSges(4,1) * t883 - rSges(4,2) * t880;
t811 = t845 * t881;
t812 = t845 * t884;
t1222 = Icges(4,4) * t883;
t840 = -Icges(4,2) * t880 + t1222;
t960 = Icges(4,1) * t880 + t1222;
t1410 = (t960 / 0.2e1 + t840 / 0.2e1) * t883 - 0.4e1 * t1328 * (t575 * t658 + t576 * t657) - 0.4e1 * t1330 * (t612 * t969 + t613 * t708) - 0.4e1 * t1331 * (t666 * t812 + t667 * t811);
t1066 = t798 * t884 + t702;
t1067 = -Icges(6,2) * t1154 + t701 + t826;
t800 = Icges(6,1) * t857 - t1219;
t1068 = -t800 * t884 + t700;
t1069 = -t800 * t881 + t699;
t1409 = -(t1068 * t881 - t1069 * t884) * t856 + (t1066 * t881 - t1067 * t884) * t857;
t732 = -Icges(5,5) * t881 + t884 * t959;
t1062 = t818 * t884 + t732;
t1063 = -Icges(5,2) * t1143 + t731 + t847;
t730 = -Icges(5,6) * t881 + t884 * t955;
t820 = Icges(5,1) * t860 - t1221;
t1064 = -t820 * t884 + t730;
t1065 = -t820 * t881 + t729;
t1408 = -(t1064 * t881 - t1065 * t884) * t859 + (t1062 * t881 - t1063 * t884) * t860;
t774 = -Icges(4,5) * t881 + t884 * t960;
t1057 = t840 * t884 + t774;
t851 = Icges(4,4) * t1127;
t773 = Icges(4,1) * t1138 + Icges(4,5) * t884 + t851;
t1058 = -Icges(4,2) * t1138 + t773 + t851;
t1223 = Icges(4,4) * t880;
t956 = Icges(4,2) * t883 + t1223;
t772 = -Icges(4,6) * t881 + t884 * t956;
t842 = Icges(4,1) * t883 - t1223;
t1059 = -t842 * t884 + t772;
t771 = Icges(4,6) * t884 + t881 * t956;
t1060 = -t842 * t881 + t771;
t1407 = -(t1059 * t881 - t1060 * t884) * t880 + (t1057 * t881 - t1058 * t884) * t883;
t332 = (-t681 * t879 + t683 * t882 + t680) * t857 + t912 * t856;
t1180 = t584 * t856;
t375 = t857 * t942 + t1180;
t421 = -t1148 * t680 + t682 * t788 + t684 * t789;
t461 = t857 * t937 + t1171;
t1405 = t1309 * t332 + t1431 * t461 + (t375 + t421) * t1154 / 0.4e1 - (-t376 + t423) * t1153 / 0.4e1;
t1404 = 0.4e1 * m(7);
t1397 = t354 * t884;
t1124 = t883 * t772;
t1395 = (t880 * t774 + t1124) * t884;
t1394 = (t730 * t860 + t732 * t859) * t884;
t696 = t884 * (-t873 + t1161);
t733 = rSges(5,3) * t884 + t881 * t967;
t581 = -t733 * t881 - t696;
t611 = -t786 * t881 - t884 * t787;
t1388 = -t1049 * t822 * t967 + t581 * t611;
t1387 = t689 * t688 - t691 * t690;
t988 = t654 * t688 - t656 * t690;
t1386 = t559 * t558 - t561 * t560;
t989 = t537 * t558 - t540 * t560;
t1384 = t876 / 0.2e1 + t877 / 0.2e1;
t1082 = -Icges(7,2) * t789 + t590 + t764;
t1084 = -Icges(7,1) * t788 + t1217 + t587;
t617 = Icges(7,5) * t788 - Icges(7,6) * t789;
t277 = t1082 * t788 - t1084 * t789 - t1148 * t617;
t1081 = Icges(7,2) * t791 - t591 + t765;
t1083 = -Icges(7,1) * t790 + t589 - t766;
t618 = Icges(7,5) * t790 + Icges(7,6) * t791;
t278 = t1081 * t788 - t1083 * t789 - t1148 * t618;
t149 = t277 * t884 + t278 * t881;
t279 = t1082 * t790 + t1084 * t791 + t1146 * t617;
t280 = t1081 * t790 + t1083 * t791 + t1146 * t618;
t150 = t279 * t884 + t280 * t881;
t1113 = t1304 * t149 + t1306 * t150;
t930 = -t749 * t884 + t751 * t881;
t1095 = (t1075 * t884 + t571 * t881) * t1427 + t930 * t1329;
t1163 = t966 * t802;
t1182 = t569 * t601;
t1379 = t391 * t430 + t594 * t593 - t596 * t595;
t975 = 0.4e1 * t1049;
t1381 = -t1379 * t1404 / 0.4e1 - m(6) * (-t1163 * t975 + 0.4e1 * t1182) / 0.4e1;
t711 = t884 * (t1242 * t883 - t815);
t738 = -t852 + t813;
t378 = t711 + (-t738 + t1075) * t881 - t1077;
t481 = t711 - t712 + (-t738 - t749) * t881;
t961 = t1006 + t1078;
t962 = t1013 + t1093;
t1380 = (t1346 * t481 + t961) * t1328 + (t1359 * t378 + t962) * t1326;
t1376 = t857 * t913 - t1178;
t1375 = t857 * t914 - t1180;
t927 = -t786 * t884 + t787 * t881;
t1003 = (-t546 * t884 + t547 * t881) * t1427 + t541 * t1328 + m(5) * t927 / 0.2e1;
t1373 = t1436 * t859 - (-t955 + t820) * t860;
t1372 = t1437 * t856 - (-t954 + t800) * t857;
t492 = t729 * t1140 + t731 * t1143 + t884 * (Icges(5,3) * t884 + t1423);
t728 = -Icges(5,3) * t881 + t1421;
t493 = -t730 * t1140 - t732 * t1143 - t884 * t728;
t495 = -t881 * t728 + t1394;
t1365 = (t492 * t884 + t493 * t881) * t1307 + ((t493 - t1438) * t881 + (t495 - t1394 + (t728 + t933) * t881 + t492) * t884) * t1306 + (t495 * t881 + t876 * t728 + (t493 + (t728 - t933) * t884 + t1438) * t884) * t1304;
t1364 = (-t820 / 0.2e1 + t955 / 0.2e1) * t859 - t1436 * t860 / 0.2e1;
t301 = t618 * t856 + (-t1081 * t879 - t1083 * t882) * t857;
t1206 = t301 * t881;
t300 = t617 * t856 + (-t1082 * t879 - t1084 * t882) * t857;
t1207 = t300 * t884;
t1072 = t684 + t743;
t1073 = t682 - t746;
t348 = t1072 * t788 - t1073 * t789 - t1148 * t740;
t349 = t1072 * t790 + t1073 * t791 + t1146 * t740;
t979 = t1207 / 0.4e1 + t1206 / 0.4e1 + t349 * t1429 + t348 * t1428;
t926 = -t811 * t884 + t812 * t881;
t978 = (-t522 * t884 + t523 * t881) * t1427 + t507 * t1328 + t568 * t1330 + m(4) * t926 / 0.2e1;
t944 = t356 * t884 - t1444;
t174 = t857 * t944 + t1425;
t1121 = t884 * t174;
t412 = t421 * t856;
t945 = -t353 * t881 + t1397;
t173 = t857 * t945 + t412;
t1135 = t881 * t173;
t68 = t412 + (-t1107 * t881 + t1397) * t857;
t69 = -t1425 + (t1104 * t884 + t1444) * t857;
t1363 = t1121 / 0.4e1 - t1135 / 0.4e1 + t68 * t1429 + t69 * t1428;
t900 = -t681 * t1149 / 0.2e1 + t683 * t1147 / 0.2e1 + t680 * t1431 + (t1169 + t800) * t1433 + t1437 * t1432 + (t1170 + t679 + t954) * t1309;
t1361 = 0.4e1 * t328;
t939 = t597 * t884 + t599 * t881;
t339 = (-t649 * t884 - t650 * t881) * t857 + t939 * t856;
t1360 = 0.2e1 * t339;
t1358 = 0.4e1 * t425;
t444 = t939 * t857;
t1357 = -0.2e1 * t444;
t623 = rSges(7,1) * t788 - rSges(7,2) * t789;
t624 = rSges(7,1) * t790 + rSges(7,2) * t791;
t938 = -t623 * t884 - t624 * t881;
t474 = t938 * t857;
t1355 = 0.2e1 * t474;
t482 = t597 * t856 + t669 * t857;
t1354 = 0.2e1 * t482;
t1353 = 0.2e1 * t537;
t1339 = -0.2e1 * t623;
t1338 = -0.2e1 * t624;
t1327 = m(7) / 0.2e1;
t203 = t299 * t1360;
t380 = (t597 + t668) * t857 + (t649 - t669) * t856;
t321 = t380 * t1352;
t381 = (t686 * t884 - t599) * t857 + (-t687 * t884 - t650) * t856;
t322 = t381 * t1353;
t1018 = t203 + t321 + t322;
t1100 = t482 * t1348 + t1350 * t1441;
t302 = t405 * t1357;
t1324 = m(7) * (t302 + t1018 + t1100);
t1017 = t380 * t1347 + t381 * t1349 + t328 * t1360;
t963 = t1017 + t1100;
t1323 = m(7) * (t1357 * t378 + t963);
t1010 = t482 * t1341 + t1343 * t1441 + t430 * t1357;
t1321 = m(7) * (t1010 + t1018);
t1320 = m(7) * (t1010 + t1017);
t74 = t1340 * t380 + t1342 * t381 + t1359 * t339 + t1010;
t1318 = m(7) * t74;
t1282 = 0.4e1 * m(6) * (t575 * t751 + t576 * t749);
t750 = (-rSges(7,1) * t879 - rSges(7,2) * t882) * t857;
t1164 = t750 * t884;
t1039 = 0.2e1 * t1164;
t1041 = t750 * t1334;
t497 = -t623 * t881 + t624 * t884;
t1204 = t328 * t497;
t1012 = t561 * t1039 + t559 * t1041 + 0.2e1 * t1204;
t1208 = t299 * t497;
t1014 = t540 * t1039 + t537 * t1041 + 0.2e1 * t1208;
t1274 = m(7) * (t1012 + t1014);
t1197 = t391 * t497;
t1009 = t596 * t1039 + t594 * t1041 + 0.2e1 * t1197;
t1273 = m(7) * (t1009 + t1014);
t1103 = 0.2e1 * t380 * t473 + 0.2e1 * t381 * t472;
t1272 = m(7) * (t1354 * t522 + 0.2e1 * t1441 * t523 + t1103);
t1271 = m(7) * (t1009 + t1012);
t1270 = (t1441 * t381 - t339 * t444 + t380 * t482) * t1404;
t1269 = m(7) * (t1351 * t1441 + t1354 * t546 + t1103);
t1268 = m(7) * (-t1075 * t1354 + t1344 * t1441 + t1103);
t1043 = 0.4e1 * t750;
t1260 = m(7) * (0.4e1 * t1204 + (t559 * t881 + t561 * t884) * t1043);
t1040 = -0.2e1 * t1164;
t1096 = t473 * t1040 + t472 * t1041;
t1259 = m(7) * (t1338 * t537 + t1339 * t540 + t1096);
t1258 = m(7) * (t1338 * t559 + t1339 * t561 + t1096);
t1257 = m(7) * (t1338 * t594 + t1339 * t596 + t1096);
t1194 = t482 * t881;
t358 = t1332 * t1441 + 0.2e1 * t1194;
t1256 = m(7) * (t1333 * t1441 - 0.2e1 * t1194 + t358);
t1255 = m(7) * t358;
t1175 = t596 * t881;
t574 = -0.2e1 * t1175;
t1247 = m(7) * (t1332 * t594 + t574);
t1246 = m(7) * (t574 + 0.2e1 * t1175);
t1245 = t938 * t1045;
t1237 = m(5) * qJD(3);
t1236 = m(5) * qJD(4);
t1235 = m(6) * qJD(3);
t1234 = m(6) * qJD(4);
t1232 = m(7) * qJD(3);
t1231 = m(7) * qJD(4);
t1200 = t381 * t881;
t1201 = t380 * t884;
t231 = (t1201 / 0.2e1 - t1200 / 0.2e1 + t1384 * t750) * m(7);
t1048 = t231 * qJD(2);
t1202 = t376 * t884;
t943 = -t375 * t881 - t1202;
t198 = t461 * t856 + t857 * t943;
t258 = -t1375 * t881 + t640 * t788 + t642 * t789;
t259 = -t1376 * t881 - t641 * t788 - t643 * t789;
t52 = (-t258 * t881 + t259 * t884 + t421) * t857 + (t313 - t945) * t856;
t260 = t1375 * t884 + t640 * t790 - t642 * t791;
t261 = t1376 * t884 - t641 * t790 + t643 * t791;
t53 = (-t260 * t881 + t261 * t884 + t423) * t857 + (t314 - t944) * t856;
t71 = (-t275 * t881 + t276 * t884 + t461) * t857 + (t332 - t943) * t856;
t28 = t1270 / 0.4e1 + (t52 * t1307 + t53 * t1304 + t198 / 0.2e1) * t857 + (t1135 / 0.2e1 - t1121 / 0.2e1 + t71 / 0.2e1) * t856;
t1224 = t28 * qJD(6) - t1048;
t1211 = t275 * t884;
t1210 = t276 * t881;
t488 = -t696 - t752 + (-t733 - t795) * t881;
t1193 = t488 * t611;
t1125 = t883 * t771;
t976 = 0.2e1 * t1049;
t230 = (t1200 - t1201) * t1427 + t750 * t976 * t1326;
t925 = t976 / 0.2e1;
t984 = t1045 / 0.2e1;
t398 = (t593 * t881 - t595 * t884) * t984 - m(6) * t966 * t925;
t1108 = t398 * qJD(5) + t230 * qJD(6);
t1094 = 0.4e1 * t1386;
t1086 = 0.4e1 * t1387;
t1042 = 0.4e1 * t966;
t1025 = -t173 / 0.2e1 + t68 / 0.2e1;
t1024 = -t69 / 0.2e1 - t174 / 0.2e1;
t1023 = t1274 / 0.4e1 + t1113;
t1022 = t1273 / 0.4e1 + t1113;
t1021 = t1271 / 0.4e1 + t1113;
t326 = (-m(6) * t690 - m(7) * t560) * t884 - m(5) * t967 * t925 + (m(6) * t688 + m(7) * t558) * t881;
t1016 = t326 * qJD(4) + t1108;
t1011 = 0.2e1 * t1379;
t1007 = t482 * t1040 + t1041 * t1441 + t497 * t1357;
t1002 = -0.2e1 * t1049 * t1163 + 0.2e1 * t1182;
t951 = Icges(4,5) * t880 + Icges(4,6) * t883;
t769 = Icges(4,3) * t884 + t881 * t951;
t516 = t881 * t1125 + t773 * t1138 + t884 * t769;
t770 = -Icges(4,3) * t881 + t884 * t951;
t517 = -t881 * t1124 - t774 * t1138 - t884 * t770;
t995 = -t1148 / 0.2e1;
t994 = -t1148 / 0.4e1;
t993 = t1148 / 0.4e1;
t992 = -t1146 / 0.4e1;
t991 = t1146 / 0.2e1;
t990 = t1146 / 0.4e1;
t131 = t258 * t884 + t259 * t881;
t132 = t260 * t884 + t261 * t881;
t948 = Icges(6,5) * t857 - Icges(6,6) * t856;
t741 = t881 * t948;
t742 = t948 * t884;
t980 = (t132 - t876 * t742 + (t881 * t741 + t1409) * t884) * t1306 + (t131 + t877 * t741 + (-t884 * t742 - t1409) * t881) * t1304;
t973 = t1049 * t1243;
t970 = t1193 - t714 * t1036 / 0.2e1 - t969 * t1161;
t952 = Icges(4,5) * t883 - Icges(4,6) * t880;
t950 = Icges(5,5) * t860 - Icges(5,6) * t859;
t934 = t700 * t856 - t702 * t857;
t929 = -t880 * t773 - t1125;
t922 = 0.4e1 * t1328 * (t575 * t675 + t576 * t674) + 0.4e1 * t1330 * (t612 * t787 + t613 * t786);
t921 = 0.4e1 * t1378 * t1326 + (0.4e1 * t1196 - (t689 * t881 + t691 * t884) * t1042) * t1328;
t920 = t980 + t1380;
t919 = t980 + t1381;
t778 = t881 * t950;
t779 = t950 * t884;
t918 = (t877 * t778 + (-t884 * t779 - t1408) * t881) * t1304 + (-t876 * t779 + (t881 * t778 + t1408) * t884) * t1306 + t980;
t917 = t971 * t881;
t389 = 0.4e1 * t1388;
t911 = t389 * t1330 + t918;
t906 = t1426 * t994 + t232 * t992 + t91 * t990 + t1363;
t902 = t52 * t1304 + t53 * t1306 + t131 * t995 + t132 * t991 + (t1210 + t1211) * t1309 + (t375 * t884 - t376 * t881) * t1431 - t1113 + t1420 * t856;
t901 = t1418 * t990 + t1419 * t994 + t1405;
t72 = t1318 / 0.4e1;
t899 = t72 + t902;
t898 = -m(6) * (t575 * t884 + t576 * t881) - m(5) * (t612 * t884 + t613 * t881) - m(4) * (t666 * t884 + t667 * t881) - m(3) * ((rSges(3,3) * t884 + t1383) * t884 + (t1238 + (rSges(3,3) + qJ(2)) * t881) * t881);
t122 = t348 * t856 + (-t277 * t881 + t278 * t884) * t857;
t123 = t349 * t856 + (-t279 * t881 + t280 * t884) * t857;
t897 = -t1270 / 0.4e1 + t123 * t1306 + t122 * t1304 + t149 * t995 + t150 * t991 + t198 * t1432 + t52 * t1148 / 0.2e1 - t53 * t1146 / 0.2e1 + (t1135 + t71) * t1433 + (t1121 + t1206 + t1207) * t1309;
t896 = t909 + t1365;
t894 = t900 + t1364;
t892 = t1210 / 0.2e1 + t1211 / 0.2e1 + (t1066 * t856 + t1068 * t857 + t1372 * t884 - t1424 + t314) * t1306 + (-t1067 * t856 - t1069 * t857 - t1372 * t881 - t1422 + t313) * t1304 - t909;
t891 = -t1202 / 0.2e1 - t900 + t934 * t1305 + (t376 + t934) * t1304;
t890 = t1426 * t993 + t232 * t990 + t91 * t992 - t1363 + t901 + t979;
t889 = t1418 * t992 + t1419 * t993 - t1405 + t906 + t979;
t888 = t901 + t906 - t979;
t887 = t892 - t1365 + (t1062 * t859 + t1064 * t860 + t1373 * t884 - t1423) * t1306 + (-t1063 * t859 - t1065 * t860 - t1373 * t881 - t1421) * t1304;
t886 = -t1364 + t891 + (t1305 + t1304) * (t730 * t859 - t732 * t860);
t806 = t952 * t884;
t805 = t881 * t952;
t735 = t881 * t769;
t713 = t1416 * t881;
t655 = t850 + (t966 + t1244) * t884;
t653 = -t776 - t917;
t566 = t611 - t973;
t539 = t850 + (t1071 + t1244) * t884;
t536 = -t917 + t593;
t519 = -t881 * t770 + t1395;
t518 = t884 * t929 + t735;
t506 = t1146 * t750 - t624 * t856;
t505 = t1148 * t750 + t623 * t856;
t489 = t1245 / 0.4e1;
t456 = t1246 / 0.4e1;
t455 = t1247 / 0.4e1;
t454 = -t973 + t481;
t395 = t518 * t884 + t519 * t881;
t394 = t516 * t884 + t517 * t881;
t384 = (t1157 + (-t1072 * t879 - t1073 * t882) * t857) * t856;
t363 = -t973 + t378;
t357 = t1255 / 0.4e1;
t347 = 0.4e1 * t472 * t884 + 0.4e1 * t473 * t881;
t333 = 0.4e1 * t1193 - 0.4e1 * (t714 * t881 + t884 * t969) * t967;
t290 = -0.4e1 * t472 * t624 + 0.4e1 * t473 * t623;
t274 = -0.4e1 * t1075 * t473 + 0.4e1 * t472 * t571;
t264 = 0.4e1 * t472 * t547 + 0.4e1 * t473 * t546;
t263 = 0.4e1 * t1199 - (t654 * t881 + t656 * t884) * t1042;
t262 = 0.4e1 * t472 * t523 + 0.4e1 * t473 * t522;
t257 = t1358 * t481 + t1086;
t252 = 0.4e1 * t383 * t538 + 0.4e1 * t988;
t243 = t455 + t456 - t1095;
t242 = t456 - t1247 / 0.4e1 + t1095;
t241 = t455 - t1246 / 0.4e1 + t1095;
t240 = 0.4e1 * t1197 + (t594 * t881 + t596 * t884) * t1043;
t239 = t1256 / 0.4e1;
t227 = t231 * qJD(6);
t214 = t1257 / 0.4e1;
t211 = t290 * t1326 - t1446;
t210 = t1258 / 0.4e1;
t209 = t876 * t770 + (t517 - t735 + (t770 - t929) * t884) * t884;
t208 = (-t518 + t735 + t517) * t881 + (t519 - t1395 + (t770 + t929) * t881 + t516) * t884;
t206 = t1259 / 0.4e1;
t189 = 0.4e1 * t1208 + (t537 * t881 + t540 * t884) * t1043;
t188 = t1326 * t347 - t898;
t144 = -0.4e1 * t1189 + 0.4e1 * t1191 + 0.4e1 * t1209;
t140 = t1097 + t1098 - t1003;
t139 = t1003 + t1413;
t138 = t1003 - t1413;
t137 = t1361 * t378 + t1094;
t136 = 0.4e1 * t299 * t405 + 0.4e1 * t989;
t135 = t357 + t239 - t1245 / 0.4e1;
t134 = t489 + t357 - t1256 / 0.4e1;
t133 = t489 + t239 - t1255 / 0.4e1;
t125 = t274 * t1326 + t1282 / 0.4e1 + t900;
t118 = t1268 / 0.4e1;
t113 = t1269 / 0.4e1;
t108 = t1272 / 0.4e1;
t107 = t1004 + t1005 - t978;
t106 = t978 + t1415;
t105 = t978 - t1415;
t101 = t1326 * t264 + t894 + t922;
t90 = t894 + t262 * t1326 + (-t842 / 0.2e1 + t956 / 0.2e1) * t880 - t1410;
t73 = -t1318 / 0.4e1;
t60 = t1320 / 0.4e1;
t58 = t1321 / 0.4e1;
t54 = t1323 / 0.4e1;
t46 = t1324 / 0.4e1;
t45 = t1326 * t240 + t1113;
t44 = t1260 / 0.4e1 + t1113;
t43 = t1326 * t189 + t1113;
t42 = t980 - t1381;
t41 = t42 * qJD(5);
t40 = t921 + t980;
t39 = t1326 * t144 + t1328 * t263 + t980;
t38 = t1326 * t137 + t1328 * t257 + t911;
t37 = t1326 * t136 + t1328 * t252 + t1330 * t333 + t918;
t36 = (t1024 * t881 + t1025 * t884) * t857;
t35 = t72 - t1320 / 0.4e1 + t1021;
t34 = t60 + t73 + t1021;
t33 = t72 - t1321 / 0.4e1 + t1022;
t32 = t58 + t73 + t1022;
t31 = t54 - t1324 / 0.4e1 + t1023;
t30 = t46 - t1323 / 0.4e1 + t1023;
t26 = t920 + t1435;
t25 = t920 - t1435;
t24 = t980 - t1380 + t1401 + t1402;
t22 = t909 + t1412;
t21 = t909 - t1412;
t20 = t909 + t1411;
t19 = t909 - t1411;
t18 = t896 + (t208 / 0.2e1 - t394 / 0.2e1) * t881 + (t395 / 0.2e1 + t209 / 0.2e1) * t884;
t17 = t892 + t1109 + t1110;
t16 = t892 + t1111 + t1112;
t15 = t60 + t899 - t1271 / 0.4e1;
t14 = t896 - t1414;
t13 = t896 + t1414;
t12 = t899 + t58 - t1273 / 0.4e1;
t11 = t46 + t54 - t1274 / 0.4e1 + t902;
t10 = t890 + t214 + t118;
t9 = t888 - t1257 / 0.4e1 + t118;
t8 = t889 + t214 - t1268 / 0.4e1;
t7 = t890 + t210 + t113;
t6 = t888 + t113 - t1258 / 0.4e1;
t5 = t889 + t210 - t1269 / 0.4e1;
t4 = t887 + t1019 + t1020;
t3 = t890 + t206 + t108;
t2 = t888 - t1259 / 0.4e1 + t108;
t1 = t889 + t206 - t1272 / 0.4e1;
t23 = [t188 * qJD(2) + t90 * qJD(3) + t101 * qJD(4) + t125 * qJD(5) + t211 * qJD(6), qJD(1) * t188 + qJD(3) * t105 + qJD(4) * t138 + qJD(5) * t241 + qJD(6) * t134, t90 * qJD(1) + t105 * qJD(2) + t4 * qJD(4) + t16 * qJD(5) + t3 * qJD(6) + (t472 * t536 + t473 * t539 - t522 * t540 + t523 * t537) * t1232 + (t575 * t653 + t576 * t655 + t654 * t658 - t656 * t657) * t1235 + (-t612 * t713 + t613 * t715 + (-t708 + t714) * t969) * t1237 + (t887 + t208 * t1307 + (-t1058 * t880 - t1060 * t883) * t1304 + (t1057 * t880 + t1059 * t883 + t394) * t1306 + (t395 + t209) * t1305 - t1384 * t951 + (t926 * t845 - (t666 * t881 - t667 * t884) * t968) * m(4)) * qJD(3), t101 * qJD(1) + t138 * qJD(2) + t4 * qJD(3) + t887 * qJD(4) + t17 * qJD(5) + t7 * qJD(6) + (-t546 * t561 + t547 * t559 + t386 / 0.2e1 + t387 / 0.2e1) * t1231 + (-t674 * t691 + t675 * t689 + t500 / 0.2e1 + t501 / 0.2e1) * t1234 + (-t553 / 0.2e1 + t554 / 0.2e1 + t927 * t822) * t1236, t125 * qJD(1) + t241 * qJD(2) + t16 * qJD(3) + t17 * qJD(4) + t10 * qJD(6) + (t892 + (t1075 * t596 + t571 * t594 + t399 / 0.2e1 + t400 / 0.2e1) * m(7) + (t514 / 0.2e1 + t515 / 0.2e1 + t930 * t802) * m(6)) * qJD(5), t211 * qJD(1) + t134 * qJD(2) + t3 * qJD(3) + t7 * qJD(4) + t10 * qJD(5) + (t384 + (-t1441 * t624 + t472 * t506 + t473 * t505 + t482 * t623) * m(7) + ((t301 / 0.2e1 + t349 / 0.2e1 - t1025) * t884 + (-t300 / 0.2e1 - t348 / 0.2e1 - t1024) * t881) * t857) * qJD(6); t898 * qJD(1) + t106 * qJD(3) + t139 * qJD(4) + t242 * qJD(5) + t133 * qJD(6) + t1440 * t347, 0, t106 * qJD(1) + (-m(4) * t968 * t925 + (-m(5) * t715 - m(6) * t655 - m(7) * t539) * t884 + (-m(5) * t713 + m(6) * t653 + m(7) * t536) * t881) * qJD(3) + t1016, qJD(1) * t139 + qJD(3) * t326 + t1016, qJD(1) * t242 + t1044 * t398 + t1108, t133 * qJD(1) + (-t505 * t884 + t506 * t881) * qJD(6) * t984 + t982 * t230; (t886 + (-t956 + t842) * t880 / 0.2e1 + t1410) * qJD(1) + t107 * qJD(2) + t18 * qJD(3) + t14 * qJD(4) + t19 * qJD(5) + t1 * qJD(6) + t262 * t1440, qJD(1) * t107 + t227, t18 * qJD(1) + ((0.4e1 * (-t884 * t907 + (-t884 * rSges(4,3) - t881 * t968) * t881) * (-t811 * t881 - t812 * t884) - t845 * t968 * t975) * t1331 + (-t876 * t806 + (t881 * t805 + t1407) * t884) * t1306 + (t877 * t805 + (-t884 * t806 - t1407) * t881) * t1304 + (t299 * t363 + t536 * t537 - t539 * t540) * m(7) + t918 + m(5) * (t488 * t566 - t713 * t714 - t715 * t969) + m(6) * (t383 * t454 + t653 * t654 - t655 * t656)) * qJD(3) + t37 * qJD(4) + t39 * qJD(5) + t43 * qJD(6), t14 * qJD(1) + t37 * qJD(3) + t918 * qJD(4) + t24 * qJD(5) + t30 * qJD(6) + (-t389 / 0.4e1 + t970 + t1388) * t1236 + (-t137 / 0.4e1 + (t299 + t328) * t405 + t989 + t1386) * t1231 + (-t257 / 0.4e1 + (t383 + t425) * t538 + t988 + t1387) * t1234, t19 * qJD(1) + t39 * qJD(3) + t24 * qJD(4) + ((t1011 + t1015) * t1327 + (t1002 + t1008) * t1329 + t919) * qJD(5) + t32 * qJD(6), t1 * qJD(1) + t1048 + t43 * qJD(3) + t30 * qJD(4) + t32 * qJD(5) + (t897 + (t1352 * t505 + t1353 * t506 + t1355 * t299 + t1007) * t1327) * qJD(6); (t886 - t922) * qJD(1) + t140 * qJD(2) + t13 * qJD(3) + t896 * qJD(4) + t21 * qJD(5) + t5 * qJD(6) + t264 * t1440, qJD(1) * t140 + t227, t13 * qJD(1) + t918 * qJD(3) + t38 * qJD(4) + t25 * qJD(5) + t31 * qJD(6) + (-t333 / 0.4e1 + t581 * t566 + (-t713 * t881 - t715 * t884) * t822 + t970) * t1237 + (-t136 / 0.4e1 + t299 * t378 + t328 * t363 + t536 * t559 - t539 * t561 + t989) * t1232 + (-t252 / 0.4e1 + t383 * t481 + t425 * t454 + t653 * t689 - t655 * t691 + t988) * t1235, t896 * qJD(1) + t38 * qJD(3) + ((t1361 * t405 + t1094) * t1326 + (t1358 * t538 + t1086) * t1328 + t911) * qJD(4) + t40 * qJD(5) + t44 * qJD(6), t21 * qJD(1) + t25 * qJD(3) + t40 * qJD(4) + ((t1011 + t1013) * t1327 + (t1002 + t1006) * t1329 + t919) * qJD(5) + t34 * qJD(6), t5 * qJD(1) + t1048 + t31 * qJD(3) + t44 * qJD(4) + t34 * qJD(5) + (t897 + (t1347 * t505 + t1349 * t506 + t1355 * t328 + t1007) * t1327) * qJD(6); (t891 - t1282 / 0.4e1) * qJD(1) + t243 * qJD(2) + t20 * qJD(3) + t22 * qJD(4) + t909 * qJD(5) + t8 * qJD(6) + t274 * t1440, qJD(1) * t243 + t227, t20 * qJD(1) + t980 * qJD(3) + t26 * qJD(4) + t41 + t33 * qJD(6) + (t363 * t391 + t536 * t594 - t539 * t596 + t238 / 0.2e1 + t451 / 0.2e1 + t452 / 0.2e1 - t144 / 0.4e1) * t1232 + (t569 * t454 + t331 / 0.2e1 + t577 / 0.2e1 + t578 / 0.2e1 - t263 / 0.4e1 + (t653 * t881 - t655 * t884) * t802) * t1235, t22 * qJD(1) + t26 * qJD(3) + ((t281 + t962) * t1327 + (t434 + t961) * t1329 - t921 + t980) * qJD(4) + t41 + t35 * qJD(6), qJD(1) * t909 + qJD(6) * t45 + t1044 * t42 + t41, t8 * qJD(1) + t1048 + t33 * qJD(3) + t35 * qJD(4) + t45 * qJD(5) + (t897 + (t1340 * t505 + t1342 * t506 + t1355 * t391 + t1007) * t1327) * qJD(6); qJD(1) * t1446 + t135 * qJD(2) + t2 * qJD(3) + t6 * qJD(4) + t9 * qJD(5) + t36 * qJD(6) + t290 * t1440, qJD(1) * t135 - t231 * t982, t2 * qJD(1) + t902 * qJD(3) + t11 * qJD(4) + t12 * qJD(5) + (-t363 * t444 + t482 * t539 + t1441 * t536 + t203 / 0.2e1 + t321 / 0.2e1 + t322 / 0.2e1 - t189 / 0.4e1) * t1232 + t1224, t6 * qJD(1) + t11 * qJD(3) + ((t302 + t963) * t1327 - t1260 / 0.4e1 + t902) * qJD(4) + t15 * qJD(5) + t1224, t9 * qJD(1) + t12 * qJD(3) + t15 * qJD(4) + ((t74 / 0.2e1 - t240 / 0.4e1) * m(7) + t902) * qJD(5) + t1224, t36 * qJD(1) + (t384 * t1309 + (t1441 * t506 - t444 * t474 + t482 * t505) * m(7) + (t122 * t1307 + t123 * t1304 + (-t300 * t881 + t301 * t884) * t1309) * t857) * qJD(6) + t982 * t28;];
Cq  = t23;