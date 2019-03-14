% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:29:43
% EndTime: 2019-03-09 08:30:55
% DurationCPUTime: 61.78s
% Computational Cost: add. (82542->1147), mult. (105540->1535), div. (0->0), fcn. (111925->8), ass. (0->647)
t780 = qJ(2) + pkin(9);
t769 = cos(t780);
t785 = sin(qJ(1));
t1055 = t769 * t785;
t788 = cos(qJ(1));
t1054 = t769 * t788;
t786 = cos(qJ(5));
t1047 = t786 * t788;
t783 = sin(qJ(5));
t1050 = t785 * t783;
t768 = sin(t780);
t682 = t768 * t1047 - t1050;
t1049 = t785 * t786;
t1053 = t783 * t788;
t683 = t1053 * t768 + t1049;
t482 = Icges(7,5) * t683 + Icges(7,6) * t682 + Icges(7,3) * t1054;
t485 = Icges(6,5) * t683 + Icges(6,6) * t682 + Icges(6,3) * t1054;
t1088 = Icges(7,4) * t683;
t488 = Icges(7,2) * t682 + Icges(7,6) * t1054 + t1088;
t1091 = Icges(6,4) * t683;
t491 = Icges(6,2) * t682 + Icges(6,6) * t1054 + t1091;
t664 = Icges(7,4) * t682;
t494 = Icges(7,1) * t683 + Icges(7,5) * t1054 + t664;
t667 = Icges(6,4) * t682;
t497 = Icges(6,1) * t683 + Icges(6,5) * t1054 + t667;
t684 = t1049 * t768 + t1053;
t949 = t768 * t1050;
t685 = -t949 + t1047;
t1367 = (-t494 - t497) * t685 + (t488 + t491) * t684 + (t482 + t485) * t1055;
t484 = -Icges(7,5) * t685 + Icges(7,6) * t684 + Icges(7,3) * t1055;
t487 = -Icges(6,5) * t685 + Icges(6,6) * t684 + Icges(6,3) * t1055;
t666 = Icges(7,4) * t685;
t490 = Icges(7,2) * t684 + Icges(7,6) * t1055 - t666;
t669 = Icges(6,4) * t685;
t493 = Icges(6,2) * t684 + Icges(6,6) * t1055 - t669;
t665 = Icges(7,4) * t684;
t495 = Icges(7,1) * t685 - Icges(7,5) * t1055 - t665;
t668 = Icges(6,4) * t684;
t498 = Icges(6,1) * t685 - Icges(6,5) * t1055 - t668;
t1366 = (t495 + t498) * t685 + (t490 + t493) * t684 + (t484 + t487) * t1055;
t1223 = m(7) / 0.2e1;
t1319 = -m(6) / 0.2e1;
t1336 = rSges(7,1) + pkin(5);
t1003 = t685 * rSges(7,2) + t1336 * t684;
t462 = -t683 * rSges(7,2) + t1336 * t682;
t340 = t1003 * t785 + t462 * t788;
t536 = rSges(6,1) * t682 - rSges(6,2) * t683;
t538 = rSges(6,1) * t684 + rSges(6,2) * t685;
t418 = t536 * t788 + t785 * t538;
t1025 = -t1223 * t340 + t418 * t1319;
t890 = rSges(6,1) * t783 + rSges(6,2) * t786;
t1240 = -rSges(6,3) * t768 + t769 * t890;
t1259 = t685 * rSges(6,1) - t684 * rSges(6,2);
t501 = -rSges(6,3) * t1055 + t1259;
t1339 = -t1240 * t1055 + t501 * t768;
t892 = t683 * rSges(6,1) + t682 * rSges(6,2);
t817 = rSges(6,3) * t1054 + t892;
t416 = -t1054 * t1240 - t768 * t817;
t1262 = -t1339 * t788 + t416 * t785;
t1318 = -m(7) / 0.2e1;
t1116 = pkin(5) * t783;
t1358 = rSges(7,3) + qJ(6);
t1100 = t786 * rSges(7,2);
t887 = rSges(7,1) * t783 + t1100;
t1328 = (t1116 + t887) * t769 - t1358 * t768;
t1111 = -qJ(6) - pkin(8);
t1115 = pkin(5) * t786;
t945 = pkin(4) + t1115;
t1355 = t685 * rSges(7,1) - t684 * rSges(7,2) + t1111 * t1055 + t788 * t945;
t752 = pkin(8) * t1055;
t716 = pkin(4) * t788 - t752;
t1343 = -rSges(7,3) * t1055 - pkin(5) * t949 + t1355 - t716;
t1359 = -t1328 * t1055 + t1343 * t768;
t924 = t788 * t1111;
t845 = -t683 * rSges(7,1) - t682 * rSges(7,2) + t769 * t924;
t959 = pkin(8) * t1054;
t1257 = rSges(7,3) * t1054 + pkin(5) * t683 - t845 - t959;
t1258 = t1257 * t768;
t301 = -t1054 * t1328 - t1258;
t869 = -t1359 * t788 + t301 * t785;
t1032 = t1262 * t1319 + t1318 * t869;
t24 = t1032 - t1025;
t1365 = qJD(1) * t24;
t875 = Icges(7,5) * t783 + Icges(7,6) * t786;
t807 = -Icges(7,3) * t768 + t769 * t875;
t1087 = Icges(7,4) * t783;
t879 = Icges(7,2) * t786 + t1087;
t809 = -Icges(7,6) * t768 + t769 * t879;
t1086 = Icges(7,4) * t786;
t882 = Icges(7,1) * t783 + t1086;
t811 = -Icges(7,5) * t768 + t769 * t882;
t351 = -t1055 * t807 - t684 * t809 + t685 * t811;
t876 = Icges(6,5) * t783 + Icges(6,6) * t786;
t808 = -Icges(6,3) * t768 + t769 * t876;
t1090 = Icges(6,4) * t783;
t880 = Icges(6,2) * t786 + t1090;
t810 = -Icges(6,6) * t768 + t769 * t880;
t1089 = Icges(6,4) * t786;
t883 = Icges(6,1) * t783 + t1089;
t812 = -Icges(6,5) * t768 + t769 * t883;
t352 = -t1055 * t808 - t684 * t810 + t685 * t812;
t1303 = t351 + t352;
t1364 = t1366 * t785 + t1367 * t788;
t1199 = t768 / 0.2e1;
t1197 = t769 / 0.2e1;
t1279 = -t785 / 0.2e1;
t1278 = t788 / 0.2e1;
t1075 = qJ(4) * t768;
t1057 = t768 * t785;
t919 = rSges(5,1) * t788 - rSges(5,3) * t1057;
t1112 = -qJ(3) - pkin(7);
t787 = cos(qJ(2));
t1119 = pkin(2) * t787;
t765 = pkin(1) + t1119;
t973 = -t788 * t1112 - t785 * t765;
t464 = (-t1075 + (rSges(5,2) - pkin(3)) * t769) * t785 + t919 + t973;
t1099 = rSges(5,3) + qJ(4);
t1117 = pkin(3) * t769;
t763 = t785 * t1112;
t918 = t785 * rSges(5,1) - rSges(5,2) * t1054;
t465 = -t763 + (t1099 * t768 + t1117 + t765) * t788 + t918;
t1016 = t464 * t1054 + t465 * t1055;
t1291 = (-t1075 + (-rSges(6,3) - pkin(3)) * t769) * t785 + t716 + t973 + t1259;
t1113 = t785 * pkin(4);
t1189 = rSges(6,3) + pkin(8);
t961 = pkin(3) + t1189;
t392 = t1113 - t763 + (t769 * t961 + t1075 + t765) * t788 + t892;
t1022 = t1054 * t1291 + t392 * t1055;
t1188 = rSges(7,3) + pkin(3);
t1236 = -t1188 * t769 - t768 * (qJ(4) + t1116);
t1292 = t1236 * t785 + t1355 + t973;
t363 = -t763 + t785 * t945 + (-t1236 + t765) * t788 - t845;
t1024 = t1054 * t1292 + t363 * t1055;
t1225 = m(6) / 0.2e1;
t1226 = m(5) / 0.2e1;
t784 = sin(qJ(2));
t1114 = t784 * pkin(2);
t751 = pkin(3) * t1057;
t433 = t751 + (t1114 + t1189 * t768 + (-qJ(4) - t890) * t769) * t785;
t737 = qJ(4) * t1054;
t950 = t769 * t1047;
t951 = t769 * t1053;
t977 = rSges(6,1) * t951 + rSges(6,2) * t950;
t434 = t737 + (-t768 * t961 - t1114) * t788 + t977;
t1104 = rSges(5,2) * t768;
t512 = t751 + (-t1099 * t769 - t1104 + t1114) * t785;
t1251 = -pkin(3) * t768 - t1114;
t1056 = t768 * t788;
t974 = rSges(5,2) * t1056 + rSges(5,3) * t1054;
t513 = t1251 * t788 + t737 + t974;
t422 = t751 + (t1114 + (rSges(7,3) - t1111) * t768 + (-t1336 * t783 - qJ(4) - t1100) * t769) * t785;
t1299 = rSges(7,2) * t950 + t1336 * t951 + t768 * t924;
t423 = t737 + (-t1188 * t768 - t1114) * t788 + t1299;
t864 = t422 * t788 + t423 * t785;
t946 = (t768 * t864 + t1024) * t1223 + ((t433 * t788 + t434 * t785) * t768 + t1022) * t1225 + ((t512 * t788 + t513 * t785) * t768 + t1016) * t1226;
t923 = -qJ(4) * t769 - t1251;
t850 = t768 * pkin(8) + t923;
t818 = t850 - t1328;
t409 = t818 * t785;
t411 = t818 * t788;
t825 = -t1240 + t850;
t454 = t825 * t785;
t456 = t825 * t788;
t886 = rSges(5,3) * t769 + t1104;
t895 = -t886 + t923;
t516 = t895 * t785;
t518 = t895 * t788;
t947 = (-t409 * t1056 + t1057 * t411 + t1024) * t1223 + (-t454 * t1056 + t1057 * t456 + t1022) * t1225 + (-t516 * t1056 + t1057 * t518 + t1016) * t1226;
t13 = t947 - t946;
t1363 = t13 * qJD(1);
t1362 = Icges(5,4) - Icges(4,5);
t1361 = Icges(5,5) - Icges(4,6);
t1360 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t937 = t1364 * t1197 + t1199 * t1303;
t523 = Icges(7,5) * t682 - Icges(7,6) * t683;
t525 = Icges(6,5) * t682 - Icges(6,6) * t683;
t1357 = t523 + t525;
t524 = Icges(7,5) * t684 + Icges(7,6) * t685;
t526 = Icges(6,5) * t684 + Icges(6,6) * t685;
t1356 = t524 + t526;
t1001 = t1328 * t785;
t1007 = Icges(6,2) * t685 - t498 + t668;
t1009 = Icges(7,2) * t685 - t495 + t665;
t1353 = t1007 + t1009;
t1008 = -Icges(6,2) * t683 + t497 + t667;
t1010 = -Icges(7,2) * t683 + t494 + t664;
t1352 = t1008 + t1010;
t1011 = -Icges(6,1) * t684 + t493 - t669;
t1013 = -Icges(7,1) * t684 + t490 - t666;
t1351 = t1011 + t1013;
t1012 = -Icges(6,1) * t682 + t1091 + t491;
t1014 = -Icges(7,1) * t682 + t1088 + t488;
t1350 = t1012 + t1014;
t704 = Icges(4,5) * t769 - Icges(4,6) * t768;
t705 = -Icges(5,4) * t769 + Icges(5,5) * t768;
t728 = Icges(3,5) * t787 - Icges(3,6) * t784;
t1347 = (t704 + t705 + t728) * t788 + t1360 * t785;
t1048 = t785 * t787;
t1052 = t784 * t785;
t1327 = -Icges(3,5) * t1048 + Icges(3,6) * t1052 + t1055 * t1362 - t1361 * t1057 + t1360 * t788;
t348 = -t1054 * t808 - t682 * t810 - t683 * t812;
t278 = t485 * t1054 + t682 * t491 + t683 * t497;
t279 = t487 * t1054 + t682 * t493 - t683 * t498;
t872 = t278 * t788 + t785 * t279;
t1340 = t348 * t768 + t769 * t872;
t347 = -t1054 * t807 - t682 * t809 - t683 * t811;
t276 = t482 * t1054 + t682 * t488 + t683 * t494;
t277 = t484 * t1054 + t682 * t490 - t683 * t495;
t873 = t276 * t788 + t785 * t277;
t1341 = t347 * t768 + t769 * t873;
t938 = t1341 / 0.2e1 + t1340 / 0.2e1;
t1067 = t487 * t768;
t1320 = t493 * t786 - t498 * t783;
t307 = t1320 * t769 - t1067;
t1069 = t484 * t768;
t1321 = t490 * t786 - t495 * t783;
t304 = t1321 * t769 - t1069;
t1227 = m(4) / 0.2e1;
t712 = rSges(4,1) * t768 + rSges(4,2) * t769;
t820 = t712 + t1114;
t1267 = t820 * t788;
t1268 = t820 * t785;
t1284 = t1267 * t788 + t1268 * t785;
t899 = (t785 * t422 - t423 * t788) * t1223 + (t785 * t433 - t434 * t788) * t1225 + (t785 * t512 - t513 * t788) * t1226 + t1284 * t1227;
t447 = t785 * t454;
t1252 = t456 * t788 + t447;
t507 = t785 * t516;
t396 = t785 * t409;
t865 = t411 * t788 + t396;
t900 = t865 * t1318 + t1252 * t1319 + (-t518 * t788 - t507) * t1226 - m(4) * t1284 / 0.2e1;
t33 = t900 - t899;
t1346 = qJD(1) * t33;
t1070 = t482 * t768;
t862 = -t488 * t786 - t494 * t783;
t303 = t769 * t862 + t1070;
t1068 = t485 * t768;
t860 = -t491 * t786 - t497 * t783;
t306 = t769 * t860 + t1068;
t1345 = t303 + t306;
t1344 = -t304 - t307;
t1092 = Icges(4,4) * t768;
t709 = Icges(4,1) * t769 - t1092;
t615 = Icges(4,5) * t785 + t709 * t788;
t1080 = Icges(5,6) * t769;
t701 = Icges(5,3) * t768 - t1080;
t616 = Icges(5,5) * t785 + t701 * t788;
t1093 = Icges(3,4) * t784;
t732 = Icges(3,1) * t787 - t1093;
t645 = Icges(3,5) * t785 + t732 * t788;
t1342 = -t645 * t1048 - t615 * t1055 - t616 * t1057;
t1338 = t278 * t785 - t279 * t788;
t1337 = t276 * t785 - t277 * t788;
t1335 = -rSges(7,3) + pkin(8);
t1334 = t1054 * t1357 - t1350 * t683 + t1352 * t682;
t1333 = t1054 * t1356 - t1351 * t683 + t1353 * t682;
t1332 = t1055 * t1357 + t1350 * t685 + t1352 * t684;
t1331 = t1055 * t1356 + t1351 * t685 + t1353 * t684;
t650 = (-Icges(7,5) * t786 + Icges(7,6) * t783) * t769;
t654 = (Icges(7,2) * t783 - t1086) * t769;
t992 = -t811 + t654;
t660 = (-Icges(7,1) * t786 + t1087) * t769;
t994 = -t809 - t660;
t268 = t1054 * t650 + t682 * t992 - t683 * t994;
t651 = (-Icges(6,5) * t786 + Icges(6,6) * t783) * t769;
t655 = (Icges(6,2) * t783 - t1089) * t769;
t991 = -t812 + t655;
t661 = (-Icges(6,1) * t786 + t1090) * t769;
t993 = -t810 - t661;
t269 = t1054 * t651 + t682 * t991 - t683 * t993;
t1330 = t268 + t269;
t270 = t1055 * t650 + t684 * t992 + t685 * t994;
t271 = t1055 * t651 + t684 * t991 + t685 * t993;
t1329 = t270 + t271;
t781 = t785 ^ 2;
t782 = t788 ^ 2;
t970 = t781 + t782;
t559 = t1240 * t785;
t1326 = -Icges(3,5) * t784 - Icges(3,6) * t787 + t1361 * t769 + t1362 * t768;
t1325 = t1347 * t788 + t1342;
t1046 = t787 * t788;
t744 = Icges(4,4) * t1057;
t614 = Icges(4,1) * t1055 - Icges(4,5) * t788 - t744;
t617 = Icges(5,5) * t788 + Icges(5,6) * t1055 - Icges(5,3) * t1057;
t759 = Icges(3,4) * t1052;
t644 = Icges(3,1) * t1048 - Icges(3,5) * t788 - t759;
t1324 = -t644 * t1046 - t614 * t1054 + t617 * t1056 + t1327 * t785;
t1285 = t645 * t1046 + t615 * t1054 + t616 * t1056 + t1347 * t785;
t612 = Icges(4,4) * t1055 - Icges(4,2) * t1057 - Icges(4,6) * t788;
t739 = Icges(5,6) * t1057;
t619 = Icges(5,4) * t788 + Icges(5,2) * t1055 - t739;
t642 = Icges(3,4) * t1048 - Icges(3,2) * t1052 - Icges(3,6) * t788;
t1323 = t612 * t768 - t619 * t769 + t642 * t784;
t1263 = -t1291 * t788 - t392 * t785;
t1264 = -t1292 * t788 - t363 * t785;
t1322 = t1003 * t788 - t462 * t785;
t1280 = -t769 / 0.2e1;
t549 = t809 * t788;
t553 = t811 * t788;
t833 = t807 * t788 - t862;
t212 = (-t549 * t786 - t553 * t783 + t482) * t769 + t833 * t768;
t551 = t810 * t788;
t555 = t812 * t788;
t831 = t808 * t788 - t860;
t214 = (-t551 * t786 - t555 * t783 + t485) * t769 + t831 * t768;
t1311 = t212 + t214;
t548 = t809 * t785;
t552 = t811 * t785;
t832 = t807 * t785 + t1321;
t213 = (-t548 * t786 - t552 * t783 + t484) * t769 + t832 * t768;
t550 = t810 * t785;
t554 = t812 * t785;
t830 = t808 * t785 + t1320;
t215 = (-t550 * t786 - t554 * t783 + t487) * t769 + t830 * t768;
t1310 = t213 + t215;
t224 = t523 * t768 + (-t1010 * t786 + t1014 * t783) * t769;
t226 = t525 * t768 + (-t1008 * t786 + t1012 * t783) * t769;
t1309 = t224 + t226;
t225 = t524 * t768 + (-t1009 * t786 + t1013 * t783) * t769;
t227 = t526 * t768 + (-t1007 * t786 + t1011 * t783) * t769;
t1308 = t225 + t227;
t586 = Icges(7,6) * t769 + t768 * t879;
t590 = Icges(7,5) * t769 + t768 * t882;
t1066 = t807 * t768;
t582 = Icges(7,3) * t769 + t768 * t875;
t856 = t783 * t811 + t786 * t809;
t829 = t582 - t856;
t801 = t769 * t829 + t1066;
t240 = t586 * t684 - t590 * t685 + t785 * t801;
t588 = Icges(6,6) * t769 + t768 * t880;
t592 = Icges(6,5) * t769 + t768 * t883;
t1065 = t808 * t768;
t584 = Icges(6,3) * t769 + t768 * t876;
t855 = t783 * t812 + t786 * t810;
t828 = t584 - t855;
t800 = t769 * t828 + t1065;
t241 = t588 * t684 - t592 * t685 + t785 * t800;
t1307 = t240 + t241;
t242 = t682 * t586 + t683 * t590 + t788 * t801;
t243 = t682 * t588 + t683 * t592 + t788 * t800;
t1306 = t242 + t243;
t256 = (-t786 * t586 - t783 * t590 - t807) * t769 + t829 * t768;
t257 = (-t786 * t588 - t783 * t592 - t808) * t769 + t828 * t768;
t1305 = t256 + t257;
t1304 = t347 + t348;
t375 = t769 * t856 - t1066;
t376 = t769 * t855 - t1065;
t1302 = t375 + t376;
t1301 = t1344 * t785 + t1345 * t788;
t729 = Icges(3,2) * t787 + t1093;
t1296 = (t729 / 0.2e1 - t732 / 0.2e1) * t784;
t1107 = m(7) * qJD(6);
t691 = t970 * t769;
t541 = (t691 - t769) * t768;
t1293 = t541 * t1107;
t290 = t301 * t1056;
t1033 = (-t1057 * t1359 - t290) * t1223 + (-t416 * t1056 - t1057 * t1339) * t1225;
t823 = t301 * t788;
t1098 = (t768 * t823 - t290) * t1223;
t1289 = t1033 - t1098;
t1020 = t1054 * t1339 - t416 * t1055;
t1030 = t1054 * t1359 - t301 * t1055;
t1265 = t970 * t768;
t960 = t768 * t1116;
t251 = (-t1343 * t788 + (-pkin(5) * t1049 + (t1335 * t769 - t960) * t788 + t845) * t785) * t769;
t1242 = t501 * t788 + t785 * t817;
t365 = t1242 * t769;
t1096 = (-t1265 * t365 + t1020) * t1225 + (t1265 * t251 + t1030) * t1223;
t1261 = t1056 * t1335 + t1299;
t155 = t1001 * t1054 - t1055 * t1261 + t1056 * t1343 + t1258 * t785;
t561 = -rSges(6,3) * t1056 + t977;
t260 = (t559 * t788 - t561 * t785) * t769 + t1242 * t768;
t1102 = rSges(6,3) * t769;
t597 = t768 * t890 + t1102;
t309 = (t597 * t785 + t501) * t769;
t310 = ((-t597 + t1102) * t788 + t892) * t769 + (-t1240 * t788 + t561) * t768;
t990 = t1358 * t769 + t768 * t887 + t960;
t198 = (t785 * t990 + t1343) * t769;
t199 = -t1054 * t990 - t1056 * t1328 + t1257 * t769 + t1261 * t768;
t834 = t198 * t788 + t199 * t785 + t251;
t1110 = (-t155 * t769 + t768 * t834 + t1030) * t1223 + (-t260 * t769 + (t309 * t788 + t310 * t785 - t365) * t768 + t1020) * t1225;
t1288 = t1096 - t1110;
t1051 = t784 * t788;
t1082 = Icges(4,2) * t768;
t762 = Icges(4,4) * t769;
t613 = Icges(4,6) * t785 + (t762 - t1082) * t788;
t740 = Icges(5,6) * t1056;
t618 = Icges(5,4) * t785 - Icges(5,2) * t1054 + t740;
t776 = Icges(3,4) * t787;
t730 = -Icges(3,2) * t784 + t776;
t643 = Icges(3,6) * t785 + t730 * t788;
t1287 = -t1051 * t643 - t1054 * t618 - t1056 * t613 + t1285;
t1286 = t613 * t768 + t618 * t769 + t643 * t784 + t1327;
t1283 = -t1051 * t642 - t1052 * t643 + t1054 * t619 - t1055 * t618 - t1056 * t612 - t1057 * t613 - t1324 - t1325;
t1282 = -0.2e1 * t1265;
t1281 = -t768 / 0.2e1;
t1196 = t785 / 0.2e1;
t1193 = -t788 / 0.2e1;
t1275 = m(7) * t1322;
t1274 = (t1333 * t785 + t1334 * t788) * t769 + t1330 * t768;
t1273 = (t1331 * t785 + t1332 * t788) * t769 + t1329 * t768;
t1222 = m(7) / 0.4e1;
t1224 = m(6) / 0.4e1;
t903 = m(5) / 0.4e1 + t1224 + t1222;
t1266 = t903 * t541;
t647 = Icges(5,3) * t1054 + t740;
t706 = Icges(4,2) * t769 + t1092;
t659 = t706 * t788;
t1256 = (t618 + t647 - t615 + t659) * t785;
t874 = Icges(5,2) * t768 + t1080;
t649 = t874 * t788;
t1094 = Icges(4,1) * t768;
t884 = -t762 - t1094;
t663 = t884 * t788;
t1255 = (t616 - t649 - t613 + t663) * t785;
t646 = Icges(5,3) * t1055 + t739;
t658 = -Icges(4,2) * t1055 - t744;
t1254 = t619 - t646 + t614 + t658;
t648 = t874 * t785;
t662 = t884 * t785;
t1253 = t617 + t648 + t612 - t662;
t731 = Icges(3,1) * t784 + t776;
t1250 = -t215 / 0.2e1 - t213 / 0.2e1;
t1249 = t214 / 0.2e1 + t212 / 0.2e1;
t1248 = t1326 * t788;
t1247 = t1326 * t785;
t858 = t536 * t785 - t538 * t788;
t1026 = (m(6) * t858 - t1275) * t1199;
t673 = (-rSges(6,1) * t786 + rSges(6,2) * t783) * t769;
t672 = (-rSges(7,1) * t786 + rSges(7,2) * t783) * t769;
t896 = t1115 * t769 - t672;
t542 = t896 * t785;
t543 = t896 * t788;
t857 = t542 * t785 + t543 * t788;
t1031 = (-t769 * t340 + t768 * t857) * t1223 + (-t1265 * t673 - t418 * t769) * t1225;
t1245 = (t1003 * t411 + t1292 * t543 + t363 * t542 - t409 * t462) * t1318 + (t1263 * t673 - t454 * t536 + t456 * t538) * t1319;
t1244 = (t1292 * t198 + t1359 * t422 + t199 * t363 - t301 * t423) * t1318 + (t1291 * t309 + t1339 * t433 + t310 * t392 - t416 * t434) * t1319;
t927 = -t811 / 0.2e1 - t812 / 0.2e1;
t929 = -t810 / 0.2e1 - t809 / 0.2e1;
t1239 = -t783 * (t661 / 0.2e1 + t660 / 0.2e1 - t929) - t786 * (t655 / 0.2e1 + t654 / 0.2e1 + t927);
t695 = t729 * t788;
t697 = t731 * t788;
t1238 = (-t645 + t695) * t1052 + (-t643 - t697) * t1048;
t694 = -Icges(3,2) * t1048 - t759;
t696 = t731 * t785;
t1237 = (t644 + t694) * t784 + (t642 + t696) * t787;
t1235 = t783 * (t590 / 0.2e1 + t592 / 0.2e1) + t786 * (t588 / 0.2e1 + t586 / 0.2e1) + t807 / 0.2e1 + t808 / 0.2e1 + t701 / 0.2e1 - t874 / 0.2e1 - t762 + t1082 / 0.2e1 - t1094 / 0.2e1;
t1234 = t783 * t927 + t786 * t929 + t582 / 0.2e1 + t584 / 0.2e1 - Icges(5,3) * t1197 - Icges(5,6) * t768 - Icges(5,2) * t1280 - t706 / 0.2e1 + t709 / 0.2e1;
t1232 = 0.4e1 * qJD(1);
t1231 = 2 * qJD(2);
t1230 = 4 * qJD(2);
t1229 = 2 * qJD(5);
t1228 = 4 * qJD(5);
t1218 = m(6) * (t1339 * t309 - t260 * t365 - t310 * t416);
t1214 = m(7) * (t834 * t769 + (t155 + t869) * t768);
t1213 = m(7) * (t1359 * t198 + t155 * t251 - t199 * t301);
t291 = t301 * t1054;
t1209 = m(7) * (t769 * t823 - t291);
t713 = t1075 + t1117;
t779 = t788 * pkin(7);
t995 = -t785 * (pkin(1) * t785 - t779 + t973) + t788 * (-t785 * pkin(7) - t763 + (-pkin(1) + t765) * t788);
t898 = t713 * t970 + t995;
t840 = -t785 * t716 + t788 * (t959 + t1113) + t898;
t183 = t1257 * t788 - t1343 * t785 + t840;
t1207 = m(7) * (t183 * t340 - t409 * t542 - t411 * t543);
t1192 = -t788 / 0.4e1;
t1106 = rSges(3,1) * t787;
t925 = pkin(1) + t1106;
t971 = rSges(3,2) * t1052 + t788 * rSges(3,3);
t576 = -t785 * t925 + t779 + t971;
t761 = rSges(3,2) * t1051;
t577 = -t761 + t925 * t788 + (rSges(3,3) + pkin(7)) * t785;
t733 = rSges(3,1) * t784 + rSges(3,2) * t787;
t698 = t733 * t785;
t699 = t733 * t788;
t1187 = m(3) * (t576 * t698 - t577 * t699);
t851 = rSges(4,1) * t1055 - rSges(4,2) * t1057 - t788 * rSges(4,3);
t539 = -t851 + t973;
t1105 = rSges(4,1) * t769;
t917 = -rSges(4,2) * t1056 + t785 * rSges(4,3);
t540 = -t763 + (t765 + t1105) * t788 + t917;
t1185 = m(4) * (-t1267 * t540 + t1268 * t539);
t1184 = m(4) * (t539 * t788 + t540 * t785);
t1015 = -t518 * t1054 - t769 * t507;
t315 = -t785 * (rSges(5,2) * t1055 + t919) + t788 * (rSges(5,3) * t1056 + t918) + t898;
t1179 = m(5) * (t1265 * t315 + t1015);
t1175 = m(5) * (t464 * t512 + t465 * t513);
t1174 = m(5) * (t465 * t1056 - t1057 * t464);
t1173 = m(5) * (t464 * t788 + t465 * t785);
t235 = -t501 * t785 + t788 * t817 + t840;
t1167 = m(6) * (t1252 * t673 + t235 * t418);
t1017 = -t456 * t1054 - t769 * t447;
t1165 = m(6) * (t1265 * t235 + t1017);
t1160 = m(6) * (t1291 * t433 + t392 * t434);
t1158 = m(6) * (-t1291 * t538 + t392 * t536);
t1157 = m(6) * (t392 * t1056 - t1057 * t1291);
t1156 = m(6) * t1263;
t1146 = m(7) * (t251 * t691 + t768 * t869);
t1021 = -t411 * t1054 - t769 * t396;
t1144 = m(7) * (t1265 * t183 + t1021);
t824 = t1264 * t768;
t1141 = m(7) * (-t409 * t1054 + t1055 * t411 + t824);
t1139 = m(7) * (t769 * t864 + t824);
t1137 = m(7) * (t1292 * t422 + t363 * t423);
t1135 = m(7) * (-t1055 * t1359 - t291);
t1133 = m(7) * (-t1003 * t1292 + t363 * t462);
t356 = t363 * t1056;
t1132 = m(7) * (-t1057 * t1292 + t356);
t1131 = m(7) * t1264;
t1129 = m(7) * (t768 * t340 + t769 * t857);
t1124 = t769 * t1275;
t766 = t768 ^ 2;
t767 = t769 ^ 2;
t972 = t970 * t767;
t1122 = m(7) * (-t767 + (0.1e1 - t970) * t766 + t972);
t1121 = m(7) * (-t691 * t769 - t766 * t970);
t1120 = m(7) * (t1265 * t768 + t972);
t1109 = m(7) * qJD(1);
t1108 = m(7) * qJD(2);
t794 = (t309 * t785 - t310 * t788) * t1319 + (t198 * t785 - t199 * t788) * t1318;
t819 = (-t542 * t788 + t785 * t543) * t1223;
t69 = t819 + t794;
t1058 = t69 * qJD(3);
t969 = qJD(1) * t768;
t968 = qJD(1) * t769;
t967 = qJD(2) * t768;
t966 = qJD(2) * t769;
t965 = qJD(5) * t768;
t964 = qJD(5) * t769;
t344 = (t1223 + t1225 + t1226) * t1282;
t963 = t344 * qJD(1);
t514 = m(7) * t691;
t962 = t514 * qJD(1);
t805 = t769 * t833 - t1070;
t189 = t549 * t684 - t553 * t685 + t785 * t805;
t804 = t769 * t832 - t1069;
t190 = t548 * t684 - t552 * t685 + t785 * t804;
t803 = t769 * t831 - t1068;
t191 = t551 * t684 - t555 * t685 + t785 * t803;
t802 = t769 * t830 - t1067;
t192 = t550 * t684 - t554 * t685 + t785 * t802;
t956 = (t1307 - t1364) * t1199 + ((t189 + t191) * t788 + (t190 + t192) * t785 + t1303) * t1197;
t193 = t682 * t549 + t683 * t553 + t788 * t805;
t194 = t682 * t548 + t683 * t552 + t788 * t804;
t195 = t682 * t551 + t683 * t555 + t788 * t803;
t196 = t682 * t550 + t683 * t554 + t788 * t802;
t955 = (-t873 - t872 + t1306) * t1199 + ((t193 + t195) * t788 + (t194 + t196) * t785 + t1304) * t1197;
t954 = (-t1301 + t1305) * t1281 + (t1310 * t785 + t1311 * t788 + t1302) * t1280;
t953 = t1193 * t1333 + t1196 * t1334;
t941 = t1055 / 0.4e1;
t939 = t1278 * t1331 + t1279 * t1332;
t936 = t1280 * t1301 + t1281 * t1302;
t933 = t224 / 0.2e1 + t226 / 0.2e1;
t932 = -t225 / 0.2e1 - t227 / 0.2e1;
t926 = -t651 / 0.2e1 - t650 / 0.2e1;
t922 = -t713 - t1119;
t921 = rSges(4,2) * t768 - t1105 - t1119;
t897 = t970 * t1114;
t894 = rSges(5,2) * t769 - rSges(5,3) * t768 + t922;
t893 = t728 / 0.2e1 + t705 / 0.2e1 + t704 / 0.2e1;
t846 = t1115 * t767 - t672 * t769;
t367 = -t1003 * t768 - t785 * t846;
t368 = t462 * t768 + t788 * t846;
t866 = t367 * t788 + t368 * t785;
t849 = -t769 * pkin(8) + t922;
t848 = -t953 + t955;
t844 = -t939 - t956;
t410 = -t752 + (t922 - t990) * t785;
t412 = (t849 - t990) * t788;
t835 = t410 * t785 + t412 * t788 + t183;
t822 = t363 * t788;
t821 = t785 * (qJ(4) * t1055 - t751) + t788 * (-pkin(3) * t1056 + t737) - t897;
t806 = -t1245 + (t1309 + t1330) * t785 / 0.4e1 + (t1308 + t1329) * t1192;
t799 = -pkin(8) * t1265 + t821;
t792 = (-t1285 + t1287) * t1278 + (t1283 + t1325) * t1196 + (t1193 * t1327 + t1286 * t1278) * t788 + ((t614 * t769 - t617 * t768 + t644 * t787 - t1323) * t1193 + t1286 * t1196) * t785;
t791 = t1285 * t1279 + t1287 * t1196 + ((t1323 + t1347) * t788 + t1283 + t1324 + t1342) * t1193;
t790 = (-t1055 / 0.4e1 + t941) * (t1337 + t1338) + (t1192 + t788 / 0.4e1) * (t1341 + t1340);
t789 = -t1244 + t1305 * t1199 + t1302 * t1197 - (t1303 + t1344) * t1057 / 0.4e1 - (t1304 + t1345) * t1056 / 0.4e1 + (t1307 + t1310) * t941 + (t1306 + t1311) * t1054 / 0.4e1;
t735 = -rSges(3,2) * t784 + t1106;
t606 = t921 * t788;
t604 = t921 * t785;
t519 = t894 * t788;
t517 = t894 * t785;
t511 = t1120 / 0.2e1;
t510 = t1121 / 0.2e1;
t481 = t1122 / 0.2e1;
t457 = (-t597 + t849) * t788;
t455 = -t752 + (-t597 + t922) * t785;
t429 = -t1054 * t673 + t768 * t536;
t428 = t1055 * t673 - t538 * t768;
t395 = t858 * t769;
t371 = t781 * t886 + t788 * t974 + t821;
t357 = t363 * t1054;
t346 = 0.4e1 * t1266;
t343 = t903 * t1282 + (m(5) + m(6) + m(7)) * t1265 / 0.2e1;
t329 = -t1124 / 0.2e1;
t327 = t511 + t481 - t1121 / 0.2e1;
t326 = t510 + t511 - t1122 / 0.2e1;
t325 = t510 + t481 - t1120 / 0.2e1;
t319 = t1322 * t769;
t317 = t785 * t559 + t561 * t788 + t799;
t314 = (t768 * t651 + (t783 * t993 - t786 * t991) * t769) * t768;
t313 = (t768 * t650 + (t783 * t994 - t786 * t992) * t769) * t768;
t247 = t1001 * t785 + t1261 * t788 + t799;
t246 = t1129 / 0.2e1;
t229 = -t1055 * t1292 + t357;
t167 = t1135 / 0.2e1;
t139 = t1139 / 0.2e1;
t133 = t1141 / 0.2e1;
t111 = t183 * t691 + t768 * t865;
t105 = t1132 + t1157 + t1174;
t96 = t1146 / 0.2e1;
t93 = t195 * t785 - t196 * t788;
t92 = t193 * t785 - t194 * t788;
t91 = t191 * t785 - t192 * t788;
t90 = t189 * t785 - t190 * t788;
t89 = -t1131 - t1156 + t1173 + t1184;
t75 = t1239 * t769 - t926 * t768 + t1133 + t1158;
t71 = t1209 / 0.2e1;
t68 = t819 - t794;
t58 = t1144 + t1165 + t1179;
t43 = t133 - t1139 / 0.2e1;
t42 = t133 + t139;
t41 = t139 - t1141 / 0.2e1;
t36 = t167 + t329 - t1209 / 0.2e1;
t35 = t167 + t71 + t1124 / 0.2e1;
t34 = t329 + t71 - t1135 / 0.2e1;
t32 = t899 + t900;
t27 = t1214 / 0.2e1;
t25 = (t730 / 0.2e1 + t731 / 0.2e1) * t787 - t1296 + t1187 + t1185 + t1175 + t1160 + t1137 - t1235 * t769 + t1234 * t768;
t22 = t1025 + t1032;
t19 = t1033 + t1098 - t1026;
t18 = t1026 - t1289;
t17 = t1026 + t1289;
t16 = t96 + t246 - t1214 / 0.2e1;
t15 = t96 + t27 - t1129 / 0.2e1;
t14 = t246 + t27 - t1146 / 0.2e1;
t12 = t946 + t947;
t10 = t785 * t953 + t788 * t939 + t1167 + t1207;
t9 = t1096 + t1110 - t1031;
t8 = t1031 - t1288;
t7 = t1031 + t1288;
t5 = t785 * t792 + t788 * t791;
t4 = t1218 + t1213 + (t785 * t956 + t788 * t955 - t936) * t769 + (-t785 * t937 - t788 * t938 - t954) * t768;
t3 = t806 + (-t376 / 0.2e1 - t375 / 0.2e1 + (-t243 / 0.4e1 - t242 / 0.4e1 - t214 / 0.4e1 - t212 / 0.4e1) * t788 + (-t241 / 0.4e1 - t240 / 0.4e1 - t215 / 0.4e1 - t213 / 0.4e1) * t785) * t769 + (-t257 / 0.2e1 - t256 / 0.2e1 + (t348 / 0.4e1 + t347 / 0.4e1 + t306 / 0.4e1 + t303 / 0.4e1) * t788 + (t352 / 0.4e1 + t351 / 0.4e1 - t307 / 0.4e1 - t304 / 0.4e1) * t785) * t768 + t790 + t1244;
t2 = (t271 / 0.4e1 + t270 / 0.4e1 + t227 / 0.4e1 + t225 / 0.4e1) * t788 + (-t269 / 0.4e1 - t268 / 0.4e1 - t226 / 0.4e1 - t224 / 0.4e1) * t785 + t789 + t790 + t1245;
t1 = t806 + t789;
t6 = [t25 * qJD(2) + t89 * qJD(3) + t105 * qJD(4) + t75 * qJD(5) + t229 * t1107, t25 * qJD(1) + t32 * qJD(3) + t12 * qJD(4) + t1 * qJD(5) + t42 * qJD(6) + ((t539 * t606 + t540 * t604) * t1227 + (t464 * t519 + t465 * t517 - t512 * t518 - t513 * t516) * t1226 + (t1292 * t412 + t363 * t410 - t409 * t423 - t411 * t422) * t1223 + (t1291 * t457 + t392 * t455 - t433 * t456 - t434 * t454) * t1225) * t1231 + ((-t241 / 0.2e1 + m(3) * (-t576 * t735 - t698 * t733) - t791 + t893 * t788 + (-t644 / 0.2e1 - t694 / 0.2e1) * t787 + (t642 / 0.2e1 + t696 / 0.2e1) * t784 - t240 / 0.2e1 + t1250) * qJD(2) + (-t614 / 0.2e1 - t658 / 0.2e1 - t619 / 0.2e1 + t646 / 0.2e1) * t966 + (t612 / 0.2e1 - t662 / 0.2e1 + t617 / 0.2e1 + t648 / 0.2e1) * t967) * t788 + ((t242 / 0.2e1 + m(3) * (-t577 * t735 + t699 * t733) - t792 + t243 / 0.2e1 + (-t643 / 0.2e1 - t697 / 0.2e1) * t784 + (t645 / 0.2e1 - t695 / 0.2e1) * t787 + t893 * t785 + t1249) * qJD(2) + (t615 / 0.2e1 - t659 / 0.2e1 - t618 / 0.2e1 - t647 / 0.2e1) * t966 + (-t613 / 0.2e1 + t663 / 0.2e1 + t616 / 0.2e1 - t649 / 0.2e1) * t967) * t785, qJD(1) * t89 + qJD(2) * t32 + qJD(4) * t343 + qJD(5) * t22, qJD(1) * t105 + qJD(2) * t12 + qJD(3) * t343 + qJD(5) * t17, t75 * qJD(1) + t1 * qJD(2) + t22 * qJD(3) + t17 * qJD(4) + (t313 + t314) * qJD(5) + t36 * qJD(6) + ((-t1003 * t1359 + t1292 * t367 - t301 * t462 + t363 * t368) * t1223 + (t1291 * t428 - t1339 * t538 + t392 * t429 - t416 * t536) * t1225) * t1229 + ((t268 / 0.2e1 + t269 / 0.2e1 + t933) * t788 + (t271 / 0.2e1 + t270 / 0.2e1 - t932) * t785) * t964, t42 * qJD(2) + t36 * qJD(5) + t1109 * t229; t5 * qJD(2) + t33 * qJD(3) + t13 * qJD(4) + t3 * qJD(5) + t43 * qJD(6) + (-t1137 / 0.4e1 - t1160 / 0.4e1 - t1175 / 0.4e1 - t1185 / 0.4e1 - t1187 / 0.4e1) * t1232 + t1235 * t968 - t1234 * t969 + (t1296 - (t730 + t731) * t787 / 0.2e1) * qJD(1), t5 * qJD(1) + t58 * qJD(4) + t10 * qJD(5) + t111 * t1107 + (m(3) * ((t785 * (rSges(3,1) * t1048 - t971) + t788 * (rSges(3,1) * t1046 + t785 * rSges(3,3) - t761)) * (-t785 * t698 - t699 * t788) + t970 * t735 * t733) + m(7) * (t183 * t247 - t409 * t410 - t411 * t412) + m(6) * (t235 * t317 - t454 * t455 - t456 * t457) + m(5) * (t315 * t371 - t516 * t517 - t518 * t519) + m(4) * (-t1267 * t606 - t1268 * t604 + (t785 * t851 + t788 * (rSges(4,1) * t1054 + t917) + t995) * (-t712 * t970 - t897)) + (t92 + t93 + (t1237 * t788 - t1247 * t785 + (t1253 * t788 + t1255) * t769 + (t1254 * t788 + t1256) * t768 + t1238) * t788 + t1248 * t781) * t1196 + (t91 + t90 + (t1255 * t769 + t1256 * t768 + (t1253 * t769 + t1254 * t768 + t1237 - t1248) * t788 + t1238) * t785 + t1247 * t782) * t1193) * qJD(2), qJD(5) * t69 + t1346, t58 * qJD(2) - 0.4e1 * qJD(4) * t1266 + t7 * qJD(5) + t326 * qJD(6) + t1363, t3 * qJD(1) + t10 * qJD(2) + t1058 + t7 * qJD(4) + (t1193 * t1273 + t1196 * t1274) * qJD(5) + t16 * qJD(6) + (-t1213 / 0.4e1 - t1218 / 0.4e1) * t1228 + ((t1262 * t673 - t395 * t235 - t365 * t418 - t428 * t456 - t429 * t454) * t1225 + (t1359 * t543 + t183 * t319 + t251 * t340 - t301 * t542 - t367 * t411 - t368 * t409) * t1223) * t1229 + ((t932 + t938) * t788 + (t933 + t937) * t785 + t954) * t965 + (t785 * t844 - t788 * t848 + t936) * t964, t43 * qJD(1) + t326 * qJD(4) + t16 * qJD(5) + t111 * t1108 + t1293; -t33 * qJD(2) + t344 * qJD(4) - t24 * qJD(5) - t514 * qJD(6) + (t1131 / 0.4e1 + t1156 / 0.4e1 - t1184 / 0.4e1 - t1173 / 0.4e1) * t1232, -t1346 + t68 * qJD(5) + ((-t410 * t788 + t785 * t412) * t1223 + (-t455 * t788 + t785 * t457) * t1225 + (-t517 * t788 + t785 * t519) * t1226 + (-t604 * t788 + t785 * t606) * t1227) * t1231, 0, t963, -t1365 + t68 * qJD(2) + ((t428 * t785 - t429 * t788) * t1225 + (t367 * t785 - t368 * t788) * t1223) * t1229, -t962; -t13 * qJD(2) - t344 * qJD(3) + t18 * qJD(5) + (-t1132 / 0.4e1 - t1157 / 0.4e1 - t1174 / 0.4e1) * t1232 + 0.2e1 * (-t768 * t822 + t356) * t1223 * qJD(1), -t1363 + t346 * qJD(4) + t8 * qJD(5) + t325 * qJD(6) + (-t1144 / 0.4e1 - t1165 / 0.4e1 - t1179 / 0.4e1) * t1230 + ((-t769 * t247 + t1021) * t1223 + (-t769 * t317 + t1017) * t1225 + (-t769 * t371 + t1015) * t1226 + (t835 * t1223 + (t455 * t785 + t457 * t788 + t235) * t1225 + (t517 * t785 + t519 * t788 + t315) * t1226) * t768) * t1231, -t963, t346 * qJD(2), t18 * qJD(1) + t8 * qJD(2) + ((t395 * t769 + (t428 * t788 + t429 * t785) * t768) * t1225 + (-t319 * t769 + t768 * t866) * t1223) * t1229, t325 * qJD(2); t2 * qJD(2) + t24 * qJD(3) + t19 * qJD(4) + t35 * qJD(6) + t926 * t969 + (-t1133 / 0.4e1 - t1158 / 0.4e1) * t1232 - t1239 * t968, t2 * qJD(1) - t1058 + t9 * qJD(4) + t4 * qJD(5) + t15 * qJD(6) + (-t1207 / 0.4e1 - t1167 / 0.4e1) * t1230 + ((t1359 * t412 + t155 * t183 - t198 * t411 - t199 * t409 + t247 * t251 - t301 * t410) * t1223 + (t1339 * t457 + t235 * t260 - t309 * t456 - t310 * t454 - t317 * t365 - t416 * t455) * t1225) * t1231 + ((t93 / 0.2e1 + t92 / 0.2e1 + t304 / 0.2e1 + t307 / 0.2e1) * t788 + (t306 / 0.2e1 + t303 / 0.2e1 + t91 / 0.2e1 + t90 / 0.2e1) * t785) * t966 + ((-t1337 / 0.2e1 - t1338 / 0.2e1 + t1250) * t788 + (t1366 * t1278 + t1279 * t1367 + t1249) * t785) * t967 + (t785 * t848 + t788 * t844) * qJD(2), -qJD(2) * t69 + t1365, qJD(1) * t19 + qJD(2) * t9, t4 * qJD(2) + (t314 / 0.2e1 + t313 / 0.2e1) * t965 + ((t1359 * t367 + t251 * t319 - t301 * t368) * t1222 + (t1339 * t428 + t365 * t395 - t416 * t429) * t1224) * t1228 + ((t1308 * t785 + t1309 * t788) * t1199 + t1273 * t1196 + t1274 * t1278) * t964, qJD(1) * t35 + qJD(2) * t15; (-t769 * t822 - t229 + t357) * t1109 + t41 * qJD(2) + t514 * qJD(3) + t34 * qJD(5), t41 * qJD(1) + (t835 * t769 + (t247 + t865) * t768 - t111) * t1108 + t327 * qJD(4) + t14 * qJD(5) - t1293, t962, t327 * qJD(2), t34 * qJD(1) + t14 * qJD(2) + m(7) * (t768 * t319 + t769 * t866) * qJD(5), -t541 * t1108;];
Cq  = t6;