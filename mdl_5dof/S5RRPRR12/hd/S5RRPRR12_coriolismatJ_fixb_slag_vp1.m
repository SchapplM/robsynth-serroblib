% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR12_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR12_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:06
% EndTime: 2019-12-31 20:30:16
% DurationCPUTime: 61.61s
% Computational Cost: add. (94028->1069), mult. (245285->1472), div. (0->0), fcn. (297074->8), ass. (0->659)
t682 = sin(qJ(1));
t1142 = sin(qJ(4));
t1144 = cos(qJ(4));
t681 = sin(qJ(2));
t683 = cos(qJ(2));
t793 = -t1142 * t681 - t1144 * t683;
t1237 = t682 * t793;
t684 = cos(qJ(1));
t902 = t684 * t1142;
t904 = t684 * t1144;
t606 = -t681 * t902 - t683 * t904;
t573 = Icges(5,4) * t606;
t607 = -t681 * t904 + t683 * t902;
t495 = -Icges(5,2) * t607 - Icges(5,6) * t682 - t573;
t1270 = Icges(5,1) * t607 + t495 - t573;
t630 = -t683 * t1142 + t681 * t1144;
t605 = t630 * t682;
t570 = Icges(5,4) * t605;
t1257 = -Icges(5,2) * t1237 - t570;
t833 = Icges(5,5) * t684 + t570;
t496 = Icges(5,1) * t1237 - t833;
t1346 = t1257 + t496;
t1265 = Icges(5,4) * t1237;
t830 = Icges(5,2) * t605 + Icges(5,6) * t684;
t494 = t1265 - t830;
t1399 = Icges(5,1) * t605 + t1265 + t494;
t503 = Icges(5,5) * t607 - Icges(5,6) * t606;
t742 = Icges(5,5) * t605 + Icges(5,6) * t1237;
t572 = Icges(5,4) * t607;
t497 = -Icges(5,1) * t606 - Icges(5,5) * t682 - t572;
t962 = Icges(5,2) * t606 + t497 - t572;
t178 = (t1399 * t1237 + t1346 * t605 - t684 * t742) * t684 + (t1237 * t1270 - t684 * t503 + t962 * t605) * t682;
t180 = (-t1346 * t607 + t1399 * t606 + t682 * t742) * t684 + (t1270 * t606 + t503 * t682 - t607 * t962) * t682;
t1141 = sin(qJ(5));
t1143 = cos(qJ(5));
t794 = t1141 * t606 - t1143 * t682;
t842 = t682 * t1141 + t1143 * t606;
t433 = -Icges(6,5) * t842 + Icges(6,6) * t794 + Icges(6,3) * t607;
t912 = Icges(6,4) * t1143;
t796 = Icges(6,2) * t1141 - t912;
t456 = -Icges(6,6) * t606 + t607 * t796;
t911 = Icges(6,4) * t1141;
t797 = -Icges(6,1) * t1143 + t911;
t458 = -Icges(6,5) * t606 + t607 * t797;
t765 = t793 * t1143;
t901 = t684 * t1141;
t735 = -t682 * t765 + t901;
t764 = t793 * t1141;
t903 = t684 * t1143;
t736 = t682 * t764 + t903;
t1239 = Icges(6,4) * t842;
t435 = Icges(6,2) * t794 + Icges(6,6) * t607 - t1239;
t535 = Icges(6,4) * t794;
t437 = -Icges(6,1) * t842 + Icges(6,5) * t607 + t535;
t1210 = t1141 * t435 - t1143 * t437;
t795 = -Icges(6,5) * t1143 + Icges(6,6) * t1141;
t783 = -Icges(6,3) * t606 + t607 * t795 + t1210;
t189 = t1237 * t433 - t456 * t736 - t458 * t735 + t783 * t605;
t1050 = Icges(6,3) * t605;
t706 = Icges(6,5) * t735 + Icges(6,6) * t736 - t1050;
t1276 = t706 * t793;
t1054 = Icges(6,6) * t605;
t707 = Icges(6,4) * t735 + Icges(6,2) * t736 - t1054;
t1060 = Icges(6,5) * t605;
t708 = Icges(6,1) * t735 + Icges(6,4) * t736 - t1060;
t1217 = t707 * t1141 - t708 * t1143;
t909 = t605 * t1141;
t910 = t605 * t1143;
t1377 = -Icges(6,5) * t910 + Icges(6,6) * t909 + Icges(6,3) * t1237 - t1217;
t719 = -Icges(6,4) * t910 + Icges(6,2) * t909 + Icges(6,6) * t1237;
t721 = -Icges(6,1) * t910 + Icges(6,4) * t909 + Icges(6,5) * t1237;
t687 = -t682 * t1276 + t1377 * t605 - t719 * t736 - t721 * t735;
t95 = t189 * t682 + t684 * t687;
t192 = -t606 * t433 + t456 * t794 - t458 * t842 + t607 * t783;
t689 = t1377 * t607 + t606 * t706 + t719 * t794 - t721 * t842;
t98 = t192 * t682 + t684 * t689;
t757 = -(-t98 / 0.2e1 - t180 / 0.2e1) * t682 - (-t95 / 0.2e1 + t178 / 0.2e1) * t684;
t1147 = t684 / 0.2e1;
t1153 = -t682 / 0.2e1;
t1200 = (t180 + t98) * t1153 + (t178 - t95) * t1147;
t306 = -t1217 * t630 - t1276;
t308 = -t1210 * t630 - t433 * t793;
t1005 = t630 * (-t306 * t684 + t308 * t682);
t696 = t607 * t706 + t707 * t794 - t708 * t842;
t1321 = t684 * t696;
t283 = t607 * t433 + t435 * t794 - t437 * t842;
t162 = t283 * t682 - t1321;
t1019 = t606 * t162;
t477 = -Icges(6,3) * t793 - t630 * t795;
t478 = -Icges(6,6) * t630 + t793 * t796;
t481 = -Icges(6,5) * t630 + t793 * t797;
t480 = -Icges(6,6) * t793 - t630 * t796;
t891 = t1141 * t480;
t483 = -Icges(6,5) * t793 - t630 * t797;
t896 = t1143 * t483;
t1211 = t891 - t896;
t475 = -Icges(6,3) * t630 + t793 * t795;
t782 = t475 - t1211;
t237 = -t1237 * t477 - t478 * t736 - t481 * t735 + t782 * t605;
t1022 = t605 * t433;
t282 = t435 * t736 + t437 * t735 - t1022;
t1021 = t605 * t477;
t333 = t480 * t736 + t483 * t735 - t1021;
t705 = t605 * t706;
t693 = t707 * t736 + t708 * t735 - t705;
t692 = t693 * t682;
t29 = t605 * t687 + t189 * t607 + t282 * t606 - t333 * t630 + (t237 + t692) * t793;
t1408 = t684 * t29;
t1409 = t605 * t95;
t907 = t630 * t1141;
t908 = t630 * t1143;
t213 = t1377 * t793 + t630 * t706 + t719 * t907 - t721 * t908;
t1383 = t213 * t684;
t212 = t783 * t793 + (t1141 * t456 - t1143 * t458 - t433) * t630;
t1400 = t212 * t682;
t1419 = t793 * (-t1400 - t1383);
t238 = t606 * t477 + t478 * t794 - t481 * t842 + t607 * t782;
t335 = t607 * t477 + t480 * t794 - t483 * t842;
t695 = t682 * t696;
t32 = t192 * t607 - t283 * t606 + t335 * t630 + t605 * t689 + (t238 - t695) * t793;
t1422 = t682 * t32;
t1423 = t607 * t98;
t1151 = t682 / 0.2e1;
t1224 = Icges(6,2) * t842 + t437 + t535;
t1225 = -Icges(6,1) * t794 - t1239 + t435;
t465 = Icges(6,5) * t794 + Icges(6,6) * t842;
t218 = t1224 * t736 - t1225 * t735 - t605 * t465;
t1215 = t708 - t796 * t684 + (Icges(6,4) * t764 + Icges(6,2) * t765) * t682;
t1216 = t707 + t797 * t684 - (Icges(6,1) * t764 + Icges(6,4) * t765) * t682;
t712 = -t795 * t684 + (Icges(6,5) * t764 + Icges(6,6) * t765) * t682;
t685 = t1215 * t736 - t1216 * t735 - t605 * t712;
t120 = t218 * t682 - t684 * t685;
t219 = t1224 * t794 + t1225 * t842 + t607 * t465;
t688 = t1215 * t794 + t1216 * t842 + t607 * t712;
t121 = t219 * t682 - t684 * t688;
t1241 = -t684 / 0.2e1;
t805 = t121 * t1151 + t120 * t1241;
t1425 = t1019 / 0.2e1 - t1005 / 0.2e1 - t805 - t1408 / 0.2e1 - t1409 / 0.2e1 + t1419 / 0.2e1 - t1422 / 0.2e1 - t1423 / 0.2e1;
t1424 = t1408 / 0.4e1 + t1409 / 0.4e1 - t1419 / 0.4e1 + t1422 / 0.4e1 + t1423 / 0.4e1;
t351 = -t1211 * t630 - t477 * t793;
t1032 = t351 * t630;
t1034 = t308 * t606;
t1154 = -t793 / 0.2e1;
t1155 = t607 / 0.2e1;
t1156 = -t605 / 0.2e1;
t1168 = m(6) / 0.4e1;
t1314 = t1143 * t481;
t1315 = t1141 * t478;
t1379 = ((t1315 - t1314 + t477) * t630 + t782 * t793) * t793;
t1384 = t213 * t605;
t1401 = t212 * t607;
t1411 = t630 / 0.2e1;
t1307 = -t335 * t793 - t605 * t696;
t128 = t283 * t607 + t1307;
t1413 = -t128 / 0.2e1;
t724 = -rSges(6,1) * t842 + rSges(6,2) * t794 + t607 * rSges(6,3);
t1275 = t724 * t793;
t587 = t605 * rSges(6,3);
t438 = rSges(6,1) * t735 + rSges(6,2) * t736 - t587;
t459 = -rSges(6,1) * t910 + rSges(6,2) * t909 + rSges(6,3) * t1237;
t799 = -t1143 * rSges(6,1) + t1141 * rSges(6,2);
t461 = t606 * rSges(6,3) - t799 * t607;
t258 = -t682 * t1275 + t606 * t438 + t607 * t459 + t605 * t461;
t484 = -t630 * rSges(6,3) + t793 * t799;
t486 = -rSges(6,3) * t793 - t630 * t799;
t287 = t1237 * t486 + t630 * t438 + t459 * t793 - t605 * t484;
t289 = t461 * t793 + t607 * t484 + t606 * t486 + t630 * t724;
t350 = t607 * t438 + t605 * t724;
t1030 = t486 * t605;
t372 = t438 * t793 - t1030;
t374 = t486 * t607 + t1275;
t1414 = -t1154 * (t1237 * t306 - t1032 + t1034 + t1379 + t1384 + t1401) + t1155 * t32 - t1156 * t29 - 0.4e1 * (t258 * t350 + t287 * t372 + t289 * t374) * t1168 + t606 * t1413 + (-t306 * t605 + t308 * t607 - t351 * t793) * t1411;
t1412 = t607 / 0.4e1;
t1410 = -t1400 / 0.4e1;
t1290 = t605 / 0.4e1;
t1404 = -t237 * t1290 - t1034 / 0.4e1 + t1032 / 0.2e1 - t606 * t335 / 0.4e1 - t1379 / 0.2e1 - t1384 / 0.4e1 - t238 * t1412 - t1401 / 0.4e1;
t723 = -t606 * pkin(4) + t607 * pkin(8) + t724;
t1259 = -pkin(4) * t1237 - t605 * pkin(8);
t973 = t1259 + t438;
t313 = -t682 * t973 - t684 * t723;
t967 = pkin(4) * t630 - pkin(8) * t793 + t486;
t424 = t967 * t682;
t427 = t967 * t684;
t969 = -pkin(4) * t793 - pkin(8) * t630 + t484;
t1351 = t969 * t684;
t1352 = t969 * t682;
t1263 = t607 * pkin(4) + t606 * pkin(8) + t461;
t971 = -t605 * pkin(4) + pkin(8) * t1237 + t459;
t1371 = -t684 * t1263 - t682 * t971;
t788 = t1351 * t372 - t1352 * t374 + t1371 * t350;
t1398 = -t258 * t313 - t287 * t427 + t289 * t424 - t788;
t628 = Icges(5,4) * t630;
t531 = Icges(5,2) * t793 + t628;
t532 = -Icges(5,1) * t793 + t628;
t1220 = t531 + t532;
t627 = Icges(5,4) * t793;
t530 = -Icges(5,2) * t630 + t627;
t534 = Icges(5,1) * t630 + t627;
t1350 = t530 + t534;
t526 = -Icges(5,5) * t793 + Icges(5,6) * t630;
t332 = t1220 * t606 - t1350 * t607 + t526 * t682;
t337 = t1270 * t630 - t793 * t962;
t1397 = t238 + t337 - t332;
t330 = t1220 * t1237 + t1350 * t605 - t684 * t526;
t1274 = t496 * t793;
t338 = -t1257 * t793 + t1399 * t630 - t1274;
t1396 = t338 + t237 + t330;
t1173 = -0.2e1 * t683;
t1395 = t1371 * t1173;
t1035 = t289 * t682;
t1036 = t287 * t684;
t1394 = t1035 - t1036;
t1393 = t237 / 0.2e1 + t330 / 0.2e1 + t338 / 0.2e1;
t1392 = t338 / 0.4e1 + t237 / 0.4e1 + t330 / 0.4e1;
t1389 = t1351 * t427 + t1352 * t424 + t1371 * t313;
t1100 = t684 * pkin(7);
t997 = t682 * t683;
t632 = pkin(3) * t997 + t1100;
t677 = t682 * pkin(7);
t1003 = t681 * qJ(3);
t647 = t683 * pkin(2) + t1003;
t679 = t682 ^ 2;
t680 = t684 ^ 2;
t945 = t679 + t680;
t949 = t945 * t647;
t995 = t683 * t684;
t850 = t682 * t632 + t684 * (pkin(3) * t995 - t677) + t949;
t280 = -t313 + t850;
t996 = t683 * qJ(3);
t643 = pkin(2) * t681 - t996;
t874 = t681 * pkin(3) + t643;
t812 = t874 + t967;
t403 = t812 * t682;
t405 = t812 * t684;
t1388 = -t405 * t1351 - t403 * t1352 + t280 * t1371;
t1387 = t337 / 0.2e1 - t332 / 0.2e1 + t212 / 0.2e1 + t238 / 0.2e1;
t1235 = t684 * (t682 * (Icges(5,5) * t1237 - Icges(5,6) * t605 - Icges(5,3) * t684) + t607 * t494 + t606 * t496);
t359 = -t607 * t495 - t606 * t497 + t682 * (Icges(5,5) * t606 + Icges(5,6) * t607 + Icges(5,3) * t682);
t106 = -t1235 + (-t1237 * t496 + t682 * t1274 + t359) * t682;
t1146 = t684 / 0.4e1;
t1148 = -t684 / 0.4e1;
t1152 = -t682 / 0.4e1;
t546 = -t1141 * t1237 - t903;
t547 = t1143 * t1237 - t901;
t439 = t547 * rSges(6,1) + t546 * rSges(6,2) + t587;
t1223 = -t1259 + t439;
t1233 = (t973 + t1223) * t684;
t1250 = 0.2e1 * m(6);
t1284 = -t1250 / 0.4e1;
t984 = -t280 + t313;
t1206 = t984 * t1284 * t1233;
t161 = t282 * t682 - t684 * t693;
t1227 = -t546 * t707 - t547 * t708 - t705;
t432 = Icges(6,5) * t547 + Icges(6,6) * t546 + t1050;
t434 = Icges(6,4) * t547 + Icges(6,2) * t546 + t1054;
t436 = Icges(6,1) * t547 + Icges(6,4) * t546 + t1060;
t165 = t607 * t432 + t434 * t794 + t546 * t435 - t436 * t842 + t547 * t437 + t1022;
t61 = t1227 * t684 + t165 * t682 + t695;
t1323 = t161 + t61;
t265 = t359 * t682 - t1235;
t166 = -t605 * t432 + t434 * t736 + t436 * t735 + t283;
t62 = t166 * t682 - t1321 + t692;
t1385 = -t1206 + 0.2e1 * (t106 + t62) * t1146 + 0.2e1 * (t162 + t265) * t1148 + 0.2e1 * t1323 * t1152 + (t212 / 0.4e1 + t238 / 0.4e1 + t337 / 0.4e1 - t332 / 0.4e1) * t682;
t1378 = t258 * t1173;
t1298 = t530 / 0.2e1 + t475 / 0.2e1 + t534 / 0.2e1 + t896 / 0.2e1 - t891 / 0.2e1;
t740 = t1314 / 0.2e1 - t1315 / 0.2e1 - t477 / 0.2e1 + t531 / 0.2e1 + t532 / 0.2e1;
t943 = qJD(1) * t630;
t944 = qJD(1) * t793;
t1370 = -t1298 * t944 + t740 * t943;
t1369 = t1298 * t793 - t630 * t740;
t1002 = t681 * t682;
t678 = t684 * pkin(6);
t731 = -t682 * pkin(1) - pkin(2) * t997 - qJ(3) * t1002 - t632 + t678;
t385 = t731 - t973;
t1167 = pkin(2) + pkin(3);
t1197 = t1167 * t683 + pkin(1) + t1003;
t752 = t682 * pkin(6) + t1197 * t684 - t677;
t386 = t752 + t723;
t1368 = -t287 * t385 + t289 * t386;
t824 = -t1351 * t385 - t1352 * t386;
t821 = -t1351 * t684 - t1352 * t682;
t1363 = -t258 * t280 + t287 * t405 - t289 * t403;
t1069 = Icges(4,1) * t683;
t670 = Icges(4,5) * t681;
t835 = t670 + t1069;
t583 = Icges(4,4) * t682 + t684 * t835;
t1064 = Icges(3,4) * t681;
t642 = Icges(3,1) * t683 - t1064;
t585 = Icges(3,5) * t682 + t642 * t684;
t1349 = t583 + t585;
t510 = t605 * rSges(5,1) + rSges(5,2) * t1237;
t511 = t607 * rSges(5,1) - t606 * rSges(5,2);
t1327 = t682 * t510 - t511 * t684;
t951 = rSges(5,1) * t1237 - t605 * rSges(5,2);
t498 = t684 * rSges(5,3) - t951;
t500 = -t606 * rSges(5,1) - t607 * rSges(5,2) - t682 * rSges(5,3);
t410 = -t682 * t498 - t500 * t684;
t369 = -t410 + t850;
t1347 = t1327 * t369;
t1345 = t1327 * t1173;
t349 = (-t438 - t439) * t605;
t1184 = -0.2e1 * t349;
t307 = -t793 * t432 + (-t1141 * t434 + t1143 * t436) * t630;
t1238 = t793 * (t306 + t307);
t127 = t282 * t607 - t333 * t793 - t605 * t693;
t1287 = t682 / 0.4e1;
t334 = t480 * t546 + t483 * t547 + t1021;
t49 = -t334 * t793 + (t165 + t696) * t607 + t1227 * t605;
t50 = (t166 + t693) * t607 + t1307;
t1229 = -t128 * t1146 - t50 * t1148 + (t1238 - t127 - t49) * t1287;
t1261 = 0.2e1 * t350 * t1233;
t373 = -t439 * t793 - t1030;
t980 = t372 - t373;
t1332 = (t280 * t1184 + 0.2e1 * t403 * t980 + t1261) * t1168 - t1229;
t1331 = t1229 + (t313 * t1184 - 0.2e1 * t424 * t980 - t1261) * t1168;
t442 = t731 - t498;
t768 = -t1197 * t682 + t678;
t972 = (-pkin(7) - rSges(5,3)) * t684 + t768 + t951 - t442;
t1325 = m(5) * t972;
t1324 = t1250 / 0.4e1;
t1001 = t681 * t684;
t658 = Icges(4,5) * t995;
t575 = Icges(4,6) * t682 + Icges(4,3) * t1001 + t658;
t635 = Icges(3,5) * t683 - Icges(3,6) * t681;
t577 = Icges(3,3) * t682 + t635 * t684;
t636 = Icges(4,4) * t683 + Icges(4,6) * t681;
t579 = Icges(4,2) * t682 + t636 * t684;
t1319 = t575 * t1001 + t1349 * t995 + (t577 + t579) * t682;
t1317 = (-Icges(3,6) + Icges(4,6)) * t683 + (-Icges(4,4) - Icges(3,5)) * t681;
t637 = Icges(3,2) * t683 + t1064;
t1051 = Icges(4,3) * t683;
t831 = t1051 - t670;
t1316 = (-t637 - t831) * t684 + t1349;
t1311 = -(t1263 * t374 + t372 * t971 + t1368) * t1284 + t1404;
t1226 = t333 + t306;
t1285 = t1237 / 0.4e1;
t1310 = t1226 * t1285;
t1286 = -t1237 / 0.4e1;
t1309 = t1226 * t1286;
t468 = -t799 * t684 + (rSges(6,1) * t764 + rSges(6,2) * t765) * t682;
t1176 = -0.2e1 * t468;
t1177 = 0.2e1 * t424;
t523 = (-rSges(6,1) * t1141 - rSges(6,2) * t1143) * t630;
t938 = 0.2e1 * t523;
t1234 = (t385 * t684 + t386 * t682) * t938;
t469 = rSges(6,1) * t794 + rSges(6,2) * t842;
t247 = -t793 * t465 + (-t1141 * t1224 - t1143 * t1225) * t630;
t1039 = t247 * t682;
t246 = -t1215 * t907 - t1216 * t908 - t712 * t793;
t1040 = t246 * t684;
t521 = (-Icges(6,2) * t1143 - t911) * t630;
t1221 = t483 + t521;
t522 = (-Icges(6,1) * t1141 - t912) * t630;
t1222 = t480 - t522;
t520 = (-Icges(6,5) * t1141 - Icges(6,6) * t1143) * t630;
t273 = t1221 * t736 - t1222 * t735 - t605 * t520;
t274 = t1221 * t794 + t1222 * t842 + t607 * t520;
t853 = t1040 / 0.4e1 - t1039 / 0.4e1 + t274 * t1152 + t273 * t1146;
t840 = (t1176 * t427 + t1177 * t469 + t1234) * t1168 + t853;
t1178 = -0.2e1 * t403;
t841 = (-t1176 * t405 + t1178 * t469 - t1234) * t1168 - t853;
t1212 = t1005 / 0.4e1 - t1019 / 0.4e1;
t395 = t682 * t468 + t469 * t684;
t976 = -t405 - t427;
t977 = -t403 - t424;
t1260 = (t984 * t395 + (t682 * t977 + t684 * t976) * t523) * t1324 - t805;
t1301 = t1285 * t161 - t1212 + t1260;
t1299 = t1286 * t161 + t1212;
t1232 = t945 * t681;
t838 = -0.2e1 * t1232;
t1172 = m(4) / 0.4e1;
t1170 = m(5) / 0.4e1;
t1291 = -t62 / 0.4e1;
t443 = t500 + t752;
t538 = -rSges(5,1) * t793 + rSges(5,2) * t630;
t1281 = t538 * (t442 * t684 + t443 * t682);
t540 = rSges(5,1) * t630 + rSges(5,2) * t793;
t847 = 0.4e1 * t945;
t1280 = t538 * t540 * t847;
t848 = 0.2e1 * t945;
t1279 = t538 * t848;
t1278 = t540 * t945;
t845 = t540 + t874;
t471 = t845 * t682;
t473 = t845 * t684;
t819 = t471 * t682 + t473 * t684;
t1273 = t819 * t538;
t580 = Icges(3,4) * t997 - Icges(3,2) * t1002 - Icges(3,6) * t684;
t1057 = Icges(3,2) * t681;
t673 = Icges(3,4) * t683;
t581 = Icges(3,6) * t682 + (t673 - t1057) * t684;
t556 = t585 * t997;
t860 = t577 * t684 - t556;
t576 = Icges(3,5) * t997 - Icges(3,6) * t1002 - Icges(3,3) * t684;
t659 = Icges(3,4) * t1002;
t584 = Icges(3,1) * t997 - Icges(3,5) * t684 - t659;
t959 = -t682 * t576 - t584 * t995;
t1272 = -t1001 * t580 - t1002 * t581 - t860 - t959;
t1271 = -t1001 * t581 + t1319;
t1028 = (-Icges(4,2) * t684 + t682 * t636) * t684;
t1269 = t1028 + t1319;
t1262 = t1173 * t1233;
t924 = t403 * t1001;
t402 = -0.2e1 * t924;
t923 = t471 * t1001;
t470 = -0.2e1 * t923;
t1091 = rSges(4,1) * t681;
t644 = -rSges(4,3) * t683 + t1091;
t948 = t643 + t644;
t550 = t948 * t682;
t922 = t550 * t1001;
t525 = -0.2e1 * t922;
t914 = (t402 + 0.2e1 * t924 + t1262) * t1168 + (t470 + 0.2e1 * t923) * t1170 + (t525 + 0.2e1 * t922) * t1172;
t663 = pkin(2) * t1002;
t844 = t682 * t996 - t663;
t787 = pkin(3) * t1002 - t844;
t398 = t787 - t971;
t655 = qJ(3) * t995;
t811 = -t1001 * t1167 + t655;
t399 = t811 + t1263;
t1071 = rSges(4,3) + qJ(3);
t544 = t663 + (-t1071 * t683 + t1091) * t682;
t1145 = rSges(4,1) + pkin(2);
t662 = rSges(4,3) * t995;
t545 = -t1001 * t1145 + t655 + t662;
t463 = t510 + t787;
t464 = t811 + t511;
t798 = 0.2e1 * t463 * t684 + 0.2e1 * t464 * t682;
t936 = 0.2e1 * t681;
t1194 = t1071 * t681 + t1145 * t683 + pkin(1);
t676 = t684 * rSges(4,2);
t501 = -t1194 * t682 + t676 + t678;
t502 = (rSges(4,2) + pkin(6)) * t682 + t1194 * t684;
t930 = 0.2e1 * t995;
t931 = 0.2e1 * t997;
t966 = t501 * t930 + t502 * t931;
t974 = t442 * t930 + t443 * t931;
t979 = t385 * t930 + t386 * t931;
t915 = ((t398 * t684 + t399 * t682) * t936 + t979) * t1168 + (t681 * t798 + t974) * t1170 + ((t544 * t684 + t545 * t682) * t936 + t966) * t1172;
t1268 = t914 - t915;
t1240 = m(6) * t936;
t982 = (-t1263 * t682 + t684 * t971) * t1240 / 0.4e1 + (-t510 * t684 - t511 * t682) * t936 * t1170;
t925 = t424 * t1001;
t417 = 0.2e1 * t925;
t986 = (t417 - 0.2e1 * t925 - t1262) * t1168;
t1267 = t982 - t986;
t1175 = 0.2e1 * t1232;
t960 = 0.2e1 * t683 * t1278;
t975 = t424 * t931 + t427 * t930;
t987 = (t1175 * t313 + t975) * t1168 + (t1175 * t410 + t960) * t1170;
t988 = (t1395 + (t313 - t821) * t936 + t975) * t1168 + (t1345 + (0.2e1 * t410 + t1279) * t681 + t960) * t1170;
t1266 = t987 - t988;
t1264 = -(t372 * t398 - t374 * t399 - t1368) * t1284 - t1404;
t1249 = -0.4e1 * t313;
t1244 = t161 / 0.2e1;
t760 = t1237 / 0.2e1;
t1063 = Icges(4,5) * t683;
t634 = Icges(4,3) * t681 + t1063;
t574 = -Icges(4,6) * t684 + t634 * t682;
t582 = -Icges(4,4) * t684 + t682 * t835;
t1236 = (t574 * t681 + t582 * t683) * t682;
t934 = 0.4e1 * t683;
t440 = (t1172 + t1170 + t1168) * (-0.1e1 + t945) * t681 * t934;
t1219 = t1317 * t682;
t1218 = t1317 * t684;
t1213 = t1316 * t682;
t985 = (-t821 * t936 + t1395) * t1168 + (-t538 * t838 + t1345) * t1170;
t1207 = (t1263 * t403 - t405 * t971 - t824) * t1284 - (t471 * t511 + t473 * t510 + t1281) * m(5) / 0.2e1;
t1198 = (t372 * t684 - t374 * t682) * t938 - 0.2e1 * t350 * t395;
t1196 = (t398 * t427 + t399 * t424 - t824) * t1284 - m(5) * (t540 * t798 + 0.2e1 * t1281) / 0.4e1 + t1383 / 0.4e1;
t552 = t948 * t684;
t932 = -0.2e1 * t1002;
t913 = (-t405 * t932 + t402 + t979) * t1168 + (-t473 * t932 + t470 + t974) * t1170 + (-t552 * t932 + t525 + t966) * t1172;
t1195 = -(t162 / 0.4e1 + t1291) * t605 - (t161 / 0.4e1 + t61 / 0.4e1) * t607;
t1070 = Icges(3,1) * t681;
t1092 = rSges(3,1) * t683;
t876 = pkin(1) + t1092;
t946 = rSges(3,2) * t1002 + t684 * rSges(3,3);
t548 = -t682 * t876 + t678 + t946;
t661 = rSges(3,2) * t1001;
t549 = -t661 + t876 * t684 + (rSges(3,3) + pkin(6)) * t682;
t645 = rSges(3,1) * t681 + rSges(3,2) * t683;
t624 = t645 * t682;
t626 = t645 * t684;
t639 = Icges(4,1) * t681 - t1063;
t1193 = -t681 * (t642 / 0.2e1 - t637 / 0.2e1 + t670 + t1069 / 0.2e1 - t1051 / 0.2e1) - t683 * (t673 + t1070 / 0.2e1 - t1057 / 0.2e1 + t639 / 0.2e1 - t634 / 0.2e1) - m(3) * (t548 * t624 - t549 * t626) - m(4) * (t501 * t544 + t502 * t545);
t1192 = t162 * t1290 + t605 * t1291 + t1323 * t1412;
t616 = t639 * t682;
t617 = -Icges(4,1) * t1001 + t658;
t836 = -t673 - t1070;
t618 = t836 * t682;
t619 = t836 * t684;
t1191 = ((t580 - t618 - t574 + t616) * t684 + (t575 + t617 - t581 + t619) * t682) * t683;
t81 = t218 * t607 - t273 * t793 - t605 * t685;
t82 = t219 * t607 - t274 * t793 - t605 * t688;
t1190 = (t1039 - t1040) * t1154 + t121 * t1155 + t120 * t1156 + t81 * t1241 + t82 * t1151;
t1189 = (t1363 + t788) * t1324 + t1424;
t1188 = -t1398 * t1284 - t1424;
t1185 = 0.4e1 * t280;
t1181 = -0.2e1 * t374;
t379 = t468 * t607 + t469 * t605;
t1180 = 0.2e1 * t379;
t396 = t468 * t793 - t523 * t605;
t1179 = 0.2e1 * t396;
t1174 = 0.4e1 * t1232;
t1171 = m(5) / 0.2e1;
t1169 = m(6) / 0.2e1;
t981 = t997 * t1181 + t372 * t930;
t1161 = m(6) * (t1378 + (t350 - t1394) * t936 + t981);
t1157 = t127 / 0.2e1;
t935 = 0.4e1 * t681;
t1137 = m(4) * (-t501 * t682 + t502 * t684) * t935;
t1125 = m(5) * t1327;
t276 = t1001 * t1181 + t372 * t932;
t1122 = m(6) * (-t349 * t1173 + (t373 * t682 + t374 * t684) * t936 + t276);
t822 = t403 * t682 + t405 * t684;
t937 = 0.4e1 * t523;
t1120 = m(6) * (t395 * t1185 + t822 * t937);
t1117 = m(6) * (t395 * t1249 + (t424 * t682 + t427 * t684) * t937);
t1115 = m(6) * (t1175 * t350 + t981);
t1109 = m(6) * t276;
t1031 = t395 * t683;
t1107 = m(6) * (t523 * t838 - 0.2e1 * t1031);
t1106 = m(6) * (t427 * t932 + t417);
t1105 = (-t468 * t684 + t469 * t682) * t1240;
t1099 = m(4) * qJD(2);
t1098 = m(5) * qJD(1);
t1097 = m(5) * qJD(2);
t1096 = m(6) * qJD(1);
t1095 = m(6) * qJD(2);
t1094 = m(6) * qJD(4);
t1093 = m(6) * qJD(5);
t1082 = t607 * t1238;
t1027 = t580 * t681;
t1006 = t793 * t520;
t86 = ((t395 / 0.2e1 - t258 / 0.2e1) * t683 + (t1036 / 0.2e1 - t1035 / 0.2e1 + (t679 / 0.2e1 + t680 / 0.2e1) * t523) * t681) * m(6);
t989 = t86 * qJD(3);
t978 = t385 - t768 + t1100 - t1223;
t955 = -t831 * t682 + t582;
t953 = -Icges(3,2) * t997 + t584 - t659;
t950 = t682 * t844 + t684 * (-pkin(2) * t1001 + t655);
t648 = rSges(4,1) * t683 + rSges(4,3) * t681;
t947 = -t647 - t648;
t942 = qJD(1) * t682;
t919 = t1413 + t50 / 0.2e1;
t916 = t49 / 0.2e1 + t1157;
t899 = t1143 * t480;
t895 = t1143 * t522;
t890 = t1141 * t483;
t889 = t1141 * t521;
t877 = t636 / 0.2e1 + t635 / 0.2e1;
t873 = -t683 * pkin(3) - t647;
t866 = t955 * t684;
t864 = t953 * t684;
t859 = t581 * t681 - t576;
t846 = -t538 + t873;
t843 = -t575 * t1002 + t579 * t684 - t583 * t997;
t837 = t681 * t848;
t472 = t846 * t682;
t474 = t846 * t684;
t818 = t472 * t682 + t474 * t684;
t816 = -t550 * t682 - t552 * t684;
t813 = t873 - t969;
t808 = t106 / 0.2e1 - t265 / 0.2e1 - t162 / 0.2e1 + t62 / 0.2e1;
t807 = t61 / 0.2e1 + t1244;
t804 = t333 / 0.2e1 + t307 / 0.2e1 + t306 / 0.2e1 + t334 / 0.2e1;
t800 = -pkin(3) * t1232 + t950;
t786 = -t1188 + t1299;
t785 = t1189 + t1299;
t777 = t1192 + t1332;
t776 = -t1192 + t1331;
t767 = -t1396 * t1146 + t1397 * t1152 - t1196 + t1410;
t766 = -t1383 / 0.4e1 + t1410 - t1207 - t1397 * t1287 + t1396 * t1148;
t756 = t1264 + t1310;
t755 = t1309 + t1311;
t444 = -t1028 + t1236;
t749 = t843 * t1153 + t807 + (-t682 * (-t584 * t683 + t1027) - t576 * t684 + t444) * t1241 + (t684 * t859 - t1269 + t1271) * t1147 + (t574 * t1001 + t582 * t995 + t682 * t859 + t1272 + t843 + t860) * t1151;
t748 = -t808 + (t444 - t1236 + t1269) * t1153 + t1271 * t1151 + (-t556 + (t577 + t1027) * t684 + t959 + t1272) * t1241;
t745 = t804 + (t830 + t494) * t793 / 0.2e1 + (t833 + t496) * t1411 + (-t1237 + t760) * t534;
t730 = t682 * t807 - t684 * t808;
t729 = t890 / 0.2e1 - t895 / 0.2e1 + t899 / 0.2e1 + t889 / 0.2e1;
t14 = -t1157 * t1237 + t1414;
t15 = t127 * t760 - t1414;
t649 = -rSges(3,2) * t681 + t1092;
t553 = t947 * t684;
t551 = t947 * t682;
t452 = t684 * (-rSges(4,1) * t1001 + t662) - t644 * t679 + t950;
t430 = t682 * (t648 * t682 - t676) + (t682 * rSges(4,2) + t648 * t684) * t684 + t949;
t406 = t813 * t684;
t404 = t813 * t682;
t397 = -t469 * t793 - t523 * t607;
t391 = t1105 / 0.4e1;
t390 = -t1327 + t800;
t376 = t1174 * t430 + t816 * t934;
t360 = (-t442 * t682 + t443 * t684) * t935;
t352 = t1106 / 0.4e1;
t345 = t1107 / 0.4e1;
t320 = -0.4e1 * t442 * t510 - 0.4e1 * t443 * t511;
t311 = -t1371 + t800;
t310 = 0.4e1 * t442 * t463 + 0.4e1 * t443 * t464;
t297 = t1174 * t369 - t819 * t934;
t296 = (-t1006 + (-t890 - t889 - t899 + t895) * t630) * t793;
t285 = 0.4e1 * t1327 * t410 + t1280;
t284 = (-t385 * t682 + t386 * t684) * t935;
t272 = t1109 / 0.4e1;
t245 = -0.4e1 * t1273 + 0.4e1 * t1347;
t240 = -0.4e1 * t385 * t468 + 0.4e1 * t386 * t469;
t215 = -0.4e1 * t1263 * t386 + 0.4e1 * t385 * t971;
t195 = t1115 / 0.4e1;
t194 = t1174 * t280 - t822 * t934;
t193 = 0.4e1 * t385 * t398 + 0.4e1 * t386 * t399;
t176 = t284 * t1168 + t360 * t1170 + t1137 / 0.4e1;
t143 = t240 * t1168 - t1006 / 0.2e1 - t729 * t630;
t139 = 0.4e1 * t1389;
t136 = t1233 * t1249;
t133 = t1122 / 0.4e1;
t122 = 0.4e1 * t1388;
t114 = t1168 * t194 + t1170 * t297 + t1172 * t376;
t113 = -0.4e1 * t349 * t350 + 0.4e1 * t374 * t980;
t112 = t352 - t1267;
t111 = t352 + t1267;
t110 = -t1106 / 0.4e1 + t982 + t986;
t109 = t1233 * t1185;
t85 = (t1394 * t936 + t523 * t837 + 0.2e1 * t1031 - t1378) * t1168;
t83 = t1168 * t215 + t1170 * t320 + t1369;
t80 = t272 + t133 - t1105 / 0.4e1;
t79 = t391 + t272 - t1122 / 0.4e1;
t78 = t391 + t133 - t1109 / 0.4e1;
t68 = t1161 / 0.4e1;
t67 = t1168 * t193 + t1170 * t310 - t1193 - t1369;
t44 = t195 + t68 - t1107 / 0.4e1;
t43 = t345 + t195 - t1161 / 0.4e1;
t42 = t345 + t68 - t1115 / 0.4e1;
t37 = t987 + t988 - t985;
t36 = t985 - t1266;
t35 = t985 + t1266;
t24 = t1117 / 0.4e1 + t805;
t23 = t1120 / 0.4e1 + t805;
t22 = t913 + t1268;
t21 = t913 - t1268;
t20 = t914 + t915 - t913;
t19 = t139 * t1168 + t285 * t1170 - t757;
t18 = t1168 * t122 + t1170 * t245 + t757;
t17 = t113 * t1168 - t1082 / 0.2e1 + t916 * t607 - t919 * t605;
t16 = t1168 * t136 + t730;
t13 = t109 * t1168 + t682 * t749 + t684 * t748;
t12 = t755 + t840 - t1195 - t1331;
t11 = t776 + t840 + t1310 - t1311;
t10 = t755 + t776 - t840;
t9 = t756 + t841 + t1195 - t1332;
t8 = -t1264 + t777 + t841 + t1309;
t7 = t756 + t777 - t841;
t6 = t766 + t767 + t730 + t1206;
t5 = t1392 * t684 + t1196 + t1385 + t766;
t4 = (t213 / 0.4e1 + t1392) * t684 + t767 + t1207 + t1385;
t3 = t785 + t786 - t1260;
t2 = -t1189 + t1301 + t786;
t1 = t1188 + t1301 + t785;
t25 = [(-m(6) * t386 * t978 + t443 * t1325) * qJD(1) + t67 * qJD(2) + t176 * qJD(3) + t83 * qJD(4) + t143 * qJD(5), t67 * qJD(1) + t21 * qJD(3) + t6 * qJD(4) + t9 * qJD(5) + (t385 * t406 + t386 * t404 - t398 * t405 - t399 * t403 - t109 / 0.4e1) * t1095 + (t442 * t474 + t443 * t472 - t463 * t473 - t464 * t471) * t1097 + (t501 * t553 + t502 * t551 - t544 * t552 - t545 * t550) * t1099 + (((t580 / 0.2e1 - t618 / 0.2e1 - t574 / 0.2e1 + t616 / 0.2e1) * t684 + (t575 / 0.2e1 + t617 / 0.2e1 - t581 / 0.2e1 + t619 / 0.2e1) * t682) * t681 + t1383 / 0.2e1 + (t877 * t684 + t1393 - t748) * t684 + (t877 * t682 + t1387 - t749) * t682 + (-t864 / 0.2e1 - t866 / 0.2e1 + t1316 * t1151) * t683 + ((-t548 * t684 - t549 * t682) * t649 + (-t624 * t684 + t626 * t682) * t645) * m(3)) * qJD(2), qJD(1) * t176 + qJD(2) * t21 + qJD(4) * t111 + qJD(5) * t79, t83 * qJD(1) + t6 * qJD(2) + t111 * qJD(3) + t12 * qJD(5) + (-t136 / 0.4e1 + t971 * t427 - t1263 * t424 + t824) * t1094 + (((-t442 * t538 - t510 * t540) * m(5) + t213 / 0.2e1 + t808 + t1393) * t684 + ((-t443 * t538 - t511 * t540) * m(5) - t807 + t1387) * t682) * qJD(4), t143 * qJD(1) + t9 * qJD(2) + t79 * qJD(3) + t12 * qJD(4) + (-t372 * t468 - t374 * t469 + t385 * t396 + t386 * t397 - t113 / 0.4e1) * t1093 + (-t296 + t1082 / 0.2e1 + (t247 / 0.2e1 + t274 / 0.2e1 - t916) * t607 - (t273 / 0.2e1 + t246 / 0.2e1 - t919) * t605) * qJD(5); t13 * qJD(2) + t22 * qJD(3) + t5 * qJD(4) + t8 * qJD(5) + (-t193 / 0.4e1 + t978 * t403) * t1096 + (-t310 / 0.4e1 - t972 * t471) * t1098 + t745 * t942 + t1193 * qJD(1) - t1370, t13 * qJD(1) + (m(4) * (t430 * t452 - t550 * t551 - t552 * t553) + m(3) * (0.4e1 * (t682 * (rSges(3,1) * t997 - t946) + t684 * (rSges(3,1) * t995 + t682 * rSges(3,3) - t661)) * (-t682 * t624 - t626 * t684) + t649 * t645 * t847) / 0.4e1 + (t280 * t311 - t403 * t404 - t405 * t406) * m(6) + (t369 * t390 - t471 * t472 - t473 * t474) * m(5) + t1200 + ((-t1219 * t682 + (t866 + t864 - t1213) * t681 + t1191) * t684 + t1218 * t679) * t1151 + ((-t1218 * t684 + t1191 + ((t953 + t955) * t684 - t1213) * t681) * t682 + t1219 * t680) * t1241) * qJD(2) + t114 * qJD(3) + t18 * qJD(4) + t23 * qJD(5), t22 * qJD(1) + t114 * qJD(2) + t35 * qJD(4) + t43 * qJD(5) + (-t440 + (t1169 + t1171 + m(4) / 0.2e1) * (t838 + t837) * t683) * qJD(3), t5 * qJD(1) + t18 * qJD(2) + t35 * qJD(3) + t1 * qJD(5) + (-t139 / 0.4e1 - t976 * t1351 - t977 * t1352 + t984 * t1371) * t1094 + ((-t285 / 0.4e1 + (-t369 + t410) * t1327 - (-t819 - t1278) * t538) * m(5) + t1200) * qJD(4), t8 * qJD(1) + t23 * qJD(2) + t43 * qJD(3) + t1 * qJD(4) + ((t1178 * t397 - t1179 * t405 + t1180 * t280 - t1198) * t1169 - t15 + t1190) * qJD(5); -qJD(1) * t1137 / 0.4e1 + t20 * qJD(2) + t110 * qJD(4) + t78 * qJD(5) + (-t284 / 0.4e1 - t978 * t1002) * t1096 + (-t360 / 0.4e1 + t972 * t1002) * t1098, t20 * qJD(1) + t440 * qJD(3) + t36 * qJD(4) + t42 * qJD(5) + (-t194 / 0.4e1 + (-t311 - t822) * t683 + (t404 * t682 + t406 * t684 + t280) * t681) * t1095 + (-t297 / 0.4e1 + (-t390 - t819) * t683 + (t369 + t818) * t681) * t1097 + (-t376 / 0.4e1 + (-t452 + t816) * t683 + (t551 * t682 + t553 * t684 + t430) * t681) * t1099, t440 * qJD(2), t110 * qJD(1) + t36 * qJD(2) + ((m(6) * t1371 + t1125) * t683 + (m(6) * t821 - t1171 * t1279) * t681) * qJD(4) + t85 * qJD(5), t78 * qJD(1) + t42 * qJD(2) + t85 * qJD(4) + (t379 * t1173 + (t396 * t684 + t397 * t682) * t936) * t1093 / 0.2e1; -t320 * t1098 / 0.4e1 + t4 * qJD(2) + t112 * qJD(3) + t16 * qJD(4) + t11 * qJD(5) + (-t215 / 0.4e1 - t978 * t424) * t1096 + (t540 * t1325 - t745) * t942 + t1370, t4 * qJD(1) - t1200 * qJD(2) + t37 * qJD(3) + t19 * qJD(4) + t2 * qJD(5) + (t311 * t313 + t404 * t424 + t406 * t427 - t122 / 0.4e1 + t1388) * t1095 + (t1347 + t410 * t390 - t245 / 0.4e1 + t818 * t540 - t1273) * t1097, qJD(1) * t112 + qJD(2) * t37 + qJD(5) * t86, t16 * qJD(1) + t19 * qJD(2) + (-t1389 * m(6) - t410 * t1125 - t1170 * t1280 + t757) * qJD(4) + t24 * qJD(5), t11 * qJD(1) + t2 * qJD(2) + t989 + t24 * qJD(4) + ((t1177 * t397 + t1179 * t427 + t1180 * t313 + t1198) * t1169 - t14 - t1190) * qJD(5); t520 * t944 / 0.2e1 + t7 * qJD(2) + t80 * qJD(3) + t10 * qJD(4) + t17 * qJD(5) + t729 * t943 + t804 * qJD(1) * t607 + (-t240 / 0.4e1 - t980 * t386 + t978 * t374) * t1096, t7 * qJD(1) + (t161 * t760 + (t311 * t350 + t372 * t406 - t374 * t404 - t1363) * m(6) - t1120 / 0.4e1 + t1425) * qJD(2) + t44 * qJD(3) + t3 * qJD(4) + t15 * qJD(5), qJD(1) * t80 + qJD(2) * t44 - qJD(4) * t86, t10 * qJD(1) + t3 * qJD(2) - t989 + (t1398 * m(6) + t1237 * t1244 - t1117 / 0.4e1 + t1425) * qJD(4) + t14 * qJD(5), t17 * qJD(1) + t15 * qJD(2) + t14 * qJD(4) + ((t350 * t379 + t372 * t396 - t374 * t397) * m(6) + t82 * t1155 + t81 * t1156 + (-t246 * t605 + t247 * t607 - t296) * t1154) * qJD(5);];
Cq = t25;