% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:53
% EndTime: 2019-03-09 08:47:45
% DurationCPUTime: 97.57s
% Computational Cost: add. (48247->1372), mult. (70225->1730), div. (0->0), fcn. (70705->10), ass. (0->626)
t1254 = Icges(3,3) + Icges(4,3);
t643 = qJ(2) + pkin(10);
t618 = sin(t643);
t619 = cos(t643);
t648 = sin(qJ(2));
t651 = cos(qJ(2));
t1253 = Icges(3,5) * t651 + Icges(4,5) * t619 - Icges(3,6) * t648 - Icges(4,6) * t618;
t1026 = Icges(4,4) * t618;
t515 = Icges(4,2) * t619 + t1026;
t604 = Icges(5,5) * t618;
t1252 = Icges(5,3) * t619 + t515 - t604;
t1016 = Icges(5,5) * t619;
t517 = Icges(5,1) * t618 - t1016;
t605 = Icges(4,4) * t619;
t1251 = Icges(4,1) * t618 + t517 + t605;
t652 = cos(qJ(1));
t1250 = t1254 * t652;
t649 = sin(qJ(1));
t1235 = t1253 * t652 + t1254 * t649;
t514 = Icges(5,4) * t619 + Icges(5,6) * t618;
t412 = Icges(5,2) * t649 + t514 * t652;
t1249 = t412 + t1235;
t956 = t649 * t651;
t958 = t648 * t649;
t960 = t619 * t649;
t962 = t618 * t649;
t1218 = -Icges(3,5) * t956 - Icges(4,5) * t960 + Icges(3,6) * t958 + Icges(4,6) * t962 + t1250;
t787 = Icges(5,1) * t619 + t604;
t415 = -Icges(5,4) * t652 + t649 * t787;
t1017 = Icges(4,5) * t652;
t581 = Icges(4,4) * t962;
t417 = Icges(4,1) * t960 - t1017 - t581;
t1242 = t415 + t417;
t416 = Icges(5,4) * t649 + t652 * t787;
t520 = Icges(4,1) * t619 - t1026;
t418 = Icges(4,5) * t649 + t520 * t652;
t1241 = t416 + t418;
t510 = Icges(5,3) * t618 + t1016;
t783 = -Icges(4,2) * t618 + t605;
t1248 = t510 - t783;
t1247 = t520 + t787;
t1004 = Icges(4,6) * t652;
t413 = Icges(4,4) * t960 - Icges(4,2) * t962 - t1004;
t1005 = Icges(3,6) * t652;
t453 = Icges(3,4) * t956 - Icges(3,2) * t958 - t1005;
t1246 = t413 * t618 + t453 * t648;
t1236 = Icges(3,5) * t648 + Icges(3,6) * t651 + (Icges(4,6) - Icges(5,6)) * t619 + (Icges(5,4) + Icges(4,5)) * t618;
t1027 = Icges(3,4) * t648;
t566 = Icges(3,1) * t651 - t1027;
t456 = Icges(3,5) * t649 + t566 * t652;
t1245 = -t418 * t960 - t456 * t956;
t407 = -Icges(5,6) * t652 + t510 * t649;
t1244 = -t407 + t413;
t1003 = Icges(5,6) * t649;
t959 = t619 * t652;
t580 = Icges(5,5) * t959;
t961 = t618 * t652;
t408 = Icges(5,3) * t961 + t1003 + t580;
t414 = Icges(4,6) * t649 + t652 * t783;
t1243 = -t408 + t414;
t1240 = t1252 * qJD(2);
t1239 = t1251 * qJD(2);
t1018 = Icges(3,5) * t652;
t600 = Icges(3,4) * t958;
t455 = Icges(3,1) * t956 - t1018 - t600;
t1238 = t417 * t619 + t455 * t651 - t1246;
t1237 = t514 + t1253;
t563 = Icges(3,2) * t651 + t1027;
t633 = Icges(3,4) * t651;
t565 = Icges(3,1) * t648 + t633;
t1226 = -t1251 * t619 + t1252 * t618 + t563 * t648 - t565 * t651;
t1234 = t1248 * qJD(2);
t1233 = t1247 * qJD(2);
t1230 = -t1235 * t652 - t1245;
t784 = -Icges(3,2) * t648 + t633;
t454 = Icges(3,6) * t649 + t652 * t784;
t1229 = t414 * t618 + t454 * t648;
t1188 = -t414 * t962 - t454 * t958 + t1230;
t799 = -t408 * t962 + t412 * t652 - t416 * t960;
t1228 = -t799 + t1188;
t955 = t651 * t652;
t1227 = t1241 * t959 + t1249 * t649 + t408 * t961 + t456 * t955;
t411 = -Icges(5,2) * t652 + t514 * t649;
t390 = t649 * t411;
t1225 = t1218 * t649 - t1242 * t959 - t407 * t961 - t455 * t955 - t390;
t1135 = t1236 * t652;
t1134 = t1236 * t649;
t1059 = sin(qJ(5));
t1060 = cos(qJ(5));
t497 = -t619 * t1059 + t618 * t1060;
t661 = (-qJD(2) + qJD(5)) * t497;
t714 = -t1059 * t618 - t1060 * t619;
t892 = qJD(1) * t649;
t217 = t652 * t661 + t714 * t892;
t647 = sin(qJ(6));
t650 = cos(qJ(6));
t854 = t652 * t1059;
t855 = t652 * t1060;
t441 = -t618 * t854 - t619 * t855;
t764 = t441 * t650 + t649 * t647;
t891 = qJD(1) * t652;
t145 = qJD(6) * t764 - t217 * t647 - t650 * t891;
t364 = t441 * t647 - t649 * t650;
t146 = qJD(6) * t364 + t217 * t650 - t647 * t891;
t440 = t714 * t649;
t361 = -t440 * t650 + t647 * t652;
t362 = -t440 * t647 - t650 * t652;
t439 = t497 * t649;
t157 = Icges(7,5) * t361 - Icges(7,6) * t362 - Icges(7,3) * t439;
t1022 = Icges(7,4) * t361;
t161 = Icges(7,2) * t362 + Icges(7,6) * t439 - t1022;
t356 = Icges(7,4) * t362;
t163 = Icges(7,1) * t361 - Icges(7,5) * t439 - t356;
t1115 = t497 * qJD(1);
t695 = t714 * qJD(2);
t216 = qJD(5) * t441 - t1115 * t649 - t652 * t695;
t442 = -t618 * t855 + t619 * t854;
t219 = t649 * t661 - t714 * t891;
t147 = -qJD(6) * t361 - t219 * t647 - t650 * t892;
t148 = -qJD(6) * t362 + t219 * t650 - t647 * t892;
t837 = qJD(2) * t1059;
t838 = qJD(2) * t1060;
t218 = -qJD(5) * t440 - t1115 * t652 - t837 * t962 - t838 * t960;
t66 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t218;
t68 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t218;
t70 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t218;
t13 = -t145 * t161 + t146 * t163 - t157 * t216 + t364 * t68 + t442 * t66 - t70 * t764;
t159 = -Icges(7,5) * t764 + Icges(7,6) * t364 + Icges(7,3) * t442;
t1021 = Icges(7,4) * t764;
t162 = Icges(7,2) * t364 + Icges(7,6) * t442 - t1021;
t357 = Icges(7,4) * t364;
t165 = -Icges(7,1) * t764 + Icges(7,5) * t442 + t357;
t65 = Icges(7,5) * t146 + Icges(7,6) * t145 - Icges(7,3) * t216;
t67 = Icges(7,4) * t146 + Icges(7,2) * t145 - Icges(7,6) * t216;
t69 = Icges(7,1) * t146 + Icges(7,4) * t145 - Icges(7,5) * t216;
t14 = t145 * t162 + t146 * t165 - t159 * t216 + t364 * t67 + t442 * t65 - t69 * t764;
t880 = qJD(1) * qJD(2);
t608 = t649 * t880;
t878 = qJD(1) * qJD(5);
t530 = -t649 * t878 + t608;
t179 = qJD(6) * t218 + t530;
t609 = t652 * t880;
t531 = -t652 * t878 + t609;
t180 = -qJD(6) * t216 + t531;
t336 = qJD(5) * t497 - t618 * t838 + t619 * t837;
t337 = qJD(5) * t714 - t695;
t884 = qJD(6) * t497;
t743 = t337 * t650 - t647 * t884;
t744 = -t337 * t647 - t650 * t884;
t105 = Icges(7,5) * t743 + Icges(7,6) * t744 + Icges(7,3) * t336;
t106 = Icges(7,4) * t743 + Icges(7,2) * t744 + Icges(7,6) * t336;
t107 = Icges(7,1) * t743 + Icges(7,4) * t744 + Icges(7,5) * t336;
t780 = Icges(7,5) * t650 - Icges(7,6) * t647;
t236 = -Icges(7,3) * t714 + t497 * t780;
t1019 = Icges(7,4) * t650;
t782 = -Icges(7,2) * t647 + t1019;
t239 = -Icges(7,6) * t714 + t497 * t782;
t1020 = Icges(7,4) * t647;
t785 = Icges(7,1) * t650 - t1020;
t242 = -Icges(7,5) * t714 + t497 * t785;
t20 = t105 * t442 + t106 * t364 - t107 * t764 + t145 * t239 + t146 * t242 - t216 * t236;
t623 = qJD(2) * t649;
t559 = -qJD(5) * t649 + t623;
t353 = qJD(6) * t442 + t559;
t889 = qJD(2) * t652;
t560 = -qJD(5) * t652 + t889;
t354 = qJD(6) * t439 + t560;
t443 = -qJD(6) * t714 + qJD(1);
t53 = t442 * t157 - t161 * t364 - t163 * t764;
t54 = t442 * t159 + t364 * t162 - t165 * t764;
t79 = t236 * t442 + t239 * t364 - t242 * t764;
t885 = qJD(6) * t336;
t1 = -t13 * t354 + t14 * t353 + t179 * t53 + t180 * t54 + t20 * t443 + t79 * t885;
t1053 = -qJD(1) / 0.2e1;
t1066 = t560 / 0.2e1;
t1072 = -t443 / 0.2e1;
t1074 = t354 / 0.2e1;
t1076 = -t353 / 0.2e1;
t778 = t161 * t647 + t163 * t650;
t11 = t157 * t336 - t714 * t66 + t778 * t337 + (-t647 * t68 + t650 * t70 + (t161 * t650 - t163 * t647) * qJD(6)) * t497;
t1125 = t236 * t439 + t239 * t362 - t242 * t361;
t51 = -t157 * t439 + t161 * t362 + t163 * t361;
t52 = -t439 * t159 - t162 * t362 + t361 * t165;
t22 = -t1125 * t443 + t353 * t52 - t354 * t51;
t1123 = t22 / 0.2e1;
t1023 = Icges(6,4) * t440;
t260 = -Icges(6,2) * t439 - Icges(6,6) * t652 + t1023;
t406 = Icges(6,4) * t441;
t261 = -Icges(6,2) * t442 - Icges(6,6) * t649 - t406;
t480 = Icges(6,4) * t497;
t343 = Icges(6,2) * t714 + t480;
t1124 = qJD(1) * (-Icges(6,1) * t714 + t343 + t480) - (-Icges(6,1) * t439 - t1023 - t260) * t560 + (Icges(6,1) * t442 + t261 - t406) * t559;
t258 = -Icges(6,5) * t441 - Icges(6,6) * t442 - Icges(6,3) * t649;
t405 = Icges(6,4) * t442;
t264 = -Icges(6,1) * t441 - Icges(6,5) * t649 - t405;
t94 = t652 * t258 + t439 * t261 - t440 * t264;
t256 = -Icges(6,5) * t440 + Icges(6,6) * t439 + Icges(6,3) * t652;
t1024 = Icges(6,4) * t439;
t263 = Icges(6,1) * t440 - Icges(6,5) * t652 - t1024;
t95 = -t256 * t649 + t260 * t442 + t263 * t441;
t1146 = t95 / 0.2e1 + t94 / 0.2e1;
t340 = Icges(6,5) * t497 + Icges(6,6) * t714;
t479 = Icges(6,4) * t714;
t346 = Icges(6,1) * t497 + t479;
t115 = t340 * t652 + t343 * t439 - t346 * t440;
t777 = -t162 * t647 + t165 * t650;
t12 = t159 * t336 - t714 * t65 + t777 * t337 + (-t647 * t67 + t650 * t69 + (-t162 * t650 - t165 * t647) * qJD(6)) * t497;
t131 = -t260 * t714 - t263 * t497;
t132 = t261 * t714 + t264 * t497;
t15 = -t147 * t161 + t148 * t163 + t157 * t218 + t361 * t70 - t362 * t68 - t439 * t66;
t16 = t147 * t162 + t148 * t165 + t159 * t218 + t361 * t69 - t362 * t67 - t439 * t65;
t190 = Icges(7,6) * t440 - t439 * t782;
t192 = Icges(7,6) * t441 + t442 * t782;
t194 = Icges(7,5) * t440 - t439 * t785;
t196 = Icges(7,5) * t441 + t442 * t785;
t21 = -t105 * t439 - t106 * t362 + t107 * t361 + t147 * t239 + t148 * t242 + t218 * t236;
t2 = -t1125 * t885 - t15 * t354 + t16 * t353 + t179 * t51 + t180 * t52 + t21 * t443;
t23 = t353 * t54 - t354 * t53 + t79 * t443;
t238 = Icges(7,6) * t497 + t714 * t782;
t241 = Icges(7,5) * t497 + t714 * t785;
t55 = -t157 * t714 + t497 * t778;
t56 = -t159 * t714 + t497 * t777;
t773 = -t239 * t647 + t242 * t650;
t88 = -t236 * t714 + t497 * t773;
t26 = t353 * t56 - t354 * t55 + t443 * t88;
t118 = Icges(6,5) * t219 - Icges(6,6) * t218 - Icges(6,3) * t892;
t120 = Icges(6,4) * t219 - Icges(6,2) * t218 - Icges(6,6) * t892;
t122 = Icges(6,1) * t219 - Icges(6,4) * t218 - Icges(6,5) * t892;
t28 = -t118 * t649 - t120 * t442 - t122 * t441 - t216 * t260 - t217 * t263 - t256 * t891;
t117 = Icges(6,5) * t217 + Icges(6,6) * t216 - Icges(6,3) * t891;
t119 = Icges(6,4) * t217 + Icges(6,2) * t216 - Icges(6,6) * t891;
t121 = Icges(6,1) * t217 + Icges(6,4) * t216 - Icges(6,5) * t891;
t29 = -t117 * t649 - t119 * t442 - t121 * t441 + t216 * t261 + t217 * t264 - t258 * t891;
t30 = t118 * t652 + t120 * t439 - t122 * t440 + t218 * t260 - t219 * t263 - t256 * t892;
t31 = t117 * t652 + t119 * t439 - t121 * t440 - t218 * t261 + t219 * t264 - t258 * t892;
t116 = -t340 * t649 - t343 * t442 - t346 * t441;
t800 = t258 * t649 + t442 * t261 + t441 * t264;
t36 = t116 * qJD(1) - t559 * t800 - t560 * t95;
t181 = Icges(6,5) * t337 - Icges(6,6) * t336;
t182 = Icges(6,4) * t337 - Icges(6,2) * t336;
t183 = Icges(6,1) * t337 - Icges(6,4) * t336;
t37 = -t181 * t649 - t182 * t442 - t183 * t441 + t216 * t343 + t217 * t346 - t340 * t891;
t38 = t181 * t652 + t182 * t439 - t183 * t440 - t218 * t343 + t219 * t346 - t340 * t892;
t41 = t120 * t714 + t122 * t497 + t260 * t336 - t263 * t337;
t42 = t119 * t714 + t121 * t497 - t261 * t336 + t264 * t337;
t656 = (Icges(7,3) * t497 + t714 * t780 - t773) * t443 - (Icges(7,3) * t441 + t442 * t780 + t777) * t353 + (Icges(7,3) * t440 - t439 * t780 + t778) * t354;
t681 = qJD(1) * (Icges(6,2) * t497 - t346 - t479) + (Icges(6,2) * t440 + t1024 - t263) * t560 - (Icges(6,2) * t441 + t264 - t405) * t559;
t705 = qJD(1) * (-Icges(6,5) * t714 + Icges(6,6) * t497) - (-Icges(6,5) * t439 - Icges(6,6) * t440) * t560 + (Icges(6,5) * t442 - Icges(6,6) * t441) * t559;
t93 = t256 * t652 - t260 * t439 + t263 * t440;
t1224 = (-t41 * t652 + t42 * t649 + (t131 * t649 + t132 * t652) * qJD(1)) * t1053 - t179 * (-t51 * t652 + t52 * t649) / 0.2e1 - t180 * (-t53 * t652 + t54 * t649) / 0.2e1 - (-t55 * t652 + t56 * t649) * t885 / 0.2e1 - t26 * t884 / 0.2e1 - (qJD(1) * t115 + t559 * t94 - t560 * t93 + t22) * t892 / 0.2e1 - (t23 + t36) * t891 / 0.2e1 - (qJD(1) * t37 - t28 * t560 + t29 * t559 + t1) * t649 / 0.2e1 + (qJD(1) * t38 - t30 * t560 + t31 * t559 + t2) * t652 / 0.2e1 + (t440 * t1123 + t441 * t23 / 0.2e1) * qJD(6) + (t1146 * t652 + t649 * t800) * t531 + (-t1146 * t649 + t652 * t93) * t530 + ((-t236 * t441 + t238 * t364 - t241 * t764) * t443 - (t159 * t441 + t192 * t364 - t196 * t764) * t353 + (t157 * t441 + t190 * t364 - t194 * t764) * t354 + t656 * t442 + (-t440 * t53 - t441 * t54 + t497 * t79) * qJD(6) - t13 * t652 + t14 * t649 + (t53 * t649 + t54 * t652) * qJD(1)) * t1076 + (-(t236 * t440 + t238 * t362 - t241 * t361) * t443 + (-t159 * t440 + t192 * t362 - t196 * t361) * t353 - (-t157 * t440 + t190 * t362 - t194 * t361) * t354 - t656 * t439 + (-t1125 * t497 - t440 * t51 - t441 * t52) * qJD(6) - t15 * t652 + t16 * t649 + (t51 * t649 + t52 * t652) * qJD(1)) * t1074 - (t1124 * t441 + t681 * t442 + (-qJD(1) * t800 - t28) * t652 + (qJD(1) * t95 + t29 + t705) * t649) * t559 / 0.2e1 + (t1124 * t440 - t439 * t681 + (qJD(1) * t93 + t31) * t649 + (qJD(1) * t94 - t30 - t705) * t652) * t1066 + (-(-t88 * qJD(6) + (t238 * t647 - t241 * t650 - t236) * t443 - (t192 * t647 - t196 * t650 + t159) * t353 + (t190 * t647 - t194 * t650 + t157) * t354) * t497 - t656 * t714 + (-t440 * t55 - t441 * t56) * qJD(6) - t11 * t652 + t12 * t649 + (t55 * t649 + t56 * t652) * qJD(1)) * t1072;
t1223 = rSges(3,2) * t648;
t1222 = t1240 * t652 + (t649 * t783 - t1004 - t407) * qJD(1);
t1221 = t1240 * t649 + (t510 * t652 + t1003 - t414) * qJD(1);
t1220 = -t1239 * t652 + (-t520 * t649 + t1017 - t415) * qJD(1);
t1219 = -qJD(1) * t1241 + t1239 * t649;
t770 = t407 * t618 + t415 * t619;
t1120 = t649 * t770;
t974 = t411 * t652;
t169 = -t974 + t1120;
t1216 = t1218 * t652 + t1238 * t649 + t169;
t957 = t648 * t652;
t1185 = -t413 * t961 - t453 * t957 - t1225;
t1184 = -t414 * t961 - t454 * t957 + t1227;
t1215 = t1226 * t649 + t1135;
t1214 = -t1226 * t652 + t1134;
t1213 = t1236 * qJD(2);
t1212 = t1241 * t619 + t408 * t618 + t456 * t651 - t1229;
t1211 = -t770 - t1238;
t1138 = t1241 * t618 + t1243 * t619 + t454 * t651 + t456 * t648;
t1137 = t1242 * t618 + t1244 * t619 + t453 * t651 + t455 * t648;
t1205 = t1244 * t652 + (-Icges(5,1) * t961 + t517 * t652 - t1243 + t580) * t649;
t1204 = t1251 - t1248;
t1203 = -t1252 + t1247;
t1202 = (Icges(4,2) * t960 - t1242 + t581) * t652 + (-t515 * t652 + t1241) * t649;
t533 = t784 * qJD(2);
t534 = t566 * qJD(2);
t1201 = -t533 * t648 + t534 * t651 + t1233 * t619 + t1234 * t618 + (-t1251 * t618 - t1252 * t619 - t563 * t651 - t565 * t648) * qJD(2) + t1236 * qJD(1);
t1200 = t1249 * qJD(1);
t1199 = t1226 * qJD(1) + qJD(2) * t1237;
t166 = rSges(7,1) * t361 - rSges(7,2) * t362 - rSges(7,3) * t439;
t793 = rSges(7,1) * t650 - rSges(7,2) * t647;
t245 = -rSges(7,3) * t714 + t497 * t793;
t352 = pkin(5) * t497 - pkin(9) * t714;
t1197 = -t166 * t443 - t245 * t354 - t352 * t560;
t1194 = -t497 * t1124 - t681 * t714;
t523 = pkin(3) * t618 - qJ(4) * t619;
t886 = qJD(4) * t652;
t1187 = t523 * t892 - t619 * t886;
t1183 = t1214 * qJD(1);
t734 = qJD(2) * t563;
t330 = -t652 * t734 + (-t649 * t784 + t1005) * qJD(1);
t737 = qJD(2) * t565;
t332 = -t652 * t737 + (-t566 * t649 + t1018) * qJD(1);
t1182 = -qJD(2) * t1138 + t1220 * t619 + t1222 * t618 - t330 * t648 + t332 * t651 + t1200;
t1114 = qJD(1) * t411;
t331 = qJD(1) * t454 - t649 * t734;
t333 = qJD(1) * t456 - t649 * t737;
t1181 = qJD(1) * t1218 + qJD(2) * t1137 + t1219 * t619 - t1221 * t618 + t331 * t648 - t333 * t651 - t1114;
t1180 = (t1184 * t649 - t1185 * t652) * qJD(2);
t1179 = (-t1216 * t652 + t1228 * t649) * qJD(2);
t1178 = t1215 * qJD(1);
t1176 = qJD(1) * t1211 - t1213 * t649 + t1200;
t1175 = -t1114 - t1213 * t652 + (-t1253 * t649 - t1212 + t1250) * qJD(1);
t198 = rSges(7,3) * t440 - t439 * t793;
t200 = rSges(7,3) * t441 + t442 * t793;
t302 = -pkin(5) * t439 + pkin(9) * t440;
t306 = pkin(5) * t442 + pkin(9) * t441;
t1163 = -t198 * t353 - t200 * t354 - t302 * t559 - t306 * t560;
t244 = rSges(7,3) * t497 + t714 * t793;
t350 = -pkin(5) * t714 - pkin(9) * t497;
t1162 = t198 * t443 - t244 * t354 + t350 * t560;
t1161 = t200 * t443 + t244 * t353 - t350 * t559;
t1160 = 0.2e1 * qJD(2);
t347 = -rSges(6,1) * t714 + rSges(6,2) * t497;
t1157 = t347 * t559;
t1156 = t347 * t560;
t1046 = rSges(5,1) * t618;
t524 = -rSges(5,3) * t619 + t1046;
t470 = t524 * t649;
t1058 = pkin(2) * t648;
t827 = -t523 - t1058;
t806 = -t524 + t827;
t753 = t652 * t806;
t1153 = t1218 + t1229;
t794 = -rSges(6,1) * t440 + rSges(6,2) * t439;
t265 = rSges(6,3) * t652 + t794;
t1152 = qJD(1) * t265;
t527 = rSges(5,1) * t619 + rSges(5,3) * t618;
t639 = t652 * rSges(5,2);
t431 = t527 * t649 - t639;
t1151 = qJD(1) * t431;
t432 = rSges(4,1) * t960 - rSges(4,2) * t962 - t652 * rSges(4,3);
t634 = t649 * rSges(4,3);
t434 = rSges(4,1) * t959 - rSges(4,2) * t961 + t634;
t765 = t432 * t649 + t434 * t652;
t641 = t652 * pkin(7);
t586 = pkin(1) * t649 - t641;
t646 = -qJ(3) - pkin(7);
t610 = t652 * t646;
t1057 = pkin(2) * t651;
t612 = pkin(1) + t1057;
t904 = -t649 * t612 - t610;
t399 = t586 + t904;
t640 = t649 * pkin(7);
t587 = t652 * pkin(1) + t640;
t591 = t652 * t612;
t819 = -t646 * t649 + t591;
t400 = t819 - t587;
t935 = -t399 * t623 + t400 * t889;
t152 = qJD(2) * t765 + t935;
t934 = -t649 * t399 + t652 * t400;
t717 = t765 + t934;
t1150 = qJD(2) * t717 + t152;
t168 = -rSges(7,1) * t764 + t364 * rSges(7,2) + t442 * rSges(7,3);
t304 = -pkin(5) * t440 - pkin(9) * t439;
t308 = -t441 * pkin(5) + pkin(9) * t442;
t1054 = pkin(8) * t652;
t521 = pkin(4) * t960 + t1054;
t588 = pkin(4) * t959;
t522 = -pkin(8) * t649 + t588;
t990 = qJ(4) * t618;
t526 = pkin(3) * t619 + t990;
t472 = t526 * t649;
t476 = pkin(3) * t959 + qJ(4) * t961;
t888 = qJD(4) * t619;
t738 = t472 * t623 + t476 * t889 - t888 + t935;
t698 = t521 * t623 + t522 * t889 + t738;
t43 = t166 * t353 + t168 * t354 + t304 * t559 + t308 * t560 + t698;
t1056 = pkin(4) * t618;
t756 = t827 - t1056;
t702 = t756 * t889;
t570 = t618 * t886;
t621 = qJD(3) * t649;
t906 = t570 + t621;
t687 = t702 + t906;
t931 = t399 - t586;
t868 = -t472 + t931;
t816 = -t521 + t868;
t47 = (-t304 + t816) * qJD(1) + t687 + t1197;
t848 = t618 * t623;
t556 = pkin(4) * t848;
t887 = qJD(4) * t649;
t840 = t618 * t887;
t594 = t623 * t1058;
t902 = qJD(3) * t652 + t594;
t863 = -t523 * t623 - t902;
t707 = -t556 + t840 + t863;
t928 = t400 + t587;
t865 = t476 + t928;
t815 = t522 + t865;
t48 = t168 * t443 - t245 * t353 - t352 * t559 + (t308 + t815) * qJD(1) + t707;
t1149 = (t43 * (t166 * t441 - t168 * t440) + t47 * (t166 * t497 + t245 * t440) - t48 * (t168 * t497 + t245 * t441)) * qJD(6);
t1148 = -t1178 + t1179;
t1147 = t1180 + t1183;
t1145 = t1199 * t649 + t1201 * t652;
t1144 = -t1199 * t652 + t1201 * t649;
t1143 = qJD(2) * t1211 + t1219 * t618 + t1221 * t619 - t331 * t651 - t333 * t648;
t1142 = qJD(2) * t1212 + t1220 * t618 - t1222 * t619 + t330 * t651 + t332 * t648;
t1088 = t649 * (-t563 * t652 + t456) - t652 * (-Icges(3,2) * t956 + t455 - t600);
t699 = t453 * t652 - t454 * t649;
t1141 = -t1088 * t648 - t1202 * t618 + t1205 * t619 + t699 * t651;
t908 = t565 + t784;
t909 = -t563 + t566;
t1140 = (t1203 * t619 - t1204 * t618 - t648 * t908 + t651 * t909) * qJD(1);
t1139 = t974 + t1227;
t1136 = t1237 * qJD(1);
t1117 = -t526 - t527;
t498 = t527 * qJD(2);
t430 = qJD(2) * t526 - t888;
t873 = qJD(2) * t1057;
t802 = -t430 - t873;
t1129 = -qJD(2) * (-t1057 + t1117) - t498 + t802;
t1055 = pkin(4) * t619;
t801 = -t1055 - t1057;
t755 = -t526 + t801;
t1128 = -t755 * t889 + t1187;
t298 = -rSges(6,1) * t439 - rSges(6,2) * t440;
t300 = rSges(6,1) * t442 - rSges(6,2) * t441;
t1127 = t298 * t559 + t300 * t560;
t950 = t166 + t304;
t1122 = t43 * t950;
t946 = t245 + t352;
t1121 = t47 * t946;
t603 = qJD(4) * t618;
t1119 = (qJD(2) * t524 - t603) * t649;
t545 = qJD(1) * t586;
t1118 = qJD(1) * t399 - t545;
t790 = -qJD(1) * t472 + t1118 + t906;
t1116 = -qJD(1) * t521 + t702 + t790;
t125 = t217 * pkin(5) - pkin(9) * t216;
t126 = pkin(5) * t219 + pkin(9) * t218;
t713 = -t618 * t889 - t619 * t892;
t392 = pkin(4) * t713 - pkin(8) * t891;
t852 = t619 * t891;
t910 = -pkin(8) * t892 - t556;
t393 = pkin(4) * t852 + t910;
t847 = t619 * t889;
t549 = qJ(4) * t847;
t853 = t618 * t892;
t268 = pkin(3) * t713 - qJ(4) * t853 + t549 + t570;
t449 = t618 * t891 + t619 * t623;
t557 = pkin(3) * t848;
t269 = pkin(3) * t852 + qJ(4) * t449 - t557 + t840;
t1051 = pkin(1) - t612;
t617 = pkin(7) * t891;
t846 = t648 * t889;
t326 = -pkin(2) * t846 - t617 + t621 + (t1051 * t649 - t610) * qJD(1);
t858 = t646 * t892 + t902;
t327 = (-t1051 * t652 - t640) * qJD(1) - t858;
t871 = t326 * t889 + t327 * t623 - t399 * t609;
t877 = qJD(2) * qJD(4);
t721 = t268 * t889 + t269 * t623 + t472 * t609 + t618 * t877 + t871;
t929 = t400 + t476;
t866 = -t522 - t929;
t663 = t392 * t889 + t393 * t623 + t521 * t609 + t608 * t866 + t721;
t71 = t146 * rSges(7,1) + t145 * rSges(7,2) - t216 * rSges(7,3);
t72 = rSges(7,1) * t148 + rSges(7,2) * t147 + rSges(7,3) * t218;
t10 = t560 * t125 + t559 * t126 + t180 * t166 - t179 * t168 + t531 * t304 - t530 * t308 + t353 * t72 + t354 * t71 + t663;
t1034 = t126 + t72;
t1090 = t10 * t950 + t43 * t1034;
t1035 = t125 + t71;
t949 = t168 + t308;
t1089 = t10 * t949 + t43 * t1035;
t685 = t353 * (Icges(7,2) * t764 + t165 + t357) - t354 * (-Icges(7,2) * t361 + t163 - t356) + t443 * (t242 + (-Icges(7,2) * t650 - t1020) * t497);
t1083 = m(5) / 0.2e1;
t1082 = m(6) / 0.2e1;
t1081 = m(7) / 0.2e1;
t1080 = -t23 / 0.2e1;
t1079 = -pkin(3) - pkin(4);
t1078 = t179 / 0.2e1;
t1077 = t180 / 0.2e1;
t1075 = t353 / 0.2e1;
t1073 = -t354 / 0.2e1;
t1064 = t649 / 0.2e1;
t1063 = -t652 / 0.2e1;
t1061 = -rSges(5,1) - pkin(3);
t1050 = -pkin(8) - t646;
t25 = -t105 * t714 + t236 * t336 + t773 * t337 + (-t106 * t647 + t107 * t650 + (-t239 * t650 - t242 * t647) * qJD(6)) * t497;
t1049 = t25 * t443 + t88 * t885;
t1048 = rSges(3,1) * t651;
t1047 = rSges(4,1) * t619;
t1041 = t11 * t354;
t1040 = t12 * t353;
t1039 = t55 * t179;
t1038 = t56 * t180;
t636 = t649 * rSges(5,2);
t635 = t649 * rSges(3,3);
t1036 = -rSges(5,3) - qJ(4);
t525 = rSges(4,1) * t618 + rSges(4,2) * t619;
t826 = -t525 - t1058;
t798 = t652 * t826;
t752 = qJD(2) * t798;
t711 = t621 + t752;
t155 = (-t432 + t931) * qJD(1) + t711;
t989 = t155 * t525;
t901 = rSges(3,2) * t958 + t652 * rSges(3,3);
t477 = rSges(3,1) * t956 - t901;
t573 = rSges(3,1) * t648 + rSges(3,2) * t651;
t841 = t573 * t889;
t313 = -t841 + (-t477 - t586) * qJD(1);
t980 = t313 * t649;
t979 = t313 * t652;
t849 = t573 * t623;
t478 = rSges(3,1) * t955 - rSges(3,2) * t957 + t635;
t918 = t478 + t587;
t314 = qJD(1) * t918 - t849;
t508 = t573 * t652;
t978 = t314 * t508;
t349 = rSges(6,1) * t497 + rSges(6,2) * t714;
t977 = t349 * t560;
t108 = rSges(7,1) * t743 + rSges(7,2) * t744 + rSges(7,3) * t336;
t185 = pkin(5) * t337 + pkin(9) * t336;
t953 = t108 + t185;
t948 = t217 * rSges(6,1) + t216 * rSges(6,2);
t543 = t587 * qJD(1);
t941 = -t327 - t543;
t923 = -t441 * rSges(6,1) - t442 * rSges(6,2);
t922 = -t430 - t498;
t473 = t523 * t652;
t921 = -qJD(1) * t473 + t619 * t887;
t912 = rSges(5,2) * t891 + rSges(5,3) * t847;
t911 = rSges(4,2) * t853 + rSges(4,3) * t891;
t879 = qJD(1) * qJD(3);
t907 = qJD(1) * t594 + t652 * t879;
t903 = rSges(3,3) * t891 + t1223 * t892;
t900 = t649 ^ 2 + t652 ^ 2;
t890 = qJD(2) * t619;
t876 = -rSges(6,3) + t1050;
t875 = pkin(2) * t957;
t653 = qJD(2) ^ 2;
t874 = t653 * t1057;
t872 = -t269 + t941;
t870 = t652 * t326 + t649 * t327 - t399 * t891;
t529 = qJD(1) * (-pkin(1) * t892 + t617);
t869 = qJD(1) * t326 + t649 * t879 + t529;
t433 = rSges(5,1) * t959 + rSges(5,3) * t961 + t636;
t867 = -t433 - t929;
t469 = t523 * t649;
t864 = -t469 * t623 - t473 * t889 + t603;
t861 = t549 + t906;
t859 = t591 + t476;
t850 = t525 * t623;
t845 = t651 * t889;
t839 = -pkin(1) - t1048;
t836 = t619 * t877;
t831 = -t623 / 0.2e1;
t828 = t889 / 0.2e1;
t528 = -rSges(4,2) * t618 + t1047;
t825 = -t528 - t1057;
t817 = t649 * t472 + t652 * t476 + t934;
t814 = t523 * t608 + t652 * t836 + t907;
t813 = t557 + t858;
t812 = t588 + t859;
t807 = t900 * t1058;
t796 = t1048 - t1223;
t795 = rSges(6,1) * t219 - rSges(6,2) * t218;
t123 = -rSges(6,3) * t891 + t948;
t184 = rSges(6,1) * t337 - rSges(6,2) * t336;
t694 = -qJD(2) * t430 + t653 * t801;
t750 = t649 * t836 + t869 + (t268 + t570) * qJD(1);
t662 = qJD(1) * t392 + t649 * t694 + t750;
t39 = -t184 * t559 - t349 * t531 + (t123 + t702) * qJD(1) + t662;
t124 = -rSges(6,3) * t892 + t795;
t664 = qJD(1) * t556 + t652 * t694 + t814;
t708 = -t393 - t840 + t872;
t40 = -t184 * t560 + t349 * t530 + (-t124 + t708) * qJD(1) + t664;
t792 = t39 * t649 + t40 * t652;
t772 = -t314 * t649 - t979;
t751 = t652 * t268 + t649 * t269 + t472 * t891 + t870;
t748 = t649 * t521 + t652 * t522 + t817;
t746 = t813 - t910;
t499 = t528 * qJD(2);
t745 = -qJD(2) * t499 - t874;
t507 = t573 * t649;
t471 = t525 * t649;
t742 = -qJ(4) * t890 - t603;
t740 = -pkin(4) * t961 - t875;
t739 = -t349 + t756;
t309 = (t477 * t649 + t478 * t652) * qJD(2);
t722 = qJD(2) * t753;
t719 = t755 * t649;
t718 = -t526 - t1055;
t712 = qJD(2) * (t1079 * t618 - t1058);
t710 = t756 - t946;
t709 = -(-Icges(7,5) * t362 - Icges(7,6) * t361) * t354 + (Icges(7,5) * t364 + Icges(7,6) * t764) * t353 + (-Icges(7,5) * t647 - Icges(7,6) * t650) * t497 * t443;
t706 = -pkin(4) * t890 + t802;
t703 = -t612 + t718;
t697 = -t184 + t706;
t696 = t652 * t392 + t649 * t393 + t521 * t891 + t751;
t693 = -t1056 * t900 - t807;
t692 = t706 - t953;
t688 = -t612 + t1117;
t686 = (Icges(7,1) * t364 + t1021 - t162) * t353 - (-Icges(7,1) * t362 - t1022 + t161) * t354 + (-t239 + (-Icges(7,1) * t647 - t1019) * t497) * t443;
t535 = t796 * qJD(2);
t475 = t525 * t652;
t474 = t524 * t652;
t450 = t847 - t853;
t448 = t900 * t618 * qJD(2);
t335 = -qJD(2) * t507 + (t652 * t796 + t635) * qJD(1);
t334 = -rSges(3,2) * t845 + (-t651 * t892 - t846) * rSges(3,1) + t903;
t325 = (-rSges(7,1) * t647 - rSges(7,2) * t650) * t497;
t297 = -qJD(2) * t471 + (t528 * t652 + t634) * qJD(1);
t296 = -qJD(2) * t470 + (t527 * t652 + t636) * qJD(1);
t295 = rSges(4,1) * t713 - rSges(4,2) * t847 + t911;
t294 = rSges(5,1) * t713 - rSges(5,3) * t853 + t912;
t267 = -rSges(6,3) * t649 + t923;
t227 = rSges(7,1) * t364 + rSges(7,2) * t764;
t226 = -rSges(7,1) * t362 - rSges(7,2) * t361;
t178 = -t535 * t889 + (-t335 - t543 + t849) * qJD(1);
t177 = -t535 * t623 + t529 + (t334 - t841) * qJD(1);
t156 = -t850 + (t434 + t928) * qJD(1) - t902;
t130 = -t1119 + (t433 + t865) * qJD(1) + t863;
t129 = t722 + (-t431 + t868) * qJD(1) + t906;
t113 = (t431 * t649 + t433 * t652) * qJD(2) + t738;
t104 = t745 * t652 + (-t297 + t850 + t941) * qJD(1) + t907;
t103 = t745 * t649 + (t295 + t752) * qJD(1) + t869;
t92 = -t349 * t559 + (t267 + t815) * qJD(1) + t707;
t91 = -t977 + (-t265 + t816) * qJD(1) + t687;
t77 = t265 * t559 + t267 * t560 + t698;
t74 = (qJD(2) * t922 - t874) * t652 + (-t296 + t872 + t1119) * qJD(1) + t814;
t73 = -t649 * t874 + qJD(1) * t294 + (qJD(1) * t753 + t649 * t922) * qJD(2) + t750;
t44 = (t294 * t652 + t296 * t649 + (t431 * t652 + t649 * t867) * qJD(1)) * qJD(2) + t721;
t27 = t123 * t560 + t124 * t559 + t265 * t531 - t267 * t530 + t663;
t18 = -t166 * t885 - t354 * t108 + t179 * t245 - t560 * t185 + t530 * t352 - t443 * t72 + (-t126 + t708) * qJD(1) + t664;
t17 = t168 * t885 - t353 * t108 - t180 * t245 - t559 * t185 - t531 * t352 + t443 * t71 + (t125 + t702) * qJD(1) + t662;
t3 = [(-(-qJD(1) * t477 - t313 - t545 - t841) * t314 + t178 * (t649 * t839 + t641 + t901) + t177 * t918 + t314 * (t617 + t903) + (t573 * t980 - t978) * qJD(2) + ((-pkin(1) - t796) * t979 + (t313 * (-rSges(3,3) - pkin(7)) + t314 * t839) * t649) * qJD(1)) * m(3) + t23 * t1074 + (((t1153 * t652 - t1139 + t1184) * t652 + (t1153 * t649 + t1185 - t1230 - t390 + t799) * t649) * qJD(2) + t1148 + t1178) * t831 - t1125 * t1078 + (-(t1116 - t91 - t977 - t1152) * t92 + t40 * (-t794 + t904) + t91 * (t746 - t795) + t39 * (t812 + t923) + t92 * (t861 + t948) + (t40 * (-rSges(6,3) - pkin(8)) + t92 * t712) * t652 + (t39 * t876 + t40 * t718 + t742 * t91) * t649 + ((t91 * rSges(6,3) + t703 * t92) * t649 + (t703 * t91 + t876 * t92) * t652) * qJD(1)) * m(6) + t1038 / 0.2e1 - t1041 / 0.2e1 + t1040 / 0.2e1 + (-qJD(2) * t1226 + t1233 * t618 - t1234 * t619 + t182 * t714 + t183 * t497 - t336 * t343 + t337 * t346 + t533 * t651 + t534 * t648) * qJD(1) + t1049 + (((t169 - t1120 + t1139) * t649 + ((t1235 + t1246) * t652 + t1188 + t1225 + t1245) * t652) * qJD(2) + t1183) * t828 + (-(-qJD(1) * t432 + t1118 - t155 + t711) * t156 + t104 * (-t432 + t904) + t155 * t858 + t103 * (t434 + t819) + t156 * (t621 + t911) + (t156 * t798 + t649 * t989) * qJD(2) + ((-t155 * rSges(4,3) + t156 * (-t612 - t1047)) * t649 + (t155 * (-t528 - t612) - t156 * t646) * t652) * qJD(1)) * m(4) + ((t1137 - t1215) * t649 + (t1138 + t1214) * t652) * t880 / 0.2e1 + t36 * t1066 + t1039 / 0.2e1 + (-(-t129 + t722 + t790 - t1151) * t130 + t74 * (t639 + t904) + t129 * t813 + t73 * (t433 + t859) + t130 * (t861 + t912) + (t130 * (t1061 * t618 - t1058) * qJD(2) + (t129 * t688 - t130 * t646) * qJD(1)) * t652 + (-t73 * t646 + t74 * t1061 * t619 + (-t129 * qJD(4) + t1036 * t74) * t618 + t129 * (t1036 * t619 + t1046) * qJD(2) + (-t129 * rSges(5,2) + t130 * t688) * qJD(1)) * t649) * m(5) + (t18 * (t904 - t950 - t1054) + t17 * (t812 + t949) + (t1050 * t17 + t18 * t718) * t649 + (t742 * t649 + t703 * t891 - t1034 + t746) * t47 + (t652 * t712 + t1035 - t1116 + t47 + t861 + (t304 + (t1079 * t619 - t612 - t990) * t649 + t1050 * t652) * qJD(1) - t1197) * t48) * m(7) - (t38 + t36 + t41) * t560 / 0.2e1 + (t37 + t42) * t559 / 0.2e1 + (t116 + t132) * t531 / 0.2e1 + (t115 + t131) * t530 / 0.2e1 + (t1142 + t1145) * t623 / 0.2e1 - (-t1143 + t1144 + t1147) * t889 / 0.2e1 + t79 * t1077 + t354 * t1080 + t21 * t1073 + t20 * t1075; (0.2e1 * t309 * (t334 * t652 + t335 * t649 + (t477 * t652 - t478 * t649) * qJD(1)) + t772 * t535 + (-t177 * t649 - t178 * t652 + (-t314 * t652 + t980) * qJD(1)) * t573 - (t313 * t507 - t978) * qJD(1) - (t309 * (-t507 * t649 - t508 * t652) + t772 * t796) * qJD(2)) * m(3) + (t27 * t748 + (t27 * t267 + t40 * t739) * t652 + (t27 * t265 + t39 * t739) * t649 + (-t719 * qJD(2) + t1157 + t697 * t649 - t921 + (t652 * t739 - t300 - t740) * qJD(1)) * t92 + (t1156 + t697 * t652 + (t349 * t649 + t298 - t469) * qJD(1) + t1128) * t91 + (-t864 - t693 * qJD(2) + t696 + (t123 + t1152) * t652 + (t124 + (-t267 + t866) * qJD(1)) * t649 - t1127) * t77) * m(6) + (t1144 * qJD(1) + ((t1228 * qJD(1) + t1176 * t652) * t652 + (t1182 * t649 + t1216 * qJD(1) + (-t1175 + t1181) * t652) * t649) * t1160) * t1063 + (t44 * t817 + (t44 * t433 + t74 * t806) * t652 + (t44 * t431 + t73 * t806) * t649 + (-t921 + t1129 * t649 + (t474 + t875 + t753) * qJD(1)) * t130 + (-qJD(1) * t469 + t1129 * t652 + t1187) * t129 + (t751 + (t294 + t1151) * t652 + (qJD(1) * t867 + t296) * t649 - t864 - (-t470 * t649 - t474 * t652 - t807) * qJD(2)) * t113) * m(5) + (t1148 + t1179) * t892 / 0.2e1 + (t1147 + t1180) * t891 / 0.2e1 + (-t48 * (t1161 + t921) - t43 * (-t1163 + t864) - t1149 - t48 * (t306 + t740) * qJD(1) - (t43 * t693 + t48 * t719) * qJD(2) + t10 * t748 + t43 * t696 + (t18 * t710 + (t48 * t710 + t1122) * qJD(1) + t1089) * t652 + (t17 * t710 + t48 * t692 + (t1121 + t43 * (t866 - t949)) * qJD(1) + t1090) * t649 + (-(-t302 + t469) * qJD(1) + t692 * t652 + t1128 + t1162) * t47) * m(7) + (t104 * t798 - t155 * pkin(2) * t845 + (t295 * t889 + t297 * t623 + t871) * t717 + t152 * t870 + (-t155 * t499 + t152 * t295 + (t1150 * t432 + t156 * t826) * qJD(1)) * t652 + (t103 * t826 + t156 * (-t499 - t873) + t152 * t297 + (t989 + t1150 * (-t400 - t434)) * qJD(1)) * t649 - (t155 * t471 + t156 * (-t475 - t875)) * qJD(1) - (-t152 * t807 + (-t152 * t475 + t155 * t825) * t652 + (-t152 * t471 + t156 * t825) * t649) * qJD(2)) * m(4) + ((-t1134 * t889 - t1136) * t652 + ((t1135 * t652 + t1141) * qJD(2) + t1140) * t649) * t828 + ((-t1135 * t623 + t1136) * t649 + ((t1134 * t649 + t1141) * qJD(2) + t1140) * t652) * t831 + (t1143 * t652 + t1142 * t649 + (t1137 * t649 + t1138 * t652) * qJD(1)) * qJD(1) / 0.2e1 + ((t1088 * t651 + t1202 * t619 + t1205 * t618 + t699 * t648) * qJD(2) + (t1203 * t618 + t1204 * t619 + t648 * t909 + t651 * t908) * qJD(1) - t1194) * t1053 + (t1145 * qJD(1) + ((qJD(1) * t1184 + t1181 * t652) * t652 + (t1175 * t649 + t1185 * qJD(1) + (-t1176 + t1182) * t652) * t649) * t1160) * t1064 - t1224; 0.2e1 * (t1063 * t17 + t1064 * t18) * m(7) + 0.2e1 * (t1063 * t39 + t1064 * t40) * m(6) + 0.2e1 * (t1063 * t73 + t1064 * t74) * m(5) + 0.2e1 * (t103 * t1063 + t104 * t1064) * m(4); -m(5) * (t113 * t448 + t129 * t450 + t130 * t449) - m(6) * (t448 * t77 + t449 * t92 + t450 * t91) - m(7) * (t43 * t448 + t449 * t48 + t450 * t47) + 0.2e1 * ((t129 * t889 + t130 * t623 - t44) * t1083 + (t623 * t92 + t889 * t91 - t27) * t1082 + (t47 * t889 + t48 * t623 - t10) * t1081) * t619 + 0.2e1 * ((qJD(2) * t113 - t129 * t892 + t130 * t891 + t649 * t73 + t652 * t74) * t1083 + (qJD(2) * t77 + t891 * t92 - t892 * t91 + t792) * t1082 + (qJD(2) * t43 + t17 * t649 + t18 * t652 - t47 * t892 + t48 * t891) * t1081) * t618; t1194 * t1053 + (-t47 * (qJD(1) * t302 + t1162) - t48 * (-qJD(1) * t306 - t1161) - t43 * t1163 + t1149 + (t18 * t946 + t47 * t953 + (t48 * t946 - t1122) * qJD(1) - t1089) * t652 + (t17 * t946 + t48 * t953 + (t43 * t949 - t1121) * qJD(1) - t1090) * t649) * m(7) + (t27 * (-t265 * t649 - t267 * t652) + (t649 * t92 + t652 * t91) * t184 + ((-t649 * t91 + t652 * t92) * qJD(1) + t792) * t349 - t91 * (qJD(1) * t298 + t1156) - t92 * (-qJD(1) * t300 + t1157) + (-t123 * t652 - t124 * t649 + (-t265 * t652 + t267 * t649) * qJD(1) + t1127) * t77) * m(6) + t1224; t216 * t1080 + t442 * t1 / 0.2e1 + (-t439 * t53 + t442 * t54 - t714 * t79) * t1077 + (-t13 * t439 + t14 * t442 - t20 * t714 - t216 * t54 + t218 * t53 + t336 * t79) * t1075 + t218 * t1123 - t439 * t2 / 0.2e1 + (t1125 * t714 - t439 * t51 + t442 * t52) * t1078 + (-t1125 * t336 - t15 * t439 + t16 * t442 - t21 * t714 - t216 * t52 + t218 * t51) * t1073 + t336 * t26 / 0.2e1 - t714 * (t1038 + t1039 + t1040 - t1041 + t1049) / 0.2e1 + (-t439 * t55 + t442 * t56 - t714 * t88) * t885 / 0.2e1 + t443 * (-t11 * t439 + t12 * t442 - t216 * t56 + t218 * t55 - t25 * t714 + t336 * t88) / 0.2e1 + (t364 * t685 + t442 * t709 - t686 * t764) * t1076 + (t361 * t686 - t362 * t685 - t439 * t709) * t1074 + (-t709 * t714 + (-t647 * t685 + t686 * t650) * t497) * t1072 + (t18 * (t166 * t714 - t245 * t439) + t17 * (-t168 * t714 - t245 * t442) + t10 * (t166 * t442 + t168 * t439) + (-t108 * t442 + t168 * t336 + t216 * t245 - t227 * t443 + t325 * t353 - t71 * t714) * t48 + (-t108 * t439 - t166 * t336 + t218 * t245 + t226 * t443 + t325 * t354 + t714 * t72) * t47 + (-t166 * t216 - t168 * t218 - t226 * t353 - t227 * t354 + t439 * t71 + t442 * t72) * t43) * m(7);];
tauc  = t3(:);