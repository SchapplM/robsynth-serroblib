% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:58
% DurationCPUTime: 10.06s
% Computational Cost: add. (5620->389), mult. (7708->506), div. (0->0), fcn. (5944->6), ass. (0->225)
t460 = Icges(5,4) + Icges(4,5);
t221 = sin(qJ(1));
t222 = cos(qJ(1));
t217 = pkin(6) + qJ(3);
t198 = sin(t217);
t188 = Icges(5,5) * t198;
t199 = cos(t217);
t265 = Icges(5,1) * t199 + t188;
t459 = -t221 * t265 + t460 * t222;
t335 = t198 * t221;
t176 = Icges(4,4) * t335;
t333 = t199 * t221;
t449 = Icges(4,1) * t333 - t176 - t459;
t263 = Icges(5,3) * t199 - t188;
t353 = Icges(4,4) * t198;
t458 = Icges(4,2) * t199 + t263 + t353;
t350 = Icges(5,5) * t199;
t152 = Icges(5,1) * t198 - t350;
t189 = Icges(4,4) * t199;
t457 = Icges(4,1) * t198 + t152 + t189;
t149 = Icges(5,4) * t199 + Icges(5,6) * t198;
t101 = Icges(5,2) * t221 + t149 * t222;
t147 = Icges(4,5) * t199 - Icges(4,6) * t198;
t99 = Icges(4,3) * t221 + t147 * t222;
t456 = t101 + t99;
t105 = Icges(5,4) * t221 + t222 * t265;
t155 = Icges(4,1) * t199 - t353;
t107 = Icges(4,5) * t221 + t155 * t222;
t448 = t105 + t107;
t145 = Icges(5,3) * t198 + t350;
t264 = -Icges(4,2) * t198 + t189;
t455 = t145 - t264;
t447 = (Icges(4,6) - Icges(5,6)) * t199 + t460 * t198;
t454 = t155 + t265;
t347 = Icges(4,6) * t222;
t102 = Icges(4,4) * t333 - Icges(4,2) * t335 - t347;
t341 = t102 * t198;
t96 = -Icges(5,6) * t222 + t145 * t221;
t424 = -t198 * t96 - t449 * t199 + t341;
t344 = Icges(4,3) * t222;
t98 = Icges(4,5) * t333 - Icges(4,6) * t335 - t344;
t453 = -t424 * t221 - t222 * t98;
t452 = t102 - t96;
t103 = Icges(4,6) * t221 + t222 * t264;
t332 = t199 * t222;
t175 = Icges(5,5) * t332;
t334 = t198 * t222;
t346 = Icges(5,6) * t221;
t97 = Icges(5,3) * t334 + t175 + t346;
t451 = t103 - t97;
t446 = t458 * qJD(3);
t445 = t457 * qJD(3);
t437 = -t198 * t458 + t457 * t199;
t444 = t456 * t221 + t448 * t332 + t97 * t334;
t100 = -Icges(5,2) * t222 + t149 * t221;
t91 = t221 * t100;
t443 = -t221 * t98 - t449 * t332 - t96 * t334 - t91;
t442 = t455 * qJD(3);
t441 = t454 * qJD(3);
t440 = -t147 - t149;
t389 = t447 * t222;
t390 = t447 * t221;
t331 = t222 * t100;
t417 = -t331 + t453;
t416 = -t102 * t334 - t443;
t415 = -t103 * t334 + t444;
t436 = t221 * t437 - t389;
t435 = t222 * t437 + t390;
t397 = t449 * t198 + t199 * t452;
t396 = t198 * t448 + t199 * t451;
t434 = t446 * t222 + (t221 * t264 - t347 - t96) * qJD(1);
t433 = t446 * t221 + (t145 * t222 - t103 + t346) * qJD(1);
t432 = -t445 * t222 + (-t155 * t221 + t459) * qJD(1);
t431 = -qJD(1) * t448 + t221 * t445;
t275 = t101 * t222 - t105 * t333 - t97 * t335;
t81 = t107 * t333;
t280 = t222 * t99 - t81;
t36 = -t103 * t335 - t280;
t430 = -t275 + t36;
t429 = t447 * qJD(3);
t340 = t103 * t198;
t428 = t198 * t97 + t199 * t448 - t340;
t427 = (-t222 * t458 + t448) * t221 - (-Icges(4,2) * t333 - t263 * t221 - t176 + t449) * t222;
t413 = rSges(5,3) + qJ(4);
t426 = t441 * t199 + t442 * t198 + (-t198 * t457 - t199 * t458) * qJD(3) + t447 * qJD(1);
t425 = t452 * t222 + (-Icges(5,1) * t334 + t152 * t222 + t175 - t451) * t221;
t423 = t456 * qJD(1);
t422 = t457 - t455;
t421 = -t458 + t454;
t418 = t437 * qJD(1) + t440 * qJD(3);
t373 = rSges(5,1) + pkin(3);
t414 = t435 * qJD(1);
t156 = pkin(3) * t198 - qJ(4) * t199;
t157 = rSges(5,1) * t198 - rSges(5,3) * t199;
t314 = t156 + t157;
t160 = rSges(5,1) * t199 + rSges(5,3) * t198;
t412 = pkin(3) * t199 + qJ(4) * t198 + t160;
t411 = -t396 * qJD(3) + t434 * t198 + t432 * t199 + t423;
t382 = qJD(1) * t100;
t410 = -qJD(1) * t98 + t397 * qJD(3) - t433 * t198 + t431 * t199 - t382;
t409 = (t415 * t221 - t416 * t222) * qJD(3);
t408 = (t430 * t221 - t417 * t222) * qJD(3);
t407 = t436 * qJD(1);
t406 = t424 * qJD(1) - t429 * t221 + t423;
t405 = -t382 - t429 * t222 + (-t147 * t221 + t344 - t428) * qJD(1);
t404 = 0.2e1 * qJD(3);
t375 = t221 / 0.2e1;
t403 = t407 + t408;
t402 = t409 + t414;
t401 = -t418 * t221 + t426 * t222;
t400 = t426 * t221 + t418 * t222;
t399 = t424 * qJD(3) + t431 * t198 + t433 * t199;
t398 = t428 * qJD(3) + t432 * t198 - t434 * t199;
t203 = t222 * qJ(2);
t168 = pkin(1) * t221 - t203;
t163 = qJD(1) * t168;
t219 = cos(pkin(6));
t190 = pkin(2) * t219 + pkin(1);
t220 = -pkin(5) - qJ(2);
t193 = t222 * t220;
t309 = -t221 * t190 - t193;
t94 = t168 + t309;
t395 = qJD(1) * t94 - t163;
t296 = qJD(4) * t199;
t212 = t221 * rSges(5,2);
t322 = t332 * t373 + t413 * t334 + t212;
t215 = t222 * rSges(5,2);
t323 = t221 * t412 - t215;
t30 = -t296 + (t221 * t323 + t222 * t322) * qJD(3);
t394 = qJD(3) * t30;
t393 = t198 * t373;
t276 = t314 * qJD(3);
t297 = qJD(4) * t198;
t244 = -t276 + t297;
t392 = t244 * t221;
t210 = t221 * rSges(4,3);
t112 = rSges(4,1) * t332 - rSges(4,2) * t334 + t210;
t178 = t222 * t190;
t277 = -t220 * t221 + t178;
t391 = t112 + t277;
t298 = qJD(3) * t222;
t288 = t199 * t298;
t300 = qJD(1) * t222;
t388 = rSges(5,2) * t300 + t413 * t288;
t202 = t221 * qJ(2);
t170 = t222 * pkin(1) + t202;
t363 = rSges(3,2) * sin(pkin(6));
t365 = rSges(3,1) * t219;
t257 = t221 * rSges(3,3) + (-t363 + t365) * t222;
t387 = t170 + t257;
t386 = -t427 * t198 + t425 * t199;
t385 = (-t422 * t198 + t421 * t199) * qJD(1);
t384 = t331 + t444;
t383 = t440 * qJD(1);
t374 = -t222 / 0.2e1;
t371 = qJD(1) / 0.2e1;
t370 = pkin(1) - t190;
t167 = t222 * t297;
t301 = qJD(1) * t221;
t242 = -t198 * t298 - t199 * t301;
t291 = t198 * t301;
t369 = t242 * t373 - t413 * t291 + t167 + t388;
t127 = t157 * t221;
t299 = qJD(3) * t221;
t368 = (pkin(3) * t300 + qJ(4) * t299) * t199 + (qJ(4) * t300 + (-pkin(3) * qJD(3) + qJD(4)) * t221) * t198 - qJD(3) * t127 + (t160 * t222 + t212) * qJD(1);
t364 = rSges(4,1) * t199;
t158 = rSges(4,1) * t198 + rSges(4,2) * t199;
t132 = t158 * t222;
t201 = qJD(2) * t222;
t290 = t158 * t299;
t42 = qJD(1) * t391 - t201 - t290;
t362 = t132 * t42;
t110 = rSges(4,1) * t333 - rSges(4,2) * t335 - t222 * rSges(4,3);
t200 = qJD(2) * t221;
t289 = t158 * t298;
t272 = t200 - t289;
t356 = t94 - t168;
t41 = (-t110 + t356) * qJD(1) + t272;
t361 = t221 * t41;
t141 = qJD(1) * t170 - t201;
t185 = t220 * t301;
t358 = -t141 + t185 - (-t222 * t370 - t202) * qJD(1);
t324 = -qJD(3) * t412 + t296;
t321 = -t156 * t221 - t127;
t320 = t314 * t222;
t295 = qJD(1) * qJD(2);
t194 = qJ(2) * t300;
t305 = t194 + t200;
t319 = qJD(1) * (-pkin(1) * t301 + t305) + t221 * t295;
t311 = rSges(4,2) * t291 + rSges(4,3) * t300;
t310 = t167 + t200;
t186 = t221 * t363;
t308 = rSges(3,3) * t300 + qJD(1) * t186;
t307 = t185 + t201;
t306 = t222 * rSges(3,3) + t186;
t294 = t221 * t365;
t292 = qJD(1) * (-t194 + (t221 * t370 - t193) * qJD(1)) + t319;
t287 = -pkin(1) - t365;
t284 = -t299 / 0.2e1;
t281 = t298 / 0.2e1;
t278 = -t98 + t340;
t270 = -rSges(4,2) * t198 + t364;
t269 = -t221 * t42 - t222 * t41;
t256 = qJD(3) * (t296 + t324);
t128 = t158 * t221;
t53 = (t110 * t221 + t112 * t222) * qJD(3);
t243 = -t222 * t276 + t310;
t237 = -t190 - t412;
t192 = t222 * t295;
t143 = t270 * qJD(3);
t113 = t294 - t306;
t77 = qJD(1) * t387 - t201;
t76 = t200 + (-t113 - t168) * qJD(1);
t73 = -qJD(3) * t128 + (t222 * t270 + t210) * qJD(1);
t71 = rSges(4,1) * t242 - rSges(4,2) * t288 + t311;
t57 = t192 + (-qJD(1) * t257 - t141) * qJD(1);
t56 = qJD(1) * (-qJD(1) * t294 + t308) + t319;
t29 = -t201 + t392 + (t322 + t277) * qJD(1);
t28 = (-t323 + t356) * qJD(1) + t243;
t25 = -t143 * t298 + t192 + (-t73 + t290 + t358) * qJD(1);
t24 = -t143 * t299 + (t71 - t289) * qJD(1) + t292;
t11 = t192 + t222 * t256 + (t358 - t368 - t392) * qJD(1);
t10 = t221 * t256 + (t222 * t244 + t369) * qJD(1) + t292;
t1 = (t297 + t369 * t222 + t368 * t221 + (-t221 * t322 + t222 * t323) * qJD(1)) * qJD(3);
t2 = [(((t36 - t81 + (t99 + t341) * t222 + t443) * t222 + (t384 + t417 - t453) * t221) * qJD(3) + t414) * t281 + (t437 * qJD(3) + t441 * t198 - t442 * t199) * qJD(1) + (t11 * (t215 + t309) + t28 * t307 + t10 * (t178 + t322) + t29 * (t310 + t388) + (-t29 * qJD(3) * t393 + (-t29 * t220 + t237 * t28) * qJD(1)) * t222 + (-t10 * t220 - t11 * t373 * t199 + (-t28 * qJD(4) - t11 * t413) * t198 + t28 * (-t199 * t413 + t393) * qJD(3) + (-t28 * rSges(5,2) + t237 * t29) * qJD(1)) * t221 - (-qJD(1) * t323 + t243 - t28 + t395) * t29) * m(5) + (t25 * (-t110 + t309) + t41 * t307 + t24 * t391 + t42 * (t200 + t311) + (t158 * t361 - t362) * qJD(3) + ((-t41 * rSges(4,3) + t42 * (-t190 - t364)) * t221 + (t41 * (-t190 - t270) - t42 * t220) * t222) * qJD(1) - (-qJD(1) * t110 + t272 + t395 - t41) * t42) * m(4) + (t57 * (t221 * t287 + t203 + t306) + t76 * t201 + t56 * t387 + t77 * (t305 + t308) + (t76 * (t287 + t363) * t222 + (t76 * (-rSges(3,3) - qJ(2)) + t77 * t287) * t221) * qJD(1) - (-qJD(1) * t113 - t163 + t200 - t76) * t77) * m(3) + (((t222 * t278 - t384 + t415) * t222 + (t221 * t278 + t275 + t280 + t416 - t91) * t221) * qJD(3) + t403 - t407) * t284 + (t398 + t401) * t299 / 0.2e1 - (-t399 + t400 + t402) * t298 / 0.2e1 + ((t397 + t436) * t221 + (t396 + t435) * t222) * qJD(3) * t371; 0.2e1 * (t10 * t374 + t11 * t375) * m(5) + 0.2e1 * (t24 * t374 + t25 * t375) * m(4) + 0.2e1 * (t56 * t374 + t375 * t57) * m(3); -((t425 * t198 + t427 * t199) * qJD(3) + (t421 * t198 + t422 * t199) * qJD(1)) * qJD(1) / 0.2e1 + (t399 * t222 + t398 * t221 + (t221 * t397 + t222 * t396) * qJD(1)) * t371 + ((-t299 * t389 - t383) * t221 + ((t221 * t390 + t386) * qJD(3) + t385) * t222) * t284 + ((-t298 * t390 + t383) * t222 + ((t222 * t389 + t386) * qJD(3) + t385) * t221) * t281 + ((-t11 * t314 + t28 * t324 + t1 * t322 + t30 * t369 + (-t29 * t314 + t30 * t323) * qJD(1)) * t222 + (-t10 * t314 + t29 * t324 + t1 * t323 + t30 * t368 + (t28 * t314 - t30 * t322) * qJD(1)) * t221 - (t198 * t30 + (t221 * t29 + t222 * t28) * t199) * qJD(4) - (-t28 * t321 - t29 * t320) * qJD(1) - ((-t28 * t412 - t30 * t320) * t222 + (-t29 * t412 + t30 * t321) * t221) * qJD(3)) * m(5) + (-(t128 * t41 - t362) * qJD(1) - (t53 * (-t128 * t221 - t132 * t222) + t269 * t270) * qJD(3) + 0.2e1 * t53 * (t221 * t73 + t222 * t71 + (t110 * t222 - t112 * t221) * qJD(1)) + t269 * t143 + (-t24 * t221 - t25 * t222 + (-t222 * t42 + t361) * qJD(1)) * t158) * m(4) + (t401 * qJD(1) + ((t415 * qJD(1) + t410 * t222) * t222 + (t405 * t221 + t416 * qJD(1) + (-t406 + t411) * t222) * t221) * t404) * t375 + (t400 * qJD(1) + ((t430 * qJD(1) + t406 * t222) * t222 + (t411 * t221 + t417 * qJD(1) + (-t405 + t410) * t222) * t221) * t404) * t374 + (t403 + t408) * t301 / 0.2e1 + (t402 + t409) * t300 / 0.2e1; (-t1 * t199 + 0.2e1 * (t394 / 0.2e1 + t10 * t375 + t11 * t222 / 0.2e1 - (t221 ^ 2 + t222 ^ 2) * t394 / 0.2e1) * t198) * m(5);];
tauc = t2(:);
