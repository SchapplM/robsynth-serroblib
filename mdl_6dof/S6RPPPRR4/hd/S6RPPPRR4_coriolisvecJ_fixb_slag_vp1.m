% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:28
% EndTime: 2019-03-09 01:35:50
% DurationCPUTime: 22.98s
% Computational Cost: add. (14708->804), mult. (34991->1048), div. (0->0), fcn. (39159->8), ass. (0->380)
t296 = sin(qJ(5));
t298 = cos(qJ(5));
t295 = sin(qJ(6));
t297 = cos(qJ(6));
t484 = sin(pkin(9));
t485 = cos(pkin(9));
t509 = sin(qJ(1));
t510 = cos(qJ(1));
t249 = -t484 * t509 - t485 * t510;
t250 = t484 * t510 - t485 * t509;
t457 = t295 * t296;
t170 = t249 * t457 - t250 * t297;
t462 = t249 * t298;
t456 = t296 * t297;
t167 = t249 * t456 + t250 * t295;
t480 = Icges(7,4) * t167;
t94 = -Icges(7,2) * t170 - Icges(7,6) * t462 + t480;
t163 = Icges(7,4) * t170;
t97 = Icges(7,1) * t167 - Icges(7,5) * t462 - t163;
t566 = t295 * t94 - t297 * t97;
t91 = Icges(7,5) * t167 - Icges(7,6) * t170 - Icges(7,3) * t462;
t36 = t296 * t91 - t298 * t566;
t168 = -t249 * t297 - t250 * t457;
t458 = t250 * t298;
t169 = -t249 * t295 + t250 * t456;
t481 = Icges(7,4) * t169;
t95 = Icges(7,2) * t168 - Icges(7,6) * t458 + t481;
t162 = Icges(7,4) * t168;
t98 = Icges(7,1) * t169 - Icges(7,5) * t458 + t162;
t501 = -t167 * t98 + t170 * t95;
t502 = -t168 * t94 - t169 * t97;
t92 = Icges(7,5) * t169 + Icges(7,6) * t168 - Icges(7,3) * t458;
t569 = -t502 - t298 * (t249 * t92 + t250 * t91) - t501;
t503 = t168 * t95 + t169 * t98;
t546 = -t167 * t97 + t170 * t94;
t568 = t546 + t503 + (t249 * t91 - t250 * t92) * t298;
t425 = qJD(6) * t296;
t276 = qJD(1) - t425;
t479 = Icges(7,4) * t295;
t275 = t298 * t479;
t455 = t297 * t298;
t476 = Icges(7,5) * t296;
t205 = -Icges(7,1) * t455 + t275 - t476;
t478 = Icges(7,4) * t297;
t361 = -Icges(7,2) * t295 + t478;
t531 = Icges(7,6) * t296 + t298 * t361;
t359 = Icges(7,5) * t297 - Icges(7,6) * t295;
t532 = Icges(7,3) * t296 + t298 * t359;
t545 = t167 * t205 + t170 * t531 + t462 * t532;
t564 = t545 * t276;
t377 = rSges(7,1) * t167 - rSges(7,2) * t170;
t102 = rSges(7,3) * t462 - t377;
t424 = qJD(6) * t298;
t430 = qJD(5) * t250;
t186 = t249 * t424 - t430;
t555 = pkin(5) * t298;
t383 = pkin(8) * t296 + t555;
t428 = qJD(5) * t383;
t235 = qJD(4) * t250;
t285 = qJD(2) * t509;
t440 = t235 + t285;
t376 = rSges(7,1) * t297 - rSges(7,2) * t295;
t529 = rSges(7,3) * t296 + t298 * t376;
t562 = -t102 * t276 - t186 * t529 + t250 * t428 + t440;
t259 = Icges(6,5) * t296 + Icges(6,6) * t298;
t144 = Icges(6,3) * t250 + t249 * t259;
t469 = t144 * t250;
t289 = Icges(6,4) * t296;
t362 = Icges(6,2) * t298 + t289;
t147 = Icges(6,6) * t250 + t249 * t362;
t560 = t147 * t296;
t559 = t147 * t298;
t431 = qJD(5) * t249;
t187 = -t250 * t424 - t431;
t200 = -Icges(7,3) * t298 + t296 * t359;
t354 = -t205 * t297 - t295 * t531;
t369 = t295 * t95 - t297 * t98;
t525 = t186 * (-t249 * t532 - t566) + t187 * (t250 * t532 + t369) + t276 * (t200 + t354);
t558 = t525 * t298;
t247 = t249 * pkin(3);
t326 = pkin(1) * t510 + qJ(2) * t509;
t258 = qJD(1) * t326;
t286 = qJD(2) * t510;
t402 = qJD(1) * t510;
t434 = -pkin(2) * t402 + t286;
t408 = t258 - t434;
t423 = t249 * qJD(4);
t470 = qJ(4) * t250;
t557 = t423 - qJD(1) * (-t247 + t470) - t408;
t360 = Icges(6,5) * t298 - Icges(6,6) * t296;
t172 = t360 * t250;
t482 = Icges(6,4) * t298;
t262 = Icges(6,2) * t296 - t482;
t264 = -Icges(6,1) * t298 + t289;
t352 = t262 * t298 + t264 * t296;
t533 = t249 * t352 - t172;
t556 = t533 * qJD(1);
t433 = qJD(1) * t249;
t554 = pkin(7) * t433;
t379 = rSges(6,1) * t298 - rSges(6,2) * t296;
t429 = qJD(5) * t379;
t553 = t249 * t429;
t173 = t360 * t249;
t229 = -t286 + t258;
t426 = qJD(5) * t298;
t341 = t250 * t426 + t296 * t433;
t552 = t554 + t557;
t288 = t510 * qJ(2);
t419 = t509 * pkin(1);
t266 = t419 - t288;
t257 = qJD(1) * t266;
t435 = t285 - t257;
t436 = qJ(2) * t402 + t285;
t551 = t436 - t435;
t263 = Icges(6,1) * t296 + t482;
t437 = t262 - t263;
t438 = t264 + t362;
t550 = (t296 * t437 - t298 * t438) * qJD(1);
t544 = 0.2e1 * qJD(5);
t272 = pkin(5) * t296 - pkin(8) * t298;
t543 = t272 * t249;
t237 = t249 * qJ(4);
t507 = t250 * pkin(3);
t541 = t237 + t507;
t542 = t541 - t266;
t418 = t509 * pkin(2);
t318 = rSges(4,1) * t250 - rSges(4,2) * t249 - t418;
t540 = -t266 + t318;
t327 = rSges(3,1) * t510 + rSges(3,3) * t509;
t539 = t326 + t327;
t446 = rSges(7,1) * t169 + rSges(7,2) * t168;
t101 = -rSges(7,3) * t458 + t446;
t538 = -t101 * t276 - t187 * t529 + t249 * t428 + t552;
t459 = t250 * t296;
t154 = rSges(6,1) * t459 + rSges(6,2) * t458 - rSges(6,3) * t249;
t537 = -qJD(1) * t154 + t552 + t553;
t192 = rSges(5,2) * t249 + rSges(5,3) * t250;
t536 = -qJD(1) * t192 + t557;
t378 = rSges(6,1) * t296 + rSges(6,2) * t298;
t153 = rSges(6,3) * t250 + t378 * t249;
t363 = Icges(7,1) * t297 - t479;
t530 = t298 * t363 + t476;
t252 = t362 * qJD(5);
t253 = t263 * qJD(5);
t528 = qJD(5) * (t262 * t296 - t264 * t298) - t252 * t298 - t253 * t296;
t427 = qJD(5) * t296;
t406 = t250 * t427;
t409 = pkin(5) * t341 + pkin(8) * t406;
t467 = t433 * t298;
t104 = -pkin(8) * t467 + t409;
t233 = t250 * qJD(1);
t403 = t250 * t425;
t137 = -t433 * t424 + (t233 + t403) * qJD(5);
t206 = -rSges(7,3) * t298 + t296 * t376;
t248 = (rSges(7,1) * t295 + rSges(7,2) * t297) * t298;
t161 = qJD(5) * t206 + qJD(6) * t248;
t299 = qJD(1) ^ 2;
t401 = qJD(1) * t509;
t443 = (-pkin(1) * t401 + t285 + t436) * qJD(1);
t335 = -t299 * t418 + t443;
t392 = pkin(3) * t233 + qJ(4) * t433;
t337 = t392 + t235;
t315 = qJD(1) * t337 + qJD(4) * t233 + t335;
t256 = qJD(5) * t272;
t464 = t249 * t256;
t312 = -qJD(6) * t249 + t341;
t381 = t233 - t403;
t79 = -t295 * t312 + t297 * t381;
t80 = t295 * t381 + t297 * t312;
t420 = rSges(7,1) * t80 + rSges(7,2) * t79 + rSges(7,3) * t406;
t47 = -rSges(7,3) * t467 + t420;
t508 = t233 * pkin(7);
t12 = t137 * t529 - t187 * t161 + t276 * t47 + (t104 + t508) * qJD(1) + (-t101 * t424 + t233 * t383 + t464) * qJD(5) + t315;
t465 = t233 * t298;
t338 = t249 * t427 + t465;
t466 = t233 * t296;
t339 = -t249 * t426 + t466;
t105 = pkin(5) * t339 - pkin(8) * t338;
t404 = t249 * t425;
t136 = -t233 * t424 + (-t433 - t404) * qJD(5);
t281 = qJD(1) * t286;
t293 = t510 * pkin(2);
t350 = -t293 * t299 + t281;
t334 = qJD(4) * t433 + t350;
t499 = pkin(7) * qJD(1);
t323 = t433 * t499 + t334;
t143 = -pkin(3) * t433 + qJ(4) * t233 - t423;
t452 = -t143 - t229;
t461 = t250 * t256;
t311 = qJD(6) * t250 - t339;
t382 = -t433 + t404;
t81 = t295 * t311 + t297 * t382;
t82 = t295 * t382 - t297 * t311;
t384 = rSges(7,1) * t82 + rSges(7,2) * t81;
t48 = -rSges(7,3) * t338 + t384;
t13 = -t136 * t529 + t186 * t161 - t276 * t48 + (t102 * t424 + t383 * t433 - t461) * qJD(5) + (-t105 + t452) * qJD(1) + t323;
t228 = pkin(5) * t459;
t181 = -pkin(8) * t458 + t228;
t34 = qJD(1) * t181 - t538;
t527 = -t12 * t250 - t13 * t249 - t34 * t433;
t240 = (Icges(7,1) * t295 + t478) * t298;
t303 = t186 * (Icges(7,1) * t170 + t480 + t94) + t187 * (Icges(7,1) * t168 - t481 - t95) - t276 * (-t531 - t240);
t23 = -t458 * t92 + t503;
t24 = t458 * t91 + t502;
t60 = -t168 * t531 + t169 * t205 + t458 * t532;
t59 = t60 * t276;
t10 = t186 * t24 + t187 * t23 + t59;
t524 = -t10 / 0.2e1;
t523 = t136 / 0.2e1;
t522 = t137 / 0.2e1;
t521 = -t186 / 0.2e1;
t520 = t186 / 0.2e1;
t519 = -t187 / 0.2e1;
t518 = t187 / 0.2e1;
t513 = -t276 / 0.2e1;
t512 = t276 / 0.2e1;
t511 = rSges(7,3) + pkin(8);
t246 = t250 * pkin(7);
t340 = t406 - t467;
t41 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t340;
t43 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t340;
t45 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t340;
t7 = (-qJD(5) * t369 - t41) * t296 + (-qJD(5) * t92 + t295 * t43 - t297 * t45 + (t295 * t98 + t297 * t95) * qJD(6)) * t298;
t506 = t7 * t187;
t42 = Icges(7,5) * t82 + Icges(7,6) * t81 - Icges(7,3) * t338;
t44 = Icges(7,4) * t82 + Icges(7,2) * t81 - Icges(7,6) * t338;
t46 = Icges(7,1) * t82 + Icges(7,4) * t81 - Icges(7,5) * t338;
t8 = (qJD(5) * t566 - t42) * t296 + (qJD(5) * t91 + t295 * t44 - t297 * t46 + (-t295 * t97 - t297 * t94) * qJD(6)) * t298;
t505 = t8 * t186;
t497 = t433 * rSges(6,3);
t231 = t250 * t499;
t349 = -t418 + t542;
t33 = t231 + (t543 + t349) * qJD(1) + t562;
t495 = t233 * t33;
t494 = t249 * rSges(5,3);
t218 = Icges(6,4) * t462;
t463 = t249 * t296;
t477 = Icges(6,5) * t250;
t152 = -Icges(6,1) * t463 - t218 - t477;
t52 = t144 * t249 - t147 * t458 + t152 * t459;
t491 = t250 * t52;
t35 = -t296 * t92 + t298 * t369;
t489 = t35 * t137;
t488 = t36 * t136;
t454 = t101 + t181;
t453 = t102 - t543;
t148 = -Icges(6,6) * t249 + t250 * t362;
t451 = t250 * t264 + t148;
t450 = -t249 * t264 - t147;
t217 = Icges(6,4) * t458;
t151 = Icges(6,1) * t459 - Icges(6,5) * t249 + t217;
t449 = -Icges(6,2) * t459 + t151 + t217;
t448 = Icges(6,2) * t463 + t152 - t218;
t447 = t161 + t256;
t444 = -t529 - t383;
t442 = rSges(4,1) * t233 - rSges(4,2) * t433;
t441 = -rSges(5,2) * t233 + rSges(5,3) * t433;
t439 = -rSges(4,1) * t249 - rSges(4,2) * t250;
t432 = qJD(1) * t259;
t422 = -t510 / 0.2e1;
t421 = t509 / 0.2e1;
t416 = t33 * t511;
t415 = t509 * rSges(3,1);
t145 = -Icges(6,3) * t249 + t250 * t259;
t53 = -t250 * t145 - t148 * t462 - t151 * t463;
t410 = t145 * t249 - t148 * t458 - t151 * t459;
t407 = t293 + t326;
t400 = qJD(5) * t424;
t397 = -t431 / 0.2e1;
t396 = t431 / 0.2e1;
t394 = t430 / 0.2e1;
t393 = -t426 / 0.2e1;
t190 = qJD(1) * t541;
t391 = t190 + t235 + t435;
t389 = rSges(6,1) * t341 + rSges(6,2) * t467 + rSges(6,3) * t233;
t388 = -t247 + t407;
t387 = qJD(6) * t393;
t380 = rSges(4,1) * t433 + rSges(4,2) * t233;
t375 = -rSges(5,2) * t433 - rSges(5,3) * t233;
t374 = t23 * t250 - t24 * t249;
t25 = t462 * t92 + t501;
t26 = -t462 * t91 - t546;
t373 = t249 * t26 - t25 * t250;
t372 = t249 * t33 + t250 * t34;
t371 = t249 * t36 - t250 * t35;
t57 = t250 * t429 + t231 + (t153 + t349) * qJD(1) + t440;
t370 = -t249 * t537 - t250 * t57;
t358 = t101 * t249 + t102 * t250;
t357 = t148 * t298 + t151 * t296;
t71 = t148 * t296 - t151 * t298;
t356 = t152 * t296 - t559;
t72 = -t152 * t298 - t560;
t355 = t153 * t249 + t154 * t250;
t116 = t296 * t532 + t298 * t354;
t238 = (Icges(7,5) * t295 + Icges(7,6) * t297) * t298;
t156 = qJD(5) * t200 + qJD(6) * t238;
t202 = -Icges(7,6) * t298 + t296 * t361;
t157 = (Icges(7,2) * t297 + t479) * t424 + t202 * qJD(5);
t204 = -Icges(7,5) * t298 + t296 * t363;
t158 = qJD(5) * t204 + qJD(6) * t240;
t40 = (-qJD(5) * t354 - t156) * t296 + (qJD(5) * t532 + t157 * t295 - t158 * t297 + (t205 * t295 - t297 * t531) * qJD(6)) * t298;
t348 = -t116 * t400 + t276 * t40;
t344 = -t298 * t41 + t427 * t92;
t343 = -t298 * t42 - t427 * t91;
t336 = -t143 + t434;
t333 = (t249 * t410 - t491) * qJD(5);
t54 = -t249 * t356 + t469;
t332 = (-t249 * t53 - t250 * t54) * qJD(5);
t331 = t388 + t470;
t329 = -t419 - t418;
t328 = -t415 - t419;
t324 = t336 + t554;
t322 = -t186 * t91 + t187 * t92 - t276 * t532;
t321 = (Icges(7,5) * t168 - Icges(7,6) * t169) * t187 + (Icges(7,5) * t170 + Icges(7,6) * t167) * t186 + t238 * t276;
t320 = t296 * t451 - t298 * t449;
t319 = -t296 * t450 + t298 * t448;
t317 = -rSges(5,2) * t250 - t418 + t494;
t316 = t337 + t436;
t313 = t237 + t288 + t329;
t254 = t378 * qJD(5);
t89 = -rSges(6,2) * t406 + t389;
t38 = (t233 * t379 + t249 * t254) * qJD(5) + (t89 + t508) * qJD(1) + t315;
t90 = rSges(6,1) * t339 + rSges(6,2) * t338 - t497;
t39 = (-t250 * t254 + t379 * t433) * qJD(5) + (-t90 + t452) * qJD(1) + t323;
t310 = t233 * t537 + t249 * t38 - t250 * t39 - t433 * t57;
t309 = t246 + t313;
t308 = -t153 * t233 + t154 * t433 - t249 * t90 + t250 * t89;
t304 = (-Icges(7,2) * t169 + t162 + t98) * t187 + (Icges(7,2) * t167 + t163 - t97) * t186 + (Icges(7,2) * t455 + t205 + t275) * t276;
t28 = -t101 * t186 + t102 * t187 - qJD(3) + (t181 * t250 + t249 * t543) * qJD(5);
t301 = t28 * t358 + t372 * t529;
t291 = t510 * rSges(3,3);
t283 = rSges(3,3) * t402;
t267 = t415 - t291;
t251 = t259 * qJD(5);
t195 = t285 + (-t266 - t267) * qJD(1);
t184 = t383 * t249;
t182 = t383 * t250;
t179 = t379 * t249;
t178 = t379 * t250;
t160 = -qJD(1) * t229 - t299 * t327 + t281;
t159 = qJD(1) * (-rSges(3,1) * t401 + t283) + t443;
t150 = t249 * t263 + t477;
t140 = qJD(1) * t540 + t285;
t139 = t529 * t249;
t138 = t529 * t250;
t135 = t530 * t249;
t134 = t530 * t250;
t133 = t531 * t249;
t132 = t531 * t250;
t118 = (-t229 + t380) * qJD(1) + t350;
t117 = qJD(1) * t442 + t335;
t115 = rSges(7,1) * t170 + rSges(7,2) * t167;
t114 = rSges(7,1) * t168 - rSges(7,2) * t169;
t112 = t250 * t352 + t173;
t103 = t112 * qJD(1);
t84 = Icges(6,5) * t339 + Icges(6,6) * t338 - Icges(6,3) * t433;
t83 = Icges(6,5) * t341 - Icges(6,6) * t340 + Icges(6,3) * t233;
t75 = (t317 + t542) * qJD(1) + t440;
t64 = qJD(5) * t355 - qJD(3);
t63 = (t375 + t452) * qJD(1) + t334;
t62 = qJD(1) * t441 + t315;
t50 = t233 * t352 + t249 * t528 - t250 * t251 + t360 * t433;
t49 = -t233 * t360 - t249 * t251 - t250 * t528 + t352 * t433;
t30 = qJD(5) * t356 + t296 * (Icges(6,4) * t339 + Icges(6,2) * t338 - Icges(6,6) * t433) - t298 * (Icges(6,1) * t339 + Icges(6,4) * t338 - Icges(6,5) * t433);
t29 = qJD(5) * t357 + t296 * (Icges(6,4) * t341 - Icges(6,2) * t340 + Icges(6,6) * t233) - t298 * (Icges(6,1) * t341 - Icges(6,4) * t340 + Icges(6,5) * t233);
t27 = t308 * qJD(5);
t22 = t332 - t556;
t21 = t103 + t333;
t20 = t156 * t462 + t157 * t170 - t158 * t167 + t205 * t82 + t338 * t532 - t531 * t81;
t19 = -t156 * t458 + t157 * t168 + t158 * t169 + t205 * t80 - t340 * t532 - t531 * t79;
t14 = t116 * t276 + t186 * t36 + t187 * t35;
t11 = t186 * t26 + t187 * t25 - t564;
t9 = -t101 * t136 + t102 * t137 - t186 * t47 + t187 * t48 + (t104 * t250 - t105 * t249 + t181 * t433 - t233 * t543) * qJD(5);
t6 = -t167 * t46 + t170 * t44 - t249 * t343 + t465 * t91 - t81 * t94 - t82 * t97;
t5 = -t167 * t45 + t170 * t43 - t249 * t344 - t465 * t92 + t81 * t95 + t82 * t98;
t4 = t168 * t44 + t169 * t46 + t250 * t343 + t467 * t91 - t79 * t94 - t80 * t97;
t3 = t168 * t43 + t169 * t45 + t250 * t344 - t467 * t92 + t79 * t95 + t80 * t98;
t2 = t136 * t26 + t137 * t25 + t186 * t6 + t187 * t5 + t20 * t276 + t400 * t545;
t1 = t136 * t24 + t137 * t23 + t186 * t4 + t187 * t3 + t19 * t276 - t400 * t60;
t15 = [(t63 * (t494 + (-rSges(5,2) + pkin(3)) * t250 + t313) + t62 * (t331 + t192) - (t316 - t391 + t441 + (-t317 + t329) * qJD(1)) * t536 + (t336 + t375 - t536 - t258) * t75) * m(5) + (t39 * t309 + t38 * (t388 + t154) + (-pkin(7) * t38 + t378 * t39) * t249 - (pkin(2) * t401 + t389 - t391 + t392 + t436 + t508 + (-t153 + t329) * qJD(1)) * t537 + (-rSges(6,1) * t466 - rSges(6,2) * t465 - t258 + t324 + t497 - t537 + t553) * t57 + (t39 * (rSges(6,3) + pkin(3)) + t38 * qJ(4) - (-rSges(6,2) * t427 + qJD(4) - t429 - t499) * t537) * t250) * m(6) + (t564 + (t24 + t569) * t187 + t11 + (-t23 + t568) * t186) * t519 + (t59 + (t26 + t568) * t187 + (-t25 - t569) * t186) * t521 + (t13 * (t309 + t377 + t507) + t12 * (-pkin(7) * t249 + t228 + t331 + t446) + (-pkin(5) * t495 + (pkin(5) * t13 + qJD(5) * t416) * t249) * t296 + (t233 * t416 + t511 * t527) * t298 + (-t190 + t257 + t316 + t409 + t420 + t508 - t562) * t34 + (t431 * t555 + t324 - t384 - t538) * t33 + ((t329 - t543 + t418 - t246) * t34 + (-t326 + t181) * t33) * qJD(1)) * m(7) + t505 / 0.2e1 + t506 / 0.2e1 + t19 * t518 + t489 / 0.2e1 + t488 / 0.2e1 - (t30 + t50 + t21) * t430 / 0.2e1 + t20 * t520 + t60 * t522 - t545 * t523 + (qJD(5) * t352 + t252 * t296 - t253 * t298) * qJD(1) + t348 + (t22 + t556 + (t469 * t250 + (-t52 + (t144 - t357) * t249 + (-t145 + (-t150 - t152) * t296) * t250) * t249) * qJD(5)) * t396 + (t160 * (t288 + t291 + t328) + t159 * t539 + (-qJD(1) * t539 + t286) * t195 + (t195 + t283 + (t267 + t328) * qJD(1) + t551) * (qJD(1) * t327 + t229)) * m(3) + (t118 * t540 + t117 * (t407 + t439) + (t380 - t408) * t140 + (t140 + t442 + (-t318 + t329) * qJD(1) + t551) * (qJD(1) * t439 + t408)) * m(4) + (t29 + t49 + qJD(1) * (t150 * t298 - t560 - t72)) * t397 + (t103 + (-t491 + (t469 - t54 + (t150 * t296 + t559) * t249 + t410) * t249) * qJD(5)) * t394 - t186 * t524 + ((t533 - t72) * t433 + (t112 + t71) * t233) * qJD(5) / 0.2e1; 0.2e1 * (t12 * t422 + t13 * t421) * m(7) + 0.2e1 * (t38 * t422 + t39 * t421) * m(6) + 0.2e1 * (t421 * t63 + t422 * t62) * m(5) + 0.2e1 * (t117 * t422 + t118 * t421) * m(4) + 0.2e1 * (t159 * t422 + t160 * t421) * m(3); -m(6) * t27 - m(7) * t9; m(5) * (-t233 * t536 - t249 * t62 + t250 * t63 + t433 * t75) - m(6) * t310 + (-t12 * t249 + t13 * t250 + t233 * t34 + t33 * t433) * m(7) + (-m(5) * (t249 * t75 - t250 * t536) - m(6) * (t249 * t57 - t250 * t537) - m(7) * t372) * qJD(1); qJD(1) * (t233 * t71 - t249 * t29 - t250 * t30 - t433 * t72) / 0.2e1 + (t233 * t35 - t249 * t7 - t250 * t8 - t36 * t433) * t512 + (((t132 * t295 - t134 * t297 - t92) * t187 + (-t133 * t295 + t135 * t297 + t91) * t186 + (t202 * t295 - t204 * t297 + t532) * t276 - t116 * qJD(6)) * t298 + (-qJD(6) * t371 - t525) * t296) * t513 + (t23 * t233 - t24 * t433 - t249 * t3 - t250 * t4) * t518 + ((t132 * t168 + t134 * t169) * t187 + (-t133 * t168 - t135 * t169) * t186 + (t168 * t202 + t169 * t204) * t276 + (-t24 * t463 - t298 * t60) * qJD(6) + ((qJD(6) * t23 + t322) * t296 - t558) * t250) * t519 + (t233 * t25 - t249 * t5 - t250 * t6 - t26 * t433) * t520 + ((t132 * t170 - t134 * t167) * t187 + (-t133 * t170 + t135 * t167) * t186 + (-t167 * t204 + t170 * t202) * t276 + (t25 * t459 + t298 * t545) * qJD(6) + ((-qJD(6) * t26 - t322) * t296 + t558) * t249) * t521 + (-t23 * t249 - t24 * t250) * t522 + (-t249 * t25 - t250 * t26) * t523 + t403 * t524 + (-t249 * t35 - t250 * t36) * t387 + ((t172 * t431 - t432) * t249 + (-t550 + (-t319 * t250 + (-t173 + t320) * t249) * qJD(5)) * t250) * t396 + ((-t173 * t430 - t432) * t250 + (t550 + (-t320 * t249 + (t172 + t319) * t250) * qJD(5)) * t249) * t394 + t11 * t404 / 0.2e1 + t14 * t424 / 0.2e1 - qJD(1) * ((t438 * t296 + t437 * t298) * qJD(1) + ((-t249 * t451 - t250 * t450) * t298 + (-t249 * t449 - t250 * t448) * t296) * qJD(5)) / 0.2e1 + (-t33 * (qJD(1) * t184 + t139 * t276 + t186 * t206 - t461) - t34 * (qJD(1) * t182 + t138 * t276 - t187 * t206 + t464) - t28 * (-t138 * t186 - t139 * t187 + t182 * t430 + t184 * t431) - ((-t101 * t34 + t102 * t33) * t298 + t301 * t296) * qJD(6) + (t28 * t453 - t34 * t444) * t233 - (-t28 * t454 + t33 * t444) * t433 + (-t13 * t444 - t33 * t447 + t9 * t454 + t28 * (t104 + t47)) * t250 + (t12 * t444 + t34 * t447 - t9 * t453 + t28 * (-t105 - t48)) * t249) * m(7) + (-(-t178 * t537 + t179 * t57) * qJD(1) - (t64 * (t178 * t250 + t179 * t249) + t370 * t378) * qJD(5) + t370 * t254 - t310 * t379 + t27 * t355 + t64 * t308) * m(6) - (qJD(1) * t49 + t1 + (-(t145 * t233 - t249 * t83 + t357 * t433) * t249 - (-t144 * t233 - t249 * t84 + t356 * t433) * t250 - t433 * t52 - t233 * t410) * t544) * t249 / 0.2e1 - (qJD(1) * t50 + t2 + (-(-t145 * t433 + t233 * t357 - t250 * t83) * t249 - (t144 * t433 + t233 * t356 - t250 * t84) * t250 - t433 * t54 + t233 * t53) * t544) * t250 / 0.2e1 - (t332 + t22 + t11) * t433 / 0.2e1 + (t333 + t21 + t10) * t233 / 0.2e1; -t1 * t458 / 0.2e1 + (-t296 * t60 - t298 * t374) * t522 + ((qJD(5) * t374 - t19) * t296 + (-qJD(5) * t60 - t23 * t433 - t233 * t24 + t249 * t4 - t250 * t3) * t298) * t518 + t2 * t462 / 0.2e1 + (t296 * t545 + t298 * t373) * t523 + ((-qJD(5) * t373 - t20) * t296 + (qJD(5) * t545 - t233 * t26 + t249 * t6 - t25 * t433 - t250 * t5) * t298) * t520 + t14 * t393 - t296 * (t348 + t488 + t489 + t505 + t506) / 0.2e1 + (-t116 * t296 + t298 * t371) * t387 + ((-qJD(5) * t371 - t40) * t296 + (-qJD(5) * t116 - t233 * t36 + t249 * t8 - t250 * t7 - t35 * t433) * t298) * t512 + (t168 * t304 + t169 * t303 - t321 * t458) * t519 + (-t167 * t303 + t170 * t304 + t321 * t462) * t521 + (-t321 * t296 + (t304 * t295 - t297 * t303) * t298) * t513 + (-t465 / 0.2e1 + t296 * t397) * t11 + (-t467 / 0.2e1 + t296 * t394) * t10 + ((qJD(5) * t301 - t101 * t12 + t102 * t13 + t33 * t48 - t34 * t47) * t296 + (t33 * (qJD(5) * t102 + t161 * t249) + t34 * (-qJD(5) * t101 + t161 * t250) - t9 * t358 + t28 * (t101 * t233 - t102 * t433 - t249 * t47 - t250 * t48) - (-t495 - t527) * t529) * t298 - t33 * (-t115 * t276 + t186 * t248) - t34 * (t114 * t276 - t187 * t248) - t28 * (-t114 * t186 + t115 * t187)) * m(7);];
tauc  = t15(:);
