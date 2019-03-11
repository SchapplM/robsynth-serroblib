% Calculate time derivative of joint inertia matrix for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:56
% EndTime: 2019-03-09 02:23:14
% DurationCPUTime: 11.31s
% Computational Cost: add. (36634->798), mult. (38218->1118), div. (0->0), fcn. (36756->10), ass. (0->404)
t308 = cos(qJ(5));
t305 = sin(qJ(5));
t306 = sin(qJ(4));
t378 = qJD(5) * t306 + qJD(1);
t331 = t378 * t305;
t309 = cos(qJ(4));
t429 = qJD(4) * t309;
t527 = -t308 * t429 + t331;
t303 = qJ(1) + pkin(10);
t297 = sin(t303);
t298 = cos(t303);
t430 = qJD(4) * t306;
t390 = -t430 / 0.2e1;
t433 = qJD(1) * t309;
t391 = -t433 / 0.2e1;
t526 = t297 * t391 + t298 * t390;
t389 = t430 / 0.2e1;
t525 = t297 * t389 + t298 * t391;
t401 = t298 * t430;
t406 = t297 * t433;
t524 = t401 + t406;
t434 = qJD(1) * t306;
t523 = t297 * t429 + t298 * t434;
t484 = Icges(6,4) * t308;
t349 = -Icges(6,2) * t305 + t484;
t245 = Icges(6,6) * t306 + t309 * t349;
t485 = Icges(6,4) * t305;
t353 = Icges(6,1) * t308 - t485;
t246 = Icges(6,5) * t306 + t309 * t353;
t522 = -t245 * t305 + t246 * t308;
t304 = qJ(5) + qJ(6);
t299 = sin(t304);
t300 = cos(t304);
t482 = Icges(7,4) * t300;
t348 = -Icges(7,2) * t299 + t482;
t235 = Icges(7,6) * t306 + t309 * t348;
t483 = Icges(7,4) * t299;
t352 = Icges(7,1) * t300 - t483;
t236 = Icges(7,5) * t306 + t309 * t352;
t521 = -t235 * t299 + t236 * t300;
t509 = 2 * m(5);
t291 = t298 * rSges(5,3);
t459 = t297 * t309;
t460 = t297 * t306;
t223 = rSges(5,1) * t460 + rSges(5,2) * t459 + t291;
t368 = rSges(5,1) * t306 + rSges(5,2) * t309;
t328 = t298 * t368;
t491 = t297 * rSges(5,3);
t224 = t491 - t328;
t494 = rSges(5,2) * t306;
t279 = rSges(5,1) * t309 - t494;
t295 = t297 ^ 2;
t296 = t298 ^ 2;
t404 = t298 * t433;
t412 = rSges(5,1) * t523 + rSges(5,2) * t404;
t90 = -t297 * t412 + (-t279 * t296 + t295 * t494) * qJD(4) + ((-t223 + t291) * t298 + (-t224 + t328 + t491) * t297) * qJD(1);
t520 = t509 * t90;
t294 = pkin(5) * t308 + pkin(4);
t496 = pkin(4) - t294;
t519 = t306 * t496;
t486 = Icges(5,4) * t309;
t354 = Icges(5,1) * t306 + t486;
t221 = Icges(5,5) * t298 + t297 * t354;
t471 = t221 * t306;
t487 = Icges(5,4) * t306;
t350 = Icges(5,2) * t309 + t487;
t219 = Icges(5,6) * t298 + t297 * t350;
t474 = t219 * t309;
t336 = t471 + t474;
t518 = t336 * t297;
t311 = -pkin(9) - pkin(8);
t495 = -pkin(8) - t311;
t229 = t306 * t495 - t309 * t496;
t362 = rSges(7,1) * t300 - rSges(7,2) * t299;
t237 = rSges(7,3) * t306 + t309 * t362;
t441 = t229 + t237;
t386 = t441 * t297;
t506 = -pkin(2) - pkin(7);
t426 = -rSges(5,3) + t506;
t497 = sin(qJ(1)) * pkin(1);
t326 = t297 * t426 - t497;
t403 = t297 * t430;
t435 = qJD(1) * t298;
t438 = qJ(3) * t435 + qJD(3) * t297;
t137 = -rSges(5,2) * t403 + qJD(1) * t326 + t412 + t438;
t288 = qJD(3) * t298;
t301 = cos(qJ(1)) * pkin(1);
t431 = qJD(4) * t298;
t138 = t288 + t279 * t431 + (-t301 + t426 * t298 + (-qJ(3) - t368) * t297) * qJD(1);
t517 = -t137 * t298 + t138 * t297;
t346 = Icges(5,5) * t306 + Icges(5,6) * t309;
t516 = -Icges(5,3) * t297 + t298 * t346;
t515 = -Icges(5,6) * t297 + t298 * t350;
t514 = -Icges(5,5) * t297 + t298 * t354;
t302 = qJD(5) + qJD(6);
t384 = t302 * t306 + qJD(1);
t513 = t299 * t384 - t300 * t429;
t512 = t299 * t429 + t300 * t384;
t511 = t305 * t429 + t308 * t378;
t462 = t294 * t306;
t490 = rSges(7,3) - t311;
t510 = -t309 * t490 + t462;
t508 = 2 * m(6);
t507 = 2 * m(7);
t505 = -t297 / 0.2e1;
t504 = t298 / 0.2e1;
t503 = -t306 / 0.2e1;
t502 = t306 / 0.2e1;
t501 = t309 / 0.2e1;
t500 = rSges(4,2) - pkin(2);
t499 = m(5) * t279;
t498 = pkin(4) * t306;
t493 = rSges(6,3) * t309;
t492 = pkin(5) * qJD(5);
t344 = Icges(7,5) * t300 - Icges(7,6) * t299;
t234 = Icges(7,3) * t306 + t309 * t344;
t136 = t234 * t306 + t309 * t521;
t453 = t302 * t309;
t183 = (-Icges(7,2) * t300 - t483) * t453 + (Icges(7,6) * t309 - t306 * t348) * qJD(4);
t182 = (-Icges(7,5) * t299 - Icges(7,6) * t300) * t453 + (Icges(7,3) * t309 - t306 * t344) * qJD(4);
t184 = (-Icges(7,1) * t299 - t482) * t453 + (Icges(7,5) * t309 - t306 * t352) * qJD(4);
t317 = t309 * t300 * t184 + t306 * t182 + t234 * t429 - t430 * t521;
t420 = t302 * t300 * t235;
t489 = t136 * t429 + ((-t420 + (-t236 * t302 - t183) * t299) * t309 + t317) * t306;
t345 = Icges(6,5) * t308 - Icges(6,6) * t305;
t244 = Icges(6,3) * t306 + t309 * t345;
t142 = t244 * t306 + t309 * t522;
t427 = qJD(5) * t309;
t202 = (-Icges(6,5) * t305 - Icges(6,6) * t308) * t427 + (Icges(6,3) * t309 - t306 * t345) * qJD(4);
t204 = (-Icges(6,1) * t305 - t484) * t427 + (Icges(6,5) * t309 - t306 * t353) * qJD(4);
t316 = t309 * t308 * t204 + t306 * t202 + t244 * t429 - t430 * t522;
t203 = (-Icges(6,2) * t308 - t485) * t427 + (Icges(6,6) * t309 - t306 * t349) * qJD(4);
t452 = t305 * t203;
t464 = t245 * t308;
t488 = t142 * t429 + ((-t452 + (-t246 * t305 - t464) * qJD(5)) * t309 + t316) * t306;
t451 = t305 * t306;
t242 = t297 * t308 + t298 * t451;
t450 = t306 * t308;
t461 = t297 * t305;
t243 = -t298 * t450 + t461;
t366 = -rSges(6,1) * t243 - rSges(6,2) * t242;
t456 = t298 * t309;
t181 = rSges(6,3) * t456 - t366;
t476 = t181 * t298;
t475 = t219 * t306;
t473 = t515 * t306;
t472 = t515 * t309;
t470 = t221 * t309;
t469 = t514 * t306;
t468 = t514 * t309;
t458 = t298 * t305;
t457 = t298 * t306;
t455 = t299 * t306;
t454 = t300 * t306;
t449 = t309 * t294;
t448 = t309 * t311;
t383 = t302 + t434;
t334 = t299 * t383;
t147 = -t297 * t334 + t298 * t512;
t333 = t300 * t383;
t148 = t297 * t333 + t298 * t513;
t364 = t148 * rSges(7,1) + t147 * rSges(7,2);
t101 = -rSges(7,3) * t524 + t364;
t284 = pkin(5) * t458;
t424 = t308 * t492;
t380 = t297 * t424;
t425 = t305 * t492;
t381 = t306 * t425;
t410 = t298 * pkin(4) * t429 + pkin(8) * t524;
t447 = t101 + t380 + (t381 + (t306 * t311 - t449) * qJD(4)) * t298 + (t284 + (t448 - t519) * t297) * qJD(1) + t410;
t149 = -t297 * t512 - t298 * t334;
t150 = -t297 * t513 + t298 * t333;
t417 = t150 * rSges(7,1) + t149 * rSges(7,2) + rSges(7,3) * t403;
t102 = -rSges(7,3) * t404 + t417;
t411 = pkin(4) * t523 + pkin(8) * t403;
t320 = pkin(8) * t404 - t411;
t376 = t294 * t523 + t298 * t424 + t311 * t404;
t428 = qJD(4) * t311;
t117 = (-pkin(5) * t331 - t306 * t428) * t297 + t320 + t376;
t446 = -t102 - t117;
t230 = -t297 * t455 + t298 * t300;
t231 = t297 * t454 + t298 * t299;
t162 = t231 * rSges(7,1) + t230 * rSges(7,2) - rSges(7,3) * t459;
t283 = pkin(4) * t460;
t248 = -pkin(8) * t459 + t283;
t413 = t294 * t460 + t297 * t448 + t284;
t185 = -t248 + t413;
t445 = -t162 - t185;
t232 = t297 * t300 + t298 * t455;
t233 = t297 * t299 - t298 * t454;
t363 = -t233 * rSges(7,1) - t232 * rSges(7,2);
t423 = rSges(7,3) * t456;
t163 = -t363 + t423;
t285 = pkin(4) * t457;
t392 = t495 * t309;
t186 = pkin(5) * t461 + t285 + (t392 - t462) * t298;
t444 = t163 + t186;
t240 = -t297 * t451 + t298 * t308;
t241 = t297 * t450 + t458;
t440 = t241 * rSges(6,1) + t240 * rSges(6,2);
t180 = -rSges(6,3) * t459 + t440;
t443 = -t180 - t248;
t193 = (-rSges(7,1) * t299 - rSges(7,2) * t300) * t453 + (rSges(7,3) * t309 - t306 * t362) * qJD(4);
t396 = t305 * t427;
t208 = -pkin(5) * t396 + (t392 + t519) * qJD(4);
t442 = t193 + t208;
t131 = t306 * t162 + t237 * t459;
t262 = (pkin(8) * t309 - t498) * qJD(4);
t282 = t309 * pkin(4) + t306 * pkin(8);
t439 = t297 * t262 + t282 * t435;
t217 = Icges(5,3) * t298 + t297 * t346;
t437 = qJD(1) * t217;
t436 = qJD(1) * t297;
t432 = qJD(4) * t297;
t377 = qJD(5) + t434;
t332 = t305 * t377;
t171 = -t297 * t511 - t298 * t332;
t330 = t377 * t308;
t172 = -t297 * t527 + t298 * t330;
t322 = t403 - t404;
t104 = Icges(6,5) * t172 + Icges(6,6) * t171 + Icges(6,3) * t322;
t106 = Icges(6,4) * t172 + Icges(6,2) * t171 + Icges(6,6) * t322;
t108 = Icges(6,1) * t172 + Icges(6,4) * t171 + Icges(6,5) * t322;
t174 = Icges(6,5) * t241 + Icges(6,6) * t240 - Icges(6,3) * t459;
t176 = Icges(6,4) * t241 + Icges(6,2) * t240 - Icges(6,6) * t459;
t178 = Icges(6,1) * t241 + Icges(6,4) * t240 - Icges(6,5) * t459;
t341 = t176 * t305 - t178 * t308;
t30 = (qJD(4) * t341 + t104) * t306 + (qJD(4) * t174 - t106 * t305 + t108 * t308 + (-t176 * t308 - t178 * t305) * qJD(5)) * t309;
t59 = t171 * t245 + t172 * t246 - t202 * t459 + t203 * t240 + t204 * t241 + t244 * t322;
t422 = -t30 / 0.2e1 - t59 / 0.2e1;
t169 = -t297 * t332 + t298 * t511;
t170 = t297 * t330 + t298 * t527;
t103 = Icges(6,5) * t170 + Icges(6,6) * t169 - Icges(6,3) * t524;
t105 = Icges(6,4) * t170 + Icges(6,2) * t169 - Icges(6,6) * t524;
t107 = Icges(6,1) * t170 + Icges(6,4) * t169 - Icges(6,5) * t524;
t175 = Icges(6,5) * t243 + Icges(6,6) * t242 + Icges(6,3) * t456;
t177 = Icges(6,4) * t243 + Icges(6,2) * t242 + Icges(6,6) * t456;
t179 = Icges(6,1) * t243 + Icges(6,4) * t242 + Icges(6,5) * t456;
t340 = t177 * t305 - t179 * t308;
t31 = (qJD(4) * t340 + t103) * t306 + (qJD(4) * t175 - t105 * t305 + t107 * t308 + (-t177 * t308 - t179 * t305) * qJD(5)) * t309;
t58 = t169 * t245 + t170 * t246 + t202 * t456 + t203 * t242 + t204 * t243 - t244 * t524;
t421 = t31 / 0.2e1 + t58 / 0.2e1;
t123 = t240 * t245 + t241 * t246 - t244 * t459;
t91 = t174 * t306 - t309 * t341;
t419 = t91 / 0.2e1 + t123 / 0.2e1;
t124 = t242 * t245 + t243 * t246 + t244 * t456;
t92 = t175 * t306 - t309 * t340;
t418 = -t92 / 0.2e1 - t124 / 0.2e1;
t416 = t162 * t524 + t163 * t403;
t415 = -t248 + t445;
t414 = t172 * rSges(6,1) + t171 * rSges(6,2) + rSges(6,3) * t403;
t409 = t298 * pkin(2) + t297 * qJ(3) + t301;
t408 = (-rSges(6,3) - pkin(8)) * t309;
t365 = rSges(6,1) * t308 - rSges(6,2) * t305;
t247 = rSges(6,3) * t306 + t309 * t365;
t407 = t247 * t436;
t395 = -t459 / 0.2e1;
t394 = t456 / 0.2e1;
t156 = Icges(7,5) * t231 + Icges(7,6) * t230 - Icges(7,3) * t459;
t158 = Icges(7,4) * t231 + Icges(7,2) * t230 - Icges(7,6) * t459;
t160 = Icges(7,1) * t231 + Icges(7,4) * t230 - Icges(7,5) * t459;
t343 = t158 * t299 - t160 * t300;
t95 = Icges(7,5) * t150 + Icges(7,6) * t149 + Icges(7,3) * t322;
t97 = Icges(7,4) * t150 + Icges(7,2) * t149 + Icges(7,6) * t322;
t99 = Icges(7,1) * t150 + Icges(7,4) * t149 + Icges(7,5) * t322;
t24 = (qJD(4) * t343 + t95) * t306 + (qJD(4) * t156 + (-t158 * t302 + t99) * t300 + (-t160 * t302 - t97) * t299) * t309;
t157 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t456;
t159 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t456;
t161 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t456;
t342 = t159 * t299 - t161 * t300;
t94 = Icges(7,5) * t148 + Icges(7,6) * t147 - Icges(7,3) * t524;
t96 = Icges(7,4) * t148 + Icges(7,2) * t147 - Icges(7,6) * t524;
t98 = Icges(7,1) * t148 + Icges(7,4) * t147 - Icges(7,5) * t524;
t25 = (qJD(4) * t342 + t94) * t306 + (qJD(4) * t157 + (-t159 * t302 + t98) * t300 + (-t161 * t302 - t96) * t299) * t309;
t84 = t156 * t306 - t309 * t343;
t85 = t157 * t306 - t309 * t342;
t356 = t84 * t297 - t85 * t298;
t357 = t297 * t85 + t298 * t84;
t118 = t230 * t235 + t231 * t236 - t234 * t459;
t70 = -t156 * t459 + t158 * t230 + t160 * t231;
t71 = -t157 * t459 + t159 * t230 + t161 * t231;
t361 = t297 * t70 - t298 * t71;
t38 = t118 * t306 - t309 * t361;
t119 = t232 * t235 + t233 * t236 + t234 * t456;
t18 = t147 * t158 + t148 * t160 - t156 * t524 + t232 * t97 + t233 * t99 + t456 * t95;
t19 = t147 * t159 + t148 * t161 - t157 * t524 + t232 * t96 + t233 * t98 + t456 * t94;
t72 = t156 * t456 + t158 * t232 + t160 * t233;
t73 = t157 * t456 + t159 * t232 + t161 * t233;
t360 = t297 * t72 - t298 * t73;
t50 = t297 * t73 + t298 * t72;
t53 = t147 * t235 + t148 * t236 + t182 * t456 + t183 * t232 + t184 * t233 - t234 * t524;
t4 = (qJD(4) * t360 + t53) * t306 + (-qJD(1) * t50 + qJD(4) * t119 - t18 * t297 + t19 * t298) * t309;
t393 = t4 * t456 + t38 * t403 + (t136 * t306 - t309 * t356) * t429 + t306 * (t356 * t430 + (-qJD(1) * t357 - t24 * t297 + t25 * t298) * t309 + t489);
t80 = -t174 * t459 + t176 * t240 + t178 * t241;
t81 = -t175 * t459 + t177 * t240 + t179 * t241;
t359 = t297 * t80 - t298 * t81;
t41 = t306 * t123 - t309 * t359;
t388 = t306 * t91 + t41;
t260 = t368 * qJD(4);
t387 = t260 * (t295 + t296);
t385 = qJD(1) * t441;
t382 = -pkin(5) * t305 + t506;
t379 = t306 * t102 + t162 * t429 + t193 * t459 + t237 * t404;
t375 = t298 * pkin(7) + t409;
t39 = t119 * t306 - t309 * t360;
t82 = t174 * t456 + t176 * t242 + t178 * t243;
t83 = t175 * t456 + t177 * t242 + t179 * t243;
t358 = t297 * t82 - t298 * t83;
t42 = t306 * t124 - t309 * t358;
t374 = -t306 * t92 - t39 - t42;
t367 = rSges(6,1) * t170 + rSges(6,2) * t169;
t49 = t297 * t71 + t298 * t70;
t55 = t297 * t81 + t298 * t80;
t56 = t297 * t83 + t298 * t82;
t355 = Icges(5,1) * t309 - t487;
t351 = -Icges(5,2) * t306 + t486;
t347 = Icges(5,5) * t309 - Icges(5,6) * t306;
t339 = t180 * t298 + t181 * t297;
t335 = -t469 - t472;
t205 = (-rSges(6,1) * t305 - rSges(6,2) * t308) * t427 + (-t306 * t365 + t493) * qJD(4);
t329 = t205 * t297 + t247 * t435;
t327 = t335 * t297;
t325 = qJD(4) * t355;
t324 = qJD(4) * t351;
t319 = t297 * t382 - t497;
t318 = rSges(4,3) * t298 + t297 * t500 - t497;
t11 = -qJD(1) * t360 + t18 * t298 + t19 * t297;
t20 = t149 * t158 + t150 * t160 + t156 * t322 + t230 * t97 + t231 * t99 - t459 * t95;
t21 = t149 * t159 + t150 * t161 + t157 * t322 + t230 * t96 + t231 * t98 - t459 * t94;
t12 = -qJD(1) * t361 + t20 * t298 + t21 * t297;
t54 = t149 * t235 + t150 * t236 - t182 * t459 + t183 * t230 + t184 * t231 + t234 * t322;
t5 = (qJD(4) * t361 + t54) * t306 + (-qJD(1) * t49 + qJD(4) * t118 - t20 * t297 + t21 * t298) * t309;
t315 = t11 * t394 + (-qJD(1) * t356 + t24 * t298 + t25 * t297) * t502 + t297 * t4 / 0.2e1 + t5 * t504 - t38 * t436 / 0.2e1 + t39 * t435 / 0.2e1 + t357 * t429 / 0.2e1 + t12 * t395 + t526 * t50 + t525 * t49;
t314 = t297 * t506 + t298 * t408 - t497;
t313 = t489 + (t24 + t54) * t395 + (t25 + t53) * t394 + (t119 + t85) * t526 + (t118 + t84) * t525;
t312 = (-t297 * t5 + (-t297 * t39 - t298 * t38) * qJD(1)) * t309 - t39 * t401 + t393;
t290 = t298 * qJ(3);
t255 = t297 * t282;
t251 = t282 * t436;
t249 = pkin(8) * t456 - t285;
t227 = t298 * t249;
t214 = t237 * t456;
t212 = -rSges(4,2) * t298 + rSges(4,3) * t297 + t409;
t211 = t290 + t318;
t207 = (-t247 - t282) * t298;
t206 = t247 * t297 + t255;
t199 = t288 + (-t301 + t500 * t298 + (-rSges(4,3) - qJ(3)) * t297) * qJD(1);
t198 = qJD(1) * t318 + t438;
t196 = t298 * (qJD(1) * t283 - t410);
t195 = t375 + t223;
t194 = t290 + t328 + t326;
t188 = qJD(1) * t516 + t347 * t432;
t187 = -t347 * t431 + t437;
t166 = t193 * t456;
t144 = (-t282 - t441) * t298;
t143 = t255 + t386;
t140 = -t181 * t306 + t247 * t456;
t139 = t180 * t306 + t247 * t459;
t134 = t297 * t408 + t283 + t375 + t440;
t133 = t285 + t290 + t314 + t366;
t132 = -t163 * t306 + t214;
t130 = -t297 * t516 - t298 * t335;
t129 = t217 * t297 - t298 * t336;
t128 = -t298 * t516 + t327;
t127 = t217 * t298 + t518;
t126 = t329 + t439;
t125 = t407 + t251 + (-t205 - t262) * t298;
t122 = t162 + t375 + t413;
t121 = t298 * t510 + t290 + t319 + t363;
t120 = t339 * t309;
t115 = (-t162 * t298 - t163 * t297) * t309;
t110 = -rSges(6,3) * t404 + t414;
t109 = -rSges(6,3) * t524 + t367;
t100 = t297 * t443 + t227 + t476;
t89 = t229 * t456 - t306 * t444 + t214;
t88 = t185 * t306 + t229 * t459 + t131;
t87 = t297 * t442 + t298 * t385 + t439;
t86 = t251 + t297 * t385 + (-t262 - t442) * t298;
t75 = rSges(6,3) * t401 + t288 + (-t301 + t506 * t298 + (-qJ(3) + t493 - t498) * t297) * qJD(1) - t367 + t410;
t74 = qJD(1) * t314 + t411 + t414 + t438;
t69 = (-t297 * t444 + t298 * t445) * t309;
t67 = -t380 + t288 + (-t381 + (t306 * t490 + t449) * qJD(4)) * t298 + (-t301 + t382 * t298 + (-qJ(3) - t510) * t297) * qJD(1) - t364;
t66 = (-t425 - t428) * t460 + (t319 - t423) * qJD(1) + t376 + t417 + t438;
t65 = t297 * t415 + t298 * t444 + t227;
t64 = (-t247 * t432 + t110) * t306 + (qJD(4) * t180 + t329) * t309;
t63 = (-t247 * t431 - t109) * t306 + (-qJD(4) * t181 + t205 * t298 - t407) * t309;
t61 = -t237 * t403 + t379;
t60 = -t237 * t406 - t101 * t306 + t166 + (-t163 * t309 - t237 * t457) * qJD(4);
t43 = t339 * t430 + (-t109 * t297 - t110 * t298 + (t180 * t297 - t476) * qJD(1)) * t309;
t40 = t109 * t298 + t196 + (-t110 + t320) * t297 + (t443 * t298 + (-t181 - t249) * t297) * qJD(1);
t37 = (-t101 * t297 + (-qJD(1) * t163 - t102) * t298) * t309 + t416;
t33 = t117 * t306 + (t208 * t297 + t229 * t435) * t309 + (t185 * t309 - t306 * t386) * qJD(4) + t379;
t32 = t166 + (-t431 * t441 - t447) * t306 + (-qJD(1) * t386 - qJD(4) * t444 + t208 * t298) * t309;
t29 = -t103 * t459 + t105 * t240 + t107 * t241 + t171 * t177 + t172 * t179 + t175 * t322;
t28 = -t104 * t459 + t106 * t240 + t108 * t241 + t171 * t176 + t172 * t178 + t174 * t322;
t27 = t103 * t456 + t105 * t242 + t107 * t243 + t169 * t177 + t170 * t179 - t175 * t524;
t26 = t104 * t456 + t106 * t242 + t108 * t243 + t169 * t176 + t170 * t178 - t174 * t524;
t17 = t196 + t447 * t298 + (t320 + t446) * t297 + (t415 * t298 + (-t249 - t444) * t297) * qJD(1);
t16 = (t185 * t298 + t186 * t297) * t430 + ((qJD(1) * t185 - t447) * t297 + (-qJD(1) * t444 + t446) * t298) * t309 + t416;
t15 = -qJD(1) * t359 + t28 * t298 + t29 * t297;
t14 = -qJD(1) * t358 + t26 * t298 + t27 * t297;
t8 = (qJD(4) * t359 + t59) * t306 + (-qJD(1) * t55 + qJD(4) * t123 - t28 * t297 + t29 * t298) * t309;
t7 = (qJD(4) * t358 + t58) * t306 + (-qJD(1) * t56 + qJD(4) * t124 - t26 * t297 + t27 * t298) * t309;
t1 = [0.2e1 * m(4) * (t198 * t212 + t199 * t211) + (t137 * t195 + t138 * t194) * t509 + (t133 * t75 + t134 * t74) * t508 + (t121 * t67 + t122 * t66) * t507 - t354 * t429 + t350 * t430 - t246 * t396 - t427 * t464 - t299 * t236 * t453 + t316 + t317 - t306 * t325 + (-t183 * t299 - t324 - t420 - t452) * t309; 0; 0; m(7) * (t297 * t67 - t298 * t66 + (t121 * t298 + t122 * t297) * qJD(1)) + m(6) * (t297 * t75 - t298 * t74 + (t133 * t298 + t134 * t297) * qJD(1)) + m(5) * ((t194 * t298 + t195 * t297) * qJD(1) + t517) + m(4) * (-t198 * t298 + t199 * t297 + (t211 * t298 + t212 * t297) * qJD(1)); 0; 0; ((qJD(1) * t515 + t297 * t324) * t503 + (qJD(1) * t514 + t297 * t325) * t501 + t54 / 0.2e1 + t24 / 0.2e1 + (-t474 / 0.2e1 - t471 / 0.2e1) * qJD(4) - t422) * t298 + ((qJD(1) * t219 - t351 * t431) * t503 + (qJD(1) * t221 - t355 * t431) * t501 + t53 / 0.2e1 + t25 / 0.2e1 + (t472 / 0.2e1 + t469 / 0.2e1) * qJD(4) + t421) * t297 + m(5) * (t517 * t279 - (t194 * t297 - t195 * t298) * t260) + m(6) * (t125 * t134 + t126 * t133 + t206 * t75 + t207 * t74) + m(7) * (t121 * t87 + t122 * t86 + t143 * t67 + t144 * t66) - (t295 / 0.2e1 + t296 / 0.2e1) * t346 * qJD(4) + ((t475 / 0.2e1 - t470 / 0.2e1 + t195 * t499 - t118 / 0.2e1 - t84 / 0.2e1 - t419) * t297 + (t473 / 0.2e1 - t468 / 0.2e1 + t194 * t499 + t119 / 0.2e1 + t85 / 0.2e1 - t418) * t298) * qJD(1); m(5) * t90 + m(6) * t40 + m(7) * t17; m(6) * (-t125 * t298 + t126 * t297 + (t206 * t298 + t207 * t297) * qJD(1)) + m(7) * (t297 * t87 - t298 * t86 + (t143 * t298 + t144 * t297) * qJD(1)) - m(5) * t387; (t143 * t87 + t144 * t86 + t65 * t17) * t507 + (t100 * t40 + t125 * t207 + t126 * t206) * t508 - t279 * t387 * t509 + (t224 * t520 + t12 + t15 + (t128 * qJD(1) + (t336 * qJD(1) + t188) * t298) * t298) * t298 + (t11 + t14 - t223 * t520 + (t297 * t187 + (-t129 + t327) * qJD(1)) * t297 + ((t515 * t430 - t514 * t429 + t188 + (t468 - t473) * qJD(4)) * t297 + (t187 + (t470 - t475) * qJD(4) + t219 * t430 - t221 * t429 + t437) * t298 + (-t127 + t130 + t518 + (-t217 + t335) * t298) * qJD(1)) * t298) * t297 + (-t127 * t298 - t128 * t297 - t49 - t55) * t436 + (t129 * t298 + t130 * t297 + t50 + t56) * t435; (t421 * t298 + t422 * t297 + (t297 * t418 - t298 * t419) * qJD(1)) * t309 + m(6) * (t133 * t63 + t134 * t64 + t139 * t74 + t140 * t75) + m(7) * (t121 * t32 + t122 * t33 + t66 * t88 + t67 * t89) + t313 + (t297 * t419 + t298 * t418) * t430 + t488; m(6) * t43 + m(7) * t16; m(6) * (t297 * t63 - t298 * t64 + (t139 * t297 + t140 * t298) * qJD(1)) + m(7) * (t297 * t32 - t298 * t33 + (t297 * t88 + t298 * t89) * qJD(1)); (qJD(4) * (t92 * t297 + t91 * t298) / 0.2e1 + t14 * t504 + t15 * t505 + (-t298 * t55 / 0.2e1 + t56 * t505) * qJD(1)) * t309 + (t56 * t390 + (qJD(1) * t92 + t30) * t502 + t8 / 0.2e1 + qJD(1) * t42 / 0.2e1) * t298 + (t55 * t389 + (-qJD(1) * t91 + t31) * t502 + t7 / 0.2e1 - qJD(1) * t41 / 0.2e1) * t297 + t315 + m(6) * (t100 * t43 - t120 * t40 + t125 * t139 + t126 * t140 + t206 * t63 + t207 * t64) + m(7) * (t143 * t32 + t144 * t33 + t16 * t65 + t69 * t17 + t86 * t88 + t87 * t89); (t69 * t16 + t32 * t89 + t33 * t88) * t507 + (-t120 * t43 + t139 * t64 + t140 * t63) * t508 + ((t297 * t388 + t298 * t374) * qJD(4) + t488) * t306 + ((t306 * t31 + t7) * t298 + (-t30 * t306 - t5 - t8) * t297 + (t142 * t306 + (-t91 * t297 + t92 * t298) * t309) * qJD(4) + ((-t38 - t388) * t298 + t374 * t297) * qJD(1)) * t309 + t393; m(7) * (t121 * t60 + t122 * t61 + t131 * t66 + t132 * t67) + t313; m(7) * t37; m(7) * (t297 * t60 - t298 * t61 + (t131 * t297 + t132 * t298) * qJD(1)); m(7) * (t115 * t17 + t131 * t86 + t132 * t87 + t143 * t60 + t144 * t61 + t37 * t65) + t315; m(7) * (t115 * t16 + t131 * t33 + t132 * t32 + t37 * t69 + t60 * t89 + t61 * t88) + t312; (t115 * t37 + t131 * t61 + t132 * t60) * t507 + t312;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
