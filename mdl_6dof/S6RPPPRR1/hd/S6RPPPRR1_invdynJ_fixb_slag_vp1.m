% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:51
% EndTime: 2019-03-09 01:30:13
% DurationCPUTime: 19.34s
% Computational Cost: add. (15789->849), mult. (20341->1101), div. (0->0), fcn. (18404->8), ass. (0->398)
t306 = cos(qJ(1));
t299 = t306 * pkin(1);
t300 = qJ(1) + pkin(9);
t295 = sin(t300);
t296 = cos(t300);
t215 = rSges(3,1) * t295 + rSges(3,2) * t296;
t303 = sin(qJ(1));
t516 = pkin(1) * t303;
t198 = -t215 - t516;
t302 = sin(qJ(5));
t305 = cos(qJ(5));
t301 = sin(qJ(6));
t465 = t301 * t302;
t304 = cos(qJ(6));
t467 = t296 * t304;
t174 = t295 * t465 - t467;
t167 = Icges(7,4) * t174;
t464 = t302 * t304;
t469 = t296 * t301;
t175 = t295 * t464 + t469;
t470 = t295 * t305;
t101 = -Icges(7,1) * t175 + Icges(7,5) * t470 + t167;
t168 = Icges(7,4) * t175;
t97 = -Icges(7,2) * t174 - Icges(7,6) * t470 + t168;
t555 = t101 * t304 + t301 * t97;
t95 = -Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t470;
t38 = -t302 * t95 - t305 * t555;
t307 = qJD(1) ^ 2;
t560 = t307 * t299;
t502 = rSges(7,1) * t304;
t385 = -rSges(7,2) * t301 + t502;
t552 = t305 * t385;
t368 = t175 * t101 + t174 * t97;
t471 = t295 * t304;
t176 = -t296 * t465 - t471;
t169 = Icges(7,4) * t176;
t473 = t295 * t301;
t177 = t296 * t464 - t473;
t466 = t296 * t305;
t102 = Icges(7,1) * t177 - Icges(7,5) * t466 + t169;
t484 = Icges(7,4) * t177;
t99 = Icges(7,2) * t176 - Icges(7,6) * t466 + t484;
t503 = t177 * t102 + t176 * t99;
t96 = Icges(7,5) * t177 + Icges(7,6) * t176 - Icges(7,3) * t466;
t559 = t368 + t503 + (-t295 * t95 - t296 * t96) * t305;
t504 = -t177 * t101 + t176 * t97;
t505 = t175 * t102 - t174 * t99;
t558 = t504 + t305 * (-t295 * t96 + t296 * t95) + t505;
t218 = -rSges(4,2) * t296 + t295 * rSges(4,3);
t547 = t296 * pkin(2) + t295 * qJ(3);
t412 = t299 + t547;
t154 = t218 + t412;
t103 = t175 * rSges(7,1) - t174 * rSges(7,2) - rSges(7,3) * t470;
t193 = rSges(7,3) * t302 + t552;
t427 = qJD(6) * t305;
t431 = qJD(5) * t296;
t201 = -t295 * t427 + t431;
t510 = t305 * pkin(5);
t256 = pkin(8) * t302 + t510;
t428 = qJD(6) * t302;
t284 = qJD(1) + t428;
t557 = -t103 * t284 + t193 * t201 + t256 * t431;
t360 = Icges(7,5) * t304 - Icges(7,6) * t301;
t179 = Icges(7,3) * t302 + t305 * t360;
t482 = Icges(7,4) * t304;
t362 = -Icges(7,2) * t301 + t482;
t181 = Icges(7,6) * t302 + t305 * t362;
t483 = Icges(7,4) * t301;
t364 = Icges(7,1) * t304 - t483;
t183 = Icges(7,5) * t302 + t305 * t364;
t323 = t174 * t181 - t175 * t183 + t179 * t470;
t553 = t323 * t284;
t430 = qJD(5) * t302;
t409 = t296 * t430;
t433 = qJD(1) * t305;
t411 = t295 * t433;
t330 = t409 + t411;
t219 = t296 * rSges(3,1) - rSges(3,2) * t295;
t199 = t219 + t299;
t423 = qJD(1) * qJD(4);
t548 = qJDD(4) * t296 - 0.2e1 * t295 * t423;
t213 = t295 * rSges(5,2) + t296 * rSges(5,3);
t298 = t305 * pkin(8);
t514 = pkin(5) * t302;
t546 = t298 - t514;
t432 = qJD(5) * t295;
t200 = -t296 * t427 - t432;
t27 = t470 * t95 - t368;
t28 = -t470 * t96 + t505;
t10 = t200 * t28 + t201 * t27 - t553;
t29 = t466 * t95 + t504;
t30 = -t466 * t96 + t503;
t56 = t176 * t181 + t177 * t183 - t179 * t466;
t54 = t56 * t284;
t11 = t200 * t30 + t201 * t29 + t54;
t545 = t10 * t295 + t11 * t296;
t361 = Icges(6,5) * t302 + Icges(6,6) * t305;
t155 = Icges(6,3) * t296 + t295 * t361;
t544 = qJD(1) * t155;
t424 = qJD(1) * qJD(3);
t436 = qJD(1) * t295;
t435 = qJD(1) * t296;
t268 = qJ(3) * t435;
t282 = qJD(3) * t295;
t444 = t268 + t282;
t542 = qJD(1) * (-pkin(2) * t436 + t444) + qJDD(1) * t547 + t295 * t424;
t286 = t296 * qJ(4);
t540 = qJDD(1) * t286 + qJDD(4) * t295 + 0.2e1 * t296 * t423;
t486 = Icges(6,4) * t302;
t363 = Icges(6,2) * t305 + t486;
t157 = Icges(6,6) * t296 + t295 * t363;
t485 = Icges(6,4) * t305;
t365 = Icges(6,1) * t302 + t485;
t159 = Icges(6,5) * t296 + t295 * t365;
t386 = rSges(6,1) * t302 + rSges(6,2) * t305;
t161 = rSges(6,3) * t296 + t295 * t386;
t239 = Icges(6,5) * t305 - Icges(6,6) * t302;
t335 = qJD(5) * t239;
t158 = -Icges(6,6) * t295 + t296 * t363;
t249 = Icges(6,4) * t466;
t468 = t296 * t302;
t481 = Icges(6,5) * t295;
t160 = Icges(6,1) * t468 + t249 - t481;
t356 = t158 * t305 + t160 * t302;
t538 = qJD(1) * t356 + t296 * t335 - t544;
t357 = t157 * t305 + t159 * t302;
t156 = -Icges(6,3) * t295 + t296 * t361;
t438 = qJD(1) * t156;
t537 = qJD(1) * t357 + t295 * t335 + t438;
t241 = -Icges(6,2) * t302 + t485;
t243 = Icges(6,1) * t305 - t486;
t353 = t241 * t305 + t243 * t302;
t536 = t353 * qJD(1) - t361 * qJD(5);
t452 = -Icges(6,2) * t468 + t160 + t249;
t454 = -t243 * t296 + t158;
t535 = -t302 * t454 + t305 * t452;
t178 = Icges(7,3) * t305 - t302 * t360;
t354 = t181 * t301 - t183 * t304;
t366 = t102 * t304 - t301 * t99;
t534 = t200 * (-t179 * t296 + t366) + t201 * (-t179 * t295 - t555) - t284 * (t178 + t354);
t209 = (-Icges(7,2) * t304 - t483) * t305;
t533 = t200 * (-Icges(7,2) * t177 + t102 + t169) + t201 * (-Icges(7,2) * t175 - t101 - t167) + t284 * (t183 + t209);
t322 = -qJDD(6) * t305 + (-qJD(1) + t428) * qJD(5);
t407 = qJD(1) * t427;
t118 = (-qJDD(5) + t407) * t295 + t322 * t296;
t532 = t118 / 0.2e1;
t276 = qJDD(5) * t296;
t119 = t295 * t322 - t296 * t407 + t276;
t531 = t119 / 0.2e1;
t530 = -t200 / 0.2e1;
t529 = t200 / 0.2e1;
t528 = -t201 / 0.2e1;
t527 = t201 / 0.2e1;
t422 = qJD(1) * qJD(5);
t203 = -qJDD(5) * t295 - t296 * t422;
t526 = t203 / 0.2e1;
t204 = -t295 * t422 + t276;
t525 = t204 / 0.2e1;
t524 = -t284 / 0.2e1;
t523 = t284 / 0.2e1;
t521 = t295 / 0.2e1;
t520 = -t296 / 0.2e1;
t518 = -rSges(6,3) - pkin(7);
t517 = -rSges(7,3) - pkin(8);
t515 = pkin(2) * t295;
t513 = pkin(7) * t296;
t512 = pkin(7) * t307;
t511 = g(1) * t295;
t410 = t296 * t433;
t331 = t295 * t430 - t410;
t429 = qJD(5) * t305;
t318 = -t284 * t304 - t301 * t429;
t434 = qJD(1) * t302;
t395 = qJD(6) + t434;
t91 = t295 * t318 - t395 * t469;
t317 = -t284 * t301 + t304 * t429;
t92 = t295 * t317 + t395 * t467;
t45 = Icges(7,5) * t92 + Icges(7,6) * t91 + Icges(7,3) * t331;
t47 = Icges(7,4) * t92 + Icges(7,2) * t91 + Icges(7,6) * t331;
t49 = Icges(7,1) * t92 + Icges(7,4) * t91 + Icges(7,5) * t331;
t7 = (qJD(5) * t555 + t45) * t302 + (-qJD(5) * t95 - t301 * t47 + t304 * t49 + (t101 * t301 - t304 * t97) * qJD(6)) * t305;
t509 = t7 * t201;
t89 = t296 * t318 + t395 * t473;
t90 = t296 * t317 - t395 * t471;
t44 = Icges(7,5) * t90 + Icges(7,6) * t89 + Icges(7,3) * t330;
t46 = Icges(7,4) * t90 + Icges(7,2) * t89 + Icges(7,6) * t330;
t48 = Icges(7,1) * t90 + Icges(7,4) * t89 + Icges(7,5) * t330;
t8 = (-qJD(5) * t366 + t44) * t302 + (qJD(5) * t96 - t301 * t46 + t304 * t48 + (-t102 * t301 - t304 * t99) * qJD(6)) * t305;
t508 = t8 * t200;
t507 = -pkin(2) - qJ(4);
t222 = qJD(5) * t427 + qJDD(6) * t302 + qJDD(1);
t208 = (-Icges(7,5) * t301 - Icges(7,6) * t304) * t305;
t132 = qJD(5) * t178 + qJD(6) * t208;
t180 = Icges(7,6) * t305 - t302 * t362;
t133 = qJD(5) * t180 + qJD(6) * t209;
t182 = Icges(7,5) * t305 - t302 * t364;
t210 = (-Icges(7,1) * t301 - t482) * t305;
t134 = qJD(5) * t182 + qJD(6) * t210;
t24 = (qJD(5) * t354 + t132) * t302 + (qJD(5) * t179 - t133 * t301 + t134 * t304 + (-t181 * t304 - t183 * t301) * qJD(6)) * t305;
t70 = t179 * t302 - t305 * t354;
t506 = t70 * t222 + t24 * t284;
t498 = rSges(4,3) * t296;
t497 = rSges(6,3) * t295;
t451 = t177 * rSges(7,1) + t176 * rSges(7,2);
t105 = -rSges(7,3) * t466 + t451;
t408 = t296 * t429;
t415 = pkin(5) * t408 + pkin(8) * t330;
t472 = t295 * t302;
t419 = pkin(5) * t472;
t130 = -qJD(1) * t419 + t415;
t211 = (-rSges(7,1) * t301 - rSges(7,2) * t304) * t305;
t297 = t305 * rSges(7,3);
t135 = qJD(6) * t211 + (-t302 * t385 + t297) * qJD(5);
t263 = pkin(5) * t468;
t196 = -pkin(8) * t466 + t263;
t294 = qJDD(1) * t299;
t391 = -t307 * t516 + t294;
t332 = t391 + t542;
t309 = (-qJDD(3) - t512) * t296 + t332 + t540;
t229 = t546 * qJD(5);
t476 = qJ(4) * t307;
t350 = qJD(5) * t229 - t476;
t477 = pkin(7) * qJDD(1);
t50 = t90 * rSges(7,1) + t89 * rSges(7,2) + rSges(7,3) * t330;
t12 = t309 + (t350 - t477) * t295 + qJD(1) * t130 + qJDD(1) * t196 + t105 * t222 - t118 * t193 - t135 * t200 - t203 * t256 + t284 * t50;
t493 = t12 * t296;
t131 = t331 * pkin(8) + (t295 * t429 + t296 * t434) * pkin(5);
t283 = qJD(3) * t296;
t164 = qJD(1) * t547 - t283;
t260 = pkin(8) * t470;
t194 = -t260 + t419;
t287 = t296 * qJ(3);
t212 = -t287 + t515;
t383 = -qJ(4) * t295 - t516;
t333 = t383 - t513;
t326 = -t212 + t333;
t320 = -t194 + t326;
t445 = qJDD(3) * t295 + t296 * t424;
t349 = t445 - t560;
t324 = t295 * t512 + t349 + t548;
t389 = rSges(7,1) * t92 + rSges(7,2) * t91;
t51 = rSges(7,3) * t331 + t389;
t13 = -t103 * t222 + t119 * t193 + t135 * t201 + t204 * t256 - t284 * t51 + t350 * t296 + (-t131 - t164) * qJD(1) + t320 * qJDD(1) + t324;
t492 = t13 * t295;
t247 = rSges(6,1) * t305 - rSges(6,2) * t302;
t202 = t247 * t431;
t321 = -t161 + t326;
t281 = qJD(4) * t296;
t440 = t281 + t282;
t67 = qJD(1) * t321 + t202 + t440;
t491 = t295 * t67;
t488 = t38 * t119;
t39 = t302 * t96 + t305 * t366;
t487 = t39 * t118;
t475 = t155 * t295;
t474 = t155 * t296;
t184 = t239 * t295;
t185 = t239 * t296;
t152 = t296 * t156;
t460 = -t103 - t194;
t459 = t105 + t196;
t458 = t135 + t229;
t457 = t157 * t466 + t159 * t468;
t456 = t158 * t466 + t160 * t468;
t455 = t243 * t295 - t157;
t453 = t241 * t295 + t159;
t449 = t193 + t256;
t448 = -t363 + t243;
t447 = -t241 - t365;
t446 = rSges(6,1) * t468 + rSges(6,2) * t466;
t195 = pkin(5) * t470 + pkin(8) * t472;
t197 = pkin(5) * t466 + pkin(8) * t468;
t443 = rSges(4,2) * t436 + rSges(4,3) * t435;
t442 = pkin(7) * t436 + t283;
t441 = rSges(7,2) * t465 + t297;
t439 = -qJD(1) * t212 + t282;
t437 = qJD(1) * t361;
t114 = t295 * t353 + t185;
t426 = t114 * qJD(1);
t425 = -m(5) - m(6) - m(7);
t421 = qJDD(3) * t296;
t420 = -rSges(5,3) + t507;
t416 = -m(4) + t425;
t64 = t158 * t470 + t160 * t472 + t152;
t414 = t268 + t440;
t413 = t281 + t439;
t405 = -t435 / 0.2e1;
t404 = -t432 / 0.2e1;
t403 = t432 / 0.2e1;
t402 = -t431 / 0.2e1;
t401 = t431 / 0.2e1;
t398 = rSges(4,2) * t295 + t498 - t516;
t397 = t287 - t516;
t217 = rSges(5,2) * t296 - t295 * rSges(5,3);
t396 = -qJD(4) * t295 + t283;
t394 = qJD(5) * t247 + qJD(4);
t393 = t286 + t412;
t392 = t507 - t514;
t248 = rSges(2,1) * t306 - rSges(2,2) * t303;
t246 = rSges(2,1) * t303 + rSges(2,2) * t306;
t382 = -t286 - t299;
t336 = qJD(5) * t241;
t110 = -qJD(1) * t157 + t296 * t336;
t337 = qJD(5) * t243;
t112 = -qJD(1) * t159 + t296 * t337;
t84 = -t158 * t302 + t160 * t305;
t311 = qJD(5) * t84 + t110 * t305 + t112 * t302 - t438;
t111 = qJD(1) * t158 + t295 * t336;
t113 = t295 * t337 + (t296 * t365 - t481) * qJD(1);
t83 = -t157 * t302 + t159 * t305;
t312 = qJD(5) * t83 + t111 * t305 + t113 * t302 - t544;
t380 = (-t295 * t537 + t312 * t296) * t296 - (-t295 * t538 + t311 * t296) * t295;
t379 = (t312 * t295 + t296 * t537) * t296 - (t311 * t295 + t296 * t538) * t295;
t190 = t247 * t295;
t117 = qJD(5) * t190 + (t296 * t386 - t497) * qJD(1);
t226 = t386 * qJD(5);
t351 = -qJD(5) * t226 - t476;
t25 = t204 * t247 + t351 * t296 + (-t117 - t164) * qJD(1) + t321 * qJDD(1) + t324;
t348 = rSges(6,1) * t408 - rSges(6,2) * t409;
t116 = -qJD(1) * t161 + t348;
t162 = t446 - t497;
t26 = qJD(1) * t116 + qJDD(1) * t162 - t203 * t247 + (t351 - t477) * t295 + t309;
t378 = t25 * t296 + t26 * t295;
t377 = t27 * t296 - t28 * t295;
t376 = t27 * t295 + t28 * t296;
t375 = t29 * t296 - t295 * t30;
t374 = t29 * t295 + t296 * t30;
t373 = t295 * t39 - t296 * t38;
t372 = t295 * t38 + t296 * t39;
t63 = t295 * t357 + t474;
t371 = -t295 * t64 + t296 * t63;
t65 = t457 - t475;
t66 = -t156 * t295 + t456;
t370 = -t295 * t66 + t296 * t65;
t347 = t547 - t382;
t68 = t394 * t295 + (t162 + t347) * qJD(1) - t442;
t369 = t295 * t68 + t296 * t67;
t359 = t103 * t296 - t105 * t295;
t358 = -t116 * t296 - t117 * t295;
t355 = -t161 * t295 - t162 * t296;
t352 = -t241 * t302 + t243 * t305;
t346 = t217 + t383;
t345 = -t305 * t44 + t430 * t96;
t344 = -t305 * t45 - t430 * t95;
t150 = rSges(7,3) * t472 + t295 * t552;
t151 = rSges(7,3) * t468 + t296 * t552;
t338 = -t386 + t507;
t334 = -t212 + t346;
t329 = t179 * t284 + t200 * t96 - t201 * t95;
t327 = (-Icges(7,5) * t174 - Icges(7,6) * t175) * t201 + (Icges(7,5) * t176 - Icges(7,6) * t177) * t200 + t208 * t284;
t325 = t302 * t455 + t305 * t453;
t319 = t305 * t327;
t316 = (t302 * t447 + t305 * t448) * qJD(1);
t315 = (Icges(7,1) * t176 - t484 - t99) * t200 + (-Icges(7,1) * t174 - t168 - t97) * t201 + (-t181 + t210) * t284;
t32 = t103 * t200 - t105 * t201 + qJD(2) + (-t194 * t295 - t196 * t296) * qJD(5);
t36 = qJD(1) * t320 + t440 + t557;
t37 = t105 * t284 - t193 * t200 + (qJD(5) * t256 + qJD(4)) * t295 + (t196 + t347) * qJD(1) - t442;
t313 = t32 * t359 + (t295 * t36 - t296 * t37) * t193;
t224 = t363 * qJD(5);
t225 = t365 * qJD(5);
t310 = -qJD(1) * t239 + qJD(5) * t352 - t224 * t305 - t225 * t302;
t308 = t534 * t305;
t273 = rSges(5,2) * t435;
t192 = -rSges(7,1) * t464 + t441;
t191 = t247 * t296;
t149 = t183 * t296;
t148 = t183 * t295;
t147 = t181 * t296;
t146 = t181 * t295;
t129 = rSges(7,1) * t176 - rSges(7,2) * t177;
t128 = -rSges(7,1) * t174 - rSges(7,2) * t175;
t115 = t296 * t353 - t184;
t93 = t115 * qJD(1);
t82 = (t213 + t347) * qJD(1) - t396;
t81 = qJD(1) * t334 + t440;
t77 = qJD(5) * t355 + qJD(2);
t58 = qJD(1) * t443 + qJDD(1) * t218 + t332 - t421;
t57 = (-t212 + t398) * qJDD(1) + (-qJD(1) * t218 - t164) * qJD(1) + t349;
t53 = t294 - t421 + qJDD(1) * t213 + qJD(1) * (-rSges(5,3) * t436 + t273) + t383 * t307 + t540 + t542;
t52 = t382 * t307 + (-qJD(1) * t213 - t164) * qJD(1) + t334 * qJDD(1) + t445 + t548;
t43 = t310 * t295 + t296 * t536;
t42 = -t295 * t536 + t310 * t296;
t41 = -qJD(5) * t356 - t110 * t302 + t112 * t305;
t40 = -qJD(5) * t357 - t111 * t302 + t113 * t305;
t31 = qJD(5) * t358 + t161 * t203 - t162 * t204 + qJDD(2);
t23 = qJD(5) * t370 + t93;
t22 = qJD(5) * t371 + t426;
t16 = -t132 * t470 - t133 * t174 + t134 * t175 + t179 * t331 + t181 * t91 + t183 * t92;
t15 = -t132 * t466 + t133 * t176 + t134 * t177 + t179 * t330 + t181 * t89 + t183 * t90;
t14 = t200 * t39 + t201 * t38 + t284 * t70;
t9 = t103 * t118 - t105 * t119 + t194 * t203 - t196 * t204 + t200 * t51 - t201 * t50 + qJDD(2) + (-t130 * t296 - t131 * t295) * qJD(5);
t6 = t102 * t92 - t174 * t46 + t175 * t48 + t295 * t345 - t410 * t96 + t91 * t99;
t5 = -t101 * t92 - t174 * t47 + t175 * t49 + t295 * t344 + t410 * t95 + t91 * t97;
t4 = t102 * t90 + t176 * t46 + t177 * t48 + t296 * t345 + t411 * t96 + t89 * t99;
t3 = -t101 * t90 + t176 * t47 + t177 * t49 + t296 * t344 - t411 * t95 + t89 * t97;
t2 = t118 * t28 + t119 * t27 + t16 * t284 + t200 * t6 + t201 * t5 - t222 * t323;
t1 = t118 * t30 + t119 * t29 + t15 * t284 + t200 * t4 + t201 * t3 + t222 * t56;
t17 = [-t323 * t531 + (t67 * t442 + t68 * (t348 + t414) - t394 * t491 + ((-t303 * t68 - t306 * t67) * pkin(1) + (t338 * t67 + t518 * t68) * t296 + (t67 * (rSges(6,3) - qJ(3)) + t68 * t338) * t295) * qJD(1) - (-t67 + t202 + t413 + (-t161 + t333) * qJD(1)) * t68 + (t26 - g(2)) * (t295 * t518 + t393 + t446) + (t25 - g(1)) * (t338 * t295 + t518 * t296 + t397)) * m(6) + (t81 * t396 + t82 * (t273 + t414) + ((-t303 * t82 - t306 * t81) * pkin(1) + t81 * t420 * t296 + (t81 * (-rSges(5,2) - qJ(3)) + t82 * t420) * t295) * qJD(1) - (qJD(1) * t346 + t413 - t81) * t82 + (-g(2) + t53) * (t393 + t213) + (-g(1) + t52) * (t295 * t507 + t217 + t397)) * m(5) + (t36 * (-t389 + t442) + t37 * (t50 + t414 + t415) + (t13 * t392 + (-qJD(4) + (t302 * t517 - t510) * qJD(5)) * t36) * t295 + ((-t303 * t37 - t306 * t36) * pkin(1) + (-t36 * qJ(3) + t37 * t392) * t295 + (t36 * (t546 + t297 + t507) - t37 * pkin(7)) * t296) * qJD(1) - t392 * t511 - (-t36 + t413 + (-t194 + t333) * qJD(1) + t557) * t37 + (t12 - g(2)) * (-pkin(7) * t295 + t466 * t517 + t263 + t393 + t451) + (t13 - g(1)) * (t260 + t397 - t513 - t103)) * m(7) + (t114 + t83) * t525 + t487 / 0.2e1 + t488 / 0.2e1 + ((t57 - g(1)) * (t498 + (rSges(4,2) - pkin(2)) * t295 + t397) + (t58 - g(2)) * t154 + (t443 + t444 - t439 + (-t515 - t516 - t398) * qJD(1)) * (qJD(1) * t154 - t283)) * m(4) + (t41 + t42) * t404 + (t22 - t426 + ((t456 - t66 - t474) * t296 + (-t152 + t64 - t65 - t475) * t295) * qJD(5)) * t403 + (t16 + t11) * t527 + t508 / 0.2e1 + t509 / 0.2e1 + (-qJD(5) * t353 + t224 * t302 - t225 * t305) * qJD(1) + (t93 + (t457 * t296 + (-t63 + (t156 + t357) * t295 - t456) * t295) * qJD(5)) * t402 + (t352 + m(3) * (t198 ^ 2 + t219 * t199) + m(2) * (t246 ^ 2 + t248 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,1) + Icges(5,1)) * qJDD(1) + (t43 + t23 + t40) * t401 + (t553 + (-t30 + t559) * t201 + (t29 - t558) * t200 + t10) * t530 + (t54 + (-t28 + t558) * t201 + (t27 + t559) * t200) * t528 + ((-t215 * t307 - g(2) + t391) * t199 + (-t560 + (-0.2e1 * t219 - t299 + t199) * t307 - g(1)) * t198) * m(3) + (t115 + t84) * t526 + t506 + t15 * t529 + t56 * t532 - m(2) * (-g(1) * t246 + g(2) * t248); m(6) * t31 + m(7) * t9 + (m(3) + m(4) + m(5)) * qJDD(2) + (-m(3) + t416) * g(3); t416 * (-g(2) * t296 + t511) + 0.2e1 * (-t493 / 0.2e1 + t492 / 0.2e1) * m(7) + 0.2e1 * (t25 * t521 + t26 * t520) * m(6) + 0.2e1 * (t52 * t521 + t520 * t53) * m(5) + 0.2e1 * (t520 * t58 + t521 * t57) * m(4); t425 * (g(1) * t296 + g(2) * t295) + m(5) * (t295 * t53 + t296 * t52) + m(6) * t378 + m(7) * (t12 * t295 + t13 * t296); ((-t147 * t174 + t149 * t175) * t200 + (-t146 * t174 + t148 * t175) * t201 + (-t174 * t180 + t175 * t182) * t284 + (t28 * t468 - t305 * t323) * qJD(6) + ((qJD(6) * t27 + t329) * t302 + t308) * t295) * t528 + (((-t147 * t301 + t149 * t304 + t96) * t200 + (-t146 * t301 + t148 * t304 - t95) * t201 + (-t180 * t301 + t182 * t304 + t179) * t284 + t70 * qJD(6)) * t305 + (qJD(6) * t372 - t534) * t302) * t524 + (-(-t190 * t67 + t191 * t68) * qJD(1) - (t77 * (-t190 * t295 - t191 * t296) - t369 * t386) * qJD(5) + t31 * t355 + t77 * ((-t161 * t296 + t162 * t295) * qJD(1) + t358) - t369 * t226 + ((t296 * t68 - t491) * qJD(1) + t378) * t247 - g(1) * t191 - g(2) * t190 + g(3) * t386) * m(6) + ((t13 * t449 + t36 * t458 - t9 * t459 + t32 * (-t130 - t50) + (t32 * t460 + t37 * t449) * qJD(1)) * t296 + (t12 * t449 + t37 * t458 + t9 * t460 + t32 * (-t131 - t51) + (t32 * t459 - t36 * t449) * qJD(1)) * t295 - t36 * (-qJD(1) * t195 - t150 * t284 + t192 * t201 + t431 * t546) - t37 * (qJD(1) * t197 + t151 * t284 - t192 * t200 + t432 * t546) - t32 * (t150 * t200 - t151 * t201 - t195 * t432 - t197 * t431) - ((-t103 * t36 + t105 * t37) * t305 + t313 * t302) * qJD(6) - g(1) * (t151 + t197) - g(2) * (t150 + t195) - g(3) * (t298 + (-pkin(5) - t502) * t302 + t441)) * m(7) + qJD(1) * (-t295 * t41 + t296 * t40 + (-t295 * t83 - t296 * t84) * qJD(1)) / 0.2e1 + (t23 + t11) * t405 - t545 * t428 / 0.2e1 + (-qJD(1) * t372 - t295 * t8 + t296 * t7) * t523 + t371 * t525 - (t22 + t10) * t436 / 0.2e1 + (qJD(1) * t43 + qJD(5) * t379 + qJDD(1) * t114 + t203 * t64 + t204 * t63 + t2) * t296 / 0.2e1 - (qJD(1) * t42 + qJD(5) * t380 + qJDD(1) * t115 + t203 * t66 + t204 * t65 + t1) * t295 / 0.2e1 - qJD(1) * ((-t448 * t302 + t447 * t305) * qJD(1) + ((t295 * t454 + t296 * t455) * t305 + (t295 * t452 - t296 * t453) * t302) * qJD(5)) / 0.2e1 + ((t185 * t432 + t437) * t295 + (t316 + (t325 * t296 + (-t184 - t535) * t295) * qJD(5)) * t296) * t403 + ((t184 * t431 - t437) * t296 + (t316 + (-t535 * t295 + (-t185 + t325) * t296) * qJD(5)) * t295) * t402 + ((-t65 * t295 - t66 * t296) * qJD(1) + t380) * t404 - t222 * t373 / 0.2e1 + ((-t63 * t295 - t64 * t296) * qJD(1) + t379) * t401 - t14 * t427 / 0.2e1 + qJDD(1) * (-t295 * t84 + t296 * t83) / 0.2e1 + t370 * t526 + (-qJD(1) * t376 - t295 * t6 + t296 * t5) * t527 + (-qJD(1) * t374 - t295 * t4 + t296 * t3) * t529 + ((t147 * t176 + t149 * t177) * t200 + (t146 * t176 + t148 * t177) * t201 + (t176 * t180 + t177 * t182) * t284 + (t29 * t472 + t305 * t56) * qJD(6) + ((qJD(6) * t30 + t329) * t302 + t308) * t296) * t530 + t377 * t531 + t375 * t532; t11 * t411 / 0.2e1 - t1 * t466 / 0.2e1 + (t302 * t56 - t305 * t374) * t532 + ((qJD(5) * t374 + t15) * t302 + (-qJD(1) * t375 + qJD(5) * t56 - t295 * t3 - t296 * t4) * t305) * t529 + t305 * t10 * t405 - t2 * t470 / 0.2e1 + (-t302 * t323 - t305 * t376) * t531 + ((qJD(5) * t376 + t16) * t302 + (-qJD(1) * t377 - qJD(5) * t323 - t295 * t5 - t296 * t6) * t305) * t527 + t14 * t429 / 0.2e1 + t302 * (t487 + t488 + t506 + t508 + t509) / 0.2e1 + t222 * (t302 * t70 - t305 * t372) / 0.2e1 + ((qJD(5) * t372 + t24) * t302 + (qJD(1) * t373 + qJD(5) * t70 - t295 * t7 - t296 * t8) * t305) * t523 + (t176 * t533 + t315 * t177 - t296 * t319) * t530 + (-t174 * t533 + t175 * t315 - t295 * t319) * t528 + (t327 * t302 + (-t301 * t533 + t304 * t315) * t305) * t524 + t545 * t430 / 0.2e1 + ((qJD(5) * t313 - t13 * t103 + t12 * t105 - t36 * t51 + t37 * t50) * t302 + (t36 * (-qJD(5) * t103 - t135 * t295) + t37 * (qJD(5) * t105 + t135 * t296) - t9 * t359 + t32 * (t103 * t436 + t105 * t435 + t295 * t50 - t296 * t51) + (t493 - t492 + (-t295 * t37 - t296 * t36) * qJD(1)) * t193) * t305 - t36 * (-t128 * t284 + t201 * t211) - t37 * (t129 * t284 - t200 * t211) - t32 * (t128 * t200 - t129 * t201) - g(1) * t129 - g(2) * t128 - g(3) * t211) * m(7);];
tau  = t17;
