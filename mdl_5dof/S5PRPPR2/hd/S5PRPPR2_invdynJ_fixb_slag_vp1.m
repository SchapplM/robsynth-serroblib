% Calculate vector of inverse dynamics joint torques for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:58
% EndTime: 2019-12-05 15:24:40
% DurationCPUTime: 31.43s
% Computational Cost: add. (13228->646), mult. (15385->960), div. (0->0), fcn. (14289->10), ass. (0->343)
t300 = qJ(2) + pkin(8);
t293 = sin(t300);
t295 = cos(t300);
t390 = -Icges(4,5) * t293 - Icges(4,6) * t295;
t307 = sin(qJ(2));
t308 = cos(qJ(2));
t586 = Icges(3,5) * t307 + Icges(3,6) * t308;
t581 = t390 - t586;
t589 = -0.2e1 * t307 * (Icges(3,1) - Icges(3,2)) * t308 + (0.2e1 * t307 ^ 2 - 0.2e1 * t308 ^ 2) * Icges(3,4);
t585 = Icges(3,3) + Icges(4,3);
t584 = t581 * qJD(2);
t583 = Icges(3,5) * t308 + Icges(4,5) * t295 - Icges(3,6) * t307 - Icges(4,6) * t293;
t302 = sin(pkin(7));
t304 = cos(pkin(7));
t582 = (t583 * t302 - t585 * t304) * t304;
t297 = t302 ^ 2;
t298 = t304 ^ 2;
t547 = t297 + t298;
t278 = rSges(3,1) * t307 + rSges(3,2) * t308;
t580 = t278 * t547;
t579 = t585 * t302 + t583 * t304;
t578 = t584 * t304;
t576 = t584 * t302;
t299 = pkin(9) + qJ(5);
t292 = sin(t299);
t518 = rSges(6,2) * t292;
t294 = cos(t299);
t521 = rSges(6,1) * t294;
t418 = -t518 + t521;
t159 = t293 * rSges(6,3) + t295 * t418;
t301 = sin(pkin(9));
t519 = rSges(5,2) * t301;
t303 = cos(pkin(9));
t522 = rSges(5,1) * t303;
t419 = -t519 + t522;
t575 = t293 * rSges(5,3) + t295 * t419;
t570 = 0.2e1 * qJD(2);
t569 = 2 * qJDD(2);
t566 = t581 * t302;
t565 = t581 * t304;
t489 = t303 * t304;
t492 = t301 * t302;
t223 = -t295 * t492 - t489;
t490 = t302 * t303;
t491 = t301 * t304;
t224 = t295 * t490 - t491;
t497 = t293 * t302;
t92 = Icges(5,5) * t224 + Icges(5,6) * t223 + Icges(5,3) * t497;
t225 = -t295 * t491 + t490;
t226 = t295 * t489 + t492;
t496 = t293 * t304;
t93 = Icges(5,5) * t226 + Icges(5,6) * t225 + Icges(5,3) * t496;
t94 = Icges(5,4) * t224 + Icges(5,2) * t223 + Icges(5,6) * t497;
t95 = Icges(5,4) * t226 + Icges(5,2) * t225 + Icges(5,6) * t496;
t96 = Icges(5,1) * t224 + Icges(5,4) * t223 + Icges(5,5) * t497;
t97 = Icges(5,1) * t226 + Icges(5,4) * t225 + Icges(5,5) * t496;
t564 = (-(t301 * t94 - t303 * t96) * t293 - t92 * t295 + t390 * t304) * t304 + ((t301 * t95 - t303 * t97) * t293 + t93 * t295 + t390 * t302) * t302;
t466 = qJD(4) * t302;
t266 = t293 * t466;
t528 = pkin(3) * t293;
t253 = -qJ(4) * t295 + t528;
t363 = qJD(2) * t253;
t145 = -t302 * t363 + t266;
t465 = qJD(4) * t304;
t268 = t293 * t465;
t146 = -t304 * t363 + t268;
t495 = t295 * t302;
t271 = qJ(4) * t495;
t494 = t295 * t304;
t272 = qJ(4) * t494;
t470 = qJD(2) * t304;
t472 = qJD(2) * t302;
t516 = pkin(2) * qJD(2);
t450 = t307 * t516;
t468 = qJD(3) * t304;
t251 = -t302 * t450 - t468;
t291 = qJD(3) * t302;
t252 = -t304 * t450 + t291;
t481 = t302 * t251 + t304 * t252;
t563 = t302 * t145 + t304 * t146 - (-pkin(3) * t497 + t271) * t472 - (-pkin(3) * t496 + t272) * t470 - qJD(4) * t293 + t481;
t461 = qJD(2) * qJD(4);
t562 = qJDD(4) * t293 + t295 * t461;
t474 = qJD(2) * t293;
t473 = qJD(2) * t295;
t545 = rSges(5,3) * t295 - t293 * t419;
t559 = t547 * qJD(2) * t545;
t296 = t308 * pkin(2);
t548 = t295 * pkin(3) + t293 * qJ(4);
t433 = -t548 - t296;
t549 = t295 * rSges(4,1) - rSges(4,2) * t293;
t550 = t549 + t296;
t546 = g(1) * t304 + g(2) * t302;
t389 = Icges(6,5) * t294 - Icges(6,6) * t292;
t147 = -Icges(6,3) * t295 + t293 * t389;
t502 = Icges(6,4) * t294;
t394 = -Icges(6,2) * t292 + t502;
t149 = -Icges(6,6) * t295 + t293 * t394;
t503 = Icges(6,4) * t292;
t399 = Icges(6,1) * t294 - t503;
t151 = -Icges(6,5) * t295 + t293 * t399;
t464 = qJD(5) * t293;
t245 = t304 * t464 + t472;
t246 = t302 * t464 - t470;
t488 = t304 * t292;
t202 = t294 * t302 - t295 * t488;
t203 = t292 * t302 + t294 * t494;
t504 = Icges(6,4) * t203;
t84 = Icges(6,2) * t202 + Icges(6,6) * t496 + t504;
t184 = Icges(6,4) * t202;
t86 = Icges(6,1) * t203 + Icges(6,5) * t496 + t184;
t411 = -t292 * t84 + t294 * t86;
t200 = -t292 * t495 - t294 * t304;
t201 = t294 * t495 - t488;
t505 = Icges(6,4) * t201;
t83 = Icges(6,2) * t200 + Icges(6,6) * t497 + t505;
t183 = Icges(6,4) * t200;
t85 = Icges(6,1) * t201 + Icges(6,5) * t497 + t183;
t412 = -t292 * t83 + t294 * t85;
t544 = -(-t147 * t304 - t411) * t245 - (-t147 * t302 - t412) * t246;
t306 = -pkin(6) - qJ(4);
t486 = qJ(4) + t306;
t288 = pkin(4) * t303 + pkin(3);
t524 = pkin(3) - t288;
t543 = t293 * t524 - t295 * t486;
t542 = t589 * t302 + t586 * t304;
t186 = (-Icges(6,2) * t294 - t503) * t293;
t463 = qJD(5) * t295;
t318 = t245 * (-Icges(6,2) * t203 + t184 + t86) + t246 * (-Icges(6,2) * t201 + t183 + t85) - t463 * (t151 + t186);
t309 = qJD(2) ^ 2;
t541 = -m(5) - m(6);
t460 = qJD(2) * qJD(5);
t350 = qJDD(5) * t293 + t295 * t460;
t459 = qJDD(2) * t302;
t162 = t304 * t350 + t459;
t540 = t162 / 0.2e1;
t458 = qJDD(2) * t304;
t163 = t302 * t350 - t458;
t539 = t163 / 0.2e1;
t538 = -t547 * t474 / 0.2e1;
t238 = -qJDD(5) * t295 + t293 * t460;
t537 = t238 / 0.2e1;
t536 = -t245 / 0.2e1;
t535 = t245 / 0.2e1;
t534 = -t246 / 0.2e1;
t533 = t246 / 0.2e1;
t532 = -t295 / 0.2e1;
t529 = pkin(2) * t307;
t82 = Icges(6,5) * t203 + Icges(6,6) * t202 + Icges(6,3) * t496;
t25 = t200 * t84 + t201 * t86 + t497 * t82;
t515 = t25 * t304;
t81 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t497;
t26 = t202 * t83 + t203 * t85 + t496 * t81;
t514 = t26 * t302;
t46 = t147 * t497 + t149 * t200 + t151 * t201;
t513 = t46 * t293;
t47 = t147 * t496 + t149 * t202 + t151 * t203;
t512 = t47 * t293;
t498 = t147 * t295;
t493 = t295 * t306;
t487 = t308 * t309;
t167 = -qJ(3) * t304 + t296 * t302;
t168 = qJ(3) * t302 + t296 * t304;
t484 = t302 * t167 + t304 * t168;
t451 = t293 * t518;
t480 = rSges(6,3) * t495 + t302 * t451;
t479 = rSges(6,3) * t494 + t304 * t451;
t452 = t293 * t519;
t478 = rSges(5,3) * t495 + t302 * t452;
t477 = rSges(5,3) * t494 + t304 * t452;
t476 = t268 + t291;
t467 = qJD(4) * t295;
t462 = -m(4) + t541;
t457 = qJDD(3) * t304;
t455 = pkin(2) * t487;
t454 = t293 * t522;
t453 = t293 * t521;
t449 = t308 * t516;
t290 = qJDD(3) * t302;
t447 = t304 * t562 + t290;
t446 = t293 * t472;
t445 = t293 * t470;
t444 = t295 * t472;
t443 = t295 * t470;
t442 = t167 * t472 + t168 * t470 + qJD(1);
t437 = t470 / 0.2e1;
t436 = -t463 / 0.2e1;
t435 = t463 / 0.2e1;
t434 = -t253 - t529;
t254 = rSges(4,1) * t293 + rSges(4,2) * t295;
t343 = -t254 - t529;
t431 = t295 * t288 - t293 * t306;
t430 = t266 - t468;
t219 = t548 * t302;
t220 = t548 * t304;
t428 = t302 * t219 + t304 * t220 + t484;
t427 = t547 * t529;
t426 = t434 + t543;
t424 = t434 + t545;
t188 = qJD(2) * t548 - t467;
t422 = -t188 - t449;
t237 = t549 * qJD(2);
t421 = -t237 - t449;
t279 = rSges(3,1) * t308 - rSges(3,2) * t307;
t158 = -rSges(6,3) * t295 + t293 * t418;
t340 = -t293 * t486 - t295 * t524;
t126 = t340 * qJD(2);
t314 = -t455 + (-t126 - t188) * qJD(2) + t426 * qJDD(2);
t405 = t302 * t562 - t457;
t112 = -qJD(5) * t203 + t292 * t445;
t113 = qJD(5) * t202 - t294 * t445;
t66 = rSges(6,1) * t113 + rSges(6,2) * t112 + rSges(6,3) * t443;
t189 = (-rSges(6,1) * t292 - rSges(6,2) * t294) * t293;
t80 = qJD(2) * t159 + qJD(5) * t189;
t88 = rSges(6,1) * t203 + rSges(6,2) * t202 + rSges(6,3) * t496;
t13 = -t158 * t162 + t238 * t88 - t245 * t80 + t302 * t314 - t463 * t66 + t405;
t110 = -qJD(5) * t201 + t292 * t446;
t111 = qJD(5) * t200 - t294 * t446;
t65 = rSges(6,1) * t111 + rSges(6,2) * t110 + rSges(6,3) * t444;
t87 = rSges(6,1) * t201 + rSges(6,2) * t200 + rSges(6,3) * t497;
t14 = t158 * t163 - t238 * t87 + t246 * t80 + t304 * t314 + t463 * t65 + t447;
t416 = -t13 * t304 + t14 * t302;
t24 = t200 * t83 + t201 * t85 + t497 * t81;
t415 = t24 * t302 + t515;
t414 = t82 * t245 + t81 * t246;
t27 = t202 * t84 + t203 * t86 + t496 * t82;
t413 = t27 * t304 + t514;
t30 = t293 * t412 - t295 * t81;
t31 = t293 * t411 - t295 * t82;
t410 = t30 * t302 + t31 * t304;
t407 = -t302 * t88 + t304 * t87;
t89 = -pkin(4) * t491 + t302 * t340;
t90 = pkin(4) * t492 + t304 * t340;
t406 = t302 * t89 + t304 * t90;
t404 = qJD(2) * t343;
t100 = rSges(5,1) * t224 + rSges(5,2) * t223 + rSges(5,3) * t497;
t101 = rSges(5,1) * t226 + rSges(5,2) * t225 + rSges(5,3) * t496;
t388 = t100 * t302 + t101 * t304;
t387 = -t149 * t292 + t151 * t294;
t181 = -rSges(4,3) * t304 + t302 * t549;
t182 = rSges(4,3) * t302 + t304 * t549;
t384 = t181 * t302 + t182 * t304;
t381 = t547 * t279;
t380 = qJD(2) * t580;
t379 = t167 * t459 + t168 * t458 + t251 * t472 + t252 * t470 + qJDD(1);
t378 = -t158 + t426;
t155 = t575 * qJD(2);
t377 = -t155 + t422;
t374 = qJD(2) * t426;
t373 = qJD(2) * t424;
t59 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t444;
t372 = t293 * t59 + t473 * t81;
t60 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t443;
t371 = t293 * t60 + t473 * t82;
t148 = Icges(6,3) * t293 + t295 * t389;
t185 = (-Icges(6,5) * t292 - Icges(6,6) * t294) * t293;
t77 = qJD(2) * t148 + qJD(5) * t185;
t370 = t147 * t473 + t293 * t77;
t369 = t148 - t387;
t341 = t219 * t472 + t220 * t470 + t442 - t467;
t19 = qJD(2) * t406 + t245 * t87 - t246 * t88 + t341;
t368 = t19 * t407;
t52 = qJD(2) * t384 + t442;
t367 = t52 * t254;
t366 = -t126 + t422 - t80;
t364 = qJD(2) * t254;
t187 = (-Icges(6,1) * t292 - t502) * t293;
t342 = -t288 * t293 - t493 + t528;
t339 = -t586 * t302 + t589 * t304;
t338 = -(Icges(6,5) * t200 - Icges(6,6) * t201) * t246 - (Icges(6,5) * t202 - Icges(6,6) * t203) * t245 + t185 * t463;
t337 = qJD(2) * t543;
t336 = Icges(5,5) * t295 + (-Icges(5,1) * t303 + Icges(5,4) * t301) * t293;
t152 = Icges(6,5) * t293 + t295 * t399;
t334 = Icges(5,6) * t295 + (-Icges(5,4) * t303 + Icges(5,2) * t301) * t293;
t150 = Icges(6,6) * t293 + t295 * t394;
t330 = t293 * t338;
t327 = qJD(2) * t336;
t326 = qJD(2) * t334;
t320 = -qJD(2) * t237 - qJDD(2) * t254 + (-qJDD(2) * t307 - t487) * pkin(2);
t319 = -qJDD(4) * t295 + t145 * t472 + t146 * t470 + t219 * t459 + t220 * t458 + t293 * t461 + t379;
t317 = (Icges(6,1) * t202 - t504 - t84) * t245 + (Icges(6,1) * t200 - t505 - t83) * t246 - (-t149 + t187) * t463;
t313 = -t455 + (-t155 - t188) * qJD(2) + t424 * qJDD(2);
t310 = (-t369 * t463 - t544) * t293;
t269 = t295 * t465;
t267 = t295 * t466;
t248 = t278 * t304;
t247 = t278 * t302;
t199 = t304 * t364;
t198 = t302 * t364;
t144 = t304 * t404 + t291;
t143 = t302 * t404 - t468;
t140 = t336 * t304;
t139 = t336 * t302;
t138 = t334 * t304;
t137 = t334 * t302;
t130 = -t304 * t453 + t479;
t129 = -t302 * t453 + t480;
t125 = t151 * t304;
t124 = t151 * t302;
t123 = t149 * t304;
t122 = t149 * t302;
t119 = t304 * t327;
t118 = t302 * t327;
t117 = t304 * t326;
t116 = t302 * t326;
t109 = rSges(6,1) * t202 - rSges(6,2) * t203;
t108 = rSges(6,1) * t200 - rSges(6,2) * t201;
t99 = t304 * t337;
t98 = t302 * t337;
t79 = qJD(2) * t152 + qJD(5) * t187;
t78 = qJD(2) * t150 + qJD(5) * t186;
t76 = t304 * t320 + t290;
t75 = t302 * t320 - t457;
t74 = t304 * t373 + t476;
t73 = t302 * t373 + t430;
t64 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t443;
t63 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t444;
t62 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t443;
t61 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t444;
t54 = -qJD(2) * t380 + qJDD(2) * t381 + qJDD(1);
t53 = t293 * t387 - t498;
t45 = t304 * t313 + t447;
t44 = t302 * t313 + t405;
t39 = t158 * t246 + t304 * t374 + t463 * t87 + t476;
t38 = -t158 * t245 + t302 * t374 - t463 * t88 + t430;
t29 = t384 * qJDD(2) + (-t198 * t302 - t199 * t304) * qJD(2) + t379;
t28 = qJD(2) * t388 + t341;
t18 = (qJD(2) * t387 - t77) * t295 + (qJD(2) * t147 - t292 * t78 + t294 * t79 + (-t149 * t294 - t151 * t292) * qJD(5)) * t293;
t17 = t112 * t149 + t113 * t151 + t202 * t78 + t203 * t79 + t304 * t370;
t16 = t110 * t149 + t111 * t151 + t200 * t78 + t201 * t79 + t302 * t370;
t15 = qJD(2) * t559 + t388 * qJDD(2) + t319;
t12 = t245 * t31 + t246 * t30 - t463 * t53;
t11 = t112 * t84 + t113 * t86 + t202 * t62 + t203 * t64 + t304 * t371;
t10 = t112 * t83 + t113 * t85 + t202 * t61 + t203 * t63 + t304 * t372;
t9 = t110 * t84 + t111 * t86 + t200 * t62 + t201 * t64 + t302 * t371;
t8 = t110 * t83 + t111 * t85 + t200 * t61 + t201 * t63 + t302 * t372;
t7 = t245 * t27 + t246 * t26 - t463 * t47;
t6 = t24 * t246 + t245 * t25 - t46 * t463;
t5 = (qJD(2) * t411 - t60) * t295 + (qJD(2) * t82 - t292 * t62 + t294 * t64 + (-t292 * t86 - t294 * t84) * qJD(5)) * t293;
t4 = (qJD(2) * t412 - t59) * t295 + (qJD(2) * t81 - t292 * t61 + t294 * t63 + (-t292 * t85 - t294 * t83) * qJD(5)) * t293;
t3 = t319 + t162 * t87 - t163 * t88 + t245 * t65 - t246 * t66 + (t302 * t98 + t304 * t99) * qJD(2) + t406 * qJDD(2);
t2 = t10 * t246 + t11 * t245 + t162 * t27 + t163 * t26 - t17 * t463 + t238 * t47;
t1 = -t16 * t463 + t162 * t25 + t163 * t24 + t238 * t46 + t245 * t9 + t246 * t8;
t20 = [m(2) * qJDD(1) + (-m(2) - m(3) + t462) * g(3) + m(3) * t54 + m(4) * t29 + m(5) * t15 + m(6) * t3; (-t30 * t304 + t302 * t31) * t537 + (-t24 * t304 + t25 * t302) * t539 + (-t26 * t304 + t27 * t302) * t540 + (t302 * t9 - t304 * t8) * t533 + ((-t123 * t200 - t125 * t201) * t245 + (-t122 * t200 - t124 * t201) * t246 + (t513 + (-t150 * t200 - t152 * t201 + t515) * t295) * qJD(5) + (((t24 - t498) * qJD(5) + t414) * t295 + t310) * t302) * t534 + (-t10 * t304 + t11 * t302) * t535 + ((-t123 * t202 - t125 * t203) * t245 + (-t122 * t202 - t124 * t203) * t246 + (t512 + (-t150 * t202 - t152 * t203 + t514) * t295) * qJD(5) + (((t27 - t498) * qJD(5) + t414) * t295 + t310) * t304) * t536 + (((t123 * t292 - t125 * t294 + t82) * t245 + (t122 * t292 - t124 * t294 + t81) * t246 + t53 * qJD(5)) * t293 + ((t369 * t295 + (t150 * t292 - t152 * t294 - t147) * t293 + t410) * qJD(5) + t544) * t295) * t435 - t12 * t464 / 0.2e1 - ((t138 * t225 + t140 * t226) * t472 + t565 * qJD(2) * t297 + (-t225 * t137 - t226 * t139 - t542 * t304 + (t339 - t566) * t302 + t564) * t470) * t472 / 0.2e1 + (-(t137 * t223 + t139 * t224) * t470 + t566 * qJD(2) * t298 + (t223 * t138 + t224 * t140 + t339 * t302 + (-t542 - t565) * t304 + t564) * t472) * t437 + ((-t4 + t7) * t304 + (t5 + t6) * t302) * t436 + (t3 * t428 + (t14 * t378 + t39 * t366 + t3 * (t88 + t90)) * t304 + (t13 * t378 + t38 * t366 + t3 * (t87 + t89)) * t302 - t39 * (t159 * t246 + t269) - t38 * (-t159 * t245 + t267) - ((t38 * t88 - t39 * t87) * t293 + (t39 * (t158 * t302 + t129) + t38 * (-t158 * t304 - t130) + t368) * t295) * qJD(5) - g(1) * t479 - g(2) * t480 - g(3) * (t159 + t296 + t431) - t546 * (-t529 - t493 + (-t288 - t521) * t293) - (t38 * t302 + t39 * t304) * qJD(2) * (-t431 + t548 + t433) + ((t66 + t99) * t304 + (t65 + t98) * t302 - t245 * t129 + t246 * t130 - (-t427 + (t304 * t342 - t272) * t304 + (t302 * t342 - t271) * t302) * qJD(2) + t563) * t19) * m(6) + (t15 * t428 + (t15 * t101 + t377 * t74 + t424 * t45) * t304 + (t15 * t100 + t377 * t73 + t424 * t44) * t302 - t74 * t269 - t73 * t267 - g(1) * (-t304 * t529 + t272 + t477) - g(2) * (-t302 * t529 + t271 + t478) - t546 * t293 * (-pkin(3) - t522) + (-(t73 * t302 + t74 * t304) * qJD(2) + g(3)) * (t433 - t575) + (-(-t427 + (-t304 * t454 + t477) * t304 + (-t302 * t454 + t478) * t302) * qJD(2) + t559 + t563) * t28) * m(5) + (-(-t52 * t427 + (-t144 * t550 - t304 * t367) * t304 + (-t143 * t550 - t302 * t367) * t302) * qJD(2) - g(3) * t550 - t546 * t343 + t29 * t484 + t52 * t481 + (t144 * t421 + t29 * t182 - t52 * t199 + t343 * t76) * t304 + (t143 * t421 + t29 * t181 - t52 * t198 + t343 * t75) * t302) * m(4) + (g(1) * t248 + g(2) * t247 - g(3) * t279 + t54 * t381 + (-t380 - (-t247 * t302 - t248 * t304) * qJD(2)) * (qJD(2) * t381 + qJD(1)) + (qJDD(2) * t278 + t279 * t309) * t580) * m(3) + (t2 + ((t117 * t225 + t119 * t226 + t578 * t302) * t302 + (-t116 * t225 - t118 * t226 - t302 * t576) * t304) * t570 + ((-t225 * t94 - t226 * t96 - t496 * t92) * t304 + (t225 * t95 + t226 * t97 + t579 * t302 + t496 * t93 - t582) * t302) * t569) * t302 / 0.2e1 - (t1 + ((-t116 * t223 - t118 * t224 + t576 * t304) * t304 + (t117 * t223 + t119 * t224 - t304 * t578) * t302) * t570 + ((-t223 * t94 - t224 * t96 - t497 * t92 + t582) * t304 + (t223 * t95 + t224 * t97 - t304 * t579 + t497 * t93) * t302) * t569) * t304 / 0.2e1; t462 * (g(1) * t302 - g(2) * t304) + m(4) * (t302 * t76 - t304 * t75) + m(5) * (t302 * t45 - t304 * t44) + m(6) * t416; -t541 * g(3) * t295 + 0.2e1 * (t19 * t538 + t3 * t532) * m(6) + 0.2e1 * (t15 * t532 + t28 * t538) * m(5) + (t541 * t546 + m(5) * (qJD(2) * t28 + t302 * t44 + t304 * t45) + m(6) * (qJD(2) * t19 + t13 * t302 + t14 * t304)) * t293; t295 * t7 * t437 + t2 * t496 / 0.2e1 + (t293 * t413 - t295 * t47) * t540 + (-t17 * t295 + (t10 * t302 + t11 * t304) * t293 + (t295 * t413 + t512) * qJD(2)) * t535 + t6 * t444 / 0.2e1 + t1 * t497 / 0.2e1 + (t293 * t415 - t295 * t46) * t539 + (-t16 * t295 + (t302 * t8 + t304 * t9) * t293 + (t295 * t415 + t513) * qJD(2)) * t533 + t12 * t474 / 0.2e1 + (t162 * t31 + t163 * t30 - t18 * t463 + t238 * t53 + t245 * t5 + t246 * t4) * t532 + (t293 * t410 - t295 * t53) * t537 + (-t18 * t295 + (t302 * t4 + t304 * t5) * t293 + (t53 * t293 + t295 * t410) * qJD(2)) * t436 + (t202 * t318 + t203 * t317 - t304 * t330) * t536 + (t200 * t318 + t201 * t317 - t302 * t330) * t534 + (t338 * t295 + (-t292 * t318 + t317 * t294) * t293) * t435 + ((-t13 * t88 + t14 * t87 - t38 * t66 + t39 * t65 + (t368 + (t302 * t39 - t304 * t38) * t158) * qJD(2)) * t295 + (t39 * (-qJD(2) * t87 + t302 * t80) + t38 * (qJD(2) * t88 - t304 * t80) + t3 * t407 + t19 * (-t302 * t66 + t304 * t65) + t416 * t158) * t293 - t39 * (t108 * t463 + t189 * t246) - t38 * (-t109 * t463 - t189 * t245) - t19 * (t108 * t245 - t109 * t246) - g(1) * t109 - g(2) * t108 - g(3) * t189) * m(6);];
tau = t20;
