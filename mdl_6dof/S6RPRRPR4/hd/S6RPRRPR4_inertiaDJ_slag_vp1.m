% Calculate time derivative of joint inertia matrix for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:09:04
% EndTime: 2019-03-09 05:09:29
% DurationCPUTime: 14.40s
% Computational Cost: add. (40984->911), mult. (33115->1239), div. (0->0), fcn. (30817->12), ass. (0->447)
t348 = pkin(10) + qJ(3);
t340 = qJ(4) + t348;
t331 = cos(t340);
t354 = cos(pkin(11));
t358 = sin(qJ(1));
t503 = t358 * t354;
t352 = sin(pkin(11));
t359 = cos(qJ(1));
t509 = t352 * t359;
t291 = t331 * t503 - t509;
t504 = t358 * t352;
t508 = t354 * t359;
t384 = t331 * t504 + t508;
t330 = sin(t340);
t519 = t330 * t358;
t198 = Icges(6,4) * t291 - Icges(6,2) * t384 + Icges(6,6) * t519;
t200 = Icges(6,1) * t291 - Icges(6,4) * t384 + Icges(6,5) * t519;
t413 = Icges(5,5) * t331 - Icges(5,6) * t330;
t256 = -Icges(5,3) * t359 + t358 * t413;
t292 = -t331 * t509 + t503;
t293 = t331 * t508 + t504;
t196 = Icges(6,5) * t291 - Icges(6,6) * t384 + Icges(6,3) * t519;
t518 = t330 * t359;
t473 = t196 * t518;
t545 = Icges(5,4) * t331;
t417 = -Icges(5,2) * t330 + t545;
t258 = -Icges(5,6) * t359 + t358 * t417;
t546 = Icges(5,4) * t330;
t422 = Icges(5,1) * t331 - t546;
t260 = -Icges(5,5) * t359 + t358 * t422;
t399 = t258 * t330 - t260 * t331;
t581 = t359 * t399;
t590 = -t292 * t198 - t293 * t200 - t358 * t256 - t473 + t581;
t199 = Icges(6,4) * t293 + Icges(6,2) * t292 + Icges(6,6) * t518;
t201 = Icges(6,1) * t293 + Icges(6,4) * t292 + Icges(6,5) * t518;
t257 = Icges(5,3) * t358 + t359 * t413;
t259 = Icges(5,6) * t358 + t359 * t417;
t261 = Icges(5,5) * t358 + t359 * t422;
t398 = t259 * t330 - t261 * t331;
t197 = Icges(6,5) * t293 + Icges(6,6) * t292 + Icges(6,3) * t518;
t471 = t197 * t518;
t589 = t292 * t199 + t293 * t201 + t358 * t257 - t359 * t398 + t471;
t347 = pkin(11) + qJ(6);
t338 = cos(t347);
t349 = qJD(3) + qJD(4);
t336 = sin(t347);
t543 = Icges(7,4) * t338;
t415 = -Icges(7,2) * t336 + t543;
t517 = t331 * t349;
t544 = Icges(7,4) * t336;
t146 = t415 * t517 + (Icges(7,6) * t349 + (-Icges(7,2) * t338 - t544) * qJD(6)) * t330;
t241 = -Icges(7,6) * t331 + t330 * t415;
t420 = Icges(7,1) * t338 - t544;
t242 = -Icges(7,5) * t331 + t330 * t420;
t588 = -t146 * t336 + (-t241 * t338 - t242 * t336) * qJD(6);
t332 = pkin(5) * t354 + pkin(4);
t557 = pkin(4) - t332;
t587 = t330 * t557;
t356 = -pkin(9) - qJ(5);
t502 = qJ(5) + t356;
t586 = t331 * t502;
t585 = t331 * t557;
t337 = sin(t348);
t339 = cos(t348);
t547 = Icges(4,4) * t339;
t419 = -Icges(4,2) * t337 + t547;
t272 = Icges(4,6) * t358 + t359 * t419;
t548 = Icges(4,4) * t337;
t424 = Icges(4,1) * t339 - t548;
t274 = Icges(4,5) * t358 + t359 * t424;
t396 = t272 * t337 - t274 * t339;
t584 = t358 * t396;
t583 = t358 * t398;
t271 = -Icges(4,6) * t359 + t358 * t419;
t273 = -Icges(4,5) * t359 + t358 * t424;
t397 = t271 * t337 - t273 * t339;
t582 = t359 * t397;
t411 = Icges(7,5) * t338 - Icges(7,6) * t336;
t145 = t411 * t517 + (Icges(7,3) * t349 + (-Icges(7,5) * t336 - Icges(7,6) * t338) * qJD(6)) * t330;
t532 = t241 * t336;
t580 = -t349 * t532 - t145;
t429 = rSges(7,1) * t338 - rSges(7,2) * t336;
t150 = t429 * t517 + (rSges(7,3) * t349 + (-rSges(7,1) * t336 - rSges(7,2) * t338) * qJD(6)) * t330;
t367 = -t330 * t502 - t585;
t579 = -t367 * t349 - t150;
t506 = t358 * t336;
t513 = t338 * t359;
t281 = -t331 * t506 - t513;
t505 = t358 * t338;
t514 = t336 * t359;
t282 = t331 * t505 - t514;
t430 = -t282 * rSges(7,1) - t281 * rSges(7,2);
t177 = rSges(7,3) * t519 - t430;
t283 = -t331 * t514 + t505;
t284 = t331 * t513 + t506;
t178 = t284 * rSges(7,1) + t283 * rSges(7,2) + rSges(7,3) * t518;
t578 = -t358 * t177 - t359 * t178;
t233 = t586 - t587;
t243 = -rSges(7,3) * t331 + t330 * t429;
t577 = -t233 - t243;
t357 = -pkin(7) - qJ(2);
t355 = cos(pkin(10));
t333 = t355 * pkin(2) + pkin(1);
t554 = rSges(4,2) * t337;
t556 = rSges(4,1) * t339;
t436 = -t554 + t556;
t388 = -t333 - t436;
t229 = (rSges(4,3) - t357) * t359 + t388 * t358;
t345 = t358 * rSges(4,3);
t490 = t359 * t556 + t345;
t230 = -t358 * t357 + (t333 - t554) * t359 + t490;
t576 = t229 * t359 + t230 * t358;
t575 = qJD(1) * t256;
t414 = Icges(4,5) * t339 - Icges(4,6) * t337;
t269 = -Icges(4,3) * t359 + t358 * t414;
t439 = qJD(1) * t331 - qJD(6);
t510 = t349 * t359;
t468 = t330 * t510;
t574 = t358 * t439 + t468;
t235 = qJD(1) * t384 + t352 * t468;
t236 = -qJD(1) * t291 - t354 * t468;
t485 = qJD(1) * t358;
t453 = t330 * t485;
t467 = t331 * t510;
t371 = -t453 + t467;
t122 = Icges(6,5) * t236 + Icges(6,6) * t235 + Icges(6,3) * t371;
t511 = t349 * t358;
t469 = t330 * t511;
t237 = qJD(1) * t292 + t352 * t469;
t238 = qJD(1) * t293 - t354 * t469;
t484 = qJD(1) * t359;
t372 = t330 * t484 + t331 * t511;
t123 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t372;
t573 = -(t123 * t330 + t196 * t517) * t359 + (t122 * t330 + t197 * t517) * t358;
t572 = t359 * t439 - t469;
t299 = Icges(5,2) * t331 + t546;
t300 = Icges(5,1) * t330 + t545;
t395 = t299 * t330 - t300 * t331;
t570 = qJD(1) * t395 + t413 * t349;
t390 = rSges(3,1) * t355 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t551 = rSges(3,3) + qJ(2);
t255 = t358 * t551 + t359 * t390;
t346 = -pkin(8) + t357;
t489 = t346 - t357;
t314 = pkin(3) * t339 + t333;
t494 = t314 - t333;
t231 = t358 * t494 + t359 * t489;
t569 = 2 * m(4);
t568 = 2 * m(5);
t567 = 2 * m(6);
t566 = 2 * m(7);
t350 = t358 ^ 2;
t351 = t359 ^ 2;
t565 = m(6) / 0.2e1;
t564 = m(7) / 0.2e1;
t563 = t358 / 0.2e1;
t562 = -t359 / 0.2e1;
t309 = rSges(4,1) * t337 + rSges(4,2) * t339;
t561 = m(4) * t309;
t302 = rSges(5,1) * t330 + rSges(5,2) * t331;
t560 = m(5) * t302;
t559 = pkin(3) * t337;
t558 = pkin(4) * t331;
t555 = rSges(5,1) * t331;
t168 = Icges(7,5) * t282 + Icges(7,6) * t281 + Icges(7,3) * t519;
t170 = Icges(7,4) * t282 + Icges(7,2) * t281 + Icges(7,6) * t519;
t172 = Icges(7,1) * t282 + Icges(7,4) * t281 + Icges(7,5) * t519;
t407 = -t170 * t336 + t172 * t338;
t440 = -qJD(6) * t331 + qJD(1);
t393 = t358 * t440;
t159 = -t336 * t572 + t338 * t393;
t160 = t336 * t393 + t338 * t572;
t90 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t372;
t92 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t372;
t94 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t372;
t23 = (t349 * t407 - t90) * t331 + (t168 * t349 - t336 * t92 + t338 * t94 + (-t170 * t338 - t172 * t336) * qJD(6)) * t330;
t553 = t23 * t359;
t169 = Icges(7,5) * t284 + Icges(7,6) * t283 + Icges(7,3) * t518;
t171 = Icges(7,4) * t284 + Icges(7,2) * t283 + Icges(7,6) * t518;
t173 = Icges(7,1) * t284 + Icges(7,4) * t283 + Icges(7,5) * t518;
t406 = -t171 * t336 + t173 * t338;
t392 = t359 * t440;
t157 = t336 * t574 + t338 * t392;
t158 = t336 * t392 - t338 * t574;
t89 = Icges(7,5) * t158 + Icges(7,6) * t157 + Icges(7,3) * t371;
t91 = Icges(7,4) * t158 + Icges(7,2) * t157 + Icges(7,6) * t371;
t93 = Icges(7,1) * t158 + Icges(7,4) * t157 + Icges(7,5) * t371;
t24 = (t349 * t406 - t89) * t331 + (t169 * t349 - t336 * t91 + t338 * t93 + (-t171 * t338 - t173 * t336) * qJD(6)) * t330;
t552 = t24 * t358;
t344 = t358 * rSges(5,3);
t550 = -rSges(6,3) - qJ(5);
t549 = -rSges(7,3) + t356;
t536 = qJ(5) * t330;
t529 = t242 * t338;
t447 = -t302 - t559;
t251 = t447 * t359;
t528 = t251 * t359;
t527 = t271 * t339;
t526 = t272 * t339;
t525 = t273 * t337;
t524 = t274 * t337;
t523 = t299 * t349;
t522 = t300 * t349;
t521 = t330 * t332;
t520 = t330 * t349;
t516 = t331 * t356;
t515 = t331 * t359;
t512 = t346 * t359;
t435 = -rSges(5,2) * t330 + t555;
t280 = t435 * t349;
t507 = t358 * t280;
t205 = t293 * rSges(6,1) + t292 * rSges(6,2) + rSges(6,3) * t518;
t324 = pkin(4) * t515;
t286 = qJ(5) * t518 + t324;
t501 = -t205 - t286;
t303 = t359 * t314;
t232 = -t333 * t359 - t358 * t489 + t303;
t500 = t358 * t231 + t359 * t232;
t432 = rSges(6,1) * t354 - rSges(6,2) * t352;
t224 = (rSges(6,3) * t330 + t331 * t432) * t349;
t244 = qJ(5) * t520 + (pkin(4) * t349 - qJD(5)) * t331;
t499 = -t224 - t244;
t249 = -rSges(6,3) * t331 + t330 * t432;
t301 = pkin(4) * t330 - qJ(5) * t331;
t287 = t301 * t485;
t498 = t249 * t485 + t287;
t497 = -t249 - t301;
t263 = -rSges(5,3) * t359 + t358 * t435;
t264 = rSges(5,1) * t515 - rSges(5,2) * t518 + t344;
t181 = t358 * t263 + t359 * t264;
t285 = (t536 + t558) * t358;
t496 = t358 * t285 + t359 * t286;
t478 = pkin(5) * t509;
t495 = qJD(1) * t478 + t356 * t453;
t493 = rSges(5,2) * t453 + rSges(5,3) * t484;
t480 = qJD(5) * t330;
t316 = t359 * t480;
t341 = qJD(2) * t358;
t492 = t316 + t341;
t483 = qJD(3) * t337;
t476 = pkin(3) * t483;
t491 = -t346 * t485 - t358 * t476;
t488 = t350 + t351;
t487 = qJD(1) * t257;
t270 = Icges(4,3) * t358 + t359 * t414;
t486 = qJD(1) * t270;
t482 = qJD(3) * t339;
t481 = qJD(3) * t358;
t329 = pkin(5) * t504;
t477 = t359 * t554;
t475 = pkin(3) * t482;
t474 = t196 * t519;
t472 = t197 * t519;
t240 = -Icges(7,3) * t331 + t330 * t411;
t100 = t240 * t519 + t241 * t281 + t242 * t282;
t74 = -t168 * t331 + t330 * t407;
t466 = t74 / 0.2e1 + t100 / 0.2e1;
t101 = t240 * t518 + t283 * t241 + t284 * t242;
t75 = -t169 * t331 + t330 * t406;
t465 = t75 / 0.2e1 + t101 / 0.2e1;
t147 = t420 * t517 + (Icges(7,5) * t349 + (-Icges(7,1) * t336 - t543) * qJD(6)) * t330;
t464 = t330 * t338 * t147 + t240 * t520 + t517 * t529;
t463 = -t244 + t579;
t462 = t158 * rSges(7,1) + t157 * rSges(7,2) + rSges(7,3) * t467;
t308 = qJ(5) * t467;
t312 = pkin(4) * t469;
t373 = -t331 * t485 - t468;
t461 = t358 * (qJ(5) * t372 + qJD(1) * t324 + t358 * t480 - t312) + t359 * (pkin(4) * t373 - qJ(5) * t453 + t308 + t316) + t285 * t484;
t383 = t302 * t349;
t460 = t358 * (-t358 * t383 + (t359 * t435 + t344) * qJD(1)) + t359 * (rSges(5,1) * t373 - rSges(5,2) * t467 + t493) + t263 * t484;
t389 = t332 * t515 - t356 * t518 + t329;
t194 = t389 - t286;
t459 = -t178 - t194 - t286;
t328 = t357 * t485;
t458 = t358 * (t484 * t494 + t328 + t491) + t359 * (-qJD(1) * t231 - t359 * t476) + t231 * t484;
t234 = t243 * t485;
t457 = t233 * t485 + t234 + t287;
t456 = t236 * rSges(6,1) + t235 * rSges(6,2) + rSges(6,3) * t467;
t455 = -t301 + t577;
t342 = qJD(2) * t359;
t454 = t342 - t491;
t452 = t337 * t485;
t451 = t517 / 0.2e1;
t450 = t485 / 0.2e1;
t449 = t484 / 0.2e1;
t448 = -t301 - t559;
t203 = t497 * t359;
t184 = -qJD(1) * t258 - t359 * t523;
t446 = t261 * t349 + t184;
t185 = qJD(1) * t259 - t358 * t523;
t445 = t260 * t349 + t185;
t186 = -qJD(1) * t260 - t359 * t522;
t444 = -t259 * t349 + t186;
t187 = qJD(1) * t261 - t358 * t522;
t443 = t258 * t349 - t187;
t442 = -t358 * t346 + t303;
t441 = -t331 * t332 - t314;
t433 = -t291 * rSges(6,1) + rSges(6,2) * t384;
t204 = rSges(6,3) * t519 - t433;
t99 = t358 * t204 + t359 * t205 + t496;
t438 = -t249 + t448;
t437 = -t244 - t475;
t133 = t455 * t359;
t434 = t238 * rSges(6,1) + t237 * rSges(6,2);
t431 = t160 * rSges(7,1) + t159 * rSges(7,2);
t64 = t168 * t519 + t170 * t281 + t172 * t282;
t65 = t169 * t519 + t171 * t281 + t173 * t282;
t44 = t65 * t358 - t359 * t64;
t428 = t358 * t64 + t359 * t65;
t66 = t168 * t518 + t283 * t170 + t284 * t172;
t67 = t169 * t518 + t283 * t171 + t284 * t173;
t45 = t67 * t358 - t359 * t66;
t427 = t358 * t66 + t359 * t67;
t426 = t75 * t358 - t359 * t74;
t425 = t358 * t74 + t359 * t75;
t423 = Icges(4,1) * t337 + t547;
t421 = Icges(6,1) * t354 - Icges(6,4) * t352;
t418 = Icges(4,2) * t339 + t548;
t416 = Icges(6,4) * t354 - Icges(6,2) * t352;
t298 = Icges(5,5) * t330 + Icges(5,6) * t331;
t412 = Icges(6,5) * t354 - Icges(6,6) * t352;
t370 = t330 * t549 + t441;
t114 = (pkin(5) * t352 - t346) * t359 + t370 * t358 + t430;
t115 = t389 + t442 + t178;
t410 = t114 * t359 + t115 * t358;
t116 = t177 * t331 + t243 * t519;
t117 = -t331 * t178 - t243 * t518;
t409 = t116 * t359 + t117 * t358;
t374 = t330 * t550 - t314 - t558;
t363 = t358 * t374 - t512;
t136 = t363 + t433;
t137 = t442 - t501;
t408 = t136 * t359 + t137 * t358;
t405 = t177 * t359 - t178 * t358;
t404 = -t198 * t352 + t200 * t354;
t403 = -t199 * t352 + t201 * t354;
t387 = -t314 - t435;
t214 = (rSges(5,3) - t346) * t359 + t387 * t358;
t215 = t264 + t442;
t400 = t214 * t359 + t215 * t358;
t394 = t448 + t577;
t391 = -t224 + t437;
t167 = t438 * t359;
t386 = t358 * (rSges(6,3) * t372 + t434) + t359 * (-rSges(6,3) * t453 + t456) + t204 * t484 + t461;
t193 = t358 * t367 - t478;
t56 = t358 * t193 + t359 * t194 + t496 - t578;
t382 = t437 + t579;
t381 = qJD(3) * t309;
t378 = t349 * t298;
t377 = qJD(3) * t423;
t376 = qJD(3) * t418;
t375 = qJD(3) * (-Icges(4,5) * t337 - Icges(4,6) * t339);
t121 = t394 * t359;
t124 = Icges(6,4) * t236 + Icges(6,2) * t235 + Icges(6,6) * t371;
t125 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t372;
t126 = Icges(6,1) * t236 + Icges(6,4) * t235 + Icges(6,5) * t371;
t127 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t372;
t128 = -t256 * t359 - t358 * t399;
t129 = -t257 * t359 - t583;
t182 = -t359 * t378 - t575;
t183 = -t358 * t378 + t487;
t76 = -t198 * t384 + t200 * t291 + t474;
t77 = -t199 * t384 + t201 * t291 + t472;
t17 = t157 * t170 + t158 * t172 + t168 * t371 + t283 * t92 + t284 * t94 + t518 * t90;
t18 = t157 * t171 + t158 * t173 + t169 * t371 + t283 * t91 + t284 * t93 + t518 * t89;
t8 = qJD(1) * t427 - t17 * t359 + t18 * t358;
t369 = t44 * t485 + t45 * t484 + ((-t128 - t76) * t485 + t590 * t484) * t359 + (t8 + (t292 * t124 + t293 * t126 + t358 * t182 + t235 * t199 + t236 * t201 + (-t472 + t583 - t590) * qJD(1)) * t358 + (t77 + t129) * t485 + t589 * t484 + (-t292 * t125 - t293 * t127 - t235 * t198 - t236 * t200 + t573 + (-t446 * t330 + t444 * t331 - t183) * t358 + (t185 * t330 - t187 * t331 + t258 * t517 + t260 * t520 - t575) * t359 + (t474 + (t257 - t399) * t358 + t589) * qJD(1)) * t359) * t358;
t95 = -rSges(7,3) * t453 + t462;
t96 = rSges(7,3) * t372 + t431;
t366 = t461 + (t177 + t193) * t484 + (t95 - t308 + (-t516 + t587) * t510 + (t536 + t585) * t485 + t495) * t359 + (t96 + t312 + (-t521 - t586) * t511 + (t359 * t367 + t329) * qJD(1)) * t358;
t28 = -t100 * t331 + t330 * t428;
t29 = -t101 * t331 + t330 * t427;
t33 = t145 * t518 + t283 * t146 + t284 * t147 + t157 * t241 + t158 * t242 + t240 * t371;
t3 = (t349 * t427 - t33) * t331 + (-qJD(1) * t45 + t349 * t101 + t17 * t358 + t18 * t359) * t330;
t19 = t159 * t170 + t160 * t172 + t168 * t372 + t281 * t92 + t282 * t94 + t519 * t90;
t20 = t159 * t171 + t160 * t173 + t169 * t372 + t281 * t91 + t282 * t93 + t519 * t89;
t34 = t145 * t519 + t281 * t146 + t282 * t147 + t159 * t241 + t160 * t242 + t240 * t372;
t4 = (t349 * t428 - t34) * t331 + (-qJD(1) * t44 + t100 * t349 + t19 * t358 + t20 * t359) * t330;
t9 = qJD(1) * t428 - t19 * t359 + t20 * t358;
t365 = t3 * t563 + t4 * t562 - t331 * (qJD(1) * t425 + t552 - t553) / 0.2e1 + t28 * t450 + t29 * t449 + t426 * t520 / 0.2e1 + t9 * t519 / 0.2e1 + t8 * t518 / 0.2e1 + (t359 * t451 - t453 / 0.2e1) * t45 + (t330 * t449 + t358 * t451) * t44;
t364 = rSges(4,2) * t452 + rSges(4,3) * t484 - t359 * t381;
t254 = -t358 * t390 + t359 * t551;
t12 = (t384 * t125 - t291 * t127 - t237 * t198 - t238 * t200 + (t77 - t473) * qJD(1)) * t359 + (-t384 * t124 + t291 * t126 + t237 * t199 + t238 * t201 + (t76 + t471) * qJD(1) + t573) * t358;
t16 = (t183 * t359 + (t129 + t581) * qJD(1)) * t359 + (t128 * qJD(1) + (-t184 * t330 + t186 * t331 - t259 * t517 - t261 * t520 + t487) * t358 + (-t182 + t443 * t331 + t445 * t330 + (-t256 - t398) * qJD(1)) * t359) * t358;
t362 = (-t9 - t16 - t12) * t359 + t369;
t276 = t417 * t349;
t277 = t422 * t349;
t361 = qJD(1) * t298 + (t277 - t523) * t331 + (-t276 - t522) * t330;
t221 = (Icges(6,3) * t330 + t331 * t412) * t349;
t222 = (Icges(6,6) * t330 + t331 * t416) * t349;
t223 = (Icges(6,5) * t330 + t331 * t421) * t349;
t245 = -Icges(6,3) * t331 + t330 * t412;
t246 = -Icges(6,6) * t331 + t330 * t416;
t247 = -Icges(6,5) * t331 + t330 * t421;
t360 = -t553 / 0.2e1 + t552 / 0.2e1 + (t221 * t518 + t292 * t222 + t293 * t223 + t235 * t246 + t236 * t247 + t245 * t371 + t358 * t570 + t361 * t359 + t33 + (t349 * t403 - t122 + t446) * t331 + (-t124 * t352 + t126 * t354 + t197 * t349 + t444) * t330) * t563 + (t221 * t519 - t384 * t222 + t291 * t223 + t237 * t246 + t238 * t247 + t245 * t372 + t361 * t358 - t359 * t570 + t34 + (t349 * t404 - t123 + t445) * t331 + (-t125 * t352 + t127 * t354 + t196 * t349 - t443) * t330) * t562 + (t245 * t519 - t246 * t384 + t247 * t291 - t298 * t359 - t358 * t395 + t100 + t74 + (-t196 + t258) * t331 + (t260 + t404) * t330) * t450 + (t245 * t518 + t292 * t246 + t293 * t247 + t358 * t298 - t359 * t395 + t101 + t75 + (-t197 + t259) * t331 + (t261 + t403) * t330) * t449;
t320 = pkin(3) * t452;
t297 = t436 * qJD(3);
t279 = -t477 + t490;
t278 = -rSges(4,3) * t359 + t358 * t436;
t250 = t447 * t358;
t226 = -qJD(1) * t255 + t342;
t225 = qJD(1) * t254 + t341;
t208 = t358 * t375 + t486;
t207 = -qJD(1) * t269 + t359 * t375;
t202 = t497 * t358;
t166 = t438 * t358;
t156 = -t302 * t484 - t507 + (-t337 * t484 - t339 * t481) * pkin(3);
t155 = t302 * t485 + t320 + (-t280 - t475) * t359;
t152 = t328 + t342 + t309 * t481 + (t359 * t388 - t345) * qJD(1);
t151 = t341 + (-t357 * t359 + (-t333 - t556) * t358) * qJD(1) + t364;
t141 = t358 * t270 - t359 * t396;
t140 = t358 * t269 - t582;
t139 = -t270 * t359 - t584;
t138 = -t269 * t359 - t358 * t397;
t135 = t302 * t511 + (t359 * t387 - t344) * qJD(1) + t454;
t134 = t341 + (-t314 - t555) * t485 + (-qJD(1) * t346 - t383 - t476) * t359 + t493;
t132 = t455 * t358;
t120 = t394 * t358;
t111 = -t240 * t331 + (t529 - t532) * t330;
t110 = qJD(1) * t203 + t358 * t499;
t109 = t359 * t499 + t498;
t108 = t405 * t330;
t107 = t181 + t500;
t106 = t111 * t520;
t103 = qJD(1) * t167 + t358 * t391;
t102 = t359 * t391 + t320 + t498;
t86 = -t264 * t485 + t460;
t83 = t312 + (t517 * t550 - t480) * t358 + t374 * t484 - t434 + t454;
t82 = t308 + (-pkin(4) * t520 - t476) * t359 + t363 * qJD(1) + t456 + t492;
t63 = t99 + t500;
t62 = qJD(1) * t133 + t358 * t463;
t61 = t359 * t463 + t457;
t60 = qJD(1) * t121 + t358 * t382;
t59 = t359 * t382 + t320 + t457;
t58 = (-t480 + (t331 * t549 + t521) * t349) * t358 + (t359 * t370 - t329) * qJD(1) - t431 + t454;
t57 = (-t476 + (-t516 - t521) * t349) * t359 + (-t512 + (-rSges(7,3) * t330 + t441) * t358) * qJD(1) + t462 + t492 + t495;
t55 = t56 + t500;
t54 = (-t232 - t264) * t485 + t458 + t460;
t51 = (t243 * t511 + t96) * t331 + (t358 * t150 - t177 * t349 + t243 * t484) * t330;
t50 = (-t243 * t510 - t95) * t331 + (-t150 * t359 + t178 * t349 + t234) * t330;
t46 = t330 * t588 + t580 * t331 + t464;
t37 = t485 * t501 + t386;
t30 = t405 * t517 + (qJD(1) * t578 - t358 * t95 + t359 * t96) * t330;
t25 = (-t232 + t501) * t485 + t386 + t458;
t14 = t459 * t485 + t366;
t13 = (-t232 + t459) * t485 + t366 + t458;
t1 = [(t114 * t58 + t115 * t57) * t566 + (t136 * t83 + t137 * t82) * t567 + (t134 * t215 + t135 * t214) * t568 + (t151 * t230 + t152 * t229) * t569 + 0.2e1 * m(3) * (t225 * t255 + t226 * t254) + t464 + (t245 - t299) * t520 + (-t418 + t424) * t483 + (t423 + t419) * t482 + (-t246 * t352 + t247 * t354 + t300) * t517 + (t276 - t221 + t580) * t331 + (-t222 * t352 + t223 * t354 + t277 + t588) * t330; m(7) * (qJD(1) * t410 + t358 * t58 - t359 * t57) + m(6) * (qJD(1) * t408 + t358 * t83 - t359 * t82) + m(5) * (qJD(1) * t400 - t134 * t359 + t358 * t135) + m(4) * (qJD(1) * t576 - t151 * t359 + t358 * t152) + m(3) * (-t225 * t359 + t358 * t226 + (t254 * t359 + t255 * t358) * qJD(1)); 0; t360 + (t350 / 0.2e1 + t351 / 0.2e1) * t414 * qJD(3) + (-qJD(3) * t397 + (qJD(1) * t272 - t358 * t376) * t339 + (qJD(1) * t274 - t358 * t377) * t337) * t562 + (-qJD(3) * t396 + (-qJD(1) * t271 - t359 * t376) * t339 + (-qJD(1) * t273 - t359 * t377) * t337) * t563 + m(7) * (t114 * t59 + t115 * t60 + t120 * t57 + t121 * t58) + m(6) * (t102 * t136 + t103 * t137 + t166 * t82 + t167 * t83) + m(5) * (t134 * t250 + t135 * t251 + t155 * t214 + t156 * t215) + ((-t230 * t561 + t526 / 0.2e1 + t524 / 0.2e1) * t359 + (t229 * t561 + t527 / 0.2e1 + t525 / 0.2e1) * t358) * qJD(1) + m(4) * ((-t151 * t358 - t152 * t359) * t309 - t576 * t297); m(5) * (t155 * t358 - t156 * t359 + (t250 * t358 + t528) * qJD(1)) + m(6) * (t102 * t358 - t103 * t359 + (t166 * t358 + t167 * t359) * qJD(1)) + m(7) * (t59 * t358 - t359 * t60 + (t120 * t358 + t121 * t359) * qJD(1)); t369 + t358 * ((t358 * t207 + (t140 + t584) * qJD(1)) * t358 + (t141 * qJD(1) + (t271 * t482 + t273 * t483) * t359 + (-t208 + (-t524 - t526) * qJD(3) + (t270 - t397) * qJD(1)) * t358) * t359) - t359 * t9 - t359 * t16 - t359 * t12 - t359 * ((t208 * t359 + (t139 + t582) * qJD(1)) * t359 + (t138 * qJD(1) + (-t272 * t482 - t274 * t483 + t486) * t358 + (-t207 + (t525 + t527) * qJD(3) - t396 * qJD(1)) * t359) * t358) + (t120 * t60 + t121 * t59 + t13 * t55) * t566 + (t102 * t167 + t103 * t166 + t25 * t63) * t567 + (t107 * t54 + t155 * t251 + t156 * t250) * t568 + (-t138 * t359 + t139 * t358) * t485 + (-t140 * t359 + t141 * t358) * t484 + ((t358 * t278 + t279 * t359) * ((qJD(1) * t278 + t364) * t359 + (-t358 * t381 + (-t279 - t477 + t345) * qJD(1)) * t358) + t488 * t309 * t297) * t569; t360 + (-t134 * t358 - t135 * t359 + (t214 * t358 - t215 * t359) * qJD(1)) * t560 - m(5) * t400 * t280 + m(7) * (t114 * t61 + t115 * t62 + t132 * t57 + t133 * t58) + m(6) * (t109 * t136 + t110 * t137 + t202 * t82 + t203 * t83); m(6) * (t109 * t358 - t110 * t359 + (t202 * t358 + t203 * t359) * qJD(1)) + m(7) * (t61 * t358 - t359 * t62 + (t132 * t358 + t133 * t359) * qJD(1)); m(5) * (t86 * t107 + t181 * t54 - t250 * t507 - t280 * t528) + m(7) * (t120 * t62 + t121 * t61 + t13 * t56 + t132 * t60 + t133 * t59 + t14 * t55) + m(6) * (t102 * t203 + t103 * t202 + t109 * t167 + t110 * t166 + t25 * t99 + t37 * t63) + t362 + (-t155 * t359 - t156 * t358 + (-t250 * t359 + t251 * t358) * qJD(1)) * t560; (t132 * t62 + t133 * t61 + t14 * t56) * t566 + (t109 * t203 + t110 * t202 + t37 * t99) * t567 + (t280 * t302 * t488 + t181 * t86) * t568 + t362; 0.2e1 * (t408 * t565 + t410 * t564) * t517 + 0.2e1 * ((-t114 * t485 + t115 * t484 + t358 * t57 + t359 * t58) * t564 + (-t136 * t485 + t137 * t484 + t358 * t82 + t359 * t83) * t565) * t330; 0; 0.2e1 * ((t120 * t511 + t121 * t510 - t13) * t564 + (t166 * t511 + t167 * t510 - t25) * t565) * t331 + 0.2e1 * ((t120 * t484 - t121 * t485 + t349 * t55 + t358 * t60 + t359 * t59) * t564 + (t102 * t359 + t103 * t358 + t166 * t484 - t167 * t485 + t349 * t63) * t565) * t330; 0.2e1 * ((t132 * t511 + t133 * t510 - t14) * t564 + (t202 * t511 + t203 * t510 - t37) * t565) * t331 + 0.2e1 * ((t132 * t484 - t133 * t485 + t349 * t56 + t358 * t62 + t359 * t61) * t564 + (t109 * t359 + t110 * t358 + t202 * t484 - t203 * t485 + t349 * t99) * t565) * t330; 0.4e1 * (t565 + t564) * (-0.1e1 + t488) * t330 * t517; m(7) * (t114 * t51 + t115 * t50 + t116 * t58 + t117 * t57) + t106 + (-t46 + (t358 * t466 + t359 * t465) * t349) * t331 + ((t24 / 0.2e1 + t33 / 0.2e1) * t359 + (t23 / 0.2e1 + t34 / 0.2e1) * t358 + (-t358 * t465 + t359 * t466) * qJD(1)) * t330; m(7) * (qJD(1) * t409 + t51 * t358 - t359 * t50); t365 + m(7) * (t108 * t13 + t116 * t59 + t117 * t60 + t120 * t50 + t121 * t51 + t30 * t55); t365 + m(7) * (t108 * t14 + t116 * t61 + t117 * t62 + t132 * t50 + t133 * t51 + t30 * t56); m(7) * ((t349 * t409 - t30) * t331 + (t108 * t349 + t358 * t50 + t359 * t51 + (-t116 * t358 + t117 * t359) * qJD(1)) * t330); (t108 * t30 + t116 * t51 + t117 * t50) * t566 + (t46 * t331 - t106 + (t358 * t28 + t359 * t29 - t331 * t425) * t349) * t331 + (t359 * t3 + t358 * t4 + t425 * t520 + (-t111 * t349 - t23 * t358 - t24 * t359) * t331 + (t359 * t28 - t358 * t29 + t331 * t426) * qJD(1)) * t330;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
