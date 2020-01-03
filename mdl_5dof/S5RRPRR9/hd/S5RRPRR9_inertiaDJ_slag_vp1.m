% Calculate time derivative of joint inertia matrix for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:20:03
% EndTime: 2019-12-31 20:20:33
% DurationCPUTime: 14.73s
% Computational Cost: add. (34708->903), mult. (38656->1250), div. (0->0), fcn. (37080->10), ass. (0->463)
t343 = cos(qJ(1));
t335 = qJ(2) + pkin(9);
t326 = cos(t335);
t482 = qJD(2) * t326;
t442 = t482 / 0.2e1;
t325 = sin(t335);
t340 = sin(qJ(1));
t485 = qJD(1) * t340;
t455 = t325 * t485;
t586 = t343 * t442 - t455 / 0.2e1;
t484 = qJD(1) * t343;
t443 = t484 / 0.2e1;
t585 = t325 * t443 + t340 * t442;
t480 = qJD(2) * t340;
t448 = t326 * t480;
t353 = t325 * t484 + t448;
t344 = -pkin(8) - pkin(7);
t551 = pkin(7) + t344;
t584 = t326 * t551;
t339 = sin(qJ(2));
t342 = cos(qJ(2));
t540 = Icges(3,4) * t342;
t399 = -Icges(3,2) * t339 + t540;
t268 = Icges(3,6) * t340 + t343 * t399;
t541 = Icges(3,4) * t339;
t405 = Icges(3,1) * t342 - t541;
t270 = Icges(3,5) * t340 + t343 * t405;
t376 = t268 * t339 - t270 * t342;
t583 = t340 * t376;
t538 = Icges(4,4) * t326;
t397 = -Icges(4,2) * t325 + t538;
t251 = Icges(4,6) * t340 + t343 * t397;
t539 = Icges(4,4) * t325;
t403 = Icges(4,1) * t326 - t539;
t253 = Icges(4,5) * t340 + t343 * t403;
t378 = t251 * t325 - t253 * t326;
t582 = t340 * t378;
t322 = pkin(2) * t342 + pkin(1);
t554 = pkin(1) - t322;
t581 = t340 * t554;
t267 = -Icges(3,6) * t343 + t340 * t399;
t269 = -Icges(3,5) * t343 + t340 * t405;
t377 = t267 * t339 - t269 * t342;
t580 = t343 * t377;
t250 = -Icges(4,6) * t343 + t340 * t397;
t252 = -Icges(4,5) * t343 + t340 * t403;
t379 = t250 * t325 - t252 * t326;
t579 = t343 * t379;
t336 = qJ(4) + qJ(5);
t329 = sin(t336);
t330 = cos(t336);
t390 = Icges(6,5) * t330 - Icges(6,6) * t329;
t334 = qJD(4) + qJD(5);
t514 = t325 * t334;
t162 = (-Icges(6,5) * t329 - Icges(6,6) * t330) * t514 + (Icges(6,3) * t325 + t326 * t390) * qJD(2);
t534 = Icges(6,4) * t330;
t394 = -Icges(6,2) * t329 + t534;
t232 = -Icges(6,6) * t326 + t325 * t394;
t526 = t232 * t329;
t578 = -qJD(2) * t526 - t162;
t338 = sin(qJ(4));
t341 = cos(qJ(4));
t391 = Icges(5,5) * t341 - Icges(5,6) * t338;
t476 = qJD(4) * t325;
t183 = (-Icges(5,5) * t338 - Icges(5,6) * t341) * t476 + (Icges(5,3) * t325 + t326 * t391) * qJD(2);
t536 = Icges(5,4) * t341;
t395 = -Icges(5,2) * t338 + t536;
t241 = -Icges(5,6) * t326 + t325 * t395;
t524 = t241 * t338;
t577 = -qJD(2) * t524 - t183;
t331 = t340 * rSges(4,3);
t512 = t325 * t343;
t576 = -rSges(4,2) * t512 + t331;
t392 = Icges(4,5) * t326 - Icges(4,6) * t325;
t248 = -Icges(4,3) * t343 + t340 * t392;
t393 = Icges(3,5) * t342 - Icges(3,6) * t339;
t265 = -Icges(3,3) * t343 + t340 * t393;
t487 = qJD(1) * t326;
t430 = -qJD(4) + t487;
t477 = qJD(2) * t343;
t449 = t325 * t477;
t575 = t340 * t430 + t449;
t435 = -t334 + t487;
t574 = t340 * t435 + t449;
t450 = t325 * t480;
t573 = t343 * t435 - t450;
t508 = t330 * t340;
t509 = t329 * t343;
t263 = -t326 * t509 + t508;
t507 = t330 * t343;
t510 = t329 * t340;
t264 = t326 * t507 + t510;
t176 = t264 * rSges(6,1) + t263 * rSges(6,2) + rSges(6,3) * t512;
t511 = t326 * t343;
t311 = pkin(3) * t511;
t276 = pkin(7) * t512 + t311;
t505 = t338 * t340;
t320 = pkin(4) * t505;
t321 = pkin(4) * t341 + pkin(3);
t368 = t321 * t511 - t344 * t512 + t320;
t194 = t368 - t276;
t496 = t176 + t194;
t261 = -t326 * t510 - t507;
t262 = t326 * t508 - t509;
t414 = -t262 * rSges(6,1) - t261 * rSges(6,2);
t513 = t325 * t340;
t175 = rSges(6,3) * t513 - t414;
t553 = pkin(3) - t321;
t350 = -t325 * t551 - t326 * t553;
t504 = t338 * t343;
t473 = pkin(4) * t504;
t193 = t340 * t350 - t473;
t497 = t175 + t193;
t572 = -t340 * t497 - t343 * t496;
t571 = 2 * m(3);
t570 = 2 * m(4);
t569 = 2 * m(5);
t568 = 2 * m(6);
t567 = t340 ^ 2;
t566 = t343 ^ 2;
t565 = -t326 / 0.2e1;
t564 = t340 / 0.2e1;
t563 = -t343 / 0.2e1;
t562 = t343 / 0.2e1;
t561 = -rSges(5,3) - pkin(7);
t306 = rSges(3,1) * t339 + rSges(3,2) * t342;
t560 = m(3) * t306;
t559 = pkin(2) * t339;
t558 = pkin(3) * t326;
t557 = pkin(6) * t340;
t556 = t325 * pkin(3);
t333 = t343 * pkin(6);
t555 = qJD(1) / 0.2e1;
t337 = -qJ(3) - pkin(6);
t552 = -pkin(6) - t337;
t550 = rSges(3,1) * t342;
t549 = rSges(4,1) * t326;
t548 = rSges(3,2) * t339;
t547 = rSges(3,3) * t343;
t546 = rSges(6,3) * t325;
t436 = -t326 * t334 + qJD(1);
t375 = t340 * t436;
t158 = -t573 * t329 + t330 * t375;
t159 = t329 * t375 + t573 * t330;
t101 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t353;
t103 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t353;
t169 = Icges(6,5) * t262 + Icges(6,6) * t261 + Icges(6,3) * t513;
t171 = Icges(6,4) * t262 + Icges(6,2) * t261 + Icges(6,6) * t513;
t173 = Icges(6,1) * t262 + Icges(6,4) * t261 + Icges(6,5) * t513;
t389 = -t171 * t329 + t173 * t330;
t99 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t353;
t25 = (qJD(2) * t389 - t99) * t326 + (qJD(2) * t169 + (-t171 * t334 + t103) * t330 + (-t173 * t334 - t101) * t329) * t325;
t545 = t25 * t343;
t374 = t343 * t436;
t156 = t574 * t329 + t330 * t374;
t157 = t329 * t374 - t574 * t330;
t447 = t326 * t477;
t352 = t447 - t455;
t100 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t352;
t102 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t352;
t170 = Icges(6,5) * t264 + Icges(6,6) * t263 + Icges(6,3) * t512;
t172 = Icges(6,4) * t264 + Icges(6,2) * t263 + Icges(6,6) * t512;
t174 = Icges(6,1) * t264 + Icges(6,4) * t263 + Icges(6,5) * t512;
t388 = -t172 * t329 + t174 * t330;
t98 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t352;
t26 = (qJD(2) * t388 - t98) * t326 + (qJD(2) * t170 + (-t172 * t334 + t102) * t330 + (-t174 * t334 - t100) * t329) * t325;
t544 = t26 * t340;
t332 = t340 * rSges(3,3);
t543 = -rSges(6,3) + t344;
t415 = t159 * rSges(6,1) + t158 * rSges(6,2);
t105 = rSges(6,3) * t353 + t415;
t542 = t105 * t512 + t175 * t447;
t537 = Icges(5,4) * t338;
t535 = Icges(6,4) * t329;
t502 = t341 * t343;
t278 = -t326 * t505 - t502;
t503 = t340 * t341;
t279 = t326 * t503 - t504;
t417 = -rSges(5,1) * t279 - rSges(5,2) * t278;
t208 = rSges(5,3) * t513 - t417;
t527 = t208 * t343;
t400 = Icges(6,1) * t330 - t535;
t233 = -Icges(6,5) * t326 + t325 * t400;
t525 = t233 * t330;
t523 = t250 * t326;
t522 = t251 * t326;
t521 = t252 * t325;
t520 = t253 * t325;
t519 = t267 * t342;
t518 = t268 * t342;
t517 = t269 * t339;
t516 = t270 * t339;
t515 = t325 * t321;
t184 = (-Icges(5,2) * t341 - t537) * t476 + (Icges(5,6) * t325 + t326 * t395) * qJD(2);
t506 = t338 * t184;
t501 = t343 * t337;
t462 = t157 * rSges(6,1) + t156 * rSges(6,2) + rSges(6,3) * t447;
t104 = -rSges(6,3) * t455 + t462;
t302 = pkin(7) * t447;
t474 = qJD(4) * t341;
t467 = pkin(4) * t474;
t458 = qJD(1) * t473 + t340 * t467 + t344 * t455;
t475 = qJD(4) * t338;
t468 = pkin(4) * t475;
t500 = t104 - t302 + (pkin(7) * t485 + t477 * t553) * t325 + ((-qJD(2) * t344 - t468) * t343 + t553 * t485) * t326 + t458;
t413 = rSges(6,1) * t330 - rSges(6,2) * t329;
t165 = (-rSges(6,1) * t329 - rSges(6,2) * t330) * t514 + (t326 * t413 + t546) * qJD(2);
t446 = t325 * t475;
t201 = -pkin(4) * t446 + qJD(2) * t350;
t499 = -t165 - t201;
t234 = -rSges(6,3) * t326 + t325 * t413;
t483 = qJD(2) * t325;
t498 = t176 * t483 + t234 * t455;
t134 = t326 * t175 + t234 * t513;
t230 = -t325 * t553 + t584;
t495 = t230 + t234;
t246 = t333 + t501 - t581;
t312 = t343 * t322;
t247 = -pkin(1) * t343 + t340 * t552 + t312;
t494 = t340 * t246 + t343 * t247;
t493 = -t247 - t276;
t292 = -t326 * pkin(7) + t556;
t453 = t339 * t485;
t316 = pkin(2) * t453;
t492 = t292 * t485 + t316;
t491 = rSges(4,2) * t455 + rSges(4,3) * t484;
t490 = t343 * t550 + t332;
t249 = Icges(4,3) * t340 + t343 * t392;
t489 = qJD(1) * t249;
t266 = Icges(3,3) * t340 + t343 * t393;
t488 = qJD(1) * t266;
t486 = qJD(1) * t337;
t481 = qJD(2) * t339;
t479 = qJD(2) * t341;
t478 = qJD(2) * t342;
t472 = t343 * t548;
t470 = pkin(2) * t481;
t469 = pkin(2) * t478;
t466 = t334 * t330 * t232;
t240 = -Icges(5,3) * t326 + t325 * t391;
t401 = Icges(5,1) * t341 - t537;
t242 = -Icges(5,5) * t326 + t325 * t401;
t125 = t240 * t513 + t241 * t278 + t242 * t279;
t202 = Icges(5,5) * t279 + Icges(5,6) * t278 + Icges(5,3) * t513;
t204 = Icges(5,4) * t279 + Icges(5,2) * t278 + Icges(5,6) * t513;
t206 = Icges(5,1) * t279 + Icges(5,4) * t278 + Icges(5,5) * t513;
t385 = -t204 * t338 + t206 * t341;
t94 = -t202 * t326 + t325 * t385;
t465 = t94 / 0.2e1 + t125 / 0.2e1;
t280 = -t326 * t504 + t503;
t281 = t326 * t502 + t505;
t126 = t240 * t512 + t241 * t280 + t242 * t281;
t203 = Icges(5,5) * t281 + Icges(5,6) * t280 + Icges(5,3) * t512;
t205 = Icges(5,4) * t281 + Icges(5,2) * t280 + Icges(5,6) * t512;
t207 = Icges(5,1) * t281 + Icges(5,4) * t280 + Icges(5,5) * t512;
t384 = -t205 * t338 + t207 * t341;
t95 = -t203 * t326 + t325 * t384;
t464 = t95 / 0.2e1 + t126 / 0.2e1;
t164 = (-Icges(6,1) * t329 - t534) * t514 + (Icges(6,5) * t325 + t326 * t400) * qJD(2);
t231 = -Icges(6,3) * t326 + t325 * t390;
t463 = t325 * t330 * t164 + t231 * t483 + t482 * t525;
t185 = (-Icges(5,1) * t338 - t536) * t476 + (Icges(5,5) * t325 + t326 * t401) * qJD(2);
t461 = t325 * t341 * t185 + t326 * t242 * t479 + t240 * t483;
t431 = -qJD(4) * t326 + qJD(1);
t372 = t341 * t431;
t186 = t575 * t338 + t343 * t372;
t371 = t431 * t338;
t187 = -t575 * t341 + t343 * t371;
t460 = t187 * rSges(5,1) + t186 * rSges(5,2) + rSges(5,3) * t447;
t327 = qJD(3) * t340;
t457 = qJD(3) * t343 + t337 * t485 + t340 * t470;
t459 = t340 * ((-t343 * t554 - t557) * qJD(1) - t457) + t343 * (-t343 * t470 + t327 + (t343 * t552 + t581) * qJD(1)) + t246 * t484;
t209 = t281 * rSges(5,1) + t280 * rSges(5,2) + rSges(5,3) * t512;
t416 = rSges(5,1) * t341 - rSges(5,2) * t338;
t243 = -rSges(5,3) * t326 + t325 * t416;
t456 = t243 * t485;
t445 = t513 / 0.2e1;
t444 = t512 / 0.2e1;
t290 = rSges(4,1) * t325 + rSges(4,2) * t326;
t441 = -t290 - t559;
t440 = -t292 - t559;
t439 = t343 * t495;
t438 = -t340 * t337 + t312;
t437 = -t321 * t326 - t322;
t434 = t326 * t468;
t433 = t343 * t467;
t432 = t326 * t105 + t165 * t513 + t353 * t234;
t422 = pkin(7) * t325 + t558;
t275 = t422 * t340;
t429 = t340 * t275 + t343 * t276 + t494;
t424 = -t243 + t440;
t423 = -t422 * qJD(2) - t469;
t128 = -t231 * t326 + (t525 - t526) * t325;
t88 = -t169 * t326 + t325 * t389;
t89 = -t170 * t326 + t325 * t388;
t409 = t88 * t340 + t89 * t343;
t110 = t231 * t513 + t232 * t261 + t233 * t262;
t75 = t169 * t513 + t171 * t261 + t173 * t262;
t76 = t170 * t513 + t172 * t261 + t174 * t262;
t412 = t340 * t75 + t343 * t76;
t41 = -t110 * t326 + t325 * t412;
t111 = t231 * t512 + t232 * t263 + t233 * t264;
t77 = t169 * t512 + t171 * t263 + t173 * t264;
t78 = t170 * t512 + t172 * t263 + t174 * t264;
t411 = t340 * t77 + t343 * t78;
t42 = -t111 * t326 + t325 * t411;
t19 = t101 * t263 + t103 * t264 + t156 * t171 + t157 * t173 + t169 * t352 + t512 * t99;
t20 = t100 * t263 + t102 * t264 + t156 * t172 + t157 * t174 + t170 * t352 + t512 * t98;
t163 = (-Icges(6,2) * t330 - t535) * t514 + (Icges(6,6) * t325 + t326 * t394) * qJD(2);
t49 = t156 * t232 + t157 * t233 + t162 * t512 + t163 * t263 + t164 * t264 + t231 * t352;
t57 = t340 * t78 - t343 * t77;
t5 = (qJD(2) * t411 - t49) * t326 + (-qJD(1) * t57 + qJD(2) * t111 + t19 * t340 + t20 * t343) * t325;
t21 = t101 * t261 + t103 * t262 + t158 * t171 + t159 * t173 + t169 * t353 + t513 * t99;
t22 = t100 * t261 + t102 * t262 + t158 * t172 + t159 * t174 + t170 * t353 + t513 * t98;
t50 = t158 * t232 + t159 * t233 + t162 * t513 + t163 * t261 + t164 * t262 + t231 * t353;
t56 = t340 * t76 - t343 * t75;
t6 = (qJD(2) * t412 - t50) * t326 + (-qJD(1) * t56 + qJD(2) * t110 + t21 * t340 + t22 * t343) * t325;
t421 = t42 * t447 + t5 * t512 + t6 * t513 + (-t128 * t326 + t325 * t409) * t483 + t353 * t41;
t420 = -t548 + t550;
t419 = -rSges(4,2) * t325 + t549;
t188 = t340 * t372 + (-t343 * t430 + t450) * t338;
t189 = t430 * t502 + (-t325 * t479 + t371) * t340;
t418 = rSges(5,1) * t189 + rSges(5,2) * t188;
t410 = t340 * t89 - t343 * t88;
t90 = t202 * t513 + t204 * t278 + t206 * t279;
t91 = t203 * t513 + t205 * t278 + t207 * t279;
t62 = t340 * t91 - t343 * t90;
t408 = t340 * t90 + t343 * t91;
t92 = t202 * t512 + t204 * t280 + t206 * t281;
t93 = t203 * t512 + t205 * t280 + t207 * t281;
t63 = t340 * t93 - t343 * t92;
t407 = t340 * t92 + t343 * t93;
t406 = t340 * t94 + t343 * t95;
t404 = Icges(3,1) * t339 + t540;
t402 = Icges(4,1) * t325 + t538;
t398 = Icges(3,2) * t342 + t541;
t396 = Icges(4,2) * t326 + t539;
t383 = -t209 * t340 + t527;
t382 = -t208 * t340 - t209 * t343;
t373 = t440 - t495;
t255 = rSges(4,1) * t511 + t576;
t190 = (-rSges(5,1) * t338 - rSges(5,2) * t341) * t476 + (rSges(5,3) * t325 + t326 * t416) * qJD(2);
t370 = -t190 + t423;
t369 = -pkin(1) - t420;
t182 = t424 * t343;
t367 = -t322 - t419;
t301 = pkin(3) * t450;
t365 = t340 * (pkin(7) * t353 + qJD(1) * t311 - t301) + t343 * (-pkin(7) * t455 + t302 + (-t326 * t485 - t449) * pkin(3)) + t275 * t484 + t459;
t364 = t423 + t499;
t363 = qJD(2) * t306;
t362 = qJD(2) * t290;
t361 = t193 * t343 - t340 * t496;
t360 = qJD(2) * t404;
t359 = qJD(2) * t402;
t358 = qJD(2) * t398;
t357 = qJD(2) * t396;
t356 = qJD(2) * (-Icges(3,5) * t339 - Icges(3,6) * t342);
t355 = qJD(2) * (-Icges(4,5) * t325 - Icges(4,6) * t326);
t137 = t373 * t343;
t354 = t325 * t561 - t322 - t558;
t351 = t325 * t543 + t437;
t12 = qJD(1) * t411 - t19 * t343 + t20 * t340;
t13 = qJD(1) * t412 - t21 * t343 + t22 * t340;
t349 = t12 * t444 + t13 * t445 + t5 * t564 + t6 * t563 + (qJD(1) * t409 + t544 - t545) * t565 + t41 * t485 / 0.2e1 + t42 * t443 + t410 * t483 / 0.2e1 + t586 * t57 + t585 * t56;
t348 = rSges(3,2) * t453 + rSges(3,3) * t484 - t343 * t363;
t347 = t340 * t354 - t501;
t127 = t128 * t483;
t60 = t578 * t326 + (-t466 + (-t233 * t334 - t163) * t329) * t325 + t463;
t7 = t127 + (qJD(2) * t409 - t60) * t326 + (-qJD(1) * t410 + t25 * t340 + t26 * t343) * t325;
t346 = -t326 * t7 - t42 * t455 + t421;
t345 = t127 + (t25 + t50) * t445 + (t26 + t49) * t444 + (t111 + t89) * t586 + (t110 + t88) * t585;
t296 = t420 * qJD(2);
t285 = t419 * qJD(2);
t274 = -t472 + t490;
t273 = t340 * t420 - t547;
t254 = -rSges(4,3) * t343 + t340 * t419;
t245 = t441 * t343;
t244 = t441 * t340;
t239 = t557 + (pkin(1) - t548) * t343 + t490;
t238 = t340 * t369 + t333 + t547;
t225 = t255 + t438;
t224 = (rSges(4,3) - t337) * t343 + t367 * t340;
t216 = t340 * t356 + t488;
t215 = -qJD(1) * t265 + t343 * t356;
t196 = t340 * t355 + t489;
t195 = -qJD(1) * t248 + t343 * t355;
t192 = t306 * t480 + ((-rSges(3,3) - pkin(6)) * t340 + t369 * t343) * qJD(1);
t191 = (t333 + (-pkin(1) - t550) * t340) * qJD(1) + t348;
t181 = t424 * t340;
t178 = -t290 * t484 - t285 * t340 + (-t339 * t484 - t340 * t478) * pkin(2);
t177 = t290 * t485 + t316 + (-t285 - t469) * t343;
t161 = t175 * t512;
t151 = t266 * t340 - t376 * t343;
t150 = t265 * t340 - t580;
t149 = -t266 * t343 - t583;
t148 = -t265 * t343 - t340 * t377;
t147 = t290 * t480 + (t343 * t367 - t331) * qJD(1) + t457;
t146 = t327 + (-t322 - t549) * t485 + (qJD(2) * t441 - t486) * t343 + t491;
t145 = t438 + t209 + t276;
t144 = t347 + t417;
t143 = -t209 * t326 - t243 * t512;
t142 = t208 * t326 + t243 * t513;
t141 = t249 * t340 - t378 * t343;
t140 = t248 * t340 - t579;
t139 = -t249 * t343 - t582;
t138 = -t248 * t343 - t340 * t379;
t136 = t373 * t340;
t135 = -t176 * t326 - t234 * t512;
t133 = -t240 * t326 + (t242 * t341 - t524) * t325;
t132 = t368 + t438 + t176;
t131 = (pkin(4) * t338 - t337) * t343 + t351 * t340 + t414;
t130 = t133 * t483;
t129 = t383 * t325;
t124 = -t433 + t301 + (-t434 + (-t515 - t584) * qJD(2)) * t340 + (t343 * t350 + t320) * qJD(1);
t122 = -t176 * t513 + t161;
t121 = rSges(5,3) * t353 + t418;
t120 = -rSges(5,3) * t455 + t460;
t119 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t353;
t118 = Icges(5,1) * t187 + Icges(5,4) * t186 + Icges(5,5) * t352;
t117 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t353;
t116 = Icges(5,4) * t187 + Icges(5,2) * t186 + Icges(5,6) * t352;
t115 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t353;
t114 = Icges(5,5) * t187 + Icges(5,6) * t186 + Icges(5,3) * t352;
t113 = qJD(1) * t182 + t340 * t370;
t112 = t343 * t370 + t456 + t492;
t87 = -t325 * t439 - t326 * t496;
t86 = t193 * t326 + t230 * t513 + t134;
t85 = t354 * t484 + t448 * t561 + t301 - t418 + t457;
t84 = t302 + t327 + (-t556 - t559) * t477 + t347 * qJD(1) + t460;
t83 = -t382 + t429;
t74 = t325 * t361 + t161;
t73 = qJD(1) * t137 + t340 * t364;
t72 = t343 * t364 + t485 * t495 + t492;
t71 = t433 + (t434 + (t326 * t543 + t515) * qJD(2)) * t340 + (t343 * t351 - t320) * qJD(1) - t415 + t457;
t70 = t327 + (t437 - t546) * t485 + (-t434 - t486 + (-t326 * t344 - t515 - t559) * qJD(2)) * t343 + t458 + t462;
t69 = (t243 * t480 + t121) * t326 + (-qJD(2) * t208 + t190 * t340 + t243 * t484) * t325;
t68 = (-t243 * t477 - t120) * t326 + (qJD(2) * t209 - t190 * t343 + t456) * t325;
t67 = t429 - t572;
t66 = t577 * t326 + (-t506 + (-t241 * t341 - t242 * t338) * qJD(4)) * t325 + t461;
t65 = -t175 * t483 + t432;
t64 = -t165 * t512 + (-t234 * t477 - t104) * t326 + t498;
t59 = t183 * t513 + t184 * t278 + t185 * t279 + t188 * t241 + t189 * t242 + t240 * t353;
t58 = t183 * t512 + t184 * t280 + t185 * t281 + t186 * t241 + t187 * t242 + t240 * t352;
t51 = t383 * t482 + (qJD(1) * t382 - t120 * t340 + t121 * t343) * t325;
t46 = -t126 * t326 + t325 * t407;
t45 = -t125 * t326 + t325 * t408;
t43 = -t176 * t448 + (-t104 * t340 + (-t175 * t340 - t176 * t343) * qJD(1)) * t325 + t542;
t35 = t120 * t343 + t121 * t340 + (t527 + (-t209 + t493) * t340) * qJD(1) + t365;
t34 = (qJD(2) * t384 - t114) * t326 + (qJD(2) * t203 - t116 * t338 + t118 * t341 + (-t205 * t341 - t207 * t338) * qJD(4)) * t325;
t33 = (qJD(2) * t385 - t115) * t326 + (qJD(2) * t202 - t117 * t338 + t119 * t341 + (-t204 * t341 - t206 * t338) * qJD(4)) * t325;
t32 = (t230 * t480 + t124) * t326 + (-qJD(2) * t497 + t201 * t340 + t230 * t484) * t325 + t432;
t31 = (-qJD(2) * t439 - t500) * t326 + (qJD(2) * t194 + t230 * t485 + t343 * t499) * t325 + t498;
t30 = t114 * t513 + t116 * t278 + t118 * t279 + t188 * t205 + t189 * t207 + t203 * t353;
t29 = t115 * t513 + t117 * t278 + t119 * t279 + t188 * t204 + t189 * t206 + t202 * t353;
t28 = t114 * t512 + t116 * t280 + t118 * t281 + t186 * t205 + t187 * t207 + t203 * t352;
t27 = t115 * t512 + t117 * t280 + t119 * t281 + t186 * t204 + t187 * t206 + t202 * t352;
t18 = t500 * t343 + (t105 + t124) * t340 + (t497 * t343 + (t493 - t496) * t340) * qJD(1) + t365;
t17 = t361 * t482 + (qJD(1) * t572 + t124 * t343 - t500 * t340) * t325 + t542;
t16 = qJD(1) * t408 - t29 * t343 + t30 * t340;
t15 = qJD(1) * t407 - t27 * t343 + t28 * t340;
t9 = (qJD(2) * t408 - t59) * t326 + (-qJD(1) * t62 + qJD(2) * t125 + t29 * t340 + t30 * t343) * t325;
t8 = (qJD(2) * t407 - t58) * t326 + (-qJD(1) * t63 + qJD(2) * t126 + t27 * t340 + t28 * t343) * t325;
t1 = [t463 - t329 * t233 * t514 + t461 - t242 * t446 + (t191 * t239 + t192 * t238) * t571 + (t146 * t225 + t147 * t224) * t570 + (t144 * t85 + t145 * t84) * t569 + (t131 * t71 + t132 * t70) * t568 + (-t396 + t403) * t483 + (t397 + t402) * t482 + (-t398 + t405) * t481 + (t399 + t404) * t478 + (t577 + t578) * t326 + (-t329 * t163 - t241 * t474 - t466 - t506) * t325; -t545 / 0.2e1 + t544 / 0.2e1 + m(3) * ((-t191 * t340 - t192 * t343) * t306 + (-t238 * t343 - t239 * t340) * t296) + m(4) * (t146 * t244 + t147 * t245 + t177 * t224 + t178 * t225) + m(5) * (t112 * t144 + t113 * t145 + t181 * t84 + t182 * t85) + m(6) * (t131 * t72 + t132 * t73 + t136 * t70 + t137 * t71) + ((t518 / 0.2e1 + t516 / 0.2e1 - t239 * t560 + t522 / 0.2e1 + t520 / 0.2e1 + t111 / 0.2e1 + t89 / 0.2e1 + t464) * t343 + (t238 * t560 + t519 / 0.2e1 + t517 / 0.2e1 + t523 / 0.2e1 + t521 / 0.2e1 + t110 / 0.2e1 + t88 / 0.2e1 + t465) * t340) * qJD(1) + (t393 + t392) * qJD(2) * (t567 / 0.2e1 + t566 / 0.2e1) + ((-t250 * qJD(1) - t343 * t357) * t326 + (-t252 * qJD(1) - t343 * t359) * t325 + (-t267 * qJD(1) - t343 * t358) * t342 + (-t269 * qJD(1) - t343 * t360) * t339 + t34 + t49 + t58 + (-t376 - t378) * qJD(2)) * t564 + ((qJD(1) * t251 - t340 * t357) * t326 + (qJD(1) * t253 - t340 * t359) * t325 + (qJD(1) * t268 - t340 * t358) * t342 + (qJD(1) * t270 - t340 * t360) * t339 + t33 + t50 + t59 + (-t377 - t379) * qJD(2)) * t563; (t136 * t73 + t137 * t72 + t18 * t67) * t568 + t340 * t12 - t343 * t13 + (t112 * t182 + t113 * t181 + t35 * t83) * t569 + t340 * t15 - t343 * t16 + (t245 * t177 + t244 * t178 + (t254 * t340 + t255 * t343 + t494) * ((qJD(1) * t254 - t343 * t362 + t491) * t343 + (-t340 * t362 + (-t247 - t255 + t576) * qJD(1)) * t340 + t459)) * t570 + ((t273 * t340 + t274 * t343) * ((qJD(1) * t273 + t348) * t343 + (-t340 * t363 + (-t274 - t472 + t332) * qJD(1)) * t340) + (t566 + t567) * t306 * t296) * t571 - t343 * ((t343 * t196 + (t139 + t579) * qJD(1)) * t343 + (t138 * qJD(1) + (-t251 * t482 - t253 * t483 + t489) * t340 + (-t195 + (t521 + t523) * qJD(2) - t378 * qJD(1)) * t343) * t340) - t343 * ((t343 * t216 + (t149 + t580) * qJD(1)) * t343 + (t148 * qJD(1) + (-t268 * t478 - t270 * t481 + t488) * t340 + (-t215 + (t517 + t519) * qJD(2) - t376 * qJD(1)) * t343) * t340) + t340 * ((t340 * t215 + (t150 + t583) * qJD(1)) * t340 + (t151 * qJD(1) + (t267 * t478 + t269 * t481) * t343 + (-t216 + (-t516 - t518) * qJD(2) + (t266 - t377) * qJD(1)) * t340) * t343) + t340 * ((t340 * t195 + (t140 + t582) * qJD(1)) * t340 + (t141 * qJD(1) + (t250 * t482 + t252 * t483) * t343 + (-t196 + (-t520 - t522) * qJD(2) + (t249 - t379) * qJD(1)) * t340) * t343) + (t56 + t62 + (-t138 - t148) * t343 + (t139 + t149) * t340) * t485 + (t57 + t63 + (-t140 - t150) * t343 + (t141 + t151) * t340) * t484; m(6) * (t340 * t71 - t343 * t70 + (t131 * t343 + t132 * t340) * qJD(1)) + m(5) * (t340 * t85 - t343 * t84 + (t144 * t343 + t145 * t340) * qJD(1)) + m(4) * (-t146 * t343 + t147 * t340 + (t224 * t343 + t225 * t340) * qJD(1)); m(6) * (t340 * t72 - t343 * t73 + (t136 * t340 + t137 * t343) * qJD(1)) + m(5) * (t112 * t340 - t113 * t343 + (t181 * t340 + t182 * t343) * qJD(1)) + m(4) * (t177 * t340 - t178 * t343 + (t244 * t340 + t245 * t343) * qJD(1)); 0; t345 + ((t34 / 0.2e1 + t58 / 0.2e1) * t343 + (t59 / 0.2e1 + t33 / 0.2e1) * t340 + (-t340 * t464 + t343 * t465) * qJD(1)) * t325 + (-t60 - t66 + (t340 * t465 + t343 * t464) * qJD(2)) * t326 + m(5) * (t142 * t85 + t143 * t84 + t144 * t69 + t145 * t68) + m(6) * (t131 * t32 + t132 * t31 + t70 * t87 + t71 * t86) + t130; (qJD(2) * (t340 * t95 - t343 * t94) / 0.2e1 + t15 * t562 + t16 * t564 + (t62 * t562 - t340 * t63 / 0.2e1) * qJD(1)) * t325 + (t63 * t442 + t46 * t555 + (qJD(1) * t95 - t33) * t565 - t9 / 0.2e1) * t343 + (t62 * t442 + (qJD(1) * t94 + t34) * t565 + t8 / 0.2e1 + t45 * t555) * t340 + m(5) * (t112 * t142 + t113 * t143 + t129 * t35 + t181 * t68 + t182 * t69 + t51 * t83) + m(6) * (t136 * t31 + t137 * t32 + t17 * t67 + t18 * t74 + t72 * t86 + t73 * t87) + t349; m(5) * (t340 * t69 - t343 * t68 + (t142 * t343 + t143 * t340) * qJD(1)) + m(6) * (-t31 * t343 + t32 * t340 + (t340 * t87 + t343 * t86) * qJD(1)); (t17 * t74 + t31 * t87 + t32 * t86) * t568 + (t129 * t51 + t142 * t69 + t143 * t68) * t569 + (t66 * t326 - t130 - t7 + (-t326 * t406 + t340 * t45 + t343 * t46) * qJD(2)) * t326 + (t343 * t8 + t340 * t9 - t326 * (t33 * t340 + t34 * t343) + (-t133 * t326 + t325 * t406) * qJD(2) + ((-t326 * t94 + t45) * t343 + (t326 * t95 - t42 - t46) * t340) * qJD(1)) * t325 + t421; t345 - t60 * t326 + m(6) * (t131 * t65 + t132 * t64 + t134 * t71 + t135 * t70); t349 + m(6) * (t122 * t18 + t134 * t72 + t135 * t73 + t136 * t64 + t137 * t65 + t43 * t67); m(6) * (t340 * t65 - t343 * t64 + (t134 * t343 + t135 * t340) * qJD(1)); m(6) * (t122 * t17 + t134 * t32 + t135 * t31 + t43 * t74 + t64 * t87 + t65 * t86) + t346; (t122 * t43 + t134 * t65 + t135 * t64) * t568 + t346;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
