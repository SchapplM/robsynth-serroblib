% Calculate time derivative of joint inertia matrix for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:37:00
% EndTime: 2019-03-09 03:37:27
% DurationCPUTime: 19.64s
% Computational Cost: add. (49939->892), mult. (39926->1215), div. (0->0), fcn. (38078->12), ass. (0->452)
t608 = Icges(4,3) + Icges(5,3);
t340 = qJ(3) + pkin(11);
t332 = sin(t340);
t334 = cos(t340);
t345 = sin(qJ(3));
t348 = cos(qJ(3));
t607 = Icges(4,5) * t348 + Icges(5,5) * t334 - Icges(4,6) * t345 - Icges(5,6) * t332;
t341 = qJ(1) + pkin(10);
t333 = sin(t341);
t335 = cos(t341);
t603 = t333 * t608 + t335 * t607;
t545 = Icges(5,4) * t334;
t404 = -Icges(5,2) * t332 + t545;
t243 = Icges(5,6) * t333 + t335 * t404;
t546 = Icges(5,4) * t332;
t410 = Icges(5,1) * t334 - t546;
t245 = Icges(5,5) * t333 + t335 * t410;
t547 = Icges(4,4) * t348;
t406 = -Icges(4,2) * t345 + t547;
t264 = Icges(4,6) * t333 + t335 * t406;
t548 = Icges(4,4) * t345;
t412 = Icges(4,1) * t348 - t548;
t267 = Icges(4,5) * t333 + t335 * t412;
t582 = t243 * t332 - t245 * t334 + t264 * t345 - t267 * t348;
t606 = t333 * t582;
t605 = t582 * t335;
t604 = t333 * t607 - t335 * t608;
t602 = (-Icges(4,5) * t345 - Icges(5,5) * t332 - Icges(4,6) * t348 - Icges(5,6) * t334) * qJD(3);
t242 = -Icges(5,6) * t335 + t333 * t404;
t244 = -Icges(5,5) * t335 + t333 * t410;
t263 = -Icges(4,6) * t335 + t333 * t406;
t266 = -Icges(4,5) * t335 + t333 * t412;
t601 = t242 * t332 - t244 * t334 + t263 * t345 - t266 * t348;
t342 = qJ(5) + qJ(6);
t336 = sin(t342);
t337 = cos(t342);
t397 = Icges(7,5) * t337 - Icges(7,6) * t336;
t339 = qJD(5) + qJD(6);
t521 = t332 * t339;
t180 = (-Icges(7,5) * t336 - Icges(7,6) * t337) * t521 + (Icges(7,3) * t332 + t334 * t397) * qJD(3);
t541 = Icges(7,4) * t337;
t401 = -Icges(7,2) * t336 + t541;
t250 = -Icges(7,6) * t334 + t332 * t401;
t489 = qJD(3) * t336;
t600 = -t250 * t489 - t180;
t344 = sin(qJ(5));
t347 = cos(qJ(5));
t398 = Icges(6,5) * t347 - Icges(6,6) * t344;
t483 = qJD(5) * t332;
t213 = (-Icges(6,5) * t344 - Icges(6,6) * t347) * t483 + (Icges(6,3) * t332 + t334 * t398) * qJD(3);
t543 = Icges(6,4) * t347;
t402 = -Icges(6,2) * t344 + t543;
t262 = -Icges(6,6) * t334 + t332 * t402;
t487 = qJD(3) * t344;
t599 = -t262 * t487 - t213;
t491 = qJD(3) * t334;
t451 = t491 / 0.2e1;
t497 = qJD(1) * t333;
t462 = t332 * t497;
t598 = t335 * t451 - t462 / 0.2e1;
t495 = qJD(1) * t335;
t452 = t495 / 0.2e1;
t597 = t332 * t452 + t333 * t451;
t457 = t333 * t491;
t363 = t332 * t495 + t457;
t596 = t603 * qJD(1);
t595 = t604 * t333 - t606 + (-t601 - t603) * t335;
t575 = t333 ^ 2;
t574 = t335 ^ 2;
t578 = 2 * m(5);
t331 = pkin(3) * t348 + pkin(2);
t310 = t335 * t331;
t343 = -qJ(4) - pkin(7);
t559 = -pkin(7) - t343;
t237 = -pkin(2) * t335 + t333 * t559 + t310;
t556 = rSges(5,1) * t334;
t427 = -rSges(5,2) * t332 + t556;
t246 = -rSges(5,3) * t335 + t333 * t427;
t516 = t334 * t335;
t327 = t333 * rSges(5,3);
t522 = t332 * t335;
t584 = -rSges(5,2) * t522 + t327;
t247 = rSges(5,1) * t516 + t584;
t294 = rSges(5,1) * t332 + rSges(5,2) * t334;
t373 = qJD(3) * t294;
t329 = t335 * pkin(7);
t513 = t335 * t343;
t561 = pkin(2) - t331;
t592 = t333 * t561;
t236 = t329 + t513 - t592;
t325 = qJD(4) * t333;
t486 = qJD(3) * t345;
t477 = pkin(3) * t486;
t464 = qJD(4) * t335 + t333 * t477 + t343 * t497;
t565 = pkin(7) * t333;
t467 = t333 * ((-t335 * t561 - t565) * qJD(1) - t464) + t335 * (-t335 * t477 + t325 + (t335 * t559 + t592) * qJD(1)) + t236 * t495;
t501 = rSges(5,2) * t462 + rSges(5,3) * t495;
t67 = (qJD(1) * t246 - t335 * t373 + t501) * t335 + (-t333 * t373 + (-t237 - t247 + t584) * qJD(1)) * t333 + t467;
t594 = t578 * t67;
t555 = rSges(4,2) * t345;
t557 = rSges(4,1) * t348;
t428 = -t555 + t557;
t554 = rSges(4,3) * t335;
t269 = t333 * t428 - t554;
t479 = t335 * t555;
t328 = t333 * rSges(4,3);
t500 = t335 * t557 + t328;
t270 = -t479 + t500;
t320 = rSges(4,1) * t345 + rSges(4,2) * t348;
t374 = qJD(3) * t320;
t494 = qJD(1) * t345;
t460 = t333 * t494;
t354 = rSges(4,2) * t460 + rSges(4,3) * t495 - t335 * t374;
t107 = (qJD(1) * t269 + t354) * t335 + (-t333 * t374 + (-t270 - t479 + t328) * qJD(1)) * t333;
t579 = 2 * m(4);
t593 = t107 * t579;
t350 = -pkin(9) - pkin(8);
t558 = pkin(8) + t350;
t591 = t334 * t558;
t590 = t333 * t601 + t335 * t604;
t587 = t333 * t603 - t605;
t586 = -qJD(1) * t604 + t335 * t602;
t585 = -t333 * t602 - t596;
t563 = sin(qJ(1)) * pkin(1);
t583 = t329 - t563;
t515 = t335 * t336;
t519 = t333 * t337;
t273 = -t334 * t515 + t519;
t514 = t335 * t337;
t520 = t333 * t336;
t274 = t334 * t514 + t520;
t175 = rSges(7,1) * t274 + rSges(7,2) * t273 + rSges(7,3) * t522;
t309 = pkin(4) * t516;
t278 = pkin(8) * t522 + t309;
t518 = t333 * t344;
t322 = pkin(5) * t518;
t330 = pkin(5) * t347 + pkin(4);
t379 = t330 * t516 - t350 * t522 + t322;
t191 = t379 - t278;
t507 = t175 + t191;
t271 = -t334 * t520 - t514;
t272 = t334 * t519 - t515;
t422 = -rSges(7,1) * t272 - rSges(7,2) * t271;
t523 = t332 * t333;
t174 = rSges(7,3) * t523 - t422;
t560 = pkin(4) - t330;
t360 = -t332 * t558 - t334 * t560;
t512 = t335 * t344;
t480 = pkin(5) * t512;
t190 = t333 * t360 - t480;
t508 = t174 + t190;
t580 = -t333 * t508 - t335 * t507;
t577 = 2 * m(6);
t576 = 2 * m(7);
t573 = t333 / 0.2e1;
t572 = -t334 / 0.2e1;
t571 = -t335 / 0.2e1;
t570 = t335 / 0.2e1;
t569 = -rSges(6,3) - pkin(8);
t568 = m(4) * t320;
t567 = pkin(3) * t345;
t566 = pkin(4) * t334;
t564 = t332 * pkin(4);
t338 = cos(qJ(1)) * pkin(1);
t562 = qJD(1) / 0.2e1;
t553 = rSges(7,3) * t332;
t446 = -t334 * t339 + qJD(1);
t359 = t332 * t489 + t337 * t446;
t496 = qJD(1) * t334;
t445 = -t339 + t496;
t160 = t333 * t359 - t445 * t515;
t488 = qJD(3) * t337;
t358 = -t332 * t488 + t336 * t446;
t161 = t333 * t358 + t445 * t514;
t100 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t363;
t102 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t363;
t104 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t363;
t168 = Icges(7,5) * t272 + Icges(7,6) * t271 + Icges(7,3) * t523;
t170 = Icges(7,4) * t272 + Icges(7,2) * t271 + Icges(7,6) * t523;
t172 = Icges(7,1) * t272 + Icges(7,4) * t271 + Icges(7,5) * t523;
t396 = -t170 * t336 + t172 * t337;
t25 = (qJD(3) * t396 - t100) * t334 + (qJD(3) * t168 + (-t170 * t339 + t104) * t337 + (-t172 * t339 - t102) * t336) * t332;
t552 = t25 * t335;
t158 = t335 * t359 + t445 * t520;
t159 = t335 * t358 - t445 * t519;
t490 = qJD(3) * t335;
t456 = t334 * t490;
t362 = t456 - t462;
t101 = Icges(7,4) * t159 + Icges(7,2) * t158 + Icges(7,6) * t362;
t103 = Icges(7,1) * t159 + Icges(7,4) * t158 + Icges(7,5) * t362;
t169 = Icges(7,5) * t274 + Icges(7,6) * t273 + Icges(7,3) * t522;
t171 = Icges(7,4) * t274 + Icges(7,2) * t273 + Icges(7,6) * t522;
t173 = Icges(7,1) * t274 + Icges(7,4) * t273 + Icges(7,5) * t522;
t395 = -t171 * t336 + t173 * t337;
t99 = Icges(7,5) * t159 + Icges(7,6) * t158 + Icges(7,3) * t362;
t26 = (qJD(3) * t395 - t99) * t334 + (qJD(3) * t169 + (-t171 * t339 + t103) * t337 + (-t173 * t339 - t101) * t336) * t332;
t551 = t26 * t333;
t550 = -rSges(7,3) + t350;
t423 = t161 * rSges(7,1) + t160 * rSges(7,2);
t106 = rSges(7,3) * t363 + t423;
t549 = t106 * t522 + t174 * t456;
t544 = Icges(6,4) * t344;
t542 = Icges(7,4) * t336;
t511 = t335 * t347;
t279 = -t334 * t518 - t511;
t517 = t333 * t347;
t280 = t334 * t517 - t512;
t425 = -rSges(6,1) * t280 - rSges(6,2) * t279;
t207 = rSges(6,3) * t523 - t425;
t534 = t207 * t335;
t216 = (-Icges(6,2) * t347 - t544) * t483 + (Icges(6,6) * t332 + t334 * t402) * qJD(3);
t533 = t216 * t344;
t532 = t242 * t334;
t531 = t243 * t334;
t530 = t244 * t332;
t529 = t245 * t332;
t528 = t263 * t348;
t527 = t264 * t348;
t526 = t266 * t345;
t525 = t267 * t345;
t524 = t332 * t330;
t470 = rSges(7,1) * t159 + rSges(7,2) * t158 + rSges(7,3) * t456;
t105 = -rSges(7,3) * t462 + t470;
t300 = pkin(8) * t456;
t481 = qJD(5) * t347;
t474 = pkin(5) * t481;
t465 = qJD(1) * t480 + t333 * t474 + t350 * t462;
t482 = qJD(5) * t344;
t475 = pkin(5) * t482;
t510 = t105 - t300 + (pkin(8) * t497 + t490 * t560) * t332 + ((-qJD(3) * t350 - t475) * t335 + t560 * t497) * t334 + t465;
t421 = rSges(7,1) * t337 - rSges(7,2) * t336;
t252 = -rSges(7,3) * t334 + t332 * t421;
t493 = qJD(3) * t332;
t509 = t175 * t493 + t252 * t462;
t187 = (-rSges(7,1) * t336 - rSges(7,2) * t337) * t521 + (t334 * t421 + t553) * qJD(3);
t455 = t332 * t482;
t223 = -pkin(5) * t455 + qJD(3) * t360;
t506 = -t187 - t223;
t138 = t174 * t334 + t252 * t523;
t505 = t236 * t333 + t237 * t335;
t504 = -t237 - t278;
t239 = -t332 * t560 + t591;
t503 = t239 + t252;
t295 = -pkin(8) * t334 + t564;
t314 = pkin(3) * t460;
t502 = t295 * t497 + t314;
t492 = qJD(3) * t333;
t485 = qJD(3) * t347;
t484 = qJD(3) * t348;
t476 = pkin(3) * t484;
t473 = t250 * t337 * t339;
t259 = -Icges(6,3) * t334 + t332 * t398;
t408 = Icges(6,1) * t347 - t544;
t265 = -Icges(6,5) * t334 + t332 * t408;
t128 = t259 * t523 + t262 * t279 + t265 * t280;
t201 = Icges(6,5) * t280 + Icges(6,6) * t279 + Icges(6,3) * t523;
t203 = Icges(6,4) * t280 + Icges(6,2) * t279 + Icges(6,6) * t523;
t205 = Icges(6,1) * t280 + Icges(6,4) * t279 + Icges(6,5) * t523;
t392 = -t203 * t344 + t205 * t347;
t95 = -t201 * t334 + t332 * t392;
t472 = t95 / 0.2e1 + t128 / 0.2e1;
t281 = -t334 * t512 + t517;
t282 = t334 * t511 + t518;
t129 = t259 * t522 + t262 * t281 + t265 * t282;
t202 = Icges(6,5) * t282 + Icges(6,6) * t281 + Icges(6,3) * t522;
t204 = Icges(6,4) * t282 + Icges(6,2) * t281 + Icges(6,6) * t522;
t206 = Icges(6,1) * t282 + Icges(6,4) * t281 + Icges(6,5) * t522;
t391 = -t204 * t344 + t206 * t347;
t96 = -t202 * t334 + t332 * t391;
t471 = t96 / 0.2e1 + t129 / 0.2e1;
t407 = Icges(7,1) * t337 - t542;
t182 = (-Icges(7,1) * t336 - t541) * t521 + (Icges(7,5) * t332 + t334 * t407) * qJD(3);
t249 = -Icges(7,3) * t334 + t332 * t397;
t251 = -Icges(7,5) * t334 + t332 * t407;
t469 = t182 * t332 * t337 + t251 * t334 * t488 + t249 * t493;
t441 = -qJD(5) * t334 + qJD(1);
t357 = t332 * t487 + t347 * t441;
t440 = -qJD(5) + t496;
t183 = t335 * t357 + t440 * t518;
t356 = -t332 * t485 + t344 * t441;
t184 = t335 * t356 - t440 * t517;
t468 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t456;
t219 = (-Icges(6,1) * t344 - t543) * t483 + (Icges(6,5) * t332 + t334 * t408) * qJD(3);
t466 = t219 * t332 * t347 + t265 * t334 * t485 + t259 * t493;
t208 = rSges(6,1) * t282 + rSges(6,2) * t281 + rSges(6,3) * t522;
t424 = rSges(6,1) * t347 - rSges(6,2) * t344;
t268 = -rSges(6,3) * t334 + t332 * t424;
t463 = t268 * t497;
t454 = t523 / 0.2e1;
t453 = t522 / 0.2e1;
t450 = -t294 - t567;
t449 = -t295 - t567;
t448 = t335 * t503;
t447 = -t330 * t334 - t331;
t444 = t334 * t475;
t443 = t335 * t474;
t442 = t334 * t106 + t187 * t523 + t252 * t363;
t432 = pkin(8) * t332 + t566;
t277 = t432 * t333;
t439 = t277 * t333 + t278 * t335 + t505;
t434 = -t268 + t449;
t433 = -qJD(3) * t432 - t476;
t431 = -t333 * t343 + t310 + t338;
t254 = t450 * t335;
t137 = -t249 * t334 + (-t250 * t336 + t251 * t337) * t332;
t121 = t249 * t523 + t250 * t271 + t251 * t272;
t74 = t168 * t523 + t170 * t271 + t172 * t272;
t75 = t169 * t523 + t171 * t271 + t173 * t272;
t419 = t333 * t74 + t335 * t75;
t41 = -t121 * t334 + t332 * t419;
t91 = -t168 * t334 + t332 * t396;
t92 = -t169 * t334 + t332 * t395;
t414 = t333 * t91 + t335 * t92;
t122 = t249 * t522 + t250 * t273 + t251 * t274;
t76 = t168 * t522 + t170 * t273 + t172 * t274;
t77 = t169 * t522 + t171 * t273 + t173 * t274;
t418 = t333 * t76 + t335 * t77;
t42 = -t122 * t334 + t332 * t418;
t19 = t100 * t522 + t102 * t273 + t104 * t274 + t158 * t170 + t159 * t172 + t168 * t362;
t20 = t101 * t273 + t103 * t274 + t158 * t171 + t159 * t173 + t169 * t362 + t522 * t99;
t181 = (-Icges(7,2) * t337 - t542) * t521 + (Icges(7,6) * t332 + t334 * t401) * qJD(3);
t54 = t158 * t250 + t159 * t251 + t180 * t522 + t181 * t273 + t182 * t274 + t249 * t362;
t57 = t333 * t77 - t335 * t76;
t5 = (qJD(3) * t418 - t54) * t334 + (-qJD(1) * t57 + qJD(3) * t122 + t19 * t333 + t20 * t335) * t332;
t21 = t100 * t523 + t102 * t271 + t104 * t272 + t160 * t170 + t161 * t172 + t168 * t363;
t22 = t101 * t271 + t103 * t272 + t160 * t171 + t161 * t173 + t169 * t363 + t523 * t99;
t55 = t160 * t250 + t161 * t251 + t180 * t523 + t181 * t271 + t182 * t272 + t249 * t363;
t56 = t333 * t75 - t335 * t74;
t6 = (qJD(3) * t419 - t55) * t334 + (-qJD(1) * t56 + qJD(3) * t121 + t21 * t333 + t22 * t335) * t332;
t430 = t42 * t456 + t5 * t522 + t6 * t523 + (-t137 * t334 + t332 * t414) * t493 + t363 * t41;
t185 = t333 * t357 - t440 * t512;
t186 = t333 * t356 + t440 * t511;
t426 = rSges(6,1) * t186 + rSges(6,2) * t185;
t420 = -t513 - t563;
t87 = t201 * t523 + t203 * t279 + t205 * t280;
t88 = t202 * t523 + t204 * t279 + t206 * t280;
t58 = t333 * t88 - t335 * t87;
t417 = t333 * t87 + t335 * t88;
t89 = t201 * t522 + t203 * t281 + t205 * t282;
t90 = t202 * t522 + t204 * t281 + t206 * t282;
t59 = t333 * t90 - t335 * t89;
t416 = t333 * t89 + t335 * t90;
t415 = t333 * t92 - t335 * t91;
t413 = t333 * t95 + t335 * t96;
t411 = Icges(4,1) * t345 + t547;
t409 = Icges(5,1) * t332 + t545;
t405 = Icges(4,2) * t348 + t548;
t403 = Icges(5,2) * t334 + t546;
t390 = -t208 * t333 + t534;
t389 = -t207 * t333 - t208 * t335;
t382 = t449 - t503;
t222 = (-rSges(6,1) * t344 - rSges(6,2) * t347) * t483 + (rSges(6,3) * t332 + t334 * t424) * qJD(3);
t381 = -t222 + t433;
t380 = -pkin(2) - t428;
t199 = t434 * t335;
t378 = -t331 - t427;
t299 = t492 * t564;
t376 = t333 * (pkin(8) * t363 + qJD(1) * t309 - t299) + t335 * (-pkin(8) * t462 + t300 + (-t332 * t490 - t333 * t496) * pkin(4)) + t277 * t495 + t467;
t375 = t433 + t506;
t372 = t190 * t335 - t333 * t507;
t371 = qJD(3) * t411;
t370 = qJD(3) * t409;
t369 = qJD(3) * t405;
t368 = qJD(3) * t403;
t146 = t382 * t335;
t365 = t332 * t569 - t331 - t566;
t361 = t332 * t550 + t447;
t12 = qJD(1) * t418 - t19 * t335 + t20 * t333;
t13 = qJD(1) * t419 - t21 * t335 + t22 * t333;
t355 = t12 * t453 + t13 * t454 + t5 * t573 + t6 * t571 + (qJD(1) * t414 + t551 - t552) * t572 + t41 * t497 / 0.2e1 + t42 * t452 + t415 * t493 / 0.2e1 + t598 * t57 + t597 * t56;
t132 = t137 * t493;
t63 = t600 * t334 + (-t473 + (-t251 * t339 - t181) * t336) * t332 + t469;
t7 = t132 + (qJD(3) * t414 - t63) * t334 + (-qJD(1) * t415 + t25 * t333 + t26 * t335) * t332;
t353 = -t334 * t7 - t42 * t462 + t430;
t352 = t132 + (t25 + t55) * t454 + (t26 + t54) * t453 + (t122 + t92) * t598 + (t121 + t91) * t597;
t351 = t333 * t365 + t420;
t304 = t428 * qJD(3);
t287 = t427 * qJD(3);
t253 = t450 * t333;
t225 = t565 + t338 + (pkin(2) - t555) * t335 + t500;
t224 = t333 * t380 + t554 + t583;
t212 = t247 + t431;
t211 = -t563 + (rSges(5,3) - t343) * t335 + t378 * t333;
t198 = t434 * t333;
t179 = -t294 * t495 - t287 * t333 + (-t333 * t484 - t335 * t494) * pkin(3);
t178 = t294 * t497 + t314 + (-t287 - t476) * t335;
t167 = t320 * t492 + (-t338 + (-rSges(4,3) - pkin(7)) * t333 + t380 * t335) * qJD(1);
t166 = ((-pkin(2) - t557) * t333 + t583) * qJD(1) + t354;
t157 = t174 * t522;
t149 = -t259 * t334 + (-t262 * t344 + t265 * t347) * t332;
t148 = -t208 * t334 - t268 * t522;
t147 = t207 * t334 + t268 * t523;
t145 = t382 * t333;
t144 = t294 * t492 + (t335 * t378 - t327 - t338) * qJD(1) + t464;
t143 = t325 + qJD(3) * t254 + ((-t331 - t556) * t333 + t420) * qJD(1) + t501;
t142 = t431 + t208 + t278;
t141 = t351 + t425;
t140 = t149 * t493;
t139 = -t175 * t334 - t252 * t522;
t131 = t379 + t431 + t175;
t130 = -t563 + (pkin(5) * t344 - t343) * t335 + t361 * t333 + t422;
t127 = t390 * t332;
t126 = -t443 + t299 + (-t444 + (-t524 - t591) * qJD(3)) * t333 + (t335 * t360 + t322) * qJD(1);
t124 = qJD(1) * t199 + t333 * t381;
t123 = t335 * t381 + t463 + t502;
t120 = -t175 * t523 + t157;
t119 = rSges(6,3) * t363 + t426;
t118 = -rSges(6,3) * t462 + t468;
t117 = Icges(6,1) * t186 + Icges(6,4) * t185 + Icges(6,5) * t363;
t116 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t362;
t115 = Icges(6,4) * t186 + Icges(6,2) * t185 + Icges(6,6) * t363;
t114 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t362;
t113 = Icges(6,5) * t186 + Icges(6,6) * t185 + Icges(6,3) * t363;
t112 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t362;
t94 = -t332 * t448 - t334 * t507;
t93 = t190 * t334 + t239 * t523 + t138;
t82 = t299 + t569 * t457 + (t335 * t365 - t338) * qJD(1) - t426 + t464;
t81 = t300 + t325 + (-t564 - t567) * t490 + t351 * qJD(1) + t468;
t80 = qJD(1) * t146 + t333 * t375;
t79 = t335 * t375 + t497 * t503 + t502;
t78 = -t389 + t439;
t73 = t332 * t372 + t157;
t72 = t443 + (t444 + (t334 * t550 + t524) * qJD(3)) * t333 + (t335 * t361 - t322 - t338) * qJD(1) - t423 + t464;
t71 = t325 + (-t444 + (-t334 * t350 - t524 - t567) * qJD(3)) * t335 + ((t447 - t553) * t333 + t420) * qJD(1) + t465 + t470;
t70 = t599 * t334 + (-t533 + (-t262 * t347 - t265 * t344) * qJD(5)) * t332 + t466;
t69 = (t268 * t492 + t119) * t334 + (-qJD(3) * t207 + t222 * t333 + t268 * t495) * t332;
t68 = (-t268 * t490 - t118) * t334 + (qJD(3) * t208 - t222 * t335 + t463) * t332;
t66 = t439 - t580;
t65 = -t174 * t493 + t442;
t64 = -t187 * t522 + (-t252 * t490 - t105) * t334 + t509;
t62 = t185 * t262 + t186 * t265 + t213 * t523 + t216 * t279 + t219 * t280 + t259 * t363;
t61 = t183 * t262 + t184 * t265 + t213 * t522 + t216 * t281 + t219 * t282 + t259 * t362;
t46 = t390 * t491 + (qJD(1) * t389 - t118 * t333 + t119 * t335) * t332;
t45 = -t129 * t334 + t332 * t416;
t44 = -t128 * t334 + t332 * t417;
t43 = -t175 * t457 + (-t105 * t333 + (-t174 * t333 - t175 * t335) * qJD(1)) * t332 + t549;
t35 = (t239 * t492 + t126) * t334 + (-qJD(3) * t508 + t223 * t333 + t239 * t495) * t332 + t442;
t34 = (-qJD(3) * t448 - t510) * t334 + (qJD(3) * t191 + t239 * t497 + t335 * t506) * t332 + t509;
t33 = t118 * t335 + t119 * t333 + (t534 + (-t208 + t504) * t333) * qJD(1) + t376;
t32 = (qJD(3) * t391 - t112) * t334 + (qJD(3) * t202 - t114 * t344 + t116 * t347 + (-t204 * t347 - t206 * t344) * qJD(5)) * t332;
t31 = (qJD(3) * t392 - t113) * t334 + (qJD(3) * t201 - t115 * t344 + t117 * t347 + (-t203 * t347 - t205 * t344) * qJD(5)) * t332;
t30 = t112 * t523 + t114 * t279 + t116 * t280 + t185 * t204 + t186 * t206 + t202 * t363;
t29 = t113 * t523 + t115 * t279 + t117 * t280 + t185 * t203 + t186 * t205 + t201 * t363;
t28 = t112 * t522 + t114 * t281 + t116 * t282 + t183 * t204 + t184 * t206 + t202 * t362;
t27 = t113 * t522 + t115 * t281 + t117 * t282 + t183 * t203 + t184 * t205 + t201 * t362;
t18 = t372 * t491 + (qJD(1) * t580 + t126 * t335 - t510 * t333) * t332 + t549;
t17 = t510 * t335 + (t106 + t126) * t333 + (t508 * t335 + (t504 - t507) * t333) * qJD(1) + t376;
t16 = qJD(1) * t417 - t29 * t335 + t30 * t333;
t15 = qJD(1) * t416 - t27 * t335 + t28 * t333;
t9 = (qJD(3) * t417 - t62) * t334 + (-qJD(1) * t58 + qJD(3) * t128 + t29 * t333 + t30 * t335) * t332;
t8 = (qJD(3) * t416 - t61) * t334 + (-qJD(1) * t59 + qJD(3) * t129 + t27 * t333 + t28 * t335) * t332;
t1 = [t469 + t466 + (t130 * t72 + t131 * t71) * t576 + (t141 * t82 + t142 * t81) * t577 + (t143 * t212 + t144 * t211) * t578 + (t166 * t225 + t167 * t224) * t579 - t336 * t251 * t521 - t265 * t455 + (-t403 + t410) * t493 + (t404 + t409) * t491 + (-t405 + t412) * t486 + (t406 + t411) * t484 + (t599 + t600) * t334 + (-t181 * t336 - t262 * t481 - t473 - t533) * t332; 0; 0; -t552 / 0.2e1 + t551 / 0.2e1 + m(5) * (t143 * t253 + t144 * t254 + t178 * t211 + t179 * t212) + m(6) * (t123 * t141 + t124 * t142 + t198 * t81 + t199 * t82) + m(7) * (t130 * t79 + t131 * t80 + t145 * t71 + t146 * t72) + m(4) * ((-t166 * t333 - t167 * t335) * t320 + (-t224 * t335 - t225 * t333) * t304) + ((t92 / 0.2e1 + t122 / 0.2e1 - t225 * t568 + t531 / 0.2e1 + t529 / 0.2e1 + t527 / 0.2e1 + t525 / 0.2e1 + t471) * t335 + (t91 / 0.2e1 + t121 / 0.2e1 + t532 / 0.2e1 + t530 / 0.2e1 + t224 * t568 + t528 / 0.2e1 + t526 / 0.2e1 + t472) * t333) * qJD(1) + t607 * qJD(3) * (t574 / 0.2e1 + t575 / 0.2e1) + (-qJD(3) * t582 + (-qJD(1) * t242 - t335 * t368) * t334 + (-qJD(1) * t244 - t335 * t370) * t332 + (-qJD(1) * t263 - t335 * t369) * t348 + (-qJD(1) * t266 - t335 * t371) * t345 + t32 + t54 + t61) * t573 + (-qJD(3) * t601 + (qJD(1) * t243 - t333 * t368) * t334 + (qJD(1) * t245 - t333 * t370) * t332 + (qJD(1) * t264 - t333 * t369) * t348 + (qJD(1) * t267 - t333 * t371) * t345 + t31 + t55 + t62) * t571; m(4) * t107 + m(5) * t67 + m(6) * t33 + m(7) * t17; (t145 * t80 + t146 * t79 + t66 * t17) * t576 + (t123 * t199 + t124 * t198 + t33 * t78) * t577 + (t254 * t178 + t253 * t179 + t505 * t67) * t578 + (t574 + t575) * t320 * t304 * t579 + (t56 + t58) * t497 + (t57 + t59) * t495 + (t247 * t594 + t270 * t593 + t590 * t497 + t585 * t574 - t13 - t16 + (-t335 * t601 - t595) * t495) * t335 + (t12 + t15 + t246 * t594 + t269 * t593 + t586 * t575 + t587 * t495 + ((t243 * t491 + t245 * t493 + t264 * t484 + t267 * t486 + t585 - t596) * t333 + (t242 * t491 + t244 * t493 + t263 * t484 + t266 * t486 + t586) * t335 + ((-t525 - t527 - t529 - t531) * t333 + (-t526 - t528 - t530 - t532) * t335) * qJD(3) + ((-t601 + t603) * t333 + t605 + t587 + t590) * qJD(1)) * t335 + (t595 + t606) * t497) * t333; m(7) * (t333 * t72 - t335 * t71 + (t130 * t335 + t131 * t333) * qJD(1)) + m(6) * (t333 * t82 - t335 * t81 + (t141 * t335 + t142 * t333) * qJD(1)) + m(5) * (-t143 * t335 + t144 * t333 + (t211 * t335 + t212 * t333) * qJD(1)); 0; m(7) * (t333 * t79 - t335 * t80 + (t145 * t333 + t146 * t335) * qJD(1)) + m(6) * (t123 * t333 - t124 * t335 + (t198 * t333 + t199 * t335) * qJD(1)) + m(5) * (t178 * t333 - t179 * t335 + (t253 * t333 + t254 * t335) * qJD(1)); 0; (-t70 - t63 + (t333 * t472 + t335 * t471) * qJD(3)) * t334 + t352 + m(6) * (t141 * t69 + t142 * t68 + t147 * t82 + t148 * t81) + m(7) * (t130 * t35 + t131 * t34 + t71 * t94 + t72 * t93) + ((t32 / 0.2e1 + t61 / 0.2e1) * t335 + (t31 / 0.2e1 + t62 / 0.2e1) * t333 + (-t333 * t471 + t335 * t472) * qJD(1)) * t332 + t140; m(6) * t46 + m(7) * t18; t355 + (t58 * t451 + (qJD(1) * t95 + t32) * t572 + t8 / 0.2e1 + t44 * t562) * t333 + (t59 * t451 + (qJD(1) * t96 - t31) * t572 - t9 / 0.2e1 + t45 * t562) * t335 + (t16 * t573 + t15 * t570 + qJD(3) * (t333 * t96 - t335 * t95) / 0.2e1 + (t58 * t570 - t333 * t59 / 0.2e1) * qJD(1)) * t332 + m(6) * (t123 * t147 + t124 * t148 + t127 * t33 + t198 * t68 + t199 * t69 + t46 * t78) + m(7) * (t145 * t34 + t146 * t35 + t17 * t73 + t18 * t66 + t79 * t93 + t80 * t94); m(6) * (t333 * t69 - t335 * t68 + (t147 * t335 + t148 * t333) * qJD(1)) + m(7) * (t333 * t35 - t335 * t34 + (t333 * t94 + t335 * t93) * qJD(1)); (t18 * t73 + t34 * t94 + t35 * t93) * t576 + (t127 * t46 + t147 * t69 + t148 * t68) * t577 + (t70 * t334 - t140 - t7 + (t333 * t44 - t334 * t413 + t335 * t45) * qJD(3)) * t334 + (t335 * t8 + t333 * t9 - t334 * (t31 * t333 + t32 * t335) + (-t149 * t334 + t332 * t413) * qJD(3) + ((-t334 * t95 + t44) * t335 + (t334 * t96 - t42 - t45) * t333) * qJD(1)) * t332 + t430; m(7) * (t130 * t65 + t131 * t64 + t138 * t72 + t139 * t71) + t352 - t63 * t334; m(7) * t43; t355 + m(7) * (t120 * t17 + t138 * t79 + t139 * t80 + t145 * t64 + t146 * t65 + t43 * t66); m(7) * (t333 * t65 - t335 * t64 + (t138 * t335 + t139 * t333) * qJD(1)); m(7) * (t120 * t18 + t138 * t35 + t139 * t34 + t43 * t73 + t64 * t94 + t65 * t93) + t353; (t120 * t43 + t138 * t65 + t139 * t64) * t576 + t353;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
