% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:23
% EndTime: 2019-12-05 18:01:43
% DurationCPUTime: 14.06s
% Computational Cost: add. (13128->495), mult. (9630->586), div. (0->0), fcn. (7393->8), ass. (0->287)
t564 = Icges(5,3) + Icges(6,3);
t274 = qJ(1) + pkin(8);
t266 = qJ(3) + t274;
t259 = sin(t266);
t260 = cos(t266);
t278 = cos(qJ(4));
t404 = t260 * t278;
t276 = sin(qJ(4));
t405 = t260 * t276;
t138 = Icges(6,4) * t404 - Icges(6,2) * t405 + Icges(6,6) * t259;
t140 = Icges(5,4) * t404 - Icges(5,2) * t405 + Icges(5,6) * t259;
t557 = t138 + t140;
t223 = Icges(6,4) * t405;
t142 = Icges(6,1) * t404 + Icges(6,5) * t259 - t223;
t225 = Icges(5,4) * t405;
t144 = Icges(5,1) * t404 + Icges(5,5) * t259 - t225;
t555 = t142 + t144;
t234 = Icges(6,5) * t278 - Icges(6,6) * t276;
t236 = Icges(5,5) * t278 - Icges(5,6) * t276;
t567 = t234 + t236;
t566 = Icges(5,5) + Icges(6,5);
t565 = -Icges(5,6) - Icges(6,6);
t429 = Icges(6,4) * t276;
t237 = Icges(6,2) * t278 + t429;
t430 = Icges(5,4) * t276;
t239 = Icges(5,2) * t278 + t430;
t561 = t237 + t239;
t242 = Icges(6,1) * t278 - t429;
t244 = Icges(5,1) * t278 - t430;
t563 = t242 + t244;
t268 = Icges(5,4) * t278;
t483 = Icges(5,1) * t276 + t268;
t267 = Icges(6,4) * t278;
t484 = Icges(6,1) * t276 + t267;
t560 = t483 + t484;
t562 = -t567 * t259 + t564 * t260;
t329 = -Icges(6,2) * t276 + t267;
t309 = t329 * t259;
t137 = Icges(6,6) * t260 - t309;
t330 = -Icges(5,2) * t276 + t268;
t139 = Icges(5,6) * t260 - t259 * t330;
t558 = t137 + t139;
t408 = t259 * t276;
t222 = Icges(6,4) * t408;
t407 = t259 * t278;
t141 = -Icges(6,1) * t407 + Icges(6,5) * t260 + t222;
t224 = Icges(5,4) * t408;
t143 = -Icges(5,1) * t407 + Icges(5,5) * t260 + t224;
t556 = t141 + t143;
t540 = t557 * t276 - t555 * t278;
t559 = t564 * t259 + t566 * t404 + t565 * t405;
t233 = Icges(6,5) * t276 + Icges(6,6) * t278;
t235 = Icges(5,5) * t276 + Icges(5,6) * t278;
t554 = t233 + t235;
t273 = qJD(1) + qJD(3);
t553 = t563 * t273;
t552 = t560 * qJD(4) - t566 * t273;
t551 = t561 * qJD(4) + t565 * t273;
t542 = t561 * t276 - t560 * t278;
t495 = -t562 * t259 - t556 * t404;
t494 = t562 * t260 + t558 * t408;
t550 = t556 * t278;
t493 = t558 * t276;
t549 = t540 * t260;
t508 = -t556 * t407 + t494;
t507 = t559 * t260 - t555 * t407 + t557 * t408;
t506 = -t558 * t405 - t495;
t505 = t559 * t259 - t549;
t504 = t556 * t276 + t558 * t278;
t503 = t555 * t276 + t557 * t278;
t310 = t330 * t273;
t548 = -t259 * t310 - t551 * t260 - t273 * t309;
t406 = t260 * t273;
t547 = t551 * t259 - t260 * t310 - t329 * t406;
t546 = -t553 * t259 - t552 * t260;
t545 = t552 * t259 - t553 * t260;
t544 = (t329 + t330) * qJD(4);
t543 = t563 * qJD(4);
t541 = t560 * t276 + t561 * t278;
t539 = t493 - t550;
t411 = t236 * t273;
t413 = t234 * t273;
t538 = -t411 - t413;
t537 = t554 * qJD(4) - t564 * t273;
t281 = qJD(1) ^ 2;
t412 = t235 * t260;
t414 = t233 * t260;
t500 = t259 * t542 + t412 + t414;
t155 = t233 * t259;
t157 = t235 * t259;
t499 = -t260 * t542 + t155 + t157;
t536 = t567 * qJD(4) + t542 * t273;
t535 = t537 * t259 + t538 * t260 + t539 * t273;
t534 = t538 * t259 - t537 * t260 + t540 * t273;
t271 = t278 * pkin(4);
t261 = t271 + pkin(3);
t201 = t260 * t261;
t372 = rSges(6,1) * t404;
t317 = rSges(6,3) * t259 + t372;
t533 = -t317 - t201;
t532 = t505 * t259 + t506 * t260;
t531 = t507 * t259 + t508 * t260;
t530 = t541 * qJD(4) - t554 * t273 + t544 * t276 - t543 * t278;
t529 = t504 * qJD(4) - t562 * t273 + t547 * t276 - t545 * t278;
t528 = t503 * qJD(4) - t273 * t559 + t548 * t276 - t546 * t278;
t527 = t500 * t273;
t230 = rSges(6,2) * t405;
t275 = -qJ(5) - pkin(7);
t443 = pkin(7) + t275;
t447 = pkin(3) * t260;
t395 = t259 * t443 + t230 + t447 + t533;
t526 = t273 * t395;
t337 = rSges(4,1) * t259 + rSges(4,2) * t260;
t150 = t337 * t273;
t264 = sin(t274);
t448 = pkin(2) * t264;
t277 = sin(qJ(1));
t450 = pkin(1) * t277;
t343 = t448 + t450;
t315 = t343 * qJD(1);
t125 = t315 + t150;
t487 = rSges(6,2) * t408 + t260 * rSges(6,3);
t524 = t443 * t260 - t487;
t523 = t499 * t273;
t522 = t259 * t535 - t260 * t529;
t521 = t259 * t534 - t260 * t528;
t409 = t259 * t273;
t217 = pkin(3) * t409;
t152 = pkin(7) * t406 - t217;
t375 = qJD(4) * t273;
t171 = qJDD(4) * t259 + t260 * t375;
t270 = t278 * rSges(6,1);
t248 = -rSges(6,2) * t276 + t270;
t202 = t248 * qJD(4);
t433 = t278 * rSges(6,2);
t245 = rSges(6,1) * t276 + t433;
t251 = qJD(5) * t259;
t272 = qJDD(1) + qJDD(3);
t263 = t281 * t450;
t265 = cos(t274);
t279 = cos(qJ(1));
t449 = pkin(1) * t279;
t342 = pkin(2) * t265 + t449;
t291 = -qJDD(1) * t342 + t281 * t448 + t263;
t446 = pkin(7) * t259;
t183 = t446 + t447;
t369 = -t183 + t395;
t377 = qJD(4) * t259;
t401 = t278 * qJD(4) ^ 2;
t374 = qJD(4) * t276;
t373 = qJD(4) * t278;
t358 = t260 * t373;
t402 = t273 * t275;
t359 = t260 * t374;
t497 = -t273 * t407 - t359;
t482 = rSges(6,1) * t497 - rSges(6,2) * t358 - t260 * t402 - t261 * t409 + t251;
t442 = t217 + (-pkin(4) * t374 - pkin(7) * t273) * t260 + t273 * t487 + t482;
t15 = t202 * t377 + qJDD(5) * t260 + t171 * t245 + (t171 * t276 + t259 * t401) * pkin(4) + t369 * t272 + (-t152 - t251 - t442) * t273 + t291;
t492 = -g(2) + t15;
t172 = qJDD(4) * t260 - t259 * t375;
t400 = t279 * t281;
t287 = (-qJDD(1) * t264 - t265 * t281) * pkin(2) + (-qJDD(1) * t277 - t400) * pkin(1);
t445 = t260 * pkin(7);
t341 = -pkin(3) * t259 + t445;
t286 = -t273 ^ 2 * t183 + t272 * t341 + t287;
t347 = -t261 - t270;
t344 = -pkin(3) - t347;
t348 = pkin(4) * t276 + t245;
t444 = pkin(3) - t261;
t361 = t259 * t374;
t384 = pkin(4) * t361 + qJD(5) * t260;
t496 = t259 * t373 + t273 * t405;
t481 = rSges(6,1) * t361 + rSges(6,2) * t496 + t259 * t402 + t384;
t441 = t481 + (t260 * t444 - t317 + t446) * t273;
t14 = t272 * t487 + t441 * t273 - t348 * t172 + (-t272 * t344 + qJDD(5)) * t259 + (-pkin(4) * t401 - qJD(4) * t202 + qJD(5) * t273 - t272 * t443) * t260 + t286;
t491 = -g(3) + t14;
t216 = rSges(4,2) * t409;
t151 = -rSges(4,1) * t406 + t216;
t520 = t151 * t273 - t272 * t337 - g(3) + t287;
t437 = rSges(5,2) * t276;
t438 = rSges(5,1) * t278;
t249 = -t437 + t438;
t203 = t249 * qJD(4);
t246 = rSges(5,1) * t276 + rSges(5,2) * t278;
t382 = rSges(5,2) * t408 + t260 * rSges(5,3);
t320 = rSges(5,1) * t407 - t382;
t376 = qJD(4) * t260;
t318 = -rSges(5,1) * t404 - rSges(5,3) * t259;
t364 = rSges(5,1) * t361 + rSges(5,2) * t496;
t95 = t273 * t318 + t364;
t26 = -t172 * t246 - t203 * t376 - t272 * t320 + t273 * t95 + t286;
t519 = t26 - g(3);
t231 = rSges(5,2) * t405;
t146 = -t231 - t318;
t386 = -t146 - t183;
t366 = rSges(5,1) * t497 - rSges(5,2) * t358;
t93 = t273 * t382 + t366;
t27 = t203 * t377 + t171 * t246 + (-t152 - t93) * t273 + t386 * t272 + t291;
t518 = t27 - g(2);
t439 = rSges(4,1) * t260;
t182 = -t259 * rSges(4,2) + t439;
t517 = t150 * t273 - t182 * t272 - g(2) + t291;
t516 = t259 * t529 + t260 * t535;
t515 = t259 * t528 + t260 * t534;
t514 = qJD(4) * t531 + t527;
t513 = qJD(4) * t532 + t523;
t512 = -t539 * qJD(4) + t545 * t276 + t547 * t278;
t511 = -t540 * qJD(4) + t546 * t276 + t548 * t278;
t510 = t259 * t536 - t530 * t260;
t509 = t530 * t259 + t260 * t536;
t380 = t484 + t329;
t381 = t237 - t242;
t502 = (t276 * t380 + t278 * t381) * t273;
t378 = t483 + t330;
t379 = t239 - t244;
t501 = (t276 * t378 + t278 * t379) * t273;
t498 = -t559 - t550;
t288 = t259 * t344 + t524;
t479 = t342 * qJD(1);
t130 = t273 * t146;
t180 = t246 * t377;
t488 = t180 - t130;
t167 = rSges(6,1) * t408 + rSges(6,2) * t407;
t232 = pkin(4) * t408;
t486 = t232 + t167;
t485 = t248 + t271;
t258 = t264 * rSges(3,2);
t440 = rSges(3,1) * t265;
t340 = -t440 - t449;
t178 = t258 + t340;
t338 = rSges(3,1) * t264 + rSges(3,2) * t265;
t306 = t338 + t450;
t128 = t273 * t320;
t175 = t273 * t341;
t319 = t246 * t376 + t128 - t175;
t63 = t315 + t319;
t64 = t273 * t386 + t180 - t479;
t480 = t259 * t64 + t260 * t63;
t368 = t245 * t377 + t384;
t478 = t368 + t526;
t387 = -Icges(5,2) * t404 + t144 - t225;
t391 = t260 * t483 + t140;
t463 = t276 * t387 + t278 * t391;
t388 = Icges(5,2) * t407 + t143 + t224;
t392 = -t259 * t483 + t139;
t462 = -t276 * t388 - t278 * t392;
t389 = -Icges(6,2) * t404 + t142 - t223;
t393 = t260 * t484 + t138;
t461 = t276 * t389 + t278 * t393;
t390 = Icges(6,2) * t407 + t141 + t222;
t394 = -t259 * t484 + t137;
t460 = -t276 * t390 - t278 * t394;
t459 = m(3) + m(4);
t458 = t171 / 0.2e1;
t457 = t172 / 0.2e1;
t451 = -rSges(5,3) - pkin(7);
t168 = t246 * t259;
t410 = t246 * t260;
t403 = t273 * t182;
t355 = -pkin(3) - t438;
t353 = -t377 / 0.2e1;
t352 = t377 / 0.2e1;
t351 = -t376 / 0.2e1;
t350 = t376 / 0.2e1;
t250 = rSges(2,1) * t279 - t277 * rSges(2,2);
t339 = rSges(2,1) * t277 + rSges(2,2) * t279;
t176 = t273 * t183;
t299 = t176 + t479;
t107 = -t372 - t201 + t230 + (-rSges(6,3) + t275) * t259;
t108 = t259 * t355 + t382 + t445;
t109 = t259 * t451 + t260 * t355 + t231;
t290 = t260 * t146 + t259 * t320;
t106 = t259 * t347 - t260 * t275 + t487;
t289 = t348 * t376 - t175 - t251 + (rSges(6,1) * t407 - t259 * t444 + t524) * t273;
t285 = t259 * t288 - t260 * t395;
t284 = (((t494 + t505 + t549) * t260 + ((-t493 + t498) * t260 - t495 - t506 + t507) * t259) * qJD(4) + t527) * t353 + (-t542 * qJD(4) + t543 * t276 + t544 * t278) * t273 + (Icges(4,3) + t541) * t272 + (t499 + t503) * t458 + (t500 + t504) * t457 + ((((t493 - t559) * t260 + t495 + t507) * t260 + (t259 * t498 + t494 - t508) * t259) * qJD(4) + t513 - t523) * t351 + (t509 + t512) * t350 + (t510 + t511 + t514) * t352;
t283 = t64 * (t217 - t366) + ((-t437 * t64 - t451 * t63) * t259 + (-t355 * t63 + t451 * t64) * t260) * t273 - t63 * t364;
t49 = t315 + t289;
t50 = t273 * t369 + t368 - t479;
t282 = t50 * (pkin(4) * t359 - t482) + (-t50 * t487 - t49 * t533) * t273 - t49 * t481;
t169 = t245 * t260;
t126 = -t479 - t403;
t69 = qJD(4) * t290 + qJD(2);
t44 = qJD(4) * t285 + qJD(2);
t25 = qJDD(2) + t172 * t146 + t171 * t320 + (-t259 * t95 + t260 * t93) * qJD(4);
t5 = qJDD(2) - t395 * t172 + t288 * t171 + (-t441 * t259 + t442 * t260) * qJD(4);
t1 = [t284 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + ((qJDD(1) * t178 + t281 * t338 - g(2) + t263) * t178 + (t400 * pkin(1) + t306 * qJDD(1) + g(3) + (-0.2e1 * t258 + t440 - t340 + t178) * t281) * t306) * m(3) + ((qJDD(1) * t339 + g(3)) * t339 + (qJDD(1) * t250 + g(2)) * t250) * m(2) + (-(t50 + t299 - t478) * t49 + (t342 * t49 + t343 * t50) * qJD(1) + t282 + t492 * (t107 - t342) + t491 * (t106 - t343)) * m(6) + (-(t64 + t299 - t488) * t63 + (t342 * t63 + t343 * t64) * qJD(1) + t283 + t518 * (t109 - t342) + t519 * (t108 - t343)) * m(5) + (t517 * (-t182 - t342) + t520 * (-t337 - t343) + (t273 * t439 - t216 - t403) * t125) * m(4); t459 * qJDD(2) + m(5) * t25 + m(6) * t5 + (-m(5) - m(6) - t459) * g(1); t284 + (-t50 * t289 + t49 * (-t176 + t478) + t282 + t492 * t107 + t491 * t106) * m(6) + (-t64 * t319 + t63 * (-t176 + t488) + t283 + t518 * t109 + t519 * t108) * m(5) + (-t125 * t151 + t126 * t150 + (-t125 * t273 - t517) * t182 - (t126 * t273 + t520) * t337) * m(4); t532 * t458 + t531 * t457 + (t510 * t273 + t499 * t272 + t506 * t172 + t505 * t171 + (t521 * t259 + t522 * t260) * qJD(4)) * t259 / 0.2e1 + (t509 * t273 + t500 * t272 + t508 * t172 + t507 * t171 + (t515 * t259 + t516 * t260) * qJD(4)) * t260 / 0.2e1 + (t259 * t503 + t504 * t260) * t272 / 0.2e1 - (((t378 + t380) * t278 + (-t379 - t381) * t276) * t273 + (((t388 + t390) * t260 + (t387 + t389) * t259) * t278 + ((-t392 - t394) * t260 + (-t391 - t393) * t259) * t276) * qJD(4)) * t273 / 0.2e1 + ((t273 * t503 + t512) * t260 + (-t273 * t504 + t511) * t259) * t273 / 0.2e1 - t514 * t409 / 0.2e1 + t513 * t406 / 0.2e1 + ((-t377 * t414 + t413) * t259 + (-t502 + (t460 * t260 + (t155 - t461) * t259) * qJD(4)) * t260 + (-t377 * t412 + t411) * t259 + (-t501 + (t462 * t260 + (t157 - t463) * t259) * qJD(4)) * t260) * t353 + ((t273 * t505 + t522) * t260 + (-t273 * t506 + t521) * t259) * t352 + ((t155 * t376 + t413) * t260 + (t502 + (t461 * t259 + (-t414 - t460) * t260) * qJD(4)) * t259 + (t157 * t376 + t411) * t260 + (t501 + (t463 * t259 + (-t412 - t462) * t260) * qJD(4)) * t259) * t351 + ((t273 * t507 + t516) * t260 + (-t508 * t273 + t515) * t259) * t350 + (-g(1) * t485 - g(2) * t486 + t15 * t232 + t5 * t285 + (t167 * t273 - t245 * t409) * t49 + (-t169 * t273 + t245 * t406) * t50 + (t44 * (-t441 + t526) + t15 * t245 + t50 * t202 - (t50 * t248 - t44 * t486) * qJD(4)) * t259 + (t44 * (t288 * t273 + t442) - t14 * t348 - t49 * (-pkin(4) * t373 - t202) - g(3) * (-t433 + (-rSges(6,1) - pkin(4)) * t276) - (t49 * t485 + (-pkin(4) * t405 - t169) * t44) * qJD(4)) * t260) * m(6) + (t25 * t290 + t69 * ((-t95 - t130) * t259 + (t93 + t128) * t260) + t27 * t168 - t26 * t410 + (t64 * t406 - t63 * t409) * t246 + t480 * t203 - (-t168 * t63 + t410 * t64) * t273 - (t69 * (-t168 * t259 - t260 * t410) + t480 * t249) * qJD(4) - g(1) * t249 - g(2) * t168 + g(3) * t410) * m(5); (t259 * t491 + t260 * t492) * m(6);];
tau = t1;
