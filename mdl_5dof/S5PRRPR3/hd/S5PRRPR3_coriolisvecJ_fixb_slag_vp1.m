% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:35
% DurationCPUTime: 23.36s
% Computational Cost: add. (16536->680), mult. (13136->887), div. (0->0), fcn. (10228->8), ass. (0->395)
t642 = Icges(4,3) + Icges(5,3);
t334 = qJ(3) + pkin(9);
t325 = sin(t334);
t327 = cos(t334);
t239 = Icges(5,5) * t327 - Icges(5,6) * t325;
t336 = sin(qJ(3));
t337 = cos(qJ(3));
t276 = Icges(4,5) * t337 - Icges(4,6) * t336;
t635 = t239 + t276;
t332 = pkin(8) + qJ(2);
t326 = cos(t332);
t641 = t642 * t326;
t324 = sin(t332);
t518 = t324 * t337;
t519 = t324 * t336;
t520 = t324 * t327;
t521 = t324 * t325;
t627 = -Icges(4,5) * t518 - Icges(5,5) * t520 + Icges(4,6) * t519 + Icges(5,6) * t521 + t641;
t636 = t642 * t324 + t635 * t326;
t545 = Icges(5,6) * t326;
t170 = Icges(5,4) * t520 - Icges(5,2) * t521 - t545;
t546 = Icges(4,6) * t326;
t187 = Icges(4,4) * t518 - Icges(4,2) * t519 - t546;
t640 = t170 * t325 + t187 * t336;
t266 = Icges(5,4) * t521;
t551 = Icges(5,5) * t326;
t172 = Icges(5,1) * t520 - t266 - t551;
t286 = Icges(4,4) * t519;
t552 = Icges(4,5) * t326;
t189 = Icges(4,1) * t518 - t286 - t552;
t623 = -t172 * t327 - t189 * t337 + t640;
t620 = -t623 * t324 + t627 * t326;
t554 = Icges(5,4) * t325;
t243 = Icges(5,1) * t327 - t554;
t173 = Icges(5,5) * t324 + t243 * t326;
t555 = Icges(4,4) * t336;
t280 = Icges(4,1) * t337 - t555;
t190 = Icges(4,5) * t324 + t280 * t326;
t638 = -t173 * t520 - t190 * t518;
t238 = Icges(5,5) * t325 + Icges(5,6) * t327;
t275 = Icges(4,5) * t336 + Icges(4,6) * t337;
t637 = t275 + t238;
t240 = Icges(5,2) * t327 + t554;
t309 = Icges(5,4) * t327;
t242 = Icges(5,1) * t325 + t309;
t277 = Icges(4,2) * t337 + t555;
t329 = Icges(4,4) * t337;
t279 = Icges(4,1) * t336 + t329;
t634 = t240 * t325 - t242 * t327 + t277 * t336 - t279 * t337;
t633 = t326 * t636 + t638;
t510 = t326 * t337;
t512 = t326 * t327;
t592 = -t173 * t512 - t190 * t510 - t324 * t636;
t632 = -t172 * t512 - t189 * t510 + t324 * t627;
t397 = -Icges(5,2) * t325 + t309;
t171 = Icges(5,6) * t324 + t326 * t397;
t398 = -Icges(4,2) * t336 + t329;
t188 = Icges(4,6) * t324 + t326 * t398;
t630 = t171 * t325 + t188 * t336;
t619 = -t171 * t521 - t188 * t519 - t633;
t511 = t326 * t336;
t517 = t325 * t326;
t618 = -t170 * t517 - t187 * t511 - t632;
t617 = -t171 * t517 - t188 * t511 - t592;
t597 = t170 * t327 + t172 * t325 + t187 * t337 + t189 * t336;
t596 = t171 * t327 + t173 * t325 + t188 * t337 + t190 * t336;
t526 = t275 * t326;
t528 = t238 * t326;
t629 = -t324 * t634 - t526 - t528;
t527 = t275 * t324;
t529 = t238 * t324;
t628 = -t326 * t634 + t527 + t529;
t626 = t637 * qJD(3);
t625 = t173 * t327 + t190 * t337 - t630;
t222 = t397 * qJD(3);
t223 = t243 * qJD(3);
t259 = t398 * qJD(3);
t260 = t280 * qJD(3);
t624 = -t222 * t325 + t223 * t327 - t259 * t336 + t260 * t337 + (-t240 * t327 - t242 * t325 - t277 * t337 - t279 * t336) * qJD(3) + t637 * qJD(2);
t622 = t636 * qJD(2);
t621 = t634 * qJD(2) + t635 * qJD(3);
t616 = rSges(4,2) * t336;
t615 = t628 * qJD(2);
t328 = qJ(5) + t334;
t320 = cos(t328);
t523 = t320 * t324;
t319 = sin(t328);
t525 = t319 * t324;
t544 = Icges(6,6) * t326;
t159 = Icges(6,4) * t523 - Icges(6,2) * t525 - t544;
t296 = Icges(6,4) * t320;
t230 = Icges(6,1) * t319 + t296;
t614 = -t230 * t324 - t159;
t396 = -Icges(6,2) * t319 + t296;
t160 = Icges(6,6) * t324 + t326 * t396;
t613 = -t230 * t326 - t160;
t612 = t230 + t396;
t369 = qJD(3) * t240;
t103 = -t326 * t369 + (-t324 * t397 + t545) * qJD(2);
t371 = qJD(3) * t242;
t105 = -t326 * t371 + (-t243 * t324 + t551) * qJD(2);
t370 = qJD(3) * t277;
t119 = -t326 * t370 + (-t324 * t398 + t546) * qJD(2);
t372 = qJD(3) * t279;
t121 = -t326 * t372 + (-t280 * t324 + t552) * qJD(2);
t611 = -t596 * qJD(3) - t103 * t325 + t105 * t327 - t119 * t336 + t121 * t337 + t622;
t104 = qJD(2) * t171 - t324 * t369;
t106 = qJD(2) * t173 - t324 * t371;
t120 = qJD(2) * t188 - t324 * t370;
t122 = qJD(2) * t190 - t324 * t372;
t610 = t627 * qJD(2) + t597 * qJD(3) + t104 * t325 - t106 * t327 + t120 * t336 - t122 * t337;
t609 = (t617 * t324 - t618 * t326) * qJD(3);
t608 = (t619 * t324 - t620 * t326) * qJD(3);
t607 = t629 * qJD(2);
t606 = qJD(2) * t623 - t324 * t626 + t622;
t605 = -t626 * t326 + (-t635 * t324 - t625 + t641) * qJD(2);
t604 = 0.2e1 * qJD(3);
t603 = t607 + t608;
t602 = t609 + t615;
t601 = t623 * qJD(3) - t104 * t327 - t106 * t325 - t120 * t337 - t122 * t336;
t600 = qJD(3) * t625 + t103 * t327 + t105 * t325 + t119 * t337 + t121 * t336;
t599 = t621 * t324 + t624 * t326;
t598 = t624 * t324 - t621 * t326;
t317 = t326 * pkin(6);
t250 = pkin(2) * t324 - t317;
t335 = -qJ(4) - pkin(6);
t292 = t326 * t335;
t330 = t337 * pkin(3);
t321 = t330 + pkin(2);
t474 = t324 * t321 + t292;
t155 = t250 - t474;
t151 = qJD(2) * t155;
t235 = qJD(2) * t250;
t595 = t151 - t235;
t298 = qJD(4) * t324;
t570 = pkin(4) * t325;
t571 = pkin(3) * t336;
t269 = -t570 - t571;
t244 = t269 * qJD(3);
t513 = t326 * t244;
t594 = t298 + t513;
t316 = t324 * pkin(6);
t251 = t326 * pkin(2) + t316;
t593 = t627 + t630;
t553 = Icges(6,4) * t319;
t231 = Icges(6,1) * t320 - t553;
t162 = Icges(6,5) * t324 + t231 * t326;
t227 = Icges(6,5) * t320 - Icges(6,6) * t319;
t333 = qJD(3) + qJD(5);
t226 = Icges(6,5) * t319 + Icges(6,6) * t320;
t531 = t226 * t326;
t537 = t160 * t319;
t541 = Icges(6,3) * t326;
t591 = -t333 * t531 + (-t162 * t320 - t227 * t324 + t537 + t541) * qJD(2);
t255 = Icges(6,4) * t525;
t550 = Icges(6,5) * t326;
t161 = Icges(6,1) * t523 - t255 - t550;
t394 = t159 * t319 - t161 * t320;
t158 = Icges(6,3) * t324 + t227 * t326;
t467 = qJD(2) * t158;
t532 = t226 * t324;
t590 = qJD(2) * t394 - t333 * t532 + t467;
t228 = Icges(6,2) * t320 + t553;
t387 = t228 * t319 - t230 * t320;
t585 = qJD(2) * t387 + t227 * t333;
t488 = -t240 * t326 + t173;
t489 = -Icges(5,2) * t520 + t172 - t266;
t582 = t324 * t488 - t326 * t489;
t485 = -Icges(4,2) * t518 + t189 - t286;
t487 = t279 * t324 + t187;
t581 = -t336 * t485 - t337 * t487;
t236 = t324 * t333;
t237 = t326 * t333;
t580 = qJD(2) * t612 + t236 * (-t228 * t326 + t162) - t237 * (-Icges(6,2) * t523 + t161 - t255);
t410 = qJD(2) * t333;
t219 = t324 * t410;
t579 = t219 / 0.2e1;
t220 = t326 * t410;
t578 = t220 / 0.2e1;
t577 = -t236 / 0.2e1;
t576 = t236 / 0.2e1;
t575 = -t237 / 0.2e1;
t574 = t237 / 0.2e1;
t573 = t324 / 0.2e1;
t572 = -t326 / 0.2e1;
t569 = pkin(4) * t327;
t568 = -qJD(2) / 0.2e1;
t567 = qJD(2) / 0.2e1;
t566 = pkin(2) - t321;
t565 = rSges(4,1) * t337;
t564 = rSges(5,1) * t327;
t563 = rSges(6,1) * t319;
t562 = rSges(6,1) * t320;
t561 = rSges(6,2) * t320;
t312 = t324 * rSges(4,3);
t311 = t324 * rSges(5,3);
t310 = t324 * rSges(6,3);
t247 = rSges(5,1) * t325 + rSges(5,2) * t327;
t267 = rSges(5,2) * t521;
t174 = rSges(5,1) * t520 - t326 * rSges(5,3) - t267;
t429 = -t247 - t571;
t459 = qJD(3) * t326;
t383 = t429 * t459;
t365 = t298 + t383;
t495 = t155 - t250;
t62 = (-t174 + t495) * qJD(2) + t365;
t560 = t62 * t247;
t469 = rSges(4,2) * t519 + t326 * rSges(4,3);
t193 = rSges(4,1) * t518 - t469;
t282 = rSges(4,1) * t336 + rSges(4,2) * t337;
t440 = t282 * t459;
t107 = -t440 + (-t193 - t250) * qJD(2);
t540 = t107 * t324;
t539 = t107 * t326;
t460 = qJD(3) * t324;
t441 = t282 * t460;
t194 = rSges(4,1) * t510 - rSges(4,2) * t511 + t312;
t483 = t194 + t251;
t108 = qJD(2) * t483 - t441;
t215 = t282 * t326;
t538 = t108 * t215;
t530 = t228 * t333;
t524 = t319 * t326;
t522 = t320 * t326;
t157 = Icges(6,5) * t523 - Icges(6,6) * t525 - t541;
t516 = t326 * t157;
t79 = -t324 * t387 - t531;
t509 = t79 * qJD(2);
t462 = qJD(2) * t324;
t281 = t335 * t462;
t457 = qJD(3) * t336;
t273 = t324 * pkin(3) * t457;
t299 = qJD(4) * t326;
t472 = t273 + t299;
t444 = t281 + t472;
t115 = (-t326 * t566 - t316) * qJD(2) - t444;
t232 = t251 * qJD(2);
t507 = -t115 - t232;
t261 = t321 + t569;
t331 = -pkin(7) + t335;
t482 = -t324 * t261 - t326 * t331;
t127 = t474 + t482;
t256 = rSges(6,2) * t525;
t477 = t326 * rSges(6,3) + t256;
t163 = rSges(6,1) * t523 - t477;
t506 = t127 - t163;
t218 = t326 * t261;
t271 = t326 * t321;
t468 = -t331 + t335;
t128 = t324 * t468 + t218 - t271;
t413 = -t324 * t335 + t271;
t156 = t413 - t251;
t505 = -t128 - t156;
t504 = -t324 * t157 - t161 * t522;
t503 = t324 * t158 + t162 * t522;
t499 = -t324 * t155 + t326 * t156;
t164 = rSges(6,1) * t522 - rSges(6,2) * t524 + t310;
t498 = t324 * t163 + t326 * t164;
t175 = rSges(5,1) * t512 - rSges(5,2) * t517 + t311;
t494 = -t156 - t175;
t491 = -t242 * t324 - t170;
t490 = -t242 * t326 - t171;
t486 = -t279 * t326 - t188;
t484 = -t277 * t326 + t190;
t480 = -t240 + t243;
t479 = t242 + t397;
t461 = qJD(2) * t326;
t478 = rSges(5,3) * t461 + qJD(2) * t267;
t476 = t261 - t321;
t454 = qJD(2) * qJD(4);
t475 = qJD(2) * t273 + t326 * t454;
t473 = rSges(4,3) * t461 + t462 * t616;
t471 = -t277 + t280;
t470 = t279 + t398;
t464 = qJD(2) * t239;
t463 = qJD(2) * t276;
t458 = qJD(3) * t327;
t456 = qJD(3) * t337;
t453 = pkin(3) * t511;
t450 = t333 * t561;
t380 = rSges(6,3) * t461 + qJD(2) * t256 - t326 * t450;
t447 = t319 * t237;
t94 = (-t320 * t462 - t447) * rSges(6,1) + t380;
t233 = t561 + t563;
t191 = t233 * t324;
t234 = -rSges(6,2) * t319 + t562;
t95 = -t333 * t191 + (t234 * t326 + t310) * qJD(2);
t452 = t163 * t461 + t324 * t95 + t326 * t94;
t297 = pkin(6) * t461;
t438 = t326 * t457;
t412 = pkin(3) * t438;
t114 = -t412 - t297 + t298 + (t324 * t566 - t292) * qJD(2);
t451 = t115 * t460 + (t114 - t151) * t459;
t449 = pkin(3) * t456;
t448 = t326 * t114 + t324 * t115 - t155 * t461;
t216 = qJD(2) * (-pkin(2) * t462 + t297);
t446 = qJD(2) * t114 + t324 * t454 + t216;
t445 = -t164 + t505;
t442 = t247 * t460;
t437 = -t155 * t460 + t156 * t459 + qJD(1);
t436 = -pkin(2) - t565;
t435 = t462 / 0.2e1;
t434 = t461 / 0.2e1;
t433 = -t460 / 0.2e1;
t430 = t459 / 0.2e1;
t249 = -rSges(5,2) * t325 + t564;
t428 = -t249 - t330;
t427 = t269 + t571;
t425 = qJD(2) * t162 + t333 * t614;
t424 = (-t231 * t324 + t550) * qJD(2) + t613 * t333;
t423 = qJD(2) * t160 + t161 * t333 - t324 * t530;
t422 = t162 * t333 - t326 * t530 + (-t324 * t396 + t544) * qJD(2);
t129 = t162 * t523;
t421 = t326 * t158 - t129;
t418 = -t157 + t537;
t415 = t612 * t333;
t414 = t231 * t333 - t530;
t411 = t460 * t570;
t409 = (-t324 ^ 2 - t326 ^ 2) * t571;
t224 = t249 * qJD(3);
t406 = -t224 - t449;
t405 = -t330 - t569;
t63 = -t442 + (t251 - t494) * qJD(2) - t472;
t404 = t63 * t429;
t402 = t565 - t616;
t395 = -t108 * t324 - t539;
t77 = t159 * t320 + t161 * t319;
t388 = t193 * t324 + t194 * t326;
t384 = -t233 + t269;
t338 = qJD(3) ^ 2;
t382 = -qJD(3) * t224 - t330 * t338;
t381 = t405 * t338;
t214 = t282 * t324;
t201 = t247 * t324;
t373 = t394 * t324;
t206 = t234 * t333;
t364 = -pkin(4) * t458 - t206 - t449;
t363 = qJD(2) * t227 - t236 * t531 + t237 * t532;
t362 = t324 * t490 - t326 * t491;
t361 = -t336 * t484 + t337 * t486;
t348 = -t319 * t422 + t320 * t424 + t467;
t10 = t324 * t591 + t348 * t326;
t349 = qJD(2) * t157 - t319 * t423 + t320 * t425;
t11 = t349 * t324 - t326 * t590;
t12 = t348 * t324 - t326 * t591;
t56 = -t373 - t516;
t57 = -t160 * t525 - t421;
t21 = t236 * t57 - t237 * t56 + t509;
t58 = -t159 * t524 - t504;
t59 = -t160 * t524 + t503;
t80 = -t326 * t387 + t532;
t76 = t80 * qJD(2);
t22 = t236 * t59 - t237 * t58 + t76;
t354 = t613 * t236 - t614 * t237 + (-t228 + t231) * qJD(2);
t340 = -t319 * t580 + t354 * t320;
t347 = qJD(2) * t226 - t319 * t415 + t320 * t414;
t36 = t324 * t585 + t347 * t326;
t37 = t347 * t324 - t326 * t585;
t38 = t319 * t425 + t320 * t423;
t39 = t319 * t424 + t320 * t422;
t78 = t160 * t320 + t162 * t319;
t9 = t324 * t590 + t349 * t326;
t360 = (qJD(2) * t36 + t10 * t236 + t219 * t58 + t220 * t59 - t237 * t9) * t573 + (t354 * t319 + t320 * t580) * t568 + t21 * t435 + (qJD(2) * t37 - t11 * t237 + t12 * t236 + t219 * t56 + t220 * t57) * t572 + t22 * t434 + (t10 * t324 - t326 * t9 + (t324 * t58 + t326 * t59) * qJD(2)) * t576 + (t324 * t57 - t326 * t56) * t579 + (t324 * t59 - t326 * t58) * t578 + (-t11 * t326 + t12 * t324 + (t324 * t56 + t326 * t57) * qJD(2)) * t575 + (t324 * t39 - t326 * t38 + (t324 * t77 + t326 * t78) * qJD(2)) * t567 + (t324 * t363 + t326 * t340) * t577 + (t324 * t340 - t326 * t363) * t574;
t359 = (-t325 * t479 + t327 * t480) * qJD(2);
t358 = (-t336 * t470 + t337 * t471) * qJD(2);
t192 = t233 * t326;
t33 = t163 * t236 + t164 * t237 + (-t127 * t324 + t128 * t326) * qJD(3) + t437;
t49 = -t411 - t233 * t236 + (t251 - t445) * qJD(2) - t472;
t357 = t33 * (-t236 * t191 - t192 * t237) + t49 * (-qJD(2) * t192 - t234 * t236);
t356 = -t237 * t233 + t594;
t125 = -rSges(4,2) * t326 * t456 + (-t337 * t462 - t438) * rSges(4,1) + t473;
t126 = -qJD(3) * t214 + (t326 * t402 + t312) * qJD(2);
t355 = t125 * t326 + t126 * t324 + (t193 * t326 - t194 * t324) * qJD(2);
t262 = t402 * qJD(3);
t202 = t247 * t326;
t178 = t427 * t326;
t177 = t427 * t324;
t167 = t237 * t234;
t110 = -qJD(3) * t201 + (t249 * t326 + t311) * qJD(2);
t109 = -rSges(5,2) * t326 * t458 + (-t325 * t459 - t327 * t462) * rSges(5,1) + t478;
t96 = qJD(3) * t388 + qJD(1);
t75 = t244 * t324 + t273 + t281 + (-t324 * t331 + t326 * t476) * qJD(2);
t74 = t412 + t513 + (-t324 * t476 + t326 * t468) * qJD(2);
t69 = -t262 * t459 + (-t126 - t232 + t441) * qJD(2);
t68 = -t262 * t460 + t216 + (t125 - t440) * qJD(2);
t53 = (t174 * t324 + t175 * t326) * qJD(3) + t437;
t50 = t355 * qJD(3);
t48 = (t495 + t506) * qJD(2) + t356;
t41 = t382 * t326 + (-t110 + t442 + t507) * qJD(2) + t475;
t40 = t382 * t324 + (t109 + t383) * qJD(2) + t446;
t24 = -t206 * t237 + t219 * t233 + t326 * t381 + (-t75 - t95 + t411 + t507) * qJD(2) + t475;
t23 = -t206 * t236 - t220 * t233 + t324 * t381 + (t74 + t94 + t513) * qJD(2) + t446;
t14 = (t109 * t326 + t110 * t324 + (t174 * t326 + t324 * t494) * qJD(2)) * qJD(3) + t451;
t5 = t163 * t220 - t164 * t219 + t236 * t95 + t237 * t94 + (t324 * t75 + t326 * t74 + (-t127 * t326 + t324 * t505) * qJD(2)) * qJD(3) + t451;
t1 = [m(4) * t50 + m(5) * t14 + m(6) * t5; (t76 + (t57 + (t159 * t326 + t160 * t324) * t319 + t421 + t504) * t237 + (-t161 * t523 + t516 + t56 + (t159 * t324 - t160 * t326) * t319 + t503) * t236) * t574 + (t77 + t79) * t579 + (t78 + t80) * t578 + (-t509 + (t59 - t373 - t503) * t237 + (t324 * t418 - t129 + t58) * t236 + ((t158 + t394) * t236 + t418 * t237) * t326 + t21) * t577 + (t39 + t36) * t576 + ((((t636 + t640) * t326 + t619 + t632 + t638) * t326 - t592 * t324) * qJD(3) + t615) * t430 + (-t634 * qJD(3) + t222 * t327 + t223 * t325 + t259 * t337 + t260 * t336 + t319 * t414 + t320 * t415) * qJD(2) + (t24 * (t477 + t482) + t48 * t299 + t23 * (t164 + t218) + t49 * (-rSges(6,1) * t447 + t380 + t594) + (-t24 * t562 + t48 * (t333 * t563 - t244 + t450) - t23 * t331) * t324 + ((t48 * (-t234 - t261) - t49 * t331) * t326 + (t48 * (-rSges(6,3) + t331) + t49 * (-t261 - t562)) * t324) * qJD(2) - (qJD(2) * t506 + t356 - t48 + t595) * t49) * m(6) + (t41 * (-t174 - t474) + t62 * t444 + t40 * (t175 + t413) + t63 * (t298 + t478) + (t324 * t560 + t326 * t404) * qJD(3) + ((-t62 * rSges(5,3) + t63 * (-t321 - t564)) * t324 + (t62 * (-t249 - t321) - t63 * t335) * t326) * qJD(2) - (-qJD(2) * t174 + t365 + t595 - t62) * t63) * m(5) + (t69 * (t324 * t436 + t317 + t469) + t68 * t483 + t108 * (t297 + t473) + (t282 * t540 - t538) * qJD(3) + ((-pkin(2) - t402) * t539 + (t107 * (-rSges(4,3) - pkin(6)) + t108 * t436) * t324) * qJD(2) - (-qJD(2) * t193 - t107 - t235 - t440) * t108) * m(4) + (t38 + t37 + t22) * t575 + (((t326 * t593 + t592 + t617) * t326 + (t324 * t593 + t618 + t633) * t324) * qJD(3) + t603 - t607) * t433 + (t599 + t600) * t460 / 0.2e1 - (t598 - t601 + t602) * t459 / 0.2e1 + ((t597 + t629) * t324 + (t596 + t628) * t326) * qJD(3) * t567; t360 + (((t324 * t484 - t326 * t485) * t337 + (t324 * t486 + t326 * t487) * t336 + t582 * t327 + t362 * t325) * qJD(3) + (t325 * t480 + t327 * t479 + t336 * t471 + t337 * t470) * qJD(2)) * t568 + (t601 * t326 + t600 * t324 + (t324 * t597 + t326 * t596) * qJD(2)) * t567 + ((-t460 * t528 + t464) * t324 + (t359 + (-(-t325 * t489 + t327 * t491) * t326 + (-t325 * t488 + t327 * t490 + t529) * t324) * qJD(3)) * t326 + (-t460 * t526 + t463) * t324 + (t358 + (-t581 * t326 + (t527 + t361) * t324) * qJD(3)) * t326) * t433 + ((-t459 * t529 - t464) * t326 + (t359 + (-t325 * t582 + t326 * t528 + t362 * t327) * qJD(3)) * t324 + (-t459 * t527 - t463) * t326 + (t358 + (t361 * t324 + (t526 - t581) * t326) * qJD(3)) * t324) * t430 + (t48 * t167 - (t48 * (-t177 + t191) + t49 * (t178 - t453)) * qJD(2) - (t33 * t409 + (t33 * t178 + t405 * t48) * t326 + (t33 * t177 + t405 * t49) * t324) * qJD(3) - t357 + t5 * (t498 + t499) + t33 * (t448 + t452) + (t24 * t384 + t48 * t364 + t5 * t128 + t33 * t74 + (-t33 * t127 + t384 * t49) * qJD(2)) * t326 + (t23 * t384 + t49 * t364 - t5 * t127 + t33 * t75 + (t48 * (t233 + t570) + t33 * t445) * qJD(2)) * t324) * m(6) + (t50 * t388 + t96 * t355 + t395 * t262 + (-t68 * t324 - t69 * t326 + (-t108 * t326 + t540) * qJD(2)) * t282 - (t107 * t214 - t538) * qJD(2) - (t96 * (-t214 * t324 - t215 * t326) + t395 * t402) * qJD(3)) * m(4) + (t599 * qJD(2) + ((t617 * qJD(2) + t610 * t326) * t326 + (t605 * t324 + t618 * qJD(2) + (-t606 + t611) * t326) * t324) * t604) * t573 + (t598 * qJD(2) + ((t619 * qJD(2) + t606 * t326) * t326 + (t611 * t324 + t620 * qJD(2) + (-t605 + t610) * t326) * t324) * t604) * t572 + (-(t62 * t201 + t63 * (-t202 - t453)) * qJD(2) - (t53 * t409 + (-t53 * t202 + t428 * t62) * t326 + (-t53 * t201 + t428 * t63) * t324) * qJD(3) + t14 * t499 + t53 * t448 + (t41 * t429 + t62 * t406 + t14 * t175 + t53 * t109 + (t53 * t174 + t404) * qJD(2)) * t326 + (t40 * t429 + t63 * t406 + t14 * t174 + t53 * t110 + (t494 * t53 + t560) * qJD(2)) * t324) * m(5) + (t603 + t608) * t435 + (t602 + t609) * t434; 0.2e1 * (t23 * t572 + t24 * t573) * m(6) + 0.2e1 * (t40 * t572 + t41 * t573) * m(5); t360 + (t5 * t498 + t33 * (-t164 * t462 + t452) + (-t324 * t49 - t326 * t48) * t206 + (-t23 * t324 - t24 * t326 + (t324 * t48 - t326 * t49) * qJD(2)) * t233 - t48 * (qJD(2) * t191 - t167) - t357) * m(6);];
tauc = t1(:);
