% Calculate vector of inverse dynamics joint torques for
% S5RRPRR6
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:30
% EndTime: 2020-01-03 12:05:50
% DurationCPUTime: 15.27s
% Computational Cost: add. (26150->772), mult. (25223->1037), div. (0->0), fcn. (24420->10), ass. (0->364)
t378 = qJ(1) + qJ(2);
t367 = sin(t378);
t369 = cos(t378);
t383 = cos(qJ(4));
t380 = cos(pkin(9));
t381 = sin(qJ(4));
t517 = t380 * t381;
t284 = -t367 * t383 + t369 * t517;
t516 = t380 * t383;
t527 = t367 * t381;
t285 = t369 * t516 + t527;
t379 = sin(pkin(9));
t524 = t369 * t379;
t178 = Icges(5,5) * t285 - Icges(5,6) * t284 + Icges(5,3) * t524;
t264 = Icges(5,4) * t285;
t182 = Icges(5,2) * t284 - Icges(5,6) * t524 - t264;
t263 = Icges(5,4) * t284;
t184 = Icges(5,1) * t285 + Icges(5,5) * t524 - t263;
t82 = t380 * t178 - (t182 * t381 + t184 * t383) * t379;
t377 = qJ(4) + qJ(5);
t366 = sin(t377);
t368 = cos(t377);
t523 = t369 * t380;
t259 = t366 * t523 - t367 * t368;
t260 = t366 * t367 + t368 * t523;
t158 = Icges(6,5) * t260 - Icges(6,6) * t259 + Icges(6,3) * t524;
t239 = Icges(6,4) * t260;
t162 = Icges(6,2) * t259 - Icges(6,6) * t524 - t239;
t238 = Icges(6,4) * t259;
t164 = Icges(6,1) * t260 + Icges(6,5) * t524 - t238;
t75 = t380 * t158 - (t162 * t366 + t164 * t368) * t379;
t282 = -t367 * t517 - t369 * t383;
t376 = qJD(1) + qJD(2);
t201 = qJD(4) * t285 + t282 * t376;
t522 = t369 * t381;
t283 = t367 * t516 - t522;
t202 = qJD(4) * t284 + t283 * t376;
t441 = rSges(5,1) * t202 + rSges(5,2) * t201;
t520 = t376 * t379;
t478 = t367 * t520;
t115 = rSges(5,3) * t478 + t441;
t189 = t285 * rSges(5,1) - t284 * rSges(5,2) + rSges(5,3) * t524;
t487 = qJD(4) * t379;
t471 = t376 * t487;
t313 = t367 * t471;
t485 = qJDD(4) * t379;
t255 = -t369 * t485 + t313;
t276 = -rSges(5,3) * t380 + (rSges(5,1) * t383 - rSges(5,2) * t381) * t379;
t373 = qJDD(1) + qJDD(2);
t346 = -qJDD(4) * t380 + t373;
t350 = -qJD(4) * t380 + t376;
t353 = qJD(3) * t367;
t525 = t369 * t376;
t491 = qJ(3) * t525 + t353;
t530 = t367 * t376;
t228 = pkin(2) * t530 - t491;
t382 = sin(qJ(1));
t371 = t382 * pkin(1);
t384 = cos(qJ(1));
t387 = qJD(1) ^ 2;
t536 = pkin(1) * qJDD(1);
t449 = -t371 * t387 + t384 * t536;
t488 = qJD(3) * t376;
t426 = t367 * t488 + t449;
t558 = pkin(7) * t379;
t560 = pkin(3) * t380;
t446 = t558 + t560;
t288 = t446 * t369;
t309 = t369 * pkin(2) + t367 * qJ(3);
t495 = t309 + t288;
t392 = (-t446 * t530 - t228) * t376 + t495 * t373 + t426;
t301 = (-rSges(5,1) * t381 - rSges(5,2) * t383) * t379;
t292 = qJD(4) * t301;
t437 = -t292 * t487 - qJDD(3);
t49 = -t115 * t350 + t189 * t346 + t255 * t276 + t369 * t437 + t392;
t613 = t49 - g(2);
t429 = t182 * t284 + t184 * t285;
t73 = t178 * t524 + t429;
t612 = t369 * t73;
t169 = t260 * rSges(6,1) - t259 * rSges(6,2) + rSges(6,3) * t524;
t360 = pkin(4) * t383 + pkin(3);
t600 = pkin(4) * t527 + t360 * t523;
t402 = t169 + t309 + t600;
t580 = t367 * rSges(3,1) + t369 * rSges(3,2);
t274 = t580 * t376;
t555 = pkin(1) * qJD(1);
t481 = t382 * t555;
t252 = t481 + t274;
t248 = -Icges(6,3) * t380 + (Icges(6,5) * t368 - Icges(6,6) * t366) * t379;
t537 = Icges(6,4) * t368;
t249 = -Icges(6,6) * t380 + (-Icges(6,2) * t366 + t537) * t379;
t538 = Icges(6,4) * t366;
t250 = -Icges(6,5) * t380 + (Icges(6,1) * t368 - t538) * t379;
t528 = t367 * t380;
t257 = -t366 * t528 - t368 * t369;
t258 = -t369 * t366 + t368 * t528;
t529 = t367 * t379;
t104 = t248 * t529 + t249 * t257 + t250 * t258;
t375 = qJD(4) + qJD(5);
t459 = t375 * t379;
t277 = t367 * t459;
t278 = t369 * t459;
t306 = -t375 * t380 + t376;
t157 = Icges(6,5) * t258 + Icges(6,6) * t257 + Icges(6,3) * t529;
t539 = Icges(6,4) * t258;
t160 = Icges(6,2) * t257 + Icges(6,6) * t529 + t539;
t237 = Icges(6,4) * t257;
t163 = Icges(6,1) * t258 + Icges(6,5) * t529 + t237;
t66 = t157 * t529 + t257 * t160 + t258 * t163;
t609 = -t257 * t162 + t258 * t164;
t67 = -t158 * t529 - t609;
t28 = t104 * t306 + t277 * t66 - t278 * t67;
t385 = -pkin(8) - pkin(7);
t556 = pkin(7) + t385;
t187 = (t379 * t556 + t560) * t369 - t600;
t557 = -pkin(3) + t360;
t236 = t379 * t557 + t380 * t556;
t251 = -rSges(6,3) * t380 + (rSges(6,1) * t368 - rSges(6,2) * t366) * t379;
t611 = t306 * t169 - t350 * t187 + (-t236 * t487 - qJD(3)) * t369 - t278 * t251;
t610 = t350 * t189 + (-t276 * t487 - qJD(3)) * t369;
t607 = -t282 * t182 + t184 * t283;
t69 = t158 * t524 + t162 * t259 + t164 * t260;
t271 = -Icges(5,3) * t380 + (Icges(5,5) * t383 - Icges(5,6) * t381) * t379;
t540 = Icges(5,4) * t383;
t272 = -Icges(5,6) * t380 + (-Icges(5,2) * t381 + t540) * t379;
t541 = Icges(5,4) * t381;
t273 = -Icges(5,5) * t380 + (Icges(5,1) * t383 - t541) * t379;
t408 = t271 * t524 - t272 * t284 + t273 * t285;
t601 = t408 * t350;
t433 = t306 * t369;
t519 = t376 * t380;
t461 = -t375 + t519;
t436 = t366 * t461;
t171 = -t367 * t436 - t368 * t433;
t435 = t461 * t368;
t172 = -t366 * t433 + t367 * t435;
t440 = rSges(6,1) * t172 + rSges(6,2) * t171;
t102 = rSges(6,3) * t478 + t440;
t198 = t259 * rSges(6,1) + t260 * rSges(6,2);
t439 = -rSges(6,1) * t366 - rSges(6,2) * t368;
t286 = t439 * t379;
t598 = t380 * t102 + t198 * t306 + t251 * t478 + t278 * t286;
t434 = t306 * t367;
t173 = t368 * t434 - t369 * t436;
t174 = t366 * t434 + t369 * t435;
t477 = t369 * t520;
t103 = t174 * rSges(6,1) + t173 * rSges(6,2) + rSges(6,3) * t477;
t197 = t257 * rSges(6,1) - t258 * rSges(6,2);
t597 = t102 * t529 + t103 * t524 - t169 * t477 - t197 * t278 - t277 * t198;
t518 = t379 * t385;
t474 = t376 * t518;
t554 = pkin(4) * qJD(4);
t596 = -t383 * t554 - t474;
t409 = t248 * t524 - t249 * t259 + t250 * t260;
t595 = t278 * t69 + t409 * t306;
t414 = t282 * qJD(4);
t476 = t369 * t519;
t493 = pkin(3) * t476 + pkin(7) * t477;
t497 = t600 * t376;
t124 = pkin(4) * t414 - t369 * t474 - t493 + t497;
t168 = t258 * rSges(6,1) + t257 * rSges(6,2) + rSges(6,3) * t529;
t287 = pkin(3) * t528 + pkin(7) * t529;
t445 = t360 * t528 - t367 * t518;
t483 = pkin(4) * t522;
t186 = t445 - t483 - t287;
t254 = t367 * t485 + t369 * t471;
t486 = qJD(5) * t376;
t205 = (qJDD(5) * t367 + t369 * t486) * t379 + t254;
t232 = t375 * t286;
t484 = -qJDD(4) - qJDD(5);
t302 = t380 * t484 + t373;
t358 = t367 * pkin(2);
t307 = -qJ(3) * t369 + t358;
t489 = qJD(3) * t369;
t492 = pkin(2) * t525 + qJ(3) * t530;
t443 = -t489 + t492;
t372 = t384 * pkin(1);
t490 = t387 * t372 + t382 * t536;
t456 = t373 * t307 + t376 * t443 + t490;
t405 = t373 * t287 - t369 * t488 + t376 * t493 + t456;
t559 = pkin(4) * t381;
t425 = qJD(4) ^ 2 * t379 ^ 2 * t559 - qJDD(3);
t20 = t103 * t306 + t124 * t350 + t168 * t302 + t186 * t346 - t205 * t251 - t232 * t277 - t236 * t254 + t367 * t425 + t405;
t594 = t20 - g(3);
t203 = -qJD(4) * t283 - t284 * t376;
t204 = t285 * t376 + t414;
t116 = t204 * rSges(5,1) + t203 * rSges(5,2) + rSges(5,3) * t477;
t188 = t283 * rSges(5,1) + t282 * rSges(5,2) + rSges(5,3) * t529;
t48 = t116 * t350 + t188 * t346 - t254 * t276 + t367 * t437 + t405;
t593 = t48 - g(3);
t342 = rSges(4,1) * t528;
t448 = -rSges(4,2) * t529 + t342;
t234 = -rSges(4,3) * t369 + t448;
t410 = rSges(4,1) * t476 + rSges(4,3) * t530 + (-rSges(4,2) * t520 - qJD(3)) * t369;
t85 = -qJDD(3) * t367 + t373 * t234 + t376 * t410 + t456;
t592 = t85 - g(3);
t235 = rSges(4,1) * t523 - rSges(4,2) * t524 + rSges(4,3) * t367;
t220 = t309 + t235;
t494 = rSges(4,2) * t478 + rSges(4,3) * t525;
t86 = -qJDD(3) * t369 + (-t342 * t376 - t228 + t494) * t376 + t220 * t373 + t426;
t591 = t86 - g(2);
t275 = rSges(3,1) * t525 - rSges(3,2) * t530;
t590 = t275 * t376 + t373 * t580 - g(3) + t490;
t310 = rSges(3,1) * t369 - t367 * rSges(3,2);
t589 = -t274 * t376 + t310 * t373 - g(2) + t449;
t457 = t554 * t517;
t473 = t596 * t367 - t376 * t483;
t123 = t369 * t457 + (t380 * t557 - t558) * t530 + t473;
t206 = t313 + (t367 * t486 + t369 * t484) * t379;
t21 = -t306 * t102 - t350 * t123 + t169 * t302 - t346 * t187 + t206 * t251 - t278 * t232 + t255 * t236 + t369 * t425 + t392;
t447 = -t353 + t481;
t498 = t287 + t307;
t585 = t376 * t498;
t411 = t447 + t585;
t472 = t367 * t487;
t573 = t168 * t306 + t186 * t350 - t236 * t472 - t277 * t251;
t61 = t411 + t573;
t365 = t384 * t555;
t431 = t376 * t495 + t365;
t62 = t431 + t611;
t588 = (-t61 * t457 + t62 * (-rSges(6,3) * t379 - t360 * t380 - pkin(2)) * t376) * t367 + (-t21 * t518 - t62 * t457 + t61 * (-qJD(3) + t596)) * t369;
t88 = t431 + t610;
t587 = t88 * (-t560 - pkin(2) + (-rSges(5,3) - pkin(7)) * t379) * t530;
t584 = t376 * (t234 + t307);
t583 = -t259 * t160 + t260 * t163;
t542 = Icges(5,4) * t283;
t180 = Icges(5,2) * t282 + Icges(5,6) * t529 + t542;
t262 = Icges(5,4) * t282;
t183 = Icges(5,1) * t283 + Icges(5,5) * t529 + t262;
t515 = -t284 * t180 + t285 * t183;
t297 = t376 * t309;
t579 = -t376 * t235 - t297 + t410 + t489 + t492;
t501 = t376 * t288 + t297;
t578 = t116 + t443 + t493 - t501 - t610;
t576 = t103 + t492 + t497 - t501 - t611;
t575 = t188 * t350 - t276 * t472;
t574 = t495 + t189;
t412 = t367 * (-Icges(5,2) * t283 + t183 + t262) - t369 * (Icges(5,2) * t285 - t184 + t263);
t572 = t367 * (-Icges(5,1) * t282 + t180 + t542) - t369 * (-Icges(5,1) * t284 + t182 - t264);
t280 = (-Icges(6,2) * t368 - t538) * t379;
t393 = t277 * (-Icges(6,2) * t258 + t163 + t237) - t278 * (Icges(6,2) * t260 - t164 + t238) + t306 * (t250 + t280);
t281 = (-Icges(6,1) * t366 - t537) * t379;
t571 = t277 * (-Icges(6,1) * t257 + t160 + t539) - t278 * (-Icges(6,1) * t259 + t162 - t239) + t306 * (t249 - t281);
t570 = t205 / 0.2e1;
t569 = t206 / 0.2e1;
t568 = t254 / 0.2e1;
t567 = t255 / 0.2e1;
t566 = -t277 / 0.2e1;
t565 = t277 / 0.2e1;
t564 = -t278 / 0.2e1;
t563 = t278 / 0.2e1;
t561 = -t380 / 0.2e1;
t95 = Icges(6,5) * t174 + Icges(6,6) * t173 + Icges(6,3) * t477;
t97 = Icges(6,4) * t174 + Icges(6,2) * t173 + Icges(6,6) * t477;
t99 = Icges(6,1) * t174 + Icges(6,4) * t173 + Icges(6,5) * t477;
t36 = -t380 * t95 + ((-t160 * t375 + t99) * t368 + (-t163 * t375 - t97) * t366) * t379;
t551 = t36 * t277;
t94 = Icges(6,5) * t172 + Icges(6,6) * t171 + Icges(6,3) * t478;
t96 = Icges(6,4) * t172 + Icges(6,2) * t171 + Icges(6,6) * t478;
t98 = Icges(6,1) * t172 + Icges(6,4) * t171 + Icges(6,5) * t478;
t37 = -t380 * t94 + ((-t162 * t375 + t98) * t368 + (t164 * t375 - t96) * t366) * t379;
t550 = t37 * t278;
t74 = -t157 * t380 + (-t160 * t366 + t163 * t368) * t379;
t549 = t74 * t205;
t548 = t75 * t206;
t177 = Icges(5,5) * t283 + Icges(5,6) * t282 + Icges(5,3) * t529;
t81 = -t177 * t380 + (-t180 * t381 + t183 * t383) * t379;
t547 = t81 * t254;
t546 = t82 * t255;
t122 = -t248 * t380 + (-t249 * t366 + t250 * t368) * t379;
t279 = (-Icges(6,5) * t366 - Icges(6,6) * t368) * t379;
t229 = t375 * t279;
t230 = t375 * t280;
t231 = t375 * t281;
t83 = -t229 * t380 + ((-t249 * t375 + t231) * t368 + (-t250 * t375 - t230) * t366) * t379;
t545 = t122 * t302 + t83 * t306;
t298 = (-Icges(5,5) * t381 - Icges(5,6) * t383) * t379;
t289 = qJD(4) * t298;
t299 = (-Icges(5,2) * t383 - t541) * t379;
t290 = qJD(4) * t299;
t300 = (-Icges(5,1) * t381 - t540) * t379;
t291 = qJD(4) * t300;
t106 = -t289 * t380 + (-t290 * t381 + t291 * t383 + (-t272 * t383 - t273 * t381) * qJD(4)) * t379;
t134 = -t271 * t380 + (-t272 * t381 + t273 * t383) * t379;
t544 = t106 * t350 + t134 * t346;
t534 = t177 * t369;
t533 = t178 * t367;
t510 = -t168 - t186;
t504 = -t236 - t251;
t500 = t272 - t300;
t499 = t273 + t299;
t70 = t177 * t529 + t282 * t180 + t283 * t183;
t71 = -t178 * t529 - t607;
t470 = t529 / 0.2e1;
t469 = -t524 / 0.2e1;
t468 = t520 / 0.2e1;
t467 = -t487 / 0.2e1;
t466 = t487 / 0.2e1;
t14 = t102 * t277 + t103 * t278 - t168 * t206 - t169 * t205 - t186 * t255 + t187 * t254 + (t123 * t367 + t124 * t369) * t487;
t465 = t14 * (t168 * t524 - t169 * t529);
t463 = t306 * t197 - t277 * t286;
t253 = t310 * t376 + t365;
t455 = t367 * t468;
t454 = t369 * t468;
t453 = t367 * t467;
t452 = t367 * t466;
t451 = t369 * t467;
t450 = t369 * t466;
t266 = t284 * pkin(4);
t339 = rSges(2,1) * t384 - t382 * rSges(2,2);
t338 = rSges(2,1) * t382 + rSges(2,2) * t384;
t87 = t411 + t575;
t438 = -t367 * t87 - t369 * t88;
t432 = t353 - t585;
t428 = t188 * t369 - t189 * t367;
t427 = (Icges(5,5) * t282 - Icges(5,6) * t283) * t367 - (Icges(5,5) * t284 + Icges(5,6) * t285) * t369;
t68 = -t157 * t524 - t583;
t419 = (t367 * t70 - t369 * t71) * t379;
t72 = -t177 * t524 - t515;
t418 = (t367 * t72 - t612) * t379;
t417 = -t441 + t491;
t265 = t282 * pkin(4);
t415 = (Icges(6,5) * t257 - Icges(6,6) * t258) * t277 - (Icges(6,5) * t259 + Icges(6,6) * t260) * t278 + t279 * t306;
t107 = t428 * t487;
t133 = t188 + t498;
t219 = t358 + (-rSges(4,3) - qJ(3)) * t369 + t448;
t403 = -t440 - t473 + t491;
t399 = (-rSges(4,1) * t380 - pkin(2)) * t530 + t491 + t494;
t398 = -t369 * t518 + t402;
t16 = t160 * t171 + t163 * t172 + t259 * t97 - t260 * t99 + (t157 * t530 - t369 * t95) * t379;
t17 = t162 * t171 - t164 * t172 + t259 * t96 - t260 * t98 + (-t158 * t530 - t369 * t94) * t379;
t18 = t160 * t173 + t163 * t174 + t257 * t97 + t258 * t99 + (t157 * t525 + t367 * t95) * t379;
t19 = t162 * t173 - t164 * t174 + t257 * t96 + t258 * t98 + (-t158 * t525 + t367 * t94) * t379;
t29 = t277 * t68 - t595;
t53 = t171 * t249 + t172 * t250 + t230 * t259 - t231 * t260 + (-t229 * t369 + t248 * t530) * t379;
t54 = t173 * t249 + t174 * t250 + t230 * t257 + t231 * t258 + (t229 * t367 + t248 * t525) * t379;
t396 = (t104 * t302 + t18 * t277 - t19 * t278 + t205 * t66 + t206 * t67 + t306 * t54) * t470 + (t257 * t393 - t258 * t571 + t415 * t529) * t566 + (t393 * t259 + t260 * t571 - t415 * t524) * t563 - (-t415 * t380 + (-t366 * t393 - t368 * t571) * t379) * t306 / 0.2e1 + (t16 * t277 - t17 * t278 + t205 * t68 + t206 * t69 - t302 * t409 + t306 * t53) * t469 + t29 * t455 + t28 * t454 + (-t104 * t380 + (t367 * t66 - t369 * t67) * t379) * t570 + (t409 * t380 + (t367 * t68 - t369 * t69) * t379) * t569 + (-t380 * t54 + ((t376 * t66 - t19) * t369 + (t376 * t67 + t18) * t367) * t379) * t565 + (-t380 * t53 + ((t376 * t68 - t17) * t369 + (t376 * t69 + t16) * t367) * t379) * t564 + t302 * (-t122 * t380 + (t367 * t74 - t369 * t75) * t379) / 0.2e1 + (t545 + t548 + t549 - t550 + t551) * t561 + t306 * (-t380 * t83 + ((t376 * t74 - t37) * t369 + (t376 * t75 + t36) * t367) * t379) / 0.2e1;
t121 = t358 + (-qJ(3) - t559) * t369 + t445 + t168;
t117 = t271 * t529 + t272 * t282 + t273 * t283;
t108 = t117 * t350;
t41 = qJD(4) * t419 + t108;
t42 = qJD(4) * t418 - t601;
t110 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t477;
t112 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t477;
t114 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t477;
t46 = -t110 * t380 + (-t112 * t381 + t114 * t383 + (-t180 * t383 - t183 * t381) * qJD(4)) * t379;
t109 = Icges(5,5) * t202 + Icges(5,6) * t201 + Icges(5,3) * t478;
t111 = Icges(5,4) * t202 + Icges(5,2) * t201 + Icges(5,6) * t478;
t113 = Icges(5,1) * t202 + Icges(5,4) * t201 + Icges(5,5) * t478;
t47 = -t109 * t380 + (-t111 * t381 + t113 * t383 + (-t182 * t383 + t184 * t381) * qJD(4)) * t379;
t59 = t201 * t272 + t202 * t273 + t284 * t290 - t285 * t291 + (t271 * t530 - t289 * t369) * t379;
t60 = t203 * t272 + t204 * t273 + t282 * t290 + t283 * t291 + (t271 * t525 + t289 * t367) * t379;
t389 = (t108 + ((-t429 + t70 + t73) * t367 + (t72 + (-t533 + t534) * t379 - t71 + t515) * t369) * t487) * t450 - t409 * t569 - t408 * t567 - t550 / 0.2e1 + t104 * t570 + t54 * t565 + t117 * t568 + t548 / 0.2e1 + t549 / 0.2e1 + t546 / 0.2e1 + t547 / 0.2e1 + t544 + t545 + t551 / 0.2e1 + t28 * t563 + (t29 + (t67 + (t157 * t369 + t158 * t367) * t379 + t583 + t609) * t277 + t595) * t566 + (t53 + t28) * t564 + (t46 + t60) * t452 + (Icges(3,3) + Icges(4,2) * t380 ^ 2 + (Icges(4,1) * t379 + 0.2e1 * Icges(4,4) * t380) * t379) * t373 + (t42 + t601 + (t612 + (t515 + t71 + (t533 + t534) * t379 + t607) * t367) * t487) * t453 + (t47 + t59 + t41) * t451;
t215 = rSges(5,1) * t284 + rSges(5,2) * t285;
t214 = rSges(5,1) * t282 - rSges(5,2) * t283;
t167 = t220 * t376 + t365 - t489;
t166 = t447 + t584;
t150 = t380 * t169;
t57 = t168 * t278 - t169 * t277 + (t186 * t369 + t187 * t367) * t487;
t33 = t111 * t282 + t113 * t283 + t182 * t203 - t184 * t204 + (t109 * t367 - t178 * t525) * t379;
t32 = t112 * t282 + t114 * t283 + t180 * t203 + t183 * t204 + (t110 * t367 + t177 * t525) * t379;
t31 = t111 * t284 - t113 * t285 + t182 * t201 - t184 * t202 + (-t109 * t369 - t178 * t530) * t379;
t30 = t112 * t284 - t114 * t285 + t180 * t201 + t183 * t202 + (-t110 * t369 + t177 * t530) * t379;
t1 = [Icges(2,3) * qJDD(1) + t389 + (t589 * (t310 + t372) + t590 * (t371 + t580) + (-t253 + t275 + t365) * t252) * m(3) + ((t338 ^ 2 + t339 ^ 2) * qJDD(1) - g(2) * t339 - g(3) * t338) * m(2) + (t21 * (t372 + t402) + t62 * (t403 - t481) - g(2) * (t372 + t398) + t594 * (t371 + t121) + (t62 + t576) * t61 + t588) * m(6) + (t88 * (t417 - t481) + t593 * (t133 + t371) + (t88 + t578) * t87 + t587 + t613 * (t372 + t574)) * m(5) + (t167 * (t399 - t481) + t591 * (t220 + t372) + t592 * (t371 + t219) + (t167 + t579) * t166) * m(4); t389 + (-g(2) * t398 + t21 * t402 + (t403 - t432 + t573) * t62 + t576 * t61 + t594 * t121 + t588) * m(6) + ((t417 - t432 + t575) * t88 + t578 * t87 + t593 * t133 + t587 + t613 * t574) * m(5) + (t591 * t220 + t592 * t219 + (-t353 + t399 + t584) * t167 + t579 * t166) * m(4) + (t252 * t275 - t253 * t274 + (-t252 * t376 + t589) * t310 + (t253 * t376 + t590) * t580) * m(3); (-m(4) - m(5) - m(6)) * (-g(2) * t369 - g(3) * t367) + m(4) * (-t367 * t85 - t369 * t86) + m(5) * (-t367 * t48 - t369 * t49) + m(6) * (-t20 * t367 - t21 * t369); t42 * t455 + t41 * t454 + t346 * (-t134 * t380 + (t367 * t81 - t369 * t82) * t379) / 0.2e1 + (t380 * t408 + t418) * t567 + (-t117 * t380 + t419) * t568 + t350 * (-t106 * t380 + ((t376 * t81 - t47) * t369 + (t376 * t82 + t46) * t367) * t379) / 0.2e1 + (-t408 * t346 + t254 * t72 + t255 * t73 + t350 * t59 + (t30 * t367 - t31 * t369) * t487) * t469 + t396 + (t117 * t346 + t254 * t70 + t255 * t71 + t350 * t60 + (t32 * t367 - t33 * t369) * t487) * t470 + (-t380 * t59 + ((t376 * t72 - t31) * t369 + (t376 * t73 + t30) * t367) * t379) * t451 + (-t380 * t60 + ((t376 * t70 - t33) * t369 + (t376 * t71 + t32) * t367) * t379) * t452 + (t547 + t546 + (t367 * t46 - t369 * t47) * t487 + t544) * t561 + ((t284 * t499 + t285 * t500 - t298 * t524) * t350 + (t412 * t284 + t285 * t572 - t427 * t524) * t487) * t450 + ((t282 * t499 - t283 * t500 + t298 * t529) * t350 + (t282 * t412 - t283 * t572 + t427 * t529) * t487) * t453 - t350 * (-t380 * t298 * t350 + ((-t381 * t499 - t383 * t500) * t350 + ((-t381 * t412 - t383 * t572) * t379 - t427 * t380) * qJD(4)) * t379) / 0.2e1 + (-t61 * (t265 * t350 + t463) + t465 - t21 * t150 + (t21 * t187 + t20 * t510 + t61 * (-t103 - t124)) * t380 - g(2) * (t265 + t197) - g(3) * (t266 + t198) + (t123 * t380 + t266 * t350 + t598) * t62 + (-(t265 * t369 + t266 * t367) * t487 + t597) * t57 + ((t14 * t186 + t57 * t124 + t21 * t504 - t62 * t232 + (t57 * t187 + t504 * t61) * t376) * t369 + (t14 * t187 + t57 * t123 + t20 * t504 - t61 * t232 + (t62 * t236 + t510 * t57) * t376) * t367 - g(1) * (t439 - t559)) * t379) * m(6) + ((t115 * t88 - t116 * t87 - t188 * t48 - t189 * t49) * t380 + ((-t188 * t255 - t189 * t254) * t428 + t438 * t292 + ((-t376 * t87 - t49) * t369 + (t376 * t88 - t48) * t367) * t276 + (0.2e1 * t115 * t367 + 0.2e1 * t116 * t369 - t188 * t530 - t189 * t525) * t107) * t379 - (t214 * t87 - t215 * t88) * t350 - (t107 * (t214 * t369 + t215 * t367) + t438 * t301) * t487 - g(1) * t301 - g(2) * t214 - g(3) * t215) * m(5); (t465 + t21 * (-t251 * t524 - t150) + t20 * (-t168 * t380 - t251 * t529) - g(1) * t286 - g(2) * t197 - g(3) * t198 + (-t232 * t524 + t598) * t62 + (-t103 * t380 + (-t232 * t367 - t251 * t525) * t379 - t463) * t61 + (-t168 * t478 + t597) * t57) * m(6) + t396;];
tau = t1;
