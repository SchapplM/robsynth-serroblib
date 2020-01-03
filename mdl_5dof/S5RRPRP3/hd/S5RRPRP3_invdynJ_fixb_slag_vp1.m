% Calculate vector of inverse dynamics joint torques for
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:55
% EndTime: 2019-12-31 19:51:12
% DurationCPUTime: 13.20s
% Computational Cost: add. (14439->522), mult. (11347->613), div. (0->0), fcn. (8712->8), ass. (0->303)
t603 = Icges(6,4) + Icges(5,5);
t602 = Icges(5,6) - Icges(6,6);
t311 = qJ(1) + qJ(2);
t305 = sin(t311);
t309 = pkin(8) + qJ(4);
t303 = sin(t309);
t304 = cos(t309);
t215 = Icges(5,5) * t304 - Icges(5,6) * t303;
t306 = cos(t311);
t352 = t215 * t306;
t150 = Icges(5,3) * t305 + t352;
t217 = Icges(6,4) * t304 + Icges(6,6) * t303;
t353 = t217 * t306;
t152 = Icges(6,2) * t305 + t353;
t601 = t150 + t152;
t277 = Icges(6,5) * t303;
t374 = Icges(6,1) * t304 + t277;
t155 = -Icges(6,4) * t306 + t305 * t374;
t461 = t303 * t305;
t253 = Icges(5,4) * t461;
t459 = t304 * t305;
t157 = Icges(5,1) * t459 - Icges(5,5) * t306 - t253;
t593 = -t155 - t157;
t355 = t374 * t306;
t156 = Icges(6,4) * t305 + t355;
t474 = Icges(5,4) * t303;
t223 = Icges(5,1) * t304 - t474;
t356 = t223 * t306;
t158 = Icges(5,5) * t305 + t356;
t592 = t156 + t158;
t372 = Icges(6,3) * t304 - t277;
t576 = -Icges(5,2) * t304 - t372 - t474;
t473 = Icges(6,5) * t304;
t220 = Icges(6,1) * t303 - t473;
t278 = Icges(5,4) * t304;
t600 = Icges(5,1) * t303 + t220 + t278;
t213 = Icges(6,3) * t303 + t473;
t373 = -Icges(5,2) * t303 + t278;
t591 = t213 - t373;
t590 = t603 * t303 + t602 * t304;
t598 = t223 + t374;
t458 = t304 * t306;
t252 = Icges(6,5) * t458;
t460 = t303 * t306;
t148 = Icges(6,6) * t305 + Icges(6,3) * t460 + t252;
t597 = t148 * t460 + t601 * t305 + t592 * t458;
t151 = -Icges(6,2) * t306 + t217 * t305;
t139 = t305 * t151;
t147 = -Icges(6,6) * t306 + t213 * t305;
t149 = Icges(5,5) * t459 - Icges(5,6) * t461 - Icges(5,3) * t306;
t596 = -t147 * t460 - t305 * t149 + t593 * t458 - t139;
t153 = Icges(5,4) * t459 - Icges(5,2) * t461 - Icges(5,6) * t306;
t595 = t147 - t153;
t354 = t373 * t306;
t154 = Icges(5,6) * t305 + t354;
t594 = t148 - t154;
t310 = qJD(1) + qJD(2);
t588 = t576 * qJD(4) + t602 * t310;
t587 = -t600 * qJD(4) + t603 * t310;
t584 = t576 * t303 + t600 * t304;
t470 = t153 * t303;
t367 = -t157 * t304 + t470;
t454 = t306 * t151;
t370 = t147 * t303 + t155 * t304;
t521 = t305 * t370;
t57 = -t454 + t521;
t528 = -t306 * t149 - t305 * t367 + t57;
t526 = -t153 * t460 - t596;
t525 = -t154 * t460 + t597;
t388 = -t148 * t461 + t152 * t306 - t156 * t459;
t126 = t158 * t459;
t396 = t306 * t150 - t126;
t60 = -t154 * t461 - t396;
t527 = t60 - t388;
t586 = t591 * qJD(4);
t585 = t598 * qJD(4);
t577 = t215 + t217;
t583 = -t303 * t600 + t304 * t576;
t519 = t590 * t306;
t520 = t590 * t305;
t554 = rSges(6,1) + pkin(4);
t547 = t305 * t584 - t519;
t546 = t306 * t584 + t520;
t545 = rSges(6,3) + qJ(5);
t457 = t305 * t310;
t582 = t588 * t306 + t591 * t457;
t453 = t306 * t310;
t581 = -t213 * t453 + t305 * t588 + t310 * t354;
t580 = -t592 * t303 + t594 * t304;
t524 = t593 * t303 + t595 * t304;
t579 = t306 * t587 - t457 * t598;
t578 = (t355 + t356) * t310 + t587 * t305;
t575 = (-Icges(6,2) - Icges(5,3)) * t310 + t590 * qJD(4);
t469 = t154 * t303;
t574 = t148 * t303 + t592 * t304 - t469;
t573 = -t367 + t370;
t514 = t545 * t459;
t572 = t583 * qJD(4) + t586 * t303 + t585 * t304 + t310 * t590;
t571 = t525 * t305 - t526 * t306;
t570 = t527 * t305 - t528 * t306;
t385 = t304 * rSges(6,1) + t303 * rSges(6,3);
t569 = t304 * pkin(4) + t303 * qJ(5) + t385;
t568 = -t577 * qJD(4) + t584 * t310;
t567 = (t576 * t306 + t592) * t305 - (-Icges(5,2) * t459 - t372 * t305 - t253 - t593) * t306;
t566 = t546 * t310;
t434 = t554 * t303 - t545 * t304;
t565 = (-t149 - t151) * t310 - t578 * t304 + t581 * t303 - t524 * qJD(4);
t564 = t580 * qJD(4) - t582 * t303 + t579 * t304 + t601 * t310;
t563 = -t595 * t306 + (-Icges(6,1) * t460 + t220 * t306 + t252 + t594) * t305;
t562 = -t591 + t600;
t561 = t576 + t598;
t558 = t547 * t310;
t557 = (-t353 - t352 + t573) * t310 + t575 * t305;
t556 = t575 * t306 + t310 * t574 + t457 * t577;
t424 = qJD(4) * t305;
t543 = t305 * rSges(6,2) + pkin(4) * t458;
t442 = rSges(6,1) * t458 + t460 * t545 + t543;
t283 = t305 * qJ(3);
t232 = t306 * pkin(2) + t283;
t313 = cos(pkin(8));
t299 = pkin(3) * t313 + pkin(2);
t255 = t306 * t299;
t314 = -pkin(7) - qJ(3);
t456 = t305 * t314;
t394 = t255 - t456;
t145 = t394 - t232;
t448 = t145 + t232;
t555 = -t434 * t424 + t310 * (t442 + t448);
t553 = qJD(4) * t570 + t558;
t552 = qJD(4) * t571 + t566;
t551 = t573 * qJD(4) + t578 * t303 + t581 * t304;
t550 = t574 * qJD(4) + t579 * t303 + t582 * t304;
t549 = -t305 * t568 + t572 * t306;
t548 = t305 * t572 + t568 * t306;
t297 = t306 * rSges(6,2);
t443 = t305 * t569 - t297;
t542 = t454 + t597;
t406 = t304 * t424;
t541 = t303 * t453 + t406;
t315 = sin(qJ(1));
t485 = pkin(1) * qJD(1);
t415 = t315 * t485;
t231 = rSges(3,1) * t305 + rSges(3,2) * t306;
t462 = t231 * t310;
t169 = -t415 - t462;
t422 = qJD(4) * t310;
t198 = -qJDD(4) * t306 + t305 * t422;
t308 = qJDD(1) + qJDD(2);
t420 = qJD(5) * t304;
t441 = -qJD(4) * t569 + t420;
t341 = qJDD(5) * t303 + (t420 + t441) * qJD(4);
t316 = cos(qJ(1));
t317 = qJD(1) ^ 2;
t351 = (-qJDD(1) * t315 - t316 * t317) * pkin(1);
t425 = qJD(3) * t310;
t344 = qJDD(3) * t305 + t306 * t425 + t351;
t285 = t306 * qJ(3);
t230 = pkin(2) * t305 - t285;
t274 = t306 * t314;
t428 = -t305 * t299 - t274;
t144 = t230 + t428;
t449 = t144 - t230;
t393 = -t443 + t449;
t421 = qJD(5) * t303;
t404 = t305 * t421;
t282 = qJD(3) * t306;
t165 = t232 * t310 - t282;
t245 = t310 * t456;
t488 = pkin(2) - t299;
t452 = t245 - (-t306 * t488 - t283) * t310 - t165;
t408 = t303 * t424;
t515 = t554 * t408;
t477 = rSges(6,3) * t406 + t404 + t541 * qJ(5) - t515 + (t306 * t385 + t543) * t310;
t5 = t434 * t198 + t393 * t308 + t341 * t306 + (-t404 + t452 - t477) * t310 + t344;
t540 = t5 - g(1);
t197 = qJDD(4) * t305 + t306 * t422;
t246 = t306 * t421;
t265 = qJ(3) * t453;
t307 = t316 * pkin(1);
t490 = pkin(1) * t315;
t391 = qJDD(1) * t307 - t317 * t490;
t281 = qJD(3) * t305;
t427 = t265 + t281;
t334 = -qJDD(3) * t306 + t310 * (-pkin(2) * t457 + t427) + t308 * t232 + t305 * t425 + t391;
t330 = t310 * (-t265 + (t305 * t488 - t274) * t310) + t308 * t145 + t334;
t423 = qJD(4) * t306;
t407 = t303 * t423;
t349 = -t304 * t457 - t407;
t414 = t303 * t457;
t405 = t304 * t423;
t516 = rSges(6,2) * t453 + t545 * t405;
t478 = t349 * t554 - t545 * t414 + t246 + t516;
t6 = t442 * t308 - t434 * t197 + (t246 + t478) * t310 + t341 * t305 + t330;
t539 = t6 - g(2);
t538 = t557 * t305 + t565 * t306;
t537 = -t556 * t305 + t564 * t306;
t359 = rSges(5,1) * t458 + t305 * rSges(5,3);
t410 = -rSges(5,1) * t408 - rSges(5,2) * t541;
t109 = t310 * t359 + t410;
t229 = rSges(5,1) * t304 - rSges(5,2) * t303;
t209 = t229 * qJD(4);
t226 = rSges(5,1) * t303 + rSges(5,2) * t304;
t418 = rSges(5,1) * t459;
t161 = -rSges(5,2) * t461 - t306 * rSges(5,3) + t418;
t412 = -t161 + t449;
t20 = -t209 * t423 + t198 * t226 + (-t109 + t452) * t310 + t412 * t308 + t344;
t536 = t20 - g(1);
t357 = rSges(5,3) * t453 + (-t405 + t414) * rSges(5,2);
t107 = rSges(5,1) * t349 + t357;
t163 = -rSges(5,2) * t460 + t359;
t21 = t107 * t310 + t163 * t308 - t197 * t226 - t209 * t424 + t330;
t535 = t21 - g(2);
t312 = sin(pkin(8));
t486 = rSges(4,2) * t312;
t417 = t306 * t486;
t237 = t310 * t417;
t487 = rSges(4,1) * t313;
t419 = t305 * t487;
t272 = t305 * t486;
t426 = t306 * rSges(4,3) + t272;
t166 = t419 - t426;
t432 = -t230 - t166;
t512 = -t305 * rSges(4,3) - t306 * t487;
t49 = t432 * t308 + (t310 * t512 - t165 + t237) * t310 + t344;
t534 = t49 - g(1);
t167 = -t417 - t512;
t431 = rSges(4,3) * t453 + t310 * t272;
t50 = t308 * t167 + t310 * (-t310 * t419 + t431) + t334;
t533 = t50 - g(2);
t195 = rSges(3,1) * t453 - rSges(3,2) * t457;
t532 = -t195 * t310 - t231 * t308 - g(1) + t351;
t233 = t306 * rSges(3,1) - rSges(3,2) * t305;
t531 = t233 * t308 - t310 * t462 - g(2) + t391;
t530 = t565 * t305 - t557 * t306;
t529 = t564 * t305 + t556 * t306;
t134 = t167 + t232;
t522 = t134 * t310;
t201 = t310 * t230;
t517 = -t310 * t144 + t201;
t513 = t545 * t458;
t416 = t316 * t485;
t389 = -t282 + t416;
t511 = -t303 * t567 + t563 * t304;
t510 = (-t562 * t303 + t561 * t304) * t310;
t429 = t246 + t281;
t350 = -t423 * t434 + t429;
t509 = t310 * t443 - t350 + t429 + t516 + t517;
t508 = t577 * t310;
t507 = t310 * t166 + t201 + t431;
t506 = t245 + t515;
t143 = t310 * t161;
t499 = -rSges(5,1) * t407 + t143 + t281 + t357 + t517;
t498 = t197 / 0.2e1;
t497 = t198 / 0.2e1;
t496 = t305 / 0.2e1;
t495 = -t306 / 0.2e1;
t489 = g(2) * t305;
t387 = -t226 * t423 + t281;
t348 = t387 - t415;
t55 = t310 * t412 + t348;
t483 = t310 * t55;
t52 = -t420 + (t443 * t305 + t442 * t306) * qJD(4);
t482 = t52 * t303;
t440 = -t461 * t554 + t514;
t439 = -t460 * t554 + t513;
t411 = t163 + t448;
t409 = t554 * t306;
t401 = -pkin(2) - t487;
t400 = -t424 / 0.2e1;
t399 = t424 / 0.2e1;
t398 = -t423 / 0.2e1;
t397 = t423 / 0.2e1;
t395 = -t149 + t469;
t271 = rSges(2,1) * t316 - rSges(2,2) * t315;
t270 = rSges(2,1) * t315 + rSges(2,2) * t316;
t45 = t310 * t393 + t350 - t415;
t347 = t404 + t389;
t46 = t347 + t555;
t383 = t305 * t46 + t306 * t45;
t193 = t226 * t424;
t56 = t310 * t411 - t193 + t389;
t382 = -t305 * t56 - t306 * t55;
t377 = t245 + t282 - t410;
t365 = t161 * t305 + t163 * t306;
t358 = t383 * t304;
t117 = -t161 + t428;
t133 = t305 * t401 + t285 + t426;
t85 = t394 + t442;
t118 = t163 + t394;
t333 = -t299 - t569;
t84 = t297 + (-t303 * t545 - t304 * t554) * t305 + t428;
t110 = t310 * t432 + t281 - t415;
t111 = t389 + t522;
t323 = (t110 * t401 * t306 + (t110 * (-rSges(4,3) - qJ(3)) + t111 * t401) * t305) * t310;
t320 = (t55 * (-t359 - t255) + t56 * (-t418 + t428)) * t310;
t319 = (((t60 - t126 + (t150 + t470) * t306 + t596) * t306 + (t57 - t521 + t542) * t305) * qJD(4) + t566) * t397 + (t584 * qJD(4) + t585 * t303 - t586 * t304) * t310 + (-t580 + t546) * t498 + (-t524 + t547) * t497 + (t549 + t550) * t399 + (Icges(3,3) + Icges(4,2) * t313 ^ 2 + (Icges(4,1) * t312 + 0.2e1 * Icges(4,4) * t313) * t312 - t583) * t308 + (((t306 * t395 + t525 - t542) * t306 + (t305 * t395 - t139 + t388 + t396 + t526) * t305) * qJD(4) + t553 - t558) * t400 + (t548 + t551 + t552) * t398;
t318 = ((-t46 * t314 + t333 * t45) * t306 + (-t45 * rSges(6,2) + t333 * t46) * t305) * t310 + (-t303 * t409 * t46 - t45 * t514) * qJD(4);
t189 = t226 * t306;
t185 = t226 * t305;
t170 = t233 * t310 + t416;
t79 = t365 * qJD(4);
t7 = -qJDD(5) * t304 - t442 * t198 + t443 * t197 + (t305 * t477 + t306 * t478 + t421) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t319 + (t531 * (t233 + t307) + t532 * (-t231 - t490) + (-t195 - t416 + t170) * t169) * m(3) + ((t270 ^ 2 + t271 ^ 2) * qJDD(1) + g(1) * t270 - g(2) * t271) * m(2) + (t45 * (-t347 + t506) + t318 + t539 * (t307 + t85) + t540 * (t84 - t490) + (t45 + t509) * t46) * m(6) + (t55 * (t377 - t416) + t320 + (-t348 + t55 - t415 + t499) * t56 + t535 * (t118 + t307) + t536 * (t117 - t490)) * m(5) + (t110 * (t237 - t389) + t323 + t533 * (t134 + t307) + t534 * (t133 - t490) + (t110 + t265 + t507) * t111) * m(4); t319 + (t509 * t46 + t539 * t85 + t540 * t84 + t318 + (t506 + t555) * t45) * m(6) + (t411 * t483 + t320 + (-t387 + t499) * t56 + (t377 - t193 - t282) * t55 + t535 * t118 + t536 * t117) * m(5) + (t323 + t533 * t134 + t534 * t133 + (t427 - t281 + t507) * t111 + (t237 + t522) * t110) * m(4) + (-t169 * t195 - t170 * t462 + (t169 * t310 + t531) * t233 + (t170 * t310 - t532) * t231) * m(3); (-m(4) - m(5) - m(6)) * (g(1) * t305 - g(2) * t306) + 0.2e1 * (t495 * t6 + t496 * t5) * m(6) + 0.2e1 * (t20 * t496 + t21 * t495) * m(5) + 0.2e1 * (t49 * t496 + t495 * t50) * m(4); t571 * t498 + t570 * t497 + (t549 * t310 + t546 * t308 + t526 * t198 + t525 * t197 + (t537 * t305 + t538 * t306) * qJD(4)) * t496 + (t548 * t310 + t547 * t308 + t528 * t198 + t527 * t197 + (t529 * t305 + t530 * t306) * qJD(4)) * t495 + (-t305 * t580 + t524 * t306) * t308 / 0.2e1 - ((t561 * t303 + t562 * t304) * t310 + (t563 * t303 + t304 * t567) * qJD(4)) * t310 / 0.2e1 + ((-t310 * t580 - t551) * t306 + (-t310 * t524 + t550) * t305) * t310 / 0.2e1 + t553 * t457 / 0.2e1 + t552 * t453 / 0.2e1 + ((-t424 * t519 + t508) * t305 + ((t305 * t520 + t511) * qJD(4) + t510) * t306) * t400 + ((t310 * t525 + t538) * t306 + (t310 * t526 + t537) * t305) * t399 + ((t310 * t527 + t530) * t306 + (t310 * t528 + t529) * t305) * t398 + ((-t423 * t520 - t508) * t306 + ((t306 * t519 + t511) * qJD(4) + t510) * t305) * t397 + (-(t358 + t482) * qJD(5) - (t439 * t46 - t440 * t45) * t310 - ((t439 * t52 - t45 * t569) * t306 + (t440 * t52 - t46 * t569) * t305) * qJD(4) + (-t5 * t434 + t45 * t441 + t7 * t442 + t52 * t478 + (-t434 * t46 + t443 * t52) * t310) * t306 + (-t6 * t434 + t46 * t441 + t7 * t443 + t52 * t477 + (t434 * t45 - t442 * t52) * t310) * t305 - g(1) * t513 - g(2) * t514 - g(3) * t569 - (-g(1) * t409 - t489 * t554) * t303) * m(6) + (-(t185 * t55 - t189 * t56) * t310 - (t79 * (-t185 * t305 - t189 * t306) + t382 * t229) * qJD(4) + (t161 * t197 - t163 * t198 + (t107 * t306 + t109 * t305) * qJD(4)) * t365 + t79 * ((t107 + t143) * t306 + (-t163 * t310 + t109) * t305) + t382 * t209 + ((-t310 * t56 - t20) * t306 + (-t21 + t483) * t305) * t226 + g(1) * t189 + g(2) * t185 - g(3) * t229) * m(5); (-(t358 + (t305 ^ 2 + t306 ^ 2) * t482) * qJD(4) + (qJD(4) * t383 + g(3) - t7) * t304 + (qJD(4) * t52 + t6 * t305 + t306 * t540 - t489) * t303) * m(6);];
tau = t1;
