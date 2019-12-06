% Calculate vector of inverse dynamics joint torques for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:32
% EndTime: 2019-12-05 15:09:09
% DurationCPUTime: 28.93s
% Computational Cost: add. (14829->648), mult. (21737->904), div. (0->0), fcn. (21124->6), ass. (0->327)
t569 = Icges(6,4) + Icges(5,5);
t568 = Icges(5,6) - Icges(6,6);
t552 = Icges(5,2) + Icges(6,3);
t592 = Icges(6,2) + Icges(5,3);
t281 = pkin(8) + qJ(3);
t279 = sin(t281);
t284 = sin(qJ(4));
t285 = cos(qJ(4));
t591 = (-t284 * t569 - t285 * t568) * t279;
t282 = sin(pkin(7));
t280 = cos(t281);
t410 = qJD(4) * t285;
t393 = t280 * t410;
t283 = cos(pkin(7));
t412 = qJD(4) * t283;
t417 = qJD(3) * t282;
t150 = -t282 * t393 + (t279 * t417 + t412) * t284;
t435 = t283 * t285;
t438 = t282 * t284;
t237 = t280 * t438 + t435;
t415 = qJD(3) * t285;
t397 = t279 * t415;
t151 = -qJD(4) * t237 - t282 * t397;
t396 = t280 * t417;
t590 = t150 * t568 + t151 * t569 + t396 * t592;
t436 = t283 * t284;
t399 = t279 * t436;
t411 = qJD(4) * t284;
t152 = qJD(3) * t399 - t282 * t411 - t283 * t393;
t437 = t282 * t285;
t239 = t280 * t436 - t437;
t153 = -qJD(4) * t239 - t283 * t397;
t416 = qJD(3) * t283;
t395 = t280 * t416;
t589 = t152 * t568 + t153 * t569 + t395 * t592;
t354 = Icges(5,5) * t285 - Icges(5,6) * t284;
t161 = Icges(5,3) * t279 + t280 * t354;
t357 = Icges(6,4) * t285 + Icges(6,6) * t284;
t163 = Icges(6,2) * t279 + t280 * t357;
t588 = -t591 * qJD(4) + (-t161 - t163) * qJD(3);
t238 = t280 * t437 - t436;
t448 = t279 * t282;
t100 = Icges(6,4) * t238 + Icges(6,2) * t448 + Icges(6,6) * t237;
t98 = Icges(5,5) * t238 - Icges(5,6) * t237 + Icges(5,3) * t448;
t559 = t98 + t100;
t240 = t280 * t435 + t438;
t447 = t279 * t283;
t101 = Icges(6,4) * t240 + Icges(6,2) * t447 + Icges(6,6) * t239;
t99 = Icges(5,5) * t240 - Icges(5,6) * t239 + Icges(5,3) * t447;
t558 = t99 + t101;
t587 = Icges(5,1) + Icges(6,1);
t586 = Icges(5,4) - Icges(6,5);
t453 = Icges(6,5) * t285;
t353 = Icges(6,3) * t284 + t453;
t159 = Icges(6,6) * t279 + t280 * t353;
t459 = Icges(5,4) * t285;
t358 = -Icges(5,2) * t284 + t459;
t165 = Icges(5,6) * t279 + t280 * t358;
t527 = -t159 + t165;
t160 = -Icges(5,3) * t280 + t279 * t354;
t162 = -Icges(6,2) * t280 + t279 * t357;
t573 = -t162 - t160;
t454 = Icges(6,5) * t284;
t361 = Icges(6,1) * t285 + t454;
t167 = Icges(6,4) * t279 + t280 * t361;
t460 = Icges(5,4) * t284;
t362 = Icges(5,1) * t285 - t460;
t169 = Icges(5,5) * t279 + t280 * t362;
t526 = -t167 - t169;
t585 = (-t285 * t552 + t454 - t460) * t279;
t584 = -t150 * t552 - t151 * t586 - t396 * t568;
t583 = -t152 * t552 - t153 * t586 - t395 * t568;
t582 = t150 * t586 + t151 * t587 + t396 * t569;
t581 = t152 * t586 + t153 * t587 + t395 * t569;
t580 = qJD(3) * t527 + qJD(4) * t585;
t231 = (-Icges(5,1) * t284 - t459) * t279;
t414 = qJD(4) * t279;
t579 = -(-Icges(6,1) * t284 + t453) * t414 - qJD(4) * t231 + t526 * qJD(3);
t462 = Icges(5,4) * t238;
t102 = -Icges(5,2) * t237 + Icges(5,6) * t448 + t462;
t212 = Icges(6,5) * t238;
t96 = Icges(6,6) * t448 + Icges(6,3) * t237 + t212;
t578 = t102 - t96;
t461 = Icges(5,4) * t240;
t103 = -Icges(5,2) * t239 + Icges(5,6) * t447 + t461;
t213 = Icges(6,5) * t240;
t97 = Icges(6,6) * t447 + Icges(6,3) * t239 + t213;
t577 = t103 - t97;
t456 = Icges(6,5) * t237;
t104 = Icges(6,1) * t238 + Icges(6,4) * t448 + t456;
t214 = Icges(5,4) * t237;
t106 = Icges(5,1) * t238 + Icges(5,5) * t448 - t214;
t576 = t104 + t106;
t455 = Icges(6,5) * t239;
t105 = Icges(6,1) * t240 + Icges(6,4) * t447 + t455;
t215 = Icges(5,4) * t239;
t107 = Icges(5,1) * t240 + Icges(5,5) * t447 - t215;
t575 = t105 + t107;
t445 = t279 * t285;
t264 = Icges(6,5) * t445;
t446 = t279 * t284;
t450 = Icges(6,6) * t280;
t158 = Icges(6,3) * t446 + t264 - t450;
t164 = -Icges(5,6) * t280 + t279 * t358;
t574 = t158 - t164;
t166 = -Icges(6,4) * t280 + t279 * t361;
t168 = -Icges(5,5) * t280 + t279 * t362;
t550 = t166 + t168;
t418 = qJD(3) * t280;
t572 = t279 * t588 + t418 * t573;
t571 = t279 * t589 + t418 * t558;
t570 = t279 * t590 + t418 * t559;
t491 = t283 ^ 2;
t492 = t282 ^ 2;
t516 = t491 + t492;
t355 = -Icges(4,5) * t279 - Icges(4,6) * t280;
t541 = t150 * t578 + t151 * t576 + t237 * t584 + t238 * t582 + t282 * t570;
t540 = t150 * t577 + t151 * t575 + t237 * t583 + t238 * t581 + t282 * t571;
t539 = t152 * t578 + t153 * t576 + t239 * t584 + t240 * t582 + t283 * t570;
t538 = t152 * t577 + t153 * t575 + t239 * t583 + t240 * t581 + t283 * t571;
t567 = t150 * t574 - t151 * t550 + t237 * t580 + t238 * t579 + t282 * t572;
t566 = t152 * t574 - t153 * t550 + t239 * t580 + t240 * t579 + t283 * t572;
t555 = -t237 * t578 + t238 * t576 + t448 * t559;
t565 = -t237 * t577 + t238 * t575 + t448 * t558;
t564 = -t239 * t578 + t240 * t576 + t447 * t559;
t535 = -t239 * t577 + t240 * t575 + t447 * t558;
t366 = t104 * t285 + t284 * t96;
t47 = -t100 * t280 + t279 * t366;
t351 = -t102 * t284 + t106 * t285;
t49 = t279 * t351 - t280 * t98;
t534 = t47 + t49;
t365 = t105 * t285 + t284 * t97;
t48 = -t101 * t280 + t279 * t365;
t350 = -t103 * t284 + t107 * t285;
t50 = t279 * t350 - t280 * t99;
t533 = t48 + t50;
t563 = t237 * t574 + t238 * t550 - t448 * t573;
t562 = t239 * t574 + t240 * t550 - t447 * t573;
t348 = t158 * t284 + t166 * t285;
t443 = t280 * t162;
t62 = t279 * t348 - t443;
t347 = -t164 * t284 + t168 * t285;
t444 = t280 * t160;
t63 = t279 * t347 - t444;
t507 = t62 + t63;
t532 = rSges(6,3) + qJ(5);
t544 = rSges(6,1) + pkin(4);
t422 = t445 * t532 - t446 * t544;
t554 = t282 * t283;
t244 = t279 * t412 + t417;
t245 = t282 * t414 - t416;
t413 = qJD(4) * t280;
t549 = t591 * t413 + (t237 * t569 + t238 * t568) * t245 + (t239 * t569 + t240 * t568) * t244;
t548 = 0.2e1 * qJD(3);
t547 = 2 * qJDD(3);
t406 = qJD(3) * qJD(4);
t315 = qJDD(4) * t279 + t280 * t406;
t404 = qJDD(3) * t282;
t154 = t283 * t315 + t404;
t403 = qJDD(3) * t283;
t155 = t282 * t315 - t403;
t243 = -qJDD(4) * t280 + t279 * t406;
t546 = t154 * t565 + t155 * t555 + t243 * t563 + t244 * t540 + t245 * t541 + t413 * t567;
t545 = t154 * t535 + t155 * t564 + t243 * t562 + t244 * t538 + t245 * t539 + t413 * t566;
t510 = ((t351 + t366) * qJD(3) - t590) * t280 + (t582 * t285 + t584 * t284 + (-t284 * t576 - t285 * t578) * qJD(4) + t559 * qJD(3)) * t279;
t509 = ((t350 + t365) * qJD(3) - t589) * t280 + (t581 * t285 + t583 * t284 + (-t284 * t575 - t285 * t577) * qJD(4) + t558 * qJD(3)) * t279;
t543 = t244 * t565 + t245 * t555 - t413 * t563;
t542 = t244 * t535 + t245 * t564 - t413 * t562;
t537 = t244 * t533 + t245 * t534 - t413 * t507;
t300 = -t279 * t353 + t450;
t134 = t300 * t282;
t140 = t164 * t282;
t531 = t134 + t140;
t135 = t300 * t283;
t141 = t164 * t283;
t530 = t135 + t141;
t142 = t166 * t282;
t144 = t168 * t282;
t529 = -t142 - t144;
t143 = t166 * t283;
t145 = t168 * t283;
t528 = -t143 - t145;
t327 = t161 - t347;
t328 = -t163 + t348;
t493 = (t162 * t283 + t365) * t244 + (t162 * t282 + t366) * t245;
t494 = -(-t160 * t283 - t350) * t244 - (-t160 * t282 - t351) * t245;
t525 = (-t493 - t494 + (-t327 + t328) * t413) * t279;
t376 = rSges(6,1) * t285 + rSges(6,3) * t284;
t524 = (pkin(4) * t285 + qJ(5) * t284 + t376) * t279;
t523 = t244 * t558 + t245 * t559;
t248 = pkin(3) * t279 - pkin(6) * t280;
t326 = qJD(3) * t248;
t204 = t282 * t326;
t205 = t283 * t326;
t442 = t280 * t282;
t261 = pkin(6) * t442;
t441 = t280 * t283;
t263 = pkin(6) * t441;
t522 = -t282 * t204 - t283 * t205 - (-pkin(3) * t448 + t261) * t417 - (-pkin(3) * t447 + t263) * t416;
t521 = -t443 - t444;
t520 = t562 * t279;
t519 = t563 * t279;
t518 = t565 * t283;
t517 = t564 * t282;
t503 = pkin(3) * t280 + pkin(6) * t279;
t242 = t503 * qJD(3);
t278 = qJD(2) * t282;
t378 = -t248 * t416 + t278;
t408 = qJD(5) * t239;
t427 = -rSges(6,2) * t280 + t524;
t434 = rSges(6,2) * t448 + t237 * t532 + t238 * t544;
t36 = t245 * t427 + t413 * t434 + t378 + t408;
t419 = qJD(2) * t283;
t323 = -t248 * t417 - t419;
t409 = qJD(5) * t237;
t433 = rSges(6,2) * t447 + t239 * t532 + t240 * t544;
t37 = -t244 * t427 - t413 * t433 + t323 + t409;
t511 = t244 * t37 - t245 * t36 - g(3);
t508 = ((-t347 - t348) * qJD(3) - t588) * t280 + (t579 * t285 + t580 * t284 + (t284 * t550 - t285 * t574) * qJD(4) + t573 * qJD(3)) * t279;
t502 = (t550 + t585) * t413 + (t238 * t552 + t214 - t456 - t576) * t245 + (t240 * t552 + t215 - t455 - t575) * t244;
t501 = (Icges(6,1) * t446 - t231 - t264 - t574) * t413 + (-t237 * t587 + t212 - t462 - t578) * t245 + (-t239 * t587 + t213 - t461 - t577) * t244;
t500 = t549 * t279;
t367 = t282 * t49 + t283 * t50;
t368 = t282 * t47 + t283 * t48;
t499 = t367 + t368;
t498 = t283 * t535 + t517;
t497 = t282 * t555 + t518;
t496 = t279 * (g(1) * t283 + g(2) * t282);
t407 = qJD(5) * t284;
t259 = t279 * t407;
t224 = t503 * t282;
t225 = t503 * t283;
t390 = t224 * t417 + t225 * t416 + qJD(1);
t35 = t244 * t434 - t245 * t433 + t259 + t390;
t476 = rSges(6,2) * t395 - t152 * t532 + t153 * t544 + t408;
t394 = t284 * t418;
t311 = t279 * t410 + t394;
t341 = -t204 * t417 - t205 * t416 + t224 * t404 + t225 * t403 + qJDD(1);
t477 = rSges(6,2) * t396 - t150 * t532 + t151 * t544 + t409;
t5 = qJD(5) * t311 + qJDD(5) * t446 + t154 * t434 - t155 * t433 + t244 * t477 - t245 * t476 + t341;
t495 = t35 * t476 + t433 * t5;
t490 = t154 / 0.2e1;
t489 = t155 / 0.2e1;
t488 = t243 / 0.2e1;
t487 = -t244 / 0.2e1;
t486 = t244 / 0.2e1;
t485 = -t245 / 0.2e1;
t484 = t245 / 0.2e1;
t272 = t279 * rSges(6,2);
t475 = t259 + t311 * qJ(5) + (-t279 * t411 + t280 * t415) * pkin(4) + (-rSges(6,1) * t284 + rSges(6,3) * t285) * t414 + (t280 * t376 + t272) * qJD(3);
t474 = rSges(5,1) * t285;
t271 = t279 * rSges(5,3);
t234 = (-rSges(5,1) * t284 - rSges(5,2) * t285) * t279;
t377 = -rSges(5,2) * t284 + t474;
t95 = qJD(4) * t234 + (t280 * t377 + t271) * qJD(3);
t465 = -t242 - t95;
t440 = t280 * t284;
t439 = t280 * t285;
t432 = -t237 * t544 + t238 * t532;
t431 = t239 * t544 - t240 * t532;
t256 = rSges(6,2) * t442;
t430 = -t282 * t524 + t256;
t258 = rSges(6,2) * t441;
t429 = t283 * t524 - t258;
t173 = -rSges(5,3) * t280 + t279 * t377;
t426 = -t173 - t248;
t423 = t224 * t282 + t225 * t283;
t421 = rSges(5,2) * t279 * t438 + rSges(5,3) * t442;
t420 = rSges(5,2) * t399 + rSges(5,3) * t441;
t405 = qJDD(2) * t283;
t402 = -t242 - t475;
t401 = rSges(5,1) * t445;
t400 = -m(3) - m(4) - m(5) - m(6);
t398 = -t248 - t427;
t387 = t416 / 0.2e1;
t385 = -t413 / 0.2e1;
t384 = t413 / 0.2e1;
t247 = rSges(4,1) * t280 - rSges(4,2) * t279;
t246 = rSges(4,1) * t279 + rSges(4,2) * t280;
t109 = rSges(5,1) * t238 - rSges(5,2) * t237 + rSges(5,3) * t448;
t277 = qJDD(2) * t282;
t339 = -qJD(3) * t242 - qJDD(3) * t248;
t307 = t283 * t339 + t277;
t83 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t396;
t33 = -t109 * t243 + t155 * t173 + t245 * t95 + t413 * t83 + t307;
t111 = rSges(5,1) * t240 - rSges(5,2) * t239 + rSges(5,3) * t447;
t296 = t282 * t339 - t405;
t85 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t395;
t34 = t111 * t243 - t154 * t173 - t244 * t95 - t413 * t85 + t296;
t373 = t282 * t33 - t283 * t34;
t356 = Icges(4,5) * t280 - Icges(4,6) * t279;
t349 = t109 * t283 - t111 * t282;
t344 = -(-t246 * t416 + t278) * t283 - (-t246 * t417 - t419) * t282;
t343 = t516 * t247;
t342 = t516 * qJD(3) * t246;
t175 = rSges(5,1) * t439 - rSges(5,2) * t440 + t271;
t241 = t247 * qJD(3);
t340 = -qJD(3) * t241 - qJDD(3) * t246;
t46 = t109 * t244 - t111 * t245 + t390;
t320 = t46 * t349;
t312 = qJD(3) * t355;
t310 = t35 * t477 + t434 * t5;
t309 = -t36 * t434 + t37 * t433;
t289 = (t35 * t434 - t37 * t427) * t283 + (-t35 * t433 + t36 * t427) * t282;
t288 = t516 * t355;
t217 = t246 * t283;
t216 = t246 * t282;
t207 = t355 * t283;
t206 = t355 * t282;
t197 = t283 * t312;
t196 = t282 * t312;
t177 = Icges(4,3) * t282 + t283 * t356;
t176 = -Icges(4,3) * t283 + t282 * t356;
t149 = -t283 * t401 + t420;
t147 = -t282 * t401 + t421;
t132 = -rSges(5,1) * t239 - rSges(5,2) * t240;
t128 = -rSges(5,1) * t237 - rSges(5,2) * t238;
t113 = t282 * t340 - t405;
t112 = t283 * t340 + t277;
t86 = qJD(3) * t343 + qJD(1);
t61 = -qJD(3) * t342 + qJDD(3) * t343 + qJDD(1);
t60 = -t111 * t413 - t173 * t244 + t323;
t59 = t109 * t413 + t173 * t245 + t378;
t26 = t109 * t154 - t111 * t155 + t244 * t83 - t245 * t85 + t341;
t7 = -qJD(5) * t150 + qJDD(5) * t237 - t154 * t427 + t243 * t433 - t244 * t475 - t413 * t476 + t296;
t6 = -qJD(5) * t152 + qJDD(5) * t239 + t155 * t427 - t243 * t434 + t245 * t475 + t413 * t477 + t307;
t1 = [(m(2) + m(3)) * qJDD(1) + m(4) * t61 + m(5) * t26 + m(6) * t5 + (-m(2) + t400) * g(3); t400 * (g(1) * t282 - g(2) * t283) + m(4) * (t112 * t282 - t113 * t283) + m(5) * t373 + m(6) * (t282 * t6 - t283 * t7) + m(3) * t516 * qJDD(2); (t206 * qJD(3) * t491 + (-t207 * t283 + t288) * t417) * t387 - (t207 * qJD(3) * t492 + (-t206 * t282 + t288) * t416) * t417 / 0.2e1 + (t535 * t282 - t283 * t564) * t490 + (t282 * t565 - t283 * t555) * t489 + (t282 * t533 - t283 * t534) * t488 + (((t239 * t527 + t240 * t526 + t517) * t280 + t520) * qJD(4) + (((t521 + t535) * qJD(4) + t523) * t280 + t525) * t283 + (t239 * t531 + t240 * t529) * t245 + (t239 * t530 + t240 * t528) * t244) * t487 + (t282 * t538 - t283 * t539) * t486 + (((t237 * t527 + t238 * t526 + t518) * t280 + t519) * qJD(4) + (((t521 + t555) * qJD(4) + t523) * t280 + t525) * t282 + (t237 * t531 + t238 * t529) * t245 + (t237 * t530 + t238 * t528) * t244) * t485 + (t282 * t540 - t283 * t541) * t484 - t537 * t414 / 0.2e1 + (((t135 * t284 - t143 * t285 + t101) * t244 + (t134 * t284 - t142 * t285 + t100) * t245 + t62 * qJD(4)) * t279 + ((-t328 * t280 + (-t159 * t284 - t167 * t285 - t162) * t279 + t368) * qJD(4) + t493) * t280 + ((t141 * t284 - t145 * t285 + t99) * t244 + (t140 * t284 - t144 * t285 + t98) * t245 + t63 * qJD(4)) * t279 + ((t327 * t280 + (t165 * t284 - t169 * t285 - t160) * t279 + t367) * qJD(4) + t494) * t280) * t384 + (-(t309 * t279 + (t36 * t430 + t37 * t429 + t289) * t280) * qJD(4) - (t282 * t37 + t283 * t36) * (-t259 - t242) - g(1) * (t258 + t263) - g(2) * (t256 + t261) - g(3) * t503 - (-t284 * t532 - t285 * t544 - pkin(3)) * t496 + t5 * t423 + (t36 * t402 + t398 * t6 + t495) * t283 + (t37 * t402 + t398 * t7 + t310) * t282 + t511 * (t439 * t544 + t440 * t532 + t272) + (-t244 * t430 - t245 * t429 - t280 * t407 + t522) * t35) * m(6) + (-t59 * (t175 * t245 - t416 * t503) - t60 * (-t175 * t244 - t417 * t503) - ((-t109 * t59 + t111 * t60) * t279 + (t59 * (t173 * t282 + t147) + t60 * (-t173 * t283 - t149) + t320) * t280) * qJD(4) - g(1) * (t263 + t420) - g(2) * (t261 + t421) - g(3) * (t175 + t503) - (-pkin(3) - t474) * t496 + t26 * t423 + (t111 * t26 + t33 * t426 + t465 * t59) * t283 + (t109 * t26 + t34 * t426 + t465 * t60) * t282 + (-t147 * t244 + t149 * t245 + t282 * t83 + t283 * t85 + t522) * t46) * m(5) + (g(1) * t217 + g(2) * t216 - g(3) * t247 - (t86 * (-t216 * t282 - t217 * t283) + t344 * t247) * qJD(3) + t61 * t343 - t86 * t342 + (-t112 * t283 - t113 * t282) * t246 + t344 * t241) * m(4) + ((-t176 * t554 + t177 * t492) * t547 + (-t196 * t554 + t197 * t492) * t548 + t545) * t282 / 0.2e1 - ((t176 * t491 - t177 * t554) * t547 + (t196 * t491 - t197 * t554) * t548 + t546) * t283 / 0.2e1 + ((-t510 + t542) * t283 + (t509 + t543) * t282) * t385; (t498 * t279 - t280 * t562) * t490 + (t497 * t279 - t280 * t563) * t489 + (t279 * t499 - t280 * t507) * t488 + (t239 * t502 + t240 * t501 - t283 * t500) * t487 + (t566 * t280 + (t282 * t539 + t283 * t538) * t279 + (t280 * t498 + t520) * qJD(3)) * t486 + (t237 * t502 + t238 * t501 - t282 * t500) * t485 + (t567 * t280 + (t282 * t541 + t283 * t540) * t279 + (t280 * t497 + t519) * qJD(3)) * t484 - (t154 * t533 + t155 * t534 + t507 * t243 + t509 * t244 + t510 * t245 + t508 * t413) * t280 / 0.2e1 + t546 * t448 / 0.2e1 + t545 * t447 / 0.2e1 + t537 * qJD(3) * t279 / 0.2e1 + (t508 * t280 + (t282 * t510 + t283 * t509) * t279 + (t279 * t507 + t280 * t499) * qJD(3)) * t385 + (t549 * t280 + (t284 * t502 + t285 * t501) * t279) * t384 + t543 * t396 / 0.2e1 + t542 * t280 * t387 + (-(t238 * t37 + t240 * t36 + t35 * t445) * qJD(5) - (t35 * t431 + t36 * t422) * t245 - (t35 * t432 - t37 * t422) * t244 - (t36 * t432 + t37 * t431) * t413 + g(1) * t431 - g(2) * t432 - g(3) * t422 + (qJD(3) * t289 + t36 * t477 - t37 * t476 - t433 * t7 + t434 * t6) * t280 + (t309 * qJD(3) + (-t37 * t475 - t427 * t7 + t310) * t283 + (t36 * t475 + t427 * t6 - t495) * t282) * t279) * m(6) + (-t59 * (t128 * t413 + t234 * t245) - t60 * (-t132 * t413 - t234 * t244) - t46 * (t128 * t244 - t132 * t245) + (t33 * t109 - t34 * t111 + t59 * t83 - t60 * t85 + (t320 + (t282 * t59 - t283 * t60) * t173) * qJD(3)) * t280 + (t59 * (-qJD(3) * t109 + t282 * t95) + t60 * (qJD(3) * t111 - t283 * t95) + t26 * t349 + t46 * (-t282 * t85 + t283 * t83) + t373 * t173) * t279 - g(1) * t132 - g(2) * t128 - g(3) * t234) * m(5); (t35 * t394 - t37 * t150 - t36 * t152 + (t284 * t5 + t35 * t410) * t279 + t511 * t446 + (t245 * t35 + t37 * t413 - g(1) + t6) * t239 + (-t244 * t35 - t36 * t413 - g(2) + t7) * t237) * m(6);];
tau = t1;
