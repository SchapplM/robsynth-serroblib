% Calculate vector of inverse dynamics joint torques for
% S5RRPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:55
% DurationCPUTime: 14.60s
% Computational Cost: add. (10147->556), mult. (10823->643), div. (0->0), fcn. (8244->6), ass. (0->299)
t307 = qJ(1) + qJ(2);
t298 = sin(t307);
t299 = cos(t307);
t310 = cos(qJ(4));
t452 = t299 * t310;
t308 = sin(qJ(4));
t453 = t299 * t308;
t142 = Icges(6,5) * t453 - Icges(6,6) * t298 - Icges(6,3) * t452;
t476 = Icges(5,4) * t308;
t369 = Icges(5,2) * t310 + t476;
t344 = t369 * t299;
t148 = -Icges(5,6) * t298 + t344;
t601 = t142 - t148;
t265 = Icges(6,5) * t452;
t150 = Icges(6,1) * t453 - Icges(6,4) * t298 - t265;
t475 = Icges(5,4) * t310;
t371 = Icges(5,1) * t308 + t475;
t152 = -Icges(5,5) * t298 + t299 * t371;
t600 = t150 + t152;
t367 = Icges(5,5) * t308 + Icges(5,6) * t310;
t342 = t367 * t298;
t143 = Icges(5,3) * t299 + t342;
t368 = Icges(6,4) * t308 - Icges(6,6) * t310;
t145 = Icges(6,2) * t299 + t298 * t368;
t605 = t143 + t145;
t243 = -Icges(5,2) * t308 + t475;
t301 = Icges(6,5) * t310;
t528 = Icges(6,3) * t308 + t301;
t604 = -t243 + t528;
t472 = Icges(6,5) * t308;
t245 = Icges(6,1) * t310 + t472;
t247 = Icges(5,1) * t310 - t476;
t603 = t245 + t247;
t584 = t600 * t308 - t601 * t310;
t370 = Icges(6,1) * t308 - t301;
t602 = t370 + t371;
t341 = t367 * t299;
t144 = -Icges(5,3) * t298 + t341;
t146 = Icges(6,4) * t453 - Icges(6,2) * t298 - Icges(6,6) * t452;
t585 = -t144 - t146;
t456 = t298 * t308;
t264 = Icges(6,5) * t456;
t455 = t298 * t310;
t141 = Icges(6,6) * t299 - Icges(6,3) * t455 + t264;
t559 = -t141 * t452 - t605 * t298;
t239 = Icges(5,5) * t310 - Icges(5,6) * t308;
t241 = Icges(6,4) * t310 + Icges(6,6) * t308;
t599 = t239 + t241;
t306 = qJD(1) + qJD(2);
t598 = t602 * t306;
t355 = t243 * t310 + t247 * t308;
t591 = t245 * t308 - t310 * t528 + t355;
t597 = (Icges(5,6) - Icges(6,6)) * t306 + t604 * qJD(4);
t596 = (-Icges(6,4) - Icges(5,5)) * t306 + t603 * qJD(4);
t595 = t584 * t299;
t149 = Icges(6,4) * t299 + t298 * t370;
t449 = t299 * t145 + t149 * t456;
t53 = -t141 * t455 + t449;
t345 = t369 * t298;
t147 = Icges(5,6) * t299 + t345;
t266 = Icges(5,4) * t455;
t151 = Icges(5,1) * t456 + Icges(5,5) * t299 + t266;
t55 = t299 * t143 + t147 * t455 + t151 * t456;
t542 = t55 + t53;
t541 = t585 * t299 + t601 * t455 - t600 * t456;
t361 = t147 * t310 + t151 * t308;
t540 = -t149 * t453 - t299 * t361 - t559;
t539 = t585 * t298 + t595;
t365 = -Icges(6,3) * t310 + t472;
t594 = (t365 - t369) * qJD(4);
t593 = t602 * qJD(4);
t592 = t604 * t308 + t603 * t310;
t558 = t149 * t308 + t361;
t300 = t310 * qJ(5);
t489 = rSges(6,1) + pkin(4);
t523 = t308 * t489;
t427 = t310 * rSges(6,3) + t300 - t523;
t424 = -t298 * rSges(6,2) - rSges(6,3) * t452;
t438 = -t299 * t300 + t453 * t489 + t424;
t564 = rSges(6,3) + qJ(5);
t426 = t308 * t564 + t310 * t489;
t340 = t365 * t306;
t590 = t298 * t340 - t299 * t597 - t306 * t345;
t589 = t298 * t597 + t299 * t340 - t306 * t344;
t588 = t298 * t598 - t299 * t596;
t587 = t298 * t596 + t299 * t598;
t586 = t601 * t308 + t600 * t310;
t538 = (t149 + t151) * t310 + (t141 - t147) * t308;
t458 = t241 * t299;
t460 = t239 * t299;
t561 = t298 * t591 + t458 + t460;
t164 = t298 * t239;
t166 = t298 * t241;
t560 = t245 * t453 + t299 * t355 - t452 * t528 - t164 - t166;
t464 = t141 * t310;
t583 = t464 - t558;
t582 = (-Icges(6,2) - Icges(5,3)) * t306 + t599 * qJD(4);
t581 = t591 * t306 + (-t367 - t368) * qJD(4);
t580 = t298 * t539 + t299 * t540;
t579 = t541 * t298 + t542 * t299;
t577 = t592 * qJD(4) - t306 * t599 - t593 * t308 + t594 * t310;
t576 = t561 * t306;
t533 = t368 * t306;
t575 = t299 * t533 + (t341 - t583) * t306 + t582 * t298;
t574 = t298 * t533 + (t342 - t584) * t306 - t582 * t299;
t573 = t560 * t306;
t572 = -t538 * qJD(4) + t605 * t306 - t587 * t308 + t589 * t310;
t571 = t586 * qJD(4) + t585 * t306 - t588 * t308 + t590 * t310;
t570 = qJD(4) * t579 + t576;
t569 = qJD(4) * t580 - t573;
t568 = t583 * qJD(4) + t589 * t308 + t587 * t310;
t567 = t584 * qJD(4) + t590 * t308 + t588 * t310;
t566 = t298 * t581 - t299 * t577;
t565 = t298 * t577 + t299 * t581;
t431 = t528 - t370;
t432 = -t365 - t245;
t563 = (t308 * t431 - t310 * t432) * t306;
t429 = t243 + t371;
t430 = -t369 + t247;
t562 = (t308 * t429 - t310 * t430) * t306;
t289 = t299 * rSges(5,3);
t154 = rSges(5,1) * t456 + rSges(5,2) * t455 + t289;
t293 = t299 * pkin(2);
t202 = t298 * qJ(3) + t293;
t292 = t299 * pkin(7);
t530 = t292 + t202;
t117 = t154 + t530;
t531 = t299 * rSges(6,2) + t456 * t489;
t417 = qJD(4) * t310;
t556 = t298 * t417 + t306 * t453;
t309 = sin(qJ(1));
t482 = pkin(1) * qJD(1);
t411 = t309 * t482;
t201 = rSges(3,1) * t298 + rSges(3,2) * t299;
t462 = t201 * t306;
t157 = -t411 - t462;
t277 = qJD(3) * t299;
t138 = t202 * t306 - t277;
t419 = qJD(4) * t306;
t183 = qJDD(4) * t298 + t299 * t419;
t279 = t299 * qJ(3);
t199 = pkin(2) * t298 - t279;
t305 = qJDD(1) + qJDD(2);
t304 = t306 ^ 2;
t311 = cos(qJ(1));
t312 = qJD(1) ^ 2;
t339 = (-qJDD(1) * t309 - t311 * t312) * pkin(1);
t422 = qJD(3) * t306;
t322 = qJDD(3) * t298 + t299 * t422 + t339;
t317 = -t292 * t304 + t322;
t416 = qJD(5) * t310;
t398 = t299 * t416;
t400 = t299 * t417;
t418 = qJD(4) * t308;
t401 = t299 * t418;
t451 = t306 * t310;
t526 = -t398 + t564 * (t298 * t451 + t401) + t489 * t400;
t478 = t531 * t306 - t526;
t297 = qJD(5) * t308;
t435 = qJD(4) * t427 + t297;
t507 = -qJD(4) * (t297 + t435) + qJDD(5) * t310;
t6 = t426 * t183 + (-t199 + t438) * t305 + (-t138 - t398 - t478) * t306 + (-pkin(7) * t305 - t507) * t298 + t317;
t555 = t6 - g(1);
t184 = qJDD(4) * t299 - t298 * t419;
t303 = t311 * pkin(1);
t488 = pkin(1) * t309;
t387 = qJDD(1) * t303 - t312 * t488;
t454 = t299 * t306;
t235 = qJ(3) * t454;
t276 = qJD(3) * t298;
t433 = t235 + t276;
t457 = t298 * t306;
t335 = t306 * (-pkin(2) * t457 + t433) + t305 * t202 + t298 * t422 + t387;
t486 = pkin(7) * t298;
t319 = t305 * t292 - t304 * t486 + t335;
t399 = t298 * t416;
t439 = -rSges(6,3) * t455 - t298 * t300 + t531;
t403 = t298 * t418;
t525 = t564 * t403 + t489 * t556;
t450 = -t424 * t306 - (-qJ(5) * t454 - qJD(5) * t298) * t310 - t525;
t7 = t439 * t305 - t426 * t184 + (-t399 - t450) * t306 + (-qJDD(3) + t507) * t299 + t319;
t554 = t7 - g(2);
t553 = t298 * t575 + t299 * t572;
t552 = t298 * t574 + t299 * t571;
t378 = rSges(5,1) * t308 + rSges(5,2) * t310;
t218 = t378 * qJD(4);
t258 = rSges(5,1) * t310 - rSges(5,2) * t308;
t156 = -t298 * rSges(5,3) + t299 * t378;
t391 = -t199 - t486;
t385 = t156 + t391;
t421 = qJD(4) * t298;
t350 = -rSges(5,1) * t400 + rSges(5,2) * t401;
t99 = (t298 * t378 + t289) * t306 + t350;
t24 = -t218 * t421 + t183 * t258 + (-t138 - t99) * t306 + t385 * t305 + t317;
t551 = t24 - g(1);
t406 = t299 * rSges(5,2) * t451 + t556 * rSges(5,1);
t101 = (-rSges(5,2) * t418 - rSges(5,3) * t306) * t298 + t406;
t25 = t101 * t306 + t154 * t305 - t184 * t258 + (qJD(4) * t218 - qJDD(3)) * t299 + t319;
t550 = t25 - g(2);
t251 = rSges(4,2) * t454;
t483 = rSges(4,3) * t299;
t200 = rSges(4,2) * t298 + t483;
t434 = -t199 + t200;
t50 = (-rSges(4,3) * t457 - t138 + t251) * t306 + t434 * t305 + t322;
t549 = t50 - g(1);
t203 = -rSges(4,2) * t299 + t298 * rSges(4,3);
t428 = rSges(4,2) * t457 + rSges(4,3) * t454;
t51 = -qJDD(3) * t299 + t305 * t203 + t306 * t428 + t335;
t548 = t51 - g(2);
t547 = pkin(7) * t306;
t161 = rSges(3,1) * t454 - rSges(3,2) * t457;
t546 = -t161 * t306 - t201 * t305 - g(1) + t339;
t204 = t299 * rSges(3,1) - rSges(3,2) * t298;
t545 = t204 * t305 - t306 * t462 - g(2) + t387;
t544 = -t298 * t572 + t299 * t575;
t543 = -t298 * t571 + t299 * t574;
t140 = t202 + t203;
t536 = t140 * t306;
t532 = t426 * t421;
t137 = t306 * t156;
t527 = -rSges(5,2) * t403 - t137 + t406 + t433;
t524 = -t306 * t200 + t428;
t380 = t276 - t399;
t522 = -t438 * t306 + t235 + t380 + t525;
t420 = qJD(4) * t299;
t519 = t306 * t117 - t258 * t420;
t518 = t306 * (t530 + t439) - t277 + t398 - t426 * t420;
t440 = t243 * t299 + t152;
t444 = -t247 * t299 + t148;
t504 = t308 * t444 - t310 * t440;
t441 = -Icges(5,2) * t456 + t151 + t266;
t445 = -t247 * t298 + t147;
t503 = t308 * t445 - t310 * t441;
t442 = -Icges(6,3) * t453 + t150 - t265;
t446 = t245 * t299 + t142;
t502 = -t308 * t446 - t310 * t442;
t443 = -t298 * t528 + t149;
t447 = Icges(6,1) * t455 + t141 + t264;
t501 = -t308 * t447 - t310 * t443;
t500 = t298 ^ 2;
t499 = -pkin(2) - pkin(7);
t497 = t183 / 0.2e1;
t496 = t184 / 0.2e1;
t495 = t298 / 0.2e1;
t494 = -t299 / 0.2e1;
t485 = g(2) * t299;
t412 = t311 * t482;
t48 = t412 + t518;
t481 = t298 * t48;
t52 = t297 + (-t298 * t439 - t299 * t438) * qJD(4);
t480 = t52 * t310;
t461 = t367 * t306;
t437 = t489 * t455 + t564 * t456;
t436 = t426 * t299;
t188 = t306 * t199;
t423 = t276 - t188;
t414 = -rSges(6,2) + t499;
t413 = -rSges(5,3) + t499;
t395 = -t421 / 0.2e1;
t394 = t421 / 0.2e1;
t393 = -t420 / 0.2e1;
t392 = t420 / 0.2e1;
t389 = t310 * t564;
t388 = t146 - t464;
t384 = t276 - t411;
t383 = -t277 + t412;
t382 = t298 * t499 + t279;
t381 = g(1) * t298 - t485;
t259 = rSges(2,1) * t311 - rSges(2,2) * t309;
t255 = rSges(2,1) * t309 + rSges(2,2) * t311;
t338 = t380 + t532;
t47 = -t411 + (t391 + t438) * t306 + t338;
t377 = t298 * t47 - t299 * t48;
t193 = t258 * t421;
t351 = t193 + t384;
t63 = t306 * t385 + t351;
t64 = t383 + t519;
t372 = t298 * t63 - t299 * t64;
t358 = -t154 * t298 - t156 * t299;
t352 = -t188 + t384;
t349 = t377 * t308;
t348 = t277 - t350;
t139 = t483 + t279 + (rSges(4,2) - pkin(2)) * t298;
t320 = t277 + t526;
t79 = -t298 * t389 + t530 + t531;
t78 = t382 + t438;
t116 = t156 + t382;
t112 = t306 * t434 + t384;
t113 = t383 + t536;
t316 = (-t112 * t293 + (t112 * (-rSges(4,3) - qJ(3)) - t113 * pkin(2)) * t298) * t306;
t315 = (t63 * t413 * t299 + (t63 * (-qJ(3) - t378) + t64 * t413) * t298) * t306;
t314 = ((-t389 * t48 + t414 * t47) * t299 + (t47 * (-qJ(3) - t523) + t48 * t414) * t298) * t306;
t313 = -t560 * t183 / 0.2e1 - t586 * t497 + (((t449 + t55 + t539 - t595) * t299 + ((t388 + t144 - t558) * t299 - t559 - t540 + t541) * t298) * qJD(4) + t576) * t395 + (-qJD(4) * t591 - t308 * t594 - t310 * t593) * t306 + (t561 + t538) * t496 + (t565 + t568) * t392 + (Icges(4,1) + Icges(3,3) + t592) * t305 + (((t298 * t388 + t449 - t53) * t298 + t144 * t500 + ((t558 - t585) * t299 + t559 + t541) * t299) * qJD(4) + t569 + t573) * t393 + (t566 + t567 + t570) * t394;
t181 = t258 * t299;
t177 = t258 * t298;
t158 = t204 * t306 + t412;
t77 = t358 * qJD(4);
t5 = qJDD(5) * t308 - t438 * t184 - t439 * t183 + (t298 * t450 + t299 * t478 + t416) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t313 + (t545 * (t204 + t303) + t546 * (-t201 - t488) + (-t161 - t412 + t158) * t157) * m(3) + ((t255 ^ 2 + t259 ^ 2) * qJDD(1) + g(1) * t255 - g(2) * t259) * m(2) + (t47 * (t320 - t412) + t314 + t554 * (t303 + t79) + t555 * (t78 - t488) + (-t411 + t47 - (-t416 - t547) * t298 - t352 + t522 - t532) * t48) * m(6) + (t63 * (t348 - t412) + t315 + (pkin(7) * t457 + t188 - t351 - t411 + t527 + t63) * t64 + t550 * (t303 + t117) + t551 * (t116 - t488)) * m(5) + (t112 * (t251 - t383) + t316 + t548 * (t140 + t303) + t549 * (t139 - t488) + (t235 + t384 + t112 - t352 + t524) * t113) * m(4); t313 + (t481 * t547 + t314 + t554 * t79 + t555 * t78 + (t188 - t338 + t522) * t48 + (t320 + t518) * t47) * m(6) + (t315 + (t306 * t486 - t193 - t423 + t527) * t64 + (-t277 + t348 + t519) * t63 + t550 * t117 + t551 * t116) * m(5) + (t316 + t548 * t140 + t549 * t139 + (-t423 + t433 + t524) * t113 + (t251 + t536) * t112) * m(4) + (-t157 * t161 - t158 * t462 + (t157 * t306 + t545) * t204 + (t158 * t306 - t546) * t201) * m(3); (-m(4) - m(5) - m(6)) * t381 + 0.2e1 * (t494 * t7 + t495 * t6) * m(6) + 0.2e1 * (t24 * t495 + t25 * t494) * m(5) + 0.2e1 * (t494 * t51 + t495 * t50) * m(4); t580 * t497 + t579 * t496 + (t566 * t306 - t560 * t305 + t540 * t184 + t539 * t183 + (t552 * t298 + t553 * t299) * qJD(4)) * t495 + (t565 * t306 + t561 * t305 + t542 * t184 + t541 * t183 + (t543 * t298 + t544 * t299) * qJD(4)) * t299 / 0.2e1 + (-t298 * t586 + t538 * t299) * t305 / 0.2e1 - (((-t429 + t431) * t310 + (-t430 + t432) * t308) * t306 + (((-t445 + t447) * t299 + (t444 - t446) * t298) * t310 + ((-t441 - t443) * t299 + (t440 + t442) * t298) * t308) * qJD(4)) * t306 / 0.2e1 + ((-t306 * t586 + t568) * t299 + (-t306 * t538 + t567) * t298) * t306 / 0.2e1 - t570 * t457 / 0.2e1 + t569 * t454 / 0.2e1 + ((-t421 * t458 - t533) * t298 + (-t563 + (t501 * t299 + (t166 - t502) * t298) * qJD(4)) * t299 + (-t421 * t460 - t461) * t298 + (t562 + (t503 * t299 + (t164 - t504) * t298) * qJD(4)) * t299) * t395 + ((t306 * t539 + t553) * t299 + (-t306 * t540 + t552) * t298) * t394 + ((t166 * t420 - t533) * t299 + (t563 + (t502 * t298 + (-t458 - t501) * t299) * qJD(4)) * t298 + (t164 * t420 - t461) * t299 + (-t562 + (t504 * t298 + (-t460 - t503) * t299) * qJD(4)) * t298) * t393 + ((t306 * t541 + t544) * t299 + (-t306 * t542 + t543) * t298) * t392 + (-g(1) * t437 - g(3) * t427 + t426 * t485 - (t349 + t480) * qJD(5) - (t436 * t47 + t437 * t48) * t306 - ((-t427 * t48 - t436 * t52) * t299 + (t427 * t47 - t437 * t52) * t298) * qJD(4) + (-t7 * t426 - t48 * t435 - t5 * t438 + t52 * t478 + (t426 * t47 - t439 * t52) * t306) * t299 + (t6 * t426 + t47 * t435 - t5 * t439 + t52 * t450 + (t426 * t48 + t438 * t52) * t306) * t298) * m(6) + (-(t177 * t64 + t181 * t63) * t306 - (t77 * (-t177 * t298 - t181 * t299) - t372 * t378) * qJD(4) + (-t154 * t183 - t156 * t184 + (-t101 * t298 + t299 * t99) * qJD(4)) * t358 + t77 * ((-t154 * t306 + t99) * t299 + (-t101 + t137) * t298) - t372 * t218 + ((t306 * t63 - t25) * t299 + (t306 * t64 + t24) * t298) * t258 - g(1) * t177 + g(2) * t181 + g(3) * t378) * m(5); (-(-t299 * t47 - t481) * t451 - ((t299 ^ 2 + t500) * t480 + t349) * qJD(4) + (qJD(4) * t52 + (-t306 * t47 + t7) * t299 + (-t306 * t48 - t6) * t298 + t381) * t310 + (qJD(4) * t377 - g(3) + t5) * t308) * m(6);];
tau = t1;
