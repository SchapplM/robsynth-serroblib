% Calculate vector of inverse dynamics joint torques for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP6_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:42
% DurationCPUTime: 27.43s
% Computational Cost: add. (7687->625), mult. (21171->874), div. (0->0), fcn. (20668->6), ass. (0->306)
t557 = Icges(5,4) + Icges(4,5);
t556 = Icges(4,6) - Icges(5,6);
t539 = Icges(4,2) + Icges(5,3);
t578 = Icges(5,2) + Icges(4,3);
t276 = sin(qJ(3));
t277 = sin(qJ(2));
t278 = cos(qJ(3));
t577 = (-t557 * t276 - t556 * t278) * t277;
t274 = sin(pkin(6));
t275 = cos(pkin(6));
t279 = cos(qJ(2));
t402 = qJD(3) * t279;
t384 = t278 * t402;
t404 = qJD(2) * t277;
t132 = -t274 * t384 + (qJD(3) * t275 + t274 * t404) * t276;
t428 = t276 * t279;
t212 = t274 * t428 + t275 * t278;
t386 = t278 * t404;
t133 = -qJD(3) * t212 - t274 * t386;
t403 = qJD(2) * t279;
t388 = t274 * t403;
t576 = t556 * t132 + t557 * t133 + t578 * t388;
t429 = t276 * t277;
t393 = t275 * t429;
t434 = t274 * t276;
t134 = qJD(2) * t393 - qJD(3) * t434 - t275 * t384;
t214 = -t274 * t278 + t275 * t428;
t135 = -qJD(3) * t214 - t275 * t386;
t387 = t275 * t403;
t575 = t556 * t134 + t557 * t135 + t578 * t387;
t426 = t278 * t279;
t213 = t274 * t426 - t275 * t276;
t433 = t274 * t277;
t89 = Icges(4,5) * t213 - Icges(4,6) * t212 + Icges(4,3) * t433;
t91 = Icges(5,4) * t213 + Icges(5,2) * t433 + Icges(5,6) * t212;
t547 = t89 + t91;
t215 = t275 * t426 + t434;
t431 = t275 * t277;
t90 = Icges(4,5) * t215 - Icges(4,6) * t214 + Icges(4,3) * t431;
t92 = Icges(5,4) * t215 + Icges(5,2) * t431 + Icges(5,6) * t214;
t546 = t90 + t92;
t574 = Icges(4,1) + Icges(5,1);
t573 = Icges(4,4) - Icges(5,5);
t343 = Icges(4,5) * t278 - Icges(4,6) * t276;
t177 = Icges(4,3) * t277 + t279 * t343;
t346 = Icges(5,4) * t278 + Icges(5,6) * t276;
t179 = Icges(5,2) * t277 + t279 * t346;
t572 = -t577 * qJD(3) + (-t177 - t179) * qJD(2);
t441 = Icges(5,5) * t278;
t342 = Icges(5,3) * t276 + t441;
t175 = Icges(5,6) * t277 + t279 * t342;
t447 = Icges(4,4) * t278;
t347 = -Icges(4,2) * t276 + t447;
t181 = Icges(4,6) * t277 + t279 * t347;
t515 = -t175 + t181;
t176 = -Icges(4,3) * t279 + t277 * t343;
t178 = -Icges(5,2) * t279 + t277 * t346;
t561 = -t178 - t176;
t442 = Icges(5,5) * t276;
t350 = Icges(5,1) * t278 + t442;
t183 = Icges(5,4) * t277 + t279 * t350;
t448 = Icges(4,4) * t276;
t351 = Icges(4,1) * t278 - t448;
t185 = Icges(4,5) * t277 + t279 * t351;
t514 = -t183 - t185;
t571 = (-t539 * t278 + t442 - t448) * t277;
t570 = -t539 * t132 - t573 * t133 - t556 * t388;
t569 = -t539 * t134 - t573 * t135 - t556 * t387;
t568 = t573 * t132 + t574 * t133 + t557 * t388;
t567 = t573 * t134 + t574 * t135 + t557 * t387;
t202 = Icges(5,5) * t213;
t87 = Icges(5,6) * t433 + Icges(5,3) * t212 + t202;
t450 = Icges(4,4) * t213;
t93 = -Icges(4,2) * t212 + Icges(4,6) * t433 + t450;
t549 = t87 - t93;
t203 = Icges(5,5) * t215;
t88 = Icges(5,6) * t431 + Icges(5,3) * t214 + t203;
t449 = Icges(4,4) * t215;
t94 = -Icges(4,2) * t214 + Icges(4,6) * t431 + t449;
t548 = t88 - t94;
t444 = Icges(5,5) * t212;
t95 = Icges(5,1) * t213 + Icges(5,4) * t433 + t444;
t204 = Icges(4,4) * t212;
t97 = Icges(4,1) * t213 + Icges(4,5) * t433 - t204;
t566 = t95 + t97;
t443 = Icges(5,5) * t214;
t96 = Icges(5,1) * t215 + Icges(5,4) * t431 + t443;
t205 = Icges(4,4) * t214;
t98 = Icges(4,1) * t215 + Icges(4,5) * t431 - t205;
t565 = t96 + t98;
t564 = t515 * qJD(2) + t571 * qJD(3);
t231 = (-Icges(4,1) * t276 - t447) * t277;
t398 = t277 * qJD(3);
t563 = -(-Icges(5,1) * t276 + t441) * t398 - qJD(3) * t231 + t514 * qJD(2);
t427 = t277 * t278;
t260 = Icges(5,5) * t427;
t438 = Icges(5,6) * t279;
t174 = Icges(5,3) * t429 + t260 - t438;
t180 = -Icges(4,6) * t279 + t277 * t347;
t562 = t174 - t180;
t182 = -Icges(5,4) * t279 + t277 * t350;
t184 = -Icges(4,5) * t279 + t277 * t351;
t537 = t182 + t184;
t560 = t572 * t277 + t561 * t403;
t559 = t575 * t277 + t546 * t403;
t558 = t576 * t277 + t547 * t403;
t529 = -t549 * t132 + t566 * t133 + t570 * t212 + t568 * t213 + t558 * t274;
t528 = -t548 * t132 + t565 * t133 + t569 * t212 + t567 * t213 + t559 * t274;
t527 = -t549 * t134 + t566 * t135 + t570 * t214 + t568 * t215 + t558 * t275;
t526 = -t548 * t134 + t565 * t135 + t569 * t214 + t567 * t215 + t559 * t275;
t555 = t562 * t132 - t537 * t133 + t564 * t212 + t563 * t213 + t560 * t274;
t554 = t562 * t134 - t537 * t135 + t564 * t214 + t563 * t215 + t560 * t275;
t543 = t549 * t212 + t566 * t213 + t547 * t433;
t553 = t548 * t212 + t565 * t213 + t546 * t433;
t552 = t549 * t214 + t566 * t215 + t547 * t431;
t523 = t548 * t214 + t565 * t215 + t546 * t431;
t357 = t276 * t87 + t278 * t95;
t47 = t277 * t357 - t279 * t91;
t355 = -t276 * t93 + t278 * t97;
t49 = t277 * t355 - t279 * t89;
t522 = t47 + t49;
t356 = t276 * t88 + t278 * t96;
t48 = t277 * t356 - t279 * t92;
t354 = -t276 * t94 + t278 * t98;
t50 = t277 * t354 - t279 * t90;
t521 = t48 + t50;
t551 = t562 * t212 + t537 * t213 - t561 * t433;
t550 = t562 * t214 + t537 * t215 - t561 * t431;
t337 = t174 * t276 + t182 * t278;
t435 = t178 * t279;
t68 = t277 * t337 - t435;
t336 = -t180 * t276 + t184 * t278;
t436 = t176 * t279;
t69 = t277 * t336 - t436;
t496 = t68 + t69;
t520 = rSges(5,3) + qJ(4);
t530 = rSges(5,1) + pkin(3);
t410 = t520 * t427 - t429 * t530;
t272 = t274 ^ 2;
t273 = t275 ^ 2;
t246 = rSges(3,1) * t277 + rSges(3,2) * t279;
t513 = t272 + t273;
t542 = t246 * t513;
t541 = t274 * t275;
t406 = qJD(2) * t274;
t237 = t275 * t398 + t406;
t405 = qJD(2) * t275;
t238 = t274 * t398 - t405;
t536 = t577 * t402 + (t212 * t557 + t556 * t213) * t238 + (t214 * t557 + t556 * t215) * t237;
t535 = (t237 * t553 + t238 * t543 - t402 * t551) * t274 + (t237 * t523 + t238 * t552 - t402 * t550) * t275;
t534 = 0.2e1 * qJD(2);
t533 = 2 * qJDD(2);
t397 = qJD(2) * qJD(3);
t310 = qJDD(3) * t277 + t279 * t397;
t396 = qJDD(2) * t274;
t152 = t275 * t310 + t396;
t395 = qJDD(2) * t275;
t153 = t274 * t310 - t395;
t241 = -qJDD(3) * t279 + t277 * t397;
t532 = t553 * t152 + t543 * t153 + t528 * t237 + t238 * t529 + t551 * t241 + t555 * t402;
t531 = t152 * t523 + t153 * t552 + t237 * t526 + t238 * t527 + t241 * t550 + t402 * t554;
t499 = ((t355 + t357) * qJD(2) - t576) * t279 + (t568 * t278 + t570 * t276 + (-t276 * t566 + t278 * t549) * qJD(3) + t547 * qJD(2)) * t277;
t498 = ((t354 + t356) * qJD(2) - t575) * t279 + (t567 * t278 + t569 * t276 + (-t276 * t565 + t278 * t548) * qJD(3) + t546 * qJD(2)) * t277;
t525 = t237 * t521 + t238 * t522 - t402 * t496;
t293 = -t277 * t342 + t438;
t136 = t293 * t274;
t142 = t180 * t274;
t519 = t136 + t142;
t137 = t293 * t275;
t143 = t180 * t275;
t518 = t137 + t143;
t144 = t182 * t274;
t146 = t184 * t274;
t517 = -t144 - t146;
t145 = t182 * t275;
t147 = t184 * t275;
t516 = -t145 - t147;
t323 = t177 - t336;
t324 = -t179 + t337;
t480 = -(-t176 * t275 - t354) * t237 - (-t176 * t274 - t355) * t238;
t481 = (t178 * t275 + t356) * t237 + (t178 * t274 + t357) * t238;
t512 = (-t480 - t481 + (-t323 + t324) * t402) * t277;
t367 = rSges(5,1) * t278 + rSges(5,3) * t276;
t511 = (pkin(3) * t278 + qJ(4) * t276 + t367) * t277;
t510 = t237 * t546 + t238 * t547;
t248 = pkin(2) * t277 - pkin(5) * t279;
t322 = qJD(2) * t248;
t200 = t274 * t322;
t201 = t275 * t322;
t432 = t274 * t279;
t256 = pkin(5) * t432;
t430 = t275 * t279;
t258 = pkin(5) * t430;
t509 = -t274 * t200 - t275 * t201 - (-pkin(2) * t433 + t256) * t406 - (-pkin(2) * t431 + t258) * t405;
t508 = -t435 - t436;
t507 = t550 * t277;
t506 = t551 * t277;
t505 = t553 * t275;
t504 = t552 * t274;
t491 = t279 * pkin(2) + t277 * pkin(5);
t240 = t491 * qJD(2);
t390 = t248 * t406;
t401 = qJD(4) * t212;
t414 = -rSges(5,2) * t279 + t511;
t425 = rSges(5,2) * t431 + t214 * t520 + t215 * t530;
t36 = -t237 * t414 - t402 * t425 - t390 + t401;
t389 = t248 * t405;
t400 = qJD(4) * t214;
t453 = rSges(5,2) * t433 + t212 * t520 + t213 * t530;
t37 = t238 * t414 + t402 * t453 - t389 + t400;
t500 = t237 * t36 - t238 * t37 - g(3);
t497 = ((-t336 - t337) * qJD(2) - t572) * t279 + (t563 * t278 + t564 * t276 + (t276 * t537 - t278 * t562) * qJD(3) + t561 * qJD(2)) * t277;
t489 = (Icges(5,1) * t429 - t231 - t260 - t562) * t402 + (-t212 * t574 + t202 - t450 + t549) * t238 + (-t214 * t574 + t203 - t449 + t548) * t237;
t488 = (t537 + t571) * t402 + (t213 * t539 + t204 - t444 - t566) * t238 + (t215 * t539 + t205 - t443 - t565) * t237;
t487 = t536 * t277;
t358 = t49 * t274 + t50 * t275;
t359 = t47 * t274 + t48 * t275;
t486 = t358 + t359;
t485 = t275 * t523 + t504;
t484 = t274 * t543 + t505;
t483 = t277 * (g(1) * t275 + g(2) * t274);
t399 = qJD(4) * t276;
t255 = t277 * t399;
t224 = t491 * t274;
t225 = t491 * t275;
t381 = t224 * t406 + t225 * t405 + qJD(1);
t33 = t237 * t453 - t238 * t425 + t255 + t381;
t463 = rSges(5,2) * t387 - t134 * t520 + t135 * t530 + t400;
t385 = t278 * t398;
t303 = t276 * t403 + t385;
t334 = -t200 * t406 - t201 * t405 + t224 * t396 + t225 * t395 + qJDD(1);
t464 = rSges(5,2) * t388 - t132 * t520 + t133 * t530 + t401;
t5 = qJD(4) * t303 + qJDD(4) * t429 + t152 * t453 - t153 * t425 + t237 * t464 - t238 * t463 + t334;
t482 = t33 * t463 + t425 * t5;
t477 = t152 / 0.2e1;
t476 = t153 / 0.2e1;
t475 = -t237 / 0.2e1;
t474 = t237 / 0.2e1;
t473 = -t238 / 0.2e1;
t472 = t238 / 0.2e1;
t471 = t241 / 0.2e1;
t462 = rSges(4,1) * t278;
t269 = t277 * rSges(5,2);
t268 = t277 * rSges(4,3);
t424 = t255 + t303 * qJ(4) + (-t276 * t398 + t278 * t403) * pkin(3) + (-rSges(5,1) * t276 + rSges(5,3) * t278) * t398 + (t279 * t367 + t269) * qJD(2);
t234 = (-rSges(4,1) * t276 - rSges(4,2) * t278) * t277;
t368 = -rSges(4,2) * t276 + t462;
t111 = qJD(3) * t234 + (t279 * t368 + t268) * qJD(2);
t423 = -t111 - t240;
t422 = -t212 * t530 + t213 * t520;
t421 = t214 * t530 - t215 * t520;
t252 = rSges(5,2) * t432;
t420 = -t274 * t511 + t252;
t254 = rSges(5,2) * t430;
t419 = t275 * t511 - t254;
t187 = -rSges(4,3) * t279 + t277 * t368;
t413 = -t187 - t248;
t411 = t274 * t224 + t275 * t225;
t409 = t274 * rSges(4,2) * t429 + rSges(4,3) * t432;
t408 = rSges(4,2) * t393 + rSges(4,3) * t430;
t394 = rSges(4,1) * t427;
t392 = -t240 - t424;
t391 = -t248 - t414;
t376 = -t402 / 0.2e1;
t375 = t402 / 0.2e1;
t247 = rSges(3,1) * t279 - rSges(3,2) * t277;
t345 = Icges(3,5) * t279 - Icges(3,6) * t277;
t344 = -Icges(3,5) * t277 - Icges(3,6) * t279;
t100 = rSges(4,1) * t213 - rSges(4,2) * t212 + rSges(4,3) * t433;
t102 = rSges(4,1) * t215 - rSges(4,2) * t214 + rSges(4,3) * t431;
t341 = t100 * t275 - t102 * t274;
t338 = t513 * t247;
t335 = qJD(2) * t542;
t189 = rSges(4,1) * t426 - rSges(4,2) * t428 + t268;
t333 = -qJD(2) * t240 - qJDD(2) * t248;
t46 = t100 * t237 - t102 * t238 + t381;
t315 = t46 * t341;
t307 = qJD(2) * t344;
t306 = t333 * t274;
t305 = t333 * t275;
t304 = t33 * t464 + t453 * t5;
t302 = t36 * t425 - t37 * t453;
t283 = (t33 * t453 - t36 * t414) * t275 + (-t33 * t425 + t37 * t414) * t274;
t223 = t246 * t275;
t222 = t246 * t274;
t217 = t344 * t275;
t216 = t344 * t274;
t193 = t275 * t307;
t192 = t274 * t307;
t159 = Icges(3,3) * t274 + t275 * t345;
t158 = -Icges(3,3) * t275 + t274 * t345;
t151 = -t275 * t394 + t408;
t149 = -t274 * t394 + t409;
t130 = -rSges(4,1) * t214 - rSges(4,2) * t215;
t126 = -rSges(4,1) * t212 - rSges(4,2) * t213;
t85 = rSges(4,1) * t135 + rSges(4,2) * t134 + rSges(4,3) * t387;
t83 = rSges(4,1) * t133 + rSges(4,2) * t132 + rSges(4,3) * t388;
t61 = t100 * t402 + t187 * t238 - t389;
t60 = -t102 * t402 - t187 * t237 - t390;
t55 = -qJD(2) * t335 + qJDD(2) * t338 + qJDD(1);
t35 = -t100 * t241 + t111 * t238 + t153 * t187 + t402 * t83 + t305;
t34 = t102 * t241 - t111 * t237 - t152 * t187 - t402 * t85 + t306;
t26 = t100 * t152 - t102 * t153 + t237 * t83 - t238 * t85 + t334;
t7 = -qJD(4) * t134 + qJDD(4) * t214 + t153 * t414 + t238 * t424 - t241 * t453 + t402 * t464 + t305;
t6 = -qJD(4) * t132 + qJDD(4) * t212 - t152 * t414 - t237 * t424 + t241 * t425 - t402 * t463 + t306;
t1 = [m(2) * qJDD(1) + m(3) * t55 + m(4) * t26 + m(5) * t5 + (-m(2) - m(3) - m(4) - m(5)) * g(3); (t216 * qJD(2) * t273 - t217 * t275 * t406) * t405 / 0.2e1 - (t217 * qJD(2) * t272 - t216 * t274 * t405) * t406 / 0.2e1 + (t523 * t274 - t275 * t552) * t477 + (t274 * t553 - t275 * t543) * t476 + (((t214 * t515 + t215 * t514 + t504) * t279 + t507) * qJD(3) + (((t508 + t523) * qJD(3) + t510) * t279 + t512) * t275 + (t214 * t519 + t215 * t517) * t238 + (t214 * t518 + t215 * t516) * t237) * t475 + (t274 * t526 - t275 * t527) * t474 + (((t212 * t515 + t213 * t514 + t505) * t279 + t506) * qJD(3) + (((t508 + t543) * qJD(3) + t510) * t279 + t512) * t274 + (t212 * t519 + t213 * t517) * t238 + (t212 * t518 + t213 * t516) * t237) * t473 + (t274 * t528 - t275 * t529) * t472 + (t274 * t521 - t275 * t522) * t471 + (((t143 * t276 - t147 * t278 + t90) * t237 + (t142 * t276 - t146 * t278 + t89) * t238 + t69 * qJD(3)) * t277 + ((t323 * t279 + (t181 * t276 - t185 * t278 - t176) * t277 + t358) * qJD(3) + t480) * t279 + ((t137 * t276 - t145 * t278 + t92) * t237 + (t136 * t276 - t144 * t278 + t91) * t238 + t68 * qJD(3)) * t277 + ((-t324 * t279 + (-t175 * t276 - t183 * t278 - t178) * t277 + t359) * qJD(3) + t481) * t279) * t375 - t525 * t398 / 0.2e1 + (-g(1) * (t254 + t258) - g(2) * (t252 + t256) - g(3) * t491 - (-t276 * t520 - t278 * t530 - pkin(2)) * t483 - (t302 * t277 + (t36 * t419 + t37 * t420 + t283) * t279) * qJD(3) - (t274 * t36 + t275 * t37) * (-t255 - t240) + t5 * t411 + (t37 * t392 + t391 * t7 + t482) * t275 + (t36 * t392 + t391 * t6 + t304) * t274 + t500 * (t426 * t530 + t428 * t520 + t269) + (-t237 * t420 - t238 * t419 - t279 * t399 + t509) * t33) * m(5) + (t26 * t411 + (t26 * t102 + t35 * t413 + t423 * t61) * t275 + (t26 * t100 + t34 * t413 + t423 * t60) * t274 - t61 * (t189 * t238 - t405 * t491) - t60 * (-t189 * t237 - t406 * t491) - ((-t100 * t61 + t102 * t60) * t277 + (t61 * (t187 * t274 + t149) + t60 * (-t187 * t275 - t151) + t315) * t279) * qJD(3) - g(1) * (t258 + t408) - g(2) * (t256 + t409) - g(3) * (t189 + t491) - (-pkin(2) - t462) * t483 + (-t149 * t237 + t151 * t238 + t274 * t83 + t275 * t85 + t509) * t46) * m(4) + (g(1) * t223 + g(2) * t222 - g(3) * t247 + t55 * t338 + (-t335 - (-t222 * t274 - t223 * t275) * qJD(2)) * (qJD(2) * t338 + qJD(1)) + (t247 * qJD(2) ^ 2 + qJDD(2) * t246) * t542) * m(3) + ((-t158 * t541 + t159 * t272) * t533 + (-t192 * t541 + t193 * t272) * t534 + t531) * t274 / 0.2e1 - ((t158 * t273 - t159 * t541) * t533 + (t192 * t273 - t193 * t541) * t534 + t532) * t275 / 0.2e1 + (t274 * t498 - t275 * t499 + t535) * t376; (t485 * t277 - t279 * t550) * t477 + (t484 * t277 - t279 * t551) * t476 + (t214 * t488 + t215 * t489 - t275 * t487) * t475 + (t554 * t279 + (t274 * t527 + t275 * t526) * t277 + (t279 * t485 + t507) * qJD(2)) * t474 + (t212 * t488 + t213 * t489 - t274 * t487) * t473 + (t555 * t279 + (t274 * t529 + t275 * t528) * t277 + (t279 * t484 + t506) * qJD(2)) * t472 + (t277 * t486 - t279 * t496) * t471 - (t152 * t521 + t153 * t522 + t498 * t237 + t499 * t238 + t496 * t241 + t497 * t402) * t279 / 0.2e1 + t532 * t433 / 0.2e1 + t531 * t431 / 0.2e1 + t525 * t404 / 0.2e1 + (t497 * t279 + (t274 * t499 + t275 * t498) * t277 + (t277 * t496 + t279 * t486) * qJD(2)) * t376 + (t536 * t279 + (t276 * t488 + t278 * t489) * t277) * t375 + (-(t213 * t36 + t215 * t37 + t33 * t427) * qJD(4) - (t33 * t421 + t37 * t410) * t238 - (t33 * t422 - t36 * t410) * t237 - (t36 * t421 + t37 * t422) * t402 + (qJD(2) * t283 - t36 * t463 + t37 * t464 - t425 * t6 + t453 * t7) * t279 + (t302 * qJD(2) + (-t36 * t424 - t414 * t6 + t304) * t275 + (t37 * t424 + t414 * t7 - t482) * t274) * t277 + g(1) * t421 - g(2) * t422 - g(3) * t410) * m(5) + (-t61 * (t126 * t402 + t234 * t238) - t60 * (-t130 * t402 - t234 * t237) - t46 * (t126 * t237 - t130 * t238) + (t35 * t100 - t34 * t102 - t60 * t85 + t61 * t83 + (t315 + (t274 * t61 - t275 * t60) * t187) * qJD(2)) * t279 + (t61 * (-qJD(2) * t100 + t111 * t274) + t60 * (qJD(2) * t102 - t111 * t275) + t26 * t341 + t46 * (-t274 * t85 + t275 * t83) + (t274 * t35 - t275 * t34) * t187) * t277 - g(1) * t130 - g(2) * t126 - g(3) * t234) * m(4) + t535 * t403 / 0.2e1; (t33 * t385 - t132 * t36 - t134 * t37 + (t277 * t5 + t33 * t403) * t276 + t500 * t429 + (t238 * t33 + t36 * t402 - g(1) + t7) * t214 + (-t237 * t33 - t37 * t402 - g(2) + t6) * t212) * m(5);];
tau = t1;
