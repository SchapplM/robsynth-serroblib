% Calculate vector of inverse dynamics joint torques for
% S4PRRP5
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:29:21
% DurationCPUTime: 28.10s
% Computational Cost: add. (7927->631), mult. (21105->877), div. (0->0), fcn. (20353->6), ass. (0->322)
t580 = Icges(4,5) + Icges(5,5);
t578 = Icges(4,6) + Icges(5,6);
t561 = Icges(4,1) + Icges(5,1);
t599 = Icges(4,2) + Icges(5,2);
t604 = Icges(4,3) + Icges(5,3);
t270 = sin(qJ(3));
t271 = sin(qJ(2));
t272 = cos(qJ(3));
t603 = (-t580 * t270 - t578 * t272) * t271;
t267 = sin(pkin(6));
t273 = cos(qJ(2));
t431 = t272 * t273;
t268 = cos(pkin(6));
t437 = t268 * t270;
t213 = t267 * t431 - t437;
t405 = qJD(2) * t271;
t382 = t270 * t405;
t132 = -qJD(3) * t213 + t267 * t382;
t433 = t270 * t273;
t212 = -t267 * t433 - t268 * t272;
t381 = t272 * t405;
t133 = qJD(3) * t212 - t267 * t381;
t404 = qJD(2) * t273;
t384 = t267 * t404;
t602 = t578 * t132 + t580 * t133 + t604 * t384;
t441 = t267 * t270;
t215 = t268 * t431 + t441;
t134 = -qJD(3) * t215 + t268 * t382;
t529 = -t267 * t272 + t268 * t433;
t135 = -qJD(3) * t529 - t268 * t381;
t383 = t268 * t404;
t601 = t578 * t134 + t580 * t135 + t604 * t383;
t440 = t267 * t271;
t87 = Icges(5,5) * t213 + Icges(5,6) * t212 + Icges(5,3) * t440;
t89 = Icges(4,5) * t213 + Icges(4,6) * t212 + Icges(4,3) * t440;
t571 = t89 + t87;
t436 = t268 * t271;
t88 = Icges(5,5) * t215 - Icges(5,6) * t529 + Icges(5,3) * t436;
t90 = Icges(4,5) * t215 - Icges(4,6) * t529 + Icges(4,3) * t436;
t570 = t90 + t88;
t600 = Icges(4,4) + Icges(5,4);
t335 = Icges(5,5) * t272 - Icges(5,6) * t270;
t175 = Icges(5,3) * t271 + t273 * t335;
t336 = Icges(4,5) * t272 - Icges(4,6) * t270;
t177 = Icges(4,3) * t271 + t273 * t336;
t598 = -t603 * qJD(3) + (-t175 - t177) * qJD(2);
t174 = -Icges(5,3) * t273 + t271 * t335;
t176 = -Icges(4,3) * t273 + t271 * t336;
t584 = t174 + t176;
t450 = Icges(5,4) * t272;
t339 = -Icges(5,2) * t270 + t450;
t179 = Icges(5,6) * t271 + t273 * t339;
t454 = Icges(4,4) * t272;
t340 = -Icges(4,2) * t270 + t454;
t181 = Icges(4,6) * t271 + t273 * t340;
t597 = t179 + t181;
t451 = Icges(5,4) * t270;
t343 = Icges(5,1) * t272 - t451;
t183 = Icges(5,5) * t271 + t273 * t343;
t455 = Icges(4,4) * t270;
t344 = Icges(4,1) * t272 - t455;
t185 = Icges(4,5) * t271 + t273 * t344;
t535 = -t183 - t185;
t596 = (-t599 * t272 - t451 - t455) * t271;
t595 = (t561 * t270 + t450 + t454) * t271;
t594 = t599 * t132 + t600 * t133 + t578 * t384;
t593 = t599 * t134 + t600 * t135 + t578 * t383;
t592 = t600 * t132 + t561 * t133 + t580 * t384;
t591 = t600 * t134 + t561 * t135 + t580 * t383;
t453 = Icges(5,4) * t213;
t91 = Icges(5,2) * t212 + Icges(5,6) * t440 + t453;
t457 = Icges(4,4) * t213;
t93 = Icges(4,2) * t212 + Icges(4,6) * t440 + t457;
t590 = t91 + t93;
t452 = Icges(5,4) * t215;
t92 = -Icges(5,2) * t529 + Icges(5,6) * t436 + t452;
t456 = Icges(4,4) * t215;
t94 = -Icges(4,2) * t529 + Icges(4,6) * t436 + t456;
t589 = t92 + t94;
t202 = Icges(5,4) * t212;
t95 = Icges(5,1) * t213 + Icges(5,5) * t440 + t202;
t204 = Icges(4,4) * t212;
t97 = Icges(4,1) * t213 + Icges(4,5) * t440 + t204;
t588 = t95 + t97;
t203 = Icges(5,4) * t529;
t96 = Icges(5,1) * t215 + Icges(5,5) * t436 - t203;
t205 = Icges(4,4) * t529;
t98 = Icges(4,1) * t215 + Icges(4,5) * t436 - t205;
t587 = t96 + t98;
t586 = qJD(2) * t597 + qJD(3) * t596;
t585 = qJD(2) * t535 + qJD(3) * t595;
t178 = -Icges(5,6) * t273 + t271 * t339;
t180 = -Icges(4,6) * t273 + t271 * t340;
t560 = t178 + t180;
t182 = -Icges(5,5) * t273 + t271 * t343;
t184 = -Icges(4,5) * t273 + t271 * t344;
t565 = t182 + t184;
t583 = t598 * t271 - t584 * t404;
t582 = t601 * t271 + t570 * t404;
t581 = t602 * t271 + t571 * t404;
t550 = rSges(5,1) + pkin(3);
t549 = t590 * t132 + t588 * t133 + t594 * t212 + t592 * t213 + t581 * t267;
t548 = t589 * t132 + t587 * t133 + t593 * t212 + t591 * t213 + t582 * t267;
t547 = t590 * t134 + t588 * t135 + t592 * t215 + t581 * t268 - t594 * t529;
t546 = t589 * t134 + t587 * t135 + t591 * t215 + t582 * t268 - t593 * t529;
t577 = -t560 * t132 - t565 * t133 - t586 * t212 + t585 * t213 + t583 * t267;
t576 = -t560 * t134 - t565 * t135 + t585 * t215 + t583 * t268 + t586 * t529;
t564 = t590 * t212 + t588 * t213 + t571 * t440;
t575 = t589 * t212 + t587 * t213 + t570 * t440;
t574 = t588 * t215 + t571 * t436 - t590 * t529;
t543 = t587 * t215 + t570 * t436 - t589 * t529;
t350 = -t270 * t91 + t272 * t95;
t47 = t271 * t350 - t273 * t87;
t348 = -t270 * t93 + t272 * t97;
t49 = t271 * t348 - t273 * t89;
t542 = t47 + t49;
t349 = -t270 * t92 + t272 * t96;
t48 = t271 * t349 - t273 * t88;
t347 = -t270 * t94 + t272 * t98;
t50 = t271 * t347 - t273 * t90;
t541 = t48 + t50;
t573 = t560 * t212 + t565 * t213 + t584 * t440;
t572 = t565 * t215 + t584 * t436 - t560 * t529;
t330 = -t178 * t270 + t182 * t272;
t443 = t174 * t273;
t66 = t271 * t330 - t443;
t329 = -t180 * t270 + t184 * t272;
t442 = t176 * t273;
t67 = t271 * t329 - t442;
t516 = t66 + t67;
t403 = qJD(3) * t271;
t379 = t268 * t403;
t407 = qJD(2) * t267;
t234 = t379 + t407;
t380 = t267 * t403;
t406 = qJD(2) * t268;
t235 = t380 - t406;
t402 = qJD(3) * t273;
t507 = t234 * (-t215 * t599 - t203 - t205 + t587) + t235 * (-t213 * t599 + t202 + t204 + t588) - t402 * (t565 + t596);
t265 = t267 ^ 2;
t266 = t268 ^ 2;
t243 = rSges(3,1) * t271 + rSges(3,2) * t273;
t533 = t265 + t266;
t563 = t243 * t533;
t562 = t267 * t268;
t559 = t603 * t402 + (-t580 * t212 + t578 * t213) * t235 + (t578 * t215 + t580 * t529) * t234;
t555 = (t234 * t575 + t235 * t564 - t402 * t573) * t267 + (t234 * t543 + t235 * t574 - t402 * t572) * t268;
t554 = 0.2e1 * qJD(2);
t553 = 2 * qJDD(2);
t399 = qJD(2) * qJD(3);
t304 = qJDD(3) * t271 + t273 * t399;
t397 = qJDD(2) * t267;
t154 = t268 * t304 + t397;
t396 = qJDD(2) * t268;
t155 = t267 * t304 - t396;
t238 = -qJDD(3) * t273 + t271 * t399;
t552 = t575 * t154 + t564 * t155 + t548 * t234 + t235 * t549 + t573 * t238 + t577 * t402;
t551 = t154 * t543 + t155 * t574 + t234 * t546 + t235 * t547 + t238 * t572 + t402 * t576;
t519 = ((t348 + t350) * qJD(2) - t602) * t273 + (t592 * t272 - t594 * t270 + (-t270 * t588 - t272 * t590) * qJD(3) + t571 * qJD(2)) * t271;
t518 = ((t347 + t349) * qJD(2) - t601) * t273 + (t591 * t272 - t593 * t270 + (-t270 * t587 - t272 * t589) * qJD(3) + t570 * qJD(2)) * t271;
t545 = t234 * t541 + t235 * t542 - t402 * t516;
t140 = t178 * t267;
t142 = t180 * t267;
t540 = -t140 - t142;
t141 = t178 * t268;
t143 = t180 * t268;
t539 = -t141 - t143;
t144 = t182 * t267;
t146 = t184 * t267;
t538 = -t144 - t146;
t145 = t182 * t268;
t147 = t184 * t268;
t537 = -t145 - t147;
t510 = t273 * pkin(2) + t271 * pkin(5);
t237 = t510 * qJD(2);
t400 = qJD(4) * t273;
t534 = -t237 + t400;
t315 = t177 - t329;
t316 = t175 - t330;
t498 = -(-t174 * t268 - t349) * t234 - (-t174 * t267 - t350) * t235;
t499 = -(-t176 * t268 - t347) * t234 - (-t176 * t267 - t348) * t235;
t532 = (-t498 - t499 + (-t315 - t316) * t402) * t271;
t531 = t234 * t570 + t235 * t571;
t260 = pkin(3) * t272 + pkin(2);
t483 = pkin(2) - t260;
t373 = t483 * t271;
t269 = -qJ(4) - pkin(5);
t482 = pkin(5) + t269;
t530 = -t273 * t482 + t373;
t245 = pkin(2) * t271 - pkin(5) * t273;
t314 = qJD(2) * t245;
t200 = t267 * t314;
t201 = t268 * t314;
t438 = t267 * t273;
t255 = pkin(5) * t438;
t435 = t268 * t273;
t256 = pkin(5) * t435;
t528 = -t267 * t200 - t268 * t201 - (-pkin(2) * t440 + t255) * t407 - (-pkin(2) * t436 + t256) * t406;
t527 = -t442 - t443;
t526 = t572 * t271;
t525 = t573 * t271;
t524 = t575 * t268;
t523 = t574 * t267;
t517 = ((-t329 - t330) * qJD(2) - t598) * t273 + (t585 * t272 + t586 * t270 + (t270 * t565 + t560 * t272) * qJD(3) - t584 * qJD(2)) * t271;
t261 = t271 * rSges(5,3);
t511 = rSges(5,1) * t431 - rSges(5,2) * t433 + t273 * t260 - t269 * t271 + t261;
t508 = (t560 + t595) * t402 + (t212 * t561 - t453 - t457 - t590) * t235 + (-t529 * t561 - t452 - t456 - t589) * t234;
t506 = t559 * t271;
t351 = t49 * t267 + t50 * t268;
t352 = t47 * t267 + t48 * t268;
t505 = t351 + t352;
t504 = t268 * t543 + t523;
t503 = t267 * t564 + t524;
t401 = qJD(4) * t271;
t248 = t267 * t401;
t386 = t245 * t407;
t478 = rSges(5,1) * t272;
t359 = -rSges(5,2) * t270 + t478;
t422 = -rSges(5,3) * t273 + t271 * t359 - t530;
t299 = -t271 * t482 - t273 * t483;
t429 = rSges(5,1) * t215 - rSges(5,2) * t529 + rSges(5,3) * t436 + pkin(3) * t441 + t268 * t299;
t36 = -t234 * t422 - t402 * t429 + t248 - t386;
t249 = t268 * t401;
t385 = t245 * t406;
t464 = rSges(5,1) * t213 + rSges(5,2) * t212 + rSges(5,3) * t440 - pkin(3) * t437 + t267 * t299;
t37 = t235 * t422 + t402 * t464 + t249 - t385;
t502 = t267 * t36 + t268 * t37;
t501 = g(1) * t268 + g(2) * t267;
t224 = t510 * t267;
t225 = t510 * t268;
t376 = t224 * t407 + t225 * t406 + qJD(1);
t31 = t234 * t464 - t235 * t429 + t376 - t400;
t476 = pkin(3) * qJD(3);
t395 = t270 * t476;
t278 = qJD(2) * t530 - t273 * t395;
t394 = t272 * t476;
t480 = rSges(5,1) * t135 + rSges(5,2) * t134 + rSges(5,3) * t383 + t267 * t394 + t268 * t278 + t249;
t327 = -t200 * t407 - t201 * t406 + t224 * t397 + t225 * t396 + qJDD(1);
t481 = rSges(5,1) * t133 + rSges(5,2) * t132 + rSges(5,3) * t384 + t267 * t278 - t268 * t394 + t248;
t5 = qJD(2) * t401 - qJDD(4) * t273 + t154 * t464 - t155 * t429 + t234 * t481 - t235 * t480 + t327;
t500 = t31 * t480 + t429 * t5;
t495 = t154 / 0.2e1;
t494 = t155 / 0.2e1;
t493 = -t234 / 0.2e1;
t492 = t234 / 0.2e1;
t491 = -t235 / 0.2e1;
t490 = t235 / 0.2e1;
t489 = t238 / 0.2e1;
t479 = rSges(4,1) * t272;
t477 = rSges(5,2) * t272;
t262 = t271 * rSges(4,3);
t434 = t270 * t271;
t432 = t271 * t272;
t430 = t273 * t269;
t232 = (-rSges(5,1) * t270 - t477) * t271;
t428 = qJD(3) * t232 - t271 * t395 - t400 + (t273 * t359 + t261 + t299) * qJD(2);
t233 = (-rSges(4,1) * t270 - rSges(4,2) * t272) * t271;
t360 = -rSges(4,2) * t270 + t479;
t113 = qJD(3) * t233 + (t273 * t360 + t262) * qJD(2);
t427 = -t113 - t237;
t311 = t373 - t430;
t392 = t267 * t432;
t393 = t267 * t434;
t412 = rSges(5,2) * t393 + rSges(5,3) * t438;
t426 = -rSges(5,1) * t392 + t267 * t311 - t255 + t412;
t389 = t268 * t432;
t391 = t268 * t434;
t410 = rSges(5,2) * t391 + rSges(5,3) * t435;
t425 = rSges(5,1) * t389 - t268 * t311 + t256 - t410;
t424 = -t213 * rSges(5,2) + t212 * t550;
t423 = t215 * rSges(5,2) + t550 * t529;
t187 = -rSges(4,3) * t273 + t271 * t360;
t414 = -t187 - t245;
t413 = t267 * t224 + t268 * t225;
t411 = rSges(4,2) * t393 + rSges(4,3) * t438;
t409 = rSges(4,2) * t391 + rSges(4,3) * t435;
t398 = qJDD(2) * t245;
t388 = -t237 - t428;
t387 = -t245 - t422;
t369 = -t402 / 0.2e1;
t368 = t402 / 0.2e1;
t362 = pkin(3) * t434 - t232;
t244 = rSges(3,1) * t273 - rSges(3,2) * t271;
t338 = Icges(3,5) * t273 - Icges(3,6) * t271;
t337 = -Icges(3,5) * t271 - Icges(3,6) * t273;
t102 = rSges(4,1) * t213 + rSges(4,2) * t212 + rSges(4,3) * t440;
t104 = rSges(4,1) * t215 - rSges(4,2) * t529 + rSges(4,3) * t436;
t334 = t102 * t268 - t104 * t267;
t331 = t533 * t244;
t328 = qJD(2) * t563;
t189 = rSges(4,1) * t431 - rSges(4,2) * t433 + t262;
t325 = -qJD(2) * t237 - t398;
t46 = t102 * t234 - t104 * t235 + t376;
t309 = t46 * t334;
t301 = qJD(2) * t337;
t300 = t31 * t481 + t464 * t5;
t298 = t36 * t429 - t37 * t464;
t285 = qJD(2) * t534 + qJDD(4) * t271 - t398;
t277 = (t31 * t464 - t36 * t422) * t268 + (-t31 * t429 + t37 * t422) * t267;
t223 = t243 * t268;
t222 = t243 * t267;
t217 = t337 * t268;
t216 = t337 * t267;
t193 = t268 * t301;
t192 = t267 * t301;
t161 = Icges(3,3) * t267 + t268 * t338;
t160 = -Icges(3,3) * t268 + t267 * t338;
t151 = -rSges(4,1) * t389 + t409;
t149 = -rSges(4,1) * t392 + t411;
t131 = -rSges(4,1) * t529 - rSges(4,2) * t215;
t129 = rSges(4,1) * t212 - rSges(4,2) * t213;
t83 = rSges(4,1) * t135 + rSges(4,2) * t134 + rSges(4,3) * t383;
t81 = rSges(4,1) * t133 + rSges(4,2) * t132 + rSges(4,3) * t384;
t61 = t102 * t402 + t187 * t235 - t385;
t60 = -t104 * t402 - t187 * t234 - t386;
t55 = -qJD(2) * t328 + qJDD(2) * t331 + qJDD(1);
t35 = -t102 * t238 + t113 * t235 + t155 * t187 + t268 * t325 + t402 * t81;
t34 = t104 * t238 - t113 * t234 - t154 * t187 + t267 * t325 - t402 * t83;
t26 = t102 * t154 - t104 * t155 + t234 * t81 - t235 * t83 + t327;
t7 = t155 * t422 + t235 * t428 - t238 * t464 + t268 * t285 + t402 * t481;
t6 = -t154 * t422 - t234 * t428 + t238 * t429 + t267 * t285 - t402 * t480;
t1 = [m(2) * qJDD(1) + m(3) * t55 + m(4) * t26 + m(5) * t5 + (-m(2) - m(3) - m(4) - m(5)) * g(3); -(t217 * qJD(2) * t265 - t216 * t267 * t406) * t407 / 0.2e1 + (t216 * qJD(2) * t266 - t217 * t268 * t407) * t406 / 0.2e1 + (t543 * t267 - t268 * t574) * t495 + (t267 * t575 - t268 * t564) * t494 + (((t215 * t535 + t529 * t597 + t523) * t273 + t526) * qJD(3) + (((t527 + t543) * qJD(3) + t531) * t273 + t532) * t268 + (t215 * t538 - t529 * t540) * t235 + (t215 * t537 - t529 * t539) * t234) * t493 + (t267 * t546 - t268 * t547) * t492 + (((-t212 * t597 + t213 * t535 + t524) * t273 + t525) * qJD(3) + (((t527 + t564) * qJD(3) + t531) * t273 + t532) * t267 + (t212 * t540 + t213 * t538) * t235 + (t212 * t539 + t213 * t537) * t234) * t491 + (t267 * t548 - t268 * t549) * t490 + (t267 * t541 - t268 * t542) * t489 - t545 * t403 / 0.2e1 + (((t141 * t270 - t145 * t272 + t88) * t234 + (t140 * t270 - t144 * t272 + t87) * t235 + t66 * qJD(3)) * t271 + ((t316 * t273 + (t179 * t270 - t183 * t272 - t174) * t271 + t352) * qJD(3) + t498) * t273 + ((t143 * t270 - t147 * t272 + t90) * t234 + (t142 * t270 - t146 * t272 + t89) * t235 + t67 * qJD(3)) * t271 + ((t315 * t273 + (t181 * t270 - t185 * t272 - t176) * t271 + t351) * qJD(3) + t499) * t273) * t368 + (-g(1) * t410 - g(2) * t412 - g(3) * t511 - t501 * (-t430 + (-t260 - t478) * t271) - (t298 * t271 + (t36 * t425 + t37 * t426 + t277) * t273) * qJD(3) - t502 * t534 + t5 * t413 + (t37 * t388 + t387 * t7 + t500) * t268 + (t36 * t388 + t387 * t6 + t300) * t267 + (t36 * t234 - t37 * t235) * (-t510 + t511) + (-t426 * t234 - t425 * t235 - t401 + t528) * t31) * m(5) + (-g(1) * (t256 + t409) - g(2) * (t255 + t411) - g(3) * (t189 + t510) - t501 * t271 * (-pkin(2) - t479) - t61 * (t189 * t235 - t406 * t510) - t60 * (-t189 * t234 - t407 * t510) - ((-t102 * t61 + t104 * t60) * t271 + (t61 * (t187 * t267 + t149) + t60 * (-t187 * t268 - t151) + t309) * t273) * qJD(3) + t26 * t413 + (t26 * t104 + t35 * t414 + t427 * t61) * t268 + (t26 * t102 + t34 * t414 + t427 * t60) * t267 + (-t149 * t234 + t151 * t235 + t81 * t267 + t83 * t268 + t528) * t46) * m(4) + (g(1) * t223 + g(2) * t222 - g(3) * t244 + t55 * t331 + (-(-t222 * t267 - t223 * t268) * qJD(2) - t328) * (qJD(2) * t331 + qJD(1)) + (t244 * qJD(2) ^ 2 + qJDD(2) * t243) * t563) * m(3) + ((-t160 * t562 + t161 * t265) * t553 + (-t192 * t562 + t193 * t265) * t554 + t551) * t267 / 0.2e1 - ((t160 * t266 - t161 * t562) * t553 + (t192 * t266 - t193 * t562) * t554 + t552) * t268 / 0.2e1 + (t267 * t518 - t268 * t519 + t555) * t369; (t504 * t271 - t273 * t572) * t495 + (t503 * t271 - t273 * t573) * t494 + (t215 * t508 - t268 * t506 - t507 * t529) * t493 + (t576 * t273 + (t267 * t547 + t268 * t546) * t271 + (t273 * t504 + t526) * qJD(2)) * t492 + (t212 * t507 + t213 * t508 - t267 * t506) * t491 + (t577 * t273 + (t267 * t549 + t268 * t548) * t271 + (t273 * t503 + t525) * qJD(2)) * t490 + (t271 * t505 - t273 * t516) * t489 - (t154 * t541 + t155 * t542 + t518 * t234 + t519 * t235 + t516 * t238 + t517 * t402) * t273 / 0.2e1 + t552 * t440 / 0.2e1 + t551 * t436 / 0.2e1 + t545 * t405 / 0.2e1 + (t517 * t273 + (t267 * t519 + t268 * t518) * t271 + (t271 * t516 + t273 * t505) * qJD(2)) * t369 + (t559 * t273 + (-t270 * t507 + t508 * t272) * t271) * t368 + (g(1) * t423 - g(2) * t424 + (qJD(2) * t277 - t36 * t480 + t37 * t481 - t429 * t6 + t464 * t7) * t273 - (t31 * t423 - t362 * t37) * t235 - (t31 * t424 + t36 * t362) * t234 - (t36 * t423 + t37 * t424) * t402 + (-g(3) * (-t270 * t550 - t477) + t298 * qJD(2) + (-t36 * t428 - t422 * t6 + t300) * t268 + (t37 * t428 + t422 * t7 - t500) * t267) * t271) * m(5) + ((t35 * t102 - t34 * t104 - t60 * t83 + t61 * t81 + (t309 + (t267 * t61 - t268 * t60) * t187) * qJD(2)) * t273 + (t61 * (-qJD(2) * t102 + t113 * t267) + t60 * (qJD(2) * t104 - t113 * t268) + t26 * t334 + t46 * (-t267 * t83 + t268 * t81) + (t267 * t35 - t268 * t34) * t187) * t271 - t61 * (t129 * t402 + t233 * t235) - t60 * (-t131 * t402 - t233 * t234) - t46 * (t129 * t234 - t131 * t235) - g(1) * t131 - g(2) * t129 - g(3) * t233) * m(4) + t555 * t404 / 0.2e1; ((-t5 + t502 * qJD(2) - t37 * (-t235 + t380) - t36 * (t234 - t379) + g(3)) * t273 + ((-t234 * t267 + t235 * t268 + qJD(2)) * t31 + t267 * t6 + t268 * t7 - t501) * t271) * m(5);];
tau = t1;
