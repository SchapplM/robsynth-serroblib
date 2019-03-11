% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:29
% EndTime: 2019-03-08 23:06:06
% DurationCPUTime: 22.32s
% Computational Cost: add. (12717->768), mult. (28455->1051), div. (0->0), fcn. (21998->18), ass. (0->359)
t364 = sin(qJ(2));
t357 = sin(pkin(6));
t450 = qJD(1) * t357;
t426 = t364 * t450;
t363 = sin(qJ(3));
t444 = qJD(3) * t363;
t538 = pkin(3) * t444 - t426;
t576 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t362 = sin(qJ(4));
t366 = cos(qJ(4));
t367 = cos(qJ(3));
t305 = t362 * t363 - t366 * t367;
t381 = t305 * qJD(4);
t237 = -qJD(3) * t305 - t381;
t306 = t362 * t367 + t363 * t366;
t382 = t306 * qJD(4);
t238 = qJD(3) * t306 + t382;
t587 = pkin(4) * t238 - qJ(5) * t237 - qJD(5) * t306 + t538;
t369 = -pkin(9) - pkin(8);
t427 = qJD(3) * t369;
t310 = t363 * t427;
t311 = t367 * t427;
t368 = cos(qJ(2));
t425 = t368 * t450;
t326 = t369 * t363;
t327 = t369 * t367;
t539 = t366 * t326 + t327 * t362;
t541 = qJD(4) * t539 + t305 * t425 + t310 * t366 + t311 * t362;
t353 = pkin(12) + qJ(6);
t347 = sin(t353);
t348 = cos(t353);
t356 = sin(pkin(12));
t358 = cos(pkin(12));
t398 = -mrSges(6,1) * t358 + mrSges(6,2) * t356;
t584 = mrSges(7,1) * t348 - mrSges(7,2) * t347 + mrSges(5,1) - t398;
t550 = -t356 * t541 + t358 * t587;
t549 = t356 * t587 + t358 * t541;
t341 = pkin(5) * t358 + pkin(4);
t355 = qJ(3) + qJ(4);
t349 = sin(t355);
t350 = cos(t355);
t360 = -pkin(10) - qJ(5);
t583 = (-m(6) * pkin(4) - m(7) * t341 - t584) * t350 + (-m(6) * qJ(5) + m(7) * t360 + t576) * t349;
t476 = t237 * t358;
t582 = pkin(5) * t238 - pkin(10) * t476 + t550;
t477 = t237 * t356;
t581 = pkin(10) * t477 - t549;
t298 = t306 * qJD(2);
t354 = qJD(3) + qJD(4);
t249 = t298 * t358 + t354 * t356;
t404 = -t298 * t356 + t358 * t354;
t494 = mrSges(5,3) * t298;
t540 = mrSges(5,1) * t354 + mrSges(6,1) * t404 - mrSges(6,2) * t249 - t494;
t361 = sin(qJ(6));
t365 = cos(qJ(6));
t147 = t249 * t365 + t361 * t404;
t565 = -t249 * t361 + t365 * t404;
t76 = -mrSges(7,1) * t565 + mrSges(7,2) * t147;
t580 = t76 - t540;
t445 = qJD(2) * t367;
t447 = qJD(2) * t363;
t297 = t362 * t447 - t366 * t445;
t386 = t356 * t361 - t358 * t365;
t182 = t386 * t297;
t291 = t386 * qJD(6);
t579 = t182 + t291;
t314 = qJD(2) * pkin(8) + t426;
t408 = pkin(9) * qJD(2) + t314;
t359 = cos(pkin(6));
t449 = qJD(1) * t359;
t424 = t363 * t449;
t243 = t367 * t408 + t424;
t235 = t366 * t243;
t336 = t367 * t449;
t242 = -t363 * t408 + t336;
t236 = qJD(3) * pkin(3) + t242;
t137 = t236 * t362 + t235;
t126 = qJ(5) * t354 + t137;
t344 = pkin(3) * t367 + pkin(2);
t284 = -qJD(2) * t344 - t425;
t169 = pkin(4) * t297 - qJ(5) * t298 + t284;
t72 = -t126 * t356 + t358 * t169;
t51 = pkin(5) * t297 - pkin(10) * t249 + t72;
t73 = t358 * t126 + t356 * t169;
t59 = pkin(10) * t404 + t73;
t19 = -t361 * t59 + t365 * t51;
t574 = t19 * mrSges(7,1);
t20 = t361 * t51 + t365 * t59;
t573 = t20 * mrSges(7,2);
t286 = qJD(6) + t297;
t556 = t404 * Ifges(6,6);
t560 = t249 * Ifges(6,5);
t572 = t147 * Ifges(7,5) + Ifges(7,6) * t565 + t297 * Ifges(6,3) + t286 * Ifges(7,3) + t556 + t560;
t304 = t356 * t365 + t358 * t361;
t181 = t304 * t297;
t292 = t304 * qJD(6);
t571 = t181 + t292;
t397 = mrSges(6,1) * t356 + mrSges(6,2) * t358;
t570 = -m(7) * pkin(5) * t356 - t347 * mrSges(7,1) - t348 * mrSges(7,2) - mrSges(5,3) - t397;
t481 = sin(pkin(11));
t411 = t481 * t364;
t482 = cos(pkin(11));
t412 = t482 * t368;
t290 = -t359 * t411 + t412;
t414 = t357 * t481;
t569 = -t290 * t363 + t367 * t414;
t468 = t357 * t364;
t293 = t359 * t367 - t363 * t468;
t410 = t481 * t368;
t413 = t482 * t364;
t288 = t359 * t413 + t410;
t415 = t357 * t482;
t228 = t288 * t349 + t350 * t415;
t229 = t288 * t350 - t349 * t415;
t568 = t228 * t584 + t229 * t576;
t230 = t290 * t349 - t350 * t414;
t231 = t290 * t350 + t349 * t414;
t567 = t230 * t584 + t231 * t576;
t274 = t349 * t468 - t359 * t350;
t275 = t349 * t359 + t350 * t468;
t566 = t274 * t584 + t275 * t576;
t440 = qJD(2) * qJD(3);
t312 = qJDD(2) * t367 - t363 * t440;
t313 = qJDD(2) * t363 + t367 * t440;
t173 = -qJD(2) * t381 + t312 * t362 + t313 * t366;
t352 = qJDD(3) + qJDD(4);
t150 = -t173 * t356 + t352 * t358;
t151 = t173 * t358 + t352 * t356;
t46 = qJD(6) * t565 + t150 * t361 + t151 * t365;
t525 = t46 / 0.2e1;
t47 = -qJD(6) * t147 + t150 * t365 - t151 * t361;
t524 = t47 / 0.2e1;
t516 = t150 / 0.2e1;
t515 = t151 / 0.2e1;
t174 = qJD(2) * t382 - t366 * t312 + t313 * t362;
t172 = qJDD(6) + t174;
t514 = t172 / 0.2e1;
t513 = t174 / 0.2e1;
t564 = t312 / 0.2e1;
t563 = t313 / 0.2e1;
t226 = pkin(4) * t305 - qJ(5) * t306 - t344;
t256 = t326 * t362 - t327 * t366;
t135 = t356 * t226 + t358 * t256;
t472 = t306 * t356;
t109 = -pkin(10) * t472 + t135;
t134 = t358 * t226 - t256 * t356;
t471 = t306 * t358;
t96 = pkin(5) * t305 - pkin(10) * t471 + t134;
t43 = t109 * t365 + t361 * t96;
t562 = -qJD(6) * t43 + t581 * t361 + t365 * t582;
t42 = -t109 * t361 + t365 * t96;
t561 = qJD(6) * t42 + t361 * t582 - t581 * t365;
t559 = t354 * Ifges(5,5);
t558 = t354 * Ifges(5,6);
t557 = t367 * Ifges(4,2);
t503 = pkin(3) * t362;
t340 = qJ(5) + t503;
t300 = (-pkin(10) - t340) * t356;
t351 = t358 * pkin(10);
t469 = t340 * t358;
t301 = t351 + t469;
t207 = t300 * t365 - t301 * t361;
t441 = qJD(4) * t366;
t432 = pkin(3) * t441;
t335 = qJD(5) + t432;
t473 = t297 * t358;
t402 = t298 * pkin(5) + pkin(10) * t473;
t234 = t362 * t243;
t139 = t242 * t366 - t234;
t221 = pkin(4) * t298 + qJ(5) * t297;
t435 = pkin(3) * t447;
t187 = t221 + t435;
t82 = -t139 * t356 + t358 * t187;
t65 = t402 + t82;
t474 = t297 * t356;
t437 = pkin(10) * t474;
t83 = t358 * t139 + t356 * t187;
t71 = t437 + t83;
t554 = qJD(6) * t207 - t335 * t386 - t361 * t65 - t365 * t71;
t208 = t300 * t361 + t301 * t365;
t553 = -qJD(6) * t208 - t304 * t335 + t361 * t71 - t365 * t65;
t319 = t360 * t356;
t480 = qJ(5) * t358;
t320 = t351 + t480;
t244 = t319 * t365 - t320 * t361;
t136 = t236 * t366 - t234;
t88 = -t136 * t356 + t358 * t221;
t66 = t402 + t88;
t89 = t358 * t136 + t356 * t221;
t74 = t437 + t89;
t552 = -qJD(5) * t386 + qJD(6) * t244 - t361 * t66 - t365 * t74;
t245 = t319 * t361 + t320 * t365;
t551 = -qJD(5) * t304 - qJD(6) * t245 + t361 * t74 - t365 * t66;
t159 = mrSges(5,1) * t352 - mrSges(5,3) * t173;
t77 = -t150 * mrSges(6,1) + t151 * mrSges(6,2);
t548 = t77 - t159;
t125 = -pkin(4) * t354 + qJD(5) - t136;
t547 = t125 * t397;
t448 = qJD(2) * t357;
t417 = qJD(1) * t448;
t329 = t368 * t417;
t439 = qJDD(1) * t357;
t280 = t364 * t439 + t329;
t272 = qJDD(2) * pkin(8) + t280;
t438 = qJDD(1) * t359;
t155 = qJD(3) * t336 + t367 * t272 - t314 * t444 + t363 * t438;
t261 = t314 * t367 + t424;
t156 = -qJD(3) * t261 - t272 * t363 + t367 * t438;
t537 = t155 * t367 - t156 * t363;
t536 = m(7) + m(6) + m(5);
t124 = qJDD(3) * pkin(3) - pkin(9) * t313 + t156;
t130 = pkin(9) * t312 + t155;
t442 = qJD(4) * t362;
t36 = t362 * t124 + t366 * t130 + t236 * t441 - t243 * t442;
t31 = qJ(5) * t352 + qJD(5) * t354 + t36;
t328 = t364 * t417;
t279 = t368 * t439 - t328;
t271 = -qJDD(2) * pkin(2) - t279;
t227 = -pkin(3) * t312 + t271;
t62 = pkin(4) * t174 - qJ(5) * t173 - qJD(5) * t298 + t227;
t13 = -t31 * t356 + t358 * t62;
t6 = pkin(5) * t174 - pkin(10) * t151 + t13;
t14 = t358 * t31 + t356 * t62;
t7 = pkin(10) * t150 + t14;
t2 = qJD(6) * t19 + t361 * t6 + t365 * t7;
t3 = -qJD(6) * t20 - t361 * t7 + t365 * t6;
t532 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t431 = m(4) * pkin(8) + mrSges(4,3);
t531 = mrSges(3,2) - t431 + t570;
t400 = -mrSges(4,1) * t367 + mrSges(4,2) * t363;
t380 = m(4) * pkin(2) - t400;
t530 = mrSges(3,1) + t380 - t583;
t483 = t298 * Ifges(5,4);
t189 = -t297 * Ifges(5,2) + t483 + t558;
t283 = Ifges(5,4) * t297;
t190 = t298 * Ifges(5,1) - t283 + t559;
t37 = t124 * t366 - t362 * t130 - t236 * t442 - t243 * t441;
t32 = -pkin(4) * t352 + qJDD(5) - t37;
t22 = -pkin(5) * t150 + t32;
t390 = Ifges(6,5) * t358 - Ifges(6,6) * t356;
t490 = Ifges(6,4) * t358;
t392 = -Ifges(6,2) * t356 + t490;
t491 = Ifges(6,4) * t356;
t394 = Ifges(6,1) * t358 - t491;
t504 = t358 / 0.2e1;
t420 = (t249 * Ifges(6,1) + Ifges(6,4) * t404 + t297 * Ifges(6,5)) * t504;
t421 = -t356 * (t249 * Ifges(6,4) + Ifges(6,2) * t404 + t297 * Ifges(6,6)) / 0.2e1;
t486 = t14 * t358;
t488 = t136 * mrSges(5,3);
t506 = t298 / 0.2e1;
t508 = t297 / 0.2e1;
t509 = -t297 / 0.2e1;
t510 = t286 / 0.2e1;
t511 = -t286 / 0.2e1;
t517 = t147 / 0.2e1;
t518 = -t147 / 0.2e1;
t519 = t565 / 0.2e1;
t52 = t151 * Ifges(6,4) + t150 * Ifges(6,2) + t174 * Ifges(6,6);
t520 = -t565 / 0.2e1;
t142 = Ifges(7,4) * t565;
t69 = Ifges(7,1) * t147 + Ifges(7,5) * t286 + t142;
t521 = t69 / 0.2e1;
t485 = t147 * Ifges(7,4);
t68 = Ifges(7,2) * t565 + t286 * Ifges(7,6) + t485;
t522 = t68 / 0.2e1;
t523 = Ifges(6,1) * t515 + Ifges(6,4) * t516 + Ifges(6,5) * t513;
t526 = Ifges(7,1) * t525 + Ifges(7,4) * t524 + Ifges(7,5) * t514;
t527 = Ifges(7,4) * t525 + Ifges(7,2) * t524 + Ifges(7,6) * t514;
t92 = -pkin(5) * t404 + t125;
t529 = (-Ifges(7,1) * t291 - Ifges(7,4) * t292) * t517 + (-Ifges(7,4) * t291 - Ifges(7,2) * t292) * t519 + (-Ifges(7,5) * t291 - Ifges(7,6) * t292) * t510 - t404 * (Ifges(6,6) * t298 - t297 * t392) / 0.2e1 - t73 * (-mrSges(6,2) * t298 + mrSges(6,3) * t474) - t72 * (mrSges(6,1) * t298 + mrSges(6,3) * t473) - t284 * (mrSges(5,1) * t298 - mrSges(5,2) * t297) + (-Ifges(5,2) * t298 + t190 - t283) * t508 - t298 * t574 + t52 * t504 + t22 * (mrSges(7,1) * t386 + mrSges(7,2) * t304) + (Ifges(7,5) * t304 - Ifges(7,6) * t386) * t514 + (Ifges(7,4) * t304 - Ifges(7,2) * t386) * t524 + (Ifges(7,1) * t304 - Ifges(7,4) * t386) * t525 - t386 * t527 - (-Ifges(5,1) * t297 - t483 + t572) * t298 / 0.2e1 + mrSges(6,3) * t486 - t297 * t488 + Ifges(5,3) * t352 - t354 * (-Ifges(5,5) * t297 - Ifges(5,6) * t298) / 0.2e1 - t181 * t68 / 0.2e1 - t182 * t69 / 0.2e1 + Ifges(5,5) * t173 - Ifges(5,6) * t174 + t37 * mrSges(5,1) - t36 * mrSges(5,2) + t297 * t420 + t297 * t421 + t298 * t573 + (t19 * t579 - t2 * t386 - t571 * t20 - t3 * t304) * mrSges(7,3) + (t571 * mrSges(7,1) - mrSges(7,2) * t579) * t92 + t189 * t506 + (Ifges(6,3) * t298 - t297 * t390) * t509 + (Ifges(7,5) * t182 + Ifges(7,6) * t181 + Ifges(7,3) * t298) * t511 + (Ifges(6,5) * t356 + Ifges(6,6) * t358) * t513 + (Ifges(6,1) * t356 + t490) * t515 + (Ifges(6,2) * t358 + t491) * t516 + (Ifges(7,1) * t182 + Ifges(7,4) * t181 + Ifges(7,5) * t298) * t518 + (Ifges(7,4) * t182 + Ifges(7,2) * t181 + Ifges(7,6) * t298) * t520 - t291 * t521 - t292 * t522 + t356 * t523 + t304 * t526 - t249 * (Ifges(6,5) * t298 - t297 * t394) / 0.2e1 + t32 * t398 + t297 * t547;
t528 = m(6) * t125 + m(7) * t92 + t580;
t370 = qJD(2) ^ 2;
t502 = pkin(3) * t366;
t500 = g(3) * t357;
t493 = Ifges(4,4) * t363;
t492 = Ifges(4,4) * t367;
t489 = t13 * t356;
t487 = t137 * mrSges(5,3);
t467 = t357 * t368;
t465 = -t228 * t341 - t229 * t360;
t464 = -t230 * t341 - t231 * t360;
t457 = -t274 * t341 - t275 * t360;
t446 = qJD(2) * t364;
t443 = qJD(3) * t367;
t436 = Ifges(7,5) * t46 + Ifges(7,6) * t47 + Ifges(7,3) * t172;
t430 = mrSges(4,3) * t447;
t429 = mrSges(4,3) * t445;
t423 = t357 * t446;
t422 = t368 * t448;
t15 = -t47 * mrSges(7,1) + t46 * mrSges(7,2);
t409 = t440 / 0.2e1;
t407 = -t228 * pkin(4) + t229 * qJ(5);
t406 = -t230 * pkin(4) + qJ(5) * t231;
t405 = -t274 * pkin(4) + qJ(5) * t275;
t138 = t242 * t362 + t235;
t403 = t569 * pkin(3);
t393 = t493 + t557;
t391 = Ifges(4,5) * t367 - Ifges(4,6) * t363;
t389 = t486 - t489;
t388 = -t356 * t72 + t358 * t73;
t294 = t359 * t363 + t367 * t468;
t196 = t293 * t362 + t294 * t366;
t164 = -t196 * t356 - t358 * t467;
t165 = t196 * t358 - t356 * t467;
t84 = t164 * t365 - t165 * t361;
t85 = t164 * t361 + t165 * t365;
t387 = t366 * t293 - t294 * t362;
t385 = t293 * pkin(3);
t315 = -qJD(2) * pkin(2) - t425;
t384 = t315 * (mrSges(4,1) * t363 + mrSges(4,2) * t367);
t383 = t363 * (Ifges(4,1) * t367 - t493);
t378 = -g(1) * t230 - g(2) * t228 - g(3) * t274;
t376 = -t288 * t363 - t367 * t415;
t374 = t376 * pkin(3);
t153 = qJD(4) * t256 + t310 * t362 - t366 * t311;
t345 = Ifges(4,4) * t445;
t343 = -pkin(4) - t502;
t323 = -qJD(3) * mrSges(4,2) + t429;
t322 = qJD(3) * mrSges(4,1) - t430;
t318 = -t341 - t502;
t308 = t400 * qJD(2);
t296 = Ifges(4,1) * t447 + Ifges(4,5) * qJD(3) + t345;
t295 = Ifges(4,6) * qJD(3) + qJD(2) * t393;
t289 = t359 * t410 + t413;
t287 = -t359 * t412 + t411;
t282 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t313;
t281 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t312;
t273 = pkin(5) * t474;
t262 = -mrSges(5,2) * t354 - mrSges(5,3) * t297;
t260 = -t314 * t363 + t336;
t246 = -mrSges(4,1) * t312 + mrSges(4,2) * t313;
t241 = qJD(3) * t293 + t367 * t422;
t240 = -qJD(3) * t294 - t363 * t422;
t223 = mrSges(5,1) * t297 + mrSges(5,2) * t298;
t200 = t386 * t306;
t199 = t304 * t306;
t183 = pkin(5) * t472 - t539;
t176 = mrSges(6,1) * t297 - mrSges(6,3) * t249;
t175 = -mrSges(6,2) * t297 + mrSges(6,3) * t404;
t160 = -mrSges(5,2) * t352 - mrSges(5,3) * t174;
t115 = mrSges(7,1) * t286 - mrSges(7,3) * t147;
t114 = -mrSges(7,2) * t286 + mrSges(7,3) * t565;
t112 = t138 - t273;
t106 = t137 - t273;
t99 = pkin(5) * t477 + t153;
t93 = mrSges(5,1) * t174 + mrSges(5,2) * t173;
t91 = qJD(4) * t196 - t366 * t240 + t241 * t362;
t90 = qJD(4) * t387 + t240 * t362 + t241 * t366;
t87 = mrSges(6,1) * t174 - mrSges(6,3) * t151;
t86 = -mrSges(6,2) * t174 + mrSges(6,3) * t150;
t81 = -t237 * t304 + t291 * t306;
t80 = -t237 * t386 - t292 * t306;
t79 = t356 * t423 + t358 * t90;
t78 = -t356 * t90 + t358 * t423;
t28 = -mrSges(7,2) * t172 + mrSges(7,3) * t47;
t27 = mrSges(7,1) * t172 - mrSges(7,3) * t46;
t17 = -qJD(6) * t85 - t361 * t79 + t365 * t78;
t16 = qJD(6) * t84 + t361 * t78 + t365 * t79;
t1 = [m(2) * qJDD(1) + t16 * t114 + t17 * t115 + t196 * t160 + t164 * t87 + t165 * t86 + t79 * t175 + t78 * t176 + t240 * t322 + t241 * t323 + t90 * t262 + t84 * t27 + t85 * t28 + t294 * t281 + t293 * t282 + t580 * t91 - (t15 + t548) * t387 + (-m(2) - m(3) - m(4) - t536) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t370 - t246 - t93) * t368 + (-mrSges(3,1) * t370 - mrSges(3,2) * qJDD(2) + (t223 + t308) * qJD(2)) * t364) * t357 + m(7) * (t16 * t20 + t17 * t19 + t2 * t85 - t22 * t387 + t3 * t84 + t91 * t92) + m(6) * (t125 * t91 + t13 * t164 + t14 * t165 - t32 * t387 + t72 * t78 + t73 * t79) + m(5) * (-t136 * t91 + t137 * t90 + t387 * t37 + t196 * t36 + (-t227 * t368 + t284 * t446) * t357) + m(4) * (t155 * t294 + t156 * t293 + t240 * t260 + t241 * t261 + (-t271 * t368 + t315 * t446) * t357) + m(3) * (qJDD(1) * t359 ^ 2 + (t279 * t368 + t280 * t364) * t357); (-pkin(2) * t271 - (t315 * t364 + (-t260 * t363 + t261 * t367) * t368) * t450) * m(4) + (t364 * t500 - t280 + t329) * mrSges(3,2) - t540 * t153 + t541 * t262 + (-t368 * t500 + t279 + t328) * mrSges(3,1) + (-t487 - t558 / 0.2e1 + t72 * mrSges(6,1) - t73 * mrSges(6,2) + t284 * mrSges(5,1) + t574 - t573 - t189 / 0.2e1 - Ifges(5,4) * t506 + Ifges(6,3) * t508 - Ifges(5,2) * t509 + Ifges(7,3) * t510 + Ifges(7,5) * t517 + Ifges(7,6) * t519 + t556 / 0.2e1 + t560 / 0.2e1 + t572 / 0.2e1) * t238 + (-Ifges(7,5) * t200 - Ifges(7,6) * t199) * t514 + (-Ifges(7,1) * t200 - Ifges(7,4) * t199) * t525 + (-t19 * t80 - t199 * t2 + t20 * t81 + t200 * t3) * mrSges(7,3) + (-Ifges(7,4) * t200 - Ifges(7,2) * t199) * t524 + t22 * (mrSges(7,1) * t199 - mrSges(7,2) * t200) + (Ifges(6,5) * t151 + Ifges(6,6) * t150 + Ifges(6,3) * t174 + t436) * t305 / 0.2e1 + (t227 * mrSges(5,2) - t37 * mrSges(5,3) + Ifges(5,1) * t173 - Ifges(5,4) * t174 + Ifges(5,5) * t352 + t32 * t397 + t390 * t513 + t392 * t516 + t394 * t515 + (m(5) * t136 - t528) * t425) * t306 + (-t136 * t153 + t137 * t541 - t227 * t344 + t256 * t36 + t284 * t538 + t37 * t539) * m(5) + (t125 * t153 + t13 * t134 + t135 * t14 - t32 * t539 + t549 * t73 + t550 * t72) * m(6) - t548 * t539 + (m(4) * ((-t260 * t367 - t261 * t363) * qJD(3) + t537) - t322 * t443 - t323 * t444 - t363 * t282 + t367 * t281) * pkin(8) + (-t260 * t443 - t261 * t444 + t537) * mrSges(4,3) + t538 * t223 + t561 * t114 + (t183 * t22 + t19 * t562 + t2 * t43 + t20 * t561 + t3 * t42 + t92 * t99) * m(7) + t562 * t115 + t549 * t175 + (-t13 * t471 - t14 * t472 - t476 * t72 - t477 * t73) * mrSges(6,3) + t492 * t563 + (-t536 * (-t289 * t344 - t290 * t369) + t531 * t290 + t530 * t289) * g(1) + (-t536 * (-t287 * t344 - t288 * t369) + t531 * t288 + t530 * t287) * g(2) + t550 * t176 + (t384 + t391 * qJD(3) / 0.2e1) * qJD(3) + (Ifges(7,5) * t80 + Ifges(7,6) * t81) * t510 + (Ifges(7,1) * t80 + Ifges(7,4) * t81) * t517 + (t227 * mrSges(5,1) + t13 * mrSges(6,1) - t14 * mrSges(6,2) - t36 * mrSges(5,3) - Ifges(5,4) * t173 + Ifges(6,5) * t515 + Ifges(7,5) * t525 + Ifges(5,2) * t174 - Ifges(5,6) * t352 + Ifges(6,6) * t516 + Ifges(7,6) * t524 + Ifges(6,3) * t513 + Ifges(7,3) * t514 + t532) * t305 + (Ifges(7,4) * t80 + Ifges(7,2) * t81) * t519 + (-t488 + t559 / 0.2e1 + t284 * mrSges(5,2) + t190 / 0.2e1 + t547 + t420 + t421 + Ifges(5,1) * t506 + t390 * t508 + Ifges(5,4) * t509 + t404 * t392 / 0.2e1 + t249 * t394 / 0.2e1) * t237 + Ifges(3,3) * qJDD(2) - (t364 * t431 + t368 * t380) * t500 + qJDD(3) * (Ifges(4,5) * t363 + Ifges(4,6) * t367) + (Ifges(4,4) * t563 + Ifges(4,2) * t564 - t323 * t425 + t492 * t409) * t367 + (Ifges(4,1) * t313 + Ifges(4,4) * t564 + t322 * t425 - t409 * t557) * t363 - t344 * t93 - t52 * t472 / 0.2e1 - t295 * t444 / 0.2e1 + t296 * t443 / 0.2e1 + t256 * t160 - pkin(2) * t246 - t308 * t426 + t183 * t15 + t134 * t87 + t135 * t86 + t99 * t76 + t92 * (-mrSges(7,1) * t81 + mrSges(7,2) * t80) + t42 * t27 + t43 * t28 + ((t583 * t368 + (t369 * t536 + t570) * t364) * t357 - t536 * t344 * t467) * g(3) + t383 * t409 + t393 * t564 + t80 * t521 + t81 * t522 + t471 * t523 - t200 * t526 - t199 * t527 + t271 * t400; t540 * t138 - (-Ifges(4,2) * t447 + t296 + t345) * t445 / 0.2e1 + t529 + (-m(6) * (t385 + t405) - m(7) * (t385 + t457) - m(5) * t385 - mrSges(4,1) * t293 + mrSges(4,2) * t294 + t566) * g(3) + (-m(6) * (t374 + t407) - m(7) * (t374 + t465) - m(5) * t374 - t376 * mrSges(4,1) - (-t288 * t367 + t363 * t415) * mrSges(4,2) + t568) * g(2) + (-m(6) * (t403 + t406) - m(7) * (t403 + t464) - t569 * mrSges(4,1) - (-t290 * t367 - t363 * t414) * mrSges(4,2) - m(5) * t403 + t567) * g(1) + t159 * t502 + t160 * t503 + t553 * t115 + (-t112 * t92 + t19 * t553 + t2 * t208 + t20 * t554 + t207 * t3 + t22 * t318) * m(7) + t554 * t114 + (-t323 + t429) * t260 + (-t335 * t356 - t82) * t176 + (-t139 + t432) * t262 + (t335 * t358 - t83) * t175 + t86 * t469 + t298 * t487 + t528 * pkin(3) * t442 + ((t36 * t362 + t366 * t37 + (-t136 * t362 + t137 * t366) * qJD(4)) * pkin(3) + t136 * t138 - t137 * t139 - t284 * t435) * m(5) + (-t125 * t138 + t32 * t343 + t335 * t388 + t340 * t389 - t72 * t82 - t73 * t83) * m(6) + (t322 + t430) * t261 - t370 * t383 / 0.2e1 + Ifges(4,3) * qJDD(3) - mrSges(6,3) * t489 + t343 * t77 + Ifges(4,5) * t313 + t318 * t15 + Ifges(4,6) * t312 + t295 * t447 / 0.2e1 - t391 * t440 / 0.2e1 - t223 * t435 + t207 * t27 + t208 * t28 - t155 * mrSges(4,2) + t156 * mrSges(4,1) - t112 * t76 - qJD(2) * t384 - t356 * t340 * t87; t551 * t115 + t552 * t114 + (t540 + t494) * t137 + (-mrSges(6,3) * t13 - qJ(5) * t87 - qJD(5) * t176) * t356 + t568 * g(2) + t567 * g(1) + t566 * g(3) + (qJD(5) * t358 - t89) * t175 - t341 * t15 - t136 * t262 + t244 * t27 + t245 * t28 - t88 * t176 - t106 * t76 - pkin(4) * t77 + t86 * t480 + (-t464 * g(1) - t465 * g(2) - t457 * g(3) - t106 * t92 + t19 * t551 + t2 * t245 + t20 * t552 - t22 * t341 + t244 * t3) * m(7) + (-pkin(4) * t32 - g(1) * t406 - g(2) * t407 - g(3) * t405 + qJ(5) * t389 + qJD(5) * t388 - t125 * t137 - t72 * t88 - t73 * t89) * m(6) + t529; -t565 * t114 + t147 * t115 - t404 * t175 + t249 * t176 + t15 + t77 + (t147 * t19 - t20 * t565 + t22 + t378) * m(7) + (t249 * t72 - t404 * t73 + t32 + t378) * m(6); -t92 * (mrSges(7,1) * t147 + mrSges(7,2) * t565) + (Ifges(7,1) * t565 - t485) * t518 + t68 * t517 + (Ifges(7,5) * t565 - Ifges(7,6) * t147) * t511 - t19 * t114 + t20 * t115 - g(1) * ((-t231 * t347 + t289 * t348) * mrSges(7,1) + (-t231 * t348 - t289 * t347) * mrSges(7,2)) - g(2) * ((-t229 * t347 + t287 * t348) * mrSges(7,1) + (-t229 * t348 - t287 * t347) * mrSges(7,2)) - g(3) * ((-t275 * t347 - t348 * t467) * mrSges(7,1) + (-t275 * t348 + t347 * t467) * mrSges(7,2)) + (t147 * t20 + t19 * t565) * mrSges(7,3) + t436 + (-Ifges(7,2) * t147 + t142 + t69) * t520 + t532;];
tau  = t1;
