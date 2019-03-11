% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:52
% EndTime: 2019-03-10 01:06:29
% DurationCPUTime: 17.64s
% Computational Cost: add. (18739->682), mult. (44841->902), div. (0->0), fcn. (32554->8), ass. (0->315)
t381 = cos(qJ(2));
t510 = -pkin(8) - pkin(7);
t361 = t510 * t381;
t348 = qJD(1) * t361;
t376 = sin(qJ(3));
t326 = t376 * t348;
t377 = sin(qJ(2));
t359 = t510 * t377;
t347 = qJD(1) * t359;
t332 = qJD(2) * pkin(2) + t347;
t380 = cos(qJ(3));
t279 = t332 * t380 + t326;
t373 = qJD(2) + qJD(3);
t251 = -pkin(3) * t373 - t279;
t342 = t376 * t381 + t380 * t377;
t324 = t342 * qJD(1);
t375 = sin(qJ(4));
t379 = cos(qJ(4));
t295 = -t324 * t375 + t373 * t379;
t200 = -pkin(4) * t295 + t251;
t296 = t324 * t379 + t373 * t375;
t374 = sin(qJ(5));
t378 = cos(qJ(5));
t405 = t378 * t295 - t296 * t374;
t115 = -pkin(5) * t405 + qJD(6) + t200;
t213 = t295 * t374 + t296 * t378;
t427 = qJD(1) * t381;
t428 = qJD(1) * t377;
t323 = -t376 * t428 + t380 * t427;
t319 = qJD(4) - t323;
t309 = qJD(5) + t319;
t368 = -pkin(2) * t381 - pkin(1);
t357 = qJD(1) * t368;
t240 = -t323 * pkin(3) - t324 * pkin(9) + t357;
t431 = t380 * t348;
t280 = t376 * t332 - t431;
t252 = pkin(9) * t373 + t280;
t158 = t379 * t240 - t252 * t375;
t129 = -pkin(10) * t296 + t158;
t113 = pkin(4) * t319 + t129;
t159 = t240 * t375 + t252 * t379;
t130 = pkin(10) * t295 + t159;
t124 = t374 * t130;
t45 = t378 * t113 - t124;
t568 = qJ(6) * t213;
t32 = t45 - t568;
t29 = pkin(5) * t309 + t32;
t126 = t378 * t130;
t46 = t113 * t374 + t126;
t526 = qJ(6) * t405;
t33 = t46 + t526;
t485 = -t309 / 0.2e1;
t496 = t213 / 0.2e1;
t497 = -t213 / 0.2e1;
t288 = t373 * t342;
t270 = t288 * qJD(1);
t340 = t376 * t377 - t380 * t381;
t287 = t373 * t340;
t269 = t287 * qJD(1);
t196 = qJD(4) * t295 - t269 * t379;
t197 = -qJD(4) * t296 + t269 * t375;
t73 = qJD(5) * t405 + t196 * t378 + t197 * t374;
t426 = qJD(2) * t377;
t418 = pkin(2) * t426;
t554 = qJD(1) * t418;
t169 = pkin(3) * t270 + pkin(9) * t269 + t554;
t415 = qJD(2) * t510;
t404 = qJD(1) * t415;
t333 = t377 * t404;
t391 = t381 * t404;
t189 = qJD(3) * t279 + t380 * t333 + t376 * t391;
t51 = -qJD(4) * t159 + t379 * t169 - t189 * t375;
t26 = pkin(4) * t270 - pkin(10) * t196 + t51;
t424 = qJD(4) * t379;
t425 = qJD(4) * t375;
t50 = t375 * t169 + t379 * t189 + t240 * t424 - t252 * t425;
t31 = pkin(10) * t197 + t50;
t9 = -qJD(5) * t46 + t378 * t26 - t31 * t374;
t2 = pkin(5) * t270 - qJ(6) * t73 - qJD(6) * t213 + t9;
t74 = -qJD(5) * t213 - t196 * t374 + t197 * t378;
t422 = qJD(5) * t378;
t423 = qJD(5) * t374;
t8 = t113 * t422 - t130 * t423 + t374 * t26 + t378 * t31;
t4 = qJ(6) * t74 + qJD(6) * t405 + t8;
t515 = t9 * mrSges(6,1) + t2 * mrSges(7,1) - t8 * mrSges(6,2) - t4 * mrSges(7,2);
t544 = Ifges(6,6) + Ifges(7,6);
t546 = Ifges(6,5) + Ifges(7,5);
t572 = Ifges(6,3) + Ifges(7,3);
t520 = t270 * t572 + t544 * t74 + t546 * t73;
t545 = Ifges(6,2) + Ifges(7,2);
t548 = Ifges(6,1) + Ifges(7,1);
t549 = -t405 / 0.2e1;
t547 = Ifges(6,4) + Ifges(7,4);
t580 = t547 * t405;
t566 = t213 * t548 + t546 * t309 + t580;
t579 = t213 * t547;
t567 = t309 * t544 + t405 * t545 + t579;
t585 = t515 + t520 + (-t213 * t544 + t405 * t546) * t485 + (t213 * t33 + t29 * t405) * mrSges(7,3) + (t213 * t46 + t405 * t45) * mrSges(6,3) - t115 * (mrSges(7,1) * t213 + mrSges(7,2) * t405) - t200 * (mrSges(6,1) * t213 + mrSges(6,2) * t405) + t567 * t496 + (-t213 * t545 + t566 + t580) * t549 + (t548 * t405 - t579) * t497;
t477 = -t375 / 0.2e1;
t275 = pkin(3) * t324 - pkin(9) * t323;
t245 = pkin(2) * t428 + t275;
t283 = t347 * t380 + t326;
t179 = t379 * t245 - t283 * t375;
t438 = t323 * t379;
t403 = t324 * pkin(4) - pkin(10) * t438;
t364 = pkin(2) * t376 + pkin(9);
t469 = -pkin(10) - t364;
t408 = qJD(4) * t469;
t460 = pkin(2) * qJD(3);
t417 = t380 * t460;
t584 = -t375 * t417 + t379 * t408 - t179 - t403;
t183 = t379 * t275 - t279 * t375;
t509 = -pkin(10) - pkin(9);
t414 = qJD(4) * t509;
t583 = t379 * t414 - t183 - t403;
t180 = t375 * t245 + t379 * t283;
t439 = t323 * t375;
t421 = pkin(10) * t439;
t582 = -t375 * t408 - t379 * t417 + t180 - t421;
t184 = t375 * t275 + t379 * t279;
t581 = -t375 * t414 + t184 - t421;
t341 = t374 * t379 + t375 * t378;
t238 = t341 * t323;
t516 = qJD(4) + qJD(5);
t286 = t516 * t341;
t564 = t238 - t286;
t392 = t374 * t375 - t378 * t379;
t239 = t392 * t323;
t285 = t516 * t392;
t563 = t239 - t285;
t578 = t373 * Ifges(4,6) / 0.2e1;
t577 = -t545 * t74 / 0.2e1 - t547 * t73 / 0.2e1 - t544 * t270 / 0.2e1;
t465 = Ifges(5,4) * t296;
t186 = t295 * Ifges(5,2) + t319 * Ifges(5,6) + t465;
t401 = mrSges(5,1) * t375 + mrSges(5,2) * t379;
t576 = t186 * t477 + t251 * t401;
t453 = t373 * Ifges(4,5);
t575 = t357 * mrSges(4,2) + t453 / 0.2e1;
t454 = t324 * Ifges(4,4);
t574 = t578 + t454 / 0.2e1 + t323 * Ifges(4,2) / 0.2e1;
t570 = t270 * t546 + t547 * t74 + t548 * t73;
t336 = t469 * t375;
t372 = t379 * pkin(10);
t337 = t364 * t379 + t372;
t272 = t374 * t336 + t378 * t337;
t530 = -qJD(5) * t272 + t374 * t582 + t378 * t584;
t529 = t336 * t422 - t337 * t423 + t374 * t584 - t378 * t582;
t358 = t509 * t375;
t360 = pkin(9) * t379 + t372;
t299 = t374 * t358 + t378 * t360;
t528 = -qJD(5) * t299 + t374 * t581 + t378 * t583;
t527 = t358 * t422 - t360 * t423 + t374 * t583 - t378 * t581;
t565 = qJ(6) * t564 - qJD(6) * t392;
t370 = pkin(4) * t425;
t562 = -pkin(5) * t564 + t370;
t282 = t347 * t376 - t431;
t307 = pkin(4) * t439;
t561 = t376 * t460 - t282 - t307;
t560 = -t324 * pkin(5) - qJ(6) * t563 - qJD(6) * t341;
t559 = -t287 * t375 + t342 * t424;
t558 = -t51 * t375 + t379 * t50;
t396 = Ifges(5,5) * t379 - Ifges(5,6) * t375;
t463 = Ifges(5,4) * t379;
t398 = -Ifges(5,2) * t375 + t463;
t464 = Ifges(5,4) * t375;
t400 = Ifges(5,1) * t379 - t464;
t294 = Ifges(5,4) * t295;
t187 = Ifges(5,1) * t296 + Ifges(5,5) * t319 + t294;
t432 = t379 * t187;
t486 = t296 / 0.2e1;
t552 = t319 * t396 / 0.2e1 + t400 * t486 + t295 * t398 / 0.2e1 + t432 / 0.2e1 + t576;
t551 = -t357 * mrSges(4,1) - t158 * mrSges(5,1) - t45 * mrSges(6,1) - t29 * mrSges(7,1) + t159 * mrSges(5,2) + t46 * mrSges(6,2) + t33 * mrSges(7,2) + t574;
t543 = mrSges(4,2) * t324;
t534 = t560 + t530;
t533 = t560 + t528;
t532 = t565 + t529;
t531 = t565 + t527;
t258 = t392 * t342;
t221 = t307 + t280;
t525 = -t221 + t562;
t524 = t561 + t562;
t278 = t340 * pkin(3) - t342 * pkin(9) + t368;
t300 = t359 * t376 - t361 * t380;
t201 = t379 * t278 - t300 * t375;
t147 = pkin(4) * t340 - t342 * t372 + t201;
t290 = t379 * t300;
t202 = t375 * t278 + t290;
t435 = t342 * t375;
t167 = -pkin(10) * t435 + t202;
t83 = t374 * t147 + t378 * t167;
t523 = Ifges(5,5) * t196 + Ifges(5,6) * t197;
t522 = t370 + t561;
t521 = t380 * t359 + t361 * t376;
t519 = -t158 * t375 + t159 * t379;
t518 = t51 * mrSges(5,1) - t50 * mrSges(5,2);
t111 = -mrSges(5,1) * t197 + mrSges(5,2) * t196;
t190 = qJD(3) * t280 + t333 * t376 - t380 * t391;
t517 = m(5) * t190 + t111;
t514 = m(4) / 0.2e1;
t513 = t73 / 0.2e1;
t512 = t74 / 0.2e1;
t508 = pkin(1) * mrSges(3,1);
t507 = pkin(1) * mrSges(3,2);
t503 = t196 / 0.2e1;
t502 = t197 / 0.2e1;
t499 = t405 / 0.2e1;
t491 = t270 / 0.2e1;
t488 = -t295 / 0.2e1;
t487 = -t296 / 0.2e1;
t484 = t309 / 0.2e1;
t483 = -t319 / 0.2e1;
t482 = -t323 / 0.2e1;
t476 = t379 / 0.2e1;
t475 = m(4) * t357;
t473 = pkin(2) * t380;
t468 = mrSges(4,3) * t323;
t467 = mrSges(4,3) * t324;
t466 = Ifges(3,4) * t377;
t316 = Ifges(4,4) * t323;
t459 = t295 * Ifges(5,6);
t458 = t296 * Ifges(5,5);
t457 = t319 * Ifges(5,3);
t455 = t324 * Ifges(4,1);
t449 = Ifges(3,5) * qJD(2);
t448 = Ifges(3,6) * qJD(2);
t447 = qJ(6) * t341;
t446 = qJD(2) * mrSges(3,1);
t445 = qJD(2) * mrSges(3,2);
t442 = t190 * t521;
t58 = t378 * t129 - t124;
t430 = -mrSges(4,1) * t373 - mrSges(5,1) * t295 + mrSges(5,2) * t296 + t467;
t416 = Ifges(5,3) * t270 + t523;
t367 = -pkin(4) * t379 - pkin(3);
t412 = t449 / 0.2e1;
t411 = -t448 / 0.2e1;
t21 = -t74 * mrSges(7,1) + t73 * mrSges(7,2);
t57 = -t129 * t374 - t126;
t82 = t378 * t147 - t167 * t374;
t195 = pkin(3) * t288 + pkin(9) * t287 + t418;
t352 = t377 * t415;
t353 = t381 * t415;
t216 = qJD(3) * t521 + t352 * t380 + t353 * t376;
t406 = t379 * t195 - t216 * t375;
t271 = t378 * t336 - t337 * t374;
t297 = t378 * t358 - t360 * t374;
t241 = pkin(4) * t435 - t521;
t402 = mrSges(5,1) * t379 - mrSges(5,2) * t375;
t399 = Ifges(5,1) * t375 + t463;
t397 = Ifges(5,2) * t379 + t464;
t395 = Ifges(5,5) * t375 + Ifges(5,6) * t379;
t136 = mrSges(5,1) * t270 - mrSges(5,3) * t196;
t137 = -mrSges(5,2) * t270 + mrSges(5,3) * t197;
t394 = -t375 * t136 + t379 * t137;
t393 = -t158 * t379 - t159 * t375;
t306 = pkin(5) * t392 + t367;
t39 = t287 * t372 + pkin(4) * t288 + (-t290 + (pkin(10) * t342 - t278) * t375) * qJD(4) + t406;
t75 = t375 * t195 + t379 * t216 + t278 * t424 - t300 * t425;
t49 = -pkin(10) * t559 + t75;
t11 = t147 * t422 - t167 * t423 + t374 * t39 + t378 * t49;
t388 = t393 * mrSges(5,3);
t12 = -qJD(5) * t83 - t374 * t49 + t378 * t39;
t217 = qJD(3) * t300 + t352 * t376 - t380 * t353;
t135 = pkin(4) * t559 + t217;
t107 = -pkin(4) * t197 + t190;
t383 = m(5) * (qJD(4) * t393 + t558);
t101 = t213 * Ifges(7,5) + Ifges(7,6) * t405 + t309 * Ifges(7,3);
t102 = t213 * Ifges(6,5) + Ifges(6,6) * t405 + t309 * Ifges(6,3);
t185 = t457 + t458 + t459;
t247 = t316 + t453 + t455;
t28 = -pkin(5) * t74 + t107;
t86 = t196 * Ifges(5,4) + t197 * Ifges(5,2) + t270 * Ifges(5,6);
t87 = Ifges(5,1) * t196 + Ifges(5,4) * t197 + Ifges(5,5) * t270;
t382 = t552 * qJD(4) + (Ifges(5,5) * t487 - Ifges(4,2) * t482 + Ifges(5,6) * t488 + Ifges(5,3) * t483 + t572 * t485 + t546 * t497 + t544 * t549 + t551 + t578) * t324 + t566 * (-t285 / 0.2e1 + t239 / 0.2e1) + t567 * (-t286 / 0.2e1 + t238 / 0.2e1) + (t396 * t483 + t398 * t488 + t400 * t487 - t575 - t576) * t323 + (-t402 - mrSges(4,1)) * t190 + (t341 * t546 - t392 * t544 + t395) * t491 + (-t285 * t546 - t286 * t544) * t484 + (-t238 * t544 - t239 * t546) * t485 + (-t285 * t547 - t286 * t545) * t499 + (t341 * t547 - t392 * t545) * t512 + (-t238 * t545 - t239 * t547) * t549 + (-t285 * t548 - t286 * t547) * t496 + (t341 * t548 - t392 * t547) * t513 + (-t238 * t547 - t239 * t548) * t497 + t570 * t341 / 0.2e1 + (t158 * t438 + t159 * t439 + t558) * mrSges(5,3) + (-t341 * t9 - t392 * t8 - t45 * t563 + t46 * t564) * mrSges(6,3) + (-t2 * t341 - t29 * t563 + t33 * t564 - t392 * t4) * mrSges(7,3) + (-mrSges(6,1) * t564 + mrSges(6,2) * t563) * t200 + (-mrSges(7,1) * t564 + mrSges(7,2) * t563) * t115 - (Ifges(4,1) * t323 + t101 + t102 + t185 - t454) * t324 / 0.2e1 + t107 * (mrSges(6,1) * t392 + mrSges(6,2) * t341) + t28 * (mrSges(7,1) * t392 + mrSges(7,2) * t341) + (t316 + t432 + t247) * t482 + t279 * t468 + t86 * t476 + t397 * t502 + t399 * t503 - t189 * mrSges(4,2) - Ifges(4,5) * t269 - Ifges(4,6) * t270 + t392 * t577 + t375 * t87 / 0.2e1;
t369 = Ifges(3,4) * t427;
t365 = pkin(4) * t378 + pkin(5);
t356 = mrSges(3,3) * t427 - t445;
t355 = -mrSges(3,3) * t428 + t446;
t354 = t367 - t473;
t335 = t392 * qJ(6);
t322 = Ifges(3,1) * t428 + t369 + t449;
t321 = t448 + (t381 * Ifges(3,2) + t466) * qJD(1);
t304 = -mrSges(4,2) * t373 + t468;
t303 = t306 - t473;
t274 = -mrSges(4,1) * t323 + t543;
t257 = t341 * t342;
t250 = -t335 + t299;
t249 = t297 - t447;
t231 = mrSges(5,1) * t319 - mrSges(5,3) * t296;
t230 = -mrSges(5,2) * t319 + mrSges(5,3) * t295;
t228 = -t335 + t272;
t227 = t271 - t447;
t174 = pkin(5) * t257 + t241;
t173 = mrSges(6,1) * t309 - mrSges(6,3) * t213;
t172 = mrSges(7,1) * t309 - mrSges(7,3) * t213;
t171 = -mrSges(6,2) * t309 + mrSges(6,3) * t405;
t170 = -mrSges(7,2) * t309 + mrSges(7,3) * t405;
t157 = pkin(4) * t296 + pkin(5) * t213;
t123 = -mrSges(6,1) * t405 + mrSges(6,2) * t213;
t122 = -mrSges(7,1) * t405 + mrSges(7,2) * t213;
t100 = t258 * t516 + t341 * t287;
t99 = -t286 * t342 + t287 * t392;
t76 = -qJD(4) * t202 + t406;
t63 = -qJ(6) * t257 + t83;
t56 = pkin(5) * t340 + qJ(6) * t258 + t82;
t55 = -mrSges(6,2) * t270 + mrSges(6,3) * t74;
t54 = -mrSges(7,2) * t270 + mrSges(7,3) * t74;
t53 = mrSges(6,1) * t270 - mrSges(6,3) * t73;
t52 = mrSges(7,1) * t270 - mrSges(7,3) * t73;
t43 = -pkin(5) * t100 + t135;
t35 = t58 - t568;
t34 = t57 - t526;
t22 = -mrSges(6,1) * t74 + mrSges(6,2) * t73;
t10 = qJ(6) * t100 - qJD(6) * t257 + t11;
t5 = pkin(5) * t288 - qJ(6) * t99 + qJD(6) * t258 + t12;
t1 = [-(t247 / 0.2e1 + t316 / 0.2e1 + t455 / 0.2e1 + t388 + t552 + t575) * t287 - (mrSges(4,2) * t368 - mrSges(4,3) * t521) * t269 - t521 * t111 + (mrSges(4,1) * t554 - mrSges(4,3) * t189 + Ifges(4,4) * t269 + t491 * t572 + t544 * t512 + t546 * t513 + t515 + t518) * t340 + (t87 * t476 + t86 * t477 - Ifges(4,1) * t269 + t400 * t503 + t398 * t502 + (mrSges(4,3) + t401) * t190 + (-t375 * t50 - t379 * t51) * mrSges(5,3) + (-t379 * t186 / 0.2e1 + t187 * t477 + t251 * t402 + t397 * t488 + t399 * t487 + t395 * t483 - t519 * mrSges(5,3)) * qJD(4)) * t342 - t570 * t258 / 0.2e1 + t566 * t99 / 0.2e1 + t567 * t100 / 0.2e1 + (t279 * t287 - t280 * t288) * mrSges(4,3) + m(5) * (t158 * t76 + t159 * t75 + t201 * t51 + t202 * t50 + t217 * t251 - t442) + m(4) * (t189 * t300 + t216 * t280 - t217 * t279 - t442) + (t416 + t520 + t523) * t340 / 0.2e1 + (-t551 + (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t405 + t102 / 0.2e1 + t101 / 0.2e1 + t458 / 0.2e1 + t459 / 0.2e1 + t185 / 0.2e1 + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t309 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t213 + t457 / 0.2e1 - t574) * t288 + (t100 * t33 + t2 * t258 - t257 * t4 - t29 * t99) * mrSges(7,3) + (t100 * t46 - t257 * t8 + t258 * t9 - t45 * t99) * mrSges(6,3) + t28 * (mrSges(7,1) * t257 - mrSges(7,2) * t258) + t107 * (mrSges(6,1) * t257 - mrSges(6,2) * t258) + m(7) * (t10 * t33 + t115 * t43 + t174 * t28 + t2 * t56 + t29 * t5 + t4 * t63) + m(6) * (t107 * t241 + t11 * t46 + t12 * t45 + t135 * t200 + t8 * t83 + t82 * t9) + (-Ifges(4,4) * t342 + t368 * mrSges(4,1) + (Ifges(4,2) + Ifges(5,3) / 0.2e1) * t340 - t300 * mrSges(4,3)) * t270 + (t322 / 0.2e1 - pkin(7) * t355 + t412 + (-0.2e1 * t507 + 0.3e1 / 0.2e1 * Ifges(3,4) * t381) * qJD(1)) * t381 * qJD(2) + t56 * t52 + t63 * t54 + t82 * t53 + t83 * t55 + (-pkin(7) * t356 - t321 / 0.2e1 + t411 + (-0.2e1 * t508 - 0.3e1 / 0.2e1 * t466 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t381) * qJD(1) + (t274 + 0.2e1 * t475 + t543) * pkin(2)) * t426 + (-t257 * t544 - t258 * t546 + t342 * t396) * t491 + (t100 * t544 + t546 * t99) * t484 + (t100 * t545 + t547 * t99) * t499 + (-t257 * t545 - t258 * t547) * t512 + (t100 * t547 + t548 * t99) * t496 + (-t257 * t547 - t258 * t548) * t513 + t115 * (-mrSges(7,1) * t100 + mrSges(7,2) * t99) + t43 * t122 + t135 * t123 + t10 * t170 + t11 * t171 + t5 * t172 + t12 * t173 + t174 * t21 + t200 * (-mrSges(6,1) * t100 + mrSges(6,2) * t99) + t201 * t136 + t202 * t137 + t75 * t230 + t76 * t231 + t241 * t22 + t257 * t577 + t216 * t304 + t430 * t217; ((t269 * t380 - t270 * t376) * mrSges(4,3) + (t430 * t376 + (t230 * t379 - t231 * t375 + t304) * t380) * qJD(3)) * pkin(2) + t522 * t123 + t524 * t122 + t517 * (-pkin(3) - t473) + 0.2e1 * ((t189 * t376 - t190 * t380) * t514 + (m(5) * (t251 * t376 + t380 * t519) / 0.2e1 + (-t279 * t376 + t280 * t380) * t514) * qJD(3)) * pkin(2) - m(4) * (-t279 * t282 + t280 * t283) - m(5) * (t158 * t179 + t159 * t180 + t251 * t282) + t382 + t280 * t467 + t388 * qJD(4) + ((-t230 * t375 - t231 * t379) * qJD(4) + t383 + t394) * t364 + ((-t322 / 0.2e1 - t369 / 0.2e1 + t412 + qJD(1) * t507 + (t355 - t446) * pkin(7)) * t381 + (t411 + t321 / 0.2e1 + (t508 + t466 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t381) * qJD(1) + (t356 + t445) * pkin(7) + (-t274 - t475) * pkin(2)) * t377) * qJD(1) + t529 * t171 + t530 * t173 + (t107 * t354 + t200 * t522 + t271 * t9 + t272 * t8 + t45 * t530 + t46 * t529) * m(6) + t532 * t170 + t534 * t172 + (t115 * t524 + t2 * t227 + t228 * t4 + t28 * t303 + t29 * t534 + t33 * t532) * m(7) + t227 * t52 + t228 * t54 - t180 * t230 - t179 * t231 + t271 * t53 + t272 * t55 + t303 * t21 - t283 * t304 - t430 * t282 + t354 * t22; ((-mrSges(5,3) * t158 - pkin(9) * t231) * t379 + (-t159 * mrSges(5,3) + pkin(4) * t123 - pkin(9) * t230) * t375) * qJD(4) + t382 + t525 * t122 + pkin(9) * t383 + t533 * t172 + t527 * t171 + t531 * t170 + t528 * t173 + (-t430 + t467) * t280 + t394 * pkin(9) - t221 * t123 - t184 * t230 - t183 * t231 + t249 * t52 + t250 * t54 + t297 * t53 + t299 * t55 - t279 * t304 + t306 * t21 - m(5) * (t158 * t183 + t159 * t184 + t251 * t280) + t367 * t22 - t517 * pkin(3) + (t115 * t525 + t2 * t249 + t250 * t4 + t28 * t306 + t29 * t533 + t33 * t531) * m(7) + (t107 * t367 + t297 * t9 + t299 * t8 + t527 * t46 + t528 * t45 + (-t221 + t370) * t200) * m(6); (-Ifges(5,2) * t296 + t187 + t294) * t488 + (-t296 * t123 + t378 * t53 + (t54 + t55) * t374 + ((t170 + t171) * t378 + (-t172 - t173) * t374) * qJD(5) + (-t200 * t296 + t374 * t8 + t378 * t9 + t422 * t46 - t423 * t45) * m(6)) * pkin(4) + (t2 * t365 - t115 * t157 - t29 * t34 - t33 * t35 + (-t29 * t423 + t33 * t422 + t374 * t4) * pkin(4)) * m(7) + t518 - m(6) * (t45 * t57 + t46 * t58) + (Ifges(5,5) * t295 - Ifges(5,6) * t296) * t483 + t186 * t486 + (Ifges(5,1) * t295 - t465) * t487 + t416 - t157 * t122 - t35 * t170 - t58 * t171 - t34 * t172 - t57 * t173 - t158 * t230 + t159 * t231 - t251 * (mrSges(5,1) * t296 + mrSges(5,2) * t295) + (t158 * t295 + t159 * t296) * mrSges(5,3) + t365 * t52 + t585; (-t122 * t213 + t52) * pkin(5) + (-(-t29 + t32) * t33 + (-t115 * t213 + t2) * pkin(5)) * m(7) - t32 * t170 - t45 * t171 + t33 * t172 + t46 * t173 + t585; -t405 * t170 + t213 * t172 + 0.2e1 * (t28 / 0.2e1 + t33 * t549 + t29 * t496) * m(7) + t21;];
tauc  = t1(:);
