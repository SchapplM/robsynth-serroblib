% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:50
% EndTime: 2019-03-09 13:22:55
% DurationCPUTime: 35.86s
% Computational Cost: add. (26632->958), mult. (60563->1261), div. (0->0), fcn. (46093->18), ass. (0->436)
t350 = sin(qJ(2));
t355 = cos(qJ(2));
t462 = sin(pkin(11));
t394 = qJD(1) * t462;
t463 = cos(pkin(11));
t395 = qJD(1) * t463;
t287 = -t350 * t394 + t355 * t395;
t288 = -t350 * t395 - t355 * t394;
t430 = qJD(1) * t350;
t414 = pkin(2) * t430;
t190 = -pkin(3) * t288 - pkin(8) * t287 + t414;
t346 = -qJ(3) - pkin(7);
t320 = t346 * t355;
t308 = qJD(1) * t320;
t291 = t462 * t308;
t318 = t346 * t350;
t307 = qJD(1) * t318;
t227 = t307 * t463 + t291;
t349 = sin(qJ(4));
t354 = cos(qJ(4));
t144 = t354 * t190 - t227 * t349;
t406 = t462 * pkin(2);
t326 = t406 + pkin(8);
t486 = pkin(9) + t326;
t399 = qJD(4) * t486;
t457 = t287 * t354;
t606 = pkin(4) * t288 + pkin(9) * t457 - t354 * t399 - t144;
t145 = t349 * t190 + t354 * t227;
t458 = t287 * t349;
t605 = -pkin(9) * t458 + t349 * t399 + t145;
t348 = sin(qJ(5));
t353 = cos(qJ(5));
t306 = t348 * t354 + t349 * t353;
t185 = t306 * t287;
t567 = qJD(4) + qJD(5);
t233 = t567 * t306;
t575 = t185 - t233;
t373 = t348 * t349 - t353 * t354;
t186 = t373 * t287;
t232 = t567 * t373;
t574 = t186 - t232;
t298 = t486 * t349;
t299 = t486 * t354;
t213 = -t348 * t298 + t353 * t299;
t584 = -qJD(5) * t213 + t605 * t348 + t606 * t353;
t422 = qJD(5) * t353;
t423 = qJD(5) * t348;
t583 = -t298 * t422 - t299 * t423 + t606 * t348 - t605 * t353;
t604 = pkin(10) * t575 + t583;
t603 = pkin(5) * t288 - pkin(10) * t574 + t584;
t564 = m(4) + m(5) + m(6) + m(7);
t242 = qJD(2) * t354 + t288 * t349;
t243 = qJD(2) * t349 - t288 * t354;
t170 = t242 * t348 + t243 * t353;
t347 = sin(qJ(6));
t352 = cos(qJ(6));
t391 = t353 * t242 - t243 * t348;
t110 = t170 * t352 + t347 * t391;
t419 = qJD(1) * qJD(2);
t403 = t350 * t419;
t418 = qJDD(1) * t355;
t310 = -t403 + t418;
t311 = qJDD(1) * t350 + t355 * t419;
t229 = t310 * t463 - t462 * t311;
t224 = qJDD(4) - t229;
t221 = qJDD(5) + t224;
t211 = qJDD(6) + t221;
t393 = -t170 * t347 + t352 * t391;
t230 = t310 * t462 + t311 * t463;
t162 = qJD(4) * t242 + qJDD(2) * t349 + t230 * t354;
t163 = -qJD(4) * t243 + qJDD(2) * t354 - t230 * t349;
t78 = qJD(5) * t391 + t162 * t353 + t163 * t348;
t79 = -qJD(5) * t170 - t162 * t348 + t163 * t353;
t30 = qJD(6) * t393 + t347 * t79 + t352 * t78;
t31 = -qJD(6) * t110 - t347 * t78 + t352 * t79;
t416 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t211;
t275 = qJD(4) - t287;
t267 = qJD(5) + t275;
t342 = t355 * pkin(2);
t333 = t342 + pkin(1);
t314 = -qJD(1) * t333 + qJD(3);
t179 = -pkin(3) * t287 + pkin(8) * t288 + t314;
t297 = qJD(2) * pkin(2) + t307;
t398 = t463 * t308;
t219 = t462 * t297 - t398;
t201 = qJD(2) * pkin(8) + t219;
t133 = t179 * t349 + t201 * t354;
t113 = pkin(9) * t242 + t133;
t102 = t348 * t113;
t132 = t354 * t179 - t201 * t349;
t112 = -pkin(9) * t243 + t132;
t95 = pkin(4) * t275 + t112;
t59 = t353 * t95 - t102;
t599 = pkin(10) * t170;
t48 = t59 - t599;
t43 = pkin(5) * t267 + t48;
t592 = pkin(10) * t391;
t104 = t353 * t113;
t60 = t348 * t95 + t104;
t49 = t60 + t592;
t468 = t347 * t49;
t16 = t352 * t43 - t468;
t461 = qJDD(1) * pkin(1);
t270 = -pkin(2) * t310 + qJDD(3) - t461;
t142 = -pkin(3) * t229 - pkin(8) * t230 + t270;
t301 = t311 * pkin(7);
t427 = qJD(3) * t350;
t208 = qJDD(2) * pkin(2) - qJ(3) * t311 - qJD(1) * t427 - t301;
t334 = pkin(7) * t418;
t428 = qJD(2) * t350;
t411 = pkin(7) * t428;
t426 = qJD(3) * t355;
t222 = qJ(3) * t310 + t334 + (-t411 + t426) * qJD(1);
t154 = t462 * t208 + t463 * t222;
t148 = qJDD(2) * pkin(8) + t154;
t58 = -qJD(4) * t133 + t354 * t142 - t148 * t349;
t39 = pkin(4) * t224 - pkin(9) * t162 + t58;
t424 = qJD(4) * t354;
t425 = qJD(4) * t349;
t57 = t349 * t142 + t354 * t148 + t179 * t424 - t201 * t425;
t47 = pkin(9) * t163 + t57;
t13 = -qJD(5) * t60 - t348 * t47 + t353 * t39;
t6 = pkin(5) * t221 - pkin(10) * t78 + t13;
t12 = -t113 * t423 + t348 * t39 + t353 * t47 + t95 * t422;
t7 = pkin(10) * t79 + t12;
t2 = qJD(6) * t16 + t347 * t6 + t352 * t7;
t466 = t352 * t49;
t17 = t347 * t43 + t466;
t3 = -qJD(6) * t17 - t347 * t7 + t352 * t6;
t559 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t372 = t416 + t559;
t415 = Ifges(6,5) * t78 + Ifges(6,6) * t79 + Ifges(6,3) * t221;
t263 = qJD(6) + t267;
t514 = -t263 / 0.2e1;
t528 = -t110 / 0.2e1;
t530 = -t393 / 0.2e1;
t99 = Ifges(7,4) * t393;
t56 = Ifges(7,1) * t110 + Ifges(7,5) * t263 + t99;
t539 = -t56 / 0.2e1;
t478 = Ifges(7,4) * t110;
t55 = Ifges(7,2) * t393 + Ifges(7,6) * t263 + t478;
t541 = -t55 / 0.2e1;
t557 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t218 = t297 * t463 + t291;
t200 = -qJD(2) * pkin(3) - t218;
t164 = -t242 * pkin(4) + t200;
t98 = -pkin(5) * t391 + t164;
t562 = mrSges(7,2) * t98 - mrSges(7,3) * t16;
t563 = -mrSges(7,1) * t98 + t17 * mrSges(7,3);
t602 = t372 + t415 + t557 + (Ifges(7,1) * t528 + Ifges(7,4) * t530 + Ifges(7,5) * t514 + t539 - t562) * t393 - (Ifges(7,4) * t528 + Ifges(7,2) * t530 + Ifges(7,6) * t514 + t541 - t563) * t110;
t345 = qJ(4) + qJ(5);
t339 = cos(t345);
t341 = t354 * pkin(4);
t313 = pkin(5) * t339 + t341;
t309 = pkin(3) + t313;
t340 = qJ(6) + t345;
t329 = sin(t340);
t330 = cos(t340);
t332 = t341 + pkin(3);
t338 = sin(t345);
t386 = -mrSges(5,1) * t354 + mrSges(5,2) * t349;
t601 = -m(5) * pkin(3) - m(6) * t332 - m(7) * t309 - mrSges(6,1) * t339 - mrSges(7,1) * t330 + mrSges(6,2) * t338 + mrSges(7,2) * t329 + t386;
t600 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t496 = pkin(5) * t170;
t319 = -mrSges(3,1) * t355 + mrSges(3,2) * t350;
t343 = qJ(2) + pkin(11);
t336 = sin(t343);
t337 = cos(t343);
t596 = -mrSges(4,1) * t337 - t336 * t600 + t319;
t302 = t350 * t462 - t355 * t463;
t290 = t302 * qJD(2);
t303 = t350 * t463 + t355 * t462;
t405 = t303 * t424;
t370 = -t290 * t349 + t405;
t225 = t307 * t462 - t398;
t572 = -t225 + (t425 - t458) * pkin(4);
t547 = m(6) * pkin(4);
t545 = t30 / 0.2e1;
t544 = t31 / 0.2e1;
t537 = t78 / 0.2e1;
t536 = t79 / 0.2e1;
t526 = t162 / 0.2e1;
t525 = t163 / 0.2e1;
t520 = t211 / 0.2e1;
t519 = t221 / 0.2e1;
t518 = t224 / 0.2e1;
t595 = t310 / 0.2e1;
t498 = pkin(4) * t349;
t594 = m(6) * t498;
t495 = pkin(5) * t338;
t312 = t495 + t498;
t593 = m(7) * t312;
t488 = qJD(2) / 0.2e1;
t212 = -t353 * t298 - t299 * t348;
t174 = -pkin(10) * t306 + t212;
t175 = -pkin(10) * t373 + t213;
t117 = t174 * t352 - t175 * t347;
t591 = qJD(6) * t117 + t603 * t347 + t352 * t604;
t118 = t174 * t347 + t175 * t352;
t590 = -qJD(6) * t118 - t347 * t604 + t603 * t352;
t589 = t242 * Ifges(5,6);
t588 = t275 * Ifges(5,3);
t464 = qJDD(2) / 0.2e1;
t497 = pkin(4) * t353;
t331 = pkin(5) + t497;
t420 = qJD(6) * t352;
t421 = qJD(6) * t347;
t441 = t348 * t352;
t63 = -t112 * t348 - t104;
t50 = t63 - t592;
t64 = t353 * t112 - t102;
t51 = t64 - t599;
t587 = t347 * t51 - t352 * t50 - t331 * t421 + (-t348 * t420 + (-t347 * t353 - t441) * qJD(5)) * pkin(4);
t442 = t347 * t348;
t586 = -t347 * t50 - t352 * t51 + t331 * t420 + (-t348 * t421 + (t352 * t353 - t442) * qJD(5)) * pkin(4);
t585 = t547 + mrSges(5,1);
t582 = Ifges(4,5) * qJD(2);
t581 = Ifges(4,6) * qJD(2);
t197 = t373 * t303;
t226 = -t306 * t347 - t352 * t373;
t125 = qJD(6) * t226 - t232 * t352 - t233 * t347;
t131 = -t185 * t347 - t186 * t352;
t579 = t125 - t131;
t228 = t306 * t352 - t347 * t373;
t126 = -qJD(6) * t228 + t232 * t347 - t233 * t352;
t130 = -t185 * t352 + t186 * t347;
t578 = t126 - t130;
t577 = -pkin(5) * t575 + t572;
t217 = pkin(3) * t302 - pkin(8) * t303 - t333;
t236 = t318 * t462 - t320 * t463;
t158 = t354 * t217 - t236 * t349;
t451 = t303 * t354;
t127 = pkin(4) * t302 - pkin(9) * t451 + t158;
t231 = t354 * t236;
t159 = t349 * t217 + t231;
t452 = t303 * t349;
t137 = -pkin(9) * t452 + t159;
t81 = t348 * t127 + t353 * t137;
t576 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t242 + mrSges(5,2) * t243 - mrSges(4,3) * t288;
t383 = -mrSges(7,1) * t329 - mrSges(7,2) * t330;
t573 = mrSges(6,1) * t338 + mrSges(6,2) * t339 - t383;
t351 = sin(qJ(1));
t444 = t339 * t351;
t356 = cos(qJ(1));
t445 = t338 * t356;
t257 = -t337 * t445 + t444;
t443 = t339 * t356;
t446 = t338 * t351;
t258 = t337 * t443 + t446;
t447 = t337 * t356;
t247 = -t329 * t447 + t330 * t351;
t248 = t329 * t351 + t330 * t447;
t433 = t247 * mrSges(7,1) - t248 * mrSges(7,2);
t571 = -t257 * mrSges(6,1) + t258 * mrSges(6,2) - t433;
t255 = t337 * t446 + t443;
t256 = -t337 * t444 + t445;
t448 = t337 * t351;
t245 = t329 * t448 + t330 * t356;
t246 = t329 * t356 - t330 * t448;
t434 = -t245 * mrSges(7,1) + t246 * mrSges(7,2);
t570 = t255 * mrSges(6,1) - t256 * mrSges(6,2) - t434;
t300 = -pkin(7) * t403 + t334;
t569 = t300 * t355 + t301 * t350;
t568 = -t349 * t58 + t354 * t57;
t566 = t243 * Ifges(5,5) + t170 * Ifges(6,5) + t110 * Ifges(7,5) + Ifges(6,6) * t391 + Ifges(7,6) * t393 + t267 * Ifges(6,3) + t263 * Ifges(7,3) + t588 + t589;
t565 = 0.2e1 * t464;
t561 = -t314 * mrSges(4,2) + t218 * mrSges(4,3);
t558 = t58 * mrSges(5,1) - t57 * mrSges(5,2);
t555 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t593;
t554 = m(3) * pkin(1) + mrSges(2,1) - t596;
t551 = t314 * mrSges(4,1) + t132 * mrSges(5,1) + t59 * mrSges(6,1) + t16 * mrSges(7,1) - t133 * mrSges(5,2) - t60 * mrSges(6,2) - t17 * mrSges(7,2) - t219 * mrSges(4,3);
t549 = Ifges(7,4) * t545 + Ifges(7,2) * t544 + Ifges(7,6) * t520;
t548 = Ifges(7,1) * t545 + Ifges(7,4) * t544 + Ifges(7,5) * t520;
t546 = m(7) * pkin(5);
t543 = Ifges(6,4) * t537 + Ifges(6,2) * t536 + Ifges(6,6) * t519;
t542 = Ifges(6,1) * t537 + Ifges(6,4) * t536 + Ifges(6,5) * t519;
t540 = t55 / 0.2e1;
t538 = t56 / 0.2e1;
t535 = Ifges(5,1) * t526 + Ifges(5,4) * t525 + Ifges(5,5) * t518;
t479 = Ifges(6,4) * t170;
t92 = Ifges(6,2) * t391 + t267 * Ifges(6,6) + t479;
t534 = -t92 / 0.2e1;
t533 = t92 / 0.2e1;
t166 = Ifges(6,4) * t391;
t93 = t170 * Ifges(6,1) + t267 * Ifges(6,5) + t166;
t532 = -t93 / 0.2e1;
t531 = t93 / 0.2e1;
t357 = -pkin(9) - pkin(8);
t529 = t393 / 0.2e1;
t527 = t110 / 0.2e1;
t524 = -t391 / 0.2e1;
t523 = t391 / 0.2e1;
t522 = -t170 / 0.2e1;
t521 = t170 / 0.2e1;
t517 = -t242 / 0.2e1;
t516 = -t243 / 0.2e1;
t515 = t243 / 0.2e1;
t513 = t263 / 0.2e1;
t512 = -t267 / 0.2e1;
t511 = t267 / 0.2e1;
t510 = -t275 / 0.2e1;
t508 = t287 / 0.2e1;
t507 = -t288 / 0.2e1;
t503 = mrSges(6,3) * t59;
t502 = mrSges(6,3) * t60;
t499 = pkin(4) * t243;
t494 = pkin(7) * t355;
t493 = pkin(8) * t336;
t490 = g(3) * t336;
t485 = mrSges(6,3) * t391;
t484 = mrSges(6,3) * t170;
t483 = Ifges(3,4) * t350;
t482 = Ifges(3,4) * t355;
t481 = Ifges(5,4) * t349;
t480 = Ifges(5,4) * t354;
t477 = t132 * mrSges(5,3);
t476 = t133 * mrSges(5,3);
t473 = t243 * Ifges(5,4);
t472 = t288 * Ifges(4,4);
t455 = t290 * t354;
t344 = -pkin(10) + t357;
t450 = t336 * t344;
t449 = t336 * t357;
t156 = t242 * Ifges(5,2) + t275 * Ifges(5,6) + t473;
t440 = t349 * t156;
t439 = t349 * t351;
t438 = t349 * t356;
t437 = t351 * t354;
t237 = Ifges(5,4) * t242;
t157 = t243 * Ifges(5,1) + t275 * Ifges(5,5) + t237;
t436 = t354 * t157;
t435 = t354 * t356;
t429 = qJD(1) * t355;
t417 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t430) * t494;
t413 = pkin(2) * t428;
t408 = Ifges(5,5) * t162 + Ifges(5,6) * t163 + Ifges(5,3) * t224;
t407 = t463 * pkin(2);
t404 = t436 / 0.2e1;
t401 = -t425 / 0.2e1;
t400 = qJD(2) * t346;
t396 = -t229 * mrSges(4,1) + t230 * mrSges(4,2);
t80 = t353 * t127 - t137 * t348;
t283 = t350 * t400 + t426;
t284 = t355 * t400 - t427;
t189 = t283 * t463 + t284 * t462;
t289 = t303 * qJD(2);
t191 = pkin(3) * t289 + pkin(8) * t290 + t413;
t392 = -t189 * t349 + t354 * t191;
t327 = -t407 - pkin(3);
t389 = pkin(3) * t337 + t493;
t388 = mrSges(3,1) * t350 + mrSges(3,2) * t355;
t385 = mrSges(5,1) * t349 + mrSges(5,2) * t354;
t382 = Ifges(5,1) * t354 - t481;
t381 = t355 * Ifges(3,2) + t483;
t380 = -Ifges(5,2) * t349 + t480;
t379 = Ifges(3,5) * t355 - Ifges(3,6) * t350;
t378 = Ifges(5,5) * t354 - Ifges(5,6) * t349;
t61 = pkin(5) * t302 + pkin(10) * t197 + t80;
t196 = t306 * t303;
t67 = -pkin(10) * t196 + t81;
t32 = -t347 * t67 + t352 * t61;
t33 = t347 * t61 + t352 * t67;
t153 = t208 * t463 - t462 * t222;
t188 = t283 * t462 - t463 * t284;
t235 = -t463 * t318 - t320 * t462;
t177 = -mrSges(5,2) * t275 + mrSges(5,3) * t242;
t178 = mrSges(5,1) * t275 - mrSges(5,3) * t243;
t376 = t177 * t354 - t178 * t349;
t138 = -t196 * t352 + t197 * t347;
t139 = -t196 * t347 - t197 * t352;
t375 = t309 * t337 - t450;
t374 = t337 * t332 - t449;
t315 = -t341 + t327;
t187 = pkin(4) * t452 + t235;
t371 = pkin(1) * t388;
t278 = -t337 * t438 + t437;
t276 = t337 * t439 + t435;
t369 = t303 * t425 + t455;
t368 = t200 * t385;
t367 = t350 * (Ifges(3,1) * t355 - t483);
t72 = pkin(9) * t455 + pkin(4) * t289 + (-t231 + (pkin(9) * t303 - t217) * t349) * qJD(4) + t392;
t84 = t354 * t189 + t349 * t191 + t217 * t424 - t236 * t425;
t77 = -pkin(9) * t370 + t84;
t20 = t127 * t422 - t137 * t423 + t348 * t72 + t353 * t77;
t146 = pkin(4) * t370 + t188;
t147 = -qJDD(2) * pkin(3) - t153;
t86 = -t163 * pkin(4) + t147;
t21 = -qJD(5) * t81 - t348 * t77 + t353 * t72;
t335 = Ifges(3,4) * t429;
t317 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t429;
t295 = Ifges(3,1) * t430 + Ifges(3,5) * qJD(2) + t335;
t294 = Ifges(3,6) * qJD(2) + qJD(1) * t381;
t286 = pkin(4) * t441 + t331 * t347;
t285 = -pkin(4) * t442 + t331 * t352;
t279 = t337 * t435 + t439;
t277 = -t337 * t437 + t438;
t271 = Ifges(4,4) * t287;
t253 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t287;
t244 = pkin(5) * t373 + t315;
t209 = -mrSges(4,1) * t287 - mrSges(4,2) * t288;
t205 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t230;
t204 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t229;
t195 = -t288 * Ifges(4,1) + t271 + t582;
t194 = t287 * Ifges(4,2) - t472 + t581;
t150 = mrSges(6,1) * t267 - t484;
t149 = -mrSges(6,2) * t267 + t485;
t141 = t196 * pkin(5) + t187;
t140 = t499 + t496;
t124 = -mrSges(5,2) * t224 + mrSges(5,3) * t163;
t123 = mrSges(5,1) * t224 - mrSges(5,3) * t162;
t114 = -mrSges(6,1) * t391 + mrSges(6,2) * t170;
t97 = t197 * t567 + t306 * t290;
t96 = -t233 * t303 + t290 * t373;
t94 = -mrSges(5,1) * t163 + mrSges(5,2) * t162;
t88 = mrSges(7,1) * t263 - mrSges(7,3) * t110;
t87 = -mrSges(7,2) * t263 + mrSges(7,3) * t393;
t85 = -qJD(4) * t159 + t392;
t82 = t162 * Ifges(5,4) + t163 * Ifges(5,2) + t224 * Ifges(5,6);
t73 = -t97 * pkin(5) + t146;
t66 = -mrSges(6,2) * t221 + mrSges(6,3) * t79;
t65 = mrSges(6,1) * t221 - mrSges(6,3) * t78;
t62 = -mrSges(7,1) * t393 + mrSges(7,2) * t110;
t42 = -qJD(6) * t139 - t347 * t96 + t352 * t97;
t41 = qJD(6) * t138 + t347 * t97 + t352 * t96;
t40 = -t79 * pkin(5) + t86;
t36 = -mrSges(6,1) * t79 + mrSges(6,2) * t78;
t27 = -mrSges(7,2) * t211 + mrSges(7,3) * t31;
t26 = mrSges(7,1) * t211 - mrSges(7,3) * t30;
t19 = t352 * t48 - t468;
t18 = -t347 * t48 - t466;
t15 = pkin(10) * t97 + t20;
t14 = pkin(5) * t289 - pkin(10) * t96 + t21;
t10 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t5 = -qJD(6) * t33 + t14 * t352 - t15 * t347;
t4 = qJD(6) * t32 + t14 * t347 + t15 * t352;
t1 = [(Ifges(7,1) * t41 + Ifges(7,4) * t42) * t527 + (Ifges(7,1) * t139 + Ifges(7,4) * t138) * t545 + (Ifges(6,1) * t96 + Ifges(6,4) * t97) * t521 + (-Ifges(5,1) * t369 - Ifges(5,4) * t370) * t515 + Ifges(3,6) * t355 * t464 + t311 * t482 / 0.2e1 + (t138 * t2 - t139 * t3 - t16 * t41 + t17 * t42) * mrSges(7,3) + (t132 * t369 - t133 * t370 - t451 * t58 - t452 * t57) * mrSges(5,3) + (t551 - Ifges(4,6) * t488 + Ifges(6,3) * t511 + Ifges(7,3) * t513 - Ifges(4,4) * t507 - Ifges(4,2) * t508 + t589 / 0.2e1 + t588 / 0.2e1 + Ifges(5,5) * t515 + Ifges(6,5) * t521 + Ifges(6,6) * t523 + Ifges(7,5) * t527 + Ifges(7,6) * t529 - t194 / 0.2e1 + t566 / 0.2e1) * t289 + (Ifges(7,4) * t41 + Ifges(7,2) * t42) * t529 + (Ifges(7,4) * t139 + Ifges(7,2) * t138) * t544 + t275 * (-Ifges(5,5) * t369 - Ifges(5,6) * t370) / 0.2e1 + m(5) * (t132 * t85 + t133 * t84 + t158 * t58 + t159 * t57) + m(4) * (t154 * t236 + t189 * t219 - t270 * t333 + t314 * t413) + (-m(4) * t153 + m(5) * t147 - t205 + t94) * t235 - qJDD(2) * mrSges(3,2) * t494 + (Ifges(7,5) * t41 + Ifges(7,6) * t42) * t513 + (Ifges(7,5) * t139 + Ifges(7,6) * t138) * t520 + t242 * (-Ifges(5,4) * t369 - Ifges(5,2) * t370) / 0.2e1 + m(7) * (t141 * t40 + t16 * t5 + t17 * t4 + t2 * t33 + t3 * t32 + t73 * t98) + m(6) * (t12 * t81 + t13 * t80 + t146 * t164 + t187 * t86 + t20 * t60 + t21 * t59) + (Ifges(6,4) * t96 + Ifges(6,2) * t97) * t523 + t200 * (mrSges(5,1) * t370 - mrSges(5,2) * t369) + (Ifges(6,5) * t96 + Ifges(6,6) * t97) * t511 + t355 * t295 * t488 + t381 * t595 + Ifges(2,3) * qJDD(1) - t319 * t461 - t82 * t452 / 0.2e1 + t355 * (Ifges(3,4) * t311 + Ifges(3,2) * t310 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t294 * t428 / 0.2e1 - t371 * t419 + t98 * (-mrSges(7,1) * t42 + mrSges(7,2) * t41) - t317 * t411 + t5 * t88 + t4 * t87 + t80 * t65 + t81 * t66 + t73 * t62 - t156 * t405 / 0.2e1 + t32 * t26 + t33 * t27 - t333 * t396 + (-t12 * t196 + t13 * t197 - t59 * t96 + t60 * t97) * mrSges(6,3) + (-Ifges(6,4) * t197 - Ifges(6,2) * t196) * t536 + (-Ifges(6,5) * t197 - Ifges(6,6) * t196) * t519 + t86 * (mrSges(6,1) * t196 - mrSges(6,2) * t197) + (-Ifges(6,1) * t197 - Ifges(6,4) * t196) * t537 + (t367 + t355 * (-Ifges(3,2) * t350 + t482)) * t419 / 0.2e1 + (t416 + t415 + t408) * t302 / 0.2e1 + (-t439 * t547 - t279 * mrSges(5,1) - t258 * mrSges(6,1) - t248 * mrSges(7,1) - t278 * mrSges(5,2) - t257 * mrSges(6,2) - t247 * mrSges(7,2) - t564 * (t356 * t333 - t351 * t346) + t555 * t351 + (-m(5) * t389 - m(6) * t374 - m(7) * t375 - t554) * t356) * g(2) + (Ifges(3,1) * t311 + Ifges(3,4) * t595 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t311) + t565 * Ifges(3,5)) * t350 + (-t277 * mrSges(5,1) - t256 * mrSges(6,1) - t246 * mrSges(7,1) - t276 * mrSges(5,2) - t255 * mrSges(6,2) - t245 * mrSges(7,2) + (t346 * t564 + t555 - t594) * t356 + (m(4) * t333 - m(6) * (-t333 - t374) - m(7) * (-t333 - t375) - m(5) * (-t333 - t389) + t554) * t351) * g(1) + (t379 * t488 - t417) * qJD(2) - pkin(1) * (-mrSges(3,1) * t310 + mrSges(3,2) * t311) + t40 * (-mrSges(7,1) * t138 + mrSges(7,2) * t139) + t141 * t10 + t146 * t114 + t20 * t149 + t21 * t150 + t158 * t123 + t159 * t124 + t209 * t413 + t164 * (-mrSges(6,1) * t97 + mrSges(6,2) * t96) + t84 * t177 + t85 * t178 + t187 * t36 + t96 * t531 + t97 * t533 + t451 * t535 + t41 * t538 + t42 * t540 - t197 * t542 - t196 * t543 + t139 * t548 + t138 * t549 - (t404 + Ifges(4,5) * t488 + Ifges(4,1) * t507 + Ifges(4,4) * t508 - t440 / 0.2e1 + t195 / 0.2e1 - t561) * t290 + t236 * t204 + (t270 * mrSges(4,2) - t153 * mrSges(4,3) + Ifges(4,1) * t230 + Ifges(4,4) * t229 + Ifges(4,5) * t565 + t147 * t385 + t157 * t401 + t378 * t518 + t380 * t525 + t382 * t526) * t303 + (t270 * mrSges(4,1) - t154 * mrSges(4,3) - Ifges(4,4) * t230 + Ifges(5,5) * t526 + Ifges(6,5) * t537 + Ifges(7,5) * t545 - Ifges(4,2) * t229 - Ifges(4,6) * t565 + Ifges(5,6) * t525 + Ifges(6,6) * t536 + Ifges(7,6) * t544 + Ifges(5,3) * t518 + Ifges(6,3) * t519 + Ifges(7,3) * t520 + t557 + t558 + t559) * t302 + t189 * t253 + (t310 * t494 + t569) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t569) + (-m(4) * t218 + m(5) * t200 + t576) * t188; ((t153 * t463 + t154 * t462) * pkin(2) + t218 * t225 - t219 * t227 - t314 * t414) * m(4) + (t417 + (t371 - t367 / 0.2e1) * qJD(1)) * qJD(1) - (Ifges(6,4) * t521 + Ifges(6,2) * t523 + Ifges(6,6) * t511 + t502 + t533) * t233 - (Ifges(6,1) * t521 + Ifges(6,4) * t523 + Ifges(6,5) * t511 - t503 + t531) * t232 + (Ifges(7,1) * t131 + Ifges(7,4) * t130) * t528 + (-t130 * t17 + t131 * t16 + t2 * t226 - t228 * t3) * mrSges(7,3) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t404 + t368) * qJD(4) + (-t132 * t144 - t133 * t145 + t147 * t327 - t200 * t225) * m(5) + (t294 / 0.2e1 + pkin(7) * t317) * t430 + (-t12 * t373 - t13 * t306 + t185 * t60 - t186 * t59) * mrSges(6,3) + t86 * (mrSges(6,1) * t373 + mrSges(6,2) * t306) + t156 * t401 + t204 * t406 + t205 * t407 + t194 * t507 + t440 * t508 + t147 * t386 + (Ifges(6,5) * t306 - Ifges(6,6) * t373) * t519 + (Ifges(6,4) * t306 - Ifges(6,2) * t373) * t536 + (Ifges(6,1) * t306 - Ifges(6,4) * t373) * t537 - t373 * t543 + (Ifges(7,4) * t131 + Ifges(7,2) * t130) * t530 + t117 * t26 + t118 * t27 - t425 * t476 - t424 * t477 + t354 * t82 / 0.2e1 - t379 * t419 / 0.2e1 - t209 * t414 + (Ifges(7,5) * t131 + Ifges(7,6) * t130) * t514 + (-Ifges(6,4) * t186 - Ifges(6,2) * t185) * t524 + (-Ifges(6,1) * t186 - Ifges(6,4) * t185) * t522 + (-Ifges(6,5) * t186 - Ifges(6,6) * t185) * t512 - (-Ifges(3,2) * t430 + t295 + t335) * t429 / 0.2e1 + (t242 * t380 + t243 * t382 + t275 * t378) * qJD(4) / 0.2e1 - (Ifges(4,2) * t288 + t195 + t271 + t436) * t287 / 0.2e1 + t590 * t88 + t591 * t87 + (t117 * t3 + t118 * t2 + t16 * t590 + t17 * t591 + t244 * t40 + t577 * t98) * m(7) + t327 * t94 + t315 * t36 + Ifges(3,6) * t310 + Ifges(3,5) * t311 + t583 * t149 + t584 * t150 + (t12 * t213 + t13 * t212 + t164 * t572 + t315 * t86 + t583 * t60 + t584 * t59) * m(6) + (t551 - Ifges(6,3) * t512 - Ifges(5,3) * t510 - Ifges(7,3) * t514 - Ifges(5,5) * t516 - Ifges(5,6) * t517 - Ifges(6,5) * t522 - Ifges(6,6) * t524 - Ifges(7,5) * t528 - Ifges(7,6) * t530 - t581 / 0.2e1) * t288 + (t378 * t510 - t368 - t582 / 0.2e1 + t382 * t516 + t380 * t517 + t561) * t287 - t300 * mrSges(3,2) - t301 * mrSges(3,1) + (Ifges(5,5) * t349 + Ifges(5,6) * t354) * t518 + (Ifges(7,5) * t228 + Ifges(7,6) * t226) * t520 - t98 * (-mrSges(7,1) * t130 + mrSges(7,2) * t131) + t153 * mrSges(4,1) - t154 * mrSges(4,2) - t145 * t177 - t144 * t178 + (Ifges(5,2) * t354 + t481) * t525 + (Ifges(5,1) * t349 + t480) * t526 - t186 * t532 - t185 * t534 + t349 * t535 + t131 * t539 + t130 * t541 + t306 * t542 + (Ifges(7,4) * t228 + Ifges(7,2) * t226) * t544 + (Ifges(7,1) * t228 + Ifges(7,4) * t226) * t545 + t228 * t548 + t226 * t549 + t212 * t65 + t213 * t66 + (Ifges(7,1) * t527 + Ifges(7,4) * t529 + Ifges(7,5) * t513 + t538 + t562) * t125 + (Ifges(7,4) * t527 + Ifges(7,2) * t529 + Ifges(7,6) * t513 + t540 + t563) * t126 + t40 * (-mrSges(7,1) * t226 + mrSges(7,2) * t228) + Ifges(4,6) * t229 + Ifges(4,5) * t230 + (-m(4) * t342 - m(5) * (t342 + t493) - m(7) * (t342 - t450) - m(6) * (t342 - t449) + t601 * t337 + t596) * g(3) + (g(1) * t356 + g(2) * t351) * (t388 + t564 * pkin(2) * t350 + (-m(5) * pkin(8) + m(6) * t357 + m(7) * t344 - t600) * t337 + (mrSges(4,1) - t601) * t336) + (Ifges(4,1) * t287 + t472 + t566) * t288 / 0.2e1 + t244 * t10 - t227 * t253 + (m(5) * ((-t132 * t354 - t133 * t349) * qJD(4) + t568) + t354 * t124 - t349 * t123 - t178 * t424 - t177 * t425) * t326 + (t132 * t457 + t133 * t458 + t568) * mrSges(5,3) + t572 * t114 + (-mrSges(6,1) * t575 + mrSges(6,2) * t574) * t164 - t576 * t225 + t577 * t62; (-t253 - t376) * t287 + t376 * qJD(4) + t575 * t150 + t574 * t149 + t396 + (t114 + t62 + t576) * t288 + t354 * t123 + t349 * t124 + t578 * t88 + t579 * t87 - t373 * t65 + t306 * t66 + t226 * t26 + t228 * t27 + (-g(1) * t351 + g(2) * t356) * t564 + (t16 * t578 + t17 * t579 + t2 * t228 + t226 * t3 + t288 * t98) * m(7) + (t12 * t306 - t13 * t373 + t164 * t288 + t574 * t60 + t575 * t59) * m(6) + (t200 * t288 + t349 * t57 + t354 * t58 + t275 * (-t132 * t349 + t133 * t354)) * m(5) + (-t218 * t288 - t219 * t287 + t270) * m(4); -(mrSges(6,1) * t164 + Ifges(6,4) * t522 + Ifges(6,2) * t524 + Ifges(6,6) * t512 - t502 + t534) * t170 + (-Ifges(5,2) * t243 + t157 + t237) * t517 + (t149 * t422 - t150 * t423 + t348 * t66) * pkin(4) + t243 * t476 + t242 * t477 + t65 * t497 + (Ifges(5,5) * t242 - Ifges(5,6) * t243) * t510 + t408 - t114 * t499 - m(6) * (t164 * t499 + t59 * t63 + t60 * t64) + t602 + t558 + (-mrSges(6,2) * t164 + Ifges(6,1) * t522 + Ifges(6,4) * t524 + Ifges(6,5) * t512 + t503 + t532) * t391 + (t385 + t573 + t593 + t594) * t490 + (-m(7) * (-t312 * t448 - t313 * t356) - mrSges(5,2) * t277 + t585 * t276 + t570) * g(2) + (-m(7) * (-t312 * t447 + t313 * t351) + mrSges(5,2) * t279 - t585 * t278 + t571) * g(1) + t586 * t87 + t587 * t88 + (-t140 * t98 + t587 * t16 + t586 * t17 + t2 * t286 + t285 * t3) * m(7) + t156 * t515 + (Ifges(5,1) * t242 - t473) * t516 - t140 * t62 - t64 * t149 - t63 * t150 - t132 * t177 + t133 * t178 + (t12 * t348 + t13 * t353 + (-t348 * t59 + t353 * t60) * qJD(5)) * t547 - t200 * (mrSges(5,1) * t243 + mrSges(5,2) * t242) + t285 * t26 + t286 * t27; (Ifges(6,5) * t391 - Ifges(6,6) * t170) * t512 - t62 * t496 - m(7) * (t16 * t18 + t17 * t19 + t496 * t98) - t18 * t88 - t19 * t87 + t92 * t521 + (Ifges(6,1) * t391 - t479) * t522 - t164 * (mrSges(6,1) * t170 + mrSges(6,2) * t391) + (t2 * t347 + t3 * t352 + (-t16 * t347 + t17 * t352) * qJD(6)) * t546 + (t484 + t150) * t60 + (t485 - t149) * t59 + (-Ifges(6,2) * t170 + t166 + t93) * t524 + (m(7) * t495 + t573) * t490 + (t255 * t546 + t570) * g(2) + (-t257 * t546 + t571) * g(1) + (t26 * t352 + t27 * t347 + t420 * t87 - t421 * t88) * pkin(5) + t602; -t98 * (mrSges(7,1) * t110 + mrSges(7,2) * t393) + (Ifges(7,1) * t393 - t478) * t528 + t55 * t527 + (Ifges(7,5) * t393 - Ifges(7,6) * t110) * t514 - t16 * t87 + t17 * t88 - g(1) * t433 - g(2) * t434 - t383 * t490 + (t110 * t17 + t16 * t393) * mrSges(7,3) + t372 + (-Ifges(7,2) * t110 + t56 + t99) * t530;];
tau  = t1;
