% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:16
% EndTime: 2019-03-09 07:13:12
% DurationCPUTime: 33.04s
% Computational Cost: add. (25473->886), mult. (60696->1161), div. (0->0), fcn. (47372->18), ass. (0->392)
t328 = sin(pkin(11));
t329 = cos(pkin(11));
t334 = sin(qJ(3));
t339 = cos(qJ(3));
t557 = -t328 * t334 + t339 * t329;
t274 = t557 * qJD(1);
t284 = t328 * t339 + t329 * t334;
t275 = t284 * qJD(1);
t208 = pkin(3) * t275 - pkin(8) * t274;
t458 = pkin(7) + qJ(2);
t299 = t458 * t328;
t285 = qJD(1) * t299;
t300 = t458 * t329;
t286 = qJD(1) * t300;
t215 = -t285 * t339 - t334 * t286;
t333 = sin(qJ(4));
t338 = cos(qJ(4));
t151 = t338 * t208 - t215 * t333;
t341 = -pkin(9) - pkin(8);
t387 = qJD(4) * t341;
t434 = t274 * t338;
t577 = -pkin(4) * t275 + pkin(9) * t434 + t338 * t387 - t151;
t152 = t333 * t208 + t338 * t215;
t435 = t274 * t333;
t576 = -pkin(9) * t435 - t333 * t387 + t152;
t332 = sin(qJ(5));
t337 = cos(qJ(5));
t288 = t332 * t338 + t333 * t337;
t182 = t288 * t274;
t523 = qJD(4) + qJD(5);
t221 = t523 * t288;
t532 = t182 - t221;
t355 = t332 * t333 - t337 * t338;
t183 = t355 * t274;
t220 = t523 * t355;
t531 = t183 - t220;
t302 = t341 * t333;
t303 = t341 * t338;
t227 = t332 * t302 - t337 * t303;
t540 = -qJD(5) * t227 + t576 * t332 + t337 * t577;
t400 = qJD(5) * t337;
t401 = qJD(5) * t332;
t539 = t302 * t400 + t303 * t401 + t332 * t577 - t576 * t337;
t396 = qJD(1) * qJD(2);
t307 = qJ(2) * qJDD(1) + t396;
t575 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t574 = t532 * pkin(10) + t539;
t573 = -pkin(5) * t275 - t531 * pkin(10) + t540;
t406 = t328 ^ 2 + t329 ^ 2;
t206 = -qJD(3) * pkin(3) - t215;
t232 = qJD(3) * t338 - t275 * t333;
t233 = qJD(3) * t333 + t275 * t338;
t456 = mrSges(4,3) * t275;
t530 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t232 - mrSges(5,2) * t233 - t456;
t572 = -m(5) * t206 + t530;
t311 = pkin(2) * t329 + pkin(1);
t296 = -qJD(1) * t311 + qJD(2);
t538 = Ifges(4,5) * qJD(3);
t571 = -t296 * mrSges(4,2) + t215 * mrSges(4,3) - t538 / 0.2e1;
t277 = t284 * qJD(3);
t178 = -pkin(3) * t274 - pkin(8) * t275 + t296;
t216 = -t334 * t285 + t339 * t286;
t207 = qJD(3) * pkin(8) + t216;
t131 = t338 * t178 - t207 * t333;
t132 = t178 * t333 + t207 * t338;
t336 = cos(qJ(6));
t263 = qJD(4) - t274;
t255 = qJD(5) + t263;
t170 = t232 * t332 + t233 * t337;
t562 = pkin(10) * t170;
t113 = pkin(9) * t232 + t132;
t106 = t332 * t113;
t112 = -pkin(9) * t233 + t131;
t95 = pkin(4) * t263 + t112;
t59 = t337 * t95 - t106;
t46 = t59 - t562;
t44 = pkin(5) * t255 + t46;
t331 = sin(qJ(6));
t373 = t337 * t232 - t233 * t332;
t546 = pkin(10) * t373;
t108 = t337 * t113;
t60 = t332 * t95 + t108;
t47 = t60 + t546;
t441 = t331 * t47;
t16 = t336 * t44 - t441;
t439 = t336 * t47;
t17 = t331 * t44 + t439;
t537 = Ifges(4,6) * qJD(3);
t570 = -t296 * mrSges(4,1) - t131 * mrSges(5,1) - t59 * mrSges(6,1) - t16 * mrSges(7,1) + t132 * mrSges(5,2) + t60 * mrSges(6,2) + t17 * mrSges(7,2) + t537 / 0.2e1;
t104 = t170 * t336 + t331 * t373;
t162 = -pkin(4) * t232 + t206;
t105 = -pkin(5) * t373 + t162;
t393 = qJDD(1) * t329;
t394 = qJDD(1) * t328;
t214 = -qJD(1) * t277 - t334 * t394 + t339 * t393;
t205 = qJDD(4) - t214;
t202 = qJDD(5) + t205;
t191 = qJDD(6) + t202;
t375 = -t170 * t331 + t336 * t373;
t276 = t557 * qJD(3);
t213 = qJD(1) * t276 + qJDD(1) * t284;
t160 = qJD(4) * t232 + qJDD(3) * t333 + t213 * t338;
t161 = -qJD(4) * t233 + qJDD(3) * t338 - t213 * t333;
t76 = qJD(5) * t373 + t160 * t337 + t161 * t332;
t77 = -qJD(5) * t170 - t160 * t332 + t161 * t337;
t30 = qJD(6) * t375 + t331 * t77 + t336 * t76;
t31 = -qJD(6) * t104 - t331 * t76 + t336 * t77;
t392 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t191;
t295 = -qJDD(1) * t311 + qJDD(2);
t144 = -pkin(3) * t214 - pkin(8) * t213 + t295;
t376 = pkin(7) * qJDD(1) + t307;
t257 = t376 * t328;
t258 = t376 * t329;
t404 = qJD(3) * t339;
t405 = qJD(3) * t334;
t149 = -t334 * t257 + t339 * t258 - t285 * t404 - t286 * t405;
t145 = qJDD(3) * pkin(8) + t149;
t57 = -t132 * qJD(4) + t338 * t144 - t145 * t333;
t39 = pkin(4) * t205 - pkin(9) * t160 + t57;
t402 = qJD(4) * t338;
t403 = qJD(4) * t333;
t56 = t333 * t144 + t338 * t145 + t178 * t402 - t207 * t403;
t45 = pkin(9) * t161 + t56;
t13 = -qJD(5) * t60 - t332 * t45 + t337 * t39;
t6 = pkin(5) * t202 - pkin(10) * t76 + t13;
t12 = -t113 * t401 + t332 * t39 + t337 * t45 + t95 * t400;
t7 = pkin(10) * t77 + t12;
t2 = qJD(6) * t16 + t331 * t6 + t336 * t7;
t3 = -qJD(6) * t17 - t331 * t7 + t336 * t6;
t552 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t354 = t392 + t552;
t391 = Ifges(6,5) * t76 + Ifges(6,6) * t77 + Ifges(6,3) * t202;
t469 = mrSges(7,3) * t17;
t470 = mrSges(7,3) * t16;
t251 = qJD(6) + t255;
t484 = -t251 / 0.2e1;
t498 = -t104 / 0.2e1;
t500 = -t375 / 0.2e1;
t98 = Ifges(7,4) * t375;
t55 = Ifges(7,1) * t104 + Ifges(7,5) * t251 + t98;
t509 = -t55 / 0.2e1;
t449 = Ifges(7,4) * t104;
t54 = Ifges(7,2) * t375 + Ifges(7,6) * t251 + t449;
t511 = -t54 / 0.2e1;
t550 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t569 = t354 + t391 + t550 + (-mrSges(7,2) * t105 + Ifges(7,1) * t498 + Ifges(7,4) * t500 + Ifges(7,5) * t484 + t470 + t509) * t375 - (mrSges(7,1) * t105 + Ifges(7,4) * t498 + Ifges(7,2) * t500 + Ifges(7,6) * t484 - t469 + t511) * t104;
t327 = qJ(4) + qJ(5);
t321 = cos(t327);
t465 = pkin(4) * t338;
t298 = pkin(5) * t321 + t465;
t292 = pkin(3) + t298;
t322 = qJ(6) + t327;
t313 = sin(t322);
t314 = cos(t322);
t316 = pkin(3) + t465;
t320 = sin(t327);
t366 = -mrSges(5,1) * t338 + mrSges(5,2) * t333;
t568 = -m(5) * pkin(3) - m(6) * t316 - m(7) * t292 - mrSges(6,1) * t321 - mrSges(7,1) * t314 + mrSges(6,2) * t320 + mrSges(7,2) * t313 + t366;
t326 = -pkin(10) + t341;
t567 = -m(5) * pkin(8) + m(6) * t341 + m(7) * t326 - t575;
t517 = m(6) * pkin(4);
t467 = pkin(4) * t333;
t564 = m(6) * t467;
t463 = pkin(5) * t320;
t297 = t463 + t467;
t563 = m(7) * t297;
t464 = pkin(5) * t170;
t561 = Ifges(5,6) * t232;
t560 = Ifges(5,3) * t263;
t384 = t284 * t402;
t353 = t276 * t333 + t384;
t527 = -t216 + (t403 - t435) * pkin(4);
t556 = Ifges(5,5) * t233 + Ifges(6,5) * t170 + Ifges(7,5) * t104 + Ifges(6,6) * t373 + Ifges(7,6) * t375 + Ifges(6,3) * t255 + Ifges(7,3) * t251 + t560 + t561;
t555 = -m(6) - m(5) - m(7) - m(4);
t551 = t57 * mrSges(5,1) - t56 * mrSges(5,2);
t382 = m(3) * qJ(2) + mrSges(3,3);
t549 = mrSges(2,2) - mrSges(4,3) - t382 - t563;
t325 = pkin(11) + qJ(3);
t318 = sin(t325);
t319 = cos(t325);
t368 = mrSges(4,1) * t319 - mrSges(4,2) * t318;
t369 = -mrSges(3,1) * t329 + mrSges(3,2) * t328;
t548 = m(3) * pkin(1) + t318 * t575 + mrSges(2,1) + t368 - t369;
t515 = t30 / 0.2e1;
t514 = t31 / 0.2e1;
t507 = t76 / 0.2e1;
t506 = t77 / 0.2e1;
t496 = t160 / 0.2e1;
t495 = t161 / 0.2e1;
t490 = t191 / 0.2e1;
t489 = t202 / 0.2e1;
t488 = t205 / 0.2e1;
t226 = t337 * t302 + t303 * t332;
t186 = -pkin(10) * t288 + t226;
t187 = -pkin(10) * t355 + t227;
t136 = t186 * t336 - t187 * t331;
t545 = qJD(6) * t136 + t573 * t331 + t336 * t574;
t137 = t186 * t331 + t187 * t336;
t544 = -qJD(6) * t137 - t331 * t574 + t573 * t336;
t466 = pkin(4) * t337;
t315 = pkin(5) + t466;
t398 = qJD(6) * t336;
t399 = qJD(6) * t331;
t418 = t332 * t336;
t65 = -t112 * t332 - t108;
t48 = t65 - t546;
t66 = t337 * t112 - t106;
t49 = t66 - t562;
t543 = t331 * t49 - t336 * t48 - t315 * t399 + (-t332 * t398 + (-t331 * t337 - t418) * qJD(5)) * pkin(4);
t419 = t331 * t332;
t542 = -t331 * t48 - t336 * t49 + t315 * t398 + (-t332 * t399 + (t336 * t337 - t419) * qJD(5)) * pkin(4);
t541 = t517 + mrSges(5,1);
t197 = t355 * t284;
t217 = -t288 * t331 - t336 * t355;
t121 = qJD(6) * t217 - t220 * t336 - t221 * t331;
t127 = -t182 * t331 - t183 * t336;
t535 = t121 - t127;
t218 = t288 * t336 - t331 * t355;
t122 = -qJD(6) * t218 + t220 * t331 - t221 * t336;
t126 = -t182 * t336 + t183 * t331;
t534 = t122 - t126;
t533 = -pkin(5) * t532 + t527;
t212 = -pkin(3) * t557 - pkin(8) * t284 - t311;
t224 = -t299 * t334 + t300 * t339;
t158 = t338 * t212 - t224 * t333;
t430 = t284 * t338;
t120 = -pkin(4) * t557 - pkin(9) * t430 + t158;
t219 = t338 * t224;
t159 = t333 * t212 + t219;
t431 = t284 * t333;
t133 = -pkin(9) * t431 + t159;
t81 = t332 * t120 + t337 * t133;
t529 = -t339 * t299 - t300 * t334;
t363 = -mrSges(7,1) * t313 - mrSges(7,2) * t314;
t528 = mrSges(6,1) * t320 + mrSges(6,2) * t321 - t363;
t335 = sin(qJ(1));
t424 = t321 * t335;
t340 = cos(qJ(1));
t425 = t320 * t340;
t249 = -t319 * t425 + t424;
t423 = t321 * t340;
t426 = t320 * t335;
t250 = t319 * t423 + t426;
t427 = t319 * t340;
t236 = -t313 * t427 + t314 * t335;
t237 = t313 * t335 + t314 * t427;
t409 = t236 * mrSges(7,1) - t237 * mrSges(7,2);
t526 = -t249 * mrSges(6,1) + t250 * mrSges(6,2) - t409;
t247 = t319 * t426 + t423;
t248 = -t319 * t424 + t425;
t428 = t319 * t335;
t234 = t313 * t428 + t314 * t340;
t235 = t313 * t340 - t314 * t428;
t410 = -t234 * mrSges(7,1) + t235 * mrSges(7,2);
t525 = t247 * mrSges(6,1) - t248 * mrSges(6,2) - t410;
t524 = -t333 * t57 + t338 * t56;
t519 = Ifges(7,4) * t515 + Ifges(7,2) * t514 + Ifges(7,6) * t490;
t518 = Ifges(7,1) * t515 + Ifges(7,4) * t514 + Ifges(7,5) * t490;
t516 = m(7) * pkin(5);
t513 = Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * t489;
t512 = Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * t489;
t510 = t54 / 0.2e1;
t508 = t55 / 0.2e1;
t505 = Ifges(5,1) * t496 + Ifges(5,4) * t495 + Ifges(5,5) * t488;
t450 = Ifges(6,4) * t170;
t92 = Ifges(6,2) * t373 + Ifges(6,6) * t255 + t450;
t504 = -t92 / 0.2e1;
t503 = t92 / 0.2e1;
t166 = Ifges(6,4) * t373;
t93 = Ifges(6,1) * t170 + Ifges(6,5) * t255 + t166;
t502 = -t93 / 0.2e1;
t501 = t93 / 0.2e1;
t499 = t375 / 0.2e1;
t497 = t104 / 0.2e1;
t494 = -t373 / 0.2e1;
t493 = t373 / 0.2e1;
t492 = -t170 / 0.2e1;
t491 = t170 / 0.2e1;
t487 = -t232 / 0.2e1;
t486 = -t233 / 0.2e1;
t485 = t233 / 0.2e1;
t483 = t251 / 0.2e1;
t482 = -t255 / 0.2e1;
t481 = t255 / 0.2e1;
t480 = -t263 / 0.2e1;
t478 = t274 / 0.2e1;
t476 = t275 / 0.2e1;
t472 = mrSges(6,3) * t59;
t471 = mrSges(6,3) * t60;
t468 = pkin(4) * t233;
t460 = g(3) * t318;
t455 = mrSges(6,3) * t373;
t454 = mrSges(6,3) * t170;
t453 = Ifges(4,4) * t275;
t452 = Ifges(5,4) * t333;
t451 = Ifges(5,4) * t338;
t448 = t131 * mrSges(5,3);
t447 = t132 * mrSges(5,3);
t445 = t233 * Ifges(5,4);
t432 = t276 * t338;
t154 = t232 * Ifges(5,2) + t263 * Ifges(5,6) + t445;
t417 = t333 * t154;
t416 = t333 * t335;
t415 = t333 * t340;
t414 = t335 * t338;
t225 = Ifges(5,4) * t232;
t155 = t233 * Ifges(5,1) + t263 * Ifges(5,5) + t225;
t413 = t338 * t155;
t412 = t338 * t340;
t388 = Ifges(5,5) * t160 + Ifges(5,6) * t161 + Ifges(5,3) * t205;
t383 = t413 / 0.2e1;
t379 = -t403 / 0.2e1;
t378 = -t214 * mrSges(4,1) + t213 * mrSges(4,2);
t80 = t337 * t120 - t133 * t332;
t179 = qJD(2) * t557 + qJD(3) * t529;
t209 = pkin(3) * t277 - pkin(8) * t276;
t374 = -t179 * t333 + t338 * t209;
t371 = pkin(3) * t319 + pkin(8) * t318;
t181 = pkin(4) * t431 - t529;
t370 = -mrSges(3,1) * t393 + mrSges(3,2) * t394;
t365 = mrSges(5,1) * t333 + mrSges(5,2) * t338;
t362 = Ifges(5,1) * t338 - t452;
t361 = -Ifges(5,2) * t333 + t451;
t360 = Ifges(5,5) * t338 - Ifges(5,6) * t333;
t61 = -pkin(5) * t557 + pkin(10) * t197 + t80;
t196 = t288 * t284;
t67 = -pkin(10) * t196 + t81;
t32 = -t331 * t67 + t336 * t61;
t33 = t331 * t61 + t336 * t67;
t176 = -mrSges(5,2) * t263 + mrSges(5,3) * t232;
t177 = mrSges(5,1) * t263 - mrSges(5,3) * t233;
t358 = t176 * t338 - t177 * t333;
t139 = -t196 * t336 + t197 * t331;
t140 = -t196 * t331 - t197 * t336;
t357 = t292 * t319 - t318 * t326;
t356 = t316 * t319 - t318 * t341;
t150 = -t257 * t339 - t334 * t258 + t285 * t405 - t286 * t404;
t266 = -t319 * t415 + t414;
t264 = t319 * t416 + t412;
t352 = t284 * t403 - t432;
t351 = t206 * t365;
t70 = -pkin(9) * t432 + pkin(4) * t277 + (-t219 + (pkin(9) * t284 - t212) * t333) * qJD(4) + t374;
t84 = t338 * t179 + t333 * t209 + t212 * t402 - t224 * t403;
t79 = -pkin(9) * t353 + t84;
t20 = t120 * t400 - t133 * t401 + t332 * t70 + t337 * t79;
t146 = -qJDD(3) * pkin(3) - t150;
t86 = -pkin(4) * t161 + t146;
t21 = -qJD(5) * t81 - t332 * t79 + t337 * t70;
t180 = qJD(2) * t284 + qJD(3) * t224;
t143 = pkin(4) * t353 + t180;
t317 = -qJDD(1) * pkin(1) + qJDD(2);
t270 = pkin(4) * t418 + t315 * t331;
t269 = -pkin(4) * t419 + t315 * t336;
t267 = t319 * t412 + t416;
t265 = -t319 * t414 + t415;
t259 = Ifges(4,4) * t274;
t254 = pkin(5) * t355 - t316;
t242 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t274;
t189 = t275 * Ifges(4,1) + t259 + t538;
t188 = Ifges(4,2) * t274 + t453 + t537;
t148 = mrSges(6,1) * t255 - t454;
t147 = -mrSges(6,2) * t255 + t455;
t141 = t468 + t464;
t138 = pkin(5) * t196 + t181;
t116 = -mrSges(5,2) * t205 + mrSges(5,3) * t161;
t115 = mrSges(5,1) * t205 - mrSges(5,3) * t160;
t109 = -mrSges(6,1) * t373 + mrSges(6,2) * t170;
t97 = t197 * t523 - t288 * t276;
t96 = -t221 * t284 - t276 * t355;
t89 = -mrSges(5,1) * t161 + mrSges(5,2) * t160;
t88 = mrSges(7,1) * t251 - mrSges(7,3) * t104;
t87 = -mrSges(7,2) * t251 + mrSges(7,3) * t375;
t85 = -qJD(4) * t159 + t374;
t82 = t160 * Ifges(5,4) + t161 * Ifges(5,2) + t205 * Ifges(5,6);
t73 = -pkin(5) * t97 + t143;
t64 = -mrSges(7,1) * t375 + mrSges(7,2) * t104;
t63 = -mrSges(6,2) * t202 + mrSges(6,3) * t77;
t62 = mrSges(6,1) * t202 - mrSges(6,3) * t76;
t42 = -qJD(6) * t140 - t331 * t96 + t336 * t97;
t41 = qJD(6) * t139 + t331 * t97 + t336 * t96;
t40 = -pkin(5) * t77 + t86;
t36 = -mrSges(6,1) * t77 + mrSges(6,2) * t76;
t25 = -mrSges(7,2) * t191 + mrSges(7,3) * t31;
t24 = mrSges(7,1) * t191 - mrSges(7,3) * t30;
t19 = t336 * t46 - t441;
t18 = -t331 * t46 - t439;
t15 = pkin(10) * t97 + t20;
t14 = pkin(5) * t277 - pkin(10) * t96 + t21;
t10 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t5 = -qJD(6) * t33 + t14 * t336 - t15 * t331;
t4 = qJD(6) * t32 + t14 * t331 + t15 * t336;
t1 = [(t131 * t352 - t132 * t353 - t430 * t57 - t431 * t56) * mrSges(5,3) + (t295 * mrSges(4,2) - t150 * mrSges(4,3) + Ifges(4,1) * t213 + Ifges(4,4) * t214 + Ifges(4,5) * qJDD(3) + t146 * t365 + t155 * t379 + t360 * t488 + t361 * t495 + t362 * t496) * t284 + (-Ifges(6,5) * t197 - Ifges(6,6) * t196) * t489 - (-Ifges(4,1) * t476 - Ifges(4,4) * t478 + t417 / 0.2e1 - t189 / 0.2e1 - t383 + t571) * t276 + (-m(4) * t215 - t572) * t180 - (t392 + t391 + t388) * t557 / 0.2e1 - (t295 * mrSges(4,1) - t149 * mrSges(4,3) - Ifges(4,4) * t213 + Ifges(5,5) * t496 + Ifges(6,5) * t507 + Ifges(7,5) * t515 - Ifges(4,2) * t214 - Ifges(4,6) * qJDD(3) + Ifges(5,6) * t495 + Ifges(6,6) * t506 + Ifges(7,6) * t514 + Ifges(5,3) * t488 + Ifges(6,3) * t489 + Ifges(7,3) * t490 + t550 + t551 + t552) * t557 + m(5) * (t131 * t85 + t132 * t84 + t158 * t57 + t159 * t56) + m(4) * (t149 * t224 + t179 * t216 - t295 * t311) + (Ifges(7,1) * t41 + Ifges(7,4) * t42) * t497 + (Ifges(7,1) * t140 + Ifges(7,4) * t139) * t515 + (Ifges(7,5) * t41 + Ifges(7,6) * t42) * t483 + (Ifges(7,5) * t140 + Ifges(7,6) * t139) * t490 + (-t12 * t196 + t13 * t197 - t59 * t96 + t60 * t97) * mrSges(6,3) + (Ifges(6,5) * t96 + Ifges(6,6) * t97) * t481 + t263 * (-Ifges(5,5) * t352 - Ifges(5,6) * t353) / 0.2e1 + (Ifges(6,1) * t96 + Ifges(6,4) * t97) * t491 + t86 * (mrSges(6,1) * t196 - mrSges(6,2) * t197) + (Ifges(7,4) * t41 + Ifges(7,2) * t42) * t499 + (Ifges(7,4) * t140 + Ifges(7,2) * t139) * t514 + (Ifges(6,4) * t96 + Ifges(6,2) * t97) * t493 + 0.2e1 * t406 * t307 * mrSges(3,3) + t232 * (-Ifges(5,4) * t352 - Ifges(5,2) * t353) / 0.2e1 + (-Ifges(5,1) * t352 - Ifges(5,4) * t353) * t485 + (-Ifges(6,1) * t197 - Ifges(6,4) * t196) * t507 + (-Ifges(6,4) * t197 - Ifges(6,2) * t196) * t506 - pkin(1) * t370 + t96 * t501 + t97 * t503 + t430 * t505 + t41 * t508 + t42 * t510 - t197 * t512 - t196 * t513 + (-Ifges(4,4) * t476 - Ifges(4,2) * t478 + Ifges(6,3) * t481 + Ifges(7,3) * t483 + Ifges(7,6) * t499 + Ifges(5,5) * t485 + Ifges(6,5) * t491 + Ifges(6,6) * t493 + Ifges(7,5) * t497 - t216 * mrSges(4,3) + t561 / 0.2e1 + t560 / 0.2e1 - t188 / 0.2e1 + t556 / 0.2e1 - t570) * t277 + m(7) * (t105 * t73 + t138 * t40 + t16 * t5 + t17 * t4 + t2 * t33 + t3 * t32) + m(6) * (t12 * t81 + t13 * t80 + t143 * t162 + t181 * t86 + t20 * t60 + t21 * t59) - (-m(4) * t150 + m(5) * t146 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t213 + t89) * t529 + t206 * (mrSges(5,1) * t353 - mrSges(5,2) * t352) + t317 * t369 + t140 * t518 + t139 * t519 + Ifges(2,3) * qJDD(1) + (-t416 * t517 - t267 * mrSges(5,1) - t250 * mrSges(6,1) - t237 * mrSges(7,1) - t266 * mrSges(5,2) - t249 * mrSges(6,2) - t236 * mrSges(7,2) + t555 * (t340 * t311 + t335 * t458) + t549 * t335 + (-m(5) * t371 - m(6) * t356 - m(7) * t357 - t548) * t340) * g(2) + m(3) * (-pkin(1) * t317 + (t307 + t396) * qJ(2) * t406) + (-t265 * mrSges(5,1) - t248 * mrSges(6,1) - t235 * mrSges(7,1) - t264 * mrSges(5,2) - t247 * mrSges(6,2) - t234 * mrSges(7,2) + (t458 * t555 + t549 - t564) * t340 + (m(4) * t311 - m(6) * (-t311 - t356) - m(5) * (-t311 - t371) - m(7) * (-t311 - t357) + t548) * t335) * g(1) + t105 * (-mrSges(7,1) * t42 + mrSges(7,2) * t41) + t4 * t87 + t5 * t88 + t81 * t63 + t80 * t62 + t73 * t64 + (Ifges(3,4) * t328 + Ifges(3,2) * t329) * t393 + (Ifges(3,1) * t328 + Ifges(3,4) * t329) * t394 - t82 * t431 / 0.2e1 + t32 * t24 + t33 * t25 - t154 * t384 / 0.2e1 - t311 * t378 + t138 * t10 + t40 * (-mrSges(7,1) * t139 + mrSges(7,2) * t140) + t143 * t109 + t20 * t147 + t21 * t148 + t158 * t115 + t159 * t116 + t162 * (-mrSges(6,1) * t97 + mrSges(6,2) * t96) + t84 * t176 + t85 * t177 + t181 * t36 + t224 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t214) + t179 * t242 + (t139 * t2 - t140 * t3 - t16 * t41 + t17 * t42) * mrSges(7,3); m(3) * t317 + t534 * t88 + t535 * t87 + (-t242 - t358) * t274 + t358 * qJD(4) + t378 + t532 * t148 + t531 * t147 + (-t109 - t64 + t530) * t275 + t370 + t338 * t115 + t333 * t116 - t355 * t62 + t288 * t63 + t217 * t24 + t218 * t25 + (-g(1) * t335 + g(2) * t340) * (m(3) - t555) - t382 * t406 * qJD(1) ^ 2 + (-t105 * t275 + t16 * t534 + t17 * t535 + t2 * t218 + t217 * t3) * m(7) + (t12 * t288 - t13 * t355 - t162 * t275 + t531 * t60 + t532 * t59) * m(6) + (-t206 * t275 + t333 * t56 + t338 * t57 + t263 * (-t131 * t333 + t132 * t338)) * m(5) + (t215 * t275 - t216 * t274 + t295) * m(4); (t351 + t383) * qJD(4) + t86 * (mrSges(6,1) * t355 + mrSges(6,2) * t288) - t355 * t513 + (-t12 * t355 - t13 * t288 + t182 * t60 - t183 * t59) * mrSges(6,3) + (Ifges(6,4) * t288 - Ifges(6,2) * t355) * t506 + (Ifges(6,1) * t288 - Ifges(6,4) * t355) * t507 + (Ifges(6,5) * t288 - Ifges(6,6) * t355) * t489 + (t360 * t480 + t361 * t487 + t362 * t486 - t351 + t571) * t274 + (t456 + t572) * t216 + (-t176 * t403 - t177 * t402 + m(5) * ((-t131 * t338 - t132 * t333) * qJD(4) + t524) + t338 * t116 - t333 * t115) * pkin(8) + (t131 * t434 + t132 * t435 + t524) * mrSges(5,3) + t527 * t109 + (-mrSges(6,1) * t532 + mrSges(6,2) * t531) * t162 + t533 * t64 + (-mrSges(7,1) * t534 + mrSges(7,2) * t535) * t105 - (Ifges(6,4) * t491 + Ifges(6,2) * t493 + Ifges(6,6) * t481 + t471 + t503) * t221 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t500 + (-Ifges(6,4) * t183 - Ifges(6,2) * t182) * t494 + (-Ifges(6,1) * t183 - Ifges(6,4) * t182) * t492 + (-Ifges(6,5) * t183 - Ifges(6,6) * t182) * t482 + (t232 * t361 + t233 * t362 + t263 * t360) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t275 + t189 + t259 + t413) * t274 / 0.2e1 - (Ifges(6,1) * t491 + Ifges(6,4) * t493 + Ifges(6,5) * t481 - t472 + t501) * t220 + (Ifges(7,4) * t497 + Ifges(7,2) * t499 + Ifges(7,6) * t483 + t469 + t510) * t122 + (Ifges(7,1) * t497 + Ifges(7,4) * t499 + Ifges(7,5) * t483 - t470 + t508) * t121 + (Ifges(7,5) * t127 + Ifges(7,6) * t126) * t484 + t188 * t476 + t417 * t478 - t183 * t502 - t182 * t504 + t333 * t505 + t127 * t509 + t126 * t511 + t288 * t512 + (g(1) * t340 + g(2) * t335) * ((mrSges(4,2) + t567) * t319 + (mrSges(4,1) - t568) * t318) + (t318 * t567 + t319 * t568 - t368) * g(3) + (Ifges(5,5) * t486 + Ifges(6,5) * t492 + Ifges(7,5) * t498 + Ifges(5,6) * t487 + Ifges(6,6) * t494 + Ifges(7,6) * t500 + Ifges(5,3) * t480 + Ifges(6,3) * t482 + Ifges(7,3) * t484 + t570) * t275 + (-pkin(3) * t146 - t131 * t151 - t132 * t152) * m(5) + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t498 + (Ifges(5,5) * t333 + Ifges(5,6) * t338) * t488 + (Ifges(7,5) * t218 + Ifges(7,6) * t217) * t490 + (Ifges(5,2) * t338 + t452) * t495 + (Ifges(5,1) * t333 + t451) * t496 + (-t126 * t17 + t127 * t16 + t2 * t217 - t218 * t3) * mrSges(7,3) + t146 * t366 + (Ifges(7,4) * t218 + Ifges(7,2) * t217) * t514 + (Ifges(7,1) * t218 + Ifges(7,4) * t217) * t515 + t218 * t518 + t217 * t519 - (Ifges(4,1) * t274 - t453 + t556) * t275 / 0.2e1 - pkin(3) * t89 - t403 * t447 - t402 * t448 + Ifges(4,3) * qJDD(3) + t338 * t82 / 0.2e1 - t316 * t36 + t136 * t24 + t137 * t25 - t149 * mrSges(4,2) + t150 * mrSges(4,1) + t539 * t147 + t540 * t148 + (t12 * t227 + t13 * t226 + t162 * t527 - t316 * t86 + t539 * t60 + t540 * t59) * m(6) - t152 * t176 - t151 * t177 + t154 * t379 + Ifges(4,5) * t213 + Ifges(4,6) * t214 + t544 * t88 + t545 * t87 + (t105 * t533 + t136 * t3 + t137 * t2 + t16 * t544 + t17 * t545 + t254 * t40) * m(7) + t40 * (-mrSges(7,1) * t217 + mrSges(7,2) * t218) + t226 * t62 + t227 * t63 - t215 * t242 + t254 * t10; -(mrSges(6,1) * t162 + Ifges(6,4) * t492 + Ifges(6,2) * t494 + Ifges(6,6) * t482 - t471 + t504) * t170 + (-Ifges(5,2) * t233 + t155 + t225) * t487 + (-mrSges(6,2) * t162 + Ifges(6,1) * t492 + Ifges(6,4) * t494 + Ifges(6,5) * t482 + t472 + t502) * t373 + (t147 * t400 - t148 * t401 + t332 * t63) * pkin(4) + t551 + (Ifges(5,5) * t232 - Ifges(5,6) * t233) * t480 + t569 + t154 * t485 + (Ifges(5,1) * t232 - t445) * t486 + t388 + (t12 * t332 + t13 * t337 + (-t332 * t59 + t337 * t60) * qJD(5)) * t517 + (t365 + t528 + t563 + t564) * t460 - t109 * t468 - m(6) * (t162 * t468 + t59 * t65 + t60 * t66) + t233 * t447 + t232 * t448 + t62 * t466 - t141 * t64 - t66 * t147 - t65 * t148 + (-m(7) * (-t297 * t428 - t298 * t340) - mrSges(5,2) * t265 + t541 * t264 + t525) * g(2) + (-m(7) * (-t297 * t427 + t298 * t335) + mrSges(5,2) * t267 - t541 * t266 + t526) * g(1) + t542 * t87 + t543 * t88 + (-t105 * t141 + t16 * t543 + t17 * t542 + t2 * t270 + t269 * t3) * m(7) - t131 * t176 + t132 * t177 - t206 * (mrSges(5,1) * t233 + mrSges(5,2) * t232) + t269 * t24 + t270 * t25; (Ifges(6,5) * t373 - Ifges(6,6) * t170) * t482 + t92 * t491 + (Ifges(6,1) * t373 - t450) * t492 + (t2 * t331 + t3 * t336 + (-t16 * t331 + t17 * t336) * qJD(6)) * t516 - t64 * t464 - m(7) * (t105 * t464 + t16 * t18 + t17 * t19) - t19 * t87 - t18 * t88 - t162 * (mrSges(6,1) * t170 + mrSges(6,2) * t373) + (t454 + t148) * t60 + (t455 - t147) * t59 + (-Ifges(6,2) * t170 + t166 + t93) * t494 + (m(7) * t463 + t528) * t460 + (t247 * t516 + t525) * g(2) + (-t249 * t516 + t526) * g(1) + (t24 * t336 + t25 * t331 + t398 * t87 - t399 * t88) * pkin(5) + t569; -t105 * (mrSges(7,1) * t104 + mrSges(7,2) * t375) + (Ifges(7,1) * t375 - t449) * t498 + t54 * t497 + (Ifges(7,5) * t375 - Ifges(7,6) * t104) * t484 - t16 * t87 + t17 * t88 - g(1) * t409 - g(2) * t410 - t363 * t460 + (t104 * t17 + t16 * t375) * mrSges(7,3) + t354 + (-Ifges(7,2) * t104 + t55 + t98) * t500;];
tau  = t1;
