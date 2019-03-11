% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:52
% EndTime: 2019-03-09 10:16:49
% DurationCPUTime: 34.69s
% Computational Cost: add. (20435->915), mult. (46894->1200), div. (0->0), fcn. (35721->18), ass. (0->412)
t339 = sin(qJ(2));
t343 = cos(qJ(2));
t440 = sin(pkin(10));
t381 = qJD(1) * t440;
t441 = cos(pkin(10));
t382 = qJD(1) * t441;
t273 = -t339 * t381 + t343 * t382;
t274 = -t339 * t382 - t343 * t381;
t413 = qJD(1) * t339;
t401 = pkin(2) * t413;
t184 = -pkin(3) * t274 - pkin(8) * t273 + t401;
t336 = -qJ(3) - pkin(7);
t307 = t336 * t343;
t296 = qJD(1) * t307;
t278 = t440 * t296;
t305 = t336 * t339;
t295 = qJD(1) * t305;
t218 = t295 * t441 + t278;
t338 = sin(qJ(4));
t342 = cos(qJ(4));
t136 = t342 * t184 - t218 * t338;
t392 = t440 * pkin(2);
t316 = t392 + pkin(8);
t416 = qJ(5) + t316;
t377 = qJD(4) * t416;
t437 = t273 * t342;
t573 = pkin(4) * t274 + qJ(5) * t437 - qJD(5) * t338 - t342 * t377 - t136;
t137 = t338 * t184 + t342 * t218;
t438 = t273 * t338;
t572 = -qJ(5) * t438 - qJD(5) * t342 + t338 * t377 + t137;
t333 = sin(pkin(11));
t334 = cos(pkin(11));
t291 = t333 * t342 + t334 * t338;
t174 = t291 * t273;
t272 = t291 * qJD(4);
t538 = t272 - t174;
t360 = t333 * t338 - t334 * t342;
t175 = t360 * t273;
t276 = t360 * qJD(4);
t571 = t276 - t175;
t553 = t572 * t333 + t334 * t573;
t552 = t333 * t573 - t572 * t334;
t408 = qJD(4) * t338;
t570 = t408 - t438;
t498 = m(6) + m(7);
t530 = t498 + m(4) + m(5);
t557 = Ifges(6,3) + Ifges(5,3);
t568 = pkin(5) * t274 + t571 * pkin(9) + t553;
t567 = -t538 * pkin(9) + t552;
t331 = qJ(4) + pkin(11);
t325 = cos(t331);
t328 = t342 * pkin(4);
t300 = pkin(5) * t325 + t328;
t294 = pkin(3) + t300;
t327 = qJ(6) + t331;
t314 = sin(t327);
t315 = cos(t327);
t319 = t328 + pkin(3);
t323 = sin(t331);
t372 = -mrSges(5,1) * t342 + mrSges(5,2) * t338;
t566 = -m(5) * pkin(3) - m(6) * t319 - m(7) * t294 - mrSges(6,1) * t325 - mrSges(7,1) * t315 + mrSges(6,2) * t323 + mrSges(7,2) * t314 + t372;
t565 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t234 = qJD(2) * t342 + t274 * t338;
t235 = qJD(2) * t338 - t274 * t342;
t162 = t234 * t333 + t235 * t334;
t564 = pkin(9) * t162;
t385 = t441 * t296;
t217 = t295 * t440 - t385;
t541 = pkin(4) * t570 - t217;
t306 = -mrSges(3,1) * t343 + mrSges(3,2) * t339;
t332 = qJ(2) + pkin(10);
t324 = sin(t332);
t326 = cos(t332);
t563 = -mrSges(4,1) * t326 + t324 * t565 + t306;
t290 = t339 * t440 - t343 * t441;
t277 = t290 * qJD(2);
t292 = t339 * t441 + t343 * t440;
t407 = qJD(4) * t342;
t391 = t292 * t407;
t356 = -t277 * t338 + t391;
t340 = sin(qJ(1));
t344 = cos(qJ(1));
t562 = g(1) * t344 + g(2) * t340;
t337 = sin(qJ(6));
t341 = cos(qJ(6));
t379 = t334 * t234 - t235 * t333;
t561 = -t162 * t337 + t341 * t379;
t101 = t162 * t341 + t337 * t379;
t516 = m(6) * pkin(4);
t406 = qJD(1) * qJD(2);
t389 = t339 * t406;
t405 = qJDD(1) * t343;
t297 = -t389 + t405;
t298 = qJDD(1) * t339 + t343 * t406;
t220 = t297 * t440 + t298 * t441;
t156 = qJD(4) * t234 + qJDD(2) * t338 + t220 * t342;
t157 = -qJD(4) * t235 + qJDD(2) * t342 - t220 * t338;
t85 = -t156 * t333 + t157 * t334;
t86 = t156 * t334 + t157 * t333;
t30 = qJD(6) * t561 + t337 * t85 + t341 * t86;
t515 = t30 / 0.2e1;
t31 = -qJD(6) * t101 - t337 * t86 + t341 * t85;
t514 = t31 / 0.2e1;
t506 = t85 / 0.2e1;
t505 = t86 / 0.2e1;
t495 = t156 / 0.2e1;
t494 = t157 / 0.2e1;
t219 = t297 * t441 - t440 * t298;
t216 = qJDD(4) - t219;
t211 = qJDD(6) + t216;
t489 = t211 / 0.2e1;
t488 = t216 / 0.2e1;
t559 = t297 / 0.2e1;
t558 = pkin(9) * t379;
t460 = qJD(2) / 0.2e1;
t285 = t416 * t338;
t286 = t416 * t342;
t196 = -t334 * t285 - t286 * t333;
t167 = -pkin(9) * t291 + t196;
t197 = -t333 * t285 + t334 * t286;
t168 = -pkin(9) * t360 + t197;
t111 = t167 * t337 + t168 * t341;
t556 = -qJD(6) * t111 - t337 * t567 + t341 * t568;
t110 = t167 * t341 - t168 * t337;
t555 = qJD(6) * t110 + t337 * t568 + t341 * t567;
t329 = t343 * pkin(2);
t320 = t329 + pkin(1);
t554 = t234 * Ifges(5,6);
t442 = qJDD(2) / 0.2e1;
t467 = pkin(4) * t334;
t317 = pkin(5) + t467;
t468 = pkin(4) * t333;
t268 = t317 * t341 - t337 * t468;
t301 = -qJD(1) * t320 + qJD(3);
t173 = -pkin(3) * t273 + pkin(8) * t274 + t301;
t284 = qJD(2) * pkin(2) + t295;
t208 = t440 * t284 - t385;
t195 = qJD(2) * pkin(8) + t208;
t122 = t342 * t173 - t195 * t338;
t108 = -qJ(5) * t235 + t122;
t123 = t173 * t338 + t195 * t342;
t109 = qJ(5) * t234 + t123;
t428 = t334 * t109;
t60 = -t108 * t333 - t428;
t41 = t60 - t558;
t102 = t333 * t109;
t61 = t334 * t108 - t102;
t42 = t61 - t564;
t551 = t268 * qJD(6) - t337 * t41 - t341 * t42;
t269 = t317 * t337 + t341 * t468;
t550 = -t269 * qJD(6) + t337 * t42 - t341 * t41;
t549 = t516 + mrSges(5,1);
t548 = Ifges(4,5) * qJD(2);
t547 = Ifges(4,6) * qJD(2);
t546 = t324 * t562;
t121 = -t174 * t337 - t175 * t341;
t212 = -t291 * t337 - t341 * t360;
t145 = qJD(6) * t212 - t272 * t337 - t276 * t341;
t543 = t145 - t121;
t120 = -t174 * t341 + t175 * t337;
t213 = t291 * t341 - t337 * t360;
t146 = -qJD(6) * t213 - t272 * t341 + t276 * t337;
t542 = t146 - t120;
t206 = pkin(3) * t290 - pkin(8) * t292 - t320;
t224 = t305 * t440 - t307 * t441;
t221 = t342 * t224;
t153 = t338 * t206 + t221;
t540 = pkin(5) * t538 + t541;
t539 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t234 - mrSges(5,2) * t235 + mrSges(4,3) * t274;
t536 = Ifges(5,5) * t156 + Ifges(6,5) * t86 + Ifges(5,6) * t157 + Ifges(6,6) * t85 + t216 * t557;
t321 = pkin(7) * t405;
t287 = -pkin(7) * t389 + t321;
t288 = t298 * pkin(7);
t535 = t287 * t343 + t288 * t339;
t439 = qJDD(1) * pkin(1);
t257 = -pkin(2) * t297 + qJDD(3) - t439;
t134 = -pkin(3) * t219 - pkin(8) * t220 + t257;
t410 = qJD(3) * t339;
t202 = qJDD(2) * pkin(2) - qJ(3) * t298 - qJD(1) * t410 - t288;
t411 = qJD(2) * t339;
t398 = pkin(7) * t411;
t409 = qJD(3) * t343;
t214 = qJ(3) * t297 + t321 + (-t398 + t409) * qJD(1);
t144 = t440 * t202 + t441 * t214;
t140 = qJDD(2) * pkin(8) + t144;
t52 = t338 * t134 + t342 * t140 + t173 * t407 - t195 * t408;
t53 = -qJD(4) * t123 + t342 * t134 - t140 * t338;
t534 = -t338 * t53 + t342 * t52;
t262 = qJD(4) - t273;
t533 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t254 = qJD(6) + t262;
t532 = t235 * Ifges(5,5) + t162 * Ifges(6,5) + t101 * Ifges(7,5) + Ifges(6,6) * t379 + t561 * Ifges(7,6) + t254 * Ifges(7,3) + t262 * t557 + t554;
t531 = 0.2e1 * t442;
t93 = pkin(4) * t262 + t108;
t54 = t334 * t93 - t102;
t37 = pkin(5) * t262 + t54 - t564;
t55 = t333 * t93 + t428;
t40 = t55 + t558;
t14 = t337 * t37 + t341 * t40;
t207 = t284 * t441 + t278;
t194 = -qJD(2) * pkin(3) - t207;
t158 = -t234 * pkin(4) + qJD(5) + t194;
t94 = -pkin(5) * t379 + t158;
t529 = -mrSges(7,1) * t94 + mrSges(7,3) * t14;
t13 = -t337 * t40 + t341 * t37;
t528 = mrSges(7,2) * t94 - mrSges(7,3) * t13;
t527 = -t301 * mrSges(4,2) + t207 * mrSges(4,3);
t34 = pkin(4) * t216 - qJ(5) * t156 - qJD(5) * t235 + t53;
t36 = qJ(5) * t157 + qJD(5) * t234 + t52;
t11 = -t333 * t36 + t334 * t34;
t6 = pkin(5) * t216 - pkin(9) * t86 + t11;
t12 = t333 * t34 + t334 * t36;
t7 = pkin(9) * t85 + t12;
t2 = qJD(6) * t13 + t337 * t6 + t341 * t7;
t3 = -qJD(6) * t14 - t337 * t7 + t341 * t6;
t525 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t466 = pkin(4) * t338;
t299 = pkin(5) * t323 + t466;
t524 = -m(6) * t466 - m(7) * t299;
t522 = m(3) * pkin(1) + mrSges(2,1) - t563;
t521 = t53 * mrSges(5,1) + t11 * mrSges(6,1) - t52 * mrSges(5,2) - t12 * mrSges(6,2);
t520 = t301 * mrSges(4,1) + t122 * mrSges(5,1) + t54 * mrSges(6,1) + t13 * mrSges(7,1) - t123 * mrSges(5,2) - t55 * mrSges(6,2) - t14 * mrSges(7,2) - t208 * mrSges(4,3);
t518 = Ifges(7,4) * t515 + Ifges(7,2) * t514 + Ifges(7,6) * t489;
t517 = Ifges(7,1) * t515 + Ifges(7,4) * t514 + Ifges(7,5) * t489;
t513 = Ifges(6,4) * t505 + Ifges(6,2) * t506 + Ifges(6,6) * t488;
t512 = Ifges(6,1) * t505 + Ifges(6,4) * t506 + Ifges(6,5) * t488;
t453 = Ifges(7,4) * t101;
t50 = Ifges(7,2) * t561 + Ifges(7,6) * t254 + t453;
t511 = -t50 / 0.2e1;
t510 = t50 / 0.2e1;
t95 = Ifges(7,4) * t561;
t51 = Ifges(7,1) * t101 + Ifges(7,5) * t254 + t95;
t509 = -t51 / 0.2e1;
t508 = t51 / 0.2e1;
t507 = Ifges(5,1) * t495 + Ifges(5,4) * t494 + Ifges(5,5) * t488;
t90 = t162 * Ifges(6,4) + Ifges(6,2) * t379 + t262 * Ifges(6,6);
t504 = -t90 / 0.2e1;
t503 = t90 / 0.2e1;
t91 = t162 * Ifges(6,1) + Ifges(6,4) * t379 + t262 * Ifges(6,5);
t502 = -t91 / 0.2e1;
t501 = t91 / 0.2e1;
t500 = -t561 / 0.2e1;
t499 = t561 / 0.2e1;
t497 = -t101 / 0.2e1;
t496 = t101 / 0.2e1;
t493 = -t379 / 0.2e1;
t492 = t379 / 0.2e1;
t491 = -t162 / 0.2e1;
t490 = t162 / 0.2e1;
t487 = -t234 / 0.2e1;
t486 = -t235 / 0.2e1;
t485 = t235 / 0.2e1;
t484 = -t254 / 0.2e1;
t483 = t254 / 0.2e1;
t482 = -t262 / 0.2e1;
t481 = t262 / 0.2e1;
t479 = t273 / 0.2e1;
t478 = -t274 / 0.2e1;
t474 = mrSges(6,3) * t54;
t473 = mrSges(6,3) * t55;
t469 = pkin(4) * t235;
t465 = pkin(7) * t343;
t464 = pkin(8) * t324;
t461 = g(3) * t324;
t335 = -qJ(5) - pkin(8);
t275 = t292 * qJD(2);
t359 = qJ(5) * t277 - qJD(5) * t292;
t386 = qJD(2) * t336;
t270 = t339 * t386 + t409;
t271 = t343 * t386 - t410;
t183 = t270 * t441 + t271 * t440;
t400 = pkin(2) * t411;
t185 = pkin(3) * t275 + pkin(8) * t277 + t400;
t380 = -t183 * t338 + t342 * t185;
t65 = pkin(4) * t275 + t359 * t342 + (-t221 + (qJ(5) * t292 - t206) * t338) * qJD(4) + t380;
t394 = t342 * t183 + t338 * t185 + t206 * t407;
t69 = -qJ(5) * t391 + (-qJD(4) * t224 + t359) * t338 + t394;
t24 = t333 * t65 + t334 * t69;
t458 = mrSges(5,3) * t235;
t457 = Ifges(3,4) * t339;
t456 = Ifges(3,4) * t343;
t455 = Ifges(5,4) * t338;
t454 = Ifges(5,4) * t342;
t452 = t122 * mrSges(5,3);
t449 = t235 * Ifges(5,4);
t448 = t274 * Ifges(4,4);
t433 = t292 * t338;
t432 = t292 * t342;
t330 = -pkin(9) + t335;
t431 = t324 * t330;
t430 = t324 * t335;
t429 = t326 * t344;
t148 = t234 * Ifges(5,2) + t262 * Ifges(5,6) + t449;
t427 = t338 * t148;
t426 = t338 * t344;
t425 = t340 * t299;
t424 = t340 * t314;
t423 = t340 * t315;
t422 = t340 * t323;
t421 = t340 * t325;
t420 = t340 * t338;
t419 = t340 * t342;
t225 = Ifges(5,4) * t234;
t149 = t235 * Ifges(5,1) + t262 * Ifges(5,5) + t225;
t418 = t342 * t149;
t417 = t342 * t344;
t152 = t342 * t206 - t224 * t338;
t117 = pkin(4) * t290 - qJ(5) * t432 + t152;
t127 = -qJ(5) * t433 + t153;
t71 = t333 * t117 + t334 * t127;
t236 = t315 * t344 + t326 * t424;
t237 = t314 * t344 - t326 * t423;
t415 = -t236 * mrSges(7,1) + t237 * mrSges(7,2);
t238 = -t314 * t429 + t423;
t239 = t315 * t429 + t424;
t414 = t238 * mrSges(7,1) - t239 * mrSges(7,2);
t412 = qJD(1) * t343;
t404 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t413) * t465;
t403 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t211;
t393 = t441 * pkin(2);
t390 = t418 / 0.2e1;
t45 = -t85 * mrSges(6,1) + t86 * mrSges(6,2);
t10 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t387 = -t408 / 0.2e1;
t23 = -t333 * t69 + t334 * t65;
t383 = -t219 * mrSges(4,1) + t220 * mrSges(4,2);
t70 = t334 * t117 - t127 * t333;
t318 = -t393 - pkin(3);
t376 = pkin(3) * t326 + t464;
t374 = mrSges(3,1) * t339 + mrSges(3,2) * t343;
t371 = mrSges(5,1) * t338 + mrSges(5,2) * t342;
t370 = -mrSges(7,1) * t314 - mrSges(7,2) * t315;
t369 = Ifges(5,1) * t342 - t455;
t368 = t343 * Ifges(3,2) + t457;
t367 = -Ifges(5,2) * t338 + t454;
t366 = Ifges(3,5) * t343 - Ifges(3,6) * t339;
t365 = Ifges(5,5) * t342 - Ifges(5,6) * t338;
t189 = t360 * t292;
t58 = pkin(5) * t290 + pkin(9) * t189 + t70;
t188 = t291 * t292;
t64 = -pkin(9) * t188 + t71;
t21 = -t337 * t64 + t341 * t58;
t22 = t337 * t58 + t341 * t64;
t143 = t202 * t441 - t440 * t214;
t182 = t270 * t440 - t441 * t271;
t223 = -t441 * t305 - t307 * t440;
t171 = -mrSges(5,2) * t262 + mrSges(5,3) * t234;
t172 = mrSges(5,1) * t262 - t458;
t363 = t171 * t342 - t172 * t338;
t130 = -t188 * t341 + t189 * t337;
t131 = -t188 * t337 - t189 * t341;
t362 = t294 * t326 - t431;
t361 = t319 * t326 - t430;
t302 = -t328 + t318;
t358 = t403 + t525;
t181 = pkin(4) * t433 + t223;
t357 = pkin(1) * t374;
t265 = -t326 * t426 + t419;
t263 = t326 * t420 + t417;
t355 = t277 * t342 + t292 * t408;
t354 = t194 * t371;
t353 = t339 * (Ifges(3,1) * t343 - t457);
t138 = pkin(4) * t356 + t182;
t139 = -qJDD(2) * pkin(3) - t143;
t79 = -t157 * pkin(4) + qJDD(5) + t139;
t322 = Ifges(3,4) * t412;
t304 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t412;
t282 = Ifges(3,1) * t413 + Ifges(3,5) * qJD(2) + t322;
t281 = Ifges(3,6) * qJD(2) + qJD(1) * t368;
t266 = t326 * t417 + t420;
t264 = -t326 * t419 + t426;
t261 = Ifges(4,4) * t273;
t247 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t273;
t246 = t325 * t429 + t422;
t245 = -t323 * t429 + t421;
t244 = t323 * t344 - t326 * t421;
t243 = t325 * t344 + t326 * t422;
t233 = pkin(5) * t360 + t302;
t203 = -mrSges(4,1) * t273 - mrSges(4,2) * t274;
t199 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t220;
t198 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t219;
t191 = -t274 * Ifges(4,1) + t261 + t548;
t190 = t273 * Ifges(4,2) - t448 + t547;
t142 = mrSges(6,1) * t262 - mrSges(6,3) * t162;
t141 = -mrSges(6,2) * t262 + mrSges(6,3) * t379;
t133 = t188 * pkin(5) + t181;
t132 = pkin(5) * t162 + t469;
t129 = -t272 * t292 + t277 * t360;
t128 = t276 * t292 + t277 * t291;
t116 = -mrSges(5,2) * t216 + mrSges(5,3) * t157;
t115 = mrSges(5,1) * t216 - mrSges(5,3) * t156;
t107 = -mrSges(6,1) * t379 + mrSges(6,2) * t162;
t92 = -mrSges(5,1) * t157 + mrSges(5,2) * t156;
t81 = mrSges(7,1) * t254 - mrSges(7,3) * t101;
t80 = -mrSges(7,2) * t254 + mrSges(7,3) * t561;
t78 = -qJD(4) * t153 + t380;
t77 = -t224 * t408 + t394;
t76 = -t128 * pkin(5) + t138;
t74 = t156 * Ifges(5,4) + t157 * Ifges(5,2) + t216 * Ifges(5,6);
t73 = mrSges(6,1) * t216 - mrSges(6,3) * t86;
t72 = -mrSges(6,2) * t216 + mrSges(6,3) * t85;
t59 = -mrSges(7,1) * t561 + mrSges(7,2) * t101;
t48 = -qJD(6) * t131 + t128 * t341 - t129 * t337;
t47 = qJD(6) * t130 + t128 * t337 + t129 * t341;
t43 = -t85 * pkin(5) + t79;
t26 = -mrSges(7,2) * t211 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t211 - mrSges(7,3) * t30;
t20 = pkin(9) * t128 + t24;
t19 = pkin(5) * t275 - pkin(9) * t129 + t23;
t5 = -qJD(6) * t22 + t19 * t341 - t20 * t337;
t4 = qJD(6) * t21 + t19 * t337 + t20 * t341;
t1 = [t368 * t559 + (t257 * mrSges(4,2) - t143 * mrSges(4,3) + Ifges(4,1) * t220 + Ifges(4,4) * t219 + Ifges(4,5) * t531 + t139 * t371 + t149 * t387 + t365 * t488 + t367 * t494 + t369 * t495) * t292 + (t297 * t465 + t535) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t535) + (t403 + t536) * t290 / 0.2e1 + (-m(4) * t207 + m(5) * t194 - t539) * t182 + (t532 / 0.2e1 - t190 / 0.2e1 + t554 / 0.2e1 + t557 * t481 + Ifges(7,6) * t499 - Ifges(4,6) * t460 - Ifges(4,4) * t478 - Ifges(4,2) * t479 + Ifges(7,3) * t483 + Ifges(5,5) * t485 + Ifges(6,5) * t490 + Ifges(6,6) * t492 + Ifges(7,5) * t496 + t520) * t275 + (Ifges(6,1) * t129 + Ifges(6,4) * t128) * t490 + (-Ifges(5,5) * t355 + Ifges(6,5) * t129 - Ifges(5,6) * t356 + Ifges(6,6) * t128) * t481 + (Ifges(7,5) * t131 + Ifges(7,6) * t130) * t489 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t483 + t79 * (mrSges(6,1) * t188 - mrSges(6,2) * t189) - t320 * t383 - t148 * t391 / 0.2e1 + (-Ifges(5,1) * t355 - Ifges(5,4) * t356) * t485 - t304 * t398 + (-t420 * t516 - m(7) * t425 - t266 * mrSges(5,1) - t246 * mrSges(6,1) - t239 * mrSges(7,1) - t265 * mrSges(5,2) - t245 * mrSges(6,2) - t238 * mrSges(7,2) - t530 * (t344 * t320 - t340 * t336) + t533 * t340 + (-m(5) * t376 - m(6) * t361 - m(7) * t362 - t522) * t344) * g(2) + (-t13 * t47 + t130 * t2 - t131 * t3 + t14 * t48) * mrSges(7,3) + m(4) * (t144 * t224 + t183 * t208 - t257 * t320 + t301 * t400) + m(5) * (t122 * t78 + t123 * t77 + t152 * t53 + t153 * t52) - t357 * t406 + (t122 * t355 - t123 * t356 - t432 * t53 - t433 * t52) * mrSges(5,3) + (t353 + t343 * (-Ifges(3,2) * t339 + t456)) * t406 / 0.2e1 + t131 * t517 + t130 * t518 + t129 * t501 + t128 * t503 + t432 * t507 + t47 * t508 + t48 * t510 - t189 * t512 - t188 * t513 + (Ifges(7,1) * t131 + Ifges(7,4) * t130) * t515 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t496 + (-Ifges(6,1) * t189 - Ifges(6,4) * t188) * t505 + (t11 * t189 - t12 * t188 + t128 * t55 - t129 * t54) * mrSges(6,3) + (-Ifges(6,5) * t189 - Ifges(6,6) * t188) * t488 + (-Ifges(6,4) * t189 - Ifges(6,2) * t188) * t506 + Ifges(3,6) * t343 * t442 - t281 * t411 / 0.2e1 + (-m(4) * t143 + m(5) * t139 - t199 + t92) * t223 + t298 * t456 / 0.2e1 + m(7) * (t13 * t5 + t133 * t43 + t14 * t4 + t2 * t22 + t21 * t3 + t76 * t94) + m(6) * (t11 * t70 + t12 * t71 + t138 * t158 + t181 * t79 + t23 * t54 + t24 * t55) + (Ifges(7,4) * t131 + Ifges(7,2) * t130) * t514 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t499 + (-t264 * mrSges(5,1) - t244 * mrSges(6,1) - t237 * mrSges(7,1) - t263 * mrSges(5,2) - t243 * mrSges(6,2) - t236 * mrSges(7,2) + (t336 * t530 + t524 + t533) * t344 + (-m(6) * (-t320 - t361) - m(7) * (-t320 - t362) - m(5) * (-t320 - t376) + m(4) * t320 + t522) * t340) * g(1) + t343 * t282 * t460 - t306 * t439 - t74 * t433 / 0.2e1 + (Ifges(3,1) * t298 + Ifges(3,4) * t559 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t298) + t531 * Ifges(3,5)) * t339 - qJDD(2) * mrSges(3,2) * t465 + (t257 * mrSges(4,1) - t144 * mrSges(4,3) - Ifges(4,4) * t220 + Ifges(5,5) * t495 + Ifges(6,5) * t505 + Ifges(7,5) * t515 - Ifges(4,2) * t219 - t531 * Ifges(4,6) + Ifges(5,6) * t494 + Ifges(6,6) * t506 + Ifges(7,6) * t514 + Ifges(7,3) * t489 + t557 * t488 + t521 + t525) * t290 + (t366 * t460 - t404) * qJD(2) + Ifges(2,3) * qJDD(1) + (Ifges(6,4) * t129 + Ifges(6,2) * t128) * t492 + t234 * (-Ifges(5,4) * t355 - Ifges(5,2) * t356) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t297 + mrSges(3,2) * t298) + t183 * t247 + t224 * t198 + t194 * (mrSges(5,1) * t356 - mrSges(5,2) * t355) + t181 * t45 + t77 * t171 + t78 * t172 + t158 * (-mrSges(6,1) * t128 + mrSges(6,2) * t129) + t152 * t115 + t153 * t116 + t24 * t141 + t23 * t142 + t138 * t107 + t43 * (-mrSges(7,1) * t130 + mrSges(7,2) * t131) + t133 * t10 + t94 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + t4 * t80 + t5 * t81 + t76 * t59 + t71 * t72 + t70 * t73 + t203 * t400 + t343 * (Ifges(3,4) * t298 + Ifges(3,2) * t297 + Ifges(3,6) * qJDD(2)) / 0.2e1 - (-t427 / 0.2e1 + t390 + t191 / 0.2e1 + Ifges(4,5) * t460 + Ifges(4,1) * t478 + Ifges(4,4) * t479 - t527) * t277 + t21 * t25 + t22 * t26; (-t122 * t136 - t123 * t137 + t139 * t318 - t194 * t217) * m(5) - (Ifges(6,4) * t490 + Ifges(6,2) * t492 + Ifges(6,6) * t481 + t473 + t503) * t272 + (-t120 * t14 + t121 * t13 + t2 * t212 - t213 * t3) * mrSges(7,3) + t540 * t59 + t541 * t107 + (Ifges(7,1) * t496 + Ifges(7,4) * t499 + Ifges(7,5) * t483 + t508 + t528) * t145 + (Ifges(7,4) * t496 + Ifges(7,2) * t499 + Ifges(7,6) * t483 + t510 + t529) * t146 + (t404 + (-t353 / 0.2e1 + t357) * qJD(1)) * qJD(1) + (t281 / 0.2e1 + pkin(7) * t304) * t413 + (m(5) * ((-t122 * t342 - t123 * t338) * qJD(4) + t534) - t172 * t407 - t171 * t408 - t338 * t115 + t342 * t116) * t316 + t539 * t217 - (Ifges(6,1) * t490 + Ifges(6,4) * t492 + Ifges(6,5) * t481 - t474 + t501) * t276 + (-t11 * t291 - t12 * t360 + t174 * t55 - t175 * t54) * mrSges(6,3) + (Ifges(6,1) * t291 - Ifges(6,4) * t360) * t505 + (Ifges(6,4) * t291 - Ifges(6,2) * t360) * t506 + t79 * (mrSges(6,1) * t360 + mrSges(6,2) * t291) + (Ifges(5,5) * t338 + Ifges(6,5) * t291 + Ifges(5,6) * t342 - Ifges(6,6) * t360) * t488 - t360 * t513 + (Ifges(7,4) * t121 + Ifges(7,2) * t120) * t500 + t318 * t92 + (t390 + t354) * qJD(4) + (-Ifges(6,4) * t175 - Ifges(6,2) * t174) * t493 + (mrSges(6,1) * t538 - mrSges(6,2) * t571) * t158 + (Ifges(7,1) * t121 + Ifges(7,4) * t120) * t497 - t366 * t406 / 0.2e1 - (-Ifges(3,2) * t413 + t282 + t322) * t412 / 0.2e1 + (t234 * t367 + t235 * t369 + t262 * t365) * qJD(4) / 0.2e1 - (Ifges(4,2) * t274 + t191 + t261 + t418) * t273 / 0.2e1 - t407 * t452 + t555 * t80 + (t110 * t3 + t111 * t2 + t13 * t556 + t14 * t555 + t233 * t43 + t540 * t94) * m(7) + t213 * t517 + t212 * t518 - t175 * t502 - t174 * t504 + t338 * t507 + t121 * t509 + t120 * t511 + t291 * t512 + (Ifges(7,4) * t213 + Ifges(7,2) * t212) * t514 + (Ifges(7,1) * t213 + Ifges(7,4) * t212) * t515 + (Ifges(7,5) * t213 + Ifges(7,6) * t212) * t489 + (Ifges(5,2) * t342 + t455) * t494 + (Ifges(5,1) * t338 + t454) * t495 + (-Ifges(6,5) * t175 - Ifges(6,6) * t174) * t482 + (-Ifges(6,1) * t175 - Ifges(6,4) * t174) * t491 - t203 * t401 + (t369 * t486 + t367 * t487 - t548 / 0.2e1 - t354 + t365 * t482 + t527) * t273 + t552 * t141 + t553 * t142 + (t11 * t196 + t12 * t197 + t158 * t541 + t302 * t79 + t54 * t553 + t55 * t552) * m(6) + (-m(5) * (t329 + t464) - m(7) * (t329 - t431) - m(6) * (t329 - t430) - m(4) * t329 + t566 * t326 + t563) * g(3) + (Ifges(7,5) * t121 + Ifges(7,6) * t120) * t484 + (Ifges(4,1) * t273 + t448 + t532) * t274 / 0.2e1 + t148 * t387 + t556 * t81 + (-t557 * t482 - Ifges(6,5) * t491 - Ifges(6,6) * t493 - Ifges(7,5) * t497 - Ifges(7,6) * t500 - Ifges(5,5) * t486 - Ifges(5,6) * t487 - Ifges(7,3) * t484 - t547 / 0.2e1 + t520) * t274 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t302 * t45 + Ifges(3,6) * t297 + Ifges(3,5) * t298 + t139 * t372 - t287 * mrSges(3,2) - t288 * mrSges(3,1) - t218 * t247 + t233 * t10 + Ifges(4,6) * t219 + Ifges(4,5) * t220 + t43 * (-mrSges(7,1) * t212 + mrSges(7,2) * t213) + t197 * t72 + t196 * t73 - t137 * t171 - t136 * t172 + t143 * mrSges(4,1) - t144 * mrSges(4,2) - t94 * (-mrSges(7,1) * t120 + mrSges(7,2) * t121) + t110 * t25 + t111 * t26 + t562 * (t374 + t530 * pkin(2) * t339 + (-m(5) * pkin(8) + m(6) * t335 + m(7) * t330 + t565) * t326 + (mrSges(4,1) - t566) * t324) + (t122 * t437 - t123 * t570 + t534) * mrSges(5,3) + t342 * t74 / 0.2e1 + t198 * t392 + t199 * t393 + t190 * t478 + t427 * t479 + (t207 * t217 - t208 * t218 - t301 * t401 + (t143 * t441 + t144 * t440) * pkin(2)) * m(4); t542 * t81 + t543 * t80 - t571 * t141 + (t107 + t59 - t539) * t274 + t383 + t291 * t72 - t360 * t73 + (-t247 - t363) * t273 + t363 * qJD(4) + t212 * t25 + t213 * t26 - t538 * t142 + t338 * t116 + t342 * t115 + (-g(1) * t340 + g(2) * t344) * t530 + (t13 * t542 + t14 * t543 + t2 * t213 + t212 * t3 + t274 * t94) * m(7) + (-t11 * t360 + t12 * t291 + t158 * t274 - t538 * t54 - t55 * t571) * m(6) + (t194 * t274 + t338 * t52 + t342 * t53 + t262 * (-t122 * t338 + t123 * t342)) * m(5) + (-t207 * t274 - t208 * t273 + t257) * m(4); (mrSges(6,1) * t323 + mrSges(6,2) * t325 - t370 + t371 - t524) * t461 - (mrSges(6,1) * t158 + Ifges(6,4) * t491 + Ifges(6,2) * t493 + Ifges(6,6) * t482 - t473 + t504) * t162 + t536 - (Ifges(7,4) * t497 + Ifges(7,2) * t500 + Ifges(7,6) * t484 + t511 - t529) * t101 + t521 - t107 * t469 - m(6) * (t158 * t469 + t54 * t60 + t55 * t61) + (Ifges(7,1) * t497 + Ifges(7,4) * t500 + Ifges(7,5) * t484 + t509 - t528) * t561 + t358 + (-mrSges(6,2) * t158 + Ifges(6,1) * t491 + Ifges(6,4) * t493 + Ifges(6,5) * t482 + t474 + t502) * t379 + (t11 * t334 + t12 * t333) * t516 + (Ifges(5,1) * t234 - t449) * t486 + (t243 * mrSges(6,1) - t244 * mrSges(6,2) - m(7) * (-t300 * t344 - t326 * t425) - t415 - mrSges(5,2) * t264 + t549 * t263) * g(2) + (-m(7) * (-t299 * t429 + t340 * t300) - t414 - t245 * mrSges(6,1) + t246 * mrSges(6,2) + mrSges(5,2) * t266 - t549 * t265) * g(1) + t550 * t81 + (t13 * t550 - t132 * t94 + t14 * t551 + t2 * t269 + t268 * t3) * m(7) + t551 * t80 + t268 * t25 + t269 * t26 - t194 * (mrSges(5,1) * t235 + mrSges(5,2) * t234) - t122 * t171 - t61 * t141 - t60 * t142 - t132 * t59 + t234 * t452 + t73 * t467 + t72 * t468 + (Ifges(5,5) * t234 - Ifges(5,6) * t235) * t482 + t148 * t485 + (-Ifges(5,2) * t235 + t149 + t225) * t487 + (t172 + t458) * t123; t498 * t326 * g(3) - t561 * t80 - t379 * t141 + t162 * t142 + t101 * t81 + t10 + t45 + (t101 * t13 - t14 * t561 + t43 - t546) * m(7) + (t162 * t54 - t379 * t55 - t546 + t79) * m(6); -t94 * (mrSges(7,1) * t101 + mrSges(7,2) * t561) + (Ifges(7,1) * t561 - t453) * t497 + t50 * t496 + (Ifges(7,5) * t561 - Ifges(7,6) * t101) * t484 - t13 * t80 + t14 * t81 - g(1) * t414 - g(2) * t415 - t370 * t461 + (t101 * t14 + t13 * t561) * mrSges(7,3) + t358 + (-Ifges(7,2) * t101 + t51 + t95) * t500;];
tau  = t1;
