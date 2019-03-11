% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR4
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:51
% EndTime: 2019-03-08 23:16:38
% DurationCPUTime: 32.04s
% Computational Cost: add. (13229->868), mult. (30230->1196), div. (0->0), fcn. (23077->18), ass. (0->381)
t332 = sin(qJ(3));
t336 = cos(qJ(3));
t378 = pkin(3) * t332 - pkin(9) * t336;
t276 = t378 * qJD(3);
t331 = sin(qJ(4));
t333 = sin(qJ(2));
t335 = cos(qJ(4));
t423 = qJD(3) * t332;
t410 = pkin(8) * t423;
t326 = sin(pkin(6));
t429 = qJD(1) * t326;
t337 = cos(qJ(2));
t436 = t336 * t337;
t574 = t335 * t276 + t331 * t410 - (-t331 * t436 + t333 * t335) * t429;
t284 = -pkin(3) * t336 - pkin(9) * t332 - pkin(2);
t418 = qJD(4) * t335;
t573 = t331 * t276 + t284 * t418 - (t331 * t333 + t335 * t436) * t429;
t437 = t335 * t336;
t309 = pkin(8) * t437;
t360 = pkin(4) * t332 - qJ(5) * t437;
t417 = qJD(5) * t335;
t572 = -t332 * t417 + t360 * qJD(3) + (-t309 + (qJ(5) * t332 - t284) * t331) * qJD(4) + t574;
t439 = t332 * t335;
t571 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t439 + (-qJD(5) * t332 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t336) * t331 + t573;
t405 = t333 * t429;
t280 = qJD(2) * pkin(8) + t405;
t328 = cos(pkin(6));
t428 = qJD(1) * t328;
t211 = -t280 * t332 + t336 * t428;
t275 = t378 * qJD(2);
t156 = -t211 * t331 + t275 * t335;
t329 = -qJ(5) - pkin(9);
t388 = qJD(4) * t329;
t570 = -qJD(2) * t360 - qJD(5) * t331 + t335 * t388 - t156;
t157 = t211 * t335 + t275 * t331;
t424 = qJD(2) * t336;
t401 = t331 * t424;
t569 = -qJ(5) * t401 - t331 * t388 + t157 - t417;
t325 = sin(pkin(12));
t327 = cos(pkin(12));
t265 = t325 * t335 + t327 * t331;
t355 = t265 * t336;
t215 = qJD(2) * t355;
t245 = t265 * qJD(4);
t568 = t215 - t245;
t362 = t325 * t331 - t327 * t335;
t354 = t362 * t336;
t216 = qJD(2) * t354;
t246 = t362 * qJD(4);
t567 = t216 - t246;
t377 = mrSges(4,1) * t336 - mrSges(4,2) * t332;
t517 = m(7) * (-pkin(10) + t329) - mrSges(7,3) + m(6) * t329 - mrSges(6,3) - m(5) * pkin(9) - mrSges(5,3);
t324 = qJ(4) + pkin(12);
t320 = cos(t324);
t464 = pkin(4) * t335;
t283 = pkin(5) * t320 + t464;
t321 = qJ(6) + t324;
t312 = sin(t321);
t313 = cos(t321);
t316 = pkin(3) + t464;
t319 = sin(t324);
t376 = -mrSges(5,1) * t335 + mrSges(5,2) * t331;
t557 = m(7) * (pkin(3) + t283) + mrSges(7,1) * t313 - mrSges(7,2) * t312 + m(6) * t316 + mrSges(6,1) * t320 - mrSges(6,2) * t319 + m(5) * pkin(3) - t376;
t566 = t332 * t517 - t336 * t557 - mrSges(3,1) - t377;
t538 = -t325 * t571 + t327 * t572;
t537 = t325 * t572 + t327 * t571;
t536 = t325 * t569 + t327 * t570;
t535 = t325 * t570 - t327 * t569;
t564 = -Ifges(6,3) - Ifges(5,3);
t426 = qJD(2) * t332;
t563 = -pkin(5) * t426 - pkin(10) * t567 + t536;
t562 = pkin(10) * t568 + t535;
t159 = -qJD(3) * t354 - t245 * t332;
t561 = pkin(5) * t423 - pkin(10) * t159 + t538;
t419 = qJD(4) * t332;
t158 = -qJD(3) * t355 + t362 * t419;
t560 = -pkin(10) * t158 - t537;
t519 = m(7) + m(6) + m(5) + m(4);
t465 = pkin(4) * t331;
t282 = pkin(5) * t319 + t465;
t375 = mrSges(5,1) * t331 + mrSges(5,2) * t335;
t556 = -m(6) * t465 - m(7) * t282 - mrSges(6,1) * t319 - mrSges(7,1) * t312 - mrSges(6,2) * t320 - mrSges(7,2) * t313 + mrSges(3,2) - mrSges(4,3) - t375;
t422 = qJD(3) * t335;
t271 = -t331 * t426 + t422;
t272 = qJD(3) * t331 + t335 * t426;
t181 = t271 * t325 + t272 * t327;
t555 = pkin(10) * t181;
t380 = t271 * t327 - t272 * t325;
t113 = -mrSges(6,1) * t380 + mrSges(6,2) * t181;
t330 = sin(qJ(6));
t334 = cos(qJ(6));
t105 = t181 * t334 + t330 * t380;
t551 = -t181 * t330 + t334 * t380;
t58 = -mrSges(7,1) * t551 + mrSges(7,2) * t105;
t554 = t113 + t58;
t212 = t280 * t336 + t332 * t428;
t420 = qJD(4) * t331;
t529 = -t212 + (-t401 + t420) * pkin(4);
t427 = qJD(2) * t326;
t395 = qJD(1) * t427;
t294 = t337 * t395;
t415 = qJDD(1) * t326;
t232 = t333 * t415 + t294;
t553 = qJDD(2) * pkin(8) + qJD(3) * t428 + t232;
t421 = qJD(3) * t336;
t347 = t331 * t421 + t332 * t418;
t552 = pkin(2) * t519 - t566;
t416 = qJD(2) * qJD(3);
t278 = qJDD(2) * t332 + t336 * t416;
t174 = qJD(4) * t271 + qJDD(3) * t331 + t278 * t335;
t175 = -qJD(4) * t272 + qJDD(3) * t335 - t278 * t331;
t93 = -t174 * t325 + t175 * t327;
t94 = t174 * t327 + t175 * t325;
t28 = qJD(6) * t551 + t330 * t93 + t334 * t94;
t504 = t28 / 0.2e1;
t29 = -qJD(6) * t105 - t330 * t94 + t334 * t93;
t503 = t29 / 0.2e1;
t491 = t93 / 0.2e1;
t490 = t94 / 0.2e1;
t485 = t174 / 0.2e1;
t484 = t175 / 0.2e1;
t277 = qJDD(2) * t336 - t332 * t416;
t263 = qJDD(4) - t277;
t258 = qJDD(6) + t263;
t479 = t258 / 0.2e1;
t478 = t263 / 0.2e1;
t549 = t277 / 0.2e1;
t548 = t278 / 0.2e1;
t267 = t335 * t284;
t182 = -qJ(5) * t439 + t267 + (-pkin(8) * t331 - pkin(4)) * t336;
t220 = t284 * t331 + t309;
t441 = t331 * t332;
t195 = -qJ(5) * t441 + t220;
t116 = t182 * t327 - t195 * t325;
t228 = t362 * t332;
t74 = -pkin(5) * t336 + pkin(10) * t228 + t116;
t117 = t182 * t325 + t195 * t327;
t227 = t265 * t332;
t79 = -pkin(10) * t227 + t117;
t40 = t330 * t74 + t334 * t79;
t547 = -qJD(6) * t40 + t330 * t560 + t334 * t561;
t39 = -t330 * t79 + t334 * t74;
t546 = qJD(6) * t39 + t330 * t561 - t334 * t560;
t545 = pkin(10) * t380;
t10 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t51 = -mrSges(6,1) * t93 + mrSges(6,2) * t94;
t543 = t10 + t51;
t287 = t329 * t331;
t288 = t329 * t335;
t198 = t287 * t327 + t288 * t325;
t163 = -pkin(10) * t265 + t198;
t199 = t287 * t325 - t288 * t327;
t164 = -pkin(10) * t362 + t199;
t75 = t163 * t334 - t164 * t330;
t542 = qJD(6) * t75 + t330 * t563 + t334 * t562;
t76 = t163 * t330 + t164 * t334;
t541 = -qJD(6) * t76 - t330 * t562 + t334 * t563;
t505 = m(6) * pkin(4);
t539 = -mrSges(5,1) - t505;
t466 = pkin(4) * t327;
t314 = pkin(5) + t466;
t467 = pkin(4) * t325;
t235 = t314 * t334 - t330 * t467;
t202 = qJD(3) * pkin(9) + t212;
t404 = t337 * t429;
t214 = qJD(2) * t284 - t404;
t136 = -t202 * t331 + t214 * t335;
t114 = -qJ(5) * t272 + t136;
t138 = t202 * t335 + t214 * t331;
t115 = qJ(5) * t271 + t138;
t442 = t327 * t115;
t59 = -t114 * t325 - t442;
t41 = t59 - t545;
t108 = t325 * t115;
t60 = t114 * t327 - t108;
t42 = t60 - t555;
t534 = qJD(6) * t235 - t330 * t41 - t334 * t42;
t236 = t314 * t330 + t334 * t467;
t533 = -qJD(6) * t236 + t330 * t42 - t334 * t41;
t100 = -mrSges(5,1) * t175 + mrSges(5,2) * t174;
t530 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t278 + t100;
t528 = -qJD(4) * t220 + t574;
t527 = (-t332 * t422 - t336 * t420) * pkin(8) + t573;
t526 = -pkin(5) * t568 + t529;
t262 = Ifges(5,4) * t271;
t307 = qJD(4) - t424;
t167 = Ifges(5,1) * t272 + Ifges(5,5) * t307 + t262;
t317 = Ifges(4,4) * t424;
t525 = Ifges(4,1) * t426 + Ifges(4,5) * qJD(3) + t335 * t167 + t317;
t408 = mrSges(4,3) * t426;
t524 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t271 + mrSges(5,2) * t272 + t408;
t523 = Ifges(5,5) * t174 + Ifges(6,5) * t94 + Ifges(5,6) * t175 + Ifges(6,6) * t93 - t263 * t564;
t414 = qJDD(1) * t328;
t132 = -t280 * t423 + t332 * t414 + t336 * t553;
t133 = -t280 * t421 - t332 * t553 + t336 * t414;
t522 = t132 * t336 - t133 * t332;
t126 = qJDD(3) * pkin(9) + t132;
t293 = t333 * t395;
t231 = t337 * t415 - t293;
t217 = -qJDD(2) * pkin(2) - t231;
t149 = -pkin(3) * t277 - pkin(9) * t278 + t217;
t46 = t126 * t335 + t149 * t331 - t202 * t420 + t214 * t418;
t47 = -qJD(4) * t138 - t126 * t331 + t149 * t335;
t521 = -t331 * t47 + t335 * t46;
t297 = qJD(6) + t307;
t520 = t272 * Ifges(5,5) + t181 * Ifges(6,5) + t105 * Ifges(7,5) + t271 * Ifges(5,6) + Ifges(6,6) * t380 + Ifges(7,6) * t551 + t297 * Ifges(7,3) - t307 * t564;
t516 = mrSges(4,1) + t557;
t515 = mrSges(4,2) + t517;
t92 = pkin(4) * t307 + t114;
t52 = t327 * t92 - t108;
t36 = pkin(5) * t307 + t52 - t555;
t53 = t325 * t92 + t442;
t38 = t53 + t545;
t15 = -t330 * t38 + t334 * t36;
t30 = pkin(4) * t263 - qJ(5) * t174 - qJD(5) * t272 + t47;
t32 = qJ(5) * t175 + qJD(5) * t271 + t46;
t11 = t30 * t327 - t32 * t325;
t6 = pkin(5) * t263 - pkin(10) * t94 + t11;
t12 = t30 * t325 + t32 * t327;
t7 = pkin(10) * t93 + t12;
t2 = qJD(6) * t15 + t330 * t6 + t334 * t7;
t16 = t330 * t36 + t334 * t38;
t3 = -qJD(6) * t16 - t330 * t7 + t334 * t6;
t513 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t201 = -qJD(3) * pkin(3) - t211;
t512 = -m(5) * t201 - t524;
t510 = -pkin(8) * t519 + t556;
t509 = -t47 * mrSges(5,1) - t11 * mrSges(6,1) + t46 * mrSges(5,2) + t12 * mrSges(6,2);
t459 = Ifges(4,4) * t332;
t370 = Ifges(4,2) * t336 + t459;
t508 = t16 * mrSges(7,2) + t53 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t370 / 0.2e1 - t15 * mrSges(7,1) - t52 * mrSges(6,1);
t338 = qJD(2) ^ 2;
t507 = Ifges(7,4) * t504 + Ifges(7,2) * t503 + Ifges(7,6) * t479;
t506 = Ifges(7,1) * t504 + Ifges(7,4) * t503 + Ifges(7,5) * t479;
t502 = Ifges(6,4) * t490 + Ifges(6,2) * t491 + Ifges(6,6) * t478;
t501 = Ifges(6,1) * t490 + Ifges(6,4) * t491 + Ifges(6,5) * t478;
t455 = Ifges(7,4) * t105;
t49 = Ifges(7,2) * t551 + Ifges(7,6) * t297 + t455;
t500 = -t49 / 0.2e1;
t499 = t49 / 0.2e1;
t97 = Ifges(7,4) * t551;
t50 = Ifges(7,1) * t105 + Ifges(7,5) * t297 + t97;
t498 = -t50 / 0.2e1;
t497 = t50 / 0.2e1;
t496 = Ifges(5,1) * t485 + Ifges(5,4) * t484 + Ifges(5,5) * t478;
t90 = Ifges(6,4) * t181 + Ifges(6,2) * t380 + Ifges(6,6) * t307;
t495 = -t90 / 0.2e1;
t494 = t90 / 0.2e1;
t91 = Ifges(6,1) * t181 + Ifges(6,4) * t380 + Ifges(6,5) * t307;
t493 = -t91 / 0.2e1;
t492 = t91 / 0.2e1;
t489 = -t551 / 0.2e1;
t488 = t551 / 0.2e1;
t487 = -t105 / 0.2e1;
t486 = t105 / 0.2e1;
t483 = -t380 / 0.2e1;
t482 = t380 / 0.2e1;
t481 = -t181 / 0.2e1;
t480 = t181 / 0.2e1;
t476 = t272 / 0.2e1;
t475 = -t297 / 0.2e1;
t474 = t297 / 0.2e1;
t473 = -t307 / 0.2e1;
t472 = t307 / 0.2e1;
t470 = mrSges(6,3) * t53;
t469 = mrSges(7,3) * t15;
t468 = pkin(4) * t272;
t463 = t16 * mrSges(7,3);
t322 = t332 * pkin(8);
t462 = t52 * mrSges(6,3);
t458 = Ifges(4,4) * t336;
t457 = Ifges(5,4) * t331;
t456 = Ifges(5,4) * t335;
t454 = t136 * mrSges(5,3);
t453 = t138 * mrSges(5,3);
t452 = t272 * Ifges(5,4);
t449 = cos(pkin(11));
t448 = sin(pkin(11));
t127 = -qJDD(3) * pkin(3) - t133;
t447 = t127 * t332;
t444 = t326 * t333;
t443 = t326 * t337;
t440 = t331 * t336;
t382 = t448 * t337;
t385 = t449 * t333;
t242 = t328 * t385 + t382;
t387 = t326 * t449;
t188 = t242 * t336 - t332 * t387;
t383 = t448 * t333;
t384 = t449 * t337;
t241 = -t328 * t384 + t383;
t435 = (-t188 * t312 + t241 * t313) * mrSges(7,1) + (-t188 * t313 - t241 * t312) * mrSges(7,2);
t244 = -t328 * t383 + t384;
t386 = t326 * t448;
t190 = t244 * t336 + t332 * t386;
t243 = t328 * t382 + t385;
t434 = (-t190 * t312 + t243 * t313) * mrSges(7,1) + (-t190 * t313 - t243 * t312) * mrSges(7,2);
t249 = t328 * t332 + t336 * t444;
t433 = (-t249 * t312 - t313 * t443) * mrSges(7,1) + (-t249 * t313 + t312 * t443) * mrSges(7,2);
t279 = pkin(4) * t441 + t322;
t425 = qJD(2) * t333;
t413 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t258;
t318 = pkin(8) * t421;
t407 = mrSges(4,3) * t424;
t213 = pkin(4) * t347 + t318;
t403 = t326 * t425;
t402 = t337 * t427;
t166 = t271 * Ifges(5,2) + Ifges(5,6) * t307 + t452;
t398 = -t331 * t166 / 0.2e1;
t381 = t416 / 0.2e1;
t372 = Ifges(5,1) * t335 - t457;
t371 = Ifges(5,1) * t331 + t456;
t369 = -Ifges(5,2) * t331 + t456;
t368 = Ifges(5,2) * t335 + t457;
t367 = Ifges(4,5) * t336 - Ifges(4,6) * t332;
t366 = Ifges(5,5) * t335 - Ifges(5,6) * t331;
t365 = Ifges(5,5) * t331 + Ifges(5,6) * t335;
t193 = -t249 * t331 - t335 * t443;
t359 = -t249 * t335 + t331 * t443;
t119 = t193 * t327 + t325 * t359;
t120 = t193 * t325 - t327 * t359;
t62 = t119 * t334 - t120 * t330;
t63 = t119 * t330 + t120 * t334;
t154 = -t227 * t334 + t228 * t330;
t155 = -t227 * t330 - t228 * t334;
t176 = -t265 * t330 - t334 * t362;
t177 = t265 * t334 - t330 * t362;
t361 = t413 - t513;
t248 = -t328 * t336 + t332 * t444;
t358 = t201 * t375;
t281 = -qJD(2) * pkin(2) - t404;
t357 = t281 * (mrSges(4,1) * t332 + mrSges(4,2) * t336);
t356 = t332 * (Ifges(4,1) * t336 - t459);
t187 = t242 * t332 + t336 * t387;
t189 = t244 * t332 - t336 * t386;
t351 = -g(1) * t189 - g(2) * t187 - g(3) * t248;
t348 = -t331 * t419 + t335 * t421;
t345 = Ifges(5,5) * t332 + t336 * t372;
t344 = Ifges(5,6) * t332 + t336 * t369;
t343 = Ifges(5,3) * t332 + t336 * t366;
t162 = -pkin(4) * t271 + qJD(5) + t201;
t71 = -pkin(4) * t175 + qJDD(5) + t127;
t286 = -qJD(3) * mrSges(4,2) + t407;
t274 = t377 * qJD(2);
t233 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t277;
t221 = pkin(5) * t362 - t316;
t219 = -pkin(8) * t440 + t267;
t210 = mrSges(5,1) * t307 - mrSges(5,3) * t272;
t209 = -mrSges(5,2) * t307 + mrSges(5,3) * t271;
t196 = -mrSges(4,1) * t277 + mrSges(4,2) * t278;
t192 = -qJD(3) * t248 + t336 * t402;
t191 = qJD(3) * t249 + t332 * t402;
t183 = pkin(5) * t227 + t279;
t153 = mrSges(6,1) * t307 - mrSges(6,3) * t181;
t152 = -mrSges(6,2) * t307 + mrSges(6,3) * t380;
t148 = -t215 * t330 - t216 * t334;
t147 = -t215 * t334 + t216 * t330;
t146 = pkin(5) * t181 + t468;
t142 = -mrSges(5,2) * t263 + mrSges(5,3) * t175;
t141 = mrSges(5,1) * t263 - mrSges(5,3) * t174;
t118 = -pkin(5) * t158 + t213;
t107 = qJD(4) * t193 + t192 * t335 + t331 * t403;
t106 = qJD(4) * t359 - t192 * t331 + t335 * t403;
t99 = -qJD(6) * t177 - t245 * t334 + t246 * t330;
t98 = qJD(6) * t176 - t245 * t330 - t246 * t334;
t96 = -pkin(5) * t380 + t162;
t81 = mrSges(7,1) * t297 - mrSges(7,3) * t105;
t80 = -mrSges(7,2) * t297 + mrSges(7,3) * t551;
t77 = Ifges(5,4) * t174 + Ifges(5,2) * t175 + Ifges(5,6) * t263;
t73 = mrSges(6,1) * t263 - mrSges(6,3) * t94;
t72 = -mrSges(6,2) * t263 + mrSges(6,3) * t93;
t66 = -qJD(6) * t155 + t158 * t334 - t159 * t330;
t65 = qJD(6) * t154 + t158 * t330 + t159 * t334;
t57 = t106 * t325 + t107 * t327;
t55 = t106 * t327 - t107 * t325;
t43 = -pkin(5) * t93 + t71;
t22 = -mrSges(7,2) * t258 + mrSges(7,3) * t29;
t21 = mrSges(7,1) * t258 - mrSges(7,3) * t28;
t14 = -qJD(6) * t63 - t330 * t57 + t334 * t55;
t13 = qJD(6) * t62 + t330 * t55 + t334 * t57;
t1 = [m(2) * qJDD(1) + t106 * t210 + t107 * t209 + t119 * t73 + t120 * t72 + t13 * t80 + t14 * t81 + t193 * t141 - t359 * t142 + t57 * t152 + t55 * t153 + t192 * t286 + t62 * t21 + t63 * t22 + t249 * t233 + (t530 + t543) * t248 + (t524 + t554) * t191 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t338 - t196) * t337 + (-mrSges(3,1) * t338 - mrSges(3,2) * qJDD(2) - qJD(2) * t274) * t333) * t326 + (-m(2) - m(3) - t519) * g(3) + m(7) * (t13 * t16 + t14 * t15 + t191 * t96 + t2 * t63 + t248 * t43 + t3 * t62) + m(6) * (t11 * t119 + t12 * t120 + t162 * t191 + t248 * t71 + t52 * t55 + t53 * t57) + m(5) * (t106 * t136 + t107 * t138 + t127 * t248 + t191 * t201 + t193 * t47 - t359 * t46) + m(4) * (t132 * t249 - t133 * t248 - t191 * t211 + t192 * t212 + (-t217 * t337 + t281 * t425) * t326) + m(3) * (qJDD(1) * t328 ^ 2 + (t231 * t337 + t232 * t333) * t326); (-Ifges(6,4) * t228 - Ifges(6,2) * t227) * t491 + (t11 * t228 - t12 * t227 + t158 * t53 - t159 * t52) * mrSges(6,3) + (-Ifges(6,1) * t228 - Ifges(6,4) * t227) * t490 + (-Ifges(6,5) * t228 - Ifges(6,6) * t227) * t478 + t71 * (mrSges(6,1) * t227 - mrSges(6,2) * t228) + (Ifges(6,4) * t159 + Ifges(6,2) * t158) * t482 + t356 * t381 + ((-Ifges(4,2) * t332 + t458) * t381 - t286 * t404 + Ifges(4,4) * t548 + Ifges(4,2) * t549 - Ifges(7,3) * t479 - Ifges(5,6) * t484 - Ifges(5,5) * t485 - Ifges(6,5) * t490 - Ifges(6,6) * t491 - Ifges(7,6) * t503 - Ifges(7,5) * t504 + pkin(8) * t233 + t564 * t478 + t509 + t513) * t336 + (t293 + t231) * mrSges(3,1) - (t335 * t166 + t331 * t167) * t419 / 0.2e1 + t458 * t548 + t370 * t549 + (Ifges(7,5) * t155 + Ifges(7,6) * t154) * t479 + (Ifges(7,5) * t65 + Ifges(7,6) * t66) * t474 + (Ifges(7,4) * t65 + Ifges(7,2) * t66) * t488 + (Ifges(7,4) * t155 + Ifges(7,2) * t154) * t503 + (t241 * t552 + t242 * t510) * g(2) + qJD(3) ^ 2 * t367 / 0.2e1 + t271 * (qJD(3) * t344 - t368 * t419) / 0.2e1 + t219 * t141 + t220 * t142 + t213 * t113 - pkin(2) * t196 + t274 * t405 + t398 * t421 + qJD(3) * t357 + t201 * (mrSges(5,1) * t347 + mrSges(5,2) * t348) + (t243 * t552 + t244 * t510) * g(1) + (-t519 * (pkin(2) * t443 + pkin(8) * t444) + (t556 * t333 + t337 * t566) * t326) * g(3) - t286 * t410 + (Ifges(6,1) * t159 + Ifges(6,4) * t158) * t480 + (Ifges(6,5) * t159 + Ifges(6,6) * t158 + qJD(3) * t343 - t365 * t419) * t472 + (-m(6) * t162 - m(7) * t96 + t512 - t554) * t332 * t404 + (Ifges(7,1) * t65 + Ifges(7,4) * t66) * t486 + (Ifges(7,1) * t155 + Ifges(7,4) * t154) * t504 + t375 * t447 + (-t136 * t348 - t138 * t347 - t439 * t47 - t441 * t46) * mrSges(5,3) - t77 * t441 / 0.2e1 + (t294 - t232) * mrSges(3,2) + t183 * t10 + t162 * (-mrSges(6,1) * t158 + mrSges(6,2) * t159) + t43 * (-mrSges(7,1) * t154 + mrSges(7,2) * t155) + t116 * t73 + t117 * t72 + t118 * t58 + Ifges(3,3) * qJDD(2) + t96 * (-mrSges(7,1) * t66 + mrSges(7,2) * t65) + t39 * t21 + t40 * t22 + (Ifges(4,1) * t278 + Ifges(4,4) * t549 + t366 * t478 + t369 * t484 + t372 * t485) * t332 + t546 * t80 + t547 * t81 + (t118 * t96 + t15 * t547 + t16 * t546 + t183 * t43 + t2 * t40 + t3 * t39) * m(7) + t537 * t152 + t538 * t153 + (t11 * t116 + t117 * t12 + t162 * t213 + t279 * t71 + t52 * t538 + t53 * t537) * m(6) + t524 * t318 + t525 * t421 / 0.2e1 + t527 * t209 + t528 * t210 + (t219 * t47 + t220 * t46 + (t201 * t421 + t447) * pkin(8) + t527 * t138 + t528 * t136) * m(5) + t530 * t322 + (-(t281 * t333 + (-t211 * t332 + t212 * t336) * t337) * t429 - pkin(2) * t217 + ((-t211 * t336 - t212 * t332) * qJD(3) + t522) * pkin(8)) * m(4) - (t413 + t523) * t336 / 0.2e1 - t217 * t377 + (qJD(3) * t345 - t371 * t419) * t476 + t279 * t51 + (-t212 * mrSges(4,3) + t520 / 0.2e1 + t136 * mrSges(5,1) - t138 * mrSges(5,2) + Ifges(6,5) * t480 + Ifges(7,5) * t486 + Ifges(6,6) * t482 + Ifges(7,6) * t488 + Ifges(6,3) * t472 + Ifges(7,3) * t474 - t508) * t423 + t159 * t492 + t158 * t494 + t439 * t496 + t65 * t497 + t66 * t499 - t228 * t501 - t227 * t502 + t155 * t506 + t154 * t507 + (-t15 * t65 + t154 * t2 - t155 * t3 + t16 * t66) * mrSges(7,3) + (-t211 * t421 + t522) * mrSges(4,3) + qJDD(3) * (Ifges(4,5) * t332 + Ifges(4,6) * t336); (t398 + t358) * qJD(4) + (-Ifges(6,5) * t216 - Ifges(6,6) * t215) * t473 + (-Ifges(6,4) * t216 - Ifges(6,2) * t215) * t483 + (-Ifges(6,1) * t216 - Ifges(6,4) * t215) * t481 + (-t11 * t265 - t12 * t362 + t215 * t53 - t216 * t52) * mrSges(6,3) + (Ifges(6,5) * t265 - Ifges(6,6) * t362 + t365) * t478 + t71 * (mrSges(6,1) * t362 + mrSges(6,2) * t265) + (Ifges(6,1) * t265 - Ifges(6,4) * t362) * t490 + (Ifges(6,4) * t265 - Ifges(6,2) * t362) * t491 - t362 * t502 + (t407 - t286) * t211 + (Ifges(7,4) * t486 + Ifges(7,2) * t488 + Ifges(7,6) * t474 + t463 + t499) * t99 + (Ifges(7,1) * t486 + Ifges(7,4) * t488 + Ifges(7,5) * t474 - t469 + t497) * t98 + (t271 * t369 + t272 * t372 + t307 * t366) * qJD(4) / 0.2e1 - (t271 * t344 + t272 * t345 + t307 * t343) * qJD(2) / 0.2e1 - (Ifges(6,1) * t480 + Ifges(6,4) * t482 + Ifges(6,5) * t472 - t462 + t492) * t246 - t358 * t424 - t367 * t416 / 0.2e1 + (-mrSges(6,1) * t568 + mrSges(6,2) * t567) * t162 + (Ifges(7,1) * t148 + Ifges(7,4) * t147) * t487 + t221 * t10 - t157 * t209 - t156 * t210 + t198 * t73 + t199 * t72 + (t167 / 0.2e1 - t454) * t418 + (Ifges(7,4) * t148 + Ifges(7,2) * t147) * t489 + (Ifges(7,5) * t148 + Ifges(7,6) * t147) * t475 - t420 * t453 - (Ifges(6,4) * t480 + Ifges(6,2) * t482 + Ifges(6,6) * t472 + t470 + t494) * t245 - t338 * t356 / 0.2e1 + ((-t148 + t98) * mrSges(7,2) + (t147 - t99) * mrSges(7,1)) * t96 + (-pkin(3) * t127 - t136 * t156 - t138 * t157) * m(5) + (-t357 - t138 * (-mrSges(5,2) * t332 - mrSges(5,3) * t440) - t136 * (mrSges(5,1) * t332 - mrSges(5,3) * t437)) * qJD(2) + t166 * t401 / 0.2e1 + t43 * (-mrSges(7,1) * t176 + mrSges(7,2) * t177) - t132 * mrSges(4,2) + t133 * mrSges(4,1) + Ifges(4,3) * qJDD(3) - pkin(3) * t100 + t75 * t21 + t76 * t22 + t541 * t81 + t542 * t80 + (t15 * t541 + t16 * t542 + t2 * t76 + t221 * t43 + t3 * t75 + t526 * t96) * m(7) + t535 * t152 + t536 * t153 + (t11 * t198 + t12 * t199 + t162 * t529 - t316 * t71 + t52 * t536 + t53 * t535) * m(6) - (-Ifges(4,2) * t426 + t317 + t525) * t424 / 0.2e1 + t526 * t58 + t529 * t113 + (m(5) * ((-t136 * t335 - t138 * t331) * qJD(4) + t521) - t209 * t420 - t210 * t418 - t331 * t141 + t335 * t142) * pkin(9) + t521 * mrSges(5,3) - t520 * t426 / 0.2e1 + (t189 * t516 + t190 * t515) * g(1) + (t187 * t516 + t188 * t515) * g(2) + (t248 * t516 + t249 * t515) * g(3) + (t408 + t512) * t212 + (Ifges(7,5) * t177 + Ifges(7,6) * t176) * t479 + t368 * t484 + t371 * t485 + Ifges(4,6) * t277 + Ifges(4,5) * t278 - t216 * t493 - t215 * t495 + t331 * t496 + t148 * t498 + t147 * t500 + t265 * t501 + (Ifges(7,4) * t177 + Ifges(7,2) * t176) * t503 + (Ifges(7,1) * t177 + Ifges(7,4) * t176) * t504 + t177 * t506 + t176 * t507 - t316 * t51 + t127 * t376 + t335 * t77 / 0.2e1 + (-t147 * t16 + t148 * t15 + t176 * t2 - t177 * t3) * mrSges(7,3) + (Ifges(6,5) * t481 + Ifges(7,5) * t487 + Ifges(6,6) * t483 + Ifges(7,6) * t489 + Ifges(6,3) * t473 + Ifges(7,3) * t475 + t508) * t426; (-mrSges(7,2) * t96 + Ifges(7,1) * t487 + Ifges(7,4) * t489 + Ifges(7,5) * t475 + t469 + t498) * t551 - (mrSges(6,1) * t162 + Ifges(6,4) * t481 + Ifges(6,2) * t483 + Ifges(6,6) * t473 - t470 + t495) * t181 + (-mrSges(5,2) * t359 - (-t249 * t319 - t320 * t443) * mrSges(6,1) - (-t249 * t320 + t319 * t443) * mrSges(6,2) - m(7) * (-t249 * t282 - t283 * t443) - t433 + t539 * t193) * g(3) - (mrSges(7,1) * t96 + Ifges(7,4) * t487 + Ifges(7,2) * t489 + Ifges(7,6) * t475 - t463 + t500) * t105 - (-Ifges(5,2) * t272 + t167 + t262) * t271 / 0.2e1 + t361 - t509 - t136 * t209 + t138 * t210 + t523 + (-mrSges(6,2) * t162 + Ifges(6,1) * t481 + Ifges(6,4) * t483 + Ifges(6,5) * t473 + t462 + t493) * t380 + t272 * t453 + t271 * t454 + t73 * t466 + t72 * t467 - t272 * (Ifges(5,1) * t271 - t452) / 0.2e1 - t113 * t468 - m(6) * (t162 * t468 + t52 * t59 + t53 * t60) - t60 * t152 - t59 * t153 - t146 * t58 + t533 * t81 + t534 * t80 + (-t146 * t96 + t533 * t15 + t534 * t16 + t2 * t236 + t235 * t3) * m(7) + (-m(7) * (-t190 * t282 + t243 * t283) - t434 - (-t190 * t319 + t243 * t320) * mrSges(6,1) - (-t190 * t320 - t243 * t319) * mrSges(6,2) - (-t190 * t335 - t243 * t331) * mrSges(5,2) + t539 * (-t190 * t331 + t243 * t335)) * g(1) + (-m(7) * (-t188 * t282 + t241 * t283) - t435 - (-t188 * t319 + t241 * t320) * mrSges(6,1) - (-t188 * t320 - t241 * t319) * mrSges(6,2) - (-t188 * t335 - t241 * t331) * mrSges(5,2) + t539 * (-t188 * t331 + t241 * t335)) * g(2) + t235 * t21 + t236 * t22 + (Ifges(5,5) * t271 - Ifges(5,6) * t272) * t473 + t166 * t476 - t201 * (mrSges(5,1) * t272 + mrSges(5,2) * t271) + (t11 * t327 + t12 * t325) * t505; t105 * t81 - t551 * t80 - t380 * t152 + t181 * t153 + (t105 * t15 - t16 * t551 + t351 + t43) * m(7) + (t181 * t52 - t380 * t53 + t351 + t71) * m(6) + t543; -t96 * (mrSges(7,1) * t105 + mrSges(7,2) * t551) + (Ifges(7,1) * t551 - t455) * t487 + t49 * t486 + (Ifges(7,5) * t551 - Ifges(7,6) * t105) * t475 - t15 * t80 + t16 * t81 - g(1) * t434 - g(2) * t435 - g(3) * t433 + (t105 * t16 + t15 * t551) * mrSges(7,3) + t361 + (-Ifges(7,2) * t105 + t50 + t97) * t489;];
tau  = t1;
