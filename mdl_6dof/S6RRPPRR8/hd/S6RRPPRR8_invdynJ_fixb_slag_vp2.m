% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:54
% EndTime: 2019-03-09 09:23:44
% DurationCPUTime: 33.42s
% Computational Cost: add. (11327->876), mult. (25163->1135), div. (0->0), fcn. (17835->12), ass. (0->399)
t586 = Ifges(5,1) + Ifges(4,1);
t312 = sin(pkin(10));
t313 = cos(pkin(10));
t316 = sin(qJ(5));
t320 = cos(qJ(5));
t251 = t312 * t320 - t313 * t316;
t321 = cos(qJ(2));
t344 = t251 * t321;
t206 = qJD(1) * t344;
t224 = t251 * qJD(5);
t549 = t206 + t224;
t441 = t312 * t316;
t353 = t313 * t320 + t441;
t343 = t353 * t321;
t207 = qJD(1) * t343;
t223 = t353 * qJD(5);
t548 = t207 + t223;
t298 = t321 * qJDD(1);
t317 = sin(qJ(2));
t411 = qJD(1) * qJD(2);
t397 = t317 * t411;
t257 = -t298 + t397;
t489 = t257 / 0.2e1;
t258 = qJDD(1) * t317 + t321 * t411;
t209 = qJDD(2) * t312 + t258 * t313;
t492 = t209 / 0.2e1;
t584 = Ifges(4,5) + Ifges(5,4);
t590 = t489 * t584 + t586 * t492;
t564 = mrSges(6,3) + mrSges(7,3);
t404 = -pkin(7) * t312 - pkin(3);
t437 = t313 * t321;
t334 = -pkin(8) * t437 + (-pkin(4) + t404) * t317;
t365 = pkin(2) * t317 - qJ(3) * t321;
t256 = t365 * qJD(1);
t442 = t256 * t313;
t120 = qJD(1) * t334 - t442;
t229 = t312 * t256;
t424 = qJD(1) * t317;
t287 = qJ(4) * t424;
t438 = t313 * t317;
t439 = t312 * t321;
t345 = -pkin(7) * t438 + pkin(8) * t439;
t142 = qJD(1) * t345 + t229 + t287;
t470 = pkin(8) - qJ(3);
t262 = t470 * t312;
t263 = t470 * t313;
t413 = qJD(5) * t320;
t414 = qJD(5) * t316;
t417 = qJD(3) * t320;
t419 = qJD(3) * t316;
t554 = -t316 * t120 - t320 * t142 - t262 * t413 + t263 * t414 + t312 * t419 + t313 * t417;
t183 = -t316 * t262 - t320 * t263;
t553 = -qJD(5) * t183 - t320 * t120 + t142 * t316 + t312 * t417 - t313 * t419;
t299 = t321 * qJD(1);
t283 = t299 + qJD(5);
t270 = qJD(6) + t283;
t486 = -t283 / 0.2e1;
t488 = -t270 / 0.2e1;
t401 = t312 * t424;
t412 = t313 * qJD(2);
t244 = t401 - t412;
t399 = t313 * t424;
t245 = qJD(2) * t312 + t399;
t356 = t244 * t316 + t245 * t320;
t496 = -t356 / 0.2e1;
t163 = t244 * t320 - t245 * t316;
t498 = -t163 / 0.2e1;
t315 = sin(qJ(6));
t319 = cos(qJ(6));
t572 = t163 * t315 + t319 * t356;
t501 = -t572 / 0.2e1;
t388 = t319 * t163 - t356 * t315;
t503 = -t388 / 0.2e1;
t588 = Ifges(6,5) * t496 + Ifges(7,5) * t501 + Ifges(6,6) * t498 + Ifges(7,6) * t503 + Ifges(6,3) * t486 + Ifges(7,3) * t488;
t208 = -t313 * qJDD(2) + t258 * t312;
t493 = t208 / 0.2e1;
t587 = mrSges(5,2) + mrSges(4,3);
t585 = -Ifges(4,4) + Ifges(5,5);
t562 = Ifges(4,6) - Ifges(5,6);
t561 = -Ifges(4,3) - Ifges(5,2);
t583 = pkin(5) * t424 + pkin(9) * t548 + t553;
t582 = pkin(9) * t549 - t554;
t461 = Ifges(5,5) * t312;
t464 = Ifges(4,4) * t312;
t580 = t313 * t586 + t461 - t464;
t579 = m(6) * pkin(8) + t564 - t587;
t310 = qJ(5) + qJ(6);
t300 = sin(t310);
t301 = cos(t310);
t354 = t300 * t312 + t301 * t313;
t355 = t300 * t313 - t301 * t312;
t376 = t313 * mrSges(5,1) + t312 * mrSges(5,3);
t378 = mrSges(4,1) * t313 - mrSges(4,2) * t312;
t573 = t353 * mrSges(6,1) + mrSges(6,2) * t251;
t578 = t354 * mrSges(7,1) - t355 * mrSges(7,2) + t376 + t378 + t573;
t569 = -t257 / 0.2e1;
t577 = m(4) * qJD(3);
t291 = Ifges(3,4) * t299;
t576 = Ifges(3,1) * t424 + Ifges(3,5) * qJD(2) + t291 + t313 * (t244 * t585 + t245 * t586 - t299 * t584);
t400 = t312 * t299;
t416 = qJD(4) * t312;
t293 = pkin(7) * t299;
t402 = qJ(4) * t299;
t428 = t313 * t402 - t293;
t499 = -pkin(3) - pkin(4);
t543 = -t400 * t499 + t416 - t428;
t303 = t313 * pkin(4);
t535 = -t313 * pkin(3) - t312 * qJ(4) - pkin(2);
t575 = t535 - t303;
t418 = qJD(3) * t317;
t449 = qJDD(1) * pkin(1);
t152 = pkin(2) * t257 - qJ(3) * t258 - qJD(1) * t418 - t449;
t290 = pkin(7) * t298;
t292 = pkin(7) * t424;
t213 = qJDD(2) * qJ(3) + t290 + (qJD(3) - t292) * qJD(2);
t90 = t152 * t313 - t312 * t213;
t363 = qJDD(4) - t90;
t76 = -pkin(3) * t257 + t363;
t574 = t76 * mrSges(5,2) + t493 * t585 + t590;
t302 = t317 * qJ(3);
t307 = t321 * pkin(2);
t405 = -pkin(1) - t307;
t351 = t405 - t302;
t235 = t351 * qJD(1);
t265 = qJD(2) * qJ(3) + t293;
t170 = t312 * t235 + t313 * t265;
t145 = -t402 + t170;
t103 = pkin(8) * t244 + t145;
t169 = t235 * t313 - t312 * t265;
t141 = pkin(3) * t299 + qJD(4) - t169;
t96 = pkin(4) * t299 - pkin(8) * t245 + t141;
t46 = -t103 * t316 + t320 * t96;
t41 = -pkin(9) * t356 + t46;
t37 = pkin(5) * t283 + t41;
t47 = t103 * t320 + t316 * t96;
t42 = pkin(9) * t163 + t47;
t456 = t315 * t42;
t13 = t319 * t37 - t456;
t453 = t319 * t42;
t14 = t315 * t37 + t453;
t466 = Ifges(3,4) * t317;
t558 = t321 * Ifges(3,2);
t371 = t466 + t558;
t571 = -t46 * mrSges(6,1) - t13 * mrSges(7,1) + t47 * mrSges(6,2) + t14 * mrSges(7,2) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t371 / 0.2e1 + t561 * t299 / 0.2e1 + t584 * t245 / 0.2e1 - t562 * t244 / 0.2e1 + t588;
t516 = m(7) * pkin(5);
t69 = qJD(5) * t163 + t208 * t316 + t209 * t320;
t70 = -qJD(5) * t356 + t208 * t320 - t209 * t316;
t21 = qJD(6) * t388 + t315 * t70 + t319 * t69;
t515 = t21 / 0.2e1;
t22 = -qJD(6) * t572 - t315 * t69 + t319 * t70;
t514 = t22 / 0.2e1;
t507 = t69 / 0.2e1;
t506 = t70 / 0.2e1;
t570 = -m(6) - m(7);
t249 = qJDD(5) - t257;
t238 = qJDD(6) + t249;
t491 = t238 / 0.2e1;
t490 = t249 / 0.2e1;
t243 = t258 * pkin(7);
t567 = qJD(2) / 0.2e1;
t566 = mrSges(2,2) - mrSges(3,3);
t565 = -mrSges(4,2) + mrSges(5,3);
t182 = -t320 * t262 + t263 * t316;
t130 = -pkin(9) * t251 + t182;
t131 = -pkin(9) * t353 + t183;
t60 = t130 * t319 - t131 * t315;
t560 = qJD(6) * t60 + t315 * t583 - t319 * t582;
t61 = t130 * t315 + t131 * t319;
t559 = -qJD(6) * t61 + t315 * t582 + t319 * t583;
t380 = mrSges(3,1) * t321 - mrSges(3,2) * t317;
t557 = mrSges(2,1) + t380;
t556 = -mrSges(6,1) - t516;
t426 = t307 + t302;
t261 = -pkin(1) - t426;
t277 = pkin(7) * t439;
t306 = t321 * pkin(3);
t475 = pkin(8) * t317;
t147 = pkin(4) * t321 + t277 + t306 + (-t261 - t475) * t313;
t201 = pkin(7) * t437 + t312 * t261;
t450 = qJ(4) * t321;
t185 = t201 - t450;
t440 = t312 * t317;
t168 = pkin(8) * t440 + t185;
t79 = t316 * t147 + t320 * t168;
t550 = pkin(5) * t549 + t543;
t255 = t315 * t320 + t316 * t319;
t547 = t270 * t255;
t352 = t315 * t316 - t319 * t320;
t546 = t270 * t352;
t545 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t244 - mrSges(4,2) * t245 - mrSges(3,3) * t424;
t423 = qJD(2) * t317;
t544 = qJ(4) * t423 - qJD(4) * t321;
t242 = -pkin(7) * t397 + t290;
t542 = t242 * t321 + t243 * t317;
t541 = t587 * t317;
t91 = t312 * t152 + t313 * t213;
t540 = -t312 * t90 + t313 * t91;
t318 = sin(qJ(1));
t322 = cos(qJ(1));
t539 = g(1) * t322 + g(2) * t318;
t537 = -m(5) + t570;
t533 = m(4) - t537;
t460 = Ifges(5,5) * t313;
t366 = Ifges(5,3) * t312 + t460;
t531 = t244 * (Ifges(5,6) * t317 + t321 * t366) + (t317 * t584 + t321 * t580) * t245;
t323 = -pkin(9) - pkin(8);
t289 = pkin(5) * t320 + pkin(4);
t348 = pkin(5) * t441 + t289 * t313;
t528 = (-m(7) * t323 + t579) * t321 + (-m(7) * (-t348 + t535) - m(6) * t575 - m(5) * t535 + m(4) * pkin(2) + t578) * t317;
t51 = -pkin(8) * t209 + t257 * t499 + t363;
t71 = t257 * qJ(4) - qJD(4) * t299 + t91;
t54 = pkin(8) * t208 + t71;
t11 = -t103 * t414 + t316 * t51 + t320 * t54 + t96 * t413;
t10 = pkin(9) * t70 + t11;
t12 = -qJD(5) * t47 - t316 * t54 + t320 * t51;
t9 = pkin(5) * t249 - pkin(9) * t69 + t12;
t2 = qJD(6) * t13 + t10 * t319 + t315 * t9;
t3 = -qJD(6) * t14 - t10 * t315 + t319 * t9;
t527 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t526 = qJ(3) * t533;
t525 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t523 = Ifges(5,5) * t492 + Ifges(5,6) * t489 - Ifges(4,4) * t209 / 0.2e1 + Ifges(4,6) * t569 - t71 * mrSges(5,2) + (Ifges(5,3) + Ifges(4,2)) * t493;
t522 = m(6) * pkin(4) + m(7) * t289 + mrSges(4,1) + mrSges(5,1);
t260 = -qJD(2) * pkin(2) + qJD(3) + t292;
t135 = t244 * pkin(3) - t245 * qJ(4) + t260;
t463 = Ifges(4,4) * t313;
t370 = -Ifges(4,2) * t312 + t463;
t375 = mrSges(5,1) * t312 - mrSges(5,3) * t313;
t377 = mrSges(4,1) * t312 + mrSges(4,2) * t313;
t521 = -t260 * t321 * t377 - t135 * t321 * t375 - t170 * (-mrSges(4,2) * t317 - mrSges(4,3) * t439) - t169 * (mrSges(4,1) * t317 - mrSges(4,3) * t437) - t145 * (-mrSges(5,2) * t439 + mrSges(5,3) * t317) - t141 * (-mrSges(5,1) * t317 + mrSges(5,2) * t437) + t244 * (Ifges(4,6) * t317 + t321 * t370) / 0.2e1;
t518 = Ifges(7,4) * t515 + Ifges(7,2) * t514 + Ifges(7,6) * t491;
t517 = Ifges(7,1) * t515 + Ifges(7,4) * t514 + Ifges(7,5) * t491;
t513 = Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * t490;
t512 = Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * t490;
t479 = Ifges(7,4) * t572;
t39 = Ifges(7,2) * t388 + Ifges(7,6) * t270 + t479;
t511 = -t39 / 0.2e1;
t510 = t39 / 0.2e1;
t80 = Ifges(7,4) * t388;
t40 = Ifges(7,1) * t572 + Ifges(7,5) * t270 + t80;
t509 = -t40 / 0.2e1;
t508 = t40 / 0.2e1;
t462 = Ifges(6,4) * t356;
t74 = t163 * Ifges(6,2) + t283 * Ifges(6,6) + t462;
t505 = t74 / 0.2e1;
t160 = Ifges(6,4) * t163;
t75 = Ifges(6,1) * t356 + t283 * Ifges(6,5) + t160;
t504 = t75 / 0.2e1;
t502 = t388 / 0.2e1;
t500 = t572 / 0.2e1;
t497 = t163 / 0.2e1;
t495 = t356 / 0.2e1;
t494 = -t208 / 0.2e1;
t487 = t270 / 0.2e1;
t485 = t283 / 0.2e1;
t481 = mrSges(7,3) * t13;
t480 = mrSges(7,3) * t14;
t478 = pkin(5) * t356;
t216 = t251 * t317;
t477 = pkin(5) * t216;
t476 = pkin(7) * t317;
t468 = mrSges(6,3) * t163;
t467 = mrSges(6,3) * t356;
t465 = Ifges(3,4) * t321;
t218 = -qJDD(2) * pkin(2) + qJDD(3) + t243;
t447 = t218 * t317;
t221 = qJD(2) * t365 - t418;
t446 = t221 * t313;
t434 = t318 * t321;
t225 = t312 * t434 + t313 * t322;
t445 = t225 * t316;
t436 = t317 * t322;
t435 = t317 * t323;
t433 = t321 * t322;
t432 = t322 * t312;
t226 = t313 * t434 - t432;
t359 = -t225 * t300 - t226 * t301;
t360 = t225 * t301 - t226 * t300;
t431 = t360 * mrSges(7,1) + t359 * mrSges(7,2);
t227 = -t318 * t313 + t321 * t432;
t228 = t312 * t318 + t313 * t433;
t127 = t227 * t301 - t228 * t300;
t128 = t227 * t300 + t228 * t301;
t430 = t127 * mrSges(7,1) - t128 * mrSges(7,2);
t429 = (-mrSges(7,1) * t355 - mrSges(7,2) * t354) * t317;
t427 = -qJD(4) * t438 - t412 * t450;
t425 = t322 * pkin(1) + t318 * pkin(7);
t422 = qJD(2) * t321;
t421 = qJD(3) * t312;
t420 = qJD(3) * t313;
t410 = Ifges(7,5) * t21 + Ifges(7,6) * t22 + Ifges(7,3) * t238;
t409 = Ifges(6,5) * t69 + Ifges(6,6) * t70 + Ifges(6,3) * t249;
t408 = pkin(7) * t423;
t403 = pkin(3) * t312 + pkin(7);
t391 = -t299 / 0.2e1;
t390 = -t411 / 0.2e1;
t389 = t411 / 0.2e1;
t115 = t208 * mrSges(4,1) + t209 * mrSges(4,2);
t151 = -t257 * mrSges(5,1) + t209 * mrSges(5,2);
t114 = t208 * mrSges(5,1) - t209 * mrSges(5,3);
t78 = t320 * t147 - t168 * t316;
t200 = t261 * t313 - t277;
t387 = pkin(3) * t437 + qJ(4) * t439 + t426;
t386 = pkin(2) * t433 + qJ(3) * t436 + t425;
t383 = t312 * t499 - pkin(7);
t382 = t404 * t317;
t36 = -t70 * mrSges(6,1) + t69 * mrSges(6,2);
t8 = -t22 * mrSges(7,1) + t21 * mrSges(7,2);
t379 = mrSges(3,1) * t317 + mrSges(3,2) * t321;
t217 = t353 * t317;
t374 = mrSges(6,1) * t216 - mrSges(6,2) * t217;
t369 = Ifges(5,4) * t313 + Ifges(5,6) * t312;
t368 = Ifges(3,5) * t321 - Ifges(3,6) * t317;
t367 = Ifges(4,5) * t313 - Ifges(4,6) * t312;
t58 = pkin(5) * t321 - pkin(9) * t217 + t78;
t59 = pkin(9) * t216 + t79;
t29 = -t315 * t59 + t319 * t58;
t30 = t315 * t58 + t319 * t59;
t362 = t228 * pkin(3) + t386;
t118 = -mrSges(6,2) * t283 + t468;
t119 = mrSges(6,1) * t283 - t467;
t361 = t118 * t320 - t119 * t316;
t125 = t216 * t319 - t217 * t315;
t126 = t216 * t315 + t217 * t319;
t358 = t225 * t320 - t226 * t316;
t357 = -t226 * t320 - t445;
t143 = t227 * t320 - t228 * t316;
t166 = -t251 * t315 - t319 * t353;
t167 = t251 * t319 - t315 * t353;
t188 = -pkin(7) * t399 + t229;
t212 = t312 * t221;
t181 = -t313 * t408 + t212;
t349 = t410 + t527;
t347 = pkin(1) * t379;
t346 = t317 * (Ifges(3,1) * t321 - t466);
t109 = qJD(2) * t334 - t446;
t110 = qJD(2) * t345 + t212 + t544;
t34 = t316 * t109 + t320 * t110 + t147 * t413 - t168 * t414;
t275 = qJ(4) * t438;
t184 = t317 * t383 + t275;
t101 = -pkin(4) * t244 - t135;
t335 = -g(1) * t227 - g(2) * t225 - g(3) * t440;
t140 = t383 * t422 - t427;
t329 = t321 * (Ifges(5,2) * t317 + t321 * t369);
t328 = t321 * (Ifges(4,3) * t317 + t321 * t367);
t35 = -qJD(5) * t79 + t320 * t109 - t110 * t316;
t77 = t208 * pkin(3) - t209 * qJ(4) - t245 * qJD(4) + t218;
t57 = pkin(4) * t208 + t77;
t308 = t322 * pkin(7);
t266 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t299;
t214 = t317 * t403 - t275;
t205 = mrSges(5,1) * t299 + mrSges(5,2) * t245;
t204 = -mrSges(4,1) * t299 - mrSges(4,3) * t245;
t203 = mrSges(4,2) * t299 - mrSges(4,3) * t244;
t202 = -mrSges(5,2) * t244 - mrSges(5,3) * t299;
t193 = pkin(3) * t400 - t428;
t187 = pkin(7) * t401 + t442;
t186 = -t200 + t306;
t180 = t312 * t408 + t446;
t177 = pkin(5) * t353 - t575;
t176 = qJD(1) * t382 - t442;
t174 = mrSges(5,1) * t244 - mrSges(5,3) * t245;
t173 = t188 + t287;
t172 = t403 * t422 + t427;
t161 = qJD(2) * t382 - t446;
t156 = t245 * Ifges(4,4) - t244 * Ifges(4,2) - Ifges(4,6) * t299;
t153 = t245 * Ifges(5,5) - Ifges(5,6) * t299 + t244 * Ifges(5,3);
t150 = mrSges(4,1) * t257 - mrSges(4,3) * t209;
t149 = -mrSges(4,2) * t257 - mrSges(4,3) * t208;
t148 = -mrSges(5,2) * t208 + mrSges(5,3) * t257;
t144 = t227 * t316 + t228 * t320;
t138 = t181 + t544;
t134 = qJD(2) * t343 + qJD(5) * t216;
t133 = qJD(2) * t344 - t223 * t317;
t113 = t184 - t477;
t108 = t206 * t315 + t207 * t319;
t107 = t206 * t319 - t207 * t315;
t89 = -mrSges(6,1) * t163 + mrSges(6,2) * t356;
t88 = -qJD(6) * t167 + t223 * t315 - t224 * t319;
t87 = qJD(6) * t166 - t223 * t319 - t224 * t315;
t72 = -pkin(5) * t133 + t140;
t66 = mrSges(7,1) * t270 - mrSges(7,3) * t572;
t65 = -mrSges(7,2) * t270 + mrSges(7,3) * t388;
t62 = -pkin(5) * t163 + t101;
t56 = -mrSges(6,2) * t249 + mrSges(6,3) * t70;
t55 = mrSges(6,1) * t249 - mrSges(6,3) * t69;
t45 = -qJD(6) * t126 + t133 * t319 - t134 * t315;
t44 = qJD(6) * t125 + t133 * t315 + t134 * t319;
t43 = -mrSges(7,1) * t388 + mrSges(7,2) * t572;
t33 = pkin(5) * t70 + t57;
t26 = pkin(9) * t133 + t34;
t25 = -pkin(5) * t423 - pkin(9) * t134 + t35;
t18 = -mrSges(7,2) * t238 + mrSges(7,3) * t22;
t17 = mrSges(7,1) * t238 - mrSges(7,3) * t21;
t16 = t319 * t41 - t456;
t15 = -t315 * t41 - t453;
t5 = -qJD(6) * t30 + t25 * t319 - t26 * t315;
t4 = qJD(6) * t29 + t25 * t315 + t26 * t319;
t1 = [(-t91 * mrSges(4,3) + t523) * t440 + (Ifges(7,5) * t44 + Ifges(7,6) * t45) * t487 + (Ifges(7,5) * t126 + Ifges(7,6) * t125) * t491 + (Ifges(6,5) * t134 + Ifges(6,6) * t133) * t485 + (Ifges(6,5) * t217 + Ifges(6,6) * t216) * t490 + (Ifges(7,4) * t44 + Ifges(7,2) * t45) * t502 + (Ifges(7,4) * t126 + Ifges(7,2) * t125) * t514 + (Ifges(6,4) * t134 + Ifges(6,2) * t133) * t497 + (Ifges(6,4) * t217 + Ifges(6,2) * t216) * t506 + (t11 * t216 - t12 * t217 + t133 * t47 - t134 * t46) * mrSges(6,3) + (Ifges(7,1) * t44 + Ifges(7,4) * t45) * t500 + (Ifges(7,1) * t126 + Ifges(7,4) * t125) * t515 + (Ifges(6,1) * t134 + Ifges(6,4) * t133) * t495 + (Ifges(6,1) * t217 + Ifges(6,4) * t216) * t507 + qJDD(2) * Ifges(3,6) * t321 + (t465 * t389 - Ifges(5,6) * t493 - Ifges(4,6) * t494 + Ifges(6,3) * t490 + Ifges(7,3) * t491 + Ifges(6,6) * t506 + Ifges(6,5) * t507 + Ifges(7,6) * t514 + Ifges(7,5) * t515 - t90 * mrSges(4,1) + t76 * mrSges(5,1) + t91 * mrSges(4,2) - t71 * mrSges(5,3) + pkin(7) * (-qJDD(2) * mrSges(3,2) - mrSges(3,3) * t257) - t584 * t492 + t561 * t489 + t525 + t527) * t321 - (-t208 * t562 + t209 * t584 - t257 * t561) * t321 / 0.2e1 + t258 * t465 / 0.2e1 + (-m(3) * t425 - m(4) * t386 - m(7) * t362 - t144 * mrSges(6,1) - t128 * mrSges(7,1) - t143 * mrSges(6,2) - t127 * mrSges(7,2) + (-m(6) - m(5)) * (t227 * qJ(4) + t362) + (-m(7) * t435 - t557) * t322 + t566 * t318 - t522 * t228 + (-m(7) * (pkin(5) * t316 + qJ(4)) - t565) * t227 + t579 * t436) * g(2) + (t312 * t153 + t576) * t422 / 0.2e1 + (-mrSges(4,3) * t90 + t574) * t438 + (t329 + t328) * t390 + (-qJDD(2) * mrSges(3,1) + t115) * t476 + (Ifges(3,4) * t258 - Ifges(3,2) * t257 + t409 + t410) * t321 / 0.2e1 + t531 * t567 + t371 * t569 - t312 * t156 * t422 / 0.2e1 + (t125 * t2 - t126 * t3 - t13 * t44 + t14 * t45) * mrSges(7,3) + (-Ifges(6,5) * t495 - Ifges(7,5) * t500 - Ifges(6,6) * t497 - Ifges(7,6) * t502 - Ifges(6,3) * t485 - Ifges(7,3) * t487 + t571) * t423 - t545 * pkin(7) * t422 + t380 * t449 + t377 * t447 + t346 * t389 + t57 * t374 + Ifges(2,3) * qJDD(1) - pkin(1) * (mrSges(3,1) * t257 + mrSges(3,2) * t258) + t214 * t114 + t138 * t202 + t181 * t203 + t180 * t204 + t161 * t205 + t200 * t150 + t201 * t149 + t184 * t36 + t185 * t148 + t186 * t151 + t172 * t174 + t140 * t89 + t101 * (-mrSges(6,1) * t133 + mrSges(6,2) * t134) + t34 * t118 + t35 * t119 - t33 * (-mrSges(7,1) * t125 + mrSges(7,2) * t126) + t113 * t8 + t72 * t43 + t78 * t55 + t79 * t56 + t4 * t65 + t5 * t66 + t62 * (-mrSges(7,1) * t45 + mrSges(7,2) * t44) + t29 * t17 + t30 * t18 + (t445 * t516 - t357 * mrSges(6,1) - t359 * mrSges(7,1) + t358 * mrSges(6,2) + t360 * mrSges(7,2) + t537 * (-t226 * pkin(3) - t225 * qJ(4) + t308) + t566 * t322 + (-m(3) - m(4)) * t308 + t522 * t226 + t565 * t225 + (m(3) * pkin(1) + t570 * t405 + (-m(4) - m(5)) * t351 + (-m(6) * t470 - m(7) * (-qJ(3) - t323) - t564) * t317 + t541 + t557) * t318) * g(1) + (t368 * t567 - t521) * qJD(2) + (t77 * t375 - t389 * t558 + Ifges(3,5) * qJDD(2) + Ifges(3,1) * t258 + Ifges(3,4) * t569 + t366 * t493 + t370 * t494 + t580 * t492 + (t367 + t369) * t489) * t317 - t266 * t408 + m(7) * (-t113 * t33 + t13 * t5 + t14 * t4 + t2 * t30 + t29 * t3 + t62 * t72) + m(6) * (t101 * t140 + t11 * t79 + t12 * t78 - t184 * t57 + t34 * t47 + t35 * t46) + m(5) * (t135 * t172 + t138 * t145 + t141 * t161 + t185 * t71 + t186 * t76 + t214 * t77) - t347 * t411 + t134 * t504 + t133 * t505 + t44 * t508 + t45 * t510 + t217 * t512 + t216 * t513 + t126 * t517 + t125 * t518 + (t258 * t476 + t542) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t542 * pkin(7)) + m(4) * (t169 * t180 + t170 * t181 + t200 * t90 + t201 * t91 + (t260 * t422 + t447) * pkin(7)); (t153 * t391 + (t151 - t150) * qJ(3) + m(5) * (qJ(3) * t76 + qJD(3) * t141 - qJD(4) * t135) - t169 * t577 + t574 + t590) * t312 + (t420 - t173) * t202 + (Ifges(7,4) * t108 + Ifges(7,2) * t107) * t503 + (Ifges(7,4) * t500 + Ifges(7,2) * t502 + Ifges(7,6) * t487 + t480 + t510) * t88 + (Ifges(7,1) * t500 + Ifges(7,4) * t502 + Ifges(7,5) * t487 - t481 + t508) * t87 + (Ifges(6,4) * t207 + Ifges(6,2) * t206) * t498 + (t521 - t531 / 0.2e1 + (-t346 / 0.2e1 + t329 / 0.2e1 + t328 / 0.2e1 + t347) * qJD(1)) * qJD(1) + (-pkin(2) * t218 + t540 * qJ(3) - t169 * t187 - t170 * t188 - t260 * t293) * m(4) + (-t460 + t463) * t492 - t57 * t573 + (t291 + t576) * t391 + t156 * t400 / 0.2e1 + (-t380 - m(6) * (t387 - t475) - m(7) * (t387 + t435) - m(4) * t426 - m(5) * t387 + t564 * t317 + (-m(6) * t303 - m(7) * t348 - t578) * t321 - t541) * g(3) + ((t149 + t148) * qJ(3) + m(5) * (qJ(3) * t71 + qJD(3) * t145) - Ifges(5,3) * t493 + Ifges(4,2) * t494 + t562 * t489 - t523 + t170 * t577) * t313 + (t101 * t543 + t11 * t183 + t12 * t182 + t46 * t553 + t47 * t554 + t57 * t575) * m(6) - t575 * t36 + (-t187 - t421) * t204 + (t421 - t176) * t205 + (Ifges(7,1) * t108 + Ifges(7,4) * t107) * t501 + (Ifges(6,5) * t251 - Ifges(6,6) * t353) * t490 + (Ifges(6,4) * t251 - Ifges(6,2) * t353) * t506 + (Ifges(6,1) * t251 - Ifges(6,4) * t353) * t507 - t353 * t513 + (-t11 * t353 - t12 * t251 + t46 * t548 - t47 * t549) * mrSges(6,3) + (-t135 * t193 - t141 * t176 - t145 * t173 + t535 * t77) * m(5) + t535 * t114 + (t318 * t528 - t434 * t526) * g(2) + (t322 * t528 - t433 * t526) * g(1) + (t420 - t188) * t203 + (-t107 * t14 + t108 * t13 + t166 * t2 - t167 * t3) * mrSges(7,3) + t266 * t292 + t368 * t390 + (-t193 - t416) * t174 - t77 * t376 - t218 * t378 + (Ifges(6,5) * t207 + Ifges(6,6) * t206) * t486 + (Ifges(7,5) * t108 + Ifges(7,6) * t107) * t488 + (-Ifges(3,2) * t391 - t571 - t588) * t424 + (Ifges(6,1) * t207 + Ifges(6,4) * t206) * t496 + t461 * t493 + t464 * t494 - Ifges(3,6) * t257 + Ifges(3,5) * t258 - t242 * mrSges(3,2) - t243 * mrSges(3,1) - t207 * t75 / 0.2e1 - t206 * t74 / 0.2e1 + t182 * t55 + t183 * t56 + t177 * t8 - t33 * (-mrSges(7,1) * t166 + mrSges(7,2) * t167) + Ifges(3,3) * qJDD(2) - pkin(2) * t115 + t60 * t17 + t61 * t18 + t559 * t66 + t560 * t65 + (t13 * t559 + t14 * t560 - t177 * t33 + t2 * t61 + t3 * t60 + t550 * t62) * m(7) + (-Ifges(6,5) * t223 - Ifges(6,6) * t224) * t485 + (-Ifges(6,1) * t223 - Ifges(6,4) * t224) * t495 + (-Ifges(6,4) * t223 - Ifges(6,2) * t224) * t497 + t539 * t379 + t540 * mrSges(4,3) + (Ifges(7,5) * t167 + Ifges(7,6) * t166) * t491 - t223 * t504 - t224 * t505 + t108 * t509 + t107 * t511 + t251 * t512 + (Ifges(7,4) * t167 + Ifges(7,2) * t166) * t514 + (Ifges(7,1) * t167 + Ifges(7,4) * t166) * t515 + t167 * t517 + t166 * t518 + t543 * t89 + t545 * t293 + (mrSges(6,1) * t549 - mrSges(6,2) * t548) * t101 + t550 * t43 + t553 * t119 + t554 * t118 + ((-t108 + t87) * mrSges(7,2) + (t107 - t88) * mrSges(7,1)) * t62; -t8 - t36 + t114 + (-t205 + t204) * t245 + (t202 + t203) * t244 + t163 * t118 - t356 * t119 - t572 * t66 + t388 * t65 + t115 + (-t13 * t572 + t14 * t388 + t33) * m(7) + (t163 * t47 - t356 * t46 + t57) * m(6) + (-t141 * t245 + t145 * t244 + t77) * m(5) + (t169 * t245 + t170 * t244 + t218) * m(4) + (g(3) * t321 - t317 * t539) * t533; -t352 * t17 + t255 * t18 + t316 * t56 + t320 * t55 - t547 * t66 - t546 * t65 + t361 * qJD(5) + (-t89 - t43 + t174) * t245 + (t202 + t361) * t299 + t151 + (-t13 * t547 - t14 * t546 + t2 * t255 - t245 * t62 - t3 * t352 + t335) * m(7) + (-t101 * t245 + t11 * t316 + t12 * t320 + t335 - t283 * (t316 * t46 - t320 * t47)) * m(6) + (t135 * t245 + t145 * t299 + t335 + t76) * m(5); (mrSges(6,2) * t144 + t143 * t556 - t430) * g(1) + (-mrSges(6,2) * t357 + t358 * t556 - t431) * g(2) + (-m(7) * t477 - t374 - t429) * g(3) - t101 * (mrSges(6,1) * t356 + mrSges(6,2) * t163) + (Ifges(6,5) * t163 - Ifges(6,6) * t356) * t486 + (-Ifges(6,2) * t356 + t160 + t75) * t498 - (mrSges(7,1) * t62 + Ifges(7,4) * t501 + Ifges(7,2) * t503 + Ifges(7,6) * t488 - t480 + t511) * t572 + (-mrSges(7,2) * t62 + Ifges(7,1) * t501 + Ifges(7,4) * t503 + Ifges(7,5) * t488 + t481 + t509) * t388 + t349 + (t468 - t118) * t46 + t409 + t525 + (t467 + t119) * t47 - t16 * t65 - t15 * t66 - t43 * t478 - m(7) * (t13 * t15 + t14 * t16 + t478 * t62) + t74 * t495 + (Ifges(6,1) * t163 - t462) * t496 + (t2 * t315 + t3 * t319 + (-t13 * t315 + t14 * t319) * qJD(6)) * t516 + (t17 * t319 + t18 * t315 + (-t315 * t66 + t319 * t65) * qJD(6)) * pkin(5); -t62 * (mrSges(7,1) * t572 + mrSges(7,2) * t388) + (Ifges(7,1) * t388 - t479) * t501 + t39 * t500 + (Ifges(7,5) * t388 - Ifges(7,6) * t572) * t488 - t13 * t65 + t14 * t66 - g(1) * t430 - g(2) * t431 - g(3) * t429 + (t13 * t388 + t14 * t572) * mrSges(7,3) + t349 + (-Ifges(7,2) * t572 + t40 + t80) * t503;];
tau  = t1;
