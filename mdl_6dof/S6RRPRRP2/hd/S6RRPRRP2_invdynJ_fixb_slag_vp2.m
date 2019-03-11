% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:06
% EndTime: 2019-03-09 11:43:46
% DurationCPUTime: 24.92s
% Computational Cost: add. (16546->764), mult. (38449->972), div. (0->0), fcn. (29107->14), ass. (0->352)
t582 = -mrSges(6,3) - mrSges(7,2);
t577 = Ifges(6,1) + Ifges(7,1);
t585 = -Ifges(6,4) + Ifges(7,5);
t556 = Ifges(6,5) + Ifges(7,4);
t576 = Ifges(6,6) - Ifges(7,6);
t315 = cos(qJ(5));
t407 = qJD(5) * t315;
t308 = sin(pkin(10));
t309 = cos(pkin(10));
t313 = sin(qJ(2));
t317 = cos(qJ(2));
t243 = -t308 * t313 + t309 * t317;
t227 = t243 * qJD(1);
t415 = qJD(1) * t317;
t416 = qJD(1) * t313;
t228 = -t308 * t415 - t309 * t416;
t312 = sin(qJ(4));
t316 = cos(qJ(4));
t370 = t316 * t227 + t228 * t312;
t540 = t370 * t315;
t584 = t407 - t540;
t311 = sin(qJ(5));
t408 = qJD(5) * t311;
t583 = t370 * t311 - t408;
t406 = qJD(1) * qJD(2);
t381 = t313 * t406;
t405 = qJDD(1) * t317;
t253 = -t381 + t405;
t254 = qJDD(1) * t313 + t317 * t406;
t192 = t253 * t309 - t254 * t308;
t193 = t253 * t308 + t254 * t309;
t107 = qJD(4) * t370 + t192 * t312 + t193 * t316;
t305 = qJDD(2) + qJDD(4);
t341 = t227 * t312 - t316 * t228;
t404 = qJD(2) + qJD(4);
t157 = t311 * t341 - t315 * t404;
t409 = qJD(5) * t157;
t73 = t315 * t107 + t311 * t305 - t409;
t511 = t73 / 0.2e1;
t158 = t311 * t404 + t315 * t341;
t74 = qJD(5) * t158 + t311 * t107 - t315 * t305;
t509 = t74 / 0.2e1;
t108 = -qJD(4) * t341 + t192 * t316 - t193 * t312;
t106 = qJDD(5) - t108;
t507 = t106 / 0.2e1;
t559 = mrSges(6,1) + mrSges(7,1);
t558 = mrSges(6,2) - mrSges(7,3);
t554 = Ifges(6,3) + Ifges(7,2);
t169 = qJD(5) - t370;
t471 = mrSges(7,2) * t157;
t120 = mrSges(7,3) * t169 - t471;
t469 = mrSges(6,3) * t157;
t121 = -mrSges(6,2) * t169 - t469;
t581 = -t121 - t120;
t468 = mrSges(6,3) * t158;
t122 = mrSges(6,1) * t169 - t468;
t470 = mrSges(7,2) * t158;
t123 = -mrSges(7,1) * t169 + t470;
t580 = t123 - t122;
t349 = pkin(5) * t315 + qJ(6) * t311;
t265 = -pkin(4) - t349;
t307 = qJ(2) + pkin(10);
t303 = qJ(4) + t307;
t289 = sin(t303);
t359 = -t315 * mrSges(7,1) - t311 * mrSges(7,3);
t433 = t289 * t311;
t472 = mrSges(6,1) * t315;
t579 = -mrSges(6,2) * t433 + (-m(7) * t265 - t359 + t472) * t289;
t563 = m(6) + m(7);
t578 = m(5) + t563;
t575 = t106 * t556 + t577 * t73 + t585 * t74;
t155 = Ifges(6,4) * t157;
t459 = t157 * Ifges(7,5);
t550 = t158 * t577 + t556 * t169 - t155 + t459;
t310 = -qJ(3) - pkin(7);
t273 = t310 * t313;
t251 = qJD(1) * t273;
t275 = t310 * t317;
t252 = qJD(1) * t275;
t425 = t309 * t252;
t190 = -t251 * t308 + t425;
t487 = pkin(8) * t227;
t159 = t190 - t487;
t231 = t308 * t252;
t191 = t309 * t251 + t231;
t486 = pkin(8) * t228;
t160 = t191 + t486;
t490 = pkin(2) * t309;
t291 = pkin(3) + t490;
t491 = pkin(2) * t308;
t222 = t312 * t291 + t316 * t491;
t533 = t222 * qJD(4) + t316 * t159 - t160 * t312;
t348 = pkin(5) * t311 - qJ(6) * t315;
t574 = pkin(5) * t408 - qJ(6) * t407 - qJD(6) * t311 - t348 * t370;
t274 = -t317 * mrSges(3,1) + t313 * mrSges(3,2);
t301 = sin(t307);
t302 = cos(t307);
t573 = -t302 * mrSges(4,1) + t301 * mrSges(4,2) + t274;
t572 = t582 * t289;
t358 = t311 * mrSges(7,1) - t315 * mrSges(7,3);
t360 = mrSges(6,1) * t311 + mrSges(6,2) * t315;
t238 = qJD(2) * pkin(2) + t251;
t184 = t309 * t238 + t231;
t149 = qJD(2) * pkin(3) + t184 + t486;
t185 = t308 * t238 - t425;
t153 = t185 + t487;
t95 = t316 * t149 - t312 * t153;
t89 = -pkin(4) * t404 - t95;
t48 = t157 * pkin(5) - t158 * qJ(6) + t89;
t571 = t48 * t358 + t89 * t360;
t530 = -t311 * t576 + t315 * t556;
t461 = Ifges(7,5) * t311;
t463 = Ifges(6,4) * t311;
t529 = t315 * t577 + t461 - t463;
t304 = t317 * pkin(2);
t294 = t304 + pkin(1);
t256 = -qJD(1) * t294 + qJD(3);
t194 = -pkin(3) * t227 + t256;
t111 = -pkin(4) * t370 - pkin(9) * t341 + t194;
t242 = t254 * pkin(7);
t413 = qJD(3) * t313;
t182 = qJDD(2) * pkin(2) - qJ(3) * t254 - qJD(1) * t413 - t242;
t295 = pkin(7) * t405;
t414 = qJD(2) * t313;
t398 = pkin(7) * t414;
t412 = qJD(3) * t317;
t188 = qJ(3) * t253 + t295 + (-t398 + t412) * qJD(1);
t134 = t309 * t182 - t188 * t308;
t102 = qJDD(2) * pkin(3) - pkin(8) * t193 + t134;
t135 = t308 * t182 + t309 * t188;
t110 = pkin(8) * t192 + t135;
t410 = qJD(4) * t316;
t411 = qJD(4) * t312;
t30 = t312 * t102 + t316 * t110 + t149 * t410 - t153 * t411;
t25 = pkin(9) * t305 + t30;
t450 = qJDD(1) * pkin(1);
t213 = -pkin(2) * t253 + qJDD(3) - t450;
t156 = -pkin(3) * t192 + t213;
t40 = -pkin(4) * t108 - pkin(9) * t107 + t156;
t96 = t312 * t149 + t316 * t153;
t90 = pkin(9) * t404 + t96;
t6 = t111 * t407 + t315 * t25 + t311 * t40 - t408 * t90;
t46 = t111 * t311 + t315 * t90;
t7 = -qJD(5) * t46 - t25 * t311 + t315 * t40;
t569 = -t311 * t7 + t315 * t6;
t2 = qJ(6) * t106 + qJD(6) * t169 + t6;
t4 = -pkin(5) * t106 + qJDD(6) - t7;
t568 = t2 * t315 + t311 * t4;
t314 = sin(qJ(1));
t318 = cos(qJ(1));
t567 = g(1) * t318 + g(2) * t314;
t566 = Ifges(7,5) * t511 + Ifges(7,6) * t507 - Ifges(6,4) * t73 / 0.2e1 - Ifges(6,6) * t106 / 0.2e1 + (Ifges(7,3) + Ifges(6,2)) * t509;
t45 = t111 * t315 - t311 * t90;
t37 = -pkin(5) * t169 + qJD(6) - t45;
t38 = qJ(6) * t169 + t46;
t33 = mrSges(6,1) * t106 - mrSges(6,3) * t73;
t34 = -t106 * mrSges(7,1) + t73 * mrSges(7,2);
t552 = -t33 + t34;
t32 = -mrSges(7,2) * t74 + mrSges(7,3) * t106;
t35 = -mrSges(6,2) * t106 - mrSges(6,3) * t74;
t553 = t32 + t35;
t565 = m(6) * ((-t311 * t46 - t315 * t45) * qJD(5) + t569) + m(7) * ((-t311 * t38 + t315 * t37) * qJD(5) + t568) + t581 * t408 + t580 * t407 + t315 * t553 + t311 * t552;
t31 = t102 * t316 - t312 * t110 - t149 * t411 - t153 * t410;
t26 = -pkin(4) * t305 - t31;
t460 = Ifges(7,5) * t315;
t350 = Ifges(7,3) * t311 + t460;
t462 = Ifges(6,4) * t315;
t354 = -Ifges(6,2) * t311 + t462;
t376 = t407 / 0.2e1;
t377 = -t408 / 0.2e1;
t154 = Ifges(7,5) * t158;
t80 = t169 * Ifges(7,6) + t157 * Ifges(7,3) + t154;
t452 = t311 * t80;
t394 = t452 / 0.2e1;
t290 = cos(t303);
t428 = t290 * t318;
t430 = t290 * t314;
t492 = t311 / 0.2e1;
t510 = -t74 / 0.2e1;
t458 = t158 * Ifges(6,4);
t83 = -t157 * Ifges(6,2) + t169 * Ifges(6,6) + t458;
t9 = pkin(5) * t74 - qJ(6) * t73 - qJD(6) * t158 + t26;
t564 = (t394 + t571) * qJD(5) + t461 * t509 + t463 * t510 + (-t354 / 0.2e1 + t350 / 0.2e1) * t409 - t26 * t472 + t575 * t492 + (Ifges(6,2) * t510 - Ifges(7,3) * t509 + t507 * t576 - t566) * t315 + (t26 * mrSges(6,2) + t507 * t556 + t511 * t577) * t311 + Ifges(5,3) * t305 + (t158 * t529 + t169 * t530) * qJD(5) / 0.2e1 + t9 * t359 + Ifges(5,5) * t107 + Ifges(5,6) * t108 - t30 * mrSges(5,2) + t31 * mrSges(5,1) + (-t460 + t462) * t511 + t83 * t377 + (t318 * t579 + t428 * t582) * g(1) + (t314 * t579 + t430 * t582) * g(2) + (-t540 / 0.2e1 + t376) * t550 + (-t45 * t584 + t46 * t583 + t569) * mrSges(6,3) + (t37 * t584 + t38 * t583 + t568) * mrSges(7,2);
t562 = t253 / 0.2e1;
t561 = t254 / 0.2e1;
t497 = -t341 / 0.2e1;
t560 = pkin(5) * t341;
t474 = qJD(2) / 0.2e1;
t551 = -t157 * t576 + t158 * t556 + t554 * t169;
t244 = t308 * t317 + t309 * t313;
t549 = Ifges(4,5) * t244;
t548 = Ifges(5,5) * t404;
t547 = Ifges(4,6) * t243;
t546 = Ifges(5,6) * t404;
t451 = qJDD(2) / 0.2e1;
t545 = t574 + t533;
t544 = -t96 + t574;
t543 = qJ(6) * t341;
t187 = t243 * t312 + t244 * t316;
t202 = -pkin(3) * t243 - t294;
t340 = t316 * t243 - t244 * t312;
t126 = -pkin(4) * t340 - pkin(9) * t187 + t202;
t196 = t309 * t273 + t275 * t308;
t166 = -pkin(8) * t244 + t196;
t197 = t308 * t273 - t309 * t275;
t167 = pkin(8) * t243 + t197;
t128 = t166 * t312 + t167 * t316;
t537 = t311 * t126 + t315 * t128;
t456 = t341 * mrSges(5,3);
t536 = mrSges(5,1) * t404 - mrSges(6,1) * t157 - mrSges(6,2) * t158 - t456;
t535 = t316 * t166 - t167 * t312;
t114 = t159 * t312 + t160 * t316;
t221 = t291 * t316 - t312 * t491;
t206 = t221 * qJD(4);
t534 = t206 - t114;
t429 = t290 * t315;
t431 = t290 * t311;
t532 = pkin(5) * t429 + qJ(6) * t431;
t531 = t290 * mrSges(5,1) - t289 * mrSges(5,2);
t418 = t290 * pkin(4) + t289 * pkin(9);
t229 = t244 * qJD(2);
t230 = t243 * qJD(2);
t136 = qJD(4) * t340 - t229 * t312 + t230 * t316;
t333 = t136 * t311 + t187 * t407;
t528 = t106 * t554 + t556 * t73 - t576 * t74;
t241 = -pkin(7) * t381 + t295;
t527 = t241 * t317 + t242 * t313;
t525 = -t429 * t559 + t431 * t558 - t531 + t572;
t522 = -m(3) * pkin(7) + m(4) * t310 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t375 = qJD(2) * t310;
t224 = t313 * t375 + t412;
t225 = t317 * t375 - t413;
t164 = -t224 * t308 + t309 * t225;
t145 = -pkin(8) * t230 + t164;
t165 = t309 * t224 + t308 * t225;
t146 = -pkin(8) * t229 + t165;
t54 = qJD(4) * t535 + t145 * t312 + t146 * t316;
t137 = qJD(4) * t187 + t316 * t229 + t230 * t312;
t299 = pkin(2) * t414;
t199 = pkin(3) * t229 + t299;
t67 = pkin(4) * t137 - pkin(9) * t136 + t199;
t13 = -qJD(5) * t537 - t311 * t54 + t315 * t67;
t521 = m(7) * pkin(5) + t559;
t520 = -m(3) * pkin(1) - m(4) * t294 - mrSges(2,1) - t531 + t573;
t519 = m(7) * qJ(6) - t558;
t518 = -m(6) * t89 + t536;
t133 = pkin(4) * t341 - pkin(9) * t370;
t517 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t501 = -t169 / 0.2e1;
t503 = -t158 / 0.2e1;
t504 = t157 / 0.2e1;
t505 = -t157 / 0.2e1;
t514 = -t452 / 0.2e1 + t83 * t492 - t548 / 0.2e1 - t194 * mrSges(5,2) + Ifges(5,1) * t497 + t354 * t504 + t350 * t505 + t529 * t503 + t530 * t501 - t571;
t498 = -t370 / 0.2e1;
t513 = -t45 * mrSges(6,1) + t37 * mrSges(7,1) + t46 * mrSges(6,2) - t38 * mrSges(7,3) + t546 / 0.2e1 - t194 * mrSges(5,1) - Ifges(5,2) * t498 + Ifges(6,6) * t504 + Ifges(7,6) * t505 + t556 * t503 + t554 * t501;
t502 = t158 / 0.2e1;
t496 = t341 / 0.2e1;
t493 = -t228 / 0.2e1;
t489 = pkin(2) * t313;
t488 = pkin(7) * t317;
t481 = g(3) * t289;
t476 = t95 * mrSges(5,3);
t475 = t96 * mrSges(5,3);
t467 = Ifges(3,4) * t313;
t466 = Ifges(3,4) * t317;
t465 = Ifges(4,4) * t228;
t464 = Ifges(5,4) * t341;
t457 = t370 * mrSges(5,3);
t455 = t184 * mrSges(4,3);
t454 = t228 * mrSges(4,3);
t57 = t311 * t133 + t315 * t95;
t448 = t136 * t315;
t440 = t206 * t311;
t439 = t206 * t315;
t432 = t289 * t318;
t424 = t314 * t311;
t423 = t314 * t315;
t422 = t315 * t318;
t421 = t318 * t311;
t297 = pkin(2) * t416;
t198 = -pkin(3) * t228 + t297;
t115 = t133 + t198;
t53 = t315 * t114 + t311 * t115;
t417 = pkin(3) * t302 + t304;
t401 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t416) * t488;
t373 = -t192 * mrSges(4,1) + t193 * mrSges(4,2);
t372 = -t108 * mrSges(5,1) + t107 * mrSges(5,2);
t250 = pkin(1) + t417;
t306 = -pkin(8) + t310;
t369 = t318 * t250 - t306 * t314;
t368 = t417 + t418;
t266 = pkin(9) * t430;
t366 = -pkin(4) * t289 * t314 + t266;
t269 = pkin(9) * t428;
t365 = -pkin(4) * t432 + t269;
t364 = mrSges(3,1) * t313 + mrSges(3,2) * t317;
t361 = mrSges(5,1) * t289 + mrSges(5,2) * t290;
t355 = Ifges(3,2) * t317 + t467;
t352 = Ifges(3,5) * t317 - Ifges(3,6) * t313;
t346 = t311 * t37 + t315 * t38;
t345 = -t311 * t45 + t315 * t46;
t56 = t133 * t315 - t311 * t95;
t52 = -t114 * t311 + t115 * t315;
t60 = t126 * t315 - t128 * t311;
t336 = pkin(1) * t364;
t332 = t187 * t408 - t448;
t331 = t313 * (Ifges(3,1) * t317 - t467);
t12 = t126 * t407 - t128 * t408 + t311 * t67 + t315 * t54;
t55 = qJD(4) * t128 - t316 * t145 + t146 * t312;
t296 = Ifges(3,4) * t415;
t272 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t415;
t255 = -pkin(3) * t301 - t489;
t240 = t318 * t255;
t239 = t314 * t255;
t237 = Ifges(3,1) * t416 + Ifges(3,5) * qJD(2) + t296;
t236 = Ifges(3,6) * qJD(2) + qJD(1) * t355;
t220 = Ifges(4,4) * t227;
t218 = -pkin(4) - t221;
t217 = t290 * t422 + t424;
t216 = t290 * t421 - t423;
t215 = t290 * t423 - t421;
t214 = t290 * t424 + t422;
t201 = qJD(2) * mrSges(4,1) + t454;
t200 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t227;
t195 = -t221 + t265;
t183 = -mrSges(4,1) * t227 - mrSges(4,2) * t228;
t179 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t193;
t178 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t192;
t173 = -t228 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t220;
t172 = t227 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t465;
t168 = Ifges(5,4) * t370;
t161 = -mrSges(5,2) * t404 + t457;
t132 = -mrSges(5,1) * t370 + mrSges(5,2) * t341;
t130 = Ifges(5,1) * t341 + t168 + t548;
t129 = Ifges(5,2) * t370 + t464 + t546;
t117 = mrSges(7,1) * t157 - mrSges(7,3) * t158;
t116 = pkin(5) * t158 + qJ(6) * t157;
t92 = -mrSges(5,2) * t305 + mrSges(5,3) * t108;
t91 = mrSges(5,1) * t305 - mrSges(5,3) * t107;
t75 = t187 * t348 - t535;
t50 = pkin(5) * t340 - t60;
t49 = -qJ(6) * t340 + t537;
t44 = -t56 - t560;
t43 = t57 + t543;
t42 = -t52 - t560;
t41 = t53 + t543;
t28 = mrSges(6,1) * t74 + mrSges(6,2) * t73;
t27 = mrSges(7,1) * t74 - mrSges(7,3) * t73;
t14 = t348 * t136 + (qJD(5) * t349 - qJD(6) * t315) * t187 + t55;
t11 = -pkin(5) * t137 - t13;
t10 = qJ(6) * t137 - qJD(6) * t340 + t12;
t1 = [t256 * (mrSges(4,1) * t229 + mrSges(4,2) * t230) + t227 * (Ifges(4,4) * t230 - Ifges(4,2) * t229) / 0.2e1 + (Ifges(4,5) * t230 - Ifges(4,6) * t229) * t474 + (Ifges(4,1) * t230 - Ifges(4,4) * t229) * t493 + (-t134 * t244 + t135 * t243 - t185 * t229) * mrSges(4,3) + (t331 + t317 * (-Ifges(3,2) * t313 + t466)) * t406 / 0.2e1 + t370 * (Ifges(5,4) * t136 - Ifges(5,2) * t137) / 0.2e1 + (t521 * t215 + t519 * t214 + (m(5) * t250 - t563 * (-t250 - t418) - t520 - t572) * t314 + (t306 * t578 + t522) * t318) * g(1) + m(6) * (t12 * t46 + t13 * t45 + t537 * t6 + t60 * t7) + t537 * t35 + ((mrSges(7,2) * t4 - mrSges(6,3) * t7 + t575 / 0.2e1) * t315 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t566) * t311 + t550 * t377 + t156 * mrSges(5,2) - t31 * mrSges(5,3) + Ifges(5,1) * t107 + Ifges(5,4) * t108 + Ifges(5,5) * t305 + t26 * t360 + t350 * t509 + t354 * t510 + t358 * t9 + t376 * t80 + t507 * t530 + t511 * t529) * t187 + m(5) * (t128 * t30 + t156 * t202 + t194 * t199 + t54 * t96) + (t137 * t554 - t332 * t556 - t333 * t576) * t169 / 0.2e1 + (t352 * t474 - t401) * qJD(2) - (-m(5) * t31 + m(6) * t26 + t28 - t91) * t535 + t136 * t394 - t236 * t414 / 0.2e1 + t550 * t448 / 0.2e1 + t551 * t137 / 0.2e1 + (Ifges(3,4) * t561 + Ifges(3,2) * t562 + Ifges(3,6) * t451 + t237 * t474) * t317 + (Ifges(3,1) * t254 + Ifges(3,4) * t562 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t254) + 0.2e1 * Ifges(3,5) * t451) * t313 + (-m(5) * t95 - t518) * t55 - t336 * t406 - t230 * t455 + t404 * (Ifges(5,5) * t136 - Ifges(5,6) * t137) / 0.2e1 + (t253 * t488 + t527) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t527) - t333 * t83 / 0.2e1 - t274 * t450 + (t547 / 0.2e1 + t549 / 0.2e1 + Ifges(3,6) * t317 / 0.2e1 - mrSges(3,2) * t488) * qJDD(2) + (t547 + t549) * t451 - t137 * t475 - t136 * t476 - pkin(1) * (-mrSges(3,1) * t253 + mrSges(3,2) * t254) + (-t528 / 0.2e1 - t156 * mrSges(5,1) + t30 * mrSges(5,3) + Ifges(5,4) * t107 + Ifges(5,2) * t108 + Ifges(5,6) * t305 - Ifges(6,6) * t510 - Ifges(7,6) * t509 - t507 * t554 - t511 * t556 - t517) * t340 + t192 * (Ifges(4,4) * t244 + Ifges(4,2) * t243) + t213 * (-mrSges(4,1) * t243 + mrSges(4,2) * t244) + t193 * (Ifges(4,1) * t244 + Ifges(4,4) * t243) - t229 * t172 / 0.2e1 + t230 * t173 / 0.2e1 + t194 * (mrSges(5,1) * t137 + mrSges(5,2) * t136) + t196 * t179 + t197 * t178 + t199 * t132 + t165 * t200 + t164 * t201 + Ifges(2,3) * qJDD(1) + t54 * t161 - t137 * t129 / 0.2e1 + t136 * t130 / 0.2e1 + t10 * t120 + t12 * t121 + t13 * t122 + t11 * t123 + t128 * t92 + t14 * t117 + t37 * (-mrSges(7,1) * t137 - mrSges(7,2) * t332) + t48 * (mrSges(7,1) * t333 + mrSges(7,3) * t332) + t89 * (mrSges(6,1) * t333 - mrSges(6,2) * t332) + t38 * (-mrSges(7,2) * t333 + mrSges(7,3) * t137) + t46 * (-mrSges(6,2) * t137 - mrSges(6,3) * t333) + t45 * (mrSges(6,1) * t137 + mrSges(6,3) * t332) + t75 * t27 + t60 * t33 + t49 * t32 + t50 * t34 + t202 * t372 - t294 * t373 + (t556 * t137 - t577 * t332 + t333 * t585) * t502 + t183 * t299 + m(4) * (t134 * t196 + t135 * t197 + t164 * t184 + t165 * t185 - t213 * t294 + t256 * t299) - t272 * t398 + m(7) * (t10 * t38 + t11 * t37 + t14 * t48 + t2 * t49 + t4 * t50 + t75 * t9) + t466 * t561 + t355 * t562 + (Ifges(5,1) * t136 - Ifges(5,4) * t137) * t496 + (-Ifges(7,5) * t332 + Ifges(7,6) * t137 + Ifges(7,3) * t333) * t504 + (-Ifges(6,4) * t332 - Ifges(6,2) * t333 + Ifges(6,6) * t137) * t505 + (-m(5) * t369 + t582 * t432 - t563 * (pkin(4) * t428 + pkin(9) * t432 + t369) - t521 * t217 - t519 * t216 + t520 * t318 + t522 * t314) * g(2); -(-Ifges(3,2) * t416 + t237 + t296) * t415 / 0.2e1 - (Ifges(4,2) * t228 + t173 + t220) * t227 / 0.2e1 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t567 * (m(4) * t489 - m(5) * t255 + mrSges(4,1) * t301 + mrSges(4,2) * t302 + t361 + t364) + (-m(7) * (t368 + t532) - m(5) * t417 - m(6) * t368 - m(4) * t304 + t525 + t573) * g(3) - t533 * t536 + t564 + (-t42 + t440) * t123 + (-t440 - t52) * t122 + (-t53 + t439) * t121 + (-t41 + t439) * t120 + (t195 * t9 + t206 * t346 - t37 * t42 - t38 * t41 - g(1) * (t240 + t269) - g(2) * (t239 + t266) + t545 * t48) * m(7) + t545 * t117 + (-Ifges(5,4) * t497 + t513 + t129 / 0.2e1 + t475) * t341 + (Ifges(5,4) * t498 + t514 + t476 - t130 / 0.2e1) * t370 + t565 * (pkin(9) + t222) - t183 * t297 - t352 * t406 / 0.2e1 + (-t184 * t190 - t185 * t191 - t256 * t297 + (t134 * t309 + t135 * t308) * pkin(2)) * m(4) - t185 * t454 + (t206 * t345 + t218 * t26 - g(1) * (t240 + t365) - g(2) * (t239 + t366) - t45 * t52 - t46 * t53 + t533 * t89) * m(6) + t534 * t161 + (-t194 * t198 + t221 * t31 + t222 * t30 - t533 * t95 + t534 * t96) * m(5) + (t401 + (-t331 / 0.2e1 + t336) * qJD(1)) * qJD(1) + Ifges(3,5) * t254 - t256 * (-mrSges(4,1) * t228 + mrSges(4,2) * t227) + Ifges(3,6) * t253 - t241 * mrSges(3,2) - t242 * mrSges(3,1) - qJD(2) * (Ifges(4,5) * t227 + Ifges(4,6) * t228) / 0.2e1 + t221 * t91 + t222 * t92 + t218 * t28 + Ifges(4,5) * t193 + t195 * t27 - t198 * t132 - t191 * t200 - t190 * t201 + Ifges(4,6) * t192 - t135 * mrSges(4,2) + t134 * mrSges(4,1) + (pkin(7) * t272 + t236 / 0.2e1) * t416 + t551 * t497 + t227 * t455 + t179 * t490 + t178 * t491 + t172 * t493 + t228 * (Ifges(4,1) * t227 + t465) / 0.2e1; -t370 * t161 - t227 * t200 - t228 * t201 + (-t117 + t536) * t341 + (-t169 * t581 - t552) * t315 + (t169 * t580 + t553) * t311 + t372 + t373 + (-g(1) * t314 + g(2) * t318) * (m(4) + t578) + (t169 * t346 + t2 * t311 - t315 * t4 - t341 * t48) * m(7) + (t169 * t345 + t311 * t6 + t315 * t7 - t341 * t89) * m(6) + (t341 * t95 - t370 * t96 + t156) * m(5) + (-t184 * t228 - t185 * t227 + t213) * m(4); t565 * pkin(9) + t564 + (-t161 + t457) * t95 + t567 * t361 + (-g(1) * t269 - g(2) * t266 + t265 * t9 - t37 * t44 - t38 * t43 + t544 * t48) * m(7) + t544 * t117 + (-pkin(4) * t26 - g(1) * t365 - g(2) * t366 - t45 * t56 - t46 * t57) * m(6) + t513 * t341 + t514 * t370 + (-t464 + t551) * t497 + (t456 + t518) * t96 + (-m(7) * (t418 + t532) - m(6) * t418 + t525) * g(3) + t265 * t27 - t43 * t120 - t57 * t121 - t56 * t122 - t44 * t123 + (t168 + t130) * t498 - pkin(4) * t28 + t129 * t496; t528 + (-t157 * t577 + t154 - t458 + t80) * t503 + (t216 * t559 + t217 * t558) * g(1) + (t559 * t214 + t558 * t215) * g(2) + (-t157 * t556 - t158 * t576) * t501 + t517 + (t358 + t360) * t481 + (-m(7) * t37 + t468 - t580) * t46 + (-m(7) * t38 - t469 + t581) * t45 + (-Ifges(6,2) * t158 - t155 + t550) * t504 - t48 * (mrSges(7,1) * t158 + mrSges(7,3) * t157) - t89 * (mrSges(6,1) * t158 - mrSges(6,2) * t157) + qJD(6) * t120 - t116 * t117 + qJ(6) * t32 - pkin(5) * t34 + (t348 * t481 - t116 * t48 - pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t38 - g(1) * (-pkin(5) * t216 + qJ(6) * t217) - g(2) * (-pkin(5) * t214 + qJ(6) * t215)) * m(7) + t38 * t470 + t37 * t471 + t83 * t502 + (Ifges(7,3) * t158 - t459) * t505; t158 * t117 - t169 * t120 + (-g(1) * t216 - g(2) * t214 - g(3) * t433 + t48 * t158 - t38 * t169 + t4) * m(7) + t34;];
tau  = t1;
