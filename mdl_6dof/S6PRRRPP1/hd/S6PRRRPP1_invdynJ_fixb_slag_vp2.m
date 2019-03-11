% Calculate vector of inverse dynamics joint torques for
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:52
% EndTime: 2019-03-08 22:43:38
% DurationCPUTime: 29.29s
% Computational Cost: add. (7933->703), mult. (18305->936), div. (0->0), fcn. (13607->14), ass. (0->328)
t294 = sin(qJ(2));
t289 = sin(pkin(6));
t393 = qJD(1) * t289;
t366 = t294 * t393;
t245 = qJD(2) * pkin(8) + t366;
t293 = sin(qJ(3));
t296 = cos(qJ(3));
t290 = cos(pkin(6));
t392 = qJD(1) * t290;
t176 = -t293 * t245 + t296 * t392;
t163 = -qJD(3) * pkin(3) - t176;
t292 = sin(qJ(4));
t295 = cos(qJ(4));
t386 = qJD(3) * t295;
t390 = qJD(2) * t293;
t237 = -t292 * t390 + t386;
t122 = -pkin(4) * t237 + qJD(5) + t163;
t423 = cos(pkin(11));
t288 = sin(pkin(11));
t177 = t296 * t245 + t293 * t392;
t164 = qJD(3) * pkin(9) + t177;
t252 = -pkin(3) * t296 - pkin(9) * t293 - pkin(2);
t297 = cos(qJ(2));
t365 = t297 * t393;
t179 = qJD(2) * t252 - t365;
t97 = t164 * t295 + t179 * t292;
t77 = qJ(5) * t237 + t97;
t428 = t288 * t77;
t388 = qJD(2) * t296;
t274 = qJD(4) - t388;
t238 = qJD(3) * t292 + t295 * t390;
t95 = -t164 * t292 + t295 * t179;
t76 = -qJ(5) * t238 + t95;
t59 = pkin(4) * t274 + t76;
t22 = t423 * t59 - t428;
t16 = -t274 * pkin(5) + qJD(6) - t22;
t136 = -t423 * t237 + t238 * t288;
t316 = t288 * t237 + t238 * t423;
t35 = pkin(5) * t136 - qJ(6) * t316 + t122;
t461 = -t136 / 0.2e1;
t520 = Ifges(6,4) - Ifges(7,5);
t448 = t274 / 0.2e1;
t457 = t316 / 0.2e1;
t519 = Ifges(7,4) + Ifges(6,5);
t521 = Ifges(6,1) + Ifges(7,1);
t562 = t448 * t519 + t457 * t521;
t565 = mrSges(6,2) * t122 + mrSges(7,2) * t16 - t22 * mrSges(6,3) - mrSges(7,3) * t35 + t520 * t461 + t562;
t449 = -t274 / 0.2e1;
t458 = -t316 / 0.2e1;
t460 = t136 / 0.2e1;
t564 = Ifges(6,4) * t460 + Ifges(7,5) * t461 + t449 * t519 + t458 * t521 - t565;
t563 = Ifges(6,4) * t461 + Ifges(7,5) * t460 + t562 + t565;
t334 = pkin(3) * t293 - pkin(9) * t296;
t241 = t334 * qJD(3);
t387 = qJD(3) * t293;
t374 = pkin(8) * t387;
t402 = t296 * t297;
t561 = t295 * t241 + t292 * t374 - (-t292 * t402 + t294 * t295) * t393;
t382 = qJD(4) * t295;
t407 = t292 * t294;
t560 = t292 * t241 + t252 * t382 - (t295 * t402 + t407) * t393;
t287 = qJ(4) + pkin(11);
t284 = sin(t287);
t285 = cos(t287);
t486 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t546 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t559 = t284 * t546 + t285 * t486;
t517 = Ifges(6,3) + Ifges(7,2);
t548 = Ifges(5,3) + t517;
t403 = t295 * t296;
t276 = pkin(8) * t403;
t321 = pkin(4) * t293 - qJ(5) * t403;
t381 = qJD(5) * t295;
t558 = -t293 * t381 + t321 * qJD(3) + (-t276 + (qJ(5) * t293 - t252) * t292) * qJD(4) + t561;
t405 = t293 * t295;
t557 = -(-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t405 - (-qJD(5) * t293 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t296) * t292 - t560;
t240 = t334 * qJD(2);
t116 = -t176 * t292 + t295 * t240;
t291 = -qJ(5) - pkin(9);
t347 = qJD(4) * t291;
t556 = -qJD(2) * t321 - qJD(5) * t292 + t295 * t347 - t116;
t117 = t295 * t176 + t292 * t240;
t362 = t292 * t388;
t554 = -qJ(5) * t362 - t292 * t347 + t117 - t381;
t69 = t423 * t77;
t23 = t288 * t59 + t69;
t18 = qJ(6) * t274 + t23;
t553 = Ifges(6,4) * t457 + Ifges(7,5) * t458 + Ifges(6,6) * t448 + Ifges(7,6) * t449 + (Ifges(6,2) + Ifges(7,3)) * t461 - t122 * mrSges(6,1) - t35 * mrSges(7,1) + t18 * mrSges(7,2) + t23 * mrSges(6,3);
t380 = qJD(2) * qJD(3);
t243 = qJDD(2) * t293 + t296 * t380;
t133 = qJD(4) * t237 + qJDD(3) * t292 + t243 * t295;
t134 = -qJD(4) * t238 + qJDD(3) * t295 - t243 * t292;
t60 = t133 * t288 - t134 * t423;
t467 = -t60 / 0.2e1;
t61 = t133 * t423 + t288 * t134;
t465 = t61 / 0.2e1;
t242 = qJDD(2) * t296 - t293 * t380;
t231 = qJDD(4) - t242;
t453 = t231 / 0.2e1;
t526 = -m(6) - m(7);
t542 = pkin(4) * t526;
t518 = Ifges(6,6) - Ifges(7,6);
t551 = -mrSges(5,1) + t542;
t549 = Ifges(6,2) * t461 - Ifges(7,3) * t460 + t520 * t457 + t553;
t391 = qJD(2) * t289;
t356 = qJD(1) * t391;
t261 = t294 * t356;
t379 = qJDD(1) * t289;
t201 = t297 * t379 - t261;
t187 = -qJDD(2) * pkin(2) - t201;
t105 = -pkin(3) * t242 - pkin(9) * t243 + t187;
t378 = qJDD(1) * t290;
t262 = t297 * t356;
t202 = t294 * t379 + t262;
t539 = qJDD(2) * pkin(8) + qJD(3) * t392 + t202;
t86 = -t245 * t387 + t293 * t378 + t296 * t539;
t83 = qJDD(3) * pkin(9) + t86;
t15 = -qJD(4) * t97 + t295 * t105 - t292 * t83;
t7 = pkin(4) * t231 - qJ(5) * t133 - qJD(5) * t238 + t15;
t384 = qJD(4) * t292;
t14 = t292 * t105 - t164 * t384 + t179 * t382 + t295 * t83;
t9 = qJ(5) * t134 + qJD(5) * t237 + t14;
t3 = -t288 * t9 + t423 * t7;
t2 = -t231 * pkin(5) + qJDD(6) - t3;
t385 = qJD(3) * t296;
t87 = -t245 * t385 - t293 * t539 + t296 * t378;
t84 = -qJDD(3) * pkin(3) - t87;
t36 = -pkin(4) * t134 + qJDD(5) + t84;
t466 = t60 / 0.2e1;
t5 = pkin(5) * t60 - qJ(6) * t61 - qJD(6) * t316 + t36;
t545 = mrSges(6,2) * t36 + mrSges(7,2) * t2 - mrSges(6,3) * t3 - mrSges(7,3) * t5 + Ifges(7,5) * t466 + 0.2e1 * t453 * t519 + 0.2e1 * t465 * t521 + (Ifges(6,4) + t520) * t467;
t541 = -m(4) + t526;
t512 = t288 * t557 + t423 * t558;
t511 = t288 * t558 - t423 * t557;
t74 = mrSges(7,1) * t136 - mrSges(7,3) * t316;
t75 = mrSges(6,1) * t136 + mrSges(6,2) * t316;
t540 = t74 + t75;
t506 = t288 * t554 + t423 * t556;
t503 = t288 * t556 - t423 * t554;
t499 = -t177 + (-t362 + t384) * pkin(4);
t331 = t292 * mrSges(5,1) + t295 * mrSges(5,2);
t538 = -t331 + mrSges(3,2) - mrSges(4,3);
t309 = t292 * t385 + t293 * t382;
t532 = -Ifges(6,2) * t460 + Ifges(7,3) * t461 - t449 * t518 - t458 * t520 + t553;
t529 = -m(5) * pkin(8) - t284 * t486 + t285 * t546 + t292 * t542 + t538;
t281 = pkin(4) * t295 + pkin(3);
t332 = -mrSges(5,1) * t295 + mrSges(5,2) * t292;
t314 = m(5) * pkin(3) - t332;
t333 = mrSges(4,1) * t296 - mrSges(4,2) * t293;
t373 = m(5) * pkin(9) + mrSges(5,3);
t490 = t293 * t373 + t296 * t314 + mrSges(3,1) + t333;
t522 = mrSges(6,3) + mrSges(7,2);
t528 = t490 + (t291 * t526 + t522) * t293 + (-t281 * t526 + t559) * t296;
t527 = -m(4) - m(5);
t463 = t133 / 0.2e1;
t462 = t134 / 0.2e1;
t525 = t242 / 0.2e1;
t524 = t243 / 0.2e1;
t515 = qJ(6) * t387 - qJD(6) * t296 + t511;
t514 = -pkin(5) * t387 - t512;
t20 = t60 * mrSges(7,1) - t61 * mrSges(7,3);
t21 = t60 * mrSges(6,1) + t61 * mrSges(6,2);
t513 = t20 + t21;
t37 = -mrSges(7,2) * t60 + mrSges(7,3) * t231;
t38 = -mrSges(6,2) * t231 - mrSges(6,3) * t60;
t510 = t37 + t38;
t39 = mrSges(6,1) * t231 - mrSges(6,3) * t61;
t40 = -t231 * mrSges(7,1) + t61 * mrSges(7,2);
t509 = t40 - t39;
t342 = t423 * t292;
t233 = t288 * t295 + t342;
t185 = t233 * t388;
t341 = t423 * t295;
t412 = t288 * t292;
t315 = t341 - t412;
t307 = t296 * t315;
t186 = qJD(2) * t307;
t212 = t233 * qJD(4);
t213 = t315 * qJD(4);
t507 = -qJD(6) * t233 + t499 + (t186 - t213) * qJ(6) + (-t185 + t212) * pkin(5);
t505 = pkin(5) * t390 - t506;
t504 = -qJ(6) * t390 + t503;
t65 = -mrSges(5,1) * t134 + mrSges(5,2) * t133;
t502 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t243 + t65;
t501 = (-t293 * t386 - t296 * t384) * pkin(8) + t560;
t190 = t292 * t252 + t276;
t500 = -qJD(4) * t190 + t561;
t113 = mrSges(6,1) * t274 - mrSges(6,3) * t316;
t114 = -mrSges(7,1) * t274 + mrSges(7,2) * t316;
t400 = t113 - t114;
t372 = mrSges(4,3) * t390;
t498 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t237 + mrSges(5,2) * t238 + t372;
t230 = Ifges(5,4) * t237;
t126 = t238 * Ifges(5,1) + t274 * Ifges(5,5) + t230;
t282 = Ifges(4,4) * t388;
t497 = Ifges(4,1) * t390 + Ifges(4,5) * qJD(3) + t295 * t126 + t282;
t29 = t423 * t76 - t428;
t496 = -t29 + qJD(6);
t383 = qJD(4) * t293;
t495 = t292 * t383 - t295 * t385;
t494 = -t293 * t87 + t296 * t86;
t493 = t14 * t295 - t15 * t292;
t491 = t238 * Ifges(5,5) + t237 * Ifges(5,6) - t136 * t518 + t274 * t548 + t316 * t519;
t489 = Ifges(5,5) * t133 + Ifges(5,6) * t134 + t231 * t548 - t518 * t60 + t519 * t61;
t487 = mrSges(4,2) - t373 - t522;
t485 = mrSges(4,1) + t314 + t559;
t483 = -m(5) * t163 - t498;
t435 = Ifges(4,4) * t293;
t327 = t296 * Ifges(4,2) + t435;
t479 = t16 * mrSges(7,1) + t23 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t327 / 0.2e1 - t18 * mrSges(7,3) - t22 * mrSges(6,1);
t4 = t288 * t7 + t423 * t9;
t1 = qJ(6) * t231 + qJD(6) * t274 + t4;
t478 = -t15 * mrSges(5,1) - t3 * mrSges(6,1) + t2 * mrSges(7,1) + t14 * mrSges(5,2) + t4 * mrSges(6,2) - t1 * mrSges(7,3);
t473 = mrSges(6,1) * t36 + mrSges(7,1) * t5 - mrSges(7,2) * t1 - mrSges(6,3) * t4 + 0.2e1 * Ifges(7,3) * t466 - t61 * Ifges(6,4) / 0.2e1 - t231 * Ifges(6,6) / 0.2e1 + Ifges(7,6) * t453 + (-t520 + Ifges(7,5)) * t465 + (-t467 + t466) * Ifges(6,2);
t298 = qJD(2) ^ 2;
t470 = Ifges(5,1) * t463 + Ifges(5,4) * t462 + Ifges(5,5) * t453;
t450 = t238 / 0.2e1;
t445 = pkin(4) * t238;
t444 = pkin(4) * t288;
t286 = t293 * pkin(8);
t437 = mrSges(5,3) * t237;
t436 = mrSges(5,3) * t238;
t434 = Ifges(4,4) * t296;
t433 = Ifges(5,4) * t292;
t432 = Ifges(5,4) * t295;
t429 = t238 * Ifges(5,4);
t427 = t293 * t84;
t424 = cos(pkin(10));
t422 = sin(pkin(10));
t411 = t289 * t294;
t410 = t289 * t297;
t408 = t292 * t293;
t406 = t292 * t296;
t112 = -mrSges(6,2) * t274 - mrSges(6,3) * t136;
t115 = -mrSges(7,2) * t136 + mrSges(7,3) * t274;
t401 = t112 + t115;
t235 = t295 * t252;
t144 = -qJ(5) * t405 + t235 + (-pkin(8) * t292 - pkin(4)) * t296;
t156 = -qJ(5) * t408 + t190;
t79 = t288 * t144 + t423 * t156;
t394 = pkin(2) * t410 + pkin(8) * t411;
t244 = pkin(4) * t408 + t286;
t389 = qJD(2) * t294;
t283 = pkin(8) * t385;
t371 = mrSges(4,3) * t388;
t370 = t293 * t410;
t369 = t284 * t410;
t178 = pkin(4) * t309 + t283;
t367 = t423 * pkin(4);
t364 = t289 * t389;
t363 = t297 * t391;
t125 = t237 * Ifges(5,2) + t274 * Ifges(5,6) + t429;
t357 = -t292 * t125 / 0.2e1;
t346 = t289 * t424;
t345 = t289 * t422;
t344 = t424 * t294;
t343 = t424 * t297;
t340 = t422 * t294;
t339 = t422 * t297;
t338 = t380 / 0.2e1;
t329 = Ifges(5,1) * t295 - t433;
t328 = Ifges(5,1) * t292 + t432;
t326 = -Ifges(5,2) * t292 + t432;
t325 = Ifges(5,2) * t295 + t433;
t324 = Ifges(4,5) * t296 - Ifges(4,6) * t293;
t323 = Ifges(5,5) * t295 - Ifges(5,6) * t292;
t322 = Ifges(5,5) * t292 + Ifges(5,6) * t295;
t216 = t290 * t293 + t296 * t411;
t154 = -t216 * t292 - t295 * t410;
t320 = -t216 * t295 + t292 * t410;
t215 = -t290 * t296 + t293 * t411;
t319 = t163 * t331;
t246 = -qJD(2) * pkin(2) - t365;
t318 = t246 * (mrSges(4,1) * t293 + mrSges(4,2) * t296);
t317 = t293 * (Ifges(4,1) * t296 - t435);
t78 = t144 * t423 - t288 * t156;
t209 = t290 * t344 + t339;
t148 = t209 * t293 + t296 * t346;
t211 = -t290 * t340 + t343;
t150 = t211 * t293 - t296 * t345;
t312 = -g(1) * t150 - g(2) * t148 - g(3) * t215;
t303 = Ifges(5,5) * t293 + t296 * t329;
t302 = Ifges(5,6) * t293 + t296 * t326;
t301 = Ifges(5,3) * t293 + t296 * t323;
t280 = -t367 - pkin(5);
t277 = qJ(6) + t444;
t256 = t291 * t295;
t255 = -qJD(3) * mrSges(4,2) + t371;
t239 = t333 * qJD(2);
t210 = t290 * t339 + t344;
t208 = -t290 * t343 + t340;
t203 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t242;
t200 = t210 * pkin(2);
t199 = t208 * pkin(2);
t198 = t315 * t293;
t197 = t233 * t293;
t189 = -pkin(8) * t406 + t235;
t175 = mrSges(5,1) * t274 - t436;
t174 = -mrSges(5,2) * t274 + t437;
t161 = -t256 * t423 + t291 * t412;
t160 = -t256 * t288 - t291 * t342;
t157 = -mrSges(4,1) * t242 + mrSges(4,2) * t243;
t153 = -qJD(3) * t215 + t296 * t363;
t152 = qJD(3) * t216 + t293 * t363;
t151 = t211 * t296 + t293 * t345;
t149 = t209 * t296 - t293 * t346;
t139 = t216 * t284 + t285 * t410;
t131 = -pkin(5) * t315 - qJ(6) * t233 - t281;
t119 = qJD(3) * t307 - t233 * t383;
t118 = t288 * t495 - t341 * t383 - t342 * t385;
t101 = pkin(5) * t197 - qJ(6) * t198 + t244;
t100 = -mrSges(5,2) * t231 + mrSges(5,3) * t134;
t99 = mrSges(5,1) * t231 - mrSges(5,3) * t133;
t90 = t151 * t284 - t210 * t285;
t88 = t149 * t284 - t208 * t285;
t81 = t288 * t154 - t320 * t423;
t80 = -t154 * t423 - t288 * t320;
t72 = t296 * pkin(5) - t78;
t71 = -qJ(6) * t296 + t79;
t68 = qJD(4) * t154 + t153 * t295 + t292 * t364;
t67 = qJD(4) * t320 - t153 * t292 + t295 * t364;
t43 = pkin(5) * t316 + qJ(6) * t136 + t445;
t41 = t133 * Ifges(5,4) + t134 * Ifges(5,2) + t231 * Ifges(5,6);
t30 = -pkin(5) * t118 - qJ(6) * t119 - qJD(6) * t198 + t178;
t28 = t288 * t76 + t69;
t27 = t288 * t67 + t423 * t68;
t25 = t288 * t68 - t423 * t67;
t6 = [m(2) * qJDD(1) - t320 * t100 + t153 * t255 + t154 * t99 + t68 * t174 + t67 * t175 + t216 * t203 + t510 * t81 + t509 * t80 + t401 * t27 - t400 * t25 + (t502 + t513) * t215 + (t498 + t540) * t152 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t298 - t157) * t297 + (-mrSges(3,1) * t298 - mrSges(3,2) * qJDD(2) - qJD(2) * t239) * t294) * t289 + (-m(2) - m(3) + t526 + t527) * g(3) + m(4) * (-t152 * t176 + t153 * t177 - t215 * t87 + t216 * t86 + (-t187 * t297 + t246 * t389) * t289) + m(3) * (qJDD(1) * t290 ^ 2 + (t201 * t297 + t202 * t294) * t289) + m(7) * (t1 * t81 + t152 * t35 + t16 * t25 + t18 * t27 + t2 * t80 + t215 * t5) + m(6) * (t122 * t152 + t215 * t36 - t22 * t25 + t23 * t27 - t3 * t80 + t4 * t81) + m(5) * (-t14 * t320 + t15 * t154 + t152 * t163 + t215 * t84 + t67 * t95 + t68 * t97); (-t176 * t385 + t494) * mrSges(4,3) + (-t202 + t262) * mrSges(3,2) + t473 * t197 + (-m(6) * t122 - m(7) * t35 + t483 - t540) * t293 * t365 + (-t518 * t197 + t293 * t323 - t296 * t548) * t453 + t549 * t118 + (qJD(3) * t301 + t118 * t518 - t322 * t383) * t448 + t545 * t198 + t237 * (qJD(3) * t302 - t325 * t383) / 0.2e1 - (t295 * t125 + t292 * t126) * t383 / 0.2e1 - pkin(2) * t157 + t101 * t20 + t78 * t39 + t79 * t38 + t71 * t37 + t72 * t40 + t30 * t74 + t434 * t524 + t327 * t525 + (t201 + t261) * mrSges(3,1) - t41 * t408 / 0.2e1 - t187 * t333 + qJD(3) ^ 2 * t324 / 0.2e1 + (-t14 * t408 - t15 * t405 - t309 * t97 + t495 * t95) * mrSges(5,3) + t163 * (mrSges(5,1) * t309 - mrSges(5,2) * t495) + t563 * t119 + t239 * t366 + (t527 * t394 - t522 * t370 + t526 * (-t291 * t370 + t394) - t546 * (-t285 * t411 + t296 * t369) + (t526 * (pkin(4) * t407 + t281 * t402) - t486 * (t284 * t294 + t285 * t402) - t490 * t297 + t538 * t294) * t289) * g(3) - t255 * t374 + Ifges(3,3) * qJDD(2) + t405 * t470 + (t491 / 0.2e1 - t177 * mrSges(4,3) + t95 * mrSges(5,1) - t97 * mrSges(5,2) + Ifges(6,6) * t461 + Ifges(7,6) * t460 + t457 * t519 - t479 + t517 * t448) * t387 + qJD(3) * t318 + (m(5) * t200 + t541 * (pkin(8) * t211 - t200) + t529 * t211 + t528 * t210) * g(1) + (m(5) * t199 + t541 * (pkin(8) * t209 - t199) + t529 * t209 + t528 * t208) * g(2) + t178 * t75 + t189 * t99 + t190 * t100 - t489 * t296 / 0.2e1 + (-pkin(2) * t187 + ((-t176 * t296 - t177 * t293) * qJD(3) + t494) * pkin(8) - (t246 * t294 + (-t176 * t293 + t177 * t296) * t297) * t393) * m(4) + t497 * t385 / 0.2e1 + t498 * t283 + t500 * t175 + t501 * t174 + (t14 * t190 + t15 * t189 + (t163 * t385 + t427) * pkin(8) + t501 * t97 + t500 * t95) * m(5) + t317 * t338 + t357 * t385 + t331 * t427 + (qJD(3) * t303 - t328 * t383) * t450 + t244 * t21 + t502 * t286 + t511 * t112 + t512 * t113 + (t122 * t178 + t22 * t512 + t23 * t511 + t244 * t36 + t3 * t78 + t4 * t79) * m(6) + t514 * t114 + t515 * t115 + (t1 * t71 + t101 * t5 + t16 * t514 + t18 * t515 + t2 * t72 + t30 * t35) * m(7) + (Ifges(4,4) * t524 + Ifges(4,2) * t525 - t255 * t365 + pkin(8) * t203 - Ifges(7,6) * t466 - Ifges(6,6) * t467 - Ifges(5,6) * t462 - Ifges(5,5) * t463 + (-Ifges(4,2) * t293 + t434) * t338 - t519 * t465 + t478) * t296 + (Ifges(4,1) * t243 + Ifges(4,4) * t525 + t326 * t462 + t329 * t463) * t293 + qJDD(3) * (Ifges(4,5) * t293 + Ifges(4,6) * t296); t532 * t185 + (t371 - t255) * t176 + (-t518 * t448 - t549) * t212 + t545 * t233 + (t237 * t326 + t238 * t329 + t274 * t323) * qJD(4) / 0.2e1 - (t237 * t302 + t238 * t303 + t274 * t301) * qJD(2) / 0.2e1 + t125 * t362 / 0.2e1 + t131 * t20 - t86 * mrSges(4,2) + t87 * mrSges(4,1) - pkin(3) * t65 - t324 * t380 / 0.2e1 + t126 * t382 / 0.2e1 - (-t518 * t453 + t473) * t315 + (-t318 - t95 * (mrSges(5,1) * t293 - mrSges(5,3) * t403) - t97 * (-mrSges(5,2) * t293 - mrSges(5,3) * t406)) * qJD(2) + (-pkin(3) * t84 - t116 * t95 - t117 * t97) * m(5) + (t357 + t319) * qJD(4) + t84 * t332 - t319 * t388 - t298 * t317 / 0.2e1 + t563 * t213 + t564 * t186 + Ifges(4,3) * qJDD(3) + t292 * t470 + t325 * t462 + t328 * t463 + t322 * t453 - t117 * t174 - t116 * t175 + (t372 + t483) * t177 - t491 * t390 / 0.2e1 + (-t175 * t382 - t174 * t384 + m(5) * ((-t292 * t97 - t295 * t95) * qJD(4) + t493) + t295 * t100 - t292 * t99) * pkin(9) + (-t382 * t95 - t384 * t97 + t493) * mrSges(5,3) - (-Ifges(4,2) * t390 + t282 + t497) * t388 / 0.2e1 + t499 * t75 + Ifges(4,6) * t242 + Ifges(4,5) * t243 + t503 * t112 + t504 * t115 + t505 * t114 + t506 * t113 + (t122 * t499 - t160 * t3 + t161 * t4 + t22 * t506 + t23 * t503 - t281 * t36) * m(6) - t281 * t21 + (t1 * t161 + t131 * t5 + t16 * t505 + t160 * t2 + t18 * t504 + t35 * t507) * m(7) + t507 * t74 + t509 * t160 + t510 * t161 + (Ifges(6,6) * t460 + Ifges(7,6) * t461 + t449 * t517 + t458 * t519 + t479) * t390 + t295 * t41 / 0.2e1 + (t526 * (-t215 * t281 - t216 * t291) + t487 * t216 + t485 * t215) * g(3) + (t526 * (-t150 * t281 - t151 * t291) + t487 * t151 + t485 * t150) * g(1) + (t526 * (-t148 * t281 - t149 * t291) + t487 * t149 + t485 * t148) * g(2); t532 * t316 + t489 - (-Ifges(5,2) * t238 + t126 + t230) * t237 / 0.2e1 - t29 * t112 - t43 * t74 - t238 * (Ifges(5,1) * t237 - t429) / 0.2e1 + (t175 + t436) * t97 - t564 * t136 + (-(-t149 * t295 - t208 * t292) * mrSges(5,2) - t546 * (t149 * t285 + t208 * t284) + t486 * t88 + t551 * (-t149 * t292 + t208 * t295)) * g(2) + (-(-t151 * t295 - t210 * t292) * mrSges(5,2) - t546 * (t151 * t285 + t210 * t284) + t486 * t90 + t551 * (-t151 * t292 + t210 * t295)) * g(1) + (-t320 * mrSges(5,2) - t546 * (t216 * t285 - t369) + t486 * t139 + t551 * t154) * g(3) + (-t174 + t437) * t95 + ((t288 * t4 + t3 * t423) * pkin(4) - t122 * t445 + t22 * t28 - t23 * t29) * m(6) - t478 + (t1 * t277 - t16 * t28 + t18 * t496 + t2 * t280 - t35 * t43) * m(7) + t496 * t115 + t400 * t28 + t39 * t367 + t38 * t444 + (Ifges(5,5) * t237 - Ifges(5,6) * t238) * t449 + t125 * t450 - t163 * (mrSges(5,1) * t238 + mrSges(5,2) * t237) + t277 * t37 + t280 * t40 - t75 * t445; t400 * t316 + t401 * t136 + (t136 * t18 - t16 * t316 + t312 + t5) * m(7) + (t136 * t23 + t22 * t316 + t312 + t36) * m(6) + t513; -t274 * t115 + t316 * t74 + (-g(1) * t90 - g(2) * t88 - g(3) * t139 - t18 * t274 + t316 * t35 + t2) * m(7) + t40;];
tau  = t6;
