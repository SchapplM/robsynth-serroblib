% Calculate vector of inverse dynamics joint torques for
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:52
% EndTime: 2019-03-08 22:49:29
% DurationCPUTime: 24.79s
% Computational Cost: add. (5160->720), mult. (11933->910), div. (0->0), fcn. (8583->10), ass. (0->329)
t533 = Ifges(7,4) + Ifges(6,5);
t534 = Ifges(6,4) + Ifges(5,5);
t532 = Ifges(7,2) + Ifges(6,3);
t531 = Ifges(5,6) - Ifges(6,6);
t540 = Ifges(7,5) - t534;
t526 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t245 = sin(qJ(3));
t434 = Ifges(4,4) * t245;
t248 = cos(qJ(3));
t448 = t248 / 0.2e1;
t536 = Ifges(4,2) * t448 + t434 / 0.2e1;
t530 = Ifges(7,6) - Ifges(6,6);
t247 = cos(qJ(4));
t538 = t533 * t247;
t244 = sin(qJ(4));
t537 = t533 * t244;
t376 = qJD(2) * qJD(3);
t198 = qJDD(2) * t245 + t248 * t376;
t377 = t247 * qJD(3);
t391 = qJD(2) * t245;
t192 = t244 * t391 - t377;
t386 = qJD(4) * t192;
t89 = qJDD(3) * t244 + t198 * t247 - t386;
t464 = t89 / 0.2e1;
t193 = qJD(3) * t244 + t247 * t391;
t90 = qJD(4) * t193 - t247 * qJDD(3) + t198 * t244;
t463 = -t90 / 0.2e1;
t462 = t90 / 0.2e1;
t197 = t248 * qJDD(2) - t245 * t376;
t190 = qJDD(4) - t197;
t458 = t190 / 0.2e1;
t456 = t192 / 0.2e1;
t455 = -t193 / 0.2e1;
t454 = t193 / 0.2e1;
t389 = qJD(2) * t248;
t230 = -qJD(4) + t389;
t453 = -t230 / 0.2e1;
t452 = t230 / 0.2e1;
t535 = qJD(3) / 0.2e1;
t529 = -Ifges(5,3) - Ifges(6,2);
t528 = t533 * t193;
t490 = t244 * t532 + t538;
t488 = -t531 * t244 + t247 * t534;
t417 = qJ(5) * t244;
t285 = pkin(4) * t247 + t417;
t527 = t533 * t192;
t431 = Ifges(5,4) * t244;
t475 = t247 * t526 - t431 + t537;
t524 = (-Ifges(5,4) + t533) * t462 + t526 * t464 - t540 * t458;
t510 = t192 * t532 + t230 * t530 + t528;
t246 = sin(qJ(2));
t242 = sin(pkin(6));
t396 = qJD(1) * t242;
t360 = t246 * t396;
t199 = qJD(2) * pkin(8) + t360;
t243 = cos(pkin(6));
t395 = qJD(1) * t243;
t358 = t245 * t395;
t137 = t199 * t248 + t358;
t120 = qJD(3) * pkin(9) + t137;
t447 = pkin(3) * t248;
t205 = -pkin(9) * t245 - pkin(2) - t447;
t249 = cos(qJ(2));
t359 = t249 * t396;
t139 = qJD(2) * t205 - t359;
t39 = -t244 * t120 + t247 * t139;
t25 = qJ(6) * t193 + t39;
t522 = qJD(5) - t25;
t521 = qJD(5) - t39;
t392 = qJD(2) * t242;
t349 = qJD(1) * t392;
t216 = t249 * t349;
t375 = qJDD(1) * t242;
t161 = t246 * t375 + t216;
t520 = qJDD(2) * pkin(8) + qJD(3) * t395 + t161;
t308 = mrSges(6,1) * t247 + mrSges(6,3) * t244;
t310 = mrSges(5,1) * t247 - mrSges(5,2) * t244;
t366 = m(7) * pkin(5) + mrSges(7,1);
t423 = t244 * mrSges(7,2);
t519 = t247 * t366 + mrSges(4,1) + t308 + t310 + t423;
t518 = -qJD(5) * t244 - t358;
t189 = Ifges(5,4) * t192;
t480 = t526 * t193 + t230 * t540 - t189 + t527;
t236 = Ifges(4,4) * t389;
t517 = Ifges(4,1) * t391 + Ifges(4,5) * qJD(3) + t244 * t510 + t247 * t480 + t236;
t516 = Ifges(7,5) * t454 + Ifges(4,6) * t535 + Ifges(7,3) * t452 + qJD(2) * t536 + t453 * t529 + t455 * t534 + (Ifges(7,6) + t531) * t456;
t374 = qJDD(1) * t243;
t388 = qJD(3) * t245;
t35 = -t199 * t388 + t245 * t374 + t248 * t520;
t33 = qJDD(3) * pkin(9) + t35;
t382 = qJD(4) * t247;
t384 = qJD(4) * t244;
t215 = t246 * t349;
t160 = t249 * t375 - t215;
t143 = -qJDD(2) * pkin(2) - t160;
t59 = -pkin(3) * t197 - pkin(9) * t198 + t143;
t6 = -t120 * t384 + t139 * t382 + t244 * t59 + t247 * t33;
t3 = t190 * qJ(5) - t230 * qJD(5) + t6;
t2 = qJ(6) * t90 + qJD(6) * t192 + t3;
t515 = Ifges(5,4) * t464 - t533 * t89 / 0.2e1 - t2 * mrSges(7,3) + (Ifges(5,2) + t532) * t463 + (Ifges(5,6) + t530) * t458;
t222 = t230 * qJ(5);
t40 = t247 * t120 + t244 * t139;
t26 = qJ(6) * t192 + t40;
t18 = -t222 + t26;
t28 = -t222 + t40;
t514 = -t28 * mrSges(6,2) - t40 * mrSges(5,3) + t18 * mrSges(7,3);
t461 = -m(6) - m(7);
t512 = mrSges(3,2) - mrSges(4,3);
t49 = mrSges(5,1) * t190 - mrSges(5,3) * t89;
t50 = -t190 * mrSges(6,1) + t89 * mrSges(6,2);
t511 = -t49 + t50;
t419 = cos(pkin(10));
t328 = t419 * t246;
t241 = sin(pkin(10));
t412 = t241 * t249;
t169 = t243 * t328 + t412;
t329 = t242 * t419;
t104 = -t169 * t245 - t248 * t329;
t508 = t285 * t104;
t327 = t419 * t249;
t413 = t241 * t246;
t171 = -t243 * t413 + t327;
t410 = t242 * t248;
t106 = -t171 * t245 + t241 * t410;
t507 = t285 * t106;
t21 = mrSges(5,1) * t90 + mrSges(5,2) * t89;
t506 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t198 + t21;
t311 = mrSges(4,1) * t248 - mrSges(4,2) * t245;
t504 = t311 + mrSges(3,1);
t438 = pkin(9) - qJ(6);
t503 = -m(7) * t438 + mrSges(4,2) + mrSges(7,3);
t211 = t438 * t247;
t405 = t247 * t248;
t460 = pkin(4) + pkin(5);
t136 = -t245 * t199 + t248 * t395;
t314 = pkin(3) * t245 - pkin(9) * t248;
t195 = t314 * qJD(2);
t70 = -t244 * t136 + t195 * t247;
t502 = -(-qJ(6) * t405 - t245 * t460) * qJD(2) + t70 + qJD(4) * t211 - qJD(6) * t244;
t355 = t244 * t389;
t378 = qJD(6) * t247;
t71 = t247 * t136 + t244 * t195;
t56 = qJ(5) * t391 + t71;
t501 = -qJ(6) * t355 - t384 * t438 - t378 - t56;
t416 = qJ(5) * t247;
t280 = -t244 * t460 + t416;
t500 = -(qJD(2) * t280 - t199) * t248 + qJD(4) * t280 - t518;
t284 = pkin(4) * t244 - t416;
t499 = -(qJD(2) * t284 + t199) * t248 + qJD(4) * t284 + t518;
t101 = mrSges(6,1) * t192 - mrSges(6,3) * t193;
t102 = -mrSges(7,1) * t192 + mrSges(7,2) * t193;
t403 = t101 - t102;
t130 = -mrSges(7,2) * t230 + mrSges(7,3) * t192;
t436 = mrSges(5,3) * t192;
t131 = mrSges(5,2) * t230 - t436;
t135 = -mrSges(6,2) * t192 - mrSges(6,3) * t230;
t401 = t131 + t135;
t498 = t130 + t401;
t132 = mrSges(7,1) * t230 - mrSges(7,3) * t193;
t435 = mrSges(5,3) * t193;
t133 = -mrSges(5,1) * t230 - t435;
t134 = mrSges(6,1) * t230 + mrSges(6,2) * t193;
t400 = t134 - t133;
t497 = t132 + t400;
t411 = t242 * t246;
t173 = -t243 * t248 + t245 * t411;
t496 = t285 * t173;
t365 = mrSges(4,3) * t391;
t495 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t192 + mrSges(5,2) * t193 + t365;
t494 = qJ(5) * t388 - qJD(5) * t248;
t493 = -t245 * t530 + t248 * t490;
t492 = -t245 * t529 + t248 * t488;
t491 = -t247 * t532 + t537;
t489 = t244 * t534 + t247 * t531;
t487 = -t190 * t529 - t531 * t90 + t534 * t89;
t312 = mrSges(5,1) + mrSges(6,1) + t366;
t486 = -pkin(4) * t461 + t312;
t387 = qJD(3) * t248;
t36 = -t199 * t387 - t245 * t520 + t248 * t374;
t485 = -t245 * t36 + t248 * t35;
t7 = -t120 * t382 - t139 * t384 - t244 * t33 + t247 * t59;
t484 = -t244 * t7 + t247 * t6;
t273 = qJDD(5) - t7;
t4 = -pkin(4) * t190 + t273;
t483 = t244 * t4 + t247 * t3;
t482 = -mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t478 = -t245 * t540 + t248 * t475;
t315 = qJD(3) * pkin(3) + t136;
t268 = qJ(5) * t193 + t315;
t24 = -t192 * t460 + qJD(6) + t268;
t306 = -mrSges(7,1) * t244 + mrSges(7,2) * t247;
t307 = mrSges(6,1) * t244 - mrSges(6,3) * t247;
t309 = mrSges(5,1) * t244 + mrSges(5,2) * t247;
t41 = pkin(4) * t192 - t268;
t477 = -t24 * t306 - t41 * t307 + t309 * t315;
t430 = Ifges(5,4) * t247;
t476 = t244 * t526 + t430 - t538;
t473 = -t247 * t460 - t417;
t470 = m(5) * t315 - t495;
t468 = m(7) * qJ(6) - mrSges(6,2) - mrSges(5,3) + mrSges(7,3);
t105 = t169 * t248 - t245 * t329;
t107 = t241 * t242 * t245 + t171 * t248;
t174 = t243 * t245 + t246 * t410;
t467 = -g(1) * t107 - g(2) * t105 - g(3) * t174;
t1 = -qJ(6) * t89 - qJD(6) * t193 - t190 * t460 + t273;
t466 = -t7 * mrSges(5,1) + t4 * mrSges(6,1) + t1 * mrSges(7,1) + t6 * mrSges(5,2) - t2 * mrSges(7,2) - t3 * mrSges(6,3);
t251 = qJD(2) ^ 2;
t432 = Ifges(5,4) * t193;
t78 = -t192 * Ifges(5,2) - t230 * Ifges(5,6) + t432;
t465 = -t78 / 0.2e1;
t459 = -t190 / 0.2e1;
t457 = -t192 / 0.2e1;
t51 = mrSges(7,2) * t190 + mrSges(7,3) * t90;
t53 = -mrSges(6,2) * t90 + mrSges(6,3) * t190;
t437 = t51 + t53;
t433 = Ifges(4,4) * t248;
t271 = qJDD(3) * pkin(3) + t36;
t422 = t245 * t271;
t418 = qJ(5) * t192;
t168 = -t243 * t327 + t413;
t415 = t168 * t245;
t170 = t243 * t412 + t328;
t414 = t170 * t245;
t409 = t242 * t249;
t408 = t244 * t245;
t407 = t244 * t248;
t406 = t245 * t247;
t404 = t248 * t249;
t402 = t130 + t135;
t196 = t314 * qJD(3);
t399 = t244 * t196 + t205 * t382;
t232 = pkin(8) * t405;
t398 = qJD(4) * t232 + t205 * t384;
t397 = pkin(2) * t409 + pkin(8) * t411;
t146 = t244 * t205 + t232;
t390 = qJD(2) * t246;
t383 = qJD(4) * t245;
t380 = qJD(5) * t247;
t370 = pkin(8) * t388;
t368 = pkin(9) * t384;
t367 = pkin(9) * t382;
t364 = mrSges(4,3) * t389;
t363 = t245 * t409;
t362 = t244 * t409;
t361 = -pkin(8) * t244 - pkin(4);
t357 = t242 * t390;
t356 = t249 * t392;
t354 = t244 * t387;
t20 = -t90 * mrSges(7,1) + t89 * mrSges(7,2);
t96 = t104 * pkin(3);
t351 = pkin(9) * t105 + t96;
t97 = t106 * pkin(3);
t350 = pkin(9) * t107 + t97;
t333 = t382 / 0.2e1;
t48 = -t190 * mrSges(7,1) - t89 * mrSges(7,3);
t332 = -t168 * pkin(2) + pkin(8) * t169;
t331 = -t170 * pkin(2) + pkin(8) * t171;
t164 = t173 * pkin(3);
t330 = pkin(9) * t174 - t164;
t326 = t376 / 0.2e1;
t231 = pkin(8) * t407;
t145 = t205 * t247 - t231;
t325 = t242 * pkin(3) * t404 + pkin(9) * t363 + t397;
t323 = t248 * t359;
t124 = -qJ(5) * t248 + t146;
t313 = -t196 * t247 + t398;
t298 = -Ifges(5,2) * t244 + t430;
t297 = Ifges(5,2) * t247 + t431;
t292 = Ifges(4,5) * t248 - Ifges(4,6) * t245;
t287 = Ifges(7,5) * t247 + Ifges(7,6) * t244;
t286 = Ifges(7,5) * t244 - Ifges(7,6) * t247;
t283 = -pkin(9) * t415 - t168 * t447 + t332;
t282 = -pkin(9) * t414 - t170 * t447 + t331;
t281 = pkin(8) + t284;
t279 = qJ(5) * t461 + t482;
t112 = t174 * t244 + t247 * t409;
t200 = -qJD(2) * pkin(2) - t359;
t275 = t200 * (mrSges(4,1) * t245 + mrSges(4,2) * t248);
t274 = t245 * (Ifges(4,1) * t248 - t434);
t44 = t105 * t244 - t168 * t247;
t46 = t107 * t244 - t170 * t247;
t272 = -g(1) * t46 - g(2) * t44 - g(3) * t112;
t270 = -pkin(8) + t280;
t269 = Ifges(7,5) * t89 + Ifges(7,6) * t90 - Ifges(7,3) * t190;
t142 = (t244 * t246 + t247 * t404) * t242;
t257 = Ifges(5,6) * t245 + t248 * t298;
t253 = -Ifges(7,3) * t245 + t248 * t287;
t61 = (-t245 * t377 - t248 * t384) * pkin(8) + t399;
t252 = qJ(5) * t89 + qJD(5) * t193 + t271;
t240 = t248 * pkin(4);
t210 = t438 * t244;
t209 = -qJD(3) * mrSges(4,2) + t364;
t201 = -pkin(3) - t285;
t194 = t311 * qJD(2);
t186 = pkin(3) - t473;
t162 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t197;
t157 = t281 * t245;
t141 = -t247 * t411 + t248 * t362;
t125 = -t145 + t240;
t123 = t270 * t245;
t114 = -mrSges(4,1) * t197 + mrSges(4,2) * t198;
t113 = t174 * t247 - t362;
t111 = -qJD(3) * t173 + t248 * t356;
t110 = qJD(3) * t174 + t245 * t356;
t100 = pkin(4) * t193 + t418;
t95 = qJ(6) * t408 + t124;
t88 = pkin(5) * t248 + t231 + t240 + (-qJ(6) * t245 - t205) * t247;
t69 = -t170 * t405 + t171 * t244;
t68 = -t170 * t407 - t171 * t247;
t67 = -t168 * t405 + t169 * t244;
t66 = -t168 * t407 - t169 * t247;
t63 = -t193 * t460 - t418;
t62 = t244 * t370 - t313;
t60 = (qJD(4) * t285 - t380) * t245 + t281 * t387;
t58 = -pkin(4) * t391 - t70;
t54 = t361 * t388 + t313;
t52 = -mrSges(5,2) * t190 - mrSges(5,3) * t90;
t37 = t61 + t494;
t31 = (qJD(4) * t473 + t380) * t245 + t270 * t387;
t27 = pkin(4) * t230 + t521;
t23 = -qJD(4) * t112 + t111 * t247 + t244 * t357;
t22 = -qJD(4) * t362 + t111 * t244 + t174 * t382 - t247 * t357;
t19 = mrSges(6,1) * t90 - mrSges(6,3) * t89;
t17 = (-pkin(8) * qJD(3) + qJ(6) * qJD(4)) * t406 + (qJD(6) * t245 + (-pkin(8) * qJD(4) + qJ(6) * qJD(3)) * t248) * t244 + t399 + t494;
t16 = (-qJ(6) * t387 - t196) * t247 + (qJ(6) * t384 - t378 + (-pkin(5) + t361) * qJD(3)) * t245 + t398;
t15 = t230 * t460 + t522;
t8 = pkin(4) * t90 - t252;
t5 = -t460 * t90 + qJDD(6) + t252;
t9 = [m(2) * qJDD(1) + t111 * t209 + t174 * t162 + t498 * t23 + t497 * t22 + (t52 + t437) * t113 + (t48 + t511) * t112 + (t19 - t20 + t506) * t173 + (t403 + t495) * t110 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t251 - t114) * t249 + (-mrSges(3,1) * t251 - mrSges(3,2) * qJDD(2) - qJD(2) * t194) * t246) * t242 + (-m(2) - m(3) - m(4) - m(5) + t461) * g(3) + m(6) * (t110 * t41 + t112 * t4 + t113 * t3 + t173 * t8 + t22 * t27 + t23 * t28) + m(5) * (-t110 * t315 - t112 * t7 + t113 * t6 - t173 * t271 - t22 * t39 + t23 * t40) + m(7) * (t1 * t112 - t110 * t24 + t113 * t2 + t15 * t22 - t173 * t5 + t18 * t23) + m(4) * (-t110 * t136 + t111 * t137 - t173 * t36 + t174 * t35 + (-t143 * t249 + t200 * t390) * t242) + m(3) * (qJDD(1) * t243 ^ 2 + (t160 * t249 + t161 * t246) * t242); (-t3 * mrSges(6,2) - t6 * mrSges(5,3) - t515) * t408 + (Ifges(4,5) * t245 + 0.2e1 * Ifges(4,6) * t448) * qJDD(3) + (-(t200 * t246 + (-t136 * t245 + t137 * t248) * t249) * t396 - pkin(2) * t143) * m(4) - (t244 * t480 + t247 * t78) * t383 / 0.2e1 + (t162 * t248 + t506 * t245 + m(5) * (-t315 * t387 - t422) + ((-t136 * t248 - t137 * t245) * qJD(3) + t485) * m(4) + t495 * t387) * pkin(8) + (-mrSges(5,2) * t315 + mrSges(6,2) * t27 + mrSges(7,2) * t24 - mrSges(5,3) * t39 - mrSges(6,3) * t41 - mrSges(7,3) * t15) * (-t244 * t383 + t248 * t377) + (-t136 * t387 + t485) * mrSges(4,3) + (mrSges(6,2) * t4 - mrSges(5,3) * t7 - mrSges(7,3) * t1 + t524) * t406 + (-m(4) * t331 - m(5) * t282 + t461 * (t69 * pkin(4) + qJ(5) * t68 + t282) + t512 * t171 + t504 * t170 - t312 * t69 + t482 * t68 - t468 * t414) * g(1) + (-m(4) * t332 - m(5) * t283 + t461 * (t67 * pkin(4) + qJ(5) * t66 + t283) + t512 * t169 + t504 * t168 - t312 * t67 + t482 * t66 - t468 * t415) * g(2) + (-m(4) * t397 - m(5) * t325 + t461 * (t142 * pkin(4) + qJ(5) * t141 + t325) + (t246 * t512 - t249 * t504) * t242 + t468 * t363 - t312 * t142 + t482 * t141) * g(3) + (-m(5) * t40 - m(6) * t28 - m(7) * t18 - t498) * qJD(1) * t142 + (-m(6) * t41 + m(7) * t24 - t403 + t470) * t245 * t359 + (t245 * Ifges(4,1) + t433 / 0.2e1 + Ifges(4,4) * t448) * t198 + (t39 * mrSges(5,1) - t27 * mrSges(6,1) - t15 * mrSges(7,1) - t40 * mrSges(5,2) + t18 * mrSges(7,2) - t137 * mrSges(4,3) + t28 * mrSges(6,3) - t516) * t388 + ((-Ifges(4,2) * t245 + t433) * t326 + Ifges(7,3) * t459 - Ifges(5,6) * t463 + t530 * t462 + t529 * t458 + t466) * t248 + (t253 * t452 + t257 * t457 + t292 * t535 + t492 * t453 + t493 * t456 + t275) * qJD(3) + (-t370 - t323) * t209 - t487 * t248 / 0.2e1 + (t287 * t459 + t298 * t463 + t306 * t5 + t307 * t8 + t458 * t488 + t462 * t490) * t245 + (-t286 * t452 - t297 * t457 - t453 * t489 - t456 * t491) * t383 + (-mrSges(5,1) * t315 + mrSges(6,1) * t41 - mrSges(7,1) * t24 + t514) * (t245 * t382 + t354) + (t216 - t161) * mrSges(3,2) + m(5) * (t145 * t7 + t146 * t6 + t39 * t62 + t40 * t61) + m(7) * (t1 * t88 + t123 * t5 + t15 * t16 + t17 * t18 + t2 * t95 + t24 * t31) + m(6) * (t124 * t3 + t125 * t4 + t157 * t8 + t27 * t54 + t28 * t37 + t41 * t60) + t194 * t360 + (t215 + t160) * mrSges(3,1) + Ifges(3,3) * qJDD(2) + (m(5) * t39 - m(6) * t27 - m(7) * t15 - t497) * (t244 * t323 - t247 * t360) + (qJD(3) * t478 - t383 * t476) * t454 + t517 * t387 / 0.2e1 - t143 * t311 + t510 * t245 * t333 + t274 * t326 - t309 * t422 + 0.2e1 * t536 * t197 + t88 * t48 + t95 * t51 + t60 * t101 + t31 * t102 + (t475 * t245 + t248 * t540) * t464 - pkin(2) * t114 + t123 * t20 + t124 * t53 + t125 * t50 + t17 * t130 + t61 * t131 + t16 * t132 + t62 * t133 + t54 * t134 + t37 * t135 + t145 * t49 + t146 * t52 + t157 * t19 + t269 * t448 + t354 * t465; (t5 * mrSges(7,1) + t515) * t247 + (-t382 * t39 + t467 + t484) * mrSges(5,3) + t499 * t101 + t500 * t102 + t501 * t130 + (t1 * t210 + t15 * t502 + t18 * t501 + t186 * t5 + t2 * t211 + t24 * t500) * m(7) + t502 * t132 + (pkin(3) * t271 - t39 * t70 - t40 * t71) * m(5) + t271 * t310 + (t201 * t8 - t27 * t58 - t28 * t56 + t499 * t41) * m(6) + (((-t244 * t28 + t247 * t27) * qJD(4) + t483) * m(6) + ((-t244 * t40 - t247 * t39) * qJD(4) + t484) * m(5) + (t53 + t52) * t247 + t511 * t244) * pkin(9) + (t287 / 0.2e1 - t488 / 0.2e1) * qJD(4) * t230 + (-t298 / 0.2e1 + t490 / 0.2e1) * t386 + (t27 * t382 + t467 + t483) * mrSges(6,2) + t516 * t391 + (-t1 * t244 - t15 * t382) * mrSges(7,3) + (t367 - t58) * t134 - t251 * t274 / 0.2e1 + t489 * t458 + t491 * t462 + t5 * t423 + ((t492 / 0.2e1 - t253 / 0.2e1) * t230 + t478 * t455 - t40 * (-mrSges(5,2) * t245 - mrSges(5,3) * t407) - t28 * (-mrSges(6,2) * t407 + mrSges(6,3) * t245) - t18 * (mrSges(7,2) * t245 + mrSges(7,3) * t407) - t15 * (-mrSges(7,1) * t245 - mrSges(7,3) * t405) - t39 * (mrSges(5,1) * t245 - mrSges(5,3) * t405) - t27 * (-mrSges(6,1) * t245 + mrSges(6,2) * t405) - t275 + (t257 / 0.2e1 - t493 / 0.2e1) * t192) * qJD(2) + (t510 / 0.2e1 + t465 + t514) * t384 + (-m(5) * t350 - m(7) * (t97 + t507) - m(6) * (t350 + t507) + t503 * t107 - t519 * t106) * g(1) + (-m(5) * t351 - m(7) * (t96 + t508) - m(6) * (t351 + t508) + t503 * t105 - t519 * t104) * g(2) + (-m(5) * t330 - m(7) * (-t164 - t496) - m(6) * (t330 - t496) + t503 * t174 + t519 * t173) * g(3) + (-t368 - t71) * t131 + (t364 - t209) * t136 + (t454 * t475 - t477) * qJD(4) + (-t368 - t56) * t135 + t78 * t355 / 0.2e1 + (t365 + t470) * t137 + (-t367 - t70) * t133 + Ifges(4,3) * qJDD(3) + t476 * t464 + t477 * t389 - (-Ifges(4,2) * t391 + t236 + t517) * t389 / 0.2e1 - t8 * t308 - pkin(3) * t21 - t35 * mrSges(4,2) + t36 * mrSges(4,1) + t480 * t333 - t292 * t376 / 0.2e1 + t186 * t20 + Ifges(4,6) * t197 + Ifges(4,5) * t198 + t201 * t19 + t210 * t48 + t211 * t51 + t286 * t459 + t297 * t463 + t244 * t524; -t466 + (t193 * t532 - t527) * t457 + t315 * (mrSges(5,1) * t193 - mrSges(5,2) * t192) - t460 * t48 + (-t192 * t526 - t432 + t510 + t528) * t455 + t487 + (-t192 * t534 - t193 * t531) * t452 + t402 * qJD(5) + (-t401 - t436) * t39 + (t279 * (t105 * t247 + t168 * t244) + t486 * t44) * g(2) + (t279 * (t107 * t247 + t170 * t244) + t486 * t46) * g(1) + (t112 * t486 + t113 * t279) * g(3) + (-t400 + t435) * t40 + (-pkin(4) * t4 + qJ(5) * t3 - t100 * t41 - t27 * t40 + t28 * t521) * m(6) + (qJ(5) * t2 - t1 * t460 - t15 * t26 + t18 * t522 - t24 * t63) * m(7) + t437 * qJ(5) - t269 + (-Ifges(5,2) * t193 - t189 + t480) * t456 - pkin(4) * t50 + (-t15 * t192 - t18 * t193) * mrSges(7,3) + (t192 * t27 + t193 * t28) * mrSges(6,2) - t100 * t101 - t63 * t102 - t25 * t130 - t26 * t132 - t24 * (-mrSges(7,1) * t193 - mrSges(7,2) * t192) - t41 * (mrSges(6,1) * t193 + mrSges(6,3) * t192) + (-Ifges(7,5) * t192 + Ifges(7,6) * t193) * t453 + t78 * t454; t403 * t193 + t402 * t230 + t48 + t50 + (t18 * t230 - t193 * t24 + t1 + t272) * m(7) + (t193 * t41 + t230 * t28 + t272 + t4) * m(6); -t192 * t130 + t193 * t132 + (-g(1) * t106 - g(2) * t104 + g(3) * t173 + t15 * t193 - t18 * t192 + t5) * m(7) + t20;];
tau  = t9;
