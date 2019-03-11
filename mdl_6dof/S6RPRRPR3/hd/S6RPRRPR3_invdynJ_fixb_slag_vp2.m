% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:42
% EndTime: 2019-03-09 05:05:21
% DurationCPUTime: 25.45s
% Computational Cost: add. (7470->768), mult. (15099->972), div. (0->0), fcn. (9500->12), ass. (0->344)
t533 = Ifges(5,1) + Ifges(6,1);
t532 = Ifges(6,4) + Ifges(5,5);
t530 = Ifges(5,6) - Ifges(6,6);
t250 = cos(qJ(3));
t541 = t250 / 0.2e1;
t246 = sin(qJ(3));
t414 = Ifges(4,4) * t246;
t540 = t414 / 0.2e1;
t242 = sin(pkin(10));
t222 = pkin(1) * t242 + pkin(7);
t198 = t222 * qJD(1);
t373 = qJD(2) * t246;
t154 = t198 * t250 + t373;
t136 = qJD(3) * pkin(8) + t154;
t237 = t246 * pkin(8);
t239 = t250 * pkin(3);
t345 = -pkin(2) - t239;
t280 = t345 - t237;
t243 = cos(pkin(10));
t428 = pkin(1) * t243;
t141 = (t280 - t428) * qJD(1);
t245 = sin(qJ(4));
t249 = cos(qJ(4));
t66 = -t245 * t136 + t249 * t141;
t539 = qJD(5) - t66;
t374 = qJD(1) * t250;
t514 = t374 - qJD(4);
t212 = qJD(6) + t514;
t436 = t212 / 0.2e1;
t363 = t249 * qJD(3);
t375 = qJD(1) * t246;
t181 = t245 * t375 - t363;
t372 = qJD(3) * t245;
t182 = t249 * t375 + t372;
t244 = sin(qJ(6));
t248 = cos(qJ(6));
t283 = t181 * t244 + t182 * t248;
t444 = t283 / 0.2e1;
t106 = t181 * t248 - t182 * t244;
t446 = t106 / 0.2e1;
t538 = Ifges(7,5) * t444 + Ifges(7,6) * t446 + Ifges(7,3) * t436;
t362 = qJD(1) * qJD(3);
t192 = qJDD(1) * t246 + t250 * t362;
t369 = qJD(4) * t181;
t99 = qJDD(3) * t245 + t192 * t249 - t369;
t452 = t99 / 0.2e1;
t100 = qJD(4) * t182 - t249 * qJDD(3) + t192 * t245;
t448 = t100 / 0.2e1;
t191 = t250 * qJDD(1) - t246 * t362;
t177 = qJDD(4) - t191;
t442 = t177 / 0.2e1;
t440 = t181 / 0.2e1;
t537 = -t182 / 0.2e1;
t536 = t514 / 0.2e1;
t535 = Ifges(4,2) * t541 + t540;
t534 = mrSges(6,2) + mrSges(5,3);
t531 = Ifges(6,2) + Ifges(5,3);
t344 = t245 * t374;
t368 = qJD(4) * t245;
t450 = pkin(8) - pkin(9);
t153 = t250 * qJD(2) - t246 * t198;
t422 = pkin(8) * t250;
t313 = pkin(3) * t246 - t422;
t187 = t313 * qJD(1);
t92 = t249 * t153 + t245 * t187;
t76 = qJ(5) * t375 + t92;
t529 = pkin(9) * t344 + t450 * t368 + t76;
t204 = t450 * t249;
t451 = pkin(4) + pkin(5);
t353 = t451 * t246;
t381 = t249 * t250;
t360 = pkin(9) * t381;
t91 = -t245 * t153 + t187 * t249;
t528 = qJD(4) * t204 - (-t353 - t360) * qJD(1) + t91;
t527 = -pkin(9) * t182 + t539;
t485 = -t245 * t530 + t249 * t532;
t409 = Ifges(6,5) * t245;
t412 = Ifges(5,4) * t245;
t483 = t249 * t533 + t409 - t412;
t526 = m(7) * pkin(9) + mrSges(7,3) - t534;
t307 = t249 * mrSges(6,1) + t245 * mrSges(6,3);
t309 = mrSges(5,1) * t249 - mrSges(5,2) * t245;
t281 = t244 * t245 + t248 * t249;
t384 = t245 * t248;
t282 = t244 * t249 - t384;
t512 = t281 * mrSges(7,1) - t282 * mrSges(7,2);
t525 = t307 + t309 + t512;
t241 = qJ(1) + pkin(10);
t229 = sin(t241);
t230 = cos(t241);
t383 = t245 * t250;
t144 = t229 * t383 + t230 * t249;
t145 = t229 * t381 - t230 * t245;
t284 = t144 * t244 + t145 * t248;
t478 = t144 * t248 - t145 * t244;
t524 = mrSges(7,1) * t478 - t284 * mrSges(7,2);
t449 = -t100 / 0.2e1;
t523 = t362 / 0.2e1;
t522 = t533 * t452 + t532 * t442 + (-Ifges(5,4) + Ifges(6,5)) * t448;
t473 = qJD(4) - qJD(6);
t521 = t282 * t473;
t520 = qJD(2) * qJD(3) + t222 * qJDD(1);
t366 = qJD(4) * t249;
t371 = qJD(3) * t246;
t83 = t246 * qJDD(2) - t198 * t371 + t250 * t520;
t79 = qJDD(3) * pkin(8) + t83;
t223 = -pkin(2) - t428;
t197 = t223 * qJDD(1);
t98 = -pkin(3) * t191 - pkin(8) * t192 + t197;
t15 = -t136 * t368 + t141 * t366 + t245 * t98 + t249 * t79;
t16 = -t136 * t366 - t141 * t368 - t245 * t79 + t249 * t98;
t287 = t15 * t249 - t16 * t245;
t67 = t249 * t136 + t245 * t141;
t519 = -t66 * t366 - t67 * t368 + t287;
t11 = t177 * qJ(5) - qJD(5) * t514 + t15;
t273 = qJDD(5) - t16;
t14 = -pkin(4) * t177 + t273;
t289 = t11 * t249 + t14 * t245;
t54 = pkin(4) * t514 + t539;
t210 = t514 * qJ(5);
t57 = -t210 + t67;
t518 = t54 * t366 - t57 * t368 + t289;
t517 = qJD(5) * t245 + t373;
t59 = mrSges(5,1) * t177 - mrSges(5,3) * t99;
t60 = -t177 * mrSges(6,1) + t99 * mrSges(6,2);
t61 = -mrSges(5,2) * t177 - mrSges(5,3) * t100;
t62 = -mrSges(6,2) * t100 + mrSges(6,3) * t177;
t516 = (t61 + t62) * t249 + (-t59 + t60) * t245;
t228 = Ifges(4,4) * t374;
t174 = Ifges(5,4) * t181;
t402 = t181 * Ifges(6,5);
t498 = t182 * t533 - t514 * t532 - t174 + t402;
t173 = Ifges(6,5) * t182;
t85 = -Ifges(6,6) * t514 + t181 * Ifges(6,3) + t173;
t515 = Ifges(4,1) * t375 + Ifges(4,5) * qJD(3) + t245 * t85 + t249 * t498 + t228;
t36 = t451 * t514 + t527;
t49 = pkin(9) * t181 + t67;
t43 = -t210 + t49;
t12 = -t244 * t43 + t248 * t36;
t8 = -pkin(9) * t99 - t177 * t451 + t273;
t9 = pkin(9) * t100 + t11;
t1 = qJD(6) * t12 + t244 * t8 + t248 * t9;
t13 = t244 * t36 + t248 * t43;
t2 = -qJD(6) * t13 - t244 * t9 + t248 * t8;
t513 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t511 = -t16 * mrSges(5,1) + t14 * mrSges(6,1) + t15 * mrSges(5,2) - t11 * mrSges(6,3);
t510 = t12 * mrSges(7,1) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t535 + t531 * t536 + t532 * t537 + t530 * t440 - t13 * mrSges(7,2) + t538;
t169 = qJDD(6) - t177;
t26 = qJD(6) * t106 + t100 * t244 + t248 * t99;
t27 = -qJD(6) * t283 + t100 * t248 - t244 * t99;
t509 = -t26 * Ifges(7,4) / 0.2e1 - t27 * Ifges(7,2) / 0.2e1 - t169 * Ifges(7,6) / 0.2e1;
t459 = t26 / 0.2e1;
t458 = t27 / 0.2e1;
t508 = -m(6) - m(7);
t443 = t169 / 0.2e1;
t420 = -qJD(1) / 0.2e1;
t505 = -mrSges(4,3) + mrSges(3,2);
t504 = mrSges(6,3) - mrSges(5,2);
t203 = t450 * t245;
t119 = t203 * t244 + t204 * t248;
t502 = -qJD(6) * t119 + t244 * t529 + t248 * t528;
t118 = t203 * t248 - t204 * t244;
t501 = qJD(6) * t118 + t244 * t528 - t248 * t529;
t354 = t451 * t245;
t390 = qJ(5) * t249;
t279 = -t354 + t390;
t496 = qJD(4) * t279 - (qJD(1) * t279 - t198) * t250 + t517;
t193 = -qJ(5) * t244 - t248 * t451;
t495 = qJD(6) * t193 - t244 * t49 + t248 * t527;
t194 = qJ(5) * t248 - t244 * t451;
t494 = -qJD(6) * t194 - t244 * t527 - t248 * t49;
t41 = mrSges(5,1) * t100 + mrSges(5,2) * t99;
t493 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t192 + t41;
t311 = mrSges(4,1) * t250 - mrSges(4,2) * t246;
t491 = -t311 - mrSges(3,1);
t352 = mrSges(4,3) * t375;
t490 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t181 + mrSges(5,2) * t182 + t352;
t416 = mrSges(5,3) * t181;
t127 = mrSges(5,2) * t514 - t416;
t418 = mrSges(6,2) * t181;
t130 = -mrSges(6,3) * t514 - t418;
t380 = t127 + t130;
t415 = mrSges(5,3) * t182;
t128 = -mrSges(5,1) * t514 - t415;
t417 = mrSges(6,2) * t182;
t129 = mrSges(6,1) * t514 + t417;
t379 = -t128 + t129;
t382 = t246 * t249;
t155 = t244 * t382 - t246 * t384;
t156 = t281 * t246;
t323 = t155 * mrSges(7,1) + t156 * mrSges(7,2);
t426 = pkin(4) * t245;
t290 = -t390 + t426;
t489 = qJD(4) * t290 - (qJD(1) * t290 + t198) * t250 - t517;
t488 = t246 * t531 + t250 * t485;
t487 = t246 * t532 + t250 * t483;
t484 = t245 * t532 + t249 * t530;
t408 = Ifges(6,5) * t249;
t411 = Ifges(5,4) * t249;
t482 = t245 * t533 - t408 + t411;
t479 = -t100 * t530 + t177 * t531 + t532 * t99;
t370 = qJD(3) * t250;
t84 = t250 * qJDD(2) - t198 * t370 - t246 * t520;
t477 = -t246 * t84 + t250 * t83;
t476 = t534 * t246;
t474 = -m(5) + t508;
t472 = -t99 * Ifges(6,5) / 0.2e1 - t177 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t452 + Ifges(5,6) * t442 + (Ifges(6,3) + Ifges(5,2)) * t449;
t315 = qJD(3) * pkin(3) + t153;
t269 = qJ(5) * t182 + t315;
t47 = -t181 * t451 + t269;
t467 = -t47 * mrSges(7,1) + mrSges(7,3) * t13;
t466 = t47 * mrSges(7,2) - mrSges(7,3) * t12;
t391 = qJ(5) * t245;
t465 = -t249 * t451 - t391;
t464 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1);
t291 = pkin(4) * t249 + t391;
t195 = -pkin(3) - t291;
t310 = mrSges(4,1) * t246 + mrSges(4,2) * t250;
t425 = pkin(5) * t249;
t463 = t310 + t526 * t250 + (-m(7) * (t195 - t425) - m(6) * t195 + m(5) * pkin(3) + t525) * t246;
t462 = m(5) * t315 - t490;
t460 = Ifges(7,1) * t459 + Ifges(7,4) * t458 + Ifges(7,5) * t443;
t410 = Ifges(7,4) * t283;
t38 = Ifges(7,2) * t106 + Ifges(7,6) * t212 + t410;
t457 = -t38 / 0.2e1;
t456 = t38 / 0.2e1;
t104 = Ifges(7,4) * t106;
t39 = Ifges(7,1) * t283 + Ifges(7,5) * t212 + t104;
t455 = -t39 / 0.2e1;
t454 = t39 / 0.2e1;
t401 = t182 * Ifges(5,4);
t88 = -t181 * Ifges(5,2) - Ifges(5,6) * t514 + t401;
t453 = -t88 / 0.2e1;
t447 = -t106 / 0.2e1;
t445 = -t283 / 0.2e1;
t441 = -t181 / 0.2e1;
t438 = t182 / 0.2e1;
t437 = -t212 / 0.2e1;
t247 = sin(qJ(1));
t427 = pkin(1) * t247;
t421 = pkin(9) * t246;
t251 = cos(qJ(1));
t240 = t251 * pkin(1);
t419 = qJD(4) / 0.2e1;
t413 = Ifges(4,4) * t250;
t396 = t249 * mrSges(6,3);
t115 = mrSges(6,1) * t181 - mrSges(6,3) * t182;
t44 = -mrSges(7,1) * t106 + mrSges(7,2) * t283;
t394 = t115 - t44;
t392 = qJ(5) * t181;
t387 = t230 * t246;
t386 = t230 * t250;
t385 = t245 * t246;
t377 = t239 + t237;
t171 = t223 - t377;
t190 = t313 * qJD(3);
t378 = t171 * t366 + t245 * t190;
t186 = t222 * t381;
t111 = t245 * t171 + t186;
t367 = qJD(4) * t246;
t364 = qJD(5) * t249;
t359 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t169;
t356 = pkin(8) * t368;
t355 = pkin(8) * t366;
t351 = mrSges(4,3) * t374;
t346 = t230 * pkin(2) + t229 * pkin(7) + t240;
t343 = t222 * t371;
t341 = t245 * t370;
t340 = t245 * t367;
t334 = -t374 / 0.2e1;
t327 = t366 / 0.2e1;
t326 = t230 * pkin(7) - t427;
t325 = -t222 * t245 - pkin(4);
t185 = t222 * t383;
t110 = t171 * t249 - t185;
t319 = pkin(4) * t381 + qJ(5) * t383 + t377;
t146 = -t229 * t249 + t230 * t383;
t147 = t229 * t245 + t230 * t381;
t68 = t146 * t248 - t147 * t244;
t69 = t146 * t244 + t147 * t248;
t314 = mrSges(7,1) * t68 - mrSges(7,2) * t69;
t102 = -qJ(5) * t250 + t111;
t312 = qJD(4) * t186 + t171 * t368 - t190 * t249;
t308 = mrSges(5,1) * t245 + mrSges(5,2) * t249;
t306 = mrSges(6,1) * t245 - t396;
t300 = -Ifges(5,2) * t245 + t411;
t299 = Ifges(5,2) * t249 + t412;
t296 = Ifges(4,5) * t250 - Ifges(4,6) * t246;
t293 = Ifges(6,3) * t245 + t408;
t292 = -Ifges(6,3) * t249 + t409;
t238 = t250 * pkin(4);
t70 = pkin(5) * t250 + t185 + t238 + (-t171 - t421) * t249;
t81 = pkin(9) * t385 + t102;
t28 = -t244 * t81 + t248 * t70;
t29 = t244 * t70 + t248 * t81;
t73 = -mrSges(7,2) * t212 + mrSges(7,3) * t106;
t74 = mrSges(7,1) * t212 - mrSges(7,3) * t283;
t286 = -t244 * t74 + t248 * t73;
t285 = pkin(3) * t386 + pkin(8) * t387 + t346;
t276 = t223 * qJD(1) * t310;
t275 = t246 * (Ifges(4,1) * t250 - t414);
t274 = qJDD(3) * pkin(3) + t84;
t272 = t250 * t281;
t271 = t282 * t250;
t268 = t359 + t513;
t265 = -g(1) * t146 - g(2) * t144 - g(3) * t385;
t260 = Ifges(5,6) * t246 + t250 * t300;
t259 = Ifges(6,6) * t246 + t250 * t293;
t50 = (-t246 * t363 - t250 * t368) * t222 + t378;
t256 = qJ(5) * t99 + qJD(5) * t182 + t274;
t112 = t473 * t281;
t227 = qJ(5) * t371;
t215 = qJ(5) * t382;
t201 = -qJD(3) * mrSges(4,2) + t351;
t172 = pkin(3) - t465;
t170 = t308 * t246;
t158 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t191;
t143 = qJD(1) * t272;
t142 = qJD(1) * t271;
t124 = -t215 + (t222 + t426) * t246;
t114 = pkin(4) * t182 + t392;
t109 = t215 + (-t222 - t354) * t246;
t103 = -t110 + t238;
t78 = -pkin(4) * t375 - t91;
t75 = -t182 * t451 - t392;
t65 = pkin(4) * t181 - t269;
t63 = (qJD(4) * t291 - t364) * t246 + (t222 + t290) * t370;
t53 = qJD(3) * t272 + t246 * t521;
t52 = -qJD(3) * t271 + t112 * t246;
t51 = t245 * t343 - t312;
t46 = t325 * t371 + t312;
t45 = (qJD(4) * t465 + t364) * t246 + (-t222 + t279) * t370;
t42 = -qJD(5) * t250 + t227 + t50;
t40 = mrSges(6,1) * t100 - mrSges(6,3) * t99;
t35 = t227 + (pkin(9) * qJD(4) - qJD(3) * t222) * t382 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t222) * t245) * t250 + t378;
t34 = pkin(9) * t340 + (-t360 + (-pkin(5) + t325) * t246) * qJD(3) + t312;
t21 = pkin(4) * t100 - t256;
t20 = -mrSges(7,2) * t169 + mrSges(7,3) * t27;
t19 = mrSges(7,1) * t169 - mrSges(7,3) * t26;
t10 = -t100 * t451 + t256;
t7 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t4 = -qJD(6) * t29 - t244 * t35 + t248 * t34;
t3 = qJD(6) * t28 + t244 * t34 + t248 * t35;
t5 = [t515 * t370 / 0.2e1 + (-t11 * mrSges(6,2) - t15 * mrSges(5,3) - t472) * t385 + (Ifges(7,4) * t156 - Ifges(7,2) * t155) * t458 + (Ifges(7,4) * t53 + Ifges(7,2) * t52) * t446 + (t413 / 0.2e1 + t223 * mrSges(4,2)) * t192 + (Ifges(7,5) * t53 + Ifges(7,6) * t52) * t436 + (Ifges(7,5) * t156 - Ifges(7,6) * t155) * t443 + (-t479 / 0.2e1 - t532 * t452 - t531 * t442 - Ifges(6,6) * t448 - Ifges(5,6) * t449 + Ifges(7,6) * t458 + Ifges(7,5) * t459 + (-Ifges(4,2) * t246 + t413) * t523 + Ifges(7,3) * t443 + t511 + t513) * t250 + t275 * t523 + (-t153 * t370 + t477) * mrSges(4,3) + (qJD(3) * t487 - t367 * t482) * t438 - (qJD(3) * t488 - t367 * t484) * t514 / 0.2e1 + (m(3) * t427 + mrSges(2,1) * t247 + t284 * mrSges(7,1) + mrSges(2,2) * t251 + t478 * mrSges(7,2) + (-m(5) - m(4)) * t326 + t508 * (-t145 * pkin(4) - qJ(5) * t144 + t326) + t505 * t230 + t464 * t145 + t504 * t144 + (-m(7) * t345 + m(4) * pkin(2) + (-m(6) - m(5)) * t280 + t476 - t491) * t229) * g(1) + (Ifges(7,1) * t156 - Ifges(7,4) * t155) * t459 + (Ifges(7,1) * t53 + Ifges(7,4) * t52) * t444 + t10 * t323 + (t66 * mrSges(5,1) - t54 * mrSges(6,1) - t67 * mrSges(5,2) - t154 * mrSges(4,3) + t57 * mrSges(6,3) - t510 - t538) * t371 - (t245 * t498 + t249 * t88) * t367 / 0.2e1 - t274 * t170 + (-mrSges(5,1) * t315 + mrSges(6,1) * t65 - mrSges(6,2) * t57 - mrSges(5,3) * t67) * (t246 * t366 + t341) + (-mrSges(5,2) * t315 + mrSges(6,2) * t54 - mrSges(5,3) * t66 - mrSges(6,3) * t65) * (t250 * t363 - t340) + (m(4) * ((-t153 * t250 - t154 * t246) * qJD(3) + t477) + t490 * t370 + t493 * t246 + m(5) * (-t246 * t274 - t315 * t370) + t250 * t158) * t222 + (-t223 * mrSges(4,1) + t535 + t540) * t191 + (Ifges(4,4) * t192 + Ifges(4,2) * t191 + 0.2e1 * Ifges(4,6) * qJDD(3) + t359) * t541 + (-t1 * t155 - t12 * t53 + t13 * t52 - t156 * t2) * mrSges(7,3) + qJD(3) ^ 2 * t296 / 0.2e1 + m(7) * (t1 * t29 + t10 * t109 + t12 * t4 + t13 * t3 + t2 * t28 + t45 * t47) + m(6) * (t102 * t11 + t103 * t14 + t124 * t21 + t42 * t57 + t46 * t54 + t63 * t65) + (t14 * mrSges(6,2) - t16 * mrSges(5,3) + t522) * t382 + t341 * t453 + t53 * t454 + t52 * t456 + t156 * t460 + qJD(3) * t276 + (t483 * t452 + t485 * t442 + t293 * t448 + t300 * t449 + Ifges(4,1) * t192 + (m(7) * t450 - mrSges(7,3)) * t229 * g(1) + t85 * t327 + Ifges(4,5) * qJDD(3) + t21 * t306) * t246 + (-m(3) * t240 - m(4) * t346 - m(5) * t285 - mrSges(2,1) * t251 - t69 * mrSges(7,1) + mrSges(2,2) * t247 - t68 * mrSges(7,2) + t508 * (t147 * pkin(4) + qJ(5) * t146 + t285) + t491 * t230 + t505 * t229 - t464 * t147 - t504 * t146 + t526 * t387) * g(2) - t201 * t343 + m(5) * (t110 * t16 + t111 * t15 + t50 * t67 + t51 * t66) + (qJD(3) * t259 - t292 * t367) * t440 + (qJD(3) * t260 - t299 * t367) * t441 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t243 - 0.2e1 * mrSges(3,2) * t242 + m(3) * (t242 ^ 2 + t243 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t124 * t40 + t50 * t127 + t51 * t128 + t46 * t129 + t42 * t130 + t111 * t61 + t63 * t115 + t109 * t7 + t110 * t59 + t103 * t60 + t102 * t62 + t3 * t73 + t4 * t74 + (m(4) * t223 - t311) * t197 + t155 * t509 + t28 * t19 + t29 * t20 + t45 * t44 + t47 * (-mrSges(7,1) * t52 + mrSges(7,2) * t53); m(3) * qJDD(2) - t155 * t19 + t156 * t20 + t53 * t73 + t52 * t74 + m(7) * (t1 * t156 + t12 * t52 + t13 * t53 - t155 * t2) + (-m(3) - m(4) + t474) * g(3) + (-t40 + t7 + (t245 * t379 + t249 * t380 + t201) * qJD(3) + m(4) * (qJD(3) * t154 + t84) + m(7) * t10 + m(6) * (t363 * t57 + t372 * t54 - t21) + m(5) * (t363 * t67 - t372 * t66 + t274) - t493) * t250 + (t158 + (-t245 * t380 + t249 * t379) * qJD(4) + m(4) * t83 + m(6) * t518 + m(5) * t519 + (-m(4) * t153 + m(6) * t65 - m(7) * t47 + t394 - t462) * qJD(3) + t516) * t246; (t228 + t515) * t334 + (t386 * t474 * g(1) + ((-t245 * t57 + t249 * t54) * qJD(4) + t289) * m(6) + ((-t245 * t67 - t249 * t66) * qJD(4) + t287) * m(5) + t516) * pkin(8) + t518 * mrSges(6,2) + t472 * t249 + (Ifges(7,1) * t444 + Ifges(7,4) * t446 + Ifges(7,5) * t436 + t454 + t466) * t112 + (-t76 - t356) * t130 + t482 * t452 + (t195 * t21 + t489 * t65 - t54 * t78 - t57 * t76) * m(6) + (t355 - t78) * t129 + (-t91 - t355) * t128 + (-t92 - t356) * t127 + t10 * t512 + (-Ifges(7,5) * t445 - Ifges(4,2) * t334 - Ifges(7,6) * t447 - Ifges(7,3) * t437 + t510) * t375 + t245 * t522 + (t85 / 0.2e1 + t453) * t368 + t498 * t327 + t501 * t73 + (t1 * t119 + t10 * t172 + t118 * t2 + t12 * t502 + t13 * t501 + t47 * t496) * m(7) + t502 * t74 + (t351 - t201) * t153 + t484 * t442 + t489 * t115 + (-t300 / 0.2e1 + t293 / 0.2e1) * t369 + (t422 * t474 + t463) * g(2) * t229 + t496 * t44 + t463 * t230 * g(1) + t88 * t344 / 0.2e1 + (pkin(3) * t274 - t66 * t91 - t67 * t92) * m(5) + t274 * t309 + t281 * t509 + (-Ifges(7,4) * t282 - Ifges(7,2) * t281) * t458 + (-Ifges(7,1) * t282 - Ifges(7,4) * t281) * t459 + (-Ifges(7,5) * t282 - Ifges(7,6) * t281) * t443 + ((t260 / 0.2e1 - t259 / 0.2e1) * t181 - t66 * (mrSges(5,1) * t246 - mrSges(5,3) * t381) - t54 * (-mrSges(6,1) * t246 + mrSges(6,2) * t381) - t67 * (-mrSges(5,2) * t246 - mrSges(5,3) * t383) - t57 * (-mrSges(6,2) * t383 + mrSges(6,3) * t246) - t276 + t275 * t420) * qJD(1) + (t352 + t462) * t154 + t292 * t448 + t299 * t449 + t143 * t455 - t142 * t457 + (-m(6) * t319 - m(5) * t377 - m(7) * (t319 - t421) + t246 * mrSges(7,3) - t311 + (-m(7) * t425 - t525) * t250 - t476) * g(3) + (-t65 * t306 + t308 * t315 - t419 * t485 - t420 * t488) * t514 - (Ifges(7,4) * t444 + Ifges(7,2) * t446 + Ifges(7,6) * t436 + t456 + t467) * t521 + t519 * mrSges(5,3) - t21 * t307 - t296 * t362 / 0.2e1 + t195 * t40 + Ifges(4,6) * t191 + Ifges(4,5) * t192 + Ifges(4,3) * qJDD(3) + t172 * t7 + t118 * t19 + t119 * t20 - t83 * mrSges(4,2) + t84 * mrSges(4,1) - t282 * t460 + (Ifges(7,4) * t143 - Ifges(7,2) * t142) * t447 + (Ifges(7,5) * t143 - Ifges(7,6) * t142) * t437 + (-t1 * t281 + t12 * t143 + t13 * t142 + t2 * t282) * mrSges(7,3) + (Ifges(7,1) * t143 - Ifges(7,4) * t142) * t445 - t47 * (mrSges(7,1) * t142 + mrSges(7,2) * t143) + (t419 * t483 + t420 * t487) * t182 - pkin(3) * t41; -(Ifges(7,1) * t445 + Ifges(7,4) * t447 + Ifges(7,5) * t437 + t455 - t466) * t106 + (Ifges(7,4) * t445 + Ifges(7,2) * t447 + Ifges(7,6) * t437 + t457 - t467) * t283 + t479 - t511 + (-t181 * t532 - t182 * t530) * t536 + (t508 * (-t144 * pkin(4) + qJ(5) * t145) - t504 * t145 + t464 * t144 + t524) * g(2) + (-Ifges(5,2) * t182 - t174 + t498) * t440 + (-m(6) * t54 - t379 + t415) * t67 + (-m(6) * t57 - t380 - t416) * t66 - t268 + t57 * t417 + t54 * t418 + t494 * t74 + t495 * t73 + (t1 * t194 + t12 * t494 + t13 * t495 + t193 * t2 - t47 * t75) * m(7) + (t314 + t508 * (-t146 * pkin(4) + qJ(5) * t147) - t504 * t147 + t464 * t146) * g(1) + t315 * (mrSges(5,1) * t182 - mrSges(5,2) * t181) + (-t181 * t533 + t173 - t401 + t85) * t537 + (-pkin(4) * t14 + qJ(5) * t11 + qJD(5) * t57 - t114 * t65) * m(6) + (-m(6) * t215 - (t396 + (-m(6) * pkin(4) - mrSges(6,1)) * t245) * t246 - m(7) * (-t245 * t353 + t215) - t323 + t170) * g(3) + t88 * t438 + (Ifges(6,3) * t182 - t402) * t441 + t194 * t20 + t193 * t19 - t65 * (mrSges(6,1) * t182 + mrSges(6,3) * t181) + qJD(5) * t130 - t114 * t115 - t75 * t44 - pkin(4) * t60 + qJ(5) * t62; t248 * t19 + t244 * t20 + t394 * t182 + t286 * qJD(6) - (-t130 - t286) * t514 + t60 + (t1 * t244 - t182 * t47 + t2 * t248 + t265 + t212 * (-t12 * t244 + t13 * t248)) * m(7) + (t182 * t65 + t514 * t57 + t14 + t265) * m(6); -t47 * (mrSges(7,1) * t283 + mrSges(7,2) * t106) + (Ifges(7,1) * t106 - t410) * t445 + t38 * t444 + (Ifges(7,5) * t106 - Ifges(7,6) * t283) * t437 - t12 * t73 + t13 * t74 - g(1) * t314 - g(2) * t524 + g(3) * t323 + (t106 * t12 + t13 * t283) * mrSges(7,3) + t268 + (-Ifges(7,2) * t283 + t104 + t39) * t447;];
tau  = t5;
