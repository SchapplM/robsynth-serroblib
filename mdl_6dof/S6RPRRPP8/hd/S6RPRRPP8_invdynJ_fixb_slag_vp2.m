% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:54
% EndTime: 2019-03-09 04:54:28
% DurationCPUTime: 21.82s
% Computational Cost: add. (4535->684), mult. (8824->821), div. (0->0), fcn. (4923->6), ass. (0->308)
t485 = Ifges(5,1) + Ifges(7,3);
t458 = Ifges(5,5) + Ifges(7,5);
t193 = cos(qJ(4));
t190 = sin(qJ(4));
t335 = qJD(3) * t190;
t194 = cos(qJ(3));
t336 = qJD(1) * t194;
t136 = t193 * t336 + t335;
t127 = Ifges(7,6) * t136;
t324 = t193 * qJD(3);
t135 = t190 * t336 - t324;
t191 = sin(qJ(3));
t337 = qJD(1) * t191;
t169 = qJD(4) + t337;
t369 = t136 * Ifges(6,6);
t459 = Ifges(7,4) + Ifges(6,5);
t484 = Ifges(7,2) + Ifges(6,3);
t454 = t135 * t484 + t169 * t459 + t127 - t369;
t370 = t136 * Ifges(5,4);
t49 = -t135 * Ifges(5,2) + t169 * Ifges(5,6) + t370;
t494 = -t49 / 0.2e1 + t454 / 0.2e1;
t492 = -t193 / 0.2e1;
t129 = Ifges(5,4) * t135;
t373 = Ifges(7,6) * t135;
t453 = t136 * t485 + t169 * t458 - t129 + t373;
t491 = t453 / 0.2e1;
t489 = qJD(1) / 0.2e1;
t488 = qJD(3) / 0.2e1;
t481 = Ifges(6,4) - t458;
t480 = -Ifges(5,6) + t459;
t244 = Ifges(6,4) * t193 - Ifges(6,5) * t190;
t439 = t458 * t193 + (Ifges(7,4) - Ifges(5,6)) * t190;
t487 = t244 - t439;
t322 = qJD(1) * qJD(3);
t221 = -qJDD(1) * t194 + t191 * t322;
t327 = qJD(4) * t194;
t301 = t190 * t327;
t65 = qJD(1) * t301 - qJD(4) * t324 - t190 * qJDD(3) + t193 * t221;
t416 = -t65 / 0.2e1;
t330 = qJD(4) * t136;
t66 = -t193 * qJDD(3) - t190 * t221 + t330;
t478 = t66 / 0.2e1;
t144 = -qJDD(1) * t191 - t194 * t322;
t133 = qJDD(4) - t144;
t410 = t133 / 0.2e1;
t409 = -t135 / 0.2e1;
t406 = t136 / 0.2e1;
t405 = -t169 / 0.2e1;
t196 = -pkin(1) - pkin(7);
t486 = m(5) * t196;
t372 = Ifges(7,6) * t190;
t377 = Ifges(5,4) * t190;
t441 = t485 * t193 + t372 - t377;
t374 = Ifges(6,6) * t193;
t233 = -Ifges(6,3) * t190 + t374;
t371 = Ifges(7,6) * t193;
t234 = Ifges(7,2) * t190 + t371;
t483 = t233 - t234;
t482 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t379 = Ifges(4,4) * t191;
t250 = t194 * Ifges(4,1) - t379;
t128 = Ifges(6,6) * t135;
t48 = t169 * Ifges(6,4) - t136 * Ifges(6,2) + t128;
t479 = Ifges(4,5) * t488 + t494 * t190 + t193 * t491 + t250 * t489 + t48 * t492;
t477 = t144 / 0.2e1;
t476 = -t221 / 0.2e1;
t186 = t194 * pkin(8);
t146 = pkin(3) * t191 + qJ(2) - t186;
t218 = t146 * qJD(1);
t158 = qJDD(1) * t196 + qJDD(2);
t165 = qJD(1) * t196 + qJD(2);
t333 = qJD(3) * t194;
t86 = t191 * t158 + t165 * t333;
t475 = qJDD(3) * pkin(8) + qJD(4) * t218 + t86;
t145 = t191 * t165;
t115 = qJD(3) * pkin(8) + t145;
t329 = qJD(4) * t190;
t323 = qJD(1) * qJD(2);
t166 = qJDD(1) * qJ(2) + t323;
t68 = -t144 * pkin(3) + pkin(8) * t221 + t166;
t7 = -t115 * t329 + t190 * t68 + t193 * t475;
t328 = qJD(4) * t193;
t8 = -t115 * t328 - t190 * t475 + t193 * t68;
t258 = -t190 * t8 + t193 * t7;
t53 = t115 * t190 - t193 * t218;
t54 = t193 * t115 + t190 * t218;
t474 = t53 * t328 - t54 * t329 + t258;
t4 = -qJ(5) * t133 - qJD(5) * t169 - t7;
t225 = qJDD(5) - t8;
t5 = -pkin(4) * t133 + t225;
t259 = t190 * t5 - t193 * t4;
t38 = -pkin(4) * t169 + qJD(5) + t53;
t39 = -t169 * qJ(5) - t54;
t473 = t38 * t328 + t39 * t329 + t259;
t359 = qJ(5) * t193;
t472 = pkin(4) * t329 - qJD(5) * t190 - t337 * t359 - t145;
t471 = t191 * t324 + t301;
t470 = -t135 * pkin(5) + qJD(6);
t226 = pkin(5) * t136 + t53;
t469 = t226 + qJD(5);
t414 = -t66 / 0.2e1;
t468 = -t133 * Ifges(6,4) / 0.2e1 + Ifges(6,6) * t414 + (-Ifges(5,4) + Ifges(7,6)) * t478 + t458 * t410 + (Ifges(6,2) + t485) * t416;
t253 = t193 * mrSges(6,2) - t190 * mrSges(6,3);
t255 = mrSges(5,1) * t193 - mrSges(5,2) * t190;
t368 = t190 * mrSges(7,2);
t466 = t255 + t368 - t253;
t378 = Ifges(4,4) * t194;
t247 = -t191 * Ifges(4,2) + t378;
t465 = Ifges(4,6) * t488 + t247 * t489 + t482 * t405 + t481 * t406 + t480 * t409;
t388 = pkin(4) + qJ(6);
t1 = -pkin(5) * t65 - qJD(6) * t169 - t133 * t388 + t225;
t2 = -pkin(5) * t66 + qJDD(6) - t4;
t464 = t8 * mrSges(5,1) - t7 * mrSges(5,2) + t5 * mrSges(6,2) + t2 * mrSges(7,2) - t4 * mrSges(6,3) - t1 * mrSges(7,3);
t463 = -m(7) - m(6);
t461 = m(6) * (t190 * t39 + t193 * t38);
t460 = mrSges(5,3) + mrSges(6,1);
t457 = t484 * t66 + (Ifges(6,6) - Ifges(7,6)) * t65 + t459 * t133;
t28 = -t65 * mrSges(7,1) - t133 * mrSges(7,3);
t31 = -t65 * mrSges(6,1) + t133 * mrSges(6,2);
t456 = t28 + t31;
t29 = mrSges(6,1) * t66 - mrSges(6,3) * t133;
t30 = -t66 * mrSges(7,1) + t133 * mrSges(7,2);
t455 = -t29 + t30;
t295 = t388 * t191;
t358 = qJ(6) * t190;
t452 = qJD(1) * t190 * t295 - qJD(6) * t193 + (t358 - t359) * qJD(4) + t472;
t381 = mrSges(5,3) * t135;
t87 = -mrSges(5,2) * t169 - t381;
t385 = mrSges(6,1) * t135;
t90 = -mrSges(6,3) * t169 + t385;
t451 = t90 - t87;
t383 = mrSges(7,1) * t135;
t91 = mrSges(7,2) * t169 - t383;
t386 = t90 - t91;
t380 = mrSges(5,3) * t136;
t88 = mrSges(5,1) * t169 - t380;
t384 = mrSges(6,1) * t136;
t92 = mrSges(6,2) * t169 + t384;
t450 = t92 - t88;
t18 = mrSges(5,1) * t66 - mrSges(5,2) * t65;
t449 = qJDD(3) * mrSges(4,1) + mrSges(4,3) * t221 - t18;
t411 = pkin(5) + pkin(8);
t157 = t411 * t193;
t353 = t191 * t193;
t219 = -pkin(5) * t353 - t194 * t388;
t260 = pkin(3) * t194 + pkin(8) * t191;
t140 = t260 * qJD(1);
t355 = t190 * t194;
t72 = t140 * t193 - t165 * t355;
t448 = -qJD(1) * t219 + qJD(4) * t157 + t72;
t356 = t190 * t191;
t347 = t193 * t194;
t73 = t190 * t140 + t165 * t347;
t447 = -(pkin(5) * t356 + qJ(5) * t194) * qJD(1) - t73 - t411 * t329;
t446 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t135 + mrSges(5,2) * t136 + mrSges(4,3) * t336;
t445 = pkin(4) * t190 * t337 - qJ(5) * t328 + t472;
t192 = sin(qJ(1));
t354 = t191 * t192;
t349 = t192 * t193;
t195 = cos(qJ(1));
t352 = t191 * t195;
t112 = t190 * t352 + t349;
t344 = t195 * t193;
t350 = t192 * t190;
t113 = t191 * t344 - t350;
t444 = t113 * pkin(4) + qJ(5) * t112;
t443 = -t191 * t441 + t194 * t458;
t442 = t191 * t483 + t194 * t459;
t375 = Ifges(6,6) * t190;
t440 = t193 * t484 - t372 + t375;
t436 = -t54 - t470;
t334 = qJD(3) * t191;
t85 = t158 * t194 - t165 * t334;
t434 = -t191 * t86 - t194 * t85;
t433 = -m(5) + t463;
t257 = mrSges(4,1) * t194 - mrSges(4,2) * t191;
t432 = qJ(2) * t257 + t194 * (-Ifges(4,1) * t191 - t378) / 0.2e1;
t431 = mrSges(7,2) + mrSges(6,3) - mrSges(5,2);
t430 = -mrSges(4,3) + mrSges(3,2) - mrSges(2,1);
t428 = g(1) * t354;
t360 = qJ(5) * t190;
t277 = -pkin(3) - t360;
t147 = -pkin(4) * t193 + t277;
t427 = -m(6) * t147 + m(5) * pkin(3) - m(7) * t277 - (-m(7) * t388 - mrSges(7,3)) * t193 + t466;
t426 = t191 * t487 + t482 * t194;
t425 = t133 * t482 + t480 * t66 + t481 * t65;
t296 = m(7) * qJ(6) + mrSges(7,3);
t256 = mrSges(4,1) * t191 + mrSges(4,2) * t194;
t267 = -m(7) * t411 - mrSges(7,1);
t422 = -t194 * t267 + mrSges(2,2) - mrSges(3,3) - t256;
t21 = -t169 * t388 + t469;
t23 = -t39 + t470;
t421 = t1 * t190 + t193 * t2 + t21 * t328 - t23 * t329;
t420 = mrSges(5,1) - mrSges(6,2) + t296;
t357 = t165 * t194;
t116 = -qJD(3) * pkin(3) - t357;
t215 = -qJ(5) * t136 + t116;
t24 = t135 * t388 + t215;
t365 = t193 * mrSges(7,2);
t251 = -mrSges(7,3) * t190 + t365;
t252 = mrSges(6,2) * t190 + mrSges(6,3) * t193;
t254 = mrSges(5,1) * t190 + mrSges(5,2) * t193;
t42 = pkin(4) * t135 + t215;
t419 = t116 * t254 - t24 * t251 - t252 * t42;
t418 = qJD(1) ^ 2;
t415 = t65 / 0.2e1;
t412 = -m(3) - m(4);
t408 = t135 / 0.2e1;
t407 = -t136 / 0.2e1;
t404 = t169 / 0.2e1;
t400 = pkin(5) * t194;
t397 = g(1) * t192;
t396 = g(2) * t195;
t74 = -mrSges(7,2) * t136 + mrSges(7,3) * t135;
t77 = -mrSges(6,2) * t135 - mrSges(6,3) * t136;
t387 = t74 + t77;
t382 = mrSges(7,1) * t136;
t376 = Ifges(5,4) * t193;
t361 = qJ(5) * t135;
t351 = t191 * t196;
t348 = t192 * t194;
t346 = t194 * t195;
t345 = t194 * t196;
t94 = t190 * t146 + t193 * t351;
t342 = -pkin(4) * t355 + qJ(5) * t347;
t341 = pkin(3) * t348 + pkin(8) * t354;
t340 = t195 * pkin(1) + t192 * qJ(2);
t332 = qJD(3) * t196;
t326 = qJD(4) * t196;
t325 = qJDD(1) * mrSges(3,2);
t320 = t87 - t386;
t89 = -mrSges(7,3) * t169 + t382;
t319 = t89 + t450;
t318 = pkin(8) * t348;
t134 = qJD(3) * t260 + qJD(2);
t302 = t194 * t332;
t306 = t190 * t134 + t146 * t328 + t193 * t302;
t305 = t195 * pkin(7) + t340;
t304 = t190 * t334;
t168 = t191 * t332;
t300 = t191 * t326;
t299 = t193 * t327;
t281 = -t327 / 0.2e1;
t280 = t327 / 0.2e1;
t185 = t195 * qJ(2);
t279 = -pkin(1) * t192 + t185;
t276 = -t322 / 0.2e1;
t274 = (t166 + t323) * qJ(2);
t162 = t190 * t351;
t93 = t146 * t193 - t162;
t269 = pkin(4) * t299 + qJ(5) * t471 + t168;
t268 = pkin(3) * t354 + t305;
t83 = -qJ(5) * t191 - t94;
t248 = Ifges(5,1) * t190 + t376;
t246 = -Ifges(5,2) * t190 + t376;
t245 = Ifges(5,2) * t193 + t377;
t243 = Ifges(6,4) * t190 + Ifges(6,5) * t193;
t242 = Ifges(7,4) * t193 - Ifges(7,5) * t190;
t240 = -Ifges(4,5) * t191 - Ifges(4,6) * t194;
t238 = Ifges(5,5) * t190 + Ifges(5,6) * t193;
t237 = Ifges(6,2) * t193 - t375;
t236 = Ifges(6,2) * t190 + t374;
t231 = -Ifges(7,3) * t190 + t371;
t228 = t190 * t54 - t193 * t53;
t33 = -t146 * t329 - t190 * t302 + (t134 - t300) * t193;
t79 = -qJDD(3) * pkin(3) - t85;
t223 = t191 * (-Ifges(4,2) * t194 - t379);
t110 = t191 * t350 - t344;
t111 = t190 * t195 + t191 * t349;
t214 = t111 * pkin(4) + qJ(5) * t110 + t268;
t213 = -g(1) * t110 + g(2) * t112 - g(3) * t355;
t174 = pkin(3) * t352;
t211 = -pkin(8) * t346 + t192 * t196 + t174 + t185;
t208 = Ifges(6,4) * t194 + t191 * t237;
t203 = Ifges(5,6) * t194 - t191 * t246;
t201 = -t190 * t320 + t193 * t319;
t200 = qJ(5) * t65 - qJD(5) * t136 + t79;
t179 = -qJDD(1) * pkin(1) + qJDD(2);
t156 = t411 * t190;
t152 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t337;
t139 = t256 * qJD(1);
t126 = t254 * t194;
t125 = -t193 * t388 + t277;
t101 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t144;
t95 = -t342 - t345;
t84 = -pkin(4) * t191 - t93;
t82 = (-t196 + t358) * t194 - t342;
t75 = pkin(4) * t136 + t361;
t69 = -pkin(5) * t355 - t83;
t67 = -pkin(4) * t336 - t72;
t64 = -qJ(5) * t336 - t73;
t43 = t162 + (-t146 + t400) * t193 - t295;
t40 = t136 * t388 + t361;
t37 = -pkin(4) * t304 - qJD(5) * t347 + t269;
t32 = -t190 * t300 + t306;
t27 = -mrSges(5,2) * t133 - mrSges(5,3) * t66;
t26 = mrSges(5,1) * t133 + mrSges(5,3) * t65;
t25 = -pkin(4) * t333 - t33;
t22 = -qJ(5) * t333 + (t190 * t326 - qJD(5)) * t191 - t306;
t20 = (qJ(6) * qJD(4) - qJD(5)) * t347 + (-qJD(3) * t295 + qJD(6) * t194) * t190 + t269;
t19 = mrSges(7,2) * t65 + mrSges(7,3) * t66;
t17 = -mrSges(6,2) * t66 + mrSges(6,3) * t65;
t16 = (-pkin(5) * t328 + qJ(5) * qJD(3)) * t194 + (qJD(5) + (pkin(5) * qJD(3) - t326) * t190) * t191 + t306;
t15 = -pkin(5) * t301 + qJD(3) * t219 - qJD(6) * t191 - t33;
t13 = -t65 * Ifges(5,4) - t66 * Ifges(5,2) + t133 * Ifges(5,6);
t6 = pkin(4) * t66 + t200;
t3 = qJD(6) * t135 + t388 * t66 + t200;
t9 = [(Ifges(4,1) * t476 + Ifges(4,4) * t477 + Ifges(4,5) * qJDD(3) - t237 * t415 + t246 * t414 - t3 * t251 - t6 * t252 - t410 * t487 + t441 * t416 - t483 * t478 - t79 * t486) * t194 + (Ifges(6,4) * t415 + Ifges(5,6) * t414 + t425 / 0.2e1 + Ifges(4,4) * t221 / 0.2e1 - Ifges(4,2) * t144 / 0.2e1 + t482 * t410 + t458 * t416 + t459 * t478 - Ifges(4,6) * qJDD(3) + t464) * t191 + (t256 + 0.2e1 * mrSges(3,3)) * t166 + qJ(2) * (-t144 * mrSges(4,1) - mrSges(4,2) * t221) - pkin(1) * t325 + m(7) * (t1 * t43 + t15 * t21 + t16 * t23 + t2 * t69 + t20 * t24 + t3 * t82) + m(6) * (t22 * t39 + t25 * t38 + t37 * t42 + t4 * t83 + t5 * t84 + t6 * t95) + (t116 * t486 - t479) * t334 + t223 * t276 + qJD(3) ^ 2 * t240 / 0.2e1 + t250 * t476 + t247 * t477 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t152 * t302 + (qJD(3) * t208 + t236 * t327) * t407 + (qJD(3) * t203 - t245 * t327) * t409 + m(3) * (-pkin(1) * t179 + t274) + (-mrSges(5,1) * t116 - mrSges(6,1) * t39 + mrSges(7,1) * t23 + mrSges(6,2) * t42 + mrSges(5,3) * t54 - mrSges(7,3) * t24) * (-t299 + t304) + ((-t238 + t242 + t243) * t327 + t426 * qJD(3)) * t404 + (t280 * t454 + t281 * t49) * t193 + (-t13 / 0.2e1 - t7 * mrSges(5,3) - t2 * mrSges(7,1) + t4 * mrSges(6,1) + t457 / 0.2e1) * t355 + m(5) * (t32 * t54 - t33 * t53 + t7 * t94 + t8 * t93) + (-t53 * mrSges(5,1) - t54 * mrSges(5,2) + t38 * mrSges(6,2) + t23 * mrSges(7,2) - t39 * mrSges(6,3) - t21 * mrSges(7,3) - t465) * t333 + (t280 * t48 + t281 * t453) * t190 + t101 * t351 + (t5 * mrSges(6,1) + t1 * mrSges(7,1) - t8 * mrSges(5,3) + t468) * t347 + t43 * t28 + t434 * mrSges(4,3) + m(4) * (-t196 * t434 + t274) + (-mrSges(6,1) * t38 - mrSges(7,1) * t21 - mrSges(5,2) * t116 + mrSges(7,2) * t24 - mrSges(5,3) * t53 + mrSges(6,3) * t42) * t471 + (qJD(3) * t442 + t327 * t440) * t408 + ((t231 - t248) * t327 + t443 * qJD(3)) * t406 + t69 * t30 + t20 * t74 + t37 * t77 + t82 * t19 + t83 * t29 + t84 * t31 + t32 * t87 + t33 * t88 + t15 * t89 + t22 * t90 + t16 * t91 + t25 * t92 + t93 * t26 + t94 * t27 + t95 * t17 + t446 * t168 + t449 * t345 + t79 * t126 + (-m(7) * (t174 + t279 + t444) - m(6) * (t211 + t444) - m(3) * t279 - m(4) * t185 - m(5) * t211 + t460 * t346 + (-m(4) * t196 + m(7) * pkin(7) - t430) * t192 - t420 * t113 - t431 * t112 + t422 * t195) * g(1) + (-m(4) * t305 - m(5) * (t268 - t318) - m(7) * t214 - m(3) * t340 - m(6) * (t214 - t318) + t460 * t348 + t430 * t195 - t420 * t111 - t431 * t110 + t422 * t192) * g(2) + qJD(2) * t139 + t179 * mrSges(3,2) + t432 * t322; t325 + m(3) * t179 + (qJ(2) * t412 - mrSges(3,3)) * t418 + (-t139 - m(5) * t228 + t461 - m(7) * (t190 * t23 - t193 * t21) + t201) * qJD(1) + (-t17 - t19 + (t190 * t319 + t193 * t320 + t152) * qJD(3) + m(4) * t85 + m(6) * (-t324 * t39 + t335 * t38 - t6) + m(7) * (t21 * t335 + t23 * t324 - t3) + m(5) * (t324 * t54 + t335 * t53 - t79) + t449) * t194 + (t101 + (t27 + t455) * t193 + (-t26 + t456) * t190 + (t387 + t446) * qJD(3) + t201 * qJD(4) + m(4) * t86 + m(6) * (qJD(3) * t42 + t473) + m(7) * (qJD(3) * t24 + t421) + m(5) * (qJD(3) * t116 + t474)) * t191 + (t396 - t397) * (-t412 - t433); (-t48 / 0.2e1 + t491) * t328 + t494 * t329 + ((t203 / 0.2e1 - t442 / 0.2e1) * t135 + t426 * t405 - t21 * (-mrSges(7,1) * t353 - t194 * mrSges(7,3)) - t38 * (-mrSges(6,1) * t353 + mrSges(6,2) * t194) + t53 * (mrSges(5,1) * t194 + mrSges(5,3) * t353) - t39 * (-mrSges(6,1) * t356 - mrSges(6,3) * t194) - t54 * (-mrSges(5,2) * t194 + mrSges(5,3) * t356) - t23 * (mrSges(7,1) * t356 + mrSges(7,2) * t194) + (t208 / 0.2e1 - t443 / 0.2e1) * t136) * qJD(1) + ((-t246 / 0.2e1 + t234 / 0.2e1 - t233 / 0.2e1) * t135 + t244 * t405 + t419 + t439 * t404) * qJD(4) + (t419 + t479) * t337 + (-m(7) * t400 + t256 + t433 * t186 + (-mrSges(7,1) - t460) * t194 + t427 * t191) * g(3) - Ifges(4,5) * t221 - pkin(3) * t18 + t3 * (-t193 * mrSges(7,3) - t368) + (-pkin(3) * t79 - t116 * t145 + t53 * t72 - t54 * t73) * m(5) + (t147 * t6 - t38 * t67 - t39 * t64 + t42 * t445) * m(6) + (t257 + t427 * t194 + (-t267 + t460) * t191) * t396 + (m(5) * (-qJD(4) * t228 + t258) + (-t29 + t27) * t193 + t450 * t328 + t451 * t329 + (-t26 + t31) * t190 + t259 * m(6) + (m(5) + m(6)) * t191 * t396 + t461 * qJD(4)) * pkin(8) + (t248 + t236) * t416 + t240 * t276 + t231 * t415 + (-t242 / 0.2e1 - t243 / 0.2e1) * t133 - t152 * t357 + t238 * t410 - t257 * t397 + t465 * t336 + (-m(5) * t341 + t463 * (t192 * pkin(4) * t347 + t348 * t360 + t341) + (-(m(7) * pkin(5) + mrSges(7,1)) * t191 + (-t193 * t296 - t466) * t194) * t192) * g(1) + t421 * mrSges(7,1) + t468 * t190 + t6 * t253 - t79 * t255 + (-t428 + t473) * mrSges(6,1) + (-t428 + t474) * mrSges(5,3) + (t245 + t440) * t414 + (t237 + t441) * t330 / 0.2e1 + t85 * mrSges(4,1) - t86 * mrSges(4,2) - t73 * t87 - t72 * t88 - t64 * t90 - t67 * t92 + t445 * t77 - t446 * t145 + t447 * t91 + t448 * t89 + t125 * t19 + t457 * t492 + t452 * t74 + (t1 * t156 + t125 * t3 + t157 * t2 + t21 * t448 + t23 * t447 + t24 * t452) * m(7) + Ifges(4,6) * t144 + t147 * t17 + t156 * t28 + t157 * t30 + t193 * t13 / 0.2e1 + Ifges(4,3) * qJDD(3) + (t223 / 0.2e1 - t432) * t418; (t135 * t481 + t136 * t480) * t405 + t226 * t91 - t388 * t28 - pkin(4) * t31 - t39 * t384 + (Ifges(6,2) * t135 + t369 + t49) * t406 + t464 + t425 + (t136 * t484 + t128 - t373 + t48) * t409 + (-t135 * t485 + t127 - t370 + t454) * t407 + (-pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t39 - t42 * t75) * m(6) + t23 * t382 + t21 * t383 + t38 * t385 + (qJ(5) * t2 - t1 * t388 + t436 * t21 + t23 * t469 - t24 * t40) * m(7) + t436 * t89 - t40 * t74 - t75 * t77 + (-m(6) * t38 + t380 - t450) * t54 - t386 * qJD(5) + (-m(6) * t39 + t381 - t451) * t53 + (-Ifges(5,2) * t136 - t129 + t453) * t408 + t455 * qJ(5) - t42 * (-mrSges(6,2) * t136 + mrSges(6,3) * t135) - t24 * (mrSges(7,2) * t135 + mrSges(7,3) * t136) - t116 * (mrSges(5,1) * t136 - mrSges(5,2) * t135) + (t126 + t463 * t342 + (t190 * t296 - t252 - t365) * t194) * g(3) + (t463 * (t112 * pkin(4) - qJ(5) * t113) + t431 * t113 - t420 * t112) * g(2) + (t463 * (-t110 * pkin(4) + qJ(5) * t111) - t431 * t111 + t420 * t110) * g(1); t387 * t136 + t386 * t169 + (t136 * t24 - t169 * t23 + t1 + t213) * m(7) + (t136 * t42 + t169 * t39 + t213 + t5) * m(6) + t456; -t135 * t74 + t169 * t89 + (-g(1) * t111 + g(2) * t113 - g(3) * t347 - t24 * t135 + t21 * t169 + t2) * m(7) + t30;];
tau  = t9;
