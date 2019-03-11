% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:34
% EndTime: 2019-03-08 21:13:09
% DurationCPUTime: 22.55s
% Computational Cost: add. (5543->724), mult. (13037->958), div. (0->0), fcn. (9785->12), ass. (0->349)
t486 = Ifges(6,1) + Ifges(5,1);
t259 = cos(qJ(3));
t413 = -t259 / 0.2e1;
t256 = sin(qJ(3));
t404 = Ifges(4,4) * t256;
t490 = -t404 / 0.2e1;
t359 = qJD(2) * qJD(3);
t210 = -t259 * qJDD(2) + t256 * t359;
t417 = t210 / 0.2e1;
t211 = qJDD(2) * t256 + t259 * t359;
t251 = sin(pkin(11));
t253 = cos(pkin(11));
t163 = qJDD(3) * t251 + t211 * t253;
t419 = t163 / 0.2e1;
t466 = Ifges(5,5) + Ifges(6,4);
t489 = t466 * t417 + t486 * t419;
t367 = qJD(2) * t259;
t240 = qJD(6) + t367;
t416 = -t240 / 0.2e1;
t369 = qJD(2) * t256;
t203 = -t253 * qJD(3) + t251 * t369;
t204 = qJD(3) * t251 + t253 * t369;
t255 = sin(qJ(6));
t258 = cos(qJ(6));
t292 = t203 * t255 + t204 * t258;
t423 = -t292 / 0.2e1;
t106 = t203 * t258 - t204 * t255;
t425 = -t106 / 0.2e1;
t488 = Ifges(7,5) * t423 + Ifges(7,6) * t425 + Ifges(7,3) * t416;
t162 = -t253 * qJDD(3) + t211 * t251;
t420 = t162 / 0.2e1;
t487 = Ifges(4,2) * t413 + t490;
t335 = -t367 / 0.2e1;
t485 = -Ifges(5,4) + Ifges(6,5);
t465 = Ifges(5,6) - Ifges(6,6);
t464 = Ifges(5,3) + Ifges(6,2);
t260 = cos(qJ(2));
t252 = sin(pkin(6));
t372 = qJD(1) * t252;
t345 = t260 * t372;
t314 = t259 * t345;
t257 = sin(qJ(2));
t346 = t257 * t372;
t140 = t251 * t314 - t253 * t346;
t347 = -pkin(8) * t251 - pkin(4);
t377 = t253 * t259;
t356 = pkin(9) * t377;
t295 = pkin(3) * t256 - qJ(4) * t259;
t362 = qJD(4) * t256;
t185 = qJD(3) * t295 - t362;
t385 = t185 * t253;
t484 = -t140 - t385 + (-t356 + (-pkin(5) + t347) * t256) * qJD(3);
t376 = t259 * t260;
t153 = (t251 * t257 + t253 * t376) * t252;
t141 = qJD(1) * t153;
t169 = t251 * t185;
t378 = t253 * t256;
t381 = t251 * t259;
t360 = qJD(5) * t259;
t366 = qJD(3) * t256;
t448 = qJ(5) * t366 - t360;
t483 = t141 - t169 - (-pkin(8) * t378 + pkin(9) * t381) * qJD(3) - t448;
t399 = Ifges(6,5) * t251;
t402 = Ifges(5,4) * t251;
t482 = t253 * t486 + t399 - t402;
t387 = qJ(5) * t251;
t410 = pkin(4) * t253;
t481 = -t387 - t410;
t480 = m(5) * qJD(4);
t479 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t244 = Ifges(4,4) * t367;
t478 = Ifges(4,1) * t369 + Ifges(4,5) * qJD(3) + t244 + t253 * (t203 * t485 + t204 * t486 - t367 * t466);
t212 = qJD(2) * pkin(8) + t346;
t200 = t256 * t212;
t370 = qJD(2) * t252;
t338 = qJD(1) * t370;
t230 = t260 * t338;
t358 = qJDD(1) * t252;
t181 = t257 * t358 + t230;
t165 = qJDD(2) * pkin(8) + t181;
t254 = cos(pkin(6));
t371 = qJD(1) * t254;
t235 = t259 * t371;
t357 = qJDD(1) * t254;
t348 = qJD(3) * t235 + t259 * t165 + t256 * t357;
t41 = qJDD(3) * qJ(4) + (qJD(4) - t200) * qJD(3) + t348;
t229 = t257 * t338;
t180 = t260 * t358 - t229;
t164 = -qJDD(2) * pkin(2) - t180;
t59 = pkin(3) * t210 - qJ(4) * t211 - qJD(2) * t362 + t164;
t19 = -t251 * t41 + t253 * t59;
t310 = qJDD(5) - t19;
t16 = -pkin(4) * t210 + t310;
t477 = t16 * mrSges(6,2) + t420 * t485 + t489;
t207 = t251 * t258 - t253 * t255;
t291 = t251 * t255 + t253 * t258;
t476 = mrSges(7,1) * t291 + mrSges(7,2) * t207;
t344 = t256 * t371;
t146 = t212 * t259 + t344;
t138 = qJD(3) * qJ(4) + t146;
t411 = pkin(3) * t259;
t217 = -qJ(4) * t256 - pkin(2) - t411;
t150 = qJD(2) * t217 - t345;
t55 = -t251 * t138 + t150 * t253;
t42 = pkin(4) * t367 + qJD(5) - t55;
t28 = pkin(5) * t367 - pkin(9) * t204 + t42;
t56 = t253 * t138 + t251 * t150;
t43 = -qJ(5) * t367 + t56;
t32 = pkin(9) * t203 + t43;
t10 = t255 * t28 + t258 * t32;
t9 = -t255 * t32 + t258 * t28;
t475 = -t9 * mrSges(7,1) + t10 * mrSges(7,2) - Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t487 + t464 * t335 + t466 * t204 / 0.2e1 - t465 * t203 / 0.2e1 + t488;
t474 = m(7) * pkin(5);
t26 = qJD(6) * t106 + t162 * t255 + t163 * t258;
t430 = t26 / 0.2e1;
t27 = -qJD(6) * t292 + t162 * t258 - t163 * t255;
t429 = t27 / 0.2e1;
t473 = -m(7) - m(6);
t205 = qJDD(6) - t210;
t418 = t205 / 0.2e1;
t412 = t259 / 0.2e1;
t155 = pkin(8) * t377 + t251 * t217;
t137 = -qJ(5) * t259 + t155;
t382 = t251 * t256;
t111 = pkin(9) * t382 + t137;
t236 = pkin(8) * t381;
t249 = t259 * pkin(4);
t90 = pkin(5) * t259 + t236 + t249 + (-pkin(9) * t256 - t217) * t253;
t33 = -t111 * t255 + t258 * t90;
t472 = qJD(6) * t33 + t255 * t484 - t258 * t483;
t34 = t111 * t258 + t255 * t90;
t471 = -qJD(6) * t34 + t255 * t483 + t258 * t484;
t469 = qJD(3) / 0.2e1;
t468 = mrSges(3,2) - mrSges(4,3);
t406 = -pkin(9) + qJ(4);
t218 = t406 * t251;
t219 = t406 * t253;
t130 = t218 * t255 + t219 * t258;
t426 = pkin(4) + pkin(5);
t145 = -t200 + t235;
t209 = t295 * qJD(2);
t79 = -t251 * t145 + t209 * t253;
t44 = (-t256 * t426 - t356) * qJD(2) - t79;
t343 = t251 * t367;
t80 = t253 * t145 + t251 * t209;
t65 = qJ(5) * t369 + t80;
t53 = pkin(9) * t343 + t65;
t463 = qJD(4) * t207 - qJD(6) * t130 + t255 * t53 - t258 * t44;
t129 = t218 * t258 - t219 * t255;
t462 = qJD(4) * t291 + qJD(6) * t129 - t255 * t44 - t258 * t53;
t91 = -mrSges(6,2) * t162 + mrSges(6,3) * t210;
t92 = -mrSges(5,2) * t210 - mrSges(5,3) * t162;
t460 = t91 + t92;
t93 = mrSges(5,1) * t210 - mrSges(5,3) * t163;
t94 = -t210 * mrSges(6,1) + t163 * mrSges(6,2);
t459 = t94 - t93;
t113 = mrSges(6,1) * t203 - mrSges(6,3) * t204;
t36 = -mrSges(7,1) * t106 + mrSges(7,2) * t292;
t392 = t113 - t36;
t308 = mrSges(4,1) * t259 - mrSges(4,2) * t256;
t456 = t308 + mrSges(3,1);
t386 = qJ(5) * t253;
t286 = -t251 * t426 + t386;
t361 = qJD(5) * t251;
t455 = t361 + t344 - (qJD(2) * t286 - t212) * t259;
t390 = sin(pkin(10));
t323 = t390 * t260;
t391 = cos(pkin(10));
t326 = t391 * t257;
t187 = t254 * t326 + t323;
t328 = t252 * t391;
t123 = t187 * t256 + t259 * t328;
t454 = t481 * t123;
t324 = t390 * t257;
t325 = t391 * t260;
t189 = -t254 * t324 + t325;
t327 = t252 * t390;
t125 = t189 * t256 - t259 * t327;
t453 = t481 * t125;
t156 = -mrSges(6,2) * t203 - mrSges(6,3) * t367;
t157 = mrSges(5,2) * t367 - mrSges(5,3) * t203;
t375 = t156 + t157;
t158 = -mrSges(5,1) * t367 - mrSges(5,3) * t204;
t159 = mrSges(6,1) * t367 + mrSges(6,2) * t204;
t374 = t159 - t158;
t282 = t207 * t259;
t160 = qJD(2) * t282;
t191 = t207 * qJD(6);
t452 = t160 + t191;
t281 = t291 * t259;
t161 = qJD(2) * t281;
t190 = t291 * qJD(6);
t451 = t161 + t190;
t380 = t252 * t257;
t192 = -t254 * t259 + t256 * t380;
t450 = t481 * t192;
t352 = mrSges(4,3) * t369;
t449 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t203 + mrSges(5,2) * t204 + t352;
t51 = -t212 * t366 + t348;
t52 = -qJD(3) * t146 - t165 * t256 + t259 * t357;
t447 = -t256 * t52 + t259 * t51;
t20 = t251 * t59 + t253 * t41;
t446 = -t19 * t251 + t20 * t253;
t69 = t162 * mrSges(6,1) - t163 * mrSges(6,3);
t70 = t162 * mrSges(5,1) + t163 * mrSges(5,2);
t8 = -t27 * mrSges(7,1) + t26 * mrSges(7,2);
t445 = t69 + t70 - t8;
t398 = Ifges(6,5) * t253;
t296 = Ifges(6,3) * t251 + t398;
t444 = t203 * (Ifges(6,6) * t256 + t259 * t296) + (t256 * t466 + t259 * t482) * t204;
t442 = -m(7) * t406 + mrSges(4,2) + t479;
t305 = -mrSges(6,1) * t253 - mrSges(6,3) * t251;
t307 = -mrSges(5,1) * t253 + mrSges(5,2) * t251;
t441 = t253 * t474 + mrSges(4,1) - t305 - t307 + t476;
t15 = t210 * qJ(5) - qJD(2) * t360 + t20;
t11 = pkin(9) * t162 + t15;
t7 = -pkin(9) * t163 - t210 * t426 + t310;
t1 = qJD(6) * t9 + t11 * t258 + t255 * t7;
t2 = -qJD(6) * t10 - t11 * t255 + t258 * t7;
t440 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t287 = qJD(3) * pkin(3) - qJD(4) + t145;
t439 = m(5) * t287 - t449;
t437 = m(7) * pkin(9) + t479;
t436 = Ifges(6,5) * t419 + Ifges(6,6) * t417 - t163 * Ifges(5,4) / 0.2e1 - t210 * Ifges(5,6) / 0.2e1 - t15 * mrSges(6,2) + (Ifges(6,3) + Ifges(5,2)) * t420;
t435 = -mrSges(7,1) * t255 - mrSges(7,2) * t258 + mrSges(5,2) - mrSges(6,3);
t213 = -qJD(2) * pkin(2) - t345;
t401 = Ifges(5,4) * t253;
t300 = -Ifges(5,2) * t251 + t401;
t304 = mrSges(6,1) * t251 - mrSges(6,3) * t253;
t306 = mrSges(5,1) * t251 + mrSges(5,2) * t253;
t271 = qJ(5) * t204 + t287;
t50 = pkin(4) * t203 - t271;
t434 = -t213 * (mrSges(4,1) * t256 + mrSges(4,2) * t259) + t287 * t259 * t306 - t259 * t50 * t304 + t203 * (Ifges(5,6) * t256 + t259 * t300) / 0.2e1 - t56 * (-mrSges(5,2) * t256 - mrSges(5,3) * t381) - t55 * (mrSges(5,1) * t256 - mrSges(5,3) * t377) - t43 * (-mrSges(6,2) * t381 + mrSges(6,3) * t256) - t42 * (-mrSges(6,1) * t256 + mrSges(6,2) * t377);
t433 = -mrSges(7,1) * t258 + mrSges(7,2) * t255 - mrSges(5,1) - mrSges(6,1) - t474;
t261 = qJD(2) ^ 2;
t432 = Ifges(7,4) * t430 + Ifges(7,2) * t429 + Ifges(7,6) * t418;
t431 = Ifges(7,1) * t430 + Ifges(7,4) * t429 + Ifges(7,5) * t418;
t400 = Ifges(7,4) * t292;
t30 = Ifges(7,2) * t106 + Ifges(7,6) * t240 + t400;
t428 = t30 / 0.2e1;
t101 = Ifges(7,4) * t106;
t31 = Ifges(7,1) * t292 + Ifges(7,5) * t240 + t101;
t427 = t31 / 0.2e1;
t424 = t106 / 0.2e1;
t422 = t292 / 0.2e1;
t421 = -t162 / 0.2e1;
t415 = t240 / 0.2e1;
t403 = Ifges(4,4) * t259;
t45 = -qJDD(3) * pkin(3) + qJDD(4) - t52;
t395 = t256 * t45;
t186 = -t254 * t325 + t324;
t384 = t186 * t256;
t188 = t254 * t323 + t326;
t383 = t188 * t256;
t379 = t252 * t260;
t373 = pkin(2) * t379 + pkin(8) * t380;
t368 = qJD(2) * t257;
t365 = qJD(3) * t259;
t364 = qJD(4) * t251;
t363 = qJD(4) * t253;
t355 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t205;
t354 = pkin(8) * t366;
t351 = mrSges(4,3) * t367;
t350 = t256 * t379;
t349 = t252 * t376;
t342 = t252 * t368;
t341 = t260 * t370;
t340 = qJD(5) * t378;
t331 = pkin(3) + t387;
t330 = -t186 * pkin(2) + pkin(8) * t187;
t329 = -t188 * pkin(2) + pkin(8) * t189;
t322 = -t359 / 0.2e1;
t321 = t359 / 0.2e1;
t115 = t123 * pkin(3);
t124 = t187 * t259 - t256 * t328;
t320 = qJ(4) * t124 - t115;
t116 = t125 * pkin(3);
t126 = t189 * t259 + t256 * t327;
t319 = qJ(4) * t126 - t116;
t184 = t192 * pkin(3);
t193 = t254 * t256 + t259 * t380;
t318 = qJ(4) * t193 - t184;
t154 = t217 * t253 - t236;
t316 = pkin(3) * t349 + qJ(4) * t350 + t373;
t299 = Ifges(6,4) * t253 + Ifges(6,6) * t251;
t298 = Ifges(4,5) * t259 - Ifges(4,6) * t256;
t297 = Ifges(5,5) * t253 - Ifges(5,6) * t251;
t294 = pkin(4) * t251 - t386;
t73 = -mrSges(7,2) * t240 + mrSges(7,3) * t106;
t74 = mrSges(7,1) * t240 - mrSges(7,3) * t292;
t293 = -t255 * t74 + t258 * t73;
t121 = t193 * t251 + t253 * t379;
t122 = t193 * t253 - t251 * t379;
t39 = t121 * t258 - t122 * t255;
t40 = t121 * t255 + t122 * t258;
t290 = -qJ(4) * t384 - t186 * t411 + t330;
t289 = -qJ(4) * t383 - t188 * t411 + t329;
t120 = -t253 * t354 + t169;
t288 = pkin(8) + t294;
t284 = t256 * (Ifges(4,1) * t259 - t404);
t60 = t124 * t251 - t186 * t253;
t62 = t126 * t251 - t188 * t253;
t283 = -g(1) * t62 - g(2) * t60 - g(3) * t121;
t176 = t207 * t256;
t280 = -pkin(8) + t286;
t279 = -g(1) * t125 - g(2) * t123 - g(3) * t192;
t264 = t259 * (Ifges(6,2) * t256 + t259 * t299);
t263 = t259 * (Ifges(5,3) * t256 + t259 * t297);
t21 = t162 * pkin(4) - t163 * qJ(5) - t204 * qJD(5) + t45;
t223 = -qJD(3) * mrSges(4,2) + t351;
t214 = -t331 - t410;
t208 = t308 * qJD(2);
t195 = t253 * t426 + t331;
t183 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t211;
t182 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t210;
t177 = t291 * t256;
t175 = t288 * t256;
t152 = t251 * t349 - t253 * t380;
t139 = -t154 + t249;
t133 = t280 * t256;
t131 = mrSges(4,1) * t210 + mrSges(4,2) * t211;
t128 = -qJD(3) * t192 + t259 * t341;
t127 = qJD(3) * t193 + t256 * t341;
t119 = t251 * t354 + t385;
t112 = t288 * t365 - t340;
t102 = t347 * t366 - t385;
t98 = t204 * Ifges(5,4) - t203 * Ifges(5,2) - Ifges(5,6) * t367;
t95 = t204 * Ifges(6,5) - Ifges(6,6) * t367 + t203 * Ifges(6,3);
t89 = t280 * t365 + t340;
t88 = t344 + (qJD(2) * t294 + t212) * t259;
t87 = t120 + t448;
t86 = qJD(3) * t282 - t190 * t256;
t85 = qJD(3) * t281 + qJD(6) * t176;
t84 = -t188 * t377 + t189 * t251;
t83 = -t188 * t381 - t189 * t253;
t82 = -t186 * t377 + t187 * t251;
t81 = -t186 * t381 - t187 * t253;
t78 = t128 * t253 + t251 * t342;
t77 = t128 * t251 - t253 * t342;
t68 = -pkin(4) * t369 - t79;
t63 = t126 * t253 + t188 * t251;
t61 = t124 * t253 + t186 * t251;
t35 = -t203 * t426 + t271;
t23 = -mrSges(7,2) * t205 + mrSges(7,3) * t27;
t22 = mrSges(7,1) * t205 - mrSges(7,3) * t26;
t14 = pkin(5) * t162 + t21;
t13 = qJD(6) * t39 + t255 * t77 + t258 * t78;
t12 = -qJD(6) * t40 - t255 * t78 + t258 * t77;
t3 = [m(2) * qJDD(1) + t12 * t74 + t128 * t223 + t13 * t73 + t193 * t182 + t39 * t22 + t40 * t23 + t375 * t78 + t374 * t77 + t460 * t122 + t459 * t121 + (-t183 + t445) * t192 + (t392 + t449) * t127 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t261 - t131) * t260 + (-mrSges(3,1) * t261 - mrSges(3,2) * qJDD(2) - qJD(2) * t208) * t257) * t252 + (-m(2) - m(3) - m(4) - m(5) + t473) * g(3) + m(4) * (-t127 * t145 + t128 * t146 - t192 * t52 + t193 * t51 + (-t164 * t260 + t213 * t368) * t252) + m(3) * (qJDD(1) * t254 ^ 2 + (t180 * t260 + t181 * t257) * t252) + m(6) * (t121 * t16 + t122 * t15 + t127 * t50 + t192 * t21 + t42 * t77 + t43 * t78) + m(5) * (-t121 * t19 + t122 * t20 - t127 * t287 + t192 * t45 - t55 * t77 + t56 * t78) + m(7) * (t1 * t40 + t10 * t13 + t12 * t9 - t127 * t35 + t14 * t192 + t2 * t39); -t164 * t308 + (-m(6) * t50 + m(7) * t35 - t392 + t439) * t256 * t345 + (t449 * t365 + (-t183 + t70) * t256 + (-t287 * t365 + t395) * m(5) + t182 * t259 + (t447 + (-t145 * t259 - t146 * t256) * qJD(3)) * m(4)) * pkin(8) + (-t146 * mrSges(4,3) - Ifges(7,5) * t422 - Ifges(7,6) * t424 - Ifges(7,3) * t415 + t475) * t366 + (-mrSges(5,3) * t19 + t477) * t378 + (t251 * t95 + t478) * t365 / 0.2e1 + (t112 * t50 + t137 * t15 + t139 * t16 + t175 * t21 + (-t141 + t87) * t43 + (-t140 + t102) * t42) * m(6) + (t230 - t181) * mrSges(3,2) + (-t354 - t314) * t223 + (t1 * t176 + t10 * t86 - t177 * t2 - t85 * t9) * mrSges(7,3) + (Ifges(7,4) * t177 + Ifges(7,2) * t176) * t429 + (Ifges(7,4) * t85 + Ifges(7,2) * t86) * t424 + (t229 + t180) * mrSges(3,1) + (Ifges(7,5) * t85 + Ifges(7,6) * t86) * t415 + (Ifges(7,5) * t177 + Ifges(7,6) * t176) * t418 + (Ifges(4,4) * t412 + t256 * Ifges(4,1) + t403 / 0.2e1) * t211 + t85 * t427 + t86 * t428 + t177 * t431 + t176 * t432 + t355 * t412 + t306 * t395 + (t264 + t263) * t322 + t208 * t346 + t284 * t321 + (t21 * t304 + t296 * t420 + t300 * t421 + t482 * t419 + (t297 + t299) * t417) * t256 + (-t145 * t365 + t447) * mrSges(4,3) + (-t20 * mrSges(5,3) + t436) * t382 - t251 * t98 * t365 / 0.2e1 + (t154 * t19 + t155 * t20 + (-t141 + t120) * t56 + (t140 + t119) * t55) * m(5) + (-pkin(2) * t164 - (t213 * t257 + (-t145 * t256 + t146 * t259) * t260) * t372) * m(4) + (t20 * mrSges(5,2) - t15 * mrSges(6,3) - t19 * mrSges(5,1) + t16 * mrSges(6,1) + Ifges(7,6) * t429 + Ifges(7,5) * t430 + Ifges(7,3) * t418 - Ifges(6,6) * t420 - Ifges(5,6) * t421 + (-Ifges(4,2) * t256 + t403) * t321 - t466 * t419 - t464 * t417 + t440) * t259 + t444 * t469 + Ifges(3,3) * qJDD(2) + (Ifges(7,1) * t177 + Ifges(7,4) * t176) * t430 + (Ifges(7,1) * t85 + Ifges(7,4) * t86) * t422 + (t298 * t469 - t434) * qJD(3) + t175 * t69 - t14 * (-mrSges(7,1) * t176 + mrSges(7,2) * t177) + t102 * t159 + t154 * t93 + t155 * t92 + t87 * t156 + t120 * t157 + t119 * t158 - pkin(2) * t131 + t133 * t8 + t137 * t91 + t139 * t94 - t374 * t140 - t375 * t141 + (-Ifges(4,2) * t412 + t464 * t413 + t487 + t490) * t210 + (-t162 * t465 + t163 * t466) * t413 + t471 * t74 + t472 * t73 + (t1 * t34 + t10 * t472 - t133 * t14 + t2 * t33 + t35 * t89 + t471 * t9) * m(7) + (Ifges(4,5) * t256 + 0.2e1 * Ifges(4,6) * t412) * qJDD(3) + (-m(4) * t329 - m(5) * t289 + t473 * (t84 * pkin(4) + qJ(5) * t83 + t289) + t433 * t84 + t435 * t83 + t468 * t189 + t456 * t188 - t437 * t383) * g(1) + (-m(4) * t330 - m(5) * t290 + t473 * (t82 * pkin(4) + qJ(5) * t81 + t290) + t433 * t82 + t435 * t81 + t468 * t187 + t456 * t186 - t437 * t384) * g(2) + (-m(4) * t373 - m(5) * t316 + t473 * (t153 * pkin(4) + qJ(5) * t152 + t316) + (t257 * t468 - t260 * t456) * t252 + t433 * t153 + t435 * t152 + t437 * t350) * g(3) + t33 * t22 + t34 * t23 + t35 * (-mrSges(7,1) * t86 + mrSges(7,2) * t85) + t89 * t36 + t112 * t113; t21 * t305 + t45 * t307 + t98 * t343 / 0.2e1 + (-pkin(3) * t45 + t446 * qJ(4) - t55 * t79 - t56 * t80) * m(5) + (Ifges(7,4) * t161 + Ifges(7,2) * t160) * t425 + (Ifges(7,4) * t207 - Ifges(7,2) * t291) * t429 + (Ifges(7,1) * t207 - Ifges(7,4) * t291) * t430 + (Ifges(7,5) * t207 - Ifges(7,6) * t291) * t418 - t291 * t432 + (-t1 * t291 - t10 * t452 - t2 * t207 + t451 * t9) * mrSges(7,3) + (-Ifges(4,2) * t335 - t475 - t488) * t369 + (-t398 + t401) * t419 + (t352 + t439) * t146 + (t244 + t478) * t335 + (-t68 + t364) * t159 + (-t80 + t363) * t157 + (-t361 - t88) * t113 - t14 * t476 + (t460 * qJ(4) + t56 * t480 - Ifges(6,3) * t420 + Ifges(5,2) * t421 + m(6) * (qJ(4) * t15 + qJD(4) * t43) + t465 * t417 - t436) * t253 + (t21 * t214 - t42 * t68 - t43 * t65 - t50 * t88) * m(6) + (-Ifges(7,5) * t190 - Ifges(7,6) * t191) * t415 + (-Ifges(7,1) * t190 - Ifges(7,4) * t191) * t422 + (-Ifges(7,4) * t190 - Ifges(7,2) * t191) * t424 + (Ifges(7,5) * t161 + Ifges(7,6) * t160) * t416 + (-t364 - t79) * t158 - t190 * t427 - t191 * t428 + t207 * t431 + (-t223 + t351) * t145 + (-t65 + t363) * t156 + t399 * t420 + t402 * t421 + (t95 * t335 + t459 * qJ(4) - t55 * t480 + m(6) * (qJ(4) * t16 + qJD(4) * t42 - qJD(5) * t50) + t477 + t489) * t251 + (t434 - t444 / 0.2e1) * qJD(2) + t298 * t322 + (-t284 / 0.2e1 + t264 / 0.2e1 + t263 / 0.2e1) * t261 + (Ifges(7,1) * t161 + Ifges(7,4) * t160) * t423 + Ifges(4,3) * qJDD(3) - Ifges(4,6) * t210 + Ifges(4,5) * t211 + t214 * t69 + t195 * t8 - t160 * t30 / 0.2e1 - t161 * t31 / 0.2e1 + t129 * t22 + t130 * t23 + t446 * mrSges(5,3) + (-m(7) * (-t184 + t450) - m(6) * (t318 + t450) - m(5) * t318 + t442 * t193 + t441 * t192) * g(3) + (mrSges(7,1) * t452 - mrSges(7,2) * t451) * t35 + (-m(7) * (-t116 + t453) - m(6) * (t319 + t453) - m(5) * t319 + t442 * t126 + t441 * t125) * g(1) + (-m(7) * (-t115 + t454) - m(6) * (t320 + t454) - m(5) * t320 + t442 * t124 + t441 * t123) * g(2) + t455 * t36 + t462 * t73 + t463 * t74 + (t1 * t130 + t462 * t10 + t129 * t2 - t14 * t195 + t455 * t35 + t463 * t9) * m(7) - t51 * mrSges(4,2) + t52 * mrSges(4,1) - pkin(3) * t70; -t292 * t74 + t106 * t73 + t375 * t203 - t374 * t204 + (t10 * t106 - t292 * t9 + t14 + t279) * m(7) + (t203 * t43 - t204 * t42 + t21 + t279) * m(6) + (t203 * t56 + t204 * t55 + t279 + t45) * m(5) + t445; t258 * t22 + t255 * t23 + t392 * t204 + t293 * qJD(6) + (t156 + t293) * t367 + t94 + (t1 * t255 + t2 * t258 - t204 * t35 + t283 + t240 * (t10 * t258 - t255 * t9)) * m(7) + (t204 * t50 + t367 * t43 + t16 + t283) * m(6); -t35 * (mrSges(7,1) * t292 + mrSges(7,2) * t106) + (Ifges(7,1) * t106 - t400) * t423 + t30 * t422 + (Ifges(7,5) * t106 - Ifges(7,6) * t292) * t416 - t9 * t73 + t10 * t74 - g(1) * ((-t255 * t63 + t258 * t62) * mrSges(7,1) + (-t255 * t62 - t258 * t63) * mrSges(7,2)) - g(2) * ((-t255 * t61 + t258 * t60) * mrSges(7,1) + (-t255 * t60 - t258 * t61) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t39 - mrSges(7,2) * t40) + (t10 * t292 + t106 * t9) * mrSges(7,3) + t355 + (-Ifges(7,2) * t292 + t101 + t31) * t425 + t440;];
tau  = t3;
