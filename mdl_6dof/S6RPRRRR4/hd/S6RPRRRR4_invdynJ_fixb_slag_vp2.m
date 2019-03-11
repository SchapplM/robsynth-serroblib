% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR4
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:48
% EndTime: 2019-03-09 07:05:25
% DurationCPUTime: 21.62s
% Computational Cost: add. (26725->722), mult. (68001->952), div. (0->0), fcn. (54945->18), ass. (0->332)
t367 = qJD(1) * qJD(2);
t248 = qJ(2) * qJDD(1) + t367;
t278 = sin(pkin(11));
t279 = cos(pkin(11));
t376 = t278 ^ 2 + t279 ^ 2;
t276 = pkin(11) + qJ(3);
t270 = qJ(4) + t276;
t260 = qJ(5) + t270;
t251 = sin(t260);
t252 = cos(t260);
t281 = sin(qJ(6));
t420 = mrSges(7,2) * t281;
t501 = t251 * t420 + t252 * (m(7) * pkin(10) + mrSges(7,3));
t284 = sin(qJ(3));
t289 = cos(qJ(3));
t384 = t279 * t289;
t228 = -t278 * t284 + t384;
t217 = t228 * qJD(1);
t229 = t278 * t289 + t279 * t284;
t218 = t229 * qJD(1);
t283 = sin(qJ(4));
t288 = cos(qJ(4));
t181 = t217 * t283 + t218 * t288;
t282 = sin(qJ(5));
t287 = cos(qJ(5));
t335 = t288 * t217 - t218 * t283;
t488 = -t181 * t282 + t287 * t335;
t130 = -t488 + qJD(6);
t487 = t287 * t181 + t282 * t335;
t88 = pkin(5) * t487 - pkin(10) * t488;
t219 = t228 * qJD(3);
t186 = qJD(1) * t219 + qJDD(1) * t229;
t220 = t229 * qJD(3);
t187 = -qJD(1) * t220 + qJDD(1) * t228;
t109 = qJD(4) * t335 + t186 * t288 + t187 * t283;
t110 = -qJD(4) * t181 - t186 * t283 + t187 * t288;
t286 = cos(qJ(6));
t277 = qJD(3) + qJD(4);
t269 = qJD(5) + t277;
t424 = pkin(7) + qJ(2);
t240 = t424 * t278;
t230 = qJD(1) * t240;
t241 = t424 * t279;
t231 = qJD(1) * t241;
t191 = -t230 * t284 + t231 * t289;
t160 = pkin(8) * t217 + t191;
t157 = t288 * t160;
t387 = t231 * t284;
t190 = -t289 * t230 - t387;
t159 = -pkin(8) * t218 + t190;
t158 = qJD(3) * pkin(3) + t159;
t115 = t158 * t283 + t157;
t481 = pkin(9) * t335;
t95 = t115 + t481;
t393 = t287 * t95;
t155 = t283 * t160;
t114 = t288 * t158 - t155;
t497 = pkin(9) * t181;
t94 = t114 - t497;
t91 = pkin(4) * t277 + t94;
t60 = t282 * t91 + t393;
t58 = pkin(10) * t269 + t60;
t259 = pkin(2) * t279 + pkin(1);
t239 = -qJD(1) * t259 + qJD(2);
t192 = -pkin(3) * t217 + t239;
t145 = -pkin(4) * t335 + t192;
t74 = -pkin(5) * t488 - pkin(10) * t487 + t145;
t22 = -t281 * t58 + t286 * t74;
t23 = t281 * t74 + t286 * t58;
t272 = qJDD(3) + qJDD(4);
t338 = pkin(7) * qJDD(1) + t248;
t208 = t338 * t278;
t209 = t338 * t279;
t139 = -qJD(3) * t191 - t289 * t208 - t209 * t284;
t108 = qJDD(3) * pkin(3) - pkin(8) * t186 + t139;
t375 = qJD(3) * t289;
t138 = -qJD(3) * t387 - t284 * t208 + t289 * t209 - t230 * t375;
t111 = pkin(8) * t187 + t138;
t45 = -qJD(4) * t115 + t288 * t108 - t111 * t283;
t32 = pkin(4) * t272 - pkin(9) * t109 + t45;
t373 = qJD(4) * t288;
t374 = qJD(4) * t283;
t44 = t283 * t108 + t288 * t111 + t158 * t373 - t160 * t374;
t36 = pkin(9) * t110 + t44;
t371 = qJD(5) * t287;
t372 = qJD(5) * t282;
t10 = t282 * t32 + t287 * t36 + t91 * t371 - t372 * t95;
t11 = -qJD(5) * t60 - t282 * t36 + t287 * t32;
t120 = t269 * t286 - t281 * t487;
t266 = qJDD(5) + t272;
t55 = qJD(5) * t488 + t109 * t287 + t110 * t282;
t41 = qJD(6) * t120 + t266 * t281 + t286 * t55;
t121 = t269 * t281 + t286 * t487;
t42 = -qJD(6) * t121 + t266 * t286 - t281 * t55;
t56 = -qJD(5) * t487 - t109 * t282 + t110 * t287;
t54 = qJDD(6) - t56;
t14 = t41 * Ifges(7,4) + t42 * Ifges(7,2) + t54 * Ifges(7,6);
t15 = t41 * Ifges(7,1) + t42 * Ifges(7,4) + t54 * Ifges(7,5);
t236 = -qJDD(1) * t259 + qJDD(2);
t161 = -pkin(3) * t187 + t236;
t89 = -pkin(4) * t110 + t161;
t19 = -pkin(5) * t56 - pkin(10) * t55 + t89;
t7 = pkin(10) * t266 + t10;
t2 = qJD(6) * t22 + t19 * t281 + t286 * t7;
t321 = Ifges(7,5) * t286 - Ifges(7,6) * t281;
t306 = t130 * t321;
t414 = Ifges(7,4) * t281;
t323 = Ifges(7,1) * t286 - t414;
t307 = t121 * t323;
t413 = Ifges(7,4) * t286;
t322 = -Ifges(7,2) * t281 + t413;
t308 = t120 * t322;
t324 = mrSges(7,1) * t281 + mrSges(7,2) * t286;
t397 = t282 * t95;
t59 = t287 * t91 - t397;
t57 = -pkin(5) * t269 - t59;
t309 = t57 * t324;
t370 = qJD(6) * t281;
t344 = -t370 / 0.2e1;
t119 = Ifges(7,4) * t120;
t71 = t121 * Ifges(7,1) + t130 * Ifges(7,5) + t119;
t395 = t286 * t71;
t353 = t395 / 0.2e1;
t418 = mrSges(7,3) * t286;
t442 = t281 / 0.2e1;
t460 = t54 / 0.2e1;
t461 = t42 / 0.2e1;
t462 = t41 / 0.2e1;
t423 = mrSges(7,1) * t286;
t492 = t420 - t423;
t405 = t121 * Ifges(7,4);
t70 = t120 * Ifges(7,2) + t130 * Ifges(7,6) + t405;
t8 = -pkin(5) * t266 - t11;
t293 = -t10 * mrSges(6,2) + t2 * t418 + t15 * t442 + t286 * t14 / 0.2e1 + Ifges(6,3) * t266 + (Ifges(7,1) * t281 + t413) * t462 + (Ifges(7,2) * t286 + t414) * t461 + (Ifges(7,5) * t281 + Ifges(7,6) * t286) * t460 + Ifges(6,6) * t56 + Ifges(6,5) * t55 + t8 * t492 + t70 * t344 + t11 * mrSges(6,1) + (t309 + t353) * qJD(6) + (t308 + t307 + t306) * qJD(6) / 0.2e1;
t3 = -qJD(6) * t23 + t19 * t286 - t281 * t7;
t419 = mrSges(7,3) * t281;
t444 = -t269 / 0.2e1;
t450 = -t487 / 0.2e1;
t451 = -t488 / 0.2e1;
t452 = -t130 / 0.2e1;
t454 = -t121 / 0.2e1;
t455 = -t120 / 0.2e1;
t426 = t60 * mrSges(6,3);
t495 = t23 * mrSges(7,2);
t496 = t22 * mrSges(7,1);
t408 = Ifges(7,3) * t130;
t409 = Ifges(7,6) * t120;
t411 = Ifges(7,5) * t121;
t69 = t408 + t409 + t411;
t410 = Ifges(6,6) * t269;
t415 = Ifges(6,4) * t487;
t85 = Ifges(6,2) * t488 + t410 + t415;
t466 = -t145 * mrSges(6,1) - t496 + t495 + t426 + t85 / 0.2e1 - t69 / 0.2e1;
t427 = t59 * mrSges(6,3);
t129 = Ifges(6,4) * t488;
t412 = Ifges(6,5) * t269;
t86 = Ifges(6,1) * t487 + t129 + t412;
t467 = -t145 * mrSges(6,2) - t309 - t395 / 0.2e1 + t70 * t442 + t427 - t86 / 0.2e1;
t369 = qJD(6) * t286;
t474 = -t22 * t369 - t23 * t370;
t500 = t45 * mrSges(5,1) - t44 * mrSges(5,2) + mrSges(7,3) * t474 + Ifges(5,5) * t109 + Ifges(5,6) * t110 + Ifges(5,3) * t272 - t3 * t419 + (Ifges(6,1) * t450 + Ifges(6,4) * t451 + Ifges(6,5) * t444 + t22 * t418 + t23 * t419 + t321 * t452 + t322 * t455 + t323 * t454 + t467) * t488 + (-Ifges(6,4) * t450 + Ifges(7,5) * t454 - Ifges(6,2) * t451 - Ifges(6,6) * t444 + Ifges(7,6) * t455 + Ifges(7,3) * t452 + t466) * t487 + t293;
t499 = -m(7) - m(6);
t476 = t252 * pkin(5) + t251 * pkin(10);
t498 = m(7) * t476;
t436 = pkin(4) * t181;
t392 = mrSges(6,1) * t269 + mrSges(7,1) * t120 - mrSges(7,2) * t121 - mrSges(6,3) * t487;
t494 = -t252 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t251;
t257 = sin(t270);
t258 = cos(t270);
t421 = mrSges(6,2) * t252;
t493 = mrSges(5,1) * t257 + mrSges(6,1) * t251 + mrSges(5,2) * t258 + t421;
t188 = t228 * t288 - t229 * t283;
t189 = t228 * t283 + t229 * t288;
t144 = t188 * t282 + t189 * t287;
t141 = qJD(4) * t188 + t219 * t288 - t220 * t283;
t142 = -qJD(4) * t189 - t219 * t283 - t220 * t288;
t314 = t287 * t188 - t189 * t282;
t75 = qJD(5) * t314 + t141 * t287 + t142 * t282;
t311 = t144 * t369 + t281 * t75;
t285 = sin(qJ(1));
t290 = cos(qJ(1));
t491 = g(1) * t290 + g(2) * t285;
t486 = -t3 * mrSges(7,1) + t2 * mrSges(7,2);
t273 = -pkin(8) - t424;
t348 = m(3) * qJ(2) + mrSges(3,3);
t485 = -m(4) * t424 + m(5) * t273 + mrSges(2,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - t348;
t268 = cos(t276);
t255 = pkin(3) * t268;
t232 = t255 + t259;
t267 = sin(t276);
t329 = mrSges(4,1) * t268 - mrSges(4,2) * t267;
t330 = -mrSges(3,1) * t279 + mrSges(3,2) * t278;
t341 = t258 * mrSges(5,1) - mrSges(5,2) * t257;
t484 = m(3) * pkin(1) + m(4) * t259 + m(5) * t232 + mrSges(2,1) + t329 - t330 + t341 - t494;
t483 = -m(6) * t59 + m(7) * t57 - t392;
t123 = -mrSges(6,2) * t269 + mrSges(6,3) * t488;
t83 = -mrSges(7,2) * t130 + mrSges(7,3) * t120;
t84 = mrSges(7,1) * t130 - mrSges(7,3) * t121;
t318 = -t281 * t84 + t286 * t83;
t477 = -t123 - t318;
t193 = -t289 * t240 - t241 * t284;
t173 = -pkin(8) * t229 + t193;
t194 = -t284 * t240 + t289 * t241;
t174 = pkin(8) * t228 + t194;
t126 = t283 * t173 + t288 * t174;
t475 = t252 * t492 + t494;
t471 = -t341 + t475;
t20 = mrSges(7,1) * t54 - mrSges(7,3) * t41;
t21 = -mrSges(7,2) * t54 + mrSges(7,3) * t42;
t469 = -t281 * t20 + t286 * t21 - t84 * t369 - t83 * t370;
t320 = t22 * t286 + t23 * t281;
t428 = t3 * t281;
t298 = -qJD(6) * t320 - t428;
t429 = t2 * t286;
t468 = m(7) * (t298 + t429) + t469;
t463 = m(7) * pkin(5);
t453 = t121 / 0.2e1;
t449 = -t335 / 0.2e1;
t448 = -t181 / 0.2e1;
t447 = t181 / 0.2e1;
t445 = t218 / 0.2e1;
t443 = -t277 / 0.2e1;
t440 = pkin(3) * t218;
t439 = pkin(3) * t220;
t438 = pkin(3) * t267;
t437 = pkin(3) * t288;
t435 = pkin(4) * t257;
t250 = pkin(4) * t258;
t434 = pkin(4) * t282;
t433 = pkin(4) * t287;
t432 = pkin(5) * t251;
t417 = Ifges(4,4) * t218;
t416 = Ifges(5,4) * t181;
t407 = t114 * mrSges(5,3);
t406 = t115 * mrSges(5,3);
t404 = t190 * mrSges(4,3);
t403 = t191 * mrSges(4,3);
t391 = t144 * t281;
t390 = t144 * t286;
t383 = t281 * t285;
t382 = t281 * t290;
t381 = t282 * t283;
t380 = t283 * t287;
t379 = t285 * t286;
t378 = t286 * t290;
t117 = t288 * t159 - t155;
t263 = pkin(4) + t437;
t215 = pkin(3) * t380 + t282 * t263;
t366 = qJDD(1) * t278;
t365 = qJDD(1) * t279;
t364 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t54;
t360 = t251 * t423;
t352 = t250 + t476;
t345 = -t56 * mrSges(6,1) + t55 * mrSges(6,2);
t343 = -t187 * mrSges(4,1) + t186 * mrSges(4,2);
t342 = -t110 * mrSges(5,1) + t109 * mrSges(5,2);
t116 = -t159 * t283 - t157;
t125 = t288 * t173 - t174 * t283;
t334 = t501 * t285;
t333 = t501 * t290;
t122 = -pkin(4) * t142 + t439;
t151 = t440 + t436;
t331 = -mrSges(3,1) * t365 + mrSges(3,2) * t366;
t319 = -t22 * t281 + t23 * t286;
t98 = -pkin(9) * t189 + t125;
t99 = pkin(9) * t188 + t126;
t73 = t282 * t98 + t287 * t99;
t201 = -pkin(3) * t228 - t259;
t154 = -pkin(4) * t188 + t201;
t81 = -pkin(5) * t314 - pkin(10) * t144 + t154;
t34 = t281 * t81 + t286 * t73;
t33 = -t281 * t73 + t286 * t81;
t72 = t282 * t99 - t287 * t98;
t214 = -pkin(3) * t381 + t263 * t287;
t312 = t116 - t481;
t310 = t144 * t370 - t286 * t75;
t165 = -t240 * t375 + qJD(2) * t384 + (-qJD(2) * t278 - qJD(3) * t241) * t284;
t149 = -pkin(8) * t220 + t165;
t166 = -t229 * qJD(2) - qJD(3) * t194;
t150 = -pkin(8) * t219 + t166;
t78 = t288 * t149 + t283 * t150 + t173 * t373 - t174 * t374;
t221 = -t435 - t438;
t299 = m(7) * (t221 - t432) - t360;
t297 = m(7) * (-t432 - t435) - t360;
t296 = t421 + (mrSges(6,1) + t423 + t463) * t251;
t79 = -qJD(4) * t126 - t149 * t283 + t288 * t150;
t294 = -pkin(9) * t141 + t79;
t271 = -pkin(9) + t273;
t265 = -qJDD(1) * pkin(1) + qJDD(2);
t262 = -pkin(5) - t433;
t212 = -pkin(5) - t214;
t210 = Ifges(4,4) * t217;
t206 = t252 * t378 + t383;
t205 = -t252 * t382 + t379;
t204 = -t252 * t379 + t382;
t203 = t252 * t383 + t378;
t202 = t232 + t250;
t197 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t218;
t196 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t217;
t177 = t218 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t210;
t176 = t217 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t417;
t175 = Ifges(5,4) * t335;
t163 = mrSges(5,1) * t277 - mrSges(5,3) * t181;
t162 = -mrSges(5,2) * t277 + mrSges(5,3) * t335;
t137 = -mrSges(5,1) * t335 + mrSges(5,2) * t181;
t128 = t181 * Ifges(5,1) + t277 * Ifges(5,5) + t175;
t127 = Ifges(5,2) * t335 + t277 * Ifges(5,6) + t416;
t101 = -mrSges(5,2) * t272 + mrSges(5,3) * t110;
t100 = mrSges(5,1) * t272 - mrSges(5,3) * t109;
t96 = t117 - t497;
t87 = -mrSges(6,1) * t488 + mrSges(6,2) * t487;
t80 = t436 + t88;
t77 = t151 + t88;
t76 = qJD(5) * t144 + t141 * t282 - t287 * t142;
t66 = pkin(9) * t142 + t78;
t65 = t282 * t312 + t287 * t96;
t62 = t287 * t94 - t397;
t61 = t282 * t94 + t393;
t48 = -mrSges(6,2) * t266 + mrSges(6,3) * t56;
t47 = mrSges(6,1) * t266 - mrSges(6,3) * t55;
t30 = t281 * t88 + t286 * t59;
t29 = -t281 * t59 + t286 * t88;
t28 = t281 * t80 + t286 * t62;
t27 = -t281 * t62 + t286 * t80;
t26 = t281 * t77 + t286 * t65;
t25 = -t281 * t65 + t286 * t77;
t24 = pkin(5) * t76 - pkin(10) * t75 + t122;
t17 = -qJD(5) * t72 + t282 * t294 + t287 * t66;
t16 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t5 = -qJD(6) * t34 - t17 * t281 + t24 * t286;
t4 = qJD(6) * t33 + t17 * t286 + t24 * t281;
t1 = [(-t206 * mrSges(7,1) - t205 * mrSges(7,2) + t499 * (t290 * t202 - t271 * t285) + t485 * t285 + (-t484 - t498) * t290) * g(2) + (-t204 * mrSges(7,1) - t203 * mrSges(7,2) + (-t271 * t499 + t485) * t290 + (m(6) * t202 - m(7) * (-t202 - t476) + t484) * t285) * g(1) - t76 * t495 + t487 * (Ifges(6,1) * t75 - Ifges(6,4) * t76) / 0.2e1 + t488 * (Ifges(6,4) * t75 - Ifges(6,2) * t76) / 0.2e1 + (-m(6) * t11 + m(7) * t8 + t16 - t47) * t72 + t217 * (Ifges(4,4) * t219 - Ifges(4,2) * t220) / 0.2e1 + qJD(3) * (Ifges(4,5) * t219 - Ifges(4,6) * t220) / 0.2e1 + t239 * (mrSges(4,1) * t220 + mrSges(4,2) * t219) + (Ifges(4,1) * t219 - Ifges(4,4) * t220) * t445 + m(5) * (t114 * t79 + t115 * t78 + t125 * t45 + t126 * t44 + t161 * t201 + t192 * t439) + m(7) * (t2 * t34 + t22 * t5 + t23 * t4 + t3 * t33) + m(6) * (t10 * t73 + t122 * t145 + t154 * t89 + t17 * t60) + t483 * (qJD(5) * t73 + t282 * t66 - t287 * t294) + t219 * t177 / 0.2e1 - t220 * t176 / 0.2e1 + t166 * t197 + (mrSges(5,2) * t161 - mrSges(5,3) * t45 + Ifges(5,1) * t109 + Ifges(5,4) * t110 + Ifges(5,5) * t272) * t189 + (-Ifges(7,3) * t460 - Ifges(7,6) * t461 - Ifges(7,5) * t462 + Ifges(6,2) * t56 + Ifges(6,4) * t55 - t89 * mrSges(6,1) + Ifges(6,6) * t266 - t364 / 0.2e1 + t10 * mrSges(6,3) + t486) * t314 + m(4) * (t138 * t194 + t139 * t193 + t165 * t191 + t166 * t190 - t236 * t259) + t265 * t330 + t192 * (-mrSges(5,1) * t142 + mrSges(5,2) * t141) + t193 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t186) + t194 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t187) + t165 * t196 - t220 * t403 + t142 * t406 + (-t2 * t391 + t22 * t310 - t23 * t311 - t3 * t390) * mrSges(7,3) + (mrSges(4,2) * t236 - mrSges(4,3) * t139 + Ifges(4,1) * t186 + Ifges(4,4) * t187 + Ifges(4,5) * qJDD(3)) * t229 - pkin(1) * t331 + t335 * (Ifges(5,4) * t141 + Ifges(5,2) * t142) / 0.2e1 + 0.2e1 * t376 * t248 * mrSges(3,3) + (Ifges(3,4) * t278 + Ifges(3,2) * t279) * t365 + (Ifges(3,1) * t278 + Ifges(3,4) * t279) * t366 + (t89 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t55 + Ifges(6,4) * t56 + Ifges(6,5) * t266 + t321 * t460 + t322 * t461 + t323 * t462 + t324 * t8 + t344 * t71) * t144 + (Ifges(5,1) * t141 + Ifges(5,4) * t142) * t447 + (-Ifges(7,1) * t310 - Ifges(7,4) * t311 + Ifges(7,5) * t76) * t453 + (-mrSges(5,1) * t161 + mrSges(5,3) * t44 + Ifges(5,4) * t109 + Ifges(5,2) * t110 + Ifges(5,6) * t272) * t188 + t78 * t162 + t79 * t163 + t145 * (mrSges(6,1) * t76 + mrSges(6,2) * t75) + t142 * t127 / 0.2e1 + Ifges(2,3) * qJDD(1) + t141 * t128 / 0.2e1 + t125 * t100 + t126 * t101 + t122 * t87 + t17 * t123 - t76 * t85 / 0.2e1 + t75 * t86 / 0.2e1 + t4 * t83 + t5 * t84 + t76 * t69 / 0.2e1 + m(3) * (-pkin(1) * t265 + (t248 + t367) * qJ(2) * t376) + t73 * t48 + t137 * t439 - t76 * t426 + t33 * t20 + t34 * t21 - t75 * t427 - t219 * t404 - t141 * t407 - t14 * t391 / 0.2e1 + t15 * t390 / 0.2e1 - t311 * t70 / 0.2e1 + t154 * t345 + t201 * t342 - t259 * t343 + t269 * (Ifges(6,5) * t75 - Ifges(6,6) * t76) / 0.2e1 + t277 * (Ifges(5,5) * t141 + Ifges(5,6) * t142) / 0.2e1 + t57 * (mrSges(7,1) * t311 - mrSges(7,2) * t310) + t130 * (-Ifges(7,5) * t310 - Ifges(7,6) * t311 + Ifges(7,3) * t76) / 0.2e1 + t120 * (-Ifges(7,4) * t310 - Ifges(7,2) * t311 + Ifges(7,6) * t76) / 0.2e1 + t75 * t353 + t76 * t496 + (-mrSges(4,1) * t236 + mrSges(4,3) * t138 + Ifges(4,4) * t186 + Ifges(4,2) * t187 + Ifges(4,6) * qJDD(3)) * t228; t342 + t343 + t218 * t197 - t217 * t196 + t331 + t477 * t488 + t318 * qJD(6) + t345 + m(3) * t265 + t181 * t163 - t335 * t162 + t392 * t487 + t281 * t21 + t286 * t20 + (-g(1) * t285 + g(2) * t290) * (m(3) + m(4) + m(5) - t499) - t348 * t376 * qJD(1) ^ 2 + (t130 * t319 + t2 * t281 + t3 * t286 - t487 * t57) * m(7) + (t487 * t59 - t488 * t60 + t89) * m(6) + (t114 * t181 - t115 * t335 + t161) * m(5) + (t190 * t218 - t191 * t217 + t236) * m(4); -(t192 * mrSges(5,1) + Ifges(5,6) * t443 + Ifges(5,4) * t448 + Ifges(5,2) * t449 - t127 / 0.2e1 - t406) * t181 + (t212 * t8 - t22 * t25 - t23 * t26) * m(7) + (t10 * t215 + t11 * t214 - t145 * t151 - t60 * t65) * m(6) + (m(6) * t60 + m(7) * t319 - t477) * (t263 * t371 + (-t283 * t372 + (t287 * t288 - t381) * qJD(4)) * pkin(3)) + t468 * (pkin(10) + t215) + t483 * (-t282 * t96 + t287 * t312 + t263 * t372 + (t283 * t371 + (t282 * t288 + t380) * qJD(4)) * pkin(3)) - qJD(3) * (Ifges(4,5) * t217 - Ifges(4,6) * t218) / 0.2e1 + t212 * t16 + t214 * t47 + t215 * t48 + t191 * t197 + t100 * t437 - t190 * t196 + t218 * t403 + t217 * t404 - t239 * (mrSges(4,1) * t218 + mrSges(4,2) * t217) - g(1) * (t290 * t299 + t333) - g(2) * (t285 * t299 + t334) - (-Ifges(4,2) * t218 + t177 + t210) * t217 / 0.2e1 + (-t192 * mrSges(5,2) + t407 + Ifges(5,5) * t443 + Ifges(5,1) * t448 + Ifges(5,4) * t449 - t128 / 0.2e1) * t335 + Ifges(4,5) * t186 + Ifges(4,6) * t187 + t500 + t176 * t445 - t117 * t162 - t116 * t163 - t151 * t87 - t138 * mrSges(4,2) + t139 * mrSges(4,1) + (m(5) * (t283 * t44 + t288 * t45 + (-t114 * t283 + t115 * t288) * qJD(4)) - t163 * t374 + t283 * t101 + t162 * t373) * pkin(3) - t65 * t123 - t26 * t83 - t25 * t84 - t137 * t440 - m(5) * (t114 * t116 + t115 * t117 + t192 * t440) + Ifges(4,3) * qJDD(3) - t218 * (Ifges(4,1) * t217 - t417) / 0.2e1 + t491 * (m(5) * t438 - m(6) * t221 + mrSges(4,1) * t267 + mrSges(4,2) * t268 + t493) + (-t329 - m(5) * t255 - m(7) * (t255 + t352) - m(6) * (t250 + t255) + t471) * g(3); -t477 * pkin(4) * t371 + t468 * (pkin(10) + t434) + t47 * t433 + t48 * t434 + t181 * t406 + ((t10 * t282 + t11 * t287 + (-t282 * t59 + t287 * t60) * qJD(5)) * pkin(4) - t145 * t436 + t59 * t61 - t60 * t62) * m(6) - g(1) * (t290 * t297 + t333) - g(2) * (t285 * t297 + t334) + t392 * (-pkin(4) * t372 + t61) - t192 * (mrSges(5,1) * t181 + mrSges(5,2) * t335) + t335 * t407 + (Ifges(5,5) * t335 - Ifges(5,6) * t181) * t443 + (Ifges(5,1) * t335 - t416) * t448 + t500 + (-Ifges(5,2) * t181 + t128 + t175) * t449 + t127 * t447 + (t262 * t8 + (t282 * t57 + t287 * t319) * qJD(5) * pkin(4) - t22 * t27 - t23 * t28 - t57 * t61) * m(7) - t114 * t162 + t115 * t163 - t62 * t123 - t28 * t83 - t27 * t84 - t87 * t436 + (m(6) * t435 + t493) * t491 + (-m(6) * t250 - m(7) * t352 + t471) * g(3) + t262 * t16; t293 + (t410 / 0.2e1 - t408 / 0.2e1 - t409 / 0.2e1 - t411 / 0.2e1 + t415 / 0.2e1 + t466) * t487 + t392 * t60 + t298 * mrSges(7,3) + (-t129 / 0.2e1 - t412 / 0.2e1 - t306 / 0.2e1 - t308 / 0.2e1 - t307 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t487 + t320 * mrSges(7,3) + t467) * t488 - t8 * t463 + (t290 * t296 - t333) * g(1) + (t285 * t296 - t334) * g(2) - m(7) * (t22 * t29 + t23 * t30 + t57 * t60) + (t475 - t498) * g(3) - t59 * t123 - t30 * t83 - t29 * t84 - pkin(5) * t16 + (m(7) * (-t428 + t429 + t474) + t469) * pkin(10); -t57 * (mrSges(7,1) * t121 + mrSges(7,2) * t120) + (Ifges(7,1) * t120 - t405) * t454 + t70 * t453 + (Ifges(7,5) * t120 - Ifges(7,6) * t121) * t452 - t22 * t83 + t23 * t84 - g(1) * (mrSges(7,1) * t205 - mrSges(7,2) * t206) - g(2) * (-mrSges(7,1) * t203 + mrSges(7,2) * t204) + g(3) * t324 * t251 + (t120 * t22 + t121 * t23) * mrSges(7,3) + t364 + (-Ifges(7,2) * t121 + t119 + t71) * t455 - t486;];
tau  = t1;
