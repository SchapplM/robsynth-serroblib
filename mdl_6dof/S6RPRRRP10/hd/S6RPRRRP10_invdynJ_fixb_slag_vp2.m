% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:19
% EndTime: 2019-03-09 06:31:08
% DurationCPUTime: 31.19s
% Computational Cost: add. (9225->721), mult. (18085->923), div. (0->0), fcn. (11472->10), ass. (0->341)
t520 = -Ifges(4,2) / 0.2e1;
t262 = sin(qJ(3));
t252 = t262 * qJD(1);
t354 = qJD(4) + qJD(5);
t231 = t252 + t354;
t417 = -t231 / 0.2e1;
t494 = Ifges(7,2) + Ifges(6,3);
t519 = t494 * t417;
t266 = cos(qJ(3));
t356 = qJD(1) * qJD(3);
t204 = qJDD(1) * t266 - t262 * t356;
t261 = sin(qJ(4));
t265 = cos(qJ(4));
t366 = qJD(3) * t265;
t368 = qJD(1) * t266;
t287 = t261 * t368 - t366;
t281 = t287 * qJD(4);
t107 = qJDD(3) * t261 + t204 * t265 - t281;
t194 = qJD(3) * t261 + t265 * t368;
t108 = -qJD(4) * t194 + qJDD(3) * t265 - t204 * t261;
t260 = sin(qJ(5));
t264 = cos(qJ(5));
t115 = t260 * t194 + t264 * t287;
t36 = -qJD(5) * t115 + t264 * t107 + t260 * t108;
t518 = -t36 / 0.2e1;
t276 = t264 * t194 - t260 * t287;
t37 = qJD(5) * t276 + t260 * t107 - t264 * t108;
t440 = -t37 / 0.2e1;
t429 = t115 / 0.2e1;
t420 = -t194 / 0.2e1;
t239 = t252 + qJD(4);
t517 = -t239 / 0.2e1;
t473 = t287 / 0.2e1;
t404 = Ifges(4,4) * t266;
t516 = t262 * t520 + t404 / 0.2e1;
t469 = Ifges(7,4) + Ifges(6,5);
t515 = -t469 / 0.2e1;
t416 = t231 / 0.2e1;
t426 = t276 / 0.2e1;
t430 = -t115 / 0.2e1;
t470 = Ifges(6,1) + Ifges(7,1);
t269 = -pkin(1) - pkin(7);
t233 = qJD(1) * t269 + qJD(2);
t375 = t266 * t233;
t182 = -qJD(3) * pkin(3) - t375;
t126 = pkin(4) * t287 + t182;
t316 = pkin(3) * t262 - pkin(8) * t266;
t210 = qJ(2) + t316;
t173 = t210 * qJD(1);
t209 = t262 * t233;
t181 = qJD(3) * pkin(8) + t209;
t101 = t261 * t173 + t265 * t181;
t86 = -pkin(9) * t287 + t101;
t398 = t260 * t86;
t100 = t265 * t173 - t181 * t261;
t85 = -pkin(9) * t194 + t100;
t76 = pkin(4) * t239 + t85;
t29 = t264 * t76 - t398;
t459 = qJD(6) - t29;
t27 = -pkin(5) * t231 + t459;
t42 = t115 * pkin(5) - qJ(6) * t276 + t126;
t110 = Ifges(6,4) * t115;
t399 = Ifges(7,5) * t115;
t463 = t231 * t469 + t276 * t470 - t110 + t399;
t503 = t126 * mrSges(6,2) + t27 * mrSges(7,2) - t29 * mrSges(6,3) - t42 * mrSges(7,3) + t463 / 0.2e1;
t514 = Ifges(6,4) * t430 + Ifges(7,5) * t429 + t469 * t416 + t470 * t426 + t503;
t497 = mrSges(6,1) + mrSges(7,1);
t496 = mrSges(6,2) - mrSges(7,3);
t468 = Ifges(7,5) - Ifges(6,4);
t467 = Ifges(6,6) - Ifges(7,6);
t342 = t261 * t252;
t363 = qJD(4) * t261;
t513 = t363 + t342;
t395 = t264 * t86;
t30 = t260 * t76 + t395;
t28 = qJ(6) * t231 + t30;
t109 = Ifges(7,5) * t276;
t51 = Ifges(7,6) * t231 + Ifges(7,3) * t115 + t109;
t400 = Ifges(6,4) * t276;
t54 = -Ifges(6,2) * t115 + Ifges(6,6) * t231 + t400;
t512 = t126 * mrSges(6,1) + t42 * mrSges(7,1) - t28 * mrSges(7,2) - t30 * mrSges(6,3) - t54 / 0.2e1 + t51 / 0.2e1;
t205 = -t262 * qJDD(1) - t266 * t356;
t357 = qJD(1) * qJD(2);
t235 = qJDD(1) * qJ(2) + t357;
t112 = -pkin(3) * t205 - pkin(8) * t204 + t235;
t226 = qJDD(1) * t269 + qJDD(2);
t365 = qJD(3) * t266;
t132 = t262 * t226 + t233 * t365;
t128 = qJDD(3) * pkin(8) + t132;
t362 = qJD(4) * t265;
t40 = t261 * t112 + t265 * t128 + t173 * t362 - t181 * t363;
t41 = -qJD(4) * t101 + t265 * t112 - t128 * t261;
t511 = t41 * mrSges(5,1) - t40 * mrSges(5,2);
t191 = qJDD(4) - t205;
t185 = qJDD(5) + t191;
t20 = pkin(4) * t191 - pkin(9) * t107 + t41;
t22 = pkin(9) * t108 + t40;
t359 = qJD(5) * t264;
t360 = qJD(5) * t260;
t5 = t260 * t20 + t264 * t22 + t76 * t359 - t360 * t86;
t2 = qJ(6) * t185 + qJD(6) * t231 + t5;
t6 = -qJD(5) * t30 + t20 * t264 - t22 * t260;
t3 = -pkin(5) * t185 + qJDD(6) - t6;
t510 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t508 = t27 * mrSges(7,1) + t30 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t516 + Ifges(5,5) * t420 + Ifges(5,6) * t473 + Ifges(5,3) * t517 + t467 * t429 + t519 + t276 * t515 - t28 * mrSges(7,3) - t29 * mrSges(6,1);
t422 = t185 / 0.2e1;
t439 = t37 / 0.2e1;
t441 = t36 / 0.2e1;
t367 = qJD(3) * t262;
t131 = t226 * t266 - t233 * t367;
t127 = -qJDD(3) * pkin(3) - t131;
t69 = -pkin(4) * t108 + t127;
t7 = pkin(5) * t37 - qJ(6) * t36 - qJD(6) * t276 + t69;
t507 = -mrSges(6,2) * t69 - mrSges(7,2) * t3 + mrSges(6,3) * t6 + mrSges(7,3) * t7 - Ifges(7,5) * t439 + t185 * t515 - t422 * t469 + (-t441 + t518) * t470 + (-Ifges(6,4) + t468) * t440;
t506 = m(6) + m(7);
t505 = t266 / 0.2e1;
t427 = -t276 / 0.2e1;
t504 = Ifges(6,2) * t429 - Ifges(7,3) * t430 + t467 * t417 - t468 * t427 + t512;
t502 = mrSges(6,1) * t69 + mrSges(7,1) * t7 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + Ifges(6,4) * t518 - t185 * Ifges(6,6) / 0.2e1 + 0.2e1 * Ifges(7,3) * t439 + (-t440 + t439) * Ifges(6,2) + (t468 + Ifges(7,5)) * t441 + (-t467 + Ifges(7,6)) * t422;
t501 = -Ifges(6,2) * t430 + Ifges(7,3) * t429 - t416 * t467 + t426 * t468 + t512;
t500 = -m(4) - m(5);
t495 = mrSges(6,3) + mrSges(7,2);
t59 = -mrSges(5,1) * t108 + mrSges(5,2) * t107;
t491 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t204 - t59;
t387 = t260 * t265;
t198 = t261 * t264 + t387;
t490 = t354 * t198;
t405 = Ifges(4,4) * t262;
t309 = t266 * Ifges(4,1) - t405;
t187 = Ifges(5,4) * t287;
t98 = t194 * Ifges(5,1) + t239 * Ifges(5,5) - t187;
t489 = Ifges(4,5) * qJD(3) + qJD(1) * t309 + t265 * t98;
t488 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t287 + t194 * mrSges(5,2) + mrSges(4,3) * t368;
t259 = qJ(4) + qJ(5);
t253 = sin(t259);
t254 = cos(t259);
t487 = -t496 * t253 + t254 * t497;
t341 = t261 * t367;
t361 = qJD(4) * t266;
t486 = t265 * t361 - t341;
t263 = sin(qJ(1));
t380 = t263 * t265;
t267 = cos(qJ(1));
t382 = t262 * t267;
t176 = t261 * t382 + t380;
t453 = pkin(4) * t513 - t209;
t485 = t185 * t494 + t36 * t469 - t37 * t467;
t295 = t131 * t266 + t132 * t262;
t315 = mrSges(4,1) * t266 - mrSges(4,2) * t262;
t483 = qJ(2) * t315 + (-Ifges(4,1) * t262 - t404) * t505;
t482 = -Ifges(6,4) * t429 - Ifges(7,5) * t430 - t417 * t469 - t427 * t470 + t503;
t481 = mrSges(3,2) - mrSges(4,3) - mrSges(2,1);
t65 = pkin(5) * t276 + qJ(6) * t115;
t314 = t262 * mrSges(4,1) + t266 * mrSges(4,2);
t479 = -m(5) * t316 + t266 * mrSges(5,3) + mrSges(2,2) - mrSges(3,3) - t314;
t478 = m(7) * pkin(5) + t497;
t477 = -m(7) * qJ(6) + t496;
t432 = t107 / 0.2e1;
t431 = t108 / 0.2e1;
t421 = t191 / 0.2e1;
t268 = -pkin(9) - pkin(8);
t371 = t267 * t268;
t471 = g(2) * t262 * t371;
t413 = g(2) * t267;
t23 = mrSges(6,1) * t185 - mrSges(6,3) * t36;
t24 = -t185 * mrSges(7,1) + t36 * mrSges(7,2);
t466 = t24 - t23;
t25 = -mrSges(6,2) * t185 - mrSges(6,3) * t37;
t26 = -mrSges(7,2) * t37 + mrSges(7,3) * t185;
t465 = t25 + t26;
t197 = t260 * t261 - t264 * t265;
t120 = t354 * t197;
t172 = t198 * qJD(1);
t146 = t262 * t172;
t460 = t197 * t262;
t147 = qJD(1) * t460;
t464 = -qJD(6) * t198 + t453 + (t120 + t147) * qJ(6) + (t490 + t146) * pkin(5);
t88 = -mrSges(7,2) * t115 + mrSges(7,3) * t231;
t407 = mrSges(6,3) * t115;
t89 = -mrSges(6,2) * t231 - t407;
t409 = t88 + t89;
t406 = mrSges(6,3) * t276;
t90 = mrSges(6,1) * t231 - t406;
t91 = -mrSges(7,1) * t231 + mrSges(7,2) * t276;
t408 = -t91 + t90;
t160 = t197 * t266;
t462 = -qJD(3) * t160 - t262 * t490 - t172;
t321 = t354 * t265;
t335 = t261 * t360;
t338 = t262 * t363;
t461 = -t365 * t387 + t260 * t338 + t262 * t335 - (t261 * t365 + t262 * t321) * t264 + t197 * qJD(1);
t190 = t265 * t210;
t327 = -t261 * t269 + pkin(4);
t377 = t265 * t266;
t113 = -pkin(9) * t377 + t262 * t327 + t190;
t381 = t262 * t269;
t228 = t265 * t381;
t137 = t261 * t210 + t228;
t386 = t261 * t266;
t125 = -pkin(9) * t386 + t137;
t458 = t260 * t113 + t264 * t125;
t248 = pkin(4) * t265 + pkin(3);
t300 = pkin(5) * t254 + qJ(6) * t253;
t456 = -m(7) * (-t248 - t300) + m(6) * t248 + t487;
t454 = t268 * t506 - t495;
t283 = t261 * t361 + t262 * t366;
t80 = mrSges(5,1) * t191 - mrSges(5,3) * t107;
t81 = -mrSges(5,2) * t191 + mrSges(5,3) * t108;
t452 = -t261 * t80 + t265 * t81;
t451 = -t261 * t41 + t265 * t40;
t414 = g(1) * t263;
t450 = t413 - t414;
t317 = pkin(3) * t266 + pkin(8) * t262;
t192 = qJD(3) * t317 + qJD(2);
t164 = t265 * t192;
t383 = t262 * t265;
t353 = pkin(9) * t383;
t61 = t164 + (-t228 + (pkin(9) * t266 - t210) * t261) * qJD(4) + (t266 * t327 + t353) * qJD(3);
t364 = qJD(3) * t269;
t339 = t266 * t364;
t83 = t261 * t192 + t210 * t362 + t265 * t339 - t269 * t338;
t68 = -pkin(9) * t486 + t83;
t17 = -qJD(5) * t458 - t260 * t68 + t264 * t61;
t313 = -mrSges(5,1) * t265 + mrSges(5,2) * t261;
t285 = m(5) * pkin(3) - t313;
t348 = m(5) * pkin(8) + mrSges(5,3);
t449 = t262 * t348 + t266 * t285;
t444 = qJD(1) ^ 2;
t438 = Ifges(5,1) * t432 + Ifges(5,4) * t431 + Ifges(5,5) * t421;
t433 = -m(3) - m(4);
t419 = t194 / 0.2e1;
t412 = g(3) * t266;
t411 = -qJD(1) / 0.2e1;
t410 = qJD(4) / 0.2e1;
t201 = t317 * qJD(1);
t122 = t265 * t201 - t261 * t375;
t92 = (pkin(4) * t266 + t353) * qJD(1) + t122;
t123 = t261 * t201 + t265 * t375;
t99 = pkin(9) * t342 + t123;
t48 = t260 * t92 + t264 * t99;
t403 = Ifges(5,4) * t194;
t402 = Ifges(5,4) * t261;
t401 = Ifges(5,4) * t265;
t389 = t254 * t266;
t385 = t261 * t267;
t384 = t262 * t263;
t379 = t263 * t266;
t376 = t265 * t267;
t373 = t266 * t269;
t372 = t267 * t254;
t370 = t176 * pkin(4);
t369 = t267 * pkin(1) + t263 * qJ(2);
t358 = qJDD(1) * mrSges(3,2);
t352 = g(1) * t384;
t346 = Ifges(5,5) * t107 + Ifges(5,6) * t108 + Ifges(5,3) * t191;
t344 = t267 * pkin(7) + t369;
t343 = qJD(4) * t268;
t238 = t262 * t364;
t326 = -t356 / 0.2e1;
t324 = (t235 + t357) * qJ(2);
t152 = t253 * t384 - t372;
t153 = t253 * t267 + t254 * t384;
t323 = -t152 * pkin(5) + qJ(6) * t153;
t154 = t253 * t382 + t254 * t263;
t155 = -t253 * t263 + t262 * t372;
t322 = t154 * pkin(5) - qJ(6) * t155;
t240 = pkin(4) * t386;
t193 = t240 - t373;
t320 = t152 * t497 + t153 * t496;
t319 = -t154 * t497 - t155 * t496;
t318 = t265 * t343;
t312 = mrSges(5,1) * t261 + mrSges(5,2) * t265;
t308 = Ifges(5,1) * t265 - t402;
t307 = Ifges(5,1) * t261 + t401;
t305 = -Ifges(5,2) * t261 + t401;
t304 = Ifges(5,2) * t265 + t402;
t303 = -Ifges(4,5) * t262 - Ifges(4,6) * t266;
t302 = Ifges(5,5) * t265 - Ifges(5,6) * t261;
t301 = Ifges(5,5) * t261 + Ifges(5,6) * t265;
t47 = -t260 * t99 + t264 * t92;
t297 = t100 * t265 + t101 * t261;
t63 = t113 * t264 - t125 * t260;
t134 = -t239 * mrSges(5,2) - mrSges(5,3) * t287;
t135 = mrSges(5,1) * t239 - mrSges(5,3) * t194;
t294 = -t261 * t134 - t265 * t135;
t222 = t268 * t261;
t223 = t268 * t265;
t293 = t264 * t222 + t223 * t260;
t130 = t222 * t260 - t223 * t264;
t16 = t113 * t359 - t125 * t360 + t260 * t61 + t264 * t68;
t289 = t262 * (-Ifges(4,2) * t266 - t405);
t133 = pkin(4) * t486 + t238;
t279 = Ifges(5,5) * t266 - t262 * t308;
t278 = Ifges(5,6) * t266 - t262 * t305;
t277 = Ifges(5,3) * t266 - t262 * t302;
t275 = t485 + t510;
t274 = (mrSges(6,2) * t254 + t253 * t478) * t266 - mrSges(7,3) * t389;
t273 = -qJD(4) * t297 + t451;
t256 = t267 * qJ(2);
t249 = -pkin(1) * qJDD(1) + qJDD(2);
t247 = -pkin(4) * t264 - pkin(5);
t244 = pkin(4) * t260 + qJ(6);
t243 = pkin(4) * t376;
t236 = pkin(4) * t359 + qJD(6);
t224 = qJ(6) * t389;
t218 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t252;
t207 = t248 * t379;
t202 = t261 * t343;
t200 = t314 * qJD(1);
t186 = t312 * t266;
t177 = -t261 * t263 + t262 * t376;
t175 = t262 * t380 + t385;
t174 = -t261 * t384 + t376;
t162 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t205;
t158 = t198 * t266;
t157 = t198 * t262;
t136 = -t261 * t381 + t190;
t111 = pkin(5) * t197 - qJ(6) * t198 - t248;
t97 = -Ifges(5,2) * t287 + Ifges(5,6) * t239 + t403;
t84 = -qJD(4) * t137 - t261 * t339 + t164;
t82 = pkin(5) * t158 + qJ(6) * t160 + t193;
t78 = qJD(5) * t130 + t202 * t260 - t264 * t318;
t77 = qJD(5) * t293 + t264 * t202 + t260 * t318;
t75 = -t266 * t335 + (t266 * t321 - t341) * t264 - t283 * t260;
t73 = qJD(3) * t460 - t266 * t490;
t67 = mrSges(6,1) * t115 + mrSges(6,2) * t276;
t66 = mrSges(7,1) * t115 - mrSges(7,3) * t276;
t60 = -pkin(5) * t262 - t63;
t58 = qJ(6) * t262 + t458;
t50 = pkin(4) * t194 + t65;
t45 = t107 * Ifges(5,4) + t108 * Ifges(5,2) + t191 * Ifges(5,6);
t44 = -pkin(5) * t368 - t47;
t43 = qJ(6) * t368 + t48;
t39 = t264 * t85 - t398;
t38 = t260 * t85 + t395;
t18 = pkin(5) * t75 - qJ(6) * t73 + qJD(6) * t160 + t133;
t15 = -pkin(5) * t365 - t17;
t14 = mrSges(6,1) * t37 + mrSges(6,2) * t36;
t13 = mrSges(7,1) * t37 - mrSges(7,3) * t36;
t12 = qJ(6) * t365 + qJD(6) * t262 + t16;
t1 = [((t495 * t267 - t371 * t506) * g(1) + Ifges(4,5) * qJDD(3) - m(5) * t127 * t269 + t302 * t421 + t305 * t431 + t308 * t432) * t266 + (-t177 * mrSges(5,1) + t176 * mrSges(5,2) - t506 * ((-pkin(4) * t261 + t269) * t263 + t248 * t382 + t256) - t478 * t155 + t477 * t154 + (-m(3) + t500) * t256 + (m(3) * pkin(1) + t269 * t500 - t481) * t263 + t479 * t267) * g(1) + (-m(3) * t369 - t175 * mrSges(5,1) - t174 * mrSges(5,2) + t495 * t379 + t500 * t344 - t506 * (pkin(4) * t385 + t248 * t384 + t268 * t379 + t344) - t478 * t153 + t477 * t152 + t481 * t267 + t479 * t263) * g(2) + t182 * (mrSges(5,1) * t486 - mrSges(5,2) * t283) + (t100 * t283 - t101 * t486 - t377 * t41 - t386 * t40) * mrSges(5,3) + (Ifges(4,1) * t204 + Ifges(4,4) * t205) * t505 + t507 * t160 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(6,6) * t430 + Ifges(7,6) * t429 + t494 * t416 + t469 * t426 - t508) * t365 + t205 * t516 + (t314 + 0.2e1 * mrSges(3,3)) * t235 + t514 * t73 + t502 * t158 + t501 * t75 + m(5) * (t182 * t367 * t269 + t100 * t84 + t101 * t83 + t136 * t41 + t137 * t40) + (Ifges(7,6) * t439 + Ifges(6,6) * t440 - Ifges(4,4) * t204 / 0.2e1 + t205 * t520 - Ifges(4,6) * qJDD(3) + t346 / 0.2e1 + t485 / 0.2e1 + t494 * t422 + t469 * t441 + Ifges(5,3) * t421 + Ifges(5,6) * t431 + Ifges(5,5) * t432 + t510 + t511) * t262 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t58 * t26 + t60 * t24 + t63 * t23 - pkin(1) * t358 - t287 * (qJD(3) * t278 - t304 * t361) / 0.2e1 + t239 * (qJD(3) * t277 - t301 * t361) / 0.2e1 + qJD(3) ^ 2 * t303 / 0.2e1 - t45 * t386 / 0.2e1 + t204 * t309 / 0.2e1 + m(7) * (t12 * t28 + t15 * t27 + t18 * t42 + t2 * t58 + t3 * t60 + t7 * t82) + t97 * t341 / 0.2e1 - (t261 * t98 + t265 * t97) * t361 / 0.2e1 + t458 * t25 + m(6) * (t126 * t133 + t16 * t30 + t17 * t29 + t193 * t69 + t458 * t5 + t6 * t63) + t249 * mrSges(3,2) + t162 * t381 + t18 * t66 + t483 * t356 - t295 * mrSges(4,3) + t488 * t238 - t489 * t367 / 0.2e1 + t82 * t13 + t12 * t88 + t16 * t89 + t17 * t90 + t15 * t91 + t218 * t339 + t491 * t373 + t133 * t67 + t83 * t134 + t84 * t135 + t136 * t80 + t137 * t81 + t289 * t326 + t127 * t186 + t193 * t14 + qJD(2) * t200 + qJ(2) * (-mrSges(4,1) * t205 + mrSges(4,2) * t204) + m(4) * (t269 * t295 + t324) + m(3) * (-pkin(1) * t249 + t324) + (qJD(3) * t279 - t307 * t361) * t419 + t377 * t438; t358 - t465 * t460 + t466 * t157 + (qJ(2) * t433 - mrSges(3,3)) * t444 + (-t200 + t294) * qJD(1) + (-t13 - t14 + (t134 * t265 - t135 * t261 + t218) * qJD(3) + t491) * t266 + (t162 + t294 * qJD(4) + (t66 + t67 + t488) * qJD(3) + t452) * t262 + m(4) * t295 + m(3) * t249 + t462 * t409 + t461 * t408 + t450 * (m(5) - t433 + t506) + (t157 * t3 - t2 * t460 - t266 * t7 - t27 * t461 + t28 * t462 + t367 * t42) * m(7) + (t126 * t367 - t157 * t6 - t266 * t69 + t29 * t461 + t30 * t462 - t460 * t5) * m(6) + ((-t127 + (-t100 * t261 + t101 * t265) * qJD(3)) * t266 + (qJD(3) * t182 + t273) * t262 - t297 * qJD(1)) * m(5); -t488 * t209 - t513 * t97 / 0.2e1 - t507 * t198 + (t279 * t411 + t308 * t410) * t194 - t514 * t120 + (-t101 * (mrSges(5,3) * t261 * t262 - mrSges(5,2) * t266) - t100 * (mrSges(5,1) * t266 + mrSges(5,3) * t383) + t278 * t473) * qJD(1) + t502 * t197 + (t182 * t312 + t277 * t411 + t302 * t410) * t239 + t501 * t490 + t127 * t313 - pkin(3) * t59 + (Ifges(6,6) * t429 + Ifges(7,6) * t430 + t469 * t427 + t508 + t519) * t368 - t449 * t414 - t408 * t78 + t409 * t77 + t464 * t66 + t465 * t130 + (-t100 * t362 - t101 * t363 + t451) * mrSges(5,3) + (m(5) * t273 - t134 * t363 - t135 * t362 + t452) * pkin(8) + t453 * t67 - t218 * t375 + t450 * t315 + t98 * t362 / 0.2e1 - t305 * t281 / 0.2e1 + t504 * t146 + (-pkin(3) * t127 - t100 * t122 - t101 * t123 - t182 * t209) * m(5) + (-t471 + t111 * t7 - t293 * t3 + t130 * t2 + t464 * t42 + (-t43 + t77) * t28 + (-t44 + t78) * t27) * m(7) + (-g(1) * t207 - t471 + t293 * t6 + t130 * t5 - t248 * t69 + (-t48 + t77) * t30 + (-t47 - t78) * t29 + t453 * t126) * m(6) - t466 * t293 + (t314 + (-t348 + t454) * t266 + (t285 + t456) * t262) * g(3) + Ifges(4,3) * qJDD(3) + t265 * t45 / 0.2e1 - t248 * t14 - t482 * t147 + (t289 / 0.2e1 - t483) * t444 + (-m(7) * t207 + ((-m(7) * t300 - t487) * t266 + t454 * t262) * t263) * g(1) + t489 * t252 / 0.2e1 - t43 * t88 - t48 * t89 - t47 * t90 - t44 * t91 + t111 * t13 + (t262 * t495 + t456 * t266 + t449) * t413 + t131 * mrSges(4,1) - t132 * mrSges(4,2) - t123 * t134 - t122 * t135 + t303 * t326 + Ifges(4,5) * t204 + Ifges(4,6) * t205 + t301 * t421 + t304 * t431 + t307 * t432 + t261 * t438; (-Ifges(5,5) * t287 - Ifges(5,6) * t194) * t517 + ((-t322 - t370) * g(2) - t27 * t38 - t42 * t50 + (-t224 + t240) * g(3) + t2 * t244 + t247 * t3 + (-t243 - t323) * g(1) + (-t39 + t236) * t28) * m(7) + (t186 + t274) * g(3) + (-m(6) * t370 - mrSges(5,1) * t176 - mrSges(5,2) * t177 + t319) * g(2) + t275 + t408 * t38 - t409 * t39 + t346 + (-Ifges(5,2) * t194 - t187 + t98) * t473 + (-m(6) * t243 - mrSges(5,1) * t174 + mrSges(5,2) * t175 + t320) * g(1) + (-t194 * t67 + t264 * t23 + t260 * t25 + (-t260 * t408 + t264 * t89) * qJD(5) + (t261 * t352 + t27 * t360) * m(7) + ((t352 + t412) * t261 + t260 * t5 + t264 * t6 - t29 * t360 + t30 * t359 + 0.2e1 * t126 * t420) * m(6)) * pkin(4) - t504 * t276 - t182 * (t194 * mrSges(5,1) - mrSges(5,2) * t287) + (-t100 * t287 + t101 * t194) * mrSges(5,3) - m(6) * (-t29 * t38 + t30 * t39) + t247 * t24 + t236 * t88 + t244 * t26 - t50 * t66 + t482 * t115 - t100 * t134 + t101 * t135 + t511 + t97 * t419 + (-Ifges(5,1) * t287 - t403) * t420; qJ(6) * t26 - pkin(5) * t24 + t275 + (t406 + t408) * t30 + (-t407 - t409) * t29 + t274 * g(3) + (t115 * t27 + t276 * t28) * mrSges(7,2) - t65 * t66 + qJD(6) * t88 - t42 * (mrSges(7,1) * t276 + mrSges(7,3) * t115) - t126 * (mrSges(6,1) * t276 - mrSges(6,2) * t115) + t319 * g(2) + t320 * g(1) + t54 * t426 + (Ifges(7,3) * t276 - t399) * t430 + (-t115 * t469 - t276 * t467) * t417 + (-pkin(5) * t3 - t323 * g(1) - t322 * g(2) - t224 * g(3) + qJ(6) * t2 - t27 * t30 + t28 * t459 - t42 * t65) * m(7) + (-Ifges(6,2) * t276 - t110 + t463) * t429 + (-t115 * t470 + t109 - t400 + t51) * t427; t276 * t66 - t231 * t88 + (-g(1) * t152 + g(2) * t154 - t28 * t231 - t253 * t412 + t276 * t42 + t3) * m(7) + t24;];
tau  = t1;
