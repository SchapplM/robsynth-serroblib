% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:51:06
% EndTime: 2019-03-08 21:51:45
% DurationCPUTime: 23.06s
% Computational Cost: add. (12500->723), mult. (29148->995), div. (0->0), fcn. (23081->18), ass. (0->337)
t295 = sin(qJ(3));
t299 = cos(qJ(3));
t292 = -qJ(4) - pkin(8);
t351 = qJD(3) * t292;
t224 = qJD(4) * t299 + t295 * t351;
t225 = -qJD(4) * t295 + t299 * t351;
t287 = sin(pkin(12));
t290 = cos(pkin(12));
t241 = t287 * t299 + t290 * t295;
t300 = cos(qJ(2));
t289 = sin(pkin(6));
t387 = qJD(1) * t289;
t362 = t300 * t387;
t476 = -t224 * t287 + t290 * t225 + t241 * t362;
t240 = -t287 * t295 + t290 * t299;
t475 = t290 * t224 + t287 * t225 - t240 * t362;
t233 = t241 * qJD(3);
t510 = pkin(9) * t233 - t475;
t234 = t240 * qJD(3);
t509 = -pkin(9) * t234 + t476;
t508 = mrSges(6,2) - mrSges(7,3);
t293 = sin(qJ(6));
t297 = cos(qJ(6));
t342 = -mrSges(7,1) * t297 + mrSges(7,2) * t293;
t507 = -mrSges(6,1) + t342;
t294 = sin(qJ(5));
t298 = cos(qJ(5));
t256 = t292 * t295;
t257 = t292 * t299;
t183 = t290 * t256 + t257 * t287;
t127 = -pkin(9) * t241 + t183;
t184 = t287 * t256 - t290 * t257;
t128 = pkin(9) * t240 + t184;
t331 = t298 * t127 - t128 * t294;
t506 = -qJD(5) * t331 - t509 * t294 + t298 * t510;
t381 = qJD(3) * t295;
t278 = pkin(3) * t381;
t191 = pkin(4) * t233 + t278;
t296 = sin(qJ(2));
t363 = t296 * t387;
t327 = t298 * t240 - t241 * t294;
t98 = qJD(5) * t327 - t233 * t294 + t234 * t298;
t161 = t240 * t294 + t241 * t298;
t99 = qJD(5) * t161 + t298 * t233 + t234 * t294;
t505 = pkin(5) * t99 - pkin(10) * t98 + t191 - t363;
t80 = t127 * t294 + t128 * t298;
t481 = -qJD(5) * t80 + t294 * t510 + t509 * t298;
t285 = qJD(3) + qJD(5);
t231 = t240 * qJD(2);
t382 = qJD(2) * t299;
t384 = qJD(2) * t295;
t232 = -t287 * t382 - t290 * t384;
t328 = t231 * t294 - t298 * t232;
t119 = t285 * t297 - t293 * t328;
t120 = t285 * t293 + t297 * t328;
t408 = mrSges(6,1) * t285 + mrSges(7,1) * t119 - mrSges(7,2) * t120 - mrSges(6,3) * t328;
t385 = qJD(2) * t289;
t356 = qJD(1) * t385;
t262 = t300 * t356;
t375 = qJDD(1) * t289;
t218 = t296 * t375 + t262;
t199 = qJDD(2) * pkin(8) + t218;
t291 = cos(pkin(6));
t386 = qJD(1) * t291;
t504 = qJD(3) * t386 + t199;
t250 = qJD(2) * pkin(8) + t363;
t345 = qJ(4) * qJD(2) + t250;
t361 = t295 * t386;
t176 = t299 * t345 + t361;
t163 = t287 * t176;
t269 = t299 * t386;
t175 = -t295 * t345 + t269;
t167 = qJD(3) * pkin(3) + t175;
t101 = t290 * t167 - t163;
t440 = pkin(9) * t232;
t81 = qJD(3) * pkin(4) + t101 + t440;
t395 = t290 * t176;
t102 = t287 * t167 + t395;
t441 = pkin(9) * t231;
t87 = t102 + t441;
t44 = -t294 * t87 + t298 * t81;
t341 = mrSges(7,1) * t293 + mrSges(7,2) * t297;
t465 = -m(4) * pkin(8) + m(5) * t292 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - t341;
t442 = pkin(3) * t299;
t275 = pkin(2) + t442;
t212 = -qJD(2) * t275 + qJD(4) - t362;
t156 = -pkin(4) * t231 + t212;
t418 = t285 * Ifges(6,6);
t427 = Ifges(6,4) * t328;
t45 = t294 * t81 + t298 * t87;
t437 = t45 * mrSges(6,3);
t41 = pkin(10) * t285 + t45;
t346 = t298 * t231 + t232 * t294;
t63 = -pkin(5) * t346 - pkin(10) * t328 + t156;
t20 = t293 * t63 + t297 * t41;
t489 = t20 * mrSges(7,2);
t19 = -t293 * t41 + t297 * t63;
t490 = t19 * mrSges(7,1);
t130 = qJD(6) - t346;
t420 = t130 * Ifges(7,3);
t421 = t120 * Ifges(7,5);
t423 = t119 * Ifges(7,6);
t57 = t420 + t421 + t423;
t484 = t346 * Ifges(6,2);
t82 = t418 + t427 + t484;
t503 = -t156 * mrSges(6,1) + t437 + t82 / 0.2e1 - t57 / 0.2e1 + t489 - t490 + t427 / 0.2e1 + t418 / 0.2e1;
t129 = Ifges(6,4) * t346;
t40 = -pkin(5) * t285 - t44;
t316 = t40 * t341;
t118 = Ifges(7,4) * t119;
t59 = t120 * Ifges(7,1) + t130 * Ifges(7,5) + t118;
t412 = t297 * t59;
t419 = t285 * Ifges(6,5);
t438 = t44 * mrSges(6,3);
t445 = t293 / 0.2e1;
t422 = t120 * Ifges(7,4);
t58 = t119 * Ifges(7,2) + t130 * Ifges(7,6) + t422;
t485 = t328 * Ifges(6,1);
t83 = t129 + t419 + t485;
t502 = -t156 * mrSges(6,2) - t316 - t412 / 0.2e1 + t58 * t445 + t438 - t83 / 0.2e1 - t419 / 0.2e1 - t129 / 0.2e1;
t286 = qJ(3) + pkin(12);
t282 = qJ(5) + t286;
t272 = sin(t282);
t273 = cos(t282);
t280 = sin(t286);
t281 = cos(t286);
t344 = -mrSges(4,1) * t299 + mrSges(4,2) * t295;
t459 = m(7) * pkin(10);
t460 = m(7) * pkin(5);
t501 = -m(4) * pkin(2) - m(5) * t275 - t281 * mrSges(5,1) + t280 * mrSges(5,2) - mrSges(3,1) + t344 + (-t460 + t507) * t273 + (-t459 + t508) * t272;
t500 = -m(7) - m(6);
t498 = -t233 / 0.2e1;
t497 = t234 / 0.2e1;
t377 = qJD(2) * qJD(3);
t248 = qJDD(2) * t299 - t295 * t377;
t496 = t248 / 0.2e1;
t200 = -pkin(4) * t240 - t275;
t78 = -pkin(5) * t327 - pkin(10) * t161 + t200;
t39 = t293 * t78 + t297 * t80;
t492 = -qJD(6) * t39 + t293 * t506 + t297 * t505;
t38 = -t293 * t80 + t297 * t78;
t491 = qJD(6) * t38 + t293 * t505 - t297 * t506;
t283 = qJDD(3) + qJDD(5);
t249 = qJDD(2) * t295 + t299 * t377;
t168 = t248 * t290 - t249 * t287;
t169 = t248 * t287 + t249 * t290;
t69 = qJD(5) * t346 + t168 * t294 + t169 * t298;
t49 = qJD(6) * t119 + t283 * t293 + t297 * t69;
t50 = -qJD(6) * t120 + t283 * t297 - t293 * t69;
t16 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t61 = mrSges(6,1) * t283 - mrSges(6,3) * t69;
t488 = t16 - t61;
t104 = -t168 * mrSges(5,1) + t169 * mrSges(5,2);
t70 = -qJD(5) * t328 + t168 * t298 - t169 * t294;
t30 = -t70 * mrSges(6,1) + t69 * mrSges(6,2);
t482 = t104 + t30;
t461 = m(5) * pkin(3);
t478 = -t461 - mrSges(4,1);
t121 = -mrSges(6,2) * t285 + mrSges(6,3) * t346;
t75 = -mrSges(7,2) * t130 + mrSges(7,3) * t119;
t76 = mrSges(7,1) * t130 - mrSges(7,3) * t120;
t333 = -t293 * t76 + t297 * t75;
t477 = -t121 - t333;
t443 = pkin(3) * t290;
t274 = pkin(4) + t443;
t444 = pkin(3) * t287;
t223 = t294 * t274 + t298 * t444;
t398 = t289 * t296;
t194 = -t272 * t398 + t273 * t291;
t195 = t272 * t291 + t273 * t398;
t474 = t507 * t194 + t195 * t508;
t407 = cos(pkin(11));
t348 = t407 * t300;
t288 = sin(pkin(11));
t400 = t288 * t296;
t230 = -t291 * t400 + t348;
t401 = t288 * t289;
t153 = -t230 * t272 + t273 * t401;
t154 = t230 * t273 + t272 * t401;
t473 = t507 * t153 + t154 * t508;
t349 = t407 * t296;
t399 = t288 * t300;
t228 = t291 * t349 + t399;
t350 = t289 * t407;
t151 = -t228 * t272 - t273 * t350;
t152 = t228 * t273 - t272 * t350;
t472 = t507 * t151 + t152 * t508;
t378 = qJD(6) * t297;
t319 = t161 * t378 + t293 * t98;
t374 = qJDD(1) * t291;
t115 = -t250 * t381 + t295 * t374 + t299 * t504;
t193 = t250 * t299 + t361;
t266 = t299 * t374;
t116 = -t193 * qJD(3) - t199 * t295 + t266;
t471 = t115 * t299 - t116 * t295;
t68 = qJDD(6) - t70;
t22 = mrSges(7,1) * t68 - mrSges(7,3) * t49;
t23 = -mrSges(7,2) * t68 + mrSges(7,3) * t50;
t470 = -t293 * t22 + t297 * t23;
t155 = -mrSges(5,1) * t231 - mrSges(5,2) * t232;
t91 = -mrSges(6,1) * t346 + mrSges(6,2) * t328;
t467 = t344 * qJD(2) + t155 + t91;
t376 = qJD(2) * qJD(4);
t380 = qJD(3) * t299;
t94 = -t250 * t380 + qJDD(3) * pkin(3) - qJ(4) * t249 + t266 + (-t376 - t504) * t295;
t95 = qJ(4) * t248 + t299 * t376 + t115;
t53 = -t287 * t95 + t290 * t94;
t35 = qJDD(3) * pkin(4) - pkin(9) * t169 + t53;
t54 = t287 * t94 + t290 * t95;
t37 = pkin(9) * t168 + t54;
t11 = -qJD(5) * t45 - t294 * t37 + t298 * t35;
t261 = t296 * t356;
t217 = t300 * t375 - t261;
t198 = -qJDD(2) * pkin(2) - t217;
t157 = -pkin(3) * t248 + qJDD(4) + t198;
t103 = -pkin(4) * t168 + t157;
t21 = -pkin(5) * t70 - pkin(10) * t69 + t103;
t10 = qJD(5) * t44 + t294 * t35 + t298 * t37;
t5 = pkin(10) * t283 + t10;
t2 = qJD(6) * t19 + t21 * t293 + t297 * t5;
t3 = -qJD(6) * t20 + t21 * t297 - t293 * t5;
t466 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t92 = pkin(5) * t328 - pkin(10) * t346;
t301 = qJD(2) ^ 2;
t458 = t49 / 0.2e1;
t457 = t50 / 0.2e1;
t454 = t68 / 0.2e1;
t451 = -t119 / 0.2e1;
t450 = -t120 / 0.2e1;
t449 = t120 / 0.2e1;
t448 = -t130 / 0.2e1;
t446 = -t232 / 0.2e1;
t433 = mrSges(5,3) * t232;
t432 = mrSges(7,3) * t293;
t431 = mrSges(7,3) * t297;
t430 = Ifges(4,4) * t295;
t429 = Ifges(4,4) * t299;
t428 = Ifges(5,4) * t232;
t426 = Ifges(7,4) * t293;
t425 = Ifges(7,4) * t297;
t424 = t101 * mrSges(5,3);
t404 = t161 * t293;
t403 = t161 * t297;
t397 = t289 * t299;
t396 = t289 * t300;
t108 = t290 * t175 - t163;
t252 = -pkin(3) * t295 - pkin(4) * t280;
t253 = pkin(4) * t281 + t442;
t390 = t230 * t252 + t253 * t401;
t388 = t252 * t398 + t291 * t253;
t383 = qJD(2) * t296;
t379 = qJD(6) * t293;
t373 = Ifges(7,5) * t49 + Ifges(7,6) * t50 + Ifges(7,3) * t68;
t277 = pkin(3) * t384;
t370 = mrSges(4,3) * t384;
t369 = mrSges(4,3) * t382;
t368 = t293 * t396;
t367 = t297 * t396;
t364 = t412 / 0.2e1;
t360 = t289 * t383;
t359 = t300 * t385;
t353 = -t379 / 0.2e1;
t190 = -pkin(4) * t232 + t277;
t352 = t153 * pkin(5) + pkin(10) * t154;
t106 = -t175 * t287 - t395;
t340 = Ifges(7,1) * t297 - t426;
t339 = t299 * Ifges(4,2) + t430;
t338 = -Ifges(7,2) * t293 + t425;
t337 = Ifges(4,5) * t299 - Ifges(4,6) * t295;
t336 = Ifges(7,5) * t297 - Ifges(7,6) * t293;
t335 = t19 * t297 + t20 * t293;
t334 = -t19 * t293 + t20 * t297;
t235 = t291 * t299 - t295 * t398;
t236 = t291 * t295 + t296 * t397;
t133 = t235 * t290 - t236 * t287;
t134 = t235 * t287 + t236 * t290;
t330 = t298 * t133 - t134 * t294;
t85 = t133 * t294 + t134 * t298;
t222 = t274 * t298 - t294 * t444;
t323 = t106 - t441;
t321 = t228 * t252 - t253 * t350;
t73 = -t293 * t85 - t367;
t320 = -t297 * t85 + t368;
t318 = t161 * t379 - t297 * t98;
t315 = t119 * t338;
t314 = t120 * t340;
t313 = t130 * t336;
t251 = -qJD(2) * pkin(2) - t362;
t312 = t251 * (mrSges(4,1) * t295 + mrSges(4,2) * t299);
t311 = t295 * (Ifges(4,1) * t299 - t430);
t227 = -t291 * t348 + t400;
t229 = t291 * t399 + t349;
t306 = -g(1) * t229 - g(2) * t227 + g(3) * t396;
t305 = -qJD(6) * t335 - t3 * t293;
t304 = t2 * t297 + t305;
t14 = t49 * Ifges(7,4) + t50 * Ifges(7,2) + t68 * Ifges(7,6);
t15 = t49 * Ifges(7,1) + t50 * Ifges(7,4) + t68 * Ifges(7,5);
t6 = -pkin(5) * t283 - t11;
t302 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + t15 * t445 + t2 * t431 + t6 * t342 + t297 * t14 / 0.2e1 + Ifges(6,3) * t283 + (Ifges(7,1) * t293 + t425) * t458 + (Ifges(7,2) * t297 + t426) * t457 + t58 * t353 + (Ifges(7,5) * t293 + Ifges(7,6) * t297) * t454 + Ifges(6,6) * t70 + Ifges(6,5) * t69 + (t316 + t364) * qJD(6) + (t315 + t314 + t313) * qJD(6) / 0.2e1;
t284 = -pkin(9) + t292;
t276 = Ifges(4,4) * t382;
t255 = -qJD(3) * mrSges(4,2) + t369;
t254 = qJD(3) * mrSges(4,1) - t370;
t246 = pkin(2) + t253;
t238 = Ifges(4,1) * t384 + Ifges(4,5) * qJD(3) + t276;
t237 = Ifges(4,6) * qJD(3) + qJD(2) * t339;
t221 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t249;
t220 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t248;
t216 = Ifges(5,4) * t231;
t213 = -pkin(5) - t222;
t197 = qJD(3) * mrSges(5,1) + t433;
t196 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t231;
t192 = -t250 * t295 + t269;
t188 = t194 * pkin(5);
t179 = -mrSges(4,1) * t248 + mrSges(4,2) * t249;
t174 = qJD(3) * t235 + t299 * t359;
t173 = -qJD(3) * t236 - t295 * t359;
t150 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t169;
t149 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t168;
t147 = t151 * pkin(5);
t136 = -t232 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t216;
t135 = t231 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t428;
t107 = t173 * t287 + t174 * t290;
t105 = t173 * t290 - t174 * t287;
t90 = t108 + t440;
t71 = t190 + t92;
t62 = -mrSges(6,2) * t283 + mrSges(6,3) * t70;
t52 = t294 * t323 + t298 * t90;
t29 = qJD(5) * t85 - t298 * t105 + t107 * t294;
t28 = qJD(5) * t330 + t105 * t294 + t107 * t298;
t27 = t293 * t92 + t297 * t44;
t26 = -t293 * t44 + t297 * t92;
t25 = t293 * t71 + t297 * t52;
t24 = -t293 * t52 + t297 * t71;
t18 = qJD(6) * t320 - t28 * t293 + t297 * t360;
t17 = qJD(6) * t73 + t28 * t297 + t293 * t360;
t1 = [m(2) * qJDD(1) + t105 * t197 + t107 * t196 + t28 * t121 + t133 * t150 + t134 * t149 + t17 * t75 + t173 * t254 + t174 * t255 + t18 * t76 + t73 * t22 + t236 * t220 + t235 * t221 - t320 * t23 + t85 * t62 - t488 * t330 - t408 * t29 + (-m(2) - m(3) - m(4) - m(5) + t500) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t301 - t179 - t482) * t300 + (-mrSges(3,1) * t301 - mrSges(3,2) * qJDD(2) + qJD(2) * t467) * t296) * t289 + m(5) * (t101 * t105 + t102 * t107 + t133 * t53 + t134 * t54 + (-t157 * t300 + t212 * t383) * t289) + m(4) * (t115 * t236 + t116 * t235 + t173 * t192 + t174 * t193 + (-t198 * t300 + t251 * t383) * t289) + m(3) * (qJDD(1) * t291 ^ 2 + (t217 * t300 + t218 * t296) * t289) + m(6) * (t10 * t85 + t11 * t330 + t28 * t45 - t29 * t44 + (-t103 * t300 + t156 * t383) * t289) + m(7) * (t17 * t20 + t18 * t19 - t2 * t320 + t29 * t40 + t3 * t73 - t330 * t6); t475 * t196 + (t101 * t476 + t102 * t475 - t157 * t275 + t183 * t53 + t184 * t54 + t212 * t278) * m(5) + t476 * t197 + (t500 * t246 * t396 + (t501 * t300 + (-t284 * t500 + t465) * t296) * t289) * g(3) + (t262 - t218) * mrSges(3,2) + t238 * t380 / 0.2e1 - t237 * t381 / 0.2e1 + (t500 * (-t229 * t246 - t230 * t284) + t465 * t230 - t501 * t229) * g(1) + (t500 * (-t227 * t246 - t228 * t284) + t465 * t228 - t501 * t227) * g(2) + t249 * t429 / 0.2e1 + (t10 * t80 + t103 * t200 + t11 * t331 + t156 * t191 + t44 * t481 - t45 * t506) * m(6) - t506 * t121 + (-pkin(2) * t198 - (t251 * t296 + (-t192 * t295 + t193 * t299) * t300) * t387) * m(4) + (m(4) * ((-t192 * t299 - t193 * t295) * qJD(3) + t471) - t254 * t380 - t255 * t381 + t299 * t220 - t295 * t221) * pkin(8) + (-t192 * t380 - t193 * t381 + t471) * mrSges(4,3) - t319 * t58 / 0.2e1 + (-t102 * t233 + t240 * t54 - t241 * t53) * mrSges(5,3) + t231 * (Ifges(5,4) * t234 - Ifges(5,2) * t233) / 0.2e1 + t212 * (mrSges(5,1) * t233 + mrSges(5,2) * t234) + (Ifges(5,1) * t234 - Ifges(5,4) * t233) * t446 + (t103 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t69 + Ifges(6,4) * t70 + Ifges(6,5) * t283 + t336 * t454 + t338 * t457 + t340 * t458 + t341 * t6 + t353 * t59) * t161 + (Ifges(4,5) * t295 + Ifges(5,5) * t241 + Ifges(4,6) * t299 + Ifges(5,6) * t240) * qJDD(3) + t40 * (mrSges(7,1) * t319 - mrSges(7,2) * t318) + t130 * (-Ifges(7,5) * t318 - Ifges(7,6) * t319 + Ifges(7,3) * t99) / 0.2e1 + t119 * (-Ifges(7,4) * t318 - Ifges(7,2) * t319 + Ifges(7,6) * t99) / 0.2e1 + t299 * (Ifges(4,4) * t249 + Ifges(4,2) * t248) / 0.2e1 - t98 * t438 + (t299 * (-Ifges(4,2) * t295 + t429) + t311) * t377 / 0.2e1 + t346 * (Ifges(6,4) * t98 - Ifges(6,2) * t99) / 0.2e1 + t328 * (Ifges(6,1) * t98 - Ifges(6,4) * t99) / 0.2e1 - t99 * t437 + t15 * t403 / 0.2e1 - t14 * t404 / 0.2e1 - t234 * t424 - t299 * t255 * t362 + (t261 + t217) * mrSges(3,1) + (-m(5) * t212 - m(6) * t156 - t467) * t363 + t136 * t497 + t135 * t498 + t99 * t490 + t339 * t496 + (t19 * t318 - t2 * t404 - t20 * t319 - t3 * t403) * mrSges(7,3) + t198 * t344 + Ifges(3,3) * qJDD(2) + t285 * (Ifges(6,5) * t98 - Ifges(6,6) * t99) / 0.2e1 - t275 * t104 + t169 * (Ifges(5,1) * t241 + Ifges(5,4) * t240) + t168 * (Ifges(5,4) * t241 + Ifges(5,2) * t240) + t157 * (-mrSges(5,1) * t240 + mrSges(5,2) * t241) + t200 * t30 + t191 * t91 - pkin(2) * t179 + t183 * t150 + t184 * t149 + t156 * (mrSges(6,1) * t99 + mrSges(6,2) * t98) + t98 * t83 / 0.2e1 + t99 * t57 / 0.2e1 - t99 * t82 / 0.2e1 + t80 * t62 + t38 * t22 + t39 * t23 + t481 * t408 + t155 * t278 - t99 * t489 + t491 * t75 + t492 * t76 + (t312 + Ifges(5,5) * t497 + Ifges(5,6) * t498 + t337 * qJD(3) / 0.2e1) * qJD(3) + t98 * t364 + (-Ifges(7,1) * t318 - Ifges(7,4) * t319 + Ifges(7,5) * t99) * t449 + (Ifges(4,1) * t249 + Ifges(4,4) * t496 + t254 * t362) * t295 - t488 * t331 + (t19 * t492 + t2 * t39 + t20 * t491 + t3 * t38 - t331 * t6 - t40 * t481) * m(7) - (t373 / 0.2e1 - Ifges(6,4) * t69 - Ifges(6,2) * t70 - Ifges(6,6) * t283 + t103 * mrSges(6,1) - t10 * mrSges(6,3) + Ifges(7,3) * t454 + Ifges(7,6) * t457 + Ifges(7,5) * t458 + t466) * t327; -t3 * t432 - t102 * t433 + t237 * t384 / 0.2e1 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) - t301 * t311 / 0.2e1 + (m(7) * t304 - t378 * t76 - t379 * t75 + t470) * (pkin(10) + t223) - qJD(2) * t312 + (-t255 + t369) * t192 + ((m(6) * t45 + m(7) * t334 - t477) * qJD(5) + t61 + m(6) * t11) * t222 + (t10 * t223 - t156 * t190 - t45 * t52) * m(6) + (t20 * t432 + t19 * t431 - t485 / 0.2e1 + t336 * t448 + t340 * t450 + t338 * t451 + t502) * t346 + (t484 / 0.2e1 + Ifges(7,3) * t448 + Ifges(7,5) * t450 + Ifges(7,6) * t451 + t503) * t328 + (-(-t228 * t280 - t281 * t350) * mrSges(5,1) - (-t228 * t281 + t280 * t350) * mrSges(5,2) - (-t228 * t299 + t295 * t350) * mrSges(4,2) - m(7) * (t152 * pkin(10) + t147 + t321) - m(6) * t321 + t478 * (-t228 * t295 - t299 * t350) + t472) * g(2) + (-m(7) * (t352 + t390) - (-t230 * t280 + t281 * t401) * mrSges(5,1) - (-t230 * t281 - t280 * t401) * mrSges(5,2) - (-t230 * t299 - t295 * t401) * mrSges(4,2) - m(6) * t390 + t478 * (-t230 * t295 + t288 * t397) + t473) * g(1) + (-m(7) * (pkin(10) * t195 + t188 + t388) - m(6) * t388 - (-t280 * t398 + t281 * t291) * mrSges(5,1) - (-t280 * t291 - t281 * t398) * mrSges(5,2) + mrSges(4,2) * t236 + t478 * t235 + t474) * g(3) + (-t19 * t378 - t20 * t379) * mrSges(7,3) + t302 + t232 * (Ifges(5,1) * t231 + t428) / 0.2e1 - (-Ifges(4,2) * t384 + t238 + t276) * t382 / 0.2e1 - (Ifges(5,2) * t232 + t136 + t216) * t231 / 0.2e1 - t155 * t277 - m(5) * (t101 * t106 + t102 * t108 + t212 * t277) + t231 * t424 - t337 * t377 / 0.2e1 + (t254 + t370) * t193 + (-t19 * t24 - t20 * t25 + t213 * t6) * m(7) + Ifges(4,6) * t248 + Ifges(4,5) * t249 - qJD(3) * (Ifges(5,5) * t231 + Ifges(5,6) * t232) / 0.2e1 - t212 * (-mrSges(5,1) * t232 + mrSges(5,2) * t231) + t223 * t62 + t213 * t16 - t190 * t91 - t108 * t196 - t106 * t197 + Ifges(5,5) * t169 + Ifges(5,6) * t168 - t52 * t121 - t115 * mrSges(4,2) + t116 * mrSges(4,1) - t25 * t75 - t24 * t76 + t53 * mrSges(5,1) - t54 * mrSges(5,2) + t150 * t443 + t149 * t444 + (-m(6) * t44 + m(7) * t40 - t408) * (t223 * qJD(5) - t294 * t90 + t298 * t323) + t135 * t446 + (t287 * t54 + t290 * t53) * t461; -t231 * t196 - t232 * t197 + t297 * t22 + t293 * t23 + t408 * t328 + t333 * qJD(6) + t477 * t346 + (t130 * t334 + t2 * t293 + t3 * t297 - t328 * t40 + t306) * m(7) + (t328 * t44 - t346 * t45 + t103 + t306) * m(6) + (-t101 * t232 - t102 * t231 + t157 + t306) * m(5) + t482; (-m(7) * t147 + t472) * g(2) + (-m(7) * t352 + t473) * g(1) + (-t313 / 0.2e1 - t315 / 0.2e1 - t314 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t328 + t335 * mrSges(7,3) + t502) * t346 + t305 * mrSges(7,3) + ((-t293 * t75 - t297 * t76) * qJD(6) + (-g(2) * t152 - g(3) * t195) * m(7) + t470) * pkin(10) + t304 * t459 - t6 * t460 + t408 * t45 + (-t420 / 0.2e1 - t423 / 0.2e1 - t421 / 0.2e1 + t503) * t328 - m(7) * (t19 * t26 + t20 * t27 + t40 * t45) + t302 + (-m(7) * t188 + t474) * g(3) - t44 * t121 - t27 * t75 - t26 * t76 - pkin(5) * t16; -t40 * (mrSges(7,1) * t120 + mrSges(7,2) * t119) + (Ifges(7,1) * t119 - t422) * t450 + t58 * t449 + (Ifges(7,5) * t119 - Ifges(7,6) * t120) * t448 - t19 * t75 + t20 * t76 - g(1) * ((-t154 * t293 + t229 * t297) * mrSges(7,1) + (-t154 * t297 - t229 * t293) * mrSges(7,2)) - g(2) * ((-t152 * t293 + t227 * t297) * mrSges(7,1) + (-t152 * t297 - t227 * t293) * mrSges(7,2)) - g(3) * ((-t195 * t293 - t367) * mrSges(7,1) + (-t195 * t297 + t368) * mrSges(7,2)) + (t119 * t19 + t120 * t20) * mrSges(7,3) + t373 + (-Ifges(7,2) * t120 + t118 + t59) * t451 + t466;];
tau  = t1;
