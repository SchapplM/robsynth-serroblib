% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR1
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
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:54:11
% EndTime: 2019-03-09 06:54:35
% DurationCPUTime: 14.18s
% Computational Cost: add. (16245->678), mult. (34981->904), div. (0->0), fcn. (24922->18), ass. (0->321)
t274 = qJ(3) + qJ(4);
t267 = qJ(5) + t274;
t248 = sin(t267);
t249 = cos(t267);
t277 = sin(qJ(6));
t415 = mrSges(7,2) * t277;
t488 = t248 * t415 + t249 * (m(7) * pkin(10) + mrSges(7,3));
t271 = qJD(3) + qJD(4);
t261 = qJD(5) + t271;
t282 = cos(qJ(6));
t279 = sin(qJ(4));
t280 = sin(qJ(3));
t284 = cos(qJ(4));
t285 = cos(qJ(3));
t217 = -t279 * t280 + t284 * t285;
t209 = t217 * qJD(1);
t218 = t279 * t285 + t280 * t284;
t210 = t218 * qJD(1);
t278 = sin(qJ(5));
t283 = cos(qJ(5));
t315 = t209 * t278 + t283 * t210;
t133 = t261 * t282 - t277 * t315;
t134 = t261 * t277 + t282 * t315;
t391 = mrSges(6,1) * t261 + mrSges(7,1) * t133 - mrSges(7,2) * t134 - mrSges(6,3) * t315;
t275 = sin(pkin(11));
t244 = pkin(1) * t275 + pkin(7);
t231 = t244 * qJD(1);
t338 = pkin(8) * qJD(1) + t231;
t370 = qJD(2) * t280;
t175 = t285 * t338 + t370;
t169 = t279 * t175;
t264 = t285 * qJD(2);
t174 = -t280 * t338 + t264;
t172 = qJD(3) * pkin(3) + t174;
t117 = t284 * t172 - t169;
t203 = t210 * pkin(9);
t109 = t117 - t203;
t102 = pkin(4) * t271 + t109;
t171 = t284 * t175;
t118 = t172 * t279 + t171;
t427 = pkin(9) * t209;
t110 = t118 + t427;
t389 = t110 * t278;
t56 = t102 * t283 - t389;
t53 = -pkin(5) * t261 - t56;
t487 = -m(6) * t56 + m(7) * t53 - t391;
t337 = t283 * t209 - t210 * t278;
t146 = Ifges(6,4) * t337;
t276 = cos(pkin(11));
t245 = -pkin(1) * t276 - pkin(2);
t268 = t285 * pkin(3);
t227 = t245 - t268;
t211 = t227 * qJD(1);
t164 = -pkin(4) * t209 + t211;
t328 = mrSges(7,1) * t277 + mrSges(7,2) * t282;
t309 = t53 * t328;
t132 = Ifges(7,4) * t133;
t147 = qJD(6) - t337;
t70 = t134 * Ifges(7,1) + t147 * Ifges(7,5) + t132;
t393 = t282 * t70;
t399 = t261 * Ifges(6,5);
t422 = t56 * mrSges(6,3);
t436 = t277 / 0.2e1;
t404 = t134 * Ifges(7,4);
t69 = t133 * Ifges(7,2) + t147 * Ifges(7,6) + t404;
t95 = Ifges(6,1) * t315 + t146 + t399;
t486 = -t393 / 0.2e1 + t69 * t436 + t422 - t95 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t315 - t146 / 0.2e1 - t309 - t164 * mrSges(6,2) - t399 / 0.2e1;
t398 = t261 * Ifges(6,6);
t410 = Ifges(6,4) * t315;
t375 = t283 * t110;
t57 = t102 * t278 + t375;
t421 = t57 * mrSges(6,3);
t54 = pkin(10) * t261 + t57;
t81 = -pkin(5) * t337 - pkin(10) * t315 + t164;
t23 = t277 * t81 + t282 * t54;
t468 = t23 * mrSges(7,2);
t22 = -t277 * t54 + t282 * t81;
t469 = t22 * mrSges(7,1);
t402 = t147 * Ifges(7,3);
t403 = t134 * Ifges(7,5);
t405 = t133 * Ifges(7,6);
t68 = t402 + t403 + t405;
t94 = Ifges(6,2) * t337 + t398 + t410;
t485 = t468 + t421 + t94 / 0.2e1 - t68 / 0.2e1 + t410 / 0.2e1 - t164 * mrSges(6,1) + t398 / 0.2e1 - t469;
t464 = t249 * pkin(5) + t248 * pkin(10);
t482 = m(7) * t464;
t481 = -t249 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t248;
t270 = qJDD(3) + qJDD(4);
t258 = qJDD(5) + t270;
t361 = qJD(1) * qJD(3);
t224 = qJDD(1) * t285 - t280 * t361;
t225 = qJDD(1) * t280 + t285 * t361;
t302 = t217 * qJD(4);
t139 = qJD(1) * t302 + t224 * t279 + t225 * t284;
t303 = t218 * qJD(4);
t140 = -qJD(1) * t303 + t224 * t284 - t225 * t279;
t74 = qJD(5) * t337 + t139 * t283 + t140 * t278;
t39 = qJD(6) * t133 + t258 * t277 + t282 * t74;
t75 = -qJD(5) * t315 - t139 * t278 + t140 * t283;
t73 = qJDD(6) - t75;
t18 = mrSges(7,1) * t73 - mrSges(7,3) * t39;
t40 = -qJD(6) * t134 + t258 * t282 - t277 * t74;
t19 = -mrSges(7,2) * t73 + mrSges(7,3) * t40;
t321 = -t277 * t18 + t282 * t19;
t362 = qJD(6) * t282;
t363 = qJD(6) * t277;
t90 = -mrSges(7,2) * t147 + mrSges(7,3) * t133;
t91 = mrSges(7,1) * t147 - mrSges(7,3) * t134;
t480 = -t91 * t362 - t90 * t363 + t321;
t265 = sin(t274);
t266 = cos(t274);
t416 = mrSges(6,2) * t249;
t479 = mrSges(5,1) * t265 + mrSges(6,1) * t248 + mrSges(5,2) * t266 + t416;
t418 = mrSges(7,1) * t282;
t478 = t415 - t418;
t272 = qJ(1) + pkin(11);
t259 = sin(t272);
t260 = cos(t272);
t477 = g(1) * t260 + g(2) * t259;
t471 = -m(7) - m(6);
t470 = t224 / 0.2e1;
t16 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t61 = mrSges(6,1) * t258 - mrSges(6,3) * t74;
t467 = t16 - t61;
t419 = pkin(8) + t244;
t214 = t419 * t280;
t215 = t419 * t285;
t156 = -t279 * t214 + t284 * t215;
t340 = t266 * mrSges(5,1) - mrSges(5,2) * t265;
t463 = t249 * t478 + t481;
t163 = t217 * t278 + t218 * t283;
t167 = qJD(3) * t217 + t302;
t168 = -qJD(3) * t218 - t303;
t314 = t283 * t217 - t218 * t278;
t87 = qJD(5) * t314 + t167 * t283 + t168 * t278;
t311 = t163 * t362 + t277 * t87;
t461 = -t22 * t362 - t23 * t363;
t229 = t244 * qJDD(1);
t369 = qJD(3) * t280;
t148 = qJD(3) * t264 + t280 * qJDD(2) + t285 * t229 - t231 * t369;
t191 = t231 * t285 + t370;
t149 = -qJD(3) * t191 + t285 * qJDD(2) - t229 * t280;
t460 = t148 * t285 - t149 * t280;
t458 = -m(5) - m(4) - m(3);
t141 = -mrSges(6,2) * t261 + mrSges(6,3) * t337;
t457 = -t277 * t91 + t282 * t90 + t141;
t455 = -t340 + t463;
t124 = qJDD(3) * pkin(3) - pkin(8) * t225 + t149;
t129 = pkin(8) * t224 + t148;
t52 = -qJD(4) * t118 + t284 * t124 - t129 * t279;
t32 = pkin(4) * t270 - pkin(9) * t139 + t52;
t366 = qJD(4) * t284;
t367 = qJD(4) * t279;
t51 = t279 * t124 + t284 * t129 + t172 * t366 - t175 * t367;
t34 = pkin(9) * t140 + t51;
t11 = -qJD(5) * t57 - t278 * t34 + t283 * t32;
t230 = t245 * qJDD(1);
t176 = -pkin(3) * t224 + t230;
t111 = -pkin(4) * t140 + t176;
t17 = -pkin(5) * t75 - pkin(10) * t74 + t111;
t364 = qJD(5) * t283;
t365 = qJD(5) * t278;
t10 = t102 * t364 - t110 * t365 + t278 * t32 + t283 * t34;
t7 = pkin(10) * t258 + t10;
t2 = qJD(6) * t22 + t17 * t277 + t282 * t7;
t3 = -qJD(6) * t23 + t17 * t282 - t277 * t7;
t454 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t103 = pkin(5) * t315 - pkin(10) * t337;
t235 = -mrSges(4,1) * t285 + mrSges(4,2) * t280;
t453 = mrSges(3,1) + m(5) * (t268 + pkin(2)) + t340 + m(4) * pkin(2) - t235 - t481;
t320 = t22 * t282 + t23 * t277;
t423 = t277 * t3;
t295 = -qJD(6) * t320 - t423;
t424 = t2 * t282;
t292 = t295 + t424;
t452 = m(7) * t292 + t480;
t287 = -pkin(8) - pkin(7);
t451 = -m(4) * pkin(7) + m(5) * t287 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t450 = m(7) * pkin(5);
t449 = t39 / 0.2e1;
t448 = t40 / 0.2e1;
t445 = t73 / 0.2e1;
t442 = -t133 / 0.2e1;
t441 = -t134 / 0.2e1;
t440 = t134 / 0.2e1;
t439 = -t147 / 0.2e1;
t437 = t210 / 0.2e1;
t281 = sin(qJ(1));
t435 = pkin(1) * t281;
t434 = pkin(3) * t280;
t433 = pkin(3) * t284;
t432 = pkin(4) * t210;
t431 = pkin(4) * t265;
t247 = pkin(4) * t266;
t430 = pkin(4) * t278;
t429 = pkin(4) * t283;
t428 = pkin(5) * t248;
t286 = cos(qJ(1));
t269 = t286 * pkin(1);
t414 = mrSges(7,3) * t277;
t413 = mrSges(7,3) * t282;
t412 = Ifges(4,4) * t280;
t411 = Ifges(4,4) * t285;
t409 = Ifges(7,4) * t277;
t408 = Ifges(7,4) * t282;
t407 = t117 * mrSges(5,3);
t406 = t118 * mrSges(5,3);
t401 = t210 * Ifges(5,4);
t386 = t163 * t277;
t385 = t163 * t282;
t381 = t259 * t277;
t380 = t259 * t282;
t379 = t260 * t277;
t378 = t260 * t282;
t377 = t278 * t279;
t376 = t279 * t283;
t126 = t284 * t174 - t169;
t252 = pkin(4) + t433;
t205 = pkin(3) * t376 + t278 * t252;
t373 = t247 + t268;
t372 = qJD(1) * t280;
t371 = qJD(1) * t285;
t368 = qJD(3) * t285;
t360 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t73;
t256 = pkin(3) * t369;
t356 = t248 * t418;
t255 = pkin(3) * t372;
t355 = mrSges(4,3) * t372;
t354 = mrSges(4,3) * t371;
t347 = t393 / 0.2e1;
t346 = t247 + t464;
t343 = -t363 / 0.2e1;
t154 = -pkin(4) * t168 + t256;
t342 = qJD(3) * t419;
t125 = -t174 * t279 - t171;
t155 = -t284 * t214 - t215 * t279;
t336 = t488 * t259;
t335 = t488 * t260;
t332 = mrSges(4,1) * t280 + mrSges(4,2) * t285;
t327 = Ifges(7,1) * t282 - t409;
t326 = t285 * Ifges(4,2) + t412;
t325 = -Ifges(7,2) * t277 + t408;
t324 = Ifges(4,5) * t285 - Ifges(4,6) * t280;
t323 = Ifges(7,5) * t282 - Ifges(7,6) * t277;
t319 = -t22 * t277 + t23 * t282;
t135 = -pkin(9) * t218 + t155;
t136 = pkin(9) * t217 + t156;
t84 = t135 * t278 + t136 * t283;
t173 = -pkin(4) * t217 + t227;
t92 = -pkin(5) * t314 - pkin(10) * t163 + t173;
t42 = t277 * t92 + t282 * t84;
t41 = -t277 * t84 + t282 * t92;
t317 = t283 * t135 - t136 * t278;
t204 = -pkin(3) * t377 + t252 * t283;
t312 = t125 - t427;
t310 = t163 * t363 - t282 * t87;
t308 = t133 * t325;
t307 = t134 * t327;
t306 = t147 * t323;
t305 = t245 * qJD(1) * t332;
t304 = t280 * (Ifges(4,1) * t285 - t412);
t200 = t280 * t342;
t201 = t285 * t342;
t104 = -t284 * t200 - t279 * t201 - t214 * t366 - t215 * t367;
t89 = t103 + t432;
t226 = -t431 - t434;
t296 = m(7) * (t226 - t428) - t356;
t294 = m(7) * (-t428 - t431) - t356;
t293 = t416 + (mrSges(6,1) + t418 + t450) * t248;
t105 = -qJD(4) * t156 + t200 * t279 - t284 * t201;
t291 = -pkin(9) * t167 + t105;
t14 = t39 * Ifges(7,4) + t40 * Ifges(7,2) + t73 * Ifges(7,6);
t15 = t39 * Ifges(7,1) + t40 * Ifges(7,4) + t73 * Ifges(7,5);
t8 = -pkin(5) * t258 - t11;
t290 = -t10 * mrSges(6,2) + t2 * t413 + t15 * t436 + t282 * t14 / 0.2e1 + Ifges(6,3) * t258 + (Ifges(7,1) * t277 + t408) * t449 + (Ifges(7,2) * t282 + t409) * t448 + t8 * t478 + (Ifges(7,5) * t277 + Ifges(7,6) * t282) * t445 + t69 * t343 + Ifges(6,6) * t75 + Ifges(6,5) * t74 + t11 * mrSges(6,1) + (t309 + t347) * qJD(6) + (t308 + t307 + t306) * qJD(6) / 0.2e1;
t144 = t209 * Ifges(5,2) + t271 * Ifges(5,6) + t401;
t202 = Ifges(5,4) * t209;
t145 = t210 * Ifges(5,1) + t271 * Ifges(5,5) + t202;
t289 = t461 * mrSges(7,3) - t211 * (mrSges(5,1) * t210 + mrSges(5,2) * t209) + t144 * t437 + t290 + t210 * t406 + t209 * t407 - t3 * t414 - t210 * (Ifges(5,1) * t209 - t401) / 0.2e1 + Ifges(5,3) * t270 - t271 * (Ifges(5,5) * t209 - Ifges(5,6) * t210) / 0.2e1 + Ifges(5,6) * t140 + Ifges(5,5) * t139 - t51 * mrSges(5,2) + t52 * mrSges(5,1) - (-Ifges(5,2) * t210 + t145 + t202) * t209 / 0.2e1 + (Ifges(7,5) * t441 + Ifges(7,6) * t442 + Ifges(7,3) * t439 + t485) * t315 + (t22 * t413 + t23 * t414 + t323 * t439 + t325 * t442 + t327 * t441 + t486) * t337;
t273 = -pkin(9) + t287;
t254 = Ifges(4,4) * t371;
t251 = -pkin(5) - t429;
t234 = -qJD(3) * mrSges(4,2) + t354;
t232 = qJD(3) * mrSges(4,1) - t355;
t223 = pkin(2) + t373;
t208 = Ifges(4,1) * t372 + Ifges(4,5) * qJD(3) + t254;
t207 = Ifges(4,6) * qJD(3) + qJD(1) * t326;
t198 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t225;
t197 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t224;
t195 = -pkin(5) - t204;
t190 = -t231 * t280 + t264;
t183 = t249 * t378 + t381;
t182 = -t249 * t379 + t380;
t181 = -t249 * t380 + t379;
t180 = t249 * t381 + t378;
t179 = mrSges(5,1) * t271 - mrSges(5,3) * t210;
t178 = -mrSges(5,2) * t271 + mrSges(5,3) * t209;
t177 = t255 + t432;
t160 = -mrSges(5,1) * t209 + mrSges(5,2) * t210;
t122 = -mrSges(5,2) * t270 + mrSges(5,3) * t140;
t121 = mrSges(5,1) * t270 - mrSges(5,3) * t139;
t112 = -t203 + t126;
t101 = -mrSges(6,1) * t337 + mrSges(6,2) * t315;
t88 = qJD(5) * t163 + t167 * t278 - t283 * t168;
t86 = t255 + t89;
t82 = pkin(9) * t168 + t104;
t66 = t283 * t112 + t278 * t312;
t62 = -mrSges(6,2) * t258 + mrSges(6,3) * t75;
t59 = t109 * t283 - t389;
t58 = t109 * t278 + t375;
t30 = pkin(5) * t88 - pkin(10) * t87 + t154;
t29 = t103 * t277 + t282 * t56;
t28 = t103 * t282 - t277 * t56;
t27 = t277 * t86 + t282 * t66;
t26 = -t277 * t66 + t282 * t86;
t25 = t277 * t89 + t282 * t59;
t24 = -t277 * t59 + t282 * t89;
t20 = qJD(5) * t317 + t278 * t291 + t283 * t82;
t5 = -qJD(6) * t42 - t20 * t277 + t282 * t30;
t4 = qJD(6) * t41 + t20 * t282 + t277 * t30;
t1 = [(t285 * t197 - t280 * t198 + m(4) * ((-t190 * t285 - t191 * t280) * qJD(3) + t460) - t232 * t368 - t234 * t369) * t244 + (-t190 * t368 - t191 * t369 + t460) * mrSges(4,3) - t311 * t69 / 0.2e1 - (-m(6) * t11 + m(7) * t8 + t467) * t317 + (m(4) * t245 + t235) * t230 + (mrSges(5,2) * t176 - mrSges(5,3) * t52 + Ifges(5,1) * t139 + Ifges(5,4) * t140 + Ifges(5,5) * t270) * t218 + m(7) * (t2 * t42 + t22 * t5 + t23 * t4 + t3 * t41) + m(6) * (t10 * t84 + t111 * t173 + t154 * t164 + t20 * t57) + (-mrSges(5,1) * t176 + mrSges(5,3) * t51 + Ifges(5,4) * t139 + Ifges(5,2) * t140 + Ifges(5,6) * t270) * t217 + (-mrSges(2,1) * t286 - t183 * mrSges(7,1) + mrSges(2,2) * t281 - t182 * mrSges(7,2) + t471 * (t260 * t223 - t259 * t273 + t269) + t458 * t269 + t451 * t259 + (-t453 - t482) * t260) * g(2) + (t111 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t74 + Ifges(6,4) * t75 + Ifges(6,5) * t258 + t323 * t445 + t325 * t448 + t327 * t449 + t328 * t8 + t343 * t70) * t163 + t87 * t347 + (t305 + t324 * qJD(3) / 0.2e1) * qJD(3) + (t304 + t285 * (-Ifges(4,2) * t280 + t411)) * t361 / 0.2e1 + t337 * (Ifges(6,4) * t87 - Ifges(6,2) * t88) / 0.2e1 + t315 * (Ifges(6,1) * t87 - Ifges(6,4) * t88) / 0.2e1 + t160 * t256 + (Ifges(5,1) * t167 + Ifges(5,4) * t168) * t437 + (-Ifges(7,1) * t310 - Ifges(7,4) * t311 + Ifges(7,5) * t88) * t440 + t225 * t411 / 0.2e1 + (-t2 * t386 + t22 * t310 - t23 * t311 - t3 * t385) * mrSges(7,3) + t487 * (qJD(5) * t84 + t278 * t82 - t283 * t291) - (Ifges(7,3) * t445 + Ifges(7,6) * t448 + Ifges(7,5) * t449 - Ifges(6,4) * t74 - Ifges(6,2) * t75 - Ifges(6,6) * t258 + t360 / 0.2e1 + t111 * mrSges(6,1) - t10 * mrSges(6,3) + t454) * t314 + (mrSges(2,1) * t281 - t181 * mrSges(7,1) + mrSges(2,2) * t286 - t180 * mrSges(7,2) + t471 * (-t260 * t273 - t435) - t458 * t435 + t451 * t260 + (-m(7) * (-t223 - t464) + m(6) * t223 + t453) * t259) * g(1) + t285 * (Ifges(4,4) * t225 + Ifges(4,2) * t224) / 0.2e1 + t88 * t469 + t326 * t470 + t168 * t406 - t88 * t421 - t87 * t422 - t167 * t407 + qJDD(3) * (Ifges(4,5) * t280 + Ifges(4,6) * t285) + t271 * (Ifges(5,5) * t167 + Ifges(5,6) * t168) / 0.2e1 + t261 * (Ifges(6,5) * t87 - Ifges(6,6) * t88) / 0.2e1 + t15 * t385 / 0.2e1 - t14 * t386 / 0.2e1 + t245 * (-mrSges(4,1) * t224 + mrSges(4,2) * t225) + t227 * (-mrSges(5,1) * t140 + mrSges(5,2) * t139) + t211 * (-mrSges(5,1) * t168 + mrSges(5,2) * t167) - t207 * t369 / 0.2e1 + t208 * t368 / 0.2e1 + t209 * (Ifges(5,4) * t167 + Ifges(5,2) * t168) / 0.2e1 + t173 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) + m(5) * (t104 * t118 + t105 * t117 + t155 * t52 + t156 * t51 + t176 * t227 + t211 * t256) + t104 * t178 + t105 * t179 + t168 * t144 / 0.2e1 + t164 * (mrSges(6,1) * t88 + mrSges(6,2) * t87) + t167 * t145 / 0.2e1 + t154 * t101 + t155 * t121 + t156 * t122 + t20 * t141 - t88 * t94 / 0.2e1 + t87 * t95 / 0.2e1 + t88 * t68 / 0.2e1 + t4 * t90 + t5 * t91 + t84 * t62 + t42 * t19 + t41 * t18 - t88 * t468 + (Ifges(4,1) * t225 + Ifges(4,4) * t470) * t280 + t53 * (mrSges(7,1) * t311 - mrSges(7,2) * t310) + t133 * (-Ifges(7,4) * t310 - Ifges(7,2) * t311 + Ifges(7,6) * t88) / 0.2e1 + t147 * (-Ifges(7,5) * t310 - Ifges(7,6) * t311 + Ifges(7,3) * t88) / 0.2e1 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t276 - 0.2e1 * mrSges(3,2) * t275 + m(3) * (t275 ^ 2 + t276 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1); m(3) * qJDD(2) + t217 * t121 + t218 * t122 + t167 * t178 + t168 * t179 + t280 * t197 + t285 * t198 - t391 * t88 - t467 * t314 + (-t232 * t280 + t234 * t285) * qJD(3) + t457 * t87 + (t62 + (-t277 * t90 - t282 * t91) * qJD(6) + t321) * t163 + (t458 + t471) * g(3) + m(4) * (t148 * t280 + t149 * t285 + (-t190 * t280 + t191 * t285) * qJD(3)) + m(7) * (t163 * t292 - t314 * t8 + t319 * t87 + t53 * t88) + m(6) * (t10 * t163 + t11 * t314 - t56 * t88 + t57 * t87) + m(5) * (t117 * t168 + t118 * t167 + t217 * t52 + t218 * t51); t487 * (-t112 * t278 + t283 * t312 + t252 * t365 + (t279 * t364 + (t278 * t284 + t376) * qJD(4)) * pkin(3)) + t477 * (m(5) * t434 - m(6) * t226 + t332 + t479) + (t195 * t8 - t22 * t26 - t23 * t27) * m(7) + (t10 * t205 + t11 * t204 - t164 * t177 - t57 * t66) * m(6) - (-Ifges(4,2) * t372 + t208 + t254) * t371 / 0.2e1 + (t355 + t232) * t191 + t452 * (pkin(10) + t205) + t121 * t433 + t289 + (-m(5) * t268 + t235 - m(6) * t373 - m(7) * (t268 + t346) + t455) * g(3) + (m(6) * t57 + m(7) * t319 + t457) * (t252 * t364 + (-t279 * t365 + (t283 * t284 - t377) * qJD(4)) * pkin(3)) + (t354 - t234) * t190 + Ifges(4,6) * t224 + Ifges(4,5) * t225 + t207 * t372 / 0.2e1 - t324 * t361 / 0.2e1 + t204 * t61 + t205 * t62 + t195 * t16 - t177 * t101 - t126 * t178 - t125 * t179 - m(5) * (t117 * t125 + t118 * t126 + t211 * t255) - t160 * t255 - t148 * mrSges(4,2) + t149 * mrSges(4,1) - t66 * t141 + Ifges(4,3) * qJDD(3) - t27 * t90 - t26 * t91 + (m(5) * (t279 * t51 + t284 * t52 + (-t117 * t279 + t118 * t284) * qJD(4)) + t279 * t122 - t179 * t367 + t178 * t366) * pkin(3) - g(1) * (t260 * t296 + t335) - g(2) * (t259 * t296 + t336) + (-t305 - t304 * qJD(1) / 0.2e1) * qJD(1); t61 * t429 + t62 * t430 + t289 - t101 * t432 + t251 * t16 - t117 * t178 + t118 * t179 - t59 * t141 - t25 * t90 - t24 * t91 - g(1) * (t260 * t294 + t335) - g(2) * (t259 * t294 + t336) + (t251 * t8 + (t278 * t53 + t283 * t319) * qJD(5) * pkin(4) - t22 * t24 - t23 * t25 - t53 * t58) * m(7) + ((t10 * t278 + t11 * t283 + (-t278 * t56 + t283 * t57) * qJD(5)) * pkin(4) - t164 * t432 + t56 * t58 - t57 * t59) * m(6) + t457 * pkin(4) * t364 + (-m(6) * t247 - m(7) * t346 + t455) * g(3) + t452 * (pkin(10) + t430) - t391 * (pkin(4) * t365 - t58) + (m(6) * t431 + t479) * t477; (t259 * t293 - t336) * g(2) + t391 * t57 + (-t402 / 0.2e1 - t405 / 0.2e1 - t403 / 0.2e1 + t485) * t315 - m(7) * (t22 * t28 + t23 * t29 + t53 * t57) - t8 * t450 + t290 + (t463 - t482) * g(3) + (-t306 / 0.2e1 - t308 / 0.2e1 - t307 / 0.2e1 + t320 * mrSges(7,3) + t486) * t337 + t295 * mrSges(7,3) + (m(7) * (-t423 + t424 + t461) + t480) * pkin(10) + (t260 * t293 - t335) * g(1) - t56 * t141 - t29 * t90 - t28 * t91 - pkin(5) * t16; -t53 * (mrSges(7,1) * t134 + mrSges(7,2) * t133) + (Ifges(7,1) * t133 - t404) * t441 + t69 * t440 + (Ifges(7,5) * t133 - Ifges(7,6) * t134) * t439 - t22 * t90 + t23 * t91 - g(1) * (mrSges(7,1) * t182 - mrSges(7,2) * t183) - g(2) * (-mrSges(7,1) * t180 + mrSges(7,2) * t181) + g(3) * t328 * t248 + (t133 * t22 + t134 * t23) * mrSges(7,3) + t360 + (-Ifges(7,2) * t134 + t132 + t70) * t442 + t454;];
tau  = t1;
