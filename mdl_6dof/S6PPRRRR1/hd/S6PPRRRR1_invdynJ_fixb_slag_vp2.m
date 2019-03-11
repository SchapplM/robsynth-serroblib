% Calculate vector of inverse dynamics joint torques for
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:45
% EndTime: 2019-03-08 19:00:05
% DurationCPUTime: 10.96s
% Computational Cost: add. (9990->616), mult. (24469->887), div. (0->0), fcn. (21501->18), ass. (0->293)
t447 = mrSges(6,2) - mrSges(7,3);
t258 = sin(qJ(6));
t262 = cos(qJ(6));
t304 = -mrSges(7,1) * t262 + mrSges(7,2) * t258;
t443 = -mrSges(6,1) + t304;
t259 = sin(qJ(5));
t260 = sin(qJ(4));
t263 = cos(qJ(5));
t264 = cos(qJ(4));
t221 = t259 * t260 - t263 * t264;
t265 = -pkin(10) - pkin(9);
t333 = qJD(4) * t265;
t229 = t260 * t333;
t375 = cos(pkin(6));
t241 = qJD(1) * t375 + qJD(2);
t261 = sin(qJ(3));
t256 = sin(pkin(6));
t257 = cos(pkin(13));
t374 = cos(pkin(7));
t322 = t257 * t374;
t397 = cos(qJ(3));
t288 = t397 * t322;
t285 = t256 * t288;
t277 = qJD(1) * t285;
t254 = sin(pkin(13));
t357 = t254 * t256;
t331 = qJD(1) * t357;
t255 = sin(pkin(7));
t332 = t255 * t397;
t268 = -t241 * t332 + t261 * t331 - t277;
t238 = t265 * t260;
t239 = t265 * t264;
t291 = t263 * t238 + t239 * t259;
t317 = t264 * t333;
t429 = qJD(5) * t291 - t221 * t268 + t263 * t229 + t259 * t317;
t310 = t261 * t322;
t275 = (t254 * t397 + t310) * t256;
t355 = t255 * t261;
t159 = qJD(1) * t275 + t241 * t355;
t281 = t221 * qJD(5);
t179 = -qJD(4) * t221 - t281;
t222 = t259 * t264 + t260 * t263;
t282 = t222 * qJD(5);
t180 = qJD(4) * t222 + t282;
t343 = qJD(4) * t260;
t338 = pkin(4) * t343;
t448 = pkin(5) * t180 - pkin(11) * t179 - t159 + t338;
t192 = t238 * t259 - t239 * t263;
t432 = qJD(5) * t192 + t222 * t268 + t229 * t259 - t263 * t317;
t340 = qJD(3) * qJD(4);
t230 = qJDD(3) * t264 - t260 * t340;
t231 = qJDD(3) * t260 + t264 * t340;
t142 = -qJD(3) * t281 + t230 * t259 + t231 * t263;
t143 = -qJD(3) * t282 + t230 * t263 - t231 * t259;
t57 = -mrSges(6,1) * t143 + mrSges(6,2) * t142;
t446 = mrSges(4,1) * qJDD(3) - t57;
t325 = t255 * t375;
t175 = t261 * t325 + t275;
t336 = t256 * t257 * t255;
t204 = t374 * t375 - t336;
t132 = -t175 * t260 + t204 * t264;
t372 = sin(pkin(12));
t307 = t375 * t372;
t373 = cos(pkin(12));
t206 = -t254 * t307 + t257 * t373;
t274 = t254 * t373 + t257 * t307;
t323 = t256 * t372;
t427 = -t255 * t323 + t274 * t374;
t129 = t206 * t397 - t261 * t427;
t178 = t255 * t274 + t323 * t374;
t445 = -t129 * t260 + t178 * t264;
t308 = t375 * t373;
t205 = t254 * t308 + t257 * t372;
t273 = t254 * t372 - t257 * t308;
t324 = t256 * t373;
t428 = t255 * t324 + t273 * t374;
t127 = t205 * t397 - t261 * t428;
t177 = t255 * t273 - t324 * t374;
t444 = -t127 * t260 + t177 * t264;
t148 = qJD(3) * pkin(9) + t159;
t320 = pkin(10) * qJD(3) + t148;
t198 = -qJD(1) * t336 + t241 * t374;
t361 = t198 * t260;
t98 = t264 * t320 + t361;
t380 = t259 * t98;
t190 = t264 * t198;
t97 = -t260 * t320 + t190;
t96 = qJD(4) * pkin(4) + t97;
t40 = t263 * t96 - t380;
t344 = qJD(3) * t264;
t442 = Ifges(5,4) * t344 / 0.2e1 + Ifges(5,5) * qJD(4) / 0.2e1;
t250 = qJDD(4) + qJDD(5);
t104 = t148 * t264 + t361;
t240 = qJDD(1) * t375 + qJDD(2);
t195 = -qJDD(1) * t336 + t240 * t374;
t290 = t256 * t310;
t315 = qJD(3) * t331;
t316 = qJD(3) * t332;
t327 = qJDD(1) * t357;
t92 = qJD(3) * t277 + qJDD(1) * t290 + t240 * t355 + t241 * t316 - t261 * t315 + t397 * t327;
t88 = qJDD(3) * pkin(9) + t92;
t32 = -qJD(4) * t104 + t264 * t195 - t260 * t88;
t28 = qJDD(4) * pkin(4) - pkin(10) * t231 + t32;
t342 = qJD(4) * t264;
t31 = -t148 * t343 + t260 * t195 + t198 * t342 + t264 * t88;
t29 = pkin(10) * t230 + t31;
t377 = t263 * t98;
t41 = t259 * t96 + t377;
t9 = -qJD(5) * t41 - t259 * t29 + t263 * t28;
t6 = -pkin(5) * t250 - t9;
t410 = m(7) * t6;
t213 = t222 * qJD(3);
t252 = qJD(4) + qJD(5);
t188 = -t213 * t258 + t252 * t262;
t82 = qJD(6) * t188 + t142 * t262 + t250 * t258;
t408 = t82 / 0.2e1;
t189 = t213 * t262 + t252 * t258;
t83 = -qJD(6) * t189 - t142 * t258 + t250 * t262;
t407 = t83 / 0.2e1;
t441 = -m(6) - m(7);
t141 = qJDD(6) - t143;
t406 = t141 / 0.2e1;
t389 = Ifges(5,4) * t260;
t440 = t389 / 0.2e1;
t245 = pkin(4) * t264 + pkin(3);
t172 = pkin(5) * t221 - pkin(11) * t222 - t245;
t107 = t172 * t258 + t192 * t262;
t439 = -qJD(6) * t107 - t429 * t258 + t448 * t262;
t106 = t172 * t262 - t192 * t258;
t438 = qJD(6) * t106 + t448 * t258 + t429 * t262;
t437 = t188 * Ifges(7,6);
t212 = t221 * qJD(3);
t207 = qJD(6) + t212;
t436 = t207 * Ifges(7,3);
t435 = t252 * Ifges(6,5);
t434 = t252 * Ifges(6,6);
t303 = mrSges(7,1) * t258 + mrSges(7,2) * t262;
t38 = -pkin(5) * t252 - t40;
t433 = t303 * t38;
t130 = mrSges(6,1) * t250 - mrSges(6,3) * t142;
t35 = -mrSges(7,1) * t83 + mrSges(7,2) * t82;
t431 = t35 - t130;
t430 = t410 + t35;
t289 = t397 * t325;
t356 = t254 * t261;
t174 = t256 * t356 - t285 - t289;
t169 = t174 * t262;
t133 = t175 * t264 + t204 * t260;
t51 = t132 * t259 + t133 * t263;
t44 = -t258 * t51 + t169;
t392 = mrSges(6,3) * t213;
t349 = mrSges(6,1) * t252 + mrSges(7,1) * t188 - mrSges(7,2) * t189 - t392;
t426 = -t260 * t32 + t264 * t31;
t39 = pkin(11) * t252 + t41;
t146 = -qJD(3) * t245 + t268;
t94 = t212 * pkin(5) - t213 * pkin(11) + t146;
t14 = -t258 * t39 + t262 * t94;
t15 = t258 * t94 + t262 * t39;
t424 = qJD(6) * (-t14 * t262 - t15 * t258);
t301 = t264 * Ifges(5,2) + t389;
t370 = Ifges(5,6) * qJD(4);
t423 = -t104 * mrSges(5,3) - qJD(3) * t301 / 0.2e1 - t370 / 0.2e1;
t103 = -t148 * t260 + t190;
t346 = qJD(3) * t260;
t422 = -t103 * mrSges(5,3) + Ifges(5,1) * t346 / 0.2e1 + t442;
t184 = -mrSges(5,1) * t230 + mrSges(5,2) * t231;
t345 = qJD(3) * t261;
t330 = t255 * t345;
t93 = -qJD(3) * qJD(1) * t290 + qJDD(1) * t285 + t240 * t332 - t241 * t330 - t261 * t327 - t397 * t315;
t89 = -qJDD(3) * pkin(3) - t93;
t421 = -m(5) * t89 - t184;
t168 = mrSges(6,1) * t212 + mrSges(6,2) * t213;
t420 = -m(6) * t146 - t168;
t419 = m(5) + m(4) + m(3) - t441;
t418 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - t303;
t253 = qJ(4) + qJ(5);
t248 = sin(t253);
t249 = cos(t253);
t306 = -mrSges(5,1) * t264 + mrSges(5,2) * t260;
t417 = m(5) * pkin(3) + mrSges(4,1) - t306 + (m(7) * pkin(5) - t443) * t249 + (m(7) * pkin(11) - t447) * t248;
t72 = -t230 * pkin(4) + t89;
t30 = -t143 * pkin(5) - t142 * pkin(11) + t72;
t8 = qJD(5) * t40 + t259 * t28 + t263 * t29;
t5 = pkin(11) * t250 + t8;
t2 = qJD(6) * t14 + t258 * t30 + t262 * t5;
t3 = -qJD(6) * t15 - t258 * t5 + t262 * t30;
t416 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t144 = -mrSges(7,2) * t207 + mrSges(7,3) * t188;
t145 = mrSges(7,1) * t207 - mrSges(7,3) * t189;
t276 = -t3 * t258 + t424;
t46 = mrSges(7,1) * t141 - mrSges(7,3) * t82;
t47 = -mrSges(7,2) * t141 + mrSges(7,3) * t83;
t415 = (-t144 * t258 - t145 * t262) * qJD(6) + m(7) * (t2 * t262 + t276) - t258 * t46 + t262 * t47;
t414 = m(6) * t9 - t410 - t431;
t413 = -m(6) * t40 + m(7) * t38 - t349;
t147 = -qJD(3) * pkin(3) + t268;
t224 = t306 * qJD(3);
t412 = -m(5) * t147 + qJD(3) * mrSges(4,1) - t224 + t420;
t411 = m(6) / 0.2e1;
t409 = Ifges(7,1) * t408 + Ifges(7,4) * t407 + Ifges(7,5) * t406;
t405 = -t188 / 0.2e1;
t404 = -t189 / 0.2e1;
t403 = t189 / 0.2e1;
t402 = -t207 / 0.2e1;
t400 = -t212 / 0.2e1;
t398 = t213 / 0.2e1;
t393 = mrSges(6,3) * t212;
t391 = mrSges(7,3) * t258;
t390 = mrSges(7,3) * t262;
t388 = Ifges(5,4) * t264;
t387 = Ifges(7,4) * t258;
t386 = Ifges(7,4) * t262;
t383 = t189 * Ifges(7,4);
t382 = t213 * Ifges(6,4);
t368 = qJD(3) * mrSges(4,2);
t365 = t174 * t258;
t359 = t222 * t258;
t358 = t222 * t262;
t101 = t188 * Ifges(7,2) + t207 * Ifges(7,6) + t383;
t354 = t258 * t101;
t183 = Ifges(7,4) * t188;
t102 = t189 * Ifges(7,1) + t207 * Ifges(7,5) + t183;
t353 = t262 * t102;
t347 = mrSges(4,2) * qJDD(3);
t341 = qJD(6) * t258;
t339 = Ifges(7,5) * t82 + Ifges(7,6) * t83 + Ifges(7,3) * t141;
t70 = -t129 * t248 + t178 * t249;
t71 = t129 * t249 + t178 * t248;
t335 = t70 * pkin(5) + pkin(11) * t71;
t329 = qJD(6) * t358;
t328 = t353 / 0.2e1;
t326 = -t341 / 0.2e1;
t314 = t444 * pkin(4);
t313 = t445 * pkin(4);
t312 = t132 * pkin(4);
t170 = pkin(5) * t213 + pkin(11) * t212;
t302 = Ifges(7,1) * t262 - t387;
t300 = -Ifges(7,2) * t258 + t386;
t299 = Ifges(7,5) * t262 - Ifges(7,6) * t258;
t45 = t262 * t51 + t365;
t295 = t263 * t132 - t133 * t259;
t208 = -t260 * t355 + t264 * t374;
t209 = t260 * t374 + t264 * t355;
t293 = t263 * t208 - t209 * t259;
t157 = t208 * t259 + t209 * t263;
t287 = t179 * t258 + t329;
t286 = -t179 * t262 + t222 * t341;
t68 = -t127 * t248 + t177 * t249;
t69 = t127 * t249 + t177 * t248;
t284 = t443 * t68 + t447 * t69;
t283 = t443 * t70 + t447 * t71;
t136 = -t258 * t157 - t262 * t332;
t279 = -t262 * t157 + t258 * t332;
t116 = -t175 * t248 + t204 * t249;
t117 = t175 * t249 + t204 * t248;
t278 = t443 * t116 + t117 * t447;
t100 = t189 * Ifges(7,5) + t436 + t437;
t153 = -t212 * Ifges(6,2) + t382 + t434;
t203 = Ifges(6,4) * t212;
t154 = t213 * Ifges(6,1) - t203 + t435;
t18 = t82 * Ifges(7,4) + t83 * Ifges(7,2) + t141 * Ifges(7,6);
t267 = (t188 * t300 + t189 * t302 + t207 * t299) * qJD(6) / 0.2e1 - (-Ifges(6,1) * t212 + t100 - t382) * t213 / 0.2e1 + (-Ifges(6,2) * t213 + t154 - t203 + t353) * t212 / 0.2e1 + (t328 + t433) * qJD(6) - t15 * (-mrSges(7,2) * t213 + t212 * t391) - t14 * (mrSges(7,1) * t213 + t212 * t390) - t146 * (mrSges(6,1) * t213 - mrSges(6,2) * t212) + (Ifges(7,3) * t213 - t212 * t299) * t402 + (Ifges(7,5) * t213 - t212 * t302) * t404 + (Ifges(7,6) * t213 - t212 * t300) * t405 - t252 * (-Ifges(6,5) * t212 - Ifges(6,6) * t213) / 0.2e1 + (Ifges(7,5) * t258 + Ifges(7,6) * t262) * t406 + (Ifges(7,2) * t262 + t387) * t407 + (Ifges(7,1) * t258 + t386) * t408 + t258 * t409 + t2 * t390 - t40 * t393 + t153 * t398 + t354 * t400 + t262 * t18 / 0.2e1 + Ifges(6,3) * t250 + Ifges(6,5) * t142 + Ifges(6,6) * t143 - t8 * mrSges(6,2) + t9 * mrSges(6,1) + t101 * t326 + t6 * t304 + t212 * t433;
t266 = qJD(3) ^ 2;
t236 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t344;
t235 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t346;
t202 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t231;
t201 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t230;
t196 = -mrSges(6,2) * t252 - t393;
t182 = -qJD(4) * t209 - t260 * t316;
t181 = qJD(4) * t208 + t264 * t316;
t164 = t175 * qJD(3);
t163 = (t289 + (t288 - t356) * t256) * qJD(3);
t152 = pkin(4) * t346 + t170;
t131 = -mrSges(6,2) * t250 + mrSges(6,3) * t143;
t128 = t206 * t261 + t397 * t427;
t126 = t205 * t261 + t397 * t428;
t115 = t116 * pkin(5);
t74 = qJD(4) * t132 + t163 * t264;
t73 = -qJD(4) * t133 - t163 * t260;
t64 = t68 * pkin(5);
t52 = qJD(5) * t293 + t181 * t263 + t182 * t259;
t43 = t263 * t97 - t380;
t42 = t259 * t97 + t377;
t37 = qJD(6) * t279 - t258 * t52 + t262 * t330;
t36 = qJD(6) * t136 + t258 * t330 + t262 * t52;
t24 = t170 * t258 + t262 * t40;
t23 = t170 * t262 - t258 * t40;
t21 = t152 * t258 + t262 * t43;
t20 = t152 * t262 - t258 * t43;
t12 = qJD(5) * t295 + t259 * t73 + t263 * t74;
t11 = -qJD(6) * t45 - t12 * t258 + t164 * t262;
t10 = qJD(6) * t44 + t12 * t262 + t164 * t258;
t1 = [m(6) * (t12 * t41 + t51 * t8) + m(5) * (t103 * t73 + t104 * t74 + t132 * t32 + t133 * t31) + m(4) * (t159 * t163 + t175 * t92 + t195 * t204) + m(3) * (t240 * t375 + (t254 ^ 2 + t257 ^ 2) * t256 ^ 2 * qJDD(1)) - t163 * t368 - t175 * t347 + t74 * t236 + t73 * t235 + t133 * t201 + t132 * t202 + t12 * t196 + t10 * t144 + t11 * t145 + t51 * t131 + t44 * t46 + t45 * t47 + m(7) * (t10 * t15 + t11 * t14 + t2 * t45 + t3 * t44) + m(2) * qJDD(1) + t414 * t295 + t413 * (qJD(5) * t51 + t259 * t74 - t263 * t73) + (-m(4) * t93 + m(6) * t72 - t421 - t446) * t174 + (m(4) * t268 - t412) * t164 + (-m(2) - t419) * g(3); m(7) * (t136 * t3 + t14 * t37 + t15 * t36 - t2 * t279) + m(3) * t240 + m(4) * (t195 * t374 + (t397 * t93 + t261 * t92 + (t159 * t397 + t261 * t268) * qJD(3)) * t255) + m(5) * (t103 * t182 + t104 * t181 + t32 * t208 + t31 * t209 + (t147 * t345 - t397 * t89) * t255) + m(6) * (t8 * t157 + t41 * t52 + (t146 * t345 - t397 * t72) * t255) + t181 * t236 + t182 * t235 + t208 * t202 + t209 * t201 + t52 * t196 + t157 * t131 + t36 * t144 + t37 * t145 + t136 * t46 - t279 * t47 + t413 * (qJD(5) * t157 + t181 * t259 - t263 * t182) + (-mrSges(4,1) * t266 - t347) * t355 + (t224 + t168) * t330 + t414 * t293 + (-mrSges(4,2) * t266 - t184 + t446) * t332 + (-g(1) * t323 + g(2) * t324 - g(3) * t375) * t419; -(t368 - m(5) * (-t103 * t260 + t104 * t264) + t260 * t235 - t264 * t236) * t268 + (t106 * t3 + t107 * t2 + t14 * t439 + t15 * t438 - t291 * t6 + t38 * t432) * m(7) - t431 * t291 + (t146 * t338 + t192 * t8 - t245 * t72 + t291 * t9 - t40 * t432 + t41 * t429) * m(6) - t432 * t349 + (m(5) * ((-t103 * t264 - t104 * t260) * qJD(4) + t426) - t235 * t342 - t236 * t343 - t260 * t202 + t264 * t201) * pkin(9) + t426 * mrSges(5,3) + (-Ifges(7,1) * t286 - Ifges(7,4) * t287) * t403 + t38 * (mrSges(7,1) * t287 - mrSges(7,2) * t286) + t422 * t342 + t423 * t343 + t421 * pkin(3) + t412 * t159 + t231 * t260 * Ifges(5,1) + (t264 * (-Ifges(5,2) * t260 + t388) + t260 * (Ifges(5,1) * t264 - t389)) * t340 / 0.2e1 + t231 * t388 / 0.2e1 + t358 * t409 + t188 * (-Ifges(7,4) * t286 - Ifges(7,2) * t287) / 0.2e1 + t207 * (-Ifges(7,5) * t286 - Ifges(7,6) * t287) / 0.2e1 + qJD(4) ^ 2 * (Ifges(5,5) * t264 - Ifges(5,6) * t260) / 0.2e1 + (t72 * mrSges(6,2) - t9 * mrSges(6,3) + Ifges(6,1) * t142 + Ifges(6,4) * t143 + Ifges(6,5) * t250 + t102 * t326 + t299 * t406 + t300 * t407 + t302 * t408 + t303 * t6) * t222 + t168 * t338 + (t441 * (-t174 * t245 - t175 * t265) + t418 * t175 + t417 * t174) * g(3) + (t441 * (-t128 * t245 - t129 * t265) + t418 * t129 + t417 * t128) * g(1) + (t441 * (-t126 * t245 - t127 * t265) + t418 * t127 + t417 * t126) * g(2) + t439 * t145 + t438 * t144 + (t14 * mrSges(7,1) + t436 / 0.2e1 + t437 / 0.2e1 - t15 * mrSges(7,2) - Ifges(6,4) * t398 - Ifges(6,2) * t400 + Ifges(7,5) * t403 - t434 / 0.2e1 + t146 * mrSges(6,1) - t153 / 0.2e1 + t100 / 0.2e1 - t41 * mrSges(6,3)) * t180 + (Ifges(6,1) * t398 + Ifges(6,4) * t400 - t354 / 0.2e1 + t435 / 0.2e1 + t146 * mrSges(6,2) + t154 / 0.2e1 - t40 * mrSges(6,3) + t328) * t179 + t147 * (mrSges(5,1) * t260 + mrSges(5,2) * t264) * qJD(4) + (t14 * t286 - t15 * t287 - t2 * t359 - t3 * t358) * mrSges(7,3) + t429 * t196 + (Ifges(7,3) * t406 + Ifges(7,6) * t407 + Ifges(7,5) * t408 + t339 / 0.2e1 - Ifges(6,6) * t250 - Ifges(6,4) * t142 - Ifges(6,2) * t143 + t72 * mrSges(6,1) - t8 * mrSges(6,3) + t416) * t221 + t264 * (Ifges(5,4) * t231 + Ifges(5,2) * t230) / 0.2e1 + Ifges(4,3) * qJDD(3) - t18 * t359 / 0.2e1 + qJDD(4) * (Ifges(5,5) * t260 + Ifges(5,6) * t264) - t245 * t57 - t101 * t329 / 0.2e1 + t192 * t131 + t106 * t46 + t107 * t47 - t92 * mrSges(4,2) + t93 * mrSges(4,1) + t230 * t301 / 0.2e1 + t89 * t306 + t230 * t440; ((-t147 * mrSges(5,2) - t422 - t442) * t264 + (-t147 * mrSges(5,1) + t370 / 0.2e1 + (t440 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t264) * qJD(3) + t420 * pkin(4) - t423) * t260) * qJD(3) + (-t444 * mrSges(5,1) - (-t127 * t264 - t177 * t260) * mrSges(5,2) - m(7) * (pkin(11) * t69 + t314 + t64) - m(6) * t314 + t284) * g(2) + (-t445 * mrSges(5,1) - (-t129 * t264 - t178 * t260) * mrSges(5,2) - m(7) * (t313 + t335) - m(6) * t313 + t283) * g(1) + (t130 * t263 + t131 * t259) * pkin(4) - m(7) * (t14 * t20 + t15 * t21 + t38 * t42) + mrSges(7,3) * t424 + t415 * (pkin(4) * t259 + pkin(11)) + t41 * t392 + 0.2e1 * ((t259 * t8 + t263 * t9) * t411 + ((-t259 * t40 + t263 * t41) * t411 + m(7) * (t259 * t38 + (-t14 * t258 + t15 * t262) * t263) / 0.2e1) * qJD(5)) * pkin(4) + t267 + t430 * (-pkin(4) * t263 - pkin(5)) - m(6) * (-t40 * t42 + t41 * t43) + Ifges(5,3) * qJDD(4) - t3 * t391 + t349 * t42 - t103 * t236 + Ifges(5,6) * t230 + Ifges(5,5) * t231 + t104 * t235 - t43 * t196 - t21 * t144 - t20 * t145 - t31 * mrSges(5,2) + t32 * mrSges(5,1) + (-t349 * t259 + (t144 * t262 - t145 * t258 + t196) * t263) * pkin(4) * qJD(5) + (-m(6) * t312 - m(7) * (pkin(11) * t117 + t115 + t312) - mrSges(5,1) * t132 + mrSges(5,2) * t133 + t278) * g(3); t276 * mrSges(7,3) - m(7) * (t14 * t23 + t15 * t24 + t38 * t41) + t267 + (t349 + t392) * t41 + (-m(7) * t64 + t284) * g(2) + (-m(7) * t335 + t283) * g(1) + (-m(7) * t115 + t278) * g(3) - t40 * t196 - t24 * t144 - t23 * t145 + ((-g(2) * t69 - g(3) * t117) * m(7) + t415) * pkin(11) - t430 * pkin(5); -t38 * (mrSges(7,1) * t189 + mrSges(7,2) * t188) + (Ifges(7,1) * t188 - t383) * t404 + t101 * t403 + (Ifges(7,5) * t188 - Ifges(7,6) * t189) * t402 - t14 * t144 + t15 * t145 - g(1) * ((t128 * t262 - t258 * t71) * mrSges(7,1) + (-t128 * t258 - t262 * t71) * mrSges(7,2)) - g(2) * ((t126 * t262 - t258 * t69) * mrSges(7,1) + (-t126 * t258 - t262 * t69) * mrSges(7,2)) - g(3) * ((-t117 * t258 + t169) * mrSges(7,1) + (-t117 * t262 - t365) * mrSges(7,2)) + (t14 * t188 + t15 * t189) * mrSges(7,3) + t339 + (-Ifges(7,2) * t189 + t102 + t183) * t405 + t416;];
tau  = t1;
