% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:56
% EndTime: 2019-03-08 20:41:25
% DurationCPUTime: 14.94s
% Computational Cost: add. (7108->595), mult. (14323->828), div. (0->0), fcn. (10253->14), ass. (0->284)
t437 = mrSges(6,2) - mrSges(7,3);
t241 = cos(qJ(2));
t231 = sin(pkin(6));
t332 = qJD(1) * t231;
t298 = t241 * t332;
t274 = qJD(3) - t298;
t229 = qJ(4) + qJ(5);
t224 = sin(t229);
t225 = cos(t229);
t236 = sin(qJ(4));
t240 = cos(qJ(4));
t284 = mrSges(5,1) * t236 + mrSges(5,2) * t240;
t405 = m(7) * pkin(10);
t406 = m(7) * pkin(5);
t426 = mrSges(3,2) - mrSges(4,3);
t234 = sin(qJ(6));
t238 = cos(qJ(6));
t282 = -mrSges(7,1) * t238 + mrSges(7,2) * t234;
t436 = -mrSges(6,1) + t282;
t440 = -t284 + t426 + (-t406 + t436) * t224 + (t405 - t437) * t225;
t325 = qJD(4) * t240;
t416 = pkin(4) * t325 + t274;
t401 = -m(5) - m(4);
t427 = m(6) + m(7);
t438 = t401 - t427;
t243 = -pkin(2) - pkin(8);
t383 = pkin(9) - t243;
t195 = t383 * t236;
t196 = t383 * t240;
t235 = sin(qJ(5));
t239 = cos(qJ(5));
t131 = -t195 * t235 + t239 * t196;
t185 = qJD(4) * t196;
t189 = t235 * t240 + t236 * t239;
t326 = qJD(4) * t236;
t285 = t383 * t326;
t237 = sin(qJ(2));
t299 = t237 * t332;
t435 = qJD(5) * t131 + t239 * t185 + t189 * t299 - t235 * t285;
t322 = qJD(5) * t239;
t323 = qJD(5) * t235;
t125 = -t235 * t325 - t236 * t322 - t239 * t326 - t240 * t323;
t228 = qJD(4) + qJD(5);
t341 = t239 * t240;
t126 = t228 * t341 - t235 * t326 - t236 * t323;
t434 = pkin(5) * t126 - pkin(10) * t125 + t416;
t132 = -t195 * t239 - t196 * t235;
t287 = t240 * t299;
t288 = t236 * t299;
t422 = -qJD(5) * t132 + (t285 - t287) * t239 + (t185 + t288) * t235;
t328 = qJD(2) * t240;
t330 = qJD(2) * t236;
t183 = -t235 * t330 + t239 * t328;
t139 = -t183 * t234 + t228 * t238;
t140 = t183 * t238 + t228 * t234;
t378 = mrSges(6,3) * t183;
t362 = mrSges(6,1) * t228 + mrSges(7,1) * t139 - mrSges(7,2) * t140 - t378;
t230 = sin(pkin(11));
t232 = cos(pkin(11));
t233 = cos(pkin(6));
t346 = t233 * t241;
t173 = t230 * t237 - t232 * t346;
t351 = t231 * t236;
t433 = t173 * t240 + t232 * t351;
t175 = t230 * t346 + t232 * t237;
t432 = t175 * t240 - t230 * t351;
t179 = qJD(2) * t243 + t274;
t331 = qJD(1) * t233;
t297 = t236 * t331;
t135 = t240 * t179 - t297;
t296 = t240 * t331;
t136 = t179 * t236 + t296;
t318 = qJD(2) * qJD(4);
t192 = qJDD(2) * t240 - t236 * t318;
t167 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t192;
t193 = -qJDD(2) * t236 - t240 * t318;
t168 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t193;
t308 = mrSges(5,3) * t330;
t199 = -qJD(4) * mrSges(5,2) - t308;
t307 = mrSges(5,3) * t328;
t200 = qJD(4) * mrSges(5,1) - t307;
t329 = qJD(2) * t237;
t295 = t231 * t329;
t203 = qJD(1) * t295;
t348 = t231 * t241;
t165 = qJDD(1) * t348 - t203;
t261 = qJDD(3) - t165;
t141 = qJDD(2) * t243 + t261;
t316 = qJDD(1) * t233;
t63 = -qJD(4) * t297 + t236 * t141 + t179 * t325 + t240 * t316;
t64 = -qJD(4) * t136 + t240 * t141 - t236 * t316;
t414 = t236 * t63 + t240 * t64;
t430 = m(5) * ((-t135 * t236 + t136 * t240) * qJD(4) + t414) + (t240 * t199 - t236 * t200) * qJD(4) + t240 * t167 + t236 * t168;
t182 = t189 * qJD(2);
t171 = qJD(6) + t182;
t371 = t140 * Ifges(7,4);
t56 = t139 * Ifges(7,2) + t171 * Ifges(7,6) + t371;
t429 = -t56 / 0.2e1;
t133 = Ifges(7,4) * t139;
t57 = t140 * Ifges(7,1) + t171 * Ifges(7,5) + t133;
t428 = t57 / 0.2e1;
t385 = mrSges(3,1) - mrSges(4,2);
t227 = qJDD(4) + qJDD(5);
t82 = -qJD(5) * t182 + t192 * t239 + t193 * t235;
t40 = qJD(6) * t139 + t227 * t234 + t238 * t82;
t41 = -qJD(6) * t140 + t227 * t238 - t234 * t82;
t16 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t71 = mrSges(6,1) * t227 - mrSges(6,3) * t82;
t382 = t16 - t71;
t188 = t235 * t236 - t341;
t226 = t236 * pkin(4);
t216 = qJ(3) + t226;
t116 = pkin(5) * t189 + pkin(10) * t188 + t216;
t62 = t116 * t234 + t132 * t238;
t425 = -qJD(6) * t62 + t234 * t435 + t238 * t434;
t61 = t116 * t238 - t132 * t234;
t424 = qJD(6) * t61 + t234 * t434 - t238 * t435;
t281 = mrSges(7,1) * t234 + mrSges(7,2) * t238;
t121 = -pkin(9) * t328 + t135;
t104 = qJD(4) * pkin(4) + t121;
t122 = t296 + (-pkin(9) * qJD(2) + t179) * t236;
t361 = t122 * t235;
t52 = t104 * t239 - t361;
t47 = -pkin(5) * t228 - t52;
t423 = t281 * t47;
t350 = t231 * t237;
t204 = t238 * t350;
t177 = -t233 * t236 - t240 * t348;
t262 = -t233 * t240 + t236 * t348;
t97 = t177 * t235 - t239 * t262;
t75 = -t234 * t97 + t204;
t157 = -t224 * t233 - t225 * t348;
t158 = -t224 * t348 + t225 * t233;
t419 = t157 * t436 + t158 * t437;
t352 = t231 * t232;
t119 = t173 * t225 + t224 * t352;
t120 = -t173 * t224 + t225 * t352;
t418 = t119 * t436 - t120 * t437;
t353 = t230 * t231;
t117 = t175 * t225 - t224 * t353;
t118 = t175 * t224 + t225 * t353;
t417 = t117 * t436 + t118 * t437;
t83 = qJD(2) * qJD(5) * t188 - t192 * t235 + t193 * t239;
t81 = qJDD(6) - t83;
t19 = mrSges(7,1) * t81 - mrSges(7,3) * t40;
t20 = -mrSges(7,2) * t81 + mrSges(7,3) * t41;
t413 = -t234 * t19 + t238 * t20;
t379 = mrSges(6,3) * t182;
t144 = -mrSges(6,2) * t228 - t379;
t84 = -mrSges(7,2) * t171 + mrSges(7,3) * t139;
t85 = mrSges(7,1) * t171 - mrSges(7,3) * t140;
t412 = -t234 * t85 + t238 * t84 + t144;
t46 = qJDD(4) * pkin(4) - pkin(9) * t192 + t64;
t51 = pkin(9) * t193 + t63;
t342 = t239 * t122;
t53 = t104 * t235 + t342;
t9 = -qJD(5) * t53 - t235 * t51 + t239 * t46;
t8 = t104 * t322 - t122 * t323 + t235 * t46 + t239 * t51;
t410 = -t125 * t52 - t126 * t53 + t188 * t9 - t189 * t8;
t48 = pkin(10) * t228 + t53;
t191 = qJD(2) * qJ(3) + t299;
t169 = pkin(4) * t330 + t191;
t78 = pkin(5) * t182 - pkin(10) * t183 + t169;
t21 = -t234 * t48 + t238 * t78;
t314 = qJDD(2) * qJ(3);
t315 = qJDD(1) * t237;
t142 = t231 * t315 + t314 + (qJD(3) + t298) * qJD(2);
t98 = -pkin(4) * t193 + t142;
t23 = -pkin(5) * t83 - pkin(10) * t82 + t98;
t5 = pkin(10) * t227 + t8;
t2 = qJD(6) * t21 + t23 * t234 + t238 * t5;
t22 = t234 * t78 + t238 * t48;
t3 = -qJD(6) * t22 + t23 * t238 - t234 * t5;
t409 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t408 = m(5) * pkin(8) + mrSges(5,3) + mrSges(6,3) + t281 + t385;
t407 = qJ(3) * t438 + t440;
t244 = qJD(2) ^ 2;
t404 = t40 / 0.2e1;
t403 = t41 / 0.2e1;
t402 = t81 / 0.2e1;
t400 = -t139 / 0.2e1;
t399 = -t140 / 0.2e1;
t398 = t140 / 0.2e1;
t397 = -t171 / 0.2e1;
t394 = t183 / 0.2e1;
t393 = pkin(4) * t235;
t392 = pkin(4) * t239;
t6 = -pkin(5) * t227 - t9;
t391 = t188 * t6;
t388 = t2 * t238;
t387 = t234 * t3;
t377 = Ifges(5,4) * t236;
t376 = Ifges(5,4) * t240;
t375 = Ifges(7,4) * t234;
t374 = Ifges(7,4) * t238;
t370 = t183 * Ifges(6,4);
t358 = t182 * t234;
t357 = t182 * t238;
t356 = t188 * t234;
t355 = t188 * t238;
t354 = t191 * t241;
t349 = t231 * t240;
t347 = t233 * t237;
t343 = t238 * t125;
t112 = mrSges(6,1) * t182 + mrSges(6,2) * t183;
t190 = t284 * qJD(2);
t336 = t112 + t190;
t334 = t433 * pkin(4);
t333 = pkin(2) * t348 + qJ(3) * t350;
t327 = qJD(2) * t241;
t321 = qJD(6) * t234;
t320 = qJD(6) * t238;
t319 = qJDD(2) * mrSges(4,2);
t313 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t81;
t312 = pkin(4) * t328;
t304 = t234 * t350;
t303 = t234 * t429;
t294 = t231 * t327;
t292 = qJD(1) * t327;
t291 = t320 / 0.2e1;
t290 = t117 * pkin(5) + pkin(10) * t118;
t289 = -t318 / 0.2e1;
t286 = t231 * t292;
t115 = pkin(5) * t183 + pkin(10) * t182;
t280 = t240 * Ifges(5,1) - t377;
t279 = Ifges(7,1) * t238 - t375;
t278 = -t236 * Ifges(5,2) + t376;
t277 = -Ifges(7,2) * t234 + t374;
t276 = -Ifges(5,5) * t236 - Ifges(5,6) * t240;
t275 = Ifges(7,5) * t238 - Ifges(7,6) * t234;
t273 = t21 * t238 + t22 * t234;
t272 = -t21 * t234 + t22 * t238;
t271 = -t234 * t84 - t238 * t85;
t269 = t239 * t177 + t235 * t262;
t267 = qJ(3) * t142 + qJD(3) * t191;
t266 = t432 * pkin(4);
t76 = t238 * t97 + t304;
t260 = t142 * t237 + t191 * t327;
t259 = -t125 * t234 + t188 * t320;
t258 = t188 * t321 + t343;
t257 = t191 * (mrSges(5,1) * t240 - mrSges(5,2) * t236);
t256 = t236 * (-Ifges(5,2) * t240 - t377);
t255 = t240 * (-Ifges(5,1) * t236 - t376);
t252 = t177 * pkin(4);
t250 = -qJD(6) * t273 - t387;
t249 = t250 + t388;
t248 = qJD(6) * t271 + t413;
t12 = t40 * Ifges(7,4) + t41 * Ifges(7,2) + t81 * Ifges(7,6);
t13 = t40 * Ifges(7,1) + t41 * Ifges(7,4) + t81 * Ifges(7,5);
t170 = Ifges(6,4) * t182;
t55 = t140 * Ifges(7,5) + t139 * Ifges(7,6) + t171 * Ifges(7,3);
t93 = -t182 * Ifges(6,2) + t228 * Ifges(6,6) + t370;
t94 = t183 * Ifges(6,1) + t228 * Ifges(6,5) - t170;
t245 = (t139 * t277 + t140 * t279 + t171 * t275) * qJD(6) / 0.2e1 - (-Ifges(6,1) * t182 - t370 + t55) * t183 / 0.2e1 - t52 * t379 + mrSges(7,3) * t388 + t6 * t282 + t93 * t394 + (Ifges(7,5) * t234 + Ifges(7,6) * t238) * t402 + (Ifges(7,2) * t238 + t375) * t403 + (Ifges(7,1) * t234 + t374) * t404 + (-Ifges(6,2) * t183 - t170 + t94) * t182 / 0.2e1 + t238 * t12 / 0.2e1 + t234 * t13 / 0.2e1 + Ifges(6,3) * t227 + Ifges(6,5) * t82 + Ifges(6,6) * t83 + (t303 + t423) * qJD(6) + t9 * mrSges(6,1) - t8 * mrSges(6,2) - t228 * (-Ifges(6,5) * t182 - Ifges(6,6) * t183) / 0.2e1 - t169 * (mrSges(6,1) * t183 - mrSges(6,2) * t182) - t21 * (mrSges(7,1) * t183 + mrSges(7,3) * t357) - t22 * (-mrSges(7,2) * t183 + mrSges(7,3) * t358) + (Ifges(7,3) * t183 - t182 * t275) * t397 + (Ifges(7,5) * t183 - t182 * t279) * t399 + (Ifges(7,6) * t183 - t182 * t277) * t400 + t57 * t291 + t182 * t423 + t357 * t428 + t358 * t429;
t242 = -pkin(9) - pkin(8);
t219 = -pkin(5) - t392;
t187 = -qJD(2) * pkin(2) + t274;
t181 = Ifges(5,5) * qJD(4) + qJD(2) * t280;
t180 = Ifges(5,6) * qJD(4) + qJD(2) * t278;
t176 = -t230 * t347 + t232 * t241;
t174 = t230 * t241 + t232 * t347;
t166 = (t292 + t315) * t231;
t164 = t175 * pkin(2);
t163 = t173 * pkin(2);
t151 = t157 * pkin(5);
t146 = -qJDD(2) * pkin(2) + t261;
t134 = -mrSges(5,1) * t193 + mrSges(5,2) * t192;
t128 = qJD(4) * t262 + t240 * t295;
t127 = qJD(4) * t177 + t236 * t295;
t114 = t119 * pkin(5);
t92 = t115 + t312;
t72 = -mrSges(6,2) * t227 + mrSges(6,3) * t83;
t59 = t121 * t239 - t361;
t58 = t121 * t235 + t342;
t30 = -mrSges(6,1) * t83 + mrSges(6,2) * t82;
t29 = qJD(5) * t97 + t127 * t235 - t239 * t128;
t28 = qJD(5) * t269 + t127 * t239 + t128 * t235;
t27 = t115 * t234 + t238 * t52;
t26 = t115 * t238 - t234 * t52;
t25 = t234 * t92 + t238 * t59;
t24 = -t234 * t59 + t238 * t92;
t18 = -qJD(6) * t76 - t234 * t28 + t238 * t294;
t17 = qJD(6) * t75 + t234 * t294 + t238 * t28;
t1 = [t127 * t199 + t128 * t200 + t28 * t144 + t177 * t167 - t262 * t168 + t17 * t84 + t18 * t85 + t75 * t19 + t76 * t20 + t97 * t72 - t382 * t269 - t362 * t29 + (-m(2) - m(3) + t438) * g(3) + m(5) * (t127 * t136 + t128 * t135 + t177 * t64 - t262 * t63) + m(6) * (t269 * t9 + t28 * t53 - t29 * t52 + t8 * t97) + m(7) * (t17 * t22 + t18 * t21 + t2 * t76 - t269 * t6 + t29 * t47 + t3 * t75) + ((qJD(2) * t336 + qJDD(2) * t385 - t244 * t426) * t241 + (-qJDD(2) * t426 - t244 * t385 + t134 + t30) * t237 + m(4) * (-t146 * t241 + t187 * t329 + t260) + m(3) * (t165 * t241 + t166 * t237) + m(5) * t260 + m(6) * (t169 * t327 + t237 * t98)) * t231 + (m(2) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * t233 ^ 2) * qJDD(1); (-t98 * mrSges(6,2) - Ifges(6,1) * t82 - Ifges(6,4) * t83 - Ifges(6,5) * t227 - t275 * t402 - t277 * t403 - t279 * t404 + t291 * t56) * t188 + (t135 * t326 - t136 * t325 - t414) * mrSges(5,3) + t416 * t112 - t236 * (Ifges(5,4) * t192 + Ifges(5,2) * t193) / 0.2e1 + (qJD(2) * qJD(3) + t142 - t286 + t314) * mrSges(4,3) + (t2 * t356 - t21 * t258 + t22 * t259 + t3 * t355) * mrSges(7,3) - t13 * t355 / 0.2e1 + t125 * t303 + t193 * t278 / 0.2e1 + t192 * t280 / 0.2e1 + t142 * t284 + t240 * (Ifges(5,1) * t192 + Ifges(5,4) * t193) / 0.2e1 + t274 * t190 + (Ifges(7,3) * t402 + Ifges(7,6) * t403 + Ifges(7,5) * t404 - Ifges(6,4) * t82 - Ifges(6,2) * t83 - Ifges(6,6) * t227 + t98 * mrSges(6,1) + t313 / 0.2e1 + t409) * t189 + t410 * mrSges(6,3) + (-t166 + t286) * mrSges(3,2) - t199 * t288 + (-(t354 + (t135 * t240 + t136 * t236) * t237) * t332 + t267) * m(5) + t430 * t243 - t200 * t287 + (-t203 + t146) * mrSges(4,2) - t180 * t325 / 0.2e1 + (Ifges(6,1) * t125 - Ifges(6,4) * t126) * t394 + (Ifges(7,1) * t258 + Ifges(7,4) * t259 + Ifges(7,5) * t126) * t398 + (-t427 * (t175 * t242 + t176 * t226 - t164) - t401 * t164 + t407 * t176 + t408 * t175) * g(1) + (-t427 * (t173 * t242 + t174 * t226 - t163) - t401 * t163 + t407 * t174 + t408 * t173) * g(2) + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) - t181 * t326 / 0.2e1 + (t401 * t333 - t427 * (t350 * t226 + t333) + ((t242 * t427 - t408) * t241 + t440 * t237) * t231) * g(3) + t255 * t318 / 0.2e1 - pkin(2) * t319 + qJDD(4) * (Ifges(5,5) * t240 - Ifges(5,6) * t236) + (t257 + t276 * qJD(4) / 0.2e1) * qJD(4) + t228 * (Ifges(6,5) * t125 - Ifges(6,6) * t126) / 0.2e1 + t216 * t30 - t182 * (Ifges(6,4) * t125 - Ifges(6,2) * t126) / 0.2e1 + t169 * (mrSges(6,1) * t126 + mrSges(6,2) * t125) + qJ(3) * t134 + t126 * t55 / 0.2e1 - t126 * t93 / 0.2e1 + t132 * t72 + t125 * t94 / 0.2e1 + t61 * t19 + t62 * t20 + (-(t187 * t237 + t354) * t332 - pkin(2) * t146 + t267) * m(4) - t435 * t144 + (-t131 * t9 + t132 * t8 + t169 * t416 + t216 * t98 + t422 * t52 - t435 * t53) * m(6) - t22 * mrSges(7,2) * t126 + t21 * mrSges(7,1) * t126 + (qJD(6) * t57 + t12) * t356 / 0.2e1 + t422 * t362 - t281 * t391 + t256 * t289 + t343 * t428 + t171 * (Ifges(7,5) * t258 + Ifges(7,6) * t259 + Ifges(7,3) * t126) / 0.2e1 + t139 * (Ifges(7,4) * t258 + Ifges(7,2) * t259 + Ifges(7,6) * t126) / 0.2e1 + t47 * (-mrSges(7,1) * t259 + mrSges(7,2) * t258) + t424 * t84 + t425 * t85 + (t131 * t6 + t2 * t62 + t21 * t425 + t22 * t424 + t3 * t61 - t422 * t47) * m(7) + t382 * t131 + (t165 + t203) * mrSges(3,1); t319 - t244 * mrSges(4,3) + t382 * t188 + t362 * t125 + t412 * t126 + (t72 + t248) * t189 + m(7) * (-t125 * t47 + t126 * t272 + t189 * t249 + t391) - m(6) * t410 + m(4) * t146 + (-m(6) * t169 - m(7) * t273 + t191 * t401 + t271 - t336) * qJD(2) - (-g(1) * t175 - g(2) * t173 + g(3) * t348) * t438 + t430; (m(7) * t249 - t320 * t85 - t321 * t84 + t413) * (pkin(10) + t393) + (t307 + t200) * t136 + t53 * t378 + (-t255 / 0.2e1 + t256 / 0.2e1) * t244 + t245 + t181 * t330 / 0.2e1 + (-t21 * t24 - t22 * t25 - t47 * t58 + t219 * t6 + (t235 * t47 + t239 * t272) * qJD(5) * pkin(4)) * m(7) + ((t235 * t8 + t239 * t9 + (-t235 * t52 + t239 * t53) * qJD(5)) * pkin(4) - t169 * t312 + t52 * t58 - t53 * t59) * m(6) + t412 * pkin(4) * t322 + (-t308 - t199) * t135 + t71 * t392 + t72 * t393 + (-t432 * mrSges(5,1) - (-t175 * t236 - t230 * t349) * mrSges(5,2) - m(6) * t266 - m(7) * (t266 + t290) + t417) * g(1) + (-m(7) * (-pkin(10) * t120 + t114 + t334) - m(6) * t334 - t433 * mrSges(5,1) - (-t173 * t236 + t232 * t349) * mrSges(5,2) + t418) * g(2) + (-m(6) * t252 - m(7) * (pkin(10) * t158 + t151 + t252) - mrSges(5,1) * t177 - mrSges(5,2) * t262 + t419) * g(3) + t180 * t328 / 0.2e1 + (-t21 * t320 - t22 * t321 - t387) * mrSges(7,3) + t219 * t16 + Ifges(5,5) * t192 + Ifges(5,6) * t193 + Ifges(5,3) * qJDD(4) - t59 * t144 - t25 * t84 - t24 * t85 - t63 * mrSges(5,2) + t64 * mrSges(5,1) - t362 * (pkin(4) * t323 - t58) + t276 * t289 - t112 * t312 - qJD(2) * t257; (t362 + t378) * t53 + t250 * mrSges(7,3) + (-m(7) * t290 + t417) * g(1) + (-m(7) * t114 + t418) * g(2) + (-m(7) * t151 + t419) * g(3) - m(7) * (t21 * t26 + t22 * t27 + t47 * t53) + t245 + ((g(2) * t120 - g(3) * t158) * m(7) + t248) * pkin(10) + t249 * t405 - t6 * t406 - t52 * t144 - t27 * t84 - t26 * t85 - pkin(5) * t16; -t47 * (mrSges(7,1) * t140 + mrSges(7,2) * t139) + (Ifges(7,1) * t139 - t371) * t399 + t56 * t398 + (Ifges(7,5) * t139 - Ifges(7,6) * t140) * t397 - t21 * t84 + t22 * t85 - g(1) * ((-t118 * t234 + t176 * t238) * mrSges(7,1) + (-t118 * t238 - t176 * t234) * mrSges(7,2)) - g(2) * ((t120 * t234 + t174 * t238) * mrSges(7,1) + (t120 * t238 - t174 * t234) * mrSges(7,2)) - g(3) * ((-t158 * t234 + t204) * mrSges(7,1) + (-t158 * t238 - t304) * mrSges(7,2)) + (t139 * t21 + t140 * t22) * mrSges(7,3) + t313 + (-Ifges(7,2) * t140 + t133 + t57) * t400 + t409;];
tau  = t1;
