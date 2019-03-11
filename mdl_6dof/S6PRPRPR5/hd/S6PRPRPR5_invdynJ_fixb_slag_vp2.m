% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:10
% EndTime: 2019-03-08 19:43:28
% DurationCPUTime: 13.56s
% Computational Cost: add. (5310->583), mult. (12589->752), div. (0->0), fcn. (9719->14), ass. (0->286)
t194 = sin(pkin(11));
t331 = pkin(8) + qJ(3);
t165 = t331 * t194;
t197 = cos(pkin(11));
t166 = t331 * t197;
t201 = sin(qJ(4));
t338 = cos(qJ(4));
t106 = -t201 * t165 + t166 * t338;
t154 = t194 * t338 + t201 * t197;
t196 = sin(pkin(6));
t204 = cos(qJ(2));
t290 = t196 * t204;
t210 = t154 * t290;
t400 = -qJD(1) * t210 + qJD(3) * t154 + qJD(4) * t106;
t388 = -m(7) - m(6);
t399 = mrSges(5,2) - mrSges(6,3);
t385 = -mrSges(6,2) + mrSges(5,1);
t250 = qJD(4) * t338;
t278 = qJD(4) * t201;
t253 = t194 * t278;
t147 = -t197 * t250 + t253;
t398 = -t147 * pkin(5) + t400;
t148 = t154 * qJD(4);
t223 = qJ(5) * t147 - qJD(5) * t154;
t202 = sin(qJ(2));
t282 = qJD(1) * t196;
t257 = t202 * t282;
t349 = pkin(4) + pkin(9);
t397 = -t148 * t349 - t223 + t257;
t396 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t200 = sin(qJ(6));
t203 = cos(qJ(6));
t260 = t338 * t197;
t240 = qJD(2) * t260;
t295 = t194 * t201;
t145 = qJD(2) * t295 - t240;
t110 = -qJD(4) * t200 + t145 * t203;
t247 = qJDD(2) * t338;
t268 = qJDD(2) * t201;
t95 = qJD(2) * t148 + t194 * t268 - t197 * t247;
t40 = qJD(6) * t110 + qJDD(4) * t203 + t200 * t95;
t94 = qJD(2) * t253 - qJD(4) * t240 - t194 * t247 - t197 * t268;
t87 = qJDD(6) - t94;
t22 = mrSges(7,1) * t87 - mrSges(7,3) * t40;
t146 = t154 * qJD(2);
t138 = qJD(6) + t146;
t59 = -mrSges(7,2) * t138 + mrSges(7,3) * t110;
t111 = qJD(4) * t203 + t145 * t200;
t60 = mrSges(7,1) * t138 - mrSges(7,3) * t111;
t228 = -t200 * t60 + t203 * t59;
t41 = -qJD(6) * t111 - qJDD(4) * t200 + t203 * t95;
t23 = -mrSges(7,2) * t87 + mrSges(7,3) * t41;
t395 = t228 * qJD(6) + t200 * t23 + t203 * t22;
t234 = mrSges(7,1) * t200 + mrSges(7,2) * t203;
t379 = mrSges(4,3) * (t194 ^ 2 + t197 ^ 2);
t383 = -Ifges(5,6) + Ifges(6,5);
t389 = m(7) * pkin(9);
t392 = -pkin(4) * t388 + t385 + t389;
t162 = qJD(2) * qJ(3) + t257;
t198 = cos(pkin(6));
t281 = qJD(1) * t198;
t177 = t197 * t281;
t322 = pkin(8) * qJD(2);
t101 = t177 + (-t162 - t322) * t194;
t114 = t197 * t162 + t194 * t281;
t102 = t197 * t322 + t114;
t307 = t102 * t201;
t45 = -t338 * t101 + t307;
t219 = pkin(5) * t146 + t45;
t390 = qJD(5) + t219;
t352 = t40 / 0.2e1;
t351 = t41 / 0.2e1;
t350 = t87 / 0.2e1;
t153 = -t260 + t295;
t187 = pkin(3) * t197 + pkin(2);
t220 = -qJ(5) * t154 - t187;
t55 = t153 * t349 + t220;
t105 = t338 * t165 + t166 * t201;
t65 = pkin(5) * t154 + t105;
t25 = t200 * t65 + t203 * t55;
t387 = -qJD(6) * t25 + t200 * t397 + t203 * t398;
t24 = -t200 * t55 + t203 * t65;
t386 = qJD(6) * t24 + t200 * t398 - t203 * t397;
t337 = g(3) * t196;
t384 = Ifges(6,4) - Ifges(5,5);
t18 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t73 = mrSges(6,1) * t95 - qJDD(4) * mrSges(6,3);
t382 = -t73 + t18;
t74 = -t94 * mrSges(6,1) + qJDD(4) * mrSges(6,2);
t381 = -qJDD(4) * mrSges(5,1) - mrSges(5,3) * t94 + t74;
t135 = Ifges(5,4) * t145;
t377 = t138 * Ifges(7,3);
t378 = t110 * Ifges(7,6);
t380 = t146 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t111 * Ifges(7,5) - t135 + t377 + t378;
t312 = qJDD(4) / 0.2e1;
t330 = mrSges(6,1) * t145;
t125 = -qJD(4) * mrSges(6,3) + t330;
t51 = -mrSges(7,1) * t110 + mrSges(7,2) * t111;
t311 = t125 - t51;
t327 = mrSges(5,3) * t146;
t329 = mrSges(6,1) * t146;
t285 = -qJD(4) * t385 + t327 + t329;
t289 = t200 * t204;
t172 = t196 * t289;
t291 = t196 * t202;
t139 = -t194 * t291 + t197 * t198;
t140 = t194 * t198 + t197 * t291;
t69 = -t139 * t338 + t140 * t201;
t52 = t203 * t69 + t172;
t192 = pkin(11) + qJ(4);
t188 = sin(t192);
t189 = cos(t192);
t376 = t188 * t399 - t189 * t385;
t249 = qJD(2) * t282;
t170 = t204 * t249;
t273 = qJDD(1) * t196;
t137 = t202 * t273 + t170;
t112 = t137 + t396;
t272 = qJDD(1) * t198;
t175 = t197 * t272;
t80 = -t112 * t194 + t175;
t81 = t197 * t112 + t194 * t272;
t374 = -t194 * t80 + t197 * t81;
t373 = m(5) - t388;
t237 = -mrSges(4,1) * t197 + mrSges(4,2) * t194;
t84 = -mrSges(6,2) * t145 - mrSges(6,3) * t146;
t371 = mrSges(5,1) * t145 + mrSges(5,2) * t146 + t237 * qJD(2) + t84;
t370 = 0.2e1 * t312;
t235 = mrSges(7,1) * t203 - mrSges(7,2) * t200;
t335 = t145 * pkin(5);
t46 = t201 * t101 + t102 * t338;
t39 = -qJD(4) * qJ(5) - t46;
t27 = -t39 - t335;
t325 = Ifges(7,4) * t111;
t33 = Ifges(7,2) * t110 + Ifges(7,6) * t138 + t325;
t368 = t235 * t27 - t203 * t33 / 0.2e1;
t367 = t201 * (qJD(3) * t194 + qJD(4) * t166) - qJD(3) * t260 + t165 * t250;
t63 = t175 + (-pkin(8) * qJDD(2) - t112) * t194;
t269 = qJDD(2) * t197;
t64 = pkin(8) * t269 + t81;
t266 = t101 * t250 + t201 * t63 + t338 * t64;
t11 = -qJDD(4) * qJ(5) + qJD(4) * (-qJD(5) + t307) - t266;
t328 = mrSges(5,3) * t145;
t123 = -qJD(4) * mrSges(5,2) - t328;
t365 = -m(6) * t39 + t123 - t125;
t38 = -qJD(4) * pkin(4) + qJD(5) + t45;
t364 = -m(6) * t38 - t285;
t169 = t202 * t249;
t136 = t204 * t273 - t169;
t217 = qJDD(3) - t136;
t103 = -qJDD(2) * t187 + t217;
t207 = qJ(5) * t94 - qJD(5) * t146 + t103;
t15 = t349 * t95 + t207;
t14 = -t101 * t278 - t102 * t250 - t201 * t64 + t338 * t63;
t211 = qJDD(5) - t14;
t5 = -t94 * pkin(5) - qJDD(4) * t349 + t211;
t26 = -qJD(4) * t349 + t390;
t256 = t204 * t282;
t229 = qJD(3) - t256;
t133 = -qJD(2) * t187 + t229;
t208 = -qJ(5) * t146 + t133;
t37 = t145 * t349 + t208;
t9 = -t200 * t37 + t203 * t26;
t1 = qJD(6) * t9 + t15 * t203 + t200 * t5;
t10 = t200 * t26 + t203 * t37;
t2 = -qJD(6) * t10 - t15 * t200 + t203 * t5;
t363 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t49 = pkin(4) * t145 + t208;
t362 = -t133 * mrSges(5,1) + t49 * mrSges(6,2);
t361 = qJ(5) * t388 - t234 + t399;
t360 = -m(5) * t46 - t365;
t128 = t188 * t291 - t198 * t189;
t310 = cos(pkin(10));
t245 = t310 * t202;
t195 = sin(pkin(10));
t292 = t195 * t204;
t142 = t198 * t245 + t292;
t246 = t196 * t310;
t96 = t142 * t188 + t189 * t246;
t244 = t310 * t204;
t293 = t195 * t202;
t144 = -t198 * t293 + t244;
t294 = t195 * t196;
t98 = t144 * t188 - t189 * t294;
t358 = g(1) * t98 + g(2) * t96 + g(3) * t128;
t357 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(5,3) - t235;
t213 = m(4) * pkin(2) - t237;
t356 = t188 * t234 + mrSges(3,1) + t213 - t376;
t355 = t9 * mrSges(7,1) + t133 * mrSges(5,2) - t10 * mrSges(7,2) - t49 * mrSges(6,3);
t251 = m(4) * qJ(3) + mrSges(4,3);
t354 = mrSges(3,2) - t251 + t357;
t353 = Ifges(7,1) * t352 + Ifges(7,4) * t351 + Ifges(7,5) * t350;
t348 = -t110 / 0.2e1;
t347 = -t111 / 0.2e1;
t346 = t111 / 0.2e1;
t345 = -t138 / 0.2e1;
t344 = -t145 / 0.2e1;
t343 = t145 / 0.2e1;
t342 = -t146 / 0.2e1;
t341 = t146 / 0.2e1;
t336 = t1 * t200;
t334 = -qJD(4) / 0.2e1;
t333 = qJD(4) / 0.2e1;
t326 = mrSges(7,3) * t203;
t324 = Ifges(7,4) * t200;
t323 = Ifges(7,4) * t203;
t321 = t146 * Ifges(5,4);
t320 = t146 * Ifges(6,6);
t309 = qJ(5) * t145;
t308 = qJ(5) * t188;
t306 = t114 * t197;
t141 = -t198 * t244 + t293;
t305 = t141 * t189;
t143 = t198 * t292 + t245;
t304 = t143 * t189;
t303 = t146 * t200;
t302 = t148 * t200;
t301 = t148 * t203;
t300 = t153 * t200;
t299 = t153 * t203;
t296 = t189 * t204;
t288 = t203 * t204;
t287 = -t141 * t187 + t142 * t331;
t286 = -t143 * t187 + t144 * t331;
t280 = qJD(2) * t202;
t279 = qJD(2) * t204;
t277 = qJD(6) * t200;
t276 = qJD(6) * t203;
t270 = qJDD(2) * t194;
t267 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t87;
t265 = t123 - t311;
t264 = m(4) + t373;
t263 = t196 * t288;
t255 = t194 * t279;
t254 = t196 * t280;
t248 = -t277 / 0.2e1;
t242 = -pkin(4) * t305 - t141 * t308 + t287;
t241 = -pkin(4) * t304 - t143 * t308 + t286;
t152 = -mrSges(4,1) * t269 + mrSges(4,2) * t270;
t239 = t10 * t203 - t200 * t9;
t232 = Ifges(7,1) * t200 + t323;
t231 = Ifges(7,2) * t203 + t324;
t230 = Ifges(7,5) * t200 + Ifges(7,6) * t203;
t227 = -t200 * t59 - t203 * t60;
t226 = t260 * t290;
t113 = -t162 * t194 + t177;
t224 = -t113 * t194 + t306;
t218 = -t200 * t69 + t263;
t70 = t201 * t139 + t140 * t338;
t216 = t153 * t276 + t302;
t215 = t153 * t277 - t301;
t209 = qJD(6) * t239 + t2 * t203 + t336;
t206 = qJD(2) ^ 2;
t160 = t187 * t290;
t158 = -qJD(2) * pkin(2) + t229;
t134 = Ifges(6,6) * t145;
t121 = -qJDD(2) * pkin(2) + t217;
t109 = Ifges(7,4) * t110;
t88 = pkin(4) * t153 + t220;
t86 = t94 * mrSges(5,2);
t85 = t94 * mrSges(6,3);
t82 = pkin(4) * t146 + t309;
t77 = -t145 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t321;
t76 = Ifges(6,4) * qJD(4) - t146 * Ifges(6,2) + t134;
t75 = Ifges(6,5) * qJD(4) + t145 * Ifges(6,3) - t320;
t72 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t95;
t66 = -t153 * pkin(5) + t106;
t54 = pkin(4) * t148 + t223;
t50 = t146 * t349 + t309;
t47 = -pkin(5) * t148 - t367;
t43 = qJD(2) * t210 + qJD(4) * t70;
t42 = -t139 * t250 - qJD(2) * t226 + (qJD(4) * t140 + t196 * t255) * t201;
t34 = t111 * Ifges(7,1) + t138 * Ifges(7,5) + t109;
t31 = t95 * mrSges(5,1) - t86;
t30 = -t95 * mrSges(6,2) + t85;
t29 = t46 - t335;
t21 = pkin(4) * t95 + t207;
t20 = qJD(6) * t52 + t200 * t43 + t203 * t254;
t19 = qJD(6) * t218 - t200 * t254 + t203 * t43;
t17 = t200 * t29 + t203 * t50;
t16 = -t200 * t50 + t203 * t29;
t13 = -t102 * t278 + t266;
t12 = -qJDD(4) * pkin(4) + t211;
t7 = t40 * Ifges(7,4) + t41 * Ifges(7,2) + t87 * Ifges(7,6);
t6 = -pkin(5) * t95 - t11;
t3 = [t19 * t60 + t20 * t59 + t52 * t22 - t218 * t23 + t381 * t69 + t285 * t43 + (-t139 * t194 + t140 * t197) * qJDD(2) * mrSges(4,3) + (t72 + t382) * t70 - t265 * t42 + (-m(2) - m(3) - t264) * g(3) + m(5) * (t13 * t70 - t14 * t69 - t42 * t46 + t43 * t45) + m(6) * (-t11 * t70 + t12 * t69 + t38 * t43 + t39 * t42) + m(4) * (t139 * t80 + t140 * t81) + m(7) * (-t1 * t218 + t10 * t20 + t19 * t9 + t2 * t52 - t27 * t42 + t6 * t70) + ((-t206 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + qJD(2) * t371) * t202 + (qJDD(2) * mrSges(3,1) - t152 - t30 - t31 + (-mrSges(3,2) + t379) * t206) * t204 + m(5) * (-t103 * t204 + t133 * t280) + m(6) * (-t204 * t21 + t280 * t49) + m(4) * (-t113 * t255 - t121 * t204 + t158 * t280 + t279 * t306) + m(3) * (t136 * t204 + t137 * t202)) * t196 + (m(3) * t198 ^ 2 + m(2)) * qJDD(1); m(6) * (t21 * t88 + t49 * t54) + t138 * (Ifges(7,5) * t216 - Ifges(7,6) * t215) / 0.2e1 + (-t204 * t337 + t136 + t169) * mrSges(3,1) + (t202 * t337 - t137 + t170) * mrSges(3,2) + (Ifges(7,1) * t216 - Ifges(7,4) * t215) * t346 - (t202 * t251 + t204 * t213) * t337 + t33 * t301 / 0.2e1 + t34 * t302 / 0.2e1 + (-m(5) * t103 - t31) * t187 + t400 * (m(5) * t45 - t364) + (-m(5) * t160 + t388 * (t160 + (pkin(4) * t189 + t308) * t290) + (-t296 * t389 + t376 * t204 + (-mrSges(7,1) * t289 - mrSges(7,2) * t288) * t188 + (-t331 * t373 + t357) * t202) * t196) * g(3) + (Ifges(4,4) * t194 + Ifges(4,2) * t197) * t269 + (Ifges(4,1) * t194 + Ifges(4,4) * t197) * t270 + t110 * (Ifges(7,4) * t216 - Ifges(7,2) * t215) / 0.2e1 + (-t377 / 0.2e1 - t378 / 0.2e1 + t76 / 0.2e1 - t38 * mrSges(6,1) - t45 * mrSges(5,3) - Ifges(5,1) * t341 + Ifges(6,2) * t342 + Ifges(6,6) * t343 - Ifges(5,4) * t344 - Ifges(7,5) * t346 + t384 * t333 - t355 - t380 / 0.2e1) * t147 + (t12 * mrSges(6,1) + t103 * mrSges(5,2) - t14 * mrSges(5,3) - t21 * mrSges(6,3) + Ifges(5,5) * t312 + Ifges(7,5) * t352 + Ifges(7,6) * t351 + Ifges(7,3) * t350 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t95 - Ifges(6,2) * t94 - t370 * Ifges(6,4) + t363) * t154 + (-Ifges(5,4) * t95 + Ifges(5,5) * qJDD(4) + t267) * t154 / 0.2e1 + t27 * (mrSges(7,1) * t215 + mrSges(7,2) * t216) + (-t170 + t396) * t379 + t121 * t237 + Ifges(3,3) * qJDD(2) + (qJD(6) * t34 + t7) * t299 / 0.2e1 + t386 * t59 + t387 * t60 + (t1 * t25 + t10 * t386 + t2 * t24 + t27 * t47 + t387 * t9 + t6 * t66) * m(7) + (t1 * t299 - t10 * t215 - t2 * t300 - t216 * t9 - t296 * t337) * mrSges(7,3) + (-t77 / 0.2e1 + t75 / 0.2e1 + t39 * mrSges(6,1) - t46 * mrSges(5,3) - Ifges(5,4) * t341 + Ifges(6,6) * t342 + Ifges(6,3) * t343 - Ifges(5,2) * t344 + t383 * t333 - t362) * t148 + (-m(5) * t14 + m(6) * t12 + t381) * t105 + (-pkin(2) * t121 + t224 * qJD(3) + t374 * qJ(3) - (t158 * t202 + t204 * t224) * t282) * m(4) + t374 * mrSges(4,3) + (-m(5) * t133 - m(6) * t49 - t371) * t257 + t360 * t367 + (-m(7) * t27 + t360 - t51) * (qJD(1) * t226 - t256 * t295) + (-m(6) * t242 - m(7) * (-pkin(9) * t305 + t242) + mrSges(7,3) * t305 - m(5) * t287 + t354 * t142 + t356 * t141) * g(2) + (-m(7) * (-pkin(9) * t304 + t241) + mrSges(7,3) * t304 - m(6) * t241 - m(5) * t286 + t354 * t144 + t356 * t143) * g(1) - pkin(2) * t152 + (t103 * mrSges(5,1) + t11 * mrSges(6,1) - t21 * mrSges(6,2) - t13 * mrSges(5,3) + t230 * t350 + t231 * t351 + t232 * t352 - t6 * t235 + t33 * t248 + (Ifges(5,2) + Ifges(6,3)) * t95 + (Ifges(5,4) + Ifges(6,6)) * t94 + t383 * t370) * t153 + t54 * t84 + t88 * t30 + t66 * t18 + t47 * t51 + t24 * t22 + t25 * t23 - t94 * Ifges(5,1) * t154 + t300 * t353 + (m(5) * t13 - m(6) * t11 + t72 - t73) * t106; -t206 * t379 + t385 * t95 + t265 * t145 + t152 + t85 - t86 + t203 * t23 - t200 * t22 + (t227 - t285) * t146 + t227 * qJD(6) + (t1 * t203 + t145 * t27 - t2 * t200 - t138 * (t10 * t200 + t203 * t9)) * m(7) + (-t145 * t39 - t146 * t38 + t21) * m(6) + (t145 * t46 - t146 * t45 + t103) * m(5) + (-qJD(2) * t224 + t121) * m(4) + (-g(1) * t143 - g(2) * t141 + g(3) * t290) * t264; -t2 * t326 + (-t303 / 0.2e1 + t248) * t34 + (-pkin(4) * t12 - qJ(5) * t11 - qJD(5) * t39 - t49 * t82) * m(6) + (qJ(5) * t6 - t10 * t17 - t16 * t9 + t27 * t390) * m(7) - t39 * t329 + t219 * t51 + (Ifges(6,1) + Ifges(5,3)) * qJDD(4) + (-m(7) * t209 - t395) * t349 + t6 * t234 - (t110 * t231 + t111 * t232 + t138 * t230) * qJD(6) / 0.2e1 + (-t321 + t75) * t342 + (t320 + t77) * t341 + t383 * t95 + (Ifges(6,3) * t344 - t10 * t326 + t230 * t345 + t231 * t348 + t232 * t347 + t334 * t383 + t362 + t368) * t146 + t384 * t94 + (-Ifges(5,1) * t342 - Ifges(7,5) * t347 + Ifges(6,2) * t341 - Ifges(7,6) * t348 - Ifges(7,3) * t345 + t334 * t384 + t355) * t145 + (-Ifges(5,2) * t146 - t135 + t380) * t343 + t382 * qJ(5) - t311 * qJD(5) + t368 * qJD(6) + (t327 + t364) * t46 + (t328 + t365) * t45 - t200 * t7 / 0.2e1 + (t134 + t76) * t344 + (-t10 * t276 - t336 + (t277 + t303) * t9 + t358) * mrSges(7,3) - t82 * t84 - pkin(4) * t74 + (t361 * (t142 * t189 - t188 * t246) + t392 * t96) * g(2) + (t361 * (t188 * t198 + t189 * t291) + t392 * t128) * g(3) - t17 * t59 - t16 * t60 + (t361 * (t144 * t189 + t188 * t294) + t392 * t98) * g(1) + t38 * t330 + (Ifges(7,5) * t203 - Ifges(7,6) * t200) * t350 + (-Ifges(7,2) * t200 + t323) * t351 + (Ifges(7,1) * t203 - t324) * t352 + t203 * t353 - t11 * mrSges(6,3) + t12 * mrSges(6,2) - t13 * mrSges(5,2) + t14 * mrSges(5,1); t311 * qJD(4) + (t228 + t84) * t146 + t74 + (-qJD(4) * t27 + t146 * t239 + t209 - t358) * m(7) + (qJD(4) * t39 + t146 * t49 + t12 - t358) * m(6) + t395; -t27 * (mrSges(7,1) * t111 + mrSges(7,2) * t110) + (Ifges(7,1) * t110 - t325) * t347 + t33 * t346 + (Ifges(7,5) * t110 - Ifges(7,6) * t111) * t345 - t9 * t59 + t10 * t60 - g(1) * ((-t143 * t200 + t203 * t98) * mrSges(7,1) + (-t143 * t203 - t200 * t98) * mrSges(7,2)) - g(2) * ((-t141 * t200 + t203 * t96) * mrSges(7,1) + (-t141 * t203 - t200 * t96) * mrSges(7,2)) - g(3) * ((t128 * t203 + t172) * mrSges(7,1) + (-t128 * t200 + t263) * mrSges(7,2)) + (t10 * t111 + t110 * t9) * mrSges(7,3) + t267 + (-Ifges(7,2) * t111 + t109 + t34) * t348 + t363;];
tau  = t3;
