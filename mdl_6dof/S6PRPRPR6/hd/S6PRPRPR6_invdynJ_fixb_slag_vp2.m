% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:42
% EndTime: 2019-03-08 19:47:07
% DurationCPUTime: 14.38s
% Computational Cost: add. (5002->595), mult. (10639->827), div. (0->0), fcn. (7632->14), ass. (0->285)
t391 = m(6) + m(7);
t406 = m(5) + t391;
t201 = sin(qJ(4));
t204 = cos(qJ(4));
t239 = pkin(4) * t204 + qJ(5) * t201;
t295 = qJD(5) * t204;
t129 = qJD(4) * t239 + qJD(3) - t295;
t194 = sin(pkin(11));
t197 = cos(pkin(11));
t205 = cos(qJ(2));
t206 = -pkin(2) - pkin(8);
t296 = qJD(4) * t206;
t269 = t204 * t296;
t196 = sin(pkin(6));
t303 = qJD(1) * t196;
t202 = sin(qJ(2));
t307 = t201 * t202;
t379 = t194 * t129 + t197 * t269 - (t194 * t205 + t197 * t307) * t303;
t405 = t197 * t129 - (-t194 * t307 + t197 * t205) * t303;
t189 = pkin(5) * t197 + pkin(4);
t193 = pkin(11) + qJ(6);
t191 = sin(t193);
t192 = cos(t193);
t248 = -t197 * mrSges(6,1) + t194 * mrSges(6,2);
t404 = -m(6) * pkin(4) - m(7) * t189 - t192 * mrSges(7,1) + t191 * mrSges(7,2) + t248;
t293 = t201 * qJD(2);
t188 = qJD(6) + t293;
t342 = -t188 / 0.2e1;
t300 = qJD(2) * t204;
t157 = qJD(4) * t197 - t194 * t300;
t158 = qJD(4) * t194 + t197 * t300;
t200 = sin(qJ(6));
t203 = cos(qJ(6));
t80 = t157 * t200 + t158 * t203;
t349 = -t80 / 0.2e1;
t252 = t203 * t157 - t158 * t200;
t351 = -t252 / 0.2e1;
t403 = Ifges(7,5) * t349 + Ifges(7,6) * t351 + Ifges(7,3) * t342;
t347 = -m(4) - m(5);
t388 = qJD(4) / 0.2e1;
t260 = -t194 * t206 + pkin(5);
t309 = t197 * t201;
t286 = pkin(9) * t309;
t402 = (t204 * t260 + t286) * qJD(4) + t405;
t298 = qJD(4) * t201;
t271 = t194 * t298;
t401 = -pkin(9) * t271 - t379;
t247 = t194 * mrSges(6,1) + t197 * mrSges(6,2);
t335 = mrSges(3,1) - mrSges(4,2);
t338 = pkin(5) * t194;
t400 = -m(7) * t338 - t191 * mrSges(7,1) - t192 * mrSges(7,2) - mrSges(5,3) - t247 - t335;
t333 = pkin(9) + qJ(5);
t399 = m(7) * t333 + mrSges(6,3) + mrSges(7,3);
t292 = qJD(2) * qJD(4);
t166 = qJDD(2) * t204 - t201 * t292;
t119 = qJDD(4) * t197 - t166 * t194;
t198 = cos(pkin(6));
t302 = qJD(1) * t198;
t184 = t201 * t302;
t301 = qJD(2) * t202;
t273 = t196 * t301;
t178 = qJD(1) * t273;
t310 = t196 * t205;
t132 = qJDD(1) * t310 - t178;
t232 = qJDD(3) - t132;
t104 = qJDD(2) * t206 + t232;
t276 = t205 * t303;
t238 = qJD(3) - t276;
t150 = qJD(2) * t206 + t238;
t291 = qJDD(1) * t198;
t297 = qJD(4) * t204;
t283 = t201 * t104 + t150 * t297 + t204 * t291;
t32 = qJDD(4) * qJ(5) + (qJD(5) - t184) * qJD(4) + t283;
t167 = qJDD(2) * t201 + t204 * t292;
t289 = qJDD(2) * qJ(3);
t290 = qJDD(1) * t202;
t225 = t196 * t290 + t289;
t237 = qJD(3) + t276;
t48 = pkin(4) * t167 - qJ(5) * t166 + (t237 - t295) * qJD(2) + t225;
t15 = t194 * t48 + t197 * t32;
t13 = pkin(9) * t119 + t15;
t120 = qJDD(4) * t194 + t166 * t197;
t14 = -t194 * t32 + t197 * t48;
t6 = pkin(5) * t167 - pkin(9) * t120 + t14;
t253 = -qJ(5) * t204 + qJ(3);
t168 = pkin(4) * t201 + t253;
t277 = t202 * t303;
t110 = qJD(2) * t168 + t277;
t275 = t204 * t302;
t101 = t150 * t201 + t275;
t87 = qJD(4) * qJ(5) + t101;
t41 = t197 * t110 - t194 * t87;
t26 = pkin(5) * t293 - pkin(9) * t158 + t41;
t42 = t194 * t110 + t197 * t87;
t31 = pkin(9) * t157 + t42;
t9 = -t200 * t31 + t203 * t26;
t1 = qJD(6) * t9 + t13 * t203 + t200 * t6;
t398 = t1 * mrSges(7,2);
t10 = t200 * t26 + t203 * t31;
t2 = -qJD(6) * t10 - t13 * t200 + t203 * t6;
t397 = t2 * mrSges(7,1);
t396 = t406 * pkin(8) + pkin(2) * (m(4) + t406) - t400;
t249 = mrSges(5,1) * t201 + mrSges(5,2) * t204;
t387 = mrSges(3,2) - mrSges(4,3);
t394 = t201 * t404 - t249 + t387;
t343 = t167 / 0.2e1;
t345 = t120 / 0.2e1;
t393 = Ifges(6,1) * t345 + Ifges(6,5) * t343;
t331 = Ifges(5,4) * t204;
t243 = -t201 * Ifges(5,2) + t331;
t392 = t10 * mrSges(7,2) + Ifges(5,6) * t388 + qJD(2) * t243 / 0.2e1 - t158 * Ifges(6,5) / 0.2e1 - t157 * Ifges(6,6) / 0.2e1 - Ifges(6,3) * t293 / 0.2e1 - t9 * mrSges(7,1) + t403;
t24 = qJD(6) * t252 + t119 * t200 + t120 * t203;
t356 = t24 / 0.2e1;
t25 = -qJD(6) * t80 + t119 * t203 - t120 * t200;
t355 = t25 / 0.2e1;
t346 = t119 / 0.2e1;
t159 = qJDD(6) + t167;
t344 = t159 / 0.2e1;
t155 = t197 * t168;
t308 = t197 * t204;
t75 = -pkin(9) * t308 + t201 * t260 + t155;
t306 = t201 * t206;
t109 = t194 * t168 + t197 * t306;
t315 = t194 * t204;
t86 = -pkin(9) * t315 + t109;
t33 = -t200 * t86 + t203 * t75;
t390 = qJD(6) * t33 + t200 * t402 - t203 * t401;
t34 = t200 * t75 + t203 * t86;
t389 = -qJD(6) * t34 + t200 * t401 + t203 * t402;
t161 = t194 * t203 + t197 * t200;
t100 = t150 * t204 - t184;
t163 = t239 * qJD(2);
t52 = -t100 * t194 + t197 * t163;
t43 = (pkin(5) * t204 + t286) * qJD(2) + t52;
t274 = t194 * t293;
t53 = t197 * t100 + t194 * t163;
t51 = pkin(9) * t274 + t53;
t171 = t333 * t194;
t172 = t333 * t197;
t98 = -t171 * t200 + t172 * t203;
t386 = -qJD(5) * t161 - qJD(6) * t98 + t200 * t51 - t203 * t43;
t234 = t194 * t200 - t197 * t203;
t97 = -t171 * t203 - t172 * t200;
t385 = -qJD(5) * t234 + qJD(6) * t97 - t200 * t43 - t203 * t51;
t56 = -t119 * mrSges(6,1) + t120 * mrSges(6,2);
t383 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t166 - t56;
t128 = t234 * t204;
t143 = t161 * qJD(2);
t145 = t161 * qJD(6);
t382 = -qJD(4) * t128 - t145 * t201 - t143;
t126 = t161 * t204;
t368 = qJD(6) * t234;
t381 = t234 * qJD(2) - qJD(4) * t126 + t201 * t368;
t380 = -t194 * t269 + t405;
t280 = mrSges(5,3) * t300;
t378 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t157 + mrSges(6,2) * t158 + t280;
t377 = t234 * t201;
t117 = t201 * t143;
t376 = t117 + t145;
t118 = qJD(2) * t377;
t375 = t118 + t368;
t329 = Ifges(6,4) * t197;
t242 = -Ifges(6,2) * t194 + t329;
t330 = Ifges(6,4) * t194;
t244 = Ifges(6,1) * t197 - t330;
t374 = t157 * (Ifges(6,6) * t204 - t201 * t242) + t158 * (Ifges(6,5) * t204 - t201 * t244);
t240 = Ifges(6,5) * t197 - Ifges(6,6) * t194;
t373 = t204 * (-Ifges(5,1) * t201 - t331) + t201 * (Ifges(6,3) * t204 - t201 * t240);
t372 = m(6) * qJ(5) + t399;
t115 = -mrSges(6,2) * t293 + mrSges(6,3) * t157;
t116 = mrSges(6,1) * t293 - mrSges(6,3) * t158;
t371 = t115 * t197 - t116 * t194;
t266 = qJD(4) * t302;
t37 = -t201 * t266 + t283;
t38 = -t150 * t298 - t201 * t291 + (t104 - t266) * t204;
t370 = t201 * t37 + t204 * t38;
t68 = -mrSges(6,2) * t167 + mrSges(6,3) * t119;
t69 = mrSges(6,1) * t167 - mrSges(6,3) * t120;
t369 = -t194 * t69 + t197 * t68;
t236 = -t14 * t194 + t15 * t197;
t365 = mrSges(5,1) - t404;
t364 = mrSges(5,2) - t372;
t165 = qJD(2) * qJ(3) + t277;
t83 = -qJD(4) * pkin(4) + qJD(5) - t100;
t323 = t201 * t83;
t363 = t247 * t323 - t42 * (mrSges(6,3) * t194 * t201 - mrSges(6,2) * t204) - t41 * (mrSges(6,1) * t204 + mrSges(6,3) * t309) - t165 * (mrSges(5,1) * t204 - mrSges(5,2) * t201);
t362 = -m(6) * t83 - t378;
t359 = -m(6) * t253 + t394 + t399 * t204 + (-m(7) + t347) * qJ(3);
t207 = qJD(2) ^ 2;
t358 = Ifges(7,4) * t356 + Ifges(7,2) * t355 + Ifges(7,6) * t344;
t357 = Ifges(7,1) * t356 + Ifges(7,4) * t355 + Ifges(7,5) * t344;
t339 = Ifges(7,4) * t80;
t28 = Ifges(7,2) * t252 + Ifges(7,6) * t188 + t339;
t354 = t28 / 0.2e1;
t74 = Ifges(7,4) * t252;
t29 = Ifges(7,1) * t80 + Ifges(7,5) * t188 + t74;
t353 = t29 / 0.2e1;
t352 = Ifges(6,4) * t346 + t393;
t350 = t252 / 0.2e1;
t348 = t80 / 0.2e1;
t341 = t188 / 0.2e1;
t332 = Ifges(5,4) * t201;
t36 = -qJDD(4) * pkin(4) + qJDD(5) - t38;
t322 = t204 * t36;
t320 = cos(pkin(10));
t316 = t165 * t205;
t195 = sin(pkin(10));
t314 = t195 * t196;
t313 = t195 * t202;
t312 = t195 * t205;
t311 = t196 * t202;
t304 = pkin(2) * t310 + qJ(3) * t311;
t299 = qJD(2) * t205;
t294 = qJDD(2) * mrSges(4,2);
t5 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t287 = t5 - t383;
t35 = -mrSges(7,1) * t252 + mrSges(7,2) * t80;
t285 = t35 + t378;
t284 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t159;
t282 = -t347 + t391;
t281 = mrSges(5,3) * t293;
t279 = t201 * t310;
t272 = t196 * t299;
t265 = qJD(1) * t299;
t261 = t293 / 0.2e1;
t259 = -t206 + t338;
t258 = t196 * t320;
t257 = t320 * t202;
t256 = t320 * t205;
t255 = -t292 / 0.2e1;
t250 = t196 * t265;
t245 = t204 * Ifges(5,1) - t332;
t241 = -Ifges(5,5) * t201 - Ifges(5,6) * t204;
t235 = -t194 * t41 + t197 * t42;
t147 = t198 * t204 - t279;
t92 = -t147 * t194 + t197 * t311;
t93 = t147 * t197 + t194 * t311;
t39 = -t200 * t93 + t203 * t92;
t40 = t200 * t92 + t203 * t93;
t105 = qJD(2) * t237 + t225;
t233 = qJ(3) * t105 + qJD(3) * t165;
t146 = t198 * t201 + t204 * t310;
t231 = t105 * t202 + t165 * t299;
t229 = t201 * (-Ifges(5,2) * t204 - t332);
t140 = t198 * t312 + t257;
t88 = -t140 * t204 + t201 * t314;
t138 = -t198 * t256 + t313;
t90 = t138 * t204 + t201 * t258;
t227 = -g(1) * t88 + g(2) * t90 - g(3) * t146;
t210 = (-t100 * t201 + t101 * t204) * qJD(4) + t370;
t174 = -qJD(4) * mrSges(5,2) - t281;
t164 = t249 * qJD(2);
t162 = -qJD(2) * pkin(2) + t238;
t156 = t259 * t204;
t152 = Ifges(5,5) * qJD(4) + qJD(2) * t245;
t141 = -t198 * t313 + t256;
t139 = t198 * t257 + t312;
t137 = t259 * t298;
t135 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t167;
t133 = (t265 + t290) * t196;
t125 = t161 * t201;
t112 = -qJDD(2) * pkin(2) + t232;
t108 = -t194 * t306 + t155;
t99 = mrSges(5,1) * t167 + mrSges(5,2) * t166;
t95 = -qJD(4) * t279 + t198 * t297 - t204 * t273;
t94 = -qJD(4) * t146 + t201 * t273;
t91 = -t138 * t201 + t204 * t258;
t89 = t140 * t201 + t204 * t314;
t73 = t275 + (-qJD(2) * t338 + t150) * t201;
t72 = t158 * Ifges(6,1) + t157 * Ifges(6,4) + Ifges(6,5) * t293;
t71 = t158 * Ifges(6,4) + t157 * Ifges(6,2) + Ifges(6,6) * t293;
t66 = t161 * t298 + t204 * t368;
t64 = qJD(4) * t377 - t145 * t204;
t62 = t194 * t272 + t197 * t94;
t61 = -t194 * t94 + t197 * t272;
t60 = mrSges(7,1) * t188 - mrSges(7,3) * t80;
t59 = -mrSges(7,2) * t188 + mrSges(7,3) * t252;
t55 = -pkin(5) * t157 + t83;
t46 = t120 * Ifges(6,4) + t119 * Ifges(6,2) + t167 * Ifges(6,6);
t20 = -pkin(5) * t119 + t36;
t19 = -mrSges(7,2) * t159 + mrSges(7,3) * t25;
t18 = mrSges(7,1) * t159 - mrSges(7,3) * t24;
t12 = -qJD(6) * t40 - t200 * t62 + t203 * t61;
t11 = qJD(6) * t39 + t200 * t61 + t203 * t62;
t3 = [t11 * t59 + t62 * t115 + t61 * t116 + t12 * t60 + t147 * t135 + t94 * t174 + t39 * t18 + t40 * t19 + t93 * t68 + t92 * t69 + t285 * t95 + t287 * t146 + (-m(2) - m(3) - t282) * g(3) + m(6) * (t14 * t92 + t146 * t36 + t15 * t93 + t41 * t61 + t42 * t62 + t83 * t95) + m(7) * (t1 * t40 + t10 * t11 + t12 * t9 + t146 * t20 + t2 * t39 + t55 * t95) + m(5) * (-t100 * t95 + t101 * t94 - t146 * t38 + t147 * t37) + (t164 * t299 + t202 * t99 + m(5) * t231 + m(3) * (t132 * t205 + t133 * t202) + m(4) * (-t112 * t205 + t162 * t301 + t231) + (-t202 * t335 - t205 * t387) * t207 + (-t202 * t387 + t205 * t335) * qJDD(2)) * t196 + (m(2) + 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1) * t198 ^ 2) * qJDD(1); (t132 + t178) * mrSges(3,1) + t389 * t60 + (t1 * t34 + t10 * t390 - t137 * t55 + t156 * t20 + t2 * t33 + t389 * t9) * m(7) + t390 * t59 + (-Ifges(7,4) * t128 - Ifges(7,2) * t126 + Ifges(7,6) * t201) * t355 + (-Ifges(7,1) * t128 - Ifges(7,4) * t126 + Ifges(7,5) * t201) * t356 + (-Ifges(7,5) * t128 - Ifges(7,6) * t126 + Ifges(7,3) * t201) * t344 + t20 * (mrSges(7,1) * t126 - mrSges(7,2) * t128) - (t197 * t72 + t152) * t298 / 0.2e1 + (Ifges(6,5) * t120 + Ifges(6,6) * t119 + Ifges(6,3) * t167 + t284) * t201 / 0.2e1 - pkin(2) * t294 + (qJD(2) * qJD(3) + t105 - t250 + t289) * mrSges(4,3) + (Ifges(7,4) * t64 + Ifges(7,2) * t66) * t350 + (Ifges(7,1) * t64 + Ifges(7,4) * t66) * t348 + (-t178 + t112) * mrSges(4,2) + t71 * t271 / 0.2e1 - t46 * t315 / 0.2e1 + t15 * (-mrSges(6,2) * t201 - mrSges(6,3) * t315) + (-(t316 + (t100 * t204 + t101 * t201) * t202) * t303 + t206 * t210 + t233) * m(5) + (-(t162 * t202 + t316) * t303 - pkin(2) * t112 + t233) * m(4) - t167 * t243 / 0.2e1 + t166 * t245 / 0.2e1 + t105 * t249 + (Ifges(7,5) * t64 + Ifges(7,6) * t66) * t341 + t14 * (mrSges(6,1) * t201 - mrSges(6,3) * t308) + t238 * t164 + (-t133 + t250) * mrSges(3,2) + (-m(4) * t304 - t406 * (pkin(8) * t310 + t304) + (t400 * t205 + (t204 * t372 + t394) * t202) * t196) * g(3) + (-t101 * mrSges(5,3) + Ifges(7,5) * t348 + Ifges(7,6) * t350 + Ifges(7,3) * t341 - t392) * t297 + (t100 * t298 - t370) * mrSges(5,3) + t204 * (Ifges(5,1) * t166 - Ifges(5,4) * t167) / 0.2e1 + t379 * t115 + t380 * t116 + (t108 * t14 + t109 * t15 + (t298 * t83 - t322) * t206 + t379 * t42 + t380 * t41) * m(6) - t201 * (Ifges(5,4) * t166 - Ifges(5,2) * t167) / 0.2e1 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (m(7) * t55 + t35 - t362) * t204 * t277 + (-t201 * t277 + t269) * t174 + qJDD(4) * (Ifges(5,5) * t204 - Ifges(5,6) * t201) + t64 * t353 + t66 * t354 - t128 * t357 - t126 * t358 + (t201 * Ifges(6,3) + t204 * t240) * t343 + (Ifges(6,5) * t201 + t204 * t244) * t345 + (Ifges(6,6) * t201 + t204 * t242) * t346 + t308 * t352 + t374 * t388 + t247 * t322 + t201 * t397 + t135 * t306 + t373 * t292 / 0.2e1 + t378 * t201 * t296 + t383 * t204 * t206 + (-t1 * t126 + t10 * t66 + t128 * t2 - t64 * t9) * mrSges(7,3) + (t138 * t396 + t139 * t359) * g(2) + (t140 * t396 + t141 * t359) * g(1) + t33 * t18 + t34 * t19 - t201 * t398 + t55 * (-mrSges(7,1) * t66 + mrSges(7,2) * t64) + qJ(3) * t99 + t108 * t69 + t109 * t68 - t137 * t35 + (t241 * t388 - t363) * qJD(4) + t156 * t5 + t229 * t255; t294 - t207 * mrSges(4,3) - t125 * t18 - t377 * t19 + t381 * t60 + t382 * t59 - t287 * t204 + (t135 + t369) * t201 + (-t194 * t115 - t197 * t116 + t165 * t347 - t164) * qJD(2) + ((t174 + t371) * t204 + t285 * t201) * qJD(4) + m(4) * t112 + m(5) * t210 + (-t1 * t377 + t10 * t382 - t125 * t2 - t20 * t204 + t55 * t298 + t381 * t9) * m(7) + (-t322 + t236 * t201 + (t204 * t235 + t323) * qJD(4) - (t42 * t194 + t41 * t197) * qJD(2)) * m(6) + (-g(1) * t140 - g(2) * t138 + g(3) * t310) * t282; (-pkin(4) * t36 + qJ(5) * t236 + qJD(5) * t235 - t41 * t52 - t42 * t53) * m(6) + (t280 + t362) * t101 + (t352 + t393) * t194 - t71 * t274 / 0.2e1 + t329 * t345 + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t351 + (t229 / 0.2e1 - t373 / 0.2e1) * t207 + (t364 * t89 + t365 * t88) * g(1) + (t146 * t365 + t147 * t364) * g(3) + (-t364 * t91 - t365 * t90) * g(2) + t330 * t346 + (t72 * t261 + Ifges(6,6) * t343 + Ifges(6,2) * t346 + t46 / 0.2e1) * t197 + (-t281 - t174) * t100 + t36 * t248 + t236 * mrSges(6,3) + t369 * qJ(5) + (-Ifges(7,5) * t368 - Ifges(7,6) * t145) * t341 + (-Ifges(7,1) * t368 - Ifges(7,4) * t145) * t348 + (-Ifges(7,4) * t368 - Ifges(7,2) * t145) * t350 - t368 * t353 + (Ifges(7,4) * t161 - Ifges(7,2) * t234) * t355 + (Ifges(7,1) * t161 - Ifges(7,4) * t234) * t356 + (Ifges(7,5) * t161 - Ifges(7,6) * t234) * t344 + t20 * (mrSges(7,1) * t234 + mrSges(7,2) * t161) + (-t1 * t234 - t10 * t376 - t161 * t2 + t375 * t9) * mrSges(7,3) - t234 * t358 + (Ifges(7,5) * t118 + Ifges(7,6) * t117) * t342 + (t392 + t403) * t300 + Ifges(5,3) * qJDD(4) - t145 * t354 + t161 * t357 + t371 * qJD(5) + (mrSges(7,1) * t376 - mrSges(7,2) * t375) * t55 - t189 * t5 + Ifges(5,5) * t166 - Ifges(5,6) * t167 - t37 * mrSges(5,2) + t38 * mrSges(5,1) - pkin(4) * t56 - t73 * t35 + t97 * t18 + t98 * t19 - t53 * t115 - t52 * t116 - t117 * t28 / 0.2e1 - t118 * t29 / 0.2e1 + (Ifges(7,1) * t118 + Ifges(7,4) * t117) * t349 + t385 * t59 + t386 * t60 + (t1 * t98 + t10 * t385 - t189 * t20 + t2 * t97 + t386 * t9 - t55 * t73) * m(7) + (t363 - t374 / 0.2e1) * qJD(2) + t241 * t255 + t152 * t261; -t157 * t115 + t158 * t116 - t252 * t59 + t80 * t60 + t5 + t56 + (-t10 * t252 + t80 * t9 + t20 + t227) * m(7) + (-t157 * t42 + t158 * t41 + t227 + t36) * m(6); -t398 + t397 - t55 * (mrSges(7,1) * t80 + mrSges(7,2) * t252) + (Ifges(7,1) * t252 - t339) * t349 + t28 * t348 + (Ifges(7,5) * t252 - Ifges(7,6) * t80) * t342 - t9 * t59 + t10 * t60 - g(1) * ((t141 * t192 - t191 * t89) * mrSges(7,1) + (-t141 * t191 - t192 * t89) * mrSges(7,2)) - g(2) * ((t139 * t192 + t191 * t91) * mrSges(7,1) + (-t139 * t191 + t192 * t91) * mrSges(7,2)) - g(3) * ((-t147 * t191 + t192 * t311) * mrSges(7,1) + (-t147 * t192 - t191 * t311) * mrSges(7,2)) + (t10 * t80 + t252 * t9) * mrSges(7,3) + t284 + (-Ifges(7,2) * t80 + t29 + t74) * t351;];
tau  = t3;
