% Calculate vector of inverse dynamics joint torques for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:56
% EndTime: 2022-01-23 09:25:17
% DurationCPUTime: 11.00s
% Computational Cost: add. (4906->502), mult. (12058->688), div. (0->0), fcn. (8646->14), ass. (0->264)
t216 = cos(pkin(8));
t277 = qJD(1) * t216;
t189 = qJD(3) - t277;
t270 = qJDD(1) * qJ(2);
t272 = qJD(1) * qJD(2);
t184 = t270 + t272;
t360 = t184 + t272;
t359 = Ifges(4,3) + Ifges(5,3);
t213 = sin(pkin(9));
t215 = cos(pkin(9));
t218 = sin(qJ(3));
t221 = cos(qJ(3));
t167 = t213 * t221 + t215 * t218;
t154 = t167 * qJD(3);
t230 = qJD(1) * t167;
t338 = t216 * t230 - t154;
t234 = t213 * t218 - t215 * t221;
t337 = t189 * t234;
t212 = qJ(3) + pkin(9);
t201 = sin(t212);
t311 = pkin(3) * t218;
t171 = pkin(4) * t201 + t311;
t353 = -m(3) - m(6);
t358 = (-m(4) + t353) * qJ(2) - m(5) * (qJ(2) + t311) - m(6) * t171 + mrSges(2,2) - mrSges(3,3);
t214 = sin(pkin(8));
t210 = t214 ^ 2;
t211 = t216 ^ 2;
t357 = mrSges(3,3) * (t211 + t210);
t356 = t360 * t210;
t175 = pkin(2) * t216 + pkin(6) * t214 + pkin(1);
t207 = t221 * pkin(3);
t240 = mrSges(3,1) * t216 - mrSges(3,2) * t214;
t307 = qJ(4) + pkin(6);
t354 = m(4) * t175 + m(5) * (pkin(1) + (pkin(2) + t207) * t216) + mrSges(2,1) + t240 + (m(5) * t307 + mrSges(4,3) + mrSges(5,3) + mrSges(6,3)) * t214;
t268 = qJDD(1) * t216;
t187 = qJDD(3) - t268;
t316 = t187 / 0.2e1;
t271 = qJD(1) * qJD(3);
t148 = (qJDD(1) * t221 - t218 * t271) * t214;
t149 = (-qJDD(1) * t218 - t221 * t271) * t214;
t82 = t148 * t215 + t149 * t213;
t352 = Ifges(5,5) * t82;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t278 = qJD(1) * t214;
t255 = t221 * t278;
t256 = t218 * t278;
t123 = t213 * t256 - t215 * t255;
t124 = t214 * t230;
t247 = t123 * t217 - t220 * t124;
t81 = -t148 * t213 + t149 * t215;
t21 = qJD(5) * t247 + t217 * t81 + t220 * t82;
t351 = Ifges(6,5) * t21;
t350 = Ifges(5,6) * t81;
t64 = t123 * t220 + t124 * t217;
t22 = qJD(5) * t64 - t217 * t82 + t220 * t81;
t349 = Ifges(6,6) * t22;
t94 = -t167 * t217 - t220 * t234;
t347 = qJD(5) * t94 + t217 * t338 - t220 * t337;
t95 = t167 * t220 - t217 * t234;
t346 = -qJD(5) * t95 + t217 * t337 + t220 * t338;
t345 = Ifges(4,5) * t148;
t344 = Ifges(4,6) * t149;
t179 = qJDD(5) + t187;
t343 = Ifges(6,3) * t179;
t312 = pkin(3) * t215;
t198 = pkin(4) + t312;
t313 = pkin(3) * t213;
t151 = t198 * t217 + t220 * t313;
t309 = pkin(7) * t124;
t152 = -qJD(1) * t175 + qJD(2);
t259 = qJ(2) * t277;
t105 = t152 * t218 + t221 * t259;
t89 = -qJ(4) * t256 + t105;
t299 = t215 * t89;
t132 = t221 * t152;
t293 = t214 * t221;
t262 = qJ(4) * t293;
t291 = t216 * t218;
t263 = qJ(2) * t291;
t229 = -t262 - t263;
t88 = qJD(1) * t229 + t132;
t45 = -t213 * t88 - t299;
t29 = t45 + t309;
t310 = pkin(7) * t123;
t83 = t213 * t89;
t46 = t215 * t88 - t83;
t30 = t46 + t310;
t342 = -t151 * qJD(5) + t217 * t30 - t220 * t29;
t150 = t198 * t220 - t217 * t313;
t341 = t150 * qJD(5) - t217 * t29 - t220 * t30;
t326 = m(5) * pkin(3);
t340 = t326 + mrSges(4,1);
t219 = sin(qJ(1));
t222 = cos(qJ(1));
t339 = t214 * (-g(1) * t222 - g(2) * t219);
t289 = t216 * t221;
t188 = qJ(2) * t289;
t122 = -t218 * t175 + t188;
t335 = t187 * t359 + t344 + t345 + t350 + t352;
t334 = 0.2e1 * t316;
t333 = t211 * t360;
t303 = Ifges(4,4) * t221;
t304 = Ifges(4,4) * t218;
t332 = (-qJ(2) * (mrSges(4,1) * t221 - mrSges(4,2) * t218) - t221 * (-Ifges(4,1) * t218 - t303) / 0.2e1 + t218 * (-Ifges(4,2) * t221 - t304) / 0.2e1) * t210;
t147 = -qJDD(1) * t175 + qJDD(2);
t131 = t221 * t147;
t258 = qJ(2) * qJD(3) * t216;
t273 = qJD(4) * t214;
t224 = qJD(1) * (-t258 - t273);
t36 = pkin(3) * t187 - qJ(4) * t148 + t131 + (-qJD(3) * t152 - t184 * t216) * t218 + t221 * t224;
t274 = qJD(3) * t221;
t261 = t218 * t147 + t152 * t274 + t184 * t289;
t41 = qJ(4) * t149 + t218 * t224 + t261;
t12 = -t213 * t41 + t215 * t36;
t6 = pkin(4) * t187 - pkin(7) * t82 + t12;
t13 = t213 * t36 + t215 * t41;
t7 = pkin(7) * t81 + t13;
t75 = pkin(3) * t189 + t88;
t42 = t215 * t75 - t83;
t25 = pkin(4) * t189 + t310 + t42;
t43 = t213 * t75 + t299;
t26 = t43 - t309;
t8 = -t217 * t26 + t220 * t25;
t2 = qJD(5) * t8 + t217 * t6 + t220 * t7;
t9 = t217 * t25 + t220 * t26;
t3 = -qJD(5) * t9 - t217 * t7 + t220 * t6;
t331 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t231 = (-t218 * Ifges(4,2) + t303) * t214;
t232 = (t221 * Ifges(4,1) - t304) * t214;
t330 = t221 * (t189 * Ifges(4,6) + qJD(1) * t231) + t218 * (t189 * Ifges(4,5) + qJD(1) * t232);
t243 = t218 * t258;
t52 = -qJD(1) * t243 + t261;
t53 = -qJD(3) * t105 - t184 * t291 + t131;
t327 = -t53 * mrSges(4,1) - t12 * mrSges(5,1) + t52 * mrSges(4,2) + t13 * mrSges(5,2);
t223 = qJD(1) ^ 2;
t325 = -t247 / 0.2e1;
t324 = t64 / 0.2e1;
t323 = -t64 / 0.2e1;
t322 = m(5) + m(6);
t321 = t8 * mrSges(6,3);
t320 = t9 * mrSges(6,3);
t319 = -t123 / 0.2e1;
t181 = qJD(5) + t189;
t317 = -t181 / 0.2e1;
t314 = Ifges(6,4) * t64;
t308 = g(3) * t214;
t276 = qJD(2) * t216;
t280 = -t175 * t274 + t221 * t276;
t72 = qJD(3) * t229 - t218 * t273 + t280;
t254 = t218 * t276;
t73 = -t254 - t221 * t273 + (-t188 + (qJ(4) * t214 + t175) * t218) * qJD(3);
t35 = t213 * t73 + t215 * t72;
t295 = t214 * t218;
t106 = -qJ(4) * t295 + t122;
t165 = t221 * t175;
t93 = -t262 - t165 + (-qJ(2) * t218 - pkin(3)) * t216;
t51 = t215 * t106 + t213 * t93;
t306 = mrSges(5,3) * t123;
t305 = mrSges(5,3) * t124;
t302 = Ifges(5,4) * t123;
t298 = qJ(2) * t223;
t203 = t214 * qJ(2);
t173 = t214 * t184;
t290 = t216 * t219;
t288 = t216 * t222;
t287 = t218 * t219;
t286 = t218 * t222;
t285 = t219 * t221;
t284 = t221 * t222;
t204 = qJ(5) + t212;
t195 = sin(t204);
t196 = cos(t204);
t117 = t195 * t290 + t196 * t222;
t118 = t195 * t222 - t196 * t290;
t283 = -t117 * mrSges(6,1) + t118 * mrSges(6,2);
t119 = -t195 * t288 + t196 * t219;
t120 = t195 * t219 + t196 * t288;
t282 = t119 * mrSges(6,1) - t120 * mrSges(6,2);
t281 = t356 * qJ(2);
t202 = cos(t212);
t172 = pkin(4) * t202 + t207;
t200 = t214 * qJD(2);
t253 = t214 * t274;
t162 = pkin(3) * t253 + t200;
t168 = pkin(3) * t295 + t203;
t275 = qJD(3) * t214;
t269 = qJDD(1) * t214;
t267 = t343 + t349 + t351;
t265 = mrSges(4,3) * t295;
t153 = pkin(3) * t256 + qJ(2) * t278 + qJD(4);
t252 = -t81 * mrSges(5,1) + t82 * mrSges(5,2);
t251 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t34 = -t213 * t72 + t215 * t73;
t50 = -t106 * t213 + t215 * t93;
t246 = pkin(3) * t255;
t245 = mrSges(4,3) * t256;
t244 = mrSges(4,3) * t255;
t100 = -pkin(3) * t149 + qJDD(4) + t173;
t241 = -mrSges(3,1) * t268 + mrSges(3,2) * t269;
t239 = mrSges(4,1) * t218 + mrSges(4,2) * t221;
t238 = -mrSges(6,1) * t195 - mrSges(6,2) * t196;
t143 = t234 * t214;
t40 = -pkin(4) * t216 + pkin(7) * t143 + t50;
t142 = t167 * t214;
t44 = -pkin(7) * t142 + t51;
t14 = -t217 * t44 + t220 * t40;
t15 = t217 * t40 + t220 * t44;
t76 = -t142 * t220 + t143 * t217;
t77 = -t142 * t217 - t143 * t220;
t144 = -mrSges(4,2) * t189 - t245;
t145 = mrSges(4,1) * t189 - t244;
t236 = t144 * t221 - t145 * t218;
t235 = (pkin(2) + t172) * t216 - (-pkin(7) - t307) * t214;
t233 = t267 - t331;
t158 = -t216 * t286 + t285;
t156 = t216 * t287 + t284;
t227 = t189 * t214 * (-Ifges(4,5) * t218 - Ifges(4,6) * t221);
t199 = -qJDD(1) * pkin(1) + qJDD(2);
t192 = t210 * t298;
t160 = t239 * t214;
t159 = t216 * t284 + t287;
t157 = -t216 * t285 + t286;
t146 = t239 * t278;
t138 = t201 * t219 + t202 * t288;
t137 = -t201 * t288 + t202 * t219;
t136 = t201 * t222 - t202 * t290;
t135 = t201 * t290 + t202 * t222;
t128 = t214 * t154;
t126 = t234 * t275;
t121 = -t165 - t263;
t116 = Ifges(5,4) * t124;
t108 = -mrSges(4,2) * t187 + mrSges(4,3) * t149;
t107 = mrSges(4,1) * t187 - mrSges(4,3) * t148;
t104 = -t218 * t259 + t132;
t103 = -qJD(3) * t122 - t254;
t102 = -t243 + t280;
t101 = -pkin(4) * t123 + t246;
t98 = mrSges(5,1) * t189 + t306;
t97 = -mrSges(5,2) * t189 - t305;
t96 = pkin(4) * t142 + t168;
t90 = -pkin(4) * t126 + t162;
t87 = pkin(4) * t124 + t153;
t74 = mrSges(5,1) * t124 - mrSges(5,2) * t123;
t60 = Ifges(6,4) * t247;
t59 = mrSges(5,1) * t187 - mrSges(5,3) * t82;
t58 = -mrSges(5,2) * t187 + mrSges(5,3) * t81;
t57 = -t123 * Ifges(5,1) + t189 * Ifges(5,5) - t116;
t56 = -t124 * Ifges(5,2) + t189 * Ifges(5,6) - t302;
t55 = mrSges(6,1) * t181 + mrSges(6,3) * t64;
t54 = -mrSges(6,2) * t181 + mrSges(6,3) * t247;
t47 = -pkin(4) * t81 + t100;
t39 = -qJD(5) * t77 + t126 * t220 + t128 * t217;
t38 = qJD(5) * t76 + t126 * t217 - t128 * t220;
t33 = -mrSges(6,1) * t247 - mrSges(6,2) * t64;
t28 = -Ifges(6,1) * t64 + Ifges(6,5) * t181 + t60;
t27 = Ifges(6,2) * t247 + Ifges(6,6) * t181 - t314;
t24 = pkin(7) * t126 + t35;
t23 = pkin(7) * t128 + t34;
t17 = -mrSges(6,2) * t179 + mrSges(6,3) * t22;
t16 = mrSges(6,1) * t179 - mrSges(6,3) * t21;
t5 = -qJD(5) * t15 - t217 * t24 + t220 * t23;
t4 = qJD(5) * t14 + t217 * t23 + t220 * t24;
t1 = [-t330 * t275 / 0.2e1 + (t227 / 0.2e1 + t104 * t265) * qJD(3) + t148 * t232 / 0.2e1 + (-t157 * mrSges(4,1) - t136 * mrSges(5,1) - t118 * mrSges(6,1) - t156 * mrSges(4,2) - t135 * mrSges(5,2) - t117 * mrSges(6,2) + t358 * t222 + (m(3) * pkin(1) - m(6) * (-pkin(1) - t235) + t354) * t219) * g(1) + (-t159 * mrSges(4,1) - t138 * mrSges(5,1) - t120 * mrSges(6,1) - t158 * mrSges(4,2) - t137 * mrSges(5,2) - t119 * mrSges(6,2) + (-m(6) * t235 + pkin(1) * t353 - t354) * t222 + t358 * t219) * g(2) + t39 * t320 + (Ifges(6,1) * t38 + Ifges(6,4) * t39) * t323 - pkin(1) * t241 + (-t105 * t253 - t293 * t53) * mrSges(4,3) + t96 * t251 + t168 * t252 + t247 * (Ifges(6,4) * t38 + Ifges(6,2) * t39) / 0.2e1 - t52 * t265 + ((Ifges(4,5) * t221 - Ifges(4,6) * t218) * t316 + Ifges(3,4) * t268 + Ifges(3,1) * t269) * t214 + m(4) * (t102 * t105 + t103 * t104 + t121 * t53 + t122 * t52 + t281) + (mrSges(6,2) * t47 - mrSges(6,3) * t3 + Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t179) * t77 - t199 * t240 - (-t100 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t82 + Ifges(5,2) * t81 + Ifges(5,6) * t334) * t142 - (t100 * mrSges(5,2) - t12 * mrSges(5,3) + Ifges(5,1) * t82 + Ifges(5,4) * t81 + Ifges(5,5) * t334) * t143 - (t267 + t335) * t216 / 0.2e1 - (Ifges(4,4) * t148 + Ifges(4,2) * t149 + Ifges(4,6) * t187) * t295 / 0.2e1 + t160 * t173 + (Ifges(4,1) * t148 + Ifges(4,4) * t149 + Ifges(4,5) * t187) * t293 / 0.2e1 + t149 * t231 / 0.2e1 + (t126 * t43 + t128 * t42) * mrSges(5,3) + (-Ifges(5,1) * t128 + Ifges(5,4) * t126) * t319 - t124 * (-Ifges(5,4) * t128 + Ifges(5,2) * t126) / 0.2e1 + t153 * (-mrSges(5,1) * t126 - mrSges(5,2) * t128) + t189 * (-Ifges(5,5) * t128 + Ifges(5,6) * t126) / 0.2e1 + (-t344 / 0.2e1 - t345 / 0.2e1 + Ifges(3,2) * t268 + Ifges(3,4) * t269 - t351 / 0.2e1 - t349 / 0.2e1 - t350 / 0.2e1 - t352 / 0.2e1 - t343 / 0.2e1 - t359 * t316 + t327 + t331) * t216 + (-mrSges(4,1) * t149 + mrSges(4,2) * t148) * t203 + t146 * t200 + (t333 + t356) * mrSges(3,3) + (-mrSges(6,1) * t47 + mrSges(6,3) * t2 + Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t179) * t76 - t332 * t271 + m(3) * (-pkin(1) * t199 + qJ(2) * t333 + t281) - t38 * t321 + m(5) * (t100 * t168 + t12 * t50 + t13 * t51 + t153 * t162 + t34 * t42 + t35 * t43) + m(6) * (t14 * t3 + t15 * t2 + t4 * t9 + t47 * t96 + t5 * t8 + t87 * t90) + t270 * t357 + t14 * t16 + t15 * t17 + t38 * t28 / 0.2e1 + t39 * t27 / 0.2e1 + t4 * t54 + t5 * t55 + t51 * t58 + t50 * t59 + t87 * (-mrSges(6,1) * t39 + mrSges(6,2) * t38) + t90 * t33 + t35 * t97 + t34 * t98 + t121 * t107 + t122 * t108 + t126 * t56 / 0.2e1 - t128 * t57 / 0.2e1 + t102 * t144 + t103 * t145 + t162 * t74 + t181 * (Ifges(6,5) * t38 + Ifges(6,6) * t39) / 0.2e1 + Ifges(2,3) * qJDD(1); -t337 * t97 + t346 * t55 + t347 * t54 + t241 + t236 * qJD(3) + (-t236 * t216 + (-t146 - t33 - t74) * t214) * qJD(1) - t223 * t357 + t338 * t98 + t94 * t16 + t95 * t17 - t234 * t59 + t167 * t58 + t218 * t108 + t221 * t107 + (-g(1) * t219 + g(2) * t222) * (m(3) + m(4) + t322) + (t2 * t95 - t87 * t278 + t3 * t94 + t346 * t8 + t347 * t9) * m(6) + (-t12 * t234 + t13 * t167 - t153 * t278 - t337 * t43 + t338 * t42) * m(5) + (t218 * t52 + t221 * t53 - t192 + t189 * (-t104 * t218 + t105 * t221)) * m(4) + (-t211 * t298 - t192 + t199) * m(3); -t74 * t246 - m(5) * (t153 * t246 + t42 * t45 + t43 * t46) + (Ifges(6,6) * t317 + Ifges(6,4) * t324 + Ifges(6,2) * t325 - t320 - t27 / 0.2e1 + t87 * mrSges(6,1)) * t64 + t330 * t278 / 0.2e1 + (-m(6) * (-t171 * t290 - t172 * t222) - t283 + t135 * mrSges(5,1) - t136 * mrSges(5,2) - mrSges(4,2) * t157 + t340 * t156) * g(2) - qJD(1) * t227 / 0.2e1 + (m(5) * t311 + mrSges(5,1) * t201 + mrSges(5,2) * t202 - t238) * t308 + t56 * t319 + (t12 * t215 + t13 * t213) * t326 + (Ifges(6,5) * t317 + t321 + Ifges(6,1) * t324 + Ifges(6,4) * t325 - t28 / 0.2e1 - t87 * mrSges(6,2)) * t247 + (Ifges(5,2) * t123 - t116 + t57) * t124 / 0.2e1 + t233 - t327 + t335 + (t244 + t145) * t105 + (-m(6) * (-t171 * t288 + t172 * t219) - t282 - t137 * mrSges(5,1) + t138 * mrSges(5,2) + mrSges(4,2) * t159 - t340 * t158) * g(1) + t341 * t54 + t342 * t55 + (-t101 * t87 + t150 * t3 + t151 * t2 + t171 * t308 + t341 * t9 + t342 * t8) * m(6) + (-t245 - t144) * t104 + t123 * (-Ifges(5,1) * t124 + t302) / 0.2e1 - t153 * (-mrSges(5,1) * t123 - mrSges(5,2) * t124) - t189 * (-Ifges(5,5) * t124 + Ifges(5,6) * t123) / 0.2e1 - t42 * t305 + t59 * t312 + t58 * t313 + t332 * t223 - t43 * t306 - t46 * t97 - t45 * t98 - t101 * t33 + t150 * t16 + t151 * t17 + g(3) * t160; t322 * t216 * g(3) - t123 * t98 + t124 * t97 - t247 * t54 - t64 * t55 + t251 + t252 + (-t247 * t9 - t64 * t8 + t339 + t47) * m(6) + (-t123 * t42 + t124 * t43 + t100 + t339) * m(5); -t87 * (-mrSges(6,1) * t64 + mrSges(6,2) * t247) + (Ifges(6,1) * t247 + t314) * t324 + t27 * t323 + (Ifges(6,5) * t247 + Ifges(6,6) * t64) * t317 - t8 * t54 + t9 * t55 - g(1) * t282 - g(2) * t283 - t238 * t308 + (t247 * t8 - t64 * t9) * mrSges(6,3) + t233 + (Ifges(6,2) * t64 + t28 + t60) * t325;];
tau = t1;
