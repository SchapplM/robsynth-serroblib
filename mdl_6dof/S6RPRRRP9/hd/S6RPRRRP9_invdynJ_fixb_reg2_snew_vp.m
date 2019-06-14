% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP9
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
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:53:26
% EndTime: 2019-05-06 01:53:36
% DurationCPUTime: 4.40s
% Computational Cost: add. (19153->367), mult. (37852->461), div. (0->0), fcn. (25488->8), ass. (0->260)
t222 = sin(qJ(4));
t226 = cos(qJ(4));
t227 = cos(qJ(3));
t277 = qJD(1) * t227;
t200 = qJD(3) * t226 - t222 * t277;
t201 = qJD(3) * t222 + t226 * t277;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t176 = -t225 * t200 + t201 * t221;
t178 = t200 * t221 + t201 * t225;
t137 = t178 * t176;
t271 = qJD(1) * qJD(3);
t259 = t227 * t271;
t223 = sin(qJ(3));
t269 = t223 * qJDD(1);
t203 = -t259 - t269;
t199 = qJDD(4) - t203;
t196 = qJDD(5) + t199;
t320 = t137 - t196;
t326 = t320 * pkin(5);
t215 = t227 * qJDD(1);
t260 = t223 * t271;
t204 = t215 - t260;
t255 = -t226 * qJDD(3) + t222 * t204;
t170 = -qJD(4) * t201 - t255;
t242 = -t222 * qJDD(3) - t226 * t204;
t171 = qJD(4) * t200 - t242;
t120 = -qJD(5) * t176 + t170 * t221 + t171 * t225;
t213 = qJD(1) * t223 + qJD(4);
t210 = qJD(5) + t213;
t158 = t210 * t176;
t99 = t120 + t158;
t323 = qJ(6) * t99;
t270 = qJD(2) * qJD(1);
t216 = 0.2e1 * t270;
t218 = qJDD(1) * qJ(2);
t224 = sin(qJ(1));
t228 = cos(qJ(1));
t248 = t228 * g(1) + t224 * g(2);
t240 = -t218 + t248;
t238 = t216 - t240;
t244 = -t204 + t260;
t245 = -t203 + t259;
t230 = qJD(1) ^ 2;
t315 = pkin(7) + pkin(1);
t319 = t315 * t230;
t141 = t245 * pkin(3) + t244 * pkin(8) + t238 - t319;
t257 = t224 * g(1) - t228 * g(2);
t247 = qJDD(2) - t257;
t279 = t230 * qJ(2);
t236 = t247 - t279;
t190 = -t315 * qJDD(1) + t236;
t180 = t227 * g(3) - t223 * t190;
t229 = qJD(3) ^ 2;
t250 = pkin(3) * t223 - pkin(8) * t227;
t239 = t230 * t250;
t153 = -t229 * pkin(3) + qJDD(3) * pkin(8) - t223 * t239 - t180;
t113 = -t226 * t141 + t222 * t153;
t191 = t213 * t200;
t146 = t171 - t191;
t183 = t200 * t201;
t316 = t183 + t199;
t78 = pkin(4) * t316 - pkin(9) * t146 - t113;
t114 = t222 * t141 + t226 * t153;
t197 = t200 ^ 2;
t249 = pkin(4) * t213 - pkin(9) * t201;
t81 = -t197 * pkin(4) + t170 * pkin(9) - t213 * t249 + t114;
t46 = t221 * t81 - t225 * t78;
t234 = 0.2e1 * qJD(6) * t178 + t323 + t326 + t46;
t233 = -t234 - t326;
t174 = t176 ^ 2;
t209 = t210 ^ 2;
t132 = -t209 - t174;
t292 = t320 * t225;
t86 = t132 * t221 - t292;
t85 = pkin(4) * t86;
t325 = t233 + t85;
t175 = t178 ^ 2;
t151 = -t175 - t209;
t127 = t137 + t196;
t295 = t127 * t221;
t104 = t151 * t225 - t295;
t103 = pkin(4) * t104;
t256 = -t225 * t170 + t221 * t171;
t119 = -qJD(5) * t178 - t256;
t154 = pkin(5) * t210 - qJ(6) * t178;
t47 = t221 * t78 + t225 * t81;
t30 = -t174 * pkin(5) + t119 * qJ(6) - 0.2e1 * qJD(6) * t176 - t210 * t154 + t47;
t237 = pkin(5) * t151 - t30;
t324 = t103 + t237;
t322 = t222 * t316;
t321 = t226 * t316;
t293 = t320 * t221;
t318 = t120 - t158;
t96 = (qJD(5) - t210) * t178 + t256;
t142 = (qJD(4) - t213) * t201 + t255;
t198 = t201 ^ 2;
t211 = t213 ^ 2;
t87 = t132 * t225 + t293;
t54 = t222 * t87 + t226 * t86;
t314 = pkin(3) * t54;
t294 = t127 * t225;
t105 = -t151 * t221 - t294;
t68 = t104 * t226 + t105 * t222;
t313 = pkin(3) * t68;
t20 = t221 * t47 - t225 * t46;
t312 = pkin(4) * t20;
t61 = -t221 * t96 - t225 * t99;
t63 = t221 * t99 - t225 * t96;
t34 = t222 * t63 + t226 * t61;
t311 = pkin(8) * t34;
t310 = pkin(8) * t54;
t309 = pkin(8) * t68;
t308 = pkin(9) * t61;
t307 = pkin(9) * t86;
t306 = pkin(9) * t104;
t55 = -t222 * t86 + t226 * t87;
t95 = (qJD(5) + t210) * t178 + t256;
t305 = -pkin(3) * t95 + pkin(8) * t55;
t69 = -t104 * t222 + t105 * t226;
t304 = -pkin(3) * t318 + pkin(8) * t69;
t303 = t20 * t222;
t302 = t20 * t226;
t301 = t221 * t234;
t300 = t225 * t234;
t122 = -t174 - t175;
t35 = -t222 * t61 + t226 * t63;
t299 = -pkin(3) * t122 + pkin(8) * t35;
t298 = qJDD(1) * pkin(1);
t179 = t223 * g(3) + t227 * t190;
t152 = qJDD(3) * pkin(3) + t229 * pkin(8) - t227 * t239 + t179;
t106 = t170 * pkin(4) + t197 * pkin(9) - t201 * t249 + t152;
t297 = t106 * t221;
t296 = t106 * t225;
t291 = t152 * t222;
t290 = t152 * t226;
t162 = -t183 + t199;
t289 = t162 * t222;
t288 = t162 * t226;
t287 = t210 * t221;
t286 = t210 * t225;
t285 = t213 * t222;
t284 = t213 * t226;
t219 = t223 ^ 2;
t283 = t219 * t230;
t220 = t227 ^ 2;
t282 = t220 * t230;
t265 = t223 * t230 * t227;
t281 = t223 * (qJDD(3) + t265);
t280 = t227 * (qJDD(3) - t265);
t278 = t219 + t220;
t274 = -qJD(4) - t213;
t268 = t103 - t47;
t267 = t223 * t137;
t266 = t223 * t183;
t59 = pkin(4) * t61;
t264 = -pkin(3) * t34 - t59;
t11 = t221 * t30 - t300;
t28 = pkin(5) * t234;
t263 = pkin(4) * t11 - t28;
t262 = -pkin(4) * t95 + pkin(9) * t87;
t261 = -pkin(4) * t122 + pkin(9) * t63;
t258 = -pkin(4) * t318 + pkin(9) * t105;
t21 = t221 * t46 + t225 * t47;
t74 = t113 * t222 + t226 * t114;
t26 = -t122 * t227 + t223 * t35;
t254 = qJ(2) * t34 - t315 * t26;
t38 = t223 * t55 - t227 * t95;
t253 = qJ(2) * t54 - t315 * t38;
t41 = t223 * t69 - t227 * t318;
t252 = qJ(2) * t68 - t315 * t41;
t251 = t46 - t85;
t243 = t113 * t226 - t114 * t222;
t140 = t227 * t179 - t223 * t180;
t241 = qJ(2) + t250;
t232 = -t178 * t154 - qJDD(6) + t106;
t231 = -t119 * pkin(5) - t232;
t206 = t278 * qJDD(1);
t205 = t215 - 0.2e1 * t260;
t202 = 0.2e1 * t259 + t269;
t192 = -t236 + t298;
t189 = -t198 + t211;
t188 = t197 - t211;
t187 = t240 - 0.2e1 * t270 + t319;
t185 = -t281 + t227 * (-t229 - t282);
t184 = t223 * (-t229 - t283) + t280;
t182 = t198 - t197;
t181 = -t198 - t211;
t172 = -t211 - t197;
t160 = t197 + t198;
t156 = -t175 + t209;
t155 = t174 - t209;
t147 = t274 * t200 + t242;
t145 = t171 + t191;
t143 = t274 * t201 - t255;
t135 = t175 - t174;
t134 = -t181 * t222 - t288;
t133 = t181 * t226 - t289;
t130 = t172 * t226 - t322;
t129 = t172 * t222 + t321;
t124 = (-t176 * t225 + t178 * t221) * t210;
t123 = (-t176 * t221 - t178 * t225) * t210;
t116 = -t142 * t226 + t146 * t222;
t111 = t155 * t225 - t295;
t110 = -t156 * t221 - t292;
t109 = t155 * t221 + t294;
t108 = t156 * t225 - t293;
t107 = t134 * t223 + t147 * t227;
t101 = t130 * t223 + t143 * t227;
t92 = pkin(5) * t99;
t91 = t120 * t225 - t178 * t287;
t90 = t120 * t221 + t178 * t286;
t89 = -t119 * t221 + t176 * t286;
t88 = t119 * t225 + t176 * t287;
t83 = t116 * t223 + t160 * t227;
t82 = t123 * t226 + t124 * t222;
t79 = t227 * (-t123 * t222 + t124 * t226) + t223 * t196;
t75 = -pkin(5) * t318 - qJ(6) * t127;
t72 = t109 * t226 + t111 * t222;
t71 = t108 * t226 + t110 * t222;
t70 = -t296 - t306;
t65 = t152 * t227 + t223 * t74;
t64 = -t297 - t307;
t62 = -t221 * t318 - t225 * t95;
t60 = -t221 * t95 + t225 * t318;
t57 = t222 * t91 + t226 * t90;
t56 = t222 * t89 + t226 * t88;
t51 = t174 * qJ(6) - t231;
t50 = t227 * (-t222 * t90 + t226 * t91) + t267;
t49 = t227 * (-t222 * t88 + t226 * t89) - t267;
t48 = (-t151 - t174) * qJ(6) + t231;
t45 = t227 * (-t109 * t222 + t111 * t226) - t223 * t96;
t44 = t227 * (-t108 * t222 + t110 * t226) + t223 * t99;
t42 = t258 - t297;
t39 = t262 + t296;
t36 = (t132 + t174) * qJ(6) + (t119 - t95) * pkin(5) + t232;
t33 = t222 * t62 + t226 * t60;
t27 = t227 * (-t222 * t60 + t226 * t62) + t223 * t135;
t24 = -t221 * t75 + t225 * t48 - t306;
t23 = qJ(6) * t292 - t221 * t36 - t307;
t22 = t234 + t323;
t19 = -pkin(5) * t122 - qJ(6) * t96 + t30;
t18 = t221 * t48 + t225 * t75 + t258;
t17 = pkin(4) * t106 + pkin(9) * t21;
t16 = qJ(6) * t293 + t225 * t36 + t262;
t15 = pkin(5) * t51 + qJ(6) * t30;
t14 = -t20 - t308;
t13 = t21 + t261;
t12 = t225 * t30 + t301;
t10 = t21 * t226 - t303;
t9 = t21 * t222 + t302;
t8 = t10 * t223 + t106 * t227;
t7 = -t19 * t221 + t22 * t225 - t308;
t6 = t19 * t225 + t22 * t221 + t261;
t5 = -t11 * t222 + t12 * t226;
t4 = t11 * t226 + t12 * t222;
t3 = -pkin(9) * t11 + qJ(6) * t300 - t15 * t221;
t2 = t223 * t5 + t227 * t51;
t1 = pkin(4) * t51 + pkin(9) * t12 + qJ(6) * t301 + t15 * t225;
t25 = [0, 0, 0, 0, 0, qJDD(1), t257, t248, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t247 - 0.2e1 * t298, t216 + 0.2e1 * t218 - t248, pkin(1) * t192 + qJ(2) * (-t230 * pkin(1) + t238), -t244 * t227, -t202 * t227 - t205 * t223, t280 - t223 * (t229 - t282), t245 * t223, t227 * (-t229 + t283) - t281, 0, qJ(2) * t202 - t315 * t184 - t223 * t187, qJ(2) * t205 - t315 * t185 - t227 * t187, t315 * t206 - t278 * t279 - t140, -qJ(2) * t187 - t315 * t140, t227 * (t171 * t226 - t201 * t285) - t266, t227 * (t143 * t226 - t145 * t222) + t223 * t182, t227 * (-t189 * t222 + t321) + t223 * t146, t227 * (-t170 * t222 - t200 * t284) + t266, t227 * (t188 * t226 - t289) - t223 * t142, t223 * t199 + t227 * (t200 * t226 + t201 * t222) * t213, t227 * (-pkin(8) * t129 - t291) - t223 * (-pkin(3) * t129 + t113) + qJ(2) * t129 - t315 * t101, t227 * (-pkin(8) * t133 - t290) - t223 * (-pkin(3) * t133 + t114) + qJ(2) * t133 - t315 * t107, t227 * t243 - t315 * t83 + t241 * (-t142 * t222 - t146 * t226), -t241 * t243 - t315 * t65, t50, t27, t44, t49, t45, t79, t227 * (-t222 * t39 + t226 * t64 - t310) - t223 * (t251 - t314) + t253, t227 * (-t222 * t42 + t226 * t70 - t309) - t223 * (-t268 - t313) + t252, t227 * (-t13 * t222 + t14 * t226 - t311) - t223 * t264 + t254, t227 * (-pkin(8) * t9 - pkin(9) * t302 - t17 * t222) - t223 * (-pkin(3) * t9 - t312) + qJ(2) * t9 - t315 * t8, t50, t27, t44, t49, t45, t79, t227 * (-t16 * t222 + t226 * t23 - t310) - t223 * (-t314 - t325) + t253, t227 * (-t18 * t222 + t226 * t24 - t309) - t223 * (-t313 - t324) + t252, t227 * (-t222 * t6 + t226 * t7 - t311) - t223 * (t264 + t92) + t254, t227 * (-pkin(8) * t4 - t1 * t222 + t226 * t3) - t223 * (-pkin(3) * t4 - t263) + qJ(2) * t4 - t315 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t230, -t192, 0, 0, 0, 0, 0, 0, t184, t185, -t206, t140, 0, 0, 0, 0, 0, 0, t101, t107, t83, t65, 0, 0, 0, 0, 0, 0, t38, t41, t26, t8, 0, 0, 0, 0, 0, 0, t38, t41, t26, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, (-t219 + t220) * t230, t215, -t265, -t269, qJDD(3), t179, t180, 0, 0, t171 * t222 + t201 * t284, t143 * t222 + t145 * t226, t189 * t226 + t322, t170 * t226 - t200 * t285, t188 * t222 + t288, (t200 * t222 - t201 * t226) * t213, pkin(3) * t143 + pkin(8) * t130 + t290, pkin(3) * t147 + pkin(8) * t134 - t291, pkin(3) * t160 + pkin(8) * t116 + t74, pkin(3) * t152 + pkin(8) * t74, t57, t33, t71, t56, t72, t82, t222 * t64 + t226 * t39 + t305, t222 * t70 + t226 * t42 + t304, t13 * t226 + t14 * t222 + t299, pkin(3) * t106 + pkin(8) * t10 - pkin(9) * t303 + t17 * t226, t57, t33, t71, t56, t72, t82, t16 * t226 + t222 * t23 + t305, t18 * t226 + t222 * t24 + t304, t222 * t7 + t226 * t6 + t299, pkin(3) * t51 + pkin(8) * t5 + t1 * t226 + t222 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, t182, t146, t183, -t142, t199, -t113, -t114, 0, 0, t137, t135, t99, -t137, -t96, t196, -t251, t268, t59, t312, t137, t135, t99, -t137, -t96, t196, t325, t324, -t92 + t59, t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t135, t99, -t137, -t96, t196, -t46, -t47, 0, 0, t137, t135, t99, -t137, -t96, t196, t233, t237, -t92, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t318, t122, -t51;];
tauJ_reg  = t25;
