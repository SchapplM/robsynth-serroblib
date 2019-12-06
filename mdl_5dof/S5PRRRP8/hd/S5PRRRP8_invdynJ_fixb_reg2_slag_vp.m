% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:42
% EndTime: 2019-12-05 17:00:53
% DurationCPUTime: 4.39s
% Computational Cost: add. (3235->449), mult. (7749->596), div. (0->0), fcn. (5932->10), ass. (0->219)
t155 = sin(qJ(4));
t310 = pkin(7) * t155;
t158 = cos(qJ(4));
t250 = qJD(3) * t155;
t156 = sin(qJ(3));
t253 = qJD(2) * t156;
t113 = t158 * t253 + t250;
t159 = cos(qJ(3));
t251 = qJD(2) * t159;
t140 = -qJD(4) + t251;
t270 = t113 * t140;
t243 = t158 * qJD(3);
t111 = t155 * t253 - t243;
t272 = t111 * t140;
t241 = qJD(2) * qJD(3);
t223 = t159 * t241;
t239 = t156 * qJDD(2);
t246 = qJD(4) * t156;
t303 = qJD(2) * t246 - qJDD(3);
t50 = -qJD(4) * t243 + (-t223 - t239) * t158 + t303 * t155;
t51 = ((qJD(4) + t251) * qJD(3) + t239) * t155 + t303 * t158;
t309 = (t50 - t272) * t158 + (t51 - t270) * t155;
t153 = sin(pkin(5));
t157 = sin(qJ(2));
t160 = cos(qJ(2));
t242 = qJD(1) * qJD(2);
t225 = t160 * t242;
t154 = cos(pkin(5));
t254 = qJD(1) * t154;
t274 = qJDD(2) * pkin(7);
t308 = qJD(3) * t254 + t274 + (qJDD(1) * t157 + t225) * t153;
t305 = t51 + t270;
t288 = pkin(8) * t156;
t304 = -pkin(3) * t159 - t288;
t248 = qJD(3) * t159;
t278 = t51 * t158;
t279 = t50 * t155;
t302 = ((t111 * t155 - t113 * t158) * qJD(4) - t278 + t279) * t156 - (t111 * t158 + t113 * t155) * t248;
t206 = pkin(3) * t156 - pkin(8) * t159;
t116 = t206 * qJD(3);
t123 = -pkin(2) + t304;
t244 = qJD(4) * t159;
t245 = qJD(4) * t158;
t261 = t159 * t160;
t80 = (t155 * t157 + t158 * t261) * t153;
t285 = t155 * t116 + t123 * t245 + (-t155 * t244 - t156 * t243) * pkin(7) - qJD(1) * t80;
t203 = pkin(4) * t158 + qJ(5) * t155;
t255 = qJD(1) * t153;
t233 = t157 * t255;
t117 = qJD(2) * pkin(7) + t233;
t108 = t156 * t117;
t76 = t159 * t254 - t108;
t68 = -qJD(3) * pkin(3) - t76;
t26 = pkin(4) * t111 - qJ(5) * t113 + t68;
t147 = t159 * qJDD(2);
t109 = t156 * t241 + qJDD(4) - t147;
t290 = pkin(8) * t109;
t300 = t140 * t26 + t290;
t161 = qJD(3) ^ 2;
t277 = cos(pkin(9));
t219 = t277 * t157;
t152 = sin(pkin(9));
t268 = t152 * t160;
t100 = t154 * t268 + t219;
t218 = t277 * t160;
t269 = t152 * t157;
t98 = -t154 * t218 + t269;
t208 = g(1) * t100 + g(2) * t98;
t226 = t157 * t242;
t275 = qJDD(2) * pkin(2);
t265 = t153 * t160;
t201 = -qJDD(1) * t265 + t153 * t226;
t81 = t201 - t275;
t299 = -pkin(7) * t161 + t153 * (-g(3) * t160 + t226) + t208 + t275 - t81;
t247 = qJD(4) * t155;
t187 = t158 * t109 + t140 * t247;
t264 = t155 * t159;
t298 = qJD(2) * (-t111 * t156 + t140 * t264) - t187;
t235 = t155 * t265;
t214 = t159 * t235;
t249 = qJD(3) * t156;
t284 = t249 * t310 + qJD(1) * t214 - t123 * t247 + (-t244 * pkin(7) + t116 - t233) * t158;
t188 = t155 * t109 - t140 * t245;
t295 = t156 * (-qJD(3) * t111 - t188) + t159 * (t140 * t250 + t51);
t294 = t113 ^ 2;
t292 = pkin(4) * t109;
t289 = pkin(8) * t113;
t287 = qJ(5) * t249 - t159 * qJD(5) + t285;
t286 = -pkin(4) * t249 - t284;
t202 = pkin(4) * t155 - qJ(5) * t158;
t228 = t156 * t254;
t283 = t202 * qJD(4) - t155 * qJD(5) - t228 - (t202 * qJD(2) + t117) * t159;
t282 = qJD(2) * pkin(2);
t77 = t159 * t117 + t228;
t69 = qJD(3) * pkin(8) + t77;
t229 = t160 * t255;
t78 = t123 * qJD(2) - t229;
t25 = t155 * t78 + t158 * t69;
t14 = -qJ(5) * t140 + t25;
t281 = t14 * t140;
t280 = t140 * t25;
t115 = t206 * qJD(2);
t47 = t155 * t115 + t158 * t76;
t273 = t109 * qJ(5);
t271 = t113 * t111;
t267 = t153 * t157;
t266 = t153 * t159;
t263 = t158 * t123;
t262 = t158 * t159;
t24 = -t155 * t69 + t158 * t78;
t260 = qJD(5) - t24;
t259 = qJDD(1) - g(3);
t258 = pkin(2) * t265 + pkin(7) * t267;
t84 = pkin(7) * t262 + t155 * t123;
t150 = t156 ^ 2;
t151 = t159 ^ 2;
t257 = t150 - t151;
t256 = t150 + t151;
t252 = qJD(2) * t157;
t240 = qJDD(1) * t154;
t94 = t98 * pkin(2);
t238 = t304 * t98 - t94;
t95 = t100 * pkin(2);
t237 = t304 * t100 - t95;
t236 = -t156 * t240 - t308 * t159;
t162 = qJD(2) ^ 2;
t234 = t156 * t162 * t159;
t232 = t153 * t252;
t231 = qJD(2) * t265;
t230 = t140 * t253;
t227 = t111 ^ 2 - t294;
t222 = t160 * t241;
t220 = t153 * t277;
t20 = -t117 * t249 - t236;
t17 = qJDD(3) * pkin(8) + t20;
t40 = qJD(2) * t116 + t123 * qJDD(2) + t201;
t217 = t155 * t17 - t158 * t40 + t69 * t245 + t78 * t247;
t216 = t117 * t248 + t308 * t156 - t159 * t240;
t213 = t153 * pkin(3) * t261 + t265 * t288 + t258;
t212 = t111 * t229;
t211 = t113 * t229;
t210 = t156 * t229;
t209 = t156 * t223;
t101 = -t154 * t269 + t218;
t99 = t154 * t219 + t268;
t207 = g(1) * t101 + g(2) * t99;
t204 = (qJD(4) * t111 - t50) * pkin(8);
t13 = pkin(4) * t140 + t260;
t200 = t13 * t158 - t14 * t155;
t199 = -t155 * t25 - t158 * t24;
t46 = t115 * t158 - t155 * t76;
t193 = qJDD(2) * t160 - t157 * t162;
t192 = pkin(7) + t202;
t104 = t154 * t156 + t157 * t266;
t63 = t104 * t155 + t158 * t265;
t103 = -t154 * t159 + t156 * t267;
t3 = t155 * t40 + t158 * t17 + t78 * t245 - t69 * t247;
t60 = -t103 * qJD(3) + t159 * t231;
t10 = -qJD(4) * t235 + t104 * t245 + t155 * t60 - t158 * t232;
t11 = -t63 * qJD(4) + t155 * t232 + t60 * t158;
t64 = t104 * t158 - t235;
t185 = t10 * t113 - t11 * t111 - t50 * t63 - t64 * t51;
t41 = -t99 * t158 - t98 * t264;
t43 = -t100 * t264 - t101 * t158;
t79 = -t158 * t267 + t214;
t184 = g(1) * t43 + g(2) * t41 + g(3) * t79;
t42 = t155 * t99 - t98 * t262;
t44 = -t100 * t262 + t101 * t155;
t183 = -g(1) * t44 - g(2) * t42 - g(3) * t80;
t18 = -qJDD(3) * pkin(3) + t216;
t56 = -t99 * t156 - t159 * t220;
t58 = -t101 * t156 + t152 * t266;
t182 = g(1) * t58 + g(2) * t56 - g(3) * t103;
t57 = -t156 * t220 + t159 * t99;
t59 = t152 * t153 * t156 + t101 * t159;
t181 = -g(1) * t59 - g(2) * t57 - g(3) * t104;
t61 = t104 * qJD(3) + t156 * t231;
t180 = t10 * t140 + t103 * t51 - t109 * t63 + t61 * t111;
t178 = -pkin(8) * t278 + t181;
t177 = -g(3) * t265 + t208;
t176 = -t140 * t68 - t290;
t174 = -t155 * t272 - t278;
t173 = t103 * t50 + t109 * t64 - t11 * t140 - t113 * t61;
t28 = t155 * t57 - t98 * t158;
t30 = -t100 * t158 + t155 * t59;
t172 = g(1) * t30 + g(2) * t28 + g(3) * t63 - t217;
t171 = pkin(8) * qJD(4) * t140 - t182;
t5 = pkin(4) * t51 + qJ(5) * t50 - qJD(5) * t113 + t18;
t169 = t171 - t5;
t168 = t171 - t18;
t118 = -t229 - t282;
t167 = -pkin(7) * qJDD(3) + (t118 + t229 - t282) * qJD(3);
t29 = t155 * t98 + t158 * t57;
t31 = t100 * t155 + t158 * t59;
t166 = -g(1) * t31 - g(2) * t29 - g(3) * t64 + t3;
t165 = t51 * t155 * t156 + (t155 * t248 + t156 * t245) * t111;
t164 = t113 * t26 + qJDD(5) - t172;
t163 = t216 * t156 + t20 * t159 + (-t156 * t77 - t159 * t76) * qJD(3) - t207;
t119 = -pkin(3) - t203;
t96 = t103 * pkin(3);
t93 = t192 * t156;
t83 = -pkin(7) * t264 + t263;
t71 = -t263 + (pkin(4) + t310) * t159;
t70 = -qJ(5) * t159 + t84;
t62 = -t109 * t159 - t140 * t249;
t55 = pkin(4) * t113 + qJ(5) * t111;
t53 = t58 * pkin(3);
t52 = t56 * pkin(3);
t37 = (t203 * qJD(4) - qJD(5) * t158) * t156 + t192 * t248;
t36 = -pkin(4) * t253 - t46;
t35 = qJ(5) * t253 + t47;
t23 = -t50 - t272;
t19 = (-t113 * t156 + t140 * t262) * qJD(2) + t188;
t9 = -t158 * t270 - t279;
t8 = -t50 * t158 * t156 + (-t155 * t246 + t159 * t243) * t113;
t6 = (-t140 * t243 + t50) * t159 + (qJD(3) * t113 + t187) * t156;
t2 = qJDD(5) + t217 - t292;
t1 = -qJD(5) * t140 + t273 + t3;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t259, 0, 0, 0, 0, 0, 0, t193 * t153, (-qJDD(2) * t157 - t160 * t162) * t153, 0, -g(3) + (t154 ^ 2 + (t157 ^ 2 + t160 ^ 2) * t153 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t61 * qJD(3) - t103 * qJDD(3) + (-t156 * t222 + t159 * t193) * t153, -t60 * qJD(3) - t104 * qJDD(3) + (-t156 * t193 - t159 * t222) * t153, (t103 * t156 + t104 * t159) * qJDD(2) + (t156 * t61 + t159 * t60 + (t103 * t159 - t104 * t156) * qJD(3)) * qJD(2), t216 * t103 + t20 * t104 + t77 * t60 - t76 * t61 - g(3) + (t118 * t252 - t160 * t81) * t153, 0, 0, 0, 0, 0, 0, t180, -t173, t185, -t10 * t24 + t103 * t18 + t11 * t25 + t217 * t63 + t3 * t64 + t61 * t68 - g(3), 0, 0, 0, 0, 0, 0, t180, t185, t173, t1 * t64 + t10 * t13 + t103 * t5 + t11 * t14 + t2 * t63 + t26 * t61 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t259 * t265 + t208, -t259 * t267 + t207, 0, 0, qJDD(2) * t150 + 0.2e1 * t209, 0.2e1 * t156 * t147 - 0.2e1 * t257 * t241, qJDD(3) * t156 + t159 * t161, qJDD(2) * t151 - 0.2e1 * t209, qJDD(3) * t159 - t156 * t161, 0, t167 * t156 + t299 * t159, -t299 * t156 + t167 * t159, t256 * t274 + (-g(3) * t157 - t256 * t225) * t153 + t163, -t81 * pkin(2) + g(1) * t95 + g(2) * t94 - g(3) * t258 + (-t118 * t157 + (t156 * t76 - t159 * t77) * t160) * t255 + t163 * pkin(7), t8, t302, t6, t165, t295, t62, t83 * t109 - t284 * t140 + (t217 + (pkin(7) * t111 + t155 * t68) * qJD(3)) * t159 + (pkin(7) * t51 + qJD(3) * t24 + t18 * t155 + t245 * t68 - t212) * t156 + t183, -t84 * t109 + t285 * t140 + (t3 + (pkin(7) * t113 + t158 * t68) * qJD(3)) * t159 + (-pkin(7) * t50 - qJD(3) * t25 + t18 * t158 - t247 * t68 - t211) * t156 + t184, t83 * t50 - t84 * t51 - t284 * t113 - t285 * t111 + t199 * t248 + (-t155 * t3 + t158 * t217 + (t155 * t24 - t158 * t25) * qJD(4) + t177) * t156, t3 * t84 - t217 * t83 - t68 * t210 - g(1) * t237 - g(2) * t238 - g(3) * t213 + t285 * t25 + t284 * t24 + (t156 * t18 + t248 * t68 - t207) * pkin(7), t8, t6, -t302, t62, -t295, t165, -t71 * t109 + t37 * t111 + t93 * t51 + (t250 * t26 + t2) * t159 + t286 * t140 + (-qJD(3) * t13 + t5 * t155 + t245 * t26 - t212) * t156 + t183, -t71 * t50 - t70 * t51 + t286 * t113 - t287 * t111 + t200 * t248 + (-t1 * t155 + t158 * t2 + (-t13 * t155 - t14 * t158) * qJD(4) + t177) * t156, t70 * t109 - t37 * t113 + t93 * t50 + (-t243 * t26 - t1) * t159 - t287 * t140 + (qJD(3) * t14 - t5 * t158 + t247 * t26 + t211) * t156 - t184, t1 * t70 + t5 * t93 + t2 * t71 - g(1) * (pkin(4) * t44 + pkin(7) * t101 + qJ(5) * t43 + t237) - g(2) * (pkin(4) * t42 + pkin(7) * t99 + qJ(5) * t41 + t238) - g(3) * (pkin(4) * t80 + qJ(5) * t79 + t213) + (t37 - t210) * t26 + t287 * t14 + t286 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t257 * t162, t239, t234, t147, qJDD(3), qJD(3) * t77 - t118 * t253 - t182 - t216, -t118 * t251 + (t76 + t108) * qJD(3) - t181 + t236, 0, 0, t9, -t309, t19, t174, -t298, t230, -pkin(3) * t51 - t77 * t111 + t46 * t140 + t155 * t176 + t158 * t168 - t24 * t253, pkin(3) * t50 - t77 * t113 - t47 * t140 - t155 * t168 + t158 * t176 + t25 * t253, t47 * t111 + t46 * t113 + (t24 * t251 + t3 + (-t24 + t289) * qJD(4)) * t158 + (t204 + t217 + t280) * t155 + t178, -t18 * pkin(3) - g(1) * t53 - g(2) * t52 + g(3) * t96 - t24 * t46 - t25 * t47 - t68 * t77 + (qJD(4) * t199 + t155 * t217 + t3 * t158 + t181) * pkin(8), t9, t19, t309, t230, t298, t174, t283 * t111 + t119 * t51 + t13 * t253 - t36 * t140 - t300 * t155 + t169 * t158, t35 * t111 - t36 * t113 + (-t13 * t251 + t1 + (t13 + t289) * qJD(4)) * t158 + (t2 + t204 + t281) * t155 + t178, -t283 * t113 + t119 * t50 - t14 * t253 + t35 * t140 + t169 * t155 + t300 * t158, t5 * t119 - t14 * t35 - t13 * t36 - g(1) * (t203 * t58 + t53) - g(2) * (t203 * t56 + t52) - g(3) * (-t203 * t103 - t96) + t283 * t26 + (qJD(4) * t200 + t1 * t158 + t2 * t155 + t181) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, -t227, t23, -t271, -t305, t109, -t113 * t68 + t172 - t280, t111 * t68 - t140 * t24 - t166, 0, 0, t271, t23, t227, t109, t305, -t271, -t111 * t55 - t164 - t280 + 0.2e1 * t292, pkin(4) * t50 - t51 * qJ(5) + (t14 - t25) * t113 + (t13 - t260) * t111, 0.2e1 * t273 - t26 * t111 + t55 * t113 + (-0.2e1 * qJD(5) + t24) * t140 + t166, t1 * qJ(5) - t2 * pkin(4) - t26 * t55 - t13 * t25 - g(1) * (-pkin(4) * t30 + qJ(5) * t31) - g(2) * (-pkin(4) * t28 + qJ(5) * t29) - g(3) * (-pkin(4) * t63 + qJ(5) * t64) + t260 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 + t271, t23, -t140 ^ 2 - t294, t164 + t281 - t292;];
tau_reg = t4;
