% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:30:16
% EndTime: 2019-05-05 17:30:26
% DurationCPUTime: 3.70s
% Computational Cost: add. (15577->395), mult. (34134->538), div. (0->0), fcn. (23234->10), ass. (0->247)
t217 = sin(pkin(10));
t219 = cos(pkin(10));
t227 = cos(qJ(3));
t275 = qJD(1) * t227;
t224 = sin(qJ(3));
t276 = qJD(1) * t224;
t187 = t217 * t275 + t219 * t276;
t223 = sin(qJ(5));
t226 = cos(qJ(5));
t171 = -t226 * qJD(3) + t187 * t223;
t173 = qJD(3) * t223 + t187 * t226;
t141 = t173 * t171;
t265 = qJD(1) * qJD(3);
t251 = t227 * t265;
t264 = t224 * qJDD(1);
t193 = t251 + t264;
t209 = t227 * qJDD(1);
t252 = t224 * t265;
t242 = t209 - t252;
t247 = t193 * t217 - t219 * t242;
t161 = qJDD(5) + t247;
t311 = -t141 + t161;
t318 = pkin(5) * t311;
t164 = t219 * t193 + t217 * t242;
t133 = -qJD(5) * t171 + qJDD(3) * t223 + t164 * t226;
t185 = t217 * t276 - t219 * t275;
t182 = qJD(5) + t185;
t152 = t182 * t171;
t109 = t133 + t152;
t317 = qJ(6) * t109;
t162 = t187 * t185;
t310 = qJDD(3) - t162;
t316 = t217 * t310;
t315 = t219 * t310;
t291 = t311 * t223;
t290 = t311 * t226;
t170 = t173 ^ 2;
t181 = t182 ^ 2;
t136 = -t170 - t181;
t169 = t171 ^ 2;
t249 = -t226 * qJDD(3) + t164 * t223;
t132 = -qJD(5) * t173 - t249;
t148 = pkin(5) * t182 - qJ(6) * t173;
t155 = pkin(4) * t185 - pkin(8) * t187;
t229 = qJD(3) ^ 2;
t218 = sin(pkin(9));
t230 = qJD(1) ^ 2;
t225 = sin(qJ(1));
t228 = cos(qJ(1));
t243 = g(1) * t225 - g(2) * t228;
t237 = qJDD(1) * pkin(1) + t243;
t244 = g(1) * t228 + g(2) * t225;
t191 = -pkin(1) * t230 - t244;
t220 = cos(pkin(9));
t283 = t220 * t191;
t233 = -t230 * pkin(2) + qJDD(1) * pkin(7) + t218 * t237 + t283;
t279 = -g(3) + qJDD(2);
t147 = t224 * t279 + t227 * t233;
t197 = qJD(3) * pkin(3) - qJ(4) * t276;
t214 = t227 ^ 2;
t212 = t214 * t230;
t122 = -pkin(3) * t212 + t242 * qJ(4) - qJD(3) * t197 + t147;
t232 = t224 * t233;
t281 = t224 * t230;
t231 = -t232 - t193 * qJ(4) + qJDD(3) * pkin(3) + (pkin(3) * t281 + qJ(4) * t265 + t279) * t227;
t81 = -0.2e1 * qJD(4) * t185 + t219 * t122 + t217 * t231;
t70 = -pkin(4) * t229 + qJDD(3) * pkin(8) - t155 * t185 + t81;
t188 = t220 * t237;
t248 = -t218 * t191 + t188;
t154 = -qJDD(1) * pkin(2) - t230 * pkin(7) - t248;
t127 = -t242 * pkin(3) - qJ(4) * t212 + t197 * t276 + qJDD(4) + t154;
t273 = qJD(3) * t187;
t142 = t247 + t273;
t274 = qJD(3) * t185;
t245 = -t164 + t274;
t84 = pkin(4) * t142 + t245 * pkin(8) + t127;
t46 = t223 * t84 + t226 * t70;
t238 = t132 * qJ(6) - 0.2e1 * qJD(6) * t171 - t148 * t182 + t46;
t314 = -t238 + (t136 + t169) * pkin(5);
t312 = t133 - t152;
t106 = (qJD(5) - t182) * t173 + t249;
t183 = t185 ^ 2;
t184 = t187 ^ 2;
t128 = -t181 - t169;
t88 = t128 * t223 + t290;
t309 = pkin(4) * t88;
t120 = t141 + t161;
t293 = t120 * t223;
t92 = t136 * t226 - t293;
t308 = pkin(4) * t92;
t267 = qJD(6) * t173;
t166 = -0.2e1 * t267;
t45 = t223 * t70 - t226 * t84;
t236 = -t317 - t45 + t318;
t26 = t166 + t236;
t307 = pkin(5) * t26;
t74 = -t106 * t223 - t109 * t226;
t306 = pkin(8) * t74;
t305 = pkin(8) * t88;
t304 = pkin(8) * t92;
t303 = pkin(4) * t217;
t302 = pkin(5) * t109;
t126 = -t169 - t170;
t76 = -t106 * t226 + t109 * t223;
t57 = -t126 * t219 + t217 * t76;
t301 = qJ(4) * t57;
t105 = (qJD(5) + t182) * t173 + t249;
t89 = t128 * t226 - t291;
t61 = -t105 * t219 + t217 * t89;
t300 = qJ(4) * t61;
t292 = t120 * t226;
t93 = -t136 * t223 - t292;
t66 = t217 * t93 - t219 * t312;
t299 = qJ(4) * t66;
t298 = t223 * t26;
t250 = t122 * t217 - t219 * t231;
t239 = -qJDD(3) * pkin(4) - t229 * pkin(8) + t250;
t246 = (0.2e1 * qJD(4) + t155) * t187;
t69 = t246 + t239;
t297 = t223 * t69;
t269 = qJD(4) * t187;
t80 = t250 + 0.2e1 * t269;
t50 = t217 * t81 - t219 * t80;
t296 = t224 * t50;
t295 = t226 * t26;
t294 = t226 * t69;
t289 = t127 * t217;
t288 = t127 * t219;
t158 = qJDD(3) + t162;
t287 = t158 * t217;
t286 = t158 * t219;
t285 = t182 * t223;
t284 = t182 * t226;
t205 = t227 * t281;
t198 = qJDD(3) + t205;
t282 = t224 * t198;
t199 = qJDD(3) - t205;
t280 = t227 * t199;
t272 = qJD(3) * t217;
t271 = qJD(3) * t219;
t58 = t126 * t217 + t219 * t76;
t29 = -t224 * t57 + t227 * t58;
t263 = pkin(1) * (t218 * t29 - t220 * t74) + pkin(7) * t29 - pkin(2) * t74;
t62 = t105 * t217 + t219 * t89;
t37 = -t224 * t61 + t227 * t62;
t262 = pkin(1) * (t218 * t37 - t220 * t88) + pkin(7) * t37 - pkin(2) * t88;
t67 = t217 * t312 + t219 * t93;
t40 = -t224 * t66 + t227 * t67;
t261 = pkin(1) * (t218 * t40 - t220 * t92) + pkin(7) * t40 - pkin(2) * t92;
t260 = pkin(3) * t57 - pkin(4) * t126 + pkin(8) * t76;
t259 = pkin(3) * t61 - pkin(4) * t105 + pkin(8) * t89;
t258 = pkin(3) * t66 - pkin(4) * t312 + pkin(8) * t93;
t257 = t217 * t141;
t256 = t219 * t141;
t255 = -pkin(4) * t219 - pkin(3);
t254 = -pkin(3) * t88 + qJ(4) * t62;
t253 = -pkin(3) * t92 + qJ(4) * t67;
t51 = t217 * t80 + t219 * t81;
t16 = t223 * t45 + t226 * t46;
t146 = -t227 * t279 + t232;
t113 = t224 * t146 + t227 * t147;
t194 = t209 - 0.2e1 * t252;
t15 = t223 * t46 - t226 * t45;
t143 = -t247 + t273;
t235 = t236 + t318;
t234 = -t132 * pkin(5) - t169 * qJ(6) + t148 * t173 + qJDD(6) + t239;
t47 = t246 + t234;
t213 = t224 ^ 2;
t210 = t213 * t230;
t204 = -t212 - t229;
t203 = -t210 - t229;
t196 = t210 + t212;
t195 = (t213 + t214) * qJDD(1);
t192 = 0.2e1 * t251 + t264;
t180 = -0.2e1 * t269;
t178 = -t184 - t229;
t177 = -t184 + t229;
t176 = t183 - t229;
t175 = -t203 * t224 - t280;
t174 = t204 * t227 - t282;
t167 = 0.2e1 * t267;
t156 = -t229 - t183;
t150 = -t170 + t181;
t149 = t169 - t181;
t145 = t164 + t274;
t140 = -t183 - t184;
t138 = t170 - t169;
t135 = -t178 * t217 - t286;
t134 = t178 * t219 - t287;
t124 = t156 * t219 - t316;
t123 = t156 * t217 + t315;
t115 = (-t171 * t226 + t173 * t223) * t182;
t114 = (-t171 * t223 - t173 * t226) * t182;
t112 = t143 * t219 + t145 * t217;
t111 = t143 * t217 - t145 * t219;
t102 = t133 * t226 - t173 * t285;
t101 = t133 * t223 + t173 * t284;
t100 = -t132 * t223 + t171 * t284;
t99 = t132 * t226 + t171 * t285;
t98 = -t224 * t134 + t227 * t135;
t97 = t149 * t226 - t293;
t96 = -t150 * t223 + t290;
t95 = t149 * t223 + t292;
t94 = t150 * t226 + t291;
t85 = -t224 * t123 + t227 * t124;
t78 = -pkin(5) * t312 - qJ(6) * t120;
t77 = -t224 * t111 + t227 * t112;
t75 = -t105 * t226 - t223 * t312;
t73 = -t105 * t223 + t226 * t312;
t63 = t224 * (t115 * t219 + t161 * t217) + t227 * (t115 * t217 - t161 * t219);
t55 = qJ(4) * t58;
t54 = -pkin(4) * t74 + t302;
t53 = t224 * (t102 * t219 + t257) + t227 * (t102 * t217 - t256);
t52 = t224 * (t100 * t219 - t257) + t227 * (t100 * t217 + t256);
t49 = t294 - t304;
t48 = t297 - t305;
t43 = -qJ(6) * t136 + t47;
t42 = t224 * (-t106 * t217 + t219 * t97) + t227 * (t106 * t219 + t217 * t97);
t41 = t224 * (t109 * t217 + t219 * t96) + t227 * (-t109 * t219 + t217 * t96);
t39 = t224 * t67 + t227 * t66;
t36 = t224 * t62 + t227 * t61;
t34 = t224 * (t138 * t217 + t219 * t75) + t227 * (-t138 * t219 + t217 * t75);
t33 = t46 - t308;
t32 = t45 - t309;
t31 = -pkin(5) * t105 + qJ(6) * t128 - t155 * t187 + t180 - t234;
t30 = -pkin(5) * t169 + t238;
t28 = t224 * t58 + t227 * t57;
t23 = t167 - t236 + t317;
t22 = t227 * t51 - t296;
t21 = -qJ(6) * t106 + (-t126 - t169) * pkin(5) + t238;
t20 = -t308 - t314;
t19 = -t223 * t78 + t226 * t43 - t304;
t18 = -qJ(6) * t290 - t223 * t31 - t305;
t14 = t167 - t235 - t309;
t13 = -pkin(5) * t47 + qJ(6) * t30;
t12 = -t15 - t306;
t11 = t16 * t219 + t217 * t69;
t10 = t16 * t217 - t219 * t69;
t9 = t226 * t30 - t298;
t8 = t223 * t30 + t295;
t7 = t217 * t47 + t219 * t9;
t6 = t217 * t9 - t219 * t47;
t5 = -t21 * t223 + t226 * t23 - t306;
t4 = -pkin(4) * t8 - t307;
t2 = -pkin(8) * t8 - qJ(6) * t295 - t13 * t223;
t1 = -t224 * t6 + t227 * t7;
t3 = [0, 0, 0, 0, 0, qJDD(1), t243, t244, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t220 - t218 * t230) + t248, -t283 - t218 * t243 + (-0.2e1 * qJDD(1) * t218 - t220 * t230) * pkin(1), 0, pkin(1) * (t218 ^ 2 * t237 + t220 * t188), (t193 + t251) * t224, t192 * t227 + t194 * t224, t282 + t227 * (-t210 + t229), t194 * t227, t224 * (t212 - t229) + t280, 0, -t227 * t154 + pkin(2) * t194 + pkin(7) * t174 + pkin(1) * (t174 * t218 + t194 * t220), t224 * t154 - pkin(2) * t192 + pkin(7) * t175 + pkin(1) * (t175 * t218 - t192 * t220), pkin(2) * t196 + pkin(7) * t195 + pkin(1) * (t195 * t218 + t196 * t220) + t113, -pkin(2) * t154 + pkin(7) * t113 + pkin(1) * (t113 * t218 - t154 * t220), t224 * (t164 * t219 - t187 * t272) + t227 * (t164 * t217 + t187 * t271), t224 * (-t142 * t219 + t217 * t245) + t227 * (-t142 * t217 - t219 * t245), t224 * (-t177 * t217 + t315) + t227 * (t177 * t219 + t316), t224 * (t185 * t271 + t217 * t247) + t227 * (t185 * t272 - t219 * t247), t224 * (t176 * t219 - t287) + t227 * (t176 * t217 + t286), (t224 * (-t185 * t219 + t187 * t217) + t227 * (-t185 * t217 - t187 * t219)) * qJD(3), t224 * (-qJ(4) * t123 + t289) + t227 * (-pkin(3) * t142 + qJ(4) * t124 - t288) - pkin(2) * t142 + pkin(7) * t85 + pkin(1) * (-t142 * t220 + t218 * t85), t224 * (-qJ(4) * t134 + t288) + t227 * (pkin(3) * t245 + qJ(4) * t135 + t289) + pkin(2) * t245 + pkin(7) * t98 + pkin(1) * (t218 * t98 + t220 * t245), t224 * (-qJ(4) * t111 - t50) + t227 * (-pkin(3) * t140 + qJ(4) * t112 + t51) - pkin(2) * t140 + pkin(7) * t77 + pkin(1) * (-t140 * t220 + t218 * t77), -qJ(4) * t296 + t227 * (-pkin(3) * t127 + qJ(4) * t51) - pkin(2) * t127 + pkin(7) * t22 + pkin(1) * (-t127 * t220 + t218 * t22), t53, t34, t41, t52, t42, t63, t224 * (-t217 * t32 + t219 * t48 - t300) + t227 * (t217 * t48 + t219 * t32 + t254) + t262, t224 * (-t217 * t33 + t219 * t49 - t299) + t227 * (t217 * t49 + t219 * t33 + t253) + t261, t224 * (t12 * t219 + t74 * t303 - t301) + t227 * (t12 * t217 + t255 * t74 + t55) + t263, (t224 * (-pkin(8) * t219 + t303) + t227 * (-pkin(8) * t217 + t255) - pkin(2) - pkin(1) * t220) * t15 + (t218 * pkin(1) + pkin(7) + qJ(4)) * (-t224 * t10 + t227 * t11), t53, t34, t41, t52, t42, t63, t224 * (-t14 * t217 + t18 * t219 - t300) + t227 * (t14 * t219 + t18 * t217 + t254) + t262, t224 * (t19 * t219 - t20 * t217 - t299) + t227 * (t19 * t217 + t20 * t219 + t253) + t261, t224 * (-t217 * t54 + t219 * t5 - t301) + t227 * (-pkin(3) * t74 + t217 * t5 + t219 * t54 + t55) + t263, t224 * (-qJ(4) * t6 + t2 * t219 - t217 * t4) + t227 * (-pkin(3) * t8 + qJ(4) * t7 + t2 * t217 + t219 * t4) - pkin(2) * t8 + pkin(7) * t1 + pkin(1) * (t1 * t218 - t220 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t279, 0, 0, 0, 0, 0, 0, t198 * t227 + t204 * t224, -t199 * t224 + t203 * t227, 0, -t146 * t227 + t147 * t224, 0, 0, 0, 0, 0, 0, t123 * t227 + t124 * t224, t134 * t227 + t135 * t224, t111 * t227 + t112 * t224, t224 * t51 + t227 * t50, 0, 0, 0, 0, 0, 0, t36, t39, t28, t10 * t227 + t11 * t224, 0, 0, 0, 0, 0, 0, t36, t39, t28, t224 * t7 + t227 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t210 - t212, t264, t205, t209, qJDD(3), -t146, -t147, 0, 0, t162, t184 - t183, t145, -t162, t143, qJDD(3), pkin(3) * t123 + t180 - t250, pkin(3) * t134 - t81, pkin(3) * t111, pkin(3) * t50, t101, t73, t94, t99, t95, t114, t259 - t294, t258 + t297, t16 + t260, pkin(3) * t10 - pkin(4) * t69 + pkin(8) * t16, t101, t73, t94, t99, t95, t114, -qJ(6) * t291 + t226 * t31 + t259, t223 * t43 + t226 * t78 + t258, t21 * t226 + t223 * t23 + t260, pkin(3) * t6 - pkin(4) * t47 + pkin(8) * t9 - qJ(6) * t298 + t13 * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, -t245, t140, t127, 0, 0, 0, 0, 0, 0, t88, t92, t74, t15, 0, 0, 0, 0, 0, 0, t88, t92, t74, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, t138, t109, -t141, -t106, t161, -t45, -t46, 0, 0, t141, t138, t109, -t141, -t106, t161, t166 + t235, t314, -t302, t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t312, t126, t47;];
tauJ_reg  = t3;
