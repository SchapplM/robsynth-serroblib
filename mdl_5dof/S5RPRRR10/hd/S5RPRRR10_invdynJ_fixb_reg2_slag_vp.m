% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:50
% EndTime: 2019-12-31 19:11:00
% DurationCPUTime: 5.22s
% Computational Cost: add. (8340->473), mult. (19918->612), div. (0->0), fcn. (15116->14), ass. (0->232)
t192 = cos(pkin(9));
t317 = cos(qJ(3));
t253 = t317 * t192;
t169 = qJD(1) * t253;
t191 = sin(pkin(9));
t196 = sin(qJ(3));
t271 = t196 * t191;
t250 = qJD(1) * t271;
t140 = -t169 + t250;
t133 = qJD(4) + t140;
t189 = pkin(9) + qJ(3);
t180 = sin(t189);
t181 = cos(t189);
t197 = sin(qJ(1));
t199 = cos(qJ(1));
t229 = g(1) * t199 + g(2) * t197;
t210 = -g(3) * t181 + t180 * t229;
t263 = qJD(1) * qJD(2);
t306 = pkin(6) + qJ(2);
t319 = qJDD(1) * t306 + t263;
t127 = t319 * t191;
t128 = t319 * t192;
t159 = t306 * t191;
t149 = qJD(1) * t159;
t160 = t306 * t192;
t150 = qJD(1) * t160;
t247 = qJD(3) * t317;
t267 = qJD(3) * t196;
t236 = t127 * t317 + t128 * t196 - t149 * t267 + t150 * t247;
t290 = qJDD(3) * pkin(3);
t43 = t236 - t290;
t207 = -t43 + t210;
t341 = qJD(4) * pkin(7) * t133 - t207;
t148 = t191 * t317 + t192 * t196;
t142 = t148 * qJD(1);
t195 = sin(qJ(4));
t198 = cos(qJ(4));
t112 = -qJD(3) * t198 + t142 * t195;
t114 = qJD(3) * t195 + t142 * t198;
t194 = sin(qJ(5));
t316 = cos(qJ(5));
t217 = -t112 * t194 + t114 * t316;
t61 = t112 * t316 + t114 * t194;
t307 = t61 * t217;
t276 = t194 * t195;
t216 = t198 * t316 - t276;
t324 = qJD(4) + qJD(5);
t245 = t316 * qJD(5);
t327 = qJD(4) * t316 + t245;
t294 = -t140 * t216 - t198 * t327 + t276 * t324;
t252 = t316 * t195;
t152 = t194 * t198 + t252;
t103 = t324 * t152;
t293 = t140 * t152 + t103;
t200 = -pkin(8) - pkin(7);
t254 = qJD(4) * t200;
t275 = t195 * t140;
t95 = pkin(3) * t142 + pkin(7) * t140;
t139 = t196 * t150;
t99 = -t149 * t317 - t139;
t49 = t195 * t95 + t198 * t99;
t340 = -pkin(8) * t275 + t195 * t254 - t49;
t48 = -t195 * t99 + t198 * t95;
t339 = -pkin(4) * t142 - t48 + (-pkin(8) * t140 + t254) * t198;
t291 = qJDD(1) * pkin(1);
t338 = -g(1) * t197 + g(2) * t199;
t221 = -qJDD(2) + t291 - t338;
t249 = t191 * t267;
t240 = qJDD(1) * t317;
t260 = t192 * qJDD(1);
t255 = qJD(3) * t169 + t191 * t240 + t196 * t260;
t212 = qJD(1) * t249 - t255;
t337 = -qJD(3) * qJD(4) + t212;
t266 = qJD(4) * t195;
t336 = t266 + t275;
t335 = t217 ^ 2 - t61 ^ 2;
t129 = qJD(5) + t133;
t265 = qJD(4) * t198;
t258 = t142 * t265 - t195 * t337;
t220 = qJDD(3) * t198 - t258;
t264 = qJD(5) * t194;
t56 = -t195 * qJDD(3) + t142 * t266 + t198 * t337;
t18 = t112 * t245 + t114 * t264 - t194 * t220 + t316 * t56;
t334 = t129 * t61 - t18;
t174 = pkin(2) * t192 + pkin(1);
t156 = -qJD(1) * t174 + qJD(2);
t71 = pkin(3) * t140 - pkin(7) * t142 + t156;
t100 = -t149 * t196 + t150 * t317;
t94 = qJD(3) * pkin(7) + t100;
t35 = -t195 * t94 + t198 * t71;
t30 = -pkin(8) * t114 + t35;
t26 = pkin(4) * t133 + t30;
t36 = t195 * t71 + t198 * t94;
t31 = -pkin(8) * t112 + t36;
t155 = -qJDD(1) * t174 + qJDD(2);
t145 = t148 * qJD(3);
t261 = t191 * qJDD(1);
t226 = -t192 * t240 + t196 * t261;
t98 = qJD(1) * t145 + t226;
t45 = t98 * pkin(3) + pkin(7) * t212 + t155;
t41 = t198 * t45;
t256 = t127 * t196 - t128 * t317 + t149 * t247;
t46 = -t150 * t267 - t256;
t42 = qJDD(3) * pkin(7) + t46;
t9 = -t36 * qJD(4) - t195 * t42 + t41;
t92 = qJDD(4) + t98;
t6 = pkin(4) * t92 + pkin(8) * t56 + t9;
t8 = t195 * t45 + t198 * t42 + t265 * t71 - t266 * t94;
t7 = pkin(8) * t220 + t8;
t1 = t194 * t6 + t245 * t26 - t264 * t31 + t316 * t7;
t190 = qJ(4) + qJ(5);
t185 = cos(t190);
t278 = t185 * t197;
t184 = sin(t190);
t279 = t184 * t199;
t121 = -t181 * t278 + t279;
t277 = t185 * t199;
t280 = t184 * t197;
t123 = t181 * t277 + t280;
t309 = g(3) * t180;
t93 = -qJD(3) * pkin(3) - t99;
t57 = pkin(4) * t112 + t93;
t333 = g(1) * t123 - g(2) * t121 + t185 * t309 + t57 * t61 - t1;
t227 = -t133 * t35 + t8;
t302 = t133 * t36;
t331 = t9 + t302;
t109 = -t159 * t196 + t160 * t317;
t101 = t198 * t109;
t218 = t253 - t271;
t97 = -pkin(3) * t218 - pkin(7) * t148 - t174;
t55 = t195 * t97 + t101;
t238 = t133 * t195;
t330 = t114 * t238;
t329 = -t159 * t317 - t160 * t196;
t269 = t198 * t199;
t274 = t195 * t197;
t134 = t181 * t274 + t269;
t270 = t197 * t198;
t273 = t195 * t199;
t136 = -t181 * t273 + t270;
t326 = -g(1) * t136 + g(2) * t134;
t325 = qJ(2) * qJDD(1);
t120 = t181 * t280 + t277;
t122 = -t181 * t279 + t278;
t257 = t316 * t31;
t11 = t194 * t26 + t257;
t2 = -qJD(5) * t11 - t194 * t7 + t316 * t6;
t323 = -g(1) * t122 + g(2) * t120 + t184 * t309 - t217 * t57 + t2;
t19 = qJD(5) * t217 - t194 * t56 - t220 * t316;
t322 = t129 * t217 - t19;
t89 = qJDD(5) + t92;
t321 = -t129 * t294 + t152 * t89;
t208 = -t181 * t229 - t309;
t320 = -t18 * t216 - t217 * t293;
t318 = t142 ^ 2;
t161 = t199 * t174;
t311 = g(2) * t161;
t162 = t200 * t195;
t163 = t200 * t198;
t110 = t162 * t316 + t163 * t194;
t305 = t110 * qJD(5) + t194 * t339 + t316 * t340;
t111 = t162 * t194 - t163 * t316;
t304 = -t111 * qJD(5) - t194 * t340 + t316 * t339;
t301 = t142 * t61;
t300 = t142 * t217;
t297 = t194 * t31;
t296 = t195 * t92;
t295 = t56 * t195;
t51 = t195 * t220;
t292 = -t112 * t265 + t51;
t289 = t112 * t140;
t288 = t112 * t142;
t287 = t114 * t112;
t286 = t114 * t142;
t285 = t142 * t140;
t144 = -t192 * t247 + t249;
t284 = t144 * t195;
t283 = t144 * t198;
t282 = t148 * t195;
t281 = t148 * t198;
t187 = t191 ^ 2;
t188 = t192 ^ 2;
t268 = t187 + t188;
t248 = t148 * t266;
t244 = pkin(4) * t195 + t306;
t72 = qJD(2) * t218 + qJD(3) * t329;
t96 = pkin(3) * t145 + pkin(7) * t144;
t241 = -t195 * t72 + t198 * t96;
t54 = -t109 * t195 + t198 * t97;
t239 = t268 * qJD(1) ^ 2;
t237 = t133 * t198;
t235 = -t152 * t19 + t294 * t61;
t234 = 0.2e1 * t268;
t233 = -t129 * t293 + t216 * t89;
t232 = pkin(4) * t336 - t100;
t231 = t338 * t180;
t230 = pkin(3) * t181 + pkin(7) * t180;
t225 = t195 * t36 + t198 * t35;
t178 = pkin(4) * t198 + pkin(3);
t224 = t178 * t181 - t180 * t200;
t222 = -t133 * t336 + t198 * t92;
t33 = -pkin(4) * t218 - pkin(8) * t281 + t54;
t38 = -pkin(8) * t282 + t55;
t21 = -t194 * t38 + t316 * t33;
t22 = t194 * t33 + t316 * t38;
t219 = -pkin(7) * t92 + t133 * t93;
t215 = t148 * t265 - t284;
t214 = -t248 - t283;
t23 = -t109 * t266 + t195 * t96 + t198 * t72 + t265 * t97;
t213 = t220 * t198;
t211 = t221 + t291;
t206 = t234 * t263 - t229;
t73 = qJD(2) * t148 + qJD(3) * t109;
t138 = t140 ^ 2;
t137 = t181 * t269 + t274;
t135 = -t181 * t270 + t273;
t84 = t216 * t148;
t83 = t152 * t148;
t74 = pkin(4) * t282 - t329;
t44 = pkin(4) * t215 + t73;
t28 = -t144 * t252 - t194 * t248 - t264 * t282 + (-t144 * t194 + t148 * t327) * t198;
t27 = t103 * t148 + t144 * t216;
t25 = -pkin(4) * t220 + t43;
t24 = -qJD(4) * t55 + t241;
t20 = -pkin(8) * t215 + t23;
t15 = pkin(8) * t283 + pkin(4) * t145 + (-t101 + (pkin(8) * t148 - t97) * t195) * qJD(4) + t241;
t13 = t30 * t316 - t297;
t12 = -t194 * t30 - t257;
t10 = t26 * t316 - t297;
t4 = -qJD(5) * t22 + t15 * t316 - t194 * t20;
t3 = qJD(5) * t21 + t194 * t15 + t20 * t316;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t338, t229, 0, 0, t187 * qJDD(1), 0.2e1 * t191 * t260, 0, t188 * qJDD(1), 0, 0, t211 * t192, -t211 * t191, t234 * t325 + t206, t221 * pkin(1) + (t268 * t325 + t206) * qJ(2), -t142 * t144 - t148 * t212, t140 * t144 - t142 * t145 - t148 * t98 - t212 * t218, -qJD(3) * t144 + qJDD(3) * t148, t140 * t145 - t218 * t98, -qJD(3) * t145 + qJDD(3) * t218, 0, -qJD(3) * t73 + qJDD(3) * t329 + t145 * t156 - t155 * t218 - t174 * t98 - t181 * t338, -t72 * qJD(3) - t109 * qJDD(3) - t156 * t144 + t155 * t148 + t174 * t212 + t231, -t100 * t145 - t109 * t98 - t140 * t72 + t142 * t73 + t144 * t99 + t148 * t236 + t212 * t329 + t218 * t46 - t229, t46 * t109 + t100 * t72 - t236 * t329 - t99 * t73 - t155 * t174 - g(1) * (-t174 * t197 + t199 * t306) - g(2) * (t197 * t306 + t161), t114 * t214 - t281 * t56, (t112 * t198 + t114 * t195) * t144 + (t213 + t295 + (t112 * t195 - t114 * t198) * qJD(4)) * t148, t114 * t145 + t133 * t214 + t218 * t56 + t281 * t92, t112 * t215 - t148 * t51, -t112 * t145 - t133 * t215 - t218 * t220 - t282 * t92, t133 * t145 - t218 * t92, t24 * t133 + t54 * t92 - t9 * t218 + t35 * t145 + t73 * t112 + t329 * t220 - t93 * t284 - g(1) * t135 - g(2) * t137 + (t195 * t43 + t265 * t93) * t148, -t93 * t283 - g(1) * t134 - g(2) * t136 + t329 * t56 + t114 * t73 - t133 * t23 - t145 * t36 + t218 * t8 - t55 * t92 + (t43 * t198 - t266 * t93) * t148, -t23 * t112 + t55 * t220 - t24 * t114 + t54 * t56 + t225 * t144 + (-t8 * t195 - t9 * t198 + (t195 * t35 - t198 * t36) * qJD(4)) * t148 - t231, -t311 - t43 * t329 + t36 * t23 + t35 * t24 + t9 * t54 + t8 * t55 + t93 * t73 + (-g(1) * t306 - g(2) * t230) * t199 + (-g(1) * (-t174 - t230) - g(2) * t306) * t197, -t18 * t84 - t217 * t27, t18 * t83 - t19 * t84 - t217 * t28 + t27 * t61, -t129 * t27 + t145 * t217 + t18 * t218 + t84 * t89, t19 * t83 + t28 * t61, -t129 * t28 - t145 * t61 + t19 * t218 - t83 * t89, t129 * t145 - t218 * t89, -g(1) * t121 - g(2) * t123 + t10 * t145 + t129 * t4 + t19 * t74 - t2 * t218 + t21 * t89 + t25 * t83 + t28 * t57 + t44 * t61, -g(1) * t120 - g(2) * t122 + t1 * t218 - t11 * t145 - t129 * t3 - t18 * t74 + t217 * t44 - t22 * t89 + t25 * t84 - t27 * t57, -t1 * t83 + t10 * t27 - t11 * t28 + t18 * t21 - t19 * t22 - t2 * t84 - t217 * t4 - t3 * t61 - t231, -t311 + t1 * t22 + t10 * t4 + t11 * t3 + t2 * t21 + t25 * t74 + t57 * t44 + (-g(1) * t244 - g(2) * t224) * t199 + (-g(1) * (-t174 - t224) - g(2) * t244) * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t260, t261, -t239, -qJ(2) * t239 - t221, 0, 0, 0, 0, 0, 0, 0.2e1 * t142 * qJD(3) + t226, (-t140 - t250) * qJD(3) + t255, -t138 - t318, t100 * t140 + t142 * t99 + t155 + t338, 0, 0, 0, 0, 0, 0, t222 - t288, -t133 ^ 2 * t198 - t286 - t296, (t56 - t289) * t198 + t330 + t292, -t142 * t93 + t195 * t227 + t198 * t331 + t338, 0, 0, 0, 0, 0, 0, t233 - t301, -t300 - t321, t235 - t320, t1 * t152 - t10 * t293 - t11 * t294 - t142 * t57 + t2 * t216 + t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t138 + t318, (t140 - t250) * qJD(3) + t255, -t285, -t226, qJDD(3), qJD(3) * t100 - t142 * t156 + t210 - t236, t140 * t156 + (t99 + t139) * qJD(3) + t256 - t208, 0, 0, t114 * t237 - t295, (-t56 - t289) * t198 - t330 + t292, t133 * t237 - t286 + t296, t112 * t238 + t213, t222 + t288, -t133 * t142, -pkin(3) * t258 - t48 * t133 - t35 * t142 - t100 * t112 + t219 * t195 + (t290 - t341) * t198, pkin(3) * t56 - t100 * t114 + t133 * t49 + t142 * t36 + t195 * t341 + t219 * t198, t49 * t112 + t48 * t114 + ((qJD(4) * t114 + t220) * pkin(7) + t227) * t198 + ((qJD(4) * t112 - t56) * pkin(7) - t331) * t195 + t208, -t93 * t100 - t35 * t48 - t36 * t49 + t207 * pkin(3) + (-qJD(4) * t225 - t9 * t195 + t8 * t198 + t208) * pkin(7), -t152 * t18 - t217 * t294, t235 + t320, -t300 + t321, -t19 * t216 + t293 * t61, t233 + t301, -t129 * t142, -t10 * t142 + t110 * t89 + t129 * t304 - t178 * t19 + t185 * t210 - t216 * t25 + t232 * t61 + t293 * t57, t11 * t142 - t111 * t89 - t129 * t305 + t152 * t25 + t178 * t18 - t184 * t210 + t217 * t232 - t294 * t57, t1 * t216 + t10 * t294 - t11 * t293 + t110 * t18 - t111 * t19 - t152 * t2 - t217 * t304 - t305 * t61 + t208, -g(3) * t224 + t1 * t111 + t10 * t304 + t11 * t305 + t2 * t110 - t25 * t178 + t232 * t57 + t229 * (t178 * t180 + t181 * t200); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, -t112 ^ 2 + t114 ^ 2, t112 * t133 - t56, -t287, t114 * t133 + t220, t92, -t94 * t265 - t114 * t93 + t302 + t41 + (-qJD(4) * t71 + t309 - t42) * t195 + t326, g(1) * t137 - g(2) * t135 + t112 * t93 + t198 * t309 - t227, 0, 0, t307, t335, t334, -t307, t322, t89, -t12 * t129 + (-t114 * t61 - t129 * t264 + t316 * t89) * pkin(4) + t323, t13 * t129 + (-t114 * t217 - t129 * t245 - t194 * t89) * pkin(4) + t333, -t10 * t61 + t11 * t217 + t12 * t217 + t13 * t61 + (t316 * t18 - t19 * t194 + (t194 * t217 - t316 * t61) * qJD(5)) * pkin(4), -t10 * t12 - t11 * t13 + (t1 * t194 + t2 * t316 - t57 * t114 + t195 * t309 + (-t10 * t194 + t11 * t316) * qJD(5) + t326) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, t335, t334, -t307, t322, t89, t11 * t129 + t323, t10 * t129 + t333, 0, 0;];
tau_reg = t5;
