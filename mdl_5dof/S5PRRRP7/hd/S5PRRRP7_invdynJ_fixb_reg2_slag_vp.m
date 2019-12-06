% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRP7
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:33
% EndTime: 2019-12-05 16:56:41
% DurationCPUTime: 3.60s
% Computational Cost: add. (3178->438), mult. (7668->592), div. (0->0), fcn. (5834->10), ass. (0->221)
t151 = sin(pkin(5));
t156 = sin(qJ(2));
t159 = cos(qJ(2));
t222 = qJD(1) * qJD(2);
t207 = t159 * t222;
t152 = cos(pkin(5));
t235 = qJD(1) * t152;
t257 = qJDD(2) * pkin(7);
t302 = qJD(3) * t235 + t257 + (qJDD(1) * t156 + t207) * t151;
t155 = sin(qJ(3));
t226 = qJD(4) * t155;
t301 = qJD(2) * t226 - qJDD(3);
t158 = cos(qJ(3));
t249 = t151 * t156;
t100 = -t152 * t158 + t155 * t249;
t260 = cos(pkin(9));
t201 = t151 * t260;
t259 = sin(pkin(9));
t196 = t259 * t159;
t199 = t260 * t156;
t96 = t152 * t199 + t196;
t61 = t96 * t155 + t158 * t201;
t200 = t151 * t259;
t197 = t259 * t156;
t198 = t260 * t159;
t98 = -t152 * t197 + t198;
t63 = t155 * t98 - t158 * t200;
t173 = g(1) * t63 + g(2) * t61 + g(3) * t100;
t154 = sin(qJ(4));
t157 = cos(qJ(4));
t219 = t155 * qJDD(2);
t231 = qJD(2) * t158;
t55 = ((qJD(4) + t231) * qJD(3) + t219) * t154 + t301 * t157;
t186 = pkin(3) * t155 - pkin(8) * t158;
t117 = t186 * qJD(3);
t229 = qJD(3) * t155;
t282 = pkin(7) * t154;
t242 = t158 * t159;
t84 = (-t154 * t242 + t156 * t157) * t151;
t300 = -qJD(1) * t84 + t157 * t117 + t229 * t282;
t187 = pkin(3) * t158 + pkin(8) * t155;
t121 = -pkin(2) - t187;
t225 = qJD(4) * t157;
t246 = t154 * t156;
t85 = (t157 * t242 + t246) * t151;
t299 = -qJD(1) * t85 + t154 * t117 + t121 * t225;
t224 = t157 * qJD(3);
t233 = qJD(2) * t155;
t112 = t154 * t233 - t224;
t139 = -qJD(4) + t231;
t255 = t112 * t139;
t221 = qJD(2) * qJD(3);
t205 = t158 * t221;
t54 = -qJD(4) * t224 + (-t205 - t219) * t157 + t301 * t154;
t298 = -t54 + t255;
t230 = qJD(3) * t154;
t114 = t157 * t233 + t230;
t252 = t114 * t139;
t297 = t55 - t252;
t243 = t157 * t158;
t140 = pkin(7) * t243;
t89 = t154 * t121 + t140;
t296 = pkin(4) * t55 + qJDD(5);
t62 = -t155 * t201 + t158 * t96;
t64 = t155 * t200 + t98 * t158;
t101 = t152 * t155 + t158 * t249;
t248 = t151 * t159;
t68 = -t101 * t154 - t157 * t248;
t95 = -t152 * t198 + t197;
t97 = t152 * t196 + t199;
t295 = -g(3) * t68 - g(2) * (-t154 * t62 + t157 * t95) - g(1) * (-t154 * t64 + t157 * t97);
t160 = qJD(3) ^ 2;
t189 = g(1) * t97 + g(2) * t95;
t208 = t156 * t222;
t258 = qJDD(2) * pkin(2);
t185 = -qJDD(1) * t248 + t151 * t208;
t86 = t185 - t258;
t294 = -pkin(7) * t160 + t151 * (-g(3) * t159 + t208) + t189 + t258 - t86;
t292 = t114 ^ 2;
t284 = pkin(4) * t112;
t283 = pkin(4) * t154;
t146 = t158 * qJDD(2);
t109 = t155 * t221 + qJDD(4) - t146;
t279 = t109 * pkin(4);
t278 = qJ(5) + pkin(8);
t181 = pkin(4) * t155 - qJ(5) * t243;
t223 = t157 * qJD(5);
t277 = -t155 * t223 + t181 * qJD(3) + (-t140 + (qJ(5) * t155 - t121) * t154) * qJD(4) + t300;
t244 = t155 * t157;
t276 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t244 + (-qJD(5) * t155 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t158) * t154 + t299;
t236 = qJD(1) * t151;
t119 = qJD(2) * pkin(7) + t156 * t236;
t234 = qJD(1) * t155;
t210 = t152 * t234;
t81 = t119 * t158 + t210;
t74 = qJD(3) * pkin(8) + t81;
t211 = t159 * t236;
t83 = t121 * qJD(2) - t211;
t34 = -t154 * t74 + t157 * t83;
t22 = -qJ(5) * t114 + t34;
t15 = -pkin(4) * t139 + t22;
t275 = -t22 + t15;
t227 = qJD(4) * t154;
t274 = (-t155 * t224 - t158 * t227) * pkin(7) + t299;
t273 = -t89 * qJD(4) + t300;
t202 = qJD(4) * t278;
t116 = t186 * qJD(2);
t107 = t155 * t119;
t80 = t158 * t235 - t107;
t51 = t154 * t116 + t157 * t80;
t272 = t223 - t51 + (qJ(5) * t231 - t202) * t154;
t50 = t157 * t116 - t154 * t80;
t271 = -t181 * qJD(2) - t154 * qJD(5) - t157 * t202 - t50;
t270 = qJD(2) * pkin(2);
t220 = qJDD(1) * t152;
t228 = qJD(3) * t158;
t194 = t119 * t228 + t302 * t155 - t158 * t220;
t256 = qJDD(3) * pkin(3);
t26 = t194 - t256;
t12 = t26 + t296;
t269 = t12 * t154;
t268 = t139 * t34;
t35 = t154 * t83 + t157 * t74;
t23 = -qJ(5) * t112 + t35;
t267 = t23 * t139;
t266 = t26 * t154;
t265 = t35 * t139;
t264 = t54 * qJ(5);
t263 = t54 * t154;
t262 = t55 * qJ(5);
t261 = t55 * t157;
t254 = t112 * t154;
t253 = t114 * t112;
t251 = t114 * t154;
t250 = t114 * t157;
t247 = t154 * t155;
t245 = t154 * t158;
t241 = qJDD(1) - g(3);
t149 = t155 ^ 2;
t150 = t158 ^ 2;
t238 = t149 - t150;
t237 = t149 + t150;
t232 = qJD(2) * t156;
t218 = t173 * t154;
t217 = -t155 * t220 - t302 * t158;
t161 = qJD(2) ^ 2;
t216 = t155 * t161 * t158;
t215 = pkin(7) + t283;
t214 = t151 * t232;
t213 = qJD(2) * t248;
t212 = t139 * t233;
t209 = g(3) * (pkin(2) * t248 + pkin(7) * t249);
t204 = t159 * t221;
t29 = -t119 * t229 - t217;
t25 = qJDD(3) * pkin(8) + t29;
t45 = qJD(2) * t117 + t121 * qJDD(2) + t185;
t6 = t154 * t45 + t157 * t25 + t83 * t225 - t74 * t227;
t195 = -qJD(5) - t284;
t192 = t112 * t211;
t191 = t114 * t211;
t190 = t155 * t205;
t188 = g(1) * t98 + g(2) * t96;
t184 = -t154 * t35 - t157 * t34;
t143 = pkin(4) * t157 + pkin(3);
t183 = t143 * t158 + t155 * t278;
t182 = qJDD(2) * t159 - t156 * t161;
t178 = -t101 * t157 + t154 * t248;
t177 = t154 * t109 - t139 * t225;
t176 = t157 * t109 + t139 * t227;
t175 = -g(1) * (t157 * t98 + t97 * t245) - g(2) * (t157 * t96 + t95 * t245) - g(3) * t84;
t174 = -g(1) * (t154 * t98 - t97 * t243) - g(2) * (t154 * t96 - t95 * t243) - g(3) * t85;
t172 = g(1) * t64 + g(2) * t62 + g(3) * t101;
t73 = -qJD(3) * pkin(3) - t80;
t171 = t154 * t228 + t155 * t225;
t170 = t173 - t26;
t169 = -g(3) * t248 + t189;
t168 = -pkin(8) * t109 - t139 * t73;
t167 = -g(1) * (-t154 * t97 - t157 * t64) - g(2) * (-t154 * t95 - t157 * t62) - g(3) * t178 - t6;
t165 = t173 - t194;
t7 = -qJD(4) * t35 - t154 * t25 + t157 * t45;
t120 = -t211 - t270;
t164 = -pkin(7) * qJDD(3) + (t120 + t211 - t270) * qJD(3);
t163 = t194 * t155 + t29 * t158 + (-t155 * t81 - t158 * t80) * qJD(3) - t188;
t162 = t7 + t295;
t123 = t278 * t157;
t122 = t278 * t154;
t118 = t215 * t155;
t111 = t157 * t121;
t108 = t112 ^ 2;
t92 = t97 * pkin(2);
t91 = t95 * pkin(2);
t88 = -pkin(7) * t245 + t111;
t82 = t171 * pkin(4) + pkin(7) * t228;
t70 = -qJ(5) * t247 + t89;
t67 = -t109 * t158 - t139 * t229;
t66 = t101 * qJD(3) + t155 * t213;
t65 = -t100 * qJD(3) + t158 * t213;
t59 = t210 + (qJD(2) * t283 + t119) * t158;
t58 = -qJ(5) * t244 + t111 + (-pkin(4) - t282) * t158;
t53 = -t108 + t292;
t52 = -t195 + t73;
t32 = -t252 - t55;
t31 = -t54 - t255;
t28 = (-t114 * t155 + t139 * t243) * qJD(2) + t177;
t27 = (t112 * t155 - t139 * t245) * qJD(2) + t176;
t20 = t68 * qJD(4) + t154 * t214 + t65 * t157;
t19 = t178 * qJD(4) - t65 * t154 + t157 * t214;
t17 = -t139 * t254 - t261;
t16 = -t139 * t250 - t263;
t14 = t171 * t112 + t55 * t247;
t13 = -t54 * t244 + (-t154 * t226 + t158 * t224) * t114;
t11 = (t139 * t230 + t55) * t158 + (-qJD(3) * t112 - t177) * t155;
t10 = (-t139 * t224 + t54) * t158 + (qJD(3) * t114 + t176) * t155;
t9 = -t297 * t154 + t298 * t157;
t8 = (-t112 * t157 - t251) * t228 + (t263 - t261 + (-t250 + t254) * qJD(4)) * t155;
t5 = -t100 * t54 + t109 * t178 + t114 * t66 + t139 * t20;
t4 = t100 * t55 + t109 * t68 + t112 * t66 - t139 * t19;
t3 = -qJD(5) * t112 - t262 + t6;
t2 = -t114 * qJD(5) + t264 + t279 + t7;
t1 = -t112 * t20 - t114 * t19 + t178 * t55 + t54 * t68;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t241, 0, 0, 0, 0, 0, 0, t182 * t151, (-qJDD(2) * t156 - t159 * t161) * t151, 0, -g(3) + (t152 ^ 2 + (t156 ^ 2 + t159 ^ 2) * t151 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t66 * qJD(3) - t100 * qJDD(3) + (-t155 * t204 + t182 * t158) * t151, -t65 * qJD(3) - t101 * qJDD(3) + (-t182 * t155 - t158 * t204) * t151, (t100 * t155 + t101 * t158) * qJDD(2) + (t155 * t66 + t158 * t65 + (t100 * t158 - t101 * t155) * qJD(3)) * qJD(2), t194 * t100 + t29 * t101 + t81 * t65 - t80 * t66 - g(3) + (t120 * t232 - t159 * t86) * t151, 0, 0, 0, 0, 0, 0, t4, t5, t1, t100 * t26 - t178 * t6 + t19 * t34 + t20 * t35 + t66 * t73 + t68 * t7 - g(3), 0, 0, 0, 0, 0, 0, t4, t5, t1, t100 * t12 + t15 * t19 - t178 * t3 + t2 * t68 + t20 * t23 + t52 * t66 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t241 * t248 + t189, -t241 * t249 + t188, 0, 0, qJDD(2) * t149 + 0.2e1 * t190, 0.2e1 * t155 * t146 - 0.2e1 * t238 * t221, qJDD(3) * t155 + t158 * t160, qJDD(2) * t150 - 0.2e1 * t190, qJDD(3) * t158 - t155 * t160, 0, t164 * t155 + t294 * t158, -t294 * t155 + t164 * t158, t237 * t257 + (-g(3) * t156 - t237 * t207) * t151 + t163, -t86 * pkin(2) + g(1) * t92 + g(2) * t91 - t209 + (-t120 * t156 + (t155 * t80 - t158 * t81) * t159) * t236 + t163 * pkin(7), t13, t8, t10, t14, t11, t67, t88 * t109 - t273 * t139 + (-t7 + (pkin(7) * t112 + t154 * t73) * qJD(3)) * t158 + (pkin(7) * t55 + qJD(3) * t34 + t225 * t73 - t192 + t266) * t155 + t174, -t89 * t109 + t274 * t139 + (t6 + (pkin(7) * t114 + t157 * t73) * qJD(3)) * t158 + (-pkin(7) * t54 - qJD(3) * t35 + t26 * t157 - t227 * t73 - t191) * t155 + t175, t54 * t88 - t55 * t89 - t273 * t114 - t274 * t112 + t184 * t228 + (-t154 * t6 - t157 * t7 + (t154 * t34 - t157 * t35) * qJD(4) + t169) * t155, t6 * t89 + t7 * t88 - g(1) * (-t187 * t97 - t92) - g(2) * (-t187 * t95 - t91) - t209 + t274 * t35 + t273 * t34 + (-g(3) * t187 - t234 * t73) * t248 + (t155 * t26 + t228 * t73 - t188) * pkin(7), t13, t8, t10, t14, t11, t67, t58 * t109 + t82 * t112 + t118 * t55 + (t230 * t52 - t2) * t158 - t277 * t139 + (qJD(3) * t15 + t225 * t52 - t192 + t269) * t155 + t174, -t70 * t109 + t82 * t114 - t118 * t54 + (t224 * t52 + t3) * t158 + t276 * t139 + (-qJD(3) * t23 + t12 * t157 - t227 * t52 - t191) * t155 + t175, t54 * t58 - t55 * t70 - t277 * t114 - t276 * t112 + (-t15 * t157 - t154 * t23) * t228 + (-t154 * t3 - t157 * t2 + (t15 * t154 - t157 * t23) * qJD(4) + t169) * t155, t3 * t70 + t2 * t58 + t12 * t118 + t52 * t82 - g(1) * (-t183 * t97 + t215 * t98 - t92) - g(2) * (-t183 * t95 + t215 * t96 - t91) - t209 + t276 * t23 + t277 * t15 + (-g(3) * pkin(4) * t246 + (-g(3) * t183 - t234 * t52) * t159) * t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216, t238 * t161, t219, t216, t146, qJDD(3), qJD(3) * t81 - t120 * t233 + t165, -t120 * t231 + (t80 + t107) * qJD(3) + t172 + t217, 0, 0, t16, t9, t28, t17, t27, t212, -t34 * t233 - pkin(3) * t55 - t81 * t112 + t50 * t139 + t168 * t154 + (pkin(8) * qJD(4) * t139 + t170) * t157, t35 * t233 + pkin(3) * t54 - t81 * t114 + t266 + (-pkin(8) * t227 - t51) * t139 + t168 * t157 - t218, t51 * t112 + t50 * t114 + (t6 + t268 + (qJD(4) * t114 - t55) * pkin(8)) * t157 + (-t7 + t265 + (qJD(4) * t112 - t54) * pkin(8)) * t154 - t172, -t34 * t50 - t35 * t51 - t73 * t81 + t170 * pkin(3) + (qJD(4) * t184 - t7 * t154 + t6 * t157 - t172) * pkin(8), t16, t9, t28, t17, t27, t212, -t15 * t233 - t122 * t109 - t59 * t112 - t143 * t55 - t271 * t139 + (-t52 * t231 + (t52 + t284) * qJD(4)) * t154 + (-t12 + t173) * t157, -t123 * t109 - t59 * t114 + t269 + t143 * t54 + t272 * t139 + (pkin(4) * t251 + t157 * t52) * qJD(4) + (t155 * t23 - t243 * t52) * qJD(2) - t218, -t122 * t54 - t123 * t55 - t271 * t114 - t272 * t112 + (t139 * t15 + t3) * t157 + (-t2 + t267) * t154 - t172, t3 * t123 - t2 * t122 - t12 * t143 - g(1) * (-t143 * t63 + t278 * t64) - g(2) * (-t143 * t61 + t278 * t62) - g(3) * (-t100 * t143 + t101 * t278) + (pkin(4) * t227 - t59) * t52 + t272 * t23 + t271 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t53, t31, -t253, t32, t109, -t73 * t114 + t162 - t265, t112 * t73 + t167 - t268, 0, 0, t253, t53, t31, -t253, t32, t109, 0.2e1 * t279 + t264 - t267 + (t195 - t52) * t114 + t162, -t292 * pkin(4) + t262 - t22 * t139 + (qJD(5) + t52) * t112 + t167, t54 * pkin(4) - t275 * t112, t275 * t23 + (-t52 * t114 + t2 + t295) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, t298, -t108 - t292, t112 * t23 + t114 * t15 - t165 - t256 + t296;];
tau_reg = t18;
