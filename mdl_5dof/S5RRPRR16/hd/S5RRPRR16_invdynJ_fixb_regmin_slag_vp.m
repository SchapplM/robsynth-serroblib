% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR16_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:28
% EndTime: 2019-12-31 20:47:40
% DurationCPUTime: 4.55s
% Computational Cost: add. (3481->449), mult. (8639->625), div. (0->0), fcn. (6627->10), ass. (0->233)
t156 = sin(pkin(5));
t160 = sin(qJ(2));
t247 = qJD(1) * qJD(2);
t231 = t160 * t247;
t217 = t156 * t231;
t164 = cos(qJ(2));
t245 = qJDD(1) * t164;
t311 = -t156 * t245 + t217;
t260 = qJD(1) * t160;
t232 = t156 * t260;
t157 = cos(pkin(5));
t249 = t157 * qJD(1);
t241 = pkin(1) * t249;
t96 = pkin(7) * t232 - t164 * t241;
t264 = qJD(3) + t96;
t144 = qJD(2) + t249;
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t259 = qJD(1) * t164;
t235 = t156 * t259;
t87 = t159 * t144 + t163 * t235;
t83 = qJD(5) + t87;
t153 = t156 ^ 2;
t243 = 0.2e1 * t153;
t123 = qJD(4) + t232;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t218 = t159 * t235;
t89 = t163 * t144 - t218;
t45 = -t162 * t123 + t158 * t89;
t310 = t123 * t45;
t165 = cos(qJ(1));
t266 = t165 * t164;
t161 = sin(qJ(1));
t270 = t161 * t160;
t101 = -t157 * t266 + t270;
t267 = t165 * t160;
t269 = t161 * t164;
t103 = t157 * t269 + t267;
t275 = t156 * t164;
t309 = g(1) * t103 + g(2) * t101 - g(3) * t275;
t265 = pkin(3) * t232 + t264;
t286 = qJ(3) * t160;
t209 = pkin(2) * t164 + t286;
t198 = -pkin(1) - t209;
t91 = t198 * t156;
t79 = qJD(1) * t91;
t308 = t79 * t232 + qJDD(3);
t272 = t159 * t165;
t284 = t101 * t163;
t276 = t156 * t161;
t65 = t103 * t163 - t159 * t276;
t99 = t157 * t159 + t163 * t275;
t185 = g(1) * t65 + g(2) * (t156 * t272 + t284) - g(3) * t99;
t166 = -pkin(2) - pkin(8);
t44 = t166 * t144 + t265;
t182 = t166 * t164 - pkin(1) - t286;
t72 = t182 * t156;
t53 = qJD(1) * t72;
t17 = t159 * t44 + t163 * t53;
t244 = t157 * qJDD(1);
t143 = qJDD(2) + t244;
t230 = t164 * t247;
t246 = qJDD(1) * t160;
t183 = t230 + t246;
t178 = t183 * t156;
t138 = pkin(7) * t235;
t277 = t156 * t160;
t146 = pkin(7) * t277;
t240 = pkin(1) * qJD(2) * t157;
t222 = qJD(1) * t240;
t239 = pkin(1) * t244;
t219 = -qJD(2) * t138 - qJDD(1) * t146 - t160 * t222 + t164 * t239;
t200 = qJDD(3) - t219;
t21 = pkin(3) * t178 + t166 * t143 + t200;
t118 = pkin(2) * t217;
t285 = qJ(3) * t164;
t208 = pkin(8) * t160 - t285;
t256 = qJD(3) * t160;
t172 = t208 * qJD(2) - t256;
t28 = t118 + (t172 * qJD(1) + t182 * qJDD(1)) * t156;
t171 = -t17 * qJD(4) - t159 * t28 + t163 * t21;
t93 = qJDD(4) + t178;
t4 = -t93 * pkin(4) - t171;
t307 = (pkin(4) * t89 + pkin(9) * t83) * t83 + t185 + t4;
t112 = t159 * pkin(4) - t163 * pkin(9) + qJ(3);
t216 = pkin(4) * t163 + pkin(9) * t159;
t35 = -qJD(4) * t218 + t159 * t143 + (qJD(4) * t144 - t311) * t163;
t33 = qJDD(5) + t35;
t306 = (t216 * qJD(4) - (-pkin(3) - t216) * t232 + t264) * t83 + t112 * t33;
t34 = -t87 * qJD(4) + t163 * t143 + t159 * t311;
t47 = t158 * t123 + t162 * t89;
t12 = t47 * qJD(5) + t158 * t34 - t162 * t93;
t102 = t157 * t267 + t269;
t274 = t156 * t165;
t194 = -t101 * t159 + t163 * t274;
t305 = t102 * t162 + t158 * t194;
t304 = -t102 * t158 + t162 * t194;
t236 = -pkin(1) * t164 - pkin(2);
t56 = pkin(3) * t277 + t146 + (-pkin(8) + t236) * t157;
t202 = t159 * t56 + t163 * t72;
t258 = qJD(2) * t160;
t234 = t156 * t258;
t137 = pkin(2) * t234;
t52 = t172 * t156 + t137;
t301 = pkin(1) * t160;
t149 = t157 * t301;
t302 = pkin(3) + pkin(7);
t77 = (t302 * t275 + t149) * qJD(2);
t303 = -t202 * qJD(4) - t159 * t52 + t163 * t77;
t14 = t123 * pkin(9) + t17;
t126 = t144 * qJ(3);
t97 = t160 * t241 + t138;
t76 = pkin(3) * t235 + t97;
t50 = t126 + t76;
t18 = t87 * pkin(4) - t89 * pkin(9) + t50;
t207 = t158 * t14 - t162 * t18;
t253 = qJD(4) * t163;
t255 = qJD(4) * t159;
t192 = t159 * t21 + t163 * t28 + t44 * t253 - t53 * t255;
t3 = t93 * pkin(9) + t192;
t124 = t143 * qJ(3);
t125 = t144 * qJD(3);
t220 = pkin(7) * t311 - t160 * t239 - t164 * t222;
t36 = -t124 - t125 + t220;
t23 = -pkin(3) * t311 - t36;
t6 = t35 * pkin(4) - t34 * pkin(9) + t23;
t1 = -t207 * qJD(5) + t158 * t6 + t162 * t3;
t300 = t143 * pkin(2);
t299 = t45 * t83;
t298 = t47 * t83;
t142 = pkin(2) * t232;
t261 = qJD(1) * t156;
t74 = t208 * t261 + t142;
t297 = t159 * t76 + t163 * t74;
t250 = qJD(5) * t162;
t251 = qJD(5) * t158;
t11 = t123 * t250 + t158 * t93 + t162 * t34 - t89 * t251;
t296 = t11 * t158;
t294 = t158 * t33;
t293 = t158 * t83;
t292 = t159 * t93;
t291 = t162 * t33;
t225 = t162 * t83;
t290 = t166 * t93;
t289 = t87 * t123;
t288 = t89 * t123;
t281 = t123 * t159;
t280 = t123 * t163;
t279 = t143 * t157;
t278 = t153 * qJD(1) ^ 2;
t273 = t158 * t166;
t271 = t160 * t162;
t268 = t162 * t166;
t263 = pkin(7) * t275 + t149;
t154 = t160 ^ 2;
t262 = -t164 ^ 2 + t154;
t257 = qJD(2) * t164;
t254 = qJD(4) * t162;
t252 = qJD(4) * t166;
t248 = qJD(2) - t144;
t242 = g(3) * t277;
t238 = t164 * t278;
t237 = t159 * t275;
t90 = -t157 * qJ(3) - t263;
t233 = t156 * t257;
t224 = t144 + t249;
t223 = t143 + t244;
t221 = t160 * t238;
t71 = pkin(3) * t275 - t90;
t214 = g(1) * t101 - g(2) * t103;
t104 = -t157 * t270 + t266;
t213 = -g(1) * t104 - g(2) * t102;
t212 = g(1) * t165 + g(2) * t161;
t8 = t162 * t14 + t158 * t18;
t27 = pkin(9) * t277 + t202;
t100 = t157 * t163 - t237;
t31 = t99 * pkin(4) - t100 * pkin(9) + t71;
t206 = t158 * t31 + t162 * t27;
t205 = -t158 * t27 + t162 * t31;
t16 = -t159 * t53 + t163 * t44;
t203 = -t159 * t72 + t163 * t56;
t201 = -0.2e1 * qJD(2) * t79;
t141 = t164 * t240;
t199 = -pkin(7) * t234 + t141;
t197 = -t83 * t250 - t294;
t196 = -t83 * t251 + t291;
t195 = -t100 * t158 + t156 * t271;
t64 = t100 * t162 + t158 * t277;
t187 = -qJ(3) * t257 - t256;
t37 = t118 + (t187 * qJD(1) + t198 * qJDD(1)) * t156;
t78 = t187 * t156 + t137;
t193 = qJD(1) * t78 + qJDD(1) * t91 + t37;
t191 = t159 * t77 + t163 * t52 + t56 * t253 - t72 * t255;
t190 = t83 * t123;
t189 = t123 * t47;
t98 = t263 * qJD(2);
t181 = -g(1) * t102 + g(2) * t104 + t98 * t144;
t151 = t157 * qJD(3);
t55 = -t302 * t234 + t141 + t151;
t179 = -t213 + t242;
t13 = -t123 * pkin(4) - t16;
t177 = -pkin(9) * t33 + (t13 + t16) * t83;
t176 = t213 - t220;
t175 = -t179 + t23;
t174 = t219 + t309;
t173 = (t248 * t259 + t246) * t156;
t2 = -t8 * qJD(5) - t158 * t3 + t162 * t6;
t170 = qJD(5) * t166 * t83 + t179;
t169 = t97 * t144 + t174;
t168 = (pkin(9) * t235 - qJD(5) * t112 + t297) * t83 + t309;
t95 = -qJ(3) * t235 + t142;
t92 = t236 * t157 + t146;
t86 = t163 * t93;
t82 = -t151 - t199;
t81 = (t158 * t164 + t159 * t271) * t261;
t80 = t158 * t159 * t232 - t162 * t235;
t73 = -t126 - t97;
t70 = -t144 * pkin(2) + t264;
t66 = t103 * t159 + t163 * t276;
t62 = -qJD(4) * t237 + t157 * t253 - t163 * t234;
t61 = -t99 * qJD(4) + t159 * t234;
t41 = t200 - t300;
t39 = t104 * t158 + t66 * t162;
t38 = t104 * t162 - t66 * t158;
t29 = -pkin(4) * t235 + t159 * t74 - t163 * t76;
t26 = -pkin(4) * t277 - t203;
t25 = t195 * qJD(5) + t158 * t233 + t61 * t162;
t24 = t64 * qJD(5) + t61 * t158 - t162 * t233;
t15 = t62 * pkin(4) - t61 * pkin(9) + t55;
t10 = -pkin(4) * t233 - t303;
t9 = pkin(9) * t233 + t191;
t5 = [qJDD(1), g(1) * t161 - g(2) * t165, t212, (qJDD(1) * t154 + 0.2e1 * t160 * t230) * t153, (t160 * t245 - t262 * t247) * t243, (t223 * t160 + t224 * t257) * t156, (t223 * t164 - t224 * t258) * t156, t279, -t146 * t143 + t219 * t157 + (t164 * t279 + (-t231 + t245) * t243) * pkin(1) - t181, -t183 * pkin(1) * t243 - t263 * t143 - t199 * t144 + t220 * t157 - t214, ((qJD(2) * t70 - qJDD(1) * t90 - t36 + (qJD(2) * t92 - t82) * qJD(1)) * t164 + (qJD(2) * t73 + qJDD(1) * t92 + t41 + (qJD(2) * t90 + t98) * qJD(1)) * t160 - t212) * t156, t92 * t143 + t41 * t157 + (t160 * t201 + t193 * t164) * t156 + t181, -t90 * t143 - t82 * t144 - t36 * t157 + (-t193 * t160 + t164 * t201) * t156 + t214, t37 * t91 + t79 * t78 + t36 * t90 + t73 * t82 + t41 * t92 + t70 * t98 - g(1) * (-t161 * pkin(1) - t102 * pkin(2) + pkin(7) * t274 - t101 * qJ(3)) - g(2) * (t165 * pkin(1) + t104 * pkin(2) + pkin(7) * t276 + t103 * qJ(3)), t34 * t100 + t89 * t61, -t100 * t35 - t34 * t99 - t61 * t87 - t89 * t62, t100 * t93 + t61 * t123 + (t160 * t34 + t89 * t257) * t156, -t62 * t123 - t99 * t93 + (-t160 * t35 - t87 * t257) * t156, (t123 * t257 + t160 * t93) * t156, t303 * t123 + t203 * t93 + t55 * t87 + t71 * t35 + t23 * t99 + t50 * t62 - g(1) * t194 - g(2) * t66 + (t16 * t257 + t160 * t171) * t156, -t191 * t123 - t202 * t93 + t55 * t89 + t71 * t34 + t23 * t100 + t50 * t61 + g(1) * t284 - g(2) * t65 + (g(1) * t272 - t160 * t192 - t17 * t257) * t156, t11 * t64 + t25 * t47, t11 * t195 - t12 * t64 - t24 * t47 - t25 * t45, t11 * t99 + t25 * t83 + t64 * t33 + t47 * t62, -t12 * t99 + t195 * t33 - t24 * t83 - t45 * t62, t33 * t99 + t83 * t62, (-qJD(5) * t206 + t162 * t15 - t158 * t9) * t83 + t205 * t33 + t2 * t99 - t207 * t62 + t10 * t45 + t26 * t12 - t4 * t195 + t13 * t24 - g(1) * t304 - g(2) * t39, -(qJD(5) * t205 + t158 * t15 + t162 * t9) * t83 - t206 * t33 - t1 * t99 - t8 * t62 + t10 * t47 + t26 * t11 + t4 * t64 + t13 * t25 + g(1) * t305 - g(2) * t38; 0, 0, 0, -t221, t262 * t278, t173, (-t248 * t260 + t245) * t156, t143, t278 * t301 + t169, pkin(1) * t238 - t96 * t144 - t176 + t242, ((-pkin(2) * t160 + t285) * qJDD(1) + ((-qJ(3) * qJD(2) - t73 - t97) * t160 + (-pkin(2) * qJD(2) + t264 - t70) * t164) * qJD(1)) * t156, -t95 * t235 - t169 - 0.2e1 * t300 + t308, 0.2e1 * t124 + t125 + t264 * t144 + (-g(3) * t160 + (t160 * t95 + t164 * t79) * qJD(1)) * t156 + t176, -t36 * qJ(3) - t41 * pkin(2) - t79 * t95 - t70 * t97 - g(1) * (-t103 * pkin(2) + t104 * qJ(3)) - g(2) * (-t101 * pkin(2) + t102 * qJ(3)) - t264 * t73 - g(3) * t209 * t156, t34 * t163 - t281 * t89, (-t35 - t288) * t163 + (-t34 + t289) * t159, -t123 * t255 + t86 + (-t160 * t281 - t164 * t89) * t261, -t123 * t253 - t292 + (-t160 * t280 + t164 * t87) * t261, -t123 * t235, -t16 * t235 + qJ(3) * t35 + t265 * t87 + (t290 + (t50 - t76) * t123) * t163 + ((t74 - t252) * t123 + t175) * t159, qJ(3) * t34 + t297 * t123 + t17 * t235 + t265 * t89 + (-t123 * t50 - t290) * t159 + (-t123 * t252 + t175) * t163, t11 * t162 * t163 + (-t159 * t254 - t163 * t251 - t81) * t47, t81 * t45 + t47 * t80 + (t158 * t47 + t162 * t45) * t255 + (-t296 - t12 * t162 + (t158 * t45 - t162 * t47) * qJD(5)) * t163, -t81 * t83 + (-t254 * t83 + t11) * t159 + (t189 + t196) * t163, t80 * t83 + (qJD(4) * t293 - t12) * t159 + (t197 - t310) * t163, t33 * t159 + t280 * t83, -t13 * t80 - t29 * t45 + t306 * t162 + t168 * t158 + (-t33 * t273 + t2 + (-t13 * t158 + t166 * t45) * qJD(4) - t170 * t162) * t159 + (-t207 * t232 + t13 * t250 - t166 * t12 + t4 * t158 + (-t273 * t83 - t207) * qJD(4)) * t163, -t13 * t81 - t29 * t47 - t306 * t158 + t168 * t162 + (-t33 * t268 - t1 + (-t13 * t162 + t166 * t47) * qJD(4) + t170 * t158) * t159 + (-t8 * t232 - t13 * t251 - t166 * t11 + t4 * t162 + (-t268 * t83 - t8) * qJD(4)) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t143 + t221, -t144 ^ 2 - t154 * t278, t73 * t144 - t174 - t300 + t308, 0, 0, 0, 0, 0, -t123 * t281 - t144 * t87 + t86, -t123 ^ 2 * t163 - t144 * t89 - t292, 0, 0, 0, 0, 0, -t144 * t225 + (-t158 * t190 - t12) * t163 + (t197 + t310) * t159, t144 * t293 + (-t162 * t190 - t11) * t163 + (t189 - t196) * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t87, -t87 ^ 2 + t89 ^ 2, t34 + t289, t288 - t35, t93, t17 * t123 - t50 * t89 + t171 - t185, g(1) * t66 - g(2) * t194 + g(3) * t100 + t16 * t123 + t50 * t87 - t192, t225 * t47 + t296, (t11 - t299) * t162 + (-t12 - t298) * t158, t225 * t83 - t47 * t89 + t294, -t158 * t83 ^ 2 + t45 * t89 + t291, -t83 * t89, -pkin(4) * t12 + t177 * t158 - t307 * t162 - t17 * t45 + t207 * t89, -pkin(4) * t11 + t307 * t158 + t177 * t162 - t17 * t47 + t8 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t45, -t45 ^ 2 + t47 ^ 2, t11 + t299, -t12 + t298, t33, -g(1) * t38 - g(2) * t305 - g(3) * t195 - t13 * t47 + t8 * t83 + t2, g(1) * t39 - g(2) * t304 + g(3) * t64 + t13 * t45 - t207 * t83 - t1;];
tau_reg = t5;
