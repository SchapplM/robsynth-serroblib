% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:51:07
% EndTime: 2019-12-05 18:51:21
% DurationCPUTime: 5.90s
% Computational Cost: add. (7276->448), mult. (17101->602), div. (0->0), fcn. (13764->14), ass. (0->222)
t166 = cos(qJ(5));
t249 = qJD(5) * t166;
t163 = sin(qJ(3));
t167 = cos(qJ(3));
t168 = cos(qJ(2));
t254 = qJD(1) * t168;
t164 = sin(qJ(2));
t255 = qJD(1) * t164;
t105 = t163 * t255 - t167 * t254;
t107 = t163 * t168 + t164 * t167;
t106 = t107 * qJD(1);
t162 = sin(qJ(4));
t297 = cos(qJ(4));
t66 = t297 * t105 + t106 * t162;
t313 = t166 * t66;
t321 = t249 - t313;
t198 = -t162 * t105 + t297 * t106;
t158 = qJ(2) + qJ(3);
t150 = qJ(4) + t158;
t140 = cos(t150);
t136 = g(3) * t140;
t139 = sin(t150);
t169 = cos(qJ(1));
t265 = t139 * t169;
t165 = sin(qJ(1));
t266 = t139 * t165;
t301 = g(1) * t265 + g(2) * t266 + t136;
t123 = qJD(1) * pkin(1) + pkin(2) * t254;
t80 = -pkin(3) * t105 + t123;
t320 = t80 * t198 + t301;
t161 = sin(qJ(5));
t236 = qJD(4) * t297;
t251 = qJD(4) * t162;
t203 = t163 * t164 - t167 * t168;
t191 = t203 * qJD(3);
t174 = t203 * qJD(2) + t191;
t55 = t174 * qJD(1) - t107 * qJDD(1);
t192 = qJD(3) * t107;
t175 = t107 * qJD(2) + t192;
t56 = t175 * qJD(1) + t203 * qJDD(1);
t195 = t105 * t236 + t106 * t251 + t162 * t56 + t297 * t55;
t155 = qJD(2) + qJD(3);
t230 = qJD(4) + t155;
t205 = t166 * t230;
t154 = qJDD(2) + qJDD(3);
t220 = qJDD(4) + t154;
t250 = qJD(5) * t161;
t16 = -qJD(5) * t205 - t161 * t220 - t166 * t195 - t198 * t250;
t59 = t161 * t230 - t166 * t198;
t271 = qJD(5) * t59;
t17 = t161 * t195 - t166 * t220 + t271;
t309 = qJD(5) - t66;
t317 = t161 * t309;
t57 = -t161 * t198 - t205;
t1 = -t16 * t166 - t161 * t17 - t317 * t59 - t321 * t57;
t231 = t162 * t55 - t297 * t56;
t28 = -qJD(4) * t198 + t231;
t26 = qJDD(5) + t28;
t7 = t166 * t26 - t198 * t57 - t309 * t317;
t35 = -pkin(4) * t66 + pkin(6) * t198 + t80;
t279 = qJD(2) * pkin(2);
t111 = pkin(3) * t155 + t167 * t279;
t218 = t297 * t279;
t204 = t163 * t218;
t82 = t162 * t111 + t204;
t79 = t230 * pkin(6) + t82;
t29 = -t161 * t79 + t166 * t35;
t319 = t29 * t309;
t30 = t161 * t35 + t166 * t79;
t318 = t309 * t30;
t148 = sin(t158);
t149 = cos(t158);
t304 = g(1) * t169 + g(2) * t165;
t181 = g(3) * t149 + t148 * t304;
t215 = t163 * t236;
t247 = qJDD(2) * t163;
t252 = qJD(3) * t167;
t296 = pkin(2) * t167;
t145 = qJDD(2) * t296;
t243 = t163 * t279;
t88 = pkin(3) * t154 - qJD(3) * t243 + t145;
t85 = t297 * t88;
t49 = -(qJD(2) * (t162 * t252 + t215) + t162 * t247) * pkin(2) - t111 * t251 + t85;
t316 = t49 + t320;
t13 = t16 * t161;
t9 = t321 * t59 - t13;
t282 = t161 * t26 + t249 * t309;
t8 = t198 * t59 - t309 * t313 + t282;
t221 = t162 * t243;
t81 = t297 * t111 - t221;
t78 = -t230 * pkin(4) - t81;
t314 = t66 * t78;
t289 = t309 * t198;
t312 = t66 * t198;
t170 = qJD(2) ^ 2;
t310 = t168 * t170;
t22 = t198 ^ 2 - t66 ^ 2;
t135 = g(3) * t139;
t263 = t140 * t169;
t264 = t140 * t165;
t238 = -g(1) * t263 - g(2) * t264 + t135;
t237 = t297 * t163;
t137 = pkin(2) * t237;
t48 = -qJD(4) * t221 + qJDD(2) * t137 + t111 * t236 + t162 * t88 + t218 * t252;
t178 = -t80 * t66 - t238 - t48;
t18 = -t66 * t230 + t195;
t199 = t301 * t166 + t198 * t29 + t78 * t250;
t44 = -t220 * pkin(4) - t49;
t216 = t44 * t161 - t30 * t198 + t78 * t249;
t19 = -t198 * t155 - t231;
t41 = -pkin(4) * t198 - pkin(6) * t66;
t244 = pkin(2) * t255;
t295 = pkin(3) * t106;
t39 = -t295 + t41;
t36 = t39 - t244;
t143 = pkin(3) + t296;
t102 = t162 * t143 + t137;
t96 = pkin(6) + t102;
t229 = qJD(5) * t96 + t36;
t306 = t229 * t309;
t45 = -t105 * t155 + t55;
t46 = -t106 * t155 + t56;
t210 = g(1) * t165 - g(2) * t169;
t305 = t210 * t139;
t303 = -t161 * t29 + t166 * t30;
t302 = t163 ^ 2 + t167 ^ 2;
t180 = -t139 * t304 - t136;
t292 = g(3) * t168;
t291 = t164 * pkin(2);
t43 = t220 * pkin(6) + t48;
t159 = qJDD(1) * pkin(1);
t248 = qJD(1) * qJD(2);
t235 = t164 * t248;
t246 = t168 * qJDD(1);
t100 = t159 + (-t235 + t246) * pkin(2);
t47 = -pkin(3) * t56 + t100;
t6 = pkin(4) * t28 - pkin(6) * t195 + t47;
t3 = qJD(5) * t29 + t161 * t6 + t166 * t43;
t2 = t3 * t166;
t290 = t59 * t57;
t151 = t168 * pkin(2);
t144 = t151 + pkin(1);
t277 = t161 * t30;
t276 = t161 * t57;
t15 = t17 * t166;
t270 = qJD(5) * t309;
t269 = t105 * t106;
t262 = t161 * t165;
t261 = t161 * t169;
t260 = t162 * t163;
t259 = t165 * t166;
t258 = t166 * t169;
t156 = t164 ^ 2;
t157 = t168 ^ 2;
t256 = t156 - t157;
t141 = pkin(3) * t162 + pkin(6);
t242 = t141 * t270;
t77 = -t297 * t107 + t162 * t203;
t241 = t77 * t250;
t240 = t77 * t249;
t171 = qJD(1) ^ 2;
t239 = t164 * t171 * t168;
t233 = pkin(3) * t149 + t151;
t112 = -pkin(3) * t148 - t291;
t232 = -pkin(4) * t139 + t112;
t225 = -0.2e1 * pkin(1) * t248;
t224 = t2 + t238;
t223 = qJD(2) * (-qJD(3) + t155);
t222 = qJD(3) * (-qJD(2) - t155);
t217 = pkin(3) * t236;
t214 = t168 * t235;
t197 = t162 * t167 + t237;
t98 = t197 * t279;
t213 = pkin(3) * t251 - t98;
t212 = pkin(4) * t140 + pkin(6) * t139;
t186 = t297 * t203;
t37 = -qJD(4) * t186 - t107 * t251 - t162 * t175 - t297 * t174;
t208 = t26 * t77 - t309 * t37;
t207 = -t198 * t82 + t81 * t66;
t206 = t166 * t29 + t277;
t202 = t66 * t277 + t29 * t313 + t224;
t200 = -qJD(5) * t35 - t135 - t43;
t196 = t297 * t167 - t260;
t194 = t123 * t106 + t145 + t181;
t193 = t123 * t107;
t190 = pkin(1) * t171 + t304;
t69 = t143 * t236 + (qJD(3) * t196 - t163 * t251) * pkin(2);
t189 = -t26 * t96 - t309 * t69 - t314;
t188 = -g(3) * t148 - t123 * t105 + t304 * t149;
t101 = -pkin(2) * t260 + t297 * t143;
t187 = t210 + 0.2e1 * t159;
t38 = t77 * qJD(4) + t162 * t174 - t297 * t175;
t68 = -pkin(3) * t192 + (-t107 * pkin(3) - t291) * qJD(2);
t11 = t38 * pkin(4) + t37 * pkin(6) + t68;
t76 = -t107 * t162 - t186;
t87 = -t203 * pkin(3) + t144;
t40 = t76 * pkin(4) - t77 * pkin(6) + t87;
t184 = qJD(5) * t77 * t78 + t11 * t309 + t26 * t40;
t183 = -t40 * t270 - t37 * t78 + t44 * t77;
t179 = -t141 * t26 - t217 * t309 - t314;
t5 = t166 * t6;
t4 = -qJD(5) * t30 - t161 * t43 + t5;
t177 = -t206 * qJD(5) - t4 * t161 + t2;
t142 = -t297 * pkin(3) - pkin(4);
t114 = pkin(6) * t263;
t113 = pkin(6) * t264;
t110 = pkin(1) + t233;
t99 = t196 * t279;
t95 = -pkin(4) - t101;
t94 = t140 * t258 - t262;
t93 = -t140 * t261 - t259;
t92 = -t140 * t259 - t261;
t91 = t140 * t262 - t258;
t83 = -t244 - t295;
t70 = t143 * t251 + (t197 * qJD(3) + t215) * pkin(2);
t60 = -t105 ^ 2 + t106 ^ 2;
t34 = t161 * t41 + t166 * t81;
t33 = -t161 * t81 + t166 * t41;
t32 = t161 * t39 + t166 * t99;
t31 = -t161 * t99 + t166 * t39;
t10 = t317 * t57 - t15;
t12 = [0, 0, 0, 0, 0, qJDD(1), t210, t304, 0, 0, qJDD(1) * t156 + 0.2e1 * t214, 0.2e1 * t164 * t246 - 0.2e1 * t256 * t248, -qJDD(2) * t164 - t310, qJDD(1) * t157 - 0.2e1 * t214, -qJDD(2) * t168 + t164 * t170, 0, t164 * t225 + t168 * t187, -t164 * t187 + t168 * t225, t304, (t210 + t159) * pkin(1), -t106 * t174 - t55 * t107, -t56 * t107 + t203 * t55 + t155 * (t105 * t203 - t107 * t106), -t154 * t107 + t155 * t174, t105 * t175 + t203 * t56, t154 * t203 + t155 * t175, 0, -t144 * t56 - t100 * t203 + t210 * t149 - qJD(3) * t193 + (t105 * t291 - t193) * qJD(2), -t100 * t107 + t144 * t55 - t210 * t148 + t123 * t191 + (t106 * t291 + t123 * t203) * qJD(2), ((t167 * t107 + t163 * t203) * qJDD(2) + t302 * t310) * pkin(2) + t304, -t123 * t164 * t279 + (t100 + t210) * t144, t195 * t77 + t198 * t37, -t195 * t76 + t198 * t38 - t28 * t77 - t37 * t66, t220 * t77 - t230 * t37, t28 * t76 - t38 * t66, -t220 * t76 - t230 * t38, 0, t140 * t210 + t87 * t28 + t80 * t38 + t47 * t76 - t66 * t68, t195 * t87 - t198 * t68 - t80 * t37 + t47 * t77 - t305, t37 * t81 - t38 * t82 - t48 * t76 - t49 * t77 + t304, t110 * t210 + t47 * t87 + t80 * t68, -t59 * t241 + (-t16 * t77 - t37 * t59) * t166, (t161 * t59 + t166 * t57) * t37 + (t13 - t15 + (-t166 * t59 + t276) * qJD(5)) * t77, -t16 * t76 + t166 * t208 - t241 * t309 + t59 * t38, t57 * t240 + (t17 * t77 - t37 * t57) * t161, -t161 * t208 - t17 * t76 - t240 * t309 - t57 * t38, t26 * t76 + t309 * t38, -g(1) * t92 - g(2) * t94 + t161 * t183 + t166 * t184 + t29 * t38 + t4 * t76, -g(1) * t91 - g(2) * t93 - t161 * t184 + t166 * t183 - t3 * t76 - t30 * t38, t305 + (-t11 * t59 + t16 * t40 + t29 * t37 - t4 * t77 + (-t30 * t77 - t40 * t57) * qJD(5)) * t166 + (-t11 * t57 - t17 * t40 - t3 * t77 + t30 * t37 + (t29 * t77 + t40 * t59) * qJD(5)) * t161, t206 * t11 + (t303 * qJD(5) + t3 * t161 + t4 * t166) * t40 + t210 * (t110 + t212); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239, t256 * t171, -t164 * qJDD(1), t239, -t246, qJDD(2), t164 * t190 + t292, -g(3) * t164 + t168 * t190, 0, 0, t269, t60, t45, -t269, t46, t154, (-t105 * t255 + t167 * t154 + t163 * t222) * pkin(2) + t194, (-t106 * t255 + (-qJDD(2) - t154) * t163 + t167 * t222) * pkin(2) + t188, (t46 * t163 - t45 * t167) * pkin(2), (t292 + (qJD(1) * t123 + t304) * t164 + t302 * qJDD(2) * pkin(2)) * pkin(2), t312, t22, t18, -t312, t19, t220, t101 * t220 - t230 * t70 + t66 * t83 + t316, -t102 * t220 + t198 * t83 - t230 * t69 + t178, -t101 * t195 - t102 * t28 - t198 * t70 + t66 * t69 + t207, g(3) * t233 + t49 * t101 + t48 * t102 - t112 * t304 + t82 * t69 - t81 * t70 - t80 * t83, t9, t1, t8, t10, t7, t289, t17 * t95 + t57 * t70 + (-t44 - t306) * t166 + t189 * t161 + t199, -t16 * t95 + t59 * t70 + t189 * t166 + (t180 + t306) * t161 + t216, (-t17 * t96 + t36 * t59 - t57 * t69 + (t59 * t96 - t29) * qJD(5)) * t166 + (-t16 * t96 + t36 * t57 + t59 * t69 - t4 + (t57 * t96 - t30) * qJD(5)) * t161 + t202, t44 * t95 + t78 * t70 - g(1) * (t169 * t232 + t114) - g(2) * (t165 * t232 + t113) - g(3) * (-t212 - t233) + (-t229 * t29 + t3 * t96 + t30 * t69) * t166 + (-t229 * t30 - t29 * t69 - t4 * t96) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, t60, t45, -t269, t46, t154, pkin(2) * t163 * t223 + t194, (t167 * t223 - t247) * pkin(2) + t188, 0, 0, t312, t22, t18, -t312, t19, t220, -qJD(4) * t204 + t85 + t98 * t230 + (-t106 * t66 + t297 * t220) * pkin(3) + ((-pkin(3) * t230 - t111) * qJD(4) + (-qJD(2) * t252 - t247) * pkin(2)) * t162 + t320, t99 * t230 + (-t106 * t198 - t162 * t220 - t230 * t236) * pkin(3) + t178, -t99 * t66 + t98 * t198 + (-t297 * t195 - t162 * t28 + (-t162 * t198 + t297 * t66) * qJD(4)) * pkin(3) + t207, t81 * t98 - t82 * t99 + (t297 * t49 + t106 * t80 + t162 * t48 + (-t162 * t81 + t297 * t82) * qJD(4) + t181) * pkin(3), t9, t1, t8, t10, t7, t289, t142 * t17 - t31 * t309 + t213 * t57 + (-t44 - t242) * t166 + t179 * t161 + t199, -t142 * t16 + t32 * t309 + t213 * t59 + t179 * t166 + (t180 + t242) * t161 + t216, t31 * t59 + t32 * t57 + (-t57 * t217 - t141 * t17 + (t141 * t59 - t29) * qJD(5)) * t166 + (t59 * t217 - t141 * t16 - t4 + (t141 * t57 - t30) * qJD(5)) * t161 + t202, t44 * t142 - t30 * t32 - t29 * t31 - t78 * t98 - g(1) * (-pkin(4) * t265 + t114) - g(2) * (-pkin(4) * t266 + t113) + g(3) * t212 + ((t162 * t78 + t303 * t297) * qJD(4) + t181) * pkin(3) + t177 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, t22, t18, -t312, t19, t220, t230 * t82 + t316, t230 * t81 + t178, 0, 0, t9, t1, t8, t276 * t309 - t15, t7, t289, -pkin(4) * t17 - pkin(6) * t282 - t161 * t314 - t44 * t166 - t309 * t33 - t57 * t82 + t199, pkin(4) * t16 + t34 * t309 - t59 * t82 + (-pkin(6) * t26 - t314) * t166 + (pkin(6) * t270 + t180) * t161 + t216, t33 * t59 + t34 * t57 + (-t319 + (-t17 + t271) * pkin(6)) * t166 + (-t4 - t318 + (qJD(5) * t57 - t16) * pkin(6)) * t161 + t224, -g(1) * t114 - g(2) * t113 - t29 * t33 - t30 * t34 - t78 * t82 + (-t180 - t44) * pkin(4) + (t177 + t135) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, -t57 ^ 2 + t59 ^ 2, t309 * t57 - t16, -t290, t309 * t59 - t17, t26, -g(1) * t93 + g(2) * t91 + t161 * t200 - t249 * t79 - t59 * t78 + t318 + t5, g(1) * t94 - g(2) * t92 + t319 + t57 * t78 + (qJD(5) * t79 - t6) * t161 + t200 * t166, 0, 0;];
tau_reg = t12;
