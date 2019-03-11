% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:03
% EndTime: 2019-03-09 09:05:13
% DurationCPUTime: 3.25s
% Computational Cost: add. (5469->412), mult. (16703->559), div. (0->0), fcn. (13419->10), ass. (0->215)
t170 = sin(pkin(11));
t172 = cos(pkin(11));
t179 = cos(qJ(2));
t171 = sin(pkin(6));
t250 = qJD(1) * t171;
t232 = t179 * t250;
t216 = t172 * t232;
t176 = sin(qJ(2));
t233 = t176 * t250;
t137 = t170 * t233 - t216;
t173 = cos(pkin(6));
t249 = qJD(1) * t173;
t160 = qJD(2) + t249;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t101 = -t178 * t137 + t160 * t175;
t100 = qJD(6) + t101;
t248 = qJD(2) * t171;
t231 = t176 * t248;
t288 = pkin(2) * t231;
t200 = t170 * t179 + t172 * t176;
t141 = t200 * t250;
t241 = qJD(5) + t141;
t103 = t137 * t175 + t160 * t178;
t144 = t200 * t171;
t132 = qJD(1) * t144 + qJD(5);
t174 = sin(qJ(6));
t177 = cos(qJ(6));
t72 = t103 * t174 - t177 * t132;
t287 = t241 * t72;
t220 = t178 * t241;
t279 = pkin(1) * t176;
t239 = t173 * t279;
t256 = t171 * t179;
t274 = pkin(8) + qJ(3);
t134 = t274 * t256 + t239;
t121 = t134 * qJD(1);
t108 = t170 * t121;
t278 = pkin(1) * t179;
t238 = t173 * t278;
t158 = qJD(1) * t238;
t229 = t274 * t176;
t215 = t171 * t229;
t120 = -qJD(1) * t215 + t158;
t79 = t120 * t172 - t108;
t252 = qJD(4) - t79;
t133 = t141 ^ 2;
t286 = -t137 ^ 2 - t133;
t285 = -qJD(5) + t132;
t213 = qJD(1) * t231;
t131 = qJD(2) * t216 - t170 * t213;
t104 = pkin(2) * t160 + t120;
t68 = t104 * t172 - t108;
t209 = qJD(4) - t68;
t275 = pkin(4) * t141;
t280 = pkin(3) + pkin(9);
t35 = -t280 * t160 + t209 + t275;
t212 = (-pkin(2) * t179 - pkin(1)) * t171;
t199 = qJD(1) * t212;
t145 = qJD(3) + t199;
t182 = -qJ(4) * t141 + t145;
t50 = t280 * t137 + t182;
t15 = t175 * t35 + t178 * t50;
t152 = qJD(2) * t158;
t183 = (-qJD(2) * t229 + qJD(3) * t179) * t171;
t94 = qJD(1) * t183 + t152;
t257 = t171 * t176;
t106 = -qJD(2) * t134 - qJD(3) * t257;
t95 = t106 * qJD(1);
t47 = t170 * t94 - t172 * t95;
t33 = pkin(4) * t131 + t47;
t139 = qJD(2) * t144;
t130 = qJD(1) * t139;
t151 = pkin(2) * t213;
t223 = -qJ(4) * t131 + t151;
t192 = -qJD(4) * t141 + t223;
t34 = t280 * t130 + t192;
t225 = t175 * t34 - t178 * t33;
t4 = -pkin(5) * t131 + qJD(5) * t15 + t225;
t284 = t100 * (pkin(5) * t103 + t100 * pkin(10)) + t4;
t164 = -pkin(2) * t172 - pkin(3);
t162 = -pkin(9) + t164;
t276 = pkin(4) * t137;
t255 = t172 * t121;
t69 = t170 * t104 + t255;
t59 = -t160 * qJ(4) - t69;
t39 = -t59 - t276;
t283 = -t162 * t131 - t241 * t39;
t245 = qJD(5) * t178;
t247 = qJD(5) * t175;
t70 = t175 * t130 + t137 * t245 - t160 * t247;
t74 = t103 * t177 + t132 * t174;
t23 = qJD(6) * t74 - t177 * t131 + t174 * t70;
t244 = qJD(6) * t174;
t246 = qJD(5) * t177;
t260 = t141 * t175;
t89 = -t137 * t174 + t177 * t260;
t185 = t175 * t246 + t178 * t244 + t89;
t254 = t177 * t178;
t117 = t178 * t130;
t71 = qJD(5) * t103 - t117;
t282 = -t100 * t185 + t71 * t254;
t193 = -t241 * t132 * t175 + t178 * t131;
t48 = t170 * t95 + t172 * t94;
t236 = t160 * qJD(4) + t48;
t28 = -pkin(4) * t130 + t236;
t10 = pkin(5) * t71 - pkin(10) * t70 + t28;
t12 = pkin(10) * t132 + t15;
t20 = pkin(5) * t101 - pkin(10) * t103 + t39;
t208 = t12 * t174 - t177 * t20;
t195 = t175 * t33 + t178 * t34 + t35 * t245 - t50 * t247;
t3 = pkin(10) * t131 + t195;
t1 = -qJD(6) * t208 + t10 * t174 + t177 * t3;
t277 = pkin(3) * t130;
t78 = t120 * t170 + t255;
t52 = t78 - t276;
t222 = pkin(2) * t233 + qJ(4) * t137;
t57 = t280 * t141 + t222;
t273 = t175 * t52 + t178 * t57;
t119 = (pkin(2) + t278) * t173 - t215;
t82 = t119 * t172 - t170 * t134;
t51 = pkin(4) * t144 - t280 * t173 - t82;
t143 = t170 * t257 - t172 * t256;
t188 = -qJ(4) * t144 + t212;
t65 = t280 * t143 + t188;
t202 = t175 * t51 + t178 * t65;
t159 = qJD(2) * t238;
t105 = t159 + t183;
t63 = t172 * t105 + t170 * t106;
t272 = t100 * t72;
t271 = t100 * t74;
t243 = qJD(6) * t177;
t22 = -t103 * t244 + t174 * t131 + t132 * t243 + t177 * t70;
t270 = t174 * t22;
t269 = t174 * t71;
t268 = t175 * t23;
t267 = t177 * t71;
t211 = pkin(5) * t178 + pkin(10) * t175;
t266 = qJD(5) * t211 - (-pkin(4) - t211) * t141 + t252;
t265 = t100 * t174;
t218 = t100 * t177;
t264 = t101 * t137;
t263 = t103 * t137;
t261 = t137 * t160;
t167 = t171 ^ 2;
t180 = qJD(1) ^ 2;
t258 = t167 * t180;
t253 = t275 + t252;
t83 = t170 * t119 + t172 * t134;
t251 = t176 ^ 2 - t179 ^ 2;
t242 = qJD(2) - t160;
t240 = t22 * t175 + (t141 * t178 + t245) * t74;
t237 = t167 * t279;
t56 = -t173 * qJD(4) - t63;
t234 = t179 * t258;
t76 = -t173 * qJ(4) - t83;
t230 = t179 * t248;
t228 = qJD(1) * qJD(2) * t167;
t163 = pkin(2) * t170 + qJ(4);
t88 = t177 * t137 + t174 * t260;
t226 = t100 * t88 + t247 * t265;
t62 = t105 * t170 - t172 * t106;
t146 = pkin(5) * t175 - pkin(10) * t178 + t163;
t221 = -pkin(10) * t137 - qJD(6) * t146 + t273;
t219 = t241 * t103;
t217 = t160 + t249;
t214 = t179 * t228;
t6 = t12 * t177 + t174 * t20;
t207 = t141 * t62 + t144 * t47;
t19 = pkin(10) * t144 + t202;
t112 = t143 * t175 + t173 * t178;
t201 = t143 * t178 - t173 * t175;
t55 = -pkin(4) * t143 - t76;
t24 = -pkin(5) * t201 - pkin(10) * t112 + t55;
t206 = t174 * t24 + t177 * t19;
t205 = -t174 * t19 + t177 * t24;
t14 = -t175 * t50 + t178 * t35;
t140 = -t170 * t231 + t172 * t230;
t191 = -qJ(4) * t140 - qJD(4) * t144 + t288;
t40 = t280 * t139 + t191;
t41 = pkin(4) * t140 + t62;
t204 = -t175 * t40 + t178 * t41;
t203 = -t175 * t65 + t178 * t51;
t36 = -pkin(4) * t139 - t56;
t86 = t112 * t177 + t144 * t174;
t85 = t112 * t174 - t144 * t177;
t75 = pkin(3) * t137 + t182;
t198 = t141 * t75 + t47;
t197 = -t100 * t243 - t269;
t196 = t100 * t244 - t267;
t194 = t175 * t41 + t178 * t40 + t51 * t245 - t65 * t247;
t190 = -pkin(8) * t256 - t239;
t189 = -pkin(8) * t213 + t152;
t187 = -t175 * t131 - t132 * t220;
t186 = t190 * t160;
t11 = -pkin(5) * t132 - t14;
t184 = -pkin(10) * t71 + (t11 + t14) * t100;
t2 = -qJD(6) * t6 + t177 * t10 - t174 * t3;
t181 = t197 + t287;
t87 = pkin(3) * t143 + t188;
t84 = pkin(3) * t141 + t222;
t81 = qJD(5) * t112 - t139 * t178;
t80 = qJD(5) * t201 + t139 * t175;
t77 = -pkin(3) * t173 - t82;
t61 = pkin(3) * t139 + t191;
t58 = -pkin(3) * t160 + t209;
t49 = t192 + t277;
t26 = -qJD(6) * t85 + t140 * t174 + t177 * t80;
t25 = qJD(6) * t86 - t140 * t177 + t174 * t80;
t18 = -pkin(5) * t144 - t203;
t16 = pkin(5) * t137 + t175 * t57 - t178 * t52;
t13 = pkin(5) * t81 - pkin(10) * t80 + t36;
t8 = -pkin(5) * t140 + qJD(5) * t202 - t204;
t7 = pkin(10) * t140 + t194;
t5 = [0, 0, 0, 0.2e1 * t176 * t214, -0.2e1 * t251 * t228, t217 * t230, -t217 * t231, 0 (t186 + (t173 * t190 - 0.2e1 * t237) * qJD(1)) * qJD(2), -0.2e1 * pkin(1) * t214 - (-pkin(8) * t231 + t159) * t160 - t189 * t173, -t130 * t83 - t131 * t82 - t137 * t63 - t139 * t69 - t140 * t68 - t143 * t48 + t207, -t47 * t82 + t48 * t83 - t68 * t62 + t69 * t63 + (t145 + t199) * t288, t130 * t76 + t131 * t77 + t137 * t56 + t139 * t59 + t140 * t58 - t143 * t236 + t207, -t130 * t87 - t137 * t61 - t139 * t75 - t143 * t49 + t160 * t62 + t173 * t47, -t131 * t87 - t140 * t75 - t141 * t61 - t144 * t49 - t160 * t56 + t173 * t236, -t236 * t76 + t47 * t77 + t49 * t87 + t56 * t59 + t58 * t62 + t61 * t75, t103 * t80 + t112 * t70, -t101 * t80 - t103 * t81 - t112 * t71 + t201 * t70, t103 * t140 + t112 * t131 + t132 * t80 + t144 * t70, -t101 * t140 + t131 * t201 - t132 * t81 - t144 * t71, t131 * t144 + t132 * t140, t204 * t132 + t203 * t131 - t225 * t144 + t14 * t140 + t36 * t101 + t55 * t71 - t28 * t201 + t39 * t81 + (-t132 * t202 - t144 * t15) * qJD(5), t36 * t103 + t28 * t112 - t202 * t131 - t194 * t132 - t15 * t140 - t195 * t144 + t39 * t80 + t55 * t70, t22 * t86 + t26 * t74, -t22 * t85 - t23 * t86 - t25 * t74 - t26 * t72, t100 * t26 - t201 * t22 + t71 * t86 + t74 * t81, -t100 * t25 + t201 * t23 - t71 * t85 - t72 * t81, t100 * t81 - t201 * t71 (-qJD(6) * t206 + t13 * t177 - t174 * t7) * t100 + t205 * t71 - t2 * t201 - t208 * t81 + t8 * t72 + t18 * t23 + t4 * t85 + t11 * t25 -(qJD(6) * t205 + t13 * t174 + t177 * t7) * t100 - t206 * t71 + t1 * t201 - t6 * t81 + t8 * t74 + t18 * t22 + t4 * t86 + t11 * t26; 0, 0, 0, -t176 * t234, t251 * t258, t242 * t232, -t242 * t233, 0, t180 * t237 + (qJD(2) * t190 - t186) * qJD(1), pkin(1) * t234 + (-pkin(8) * t233 + t158) * t160 - t189 (t69 - t78) * t141 + (-t68 + t79) * t137 + (-t130 * t170 - t131 * t172) * pkin(2), t68 * t78 - t69 * t79 + (-t145 * t233 + t170 * t48 - t172 * t47) * pkin(2), -t130 * t163 + t131 * t164 + (-t59 - t78) * t141 + (t58 - t252) * t137, t137 * t84 - t160 * t78 + t198, -t137 * t75 + t141 * t84 + t252 * t160 + t236, t163 * t236 + t164 * t47 - t252 * t59 - t58 * t78 - t75 * t84, -t175 * t219 + t178 * t70 (-t71 - t219) * t178 + (t241 * t101 - t70) * t175, t193 + t263, t187 - t264, t132 * t137, t14 * t137 + t163 * t71 + t253 * t101 + (t28 + (-qJD(5) * t162 + t57) * t132) * t175 + (-t52 * t132 - t283) * t178, -t15 * t137 + t163 * t70 + t28 * t178 + (-t162 * t245 + t273) * t132 + t253 * t103 + t283 * t175, -t185 * t74 + t22 * t254, t72 * t89 + t74 * t88 + (t174 * t74 + t177 * t72) * t247 + (-t270 - t177 * t23 + (t174 * t72 - t177 * t74) * qJD(6)) * t178, t240 + t282, -t268 + (t197 - t287) * t178 + t226, t100 * t220 + t175 * t71, t146 * t267 - t11 * t88 - t16 * t72 + (t221 * t174 + t266 * t177) * t100 + (-t11 * t174 * qJD(5) + t2 + (qJD(5) * t72 + t197) * t162) * t175 + (t11 * t243 - t208 * t141 - t162 * t23 + t4 * t174 + (-t162 * t265 - t208) * qJD(5)) * t178, -t146 * t269 - t11 * t89 - t16 * t74 + (-t266 * t174 + t177 * t221) * t100 + (-t11 * t246 - t1 + (qJD(5) * t74 + t196) * t162) * t175 + (-t11 * t244 - t6 * t141 - t162 * t22 + t4 * t177 + (-t162 * t218 - t6) * qJD(5)) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t286, t137 * t69 + t141 * t68 + t151, t286 -(qJD(2) + t160) * t141, -t131 + t261, t277 - t137 * t59 + (-qJD(4) - t58) * t141 + t223, 0, 0, 0, 0, 0, t187 + t264, -t193 + t263, 0, 0, 0, 0, 0, t178 * t181 + t226 + t268, t240 - t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131 + t261, -t141 * t137, -t160 ^ 2 - t133, t160 * t59 + t198, 0, 0, 0, 0, 0, -t101 * t160 + t193, -t103 * t160 + t187, 0, 0, 0, 0, 0, -t178 * t23 + (-t177 * t160 - t174 * t220) * t100 + t181 * t175, -t178 * t22 + (t174 * t160 - t177 * t220) * t100 + (t241 * t74 + t196) * t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t101, -t101 ^ 2 + t103 ^ 2, t101 * t132 + t70, t103 * t285 + t117, t131, -t103 * t39 + t15 * t285 - t225, t101 * t39 + t132 * t14 - t195, t218 * t74 + t270 (t22 - t272) * t177 + (-t23 - t271) * t174, t100 * t218 - t103 * t74 + t269, -t100 ^ 2 * t174 + t103 * t72 + t267, -t100 * t103, -pkin(5) * t23 + t103 * t208 - t15 * t72 + t184 * t174 - t177 * t284, -pkin(5) * t22 + t6 * t103 - t15 * t74 + t174 * t284 + t184 * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72, -t72 ^ 2 + t74 ^ 2, t22 + t272, -t23 + t271, t71, t100 * t6 - t11 * t74 + t2, -t100 * t208 + t11 * t72 - t1;];
tauc_reg  = t5;
