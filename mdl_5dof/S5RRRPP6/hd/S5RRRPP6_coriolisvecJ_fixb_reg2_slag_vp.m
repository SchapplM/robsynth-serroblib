% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:00
% EndTime: 2019-12-31 21:02:10
% DurationCPUTime: 3.63s
% Computational Cost: add. (4278->407), mult. (11048->540), div. (0->0), fcn. (7333->6), ass. (0->204)
t166 = cos(qJ(2));
t220 = t166 * qJD(1);
t145 = -qJD(3) + t220;
t163 = sin(qJ(3));
t164 = sin(qJ(2));
t228 = qJD(1) * t164;
t205 = t163 * t228;
t165 = cos(qJ(3));
t222 = t165 * qJD(2);
t119 = t205 - t222;
t204 = t165 * t228;
t223 = t163 * qJD(2);
t121 = t204 + t223;
t162 = sin(pkin(8));
t250 = cos(pkin(8));
t70 = t250 * t119 + t162 * t121;
t257 = t70 * t145;
t219 = qJD(1) * qJD(2);
t203 = t166 * t219;
t218 = qJD(2) * qJD(3);
t84 = qJD(3) * t205 + (-t203 - t218) * t165;
t212 = t166 * t223;
t224 = qJD(3) * t165;
t172 = t164 * t224 + t212;
t85 = t172 * qJD(1) + t163 * t218;
t50 = -t162 * t85 - t250 * t84;
t27 = t50 - t257;
t173 = -t162 * t119 + t250 * t121;
t279 = t70 * t173;
t200 = t250 * t163;
t115 = t162 * t165 + t200;
t104 = t115 * qJD(3);
t252 = -t115 * t220 + t104;
t199 = t250 * t165;
t225 = qJD(3) * t163;
t105 = qJD(3) * t199 - t162 * t225;
t213 = t163 * t220;
t251 = t162 * t213 - t199 * t220 + t105;
t268 = t173 ^ 2;
t278 = -0.2e1 * t219;
t266 = -qJ(4) - pkin(7);
t201 = qJD(3) * t266;
t221 = t165 * qJD(4);
t101 = t163 * t201 + t221;
t171 = -t163 * qJD(4) + t165 * t201;
t236 = t165 * t166;
t175 = pkin(3) * t164 - qJ(4) * t236;
t185 = pkin(2) * t164 - pkin(7) * t166;
t122 = t185 * qJD(1);
t86 = pkin(6) * t205 + t165 * t122;
t64 = t175 * qJD(1) + t86;
t106 = t163 * t122;
t237 = t164 * t165;
t238 = t163 * t166;
t74 = t106 + (-pkin(6) * t237 - qJ(4) * t238) * qJD(1);
t263 = (-t171 + t64) * t250 + (t101 - t74) * t162;
t155 = pkin(6) * t228;
t260 = qJD(2) * pkin(2);
t133 = t155 - t260;
t83 = t119 * pkin(3) + qJD(4) + t133;
t28 = t70 * pkin(4) - qJ(5) * t173 + t83;
t277 = t28 * t173;
t130 = -t166 * pkin(2) - t164 * pkin(7) - pkin(1);
t110 = t130 * qJD(1);
t123 = t185 * qJD(2);
t111 = qJD(1) * t123;
t156 = pkin(6) * t220;
t134 = qJD(2) * pkin(7) + t156;
t152 = t164 * t219;
t192 = pkin(6) * t152;
t44 = t110 * t224 + t163 * t111 - t134 * t225 - t165 * t192;
t77 = t165 * t110 - t163 * t134;
t276 = t77 * t145 + t44;
t78 = t163 * t110 + t165 * t134;
t45 = -qJD(3) * t78 + t165 * t111 + t163 * t192;
t275 = t78 * t145 - t45;
t188 = -t156 + (-t213 + t225) * pkin(3);
t148 = pkin(6) * t236;
t93 = t163 * t130 + t148;
t227 = qJD(2) * t164;
t49 = -t162 * t84 + t250 * t85;
t197 = qJD(2) * t250;
t187 = t166 * t197;
t211 = t166 * t222;
t59 = -t105 * t164 - t162 * t211 - t163 * t187;
t96 = t115 * t164;
t274 = t59 * t145 - t166 * t49 + (qJD(1) * t96 + t70) * t227;
t240 = t162 * t163;
t114 = -t199 + t240;
t273 = t145 * t252 - (qJD(2) * t114 - t70) * t228;
t272 = t50 * t114 + t115 * t49 + t173 * t252 + t251 * t70;
t159 = t164 * pkin(6);
t58 = -t119 * qJ(4) + t78;
t53 = t250 * t58;
t57 = -t121 * qJ(4) + t77;
t25 = t162 * t57 + t53;
t271 = t25 * t173;
t269 = t70 ^ 2;
t19 = pkin(3) * t152 + t84 * qJ(4) - t121 * qJD(4) + t45;
t24 = -t85 * qJ(4) - t119 * qJD(4) + t44;
t3 = -t162 * t24 + t250 * t19;
t4 = t162 * t19 + t250 * t24;
t35 = t162 * t64 + t250 * t74;
t30 = qJ(5) * t228 + t35;
t62 = t250 * t101 + t162 * t171;
t265 = t30 - t62;
t264 = pkin(4) * t228 + t263;
t262 = -t35 + t62;
t230 = t165 * t123 + t223 * t159;
t36 = -t164 * t221 + t175 * qJD(2) + (-t148 + (qJ(4) * t164 - t130) * t163) * qJD(3) + t230;
t231 = t163 * t123 + t130 * t224;
t41 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t237 + (-qJD(4) * t164 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t166) * t163 + t231;
t11 = t162 * t36 + t250 * t41;
t261 = -t252 * pkin(4) + t251 * qJ(5) + t115 * qJD(5) - t188;
t51 = -t145 * pkin(3) + t57;
t21 = t162 * t51 + t53;
t117 = t165 * t130;
t75 = -qJ(4) * t237 + t117 + (-pkin(6) * t163 - pkin(3)) * t166;
t239 = t163 * t164;
t79 = -qJ(4) * t239 + t93;
t43 = t162 * t75 + t250 * t79;
t259 = t162 * t58;
t258 = t173 * t145;
t254 = t84 * t163;
t253 = t85 * t165;
t132 = t266 * t165;
t81 = -t162 * t132 - t266 * t200;
t249 = qJD(2) * t81;
t82 = -t250 * t132 + t266 * t240;
t248 = qJD(2) * t82;
t247 = t119 * t145;
t246 = t121 * t119;
t245 = t121 * t145;
t244 = t133 * t163;
t243 = t133 * t165;
t242 = t145 * t163;
t241 = t145 * t165;
t168 = qJD(1) ^ 2;
t235 = t166 * t168;
t167 = qJD(2) ^ 2;
t234 = t167 * t164;
t233 = t167 * t166;
t26 = t250 * t57 - t259;
t232 = qJD(5) - t26;
t124 = pkin(3) * t239 + t159;
t160 = t164 ^ 2;
t229 = -t166 ^ 2 + t160;
t226 = qJD(2) * t166;
t217 = pkin(6) * t238;
t216 = qJ(5) * t152 + t4;
t157 = pkin(6) * t226;
t214 = t164 * t235;
t88 = t172 * pkin(3) + t157;
t154 = -t165 * pkin(3) - pkin(2);
t210 = t164 * t225;
t208 = t145 * t224;
t207 = t145 * t228;
t68 = t85 * pkin(3) + pkin(6) * t203;
t196 = pkin(1) * t278;
t194 = t119 + t222;
t193 = -t121 + t223;
t191 = -t82 * t49 + t81 * t50 - t62 * t70;
t190 = -t25 * t145 + t3;
t189 = t166 * t152;
t184 = t49 * t96 - t70 * t59;
t183 = -t268 - t269;
t182 = t268 - t269;
t179 = -t163 * t78 - t165 * t77;
t178 = qJD(1) * t160 - t145 * t166;
t177 = t49 - t258;
t176 = t49 + t258;
t10 = -t162 * t41 + t250 * t36;
t20 = t250 * t51 - t259;
t42 = -t162 * t79 + t250 * t75;
t174 = t49 * t114 + t252 * t70;
t2 = -pkin(4) * t152 - t3;
t170 = -t50 - t257;
t60 = t164 * t104 + t162 * t212 - t165 * t187;
t97 = -t162 * t239 + t164 * t199;
t169 = -t173 * t59 + t97 * t49 + t50 * t96 - t60 * t70;
t5 = t49 * pkin(4) - t50 * qJ(5) - qJD(5) * t173 + t68;
t151 = -t250 * pkin(3) - pkin(4);
t149 = t162 * pkin(3) + qJ(5);
t92 = t117 - t217;
t91 = (-t145 - t220) * t227;
t87 = -pkin(6) * t204 + t106;
t66 = t114 * pkin(4) - t115 * qJ(5) + t154;
t56 = -t93 * qJD(3) + t230;
t55 = (-t164 * t222 - t166 * t225) * pkin(6) + t231;
t52 = t96 * pkin(4) - t97 * qJ(5) + t124;
t40 = t166 * pkin(4) - t42;
t39 = -t166 * qJ(5) + t43;
t29 = t121 * pkin(3) + pkin(4) * t173 + qJ(5) * t70;
t18 = -t251 * t145 + (qJD(2) * t115 - t173) * t228;
t15 = -t145 * qJ(5) + t21;
t14 = t145 * pkin(4) + qJD(5) - t20;
t13 = -t59 * pkin(4) + t60 * qJ(5) - t97 * qJD(5) + t88;
t12 = -t173 * t60 + t50 * t97;
t9 = -pkin(4) * t227 - t10;
t8 = qJ(5) * t227 - t166 * qJD(5) + t11;
t7 = t60 * t145 - t50 * t166 + (qJD(1) * t97 + t173) * t227;
t6 = t50 * t115 + t173 * t251;
t1 = -t145 * qJD(5) + t216;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t189, t229 * t278, t233, -0.2e1 * t189, -t234, 0, -pkin(6) * t233 + t164 * t196, pkin(6) * t234 + t166 * t196, 0, 0, -t84 * t237 + (-t210 + t211) * t121, (-t119 * t165 - t121 * t163) * t226 + (t254 - t253 + (t119 * t163 - t121 * t165) * qJD(3)) * t164, t145 * t210 + t84 * t166 + (t121 * t164 + t165 * t178) * qJD(2), t119 * t172 + t239 * t85, t164 * t208 + t85 * t166 + (-t119 * t164 - t163 * t178) * qJD(2), t91, -t56 * t145 - t45 * t166 + (pkin(6) * t85 + t133 * t224) * t164 + ((pkin(6) * t119 + t244) * t166 + (t77 + (t92 + t217) * qJD(1)) * t164) * qJD(2), t55 * t145 + t44 * t166 + (-pkin(6) * t84 - t133 * t225) * t164 + ((pkin(6) * t121 + t243) * t166 + (-t78 + (-t93 + t148) * qJD(1)) * t164) * qJD(2), -t55 * t119 - t56 * t121 + t92 * t84 - t93 * t85 + t179 * t226 + (-t163 * t44 - t165 * t45 + (t163 * t77 - t165 * t78) * qJD(3)) * t164, t44 * t93 + t45 * t92 + t78 * t55 + t77 * t56 + (t133 + t155) * t157, t12, -t169, t7, t184, -t274, t91, -t10 * t145 + t124 * t49 - t3 * t166 - t83 * t59 + t68 * t96 + t88 * t70 + (qJD(1) * t42 + t20) * t227, t11 * t145 + t124 * t50 + t4 * t166 - t83 * t60 + t68 * t97 + t88 * t173 + (-qJD(1) * t43 - t21) * t227, -t10 * t173 - t11 * t70 + t20 * t60 + t21 * t59 - t3 * t97 - t4 * t96 - t42 * t50 - t43 * t49, t20 * t10 + t21 * t11 + t68 * t124 + t3 * t42 + t4 * t43 + t83 * t88, t12, t7, t169, t91, t274, t184, t13 * t70 + t9 * t145 + t2 * t166 - t28 * t59 + t52 * t49 + t5 * t96 + (-qJD(1) * t40 - t14) * t227, -t1 * t96 - t14 * t60 + t15 * t59 + t173 * t9 + t2 * t97 - t39 * t49 + t40 * t50 - t8 * t70, -t1 * t166 - t13 * t173 - t8 * t145 + t28 * t60 - t5 * t97 - t52 * t50 + (qJD(1) * t39 + t15) * t227, t1 * t39 + t28 * t13 + t14 * t9 + t15 * t8 + t2 * t40 + t5 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, t229 * t168, 0, t214, 0, 0, t168 * pkin(1) * t164, pkin(1) * t235, 0, 0, -t121 * t241 - t254, (-t84 + t247) * t165 + (-t85 + t245) * t163, -t208 + (t145 * t236 + t164 * t193) * qJD(1), -t119 * t242 - t253, t145 * t225 + (-t145 * t238 + t164 * t194) * qJD(1), t207, -pkin(2) * t85 + t86 * t145 + (pkin(7) * t241 + t244) * qJD(3) + ((-pkin(7) * t223 - t77) * t164 + (-pkin(6) * t194 - t244) * t166) * qJD(1), pkin(2) * t84 - t87 * t145 + (-pkin(7) * t242 + t243) * qJD(3) + ((-pkin(7) * t222 + t78) * t164 + (pkin(6) * t193 - t243) * t166) * qJD(1), t87 * t119 + t86 * t121 + ((qJD(3) * t121 - t85) * pkin(7) + t276) * t165 + ((qJD(3) * t119 - t84) * pkin(7) + t275) * t163, -t77 * t86 - t78 * t87 + (-t133 - t260) * t156 + (qJD(3) * t179 - t45 * t163 + t44 * t165) * pkin(7), t6, -t272, t18, t174, t273, t207, t68 * t114 + t154 * t49 + t252 * t83 + t188 * t70 + t263 * t145 + (-t20 - t249) * t228, t68 * t115 + t154 * t50 + t251 * t83 + t188 * t173 + t262 * t145 + (t21 - t248) * t228, -t4 * t114 - t3 * t115 + t173 * t263 - t251 * t20 - t252 * t21 + t35 * t70 + t191, t68 * t154 + t188 * t83 - t263 * t20 + t262 * t21 - t3 * t81 + t4 * t82, t6, t18, t272, t207, -t273, t174, t5 * t114 + t66 * t49 - t261 * t70 + t252 * t28 + t264 * t145 + (t14 - t249) * t228, -t1 * t114 + t2 * t115 + t251 * t14 - t252 * t15 + t173 * t264 + t30 * t70 + t191, -t5 * t115 - t66 * t50 + t261 * t173 - t251 * t28 + t265 * t145 + (-t15 + t248) * t228, t1 * t82 + t264 * t14 - t265 * t15 + t2 * t81 - t261 * t28 + t5 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, -t119 ^ 2 + t121 ^ 2, -t84 - t247, -t246, -t245 - t85, t152, -t133 * t121 - t275, t133 * t119 - t276, 0, 0, t279, t182, t27, -t279, -t176, t152, -t83 * t173 + (-t121 * t70 + t197 * t228) * pkin(3) + t190, -t26 * t145 + t83 * t70 + (-t121 * t173 - t152 * t162) * pkin(3) - t4, t21 * t173 - t271 + (-t162 * t49 - t250 * t50) * pkin(3) + (-t20 + t26) * t70, t20 * t25 - t21 * t26 + (-t121 * t83 + t162 * t4 + t250 * t3) * pkin(3), t279, t27, -t182, t152, t176, -t279, -t277 - t29 * t70 + (pkin(4) - t151) * t152 + t190, -t149 * t49 + t15 * t173 + t151 * t50 - t271 + (t14 - t232) * t70, t149 * t152 - t28 * t70 + t29 * t173 + (-0.2e1 * qJD(5) + t26) * t145 + t216, t1 * t149 - t14 * t25 + t15 * t232 + t2 * t151 - t28 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, -t170, t183, t173 * t20 + t21 * t70 + t68, 0, 0, 0, 0, 0, 0, t177, t183, t170, -t14 * t173 + t15 * t70 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 + t279, t27, -t145 ^ 2 - t268, t15 * t145 + t2 + t277;];
tauc_reg = t16;
