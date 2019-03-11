% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:26:32
% EndTime: 2019-03-09 21:26:45
% DurationCPUTime: 4.20s
% Computational Cost: add. (8193->413), mult. (21834->727), div. (0->0), fcn. (20520->10), ass. (0->182)
t160 = sin(qJ(3));
t159 = sin(qJ(4));
t163 = cos(qJ(3));
t219 = qJD(3) * t163;
t202 = t159 * t219;
t162 = cos(qJ(4));
t216 = qJD(4) * t162;
t247 = t160 * t216 + t202;
t246 = -0.4e1 * t160;
t180 = -pkin(3) * t163 - pkin(10) * t160;
t130 = -pkin(2) + t180;
t229 = t162 * t163;
t141 = pkin(9) * t229;
t100 = t159 * t130 + t141;
t156 = sin(pkin(11));
t235 = cos(pkin(11));
t189 = t235 * t162;
t217 = qJD(4) * t159;
t112 = qJD(4) * t189 - t156 * t217;
t154 = t162 ^ 2;
t226 = t159 ^ 2 - t154;
t187 = t226 * qJD(4);
t157 = sin(pkin(6));
t245 = 0.2e1 * t157;
t244 = 2 * qJD(6);
t158 = cos(pkin(6));
t164 = cos(qJ(2));
t232 = t157 * t164;
t207 = pkin(8) * t232;
t161 = sin(qJ(2));
t243 = pkin(1) * t161;
t109 = (t158 * t243 + t207) * qJD(2);
t233 = t157 * t161;
t114 = t158 * t160 + t163 * t233;
t223 = qJD(2) * t164;
t195 = t157 * t223;
t82 = t114 * qJD(3) + t160 * t195;
t113 = -t158 * t163 + t160 * t233;
t222 = qJD(3) * t113;
t83 = t163 * t195 - t222;
t165 = t82 * pkin(3) - t83 * pkin(10) + t109;
t224 = qJD(2) * t161;
t196 = t157 * t224;
t103 = t207 + (pkin(9) + t243) * t158;
t104 = (-pkin(2) * t164 - pkin(9) * t161 - pkin(1)) * t157;
t107 = (pkin(2) * t161 - pkin(9) * t164) * t157 * qJD(2);
t108 = -pkin(1) * t158 * t223 + pkin(8) * t196;
t221 = qJD(3) * t160;
t37 = t103 * t221 - t104 * t219 - t160 * t107 + t163 * t108;
t35 = pkin(10) * t196 - t37;
t102 = pkin(8) * t233 + (-pkin(1) * t164 - pkin(2)) * t158;
t65 = t113 * pkin(3) - t114 * pkin(10) + t102;
t238 = t163 * t103 + t160 * t104;
t67 = -pkin(10) * t232 + t238;
t15 = -t159 * t165 - t162 * t35 - t65 * t216 + t67 * t217;
t205 = t159 * t232;
t47 = -qJD(4) * t205 + t114 * t216 + t159 * t83 - t162 * t196;
t84 = t114 * t159 + t162 * t232;
t10 = -qJ(5) * t47 - qJD(5) * t84 - t15;
t32 = t159 * t65 + t162 * t67;
t16 = -t32 * qJD(4) - t159 * t35 + t162 * t165;
t48 = -t84 * qJD(4) + t159 * t196 + t83 * t162;
t85 = t114 * t162 - t205;
t8 = t82 * pkin(4) - t48 * qJ(5) - t85 * qJD(5) + t16;
t4 = t235 * t10 + t156 * t8;
t242 = pkin(9) * t157;
t241 = pkin(9) * t159;
t145 = -t235 * pkin(4) - pkin(5);
t240 = pkin(5) - t145;
t239 = -qJ(5) - pkin(10);
t31 = -t159 * t67 + t162 * t65;
t22 = pkin(4) * t113 - qJ(5) * t85 + t31;
t29 = -qJ(5) * t84 + t32;
t14 = t156 * t22 + t235 * t29;
t212 = t162 * qJD(5);
t179 = pkin(3) * t160 - pkin(10) * t163;
t173 = t179 * qJD(3);
t203 = t159 * t221;
t227 = pkin(9) * t203 + t162 * t173;
t45 = -t160 * t212 + (pkin(4) * t160 - qJ(5) * t229) * qJD(3) + (-t141 + (qJ(5) * t160 - t130) * t159) * qJD(4) + t227;
t228 = -t130 * t216 - t159 * t173;
t230 = t160 * t162;
t51 = (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t230 + (-qJD(5) * t160 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t163) * t159 - t228;
t26 = t156 * t45 + t235 * t51;
t121 = t162 * t130;
t79 = -qJ(5) * t230 + t121 + (-pkin(4) - t241) * t163;
t231 = t159 * t160;
t86 = -qJ(5) * t231 + t100;
t55 = t156 * t79 + t235 * t86;
t237 = t48 * t159;
t236 = t48 * t162;
t234 = t156 * t159;
t123 = pkin(4) * t231 + t160 * pkin(9);
t153 = t160 ^ 2;
t225 = -t163 ^ 2 + t153;
t220 = qJD(3) * t162;
t218 = qJD(3) * t164;
t215 = qJD(4) * t163;
t214 = qJD(6) * t113;
t213 = qJD(6) * t163;
t211 = t82 * qJ(6) + t4;
t210 = -0.2e1 * pkin(2) * qJD(3);
t209 = -0.2e1 * pkin(3) * qJD(4);
t208 = t163 * t241;
t206 = qJ(6) * t221 + t26;
t148 = pkin(4) * t217;
t149 = pkin(9) * t219;
t95 = t247 * pkin(4) + t149;
t147 = -pkin(4) * t162 - pkin(3);
t151 = t157 ^ 2;
t204 = t151 * t223;
t201 = t162 * t219;
t200 = t159 * t215;
t198 = t162 * t215;
t194 = t159 * t216;
t193 = t160 * t219;
t3 = -t156 * t10 + t235 * t8;
t191 = qJD(4) * t239;
t110 = t159 * t191 + t212;
t167 = -t159 * qJD(5) + t162 * t191;
t74 = t156 * t110 - t235 * t167;
t75 = t235 * t110 + t156 * t167;
t131 = t239 * t162;
t190 = t235 * t159;
t87 = -t156 * t131 - t239 * t190;
t88 = -t235 * t131 + t239 * t234;
t192 = t74 * t87 + t88 * t75;
t24 = -t156 * t51 + t235 * t45;
t188 = -t160 * t103 + t104 * t163;
t186 = t225 * qJD(3);
t185 = 0.2e1 * t193;
t184 = t161 * t204;
t183 = t159 * t201;
t66 = pkin(3) * t232 - t188;
t182 = t235 * t219;
t178 = -t159 * t85 - t162 * t84;
t38 = -t103 * t219 - t104 * t221 + t107 * t163 + t160 * t108;
t36 = -pkin(3) * t196 - t38;
t177 = t36 * t159 + t66 * t216;
t176 = -t36 * t162 + t66 * t217;
t175 = t113 * t216 + t159 * t82;
t174 = t113 * t217 - t162 * t82;
t13 = -t156 * t29 + t235 * t22;
t54 = -t156 * t86 + t235 * t79;
t25 = t156 * t48 + t235 * t47;
t27 = -t156 * t47 + t235 * t48;
t56 = t156 * t85 + t235 * t84;
t57 = -t156 * t84 + t235 * t85;
t172 = -t88 * t25 + t27 * t87 - t75 * t56 + t57 * t74;
t119 = t156 * t162 + t190;
t105 = t119 * t160;
t106 = -t156 * t231 + t160 * t189;
t72 = -t112 * t160 - t156 * t201 - t159 * t182;
t111 = t119 * qJD(4);
t73 = t160 * t111 + t156 * t202 - t162 * t182;
t171 = -t75 * t105 + t106 * t74 + t88 * t72 - t73 * t87;
t39 = pkin(4) * t84 + t66;
t170 = t160 * t218 + t163 * t224;
t169 = t160 * t224 - t163 * t218;
t168 = t160 * t220 + t200;
t118 = -t189 + t234;
t166 = -0.2e1 * t88 * t111 + 0.2e1 * t112 * t87 - 0.2e1 * t75 * t118 + 0.2e1 * t119 * t74;
t19 = pkin(4) * t47 + t36;
t142 = pkin(4) * t156 + qJ(6);
t99 = t121 - t208;
t76 = pkin(5) * t118 - qJ(6) * t119 + t147;
t70 = -t100 * qJD(4) + t227;
t69 = pkin(9) * t168 + t228;
t68 = pkin(5) * t105 - qJ(6) * t106 + t123;
t62 = pkin(5) * t111 - qJ(6) * t112 - qJD(6) * t119 + t148;
t50 = t163 * pkin(5) - t54;
t49 = -qJ(6) * t163 + t55;
t30 = -pkin(5) * t72 + qJ(6) * t73 - qJD(6) * t106 + t95;
t23 = -pkin(5) * t221 - t24;
t21 = t206 - t213;
t17 = pkin(5) * t56 - qJ(6) * t57 + t39;
t12 = -t113 * pkin(5) - t13;
t11 = qJ(6) * t113 + t14;
t5 = pkin(5) * t25 - qJ(6) * t27 - qJD(6) * t57 + t19;
t2 = -pkin(5) * t82 - t3;
t1 = t211 + t214;
t6 = [0, 0, 0, 0.2e1 * t184, 0.2e1 * (-t161 ^ 2 + t164 ^ 2) * t151 * qJD(2), 0.2e1 * t158 * t195, -0.2e1 * t158 * t196, 0, -0.2e1 * pkin(1) * t151 * t224 - 0.2e1 * t109 * t158, -0.2e1 * pkin(1) * t204 + 0.2e1 * t108 * t158, 0.2e1 * t114 * t83, -0.2e1 * t113 * t83 - 0.2e1 * t114 * t82 (t114 * t224 - t164 * t83) * t245 (-t113 * t224 + t164 * t82) * t245, -0.2e1 * t184, 0.2e1 * t102 * t82 + 0.2e1 * t109 * t113 + 0.2e1 * (-t38 * t164 + t188 * t224) * t157, 0.2e1 * t102 * t83 + 0.2e1 * t109 * t114 + 0.2e1 * (-t37 * t164 - t238 * t224) * t157, 0.2e1 * t85 * t48, -0.2e1 * t47 * t85 - 0.2e1 * t48 * t84, 0.2e1 * t113 * t48 + 0.2e1 * t82 * t85, -0.2e1 * t113 * t47 - 0.2e1 * t82 * t84, 0.2e1 * t113 * t82, 0.2e1 * t113 * t16 + 0.2e1 * t31 * t82 + 0.2e1 * t36 * t84 + 0.2e1 * t47 * t66, 0.2e1 * t113 * t15 - 0.2e1 * t32 * t82 + 0.2e1 * t36 * t85 + 0.2e1 * t48 * t66, -0.2e1 * t13 * t27 - 0.2e1 * t14 * t25 - 0.2e1 * t3 * t57 - 0.2e1 * t4 * t56, 0.2e1 * t13 * t3 + 0.2e1 * t14 * t4 + 0.2e1 * t19 * t39, -0.2e1 * t113 * t2 - 0.2e1 * t12 * t82 + 0.2e1 * t17 * t25 + 0.2e1 * t5 * t56, -0.2e1 * t1 * t56 - 0.2e1 * t11 * t25 + 0.2e1 * t12 * t27 + 0.2e1 * t2 * t57, 0.2e1 * t1 * t113 + 0.2e1 * t11 * t82 - 0.2e1 * t17 * t27 - 0.2e1 * t5 * t57, 0.2e1 * t1 * t11 + 0.2e1 * t12 * t2 + 0.2e1 * t17 * t5; 0, 0, 0, 0, 0, t195, -t196, 0, -t109, t108, t114 * t219 + t160 * t83, -t160 * t82 + t83 * t163 + (-t113 * t163 - t114 * t160) * qJD(3), t169 * t157, t170 * t157, 0, -pkin(2) * t82 + t102 * t221 - t109 * t163 - t169 * t242, -pkin(2) * t83 + t102 * t219 + t109 * t160 - t170 * t242, t85 * t201 + (-t217 * t85 + t236) * t160, t178 * t219 + (-t237 - t162 * t47 + (t159 * t84 - t162 * t85) * qJD(4)) * t160 (t113 * t220 - t48) * t163 + (qJD(3) * t85 - t174) * t160 (-t159 * t222 + t47) * t163 + (-qJD(3) * t84 - t175) * t160, t113 * t221 - t163 * t82, t70 * t113 + t99 * t82 + (-t16 + (pkin(9) * t84 + t159 * t66) * qJD(3)) * t163 + (pkin(9) * t47 + qJD(3) * t31 + t177) * t160, -t100 * t82 + t69 * t113 + (-t15 + (pkin(9) * t85 + t162 * t66) * qJD(3)) * t163 + (pkin(9) * t48 - qJD(3) * t32 - t176) * t160, -t105 * t4 - t106 * t3 + t13 * t73 + t14 * t72 - t24 * t57 - t25 * t55 - t26 * t56 - t27 * t54, t123 * t19 + t13 * t24 + t14 * t26 + t3 * t54 + t39 * t95 + t4 * t55, t105 * t5 - t113 * t23 - t12 * t221 + t163 * t2 - t17 * t72 + t25 * t68 + t30 * t56 - t50 * t82, -t1 * t105 + t106 * t2 + t11 * t72 - t12 * t73 - t21 * t56 + t23 * t57 - t25 * t49 + t27 * t50, -t1 * t163 - t106 * t5 + t11 * t221 + t113 * t21 + t17 * t73 - t27 * t68 - t30 * t57 + t49 * t82, t1 * t49 + t11 * t21 + t12 * t23 + t17 * t30 + t2 * t50 + t5 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, -0.2e1 * t186, 0, 0, 0, t160 * t210, t163 * t210, -0.2e1 * t153 * t194 + 0.2e1 * t154 * t193, 0.2e1 * t153 * t187 + t183 * t246, 0.2e1 * t160 * t200 + 0.2e1 * t225 * t220, -0.2e1 * t159 * t186 + 0.2e1 * t160 * t198, -0.2e1 * t193, 0.2e1 * t99 * t221 - 0.2e1 * t70 * t163 + 0.2e1 * (t153 * t216 + t159 * t185) * pkin(9), -0.2e1 * t100 * t221 - 0.2e1 * t69 * t163 + 0.2e1 * (-t153 * t217 + t162 * t185) * pkin(9), -0.2e1 * t105 * t26 - 0.2e1 * t106 * t24 + 0.2e1 * t54 * t73 + 0.2e1 * t55 * t72, 0.2e1 * t123 * t95 + 0.2e1 * t24 * t54 + 0.2e1 * t26 * t55, 0.2e1 * t105 * t30 + 0.2e1 * t163 * t23 - 0.2e1 * t221 * t50 - 0.2e1 * t68 * t72, -0.2e1 * t105 * t21 + 0.2e1 * t106 * t23 + 0.2e1 * t49 * t72 - 0.2e1 * t50 * t73, -0.2e1 * t106 * t30 - 0.2e1 * t163 * t21 + 0.2e1 * t221 * t49 + 0.2e1 * t68 * t73, 0.2e1 * t21 * t49 + 0.2e1 * t23 * t50 + 0.2e1 * t30 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t82, t196, t38, t37, t216 * t85 + t237, qJD(4) * t178 - t159 * t47 + t236, t175, -t174, 0, -pkin(3) * t47 - pkin(10) * t175 + t176, -pkin(3) * t48 + pkin(10) * t174 + t177, -t111 * t14 - t112 * t13 - t118 * t4 - t119 * t3 + t172, -t13 * t74 + t14 * t75 + t147 * t19 + t148 * t39 - t3 * t87 + t4 * t88, t111 * t17 - t113 * t74 + t118 * t5 + t25 * t76 + t56 * t62 - t82 * t87, -t1 * t118 - t11 * t111 + t112 * t12 + t119 * t2 + t172, -t112 * t17 + t113 * t75 - t119 * t5 - t27 * t76 - t57 * t62 + t82 * t88, t1 * t88 + t11 * t75 + t12 * t74 + t17 * t62 + t2 * t87 + t5 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, -t221, 0, -t149, pkin(9) * t221, -t160 * t187 + t183, t194 * t246 - t226 * t219, -t198 + t203, t168, 0 (pkin(10) * t229 + (-pkin(3) * t162 + t241) * t160) * qJD(4) + (t159 * t180 - t141) * qJD(3) (pkin(9) * t230 + t159 * t179) * qJD(4) + (t162 * t180 + t208) * qJD(3), -t111 * t55 - t112 * t54 - t118 * t26 - t119 * t24 + t171, t123 * t148 + t147 * t95 - t24 * t87 + t26 * t88 - t54 * t74 + t55 * t75, t105 * t62 + t111 * t68 + t118 * t30 + t163 * t74 - t221 * t87 - t72 * t76, -t111 * t49 + t112 * t50 - t118 * t21 + t119 * t23 + t171, -t106 * t62 - t112 * t68 - t119 * t30 - t163 * t75 + t221 * t88 + t73 * t76, t21 * t88 + t23 * t87 + t30 * t76 + t49 * t75 + t50 * t74 + t62 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t194, -0.2e1 * t187, 0, 0, 0, t159 * t209, t162 * t209, t166, 0.2e1 * t147 * t148 + 0.2e1 * t192, 0.2e1 * t111 * t76 + 0.2e1 * t118 * t62, t166, -0.2e1 * t112 * t76 - 0.2e1 * t119 * t62, 0.2e1 * t62 * t76 + 0.2e1 * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t47, t82, t16, t15 (-t156 * t25 - t235 * t27) * pkin(4) (t156 * t4 + t235 * t3) * pkin(4), t240 * t82 + t3, -qJD(6) * t56 - t142 * t25 + t145 * t27, t142 * t82 + t211 + 0.2e1 * t214, qJD(6) * t11 + t1 * t142 + t145 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t217 + t201, -t247, t221, t70, t69 (t156 * t72 + t235 * t73) * pkin(4) (t156 * t26 + t235 * t24) * pkin(4), t221 * t240 + t24, -qJD(6) * t105 + t142 * t72 - t145 * t73, t142 * t221 + t206 - 0.2e1 * t213, qJD(6) * t49 + t142 * t21 + t145 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, -t217, 0, -pkin(10) * t216, pkin(10) * t217 (-t111 * t156 - t235 * t112) * pkin(4) (t156 * t75 - t235 * t74) * pkin(4), -t74, -qJD(6) * t118 - t111 * t142 + t112 * t145, t75, qJD(6) * t88 + t142 * t75 + t145 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t142 * t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t25, 0, -t27, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t72, 0, t73, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t111, 0, -t112, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t27, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, -t73, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;
