% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:40
% EndTime: 2019-03-09 17:19:51
% DurationCPUTime: 4.29s
% Computational Cost: add. (3921->366), mult. (9265->609), div. (0->0), fcn. (7605->6), ass. (0->191)
t156 = sin(qJ(3));
t158 = cos(qJ(3));
t247 = sin(qJ(5));
t196 = qJD(3) * t247;
t248 = cos(qJ(5));
t197 = qJD(3) * t248;
t256 = t156 * t197 - t158 * t196;
t157 = sin(qJ(2));
t198 = qJD(2) * t247;
t199 = qJD(2) * t248;
t159 = cos(qJ(2));
t237 = t158 * t159;
t239 = t156 * t159;
t107 = t247 * t156 + t248 * t158;
t187 = t158 * t197;
t61 = t107 * qJD(5) - t156 * t196 - t187;
t30 = t61 * t157 + t198 * t237 - t199 * t239;
t209 = t248 * t156;
t108 = -t247 * t158 + t209;
t88 = t108 * t157;
t255 = -t30 * qJ(6) + t88 * qJD(6);
t149 = t156 * qJ(4);
t254 = -t158 * pkin(3) - t149;
t152 = t156 ^ 2;
t154 = t158 ^ 2;
t232 = t152 - t154;
t240 = t159 * pkin(2);
t241 = t157 * pkin(8);
t182 = -t240 - t241;
t176 = pkin(1) - t182;
t217 = pkin(7) * t239;
t167 = t158 * t176 + t217;
t163 = t167 * qJD(3);
t116 = t232 * qJD(3);
t146 = qJD(5) * t247;
t147 = qJD(5) * t248;
t62 = -t158 * t146 + t156 * t147 - t256;
t249 = pkin(8) - pkin(9);
t253 = t249 * t157 + pkin(1);
t115 = -pkin(2) + t254;
t238 = t157 * t158;
t137 = qJ(4) * t238;
t85 = -t137 + (pkin(3) * t156 + pkin(7)) * t157;
t229 = qJD(3) * t156;
t148 = qJD(3) * t158;
t233 = qJ(4) * t148 + t156 * qJD(4);
t90 = pkin(3) * t229 - t233;
t252 = qJD(2) * (-t115 * t159 + t241) - qJD(3) * t85 - t157 * t90;
t139 = pkin(7) * t237;
t234 = -t156 * t176 + t139;
t242 = pkin(8) * t159;
t246 = pkin(2) * t157;
t181 = -t242 + t246;
t144 = t157 * qJD(2);
t203 = t156 * t144;
t230 = qJD(2) * t158;
t236 = pkin(7) * t203 + t181 * t230;
t41 = -qJD(3) * t234 + t236;
t251 = 0.2e1 * qJD(4);
t250 = -pkin(3) - pkin(4);
t245 = pkin(2) * t158;
t244 = pkin(7) * t156;
t243 = pkin(7) * t158;
t175 = pkin(3) + t244 + t245;
t162 = t253 * t158 + (pkin(4) + t175) * t159;
t50 = t247 * t162;
t71 = -qJ(4) * t159 + t234;
t56 = t156 * t157 * pkin(9) + t71;
t26 = t248 * t56 + t50;
t123 = t249 * t158;
t192 = t249 * t247;
t68 = t248 * t123 + t156 * t192;
t225 = t159 * qJD(2);
t207 = t158 * t225;
t226 = t158 * qJD(4);
t235 = qJ(4) * t207 + t157 * t226;
t153 = t157 ^ 2;
t231 = -t159 ^ 2 + t153;
t228 = qJD(3) * t157;
t227 = qJD(3) * t159;
t224 = t159 * qJD(4);
t223 = -0.2e1 * pkin(1) * qJD(2);
t222 = -0.2e1 * pkin(2) * qJD(3);
t114 = t248 * qJ(4) + t247 * t250;
t142 = t248 * t250;
t86 = qJ(4) * t146 - t248 * qJD(4) - qJD(5) * t142;
t87 = t247 * qJD(4) + t114 * qJD(5);
t89 = t107 * t157;
t221 = -t114 * t30 - t86 * t88 + t87 * t89;
t220 = t86 * t107 + t87 * t108 - t114 * t62;
t219 = t114 * t147 - t86 * t247 - t87 * t248;
t218 = pkin(8) * t239;
t216 = pkin(7) * t238;
t215 = pkin(5) * t144;
t214 = pkin(3) * t144;
t213 = pkin(8) * t229;
t212 = pkin(8) * t148;
t211 = pkin(7) * t225;
t210 = t250 * t156;
t206 = t156 * t228;
t205 = t156 * t227;
t204 = t158 * t227;
t202 = t156 * t148;
t201 = t157 * t225;
t200 = t158 * t144;
t173 = -t240 - t253;
t160 = -t224 + (t158 * t173 - t217) * qJD(3) + ((qJ(4) - t243) * t157 + (-t249 * t159 + t246) * t156) * qJD(2);
t161 = (t156 * t173 + t139) * qJD(3) + (-pkin(9) * t237 + t250 * t157) * qJD(2) - t236;
t195 = qJD(5) * t50 + t56 * t147 + t247 * t160 - t248 * t161;
t194 = t231 * qJD(2);
t193 = 0.2e1 * t201;
t104 = t158 * pkin(4) - t115;
t191 = t159 * t200;
t190 = t153 * t202;
t184 = -pkin(7) + t210;
t51 = t248 * t162;
t25 = -t247 * t56 + t51;
t183 = pkin(2) * t156 - t243;
t113 = -t247 * qJ(4) + t142;
t110 = t249 * t209;
t67 = -t247 * t123 + t110;
t73 = -t158 * (-pkin(1) - t241) + t175 * t159;
t179 = -t156 * t71 + t158 * t73;
t178 = -t156 * t234 + t158 * t167;
t174 = t156 * t181;
t76 = qJD(3) * t210 + t233;
t4 = -qJD(5) * t51 + t56 * t146 - t248 * t160 - t247 * t161;
t66 = t157 * t184 + t137;
t97 = t157 * t148 + t156 * t225;
t32 = -qJD(5) * t110 + t123 * t146 + t249 * t256;
t31 = t107 * t225 + t62 * t157;
t172 = t31 * qJ(6) + t89 * qJD(6) + t195;
t169 = t254 * qJD(3) + t226;
t168 = t114 * t144 + t86 * t159 - t4;
t34 = -t224 - t163 + (-t218 + (qJ(4) + t183) * t157) * qJD(2);
t36 = -t214 - t41;
t165 = qJD(3) * t179 + t36 * t156 + t34 * t158;
t40 = t163 + (-t174 + t216) * qJD(2);
t164 = qJD(3) * t178 - t41 * t156 - t40 * t158;
t28 = (t250 * t158 - t149) * t228 + t184 * t225 + t235;
t33 = t123 * t147 - t249 * t187 + (qJD(5) * t192 - t249 * t196) * t156;
t130 = -0.2e1 * t201;
t129 = -0.2e1 * t202;
t128 = 0.2e1 * t202;
t127 = pkin(8) * t204;
t112 = -pkin(5) + t113;
t98 = t203 - t204;
t96 = t200 + t205;
t95 = -t206 + t207;
t94 = -t159 * t147 + t157 * t198;
t93 = -t159 * t146 - t157 * t199;
t81 = 0.2e1 * t87;
t80 = 0.2e1 * t86;
t79 = t87 * t159;
t75 = 0.2e1 * t154 * t201 - 0.2e1 * t190;
t74 = 0.2e1 * t152 * t201 + 0.2e1 * t190;
t72 = -t156 * t207 + t232 * t228;
t70 = -t156 * t194 + t157 * t204;
t69 = 0.4e1 * t157 * t202 + t232 * t225;
t65 = 0.2e1 * t157 * t205 + 0.2e1 * t231 * t230;
t64 = t153 * t116 - 0.2e1 * t156 * t191;
t63 = pkin(5) * t107 + t104;
t60 = t114 * t86;
t49 = -qJ(6) * t107 + t68;
t48 = -t108 * qJ(6) + t67;
t44 = -0.2e1 * t108 * t61;
t43 = 0.2e1 * t107 * t62;
t42 = -t88 * pkin(5) + t66;
t39 = t97 * pkin(3) + qJ(4) * t206 + t211 - t235;
t38 = -t108 * t144 - t159 * t61;
t37 = t107 * t144 - t159 * t62;
t35 = t62 * pkin(5) + t76;
t20 = 0.2e1 * t89 * t31;
t19 = -0.2e1 * t88 * t30;
t18 = -0.2e1 * t89 * t144 + 0.2e1 * t159 * t31;
t17 = -0.2e1 * t88 * t144 - 0.2e1 * t159 * t30;
t16 = 0.2e1 * t107 * t61 - 0.2e1 * t108 * t62;
t15 = qJ(6) * t88 + t26;
t14 = t159 * pkin(5) - t89 * qJ(6) + t25;
t13 = t248 * t61 - t247 * t62 + (-t248 * t107 + t247 * t108) * qJD(5);
t12 = -t61 * qJ(6) + t108 * qJD(6) + t33;
t11 = t62 * qJ(6) + t107 * qJD(6) + t32;
t10 = t108 * t31 - t61 * t89;
t9 = t107 * t30 - t62 * t88;
t8 = t30 * pkin(5) + t28;
t7 = -0.2e1 * t30 * t89 + 0.2e1 * t31 * t88;
t6 = -t248 * t31 - t247 * t30 + (t247 * t89 + t248 * t88) * qJD(5);
t3 = -t107 * t31 - t108 * t30 - t61 * t88 - t62 * t89;
t2 = t4 - t255;
t1 = -t172 - t215;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -0.2e1 * t194, 0, t130, 0, 0, t157 * t223, t159 * t223, 0, 0, t75, 0.2e1 * t64, t65, t74, 0.2e1 * t70, t130, -0.2e1 * t167 * t144 - 0.2e1 * t41 * t159 + 0.2e1 * (t153 * t148 + t156 * t193) * pkin(7), -0.2e1 * t234 * t144 - 0.2e1 * t40 * t159 + 0.2e1 * (-t153 * t229 + 0.2e1 * t191) * pkin(7), 0.2e1 * t178 * t225 + 0.2e1 * (t156 * t40 - t158 * t41 + (-t156 * t167 - t158 * t234) * qJD(3)) * t157, 0.2e1 * pkin(7) ^ 2 * t201 - 0.2e1 * t167 * t41 - 0.2e1 * t234 * t40, t75, t65, -0.2e1 * t64, t130, -0.2e1 * t70, t74, 0.2e1 * (qJD(2) * t156 * t85 + t36) * t159 + 0.2e1 * (-qJD(2) * t73 + t85 * t148 + t39 * t156) * t157, 0.2e1 * t179 * t225 + 0.2e1 * (-t156 * t34 + t158 * t36 + (-t156 * t73 - t158 * t71) * qJD(3)) * t157, 0.2e1 * (-t85 * t230 - t34) * t159 + 0.2e1 * (qJD(2) * t71 - t39 * t158 + t85 * t229) * t157, 0.2e1 * t34 * t71 + 0.2e1 * t36 * t73 + 0.2e1 * t39 * t85, t20, t7, t18, t19, t17, t130, -0.2e1 * t25 * t144 - 0.2e1 * t159 * t195 - 0.2e1 * t28 * t88 + 0.2e1 * t66 * t30, 0.2e1 * t26 * t144 + 0.2e1 * t159 * t4 + 0.2e1 * t28 * t89 + 0.2e1 * t66 * t31, 0.2e1 * t195 * t89 - 0.2e1 * t25 * t31 - 0.2e1 * t26 * t30 - 0.2e1 * t4 * t88, -0.2e1 * t195 * t25 - 0.2e1 * t26 * t4 + 0.2e1 * t28 * t66, t20, t7, t18, t19, t17, t130, 0.2e1 * t1 * t159 - 0.2e1 * t14 * t144 + 0.2e1 * t42 * t30 - 0.2e1 * t8 * t88, 0.2e1 * t144 * t15 + 0.2e1 * t159 * t2 + 0.2e1 * t42 * t31 + 0.2e1 * t8 * t89, -0.2e1 * t1 * t89 - 0.2e1 * t14 * t31 - 0.2e1 * t15 * t30 - 0.2e1 * t2 * t88, 0.2e1 * t1 * t14 - 0.2e1 * t15 * t2 + 0.2e1 * t42 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, 0, -t144, 0, -t211, pkin(7) * t144, 0, 0, -t72, -t69, t98, t72, t96, 0, t127 + (t244 - t245) * t228 + (t156 * t182 - t139) * qJD(2) (t174 + t216) * qJD(3) + (t158 * t182 + t217) * qJD(2), t164, -pkin(2) * t211 + pkin(8) * t164, -t72, t98, t69, 0, -t96, t72, t127 + (t115 * t228 - t39) * t158 - t252 * t156, t165 (-t39 + (t115 * t157 + t242) * qJD(3)) * t156 + t252 * t158, pkin(8) * t165 + t39 * t115 + t85 * t90, t10, t3, t38, t9, t37, 0, t104 * t30 + t28 * t107 - t67 * t144 - t159 * t33 + t66 * t62 - t76 * t88, t104 * t31 + t28 * t108 + t68 * t144 + t159 * t32 - t66 * t61 + t76 * t89, t107 * t4 + t108 * t195 + t25 * t61 - t26 * t62 - t30 * t68 - t31 * t67 - t32 * t88 + t33 * t89, t104 * t28 - t195 * t67 - t25 * t33 - t26 * t32 - t4 * t68 + t66 * t76, t10, t3, t38, t9, t37, 0, t8 * t107 - t12 * t159 - t144 * t48 + t63 * t30 - t35 * t88 + t42 * t62, t8 * t108 + t11 * t159 + t144 * t49 + t63 * t31 + t35 * t89 - t42 * t61, -t1 * t108 + t107 * t2 - t11 * t88 + t12 * t89 + t14 * t61 - t15 * t62 - t30 * t49 - t31 * t48, t1 * t48 - t11 * t15 - t12 * t14 - t2 * t49 + t35 * t42 + t63 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -0.2e1 * t116, 0, t129, 0, 0, t156 * t222, t158 * t222, 0, 0, t128, 0, 0.2e1 * t116, 0, 0, t129, 0.2e1 * t115 * t229 - 0.2e1 * t158 * t90, 0, -0.2e1 * t115 * t148 - 0.2e1 * t156 * t90, 0.2e1 * t115 * t90, t44, t16, 0, t43, 0, 0, 0.2e1 * t104 * t62 + 0.2e1 * t107 * t76, -0.2e1 * t104 * t61 + 0.2e1 * t108 * t76, 0.2e1 * t107 * t32 + 0.2e1 * t108 * t33 + 0.2e1 * t61 * t67 - 0.2e1 * t62 * t68, 0.2e1 * t104 * t76 - 0.2e1 * t32 * t68 - 0.2e1 * t33 * t67, t44, t16, 0, t43, 0, 0, 0.2e1 * t107 * t35 + 0.2e1 * t62 * t63, 0.2e1 * t108 * t35 - 0.2e1 * t61 * t63, 0.2e1 * t107 * t11 + 0.2e1 * t108 * t12 + 0.2e1 * t48 * t61 - 0.2e1 * t49 * t62, -0.2e1 * t11 * t49 - 0.2e1 * t12 * t48 + 0.2e1 * t35 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, -t97, t144, t41, t40, 0, 0, 0, t95, 0, t144, t97, 0, t41 + 0.2e1 * t214 (-pkin(3) * t225 - qJ(4) * t228) * t158 + (-qJ(4) * t225 + (pkin(3) * qJD(3) - qJD(4)) * t157) * t156, -0.2e1 * t224 - t163 + (-t218 + (0.2e1 * qJ(4) + t183) * t157) * qJD(2), -pkin(3) * t36 + qJ(4) * t34 + qJD(4) * t71, 0, 0, -t31, 0, t30, t144, -t113 * t144 + t195 - t79, t168, -t113 * t31 + t221, -t113 * t195 - t114 * t4 - t25 * t87 - t26 * t86, 0, 0, -t31, 0, t30, t144, -t79 + (pkin(5) - t112) * t144 + t172, t168 + t255, -t112 * t31 + t221, t1 * t112 - t114 * t2 - t14 * t87 - t15 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, -t229, 0, -t212, t213, 0, 0, 0, t148, 0, 0, t229, 0, -t212, t169, -t213, t169 * pkin(8), 0, 0, t61, 0, t62, 0, t33, -t32, t113 * t61 + t220, -t113 * t33 - t114 * t32 - t67 * t87 - t68 * t86, 0, 0, t61, 0, t62, 0, t12, -t11, t112 * t61 + t220, -t11 * t114 - t112 * t12 - t48 * t87 - t49 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, qJ(4) * t251, 0, 0, 0, 0, 0, 0, t81, -t80, 0, -0.2e1 * t113 * t87 - 0.2e1 * t60, 0, 0, 0, 0, 0, 0, t81, -t80, 0, -0.2e1 * t112 * t87 - 0.2e1 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t95, 0, t36, 0, 0, 0, 0, 0, 0, t93, t94, t6, -t195 * t248 - t4 * t247 + (-t247 * t25 + t248 * t26) * qJD(5), 0, 0, 0, 0, 0, 0, t93, t94, t6, t1 * t248 - t2 * t247 + (-t247 * t14 + t248 * t15) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, t212, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t33 * t248 - t32 * t247 + (-t247 * t67 + t248 * t68) * qJD(5), 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12 * t248 - t11 * t247 + (-t247 * t48 + t248 * t49) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t147, 0, -t113 * t146 + t219, 0, 0, 0, 0, 0, 0, t146, t147, 0, -t112 * t146 + t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t30, -t144, -t195, t4, 0, 0, 0, 0, t31, 0, -t30, -t144, -t172 - 0.2e1 * t215, t2, -t31 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, 0, -t62, 0, -t33, t32, 0, 0, 0, 0, -t61, 0, -t62, 0, -t12, t11, t61 * pkin(5), -t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t86, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t86, 0, -t87 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t146, -t147, 0, 0, 0, 0, 0, 0, 0, 0, -t146, -t147, 0, -pkin(5) * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t31, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t61, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
