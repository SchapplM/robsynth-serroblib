% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:10:59
% EndTime: 2019-03-09 09:11:12
% DurationCPUTime: 4.97s
% Computational Cost: add. (3810->372), mult. (9415->636), div. (0->0), fcn. (8387->8), ass. (0->218)
t235 = pkin(2) + pkin(3);
t84 = pkin(4) + t235;
t95 = cos(qJ(5));
t244 = t95 * t84;
t92 = sin(qJ(5));
t125 = t92 * qJ(4) + t244;
t96 = cos(qJ(2));
t106 = (-t95 * qJD(4) + qJD(5) * t125) * t96;
t93 = sin(qJ(2));
t211 = qJ(3) * t93;
t164 = pkin(1) + t211;
t121 = qJD(5) * (-pkin(9) * t93 + t164);
t190 = t92 * qJD(3);
t248 = -t121 * t95 - t190 * t93 - t106;
t247 = t235 * t96;
t188 = t95 * qJD(3);
t245 = t92 * t84;
t246 = t121 * t92 - t188 * t93 - (qJD(5) * (t95 * qJ(4) - t245) + t92 * qJD(4)) * t96;
t214 = pkin(8) - qJ(4);
t89 = sin(pkin(6));
t220 = t89 * t96;
t152 = t214 * t220;
t201 = qJD(4) * t93;
t208 = cos(pkin(6));
t166 = pkin(1) * t208;
t77 = t93 * t166;
t24 = -(-t77 - t152) * qJD(2) - t89 * t201;
t91 = sin(qJ(6));
t217 = t91 * t92;
t90 = qJ(3) - pkin(9);
t243 = t90 * t217;
t200 = qJD(4) * t96;
t207 = qJD(2) * t93;
t117 = (-t214 * t207 - t200) * t89;
t158 = t208 * qJD(3);
t160 = qJD(2) * t208;
t144 = t96 * t160;
t67 = pkin(1) * t144;
t135 = t158 + t67;
t241 = t117 + t135;
t85 = t91 ^ 2;
t94 = cos(qJ(6));
t87 = t94 ^ 2;
t162 = (t85 - t87) * qJD(6);
t195 = qJD(5) * t95;
t127 = t195 * t90 + t190;
t150 = pkin(5) * t95 + pkin(10) * t92;
t126 = t150 + t84;
t120 = t94 * t126;
t219 = t90 * t95;
t183 = t91 * t219;
t34 = t120 - t183;
t182 = t94 * t219;
t35 = t126 * t91 + t182;
t136 = t34 * t91 - t35 * t94;
t233 = pkin(10) * t95;
t234 = pkin(5) * t92;
t149 = t233 - t234;
t197 = qJD(5) * t91;
t191 = qJD(6) * t95;
t172 = t91 * t191;
t196 = qJD(5) * t94;
t52 = t196 * t92 + t172;
t17 = -qJD(6) * t120 - t149 * t197 - t188 * t94 + t52 * t90;
t18 = -t91 * t188 - t35 * qJD(6) + (t149 * t94 + t243) * qJD(5);
t238 = qJD(6) * t136 + t17 * t91 - t18 * t94;
t57 = pkin(8) * t220 + t77;
t40 = qJ(3) * t208 + t57;
t116 = pkin(9) * t208 - t40;
t111 = qJD(5) * t116;
t102 = t111 * t92 + t135 * t95;
t143 = t90 * t93 + pkin(1);
t209 = qJ(4) * t89;
t178 = t96 * t209;
t13 = t95 * (-t116 - t178) + t92 * (t84 * t96 + t143) * t89;
t11 = pkin(10) * t220 + t13;
t118 = t95 * t214 + t245;
t115 = (-pkin(10) - t118) * t93;
t159 = qJD(5) * t208;
t206 = qJD(2) * t96;
t26 = -t92 * t159 + (t195 * t93 + t206 * t92) * t89;
t221 = t89 * t93;
t184 = t92 * t221;
t69 = t89 * t206;
t27 = -qJD(5) * t184 + (-t159 + t69) * t95;
t151 = t26 * pkin(5) - t27 * pkin(10);
t165 = t90 * t96;
t154 = t92 * t165;
t193 = qJD(6) * t91;
t56 = -pkin(8) * t221 + t166 * t96;
t42 = -pkin(2) * t208 - t56;
t105 = -pkin(3) * t208 + t42;
t104 = pkin(4) * t208 - t105;
t48 = t208 * t95 + t184;
t49 = -t208 * t92 + t95 * t221;
t101 = t48 * pkin(5) - t49 * pkin(10) + t104;
t179 = t93 * t209;
t100 = t101 + t179;
t99 = t94 * t100;
t1 = t11 * t193 - qJD(6) * t99 - t94 * (((t154 + t115) * qJD(2) - t248) * t89 + t102) - t91 * (t151 - t24);
t7 = -t91 * t11 + t99;
t215 = t94 * t11;
t8 = t100 * t91 + t215;
t146 = t7 * t91 - t8 * t94;
t145 = t93 * t160;
t2 = -t91 * t102 + t94 * (-pkin(1) * t145 + t151) + (-t101 * t91 - t215) * qJD(6) + (t94 * t201 + (-pkin(1) * t195 + (-qJD(6) * qJ(4) - t127) * t93 - t106) * t91 + ((-t94 * t214 - t243) * t96 - t91 * t115) * qJD(2)) * t89;
t237 = qJD(6) * t146 + t1 * t91 - t2 * t94;
t83 = t89 ^ 2;
t43 = 0.2e1 * (-t93 ^ 2 + t96 ^ 2) * t83 * qJD(2);
t236 = 0.2e1 * t89;
t98 = 0.2e1 * qJD(3);
t170 = t89 * t207;
t29 = t91 * t220 + t49 * t94;
t15 = qJD(6) * t29 + t170 * t94 + t27 * t91;
t232 = t15 * t94;
t16 = -t49 * t193 - t91 * t170 + (qJD(6) * t220 + t27) * t94;
t231 = t16 * t91;
t230 = t16 * t94;
t229 = t27 * t95;
t28 = -t94 * t220 + t49 * t91;
t228 = t28 * t91;
t227 = t28 * t94;
t226 = t29 * t91;
t225 = t29 * t94;
t45 = t57 * qJD(2);
t224 = t45 * t93;
t223 = t48 * t95;
t222 = t49 * t92;
t218 = t91 * t15;
t216 = t92 * t26;
t86 = t92 ^ 2;
t88 = t95 ^ 2;
t212 = t86 - t88;
t210 = qJ(3) * t96;
t205 = qJD(3) * t90;
t204 = qJD(3) * t93;
t203 = qJD(3) * t94;
t202 = qJD(3) * t96;
t199 = qJD(5) * t28;
t198 = qJD(5) * t29;
t81 = qJD(5) * t92;
t194 = qJD(5) * t96;
t192 = qJD(6) * t92;
t80 = qJD(6) * t94;
t186 = 0.2e1 * t48 * t26;
t185 = -0.2e1 * pkin(5) * qJD(6);
t181 = -0.2e1 * qJD(5) * t84;
t180 = pkin(10) * t80;
t177 = t83 * t206;
t176 = t48 * t197;
t175 = t48 * t196;
t174 = t87 * t195;
t173 = t94 * t195;
t171 = t94 * t191;
t168 = t91 * t80;
t167 = t92 * t195;
t163 = t45 * t208;
t161 = t212 * qJD(5);
t157 = t86 * t168;
t156 = t93 * t177;
t155 = t91 * t173;
t153 = t95 * t165;
t147 = t7 * t94 + t8 * t91;
t142 = -pkin(2) * t96 - t211;
t113 = t92 * t116;
t128 = t95 * t143;
t12 = t113 + (t125 * t96 + t128) * t89;
t141 = -t12 * t92 + t13 * t95;
t139 = -t226 + t227;
t138 = t226 + t227;
t137 = t34 * t94 + t35 * t91;
t44 = pkin(8) * t170 - t67;
t134 = t89 * t145;
t133 = t89 * t144;
t10 = -t113 + (-t128 + (-pkin(5) - t125) * t96) * t89;
t103 = t111 * t95 - t135 * t92;
t119 = t92 * t214 - t244;
t4 = ((-t153 + (pkin(5) - t119) * t93) * qJD(2) + t246) * t89 - t103;
t132 = t10 * t80 + t4 * t91;
t131 = t10 * t193 - t4 * t94;
t130 = t26 * t95 - t48 * t81;
t129 = t26 * t91 + t48 * t80;
t19 = t193 * t48 - t26 * t94;
t123 = t194 * t92 + t207 * t95;
t122 = -t194 * t95 + t207 * t92;
t51 = t192 * t91 - t173;
t53 = t195 * t91 + t80 * t92;
t114 = t67 + t117;
t112 = -qJD(6) * t147 - t1 * t94 - t2 * t91;
t110 = -t232 + t231 + (t225 + t228) * qJD(6);
t109 = -qJD(6) * t137 - t17 * t94 - t18 * t91;
t82 = qJ(3) * t98;
t79 = 0.2e1 * t158;
t73 = t85 * t195;
t72 = t86 * t205;
t68 = -0.2e1 * t167;
t65 = t87 * t167;
t64 = t85 * t167;
t61 = -0.2e1 * t156;
t60 = 0.2e1 * t156;
t59 = 0.2e1 * t133;
t58 = -0.2e1 * t134;
t54 = -t81 * t91 + t171;
t50 = -t73 - t174;
t41 = (-pkin(1) + t142) * t89;
t38 = t123 * t89;
t37 = t122 * t89;
t36 = t158 - t44;
t33 = -t162 * t92 + t155;
t32 = (-t204 + (pkin(2) * t93 - t210) * qJD(2)) * t89;
t31 = (t164 + t247) * t89;
t30 = t40 - t178;
t25 = t105 - t179;
t23 = (t204 + (-t235 * t93 + t210) * qJD(2)) * t89;
t21 = t104 + t179;
t6 = ((t119 * t93 + t153) * qJD(2) - t246) * t89 + t103;
t5 = ((t118 * t93 - t154) * qJD(2) + t248) * t89 - t102;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t43, t59, t61, t58, 0, -0.2e1 * pkin(1) * t207 * t83 - 0.2e1 * t163, -0.2e1 * pkin(1) * t177 + 0.2e1 * t208 * t44 (-t44 * t96 + t224 + (-t56 * t96 - t57 * t93) * qJD(2)) * t236, -0.2e1 * t44 * t57 - 0.2e1 * t45 * t56, t60, t59, -t43, 0, 0.2e1 * t134, t61, -0.2e1 * t163 + 0.2e1 * (t207 * t41 - t32 * t96) * t89 (t36 * t96 + t224 + (-t40 * t93 + t42 * t96) * qJD(2)) * t236, 0.2e1 * t36 * t208 + 0.2e1 * (-t206 * t41 - t32 * t93) * t89, 0.2e1 * t32 * t41 + 0.2e1 * t36 * t40 + 0.2e1 * t42 * t45, t60, -t43, -0.2e1 * t133, t61, t58, 0, -0.2e1 * t24 * t208 + 0.2e1 * (-t207 * t31 + t23 * t96) * t89, 0.2e1 * (t158 + t114) * t208 + 0.2e1 * t23 * t221 + 0.2e1 * t31 * t69 (-t24 * t93 + (t200 * t89 - t135) * t96 + (-t25 * t96 + (t30 + t152) * t93) * qJD(2)) * t236, 0.2e1 * t31 * t23 + 0.2e1 * t25 * t24 + 0.2e1 * t241 * t30, 0.2e1 * t49 * t27, -0.2e1 * t26 * t49 - 0.2e1 * t27 * t48 (-t207 * t49 + t27 * t96) * t236, t186 (t207 * t48 - t26 * t96) * t236, t61, 0.2e1 * t21 * t26 - 0.2e1 * t24 * t48 + 0.2e1 * (-t12 * t207 + t6 * t96) * t89, 0.2e1 * t21 * t27 - 0.2e1 * t24 * t49 + 0.2e1 * (t13 * t207 + t5 * t96) * t89, -0.2e1 * t12 * t27 - 0.2e1 * t13 * t26 + 0.2e1 * t48 * t5 - 0.2e1 * t49 * t6, 0.2e1 * t12 * t6 - 0.2e1 * t13 * t5 - 0.2e1 * t21 * t24, 0.2e1 * t29 * t16, -0.2e1 * t15 * t29 - 0.2e1 * t16 * t28, 0.2e1 * t16 * t48 + 0.2e1 * t26 * t29, 0.2e1 * t28 * t15, -0.2e1 * t15 * t48 - 0.2e1 * t26 * t28, t186, 0.2e1 * t10 * t15 + 0.2e1 * t2 * t48 + 0.2e1 * t26 * t7 + 0.2e1 * t28 * t4, 0.2e1 * t1 * t48 + 0.2e1 * t10 * t16 - 0.2e1 * t26 * t8 + 0.2e1 * t29 * t4, 0.2e1 * t1 * t28 - 0.2e1 * t15 * t8 - 0.2e1 * t16 * t7 - 0.2e1 * t2 * t29, -0.2e1 * t1 * t8 + 0.2e1 * t10 * t4 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, -t170, 0, -t45, t44, 0, 0, 0, t69, 0, 0, t170, 0, -t45 (qJD(2) * t142 + t202) * t89, -t44 + t79, -pkin(2) * t45 + qJ(3) * t36 + qJD(3) * t40, 0, 0, -t69, 0, -t170, 0, -t24, t79 + t114 (-t202 + (t211 + t247) * qJD(2)) * t89, t241 * qJ(3) + t30 * qJD(3) - t235 * t24, -t195 * t49 - t27 * t92, t216 - t229 + (t222 + t223) * qJD(5), t37, t130, t38, 0, -t21 * t81 - t24 * t95 + t26 * t84 + (t122 * t90 - t190 * t96) * t89, -t21 * t195 + t24 * t92 + t27 * t84 + (t123 * t90 - t188 * t96) * t89 (-qJD(3) * t48 - t26 * t90 + t5 + (t49 * t90 + t12) * qJD(5)) * t95 + (qJD(3) * t49 + t27 * t90 + t6 + (t48 * t90 + t13) * qJD(5)) * t92, -t24 * t84 + t141 * qJD(3) + (-t5 * t95 - t6 * t92 + (-t12 * t95 - t13 * t92) * qJD(5)) * t90, -t92 * t230 + t29 * t51, t138 * t195 + (t232 + t231 + (t225 - t228) * qJD(6)) * t92 (t16 - t175) * t95 + (t19 - t198) * t92, -t15 * t217 - t28 * t53 (-t15 + t176) * t95 + (t129 + t199) * t92, t130, t18 * t48 + t26 * t34 + (t2 + (-t10 * t91 + t28 * t90) * qJD(5)) * t95 + (qJD(3) * t28 - qJD(5) * t7 + t15 * t90 - t132) * t92, t17 * t48 - t26 * t35 + (t1 + (-t10 * t94 + t29 * t90) * qJD(5)) * t95 + (qJD(3) * t29 + qJD(5) * t8 + t16 * t90 + t131) * t92, t147 * t195 - t15 * t35 - t16 * t34 + t17 * t28 - t18 * t29 - t237 * t92, t4 * t90 * t92 - t1 * t35 + t10 * t127 - t17 * t8 + t18 * t7 + t2 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t82, 0, 0, 0, 0, 0, 0, 0, t98, 0, t82, 0.2e1 * t167, -0.2e1 * t161, 0, t68, 0, 0, t92 * t181, t95 * t181 (-t86 - t88) * t98, 0.2e1 * t205 * t88 + 0.2e1 * t72, 0.2e1 * t65 - 0.2e1 * t157, -0.4e1 * t155 * t92 + 0.2e1 * t162 * t86, 0.2e1 * t172 * t92 + 0.2e1 * t196 * t212, 0.2e1 * t64 + 0.2e1 * t157, -0.2e1 * t161 * t91 + 0.2e1 * t171 * t92, t68, 0.2e1 * t18 * t95 + 0.2e1 * (-qJD(3) * t91 - t80 * t90) * t86 + 0.2e1 * (-t34 - 0.2e1 * t183) * t81, 0.2e1 * t17 * t95 + 0.2e1 * (t193 * t90 - t203) * t86 + 0.2e1 * (t35 - 0.2e1 * t182) * t81, 0.2e1 * t137 * t195 - 0.2e1 * t238 * t92, 0.2e1 * t167 * t90 ^ 2 - 0.2e1 * t17 * t35 + 0.2e1 * t18 * t34 + 0.2e1 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t24, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, t24, 0, 0, 0, 0, 0, 0, t19, t129, qJD(6) * t139 + t218 + t230, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t195, 0, 0, 0, 0, 0, 0, 0, 0, t52, t54, t50, t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, t69, 0, t23, 0, 0, 0, 0, 0, 0, -t38, t37, -t216 - t229 + (t222 - t223) * qJD(5), qJD(5) * t141 - t5 * t92 + t6 * t95, 0, 0, 0, 0, 0, 0 (-t15 - t176) * t95 + (-t129 + t199) * t92 (-t16 - t175) * t95 + (t19 + t198) * t92, t110 * t92 - t139 * t195 (-qJD(5) * t146 - t4) * t95 + (qJD(5) * t10 + t112) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t136 * t95 + t212 * t90) * qJD(5) + (t109 - t188) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t64 + 0.2e1 * t65 - 0.2e1 * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t26, -t170, t6, t5, 0, 0, t29 * t80 + t231, -qJD(6) * t138 - t218 + t230, t129, t193 * t28 - t232, -t19, 0, -pkin(5) * t15 - pkin(10) * t129 + t131, -pkin(5) * t16 + pkin(10) * t19 + t132, pkin(10) * t110 + t112, -pkin(5) * t4 + pkin(10) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, 0, t81, 0, -t127, t81 * t90 - t188, 0, 0, -t33, 0.4e1 * t168 * t92 - t174 + t73, t54, t33, -t52, 0 (-t180 + (pkin(5) * t91 - t90 * t94) * qJD(5)) * t95 + (pkin(10) * t197 - t203 + (pkin(5) * t94 + t90 * t91) * qJD(6)) * t92 (qJD(5) * t150 + t192 * t90) * t94 + (qJD(6) * t149 + t127) * t91, t109, -pkin(5) * t127 + pkin(10) * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t195, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t54, -t50 (-t234 + (t85 + t87) * t233) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t168, -0.2e1 * t162, 0, -0.2e1 * t168, 0, 0, t91 * t185, t94 * t185, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t15, t26, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t53, -t81, t18, t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t80, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, -t193, 0, -t180, pkin(10) * t193, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
