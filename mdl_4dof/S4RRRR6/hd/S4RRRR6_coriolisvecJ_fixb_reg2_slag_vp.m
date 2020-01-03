% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:51
% EndTime: 2019-12-31 17:31:00
% DurationCPUTime: 3.33s
% Computational Cost: add. (3960->357), mult. (11162->539), div. (0->0), fcn. (8225->8), ass. (0->163)
t121 = cos(qJ(3));
t213 = pkin(7) * t121;
t117 = sin(qJ(4));
t120 = cos(qJ(4));
t116 = cos(pkin(4));
t165 = t116 * qJD(1);
t146 = qJD(2) + t165;
t115 = sin(pkin(4));
t119 = sin(qJ(2));
t173 = qJD(1) * t119;
t159 = t115 * t173;
t122 = cos(qJ(2));
t198 = pkin(1) * t122;
t161 = t116 * t198;
t84 = -pkin(6) * t159 + qJD(1) * t161;
t57 = -t146 * pkin(2) - t84;
t118 = sin(qJ(3));
t69 = t118 * t159 - t121 * t146;
t71 = t118 * t146 + t121 * t159;
t28 = t69 * pkin(3) - t71 * pkin(8) + t57;
t172 = qJD(1) * t122;
t107 = t115 * t172;
t135 = t107 - qJD(3);
t179 = t115 * t122;
t109 = pkin(6) * t179;
t199 = pkin(1) * t119;
t58 = qJD(2) * pkin(7) + (t109 + (pkin(7) + t199) * t116) * qJD(1);
t83 = (-pkin(2) * t122 - pkin(7) * t119 - pkin(1)) * t115;
t65 = qJD(1) * t83;
t32 = t118 * t65 + t121 * t58;
t30 = -t135 * pkin(8) + t32;
t134 = t117 * t30 - t120 * t28;
t162 = qJD(1) * qJD(2);
t154 = t122 * t162;
t140 = t115 * t154;
t210 = qJD(3) * t69;
t50 = -t121 * t140 + t210;
t130 = t118 * t140;
t183 = qJD(3) * t71;
t51 = t130 + t183;
t110 = t116 * t199;
t176 = t109 + t110;
t89 = t176 * qJD(2);
t80 = qJD(1) * t89;
t20 = t51 * pkin(3) + t50 * pkin(8) + t80;
t163 = t121 * qJD(3);
t170 = qJD(3) * t118;
t129 = t115 * (pkin(2) * t119 - pkin(7) * t122);
t86 = qJD(2) * t129;
t78 = qJD(1) * t86;
t180 = t115 * t119;
t108 = pkin(6) * t180;
t92 = -t108 + t161;
t88 = t92 * qJD(2);
t79 = qJD(1) * t88;
t127 = -t118 * t78 - t121 * t79 - t65 * t163 + t58 * t170;
t171 = qJD(2) * t119;
t158 = t115 * t171;
t141 = qJD(1) * t158;
t7 = pkin(8) * t141 - t127;
t1 = -t134 * qJD(4) + t117 * t20 + t120 * t7;
t66 = qJD(4) + t69;
t212 = t134 * t66 + t1;
t147 = qJD(3) * t135;
t211 = -pkin(7) * t147 + t80;
t31 = -t118 * t58 + t121 * t65;
t209 = t135 * t31 - t127;
t138 = pkin(3) * t118 - pkin(8) * t121;
t208 = qJD(4) * t213 + (t110 + (pkin(6) + t138) * t179) * qJD(1) - t138 * qJD(3);
t112 = t115 ^ 2;
t207 = -0.2e1 * t112 * t162;
t6 = t117 * t28 + t120 * t30;
t2 = -qJD(4) * t6 - t117 * t7 + t120 * t20;
t205 = -t6 * t66 - t2;
t203 = t135 * t69;
t202 = t135 * t71;
t49 = -t117 * t135 + t120 * t71;
t182 = qJD(4) * t49;
t22 = -t117 * t50 - t120 * t141 + t182;
t152 = t118 * t79 - t121 * t78;
t14 = -qJD(3) * t32 - t152;
t82 = t116 * pkin(7) + t176;
t190 = t118 * t83 + t121 * t82;
t24 = -qJD(3) * t190 - t118 * t88 + t121 * t86;
t123 = qJD(1) ^ 2;
t150 = t120 * t135;
t47 = t117 * t71 + t150;
t197 = t47 * t66;
t196 = t49 * t47;
t195 = t49 * t66;
t90 = -t116 * t121 + t118 * t180;
t194 = t51 * t90;
t193 = t71 * t69;
t105 = -t121 * pkin(3) - t118 * pkin(8) - pkin(2);
t168 = qJD(4) * t117;
t85 = qJD(1) * t129;
t44 = t118 * t85 + t121 * t84;
t37 = pkin(8) * t159 + t44;
t192 = t105 * t168 + (-t170 * pkin(7) - t37) * t117 + t208 * t120;
t167 = qJD(4) * t120;
t169 = qJD(3) * t120;
t191 = t118 * t169 * pkin(7) - t105 * t167 + t208 * t117 + t120 * t37;
t189 = t117 * t51;
t188 = t117 * t66;
t187 = t120 * t51;
t21 = qJD(4) * t150 - t117 * t141 + t120 * t50 + t71 * t168;
t186 = t21 * t117;
t185 = t22 * t120;
t184 = t51 * t121;
t181 = t112 * t123;
t178 = t117 * t121;
t177 = t121 * t122;
t175 = t119 ^ 2 - t122 ^ 2;
t174 = qJD(1) * t115;
t164 = t121 * qJD(2);
t160 = t117 * t179;
t157 = qJD(2) * t179;
t156 = t115 * t116 * t123;
t151 = t120 * t66;
t149 = t122 * t135;
t148 = t135 * t115;
t145 = qJD(2) + 0.2e1 * t165;
t144 = t119 * t122 * t181;
t8 = -pkin(3) * t141 - t14;
t143 = pkin(8) * qJD(4) * t66 + t8;
t139 = pkin(1) * t207;
t137 = -t117 * t6 + t120 * t134;
t81 = t108 + (-pkin(2) - t198) * t116;
t91 = t116 * t118 + t121 * t180;
t33 = t90 * pkin(3) - t91 * pkin(8) + t81;
t35 = -pkin(8) * t179 + t190;
t9 = -t117 * t35 + t120 * t33;
t10 = t117 * t33 + t120 * t35;
t41 = -t118 * t82 + t121 * t83;
t43 = -t118 * t84 + t121 * t85;
t131 = t112 * t119 * t154;
t54 = t91 * t117 + t120 * t179;
t29 = t135 * pkin(3) - t31;
t128 = -pkin(8) * t51 + t66 * t29;
t23 = t118 * t86 + t121 * t88 + t83 * t163 - t82 * t170;
t126 = t118 * t135;
t124 = pkin(1) * (-t116 * t162 + t181);
t87 = t176 * qJD(1);
t77 = t117 * t105 + t120 * t213;
t76 = -pkin(7) * t178 + t120 * t105;
t61 = (t117 * t119 + t120 * t177) * t174;
t60 = t107 * t178 - t120 * t159;
t55 = t91 * t120 - t160;
t53 = -t90 * qJD(3) + t121 * t157;
t52 = t91 * qJD(3) + t118 * t157;
t40 = t71 * pkin(3) + t69 * pkin(8);
t36 = -pkin(3) * t159 - t43;
t34 = pkin(3) * t179 - t41;
t27 = -t54 * qJD(4) + t117 * t158 + t53 * t120;
t26 = -qJD(4) * t160 + t53 * t117 - t120 * t158 + t91 * t167;
t25 = t52 * pkin(3) - t53 * pkin(8) + t89;
t16 = -pkin(3) * t158 - t24;
t15 = pkin(8) * t158 + t23;
t12 = t117 * t40 + t120 * t31;
t11 = -t117 * t31 + t120 * t40;
t4 = -t10 * qJD(4) - t117 * t15 + t120 * t25;
t3 = t9 * qJD(4) + t117 * t25 + t120 * t15;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t131, t175 * t207, t145 * t157, -0.2e1 * t131, -t145 * t158, 0, -t80 * t116 + t119 * t139 - t89 * t146, -t79 * t116 + t122 * t139 - t88 * t146, (t119 * t80 + t122 * t79 + (-t119 * t87 - t122 * t84) * qJD(2) + (t119 * t89 + t122 * t88 + (-t119 * t176 - t122 * t92) * qJD(2)) * qJD(1)) * t115, t176 * t79 - t80 * t92 - t84 * t89 + t87 * t88, -t50 * t91 + t71 * t53, t50 * t90 - t91 * t51 - t71 * t52 - t53 * t69, -t53 * t135 + (t50 * t122 + (qJD(1) * t91 + t71) * t171) * t115, t69 * t52 + t194, t52 * t135 + (t51 * t122 + (-qJD(1) * t90 - t69) * t171) * t115, (-t112 * t172 - t148) * t171, -t24 * t135 + t89 * t69 + t81 * t51 + t80 * t90 + t57 * t52 + (-t14 * t122 + (qJD(1) * t41 + t31) * t171) * t115, t23 * t135 + t89 * t71 - t81 * t50 + t80 * t91 + t57 * t53 + (-t127 * t122 + (-qJD(1) * t190 - t32) * t171) * t115, t127 * t90 - t14 * t91 - t190 * t51 - t23 * t69 - t24 * t71 - t31 * t53 - t32 * t52 + t41 * t50, -t127 * t190 + t14 * t41 + t32 * t23 + t31 * t24 + t57 * t89 + t80 * t81, -t21 * t55 + t49 * t27, t21 * t54 - t55 * t22 - t49 * t26 - t27 * t47, -t21 * t90 + t27 * t66 + t49 * t52 + t55 * t51, t22 * t54 + t47 * t26, -t22 * t90 - t26 * t66 - t47 * t52 - t54 * t51, t66 * t52 + t194, -t134 * t52 + t16 * t47 + t2 * t90 + t34 * t22 + t29 * t26 + t4 * t66 + t9 * t51 + t8 * t54, -t1 * t90 - t10 * t51 + t16 * t49 - t34 * t21 + t29 * t27 - t3 * t66 - t6 * t52 + t8 * t55, -t1 * t54 - t10 * t22 + t134 * t27 - t2 * t55 + t9 * t21 - t6 * t26 - t3 * t47 - t4 * t49, t1 * t10 - t134 * t4 + t29 * t16 + t2 * t9 + t6 * t3 + t8 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t175 * t181, -t122 * t156, t144, t119 * t156, 0, -pkin(6) * t140 + t119 * t124 + t87 * t146, pkin(6) * t141 + t122 * t124 + t84 * t146, 0, 0, -t50 * t118 - t121 * t202, (-t50 + t203) * t121 + (-t51 + t202) * t118, -t121 * t147 + (t121 * t149 + (t118 * qJD(2) - t71) * t119) * t174, -t126 * t69 - t184, t118 * t147 + (-t118 * t149 + (t69 + t164) * t119) * t174, t148 * t173, -pkin(2) * t51 + t57 * t170 + t43 * t135 - t87 * t69 - t211 * t121 + (-t119 * t31 + (-pkin(7) * t171 - t122 * t57) * t118) * t174, pkin(2) * t50 + t57 * t163 - t44 * t135 - t87 * t71 + t211 * t118 + (-t57 * t177 + (-pkin(7) * t164 + t32) * t119) * t174, t43 * t71 + t44 * t69 + ((-t51 + t183) * pkin(7) + t209) * t121 + (-t14 + t135 * t32 + (-t50 + t210) * pkin(7)) * t118, -t80 * pkin(2) - t31 * t43 - t32 * t44 - t57 * t87 + (-t14 * t118 - t127 * t121 + (-t118 * t32 - t121 * t31) * qJD(3)) * pkin(7), -t21 * t120 * t118 + (-t118 * t168 + t120 * t163 - t61) * t49, t61 * t47 + t49 * t60 + (-t117 * t49 - t120 * t47) * t163 + (t186 - t185 + (t117 * t47 - t120 * t49) * qJD(4)) * t118, -t61 * t66 + (t169 * t66 + t21) * t121 + (-t135 * t49 - t168 * t66 + t187) * t118, t22 * t117 * t118 + (t117 * t163 + t118 * t167 - t60) * t47, t60 * t66 + (-qJD(3) * t188 + t22) * t121 + (t135 * t47 - t167 * t66 - t189) * t118, -t126 * t66 - t184, -t29 * t60 - t36 * t47 + t76 * t51 - t192 * t66 + (-t2 + (pkin(7) * t47 + t117 * t29) * qJD(3)) * t121 + (pkin(7) * t22 + t8 * t117 + t134 * t135 + t167 * t29) * t118, -t29 * t61 - t36 * t49 - t77 * t51 + t191 * t66 + (t1 + (pkin(7) * t49 + t120 * t29) * qJD(3)) * t121 + (-pkin(7) * t21 + t8 * t120 + t135 * t6 - t168 * t29) * t118, t76 * t21 - t77 * t22 - t134 * t61 + t6 * t60 + t192 * t49 + t191 * t47 + t137 * t163 + (-t1 * t117 - t120 * t2 + (-t117 * t134 - t120 * t6) * qJD(4)) * t118, t1 * t77 + t2 * t76 - t29 * t36 - t191 * t6 + t192 * t134 + (t118 * t8 + t163 * t29) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t69 ^ 2 + t71 ^ 2, -t50 - t203, -t193, -t71 * t107 - t130, t141, -t32 * t107 - t57 * t71 - t152, t57 * t69 - t209, 0, 0, t151 * t49 - t186, (-t21 - t197) * t120 + (-t22 - t195) * t117, t151 * t66 - t49 * t71 + t189, t188 * t47 - t185, -t117 * t66 ^ 2 + t47 * t71 + t187, -t66 * t71, -pkin(3) * t22 - t11 * t66 + t117 * t128 - t120 * t143 + t134 * t71 - t32 * t47, pkin(3) * t21 + t117 * t143 + t12 * t66 + t120 * t128 - t32 * t49 + t6 * t71, t11 * t49 + t12 * t47 + ((-t22 + t182) * pkin(8) + t212) * t120 + ((qJD(4) * t47 - t21) * pkin(8) + t205) * t117, -t8 * pkin(3) + t134 * t11 - t6 * t12 - t29 * t32 + (qJD(4) * t137 + t1 * t120 - t2 * t117) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, -t47 ^ 2 + t49 ^ 2, -t21 + t197, -t196, t195 - t22, t51, -t29 * t49 - t205, t29 * t47 - t212, 0, 0;];
tauc_reg = t5;
