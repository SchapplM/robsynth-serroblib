% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:46
% DurationCPUTime: 1.97s
% Computational Cost: add. (4118->291), mult. (9565->397), div. (0->0), fcn. (6377->8), ass. (0->173)
t109 = cos(qJ(3));
t200 = cos(qJ(4));
t148 = t200 * t109;
t106 = sin(qJ(4));
t107 = sin(qJ(3));
t170 = t106 * t107;
t120 = t148 - t170;
t207 = qJD(1) * t120;
t163 = qJD(3) * t107;
t155 = pkin(3) * t163;
t206 = 0.2e1 * t155;
t157 = qJD(3) + qJD(4);
t135 = qJD(3) * t148;
t145 = t200 * qJD(4);
t59 = -t109 * t145 + t157 * t170 - t135;
t205 = t59 * t157;
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t94 = sin(pkin(9)) * pkin(1) + pkin(6);
t87 = t94 * qJD(1);
t143 = pkin(7) * qJD(1) + t87;
t159 = t107 * qJD(2);
t69 = t143 * t109 + t159;
t152 = t200 * t69;
t100 = t109 * qJD(2);
t127 = t143 * t107;
t68 = t100 - t127;
t64 = qJD(3) * pkin(3) + t68;
t39 = t106 * t64 + t152;
t36 = t157 * pkin(8) + t39;
t165 = qJD(1) * t107;
t78 = -qJD(1) * t148 + t106 * t165;
t169 = t106 * t109;
t85 = t200 * t107 + t169;
t80 = qJD(1) * t85;
t95 = -cos(pkin(9)) * pkin(1) - pkin(2);
t86 = -t109 * pkin(3) + t95;
t81 = qJD(1) * t86;
t45 = t78 * pkin(4) - t80 * pkin(8) + t81;
t123 = t105 * t36 - t108 * t45;
t17 = t105 * t45 + t108 * t36;
t204 = t105 * t123 + t108 * t17;
t202 = t157 * qJD(1);
t54 = t120 * t202;
t67 = t105 * t157 + t108 * t80;
t32 = t67 * qJD(5) + t105 * t54;
t60 = t157 * t85;
t55 = t60 * qJD(1);
t46 = t85 * t55;
t77 = qJD(5) + t78;
t128 = -t59 * t77 + t46;
t161 = qJD(5) * t105;
t151 = t85 * t161;
t203 = -t128 * t108 + t77 * t151;
t117 = t69 * qJD(3);
t162 = qJD(4) * t106;
t96 = qJD(3) * t100;
t62 = -qJD(3) * t127 + t96;
t132 = -t69 * t162 + t200 * t62;
t14 = -t106 * t117 + t64 * t145 + t132;
t158 = qJD(1) * qJD(3);
t144 = t107 * t158;
t24 = pkin(3) * t144 + t55 * pkin(4) - t54 * pkin(8);
t2 = -t123 * qJD(5) + t105 * t24 + t108 * t14;
t201 = pkin(7) + t94;
t82 = t201 * t107;
t83 = t201 * t109;
t121 = -t106 * t83 - t200 * t82;
t58 = t200 * t117;
t142 = t106 * t62 + t58;
t15 = t39 * qJD(4) + t142;
t199 = t15 * t121;
t198 = t15 * t120;
t197 = t15 * t85;
t1 = t2 * t108;
t196 = t55 * t120;
t140 = t108 * t157;
t65 = t105 * t80 - t140;
t195 = t59 * t65;
t194 = t59 * t67;
t193 = t65 * t77;
t192 = t67 * t65;
t191 = t67 * t77;
t190 = t77 * t80;
t189 = t80 * t78;
t188 = t81 * t80;
t172 = t32 * t108;
t187 = t108 * t195 - t85 * t172;
t31 = -qJD(5) * t140 - t108 * t54 + t80 * t161;
t186 = t120 * t31 + t67 * t60;
t185 = t59 * t78 - t46;
t184 = pkin(3) * qJD(4);
t181 = t105 * t55;
t180 = t105 * t65;
t179 = t105 * t78;
t178 = t106 * t69;
t177 = t107 * t87;
t175 = t108 * t55;
t174 = t108 * t78;
t173 = t31 * t105;
t88 = qJD(1) * t95;
t171 = qJD(5) * t65;
t110 = qJD(3) ^ 2;
t168 = t110 * t107;
t167 = t110 * t109;
t166 = t107 ^ 2 - t109 ^ 2;
t164 = qJD(1) * t109;
t160 = qJD(5) * t108;
t156 = t123 * t174 - t17 * t179 + t1;
t154 = pkin(3) * t165;
t153 = t200 * t64;
t150 = t85 * t160;
t111 = qJD(1) ^ 2;
t149 = t107 * t111 * t109;
t38 = t153 - t178;
t35 = -t157 * pkin(4) - t38;
t33 = t35 * t161;
t34 = t35 * t160;
t147 = t123 * t80 + t33;
t146 = qJD(3) * t201;
t141 = t108 * t77;
t138 = t67 * t150;
t137 = pkin(3) * t145;
t136 = t15 * t105 + t17 * t80 + t34;
t134 = t109 * t144;
t43 = t106 * t68 + t152;
t133 = pkin(3) * t162 - t43;
t56 = t80 * pkin(4) + t78 * pkin(8);
t131 = -t35 * t59 + t197;
t130 = t120 * t32 - t60 * t65;
t129 = -t120 * t54 + t80 * t60;
t125 = t105 * t17 - t108 * t123;
t51 = -pkin(4) * t120 - t85 * pkin(8) + t86;
t53 = -t106 * t82 + t200 * t83;
t25 = -t105 * t53 + t108 * t51;
t26 = t105 * t51 + t108 * t53;
t122 = 0.2e1 * qJD(3) * t88;
t74 = t109 * t87 + t159;
t119 = t81 * t78 - t132;
t3 = -qJD(5) * t17 - t105 * t14 + t108 * t24;
t118 = -t125 * qJD(5) - t3 * t105;
t98 = t106 * pkin(3) + pkin(8);
t116 = -t77 * t137 + t35 * t78 - t55 * t98;
t115 = t106 * (-pkin(7) * t164 - t74);
t114 = -t128 * t105 - t77 * t150;
t113 = t118 + t1;
t70 = -t87 * t163 + t96;
t71 = t74 * qJD(3);
t73 = t100 - t177;
t112 = t71 * t107 + t70 * t109 + (-t107 * t74 - t109 * t73) * qJD(3);
t99 = -t200 * pkin(3) - pkin(4);
t76 = t107 * t146;
t57 = t60 * t157;
t49 = t56 + t154;
t47 = -t78 ^ 2 + t80 ^ 2;
t44 = t200 * t68 - t178;
t42 = t80 * t157 - t85 * t202;
t41 = (t78 + t207) * t157;
t30 = t60 * pkin(4) + t59 * pkin(8) + t155;
t29 = t53 * qJD(4) - t106 * t76 + t201 * t135;
t28 = t121 * qJD(4) - t146 * t169 - t200 * t76;
t21 = t105 * t56 + t108 * t38;
t20 = -t105 * t38 + t108 * t56;
t19 = t105 * t49 + t108 * t44;
t18 = -t105 * t44 + t108 * t49;
t10 = t77 * t141 - t67 * t80 + t181;
t9 = -t77 ^ 2 * t105 + t65 * t80 + t175;
t8 = t77 * t180 - t172;
t7 = t67 * t141 - t173;
t6 = -t26 * qJD(5) - t105 * t28 + t108 * t30;
t5 = t25 * qJD(5) + t105 * t30 + t108 * t28;
t4 = (-t31 - t193) * t108 + (-t32 - t191) * t105;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t134, -0.2e1 * t166 * t158, t167, -0.2e1 * t134, -t168, 0, t107 * t122 - t94 * t167, t109 * t122 + t94 * t168, t112, t112 * t94, t54 * t85 - t80 * t59, -t129 + t185, -t205, t78 * t60 - t196, -t57, 0, t86 * t55 + t81 * t60 - t29 * t157 + (t78 - t207) * t155, -t28 * t157 + t80 * t206 + t86 * t54 - t81 * t59, t120 * t14 - t121 * t54 - t28 * t78 + t29 * t80 + t38 * t59 - t39 * t60 - t53 * t55 + t197, t14 * t53 + t81 * t206 + t39 * t28 - t38 * t29 - t199, -t67 * t151 + (-t31 * t85 - t194) * t108, -t138 + (t194 + (t31 + t171) * t85) * t105 + t187, t186 - t203, t65 * t150 + (t32 * t85 - t195) * t105, t114 + t130, t77 * t60 - t196, t131 * t105 - t120 * t3 - t121 * t32 - t123 * t60 + t25 * t55 + t29 * t65 + t85 * t34 + t6 * t77, t108 * t131 + t120 * t2 + t121 * t31 - t17 * t60 - t26 * t55 + t29 * t67 - t85 * t33 - t5 * t77, t25 * t31 - t26 * t32 - t5 * t65 - t6 * t67 + t125 * t59 + (-qJD(5) * t204 - t105 * t2 - t108 * t3) * t85, -t123 * t6 + t17 * t5 + t2 * t26 + t3 * t25 + t35 * t29 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, -t167, 0, t70 * t107 - t71 * t109 + (-t107 * t73 + t109 * t74) * qJD(3), 0, 0, 0, 0, 0, 0, -t57, t205, t129 + t185, t14 * t85 - t38 * t60 - t39 * t59 - t198, 0, 0, 0, 0, 0, 0, t114 - t130, t186 + t203, t138 + (-t194 + (-t31 + t171) * t85) * t105 + t187, t113 * t85 - t204 * t59 + t35 * t60 - t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t166 * t111, 0, t149, 0, 0, -t88 * t165, -t88 * t164 - t96 + (t73 + t177) * qJD(3), 0, 0, t189, t47, t41, -t189, t42, 0, -t69 * t145 - t58 - t78 * t154 - t188 + t43 * t157 + (-qJD(4) * t64 - t157 * t184 - t62) * t106, (-t153 + t44) * qJD(4) + (-t115 + t44) * qJD(3) + (-t157 * t145 - t80 * t165) * pkin(3) + t119, (t39 - t43) * t80 + (-t38 + t44) * t78 + (-t200 * t54 - t106 * t55 + (t106 * t80 - t200 * t78) * qJD(4)) * pkin(3), t38 * t43 - t39 * t44 + (-t81 * t165 - t200 * t15 + t106 * t14 + (-t106 * t38 + t200 * t39) * qJD(4)) * pkin(3), t7, t4, t10, t8, t9, -t190, -t18 * t77 + t99 * t32 + t133 * t65 + (-qJD(5) * t77 * t98 - t15) * t108 + t116 * t105 + t147, -t99 * t31 + (t98 * t161 + t19) * t77 + t133 * t67 + t116 * t108 + t136, t18 * t67 + t19 * t65 + (-t65 * t137 - t32 * t98 + (t67 * t98 + t123) * qJD(5)) * t108 + (t67 * t137 - t31 * t98 - t3 + (t65 * t98 - t17) * qJD(5)) * t105 + t156, t15 * t99 + t123 * t18 - t17 * t19 - t35 * t43 + (t106 * t35 + t204 * t200) * t184 + t113 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t47, t41, -t189, t42, 0, t39 * qJD(3) - t142 - t188, (-t153 + t38) * qJD(4) + (-t115 + t38) * qJD(3) + t119, 0, 0, t7, t4, t10, t8, t9, -t190, t35 * t179 - pkin(4) * t32 - t15 * t108 - t20 * t77 - t39 * t65 + (-t77 * t160 - t181) * pkin(8) + t147, t35 * t174 + pkin(4) * t31 + t21 * t77 - t39 * t67 + (t77 * t161 - t175) * pkin(8) + t136, t20 * t67 + t21 * t65 + (-t173 - t172 + (t108 * t67 + t180) * qJD(5)) * pkin(8) + t118 + t156, -t15 * pkin(4) + pkin(8) * t113 + t123 * t20 - t17 * t21 - t35 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, -t65 ^ 2 + t67 ^ 2, -t31 + t193, -t192, t191 - t32, t55, t17 * t77 - t35 * t67 + t3, -t123 * t77 + t35 * t65 - t2, 0, 0;];
tauc_reg = t11;
