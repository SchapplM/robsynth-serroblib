% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:35
% DurationCPUTime: 2.29s
% Computational Cost: add. (2292->314), mult. (6073->445), div. (0->0), fcn. (4279->8), ass. (0->176)
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t154 = qJD(4) * t109;
t155 = qJD(4) * t106;
t110 = cos(qJ(3));
t105 = cos(pkin(5));
t111 = cos(qJ(2));
t104 = sin(pkin(5));
t161 = qJD(2) * t104;
t144 = t111 * t161;
t117 = qJD(1) * (qJD(3) * t105 + t144);
t107 = sin(qJ(3));
t157 = qJD(3) * t107;
t108 = sin(qJ(2));
t163 = qJD(1) * t104;
t140 = t108 * t163;
t88 = qJD(2) * pkin(7) + t140;
t39 = t110 * t117 - t157 * t88;
t162 = qJD(1) * t105;
t138 = t107 * t162;
t62 = t110 * t88 + t138;
t54 = qJD(3) * pkin(8) + t62;
t126 = pkin(3) * t107 - pkin(8) * t110;
t86 = t126 * qJD(3);
t60 = (t86 + t140) * qJD(2);
t139 = t111 * t163;
t90 = -pkin(3) * t110 - pkin(8) * t107 - pkin(2);
t64 = qJD(2) * t90 - t139;
t134 = -t106 * t60 - t109 * t39 - t64 * t154 + t54 * t155;
t25 = -t106 * t54 + t109 * t64;
t159 = qJD(2) * t110;
t95 = -qJD(4) + t159;
t214 = t25 * t95 - t134;
t26 = t106 * t64 + t109 * t54;
t9 = -qJD(4) * t26 - t106 * t39 + t109 * t60;
t213 = t26 * t95 - t9;
t167 = t110 * t111;
t205 = pkin(7) * t106;
t212 = t109 * t86 + t157 * t205 - (-t106 * t167 + t108 * t109) * t163;
t211 = t106 * t86 + t90 * t154 - (t106 * t108 + t109 * t167) * t163;
t153 = t109 * qJD(3);
t160 = qJD(2) * t107;
t81 = t106 * t160 - t153;
t200 = t81 * t95;
t142 = t107 * t155;
t143 = t110 * t153;
t150 = qJD(3) * qJD(4);
t56 = -t109 * t150 + (t142 - t143) * qJD(2);
t210 = -t56 + t200;
t158 = qJD(3) * t106;
t83 = t109 * t160 + t158;
t198 = t83 * t95;
t141 = t107 * t154;
t156 = qJD(3) * t110;
t118 = t106 * t156 + t141;
t57 = qJD(2) * t118 + t106 * t150;
t209 = t57 - t198;
t61 = -t107 * t88 + t110 * t162;
t208 = t83 ^ 2;
t207 = pkin(4) * t81;
t206 = pkin(4) * t106;
t20 = -qJ(5) * t81 + t26;
t204 = t20 * t95;
t40 = t107 * t117 + t156 * t88;
t173 = t104 * t108;
t70 = -t105 * t110 + t107 * t173;
t201 = t40 * t70;
t199 = t83 * t81;
t197 = -qJ(5) - pkin(8);
t19 = -qJ(5) * t83 + t25;
t12 = -pkin(4) * t95 + t19;
t196 = t12 - t19;
t168 = t109 * t110;
t122 = pkin(4) * t107 - qJ(5) * t168;
t152 = t109 * qJD(5);
t174 = qJ(5) * t107;
t96 = pkin(7) * t168;
t195 = -t107 * t152 + t122 * qJD(3) + (-t96 + (-t90 + t174) * t106) * qJD(4) + t212;
t169 = t107 * t109;
t194 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t169 + (-qJD(5) * t107 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t110) * t106 + t211;
t135 = qJD(4) * t197;
t85 = t126 * qJD(2);
t34 = -t106 * t61 + t109 * t85;
t193 = qJD(2) * t122 + t106 * qJD(5) - t109 * t135 + t34;
t35 = t106 * t85 + t109 * t61;
t192 = -t152 + t35 + (-qJ(5) * t159 - t135) * t106;
t191 = (-t107 * t153 - t110 * t155) * pkin(7) + t211;
t67 = t106 * t90 + t96;
t190 = -t67 * qJD(4) + t212;
t187 = qJD(2) * pkin(2);
t23 = t57 * pkin(4) + t40;
t186 = t106 * t23;
t185 = t106 * t40;
t53 = -qJD(3) * pkin(3) - t61;
t184 = t106 * t53;
t183 = t106 * t83;
t182 = t106 * t95;
t181 = t109 * t23;
t180 = t109 * t40;
t179 = t109 * t53;
t178 = t109 * t95;
t177 = t40 * t107;
t176 = t56 * t106;
t175 = t57 * t109;
t172 = t104 * t111;
t113 = qJD(2) ^ 2;
t171 = t104 * t113;
t170 = t106 * t110;
t112 = qJD(3) ^ 2;
t166 = t112 * t107;
t165 = t112 * t110;
t102 = t107 ^ 2;
t103 = t110 ^ 2;
t164 = t102 - t103;
t151 = qJD(2) * qJD(3);
t149 = t108 * t171;
t148 = t107 * t113 * t110;
t147 = t95 * t160;
t145 = t108 * t161;
t98 = t107 * t151;
t136 = -qJD(5) - t207;
t133 = pkin(4) * t98;
t132 = t81 * t139;
t131 = t83 * t139;
t130 = t107 * t139;
t129 = t107 * t144;
t128 = t110 * t144;
t127 = t110 * t98;
t89 = -t139 - t187;
t125 = -t89 - t139;
t124 = -t106 * t26 - t109 * t25;
t123 = qJD(2) * t102 - t110 * t95;
t121 = qJ(5) * t57 + t134;
t71 = t105 * t107 + t110 * t173;
t46 = -t106 * t71 - t109 * t172;
t120 = t106 * t172 - t109 * t71;
t119 = qJD(2) * t125;
t116 = qJD(3) * (-t125 - t187);
t115 = t177 + t39 * t110 + (-t107 * t62 - t110 * t61) * qJD(3);
t114 = qJ(5) * t56 + t9;
t100 = -pkin(4) * t109 - pkin(3);
t92 = t197 * t109;
t91 = t197 * t106;
t87 = (pkin(7) + t206) * t107;
t80 = t109 * t90;
t78 = t81 ^ 2;
t66 = -pkin(7) * t170 + t80;
t65 = (-t95 - t159) * t157;
t63 = pkin(4) * t118 + pkin(7) * t156;
t48 = -t106 * t174 + t67;
t45 = qJD(3) * t71 + t129;
t44 = -qJD(3) * t70 + t128;
t42 = t138 + (qJD(2) * t206 + t88) * t110;
t41 = -qJ(5) * t169 + t80 + (-pkin(4) - t205) * t110;
t38 = -t78 + t208;
t36 = -t136 + t53;
t33 = -t198 - t57;
t32 = -t56 - t200;
t28 = -t95 * t154 + (t95 * t168 + (-t83 + t158) * t107) * qJD(2);
t27 = t95 * t155 + (-t95 * t170 + (t81 + t153) * t107) * qJD(2);
t22 = -t182 * t81 - t175;
t21 = -t178 * t83 - t176;
t18 = t81 * t141 + (t107 * t57 + t156 * t81) * t106;
t17 = t83 * t143 + (-t56 * t109 - t155 * t83) * t107;
t15 = qJD(4) * t46 + t106 * t145 + t44 * t109;
t14 = qJD(4) * t120 - t44 * t106 + t109 * t145;
t11 = t95 * t141 + t57 * t110 + (-t123 * t106 - t107 * t81) * qJD(3);
t10 = t95 * t142 + t56 * t110 + (t107 * t83 + t109 * t123) * qJD(3);
t7 = -t209 * t106 + t210 * t109;
t6 = (-t109 * t81 - t183) * t156 + (t176 - t175 + (t106 * t81 - t109 * t83) * qJD(4)) * t107;
t5 = -qJD(5) * t81 - t121;
t4 = -qJD(5) * t83 + t114 + t133;
t3 = t120 * t98 + t15 * t95 + t45 * t83 - t56 * t70;
t2 = -t14 * t95 + t45 * t81 + t46 * t98 + t57 * t70;
t1 = t120 * t57 - t14 * t83 - t15 * t81 + t46 * t56;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t111 * t171, 0, 0, 0, 0, 0, 0, 0, 0, -t110 * t149 + (-t45 - t129) * qJD(3), t107 * t149 + (-t44 - t128) * qJD(3), (t107 * t45 + t110 * t44 + (-t107 * t71 + t110 * t70) * qJD(3)) * qJD(2), t39 * t71 + t201 + t44 * t62 - t45 * t61 + (t89 - t139) * t145, 0, 0, 0, 0, 0, 0, t2, t3, t1, t120 * t134 + t14 * t25 + t15 * t26 + t45 * t53 + t46 * t9 + t201, 0, 0, 0, 0, 0, 0, t2, t3, t1, t12 * t14 - t120 * t5 + t15 * t20 + t23 * t70 + t36 * t45 + t4 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t127, -0.2e1 * t164 * t151, t165, -0.2e1 * t127, -t166, 0, -pkin(7) * t165 + t107 * t116, pkin(7) * t166 + t110 * t116, (-t102 - t103) * qJD(2) * t139 + t115, ((t107 * t61 - t110 * t62) * t111 + (-t89 - t187) * t108) * t163 + t115 * pkin(7), t17, t6, t10, t18, t11, t65, -t190 * t95 + (-t9 + (pkin(7) * t81 + t184) * qJD(3)) * t110 + (-t132 + t53 * t154 + pkin(7) * t57 + t185 + (qJD(2) * t66 + t25) * qJD(3)) * t107, t191 * t95 + (-t134 + (pkin(7) * t83 + t179) * qJD(3)) * t110 + (-t131 - t53 * t155 - pkin(7) * t56 + t180 + (-qJD(2) * t67 - t26) * qJD(3)) * t107, t56 * t66 - t57 * t67 - t190 * t83 - t191 * t81 + t124 * t156 + (t106 * t134 - t109 * t9 + (t106 * t25 - t109 * t26) * qJD(4)) * t107, -t53 * t130 + t66 * t9 - t67 * t134 + t191 * t26 + t190 * t25 + (t156 * t53 + t177) * pkin(7), t17, t6, t10, t18, t11, t65, t57 * t87 + t63 * t81 - t195 * t95 + (t158 * t36 - t4) * t110 + (-t132 + t36 * t154 + t186 + (qJD(2) * t41 + t12) * qJD(3)) * t107, -t56 * t87 + t63 * t83 + t194 * t95 + (t153 * t36 + t5) * t110 + (-t131 - t36 * t155 + t181 + (-qJD(2) * t48 - t20) * qJD(3)) * t107, t41 * t56 - t48 * t57 - t195 * t83 - t194 * t81 + (-t106 * t20 - t109 * t12) * t156 + (-t106 * t5 - t109 * t4 + (t106 * t12 - t109 * t20) * qJD(4)) * t107, t23 * t87 + t4 * t41 + t48 * t5 + (t63 - t130) * t36 + t194 * t20 + t195 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t164 * t113, 0, t148, 0, 0, t107 * t119, t110 * t119, 0, 0, t21, t7, t28, t22, t27, t147, -pkin(3) * t57 - t180 + t34 * t95 - t62 * t81 + (pkin(8) * t178 + t184) * qJD(4) + (-t107 * t25 + (-pkin(8) * t157 - t110 * t53) * t106) * qJD(2), pkin(3) * t56 + t185 - t35 * t95 - t62 * t83 + (-pkin(8) * t182 + t179) * qJD(4) + (-t53 * t168 + (-pkin(8) * t153 + t26) * t107) * qJD(2), t34 * t83 + t35 * t81 + ((qJD(4) * t83 - t57) * pkin(8) + t214) * t109 + ((qJD(4) * t81 - t56) * pkin(8) + t213) * t106, -pkin(3) * t40 - t25 * t34 - t26 * t35 - t53 * t62 + (qJD(4) * t124 - t106 * t9 - t109 * t134) * pkin(8), t21, t7, t28, t22, t27, t147, t100 * t57 - t181 - t42 * t81 + t193 * t95 + (t36 + t207) * t155 + (-t36 * t170 + (qJD(3) * t91 - t12) * t107) * qJD(2), -t100 * t56 + t186 - t42 * t83 - t192 * t95 + (pkin(4) * t183 + t109 * t36) * qJD(4) + (-t36 * t168 + (qJD(3) * t92 + t20) * t107) * qJD(2), t56 * t91 + t57 * t92 + t193 * t83 + t192 * t81 + (t12 * t95 + t5) * t109 + (-t4 + t204) * t106, t100 * t23 + t4 * t91 - t5 * t92 + (pkin(4) * t155 - t42) * t36 - t192 * t20 - t193 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t38, t32, -t199, t33, t98, -t53 * t83 - t213, t53 * t81 - t214, 0, 0, t199, t38, t32, -t199, t33, t98, 0.2e1 * t133 - t204 + (t136 - t36) * t83 + t114, -pkin(4) * t208 - t19 * t95 + (qJD(5) + t36) * t81 + t121, t56 * pkin(4) - t196 * t81, t196 * t20 + (-t36 * t83 + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, t210, -t78 - t208, t12 * t83 + t20 * t81 + t23;];
tauc_reg = t8;
