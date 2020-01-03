% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:08
% EndTime: 2019-12-31 22:02:16
% DurationCPUTime: 2.46s
% Computational Cost: add. (2741->258), mult. (6918->460), div. (0->0), fcn. (5813->6), ass. (0->143)
t102 = sin(qJ(4));
t105 = cos(qJ(3));
t175 = -pkin(8) - pkin(7);
t180 = t175 * t105;
t184 = t102 * t180;
t104 = sin(qJ(2));
t103 = sin(qJ(3));
t106 = cos(qJ(2));
t157 = t106 * qJD(2);
t138 = t103 * t157;
t160 = qJD(3) * t105;
t183 = t104 * t160 + t138;
t182 = -0.4e1 * t104;
t174 = cos(qJ(4));
t134 = qJD(2) * t174;
t124 = t106 * t134;
t161 = qJD(3) * t103;
t142 = t104 * t161;
t163 = t103 * t104;
t146 = t102 * t163;
t132 = t174 * qJD(4);
t178 = t174 * qJD(3) + t132;
t27 = t103 * t124 - t102 * t142 - qJD(4) * t146 + (t102 * t157 + t178 * t104) * t105;
t70 = t102 * t105 + t174 * t103;
t59 = t70 * t104;
t181 = t27 * qJ(5) + t59 * qJD(5);
t164 = t102 * t103;
t179 = t175 * t164;
t99 = t104 ^ 2;
t130 = (t106 ^ 2 - t99) * qJD(2);
t100 = t105 ^ 2;
t98 = t103 ^ 2;
t166 = t100 - t98;
t131 = t166 * qJD(3);
t177 = qJD(3) + qJD(4);
t159 = qJD(3) * t106;
t141 = t103 * t159;
t96 = t104 * qJD(2);
t116 = t105 * t96 + t141;
t170 = t106 * pkin(2);
t123 = -t104 * pkin(7) - t170;
t119 = -pkin(1) + t123;
t117 = t105 * t119;
t122 = pkin(2) * t104 - pkin(7) * t106;
t118 = t122 * t103;
t32 = t116 * pkin(6) - qJD(2) * t118 - qJD(3) * t117;
t176 = t175 * t104 - pkin(1) - t170;
t173 = pkin(3) * t102;
t172 = t103 * pkin(6);
t97 = t104 * pkin(6);
t171 = t105 * pkin(2);
t128 = pkin(3) * t132;
t169 = -t59 * t128 - t27 * t173;
t47 = t177 * t70;
t143 = t174 * t105;
t69 = -t143 + t164;
t168 = -t69 * t128 - t47 * t173;
t144 = -pkin(3) - t172;
t112 = t176 * t105 + t144 * t106;
t40 = t102 * t112;
t162 = t105 * t106;
t88 = pkin(6) * t162;
t55 = t103 * t119 + t88;
t48 = -pkin(8) * t163 + t55;
t22 = t174 * t48 + t40;
t50 = -t174 * t180 + t179;
t74 = pkin(3) * t163 + t97;
t158 = qJD(4) * t102;
t156 = -0.2e1 * pkin(1) * qJD(2);
t155 = -0.2e1 * pkin(2) * qJD(3);
t109 = (-t176 * t103 - t88) * qJD(3) + (t106 * t180 + (-t144 + t171) * t104) * qJD(2);
t111 = -pkin(8) * t183 - t32;
t41 = t174 * t112;
t154 = -qJD(4) * t41 - t102 * t109 - t174 * t111;
t94 = pkin(6) * t157;
t53 = t183 * pkin(3) + t94;
t153 = t106 * t172;
t152 = t174 * pkin(3);
t151 = pkin(3) * t161;
t150 = pkin(3) * t158;
t60 = t104 * t143 - t146;
t148 = t60 * t158;
t147 = t70 * t158;
t93 = -t105 * pkin(3) - pkin(2);
t139 = t105 * t159;
t137 = t103 * t160;
t136 = t104 * t157;
t135 = t105 * t157;
t21 = -t102 * t48 + t41;
t126 = t175 * t174;
t72 = t103 * t126;
t49 = t72 + t184;
t129 = 0.2e1 * t136;
t127 = t99 * t137;
t125 = t103 * t135;
t54 = t117 - t153;
t121 = -t103 * t55 - t105 * t54;
t120 = qJD(3) * t126;
t4 = t48 * t158 + t154;
t28 = -qJD(3) * t184 - qJD(4) * t72 - t103 * t120 - t158 * t180;
t114 = t106 * t128 + (-pkin(3) * t96 + qJD(4) * t48) * t102 + t154;
t33 = -t55 * qJD(3) + (pkin(6) * t163 + t105 * t122) * qJD(2);
t113 = t121 * qJD(3) - t33 * t103 - t32 * t105;
t29 = t105 * t120 + t132 * t180 - t177 * t179;
t110 = -t102 * t111 + t174 * t109;
t108 = -qJD(4) * t40 - t48 * t132 + t110;
t26 = t102 * t138 + t47 * t104 - t105 * t124;
t107 = t26 * qJ(5) - t60 * qJD(5) + t108;
t95 = pkin(4) * t96;
t1 = t95 + t107;
t92 = t152 + pkin(4);
t90 = -0.2e1 * t128;
t89 = -0.2e1 * t150;
t84 = -0.2e1 * t136;
t81 = t106 * t150;
t56 = t69 * pkin(4) + t93;
t51 = -t104 * t131 - t125;
t46 = -t178 * t105 + t177 * t164;
t44 = t59 * pkin(4) + t74;
t38 = t47 * pkin(4) + t151;
t37 = -t69 * qJ(5) + t50;
t36 = -t70 * qJ(5) + t49;
t35 = -0.2e1 * t70 * t46;
t34 = 0.2e1 * t69 * t47;
t31 = t47 * t106 - t69 * t96;
t30 = t46 * t106 + t70 * t96;
t18 = -0.2e1 * t60 * t26;
t17 = 0.2e1 * t59 * t27;
t16 = t27 * pkin(4) + t53;
t15 = -t59 * qJ(5) + t22;
t14 = 0.2e1 * t27 * t106 - 0.2e1 * t59 * t96;
t13 = 0.2e1 * t26 * t106 + 0.2e1 * t60 * t96;
t12 = 0.2e1 * t46 * t69 - 0.2e1 * t70 * t47;
t11 = -t106 * pkin(4) - t60 * qJ(5) + t21;
t10 = t46 * qJ(5) - t70 * qJD(5) + t29;
t9 = t47 * qJ(5) + t69 * qJD(5) + t28;
t8 = t27 * t69 + t59 * t47;
t7 = -t26 * t70 - t60 * t46;
t6 = 0.2e1 * t26 * t59 - 0.2e1 * t60 * t27;
t5 = -t22 * qJD(4) + t110;
t3 = t26 * t69 - t70 * t27 + t46 * t59 - t60 * t47;
t2 = t4 + t181;
t19 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0.2e1 * t130, 0, t84, 0, 0, t104 * t156, t106 * t156, 0, 0, 0.2e1 * t100 * t136 - 0.2e1 * t127, t125 * t182 - 0.2e1 * t99 * t131, 0.2e1 * t104 * t141 - 0.2e1 * t105 * t130, 0.2e1 * t98 * t136 + 0.2e1 * t127, 0.2e1 * t103 * t130 + 0.2e1 * t104 * t139, t84, 0.2e1 * t54 * t96 - 0.2e1 * t33 * t106 + 0.2e1 * (t103 * t129 + t99 * t160) * pkin(6), -0.2e1 * t55 * t96 - 0.2e1 * t32 * t106 + 0.2e1 * (t105 * t129 - t99 * t161) * pkin(6), 0.2e1 * t121 * t157 + 0.2e1 * (t103 * t32 - t105 * t33 + (t103 * t54 - t105 * t55) * qJD(3)) * t104, 0.2e1 * pkin(6) ^ 2 * t136 - 0.2e1 * t55 * t32 + 0.2e1 * t54 * t33, t18, t6, t13, t17, t14, t84, -0.2e1 * t5 * t106 + 0.2e1 * t21 * t96 + 0.2e1 * t74 * t27 + 0.2e1 * t53 * t59, -0.2e1 * t4 * t106 - 0.2e1 * t22 * t96 - 0.2e1 * t74 * t26 + 0.2e1 * t53 * t60, 0.2e1 * t21 * t26 - 0.2e1 * t22 * t27 + 0.2e1 * t4 * t59 - 0.2e1 * t5 * t60, 0.2e1 * t21 * t5 - 0.2e1 * t22 * t4 + 0.2e1 * t74 * t53, t18, t6, t13, t17, t14, t84, -0.2e1 * t1 * t106 + 0.2e1 * t11 * t96 + 0.2e1 * t16 * t59 + 0.2e1 * t44 * t27, -0.2e1 * t2 * t106 - 0.2e1 * t15 * t96 + 0.2e1 * t16 * t60 - 0.2e1 * t44 * t26, -0.2e1 * t1 * t60 + 0.2e1 * t11 * t26 - 0.2e1 * t15 * t27 + 0.2e1 * t2 * t59, 0.2e1 * t11 * t1 - 0.2e1 * t15 * t2 + 0.2e1 * t44 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, 0, -t96, 0, -t94, pkin(6) * t96, 0, 0, -t51, t137 * t182 + t166 * t157, t103 * t96 - t139, t51, t116, 0, (pkin(7) * t162 + (-t171 + t172) * t104) * qJD(3) + (t103 * t123 - t88) * qJD(2), (t105 * t97 + t118) * qJD(3) + (t105 * t123 + t153) * qJD(2), t113, -pkin(2) * t94 + pkin(7) * t113, t7, t3, t30, t8, t31, 0, -t29 * t106 + t59 * t151 + t93 * t27 + t74 * t47 + t49 * t96 + t53 * t69, -t28 * t106 + t60 * t151 - t93 * t26 - t74 * t46 - t50 * t96 + t53 * t70, t21 * t46 - t22 * t47 + t49 * t26 - t50 * t27 + t28 * t59 - t29 * t60 + t4 * t69 - t5 * t70, t74 * t151 + t21 * t29 - t22 * t28 - t4 * t50 + t5 * t49 + t53 * t93, t7, t3, t30, t8, t31, 0, -t10 * t106 + t16 * t69 + t56 * t27 + t36 * t96 + t38 * t59 + t44 * t47, -t9 * t106 + t16 * t70 - t56 * t26 - t37 * t96 + t38 * t60 - t44 * t46, -t1 * t70 - t10 * t60 + t11 * t46 - t15 * t47 + t2 * t69 + t36 * t26 - t37 * t27 + t9 * t59, t1 * t36 + t11 * t10 - t15 * t9 + t16 * t56 - t2 * t37 + t44 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t137, 0.2e1 * t131, 0, -0.2e1 * t137, 0, 0, t103 * t155, t105 * t155, 0, 0, t35, t12, 0, t34, 0, 0, 0.2e1 * t69 * t151 + 0.2e1 * t93 * t47, 0.2e1 * t70 * t151 - 0.2e1 * t93 * t46, 0.2e1 * t28 * t69 - 0.2e1 * t29 * t70 + 0.2e1 * t49 * t46 - 0.2e1 * t50 * t47, 0.2e1 * t93 * t151 - 0.2e1 * t50 * t28 + 0.2e1 * t49 * t29, t35, t12, 0, t34, 0, 0, 0.2e1 * t38 * t69 + 0.2e1 * t56 * t47, 0.2e1 * t38 * t70 - 0.2e1 * t56 * t46, -0.2e1 * t10 * t70 + 0.2e1 * t36 * t46 - 0.2e1 * t37 * t47 + 0.2e1 * t9 * t69, 0.2e1 * t36 * t10 - 0.2e1 * t37 * t9 + 0.2e1 * t56 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135 - t142, 0, -t183, t96, t33, t32, 0, 0, 0, 0, -t26, 0, -t27, t96, pkin(3) * t104 * t134 + t108 + t81, t114, (t174 * t26 + t148) * pkin(3) + t169, (t174 * t5 - t102 * t4 + (-t102 * t21 + t174 * t22) * qJD(4)) * pkin(3), 0, 0, -t26, 0, -t27, t96, t92 * t96 + t1 + t81, t114 + t181, pkin(3) * t148 + t92 * t26 + t169, t1 * t92 + (-t102 * t2 + (-t102 * t11 + t174 * t15) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, 0, -t161, 0, -pkin(7) * t160, pkin(7) * t161, 0, 0, 0, 0, -t46, 0, -t47, 0, t29, t28, (t174 * t46 + t147) * pkin(3) + t168, (t174 * t29 - t102 * t28 + (-t102 * t49 + t174 * t50) * qJD(4)) * pkin(3), 0, 0, -t46, 0, -t47, 0, t10, t9, pkin(3) * t147 + t92 * t46 + t168, t10 * t92 + (-t102 * t9 + (-t102 * t36 + t174 * t37) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t90, 0, 0, 0, 0, 0, 0, 0, 0, t89, t90, 0, 0.2e1 * (t152 - t92) * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, -t27, t96, t5, t4, 0, 0, 0, 0, -t26, 0, -t27, t96, 0.2e1 * t95 + t107, t2, t26 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, -t47, 0, t29, t28, 0, 0, 0, 0, -t46, 0, -t47, 0, t10, t9, t46 * pkin(4), t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t128, 0, 0, 0, 0, 0, 0, 0, 0, -t150, -t128, 0, -pkin(4) * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t19;