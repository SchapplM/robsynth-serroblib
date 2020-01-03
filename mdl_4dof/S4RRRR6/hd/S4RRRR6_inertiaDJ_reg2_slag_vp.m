% Calculate inertial parameters regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:51
% EndTime: 2019-12-31 17:30:57
% DurationCPUTime: 1.81s
% Computational Cost: add. (1840->227), mult. (5330->469), div. (0->0), fcn. (4761->8), ass. (0->138)
t55 = sin(qJ(3));
t156 = -0.4e1 * t55;
t53 = sin(pkin(4));
t59 = cos(qJ(2));
t144 = t53 * t59;
t134 = qJD(2) * t59;
t112 = t53 * t134;
t56 = sin(qJ(2));
t145 = t53 * t56;
t120 = t55 * t145;
t136 = cos(pkin(4));
t58 = cos(qJ(3));
t20 = -qJD(3) * t120 + (t136 * qJD(3) + t112) * t58;
t155 = -qJD(4) * t144 + t20;
t54 = sin(qJ(4));
t49 = t54 ^ 2;
t57 = cos(qJ(4));
t51 = t57 ^ 2;
t138 = t49 - t51;
t108 = qJD(4) * t138;
t154 = 0.2e1 * t53;
t153 = pkin(7) * t53;
t152 = pkin(7) * t54;
t135 = qJD(2) * t56;
t113 = t53 * t135;
t128 = qJD(4) * t57;
t33 = t136 * t55 + t58 * t145;
t10 = -t57 * t113 + t33 * t128 + t155 * t54;
t151 = t10 * t57;
t129 = qJD(4) * t54;
t11 = t54 * t113 - t33 * t129 + t155 * t57;
t150 = t11 * t54;
t19 = t33 * qJD(3) + t55 * t112;
t149 = t19 * t58;
t148 = t20 * t55;
t139 = t57 * t59;
t21 = t53 * t139 + t33 * t54;
t147 = t21 * t54;
t22 = -t54 * t144 + t33 * t57;
t146 = t22 * t57;
t143 = t54 * t55;
t142 = t55 * t57;
t92 = pkin(2) * t59 + pkin(7) * t56 + pkin(1);
t81 = t92 * t53;
t109 = pkin(1) * t136;
t36 = pkin(6) * t144 + t56 * t109;
t90 = t136 * pkin(7) + t36;
t18 = -t55 * t81 + t58 * t90;
t14 = -pkin(8) * t144 + t18;
t141 = t57 * t14;
t140 = t57 * t58;
t50 = t55 ^ 2;
t137 = -t58 ^ 2 + t50;
t133 = qJD(3) * t55;
t132 = qJD(3) * t57;
t131 = qJD(3) * t58;
t130 = qJD(3) * t59;
t127 = qJD(4) * t58;
t32 = -t136 * t58 + t120;
t126 = 0.2e1 * t32 * t19;
t125 = t58 * t152;
t124 = pkin(6) * t145;
t123 = pkin(7) * t140;
t122 = -0.2e1 * pkin(2) * qJD(3);
t121 = -0.2e1 * pkin(3) * qJD(4);
t119 = pkin(7) * t131;
t48 = t53 ^ 2;
t118 = t48 * t134;
t117 = t57 * t131;
t115 = t54 * t127;
t114 = t57 * t127;
t111 = t54 * t128;
t110 = t55 * t131;
t107 = t137 * qJD(3);
t106 = 0.2e1 * t110;
t105 = qJD(2) * t136;
t104 = t50 * t111;
t103 = t56 * t118;
t102 = t54 * t117;
t101 = t59 * t109;
t100 = pkin(2) * t56 - pkin(7) * t59;
t99 = t19 * pkin(3) - t20 * pkin(8);
t98 = -pkin(3) * t58 - pkin(8) * t55;
t97 = pkin(3) * t55 - pkin(8) * t58;
t69 = -t136 * pkin(2) - t101;
t62 = t32 * pkin(3) - t33 * pkin(8) + t69;
t61 = t62 + t124;
t60 = t57 * t61;
t3 = -t54 * t14 + t60;
t4 = t54 * t61 + t141;
t96 = -t3 * t57 - t4 * t54;
t95 = t56 * t105;
t94 = -t21 * t57 - t22 * t54;
t91 = pkin(2) - t98;
t82 = t57 * t91;
t26 = -t82 - t125;
t27 = -t54 * t91 + t123;
t93 = -t26 * t57 - t27 * t54;
t79 = t55 * t90;
t13 = t79 + (t59 * pkin(3) + t58 * t92) * t53;
t70 = -t101 + t124;
t30 = t70 * qJD(2);
t73 = qJD(3) * t90;
t67 = t55 * t30 - t58 * t73;
t80 = qJD(3) * t92;
t72 = t55 * t80;
t6 = (-t72 + (-t56 * pkin(3) - t58 * t100) * qJD(2)) * t53 - t67;
t89 = t13 * t128 + t54 * t6;
t88 = t13 * t129 - t57 * t6;
t87 = t97 * t54;
t86 = t32 * t133 - t149;
t85 = t32 * t128 + t19 * t54;
t84 = t32 * t129 - t19 * t57;
t83 = qJD(2) * t100;
t78 = t55 * t130 + t58 * t135;
t77 = -t58 * t130 + t55 * t135;
t76 = -t55 * t129 + t117;
t75 = t55 * t132 + t115;
t74 = t55 * t128 + t54 * t131;
t71 = t58 * t80;
t68 = t56 * pkin(8) + t55 * t100;
t31 = t36 * qJD(2);
t66 = -t58 * t30 - t55 * t73;
t1 = t14 * t129 - qJD(4) * t60 - t57 * ((t68 * qJD(2) - t71) * t53 + t66) - t54 * (t31 + t99);
t2 = -t54 * t66 + t57 * (pkin(1) * t95 + t99) + (-t54 * t62 - t141) * qJD(4) + ((-qJD(4) * pkin(6) * t56 + t71) * t54 + (pkin(6) * t139 - t54 * t68) * qJD(2)) * t53;
t65 = t96 * qJD(4) - t1 * t57 - t2 * t54;
t17 = -t58 * t81 - t79;
t7 = (-t55 * t83 + t71) * t53 - t66;
t8 = (t58 * t83 + t72) * t53 + t67;
t64 = -t8 * t55 - t7 * t58 + (-t17 * t58 - t18 * t55) * qJD(3);
t15 = t75 * pkin(7) - qJD(3) * t87 + qJD(4) * t82;
t16 = -t27 * qJD(4) + (pkin(7) * t143 + t57 * t97) * qJD(3);
t63 = t93 * qJD(4) - t15 * t57 - t16 * t54;
t44 = -0.2e1 * t110;
t38 = -0.2e1 * t103;
t29 = t69 + t124;
t23 = t55 * t108 - t102;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t103, 0.2e1 * (-t56 ^ 2 + t59 ^ 2) * t48 * qJD(2), 0.2e1 * t105 * t144, t38, -0.2e1 * t53 * t95, 0, -0.2e1 * t48 * pkin(1) * t135 - 0.2e1 * t31 * t136, -0.2e1 * pkin(1) * t118 + 0.2e1 * t30 * t136, (-t30 * t59 + t31 * t56 + (-t36 * t56 + t59 * t70) * qJD(2)) * t154, -0.2e1 * t30 * t36 + 0.2e1 * t31 * t70, 0.2e1 * t33 * t20, -0.2e1 * t19 * t33 - 0.2e1 * t20 * t32, (t33 * t135 - t20 * t59) * t154, t126, (-t32 * t135 + t19 * t59) * t154, t38, 0.2e1 * t29 * t19 + 0.2e1 * t31 * t32 + 0.2e1 * (t17 * t135 - t59 * t8) * t53, 0.2e1 * t29 * t20 + 0.2e1 * t31 * t33 + 0.2e1 * (-t18 * t135 - t59 * t7) * t53, -0.2e1 * t17 * t20 - 0.2e1 * t18 * t19 + 0.2e1 * t32 * t7 - 0.2e1 * t33 * t8, 0.2e1 * t17 * t8 - 0.2e1 * t18 * t7 + 0.2e1 * t29 * t31, 0.2e1 * t22 * t11, -0.2e1 * t10 * t22 - 0.2e1 * t11 * t21, 0.2e1 * t11 * t32 + 0.2e1 * t19 * t22, 0.2e1 * t21 * t10, -0.2e1 * t10 * t32 - 0.2e1 * t19 * t21, t126, 0.2e1 * t10 * t13 + 0.2e1 * t19 * t3 + 0.2e1 * t2 * t32 + 0.2e1 * t21 * t6, 0.2e1 * t1 * t32 + 0.2e1 * t11 * t13 - 0.2e1 * t19 * t4 + 0.2e1 * t22 * t6, 0.2e1 * t1 * t21 - 0.2e1 * t10 * t4 - 0.2e1 * t11 * t3 - 0.2e1 * t2 * t22, -0.2e1 * t1 * t4 + 0.2e1 * t13 * t6 + 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, 0, -t113, 0, -t31, t30, 0, 0, t33 * t131 + t148, -t55 * t19 + t20 * t58 + (-t32 * t58 - t33 * t55) * qJD(3), t77 * t53, t86, t78 * t53, 0, -pkin(2) * t19 + t29 * t133 - t77 * t153 - t31 * t58, -pkin(2) * t20 + t29 * t131 - t78 * t153 + t31 * t55, (-t149 + t148 + (t32 * t55 + t33 * t58) * qJD(3)) * pkin(7) + t64, -t31 * pkin(2) + t64 * pkin(7), t11 * t142 + t76 * t22, t94 * t131 + (-t151 - t150 + (-t146 + t147) * qJD(4)) * t55, (t32 * t132 - t11) * t58 + (qJD(3) * t22 - t84) * t55, t10 * t143 + t74 * t21, (-qJD(3) * t32 * t54 + t10) * t58 + (-qJD(3) * t21 - t85) * t55, t86, t16 * t32 + t19 * t26 + (-t2 + (pkin(7) * t21 + t13 * t54) * qJD(3)) * t58 + (pkin(7) * t10 + qJD(3) * t3 + t89) * t55, t15 * t32 - t19 * t27 + (-t1 + (pkin(7) * t22 + t13 * t57) * qJD(3)) * t58 + (pkin(7) * t11 - qJD(3) * t4 - t88) * t55, -t10 * t27 - t11 * t26 + t15 * t21 - t16 * t22 + t96 * t131 + (t1 * t54 - t2 * t57 + (t3 * t54 - t4 * t57) * qJD(4)) * t55, -t1 * t27 - t15 * t4 + t16 * t3 + t2 * t26 + (t13 * t131 + t55 * t6) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -0.2e1 * t107, 0, t44, 0, 0, t55 * t122, t58 * t122, 0, 0, 0.2e1 * t51 * t110 - 0.2e1 * t104, t102 * t156 + 0.2e1 * t50 * t108, 0.2e1 * t55 * t115 + 0.2e1 * t137 * t132, 0.2e1 * t49 * t110 + 0.2e1 * t104, -0.2e1 * t54 * t107 + 0.2e1 * t55 * t114, t44, 0.2e1 * t26 * t133 - 0.2e1 * t16 * t58 + 0.2e1 * (t54 * t106 + t50 * t128) * pkin(7), -0.2e1 * t27 * t133 - 0.2e1 * t15 * t58 + 0.2e1 * (t57 * t106 - t50 * t129) * pkin(7), 0.2e1 * t93 * t131 + 0.2e1 * (t15 * t54 - t16 * t57 + (t26 * t54 - t27 * t57) * qJD(4)) * t55, 0.2e1 * pkin(7) ^ 2 * t110 - 0.2e1 * t15 * t27 + 0.2e1 * t16 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t19, t113, t8, t7, 0, 0, t22 * t128 + t150, t94 * qJD(4) - t54 * t10 + t11 * t57, t85, t21 * t129 - t151, -t84, 0, -pkin(3) * t10 - t85 * pkin(8) + t88, -pkin(3) * t11 + t84 * pkin(8) + t89, (-t151 + t150 + (t146 + t147) * qJD(4)) * pkin(8) + t65, -t6 * pkin(3) + t65 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, -t133, 0, -t119, pkin(7) * t133, 0, 0, -t23, t111 * t156 - t138 * t131, t54 * t133 - t114, t23, t75, 0, (pkin(8) * t140 + (-pkin(3) * t57 + t152) * t55) * qJD(4) + (t98 * t54 - t123) * qJD(3), (pkin(7) * t142 + t87) * qJD(4) + (t98 * t57 + t125) * qJD(3), t63, -pkin(3) * t119 + pkin(8) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t111, -0.2e1 * t108, 0, -0.2e1 * t111, 0, 0, t54 * t121, t57 * t121, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t10, t19, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, -t74, t133, t16, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, 0, -t129, 0, -pkin(8) * t128, pkin(8) * t129, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
