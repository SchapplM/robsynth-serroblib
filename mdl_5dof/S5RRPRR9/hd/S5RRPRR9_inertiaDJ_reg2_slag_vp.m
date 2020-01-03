% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:58
% EndTime: 2019-12-31 20:22:05
% DurationCPUTime: 2.10s
% Computational Cost: add. (3432->220), mult. (7583->421), div. (0->0), fcn. (7220->8), ass. (0->125)
t138 = -qJ(3) - pkin(6);
t70 = sin(pkin(9));
t76 = cos(qJ(2));
t144 = t70 * t76;
t71 = cos(pkin(9));
t74 = sin(qJ(2));
t53 = t70 * t74 - t71 * t76;
t59 = t138 * t74;
t55 = t71 * t59;
t81 = (t138 * t144 + t55) * qJD(2) - t53 * qJD(3);
t54 = t71 * t74 + t144;
t67 = -t76 * pkin(2) - pkin(1);
t99 = t53 * pkin(3) + t67;
t87 = -pkin(7) * t54 + t99;
t159 = -qJD(4) * t87 - t81;
t60 = t138 * t76;
t101 = t59 * t70 - t60 * t71;
t131 = t74 * qJD(2);
t125 = pkin(2) * t131;
t49 = t54 * qJD(2);
t130 = t76 * qJD(2);
t50 = t130 * t71 - t131 * t70;
t158 = t49 * pkin(3) - t50 * pkin(7) - qJD(4) * t101 + t125;
t150 = cos(qJ(5));
t115 = t150 * qJD(5);
t157 = t150 * qJD(4) + t115;
t73 = sin(qJ(4));
t134 = qJD(4) * t73;
t123 = t54 * t134;
t75 = cos(qJ(4));
t140 = t75 * t50;
t156 = t123 - t140;
t68 = t73 ^ 2;
t69 = t75 ^ 2;
t135 = t68 - t69;
t113 = qJD(4) * t135;
t155 = qJD(4) + qJD(5);
t94 = t73 * t101;
t20 = t75 * t87 - t94;
t21 = t101 * t75 + t73 * t87;
t8 = -t158 * t73 + t159 * t75;
t9 = t158 * t75 + t159 * t73;
t154 = t8 * t73 - t9 * t75 + (t20 * t73 - t21 * t75) * qJD(4);
t153 = t49 * pkin(4);
t152 = t70 * pkin(2);
t151 = t71 * pkin(2);
t114 = qJD(2) * t138;
t28 = t70 * (t76 * qJD(3) + t114 * t74) - t71 * (-t74 * qJD(3) + t114 * t76);
t40 = -t60 * t70 - t55;
t149 = t40 * t28;
t148 = t54 * t50;
t147 = t54 * t73;
t119 = t150 * t73;
t72 = sin(qJ(5));
t57 = t72 * t75 + t119;
t39 = t155 * t57;
t118 = t150 * t75;
t143 = t72 * t73;
t56 = -t118 + t143;
t146 = t56 * t39;
t38 = t143 * t155 - t157 * t75;
t145 = t57 * t38;
t142 = t73 * t49;
t141 = t73 * t50;
t139 = t75 * t54;
t128 = t54 * t143;
t12 = t50 * t119 - t72 * t123 - qJD(5) * t128 + (t157 * t54 + t50 * t72) * t75;
t30 = t57 * t54;
t137 = -t12 * t57 + t30 * t38;
t133 = qJD(4) * t75;
t132 = qJD(5) * t72;
t36 = 0.2e1 * t53 * t49;
t129 = -0.2e1 * pkin(1) * qJD(2);
t65 = -pkin(3) - t151;
t127 = 0.2e1 * qJD(4) * t65;
t126 = pkin(4) * t134;
t124 = pkin(4) * t132;
t122 = t73 * t133;
t121 = t74 * t130;
t120 = pkin(7) + t152;
t117 = -0.4e1 * t73 * t139;
t112 = pkin(8) + t120;
t52 = t54 ^ 2;
t111 = t52 * t122;
t110 = pkin(4) * t115;
t109 = t73 * t120;
t108 = t75 * t120;
t106 = qJD(4) * t120;
t11 = -t118 * t50 + t141 * t72 + t39 * t54;
t31 = t118 * t54 - t128;
t105 = -t11 * t56 + t31 * t39;
t104 = -t20 * t75 - t21 * t73;
t102 = t38 * t53 - t49 * t57;
t100 = t72 * t112;
t98 = t133 * t53 + t142;
t97 = t133 * t54 + t141;
t14 = -pkin(8) * t147 + t21;
t77 = pkin(8) * t156 + t153 + t9;
t80 = -t94 + t53 * pkin(4) + ((-pkin(8) - pkin(7)) * t54 + t99) * t75;
t78 = t150 * t80;
t82 = -pkin(8) * t97 - t8;
t1 = -qJD(5) * t78 + t132 * t14 - t150 * t82 - t72 * t77;
t95 = t73 * t100;
t93 = qJD(4) * t100;
t91 = t112 * t150;
t90 = t73 * t91;
t89 = -t120 * t49 + t65 * t50;
t88 = t120 * t53 - t65 * t54;
t86 = qJD(4) * t91;
t79 = t72 * t80;
t2 = -qJD(5) * t79 - t115 * t14 + t150 * t77 - t72 * t82;
t58 = -pkin(4) * t75 + t65;
t51 = t112 * t75;
t35 = t150 * t51 - t95;
t34 = -t72 * t51 - t90;
t33 = -t134 * t53 + t49 * t75;
t25 = pkin(4) * t147 + t40;
t22 = t113 * t54 - t140 * t73;
t18 = -t39 * t53 - t49 * t56;
t17 = qJD(5) * t95 - t115 * t51 + t73 * t93 - t75 * t86;
t16 = qJD(5) * t90 + t132 * t51 + t73 * t86 + t75 * t93;
t15 = pkin(4) * t97 + t28;
t7 = t14 * t150 + t79;
t6 = -t72 * t14 + t78;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t121, 0.2e1 * (-t74 ^ 2 + t76 ^ 2) * qJD(2), 0, -0.2e1 * t121, 0, 0, t74 * t129, t76 * t129, 0, 0, 0.2e1 * t148, -0.2e1 * t49 * t54 - 0.2e1 * t50 * t53, 0, t36, 0, 0, 0.2e1 * t125 * t53 + 0.2e1 * t49 * t67, 0.2e1 * t125 * t54 + 0.2e1 * t50 * t67, -0.2e1 * t101 * t49 + 0.2e1 * t28 * t54 + 0.2e1 * t40 * t50 - 0.2e1 * t53 * t81, 0.2e1 * t101 * t81 + 0.2e1 * t125 * t67 + 0.2e1 * t149, 0.2e1 * t148 * t69 - 0.2e1 * t111, 0.2e1 * t113 * t52 + t117 * t50, 0.2e1 * t139 * t49 - 0.2e1 * t156 * t53, 0.2e1 * t148 * t68 + 0.2e1 * t111, -0.2e1 * t142 * t54 - 0.2e1 * t53 * t97, t36, 0.2e1 * t147 * t28 + 0.2e1 * t20 * t49 + 0.2e1 * t40 * t97 + 0.2e1 * t9 * t53, 0.2e1 * t139 * t28 - 0.2e1 * t156 * t40 - 0.2e1 * t21 * t49 + 0.2e1 * t8 * t53, 0.2e1 * t104 * t50 + 0.2e1 * t154 * t54, 0.2e1 * t20 * t9 - 0.2e1 * t21 * t8 + 0.2e1 * t149, -0.2e1 * t31 * t11, 0.2e1 * t11 * t30 - 0.2e1 * t12 * t31, -0.2e1 * t11 * t53 + 0.2e1 * t31 * t49, 0.2e1 * t30 * t12, -0.2e1 * t12 * t53 - 0.2e1 * t30 * t49, t36, 0.2e1 * t12 * t25 + 0.2e1 * t15 * t30 + 0.2e1 * t2 * t53 + 0.2e1 * t49 * t6, 0.2e1 * t1 * t53 - 0.2e1 * t11 * t25 + 0.2e1 * t15 * t31 - 0.2e1 * t49 * t7, 0.2e1 * t1 * t30 + 0.2e1 * t11 * t6 - 0.2e1 * t12 * t7 - 0.2e1 * t2 * t31, -0.2e1 * t1 * t7 + 0.2e1 * t15 * t25 + 0.2e1 * t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, -t131, 0, -pkin(6) * t130, pkin(6) * t131, 0, 0, 0, 0, t50, 0, -t49, 0, -t28, -t81, (-t49 * t70 - t50 * t71) * pkin(2), -t151 * t28 + t152 * t81, -t22, qJD(4) * t117 - t135 * t50, t98, t22, t33, 0, -t28 * t75 + t89 * t73 + (t40 * t73 - t75 * t88) * qJD(4), t28 * t73 + t89 * t75 + (t40 * t75 + t73 * t88) * qJD(4), qJD(4) * t104 - t9 * t73 - t8 * t75, -t8 * t108 - t9 * t109 + t28 * t65 + (-t108 * t20 - t109 * t21) * qJD(4), -t11 * t57 - t31 * t38, -t105 + t137, -t102, t12 * t56 + t30 * t39, t18, 0, t12 * t58 + t126 * t30 + t15 * t56 + t17 * t53 + t25 * t39 + t34 * t49, -t11 * t58 + t126 * t31 + t15 * t57 + t16 * t53 - t25 * t38 - t35 * t49, t1 * t56 + t11 * t34 - t12 * t35 + t16 * t30 - t17 * t31 - t2 * t57 + t38 * t6 - t39 * t7, -t1 * t35 + t126 * t25 + t15 * t58 - t16 * t7 + t17 * t6 + t2 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t122, -0.2e1 * t113, 0, -0.2e1 * t122, 0, 0, t73 * t127, t75 * t127, 0, 0, -0.2e1 * t145, 0.2e1 * t38 * t56 - 0.2e1 * t39 * t57, 0, 0.2e1 * t146, 0, 0, 0.2e1 * t126 * t56 + 0.2e1 * t39 * t58, 0.2e1 * t126 * t57 - 0.2e1 * t38 * t58, 0.2e1 * t16 * t56 - 0.2e1 * t17 * t57 + 0.2e1 * t34 * t38 - 0.2e1 * t35 * t39, 0.2e1 * t126 * t58 - 0.2e1 * t16 * t35 + 0.2e1 * t17 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t50, 0, t125, 0, 0, 0, 0, 0, 0, t33, -t98, (-t68 - t69) * t50, -t154, 0, 0, 0, 0, 0, 0, t18, t102, t105 + t137, -t1 * t57 - t2 * t56 - t38 * t7 - t39 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t57 - t17 * t56 - t34 * t39 - t35 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t145 + 0.2e1 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, 0, -t97, t49, t9, t8, 0, 0, 0, 0, -t11, 0, -t12, t49, -t124 * t53 + t150 * t153 + t2, (-t115 * t53 - t49 * t72) * pkin(4) + t1, (t150 * t11 - t12 * t72 + (-t150 * t30 + t31 * t72) * qJD(5)) * pkin(4), (t150 * t2 - t1 * t72 + (t150 * t7 - t6 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0, -t134, 0, -t75 * t106, t73 * t106, 0, 0, 0, 0, -t38, 0, -t39, 0, t17, t16, (t150 * t38 - t39 * t72 + (-t150 * t56 + t57 * t72) * qJD(5)) * pkin(4), (t150 * t17 - t16 * t72 + (t150 * t35 - t34 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t133, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, (-t150 * t39 - t38 * t72 + (t150 * t57 + t56 * t72) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t124, -0.2e1 * t110, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, t49, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, -t39, 0, t17, t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, -t110, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;