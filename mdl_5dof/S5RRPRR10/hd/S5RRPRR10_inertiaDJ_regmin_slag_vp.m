% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:43
% EndTime: 2021-01-15 22:02:53
% DurationCPUTime: 1.44s
% Computational Cost: add. (1949->194), mult. (5825->399), div. (0->0), fcn. (5692->10), ass. (0->122)
t81 = sin(qJ(4));
t150 = -0.4e1 * t81;
t138 = pkin(7) + qJ(3);
t82 = sin(qJ(2));
t106 = t138 * t82;
t77 = sin(pkin(5));
t79 = cos(pkin(5));
t85 = cos(qJ(2));
t44 = (pkin(1) * t85 + pkin(2)) * t79 - t77 * t106;
t146 = pkin(1) * t79;
t123 = t82 * t146;
t141 = t77 * t85;
t48 = t138 * t141 + t123;
t76 = sin(pkin(10));
t78 = cos(pkin(10));
t30 = t76 * t44 + t78 * t48;
t21 = pkin(8) * t79 + t30;
t142 = t77 * t82;
t52 = -t78 * t141 + t76 * t142;
t53 = (t76 * t85 + t78 * t82) * t77;
t63 = (-pkin(2) * t85 - pkin(1)) * t77;
t35 = t52 * pkin(3) - t53 * pkin(8) + t63;
t84 = cos(qJ(4));
t149 = t84 * t21 + t81 * t35;
t148 = -t48 * qJD(2) - qJD(3) * t142;
t83 = cos(qJ(5));
t140 = t81 * t83;
t43 = t53 * t84 + t79 * t81;
t133 = qJD(2) * t85;
t110 = t77 * t133;
t134 = qJD(2) * t82;
t111 = t77 * t134;
t51 = t78 * t110 - t76 * t111;
t28 = t43 * qJD(4) + t81 * t51;
t42 = t53 * t81 - t79 * t84;
t128 = qJD(4) * t84;
t115 = t83 * t128;
t80 = sin(qJ(5));
t126 = qJD(5) * t80;
t89 = -t81 * t126 + t115;
t147 = -t28 * t140 - t89 * t42;
t67 = t133 * t146;
t38 = t67 + (-qJD(2) * t106 + qJD(3) * t85) * t77;
t18 = t148 * t76 + t78 * t38;
t50 = qJD(2) * t53;
t66 = pkin(2) * t111;
t32 = pkin(3) * t50 - pkin(8) * t51 + t66;
t8 = -qJD(4) * t149 - t81 * t18 + t84 * t32;
t131 = qJD(4) * t42;
t27 = t84 * t51 - t131;
t33 = t43 * t80 - t83 * t52;
t10 = -qJD(5) * t33 + t83 * t27 + t80 * t50;
t145 = t10 * t80;
t69 = pkin(2) * t76 + pkin(8);
t144 = t69 * t81;
t143 = t69 * t84;
t139 = t83 * t84;
t74 = t83 ^ 2;
t136 = t80 ^ 2 - t74;
t73 = t81 ^ 2;
t135 = -t84 ^ 2 + t73;
t132 = qJD(4) * t33;
t130 = qJD(4) * t81;
t129 = qJD(4) * t83;
t127 = qJD(5) * t73;
t125 = qJD(5) * t83;
t124 = qJD(5) * t84;
t122 = -0.2e1 * pkin(4) * qJD(5);
t121 = t80 * t143;
t119 = t69 * t139;
t70 = -pkin(2) * t78 - pkin(3);
t118 = 0.2e1 * qJD(4) * t70;
t71 = t77 ^ 2;
t117 = t71 * t133;
t116 = t80 * t131;
t114 = t69 * t127;
t113 = t80 * t124;
t112 = t83 * t124;
t108 = t80 * t125;
t107 = t81 * t128;
t34 = t43 * t83 + t52 * t80;
t105 = -t10 * t84 + t34 * t130;
t17 = -t78 * t148 + t38 * t76;
t29 = t44 * t78 - t76 * t48;
t104 = t136 * qJD(5);
t103 = t135 * qJD(4);
t102 = t80 * t115;
t101 = -pkin(4) * t84 - pkin(9) * t81;
t100 = pkin(4) * t81 - pkin(9) * t84;
t13 = pkin(9) * t52 + t149;
t20 = -pkin(3) * t79 - t29;
t15 = pkin(4) * t42 - pkin(9) * t43 + t20;
t6 = t13 * t83 + t15 * t80;
t97 = -t21 * t81 + t35 * t84;
t96 = -t33 * t83 - t34 * t80;
t62 = t101 + t70;
t41 = t62 * t80 + t119;
t12 = -pkin(4) * t52 - t97;
t4 = -t50 * pkin(4) - t8;
t95 = t12 * t125 + t4 * t80;
t94 = t12 * t126 - t4 * t83;
t93 = t100 * t80;
t92 = t52 * t128 + t50 * t81;
t91 = t42 * t125 + t28 * t80;
t90 = t42 * t126 - t28 * t83;
t7 = -t35 * t128 + t21 * t130 - t84 * t18 - t81 * t32;
t58 = t81 * t129 + t113;
t88 = pkin(4) * t28 - pkin(9) * t27 + t17;
t87 = -pkin(9) * t50 + t7;
t60 = t80 * t130 - t112;
t59 = -t81 * t125 - t80 * t128;
t56 = (-pkin(7) * t141 - t123) * qJD(2);
t55 = pkin(7) * t111 - t67;
t40 = t62 * t83 - t121;
t36 = -t52 * t130 + t50 * t84;
t23 = -t41 * qJD(5) + (t83 * t100 + t80 * t144) * qJD(4);
t22 = -qJD(4) * t93 - t62 * t125 + t58 * t69;
t11 = qJD(5) * t34 + t80 * t27 - t83 * t50;
t5 = -t13 * t80 + t15 * t83;
t2 = -t6 * qJD(5) + t80 * t87 + t83 * t88;
t1 = -t15 * t125 + t13 * t126 - t80 * t88 + t83 * t87;
t3 = [0, 0, 0, 0.2e1 * t82 * t117, 0.2e1 * (-t82 ^ 2 + t85 ^ 2) * t71 * qJD(2), 0.2e1 * t79 * t110, -0.2e1 * t79 * t111, 0, -0.2e1 * pkin(1) * t71 * t134 + 0.2e1 * t56 * t79, -0.2e1 * pkin(1) * t117 + 0.2e1 * t55 * t79, -0.2e1 * t17 * t79 + 0.2e1 * t50 * t63 + 0.2e1 * t52 * t66, -0.2e1 * t18 * t79 + 0.2e1 * t51 * t63 + 0.2e1 * t53 * t66, 0.2e1 * t17 * t53 - 0.2e1 * t18 * t52 - 0.2e1 * t29 * t51 - 0.2e1 * t30 * t50, -0.2e1 * t17 * t29 + 0.2e1 * t18 * t30 + 0.2e1 * t63 * t66, 0.2e1 * t43 * t27, -0.2e1 * t27 * t42 - 0.2e1 * t28 * t43, 0.2e1 * t27 * t52 + 0.2e1 * t43 * t50, -0.2e1 * t28 * t52 - 0.2e1 * t42 * t50, 0.2e1 * t52 * t50, 0.2e1 * t17 * t42 + 0.2e1 * t20 * t28 + 0.2e1 * t97 * t50 + 0.2e1 * t8 * t52, -0.2e1 * t149 * t50 + 0.2e1 * t17 * t43 + 0.2e1 * t20 * t27 + 0.2e1 * t7 * t52, 0.2e1 * t34 * t10, -0.2e1 * t10 * t33 - 0.2e1 * t11 * t34, 0.2e1 * t10 * t42 + 0.2e1 * t28 * t34, -0.2e1 * t11 * t42 - 0.2e1 * t28 * t33, 0.2e1 * t42 * t28, 0.2e1 * t11 * t12 + 0.2e1 * t2 * t42 + 0.2e1 * t28 * t5 + 0.2e1 * t33 * t4, 0.2e1 * t1 * t42 + 0.2e1 * t10 * t12 - 0.2e1 * t28 * t6 + 0.2e1 * t34 * t4; 0, 0, 0, 0, 0, t110, -t111, 0, t56, t55, -t17, -t18, (-t50 * t76 - t51 * t78) * pkin(2), (-t17 * t78 + t18 * t76) * pkin(2), t43 * t128 + t27 * t81, t27 * t84 - t81 * t28 + (-t42 * t84 - t43 * t81) * qJD(4), t92, t36, 0, -t50 * t144 - t17 * t84 + t70 * t28 + (-t52 * t143 + t20 * t81) * qJD(4), -t50 * t143 + t17 * t81 + t70 * t27 + (t52 * t144 + t20 * t84) * qJD(4), t10 * t140 + t89 * t34, t96 * t128 + (-t145 - t11 * t83 + (t33 * t80 - t34 * t83) * qJD(5)) * t81, t105 - t147, (t11 - t116) * t84 + (-t91 - t132) * t81, t130 * t42 - t28 * t84, t23 * t42 + t28 * t40 + (-t2 + (t12 * t80 + t33 * t69) * qJD(4)) * t84 + (qJD(4) * t5 + t11 * t69 + t95) * t81, t22 * t42 - t28 * t41 + (-t1 + (t12 * t83 + t34 * t69) * qJD(4)) * t84 + (-qJD(4) * t6 + t10 * t69 - t94) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t107, -0.2e1 * t103, 0, 0, 0, t81 * t118, t84 * t118, 0.2e1 * t74 * t107 - 0.2e1 * t73 * t108, t102 * t150 + 0.2e1 * t136 * t127, 0.2e1 * t81 * t113 + 0.2e1 * t135 * t129, -0.2e1 * t80 * t103 + 0.2e1 * t81 * t112, -0.2e1 * t107, 0.2e1 * t83 * t114 - 0.2e1 * t23 * t84 + 0.2e1 * (t40 + 0.2e1 * t121) * t130, -0.2e1 * t80 * t114 - 0.2e1 * t22 * t84 + 0.2e1 * (-t41 + 0.2e1 * t119) * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, 0, t66, 0, 0, 0, 0, 0, t36, -t92, 0, 0, 0, 0, 0, (-t11 - t116) * t84 + (-t91 + t132) * t81, t105 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t28, t50, t8, t7, t34 * t125 + t145, t96 * qJD(5) + t10 * t83 - t80 * t11, t91, -t90, 0, -pkin(4) * t11 - pkin(9) * t91 + t94, -pkin(4) * t10 + pkin(9) * t90 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t130, 0, -t69 * t128, t69 * t130, -t81 * t104 + t102, t108 * t150 - t136 * t128, t60, t58, 0, (pkin(9) * t139 + (-pkin(4) * t83 + t69 * t80) * t81) * qJD(5) + (t101 * t80 - t119) * qJD(4), (t69 * t140 + t93) * qJD(5) + (t101 * t83 + t121) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t128, 0, 0, 0, 0, 0, -t58, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t108, -0.2e1 * t104, 0, 0, 0, t80 * t122, t83 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t28, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t59, t130, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, -t126, 0, -pkin(9) * t125, pkin(9) * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
