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
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 20:26:39
% EndTime: 2019-12-31 20:26:44
% DurationCPUTime: 1.32s
% Computational Cost: add. (1909->188), mult. (5673->388), div. (0->0), fcn. (5572->10), ass. (0->124)
t80 = sin(qJ(4));
t152 = -0.4e1 * t80;
t138 = pkin(7) + qJ(3);
t81 = sin(qJ(2));
t105 = t138 * t81;
t76 = sin(pkin(5));
t78 = cos(pkin(5));
t84 = cos(qJ(2));
t44 = (pkin(1) * t84 + pkin(2)) * t78 - t76 * t105;
t147 = pkin(1) * t78;
t123 = t81 * t147;
t142 = t76 * t84;
t48 = t138 * t142 + t123;
t75 = sin(pkin(10));
t77 = cos(pkin(10));
t30 = t75 * t44 + t77 * t48;
t21 = t78 * pkin(8) + t30;
t106 = -pkin(2) * t84 - pkin(1);
t143 = t76 * t81;
t52 = -t77 * t142 + t75 * t143;
t53 = (t75 * t84 + t77 * t81) * t76;
t35 = t52 * pkin(3) - t53 * pkin(8) + t106 * t76;
t83 = cos(qJ(4));
t151 = t83 * t21 + t80 * t35;
t134 = qJD(2) * t81;
t70 = t76 ^ 2;
t150 = t134 * t70;
t149 = -t48 * qJD(2) - qJD(3) * t143;
t82 = cos(qJ(5));
t141 = t80 * t82;
t43 = t53 * t83 + t78 * t80;
t133 = qJD(2) * t84;
t111 = t76 * t133;
t112 = t76 * t134;
t51 = t77 * t111 - t75 * t112;
t27 = t43 * qJD(4) + t51 * t80;
t42 = t53 * t80 - t78 * t83;
t128 = qJD(4) * t83;
t107 = t82 * t128;
t79 = sin(qJ(5));
t126 = qJD(5) * t79;
t88 = -t80 * t126 + t107;
t148 = -t27 * t141 - t42 * t88;
t66 = t133 * t147;
t38 = t66 + (-qJD(2) * t105 + qJD(3) * t84) * t76;
t18 = t149 * t75 + t77 * t38;
t50 = qJD(2) * t53;
t65 = pkin(2) * t112;
t32 = t50 * pkin(3) - t51 * pkin(8) + t65;
t8 = -qJD(4) * t151 - t80 * t18 + t83 * t32;
t131 = qJD(4) * t42;
t28 = t51 * t83 - t131;
t33 = t43 * t79 - t52 * t82;
t11 = -t33 * qJD(5) + t28 * t82 + t50 * t79;
t146 = t11 * t79;
t68 = t75 * pkin(2) + pkin(8);
t145 = t68 * t80;
t144 = t68 * t83;
t140 = t82 * t83;
t139 = t83 * t50;
t73 = t82 ^ 2;
t136 = t79 ^ 2 - t73;
t72 = t80 ^ 2;
t135 = -t83 ^ 2 + t72;
t132 = qJD(4) * t33;
t130 = qJD(4) * t80;
t129 = qJD(4) * t82;
t127 = qJD(5) * t72;
t125 = qJD(5) * t82;
t124 = qJD(5) * t83;
t122 = -0.2e1 * pkin(4) * qJD(5);
t121 = t79 * t144;
t119 = t68 * t140;
t69 = -t77 * pkin(2) - pkin(3);
t118 = 0.2e1 * qJD(4) * t69;
t117 = t70 * t133;
t116 = t79 * t131;
t115 = t68 * t127;
t114 = t79 * t124;
t113 = t82 * t124;
t109 = t79 * t125;
t108 = t80 * t128;
t34 = t43 * t82 + t52 * t79;
t104 = -t11 * t83 + t34 * t130;
t17 = -t77 * t149 + t75 * t38;
t29 = t77 * t44 - t75 * t48;
t103 = t136 * qJD(5);
t102 = t135 * qJD(4);
t101 = t79 * t107;
t100 = -t83 * pkin(4) - t80 * pkin(9);
t99 = pkin(4) * t80 - pkin(9) * t83;
t13 = t52 * pkin(9) + t151;
t20 = -t78 * pkin(3) - t29;
t15 = t42 * pkin(4) - t43 * pkin(9) + t20;
t6 = t13 * t82 + t15 * t79;
t96 = -t80 * t21 + t83 * t35;
t95 = -t33 * t82 - t34 * t79;
t62 = t100 + t69;
t41 = t79 * t62 + t119;
t12 = -t52 * pkin(4) - t96;
t4 = -t50 * pkin(4) - t8;
t94 = t12 * t125 + t4 * t79;
t93 = t12 * t126 - t4 * t82;
t92 = t99 * t79;
t91 = t52 * t128 + t50 * t80;
t90 = t42 * t125 + t27 * t79;
t89 = t42 * t126 - t27 * t82;
t7 = -t35 * t128 + t21 * t130 - t83 * t18 - t80 * t32;
t58 = t80 * t129 + t114;
t87 = pkin(4) * t27 - pkin(9) * t28 + t17;
t86 = -pkin(9) * t50 + t7;
t60 = t79 * t130 - t113;
t59 = -t80 * t125 - t79 * t128;
t56 = (-pkin(7) * t142 - t123) * qJD(2);
t55 = pkin(7) * t112 - t66;
t40 = t82 * t62 - t121;
t36 = -t52 * t130 + t139;
t23 = -t41 * qJD(5) + (t79 * t145 + t82 * t99) * qJD(4);
t22 = -qJD(4) * t92 - t62 * t125 + t58 * t68;
t10 = t34 * qJD(5) + t28 * t79 - t50 * t82;
t5 = -t13 * t79 + t15 * t82;
t2 = -t6 * qJD(5) + t79 * t86 + t82 * t87;
t1 = -t15 * t125 + t13 * t126 - t79 * t87 + t82 * t86;
t3 = [0, 0, 0, 0.2e1 * t81 * t117, 0.2e1 * (-t81 ^ 2 + t84 ^ 2) * t70 * qJD(2), 0.2e1 * t78 * t111, -0.2e1 * t78 * t112, 0, -0.2e1 * pkin(1) * t150 + 0.2e1 * t56 * t78, -0.2e1 * pkin(1) * t117 + 0.2e1 * t55 * t78, 0.2e1 * t17 * t53 - 0.2e1 * t18 * t52 - 0.2e1 * t29 * t51 - 0.2e1 * t30 * t50, 0.2e1 * t106 * pkin(2) * t150 - 0.2e1 * t29 * t17 + 0.2e1 * t30 * t18, 0.2e1 * t43 * t28, -0.2e1 * t27 * t43 - 0.2e1 * t28 * t42, 0.2e1 * t28 * t52 + 0.2e1 * t43 * t50, -0.2e1 * t27 * t52 - 0.2e1 * t42 * t50, 0.2e1 * t52 * t50, 0.2e1 * t17 * t42 + 0.2e1 * t20 * t27 + 0.2e1 * t50 * t96 + 0.2e1 * t8 * t52, -0.2e1 * t151 * t50 + 0.2e1 * t17 * t43 + 0.2e1 * t20 * t28 + 0.2e1 * t7 * t52, 0.2e1 * t34 * t11, -0.2e1 * t10 * t34 - 0.2e1 * t11 * t33, 0.2e1 * t11 * t42 + 0.2e1 * t27 * t34, -0.2e1 * t10 * t42 - 0.2e1 * t27 * t33, 0.2e1 * t42 * t27, 0.2e1 * t10 * t12 + 0.2e1 * t2 * t42 + 0.2e1 * t27 * t5 + 0.2e1 * t33 * t4, 0.2e1 * t1 * t42 + 0.2e1 * t11 * t12 - 0.2e1 * t27 * t6 + 0.2e1 * t34 * t4; 0, 0, 0, 0, 0, t111, -t112, 0, t56, t55, (-t50 * t75 - t51 * t77) * pkin(2), (-t17 * t77 + t18 * t75) * pkin(2), t43 * t128 + t28 * t80, -t80 * t27 + t28 * t83 + (-t42 * t83 - t43 * t80) * qJD(4), t91, t36, 0, -t50 * t145 - t17 * t83 + t69 * t27 + (-t144 * t52 + t20 * t80) * qJD(4), -t68 * t139 + t17 * t80 + t69 * t28 + (t145 * t52 + t20 * t83) * qJD(4), t11 * t141 + t34 * t88, t95 * t128 + (-t10 * t82 - t146 + (t33 * t79 - t34 * t82) * qJD(5)) * t80, t104 - t148, (t10 - t116) * t83 + (-t90 - t132) * t80, t130 * t42 - t27 * t83, t23 * t42 + t40 * t27 + (-t2 + (t12 * t79 + t33 * t68) * qJD(4)) * t83 + (qJD(4) * t5 + t10 * t68 + t94) * t80, t22 * t42 - t41 * t27 + (-t1 + (t12 * t82 + t34 * t68) * qJD(4)) * t83 + (-qJD(4) * t6 + t11 * t68 - t93) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t108, -0.2e1 * t102, 0, 0, 0, t80 * t118, t83 * t118, 0.2e1 * t108 * t73 - 0.2e1 * t109 * t72, t101 * t152 + 0.2e1 * t127 * t136, 0.2e1 * t114 * t80 + 0.2e1 * t129 * t135, -0.2e1 * t102 * t79 + 0.2e1 * t113 * t80, -0.2e1 * t108, 0.2e1 * t82 * t115 - 0.2e1 * t23 * t83 + 0.2e1 * (t40 + 0.2e1 * t121) * t130, -0.2e1 * t79 * t115 - 0.2e1 * t22 * t83 + 0.2e1 * (-t41 + 0.2e1 * t119) * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, t36, -t91, 0, 0, 0, 0, 0, (-t10 - t116) * t83 + (-t90 + t132) * t80, t104 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, t50, t8, t7, t125 * t34 + t146, qJD(5) * t95 - t79 * t10 + t11 * t82, t90, -t89, 0, -pkin(4) * t10 - pkin(9) * t90 + t93, -pkin(4) * t11 + pkin(9) * t89 + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t130, 0, -t68 * t128, t68 * t130, -t103 * t80 + t101, t109 * t152 - t128 * t136, t60, t58, 0, (pkin(9) * t140 + (-pkin(4) * t82 + t68 * t79) * t80) * qJD(5) + (t100 * t79 - t119) * qJD(4), (t141 * t68 + t92) * qJD(5) + (t100 * t82 + t121) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, -t128, 0, 0, 0, 0, 0, -t58, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t109, -0.2e1 * t103, 0, 0, 0, t79 * t122, t82 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, t27, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t59, t130, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, -t126, 0, -pkin(9) * t125, pkin(9) * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
