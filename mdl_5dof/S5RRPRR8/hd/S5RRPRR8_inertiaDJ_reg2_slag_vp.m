% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:18
% DurationCPUTime: 1.44s
% Computational Cost: add. (3382->157), mult. (7345->300), div. (0->0), fcn. (7284->8), ass. (0->106)
t72 = sin(qJ(5));
t70 = t72 ^ 2;
t74 = cos(qJ(5));
t71 = t74 ^ 2;
t121 = t70 - t71;
t136 = t121 * qJD(5);
t118 = sin(pkin(9));
t104 = t118 * pkin(2);
t130 = sin(qJ(4));
t131 = cos(qJ(4));
t119 = cos(pkin(9));
t94 = pkin(2) * t119 + pkin(3);
t50 = -t130 * t104 + t131 * t94;
t116 = qJD(5) * t72;
t124 = -qJ(3) - pkin(6);
t73 = sin(qJ(2));
t59 = t124 * t73;
t75 = cos(qJ(2));
t60 = t124 * t75;
t35 = t118 * t60 + t119 * t59;
t53 = t118 * t75 + t119 * t73;
t29 = -pkin(7) * t53 + t35;
t36 = t118 * t59 - t119 * t60;
t83 = t118 * t73 - t119 * t75;
t30 = -pkin(7) * t83 + t36;
t18 = t130 * t29 + t131 * t30;
t102 = qJD(4) * t130;
t105 = t131 * t29;
t100 = qJD(2) * t119;
t99 = qJD(2) * t118;
t122 = t100 * t73 + t75 * t99;
t101 = qJD(2) * t124;
t84 = t75 * qJD(3) + t101 * t73;
t85 = -t73 * qJD(3) + t101 * t75;
t28 = t118 * t85 + t119 * t84;
t23 = -pkin(7) * t122 + t28;
t27 = -t118 * t84 + t119 * t85;
t87 = t100 * t75 - t73 * t99;
t77 = pkin(7) * t87 - t27;
t76 = -qJD(4) * t105 + t102 * t30 + t130 * t77 - t131 * t23;
t81 = t131 * t83;
t33 = t130 * t53 + t81;
t34 = -t130 * t83 + t131 * t53;
t68 = -t75 * pkin(2) - pkin(1);
t41 = pkin(3) * t83 + t68;
t79 = -t33 * pkin(4) + t34 * pkin(8) - t41;
t78 = t74 * t79;
t20 = qJD(4) * t81 + t102 * t53 + t122 * t130 - t131 * t87;
t21 = qJD(4) * t34 + t122 * t131 + t130 * t87;
t115 = t73 * qJD(2);
t108 = pkin(2) * t115;
t37 = pkin(3) * t122 + t108;
t80 = t21 * pkin(4) + t20 * pkin(8) + t37;
t2 = qJD(5) * t78 + t116 * t18 - t72 * t80 + t74 * t76;
t9 = t74 * t18 - t72 * t79;
t3 = -qJD(5) * t9 + t72 * t76 + t74 * t80;
t8 = -t72 * t18 - t78;
t95 = t72 * t8 - t74 * t9;
t135 = qJD(5) * t95 + t2 * t72 - t3 * t74;
t17 = t130 * t30 - t105;
t6 = qJD(4) * t18 + t130 * t23 + t131 * t77;
t134 = t17 * t6;
t133 = t6 * t34;
t5 = t6 * t72;
t69 = qJD(5) * t74;
t132 = t17 * t69 + t5;
t51 = t104 * t131 + t130 * t94;
t43 = t51 * qJD(4);
t129 = t17 * t43;
t128 = t34 * t20;
t127 = t72 * t21;
t126 = t74 * t20;
t125 = t74 * t21;
t47 = -pkin(4) - t50;
t123 = t43 * t72 + t47 * t69;
t120 = t70 + t71;
t114 = t75 * qJD(2);
t113 = 0.2e1 * t33 * t21;
t112 = -0.2e1 * pkin(1) * qJD(2);
t111 = t72 * t126;
t110 = pkin(4) * t116;
t109 = pkin(4) * t69;
t107 = t72 * t69;
t106 = t73 * t114;
t42 = t50 * qJD(4);
t24 = t120 * t42;
t103 = t116 * t47 - t43 * t74;
t32 = t34 ^ 2;
t98 = t32 * t107;
t96 = t72 * t9 + t74 * t8;
t93 = -t33 * t42 + t34 * t43;
t48 = pkin(8) + t51;
t92 = t33 * t48 - t34 * t47;
t90 = -t20 * t72 + t34 * t69;
t89 = t116 * t34 + t126;
t13 = t33 * t69 + t127;
t88 = t116 * t33 - t125;
t82 = -t20 * t47 - t21 * t48 + t93;
t1 = -qJD(5) * t96 - t2 * t74 - t3 * t72;
t63 = -0.2e1 * t107;
t62 = 0.2e1 * t107;
t57 = -0.2e1 * t136;
t15 = t17 * t116;
t11 = t136 * t34 + t111;
t7 = -0.4e1 * t107 * t34 + t121 * t20;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t106, 0.2e1 * (-t73 ^ 2 + t75 ^ 2) * qJD(2), 0, -0.2e1 * t106, 0, 0, t73 * t112, t75 * t112, 0, 0, 0.2e1 * t53 * t87, -0.2e1 * t122 * t53 - 0.2e1 * t83 * t87, 0, 0.2e1 * t83 * t122, 0, 0, 0.2e1 * t108 * t83 + 0.2e1 * t122 * t68, 0.2e1 * t108 * t53 + 0.2e1 * t68 * t87, -0.2e1 * t122 * t36 - 0.2e1 * t27 * t53 - 0.2e1 * t28 * t83 - 0.2e1 * t35 * t87, 0.2e1 * t108 * t68 + 0.2e1 * t27 * t35 + 0.2e1 * t28 * t36, -0.2e1 * t128, 0.2e1 * t20 * t33 - 0.2e1 * t21 * t34, 0, t113, 0, 0, 0.2e1 * t21 * t41 + 0.2e1 * t33 * t37, -0.2e1 * t20 * t41 + 0.2e1 * t34 * t37, -0.2e1 * t17 * t20 - 0.2e1 * t18 * t21 + 0.2e1 * t33 * t76 + 0.2e1 * t133, -0.2e1 * t18 * t76 + 0.2e1 * t41 * t37 + 0.2e1 * t134, -0.2e1 * t128 * t71 - 0.2e1 * t98, 0.4e1 * t111 * t34 + 0.2e1 * t136 * t32, 0.2e1 * t125 * t34 - 0.2e1 * t33 * t89, -0.2e1 * t128 * t70 + 0.2e1 * t98, -0.2e1 * t127 * t34 - 0.2e1 * t33 * t90, t113, 0.2e1 * t17 * t90 + 0.2e1 * t8 * t21 + 0.2e1 * t3 * t33 + 0.2e1 * t34 * t5, 0.2e1 * t133 * t74 - 0.2e1 * t17 * t89 + 0.2e1 * t2 * t33 - 0.2e1 * t9 * t21, 0.2e1 * t135 * t34 + 0.2e1 * t96 * t20, -0.2e1 * t2 * t9 + 0.2e1 * t3 * t8 + 0.2e1 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, -t115, 0, -pkin(6) * t114, pkin(6) * t115, 0, 0, 0, 0, t87, 0, -t122, 0, t27, -t28, (-t118 * t122 - t119 * t87) * pkin(2), (t118 * t28 + t119 * t27) * pkin(2), 0, 0, -t20, 0, -t21, 0, -t6, t76, t20 * t50 - t21 * t51 + t93, t18 * t42 - t6 * t50 - t51 * t76 + t129, -t11, t7, t13, t11, -t88, 0, t15 + (-qJD(5) * t92 - t6) * t74 + t82 * t72, t116 * t92 + t74 * t82 + t132, t1, t1 * t48 - t42 * t95 + t6 * t47 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t43, -0.2e1 * t42, 0, 0.2e1 * t42 * t51 - 0.2e1 * t43 * t50, t62, t57, 0, t63, 0, 0, 0.2e1 * t103, 0.2e1 * t123, 0.2e1 * t24, 0.2e1 * t24 * t48 + 0.2e1 * t47 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t87, 0, t108, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t37, 0, 0, 0, 0, 0, 0, -t88, -t13, t120 * t20, -t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t21, 0, -t6, t76, 0, 0, -t11, t7, t13, t11, -t88, 0, t15 + (pkin(4) * t20 - pkin(8) * t21) * t72 + (-t6 + (-pkin(4) * t34 - pkin(8) * t33) * qJD(5)) * t74, pkin(4) * t89 + pkin(8) * t88 + t132, t1, -t6 * pkin(4) + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, 0, 0, t62, t57, 0, t63, 0, 0, t103 - t110, -t109 + t123, t24, -t43 * pkin(4) + pkin(8) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t57, 0, t63, 0, 0, -0.2e1 * t110, -0.2e1 * t109, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, 0, -t90, t21, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, -t116, 0, -t42 * t72 - t48 * t69, t116 * t48 - t42 * t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t69, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, -t116, 0, -pkin(8) * t69, pkin(8) * t116, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
