% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:12
% EndTime: 2019-12-05 17:10:17
% DurationCPUTime: 1.05s
% Computational Cost: add. (1160->132), mult. (2905->226), div. (0->0), fcn. (2615->8), ass. (0->100)
t127 = cos(qJ(3));
t128 = cos(qJ(2));
t81 = sin(qJ(3));
t82 = sin(qJ(2));
t55 = -t127 * t128 + t81 * t82;
t57 = t127 * t82 + t81 * t128;
t126 = cos(qJ(5));
t83 = cos(qJ(4));
t103 = t126 * t83;
t79 = sin(qJ(5));
t80 = sin(qJ(4));
t124 = t79 * t80;
t87 = t103 - t124;
t29 = t87 * t57;
t77 = t80 ^ 2;
t78 = t83 ^ 2;
t134 = t77 + t78;
t133 = qJD(2) + qJD(3);
t132 = qJD(4) + qJD(5);
t131 = pkin(8) + pkin(7);
t130 = t83 * pkin(4);
t72 = t81 * pkin(2) + pkin(7);
t129 = -pkin(8) - t72;
t38 = t133 * t57;
t25 = t55 * t38;
t125 = t55 * t81;
t37 = t55 * t133;
t123 = t80 * t37;
t104 = t126 * t80;
t56 = t79 * t83 + t104;
t36 = t132 * t56;
t114 = t80 * qJD(4);
t110 = pkin(4) * t114;
t116 = pkin(2) * qJD(3);
t111 = t81 * t116;
t58 = t110 + t111;
t107 = t127 * pkin(2);
t73 = -t107 - pkin(3);
t61 = t73 - t130;
t121 = t61 * t36 - t58 * t87;
t100 = qJD(5) * t126;
t35 = -qJD(4) * t103 - t83 * t100 + t132 * t124;
t120 = -t61 * t35 + t58 * t56;
t74 = -pkin(3) - t130;
t119 = -t110 * t87 + t74 * t36;
t118 = t56 * t110 - t74 * t35;
t75 = t83 * qJD(4);
t117 = t80 * t111 + t73 * t75;
t115 = qJD(5) * t79;
t113 = pkin(3) * t114;
t112 = pkin(3) * t75;
t109 = pkin(4) * t115;
t108 = t79 * t131;
t106 = t55 * t114;
t105 = t80 * t75;
t14 = t134 * t37;
t76 = t83 * pkin(8);
t51 = t83 * t72 + t76;
t88 = t129 * t104;
t33 = -t79 * t51 + t88;
t97 = t129 * t124;
t34 = t126 * t51 + t97;
t101 = qJD(4) * t129;
t96 = qJD(3) * t107;
t91 = t83 * t96;
t84 = t80 * t101 + t91;
t92 = t80 * t96;
t85 = t83 * t101 - t92;
t7 = -qJD(5) * t88 + t51 * t115 - t126 * t84 - t79 * t85;
t8 = -qJD(5) * t97 - t51 * t100 + t126 * t85 - t79 * t84;
t102 = t33 * t35 - t34 * t36 - t8 * t56 - t7 * t87;
t62 = t83 * pkin(7) + t76;
t94 = t131 * t126;
t90 = t80 * t94;
t17 = t108 * t75 + t62 * t115 + t132 * t90;
t98 = t80 * t108;
t42 = t126 * t62 - t98;
t18 = -t42 * qJD(5) + (-t83 * t94 + t98) * qJD(4);
t41 = -t79 * t62 - t90;
t99 = -t17 * t87 - t18 * t56 + t41 * t35 - t42 * t36;
t95 = pkin(4) * t100;
t89 = -t83 * t111 + t73 * t114;
t86 = t134 * t127;
t65 = -0.2e1 * t105;
t64 = 0.2e1 * t105;
t53 = 0.2e1 * (-t77 + t78) * qJD(4);
t48 = t86 * t116;
t28 = t56 * t57;
t24 = -0.2e1 * t56 * t35;
t23 = -0.2e1 * t87 * t36;
t22 = -t38 * t83 + t106;
t21 = t38 * t80 + t55 * t75;
t11 = -t55 * t35 + t38 * t56;
t10 = t55 * t36 - t38 * t87;
t9 = -0.2e1 * t35 * t87 - 0.2e1 * t56 * t36;
t6 = (t126 * t35 - t36 * t79 + (t126 * t87 + t56 * t79) * qJD(5)) * pkin(4);
t3 = -t132 * t29 + t56 * t37;
t2 = t37 * t103 - t79 * t123 + t36 * t57;
t1 = -t2 * t87 - t28 * t35 - t29 * t36 - t3 * t56;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t57 * t37 + 0.2e1 * t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t57 * t14 + 0.2e1 * t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t29 * t2 - 0.2e1 * t28 * t3 + 0.2e1 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 * qJD(2), -t128 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37, 0, (-t127 * t38 - t37 * t81 + (t127 * t57 + t125) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t22, t21, -t14, t38 * t73 - t72 * t14 + (t86 * t57 + t125) * t116, 0, 0, 0, 0, 0, 0, t10, t11, t1, -t2 * t34 - t28 * t8 - t29 * t7 + t3 * t33 + t38 * t61 + t55 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t111, -0.2e1 * t96, 0, 0, t64, t53, 0, t65, 0, 0, 0.2e1 * t89, 0.2e1 * t117, 0.2e1 * t48, 0.2e1 * (t86 * t72 + t73 * t81) * t116, t24, t9, 0, t23, 0, 0, 0.2e1 * t121, 0.2e1 * t120, 0.2e1 * t102, 0.2e1 * t33 * t8 - 0.2e1 * t34 * t7 + 0.2e1 * t61 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37, 0, 0, 0, 0, 0, 0, 0, 0, t22, t21, -t14, -t38 * pkin(3) - pkin(7) * t14, 0, 0, 0, 0, 0, 0, t10, t11, t1, pkin(4) * t106 - t29 * t17 - t28 * t18 - t2 * t42 + t3 * t41 + t38 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t96, 0, 0, t64, t53, 0, t65, 0, 0, t89 - t113, -t112 + t117, t48, (-pkin(3) * t81 + t86 * pkin(7)) * t116, t24, t9, 0, t23, 0, 0, t119 + t121, t118 + t120, t99 + t102, t61 * t110 - t34 * t17 + t33 * t18 + t8 * t41 - t7 * t42 + t58 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t53, 0, t65, 0, 0, -0.2e1 * t113, -0.2e1 * t112, 0, 0, t24, t9, 0, t23, 0, 0, 0.2e1 * t119, 0.2e1 * t118, 0.2e1 * t99, 0.2e1 * t74 * t110 - 0.2e1 * t42 * t17 + 0.2e1 * t41 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57 * t75 + t123, t57 * t114 + t83 * t37, 0, 0, 0, 0, 0, 0, 0, 0, t3, t2, 0, (t126 * t3 - t2 * t79 + (t126 * t29 + t28 * t79) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, -t114, 0, -t72 * t75 - t92, t72 * t114 - t91, 0, 0, 0, 0, -t35, 0, -t36, 0, t8, t7, t6, (t126 * t8 - t7 * t79 + (t126 * t34 - t33 * t79) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, -t114, 0, -pkin(7) * t75, pkin(7) * t114, 0, 0, 0, 0, -t35, 0, -t36, 0, t18, t17, t6, (t126 * t18 - t17 * t79 + (t126 * t42 - t41 * t79) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t109, -0.2e1 * t95, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, -t36, 0, t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, -t36, 0, t18, t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t95, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
