% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:16
% EndTime: 2019-03-09 06:55:20
% DurationCPUTime: 1.20s
% Computational Cost: add. (2671->170), mult. (5754->275), div. (0->0), fcn. (5641->10), ass. (0->111)
t134 = sin(qJ(4));
t135 = cos(qJ(5));
t103 = t134 * t135;
t136 = cos(qJ(4));
t74 = sin(qJ(5));
t142 = ((t136 * t74 + t103) * qJD(4) + qJD(5) * t103) * pkin(3);
t76 = cos(qJ(6));
t71 = t76 ^ 2;
t73 = sin(qJ(6));
t125 = t73 ^ 2 - t71;
t107 = t125 * qJD(6);
t108 = t134 * qJD(4);
t109 = t136 * qJD(4);
t75 = sin(qJ(3));
t113 = t134 * t75;
t77 = cos(qJ(3));
t141 = -(t136 * qJD(3) + t109) * t77 + qJD(3) * t113 + t75 * t108;
t139 = pkin(3) * t74;
t63 = sin(pkin(11)) * pkin(1) + pkin(7);
t138 = pkin(8) + t63;
t53 = t134 * t77 + t136 * t75;
t50 = t138 * t75;
t51 = t138 * t77;
t92 = t134 * t51 + t136 * t50;
t29 = -t53 * pkin(9) - t92;
t90 = -t136 * t77 + t113;
t91 = t134 * t50 - t136 * t51;
t30 = -t90 * pkin(9) - t91;
t20 = -t135 * t29 + t74 * t30;
t110 = qJD(5) * t135;
t124 = qJD(5) * t74;
t80 = (-qJD(3) - qJD(4)) * t53;
t111 = qJD(3) * t138;
t47 = t75 * t111;
t48 = t77 * t111;
t94 = t134 * t48 + t136 * t47;
t78 = t80 * pkin(9) - t51 * t108 - t50 * t109 - t94;
t93 = t134 * t47 - t136 * t48;
t79 = -pkin(9) * t141 - t50 * t108 + t51 * t109 - t93;
t5 = t30 * t110 + t29 * t124 + t135 * t79 + t74 * t78;
t69 = qJD(6) * t76;
t137 = t20 * t69 + t5 * t73;
t84 = t135 * t90;
t23 = qJD(5) * t84 + t53 * t124 + t135 * t141 - t74 * t80;
t86 = t74 * t90;
t40 = t135 * t53 - t86;
t133 = t40 * t23;
t132 = t40 * t73;
t131 = t76 * t23;
t24 = -qJD(5) * t86 + t53 * t110 - t135 * t80 - t141 * t74;
t130 = t76 * t24;
t39 = t53 * t74 + t84;
t129 = t40 * t130 - t39 * t131;
t67 = t136 * pkin(3) + pkin(4);
t34 = t67 * t124 + t142;
t44 = t134 * t139 - t135 * t67 - pkin(5);
t128 = t34 * t73 + t44 * t69;
t116 = pkin(4) * t124;
t66 = -t135 * pkin(4) - pkin(5);
t127 = t73 * t116 + t66 * t69;
t123 = qJD(6) * t73;
t122 = t75 * qJD(3);
t121 = t77 * qJD(3);
t120 = t73 * t131;
t119 = 0.2e1 * t121;
t118 = pkin(5) * t123;
t117 = pkin(5) * t69;
t68 = pkin(3) * t122;
t115 = t40 * t123;
t114 = t73 * t69;
t64 = -cos(pkin(11)) * pkin(1) - pkin(2);
t42 = t44 * t123;
t112 = -t34 * t76 + t42;
t106 = pkin(3) * t109;
t105 = pkin(3) * t108;
t104 = pkin(4) * t110;
t21 = t135 * t30 + t74 * t29;
t57 = -t77 * pkin(3) + t64;
t41 = t90 * pkin(4) + t57;
t25 = t39 * pkin(5) - t40 * pkin(10) + t41;
t102 = t21 * t76 + t25 * t73;
t101 = t21 * t73 - t25 * t76;
t100 = t23 * t39 - t24 * t40;
t45 = pkin(3) * t103 + t74 * t67 + pkin(10);
t99 = t39 * t45 - t40 * t44;
t65 = pkin(4) * t74 + pkin(10);
t98 = t39 * t65 - t40 * t66;
t54 = t66 * t123;
t97 = -t76 * t116 + t54;
t95 = -t73 * t23 + t40 * t69;
t10 = t115 + t131;
t9 = t39 * t123 - t130;
t87 = t57 * t53;
t33 = -t67 * t110 - t135 * t106 + (t134 * qJD(5) + t108) * t139;
t85 = -t23 * t44 - t24 * t45 + t33 * t39 + t34 * t40;
t81 = -t23 * t66 - t24 * t65 + (-t135 * t39 + t40 * t74) * qJD(5) * pkin(4);
t32 = -t80 * pkin(4) + t68;
t60 = 0.2e1 * t114;
t52 = -0.2e1 * t107;
t38 = t40 ^ 2;
t27 = t91 * qJD(4) + t93;
t26 = t92 * qJD(4) + t94;
t18 = t20 * t123;
t11 = t24 * t73 + t39 * t69;
t8 = -t40 * t107 - t120;
t7 = t24 * pkin(5) + t23 * pkin(10) + t32;
t6 = -0.4e1 * t40 * t114 + t125 * t23;
t4 = -t29 * t110 + t30 * t124 - t135 * t78 + t74 * t79;
t2 = -t102 * qJD(6) + t4 * t73 + t7 * t76;
t1 = t101 * qJD(6) + t4 * t76 - t7 * t73;
t3 = [0, 0, 0, 0, t75 * t119, 0.2e1 * (-t75 ^ 2 + t77 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t64 * t122, t64 * t119, -0.2e1 * t53 * t141, 0.2e1 * t141 * t90 + 0.2e1 * t53 * t80, 0, 0, 0, 0.2e1 * qJD(4) * t87 + 0.2e1 * (t75 * pkin(3) * t90 + t87) * qJD(3), -0.2e1 * t141 * t57 + 0.2e1 * t53 * t68, -0.2e1 * t133, 0.2e1 * t100, 0, 0, 0, 0.2e1 * t24 * t41 + 0.2e1 * t32 * t39, -0.2e1 * t23 * t41 + 0.2e1 * t32 * t40, -0.2e1 * t38 * t114 - 0.2e1 * t71 * t133, 0.2e1 * t38 * t107 + 0.4e1 * t40 * t120, -0.2e1 * t39 * t115 + 0.2e1 * t129, -0.2e1 * t24 * t132 - 0.2e1 * t95 * t39, 0.2e1 * t39 * t24, -0.2e1 * t101 * t24 + 0.2e1 * t5 * t132 + 0.2e1 * t2 * t39 + 0.2e1 * t95 * t20, 0.2e1 * t5 * t76 * t40 + 0.2e1 * t1 * t39 - 0.2e1 * t10 * t20 - 0.2e1 * t102 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t76 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t121, -t122, 0, -t63 * t121, t63 * t122, 0, 0, -t141, t80, 0, t27, t26, 0, 0, -t23, -t24, 0, -t5, t4, t8, t6, t11, -t9, 0, t18 + (-t99 * qJD(6) - t5) * t76 + t85 * t73, t99 * t123 + t85 * t76 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t121, 0, 0, 0, 0, 0, t80, t141, 0, 0, 0, 0, 0, -t24, t23, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t105, -0.2e1 * t106, 0, 0, 0, 0, 0, -0.2e1 * t34, 0.2e1 * t33, t60, t52, 0, 0, 0, 0.2e1 * t112, 0.2e1 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t80, 0, t27, t26, 0, 0, -t23, -t24, 0, -t5, t4, t8, t6, t11, -t9, 0, t18 + (-t98 * qJD(6) - t5) * t76 + t81 * t73, t98 * t123 + t81 * t76 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t141, 0, 0, 0, 0, 0, -t24, t23, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t106, 0, 0, 0, 0, 0 (-pkin(4) - t67) * t124 - t142, -t104 + t33, t60, t52, 0, 0, 0, t42 + t54 + (-t34 - t116) * t76, t127 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t116, -0.2e1 * t104, t60, t52, 0, 0, 0, 0.2e1 * t97, 0.2e1 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, -t5, t4, t8, t6, t11, -t9, 0, t18 + (pkin(5) * t23 - pkin(10) * t24) * t73 + (-t5 + (-pkin(5) * t40 - pkin(10) * t39) * qJD(6)) * t76, t10 * pkin(5) + t9 * pkin(10) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t23, 0, 0, 0, 0, 0, t9, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, t60, t52, 0, 0, 0, t112 - t118, -t117 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t104, t60, t52, 0, 0, 0, t97 - t118, -t117 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t52, 0, 0, 0, -0.2e1 * t118, -0.2e1 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t95, t24, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t123, 0, t33 * t73 - t45 * t69, t45 * t123 + t33 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t123, 0, -t73 * t104 - t65 * t69, -t76 * t104 + t65 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t123, 0, -pkin(10) * t69, pkin(10) * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
