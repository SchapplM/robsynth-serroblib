% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR15_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:25
% EndTime: 2024-09-27 22:26:26
% DurationCPUTime: 1.08s
% Computational Cost: add. (865->153), mult. (1079->262), div. (0->0), fcn. (1151->26), ass. (0->135)
t78 = sin(pkin(5));
t89 = cos(qJ(2));
t131 = t78 * t89;
t84 = sin(qJ(2));
t90 = cos(qJ(1));
t117 = t90 * t84;
t85 = sin(qJ(1));
t122 = t85 * t89;
t80 = cos(pkin(5));
t44 = -t80 * t122 - t117;
t115 = t90 * t89;
t123 = t85 * t84;
t94 = t80 * t115 - t123;
t145 = -g(1) * t44 - g(2) * t94 - g(3) * t131;
t77 = sin(pkin(6));
t141 = pkin(11) * t77;
t82 = sin(qJ(4));
t87 = cos(qJ(4));
t53 = -t82 * pkin(4) + t87 * t141;
t88 = cos(qJ(3));
t46 = t53 * t88;
t52 = -pkin(4) * t87 - t82 * t141;
t49 = pkin(3) - t52;
t83 = sin(qJ(3));
t22 = -t83 * t49 + t46;
t19 = t22 * t89;
t135 = t53 * t83;
t21 = -t49 * t88 - t135;
t20 = pkin(2) - t21;
t10 = -t20 * t84 + t19;
t143 = -t85 / 0.2e1;
t142 = pkin(8) + pkin(9);
t139 = g(3) * t78;
t138 = t84 * pkin(2);
t69 = t89 * pkin(2) + pkin(1);
t136 = t22 * t84;
t134 = t77 * t80;
t133 = t78 * t84;
t132 = t78 * t85;
t130 = t78 * t90;
t79 = cos(pkin(6));
t81 = sin(qJ(5));
t129 = t79 * t81;
t86 = cos(qJ(5));
t128 = t79 * t86;
t127 = t83 * t84;
t76 = qJ(2) + qJ(3);
t74 = qJ(4) + t76;
t66 = sin(t74);
t126 = t85 * t66;
t67 = cos(t74);
t125 = t85 * t67;
t124 = t85 * t81;
t121 = t86 * t85;
t120 = t90 * t66;
t119 = t90 * t67;
t118 = t90 * t81;
t116 = t90 * t86;
t114 = pkin(10) + t142;
t112 = t77 * t132;
t111 = t77 * t130;
t110 = t80 * t118;
t109 = t80 * t124;
t108 = t80 * t121;
t107 = t80 * t116;
t70 = pkin(5) + t76;
t64 = qJ(4) + t70;
t106 = sin(t64) / 0.2e1;
t71 = pkin(5) - t76;
t105 = -qJ(4) + t71;
t104 = g(1) * t90 + g(2) * t85;
t103 = g(1) * t85 - g(2) * t90;
t102 = t20 * t89 + t136;
t23 = t52 * t88 - t135;
t24 = t52 * t83 + t46;
t101 = t23 * t89 - t24 * t84;
t93 = sin(t105);
t40 = t106 - t93 / 0.2e1;
t100 = t90 * t40 + t125;
t99 = t85 * t40 - t119;
t60 = sin(t70);
t61 = sin(t71);
t50 = t60 - t61;
t73 = cos(t76);
t97 = t85 * t73 + t90 * t50 / 0.2e1;
t96 = t50 * t143 + t90 * t73;
t68 = t88 * pkin(3) + pkin(2);
t95 = -pkin(3) * t127 + t68 * t89;
t47 = t104 * t78;
t32 = t79 * t109 - t116;
t38 = t79 * t118 + t108;
t92 = t81 * t112 - t32 * t67 - t38 * t66;
t35 = t79 * t107 - t124;
t37 = t79 * t121 + t110;
t91 = t86 * t111 - t35 * t67 + t37 * t66;
t72 = sin(t76);
t63 = cos(t71);
t62 = cos(t70);
t59 = cos(t105);
t58 = cos(t64) / 0.2e1;
t55 = -pkin(3) * t72 - t138;
t54 = pkin(3) * t73 + t69;
t51 = t63 + t62;
t48 = t80 * t138 - t78 * t142;
t45 = -t80 * t123 + t115;
t43 = -t80 * t117 - t122;
t41 = t59 / 0.2e1 + t58;
t39 = -t79 * t116 + t109;
t36 = t79 * t124 - t107;
t34 = t79 * t110 + t121;
t33 = -t79 * t108 - t118;
t30 = t95 * t80;
t29 = t51 * t143 - t90 * t72;
t28 = t85 * t72 - t90 * t51 / 0.2e1;
t27 = -t78 * t114 + (t89 * t83 * pkin(3) + t84 * t68) * t80;
t26 = -t85 * t41 - t120;
t25 = -t90 * t41 + t126;
t18 = t22 * t133;
t17 = t81 * t111 - t34 * t67 + t36 * t66;
t16 = t86 * t112 + t33 * t67 + t39 * t66;
t15 = t23 * t84 + t24 * t89;
t14 = g(1) * t96 + g(2) * t97 - g(3) * (t62 / 0.2e1 - t63 / 0.2e1);
t13 = -g(1) * t29 + g(2) * t28 - g(3) * (t60 / 0.2e1 + t61 / 0.2e1);
t12 = t21 * t84 + t19;
t11 = (-g(1) * (-t80 * t126 + t119) - g(2) * (t80 * t120 + t125) - t66 * t139) * t77;
t9 = pkin(1) + t102;
t8 = t101 * t80;
t7 = (t21 * t89 - t136) * t80;
t6 = t102 * t80;
t5 = -g(1) * t99 + g(2) * t100 - g(3) * (t58 - t59 / 0.2e1);
t4 = -g(1) * t26 + g(2) * t25 - g(3) * (t106 + t93 / 0.2e1);
t3 = t78 * (t79 * pkin(11) + t114) + t10 * t80;
t2 = -g(1) * (-t33 * t66 + t39 * t67) - g(2) * (-t35 * t66 - t37 * t67) - (-t66 * t128 - t67 * t81) * t139;
t1 = -g(1) * (t32 * t66 - t38 * t67) - g(2) * (-t34 * t66 - t36 * t67) - (-t66 * t129 + t67 * t86) * t139;
t31 = [0, 0, 0, 0, 0, 0, t103, t104, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t43 - g(2) * t45, g(1) * t94 - g(2) * t44, -t47, -g(1) * (-t85 * pkin(1) + pkin(8) * t130) - g(2) * (t90 * pkin(1) + pkin(8) * t132), 0, 0, 0, 0, 0, 0, g(1) * t97 - g(2) * t96, -g(1) * t28 - g(2) * t29, -t47, -g(1) * (-t48 * t90 - t85 * t69) - g(2) * (-t85 * t48 + t90 * t69), 0, 0, 0, 0, 0, 0, g(1) * t100 + g(2) * t99, -g(1) * t25 - g(2) * t26, -t47, -g(1) * (-t90 * t27 - t85 * t54) - g(2) * (-t85 * t27 + t90 * t54), 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t92, -g(1) * t91 - g(2) * t16, -t79 * t47 + (-g(1) * (t80 * t119 - t126) - g(2) * (t80 * t125 + t120)) * t77, -g(1) * (t3 * t90 - t9 * t85) - g(2) * (t3 * t85 + t9 * t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, g(1) * t45 - g(2) * t43 + g(3) * t133, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t145 * pkin(2), 0, 0, 0, 0, 0, 0, t4, t5, 0, -g(1) * (-t85 * t30 + t90 * t55) - g(2) * (t90 * t30 + t85 * t55) - t95 * t139, 0, 0, 0, 0, 0, 0, t1, t2, t11, -g(1) * (t10 * t90 - t6 * t85) - g(2) * (t10 * t85 + t6 * t90) - g(3) * (t20 * t131 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, (t104 * t72 + (t103 * t80 - t139) * (t88 * t89 - t127)) * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, t11, -g(1) * (t12 * t90 + t7 * t85) - g(2) * (t12 * t85 - t7 * t90) - g(3) * (-t21 * t131 + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t11, -g(1) * (t15 * t90 + t8 * t85) - g(2) * (t15 * t85 - t8 * t90) + t101 * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t91 - g(3) * (t86 * t134 + (t67 * t128 - t66 * t81) * t78), g(1) * t92 - g(2) * t17 - g(3) * (-t81 * t134 + (-t67 * t129 - t66 * t86) * t78), 0, 0;];
taug_reg = t31;
