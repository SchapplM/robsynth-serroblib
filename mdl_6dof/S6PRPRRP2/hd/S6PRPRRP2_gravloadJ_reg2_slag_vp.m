% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:38:15
% EndTime: 2019-05-04 23:38:16
% DurationCPUTime: 0.51s
% Computational Cost: add. (644->111), mult. (1687->165), div. (0->0), fcn. (2171->12), ass. (0->76)
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t106 = pkin(4) * t71 + pkin(9) * t68;
t62 = sin(pkin(11));
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t94 = cos(pkin(11));
t84 = -t69 * t62 + t72 * t94;
t63 = sin(pkin(10));
t103 = t63 * t69;
t64 = sin(pkin(6));
t102 = t64 * t68;
t101 = t64 * t71;
t100 = t64 * t72;
t66 = cos(pkin(6));
t99 = t66 * t69;
t98 = t66 * t72;
t67 = sin(qJ(5));
t97 = t67 * t71;
t70 = cos(qJ(5));
t95 = t70 * t71;
t65 = cos(pkin(10));
t93 = t65 * t98;
t56 = -t72 * t62 - t69 * t94;
t54 = t56 * t66;
t37 = -t65 * t54 + t63 * t84;
t19 = -t101 * t65 - t37 * t68;
t20 = -t102 * t65 + t37 * t71;
t92 = t19 * pkin(4) + t20 * pkin(9);
t38 = -t63 * t54 - t65 * t84;
t21 = t101 * t63 + t38 * t68;
t22 = t102 * t63 - t38 * t71;
t91 = t21 * pkin(4) + t22 * pkin(9);
t52 = t56 * t64;
t42 = t52 * t68 + t66 * t71;
t43 = -t52 * t71 + t66 * t68;
t90 = t42 * pkin(4) + t43 * pkin(9);
t51 = t84 * t64;
t88 = pkin(2) * t100 + t51 * pkin(3) - t52 * pkin(8);
t87 = pkin(5) * t70 + qJ(6) * t67;
t86 = -t63 * t98 - t65 * t69;
t85 = t106 * t51 + t88;
t15 = t43 * t67 + t51 * t70;
t53 = t84 * t66;
t36 = t65 * t53 + t63 * t56;
t7 = t20 * t67 + t36 * t70;
t39 = -t63 * t53 + t65 * t56;
t9 = t22 * t67 + t39 * t70;
t1 = g(1) * t9 + g(2) * t7 + g(3) * t15;
t10 = t22 * t70 - t39 * t67;
t16 = t43 * t70 - t51 * t67;
t8 = t20 * t70 - t36 * t67;
t83 = g(1) * t10 + g(2) * t8 + g(3) * t16;
t11 = t36 * t97 - t37 * t70;
t13 = t38 * t70 + t39 * t97;
t23 = t51 * t97 + t52 * t70;
t82 = g(1) * t13 + g(2) * t11 + g(3) * t23;
t81 = g(1) * t21 + g(2) * t19 + g(3) * t42;
t80 = g(1) * t22 + g(2) * t20 + g(3) * t43;
t79 = g(1) * t38 - g(2) * t37 + g(3) * t52;
t78 = g(1) * t39 + g(2) * t36 + g(3) * t51;
t57 = pkin(2) * t93;
t77 = -pkin(2) * t103 + t36 * pkin(3) + pkin(8) * t37 + t57;
t76 = t106 * t36 + t77;
t75 = -g(1) * t86 - g(3) * t100;
t74 = pkin(2) * t86 + t39 * pkin(3) - t38 * pkin(8);
t73 = t106 * t39 + t74;
t50 = -g(3) * t66 + (-g(1) * t63 + g(2) * t65) * t64;
t24 = t51 * t95 - t52 * t67;
t14 = -t38 * t67 + t39 * t95;
t12 = t36 * t95 + t37 * t67;
t6 = t78 * t68;
t4 = t81 * t70;
t3 = t81 * t67;
t2 = -g(1) * t14 - g(2) * t12 - g(3) * t24;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t93 - t103) + t75, -g(1) * (t63 * t99 - t65 * t72) - g(2) * (-t63 * t72 - t65 * t99) + g(3) * t64 * t69, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t79, 0, -g(2) * t57 + (g(2) * t103 + t75) * pkin(2), 0, 0, 0, 0, 0, 0, -t78 * t71, t6, t79, -g(1) * t74 - g(2) * t77 - g(3) * t88, 0, 0, 0, 0, 0, 0, t2, t82, -t6, -g(1) * t73 - g(2) * t76 - g(3) * t85, 0, 0, 0, 0, 0, 0, t2, -t6, -t82, -g(1) * (t14 * pkin(5) + t13 * qJ(6) + t73) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t76) - g(3) * (t24 * pkin(5) + t23 * qJ(6) + t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t80, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t80, -g(1) * t91 - g(2) * t92 - g(3) * t90, 0, 0, 0, 0, 0, 0, -t4, -t80, -t3, -g(1) * (t21 * t87 + t91) - g(2) * (t19 * t87 + t92) - g(3) * (t42 * t87 + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t83, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t83, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (-t15 * pkin(5) + t16 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
