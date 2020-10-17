% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:35:51
% EndTime: 2019-05-07 08:35:54
% DurationCPUTime: 0.87s
% Computational Cost: add. (379->129), mult. (951->167), div. (0->0), fcn. (1071->8), ass. (0->80)
t49 = sin(qJ(3));
t53 = cos(qJ(3));
t55 = cos(qJ(1));
t82 = t55 * t53;
t51 = sin(qJ(1));
t54 = cos(qJ(2));
t86 = t51 * t54;
t25 = t49 * t86 + t82;
t83 = t55 * t49;
t85 = t53 * t54;
t26 = t51 * t85 - t83;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t108 = t25 * t52 - t26 * t48;
t27 = -t51 * t53 + t54 * t83;
t28 = t51 * t49 + t54 * t82;
t11 = t27 * t52 - t28 * t48;
t110 = g(1) * t11 + g(2) * t108;
t50 = sin(qJ(2));
t91 = t48 * t53;
t65 = -t49 * t52 + t91;
t109 = t65 * t50;
t2 = -g(3) * t109 + t110;
t30 = g(1) * t55 + g(2) * t51;
t106 = t30 * t50;
t16 = -g(3) * t54 + t106;
t97 = g(3) * t50;
t12 = t27 * t48 + t28 * t52;
t92 = t48 * t49;
t64 = t52 * t53 + t92;
t94 = t25 * t48;
t66 = t26 * t52 + t94;
t4 = g(1) * t12 + g(2) * t66 + t64 * t97;
t105 = -pkin(3) - pkin(4);
t102 = g(1) * t51;
t42 = t50 * pkin(8);
t44 = t54 * pkin(2);
t41 = t52 * pkin(5) + pkin(4);
t95 = -pkin(3) - t41;
t90 = t49 * t50;
t89 = t50 * t51;
t88 = t50 * t53;
t87 = t50 * t55;
t84 = t54 * t55;
t81 = t44 + t42;
t80 = t55 * pkin(1) + t51 * pkin(7);
t79 = qJ(4) * t49;
t78 = -pkin(1) - t44;
t77 = -pkin(2) - t79;
t76 = pkin(5) * t48 + qJ(4);
t21 = t25 * pkin(3);
t74 = t26 * qJ(4) - t21;
t23 = t27 * pkin(3);
t73 = t28 * qJ(4) - t23;
t72 = pkin(3) * t85 + t54 * t79 + t81;
t71 = pkin(2) * t84 + pkin(8) * t87 + t80;
t45 = t55 * pkin(7);
t70 = -t26 * pkin(3) - t25 * qJ(4) + t45;
t69 = t28 * pkin(3) + t71;
t68 = g(1) * t25 - g(2) * t27;
t67 = -g(2) * t55 + t102;
t61 = g(3) * t65;
t59 = t27 * qJ(4) + t69;
t58 = (t78 - t42) * t102;
t10 = g(1) * t27 + g(2) * t25 + g(3) * t90;
t57 = g(1) * t28 + g(2) * t26 + g(3) * t88;
t47 = -qJ(6) - pkin(9);
t36 = pkin(8) * t84;
t33 = pkin(8) * t86;
t31 = qJ(4) * t88;
t29 = g(1) * t89 - g(2) * t87;
t17 = t30 * t54 + t97;
t15 = t16 * t53;
t14 = t16 * t49;
t13 = g(1) * t26 - g(2) * t28;
t8 = t16 * t64;
t7 = -t30 * t109 + t54 * t61;
t6 = g(1) * t66 - g(2) * t12;
t5 = g(1) * t108 - g(2) * t11;
t1 = [0, 0, 0, 0, 0, 0, t67, t30, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t54, -t29, -t30, -g(1) * (-t51 * pkin(1) + t45) - g(2) * t80, 0, 0, 0, 0, 0, 0, t13, -t68, t29, -g(1) * t45 - g(2) * t71 - t58, 0, 0, 0, 0, 0, 0, t13, t29, t68, -g(1) * t70 - g(2) * t59 - t58, 0, 0, 0, 0, 0, 0, t6, t5, -t29, -g(1) * (-t26 * pkin(4) + t70) - g(2) * (t28 * pkin(4) - pkin(9) * t87 + t59) - ((-pkin(8) + pkin(9)) * t50 + t78) * t102, 0, 0, 0, 0, 0, 0, t6, t5, -t29, -g(1) * (-pkin(5) * t94 - t26 * t41 + t70) - g(2) * (t76 * t27 + t28 * t41 + t47 * t87 + t69) - ((-pkin(8) - t47) * t50 + t78) * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, -t17, -g(1) * (-pkin(2) * t87 + t36) - g(2) * (-pkin(2) * t89 + t33) - g(3) * t81, 0, 0, 0, 0, 0, 0, t15, -t17, t14, -g(1) * t36 - g(2) * t33 - g(3) * t72 + (pkin(3) * t53 - t77) * t106, 0, 0, 0, 0, 0, 0, t8, t7, t17, -g(1) * (-pkin(9) * t84 + t36) - g(2) * (-pkin(9) * t86 + t33) - g(3) * (pkin(4) * t85 + t72) + (g(3) * pkin(9) + t30 * (-t105 * t53 - t77)) * t50, 0, 0, 0, 0, 0, 0, t8, t7, t17, -g(1) * (t47 * t84 + t36) - g(2) * (t47 * t86 + t33) - g(3) * (t54 * pkin(5) * t92 + t41 * t85 + t72) + (-g(3) * t47 + t30 * (t76 * t49 - t95 * t53 + pkin(2))) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t57, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t57, -g(1) * t73 - g(2) * t74 - g(3) * (-pkin(3) * t90 + t31) 0, 0, 0, 0, 0, 0, t2, -t4, 0, -g(1) * (-t27 * pkin(4) + t73) - g(2) * (-t25 * pkin(4) + t74) - g(3) * (t105 * t90 + t31) 0, 0, 0, 0, 0, 0, t2, -t4, 0, -g(1) * (-t27 * t41 + t76 * t28 - t23) - g(2) * (-t25 * t41 + t76 * t26 - t21) - g(3) * t31 - (pkin(5) * t91 + t95 * t49) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t4, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t4, 0 (t50 * t61 - t110) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16;];
taug_reg  = t1;
