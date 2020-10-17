% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:05:28
% EndTime: 2019-05-07 11:05:31
% DurationCPUTime: 0.65s
% Computational Cost: add. (559->118), mult. (580->173), div. (0->0), fcn. (590->12), ass. (0->75)
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t30 = g(1) * t56 + g(2) * t53;
t51 = sin(qJ(3));
t54 = cos(qJ(3));
t70 = t56 * t54;
t55 = cos(qJ(2));
t81 = t53 * t55;
t22 = t51 * t81 + t70;
t71 = t56 * t51;
t24 = t53 * t54 - t55 * t71;
t52 = sin(qJ(2));
t85 = g(3) * t52;
t93 = -g(1) * t24 + g(2) * t22 + t51 * t85;
t49 = qJ(3) + pkin(11);
t41 = qJ(5) + t49;
t36 = cos(t41);
t35 = sin(t41);
t75 = t56 * t35;
t11 = t53 * t36 - t55 * t75;
t74 = t56 * t36;
t9 = t35 * t81 + t74;
t3 = -g(1) * t11 + g(2) * t9 + t35 * t85;
t20 = -g(3) * t55 + t30 * t52;
t89 = g(1) * t53;
t83 = t51 * pkin(3);
t82 = t52 * t56;
t80 = t55 * t56;
t39 = sin(t49);
t28 = pkin(4) * t39 + t83;
t18 = pkin(5) * t35 + t28;
t79 = t56 * t18;
t78 = t56 * t28;
t37 = qJ(6) + t41;
t32 = sin(t37);
t77 = t56 * t32;
t33 = cos(t37);
t76 = t56 * t33;
t73 = t56 * t39;
t40 = cos(t49);
t72 = t56 * t40;
t50 = -qJ(4) - pkin(8);
t44 = t54 * pkin(3);
t29 = pkin(4) * t40 + t44;
t69 = t56 * pkin(1) + t53 * pkin(7);
t48 = -pkin(9) + t50;
t19 = pkin(5) * t36 + t29;
t66 = t55 * pkin(2) + t52 * pkin(8);
t64 = -g(2) * t56 + t89;
t13 = pkin(2) + t19;
t42 = -pkin(10) + t48;
t63 = t55 * t13 - t52 * t42;
t27 = pkin(2) + t29;
t61 = t55 * t27 - t52 * t48;
t38 = t44 + pkin(2);
t59 = t55 * t38 - t52 * t50;
t45 = t56 * pkin(7);
t26 = t64 * t52;
t25 = t53 * t51 + t55 * t70;
t23 = -t54 * t81 + t71;
t21 = t30 * t55 + t85;
t17 = t53 * t39 + t55 * t72;
t16 = t53 * t40 - t55 * t73;
t15 = -t40 * t81 + t73;
t14 = t39 * t81 + t72;
t12 = t53 * t35 + t55 * t74;
t10 = -t36 * t81 + t75;
t8 = t53 * t32 + t55 * t76;
t7 = t53 * t33 - t55 * t77;
t6 = -t33 * t81 + t77;
t5 = t32 * t81 + t76;
t4 = g(1) * t12 - g(2) * t10 + t36 * t85;
t2 = g(1) * t8 - g(2) * t6 + t33 * t85;
t1 = -g(1) * t7 + g(2) * t5 + t32 * t85;
t31 = [0, 0, 0, 0, 0, 0, t64, t30, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t55, -t26, -t30, -g(1) * (-t53 * pkin(1) + t45) - g(2) * t69, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t25, -g(1) * t22 - g(2) * t24, t26, -g(1) * t45 - g(2) * (t66 * t56 + t69) - (-pkin(1) - t66) * t89, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t26, -g(1) * (pkin(3) * t71 + t45) - g(2) * (t38 * t80 - t50 * t82 + t69) + (-g(1) * (-pkin(1) - t59) - g(2) * t83) * t53, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t26, -g(1) * (t45 + t78) - g(2) * (t27 * t80 - t48 * t82 + t69) + (-g(1) * (-pkin(1) - t61) - g(2) * t28) * t53, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t26, -g(1) * (t45 + t79) - g(2) * (t13 * t80 - t42 * t82 + t69) + (-g(1) * (-pkin(1) - t63) - g(2) * t18) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t54, -t20 * t51, -t21, -g(3) * t66 + t30 * (pkin(2) * t52 - pkin(8) * t55) 0, 0, 0, 0, 0, 0, t20 * t40, -t20 * t39, -t21, -g(3) * t59 + t30 * (t38 * t52 + t50 * t55) 0, 0, 0, 0, 0, 0, t20 * t36, -t20 * t35, -t21, -g(3) * t61 + t30 * (t27 * t52 + t48 * t55) 0, 0, 0, 0, 0, 0, t20 * t33, -t20 * t32, -t21, -g(3) * t63 + t30 * (t13 * t52 + t42 * t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, g(1) * t25 - g(2) * t23 + t54 * t85, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t14 + t39 * t85, g(1) * t17 - g(2) * t15 + t40 * t85, 0, t93 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (t53 * t29 - t55 * t78) - g(2) * (-t28 * t81 - t56 * t29) + t28 * t85, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t53 * t19 - t55 * t79) - g(2) * (-t18 * t81 - t56 * t19) + t18 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t31;
