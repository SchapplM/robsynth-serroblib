% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:40:28
% EndTime: 2019-05-05 03:40:30
% DurationCPUTime: 0.59s
% Computational Cost: add. (519->116), mult. (939->171), div. (0->0), fcn. (1123->12), ass. (0->66)
t43 = sin(qJ(3));
t46 = cos(qJ(3));
t64 = cos(pkin(6));
t39 = sin(pkin(6));
t44 = sin(qJ(2));
t72 = t39 * t44;
t82 = -t43 * t72 + t46 * t64;
t47 = cos(qJ(2));
t38 = sin(pkin(10));
t60 = t38 * t64;
t63 = cos(pkin(10));
t24 = -t44 * t60 + t47 * t63;
t71 = t39 * t46;
t81 = -t24 * t43 + t38 * t71;
t34 = pkin(3) * t46 + pkin(2);
t70 = t39 * t47;
t25 = t34 * t70;
t80 = g(3) * t25;
t79 = g(3) * t39;
t53 = t64 * t63;
t22 = t38 * t47 + t44 * t53;
t42 = sin(qJ(5));
t78 = t22 * t42;
t77 = t24 * t42;
t37 = qJ(3) + pkin(11);
t36 = cos(t37);
t75 = t36 * t42;
t45 = cos(qJ(5));
t74 = t36 * t45;
t73 = t38 * t39;
t41 = -qJ(4) - pkin(8);
t69 = t41 * t44;
t68 = t42 * t47;
t67 = t45 * t47;
t21 = t38 * t44 - t47 * t53;
t66 = -t21 * t34 - t22 * t41;
t23 = t44 * t63 + t47 * t60;
t65 = -t23 * t34 - t24 * t41;
t59 = t39 * t63;
t57 = t81 * pkin(3);
t35 = sin(t37);
t56 = pkin(4) * t36 + pkin(9) * t35;
t33 = pkin(5) * t45 + pkin(4);
t40 = -qJ(6) - pkin(9);
t55 = t33 * t36 - t35 * t40;
t54 = t82 * pkin(3);
t11 = t22 * t35 + t36 * t59;
t13 = t24 * t35 - t36 * t73;
t17 = t35 * t72 - t36 * t64;
t52 = g(1) * t13 + g(2) * t11 + g(3) * t17;
t12 = t22 * t36 - t35 * t59;
t14 = t24 * t36 + t35 * t73;
t18 = t35 * t64 + t36 * t72;
t51 = g(1) * t14 + g(2) * t12 + g(3) * t18;
t50 = -t22 * t43 - t46 * t59;
t9 = -g(1) * t23 - g(2) * t21 + g(3) * t70;
t49 = g(1) * t24 + g(2) * t22 + g(3) * t72;
t48 = t50 * pkin(3);
t1 = -g(1) * (-t14 * t42 + t23 * t45) - g(2) * (-t12 * t42 + t21 * t45) - g(3) * (-t18 * t42 - t39 * t67);
t8 = t9 * t35;
t6 = t52 * t45;
t5 = t52 * t42;
t4 = -g(1) * (-t23 * t74 + t77) - g(2) * (-t21 * t74 + t78) - (t36 * t67 + t42 * t44) * t79;
t3 = -g(1) * (t23 * t75 + t24 * t45) - g(2) * (t21 * t75 + t22 * t45) - (-t36 * t68 + t44 * t45) * t79;
t2 = -g(1) * (-t14 * t45 - t23 * t42) - g(2) * (-t12 * t45 - t21 * t42) - g(3) * (-t18 * t45 + t39 * t68);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t49, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t46, t9 * t43, -t49, -g(1) * (-pkin(2) * t23 + pkin(8) * t24) - g(2) * (-pkin(2) * t21 + pkin(8) * t22) - (pkin(2) * t47 + pkin(8) * t44) * t79, 0, 0, 0, 0, 0, 0, -t9 * t36, t8, -t49, -g(1) * t65 - g(2) * t66 - g(3) * (-t39 * t69 + t25) 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (-t23 * t56 + t65) - g(2) * (-t21 * t56 + t66) - t80 - (t47 * t56 - t69) * t79, 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (pkin(5) * t77 - t23 * t55 + t65) - g(2) * (pkin(5) * t78 - t21 * t55 + t66) - t80 - (t55 * t47 + (pkin(5) * t42 - t41) * t44) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t50 - g(3) * t82, -g(1) * (-t24 * t46 - t43 * t73) - g(2) * (-t22 * t46 + t43 * t59) - g(3) * (-t43 * t64 - t44 * t71) 0, 0, 0, 0, 0, 0, 0, 0, t52, t51, 0, -g(1) * t57 - g(2) * t48 - g(3) * t54, 0, 0, 0, 0, 0, 0, t6, -t5, -t51, -g(1) * (-pkin(4) * t13 + pkin(9) * t14 + t57) - g(2) * (-t11 * pkin(4) + t12 * pkin(9) + t48) - g(3) * (-pkin(4) * t17 + pkin(9) * t18 + t54) 0, 0, 0, 0, 0, 0, t6, -t5, -t51, -g(1) * (-t13 * t33 - t14 * t40 + t57) - g(2) * (-t11 * t33 - t12 * t40 + t48) - g(3) * (-t17 * t33 - t18 * t40 + t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52;];
taug_reg  = t7;
