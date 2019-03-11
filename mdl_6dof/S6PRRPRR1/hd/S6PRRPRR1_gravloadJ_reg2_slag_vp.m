% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t48 = sin(pkin(11));
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t70 = cos(pkin(11));
t71 = cos(pkin(6));
t63 = t71 * t70;
t28 = t48 * t56 + t53 * t63;
t66 = t48 * t71;
t30 = -t53 * t66 + t56 * t70;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t49 = sin(pkin(6));
t65 = t49 * t70;
t79 = t49 * t55;
t80 = t49 * t53;
t85 = -g(1) * (-t30 * t52 + t48 * t79) - g(2) * (-t28 * t52 - t55 * t65) - g(3) * (-t52 * t80 + t55 * t71);
t84 = g(3) * t49;
t47 = qJ(3) + pkin(12);
t44 = qJ(5) + t47;
t40 = cos(t44);
t51 = sin(qJ(6));
t83 = t40 * t51;
t54 = cos(qJ(6));
t82 = t40 * t54;
t81 = t48 * t49;
t78 = t49 * t56;
t77 = t51 * t56;
t76 = t54 * t56;
t50 = -qJ(4) - pkin(8);
t27 = t48 * t53 - t56 * t63;
t43 = cos(t47);
t45 = t55 * pkin(3);
t35 = pkin(4) * t43 + t45;
t33 = pkin(2) + t35;
t46 = -pkin(9) + t50;
t75 = -t27 * t33 - t28 * t46;
t29 = t53 * t70 + t56 * t66;
t74 = -t29 * t33 - t30 * t46;
t42 = sin(t47);
t34 = -pkin(3) * t52 - pkin(4) * t42;
t73 = t30 * t34 + t35 * t81;
t72 = t34 * t80 + t35 * t71;
t39 = sin(t44);
t11 = -t28 * t39 - t40 * t65;
t12 = t28 * t40 - t39 * t65;
t69 = pkin(5) * t11 + pkin(10) * t12;
t13 = -t30 * t39 + t40 * t81;
t14 = t30 * t40 + t39 * t81;
t68 = pkin(5) * t13 + pkin(10) * t14;
t20 = -t39 * t80 + t40 * t71;
t21 = t39 * t71 + t40 * t80;
t67 = pkin(5) * t20 + pkin(10) * t21;
t64 = pkin(5) * t40 + pkin(10) * t39;
t62 = t28 * t34 - t35 * t65;
t61 = g(1) * t13 + g(2) * t11 + g(3) * t20;
t5 = g(1) * t14 + g(2) * t12 + g(3) * t21;
t7 = -g(1) * t29 - g(2) * t27 + g(3) * t78;
t59 = g(1) * t30 + g(2) * t28 + g(3) * t80;
t41 = t45 + pkin(2);
t24 = t33 * t78;
t6 = t7 * t39;
t2 = t61 * t54;
t1 = t61 * t51;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t55, t7 * t52, -t59, -g(1) * (-pkin(2) * t29 + pkin(8) * t30) - g(2) * (-pkin(2) * t27 + pkin(8) * t28) - (pkin(2) * t56 + pkin(8) * t53) * t84, 0, 0, 0, 0, 0, 0, -t7 * t43, t7 * t42, -t59, -g(1) * (-t29 * t41 - t30 * t50) - g(2) * (-t27 * t41 - t28 * t50) - (t41 * t56 - t50 * t53) * t84, 0, 0, 0, 0, 0, 0, -t7 * t40, t6, -t59, -g(1) * t74 - g(2) * t75 - g(3) * (-t46 * t80 + t24) 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t82 + t30 * t51) - g(2) * (-t27 * t82 + t28 * t51) - (t40 * t76 + t51 * t53) * t84, -g(1) * (t29 * t83 + t30 * t54) - g(2) * (t27 * t83 + t28 * t54) - (-t40 * t77 + t53 * t54) * t84, -t6, -g(1) * (-t29 * t64 + t74) - g(2) * (-t27 * t64 + t75) - g(3) * t24 - (-t46 * t53 + t56 * t64) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -g(1) * (-t30 * t55 - t52 * t81) - g(2) * (-t28 * t55 + t52 * t65) - g(3) * (-t52 * t71 - t53 * t79) 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t42 + t43 * t81) - g(2) * (-t28 * t42 - t43 * t65) - g(3) * (-t42 * t80 + t43 * t71) -g(1) * (-t30 * t43 - t42 * t81) - g(2) * (-t28 * t43 + t42 * t65) - g(3) * (-t42 * t71 - t43 * t80) 0, t85 * pkin(3), 0, 0, 0, 0, 0, 0, -t61, t5, 0, -g(1) * t73 - g(2) * t62 - g(3) * t72, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t68 + t73) - g(2) * (t62 + t69) - g(3) * (t67 + t72); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t68 - g(2) * t69 - g(3) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t51 + t29 * t54) - g(2) * (-t12 * t51 + t27 * t54) - g(3) * (-t21 * t51 - t49 * t76) -g(1) * (-t14 * t54 - t29 * t51) - g(2) * (-t12 * t54 - t27 * t51) - g(3) * (-t21 * t54 + t49 * t77) 0, 0;];
taug_reg  = t3;
