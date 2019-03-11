% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t68 = cos(pkin(11));
t69 = cos(pkin(6));
t54 = t69 * t68;
t67 = sin(pkin(11));
t21 = t43 * t67 - t46 * t54;
t42 = sin(qJ(3));
t70 = qJ(4) * t42;
t45 = cos(qJ(3));
t82 = t21 * t45;
t87 = -pkin(3) * t82 - t21 * t70;
t53 = t69 * t67;
t23 = t43 * t68 + t46 * t53;
t81 = t23 * t45;
t86 = -pkin(3) * t81 - t23 * t70;
t85 = pkin(4) + pkin(8);
t40 = sin(pkin(6));
t84 = g(3) * t40;
t44 = cos(qJ(5));
t36 = pkin(5) * t44 + pkin(4);
t83 = pkin(8) + t36;
t39 = qJ(5) + qJ(6);
t37 = sin(t39);
t80 = t37 * t42;
t38 = cos(t39);
t79 = t38 * t42;
t78 = t40 * t43;
t77 = t40 * t46;
t41 = sin(qJ(5));
t76 = t41 * t42;
t75 = t42 * t44;
t74 = t42 * t46;
t73 = t44 * t46;
t72 = t45 * t46;
t71 = pkin(2) * t77 + pkin(8) * t78;
t18 = t21 * pkin(2);
t66 = -t18 + t87;
t19 = t23 * pkin(2);
t65 = -t19 + t86;
t22 = t43 * t54 + t46 * t67;
t64 = t22 * pkin(8) - t18;
t24 = -t43 * t53 + t46 * t68;
t63 = t24 * pkin(8) - t19;
t62 = pkin(5) * t41 + qJ(4);
t59 = t40 * t68;
t11 = t22 * t45 - t42 * t59;
t10 = t22 * t42 + t45 * t59;
t8 = t10 * pkin(3);
t61 = qJ(4) * t11 - t8;
t58 = t40 * t67;
t13 = t24 * t45 + t42 * t58;
t12 = t24 * t42 - t45 * t58;
t9 = t12 * pkin(3);
t60 = qJ(4) * t13 - t9;
t25 = t42 * t78 - t45 * t69;
t20 = t25 * pkin(3);
t26 = t42 * t69 + t45 * t78;
t57 = qJ(4) * t26 - t20;
t56 = pkin(3) * t40 * t72 + t70 * t77 + t71;
t55 = g(3) * t56;
t47 = -pkin(10) - pkin(9);
t52 = pkin(5) * t76 - t45 * t47;
t4 = g(1) * t12 + g(2) * t10 + g(3) * t25;
t51 = g(1) * t13 + g(2) * t11 + g(3) * t26;
t50 = -g(1) * t23 - g(2) * t21 + g(3) * t77;
t49 = g(1) * t24 + g(2) * t22 + g(3) * t78;
t48 = -g(1) * (t12 * t44 - t23 * t41) - g(2) * (t10 * t44 - t21 * t41) - g(3) * (t25 * t44 + t41 * t77);
t6 = t50 * t45;
t5 = t50 * t42;
t2 = -g(1) * (-t12 * t37 - t23 * t38) - g(2) * (-t10 * t37 - t21 * t38) - g(3) * (-t25 * t37 + t38 * t77);
t1 = -g(1) * (t12 * t38 - t23 * t37) - g(2) * (t10 * t38 - t21 * t37) - g(3) * (t25 * t38 + t37 * t77);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t49, -g(1) * t63 - g(2) * t64 - g(3) * t71, 0, 0, 0, 0, 0, 0, -t49, t6, -t5, -g(1) * (t63 + t86) - g(2) * (t64 + t87) - t55, 0, 0, 0, 0, 0, 0, -g(1) * (-t23 * t76 + t24 * t44) - g(2) * (-t21 * t76 + t22 * t44) - (t41 * t74 + t43 * t44) * t84, -g(1) * (-t23 * t75 - t24 * t41) - g(2) * (-t21 * t75 - t22 * t41) - (-t41 * t43 + t42 * t73) * t84, -t6, -g(1) * (-pkin(9) * t81 + t24 * t85 + t65) - g(2) * (-pkin(9) * t82 + t22 * t85 + t66) - g(3) * ((pkin(4) * t43 + pkin(9) * t72) * t40 + t56) 0, 0, 0, 0, 0, 0, -g(1) * (-t23 * t80 + t24 * t38) - g(2) * (-t21 * t80 + t22 * t38) - (t37 * t74 + t38 * t43) * t84, -g(1) * (-t23 * t79 - t24 * t37) - g(2) * (-t21 * t79 - t22 * t37) - (-t37 * t43 + t38 * t74) * t84, -t6, -g(1) * (-t23 * t52 + t24 * t83 + t65) - g(2) * (-t21 * t52 + t22 * t83 + t66) - t55 - (t36 * t43 + t46 * t52) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t51, -g(1) * t60 - g(2) * t61 - g(3) * t57, 0, 0, 0, 0, 0, 0, -t51 * t41, -t51 * t44, t4, -g(1) * (-pkin(9) * t12 + t60) - g(2) * (-pkin(9) * t10 + t61) - g(3) * (-pkin(9) * t25 + t57) 0, 0, 0, 0, 0, 0, -t51 * t37, -t51 * t38, t4, -g(1) * (t12 * t47 + t13 * t62 - t9) - g(2) * (t10 * t47 + t11 * t62 - t8) - g(3) * (t25 * t47 + t26 * t62 - t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -g(1) * (-t12 * t41 - t23 * t44) - g(2) * (-t10 * t41 - t21 * t44) - g(3) * (-t25 * t41 + t40 * t73) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t48 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
