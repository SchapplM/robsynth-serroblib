% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = qJ(4) + pkin(10);
t38 = sin(t43);
t39 = cos(t43);
t94 = pkin(5) * t39 + qJ(6) * t38;
t48 = sin(qJ(1));
t51 = cos(qJ(1));
t24 = g(1) * t51 + g(2) * t48;
t44 = qJ(2) + qJ(3);
t40 = sin(t44);
t41 = cos(t44);
t7 = -g(3) * t41 + t24 * t40;
t93 = t94 * t41;
t92 = t41 * pkin(3) + t40 * pkin(9);
t47 = sin(qJ(2));
t91 = pkin(2) * t47;
t87 = g(3) * t40;
t85 = t40 * t48;
t84 = t40 * t51;
t49 = cos(qJ(4));
t36 = t49 * pkin(4) + pkin(3);
t22 = t41 * t36;
t45 = -qJ(5) - pkin(9);
t83 = t41 * t45;
t82 = t41 * t51;
t81 = t48 * t38;
t80 = t48 * t39;
t46 = sin(qJ(4));
t79 = t48 * t46;
t78 = t48 * t49;
t77 = t51 * t38;
t76 = t51 * t39;
t75 = t51 * t46;
t74 = t51 * t49;
t52 = -pkin(8) - pkin(7);
t73 = t51 * t52;
t70 = t41 * t75;
t69 = -t40 * t45 + t22;
t50 = cos(qJ(2));
t42 = t50 * pkin(2);
t37 = t42 + pkin(1);
t25 = t51 * t37;
t68 = -t48 * t52 + t25;
t67 = t42 + t69;
t11 = t41 * t77 - t80;
t9 = t41 * t81 + t76;
t66 = g(1) * t9 - g(2) * t11;
t65 = -pkin(3) * t40 - t91;
t63 = g(1) * t48 - g(2) * t51;
t62 = t36 * t40 + t83;
t14 = t41 * t79 + t74;
t61 = t36 + t94;
t59 = t24 * t41;
t1 = g(1) * t11 + g(2) * t9 + t38 * t87;
t10 = t41 * t80 - t77;
t12 = t41 * t76 + t81;
t57 = g(1) * t12 + g(2) * t10 + t39 * t87;
t56 = pkin(4) * t79 + t36 * t82 - t45 * t84 + t68;
t8 = t59 + t87;
t55 = -g(3) * t50 + t24 * t47;
t54 = -t73 + t45 * t85 + pkin(4) * t75 + (-t37 - t22) * t48;
t31 = pkin(4) * t78;
t27 = pkin(9) * t82;
t26 = t48 * t41 * pkin(9);
t17 = t41 * t74 + t79;
t16 = -t70 + t78;
t15 = -t41 * t78 + t75;
t13 = t63 * t40;
t6 = t7 * t49;
t5 = t7 * t46;
t4 = t7 * t39;
t3 = t7 * t38;
t2 = g(1) * t10 - g(2) * t12;
t18 = [0, 0, 0, 0, 0, 0, t63, t24, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t50, -t63 * t47, -t24, -g(1) * (-t48 * pkin(1) + t51 * pkin(7)) - g(2) * (t51 * pkin(1) + t48 * pkin(7)) 0, 0, 0, 0, 0, 0, t63 * t41, -t13, -t24, -g(1) * (-t48 * t37 - t73) - g(2) * t68, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -g(2) * t25 + (g(1) * t52 - g(2) * t92) * t51 + (-g(1) * (-t37 - t92) + g(2) * t52) * t48, 0, 0, 0, 0, 0, 0, t2, -t66, t13, -g(1) * t54 - g(2) * t56, 0, 0, 0, 0, 0, 0, t2, t13, t66, -g(1) * (-t10 * pkin(5) - t9 * qJ(6) + t54) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, g(3) * t47 + t24 * t50, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t55 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t65 * t51 + t27) - g(2) * (t65 * t48 + t26) - g(3) * (t42 + t92) 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t67 + t24 * (t62 + t91) 0, 0, 0, 0, 0, 0, t4, -t8, t3, -g(3) * (t67 + t93) + t24 * (t61 * t40 + t83 + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-pkin(3) * t84 + t27) - g(2) * (-pkin(3) * t85 + t26) - g(3) * t92, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t69 + t24 * t62, 0, 0, 0, 0, 0, 0, t4, -t8, t3, -g(3) * (t22 + t93) + t45 * t59 + (g(3) * t45 + t24 * t61) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t14 + t46 * t87, g(1) * t17 - g(2) * t15 + t49 * t87, 0, 0, 0, 0, 0, 0, 0, 0, t1, t57, 0, -g(1) * t31 + (g(2) * t74 + t46 * t8) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t57, -g(1) * (-pkin(4) * t70 - t11 * pkin(5) + t12 * qJ(6) + t31) - g(2) * (-t14 * pkin(4) - t9 * pkin(5) + t10 * qJ(6)) - (-pkin(4) * t46 - pkin(5) * t38 + qJ(6) * t39) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t18;
