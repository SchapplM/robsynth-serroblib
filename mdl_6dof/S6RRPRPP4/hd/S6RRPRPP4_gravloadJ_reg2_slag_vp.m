% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t48 = sin(qJ(1));
t51 = cos(qJ(1));
t23 = g(1) * t51 + g(2) * t48;
t47 = sin(qJ(2));
t100 = t23 * t47;
t45 = -qJ(5) - pkin(8);
t50 = cos(qJ(2));
t46 = sin(qJ(4));
t84 = t48 * t46;
t71 = pkin(4) * t84;
t88 = t47 * t48;
t99 = t45 * t88 + t50 * t71;
t78 = t51 * t46;
t33 = pkin(4) * t78;
t87 = t47 * t51;
t98 = t50 * t33 + t45 * t87;
t10 = g(3) * t47 + t23 * t50;
t96 = pkin(4) * t46;
t95 = g(1) * t48;
t91 = g(3) * t50;
t49 = cos(qJ(4));
t90 = t49 * pkin(4);
t40 = t50 * pkin(2);
t89 = t50 * pkin(8);
t44 = qJ(4) + pkin(9);
t36 = sin(t44);
t86 = t48 * t36;
t37 = cos(t44);
t85 = t48 * t37;
t83 = t48 * t49;
t82 = t50 * t45;
t81 = t50 * t51;
t80 = t51 * t36;
t79 = t51 * t37;
t77 = t51 * t49;
t70 = t47 * t83;
t76 = pkin(4) * t70 + t33;
t38 = t47 * qJ(3);
t75 = t40 + t38;
t74 = t51 * pkin(1) + t48 * pkin(7);
t73 = qJ(3) * t50;
t72 = t49 * t91;
t69 = t47 * t78;
t68 = t47 * t77;
t35 = pkin(3) + t90;
t41 = t51 * pkin(7);
t67 = t51 * t35 + t48 * t82 + t41;
t66 = -pkin(1) - t40;
t65 = pkin(2) * t81 + t51 * t38 + t74;
t29 = t48 * t73;
t64 = -pkin(2) * t88 + t29;
t31 = t51 * t73;
t63 = -pkin(2) * t87 + t31;
t62 = pkin(4) * t68 - t71;
t5 = -t47 * t79 + t86;
t7 = t47 * t85 + t80;
t61 = g(1) * t7 + g(2) * t5;
t60 = -g(2) * t51 + t95;
t59 = pkin(5) * t36 - qJ(6) * t37;
t58 = t47 * t96 + t75 - t82;
t57 = t66 - t38;
t6 = t47 * t80 + t85;
t8 = -t47 * t86 + t79;
t56 = -g(1) * t6 + g(2) * t8 + t36 * t91;
t1 = g(1) * t5 - g(2) * t7 + t37 * t91;
t55 = pkin(4) * t69 + t48 * t35 - t45 * t81 + t65;
t52 = ((-qJ(3) - t96) * t47 + t66) * t95;
t16 = t60 * t50;
t15 = t60 * t47;
t14 = -t47 * t84 + t77;
t13 = t70 + t78;
t12 = t69 + t83;
t11 = t68 - t84;
t9 = -t91 + t100;
t4 = t10 * t37;
t3 = t10 * t36;
t2 = -g(1) * t8 - g(2) * t6;
t17 = [0, 0, 0, 0, 0, 0, t60, t23, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, -t23, -g(1) * (-t48 * pkin(1) + t41) - g(2) * t74, 0, 0, 0, 0, 0, 0, -t23, -t16, t15, -g(1) * t41 - g(2) * t65 - t57 * t95, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t12, g(1) * t13 - g(2) * t11, t16, -g(1) * (t51 * pkin(3) + t41) - g(2) * (pkin(8) * t81 + t65) + (-g(1) * (t57 - t89) - g(2) * pkin(3)) * t48, 0, 0, 0, 0, 0, 0, t2, t61, t16, -g(1) * t67 - g(2) * t55 - t52, 0, 0, 0, 0, 0, 0, t2, t16, -t61, -g(1) * (t8 * pkin(5) + t7 * qJ(6) + t67) - g(2) * (t6 * pkin(5) + t5 * qJ(6) + t55) - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -g(1) * t63 - g(2) * t64 - g(3) * t75, 0, 0, 0, 0, 0, 0, -t10 * t46, -t10 * t49, t9, -g(1) * t31 - g(2) * t29 - g(3) * (t75 + t89) + (pkin(2) + pkin(8)) * t100, 0, 0, 0, 0, 0, 0, -t3, -t4, t9, -g(1) * (t63 + t98) - g(2) * (t64 + t99) - g(3) * t58, 0, 0, 0, 0, 0, 0, -t3, t9, t4, -g(1) * (t31 + t98) - g(2) * (t29 + t99) - g(3) * (t47 * t59 + t58) + t23 * (pkin(2) * t47 - t59 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13 + t72, g(1) * t12 - g(2) * t14 - t46 * t91, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t56, 0, pkin(4) * t72 - g(1) * t62 - g(2) * t76, 0, 0, 0, 0, 0, 0, t1, 0, t56, -g(1) * (-t5 * pkin(5) + t6 * qJ(6) + t62) - g(2) * (t7 * pkin(5) - t8 * qJ(6) + t76) - (-pkin(5) * t37 - qJ(6) * t36 - t90) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t17;
