% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t47 = cos(qJ(3));
t40 = qJ(1) + pkin(10);
t35 = sin(t40);
t36 = cos(t40);
t21 = g(1) * t36 + g(2) * t35;
t43 = sin(qJ(3));
t92 = t21 * t43;
t97 = g(3) * t47 - t92;
t46 = cos(qJ(4));
t42 = sin(qJ(4));
t75 = t42 * t47;
t16 = -t35 * t46 + t36 * t75;
t73 = t46 * t47;
t17 = t35 * t42 + t36 * t73;
t41 = sin(qJ(6));
t45 = cos(qJ(6));
t3 = t16 * t45 - t17 * t41;
t57 = t41 * t46 - t42 * t45;
t14 = t35 * t75 + t36 * t46;
t15 = t35 * t73 - t36 * t42;
t67 = -t14 * t45 + t15 * t41;
t94 = g(3) * t43;
t96 = -g(1) * t3 + g(2) * t67 + t57 * t94;
t4 = t16 * t41 + t17 * t45;
t56 = t41 * t42 + t45 * t46;
t58 = t14 * t41 + t15 * t45;
t91 = g(1) * t4 + g(2) * t58 + t56 * t94;
t88 = -pkin(4) - pkin(5);
t87 = g(1) * t35;
t37 = t43 * pkin(8);
t38 = t47 * pkin(3);
t80 = t35 * t43;
t79 = t35 * t47;
t78 = t36 * t43;
t77 = t36 * t47;
t76 = t42 * t43;
t74 = t43 * t46;
t72 = t38 + t37;
t71 = qJ(5) * t42;
t48 = cos(qJ(1));
t70 = pkin(1) * t48 + pkin(2) * t36 + pkin(7) * t35;
t69 = -pkin(2) - t38;
t44 = sin(qJ(1));
t68 = -t44 * pkin(1) + pkin(7) * t36;
t66 = -pkin(3) - t71;
t65 = -pkin(4) * t14 + qJ(5) * t15;
t64 = -pkin(4) * t16 + qJ(5) * t17;
t63 = pkin(4) * t73 + t47 * t71 + t72;
t62 = pkin(3) * t77 + pkin(8) * t78 + t70;
t61 = g(1) * t14 - g(2) * t16;
t60 = -g(2) * t36 + t87;
t59 = g(1) * t44 - g(2) * t48;
t54 = -pkin(4) * t15 - t14 * qJ(5) + t68;
t52 = (t69 - t37) * t87;
t2 = g(1) * t16 + g(2) * t14 + g(3) * t76;
t51 = g(1) * t17 + g(2) * t15 + g(3) * t74;
t50 = pkin(4) * t17 + qJ(5) * t16 + t62;
t29 = qJ(5) * t74;
t24 = pkin(8) * t77;
t22 = pkin(8) * t79;
t18 = g(1) * t80 - g(2) * t78;
t8 = t21 * t47 + t94;
t7 = t97 * t46;
t6 = t97 * t42;
t5 = g(1) * t15 - g(2) * t17;
t1 = [0, 0, 0, 0, 0, 0, t59, g(1) * t48 + g(2) * t44, 0, 0, 0, 0, 0, 0, 0, 0, t60, t21, 0, t59 * pkin(1), 0, 0, 0, 0, 0, 0, t60 * t47, -t18, -t21, -g(1) * (-t35 * pkin(2) + t68) - g(2) * t70, 0, 0, 0, 0, 0, 0, t5, -t61, t18, -g(1) * t68 - g(2) * t62 - t52, 0, 0, 0, 0, 0, 0, t5, t18, t61, -g(1) * t54 - g(2) * t50 - t52, 0, 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t4, -g(1) * t67 - g(2) * t3, -t18, -g(1) * (-t15 * pkin(5) + t54) - g(2) * (pkin(5) * t17 - pkin(9) * t78 + t50) - ((-pkin(8) + pkin(9)) * t43 + t69) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t8, -g(1) * (-pkin(3) * t78 + t24) - g(2) * (-pkin(3) * t80 + t22) - g(3) * t72, 0, 0, 0, 0, 0, 0, -t7, -t8, -t6, -g(1) * t24 - g(2) * t22 - g(3) * t63 + (pkin(4) * t46 - t66) * t92, 0, 0, 0, 0, 0, 0, -t97 * t56, t97 * t57, t8, -g(1) * (-pkin(9) * t77 + t24) - g(2) * (-pkin(9) * t79 + t22) - g(3) * (pkin(5) * t73 + t63) + (g(3) * pkin(9) + t21 * (-t46 * t88 - t66)) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t51, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t51, -g(1) * t64 - g(2) * t65 - g(3) * (-pkin(4) * t76 + t29) 0, 0, 0, 0, 0, 0, -t96, -t91, 0, -g(1) * (-pkin(5) * t16 + t64) - g(2) * (-pkin(5) * t14 + t65) - g(3) * (t76 * t88 + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t91, 0, 0;];
taug_reg  = t1;
