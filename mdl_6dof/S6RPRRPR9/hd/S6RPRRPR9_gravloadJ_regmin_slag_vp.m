% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t48 = sin(qJ(1));
t51 = cos(qJ(1));
t70 = sin(pkin(12));
t74 = cos(pkin(6));
t61 = t74 * t70;
t72 = cos(pkin(12));
t29 = t48 * t72 + t51 * t61;
t47 = sin(qJ(3));
t43 = sin(pkin(6));
t71 = sin(pkin(7));
t66 = t43 * t71;
t81 = cos(qJ(3));
t58 = t81 * t66;
t62 = t74 * t72;
t28 = t48 * t70 - t51 * t62;
t73 = cos(pkin(7));
t68 = t28 * t73;
t12 = t29 * t47 + t51 * t58 + t81 * t68;
t13 = t29 * t81 + (-t51 * t66 - t68) * t47;
t67 = t43 * t73;
t21 = t28 * t71 - t51 * t67;
t42 = qJ(4) + pkin(13);
t39 = sin(t42);
t40 = cos(t42);
t4 = t13 * t40 + t21 * t39;
t45 = sin(qJ(6));
t49 = cos(qJ(6));
t90 = -t12 * t49 + t4 * t45;
t89 = t12 * t45 + t4 * t49;
t50 = cos(qJ(4));
t46 = sin(qJ(4));
t80 = t21 * t46;
t88 = t13 * t50 + t80;
t84 = t13 * t46 - t21 * t50;
t59 = t71 * t74;
t60 = t73 * t72;
t19 = t47 * t59 + (t47 * t60 + t81 * t70) * t43;
t27 = -t72 * t66 + t74 * t73;
t30 = -t48 * t61 + t51 * t72;
t54 = t48 * t62 + t51 * t70;
t52 = t54 * t73;
t17 = t30 * t81 + (t48 * t66 - t52) * t47;
t22 = t48 * t67 + t54 * t71;
t8 = -t17 * t46 + t22 * t50;
t83 = -g(1) * t8 + g(2) * t84 - g(3) * (-t19 * t46 + t27 * t50);
t79 = t22 * t46;
t78 = t40 * t45;
t77 = t40 * t49;
t75 = qJ(2) * t43;
t76 = t51 * pkin(1) + t48 * t75;
t69 = -t48 * pkin(1) + t51 * t75;
t16 = t30 * t47 - t48 * t58 + t81 * t52;
t65 = -g(1) * t12 + g(2) * t16;
t64 = g(1) * t51 + g(2) * t48;
t63 = g(1) * t48 - g(2) * t51;
t57 = g(1) * (-t17 * t39 + t22 * t40) + g(2) * (-t13 * t39 + t21 * t40) + g(3) * (-t19 * t39 + t27 * t40);
t18 = -t81 * t59 + (t47 * t70 - t60 * t81) * t43;
t56 = g(1) * t16 + g(2) * t12 + g(3) * t18;
t55 = g(1) * t17 + g(2) * t13 + g(3) * t19;
t44 = -qJ(5) - pkin(10);
t38 = t50 * pkin(4) + pkin(3);
t26 = -g(3) * t74 - t63 * t43;
t11 = t19 * t40 + t27 * t39;
t9 = t17 * t50 + t79;
t7 = t17 * t40 + t22 * t39;
t2 = t16 * t45 + t7 * t49;
t1 = t16 * t49 - t7 * t45;
t3 = [0, t63, t64, g(1) * t29 - g(2) * t30, -g(1) * t28 + g(2) * t54, -t64 * t43, -g(1) * t69 - g(2) * t76, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t17, t65, 0, 0, 0, 0, 0, g(1) * t88 - g(2) * t9, -g(1) * t84 - g(2) * t8, -t65, -g(1) * (-t29 * pkin(2) - pkin(4) * t80 + t12 * t44 - t13 * t38 + t69) - g(2) * (t30 * pkin(2) + pkin(4) * t79 - t16 * t44 + t17 * t38 + t76) + (g(1) * t21 - g(2) * t22) * pkin(9), 0, 0, 0, 0, 0, g(1) * t89 - g(2) * t2, -g(1) * t90 - g(2) * t1; 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t55, 0, 0, 0, 0, 0, t56 * t50, -t56 * t46, -t55, -g(1) * (-t16 * t38 - t17 * t44) - g(2) * (-t12 * t38 - t13 * t44) - g(3) * (-t18 * t38 - t19 * t44) 0, 0, 0, 0, 0, -g(1) * (-t16 * t77 + t17 * t45) - g(2) * (-t12 * t77 + t13 * t45) - g(3) * (-t18 * t77 + t19 * t45) -g(1) * (t16 * t78 + t17 * t49) - g(2) * (t12 * t78 + t13 * t49) - g(3) * (t18 * t78 + t19 * t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, g(1) * t9 + g(2) * t88 - g(3) * (-t19 * t50 - t27 * t46) 0, t83 * pkin(4), 0, 0, 0, 0, 0, -t57 * t49, t57 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t90 - g(3) * (-t11 * t45 + t18 * t49) g(1) * t2 + g(2) * t89 - g(3) * (-t11 * t49 - t18 * t45);];
taug_reg  = t3;
