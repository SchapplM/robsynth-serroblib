% Calculate minimal parameter regressor of gravitation load for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t45 = cos(pkin(13));
t47 = cos(pkin(7));
t40 = sin(pkin(14));
t41 = sin(pkin(13));
t44 = cos(pkin(14));
t48 = cos(pkin(6));
t76 = t45 * t48;
t67 = -t41 * t40 + t44 * t76;
t43 = sin(pkin(6));
t72 = sin(pkin(7));
t70 = t43 * t72;
t83 = -t45 * t70 + t67 * t47;
t82 = cos(qJ(4));
t81 = t41 * t48;
t42 = sin(pkin(8));
t50 = sin(qJ(5));
t80 = t42 * t50;
t54 = cos(qJ(5));
t79 = t42 * t54;
t78 = t43 * t47;
t77 = t44 * t47;
t46 = cos(pkin(8));
t51 = sin(qJ(4));
t75 = t46 * t51;
t49 = sin(qJ(6));
t74 = t49 * t54;
t53 = cos(qJ(6));
t73 = t53 * t54;
t71 = t46 * t82;
t69 = t72 * t48;
t66 = -t45 * t40 - t44 * t81;
t38 = -t40 * t81 + t45 * t44;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t59 = t41 * t70 + t66 * t47;
t28 = -t38 * t52 + t59 * t55;
t29 = t38 * t55 + t59 * t52;
t58 = t41 * t78 - t66 * t72;
t56 = t58 * t42;
t10 = t29 * t82 + (t28 * t46 + t56) * t51;
t34 = t55 * t69 + (-t40 * t52 + t55 * t77) * t43;
t35 = t52 * t69 + (t40 * t55 + t52 * t77) * t43;
t62 = -t44 * t70 + t48 * t47;
t61 = t62 * t42;
t19 = t35 * t82 + (t34 * t46 + t61) * t51;
t37 = t40 * t76 + t41 * t44;
t26 = -t37 * t52 + t83 * t55;
t60 = -t45 * t78 - t67 * t72;
t20 = -t26 * t42 + t60 * t46;
t21 = -t28 * t42 + t58 * t46;
t30 = -t34 * t42 + t62 * t46;
t27 = t37 * t55 + t83 * t52;
t57 = t60 * t42;
t8 = t27 * t82 + (t26 * t46 + t57) * t51;
t65 = g(1) * (-t10 * t50 + t21 * t54) + g(2) * (t20 * t54 - t8 * t50) + g(3) * (-t19 * t50 + t30 * t54);
t18 = -t34 * t71 + t35 * t51 - t82 * t61;
t7 = -t26 * t71 + t27 * t51 - t82 * t57;
t9 = -t28 * t71 + t29 * t51 - t82 * t56;
t64 = g(1) * t9 + g(2) * t7 + g(3) * t18;
t23 = t34 * t82 - t35 * t75;
t22 = t34 * t51 + t35 * t71;
t17 = t23 * t54 + t35 * t80;
t16 = t28 * t82 - t29 * t75;
t15 = t28 * t51 + t29 * t71;
t14 = t26 * t82 - t27 * t75;
t13 = t26 * t51 + t27 * t71;
t12 = t19 * t54 + t30 * t50;
t6 = t16 * t54 + t29 * t80;
t5 = t14 * t54 + t27 * t80;
t4 = t10 * t54 + t21 * t50;
t2 = t20 * t50 + t8 * t54;
t1 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(3) * t48 + (-g(1) * t41 + g(2) * t45) * t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -g(1) * t28 - g(2) * t26 - g(3) * t34, g(1) * t29 + g(2) * t27 + g(3) * t35, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14 - g(3) * t23, g(1) * t15 + g(2) * t13 + g(3) * t22, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t5 - g(3) * t17, -g(1) * (-t16 * t50 + t29 * t79) - g(2) * (-t14 * t50 + t27 * t79) - g(3) * (-t23 * t50 + t35 * t79) 0, 0, 0, 0, 0, -g(1) * (t15 * t49 + t6 * t53) - g(2) * (t13 * t49 + t5 * t53) - g(3) * (t17 * t53 + t22 * t49) -g(1) * (t15 * t53 - t6 * t49) - g(2) * (t13 * t53 - t5 * t49) - g(3) * (-t17 * t49 + t22 * t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, g(1) * t10 + g(2) * t8 + g(3) * t19, 0, 0, 0, 0, 0, t64 * t54, -t64 * t50, 0, 0, 0, 0, 0, -g(1) * (t10 * t49 - t9 * t73) - g(2) * (t8 * t49 - t7 * t73) - g(3) * (-t18 * t73 + t19 * t49) -g(1) * (t10 * t53 + t9 * t74) - g(2) * (t8 * t53 + t7 * t74) - g(3) * (t18 * t74 + t19 * t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, g(1) * t4 + g(2) * t2 + g(3) * t12, 0, 0, 0, 0, 0, -t65 * t53, t65 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t49 + t9 * t53) - g(2) * (-t2 * t49 + t7 * t53) - g(3) * (-t12 * t49 + t18 * t53) -g(1) * (-t4 * t53 - t9 * t49) - g(2) * (-t2 * t53 - t7 * t49) - g(3) * (-t12 * t53 - t18 * t49);];
taug_reg  = t1;
