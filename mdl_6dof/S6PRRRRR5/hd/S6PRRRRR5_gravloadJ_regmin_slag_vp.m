% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t38 = sin(pkin(7));
t44 = sin(qJ(2));
t47 = cos(qJ(2));
t40 = cos(pkin(13));
t65 = cos(pkin(6));
t59 = t40 * t65;
t63 = sin(pkin(13));
t50 = t63 * t44 - t47 * t59;
t64 = cos(pkin(7));
t39 = sin(pkin(6));
t70 = t39 * t40;
t77 = t38 * t70 + t50 * t64;
t54 = t65 * t63;
t51 = t40 * t44 + t47 * t54;
t60 = t39 * t63;
t76 = -t38 * t60 + t51 * t64;
t75 = cos(qJ(3));
t37 = qJ(5) + qJ(6);
t35 = sin(t37);
t46 = cos(qJ(4));
t74 = t35 * t46;
t36 = cos(t37);
t73 = t36 * t46;
t42 = sin(qJ(4));
t72 = t38 * t42;
t71 = t38 * t46;
t69 = t39 * t44;
t68 = t39 * t47;
t41 = sin(qJ(5));
t67 = t41 * t46;
t45 = cos(qJ(5));
t66 = t45 * t46;
t62 = t38 * t69;
t43 = sin(qJ(3));
t58 = t43 * t64;
t57 = t65 * t38;
t55 = t64 * t75;
t29 = t44 * t59 + t63 * t47;
t10 = t29 * t75 - t77 * t43;
t30 = t40 * t47 - t44 * t54;
t12 = t30 * t75 - t76 * t43;
t21 = t43 * t57 + (t75 * t44 + t47 * t58) * t39;
t22 = t50 * t38 - t64 * t70;
t23 = t51 * t38 + t64 * t60;
t28 = -t38 * t68 + t65 * t64;
t53 = g(1) * (-t12 * t42 + t23 * t46) + g(2) * (-t10 * t42 + t22 * t46) + g(3) * (-t21 * t42 + t28 * t46);
t11 = t30 * t43 + t76 * t75;
t20 = t43 * t69 - t55 * t68 - t75 * t57;
t9 = t29 * t43 + t77 * t75;
t52 = g(1) * t11 + g(2) * t9 + g(3) * t20;
t27 = (-t44 * t58 + t75 * t47) * t39;
t26 = (t43 * t47 + t44 * t55) * t39;
t19 = t27 * t46 + t42 * t62;
t18 = -t30 * t58 - t51 * t75;
t17 = t30 * t55 - t51 * t43;
t16 = -t29 * t58 - t50 * t75;
t15 = t29 * t55 - t50 * t43;
t14 = t21 * t46 + t28 * t42;
t8 = t18 * t46 + t30 * t72;
t7 = t16 * t46 + t29 * t72;
t6 = t12 * t46 + t23 * t42;
t4 = t10 * t46 + t22 * t42;
t2 = -g(1) * (-t11 * t35 - t6 * t36) - g(2) * (-t9 * t35 - t4 * t36) - g(3) * (-t14 * t36 - t20 * t35);
t1 = -g(1) * (t11 * t36 - t6 * t35) - g(2) * (-t4 * t35 + t9 * t36) - g(3) * (-t14 * t35 + t20 * t36);
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t51 + g(2) * t50 - g(3) * t68, g(1) * t30 + g(2) * t29 + g(3) * t69, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16 - g(3) * t27, g(1) * t17 + g(2) * t15 + g(3) * t26, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t7 - g(3) * t19, -g(1) * (-t18 * t42 + t30 * t71) - g(2) * (-t16 * t42 + t29 * t71) - g(3) * (-t27 * t42 + t46 * t62) 0, 0, 0, 0, 0, -g(1) * (t17 * t41 + t8 * t45) - g(2) * (t15 * t41 + t7 * t45) - g(3) * (t19 * t45 + t26 * t41) -g(1) * (t17 * t45 - t8 * t41) - g(2) * (t15 * t45 - t7 * t41) - g(3) * (-t19 * t41 + t26 * t45) 0, 0, 0, 0, 0, -g(1) * (t17 * t35 + t8 * t36) - g(2) * (t15 * t35 + t7 * t36) - g(3) * (t19 * t36 + t26 * t35) -g(1) * (t17 * t36 - t8 * t35) - g(2) * (t15 * t36 - t7 * t35) - g(3) * (-t19 * t35 + t26 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, g(1) * t12 + g(2) * t10 + g(3) * t21, 0, 0, 0, 0, 0, t52 * t46, -t52 * t42, 0, 0, 0, 0, 0, -g(1) * (-t11 * t66 + t12 * t41) - g(2) * (t10 * t41 - t9 * t66) - g(3) * (-t20 * t66 + t21 * t41) -g(1) * (t11 * t67 + t12 * t45) - g(2) * (t10 * t45 + t9 * t67) - g(3) * (t20 * t67 + t21 * t45) 0, 0, 0, 0, 0, -g(1) * (-t11 * t73 + t12 * t35) - g(2) * (t10 * t35 - t9 * t73) - g(3) * (-t20 * t73 + t21 * t35) -g(1) * (t11 * t74 + t12 * t36) - g(2) * (t10 * t36 + t9 * t74) - g(3) * (t20 * t74 + t21 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, g(1) * t6 + g(2) * t4 + g(3) * t14, 0, 0, 0, 0, 0, -t53 * t45, t53 * t41, 0, 0, 0, 0, 0, -t53 * t36, t53 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t45 - t6 * t41) - g(2) * (-t4 * t41 + t9 * t45) - g(3) * (-t14 * t41 + t20 * t45) -g(1) * (-t11 * t41 - t6 * t45) - g(2) * (-t4 * t45 - t9 * t41) - g(3) * (-t14 * t45 - t20 * t41) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t3;
