% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t44 = cos(pkin(12));
t71 = cos(pkin(6));
t63 = t44 * t71;
t69 = sin(pkin(12));
t30 = t49 * t63 + t69 * t52;
t58 = t71 * t69;
t32 = t44 * t52 - t49 * t58;
t43 = sin(pkin(6));
t73 = t43 * t49;
t83 = g(1) * t32 + g(2) * t30 + g(3) * t73;
t80 = cos(qJ(3));
t41 = qJ(4) + pkin(13);
t39 = sin(t41);
t42 = sin(pkin(7));
t79 = t39 * t42;
t40 = cos(t41);
t46 = sin(qJ(6));
t78 = t40 * t46;
t50 = cos(qJ(6));
t77 = t40 * t50;
t47 = sin(qJ(4));
t76 = t42 * t47;
t51 = cos(qJ(4));
t75 = t42 * t51;
t74 = t43 * t44;
t72 = t43 * t52;
t70 = cos(pkin(7));
t67 = t42 * t74;
t66 = t42 * t73;
t64 = t43 * t69;
t48 = sin(qJ(3));
t62 = t48 * t70;
t61 = t71 * t42;
t60 = t42 * t64;
t59 = t70 * t80;
t29 = -t69 * t49 + t52 * t63;
t10 = t29 * t62 + t30 * t80 - t48 * t67;
t31 = -t44 * t49 - t52 * t58;
t12 = t32 * t80 + (t70 * t31 + t60) * t48;
t19 = t48 * t61 + (t80 * t49 + t52 * t62) * t43;
t20 = -t29 * t42 - t70 * t74;
t21 = -t31 * t42 + t70 * t64;
t28 = -t42 * t72 + t71 * t70;
t57 = g(1) * (-t12 * t39 + t21 * t40) + g(2) * (-t10 * t39 + t20 * t40) + g(3) * (-t19 * t39 + t28 * t40);
t11 = -t31 * t59 + t32 * t48 - t80 * t60;
t18 = t48 * t73 - t59 * t72 - t80 * t61;
t9 = -t29 * t59 + t30 * t48 + t80 * t67;
t56 = g(1) * t11 + g(2) * t9 + g(3) * t18;
t55 = g(1) * t12 + g(2) * t10 + g(3) * t19;
t13 = t29 * t48 + t30 * t59;
t15 = t31 * t48 + t32 * t59;
t26 = (t48 * t52 + t49 * t59) * t43;
t54 = g(1) * t15 + g(2) * t13 + g(3) * t26;
t53 = -g(1) * (-t12 * t47 + t21 * t51) - g(2) * (-t10 * t47 + t20 * t51) - g(3) * (-t19 * t47 + t28 * t51);
t45 = -qJ(5) - pkin(10);
t38 = t51 * pkin(4) + pkin(3);
t27 = (-t49 * t62 + t80 * t52) * t43;
t17 = t27 * t40 + t39 * t66;
t16 = t31 * t80 - t32 * t62;
t14 = t29 * t80 - t30 * t62;
t8 = t19 * t40 + t28 * t39;
t6 = t16 * t40 + t32 * t79;
t5 = t14 * t40 + t30 * t79;
t4 = t12 * t40 + t21 * t39;
t2 = t10 * t40 + t20 * t39;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(1) * t31 - g(2) * t29 - g(3) * t72, t83, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14 - g(3) * t27, t54, 0, 0, 0, 0, 0, -g(1) * (t16 * t51 + t32 * t76) - g(2) * (t14 * t51 + t30 * t76) - g(3) * (t27 * t51 + t47 * t66) -g(1) * (-t16 * t47 + t32 * t75) - g(2) * (-t14 * t47 + t30 * t75) - g(3) * (-t27 * t47 + t51 * t66) -t54, -g(1) * (t31 * pkin(2) - t15 * t45 + t16 * t38) - g(2) * (t29 * pkin(2) - t13 * t45 + t14 * t38) - g(3) * (pkin(2) * t72 - t26 * t45 + t27 * t38) - t83 * t42 * (pkin(4) * t47 + pkin(9)) 0, 0, 0, 0, 0, -g(1) * (t15 * t46 + t6 * t50) - g(2) * (t13 * t46 + t5 * t50) - g(3) * (t17 * t50 + t26 * t46) -g(1) * (t15 * t50 - t6 * t46) - g(2) * (t13 * t50 - t5 * t46) - g(3) * (-t17 * t46 + t26 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t55, 0, 0, 0, 0, 0, t56 * t51, -t56 * t47, -t55, -g(1) * (-t11 * t38 - t12 * t45) - g(2) * (-t10 * t45 - t9 * t38) - g(3) * (-t18 * t38 - t19 * t45) 0, 0, 0, 0, 0, -g(1) * (-t11 * t77 + t12 * t46) - g(2) * (t10 * t46 - t9 * t77) - g(3) * (-t18 * t77 + t19 * t46) -g(1) * (t11 * t78 + t12 * t50) - g(2) * (t10 * t50 + t9 * t78) - g(3) * (t18 * t78 + t19 * t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -g(1) * (-t12 * t51 - t21 * t47) - g(2) * (-t10 * t51 - t20 * t47) - g(3) * (-t19 * t51 - t28 * t47) 0, t53 * pkin(4), 0, 0, 0, 0, 0, -t57 * t50, t57 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t50 - t4 * t46) - g(2) * (-t2 * t46 + t9 * t50) - g(3) * (t18 * t50 - t8 * t46) -g(1) * (-t11 * t46 - t4 * t50) - g(2) * (-t2 * t50 - t9 * t46) - g(3) * (-t18 * t46 - t8 * t50);];
taug_reg  = t1;
