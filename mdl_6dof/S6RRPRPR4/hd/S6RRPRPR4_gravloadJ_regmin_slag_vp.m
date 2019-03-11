% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(pkin(11));
t47 = sin(qJ(2));
t51 = cos(qJ(2));
t66 = cos(pkin(11));
t30 = -t41 * t51 - t47 * t66;
t48 = sin(qJ(1));
t52 = cos(qJ(1));
t43 = cos(pkin(6));
t57 = -t41 * t47 + t51 * t66;
t54 = t57 * t43;
t11 = t48 * t30 + t52 * t54;
t45 = sin(qJ(6));
t49 = cos(qJ(6));
t23 = t30 * t43;
t12 = -t23 * t52 + t48 * t57;
t40 = qJ(4) + pkin(12);
t38 = sin(t40);
t39 = cos(t40);
t42 = sin(pkin(6));
t73 = t42 * t52;
t5 = -t12 * t39 + t38 * t73;
t85 = -t11 * t49 + t45 * t5;
t84 = t11 * t45 + t49 * t5;
t22 = t30 * t42;
t46 = sin(qJ(4));
t50 = cos(qJ(4));
t58 = t12 * t46 + t50 * t73;
t13 = -t23 * t48 - t52 * t57;
t75 = t42 * t48;
t8 = t13 * t46 + t50 * t75;
t83 = g(2) * t58 - g(3) * (t22 * t46 + t43 * t50) - g(1) * t8;
t68 = t52 * t47;
t69 = t48 * t51;
t27 = -t43 * t69 - t68;
t74 = t42 * t51;
t82 = -g(1) * t27 - g(3) * t74;
t77 = t39 * t45;
t76 = t39 * t49;
t37 = pkin(2) * t51 + pkin(1);
t71 = t48 * t37;
t70 = t48 * t47;
t67 = t52 * t51;
t64 = t43 * t67;
t63 = t12 * t50 - t46 * t73;
t24 = t43 * t47 * pkin(2) + (-pkin(8) - qJ(3)) * t42;
t61 = pkin(4) * t42 * t46 - t24;
t60 = g(1) * t52 + g(2) * t48;
t59 = g(1) * t48 - g(2) * t52;
t56 = g(1) * (t13 * t38 + t39 * t75) + g(2) * (-t12 * t38 - t39 * t73) + g(3) * (t22 * t38 + t39 * t43);
t14 = t30 * t52 - t48 * t54;
t21 = t57 * t42;
t55 = g(1) * t14 + g(2) * t11 + g(3) * t21;
t44 = -qJ(5) - pkin(9);
t36 = pkin(4) * t50 + pkin(3);
t33 = t52 * t37;
t31 = pkin(2) * t64;
t28 = -t43 * t70 + t67;
t26 = -t43 * t68 - t69;
t25 = -t64 + t70;
t20 = -g(3) * t43 - t42 * t59;
t17 = -t22 * t39 + t38 * t43;
t9 = -t13 * t50 + t46 * t75;
t7 = -t13 * t39 + t38 * t75;
t2 = -t14 * t45 + t49 * t7;
t1 = -t14 * t49 - t45 * t7;
t3 = [0, t59, t60, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t28, -g(1) * t25 - g(2) * t27, -t60 * t42, -g(1) * (-t24 * t52 - t71) - g(2) * (-t24 * t48 + t33) 0, 0, 0, 0, 0, g(1) * t63 - g(2) * t9, -g(1) * t58 - g(2) * t8, -g(1) * t11 + g(2) * t14, -g(1) * (-t11 * t44 - t12 * t36 + t52 * t61 - t71) - g(2) * (-t13 * t36 + t14 * t44 + t48 * t61 + t33) 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t2, g(1) * t85 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t25 + t82, g(3) * t42 * t47 + g(1) * t28 - g(2) * t26, 0, -g(2) * t31 + (g(2) * t70 + t82) * pkin(2), 0, 0, 0, 0, 0, -t55 * t50, t55 * t46, g(1) * t13 - g(2) * t12 + g(3) * t22, -g(1) * (pkin(2) * t27 + t13 * t44 + t14 * t36) - g(2) * (-pkin(2) * t70 + t11 * t36 - t12 * t44 + t31) - g(3) * (pkin(2) * t74 + t21 * t36 + t22 * t44) 0, 0, 0, 0, 0, -g(1) * (-t13 * t45 + t14 * t76) - g(2) * (t11 * t76 + t12 * t45) - g(3) * (t21 * t76 - t22 * t45) -g(1) * (-t13 * t49 - t14 * t77) - g(2) * (-t11 * t77 + t12 * t49) - g(3) * (-t21 * t77 - t22 * t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, g(1) * t9 + g(2) * t63 - g(3) * (t22 * t50 - t43 * t46) 0, t83 * pkin(4), 0, 0, 0, 0, 0, -t56 * t49, t56 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t85 - g(3) * (-t17 * t45 - t21 * t49) g(1) * t2 - g(2) * t84 - g(3) * (-t17 * t49 + t21 * t45);];
taug_reg  = t3;
