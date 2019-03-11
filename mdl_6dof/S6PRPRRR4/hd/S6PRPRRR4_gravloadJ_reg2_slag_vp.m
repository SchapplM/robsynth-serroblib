% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t36 = cos(pkin(12));
t25 = t36 * pkin(3) + pkin(2);
t35 = sin(pkin(6));
t41 = cos(qJ(2));
t59 = t35 * t41;
t21 = t25 * t59;
t70 = g(3) * t21;
t69 = g(3) * t35;
t34 = sin(pkin(11));
t39 = sin(qJ(2));
t52 = cos(pkin(11));
t53 = cos(pkin(6));
t47 = t53 * t52;
t18 = t34 * t41 + t39 * t47;
t38 = sin(qJ(5));
t68 = t18 * t38;
t51 = t34 * t53;
t20 = -t39 * t51 + t52 * t41;
t67 = t20 * t38;
t31 = pkin(12) + qJ(4);
t28 = cos(t31);
t32 = qJ(5) + qJ(6);
t29 = sin(t32);
t66 = t28 * t29;
t30 = cos(t32);
t65 = t28 * t30;
t64 = t28 * t38;
t40 = cos(qJ(5));
t63 = t28 * t40;
t62 = t28 * t41;
t61 = t34 * t35;
t60 = t35 * t39;
t37 = -pkin(8) - qJ(3);
t58 = t37 * t39;
t57 = t38 * t41;
t56 = t40 * t41;
t17 = t34 * t39 - t41 * t47;
t55 = -t17 * t25 - t18 * t37;
t19 = t52 * t39 + t41 * t51;
t54 = -t19 * t25 - t20 * t37;
t50 = t35 * t52;
t27 = sin(t31);
t49 = pkin(4) * t28 + pkin(9) * t27;
t26 = t40 * pkin(5) + pkin(4);
t42 = -pkin(10) - pkin(9);
t48 = t26 * t28 - t27 * t42;
t13 = -t27 * t60 + t53 * t28;
t7 = -t18 * t27 - t28 * t50;
t9 = -t20 * t27 + t28 * t61;
t46 = g(1) * t9 + g(2) * t7 + g(3) * t13;
t10 = t20 * t28 + t27 * t61;
t14 = t53 * t27 + t28 * t60;
t8 = t18 * t28 - t27 * t50;
t45 = g(1) * t10 + g(2) * t8 + g(3) * t14;
t5 = -g(1) * t19 - g(2) * t17 + g(3) * t59;
t44 = g(1) * t20 + g(2) * t18 + g(3) * t60;
t43 = -g(1) * (-t10 * t38 + t19 * t40) - g(2) * (t17 * t40 - t8 * t38) - g(3) * (-t14 * t38 - t35 * t56);
t4 = t5 * t27;
t2 = -g(1) * (-t10 * t30 - t19 * t29) - g(2) * (-t17 * t29 - t8 * t30) - g(3) * (-t14 * t30 + t29 * t59);
t1 = -g(1) * (-t10 * t29 + t19 * t30) - g(2) * (t17 * t30 - t8 * t29) - g(3) * (-t14 * t29 - t30 * t59);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t44, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t36, t5 * sin(pkin(12)) -t44, -g(1) * (-t19 * pkin(2) + t20 * qJ(3)) - g(2) * (-t17 * pkin(2) + t18 * qJ(3)) - (pkin(2) * t41 + qJ(3) * t39) * t69, 0, 0, 0, 0, 0, 0, -t5 * t28, t4, -t44, -g(1) * t54 - g(2) * t55 - g(3) * (-t35 * t58 + t21) 0, 0, 0, 0, 0, 0, -g(1) * (-t19 * t63 + t67) - g(2) * (-t17 * t63 + t68) - (t28 * t56 + t38 * t39) * t69, -g(1) * (t19 * t64 + t20 * t40) - g(2) * (t17 * t64 + t18 * t40) - (-t28 * t57 + t39 * t40) * t69, -t4, -g(1) * (-t19 * t49 + t54) - g(2) * (-t17 * t49 + t55) - t70 - (t41 * t49 - t58) * t69, 0, 0, 0, 0, 0, 0, -g(1) * (-t19 * t65 + t20 * t29) - g(2) * (-t17 * t65 + t18 * t29) - (t29 * t39 + t30 * t62) * t69, -g(1) * (t19 * t66 + t20 * t30) - g(2) * (t17 * t66 + t18 * t30) - (-t29 * t62 + t30 * t39) * t69, -t4, -g(1) * (pkin(5) * t67 - t19 * t48 + t54) - g(2) * (pkin(5) * t68 - t17 * t48 + t55) - t70 - (t48 * t41 + (pkin(5) * t38 - t37) * t39) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t40, t46 * t38, -t45, -g(1) * (t9 * pkin(4) + t10 * pkin(9)) - g(2) * (t7 * pkin(4) + t8 * pkin(9)) - g(3) * (t13 * pkin(4) + t14 * pkin(9)) 0, 0, 0, 0, 0, 0, -t46 * t30, t46 * t29, -t45, -g(1) * (-t10 * t42 + t9 * t26) - g(2) * (t7 * t26 - t8 * t42) - g(3) * (t13 * t26 - t14 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -g(1) * (-t10 * t40 - t19 * t38) - g(2) * (-t17 * t38 - t8 * t40) - g(3) * (-t14 * t40 + t35 * t57) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t43 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
