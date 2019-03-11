% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t33 = qJ(3) + qJ(4);
t28 = qJ(5) + t33;
t21 = sin(t28);
t22 = cos(t28);
t49 = pkin(5) * t22 + pkin(10) * t21;
t31 = qJ(1) + pkin(11);
t24 = sin(t31);
t25 = cos(t31);
t13 = g(1) * t25 + g(2) * t24;
t3 = -g(3) * t22 + t13 * t21;
t40 = -pkin(8) - pkin(7);
t26 = sin(t33);
t61 = pkin(4) * t26;
t60 = pkin(5) * t21;
t59 = g(3) * t21;
t36 = sin(qJ(1));
t57 = t36 * pkin(1);
t56 = t21 * t25;
t55 = t22 * t25;
t34 = sin(qJ(6));
t54 = t24 * t34;
t37 = cos(qJ(6));
t53 = t24 * t37;
t52 = t25 * t34;
t51 = t25 * t37;
t27 = cos(t33);
t20 = pkin(4) * t27;
t38 = cos(qJ(3));
t29 = t38 * pkin(3);
t48 = t20 + t29;
t16 = pkin(2) + t48;
t39 = cos(qJ(1));
t30 = t39 * pkin(1);
t50 = t16 * t25 + t30;
t47 = t20 + t49;
t35 = sin(qJ(3));
t17 = -pkin(3) * t35 - t61;
t46 = t17 - t60;
t45 = -t60 - t61;
t44 = g(1) * t24 - g(2) * t25;
t43 = g(1) * t36 - g(2) * t39;
t32 = -pkin(9) + t40;
t42 = -t25 * t32 - t57;
t5 = -g(3) * t27 + t13 * t26;
t41 = -g(3) * t38 + t13 * t35;
t23 = t29 + pkin(2);
t15 = pkin(10) * t55;
t14 = t24 * t22 * pkin(10);
t11 = t22 * t51 + t54;
t10 = -t22 * t52 + t53;
t9 = -t22 * t53 + t52;
t8 = t22 * t54 + t51;
t7 = t44 * t21;
t6 = g(3) * t26 + t13 * t27;
t4 = t13 * t22 + t59;
t2 = t3 * t37;
t1 = t3 * t34;
t12 = [0, 0, 0, 0, 0, 0, t43, g(1) * t39 + g(2) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t44, t13, 0, t43 * pkin(1), 0, 0, 0, 0, 0, 0, t44 * t38, -t44 * t35, -t13, -g(1) * (-pkin(2) * t24 + pkin(7) * t25 - t57) - g(2) * (pkin(2) * t25 + pkin(7) * t24 + t30) 0, 0, 0, 0, 0, 0, t44 * t27, -t44 * t26, -t13, -g(1) * (-t23 * t24 - t25 * t40 - t57) - g(2) * (t23 * t25 - t24 * t40 + t30) 0, 0, 0, 0, 0, 0, t44 * t22, -t7, -t13, -g(1) * (-t16 * t24 + t42) - g(2) * (-t24 * t32 + t50) 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, -g(1) * t42 - g(2) * (pkin(5) * t55 + pkin(10) * t56 + t50) + (-g(1) * (-t16 - t49) + g(2) * t32) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, g(3) * t35 + t13 * t38, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t41 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t48 - t13 * t17, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t25 * t46 + t15) - g(2) * (t24 * t46 + t14) - g(3) * (t29 + t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t25 * t45 + t15) - g(2) * (t24 * t45 + t14) - g(3) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-pkin(5) * t56 + t15) - g(2) * (-t24 * t60 + t14) - g(3) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t34 * t59, g(1) * t11 - g(2) * t9 + t37 * t59, 0, 0;];
taug_reg  = t12;
