% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t32 = qJ(3) + pkin(11);
t28 = qJ(5) + t32;
t21 = sin(t28);
t22 = cos(t28);
t46 = t22 * pkin(5) + t21 * pkin(9);
t33 = qJ(1) + pkin(10);
t25 = sin(t33);
t27 = cos(t33);
t12 = g(1) * t27 + g(2) * t25;
t3 = -g(3) * t22 + t12 * t21;
t57 = pkin(5) * t21;
t56 = g(3) * t21;
t37 = sin(qJ(1));
t54 = t37 * pkin(1);
t53 = t21 * t27;
t52 = t22 * t27;
t35 = sin(qJ(6));
t51 = t25 * t35;
t38 = cos(qJ(6));
t50 = t25 * t38;
t49 = t27 * t35;
t48 = t27 * t38;
t34 = -qJ(4) - pkin(7);
t26 = cos(t32);
t39 = cos(qJ(3));
t29 = t39 * pkin(3);
t45 = pkin(4) * t26 + t29;
t15 = pkin(2) + t45;
t40 = cos(qJ(1));
t30 = t40 * pkin(1);
t47 = t27 * t15 + t30;
t24 = sin(t32);
t36 = sin(qJ(3));
t16 = -t36 * pkin(3) - pkin(4) * t24;
t44 = t16 - t57;
t11 = g(1) * t25 - g(2) * t27;
t43 = g(1) * t37 - g(2) * t40;
t31 = -pkin(8) + t34;
t42 = -t27 * t31 - t54;
t41 = -g(3) * t39 + t12 * t36;
t23 = t29 + pkin(2);
t14 = pkin(9) * t52;
t13 = t25 * t22 * pkin(9);
t9 = t22 * t48 + t51;
t8 = -t22 * t49 + t50;
t7 = -t22 * t50 + t49;
t6 = t22 * t51 + t48;
t5 = t11 * t21;
t4 = t12 * t22 + t56;
t2 = t3 * t38;
t1 = t3 * t35;
t10 = [0, 0, 0, 0, 0, 0, t43, g(1) * t40 + g(2) * t37, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t43 * pkin(1), 0, 0, 0, 0, 0, 0, t11 * t39, -t11 * t36, -t12, -g(1) * (-t25 * pkin(2) + t27 * pkin(7) - t54) - g(2) * (t27 * pkin(2) + t25 * pkin(7) + t30) 0, 0, 0, 0, 0, 0, t11 * t26, -t11 * t24, -t12, -g(1) * (-t25 * t23 - t27 * t34 - t54) - g(2) * (t27 * t23 - t25 * t34 + t30) 0, 0, 0, 0, 0, 0, t11 * t22, -t5, -t12, -g(1) * (-t25 * t15 + t42) - g(2) * (-t25 * t31 + t47) 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5, -g(1) * t42 - g(2) * (pkin(5) * t52 + pkin(9) * t53 + t47) + (-g(1) * (-t15 - t46) + g(2) * t31) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, g(3) * t36 + t12 * t39, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t12 * t24, g(3) * t24 + t12 * t26, 0, t41 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t45 - t12 * t16, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t44 * t27 + t14) - g(2) * (t44 * t25 + t13) - g(3) * (t45 + t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-pkin(5) * t53 + t14) - g(2) * (-t25 * t57 + t13) - g(3) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t6 + t35 * t56, g(1) * t9 - g(2) * t7 + t38 * t56, 0, 0;];
taug_reg  = t10;
