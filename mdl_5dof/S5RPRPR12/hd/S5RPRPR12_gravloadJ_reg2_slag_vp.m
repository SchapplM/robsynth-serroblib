% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = sin(qJ(1));
t28 = cos(qJ(1));
t10 = g(1) * t28 + g(2) * t27;
t20 = pkin(8) + qJ(3);
t15 = sin(t20);
t17 = cos(t20);
t1 = -g(3) * t17 + t10 * t15;
t24 = cos(pkin(8));
t13 = t24 * pkin(2) + pkin(1);
t8 = t28 * t13;
t47 = g(2) * t8;
t44 = g(3) * t15;
t19 = pkin(9) + qJ(5);
t14 = sin(t19);
t42 = t27 * t14;
t16 = cos(t19);
t41 = t27 * t16;
t21 = sin(pkin(9));
t40 = t27 * t21;
t23 = cos(pkin(9));
t39 = t27 * t23;
t38 = t28 * t14;
t37 = t28 * t16;
t36 = t28 * t21;
t35 = t28 * t23;
t26 = -pkin(6) - qJ(2);
t34 = pkin(4) * t21 - t26;
t9 = g(1) * t27 - g(2) * t28;
t33 = t17 * pkin(3) + t15 * qJ(4);
t12 = t23 * pkin(4) + pkin(3);
t25 = -pkin(7) - qJ(4);
t31 = t17 * t12 - t15 * t25;
t7 = t9 * t15;
t6 = t17 * t37 + t42;
t5 = -t17 * t38 + t41;
t4 = -t17 * t41 + t38;
t3 = t17 * t42 + t37;
t2 = t10 * t17 + t44;
t11 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t24, -t9 * sin(pkin(8)), -t10, -g(1) * (-t27 * pkin(1) + t28 * qJ(2)) - g(2) * (t28 * pkin(1) + t27 * qJ(2)), 0, 0, 0, 0, 0, 0, t9 * t17, -t7, -t10, -g(1) * (-t27 * t13 - t28 * t26) - g(2) * (-t27 * t26 + t8), 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t39 + t36) - g(2) * (t17 * t35 + t40), -g(1) * (t17 * t40 + t35) - g(2) * (-t17 * t36 + t39), t7, -t47 + (g(1) * t26 - g(2) * t33) * t28 + (-g(1) * (-t13 - t33) + g(2) * t26) * t27, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t7, -t47 + (-g(1) * t34 - g(2) * t31) * t28 + (-g(1) * (-t13 - t31) - g(2) * t34) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t23, -t1 * t21, -t2, -g(3) * t33 + t10 * (pkin(3) * t15 - qJ(4) * t17), 0, 0, 0, 0, 0, 0, t1 * t16, -t1 * t14, -t2, -g(3) * t31 + t10 * (t12 * t15 + t17 * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t14 * t44, g(1) * t6 - g(2) * t4 + t16 * t44, 0, 0;];
taug_reg = t11;
