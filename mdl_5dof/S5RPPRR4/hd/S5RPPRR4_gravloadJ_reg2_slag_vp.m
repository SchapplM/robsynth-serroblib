% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = pkin(9) + qJ(4);
t20 = sin(t28);
t30 = sin(pkin(8));
t55 = g(1) * t30;
t32 = cos(pkin(8));
t21 = cos(t28);
t35 = cos(qJ(1));
t42 = t35 * t21;
t34 = sin(qJ(1));
t49 = t34 * t20;
t7 = -t32 * t49 - t42;
t43 = t35 * t20;
t48 = t34 * t21;
t9 = t32 * t43 - t48;
t58 = -g(2) * t7 - g(3) * t9 + t20 * t55;
t29 = sin(pkin(9));
t54 = t29 * pkin(3);
t31 = cos(pkin(9));
t19 = t31 * pkin(3) + pkin(2);
t53 = t30 * t34;
t52 = t32 * t34;
t22 = qJ(5) + t28;
t17 = sin(t22);
t51 = t34 * t17;
t18 = cos(t22);
t50 = t34 * t18;
t47 = t34 * t29;
t46 = t34 * t31;
t45 = t35 * t17;
t44 = t35 * t18;
t41 = t35 * t29;
t40 = t35 * t31;
t33 = -pkin(6) - qJ(3);
t39 = t35 * pkin(1) + t34 * qJ(2);
t25 = t34 * pkin(1);
t37 = -t35 * qJ(2) + t25;
t16 = g(2) * t35 + g(3) * t34;
t15 = g(2) * t34 - g(3) * t35;
t36 = pkin(2) * t32 + qJ(3) * t30;
t27 = -pkin(7) + t33;
t14 = pkin(4) * t20 + t54;
t13 = pkin(4) * t21 + t19;
t12 = t16 * t30;
t11 = g(1) * t32 - t15 * t30;
t10 = t32 * t42 + t49;
t8 = t32 * t48 - t43;
t6 = t32 * t44 + t51;
t5 = t32 * t45 - t50;
t4 = t32 * t50 - t45;
t3 = -t32 * t51 - t44;
t2 = g(2) * t4 - g(3) * t6 + t18 * t55;
t1 = -g(2) * t3 - g(3) * t5 + t17 * t55;
t23 = [0, 0, 0, 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t32, t12, -t15, -g(2) * t39 - g(3) * t37, 0, 0, 0, 0, 0, 0, -g(2) * (t32 * t40 + t47) - g(3) * (t32 * t46 - t41), -g(2) * (-t32 * t41 + t46) - g(3) * (-t32 * t47 - t40), -t12, -g(2) * (t35 * t36 + t39) - g(3) * (t34 * t36 + t37), 0, 0, 0, 0, 0, 0, -g(2) * t10 - g(3) * t8, g(2) * t9 - g(3) * t7, -t12, -g(2) * (pkin(3) * t47 + t39) - g(3) * (t19 * t52 - t33 * t53 + t25) + (-g(2) * (t19 * t32 - t30 * t33) - g(3) * (-qJ(2) - t54)) * t35, 0, 0, 0, 0, 0, 0, -g(2) * t6 - g(3) * t4, g(2) * t5 - g(3) * t3, -t12, -g(2) * (t34 * t14 + t39) - g(3) * (t13 * t52 - t27 * t53 + t25) + (-g(2) * (t13 * t32 - t27 * t30) - g(3) * (-qJ(2) - t14)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, g(2) * t8 - g(3) * t10 + t21 * t55, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t58 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t23;
