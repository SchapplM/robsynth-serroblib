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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
t26 = pkin(9) + qJ(4);
t20 = sin(t26);
t28 = sin(pkin(8));
t53 = g(1) * t28;
t30 = cos(pkin(8));
t21 = cos(t26);
t33 = cos(qJ(1));
t40 = t33 * t21;
t32 = sin(qJ(1));
t47 = t32 * t20;
t7 = t30 * t47 + t40;
t41 = t33 * t20;
t46 = t32 * t21;
t9 = t30 * t41 - t46;
t56 = -g(2) * t7 + g(3) * t9 + t20 * t53;
t52 = g(2) * t33;
t24 = t33 * qJ(2);
t51 = g(3) * t24;
t27 = sin(pkin(9));
t50 = t27 * pkin(3);
t29 = cos(pkin(9));
t19 = t29 * pkin(3) + pkin(2);
t22 = qJ(5) + t26;
t17 = sin(t22);
t49 = t32 * t17;
t18 = cos(t22);
t48 = t32 * t18;
t45 = t32 * t27;
t44 = t32 * t29;
t43 = t33 * t17;
t42 = t33 * t18;
t39 = t33 * t27;
t38 = t33 * t29;
t31 = -pkin(6) - qJ(3);
t16 = g(3) * t32 + t52;
t15 = g(2) * t32 - g(3) * t33;
t36 = pkin(2) * t30 + qJ(3) * t28 + pkin(1);
t35 = (pkin(4) * t21 + t19) * t30 - (-pkin(7) + t31) * t28 + pkin(1);
t34 = t19 * t30 - t28 * t31 + pkin(1);
t14 = pkin(4) * t20 + t50;
t12 = t16 * t28;
t11 = g(1) * t30 + t15 * t28;
t10 = -t30 * t40 - t47;
t8 = t30 * t46 - t41;
t6 = -t30 * t42 - t49;
t5 = t30 * t43 - t48;
t4 = t30 * t48 - t43;
t3 = t30 * t49 + t42;
t2 = -g(2) * t4 - g(3) * t6 + t18 * t53;
t1 = -g(2) * t3 + g(3) * t5 + t17 * t53;
t13 = [0, 0, 0, 0, 0, 0, t16, -t15, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t30, -t12, t15, -g(2) * (-t33 * pkin(1) - t32 * qJ(2)) - g(3) * (-t32 * pkin(1) + t24), 0, 0, 0, 0, 0, 0, -g(2) * (-t30 * t38 - t45) - g(3) * (-t30 * t44 + t39), -g(2) * (t30 * t39 - t44) - g(3) * (t30 * t45 + t38), t12, -t51 + t36 * t52 + (g(2) * qJ(2) + g(3) * t36) * t32, 0, 0, 0, 0, 0, 0, -g(2) * t10 + g(3) * t8, -g(2) * t9 - g(3) * t7, t12, -t51 + (g(2) * t34 - g(3) * t50) * t33 + (-g(2) * (-qJ(2) - t50) + g(3) * t34) * t32, 0, 0, 0, 0, 0, 0, -g(2) * t6 + g(3) * t4, -g(2) * t5 - g(3) * t3, t12, -t51 + (g(2) * t35 - g(3) * t14) * t33 + (-g(2) * (-qJ(2) - t14) + g(3) * t35) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -g(2) * t8 - g(3) * t10 + t21 * t53, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t56 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t13;
