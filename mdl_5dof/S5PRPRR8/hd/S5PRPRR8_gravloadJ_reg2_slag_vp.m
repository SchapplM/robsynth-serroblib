% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = sin(pkin(9));
t25 = cos(pkin(9));
t31 = cos(qJ(2));
t28 = sin(qJ(2));
t42 = cos(pkin(5));
t39 = t28 * t42;
t12 = t23 * t31 + t25 * t39;
t14 = -t23 * t39 + t25 * t31;
t57 = -g(1) * t14 - g(2) * t12;
t24 = sin(pkin(5));
t54 = g(3) * t24;
t38 = t31 * t42;
t11 = t23 * t28 - t25 * t38;
t53 = t11 * pkin(7);
t13 = t23 * t38 + t25 * t28;
t52 = t13 * pkin(7);
t27 = sin(qJ(4));
t51 = t24 * t27;
t50 = t24 * t28;
t30 = cos(qJ(4));
t49 = t24 * t30;
t48 = t24 * t31;
t26 = sin(qJ(5));
t47 = t26 * t27;
t46 = t26 * t28;
t29 = cos(qJ(5));
t45 = t27 * t29;
t44 = t28 * t29;
t43 = pkin(2) * t48 + qJ(3) * t50;
t41 = pkin(7) * t48 + t43;
t9 = t11 * pkin(2);
t40 = t12 * qJ(3) - t9;
t10 = t13 * pkin(2);
t37 = t14 * qJ(3) - t10;
t36 = pkin(4) * t27 - pkin(8) * t30;
t15 = -t42 * t27 - t30 * t48;
t4 = t13 * t30 - t23 * t51;
t6 = t11 * t30 + t25 * t51;
t34 = g(1) * t4 + g(2) * t6 + g(3) * t15;
t16 = -t27 * t48 + t42 * t30;
t5 = t13 * t27 + t23 * t49;
t7 = -t11 * t27 + t25 * t49;
t33 = g(1) * t5 - g(2) * t7 + g(3) * t16;
t2 = -g(1) * t13 - g(2) * t11 + g(3) * t48;
t32 = g(3) * t50 - t57;
t1 = t32 * t30;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t32, -g(1) * t37 - g(2) * t40 - g(3) * t43, 0, 0, 0, 0, 0, 0, -t32 * t27, -t1, -t2, -g(1) * (t37 - t52) - g(2) * (t40 - t53) - g(3) * t41, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t26 + t14 * t45) - g(2) * (-t11 * t26 + t12 * t45) - (t26 * t31 + t27 * t44) * t54, -g(1) * (-t13 * t29 - t14 * t47) - g(2) * (-t11 * t29 - t12 * t47) - (-t27 * t46 + t29 * t31) * t54, t1, -g(1) * (-t10 - t52) - g(2) * (-t9 - t53) - g(3) * (t36 * t50 + t41) + t57 * (qJ(3) + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t29, t34 * t26, -t33, -g(1) * (t4 * pkin(4) + t5 * pkin(8)) - g(2) * (t6 * pkin(4) - t7 * pkin(8)) - g(3) * (t15 * pkin(4) + t16 * pkin(8)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t29 - t5 * t26) - g(2) * (t12 * t29 + t7 * t26) - g(3) * (-t16 * t26 + t24 * t44), -g(1) * (-t14 * t26 - t5 * t29) - g(2) * (-t12 * t26 + t7 * t29) - g(3) * (-t16 * t29 - t24 * t46), 0, 0;];
taug_reg = t3;
