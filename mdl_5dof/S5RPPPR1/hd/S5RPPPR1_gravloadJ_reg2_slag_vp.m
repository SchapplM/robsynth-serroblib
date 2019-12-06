% Calculate inertial parameters regressor of gravitation load for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t17 = sin(pkin(9));
t36 = pkin(4) * t17;
t18 = sin(pkin(8));
t35 = g(1) * t18;
t16 = qJ(1) + pkin(7);
t14 = cos(t16);
t34 = g(2) * t14;
t23 = cos(qJ(1));
t33 = t23 * pkin(1);
t12 = sin(t16);
t20 = cos(pkin(8));
t32 = t12 * t20;
t31 = t14 * t20;
t30 = t17 * t20;
t19 = cos(pkin(9));
t29 = t19 * t20;
t22 = sin(qJ(1));
t28 = -t22 * pkin(1) + t14 * qJ(3);
t8 = g(3) * t12 + t34;
t7 = g(2) * t12 - g(3) * t14;
t27 = g(2) * t23 + g(3) * t22;
t26 = pkin(3) * t20 + qJ(4) * t18 + pkin(2);
t25 = (t19 * pkin(4) + pkin(3)) * t20 - t18 * (-pkin(6) - qJ(4)) + pkin(2);
t24 = g(2) * t33 - g(3) * t28;
t15 = pkin(9) + qJ(5);
t13 = cos(t15);
t11 = sin(t15);
t6 = t8 * t18;
t5 = g(1) * t20 + t7 * t18;
t4 = -t12 * t11 - t13 * t31;
t3 = t11 * t31 - t12 * t13;
t2 = -t14 * t11 + t13 * t32;
t1 = t11 * t32 + t14 * t13;
t9 = [0, 0, 0, 0, 0, 0, t27, -g(2) * t22 + g(3) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t27 * pkin(1), 0, 0, 0, 0, 0, 0, t8 * t20, -t6, t7, -g(2) * (-t14 * pkin(2) - t12 * qJ(3) - t33) - g(3) * (-t12 * pkin(2) + t28), 0, 0, 0, 0, 0, 0, -g(2) * (-t12 * t17 - t14 * t29) - g(3) * (-t12 * t29 + t14 * t17), -g(2) * (-t12 * t19 + t14 * t30) - g(3) * (t12 * t30 + t14 * t19), t6, t26 * t34 + (g(2) * qJ(3) + g(3) * t26) * t12 + t24, 0, 0, 0, 0, 0, 0, -g(2) * t4 + g(3) * t2, -g(2) * t3 - g(3) * t1, t6, (g(2) * t25 - g(3) * t36) * t14 + (-g(2) * (-qJ(3) - t36) + g(3) * t25) * t12 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t1 + g(3) * t3 + t11 * t35, -g(2) * t2 - g(3) * t4 + t13 * t35, 0, 0;];
taug_reg = t9;
