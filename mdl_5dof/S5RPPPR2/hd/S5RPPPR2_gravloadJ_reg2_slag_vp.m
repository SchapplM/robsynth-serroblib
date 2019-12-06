% Calculate inertial parameters regressor of gravitation load for
% S5RPPPR2
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = cos(pkin(7));
t23 = sin(pkin(8));
t31 = cos(qJ(1));
t39 = t31 * t23;
t26 = cos(pkin(8));
t29 = sin(qJ(1));
t40 = t29 * t26;
t15 = t27 * t39 - t40;
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t38 = t31 * t26;
t41 = t29 * t23;
t16 = t27 * t38 + t41;
t22 = sin(pkin(9));
t25 = cos(pkin(9));
t24 = sin(pkin(7));
t42 = t24 * t31;
t6 = t16 * t25 + t22 * t42;
t49 = -t15 * t30 + t6 * t28;
t48 = t15 * t28 + t6 * t30;
t47 = g(2) * t31;
t44 = t23 * t24;
t43 = t24 * t29;
t37 = -t16 * pkin(3) - t15 * qJ(4);
t14 = -t27 * t40 + t39;
t4 = t14 * t22 + t25 * t43;
t7 = -t16 * t22 + t25 * t42;
t36 = g(2) * t7 + g(3) * t4;
t13 = t27 * t41 + t38;
t21 = t31 * qJ(2);
t35 = t14 * pkin(3) - t13 * qJ(4) + t21;
t34 = g(2) * t15 + g(3) * t13;
t19 = g(3) * t29 + t47;
t18 = g(2) * t29 - g(3) * t31;
t33 = pkin(2) * t27 + qJ(3) * t24 + pkin(1);
t32 = (g(2) * qJ(2) + g(3) * t33) * t29 + t33 * t47;
t17 = t19 * t24;
t12 = t24 * t26 * t25 - t27 * t22;
t9 = g(1) * t27 + t18 * t24;
t5 = t14 * t25 - t22 * t43;
t3 = -g(1) * t44 + g(2) * t13 - g(3) * t15;
t2 = -t13 * t28 + t5 * t30;
t1 = -t13 * t30 - t5 * t28;
t8 = [0, 0, 0, 0, 0, 0, t19, -t18, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t27, -t17, t18, -g(2) * (-t31 * pkin(1) - t29 * qJ(2)) - g(3) * (-t29 * pkin(1) + t21), 0, 0, 0, 0, 0, 0, g(2) * t16 - g(3) * t14, -t34, t17, -g(3) * t21 + t32, 0, 0, 0, 0, 0, 0, g(2) * t6 - g(3) * t5, t36, t34, -g(2) * t37 - g(3) * t35 + t32, 0, 0, 0, 0, 0, 0, g(2) * t48 - g(3) * t2, -g(2) * t49 - g(3) * t1, -t36, -g(2) * (-pkin(4) * t6 + t7 * pkin(6) + t37) - g(3) * (t5 * pkin(4) + t4 * pkin(6) + t35) + t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t28 + t30 * t44) - g(2) * t1 + g(3) * t49, -g(1) * (-t12 * t30 - t28 * t44) + g(2) * t2 + g(3) * t48, 0, 0;];
taug_reg = t8;
