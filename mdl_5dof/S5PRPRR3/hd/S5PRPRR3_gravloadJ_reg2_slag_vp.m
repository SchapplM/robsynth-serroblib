% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t26 = g(1) * t14 + g(2) * t13;
t11 = qJ(2) + pkin(9);
t6 = sin(t11);
t7 = cos(t11);
t23 = -g(3) * t7 + t26 * t6;
t39 = g(3) * t6;
t16 = sin(qJ(2));
t37 = pkin(2) * t16;
t12 = qJ(4) + qJ(5);
t8 = sin(t12);
t34 = t13 * t8;
t9 = cos(t12);
t33 = t13 * t9;
t32 = t14 * t8;
t31 = t14 * t9;
t15 = sin(qJ(4));
t30 = t13 * t15;
t17 = cos(qJ(4));
t29 = t13 * t17;
t28 = t14 * t15;
t27 = t14 * t17;
t18 = cos(qJ(2));
t21 = -g(3) * t18 + t26 * t16;
t20 = -g(1) * (-t7 * t28 + t29) - g(2) * (-t7 * t30 - t27) + t15 * t39;
t19 = -pkin(7) - pkin(6);
t10 = t18 * pkin(2);
t5 = t17 * pkin(4) + pkin(3);
t4 = -g(1) * t13 + g(2) * t14;
t3 = t26 * t7 + t39;
t2 = -g(1) * (-t7 * t31 - t34) - g(2) * (-t7 * t33 + t32) + t9 * t39;
t1 = -g(1) * (-t7 * t32 + t33) - g(2) * (-t7 * t34 - t31) + t8 * t39;
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, g(3) * t16 + t26 * t18, 0, 0, 0, 0, 0, 0, 0, 0, t23, t3, 0, t21 * pkin(2), 0, 0, 0, 0, 0, 0, t23 * t17, -t23 * t15, -t3, -g(3) * (t7 * pkin(3) + t6 * pkin(6) + t10) + t26 * (pkin(3) * t6 - pkin(6) * t7 + t37), 0, 0, 0, 0, 0, 0, t23 * t9, -t23 * t8, -t3, -g(3) * (-t6 * t19 + t7 * t5 + t10) + t26 * (t19 * t7 + t5 * t6 + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -g(1) * (-t7 * t27 - t30) - g(2) * (-t7 * t29 + t28) + t17 * t39, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t22;
