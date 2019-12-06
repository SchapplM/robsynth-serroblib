% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = -pkin(7) - pkin(6);
t29 = qJ(2) + qJ(3);
t23 = sin(t29);
t40 = pkin(3) * t23;
t30 = sin(qJ(2));
t39 = t30 * pkin(2);
t24 = cos(t29);
t19 = pkin(3) * t24;
t32 = cos(qJ(2));
t26 = t32 * pkin(2);
t38 = t19 + t26;
t28 = -qJ(4) + t34;
t22 = pkin(9) + t29;
t18 = cos(t22);
t14 = pkin(4) * t18;
t37 = t14 + t38;
t17 = sin(t22);
t36 = -pkin(4) * t17 - t40;
t31 = sin(qJ(1));
t33 = cos(qJ(1));
t13 = g(1) * t33 + g(2) * t31;
t12 = g(1) * t31 - g(2) * t33;
t5 = -g(3) * t24 + t13 * t23;
t35 = -g(3) * t32 + t13 * t30;
t25 = -pkin(8) + t28;
t21 = t26 + pkin(1);
t20 = qJ(5) + t22;
t16 = cos(t20);
t15 = sin(t20);
t10 = pkin(1) + t38;
t7 = pkin(1) + t37;
t6 = g(3) * t23 + t13 * t24;
t4 = g(3) * t17 + t13 * t18;
t3 = -g(3) * t18 + t13 * t17;
t2 = g(3) * t15 + t13 * t16;
t1 = -g(3) * t16 + t13 * t15;
t8 = [0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t32, -t12 * t30, -t13, -g(1) * (-t31 * pkin(1) + t33 * pkin(6)) - g(2) * (t33 * pkin(1) + t31 * pkin(6)), 0, 0, 0, 0, 0, 0, t12 * t24, -t12 * t23, -t13, -g(1) * (-t31 * t21 - t33 * t34) - g(2) * (t33 * t21 - t31 * t34), 0, 0, 0, 0, 0, 0, t12 * t18, -t12 * t17, -t13, -g(1) * (-t31 * t10 - t33 * t28) - g(2) * (t33 * t10 - t31 * t28), 0, 0, 0, 0, 0, 0, t12 * t16, -t12 * t15, -t13, -g(1) * (-t33 * t25 - t31 * t7) - g(2) * (-t31 * t25 + t33 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, g(3) * t30 + t13 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t35 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t38 - t13 * (-t39 - t40), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t37 - t13 * (t36 - t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * (t14 + t19) - t13 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t8;
