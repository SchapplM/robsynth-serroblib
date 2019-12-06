% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = qJ(1) + qJ(2);
t18 = sin(t20);
t43 = pkin(2) * t18;
t19 = cos(t20);
t42 = pkin(2) * t19;
t21 = sin(pkin(9));
t41 = g(1) * t21;
t17 = pkin(8) + t20;
t16 = cos(t17);
t40 = g(2) * t16;
t39 = g(2) * t19;
t24 = sin(qJ(1));
t38 = t24 * pkin(1);
t26 = cos(qJ(1));
t37 = t26 * pkin(1);
t22 = cos(pkin(9));
t23 = sin(qJ(5));
t36 = t22 * t23;
t25 = cos(qJ(5));
t35 = t22 * t25;
t13 = t16 * qJ(4);
t34 = t13 - t43;
t33 = -t38 - t43;
t15 = sin(t17);
t10 = g(3) * t15 + t40;
t12 = g(3) * t18 + t39;
t32 = g(2) * t26 + g(3) * t24;
t31 = pkin(4) * t22 + pkin(7) * t21 + pkin(3);
t30 = -t15 * pkin(3) + t34;
t29 = g(2) * (-t37 - t42);
t28 = -t16 * pkin(3) - t15 * qJ(4) - t42;
t27 = (g(2) * qJ(4) + g(3) * t31) * t15 + t31 * t40;
t11 = -g(2) * t18 + g(3) * t19;
t9 = g(2) * t15 - g(3) * t16;
t8 = t10 * t22;
t7 = t10 * t21;
t6 = -t15 * t23 - t16 * t35;
t5 = -t15 * t25 + t16 * t36;
t4 = t15 * t35 - t16 * t23;
t3 = t15 * t36 + t16 * t25;
t2 = -g(2) * t6 + g(3) * t4;
t1 = -g(2) * t5 - g(3) * t3;
t14 = [0, 0, 0, 0, 0, 0, t32, -g(2) * t24 + g(3) * t26, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, t32 * pkin(1), 0, 0, 0, 0, 0, 0, t10, -t9, 0, -g(3) * t33 - t29, 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * (t28 - t37) - g(3) * (t30 - t38), 0, 0, 0, 0, 0, 0, t2, t1, t7, -t29 - g(3) * (t13 + t33) + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t12 * pkin(2), 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * t28 - g(3) * t30, 0, 0, 0, 0, 0, 0, t2, t1, t7, pkin(2) * t39 - g(3) * t34 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 + g(3) * t5 + t23 * t41, -g(2) * t4 - g(3) * t6 + t25 * t41, 0, 0;];
taug_reg = t14;
