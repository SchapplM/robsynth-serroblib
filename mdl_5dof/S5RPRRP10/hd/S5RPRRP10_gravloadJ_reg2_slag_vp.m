% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t16 = g(1) * t30 + g(2) * t28;
t22 = pkin(8) + qJ(3);
t20 = cos(t22);
t29 = cos(qJ(4));
t37 = t30 * t29;
t27 = sin(qJ(4));
t40 = t28 * t27;
t10 = t20 * t40 + t37;
t38 = t30 * t27;
t39 = t28 * t29;
t12 = -t20 * t38 + t39;
t19 = sin(t22);
t42 = g(3) * t19;
t1 = -g(1) * t12 + g(2) * t10 + t27 * t42;
t7 = -g(3) * t20 + t16 * t19;
t24 = cos(pkin(8));
t17 = t24 * pkin(2) + pkin(1);
t14 = t30 * t17;
t44 = g(2) * t14;
t26 = -pkin(6) - qJ(2);
t35 = pkin(4) * t27 - t26;
t34 = t20 * pkin(3) + t19 * pkin(7);
t15 = g(1) * t28 - g(2) * t30;
t18 = t29 * pkin(4) + pkin(3);
t25 = -qJ(5) - pkin(7);
t32 = t20 * t18 - t19 * t25;
t13 = t20 * t37 + t40;
t11 = -t20 * t39 + t38;
t9 = t15 * t19;
t8 = t16 * t20 + t42;
t6 = t7 * t29;
t5 = t7 * t27;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t29 * t42;
t21 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t24, -t15 * sin(pkin(8)), -t16, -g(1) * (-t28 * pkin(1) + t30 * qJ(2)) - g(2) * (t30 * pkin(1) + t28 * qJ(2)), 0, 0, 0, 0, 0, 0, t15 * t20, -t9, -t16, -g(1) * (-t28 * t17 - t30 * t26) - g(2) * (-t28 * t26 + t14), 0, 0, 0, 0, 0, 0, t4, t3, t9, -t44 + (g(1) * t26 - g(2) * t34) * t30 + (-g(1) * (-t17 - t34) + g(2) * t26) * t28, 0, 0, 0, 0, 0, 0, t4, t3, t9, -t44 + (-g(1) * t35 - g(2) * t32) * t30 + (-g(1) * (-t17 - t32) - g(2) * t35) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t34 + t16 * (pkin(3) * t19 - pkin(7) * t20), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t32 + t16 * (t18 * t19 + t20 * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t21;
