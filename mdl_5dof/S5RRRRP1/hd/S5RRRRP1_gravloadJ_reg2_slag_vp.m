% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = -pkin(7) - pkin(6);
t25 = qJ(2) + qJ(3);
t19 = sin(t25);
t36 = pkin(3) * t19;
t26 = sin(qJ(2));
t35 = t26 * pkin(2);
t20 = cos(t25);
t15 = pkin(3) * t20;
t28 = cos(qJ(2));
t23 = t28 * pkin(2);
t34 = t15 + t23;
t24 = -pkin(8) + t30;
t22 = qJ(4) + t25;
t17 = cos(t22);
t14 = pkin(4) * t17;
t33 = t14 + t34;
t16 = sin(t22);
t32 = -pkin(4) * t16 - t36;
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t13 = g(1) * t29 + g(2) * t27;
t12 = g(1) * t27 - g(2) * t29;
t1 = -g(3) * t17 + t13 * t16;
t3 = -g(3) * t20 + t13 * t19;
t31 = -g(3) * t28 + t13 * t26;
t21 = -qJ(5) + t24;
t18 = t23 + pkin(1);
t10 = pkin(1) + t34;
t7 = pkin(1) + t33;
t6 = t12 * t17;
t5 = t12 * t16;
t4 = g(3) * t19 + t13 * t20;
t2 = g(3) * t16 + t13 * t17;
t8 = [0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t28, -t12 * t26, -t13, -g(1) * (-t27 * pkin(1) + t29 * pkin(6)) - g(2) * (t29 * pkin(1) + t27 * pkin(6)), 0, 0, 0, 0, 0, 0, t12 * t20, -t12 * t19, -t13, -g(1) * (-t27 * t18 - t29 * t30) - g(2) * (t29 * t18 - t27 * t30), 0, 0, 0, 0, 0, 0, t6, -t5, -t13, -g(1) * (-t27 * t10 - t29 * t24) - g(2) * (t29 * t10 - t27 * t24), 0, 0, 0, 0, 0, 0, t6, -t5, -t13, -g(1) * (-t29 * t21 - t27 * t7) - g(2) * (-t27 * t21 + t29 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, g(3) * t26 + t13 * t28, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t31 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t34 - t13 * (-t35 - t36), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t33 - t13 * (t32 - t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * (t14 + t15) - t13 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12;];
taug_reg = t8;
