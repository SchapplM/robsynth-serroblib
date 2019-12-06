% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = sin(pkin(9));
t37 = g(1) * t19;
t18 = qJ(1) + pkin(8);
t17 = qJ(3) + t18;
t14 = cos(t17);
t36 = g(2) * t14;
t13 = sin(t17);
t35 = t13 * pkin(3);
t20 = cos(pkin(9));
t21 = sin(qJ(5));
t34 = t20 * t21;
t23 = cos(qJ(5));
t33 = t20 * t23;
t15 = sin(t18);
t22 = sin(qJ(1));
t32 = -t22 * pkin(1) - pkin(2) * t15;
t16 = cos(t18);
t24 = cos(qJ(1));
t31 = -t24 * pkin(1) - pkin(2) * t16;
t10 = g(3) * t13 + t36;
t30 = g(2) * t24 + g(3) * t22;
t29 = -t14 * pkin(3) - t13 * qJ(4);
t28 = pkin(4) * t20 + pkin(7) * t19 + pkin(3);
t11 = t14 * qJ(4);
t27 = t11 + t32;
t26 = g(2) * t31;
t25 = (g(2) * qJ(4) + g(3) * t28) * t13 + t28 * t36;
t9 = g(2) * t13 - g(3) * t14;
t8 = t10 * t20;
t7 = t10 * t19;
t6 = -t13 * t21 - t14 * t33;
t5 = -t13 * t23 + t14 * t34;
t4 = t13 * t33 - t14 * t21;
t3 = t13 * t34 + t14 * t23;
t2 = -g(2) * t6 + g(3) * t4;
t1 = -g(2) * t5 - g(3) * t3;
t12 = [0, 0, 0, 0, 0, 0, t30, -g(2) * t22 + g(3) * t24, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t16 + g(3) * t15, -g(2) * t15 + g(3) * t16, 0, t30 * pkin(1), 0, 0, 0, 0, 0, 0, t10, -t9, 0, -g(3) * t32 - t26, 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * (t29 + t31) - g(3) * (t27 - t35), 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(3) * t27 + t25 - t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * t29 - g(3) * (t11 - t35), 0, 0, 0, 0, 0, 0, t2, t1, t7, -g(3) * t11 + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 + g(3) * t5 + t21 * t37, -g(2) * t4 - g(3) * t6 + t23 * t37, 0, 0;];
taug_reg = t12;
