% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR4
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = sin(qJ(1));
t29 = t22 * pkin(1);
t24 = cos(qJ(1));
t28 = t24 * pkin(1);
t18 = qJ(3) + pkin(9);
t13 = cos(t18);
t23 = cos(qJ(3));
t16 = t23 * pkin(3);
t27 = pkin(4) * t13 + t16;
t20 = -qJ(4) - pkin(6);
t19 = qJ(1) + pkin(8);
t12 = sin(t19);
t14 = cos(t19);
t4 = g(2) * t14 + g(3) * t12;
t3 = g(2) * t12 - g(3) * t14;
t26 = g(2) * t24 + g(3) * t22;
t21 = sin(qJ(3));
t25 = -g(1) * t23 - t3 * t21;
t17 = -pkin(7) + t20;
t15 = qJ(5) + t18;
t11 = sin(t18);
t10 = t16 + pkin(2);
t9 = cos(t15);
t8 = sin(t15);
t5 = pkin(2) + t27;
t2 = -g(1) * t9 - t3 * t8;
t1 = g(1) * t8 - t3 * t9;
t6 = [0, 0, 0, 0, 0, 0, t26, -g(2) * t22 + g(3) * t24, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, t4 * t23, -t4 * t21, t3, -g(2) * (-t14 * pkin(2) - t12 * pkin(6) - t28) - g(3) * (-t12 * pkin(2) + t14 * pkin(6) - t29), 0, 0, 0, 0, 0, 0, t4 * t13, -t4 * t11, t3, -g(2) * (-t14 * t10 + t12 * t20 - t28) - g(3) * (-t12 * t10 - t14 * t20 - t29), 0, 0, 0, 0, 0, 0, t4 * t9, -t4 * t8, t3, -g(2) * (t12 * t17 - t14 * t5 - t28) - g(3) * (-t12 * t5 - t14 * t17 - t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(1) * t21 - t3 * t23, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - t3 * t11, g(1) * t11 - t3 * t13, 0, t25 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t27 + t3 * (-t21 * pkin(3) - pkin(4) * t11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t6;
