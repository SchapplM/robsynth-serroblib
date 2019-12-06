% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR2
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = qJ(1) + pkin(8);
t17 = qJ(3) + t19;
t10 = sin(t17);
t11 = cos(t17);
t31 = -t10 * pkin(3) + t11 * qJ(4);
t21 = cos(pkin(9));
t12 = t21 * pkin(4) + pkin(3);
t22 = -pkin(7) - qJ(4);
t30 = t10 * t22 - t11 * t12;
t14 = sin(t19);
t23 = sin(qJ(1));
t29 = -t23 * pkin(1) - pkin(2) * t14;
t16 = cos(t19);
t24 = cos(qJ(1));
t28 = -t24 * pkin(1) - pkin(2) * t16;
t6 = g(2) * t11 + g(3) * t10;
t5 = g(2) * t10 - g(3) * t11;
t27 = g(2) * t24 + g(3) * t23;
t26 = -t11 * pkin(3) - t10 * qJ(4);
t25 = -t10 * t12 - t11 * t22;
t18 = pkin(9) + qJ(5);
t15 = cos(t18);
t13 = sin(t18);
t4 = t6 * t21;
t3 = t6 * sin(pkin(9));
t2 = t6 * t15;
t1 = t6 * t13;
t7 = [0, 0, 0, 0, 0, 0, t27, -g(2) * t23 + g(3) * t24, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t16 + g(3) * t14, -g(2) * t14 + g(3) * t16, 0, t27 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, 0, -g(2) * t28 - g(3) * t29, 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * (t26 + t28) - g(3) * (t29 + t31), 0, 0, 0, 0, 0, 0, t2, -t1, t5, -g(2) * (t28 + t30) - g(3) * (t25 + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, t5, -g(2) * t26 - g(3) * t31, 0, 0, 0, 0, 0, 0, t2, -t1, t5, -g(2) * t30 - g(3) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 - t5 * t13, g(1) * t13 - t5 * t15, 0, 0;];
taug_reg = t7;
