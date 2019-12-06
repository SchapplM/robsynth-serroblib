% Calculate inertial parameters regressor of gravitation load for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t11 = sin(pkin(8));
t23 = g(3) * t11;
t15 = sin(qJ(5));
t22 = t11 * t15;
t16 = cos(qJ(5));
t21 = t11 * t16;
t12 = sin(pkin(7));
t13 = cos(pkin(8));
t20 = t12 * t13;
t14 = cos(pkin(7));
t19 = t13 * t14;
t10 = pkin(9) + qJ(4);
t8 = sin(t10);
t9 = cos(t10);
t1 = -t14 * t9 - t8 * t20;
t3 = t12 * t9 - t8 * t19;
t18 = -g(1) * t3 - g(2) * t1 + t8 * t23;
t2 = -t14 * t8 + t9 * t20;
t4 = t12 * t8 + t9 * t19;
t17 = g(1) * t4 + g(2) * t2 + t9 * t23;
t6 = -g(1) * t12 + g(2) * t14;
t5 = g(3) * t13 + (-g(1) * t14 - g(2) * t12) * t11;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t16, -t18 * t15, -t17, -g(1) * (t3 * pkin(4) + t4 * pkin(6)) - g(2) * (t1 * pkin(4) + t2 * pkin(6)) - (-pkin(4) * t8 + pkin(6) * t9) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t21 - t4 * t15) - g(2) * (t12 * t21 - t2 * t15) - g(3) * (-t13 * t16 - t9 * t22), -g(1) * (-t14 * t22 - t4 * t16) - g(2) * (-t12 * t22 - t2 * t16) - g(3) * (t13 * t15 - t9 * t21), 0, 0;];
taug_reg = t7;
