% Calculate inertial parameters regressor of gravitation load for
% S5RPPPR1
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:44
% EndTime: 2022-01-20 09:12:44
% DurationCPUTime: 0.20s
% Computational Cost: add. (158->52), mult. (144->74), div. (0->0), fcn. (148->10), ass. (0->34)
t20 = qJ(1) + pkin(7);
t15 = sin(t20);
t39 = g(1) * t15;
t22 = sin(pkin(8));
t38 = g(3) * t22;
t24 = cos(pkin(8));
t37 = t15 * t24;
t17 = cos(t20);
t21 = sin(pkin(9));
t36 = t17 * t21;
t35 = t17 * t24;
t34 = t21 * t24;
t33 = t22 * (-pkin(6) - qJ(4));
t23 = cos(pkin(9));
t32 = t23 * t24;
t27 = cos(qJ(1));
t31 = t27 * pkin(1) + t17 * pkin(2) + t15 * qJ(3);
t26 = sin(qJ(1));
t30 = -t26 * pkin(1) + t17 * qJ(3);
t8 = g(1) * t17 + g(2) * t15;
t7 = -g(2) * t17 + t39;
t29 = g(1) * t26 - g(2) * t27;
t28 = pkin(3) * t24 + qJ(4) * t22;
t19 = pkin(9) + qJ(5);
t16 = cos(t19);
t14 = sin(t19);
t13 = t23 * pkin(4) + pkin(3);
t6 = t7 * t22;
t5 = g(3) * t24 - t8 * t22;
t4 = t15 * t14 + t16 * t35;
t3 = -t14 * t35 + t15 * t16;
t2 = t17 * t14 - t16 * t37;
t1 = t14 * t37 + t17 * t16;
t9 = [0, 0, 0, 0, 0, 0, t29, g(1) * t27 + g(2) * t26, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t29 * pkin(1), 0, 0, 0, 0, 0, 0, t7 * t24, -t6, -t8, -g(1) * (-t15 * pkin(2) + t30) - g(2) * t31, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t32 + t36) - g(2) * (t15 * t21 + t17 * t32), -g(1) * (t15 * t34 + t17 * t23) - g(2) * (t15 * t23 - t17 * t34), t6, -g(1) * t30 - g(2) * (t28 * t17 + t31) - (-pkin(2) - t28) * t39, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t6, -g(1) * (pkin(4) * t36 + t30) - g(2) * (t13 * t35 - t17 * t33 + t31) + (-g(1) * (-t13 * t24 - pkin(2) + t33) - g(2) * pkin(4) * t21) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t14 * t38, g(1) * t4 - g(2) * t2 + t16 * t38, 0, 0;];
taug_reg = t9;
