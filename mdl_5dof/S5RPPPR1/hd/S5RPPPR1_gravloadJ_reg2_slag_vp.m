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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:20:33
% EndTime: 2020-01-03 11:20:34
% DurationCPUTime: 0.19s
% Computational Cost: add. (158->51), mult. (144->74), div. (0->0), fcn. (148->10), ass. (0->34)
t22 = sin(pkin(8));
t39 = g(1) * t22;
t20 = qJ(1) + pkin(7);
t14 = sin(t20);
t21 = sin(pkin(9));
t38 = t14 * t21;
t24 = cos(pkin(8));
t37 = t14 * t24;
t16 = cos(t20);
t36 = t16 * t24;
t35 = t21 * t24;
t34 = t22 * (-pkin(6) - qJ(4));
t23 = cos(pkin(9));
t33 = t23 * t24;
t26 = sin(qJ(1));
t32 = t26 * pkin(1) + t14 * pkin(2);
t27 = cos(qJ(1));
t31 = t27 * pkin(1) + t16 * pkin(2) + t14 * qJ(3);
t30 = -t16 * qJ(3) + t32;
t8 = g(2) * t16 + g(3) * t14;
t7 = g(2) * t14 - g(3) * t16;
t29 = -g(2) * t27 - g(3) * t26;
t28 = pkin(3) * t24 + qJ(4) * t22;
t19 = pkin(9) + qJ(5);
t15 = cos(t19);
t13 = sin(t19);
t12 = t23 * pkin(4) + pkin(3);
t6 = t8 * t22;
t5 = g(1) * t24 - t7 * t22;
t4 = t14 * t13 + t15 * t36;
t3 = t13 * t36 - t14 * t15;
t2 = -t16 * t13 + t15 * t37;
t1 = -t13 * t37 - t16 * t15;
t9 = [0, 0, 0, 0, 0, 0, t29, g(2) * t26 - g(3) * t27, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, t29 * pkin(1), 0, 0, 0, 0, 0, 0, -t8 * t24, t6, -t7, -g(2) * t31 - g(3) * t30, 0, 0, 0, 0, 0, 0, -g(2) * (t16 * t33 + t38) - g(3) * (t14 * t33 - t16 * t21), -g(2) * (t14 * t23 - t16 * t35) - g(3) * (-t14 * t35 - t16 * t23), -t6, -g(2) * (t28 * t16 + t31) - g(3) * (t28 * t14 + t30), 0, 0, 0, 0, 0, 0, -g(2) * t4 - g(3) * t2, g(2) * t3 - g(3) * t1, -t6, -g(2) * (pkin(4) * t38 + t31) - g(3) * (t12 * t37 - t14 * t34 + t32) + (-g(2) * (t12 * t24 - t34) - g(3) * (-pkin(4) * t21 - qJ(3))) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t1 - g(3) * t3 + t13 * t39, g(2) * t2 - g(3) * t4 + t15 * t39, 0, 0;];
taug_reg = t9;
