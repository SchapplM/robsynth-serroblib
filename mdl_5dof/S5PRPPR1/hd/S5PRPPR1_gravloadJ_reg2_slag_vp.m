% Calculate inertial parameters regressor of gravitation load for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:10
% EndTime: 2019-12-05 15:22:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (150->49), mult. (130->67), div. (0->0), fcn. (136->8), ass. (0->31)
t19 = pkin(7) + qJ(2);
t15 = sin(t19);
t34 = g(1) * t15;
t21 = sin(pkin(8));
t33 = g(3) * t21;
t17 = cos(t19);
t32 = t17 * pkin(2) + t15 * qJ(3);
t23 = cos(pkin(8));
t31 = t15 * t23;
t20 = sin(pkin(9));
t30 = t17 * t20;
t29 = t17 * t23;
t28 = t20 * t23;
t27 = t21 * (-pkin(6) - qJ(4));
t22 = cos(pkin(9));
t26 = t22 * t23;
t8 = g(1) * t17 + g(2) * t15;
t7 = -g(2) * t17 + t34;
t25 = pkin(3) * t23 + qJ(4) * t21;
t18 = pkin(9) + qJ(5);
t16 = cos(t18);
t14 = sin(t18);
t13 = t22 * pkin(4) + pkin(3);
t10 = t17 * qJ(3);
t6 = t7 * t21;
t5 = g(3) * t23 - t8 * t21;
t4 = t15 * t14 + t16 * t29;
t3 = -t14 * t29 + t15 * t16;
t2 = t17 * t14 - t16 * t31;
t1 = t14 * t31 + t17 * t16;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t23, -t6, -t8, -g(1) * (-t15 * pkin(2) + t10) - g(2) * t32, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t26 + t30) - g(2) * (t15 * t20 + t17 * t26), -g(1) * (t15 * t28 + t17 * t22) - g(2) * (t15 * t22 - t17 * t28), t6, -g(1) * t10 - g(2) * (t25 * t17 + t32) - (-pkin(2) - t25) * t34, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t6, -g(1) * (pkin(4) * t30 + t10) - g(2) * (t13 * t29 - t17 * t27 + t32) + (-g(1) * (-t13 * t23 - pkin(2) + t27) - g(2) * pkin(4) * t20) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t14 * t33, g(1) * t4 - g(2) * t2 + t16 * t33, 0, 0;];
taug_reg = t9;
