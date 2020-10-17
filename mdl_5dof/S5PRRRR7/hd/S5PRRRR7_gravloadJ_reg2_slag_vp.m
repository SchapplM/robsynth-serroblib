% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:00
% DurationCPUTime: 0.29s
% Computational Cost: add. (236->59), mult. (291->89), div. (0->0), fcn. (297->10), ass. (0->33)
t19 = sin(pkin(9));
t20 = cos(pkin(9));
t31 = g(1) * t20 + g(2) * t19;
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t28 = -g(3) * t24 + t31 * t22;
t25 = -pkin(7) - pkin(6);
t38 = g(3) * t22;
t18 = qJ(3) + qJ(4);
t14 = cos(t18);
t23 = cos(qJ(3));
t16 = t23 * pkin(3);
t8 = pkin(4) * t14 + t16;
t36 = t19 * t24;
t35 = t20 * t24;
t21 = sin(qJ(3));
t34 = t21 * t24;
t33 = t23 * t24;
t13 = sin(t18);
t3 = -g(1) * (-t13 * t35 + t19 * t14) - g(2) * (-t13 * t36 - t20 * t14) + t13 * t38;
t26 = -g(1) * (t19 * t23 - t20 * t34) - g(2) * (-t19 * t34 - t20 * t23) + t21 * t38;
t17 = -pkin(8) + t25;
t15 = qJ(5) + t18;
t12 = t16 + pkin(2);
t11 = cos(t15);
t10 = sin(t15);
t7 = -t21 * pkin(3) - pkin(4) * t13;
t6 = pkin(2) + t8;
t5 = t31 * t24 + t38;
t4 = -g(1) * (-t19 * t13 - t14 * t35) - g(2) * (t20 * t13 - t14 * t36) + t14 * t38;
t2 = -g(1) * (-t19 * t10 - t11 * t35) - g(2) * (t20 * t10 - t11 * t36) + t11 * t38;
t1 = -g(1) * (-t10 * t35 + t19 * t11) - g(2) * (-t10 * t36 - t20 * t11) + t10 * t38;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t5, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t23, -t28 * t21, -t5, -g(3) * (t24 * pkin(2) + t22 * pkin(6)) + t31 * (pkin(2) * t22 - pkin(6) * t24), 0, 0, 0, 0, 0, 0, t28 * t14, -t28 * t13, -t5, -g(3) * (t24 * t12 - t22 * t25) + t31 * (t12 * t22 + t24 * t25), 0, 0, 0, 0, 0, 0, t28 * t11, -t28 * t10, -t5, -g(3) * (-t22 * t17 + t24 * t6) + t31 * (t17 * t24 + t22 * t6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -g(1) * (-t19 * t21 - t20 * t33) - g(2) * (-t19 * t33 + t20 * t21) + t23 * t38, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t26 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t19 * t8 + t35 * t7) - g(2) * (-t20 * t8 + t36 * t7) - t7 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
