% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:41
% EndTime: 2019-12-05 15:54:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (170->48), mult. (211->63), div. (0->0), fcn. (211->10), ass. (0->25)
t17 = sin(pkin(8));
t19 = cos(pkin(8));
t28 = g(1) * t19 + g(2) * t17;
t21 = sin(qJ(2));
t22 = cos(qJ(2));
t3 = -g(3) * t22 + t28 * t21;
t32 = g(3) * t21;
t18 = cos(pkin(9));
t8 = t18 * pkin(3) + pkin(2);
t30 = t17 * t22;
t29 = t19 * t22;
t20 = -pkin(6) - qJ(3);
t15 = pkin(9) + qJ(4);
t10 = cos(t15);
t9 = sin(t15);
t23 = -g(1) * (t17 * t10 - t9 * t29) - g(2) * (-t19 * t10 - t9 * t30) + t9 * t32;
t14 = -pkin(7) + t20;
t11 = qJ(5) + t15;
t7 = cos(t11);
t6 = sin(t11);
t5 = pkin(4) * t10 + t8;
t4 = t28 * t22 + t32;
t2 = -g(1) * (-t17 * t6 - t7 * t29) - g(2) * (t19 * t6 - t7 * t30) + t7 * t32;
t1 = -g(1) * (t17 * t7 - t6 * t29) - g(2) * (-t19 * t7 - t6 * t30) + t6 * t32;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t18, -t3 * sin(pkin(9)), -t4, -g(3) * (t22 * pkin(2) + t21 * qJ(3)) + t28 * (pkin(2) * t21 - qJ(3) * t22), 0, 0, 0, 0, 0, 0, t3 * t10, -t3 * t9, -t4, -g(3) * (-t21 * t20 + t22 * t8) + t28 * (t20 * t22 + t21 * t8), 0, 0, 0, 0, 0, 0, t3 * t7, -t3 * t6, -t4, -g(3) * (-t21 * t14 + t22 * t5) + t28 * (t14 * t22 + t21 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -g(1) * (-t10 * t29 - t17 * t9) - g(2) * (-t10 * t30 + t19 * t9) + t10 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t12;
