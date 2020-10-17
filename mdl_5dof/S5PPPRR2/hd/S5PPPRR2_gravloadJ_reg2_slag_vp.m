% Calculate inertial parameters regressor of gravitation load for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:38
% EndTime: 2019-12-05 14:59:38
% DurationCPUTime: 0.14s
% Computational Cost: add. (89->40), mult. (220->65), div. (0->0), fcn. (272->10), ass. (0->31)
t14 = sin(pkin(9));
t15 = sin(pkin(8));
t31 = t14 * t15;
t21 = sin(qJ(4));
t30 = t15 * t21;
t23 = cos(qJ(4));
t29 = t15 * t23;
t16 = sin(pkin(7));
t18 = cos(pkin(8));
t28 = t16 * t18;
t19 = cos(pkin(7));
t27 = t19 * t14;
t17 = cos(pkin(9));
t26 = t19 * t17;
t7 = t17 * t28 - t27;
t1 = t16 * t29 - t7 * t21;
t10 = -t17 * t30 - t18 * t23;
t9 = t16 * t14 + t18 * t26;
t3 = t19 * t29 - t9 * t21;
t25 = g(1) * t3 + g(2) * t1 + g(3) * t10;
t11 = t17 * t29 - t18 * t21;
t2 = t16 * t30 + t7 * t23;
t4 = t19 * t30 + t9 * t23;
t24 = g(1) * t4 + g(2) * t2 + g(3) * t11;
t22 = cos(qJ(5));
t20 = sin(qJ(5));
t12 = -g(1) * t16 + g(2) * t19;
t8 = -t16 * t17 + t18 * t27;
t6 = t14 * t28 + t26;
t5 = g(3) * t18 + (-g(1) * t19 - g(2) * t16) * t15;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t22, t25 * t20, -t24, -g(1) * (t3 * pkin(4) + t4 * pkin(6)) - g(2) * (t1 * pkin(4) + t2 * pkin(6)) - g(3) * (t10 * pkin(4) + t11 * pkin(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t20 + t8 * t22) - g(2) * (-t2 * t20 + t6 * t22) - g(3) * (-t11 * t20 + t22 * t31), -g(1) * (-t8 * t20 - t4 * t22) - g(2) * (-t2 * t22 - t6 * t20) - g(3) * (-t11 * t22 - t20 * t31), 0, 0;];
taug_reg = t13;
