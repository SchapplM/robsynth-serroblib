% Calculate inertial parameters regressor of gravitation load for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:59
% EndTime: 2019-12-05 15:17:00
% DurationCPUTime: 0.22s
% Computational Cost: add. (144->51), mult. (281->82), div. (0->0), fcn. (328->10), ass. (0->34)
t14 = sin(pkin(9));
t35 = g(3) * t14;
t15 = sin(pkin(8));
t34 = t14 * t15;
t17 = cos(pkin(8));
t33 = t14 * t17;
t18 = sin(qJ(4));
t32 = t14 * t18;
t20 = cos(qJ(4));
t31 = t14 * t20;
t21 = cos(qJ(3));
t30 = t14 * t21;
t19 = sin(qJ(3));
t29 = t15 * t19;
t28 = t15 * t21;
t27 = t17 * t19;
t26 = t17 * t21;
t16 = cos(pkin(9));
t4 = -t16 * t29 - t26;
t6 = -t16 * t27 + t28;
t25 = -g(1) * t6 - g(2) * t4 + t19 * t35;
t5 = t16 * t28 - t27;
t7 = t16 * t26 + t29;
t24 = g(1) * t7 + g(2) * t5 + g(3) * t30;
t23 = -g(1) * (t17 * t31 - t7 * t18) - g(2) * (t15 * t31 - t5 * t18) - g(3) * (-t16 * t20 - t18 * t30);
t22 = -pkin(7) - pkin(6);
t13 = qJ(4) + qJ(5);
t12 = cos(t13);
t11 = sin(t13);
t10 = t20 * pkin(4) + pkin(3);
t8 = -g(1) * t15 + g(2) * t17;
t2 = -g(1) * (-t11 * t33 - t7 * t12) - g(2) * (-t11 * t34 - t5 * t12) - g(3) * (t16 * t11 - t12 * t30);
t1 = -g(1) * (-t7 * t11 + t12 * t33) - g(2) * (-t5 * t11 + t12 * t34) - g(3) * (-t11 * t30 - t16 * t12);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t24, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t20, -t25 * t18, -t24, -g(1) * (t6 * pkin(3) + t7 * pkin(6)) - g(2) * (t4 * pkin(3) + t5 * pkin(6)) - (-pkin(3) * t19 + pkin(6) * t21) * t35, 0, 0, 0, 0, 0, 0, t25 * t12, -t25 * t11, -t24, -g(1) * (t6 * t10 - t7 * t22) - g(2) * (t4 * t10 - t5 * t22) - (-t10 * t19 - t21 * t22) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -g(1) * (-t17 * t32 - t7 * t20) - g(2) * (-t15 * t32 - t5 * t20) - g(3) * (t16 * t18 - t20 * t30), 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t3;
