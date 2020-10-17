% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:07
% EndTime: 2019-12-31 17:58:08
% DurationCPUTime: 0.18s
% Computational Cost: add. (186->51), mult. (158->64), div. (0->0), fcn. (155->10), ass. (0->34)
t17 = pkin(9) + qJ(4);
t12 = sin(t17);
t14 = cos(t17);
t30 = t14 * pkin(4) + t12 * pkin(7);
t18 = qJ(1) + pkin(8);
t13 = sin(t18);
t15 = cos(t18);
t8 = g(1) * t15 + g(2) * t13;
t26 = -g(3) * t14 + t8 * t12;
t40 = g(3) * t12;
t23 = sin(qJ(1));
t36 = t23 * pkin(1);
t20 = cos(pkin(9));
t11 = t20 * pkin(3) + pkin(2);
t25 = cos(qJ(1));
t16 = t25 * pkin(1);
t35 = t15 * t11 + t16;
t22 = sin(qJ(5));
t34 = t13 * t22;
t24 = cos(qJ(5));
t33 = t13 * t24;
t32 = t15 * t22;
t31 = t15 * t24;
t7 = g(1) * t13 - g(2) * t15;
t28 = g(1) * t23 - g(2) * t25;
t21 = -pkin(6) - qJ(3);
t27 = -t15 * t21 - t36;
t6 = t14 * t31 + t34;
t5 = -t14 * t32 + t33;
t4 = -t14 * t33 + t32;
t3 = t14 * t34 + t31;
t2 = t7 * t12;
t1 = t8 * t14 + t40;
t9 = [0, 0, 0, 0, 0, 0, t28, g(1) * t25 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, t7 * t20, -t7 * sin(pkin(9)), -t8, -g(1) * (-t13 * pkin(2) + t15 * qJ(3) - t36) - g(2) * (t15 * pkin(2) + t13 * qJ(3) + t16), 0, 0, 0, 0, 0, 0, t7 * t14, -t2, -t8, -g(1) * (-t13 * t11 + t27) - g(2) * (-t13 * t21 + t35), 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t2, -g(1) * t27 - g(2) * (t30 * t15 + t35) + (-g(1) * (-t11 - t30) + g(2) * t21) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t1, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t24, -t26 * t22, -t1, -g(3) * t30 + t8 * (pkin(4) * t12 - pkin(7) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t22 * t40, g(1) * t6 - g(2) * t4 + t24 * t40, 0, 0;];
taug_reg = t9;
