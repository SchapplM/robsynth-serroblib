% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:19
% DurationCPUTime: 0.24s
% Computational Cost: add. (184->51), mult. (276->73), div. (0->0), fcn. (278->8), ass. (0->30)
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t28 = g(1) * t18 + g(2) * t17;
t20 = sin(qJ(2));
t22 = cos(qJ(2));
t5 = -g(3) * t22 + t28 * t20;
t23 = -pkin(7) - pkin(6);
t35 = g(3) * t20;
t33 = t17 * t22;
t32 = t18 * t22;
t19 = sin(qJ(3));
t31 = t19 * t22;
t21 = cos(qJ(3));
t30 = t21 * t22;
t16 = qJ(3) + qJ(4);
t13 = cos(t16);
t14 = t21 * pkin(3);
t9 = pkin(4) * t13 + t14;
t12 = sin(t16);
t1 = -g(1) * (-t12 * t32 + t17 * t13) - g(2) * (-t12 * t33 - t18 * t13) + t12 * t35;
t24 = -g(1) * (t17 * t21 - t18 * t31) - g(2) * (-t17 * t31 - t18 * t21) + t19 * t35;
t15 = -qJ(5) + t23;
t11 = t14 + pkin(2);
t8 = -t19 * pkin(3) - pkin(4) * t12;
t7 = pkin(2) + t9;
t6 = t28 * t22 + t35;
t4 = t5 * t13;
t3 = t5 * t12;
t2 = -g(1) * (-t17 * t12 - t13 * t32) - g(2) * (t18 * t12 - t13 * t33) + t13 * t35;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t21, -t5 * t19, -t6, -g(3) * (t22 * pkin(2) + t20 * pkin(6)) + t28 * (pkin(2) * t20 - pkin(6) * t22), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (t22 * t11 - t20 * t23) + t28 * (t11 * t20 + t22 * t23), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (-t20 * t15 + t22 * t7) + t28 * (t15 * t22 + t20 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -g(1) * (-t17 * t19 - t18 * t30) - g(2) * (-t17 * t30 + t18 * t19) + t21 * t35, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t24 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t17 * t9 + t8 * t32) - g(2) * (-t18 * t9 + t8 * t33) - t8 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t10;
