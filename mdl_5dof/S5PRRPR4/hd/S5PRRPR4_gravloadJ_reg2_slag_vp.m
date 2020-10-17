% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:22
% EndTime: 2019-12-05 16:23:24
% DurationCPUTime: 0.27s
% Computational Cost: add. (185->61), mult. (246->88), div. (0->0), fcn. (247->10), ass. (0->31)
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t30 = g(1) * t20 + g(2) * t19;
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t3 = -g(3) * t25 + t30 * t23;
t37 = g(3) * t23;
t18 = qJ(3) + pkin(9);
t13 = cos(t18);
t24 = cos(qJ(3));
t15 = t24 * pkin(3);
t7 = pkin(4) * t13 + t15;
t35 = t19 * t25;
t34 = t20 * t25;
t22 = sin(qJ(3));
t33 = t22 * t25;
t32 = t24 * t25;
t21 = -qJ(4) - pkin(6);
t26 = -g(1) * (t19 * t24 - t20 * t33) - g(2) * (-t19 * t33 - t20 * t24) + t22 * t37;
t17 = -pkin(7) + t21;
t14 = qJ(5) + t18;
t12 = sin(t18);
t11 = t15 + pkin(2);
t10 = cos(t14);
t9 = sin(t14);
t6 = -t22 * pkin(3) - pkin(4) * t12;
t5 = pkin(2) + t7;
t4 = t30 * t25 + t37;
t2 = -g(1) * (-t10 * t34 - t19 * t9) - g(2) * (-t10 * t35 + t20 * t9) + t10 * t37;
t1 = -g(1) * (t19 * t10 - t9 * t34) - g(2) * (-t20 * t10 - t9 * t35) + t9 * t37;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t24, -t3 * t22, -t4, -g(3) * (t25 * pkin(2) + t23 * pkin(6)) + t30 * (pkin(2) * t23 - pkin(6) * t25), 0, 0, 0, 0, 0, 0, t3 * t13, -t3 * t12, -t4, -g(3) * (t25 * t11 - t23 * t21) + t30 * (t11 * t23 + t21 * t25), 0, 0, 0, 0, 0, 0, t3 * t10, -t3 * t9, -t4, -g(3) * (-t23 * t17 + t25 * t5) + t30 * (t17 * t25 + t23 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -g(1) * (-t19 * t22 - t20 * t32) - g(2) * (-t19 * t32 + t20 * t22) + t24 * t37, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t34 + t19 * t13) - g(2) * (-t12 * t35 - t20 * t13) + t12 * t37, -g(1) * (-t19 * t12 - t13 * t34) - g(2) * (t20 * t12 - t13 * t35) + t13 * t37, 0, t26 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t19 * t7 + t6 * t34) - g(2) * (-t20 * t7 + t6 * t35) - t6 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t8;
