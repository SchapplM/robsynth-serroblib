% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR2
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:14
% DurationCPUTime: 0.16s
% Computational Cost: add. (175->41), mult. (161->55), div. (0->0), fcn. (153->10), ass. (0->33)
t18 = qJ(2) + pkin(9);
t16 = qJ(4) + t18;
t12 = sin(t16);
t13 = cos(t16);
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t26 = g(1) * t20 + g(2) * t19;
t3 = -g(3) * t13 + t26 * t12;
t37 = pkin(4) * t12;
t36 = pkin(7) * t13;
t35 = g(3) * t12;
t33 = t13 * pkin(4) + t12 * pkin(7);
t21 = sin(qJ(5));
t32 = t19 * t21;
t23 = cos(qJ(5));
t31 = t19 * t23;
t30 = t20 * t21;
t29 = t20 * t23;
t15 = cos(t18);
t24 = cos(qJ(2));
t28 = t24 * pkin(2) + pkin(3) * t15;
t14 = sin(t18);
t22 = sin(qJ(2));
t5 = -t22 * pkin(2) - pkin(3) * t14;
t27 = t5 - t37;
t25 = -g(3) * t24 + t26 * t22;
t8 = -g(1) * t19 + g(2) * t20;
t7 = t20 * t36;
t6 = t19 * t36;
t4 = t26 * t13 + t35;
t2 = t3 * t23;
t1 = t3 * t21;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t22 + t26 * t24, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t15 + t26 * t14, g(3) * t14 + t26 * t15, 0, t25 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t28 - t26 * t5, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t27 * t20 + t7) - g(2) * (t27 * t19 + t6) - g(3) * (t28 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t20 * t37 + t7) - g(2) * (-t19 * t37 + t6) - g(3) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t30 + t31) - g(2) * (-t13 * t32 - t29) + t21 * t35, -g(1) * (-t13 * t29 - t32) - g(2) * (-t13 * t31 + t30) + t23 * t35, 0, 0;];
taug_reg = t9;
