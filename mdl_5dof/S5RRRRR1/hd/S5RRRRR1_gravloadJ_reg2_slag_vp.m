% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:51:00
% EndTime: 2019-12-05 18:51:01
% DurationCPUTime: 0.24s
% Computational Cost: add. (317->52), mult. (310->75), div. (0->0), fcn. (292->10), ass. (0->46)
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t35 = g(1) * t28 - g(2) * t31;
t16 = g(1) * t31 + g(2) * t28;
t25 = qJ(2) + qJ(3);
t23 = qJ(4) + t25;
t18 = sin(t23);
t19 = cos(t23);
t3 = g(3) * t19 + t16 * t18;
t21 = sin(t25);
t51 = pkin(3) * t21;
t22 = cos(t25);
t50 = pkin(3) * t22;
t49 = pkin(4) * t18;
t48 = pkin(6) * t19;
t45 = g(3) * t18;
t26 = sin(qJ(5));
t43 = t28 * t26;
t29 = cos(qJ(5));
t42 = t28 * t29;
t41 = t31 * t26;
t40 = t31 * t29;
t30 = cos(qJ(2));
t24 = t30 * pkin(2);
t39 = t24 + t50;
t27 = sin(qJ(2));
t13 = -t27 * pkin(2) - t51;
t38 = t13 - t49;
t37 = -t49 - t51;
t36 = t19 * pkin(4) + t18 * pkin(6);
t33 = -t36 - t50;
t5 = g(3) * t22 + t16 * t21;
t32 = g(3) * t30 + t16 * t27;
t15 = t31 * t48;
t14 = t28 * t48;
t12 = pkin(1) + t39;
t11 = t19 * t40 - t43;
t10 = -t19 * t41 - t42;
t9 = -t19 * t42 - t41;
t8 = t19 * t43 - t40;
t7 = t35 * t18;
t6 = -g(3) * t21 + t16 * t22;
t4 = t16 * t19 - t45;
t2 = t3 * t29;
t1 = t3 * t26;
t17 = [0, 0, 0, 0, 0, 0, t35, t16, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t30, -t35 * t27, t16, t35 * pkin(1), 0, 0, 0, 0, 0, 0, t35 * t22, -t35 * t21, t16, t35 * (t24 + pkin(1)), 0, 0, 0, 0, 0, 0, t35 * t19, -t7, t16, t35 * t12, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, t35 * (t12 + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -g(3) * t27 + t16 * t30, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t32 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, g(3) * t39 - t16 * t13, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t38 * t31 + t15) - g(2) * (t38 * t28 + t14) - g(3) * (-t24 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t37 * t31 + t15) - g(2) * (t37 * t28 + t14) - g(3) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t31 * t49 + t15) - g(2) * (-t28 * t49 + t14) + g(3) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 - t26 * t45, g(1) * t11 - g(2) * t9 - t29 * t45, 0, 0;];
taug_reg = t17;
