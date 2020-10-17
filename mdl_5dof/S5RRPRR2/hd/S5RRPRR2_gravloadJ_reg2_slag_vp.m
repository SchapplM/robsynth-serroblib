% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:27
% EndTime: 2019-12-05 18:28:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (232->48), mult. (197->64), div. (0->0), fcn. (177->10), ass. (0->31)
t27 = -qJ(3) - pkin(6);
t26 = qJ(2) + pkin(9);
t20 = cos(t26);
t30 = cos(qJ(2));
t23 = t30 * pkin(2);
t35 = pkin(3) * t20 + t23;
t25 = -pkin(7) + t27;
t21 = qJ(4) + t26;
t16 = cos(t21);
t34 = pkin(4) * t16 + t35;
t19 = sin(t26);
t28 = sin(qJ(2));
t33 = -t28 * pkin(2) - pkin(3) * t19;
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t10 = g(1) * t31 + g(2) * t29;
t9 = g(1) * t29 - g(2) * t31;
t15 = sin(t21);
t3 = -g(3) * t16 + t10 * t15;
t32 = -g(3) * t30 + t10 * t28;
t22 = -pkin(8) + t25;
t18 = t23 + pkin(1);
t17 = qJ(5) + t21;
t13 = cos(t17);
t12 = sin(t17);
t7 = pkin(1) + t35;
t5 = pkin(1) + t34;
t4 = g(3) * t15 + t10 * t16;
t2 = g(3) * t12 + t10 * t13;
t1 = -g(3) * t13 + t10 * t12;
t6 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t30, -t9 * t28, -t10, -g(1) * (-t29 * pkin(1) + t31 * pkin(6)) - g(2) * (t31 * pkin(1) + t29 * pkin(6)), 0, 0, 0, 0, 0, 0, t9 * t20, -t9 * t19, -t10, -g(1) * (-t29 * t18 - t31 * t27) - g(2) * (t31 * t18 - t29 * t27), 0, 0, 0, 0, 0, 0, t9 * t16, -t9 * t15, -t10, -g(1) * (-t31 * t25 - t29 * t7) - g(2) * (-t29 * t25 + t31 * t7), 0, 0, 0, 0, 0, 0, t9 * t13, -t9 * t12, -t10, -g(1) * (-t31 * t22 - t29 * t5) - g(2) * (-t29 * t22 + t31 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t28 + t10 * t30, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t20 + t10 * t19, g(3) * t19 + t10 * t20, 0, t32 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t35 - t10 * t33, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t34 - t10 * (-pkin(4) * t15 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t6;
