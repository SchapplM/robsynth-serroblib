% Calculate minimal parameter regressor of gravitation load for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:53
% EndTime: 2019-12-05 15:57:54
% DurationCPUTime: 0.16s
% Computational Cost: add. (161->52), mult. (307->93), div. (0->0), fcn. (378->12), ass. (0->32)
t18 = sin(pkin(5));
t38 = g(3) * t18;
t15 = pkin(10) + qJ(4);
t14 = cos(t15);
t20 = sin(qJ(5));
t37 = t14 * t20;
t22 = cos(qJ(5));
t36 = t14 * t22;
t17 = sin(pkin(9));
t35 = t17 * t18;
t21 = sin(qJ(2));
t34 = t18 * t21;
t23 = cos(qJ(2));
t33 = t20 * t23;
t32 = t22 * t23;
t31 = cos(pkin(5));
t30 = cos(pkin(9));
t29 = t17 * t31;
t28 = t18 * t30;
t27 = t31 * t30;
t10 = -t21 * t29 + t30 * t23;
t13 = sin(t15);
t8 = t17 * t23 + t21 * t27;
t26 = g(1) * (-t10 * t13 + t14 * t35) + g(2) * (-t8 * t13 - t14 * t28) + g(3) * (-t13 * t34 + t31 * t14);
t7 = t17 * t21 - t23 * t27;
t9 = t30 * t21 + t23 * t29;
t25 = -g(1) * t9 - g(2) * t7 + t23 * t38;
t24 = g(1) * t10 + g(2) * t8 + g(3) * t34;
t6 = t31 * t13 + t14 * t34;
t4 = t10 * t14 + t13 * t35;
t2 = -t13 * t28 + t8 * t14;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t25, t24, -t25 * cos(pkin(10)), t25 * sin(pkin(10)), -t24, -g(1) * (-t9 * pkin(2) + t10 * qJ(3)) - g(2) * (-t7 * pkin(2) + t8 * qJ(3)) - (pkin(2) * t23 + qJ(3) * t21) * t38, 0, 0, 0, 0, 0, -t25 * t14, t25 * t13, 0, 0, 0, 0, 0, -g(1) * (t10 * t20 - t9 * t36) - g(2) * (t8 * t20 - t7 * t36) - (t14 * t32 + t20 * t21) * t38, -g(1) * (t10 * t22 + t9 * t37) - g(2) * (t8 * t22 + t7 * t37) - (-t14 * t33 + t21 * t22) * t38; 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, g(1) * t4 + g(2) * t2 + g(3) * t6, 0, 0, 0, 0, 0, -t26 * t22, t26 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t20 + t9 * t22) - g(2) * (-t2 * t20 + t7 * t22) - g(3) * (-t18 * t32 - t6 * t20), -g(1) * (-t9 * t20 - t4 * t22) - g(2) * (-t2 * t22 - t7 * t20) - g(3) * (t18 * t33 - t6 * t22);];
taug_reg = t1;
