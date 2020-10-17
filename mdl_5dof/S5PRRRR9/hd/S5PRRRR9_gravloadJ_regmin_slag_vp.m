% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:10
% EndTime: 2019-12-05 17:21:11
% DurationCPUTime: 0.23s
% Computational Cost: add. (209->65), mult. (456->123), div. (0->0), fcn. (582->12), ass. (0->37)
t18 = sin(pkin(5));
t40 = g(3) * t18;
t16 = qJ(4) + qJ(5);
t14 = sin(t16);
t23 = cos(qJ(3));
t39 = t14 * t23;
t15 = cos(t16);
t38 = t15 * t23;
t21 = sin(qJ(2));
t37 = t18 * t21;
t36 = t18 * t23;
t24 = cos(qJ(2));
t35 = t18 * t24;
t19 = sin(qJ(4));
t34 = t19 * t23;
t22 = cos(qJ(4));
t33 = t22 * t23;
t32 = t23 * t24;
t31 = cos(pkin(5));
t30 = cos(pkin(10));
t17 = sin(pkin(10));
t29 = t17 * t31;
t28 = t18 * t30;
t27 = t31 * t30;
t10 = -t21 * t29 + t30 * t24;
t20 = sin(qJ(3));
t8 = t17 * t24 + t21 * t27;
t26 = g(1) * (-t10 * t20 + t17 * t36) + g(2) * (-t8 * t20 - t23 * t28) + g(3) * (-t20 * t37 + t31 * t23);
t7 = t17 * t21 - t24 * t27;
t9 = t30 * t21 + t24 * t29;
t25 = -g(1) * t9 - g(2) * t7 + g(3) * t35;
t12 = t31 * t20 + t21 * t36;
t6 = t17 * t18 * t20 + t10 * t23;
t4 = -t20 * t28 + t8 * t23;
t2 = -g(1) * (-t9 * t14 - t6 * t15) - g(2) * (-t7 * t14 - t4 * t15) - g(3) * (-t12 * t15 + t14 * t35);
t1 = -g(1) * (-t6 * t14 + t9 * t15) - g(2) * (-t4 * t14 + t7 * t15) - g(3) * (-t12 * t14 - t15 * t35);
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t25, g(1) * t10 + g(2) * t8 + g(3) * t37, 0, 0, 0, 0, 0, -t25 * t23, t25 * t20, 0, 0, 0, 0, 0, -g(1) * (t10 * t19 - t9 * t33) - g(2) * (t8 * t19 - t7 * t33) - (t19 * t21 + t22 * t32) * t40, -g(1) * (t10 * t22 + t9 * t34) - g(2) * (t8 * t22 + t7 * t34) - (-t19 * t32 + t21 * t22) * t40, 0, 0, 0, 0, 0, -g(1) * (t10 * t14 - t9 * t38) - g(2) * (t8 * t14 - t7 * t38) - (t14 * t21 + t15 * t32) * t40, -g(1) * (t10 * t15 + t9 * t39) - g(2) * (t8 * t15 + t7 * t39) - (-t14 * t32 + t15 * t21) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, g(1) * t6 + g(2) * t4 + g(3) * t12, 0, 0, 0, 0, 0, -t26 * t22, t26 * t19, 0, 0, 0, 0, 0, -t26 * t15, t26 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t6 * t19 + t9 * t22) - g(2) * (-t4 * t19 + t7 * t22) - g(3) * (-t12 * t19 - t22 * t35), -g(1) * (-t9 * t19 - t6 * t22) - g(2) * (-t7 * t19 - t4 * t22) - g(3) * (-t12 * t22 + t19 * t35), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t3;
