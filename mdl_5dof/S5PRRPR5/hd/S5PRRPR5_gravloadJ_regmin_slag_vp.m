% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:27:55
% DurationCPUTime: 0.20s
% Computational Cost: add. (148->60), mult. (308->110), div. (0->0), fcn. (374->12), ass. (0->39)
t18 = sin(pkin(5));
t43 = g(3) * t18;
t16 = qJ(3) + pkin(10);
t15 = cos(t16);
t22 = sin(qJ(5));
t42 = t15 * t22;
t25 = cos(qJ(5));
t41 = t15 * t25;
t17 = sin(pkin(9));
t40 = t17 * t18;
t19 = cos(pkin(9));
t39 = t18 * t19;
t23 = sin(qJ(3));
t38 = t18 * t23;
t24 = sin(qJ(2));
t37 = t18 * t24;
t26 = cos(qJ(3));
t36 = t18 * t26;
t20 = cos(pkin(5));
t35 = t20 * t24;
t27 = cos(qJ(2));
t34 = t20 * t27;
t33 = t22 * t27;
t32 = t25 * t27;
t10 = -t17 * t35 + t19 * t27;
t14 = sin(t16);
t8 = t17 * t27 + t19 * t35;
t31 = g(1) * (-t10 * t14 + t15 * t40) + g(2) * (-t8 * t14 - t15 * t39) + g(3) * (-t14 * t37 + t20 * t15);
t7 = t17 * t24 - t19 * t34;
t9 = t17 * t34 + t19 * t24;
t30 = -g(1) * t9 - g(2) * t7 + t27 * t43;
t29 = g(1) * t10 + g(2) * t8 + g(3) * t37;
t28 = -g(1) * (-t10 * t23 + t17 * t36) - g(2) * (-t19 * t36 - t8 * t23) - g(3) * (t20 * t26 - t23 * t37);
t21 = -qJ(4) - pkin(7);
t13 = t26 * pkin(3) + pkin(2);
t6 = t20 * t14 + t15 * t37;
t4 = t10 * t15 + t14 * t40;
t2 = -t14 * t39 + t8 * t15;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t30, t29, 0, 0, 0, 0, 0, -t30 * t26, t30 * t23, -t29, -g(1) * (-t10 * t21 - t9 * t13) - g(2) * (-t7 * t13 - t8 * t21) - (t13 * t27 - t21 * t24) * t43, 0, 0, 0, 0, 0, -g(1) * (t10 * t22 - t9 * t41) - g(2) * (t8 * t22 - t7 * t41) - (t15 * t32 + t22 * t24) * t43, -g(1) * (t10 * t25 + t9 * t42) - g(2) * (t8 * t25 + t7 * t42) - (-t15 * t33 + t24 * t25) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -g(1) * (-t10 * t26 - t17 * t38) - g(2) * (t19 * t38 - t8 * t26) - g(3) * (-t20 * t23 - t24 * t36), 0, t28 * pkin(3), 0, 0, 0, 0, 0, -t31 * t25, t31 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t22 + t9 * t25) - g(2) * (-t2 * t22 + t7 * t25) - g(3) * (-t18 * t32 - t6 * t22), -g(1) * (-t9 * t22 - t4 * t25) - g(2) * (-t2 * t25 - t7 * t22) - g(3) * (t18 * t33 - t6 * t25);];
taug_reg = t1;
