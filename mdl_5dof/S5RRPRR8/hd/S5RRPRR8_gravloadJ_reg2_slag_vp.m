% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:13
% EndTime: 2019-12-31 20:18:14
% DurationCPUTime: 0.25s
% Computational Cost: add. (260->65), mult. (257->90), div. (0->0), fcn. (245->10), ass. (0->44)
t29 = qJ(2) + pkin(9);
t25 = qJ(4) + t29;
t20 = sin(t25);
t21 = cos(t25);
t50 = t21 * pkin(4) + t20 * pkin(8);
t33 = sin(qJ(1));
t36 = cos(qJ(1));
t16 = g(1) * t36 + g(2) * t33;
t3 = -g(3) * t21 + t16 * t20;
t49 = pkin(4) * t20;
t48 = pkin(8) * t21;
t47 = g(3) * t20;
t31 = sin(qJ(5));
t45 = t31 * t36;
t44 = t33 * t31;
t34 = cos(qJ(5));
t43 = t33 * t34;
t42 = t34 * t36;
t30 = -qJ(3) - pkin(6);
t24 = cos(t29);
t35 = cos(qJ(2));
t26 = t35 * pkin(2);
t40 = pkin(3) * t24 + t26;
t23 = sin(t29);
t32 = sin(qJ(2));
t12 = -pkin(2) * t32 - pkin(3) * t23;
t39 = t12 - t49;
t15 = g(1) * t33 - g(2) * t36;
t37 = -g(3) * t35 + t16 * t32;
t28 = -pkin(7) + t30;
t22 = t26 + pkin(1);
t14 = t36 * t48;
t13 = t33 * t48;
t11 = pkin(1) + t40;
t10 = t36 * t11;
t9 = t21 * t42 + t44;
t8 = -t21 * t45 + t43;
t7 = -t21 * t43 + t45;
t6 = t21 * t44 + t42;
t5 = t15 * t20;
t4 = t16 * t21 + t47;
t2 = t3 * t34;
t1 = t3 * t31;
t17 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t35, -t15 * t32, -t16, -g(1) * (-t33 * pkin(1) + pkin(6) * t36) - g(2) * (t36 * pkin(1) + t33 * pkin(6)), 0, 0, 0, 0, 0, 0, t15 * t24, -t15 * t23, -t16, -g(1) * (-t33 * t22 - t30 * t36) - g(2) * (t22 * t36 - t33 * t30), 0, 0, 0, 0, 0, 0, t15 * t21, -t5, -t16, -g(1) * (-t33 * t11 - t28 * t36) - g(2) * (-t33 * t28 + t10), 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5, -g(2) * t10 + (g(1) * t28 - g(2) * t50) * t36 + (-g(1) * (-t11 - t50) + g(2) * t28) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(3) * t32 + t16 * t35, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t24 + t16 * t23, g(3) * t23 + t16 * t24, 0, t37 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t40 - t16 * t12, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t39 * t36 + t14) - g(2) * (t39 * t33 + t13) - g(3) * (t40 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t36 * t49 + t14) - g(2) * (-t33 * t49 + t13) - g(3) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t6 + t31 * t47, g(1) * t9 - g(2) * t7 + t34 * t47, 0, 0;];
taug_reg = t17;
