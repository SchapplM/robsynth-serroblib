% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR9
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
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:21:58
% EndTime: 2019-12-31 20:21:59
% DurationCPUTime: 0.30s
% Computational Cost: add. (250->70), mult. (306->99), div. (0->0), fcn. (306->10), ass. (0->50)
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t14 = g(1) * t32 + g(2) * t29;
t24 = qJ(2) + pkin(9);
t19 = cos(t24);
t27 = sin(qJ(4));
t44 = t32 * t27;
t30 = cos(qJ(4));
t47 = t29 * t30;
t11 = -t19 * t44 + t47;
t18 = sin(t24);
t52 = g(3) * t18;
t43 = t32 * t30;
t48 = t29 * t27;
t9 = t19 * t48 + t43;
t59 = -g(1) * t11 + g(2) * t9 + t27 * t52;
t36 = -g(3) * t19 + t14 * t18;
t28 = sin(qJ(2));
t57 = pkin(2) * t28;
t31 = cos(qJ(2));
t22 = t31 * pkin(2);
t17 = t22 + pkin(1);
t15 = t32 * t17;
t54 = g(2) * t15;
t25 = qJ(4) + qJ(5);
t20 = sin(t25);
t50 = t29 * t20;
t21 = cos(t25);
t49 = t29 * t21;
t46 = t32 * t20;
t45 = t32 * t21;
t26 = -qJ(3) - pkin(6);
t41 = pkin(4) * t27 - t26;
t40 = t19 * pkin(3) + t18 * pkin(7);
t13 = g(1) * t29 - g(2) * t32;
t16 = t30 * pkin(4) + pkin(3);
t33 = -pkin(8) - pkin(7);
t39 = t19 * t16 - t18 * t33;
t34 = -g(3) * t31 + t14 * t28;
t12 = t19 * t43 + t48;
t10 = -t19 * t47 + t44;
t8 = t13 * t18;
t7 = t19 * t45 + t50;
t6 = -t19 * t46 + t49;
t5 = -t19 * t49 + t46;
t4 = t19 * t50 + t45;
t3 = t14 * t19 + t52;
t2 = g(1) * t7 - g(2) * t5 + t21 * t52;
t1 = -g(1) * t6 + g(2) * t4 + t20 * t52;
t23 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t31, -t13 * t28, -t14, -g(1) * (-t29 * pkin(1) + t32 * pkin(6)) - g(2) * (t32 * pkin(1) + t29 * pkin(6)), 0, 0, 0, 0, 0, 0, t13 * t19, -t8, -t14, -g(1) * (-t29 * t17 - t32 * t26) - g(2) * (-t29 * t26 + t15), 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t8, -t54 + (g(1) * t26 - g(2) * t40) * t32 + (-g(1) * (-t17 - t40) + g(2) * t26) * t29, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, -g(1) * t4 - g(2) * t6, t8, -t54 + (-g(1) * t41 - g(2) * t39) * t32 + (-g(1) * (-t17 - t39) - g(2) * t41) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, g(3) * t28 + t14 * t31, 0, 0, 0, 0, 0, 0, 0, 0, t36, t3, 0, t34 * pkin(2), 0, 0, 0, 0, 0, 0, t36 * t30, -t36 * t27, -t3, -g(3) * (t22 + t40) + t14 * (pkin(3) * t18 - pkin(7) * t19 + t57), 0, 0, 0, 0, 0, 0, t36 * t21, -t36 * t20, -t3, -g(3) * (t22 + t39) + t14 * (t16 * t18 + t19 * t33 + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, g(1) * t12 - g(2) * t10 + t30 * t52, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t59 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t23;
