% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:35
% DurationCPUTime: 0.28s
% Computational Cost: add. (272->68), mult. (366->91), div. (0->0), fcn. (356->8), ass. (0->48)
t31 = cos(qJ(4));
t21 = t31 * pkin(4) + pkin(3);
t26 = qJ(2) + qJ(3);
t23 = sin(t26);
t24 = cos(t26);
t27 = -qJ(5) - pkin(8);
t62 = t24 * t21 - t23 * t27;
t61 = t24 * pkin(3) + t23 * pkin(8);
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t15 = g(1) * t33 + g(2) * t30;
t46 = t33 * t31;
t28 = sin(qJ(4));
t49 = t30 * t28;
t10 = t24 * t49 + t46;
t47 = t33 * t28;
t48 = t30 * t31;
t12 = -t24 * t47 + t48;
t52 = g(3) * t23;
t1 = -g(1) * t12 + g(2) * t10 + t28 * t52;
t7 = -g(3) * t24 + t15 * t23;
t29 = sin(qJ(2));
t60 = pkin(2) * t29;
t59 = pkin(3) * t23;
t58 = pkin(8) * t24;
t32 = cos(qJ(2));
t25 = t32 * pkin(2);
t22 = t25 + pkin(1);
t16 = t33 * t22;
t54 = g(2) * t16;
t34 = -pkin(7) - pkin(6);
t43 = pkin(4) * t28 - t34;
t41 = -t59 - t60;
t39 = g(1) * t30 - g(2) * t33;
t37 = t21 * t23 + t24 * t27;
t35 = -g(3) * t32 + t15 * t29;
t18 = t33 * t58;
t17 = t30 * t58;
t13 = t24 * t46 + t49;
t11 = -t24 * t48 + t47;
t9 = t39 * t23;
t8 = t15 * t24 + t52;
t6 = t7 * t31;
t5 = t7 * t28;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t31 * t52;
t14 = [0, 0, 0, 0, 0, 0, t39, t15, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t32, -t39 * t29, -t15, -g(1) * (-t30 * pkin(1) + t33 * pkin(6)) - g(2) * (t33 * pkin(1) + t30 * pkin(6)), 0, 0, 0, 0, 0, 0, t39 * t24, -t9, -t15, -g(1) * (-t30 * t22 - t33 * t34) - g(2) * (-t30 * t34 + t16), 0, 0, 0, 0, 0, 0, t4, t3, t9, -t54 + (g(1) * t34 - g(2) * t61) * t33 + (-g(1) * (-t22 - t61) + g(2) * t34) * t30, 0, 0, 0, 0, 0, 0, t4, t3, t9, -t54 + (-g(1) * t43 - g(2) * t62) * t33 + (-g(1) * (-t22 - t62) - g(2) * t43) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, g(3) * t29 + t15 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t35 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t41 * t33 + t18) - g(2) * (t41 * t30 + t17) - g(3) * (t25 + t61), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t25 + t62) + t15 * (t37 + t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-t33 * t59 + t18) - g(2) * (-t30 * t59 + t17) - g(3) * t61, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t62 + t15 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t14;
