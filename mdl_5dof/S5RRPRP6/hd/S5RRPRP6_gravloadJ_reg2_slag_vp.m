% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t15 = g(1) * t31 + g(2) * t28;
t23 = qJ(2) + pkin(8);
t20 = cos(t23);
t29 = cos(qJ(4));
t39 = t31 * t29;
t26 = sin(qJ(4));
t42 = t28 * t26;
t10 = t20 * t42 + t39;
t40 = t31 * t26;
t41 = t28 * t29;
t12 = -t20 * t40 + t41;
t19 = sin(t23);
t44 = g(3) * t19;
t1 = -g(1) * t12 + g(2) * t10 + t26 * t44;
t7 = -g(3) * t20 + t15 * t19;
t27 = sin(qJ(2));
t50 = pkin(2) * t27;
t30 = cos(qJ(2));
t21 = t30 * pkin(2);
t18 = t21 + pkin(1);
t16 = t31 * t18;
t46 = g(2) * t16;
t25 = -qJ(3) - pkin(6);
t37 = pkin(4) * t26 - t25;
t36 = t20 * pkin(3) + t19 * pkin(7);
t14 = g(1) * t28 - g(2) * t31;
t17 = t29 * pkin(4) + pkin(3);
t24 = -qJ(5) - pkin(7);
t35 = t20 * t17 - t19 * t24;
t32 = -g(3) * t30 + t15 * t27;
t13 = t20 * t39 + t42;
t11 = -t20 * t41 + t40;
t9 = t14 * t19;
t8 = t15 * t20 + t44;
t6 = t7 * t29;
t5 = t7 * t26;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t29 * t44;
t22 = [0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t30, -t14 * t27, -t15, -g(1) * (-t28 * pkin(1) + t31 * pkin(6)) - g(2) * (t31 * pkin(1) + t28 * pkin(6)), 0, 0, 0, 0, 0, 0, t14 * t20, -t9, -t15, -g(1) * (-t28 * t18 - t31 * t25) - g(2) * (-t28 * t25 + t16), 0, 0, 0, 0, 0, 0, t4, t3, t9, -t46 + (g(1) * t25 - g(2) * t36) * t31 + (-g(1) * (-t18 - t36) + g(2) * t25) * t28, 0, 0, 0, 0, 0, 0, t4, t3, t9, -t46 + (-g(1) * t37 - g(2) * t35) * t31 + (-g(1) * (-t18 - t35) - g(2) * t37) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t27 + t15 * t30, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t32 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t21 + t36) + t15 * (pkin(3) * t19 - pkin(7) * t20 + t50), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t21 + t35) + t15 * (t17 * t19 + t20 * t24 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t22;
