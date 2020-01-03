% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR13
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t16 = g(1) * t36 + g(2) * t34;
t29 = pkin(9) + qJ(4);
t20 = sin(t29);
t33 = sin(qJ(2));
t52 = g(3) * t33;
t21 = cos(t29);
t35 = cos(qJ(2));
t48 = t34 * t35;
t7 = t20 * t48 + t21 * t36;
t47 = t35 * t36;
t9 = -t20 * t47 + t34 * t21;
t59 = -g(1) * t9 + g(2) * t7 + t20 * t52;
t11 = -g(3) * t35 + t16 * t33;
t30 = sin(pkin(9));
t56 = pkin(3) * t30;
t55 = g(1) * t34;
t31 = cos(pkin(9));
t19 = t31 * pkin(3) + pkin(2);
t50 = t30 * t36;
t49 = t33 * t36;
t32 = -pkin(7) - qJ(3);
t46 = t36 * pkin(1) + t34 * pkin(6);
t44 = -g(2) * t36 + t55;
t43 = pkin(2) * t35 + qJ(3) * t33;
t14 = pkin(4) * t21 + t19;
t28 = -pkin(8) + t32;
t41 = t14 * t35 - t28 * t33;
t39 = t19 * t35 - t32 * t33;
t25 = t36 * pkin(6);
t22 = qJ(5) + t29;
t18 = cos(t22);
t17 = sin(t22);
t15 = pkin(4) * t20 + t56;
t13 = t44 * t33;
t12 = t16 * t35 + t52;
t10 = t34 * t20 + t21 * t47;
t8 = t20 * t36 - t21 * t48;
t6 = t34 * t17 + t18 * t47;
t5 = -t17 * t47 + t34 * t18;
t4 = t17 * t36 - t18 * t48;
t3 = t17 * t48 + t18 * t36;
t2 = g(1) * t6 - g(2) * t4 + t18 * t52;
t1 = -g(1) * t5 + g(2) * t3 + t17 * t52;
t23 = [0, 0, 0, 0, 0, 0, t44, t16, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t35, -t13, -t16, -g(1) * (-t34 * pkin(1) + t25) - g(2) * t46, 0, 0, 0, 0, 0, 0, -g(1) * (-t31 * t48 + t50) - g(2) * (t34 * t30 + t31 * t47), -g(1) * (t30 * t48 + t31 * t36) - g(2) * (-t30 * t47 + t34 * t31), t13, -g(1) * t25 - g(2) * (t43 * t36 + t46) - (-pkin(1) - t43) * t55, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t13, -g(1) * (pkin(3) * t50 + t25) - g(2) * (t19 * t47 - t32 * t49 + t46) + (-g(1) * (-pkin(1) - t39) - g(2) * t56) * t34, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t13, -g(1) * (t15 * t36 + t25) - g(2) * (t14 * t47 - t28 * t49 + t46) + (-g(1) * (-pkin(1) - t41) - g(2) * t15) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t31, -t11 * t30, -t12, -g(3) * t43 + t16 * (pkin(2) * t33 - qJ(3) * t35), 0, 0, 0, 0, 0, 0, t11 * t21, -t11 * t20, -t12, -g(3) * t39 + t16 * (t19 * t33 + t32 * t35), 0, 0, 0, 0, 0, 0, t11 * t18, -t11 * t17, -t12, -g(3) * t41 + t16 * (t14 * t33 + t28 * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, g(1) * t10 - g(2) * t8 + t21 * t52, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t59 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t23;
