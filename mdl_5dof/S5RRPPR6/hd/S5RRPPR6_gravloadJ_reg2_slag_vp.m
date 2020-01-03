% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = sin(qJ(1));
t29 = cos(qJ(1));
t9 = g(1) * t29 + g(2) * t27;
t21 = qJ(2) + pkin(8);
t15 = sin(t21);
t17 = cos(t21);
t1 = -g(3) * t17 + t9 * t15;
t26 = sin(qJ(2));
t50 = pkin(2) * t26;
t28 = cos(qJ(2));
t18 = t28 * pkin(2);
t13 = t18 + pkin(1);
t10 = t29 * t13;
t48 = g(2) * t10;
t46 = g(3) * t15;
t20 = pkin(9) + qJ(5);
t14 = sin(t20);
t44 = t27 * t14;
t16 = cos(t20);
t43 = t27 * t16;
t22 = sin(pkin(9));
t42 = t27 * t22;
t23 = cos(pkin(9));
t41 = t27 * t23;
t40 = t29 * t14;
t39 = t29 * t16;
t38 = t29 * t22;
t37 = t29 * t23;
t24 = -qJ(3) - pkin(6);
t36 = pkin(4) * t22 - t24;
t8 = g(1) * t27 - g(2) * t29;
t35 = t17 * pkin(3) + t15 * qJ(4);
t12 = t23 * pkin(4) + pkin(3);
t25 = -pkin(7) - qJ(4);
t34 = t17 * t12 - t15 * t25;
t30 = -g(3) * t28 + t9 * t26;
t7 = t8 * t15;
t6 = t17 * t39 + t44;
t5 = -t17 * t40 + t43;
t4 = -t17 * t43 + t40;
t3 = t17 * t44 + t39;
t2 = t9 * t17 + t46;
t11 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t28, -t8 * t26, -t9, -g(1) * (-t27 * pkin(1) + t29 * pkin(6)) - g(2) * (t29 * pkin(1) + t27 * pkin(6)), 0, 0, 0, 0, 0, 0, t8 * t17, -t7, -t9, -g(1) * (-t27 * t13 - t29 * t24) - g(2) * (-t27 * t24 + t10), 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t41 + t38) - g(2) * (t17 * t37 + t42), -g(1) * (t17 * t42 + t37) - g(2) * (-t17 * t38 + t41), t7, -t48 + (g(1) * t24 - g(2) * t35) * t29 + (-g(1) * (-t13 - t35) + g(2) * t24) * t27, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t7, -t48 + (-g(1) * t36 - g(2) * t34) * t29 + (-g(1) * (-t13 - t34) - g(2) * t36) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t26 + t9 * t28, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t30 * pkin(2), 0, 0, 0, 0, 0, 0, t1 * t23, -t1 * t22, -t2, -g(3) * (t18 + t35) + t9 * (pkin(3) * t15 - qJ(4) * t17 + t50), 0, 0, 0, 0, 0, 0, t1 * t16, -t1 * t14, -t2, -g(3) * (t18 + t34) + t9 * (t12 * t15 + t17 * t25 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t14 * t46, g(1) * t6 - g(2) * t4 + t16 * t46, 0, 0;];
taug_reg = t11;
