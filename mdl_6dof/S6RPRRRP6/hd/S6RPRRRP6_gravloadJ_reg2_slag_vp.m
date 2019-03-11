% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t23 = g(1) * t42 + g(2) * t40;
t34 = pkin(10) + qJ(3);
t28 = cos(t34);
t41 = cos(qJ(4));
t55 = t42 * t41;
t39 = sin(qJ(4));
t60 = t40 * t39;
t14 = t28 * t60 + t55;
t56 = t42 * t39;
t59 = t40 * t41;
t16 = -t28 * t56 + t59;
t27 = sin(t34);
t66 = g(3) * t27;
t74 = -g(1) * t16 + g(2) * t14 + t39 * t66;
t35 = qJ(4) + qJ(5);
t29 = sin(t35);
t58 = t42 * t29;
t30 = cos(t35);
t61 = t40 * t30;
t11 = -t28 * t58 + t61;
t57 = t42 * t30;
t62 = t40 * t29;
t9 = t28 * t62 + t57;
t1 = -g(1) * t11 + g(2) * t9 + t29 * t66;
t7 = -g(3) * t28 + t23 * t27;
t43 = -pkin(9) - pkin(8);
t37 = cos(pkin(10));
t24 = t37 * pkin(2) + pkin(1);
t21 = t42 * t24;
t68 = g(2) * t21;
t64 = t39 * pkin(4);
t19 = pkin(5) * t29 + t64;
t63 = t19 * t28;
t38 = -pkin(7) - qJ(2);
t54 = t19 - t38;
t31 = t41 * pkin(4);
t20 = pkin(5) * t30 + t31;
t51 = -t38 + t64;
t50 = t28 * pkin(3) + t27 * pkin(8);
t22 = g(1) * t40 - g(2) * t42;
t18 = pkin(3) + t20;
t33 = -qJ(6) + t43;
t48 = t28 * t18 - t27 * t33;
t26 = t31 + pkin(3);
t46 = t28 * t26 - t27 * t43;
t17 = t28 * t55 + t60;
t15 = -t28 * t59 + t56;
t13 = t22 * t27;
t12 = t28 * t57 + t62;
t10 = -t28 * t61 + t58;
t8 = t23 * t28 + t66;
t6 = t7 * t30;
t5 = t7 * t29;
t4 = -g(1) * t10 - g(2) * t12;
t3 = -g(1) * t9 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t10 + t30 * t66;
t25 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t37, -t22 * sin(pkin(10)) -t23, -g(1) * (-t40 * pkin(1) + t42 * qJ(2)) - g(2) * (t42 * pkin(1) + t40 * qJ(2)) 0, 0, 0, 0, 0, 0, t22 * t28, -t13, -t23, -g(1) * (-t40 * t24 - t42 * t38) - g(2) * (-t40 * t38 + t21) 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -t68 + (g(1) * t38 - g(2) * t50) * t42 + (-g(1) * (-t24 - t50) + g(2) * t38) * t40, 0, 0, 0, 0, 0, 0, t4, t3, t13, -t68 + (-g(1) * t51 - g(2) * t46) * t42 + (-g(1) * (-t24 - t46) - g(2) * t51) * t40, 0, 0, 0, 0, 0, 0, t4, t3, t13, -t68 + (-g(1) * t54 - g(2) * t48) * t42 + (-g(1) * (-t24 - t48) - g(2) * t54) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t41, -t7 * t39, -t8, -g(3) * t50 + t23 * (pkin(3) * t27 - pkin(8) * t28) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t46 + t23 * (t26 * t27 + t28 * t43) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t48 + t23 * (t18 * t27 + t28 * t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, g(1) * t17 - g(2) * t15 + t41 * t66, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t74 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t40 * t20 - t42 * t63) - g(2) * (-t42 * t20 - t40 * t63) + t19 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t25;
