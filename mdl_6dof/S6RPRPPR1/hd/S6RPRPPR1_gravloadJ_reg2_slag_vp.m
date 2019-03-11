% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t24 = qJ(3) + pkin(10);
t16 = sin(t24);
t19 = cos(t24);
t59 = t19 * pkin(4) + t16 * qJ(5);
t25 = qJ(1) + pkin(9);
t17 = sin(t25);
t20 = cos(t25);
t9 = g(1) * t20 + g(2) * t17;
t1 = -g(3) * t19 + t9 * t16;
t30 = sin(qJ(3));
t58 = pkin(3) * t30;
t55 = g(3) * t16;
t31 = sin(qJ(1));
t52 = t31 * pkin(1);
t51 = t17 * t19;
t26 = sin(pkin(11));
t50 = t17 * t26;
t27 = cos(pkin(11));
t49 = t17 * t27;
t23 = pkin(11) + qJ(6);
t15 = sin(t23);
t48 = t20 * t15;
t18 = cos(t23);
t47 = t20 * t18;
t46 = t20 * t26;
t45 = t20 * t27;
t32 = cos(qJ(3));
t21 = t32 * pkin(3);
t14 = t21 + pkin(2);
t33 = cos(qJ(1));
t22 = t33 * pkin(1);
t44 = t20 * t14 + t22;
t28 = -qJ(4) - pkin(7);
t42 = pkin(5) * t26 - t28;
t8 = g(1) * t17 - g(2) * t20;
t41 = g(1) * t31 - g(2) * t33;
t40 = -t20 * t28 - t52;
t13 = t27 * pkin(5) + pkin(4);
t29 = -pkin(8) - qJ(5);
t38 = t19 * t13 - t16 * t29;
t34 = -g(3) * t32 + t9 * t30;
t7 = t8 * t16;
t6 = t17 * t15 + t19 * t47;
t5 = t17 * t18 - t19 * t48;
t4 = -t18 * t51 + t48;
t3 = t15 * t51 + t47;
t2 = t9 * t19 + t55;
t10 = [0, 0, 0, 0, 0, 0, t41, g(1) * t33 + g(2) * t31, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t41 * pkin(1), 0, 0, 0, 0, 0, 0, t8 * t32, -t8 * t30, -t9, -g(1) * (-t17 * pkin(2) + t20 * pkin(7) - t52) - g(2) * (t20 * pkin(2) + t17 * pkin(7) + t22) 0, 0, 0, 0, 0, 0, t8 * t19, -t7, -t9, -g(1) * (-t17 * t14 + t40) - g(2) * (-t17 * t28 + t44) 0, 0, 0, 0, 0, 0, -g(1) * (-t19 * t49 + t46) - g(2) * (t19 * t45 + t50) -g(1) * (t19 * t50 + t45) - g(2) * (-t19 * t46 + t49) t7, -g(1) * t40 - g(2) * (t59 * t20 + t44) + (-g(1) * (-t14 - t59) + g(2) * t28) * t17, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t7, g(1) * t52 - g(2) * t44 + (-g(1) * t42 - g(2) * t38) * t20 + (-g(1) * (-t14 - t38) - g(2) * t42) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, g(3) * t30 + t9 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t34 * pkin(3), 0, 0, 0, 0, 0, 0, t1 * t27, -t1 * t26, -t2, -g(3) * (t21 + t59) + t9 * (pkin(4) * t16 - qJ(5) * t19 + t58) 0, 0, 0, 0, 0, 0, t1 * t18, -t1 * t15, -t2, -g(3) * (t21 + t38) + t9 * (t13 * t16 + t19 * t29 + t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t15 * t55, g(1) * t6 - g(2) * t4 + t18 * t55, 0, 0;];
taug_reg  = t10;
