% Calculate minimal parameter regressor of gravitation load for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = qJ(2) + qJ(3);
t26 = pkin(10) + t33;
t21 = sin(t26);
t22 = cos(t26);
t42 = t22 * pkin(4) + t21 * qJ(5);
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t18 = g(1) * t39 + g(2) * t37;
t41 = -g(3) * t22 + t18 * t21;
t27 = sin(t33);
t59 = pkin(3) * t27;
t58 = pkin(4) * t21;
t57 = g(3) * t21;
t32 = pkin(11) + qJ(6);
t24 = sin(t32);
t55 = t37 * t24;
t25 = cos(t32);
t54 = t37 * t25;
t34 = sin(pkin(11));
t53 = t37 * t34;
t35 = cos(pkin(11));
t52 = t37 * t35;
t51 = t39 * t24;
t50 = t39 * t25;
t49 = t39 * t34;
t48 = t39 * t35;
t28 = cos(t33);
t23 = pkin(3) * t28;
t38 = cos(qJ(2));
t29 = t38 * pkin(2);
t47 = t23 + t29;
t46 = qJ(5) * t22;
t45 = t23 + t42;
t36 = sin(qJ(2));
t14 = -t36 * pkin(2) - t59;
t44 = t14 - t58;
t43 = -t58 - t59;
t17 = g(1) * t37 - g(2) * t39;
t6 = -g(3) * t28 + t18 * t27;
t31 = -qJ(4) - pkin(8) - pkin(7);
t16 = t39 * t46;
t15 = t37 * t46;
t13 = pkin(1) + t47;
t12 = t39 * t13;
t11 = t22 * t50 + t55;
t10 = -t22 * t51 + t54;
t9 = -t22 * t54 + t51;
t8 = t22 * t55 + t50;
t7 = g(3) * t27 + t18 * t28;
t5 = -t18 * t22 - t57;
t4 = t41 * t35;
t3 = t41 * t34;
t2 = t41 * t25;
t1 = t41 * t24;
t19 = [0, t17, t18, 0, 0, 0, 0, 0, t17 * t38, -t17 * t36, 0, 0, 0, 0, 0, t17 * t28, -t17 * t27, -t18, -g(1) * (-t37 * t13 - t39 * t31) - g(2) * (-t37 * t31 + t12) -g(1) * (-t22 * t52 + t49) - g(2) * (t22 * t48 + t53) -g(1) * (t22 * t53 + t48) - g(2) * (-t22 * t49 + t52) t17 * t21, -g(2) * t12 + (g(1) * t31 - g(2) * t42) * t39 + (-g(1) * (-t13 - t42) + g(2) * t31) * t37, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t38 + t18 * t36, g(3) * t36 + t18 * t38, 0, 0, 0, 0, 0, t6, t7, 0, -g(3) * t47 - t18 * t14, t4, -t3, t5, -g(1) * (t44 * t39 + t16) - g(2) * (t44 * t37 + t15) - g(3) * (t29 + t45) 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t6 * pkin(3), t4, -t3, t5, -g(1) * (t43 * t39 + t16) - g(2) * (t43 * t37 + t15) - g(3) * t45, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t24 * t57, g(1) * t11 - g(2) * t9 + t25 * t57;];
taug_reg  = t19;
