% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:06
% EndTime: 2019-12-05 17:26:07
% DurationCPUTime: 0.30s
% Computational Cost: add. (347->91), mult. (998->176), div. (0->0), fcn. (1304->14), ass. (0->62)
t70 = cos(qJ(3));
t35 = sin(pkin(6));
t40 = sin(qJ(4));
t69 = t35 * t40;
t44 = cos(qJ(4));
t68 = t35 * t44;
t36 = sin(pkin(5));
t38 = cos(pkin(6));
t67 = t36 * t38;
t42 = sin(qJ(2));
t66 = t36 * t42;
t45 = cos(qJ(2));
t65 = t36 * t45;
t41 = sin(qJ(3));
t64 = t38 * t41;
t39 = sin(qJ(5));
t63 = t39 * t44;
t62 = t41 * t42;
t61 = t41 * t45;
t43 = cos(qJ(5));
t60 = t43 * t44;
t59 = cos(pkin(5));
t58 = sin(pkin(11));
t37 = cos(pkin(11));
t57 = t35 * t36 * t37;
t56 = t35 * t66;
t55 = t38 * t70;
t54 = t70 * t42;
t53 = t70 * t45;
t52 = t36 * t58;
t51 = t37 * t59;
t50 = t59 * t35;
t49 = t35 * t52;
t48 = t59 * t58;
t29 = -t37 * t42 - t45 * t48;
t30 = t37 * t45 - t42 * t48;
t10 = t30 * t70 + (t29 * t38 + t49) * t41;
t19 = t41 * t50 + (t38 * t61 + t54) * t36;
t27 = -t58 * t42 + t45 * t51;
t20 = -t27 * t35 - t37 * t67;
t21 = -t29 * t35 + t38 * t52;
t26 = -t35 * t65 + t59 * t38;
t28 = t42 * t51 + t58 * t45;
t8 = t28 * t70 + (t27 * t38 - t57) * t41;
t47 = g(1) * (-t10 * t40 + t21 * t44) + g(2) * (t20 * t44 - t8 * t40) + g(3) * (-t19 * t40 + t26 * t44);
t18 = t36 * t62 - t70 * t50 - t53 * t67;
t7 = -t27 * t55 + t28 * t41 + t70 * t57;
t9 = -t29 * t55 + t30 * t41 - t70 * t49;
t46 = g(1) * t9 + g(2) * t7 + g(3) * t18;
t25 = (-t38 * t62 + t53) * t36;
t24 = (t38 * t54 + t61) * t36;
t17 = t25 * t44 + t40 * t56;
t16 = t29 * t70 - t30 * t64;
t15 = t29 * t41 + t30 * t55;
t14 = t27 * t70 - t28 * t64;
t13 = t27 * t41 + t28 * t55;
t12 = t19 * t44 + t26 * t40;
t6 = t16 * t44 + t30 * t69;
t5 = t14 * t44 + t28 * t69;
t4 = t10 * t44 + t21 * t40;
t2 = t20 * t40 + t8 * t44;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(1) * t29 - g(2) * t27 - g(3) * t65, g(1) * t30 + g(2) * t28 + g(3) * t66, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14 - g(3) * t25, g(1) * t15 + g(2) * t13 + g(3) * t24, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t5 - g(3) * t17, -g(1) * (-t16 * t40 + t30 * t68) - g(2) * (-t14 * t40 + t28 * t68) - g(3) * (-t25 * t40 + t44 * t56), 0, 0, 0, 0, 0, -g(1) * (t15 * t39 + t6 * t43) - g(2) * (t13 * t39 + t5 * t43) - g(3) * (t17 * t43 + t24 * t39), -g(1) * (t15 * t43 - t6 * t39) - g(2) * (t13 * t43 - t5 * t39) - g(3) * (-t17 * t39 + t24 * t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(1) * t10 + g(2) * t8 + g(3) * t19, 0, 0, 0, 0, 0, t46 * t44, -t46 * t40, 0, 0, 0, 0, 0, -g(1) * (t10 * t39 - t9 * t60) - g(2) * (t8 * t39 - t7 * t60) - g(3) * (-t18 * t60 + t19 * t39), -g(1) * (t10 * t43 + t9 * t63) - g(2) * (t8 * t43 + t7 * t63) - g(3) * (t18 * t63 + t19 * t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, g(1) * t4 + g(2) * t2 + g(3) * t12, 0, 0, 0, 0, 0, -t47 * t43, t47 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t39 + t9 * t43) - g(2) * (-t2 * t39 + t7 * t43) - g(3) * (-t12 * t39 + t18 * t43), -g(1) * (-t9 * t39 - t4 * t43) - g(2) * (-t2 * t43 - t7 * t39) - g(3) * (-t12 * t43 - t18 * t39);];
taug_reg = t1;
