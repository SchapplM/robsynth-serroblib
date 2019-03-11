% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t24 = g(1) * t46 + g(2) * t44;
t37 = pkin(11) + qJ(3);
t31 = cos(t37);
t45 = cos(qJ(4));
t60 = t46 * t45;
t43 = sin(qJ(4));
t67 = t44 * t43;
t15 = t31 * t67 + t60;
t61 = t46 * t43;
t66 = t44 * t45;
t17 = -t31 * t61 + t66;
t30 = sin(t37);
t75 = g(3) * t30;
t83 = -g(1) * t17 + g(2) * t15 + t43 * t75;
t39 = qJ(4) + qJ(5);
t33 = cos(t39);
t62 = t46 * t33;
t32 = sin(t39);
t69 = t44 * t32;
t10 = t31 * t69 + t62;
t63 = t46 * t32;
t68 = t44 * t33;
t12 = -t31 * t63 + t68;
t3 = -g(1) * t12 + g(2) * t10 + t32 * t75;
t49 = -g(3) * t31 + t24 * t30;
t47 = -pkin(9) - pkin(8);
t41 = cos(pkin(11));
t25 = t41 * pkin(2) + pkin(1);
t22 = t46 * t25;
t77 = g(2) * t22;
t73 = t43 * pkin(4);
t20 = pkin(5) * t32 + t73;
t72 = t20 * t31;
t34 = qJ(6) + t39;
t27 = sin(t34);
t71 = t44 * t27;
t28 = cos(t34);
t70 = t44 * t28;
t65 = t46 * t27;
t64 = t46 * t28;
t42 = -pkin(7) - qJ(2);
t59 = t20 - t42;
t35 = t45 * pkin(4);
t21 = pkin(5) * t33 + t35;
t56 = -t42 + t73;
t55 = t31 * pkin(3) + t30 * pkin(8);
t23 = g(1) * t44 - g(2) * t46;
t19 = pkin(3) + t21;
t38 = -pkin(10) + t47;
t53 = t31 * t19 - t30 * t38;
t29 = t35 + pkin(3);
t51 = t31 * t29 - t30 * t47;
t18 = t31 * t60 + t67;
t16 = -t31 * t66 + t61;
t14 = t23 * t30;
t13 = t31 * t62 + t69;
t11 = -t31 * t68 + t63;
t9 = t31 * t64 + t71;
t8 = -t31 * t65 + t70;
t7 = -t31 * t70 + t65;
t6 = t31 * t71 + t64;
t5 = t24 * t31 + t75;
t4 = g(1) * t13 - g(2) * t11 + t33 * t75;
t2 = g(1) * t9 - g(2) * t7 + t28 * t75;
t1 = -g(1) * t8 + g(2) * t6 + t27 * t75;
t26 = [0, 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t41, -t23 * sin(pkin(11)) -t24, -g(1) * (-t44 * pkin(1) + t46 * qJ(2)) - g(2) * (t46 * pkin(1) + t44 * qJ(2)) 0, 0, 0, 0, 0, 0, t23 * t31, -t14, -t24, -g(1) * (-t44 * t25 - t46 * t42) - g(2) * (-t44 * t42 + t22) 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t14, -t77 + (g(1) * t42 - g(2) * t55) * t46 + (-g(1) * (-t25 - t55) + g(2) * t42) * t44, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, t14, -t77 + (-g(1) * t56 - g(2) * t51) * t46 + (-g(1) * (-t25 - t51) - g(2) * t56) * t44, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t14, -t77 + (-g(1) * t59 - g(2) * t53) * t46 + (-g(1) * (-t25 - t53) - g(2) * t59) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t5, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t45, -t49 * t43, -t5, -g(3) * t55 + t24 * (pkin(3) * t30 - pkin(8) * t31) 0, 0, 0, 0, 0, 0, t49 * t33, -t49 * t32, -t5, -g(3) * t51 + t24 * (t29 * t30 + t31 * t47) 0, 0, 0, 0, 0, 0, t49 * t28, -t49 * t27, -t5, -g(3) * t53 + t24 * (t19 * t30 + t31 * t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, g(1) * t18 - g(2) * t16 + t45 * t75, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t83 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t44 * t21 - t46 * t72) - g(2) * (-t46 * t21 - t44 * t72) + t20 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t26;
