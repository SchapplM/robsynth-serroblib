% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = cos(pkin(11));
t26 = t41 * pkin(5) + pkin(4);
t39 = qJ(2) + qJ(3);
t35 = qJ(4) + t39;
t28 = sin(t35);
t29 = cos(t35);
t42 = -pkin(10) - qJ(5);
t84 = t29 * t26 - t28 * t42;
t83 = t29 * pkin(4) + t28 * qJ(5);
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t22 = g(1) * t46 + g(2) * t44;
t5 = -g(3) * t29 + t22 * t28;
t47 = -pkin(8) - pkin(7);
t33 = sin(t39);
t82 = pkin(3) * t33;
t81 = pkin(4) * t28;
t34 = cos(t39);
t27 = pkin(3) * t34;
t45 = cos(qJ(2));
t36 = t45 * pkin(2);
t65 = t27 + t36;
t18 = pkin(1) + t65;
t14 = t46 * t18;
t79 = g(2) * t14;
t77 = g(3) * t28;
t37 = pkin(11) + qJ(6);
t31 = sin(t37);
t74 = t44 * t31;
t32 = cos(t37);
t73 = t44 * t32;
t40 = sin(pkin(11));
t72 = t44 * t40;
t71 = t44 * t41;
t70 = t46 * t31;
t69 = t46 * t32;
t68 = t46 * t40;
t67 = t46 * t41;
t64 = qJ(5) * t29;
t63 = t27 + t83;
t38 = -pkin(9) + t47;
t62 = pkin(5) * t40 - t38;
t20 = t44 * t64;
t60 = -t44 * t81 + t20;
t21 = t46 * t64;
t59 = -t46 * t81 + t21;
t58 = t27 + t84;
t57 = -t81 - t82;
t56 = g(1) * t44 - g(2) * t46;
t53 = t26 * t28 + t29 * t42;
t52 = t53 * t44;
t51 = t53 * t46;
t7 = -g(3) * t34 + t22 * t33;
t43 = sin(qJ(2));
t48 = -g(3) * t45 + t22 * t43;
t30 = t36 + pkin(1);
t19 = -t43 * pkin(2) - t82;
t16 = t46 * t19;
t15 = t44 * t19;
t13 = t56 * t28;
t12 = t29 * t69 + t74;
t11 = -t29 * t70 + t73;
t10 = -t29 * t73 + t70;
t9 = t29 * t74 + t69;
t8 = g(3) * t33 + t22 * t34;
t6 = t22 * t29 + t77;
t4 = t5 * t41;
t3 = t5 * t40;
t2 = t5 * t32;
t1 = t5 * t31;
t17 = [0, 0, 0, 0, 0, 0, t56, t22, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t45, -t56 * t43, -t22, -g(1) * (-t44 * pkin(1) + t46 * pkin(7)) - g(2) * (t46 * pkin(1) + t44 * pkin(7)) 0, 0, 0, 0, 0, 0, t56 * t34, -t56 * t33, -t22, -g(1) * (-t44 * t30 - t46 * t47) - g(2) * (t46 * t30 - t44 * t47) 0, 0, 0, 0, 0, 0, t56 * t29, -t13, -t22, -g(1) * (-t44 * t18 - t46 * t38) - g(2) * (-t44 * t38 + t14) 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t71 + t68) - g(2) * (t29 * t67 + t72) -g(1) * (t29 * t72 + t67) - g(2) * (-t29 * t68 + t71) t13, -t79 + (g(1) * t38 - g(2) * t83) * t46 + (-g(1) * (-t18 - t83) + g(2) * t38) * t44, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t79 + (-g(1) * t62 - g(2) * t84) * t46 + (-g(1) * (-t18 - t84) - g(2) * t62) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(3) * t43 + t22 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t48 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t65 - t19 * t22, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t16 + t59) - g(2) * (t15 + t60) - g(3) * (t36 + t63) 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * (t16 - t51) - g(2) * (t15 - t52) - g(3) * (t36 + t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t46 * t57 + t21) - g(2) * (t44 * t57 + t20) - g(3) * t63, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(3) * t58 + t22 * (t53 + t82); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t59 - g(2) * t60 - g(3) * t83, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, g(1) * t51 + g(2) * t52 - g(3) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t31 * t77, g(1) * t12 - g(2) * t10 + t32 * t77, 0, 0;];
taug_reg  = t17;
