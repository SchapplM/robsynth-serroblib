% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR2
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = cos(qJ(5));
t25 = t41 * pkin(5) + pkin(4);
t37 = qJ(3) + qJ(4);
t30 = sin(t37);
t32 = cos(t37);
t44 = -pkin(10) - pkin(9);
t77 = t32 * t25 - t30 * t44;
t58 = t32 * pkin(4) + t30 * pkin(9);
t35 = qJ(1) + pkin(11);
t27 = sin(t35);
t28 = cos(t35);
t18 = g(1) * t28 + g(2) * t27;
t38 = sin(qJ(5));
t61 = t32 * t38;
t14 = t27 * t61 + t28 * t41;
t16 = t27 * t41 - t28 * t61;
t69 = g(3) * t30;
t76 = -g(1) * t16 + g(2) * t14 + t38 * t69;
t7 = -g(3) * t32 + t18 * t30;
t39 = sin(qJ(3));
t75 = pkin(3) * t39;
t74 = pkin(4) * t30;
t40 = sin(qJ(1));
t67 = t40 * pkin(1);
t66 = t28 * t30;
t65 = t28 * t32;
t36 = qJ(5) + qJ(6);
t29 = sin(t36);
t64 = t29 * t32;
t31 = cos(t36);
t62 = t31 * t32;
t60 = t32 * t41;
t42 = cos(qJ(3));
t33 = t42 * pkin(3);
t26 = t33 + pkin(2);
t43 = cos(qJ(1));
t34 = t43 * pkin(1);
t59 = t28 * t26 + t34;
t45 = -pkin(8) - pkin(7);
t56 = pkin(5) * t38 - t45;
t54 = -t74 - t75;
t53 = g(1) * t27 - g(2) * t28;
t52 = g(1) * t40 - g(2) * t43;
t51 = -t28 * t45 - t67;
t49 = t25 * t30 + t32 * t44;
t46 = -g(3) * t42 + t18 * t39;
t21 = pkin(9) * t65;
t20 = t27 * t32 * pkin(9);
t17 = t27 * t38 + t28 * t60;
t15 = -t27 * t60 + t28 * t38;
t13 = t53 * t30;
t12 = t27 * t29 + t28 * t62;
t11 = t27 * t31 - t28 * t64;
t10 = -t27 * t62 + t28 * t29;
t9 = t27 * t64 + t28 * t31;
t8 = t18 * t32 + t69;
t6 = t7 * t41;
t5 = t7 * t38;
t4 = t7 * t31;
t3 = t7 * t29;
t2 = g(1) * t12 - g(2) * t10 + t31 * t69;
t1 = -g(1) * t11 + g(2) * t9 + t29 * t69;
t19 = [0, 0, 0, 0, 0, 0, t52, g(1) * t43 + g(2) * t40, 0, 0, 0, 0, 0, 0, 0, 0, t53, t18, 0, t52 * pkin(1), 0, 0, 0, 0, 0, 0, t53 * t42, -t53 * t39, -t18, -g(1) * (-t27 * pkin(2) + t28 * pkin(7) - t67) - g(2) * (t28 * pkin(2) + t27 * pkin(7) + t34) 0, 0, 0, 0, 0, 0, t53 * t32, -t13, -t18, -g(1) * (-t27 * t26 + t51) - g(2) * (-t27 * t45 + t59) 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -g(1) * t51 - g(2) * (pkin(4) * t65 + pkin(9) * t66 + t59) + (-g(1) * (-t26 - t58) + g(2) * t45) * t27, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, g(1) * t67 - g(2) * t59 + (-g(1) * t56 - g(2) * t77) * t28 + (-g(1) * (-t26 - t77) - g(2) * t56) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(3) * t39 + t18 * t42, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t46 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t54 * t28 + t21) - g(2) * (t54 * t27 + t20) - g(3) * (t33 + t58) 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * (t33 + t77) + t18 * (t49 + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-pkin(4) * t66 + t21) - g(2) * (-t27 * t74 + t20) - g(3) * t58, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t77 + t18 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, g(1) * t17 - g(2) * t15 + t41 * t69, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t76 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t19;
