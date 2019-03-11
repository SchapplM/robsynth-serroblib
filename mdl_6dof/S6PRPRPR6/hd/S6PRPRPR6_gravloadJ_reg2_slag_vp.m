% Calculate inertial parameters regressor of gravitation load for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = sin(pkin(10));
t36 = sin(qJ(2));
t38 = cos(qJ(2));
t55 = cos(pkin(10));
t56 = cos(pkin(6));
t44 = t56 * t55;
t15 = t31 * t38 + t36 * t44;
t51 = t31 * t56;
t17 = -t36 * t51 + t38 * t55;
t74 = -g(1) * t17 - g(2) * t15;
t32 = sin(pkin(6));
t71 = g(3) * t32;
t14 = t31 * t36 - t38 * t44;
t70 = t14 * pkin(8);
t16 = t36 * t55 + t38 * t51;
t69 = t16 * pkin(8);
t30 = sin(pkin(11));
t68 = t14 * t30;
t67 = t16 * t30;
t29 = pkin(11) + qJ(6);
t27 = sin(t29);
t35 = sin(qJ(4));
t66 = t27 * t35;
t28 = cos(t29);
t65 = t28 * t35;
t64 = t30 * t35;
t63 = t30 * t38;
t62 = t31 * t32;
t61 = t32 * t36;
t60 = t32 * t38;
t33 = cos(pkin(11));
t59 = t33 * t35;
t58 = t35 * t36;
t57 = pkin(2) * t60 + qJ(3) * t61;
t54 = pkin(8) * t60 + t57;
t12 = t14 * pkin(2);
t53 = -t12 - t70;
t13 = t16 * pkin(2);
t52 = -t13 - t69;
t50 = t32 * t55;
t49 = t15 * qJ(3) - t12;
t48 = t17 * qJ(3) - t13;
t47 = g(3) * t54;
t37 = cos(qJ(4));
t46 = pkin(4) * t35 - qJ(5) * t37;
t26 = pkin(5) * t33 + pkin(4);
t34 = -pkin(9) - qJ(5);
t45 = t26 * t35 + t34 * t37;
t18 = t35 * t56 + t37 * t60;
t6 = -t16 * t37 + t35 * t62;
t8 = t14 * t37 + t35 * t50;
t41 = g(1) * t6 - g(2) * t8 + g(3) * t18;
t19 = -t35 * t60 + t37 * t56;
t7 = t16 * t35 + t37 * t62;
t9 = -t14 * t35 + t37 * t50;
t40 = g(1) * t7 - g(2) * t9 + g(3) * t19;
t4 = -g(1) * t16 - g(2) * t14 + g(3) * t60;
t39 = g(3) * t61 - t74;
t3 = t39 * t37;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t39, -g(1) * t48 - g(2) * t49 - g(3) * t57, 0, 0, 0, 0, 0, 0, -t39 * t35, -t3, -t4, -g(1) * (t48 - t69) - g(2) * (t49 - t70) - t47, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t59 - t67) - g(2) * (t15 * t59 - t68) - (t33 * t58 + t63) * t71, -g(1) * (-t16 * t33 - t17 * t64) - g(2) * (-t14 * t33 - t15 * t64) - (-t30 * t58 + t33 * t38) * t71, t3, -g(1) * t52 - g(2) * t53 - g(3) * (t46 * t61 + t54) + t74 * (qJ(3) + t46) 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t27 + t17 * t65) - g(2) * (-t14 * t27 + t15 * t65) - (t27 * t38 + t28 * t58) * t71, -g(1) * (-t16 * t28 - t17 * t66) - g(2) * (-t14 * t28 - t15 * t66) - (-t27 * t58 + t28 * t38) * t71, t3, -g(1) * (-pkin(5) * t67 + t52) - g(2) * (-pkin(5) * t68 + t53) - t47 - (pkin(5) * t63 + t36 * t45) * t71 + t74 * (qJ(3) + t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t33, -t41 * t30, -t40, -g(1) * (-pkin(4) * t6 + qJ(5) * t7) - g(2) * (pkin(4) * t8 - qJ(5) * t9) - g(3) * (-pkin(4) * t18 + qJ(5) * t19) 0, 0, 0, 0, 0, 0, t41 * t28, -t41 * t27, -t40, -g(1) * (-t26 * t6 - t34 * t7) - g(2) * (t26 * t8 + t34 * t9) - g(3) * (-t18 * t26 - t19 * t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t28 - t27 * t7) - g(2) * (t15 * t28 + t27 * t9) - g(3) * (-t19 * t27 + t28 * t61) -g(1) * (-t17 * t27 - t28 * t7) - g(2) * (-t15 * t27 + t28 * t9) - g(3) * (-t19 * t28 - t27 * t61) 0, 0;];
taug_reg  = t1;
