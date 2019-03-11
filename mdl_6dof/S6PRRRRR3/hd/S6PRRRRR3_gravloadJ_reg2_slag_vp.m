% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = -pkin(10) - pkin(9);
t38 = sin(pkin(6));
t77 = g(3) * t38;
t39 = sin(qJ(4));
t76 = t39 * pkin(4);
t36 = qJ(4) + qJ(5);
t31 = sin(t36);
t20 = pkin(5) * t31 + t76;
t75 = pkin(8) + t20;
t33 = qJ(6) + t36;
t28 = sin(t33);
t43 = cos(qJ(3));
t74 = t28 * t43;
t29 = cos(t33);
t73 = t29 * t43;
t72 = t31 * t43;
t32 = cos(t36);
t71 = t32 * t43;
t41 = sin(qJ(2));
t70 = t38 * t41;
t69 = t38 * t43;
t44 = cos(qJ(2));
t68 = t38 * t44;
t67 = t39 * t41;
t66 = t39 * t43;
t42 = cos(qJ(4));
t65 = t42 * t43;
t64 = t43 * t44;
t63 = pkin(2) * t68 + pkin(8) * t70;
t34 = t42 * pkin(4);
t21 = pkin(5) * t32 + t34;
t62 = cos(pkin(6));
t61 = cos(pkin(12));
t60 = pkin(8) + t76;
t59 = g(3) * t63;
t37 = sin(pkin(12));
t51 = t62 * t61;
t13 = t37 * t41 - t44 * t51;
t11 = t13 * pkin(2);
t14 = t37 * t44 + t41 * t51;
t58 = t14 * pkin(8) - t11;
t56 = t37 * t62;
t15 = t61 * t41 + t44 * t56;
t12 = t15 * pkin(2);
t16 = -t41 * t56 + t61 * t44;
t57 = t16 * pkin(8) - t12;
t55 = t38 * t61;
t40 = sin(qJ(3));
t54 = pkin(3) * t43 + pkin(9) * t40;
t19 = pkin(3) + t21;
t35 = -pkin(11) + t45;
t53 = t19 * t43 - t35 * t40;
t30 = t34 + pkin(3);
t52 = t30 * t43 - t40 * t45;
t17 = -t40 * t70 + t62 * t43;
t7 = -t14 * t40 - t43 * t55;
t9 = -t16 * t40 + t37 * t69;
t50 = g(1) * t9 + g(2) * t7 + g(3) * t17;
t10 = t37 * t38 * t40 + t16 * t43;
t18 = t62 * t40 + t41 * t69;
t8 = t14 * t43 - t40 * t55;
t49 = g(1) * t10 + g(2) * t8 + g(3) * t18;
t48 = -g(1) * t15 - g(2) * t13 + g(3) * t68;
t47 = g(1) * t16 + g(2) * t14 + g(3) * t70;
t3 = -g(1) * (-t10 * t31 + t15 * t32) - g(2) * (t13 * t32 - t8 * t31) - g(3) * (-t18 * t31 - t32 * t68);
t46 = -g(1) * (-t10 * t39 + t15 * t42) - g(2) * (t13 * t42 - t8 * t39) - g(3) * (-t18 * t39 - t42 * t68);
t6 = t48 * t40;
t4 = -g(1) * (-t10 * t32 - t15 * t31) - g(2) * (-t13 * t31 - t8 * t32) - g(3) * (-t18 * t32 + t31 * t68);
t2 = -g(1) * (-t10 * t29 - t15 * t28) - g(2) * (-t13 * t28 - t8 * t29) - g(3) * (-t18 * t29 + t28 * t68);
t1 = -g(1) * (-t10 * t28 + t15 * t29) - g(2) * (t13 * t29 - t8 * t28) - g(3) * (-t18 * t28 - t29 * t68);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t43, t6, -t47, -g(1) * t57 - g(2) * t58 - t59, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t65 + t16 * t39) - g(2) * (-t13 * t65 + t14 * t39) - (t42 * t64 + t67) * t77, -g(1) * (t15 * t66 + t16 * t42) - g(2) * (t13 * t66 + t14 * t42) - (-t39 * t64 + t41 * t42) * t77, -t6, -g(1) * (-t54 * t15 + t57) - g(2) * (-t54 * t13 + t58) - g(3) * (t54 * t68 + t63) 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t71 + t16 * t31) - g(2) * (-t13 * t71 + t14 * t31) - (t31 * t41 + t32 * t64) * t77, -g(1) * (t15 * t72 + t16 * t32) - g(2) * (t13 * t72 + t14 * t32) - (-t31 * t64 + t32 * t41) * t77, -t6, -g(1) * (-t52 * t15 + t60 * t16 - t12) - g(2) * (-t52 * t13 + t60 * t14 - t11) - t59 - (pkin(4) * t67 + t52 * t44) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t73 + t16 * t28) - g(2) * (-t13 * t73 + t14 * t28) - (t28 * t41 + t29 * t64) * t77, -g(1) * (t15 * t74 + t16 * t29) - g(2) * (t13 * t74 + t14 * t29) - (-t28 * t64 + t29 * t41) * t77, -t6, -g(1) * (-t53 * t15 + t75 * t16 - t12) - g(2) * (-t53 * t13 + t75 * t14 - t11) - t59 - (t20 * t41 + t53 * t44) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t42, t50 * t39, -t49, -g(1) * (t9 * pkin(3) + t10 * pkin(9)) - g(2) * (t7 * pkin(3) + t8 * pkin(9)) - g(3) * (t17 * pkin(3) + t18 * pkin(9)) 0, 0, 0, 0, 0, 0, -t50 * t32, t50 * t31, -t49, -g(1) * (-t10 * t45 + t9 * t30) - g(2) * (t7 * t30 - t8 * t45) - g(3) * (t17 * t30 - t18 * t45) 0, 0, 0, 0, 0, 0, -t50 * t29, t50 * t28, -t49, -g(1) * (-t10 * t35 + t9 * t19) - g(2) * (t7 * t19 - t8 * t35) - g(3) * (t17 * t19 - t18 * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -g(1) * (-t10 * t42 - t15 * t39) - g(2) * (-t13 * t39 - t8 * t42) - g(3) * (-t18 * t42 + t39 * t68) 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t46 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t10 * t20 + t15 * t21) - g(2) * (t13 * t21 - t8 * t20) - g(3) * (-t18 * t20 - t21 * t68); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t5;
