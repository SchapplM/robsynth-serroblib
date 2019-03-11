% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = -pkin(10) - pkin(9);
t38 = sin(pkin(6));
t76 = g(3) * t38;
t39 = sin(qJ(4));
t75 = t39 * pkin(4);
t37 = qJ(4) + qJ(5);
t33 = sin(t37);
t22 = pkin(5) * t33 + t75;
t74 = pkin(8) + t22;
t43 = cos(qJ(3));
t73 = t33 * t43;
t34 = cos(t37);
t72 = t34 * t43;
t41 = sin(qJ(2));
t71 = t38 * t41;
t44 = cos(qJ(2));
t70 = t38 * t44;
t69 = t39 * t41;
t68 = t39 * t43;
t42 = cos(qJ(4));
t67 = t42 * t43;
t66 = t43 * t44;
t65 = pkin(2) * t70 + pkin(8) * t71;
t35 = t42 * pkin(4);
t23 = pkin(5) * t34 + t35;
t64 = cos(pkin(6));
t63 = cos(pkin(11));
t62 = sin(pkin(11));
t61 = pkin(8) + t75;
t60 = g(3) * t65;
t52 = t64 * t63;
t15 = t41 * t62 - t44 * t52;
t13 = t15 * pkin(2);
t16 = t41 * t52 + t44 * t62;
t59 = t16 * pkin(8) - t13;
t51 = t64 * t62;
t17 = t41 * t63 + t44 * t51;
t14 = t17 * pkin(2);
t18 = -t41 * t51 + t44 * t63;
t58 = t18 * pkin(8) - t14;
t57 = t38 * t63;
t56 = t38 * t62;
t40 = sin(qJ(3));
t55 = pkin(3) * t43 + pkin(9) * t40;
t21 = pkin(3) + t23;
t36 = -qJ(6) + t45;
t54 = t21 * t43 - t36 * t40;
t32 = t35 + pkin(3);
t53 = t32 * t43 - t40 * t45;
t11 = t18 * t40 - t43 * t56;
t19 = t40 * t71 - t43 * t64;
t9 = t16 * t40 + t43 * t57;
t50 = g(1) * t11 + g(2) * t9 + g(3) * t19;
t10 = t16 * t43 - t40 * t57;
t12 = t18 * t43 + t40 * t56;
t20 = t40 * t64 + t43 * t71;
t49 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t48 = -g(1) * t17 - g(2) * t15 + g(3) * t70;
t47 = g(1) * t18 + g(2) * t16 + g(3) * t71;
t1 = -g(1) * (-t12 * t33 + t17 * t34) - g(2) * (-t10 * t33 + t15 * t34) - g(3) * (-t20 * t33 - t34 * t70);
t46 = -g(1) * (-t12 * t39 + t17 * t42) - g(2) * (-t10 * t39 + t15 * t42) - g(3) * (-t20 * t39 - t42 * t70);
t8 = t48 * t40;
t6 = t50 * t34;
t5 = t50 * t33;
t4 = -g(1) * (-t17 * t72 + t18 * t33) - g(2) * (-t15 * t72 + t16 * t33) - (t33 * t41 + t34 * t66) * t76;
t3 = -g(1) * (t17 * t73 + t18 * t34) - g(2) * (t15 * t73 + t16 * t34) - (-t33 * t66 + t34 * t41) * t76;
t2 = -g(1) * (-t12 * t34 - t17 * t33) - g(2) * (-t10 * t34 - t15 * t33) - g(3) * (-t20 * t34 + t33 * t70);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t43, t8, -t47, -g(1) * t58 - g(2) * t59 - t60, 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t67 + t18 * t39) - g(2) * (-t15 * t67 + t16 * t39) - (t42 * t66 + t69) * t76, -g(1) * (t17 * t68 + t18 * t42) - g(2) * (t15 * t68 + t16 * t42) - (-t39 * t66 + t41 * t42) * t76, -t8, -g(1) * (-t17 * t55 + t58) - g(2) * (-t15 * t55 + t59) - g(3) * (t55 * t70 + t65) 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (-t17 * t53 + t18 * t61 - t14) - g(2) * (-t15 * t53 + t16 * t61 - t13) - t60 - (pkin(4) * t69 + t44 * t53) * t76, 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (-t17 * t54 + t18 * t74 - t14) - g(2) * (-t15 * t54 + t16 * t74 - t13) - t60 - (t22 * t41 + t44 * t54) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t49, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t42, -t50 * t39, -t49, -g(1) * (-pkin(3) * t11 + pkin(9) * t12) - g(2) * (-pkin(3) * t9 + pkin(9) * t10) - g(3) * (-pkin(3) * t19 + pkin(9) * t20) 0, 0, 0, 0, 0, 0, t6, -t5, -t49, -g(1) * (-t11 * t32 - t12 * t45) - g(2) * (-t10 * t45 - t32 * t9) - g(3) * (-t19 * t32 - t20 * t45) 0, 0, 0, 0, 0, 0, t6, -t5, -t49, -g(1) * (-t11 * t21 - t12 * t36) - g(2) * (-t10 * t36 - t21 * t9) - g(3) * (-t19 * t21 - t20 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -g(1) * (-t12 * t42 - t17 * t39) - g(2) * (-t10 * t42 - t15 * t39) - g(3) * (-t20 * t42 + t39 * t70) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t46 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t12 * t22 + t17 * t23) - g(2) * (-t10 * t22 + t15 * t23) - g(3) * (-t20 * t22 - t23 * t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50;];
taug_reg  = t7;
