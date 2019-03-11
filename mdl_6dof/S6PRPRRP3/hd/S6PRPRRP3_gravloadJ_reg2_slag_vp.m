% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = cos(pkin(11));
t31 = t38 * pkin(3) + pkin(2);
t37 = sin(pkin(6));
t44 = cos(qJ(2));
t62 = t37 * t44;
t25 = t31 * t62;
t69 = g(3) * t25;
t68 = g(3) * t37;
t42 = sin(qJ(2));
t55 = cos(pkin(10));
t56 = cos(pkin(6));
t49 = t56 * t55;
t54 = sin(pkin(10));
t22 = t42 * t49 + t54 * t44;
t41 = sin(qJ(5));
t67 = t22 * t41;
t48 = t56 * t54;
t24 = -t42 * t48 + t55 * t44;
t66 = t24 * t41;
t35 = pkin(11) + qJ(4);
t34 = cos(t35);
t65 = t34 * t41;
t43 = cos(qJ(5));
t64 = t34 * t43;
t63 = t37 * t42;
t40 = -pkin(8) - qJ(3);
t61 = t40 * t42;
t60 = t41 * t44;
t59 = t43 * t44;
t21 = t54 * t42 - t44 * t49;
t58 = -t21 * t31 - t22 * t40;
t23 = t55 * t42 + t44 * t48;
t57 = -t23 * t31 - t24 * t40;
t53 = t37 * t55;
t52 = t37 * t54;
t33 = sin(t35);
t51 = pkin(4) * t34 + pkin(9) * t33;
t32 = t43 * pkin(5) + pkin(4);
t39 = -qJ(6) - pkin(9);
t50 = t32 * t34 - t33 * t39;
t11 = t22 * t33 + t34 * t53;
t13 = t24 * t33 - t34 * t52;
t17 = t33 * t63 - t56 * t34;
t47 = g(1) * t13 + g(2) * t11 + g(3) * t17;
t12 = t22 * t34 - t33 * t53;
t14 = t24 * t34 + t33 * t52;
t18 = t56 * t33 + t34 * t63;
t46 = g(1) * t14 + g(2) * t12 + g(3) * t18;
t9 = -g(1) * t23 - g(2) * t21 + g(3) * t62;
t45 = g(1) * t24 + g(2) * t22 + g(3) * t63;
t1 = -g(1) * (-t14 * t41 + t23 * t43) - g(2) * (-t12 * t41 + t21 * t43) - g(3) * (-t18 * t41 - t37 * t59);
t8 = t9 * t33;
t6 = t47 * t43;
t5 = t47 * t41;
t4 = -g(1) * (-t23 * t64 + t66) - g(2) * (-t21 * t64 + t67) - (t34 * t59 + t41 * t42) * t68;
t3 = -g(1) * (t23 * t65 + t24 * t43) - g(2) * (t21 * t65 + t22 * t43) - (-t34 * t60 + t42 * t43) * t68;
t2 = -g(1) * (-t14 * t43 - t23 * t41) - g(2) * (-t12 * t43 - t21 * t41) - g(3) * (-t18 * t43 + t37 * t60);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t45, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * t38, t9 * sin(pkin(11)) -t45, -g(1) * (-t23 * pkin(2) + t24 * qJ(3)) - g(2) * (-t21 * pkin(2) + t22 * qJ(3)) - (pkin(2) * t44 + qJ(3) * t42) * t68, 0, 0, 0, 0, 0, 0, -t9 * t34, t8, -t45, -g(1) * t57 - g(2) * t58 - g(3) * (-t37 * t61 + t25) 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (-t51 * t23 + t57) - g(2) * (-t51 * t21 + t58) - t69 - (t51 * t44 - t61) * t68, 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (pkin(5) * t66 - t50 * t23 + t57) - g(2) * (pkin(5) * t67 - t50 * t21 + t58) - t69 - (t50 * t44 + (pkin(5) * t41 - t40) * t42) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t46, -g(1) * (-t13 * pkin(4) + t14 * pkin(9)) - g(2) * (-t11 * pkin(4) + t12 * pkin(9)) - g(3) * (-t17 * pkin(4) + t18 * pkin(9)) 0, 0, 0, 0, 0, 0, t6, -t5, -t46, -g(1) * (-t13 * t32 - t14 * t39) - g(2) * (-t11 * t32 - t12 * t39) - g(3) * (-t17 * t32 - t18 * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47;];
taug_reg  = t7;
