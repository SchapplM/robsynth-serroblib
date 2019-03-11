% Calculate inertial parameters regressor of gravitation load for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t40 = sin(pkin(11));
t48 = sin(qJ(2));
t50 = cos(qJ(2));
t65 = cos(pkin(11));
t57 = -t48 * t40 + t50 * t65;
t38 = pkin(12) + qJ(5);
t37 = cos(t38);
t47 = sin(qJ(6));
t74 = t37 * t47;
t49 = cos(qJ(6));
t73 = t37 * t49;
t41 = sin(pkin(10));
t42 = sin(pkin(6));
t72 = t41 * t42;
t71 = t41 * t48;
t44 = cos(pkin(10));
t70 = t42 * t44;
t69 = t42 * t50;
t45 = cos(pkin(6));
t68 = t45 * t48;
t67 = t45 * t50;
t64 = t44 * t67;
t25 = t57 * t42;
t29 = -t50 * t40 - t48 * t65;
t26 = t29 * t42;
t33 = pkin(2) * t69;
t43 = cos(pkin(12));
t35 = t43 * pkin(4) + pkin(3);
t46 = -pkin(8) - qJ(4);
t63 = t25 * t35 + t26 * t46 + t33;
t31 = pkin(2) * t64;
t61 = -pkin(2) * t71 + t31;
t36 = sin(t38);
t60 = pkin(5) * t37 + pkin(9) * t36;
t27 = t29 * t45;
t14 = -t44 * t27 + t41 * t57;
t15 = -t41 * t27 - t44 * t57;
t59 = -t41 * t67 - t44 * t48;
t53 = t57 * t45;
t13 = t41 * t29 + t44 * t53;
t58 = t13 * t35 - t14 * t46 + t61;
t18 = t26 * t36 + t45 * t37;
t4 = -t14 * t36 - t37 * t70;
t6 = t15 * t36 + t37 * t72;
t56 = g(1) * t6 + g(2) * t4 + g(3) * t18;
t19 = -t26 * t37 + t45 * t36;
t5 = t14 * t37 - t36 * t70;
t7 = -t15 * t37 + t36 * t72;
t55 = g(1) * t7 + g(2) * t5 + g(3) * t19;
t2 = g(1) * t15 - g(2) * t14 + g(3) * t26;
t16 = t44 * t29 - t41 * t53;
t3 = g(1) * t16 + g(2) * t13 + g(3) * t25;
t54 = t59 * pkin(2);
t52 = t15 * t46 + t16 * t35 + t54;
t51 = -g(1) * t59 - g(3) * t69;
t24 = -g(3) * t45 + (-g(1) * t41 + g(2) * t44) * t42;
t1 = t3 * t36;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t64 - t71) + t51, -g(1) * (t41 * t68 - t44 * t50) - g(2) * (-t41 * t50 - t44 * t68) + g(3) * t42 * t48, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t2, 0, -g(2) * t31 + (g(2) * t71 + t51) * pkin(2), 0, 0, 0, 0, 0, 0, -t3 * t43, t3 * sin(pkin(12)) t2, -g(1) * (t16 * pkin(3) - t15 * qJ(4) + t54) - g(2) * (t13 * pkin(3) + qJ(4) * t14 + t61) - g(3) * (t25 * pkin(3) - t26 * qJ(4) + t33) 0, 0, 0, 0, 0, 0, -t3 * t37, t1, t2, -g(1) * t52 - g(2) * t58 - g(3) * t63, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t47 + t16 * t73) - g(2) * (t13 * t73 + t14 * t47) - g(3) * (t25 * t73 - t26 * t47) -g(1) * (-t15 * t49 - t16 * t74) - g(2) * (-t13 * t74 + t14 * t49) - g(3) * (-t25 * t74 - t26 * t49) -t1, -g(1) * (t60 * t16 + t52) - g(2) * (t60 * t13 + t58) - g(3) * (t60 * t25 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t55, 0, 0, 0, 0, 0, 0, 0, 0, -t56 * t49, t56 * t47, -t55, -g(1) * (t6 * pkin(5) + t7 * pkin(9)) - g(2) * (t4 * pkin(5) + t5 * pkin(9)) - g(3) * (t18 * pkin(5) + t19 * pkin(9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t49 - t7 * t47) - g(2) * (-t13 * t49 - t5 * t47) - g(3) * (-t19 * t47 - t25 * t49) -g(1) * (t16 * t47 - t7 * t49) - g(2) * (t13 * t47 - t5 * t49) - g(3) * (-t19 * t49 + t25 * t47) 0, 0;];
taug_reg  = t8;
