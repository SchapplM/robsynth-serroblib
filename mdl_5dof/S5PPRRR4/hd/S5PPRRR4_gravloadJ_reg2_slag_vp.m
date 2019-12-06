% Calculate inertial parameters regressor of gravitation load for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t58 = sin(pkin(11));
t59 = sin(pkin(10));
t45 = t59 * t58;
t62 = cos(pkin(11));
t63 = cos(pkin(10));
t52 = t63 * t62;
t65 = cos(pkin(5));
t36 = -t65 * t52 + t45;
t60 = sin(pkin(6));
t61 = sin(pkin(5));
t49 = t61 * t60;
t64 = cos(pkin(6));
t71 = t36 * t64 + t63 * t49;
t47 = t59 * t62;
t50 = t63 * t58;
t37 = t65 * t47 + t50;
t46 = t59 * t61;
t70 = t37 * t64 - t60 * t46;
t69 = t62 * t64 * t61 + t65 * t60;
t68 = cos(qJ(3));
t29 = sin(qJ(5));
t33 = cos(qJ(4));
t67 = t29 * t33;
t32 = cos(qJ(5));
t66 = t32 * t33;
t23 = t65 * t50 + t47;
t31 = sin(qJ(3));
t8 = t23 * t31 + t71 * t68;
t9 = t23 * t68 - t71 * t31;
t57 = -t8 * pkin(3) + t9 * pkin(8);
t24 = -t65 * t45 + t52;
t10 = t24 * t31 + t70 * t68;
t11 = t24 * t68 - t70 * t31;
t56 = -t10 * pkin(3) + t11 * pkin(8);
t48 = t61 * t58;
t15 = t31 * t48 - t69 * t68;
t16 = t69 * t31 + t68 * t48;
t55 = -t15 * pkin(3) + t16 * pkin(8);
t30 = sin(qJ(4));
t54 = -pkin(4) * t33 - pkin(9) * t30;
t51 = t63 * t61;
t22 = -t62 * t49 + t65 * t64;
t12 = -t16 * t30 + t22 * t33;
t17 = t36 * t60 - t64 * t51;
t2 = t17 * t33 - t9 * t30;
t18 = t37 * t60 + t64 * t46;
t4 = -t11 * t30 + t18 * t33;
t44 = g(1) * t4 + g(2) * t2 + g(3) * t12;
t13 = t16 * t33 + t22 * t30;
t3 = t17 * t30 + t9 * t33;
t5 = t11 * t33 + t18 * t30;
t43 = g(1) * t5 + g(2) * t3 + g(3) * t13;
t42 = g(1) * t10 + g(2) * t8 + g(3) * t15;
t41 = g(1) * t11 + g(2) * t9 + g(3) * t16;
t21 = -g(1) * t46 + g(2) * t51 - g(3) * t65;
t1 = t42 * t30;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t41, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t33, -t1, -t41, -g(1) * t56 - g(2) * t57 - g(3) * t55, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t66 + t11 * t29) - g(2) * (t9 * t29 - t8 * t66) - g(3) * (-t15 * t66 + t16 * t29), -g(1) * (t10 * t67 + t11 * t32) - g(2) * (t9 * t32 + t8 * t67) - g(3) * (t15 * t67 + t16 * t32), t1, -g(1) * (t54 * t10 + t56) - g(2) * (t54 * t8 + t57) - g(3) * (t54 * t15 + t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t32, t44 * t29, -t43, -g(1) * (t4 * pkin(4) + t5 * pkin(9)) - g(2) * (t2 * pkin(4) + t3 * pkin(9)) - g(3) * (t12 * pkin(4) + t13 * pkin(9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t32 - t5 * t29) - g(2) * (-t3 * t29 + t8 * t32) - g(3) * (-t13 * t29 + t15 * t32), -g(1) * (-t10 * t29 - t5 * t32) - g(2) * (-t8 * t29 - t3 * t32) - g(3) * (-t13 * t32 - t15 * t29), 0, 0;];
taug_reg = t6;
