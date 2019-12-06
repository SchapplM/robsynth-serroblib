% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = sin(pkin(5));
t61 = g(3) * t27;
t25 = qJ(4) + qJ(5);
t23 = sin(t25);
t32 = cos(qJ(3));
t60 = t23 * t32;
t24 = cos(t25);
t59 = t24 * t32;
t30 = sin(qJ(2));
t58 = t27 * t30;
t57 = t27 * t32;
t33 = cos(qJ(2));
t56 = t27 * t33;
t28 = sin(qJ(4));
t55 = t28 * t30;
t54 = t28 * t32;
t31 = cos(qJ(4));
t53 = t31 * t32;
t52 = t32 * t33;
t51 = pkin(2) * t56 + pkin(7) * t58;
t50 = cos(pkin(5));
t49 = cos(pkin(10));
t48 = pkin(4) * t28 + pkin(7);
t26 = sin(pkin(10));
t40 = t50 * t49;
t12 = t26 * t33 + t30 * t40;
t11 = t26 * t30 - t33 * t40;
t9 = t11 * pkin(2);
t47 = t12 * pkin(7) - t9;
t46 = g(3) * t51;
t44 = t26 * t50;
t13 = t30 * t49 + t33 * t44;
t10 = t13 * pkin(2);
t14 = -t30 * t44 + t33 * t49;
t45 = t14 * pkin(7) - t10;
t43 = t27 * t49;
t29 = sin(qJ(3));
t42 = pkin(3) * t32 + pkin(8) * t29;
t22 = pkin(4) * t31 + pkin(3);
t34 = -pkin(9) - pkin(8);
t41 = t22 * t32 - t29 * t34;
t15 = -t29 * t58 + t32 * t50;
t5 = -t12 * t29 - t32 * t43;
t7 = -t14 * t29 + t26 * t57;
t39 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t16 = t29 * t50 + t30 * t57;
t6 = t12 * t32 - t29 * t43;
t8 = t26 * t27 * t29 + t14 * t32;
t38 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t37 = -g(1) * t13 - g(2) * t11 + g(3) * t56;
t36 = g(1) * t14 + g(2) * t12 + g(3) * t58;
t35 = -g(1) * (t13 * t31 - t28 * t8) - g(2) * (t11 * t31 - t28 * t6) - g(3) * (-t16 * t28 - t31 * t56);
t4 = t37 * t29;
t2 = -g(1) * (-t13 * t23 - t24 * t8) - g(2) * (-t11 * t23 - t24 * t6) - g(3) * (-t16 * t24 + t23 * t56);
t1 = -g(1) * (t13 * t24 - t23 * t8) - g(2) * (t11 * t24 - t23 * t6) - g(3) * (-t16 * t23 - t24 * t56);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t36, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t32, t4, -t36, -g(1) * t45 - g(2) * t47 - t46, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t53 + t14 * t28) - g(2) * (-t11 * t53 + t12 * t28) - (t31 * t52 + t55) * t61, -g(1) * (t13 * t54 + t14 * t31) - g(2) * (t11 * t54 + t12 * t31) - (-t28 * t52 + t30 * t31) * t61, -t4, -g(1) * (-t13 * t42 + t45) - g(2) * (-t11 * t42 + t47) - g(3) * (t42 * t56 + t51), 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t59 + t14 * t23) - g(2) * (-t11 * t59 + t12 * t23) - (t23 * t30 + t24 * t52) * t61, -g(1) * (t13 * t60 + t14 * t24) - g(2) * (t11 * t60 + t12 * t24) - (-t23 * t52 + t24 * t30) * t61, -t4, -g(1) * (-t13 * t41 + t14 * t48 - t10) - g(2) * (-t11 * t41 + t12 * t48 - t9) - t46 - (pkin(4) * t55 + t33 * t41) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * t31, t39 * t28, -t38, -g(1) * (pkin(3) * t7 + pkin(8) * t8) - g(2) * (pkin(3) * t5 + pkin(8) * t6) - g(3) * (pkin(3) * t15 + pkin(8) * t16), 0, 0, 0, 0, 0, 0, -t39 * t24, t39 * t23, -t38, -g(1) * (t22 * t7 - t34 * t8) - g(2) * (t22 * t5 - t34 * t6) - g(3) * (t15 * t22 - t16 * t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -g(1) * (-t13 * t28 - t31 * t8) - g(2) * (-t11 * t28 - t31 * t6) - g(3) * (-t16 * t31 + t28 * t56), 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t35 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t3;
