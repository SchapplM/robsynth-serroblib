% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = -qJ(4) - pkin(7);
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t33 = sin(qJ(3));
t60 = t33 * pkin(3);
t70 = t34 * t31 + t37 * t60;
t63 = g(2) * t37;
t13 = g(1) * t34 - t63;
t29 = qJ(3) + pkin(10);
t21 = sin(t29);
t32 = sin(qJ(5));
t53 = t37 * t32;
t35 = cos(qJ(5));
t56 = t34 * t35;
t11 = t21 * t53 + t56;
t22 = cos(t29);
t61 = g(3) * t22;
t52 = t37 * t35;
t57 = t34 * t32;
t9 = -t21 * t57 + t52;
t69 = -g(1) * t9 - g(2) * t11 + t32 * t61;
t41 = -g(3) * t21 + t13 * t22;
t36 = cos(qJ(3));
t67 = pkin(3) * t36;
t66 = pkin(5) * t32;
t30 = qJ(5) + qJ(6);
t23 = sin(t30);
t59 = t34 * t23;
t24 = cos(t30);
t58 = t34 * t24;
t55 = t37 * t23;
t54 = t37 * t24;
t51 = t37 * pkin(1) + t34 * qJ(2);
t49 = t34 * t60 + t51;
t26 = t37 * qJ(2);
t48 = -t34 * pkin(1) + t26;
t47 = pkin(4) * t22 + pkin(8) * t21;
t46 = t21 * pkin(4) - t22 * pkin(8);
t14 = g(1) * t37 + g(2) * t34;
t20 = t35 * pkin(5) + pkin(4);
t38 = -pkin(9) - pkin(8);
t45 = t20 * t22 - t21 * t38;
t44 = t21 * t20 + t22 * t38;
t43 = t48 + t70;
t42 = -t37 * t31 + t49;
t39 = g(3) * t33 - t13 * t36;
t17 = t34 * t67;
t12 = t21 * t52 - t57;
t10 = t21 * t56 + t53;
t8 = t14 * t22;
t7 = t21 * t54 - t59;
t6 = t21 * t55 + t58;
t5 = t21 * t58 + t55;
t4 = -t21 * t59 + t54;
t3 = t13 * t21 + t61;
t2 = g(1) * t5 - g(2) * t7 + t24 * t61;
t1 = -g(1) * t4 - g(2) * t6 + t23 * t61;
t15 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -g(1) * t48 - g(2) * t51, 0, 0, 0, 0, 0, 0, -t14 * t33, -t14 * t36, t13, -g(1) * (t26 + (-pkin(1) - pkin(7)) * t34) - g(2) * (t37 * pkin(7) + t51) 0, 0, 0, 0, 0, 0, -t14 * t21, -t8, t13, -g(1) * t43 - g(2) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t8, -g(1) * (t46 * t37 + t43) - g(2) * (t46 * t34 + t42) 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t5, g(1) * t6 - g(2) * t4, t8, -g(1) * (t26 + t70) - g(2) * t49 + (-g(1) * t44 - g(2) * (-t31 + t66)) * t37 + (-g(1) * (-pkin(1) - t66) - g(2) * t44) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, g(3) * t36 + t13 * t33, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t3, 0, t39 * pkin(3), 0, 0, 0, 0, 0, 0, -t41 * t35, t41 * t32, -t3, -g(1) * (t47 * t34 + t17) - g(3) * (-t46 - t60) - (-t47 - t67) * t63, 0, 0, 0, 0, 0, 0, -t41 * t24, t41 * t23, -t3, -g(1) * (t45 * t34 + t17) - g(3) * (-t44 - t60) - (-t45 - t67) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, g(1) * t10 - g(2) * t12 + t35 * t61, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t69 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t15;
