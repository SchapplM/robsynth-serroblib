% Calculate inertial parameters regressor of gravitation load for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = cos(pkin(11));
t27 = pkin(3) * t38 + pkin(2);
t36 = sin(pkin(6));
t42 = cos(qJ(2));
t58 = t36 * t42;
t20 = t27 * t58;
t68 = g(3) * t20;
t67 = g(3) * t36;
t41 = sin(qJ(2));
t53 = cos(pkin(10));
t54 = cos(pkin(6));
t47 = t54 * t53;
t52 = sin(pkin(10));
t17 = t41 * t47 + t52 * t42;
t34 = sin(pkin(12));
t66 = t17 * t34;
t46 = t54 * t52;
t19 = -t41 * t46 + t53 * t42;
t65 = t19 * t34;
t32 = pkin(12) + qJ(6);
t28 = sin(t32);
t33 = pkin(11) + qJ(4);
t31 = cos(t33);
t64 = t28 * t31;
t30 = cos(t32);
t63 = t30 * t31;
t62 = t31 * t34;
t37 = cos(pkin(12));
t61 = t31 * t37;
t60 = t31 * t42;
t59 = t36 * t41;
t40 = -pkin(8) - qJ(3);
t57 = t40 * t41;
t16 = t52 * t41 - t42 * t47;
t56 = -t16 * t27 - t17 * t40;
t18 = t53 * t41 + t42 * t46;
t55 = -t18 * t27 - t19 * t40;
t51 = t36 * t53;
t50 = t36 * t52;
t29 = sin(t33);
t49 = pkin(4) * t31 + qJ(5) * t29;
t26 = pkin(5) * t37 + pkin(4);
t39 = -pkin(9) - qJ(5);
t48 = t26 * t31 - t29 * t39;
t12 = t29 * t59 - t54 * t31;
t6 = t17 * t29 + t31 * t51;
t8 = t19 * t29 - t31 * t50;
t45 = g(1) * t8 + g(2) * t6 + g(3) * t12;
t13 = t54 * t29 + t31 * t59;
t7 = t17 * t31 - t29 * t51;
t9 = t19 * t31 + t29 * t50;
t44 = g(1) * t9 + g(2) * t7 + g(3) * t13;
t4 = -g(1) * t18 - g(2) * t16 + g(3) * t58;
t43 = g(1) * t19 + g(2) * t17 + g(3) * t59;
t3 = t4 * t29;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t43, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t38, t4 * sin(pkin(11)) -t43, -g(1) * (-pkin(2) * t18 + qJ(3) * t19) - g(2) * (-pkin(2) * t16 + qJ(3) * t17) - (pkin(2) * t42 + qJ(3) * t41) * t67, 0, 0, 0, 0, 0, 0, -t4 * t31, t3, -t43, -g(1) * t55 - g(2) * t56 - g(3) * (-t36 * t57 + t20) 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t61 + t65) - g(2) * (-t16 * t61 + t66) - (t34 * t41 + t37 * t60) * t67, -g(1) * (t18 * t62 + t19 * t37) - g(2) * (t16 * t62 + t17 * t37) - (-t34 * t60 + t37 * t41) * t67, -t3, -g(1) * (-t49 * t18 + t55) - g(2) * (-t49 * t16 + t56) - t68 - (t49 * t42 - t57) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t63 + t19 * t28) - g(2) * (-t16 * t63 + t17 * t28) - (t28 * t41 + t30 * t60) * t67, -g(1) * (t18 * t64 + t19 * t30) - g(2) * (t16 * t64 + t17 * t30) - (-t28 * t60 + t30 * t41) * t67, -t3, -g(1) * (pkin(5) * t65 - t48 * t18 + t55) - g(2) * (pkin(5) * t66 - t48 * t16 + t56) - t68 - (t48 * t42 + (pkin(5) * t34 - t40) * t41) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t44, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t37, -t45 * t34, -t44, -g(1) * (-pkin(4) * t8 + qJ(5) * t9) - g(2) * (-pkin(4) * t6 + qJ(5) * t7) - g(3) * (-pkin(4) * t12 + qJ(5) * t13) 0, 0, 0, 0, 0, 0, t45 * t30, -t45 * t28, -t44, -g(1) * (-t26 * t8 - t39 * t9) - g(2) * (-t26 * t6 - t39 * t7) - g(3) * (-t12 * t26 - t13 * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t18 * t30 - t28 * t9) - g(2) * (t16 * t30 - t28 * t7) - g(3) * (-t13 * t28 - t30 * t58) -g(1) * (-t18 * t28 - t30 * t9) - g(2) * (-t16 * t28 - t30 * t7) - g(3) * (-t13 * t30 + t28 * t58) 0, 0;];
taug_reg  = t1;
