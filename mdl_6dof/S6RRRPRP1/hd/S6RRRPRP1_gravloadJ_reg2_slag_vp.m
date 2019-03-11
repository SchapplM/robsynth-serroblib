% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t39 = qJ(2) + qJ(3);
t33 = pkin(10) + t39;
t28 = sin(t33);
t29 = cos(t33);
t44 = cos(qJ(5));
t31 = t44 * pkin(5) + pkin(4);
t40 = -qJ(6) - pkin(9);
t51 = -t28 * t40 + t29 * t31;
t52 = t29 * pkin(4) + t28 * pkin(9);
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t25 = g(1) * t46 + g(2) * t43;
t59 = t46 * t44;
t41 = sin(qJ(5));
t62 = t43 * t41;
t12 = t29 * t62 + t59;
t60 = t46 * t41;
t61 = t43 * t44;
t14 = -t29 * t60 + t61;
t65 = g(3) * t28;
t1 = -g(1) * t14 + g(2) * t12 + t41 * t65;
t7 = -g(3) * t29 + t25 * t28;
t47 = -pkin(8) - pkin(7);
t34 = sin(t39);
t73 = pkin(3) * t34;
t72 = pkin(4) * t28;
t71 = pkin(9) * t29;
t35 = cos(t39);
t30 = pkin(3) * t35;
t45 = cos(qJ(2));
t36 = t45 * pkin(2);
t58 = t30 + t36;
t20 = pkin(1) + t58;
t16 = t46 * t20;
t67 = g(2) * t16;
t56 = t30 + t52;
t38 = -qJ(4) + t47;
t55 = pkin(5) * t41 - t38;
t54 = t30 + t51;
t53 = -t72 - t73;
t24 = g(1) * t43 - g(2) * t46;
t50 = -t28 * t31 - t29 * t40;
t9 = -g(3) * t35 + t25 * t34;
t42 = sin(qJ(2));
t48 = -g(3) * t45 + t25 * t42;
t32 = t36 + pkin(1);
t23 = t46 * t71;
t22 = t43 * t71;
t21 = -t42 * pkin(2) - t73;
t18 = t46 * t21;
t17 = t43 * t21;
t15 = t29 * t59 + t62;
t13 = -t29 * t61 + t60;
t11 = t24 * t28;
t10 = g(3) * t34 + t25 * t35;
t8 = t25 * t29 + t65;
t6 = t7 * t44;
t5 = t7 * t41;
t4 = -g(1) * t13 - g(2) * t15;
t3 = -g(1) * t12 - g(2) * t14;
t2 = g(1) * t15 - g(2) * t13 + t44 * t65;
t19 = [0, 0, 0, 0, 0, 0, t24, t25, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t45, -t24 * t42, -t25, -g(1) * (-t43 * pkin(1) + t46 * pkin(7)) - g(2) * (t46 * pkin(1) + t43 * pkin(7)) 0, 0, 0, 0, 0, 0, t24 * t35, -t24 * t34, -t25, -g(1) * (-t43 * t32 - t46 * t47) - g(2) * (t46 * t32 - t43 * t47) 0, 0, 0, 0, 0, 0, t24 * t29, -t11, -t25, -g(1) * (-t43 * t20 - t46 * t38) - g(2) * (-t43 * t38 + t16) 0, 0, 0, 0, 0, 0, t4, t3, t11, -t67 + (g(1) * t38 - g(2) * t52) * t46 + (-g(1) * (-t20 - t52) + g(2) * t38) * t43, 0, 0, 0, 0, 0, 0, t4, t3, t11, -t67 + (-g(1) * t55 - g(2) * t51) * t46 + (-g(1) * (-t20 - t51) - g(2) * t55) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(3) * t42 + t25 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t48 * pkin(2), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(3) * t58 - t25 * t21, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-t46 * t72 + t18 + t23) - g(2) * (-t43 * t72 + t17 + t22) - g(3) * (t36 + t56) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t50 * t46 + t18) - g(2) * (t50 * t43 + t17) - g(3) * (t36 + t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t9 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t53 * t46 + t23) - g(2) * (t53 * t43 + t22) - g(3) * t56, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t54 + t25 * (-t50 + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t19;
