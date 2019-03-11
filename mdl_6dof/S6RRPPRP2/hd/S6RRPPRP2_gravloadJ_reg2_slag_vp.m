% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = qJ(2) + pkin(9);
t27 = sin(t31);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t67 = -g(1) * t39 - g(2) * t36;
t70 = t27 * t67;
t28 = cos(t31);
t69 = t28 * pkin(3) + t27 * qJ(4);
t33 = -qJ(3) - pkin(7);
t57 = pkin(4) - t33;
t37 = cos(qJ(5));
t51 = t39 * t37;
t34 = sin(qJ(5));
t55 = t36 * t34;
t11 = t27 * t51 - t55;
t52 = t39 * t34;
t54 = t36 * t37;
t13 = t27 * t54 + t52;
t59 = g(3) * t28;
t1 = -g(1) * t11 - g(2) * t13 + t37 * t59;
t8 = g(3) * t27 - t28 * t67;
t35 = sin(qJ(2));
t65 = pkin(2) * t35;
t58 = t28 * pkin(8);
t56 = t28 * t34;
t53 = t39 * t33;
t50 = t37 * pkin(5) + t57;
t49 = qJ(4) * t28;
t38 = cos(qJ(2));
t29 = t38 * pkin(2);
t47 = t29 + t69;
t26 = t29 + pkin(1);
t21 = t39 * t26;
t46 = g(2) * (t69 * t39 + t21);
t45 = -pkin(3) * t27 - t65;
t18 = g(1) * t36 - g(2) * t39;
t32 = -qJ(6) - pkin(8);
t44 = t27 * t34 * pkin(5) - t28 * t32;
t43 = -t26 - t69;
t41 = -g(3) * t38 - t35 * t67;
t17 = t39 * t49;
t15 = t36 * t49;
t14 = -t27 * t55 + t51;
t12 = t27 * t52 + t54;
t10 = t18 * t28;
t9 = t18 * t27;
t7 = -t59 - t70;
t6 = t8 * t37;
t5 = t8 * t34;
t4 = -g(1) * t14 - g(2) * t12;
t3 = g(1) * t13 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t14 - g(3) * t56;
t16 = [0, 0, 0, 0, 0, 0, t18, -t67, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t38, -t18 * t35, t67, -g(1) * (-t36 * pkin(1) + t39 * pkin(7)) - g(2) * (t39 * pkin(1) + t36 * pkin(7)) 0, 0, 0, 0, 0, 0, t10, -t9, t67, -g(1) * (-t36 * t26 - t53) - g(2) * (-t36 * t33 + t21) 0, 0, 0, 0, 0, 0, t67, -t10, t9, g(1) * t53 - t46 + (-g(1) * t43 + g(2) * t33) * t36, 0, 0, 0, 0, 0, 0, t4, t3, t10, -t46 + (-g(1) * t57 - g(2) * t58) * t39 + (-g(1) * (t43 - t58) - g(2) * t57) * t36, 0, 0, 0, 0, 0, 0, t4, t3, t10, -t46 + (-g(1) * t50 - g(2) * t44) * t39 + (-g(1) * (t43 - t44) - g(2) * t50) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, g(3) * t35 - t38 * t67, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t41 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (t45 * t39 + t17) - g(2) * (t45 * t36 + t15) - g(3) * t47, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * (-t39 * t65 + t17) - g(2) * (-t36 * t65 + t15) - g(3) * (t47 + t58) - (pkin(3) + pkin(8)) * t70, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * t17 - g(2) * t15 - g(3) * (t44 + t47) + t67 * (pkin(5) * t56 - t65 + (-pkin(3) + t32) * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg  = t16;
