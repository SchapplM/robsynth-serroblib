% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t36 = qJ(2) + pkin(10);
t32 = qJ(4) + t36;
t27 = sin(t32);
t28 = cos(t32);
t54 = t28 * pkin(4) + t27 * qJ(5);
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t21 = g(1) * t43 + g(2) * t40;
t67 = t21 * t27;
t4 = g(3) * t27 + t21 * t28;
t65 = pkin(4) * t27;
t61 = g(3) * t28;
t23 = t28 * pkin(9);
t37 = -qJ(3) - pkin(7);
t35 = -pkin(8) + t37;
t60 = pkin(5) - t35;
t38 = sin(qJ(6));
t59 = t40 * t38;
t41 = cos(qJ(6));
t58 = t40 * t41;
t57 = t43 * t35;
t56 = t43 * t38;
t55 = t43 * t41;
t31 = cos(t36);
t42 = cos(qJ(2));
t33 = t42 * pkin(2);
t53 = pkin(3) * t31 + t33;
t52 = qJ(5) * t28;
t51 = t53 + t54;
t14 = pkin(1) + t53;
t11 = t43 * t14;
t50 = g(2) * (t54 * t43 + t11);
t16 = t40 * t52;
t49 = -t40 * t65 + t16;
t18 = t43 * t52;
t48 = -t43 * t65 + t18;
t20 = g(1) * t40 - g(2) * t43;
t47 = -t14 - t54;
t39 = sin(qJ(2));
t45 = -g(3) * t42 + t21 * t39;
t44 = (pkin(4) + pkin(9)) * t67;
t30 = sin(t36);
t29 = t33 + pkin(1);
t15 = -t39 * pkin(2) - pkin(3) * t30;
t13 = t43 * t15;
t12 = t40 * t15;
t10 = -t27 * t59 + t55;
t9 = t27 * t58 + t56;
t8 = t27 * t56 + t58;
t7 = t27 * t55 - t59;
t6 = t20 * t28;
t5 = t20 * t27;
t3 = -t61 + t67;
t2 = t4 * t41;
t1 = t4 * t38;
t17 = [0, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t42, -t20 * t39, -t21, -g(1) * (-t40 * pkin(1) + t43 * pkin(7)) - g(2) * (t43 * pkin(1) + t40 * pkin(7)) 0, 0, 0, 0, 0, 0, t20 * t31, -t20 * t30, -t21, -g(1) * (-t40 * t29 - t43 * t37) - g(2) * (t43 * t29 - t40 * t37) 0, 0, 0, 0, 0, 0, t6, -t5, -t21, -g(1) * (-t40 * t14 - t57) - g(2) * (-t40 * t35 + t11) 0, 0, 0, 0, 0, 0, -t21, -t6, t5, g(1) * t57 - t50 + (-g(1) * t47 + g(2) * t35) * t40, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t6, -t50 + (-g(1) * t60 - g(2) * t23) * t43 + (-g(1) * (t47 - t23) - g(2) * t60) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(3) * t39 + t21 * t42, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t31 + t21 * t30, g(3) * t30 + t21 * t31, 0, t45 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t53 - t21 * t15, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t13 + t48) - g(2) * (t12 + t49) - g(3) * t51, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (t13 + t18) - g(2) * (t12 + t16) - g(3) * (t23 + t51) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * t48 - g(2) * t49 - g(3) * t54, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * t18 - g(2) * t16 - g(3) * (t23 + t54) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t41 * t61, g(1) * t8 - g(2) * t10 - t38 * t61, 0, 0;];
taug_reg  = t17;
