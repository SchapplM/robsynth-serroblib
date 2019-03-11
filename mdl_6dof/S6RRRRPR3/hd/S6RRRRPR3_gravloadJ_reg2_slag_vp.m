% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t36 = qJ(2) + qJ(3);
t33 = qJ(4) + t36;
t28 = sin(t33);
t29 = cos(t33);
t57 = t29 * pkin(4) + t28 * qJ(5);
t39 = sin(qJ(1));
t42 = cos(qJ(1));
t22 = g(1) * t42 + g(2) * t39;
t71 = t22 * t28;
t4 = g(3) * t28 + t22 * t29;
t43 = -pkin(8) - pkin(7);
t31 = sin(t36);
t69 = pkin(3) * t31;
t68 = pkin(4) * t28;
t64 = g(3) * t29;
t24 = t29 * pkin(10);
t35 = -pkin(9) + t43;
t63 = pkin(5) - t35;
t37 = sin(qJ(6));
t62 = t39 * t37;
t40 = cos(qJ(6));
t61 = t39 * t40;
t60 = t42 * t35;
t59 = t42 * t37;
t58 = t42 * t40;
t32 = cos(t36);
t27 = pkin(3) * t32;
t41 = cos(qJ(2));
t34 = t41 * pkin(2);
t56 = t27 + t34;
t55 = qJ(5) * t29;
t54 = t27 + t57;
t53 = t34 + t54;
t16 = pkin(1) + t56;
t13 = t42 * t16;
t52 = g(2) * (t57 * t42 + t13);
t18 = t39 * t55;
t51 = -t39 * t68 + t18;
t20 = t42 * t55;
t50 = -t42 * t68 + t20;
t49 = -t68 - t69;
t48 = g(1) * t39 - g(2) * t42;
t47 = -t16 - t57;
t5 = -g(3) * t32 + t22 * t31;
t38 = sin(qJ(2));
t45 = -g(3) * t41 + t22 * t38;
t44 = (pkin(4) + pkin(10)) * t71;
t30 = t34 + pkin(1);
t17 = -t38 * pkin(2) - t69;
t15 = t42 * t17;
t14 = t39 * t17;
t12 = -t28 * t62 + t58;
t11 = t28 * t61 + t59;
t10 = t28 * t59 + t61;
t9 = t28 * t58 - t62;
t8 = t48 * t29;
t7 = t48 * t28;
t6 = g(3) * t31 + t22 * t32;
t3 = -t64 + t71;
t2 = t4 * t40;
t1 = t4 * t37;
t19 = [0, 0, 0, 0, 0, 0, t48, t22, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t41, -t48 * t38, -t22, -g(1) * (-t39 * pkin(1) + t42 * pkin(7)) - g(2) * (t42 * pkin(1) + t39 * pkin(7)) 0, 0, 0, 0, 0, 0, t48 * t32, -t48 * t31, -t22, -g(1) * (-t39 * t30 - t42 * t43) - g(2) * (t42 * t30 - t39 * t43) 0, 0, 0, 0, 0, 0, t8, -t7, -t22, -g(1) * (-t39 * t16 - t60) - g(2) * (-t39 * t35 + t13) 0, 0, 0, 0, 0, 0, -t22, -t8, t7, g(1) * t60 - t52 + (-g(1) * t47 + g(2) * t35) * t39, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t8, -t52 + (-g(1) * t63 - g(2) * t24) * t42 + (-g(1) * (t47 - t24) - g(2) * t63) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(3) * t38 + t22 * t41, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t45 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t56 - t22 * t17, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t15 + t50) - g(2) * (t14 + t51) - g(3) * t53, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (t15 + t20) - g(2) * (t14 + t18) - g(3) * (t24 + t53) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t49 * t42 + t20) - g(2) * (t49 * t39 + t18) - g(3) * t54, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (-t42 * t69 + t20) - g(2) * (-t39 * t69 + t18) - g(3) * (t24 + t54) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * t50 - g(2) * t51 - g(3) * t57, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * t20 - g(2) * t18 - g(3) * (t24 + t57) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11 + t40 * t64, g(1) * t10 - g(2) * t12 - t37 * t64, 0, 0;];
taug_reg  = t19;
