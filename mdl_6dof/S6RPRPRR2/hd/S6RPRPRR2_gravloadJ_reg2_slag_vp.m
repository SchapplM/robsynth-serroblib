% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = qJ(3) + pkin(11);
t19 = sin(t27);
t21 = cos(t27);
t69 = t21 * pkin(4) + t19 * pkin(8);
t28 = qJ(1) + pkin(10);
t20 = sin(t28);
t22 = cos(t28);
t14 = g(1) * t22 + g(2) * t20;
t31 = sin(qJ(5));
t51 = t22 * t31;
t34 = cos(qJ(5));
t54 = t20 * t34;
t11 = -t21 * t51 + t54;
t62 = g(3) * t19;
t50 = t22 * t34;
t55 = t20 * t31;
t9 = t21 * t55 + t50;
t68 = -g(1) * t11 + g(2) * t9 + t31 * t62;
t40 = -g(3) * t21 + t14 * t19;
t32 = sin(qJ(3));
t66 = pkin(3) * t32;
t33 = sin(qJ(1));
t58 = t33 * pkin(1);
t29 = qJ(5) + qJ(6);
t23 = sin(t29);
t57 = t20 * t23;
t24 = cos(t29);
t56 = t20 * t24;
t53 = t22 * t23;
t52 = t22 * t24;
t35 = cos(qJ(3));
t25 = t35 * pkin(3);
t18 = t25 + pkin(2);
t36 = cos(qJ(1));
t26 = t36 * pkin(1);
t49 = t22 * t18 + t26;
t30 = -qJ(4) - pkin(7);
t47 = pkin(5) * t31 - t30;
t13 = g(1) * t20 - g(2) * t22;
t45 = g(1) * t33 - g(2) * t36;
t44 = -t22 * t30 - t58;
t17 = t34 * pkin(5) + pkin(4);
t37 = -pkin(9) - pkin(8);
t43 = t21 * t17 - t19 * t37;
t38 = -g(3) * t35 + t14 * t32;
t12 = t21 * t50 + t55;
t10 = -t21 * t54 + t51;
t8 = t13 * t19;
t7 = t21 * t52 + t57;
t6 = -t21 * t53 + t56;
t5 = -t21 * t56 + t53;
t4 = t21 * t57 + t52;
t3 = t14 * t21 + t62;
t2 = g(1) * t7 - g(2) * t5 + t24 * t62;
t1 = -g(1) * t6 + g(2) * t4 + t23 * t62;
t15 = [0, 0, 0, 0, 0, 0, t45, g(1) * t36 + g(2) * t33, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t45 * pkin(1), 0, 0, 0, 0, 0, 0, t13 * t35, -t13 * t32, -t14, -g(1) * (-t20 * pkin(2) + t22 * pkin(7) - t58) - g(2) * (t22 * pkin(2) + t20 * pkin(7) + t26) 0, 0, 0, 0, 0, 0, t13 * t21, -t8, -t14, -g(1) * (-t20 * t18 + t44) - g(2) * (-t20 * t30 + t49) 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t8, -g(1) * t44 - g(2) * (t69 * t22 + t49) + (-g(1) * (-t18 - t69) + g(2) * t30) * t20, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, -g(1) * t4 - g(2) * t6, t8, g(1) * t58 - g(2) * t49 + (-g(1) * t47 - g(2) * t43) * t22 + (-g(1) * (-t18 - t43) - g(2) * t47) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, g(3) * t32 + t14 * t35, 0, 0, 0, 0, 0, 0, 0, 0, t40, t3, 0, t38 * pkin(3), 0, 0, 0, 0, 0, 0, t40 * t34, -t40 * t31, -t3, -g(3) * (t25 + t69) + t14 * (pkin(4) * t19 - pkin(8) * t21 + t66) 0, 0, 0, 0, 0, 0, t40 * t24, -t40 * t23, -t3, -g(3) * (t25 + t43) + t14 * (t17 * t19 + t21 * t37 + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, g(1) * t12 - g(2) * t10 + t34 * t62, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t68 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t15;
