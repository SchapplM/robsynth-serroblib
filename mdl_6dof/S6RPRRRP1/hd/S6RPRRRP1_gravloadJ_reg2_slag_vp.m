% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t32 = qJ(1) + pkin(10);
t26 = sin(t32);
t27 = cos(t32);
t12 = g(1) * t27 + g(2) * t26;
t33 = qJ(3) + qJ(4);
t28 = sin(t33);
t45 = t12 * t28;
t29 = cos(t33);
t55 = t29 * pkin(4) + t28 * pkin(9);
t35 = sin(qJ(3));
t66 = pkin(3) * t35;
t65 = pkin(4) * t28;
t62 = g(3) * t28;
t36 = sin(qJ(1));
t61 = t36 * pkin(1);
t60 = t27 * t28;
t59 = t27 * t29;
t34 = sin(qJ(5));
t58 = t29 * t34;
t37 = cos(qJ(5));
t57 = t29 * t37;
t38 = cos(qJ(3));
t30 = t38 * pkin(3);
t25 = t30 + pkin(2);
t39 = cos(qJ(1));
t31 = t39 * pkin(1);
t56 = t27 * t25 + t31;
t54 = qJ(6) * t34;
t53 = pkin(4) * t59 + pkin(9) * t60 + t56;
t52 = pkin(5) * t57 + t29 * t54 + t55;
t10 = -t26 * t37 + t27 * t58;
t8 = t26 * t58 + t27 * t37;
t51 = g(1) * t8 - g(2) * t10;
t50 = -t65 - t66;
t49 = g(1) * t26 - g(2) * t27;
t48 = g(1) * t36 - g(2) * t39;
t40 = -pkin(8) - pkin(7);
t47 = -t27 * t40 - t61;
t1 = g(1) * t10 + g(2) * t8 + t34 * t62;
t11 = t26 * t34 + t27 * t57;
t9 = t26 * t57 - t27 * t34;
t44 = g(1) * t11 + g(2) * t9 + t37 * t62;
t5 = -g(3) * t29 + t45;
t43 = -g(3) * t38 + t12 * t35;
t42 = (-g(1) * (-t25 - t55) + g(2) * t40) * t26;
t41 = (pkin(5) * t37 + pkin(4) + t54) * t45;
t16 = pkin(9) * t59;
t14 = t26 * t29 * pkin(9);
t7 = t49 * t28;
t6 = t12 * t29 + t62;
t4 = t5 * t37;
t3 = -g(3) * t58 + t34 * t45;
t2 = g(1) * t9 - g(2) * t11;
t13 = [0, 0, 0, 0, 0, 0, t48, g(1) * t39 + g(2) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t49, t12, 0, t48 * pkin(1), 0, 0, 0, 0, 0, 0, t49 * t38, -t49 * t35, -t12, -g(1) * (-t26 * pkin(2) + t27 * pkin(7) - t61) - g(2) * (t27 * pkin(2) + t26 * pkin(7) + t31) 0, 0, 0, 0, 0, 0, t49 * t29, -t7, -t12, -g(1) * (-t26 * t25 + t47) - g(2) * (-t26 * t40 + t56) 0, 0, 0, 0, 0, 0, t2, -t51, t7, -g(1) * t47 - g(2) * t53 + t42, 0, 0, 0, 0, 0, 0, t2, t7, t51, -g(1) * (-t9 * pkin(5) - t8 * qJ(6) + t47) - g(2) * (t11 * pkin(5) + t10 * qJ(6) + t53) + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(3) * t35 + t12 * t38, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t43 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t50 * t27 + t16) - g(2) * (t50 * t26 + t14) - g(3) * (t30 + t55) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (-t27 * t66 + t16) - g(2) * (-t26 * t66 + t14) - g(3) * (t30 + t52) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-pkin(4) * t60 + t16) - g(2) * (-t26 * t65 + t14) - g(3) * t55, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t16 - g(2) * t14 - g(3) * t52 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t44, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t44, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t8 * pkin(5) + t9 * qJ(6)) - (-pkin(5) * t34 + qJ(6) * t37) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
