% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = qJ(1) + qJ(2);
t26 = sin(t30);
t31 = sin(pkin(9));
t49 = t31 * (-pkin(8) - pkin(7));
t28 = cos(t30);
t33 = sin(qJ(4));
t50 = t28 * t33;
t61 = pkin(4) * t50 + t26 * t49;
t35 = cos(qJ(4));
t32 = cos(pkin(9));
t48 = t32 * t33;
t11 = t26 * t48 + t28 * t35;
t13 = -t26 * t35 + t28 * t48;
t59 = g(1) * t31;
t60 = -g(2) * t11 + g(3) * t13 + t33 * t59;
t57 = g(2) * t28;
t55 = t28 * pkin(2);
t34 = sin(qJ(1));
t54 = t34 * pkin(1);
t36 = cos(qJ(1));
t53 = t36 * pkin(1);
t52 = t26 * t32;
t51 = t28 * t32;
t47 = t32 * t35;
t22 = t28 * qJ(3);
t45 = t22 - t54;
t44 = -t26 * pkin(2) + t22;
t18 = g(3) * t26 + t57;
t43 = g(2) * t36 + g(3) * t34;
t42 = -t26 * qJ(3) - t55;
t41 = pkin(3) * t32 + pkin(7) * t31 + pkin(2);
t24 = t35 * pkin(4) + pkin(3);
t40 = -t24 * t51 + t28 * t49 - t55;
t39 = (-g(2) * (-pkin(4) * t33 - qJ(3)) - g(3) * (-t24 * t32 - pkin(2))) * t26;
t38 = (g(2) * qJ(3) + g(3) * t41) * t26 + t41 * t57;
t29 = qJ(4) + qJ(5);
t27 = cos(t29);
t25 = sin(t29);
t17 = g(2) * t26 - g(3) * t28;
t16 = t18 * t32;
t15 = t18 * t31;
t14 = -t26 * t33 - t28 * t47;
t12 = t26 * t47 - t50;
t10 = -t26 * t25 - t27 * t51;
t9 = t25 * t51 - t26 * t27;
t8 = -t28 * t25 + t27 * t52;
t7 = t25 * t52 + t28 * t27;
t6 = -g(2) * t14 + g(3) * t12;
t5 = -g(2) * t13 - g(3) * t11;
t4 = -g(2) * t10 + g(3) * t8;
t3 = -g(2) * t9 - g(3) * t7;
t2 = -g(2) * t8 - g(3) * t10 + t27 * t59;
t1 = -g(2) * t7 + g(3) * t9 + t25 * t59;
t19 = [0, 0, 0, 0, 0, 0, t43, -g(2) * t34 + g(3) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t43 * pkin(1), 0, 0, 0, 0, 0, 0, t16, -t15, t17, -g(2) * (t42 - t53) - g(3) * (t44 - t54), 0, 0, 0, 0, 0, 0, t6, t5, t15, g(2) * t53 - g(3) * t45 + t38, 0, 0, 0, 0, 0, 0, t4, t3, t15, -g(2) * (t40 - t53) - g(3) * (t45 + t61) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t17, -g(2) * t42 - g(3) * t44, 0, 0, 0, 0, 0, 0, t6, t5, t15, -g(3) * t22 + t38, 0, 0, 0, 0, 0, 0, t4, t3, t15, -g(2) * t40 - g(3) * (t22 + t61) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -g(2) * t12 - g(3) * t14 + t35 * t59, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t60 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t19;
