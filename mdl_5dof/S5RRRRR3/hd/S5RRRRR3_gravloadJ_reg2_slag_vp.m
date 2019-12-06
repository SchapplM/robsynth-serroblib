% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t33 = qJ(2) + qJ(3);
t30 = cos(t33);
t37 = cos(qJ(4));
t39 = cos(qJ(1));
t51 = t39 * t37;
t34 = sin(qJ(4));
t36 = sin(qJ(1));
t56 = t36 * t34;
t14 = t30 * t56 + t51;
t52 = t39 * t34;
t55 = t36 * t37;
t16 = -t30 * t52 + t55;
t28 = sin(t33);
t63 = g(3) * t28;
t68 = -g(1) * t16 + g(2) * t14 + t34 * t63;
t19 = g(1) * t39 + g(2) * t36;
t7 = -g(3) * t30 + t19 * t28;
t35 = sin(qJ(2));
t67 = pkin(1) * t35;
t65 = g(1) * t36;
t24 = t28 * pkin(5);
t25 = t30 * pkin(2);
t38 = cos(qJ(2));
t31 = t38 * pkin(1);
t61 = t28 * t36;
t60 = t28 * t39;
t26 = t37 * pkin(3) + pkin(2);
t18 = t30 * t26;
t59 = t30 * t39;
t32 = qJ(4) + qJ(5);
t27 = sin(t32);
t58 = t36 * t27;
t29 = cos(t32);
t57 = t36 * t29;
t54 = t39 * t27;
t53 = t39 * t29;
t50 = t18 + t24;
t49 = pkin(5) * t60 + t39 * t31;
t48 = t25 + t24;
t46 = -pkin(2) * t28 - t67;
t45 = -t24 - t31;
t44 = -g(2) * t39 + t65;
t43 = -t26 * t28 - t67;
t42 = t44 * t38;
t40 = -g(3) * t38 + t19 * t35;
t22 = pkin(5) * t59;
t20 = t36 * t30 * pkin(5);
t17 = t30 * t51 + t56;
t15 = -t30 * t55 + t52;
t13 = t44 * t28;
t12 = t30 * t53 + t58;
t11 = -t30 * t54 + t57;
t10 = -t30 * t57 + t54;
t9 = t30 * t58 + t53;
t8 = t19 * t30 + t63;
t6 = t7 * t37;
t5 = t7 * t34;
t4 = t7 * t29;
t3 = t7 * t27;
t2 = g(1) * t12 - g(2) * t10 + t29 * t63;
t1 = -g(1) * t11 + g(2) * t9 + t27 * t63;
t21 = [0, 0, 0, 0, 0, 0, t44, t19, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t44 * t35, -t19, 0, 0, 0, 0, 0, 0, 0, t44 * t30, -t13, -t19, pkin(1) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -g(2) * (pkin(2) * t59 + t49) - (t45 - t25) * t65, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -g(1) * pkin(3) * t52 - g(2) * (t26 * t59 + t49) + (-g(1) * (t45 - t18) - g(2) * pkin(3) * t34) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(3) * t35 + t19 * t38, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t40 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t46 * t39 + t22) - g(2) * (t46 * t36 + t20) - g(3) * (t31 + t48), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t43 * t39 + t22) - g(2) * (t43 * t36 + t20) - g(3) * (t31 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-pkin(2) * t60 + t22) - g(2) * (-pkin(2) * t61 + t20) - g(3) * t48, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (-t26 * t60 + t22) - g(2) * (-t26 * t61 + t20) - g(3) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, g(1) * t17 - g(2) * t15 + t37 * t63, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t68 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t21;
