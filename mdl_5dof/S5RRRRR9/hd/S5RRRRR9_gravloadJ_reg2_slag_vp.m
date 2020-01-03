% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t22 = g(1) * t41 + g(2) * t38;
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t55 = t41 * t39;
t40 = cos(qJ(2));
t63 = t38 * t40;
t14 = t36 * t63 + t55;
t56 = t41 * t36;
t16 = t38 * t39 - t40 * t56;
t37 = sin(qJ(2));
t68 = g(3) * t37;
t76 = -g(1) * t16 + g(2) * t14 + t36 * t68;
t35 = qJ(3) + qJ(4);
t28 = cos(t35);
t27 = sin(t35);
t58 = t41 * t27;
t11 = t38 * t28 - t40 * t58;
t57 = t41 * t28;
t9 = t27 * t63 + t57;
t3 = -g(1) * t11 + g(2) * t9 + t27 * t68;
t44 = -g(3) * t40 + t22 * t37;
t42 = -pkin(8) - pkin(7);
t72 = g(1) * t38;
t66 = t36 * pkin(3);
t34 = -pkin(9) + t42;
t65 = t37 * t34;
t64 = t37 * t42;
t62 = t40 * t41;
t20 = pkin(4) * t27 + t66;
t61 = t41 * t20;
t29 = qJ(5) + t35;
t24 = sin(t29);
t60 = t41 * t24;
t25 = cos(t29);
t59 = t41 * t25;
t31 = t39 * pkin(3);
t21 = pkin(4) * t28 + t31;
t54 = t41 * pkin(1) + t38 * pkin(6);
t51 = t40 * pkin(2) + t37 * pkin(7);
t49 = -g(2) * t41 + t72;
t19 = pkin(2) + t21;
t48 = t40 * t19 - t65;
t26 = t31 + pkin(2);
t46 = t40 * t26 - t64;
t32 = t41 * pkin(6);
t18 = t49 * t37;
t17 = t38 * t36 + t40 * t55;
t15 = -t39 * t63 + t56;
t13 = t22 * t40 + t68;
t12 = t38 * t27 + t40 * t57;
t10 = -t28 * t63 + t58;
t8 = t38 * t24 + t40 * t59;
t7 = t38 * t25 - t40 * t60;
t6 = -t25 * t63 + t60;
t5 = t24 * t63 + t59;
t4 = g(1) * t12 - g(2) * t10 + t28 * t68;
t2 = g(1) * t8 - g(2) * t6 + t25 * t68;
t1 = -g(1) * t7 + g(2) * t5 + t24 * t68;
t23 = [0, 0, 0, 0, 0, 0, t49, t22, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t40, -t18, -t22, -g(1) * (-t38 * pkin(1) + t32) - g(2) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t18, -g(1) * t32 - g(2) * (t51 * t41 + t54) - (-pkin(1) - t51) * t72, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t18, -g(1) * (pkin(3) * t56 + t32) - g(2) * (t26 * t62 - t41 * t64 + t54) + (-g(1) * (-pkin(1) - t46) - g(2) * t66) * t38, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t18, -g(1) * (t32 + t61) - g(2) * (t19 * t62 - t41 * t65 + t54) + (-g(1) * (-pkin(1) - t48) - g(2) * t20) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t13, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t39, -t44 * t36, -t13, -g(3) * t51 + t22 * (pkin(2) * t37 - pkin(7) * t40), 0, 0, 0, 0, 0, 0, t44 * t28, -t44 * t27, -t13, -g(3) * t46 + t22 * (t26 * t37 + t40 * t42), 0, 0, 0, 0, 0, 0, t44 * t25, -t44 * t24, -t13, -g(3) * t48 + t22 * (t19 * t37 + t34 * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, g(1) * t17 - g(2) * t15 + t39 * t68, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t76 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t38 * t21 - t40 * t61) - g(2) * (-t20 * t63 - t41 * t21) + t20 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t23;
