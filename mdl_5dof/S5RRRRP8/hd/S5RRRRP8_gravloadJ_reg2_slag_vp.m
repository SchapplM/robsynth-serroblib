% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t21 = g(1) * t37 + g(2) * t34;
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t50 = t37 * t35;
t36 = cos(qJ(2));
t56 = t34 * t36;
t13 = t32 * t56 + t50;
t51 = t37 * t32;
t15 = t34 * t35 - t36 * t51;
t33 = sin(qJ(2));
t61 = g(3) * t33;
t69 = -g(1) * t15 + g(2) * t13 + t32 * t61;
t31 = qJ(3) + qJ(4);
t24 = sin(t31);
t25 = cos(t31);
t52 = t37 * t25;
t7 = t24 * t56 + t52;
t53 = t37 * t24;
t9 = t34 * t25 - t36 * t53;
t1 = -g(1) * t9 + g(2) * t7 + t24 * t61;
t11 = -g(3) * t36 + t21 * t33;
t38 = -pkin(8) - pkin(7);
t65 = g(1) * t34;
t59 = t32 * pkin(3);
t30 = -qJ(5) + t38;
t58 = t33 * t30;
t57 = t33 * t38;
t55 = t36 * t37;
t19 = pkin(4) * t24 + t59;
t54 = t37 * t19;
t27 = t35 * pkin(3);
t20 = pkin(4) * t25 + t27;
t49 = t37 * pkin(1) + t34 * pkin(6);
t46 = t36 * pkin(2) + t33 * pkin(7);
t44 = -g(2) * t37 + t65;
t18 = pkin(2) + t20;
t43 = t36 * t18 - t58;
t23 = t27 + pkin(2);
t41 = t36 * t23 - t57;
t28 = t37 * pkin(6);
t17 = t44 * t33;
t16 = t34 * t32 + t36 * t50;
t14 = -t35 * t56 + t51;
t12 = t21 * t36 + t61;
t10 = t34 * t24 + t36 * t52;
t8 = -t25 * t56 + t53;
t6 = t11 * t25;
t5 = t11 * t24;
t4 = -g(1) * t8 - g(2) * t10;
t3 = -g(1) * t7 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t8 + t25 * t61;
t22 = [0, 0, 0, 0, 0, 0, t44, t21, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t36, -t17, -t21, -g(1) * (-t34 * pkin(1) + t28) - g(2) * t49, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t17, -g(1) * t28 - g(2) * (t46 * t37 + t49) - (-pkin(1) - t46) * t65, 0, 0, 0, 0, 0, 0, t4, t3, t17, -g(1) * (pkin(3) * t51 + t28) - g(2) * (t23 * t55 - t37 * t57 + t49) + (-g(1) * (-pkin(1) - t41) - g(2) * t59) * t34, 0, 0, 0, 0, 0, 0, t4, t3, t17, -g(1) * (t28 + t54) - g(2) * (t18 * t55 - t37 * t58 + t49) + (-g(1) * (-pkin(1) - t43) - g(2) * t19) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t35, -t11 * t32, -t12, -g(3) * t46 + t21 * (pkin(2) * t33 - pkin(7) * t36), 0, 0, 0, 0, 0, 0, t6, -t5, -t12, -g(3) * t41 + t21 * (t23 * t33 + t36 * t38), 0, 0, 0, 0, 0, 0, t6, -t5, -t12, -g(3) * t43 + t21 * (t18 * t33 + t30 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, g(1) * t16 - g(2) * t14 + t35 * t61, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t69 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t34 * t20 - t36 * t54) - g(2) * (-t19 * t56 - t37 * t20) + t19 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11;];
taug_reg = t22;
