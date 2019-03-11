% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t36 = cos(qJ(1));
t33 = sin(qJ(1));
t64 = g(2) * t33;
t16 = g(1) * t36 + t64;
t32 = sin(qJ(2));
t71 = t16 * t32;
t31 = sin(qJ(5));
t34 = cos(qJ(5));
t56 = t33 * t34;
t58 = t32 * t36;
t11 = -t31 * t58 - t56;
t35 = cos(qJ(2));
t62 = g(3) * t35;
t54 = t34 * t36;
t57 = t33 * t31;
t9 = t32 * t57 - t54;
t1 = -g(1) * t11 + g(2) * t9 - t31 * t62;
t8 = g(3) * t32 + t16 * t35;
t69 = pkin(2) + pkin(3);
t27 = t36 * pkin(7);
t67 = g(1) * t27;
t66 = g(1) * t33;
t26 = t35 * pkin(2);
t25 = t35 * pkin(3);
t61 = g(2) * qJ(4);
t22 = pkin(5) * t34 + pkin(4);
t60 = t22 * t32;
t30 = -qJ(6) - pkin(8);
t59 = t30 * t35;
t55 = t33 * t35;
t53 = t35 * t36;
t23 = t32 * qJ(3);
t52 = t26 + t23;
t51 = t36 * pkin(1) + t33 * pkin(7);
t50 = qJ(3) * t35;
t49 = pkin(8) + t69;
t48 = -t30 + t69;
t46 = t25 + t52;
t45 = -pkin(1) - t23;
t44 = pkin(5) * t31 + qJ(4);
t43 = g(1) * t49;
t42 = pkin(2) * t53 + t36 * t23 + t51;
t41 = g(1) * t48;
t40 = pkin(3) * t53 + t42;
t15 = -g(2) * t36 + t66;
t39 = g(1) * (-qJ(4) * t36 + t27);
t38 = t45 - t26;
t37 = g(2) * t40;
t19 = t36 * t50;
t17 = t33 * t50;
t14 = t15 * t35;
t13 = t15 * t32;
t12 = t32 * t54 - t57;
t10 = -t31 * t36 - t32 * t56;
t7 = -t62 + t71;
t6 = t8 * t34;
t5 = t8 * t31;
t4 = -g(1) * t10 - g(2) * t12;
t3 = -g(1) * t9 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t10 - t34 * t62;
t18 = [0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t16, -g(1) * (-t33 * pkin(1) + t27) - g(2) * t51, 0, 0, 0, 0, 0, 0, t14, -t16, t13, -g(2) * t42 - t38 * t66 - t67, 0, 0, 0, 0, 0, 0, t13, -t14, t16, -t39 - t37 + (-g(1) * (t38 - t25) + t61) * t33, 0, 0, 0, 0, 0, 0, t4, t3, t14, -t39 - g(2) * (pkin(4) * t58 + pkin(8) * t53 + t40) + (-g(1) * (-pkin(4) * t32 + t45) + t61 + t35 * t43) * t33, 0, 0, 0, 0, 0, 0, t4, t3, t14, -t67 - t37 + (g(1) * t44 - g(2) * (-t59 + t60)) * t36 + (-g(1) * (t45 - t60) + g(2) * t44 + t35 * t41) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t8, -g(1) * (-pkin(2) * t58 + t19) - g(2) * (-pkin(2) * t32 * t33 + t17) - g(3) * t52, 0, 0, 0, 0, 0, 0, -t8, -t7, 0, -g(1) * t19 - g(2) * t17 - g(3) * t46 + t69 * t71, 0, 0, 0, 0, 0, 0, -t6, t5, t7, -g(1) * (pkin(4) * t53 + t19) - g(2) * (pkin(4) * t55 + t17) - g(3) * (pkin(8) * t35 + t46) + (-g(3) * pkin(4) + t36 * t43 + t49 * t64) * t32, 0, 0, 0, 0, 0, 0, -t6, t5, t7, -g(1) * (t22 * t53 + t19) - g(2) * (t22 * t55 + t17) - g(3) * (t46 - t59) + (-g(3) * t22 + t36 * t41 + t48 * t64) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg  = t18;
