% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 08:57:19
% EndTime: 2019-05-08 08:57:22
% DurationCPUTime: 0.64s
% Computational Cost: add. (621->113), mult. (607->156), div. (0->0), fcn. (603->12), ass. (0->85)
t50 = qJ(4) + qJ(5);
t44 = cos(t50);
t55 = cos(qJ(4));
t47 = t55 * pkin(4);
t29 = pkin(5) * t44 + t47;
t27 = pkin(3) + t29;
t51 = qJ(2) + qJ(3);
t43 = sin(t51);
t45 = cos(t51);
t58 = -pkin(10) - pkin(9);
t49 = -pkin(11) + t58;
t109 = t45 * t27 - t43 * t49;
t40 = t47 + pkin(3);
t108 = t45 * t40 - t43 * t58;
t107 = t45 * pkin(3) + t43 * pkin(9);
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t31 = g(1) * t57 + g(2) * t54;
t78 = t57 * t55;
t52 = sin(qJ(4));
t85 = t54 * t52;
t22 = t45 * t85 + t78;
t79 = t57 * t52;
t84 = t54 * t55;
t24 = -t45 * t79 + t84;
t96 = g(3) * t43;
t106 = -g(1) * t24 + g(2) * t22 + t52 * t96;
t80 = t57 * t44;
t42 = sin(t50);
t87 = t54 * t42;
t17 = t45 * t87 + t80;
t81 = t57 * t42;
t86 = t54 * t44;
t19 = -t45 * t81 + t86;
t3 = -g(1) * t19 + g(2) * t17 + t42 * t96;
t11 = -g(3) * t45 + t31 * t43;
t53 = sin(qJ(2));
t105 = pkin(2) * t53;
t104 = pkin(3) * t43;
t56 = cos(qJ(2));
t48 = t56 * pkin(2);
t41 = t48 + pkin(1);
t32 = t57 * t41;
t98 = g(2) * t32;
t94 = t52 * pkin(4);
t91 = t45 * t54;
t90 = t45 * t57;
t46 = qJ(6) + t50;
t38 = sin(t46);
t89 = t54 * t38;
t39 = cos(t46);
t88 = t54 * t39;
t83 = t57 * t38;
t82 = t57 * t39;
t28 = pkin(5) * t42 + t94;
t59 = -pkin(8) - pkin(7);
t77 = t28 - t59;
t73 = -t59 + t94;
t70 = -t104 - t105;
t68 = g(1) * t54 - g(2) * t57;
t66 = t27 * t43 + t45 * t49;
t64 = t40 * t43 + t45 * t58;
t60 = -g(3) * t56 + t31 * t53;
t34 = pkin(9) * t90;
t33 = pkin(9) * t91;
t25 = t45 * t78 + t85;
t23 = -t45 * t84 + t79;
t21 = t68 * t43;
t20 = t45 * t80 + t87;
t18 = -t45 * t86 + t81;
t16 = t45 * t82 + t89;
t15 = -t45 * t83 + t88;
t14 = -t45 * t88 + t83;
t13 = t45 * t89 + t82;
t12 = t31 * t45 + t96;
t10 = t11 * t55;
t9 = t11 * t52;
t8 = t11 * t44;
t7 = t11 * t42;
t6 = t11 * t39;
t5 = t11 * t38;
t4 = g(1) * t20 - g(2) * t18 + t44 * t96;
t2 = g(1) * t16 - g(2) * t14 + t39 * t96;
t1 = -g(1) * t15 + g(2) * t13 + t38 * t96;
t26 = [0, 0, 0, 0, 0, 0, t68, t31, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t56, -t68 * t53, -t31, -g(1) * (-t54 * pkin(1) + t57 * pkin(7)) - g(2) * (t57 * pkin(1) + t54 * pkin(7)) 0, 0, 0, 0, 0, 0, t68 * t45, -t21, -t31, -g(1) * (-t54 * t41 - t57 * t59) - g(2) * (-t54 * t59 + t32) 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t25, -g(1) * t22 - g(2) * t24, t21, -t98 + (g(1) * t59 - g(2) * t107) * t57 + (-g(1) * (-t41 - t107) + g(2) * t59) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t20, -g(1) * t17 - g(2) * t19, t21, -t98 + (-g(1) * t73 - g(2) * t108) * t57 + (-g(1) * (-t41 - t108) - g(2) * t73) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t21, -t98 + (-g(1) * t77 - g(2) * t109) * t57 + (-g(1) * (-t41 - t109) - g(2) * t77) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, g(3) * t53 + t31 * t56, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t60 * pkin(2), 0, 0, 0, 0, 0, 0, t10, -t9, -t12, -g(1) * (t70 * t57 + t34) - g(2) * (t70 * t54 + t33) - g(3) * (t48 + t107) 0, 0, 0, 0, 0, 0, t8, -t7, -t12, -g(3) * (t48 + t108) + t31 * (t64 + t105) 0, 0, 0, 0, 0, 0, t6, -t5, -t12, -g(3) * (t48 + t109) + t31 * (t66 + t105); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t12, -g(1) * (-t57 * t104 + t34) - g(2) * (-t54 * t104 + t33) - g(3) * t107, 0, 0, 0, 0, 0, 0, t8, -t7, -t12, -g(3) * t108 + t31 * t64, 0, 0, 0, 0, 0, 0, t6, -t5, -t12, -g(3) * t109 + t31 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, g(1) * t25 - g(2) * t23 + t55 * t96, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t106 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t28 * t90 + t54 * t29) - g(2) * (-t28 * t91 - t57 * t29) + t28 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t26;
