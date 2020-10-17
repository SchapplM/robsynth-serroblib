% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:58:56
% EndTime: 2019-05-07 09:58:58
% DurationCPUTime: 0.49s
% Computational Cost: add. (530->101), mult. (466->133), div. (0->0), fcn. (452->12), ass. (0->74)
t46 = qJ(2) + qJ(3);
t37 = pkin(11) + t46;
t32 = sin(t37);
t33 = cos(t37);
t50 = cos(qJ(5));
t35 = t50 * pkin(5) + pkin(4);
t53 = -pkin(10) - pkin(9);
t59 = -t32 * t53 + t33 * t35;
t60 = t33 * pkin(4) + t32 * pkin(9);
t49 = sin(qJ(1));
t52 = cos(qJ(1));
t29 = g(1) * t52 + g(2) * t49;
t67 = t52 * t50;
t47 = sin(qJ(5));
t72 = t49 * t47;
t16 = t33 * t72 + t67;
t68 = t52 * t47;
t71 = t49 * t50;
t18 = -t33 * t68 + t71;
t77 = g(3) * t32;
t86 = -g(1) * t18 + g(2) * t16 + t47 * t77;
t7 = -g(3) * t33 + t29 * t32;
t54 = -pkin(8) - pkin(7);
t39 = sin(t46);
t85 = pkin(3) * t39;
t84 = pkin(4) * t32;
t83 = pkin(9) * t33;
t41 = cos(t46);
t34 = pkin(3) * t41;
t51 = cos(qJ(2));
t42 = t51 * pkin(2);
t66 = t34 + t42;
t24 = pkin(1) + t66;
t20 = t52 * t24;
t79 = g(2) * t20;
t45 = qJ(5) + qJ(6);
t38 = sin(t45);
t74 = t49 * t38;
t40 = cos(t45);
t73 = t49 * t40;
t70 = t52 * t38;
t69 = t52 * t40;
t64 = t34 + t60;
t44 = -qJ(4) + t54;
t63 = pkin(5) * t47 - t44;
t62 = t34 + t59;
t61 = -t84 - t85;
t28 = g(1) * t49 - g(2) * t52;
t58 = -t32 * t35 - t33 * t53;
t9 = -g(3) * t41 + t29 * t39;
t48 = sin(qJ(2));
t55 = -g(3) * t51 + t29 * t48;
t36 = t42 + pkin(1);
t27 = t52 * t83;
t26 = t49 * t83;
t25 = -t48 * pkin(2) - t85;
t22 = t52 * t25;
t21 = t49 * t25;
t19 = t33 * t67 + t72;
t17 = -t33 * t71 + t68;
t15 = t28 * t32;
t14 = t33 * t69 + t74;
t13 = -t33 * t70 + t73;
t12 = -t33 * t73 + t70;
t11 = t33 * t74 + t69;
t10 = g(3) * t39 + t29 * t41;
t8 = t29 * t33 + t77;
t6 = t7 * t50;
t5 = t7 * t47;
t4 = t7 * t40;
t3 = t7 * t38;
t2 = g(1) * t14 - g(2) * t12 + t40 * t77;
t1 = -g(1) * t13 + g(2) * t11 + t38 * t77;
t23 = [0, 0, 0, 0, 0, 0, t28, t29, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t51, -t28 * t48, -t29, -g(1) * (-t49 * pkin(1) + t52 * pkin(7)) - g(2) * (t52 * pkin(1) + t49 * pkin(7)) 0, 0, 0, 0, 0, 0, t28 * t41, -t28 * t39, -t29, -g(1) * (-t49 * t36 - t52 * t54) - g(2) * (t52 * t36 - t49 * t54) 0, 0, 0, 0, 0, 0, t28 * t33, -t15, -t29, -g(1) * (-t49 * t24 - t52 * t44) - g(2) * (-t49 * t44 + t20) 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, -g(1) * t16 - g(2) * t18, t15, -t79 + (g(1) * t44 - g(2) * t60) * t52 + (-g(1) * (-t24 - t60) + g(2) * t44) * t49, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t15, -t79 + (-g(1) * t63 - g(2) * t59) * t52 + (-g(1) * (-t24 - t59) - g(2) * t63) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, g(3) * t48 + t29 * t51, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t55 * pkin(2), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(3) * t66 - t25 * t29, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-t52 * t84 + t22 + t27) - g(2) * (-t49 * t84 + t21 + t26) - g(3) * (t42 + t64) 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t52 * t58 + t22) - g(2) * (t49 * t58 + t21) - g(3) * (t42 + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t9 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t52 * t61 + t27) - g(2) * (t49 * t61 + t26) - g(3) * t64, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t62 + t29 * (-t58 + t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, g(1) * t19 - g(2) * t17 + t50 * t77, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t86 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t23;
