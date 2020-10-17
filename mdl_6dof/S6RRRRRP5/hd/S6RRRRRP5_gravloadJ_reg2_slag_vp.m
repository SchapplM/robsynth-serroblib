% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:03:11
% EndTime: 2019-05-08 05:03:14
% DurationCPUTime: 0.62s
% Computational Cost: add. (590->111), mult. (645->161), div. (0->0), fcn. (657->10), ass. (0->76)
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t32 = g(1) * t53 + g(2) * t50;
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t69 = t53 * t51;
t52 = cos(qJ(2));
t78 = t50 * t52;
t22 = t48 * t78 + t69;
t70 = t53 * t48;
t24 = t50 * t51 - t52 * t70;
t49 = sin(qJ(2));
t83 = g(3) * t49;
t94 = -g(1) * t24 + g(2) * t22 + t48 * t83;
t47 = qJ(3) + qJ(4);
t38 = sin(t47);
t39 = cos(t47);
t71 = t53 * t39;
t14 = t38 * t78 + t71;
t72 = t53 * t38;
t16 = t50 * t39 - t52 * t72;
t3 = -g(1) * t16 + g(2) * t14 + t38 * t83;
t41 = qJ(5) + t47;
t36 = cos(t41);
t35 = sin(t41);
t74 = t53 * t35;
t11 = t50 * t36 - t52 * t74;
t73 = t53 * t36;
t9 = t35 * t78 + t73;
t1 = -g(1) * t11 + g(2) * t9 + t35 * t83;
t18 = -g(3) * t52 + t32 * t49;
t54 = -pkin(9) - pkin(8);
t92 = pkin(4) * t38;
t88 = g(1) * t50;
t81 = t48 * pkin(3);
t80 = t49 * t53;
t79 = t49 * t54;
t77 = t52 * t53;
t27 = -pkin(5) * t35 - t92;
t20 = -t27 + t81;
t76 = t53 * t20;
t30 = t81 + t92;
t75 = t53 * t30;
t34 = pkin(4) * t39;
t43 = t51 * pkin(3);
t31 = t34 + t43;
t68 = t53 * pkin(1) + t50 * pkin(7);
t46 = -pkin(10) + t54;
t33 = pkin(5) * t36;
t21 = t33 + t31;
t64 = t52 * pkin(2) + t49 * pkin(8);
t62 = -g(2) * t53 + t88;
t13 = pkin(2) + t21;
t40 = -qJ(6) + t46;
t61 = t52 * t13 - t49 * t40;
t29 = pkin(2) + t31;
t59 = t52 * t29 - t49 * t46;
t37 = t43 + pkin(2);
t57 = t52 * t37 - t79;
t44 = t53 * pkin(7);
t28 = t33 + t34;
t26 = t62 * t49;
t25 = t50 * t48 + t52 * t69;
t23 = -t51 * t78 + t70;
t19 = t32 * t52 + t83;
t17 = t50 * t38 + t52 * t71;
t15 = -t39 * t78 + t72;
t12 = t50 * t35 + t52 * t73;
t10 = -t36 * t78 + t74;
t8 = t18 * t36;
t7 = t18 * t35;
t6 = -g(1) * t10 - g(2) * t12;
t5 = -g(1) * t9 - g(2) * t11;
t4 = g(1) * t17 - g(2) * t15 + t39 * t83;
t2 = g(1) * t12 - g(2) * t10 + t36 * t83;
t42 = [0, 0, 0, 0, 0, 0, t62, t32, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t52, -t26, -t32, -g(1) * (-t50 * pkin(1) + t44) - g(2) * t68, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t25, -g(1) * t22 - g(2) * t24, t26, -g(1) * t44 - g(2) * (t64 * t53 + t68) - (-pkin(1) - t64) * t88, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t26, -g(1) * (pkin(3) * t70 + t44) - g(2) * (t37 * t77 - t53 * t79 + t68) + (-g(1) * (-pkin(1) - t57) - g(2) * t81) * t50, 0, 0, 0, 0, 0, 0, t6, t5, t26, -g(1) * (t44 + t75) - g(2) * (t29 * t77 - t46 * t80 + t68) + (-g(1) * (-pkin(1) - t59) - g(2) * t30) * t50, 0, 0, 0, 0, 0, 0, t6, t5, t26, -g(1) * (t44 + t76) - g(2) * (t13 * t77 - t40 * t80 + t68) + (-g(1) * (-pkin(1) - t61) - g(2) * t20) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t51, -t18 * t48, -t19, -g(3) * t64 + t32 * (pkin(2) * t49 - pkin(8) * t52) 0, 0, 0, 0, 0, 0, t18 * t39, -t18 * t38, -t19, -g(3) * t57 + t32 * (t37 * t49 + t52 * t54) 0, 0, 0, 0, 0, 0, t8, -t7, -t19, -g(3) * t59 + t32 * (t29 * t49 + t46 * t52) 0, 0, 0, 0, 0, 0, t8, -t7, -t19, -g(3) * t61 + t32 * (t13 * t49 + t40 * t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, g(1) * t25 - g(2) * t23 + t51 * t83, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t94 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t50 * t31 - t52 * t75) - g(2) * (-t30 * t78 - t53 * t31) + t30 * t83, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t50 * t21 - t52 * t76) - g(2) * (-t20 * t78 - t53 * t21) + t20 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t27 * t77 + t50 * t28) - g(2) * (t27 * t78 - t53 * t28) - t27 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18;];
taug_reg  = t42;
