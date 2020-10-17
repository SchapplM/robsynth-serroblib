% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 16:13:53
% EndTime: 2019-05-08 16:13:57
% DurationCPUTime: 0.97s
% Computational Cost: add. (834->133), mult. (2208->256), div. (0->0), fcn. (2908->16), ass. (0->75)
t53 = sin(qJ(2));
t56 = cos(qJ(2));
t57 = cos(qJ(1));
t72 = cos(pkin(6));
t66 = t57 * t72;
t82 = sin(qJ(1));
t37 = t82 * t53 - t56 * t66;
t38 = t53 * t66 + t82 * t56;
t52 = sin(qJ(3));
t71 = cos(pkin(7));
t83 = cos(qJ(3));
t62 = t71 * t83;
t48 = sin(pkin(7));
t49 = sin(pkin(6));
t79 = t48 * t49;
t69 = t57 * t79;
t16 = t37 * t62 + t38 * t52 + t83 * t69;
t47 = qJ(5) + qJ(6);
t45 = sin(t47);
t46 = cos(t47);
t67 = t52 * t71;
t17 = -t37 * t67 + t38 * t83 - t52 * t69;
t68 = t49 * t71;
t30 = t37 * t48 - t57 * t68;
t51 = sin(qJ(4));
t55 = cos(qJ(4));
t8 = t17 * t55 + t30 * t51;
t93 = -t16 * t46 + t8 * t45;
t92 = t16 * t45 + t8 * t46;
t50 = sin(qJ(5));
t54 = cos(qJ(5));
t91 = -t16 * t54 + t8 * t50;
t90 = t16 * t50 + t8 * t54;
t85 = t17 * t51 - t30 * t55;
t63 = t72 * t82;
t59 = t57 * t53 + t56 * t63;
t84 = t59 * t71 - t82 * t79;
t81 = t45 * t55;
t80 = t46 * t55;
t78 = t48 * t51;
t77 = t48 * t55;
t76 = t49 * t53;
t75 = t49 * t56;
t74 = t50 * t55;
t73 = t54 * t55;
t70 = t48 * t76;
t65 = t72 * t48;
t39 = -t53 * t63 + t57 * t56;
t21 = t39 * t83 - t52 * t84;
t31 = t59 * t48 + t82 * t68;
t10 = -t21 * t51 + t31 * t55;
t28 = t52 * t65 + (t83 * t53 + t56 * t67) * t49;
t36 = -t48 * t75 + t72 * t71;
t61 = g(1) * t10 - g(2) * t85 + g(3) * (-t28 * t51 + t36 * t55);
t20 = t39 * t52 + t83 * t84;
t27 = t52 * t76 - t62 * t75 - t83 * t65;
t60 = g(1) * t20 + g(2) * t16 + g(3) * t27;
t35 = (-t53 * t67 + t83 * t56) * t49;
t34 = (t52 * t56 + t53 * t62) * t49;
t26 = t35 * t55 + t51 * t70;
t25 = -t39 * t67 - t59 * t83;
t24 = t39 * t62 - t59 * t52;
t23 = -t37 * t83 - t38 * t67;
t22 = -t37 * t52 + t38 * t62;
t15 = t28 * t55 + t36 * t51;
t13 = t25 * t55 + t39 * t78;
t12 = t23 * t55 + t38 * t78;
t11 = t21 * t55 + t31 * t51;
t6 = t11 * t54 + t20 * t50;
t5 = -t11 * t50 + t20 * t54;
t4 = t11 * t46 + t20 * t45;
t3 = -t11 * t45 + t20 * t46;
t2 = g(1) * t4 + g(2) * t92 - g(3) * (-t15 * t46 - t27 * t45);
t1 = -g(1) * t3 + g(2) * t93 - g(3) * (-t15 * t45 + t27 * t46);
t7 = [0, g(1) * t82 - g(2) * t57, g(1) * t57 + g(2) * t82, 0, 0, 0, 0, 0, g(1) * t38 - g(2) * t39, -g(1) * t37 + g(2) * t59, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t21, -g(1) * t16 + g(2) * t20, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t85 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t90 - g(2) * t6, -g(1) * t91 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t92 - g(2) * t4, -g(1) * t93 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t59 + g(2) * t37 - g(3) * t75, g(1) * t39 + g(2) * t38 + g(3) * t76, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t23 - g(3) * t35, g(1) * t24 + g(2) * t22 + g(3) * t34, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t12 - g(3) * t26, -g(1) * (-t25 * t51 + t39 * t77) - g(2) * (-t23 * t51 + t38 * t77) - g(3) * (-t35 * t51 + t55 * t70) 0, 0, 0, 0, 0, -g(1) * (t13 * t54 + t24 * t50) - g(2) * (t12 * t54 + t22 * t50) - g(3) * (t26 * t54 + t34 * t50) -g(1) * (-t13 * t50 + t24 * t54) - g(2) * (-t12 * t50 + t22 * t54) - g(3) * (-t26 * t50 + t34 * t54) 0, 0, 0, 0, 0, -g(1) * (t13 * t46 + t24 * t45) - g(2) * (t12 * t46 + t22 * t45) - g(3) * (t26 * t46 + t34 * t45) -g(1) * (-t13 * t45 + t24 * t46) - g(2) * (-t12 * t45 + t22 * t46) - g(3) * (-t26 * t45 + t34 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, g(1) * t21 + g(2) * t17 + g(3) * t28, 0, 0, 0, 0, 0, t60 * t55, -t60 * t51, 0, 0, 0, 0, 0, -g(1) * (-t20 * t73 + t21 * t50) - g(2) * (-t16 * t73 + t17 * t50) - g(3) * (-t27 * t73 + t28 * t50) -g(1) * (t20 * t74 + t21 * t54) - g(2) * (t16 * t74 + t17 * t54) - g(3) * (t27 * t74 + t28 * t54) 0, 0, 0, 0, 0, -g(1) * (-t20 * t80 + t21 * t45) - g(2) * (-t16 * t80 + t17 * t45) - g(3) * (-t27 * t80 + t28 * t45) -g(1) * (t20 * t81 + t21 * t46) - g(2) * (t16 * t81 + t17 * t46) - g(3) * (t27 * t81 + t28 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, g(1) * t11 + g(2) * t8 + g(3) * t15, 0, 0, 0, 0, 0, -t61 * t54, t61 * t50, 0, 0, 0, 0, 0, -t61 * t46, t61 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t91 - g(3) * (-t15 * t50 + t27 * t54) g(1) * t6 + g(2) * t90 - g(3) * (-t15 * t54 - t27 * t50) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t7;
