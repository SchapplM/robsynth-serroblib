% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t111 = cos(qJ(4));
t61 = sin(pkin(8));
t62 = sin(pkin(7));
t110 = t61 * t62;
t67 = sin(qJ(5));
t109 = t61 * t67;
t72 = cos(qJ(5));
t108 = t61 * t72;
t64 = cos(pkin(8));
t107 = t62 * t64;
t63 = sin(pkin(6));
t70 = sin(qJ(2));
t106 = t63 * t70;
t74 = cos(qJ(2));
t105 = t63 * t74;
t68 = sin(qJ(4));
t104 = t64 * t68;
t65 = cos(pkin(7));
t69 = sin(qJ(3));
t103 = t65 * t69;
t73 = cos(qJ(3));
t102 = t65 * t73;
t66 = sin(qJ(6));
t101 = t66 * t72;
t100 = t69 * t70;
t99 = t69 * t74;
t98 = t70 * t73;
t71 = cos(qJ(6));
t97 = t71 * t72;
t96 = t73 * t74;
t95 = cos(pkin(6));
t94 = cos(pkin(14));
t93 = sin(pkin(14));
t92 = t62 * t106;
t91 = t64 * t111;
t90 = t63 * t94;
t89 = t63 * t93;
t88 = t95 * t62;
t87 = t111 * t110;
t86 = t62 * t90;
t85 = t95 * t94;
t84 = t95 * t93;
t55 = -t93 * t70 + t74 * t85;
t56 = t70 * t85 + t93 * t74;
t38 = -t56 * t69 + (t55 * t65 - t86) * t73;
t39 = t55 * t103 + t56 * t73 - t69 * t86;
t81 = -t55 * t62 - t65 * t90;
t77 = t81 * t61;
t12 = t39 * t111 + (t38 * t64 + t77) * t68;
t58 = -t70 * t84 + t94 * t74;
t57 = -t94 * t70 - t74 * t84;
t78 = t57 * t65 + t62 * t89;
t40 = -t58 * t69 + t78 * t73;
t41 = t58 * t73 + t78 * t69;
t80 = -t57 * t62 + t65 * t89;
t76 = t80 * t61;
t14 = t41 * t111 + (t40 * t64 + t76) * t68;
t50 = t73 * t88 + (t65 * t96 - t100) * t63;
t51 = t69 * t88 + (t65 * t99 + t98) * t63;
t79 = -t62 * t105 + t95 * t65;
t75 = t79 * t61;
t26 = t51 * t111 + (t50 * t64 + t75) * t68;
t27 = -t38 * t61 + t81 * t64;
t28 = -t40 * t61 + t80 * t64;
t37 = -t50 * t61 + t79 * t64;
t83 = g(1) * (-t14 * t67 + t28 * t72) + g(2) * (-t12 * t67 + t27 * t72) + g(3) * (-t26 * t67 + t37 * t72);
t11 = -t111 * t77 - t38 * t91 + t39 * t68;
t13 = -t111 * t76 - t40 * t91 + t41 * t68;
t25 = -t111 * t75 - t50 * t91 + t51 * t68;
t82 = g(1) * t13 + g(2) * t11 + g(3) * t25;
t54 = (-t65 * t100 + t96) * t63;
t53 = (-t65 * t98 - t99) * t63;
t46 = -t53 * t61 + t64 * t92;
t45 = -t58 * t103 + t57 * t73;
t44 = -t58 * t102 - t57 * t69;
t43 = -t56 * t103 + t55 * t73;
t42 = -t56 * t102 - t55 * t69;
t34 = t54 * t111 + (t53 * t64 + t61 * t92) * t68;
t33 = -t87 * t106 - t53 * t91 + t54 * t68;
t32 = t58 * t107 - t44 * t61;
t31 = t56 * t107 - t42 * t61;
t30 = -t51 * t104 + t50 * t111;
t29 = t50 * t68 + t51 * t91;
t24 = t51 * t109 + t30 * t72;
t23 = t34 * t72 + t46 * t67;
t22 = -t41 * t104 + t40 * t111;
t21 = t40 * t68 + t41 * t91;
t20 = -t39 * t104 + t38 * t111;
t19 = t38 * t68 + t39 * t91;
t18 = t45 * t111 + (t58 * t110 + t44 * t64) * t68;
t17 = -t44 * t91 + t45 * t68 - t58 * t87;
t16 = t43 * t111 + (t56 * t110 + t42 * t64) * t68;
t15 = -t42 * t91 + t43 * t68 - t56 * t87;
t10 = t26 * t72 + t37 * t67;
t8 = t41 * t109 + t22 * t72;
t7 = t39 * t109 + t20 * t72;
t6 = t18 * t72 + t32 * t67;
t5 = t16 * t72 + t31 * t67;
t4 = t14 * t72 + t28 * t67;
t2 = t12 * t72 + t27 * t67;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(1) * t57 - g(2) * t55 - g(3) * t105, g(1) * t58 + g(2) * t56 + g(3) * t106, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t43 - g(3) * t54, -g(1) * t44 - g(2) * t42 - g(3) * t53, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16 - g(3) * t34, g(1) * t17 + g(2) * t15 + g(3) * t33, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t5 - g(3) * t23, -g(1) * (-t18 * t67 + t32 * t72) - g(2) * (-t16 * t67 + t31 * t72) - g(3) * (-t34 * t67 + t46 * t72) 0, 0, 0, 0, 0, -g(1) * (t17 * t66 + t6 * t71) - g(2) * (t15 * t66 + t5 * t71) - g(3) * (t23 * t71 + t33 * t66) -g(1) * (t17 * t71 - t6 * t66) - g(2) * (t15 * t71 - t5 * t66) - g(3) * (-t23 * t66 + t33 * t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t40 - g(2) * t38 - g(3) * t50, g(1) * t41 + g(2) * t39 + g(3) * t51, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t20 - g(3) * t30, g(1) * t21 + g(2) * t19 + g(3) * t29, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t7 - g(3) * t24, -g(1) * (t41 * t108 - t22 * t67) - g(2) * (t39 * t108 - t20 * t67) - g(3) * (t51 * t108 - t30 * t67) 0, 0, 0, 0, 0, -g(1) * (t21 * t66 + t8 * t71) - g(2) * (t19 * t66 + t7 * t71) - g(3) * (t24 * t71 + t29 * t66) -g(1) * (t21 * t71 - t8 * t66) - g(2) * (t19 * t71 - t7 * t66) - g(3) * (-t24 * t66 + t29 * t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, g(1) * t14 + g(2) * t12 + g(3) * t26, 0, 0, 0, 0, 0, t82 * t72, -t82 * t67, 0, 0, 0, 0, 0, -g(1) * (-t13 * t97 + t14 * t66) - g(2) * (-t11 * t97 + t12 * t66) - g(3) * (-t25 * t97 + t26 * t66) -g(1) * (t13 * t101 + t14 * t71) - g(2) * (t11 * t101 + t12 * t71) - g(3) * (t25 * t101 + t26 * t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, g(1) * t4 + g(2) * t2 + g(3) * t10, 0, 0, 0, 0, 0, -t83 * t71, t83 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t71 - t4 * t66) - g(2) * (t11 * t71 - t2 * t66) - g(3) * (-t10 * t66 + t25 * t71) -g(1) * (-t13 * t66 - t4 * t71) - g(2) * (-t11 * t66 - t2 * t71) - g(3) * (-t10 * t71 - t25 * t66);];
taug_reg  = t1;
