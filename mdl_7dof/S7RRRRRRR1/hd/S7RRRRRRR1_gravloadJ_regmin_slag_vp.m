% Calculate minimal parameter regressor of gravitation load for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% taug_reg [7x45]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S7RRRRRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_gravloadJ_regmin_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_gravloadJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 23:18:08
% EndTime: 2019-05-08 23:18:16
% DurationCPUTime: 2.08s
% Computational Cost: add. (684->165), mult. (1942->319), div. (0->0), fcn. (2528->14), ass. (0->98)
t75 = cos(qJ(3));
t68 = sin(qJ(3));
t77 = cos(qJ(1));
t89 = t77 * t68;
t70 = sin(qJ(1));
t76 = cos(qJ(2));
t95 = t70 * t76;
t54 = t75 * t95 + t89;
t74 = cos(qJ(4));
t67 = sin(qJ(4));
t69 = sin(qJ(2));
t99 = t69 * t67;
t36 = t54 * t74 + t70 * t99;
t88 = t77 * t75;
t53 = t68 * t95 - t88;
t66 = sin(qJ(5));
t73 = cos(qJ(5));
t16 = t36 * t66 + t53 * t73;
t17 = t36 * t73 - t53 * t66;
t98 = t69 * t74;
t35 = t54 * t67 - t70 * t98;
t65 = sin(qJ(6));
t72 = cos(qJ(6));
t4 = t17 * t72 + t35 * t65;
t64 = sin(qJ(7));
t71 = cos(qJ(7));
t118 = t16 * t71 + t4 * t64;
t117 = -t16 * t64 + t4 * t71;
t114 = t17 * t65 - t35 * t72;
t84 = g(1) * t77 + g(2) * t70;
t78 = -g(3) * t76 + t84 * t69;
t108 = t64 * t66;
t107 = t64 * t72;
t106 = t65 * t67;
t105 = t65 * t73;
t104 = t66 * t71;
t103 = t66 * t74;
t102 = t67 * t72;
t101 = t68 * t69;
t100 = t68 * t76;
t97 = t69 * t75;
t96 = t69 * t77;
t94 = t71 * t72;
t93 = t72 * t73;
t92 = t73 * t74;
t91 = t76 * t67;
t90 = t76 * t74;
t87 = t68 * t99;
t86 = t66 * t101;
t85 = t73 * t101;
t83 = g(1) * t70 - g(2) * t77;
t52 = t74 * t97 - t91;
t51 = t67 * t97 + t90;
t34 = t52 * t73 - t86;
t58 = -t70 * t68 + t76 * t88;
t42 = t58 * t74 + t67 * t96;
t57 = -t70 * t75 - t76 * t89;
t22 = t42 * t73 + t57 * t66;
t41 = t58 * t67 - t74 * t96;
t6 = -t22 * t65 + t41 * t72;
t82 = g(1) * t6 - g(2) * t114 + g(3) * (-t34 * t65 + t51 * t72);
t21 = t42 * t66 - t57 * t73;
t33 = t52 * t66 + t85;
t81 = g(1) * t21 + g(2) * t16 + g(3) * t33;
t80 = g(1) * t41 + g(2) * t35 + g(3) * t51;
t79 = -g(1) * t57 + g(2) * t53 + g(3) * t101;
t56 = t75 * t90 + t99;
t55 = t75 * t91 - t98;
t48 = t52 * t77;
t47 = t51 * t77;
t46 = t52 * t70;
t45 = t51 * t70;
t44 = (-t66 * t75 - t68 * t92) * t69;
t43 = -t73 * t97 + t74 * t86;
t40 = -t66 * t100 + t56 * t73;
t39 = t73 * t100 + t56 * t66;
t32 = -t48 * t73 + t77 * t86;
t31 = -t48 * t66 - t77 * t85;
t30 = -t46 * t73 + t70 * t86;
t29 = -t46 * t66 - t70 * t85;
t28 = t44 * t72 - t65 * t87;
t27 = t57 * t92 - t58 * t66;
t26 = t57 * t103 + t58 * t73;
t25 = -t53 * t92 - t54 * t66;
t24 = -t53 * t103 + t54 * t73;
t23 = -t51 * t93 + t52 * t65;
t20 = t40 * t72 + t55 * t65;
t15 = t34 * t72 + t51 * t65;
t13 = t32 * t72 - t47 * t65;
t12 = t30 * t72 - t45 * t65;
t11 = t57 * t106 + t27 * t72;
t10 = -t53 * t106 + t25 * t72;
t9 = -t41 * t93 + t42 * t65;
t8 = -t35 * t93 + t36 * t65;
t7 = t22 * t72 + t41 * t65;
t2 = -t21 * t64 + t7 * t71;
t1 = -t21 * t71 - t7 * t64;
t3 = [0, t83, t84, 0, 0, 0, 0, 0, t83 * t76, -t83 * t69, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t58, -g(1) * t53 - g(2) * t57, 0, 0, 0, 0, 0, g(1) * t36 - g(2) * t42, -g(1) * t35 + g(2) * t41, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t22, -g(1) * t16 + g(2) * t21, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t114 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t117 - g(2) * t2, -g(1) * t118 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t78, g(3) * t69 + t84 * t76, 0, 0, 0, 0, 0, t78 * t75, -t78 * t68, 0, 0, 0, 0, 0, g(1) * t48 + g(2) * t46 - g(3) * t56, -g(1) * t47 - g(2) * t45 + g(3) * t55, 0, 0, 0, 0, 0, -g(1) * t32 - g(2) * t30 - g(3) * t40, g(1) * t31 + g(2) * t29 + g(3) * t39, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t12 - g(3) * t20, -g(1) * (-t32 * t65 - t47 * t72) - g(2) * (-t30 * t65 - t45 * t72) - g(3) * (-t40 * t65 + t55 * t72) 0, 0, 0, 0, 0, -g(1) * (t13 * t71 - t31 * t64) - g(2) * (t12 * t71 - t29 * t64) - g(3) * (t20 * t71 - t39 * t64) -g(1) * (-t13 * t64 - t31 * t71) - g(2) * (-t12 * t64 - t29 * t71) - g(3) * (-t20 * t64 - t39 * t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, g(1) * t58 + g(2) * t54 + g(3) * t97, 0, 0, 0, 0, 0, t79 * t74, -t79 * t67, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t25 - g(3) * t44, g(1) * t26 + g(2) * t24 - g(3) * t43, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t10 - g(3) * t28, -g(1) * (t102 * t57 - t27 * t65) - g(2) * (-t53 * t102 - t25 * t65) - g(3) * (-t44 * t65 - t72 * t87) 0, 0, 0, 0, 0, -g(1) * (t11 * t71 - t26 * t64) - g(2) * (t10 * t71 - t24 * t64) - g(3) * (t28 * t71 + t43 * t64) -g(1) * (-t11 * t64 - t26 * t71) - g(2) * (-t10 * t64 - t24 * t71) - g(3) * (-t28 * t64 + t43 * t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, g(1) * t42 + g(2) * t36 + g(3) * t52, 0, 0, 0, 0, 0, t80 * t73, -t80 * t66, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t23, -g(1) * (t105 * t41 + t42 * t72) - g(2) * (t105 * t35 + t36 * t72) - g(3) * (t105 * t51 + t52 * t72) 0, 0, 0, 0, 0, -g(1) * (t108 * t41 + t9 * t71) - g(2) * (t108 * t35 + t8 * t71) - g(3) * (t108 * t51 + t23 * t71) -g(1) * (t104 * t41 - t9 * t64) - g(2) * (t104 * t35 - t8 * t64) - g(3) * (t104 * t51 - t23 * t64); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, g(1) * t22 + g(2) * t17 + g(3) * t34, 0, 0, 0, 0, 0, t81 * t72, -t81 * t65, 0, 0, 0, 0, 0, -g(1) * (-t21 * t94 - t22 * t64) - g(2) * (-t16 * t94 - t17 * t64) - g(3) * (-t33 * t94 - t34 * t64) -g(1) * (t107 * t21 - t22 * t71) - g(2) * (t107 * t16 - t17 * t71) - g(3) * (t107 * t33 - t34 * t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, g(1) * t7 + g(2) * t4 + g(3) * t15, 0, 0, 0, 0, 0, -t82 * t71, t82 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t118 - g(3) * (-t15 * t64 - t33 * t71) g(1) * t2 + g(2) * t117 - g(3) * (-t15 * t71 + t33 * t64);];
taug_reg  = t3;
