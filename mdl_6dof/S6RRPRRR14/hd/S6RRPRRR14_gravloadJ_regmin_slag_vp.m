% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-09 11:42:45
% EndTime: 2019-05-09 11:42:58
% DurationCPUTime: 1.69s
% Computational Cost: add. (1268->150), mult. (3679->289), div. (0->0), fcn. (4867->18), ass. (0->99)
t116 = cos(qJ(4));
t68 = sin(qJ(2));
t69 = sin(qJ(1));
t72 = cos(qJ(2));
t73 = cos(qJ(1));
t98 = cos(pkin(6));
t93 = t73 * t98;
t53 = t68 * t93 + t69 * t72;
t58 = sin(pkin(14));
t62 = cos(pkin(14));
t61 = sin(pkin(6));
t105 = t61 * t73;
t52 = t69 * t68 - t72 * t93;
t60 = sin(pkin(7));
t64 = cos(pkin(7));
t90 = t60 * t105 + t52 * t64;
t118 = t53 * t58 + t90 * t62;
t48 = t64 * t105 - t52 * t60;
t59 = sin(pkin(8));
t63 = cos(pkin(8));
t120 = t118 * t63 + t48 * t59;
t37 = -t53 * t62 + t90 * t58;
t67 = sin(qJ(4));
t15 = t37 * t116 + t120 * t67;
t25 = t118 * t59 - t48 * t63;
t66 = sin(qJ(5));
t71 = cos(qJ(5));
t5 = t15 * t71 - t25 * t66;
t65 = sin(qJ(6));
t129 = t5 * t65;
t70 = cos(qJ(6));
t128 = t5 * t70;
t127 = t15 * t66 + t25 * t71;
t124 = t37 * t67;
t102 = t64 * t72;
t92 = t98 * t60;
t75 = t62 * t92 + (t62 * t102 - t58 * t68) * t61;
t106 = t61 * t72;
t82 = t60 * t106 - t98 * t64;
t122 = -t82 * t59 + t75 * t63;
t94 = t69 * t98;
t55 = -t68 * t94 + t73 * t72;
t107 = t61 * t69;
t54 = -t73 * t68 - t72 * t94;
t88 = t60 * t107 + t54 * t64;
t78 = t55 * t58 - t88 * t62;
t89 = t64 * t107 - t54 * t60;
t121 = -t89 * t59 + t78 * t63;
t119 = -g(1) * t48 - g(2) * t89;
t111 = t58 * t64;
t110 = t60 * t59;
t109 = t60 * t63;
t108 = t61 * t68;
t104 = t62 * t64;
t103 = t64 * t68;
t101 = t65 * t71;
t100 = t70 * t71;
t99 = qJ(3) * t60;
t97 = t60 * t108;
t96 = t59 * t116;
t95 = t63 * t116;
t91 = t60 * t96;
t47 = t62 * t108 + (t61 * t102 + t92) * t58;
t24 = t47 * t116 + t122 * t67;
t34 = -t75 * t59 - t82 * t63;
t38 = t55 * t62 + t88 * t58;
t17 = t38 * t116 - t121 * t67;
t27 = t78 * t59 + t89 * t63;
t6 = -t17 * t66 + t27 * t71;
t87 = g(1) * t6 + g(2) * t127 + g(3) * (-t24 * t66 + t34 * t71);
t12 = t120 * t116 - t124;
t16 = t121 * t116 + t38 * t67;
t23 = -t122 * t116 + t47 * t67;
t86 = g(1) * t16 + g(2) * t12 + g(3) * t23;
t81 = g(1) * t55 + g(2) * t53 + g(3) * t108;
t51 = (-t58 * t103 + t62 * t72) * t61;
t50 = (-t62 * t103 - t58 * t72) * t61;
t43 = -t50 * t59 + t63 * t97;
t42 = -t55 * t111 + t54 * t62;
t41 = -t55 * t104 - t54 * t58;
t40 = -t53 * t111 - t52 * t62;
t39 = -t53 * t104 + t52 * t58;
t31 = t55 * t109 - t41 * t59;
t30 = t53 * t109 - t39 * t59;
t29 = t51 * t116 + (t50 * t63 + t59 * t97) * t67;
t28 = -t91 * t108 - t50 * t95 + t51 * t67;
t22 = t29 * t71 + t43 * t66;
t21 = t42 * t116 + (t55 * t110 + t41 * t63) * t67;
t20 = -t41 * t95 + t42 * t67 - t55 * t91;
t19 = t40 * t116 + (t53 * t110 + t39 * t63) * t67;
t18 = -t39 * t95 + t40 * t67 - t53 * t91;
t14 = -t118 * t95 - t48 * t96 + t124;
t11 = t24 * t71 + t34 * t66;
t9 = t21 * t71 + t31 * t66;
t8 = t19 * t71 + t30 * t66;
t7 = t17 * t71 + t27 * t66;
t2 = t16 * t65 + t7 * t70;
t1 = t16 * t70 - t7 * t65;
t3 = [0, g(1) * t69 - g(2) * t73, g(1) * t73 + g(2) * t69, 0, 0, 0, 0, 0, g(1) * t53 - g(2) * t55, -g(1) * t52 - g(2) * t54, -g(1) * t37 - g(2) * t38, -g(1) * t118 + g(2) * t78, t119, -g(1) * (-t69 * pkin(1) - t53 * pkin(2) + pkin(10) * t105) - g(2) * (t73 * pkin(1) + t55 * pkin(2) + pkin(10) * t107) + t119 * qJ(3), 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, g(1) * t14 + g(2) * t16, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, g(1) * t127 - g(2) * t6, 0, 0, 0, 0, 0, -g(1) * (t14 * t65 + t128) - g(2) * t2, -g(1) * (t14 * t70 - t129) - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t54 + g(2) * t52 - g(3) * t106, t81, -g(1) * t42 - g(2) * t40 - g(3) * t51, -g(1) * t41 - g(2) * t39 - g(3) * t50, -t81 * t60, -g(1) * (t54 * pkin(2) + t55 * t99) - g(2) * (-t52 * pkin(2) + t53 * t99) - g(3) * (pkin(2) * t72 + t68 * t99) * t61, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t19 - g(3) * t29, g(1) * t20 + g(2) * t18 + g(3) * t28, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t22, -g(1) * (-t21 * t66 + t31 * t71) - g(2) * (-t19 * t66 + t30 * t71) - g(3) * (-t29 * t66 + t43 * t71) 0, 0, 0, 0, 0, -g(1) * (t20 * t65 + t9 * t70) - g(2) * (t18 * t65 + t8 * t70) - g(3) * (t22 * t70 + t28 * t65) -g(1) * (t20 * t70 - t9 * t65) - g(2) * (t18 * t70 - t8 * t65) - g(3) * (-t22 * t65 + t28 * t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t89 + g(2) * t48 + g(3) * t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, g(1) * t17 - g(2) * t15 + g(3) * t24, 0, 0, 0, 0, 0, t86 * t71, -t86 * t66, 0, 0, 0, 0, 0, -g(1) * (-t16 * t100 + t17 * t65) - g(2) * (-t12 * t100 - t15 * t65) - g(3) * (-t23 * t100 + t24 * t65) -g(1) * (t16 * t101 + t17 * t70) - g(2) * (t12 * t101 - t15 * t70) - g(3) * (t23 * t101 + t24 * t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, g(1) * t7 - g(2) * t5 + g(3) * t11, 0, 0, 0, 0, 0, -t87 * t70, t87 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * (t12 * t70 + t129) - g(3) * (-t11 * t65 + t23 * t70) g(1) * t2 - g(2) * (-t12 * t65 + t128) - g(3) * (-t11 * t70 - t23 * t65);];
taug_reg  = t3;
