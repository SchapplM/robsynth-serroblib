% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t124 = pkin(5) * t72 + qJ(6) * t69;
t103 = cos(pkin(6));
t68 = sin(pkin(6));
t71 = sin(qJ(2));
t110 = t68 * t71;
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t123 = t103 * t73 - t70 * t110;
t109 = t68 * t73;
t102 = cos(pkin(11));
t74 = cos(qJ(2));
t67 = sin(pkin(11));
t95 = t67 * t103;
t52 = t102 * t74 - t71 * t95;
t122 = t67 * t109 - t52 * t70;
t66 = qJ(3) + qJ(4);
t64 = sin(t66);
t65 = cos(t66);
t121 = pkin(4) * t65 + pkin(10) * t64;
t87 = t103 * t102;
t50 = t67 * t74 + t71 * t87;
t94 = t68 * t102;
t24 = -t50 * t64 - t65 * t94;
t120 = t124 * t24;
t111 = t67 * t68;
t26 = t65 * t111 - t52 * t64;
t119 = t124 * t26;
t41 = t103 * t65 - t64 * t110;
t118 = t124 * t41;
t113 = t65 * t69;
t112 = t65 * t72;
t108 = t68 * t74;
t107 = t72 * t74;
t49 = t67 * t71 - t74 * t87;
t63 = t73 * pkin(3) + pkin(2);
t75 = -pkin(9) - pkin(8);
t106 = -t49 * t63 - t50 * t75;
t51 = t102 * t71 + t74 * t95;
t105 = -t51 * t63 - t52 * t75;
t99 = t69 * t108;
t25 = t50 * t65 - t64 * t94;
t98 = t24 * pkin(4) + t25 * pkin(10);
t27 = t64 * t111 + t52 * t65;
t97 = t26 * pkin(4) + t27 * pkin(10);
t42 = t103 * t64 + t65 * t110;
t96 = t41 * pkin(4) + t42 * pkin(10);
t92 = -t121 * t49 + t106;
t91 = -t121 * t51 + t105;
t90 = t122 * pkin(3);
t89 = t63 * t108 - t75 * t110;
t88 = t123 * pkin(3);
t86 = t121 * t108 + t89;
t11 = t27 * t69 - t51 * t72;
t28 = t68 * t107 + t42 * t69;
t9 = t25 * t69 - t49 * t72;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t28;
t10 = t25 * t72 + t49 * t69;
t12 = t27 * t72 + t51 * t69;
t29 = t42 * t72 - t99;
t85 = g(1) * t12 + g(2) * t10 + g(3) * t29;
t13 = -t49 * t113 - t50 * t72;
t15 = -t51 * t113 - t52 * t72;
t32 = -t72 * t110 + t65 * t99;
t84 = g(1) * t15 + g(2) * t13 + g(3) * t32;
t83 = g(1) * t26 + g(2) * t24 + g(3) * t41;
t7 = g(1) * t27 + g(2) * t25 + g(3) * t42;
t82 = t90 + t97;
t81 = -t50 * t70 - t73 * t94;
t80 = -g(1) * t51 - g(2) * t49 + g(3) * t108;
t79 = g(1) * t52 + g(2) * t50 + g(3) * t110;
t78 = t88 + t96;
t77 = t81 * pkin(3);
t76 = t77 + t98;
t33 = (t65 * t107 + t69 * t71) * t68;
t16 = -t51 * t112 + t52 * t69;
t14 = -t49 * t112 + t50 * t69;
t8 = t80 * t64;
t4 = t83 * t72;
t3 = t83 * t69;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t33;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t79, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t73, t80 * t70, -t79, -g(1) * (-t51 * pkin(2) + t52 * pkin(8)) - g(2) * (-t49 * pkin(2) + t50 * pkin(8)) - g(3) * (pkin(2) * t74 + pkin(8) * t71) * t68, 0, 0, 0, 0, 0, 0, -t80 * t65, t8, -t79, -g(1) * t105 - g(2) * t106 - g(3) * t89, 0, 0, 0, 0, 0, 0, t2, t84, -t8, -g(1) * t91 - g(2) * t92 - g(3) * t86, 0, 0, 0, 0, 0, 0, t2, -t8, -t84, -g(1) * (t16 * pkin(5) + t15 * qJ(6) + t91) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t92) - g(3) * (t33 * pkin(5) + t32 * qJ(6) + t86); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t122 - g(2) * t81 - g(3) * t123, -g(1) * (-t70 * t111 - t52 * t73) - g(2) * (-t50 * t73 + t70 * t94) - g(3) * (-t103 * t70 - t71 * t109) 0, 0, 0, 0, 0, 0, 0, 0, -t83, t7, 0, -g(1) * t90 - g(2) * t77 - g(3) * t88, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(1) * t82 - g(2) * t76 - g(3) * t78, 0, 0, 0, 0, 0, 0, -t4, -t7, -t3, -g(1) * (t82 + t119) - g(2) * (t76 + t120) - g(3) * (t78 + t118); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(1) * t97 - g(2) * t98 - g(3) * t96, 0, 0, 0, 0, 0, 0, -t4, -t7, -t3, -g(1) * (t97 + t119) - g(2) * (t98 + t120) - g(3) * (t96 + t118); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t85, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t85, -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t28 * pkin(5) + t29 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
