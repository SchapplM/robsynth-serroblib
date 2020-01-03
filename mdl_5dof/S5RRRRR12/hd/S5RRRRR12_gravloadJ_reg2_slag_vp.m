% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t71 = sin(pkin(5));
t78 = cos(qJ(1));
t126 = t71 * t78;
t70 = sin(pkin(6));
t118 = t70 * t126;
t132 = cos(qJ(3));
t121 = cos(pkin(5));
t110 = t78 * t121;
t131 = sin(qJ(1));
t133 = cos(qJ(2));
t75 = sin(qJ(2));
t54 = -t133 * t110 + t131 * t75;
t55 = t75 * t110 + t131 * t133;
t74 = sin(qJ(3));
t120 = cos(pkin(6));
t99 = t120 * t132;
t20 = t132 * t118 + t54 * t99 + t55 * t74;
t112 = t71 * t120;
t137 = t78 * t112 - t54 * t70;
t111 = t74 * t120;
t21 = -t54 * t111 - t74 * t118 + t55 * t132;
t73 = sin(qJ(4));
t77 = cos(qJ(4));
t5 = -t137 * t73 + t21 * t77;
t72 = sin(qJ(5));
t76 = cos(qJ(5));
t142 = -t20 * t76 + t5 * t72;
t141 = t20 * t72 + t5 * t76;
t100 = t121 * t131;
t85 = t133 * t100 + t78 * t75;
t138 = t131 * t112 + t85 * t70;
t136 = -t137 * t77 - t21 * t73;
t116 = t71 * t131;
t135 = t70 * t116 - t85 * t120;
t134 = pkin(9) * t70;
t129 = t70 * t73;
t128 = t70 * t77;
t127 = t71 * t75;
t125 = t72 * t77;
t124 = t76 * t77;
t117 = t71 * t133;
t119 = t70 * t127;
t123 = pkin(2) * t117 + pkin(9) * t119;
t122 = t78 * pkin(1) + pkin(8) * t116;
t115 = -t20 * pkin(3) + t21 * pkin(10);
t56 = -t75 * t100 + t78 * t133;
t24 = -t135 * t132 + t56 * t74;
t25 = t56 * t132 + t135 * t74;
t114 = -t24 * pkin(3) + t25 * pkin(10);
t109 = t121 * t70;
t37 = -t132 * t109 - t99 * t117 + t74 * t127;
t38 = t74 * t109 + (t133 * t111 + t132 * t75) * t71;
t113 = -t37 * pkin(3) + t38 * pkin(10);
t107 = -t54 * pkin(2) + t55 * t134;
t106 = -t85 * pkin(2) + t56 * t134;
t8 = -t138 * t77 + t25 * t73;
t105 = g(1) * t136 + g(2) * t8;
t103 = -t131 * pkin(1) + pkin(8) * t126;
t102 = -pkin(4) * t77 - pkin(11) * t73;
t101 = -g(1) * t20 + g(2) * t24;
t47 = (t133 * t74 + t75 * t99) * t71;
t48 = (-t75 * t111 + t132 * t133) * t71;
t98 = t48 * pkin(3) + t47 * pkin(10) + t123;
t96 = -g(1) * t78 - g(2) * t131;
t53 = -t70 * t117 + t121 * t120;
t18 = -t38 * t73 + t53 * t77;
t95 = g(1) * t8 - g(2) * t136 - g(3) * t18;
t19 = t38 * t77 + t53 * t73;
t9 = t138 * t73 + t25 * t77;
t94 = g(1) * t9 + g(2) * t5 + g(3) * t19;
t29 = -t55 * t111 - t54 * t132;
t10 = -t55 * t128 + t29 * t73;
t31 = -t56 * t111 - t85 * t132;
t12 = -t56 * t128 + t31 * t73;
t32 = -t77 * t119 + t48 * t73;
t93 = g(1) * t12 + g(2) * t10 + g(3) * t32;
t92 = g(1) * t24 + g(2) * t20 + g(3) * t37;
t91 = g(1) * t25 + g(2) * t21 + g(3) * t38;
t28 = -t54 * t74 + t55 * t99;
t30 = t56 * t99 - t85 * t74;
t90 = g(1) * t30 + g(2) * t28 + g(3) * t47;
t89 = t29 * pkin(3) + t28 * pkin(10) + t107;
t88 = t31 * pkin(3) + t30 * pkin(10) + t106;
t87 = g(1) * t56 + g(2) * t55 + g(3) * t127;
t86 = -t55 * pkin(2) + t137 * pkin(9) + t103;
t83 = -pkin(3) * t21 - pkin(10) * t20 + t86;
t81 = t56 * pkin(2) + t138 * pkin(9) + t122;
t79 = t25 * pkin(3) + t24 * pkin(10) + t81;
t33 = t73 * t119 + t48 * t77;
t13 = t56 * t129 + t31 * t77;
t11 = t55 * t129 + t29 * t77;
t3 = t24 * t72 + t9 * t76;
t2 = t24 * t76 - t9 * t72;
t1 = t92 * t73;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t131 - g(2) * t78, -t96, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t55 - g(2) * t56, -g(1) * t54 + g(2) * t85, t96 * t71, -g(1) * t103 - g(2) * t122, 0, 0, 0, 0, 0, 0, g(1) * t21 - g(2) * t25, t101, -g(1) * t137 - g(2) * t138, -g(1) * t86 - g(2) * t81, 0, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t105, -t101, -g(1) * t83 - g(2) * t79, 0, 0, 0, 0, 0, 0, g(1) * t141 - g(2) * t3, -g(1) * t142 - g(2) * t2, -t105, -g(1) * (-pkin(4) * t5 + pkin(11) * t136 + t83) - g(2) * (t9 * pkin(4) + t8 * pkin(11) + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t85 + g(2) * t54 - g(3) * t117, t87, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t29 - g(3) * t48, t90, -t87 * t70, -g(1) * t106 - g(2) * t107 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t33, t93, -t90, -g(1) * t88 - g(2) * t89 - g(3) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t76 + t30 * t72) - g(2) * (t11 * t76 + t28 * t72) - g(3) * (t33 * t76 + t47 * t72), -g(1) * (-t13 * t72 + t30 * t76) - g(2) * (-t11 * t72 + t28 * t76) - g(3) * (-t33 * t72 + t47 * t76), -t93, -g(1) * (t13 * pkin(4) + t12 * pkin(11) + t88) - g(2) * (t11 * pkin(4) + t10 * pkin(11) + t89) - g(3) * (t33 * pkin(4) + t32 * pkin(11) + t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t91, 0, 0, 0, 0, 0, 0, 0, 0, t92 * t77, -t1, -t91, -g(1) * t114 - g(2) * t115 - g(3) * t113, 0, 0, 0, 0, 0, 0, -g(1) * (-t24 * t124 + t25 * t72) - g(2) * (-t20 * t124 + t21 * t72) - g(3) * (-t37 * t124 + t38 * t72), -g(1) * (t24 * t125 + t25 * t76) - g(2) * (t20 * t125 + t21 * t76) - g(3) * (t37 * t125 + t38 * t76), t1, -g(1) * (t102 * t24 + t114) - g(2) * (t102 * t20 + t115) - g(3) * (t102 * t37 + t113); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t94, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t76, -t95 * t72, -t94, -g(1) * (-t8 * pkin(4) + t9 * pkin(11)) - g(2) * (pkin(4) * t136 + t5 * pkin(11)) - g(3) * (t18 * pkin(4) + t19 * pkin(11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t142 - g(3) * (-t19 * t72 + t37 * t76), g(1) * t3 + g(2) * t141 - g(3) * (-t19 * t76 - t37 * t72), 0, 0;];
taug_reg = t4;
