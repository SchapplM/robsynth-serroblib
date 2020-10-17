% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:06:17
% EndTime: 2019-05-07 19:06:20
% DurationCPUTime: 0.84s
% Computational Cost: add. (773->154), mult. (2015->214), div. (0->0), fcn. (2509->10), ass. (0->98)
t90 = sin(qJ(3));
t93 = cos(qJ(3));
t155 = -pkin(3) * t93 - pkin(10) * t90;
t147 = cos(qJ(1));
t88 = sin(pkin(6));
t127 = t88 * t147;
t130 = cos(pkin(6));
t112 = t130 * t147;
t146 = sin(qJ(1));
t91 = sin(qJ(2));
t94 = cos(qJ(2));
t71 = t91 * t112 + t146 * t94;
t44 = -t90 * t127 + t71 * t93;
t70 = -t94 * t112 + t146 * t91;
t89 = sin(qJ(4));
t92 = cos(qJ(4));
t19 = t44 * t89 - t70 * t92;
t20 = t44 * t92 + t70 * t89;
t120 = -t93 * t127 - t71 * t90;
t132 = qJ(5) * t89;
t145 = t120 * t92;
t154 = pkin(4) * t145 + t120 * t132;
t126 = t88 * t146;
t111 = t130 * t146;
t73 = -t91 * t111 + t147 * t94;
t47 = -t93 * t126 + t73 * t90;
t144 = t47 * t92;
t153 = -pkin(4) * t144 - t47 * t132;
t140 = t88 * t91;
t68 = t130 * t93 - t90 * t140;
t143 = t68 * t92;
t152 = pkin(4) * t143 + t68 * t132;
t149 = t120 * pkin(10);
t148 = t47 * pkin(10);
t139 = t88 * t94;
t138 = t89 * t93;
t137 = t92 * t93;
t136 = t93 * t94;
t135 = pkin(10) - qJ(6);
t134 = pkin(2) * t139 + pkin(9) * t140;
t133 = t147 * pkin(1) + pkin(8) * t126;
t131 = qJ(6) * t90;
t129 = t90 * t139;
t128 = t89 * t139;
t125 = -t70 * pkin(2) + t71 * pkin(9);
t72 = t94 * t111 + t147 * t91;
t124 = -t72 * pkin(2) + t73 * pkin(9);
t37 = t120 * pkin(3);
t123 = t44 * pkin(10) + t37;
t39 = t47 * pkin(3);
t48 = t90 * t126 + t73 * t93;
t122 = t48 * pkin(10) - t39;
t63 = t68 * pkin(3);
t69 = t130 * t90 + t93 * t140;
t121 = t69 * pkin(10) + t63;
t119 = -t19 * pkin(4) + t20 * qJ(5);
t23 = t48 * t89 - t72 * t92;
t24 = t48 * t92 + t72 * t89;
t118 = -t23 * pkin(4) + t24 * qJ(5);
t41 = t92 * t139 + t69 * t89;
t42 = t69 * t92 - t128;
t117 = -t41 * pkin(4) + t42 * qJ(5);
t116 = t88 * pkin(3) * t136 + pkin(10) * t129 + t134;
t115 = -t146 * pkin(1) + pkin(8) * t127;
t114 = -g(1) * t19 + g(2) * t23;
t14 = g(1) * t120 + g(2) * t47;
t113 = g(1) * t70 - g(2) * t72;
t110 = t155 * t70 + t125;
t109 = t155 * t72 + t124;
t108 = t73 * pkin(2) + t72 * pkin(9) + t133;
t107 = t48 * pkin(3) + t108;
t2 = g(1) * t23 + g(2) * t19 + g(3) * t41;
t106 = g(1) * t24 + g(2) * t20 + g(3) * t42;
t28 = -t70 * t138 - t71 * t92;
t30 = -t72 * t138 - t73 * t92;
t50 = t93 * t128 - t92 * t140;
t105 = g(1) * t30 + g(2) * t28 + g(3) * t50;
t10 = g(1) * t47 - g(2) * t120 - g(3) * t68;
t12 = g(1) * t48 + g(2) * t44 + g(3) * t69;
t51 = (t92 * t136 + t89 * t91) * t88;
t104 = t51 * pkin(4) + t50 * qJ(5) + t116;
t103 = g(1) * t147 + g(2) * t146;
t102 = -t71 * pkin(2) - t70 * pkin(9) + t115;
t101 = -g(1) * t72 - g(2) * t70 + g(3) * t139;
t100 = g(1) * t73 + g(2) * t71 + g(3) * t140;
t29 = -t70 * t137 + t71 * t89;
t99 = t29 * pkin(4) + t28 * qJ(5) + t110;
t31 = -t72 * t137 + t73 * t89;
t98 = t31 * pkin(4) + t30 * qJ(5) + t109;
t97 = -pkin(3) * t44 + t102;
t96 = t24 * pkin(4) + t23 * qJ(5) + t107;
t95 = -pkin(4) * t20 - qJ(5) * t19 + t97;
t25 = t101 * t90;
t9 = t10 * t92;
t8 = t10 * t89;
t7 = g(1) * t20 - g(2) * t24;
t5 = -g(1) * t31 - g(2) * t29 - g(3) * t51;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t146 - g(2) * t147, t103, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t73, -t113, -t103 * t88, -g(1) * t115 - g(2) * t133, 0, 0, 0, 0, 0, 0, g(1) * t44 - g(2) * t48, t14, t113, -g(1) * t102 - g(2) * t108, 0, 0, 0, 0, 0, 0, t7, t114, -t14, -g(1) * (t97 + t149) - g(2) * (t107 + t148) 0, 0, 0, 0, 0, 0, t7, -t14, -t114, -g(1) * (t95 + t149) - g(2) * (t96 + t148) 0, 0, 0, 0, 0, 0, t7, -t114, t14, -g(1) * (-pkin(5) * t20 + t120 * t135 + t95) - g(2) * (t24 * pkin(5) + t135 * t47 + t96); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, t100, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * t93, t25, -t100, -g(1) * t124 - g(2) * t125 - g(3) * t134, 0, 0, 0, 0, 0, 0, t5, t105, -t25, -g(1) * t109 - g(2) * t110 - g(3) * t116, 0, 0, 0, 0, 0, 0, t5, -t25, -t105, -g(1) * t98 - g(2) * t99 - g(3) * t104, 0, 0, 0, 0, 0, 0, t5, -t105, t25, -g(1) * (t31 * pkin(5) + t72 * t131 + t98) - g(2) * (t29 * pkin(5) + t70 * t131 + t99) - g(3) * (t51 * pkin(5) - qJ(6) * t129 + t104); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, -t12, -g(1) * t122 - g(2) * t123 - g(3) * t121, 0, 0, 0, 0, 0, 0, t9, -t12, t8, -g(1) * (t122 + t153) - g(2) * (t123 + t154) - g(3) * (t121 + t152) 0, 0, 0, 0, 0, 0, t9, t8, t12, -g(1) * (-pkin(5) * t144 + t135 * t48 + t153 - t39) - g(2) * (pkin(5) * t145 + t135 * t44 + t154 + t37) - g(3) * (pkin(5) * t143 + t135 * t69 + t152 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t106, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t106, -g(1) * t118 - g(2) * t119 - g(3) * t117, 0, 0, 0, 0, 0, 0, t2, -t106, 0, -g(1) * (-t23 * pkin(5) + t118) - g(2) * (-t19 * pkin(5) + t119) - g(3) * (-t41 * pkin(5) + t117); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10;];
taug_reg  = t1;
