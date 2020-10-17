% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:05:41
% EndTime: 2019-05-06 18:05:44
% DurationCPUTime: 0.89s
% Computational Cost: add. (936->141), mult. (2405->216), div. (0->0), fcn. (3097->12), ass. (0->95)
t87 = sin(qJ(4));
t91 = cos(qJ(4));
t143 = pkin(4) * t91 + pkin(10) * t87;
t84 = sin(pkin(6));
t93 = cos(qJ(1));
t134 = t84 * t93;
t126 = cos(pkin(11));
t83 = sin(pkin(11));
t88 = sin(qJ(2));
t92 = cos(qJ(2));
t108 = t92 * t126 - t88 * t83;
t73 = -t88 * t126 - t92 * t83;
t85 = cos(pkin(6));
t65 = t73 * t85;
t89 = sin(qJ(1));
t48 = t108 * t89 - t93 * t65;
t26 = -t87 * t134 + t48 * t91;
t64 = t108 * t85;
t47 = t93 * t64 + t89 * t73;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t9 = t26 * t86 + t47 * t90;
t10 = t26 * t90 - t47 * t86;
t135 = t84 * t92;
t128 = t93 * t88;
t130 = t89 * t92;
t69 = -t85 * t130 - t128;
t142 = -g(1) * t69 - g(3) * t135;
t136 = t84 * t89;
t133 = t86 * t91;
t131 = t89 * t88;
t129 = t90 * t91;
t127 = t93 * t92;
t124 = t85 * t127;
t120 = -t91 * t134 - t48 * t87;
t123 = pkin(4) * t120 + t26 * pkin(10);
t49 = -t108 * t93 - t89 * t65;
t29 = -t91 * t136 - t49 * t87;
t30 = t87 * t136 - t49 * t91;
t122 = -t29 * pkin(4) + t30 * pkin(10);
t63 = t73 * t84;
t53 = t63 * t87 + t85 * t91;
t54 = -t63 * t91 + t85 * t87;
t121 = t53 * pkin(4) + t54 * pkin(10);
t66 = t85 * t88 * pkin(2) + (-pkin(8) - qJ(3)) * t84;
t82 = t92 * pkin(2) + pkin(1);
t119 = -t89 * t66 + t93 * t82;
t62 = t108 * t84;
t117 = pkin(2) * t135 + t62 * pkin(3) - t63 * pkin(9);
t50 = -t89 * t64 + t93 * t73;
t13 = t30 * t86 + t50 * t90;
t116 = -g(1) * t9 + g(2) * t13;
t115 = g(1) * t120 + g(2) * t29;
t114 = g(1) * t47 - g(2) * t50;
t113 = g(1) * t93 + g(2) * t89;
t112 = g(1) * t89 - g(2) * t93;
t111 = pkin(5) * t90 + qJ(6) * t86;
t110 = -t93 * t66 - t89 * t82;
t109 = t143 * t62 + t117;
t107 = -pkin(3) * t49 - t50 * pkin(9) + t119;
t19 = t54 * t86 + t62 * t90;
t1 = g(1) * t13 + g(2) * t9 + g(3) * t19;
t14 = t30 * t90 - t50 * t86;
t20 = t54 * t90 - t62 * t86;
t106 = g(1) * t14 + g(2) * t10 + g(3) * t20;
t15 = t47 * t133 - t48 * t90;
t17 = t50 * t133 + t49 * t90;
t31 = t62 * t133 + t63 * t90;
t105 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t104 = g(1) * t29 - g(2) * t120 - g(3) * t53;
t103 = g(1) * t30 + g(2) * t26 + g(3) * t54;
t102 = g(1) * t49 - g(2) * t48 + g(3) * t63;
t101 = g(1) * t50 + g(2) * t47 + g(3) * t62;
t74 = pkin(2) * t124;
t100 = -pkin(2) * t131 + t47 * pkin(3) + pkin(9) * t48 + t74;
t99 = -t48 * pkin(3) + t47 * pkin(9) + t110;
t98 = t143 * t47 + t100;
t97 = t30 * pkin(4) + t29 * pkin(10) + t107;
t96 = t69 * pkin(2) + t50 * pkin(3) - t49 * pkin(9);
t95 = -pkin(4) * t26 + pkin(10) * t120 + t99;
t94 = t143 * t50 + t96;
t71 = t113 * t84;
t70 = -t85 * t131 + t127;
t68 = -t85 * t128 - t130;
t67 = -t124 + t131;
t61 = -g(3) * t85 - t112 * t84;
t32 = t62 * t129 - t63 * t86;
t18 = t50 * t129 - t49 * t86;
t16 = t47 * t129 + t48 * t86;
t7 = t101 * t87;
t5 = t104 * t90;
t4 = t104 * t86;
t3 = g(1) * t10 - g(2) * t14;
t2 = -g(1) * t18 - g(2) * t16 - g(3) * t32;
t6 = [0, 0, 0, 0, 0, 0, t112, t113, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t70, -g(1) * t67 - g(2) * t69, -t71, -g(1) * (-t89 * pkin(1) + pkin(8) * t134) - g(2) * (t93 * pkin(1) + pkin(8) * t136) 0, 0, 0, 0, 0, 0, g(1) * t48 + g(2) * t49, t114, -t71, -g(1) * t110 - g(2) * t119, 0, 0, 0, 0, 0, 0, g(1) * t26 - g(2) * t30, t115, -t114, -g(1) * t99 - g(2) * t107, 0, 0, 0, 0, 0, 0, t3, t116, -t115, -g(1) * t95 - g(2) * t97, 0, 0, 0, 0, 0, 0, t3, -t115, -t116, -g(1) * (-pkin(5) * t10 - qJ(6) * t9 + t95) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t97); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t67 + t142, g(3) * t84 * t88 + g(1) * t70 - g(2) * t68, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t102, 0, -g(2) * t74 + (g(2) * t131 + t142) * pkin(2), 0, 0, 0, 0, 0, 0, -t101 * t91, t7, t102, -g(1) * t96 - g(2) * t100 - g(3) * t117, 0, 0, 0, 0, 0, 0, t2, t105, -t7, -g(1) * t94 - g(2) * t98 - g(3) * t109, 0, 0, 0, 0, 0, 0, t2, -t7, -t105, -g(1) * (t18 * pkin(5) + t17 * qJ(6) + t94) - g(2) * (t16 * pkin(5) + t15 * qJ(6) + t98) - g(3) * (t32 * pkin(5) + t31 * qJ(6) + t109); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t103, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t103, -g(1) * t122 - g(2) * t123 - g(3) * t121, 0, 0, 0, 0, 0, 0, t5, -t103, t4, -g(1) * (-t111 * t29 + t122) - g(2) * (t111 * t120 + t123) - g(3) * (t111 * t53 + t121); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t106, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t106, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t19 * pkin(5) + t20 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
