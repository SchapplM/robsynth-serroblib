% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:10:37
% EndTime: 2019-05-05 09:10:40
% DurationCPUTime: 0.81s
% Computational Cost: add. (1013->174), mult. (2852->257), div. (0->0), fcn. (3658->14), ass. (0->104)
t73 = sin(qJ(4));
t125 = qJ(5) * t73;
t134 = cos(qJ(3));
t123 = cos(pkin(7));
t71 = sin(pkin(7));
t122 = cos(pkin(12));
t124 = cos(pkin(6));
t102 = t124 * t122;
t120 = sin(pkin(12));
t135 = cos(qJ(2));
t75 = sin(qJ(2));
t83 = -t135 * t102 + t120 * t75;
t121 = sin(pkin(6));
t99 = t122 * t121;
t142 = t83 * t123 + t71 * t99;
t61 = t75 * t102 + t120 * t135;
t74 = sin(qJ(3));
t27 = t142 * t134 + t61 * t74;
t77 = cos(qJ(4));
t133 = t27 * t77;
t145 = -pkin(4) * t133 - t27 * t125;
t101 = t124 * t120;
t84 = t135 * t101 + t122 * t75;
t98 = t121 * t120;
t141 = t84 * t123 - t71 * t98;
t62 = -t75 * t101 + t122 * t135;
t29 = t141 * t134 + t62 * t74;
t132 = t29 * t77;
t144 = -pkin(4) * t132 - t29 * t125;
t111 = t124 * t71;
t112 = t75 * t121;
t100 = t123 * t121;
t88 = t134 * t100;
t47 = -t134 * t111 + t74 * t112 - t135 * t88;
t131 = t47 * t77;
t143 = -pkin(4) * t131 - t47 * t125;
t140 = pkin(5) + pkin(10);
t139 = pkin(9) * t71;
t106 = t123 * t134;
t35 = t61 * t106 - t83 * t74;
t138 = t35 * pkin(10);
t37 = t62 * t106 - t84 * t74;
t137 = t37 * pkin(10);
t105 = t121 * t135;
t56 = t74 * t105 + t75 * t88;
t136 = t56 * pkin(10);
t130 = t71 * t73;
t129 = t71 * t77;
t72 = sin(qJ(6));
t128 = t72 * t73;
t76 = cos(qJ(6));
t127 = t73 * t76;
t107 = t71 * t112;
t126 = pkin(2) * t105 + pkin(9) * t107;
t57 = -t75 * t74 * t100 + t134 * t105;
t119 = t57 * pkin(3) + t126;
t24 = t27 * pkin(3);
t28 = t61 * t134 - t142 * t74;
t118 = t28 * pkin(10) - t24;
t25 = t29 * pkin(3);
t30 = t62 * t134 - t141 * t74;
t117 = t30 * pkin(10) - t25;
t46 = t47 * pkin(3);
t48 = t134 * t112 + (t135 * t100 + t111) * t74;
t116 = t48 * pkin(10) - t46;
t78 = -t123 * t99 + t83 * t71;
t10 = t28 * t73 - t78 * t77;
t11 = t28 * t77 + t78 * t73;
t115 = -t10 * pkin(4) + t11 * qJ(5);
t79 = t123 * t98 + t84 * t71;
t12 = t30 * t73 - t79 * t77;
t13 = t30 * t77 + t79 * t73;
t114 = -t12 * pkin(4) + t13 * qJ(5);
t113 = t74 * t123;
t82 = -t71 * t105 + t124 * t123;
t31 = t48 * t73 - t82 * t77;
t32 = t48 * t77 + t82 * t73;
t110 = -t31 * pkin(4) + t32 * qJ(5);
t109 = -t83 * pkin(2) + t61 * t139;
t108 = -t84 * pkin(2) + t62 * t139;
t36 = -t61 * t113 - t83 * t134;
t104 = t36 * pkin(3) + t109;
t38 = -t62 * t113 - t84 * t134;
t103 = t38 * pkin(3) + t108;
t40 = -t77 * t107 + t57 * t73;
t41 = t73 * t107 + t57 * t77;
t95 = t41 * pkin(4) + t40 * qJ(5) + t119;
t2 = g(1) * t12 + g(2) * t10 + g(3) * t31;
t94 = g(1) * t13 + g(2) * t11 + g(3) * t32;
t16 = -t61 * t129 + t36 * t73;
t18 = -t62 * t129 + t38 * t73;
t93 = g(1) * t18 + g(2) * t16 + g(3) * t40;
t17 = t61 * t130 + t36 * t77;
t19 = t62 * t130 + t38 * t77;
t92 = g(1) * t19 + g(2) * t17 + g(3) * t41;
t91 = g(1) * t29 + g(2) * t27 + g(3) * t47;
t90 = g(1) * t30 + g(2) * t28 + g(3) * t48;
t89 = g(1) * t37 + g(2) * t35 + g(3) * t56;
t87 = t17 * pkin(4) + t16 * qJ(5) + t104;
t86 = t19 * pkin(4) + t18 * qJ(5) + t103;
t85 = -g(1) * t62 - g(2) * t61 - g(3) * t112;
t5 = t91 * t77;
t4 = t91 * t73;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t84 + g(2) * t83 - g(3) * t105, -t85, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t38 - g(2) * t36 - g(3) * t57, t89, t85 * t71, -g(1) * t108 - g(2) * t109 - g(3) * t126, 0, 0, 0, 0, 0, 0, -t92, t93, -t89, -g(1) * (t103 + t137) - g(2) * (t104 + t138) - g(3) * (t119 + t136) 0, 0, 0, 0, 0, 0, -t89, t92, -t93, -g(1) * (t86 + t137) - g(2) * (t87 + t138) - g(3) * (t95 + t136) 0, 0, 0, 0, 0, 0, -g(1) * (t18 * t72 + t37 * t76) - g(2) * (t16 * t72 + t35 * t76) - g(3) * (t40 * t72 + t56 * t76) -g(1) * (t18 * t76 - t37 * t72) - g(2) * (t16 * t76 - t35 * t72) - g(3) * (t40 * t76 - t56 * t72) -t92, -g(1) * (t19 * pkin(11) + t140 * t37 + t86) - g(2) * (t17 * pkin(11) + t140 * t35 + t87) - g(3) * (t41 * pkin(11) + t140 * t56 + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t90, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t90, -g(1) * t117 - g(2) * t118 - g(3) * t116, 0, 0, 0, 0, 0, 0, -t90, -t5, t4, -g(1) * (t117 + t144) - g(2) * (t118 + t145) - g(3) * (t116 + t143) 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t128 + t30 * t76) - g(2) * (-t27 * t128 + t28 * t76) - g(3) * (-t47 * t128 + t48 * t76) -g(1) * (-t29 * t127 - t30 * t72) - g(2) * (-t27 * t127 - t28 * t72) - g(3) * (-t47 * t127 - t48 * t72) t5, -g(1) * (-pkin(11) * t132 + t140 * t30 + t144 - t25) - g(2) * (-pkin(11) * t133 + t140 * t28 + t145 - t24) - g(3) * (-pkin(11) * t131 + t140 * t48 + t143 - t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t94, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t94, -g(1) * t114 - g(2) * t115 - g(3) * t110, 0, 0, 0, 0, 0, 0, -t94 * t72, -t94 * t76, t2, -g(1) * (-t12 * pkin(11) + t114) - g(2) * (-t10 * pkin(11) + t115) - g(3) * (-t31 * pkin(11) + t110); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t76 - t29 * t72) - g(2) * (t10 * t76 - t27 * t72) - g(3) * (t31 * t76 - t47 * t72) -g(1) * (-t12 * t72 - t29 * t76) - g(2) * (-t10 * t72 - t27 * t76) - g(3) * (-t31 * t72 - t47 * t76) 0, 0;];
taug_reg  = t1;
