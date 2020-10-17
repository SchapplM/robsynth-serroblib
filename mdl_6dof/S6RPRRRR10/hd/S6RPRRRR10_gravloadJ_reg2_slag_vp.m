% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:49:25
% EndTime: 2019-05-06 04:49:29
% DurationCPUTime: 1.46s
% Computational Cost: add. (1197->143), mult. (2892->233), div. (0->0), fcn. (3718->16), ass. (0->88)
t141 = cos(qJ(3));
t121 = sin(pkin(7));
t122 = sin(pkin(6));
t103 = t122 * t121;
t124 = cos(pkin(7));
t80 = cos(qJ(1));
t123 = cos(pkin(13));
t125 = cos(pkin(6));
t106 = t125 * t123;
t120 = sin(pkin(13));
t140 = sin(qJ(1));
t94 = -t80 * t106 + t140 * t120;
t148 = t80 * t103 + t94 * t124;
t105 = t125 * t120;
t59 = t80 * t105 + t140 * t123;
t77 = sin(qJ(3));
t35 = -t59 * t141 + t148 * t77;
t104 = t124 * t122;
t49 = -t80 * t104 + t94 * t121;
t74 = qJ(4) + qJ(5);
t71 = sin(t74);
t72 = cos(t74);
t16 = t35 * t72 - t49 * t71;
t32 = t148 * t141 + t59 * t77;
t75 = sin(qJ(6));
t78 = cos(qJ(6));
t154 = t16 * t75 + t32 * t78;
t153 = t16 * t78 - t32 * t75;
t76 = sin(qJ(4));
t136 = t49 * t76;
t79 = cos(qJ(4));
t152 = -t35 * t79 + t136;
t147 = t35 * t71 + t49 * t72;
t145 = t35 * t76 + t49 * t79;
t88 = t140 * t106 + t80 * t120;
t50 = t140 * t104 + t88 * t121;
t102 = t122 * t120;
t142 = t123 * t104 + t121 * t125;
t47 = t141 * t102 + t142 * t77;
t58 = -t123 * t103 + t125 * t124;
t146 = -t47 * t76 + t58 * t79;
t144 = -t140 * t103 + t88 * t124;
t60 = -t140 * t105 + t80 * t123;
t37 = t60 * t141 - t144 * t77;
t19 = -t37 * t76 + t50 * t79;
t134 = t50 * t76;
t131 = t72 * t75;
t130 = t72 * t78;
t70 = t79 * pkin(4) + pkin(3);
t81 = -pkin(11) - pkin(10);
t129 = -t32 * t70 + t35 * t81;
t36 = t144 * t141 + t60 * t77;
t128 = -t36 * t70 - t37 * t81;
t46 = t77 * t102 - t142 * t141;
t127 = -t46 * t70 - t47 * t81;
t107 = t122 * t140;
t126 = t80 * pkin(1) + qJ(2) * t107;
t119 = pkin(5) * t147 - pkin(12) * t16;
t17 = t37 * t71 - t50 * t72;
t18 = t37 * t72 + t50 * t71;
t118 = -t17 * pkin(5) + t18 * pkin(12);
t26 = -t47 * t71 + t58 * t72;
t27 = t47 * t72 + t58 * t71;
t117 = t26 * pkin(5) + t27 * pkin(12);
t116 = t80 * t122;
t115 = t145 * pkin(4);
t114 = t19 * pkin(4);
t113 = t146 * pkin(4);
t112 = -t140 * pkin(1) + qJ(2) * t116;
t111 = -pkin(5) * t72 - pkin(12) * t71;
t110 = g(1) * t147 + g(2) * t17;
t109 = -g(1) * t32 + g(2) * t36;
t3 = g(1) * t17 - g(2) * t147 - g(3) * t26;
t5 = g(1) * t18 - g(2) * t16 + g(3) * t27;
t98 = g(1) * t36 + g(2) * t32 + g(3) * t46;
t97 = g(1) * t37 - g(2) * t35 + g(3) * t47;
t85 = -t59 * pkin(2) - pkin(9) * t49 + t112;
t84 = t60 * pkin(2) + pkin(9) * t50 + t126;
t83 = -pkin(4) * t136 + t32 * t81 + t35 * t70 + t85;
t82 = pkin(4) * t134 - t36 * t81 + t37 * t70 + t84;
t55 = -g(1) * t107 + g(2) * t116 - g(3) * t125;
t20 = t37 * t79 + t134;
t8 = t18 * t78 + t36 * t75;
t7 = -t18 * t75 + t36 * t78;
t6 = t98 * t71;
t2 = t3 * t78;
t1 = t3 * t75;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t140 - g(2) * t80, g(1) * t80 + g(2) * t140, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * t60, -g(1) * t94 + g(2) * t88, -g(1) * t116 - g(2) * t107, -g(1) * t112 - g(2) * t126, 0, 0, 0, 0, 0, 0, -g(1) * t35 - g(2) * t37, t109, g(1) * t49 - g(2) * t50, -g(1) * t85 - g(2) * t84, 0, 0, 0, 0, 0, 0, g(1) * t152 - g(2) * t20, g(1) * t145 - g(2) * t19, -t109, -g(1) * (t35 * pkin(3) - pkin(10) * t32 + t85) - g(2) * (t37 * pkin(3) + t36 * pkin(10) + t84) 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, t110, -t109, -g(1) * t83 - g(2) * t82, 0, 0, 0, 0, 0, 0, -g(1) * t153 - g(2) * t8, g(1) * t154 - g(2) * t7, -t110, -g(1) * (t16 * pkin(5) + pkin(12) * t147 + t83) - g(2) * (t18 * pkin(5) + t17 * pkin(12) + t82); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t97, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t79, -t98 * t76, -t97, -g(1) * (-t36 * pkin(3) + t37 * pkin(10)) - g(2) * (-t32 * pkin(3) - pkin(10) * t35) - g(3) * (-t46 * pkin(3) + t47 * pkin(10)) 0, 0, 0, 0, 0, 0, t98 * t72, -t6, -t97, -g(1) * t128 - g(2) * t129 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (-t36 * t130 + t37 * t75) - g(2) * (-t32 * t130 - t35 * t75) - g(3) * (-t46 * t130 + t47 * t75) -g(1) * (t36 * t131 + t37 * t78) - g(2) * (t32 * t131 - t35 * t78) - g(3) * (t46 * t131 + t47 * t78) t6, -g(1) * (t111 * t36 + t128) - g(2) * (t111 * t32 + t129) - g(3) * (t111 * t46 + t127); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t145 - g(3) * t146, g(1) * t20 + g(2) * t152 - g(3) * (-t47 * t79 - t58 * t76) 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, -g(1) * t114 - g(2) * t115 - g(3) * t113, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * (t114 + t118) - g(2) * (t115 + t119) - g(3) * (t113 + t117); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t118 - g(2) * t119 - g(3) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t154 - g(3) * (-t27 * t75 + t46 * t78) g(1) * t8 - g(2) * t153 - g(3) * (-t27 * t78 - t46 * t75) 0, 0;];
taug_reg  = t4;
