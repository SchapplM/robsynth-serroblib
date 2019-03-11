% Calculate inertial parameters regressor of potential energy for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:23
% EndTime: 2019-03-08 18:54:23
% DurationCPUTime: 0.34s
% Computational Cost: add. (527->100), mult. (1376->154), div. (0->0), fcn. (1753->14), ass. (0->64)
t107 = sin(pkin(7));
t111 = cos(pkin(7));
t112 = cos(pkin(6));
t108 = sin(pkin(6));
t109 = cos(pkin(12));
t145 = t108 * t109;
t129 = t107 * t145 - t112 * t111;
t106 = sin(pkin(11));
t143 = t108 * t111;
t105 = sin(pkin(12));
t110 = cos(pkin(11));
t146 = t106 * t112;
t92 = -t110 * t105 - t109 * t146;
t131 = t106 * t143 - t92 * t107;
t152 = cos(qJ(3));
t151 = cos(qJ(4));
t147 = t106 * t108;
t149 = t110 * pkin(1) + qJ(2) * t147;
t148 = t105 * t108;
t144 = t108 * t110;
t142 = t110 * t112;
t140 = t112 * qJ(2) + qJ(1);
t114 = sin(qJ(5));
t137 = pkin(5) * t114 + pkin(9);
t136 = t107 * t152;
t135 = t111 * t152;
t134 = t108 * t136;
t133 = g(1) * t106 - g(2) * t110;
t132 = t106 * pkin(1) - qJ(2) * t144;
t90 = -t106 * t105 + t109 * t142;
t130 = t90 * t107 + t110 * t143;
t115 = sin(qJ(4));
t116 = sin(qJ(3));
t91 = t105 * t142 + t106 * t109;
t75 = t91 * t152 + (-t107 * t144 + t111 * t90) * t116;
t68 = t75 * t115 + t130 * t151;
t93 = -t105 * t146 + t110 * t109;
t77 = t93 * t152 + (t107 * t147 + t111 * t92) * t116;
t70 = t77 * t115 - t131 * t151;
t84 = t112 * t107 * t116 + (t109 * t111 * t116 + t152 * t105) * t108;
t78 = t84 * t115 + t129 * t151;
t128 = g(1) * t70 + g(2) * t68 + g(3) * t78;
t74 = t110 * t134 + t91 * t116 - t90 * t135;
t76 = -t106 * t134 + t93 * t116 - t92 * t135;
t83 = -t112 * t136 + t116 * t148 - t135 * t145;
t127 = g(1) * t76 + g(2) * t74 + g(3) * t83;
t126 = t93 * pkin(2) + t131 * pkin(8) + t149;
t125 = t77 * pkin(3) + t126;
t124 = pkin(2) * t148 - t129 * pkin(8) + t140;
t123 = t84 * pkin(3) + t124;
t122 = t76 * pkin(9) + t125;
t121 = t83 * pkin(9) + t123;
t120 = t91 * pkin(2) - t130 * pkin(8) + t132;
t119 = t75 * pkin(3) + t120;
t118 = t74 * pkin(9) + t119;
t117 = cos(qJ(5));
t113 = -qJ(6) - pkin(10);
t101 = t117 * pkin(5) + pkin(4);
t79 = -t129 * t115 + t84 * t151;
t71 = t131 * t115 + t77 * t151;
t69 = -t130 * t115 + t75 * t151;
t66 = -g(1) * (t76 * t114 + t71 * t117) - g(2) * (t74 * t114 + t69 * t117) - g(3) * (t83 * t114 + t79 * t117);
t65 = -g(1) * (-t71 * t114 + t76 * t117) - g(2) * (-t69 * t114 + t74 * t117) - g(3) * (-t79 * t114 + t83 * t117);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t110 - g(2) * t106, t133, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t148, -g(1) * t92 - g(2) * t90 - g(3) * t145, -g(3) * t112 - t133 * t108, -g(1) * t149 - g(2) * t132 - g(3) * t140, 0, 0, 0, 0, 0, 0, -g(1) * t77 - g(2) * t75 - g(3) * t84, t127, -g(1) * t131 + g(2) * t130 + g(3) * t129, -g(1) * t126 - g(2) * t120 - g(3) * t124, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t79, t128, -t127, -g(1) * t122 - g(2) * t118 - g(3) * t121, 0, 0, 0, 0, 0, 0, t66, t65, -t128, -g(1) * (t71 * pkin(4) + t70 * pkin(10) + t122) - g(2) * (t69 * pkin(4) + t68 * pkin(10) + t118) - g(3) * (t79 * pkin(4) + t78 * pkin(10) + t121) 0, 0, 0, 0, 0, 0, t66, t65, -t128, -g(1) * (t71 * t101 - t70 * t113 + t137 * t76 + t125) - g(2) * (t69 * t101 - t68 * t113 + t137 * t74 + t119) - g(3) * (t79 * t101 - t78 * t113 + t137 * t83 + t123);];
U_reg  = t1;
