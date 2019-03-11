% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:41:04
% EndTime: 2019-03-10 04:41:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (311->108), mult. (592->154), div. (0->0), fcn. (717->14), ass. (0->52)
t112 = sin(qJ(4));
t131 = t112 * pkin(4) + pkin(9);
t119 = -pkin(11) - pkin(10);
t110 = qJ(4) + qJ(5);
t102 = sin(t110);
t143 = pkin(5) * t102 + t131;
t142 = cos(qJ(3));
t138 = cos(pkin(6));
t140 = t138 * pkin(8) + pkin(7);
t116 = cos(qJ(4));
t101 = t116 * pkin(4) + pkin(3);
t118 = cos(qJ(1));
t111 = sin(pkin(6));
t115 = sin(qJ(1));
t136 = t111 * t115;
t139 = t118 * pkin(1) + pkin(8) * t136;
t114 = sin(qJ(2));
t137 = t111 * t114;
t117 = cos(qJ(2));
t135 = t111 * t117;
t134 = t111 * t118;
t133 = pkin(2) * t137 + t140;
t129 = t115 * t138;
t91 = -t114 * t129 + t118 * t117;
t132 = t91 * pkin(2) + t139;
t130 = t111 * t142;
t128 = t118 * t138;
t127 = t115 * pkin(1) - pkin(8) * t134;
t126 = g(1) * t115 - g(2) * t118;
t90 = t118 * t114 + t117 * t129;
t125 = t90 * pkin(9) + t132;
t89 = t114 * t128 + t115 * t117;
t124 = t89 * pkin(2) + t127;
t123 = -pkin(9) * t135 + t133;
t113 = sin(qJ(3));
t80 = t89 * t113 + t118 * t130;
t82 = t91 * t113 - t115 * t130;
t86 = t113 * t137 - t138 * t142;
t122 = g(1) * t82 + g(2) * t80 + g(3) * t86;
t88 = t115 * t114 - t117 * t128;
t121 = t88 * pkin(9) + t124;
t120 = -g(1) * t90 - g(2) * t88 + g(3) * t135;
t109 = -pkin(12) + t119;
t105 = qJ(6) + t110;
t103 = cos(t110);
t100 = cos(t105);
t99 = sin(t105);
t92 = pkin(5) * t103 + t101;
t87 = t113 * t138 + t114 * t130;
t83 = t113 * t136 + t142 * t91;
t81 = -t113 * t134 + t142 * t89;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t118 - g(2) * t115, t126, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89 - g(3) * t137, -t120, -g(3) * t138 - t111 * t126, -g(1) * t139 - g(2) * t127 - g(3) * t140, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t87, t122, t120, -g(1) * t125 - g(2) * t121 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t112 + t83 * t116) - g(2) * (t88 * t112 + t81 * t116) - g(3) * (-t112 * t135 + t87 * t116) -g(1) * (-t83 * t112 + t90 * t116) - g(2) * (-t81 * t112 + t88 * t116) - g(3) * (-t87 * t112 - t116 * t135) -t122, -g(1) * (t83 * pkin(3) + t82 * pkin(10) + t125) - g(2) * (t81 * pkin(3) + t80 * pkin(10) + t121) - g(3) * (t87 * pkin(3) + t86 * pkin(10) + t123) 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t102 + t83 * t103) - g(2) * (t88 * t102 + t81 * t103) - g(3) * (-t102 * t135 + t87 * t103) -g(1) * (-t83 * t102 + t90 * t103) - g(2) * (-t81 * t102 + t88 * t103) - g(3) * (-t87 * t102 - t103 * t135) -t122, -g(1) * (t83 * t101 - t82 * t119 + t131 * t90 + t132) - g(2) * (t81 * t101 - t80 * t119 + t131 * t88 + t124) - g(3) * (t87 * t101 - t86 * t119 - t131 * t135 + t133) 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t100 + t90 * t99) - g(2) * (t81 * t100 + t88 * t99) - g(3) * (t87 * t100 - t135 * t99) -g(1) * (t90 * t100 - t83 * t99) - g(2) * (t88 * t100 - t81 * t99) - g(3) * (-t100 * t135 - t87 * t99) -t122, -g(1) * (-t82 * t109 + t143 * t90 + t83 * t92 + t132) - g(2) * (-t80 * t109 + t143 * t88 + t81 * t92 + t124) - g(3) * (-t86 * t109 - t143 * t135 + t87 * t92 + t133);];
U_reg  = t1;
