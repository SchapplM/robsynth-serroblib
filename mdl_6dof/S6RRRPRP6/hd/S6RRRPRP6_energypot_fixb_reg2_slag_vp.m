% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:07
% EndTime: 2019-03-09 17:01:08
% DurationCPUTime: 0.24s
% Computational Cost: add. (319->95), mult. (533->134), div. (0->0), fcn. (633->12), ass. (0->54)
t111 = cos(pkin(6));
t143 = t111 * pkin(8) + pkin(7);
t114 = sin(qJ(5));
t120 = cos(qJ(2));
t121 = cos(qJ(1));
t131 = t121 * t120;
t116 = sin(qJ(2));
t117 = sin(qJ(1));
t134 = t117 * t116;
t90 = -t111 * t131 + t134;
t142 = t90 * t114;
t132 = t121 * t116;
t133 = t117 * t120;
t92 = t111 * t133 + t132;
t141 = t92 * t114;
t110 = sin(pkin(6));
t140 = t110 * t116;
t139 = t110 * t117;
t119 = cos(qJ(3));
t138 = t110 * t119;
t137 = t110 * t120;
t136 = t110 * t121;
t115 = sin(qJ(3));
t135 = t111 * t115;
t130 = t121 * pkin(1) + pkin(8) * t139;
t129 = t114 * t137;
t128 = t115 * t139;
t107 = t117 * pkin(1);
t127 = -pkin(8) * t136 + t107;
t103 = t119 * pkin(3) + pkin(2);
t113 = -qJ(4) - pkin(9);
t126 = pkin(3) * t135 + t103 * t140 + t113 * t137 + t143;
t93 = -t111 * t134 + t131;
t125 = pkin(3) * t128 + t93 * t103 - t92 * t113 + t130;
t124 = g(1) * t117 - g(2) * t121;
t109 = qJ(3) + pkin(11);
t104 = sin(t109);
t105 = cos(t109);
t91 = t111 * t132 + t133;
t80 = t91 * t104 + t105 * t136;
t82 = t93 * t104 - t105 * t139;
t86 = t104 * t140 - t111 * t105;
t123 = g(1) * t82 + g(2) * t80 + g(3) * t86;
t79 = -g(1) * t92 - g(2) * t90 + g(3) * t137;
t122 = t107 + t91 * t103 - t90 * t113 + (-pkin(3) * t115 - pkin(8)) * t136;
t118 = cos(qJ(5));
t112 = -qJ(6) - pkin(10);
t102 = t118 * pkin(5) + pkin(4);
t87 = t111 * t104 + t105 * t140;
t83 = t104 * t139 + t93 * t105;
t81 = -t104 * t136 + t91 * t105;
t77 = -g(1) * (t83 * t118 + t141) - g(2) * (t81 * t118 + t142) - g(3) * (t87 * t118 - t129);
t76 = -g(1) * (-t83 * t114 + t92 * t118) - g(2) * (-t81 * t114 + t90 * t118) - g(3) * (-t87 * t114 - t118 * t137);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t117, t124, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t140, -t79, -g(3) * t111 - t124 * t110, -g(1) * t130 - g(2) * t127 - g(3) * t143, 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t119 + t128) - g(2) * (-t115 * t136 + t91 * t119) - g(3) * (t116 * t138 + t135) -g(1) * (-t93 * t115 + t117 * t138) - g(2) * (-t91 * t115 - t119 * t136) - g(3) * (t111 * t119 - t115 * t140) t79, -g(1) * (t93 * pkin(2) + t92 * pkin(9) + t130) - g(2) * (t91 * pkin(2) + t90 * pkin(9) + t127) - g(3) * ((pkin(2) * t116 - pkin(9) * t120) * t110 + t143) 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t87, t123, t79, -g(1) * t125 - g(2) * t122 - g(3) * t126, 0, 0, 0, 0, 0, 0, t77, t76, -t123, -g(1) * (t83 * pkin(4) + t82 * pkin(10) + t125) - g(2) * (t81 * pkin(4) + t80 * pkin(10) + t122) - g(3) * (t87 * pkin(4) + t86 * pkin(10) + t126) 0, 0, 0, 0, 0, 0, t77, t76, -t123, -g(1) * (pkin(5) * t141 + t83 * t102 - t82 * t112 + t125) - g(2) * (pkin(5) * t142 + t81 * t102 - t80 * t112 + t122) - g(3) * (-pkin(5) * t129 + t87 * t102 - t86 * t112 + t126);];
U_reg  = t1;
