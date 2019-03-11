% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:53:53
% EndTime: 2019-03-09 15:53:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (308->91), mult. (511->130), div. (0->0), fcn. (603->12), ass. (0->51)
t112 = cos(pkin(6));
t144 = t112 * pkin(8) + pkin(7);
t111 = sin(pkin(6));
t116 = sin(qJ(2));
t143 = t111 * t116;
t117 = sin(qJ(1));
t142 = t111 * t117;
t119 = cos(qJ(3));
t141 = t111 * t119;
t120 = cos(qJ(2));
t140 = t111 * t120;
t121 = cos(qJ(1));
t139 = t111 * t121;
t115 = sin(qJ(3));
t138 = t112 * t115;
t137 = t117 * t116;
t136 = t117 * t120;
t135 = t121 * t116;
t134 = t121 * t120;
t133 = t121 * pkin(1) + pkin(8) * t142;
t132 = t115 * t142;
t108 = t117 * pkin(1);
t131 = -pkin(8) * t139 + t108;
t104 = t119 * pkin(3) + pkin(2);
t113 = -qJ(4) - pkin(9);
t130 = pkin(3) * t138 + t104 * t143 + t113 * t140 + t144;
t93 = t112 * t136 + t135;
t94 = -t112 * t137 + t134;
t129 = pkin(3) * t132 + t94 * t104 - t93 * t113 + t133;
t128 = g(1) * t117 - g(2) * t121;
t110 = qJ(3) + pkin(11);
t105 = sin(t110);
t106 = cos(t110);
t92 = t112 * t135 + t136;
t80 = t92 * t105 + t106 * t139;
t82 = t94 * t105 - t106 * t142;
t87 = t105 * t143 - t112 * t106;
t127 = g(1) * t82 + g(2) * t80 + g(3) * t87;
t81 = -t105 * t139 + t92 * t106;
t83 = t105 * t142 + t94 * t106;
t88 = t112 * t105 + t106 * t143;
t126 = g(1) * t83 + g(2) * t81 + g(3) * t88;
t91 = -t112 * t134 + t137;
t77 = -g(1) * t93 - g(2) * t91 + g(3) * t140;
t125 = t88 * pkin(4) + t87 * qJ(5) + t130;
t124 = t83 * pkin(4) + t82 * qJ(5) + t129;
t123 = t108 + t92 * t104 - t91 * t113 + (-pkin(3) * t115 - pkin(8)) * t139;
t122 = t81 * pkin(4) + t80 * qJ(5) + t123;
t118 = cos(qJ(6));
t114 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t117, t128, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t143, -t77, -g(3) * t112 - t128 * t111, -g(1) * t133 - g(2) * t131 - g(3) * t144, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t119 + t132) - g(2) * (-t115 * t139 + t92 * t119) - g(3) * (t116 * t141 + t138) -g(1) * (-t94 * t115 + t117 * t141) - g(2) * (-t92 * t115 - t119 * t139) - g(3) * (t112 * t119 - t115 * t143) t77, -g(1) * (t94 * pkin(2) + t93 * pkin(9) + t133) - g(2) * (t92 * pkin(2) + t91 * pkin(9) + t131) - g(3) * ((pkin(2) * t116 - pkin(9) * t120) * t111 + t144) 0, 0, 0, 0, 0, 0, -t126, t127, t77, -g(1) * t129 - g(2) * t123 - g(3) * t130, 0, 0, 0, 0, 0, 0, t77, t126, -t127, -g(1) * t124 - g(2) * t122 - g(3) * t125, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t114 + t93 * t118) - g(2) * (t80 * t114 + t91 * t118) - g(3) * (t87 * t114 - t118 * t140) -g(1) * (-t93 * t114 + t82 * t118) - g(2) * (-t91 * t114 + t80 * t118) - g(3) * (t114 * t140 + t87 * t118) -t126, -g(1) * (t93 * pkin(5) + t83 * pkin(10) + t124) - g(2) * (t91 * pkin(5) + t81 * pkin(10) + t122) - g(3) * (-pkin(5) * t140 + t88 * pkin(10) + t125);];
U_reg  = t1;
