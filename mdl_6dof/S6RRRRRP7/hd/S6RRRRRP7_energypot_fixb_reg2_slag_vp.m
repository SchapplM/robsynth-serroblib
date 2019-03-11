% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:43:08
% EndTime: 2019-03-10 01:43:08
% DurationCPUTime: 0.22s
% Computational Cost: add. (319->95), mult. (533->134), div. (0->0), fcn. (633->12), ass. (0->54)
t112 = cos(pkin(6));
t144 = t112 * pkin(8) + pkin(7);
t114 = sin(qJ(5));
t120 = cos(qJ(2));
t121 = cos(qJ(1));
t132 = t121 * t120;
t116 = sin(qJ(2));
t117 = sin(qJ(1));
t135 = t117 * t116;
t91 = -t112 * t132 + t135;
t143 = t91 * t114;
t133 = t121 * t116;
t134 = t117 * t120;
t93 = t112 * t134 + t133;
t142 = t93 * t114;
t111 = sin(pkin(6));
t141 = t111 * t116;
t140 = t111 * t117;
t119 = cos(qJ(3));
t139 = t111 * t119;
t138 = t111 * t120;
t137 = t111 * t121;
t115 = sin(qJ(3));
t136 = t112 * t115;
t131 = t121 * pkin(1) + pkin(8) * t140;
t130 = t114 * t138;
t129 = t115 * t140;
t108 = t117 * pkin(1);
t128 = -pkin(8) * t137 + t108;
t104 = t119 * pkin(3) + pkin(2);
t122 = -pkin(10) - pkin(9);
t127 = pkin(3) * t136 + t104 * t141 + t122 * t138 + t144;
t94 = -t112 * t135 + t132;
t126 = pkin(3) * t129 + t94 * t104 - t93 * t122 + t131;
t125 = g(1) * t117 - g(2) * t121;
t110 = qJ(3) + qJ(4);
t105 = sin(t110);
t106 = cos(t110);
t92 = t112 * t133 + t134;
t81 = t92 * t105 + t106 * t137;
t83 = t94 * t105 - t106 * t140;
t87 = t105 * t141 - t112 * t106;
t124 = g(1) * t83 + g(2) * t81 + g(3) * t87;
t80 = -g(1) * t93 - g(2) * t91 + g(3) * t138;
t123 = t108 + t92 * t104 - t91 * t122 + (-pkin(3) * t115 - pkin(8)) * t137;
t118 = cos(qJ(5));
t113 = -qJ(6) - pkin(11);
t103 = t118 * pkin(5) + pkin(4);
t88 = t112 * t105 + t106 * t141;
t84 = t105 * t140 + t94 * t106;
t82 = -t105 * t137 + t92 * t106;
t78 = -g(1) * (t84 * t118 + t142) - g(2) * (t82 * t118 + t143) - g(3) * (t88 * t118 - t130);
t77 = -g(1) * (-t84 * t114 + t93 * t118) - g(2) * (-t82 * t114 + t91 * t118) - g(3) * (-t88 * t114 - t118 * t138);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t117, t125, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t141, -t80, -g(3) * t112 - t125 * t111, -g(1) * t131 - g(2) * t128 - g(3) * t144, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t119 + t129) - g(2) * (-t115 * t137 + t92 * t119) - g(3) * (t116 * t139 + t136) -g(1) * (-t94 * t115 + t117 * t139) - g(2) * (-t92 * t115 - t119 * t137) - g(3) * (t112 * t119 - t115 * t141) t80, -g(1) * (t94 * pkin(2) + t93 * pkin(9) + t131) - g(2) * (t92 * pkin(2) + t91 * pkin(9) + t128) - g(3) * ((pkin(2) * t116 - pkin(9) * t120) * t111 + t144) 0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82 - g(3) * t88, t124, t80, -g(1) * t126 - g(2) * t123 - g(3) * t127, 0, 0, 0, 0, 0, 0, t78, t77, -t124, -g(1) * (t84 * pkin(4) + t83 * pkin(11) + t126) - g(2) * (t82 * pkin(4) + t81 * pkin(11) + t123) - g(3) * (t88 * pkin(4) + t87 * pkin(11) + t127) 0, 0, 0, 0, 0, 0, t78, t77, -t124, -g(1) * (pkin(5) * t142 + t84 * t103 - t83 * t113 + t126) - g(2) * (pkin(5) * t143 + t82 * t103 - t81 * t113 + t123) - g(3) * (-pkin(5) * t130 + t88 * t103 - t87 * t113 + t127);];
U_reg  = t1;
