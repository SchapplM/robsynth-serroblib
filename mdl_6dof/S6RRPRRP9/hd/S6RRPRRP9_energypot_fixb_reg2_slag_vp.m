% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:33:57
% EndTime: 2019-03-09 12:33:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (319->95), mult. (533->133), div. (0->0), fcn. (633->12), ass. (0->53)
t113 = cos(pkin(6));
t142 = t113 * pkin(8) + pkin(7);
t116 = sin(qJ(5));
t120 = cos(qJ(2));
t121 = cos(qJ(1));
t131 = t121 * t120;
t117 = sin(qJ(2));
t118 = sin(qJ(1));
t134 = t118 * t117;
t90 = -t113 * t131 + t134;
t141 = t90 * t116;
t132 = t121 * t117;
t133 = t118 * t120;
t92 = t113 * t133 + t132;
t140 = t92 * t116;
t111 = sin(pkin(6));
t139 = t111 * t117;
t138 = t111 * t118;
t137 = t111 * t120;
t136 = t111 * t121;
t110 = sin(pkin(11));
t135 = t113 * t110;
t130 = t121 * pkin(1) + pkin(8) * t138;
t129 = t116 * t137;
t128 = t110 * t138;
t107 = t118 * pkin(1);
t127 = -pkin(8) * t136 + t107;
t112 = cos(pkin(11));
t102 = t112 * pkin(3) + pkin(2);
t115 = -pkin(9) - qJ(3);
t126 = pkin(3) * t135 + t102 * t139 + t115 * t137 + t142;
t93 = -t113 * t134 + t131;
t125 = pkin(3) * t128 + t93 * t102 - t92 * t115 + t130;
t124 = g(1) * t118 - g(2) * t121;
t109 = pkin(11) + qJ(4);
t104 = sin(t109);
t105 = cos(t109);
t91 = t113 * t132 + t133;
t80 = t91 * t104 + t105 * t136;
t82 = t93 * t104 - t105 * t138;
t86 = t104 * t139 - t113 * t105;
t123 = g(1) * t82 + g(2) * t80 + g(3) * t86;
t79 = -g(1) * t92 - g(2) * t90 + g(3) * t137;
t122 = t107 + t91 * t102 - t90 * t115 + (-pkin(3) * t110 - pkin(8)) * t136;
t119 = cos(qJ(5));
t114 = -qJ(6) - pkin(10);
t103 = t119 * pkin(5) + pkin(4);
t87 = t113 * t104 + t105 * t139;
t83 = t104 * t138 + t93 * t105;
t81 = -t104 * t136 + t91 * t105;
t77 = -g(1) * (t83 * t119 + t140) - g(2) * (t81 * t119 + t141) - g(3) * (t87 * t119 - t129);
t76 = -g(1) * (-t83 * t116 + t92 * t119) - g(2) * (-t81 * t116 + t90 * t119) - g(3) * (-t87 * t116 - t119 * t137);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t118, t124, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t91 - g(3) * t139, -t79, -g(3) * t113 - t124 * t111, -g(1) * t130 - g(2) * t127 - g(3) * t142, 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t112 + t128) - g(2) * (-t110 * t136 + t91 * t112) - g(3) * (t112 * t139 + t135) -g(1) * (-t93 * t110 + t112 * t138) - g(2) * (-t91 * t110 - t112 * t136) - g(3) * (-t110 * t139 + t113 * t112) t79, -g(1) * (t93 * pkin(2) + t92 * qJ(3) + t130) - g(2) * (t91 * pkin(2) + t90 * qJ(3) + t127) - g(3) * ((pkin(2) * t117 - qJ(3) * t120) * t111 + t142) 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t87, t123, t79, -g(1) * t125 - g(2) * t122 - g(3) * t126, 0, 0, 0, 0, 0, 0, t77, t76, -t123, -g(1) * (t83 * pkin(4) + t82 * pkin(10) + t125) - g(2) * (t81 * pkin(4) + t80 * pkin(10) + t122) - g(3) * (t87 * pkin(4) + t86 * pkin(10) + t126) 0, 0, 0, 0, 0, 0, t77, t76, -t123, -g(1) * (pkin(5) * t140 + t83 * t103 - t82 * t114 + t125) - g(2) * (pkin(5) * t141 + t81 * t103 - t80 * t114 + t122) - g(3) * (-pkin(5) * t129 + t87 * t103 - t86 * t114 + t126);];
U_reg  = t1;
