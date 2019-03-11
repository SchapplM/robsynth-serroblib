% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:52:20
% EndTime: 2019-03-09 18:52:21
% DurationCPUTime: 0.27s
% Computational Cost: add. (331->106), mult. (533->152), div. (0->0), fcn. (633->14), ass. (0->55)
t115 = cos(pkin(6));
t147 = t115 * pkin(8) + pkin(7);
t117 = sin(qJ(5));
t123 = cos(qJ(2));
t124 = cos(qJ(1));
t135 = t124 * t123;
t119 = sin(qJ(2));
t120 = sin(qJ(1));
t138 = t120 * t119;
t91 = -t115 * t135 + t138;
t146 = t91 * t117;
t136 = t124 * t119;
t137 = t120 * t123;
t93 = t115 * t137 + t136;
t145 = t93 * t117;
t114 = sin(pkin(6));
t144 = t114 * t119;
t143 = t114 * t120;
t122 = cos(qJ(3));
t142 = t114 * t122;
t141 = t114 * t123;
t140 = t114 * t124;
t118 = sin(qJ(3));
t139 = t115 * t118;
t134 = t124 * pkin(1) + pkin(8) * t143;
t133 = t117 * t141;
t132 = t118 * t143;
t110 = t120 * pkin(1);
t131 = -pkin(8) * t140 + t110;
t104 = t122 * pkin(3) + pkin(2);
t116 = -qJ(4) - pkin(9);
t130 = pkin(3) * t139 + t104 * t144 + t116 * t141 + t147;
t94 = -t115 * t138 + t135;
t129 = pkin(3) * t132 + t94 * t104 - t93 * t116 + t134;
t128 = g(1) * t120 - g(2) * t124;
t112 = qJ(3) + pkin(12);
t105 = sin(t112);
t106 = cos(t112);
t92 = t115 * t136 + t137;
t81 = t92 * t105 + t106 * t140;
t83 = t94 * t105 - t106 * t143;
t87 = t105 * t144 - t115 * t106;
t127 = g(1) * t83 + g(2) * t81 + g(3) * t87;
t80 = -g(1) * t93 - g(2) * t91 + g(3) * t141;
t126 = t110 + t92 * t104 - t91 * t116 + (-pkin(3) * t118 - pkin(8)) * t140;
t125 = -pkin(11) - pkin(10);
t121 = cos(qJ(5));
t113 = qJ(5) + qJ(6);
t108 = cos(t113);
t107 = sin(t113);
t103 = t121 * pkin(5) + pkin(4);
t88 = t115 * t105 + t106 * t144;
t84 = t105 * t143 + t94 * t106;
t82 = -t105 * t140 + t92 * t106;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t124 - g(2) * t120, t128, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t144, -t80, -g(3) * t115 - t114 * t128, -g(1) * t134 - g(2) * t131 - g(3) * t147, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t122 + t132) - g(2) * (-t118 * t140 + t92 * t122) - g(3) * (t119 * t142 + t139) -g(1) * (-t94 * t118 + t120 * t142) - g(2) * (-t92 * t118 - t122 * t140) - g(3) * (t115 * t122 - t118 * t144) t80, -g(1) * (t94 * pkin(2) + t93 * pkin(9) + t134) - g(2) * (t92 * pkin(2) + t91 * pkin(9) + t131) - g(3) * ((pkin(2) * t119 - pkin(9) * t123) * t114 + t147) 0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82 - g(3) * t88, t127, t80, -g(1) * t129 - g(2) * t126 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t121 + t145) - g(2) * (t82 * t121 + t146) - g(3) * (t88 * t121 - t133) -g(1) * (-t84 * t117 + t93 * t121) - g(2) * (-t82 * t117 + t91 * t121) - g(3) * (-t88 * t117 - t121 * t141) -t127, -g(1) * (t84 * pkin(4) + t83 * pkin(10) + t129) - g(2) * (t82 * pkin(4) + t81 * pkin(10) + t126) - g(3) * (t88 * pkin(4) + t87 * pkin(10) + t130) 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t107 + t84 * t108) - g(2) * (t91 * t107 + t82 * t108) - g(3) * (-t107 * t141 + t88 * t108) -g(1) * (-t84 * t107 + t93 * t108) - g(2) * (-t82 * t107 + t91 * t108) - g(3) * (-t88 * t107 - t108 * t141) -t127, -g(1) * (pkin(5) * t145 + t84 * t103 - t83 * t125 + t129) - g(2) * (pkin(5) * t146 + t82 * t103 - t81 * t125 + t126) - g(3) * (-pkin(5) * t133 + t88 * t103 - t87 * t125 + t130);];
U_reg  = t1;
