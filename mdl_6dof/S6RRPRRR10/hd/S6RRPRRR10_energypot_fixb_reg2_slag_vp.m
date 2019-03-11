% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:25:22
% EndTime: 2019-03-09 14:25:22
% DurationCPUTime: 0.27s
% Computational Cost: add. (331->106), mult. (533->151), div. (0->0), fcn. (633->14), ass. (0->54)
t117 = cos(pkin(6));
t146 = t117 * pkin(8) + pkin(7);
t119 = sin(qJ(5));
t123 = cos(qJ(2));
t124 = cos(qJ(1));
t135 = t124 * t123;
t120 = sin(qJ(2));
t121 = sin(qJ(1));
t138 = t121 * t120;
t91 = -t117 * t135 + t138;
t145 = t91 * t119;
t136 = t124 * t120;
t137 = t121 * t123;
t93 = t117 * t137 + t136;
t144 = t93 * t119;
t115 = sin(pkin(6));
t143 = t115 * t120;
t142 = t115 * t121;
t141 = t115 * t123;
t140 = t115 * t124;
t114 = sin(pkin(12));
t139 = t117 * t114;
t134 = t124 * pkin(1) + pkin(8) * t142;
t133 = t119 * t141;
t132 = t114 * t142;
t110 = t121 * pkin(1);
t131 = -pkin(8) * t140 + t110;
t116 = cos(pkin(12));
t103 = t116 * pkin(3) + pkin(2);
t118 = -pkin(9) - qJ(3);
t130 = pkin(3) * t139 + t103 * t143 + t118 * t141 + t146;
t94 = -t117 * t138 + t135;
t129 = pkin(3) * t132 + t94 * t103 - t93 * t118 + t134;
t128 = g(1) * t121 - g(2) * t124;
t112 = pkin(12) + qJ(4);
t105 = sin(t112);
t106 = cos(t112);
t92 = t117 * t136 + t137;
t81 = t92 * t105 + t106 * t140;
t83 = t94 * t105 - t106 * t142;
t87 = t105 * t143 - t117 * t106;
t127 = g(1) * t83 + g(2) * t81 + g(3) * t87;
t80 = -g(1) * t93 - g(2) * t91 + g(3) * t141;
t126 = t110 + t92 * t103 - t91 * t118 + (-pkin(3) * t114 - pkin(8)) * t140;
t125 = -pkin(11) - pkin(10);
t122 = cos(qJ(5));
t113 = qJ(5) + qJ(6);
t108 = cos(t113);
t107 = sin(t113);
t104 = t122 * pkin(5) + pkin(4);
t88 = t117 * t105 + t106 * t143;
t84 = t105 * t142 + t94 * t106;
t82 = -t105 * t140 + t92 * t106;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t124 - g(2) * t121, t128, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t143, -t80, -g(3) * t117 - t115 * t128, -g(1) * t134 - g(2) * t131 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t116 + t132) - g(2) * (-t114 * t140 + t92 * t116) - g(3) * (t116 * t143 + t139) -g(1) * (-t94 * t114 + t116 * t142) - g(2) * (-t92 * t114 - t116 * t140) - g(3) * (-t114 * t143 + t117 * t116) t80, -g(1) * (t94 * pkin(2) + t93 * qJ(3) + t134) - g(2) * (t92 * pkin(2) + t91 * qJ(3) + t131) - g(3) * ((pkin(2) * t120 - qJ(3) * t123) * t115 + t146) 0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82 - g(3) * t88, t127, t80, -g(1) * t129 - g(2) * t126 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t122 + t144) - g(2) * (t82 * t122 + t145) - g(3) * (t88 * t122 - t133) -g(1) * (-t84 * t119 + t93 * t122) - g(2) * (-t82 * t119 + t91 * t122) - g(3) * (-t88 * t119 - t122 * t141) -t127, -g(1) * (t84 * pkin(4) + t83 * pkin(10) + t129) - g(2) * (t82 * pkin(4) + t81 * pkin(10) + t126) - g(3) * (t88 * pkin(4) + t87 * pkin(10) + t130) 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t107 + t84 * t108) - g(2) * (t91 * t107 + t82 * t108) - g(3) * (-t107 * t141 + t88 * t108) -g(1) * (-t84 * t107 + t93 * t108) - g(2) * (-t82 * t107 + t91 * t108) - g(3) * (-t88 * t107 - t108 * t141) -t127, -g(1) * (pkin(5) * t144 + t84 * t104 - t83 * t125 + t129) - g(2) * (pkin(5) * t145 + t82 * t104 - t81 * t125 + t126) - g(3) * (-pkin(5) * t133 + t88 * t104 - t87 * t125 + t130);];
U_reg  = t1;
