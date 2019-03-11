% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:44:12
% EndTime: 2019-03-09 15:44:12
% DurationCPUTime: 0.27s
% Computational Cost: add. (331->106), mult. (533->152), div. (0->0), fcn. (633->14), ass. (0->55)
t117 = cos(pkin(6));
t147 = t117 * pkin(8) + pkin(7);
t114 = sin(pkin(12));
t124 = cos(qJ(2));
t125 = cos(qJ(1));
t135 = t125 * t124;
t121 = sin(qJ(2));
t122 = sin(qJ(1));
t138 = t122 * t121;
t91 = -t117 * t135 + t138;
t146 = t91 * t114;
t136 = t125 * t121;
t137 = t122 * t124;
t93 = t117 * t137 + t136;
t145 = t93 * t114;
t115 = sin(pkin(6));
t144 = t115 * t121;
t143 = t115 * t122;
t123 = cos(qJ(3));
t142 = t115 * t123;
t141 = t115 * t124;
t140 = t115 * t125;
t120 = sin(qJ(3));
t139 = t117 * t120;
t134 = t125 * pkin(1) + pkin(8) * t143;
t133 = t114 * t141;
t132 = t120 * t143;
t110 = t122 * pkin(1);
t131 = -pkin(8) * t140 + t110;
t104 = t123 * pkin(3) + pkin(2);
t118 = -qJ(4) - pkin(9);
t130 = pkin(3) * t139 + t104 * t144 + t118 * t141 + t147;
t94 = -t117 * t138 + t135;
t129 = pkin(3) * t132 + t94 * t104 - t93 * t118 + t134;
t128 = g(1) * t122 - g(2) * t125;
t113 = qJ(3) + pkin(11);
t106 = sin(t113);
t108 = cos(t113);
t92 = t117 * t136 + t137;
t81 = t92 * t106 + t108 * t140;
t83 = t94 * t106 - t108 * t143;
t87 = t106 * t144 - t117 * t108;
t127 = g(1) * t83 + g(2) * t81 + g(3) * t87;
t80 = -g(1) * t93 - g(2) * t91 + g(3) * t141;
t126 = t110 + t92 * t104 - t91 * t118 + (-pkin(3) * t120 - pkin(8)) * t140;
t119 = -pkin(10) - qJ(5);
t116 = cos(pkin(12));
t112 = pkin(12) + qJ(6);
t107 = cos(t112);
t105 = sin(t112);
t103 = t116 * pkin(5) + pkin(4);
t88 = t117 * t106 + t108 * t144;
t84 = t106 * t143 + t94 * t108;
t82 = -t106 * t140 + t92 * t108;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t125 - g(2) * t122, t128, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t144, -t80, -g(3) * t117 - t115 * t128, -g(1) * t134 - g(2) * t131 - g(3) * t147, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t123 + t132) - g(2) * (-t120 * t140 + t92 * t123) - g(3) * (t121 * t142 + t139) -g(1) * (-t94 * t120 + t122 * t142) - g(2) * (-t92 * t120 - t123 * t140) - g(3) * (t117 * t123 - t120 * t144) t80, -g(1) * (t94 * pkin(2) + t93 * pkin(9) + t134) - g(2) * (t92 * pkin(2) + t91 * pkin(9) + t131) - g(3) * ((pkin(2) * t121 - pkin(9) * t124) * t115 + t147) 0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82 - g(3) * t88, t127, t80, -g(1) * t129 - g(2) * t126 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t116 + t145) - g(2) * (t82 * t116 + t146) - g(3) * (t88 * t116 - t133) -g(1) * (-t84 * t114 + t93 * t116) - g(2) * (-t82 * t114 + t91 * t116) - g(3) * (-t88 * t114 - t116 * t141) -t127, -g(1) * (t84 * pkin(4) + t83 * qJ(5) + t129) - g(2) * (t82 * pkin(4) + t81 * qJ(5) + t126) - g(3) * (t88 * pkin(4) + t87 * qJ(5) + t130) 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t105 + t84 * t107) - g(2) * (t91 * t105 + t82 * t107) - g(3) * (-t105 * t141 + t88 * t107) -g(1) * (-t84 * t105 + t93 * t107) - g(2) * (-t82 * t105 + t91 * t107) - g(3) * (-t88 * t105 - t107 * t141) -t127, -g(1) * (pkin(5) * t145 + t84 * t103 - t83 * t119 + t129) - g(2) * (pkin(5) * t146 + t82 * t103 - t81 * t119 + t126) - g(3) * (-pkin(5) * t133 + t88 * t103 - t87 * t119 + t130);];
U_reg  = t1;
