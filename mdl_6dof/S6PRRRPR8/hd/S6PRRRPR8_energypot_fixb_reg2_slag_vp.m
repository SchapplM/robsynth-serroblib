% Calculate inertial parameters regressor of potential energy for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:53:19
% EndTime: 2019-03-08 23:53:20
% DurationCPUTime: 0.33s
% Computational Cost: add. (508->100), mult. (1330->149), div. (0->0), fcn. (1691->14), ass. (0->64)
t119 = sin(pkin(7));
t122 = cos(pkin(7));
t123 = cos(pkin(6));
t120 = sin(pkin(6));
t129 = cos(qJ(2));
t157 = t120 * t129;
t142 = t119 * t157 - t123 * t122;
t118 = sin(pkin(12));
t121 = cos(pkin(12));
t127 = sin(qJ(2));
t154 = t123 * t129;
t106 = -t118 * t154 - t121 * t127;
t159 = t120 * t122;
t143 = -t106 * t119 + t118 * t159;
t168 = pkin(5) + pkin(10);
t104 = -t118 * t127 + t121 * t154;
t155 = t123 * t127;
t105 = t118 * t129 + t121 * t155;
t126 = sin(qJ(3));
t164 = cos(qJ(3));
t149 = t119 * t164;
t147 = t120 * t149;
t148 = t122 * t164;
t88 = -t104 * t148 + t105 * t126 + t121 * t147;
t167 = t88 * pkin(10);
t107 = -t118 * t155 + t121 * t129;
t90 = -t106 * t148 + t107 * t126 - t118 * t147;
t166 = t90 * pkin(10);
t158 = t120 * t127;
t97 = -t123 * t149 + t126 * t158 - t148 * t157;
t165 = t97 * pkin(10);
t163 = cos(qJ(4));
t161 = t118 * t120;
t160 = t120 * t121;
t153 = t121 * pkin(1) + pkin(8) * t161;
t152 = t123 * pkin(8) + qJ(1);
t146 = t118 * pkin(1) - pkin(8) * t160;
t145 = g(1) * t118 - g(2) * t121;
t144 = t104 * t119 + t121 * t159;
t125 = sin(qJ(4));
t89 = t105 * t164 + (t104 * t122 - t119 * t160) * t126;
t81 = t89 * t125 + t144 * t163;
t91 = t107 * t164 + (t106 * t122 + t119 * t161) * t126;
t83 = t91 * t125 - t143 * t163;
t98 = t123 * t119 * t126 + (t122 * t126 * t129 + t127 * t164) * t120;
t92 = t98 * t125 + t142 * t163;
t141 = g(1) * t83 + g(2) * t81 + g(3) * t92;
t82 = -t125 * t144 + t163 * t89;
t84 = t125 * t143 + t163 * t91;
t93 = -t125 * t142 + t163 * t98;
t140 = g(1) * t84 + g(2) * t82 + g(3) * t93;
t139 = g(1) * t90 + g(2) * t88 + g(3) * t97;
t138 = t107 * pkin(2) + t143 * pkin(9) + t153;
t137 = t91 * pkin(3) + t138;
t136 = pkin(2) * t158 - t142 * pkin(9) + t152;
t135 = t98 * pkin(3) + t136;
t134 = t84 * pkin(4) + t83 * qJ(5) + t137;
t133 = t93 * pkin(4) + t92 * qJ(5) + t135;
t132 = t105 * pkin(2) - pkin(9) * t144 + t146;
t131 = t89 * pkin(3) + t132;
t130 = t82 * pkin(4) + t81 * qJ(5) + t131;
t128 = cos(qJ(6));
t124 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t118, t145, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t107 - g(2) * t105 - g(3) * t158, -g(1) * t106 - g(2) * t104 - g(3) * t157, -g(3) * t123 - t120 * t145, -g(1) * t153 - g(2) * t146 - g(3) * t152, 0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89 - g(3) * t98, t139, -g(1) * t143 + g(2) * t144 + g(3) * t142, -g(1) * t138 - g(2) * t132 - g(3) * t136, 0, 0, 0, 0, 0, 0, -t140, t141, -t139, -g(1) * (t137 + t166) - g(2) * (t131 + t167) - g(3) * (t135 + t165) 0, 0, 0, 0, 0, 0, -t139, t140, -t141, -g(1) * (t134 + t166) - g(2) * (t130 + t167) - g(3) * (t133 + t165) 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t124 + t90 * t128) - g(2) * (t81 * t124 + t88 * t128) - g(3) * (t92 * t124 + t97 * t128) -g(1) * (-t90 * t124 + t83 * t128) - g(2) * (-t88 * t124 + t81 * t128) - g(3) * (-t97 * t124 + t92 * t128) -t140, -g(1) * (t84 * pkin(11) + t168 * t90 + t134) - g(2) * (t82 * pkin(11) + t168 * t88 + t130) - g(3) * (t93 * pkin(11) + t168 * t97 + t133);];
U_reg  = t1;
