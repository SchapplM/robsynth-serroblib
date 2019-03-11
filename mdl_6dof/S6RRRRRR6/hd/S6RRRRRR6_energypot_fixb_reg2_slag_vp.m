% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR6
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:21:34
% EndTime: 2019-03-10 04:21:34
% DurationCPUTime: 0.27s
% Computational Cost: add. (331->106), mult. (533->152), div. (0->0), fcn. (633->14), ass. (0->55)
t116 = cos(pkin(6));
t148 = t116 * pkin(8) + pkin(7);
t117 = sin(qJ(5));
t123 = cos(qJ(2));
t124 = cos(qJ(1));
t136 = t124 * t123;
t119 = sin(qJ(2));
t120 = sin(qJ(1));
t139 = t120 * t119;
t92 = -t116 * t136 + t139;
t147 = t92 * t117;
t137 = t124 * t119;
t138 = t120 * t123;
t94 = t116 * t138 + t137;
t146 = t94 * t117;
t115 = sin(pkin(6));
t145 = t115 * t119;
t144 = t115 * t120;
t122 = cos(qJ(3));
t143 = t115 * t122;
t142 = t115 * t123;
t141 = t115 * t124;
t118 = sin(qJ(3));
t140 = t116 * t118;
t135 = t124 * pkin(1) + pkin(8) * t144;
t134 = t117 * t142;
t133 = t118 * t144;
t111 = t120 * pkin(1);
t132 = -pkin(8) * t141 + t111;
t105 = t122 * pkin(3) + pkin(2);
t126 = -pkin(10) - pkin(9);
t95 = -t116 * t139 + t136;
t131 = pkin(3) * t133 + t95 * t105 - t94 * t126 + t135;
t130 = pkin(3) * t140 + t105 * t145 + t126 * t142 + t148;
t129 = g(1) * t120 - g(2) * t124;
t114 = qJ(3) + qJ(4);
t107 = sin(t114);
t109 = cos(t114);
t93 = t116 * t137 + t138;
t82 = t93 * t107 + t109 * t141;
t84 = t95 * t107 - t109 * t144;
t88 = t107 * t145 - t116 * t109;
t128 = g(1) * t84 + g(2) * t82 + g(3) * t88;
t81 = -g(1) * t94 - g(2) * t92 + g(3) * t142;
t127 = t111 + t93 * t105 - t92 * t126 + (-pkin(3) * t118 - pkin(8)) * t141;
t125 = -pkin(12) - pkin(11);
t121 = cos(qJ(5));
t113 = qJ(5) + qJ(6);
t108 = cos(t113);
t106 = sin(t113);
t104 = t121 * pkin(5) + pkin(4);
t89 = t116 * t107 + t109 * t145;
t85 = t107 * t144 + t95 * t109;
t83 = -t107 * t141 + t93 * t109;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t124 - g(2) * t120, t129, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t145, -t81, -g(3) * t116 - t115 * t129, -g(1) * t135 - g(2) * t132 - g(3) * t148, 0, 0, 0, 0, 0, 0, -g(1) * (t95 * t122 + t133) - g(2) * (-t118 * t141 + t93 * t122) - g(3) * (t119 * t143 + t140) -g(1) * (-t95 * t118 + t120 * t143) - g(2) * (-t93 * t118 - t122 * t141) - g(3) * (t116 * t122 - t118 * t145) t81, -g(1) * (t95 * pkin(2) + t94 * pkin(9) + t135) - g(2) * (t93 * pkin(2) + t92 * pkin(9) + t132) - g(3) * ((pkin(2) * t119 - pkin(9) * t123) * t115 + t148) 0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t83 - g(3) * t89, t128, t81, -g(1) * t131 - g(2) * t127 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * (t85 * t121 + t146) - g(2) * (t83 * t121 + t147) - g(3) * (t89 * t121 - t134) -g(1) * (-t85 * t117 + t94 * t121) - g(2) * (-t83 * t117 + t92 * t121) - g(3) * (-t89 * t117 - t121 * t142) -t128, -g(1) * (t85 * pkin(4) + t84 * pkin(11) + t131) - g(2) * (t83 * pkin(4) + t82 * pkin(11) + t127) - g(3) * (t89 * pkin(4) + t88 * pkin(11) + t130) 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t106 + t85 * t108) - g(2) * (t92 * t106 + t83 * t108) - g(3) * (-t106 * t142 + t89 * t108) -g(1) * (-t85 * t106 + t94 * t108) - g(2) * (-t83 * t106 + t92 * t108) - g(3) * (-t89 * t106 - t108 * t142) -t128, -g(1) * (pkin(5) * t146 + t85 * t104 - t84 * t125 + t131) - g(2) * (pkin(5) * t147 + t83 * t104 - t82 * t125 + t127) - g(3) * (-pkin(5) * t134 + t89 * t104 - t88 * t125 + t130);];
U_reg  = t1;
