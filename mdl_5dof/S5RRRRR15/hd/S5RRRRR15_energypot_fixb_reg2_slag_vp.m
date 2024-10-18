% Calculate inertial parameters regressor of potential energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR15_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:17
% EndTime: 2024-09-27 22:26:17
% DurationCPUTime: 0.06s
% Computational Cost: add. (296->114), mult. (349->180), div. (0->0), fcn. (365->26), ass. (0->81)
t124 = pkin(8) + pkin(9);
t118 = sin(qJ(1));
t154 = -t118 / 0.2e1;
t123 = cos(qJ(1));
t153 = t123 / 0.2e1;
t110 = sin(pkin(6));
t152 = pkin(11) * t110;
t111 = sin(pkin(5));
t151 = g(3) * t111;
t122 = cos(qJ(2));
t101 = t122 * pkin(2) + pkin(1);
t112 = cos(pkin(6));
t109 = qJ(2) + qJ(3);
t106 = qJ(4) + t109;
t99 = cos(t106);
t150 = t112 * t99;
t98 = sin(t106);
t149 = t118 * t98;
t148 = t118 * t99;
t147 = t123 * t98;
t146 = t123 * t99;
t113 = cos(pkin(5));
t145 = t110 * t113;
t117 = sin(qJ(2));
t144 = t111 * t117;
t143 = t111 * t118;
t142 = t111 * t123;
t114 = sin(qJ(5));
t141 = t114 * t118;
t140 = t114 * t123;
t139 = t117 * t118;
t138 = t117 * t123;
t119 = cos(qJ(5));
t137 = t118 * t119;
t136 = t118 * t122;
t135 = t119 * t123;
t134 = t122 * t123;
t108 = pkin(10) + t124;
t133 = t113 * t141;
t132 = t113 * t140;
t131 = t113 * t137;
t130 = t113 * t135;
t129 = t110 * t143;
t128 = t110 * t142;
t103 = pkin(5) - t109;
t102 = pkin(5) + t109;
t127 = g(1) * t118 - g(2) * t123;
t116 = sin(qJ(3));
t121 = cos(qJ(3));
t115 = sin(qJ(4));
t120 = cos(qJ(4));
t82 = pkin(4) * t120 + t115 * t152 + pkin(3);
t85 = -pkin(4) * t115 + t120 * t152;
t75 = t116 * t85 + t121 * t82 + pkin(2);
t76 = -t116 * t82 + t121 * t85;
t126 = t117 * t75 - t122 * t76;
t125 = pkin(3) * t116 * t122 + (pkin(3) * t121 + pkin(2)) * t117;
t78 = -g(3) * t113 - t111 * t127;
t105 = cos(t109);
t104 = sin(t109);
t97 = -qJ(4) + t103;
t96 = qJ(4) + t102;
t95 = cos(t103);
t94 = cos(t102);
t93 = sin(t103);
t92 = sin(t102);
t91 = cos(t96);
t90 = sin(t97);
t89 = pkin(11) * t112 + t108;
t88 = cos(t97) / 0.2e1;
t87 = sin(t96) / 0.2e1;
t86 = pkin(3) * t105 + t101;
t84 = t95 + t94;
t83 = t92 - t93;
t81 = pkin(2) * t113 * t117 - t111 * t124;
t80 = t88 + t91 / 0.2e1;
t79 = t87 - t90 / 0.2e1;
t77 = -t108 * t111 + t113 * t125;
t74 = t117 * t76 + t122 * t75 + pkin(1);
t73 = t111 * t89 - t113 * t126;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t123 - g(2) * t118, t127, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t113 * t139 + t134) - g(2) * (t113 * t138 + t136) - g(3) * t144, -g(1) * (-t113 * t136 - t138) - g(2) * (t113 * t134 - t139) - t122 * t151, t78, -g(1) * (pkin(1) * t123 + pkin(8) * t143) - g(2) * (pkin(1) * t118 - pkin(8) * t142) - g(3) * (pkin(8) * t113 + pkin(7)), 0, 0, 0, 0, 0, 0, -g(1) * (t123 * t105 + t154 * t83) - g(2) * (t118 * t105 + t153 * t83) - g(3) * (t95 / 0.2e1 - t94 / 0.2e1), -g(1) * (-t123 * t104 + t154 * t84) - g(2) * (-t118 * t104 + t153 * t84) - g(3) * (t92 / 0.2e1 + t93 / 0.2e1), t78, -g(1) * (t101 * t123 - t118 * t81) - g(2) * (t101 * t118 + t123 * t81) - g(3) * (pkin(2) * t144 + t113 * t124 + pkin(7)), 0, 0, 0, 0, 0, 0, -g(1) * (-t118 * t79 + t146) - g(2) * (t123 * t79 + t148) - g(3) * (t88 - t91 / 0.2e1), -g(1) * (-t118 * t80 - t147) - g(2) * (t123 * t80 - t149) - g(3) * (t87 + t90 / 0.2e1), t78, -g(1) * (-t118 * t77 + t123 * t86) - g(2) * (t118 * t86 + t123 * t77) - g(3) * (t108 * t113 + t111 * t125 + pkin(7)), 0, 0, 0, 0, 0, 0, -g(1) * ((-t112 * t133 + t135) * t99 + (-t112 * t140 - t131) * t98 + t114 * t129) - g(2) * ((t112 * t132 + t137) * t99 + (-t112 * t141 + t130) * t98 - t114 * t128) - g(3) * (t114 * t145 + (t114 * t150 + t119 * t98) * t111), -g(1) * ((-t112 * t131 - t140) * t99 + (-t112 * t135 + t133) * t98 + t119 * t129) - g(2) * ((t112 * t130 - t141) * t99 + (-t112 * t137 - t132) * t98 - t119 * t128) - g(3) * (t119 * t145 + (-t114 * t98 + t119 * t150) * t111), t78 * t112 + (-g(1) * (t113 * t148 + t147) - g(2) * (-t113 * t146 + t149) + t99 * t151) * t110, -g(1) * (t118 * t73 + t123 * t74) - g(2) * (t118 * t74 - t123 * t73) - g(3) * (t111 * t126 + t113 * t89 + pkin(7));];
U_reg = t1;
