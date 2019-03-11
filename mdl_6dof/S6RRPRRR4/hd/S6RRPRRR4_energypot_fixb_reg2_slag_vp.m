% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR4
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:34:48
% EndTime: 2019-03-09 13:34:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (359->99), mult. (712->153), div. (0->0), fcn. (880->14), ass. (0->57)
t118 = sin(pkin(12));
t120 = cos(pkin(12));
t124 = sin(qJ(2));
t128 = cos(qJ(2));
t102 = -t124 * t118 + t128 * t120;
t121 = cos(pkin(6));
t153 = t121 * pkin(8) + pkin(7);
t119 = sin(pkin(6));
t100 = t121 * t124 * pkin(2) + (-pkin(8) - qJ(3)) * t119;
t112 = t128 * pkin(2) + pkin(1);
t125 = sin(qJ(1));
t129 = cos(qJ(1));
t152 = t129 * t100 + t125 * t112;
t151 = t119 * t124;
t150 = t119 * t125;
t149 = t119 * t129;
t123 = sin(qJ(4));
t148 = t121 * t123;
t146 = t125 * t124;
t145 = t125 * t128;
t143 = t129 * t124;
t142 = t129 * t128;
t141 = t123 * t150;
t140 = t123 * t149;
t139 = -t125 * t100 + t129 * t112;
t138 = pkin(2) * t151 + t121 * qJ(3) + t153;
t137 = g(1) * t125 - g(2) * t129;
t136 = t128 * t118 + t124 * t120;
t127 = cos(qJ(4));
t111 = t127 * pkin(4) + pkin(3);
t130 = -pkin(10) - pkin(9);
t132 = t102 * t121;
t88 = -t125 * t132 - t129 * t136;
t99 = t136 * t121;
t89 = t129 * t102 - t125 * t99;
t135 = pkin(4) * t141 + t89 * t111 + t88 * t130 + t139;
t97 = t102 * t119;
t98 = t136 * t119;
t134 = pkin(4) * t148 + t98 * t111 + t97 * t130 + t138;
t117 = qJ(4) + qJ(5);
t114 = sin(t117);
t115 = cos(t117);
t87 = t125 * t102 + t129 * t99;
t78 = t87 * t114 + t115 * t149;
t80 = t89 * t114 - t115 * t150;
t90 = t98 * t114 - t121 * t115;
t133 = g(1) * t80 + g(2) * t78 + g(3) * t90;
t86 = -t125 * t136 + t129 * t132;
t77 = g(1) * t88 + g(2) * t86 + g(3) * t97;
t131 = -pkin(4) * t140 + t87 * t111 + t86 * t130 + t152;
t126 = cos(qJ(6));
t122 = sin(qJ(6));
t96 = -g(3) * t121 - t119 * t137;
t91 = t121 * t114 + t98 * t115;
t81 = t114 * t150 + t89 * t115;
t79 = -t114 * t149 + t87 * t115;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t129 - g(2) * t125, t137, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t121 * t146 + t142) - g(2) * (t121 * t143 + t145) - g(3) * t151, -g(1) * (-t121 * t145 - t143) - g(2) * (t121 * t142 - t146) - g(3) * t119 * t128, t96, -g(1) * (t129 * pkin(1) + pkin(8) * t150) - g(2) * (t125 * pkin(1) - pkin(8) * t149) - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 - g(3) * t98, -t77, t96, -g(1) * t139 - g(2) * t152 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t127 + t141) - g(2) * (t87 * t127 - t140) - g(3) * (t98 * t127 + t148) -g(1) * (-t89 * t123 + t127 * t150) - g(2) * (-t87 * t123 - t127 * t149) - g(3) * (t121 * t127 - t98 * t123) t77, -g(1) * (t89 * pkin(3) - t88 * pkin(9) + t139) - g(2) * (t87 * pkin(3) - t86 * pkin(9) + t152) - g(3) * (t98 * pkin(3) - t97 * pkin(9) + t138) 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t91, t133, t77, -g(1) * t135 - g(2) * t131 - g(3) * t134, 0, 0, 0, 0, 0, 0, -g(1) * (-t88 * t122 + t81 * t126) - g(2) * (-t86 * t122 + t79 * t126) - g(3) * (-t97 * t122 + t91 * t126) -g(1) * (-t81 * t122 - t88 * t126) - g(2) * (-t79 * t122 - t86 * t126) - g(3) * (-t91 * t122 - t97 * t126) -t133, -g(1) * (t81 * pkin(5) + t80 * pkin(11) + t135) - g(2) * (t79 * pkin(5) + t78 * pkin(11) + t131) - g(3) * (t91 * pkin(5) + t90 * pkin(11) + t134);];
U_reg  = t1;
