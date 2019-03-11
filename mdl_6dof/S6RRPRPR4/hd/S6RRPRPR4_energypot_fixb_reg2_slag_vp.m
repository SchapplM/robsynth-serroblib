% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:26:16
% EndTime: 2019-03-09 10:26:16
% DurationCPUTime: 0.28s
% Computational Cost: add. (359->99), mult. (712->153), div. (0->0), fcn. (880->14), ass. (0->57)
t118 = sin(pkin(11));
t120 = cos(pkin(11));
t125 = sin(qJ(2));
t129 = cos(qJ(2));
t102 = -t125 * t118 + t129 * t120;
t121 = cos(pkin(6));
t153 = t121 * pkin(8) + pkin(7);
t119 = sin(pkin(6));
t100 = t121 * t125 * pkin(2) + (-pkin(8) - qJ(3)) * t119;
t112 = t129 * pkin(2) + pkin(1);
t126 = sin(qJ(1));
t130 = cos(qJ(1));
t152 = t130 * t100 + t126 * t112;
t151 = t119 * t125;
t150 = t119 * t126;
t149 = t119 * t130;
t124 = sin(qJ(4));
t148 = t121 * t124;
t146 = t126 * t125;
t145 = t126 * t129;
t143 = t130 * t125;
t142 = t130 * t129;
t141 = t124 * t150;
t140 = t124 * t149;
t139 = -t126 * t100 + t130 * t112;
t138 = pkin(2) * t151 + t121 * qJ(3) + t153;
t137 = g(1) * t126 - g(2) * t130;
t136 = t129 * t118 + t125 * t120;
t128 = cos(qJ(4));
t111 = t128 * pkin(4) + pkin(3);
t122 = -qJ(5) - pkin(9);
t132 = t102 * t121;
t88 = -t126 * t132 - t130 * t136;
t99 = t136 * t121;
t89 = t130 * t102 - t126 * t99;
t135 = pkin(4) * t141 + t89 * t111 + t88 * t122 + t139;
t97 = t102 * t119;
t98 = t136 * t119;
t134 = pkin(4) * t148 + t98 * t111 + t97 * t122 + t138;
t117 = qJ(4) + pkin(12);
t113 = sin(t117);
t114 = cos(t117);
t87 = t126 * t102 + t130 * t99;
t78 = t87 * t113 + t114 * t149;
t80 = t89 * t113 - t114 * t150;
t90 = t98 * t113 - t121 * t114;
t133 = g(1) * t80 + g(2) * t78 + g(3) * t90;
t86 = -t126 * t136 + t130 * t132;
t77 = g(1) * t88 + g(2) * t86 + g(3) * t97;
t131 = -pkin(4) * t140 + t87 * t111 + t86 * t122 + t152;
t127 = cos(qJ(6));
t123 = sin(qJ(6));
t96 = -g(3) * t121 - t119 * t137;
t91 = t121 * t113 + t98 * t114;
t81 = t113 * t150 + t89 * t114;
t79 = -t113 * t149 + t87 * t114;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t130 - g(2) * t126, t137, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t121 * t146 + t142) - g(2) * (t121 * t143 + t145) - g(3) * t151, -g(1) * (-t121 * t145 - t143) - g(2) * (t121 * t142 - t146) - g(3) * t119 * t129, t96, -g(1) * (t130 * pkin(1) + pkin(8) * t150) - g(2) * (t126 * pkin(1) - pkin(8) * t149) - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 - g(3) * t98, -t77, t96, -g(1) * t139 - g(2) * t152 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t128 + t141) - g(2) * (t87 * t128 - t140) - g(3) * (t98 * t128 + t148) -g(1) * (-t89 * t124 + t128 * t150) - g(2) * (-t87 * t124 - t128 * t149) - g(3) * (t121 * t128 - t98 * t124) t77, -g(1) * (t89 * pkin(3) - t88 * pkin(9) + t139) - g(2) * (t87 * pkin(3) - t86 * pkin(9) + t152) - g(3) * (t98 * pkin(3) - t97 * pkin(9) + t138) 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t91, t133, t77, -g(1) * t135 - g(2) * t131 - g(3) * t134, 0, 0, 0, 0, 0, 0, -g(1) * (-t88 * t123 + t81 * t127) - g(2) * (-t86 * t123 + t79 * t127) - g(3) * (-t97 * t123 + t91 * t127) -g(1) * (-t81 * t123 - t88 * t127) - g(2) * (-t79 * t123 - t86 * t127) - g(3) * (-t91 * t123 - t97 * t127) -t133, -g(1) * (t81 * pkin(5) + t80 * pkin(10) + t135) - g(2) * (t79 * pkin(5) + t78 * pkin(10) + t131) - g(3) * (t91 * pkin(5) + t90 * pkin(10) + t134);];
U_reg  = t1;
