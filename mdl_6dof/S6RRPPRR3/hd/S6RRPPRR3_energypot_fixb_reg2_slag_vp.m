% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:36
% EndTime: 2019-03-09 08:58:36
% DurationCPUTime: 0.27s
% Computational Cost: add. (359->99), mult. (712->153), div. (0->0), fcn. (880->14), ass. (0->57)
t119 = sin(pkin(11));
t122 = cos(pkin(11));
t126 = sin(qJ(2));
t129 = cos(qJ(2));
t102 = -t126 * t119 + t129 * t122;
t123 = cos(pkin(6));
t153 = t123 * pkin(8) + pkin(7);
t120 = sin(pkin(6));
t100 = t123 * t126 * pkin(2) + (-pkin(8) - qJ(3)) * t120;
t112 = t129 * pkin(2) + pkin(1);
t127 = sin(qJ(1));
t130 = cos(qJ(1));
t152 = t130 * t100 + t127 * t112;
t151 = t120 * t126;
t150 = t120 * t127;
t149 = t120 * t130;
t118 = sin(pkin(12));
t148 = t123 * t118;
t146 = t127 * t126;
t145 = t127 * t129;
t143 = t130 * t126;
t142 = t130 * t129;
t141 = t118 * t150;
t140 = t118 * t149;
t139 = -t127 * t100 + t130 * t112;
t138 = pkin(2) * t151 + t123 * qJ(3) + t153;
t137 = g(1) * t127 - g(2) * t130;
t136 = t129 * t119 + t126 * t122;
t121 = cos(pkin(12));
t111 = t121 * pkin(4) + pkin(3);
t124 = -pkin(9) - qJ(4);
t132 = t102 * t123;
t88 = -t127 * t132 - t130 * t136;
t99 = t136 * t123;
t89 = t130 * t102 - t127 * t99;
t135 = pkin(4) * t141 + t89 * t111 + t88 * t124 + t139;
t97 = t102 * t120;
t98 = t136 * t120;
t134 = pkin(4) * t148 + t98 * t111 + t97 * t124 + t138;
t117 = pkin(12) + qJ(5);
t113 = sin(t117);
t114 = cos(t117);
t87 = t127 * t102 + t130 * t99;
t78 = t87 * t113 + t114 * t149;
t80 = t89 * t113 - t114 * t150;
t90 = t98 * t113 - t123 * t114;
t133 = g(1) * t80 + g(2) * t78 + g(3) * t90;
t86 = -t127 * t136 + t130 * t132;
t77 = g(1) * t88 + g(2) * t86 + g(3) * t97;
t131 = -pkin(4) * t140 + t87 * t111 + t86 * t124 + t152;
t128 = cos(qJ(6));
t125 = sin(qJ(6));
t96 = -g(3) * t123 - t120 * t137;
t91 = t123 * t113 + t98 * t114;
t81 = t113 * t150 + t89 * t114;
t79 = -t113 * t149 + t87 * t114;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t130 - g(2) * t127, t137, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * (-t123 * t146 + t142) - g(2) * (t123 * t143 + t145) - g(3) * t151, -g(1) * (-t123 * t145 - t143) - g(2) * (t123 * t142 - t146) - g(3) * t120 * t129, t96, -g(1) * (t130 * pkin(1) + pkin(8) * t150) - g(2) * (t127 * pkin(1) - pkin(8) * t149) - g(3) * t153, 0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87 - g(3) * t98, -t77, t96, -g(1) * t139 - g(2) * t152 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t121 + t141) - g(2) * (t87 * t121 - t140) - g(3) * (t98 * t121 + t148) -g(1) * (-t89 * t118 + t121 * t150) - g(2) * (-t87 * t118 - t121 * t149) - g(3) * (-t98 * t118 + t123 * t121) t77, -g(1) * (t89 * pkin(3) - t88 * qJ(4) + t139) - g(2) * (t87 * pkin(3) - t86 * qJ(4) + t152) - g(3) * (t98 * pkin(3) - t97 * qJ(4) + t138) 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t91, t133, t77, -g(1) * t135 - g(2) * t131 - g(3) * t134, 0, 0, 0, 0, 0, 0, -g(1) * (-t88 * t125 + t81 * t128) - g(2) * (-t86 * t125 + t79 * t128) - g(3) * (-t97 * t125 + t91 * t128) -g(1) * (-t81 * t125 - t88 * t128) - g(2) * (-t79 * t125 - t86 * t128) - g(3) * (-t91 * t125 - t97 * t128) -t133, -g(1) * (t81 * pkin(5) + t80 * pkin(10) + t135) - g(2) * (t79 * pkin(5) + t78 * pkin(10) + t131) - g(3) * (t91 * pkin(5) + t90 * pkin(10) + t134);];
U_reg  = t1;
