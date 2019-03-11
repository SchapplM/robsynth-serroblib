% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:32
% EndTime: 2019-03-09 11:29:33
% DurationCPUTime: 0.24s
% Computational Cost: add. (248->97), mult. (524->133), div. (0->0), fcn. (622->12), ass. (0->56)
t108 = cos(pkin(6));
t106 = sin(pkin(6));
t111 = sin(qJ(2));
t138 = t106 * t111;
t143 = t108 * pkin(3) + pkin(9) * t138;
t105 = sin(pkin(11));
t142 = pkin(5) * t105;
t141 = t108 * pkin(8) + pkin(7);
t115 = cos(qJ(1));
t132 = t115 * t111;
t112 = sin(qJ(1));
t114 = cos(qJ(2));
t133 = t112 * t114;
t88 = t108 * t132 + t133;
t140 = t88 * t105;
t137 = t106 * t112;
t139 = t115 * pkin(1) + pkin(8) * t137;
t136 = t106 * t114;
t135 = t106 * t115;
t134 = t112 * t111;
t131 = t115 * t114;
t130 = pkin(2) * t138 + t141;
t129 = pkin(8) * t135;
t128 = (-pkin(3) - pkin(8)) * t115;
t102 = t112 * pkin(1);
t87 = -t108 * t131 + t134;
t127 = t88 * pkin(2) + t87 * qJ(3) + t102;
t126 = g(1) * t112 - g(2) * t115;
t89 = t108 * t133 + t132;
t90 = -t108 * t134 + t131;
t125 = t90 * pkin(2) + t89 * qJ(3) + t139;
t124 = pkin(3) * t137 + t125;
t123 = -qJ(3) * t136 + t130;
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t77 = t110 * t137 - t89 * t113;
t79 = t110 * t135 + t87 * t113;
t85 = t108 * t110 + t113 * t136;
t122 = g(1) * t77 - g(2) * t79 + g(3) * t85;
t121 = t88 * pkin(9) + t127;
t120 = g(1) * t90 + g(2) * t88 + g(3) * t138;
t119 = -g(1) * t89 - g(2) * t87 + g(3) * t136;
t118 = t123 + t143;
t117 = t90 * pkin(9) + t124;
t116 = t106 * t128 + t121;
t109 = -pkin(10) - qJ(5);
t107 = cos(pkin(11));
t104 = pkin(11) + qJ(6);
t99 = cos(t104);
t98 = sin(t104);
t97 = t107 * pkin(5) + pkin(4);
t86 = t108 * t113 - t110 * t136;
t81 = -g(3) * t108 - t126 * t106;
t80 = t87 * t110 - t113 * t135;
t78 = t89 * t110 + t113 * t137;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t115 - g(2) * t112, t126, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t120, -t119, t81, -g(1) * t139 - g(2) * (t102 - t129) - g(3) * t141, 0, 0, 0, 0, 0, 0, t81, t120, t119, -g(1) * t125 - g(2) * (t127 - t129) - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t80 - g(3) * t86, t122, -t120, -g(1) * t117 - g(2) * t116 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t105 + t78 * t107) - g(2) * (t80 * t107 + t140) - g(3) * (t105 * t138 + t86 * t107) -g(1) * (-t78 * t105 + t90 * t107) - g(2) * (-t80 * t105 + t88 * t107) - g(3) * (-t86 * t105 + t107 * t138) -t122, -g(1) * (t78 * pkin(4) + t77 * qJ(5) + t117) - g(2) * (t80 * pkin(4) - t79 * qJ(5) + t116) - g(3) * (t86 * pkin(4) + t85 * qJ(5) + t118) 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t99 + t90 * t98) - g(2) * (t80 * t99 + t88 * t98) - g(3) * (t98 * t138 + t86 * t99) -g(1) * (-t78 * t98 + t90 * t99) - g(2) * (-t80 * t98 + t88 * t99) - g(3) * (t99 * t138 - t86 * t98) -t122, -g(1) * (-t77 * t109 + t78 * t97 + (pkin(9) + t142) * t90 + t124) - g(2) * (pkin(5) * t140 + t79 * t109 + t80 * t97 + t121) - g(3) * (-t85 * t109 + t86 * t97 + t130 + t143) + (-g(3) * (-qJ(3) * t114 + t111 * t142) - g(2) * t128) * t106;];
U_reg  = t1;
