% Calculate inertial parameters regressor of potential energy for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:19
% EndTime: 2019-03-08 21:38:20
% DurationCPUTime: 0.23s
% Computational Cost: add. (322->92), mult. (648->131), div. (0->0), fcn. (793->12), ass. (0->56)
t108 = sin(pkin(6));
t139 = pkin(7) * t108;
t106 = sin(pkin(11));
t107 = sin(pkin(10));
t110 = cos(pkin(10));
t114 = sin(qJ(2));
t111 = cos(pkin(6));
t116 = cos(qJ(2));
t130 = t111 * t116;
t88 = t107 * t114 - t110 * t130;
t138 = t88 * t106;
t90 = t107 * t130 + t110 * t114;
t137 = t90 * t106;
t136 = t110 * pkin(1) + t107 * t139;
t113 = sin(qJ(3));
t135 = t108 * t113;
t134 = t108 * t114;
t115 = cos(qJ(3));
t133 = t108 * t115;
t132 = t108 * t116;
t131 = t111 * t114;
t129 = t111 * pkin(7) + qJ(1);
t128 = pkin(2) * t134 + t129;
t127 = t107 * pkin(1) - t110 * t139;
t126 = g(1) * t107 - g(2) * t110;
t91 = -t107 * t131 + t110 * t116;
t125 = t91 * pkin(2) + t90 * pkin(8) + t136;
t124 = -pkin(8) * t132 + t128;
t105 = pkin(11) + qJ(5);
t100 = sin(t105);
t101 = cos(t105);
t89 = t107 * t116 + t110 * t131;
t77 = -t110 * t135 + t89 * t115;
t66 = t77 * t100 - t88 * t101;
t79 = t107 * t135 + t91 * t115;
t68 = t79 * t100 - t90 * t101;
t93 = t111 * t113 + t114 * t133;
t72 = t93 * t100 + t101 * t132;
t123 = g(1) * t68 + g(2) * t66 + g(3) * t72;
t76 = t110 * t133 + t89 * t113;
t78 = -t107 * t133 + t91 * t113;
t92 = -t111 * t115 + t113 * t134;
t122 = g(1) * t78 + g(2) * t76 + g(3) * t92;
t121 = t89 * pkin(2) + t88 * pkin(8) + t127;
t112 = -pkin(9) - qJ(4);
t109 = cos(pkin(11));
t99 = t109 * pkin(4) + pkin(3);
t120 = pkin(4) * t137 - t78 * t112 + t79 * t99 + t125;
t119 = -g(1) * t90 - g(2) * t88 + g(3) * t132;
t118 = pkin(4) * t138 - t76 * t112 + t77 * t99 + t121;
t117 = t93 * t99 - t92 * t112 + (-pkin(4) * t106 - pkin(8)) * t132 + t128;
t73 = -t100 * t132 + t93 * t101;
t69 = t90 * t100 + t79 * t101;
t67 = t88 * t100 + t77 * t101;
t64 = -g(1) * t69 - g(2) * t67 - g(3) * t73;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t110 - g(2) * t107, t126, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89 - g(3) * t134, -t119, -g(3) * t111 - t126 * t108, -g(1) * t136 - g(2) * t127 - g(3) * t129, 0, 0, 0, 0, 0, 0, -g(1) * t79 - g(2) * t77 - g(3) * t93, t122, t119, -g(1) * t125 - g(2) * t121 - g(3) * t124, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t109 + t137) - g(2) * (t77 * t109 + t138) - g(3) * (-t106 * t132 + t93 * t109) -g(1) * (-t79 * t106 + t90 * t109) - g(2) * (-t77 * t106 + t88 * t109) - g(3) * (-t93 * t106 - t109 * t132) -t122, -g(1) * (t79 * pkin(3) + t78 * qJ(4) + t125) - g(2) * (t77 * pkin(3) + t76 * qJ(4) + t121) - g(3) * (t93 * pkin(3) + t92 * qJ(4) + t124) 0, 0, 0, 0, 0, 0, t64, t123, -t122, -g(1) * t120 - g(2) * t118 - g(3) * t117, 0, 0, 0, 0, 0, 0, t64, -t122, -t123, -g(1) * (t69 * pkin(5) + t68 * qJ(6) + t120) - g(2) * (t67 * pkin(5) + t66 * qJ(6) + t118) - g(3) * (t73 * pkin(5) + t72 * qJ(6) + t117);];
U_reg  = t1;
