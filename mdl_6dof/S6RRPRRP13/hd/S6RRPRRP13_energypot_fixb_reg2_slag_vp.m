% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP13_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:31
% EndTime: 2019-03-09 13:01:31
% DurationCPUTime: 0.20s
% Computational Cost: add. (236->86), mult. (524->116), div. (0->0), fcn. (622->10), ass. (0->55)
t102 = cos(pkin(6));
t101 = sin(pkin(6));
t106 = sin(qJ(2));
t135 = t101 * t106;
t139 = pkin(3) * t102 + pkin(9) * t135;
t138 = pkin(8) * t102 + pkin(7);
t104 = sin(qJ(5));
t111 = cos(qJ(1));
t128 = t111 * t106;
t107 = sin(qJ(1));
t110 = cos(qJ(2));
t129 = t107 * t110;
t87 = t102 * t128 + t129;
t137 = t87 * t104;
t134 = t101 * t107;
t136 = pkin(1) * t111 + pkin(8) * t134;
t133 = t101 * t110;
t132 = t101 * t111;
t131 = t104 * t106;
t130 = t107 * t106;
t127 = t111 * t110;
t126 = pkin(2) * t135 + t138;
t125 = pkin(8) * t132;
t124 = (-pkin(3) - pkin(8)) * t111;
t86 = -t102 * t127 + t130;
t99 = t107 * pkin(1);
t123 = pkin(2) * t87 + t86 * qJ(3) + t99;
t122 = g(1) * t107 - g(2) * t111;
t88 = t102 * t129 + t128;
t89 = -t102 * t130 + t127;
t121 = pkin(2) * t89 + t88 * qJ(3) + t136;
t120 = -qJ(3) * t133 + t126;
t119 = pkin(3) * t134 + t121;
t118 = t87 * pkin(9) + t123;
t105 = sin(qJ(4));
t109 = cos(qJ(4));
t76 = t105 * t134 - t109 * t88;
t78 = t105 * t132 + t109 * t86;
t84 = t102 * t105 + t109 * t133;
t117 = g(1) * t76 - g(2) * t78 + g(3) * t84;
t116 = g(1) * t89 + g(2) * t87 + g(3) * t135;
t115 = -g(1) * t88 - g(2) * t86 + g(3) * t133;
t114 = t120 + t139;
t113 = pkin(9) * t89 + t119;
t112 = t101 * t124 + t118;
t108 = cos(qJ(5));
t103 = -qJ(6) - pkin(10);
t96 = pkin(5) * t108 + pkin(4);
t85 = t102 * t109 - t105 * t133;
t80 = -g(3) * t102 - t101 * t122;
t79 = t105 * t86 - t109 * t132;
t77 = t105 * t88 + t109 * t134;
t73 = -g(1) * (t104 * t89 + t108 * t77) - g(2) * (t108 * t79 + t137) - g(3) * (t101 * t131 + t108 * t85);
t72 = -g(1) * (-t104 * t77 + t108 * t89) - g(2) * (-t104 * t79 + t108 * t87) - g(3) * (-t104 * t85 + t108 * t135);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t111 - g(2) * t107, t122, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t116, -t115, t80, -g(1) * t136 - g(2) * (t99 - t125) - g(3) * t138, 0, 0, 0, 0, 0, 0, t80, t116, t115, -g(1) * t121 - g(2) * (t123 - t125) - g(3) * t120, 0, 0, 0, 0, 0, 0, -g(1) * t77 - g(2) * t79 - g(3) * t85, t117, -t116, -g(1) * t113 - g(2) * t112 - g(3) * t114, 0, 0, 0, 0, 0, 0, t73, t72, -t117, -g(1) * (pkin(4) * t77 + pkin(10) * t76 + t113) - g(2) * (t79 * pkin(4) - t78 * pkin(10) + t112) - g(3) * (pkin(4) * t85 + pkin(10) * t84 + t114) 0, 0, 0, 0, 0, 0, t73, t72, -t117, -g(1) * (-t76 * t103 + t77 * t96 + (pkin(5) * t104 + pkin(9)) * t89 + t119) - g(2) * (pkin(5) * t137 + t78 * t103 + t79 * t96 + t118) - g(3) * (-t84 * t103 + t85 * t96 + t126 + t139) + (-g(3) * (pkin(5) * t131 - qJ(3) * t110) - g(2) * t124) * t101;];
U_reg  = t1;
