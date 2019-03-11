% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:21:55
% EndTime: 2019-03-09 19:21:55
% DurationCPUTime: 0.21s
% Computational Cost: add. (182->85), mult. (350->114), div. (0->0), fcn. (388->10), ass. (0->47)
t102 = cos(qJ(3));
t99 = sin(qJ(2));
t123 = t102 * t99;
t98 = sin(qJ(3));
t126 = t98 * t99;
t131 = pkin(3) * t123 + qJ(4) * t126;
t130 = g(3) * pkin(6);
t129 = g(3) * t99;
t128 = t99 * pkin(2) + pkin(6);
t97 = sin(qJ(5));
t127 = t97 * t98;
t100 = sin(qJ(1));
t104 = cos(qJ(1));
t125 = t104 * pkin(1) + t100 * pkin(7);
t124 = t100 * t99;
t122 = t104 * t99;
t105 = -pkin(10) - pkin(9);
t121 = t105 * t99;
t103 = cos(qJ(2));
t120 = t100 * t103;
t119 = t103 * t104;
t118 = t104 * t102;
t117 = pkin(5) * t97 + qJ(4);
t116 = t100 * pkin(1) - t104 * pkin(7);
t115 = t128 + t131;
t114 = pkin(2) * t119 + pkin(8) * t122 + t125;
t113 = -t103 * pkin(8) + t128;
t79 = t100 * t98 + t103 * t118;
t112 = t79 * pkin(3) + t114;
t111 = g(1) * t104 + g(2) * t100;
t110 = pkin(2) * t120 + pkin(8) * t124 + t116;
t77 = t102 * t120 - t104 * t98;
t109 = t77 * pkin(3) + t110;
t78 = -t100 * t102 + t98 * t119;
t108 = t78 * qJ(4) + t112;
t76 = t98 * t120 + t118;
t107 = g(1) * t78 + g(2) * t76 + g(3) * t126;
t106 = t76 * qJ(4) + t109;
t101 = cos(qJ(5));
t96 = qJ(5) + qJ(6);
t90 = cos(t96);
t89 = sin(t96);
t88 = t101 * pkin(5) + pkin(4);
t80 = g(1) * t100 - g(2) * t104;
t73 = -g(3) * t103 + t111 * t99;
t72 = -g(1) * t79 - g(2) * t77 - g(3) * t123;
t1 = [0, 0, 0, 0, 0, 0, -t111, t80, -g(3), -t130, 0, 0, 0, 0, 0, 0, -t111 * t103 - t129, t73, -t80, -g(1) * t125 - g(2) * t116 - t130, 0, 0, 0, 0, 0, 0, t72, t107, -t73, -g(1) * t114 - g(2) * t110 - g(3) * t113, 0, 0, 0, 0, 0, 0, t72, -t73, -t107, -g(1) * t108 - g(2) * t106 - g(3) * (t113 + t131) 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t101 + t78 * t97) - g(2) * (t77 * t101 + t76 * t97) - (t101 * t102 + t127) * t129, -g(1) * (t78 * t101 - t79 * t97) - g(2) * (t76 * t101 - t77 * t97) - (t101 * t98 - t102 * t97) * t129, t73, -g(1) * (t79 * pkin(4) - pkin(9) * t122 + t108) - g(2) * (t77 * pkin(4) - pkin(9) * t124 + t106) - g(3) * (pkin(4) * t123 + (-pkin(8) + pkin(9)) * t103 + t115) 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t89 + t79 * t90) - g(2) * (t76 * t89 + t77 * t90) - (t102 * t90 + t89 * t98) * t129, -g(1) * (t78 * t90 - t79 * t89) - g(2) * (t76 * t90 - t77 * t89) - (-t102 * t89 + t90 * t98) * t129, t73, -g(1) * (t104 * t121 + t117 * t78 + t79 * t88 + t112) - g(2) * (t100 * t121 + t117 * t76 + t77 * t88 + t109) - g(3) * ((pkin(5) * t127 + t102 * t88) * t99 + (-pkin(8) - t105) * t103 + t115);];
U_reg  = t1;
