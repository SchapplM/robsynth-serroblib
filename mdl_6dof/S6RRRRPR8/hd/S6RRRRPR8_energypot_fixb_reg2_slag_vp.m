% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:44:47
% EndTime: 2019-03-09 22:44:47
% DurationCPUTime: 0.20s
% Computational Cost: add. (220->83), mult. (294->108), div. (0->0), fcn. (314->10), ass. (0->45)
t101 = sin(qJ(1));
t104 = cos(qJ(2));
t105 = cos(qJ(1));
t117 = t104 * t105;
t97 = qJ(3) + qJ(4);
t91 = sin(t97);
t92 = cos(t97);
t78 = -t101 * t92 + t91 * t117;
t79 = t101 * t91 + t92 * t117;
t133 = t79 * pkin(4) + t78 * qJ(5);
t118 = t101 * t104;
t76 = t105 * t92 + t91 * t118;
t77 = -t105 * t91 + t92 * t118;
t132 = t77 * pkin(4) + t76 * qJ(5);
t110 = g(1) * t105 + g(2) * t101;
t131 = g(3) * pkin(6);
t100 = sin(qJ(2));
t128 = g(3) * t100;
t126 = t105 * pkin(1) + t101 * pkin(7);
t125 = t100 * t91;
t124 = t100 * t92;
t99 = sin(qJ(3));
t123 = t101 * t99;
t122 = t105 * t99;
t106 = -pkin(9) - pkin(8);
t119 = t100 * t106;
t103 = cos(qJ(3));
t116 = t105 * t103;
t89 = t103 * pkin(3) + pkin(2);
t115 = t100 * t89 + t104 * t106 + pkin(6);
t94 = t101 * pkin(1);
t114 = -t105 * pkin(7) + t94;
t113 = pkin(3) * t123 + t89 * t117 + t126;
t112 = pkin(4) * t124 + qJ(5) * t125 + t115;
t111 = pkin(2) * t104 + pkin(8) * t100;
t109 = -t105 * t119 + t113;
t108 = g(1) * t78 + g(2) * t76 + g(3) * t125;
t81 = t89 * t118;
t107 = -t101 * t119 + t81 + t94 + (-pkin(3) * t99 - pkin(7)) * t105;
t102 = cos(qJ(6));
t98 = sin(qJ(6));
t83 = g(1) * t101 - g(2) * t105;
t80 = -g(3) * t104 + t110 * t100;
t73 = -g(1) * t79 - g(2) * t77 - g(3) * t124;
t1 = [0, 0, 0, 0, 0, 0, -t110, t83, -g(3), -t131, 0, 0, 0, 0, 0, 0, -t110 * t104 - t128, t80, -t83, -g(1) * t126 - g(2) * t114 - t131, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t116 + t123) - g(2) * (t103 * t118 - t122) - t103 * t128, -g(1) * (t101 * t103 - t99 * t117) - g(2) * (-t99 * t118 - t116) + t99 * t128, -t80, -g(1) * (t111 * t105 + t126) - g(2) * (t111 * t101 + t114) - g(3) * (t100 * pkin(2) - t104 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, t73, t108, -t80, -g(1) * t109 - g(2) * t107 - g(3) * t115, 0, 0, 0, 0, 0, 0, t73, -t80, -t108, -g(1) * (t109 + t133) - g(2) * (t107 + t132) - g(3) * t112, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t102 + t78 * t98) - g(2) * (t77 * t102 + t76 * t98) - (t102 * t92 + t91 * t98) * t128, -g(1) * (t78 * t102 - t79 * t98) - g(2) * (t76 * t102 - t77 * t98) - (t102 * t91 - t92 * t98) * t128, t80, -g(1) * (t79 * pkin(5) + t113 + t133) - g(2) * (-pkin(3) * t122 + t77 * pkin(5) + t114 + t132 + t81) - g(3) * (t104 * pkin(10) + t112) + (-g(3) * pkin(5) * t92 + t110 * (pkin(10) + t106)) * t100;];
U_reg  = t1;
