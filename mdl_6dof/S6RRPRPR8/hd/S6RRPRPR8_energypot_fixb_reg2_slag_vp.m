% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:19
% EndTime: 2019-03-09 10:53:19
% DurationCPUTime: 0.20s
% Computational Cost: add. (220->83), mult. (294->108), div. (0->0), fcn. (314->10), ass. (0->44)
t103 = sin(qJ(1));
t105 = cos(qJ(2));
t106 = cos(qJ(1));
t116 = t105 * t106;
t97 = pkin(10) + qJ(4);
t91 = sin(t97);
t92 = cos(t97);
t78 = -t103 * t92 + t91 * t116;
t79 = t103 * t91 + t92 * t116;
t132 = t79 * pkin(4) + t78 * qJ(5);
t117 = t103 * t105;
t76 = t106 * t92 + t91 * t117;
t77 = -t106 * t91 + t92 * t117;
t131 = t77 * pkin(4) + t76 * qJ(5);
t111 = g(1) * t106 + g(2) * t103;
t130 = g(3) * pkin(6);
t102 = sin(qJ(2));
t127 = g(3) * t102;
t125 = t106 * pkin(1) + t103 * pkin(7);
t124 = t102 * t91;
t123 = t102 * t92;
t98 = sin(pkin(10));
t122 = t103 * t98;
t121 = t106 * t98;
t100 = -pkin(8) - qJ(3);
t118 = t100 * t102;
t99 = cos(pkin(10));
t89 = t99 * pkin(3) + pkin(2);
t115 = t105 * t100 + t102 * t89 + pkin(6);
t94 = t103 * pkin(1);
t114 = -t106 * pkin(7) + t94;
t113 = pkin(3) * t122 + t89 * t116 + t125;
t112 = pkin(4) * t123 + qJ(5) * t124 + t115;
t110 = pkin(2) * t105 + qJ(3) * t102;
t109 = -t106 * t118 + t113;
t108 = g(1) * t78 + g(2) * t76 + g(3) * t124;
t81 = t89 * t117;
t107 = -t103 * t118 + t81 + t94 + (-pkin(3) * t98 - pkin(7)) * t106;
t104 = cos(qJ(6));
t101 = sin(qJ(6));
t85 = g(1) * t103 - g(2) * t106;
t80 = -g(3) * t105 + t111 * t102;
t73 = -g(1) * t79 - g(2) * t77 - g(3) * t123;
t1 = [0, 0, 0, 0, 0, 0, -t111, t85, -g(3), -t130, 0, 0, 0, 0, 0, 0, -t111 * t105 - t127, t80, -t85, -g(1) * t125 - g(2) * t114 - t130, 0, 0, 0, 0, 0, 0, -g(1) * (t99 * t116 + t122) - g(2) * (t99 * t117 - t121) - t99 * t127, -g(1) * (t103 * t99 - t98 * t116) - g(2) * (-t106 * t99 - t98 * t117) + t98 * t127, -t80, -g(1) * (t110 * t106 + t125) - g(2) * (t110 * t103 + t114) - g(3) * (t102 * pkin(2) - t105 * qJ(3) + pkin(6)) 0, 0, 0, 0, 0, 0, t73, t108, -t80, -g(1) * t109 - g(2) * t107 - g(3) * t115, 0, 0, 0, 0, 0, 0, t73, -t80, -t108, -g(1) * (t109 + t132) - g(2) * (t107 + t131) - g(3) * t112, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t101 + t79 * t104) - g(2) * (t76 * t101 + t77 * t104) - (t101 * t91 + t104 * t92) * t127, -g(1) * (-t79 * t101 + t78 * t104) - g(2) * (-t77 * t101 + t76 * t104) - (-t101 * t92 + t104 * t91) * t127, t80, -g(1) * (t79 * pkin(5) + t113 + t132) - g(2) * (-pkin(3) * t121 + t77 * pkin(5) + t114 + t131 + t81) - g(3) * (t105 * pkin(9) + t112) + (-g(3) * pkin(5) * t92 + t111 * (pkin(9) + t100)) * t102;];
U_reg  = t1;
