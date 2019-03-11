% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:03
% EndTime: 2019-03-09 21:15:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (200->77), mult. (268->90), div. (0->0), fcn. (278->8), ass. (0->43)
t92 = sin(qJ(1));
t94 = cos(qJ(2));
t112 = t92 * t94;
t89 = qJ(3) + qJ(4);
t83 = sin(t89);
t84 = cos(t89);
t95 = cos(qJ(1));
t67 = t83 * t112 + t84 * t95;
t110 = t95 * t83;
t68 = t84 * t112 - t110;
t123 = t68 * pkin(4) + t67 * qJ(5);
t101 = g(1) * t95 + g(2) * t92;
t122 = g(3) * pkin(6);
t91 = sin(qJ(2));
t119 = g(3) * t91;
t117 = t83 * t91;
t116 = t84 * t91;
t90 = sin(qJ(3));
t115 = t90 * t92;
t114 = t90 * t95;
t96 = -pkin(9) - pkin(8);
t113 = t91 * t96;
t111 = t94 * t95;
t109 = t95 * pkin(1) + t92 * pkin(7);
t93 = cos(qJ(3));
t81 = pkin(3) * t93 + pkin(2);
t107 = t91 * t81 + t94 * t96 + pkin(6);
t106 = t95 * t113;
t86 = t92 * pkin(1);
t105 = -t95 * pkin(7) + t86;
t104 = pkin(3) * t115 + t81 * t111 + t109;
t103 = pkin(4) * t116 + qJ(5) * t117 + t107;
t102 = pkin(2) * t94 + pkin(8) * t91;
t69 = t94 * t110 - t92 * t84;
t70 = t84 * t111 + t83 * t92;
t100 = t70 * pkin(4) + t69 * qJ(5) + t104;
t99 = g(1) * t69 + g(2) * t67 + g(3) * t117;
t98 = g(1) * t70 + g(2) * t68 + g(3) * t116;
t72 = t81 * t112;
t97 = -t92 * t113 + t72 + t86 + (-pkin(3) * t90 - pkin(7)) * t95;
t74 = g(1) * t92 - g(2) * t95;
t71 = -g(3) * t94 + t101 * t91;
t1 = [0, 0, 0, 0, 0, 0, -t101, t74, -g(3), -t122, 0, 0, 0, 0, 0, 0, -t101 * t94 - t119, t71, -t74, -g(1) * t109 - g(2) * t105 - t122, 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t111 + t115) - g(2) * (t93 * t112 - t114) - t93 * t119, -g(1) * (-t90 * t111 + t92 * t93) - g(2) * (-t90 * t112 - t93 * t95) + t90 * t119, -t71, -g(1) * (t102 * t95 + t109) - g(2) * (t102 * t92 + t105) - g(3) * (pkin(2) * t91 - pkin(8) * t94 + pkin(6)) 0, 0, 0, 0, 0, 0, -t98, t99, -t71, -g(1) * (t104 - t106) - g(2) * t97 - g(3) * t107, 0, 0, 0, 0, 0, 0, -t71, t98, -t99, -g(1) * (t100 - t106) - g(2) * (t97 + t123) - g(3) * t103, 0, 0, 0, 0, 0, 0, -t71, -t99, -t98, -g(1) * (t70 * qJ(6) + t100) - g(2) * (-pkin(3) * t114 + t68 * qJ(6) + t105 + t123 + t72) - g(3) * (-pkin(5) * t94 + t103) + (-g(3) * qJ(6) * t84 - t101 * (pkin(5) - t96)) * t91;];
U_reg  = t1;
