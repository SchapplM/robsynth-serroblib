% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR5
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:23:09
% EndTime: 2019-03-09 18:23:09
% DurationCPUTime: 0.16s
% Computational Cost: add. (182->76), mult. (192->91), div. (0->0), fcn. (184->10), ass. (0->46)
t92 = qJ(2) + qJ(3);
t86 = sin(t92);
t108 = qJ(4) * t86;
t88 = cos(t92);
t98 = cos(qJ(1));
t119 = t88 * t98;
t125 = pkin(3) * t119 + t98 * t108;
t124 = g(3) * pkin(6);
t93 = sin(qJ(5));
t123 = pkin(5) * t93;
t122 = g(3) * t88;
t94 = sin(qJ(2));
t121 = t94 * pkin(2) + pkin(6);
t95 = sin(qJ(1));
t120 = t88 * t95;
t99 = -pkin(10) - pkin(9);
t118 = t88 * t99;
t91 = qJ(5) + qJ(6);
t85 = sin(t91);
t117 = t95 * t85;
t87 = cos(t91);
t116 = t95 * t87;
t115 = t95 * t93;
t96 = cos(qJ(5));
t114 = t95 * t96;
t113 = t98 * t85;
t112 = t98 * t87;
t111 = t98 * t93;
t110 = t98 * t96;
t100 = -pkin(8) - pkin(7);
t97 = cos(qJ(2));
t83 = t97 * pkin(2) + pkin(1);
t109 = t98 * t100 + t95 * t83;
t107 = t86 * pkin(3) + t121;
t106 = t86 * t111;
t78 = t98 * t83;
t105 = t78 + t125;
t104 = -t95 * t100 + t78;
t103 = pkin(3) * t120 + t95 * t108 + t109;
t102 = g(1) * t98 + g(2) * t95;
t101 = -t88 * qJ(4) + t107;
t82 = t96 * pkin(5) + pkin(4);
t74 = g(1) * t95 - g(2) * t98;
t73 = g(3) * t86 + t102 * t88;
t72 = t102 * t86 - t122;
t1 = [0, 0, 0, 0, 0, 0, -t102, t74, -g(3), -t124, 0, 0, 0, 0, 0, 0, -g(3) * t94 - t102 * t97, -g(3) * t97 + t102 * t94, -t74, -g(1) * (t98 * pkin(1) + t95 * pkin(7)) - g(2) * (t95 * pkin(1) - t98 * pkin(7)) - t124, 0, 0, 0, 0, 0, 0, -t73, t72, -t74, -g(1) * t104 - g(2) * t109 - g(3) * t121, 0, 0, 0, 0, 0, 0, -t74, t73, -t72, -g(1) * (t104 + t125) - g(2) * t103 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (t106 + t114) - g(2) * (t86 * t115 - t110) + t93 * t122, -g(1) * (t86 * t110 - t115) - g(2) * (t86 * t114 + t111) + t96 * t122, -t73, -g(1) * (pkin(9) * t119 + (pkin(4) - t100) * t95 + t105) - g(2) * (-t98 * pkin(4) + pkin(9) * t120 + t103) - g(3) * (pkin(9) * t86 + t101) 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t113 + t116) - g(2) * (t86 * t117 - t112) + t85 * t122, -g(1) * (t86 * t112 - t117) - g(2) * (t86 * t116 + t113) + t87 * t122, -t73, -g(1) * (pkin(5) * t106 - t98 * t118 + t105) - g(2) * (-t98 * t82 + t103) - g(3) * (-t86 * t99 + (-qJ(4) - t123) * t88 + t107) + (-g(1) * (-t100 + t82) - g(2) * (t86 * t123 - t118)) * t95;];
U_reg  = t1;
