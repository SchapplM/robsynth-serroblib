% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:39
% EndTime: 2019-03-09 17:19:39
% DurationCPUTime: 0.18s
% Computational Cost: add. (170->74), mult. (350->95), div. (0->0), fcn. (388->8), ass. (0->44)
t96 = sin(qJ(2));
t99 = cos(qJ(3));
t118 = t96 * t99;
t95 = sin(qJ(3));
t120 = t95 * t96;
t125 = pkin(3) * t118 + qJ(4) * t120;
t124 = g(3) * pkin(6);
t123 = g(3) * t96;
t122 = t96 * pkin(2) + pkin(6);
t94 = sin(qJ(5));
t121 = t94 * t95;
t97 = sin(qJ(1));
t119 = t96 * t97;
t101 = cos(qJ(1));
t117 = t101 * pkin(1) + t97 * pkin(7);
t100 = cos(qJ(2));
t116 = t100 * t97;
t115 = t101 * t96;
t114 = t100 * t101;
t113 = pkin(5) * t94 + qJ(4);
t112 = t97 * pkin(1) - t101 * pkin(7);
t111 = t122 + t125;
t110 = pkin(2) * t114 + pkin(8) * t115 + t117;
t109 = -t100 * pkin(8) + t122;
t78 = t99 * t114 + t97 * t95;
t108 = t78 * pkin(3) + t110;
t107 = g(1) * t101 + g(2) * t97;
t106 = pkin(2) * t116 + pkin(8) * t119 + t112;
t76 = -t101 * t95 + t99 * t116;
t105 = t76 * pkin(3) + t106;
t77 = t95 * t114 - t97 * t99;
t104 = t77 * qJ(4) + t108;
t75 = t101 * t99 + t95 * t116;
t103 = g(1) * t77 + g(2) * t75 + g(3) * t120;
t102 = t75 * qJ(4) + t105;
t98 = cos(qJ(5));
t93 = -qJ(6) - pkin(9);
t87 = pkin(5) * t98 + pkin(4);
t79 = g(1) * t97 - g(2) * t101;
t72 = -g(3) * t100 + t107 * t96;
t71 = -g(1) * t78 - g(2) * t76 - g(3) * t118;
t70 = -g(1) * (t77 * t94 + t78 * t98) - g(2) * (t75 * t94 + t76 * t98) - (t98 * t99 + t121) * t123;
t69 = -g(1) * (t77 * t98 - t78 * t94) - g(2) * (t75 * t98 - t76 * t94) - (-t94 * t99 + t95 * t98) * t123;
t1 = [0, 0, 0, 0, 0, 0, -t107, t79, -g(3), -t124, 0, 0, 0, 0, 0, 0, -t107 * t100 - t123, t72, -t79, -g(1) * t117 - g(2) * t112 - t124, 0, 0, 0, 0, 0, 0, t71, t103, -t72, -g(1) * t110 - g(2) * t106 - g(3) * t109, 0, 0, 0, 0, 0, 0, t71, -t72, -t103, -g(1) * t104 - g(2) * t102 - g(3) * (t109 + t125) 0, 0, 0, 0, 0, 0, t70, t69, t72, -g(1) * (t78 * pkin(4) - pkin(9) * t115 + t104) - g(2) * (t76 * pkin(4) - pkin(9) * t119 + t102) - g(3) * (pkin(4) * t118 + (-pkin(8) + pkin(9)) * t100 + t111) 0, 0, 0, 0, 0, 0, t70, t69, t72, -g(1) * (t113 * t77 + t93 * t115 + t78 * t87 + t108) - g(2) * (t113 * t75 + t93 * t119 + t76 * t87 + t105) - g(3) * ((pkin(5) * t121 + t87 * t99) * t96 + (-pkin(8) - t93) * t100 + t111);];
U_reg  = t1;
