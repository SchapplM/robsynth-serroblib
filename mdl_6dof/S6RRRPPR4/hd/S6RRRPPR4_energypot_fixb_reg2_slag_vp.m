% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:35:59
% EndTime: 2019-03-09 15:35:59
% DurationCPUTime: 0.19s
% Computational Cost: add. (220->83), mult. (294->109), div. (0->0), fcn. (314->10), ass. (0->45)
t104 = cos(qJ(1));
t100 = sin(qJ(1));
t103 = cos(qJ(2));
t119 = t100 * t103;
t95 = qJ(3) + pkin(10);
t89 = sin(t95);
t90 = cos(t95);
t74 = t104 * t90 + t89 * t119;
t75 = -t104 * t89 + t90 * t119;
t128 = t75 * pkin(4) + t74 * qJ(5);
t127 = g(3) * pkin(6);
t99 = sin(qJ(2));
t126 = g(3) * t99;
t125 = t89 * t99;
t124 = t90 * t99;
t123 = t104 * pkin(1) + t100 * pkin(7);
t98 = sin(qJ(3));
t122 = t100 * t98;
t121 = t100 * t99;
t118 = t103 * t104;
t102 = cos(qJ(3));
t117 = t104 * t102;
t88 = t102 * pkin(3) + pkin(2);
t96 = -qJ(4) - pkin(8);
t116 = t103 * t96 + t99 * t88 + pkin(6);
t115 = t104 * t99 * t96;
t114 = -pkin(3) * t98 - pkin(7);
t92 = t100 * pkin(1);
t113 = -t104 * pkin(7) + t92;
t112 = pkin(3) * t122 + t88 * t118 + t123;
t111 = pkin(4) * t124 + qJ(5) * t125 + t116;
t110 = pkin(2) * t103 + pkin(8) * t99;
t109 = g(1) * t104 + g(2) * t100;
t108 = t88 * t119 - t96 * t121 + t92;
t76 = -t100 * t90 + t89 * t118;
t77 = t100 * t89 + t90 * t118;
t107 = t77 * pkin(4) + t76 * qJ(5) + t112;
t106 = g(1) * t76 + g(2) * t74 + g(3) * t125;
t105 = t114 * t104 + t108;
t101 = cos(qJ(6));
t97 = sin(qJ(6));
t82 = g(1) * t100 - g(2) * t104;
t78 = -g(3) * t103 + t109 * t99;
t71 = -g(1) * t77 - g(2) * t75 - g(3) * t124;
t1 = [0, 0, 0, 0, 0, 0, -t109, t82, -g(3), -t127, 0, 0, 0, 0, 0, 0, -t109 * t103 - t126, t78, -t82, -g(1) * t123 - g(2) * t113 - t127, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t117 + t122) - g(2) * (t102 * t119 - t104 * t98) - t102 * t126, -g(1) * (t100 * t102 - t98 * t118) - g(2) * (-t98 * t119 - t117) + t98 * t126, -t78, -g(1) * (t110 * t104 + t123) - g(2) * (t110 * t100 + t113) - g(3) * (t99 * pkin(2) - t103 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, t71, t106, -t78, -g(1) * (t112 - t115) - g(2) * t105 - g(3) * t116, 0, 0, 0, 0, 0, 0, t71, -t78, -t106, -g(1) * (t107 - t115) - g(2) * (t105 + t128) - g(3) * t111, 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t101 + t76 * t97) - g(2) * (t75 * t101 + t74 * t97) - (t101 * t90 + t89 * t97) * t126, -g(1) * (t76 * t101 - t77 * t97) - g(2) * (t74 * t101 - t75 * t97) - (t101 * t89 - t90 * t97) * t126, t78, -g(1) * (t77 * pkin(5) + t107) - g(2) * (t75 * pkin(5) - pkin(9) * t121 + t108 + t128) - g(3) * (pkin(5) * t124 + t103 * pkin(9) + t111) + (-g(2) * t114 - g(1) * (-pkin(9) - t96) * t99) * t104;];
U_reg  = t1;
