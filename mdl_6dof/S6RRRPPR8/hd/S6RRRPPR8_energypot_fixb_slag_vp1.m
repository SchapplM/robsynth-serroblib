% Calculate potential energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:29
% EndTime: 2019-03-09 16:03:30
% DurationCPUTime: 0.38s
% Computational Cost: add. (278->105), mult. (590->127), div. (0->0), fcn. (697->10), ass. (0->50)
t134 = rSges(6,1) + qJ(4);
t133 = rSges(5,2) + pkin(9);
t132 = rSges(4,3) + pkin(9);
t131 = rSges(6,3) - pkin(9);
t121 = cos(pkin(6));
t130 = t121 * pkin(8) + pkin(7);
t129 = pkin(10) + rSges(7,3);
t128 = cos(qJ(3));
t104 = cos(qJ(1));
t101 = sin(qJ(1));
t97 = sin(pkin(6));
t125 = t101 * t97;
t127 = t104 * pkin(1) + pkin(8) * t125;
t100 = sin(qJ(2));
t126 = t100 * t97;
t103 = cos(qJ(2));
t124 = t103 * t97;
t123 = t104 * t97;
t122 = rSges(5,3) + qJ(4);
t120 = pkin(2) * t126 + t130;
t119 = -qJ(5) - t131;
t115 = t101 * t121;
t86 = -t100 * t115 + t104 * t103;
t118 = t86 * pkin(2) + t127;
t117 = t97 * t128;
t99 = sin(qJ(3));
t82 = t100 * t117 + t121 * t99;
t116 = t82 * pkin(3) + t120;
t114 = t104 * t121;
t76 = t99 * t125 + t86 * t128;
t113 = t76 * pkin(3) + t118;
t112 = t76 * pkin(4) + t113;
t84 = t100 * t114 + t101 * t103;
t95 = t101 * pkin(1);
t111 = t84 * pkin(2) - pkin(8) * t123 + t95;
t110 = t82 * pkin(4) + qJ(5) * t124 + t116;
t102 = cos(qJ(6));
t98 = sin(qJ(6));
t109 = rSges(7,1) * t98 + rSges(7,2) * t102 - pkin(9);
t74 = -t99 * t123 + t84 * t128;
t108 = t74 * pkin(3) + t111;
t107 = rSges(7,1) * t102 - rSges(7,2) * t98 + pkin(5) + qJ(4);
t106 = -qJ(5) - t109;
t105 = t74 * pkin(4) + t108;
t85 = t104 * t100 + t103 * t115;
t83 = t100 * t101 - t103 * t114;
t81 = -t121 * t128 + t99 * t126;
t75 = -t101 * t117 + t86 * t99;
t73 = t104 * t117 + t84 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t104 - t101 * rSges(2,2)) + g(2) * (t101 * rSges(2,1) + rSges(2,2) * t104) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t86 - rSges(3,2) * t85 + t127) + g(2) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t95) + g(3) * (t121 * rSges(3,3) + t130) + (g(1) * rSges(3,3) * t101 + g(3) * (rSges(3,1) * t100 + rSges(3,2) * t103) + g(2) * (-rSges(3,3) - pkin(8)) * t104) * t97) - m(4) * (g(1) * (t76 * rSges(4,1) - t75 * rSges(4,2) + t132 * t85 + t118) + g(2) * (t74 * rSges(4,1) - t73 * rSges(4,2) + t132 * t83 + t111) + g(3) * (t82 * rSges(4,1) - t81 * rSges(4,2) - t132 * t124 + t120)) - m(5) * (g(1) * (rSges(5,1) * t76 + t122 * t75 + t133 * t85 + t113) + g(2) * (t74 * rSges(5,1) + t122 * t73 + t133 * t83 + t108) + g(3) * (rSges(5,1) * t82 + t122 * t81 - t133 * t124 + t116)) - m(6) * (g(3) * (-rSges(6,2) * t82 + t131 * t124 + t134 * t81 + t110) + (-t74 * rSges(6,2) + t119 * t83 + t134 * t73 + t105) * g(2) + (-rSges(6,2) * t76 + t119 * t85 + t134 * t75 + t112) * g(1)) - m(7) * (g(1) * (t106 * t85 + t107 * t75 + t129 * t76 + t112) + g(2) * (t106 * t83 + t107 * t73 + t129 * t74 + t105) + g(3) * (t107 * t81 + t109 * t124 + t129 * t82 + t110));
U  = t1;
