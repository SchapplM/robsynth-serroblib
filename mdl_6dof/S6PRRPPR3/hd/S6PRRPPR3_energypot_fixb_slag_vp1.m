% Calculate potential energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:03
% EndTime: 2019-03-08 21:08:03
% DurationCPUTime: 0.38s
% Computational Cost: add. (278->105), mult. (590->127), div. (0->0), fcn. (697->10), ass. (0->50)
t134 = rSges(6,1) + qJ(4);
t98 = sin(pkin(6));
t133 = pkin(7) * t98;
t132 = rSges(5,2) + pkin(8);
t131 = rSges(4,3) + pkin(8);
t130 = rSges(6,3) - pkin(8);
t129 = pkin(9) + rSges(7,3);
t128 = cos(qJ(3));
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t127 = t99 * pkin(1) + t97 * t133;
t101 = sin(qJ(3));
t126 = t101 * t98;
t102 = sin(qJ(2));
t125 = t102 * t98;
t104 = cos(qJ(2));
t124 = t104 * t98;
t123 = rSges(5,3) + qJ(4);
t121 = cos(pkin(6));
t122 = t121 * pkin(7) + qJ(1);
t120 = -qJ(5) - t130;
t116 = t102 * t121;
t84 = t99 * t104 - t97 * t116;
t119 = t84 * pkin(2) + t127;
t118 = pkin(2) * t125 + t122;
t117 = t98 * t128;
t115 = t104 * t121;
t76 = t97 * t126 + t84 * t128;
t114 = t76 * pkin(3) + t119;
t86 = t121 * t101 + t102 * t117;
t113 = t86 * pkin(3) + t118;
t112 = t76 * pkin(4) + t114;
t82 = t97 * t104 + t99 * t116;
t94 = t97 * pkin(1);
t111 = t82 * pkin(2) - t99 * t133 + t94;
t110 = t86 * pkin(4) + qJ(5) * t124 + t113;
t100 = sin(qJ(6));
t103 = cos(qJ(6));
t109 = rSges(7,1) * t100 + rSges(7,2) * t103 - pkin(8);
t74 = -t99 * t126 + t82 * t128;
t108 = t74 * pkin(3) + t111;
t107 = t74 * pkin(4) + t108;
t106 = rSges(7,1) * t103 - rSges(7,2) * t100 + pkin(5) + qJ(4);
t105 = -qJ(5) - t109;
t85 = t101 * t125 - t121 * t128;
t83 = t99 * t102 + t97 * t115;
t81 = t102 * t97 - t99 * t115;
t75 = t101 * t84 - t97 * t117;
t73 = t82 * t101 + t99 * t117;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t99 - rSges(2,2) * t97) + g(2) * (rSges(2,1) * t97 + rSges(2,2) * t99) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t127) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t94) + g(3) * (t121 * rSges(3,3) + t122) + (g(1) * rSges(3,3) * t97 + g(3) * (rSges(3,1) * t102 + rSges(3,2) * t104) + g(2) * (-rSges(3,3) - pkin(7)) * t99) * t98) - m(4) * (g(1) * (rSges(4,1) * t76 - rSges(4,2) * t75 + t131 * t83 + t119) + g(2) * (rSges(4,1) * t74 - rSges(4,2) * t73 + t131 * t81 + t111) + g(3) * (t86 * rSges(4,1) - t85 * rSges(4,2) - t131 * t124 + t118)) - m(5) * (g(1) * (rSges(5,1) * t76 + t123 * t75 + t132 * t83 + t114) + g(2) * (rSges(5,1) * t74 + t123 * t73 + t132 * t81 + t108) + g(3) * (t86 * rSges(5,1) + t123 * t85 - t132 * t124 + t113)) - m(6) * (g(3) * (-t86 * rSges(6,2) + t130 * t124 + t134 * t85 + t110) + (-rSges(6,2) * t74 + t120 * t81 + t134 * t73 + t107) * g(2) + (-rSges(6,2) * t76 + t120 * t83 + t134 * t75 + t112) * g(1)) - m(7) * (g(1) * (t105 * t83 + t106 * t75 + t129 * t76 + t112) + g(2) * (t105 * t81 + t106 * t73 + t129 * t74 + t107) + g(3) * (t106 * t85 + t109 * t124 + t129 * t86 + t110));
U  = t1;
