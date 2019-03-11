% Calculate potential energy for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:00:11
% EndTime: 2019-03-10 02:00:11
% DurationCPUTime: 0.40s
% Computational Cost: add. (323->125), mult. (604->162), div. (0->0), fcn. (717->12), ass. (0->52)
t102 = sin(qJ(4));
t114 = t102 * pkin(4) + pkin(9);
t109 = -pkin(11) - pkin(10);
t130 = rSges(4,3) + pkin(9);
t121 = cos(pkin(6));
t129 = t121 * pkin(8) + pkin(7);
t100 = qJ(4) + qJ(5);
t93 = sin(t100);
t128 = pkin(5) * t93 + t114;
t127 = pkin(10) + rSges(5,3);
t106 = cos(qJ(4));
t92 = t106 * pkin(4) + pkin(3);
t126 = cos(qJ(3));
t124 = rSges(7,3) + qJ(6) - t109;
t108 = cos(qJ(1));
t101 = sin(pkin(6));
t105 = sin(qJ(1));
t119 = t101 * t105;
t123 = t108 * pkin(1) + pkin(8) * t119;
t122 = rSges(6,3) - t109;
t104 = sin(qJ(2));
t120 = t101 * t104;
t107 = cos(qJ(2));
t118 = t101 * t107;
t117 = t101 * t108;
t116 = pkin(2) * t120 + t129;
t112 = t105 * t121;
t84 = -t104 * t112 + t108 * t107;
t115 = t84 * pkin(2) + t123;
t113 = t101 * t126;
t111 = t108 * t121;
t82 = t104 * t111 + t105 * t107;
t96 = t105 * pkin(1);
t110 = t82 * pkin(2) - pkin(8) * t117 + t96;
t103 = sin(qJ(3));
t94 = cos(t100);
t85 = pkin(5) * t94 + t92;
t83 = t108 * t104 + t107 * t112;
t81 = t105 * t104 - t107 * t111;
t80 = t121 * t103 + t104 * t113;
t79 = t103 * t120 - t121 * t126;
t76 = t103 * t119 + t84 * t126;
t75 = t84 * t103 - t105 * t113;
t74 = -t103 * t117 + t82 * t126;
t73 = t82 * t103 + t108 * t113;
t72 = -t118 * t93 + t80 * t94;
t71 = -t118 * t94 - t80 * t93;
t70 = t76 * t94 + t83 * t93;
t69 = -t76 * t93 + t83 * t94;
t68 = t74 * t94 + t81 * t93;
t67 = -t74 * t93 + t81 * t94;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t108 * rSges(2,1) - t105 * rSges(2,2)) + g(2) * (t105 * rSges(2,1) + t108 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t123) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t96) + g(3) * (t121 * rSges(3,3) + t129) + (g(1) * rSges(3,3) * t105 + g(3) * (rSges(3,1) * t104 + rSges(3,2) * t107) + g(2) * (-rSges(3,3) - pkin(8)) * t108) * t101) - m(4) * (g(1) * (t76 * rSges(4,1) - t75 * rSges(4,2) + t130 * t83 + t115) + g(2) * (t74 * rSges(4,1) - t73 * rSges(4,2) + t130 * t81 + t110) + g(3) * (t80 * rSges(4,1) - t79 * rSges(4,2) - t130 * t118 + t116)) - m(5) * (g(1) * (t76 * pkin(3) + t83 * pkin(9) + (t83 * t102 + t76 * t106) * rSges(5,1) + (-t76 * t102 + t83 * t106) * rSges(5,2) + t127 * t75 + t115) + g(2) * (t74 * pkin(3) + t81 * pkin(9) + (t81 * t102 + t74 * t106) * rSges(5,1) + (-t74 * t102 + t81 * t106) * rSges(5,2) + t127 * t73 + t110) + g(3) * (t80 * pkin(3) - pkin(9) * t118 + (-t102 * t118 + t80 * t106) * rSges(5,1) + (-t80 * t102 - t106 * t118) * rSges(5,2) + t127 * t79 + t116)) - m(6) * (g(1) * (t70 * rSges(6,1) + t69 * rSges(6,2) + t114 * t83 + t122 * t75 + t76 * t92 + t115) + g(2) * (t68 * rSges(6,1) + t67 * rSges(6,2) + t114 * t81 + t122 * t73 + t74 * t92 + t110) + g(3) * (t72 * rSges(6,1) + t71 * rSges(6,2) - t114 * t118 + t122 * t79 + t80 * t92 + t116)) - m(7) * (g(1) * (t70 * rSges(7,1) + t69 * rSges(7,2) + t124 * t75 + t128 * t83 + t76 * t85 + t115) + g(2) * (t68 * rSges(7,1) + t67 * rSges(7,2) + t124 * t73 + t128 * t81 + t74 * t85 + t110) + g(3) * (t72 * rSges(7,1) + t71 * rSges(7,2) - t128 * t118 + t124 * t79 + t80 * t85 + t116));
U  = t1;
