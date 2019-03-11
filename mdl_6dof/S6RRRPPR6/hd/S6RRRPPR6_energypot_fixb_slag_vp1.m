% Calculate potential energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:40
% EndTime: 2019-03-09 15:47:40
% DurationCPUTime: 0.41s
% Computational Cost: add. (332->122), mult. (523->157), div. (0->0), fcn. (603->12), ass. (0->48)
t102 = cos(pkin(6));
t132 = t102 * pkin(8) + pkin(7);
t131 = pkin(9) + rSges(4,3);
t130 = pkin(10) + rSges(7,3);
t105 = sin(qJ(3));
t129 = pkin(3) * t105;
t111 = cos(qJ(1));
t101 = sin(pkin(6));
t107 = sin(qJ(1));
t125 = t101 * t107;
t128 = t111 * pkin(1) + pkin(8) * t125;
t127 = rSges(6,3) + qJ(5);
t106 = sin(qJ(2));
t126 = t101 * t106;
t110 = cos(qJ(2));
t124 = t101 * t110;
t123 = t101 * t111;
t122 = t107 * t106;
t121 = t107 * t110;
t120 = t111 * t106;
t119 = t111 * t110;
t118 = t105 * t125;
t103 = -qJ(4) - pkin(9);
t109 = cos(qJ(3));
t94 = t109 * pkin(3) + pkin(2);
t117 = t102 * t129 + t103 * t124 + t94 * t126 + t132;
t83 = t102 * t121 + t120;
t84 = -t102 * t122 + t119;
t116 = pkin(3) * t118 - t83 * t103 + t84 * t94 + t128;
t100 = qJ(3) + pkin(11);
t95 = sin(t100);
t96 = cos(t100);
t78 = t102 * t95 + t96 * t126;
t115 = t78 * pkin(4) + t117;
t73 = t95 * t125 + t84 * t96;
t114 = t73 * pkin(4) + t116;
t81 = -t102 * t119 + t122;
t82 = t102 * t120 + t121;
t98 = t107 * pkin(1);
t113 = t82 * t94 + t98 + (-pkin(8) - t129) * t123 - t81 * t103;
t71 = -t95 * t123 + t82 * t96;
t112 = t71 * pkin(4) + t113;
t108 = cos(qJ(6));
t104 = sin(qJ(6));
t77 = -t102 * t96 + t95 * t126;
t72 = -t96 * t125 + t84 * t95;
t70 = t96 * t123 + t82 * t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t111 * rSges(2,1) - t107 * rSges(2,2)) + g(2) * (t107 * rSges(2,1) + t111 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t128) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t98) + g(3) * (t102 * rSges(3,3) + t132) + (g(1) * rSges(3,3) * t107 + g(3) * (rSges(3,1) * t106 + rSges(3,2) * t110) + g(2) * (-rSges(3,3) - pkin(8)) * t111) * t101) - m(4) * (g(1) * (t84 * pkin(2) + (t84 * t109 + t118) * rSges(4,1) + (-t84 * t105 + t109 * t125) * rSges(4,2) + t131 * t83 + t128) + g(2) * (t82 * pkin(2) + t98 - pkin(8) * t123 + (-t105 * t123 + t82 * t109) * rSges(4,1) + (-t82 * t105 - t109 * t123) * rSges(4,2) + t131 * t81) + g(3) * ((t105 * rSges(4,1) + t109 * rSges(4,2)) * t102 + (-t131 * t110 + (t109 * rSges(4,1) - t105 * rSges(4,2) + pkin(2)) * t106) * t101 + t132)) - m(5) * (g(1) * (t73 * rSges(5,1) - t72 * rSges(5,2) + t83 * rSges(5,3) + t116) + g(2) * (t71 * rSges(5,1) - t70 * rSges(5,2) + t81 * rSges(5,3) + t113) + g(3) * (t78 * rSges(5,1) - t77 * rSges(5,2) - rSges(5,3) * t124 + t117)) - m(6) * (g(1) * (t83 * rSges(6,1) - t73 * rSges(6,2) + t127 * t72 + t114) + g(2) * (t81 * rSges(6,1) - t71 * rSges(6,2) + t127 * t70 + t112) + g(3) * (-rSges(6,1) * t124 - t78 * rSges(6,2) + t127 * t77 + t115)) - m(7) * (g(1) * (t83 * pkin(5) + t72 * qJ(5) + (t72 * t104 + t83 * t108) * rSges(7,1) + (-t83 * t104 + t72 * t108) * rSges(7,2) + t130 * t73 + t114) + g(2) * (t81 * pkin(5) + t70 * qJ(5) + (t70 * t104 + t81 * t108) * rSges(7,1) + (-t81 * t104 + t70 * t108) * rSges(7,2) + t130 * t71 + t112) + g(3) * (-pkin(5) * t124 + t77 * qJ(5) + (t77 * t104 - t108 * t124) * rSges(7,1) + (t104 * t124 + t77 * t108) * rSges(7,2) + t130 * t78 + t115));
U  = t1;
