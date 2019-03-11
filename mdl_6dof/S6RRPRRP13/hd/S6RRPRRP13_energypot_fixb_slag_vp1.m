% Calculate potential energy for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP13_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:28
% EndTime: 2019-03-09 12:55:29
% DurationCPUTime: 0.54s
% Computational Cost: add. (260->120), mult. (536->154), div. (0->0), fcn. (622->10), ass. (0->54)
t134 = rSges(7,3) + qJ(6) + pkin(10);
t133 = rSges(6,3) + pkin(10);
t100 = cos(pkin(6));
t132 = t100 * pkin(8) + pkin(7);
t105 = sin(qJ(1));
t131 = g(1) * t105;
t109 = cos(qJ(1));
t130 = g(2) * t109;
t99 = sin(pkin(6));
t126 = t105 * t99;
t129 = t109 * pkin(1) + pkin(8) * t126;
t102 = sin(qJ(5));
t108 = cos(qJ(2));
t119 = t105 * t108;
t104 = sin(qJ(2));
t120 = t104 * t109;
t85 = t100 * t120 + t119;
t128 = t102 * t85;
t127 = t104 * t99;
t125 = t108 * t99;
t124 = t109 * t99;
t123 = qJ(3) * t108;
t122 = t102 * t104;
t121 = t104 * t105;
t118 = t108 * t109;
t117 = pkin(2) * t127 + t132;
t116 = t109 * (-pkin(3) - pkin(8));
t115 = g(2) * t116;
t114 = t100 * pkin(3) + pkin(9) * t127 + t117;
t84 = -t100 * t118 + t121;
t97 = t105 * pkin(1);
t113 = t85 * pkin(2) + t84 * qJ(3) + t97;
t86 = t100 * t119 + t120;
t87 = -t100 * t121 + t118;
t112 = t87 * pkin(2) + qJ(3) * t86 + t129;
t111 = pkin(3) * t126 + t112;
t110 = t85 * pkin(9) + t113;
t107 = cos(qJ(4));
t106 = cos(qJ(5));
t103 = sin(qJ(4));
t94 = pkin(5) * t106 + pkin(4);
t83 = t100 * t107 - t103 * t125;
t82 = t100 * t103 + t107 * t125;
t78 = t84 * t103 - t107 * t124;
t77 = t103 * t124 + t84 * t107;
t76 = t103 * t86 + t107 * t126;
t75 = t103 * t126 - t86 * t107;
t74 = t106 * t83 + t99 * t122;
t73 = -t102 * t83 + t106 * t127;
t72 = t106 * t78 + t128;
t71 = -t102 * t78 + t106 * t85;
t70 = t102 * t87 + t106 * t76;
t69 = -t102 * t76 + t106 * t87;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t109 - t105 * rSges(2,2)) + g(2) * (t105 * rSges(2,1) + rSges(2,2) * t109) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t87 - rSges(3,2) * t86 + t129) + g(2) * (t85 * rSges(3,1) - t84 * rSges(3,2) + t97) + g(3) * (rSges(3,3) * t100 + t132) + (rSges(3,3) * t131 + g(3) * (rSges(3,1) * t104 + rSges(3,2) * t108) + (-rSges(3,3) - pkin(8)) * t130) * t99) - m(4) * (g(1) * (-rSges(4,2) * t87 + rSges(4,3) * t86 + t112) + g(2) * (-t85 * rSges(4,2) + t84 * rSges(4,3) + t113) + g(3) * (rSges(4,1) * t100 + t117) + (rSges(4,1) * t131 + g(3) * (-rSges(4,2) * t104 - rSges(4,3) * t108 - t123) + (-rSges(4,1) - pkin(8)) * t130) * t99) - m(5) * (g(1) * (rSges(5,1) * t76 - rSges(5,2) * t75 + (rSges(5,3) + pkin(9)) * t87 + t111) + g(2) * (t78 * rSges(5,1) + t77 * rSges(5,2) + t85 * rSges(5,3) + t110) + g(3) * (rSges(5,1) * t83 - rSges(5,2) * t82 + t114) + (g(3) * (rSges(5,3) * t104 - t123) + t115) * t99) - m(6) * (g(1) * (rSges(6,1) * t70 + rSges(6,2) * t69 + pkin(4) * t76 + pkin(9) * t87 + t133 * t75 + t111) + g(2) * (t72 * rSges(6,1) + t71 * rSges(6,2) + t78 * pkin(4) + t99 * t116 - t133 * t77 + t110) + g(3) * (rSges(6,1) * t74 + rSges(6,2) * t73 + pkin(4) * t83 - t99 * t123 + t133 * t82 + t114)) - m(7) * (g(1) * (rSges(7,1) * t70 + rSges(7,2) * t69 + t76 * t94 + (pkin(5) * t102 + pkin(9)) * t87 + t134 * t75 + t111) + g(2) * (t72 * rSges(7,1) + t71 * rSges(7,2) + pkin(5) * t128 - t134 * t77 + t78 * t94 + t110) + g(3) * (rSges(7,1) * t74 + rSges(7,2) * t73 + t134 * t82 + t83 * t94 + t114) + (g(3) * (pkin(5) * t122 - t123) + t115) * t99);
U  = t1;
