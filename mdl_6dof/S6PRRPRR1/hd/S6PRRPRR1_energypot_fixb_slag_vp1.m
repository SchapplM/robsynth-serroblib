% Calculate potential energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:50:56
% EndTime: 2019-03-08 21:50:56
% DurationCPUTime: 0.67s
% Computational Cost: add. (357->125), mult. (471->163), div. (0->0), fcn. (534->14), ass. (0->54)
t111 = -qJ(4) - pkin(8);
t138 = -t111 + rSges(5,3);
t137 = pkin(8) + rSges(4,3);
t136 = pkin(10) + rSges(7,3);
t113 = sin(qJ(3));
t135 = pkin(3) * t113;
t107 = sin(pkin(11));
t134 = g(1) * t107;
t109 = cos(pkin(11));
t133 = g(2) * t109;
t116 = cos(qJ(3));
t97 = t116 * pkin(3) + pkin(2);
t108 = sin(pkin(6));
t131 = t107 * t108;
t132 = t109 * pkin(1) + pkin(7) * t131;
t130 = t108 * t109;
t129 = t108 * t113;
t114 = sin(qJ(2));
t128 = t108 * t114;
t127 = t108 * t116;
t117 = cos(qJ(2));
t126 = t108 * t117;
t110 = cos(pkin(6));
t125 = t110 * t114;
t124 = t110 * t117;
t123 = t110 * pkin(7) + qJ(1);
t106 = qJ(3) + pkin(12);
t105 = -pkin(9) + t111;
t85 = t107 * t124 + t109 * t114;
t86 = -t107 * t125 + t109 * t117;
t99 = cos(t106);
t89 = pkin(4) * t99 + t97;
t98 = sin(t106);
t90 = pkin(4) * t98 + t135;
t122 = -t85 * t105 + t90 * t131 + t86 * t89 + t132;
t121 = t105 * t126 + t110 * t90 + t89 * t128 + t123;
t120 = rSges(5,1) * t99 - rSges(5,2) * t98 + t97;
t119 = rSges(5,1) * t98 + rSges(5,2) * t99 + t135;
t101 = t107 * pkin(1);
t83 = t107 * t114 - t109 * t124;
t84 = t107 * t117 + t109 * t125;
t118 = t101 + t84 * t89 + (-pkin(7) - t90) * t130 - t83 * t105;
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t100 = qJ(5) + t106;
t96 = cos(t100);
t95 = sin(t100);
t78 = t110 * t95 + t96 * t128;
t77 = -t110 * t96 + t95 * t128;
t74 = t95 * t131 + t86 * t96;
t73 = -t96 * t131 + t86 * t95;
t72 = -t95 * t130 + t84 * t96;
t71 = t96 * t130 + t84 * t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t109 - rSges(2,2) * t107) + g(2) * (rSges(2,1) * t107 + rSges(2,2) * t109) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t86 - rSges(3,2) * t85 + t132) + g(2) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t101) + g(3) * (t110 * rSges(3,3) + t123) + (rSges(3,3) * t134 + g(3) * (rSges(3,1) * t114 + rSges(3,2) * t117) + (-rSges(3,3) - pkin(7)) * t133) * t108) - m(4) * (g(1) * (t86 * pkin(2) + (t107 * t129 + t116 * t86) * rSges(4,1) + (t107 * t127 - t113 * t86) * rSges(4,2) + t137 * t85 + t132) + g(2) * (t84 * pkin(2) + t101 - pkin(7) * t130 + (-t109 * t129 + t116 * t84) * rSges(4,1) + (-t109 * t127 - t113 * t84) * rSges(4,2) + t137 * t83) + g(3) * ((rSges(4,1) * t113 + rSges(4,2) * t116) * t110 + (-t137 * t117 + (rSges(4,1) * t116 - rSges(4,2) * t113 + pkin(2)) * t114) * t108 + t123)) - m(5) * ((t119 * t134 + (-pkin(7) - t119) * t133) * t108 + (t123 + t119 * t110 + (t120 * t114 - t117 * t138) * t108) * g(3) + (t120 * t84 + t138 * t83 + t101) * g(2) + (t120 * t86 + t138 * t85 + t132) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t74 - rSges(6,2) * t73 + rSges(6,3) * t85 + t122) + g(2) * (rSges(6,1) * t72 - rSges(6,2) * t71 + rSges(6,3) * t83 + t118) + g(3) * (t78 * rSges(6,1) - t77 * rSges(6,2) - rSges(6,3) * t126 + t121)) - m(7) * (g(1) * (t74 * pkin(5) + (t112 * t85 + t115 * t74) * rSges(7,1) + (-t112 * t74 + t115 * t85) * rSges(7,2) + t136 * t73 + t122) + g(2) * (t72 * pkin(5) + (t112 * t83 + t115 * t72) * rSges(7,1) + (-t112 * t72 + t115 * t83) * rSges(7,2) + t136 * t71 + t118) + g(3) * (t78 * pkin(5) + (-t112 * t126 + t78 * t115) * rSges(7,1) + (-t78 * t112 - t115 * t126) * rSges(7,2) + t136 * t77 + t121));
U  = t1;
