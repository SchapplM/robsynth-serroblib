% Calculate potential energy for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR13_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:50:58
% EndTime: 2019-03-09 23:50:58
% DurationCPUTime: 0.46s
% Computational Cost: add. (361->118), mult. (807->156), div. (0->0), fcn. (997->12), ass. (0->52)
t140 = rSges(6,2) + pkin(10);
t139 = rSges(4,3) + pkin(9);
t138 = rSges(5,3) + pkin(10);
t137 = cos(qJ(3));
t134 = cos(pkin(6));
t136 = t134 * pkin(8) + pkin(7);
t135 = rSges(6,3) + qJ(5);
t106 = sin(pkin(6));
t110 = sin(qJ(2));
t133 = t106 * t110;
t111 = sin(qJ(1));
t132 = t106 * t111;
t114 = cos(qJ(2));
t131 = t106 * t114;
t115 = cos(qJ(1));
t130 = t106 * t115;
t129 = t115 * pkin(1) + pkin(8) * t132;
t127 = pkin(2) * t133 + t136;
t124 = t111 * t134;
t97 = -t110 * t124 + t115 * t114;
t126 = t97 * pkin(2) + t129;
t125 = t106 * t137;
t123 = t115 * t134;
t104 = t111 * pkin(1);
t95 = t110 * t123 + t111 * t114;
t122 = t95 * pkin(2) - pkin(8) * t130 + t104;
t109 = sin(qJ(3));
t86 = t109 * t132 + t137 * t97;
t96 = t115 * t110 + t114 * t124;
t121 = t86 * pkin(3) + t96 * pkin(9) + t126;
t108 = sin(qJ(4));
t113 = cos(qJ(4));
t77 = t108 * t96 + t113 * t86;
t120 = t77 * pkin(4) + t121;
t93 = t109 * t134 + t110 * t125;
t119 = t93 * pkin(3) - pkin(9) * t131 + t127;
t82 = -t108 * t131 + t113 * t93;
t118 = t82 * pkin(4) + t119;
t84 = -t109 * t130 + t137 * t95;
t94 = t110 * t111 - t114 * t123;
t117 = t84 * pkin(3) + t94 * pkin(9) + t122;
t75 = t108 * t94 + t113 * t84;
t116 = t75 * pkin(4) + t117;
t112 = cos(qJ(6));
t107 = sin(qJ(6));
t92 = t109 * t133 - t134 * t137;
t85 = t109 * t97 - t111 * t125;
t83 = t95 * t109 + t115 * t125;
t81 = t108 * t93 + t113 * t131;
t76 = t108 * t86 - t96 * t113;
t74 = t108 * t84 - t94 * t113;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t115 - t111 * rSges(2,2)) + g(2) * (t111 * rSges(2,1) + rSges(2,2) * t115) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t97 - rSges(3,2) * t96 + t129) + g(2) * (t95 * rSges(3,1) - t94 * rSges(3,2) + t104) + g(3) * (t134 * rSges(3,3) + t136) + (g(1) * rSges(3,3) * t111 + g(3) * (rSges(3,1) * t110 + rSges(3,2) * t114) + g(2) * (-rSges(3,3) - pkin(8)) * t115) * t106) - m(4) * (g(1) * (rSges(4,1) * t86 - rSges(4,2) * t85 + t139 * t96 + t126) + g(2) * (t84 * rSges(4,1) - t83 * rSges(4,2) + t139 * t94 + t122) + g(3) * (rSges(4,1) * t93 - rSges(4,2) * t92 - t131 * t139 + t127)) - m(5) * (g(1) * (rSges(5,1) * t77 - rSges(5,2) * t76 + t138 * t85 + t121) + g(2) * (t75 * rSges(5,1) - t74 * rSges(5,2) + t138 * t83 + t117) + g(3) * (rSges(5,1) * t82 - rSges(5,2) * t81 + t138 * t92 + t119)) - m(6) * (g(1) * (rSges(6,1) * t77 + t135 * t76 + t140 * t85 + t120) + g(2) * (t75 * rSges(6,1) + t135 * t74 + t140 * t83 + t116) + g(3) * (rSges(6,1) * t82 + t135 * t81 + t140 * t92 + t118)) - m(7) * (g(1) * (t77 * pkin(5) + t76 * qJ(5) + (t107 * t76 + t112 * t77) * rSges(7,1) + (-t107 * t77 + t112 * t76) * rSges(7,2) + t120) + g(2) * (t75 * pkin(5) + t74 * qJ(5) + (t107 * t74 + t112 * t75) * rSges(7,1) + (-t107 * t75 + t112 * t74) * rSges(7,2) + t116) + g(3) * (t82 * pkin(5) + t81 * qJ(5) + (t107 * t81 + t112 * t82) * rSges(7,1) + (-t107 * t82 + t112 * t81) * rSges(7,2) + t118) + (g(1) * t85 + g(2) * t83 + g(3) * t92) * (-pkin(11) + pkin(10) - rSges(7,3)));
U  = t1;
