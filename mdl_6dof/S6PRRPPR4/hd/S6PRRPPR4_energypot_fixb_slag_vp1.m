% Calculate potential energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:25
% EndTime: 2019-03-08 21:12:25
% DurationCPUTime: 0.49s
% Computational Cost: add. (361->118), mult. (807->156), div. (0->0), fcn. (997->12), ass. (0->52)
t140 = rSges(4,3) + pkin(8);
t139 = cos(qJ(3));
t108 = sin(pkin(6));
t138 = pkin(7) * t108;
t137 = rSges(6,2) + qJ(4);
t136 = rSges(5,3) + qJ(4);
t135 = rSges(6,3) + qJ(5);
t134 = cos(pkin(6));
t112 = sin(qJ(3));
t133 = t108 * t112;
t113 = sin(qJ(2));
t132 = t108 * t113;
t115 = cos(qJ(2));
t131 = t108 * t115;
t130 = t134 * pkin(7) + qJ(1);
t107 = sin(pkin(10));
t110 = cos(pkin(10));
t129 = t110 * pkin(1) + t107 * t138;
t124 = t113 * t134;
t95 = -t107 * t124 + t110 * t115;
t127 = t95 * pkin(2) + t129;
t126 = pkin(2) * t132 + t130;
t125 = t108 * t139;
t123 = t115 * t134;
t103 = t107 * pkin(1);
t93 = t107 * t115 + t110 * t124;
t122 = t93 * pkin(2) - t110 * t138 + t103;
t86 = t107 * t133 + t139 * t95;
t94 = t107 * t123 + t110 * t113;
t121 = t86 * pkin(3) + t94 * pkin(8) + t127;
t106 = sin(pkin(11));
t109 = cos(pkin(11));
t77 = t106 * t94 + t109 * t86;
t120 = t77 * pkin(4) + t121;
t97 = t112 * t134 + t113 * t125;
t119 = t97 * pkin(3) - pkin(8) * t131 + t126;
t82 = -t106 * t131 + t97 * t109;
t118 = t82 * pkin(4) + t119;
t84 = -t110 * t133 + t139 * t93;
t92 = t107 * t113 - t110 * t123;
t117 = t84 * pkin(3) + t92 * pkin(8) + t122;
t75 = t106 * t92 + t109 * t84;
t116 = t75 * pkin(4) + t117;
t114 = cos(qJ(6));
t111 = sin(qJ(6));
t96 = t112 * t132 - t134 * t139;
t85 = -t107 * t125 + t112 * t95;
t83 = t110 * t125 + t93 * t112;
t81 = t97 * t106 + t109 * t131;
t76 = t106 * t86 - t94 * t109;
t74 = t106 * t84 - t92 * t109;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t110 - rSges(2,2) * t107) + g(2) * (rSges(2,1) * t107 + rSges(2,2) * t110) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t95 - rSges(3,2) * t94 + t129) + g(2) * (rSges(3,1) * t93 - rSges(3,2) * t92 + t103) + g(3) * (t134 * rSges(3,3) + t130) + (g(1) * rSges(3,3) * t107 + g(3) * (rSges(3,1) * t113 + rSges(3,2) * t115) + g(2) * (-rSges(3,3) - pkin(7)) * t110) * t108) - m(4) * (g(1) * (rSges(4,1) * t86 - rSges(4,2) * t85 + t140 * t94 + t127) + g(2) * (rSges(4,1) * t84 - rSges(4,2) * t83 + t140 * t92 + t122) + g(3) * (t97 * rSges(4,1) - t96 * rSges(4,2) - t140 * t131 + t126)) - m(5) * (g(1) * (rSges(5,1) * t77 - rSges(5,2) * t76 + t136 * t85 + t121) + g(2) * (rSges(5,1) * t75 - rSges(5,2) * t74 + t136 * t83 + t117) + g(3) * (t82 * rSges(5,1) - t81 * rSges(5,2) + t136 * t96 + t119)) - m(6) * (g(1) * (rSges(6,1) * t77 + t135 * t76 + t137 * t85 + t120) + g(2) * (rSges(6,1) * t75 + t135 * t74 + t137 * t83 + t116) + g(3) * (t82 * rSges(6,1) + t135 * t81 + t137 * t96 + t118)) - m(7) * (g(1) * (t77 * pkin(5) + t76 * qJ(5) + (t111 * t76 + t114 * t77) * rSges(7,1) + (-t111 * t77 + t114 * t76) * rSges(7,2) + t120) + g(2) * (t75 * pkin(5) + t74 * qJ(5) + (t111 * t74 + t114 * t75) * rSges(7,1) + (-t111 * t75 + t114 * t74) * rSges(7,2) + t116) + g(3) * (t82 * pkin(5) + t81 * qJ(5) + (t111 * t81 + t114 * t82) * rSges(7,1) + (-t111 * t82 + t114 * t81) * rSges(7,2) + t118) + (g(1) * t85 + g(2) * t83 + g(3) * t96) * (-pkin(9) + qJ(4) - rSges(7,3)));
U  = t1;
