% Calculate potential energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:12
% EndTime: 2019-03-09 21:30:13
% DurationCPUTime: 0.43s
% Computational Cost: add. (335->108), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->52)
t140 = rSges(7,1) + pkin(5);
t139 = rSges(7,2) + qJ(5);
t138 = rSges(6,2) + pkin(10);
t137 = rSges(4,3) + pkin(9);
t136 = rSges(5,3) + pkin(10);
t135 = cos(qJ(3));
t132 = cos(pkin(6));
t134 = t132 * pkin(8) + pkin(7);
t133 = rSges(6,3) + qJ(5);
t106 = sin(pkin(6));
t109 = sin(qJ(2));
t131 = t106 * t109;
t110 = sin(qJ(1));
t130 = t106 * t110;
t112 = cos(qJ(2));
t129 = t106 * t112;
t113 = cos(qJ(1));
t128 = t106 * t113;
t127 = t113 * pkin(1) + pkin(8) * t130;
t125 = pkin(2) * t131 + t134;
t122 = t110 * t132;
t97 = -t109 * t122 + t113 * t112;
t124 = t97 * pkin(2) + t127;
t123 = t106 * t135;
t121 = t113 * t132;
t104 = t110 * pkin(1);
t95 = t109 * t121 + t110 * t112;
t120 = t95 * pkin(2) - pkin(8) * t128 + t104;
t108 = sin(qJ(3));
t86 = t108 * t130 + t97 * t135;
t96 = t113 * t109 + t112 * t122;
t119 = t86 * pkin(3) + t96 * pkin(9) + t124;
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t77 = t107 * t96 + t111 * t86;
t118 = t77 * pkin(4) + t119;
t93 = t132 * t108 + t109 * t123;
t117 = t93 * pkin(3) - pkin(9) * t129 + t125;
t82 = -t107 * t129 + t111 * t93;
t116 = t82 * pkin(4) + t117;
t84 = -t108 * t128 + t95 * t135;
t94 = t109 * t110 - t112 * t121;
t115 = t84 * pkin(3) + t94 * pkin(9) + t120;
t75 = t107 * t94 + t111 * t84;
t114 = t75 * pkin(4) + t115;
t92 = t108 * t131 - t132 * t135;
t85 = t108 * t97 - t110 * t123;
t83 = t95 * t108 + t113 * t123;
t81 = t107 * t93 + t111 * t129;
t76 = t107 * t86 - t96 * t111;
t74 = t107 * t84 - t94 * t111;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t113 - t110 * rSges(2,2)) + g(2) * (t110 * rSges(2,1) + rSges(2,2) * t113) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t97 - rSges(3,2) * t96 + t127) + g(2) * (t95 * rSges(3,1) - t94 * rSges(3,2) + t104) + g(3) * (t132 * rSges(3,3) + t134) + (g(1) * rSges(3,3) * t110 + g(3) * (rSges(3,1) * t109 + rSges(3,2) * t112) + g(2) * (-rSges(3,3) - pkin(8)) * t113) * t106) - m(4) * (g(1) * (t86 * rSges(4,1) - t85 * rSges(4,2) + t137 * t96 + t124) + g(2) * (t84 * rSges(4,1) - t83 * rSges(4,2) + t137 * t94 + t120) + g(3) * (t93 * rSges(4,1) - t92 * rSges(4,2) - t137 * t129 + t125)) - m(5) * (g(1) * (t77 * rSges(5,1) - t76 * rSges(5,2) + t136 * t85 + t119) + g(2) * (t75 * rSges(5,1) - t74 * rSges(5,2) + t136 * t83 + t115) + g(3) * (t82 * rSges(5,1) - t81 * rSges(5,2) + t136 * t92 + t117)) - m(6) * (g(1) * (t77 * rSges(6,1) + t133 * t76 + t138 * t85 + t118) + g(2) * (t75 * rSges(6,1) + t133 * t74 + t138 * t83 + t114) + g(3) * (t82 * rSges(6,1) + t133 * t81 + t138 * t92 + t116)) - m(7) * (g(1) * (t139 * t76 + t140 * t77 + t118) + g(2) * (t139 * t74 + t140 * t75 + t114) + g(3) * (t139 * t81 + t140 * t82 + t116) + (g(1) * t85 + g(2) * t83 + g(3) * t92) * (-rSges(7,3) + pkin(10) - qJ(6)));
U  = t1;
