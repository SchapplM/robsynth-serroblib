% Calculate potential energy for
% S6RRRRPP9
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:41
% EndTime: 2019-03-09 21:41:42
% DurationCPUTime: 0.42s
% Computational Cost: add. (335->108), mult. (739->138), div. (0->0), fcn. (903->10), ass. (0->52)
t142 = rSges(7,2) + qJ(5);
t141 = rSges(7,3) + qJ(6);
t140 = rSges(6,1) + pkin(10);
t139 = rSges(4,3) + pkin(9);
t138 = rSges(5,3) + pkin(10);
t137 = cos(qJ(3));
t134 = cos(pkin(6));
t136 = t134 * pkin(8) + pkin(7);
t135 = rSges(6,3) + qJ(5);
t108 = sin(pkin(6));
t111 = sin(qJ(2));
t133 = t108 * t111;
t112 = sin(qJ(1));
t132 = t108 * t112;
t114 = cos(qJ(2));
t131 = t108 * t114;
t115 = cos(qJ(1));
t130 = t108 * t115;
t129 = t115 * pkin(1) + pkin(8) * t132;
t127 = pkin(2) * t133 + t136;
t124 = t112 * t134;
t98 = -t111 * t124 + t115 * t114;
t126 = t98 * pkin(2) + t129;
t125 = t108 * t137;
t123 = t115 * t134;
t106 = t112 * pkin(1);
t96 = t111 * t123 + t112 * t114;
t122 = t96 * pkin(2) - pkin(8) * t130 + t106;
t110 = sin(qJ(3));
t87 = t110 * t132 + t98 * t137;
t97 = t115 * t111 + t114 * t124;
t121 = t87 * pkin(3) + t97 * pkin(9) + t126;
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t78 = t109 * t97 + t113 * t87;
t120 = t78 * pkin(4) + t121;
t94 = t134 * t110 + t111 * t125;
t119 = t94 * pkin(3) - pkin(9) * t131 + t127;
t83 = -t109 * t131 + t113 * t94;
t118 = t83 * pkin(4) + t119;
t85 = -t110 * t130 + t96 * t137;
t95 = t111 * t112 - t114 * t123;
t117 = t85 * pkin(3) + t95 * pkin(9) + t122;
t76 = t109 * t95 + t113 * t85;
t116 = t76 * pkin(4) + t117;
t93 = t110 * t133 - t134 * t137;
t86 = t110 * t98 - t112 * t125;
t84 = t96 * t110 + t115 * t125;
t82 = t109 * t94 + t113 * t131;
t77 = t109 * t87 - t97 * t113;
t75 = t109 * t85 - t95 * t113;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t115 - t112 * rSges(2,2)) + g(2) * (t112 * rSges(2,1) + rSges(2,2) * t115) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t98 - rSges(3,2) * t97 + t129) + g(2) * (t96 * rSges(3,1) - t95 * rSges(3,2) + t106) + g(3) * (t134 * rSges(3,3) + t136) + (g(1) * rSges(3,3) * t112 + g(3) * (rSges(3,1) * t111 + rSges(3,2) * t114) + g(2) * (-rSges(3,3) - pkin(8)) * t115) * t108) - m(4) * (g(1) * (t87 * rSges(4,1) - t86 * rSges(4,2) + t139 * t97 + t126) + g(2) * (t85 * rSges(4,1) - t84 * rSges(4,2) + t139 * t95 + t122) + g(3) * (t94 * rSges(4,1) - t93 * rSges(4,2) - t139 * t131 + t127)) - m(5) * (g(1) * (t78 * rSges(5,1) - t77 * rSges(5,2) + t138 * t86 + t121) + g(2) * (t76 * rSges(5,1) - t75 * rSges(5,2) + t138 * t84 + t117) + g(3) * (t83 * rSges(5,1) - t82 * rSges(5,2) + t138 * t93 + t119)) - m(6) * (g(1) * (-t78 * rSges(6,2) + t135 * t77 + t140 * t86 + t120) + g(2) * (-t76 * rSges(6,2) + t135 * t75 + t140 * t84 + t116) + g(3) * (-t83 * rSges(6,2) + t135 * t82 + t140 * t93 + t118)) - m(7) * (g(1) * (t141 * t78 + t142 * t77 + t120) + g(2) * (t141 * t76 + t142 * t75 + t116) + g(3) * (t141 * t83 + t142 * t82 + t118) + (g(1) * t86 + g(2) * t84 + g(3) * t93) * (rSges(7,1) + pkin(5) + pkin(10)));
U  = t1;
