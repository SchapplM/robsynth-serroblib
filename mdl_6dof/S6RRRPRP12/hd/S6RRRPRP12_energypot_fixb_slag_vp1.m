% Calculate potential energy for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:32
% EndTime: 2019-03-09 17:51:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (301->108), mult. (650->137), div. (0->0), fcn. (780->10), ass. (0->54)
t142 = pkin(4) + pkin(9);
t141 = rSges(5,1) + pkin(9);
t140 = rSges(7,1) + pkin(5);
t139 = rSges(7,2) + pkin(10);
t138 = rSges(4,3) + pkin(9);
t137 = rSges(6,3) + pkin(10);
t136 = cos(qJ(3));
t132 = cos(pkin(6));
t135 = t132 * pkin(8) + pkin(7);
t134 = rSges(5,3) + qJ(4);
t133 = rSges(7,3) + qJ(6);
t107 = sin(pkin(6));
t110 = sin(qJ(2));
t131 = t107 * t110;
t111 = sin(qJ(1));
t130 = t107 * t111;
t113 = cos(qJ(2));
t129 = t107 * t113;
t114 = cos(qJ(1));
t128 = t107 * t114;
t127 = t114 * pkin(1) + pkin(8) * t130;
t126 = pkin(2) * t131 + t135;
t123 = t111 * t132;
t97 = -t110 * t123 + t114 * t113;
t125 = t97 * pkin(2) + t127;
t124 = t107 * t136;
t122 = t114 * t132;
t109 = sin(qJ(3));
t93 = t132 * t109 + t110 * t124;
t121 = t93 * pkin(3) + t126;
t86 = t109 * t130 + t97 * t136;
t120 = t86 * pkin(3) + t125;
t105 = t111 * pkin(1);
t95 = t110 * t122 + t111 * t113;
t119 = t95 * pkin(2) - pkin(8) * t128 + t105;
t84 = -t109 * t128 + t95 * t136;
t118 = t84 * pkin(3) + t119;
t85 = t97 * t109 - t111 * t124;
t96 = t114 * t110 + t113 * t123;
t117 = t85 * qJ(4) + t142 * t96 + t120;
t92 = t109 * t131 - t132 * t136;
t116 = t92 * qJ(4) - t142 * t129 + t121;
t83 = t95 * t109 + t114 * t124;
t94 = t111 * t110 - t113 * t122;
t115 = t83 * qJ(4) + t142 * t94 + t118;
t112 = cos(qJ(5));
t108 = sin(qJ(5));
t82 = t92 * t108 - t112 * t129;
t81 = t108 * t129 + t92 * t112;
t76 = t85 * t108 + t96 * t112;
t75 = t96 * t108 - t85 * t112;
t74 = t83 * t108 + t94 * t112;
t73 = t94 * t108 - t83 * t112;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t114 * rSges(2,1) - t111 * rSges(2,2)) + g(2) * (t111 * rSges(2,1) + t114 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t97 * rSges(3,1) - t96 * rSges(3,2) + t127) + g(2) * (t95 * rSges(3,1) - t94 * rSges(3,2) + t105) + g(3) * (t132 * rSges(3,3) + t135) + (g(1) * rSges(3,3) * t111 + g(3) * (rSges(3,1) * t110 + rSges(3,2) * t113) + g(2) * (-rSges(3,3) - pkin(8)) * t114) * t107) - m(4) * (g(1) * (t86 * rSges(4,1) - t85 * rSges(4,2) + t138 * t96 + t125) + g(2) * (t84 * rSges(4,1) - t83 * rSges(4,2) + t138 * t94 + t119) + g(3) * (t93 * rSges(4,1) - t92 * rSges(4,2) - t138 * t129 + t126)) - m(5) * (g(1) * (-t86 * rSges(5,2) + t134 * t85 + t141 * t96 + t120) + g(2) * (-t84 * rSges(5,2) + t134 * t83 + t141 * t94 + t118) + g(3) * (-t93 * rSges(5,2) - t141 * t129 + t134 * t92 + t121)) - m(6) * (g(1) * (t76 * rSges(6,1) - t75 * rSges(6,2) + t137 * t86 + t117) + g(2) * (t74 * rSges(6,1) - t73 * rSges(6,2) + t137 * t84 + t115) + g(3) * (t82 * rSges(6,1) + t81 * rSges(6,2) + t137 * t93 + t116)) - m(7) * (g(1) * (t133 * t75 + t139 * t86 + t140 * t76 + t117) + g(2) * (t133 * t73 + t139 * t84 + t140 * t74 + t115) + g(3) * (-t133 * t81 + t139 * t93 + t140 * t82 + t116));
U  = t1;
