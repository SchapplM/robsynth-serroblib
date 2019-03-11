% Calculate potential energy for
% S6RRRPRP11
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:28
% EndTime: 2019-03-09 17:41:28
% DurationCPUTime: 0.33s
% Computational Cost: add. (289->114), mult. (612->145), div. (0->0), fcn. (727->10), ass. (0->51)
t132 = pkin(4) + pkin(9);
t131 = rSges(5,1) + pkin(9);
t130 = rSges(4,3) + pkin(9);
t129 = rSges(6,3) + pkin(10);
t118 = cos(pkin(6));
t128 = t118 * pkin(8) + pkin(7);
t105 = cos(qJ(5));
t127 = t105 * pkin(5) + t132;
t126 = cos(qJ(3));
t107 = cos(qJ(1));
t104 = sin(qJ(1));
t99 = sin(pkin(6));
t123 = t104 * t99;
t125 = t107 * pkin(1) + pkin(8) * t123;
t103 = sin(qJ(2));
t124 = t103 * t99;
t106 = cos(qJ(2));
t122 = t106 * t99;
t121 = t107 * t99;
t120 = rSges(5,3) + qJ(4);
t119 = rSges(7,3) + qJ(6) + pkin(10);
t117 = pkin(2) * t124 + t128;
t112 = t104 * t118;
t88 = -t103 * t112 + t107 * t106;
t116 = t88 * pkin(2) + t125;
t115 = t99 * t126;
t102 = sin(qJ(3));
t84 = t118 * t102 + t103 * t115;
t114 = t84 * pkin(3) + t117;
t101 = sin(qJ(5));
t113 = pkin(5) * t101 + qJ(4);
t111 = t107 * t118;
t79 = t102 * t123 + t88 * t126;
t110 = t79 * pkin(3) + t116;
t86 = t103 * t111 + t104 * t106;
t97 = t104 * pkin(1);
t109 = t86 * pkin(2) - pkin(8) * t121 + t97;
t77 = -t102 * t121 + t86 * t126;
t108 = t77 * pkin(3) + t109;
t87 = t107 * t103 + t106 * t112;
t85 = t104 * t103 - t106 * t111;
t83 = t102 * t124 - t118 * t126;
t78 = t88 * t102 - t104 * t115;
t76 = t86 * t102 + t107 * t115;
t75 = t83 * t101 - t105 * t122;
t74 = t101 * t122 + t83 * t105;
t71 = t78 * t101 + t87 * t105;
t70 = -t87 * t101 + t78 * t105;
t69 = t76 * t101 + t85 * t105;
t68 = -t85 * t101 + t76 * t105;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t107 * rSges(2,1) - t104 * rSges(2,2)) + g(2) * (t104 * rSges(2,1) + t107 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t88 * rSges(3,1) - t87 * rSges(3,2) + t125) + g(2) * (t86 * rSges(3,1) - t85 * rSges(3,2) + t97) + g(3) * (t118 * rSges(3,3) + t128) + (g(1) * rSges(3,3) * t104 + g(3) * (rSges(3,1) * t103 + rSges(3,2) * t106) + g(2) * (-rSges(3,3) - pkin(8)) * t107) * t99) - m(4) * (g(1) * (t79 * rSges(4,1) - t78 * rSges(4,2) + t130 * t87 + t116) + g(2) * (t77 * rSges(4,1) - t76 * rSges(4,2) + t130 * t85 + t109) + g(3) * (t84 * rSges(4,1) - t83 * rSges(4,2) - t130 * t122 + t117)) - m(5) * (g(1) * (-t79 * rSges(5,2) + t120 * t78 + t131 * t87 + t110) + g(2) * (-t77 * rSges(5,2) + t120 * t76 + t131 * t85 + t108) + g(3) * (-t84 * rSges(5,2) + t120 * t83 - t131 * t122 + t114)) - m(6) * (g(1) * (t71 * rSges(6,1) + t70 * rSges(6,2) + t78 * qJ(4) + t129 * t79 + t132 * t87 + t110) + g(2) * (t69 * rSges(6,1) + t68 * rSges(6,2) + t76 * qJ(4) + t129 * t77 + t132 * t85 + t108) + g(3) * (t75 * rSges(6,1) + t74 * rSges(6,2) + t83 * qJ(4) - t132 * t122 + t129 * t84 + t114)) - m(7) * (g(1) * (t71 * rSges(7,1) + t70 * rSges(7,2) + t113 * t78 + t119 * t79 + t127 * t87 + t110) + g(2) * (t69 * rSges(7,1) + t68 * rSges(7,2) + t113 * t76 + t119 * t77 + t127 * t85 + t108) + g(3) * (t75 * rSges(7,1) + t74 * rSges(7,2) + t113 * t83 + t119 * t84 - t127 * t122 + t114));
U  = t1;
