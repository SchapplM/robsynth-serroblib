% Calculate potential energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:10
% EndTime: 2019-03-08 21:40:10
% DurationCPUTime: 0.33s
% Computational Cost: add. (289->114), mult. (612->145), div. (0->0), fcn. (727->10), ass. (0->51)
t132 = pkin(4) + pkin(8);
t131 = rSges(5,1) + pkin(8);
t130 = rSges(4,3) + pkin(8);
t129 = rSges(6,3) + pkin(9);
t106 = cos(qJ(5));
t128 = t106 * pkin(5) + t132;
t127 = cos(qJ(3));
t100 = sin(pkin(6));
t126 = pkin(7) * t100;
t101 = cos(pkin(10));
t99 = sin(pkin(10));
t125 = t101 * pkin(1) + t99 * t126;
t124 = rSges(5,3) + qJ(4);
t123 = rSges(7,3) + qJ(6) + pkin(9);
t121 = cos(pkin(6));
t122 = t121 * pkin(7) + qJ(1);
t104 = sin(qJ(3));
t120 = t100 * t104;
t105 = sin(qJ(2));
t119 = t100 * t105;
t107 = cos(qJ(2));
t118 = t100 * t107;
t113 = t105 * t121;
t86 = t101 * t107 - t99 * t113;
t117 = t86 * pkin(2) + t125;
t116 = pkin(2) * t119 + t122;
t115 = t100 * t127;
t103 = sin(qJ(5));
t114 = pkin(5) * t103 + qJ(4);
t112 = t107 * t121;
t77 = t99 * t120 + t86 * t127;
t111 = t77 * pkin(3) + t117;
t88 = t121 * t104 + t105 * t115;
t110 = t88 * pkin(3) + t116;
t84 = t101 * t113 + t99 * t107;
t96 = t99 * pkin(1);
t109 = t84 * pkin(2) - t101 * t126 + t96;
t75 = -t101 * t120 + t84 * t127;
t108 = t75 * pkin(3) + t109;
t87 = t104 * t119 - t121 * t127;
t85 = t101 * t105 + t99 * t112;
t83 = -t101 * t112 + t99 * t105;
t79 = t87 * t103 - t106 * t118;
t78 = t103 * t118 + t87 * t106;
t76 = t86 * t104 - t99 * t115;
t74 = t101 * t115 + t84 * t104;
t71 = t76 * t103 + t85 * t106;
t70 = -t85 * t103 + t76 * t106;
t69 = t74 * t103 + t83 * t106;
t68 = -t83 * t103 + t74 * t106;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t101 * rSges(2,1) - t99 * rSges(2,2)) + g(2) * (t99 * rSges(2,1) + t101 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t86 * rSges(3,1) - t85 * rSges(3,2) + t125) + g(2) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t96) + g(3) * (t121 * rSges(3,3) + t122) + (g(1) * rSges(3,3) * t99 + g(3) * (rSges(3,1) * t105 + rSges(3,2) * t107) + g(2) * (-rSges(3,3) - pkin(7)) * t101) * t100) - m(4) * (g(1) * (t77 * rSges(4,1) - t76 * rSges(4,2) + t130 * t85 + t117) + g(2) * (t75 * rSges(4,1) - t74 * rSges(4,2) + t130 * t83 + t109) + g(3) * (t88 * rSges(4,1) - t87 * rSges(4,2) - t130 * t118 + t116)) - m(5) * (g(1) * (-t77 * rSges(5,2) + t124 * t76 + t131 * t85 + t111) + g(2) * (-t75 * rSges(5,2) + t124 * t74 + t131 * t83 + t108) + g(3) * (-t88 * rSges(5,2) - t131 * t118 + t124 * t87 + t110)) - m(6) * (g(1) * (t71 * rSges(6,1) + t70 * rSges(6,2) + t76 * qJ(4) + t129 * t77 + t132 * t85 + t111) + g(2) * (t69 * rSges(6,1) + t68 * rSges(6,2) + t74 * qJ(4) + t129 * t75 + t132 * t83 + t108) + g(3) * (t79 * rSges(6,1) + t78 * rSges(6,2) + t87 * qJ(4) - t132 * t118 + t129 * t88 + t110)) - m(7) * (g(1) * (t71 * rSges(7,1) + t70 * rSges(7,2) + t114 * t76 + t123 * t77 + t128 * t85 + t111) + g(2) * (t69 * rSges(7,1) + t68 * rSges(7,2) + t114 * t74 + t123 * t75 + t128 * t83 + t108) + g(3) * (t79 * rSges(7,1) + t78 * rSges(7,2) + t114 * t87 - t128 * t118 + t123 * t88 + t110));
U  = t1;
