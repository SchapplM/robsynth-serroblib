% Calculate potential energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:22
% EndTime: 2019-03-09 16:21:22
% DurationCPUTime: 0.36s
% Computational Cost: add. (301->105), mult. (612->129), div. (0->0), fcn. (727->12), ass. (0->50)
t130 = pkin(4) + pkin(9);
t129 = rSges(5,1) + pkin(9);
t128 = rSges(4,3) + pkin(9);
t117 = cos(pkin(6));
t127 = t117 * pkin(8) + pkin(7);
t126 = cos(qJ(3));
t103 = cos(qJ(1));
t101 = sin(qJ(1));
t96 = sin(pkin(6));
t122 = t101 * t96;
t125 = t103 * pkin(1) + pkin(8) * t122;
t124 = pkin(10) + qJ(5) + rSges(7,3);
t100 = sin(qJ(2));
t123 = t100 * t96;
t102 = cos(qJ(2));
t121 = t102 * t96;
t120 = t103 * t96;
t119 = rSges(5,3) + qJ(4);
t118 = qJ(5) + rSges(6,3);
t116 = pkin(2) * t123 + t127;
t112 = t101 * t117;
t81 = -t100 * t112 + t103 * t102;
t115 = t81 * pkin(2) + t125;
t114 = t96 * t126;
t99 = sin(qJ(3));
t77 = t100 * t114 + t117 * t99;
t113 = t77 * pkin(3) + t116;
t111 = t103 * t117;
t72 = t122 * t99 + t126 * t81;
t110 = t72 * pkin(3) + t115;
t79 = t100 * t111 + t101 * t102;
t92 = t101 * pkin(1);
t109 = t79 * pkin(2) - pkin(8) * t120 + t92;
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t108 = t95 * rSges(6,1) + t97 * rSges(6,2) + qJ(4);
t70 = -t120 * t99 + t126 * t79;
t107 = t70 * pkin(3) + t109;
t106 = t97 * rSges(6,1) - t95 * rSges(6,2) + t130;
t94 = pkin(11) + qJ(6);
t89 = sin(t94);
t90 = cos(t94);
t105 = t90 * rSges(7,1) - t89 * rSges(7,2) + pkin(5) * t97 + t130;
t104 = t89 * rSges(7,1) + t90 * rSges(7,2) + t95 * pkin(5) + qJ(4);
t80 = t103 * t100 + t102 * t112;
t78 = t100 * t101 - t102 * t111;
t76 = -t117 * t126 + t123 * t99;
t71 = -t101 * t114 + t81 * t99;
t69 = t103 * t114 + t79 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t103 - t101 * rSges(2,2)) + g(2) * (t101 * rSges(2,1) + rSges(2,2) * t103) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t81 - rSges(3,2) * t80 + t125) + g(2) * (t79 * rSges(3,1) - t78 * rSges(3,2) + t92) + g(3) * (t117 * rSges(3,3) + t127) + (g(1) * rSges(3,3) * t101 + g(3) * (rSges(3,1) * t100 + rSges(3,2) * t102) + g(2) * (-rSges(3,3) - pkin(8)) * t103) * t96) - m(4) * (g(1) * (rSges(4,1) * t72 - rSges(4,2) * t71 + t128 * t80 + t115) + g(2) * (t70 * rSges(4,1) - t69 * rSges(4,2) + t128 * t78 + t109) + g(3) * (rSges(4,1) * t77 - rSges(4,2) * t76 - t121 * t128 + t116)) - m(5) * (g(1) * (-rSges(5,2) * t72 + t119 * t71 + t129 * t80 + t110) + g(2) * (-t70 * rSges(5,2) + t119 * t69 + t129 * t78 + t107) + g(3) * (-rSges(5,2) * t77 + t119 * t76 - t121 * t129 + t113)) - m(6) * (g(1) * (t106 * t80 + t108 * t71 + t118 * t72 + t110) + g(2) * (t106 * t78 + t108 * t69 + t118 * t70 + t107) + g(3) * (-t106 * t121 + t108 * t76 + t118 * t77 + t113)) - m(7) * (g(1) * (t104 * t71 + t105 * t80 + t124 * t72 + t110) + g(2) * (t104 * t69 + t105 * t78 + t124 * t70 + t107) + g(3) * (t104 * t76 - t105 * t121 + t124 * t77 + t113));
U  = t1;
