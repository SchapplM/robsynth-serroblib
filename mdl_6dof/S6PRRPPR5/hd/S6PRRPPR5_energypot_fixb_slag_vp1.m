% Calculate potential energy for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:05
% EndTime: 2019-03-08 21:18:05
% DurationCPUTime: 0.37s
% Computational Cost: add. (301->105), mult. (612->129), div. (0->0), fcn. (727->12), ass. (0->50)
t130 = pkin(4) + pkin(8);
t97 = sin(pkin(6));
t129 = pkin(7) * t97;
t128 = rSges(5,1) + pkin(8);
t127 = rSges(4,3) + pkin(8);
t126 = cos(qJ(3));
t96 = sin(pkin(10));
t99 = cos(pkin(10));
t125 = t99 * pkin(1) + t96 * t129;
t101 = sin(qJ(3));
t124 = t101 * t97;
t102 = sin(qJ(2));
t123 = t102 * t97;
t103 = cos(qJ(2));
t122 = t103 * t97;
t121 = rSges(5,3) + qJ(4);
t117 = cos(pkin(6));
t120 = t117 * pkin(7) + qJ(1);
t119 = qJ(5) + rSges(6,3);
t118 = pkin(9) + qJ(5) + rSges(7,3);
t113 = t102 * t117;
t79 = t99 * t103 - t96 * t113;
t116 = t79 * pkin(2) + t125;
t115 = pkin(2) * t123 + t120;
t114 = t97 * t126;
t112 = t103 * t117;
t72 = t96 * t124 + t79 * t126;
t111 = t72 * pkin(3) + t116;
t81 = t117 * t101 + t102 * t114;
t110 = t81 * pkin(3) + t115;
t77 = t96 * t103 + t99 * t113;
t91 = t96 * pkin(1);
t109 = t77 * pkin(2) - t99 * t129 + t91;
t95 = sin(pkin(11));
t98 = cos(pkin(11));
t108 = t95 * rSges(6,1) + t98 * rSges(6,2) + qJ(4);
t70 = -t99 * t124 + t77 * t126;
t107 = t70 * pkin(3) + t109;
t106 = t98 * rSges(6,1) - t95 * rSges(6,2) + t130;
t94 = pkin(11) + qJ(6);
t89 = sin(t94);
t90 = cos(t94);
t105 = t90 * rSges(7,1) - t89 * rSges(7,2) + pkin(5) * t98 + t130;
t104 = t89 * rSges(7,1) + t90 * rSges(7,2) + t95 * pkin(5) + qJ(4);
t80 = t101 * t123 - t117 * t126;
t78 = t99 * t102 + t96 * t112;
t76 = t102 * t96 - t99 * t112;
t71 = t101 * t79 - t96 * t114;
t69 = t77 * t101 + t99 * t114;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t99 - rSges(2,2) * t96) + g(2) * (rSges(2,1) * t96 + rSges(2,2) * t99) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t79 - rSges(3,2) * t78 + t125) + g(2) * (rSges(3,1) * t77 - rSges(3,2) * t76 + t91) + g(3) * (t117 * rSges(3,3) + t120) + (g(1) * rSges(3,3) * t96 + g(3) * (rSges(3,1) * t102 + rSges(3,2) * t103) + g(2) * (-rSges(3,3) - pkin(7)) * t99) * t97) - m(4) * (g(1) * (rSges(4,1) * t72 - rSges(4,2) * t71 + t127 * t78 + t116) + g(2) * (rSges(4,1) * t70 - rSges(4,2) * t69 + t127 * t76 + t109) + g(3) * (t81 * rSges(4,1) - t80 * rSges(4,2) - t127 * t122 + t115)) - m(5) * (g(1) * (-rSges(5,2) * t72 + t121 * t71 + t128 * t78 + t111) + g(2) * (-rSges(5,2) * t70 + t121 * t69 + t128 * t76 + t107) + g(3) * (-t81 * rSges(5,2) + t121 * t80 - t128 * t122 + t110)) - m(6) * (g(1) * (t106 * t78 + t108 * t71 + t119 * t72 + t111) + g(2) * (t106 * t76 + t108 * t69 + t119 * t70 + t107) + g(3) * (-t106 * t122 + t108 * t80 + t119 * t81 + t110)) - m(7) * (g(1) * (t104 * t71 + t105 * t78 + t118 * t72 + t111) + g(2) * (t104 * t69 + t105 * t76 + t118 * t70 + t107) + g(3) * (t104 * t80 - t105 * t122 + t118 * t81 + t110));
U  = t1;
