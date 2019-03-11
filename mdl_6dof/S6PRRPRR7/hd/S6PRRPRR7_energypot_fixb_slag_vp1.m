% Calculate potential energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:08
% EndTime: 2019-03-08 22:30:09
% DurationCPUTime: 0.36s
% Computational Cost: add. (301->105), mult. (612->129), div. (0->0), fcn. (727->12), ass. (0->50)
t130 = pkin(4) + pkin(8);
t96 = sin(pkin(6));
t129 = pkin(7) * t96;
t128 = rSges(5,1) + pkin(8);
t127 = rSges(4,3) + pkin(8);
t126 = pkin(9) + rSges(6,3);
t125 = cos(qJ(3));
t99 = sin(qJ(3));
t124 = t96 * t99;
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t123 = t97 * pkin(1) + t95 * t129;
t100 = sin(qJ(2));
t122 = t100 * t96;
t102 = cos(qJ(2));
t121 = t96 * t102;
t120 = rSges(5,3) + qJ(4);
t117 = cos(pkin(6));
t119 = t117 * pkin(7) + qJ(1);
t118 = pkin(10) + pkin(9) + rSges(7,3);
t113 = t100 * t117;
t79 = t97 * t102 - t113 * t95;
t116 = t79 * pkin(2) + t123;
t115 = pkin(2) * t122 + t119;
t114 = t96 * t125;
t112 = t102 * t117;
t72 = t124 * t95 + t125 * t79;
t111 = t72 * pkin(3) + t116;
t81 = t100 * t114 + t117 * t99;
t110 = t81 * pkin(3) + t115;
t77 = t95 * t102 + t113 * t97;
t91 = t95 * pkin(1);
t109 = t77 * pkin(2) - t129 * t97 + t91;
t70 = -t124 * t97 + t125 * t77;
t108 = t70 * pkin(3) + t109;
t101 = cos(qJ(5));
t98 = sin(qJ(5));
t107 = t98 * rSges(6,1) + t101 * rSges(6,2) + qJ(4);
t106 = t101 * rSges(6,1) - t98 * rSges(6,2) + t130;
t94 = qJ(5) + qJ(6);
t89 = sin(t94);
t90 = cos(t94);
t105 = t90 * rSges(7,1) - t89 * rSges(7,2) + pkin(5) * t101 + t130;
t104 = t89 * rSges(7,1) + t90 * rSges(7,2) + t98 * pkin(5) + qJ(4);
t80 = -t117 * t125 + t122 * t99;
t78 = t97 * t100 + t112 * t95;
t76 = t100 * t95 - t112 * t97;
t71 = -t114 * t95 + t79 * t99;
t69 = t114 * t97 + t77 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t97 - rSges(2,2) * t95) + g(2) * (rSges(2,1) * t95 + rSges(2,2) * t97) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t79 - rSges(3,2) * t78 + t123) + g(2) * (rSges(3,1) * t77 - rSges(3,2) * t76 + t91) + g(3) * (t117 * rSges(3,3) + t119) + (g(1) * rSges(3,3) * t95 + g(3) * (rSges(3,1) * t100 + rSges(3,2) * t102) + g(2) * (-rSges(3,3) - pkin(7)) * t97) * t96) - m(4) * (g(1) * (rSges(4,1) * t72 - rSges(4,2) * t71 + t127 * t78 + t116) + g(2) * (rSges(4,1) * t70 - rSges(4,2) * t69 + t127 * t76 + t109) + g(3) * (rSges(4,1) * t81 - rSges(4,2) * t80 - t121 * t127 + t115)) - m(5) * (g(1) * (-rSges(5,2) * t72 + t120 * t71 + t128 * t78 + t111) + g(2) * (-rSges(5,2) * t70 + t120 * t69 + t128 * t76 + t108) + g(3) * (-rSges(5,2) * t81 + t120 * t80 - t121 * t128 + t110)) - m(6) * (g(1) * (t106 * t78 + t107 * t71 + t126 * t72 + t111) + g(2) * (t106 * t76 + t107 * t69 + t126 * t70 + t108) + g(3) * (-t106 * t121 + t107 * t80 + t126 * t81 + t110)) - m(7) * (g(1) * (t104 * t71 + t105 * t78 + t118 * t72 + t111) + g(2) * (t104 * t69 + t105 * t76 + t118 * t70 + t108) + g(3) * (t104 * t80 - t105 * t121 + t118 * t81 + t110));
U  = t1;
