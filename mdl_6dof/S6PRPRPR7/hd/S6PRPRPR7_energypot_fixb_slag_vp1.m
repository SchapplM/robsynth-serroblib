% Calculate potential energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:05
% EndTime: 2019-03-08 19:51:06
% DurationCPUTime: 0.61s
% Computational Cost: add. (249->115), mult. (514->143), div. (0->0), fcn. (592->10), ass. (0->49)
t131 = pkin(9) + rSges(7,3);
t130 = rSges(6,3) + qJ(5);
t95 = sin(pkin(10));
t129 = g(1) * t95;
t97 = cos(pkin(10));
t128 = g(2) * t97;
t96 = sin(pkin(6));
t127 = t95 * t96;
t126 = t97 * pkin(1) + pkin(7) * t127;
t100 = sin(qJ(4));
t125 = t100 * t96;
t101 = sin(qJ(2));
t124 = t101 * t95;
t123 = t101 * t96;
t103 = cos(qJ(4));
t122 = t103 * t96;
t104 = cos(qJ(2));
t121 = t104 * t96;
t98 = cos(pkin(6));
t120 = t104 * t98;
t119 = t97 * t101;
t118 = t98 * pkin(7) + qJ(1);
t117 = qJ(3) * t104;
t116 = pkin(2) * t123 + t118;
t115 = (-pkin(3) - pkin(7)) * t128;
t78 = -t97 * t120 + t124;
t79 = t104 * t95 + t98 * t119;
t91 = t95 * pkin(1);
t114 = t79 * pkin(2) + t78 * qJ(3) + t91;
t113 = t98 * pkin(3) + pkin(8) * t123 + t116;
t80 = t95 * t120 + t119;
t81 = t104 * t97 - t98 * t124;
t112 = t81 * pkin(2) + t80 * qJ(3) + t126;
t83 = -t100 * t121 + t103 * t98;
t111 = t83 * pkin(4) + t113;
t102 = cos(qJ(6));
t99 = sin(qJ(6));
t110 = t102 * rSges(7,1) - t99 * rSges(7,2) + pkin(5);
t109 = t99 * rSges(7,1) + t102 * rSges(7,2) + qJ(5);
t108 = pkin(3) * t127 + t112;
t107 = t79 * pkin(8) + t114;
t70 = t100 * t80 + t95 * t122;
t106 = t70 * pkin(4) + t108;
t71 = t103 * t78 + t97 * t125;
t72 = -t78 * t100 + t97 * t122;
t105 = -t72 * pkin(4) - t71 * qJ(5) + t107;
t82 = t98 * t100 + t103 * t121;
t69 = -t80 * t103 + t95 * t125;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t97 - rSges(2,2) * t95) + g(2) * (rSges(2,1) * t95 + rSges(2,2) * t97) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t81 - rSges(3,2) * t80 + t126) + g(2) * (rSges(3,1) * t79 - rSges(3,2) * t78 + t91) + g(3) * (t98 * rSges(3,3) + t118) + (rSges(3,3) * t129 + g(3) * (rSges(3,1) * t101 + rSges(3,2) * t104) + (-rSges(3,3) - pkin(7)) * t128) * t96) - m(4) * (g(1) * (-rSges(4,2) * t81 + rSges(4,3) * t80 + t112) + g(2) * (-rSges(4,2) * t79 + rSges(4,3) * t78 + t114) + g(3) * (t98 * rSges(4,1) + t116) + (rSges(4,1) * t129 + g(3) * (-rSges(4,2) * t101 - rSges(4,3) * t104 - t117) + (-rSges(4,1) - pkin(7)) * t128) * t96) - m(5) * (g(1) * (rSges(5,1) * t70 - rSges(5,2) * t69 + (rSges(5,3) + pkin(8)) * t81 + t108) + g(2) * (-rSges(5,1) * t72 + rSges(5,2) * t71 + rSges(5,3) * t79 + t107) + g(3) * (t83 * rSges(5,1) - t82 * rSges(5,2) + t113) + (g(3) * (rSges(5,3) * t101 - t117) + t115) * t96) - m(6) * (g(1) * (-rSges(6,2) * t70 + (rSges(6,1) + pkin(8)) * t81 + t130 * t69 + t106) + g(2) * (rSges(6,1) * t79 + rSges(6,2) * t72 - rSges(6,3) * t71 + t105) + g(3) * (-t83 * rSges(6,2) + t130 * t82 + t111) + (g(3) * (rSges(6,1) * t101 - t117) + t115) * t96) - m(7) * (g(1) * (t109 * t69 + (pkin(8) + t110) * t81 + t131 * t70 + t106) + g(2) * (t79 * pkin(5) + (t102 * t79 - t71 * t99) * rSges(7,1) + (-t102 * t71 - t79 * t99) * rSges(7,2) + t105 - t131 * t72) + t115 * t96 + (t111 + t109 * t82 + (t101 * t110 - t117) * t96 + t131 * t83) * g(3));
U  = t1;
