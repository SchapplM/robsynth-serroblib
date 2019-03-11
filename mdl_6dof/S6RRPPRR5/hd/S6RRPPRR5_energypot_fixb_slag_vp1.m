% Calculate potential energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:07
% EndTime: 2019-03-09 09:07:08
% DurationCPUTime: 0.52s
% Computational Cost: add. (232->116), mult. (469->141), div. (0->0), fcn. (530->10), ass. (0->52)
t106 = cos(qJ(1));
t97 = sin(pkin(6));
t124 = t106 * t97;
t102 = sin(qJ(1));
t105 = cos(qJ(2));
t119 = t102 * t105;
t101 = sin(qJ(2));
t120 = t101 * t106;
t98 = cos(pkin(6));
t83 = t98 * t120 + t119;
t135 = t83 * pkin(3) + qJ(4) * t124;
t134 = t98 * pkin(8) + pkin(7);
t133 = pkin(10) + rSges(7,3);
t132 = g(1) * t102;
t131 = g(2) * t106;
t95 = t102 * pkin(1);
t130 = t83 * pkin(2) + t95;
t127 = t102 * t97;
t129 = t106 * pkin(1) + pkin(8) * t127;
t128 = t101 * t97;
t104 = cos(qJ(5));
t126 = t104 * t97;
t125 = t105 * t97;
t123 = rSges(6,3) - qJ(3);
t122 = qJ(3) * t105;
t121 = t101 * t102;
t118 = t105 * t106;
t117 = pkin(2) * t128 + t134;
t116 = -pkin(9) - t123;
t85 = -t98 * t121 + t118;
t115 = t85 * pkin(2) + t129;
t82 = -t98 * t118 + t121;
t114 = t82 * qJ(3) + t130;
t84 = t98 * t119 + t120;
t113 = t84 * qJ(3) + t115;
t103 = cos(qJ(6));
t99 = sin(qJ(6));
t112 = rSges(7,1) * t103 - rSges(7,2) * t99 + pkin(5);
t111 = pkin(3) * t128 - t98 * qJ(4) + t117;
t110 = -t99 * rSges(7,1) - t103 * rSges(7,2) - pkin(9) + qJ(3);
t109 = t83 * pkin(4) - pkin(8) * t124 + t130 + t135;
t108 = pkin(4) * t128 + pkin(9) * t125 + t111;
t78 = t85 * pkin(3);
t107 = t85 * pkin(4) - qJ(4) * t127 + t115 + t78;
t100 = sin(qJ(5));
t81 = -t100 * t98 + t101 * t126;
t80 = t100 * t128 + t104 * t98;
t73 = -t100 * t127 + t104 * t85;
t72 = t100 * t85 + t102 * t126;
t71 = t100 * t124 + t83 * t104;
t70 = t100 * t83 - t104 * t124;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t106 - t102 * rSges(2,2)) + g(2) * (t102 * rSges(2,1) + t106 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t85 - rSges(3,2) * t84 + t129) + g(2) * (t83 * rSges(3,1) - t82 * rSges(3,2) + t95) + g(3) * (rSges(3,3) * t98 + t134) + (rSges(3,3) * t132 + g(3) * (rSges(3,1) * t101 + rSges(3,2) * t105) + (-rSges(3,3) - pkin(8)) * t131) * t97) - m(4) * (g(1) * (rSges(4,1) * t85 + rSges(4,3) * t84 + t113) + g(2) * (t83 * rSges(4,1) + t82 * rSges(4,3) + t114) + g(3) * (rSges(4,2) * t98 + t117) + (rSges(4,2) * t132 + g(3) * (rSges(4,1) * t101 - rSges(4,3) * t105 - t122) + (-rSges(4,2) - pkin(8)) * t131) * t97) - m(5) * (g(1) * (rSges(5,1) * t85 + rSges(5,2) * t84 + t113 + t78) + g(2) * (t83 * rSges(5,1) + t82 * rSges(5,2) + t114 + t135) + g(3) * (-rSges(5,3) * t98 + t111) + (g(3) * (rSges(5,1) * t101 + (-rSges(5,2) - qJ(3)) * t105) + (rSges(5,3) - pkin(8)) * t131 + (-rSges(5,3) - qJ(4)) * t132) * t97) - m(6) * (g(3) * (rSges(6,1) * t81 - rSges(6,2) * t80 + t123 * t125 + t108) + (t71 * rSges(6,1) - t70 * rSges(6,2) + t116 * t82 + t109) * g(2) + (rSges(6,1) * t73 - rSges(6,2) * t72 + t116 * t84 + t107) * g(1)) - m(7) * (g(1) * (t110 * t84 + t112 * t73 + t133 * t72 + t107) + g(2) * (t110 * t82 + t112 * t71 + t133 * t70 + t109) + g(3) * (t81 * pkin(5) - t97 * t122 + (t103 * t81 + t99 * t125) * rSges(7,1) + (t103 * t125 - t81 * t99) * rSges(7,2) + t133 * t80 + t108));
U  = t1;
