% Calculate potential energy for
% S6RRPPRR9
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:14
% EndTime: 2019-03-09 09:28:14
% DurationCPUTime: 0.59s
% Computational Cost: add. (228->120), mult. (459->148), div. (0->0), fcn. (516->10), ass. (0->52)
t136 = pkin(10) + rSges(7,3);
t138 = -pkin(3) - pkin(8);
t100 = cos(pkin(6));
t137 = t100 * pkin(8) + pkin(7);
t104 = sin(qJ(1));
t135 = g(1) * t104;
t108 = cos(qJ(1));
t134 = g(2) * t108;
t107 = cos(qJ(2));
t120 = t104 * t107;
t103 = sin(qJ(2));
t121 = t103 * t108;
t82 = t100 * t121 + t120;
t97 = t104 * pkin(1);
t133 = t82 * pkin(2) + t97;
t99 = sin(pkin(6));
t130 = t104 * t99;
t132 = t108 * pkin(1) + pkin(8) * t130;
t131 = t103 * t99;
t106 = cos(qJ(5));
t129 = t106 * t99;
t128 = t107 * t99;
t127 = t108 * t99;
t119 = t107 * t108;
t122 = t103 * t104;
t81 = -t100 * t119 + t122;
t126 = t81 * qJ(3);
t83 = t100 * t120 + t121;
t125 = t83 * qJ(3);
t124 = rSges(6,3) - qJ(3);
t123 = qJ(3) * t107;
t118 = pkin(2) * t131 + t137;
t117 = -pkin(9) - t124;
t84 = -t100 * t122 + t119;
t116 = t84 * pkin(2) + t132;
t115 = t100 * pkin(3) + qJ(4) * t131 + t118;
t114 = t82 * qJ(4) + t133;
t113 = (-pkin(4) + t138) * t134;
t112 = t100 * pkin(4) + pkin(9) * t128 + t115;
t111 = pkin(3) * t130 + t84 * qJ(4) + t116;
t110 = t114 + t126;
t109 = pkin(4) * t130 + t111;
t105 = cos(qJ(6));
t102 = sin(qJ(5));
t101 = sin(qJ(6));
t80 = t100 * t106 + t102 * t131;
t79 = t100 * t102 - t103 * t129;
t75 = t82 * t102 - t106 * t127;
t74 = t102 * t127 + t82 * t106;
t73 = t102 * t84 + t104 * t129;
t72 = t102 * t130 - t84 * t106;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t108 - t104 * rSges(2,2)) + g(2) * (t104 * rSges(2,1) + rSges(2,2) * t108) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t132) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t97) + g(3) * (rSges(3,3) * t100 + t137) + (rSges(3,3) * t135 + g(3) * (rSges(3,1) * t103 + rSges(3,2) * t107) + (-rSges(3,3) - pkin(8)) * t134) * t99) - m(4) * (g(1) * (-rSges(4,2) * t84 + rSges(4,3) * t83 + t116 + t125) + g(2) * (-t82 * rSges(4,2) + t81 * rSges(4,3) + t126 + t133) + g(3) * (rSges(4,1) * t100 + t118) + (rSges(4,1) * t135 + g(3) * (-rSges(4,2) * t103 - rSges(4,3) * t107 - t123) + (-rSges(4,1) - pkin(8)) * t134) * t99) - m(5) * (g(1) * (rSges(5,2) * t83 + rSges(5,3) * t84 + t111 + t125) + g(2) * (t81 * rSges(5,2) + t82 * rSges(5,3) + t110) + g(3) * (rSges(5,1) * t100 + t115) + (rSges(5,1) * t135 + g(3) * (-rSges(5,2) * t107 + rSges(5,3) * t103 - t123) + (-rSges(5,1) + t138) * t134) * t99) - m(6) * (g(3) * (rSges(6,1) * t80 - rSges(6,2) * t79 + t112) + (g(3) * t124 * t107 + t113) * t99 + (t75 * rSges(6,1) + t74 * rSges(6,2) + t117 * t81 + t114) * g(2) + (rSges(6,1) * t73 - rSges(6,2) * t72 + t117 * t83 + t109) * g(1)) - m(7) * (g(1) * ((t105 * rSges(7,1) - t101 * rSges(7,2) + pkin(5)) * t73 + (-t101 * rSges(7,1) - t105 * rSges(7,2) - pkin(9) + qJ(3)) * t83 + t136 * t72 + t109) + g(2) * (t75 * pkin(5) - t81 * pkin(9) + (-t101 * t81 + t105 * t75) * rSges(7,1) + (-t101 * t75 - t105 * t81) * rSges(7,2) + t110 - t136 * t74) + g(3) * (t80 * pkin(5) - t99 * t123 + (t101 * t128 + t105 * t80) * rSges(7,1) + (-t101 * t80 + t105 * t128) * rSges(7,2) + t136 * t79 + t112) + t99 * t113);
U  = t1;
