% Calculate potential energy for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR14_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR14_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:28
% EndTime: 2019-03-09 11:32:28
% DurationCPUTime: 0.59s
% Computational Cost: add. (249->115), mult. (514->141), div. (0->0), fcn. (592->10), ass. (0->49)
t132 = pkin(10) + rSges(7,3);
t131 = rSges(6,3) + qJ(5);
t97 = cos(pkin(6));
t130 = t97 * pkin(8) + pkin(7);
t101 = sin(qJ(1));
t129 = g(1) * t101;
t105 = cos(qJ(1));
t128 = g(2) * t105;
t96 = sin(pkin(6));
t125 = t101 * t96;
t127 = t105 * pkin(1) + pkin(8) * t125;
t100 = sin(qJ(2));
t126 = t100 * t96;
t104 = cos(qJ(2));
t124 = t104 * t96;
t123 = t105 * t96;
t122 = qJ(3) * t104;
t121 = t100 * t101;
t120 = t100 * t105;
t119 = t101 * t104;
t118 = t104 * t105;
t117 = pkin(2) * t126 + t130;
t116 = (-pkin(3) - pkin(8)) * t128;
t115 = t97 * pkin(3) + pkin(9) * t126 + t117;
t81 = -t97 * t118 + t121;
t82 = t97 * t120 + t119;
t94 = t101 * pkin(1);
t114 = t82 * pkin(2) + t81 * qJ(3) + t94;
t103 = cos(qJ(4));
t99 = sin(qJ(4));
t80 = t103 * t97 - t99 * t124;
t113 = t80 * pkin(4) + t115;
t83 = t97 * t119 + t120;
t84 = -t97 * t121 + t118;
t112 = t84 * pkin(2) + t83 * qJ(3) + t127;
t102 = cos(qJ(6));
t98 = sin(qJ(6));
t111 = t102 * rSges(7,1) - t98 * rSges(7,2) + pkin(5);
t110 = t98 * rSges(7,1) + t102 * rSges(7,2) + qJ(5);
t109 = pkin(3) * t125 + t112;
t108 = t82 * pkin(9) + t114;
t71 = t103 * t125 + t83 * t99;
t107 = t71 * pkin(4) + t109;
t72 = t81 * t103 + t99 * t123;
t73 = t103 * t123 - t81 * t99;
t106 = -t73 * pkin(4) - t72 * qJ(5) + t108;
t79 = t103 * t124 + t97 * t99;
t70 = -t83 * t103 + t99 * t125;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t105 - t101 * rSges(2,2)) + g(2) * (t101 * rSges(2,1) + rSges(2,2) * t105) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t127) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t94) + g(3) * (rSges(3,3) * t97 + t130) + (rSges(3,3) * t129 + g(3) * (rSges(3,1) * t100 + rSges(3,2) * t104) + (-rSges(3,3) - pkin(8)) * t128) * t96) - m(4) * (g(1) * (-rSges(4,2) * t84 + rSges(4,3) * t83 + t112) + g(2) * (-t82 * rSges(4,2) + t81 * rSges(4,3) + t114) + g(3) * (rSges(4,1) * t97 + t117) + (rSges(4,1) * t129 + g(3) * (-rSges(4,2) * t100 - rSges(4,3) * t104 - t122) + (-rSges(4,1) - pkin(8)) * t128) * t96) - m(5) * (g(1) * (rSges(5,1) * t71 - t70 * rSges(5,2) + (rSges(5,3) + pkin(9)) * t84 + t109) + g(2) * (-t73 * rSges(5,1) + t72 * rSges(5,2) + t82 * rSges(5,3) + t108) + g(3) * (rSges(5,1) * t80 - rSges(5,2) * t79 + t115) + (g(3) * (rSges(5,3) * t100 - t122) + t116) * t96) - m(6) * (g(1) * (-rSges(6,2) * t71 + (rSges(6,1) + pkin(9)) * t84 + t131 * t70 + t107) + g(2) * (t82 * rSges(6,1) + t73 * rSges(6,2) - t72 * rSges(6,3) + t106) + g(3) * (-rSges(6,2) * t80 + t131 * t79 + t113) + (g(3) * (rSges(6,1) * t100 - t122) + t116) * t96) - m(7) * (g(1) * (t110 * t70 + (pkin(9) + t111) * t84 + t132 * t71 + t107) + g(2) * (t82 * pkin(5) + (t102 * t82 - t72 * t98) * rSges(7,1) + (-t102 * t72 - t82 * t98) * rSges(7,2) + t106 - t132 * t73) + t116 * t96 + (t113 + t110 * t79 + (t100 * t111 - t122) * t96 + t132 * t80) * g(3));
U  = t1;
