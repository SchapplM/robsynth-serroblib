% Calculate potential energy for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:09
% EndTime: 2019-03-09 11:17:10
% DurationCPUTime: 0.76s
% Computational Cost: add. (284->121), mult. (482->150), div. (0->0), fcn. (547->12), ass. (0->53)
t136 = pkin(9) + rSges(5,3);
t135 = pkin(10) + rSges(7,3);
t101 = cos(pkin(6));
t134 = t101 * pkin(8) + pkin(7);
t104 = sin(qJ(4));
t133 = pkin(4) * t104;
t106 = sin(qJ(1));
t132 = g(1) * t106;
t110 = cos(qJ(1));
t131 = g(2) * t110;
t109 = cos(qJ(2));
t123 = t106 * t109;
t105 = sin(qJ(2));
t124 = t105 * t110;
t84 = t101 * t124 + t123;
t97 = t106 * pkin(1);
t130 = t84 * pkin(2) + t97;
t100 = sin(pkin(6));
t128 = t100 * t106;
t129 = t110 * pkin(1) + pkin(8) * t128;
t127 = t100 * t109;
t126 = t100 * t110;
t125 = t105 * t106;
t122 = t109 * t110;
t121 = t100 * t105 * pkin(2) + t134;
t86 = -t101 * t125 + t122;
t120 = t86 * pkin(2) + t129;
t108 = cos(qJ(4));
t93 = pkin(4) * t108 + pkin(3);
t119 = t101 * t93 + t121;
t118 = (-pkin(8) - t93) * t131;
t83 = -t101 * t122 + t125;
t117 = t83 * qJ(3) + t130;
t116 = (-qJ(3) - t133) * t109;
t85 = t101 * t123 + t124;
t115 = t85 * qJ(3) + t120;
t114 = rSges(5,1) * t108 - rSges(5,2) * t104 + pkin(3);
t113 = rSges(5,1) * t104 + rSges(5,2) * t108 + qJ(3);
t102 = -qJ(5) - pkin(9);
t112 = -t84 * t102 + t83 * t133 + t117;
t111 = -t86 * t102 + t93 * t128 + t85 * t133 + t115;
t107 = cos(qJ(6));
t103 = sin(qJ(6));
t99 = qJ(4) + pkin(11);
t95 = cos(t99);
t94 = sin(t99);
t76 = t101 * t95 - t127 * t94;
t75 = t101 * t94 + t127 * t95;
t73 = -t126 * t95 + t83 * t94;
t72 = t126 * t94 + t83 * t95;
t71 = t128 * t95 + t85 * t94;
t70 = t128 * t94 - t85 * t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t110 - t106 * rSges(2,2)) + g(2) * (t106 * rSges(2,1) + rSges(2,2) * t110) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t86 - rSges(3,2) * t85 + t129) + g(2) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t97) + g(3) * (rSges(3,3) * t101 + t134) + (rSges(3,3) * t132 + g(3) * (rSges(3,1) * t105 + rSges(3,2) * t109) + (-rSges(3,3) - pkin(8)) * t131) * t100) - m(4) * (g(1) * (-rSges(4,2) * t86 + rSges(4,3) * t85 + t115) + g(2) * (-t84 * rSges(4,2) + t83 * rSges(4,3) + t117) + g(3) * (rSges(4,1) * t101 + t121) + (rSges(4,1) * t132 + g(3) * (-rSges(4,2) * t105 + (-rSges(4,3) - qJ(3)) * t109) + (-rSges(4,1) - pkin(8)) * t131) * t100) - m(5) * ((t114 * t132 + (-pkin(8) - t114) * t131) * t100 + (t121 + t114 * t101 + (t136 * t105 - t113 * t109) * t100) * g(3) + (t113 * t83 + t136 * t84 + t130) * g(2) + (t113 * t85 + t136 * t86 + t120) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t71 - rSges(6,2) * t70 + rSges(6,3) * t86 + t111) + g(2) * (t73 * rSges(6,1) + t72 * rSges(6,2) + t84 * rSges(6,3) + t112) + g(3) * (rSges(6,1) * t76 - rSges(6,2) * t75 + t119) + (g(3) * (t116 + (rSges(6,3) - t102) * t105) + t118) * t100) - m(7) * (g(1) * (t71 * pkin(5) + (t103 * t86 + t107 * t71) * rSges(7,1) + (-t103 * t71 + t107 * t86) * rSges(7,2) + t135 * t70 + t111) + g(2) * (t73 * pkin(5) + (t103 * t84 + t107 * t73) * rSges(7,1) + (-t103 * t73 + t107 * t84) * rSges(7,2) + t112 - t135 * t72) + t118 * t100 + (t135 * t75 + t119 + (t107 * rSges(7,1) - t103 * rSges(7,2) + pkin(5)) * t76 + (t116 + (t103 * rSges(7,1) + t107 * rSges(7,2) - t102) * t105) * t100) * g(3));
U  = t1;
