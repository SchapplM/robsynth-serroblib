% Calculate potential energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:35
% EndTime: 2019-03-10 04:11:35
% DurationCPUTime: 0.44s
% Computational Cost: add. (355->130), mult. (545->170), div. (0->0), fcn. (633->14), ass. (0->52)
t103 = cos(pkin(6));
t134 = t103 * pkin(8) + pkin(7);
t133 = pkin(9) + rSges(4,3);
t132 = pkin(11) + rSges(6,3);
t105 = sin(qJ(3));
t131 = pkin(3) * t105;
t111 = cos(qJ(1));
t102 = sin(pkin(6));
t107 = sin(qJ(1));
t125 = t102 * t107;
t130 = t111 * pkin(1) + pkin(8) * t125;
t104 = sin(qJ(5));
t110 = cos(qJ(2));
t119 = t110 * t111;
t106 = sin(qJ(2));
t122 = t106 * t107;
t79 = -t103 * t119 + t122;
t129 = t104 * t79;
t120 = t107 * t110;
t121 = t106 * t111;
t81 = t103 * t120 + t121;
t128 = t104 * t81;
t127 = pkin(12) + pkin(11) + rSges(7,3);
t126 = t102 * t106;
t124 = t102 * t110;
t123 = t102 * t111;
t118 = t104 * t124;
t117 = t105 * t125;
t113 = -pkin(10) - pkin(9);
t109 = cos(qJ(3));
t92 = pkin(3) * t109 + pkin(2);
t116 = t103 * t131 + t113 * t124 + t92 * t126 + t134;
t82 = -t103 * t122 + t119;
t115 = pkin(3) * t117 - t81 * t113 + t82 * t92 + t130;
t80 = t103 * t121 + t120;
t98 = t107 * pkin(1);
t114 = t80 * t92 + t98 + (-pkin(8) - t131) * t123 - t79 * t113;
t108 = cos(qJ(5));
t101 = qJ(3) + qJ(4);
t100 = qJ(5) + qJ(6);
t96 = cos(t101);
t95 = cos(t100);
t94 = sin(t101);
t93 = sin(t100);
t91 = pkin(5) * t108 + pkin(4);
t76 = t103 * t94 + t126 * t96;
t75 = -t103 * t96 + t126 * t94;
t72 = t125 * t94 + t82 * t96;
t71 = -t125 * t96 + t82 * t94;
t70 = -t123 * t94 + t80 * t96;
t69 = t123 * t96 + t80 * t94;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t111 - rSges(2,2) * t107) + g(2) * (rSges(2,1) * t107 + rSges(2,2) * t111) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t82 - rSges(3,2) * t81 + t130) + g(2) * (rSges(3,1) * t80 - rSges(3,2) * t79 + t98) + g(3) * (rSges(3,3) * t103 + t134) + (g(1) * rSges(3,3) * t107 + g(3) * (rSges(3,1) * t106 + rSges(3,2) * t110) + g(2) * (-rSges(3,3) - pkin(8)) * t111) * t102) - m(4) * (g(1) * (t82 * pkin(2) + (t109 * t82 + t117) * rSges(4,1) + (-t105 * t82 + t109 * t125) * rSges(4,2) + t133 * t81 + t130) + g(2) * (t80 * pkin(2) + t98 - pkin(8) * t123 + (-t105 * t123 + t109 * t80) * rSges(4,1) + (-t105 * t80 - t109 * t123) * rSges(4,2) + t133 * t79) + g(3) * ((t105 * rSges(4,1) + t109 * rSges(4,2)) * t103 + (-t133 * t110 + (rSges(4,1) * t109 - rSges(4,2) * t105 + pkin(2)) * t106) * t102 + t134)) - m(5) * (g(1) * (rSges(5,1) * t72 - rSges(5,2) * t71 + rSges(5,3) * t81 + t115) + g(2) * (rSges(5,1) * t70 - rSges(5,2) * t69 + rSges(5,3) * t79 + t114) + g(3) * (rSges(5,1) * t76 - rSges(5,2) * t75 - rSges(5,3) * t124 + t116)) - m(6) * (g(1) * (t72 * pkin(4) + (t108 * t72 + t128) * rSges(6,1) + (-t104 * t72 + t108 * t81) * rSges(6,2) + t132 * t71 + t115) + g(2) * (t70 * pkin(4) + (t108 * t70 + t129) * rSges(6,1) + (-t104 * t70 + t108 * t79) * rSges(6,2) + t132 * t69 + t114) + g(3) * (t76 * pkin(4) + (t108 * t76 - t118) * rSges(6,1) + (-t104 * t76 - t108 * t124) * rSges(6,2) + t132 * t75 + t116)) - m(7) * (g(1) * (t72 * t91 + pkin(5) * t128 + (t72 * t95 + t81 * t93) * rSges(7,1) + (-t72 * t93 + t81 * t95) * rSges(7,2) + t127 * t71 + t115) + g(2) * (t70 * t91 + pkin(5) * t129 + (t70 * t95 + t79 * t93) * rSges(7,1) + (-t70 * t93 + t79 * t95) * rSges(7,2) + t127 * t69 + t114) + g(3) * (t76 * t91 - pkin(5) * t118 + (-t124 * t93 + t76 * t95) * rSges(7,1) + (-t124 * t95 - t76 * t93) * rSges(7,2) + t127 * t75 + t116));
U  = t1;
