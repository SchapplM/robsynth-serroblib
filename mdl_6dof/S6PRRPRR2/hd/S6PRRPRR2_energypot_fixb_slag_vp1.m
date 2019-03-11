% Calculate potential energy for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:55:47
% EndTime: 2019-03-08 21:55:48
% DurationCPUTime: 0.44s
% Computational Cost: add. (355->130), mult. (545->174), div. (0->0), fcn. (633->14), ass. (0->52)
t134 = pkin(8) + rSges(4,3);
t133 = pkin(9) + rSges(6,3);
t108 = sin(qJ(3));
t132 = pkin(3) * t108;
t104 = cos(pkin(11));
t102 = sin(pkin(11));
t103 = sin(pkin(6));
t126 = t102 * t103;
t131 = t104 * pkin(1) + pkin(7) * t126;
t107 = sin(qJ(5));
t109 = sin(qJ(2));
t105 = cos(pkin(6));
t112 = cos(qJ(2));
t119 = t105 * t112;
t79 = t102 * t109 - t104 * t119;
t130 = t107 * t79;
t81 = t102 * t119 + t104 * t109;
t129 = t107 * t81;
t128 = t105 * pkin(7) + qJ(1);
t127 = pkin(10) + pkin(9) + rSges(7,3);
t125 = t103 * t104;
t124 = t103 * t108;
t123 = t103 * t109;
t111 = cos(qJ(3));
t122 = t103 * t111;
t121 = t103 * t112;
t120 = t105 * t109;
t118 = t102 * t124;
t117 = t107 * t121;
t106 = -qJ(4) - pkin(8);
t82 = -t102 * t120 + t104 * t112;
t92 = pkin(3) * t111 + pkin(2);
t116 = pkin(3) * t118 - t81 * t106 + t82 * t92 + t131;
t115 = t105 * t132 + t106 * t121 + t92 * t123 + t128;
t80 = t102 * t112 + t104 * t120;
t97 = t102 * pkin(1);
t114 = t80 * t92 + t97 + (-pkin(7) - t132) * t125 - t79 * t106;
t110 = cos(qJ(5));
t101 = qJ(5) + qJ(6);
t100 = qJ(3) + pkin(12);
t96 = cos(t101);
t95 = sin(t101);
t94 = cos(t100);
t93 = sin(t100);
t91 = pkin(5) * t110 + pkin(4);
t76 = t105 * t93 + t123 * t94;
t75 = -t105 * t94 + t123 * t93;
t72 = t126 * t93 + t82 * t94;
t71 = -t126 * t94 + t82 * t93;
t70 = -t125 * t93 + t80 * t94;
t69 = t125 * t94 + t80 * t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t104 - rSges(2,2) * t102) + g(2) * (rSges(2,1) * t102 + rSges(2,2) * t104) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t82 - rSges(3,2) * t81 + t131) + g(2) * (rSges(3,1) * t80 - rSges(3,2) * t79 + t97) + g(3) * (rSges(3,3) * t105 + t128) + (g(1) * rSges(3,3) * t102 + g(3) * (rSges(3,1) * t109 + rSges(3,2) * t112) + g(2) * (-rSges(3,3) - pkin(7)) * t104) * t103) - m(4) * (g(1) * (t82 * pkin(2) + (t111 * t82 + t118) * rSges(4,1) + (t102 * t122 - t108 * t82) * rSges(4,2) + t134 * t81 + t131) + g(2) * (t80 * pkin(2) + t97 - pkin(7) * t125 + (-t104 * t124 + t111 * t80) * rSges(4,1) + (-t104 * t122 - t108 * t80) * rSges(4,2) + t134 * t79) + g(3) * ((t108 * rSges(4,1) + t111 * rSges(4,2)) * t105 + (-t134 * t112 + (t111 * rSges(4,1) - t108 * rSges(4,2) + pkin(2)) * t109) * t103 + t128)) - m(5) * (g(1) * (rSges(5,1) * t72 - rSges(5,2) * t71 + rSges(5,3) * t81 + t116) + g(2) * (rSges(5,1) * t70 - rSges(5,2) * t69 + rSges(5,3) * t79 + t114) + g(3) * (rSges(5,1) * t76 - rSges(5,2) * t75 - rSges(5,3) * t121 + t115)) - m(6) * (g(1) * (t72 * pkin(4) + (t110 * t72 + t129) * rSges(6,1) + (-t107 * t72 + t110 * t81) * rSges(6,2) + t133 * t71 + t116) + g(2) * (t70 * pkin(4) + (t110 * t70 + t130) * rSges(6,1) + (-t107 * t70 + t110 * t79) * rSges(6,2) + t133 * t69 + t114) + g(3) * (t76 * pkin(4) + (t110 * t76 - t117) * rSges(6,1) + (-t107 * t76 - t110 * t121) * rSges(6,2) + t133 * t75 + t115)) - m(7) * (g(1) * (t72 * t91 + pkin(5) * t129 + (t72 * t96 + t81 * t95) * rSges(7,1) + (-t72 * t95 + t81 * t96) * rSges(7,2) + t127 * t71 + t116) + g(2) * (t70 * t91 + pkin(5) * t130 + (t70 * t96 + t79 * t95) * rSges(7,1) + (-t70 * t95 + t79 * t96) * rSges(7,2) + t127 * t69 + t114) + g(3) * (t76 * t91 - pkin(5) * t117 + (-t121 * t95 + t76 * t96) * rSges(7,1) + (-t121 * t96 - t76 * t95) * rSges(7,2) + t127 * t75 + t115));
U  = t1;
