% Calculate potential energy for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR13_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR13_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:51
% EndTime: 2019-03-09 14:44:51
% DurationCPUTime: 0.63s
% Computational Cost: add. (272->123), mult. (536->157), div. (0->0), fcn. (622->12), ass. (0->52)
t130 = -pkin(11) - pkin(10) - rSges(7,3);
t99 = sin(qJ(1));
t129 = g(1) * t99;
t95 = cos(pkin(6));
t128 = t95 * pkin(8) + pkin(7);
t127 = pkin(10) + rSges(6,3);
t103 = cos(qJ(1));
t126 = g(2) * t103;
t98 = sin(qJ(2));
t117 = t103 * t98;
t102 = cos(qJ(2));
t119 = t102 * t99;
t77 = t95 * t117 + t119;
t96 = sin(qJ(5));
t125 = t77 * t96;
t94 = sin(pkin(6));
t124 = t94 * t98;
t123 = t94 * t99;
t122 = t98 * t99;
t121 = t103 * pkin(1) + pkin(8) * t123;
t120 = t102 * t94;
t118 = t103 * t94;
t116 = t102 * qJ(3);
t115 = t102 * t103;
t114 = pkin(2) * t124 + t128;
t113 = t103 * (-pkin(3) - pkin(8));
t112 = g(2) * t113;
t111 = t95 * pkin(3) + pkin(9) * t124 + t114;
t76 = -t95 * t115 + t122;
t91 = t99 * pkin(1);
t110 = t77 * pkin(2) + t76 * qJ(3) + t91;
t78 = t95 * t119 + t117;
t79 = -t95 * t122 + t115;
t109 = t79 * pkin(2) + t78 * qJ(3) + t121;
t100 = cos(qJ(5));
t86 = pkin(5) * t100 + pkin(4);
t93 = qJ(5) + qJ(6);
t87 = sin(t93);
t88 = cos(t93);
t108 = t88 * rSges(7,1) - t87 * rSges(7,2) + t86;
t107 = pkin(3) * t123 + t109;
t106 = t77 * pkin(9) + t110;
t105 = rSges(7,1) * t87 + rSges(7,2) * t88 + pkin(5) * t96;
t101 = cos(qJ(4));
t97 = sin(qJ(4));
t75 = t101 * t95 - t97 * t120;
t74 = t101 * t120 + t95 * t97;
t70 = -t101 * t118 + t76 * t97;
t69 = t101 * t76 + t97 * t118;
t68 = t101 * t123 + t78 * t97;
t67 = -t78 * t101 + t97 * t123;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t103 - rSges(2,2) * t99) + g(2) * (rSges(2,1) * t99 + rSges(2,2) * t103) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t79 - rSges(3,2) * t78 + t121) + g(2) * (rSges(3,1) * t77 - rSges(3,2) * t76 + t91) + g(3) * (rSges(3,3) * t95 + t128) + (rSges(3,3) * t129 + g(3) * (rSges(3,1) * t98 + rSges(3,2) * t102) + (-rSges(3,3) - pkin(8)) * t126) * t94) - m(4) * (g(1) * (-rSges(4,2) * t79 + rSges(4,3) * t78 + t109) + g(2) * (-rSges(4,2) * t77 + rSges(4,3) * t76 + t110) + g(3) * (rSges(4,1) * t95 + t114) + (rSges(4,1) * t129 + g(3) * (-rSges(4,2) * t98 - rSges(4,3) * t102 - t116) + (-rSges(4,1) - pkin(8)) * t126) * t94) - m(5) * (g(1) * (rSges(5,1) * t68 - rSges(5,2) * t67 + (rSges(5,3) + pkin(9)) * t79 + t107) + g(2) * (rSges(5,1) * t70 + rSges(5,2) * t69 + rSges(5,3) * t77 + t106) + g(3) * (rSges(5,1) * t75 - rSges(5,2) * t74 + t111) + (g(3) * (rSges(5,3) * t98 - t116) + t112) * t94) - m(6) * (g(1) * (t68 * pkin(4) + t79 * pkin(9) + (t100 * t68 + t79 * t96) * rSges(6,1) + (t100 * t79 - t68 * t96) * rSges(6,2) + t127 * t67 + t107) + g(2) * (t70 * pkin(4) + (t100 * t70 + t125) * rSges(6,1) + (t100 * t77 - t70 * t96) * rSges(6,2) - t127 * t69 + t94 * t113 + t106) + g(3) * (t75 * pkin(4) - t94 * t116 + (t100 * t75 + t96 * t124) * rSges(6,1) + (t100 * t124 - t75 * t96) * rSges(6,2) + t127 * t74 + t111)) - m(7) * (g(1) * (t108 * t68 + (pkin(9) + t105) * t79 - t130 * t67 + t107) + g(2) * (t70 * t86 + pkin(5) * t125 + (t70 * t88 + t77 * t87) * rSges(7,1) + (-t70 * t87 + t77 * t88) * rSges(7,2) + t106 + t130 * t69) + t112 * t94 + (t111 + t108 * t75 + (t105 * t98 - t116) * t94 - t130 * t74) * g(3));
U  = t1;
