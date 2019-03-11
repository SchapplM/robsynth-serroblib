% Calculate potential energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:15
% EndTime: 2019-03-09 00:06:15
% DurationCPUTime: 0.40s
% Computational Cost: add. (323->125), mult. (604->162), div. (0->0), fcn. (717->12), ass. (0->52)
t104 = sin(qJ(4));
t114 = t104 * pkin(4) + pkin(8);
t109 = -pkin(10) - pkin(9);
t130 = rSges(4,3) + pkin(8);
t100 = qJ(4) + qJ(5);
t93 = sin(t100);
t129 = pkin(5) * t93 + t114;
t128 = pkin(9) + rSges(5,3);
t107 = cos(qJ(4));
t92 = t107 * pkin(4) + pkin(3);
t127 = cos(qJ(3));
t102 = sin(pkin(6));
t126 = pkin(7) * t102;
t124 = rSges(7,3) + qJ(6) - t109;
t101 = sin(pkin(11));
t103 = cos(pkin(11));
t123 = t103 * pkin(1) + t101 * t126;
t122 = rSges(6,3) - t109;
t120 = cos(pkin(6));
t121 = t120 * pkin(7) + qJ(1);
t105 = sin(qJ(3));
t119 = t102 * t105;
t106 = sin(qJ(2));
t118 = t102 * t106;
t108 = cos(qJ(2));
t117 = t102 * t108;
t112 = t106 * t120;
t82 = -t101 * t112 + t103 * t108;
t116 = t82 * pkin(2) + t123;
t115 = pkin(2) * t118 + t121;
t113 = t102 * t127;
t111 = t108 * t120;
t80 = t101 * t108 + t103 * t112;
t95 = t101 * pkin(1);
t110 = t80 * pkin(2) - t103 * t126 + t95;
t94 = cos(t100);
t85 = pkin(5) * t94 + t92;
t84 = t120 * t105 + t106 * t113;
t83 = t105 * t118 - t120 * t127;
t81 = t101 * t111 + t103 * t106;
t79 = t101 * t106 - t103 * t111;
t76 = t101 * t119 + t82 * t127;
t75 = -t101 * t113 + t82 * t105;
t74 = -t103 * t119 + t80 * t127;
t73 = t103 * t113 + t80 * t105;
t72 = -t93 * t117 + t84 * t94;
t71 = -t94 * t117 - t84 * t93;
t70 = t76 * t94 + t81 * t93;
t69 = -t76 * t93 + t81 * t94;
t68 = t74 * t94 + t79 * t93;
t67 = -t74 * t93 + t79 * t94;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t103 * rSges(2,1) - t101 * rSges(2,2)) + g(2) * (t101 * rSges(2,1) + t103 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t123) + g(2) * (t80 * rSges(3,1) - t79 * rSges(3,2) + t95) + g(3) * (t120 * rSges(3,3) + t121) + (g(1) * rSges(3,3) * t101 + g(3) * (rSges(3,1) * t106 + rSges(3,2) * t108) + g(2) * (-rSges(3,3) - pkin(7)) * t103) * t102) - m(4) * (g(1) * (t76 * rSges(4,1) - t75 * rSges(4,2) + t130 * t81 + t116) + g(2) * (t74 * rSges(4,1) - t73 * rSges(4,2) + t130 * t79 + t110) + g(3) * (t84 * rSges(4,1) - t83 * rSges(4,2) - t130 * t117 + t115)) - m(5) * (g(1) * (t76 * pkin(3) + t81 * pkin(8) + (t81 * t104 + t76 * t107) * rSges(5,1) + (-t76 * t104 + t81 * t107) * rSges(5,2) + t128 * t75 + t116) + g(2) * (t74 * pkin(3) + t79 * pkin(8) + (t79 * t104 + t74 * t107) * rSges(5,1) + (-t74 * t104 + t79 * t107) * rSges(5,2) + t128 * t73 + t110) + g(3) * (t84 * pkin(3) - pkin(8) * t117 + (-t104 * t117 + t84 * t107) * rSges(5,1) + (-t84 * t104 - t107 * t117) * rSges(5,2) + t128 * t83 + t115)) - m(6) * (g(1) * (t70 * rSges(6,1) + t69 * rSges(6,2) + t114 * t81 + t122 * t75 + t76 * t92 + t116) + g(2) * (t68 * rSges(6,1) + t67 * rSges(6,2) + t114 * t79 + t122 * t73 + t74 * t92 + t110) + g(3) * (t72 * rSges(6,1) + t71 * rSges(6,2) - t114 * t117 + t122 * t83 + t84 * t92 + t115)) - m(7) * (g(1) * (t70 * rSges(7,1) + t69 * rSges(7,2) + t124 * t75 + t129 * t81 + t76 * t85 + t116) + g(2) * (t68 * rSges(7,1) + t67 * rSges(7,2) + t124 * t73 + t129 * t79 + t74 * t85 + t110) + g(3) * (t72 * rSges(7,1) + t71 * rSges(7,2) - t129 * t117 + t124 * t83 + t84 * t85 + t115));
U  = t1;
