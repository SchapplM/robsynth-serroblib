% Calculate potential energy for
% S6RRRRRR7
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
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:28:55
% EndTime: 2019-03-10 04:28:55
% DurationCPUTime: 0.41s
% Computational Cost: add. (335->116), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->51)
t99 = sin(qJ(4));
t130 = pkin(4) * t99 + pkin(9);
t106 = -pkin(11) - pkin(10);
t128 = rSges(4,3) + pkin(9);
t117 = cos(pkin(6));
t127 = t117 * pkin(8) + pkin(7);
t126 = pkin(10) + rSges(5,3);
t103 = cos(qJ(4));
t88 = t103 * pkin(4) + pkin(3);
t125 = cos(qJ(3));
t105 = cos(qJ(1));
t102 = sin(qJ(1));
t98 = sin(pkin(6));
t121 = t102 * t98;
t124 = t105 * pkin(1) + pkin(8) * t121;
t123 = pkin(12) - t106 + rSges(7,3);
t101 = sin(qJ(2));
t122 = t101 * t98;
t104 = cos(qJ(2));
t120 = t104 * t98;
t119 = t105 * t98;
t118 = -t106 + rSges(6,3);
t97 = qJ(4) + qJ(5);
t116 = pkin(2) * t122 + t127;
t113 = t102 * t117;
t78 = -t101 * t113 + t105 * t104;
t115 = t78 * pkin(2) + t124;
t114 = t98 * t125;
t112 = t105 * t117;
t76 = t101 * t112 + t102 * t104;
t93 = t102 * pkin(1);
t111 = t76 * pkin(2) - pkin(8) * t119 + t93;
t89 = sin(t97);
t90 = cos(t97);
t110 = t90 * rSges(6,1) - t89 * rSges(6,2) + t88;
t92 = qJ(6) + t97;
t86 = sin(t92);
t87 = cos(t92);
t109 = rSges(7,1) * t87 - rSges(7,2) * t86 + pkin(5) * t90 + t88;
t108 = rSges(7,1) * t86 + rSges(7,2) * t87 + pkin(5) * t89 + t130;
t107 = rSges(6,1) * t89 + rSges(6,2) * t90 + t130;
t100 = sin(qJ(3));
t77 = t105 * t101 + t104 * t113;
t75 = t101 * t102 - t104 * t112;
t74 = t117 * t100 + t101 * t114;
t73 = t100 * t122 - t117 * t125;
t70 = t100 * t121 + t78 * t125;
t69 = t100 * t78 - t102 * t114;
t68 = -t100 * t119 + t76 * t125;
t67 = t76 * t100 + t105 * t114;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t105 - rSges(2,2) * t102) + g(2) * (rSges(2,1) * t102 + rSges(2,2) * t105) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t78 - rSges(3,2) * t77 + t124) + g(2) * (rSges(3,1) * t76 - rSges(3,2) * t75 + t93) + g(3) * (t117 * rSges(3,3) + t127) + (g(1) * rSges(3,3) * t102 + g(3) * (rSges(3,1) * t101 + rSges(3,2) * t104) + g(2) * (-rSges(3,3) - pkin(8)) * t105) * t98) - m(4) * (g(1) * (rSges(4,1) * t70 - rSges(4,2) * t69 + t128 * t77 + t115) + g(2) * (rSges(4,1) * t68 - rSges(4,2) * t67 + t128 * t75 + t111) + g(3) * (rSges(4,1) * t74 - rSges(4,2) * t73 - t128 * t120 + t116)) - m(5) * (g(1) * (t70 * pkin(3) + t77 * pkin(9) + (t103 * t70 + t77 * t99) * rSges(5,1) + (t103 * t77 - t70 * t99) * rSges(5,2) + t126 * t69 + t115) + g(2) * (t68 * pkin(3) + t75 * pkin(9) + (t103 * t68 + t75 * t99) * rSges(5,1) + (t103 * t75 - t68 * t99) * rSges(5,2) + t126 * t67 + t111) + g(3) * (t74 * pkin(3) - pkin(9) * t120 + (t103 * t74 - t99 * t120) * rSges(5,1) + (-t103 * t120 - t74 * t99) * rSges(5,2) + t126 * t73 + t116)) - m(6) * (g(1) * (t107 * t77 + t110 * t70 + t118 * t69 + t115) + g(2) * (t107 * t75 + t110 * t68 + t118 * t67 + t111) + g(3) * (-t107 * t120 + t110 * t74 + t118 * t73 + t116)) - m(7) * (g(1) * (t108 * t77 + t109 * t70 + t123 * t69 + t115) + g(2) * (t108 * t75 + t109 * t68 + t123 * t67 + t111) + g(3) * (-t108 * t120 + t109 * t74 + t123 * t73 + t116));
U  = t1;
