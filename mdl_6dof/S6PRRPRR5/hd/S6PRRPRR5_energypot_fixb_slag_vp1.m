% Calculate potential energy for
% S6PRRPRR5
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:12
% EndTime: 2019-03-08 22:16:12
% DurationCPUTime: 0.40s
% Computational Cost: add. (335->116), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->51)
t98 = sin(pkin(12));
t130 = pkin(4) * t98 + pkin(8);
t128 = rSges(4,3) + pkin(8);
t101 = cos(pkin(12));
t88 = t101 * pkin(4) + pkin(3);
t127 = cos(qJ(3));
t100 = sin(pkin(6));
t126 = pkin(7) * t100;
t103 = -pkin(9) - qJ(4);
t102 = cos(pkin(11));
t99 = sin(pkin(11));
t125 = t102 * pkin(1) + t99 * t126;
t124 = pkin(10) - t103 + rSges(7,3);
t120 = cos(pkin(6));
t123 = t120 * pkin(7) + qJ(1);
t122 = qJ(4) + rSges(5,3);
t121 = -t103 + rSges(6,3);
t104 = sin(qJ(3));
t119 = t100 * t104;
t105 = sin(qJ(2));
t118 = t100 * t105;
t106 = cos(qJ(2));
t117 = t100 * t106;
t97 = pkin(12) + qJ(5);
t113 = t105 * t120;
t76 = t102 * t106 - t99 * t113;
t116 = t76 * pkin(2) + t125;
t115 = pkin(2) * t118 + t123;
t114 = t100 * t127;
t112 = t106 * t120;
t74 = t102 * t113 + t99 * t106;
t92 = t99 * pkin(1);
t111 = t74 * pkin(2) - t102 * t126 + t92;
t89 = sin(t97);
t90 = cos(t97);
t110 = t90 * rSges(6,1) - t89 * rSges(6,2) + t88;
t91 = qJ(6) + t97;
t86 = sin(t91);
t87 = cos(t91);
t109 = rSges(7,1) * t87 - rSges(7,2) * t86 + pkin(5) * t90 + t88;
t108 = rSges(7,1) * t86 + rSges(7,2) * t87 + pkin(5) * t89 + t130;
t107 = rSges(6,1) * t89 + rSges(6,2) * t90 + t130;
t78 = t120 * t104 + t105 * t114;
t77 = t104 * t118 - t120 * t127;
t75 = t102 * t105 + t99 * t112;
t73 = -t102 * t112 + t105 * t99;
t70 = t99 * t119 + t76 * t127;
t69 = t104 * t76 - t99 * t114;
t68 = -t102 * t119 + t74 * t127;
t67 = t102 * t114 + t74 * t104;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t102 - rSges(2,2) * t99) + g(2) * (rSges(2,1) * t99 + rSges(2,2) * t102) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t76 - rSges(3,2) * t75 + t125) + g(2) * (rSges(3,1) * t74 - rSges(3,2) * t73 + t92) + g(3) * (t120 * rSges(3,3) + t123) + (g(1) * rSges(3,3) * t99 + g(3) * (rSges(3,1) * t105 + rSges(3,2) * t106) + g(2) * (-rSges(3,3) - pkin(7)) * t102) * t100) - m(4) * (g(1) * (rSges(4,1) * t70 - rSges(4,2) * t69 + t128 * t75 + t116) + g(2) * (rSges(4,1) * t68 - rSges(4,2) * t67 + t128 * t73 + t111) + g(3) * (t78 * rSges(4,1) - t77 * rSges(4,2) - t128 * t117 + t115)) - m(5) * (g(1) * (t70 * pkin(3) + t75 * pkin(8) + (t101 * t70 + t75 * t98) * rSges(5,1) + (t101 * t75 - t70 * t98) * rSges(5,2) + t122 * t69 + t116) + g(2) * (t68 * pkin(3) + t73 * pkin(8) + (t101 * t68 + t73 * t98) * rSges(5,1) + (t101 * t73 - t68 * t98) * rSges(5,2) + t122 * t67 + t111) + g(3) * (t78 * pkin(3) - pkin(8) * t117 + (t78 * t101 - t98 * t117) * rSges(5,1) + (-t101 * t117 - t78 * t98) * rSges(5,2) + t122 * t77 + t115)) - m(6) * (g(1) * (t107 * t75 + t110 * t70 + t121 * t69 + t116) + g(2) * (t107 * t73 + t110 * t68 + t121 * t67 + t111) + g(3) * (-t107 * t117 + t110 * t78 + t121 * t77 + t115)) - m(7) * (g(1) * (t108 * t75 + t109 * t70 + t124 * t69 + t116) + g(2) * (t108 * t73 + t109 * t68 + t124 * t67 + t111) + g(3) * (-t108 * t117 + t109 * t78 + t124 * t77 + t115));
U  = t1;
