% Calculate potential energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:46
% EndTime: 2019-03-09 00:47:46
% DurationCPUTime: 0.40s
% Computational Cost: add. (335->116), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->51)
t101 = sin(qJ(4));
t130 = t101 * pkin(4) + pkin(8);
t106 = -pkin(10) - pkin(9);
t99 = sin(pkin(6));
t129 = pkin(7) * t99;
t128 = rSges(4,3) + pkin(8);
t127 = pkin(9) + rSges(5,3);
t104 = cos(qJ(4));
t88 = t104 * pkin(4) + pkin(3);
t126 = cos(qJ(3));
t100 = cos(pkin(12));
t98 = sin(pkin(12));
t124 = t100 * pkin(1) + t98 * t129;
t123 = pkin(11) - t106 + rSges(7,3);
t102 = sin(qJ(3));
t122 = t102 * t99;
t103 = sin(qJ(2));
t121 = t103 * t99;
t105 = cos(qJ(2));
t120 = t105 * t99;
t117 = cos(pkin(6));
t119 = t117 * pkin(7) + qJ(1);
t118 = -t106 + rSges(6,3);
t97 = qJ(4) + qJ(5);
t113 = t103 * t117;
t76 = t100 * t105 - t98 * t113;
t116 = t76 * pkin(2) + t124;
t115 = pkin(2) * t121 + t119;
t114 = t99 * t126;
t112 = t105 * t117;
t74 = t100 * t113 + t98 * t105;
t91 = t98 * pkin(1);
t111 = t74 * pkin(2) - t100 * t129 + t91;
t89 = sin(t97);
t90 = cos(t97);
t110 = t90 * rSges(6,1) - t89 * rSges(6,2) + t88;
t94 = qJ(6) + t97;
t86 = sin(t94);
t87 = cos(t94);
t109 = t87 * rSges(7,1) - t86 * rSges(7,2) + pkin(5) * t90 + t88;
t108 = t86 * rSges(7,1) + t87 * rSges(7,2) + pkin(5) * t89 + t130;
t107 = t89 * rSges(6,1) + t90 * rSges(6,2) + t130;
t78 = t117 * t102 + t103 * t114;
t77 = t102 * t121 - t117 * t126;
t75 = t100 * t103 + t98 * t112;
t73 = -t100 * t112 + t103 * t98;
t70 = t98 * t122 + t76 * t126;
t69 = t102 * t76 - t98 * t114;
t68 = -t100 * t122 + t74 * t126;
t67 = t100 * t114 + t74 * t102;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t100 - rSges(2,2) * t98) + g(2) * (rSges(2,1) * t98 + rSges(2,2) * t100) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t76 - rSges(3,2) * t75 + t124) + g(2) * (rSges(3,1) * t74 - rSges(3,2) * t73 + t91) + g(3) * (t117 * rSges(3,3) + t119) + (g(1) * rSges(3,3) * t98 + g(3) * (rSges(3,1) * t103 + rSges(3,2) * t105) + g(2) * (-rSges(3,3) - pkin(7)) * t100) * t99) - m(4) * (g(1) * (rSges(4,1) * t70 - rSges(4,2) * t69 + t128 * t75 + t116) + g(2) * (rSges(4,1) * t68 - rSges(4,2) * t67 + t128 * t73 + t111) + g(3) * (rSges(4,1) * t78 - rSges(4,2) * t77 - t128 * t120 + t115)) - m(5) * (g(1) * (t70 * pkin(3) + t75 * pkin(8) + (t101 * t75 + t104 * t70) * rSges(5,1) + (-t101 * t70 + t104 * t75) * rSges(5,2) + t127 * t69 + t116) + g(2) * (t68 * pkin(3) + t73 * pkin(8) + (t101 * t73 + t104 * t68) * rSges(5,1) + (-t101 * t68 + t104 * t73) * rSges(5,2) + t127 * t67 + t111) + g(3) * (t78 * pkin(3) - pkin(8) * t120 + (-t101 * t120 + t104 * t78) * rSges(5,1) + (-t101 * t78 - t104 * t120) * rSges(5,2) + t127 * t77 + t115)) - m(6) * (g(1) * (t107 * t75 + t110 * t70 + t118 * t69 + t116) + g(2) * (t107 * t73 + t110 * t68 + t118 * t67 + t111) + g(3) * (-t107 * t120 + t110 * t78 + t118 * t77 + t115)) - m(7) * (g(1) * (t108 * t75 + t109 * t70 + t123 * t69 + t116) + g(2) * (t108 * t73 + t109 * t68 + t123 * t67 + t111) + g(3) * (-t108 * t120 + t109 * t78 + t123 * t77 + t115));
U  = t1;
