% Calculate potential energy for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:04
% EndTime: 2019-03-09 19:38:04
% DurationCPUTime: 0.40s
% Computational Cost: add. (335->116), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->51)
t98 = sin(pkin(12));
t130 = pkin(4) * t98 + pkin(9);
t128 = rSges(4,3) + pkin(9);
t117 = cos(pkin(6));
t127 = pkin(8) * t117 + pkin(7);
t100 = cos(pkin(12));
t88 = pkin(4) * t100 + pkin(3);
t126 = cos(qJ(3));
t101 = -pkin(10) - qJ(4);
t106 = cos(qJ(1));
t104 = sin(qJ(1));
t99 = sin(pkin(6));
t122 = t104 * t99;
t125 = pkin(1) * t106 + pkin(8) * t122;
t124 = pkin(11) - t101 + rSges(7,3);
t103 = sin(qJ(2));
t123 = t103 * t99;
t105 = cos(qJ(2));
t121 = t105 * t99;
t120 = t106 * t99;
t119 = qJ(4) + rSges(5,3);
t118 = -t101 + rSges(6,3);
t97 = pkin(12) + qJ(5);
t116 = pkin(2) * t123 + t127;
t113 = t104 * t117;
t78 = -t103 * t113 + t105 * t106;
t115 = pkin(2) * t78 + t125;
t114 = t99 * t126;
t112 = t106 * t117;
t76 = t103 * t112 + t104 * t105;
t94 = t104 * pkin(1);
t111 = pkin(2) * t76 - pkin(8) * t120 + t94;
t89 = sin(t97);
t90 = cos(t97);
t110 = rSges(6,1) * t90 - rSges(6,2) * t89 + t88;
t91 = qJ(6) + t97;
t86 = sin(t91);
t87 = cos(t91);
t109 = rSges(7,1) * t87 - rSges(7,2) * t86 + pkin(5) * t90 + t88;
t108 = rSges(7,1) * t86 + rSges(7,2) * t87 + pkin(5) * t89 + t130;
t107 = rSges(6,1) * t89 + rSges(6,2) * t90 + t130;
t102 = sin(qJ(3));
t77 = t103 * t106 + t105 * t113;
t75 = t103 * t104 - t105 * t112;
t74 = t102 * t117 + t103 * t114;
t73 = t102 * t123 - t117 * t126;
t70 = t102 * t122 + t126 * t78;
t69 = t102 * t78 - t104 * t114;
t68 = -t102 * t120 + t126 * t76;
t67 = t102 * t76 + t106 * t114;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t106 - rSges(2,2) * t104) + g(2) * (rSges(2,1) * t104 + rSges(2,2) * t106) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t78 - rSges(3,2) * t77 + t125) + g(2) * (t76 * rSges(3,1) - t75 * rSges(3,2) + t94) + g(3) * (t117 * rSges(3,3) + t127) + (g(1) * rSges(3,3) * t104 + g(3) * (rSges(3,1) * t103 + rSges(3,2) * t105) + g(2) * (-rSges(3,3) - pkin(8)) * t106) * t99) - m(4) * (g(1) * (rSges(4,1) * t70 - rSges(4,2) * t69 + t128 * t77 + t115) + g(2) * (t68 * rSges(4,1) - t67 * rSges(4,2) + t128 * t75 + t111) + g(3) * (rSges(4,1) * t74 - rSges(4,2) * t73 - t121 * t128 + t116)) - m(5) * (g(1) * (t70 * pkin(3) + t77 * pkin(9) + (t100 * t70 + t77 * t98) * rSges(5,1) + (t100 * t77 - t70 * t98) * rSges(5,2) + t119 * t69 + t115) + g(2) * (t68 * pkin(3) + t75 * pkin(9) + (t100 * t68 + t75 * t98) * rSges(5,1) + (t100 * t75 - t68 * t98) * rSges(5,2) + t119 * t67 + t111) + g(3) * (t74 * pkin(3) - pkin(9) * t121 + (t100 * t74 - t121 * t98) * rSges(5,1) + (-t100 * t121 - t74 * t98) * rSges(5,2) + t119 * t73 + t116)) - m(6) * (g(1) * (t107 * t77 + t110 * t70 + t118 * t69 + t115) + g(2) * (t107 * t75 + t110 * t68 + t118 * t67 + t111) + g(3) * (-t107 * t121 + t110 * t74 + t118 * t73 + t116)) - m(7) * (g(1) * (t108 * t77 + t109 * t70 + t124 * t69 + t115) + g(2) * (t108 * t75 + t109 * t68 + t124 * t67 + t111) + g(3) * (-t108 * t121 + t109 * t74 + t124 * t73 + t116));
U  = t1;
