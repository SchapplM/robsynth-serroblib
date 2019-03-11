% Calculate potential energy for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:12:14
% EndTime: 2019-03-09 23:12:15
% DurationCPUTime: 0.40s
% Computational Cost: add. (335->116), mult. (604->146), div. (0->0), fcn. (717->14), ass. (0->51)
t100 = sin(qJ(4));
t130 = pkin(4) * t100 + pkin(9);
t129 = rSges(4,3) + pkin(9);
t117 = cos(pkin(6));
t128 = pkin(8) * t117 + pkin(7);
t127 = pkin(10) + rSges(5,3);
t104 = cos(qJ(4));
t88 = pkin(4) * t104 + pkin(3);
t126 = cos(qJ(3));
t99 = -qJ(5) - pkin(10);
t106 = cos(qJ(1));
t103 = sin(qJ(1));
t98 = sin(pkin(6));
t120 = t103 * t98;
t124 = pkin(1) * t106 + pkin(8) * t120;
t123 = pkin(11) - t99 + rSges(7,3);
t122 = -t99 + rSges(6,3);
t102 = sin(qJ(2));
t121 = t102 * t98;
t105 = cos(qJ(2));
t119 = t105 * t98;
t118 = t106 * t98;
t97 = qJ(4) + pkin(12);
t116 = pkin(2) * t121 + t128;
t113 = t103 * t117;
t78 = -t102 * t113 + t105 * t106;
t115 = pkin(2) * t78 + t124;
t114 = t98 * t126;
t112 = t106 * t117;
t76 = t102 * t112 + t103 * t105;
t93 = t103 * pkin(1);
t111 = pkin(2) * t76 - pkin(8) * t118 + t93;
t89 = sin(t97);
t90 = cos(t97);
t110 = rSges(6,1) * t90 - rSges(6,2) * t89 + t88;
t91 = qJ(6) + t97;
t86 = sin(t91);
t87 = cos(t91);
t109 = rSges(7,1) * t87 - rSges(7,2) * t86 + pkin(5) * t90 + t88;
t108 = rSges(7,1) * t86 + rSges(7,2) * t87 + pkin(5) * t89 + t130;
t107 = rSges(6,1) * t89 + rSges(6,2) * t90 + t130;
t101 = sin(qJ(3));
t77 = t102 * t106 + t105 * t113;
t75 = t102 * t103 - t105 * t112;
t74 = t101 * t117 + t102 * t114;
t73 = t101 * t121 - t117 * t126;
t70 = t101 * t120 + t126 * t78;
t69 = t101 * t78 - t103 * t114;
t68 = -t101 * t118 + t126 * t76;
t67 = t101 * t76 + t106 * t114;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + rSges(1,3) * g(3)) - m(2) * (g(1) * (rSges(2,1) * t106 - rSges(2,2) * t103) + g(2) * (rSges(2,1) * t103 + rSges(2,2) * t106) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t78 - rSges(3,2) * t77 + t124) + g(2) * (t76 * rSges(3,1) - t75 * rSges(3,2) + t93) + g(3) * (t117 * rSges(3,3) + t128) + (g(1) * rSges(3,3) * t103 + g(3) * (rSges(3,1) * t102 + rSges(3,2) * t105) + g(2) * (-rSges(3,3) - pkin(8)) * t106) * t98) - m(4) * (g(1) * (rSges(4,1) * t70 - rSges(4,2) * t69 + t129 * t77 + t115) + g(2) * (t68 * rSges(4,1) - t67 * rSges(4,2) + t129 * t75 + t111) + g(3) * (rSges(4,1) * t74 - rSges(4,2) * t73 - t119 * t129 + t116)) - m(5) * (g(1) * (t70 * pkin(3) + t77 * pkin(9) + (t100 * t77 + t104 * t70) * rSges(5,1) + (-t100 * t70 + t104 * t77) * rSges(5,2) + t127 * t69 + t115) + g(2) * (t68 * pkin(3) + t75 * pkin(9) + (t100 * t75 + t104 * t68) * rSges(5,1) + (-t100 * t68 + t104 * t75) * rSges(5,2) + t127 * t67 + t111) + g(3) * (t74 * pkin(3) - pkin(9) * t119 + (-t100 * t119 + t104 * t74) * rSges(5,1) + (-t100 * t74 - t104 * t119) * rSges(5,2) + t127 * t73 + t116)) - m(6) * (g(1) * (t107 * t77 + t110 * t70 + t122 * t69 + t115) + g(2) * (t107 * t75 + t110 * t68 + t122 * t67 + t111) + g(3) * (-t107 * t119 + t110 * t74 + t122 * t73 + t116)) - m(7) * (g(1) * (t108 * t77 + t109 * t70 + t123 * t69 + t115) + g(2) * (t108 * t75 + t109 * t68 + t123 * t67 + t111) + g(3) * (-t108 * t119 + t109 * t74 + t123 * t73 + t116));
U  = t1;
