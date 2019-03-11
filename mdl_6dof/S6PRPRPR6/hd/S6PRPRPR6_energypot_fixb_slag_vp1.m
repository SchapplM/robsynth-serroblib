% Calculate potential energy for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:34
% EndTime: 2019-03-08 19:46:34
% DurationCPUTime: 0.64s
% Computational Cost: add. (272->123), mult. (536->158), div. (0->0), fcn. (622->12), ass. (0->53)
t130 = -pkin(9) - qJ(5) - rSges(7,3);
t94 = sin(pkin(10));
t129 = g(1) * t94;
t97 = cos(pkin(10));
t128 = g(2) * t97;
t101 = sin(qJ(2));
t118 = t97 * t101;
t103 = cos(qJ(2));
t120 = t94 * t103;
t98 = cos(pkin(6));
t74 = t118 * t98 + t120;
t93 = sin(pkin(11));
t127 = t74 * t93;
t95 = sin(pkin(6));
t126 = t94 * t95;
t125 = t97 * pkin(1) + pkin(7) * t126;
t100 = sin(qJ(4));
t124 = t100 * t95;
t123 = t101 * t95;
t102 = cos(qJ(4));
t122 = t102 * t95;
t121 = t94 * t101;
t119 = t95 * t103;
t117 = t97 * t103;
t116 = t98 * pkin(7) + qJ(1);
t115 = qJ(5) + rSges(6,3);
t114 = t103 * qJ(3);
t113 = pkin(2) * t123 + t116;
t112 = (-pkin(3) - pkin(7)) * t97;
t111 = g(2) * t112;
t73 = -t117 * t98 + t121;
t88 = t94 * pkin(1);
t110 = t74 * pkin(2) + t73 * qJ(3) + t88;
t109 = t98 * pkin(3) + pkin(8) * t123 + t113;
t75 = t120 * t98 + t118;
t76 = -t121 * t98 + t117;
t108 = t76 * pkin(2) + t75 * qJ(3) + t125;
t96 = cos(pkin(11));
t85 = t96 * pkin(5) + pkin(4);
t92 = pkin(11) + qJ(6);
t86 = sin(t92);
t87 = cos(t92);
t107 = t87 * rSges(7,1) - t86 * rSges(7,2) + t85;
t106 = pkin(3) * t126 + t108;
t105 = t74 * pkin(8) + t110;
t104 = t86 * rSges(7,1) + t87 * rSges(7,2) + t93 * pkin(5);
t78 = -t100 * t119 + t98 * t102;
t77 = t98 * t100 + t102 * t119;
t69 = t73 * t100 - t122 * t97;
t68 = t73 * t102 + t124 * t97;
t67 = t75 * t100 + t122 * t94;
t66 = -t75 * t102 + t124 * t94;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t97 * rSges(2,1) - t94 * rSges(2,2)) + g(2) * (t94 * rSges(2,1) + t97 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t76 * rSges(3,1) - t75 * rSges(3,2) + t125) + g(2) * (t74 * rSges(3,1) - t73 * rSges(3,2) + t88) + g(3) * (t98 * rSges(3,3) + t116) + (rSges(3,3) * t129 + g(3) * (rSges(3,1) * t101 + rSges(3,2) * t103) + (-rSges(3,3) - pkin(7)) * t128) * t95) - m(4) * (g(1) * (-t76 * rSges(4,2) + t75 * rSges(4,3) + t108) + g(2) * (-t74 * rSges(4,2) + t73 * rSges(4,3) + t110) + g(3) * (t98 * rSges(4,1) + t113) + (rSges(4,1) * t129 + g(3) * (-rSges(4,2) * t101 - rSges(4,3) * t103 - t114) + (-rSges(4,1) - pkin(7)) * t128) * t95) - m(5) * (g(1) * (t67 * rSges(5,1) - t66 * rSges(5,2) + (rSges(5,3) + pkin(8)) * t76 + t106) + g(2) * (t69 * rSges(5,1) + t68 * rSges(5,2) + t74 * rSges(5,3) + t105) + g(3) * (t78 * rSges(5,1) - t77 * rSges(5,2) + t109) + (g(3) * (rSges(5,3) * t101 - t114) + t111) * t95) - m(6) * (g(1) * (t67 * pkin(4) + t76 * pkin(8) + (t67 * t96 + t76 * t93) * rSges(6,1) + (-t67 * t93 + t76 * t96) * rSges(6,2) + t115 * t66 + t106) + g(2) * (t69 * pkin(4) + (t69 * t96 + t127) * rSges(6,1) + (-t69 * t93 + t74 * t96) * rSges(6,2) + t95 * t112 - t115 * t68 + t105) + g(3) * (t78 * pkin(4) - t95 * t114 + (t123 * t93 + t78 * t96) * rSges(6,1) + (t123 * t96 - t78 * t93) * rSges(6,2) + t115 * t77 + t109)) - m(7) * (g(1) * (t107 * t67 + (pkin(8) + t104) * t76 - t130 * t66 + t106) + g(2) * (t69 * t85 + pkin(5) * t127 + (t69 * t87 + t74 * t86) * rSges(7,1) + (-t69 * t86 + t74 * t87) * rSges(7,2) + t105 + t130 * t68) + t111 * t95 + (t109 + t107 * t78 + (t101 * t104 - t114) * t95 - t130 * t77) * g(3));
U  = t1;
