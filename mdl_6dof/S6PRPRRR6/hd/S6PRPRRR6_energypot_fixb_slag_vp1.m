% Calculate potential energy for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:05
% EndTime: 2019-03-08 20:45:05
% DurationCPUTime: 0.60s
% Computational Cost: add. (272->123), mult. (536->159), div. (0->0), fcn. (622->12), ass. (0->52)
t129 = -pkin(10) - pkin(9) - rSges(7,3);
t93 = sin(pkin(11));
t128 = g(1) * t93;
t95 = cos(pkin(11));
t127 = g(2) * t95;
t126 = pkin(9) + rSges(6,3);
t102 = cos(qJ(2));
t118 = t93 * t102;
t96 = cos(pkin(6));
t99 = sin(qJ(2));
t121 = t96 * t99;
t74 = t121 * t95 + t118;
t97 = sin(qJ(5));
t125 = t74 * t97;
t94 = sin(pkin(6));
t124 = t93 * t94;
t98 = sin(qJ(4));
t123 = t94 * t98;
t122 = t94 * t99;
t120 = t95 * pkin(1) + pkin(7) * t124;
t101 = cos(qJ(4));
t119 = t101 * t94;
t117 = t94 * t102;
t116 = t95 * t102;
t115 = t96 * pkin(7) + qJ(1);
t114 = t102 * qJ(3);
t113 = pkin(2) * t122 + t115;
t112 = (-pkin(3) - pkin(7)) * t95;
t111 = g(2) * t112;
t73 = -t116 * t96 + t93 * t99;
t88 = t93 * pkin(1);
t110 = t74 * pkin(2) + t73 * qJ(3) + t88;
t109 = t96 * pkin(3) + pkin(8) * t122 + t113;
t75 = t118 * t96 + t95 * t99;
t76 = -t121 * t93 + t116;
t108 = t76 * pkin(2) + t75 * qJ(3) + t120;
t100 = cos(qJ(5));
t85 = t100 * pkin(5) + pkin(4);
t92 = qJ(5) + qJ(6);
t86 = sin(t92);
t87 = cos(t92);
t107 = t87 * rSges(7,1) - t86 * rSges(7,2) + t85;
t106 = pkin(3) * t124 + t108;
t105 = t74 * pkin(8) + t110;
t104 = t86 * rSges(7,1) + t87 * rSges(7,2) + t97 * pkin(5);
t78 = t96 * t101 - t117 * t98;
t77 = t101 * t117 + t96 * t98;
t69 = -t119 * t95 + t73 * t98;
t68 = t73 * t101 + t123 * t95;
t67 = t119 * t93 + t75 * t98;
t66 = -t75 * t101 + t123 * t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t95 * rSges(2,1) - t93 * rSges(2,2)) + g(2) * (t93 * rSges(2,1) + t95 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t76 * rSges(3,1) - t75 * rSges(3,2) + t120) + g(2) * (t74 * rSges(3,1) - t73 * rSges(3,2) + t88) + g(3) * (t96 * rSges(3,3) + t115) + (rSges(3,3) * t128 + g(3) * (rSges(3,1) * t99 + rSges(3,2) * t102) + (-rSges(3,3) - pkin(7)) * t127) * t94) - m(4) * (g(1) * (-t76 * rSges(4,2) + t75 * rSges(4,3) + t108) + g(2) * (-t74 * rSges(4,2) + t73 * rSges(4,3) + t110) + g(3) * (t96 * rSges(4,1) + t113) + (rSges(4,1) * t128 + g(3) * (-rSges(4,2) * t99 - rSges(4,3) * t102 - t114) + (-rSges(4,1) - pkin(7)) * t127) * t94) - m(5) * (g(1) * (t67 * rSges(5,1) - t66 * rSges(5,2) + (rSges(5,3) + pkin(8)) * t76 + t106) + g(2) * (t69 * rSges(5,1) + t68 * rSges(5,2) + t74 * rSges(5,3) + t105) + g(3) * (t78 * rSges(5,1) - t77 * rSges(5,2) + t109) + (g(3) * (rSges(5,3) * t99 - t114) + t111) * t94) - m(6) * (g(1) * (t67 * pkin(4) + t76 * pkin(8) + (t67 * t100 + t76 * t97) * rSges(6,1) + (t76 * t100 - t67 * t97) * rSges(6,2) + t126 * t66 + t106) + g(2) * (t69 * pkin(4) + (t69 * t100 + t125) * rSges(6,1) + (t74 * t100 - t69 * t97) * rSges(6,2) + t94 * t112 - t126 * t68 + t105) + g(3) * (t78 * pkin(4) - t94 * t114 + (t78 * t100 + t122 * t97) * rSges(6,1) + (t100 * t122 - t78 * t97) * rSges(6,2) + t126 * t77 + t109)) - m(7) * (g(1) * (t107 * t67 + (pkin(8) + t104) * t76 - t129 * t66 + t106) + g(2) * (t69 * t85 + pkin(5) * t125 + (t69 * t87 + t74 * t86) * rSges(7,1) + (-t69 * t86 + t74 * t87) * rSges(7,2) + t105 + t129 * t68) + t111 * t94 + (t109 + t107 * t78 + (t104 * t99 - t114) * t94 - t129 * t77) * g(3));
U  = t1;
