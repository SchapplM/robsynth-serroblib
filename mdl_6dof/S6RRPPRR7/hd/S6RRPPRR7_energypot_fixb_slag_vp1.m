% Calculate potential energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:04
% EndTime: 2019-03-09 09:17:04
% DurationCPUTime: 0.46s
% Computational Cost: add. (232->119), mult. (469->153), div. (0->0), fcn. (530->10), ass. (0->46)
t96 = sin(qJ(1));
t124 = g(1) * t96;
t92 = cos(pkin(6));
t123 = t92 * pkin(8) + pkin(7);
t122 = pkin(10) + rSges(7,3);
t100 = cos(qJ(1));
t121 = g(2) * t100;
t91 = sin(pkin(6));
t95 = sin(qJ(2));
t120 = t91 * t95;
t99 = cos(qJ(2));
t119 = t91 * t99;
t118 = t96 * t91;
t117 = t96 * t95;
t116 = t96 * t99;
t115 = t100 * pkin(1) + pkin(8) * t118;
t114 = t100 * t91;
t113 = t100 * t95;
t112 = t100 * t99;
t111 = t96 * qJ(4);
t110 = t99 * qJ(3);
t109 = pkin(2) * t120 + t123;
t76 = -t92 * t112 + t117;
t77 = t92 * t113 + t116;
t89 = t96 * pkin(1);
t108 = t77 * pkin(2) + t76 * qJ(3) + t89;
t78 = t92 * t116 + t113;
t79 = -t92 * t117 + t112;
t107 = t79 * pkin(2) + t78 * qJ(3) + t115;
t106 = pkin(3) * t120 - t92 * qJ(4) + t109;
t105 = t77 * pkin(3) + qJ(4) * t114 + t108;
t104 = t79 * pkin(3) + t107;
t103 = pkin(9) * t120 + t106;
t102 = t76 * pkin(4) + t77 * pkin(9) + t105;
t101 = t78 * pkin(4) + t79 * pkin(9) + t104;
t98 = cos(qJ(5));
t97 = cos(qJ(6));
t94 = sin(qJ(5));
t93 = sin(qJ(6));
t75 = -t98 * t119 - t92 * t94;
t74 = t94 * t119 - t92 * t98;
t67 = -t94 * t118 + t78 * t98;
t66 = t98 * t118 + t78 * t94;
t65 = t94 * t114 + t76 * t98;
t64 = -t98 * t114 + t76 * t94;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t100 * rSges(2,1) - t96 * rSges(2,2)) + g(2) * (t96 * rSges(2,1) + t100 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t79 * rSges(3,1) - t78 * rSges(3,2) + t115) + g(2) * (t77 * rSges(3,1) - t76 * rSges(3,2) + t89) + g(3) * (t92 * rSges(3,3) + t123) + (rSges(3,3) * t124 + g(3) * (rSges(3,1) * t95 + rSges(3,2) * t99) + (-rSges(3,3) - pkin(8)) * t121) * t91) - m(4) * (g(1) * (t79 * rSges(4,1) + t78 * rSges(4,3) + t107) + g(2) * (t77 * rSges(4,1) + t76 * rSges(4,3) + t108) + g(3) * (t92 * rSges(4,2) + t109) + (rSges(4,2) * t124 + g(3) * (rSges(4,1) * t95 - rSges(4,3) * t99 - t110) + (-rSges(4,2) - pkin(8)) * t121) * t91) - m(5) * (g(1) * (t78 * rSges(5,1) - t79 * rSges(5,2) + t104) + g(2) * (t76 * rSges(5,1) - t77 * rSges(5,2) + t105) + g(3) * (-t92 * rSges(5,3) + t106) + (g(3) * (-rSges(5,2) * t95 + (-rSges(5,1) - qJ(3)) * t99) + (-rSges(5,3) - qJ(4)) * t124 + (rSges(5,3) - pkin(8)) * t121) * t91) - m(6) * (g(1) * (t67 * rSges(6,1) - t66 * rSges(6,2) + t79 * rSges(6,3) + t101) + g(2) * (t65 * rSges(6,1) - t64 * rSges(6,2) + t77 * rSges(6,3) + t102) + g(3) * (t75 * rSges(6,1) + t74 * rSges(6,2) + t103) + (-g(1) * t111 - pkin(8) * t121 + g(3) * (rSges(6,3) * t95 - t99 * pkin(4) - t110)) * t91) - m(7) * (g(1) * (t67 * pkin(5) - t91 * t111 + (t67 * t97 + t79 * t93) * rSges(7,1) + (-t67 * t93 + t79 * t97) * rSges(7,2) + t122 * t66 + t101) + g(2) * (t65 * pkin(5) - pkin(8) * t114 + (t65 * t97 + t77 * t93) * rSges(7,1) + (-t65 * t93 + t77 * t97) * rSges(7,2) + t122 * t64 + t102) + g(3) * ((t97 * rSges(7,1) - t93 * rSges(7,2) + pkin(5)) * t75 + ((-pkin(4) - qJ(3)) * t99 + (t93 * rSges(7,1) + t97 * rSges(7,2)) * t95) * t91 - t122 * t74 + t103));
U  = t1;
