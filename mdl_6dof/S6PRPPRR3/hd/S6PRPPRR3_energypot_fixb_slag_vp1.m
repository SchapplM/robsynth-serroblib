% Calculate potential energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:21
% EndTime: 2019-03-08 19:21:21
% DurationCPUTime: 0.44s
% Computational Cost: add. (312->122), mult. (677->165), div. (0->0), fcn. (818->12), ass. (0->53)
t140 = rSges(6,3) + pkin(8);
t139 = pkin(9) + rSges(7,3);
t107 = sin(pkin(10));
t138 = g(1) * t107;
t110 = cos(pkin(10));
t137 = g(2) * t110;
t117 = cos(qJ(2));
t136 = qJ(3) * t117;
t108 = sin(pkin(6));
t135 = t107 * t108;
t134 = t108 * t110;
t113 = sin(qJ(5));
t133 = t108 * t113;
t114 = sin(qJ(2));
t132 = t108 * t114;
t116 = cos(qJ(5));
t131 = t108 * t116;
t111 = cos(pkin(6));
t130 = t111 * t114;
t129 = t111 * t117;
t128 = t111 * pkin(7) + qJ(1);
t127 = t110 * pkin(1) + pkin(7) * t135;
t126 = pkin(2) * t132 + t128;
t103 = t107 * pkin(1);
t93 = t107 * t114 - t110 * t129;
t94 = t107 * t117 + t110 * t130;
t125 = t94 * pkin(2) + t93 * qJ(3) + t103;
t95 = t107 * t129 + t110 * t114;
t96 = -t107 * t130 + t110 * t117;
t124 = t96 * pkin(2) + t95 * qJ(3) + t127;
t123 = t94 * pkin(3) + qJ(4) * t134 + t125;
t122 = t96 * pkin(3) + t124;
t121 = pkin(3) * t132 - t111 * qJ(4) + t126;
t106 = sin(pkin(11));
t109 = cos(pkin(11));
t79 = t106 * t93 + t109 * t94;
t120 = t79 * pkin(4) - pkin(7) * t134 + t123;
t81 = t106 * t95 + t109 * t96;
t119 = t81 * pkin(4) - qJ(4) * t135 + t122;
t88 = (-t106 * t117 + t109 * t114) * t108;
t118 = t88 * pkin(4) - t108 * t136 + t121;
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t87 = (t106 * t114 + t109 * t117) * t108;
t83 = -t111 * t113 + t116 * t88;
t82 = t111 * t116 + t113 * t88;
t80 = t106 * t96 - t95 * t109;
t78 = t106 * t94 - t93 * t109;
t75 = -t107 * t133 + t116 * t81;
t74 = t107 * t131 + t113 * t81;
t73 = t110 * t133 + t116 * t79;
t72 = -t110 * t131 + t113 * t79;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t110 - rSges(2,2) * t107) + g(2) * (rSges(2,1) * t107 + rSges(2,2) * t110) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t96 - rSges(3,2) * t95 + t127) + g(2) * (rSges(3,1) * t94 - rSges(3,2) * t93 + t103) + g(3) * (t111 * rSges(3,3) + t128) + (rSges(3,3) * t138 + g(3) * (rSges(3,1) * t114 + rSges(3,2) * t117) + (-rSges(3,3) - pkin(7)) * t137) * t108) - m(4) * (g(1) * (rSges(4,1) * t96 + rSges(4,3) * t95 + t124) + g(2) * (rSges(4,1) * t94 + rSges(4,3) * t93 + t125) + g(3) * (t111 * rSges(4,2) + t126) + (rSges(4,2) * t138 + g(3) * (rSges(4,1) * t114 - rSges(4,3) * t117 - t136) + (-rSges(4,2) - pkin(7)) * t137) * t108) - m(5) * (g(1) * (rSges(5,1) * t81 - rSges(5,2) * t80 + t122) + g(2) * (rSges(5,1) * t79 - rSges(5,2) * t78 + t123) + g(3) * (t88 * rSges(5,1) - t87 * rSges(5,2) - t111 * rSges(5,3) + t121) + (-g(3) * t136 + (rSges(5,3) - pkin(7)) * t137 + (-rSges(5,3) - qJ(4)) * t138) * t108) - m(6) * (g(1) * (rSges(6,1) * t75 - rSges(6,2) * t74 + t140 * t80 + t119) + g(2) * (rSges(6,1) * t73 - rSges(6,2) * t72 + t140 * t78 + t120) + g(3) * (t83 * rSges(6,1) - t82 * rSges(6,2) + t140 * t87 + t118)) - m(7) * (g(1) * (t75 * pkin(5) + t80 * pkin(8) + (t112 * t80 + t115 * t75) * rSges(7,1) + (-t112 * t75 + t115 * t80) * rSges(7,2) + t139 * t74 + t119) + g(2) * (t73 * pkin(5) + t78 * pkin(8) + (t112 * t78 + t115 * t73) * rSges(7,1) + (-t112 * t73 + t115 * t78) * rSges(7,2) + t139 * t72 + t120) + g(3) * (t83 * pkin(5) + t87 * pkin(8) + (t112 * t87 + t115 * t83) * rSges(7,1) + (-t112 * t83 + t115 * t87) * rSges(7,2) + t139 * t82 + t118));
U  = t1;
