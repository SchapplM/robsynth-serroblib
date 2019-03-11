% Calculate potential energy for
% S6PRPRRR5
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:48
% EndTime: 2019-03-08 20:40:48
% DurationCPUTime: 0.64s
% Computational Cost: add. (284->121), mult. (482->152), div. (0->0), fcn. (547->12), ass. (0->51)
t132 = pkin(8) + rSges(5,3);
t131 = pkin(10) + rSges(7,3);
t98 = sin(pkin(11));
t130 = g(1) * t98;
t103 = sin(qJ(4));
t129 = pkin(4) * t103;
t100 = cos(pkin(11));
t128 = g(2) * t100;
t99 = sin(pkin(6));
t127 = t98 * t99;
t107 = cos(qJ(2));
t101 = cos(pkin(6));
t104 = sin(qJ(2));
t121 = t101 * t104;
t82 = t100 * t121 + t107 * t98;
t94 = t98 * pkin(1);
t126 = t82 * pkin(2) + t94;
t125 = t100 * pkin(1) + pkin(7) * t127;
t124 = t100 * t99;
t123 = t107 * t99;
t122 = t101 * pkin(7) + qJ(1);
t120 = t101 * t107;
t84 = t100 * t107 - t98 * t121;
t119 = t84 * pkin(2) + t125;
t118 = t99 * t104 * pkin(2) + t122;
t106 = cos(qJ(4));
t91 = pkin(4) * t106 + pkin(3);
t117 = t101 * t91 + t118;
t116 = (-pkin(7) - t91) * t128;
t81 = -t100 * t120 + t104 * t98;
t115 = t81 * qJ(3) + t126;
t114 = (-qJ(3) - t129) * t107;
t83 = t100 * t104 + t98 * t120;
t113 = t83 * qJ(3) + t119;
t112 = t106 * rSges(5,1) - t103 * rSges(5,2) + pkin(3);
t111 = t103 * rSges(5,1) + t106 * rSges(5,2) + qJ(3);
t108 = -pkin(9) - pkin(8);
t110 = -t82 * t108 + t81 * t129 + t115;
t109 = -t84 * t108 + t91 * t127 + t83 * t129 + t113;
t105 = cos(qJ(6));
t102 = sin(qJ(6));
t97 = qJ(4) + qJ(5);
t93 = cos(t97);
t92 = sin(t97);
t76 = t101 * t93 - t92 * t123;
t75 = t101 * t92 + t93 * t123;
t71 = -t93 * t124 + t81 * t92;
t70 = t92 * t124 + t81 * t93;
t69 = t93 * t127 + t83 * t92;
t68 = t92 * t127 - t83 * t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t100 - rSges(2,2) * t98) + g(2) * (rSges(2,1) * t98 + rSges(2,2) * t100) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t125) + g(2) * (rSges(3,1) * t82 - rSges(3,2) * t81 + t94) + g(3) * (rSges(3,3) * t101 + t122) + (rSges(3,3) * t130 + g(3) * (rSges(3,1) * t104 + rSges(3,2) * t107) + (-rSges(3,3) - pkin(7)) * t128) * t99) - m(4) * (g(1) * (-rSges(4,2) * t84 + rSges(4,3) * t83 + t113) + g(2) * (-t82 * rSges(4,2) + rSges(4,3) * t81 + t115) + g(3) * (rSges(4,1) * t101 + t118) + (rSges(4,1) * t130 + g(3) * (-rSges(4,2) * t104 + (-rSges(4,3) - qJ(3)) * t107) + (-rSges(4,1) - pkin(7)) * t128) * t99) - m(5) * ((t112 * t130 + (-pkin(7) - t112) * t128) * t99 + (t118 + t112 * t101 + (t132 * t104 - t111 * t107) * t99) * g(3) + (t111 * t81 + t132 * t82 + t126) * g(2) + (t111 * t83 + t132 * t84 + t119) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t69 - rSges(6,2) * t68 + rSges(6,3) * t84 + t109) + g(2) * (rSges(6,1) * t71 + rSges(6,2) * t70 + rSges(6,3) * t82 + t110) + g(3) * (t76 * rSges(6,1) - t75 * rSges(6,2) + t117) + (g(3) * (t114 + (rSges(6,3) - t108) * t104) + t116) * t99) - m(7) * (g(1) * (t69 * pkin(5) + (t102 * t84 + t105 * t69) * rSges(7,1) + (-t102 * t69 + t105 * t84) * rSges(7,2) + t131 * t68 + t109) + g(2) * (t71 * pkin(5) + (t102 * t82 + t105 * t71) * rSges(7,1) + (-t102 * t71 + t105 * t82) * rSges(7,2) + t110 - t131 * t70) + t116 * t99 + (t117 + (t105 * rSges(7,1) - t102 * rSges(7,2) + pkin(5)) * t76 + (t114 + (t102 * rSges(7,1) + t105 * rSges(7,2) - t108) * t104) * t99 + t131 * t75) * g(3));
U  = t1;
