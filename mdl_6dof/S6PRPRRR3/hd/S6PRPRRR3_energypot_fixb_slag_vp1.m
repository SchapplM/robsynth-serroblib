% Calculate potential energy for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:35
% EndTime: 2019-03-08 20:31:36
% DurationCPUTime: 0.63s
% Computational Cost: add. (357->125), mult. (471->161), div. (0->0), fcn. (534->14), ass. (0->52)
t111 = -pkin(8) - qJ(3);
t134 = -t111 + rSges(5,3);
t133 = pkin(10) + rSges(7,3);
t105 = sin(pkin(12));
t132 = pkin(3) * t105;
t106 = sin(pkin(11));
t131 = g(1) * t106;
t109 = cos(pkin(11));
t130 = g(2) * t109;
t108 = cos(pkin(12));
t95 = t108 * pkin(3) + pkin(2);
t129 = qJ(3) + rSges(4,3);
t107 = sin(pkin(6));
t127 = t106 * t107;
t128 = t109 * pkin(1) + pkin(7) * t127;
t126 = t107 * t109;
t113 = sin(qJ(2));
t125 = t107 * t113;
t115 = cos(qJ(2));
t124 = t107 * t115;
t110 = cos(pkin(6));
t123 = t110 * t113;
t122 = t110 * t115;
t121 = t110 * pkin(7) + qJ(1);
t104 = pkin(12) + qJ(4);
t103 = -pkin(9) + t111;
t83 = t106 * t122 + t109 * t113;
t84 = -t106 * t123 + t109 * t115;
t97 = cos(t104);
t87 = pkin(4) * t97 + t95;
t96 = sin(t104);
t88 = pkin(4) * t96 + t132;
t120 = -t83 * t103 + t88 * t127 + t84 * t87 + t128;
t119 = t103 * t124 + t110 * t88 + t87 * t125 + t121;
t118 = t97 * rSges(5,1) - t96 * rSges(5,2) + t95;
t117 = rSges(5,1) * t96 + rSges(5,2) * t97 + t132;
t81 = t106 * t113 - t109 * t122;
t82 = t106 * t115 + t109 * t123;
t99 = t106 * pkin(1);
t116 = t82 * t87 + t99 + (-pkin(7) - t88) * t126 - t81 * t103;
t114 = cos(qJ(6));
t112 = sin(qJ(6));
t98 = qJ(5) + t104;
t94 = cos(t98);
t93 = sin(t98);
t76 = t110 * t93 + t125 * t94;
t75 = -t110 * t94 + t125 * t93;
t72 = t127 * t93 + t84 * t94;
t71 = -t127 * t94 + t84 * t93;
t70 = -t126 * t93 + t82 * t94;
t69 = t126 * t94 + t82 * t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t109 - rSges(2,2) * t106) + g(2) * (rSges(2,1) * t106 + rSges(2,2) * t109) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t128) + g(2) * (rSges(3,1) * t82 - rSges(3,2) * t81 + t99) + g(3) * (t110 * rSges(3,3) + t121) + (rSges(3,3) * t131 + g(3) * (rSges(3,1) * t113 + rSges(3,2) * t115) + (-rSges(3,3) - pkin(7)) * t130) * t107) - m(4) * (g(1) * (t84 * pkin(2) + (t105 * t127 + t108 * t84) * rSges(4,1) + (-t105 * t84 + t108 * t127) * rSges(4,2) + t129 * t83 + t128) + g(2) * (t82 * pkin(2) + t99 - pkin(7) * t126 + (-t105 * t126 + t108 * t82) * rSges(4,1) + (-t105 * t82 - t108 * t126) * rSges(4,2) + t129 * t81) + g(3) * ((rSges(4,1) * t105 + rSges(4,2) * t108) * t110 + (-t129 * t115 + (t108 * rSges(4,1) - t105 * rSges(4,2) + pkin(2)) * t113) * t107 + t121)) - m(5) * ((t117 * t131 + (-pkin(7) - t117) * t130) * t107 + (t121 + t117 * t110 + (t118 * t113 - t115 * t134) * t107) * g(3) + (t118 * t82 + t134 * t81 + t99) * g(2) + (t118 * t84 + t134 * t83 + t128) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t72 - rSges(6,2) * t71 + rSges(6,3) * t83 + t120) + g(2) * (rSges(6,1) * t70 - rSges(6,2) * t69 + rSges(6,3) * t81 + t116) + g(3) * (t76 * rSges(6,1) - t75 * rSges(6,2) - rSges(6,3) * t124 + t119)) - m(7) * (g(1) * (t72 * pkin(5) + (t112 * t83 + t114 * t72) * rSges(7,1) + (-t112 * t72 + t114 * t83) * rSges(7,2) + t133 * t71 + t120) + g(2) * (t70 * pkin(5) + (t112 * t81 + t114 * t70) * rSges(7,1) + (-t112 * t70 + t114 * t81) * rSges(7,2) + t133 * t69 + t116) + g(3) * (t76 * pkin(5) + (-t112 * t124 + t76 * t114) * rSges(7,1) + (-t76 * t112 - t114 * t124) * rSges(7,2) + t133 * t75 + t119));
U  = t1;
