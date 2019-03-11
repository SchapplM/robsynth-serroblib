% Calculate potential energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:01
% EndTime: 2019-03-08 21:04:02
% DurationCPUTime: 0.37s
% Computational Cost: add. (332->122), mult. (523->161), div. (0->0), fcn. (603->12), ass. (0->48)
t132 = pkin(8) + rSges(4,3);
t131 = pkin(9) + rSges(7,3);
t107 = sin(qJ(3));
t130 = pkin(3) * t107;
t103 = cos(pkin(10));
t101 = sin(pkin(10));
t102 = sin(pkin(6));
t126 = t101 * t102;
t129 = pkin(1) * t103 + pkin(7) * t126;
t128 = rSges(6,3) + qJ(5);
t104 = cos(pkin(6));
t127 = pkin(7) * t104 + qJ(1);
t125 = t102 * t103;
t124 = t102 * t107;
t108 = sin(qJ(2));
t123 = t102 * t108;
t110 = cos(qJ(3));
t122 = t102 * t110;
t111 = cos(qJ(2));
t121 = t102 * t111;
t120 = t104 * t108;
t119 = t104 * t111;
t118 = t101 * t124;
t105 = -qJ(4) - pkin(8);
t83 = t101 * t119 + t103 * t108;
t84 = -t101 * t120 + t103 * t111;
t94 = pkin(3) * t110 + pkin(2);
t117 = pkin(3) * t118 - t105 * t83 + t84 * t94 + t129;
t116 = t104 * t130 + t105 * t121 + t123 * t94 + t127;
t100 = qJ(3) + pkin(11);
t95 = sin(t100);
t96 = cos(t100);
t73 = t126 * t95 + t84 * t96;
t115 = pkin(4) * t73 + t117;
t78 = t104 * t95 + t123 * t96;
t114 = pkin(4) * t78 + t116;
t81 = t101 * t108 - t103 * t119;
t82 = t101 * t111 + t103 * t120;
t97 = t101 * pkin(1);
t113 = t82 * t94 + t97 + (-pkin(7) - t130) * t125 - t81 * t105;
t71 = -t125 * t95 + t82 * t96;
t112 = pkin(4) * t71 + t113;
t109 = cos(qJ(6));
t106 = sin(qJ(6));
t77 = -t104 * t96 + t123 * t95;
t72 = -t126 * t96 + t84 * t95;
t70 = t125 * t96 + t82 * t95;
t1 = -m(1) * (g(1) * rSges(1,1) + rSges(1,2) * g(2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t103 - rSges(2,2) * t101) + g(2) * (rSges(2,1) * t101 + rSges(2,2) * t103) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t84 - rSges(3,2) * t83 + t129) + g(2) * (t82 * rSges(3,1) - t81 * rSges(3,2) + t97) + g(3) * (t104 * rSges(3,3) + t127) + (g(1) * rSges(3,3) * t101 + g(3) * (rSges(3,1) * t108 + rSges(3,2) * t111) + g(2) * (-rSges(3,3) - pkin(7)) * t103) * t102) - m(4) * (g(1) * (t84 * pkin(2) + (t110 * t84 + t118) * rSges(4,1) + (t101 * t122 - t107 * t84) * rSges(4,2) + t132 * t83 + t129) + g(2) * (t82 * pkin(2) + t97 - pkin(7) * t125 + (-t103 * t124 + t110 * t82) * rSges(4,1) + (-t103 * t122 - t107 * t82) * rSges(4,2) + t132 * t81) + g(3) * ((rSges(4,1) * t107 + rSges(4,2) * t110) * t104 + (-t132 * t111 + (rSges(4,1) * t110 - rSges(4,2) * t107 + pkin(2)) * t108) * t102 + t127)) - m(5) * (g(1) * (rSges(5,1) * t73 - rSges(5,2) * t72 + rSges(5,3) * t83 + t117) + g(2) * (t71 * rSges(5,1) - t70 * rSges(5,2) + t81 * rSges(5,3) + t113) + g(3) * (rSges(5,1) * t78 - rSges(5,2) * t77 - rSges(5,3) * t121 + t116)) - m(6) * (g(1) * (t83 * rSges(6,1) - t73 * rSges(6,2) + t128 * t72 + t115) + g(2) * (t81 * rSges(6,1) - t71 * rSges(6,2) + t128 * t70 + t112) + g(3) * (-rSges(6,1) * t121 - t78 * rSges(6,2) + t128 * t77 + t114)) - m(7) * (g(1) * (t83 * pkin(5) + t72 * qJ(5) + (t106 * t72 + t109 * t83) * rSges(7,1) + (-t106 * t83 + t109 * t72) * rSges(7,2) + t131 * t73 + t115) + g(2) * (t81 * pkin(5) + t70 * qJ(5) + (t106 * t70 + t109 * t81) * rSges(7,1) + (-t106 * t81 + t109 * t70) * rSges(7,2) + t131 * t71 + t112) + g(3) * (-pkin(5) * t121 + t77 * qJ(5) + (t106 * t77 - t109 * t121) * rSges(7,1) + (t106 * t121 + t109 * t77) * rSges(7,2) + t131 * t78 + t114));
U  = t1;
