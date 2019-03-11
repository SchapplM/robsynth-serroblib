% Calculate potential energy for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:13
% EndTime: 2019-03-09 14:08:14
% DurationCPUTime: 0.64s
% Computational Cost: add. (357->125), mult. (471->159), div. (0->0), fcn. (534->14), ass. (0->54)
t111 = -pkin(9) - qJ(3);
t138 = -t111 + rSges(5,3);
t137 = pkin(11) + rSges(7,3);
t107 = sin(pkin(12));
t136 = pkin(3) * t107;
t114 = sin(qJ(1));
t135 = g(1) * t114;
t117 = cos(qJ(1));
t134 = g(2) * t117;
t110 = cos(pkin(6));
t133 = t110 * pkin(8) + pkin(7);
t109 = cos(pkin(12));
t97 = t109 * pkin(3) + pkin(2);
t132 = qJ(3) + rSges(4,3);
t108 = sin(pkin(6));
t129 = t108 * t114;
t131 = t117 * pkin(1) + pkin(8) * t129;
t113 = sin(qJ(2));
t130 = t108 * t113;
t116 = cos(qJ(2));
t128 = t108 * t116;
t127 = t108 * t117;
t126 = t113 * t114;
t125 = t113 * t117;
t124 = t114 * t116;
t123 = t116 * t117;
t106 = pkin(12) + qJ(4);
t105 = -pkin(10) + t111;
t99 = cos(t106);
t88 = pkin(4) * t99 + t97;
t98 = sin(t106);
t90 = pkin(4) * t98 + t136;
t122 = t105 * t128 + t110 * t90 + t88 * t130 + t133;
t85 = t110 * t124 + t125;
t86 = -t110 * t126 + t123;
t121 = -t85 * t105 + t90 * t129 + t86 * t88 + t131;
t120 = rSges(5,1) * t99 - rSges(5,2) * t98 + t97;
t119 = rSges(5,1) * t98 + rSges(5,2) * t99 + t136;
t103 = t114 * pkin(1);
t83 = -t110 * t123 + t126;
t84 = t110 * t125 + t124;
t118 = t103 + t84 * t88 + (-pkin(8) - t90) * t127 - t83 * t105;
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t100 = qJ(5) + t106;
t96 = cos(t100);
t95 = sin(t100);
t78 = t110 * t95 + t130 * t96;
t77 = -t110 * t96 + t130 * t95;
t74 = t129 * t95 + t86 * t96;
t73 = -t129 * t96 + t86 * t95;
t72 = -t127 * t95 + t84 * t96;
t71 = t127 * t96 + t84 * t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t117 - t114 * rSges(2,2)) + g(2) * (t114 * rSges(2,1) + rSges(2,2) * t117) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t86 - rSges(3,2) * t85 + t131) + g(2) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t103) + g(3) * (rSges(3,3) * t110 + t133) + (rSges(3,3) * t135 + g(3) * (rSges(3,1) * t113 + rSges(3,2) * t116) + (-rSges(3,3) - pkin(8)) * t134) * t108) - m(4) * (g(1) * (t86 * pkin(2) + (t107 * t129 + t109 * t86) * rSges(4,1) + (-t107 * t86 + t109 * t129) * rSges(4,2) + t132 * t85 + t131) + g(2) * (t84 * pkin(2) + t103 - pkin(8) * t127 + (-t107 * t127 + t84 * t109) * rSges(4,1) + (-t84 * t107 - t109 * t127) * rSges(4,2) + t132 * t83) + g(3) * ((rSges(4,1) * t107 + rSges(4,2) * t109) * t110 + (-t132 * t116 + (rSges(4,1) * t109 - rSges(4,2) * t107 + pkin(2)) * t113) * t108 + t133)) - m(5) * ((t119 * t135 + (-pkin(8) - t119) * t134) * t108 + (t133 + t119 * t110 + (t120 * t113 - t116 * t138) * t108) * g(3) + (t120 * t84 + t138 * t83 + t103) * g(2) + (t120 * t86 + t138 * t85 + t131) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t74 - rSges(6,2) * t73 + rSges(6,3) * t85 + t121) + g(2) * (t72 * rSges(6,1) - t71 * rSges(6,2) + t83 * rSges(6,3) + t118) + g(3) * (rSges(6,1) * t78 - rSges(6,2) * t77 - rSges(6,3) * t128 + t122)) - m(7) * (g(1) * (t74 * pkin(5) + (t112 * t85 + t115 * t74) * rSges(7,1) + (-t112 * t74 + t115 * t85) * rSges(7,2) + t137 * t73 + t121) + g(2) * (t72 * pkin(5) + (t112 * t83 + t115 * t72) * rSges(7,1) + (-t112 * t72 + t115 * t83) * rSges(7,2) + t137 * t71 + t118) + g(3) * (t78 * pkin(5) + (-t112 * t128 + t115 * t78) * rSges(7,1) + (-t78 * t112 - t115 * t128) * rSges(7,2) + t137 * t77 + t122));
U  = t1;
