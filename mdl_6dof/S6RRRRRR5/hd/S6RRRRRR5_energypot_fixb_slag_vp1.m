% Calculate potential energy for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:32
% EndTime: 2019-03-10 03:56:33
% DurationCPUTime: 0.58s
% Computational Cost: add. (357->125), mult. (471->159), div. (0->0), fcn. (534->14), ass. (0->54)
t117 = -pkin(10) - pkin(9);
t138 = -t117 + rSges(5,3);
t137 = pkin(9) + rSges(4,3);
t136 = pkin(12) + rSges(7,3);
t112 = sin(qJ(1));
t135 = g(1) * t112;
t116 = cos(qJ(1));
t134 = g(2) * t116;
t110 = sin(qJ(3));
t133 = t110 * pkin(3);
t108 = cos(pkin(6));
t132 = t108 * pkin(8) + pkin(7);
t114 = cos(qJ(3));
t97 = t114 * pkin(3) + pkin(2);
t107 = sin(pkin(6));
t129 = t107 * t112;
t131 = t116 * pkin(1) + pkin(8) * t129;
t111 = sin(qJ(2));
t130 = t107 * t111;
t115 = cos(qJ(2));
t128 = t107 * t115;
t127 = t107 * t116;
t126 = t112 * t111;
t125 = t112 * t115;
t124 = t116 * t111;
t123 = t116 * t115;
t106 = qJ(3) + qJ(4);
t105 = -pkin(11) + t117;
t99 = cos(t106);
t89 = pkin(4) * t99 + t97;
t98 = sin(t106);
t90 = pkin(4) * t98 + t133;
t122 = t105 * t128 + t108 * t90 + t89 * t130 + t132;
t85 = t108 * t125 + t124;
t86 = -t108 * t126 + t123;
t121 = -t85 * t105 + t90 * t129 + t86 * t89 + t131;
t120 = t99 * rSges(5,1) - t98 * rSges(5,2) + t97;
t119 = t98 * rSges(5,1) + t99 * rSges(5,2) + t133;
t102 = t112 * pkin(1);
t83 = -t108 * t123 + t126;
t84 = t108 * t124 + t125;
t118 = t102 + t84 * t89 + (-pkin(8) - t90) * t127 - t83 * t105;
t113 = cos(qJ(6));
t109 = sin(qJ(6));
t101 = qJ(5) + t106;
t96 = cos(t101);
t95 = sin(t101);
t78 = t108 * t95 + t130 * t96;
t77 = -t108 * t96 + t130 * t95;
t74 = t129 * t95 + t86 * t96;
t73 = -t129 * t96 + t86 * t95;
t72 = -t127 * t95 + t84 * t96;
t71 = t127 * t96 + t84 * t95;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t116 * rSges(2,1) - t112 * rSges(2,2)) + g(2) * (t112 * rSges(2,1) + t116 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t86 * rSges(3,1) - t85 * rSges(3,2) + t131) + g(2) * (t84 * rSges(3,1) - t83 * rSges(3,2) + t102) + g(3) * (t108 * rSges(3,3) + t132) + (rSges(3,3) * t135 + g(3) * (rSges(3,1) * t111 + rSges(3,2) * t115) + (-rSges(3,3) - pkin(8)) * t134) * t107) - m(4) * (g(1) * (t86 * pkin(2) + (t110 * t129 + t86 * t114) * rSges(4,1) + (-t86 * t110 + t114 * t129) * rSges(4,2) + t137 * t85 + t131) + g(2) * (t84 * pkin(2) + t102 - pkin(8) * t127 + (-t110 * t127 + t84 * t114) * rSges(4,1) + (-t84 * t110 - t114 * t127) * rSges(4,2) + t137 * t83) + g(3) * ((t110 * rSges(4,1) + t114 * rSges(4,2)) * t108 + (-t137 * t115 + (t114 * rSges(4,1) - t110 * rSges(4,2) + pkin(2)) * t111) * t107 + t132)) - m(5) * ((t119 * t135 + (-pkin(8) - t119) * t134) * t107 + (t132 + t119 * t108 + (t120 * t111 - t138 * t115) * t107) * g(3) + (t120 * t84 + t138 * t83 + t102) * g(2) + (t120 * t86 + t138 * t85 + t131) * g(1)) - m(6) * (g(1) * (t74 * rSges(6,1) - t73 * rSges(6,2) + t85 * rSges(6,3) + t121) + g(2) * (t72 * rSges(6,1) - t71 * rSges(6,2) + t83 * rSges(6,3) + t118) + g(3) * (t78 * rSges(6,1) - t77 * rSges(6,2) - rSges(6,3) * t128 + t122)) - m(7) * (g(1) * (t74 * pkin(5) + (t85 * t109 + t74 * t113) * rSges(7,1) + (-t74 * t109 + t85 * t113) * rSges(7,2) + t136 * t73 + t121) + g(2) * (t72 * pkin(5) + (t83 * t109 + t72 * t113) * rSges(7,1) + (-t72 * t109 + t83 * t113) * rSges(7,2) + t136 * t71 + t118) + g(3) * (t78 * pkin(5) + (-t109 * t128 + t78 * t113) * rSges(7,1) + (-t78 * t109 - t113 * t128) * rSges(7,2) + t136 * t77 + t122));
U  = t1;
