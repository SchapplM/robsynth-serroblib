% Calculate potential energy for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:38
% EndTime: 2019-03-08 23:55:38
% DurationCPUTime: 0.38s
% Computational Cost: add. (343->123), mult. (545->162), div. (0->0), fcn. (633->12), ass. (0->55)
t140 = rSges(6,3) + pkin(10);
t139 = pkin(8) + rSges(4,3);
t114 = sin(qJ(3));
t138 = pkin(3) * t114;
t113 = sin(qJ(5));
t108 = sin(pkin(11));
t110 = cos(pkin(11));
t115 = sin(qJ(2));
t111 = cos(pkin(6));
t118 = cos(qJ(2));
t126 = t111 * t118;
t88 = t108 * t115 - t110 * t126;
t137 = t88 * t113;
t90 = t108 * t126 + t110 * t115;
t136 = t90 * t113;
t135 = rSges(7,3) + qJ(6) + pkin(10);
t109 = sin(pkin(6));
t133 = t108 * t109;
t134 = t110 * pkin(1) + pkin(7) * t133;
t132 = t109 * t110;
t131 = t109 * t114;
t130 = t109 * t115;
t117 = cos(qJ(3));
t129 = t109 * t117;
t128 = t109 * t118;
t127 = t111 * t115;
t125 = t111 * pkin(7) + qJ(1);
t124 = t108 * t131;
t123 = t113 * t128;
t101 = t117 * pkin(3) + pkin(2);
t119 = -pkin(9) - pkin(8);
t91 = -t108 * t127 + t110 * t118;
t122 = pkin(3) * t124 + t91 * t101 - t90 * t119 + t134;
t121 = t101 * t130 + t111 * t138 + t119 * t128 + t125;
t104 = t108 * pkin(1);
t89 = t108 * t118 + t110 * t127;
t120 = t104 + t89 * t101 + (-pkin(7) - t138) * t132 - t88 * t119;
t116 = cos(qJ(5));
t107 = qJ(3) + qJ(4);
t103 = cos(t107);
t102 = sin(t107);
t100 = t116 * pkin(5) + pkin(4);
t85 = t111 * t102 + t103 * t130;
t84 = t102 * t130 - t111 * t103;
t81 = t85 * t116 - t123;
t80 = -t85 * t113 - t116 * t128;
t79 = t102 * t133 + t91 * t103;
t78 = t91 * t102 - t103 * t133;
t77 = -t102 * t132 + t89 * t103;
t76 = t89 * t102 + t103 * t132;
t75 = t79 * t116 + t136;
t74 = -t79 * t113 + t90 * t116;
t73 = t77 * t116 + t137;
t72 = -t77 * t113 + t88 * t116;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t110 * rSges(2,1) - t108 * rSges(2,2)) + g(2) * (t108 * rSges(2,1) + t110 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t91 * rSges(3,1) - t90 * rSges(3,2) + t134) + g(2) * (t89 * rSges(3,1) - t88 * rSges(3,2) + t104) + g(3) * (t111 * rSges(3,3) + t125) + (g(1) * rSges(3,3) * t108 + g(3) * (rSges(3,1) * t115 + rSges(3,2) * t118) + g(2) * (-rSges(3,3) - pkin(7)) * t110) * t109) - m(4) * (g(1) * (t91 * pkin(2) + (t91 * t117 + t124) * rSges(4,1) + (t108 * t129 - t91 * t114) * rSges(4,2) + t139 * t90 + t134) + g(2) * (t89 * pkin(2) + t104 - pkin(7) * t132 + (-t110 * t131 + t89 * t117) * rSges(4,1) + (-t110 * t129 - t89 * t114) * rSges(4,2) + t139 * t88) + g(3) * ((t114 * rSges(4,1) + t117 * rSges(4,2)) * t111 + (-t139 * t118 + (t117 * rSges(4,1) - t114 * rSges(4,2) + pkin(2)) * t115) * t109 + t125)) - m(5) * (g(1) * (t79 * rSges(5,1) - t78 * rSges(5,2) + t90 * rSges(5,3) + t122) + g(2) * (t77 * rSges(5,1) - t76 * rSges(5,2) + t88 * rSges(5,3) + t120) + g(3) * (t85 * rSges(5,1) - t84 * rSges(5,2) - rSges(5,3) * t128 + t121)) - m(6) * (g(1) * (t75 * rSges(6,1) + t74 * rSges(6,2) + t79 * pkin(4) + t140 * t78 + t122) + g(2) * (t73 * rSges(6,1) + t72 * rSges(6,2) + t77 * pkin(4) + t140 * t76 + t120) + g(3) * (t81 * rSges(6,1) + t80 * rSges(6,2) + t85 * pkin(4) + t140 * t84 + t121)) - m(7) * (g(1) * (t75 * rSges(7,1) + t74 * rSges(7,2) + pkin(5) * t136 + t79 * t100 + t135 * t78 + t122) + g(2) * (t73 * rSges(7,1) + t72 * rSges(7,2) + pkin(5) * t137 + t77 * t100 + t135 * t76 + t120) + g(3) * (t81 * rSges(7,1) + t80 * rSges(7,2) - pkin(5) * t123 + t85 * t100 + t135 * t84 + t121));
U  = t1;
