% Calculate potential energy for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:23
% EndTime: 2019-03-10 01:34:23
% DurationCPUTime: 0.38s
% Computational Cost: add. (343->123), mult. (545->158), div. (0->0), fcn. (633->12), ass. (0->55)
t140 = rSges(6,3) + pkin(11);
t139 = pkin(9) + rSges(4,3);
t112 = sin(qJ(3));
t138 = pkin(3) * t112;
t109 = cos(pkin(6));
t137 = t109 * pkin(8) + pkin(7);
t111 = sin(qJ(5));
t117 = cos(qJ(2));
t118 = cos(qJ(1));
t125 = t118 * t117;
t113 = sin(qJ(2));
t114 = sin(qJ(1));
t128 = t114 * t113;
t88 = -t109 * t125 + t128;
t136 = t88 * t111;
t126 = t118 * t113;
t127 = t114 * t117;
t90 = t109 * t127 + t126;
t135 = t90 * t111;
t134 = rSges(7,3) + qJ(6) + pkin(11);
t108 = sin(pkin(6));
t131 = t108 * t114;
t133 = t118 * pkin(1) + pkin(8) * t131;
t132 = t108 * t113;
t130 = t108 * t117;
t129 = t108 * t118;
t124 = t111 * t130;
t123 = t112 * t131;
t116 = cos(qJ(3));
t101 = t116 * pkin(3) + pkin(2);
t119 = -pkin(10) - pkin(9);
t122 = t101 * t132 + t109 * t138 + t119 * t130 + t137;
t91 = -t109 * t128 + t125;
t121 = pkin(3) * t123 + t91 * t101 - t90 * t119 + t133;
t105 = t114 * pkin(1);
t89 = t109 * t126 + t127;
t120 = t105 + t89 * t101 + (-pkin(8) - t138) * t129 - t88 * t119;
t115 = cos(qJ(5));
t107 = qJ(3) + qJ(4);
t103 = cos(t107);
t102 = sin(t107);
t100 = t115 * pkin(5) + pkin(4);
t85 = t109 * t102 + t103 * t132;
t84 = t102 * t132 - t109 * t103;
t81 = t102 * t131 + t91 * t103;
t80 = t102 * t91 - t103 * t131;
t79 = -t102 * t129 + t89 * t103;
t78 = t89 * t102 + t103 * t129;
t77 = t85 * t115 - t124;
t76 = -t85 * t111 - t115 * t130;
t75 = t81 * t115 + t135;
t74 = -t81 * t111 + t90 * t115;
t73 = t79 * t115 + t136;
t72 = -t79 * t111 + t88 * t115;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t118 * rSges(2,1) - t114 * rSges(2,2)) + g(2) * (t114 * rSges(2,1) + t118 * rSges(2,2)) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t91 * rSges(3,1) - t90 * rSges(3,2) + t133) + g(2) * (t89 * rSges(3,1) - t88 * rSges(3,2) + t105) + g(3) * (t109 * rSges(3,3) + t137) + (g(1) * rSges(3,3) * t114 + g(3) * (rSges(3,1) * t113 + rSges(3,2) * t117) + g(2) * (-rSges(3,3) - pkin(8)) * t118) * t108) - m(4) * (g(1) * (t91 * pkin(2) + (t91 * t116 + t123) * rSges(4,1) + (-t91 * t112 + t116 * t131) * rSges(4,2) + t139 * t90 + t133) + g(2) * (t89 * pkin(2) + t105 - pkin(8) * t129 + (-t112 * t129 + t89 * t116) * rSges(4,1) + (-t89 * t112 - t116 * t129) * rSges(4,2) + t139 * t88) + g(3) * ((t112 * rSges(4,1) + t116 * rSges(4,2)) * t109 + (-t139 * t117 + (t116 * rSges(4,1) - t112 * rSges(4,2) + pkin(2)) * t113) * t108 + t137)) - m(5) * (g(1) * (t81 * rSges(5,1) - t80 * rSges(5,2) + t90 * rSges(5,3) + t121) + g(2) * (t79 * rSges(5,1) - t78 * rSges(5,2) + t88 * rSges(5,3) + t120) + g(3) * (t85 * rSges(5,1) - t84 * rSges(5,2) - rSges(5,3) * t130 + t122)) - m(6) * (g(1) * (t75 * rSges(6,1) + t74 * rSges(6,2) + t81 * pkin(4) + t140 * t80 + t121) + g(2) * (t73 * rSges(6,1) + t72 * rSges(6,2) + t79 * pkin(4) + t140 * t78 + t120) + g(3) * (t77 * rSges(6,1) + t76 * rSges(6,2) + t85 * pkin(4) + t140 * t84 + t122)) - m(7) * (g(1) * (t75 * rSges(7,1) + t74 * rSges(7,2) + pkin(5) * t135 + t81 * t100 + t134 * t80 + t121) + g(2) * (t73 * rSges(7,1) + t72 * rSges(7,2) + pkin(5) * t136 + t79 * t100 + t134 * t78 + t120) + g(3) * (t77 * rSges(7,1) + t76 * rSges(7,2) - pkin(5) * t124 + t85 * t100 + t134 * t84 + t122));
U  = t1;
