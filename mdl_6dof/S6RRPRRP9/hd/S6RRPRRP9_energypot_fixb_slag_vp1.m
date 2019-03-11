% Calculate potential energy for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:24
% EndTime: 2019-03-09 12:27:24
% DurationCPUTime: 0.38s
% Computational Cost: add. (343->123), mult. (545->158), div. (0->0), fcn. (633->12), ass. (0->55)
t140 = rSges(6,3) + pkin(10);
t108 = sin(pkin(11));
t139 = pkin(3) * t108;
t111 = cos(pkin(6));
t138 = t111 * pkin(8) + pkin(7);
t114 = sin(qJ(5));
t118 = cos(qJ(2));
t119 = cos(qJ(1));
t125 = t118 * t119;
t115 = sin(qJ(2));
t116 = sin(qJ(1));
t128 = t115 * t116;
t88 = -t111 * t125 + t128;
t137 = t114 * t88;
t126 = t116 * t118;
t127 = t115 * t119;
t90 = t111 * t126 + t127;
t136 = t114 * t90;
t135 = rSges(7,3) + qJ(6) + pkin(10);
t134 = qJ(3) + rSges(4,3);
t109 = sin(pkin(6));
t131 = t109 * t116;
t133 = t119 * pkin(1) + pkin(8) * t131;
t132 = t109 * t115;
t130 = t109 * t118;
t129 = t109 * t119;
t124 = t114 * t130;
t123 = t108 * t131;
t110 = cos(pkin(11));
t100 = pkin(3) * t110 + pkin(2);
t113 = -pkin(9) - qJ(3);
t122 = t100 * t132 + t111 * t139 + t113 * t130 + t138;
t91 = -t111 * t128 + t125;
t121 = pkin(3) * t123 + t91 * t100 - t90 * t113 + t133;
t105 = t116 * pkin(1);
t89 = t111 * t127 + t126;
t120 = t105 + t89 * t100 + (-pkin(8) - t139) * t129 - t88 * t113;
t117 = cos(qJ(5));
t107 = pkin(11) + qJ(4);
t103 = cos(t107);
t102 = sin(t107);
t101 = pkin(5) * t117 + pkin(4);
t85 = t111 * t102 + t103 * t132;
t84 = t102 * t132 - t111 * t103;
t81 = t102 * t131 + t91 * t103;
t80 = t102 * t91 - t103 * t131;
t79 = -t102 * t129 + t89 * t103;
t78 = t89 * t102 + t103 * t129;
t77 = t117 * t85 - t124;
t76 = -t114 * t85 - t117 * t130;
t75 = t117 * t81 + t136;
t74 = -t114 * t81 + t117 * t90;
t73 = t117 * t79 + t137;
t72 = -t114 * t79 + t117 * t88;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t119 - t116 * rSges(2,2)) + g(2) * (t116 * rSges(2,1) + rSges(2,2) * t119) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t91 * rSges(3,1) - t90 * rSges(3,2) + t133) + g(2) * (t89 * rSges(3,1) - t88 * rSges(3,2) + t105) + g(3) * (t111 * rSges(3,3) + t138) + (g(1) * rSges(3,3) * t116 + g(3) * (rSges(3,1) * t115 + rSges(3,2) * t118) + g(2) * (-rSges(3,3) - pkin(8)) * t119) * t109) - m(4) * (g(1) * (t91 * pkin(2) + (t91 * t110 + t123) * rSges(4,1) + (-t91 * t108 + t110 * t131) * rSges(4,2) + t134 * t90 + t133) + g(2) * (t89 * pkin(2) + t105 - pkin(8) * t129 + (-t108 * t129 + t89 * t110) * rSges(4,1) + (-t89 * t108 - t110 * t129) * rSges(4,2) + t134 * t88) + g(3) * ((t108 * rSges(4,1) + t110 * rSges(4,2)) * t111 + (-t134 * t118 + (t110 * rSges(4,1) - t108 * rSges(4,2) + pkin(2)) * t115) * t109 + t138)) - m(5) * (g(1) * (rSges(5,1) * t81 - rSges(5,2) * t80 + rSges(5,3) * t90 + t121) + g(2) * (t79 * rSges(5,1) - t78 * rSges(5,2) + t88 * rSges(5,3) + t120) + g(3) * (t85 * rSges(5,1) - t84 * rSges(5,2) - rSges(5,3) * t130 + t122)) - m(6) * (g(1) * (rSges(6,1) * t75 + rSges(6,2) * t74 + pkin(4) * t81 + t140 * t80 + t121) + g(2) * (t73 * rSges(6,1) + t72 * rSges(6,2) + t79 * pkin(4) + t140 * t78 + t120) + g(3) * (rSges(6,1) * t77 + rSges(6,2) * t76 + pkin(4) * t85 + t140 * t84 + t122)) - m(7) * (g(1) * (t75 * rSges(7,1) + t74 * rSges(7,2) + pkin(5) * t136 + t81 * t101 + t135 * t80 + t121) + g(2) * (t73 * rSges(7,1) + t72 * rSges(7,2) + pkin(5) * t137 + t79 * t101 + t135 * t78 + t120) + g(3) * (t77 * rSges(7,1) + t76 * rSges(7,2) - pkin(5) * t124 + t85 * t101 + t135 * t84 + t122));
U  = t1;
