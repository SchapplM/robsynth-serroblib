% Calculate potential energy for
% S6PRRRRP2
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:03
% EndTime: 2019-03-09 00:01:04
% DurationCPUTime: 0.37s
% Computational Cost: add. (370->118), mult. (591->155), div. (0->0), fcn. (697->12), ass. (0->56)
t150 = rSges(7,1) + pkin(5);
t149 = rSges(7,2) + pkin(10);
t148 = rSges(6,3) + pkin(10);
t147 = pkin(8) + rSges(4,3);
t122 = sin(qJ(3));
t146 = pkin(3) * t122;
t145 = rSges(7,3) + qJ(6);
t117 = sin(pkin(11));
t118 = sin(pkin(6));
t144 = t117 * t118;
t119 = cos(pkin(11));
t143 = t118 * t119;
t142 = t118 * t122;
t123 = sin(qJ(2));
t141 = t118 * t123;
t125 = cos(qJ(3));
t140 = t118 * t125;
t126 = cos(qJ(2));
t139 = t118 * t126;
t120 = cos(pkin(6));
t138 = t120 * t123;
t137 = t120 * t126;
t136 = t120 * pkin(7) + qJ(1);
t135 = t119 * pkin(1) + pkin(7) * t144;
t134 = t117 * t142;
t100 = t117 * t137 + t119 * t123;
t101 = -t117 * t138 + t119 * t126;
t110 = pkin(3) * t125 + pkin(2);
t127 = -pkin(9) - pkin(8);
t133 = pkin(3) * t134 - t100 * t127 + t101 * t110 + t135;
t132 = t110 * t141 + t120 * t146 + t127 * t139 + t136;
t116 = qJ(3) + qJ(4);
t111 = sin(t116);
t112 = cos(t116);
t86 = t101 * t112 + t111 * t144;
t131 = t86 * pkin(4) + t133;
t93 = t111 * t120 + t112 * t141;
t130 = t93 * pkin(4) + t132;
t113 = t117 * pkin(1);
t98 = t117 * t123 - t119 * t137;
t99 = t117 * t126 + t119 * t138;
t129 = t113 + t99 * t110 + (-pkin(7) - t146) * t143 - t98 * t127;
t84 = -t111 * t143 + t112 * t99;
t128 = t84 * pkin(4) + t129;
t124 = cos(qJ(5));
t121 = sin(qJ(5));
t92 = t111 * t141 - t120 * t112;
t88 = -t121 * t139 + t124 * t93;
t87 = t121 * t93 + t124 * t139;
t85 = t101 * t111 - t112 * t144;
t83 = t111 * t99 + t112 * t143;
t80 = t100 * t121 + t124 * t86;
t79 = -t100 * t124 + t121 * t86;
t78 = t121 * t98 + t124 * t84;
t77 = t121 * t84 - t98 * t124;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t119 - rSges(2,2) * t117) + g(2) * (rSges(2,1) * t117 + rSges(2,2) * t119) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t101 - rSges(3,2) * t100 + t135) + g(2) * (rSges(3,1) * t99 - rSges(3,2) * t98 + t113) + g(3) * (rSges(3,3) * t120 + t136) + (g(1) * rSges(3,3) * t117 + g(3) * (rSges(3,1) * t123 + rSges(3,2) * t126) + g(2) * (-rSges(3,3) - pkin(7)) * t119) * t118) - m(4) * (g(1) * (t101 * pkin(2) + (t101 * t125 + t134) * rSges(4,1) + (-t101 * t122 + t117 * t140) * rSges(4,2) + t147 * t100 + t135) + g(2) * (t99 * pkin(2) + t113 - pkin(7) * t143 + (-t119 * t142 + t125 * t99) * rSges(4,1) + (-t119 * t140 - t122 * t99) * rSges(4,2) + t147 * t98) + g(3) * ((t122 * rSges(4,1) + t125 * rSges(4,2)) * t120 + (-t147 * t126 + (t125 * rSges(4,1) - t122 * rSges(4,2) + pkin(2)) * t123) * t118 + t136)) - m(5) * (g(1) * (rSges(5,1) * t86 - rSges(5,2) * t85 + rSges(5,3) * t100 + t133) + g(2) * (rSges(5,1) * t84 - rSges(5,2) * t83 + rSges(5,3) * t98 + t129) + g(3) * (rSges(5,1) * t93 - rSges(5,2) * t92 - rSges(5,3) * t139 + t132)) - m(6) * (g(1) * (rSges(6,1) * t80 - rSges(6,2) * t79 + t148 * t85 + t131) + g(2) * (rSges(6,1) * t78 - rSges(6,2) * t77 + t148 * t83 + t128) + g(3) * (rSges(6,1) * t88 - rSges(6,2) * t87 + t148 * t92 + t130)) - m(7) * (g(1) * (t145 * t79 + t149 * t85 + t150 * t80 + t131) + g(2) * (t145 * t77 + t149 * t83 + t150 * t78 + t128) + g(3) * (t145 * t87 + t149 * t92 + t150 * t88 + t130));
U  = t1;
