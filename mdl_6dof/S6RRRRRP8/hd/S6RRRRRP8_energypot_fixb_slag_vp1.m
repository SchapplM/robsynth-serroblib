% Calculate potential energy for
% S6RRRRRP8
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:47:09
% EndTime: 2019-03-10 01:47:10
% DurationCPUTime: 0.37s
% Computational Cost: add. (370->118), mult. (591->151), div. (0->0), fcn. (697->12), ass. (0->56)
t150 = rSges(7,1) + pkin(5);
t149 = rSges(7,2) + pkin(11);
t148 = rSges(6,3) + pkin(11);
t147 = pkin(9) + rSges(4,3);
t120 = sin(qJ(3));
t146 = pkin(3) * t120;
t118 = cos(pkin(6));
t145 = t118 * pkin(8) + pkin(7);
t144 = rSges(7,3) + qJ(6);
t117 = sin(pkin(6));
t121 = sin(qJ(2));
t143 = t117 * t121;
t122 = sin(qJ(1));
t142 = t117 * t122;
t125 = cos(qJ(2));
t141 = t117 * t125;
t126 = cos(qJ(1));
t140 = t117 * t126;
t139 = t121 * t122;
t138 = t121 * t126;
t137 = t122 * t125;
t136 = t125 * t126;
t135 = t126 * pkin(1) + pkin(8) * t142;
t134 = t120 * t142;
t124 = cos(qJ(3));
t110 = pkin(3) * t124 + pkin(2);
t127 = -pkin(10) - pkin(9);
t133 = t110 * t143 + t118 * t146 + t127 * t141 + t145;
t100 = t118 * t137 + t138;
t101 = -t118 * t139 + t136;
t132 = pkin(3) * t134 - t100 * t127 + t101 * t110 + t135;
t116 = qJ(3) + qJ(4);
t111 = sin(t116);
t112 = cos(t116);
t93 = t111 * t118 + t112 * t143;
t131 = t93 * pkin(4) + t133;
t88 = t101 * t112 + t111 * t142;
t130 = t88 * pkin(4) + t132;
t114 = t122 * pkin(1);
t98 = -t118 * t136 + t139;
t99 = t118 * t138 + t137;
t129 = t114 + t99 * t110 + (-pkin(8) - t146) * t140 - t98 * t127;
t86 = -t111 * t140 + t99 * t112;
t128 = t86 * pkin(4) + t129;
t123 = cos(qJ(5));
t119 = sin(qJ(5));
t92 = t111 * t143 - t118 * t112;
t87 = t101 * t111 - t112 * t142;
t85 = t99 * t111 + t112 * t140;
t84 = -t119 * t141 + t123 * t93;
t83 = t119 * t93 + t123 * t141;
t80 = t100 * t119 + t123 * t88;
t79 = -t100 * t123 + t119 * t88;
t78 = t119 * t98 + t123 * t86;
t77 = t119 * t86 - t98 * t123;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t126 - rSges(2,2) * t122) + g(2) * (rSges(2,1) * t122 + rSges(2,2) * t126) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t101 * rSges(3,1) - t100 * rSges(3,2) + t135) + g(2) * (t99 * rSges(3,1) - t98 * rSges(3,2) + t114) + g(3) * (rSges(3,3) * t118 + t145) + (g(1) * rSges(3,3) * t122 + g(3) * (rSges(3,1) * t121 + rSges(3,2) * t125) + g(2) * (-rSges(3,3) - pkin(8)) * t126) * t117) - m(4) * (g(1) * (t101 * pkin(2) + (t101 * t124 + t134) * rSges(4,1) + (-t101 * t120 + t124 * t142) * rSges(4,2) + t147 * t100 + t135) + g(2) * (t99 * pkin(2) + t114 - pkin(8) * t140 + (-t120 * t140 + t124 * t99) * rSges(4,1) + (-t120 * t99 - t124 * t140) * rSges(4,2) + t147 * t98) + g(3) * ((t120 * rSges(4,1) + t124 * rSges(4,2)) * t118 + (-t147 * t125 + (t124 * rSges(4,1) - t120 * rSges(4,2) + pkin(2)) * t121) * t117 + t145)) - m(5) * (g(1) * (rSges(5,1) * t88 - rSges(5,2) * t87 + rSges(5,3) * t100 + t132) + g(2) * (t86 * rSges(5,1) - t85 * rSges(5,2) + t98 * rSges(5,3) + t129) + g(3) * (t93 * rSges(5,1) - t92 * rSges(5,2) - rSges(5,3) * t141 + t133)) - m(6) * (g(1) * (rSges(6,1) * t80 - rSges(6,2) * t79 + t148 * t87 + t130) + g(2) * (t78 * rSges(6,1) - t77 * rSges(6,2) + t148 * t85 + t128) + g(3) * (rSges(6,1) * t84 - rSges(6,2) * t83 + t148 * t92 + t131)) - m(7) * (g(1) * (t144 * t79 + t149 * t87 + t150 * t80 + t130) + g(2) * (t144 * t77 + t149 * t85 + t150 * t78 + t128) + g(3) * (t144 * t83 + t149 * t92 + t150 * t84 + t131));
U  = t1;
