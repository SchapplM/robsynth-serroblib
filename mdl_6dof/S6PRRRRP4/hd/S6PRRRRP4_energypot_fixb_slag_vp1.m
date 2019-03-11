% Calculate potential energy for
% S6PRRRRP4
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:01
% EndTime: 2019-03-09 00:13:01
% DurationCPUTime: 0.37s
% Computational Cost: add. (346->120), mult. (660->157), div. (0->0), fcn. (793->12), ass. (0->56)
t150 = rSges(7,1) + pkin(5);
t149 = rSges(4,3) + pkin(8);
t148 = pkin(9) + rSges(5,3);
t118 = sin(pkin(6));
t147 = pkin(7) * t118;
t121 = sin(qJ(4));
t117 = sin(pkin(11));
t119 = cos(pkin(11));
t123 = sin(qJ(2));
t120 = cos(pkin(6));
t126 = cos(qJ(2));
t138 = t120 * t126;
t99 = t117 * t123 - t119 * t138;
t146 = t121 * t99;
t145 = rSges(7,3) + qJ(6);
t101 = t117 * t138 + t119 * t123;
t144 = t101 * t121;
t122 = sin(qJ(3));
t143 = t118 * t122;
t142 = t118 * t123;
t125 = cos(qJ(3));
t141 = t118 * t125;
t140 = t118 * t126;
t139 = t120 * t123;
t137 = t120 * pkin(7) + qJ(1);
t136 = t119 * pkin(1) + t117 * t147;
t102 = -t117 * t139 + t119 * t126;
t135 = t102 * pkin(2) + t136;
t134 = pkin(2) * t142 + t137;
t100 = t117 * t126 + t119 * t139;
t113 = t117 * pkin(1);
t133 = t100 * pkin(2) - t119 * t147 + t113;
t132 = t101 * pkin(8) + t135;
t131 = t99 * pkin(8) + t133;
t124 = cos(qJ(4));
t110 = pkin(4) * t124 + pkin(3);
t127 = -pkin(10) - pkin(9);
t89 = t102 * t122 - t117 * t141;
t90 = t102 * t125 + t117 * t143;
t130 = pkin(4) * t144 + t90 * t110 - t89 * t127 + t132;
t87 = t100 * t122 + t119 * t141;
t88 = t100 * t125 - t119 * t143;
t129 = pkin(4) * t146 + t88 * t110 - t87 * t127 + t131;
t103 = -t120 * t125 + t122 * t142;
t104 = t120 * t122 + t123 * t141;
t128 = t104 * t110 + (-pkin(4) * t121 - pkin(8)) * t140 - t103 * t127 + t134;
t116 = qJ(4) + qJ(5);
t112 = cos(t116);
t111 = sin(t116);
t84 = t104 * t112 - t111 * t140;
t83 = t104 * t111 + t112 * t140;
t80 = t101 * t111 + t112 * t90;
t79 = -t101 * t112 + t111 * t90;
t78 = t111 * t99 + t112 * t88;
t77 = t111 * t88 - t99 * t112;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t119 - rSges(2,2) * t117) + g(2) * (rSges(2,1) * t117 + rSges(2,2) * t119) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t102 - t101 * rSges(3,2) + t136) + g(2) * (rSges(3,1) * t100 - rSges(3,2) * t99 + t113) + g(3) * (rSges(3,3) * t120 + t137) + (g(1) * rSges(3,3) * t117 + g(3) * (rSges(3,1) * t123 + rSges(3,2) * t126) + g(2) * (-rSges(3,3) - pkin(7)) * t119) * t118) - m(4) * (g(1) * (rSges(4,1) * t90 - rSges(4,2) * t89 + t101 * t149 + t135) + g(2) * (rSges(4,1) * t88 - rSges(4,2) * t87 + t149 * t99 + t133) + g(3) * (rSges(4,1) * t104 - rSges(4,2) * t103 - t140 * t149 + t134)) - m(5) * (g(1) * (t90 * pkin(3) + (t124 * t90 + t144) * rSges(5,1) + (t101 * t124 - t121 * t90) * rSges(5,2) + t148 * t89 + t132) + g(2) * (t88 * pkin(3) + (t124 * t88 + t146) * rSges(5,1) + (-t121 * t88 + t124 * t99) * rSges(5,2) + t148 * t87 + t131) + g(3) * (t104 * pkin(3) - pkin(8) * t140 + (t104 * t124 - t121 * t140) * rSges(5,1) + (-t104 * t121 - t124 * t140) * rSges(5,2) + t148 * t103 + t134)) - m(6) * (g(1) * (rSges(6,1) * t80 - rSges(6,2) * t79 + rSges(6,3) * t89 + t130) + g(2) * (rSges(6,1) * t78 - rSges(6,2) * t77 + rSges(6,3) * t87 + t129) + g(3) * (rSges(6,1) * t84 - rSges(6,2) * t83 + rSges(6,3) * t103 + t128)) - m(7) * (g(1) * (rSges(7,2) * t89 + t145 * t79 + t150 * t80 + t130) + g(2) * (rSges(7,2) * t87 + t145 * t77 + t150 * t78 + t129) + g(3) * (rSges(7,2) * t103 + t145 * t83 + t150 * t84 + t128));
U  = t1;
