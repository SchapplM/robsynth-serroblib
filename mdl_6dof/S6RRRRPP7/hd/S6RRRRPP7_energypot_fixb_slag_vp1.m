% Calculate potential energy for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:02
% EndTime: 2019-03-09 21:17:02
% DurationCPUTime: 0.36s
% Computational Cost: add. (346->120), mult. (660->155), div. (0->0), fcn. (793->12), ass. (0->58)
t154 = rSges(7,1) + pkin(5);
t153 = rSges(4,3) + pkin(9);
t152 = pkin(10) + rSges(5,3);
t120 = cos(pkin(6));
t151 = t120 * pkin(8) + pkin(7);
t150 = rSges(7,3) + qJ(6);
t128 = cos(qJ(2));
t129 = cos(qJ(1));
t139 = t128 * t129;
t124 = sin(qJ(2));
t125 = sin(qJ(1));
t142 = t124 * t125;
t103 = -t120 * t139 + t142;
t122 = sin(qJ(4));
t149 = t103 * t122;
t140 = t125 * t128;
t141 = t124 * t129;
t105 = t120 * t140 + t141;
t148 = t105 * t122;
t119 = sin(pkin(6));
t147 = t119 * t124;
t146 = t119 * t125;
t127 = cos(qJ(3));
t145 = t119 * t127;
t144 = t119 * t128;
t143 = t119 * t129;
t138 = t129 * pkin(1) + pkin(8) * t146;
t137 = pkin(2) * t147 + t151;
t106 = -t120 * t142 + t139;
t136 = t106 * pkin(2) + t138;
t104 = t120 * t141 + t140;
t116 = t125 * pkin(1);
t135 = t104 * pkin(2) - pkin(8) * t143 + t116;
t134 = t105 * pkin(9) + t136;
t133 = t103 * pkin(9) + t135;
t126 = cos(qJ(4));
t112 = pkin(4) * t126 + pkin(3);
t121 = -qJ(5) - pkin(10);
t123 = sin(qJ(3));
t91 = t106 * t123 - t125 * t145;
t92 = t106 * t127 + t123 * t146;
t132 = pkin(4) * t148 + t92 * t112 - t91 * t121 + t134;
t89 = t104 * t123 + t127 * t143;
t90 = t104 * t127 - t123 * t143;
t131 = pkin(4) * t149 + t90 * t112 - t89 * t121 + t133;
t101 = -t120 * t127 + t123 * t147;
t102 = t120 * t123 + t124 * t145;
t130 = t102 * t112 + (-pkin(4) * t122 - pkin(9)) * t144 - t101 * t121 + t137;
t118 = qJ(4) + pkin(11);
t114 = cos(t118);
t113 = sin(t118);
t86 = t102 * t114 - t113 * t144;
t85 = t102 * t113 + t114 * t144;
t82 = t105 * t113 + t114 * t92;
t81 = -t105 * t114 + t113 * t92;
t80 = t103 * t113 + t114 * t90;
t79 = -t103 * t114 + t113 * t90;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t129 - t125 * rSges(2,2)) + g(2) * (t125 * rSges(2,1) + rSges(2,2) * t129) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t106 - rSges(3,2) * t105 + t138) + g(2) * (t104 * rSges(3,1) - t103 * rSges(3,2) + t116) + g(3) * (rSges(3,3) * t120 + t151) + (g(1) * rSges(3,3) * t125 + g(3) * (rSges(3,1) * t124 + rSges(3,2) * t128) + g(2) * (-rSges(3,3) - pkin(8)) * t129) * t119) - m(4) * (g(1) * (rSges(4,1) * t92 - rSges(4,2) * t91 + t105 * t153 + t136) + g(2) * (t90 * rSges(4,1) - t89 * rSges(4,2) + t103 * t153 + t135) + g(3) * (rSges(4,1) * t102 - rSges(4,2) * t101 - t144 * t153 + t137)) - m(5) * (g(1) * (t92 * pkin(3) + (t126 * t92 + t148) * rSges(5,1) + (t105 * t126 - t122 * t92) * rSges(5,2) + t152 * t91 + t134) + g(2) * (t90 * pkin(3) + (t126 * t90 + t149) * rSges(5,1) + (t103 * t126 - t122 * t90) * rSges(5,2) + t152 * t89 + t133) + g(3) * (t102 * pkin(3) - pkin(9) * t144 + (t102 * t126 - t122 * t144) * rSges(5,1) + (-t102 * t122 - t126 * t144) * rSges(5,2) + t152 * t101 + t137)) - m(6) * (g(1) * (rSges(6,1) * t82 - rSges(6,2) * t81 + rSges(6,3) * t91 + t132) + g(2) * (t80 * rSges(6,1) - t79 * rSges(6,2) + t89 * rSges(6,3) + t131) + g(3) * (rSges(6,1) * t86 - rSges(6,2) * t85 + rSges(6,3) * t101 + t130)) - m(7) * (g(1) * (rSges(7,2) * t91 + t150 * t81 + t154 * t82 + t132) + g(2) * (t89 * rSges(7,2) + t150 * t79 + t154 * t80 + t131) + g(3) * (rSges(7,2) * t101 + t150 * t85 + t154 * t86 + t130));
U  = t1;
