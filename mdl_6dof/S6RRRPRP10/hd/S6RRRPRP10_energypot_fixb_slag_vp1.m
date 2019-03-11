% Calculate potential energy for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:22
% EndTime: 2019-03-09 17:29:22
% DurationCPUTime: 0.36s
% Computational Cost: add. (346->120), mult. (660->155), div. (0->0), fcn. (793->12), ass. (0->58)
t154 = rSges(7,1) + pkin(5);
t153 = rSges(4,3) + pkin(9);
t122 = cos(pkin(6));
t152 = t122 * pkin(8) + pkin(7);
t151 = rSges(7,3) + qJ(6);
t150 = qJ(4) + rSges(5,3);
t128 = cos(qJ(2));
t129 = cos(qJ(1));
t139 = t128 * t129;
t125 = sin(qJ(2));
t126 = sin(qJ(1));
t142 = t125 * t126;
t103 = -t122 * t139 + t142;
t119 = sin(pkin(11));
t149 = t103 * t119;
t140 = t126 * t128;
t141 = t125 * t129;
t105 = t122 * t140 + t141;
t148 = t105 * t119;
t120 = sin(pkin(6));
t147 = t120 * t125;
t146 = t120 * t126;
t127 = cos(qJ(3));
t145 = t120 * t127;
t144 = t120 * t128;
t143 = t120 * t129;
t138 = t129 * pkin(1) + pkin(8) * t146;
t137 = pkin(2) * t147 + t152;
t106 = -t122 * t142 + t139;
t136 = t106 * pkin(2) + t138;
t104 = t122 * t141 + t140;
t116 = t126 * pkin(1);
t135 = t104 * pkin(2) - pkin(8) * t143 + t116;
t134 = t105 * pkin(9) + t136;
t133 = t103 * pkin(9) + t135;
t121 = cos(pkin(11));
t112 = pkin(4) * t121 + pkin(3);
t123 = -pkin(10) - qJ(4);
t124 = sin(qJ(3));
t91 = t106 * t124 - t126 * t145;
t92 = t106 * t127 + t124 * t146;
t132 = pkin(4) * t148 + t92 * t112 - t91 * t123 + t134;
t89 = t104 * t124 + t127 * t143;
t90 = t104 * t127 - t124 * t143;
t131 = pkin(4) * t149 + t90 * t112 - t89 * t123 + t133;
t101 = -t122 * t127 + t124 * t147;
t102 = t122 * t124 + t125 * t145;
t130 = t102 * t112 + (-pkin(4) * t119 - pkin(9)) * t144 - t101 * t123 + t137;
t118 = pkin(11) + qJ(5);
t114 = cos(t118);
t113 = sin(t118);
t86 = t102 * t114 - t113 * t144;
t85 = t102 * t113 + t114 * t144;
t82 = t105 * t113 + t114 * t92;
t81 = -t105 * t114 + t113 * t92;
t80 = t103 * t113 + t114 * t90;
t79 = -t103 * t114 + t113 * t90;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t129 - t126 * rSges(2,2)) + g(2) * (t126 * rSges(2,1) + rSges(2,2) * t129) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t106 - rSges(3,2) * t105 + t138) + g(2) * (t104 * rSges(3,1) - t103 * rSges(3,2) + t116) + g(3) * (rSges(3,3) * t122 + t152) + (g(1) * rSges(3,3) * t126 + g(3) * (rSges(3,1) * t125 + rSges(3,2) * t128) + g(2) * (-rSges(3,3) - pkin(8)) * t129) * t120) - m(4) * (g(1) * (rSges(4,1) * t92 - rSges(4,2) * t91 + t153 * t105 + t136) + g(2) * (t90 * rSges(4,1) - t89 * rSges(4,2) + t153 * t103 + t135) + g(3) * (rSges(4,1) * t102 - rSges(4,2) * t101 - t153 * t144 + t137)) - m(5) * (g(1) * (t92 * pkin(3) + (t121 * t92 + t148) * rSges(5,1) + (t105 * t121 - t119 * t92) * rSges(5,2) + t150 * t91 + t134) + g(2) * (t90 * pkin(3) + (t121 * t90 + t149) * rSges(5,1) + (t103 * t121 - t119 * t90) * rSges(5,2) + t150 * t89 + t133) + g(3) * (t102 * pkin(3) - pkin(9) * t144 + (t102 * t121 - t119 * t144) * rSges(5,1) + (-t102 * t119 - t121 * t144) * rSges(5,2) + t150 * t101 + t137)) - m(6) * (g(1) * (rSges(6,1) * t82 - rSges(6,2) * t81 + rSges(6,3) * t91 + t132) + g(2) * (t80 * rSges(6,1) - t79 * rSges(6,2) + t89 * rSges(6,3) + t131) + g(3) * (rSges(6,1) * t86 - rSges(6,2) * t85 + rSges(6,3) * t101 + t130)) - m(7) * (g(1) * (rSges(7,2) * t91 + t151 * t81 + t154 * t82 + t132) + g(2) * (t89 * rSges(7,2) + t151 * t79 + t154 * t80 + t131) + g(3) * (rSges(7,2) * t101 + t151 * t85 + t154 * t86 + t130));
U  = t1;
