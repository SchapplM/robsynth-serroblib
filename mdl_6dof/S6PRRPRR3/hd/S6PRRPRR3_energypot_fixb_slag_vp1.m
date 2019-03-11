% Calculate potential energy for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:07
% EndTime: 2019-03-08 22:02:08
% DurationCPUTime: 0.82s
% Computational Cost: add. (598->147), mult. (1464->213), div. (0->0), fcn. (1852->16), ass. (0->72)
t143 = sin(pkin(7));
t147 = cos(pkin(7));
t151 = sin(qJ(3));
t155 = cos(qJ(3));
t157 = t147 * pkin(9) + (rSges(4,1) * t151 + rSges(4,2) * t155) * t143;
t186 = -rSges(6,3) - pkin(10);
t185 = pkin(11) + rSges(7,3);
t142 = sin(pkin(12));
t184 = g(1) * t142;
t146 = cos(pkin(12));
t183 = g(2) * t146;
t181 = pkin(9) + qJ(4);
t152 = sin(qJ(2));
t148 = cos(pkin(6));
t156 = cos(qJ(2));
t171 = t148 * t156;
t129 = -t142 * t152 + t146 * t171;
t180 = t129 * t143;
t131 = -t142 * t171 - t146 * t152;
t179 = t131 * t143;
t144 = sin(pkin(6));
t178 = t142 * t144;
t177 = t144 * t146;
t176 = t144 * t147;
t175 = t147 * t151;
t174 = t147 * t155;
t173 = t147 * t156;
t172 = t148 * t152;
t170 = t156 * t143;
t169 = t148 * pkin(8) + qJ(1);
t168 = t146 * pkin(1) + pkin(8) * t178;
t126 = pkin(3) * t143 * t151 + t181 * t147;
t127 = pkin(3) * t175 - t181 * t143;
t137 = pkin(3) * t155 + pkin(2);
t166 = t148 * t126 + t169 + (t127 * t156 + t137 * t152) * t144;
t132 = -t142 * t172 + t146 * t156;
t165 = t126 * t178 + t131 * t127 + t132 * t137 + t168;
t141 = sin(pkin(13));
t145 = cos(pkin(13));
t164 = t141 * t155 + t145 * t151;
t134 = -t141 * t151 + t145 * t155;
t123 = t164 * t143;
t125 = t164 * t147;
t109 = t148 * t123 + (t125 * t156 + t134 * t152) * t144;
t163 = t109 * pkin(4) + t166;
t104 = t123 * t178 + t125 * t131 + t132 * t134;
t162 = t104 * pkin(4) + t165;
t161 = t134 * t143;
t160 = t144 * t161;
t130 = t142 * t156 + t146 * t172;
t138 = t142 * pkin(1);
t159 = t129 * t127 + t130 * t137 + t138 + (-pkin(8) - t126) * t177;
t102 = -t123 * t177 + t125 * t129 + t130 * t134;
t158 = t102 * pkin(4) + t159;
t154 = cos(qJ(5));
t153 = cos(qJ(6));
t150 = sin(qJ(5));
t149 = sin(qJ(6));
t128 = -t144 * t170 + t148 * t147;
t124 = t134 * t147;
t115 = t142 * t176 - t179;
t114 = -t146 * t176 - t180;
t108 = (t124 * t156 - t152 * t164) * t144 + t148 * t161;
t106 = t109 * t154 + t128 * t150;
t105 = t109 * t150 - t128 * t154;
t103 = t131 * t124 - t132 * t164 + t142 * t160;
t101 = t124 * t129 - t130 * t164 - t146 * t160;
t98 = t104 * t154 + t115 * t150;
t97 = t104 * t150 - t115 * t154;
t96 = t102 * t154 + t114 * t150;
t95 = t102 * t150 - t114 * t154;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t146 - rSges(2,2) * t142) + g(2) * (rSges(2,1) * t142 + rSges(2,2) * t146) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t132 + rSges(3,2) * t131 + t168) + g(2) * (rSges(3,1) * t130 + rSges(3,2) * t129 + t138) + g(3) * (t148 * rSges(3,3) + t169) + (rSges(3,3) * t184 + g(3) * (rSges(3,1) * t152 + rSges(3,2) * t156) + (-rSges(3,3) - pkin(8)) * t183) * t144) - m(4) * (g(1) * (t132 * pkin(2) - pkin(9) * t179 + (t131 * t175 + t132 * t155) * rSges(4,1) + (t131 * t174 - t132 * t151) * rSges(4,2) + t115 * rSges(4,3) + t168) + g(2) * (t130 * pkin(2) - pkin(9) * t180 + t138 + (t129 * t175 + t130 * t155) * rSges(4,1) + (t129 * t174 - t130 * t151) * rSges(4,2) + t114 * rSges(4,3)) + (t157 * t184 + (-pkin(8) - t157) * t183) * t144 + (t128 * rSges(4,3) + t169 + t157 * t148 + (t152 * pkin(2) - pkin(9) * t170 + (t151 * t173 + t152 * t155) * rSges(4,1) + (-t151 * t152 + t155 * t173) * rSges(4,2)) * t144) * g(3)) - m(5) * (g(1) * (rSges(5,1) * t104 + rSges(5,2) * t103 + rSges(5,3) * t115 + t165) + g(2) * (rSges(5,1) * t102 + rSges(5,2) * t101 + rSges(5,3) * t114 + t159) + g(3) * (rSges(5,1) * t109 + rSges(5,2) * t108 + rSges(5,3) * t128 + t166)) - m(6) * (g(1) * (rSges(6,1) * t98 - rSges(6,2) * t97 + t186 * t103 + t162) + g(2) * (rSges(6,1) * t96 - rSges(6,2) * t95 + t186 * t101 + t158) + g(3) * (rSges(6,1) * t106 - rSges(6,2) * t105 + t186 * t108 + t163)) - m(7) * (g(1) * (t98 * pkin(5) - t103 * pkin(10) + (-t103 * t149 + t153 * t98) * rSges(7,1) + (-t103 * t153 - t149 * t98) * rSges(7,2) + t185 * t97 + t162) + g(2) * (t96 * pkin(5) - t101 * pkin(10) + (-t101 * t149 + t153 * t96) * rSges(7,1) + (-t101 * t153 - t149 * t96) * rSges(7,2) + t185 * t95 + t158) + g(3) * (t106 * pkin(5) - t108 * pkin(10) + (t106 * t153 - t108 * t149) * rSges(7,1) + (-t106 * t149 - t108 * t153) * rSges(7,2) + t185 * t105 + t163));
U  = t1;
