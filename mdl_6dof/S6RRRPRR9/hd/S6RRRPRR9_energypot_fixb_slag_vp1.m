% Calculate potential energy for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:26
% EndTime: 2019-03-09 18:56:27
% DurationCPUTime: 0.84s
% Computational Cost: add. (598->147), mult. (1464->210), div. (0->0), fcn. (1852->16), ass. (0->73)
t143 = sin(pkin(7));
t146 = cos(pkin(7));
t150 = sin(qJ(3));
t155 = cos(qJ(3));
t158 = (rSges(4,1) * t150 + rSges(4,2) * t155) * t143 + t146 * pkin(10);
t188 = -rSges(6,3) - pkin(11);
t187 = pkin(12) + rSges(7,3);
t152 = sin(qJ(1));
t186 = g(1) * t152;
t157 = cos(qJ(1));
t185 = g(2) * t157;
t147 = cos(pkin(6));
t183 = t147 * pkin(9) + pkin(8);
t182 = pkin(10) + qJ(4);
t156 = cos(qJ(2));
t170 = t156 * t157;
t151 = sin(qJ(2));
t173 = t152 * t151;
t130 = t147 * t170 - t173;
t181 = t130 * t143;
t172 = t152 * t156;
t174 = t151 * t157;
t132 = -t147 * t172 - t174;
t180 = t132 * t143;
t144 = sin(pkin(6));
t179 = t144 * t152;
t178 = t144 * t157;
t177 = t146 * t150;
t176 = t146 * t155;
t175 = t146 * t156;
t171 = t156 * t143;
t169 = t157 * pkin(1) + pkin(9) * t179;
t127 = pkin(3) * t143 * t150 + t182 * t146;
t128 = pkin(3) * t177 - t182 * t143;
t138 = pkin(3) * t155 + pkin(2);
t168 = t147 * t127 + t183 + (t128 * t156 + t138 * t151) * t144;
t133 = -t147 * t173 + t170;
t166 = t127 * t179 + t132 * t128 + t133 * t138 + t169;
t142 = sin(pkin(13));
t145 = cos(pkin(13));
t165 = t142 * t155 + t145 * t150;
t135 = -t142 * t150 + t145 * t155;
t124 = t165 * t143;
t126 = t165 * t146;
t110 = t124 * t147 + (t126 * t156 + t135 * t151) * t144;
t164 = t110 * pkin(4) + t168;
t107 = t124 * t179 + t126 * t132 + t133 * t135;
t163 = t107 * pkin(4) + t166;
t162 = t135 * t143;
t161 = t144 * t162;
t131 = t147 * t174 + t172;
t140 = t152 * pkin(1);
t160 = t130 * t128 + t131 * t138 + t140 + (-pkin(9) - t127) * t178;
t105 = -t124 * t178 + t130 * t126 + t131 * t135;
t159 = t105 * pkin(4) + t160;
t154 = cos(qJ(5));
t153 = cos(qJ(6));
t149 = sin(qJ(5));
t148 = sin(qJ(6));
t129 = -t144 * t171 + t146 * t147;
t125 = t135 * t146;
t116 = t146 * t179 - t180;
t115 = -t146 * t178 - t181;
t109 = (t125 * t156 - t151 * t165) * t144 + t147 * t162;
t106 = t132 * t125 - t133 * t165 + t152 * t161;
t104 = t125 * t130 - t131 * t165 - t157 * t161;
t103 = t110 * t154 + t129 * t149;
t102 = t110 * t149 - t129 * t154;
t99 = t107 * t154 + t116 * t149;
t98 = t107 * t149 - t116 * t154;
t97 = t105 * t154 + t115 * t149;
t96 = t105 * t149 - t115 * t154;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t157 - t152 * rSges(2,2)) + g(2) * (t152 * rSges(2,1) + rSges(2,2) * t157) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t133 + rSges(3,2) * t132 + t169) + g(2) * (t131 * rSges(3,1) + t130 * rSges(3,2) + t140) + g(3) * (rSges(3,3) * t147 + t183) + (rSges(3,3) * t186 + g(3) * (rSges(3,1) * t151 + rSges(3,2) * t156) + (-rSges(3,3) - pkin(9)) * t185) * t144) - m(4) * (g(1) * (t133 * pkin(2) - pkin(10) * t180 + (t132 * t177 + t133 * t155) * rSges(4,1) + (t132 * t176 - t133 * t150) * rSges(4,2) + t116 * rSges(4,3) + t169) + g(2) * (t131 * pkin(2) - pkin(10) * t181 + t140 + (t130 * t177 + t131 * t155) * rSges(4,1) + (t130 * t176 - t131 * t150) * rSges(4,2) + t115 * rSges(4,3)) + (t158 * t186 + (-pkin(9) - t158) * t185) * t144 + (t129 * rSges(4,3) + t183 + t158 * t147 + (t151 * pkin(2) - pkin(10) * t171 + (t150 * t175 + t151 * t155) * rSges(4,1) + (-t150 * t151 + t155 * t175) * rSges(4,2)) * t144) * g(3)) - m(5) * (g(1) * (rSges(5,1) * t107 + rSges(5,2) * t106 + rSges(5,3) * t116 + t166) + g(2) * (t105 * rSges(5,1) + t104 * rSges(5,2) + t115 * rSges(5,3) + t160) + g(3) * (rSges(5,1) * t110 + rSges(5,2) * t109 + rSges(5,3) * t129 + t168)) - m(6) * (g(1) * (rSges(6,1) * t99 - rSges(6,2) * t98 + t188 * t106 + t163) + g(2) * (t97 * rSges(6,1) - t96 * rSges(6,2) + t188 * t104 + t159) + g(3) * (rSges(6,1) * t103 - rSges(6,2) * t102 + t188 * t109 + t164)) - m(7) * (g(1) * (t99 * pkin(5) - t106 * pkin(11) + (-t106 * t148 + t153 * t99) * rSges(7,1) + (-t106 * t153 - t148 * t99) * rSges(7,2) + t187 * t98 + t163) + g(2) * (t97 * pkin(5) - t104 * pkin(11) + (-t104 * t148 + t153 * t97) * rSges(7,1) + (-t104 * t153 - t148 * t97) * rSges(7,2) + t187 * t96 + t159) + g(3) * (t103 * pkin(5) - t109 * pkin(11) + (t103 * t153 - t109 * t148) * rSges(7,1) + (-t103 * t148 - t109 * t153) * rSges(7,2) + t187 * t102 + t164));
U  = t1;
