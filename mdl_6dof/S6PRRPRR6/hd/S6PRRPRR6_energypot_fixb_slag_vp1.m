% Calculate potential energy for
% S6PRRPRR6
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:44
% EndTime: 2019-03-08 22:22:45
% DurationCPUTime: 0.54s
% Computational Cost: add. (526->135), mult. (1191->191), div. (0->0), fcn. (1486->16), ass. (0->63)
t145 = sin(pkin(7));
t149 = cos(pkin(7));
t150 = cos(pkin(6));
t146 = sin(pkin(6));
t156 = cos(qJ(2));
t173 = t146 * t156;
t124 = -t145 * t173 + t150 * t149;
t144 = sin(pkin(12));
t148 = cos(pkin(12));
t154 = sin(qJ(2));
t170 = t150 * t156;
t127 = -t144 * t170 - t148 * t154;
t175 = t146 * t149;
t117 = -t127 * t145 + t144 * t175;
t184 = pkin(11) + rSges(7,3);
t183 = cos(qJ(3));
t182 = qJ(4) + rSges(5,3);
t125 = -t144 * t154 + t148 * t170;
t116 = -t125 * t145 - t148 * t175;
t143 = sin(pkin(13));
t181 = t116 * t143;
t180 = t117 * t143;
t179 = t124 * t143;
t177 = t144 * t146;
t176 = t146 * t148;
t174 = t146 * t154;
t171 = t150 * t154;
t169 = t150 * pkin(8) + qJ(1);
t168 = t148 * pkin(1) + pkin(8) * t177;
t165 = t145 * t183;
t164 = t149 * t183;
t163 = t146 * t165;
t128 = -t144 * t171 + t148 * t156;
t162 = t128 * pkin(2) + t117 * pkin(9) + t168;
t161 = pkin(2) * t174 + t124 * pkin(9) + t169;
t153 = sin(qJ(3));
t106 = -t127 * t164 + t128 * t153 - t144 * t163;
t107 = t128 * t183 + (t127 * t149 + t145 * t177) * t153;
t147 = cos(pkin(13));
t136 = pkin(4) * t147 + pkin(3);
t151 = -pkin(10) - qJ(4);
t160 = pkin(4) * t180 - t106 * t151 + t107 * t136 + t162;
t114 = -t150 * t165 + t153 * t174 - t164 * t173;
t115 = t150 * t145 * t153 + (t149 * t153 * t156 + t183 * t154) * t146;
t159 = pkin(4) * t179 - t114 * t151 + t115 * t136 + t161;
t126 = t144 * t156 + t148 * t171;
t139 = t144 * pkin(1);
t158 = t126 * pkin(2) - pkin(8) * t176 + t116 * pkin(9) + t139;
t104 = -t125 * t164 + t126 * t153 + t148 * t163;
t105 = t126 * t183 + (t125 * t149 - t145 * t176) * t153;
t157 = pkin(4) * t181 - t104 * t151 + t105 * t136 + t158;
t155 = cos(qJ(6));
t152 = sin(qJ(6));
t142 = pkin(13) + qJ(5);
t138 = cos(t142);
t137 = sin(t142);
t101 = t115 * t138 + t124 * t137;
t100 = t115 * t137 - t124 * t138;
t97 = t107 * t138 + t117 * t137;
t96 = t107 * t137 - t117 * t138;
t95 = t105 * t138 + t116 * t137;
t94 = t105 * t137 - t116 * t138;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t148 - rSges(2,2) * t144) + g(2) * (rSges(2,1) * t144 + rSges(2,2) * t148) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t128 + rSges(3,2) * t127 + t168) + g(2) * (rSges(3,1) * t126 + rSges(3,2) * t125 + t139) + g(3) * (t150 * rSges(3,3) + t169) + (g(1) * rSges(3,3) * t144 + g(3) * (rSges(3,1) * t154 + rSges(3,2) * t156) + g(2) * (-rSges(3,3) - pkin(8)) * t148) * t146) - m(4) * (g(1) * (rSges(4,1) * t107 - rSges(4,2) * t106 + rSges(4,3) * t117 + t162) + g(2) * (rSges(4,1) * t105 - rSges(4,2) * t104 + rSges(4,3) * t116 + t158) + g(3) * (t115 * rSges(4,1) - t114 * rSges(4,2) + t124 * rSges(4,3) + t161)) - m(5) * (g(1) * (t107 * pkin(3) + (t107 * t147 + t180) * rSges(5,1) + (-t107 * t143 + t117 * t147) * rSges(5,2) + t182 * t106 + t162) + g(2) * (t105 * pkin(3) + (t105 * t147 + t181) * rSges(5,1) + (-t105 * t143 + t116 * t147) * rSges(5,2) + t182 * t104 + t158) + g(3) * (t115 * pkin(3) + (t115 * t147 + t179) * rSges(5,1) + (-t115 * t143 + t124 * t147) * rSges(5,2) + t182 * t114 + t161)) - m(6) * (g(1) * (rSges(6,1) * t97 - rSges(6,2) * t96 + rSges(6,3) * t106 + t160) + g(2) * (rSges(6,1) * t95 - rSges(6,2) * t94 + rSges(6,3) * t104 + t157) + g(3) * (t101 * rSges(6,1) - t100 * rSges(6,2) + t114 * rSges(6,3) + t159)) - m(7) * (g(1) * (t97 * pkin(5) + (t106 * t152 + t155 * t97) * rSges(7,1) + (t106 * t155 - t152 * t97) * rSges(7,2) + t184 * t96 + t160) + g(2) * (t95 * pkin(5) + (t104 * t152 + t155 * t95) * rSges(7,1) + (t104 * t155 - t152 * t95) * rSges(7,2) + t184 * t94 + t157) + g(3) * (t101 * pkin(5) + (t101 * t155 + t114 * t152) * rSges(7,1) + (-t101 * t152 + t114 * t155) * rSges(7,2) + t184 * t100 + t159));
U  = t1;
