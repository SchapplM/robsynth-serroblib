% Calculate potential energy for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:08
% EndTime: 2019-03-09 04:01:09
% DurationCPUTime: 0.84s
% Computational Cost: add. (598->147), mult. (1464->210), div. (0->0), fcn. (1852->16), ass. (0->71)
t142 = sin(pkin(7));
t146 = cos(pkin(7));
t150 = sin(qJ(3));
t154 = cos(qJ(3));
t156 = t146 * pkin(9) + (rSges(4,1) * t150 + rSges(4,2) * t154) * t142;
t184 = -rSges(6,3) - pkin(10);
t183 = pkin(11) + rSges(7,3);
t151 = sin(qJ(1));
t182 = g(1) * t151;
t155 = cos(qJ(1));
t181 = g(2) * t155;
t147 = cos(pkin(6));
t179 = t147 * qJ(2) + pkin(8);
t178 = pkin(9) + qJ(4);
t145 = cos(pkin(12));
t141 = sin(pkin(12));
t169 = t151 * t141;
t170 = t147 * t155;
t128 = t145 * t170 - t169;
t177 = t128 * t142;
t168 = t151 * t145;
t130 = -t141 * t155 - t147 * t168;
t176 = t130 * t142;
t143 = sin(pkin(6));
t175 = t143 * t151;
t174 = t143 * t155;
t173 = t145 * t142;
t172 = t146 * t150;
t171 = t146 * t154;
t167 = t155 * pkin(1) + qJ(2) * t175;
t125 = pkin(3) * t142 * t150 + t178 * t146;
t126 = pkin(3) * t172 - t178 * t142;
t136 = pkin(3) * t154 + pkin(2);
t166 = t147 * t125 + t179 + (t126 * t145 + t136 * t141) * t143;
t131 = t145 * t155 - t147 * t169;
t164 = t125 * t175 + t130 * t126 + t131 * t136 + t167;
t140 = sin(pkin(13));
t144 = cos(pkin(13));
t163 = t140 * t154 + t144 * t150;
t133 = -t140 * t150 + t144 * t154;
t122 = t163 * t142;
t124 = t163 * t146;
t108 = t122 * t147 + (t124 * t145 + t133 * t141) * t143;
t162 = t108 * pkin(4) + t166;
t105 = t122 * t175 + t124 * t130 + t131 * t133;
t161 = t105 * pkin(4) + t164;
t160 = t133 * t142;
t159 = t143 * t160;
t129 = t141 * t170 + t168;
t138 = t151 * pkin(1);
t158 = t128 * t126 + t129 * t136 + t138 + (-qJ(2) - t125) * t174;
t103 = -t122 * t174 + t128 * t124 + t129 * t133;
t157 = t103 * pkin(4) + t158;
t153 = cos(qJ(5));
t152 = cos(qJ(6));
t149 = sin(qJ(5));
t148 = sin(qJ(6));
t127 = -t143 * t173 + t146 * t147;
t123 = t133 * t146;
t114 = t146 * t175 - t176;
t113 = -t146 * t174 - t177;
t107 = (t123 * t145 - t141 * t163) * t143 + t147 * t160;
t104 = t130 * t123 - t131 * t163 + t151 * t159;
t102 = t123 * t128 - t129 * t163 - t155 * t159;
t101 = t108 * t153 + t127 * t149;
t100 = t108 * t149 - t127 * t153;
t97 = t105 * t153 + t114 * t149;
t96 = t105 * t149 - t114 * t153;
t95 = t103 * t153 + t113 * t149;
t94 = t103 * t149 - t113 * t153;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t155 - t151 * rSges(2,2)) + g(2) * (t151 * rSges(2,1) + rSges(2,2) * t155) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t131 + rSges(3,2) * t130 + t167) + g(2) * (t129 * rSges(3,1) + t128 * rSges(3,2) + t138) + g(3) * (rSges(3,3) * t147 + t179) + (rSges(3,3) * t182 + g(3) * (rSges(3,1) * t141 + rSges(3,2) * t145) + (-rSges(3,3) - qJ(2)) * t181) * t143) - m(4) * (g(1) * (t131 * pkin(2) - pkin(9) * t176 + (t130 * t172 + t131 * t154) * rSges(4,1) + (t130 * t171 - t131 * t150) * rSges(4,2) + t114 * rSges(4,3) + t167) + g(2) * (t129 * pkin(2) - pkin(9) * t177 + t138 + (t128 * t172 + t129 * t154) * rSges(4,1) + (t128 * t171 - t129 * t150) * rSges(4,2) + t113 * rSges(4,3)) + (t156 * t182 + (-qJ(2) - t156) * t181) * t143 + (t127 * rSges(4,3) + t179 + t156 * t147 + (t141 * pkin(2) - pkin(9) * t173 + (t141 * t154 + t145 * t172) * rSges(4,1) + (-t141 * t150 + t145 * t171) * rSges(4,2)) * t143) * g(3)) - m(5) * (g(1) * (rSges(5,1) * t105 + rSges(5,2) * t104 + rSges(5,3) * t114 + t164) + g(2) * (t103 * rSges(5,1) + t102 * rSges(5,2) + t113 * rSges(5,3) + t158) + g(3) * (rSges(5,1) * t108 + rSges(5,2) * t107 + rSges(5,3) * t127 + t166)) - m(6) * (g(1) * (rSges(6,1) * t97 - rSges(6,2) * t96 + t184 * t104 + t161) + g(2) * (t95 * rSges(6,1) - t94 * rSges(6,2) + t184 * t102 + t157) + g(3) * (rSges(6,1) * t101 - rSges(6,2) * t100 + t184 * t107 + t162)) - m(7) * (g(1) * (t97 * pkin(5) - t104 * pkin(10) + (-t104 * t148 + t152 * t97) * rSges(7,1) + (-t104 * t152 - t148 * t97) * rSges(7,2) + t183 * t96 + t161) + g(2) * (t95 * pkin(5) - t102 * pkin(10) + (-t102 * t148 + t152 * t95) * rSges(7,1) + (-t102 * t152 - t148 * t95) * rSges(7,2) + t183 * t94 + t157) + g(3) * (t101 * pkin(5) - t107 * pkin(10) + (t101 * t152 - t107 * t148) * rSges(7,1) + (-t101 * t148 - t107 * t152) * rSges(7,2) + t183 * t100 + t162));
U  = t1;
