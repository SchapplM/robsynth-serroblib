% Calculate potential energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRPRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRPRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:41:55
% EndTime: 2019-03-08 18:41:56
% DurationCPUTime: 0.83s
% Computational Cost: add. (598->147), mult. (1464->212), div. (0->0), fcn. (1852->16), ass. (0->71)
t143 = sin(pkin(7));
t148 = cos(pkin(7));
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t156 = (rSges(4,1) * t152 + rSges(4,2) * t155) * t143 + t148 * pkin(8);
t184 = -rSges(6,3) - pkin(9);
t183 = pkin(10) + rSges(7,3);
t142 = sin(pkin(11));
t182 = g(1) * t142;
t147 = cos(pkin(11));
t181 = g(2) * t147;
t179 = pkin(8) + qJ(4);
t141 = sin(pkin(12));
t146 = cos(pkin(12));
t149 = cos(pkin(6));
t171 = t147 * t149;
t128 = -t141 * t142 + t146 * t171;
t178 = t128 * t143;
t175 = t142 * t149;
t130 = -t141 * t147 - t146 * t175;
t177 = t130 * t143;
t144 = sin(pkin(6));
t176 = t142 * t144;
t174 = t144 * t147;
t173 = t144 * t148;
t172 = t146 * t143;
t170 = t148 * t152;
t169 = t148 * t155;
t168 = t149 * qJ(2) + qJ(1);
t167 = t147 * pkin(1) + qJ(2) * t176;
t125 = pkin(3) * t143 * t152 + t179 * t148;
t126 = pkin(3) * t170 - t179 * t143;
t136 = pkin(3) * t155 + pkin(2);
t165 = t149 * t125 + t168 + (t126 * t146 + t136 * t141) * t144;
t131 = -t141 * t175 + t146 * t147;
t164 = t125 * t176 + t130 * t126 + t131 * t136 + t167;
t140 = sin(pkin(13));
t145 = cos(pkin(13));
t163 = t140 * t155 + t152 * t145;
t133 = -t152 * t140 + t145 * t155;
t122 = t163 * t143;
t124 = t163 * t148;
t103 = t122 * t176 + t124 * t130 + t131 * t133;
t162 = t103 * pkin(4) + t164;
t108 = t122 * t149 + (t124 * t146 + t133 * t141) * t144;
t161 = t108 * pkin(4) + t165;
t160 = t133 * t143;
t159 = t144 * t160;
t129 = t141 * t171 + t142 * t146;
t138 = t142 * pkin(1);
t158 = t128 * t126 + t129 * t136 + t138 + (-qJ(2) - t125) * t174;
t101 = -t122 * t174 + t124 * t128 + t129 * t133;
t157 = t101 * pkin(4) + t158;
t154 = cos(qJ(5));
t153 = cos(qJ(6));
t151 = sin(qJ(5));
t150 = sin(qJ(6));
t127 = -t144 * t172 + t148 * t149;
t123 = t133 * t148;
t114 = t142 * t173 - t177;
t113 = -t147 * t173 - t178;
t107 = (t123 * t146 - t141 * t163) * t144 + t149 * t160;
t105 = t108 * t154 + t127 * t151;
t104 = t108 * t151 - t127 * t154;
t102 = t130 * t123 - t131 * t163 + t142 * t159;
t100 = t123 * t128 - t129 * t163 - t147 * t159;
t97 = t103 * t154 + t114 * t151;
t96 = t103 * t151 - t114 * t154;
t95 = t101 * t154 + t113 * t151;
t94 = t101 * t151 - t113 * t154;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t147 - rSges(2,2) * t142) + g(2) * (rSges(2,1) * t142 + rSges(2,2) * t147) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t131 + rSges(3,2) * t130 + t167) + g(2) * (rSges(3,1) * t129 + rSges(3,2) * t128 + t138) + g(3) * (rSges(3,3) * t149 + t168) + (rSges(3,3) * t182 + g(3) * (rSges(3,1) * t141 + rSges(3,2) * t146) + (-rSges(3,3) - qJ(2)) * t181) * t144) - m(4) * (g(1) * (t131 * pkin(2) - pkin(8) * t177 + (t130 * t170 + t131 * t155) * rSges(4,1) + (t130 * t169 - t131 * t152) * rSges(4,2) + t114 * rSges(4,3) + t167) + g(2) * (t129 * pkin(2) - pkin(8) * t178 + t138 + (t128 * t170 + t129 * t155) * rSges(4,1) + (t128 * t169 - t129 * t152) * rSges(4,2) + t113 * rSges(4,3)) + (t156 * t182 + (-qJ(2) - t156) * t181) * t144 + (t127 * rSges(4,3) + t168 + t156 * t149 + (t141 * pkin(2) - pkin(8) * t172 + (t141 * t155 + t146 * t170) * rSges(4,1) + (-t141 * t152 + t146 * t169) * rSges(4,2)) * t144) * g(3)) - m(5) * (g(1) * (rSges(5,1) * t103 + rSges(5,2) * t102 + rSges(5,3) * t114 + t164) + g(2) * (rSges(5,1) * t101 + rSges(5,2) * t100 + rSges(5,3) * t113 + t158) + g(3) * (rSges(5,1) * t108 + rSges(5,2) * t107 + rSges(5,3) * t127 + t165)) - m(6) * (g(1) * (rSges(6,1) * t97 - rSges(6,2) * t96 + t184 * t102 + t162) + g(2) * (rSges(6,1) * t95 - rSges(6,2) * t94 + t184 * t100 + t157) + g(3) * (rSges(6,1) * t105 - rSges(6,2) * t104 + t184 * t107 + t161)) - m(7) * (g(1) * (t97 * pkin(5) - t102 * pkin(9) + (-t102 * t150 + t153 * t97) * rSges(7,1) + (-t102 * t153 - t150 * t97) * rSges(7,2) + t183 * t96 + t162) + g(2) * (t95 * pkin(5) - t100 * pkin(9) + (-t100 * t150 + t95 * t153) * rSges(7,1) + (-t100 * t153 - t150 * t95) * rSges(7,2) + t183 * t94 + t157) + g(3) * (t105 * pkin(5) - t107 * pkin(9) + (t105 * t153 - t107 * t150) * rSges(7,1) + (-t105 * t150 - t107 * t153) * rSges(7,2) + t183 * t104 + t161));
U  = t1;
