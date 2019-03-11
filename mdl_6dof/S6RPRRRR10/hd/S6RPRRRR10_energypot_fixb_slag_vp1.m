% Calculate potential energy for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:34
% EndTime: 2019-03-09 07:27:35
% DurationCPUTime: 0.53s
% Computational Cost: add. (526->135), mult. (1191->190), div. (0->0), fcn. (1486->16), ass. (0->62)
t143 = sin(pkin(13));
t146 = cos(pkin(13));
t155 = cos(qJ(1));
t148 = cos(pkin(6));
t152 = sin(qJ(1));
t170 = t148 * t152;
t127 = -t143 * t155 - t146 * t170;
t144 = sin(pkin(7));
t147 = cos(pkin(7));
t145 = sin(pkin(6));
t173 = t145 * t152;
t117 = -t127 * t144 + t147 * t173;
t174 = t145 * t146;
t124 = -t144 * t174 + t147 * t148;
t183 = pkin(10) + rSges(5,3);
t182 = pkin(12) + rSges(7,3);
t181 = cos(qJ(3));
t180 = t148 * qJ(2) + pkin(8);
t169 = t148 * t155;
t125 = -t143 * t152 + t146 * t169;
t172 = t145 * t155;
t116 = -t125 * t144 - t147 * t172;
t150 = sin(qJ(4));
t179 = t116 * t150;
t178 = t117 * t150;
t177 = t124 * t150;
t175 = t143 * t145;
t168 = t155 * pkin(1) + qJ(2) * t173;
t165 = t144 * t181;
t164 = t147 * t181;
t163 = t145 * t165;
t128 = -t143 * t170 + t146 * t155;
t162 = t128 * pkin(2) + t117 * pkin(9) + t168;
t161 = pkin(2) * t175 + t124 * pkin(9) + t180;
t151 = sin(qJ(3));
t106 = -t127 * t164 + t128 * t151 - t152 * t163;
t107 = t128 * t181 + (t127 * t147 + t144 * t173) * t151;
t154 = cos(qJ(4));
t136 = pkin(4) * t154 + pkin(3);
t156 = -pkin(11) - pkin(10);
t160 = pkin(4) * t178 - t106 * t156 + t107 * t136 + t162;
t114 = -t148 * t165 + t151 * t175 - t164 * t174;
t115 = t148 * t144 * t151 + (t146 * t147 * t151 + t181 * t143) * t145;
t159 = pkin(4) * t177 - t114 * t156 + t115 * t136 + t161;
t126 = t143 * t169 + t146 * t152;
t140 = t152 * pkin(1);
t158 = t126 * pkin(2) + t116 * pkin(9) - qJ(2) * t172 + t140;
t104 = -t125 * t164 + t126 * t151 + t155 * t163;
t105 = t126 * t181 + (t125 * t147 - t144 * t172) * t151;
t157 = pkin(4) * t179 - t104 * t156 + t105 * t136 + t158;
t153 = cos(qJ(6));
t149 = sin(qJ(6));
t142 = qJ(4) + qJ(5);
t139 = cos(t142);
t138 = sin(t142);
t101 = t115 * t139 + t124 * t138;
t100 = t115 * t138 - t124 * t139;
t97 = t107 * t139 + t117 * t138;
t96 = t107 * t138 - t117 * t139;
t95 = t105 * t139 + t116 * t138;
t94 = t105 * t138 - t116 * t139;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t155 - rSges(2,2) * t152) + g(2) * (rSges(2,1) * t152 + rSges(2,2) * t155) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t128 + rSges(3,2) * t127 + t168) + g(2) * (rSges(3,1) * t126 + rSges(3,2) * t125 + t140) + g(3) * (rSges(3,3) * t148 + t180) + (g(1) * rSges(3,3) * t152 + g(3) * (rSges(3,1) * t143 + rSges(3,2) * t146) + g(2) * (-rSges(3,3) - qJ(2)) * t155) * t145) - m(4) * (g(1) * (rSges(4,1) * t107 - rSges(4,2) * t106 + rSges(4,3) * t117 + t162) + g(2) * (rSges(4,1) * t105 - rSges(4,2) * t104 + rSges(4,3) * t116 + t158) + g(3) * (rSges(4,1) * t115 - rSges(4,2) * t114 + rSges(4,3) * t124 + t161)) - m(5) * (g(1) * (t107 * pkin(3) + (t107 * t154 + t178) * rSges(5,1) + (-t107 * t150 + t117 * t154) * rSges(5,2) + t183 * t106 + t162) + g(2) * (t105 * pkin(3) + (t105 * t154 + t179) * rSges(5,1) + (-t105 * t150 + t116 * t154) * rSges(5,2) + t183 * t104 + t158) + g(3) * (t115 * pkin(3) + (t115 * t154 + t177) * rSges(5,1) + (-t115 * t150 + t124 * t154) * rSges(5,2) + t183 * t114 + t161)) - m(6) * (g(1) * (rSges(6,1) * t97 - rSges(6,2) * t96 + rSges(6,3) * t106 + t160) + g(2) * (rSges(6,1) * t95 - rSges(6,2) * t94 + rSges(6,3) * t104 + t157) + g(3) * (rSges(6,1) * t101 - rSges(6,2) * t100 + rSges(6,3) * t114 + t159)) - m(7) * (g(1) * (t97 * pkin(5) + (t106 * t149 + t153 * t97) * rSges(7,1) + (t106 * t153 - t149 * t97) * rSges(7,2) + t182 * t96 + t160) + g(2) * (t95 * pkin(5) + (t104 * t149 + t153 * t95) * rSges(7,1) + (t104 * t153 - t149 * t95) * rSges(7,2) + t182 * t94 + t157) + g(3) * (t101 * pkin(5) + (t101 * t153 + t114 * t149) * rSges(7,1) + (-t101 * t149 + t114 * t153) * rSges(7,2) + t182 * t100 + t159));
U  = t1;
