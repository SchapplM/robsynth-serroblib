% Calculate potential energy for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:04
% EndTime: 2019-03-09 23:29:04
% DurationCPUTime: 0.54s
% Computational Cost: add. (526->135), mult. (1191->188), div. (0->0), fcn. (1486->16), ass. (0->64)
t147 = cos(pkin(6));
t153 = sin(qJ(1));
t156 = cos(qJ(2));
t171 = t153 * t156;
t152 = sin(qJ(2));
t157 = cos(qJ(1));
t173 = t152 * t157;
t128 = -t147 * t171 - t173;
t144 = sin(pkin(7));
t146 = cos(pkin(7));
t145 = sin(pkin(6));
t177 = t145 * t153;
t118 = -t128 * t144 + t146 * t177;
t176 = t145 * t156;
t125 = -t144 * t176 + t146 * t147;
t186 = pkin(11) + rSges(5,3);
t185 = pkin(12) + rSges(7,3);
t184 = cos(qJ(3));
t183 = t147 * pkin(9) + pkin(8);
t170 = t156 * t157;
t172 = t153 * t152;
t126 = t147 * t170 - t172;
t175 = t145 * t157;
t117 = -t126 * t144 - t146 * t175;
t150 = sin(qJ(4));
t182 = t117 * t150;
t181 = t118 * t150;
t180 = t125 * t150;
t178 = t145 * t152;
t169 = t157 * pkin(1) + pkin(9) * t177;
t166 = t144 * t184;
t165 = t146 * t184;
t164 = t145 * t166;
t129 = -t147 * t172 + t170;
t163 = t129 * pkin(2) + t118 * pkin(10) + t169;
t162 = pkin(2) * t178 + t125 * pkin(10) + t183;
t151 = sin(qJ(3));
t107 = -t128 * t165 + t129 * t151 - t153 * t164;
t108 = t129 * t184 + (t128 * t146 + t144 * t177) * t151;
t155 = cos(qJ(4));
t137 = pkin(4) * t155 + pkin(3);
t148 = -qJ(5) - pkin(11);
t161 = pkin(4) * t181 - t107 * t148 + t108 * t137 + t163;
t115 = -t147 * t166 + t151 * t178 - t165 * t176;
t116 = t147 * t144 * t151 + (t146 * t151 * t156 + t184 * t152) * t145;
t160 = pkin(4) * t180 - t115 * t148 + t116 * t137 + t162;
t127 = t147 * t173 + t171;
t141 = t153 * pkin(1);
t159 = t127 * pkin(2) - pkin(9) * t175 + t117 * pkin(10) + t141;
t105 = -t126 * t165 + t127 * t151 + t157 * t164;
t106 = t127 * t184 + (t126 * t146 - t144 * t175) * t151;
t158 = pkin(4) * t182 - t105 * t148 + t106 * t137 + t159;
t154 = cos(qJ(6));
t149 = sin(qJ(6));
t143 = qJ(4) + pkin(13);
t139 = cos(t143);
t138 = sin(t143);
t102 = t116 * t139 + t125 * t138;
t101 = t116 * t138 - t125 * t139;
t98 = t108 * t139 + t118 * t138;
t97 = t108 * t138 - t118 * t139;
t96 = t106 * t139 + t117 * t138;
t95 = t106 * t138 - t117 * t139;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t157 - t153 * rSges(2,2)) + g(2) * (t153 * rSges(2,1) + rSges(2,2) * t157) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t129 + rSges(3,2) * t128 + t169) + g(2) * (t127 * rSges(3,1) + t126 * rSges(3,2) + t141) + g(3) * (rSges(3,3) * t147 + t183) + (g(1) * rSges(3,3) * t153 + g(3) * (rSges(3,1) * t152 + rSges(3,2) * t156) + g(2) * (-rSges(3,3) - pkin(9)) * t157) * t145) - m(4) * (g(1) * (rSges(4,1) * t108 - rSges(4,2) * t107 + rSges(4,3) * t118 + t163) + g(2) * (t106 * rSges(4,1) - t105 * rSges(4,2) + t117 * rSges(4,3) + t159) + g(3) * (rSges(4,1) * t116 - rSges(4,2) * t115 + rSges(4,3) * t125 + t162)) - m(5) * (g(1) * (t108 * pkin(3) + (t108 * t155 + t181) * rSges(5,1) + (-t108 * t150 + t118 * t155) * rSges(5,2) + t186 * t107 + t163) + g(2) * (t106 * pkin(3) + (t106 * t155 + t182) * rSges(5,1) + (-t106 * t150 + t117 * t155) * rSges(5,2) + t186 * t105 + t159) + g(3) * (t116 * pkin(3) + (t116 * t155 + t180) * rSges(5,1) + (-t116 * t150 + t125 * t155) * rSges(5,2) + t186 * t115 + t162)) - m(6) * (g(1) * (rSges(6,1) * t98 - rSges(6,2) * t97 + rSges(6,3) * t107 + t161) + g(2) * (t96 * rSges(6,1) - t95 * rSges(6,2) + t105 * rSges(6,3) + t158) + g(3) * (rSges(6,1) * t102 - rSges(6,2) * t101 + rSges(6,3) * t115 + t160)) - m(7) * (g(1) * (t98 * pkin(5) + (t107 * t149 + t154 * t98) * rSges(7,1) + (t107 * t154 - t149 * t98) * rSges(7,2) + t185 * t97 + t161) + g(2) * (t96 * pkin(5) + (t105 * t149 + t154 * t96) * rSges(7,1) + (t105 * t154 - t149 * t96) * rSges(7,2) + t185 * t95 + t158) + g(3) * (t102 * pkin(5) + (t102 * t154 + t115 * t149) * rSges(7,1) + (-t102 * t149 + t115 * t154) * rSges(7,2) + t185 * t101 + t160));
U  = t1;
