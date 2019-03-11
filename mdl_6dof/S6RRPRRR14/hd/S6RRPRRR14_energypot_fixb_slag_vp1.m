% Calculate potential energy for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:57:11
% EndTime: 2019-03-09 14:57:11
% DurationCPUTime: 0.61s
% Computational Cost: add. (951->139), mult. (2565->199), div. (0->0), fcn. (3324->18), ass. (0->74)
t174 = cos(pkin(6));
t179 = sin(qJ(1));
t182 = cos(qJ(2));
t201 = t179 * t182;
t178 = sin(qJ(2));
t183 = cos(qJ(1));
t203 = t178 * t183;
t158 = -t174 * t201 - t203;
t169 = sin(pkin(7));
t173 = cos(pkin(7));
t170 = sin(pkin(6));
t207 = t170 * t179;
t150 = -t158 * t169 + t173 * t207;
t155 = -t169 * t170 * t182 + t173 * t174;
t167 = sin(pkin(14));
t171 = cos(pkin(14));
t204 = t173 * t182;
t209 = t169 * t174;
t147 = t171 * t209 + (-t167 * t178 + t171 * t204) * t170;
t168 = sin(pkin(8));
t172 = cos(pkin(8));
t134 = -t147 * t168 + t155 * t172;
t200 = t182 * t183;
t202 = t179 * t178;
t159 = -t174 * t202 + t200;
t193 = t158 * t173 + t169 * t207;
t139 = -t159 * t167 + t193 * t171;
t130 = -t139 * t168 + t150 * t172;
t157 = t174 * t203 + t201;
t156 = t174 * t200 - t202;
t206 = t170 * t183;
t194 = t156 * t173 - t169 * t206;
t137 = -t157 * t167 + t194 * t171;
t149 = -t156 * t169 - t173 * t206;
t129 = -t137 * t168 + t149 * t172;
t220 = rSges(6,3) + pkin(12);
t219 = pkin(13) + rSges(7,3);
t218 = cos(qJ(4));
t217 = t174 * pkin(10) + pkin(9);
t208 = t170 * t178;
t199 = t183 * pkin(1) + pkin(10) * t207;
t196 = t168 * t218;
t195 = t172 * t218;
t192 = t159 * pkin(2) + t150 * qJ(3) + t199;
t191 = pkin(2) * t208 + t155 * qJ(3) + t217;
t140 = t159 * t171 + t193 * t167;
t190 = t140 * pkin(3) + t130 * pkin(11) + t192;
t165 = t179 * pkin(1);
t189 = t157 * pkin(2) - pkin(10) * t206 + t149 * qJ(3) + t165;
t148 = t171 * t208 + (t170 * t204 + t209) * t167;
t188 = t148 * pkin(3) + t134 * pkin(11) + t191;
t177 = sin(qJ(4));
t123 = t140 * t218 + (t139 * t172 + t150 * t168) * t177;
t187 = t123 * pkin(4) + t190;
t126 = t148 * t218 + (t147 * t172 + t155 * t168) * t177;
t186 = t126 * pkin(4) + t188;
t138 = t157 * t171 + t194 * t167;
t185 = t138 * pkin(3) + t129 * pkin(11) + t189;
t121 = t138 * t218 + (t137 * t172 + t149 * t168) * t177;
t184 = t121 * pkin(4) + t185;
t181 = cos(qJ(5));
t180 = cos(qJ(6));
t176 = sin(qJ(5));
t175 = sin(qJ(6));
t125 = -t147 * t195 + t148 * t177 - t155 * t196;
t122 = -t139 * t195 + t140 * t177 - t150 * t196;
t120 = -t137 * t195 + t138 * t177 - t149 * t196;
t117 = t126 * t181 + t134 * t176;
t116 = t126 * t176 - t134 * t181;
t115 = t123 * t181 + t130 * t176;
t114 = t123 * t176 - t130 * t181;
t113 = t121 * t181 + t129 * t176;
t112 = t121 * t176 - t129 * t181;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t183 - t179 * rSges(2,2)) + g(2) * (t179 * rSges(2,1) + rSges(2,2) * t183) + g(3) * (pkin(9) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t159 + rSges(3,2) * t158 + t199) + g(2) * (t157 * rSges(3,1) + t156 * rSges(3,2) + t165) + g(3) * (rSges(3,3) * t174 + t217) + (g(1) * rSges(3,3) * t179 + g(3) * (rSges(3,1) * t178 + rSges(3,2) * t182) + g(2) * (-rSges(3,3) - pkin(10)) * t183) * t170) - m(4) * (g(1) * (rSges(4,1) * t140 + rSges(4,2) * t139 + rSges(4,3) * t150 + t192) + g(2) * (t138 * rSges(4,1) + t137 * rSges(4,2) + t149 * rSges(4,3) + t189) + g(3) * (rSges(4,1) * t148 + rSges(4,2) * t147 + rSges(4,3) * t155 + t191)) - m(5) * (g(1) * (rSges(5,1) * t123 - rSges(5,2) * t122 + rSges(5,3) * t130 + t190) + g(2) * (t121 * rSges(5,1) - t120 * rSges(5,2) + t129 * rSges(5,3) + t185) + g(3) * (rSges(5,1) * t126 - rSges(5,2) * t125 + rSges(5,3) * t134 + t188)) - m(6) * (g(1) * (rSges(6,1) * t115 - rSges(6,2) * t114 + t220 * t122 + t187) + g(2) * (t113 * rSges(6,1) - t112 * rSges(6,2) + t220 * t120 + t184) + g(3) * (rSges(6,1) * t117 - rSges(6,2) * t116 + t220 * t125 + t186)) - m(7) * (g(1) * (t115 * pkin(5) + t122 * pkin(12) + (t115 * t180 + t122 * t175) * rSges(7,1) + (-t115 * t175 + t122 * t180) * rSges(7,2) + t219 * t114 + t187) + g(2) * (t113 * pkin(5) + t120 * pkin(12) + (t113 * t180 + t120 * t175) * rSges(7,1) + (-t113 * t175 + t120 * t180) * rSges(7,2) + t219 * t112 + t184) + g(3) * (t117 * pkin(5) + t125 * pkin(12) + (t117 * t180 + t125 * t175) * rSges(7,1) + (-t117 * t175 + t125 * t180) * rSges(7,2) + t219 * t116 + t186));
U  = t1;
