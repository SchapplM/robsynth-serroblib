% Calculate potential energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPPRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPPRRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:03
% EndTime: 2019-03-08 18:39:03
% DurationCPUTime: 0.61s
% Computational Cost: add. (951->139), mult. (2565->202), div. (0->0), fcn. (3324->18), ass. (0->73)
t173 = cos(pkin(13));
t176 = cos(pkin(7));
t177 = cos(pkin(6));
t170 = sin(pkin(7));
t171 = sin(pkin(6));
t206 = t170 * t171;
t154 = -t173 * t206 + t176 * t177;
t167 = sin(pkin(13));
t174 = cos(pkin(12));
t168 = sin(pkin(12));
t207 = t168 * t177;
t157 = -t167 * t174 - t173 * t207;
t203 = t171 * t176;
t149 = -t157 * t170 + t168 * t203;
t166 = sin(pkin(14));
t172 = cos(pkin(14));
t202 = t173 * t176;
t205 = t170 * t177;
t146 = t172 * t205 + (-t166 * t167 + t172 * t202) * t171;
t169 = sin(pkin(8));
t175 = cos(pkin(8));
t139 = -t146 * t169 + t154 * t175;
t158 = -t167 * t207 + t173 * t174;
t192 = t157 * t176 + t168 * t206;
t137 = -t158 * t166 + t192 * t172;
t129 = -t137 * t169 + t149 * t175;
t201 = t174 * t177;
t156 = t167 * t201 + t168 * t173;
t155 = -t167 * t168 + t173 * t201;
t204 = t171 * t174;
t193 = t155 * t176 - t170 * t204;
t135 = -t156 * t166 + t193 * t172;
t148 = -t155 * t170 - t174 * t203;
t128 = -t135 * t169 + t148 * t175;
t218 = rSges(6,3) + pkin(10);
t217 = pkin(11) + rSges(7,3);
t216 = cos(qJ(4));
t208 = t167 * t171;
t199 = t177 * qJ(2) + qJ(1);
t198 = t168 * t171 * qJ(2) + t174 * pkin(1);
t195 = t169 * t216;
t194 = t175 * t216;
t191 = t158 * pkin(2) + t149 * qJ(3) + t198;
t190 = pkin(2) * t208 + t154 * qJ(3) + t199;
t138 = t158 * t172 + t192 * t166;
t189 = t138 * pkin(3) + t129 * pkin(9) + t191;
t164 = t168 * pkin(1);
t188 = t156 * pkin(2) - qJ(2) * t204 + t148 * qJ(3) + t164;
t147 = t172 * t208 + (t171 * t202 + t205) * t166;
t187 = t147 * pkin(3) + t139 * pkin(9) + t190;
t180 = sin(qJ(4));
t122 = t138 * t216 + (t137 * t175 + t149 * t169) * t180;
t186 = t122 * pkin(4) + t189;
t127 = t147 * t216 + (t146 * t175 + t154 * t169) * t180;
t185 = t127 * pkin(4) + t187;
t136 = t156 * t172 + t193 * t166;
t184 = t136 * pkin(3) + t128 * pkin(9) + t188;
t120 = t136 * t216 + (t135 * t175 + t148 * t169) * t180;
t183 = t120 * pkin(4) + t184;
t182 = cos(qJ(5));
t181 = cos(qJ(6));
t179 = sin(qJ(5));
t178 = sin(qJ(6));
t126 = -t146 * t194 + t147 * t180 - t154 * t195;
t121 = -t137 * t194 + t138 * t180 - t149 * t195;
t119 = -t135 * t194 + t136 * t180 - t148 * t195;
t118 = t127 * t182 + t139 * t179;
t117 = t127 * t179 - t139 * t182;
t114 = t122 * t182 + t129 * t179;
t113 = t122 * t179 - t129 * t182;
t112 = t120 * t182 + t128 * t179;
t111 = t120 * t179 - t128 * t182;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t174 - rSges(2,2) * t168) + g(2) * (rSges(2,1) * t168 + rSges(2,2) * t174) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t158 + rSges(3,2) * t157 + t198) + g(2) * (rSges(3,1) * t156 + rSges(3,2) * t155 + t164) + g(3) * (rSges(3,3) * t177 + t199) + (g(1) * rSges(3,3) * t168 + g(3) * (rSges(3,1) * t167 + rSges(3,2) * t173) + g(2) * (-rSges(3,3) - qJ(2)) * t174) * t171) - m(4) * (g(1) * (rSges(4,1) * t138 + rSges(4,2) * t137 + rSges(4,3) * t149 + t191) + g(2) * (rSges(4,1) * t136 + rSges(4,2) * t135 + rSges(4,3) * t148 + t188) + g(3) * (rSges(4,1) * t147 + rSges(4,2) * t146 + rSges(4,3) * t154 + t190)) - m(5) * (g(1) * (rSges(5,1) * t122 - rSges(5,2) * t121 + rSges(5,3) * t129 + t189) + g(2) * (rSges(5,1) * t120 - rSges(5,2) * t119 + rSges(5,3) * t128 + t184) + g(3) * (rSges(5,1) * t127 - rSges(5,2) * t126 + rSges(5,3) * t139 + t187)) - m(6) * (g(1) * (rSges(6,1) * t114 - rSges(6,2) * t113 + t121 * t218 + t186) + g(2) * (rSges(6,1) * t112 - rSges(6,2) * t111 + t119 * t218 + t183) + g(3) * (rSges(6,1) * t118 - rSges(6,2) * t117 + t126 * t218 + t185)) - m(7) * (g(1) * (t114 * pkin(5) + t121 * pkin(10) + (t114 * t181 + t121 * t178) * rSges(7,1) + (-t114 * t178 + t121 * t181) * rSges(7,2) + t217 * t113 + t186) + g(2) * (t112 * pkin(5) + t119 * pkin(10) + (t112 * t181 + t119 * t178) * rSges(7,1) + (-t112 * t178 + t119 * t181) * rSges(7,2) + t217 * t111 + t183) + g(3) * (t118 * pkin(5) + t126 * pkin(10) + (t118 * t181 + t126 * t178) * rSges(7,1) + (-t118 * t178 + t126 * t181) * rSges(7,2) + t217 * t117 + t185));
U  = t1;
