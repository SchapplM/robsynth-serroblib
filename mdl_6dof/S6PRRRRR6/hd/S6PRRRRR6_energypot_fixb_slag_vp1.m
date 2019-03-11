% Calculate potential energy for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:13:58
% EndTime: 2019-03-09 01:13:59
% DurationCPUTime: 0.61s
% Computational Cost: add. (951->139), mult. (2565->203), div. (0->0), fcn. (3324->18), ass. (0->72)
t171 = cos(pkin(7));
t172 = cos(pkin(6));
t181 = cos(qJ(2));
t167 = sin(pkin(7));
t168 = sin(pkin(6));
t206 = t167 * t168;
t153 = t172 * t171 - t181 * t206;
t165 = sin(pkin(14));
t169 = cos(pkin(14));
t177 = sin(qJ(2));
t199 = t172 * t181;
t156 = -t165 * t199 - t169 * t177;
t203 = t168 * t171;
t148 = -t156 * t167 + t165 * t203;
t176 = sin(qJ(3));
t180 = cos(qJ(3));
t202 = t171 * t181;
t205 = t167 * t172;
t145 = t180 * t205 + (-t176 * t177 + t180 * t202) * t168;
t166 = sin(pkin(8));
t170 = cos(pkin(8));
t134 = -t145 * t166 + t153 * t170;
t200 = t172 * t177;
t157 = -t165 * t200 + t169 * t181;
t191 = t156 * t171 + t165 * t206;
t137 = -t157 * t176 + t191 * t180;
t128 = -t137 * t166 + t148 * t170;
t155 = t165 * t181 + t169 * t200;
t154 = -t165 * t177 + t169 * t199;
t204 = t168 * t169;
t192 = t154 * t171 - t167 * t204;
t135 = -t155 * t176 + t192 * t180;
t147 = -t154 * t167 - t169 * t203;
t127 = -t135 * t166 + t147 * t170;
t216 = rSges(6,3) + pkin(12);
t215 = pkin(13) + rSges(7,3);
t214 = cos(qJ(4));
t198 = t172 * pkin(9) + qJ(1);
t197 = t165 * t168 * pkin(9) + t169 * pkin(1);
t194 = t166 * t214;
t193 = t170 * t214;
t190 = t157 * pkin(2) + t148 * pkin(10) + t197;
t189 = t168 * t177 * pkin(2) + t153 * pkin(10) + t198;
t138 = t157 * t180 + t191 * t176;
t188 = t138 * pkin(3) + t128 * pkin(11) + t190;
t162 = t165 * pkin(1);
t187 = t155 * pkin(2) - pkin(9) * t204 + t147 * pkin(10) + t162;
t146 = t176 * t205 + (t176 * t202 + t177 * t180) * t168;
t186 = t146 * pkin(3) + t134 * pkin(11) + t189;
t175 = sin(qJ(4));
t121 = t138 * t214 + (t137 * t170 + t148 * t166) * t175;
t185 = t121 * pkin(4) + t188;
t126 = t146 * t214 + (t145 * t170 + t153 * t166) * t175;
t184 = t126 * pkin(4) + t186;
t136 = t155 * t180 + t192 * t176;
t183 = t136 * pkin(3) + t127 * pkin(11) + t187;
t119 = t136 * t214 + (t135 * t170 + t147 * t166) * t175;
t182 = t119 * pkin(4) + t183;
t179 = cos(qJ(5));
t178 = cos(qJ(6));
t174 = sin(qJ(5));
t173 = sin(qJ(6));
t125 = -t145 * t193 + t146 * t175 - t153 * t194;
t120 = -t137 * t193 + t138 * t175 - t148 * t194;
t118 = -t135 * t193 + t136 * t175 - t147 * t194;
t117 = t126 * t179 + t134 * t174;
t116 = t126 * t174 - t134 * t179;
t113 = t121 * t179 + t128 * t174;
t112 = t121 * t174 - t128 * t179;
t111 = t119 * t179 + t127 * t174;
t110 = t119 * t174 - t127 * t179;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t169 - rSges(2,2) * t165) + g(2) * (rSges(2,1) * t165 + rSges(2,2) * t169) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t157 + rSges(3,2) * t156 + t197) + g(2) * (rSges(3,1) * t155 + rSges(3,2) * t154 + t162) + g(3) * (t172 * rSges(3,3) + t198) + (g(1) * rSges(3,3) * t165 + g(3) * (rSges(3,1) * t177 + rSges(3,2) * t181) + g(2) * (-rSges(3,3) - pkin(9)) * t169) * t168) - m(4) * (g(1) * (rSges(4,1) * t138 + rSges(4,2) * t137 + rSges(4,3) * t148 + t190) + g(2) * (rSges(4,1) * t136 + rSges(4,2) * t135 + rSges(4,3) * t147 + t187) + g(3) * (t146 * rSges(4,1) + t145 * rSges(4,2) + t153 * rSges(4,3) + t189)) - m(5) * (g(1) * (rSges(5,1) * t121 - rSges(5,2) * t120 + rSges(5,3) * t128 + t188) + g(2) * (rSges(5,1) * t119 - rSges(5,2) * t118 + rSges(5,3) * t127 + t183) + g(3) * (t126 * rSges(5,1) - t125 * rSges(5,2) + t134 * rSges(5,3) + t186)) - m(6) * (g(1) * (rSges(6,1) * t113 - rSges(6,2) * t112 + t216 * t120 + t185) + g(2) * (rSges(6,1) * t111 - rSges(6,2) * t110 + t216 * t118 + t182) + g(3) * (t117 * rSges(6,1) - t116 * rSges(6,2) + t216 * t125 + t184)) - m(7) * (g(1) * (t113 * pkin(5) + t120 * pkin(12) + (t113 * t178 + t120 * t173) * rSges(7,1) + (-t113 * t173 + t120 * t178) * rSges(7,2) + t215 * t112 + t185) + g(2) * (t111 * pkin(5) + t118 * pkin(12) + (t111 * t178 + t118 * t173) * rSges(7,1) + (-t111 * t173 + t118 * t178) * rSges(7,2) + t215 * t110 + t182) + g(3) * (t117 * pkin(5) + t125 * pkin(12) + (t117 * t178 + t125 * t173) * rSges(7,1) + (-t117 * t173 + t125 * t178) * rSges(7,2) + t215 * t116 + t184));
U  = t1;
