% Calculate potential energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:11
% EndTime: 2019-03-08 19:07:12
% DurationCPUTime: 0.60s
% Computational Cost: add. (951->139), mult. (2565->203), div. (0->0), fcn. (3324->18), ass. (0->72)
t170 = cos(pkin(14));
t173 = cos(pkin(7));
t174 = cos(pkin(6));
t168 = sin(pkin(7));
t169 = sin(pkin(6));
t205 = t168 * t169;
t153 = -t170 * t205 + t173 * t174;
t165 = sin(pkin(14));
t171 = cos(pkin(13));
t166 = sin(pkin(13));
t206 = t166 * t174;
t156 = -t165 * t171 - t170 * t206;
t202 = t169 * t173;
t148 = -t156 * t168 + t166 * t202;
t178 = sin(qJ(3));
t181 = cos(qJ(3));
t201 = t170 * t173;
t204 = t168 * t174;
t145 = t181 * t204 + (-t165 * t178 + t181 * t201) * t169;
t167 = sin(pkin(8));
t172 = cos(pkin(8));
t138 = -t145 * t167 + t153 * t172;
t157 = -t165 * t206 + t170 * t171;
t191 = t156 * t173 + t166 * t205;
t136 = -t157 * t178 + t191 * t181;
t128 = -t136 * t167 + t148 * t172;
t200 = t171 * t174;
t155 = t165 * t200 + t166 * t170;
t154 = -t165 * t166 + t170 * t200;
t203 = t169 * t171;
t192 = t154 * t173 - t168 * t203;
t134 = -t155 * t178 + t192 * t181;
t147 = -t154 * t168 - t171 * t202;
t127 = -t134 * t167 + t147 * t172;
t216 = rSges(6,3) + pkin(11);
t215 = pkin(12) + rSges(7,3);
t214 = cos(qJ(4));
t198 = t174 * qJ(2) + qJ(1);
t197 = t166 * t169 * qJ(2) + t171 * pkin(1);
t194 = t167 * t214;
t193 = t172 * t214;
t190 = t157 * pkin(2) + t148 * pkin(9) + t197;
t189 = t169 * t165 * pkin(2) + t153 * pkin(9) + t198;
t137 = t157 * t181 + t191 * t178;
t188 = t137 * pkin(3) + t128 * pkin(10) + t190;
t163 = t166 * pkin(1);
t187 = t155 * pkin(2) + t147 * pkin(9) - qJ(2) * t203 + t163;
t146 = t178 * t204 + (t165 * t181 + t178 * t201) * t169;
t186 = t146 * pkin(3) + t138 * pkin(10) + t189;
t177 = sin(qJ(4));
t119 = t137 * t214 + (t136 * t172 + t148 * t167) * t177;
t185 = t119 * pkin(4) + t188;
t126 = t146 * t214 + (t145 * t172 + t153 * t167) * t177;
t184 = t126 * pkin(4) + t186;
t135 = t155 * t181 + t192 * t178;
t183 = t135 * pkin(3) + t127 * pkin(10) + t187;
t117 = t135 * t214 + (t134 * t172 + t147 * t167) * t177;
t182 = t117 * pkin(4) + t183;
t180 = cos(qJ(5));
t179 = cos(qJ(6));
t176 = sin(qJ(5));
t175 = sin(qJ(6));
t125 = -t145 * t193 + t146 * t177 - t153 * t194;
t121 = t126 * t180 + t138 * t176;
t120 = t126 * t176 - t138 * t180;
t118 = -t136 * t193 + t137 * t177 - t148 * t194;
t116 = -t134 * t193 + t135 * t177 - t147 * t194;
t113 = t119 * t180 + t128 * t176;
t112 = t119 * t176 - t128 * t180;
t111 = t117 * t180 + t127 * t176;
t110 = t117 * t176 - t127 * t180;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t171 - rSges(2,2) * t166) + g(2) * (rSges(2,1) * t166 + rSges(2,2) * t171) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t157 + rSges(3,2) * t156 + t197) + g(2) * (rSges(3,1) * t155 + rSges(3,2) * t154 + t163) + g(3) * (rSges(3,3) * t174 + t198) + (g(1) * rSges(3,3) * t166 + g(3) * (rSges(3,1) * t165 + rSges(3,2) * t170) + g(2) * (-rSges(3,3) - qJ(2)) * t171) * t169) - m(4) * (g(1) * (rSges(4,1) * t137 + rSges(4,2) * t136 + rSges(4,3) * t148 + t190) + g(2) * (rSges(4,1) * t135 + rSges(4,2) * t134 + rSges(4,3) * t147 + t187) + g(3) * (rSges(4,1) * t146 + rSges(4,2) * t145 + rSges(4,3) * t153 + t189)) - m(5) * (g(1) * (rSges(5,1) * t119 - rSges(5,2) * t118 + rSges(5,3) * t128 + t188) + g(2) * (rSges(5,1) * t117 - rSges(5,2) * t116 + rSges(5,3) * t127 + t183) + g(3) * (rSges(5,1) * t126 - rSges(5,2) * t125 + rSges(5,3) * t138 + t186)) - m(6) * (g(1) * (rSges(6,1) * t113 - rSges(6,2) * t112 + t118 * t216 + t185) + g(2) * (rSges(6,1) * t111 - rSges(6,2) * t110 + t116 * t216 + t182) + g(3) * (rSges(6,1) * t121 - rSges(6,2) * t120 + t125 * t216 + t184)) - m(7) * (g(1) * (t113 * pkin(5) + t118 * pkin(11) + (t113 * t179 + t118 * t175) * rSges(7,1) + (-t113 * t175 + t118 * t179) * rSges(7,2) + t215 * t112 + t185) + g(2) * (t111 * pkin(5) + t116 * pkin(11) + (t111 * t179 + t116 * t175) * rSges(7,1) + (-t111 * t175 + t116 * t179) * rSges(7,2) + t215 * t110 + t182) + g(3) * (t121 * pkin(5) + t125 * pkin(11) + (t121 * t179 + t125 * t175) * rSges(7,1) + (-t121 * t175 + t125 * t179) * rSges(7,2) + t215 * t120 + t184));
U  = t1;
