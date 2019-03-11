% Calculate potential energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:18
% EndTime: 2019-03-09 07:50:19
% DurationCPUTime: 0.61s
% Computational Cost: add. (951->139), mult. (2565->201), div. (0->0), fcn. (3324->18), ass. (0->72)
t166 = sin(pkin(14));
t173 = cos(pkin(6));
t182 = cos(qJ(1));
t170 = cos(pkin(14));
t178 = sin(qJ(1));
t199 = t178 * t170;
t157 = -t166 * t182 - t173 * t199;
t168 = sin(pkin(7));
t172 = cos(pkin(7));
t169 = sin(pkin(6));
t205 = t169 * t178;
t149 = -t157 * t168 + t172 * t205;
t154 = -t168 * t169 * t170 + t172 * t173;
t177 = sin(qJ(3));
t181 = cos(qJ(3));
t203 = t170 * t172;
t206 = t168 * t173;
t146 = t181 * t206 + (-t166 * t177 + t181 * t203) * t169;
t167 = sin(pkin(8));
t171 = cos(pkin(8));
t133 = -t146 * t167 + t154 * t171;
t200 = t178 * t166;
t158 = t170 * t182 - t173 * t200;
t192 = t157 * t172 + t168 * t205;
t138 = -t158 * t177 + t192 * t181;
t129 = -t138 * t167 + t149 * t171;
t201 = t173 * t182;
t156 = t166 * t201 + t199;
t155 = t170 * t201 - t200;
t204 = t169 * t182;
t193 = t155 * t172 - t168 * t204;
t136 = -t156 * t177 + t193 * t181;
t148 = -t155 * t168 - t172 * t204;
t128 = -t136 * t167 + t148 * t171;
t217 = rSges(6,3) + pkin(12);
t216 = pkin(13) + rSges(7,3);
t215 = cos(qJ(4));
t214 = t173 * qJ(2) + pkin(9);
t198 = t182 * pkin(1) + qJ(2) * t205;
t195 = t167 * t215;
t194 = t171 * t215;
t191 = t158 * pkin(2) + t149 * pkin(10) + t198;
t190 = t169 * t166 * pkin(2) + t154 * pkin(10) + t214;
t139 = t158 * t181 + t192 * t177;
t189 = t139 * pkin(3) + t129 * pkin(11) + t191;
t147 = t177 * t206 + (t166 * t181 + t177 * t203) * t169;
t188 = t147 * pkin(3) + t133 * pkin(11) + t190;
t164 = t178 * pkin(1);
t187 = t156 * pkin(2) + t148 * pkin(10) - qJ(2) * t204 + t164;
t176 = sin(qJ(4));
t122 = t139 * t215 + (t138 * t171 + t149 * t167) * t176;
t186 = t122 * pkin(4) + t189;
t125 = t147 * t215 + (t146 * t171 + t154 * t167) * t176;
t185 = t125 * pkin(4) + t188;
t137 = t156 * t181 + t193 * t177;
t184 = t137 * pkin(3) + t128 * pkin(11) + t187;
t120 = t137 * t215 + (t136 * t171 + t148 * t167) * t176;
t183 = t120 * pkin(4) + t184;
t180 = cos(qJ(5));
t179 = cos(qJ(6));
t175 = sin(qJ(5));
t174 = sin(qJ(6));
t124 = -t146 * t194 + t147 * t176 - t154 * t195;
t121 = -t138 * t194 + t139 * t176 - t149 * t195;
t119 = -t136 * t194 + t137 * t176 - t148 * t195;
t116 = t125 * t180 + t133 * t175;
t115 = t125 * t175 - t133 * t180;
t114 = t122 * t180 + t129 * t175;
t113 = t122 * t175 - t129 * t180;
t112 = t120 * t180 + t128 * t175;
t111 = t120 * t175 - t128 * t180;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t182 - t178 * rSges(2,2)) + g(2) * (t178 * rSges(2,1) + rSges(2,2) * t182) + g(3) * (pkin(9) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t158 + rSges(3,2) * t157 + t198) + g(2) * (t156 * rSges(3,1) + t155 * rSges(3,2) + t164) + g(3) * (rSges(3,3) * t173 + t214) + (g(1) * rSges(3,3) * t178 + g(3) * (rSges(3,1) * t166 + rSges(3,2) * t170) + g(2) * (-rSges(3,3) - qJ(2)) * t182) * t169) - m(4) * (g(1) * (rSges(4,1) * t139 + rSges(4,2) * t138 + rSges(4,3) * t149 + t191) + g(2) * (t137 * rSges(4,1) + t136 * rSges(4,2) + t148 * rSges(4,3) + t187) + g(3) * (rSges(4,1) * t147 + rSges(4,2) * t146 + rSges(4,3) * t154 + t190)) - m(5) * (g(1) * (rSges(5,1) * t122 - rSges(5,2) * t121 + rSges(5,3) * t129 + t189) + g(2) * (t120 * rSges(5,1) - t119 * rSges(5,2) + t128 * rSges(5,3) + t184) + g(3) * (rSges(5,1) * t125 - rSges(5,2) * t124 + rSges(5,3) * t133 + t188)) - m(6) * (g(1) * (rSges(6,1) * t114 - rSges(6,2) * t113 + t121 * t217 + t186) + g(2) * (t112 * rSges(6,1) - t111 * rSges(6,2) + t119 * t217 + t183) + g(3) * (rSges(6,1) * t116 - rSges(6,2) * t115 + t124 * t217 + t185)) - m(7) * (g(1) * (t114 * pkin(5) + t121 * pkin(12) + (t114 * t179 + t121 * t174) * rSges(7,1) + (-t114 * t174 + t121 * t179) * rSges(7,2) + t216 * t113 + t186) + g(2) * (t112 * pkin(5) + t119 * pkin(12) + (t112 * t179 + t119 * t174) * rSges(7,1) + (-t112 * t174 + t119 * t179) * rSges(7,2) + t216 * t111 + t183) + g(3) * (t116 * pkin(5) + t124 * pkin(12) + (t116 * t179 + t124 * t174) * rSges(7,1) + (-t116 * t174 + t124 * t179) * rSges(7,2) + t216 * t115 + t185));
U  = t1;
