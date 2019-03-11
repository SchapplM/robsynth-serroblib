% Calculate potential energy for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:36
% EndTime: 2019-03-10 05:49:36
% DurationCPUTime: 0.61s
% Computational Cost: add. (951->139), mult. (2565->200), div. (0->0), fcn. (3324->18), ass. (0->73)
t171 = cos(pkin(6));
t177 = sin(qJ(1));
t181 = cos(qJ(2));
t200 = t177 * t181;
t176 = sin(qJ(2));
t182 = cos(qJ(1));
t202 = t176 * t182;
t157 = -t171 * t200 - t202;
t167 = sin(pkin(7));
t170 = cos(pkin(7));
t168 = sin(pkin(6));
t206 = t168 * t177;
t149 = -t157 * t167 + t170 * t206;
t154 = -t167 * t168 * t181 + t170 * t171;
t175 = sin(qJ(3));
t180 = cos(qJ(3));
t203 = t170 * t181;
t207 = t167 * t171;
t146 = t180 * t207 + (-t175 * t176 + t180 * t203) * t168;
t166 = sin(pkin(8));
t169 = cos(pkin(8));
t133 = -t146 * t166 + t154 * t169;
t199 = t181 * t182;
t201 = t177 * t176;
t158 = -t171 * t201 + t199;
t192 = t157 * t170 + t167 * t206;
t138 = -t158 * t175 + t180 * t192;
t129 = -t138 * t166 + t149 * t169;
t156 = t171 * t202 + t200;
t155 = t171 * t199 - t201;
t205 = t168 * t182;
t193 = t155 * t170 - t167 * t205;
t136 = -t156 * t175 + t180 * t193;
t148 = -t155 * t167 - t170 * t205;
t128 = -t136 * t166 + t148 * t169;
t218 = rSges(6,3) + pkin(13);
t217 = pkin(14) + rSges(7,3);
t216 = cos(qJ(4));
t215 = t171 * pkin(10) + pkin(9);
t198 = t182 * pkin(1) + pkin(10) * t206;
t195 = t166 * t216;
t194 = t169 * t216;
t191 = t158 * pkin(2) + t149 * pkin(11) + t198;
t190 = t168 * t176 * pkin(2) + t154 * pkin(11) + t215;
t139 = t158 * t180 + t175 * t192;
t189 = t139 * pkin(3) + t129 * pkin(12) + t191;
t164 = t177 * pkin(1);
t188 = t156 * pkin(2) - pkin(10) * t205 + t148 * pkin(11) + t164;
t147 = t175 * t207 + (t175 * t203 + t176 * t180) * t168;
t187 = t147 * pkin(3) + t133 * pkin(12) + t190;
t174 = sin(qJ(4));
t122 = t139 * t216 + (t138 * t169 + t149 * t166) * t174;
t186 = t122 * pkin(4) + t189;
t125 = t147 * t216 + (t146 * t169 + t154 * t166) * t174;
t185 = t125 * pkin(4) + t187;
t137 = t156 * t180 + t175 * t193;
t184 = t137 * pkin(3) + t128 * pkin(12) + t188;
t120 = t137 * t216 + (t136 * t169 + t148 * t166) * t174;
t183 = t120 * pkin(4) + t184;
t179 = cos(qJ(5));
t178 = cos(qJ(6));
t173 = sin(qJ(5));
t172 = sin(qJ(6));
t124 = -t146 * t194 + t147 * t174 - t154 * t195;
t121 = -t138 * t194 + t139 * t174 - t149 * t195;
t119 = -t136 * t194 + t137 * t174 - t148 * t195;
t116 = t125 * t179 + t133 * t173;
t115 = t125 * t173 - t133 * t179;
t114 = t122 * t179 + t129 * t173;
t113 = t122 * t173 - t129 * t179;
t112 = t120 * t179 + t128 * t173;
t111 = t120 * t173 - t128 * t179;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t182 - t177 * rSges(2,2)) + g(2) * (t177 * rSges(2,1) + rSges(2,2) * t182) + g(3) * (pkin(9) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t158 + rSges(3,2) * t157 + t198) + g(2) * (t156 * rSges(3,1) + t155 * rSges(3,2) + t164) + g(3) * (rSges(3,3) * t171 + t215) + (g(1) * rSges(3,3) * t177 + g(3) * (rSges(3,1) * t176 + rSges(3,2) * t181) + g(2) * (-rSges(3,3) - pkin(10)) * t182) * t168) - m(4) * (g(1) * (rSges(4,1) * t139 + rSges(4,2) * t138 + rSges(4,3) * t149 + t191) + g(2) * (t137 * rSges(4,1) + t136 * rSges(4,2) + t148 * rSges(4,3) + t188) + g(3) * (rSges(4,1) * t147 + rSges(4,2) * t146 + rSges(4,3) * t154 + t190)) - m(5) * (g(1) * (rSges(5,1) * t122 - rSges(5,2) * t121 + rSges(5,3) * t129 + t189) + g(2) * (t120 * rSges(5,1) - t119 * rSges(5,2) + t128 * rSges(5,3) + t184) + g(3) * (rSges(5,1) * t125 - rSges(5,2) * t124 + rSges(5,3) * t133 + t187)) - m(6) * (g(1) * (rSges(6,1) * t114 - rSges(6,2) * t113 + t121 * t218 + t186) + g(2) * (t112 * rSges(6,1) - t111 * rSges(6,2) + t119 * t218 + t183) + g(3) * (rSges(6,1) * t116 - rSges(6,2) * t115 + t124 * t218 + t185)) - m(7) * (g(1) * (t114 * pkin(5) + t121 * pkin(13) + (t114 * t178 + t121 * t172) * rSges(7,1) + (-t114 * t172 + t121 * t178) * rSges(7,2) + t217 * t113 + t186) + g(2) * (t112 * pkin(5) + t119 * pkin(13) + (t112 * t178 + t119 * t172) * rSges(7,1) + (-t112 * t172 + t119 * t178) * rSges(7,2) + t217 * t111 + t183) + g(3) * (t116 * pkin(5) + t124 * pkin(13) + (t116 * t178 + t124 * t172) * rSges(7,1) + (-t116 * t172 + t124 * t178) * rSges(7,2) + t217 * t115 + t185));
U  = t1;
