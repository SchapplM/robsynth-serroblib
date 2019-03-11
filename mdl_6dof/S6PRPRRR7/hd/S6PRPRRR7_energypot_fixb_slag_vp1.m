% Calculate potential energy for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:29
% EndTime: 2019-03-08 20:50:30
% DurationCPUTime: 0.61s
% Computational Cost: add. (951->139), mult. (2565->202), div. (0->0), fcn. (3324->18), ass. (0->73)
t174 = cos(pkin(7));
t175 = cos(pkin(6));
t182 = cos(qJ(2));
t169 = sin(pkin(7));
t170 = sin(pkin(6));
t208 = t169 * t170;
t154 = t175 * t174 - t182 * t208;
t167 = sin(pkin(13));
t172 = cos(pkin(13));
t179 = sin(qJ(2));
t200 = t175 * t182;
t157 = -t167 * t200 - t172 * t179;
t205 = t170 * t174;
t149 = -t157 * t169 + t167 * t205;
t166 = sin(pkin(14));
t171 = cos(pkin(14));
t203 = t174 * t182;
t207 = t169 * t175;
t146 = t171 * t207 + (-t166 * t179 + t171 * t203) * t170;
t168 = sin(pkin(8));
t173 = cos(pkin(8));
t135 = -t146 * t168 + t154 * t173;
t201 = t175 * t179;
t158 = -t167 * t201 + t172 * t182;
t192 = t157 * t174 + t167 * t208;
t138 = -t158 * t166 + t171 * t192;
t129 = -t138 * t168 + t149 * t173;
t156 = t167 * t182 + t172 * t201;
t155 = -t167 * t179 + t172 * t200;
t206 = t170 * t172;
t193 = t155 * t174 - t169 * t206;
t136 = -t156 * t166 + t171 * t193;
t148 = -t155 * t169 - t172 * t205;
t128 = -t136 * t168 + t148 * t173;
t218 = rSges(6,3) + pkin(11);
t217 = pkin(12) + rSges(7,3);
t216 = cos(qJ(4));
t204 = t170 * t179;
t199 = t175 * pkin(9) + qJ(1);
t198 = t167 * t170 * pkin(9) + t172 * pkin(1);
t195 = t168 * t216;
t194 = t173 * t216;
t191 = t158 * pkin(2) + t149 * qJ(3) + t198;
t190 = pkin(2) * t204 + t154 * qJ(3) + t199;
t139 = t158 * t171 + t166 * t192;
t189 = t139 * pkin(3) + t129 * pkin(10) + t191;
t163 = t167 * pkin(1);
t188 = t156 * pkin(2) - pkin(9) * t206 + t148 * qJ(3) + t163;
t147 = t171 * t204 + (t170 * t203 + t207) * t166;
t187 = t147 * pkin(3) + t135 * pkin(10) + t190;
t178 = sin(qJ(4));
t122 = t139 * t216 + (t138 * t173 + t149 * t168) * t178;
t186 = t122 * pkin(4) + t189;
t127 = t147 * t216 + (t146 * t173 + t154 * t168) * t178;
t185 = t127 * pkin(4) + t187;
t137 = t156 * t171 + t166 * t193;
t184 = t137 * pkin(3) + t128 * pkin(10) + t188;
t120 = t137 * t216 + (t136 * t173 + t148 * t168) * t178;
t183 = t120 * pkin(4) + t184;
t181 = cos(qJ(5));
t180 = cos(qJ(6));
t177 = sin(qJ(5));
t176 = sin(qJ(6));
t126 = -t146 * t194 + t147 * t178 - t154 * t195;
t121 = -t138 * t194 + t139 * t178 - t149 * t195;
t119 = -t136 * t194 + t137 * t178 - t148 * t195;
t116 = t127 * t181 + t135 * t177;
t115 = t127 * t177 - t135 * t181;
t114 = t122 * t181 + t129 * t177;
t113 = t122 * t177 - t129 * t181;
t112 = t120 * t181 + t128 * t177;
t111 = t120 * t177 - t128 * t181;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t172 - rSges(2,2) * t167) + g(2) * (rSges(2,1) * t167 + rSges(2,2) * t172) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t158 + rSges(3,2) * t157 + t198) + g(2) * (rSges(3,1) * t156 + rSges(3,2) * t155 + t163) + g(3) * (t175 * rSges(3,3) + t199) + (g(1) * rSges(3,3) * t167 + g(3) * (rSges(3,1) * t179 + rSges(3,2) * t182) + g(2) * (-rSges(3,3) - pkin(9)) * t172) * t170) - m(4) * (g(1) * (rSges(4,1) * t139 + rSges(4,2) * t138 + rSges(4,3) * t149 + t191) + g(2) * (rSges(4,1) * t137 + rSges(4,2) * t136 + rSges(4,3) * t148 + t188) + g(3) * (t147 * rSges(4,1) + t146 * rSges(4,2) + t154 * rSges(4,3) + t190)) - m(5) * (g(1) * (rSges(5,1) * t122 - rSges(5,2) * t121 + t129 * rSges(5,3) + t189) + g(2) * (rSges(5,1) * t120 - rSges(5,2) * t119 + rSges(5,3) * t128 + t184) + g(3) * (t127 * rSges(5,1) - t126 * rSges(5,2) + t135 * rSges(5,3) + t187)) - m(6) * (g(1) * (rSges(6,1) * t114 - rSges(6,2) * t113 + t121 * t218 + t186) + g(2) * (rSges(6,1) * t112 - rSges(6,2) * t111 + t119 * t218 + t183) + g(3) * (t116 * rSges(6,1) - t115 * rSges(6,2) + t126 * t218 + t185)) - m(7) * (g(1) * (t114 * pkin(5) + t121 * pkin(11) + (t114 * t180 + t121 * t176) * rSges(7,1) + (-t114 * t176 + t121 * t180) * rSges(7,2) + t217 * t113 + t186) + g(2) * (t112 * pkin(5) + t119 * pkin(11) + (t112 * t180 + t119 * t176) * rSges(7,1) + (-t112 * t176 + t119 * t180) * rSges(7,2) + t217 * t111 + t183) + g(3) * (t116 * pkin(5) + t126 * pkin(11) + (t116 * t180 + t126 * t176) * rSges(7,1) + (-t116 * t176 + t126 * t180) * rSges(7,2) + t217 * t115 + t185));
U  = t1;
