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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
% StartTime: 2018-12-10 18:09:46
% EndTime: 2018-12-10 18:09:46
% DurationCPUTime: 0.73s
% Computational Cost: add. (3309->159), mult. (3351->203), div. (0->0), fcn. (3324->30), ass. (0->93)
t206 = pkin(6) - qJ(2);
t195 = cos(t206) / 0.2e1;
t205 = pkin(6) + qJ(2);
t199 = cos(t205);
t187 = t195 + t199 / 0.2e1;
t218 = sin(qJ(2));
t219 = sin(qJ(1));
t224 = cos(qJ(1));
t176 = -t219 * t187 - t224 * t218;
t209 = sin(pkin(7));
t213 = cos(pkin(7));
t210 = sin(pkin(6));
t246 = t210 * t219;
t168 = -t176 * t209 + t213 * t246;
t194 = sin(t205) / 0.2e1;
t198 = sin(t206);
t184 = t194 + t198 / 0.2e1;
t214 = cos(pkin(6));
t173 = -t184 * t209 + t214 * t213;
t203 = pkin(7) + pkin(14);
t192 = sin(t203) / 0.2e1;
t204 = pkin(7) - pkin(14);
t196 = sin(t204);
t179 = t192 + t196 / 0.2e1;
t193 = cos(t204) / 0.2e1;
t197 = cos(t203);
t181 = t193 + t197 / 0.2e1;
t188 = t195 - t199 / 0.2e1;
t207 = sin(pkin(14));
t161 = t214 * t179 + t184 * t181 - t188 * t207;
t208 = sin(pkin(8));
t212 = cos(pkin(8));
t152 = -t161 * t208 + t173 * t212;
t185 = t194 - t198 / 0.2e1;
t223 = cos(qJ(2));
t177 = -t219 * t185 + t224 * t223;
t157 = t176 * t181 - t177 * t207 + t179 * t246;
t148 = -t157 * t208 + t168 * t212;
t174 = t224 * t187 - t219 * t218;
t175 = t224 * t185 + t219 * t223;
t245 = t210 * t224;
t155 = t174 * t181 - t175 * t207 - t179 * t245;
t167 = -t174 * t209 - t213 * t245;
t147 = -t155 * t208 + t167 * t212;
t257 = rSges(6,3) + pkin(12);
t256 = pkin(13) + rSges(7,3);
t255 = t214 * pkin(10) + pkin(9);
t243 = t224 * pkin(1) + pkin(10) * t246;
t242 = pkin(8) - qJ(4);
t241 = pkin(8) + qJ(4);
t239 = cos(t241);
t238 = sin(t242);
t237 = cos(t242) / 0.2e1;
t236 = sin(t241) / 0.2e1;
t235 = t188 * pkin(2) + t173 * qJ(3) + t255;
t234 = t177 * pkin(2) + t168 * qJ(3) + t243;
t233 = t237 + t239 / 0.2e1;
t232 = t236 + t238 / 0.2e1;
t180 = t192 - t196 / 0.2e1;
t182 = t193 - t197 / 0.2e1;
t211 = cos(pkin(14));
t162 = t184 * t180 + t214 * t182 + t188 * t211;
t231 = t162 * pkin(3) + t152 * pkin(11) + t235;
t158 = t176 * t180 + t177 * t211 + t182 * t246;
t230 = t158 * pkin(3) + t148 * pkin(11) + t234;
t201 = t219 * pkin(1);
t229 = t175 * pkin(2) - pkin(10) * t245 + qJ(3) * t167 + t201;
t183 = t236 - t238 / 0.2e1;
t186 = t237 - t239 / 0.2e1;
t222 = cos(qJ(4));
t144 = t161 * t183 + t162 * t222 + t173 * t186;
t228 = t144 * pkin(4) + t231;
t141 = t157 * t183 + t158 * t222 + t168 * t186;
t227 = t141 * pkin(4) + t230;
t156 = t174 * t180 + t175 * t211 - t182 * t245;
t226 = t156 * pkin(3) + t147 * pkin(11) + t229;
t139 = t155 * t183 + t156 * t222 + t167 * t186;
t225 = t139 * pkin(4) + t226;
t221 = cos(qJ(5));
t220 = cos(qJ(6));
t217 = sin(qJ(4));
t216 = sin(qJ(5));
t215 = sin(qJ(6));
t143 = -t161 * t233 + t162 * t217 - t173 * t232;
t140 = -t157 * t233 + t158 * t217 - t168 * t232;
t138 = -t155 * t233 + t156 * t217 - t167 * t232;
t135 = t144 * t221 + t152 * t216;
t134 = t144 * t216 - t152 * t221;
t133 = t141 * t221 + t148 * t216;
t132 = t141 * t216 - t148 * t221;
t131 = t139 * t221 + t147 * t216;
t130 = t139 * t216 - t147 * t221;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t224 * rSges(2,1) - t219 * rSges(2,2)) + g(2) * (t219 * rSges(2,1) + t224 * rSges(2,2)) + g(3) * (pkin(9) + rSges(2,3))) - m(3) * (g(1) * (t177 * rSges(3,1) + t176 * rSges(3,2) + rSges(3,3) * t246 + t243) + g(2) * (t175 * rSges(3,1) + t174 * rSges(3,2) + t201 + (-rSges(3,3) - pkin(10)) * t245) + g(3) * (t188 * rSges(3,1) + t184 * rSges(3,2) + t214 * rSges(3,3) + t255)) - m(4) * (g(1) * (t158 * rSges(4,1) + t157 * rSges(4,2) + t168 * rSges(4,3) + t234) + g(2) * (t156 * rSges(4,1) + t155 * rSges(4,2) + t167 * rSges(4,3) + t229) + g(3) * (t162 * rSges(4,1) + t161 * rSges(4,2) + t173 * rSges(4,3) + t235)) - m(5) * (g(1) * (t141 * rSges(5,1) - t140 * rSges(5,2) + t148 * rSges(5,3) + t230) + g(2) * (t139 * rSges(5,1) - t138 * rSges(5,2) + t147 * rSges(5,3) + t226) + g(3) * (t144 * rSges(5,1) - t143 * rSges(5,2) + t152 * rSges(5,3) + t231)) - m(6) * (g(1) * (t133 * rSges(6,1) - t132 * rSges(6,2) + t140 * t257 + t227) + g(2) * (t131 * rSges(6,1) - t130 * rSges(6,2) + t138 * t257 + t225) + g(3) * (t135 * rSges(6,1) - t134 * rSges(6,2) + t257 * t143 + t228)) - m(7) * (g(1) * (t133 * pkin(5) + t140 * pkin(12) + (t133 * t220 + t140 * t215) * rSges(7,1) + (-t133 * t215 + t140 * t220) * rSges(7,2) + t256 * t132 + t227) + g(2) * (t131 * pkin(5) + t138 * pkin(12) + (t131 * t220 + t138 * t215) * rSges(7,1) + (-t131 * t215 + t138 * t220) * rSges(7,2) + t256 * t130 + t225) + g(3) * (t135 * pkin(5) + t143 * pkin(12) + (t135 * t220 + t143 * t215) * rSges(7,1) + (-t135 * t215 + t143 * t220) * rSges(7,2) + t256 * t134 + t228));
U  = t1;
