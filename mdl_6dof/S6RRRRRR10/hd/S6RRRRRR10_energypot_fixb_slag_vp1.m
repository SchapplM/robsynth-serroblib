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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:07
% EndTime: 2018-11-23 10:29:08
% DurationCPUTime: 0.77s
% Computational Cost: add. (3309->159), mult. (3351->203), div. (0->0), fcn. (3324->30), ass. (0->93)
t207 = pkin(6) - qJ(2);
t196 = cos(t207) / 0.2e1;
t206 = pkin(6) + qJ(2);
t200 = cos(t206);
t188 = t196 + t200 / 0.2e1;
t218 = sin(qJ(2));
t219 = sin(qJ(1));
t225 = cos(qJ(1));
t177 = -t219 * t188 - t218 * t225;
t209 = sin(pkin(7));
t212 = cos(pkin(7));
t210 = sin(pkin(6));
t247 = t210 * t219;
t169 = -t177 * t209 + t212 * t247;
t194 = sin(t206) / 0.2e1;
t198 = sin(t207);
t183 = t194 + t198 / 0.2e1;
t213 = cos(pkin(6));
t174 = -t183 * t209 + t212 * t213;
t204 = pkin(7) + qJ(3);
t193 = sin(t204) / 0.2e1;
t205 = pkin(7) - qJ(3);
t197 = sin(t205);
t181 = t193 + t197 / 0.2e1;
t195 = cos(t205) / 0.2e1;
t199 = cos(t204);
t186 = t195 + t199 / 0.2e1;
t189 = t196 - t200 / 0.2e1;
t217 = sin(qJ(3));
t162 = t181 * t213 + t183 * t186 - t189 * t217;
t208 = sin(pkin(8));
t211 = cos(pkin(8));
t153 = -t162 * t208 + t174 * t211;
t184 = t194 - t198 / 0.2e1;
t224 = cos(qJ(2));
t178 = -t219 * t184 + t224 * t225;
t158 = t177 * t186 - t178 * t217 + t181 * t247;
t149 = -t158 * t208 + t169 * t211;
t175 = t188 * t225 - t219 * t218;
t176 = t184 * t225 + t219 * t224;
t246 = t210 * t225;
t156 = t175 * t186 - t176 * t217 - t181 * t246;
t168 = -t175 * t209 - t212 * t246;
t148 = -t156 * t208 + t168 * t211;
t258 = rSges(6,3) + pkin(13);
t257 = pkin(14) + rSges(7,3);
t256 = t213 * pkin(10) + pkin(9);
t244 = t225 * pkin(1) + pkin(10) * t247;
t243 = pkin(8) - qJ(4);
t242 = pkin(8) + qJ(4);
t240 = cos(t242);
t239 = sin(t243);
t238 = cos(t243) / 0.2e1;
t237 = sin(t242) / 0.2e1;
t236 = t189 * pkin(2) + t174 * pkin(11) + t256;
t235 = t178 * pkin(2) + t169 * pkin(11) + t244;
t234 = t238 + t240 / 0.2e1;
t233 = t237 + t239 / 0.2e1;
t182 = t193 - t197 / 0.2e1;
t187 = t195 - t199 / 0.2e1;
t223 = cos(qJ(3));
t163 = t182 * t183 + t187 * t213 + t189 * t223;
t232 = t163 * pkin(3) + t153 * pkin(12) + t236;
t159 = t177 * t182 + t178 * t223 + t187 * t247;
t231 = t159 * pkin(3) + t149 * pkin(12) + t235;
t202 = t219 * pkin(1);
t230 = t176 * pkin(2) - pkin(10) * t246 + t168 * pkin(11) + t202;
t180 = t237 - t239 / 0.2e1;
t185 = t238 - t240 / 0.2e1;
t222 = cos(qJ(4));
t145 = t162 * t180 + t163 * t222 + t174 * t185;
t229 = t145 * pkin(4) + t232;
t142 = t158 * t180 + t159 * t222 + t169 * t185;
t228 = t142 * pkin(4) + t231;
t157 = t175 * t182 + t176 * t223 - t187 * t246;
t227 = t157 * pkin(3) + t148 * pkin(12) + t230;
t140 = t156 * t180 + t157 * t222 + t168 * t185;
t226 = t140 * pkin(4) + t227;
t221 = cos(qJ(5));
t220 = cos(qJ(6));
t216 = sin(qJ(4));
t215 = sin(qJ(5));
t214 = sin(qJ(6));
t144 = -t162 * t234 + t163 * t216 - t174 * t233;
t141 = -t158 * t234 + t159 * t216 - t169 * t233;
t139 = -t156 * t234 + t157 * t216 - t168 * t233;
t136 = t145 * t221 + t153 * t215;
t135 = t145 * t215 - t153 * t221;
t134 = t142 * t221 + t149 * t215;
t133 = t142 * t215 - t149 * t221;
t132 = t140 * t221 + t148 * t215;
t131 = t140 * t215 - t148 * t221;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t225 - t219 * rSges(2,2)) + g(2) * (t219 * rSges(2,1) + rSges(2,2) * t225) + g(3) * (pkin(9) + rSges(2,3))) - m(3) * (g(1) * (t178 * rSges(3,1) + t177 * rSges(3,2) + rSges(3,3) * t247 + t244) + g(2) * (t176 * rSges(3,1) + t175 * rSges(3,2) + t202 + (-rSges(3,3) - pkin(10)) * t246) + g(3) * (rSges(3,1) * t189 + rSges(3,2) * t183 + rSges(3,3) * t213 + t256)) - m(4) * (g(1) * (t159 * rSges(4,1) + t158 * rSges(4,2) + t169 * rSges(4,3) + t235) + g(2) * (t157 * rSges(4,1) + t156 * rSges(4,2) + t168 * rSges(4,3) + t230) + g(3) * (t163 * rSges(4,1) + t162 * rSges(4,2) + t174 * rSges(4,3) + t236)) - m(5) * (g(1) * (t142 * rSges(5,1) - t141 * rSges(5,2) + t149 * rSges(5,3) + t231) + g(2) * (t140 * rSges(5,1) - t139 * rSges(5,2) + t148 * rSges(5,3) + t227) + g(3) * (t145 * rSges(5,1) - t144 * rSges(5,2) + t153 * rSges(5,3) + t232)) - m(6) * (g(1) * (t134 * rSges(6,1) - t133 * rSges(6,2) + t258 * t141 + t228) + g(2) * (t132 * rSges(6,1) - t131 * rSges(6,2) + t258 * t139 + t226) + g(3) * (t136 * rSges(6,1) - t135 * rSges(6,2) + t258 * t144 + t229)) - m(7) * (g(1) * (t134 * pkin(5) + t141 * pkin(13) + (t134 * t220 + t141 * t214) * rSges(7,1) + (-t134 * t214 + t141 * t220) * rSges(7,2) + t257 * t133 + t228) + g(2) * (t132 * pkin(5) + t139 * pkin(13) + (t132 * t220 + t139 * t214) * rSges(7,1) + (-t132 * t214 + t139 * t220) * rSges(7,2) + t257 * t131 + t226) + g(3) * (t136 * pkin(5) + t144 * pkin(13) + (t136 * t220 + t144 * t214) * rSges(7,1) + (-t136 * t214 + t144 * t220) * rSges(7,2) + t257 * t135 + t229));
U  = t1;
