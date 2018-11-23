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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S6RRRRRR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 10:29:13
% EndTime: 2018-11-23 10:29:14
% DurationCPUTime: 1.20s
% Computational Cost: add. (3309->143), mult. (3360->170), div. (0->0), fcn. (3324->30), ass. (0->84)
t269 = m(6) + m(7);
t214 = pkin(6) - qJ(2);
t203 = cos(t214) / 0.2e1;
t213 = pkin(6) + qJ(2);
t207 = cos(t213);
t195 = t203 + t207 / 0.2e1;
t225 = sin(qJ(2));
t226 = sin(qJ(1));
t232 = cos(qJ(1));
t184 = -t226 * t195 - t225 * t232;
t216 = sin(pkin(7));
t219 = cos(pkin(7));
t217 = sin(pkin(6));
t258 = t217 * t226;
t176 = -t184 * t216 + t219 * t258;
t201 = sin(t213) / 0.2e1;
t205 = sin(t214);
t190 = t201 + t205 / 0.2e1;
t220 = cos(pkin(6));
t181 = -t190 * t216 + t219 * t220;
t211 = pkin(7) + qJ(3);
t200 = sin(t211) / 0.2e1;
t212 = pkin(7) - qJ(3);
t204 = sin(t212);
t188 = t200 + t204 / 0.2e1;
t202 = cos(t212) / 0.2e1;
t206 = cos(t211);
t193 = t202 + t206 / 0.2e1;
t196 = t203 - t207 / 0.2e1;
t224 = sin(qJ(3));
t169 = t188 * t220 + t190 * t193 - t196 * t224;
t215 = sin(pkin(8));
t218 = cos(pkin(8));
t160 = -t169 * t215 + t181 * t218;
t191 = t201 - t205 / 0.2e1;
t231 = cos(qJ(2));
t185 = -t226 * t191 + t231 * t232;
t165 = t184 * t193 - t185 * t224 + t188 * t258;
t156 = -t165 * t215 + t176 * t218;
t182 = t195 * t232 - t226 * t225;
t183 = t191 * t232 + t226 * t231;
t257 = t217 * t232;
t163 = t182 * t193 - t183 * t224 - t188 * t257;
t261 = t182 * t216;
t175 = -t219 * t257 - t261;
t155 = -t163 * t215 + t175 * t218;
t268 = t220 * pkin(10) + pkin(9);
t255 = t232 * pkin(1) + pkin(10) * t258;
t254 = pkin(8) - qJ(4);
t253 = pkin(8) + qJ(4);
t250 = cos(t253);
t249 = sin(t254);
t248 = -m(7) * pkin(14) + mrSges(6,2) - mrSges(7,3);
t247 = cos(t254) / 0.2e1;
t246 = sin(t253) / 0.2e1;
t209 = t226 * pkin(1);
t245 = t183 * pkin(2) - pkin(11) * t261 + t209;
t244 = t196 * pkin(2) + t181 * pkin(11) + t268;
t243 = t185 * pkin(2) + t176 * pkin(11) + t255;
t221 = sin(qJ(6));
t227 = cos(qJ(6));
t242 = -m(7) * pkin(5) - mrSges(7,1) * t227 + mrSges(7,2) * t221 - mrSges(6,1);
t241 = t247 + t250 / 0.2e1;
t240 = t246 + t249 / 0.2e1;
t239 = -mrSges(7,1) * t221 - mrSges(7,2) * t227 - pkin(13) * t269 + mrSges(5,2) - mrSges(6,3);
t189 = t200 - t204 / 0.2e1;
t194 = t202 - t206 / 0.2e1;
t230 = cos(qJ(3));
t164 = t182 * t189 + t183 * t230 - t194 * t257;
t238 = t164 * pkin(3) + t155 * pkin(12) + t245;
t170 = t189 * t190 + t194 * t220 + t196 * t230;
t236 = t170 * pkin(3) + t160 * pkin(12) + t244;
t166 = t184 * t189 + t185 * t230 + t194 * t258;
t235 = t166 * pkin(3) + t156 * pkin(12) + t243;
t229 = cos(qJ(4));
t228 = cos(qJ(5));
t223 = sin(qJ(4));
t222 = sin(qJ(5));
t192 = t247 - t250 / 0.2e1;
t187 = t246 - t249 / 0.2e1;
t152 = t169 * t187 + t170 * t229 + t181 * t192;
t149 = t165 * t187 + t166 * t229 + t176 * t192;
t147 = t163 * t187 + t164 * t229 + t175 * t192;
t1 = (-mrSges(1,3) - m(2) * pkin(9) - mrSges(2,3) - m(3) * t268 - t196 * mrSges(3,1) - t190 * mrSges(3,2) - t220 * mrSges(3,3) - m(4) * t244 - t170 * mrSges(4,1) - t169 * mrSges(4,2) - t181 * mrSges(4,3) - m(5) * t236 - t152 * mrSges(5,1) - t160 * mrSges(5,3) + t242 * (t152 * t228 + t160 * t222) + t239 * (-t169 * t241 + t170 * t223 - t181 * t240) + t248 * (t152 * t222 - t160 * t228) + t269 * (-t152 * pkin(4) - t236)) * g(3) + (-mrSges(1,2) - t226 * mrSges(2,1) - m(3) * t209 - t183 * mrSges(3,1) - t182 * mrSges(3,2) - m(4) * t245 - t164 * mrSges(4,1) - t163 * mrSges(4,2) - t175 * mrSges(4,3) - m(5) * t238 - t147 * mrSges(5,1) - t155 * mrSges(5,3) + t242 * (t147 * t228 + t155 * t222) + t239 * (-t163 * t241 + t164 * t223 - t175 * t240) + t248 * (t147 * t222 - t155 * t228) + (-mrSges(2,2) + (m(3) * pkin(10) + mrSges(3,3) + (m(4) + m(5) + t269) * (pkin(11) * t219 + pkin(10))) * t217) * t232 + t269 * (-t147 * pkin(4) - t238)) * g(2) + (-mrSges(1,1) - t232 * mrSges(2,1) - m(3) * t255 - t185 * mrSges(3,1) - t184 * mrSges(3,2) - m(4) * t243 - t166 * mrSges(4,1) - t165 * mrSges(4,2) - t176 * mrSges(4,3) - m(5) * t235 - t149 * mrSges(5,1) - t156 * mrSges(5,3) + (-t217 * mrSges(3,3) + mrSges(2,2)) * t226 + t242 * (t149 * t228 + t156 * t222) + t239 * (-t165 * t241 + t166 * t223 - t176 * t240) + t248 * (t149 * t222 - t156 * t228) + t269 * (-t149 * pkin(4) - t235)) * g(1);
U  = t1;
