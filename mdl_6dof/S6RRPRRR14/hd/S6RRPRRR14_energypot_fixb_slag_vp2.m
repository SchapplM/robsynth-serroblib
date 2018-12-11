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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:09:46
% EndTime: 2018-12-10 18:09:47
% DurationCPUTime: 1.24s
% Computational Cost: add. (3309->143), mult. (3360->170), div. (0->0), fcn. (3324->30), ass. (0->84)
t267 = m(6) + m(7);
t212 = pkin(6) - qJ(2);
t201 = cos(t212) / 0.2e1;
t211 = pkin(6) + qJ(2);
t205 = cos(t211);
t193 = t201 + t205 / 0.2e1;
t224 = sin(qJ(2));
t225 = sin(qJ(1));
t230 = cos(qJ(1));
t182 = -t225 * t193 - t224 * t230;
t215 = sin(pkin(7));
t219 = cos(pkin(7));
t216 = sin(pkin(6));
t256 = t216 * t225;
t174 = -t182 * t215 + t219 * t256;
t200 = sin(t211) / 0.2e1;
t204 = sin(t212);
t190 = t200 + t204 / 0.2e1;
t220 = cos(pkin(6));
t179 = -t190 * t215 + t219 * t220;
t209 = pkin(7) + pkin(14);
t198 = sin(t209) / 0.2e1;
t210 = pkin(7) - pkin(14);
t202 = sin(t210);
t185 = t198 + t202 / 0.2e1;
t199 = cos(t210) / 0.2e1;
t203 = cos(t209);
t187 = t199 + t203 / 0.2e1;
t194 = t201 - t205 / 0.2e1;
t213 = sin(pkin(14));
t167 = t185 * t220 + t187 * t190 - t194 * t213;
t214 = sin(pkin(8));
t218 = cos(pkin(8));
t158 = -t167 * t214 + t179 * t218;
t191 = t200 - t204 / 0.2e1;
t229 = cos(qJ(2));
t183 = -t225 * t191 + t229 * t230;
t163 = t182 * t187 - t183 * t213 + t185 * t256;
t154 = -t163 * t214 + t174 * t218;
t180 = t193 * t230 - t225 * t224;
t181 = t191 * t230 + t225 * t229;
t255 = t216 * t230;
t161 = t180 * t187 - t181 * t213 - t185 * t255;
t259 = t180 * t215;
t173 = -t219 * t255 - t259;
t153 = -t161 * t214 + t173 * t218;
t266 = t220 * pkin(10) + pkin(9);
t253 = t230 * pkin(1) + pkin(10) * t256;
t252 = pkin(8) - qJ(4);
t251 = pkin(8) + qJ(4);
t248 = cos(t251);
t247 = sin(t252);
t246 = -m(7) * pkin(13) + mrSges(6,2) - mrSges(7,3);
t245 = cos(t252) / 0.2e1;
t244 = sin(t251) / 0.2e1;
t207 = t225 * pkin(1);
t243 = t181 * pkin(2) - qJ(3) * t259 + t207;
t242 = t194 * pkin(2) + t179 * qJ(3) + t266;
t241 = t183 * pkin(2) + t174 * qJ(3) + t253;
t221 = sin(qJ(6));
t226 = cos(qJ(6));
t240 = -m(7) * pkin(5) - mrSges(7,1) * t226 + mrSges(7,2) * t221 - mrSges(6,1);
t239 = t245 + t248 / 0.2e1;
t238 = t244 + t247 / 0.2e1;
t237 = -mrSges(7,1) * t221 - mrSges(7,2) * t226 - pkin(12) * t267 + mrSges(5,2) - mrSges(6,3);
t186 = t198 - t202 / 0.2e1;
t188 = t199 - t203 / 0.2e1;
t217 = cos(pkin(14));
t162 = t180 * t186 + t181 * t217 - t188 * t255;
t236 = t162 * pkin(3) + t153 * pkin(11) + t243;
t168 = t186 * t190 + t188 * t220 + t194 * t217;
t234 = t168 * pkin(3) + t158 * pkin(11) + t242;
t164 = t182 * t186 + t183 * t217 + t188 * t256;
t233 = t164 * pkin(3) + t154 * pkin(11) + t241;
t228 = cos(qJ(4));
t227 = cos(qJ(5));
t223 = sin(qJ(4));
t222 = sin(qJ(5));
t192 = t245 - t248 / 0.2e1;
t189 = t244 - t247 / 0.2e1;
t150 = t167 * t189 + t168 * t228 + t179 * t192;
t147 = t163 * t189 + t164 * t228 + t174 * t192;
t145 = t161 * t189 + t162 * t228 + t173 * t192;
t1 = (-mrSges(1,3) - m(2) * pkin(9) - mrSges(2,3) - m(3) * t266 - t194 * mrSges(3,1) - t190 * mrSges(3,2) - t220 * mrSges(3,3) - m(4) * t242 - t168 * mrSges(4,1) - t167 * mrSges(4,2) - t179 * mrSges(4,3) - m(5) * t234 - t150 * mrSges(5,1) - t158 * mrSges(5,3) + t240 * (t150 * t227 + t158 * t222) + t237 * (-t167 * t239 + t168 * t223 - t179 * t238) + t246 * (t150 * t222 - t158 * t227) + t267 * (-t150 * pkin(4) - t234)) * g(3) + (-mrSges(1,2) - t225 * mrSges(2,1) - m(3) * t207 - t181 * mrSges(3,1) - t180 * mrSges(3,2) - m(4) * t243 - t162 * mrSges(4,1) - t161 * mrSges(4,2) - t173 * mrSges(4,3) - m(5) * t236 - t145 * mrSges(5,1) - t153 * mrSges(5,3) + t240 * (t145 * t227 + t153 * t222) + t237 * (-t161 * t239 + t162 * t223 - t173 * t238) + t246 * (t145 * t222 - t153 * t227) + (-mrSges(2,2) + (m(3) * pkin(10) + mrSges(3,3) + (m(4) + m(5) + t267) * (qJ(3) * t219 + pkin(10))) * t216) * t230 + t267 * (-t145 * pkin(4) - t236)) * g(2) + (-mrSges(1,1) - t230 * mrSges(2,1) - m(3) * t253 - t183 * mrSges(3,1) - t182 * mrSges(3,2) - m(4) * t241 - t164 * mrSges(4,1) - t163 * mrSges(4,2) - t174 * mrSges(4,3) - m(5) * t233 - t147 * mrSges(5,1) - t154 * mrSges(5,3) + (-t216 * mrSges(3,3) + mrSges(2,2)) * t225 + t240 * (t147 * t227 + t154 * t222) + t237 * (-t163 * t239 + t164 * t223 - t174 * t238) + t246 * (t147 * t222 - t154 * t227) + t267 * (-t147 * pkin(4) - t233)) * g(1);
U  = t1;
