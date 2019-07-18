% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:10
% EndTime: 2019-07-18 17:17:27
% DurationCPUTime: 5.76s
% Computational Cost: add. (6570->471), mult. (16507->688), div. (0->0), fcn. (12376->8), ass. (0->213)
t186 = sin(qJ(3));
t187 = sin(qJ(2));
t190 = cos(qJ(3));
t191 = cos(qJ(2));
t171 = t186 * t191 + t187 * t190;
t165 = t171 * qJD(1);
t183 = qJD(2) + qJD(3);
t185 = sin(qJ(4));
t189 = cos(qJ(4));
t148 = -t165 * t185 + t183 * t189;
t250 = pkin(1) * qJD(2);
t223 = t190 * t250;
t174 = -pkin(2) * t183 - t223;
t124 = -pkin(3) * t148 + t174;
t145 = t183 * t171;
t135 = t145 * qJD(1);
t130 = Ifges(6,3) * t135;
t149 = t165 * t189 + t183 * t185;
t184 = sin(qJ(5));
t188 = cos(qJ(5));
t215 = t188 * t148 - t149 * t184;
t98 = t148 * t184 + t149 * t188;
t261 = Ifges(6,4) * t98;
t231 = qJD(1) * t191;
t232 = qJD(1) * t187;
t164 = -t186 * t232 + t190 * t231;
t160 = qJD(4) - t164;
t157 = qJD(5) + t160;
t271 = -t157 / 0.2e1;
t279 = -t98 / 0.2e1;
t225 = pkin(1) * t231;
t120 = -pkin(2) * t164 - pkin(5) * t165 - t225;
t182 = t186 * t250;
t173 = pkin(5) * t183 + t182;
t100 = t120 * t185 + t173 * t189;
t243 = t100 * t184;
t99 = t189 * t120 - t173 * t185;
t76 = pkin(3) * t160 + t99;
t41 = t188 * t76 - t243;
t242 = t100 * t188;
t42 = t184 * t76 + t242;
t315 = t130 + (Ifges(6,5) * t215 - Ifges(6,6) * t98) * t271 + (t215 * t41 + t42 * t98) * mrSges(6,3) - t124 * (mrSges(6,1) * t98 + mrSges(6,2) * t215) + (Ifges(6,1) * t215 - t261) * t279;
t314 = t41 * mrSges(6,1);
t313 = t42 * mrSges(6,2);
t312 = t149 * Ifges(5,5) + t98 * Ifges(6,5) + t148 * Ifges(5,6) + Ifges(6,6) * t215 + t160 * Ifges(5,3) + t157 * Ifges(6,3);
t235 = t184 * t185;
t199 = -t188 * t189 + t235;
t294 = qJD(4) + qJD(5);
t311 = t294 * t199;
t234 = t185 * t188;
t170 = t184 * t189 + t234;
t113 = t170 * t164;
t143 = t294 * t170;
t310 = t113 - t143;
t114 = t199 * t164;
t309 = t114 - t311;
t169 = t186 * t187 - t190 * t191;
t144 = t183 * t169;
t134 = t144 * qJD(1);
t87 = qJD(4) * t148 - t134 * t189;
t88 = -qJD(4) * t149 + t134 * t185;
t18 = qJD(5) * t215 + t184 * t88 + t188 * t87;
t19 = -qJD(5) * t98 - t184 * t87 + t188 * t88;
t249 = pkin(1) * qJD(3);
t220 = qJD(2) * t249;
t213 = t190 * t220;
t230 = qJD(2) * t187;
t224 = pkin(1) * t230;
t69 = pkin(2) * t135 + pkin(5) * t134 + qJD(1) * t224;
t35 = -qJD(4) * t100 - t185 * t213 + t189 * t69;
t20 = pkin(3) * t135 + t35;
t34 = qJD(4) * t99 + t185 * t69 + t189 * t213;
t3 = qJD(5) * t41 + t184 * t20 + t188 * t34;
t4 = -qJD(5) * t42 - t184 * t34 + t188 * t20;
t308 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t18 + Ifges(6,6) * t19;
t93 = Ifges(6,4) * t215;
t307 = -Ifges(6,2) * t98 + t93;
t289 = t18 / 0.2e1;
t288 = t19 / 0.2e1;
t276 = t135 / 0.2e1;
t137 = pkin(2) * t165 - pkin(5) * t164;
t119 = pkin(1) * t232 + t137;
t178 = pkin(1) * t186 + pkin(5);
t222 = t190 * t249;
t159 = t165 * pkin(3);
t94 = t119 * t189 + t159;
t305 = -t119 * t234 - t143 * t178 - t184 * t94 - t199 * t222;
t304 = t119 * t235 - t170 * t222 + t178 * t311 - t188 * t94;
t111 = t137 * t185 + t189 * t223;
t110 = t189 * t137 - t185 * t223;
t90 = t110 + t159;
t299 = -t143 * pkin(5) - t111 * t188 - t184 * t90;
t298 = pkin(5) * t311 + t111 * t184 - t188 * t90;
t50 = -mrSges(6,1) * t215 + mrSges(6,2) * t98;
t297 = m(6) * t124 + t50;
t210 = mrSges(5,1) * t185 + mrSges(5,2) * t189;
t296 = t174 * t210;
t229 = qJD(4) * t185;
t221 = pkin(3) * t229;
t238 = t164 * t185;
t226 = pkin(3) * t238;
t295 = t186 * t249 + t221 - t226;
t292 = t35 * mrSges(5,1) - t34 * mrSges(5,2) + Ifges(5,5) * t87 + Ifges(5,6) * t88 + t308;
t291 = Ifges(6,4) * t289 + Ifges(6,2) * t288 + Ifges(6,6) * t276;
t290 = Ifges(6,1) * t289 + Ifges(6,4) * t288 + Ifges(6,5) * t276;
t39 = Ifges(6,2) * t215 + Ifges(6,6) * t157 + t261;
t287 = -t39 / 0.2e1;
t286 = t39 / 0.2e1;
t40 = Ifges(6,1) * t98 + Ifges(6,5) * t157 + t93;
t285 = -t40 / 0.2e1;
t284 = t40 / 0.2e1;
t283 = t87 / 0.2e1;
t282 = t88 / 0.2e1;
t281 = -t215 / 0.2e1;
t280 = t215 / 0.2e1;
t278 = t98 / 0.2e1;
t274 = -t148 / 0.2e1;
t273 = -t149 / 0.2e1;
t272 = t149 / 0.2e1;
t270 = t157 / 0.2e1;
t269 = -t160 / 0.2e1;
t266 = t165 / 0.2e1;
t265 = -t185 / 0.2e1;
t264 = t189 / 0.2e1;
t263 = m(4) * pkin(1) ^ 2;
t260 = pkin(1) * t190;
t259 = pkin(1) * t191;
t258 = pkin(3) * t189;
t256 = mrSges(4,3) * t164;
t255 = mrSges(4,3) * t165;
t254 = Ifges(3,4) * t187;
t253 = Ifges(5,4) * t149;
t252 = Ifges(5,4) * t185;
t251 = Ifges(5,4) * t189;
t248 = t165 * Ifges(4,4);
t247 = t185 * t35;
t246 = t189 * t34;
t245 = Ifges(3,5) * qJD(2);
t244 = Ifges(3,6) * qJD(2);
t241 = t144 * t185;
t240 = t144 * t189;
t237 = t164 * t189;
t233 = -mrSges(4,1) * t183 - mrSges(5,1) * t148 + mrSges(5,2) * t149 + t255;
t228 = qJD(4) * t189;
t227 = pkin(3) * t241;
t181 = Ifges(3,4) * t231;
t80 = Ifges(5,2) * t148 + Ifges(5,6) * t160 + t253;
t219 = t80 * t265;
t147 = Ifges(5,4) * t148;
t81 = Ifges(5,1) * t149 + Ifges(5,5) * t160 + t147;
t218 = t81 * t264;
t180 = -pkin(2) - t258;
t217 = t245 / 0.2e1;
t216 = -t244 / 0.2e1;
t214 = t186 * t220;
t140 = pkin(2) * t169 - pkin(5) * t171 - t259;
t86 = pkin(2) * t145 + pkin(5) * t144 + t224;
t212 = pkin(3) * t145 + t189 * t86 + (-qJD(5) * t185 - t229) * t140;
t211 = mrSges(5,1) * t189 - mrSges(5,2) * t185;
t209 = Ifges(5,1) * t189 - t252;
t208 = Ifges(5,1) * t185 + t251;
t207 = -Ifges(5,2) * t185 + t251;
t206 = Ifges(5,2) * t189 + t252;
t205 = Ifges(5,5) * t189 - Ifges(5,6) * t185;
t204 = Ifges(5,5) * t185 + Ifges(5,6) * t189;
t203 = -t185 * t34 - t189 * t35;
t55 = mrSges(5,1) * t135 - mrSges(5,3) * t87;
t56 = -mrSges(5,2) * t135 + mrSges(5,3) * t88;
t202 = -t185 * t55 + t189 * t56;
t201 = t100 * t189 - t185 * t99;
t200 = t100 * t185 + t189 * t99;
t106 = pkin(3) * t169 + t140 * t189;
t197 = qJD(5) * t106 + t140 * t228 + t185 * t86;
t196 = -qJD(4) * t200 - t247;
t107 = -mrSges(5,2) * t160 + mrSges(5,3) * t148;
t108 = mrSges(5,1) * t160 - mrSges(5,3) * t149;
t195 = m(5) * t200 + t185 * t107 + t189 * t108;
t121 = t164 * Ifges(4,2) + t183 * Ifges(4,6) + t248;
t158 = Ifges(4,4) * t164;
t122 = t165 * Ifges(4,1) + t183 * Ifges(4,5) + t158;
t25 = t87 * Ifges(5,4) + t88 * Ifges(5,2) + t135 * Ifges(5,6);
t26 = t87 * Ifges(5,1) + t88 * Ifges(5,4) + t135 * Ifges(5,5);
t77 = -pkin(3) * t88 + t214;
t193 = t185 * t26 / 0.2e1 - t183 * (Ifges(4,5) * t164 - Ifges(4,6) * t165) / 0.2e1 - t164 * t296 - t81 * t237 / 0.2e1 + t80 * t238 / 0.2e1 + (-Ifges(6,5) * t114 - Ifges(6,6) * t113 + Ifges(6,3) * t165) * t271 + (-Ifges(6,1) * t114 - Ifges(6,4) * t113 + Ifges(6,5) * t165) * t279 + (-Ifges(6,4) * t114 - Ifges(6,2) * t113 + Ifges(6,6) * t165) * t281 - t100 * (-mrSges(5,2) * t165 - mrSges(5,3) * t238) - t99 * (mrSges(5,1) * t165 - mrSges(5,3) * t237) + (-Ifges(6,5) * t311 - Ifges(6,6) * t143) * t270 + (-Ifges(6,1) * t311 - Ifges(6,4) * t143) * t278 + (-Ifges(6,4) * t311 - Ifges(6,2) * t143) * t280 - t311 * t284 - t211 * t214 + (t218 + t219 + t296) * qJD(4) - Ifges(4,5) * t134 - Ifges(4,6) * t135 + t165 * t313 + mrSges(5,3) * t246 + (mrSges(4,1) * t165 + mrSges(4,2) * t164) * t225 - (Ifges(4,1) * t164 - t248 + t312) * t165 / 0.2e1 - t165 * t314 + (-t4 * t170 - t199 * t3 - t309 * t41 + t310 * t42) * mrSges(6,3) + (-mrSges(6,1) * t310 + mrSges(6,2) * t309) * t124 + t223 * t256 + t25 * t264 + t121 * t266 + (Ifges(5,3) * t165 + t164 * t205) * t269 + (Ifges(5,5) * t165 + t164 * t209) * t273 + (Ifges(5,6) * t165 + t164 * t207) * t274 + t206 * t282 + t208 * t283 - t114 * t285 - t143 * t286 - t113 * t287 + t170 * t290 + (t148 * t207 + t149 * t209 + t160 * t205) * qJD(4) / 0.2e1 - (-Ifges(4,2) * t165 + t122 + t158) * t164 / 0.2e1 + t77 * (mrSges(6,1) * t199 + mrSges(6,2) * t170) + (Ifges(6,4) * t170 - Ifges(6,2) * t199) * t288 + (Ifges(6,1) * t170 - Ifges(6,4) * t199) * t289 + (Ifges(6,5) * t170 - Ifges(6,6) * t199 + t204) * t276 - t199 * t291;
t179 = -pkin(2) - t260;
t175 = t180 - t260;
t167 = t199 * pkin(5);
t166 = t170 * pkin(5);
t163 = Ifges(3,1) * t232 + t181 + t245;
t162 = t244 + (t191 * Ifges(3,2) + t254) * qJD(1);
t154 = -mrSges(4,2) * t183 + t256;
t153 = t199 * t178;
t152 = t170 * t178;
t150 = t182 + t226;
t136 = -mrSges(4,1) * t164 + mrSges(4,2) * t165;
t131 = Ifges(5,3) * t135;
t127 = t199 * t171;
t126 = t170 * t171;
t71 = mrSges(6,1) * t157 - mrSges(6,3) * t98;
t70 = -mrSges(6,2) * t157 + mrSges(6,3) * t215;
t62 = t106 * t184 + t140 * t234;
t61 = t106 * t188 - t140 * t235;
t49 = t188 * t99 - t243;
t48 = -t184 * t99 - t242;
t44 = -mrSges(5,1) * t88 + mrSges(5,2) * t87;
t37 = t144 * t170 + t171 * t311;
t36 = -t143 * t171 + t144 * t199;
t13 = -mrSges(6,2) * t135 + mrSges(6,3) * t19;
t12 = mrSges(6,1) * t135 - mrSges(6,3) * t18;
t11 = -t184 * t197 + t188 * t212;
t10 = t184 * t212 + t188 * t197;
t9 = -mrSges(6,1) * t19 + mrSges(6,2) * t18;
t1 = [t183 * (-Ifges(4,5) * t144 - Ifges(4,6) * t145) / 0.2e1 + ((t144 * t190 - t145 * t186) * pkin(1) * mrSges(4,3) + (t163 / 0.2e1 + t217 + 0.3e1 / 0.2e1 * t181) * t191) * qJD(2) + t164 * (-Ifges(4,4) * t144 - Ifges(4,2) * t145) / 0.2e1 + (-Ifges(4,1) * t144 - Ifges(4,4) * t145) * t266 + t99 * (mrSges(5,1) * t145 + mrSges(5,3) * t240) + (t185 * t56 + t189 * t55 + m(5) * (t100 * t228 - t229 * t99 - t203) - t108 * t229 + t107 * t228) * t140 + (-t126 * t3 + t127 * t4 - t36 * t41 + t37 * t42) * mrSges(6,3) + t77 * (mrSges(6,1) * t126 - mrSges(6,2) * t127) + (-Ifges(6,5) * t127 - Ifges(6,6) * t126) * t276 + (-Ifges(6,4) * t127 - Ifges(6,2) * t126) * t288 + (-Ifges(6,1) * t127 - Ifges(6,4) * t126) * t289 + t174 * (-mrSges(5,1) * t241 - mrSges(5,2) * t240) + t160 * (-Ifges(5,5) * t240 + Ifges(5,6) * t241 + Ifges(5,3) * t145) / 0.2e1 + t100 * (-mrSges(5,2) * t145 + mrSges(5,3) * t241) + t148 * (-Ifges(5,4) * t240 + Ifges(5,2) * t241 + Ifges(5,6) * t145) / 0.2e1 + (-Ifges(5,1) * t240 + Ifges(5,4) * t241 + Ifges(5,5) * t145) * t272 + (-mrSges(4,1) * t135 + mrSges(4,2) * t134 - qJD(1) * (mrSges(4,1) * t145 - mrSges(4,2) * t144)) * t259 + (t205 * t276 + t209 * t283 + t207 * t282 - Ifges(4,4) * t135 - Ifges(4,1) * t134 + t25 * t265 + t26 * t264 + (m(6) * t77 + t9) * t185 * pkin(3) + (mrSges(4,3) + t210) * t214 + t203 * mrSges(5,3) + (t174 * t211 + t204 * t269 + t206 * t274 + t208 * t273 + t81 * t265 - t189 * t80 / 0.2e1 + t297 * t258 - t201 * mrSges(5,3)) * qJD(4)) * t171 + (-mrSges(4,3) * t213 + t130 / 0.2e1 + t131 / 0.2e1 + Ifges(4,4) * t134 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(4,2)) * t135 + t292) * t169 + m(6) * (t10 * t42 + t11 * t41 - t124 * t227 + t3 * t62 + t4 * t61) + t195 * t86 - t145 * t121 / 0.2e1 - t144 * t122 / 0.2e1 + t124 * (-mrSges(6,1) * t37 + mrSges(6,2) * t36) + t10 * t70 + t11 * t71 + t61 * t12 + t62 * t13 + t145 * t314 + (-t162 / 0.2e1 + pkin(1) * t136 + t216 + (pkin(1) * (mrSges(4,1) * t169 + mrSges(4,2) * t171) - 0.3e1 / 0.2e1 * t254 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.2e1 * t263) * t191) * qJD(1)) * t230 - t50 * t227 + t312 * t145 / 0.2e1 - t145 * t313 - t144 * t218 - t144 * t219 + (Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t145) * t270 + (Ifges(6,1) * t36 + Ifges(6,4) * t37 + Ifges(6,5) * t145) * t278 + (Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t145) * t280 + t36 * t284 + t37 * t286 - t127 * t290 - t126 * t291; (m(5) * (-t100 * t229 - t228 * t99 + t246 - t247) - t107 * t229 - t108 * t228 + t202) * t178 - t195 * t119 + (-t136 * t232 + (t190 * t134 + (qJD(2) * t165 - t135) * t186) * mrSges(4,3) + ((-qJD(2) * mrSges(4,1) + m(5) * (qJD(2) * t179 + t174) + t233) * t186 + (m(5) * t201 - qJD(2) * mrSges(4,2) + t189 * t107 - t185 * t108 + t154) * t190) * qJD(3)) * pkin(1) + t304 * t71 + t193 + t305 * t70 + t295 * t50 + t175 * t9 + t179 * t44 - t152 * t12 - t153 * t13 + t196 * mrSges(5,3) + ((t217 - t163 / 0.2e1 - t181 / 0.2e1) * t191 + (t216 + t162 / 0.2e1 + Ifges(3,4) * t232 / 0.2e1 + (t263 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t231) * t187) * qJD(1) + (t124 * t295 - t152 * t4 - t153 * t3 + t175 * t77 + t304 * t41 + t305 * t42) * m(6); ((-mrSges(4,2) * qJD(3) - t154) * t190 + (-mrSges(4,1) * qJD(3) - t233 + t255) * t186) * t250 + t298 * t71 + t299 * t70 + t193 + t180 * t9 - t166 * t12 - t167 * t13 - t150 * t50 - t110 * t108 - t111 * t107 - pkin(2) * t44 - mrSges(5,3) * t247 + t202 * pkin(5) + ((-t99 * mrSges(5,3) - pkin(5) * t108) * t189 + (-mrSges(5,3) * t100 + pkin(3) * t50 - pkin(5) * t107) * t185) * qJD(4) + (-t166 * t4 - t167 * t3 + t180 * t77 + t299 * t42 + t298 * t41 + (-t150 + t221) * t124) * m(6) + (-pkin(2) * t214 - t100 * t111 - t110 * t99 - t174 * t182 + (t196 + t246) * pkin(5)) * m(5); t131 + (t188 * t12 + t184 * t13 + m(6) * (t184 * t3 + t188 * t4) - t297 * t149 + (-t184 * t71 + t188 * t70 + m(6) * (-t184 * t41 + t188 * t42)) * qJD(5)) * pkin(3) - m(6) * (t41 * t48 + t42 * t49) + t215 * t285 + (-Ifges(5,2) * t149 + t147 + t81) * t274 + (t100 * t149 + t148 * t99) * mrSges(5,3) - t174 * (mrSges(5,1) * t149 + mrSges(5,2) * t148) + t292 - t99 * t107 + t100 * t108 - t49 * t70 - t48 * t71 + (Ifges(5,5) * t148 - Ifges(5,6) * t149) * t269 + t80 * t272 + (Ifges(5,1) * t148 - t253) * t273 + t307 * t281 - t98 * t287 + t315; t39 * t278 - t41 * t70 + t42 * t71 + (t307 + t40) * t281 + t308 + t315;];
tauc  = t1(:);
