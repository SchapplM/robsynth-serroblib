% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:50
% EndTime: 2019-12-31 21:53:07
% DurationCPUTime: 7.04s
% Computational Cost: add. (5541->453), mult. (13946->595), div. (0->0), fcn. (9306->6), ass. (0->217)
t362 = Ifges(5,4) + Ifges(6,4);
t363 = Ifges(5,1) + Ifges(6,1);
t353 = Ifges(5,5) + Ifges(6,5);
t361 = Ifges(5,2) + Ifges(6,2);
t352 = Ifges(6,6) + Ifges(5,6);
t214 = cos(qJ(4));
t365 = t362 * t214;
t211 = sin(qJ(4));
t364 = t362 * t211;
t342 = -t352 * t211 + t353 * t214;
t340 = -t361 * t211 + t365;
t338 = t363 * t214 - t364;
t210 = qJD(2) + qJD(3);
t360 = t210 * Ifges(4,6) / 0.2e1;
t212 = sin(qJ(3));
t213 = sin(qJ(2));
t215 = cos(qJ(3));
t216 = cos(qJ(2));
t187 = t212 * t216 + t215 * t213;
t178 = t187 * qJD(1);
t162 = -t178 * t211 + t210 * t214;
t359 = t362 * t162;
t186 = t212 * t213 - t215 * t216;
t177 = t186 * qJD(1);
t282 = t177 * t211;
t358 = -qJ(5) * t282 + t214 * qJD(5);
t163 = t178 * t214 + t210 * t211;
t357 = t362 * t163;
t327 = -pkin(7) - pkin(6);
t201 = t327 * t216;
t190 = qJD(1) * t201;
t179 = t212 * t190;
t200 = t327 * t213;
t189 = qJD(1) * t200;
t182 = qJD(2) * pkin(2) + t189;
t149 = t182 * t215 + t179;
t206 = -pkin(2) * t216 - pkin(1);
t199 = qJD(1) * t206;
t291 = t210 * Ifges(4,5);
t356 = t199 * mrSges(4,2) - t149 * mrSges(4,3) + t291 / 0.2e1;
t292 = t178 * Ifges(4,4);
t355 = t360 + t292 / 0.2e1 - t177 * Ifges(4,2) / 0.2e1;
t158 = t210 * t187;
t143 = t158 * qJD(1);
t157 = t210 * t186;
t142 = t157 * qJD(1);
t87 = qJD(4) * t162 - t142 * t214;
t88 = -qJD(4) * t163 + t142 * t211;
t351 = t352 * t143 + t361 * t88 + t362 * t87;
t350 = t353 * t143 + t362 * t88 + t363 * t87;
t174 = qJD(4) + t177;
t335 = t361 * t162 + t352 * t174 + t357;
t334 = t363 * t163 + t353 * t174 + t359;
t209 = t214 * qJ(5);
t247 = t178 * pkin(4) + t177 * t209;
t203 = pkin(2) * t212 + pkin(8);
t276 = -qJ(5) - t203;
t252 = qJD(4) * t276;
t300 = pkin(2) * qJD(3);
t268 = t215 * t300;
t146 = pkin(3) * t178 + pkin(8) * t177;
t275 = qJD(1) * t213;
t121 = pkin(2) * t275 + t146;
t152 = t189 * t215 + t179;
t62 = t214 * t121 - t152 * t211;
t348 = -t247 - t62 + (-qJD(5) - t268) * t211 + t214 * t252;
t22 = -mrSges(5,1) * t88 + mrSges(5,2) * t87;
t180 = t215 * t190;
t150 = t182 * t212 - t180;
t262 = qJD(2) * t327;
t251 = qJD(1) * t262;
t183 = t213 * t251;
t229 = t216 * t251;
t77 = qJD(3) * t150 + t183 * t212 - t215 * t229;
t347 = m(5) * t77 + t22;
t63 = t211 * t121 + t214 * t152;
t346 = t211 * t252 + t214 * t268 + t358 - t63;
t151 = t189 * t212 - t180;
t168 = pkin(4) * t282;
t272 = qJD(4) * t211;
t267 = pkin(4) * t272;
t345 = t212 * t300 - t151 + t168 + t267;
t126 = -pkin(3) * t210 - t149;
t243 = mrSges(6,1) * t211 + mrSges(6,2) * t214;
t245 = mrSges(5,1) * t211 + mrSges(5,2) * t214;
t86 = -pkin(4) * t162 + qJD(5) + t126;
t344 = t126 * t245 + t86 * t243;
t343 = t353 * t211 + t352 * t214;
t341 = t361 * t214 + t364;
t339 = t363 * t211 + t365;
t319 = t163 / 0.2e1;
t337 = t344 + t340 * t162 / 0.2e1 + t338 * t319 + t342 * t174 / 0.2e1;
t114 = t177 * pkin(3) - t178 * pkin(8) + t199;
t127 = pkin(8) * t210 + t150;
t42 = t214 * t114 - t127 * t211;
t24 = -qJ(5) * t163 + t42;
t23 = pkin(4) * t174 + t24;
t43 = t114 * t211 + t127 * t214;
t25 = qJ(5) * t162 + t43;
t336 = -t199 * mrSges(4,1) - t42 * mrSges(5,1) - t23 * mrSges(6,1) + t43 * mrSges(5,2) + t25 * mrSges(6,2) + t355;
t306 = -qJ(5) - pkin(8);
t253 = qJD(4) * t306;
t67 = t214 * t146 - t149 * t211;
t333 = -qJD(5) * t211 + t214 * t253 - t247 - t67;
t68 = t211 * t146 + t214 * t149;
t332 = t211 * t253 + t358 - t68;
t148 = t186 * pkin(3) - t187 * pkin(8) + t206;
t165 = t200 * t212 - t201 * t215;
t159 = t214 * t165;
t90 = t211 * t148 + t159;
t331 = t215 * t200 + t201 * t212;
t330 = -t211 * t42 + t214 * t43;
t329 = t87 / 0.2e1;
t328 = t88 / 0.2e1;
t325 = pkin(1) * mrSges(3,1);
t324 = pkin(1) * mrSges(3,2);
t323 = t143 / 0.2e1;
t322 = -t162 / 0.2e1;
t320 = -t163 / 0.2e1;
t318 = -t174 / 0.2e1;
t316 = t177 / 0.2e1;
t314 = -t211 / 0.2e1;
t311 = t214 / 0.2e1;
t310 = m(4) * t199;
t309 = pkin(2) * t215;
t308 = pkin(8) * t214;
t271 = qJD(4) * t214;
t273 = qJD(2) * t213;
t269 = pkin(2) * t273;
t52 = pkin(3) * t143 + pkin(8) * t142 + qJD(1) * t269;
t76 = qJD(3) * t149 + t215 * t183 + t212 * t229;
t7 = t114 * t271 - t127 * t272 + t211 * t52 + t214 * t76;
t307 = t214 * t7;
t305 = Ifges(3,4) * t213;
t170 = Ifges(4,4) * t177;
t296 = t331 * t77;
t294 = t178 * mrSges(4,3);
t293 = t178 * Ifges(4,1);
t287 = -mrSges(4,1) * t210 - mrSges(5,1) * t162 + mrSges(5,2) * t163 + t294;
t286 = Ifges(3,5) * qJD(2);
t285 = Ifges(3,6) * qJD(2);
t284 = qJD(2) * mrSges(3,1);
t283 = qJD(2) * mrSges(3,2);
t281 = t177 * t214;
t279 = t187 * t211;
t277 = t203 * t214;
t274 = qJD(1) * t216;
t85 = pkin(3) * t158 + pkin(8) * t157 + t269;
t192 = t213 * t262;
t193 = t216 * t262;
t92 = qJD(3) * t331 + t192 * t215 + t193 * t212;
t270 = t148 * t271 + t211 * t85 + t214 * t92;
t266 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t265 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t264 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t205 = -pkin(4) * t214 - pkin(3);
t261 = t187 * t271;
t260 = t286 / 0.2e1;
t259 = -t285 / 0.2e1;
t21 = -t88 * mrSges(6,1) + t87 * mrSges(6,2);
t254 = -t211 * t92 + t214 * t85;
t89 = t214 * t148 - t165 * t211;
t8 = -qJD(4) * t43 - t211 * t76 + t214 * t52;
t1 = pkin(4) * t143 - qJ(5) * t87 - qJD(5) * t163 + t8;
t250 = -t8 * mrSges(5,3) - t1 * mrSges(6,3);
t249 = -t42 * mrSges(5,3) - t23 * mrSges(6,3);
t248 = -t43 * mrSges(5,3) - t25 * mrSges(6,3);
t246 = mrSges(5,1) * t214 - mrSges(5,2) * t211;
t244 = mrSges(6,1) * t214 - mrSges(6,2) * t211;
t230 = -t211 * t43 - t214 * t42;
t228 = qJ(5) * t157 - qJD(5) * t187;
t3 = qJ(5) * t88 + qJD(5) * t162 + t7;
t219 = t8 * mrSges(5,1) + t1 * mrSges(6,1) - t7 * mrSges(5,2) - t3 * mrSges(6,2);
t93 = qJD(3) * t165 + t192 * t212 - t215 * t193;
t218 = m(5) * (qJD(4) * t230 - t8 * t211 + t307);
t123 = -t170 + t291 + t293;
t20 = -pkin(4) * t88 + t77;
t69 = t163 * Ifges(6,5) + t162 * Ifges(6,6) + t174 * Ifges(6,3);
t70 = t163 * Ifges(5,5) + t162 * Ifges(5,6) + t174 * Ifges(5,3);
t217 = -t76 * mrSges(4,2) - Ifges(4,5) * t142 - Ifges(4,6) * t143 - t20 * t244 + (-t246 - mrSges(4,1)) * t77 + t339 * t329 + t341 * t328 + t343 * t323 + (t123 - t170) * t316 + t350 * t211 / 0.2e1 + t351 * t311 + (t214 * t3 - t23 * t281 - t25 * t282) * mrSges(6,3) + (-t281 * t42 - t282 * t43 + t307) * mrSges(5,3) + (-Ifges(4,2) * t316 + t360 + t352 * t322 + t353 * t320 + (Ifges(6,3) + Ifges(5,3)) * t318 + t336) * t178 + (-t342 * t318 - t338 * t320 - t340 * t322 + t344 + t356) * t177 - (-Ifges(4,1) * t177 - t292 + t69 + t70) * t178 / 0.2e1 + t337 * qJD(4) + (-t282 / 0.2e1 - t272 / 0.2e1) * t335 + (t281 / 0.2e1 + t271 / 0.2e1) * t334;
t207 = Ifges(3,4) * t274;
t198 = t209 + t308;
t197 = t306 * t211;
t196 = mrSges(3,3) * t274 - t283;
t195 = -mrSges(3,3) * t275 + t284;
t194 = t205 - t309;
t185 = t209 + t277;
t184 = t276 * t211;
t176 = Ifges(3,1) * t275 + t207 + t286;
t175 = t285 + (t216 * Ifges(3,2) + t305) * qJD(1);
t166 = -mrSges(4,2) * t210 - mrSges(4,3) * t177;
t145 = mrSges(4,1) * t177 + mrSges(4,2) * t178;
t139 = Ifges(5,3) * t143;
t138 = Ifges(6,3) * t143;
t115 = pkin(4) * t279 - t331;
t108 = mrSges(5,1) * t174 - mrSges(5,3) * t163;
t107 = mrSges(6,1) * t174 - mrSges(6,3) * t163;
t106 = -mrSges(5,2) * t174 + mrSges(5,3) * t162;
t105 = -mrSges(6,2) * t174 + mrSges(6,3) * t162;
t98 = t150 - t168;
t96 = -mrSges(6,1) * t162 + mrSges(6,2) * t163;
t84 = Ifges(5,5) * t87;
t83 = Ifges(6,5) * t87;
t82 = Ifges(5,6) * t88;
t81 = Ifges(6,6) * t88;
t50 = -qJ(5) * t279 + t90;
t34 = pkin(4) * t186 - t187 * t209 + t89;
t32 = -mrSges(5,2) * t143 + mrSges(5,3) * t88;
t31 = -mrSges(6,2) * t143 + mrSges(6,3) * t88;
t30 = mrSges(5,1) * t143 - mrSges(5,3) * t87;
t29 = mrSges(6,1) * t143 - mrSges(6,3) * t87;
t28 = (-t211 * t157 + t261) * pkin(4) + t93;
t10 = -qJD(4) * t90 + t254;
t9 = -t165 * t272 + t270;
t5 = -qJ(5) * t261 + (-qJD(4) * t165 + t228) * t211 + t270;
t4 = pkin(4) * t158 + t228 * t214 + (-t159 + (qJ(5) * t187 - t148) * t211) * qJD(4) + t254;
t2 = [t34 * t29 + t50 * t31 + t89 * t30 + t90 * t32 + t28 * t96 + t5 * t105 + t9 * t106 + t4 * t107 + t10 * t108 + t115 * t21 - t331 * t22 + t92 * t166 + t206 * (mrSges(4,1) * t143 - mrSges(4,2) * t142) + t287 * t93 + (t142 * t331 - t143 * t165) * mrSges(4,3) + m(4) * (-t149 * t93 + t150 * t92 + t165 * t76 - t296) + m(6) * (t1 * t34 + t115 * t20 + t23 * t4 + t25 * t5 + t28 * t86 + t3 * t50) + m(5) * (t10 * t42 + t126 * t93 + t43 * t9 + t7 * t90 + t8 * t89 - t296) + (t176 / 0.2e1 - pkin(6) * t195 + t260 + (-0.2e1 * t324 + 0.3e1 / 0.2e1 * Ifges(3,4) * t216) * qJD(1)) * t216 * qJD(2) + (t69 / 0.2e1 + t70 / 0.2e1 - t150 * mrSges(4,3) + t264 * t174 + t266 * t163 + t265 * t162 - t336 - t355) * t158 - (t337 + (-t211 * t25 - t214 * t23) * mrSges(6,3) + t230 * mrSges(5,3) - t170 / 0.2e1 + t293 / 0.2e1 + t123 / 0.2e1 + t334 * t311 + t335 * t314 + t356) * t157 + (-t76 * mrSges(4,3) + t83 / 0.2e1 + t81 / 0.2e1 + t138 / 0.2e1 + t84 / 0.2e1 + t82 / 0.2e1 + t139 / 0.2e1 + Ifges(4,4) * t142 + t265 * t88 + t266 * t87 + (Ifges(4,2) + t264) * t143 + t219) * t186 + (t20 * t243 - Ifges(4,1) * t142 - Ifges(4,4) * t143 + (mrSges(4,3) + t245) * t77 + (-t1 * t214 - t211 * t3) * mrSges(6,3) + (-t211 * t7 - t214 * t8) * mrSges(5,3) + (t86 * t244 + t126 * t246 + (t211 * t23 - t214 * t25) * mrSges(6,3) - t330 * mrSges(5,3) + t341 * t322 + t339 * t320 + t343 * t318 - t335 * t214 / 0.2e1) * qJD(4) + t338 * t329 + t340 * t328 + t342 * t323 + (qJD(4) * t334 + t351) * t314 + t350 * t311) * t187 + (-t175 / 0.2e1 - pkin(6) * t196 + t259 + (-0.2e1 * t325 - 0.3e1 / 0.2e1 * t305 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t216) * qJD(1) + (qJD(1) * (mrSges(4,1) * t186 + mrSges(4,2) * t187) + 0.2e1 * t310 + t145) * pkin(2)) * t273; t203 * t218 - t287 * t151 - m(5) * (t126 * t151 + t42 * t62 + t43 * t63) + t346 * t105 - m(4) * (-t149 * t151 + t150 * t152) + t345 * t96 + (m(4) * (t212 * t76 - t215 * t77) + (t142 * t215 - t143 * t212) * mrSges(4,3) + ((-m(4) * t149 + m(5) * t126 + t287) * t212 + (m(4) * t150 + m(5) * t330 + t214 * t106 - t211 * t108 + t166) * t215) * qJD(3)) * pkin(2) + ((t260 - t176 / 0.2e1 - t207 / 0.2e1 + qJD(1) * t324 + (t195 - t284) * pkin(6)) * t216 + (t259 + t175 / 0.2e1 + (t325 + t305 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t216) * qJD(1) + (t196 + t283) * pkin(6) + (-t145 - t310) * pkin(2)) * t213) * qJD(1) + t348 * t107 + (-t203 * t30 + t250) * t211 + ((-t108 * t203 + t249) * t214 + (-t106 * t203 + t248) * t211) * qJD(4) + t217 + t150 * t294 + t32 * t277 - t63 * t106 - t62 * t108 - t152 * t166 + t184 * t29 + t185 * t31 + t194 * t21 + t347 * (-pkin(3) - t309) + (t1 * t184 + t185 * t3 + t194 * t20 + t23 * t348 + t25 * t346 + t345 * t86) * m(6); (-t287 + t294) * t150 + ((-pkin(8) * t108 + t249) * t214 + (pkin(4) * t96 - pkin(8) * t106 + t248) * t211) * qJD(4) - m(5) * (t126 * t150 + t42 * t67 + t43 * t68) + (-pkin(8) * t30 + t250) * t211 + t333 * t107 + t332 * t105 + t217 + t32 * t308 + pkin(8) * t218 - t98 * t96 - t68 * t106 - t67 * t108 - t149 * t166 + t197 * t29 + t198 * t31 + t205 * t21 - t347 * pkin(3) + (t1 * t197 + t198 * t3 + t20 * t205 + (t267 - t98) * t86 + t332 * t25 + t333 * t23) * m(6); (-t163 * t96 + t29) * pkin(4) + (t162 * t23 + t163 * t25) * mrSges(6,3) + (t162 * t42 + t163 * t43) * mrSges(5,3) + t138 + t139 + t84 + t83 + t82 + t81 + t219 + (-(-t23 + t24) * t25 + (-t163 * t86 + t1) * pkin(4)) * m(6) - t24 * t105 - t42 * t106 + t25 * t107 + t43 * t108 - t86 * (mrSges(6,1) * t163 + mrSges(6,2) * t162) - t126 * (mrSges(5,1) * t163 + mrSges(5,2) * t162) + (t363 * t162 - t357) * t320 + t335 * t319 + (t162 * t353 - t163 * t352) * t318 + (-t361 * t163 + t334 + t359) * t322; -t162 * t105 + t163 * t107 + 0.2e1 * (t20 / 0.2e1 + t25 * t322 + t23 * t319) * m(6) + t21;];
tauc = t2(:);
