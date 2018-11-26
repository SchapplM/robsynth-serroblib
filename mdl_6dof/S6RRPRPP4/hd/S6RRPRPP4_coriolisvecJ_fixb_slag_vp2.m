% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:58:31
% EndTime: 2018-11-23 16:58:41
% DurationCPUTime: 9.08s
% Computational Cost: add. (5885->598), mult. (13974->772), div. (0->0), fcn. (8365->6), ass. (0->272)
t356 = -qJD(1) / 0.2e1;
t353 = Ifges(6,1) + Ifges(7,1);
t352 = Ifges(6,4) - Ifges(7,5);
t351 = Ifges(7,4) + Ifges(6,5);
t219 = sin(qJ(2));
t220 = cos(qJ(4));
t290 = cos(pkin(9));
t255 = t290 * t220;
t249 = t219 * t255;
t282 = qJD(1) * t219;
t217 = sin(pkin(9));
t218 = sin(qJ(4));
t288 = t217 * t218;
t130 = qJD(1) * t249 - t282 * t288;
t276 = qJD(4) * t218;
t148 = -qJD(4) * t255 + t217 * t276;
t285 = t130 - t148;
t230 = -t217 * t220 - t218 * t290;
t227 = t219 * t230;
t131 = qJD(1) * t227;
t149 = t230 * qJD(4);
t284 = -t131 - t149;
t221 = cos(qJ(2));
t281 = qJD(1) * t221;
t208 = pkin(7) * t281;
t175 = pkin(3) * t281 + t208;
t216 = qJD(2) * qJ(3);
t155 = t216 + t175;
t222 = -pkin(2) - pkin(8);
t256 = -qJ(3) * t219 - pkin(1);
t163 = t221 * t222 + t256;
t134 = t163 * qJD(1);
t206 = pkin(7) * t282;
t174 = -pkin(3) * t282 - t206;
t340 = qJD(3) - t174;
t136 = qJD(2) * t222 + t340;
t81 = -t134 * t218 + t220 * t136;
t82 = t134 * t220 + t136 * t218;
t237 = t218 * t81 - t220 * t82;
t303 = Ifges(5,4) * t218;
t241 = Ifges(5,2) * t220 + t303;
t302 = Ifges(5,4) * t220;
t243 = Ifges(5,1) * t218 + t302;
t246 = mrSges(5,1) * t220 - mrSges(5,2) * t218;
t312 = -t220 / 0.2e1;
t313 = -t218 / 0.2e1;
t279 = qJD(2) * t220;
t169 = -t218 * t281 + t279;
t316 = -t169 / 0.2e1;
t168 = -qJD(2) * t218 - t220 * t281;
t317 = -t168 / 0.2e1;
t199 = qJD(4) + t282;
t295 = t169 * Ifges(5,4);
t97 = t168 * Ifges(5,2) + t199 * Ifges(5,6) + t295;
t164 = Ifges(5,4) * t168;
t98 = t169 * Ifges(5,1) + t199 * Ifges(5,5) + t164;
t367 = t237 * mrSges(5,3) + t155 * t246 + t241 * t317 + t243 * t316 + t312 * t97 + t313 * t98;
t366 = qJD(3) + t206;
t365 = Ifges(6,6) / 0.2e1;
t364 = -Ifges(7,6) / 0.2e1;
t280 = qJD(2) * t219;
t262 = t218 * t280;
t271 = qJD(2) * qJD(4);
t274 = qJD(4) * t221;
t122 = -t218 * t271 + (-t220 * t274 + t262) * qJD(1);
t261 = t218 * t274;
t229 = t219 * t279 + t261;
t123 = qJD(1) * t229 - t220 * t271;
t74 = t122 * t217 - t123 * t290;
t332 = -t74 / 0.2e1;
t75 = t122 * t290 + t217 * t123;
t330 = t75 / 0.2e1;
t327 = pkin(3) + pkin(7);
t363 = -mrSges(3,1) + mrSges(4,2);
t264 = -pkin(4) * t220 - pkin(3);
t275 = qJD(4) * t220;
t342 = pkin(4) * t275 - t264 * t282 + t366;
t181 = -pkin(2) * t221 + t256;
t156 = t181 * qJD(1);
t185 = -t208 - t216;
t297 = Ifges(5,6) * t220;
t301 = Ifges(5,5) * t218;
t240 = t297 + t301;
t315 = -t199 / 0.2e1;
t354 = qJD(2) / 0.2e1;
t355 = -qJD(2) / 0.2e1;
t362 = t185 * mrSges(4,1) - t156 * mrSges(4,2) - Ifges(4,5) * t355 - Ifges(3,6) * t354 - t240 * t315 - t367 + ((Ifges(3,2) + Ifges(4,3)) * t221 + (Ifges(3,4) + Ifges(4,6)) * t219) * t356;
t63 = qJ(5) * t168 + t82;
t293 = t217 * t63;
t62 = -qJ(5) * t169 + t81;
t54 = pkin(4) * t199 + t62;
t16 = t290 * t54 - t293;
t11 = -t199 * pkin(5) + qJD(6) - t16;
t110 = -pkin(4) * t168 + qJD(5) + t155;
t106 = -t290 * t168 + t169 * t217;
t231 = t217 * t168 + t169 * t290;
t34 = pkin(5) * t106 - qJ(6) * t231 + t110;
t348 = -t106 * t352 + t199 * t351 + t231 * t353;
t361 = t110 * mrSges(6,2) + t11 * mrSges(7,2) - t16 * mrSges(6,3) - t34 * mrSges(7,3) + t348 / 0.2e1;
t166 = -t255 + t288;
t60 = t290 * t63;
t17 = t217 * t54 + t60;
t272 = qJD(1) * qJD(2);
t257 = t221 * t272;
t258 = t219 * t272;
t198 = pkin(2) * t258;
t239 = pkin(8) * t219 - qJ(3) * t221;
t277 = qJD(3) * t219;
t226 = qJD(2) * t239 - t277;
t117 = qJD(1) * t226 + t198;
t278 = qJD(2) * t221;
t177 = t327 * t278;
t161 = qJD(1) * t177;
t36 = -qJD(4) * t82 - t117 * t218 + t220 * t161;
t13 = pkin(4) * t257 - qJ(5) * t122 - qJD(5) * t169 + t36;
t35 = t220 * t117 - t134 * t276 + t136 * t275 + t218 * t161;
t15 = qJ(5) * t123 + qJD(5) * t168 + t35;
t3 = t13 * t290 - t217 * t15;
t4 = t217 * t13 + t290 * t15;
t360 = t16 * t284 + t166 * t3 - t17 * t285 + t230 * t4;
t1 = qJ(6) * t257 + qJD(6) * t199 + t4;
t12 = qJ(6) * t199 + t17;
t2 = -pkin(5) * t257 - t3;
t359 = t1 * t230 - t11 * t284 - t12 * t285 - t166 * t2;
t358 = -t75 * Ifges(7,5) / 0.2e1 + t257 * t364 + Ifges(7,3) * t332;
t357 = Ifges(6,4) * t330 + Ifges(6,2) * t332 + t257 * t365;
t324 = -t106 / 0.2e1;
t323 = t106 / 0.2e1;
t314 = t199 / 0.2e1;
t320 = t231 / 0.2e1;
t350 = Ifges(6,6) - Ifges(7,6);
t349 = t257 * t351 - t352 * t74 + t353 * t75;
t347 = pkin(5) * t285 + qJ(6) * t284 + qJD(6) * t166 + t342;
t87 = -mrSges(6,2) * t199 - mrSges(6,3) * t106;
t88 = -mrSges(7,2) * t106 + mrSges(7,3) * t199;
t307 = t88 + t87;
t89 = mrSges(6,1) * t199 - mrSges(6,3) * t231;
t90 = -mrSges(7,1) * t199 + mrSges(7,2) * t231;
t306 = t90 - t89;
t189 = t327 * t219;
t112 = t220 * t163 + t218 * t189;
t341 = t218 * t35 + t220 * t36;
t339 = t36 * mrSges(5,1) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t35 * mrSges(5,2) - t4 * mrSges(6,2) + t1 * mrSges(7,3) + Ifges(5,5) * t122 + Ifges(5,6) * t123;
t180 = -qJD(2) * pkin(2) + t366;
t205 = Ifges(3,4) * t281;
t252 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t266 = t365 + t364;
t268 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t299 = Ifges(4,6) * t221;
t338 = -t266 * t106 + t252 * t199 + t268 * t231 + t12 * mrSges(7,3) + t16 * mrSges(6,1) + t180 * mrSges(4,1) + t81 * mrSges(5,1) + Ifges(3,1) * t282 / 0.2e1 + Ifges(3,5) * t354 + t205 / 0.2e1 + Ifges(4,4) * t355 + (-Ifges(4,2) * t219 - t299) * t356 + Ifges(6,6) * t324 + Ifges(7,6) * t323 + t169 * Ifges(5,5) + t168 * Ifges(5,6) - t11 * mrSges(7,1) - t156 * mrSges(4,3) - t17 * mrSges(6,2) - t82 * mrSges(5,2) + t351 * t320 + (Ifges(6,3) + Ifges(7,2) + Ifges(5,3)) * t314;
t43 = Ifges(7,5) * t231 + t199 * Ifges(7,6) + t106 * Ifges(7,3);
t46 = Ifges(6,4) * t231 - t106 * Ifges(6,2) + t199 * Ifges(6,6);
t335 = -mrSges(6,1) * t110 - mrSges(7,1) * t34 + mrSges(7,2) * t12 + mrSges(6,3) * t17 + t46 / 0.2e1 - t43 / 0.2e1;
t331 = t74 / 0.2e1;
t328 = t97 / 0.2e1;
t326 = pkin(1) * mrSges(3,1);
t325 = pkin(1) * mrSges(3,2);
t321 = -t231 / 0.2e1;
t311 = pkin(4) * t169;
t310 = pkin(4) * t217;
t211 = pkin(2) * t280;
t127 = t211 + t226;
t160 = t220 * t177;
t253 = qJ(5) * t221 - t163;
t30 = pkin(4) * t278 + t160 + t253 * t275 + (-qJ(5) * t280 - qJD(4) * t189 + qJD(5) * t221 - t127) * t218;
t273 = t220 * qJD(5);
t49 = t220 * t127 - t163 * t276 + t218 * t177 + t189 * t275;
t37 = qJ(5) * t229 - t221 * t273 + t49;
t8 = t217 * t30 + t290 * t37;
t56 = -mrSges(7,2) * t74 + mrSges(7,3) * t257;
t57 = -mrSges(6,2) * t257 - mrSges(6,3) * t74;
t309 = t56 + t57;
t58 = mrSges(6,1) * t257 - mrSges(6,3) * t75;
t59 = -mrSges(7,1) * t257 + t75 * mrSges(7,2);
t308 = -t58 + t59;
t207 = pkin(2) * t282;
t143 = qJD(1) * t239 + t207;
t100 = -t143 * t218 + t220 * t175;
t79 = (-qJ(5) * t218 * t219 + pkin(4) * t221) * qJD(1) + t100;
t101 = t220 * t143 + t218 * t175;
t84 = qJ(5) * t220 * t282 + t101;
t32 = t217 * t79 + t290 * t84;
t171 = t220 * t189;
t92 = pkin(4) * t219 + t218 * t253 + t171;
t287 = t220 * t221;
t99 = -qJ(5) * t287 + t112;
t41 = t217 * t92 + t290 * t99;
t305 = mrSges(5,3) * t168;
t304 = mrSges(5,3) * t169;
t300 = Ifges(5,5) * t220;
t298 = Ifges(5,6) * t218;
t289 = qJD(2) * mrSges(3,2);
t204 = t218 * pkin(4) + qJ(3);
t286 = qJ(5) - t222;
t115 = -mrSges(5,1) * t168 + mrSges(5,2) * t169;
t187 = -mrSges(4,1) * t281 - qJD(2) * mrSges(4,3);
t283 = -t187 + t115;
t190 = t327 * t221;
t270 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t269 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t267 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t265 = m(4) * pkin(7) + mrSges(4,1);
t150 = pkin(4) * t287 + t190;
t263 = t290 * pkin(4);
t26 = t74 * mrSges(6,1) + t75 * mrSges(6,2);
t25 = t74 * mrSges(7,1) - t75 * mrSges(7,3);
t254 = t286 * t220;
t176 = t327 * t280;
t251 = m(4) * t180 + (mrSges(4,1) + mrSges(3,3)) * t282 + t363 * qJD(2);
t250 = m(4) * t185 - mrSges(3,3) * t281 + t187 + t289;
t245 = mrSges(5,1) * t218 + mrSges(5,2) * t220;
t244 = Ifges(5,1) * t220 - t303;
t242 = -Ifges(5,2) * t218 + t302;
t103 = mrSges(5,1) * t257 - mrSges(5,3) * t122;
t104 = -mrSges(5,2) * t257 + mrSges(5,3) * t123;
t236 = t220 * t103 + t218 * t104;
t125 = -mrSges(5,2) * t199 + t305;
t126 = mrSges(5,1) * t199 - t304;
t235 = t220 * t125 - t218 * t126;
t7 = -t217 * t37 + t290 * t30;
t31 = -t217 * t84 + t290 * t79;
t40 = -t217 * t99 + t290 * t92;
t232 = -qJ(3) * t278 - t277;
t215 = qJD(2) * qJD(3);
t140 = -qJD(1) * t176 + t215;
t228 = t276 * t286 - t273;
t91 = -pkin(4) * t123 + t140;
t116 = -pkin(4) * t261 + (-pkin(7) + t264) * t280;
t203 = -t263 - pkin(5);
t201 = qJ(6) + t310;
t197 = Ifges(7,2) * t257;
t196 = Ifges(5,3) * t257;
t195 = Ifges(6,3) * t257;
t179 = t286 * t218;
t178 = pkin(7) * t258 - t215;
t172 = (mrSges(4,2) * t221 - mrSges(4,3) * t219) * qJD(1);
t145 = t211 + t232;
t139 = t230 * t221;
t138 = t166 * t221;
t137 = -qJD(4) * t254 - t218 * qJD(5);
t132 = qJD(1) * t232 + t198;
t114 = -t179 * t290 - t217 * t254;
t113 = -t179 * t217 + t254 * t290;
t111 = -t163 * t218 + t171;
t102 = -pkin(5) * t230 + qJ(6) * t166 + t204;
t94 = -qJD(2) * t227 + t166 * t274;
t93 = qJD(2) * t249 - t149 * t221 - t217 * t262;
t86 = t137 * t290 + t217 * t228;
t85 = t137 * t217 - t228 * t290;
t78 = -mrSges(5,1) * t123 + mrSges(5,2) * t122;
t76 = -pkin(5) * t138 - qJ(6) * t139 + t150;
t73 = Ifges(7,4) * t75;
t72 = Ifges(6,5) * t75;
t71 = Ifges(6,6) * t74;
t70 = Ifges(7,6) * t74;
t65 = t122 * Ifges(5,1) + t123 * Ifges(5,4) + Ifges(5,5) * t257;
t64 = t122 * Ifges(5,4) + t123 * Ifges(5,2) + Ifges(5,6) * t257;
t53 = mrSges(6,1) * t106 + mrSges(6,2) * t231;
t52 = mrSges(7,1) * t106 - mrSges(7,3) * t231;
t50 = -qJD(4) * t112 - t127 * t218 + t160;
t42 = pkin(5) * t231 + qJ(6) * t106 + t311;
t39 = -t219 * pkin(5) - t40;
t38 = qJ(6) * t219 + t41;
t28 = -pkin(5) * t281 - t31;
t27 = qJ(6) * t281 + t32;
t20 = t290 * t62 - t293;
t19 = t217 * t62 + t60;
t18 = -pkin(5) * t93 - qJ(6) * t94 - qJD(6) * t139 + t116;
t9 = pkin(5) * t74 - qJ(6) * t75 - qJD(6) * t231 + t91;
t6 = -pkin(5) * t278 - t7;
t5 = qJ(6) * t278 + qJD(6) * t219 + t8;
t10 = [m(7) * (t1 * t38 + t11 * t6 + t12 * t5 + t18 * t34 + t2 * t39 + t76 * t9) + m(6) * (t110 * t116 + t150 * t91 + t16 * t7 + t17 * t8 + t3 * t40 + t4 * t41) + (t1 * t138 + t139 * t2) * mrSges(7,2) + m(4) * (t132 * t181 + t145 * t156) + (t65 * t313 + t64 * t312 - t122 * t243 / 0.2e1 - t123 * t241 / 0.2e1 + t140 * t246 + t132 * mrSges(4,2) - t265 * t178 + (t218 * t36 - t220 * t35) * mrSges(5,3) + ((t298 - t300) * t314 + t244 * t316 + t242 * t317 - t155 * t245 + t98 * t312 + t218 * t328 + (t218 * t82 + t220 * t81) * mrSges(5,3)) * qJD(4) + (t251 * pkin(7) + ((-t301 / 0.2e1 - t297 / 0.2e1 - t267) * t221 - 0.2e1 * t325 - t181 * mrSges(4,3) + t268 * t139 + t266 * t138 + (-0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + t265 * pkin(7) + t252) * t219) * qJD(1) + t270 * qJD(2) + t338) * qJD(2)) * t221 + (t250 * pkin(7) + (-t181 * mrSges(4,2) + t219 * t267 - 0.2e1 * t326) * qJD(1) + t269 * qJD(2) + t362) * t280 + (Ifges(6,4) * t324 + Ifges(7,5) * t323 + t314 * t351 + t353 * t320 + t361) * t94 + (Ifges(6,4) * t139 + Ifges(6,2) * t138) * t332 + (Ifges(7,5) * t139 - Ifges(7,3) * t138) * t331 + (Ifges(6,2) * t324 - Ifges(7,3) * t323 + t314 * t350 + t320 * t352 + t335) * t93 + (t138 * t4 - t139 * t3) * mrSges(6,3) + m(5) * (t111 * t36 + t112 * t35 + t140 * t190 - t155 * t176 + t49 * t82 + t50 * t81) + (t197 / 0.2e1 + t195 / 0.2e1 + t196 / 0.2e1 + t73 / 0.2e1 + t72 / 0.2e1 - t71 / 0.2e1 + t70 / 0.2e1 + t268 * t75 - t266 * t74 - t132 * mrSges(4,3) + t339) * t219 + t190 * t78 - t176 * t115 + t145 * t172 + t150 * t26 + t91 * (-mrSges(6,1) * t138 + mrSges(6,2) * t139) + t9 * (-mrSges(7,1) * t138 - mrSges(7,3) * t139) + t49 * t125 + t50 * t126 + t112 * t104 + t116 * t53 + t111 * t103 + t5 * t88 + t7 * t89 + t6 * t90 + t8 * t87 + t76 * t25 + (t138 * t352 + t139 * t353) * t330 + t349 * t139 / 0.2e1 + t138 * t357 + t138 * t358 + t18 * t52 + t38 * t56 + t41 * t57 + t40 * t58 + t39 * t59; (-Ifges(7,5) * t166 - Ifges(7,3) * t230) * t331 + (-Ifges(6,4) * t166 + Ifges(6,2) * t230) * t332 + t9 * (-mrSges(7,1) * t230 + mrSges(7,3) * t166) + t91 * (-mrSges(6,1) * t230 - mrSges(6,2) * t166) + (-t166 * t353 + t230 * t352) * t330 + t230 * t357 + t230 * t358 + (-Ifges(6,4) * t131 + Ifges(7,5) * t149 + Ifges(6,2) * t130 - Ifges(7,3) * t148) * t323 + (Ifges(6,4) * t149 - Ifges(7,5) * t131 + Ifges(6,2) * t148 - Ifges(7,3) * t130) * t324 + (qJD(4) * t240 + t130 * t350 - t131 * t351) * t315 + (t130 * t352 - t131 * t353) * t321 + t348 * (t149 / 0.2e1 + t131 / 0.2e1) + t64 * t313 + (((-t250 + t289) * pkin(7) + (t326 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t219) * qJD(1) + (-qJ(3) * mrSges(4,1) + t269) * qJD(2) - t362) * t219 + (((-m(4) * pkin(2) + t363) * qJD(2) - t251) * pkin(7) + (t325 - t299 / 0.2e1) * qJD(1) + (-Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t282 + (t300 / 0.2e1 - t298 / 0.2e1 - pkin(2) * mrSges(4,1) + t266 * t230 - t268 * t166 + t270) * qJD(2) - t205 / 0.2e1 - t338) * t221) * qJD(1) + t367 * qJD(4) + (-m(4) * t156 - t172) * (-qJ(3) * t281 + t207) + (t46 - t43) * (t148 / 0.2e1 - t130 / 0.2e1) + (t236 + m(5) * t341 + (-m(5) * t237 + t235) * qJD(4)) * t222 + t342 * t53 + (-t113 * t3 + t114 * t4 + t204 * t91 + (-t32 + t86) * t17 + (-t31 - t85) * t16 + t342 * t110) * m(6) + (t140 * qJ(3) - t100 * t81 - t101 * t82 + t155 * t340) * m(5) + m(4) * (-qJ(3) * t178 - qJD(3) * t185) + t220 * t65 / 0.2e1 + t204 * t26 - t178 * mrSges(4,3) - t174 * t115 - t101 * t125 - t100 * t126 + t102 * t25 - t27 * t88 - t31 * t89 - t28 * t90 - t32 * t87 + qJ(3) * t78 + (t148 * t350 + t149 * t351) * t314 + (t148 * t352 + t149 * t353) * t320 + t123 * t242 / 0.2e1 + t122 * t244 / 0.2e1 + t140 * t245 + t347 * t52 + (t1 * t114 + t102 * t9 + t113 * t2 + t347 * t34 + (-t27 + t86) * t12 + (-t28 + t85) * t11) * m(7) - t349 * t166 / 0.2e1 - t341 * mrSges(5,3) + t283 * qJD(3) + (mrSges(7,1) * t285 + mrSges(7,3) * t284) * t34 + (mrSges(6,1) * t285 - mrSges(6,2) * t284) * t110 + t306 * t85 + t307 * t86 + t359 * mrSges(7,2) + t360 * mrSges(6,3) + t308 * t113 + t309 * t114; -t309 * t230 + t308 * t166 + t235 * qJD(4) + (t172 + t235) * t282 + (t265 * t281 - t283 - t52 - t53) * qJD(2) - m(4) * (-qJD(2) * t185 - t156 * t282) + t236 + t285 * t307 + t284 * t306 + (-qJD(2) * t34 - t359) * m(7) + (-qJD(2) * t110 - t360) * m(6) + (-qJD(2) * t155 - t199 * t237 + t341) * m(5); ((t217 * t4 + t290 * t3) * pkin(4) - t110 * t311 + t16 * t19 - t17 * t20) * m(6) + (t1 * t201 - t11 * t19 + t2 * t203 - t34 * t42 + (qJD(6) - t20) * t12) * m(7) + (-Ifges(5,2) * t169 + t164 + t98) * t317 + t197 + t195 + t196 + (t304 + t126) * t82 + (-Ifges(6,4) * t323 - Ifges(7,5) * t324 - t351 * t315 - t353 * t321 + t361) * t106 + t169 * t328 + t57 * t310 + (Ifges(5,5) * t168 - Ifges(5,6) * t169) * t315 + (Ifges(5,1) * t168 - t295) * t316 + t73 + t72 - t71 + t70 + (t305 - t125) * t81 + t201 * t56 + t203 * t59 - t155 * (mrSges(5,1) * t169 + mrSges(5,2) * t168) + qJD(6) * t88 + (-Ifges(6,2) * t323 + Ifges(7,3) * t324 - t350 * t315 - t352 * t321 + t335) * t231 + t58 * t263 + t339 - t306 * t19 - t307 * t20 - t53 * t311 - t42 * t52; -t306 * t231 + t307 * t106 + t25 + t26 + (t106 * t12 - t11 * t231 + t9) * m(7) + (t106 * t17 + t16 * t231 + t91) * m(6); t231 * t52 - t199 * t88 + 0.2e1 * (t2 / 0.2e1 + t34 * t320 + t12 * t315) * m(7) + t59;];
tauc  = t10(:);
