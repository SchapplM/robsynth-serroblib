% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:17:23
% EndTime: 2018-11-23 16:17:34
% DurationCPUTime: 11.34s
% Computational Cost: add. (15998->587), mult. (42951->812), div. (0->0), fcn. (34277->10), ass. (0->277)
t252 = sin(pkin(11));
t254 = cos(pkin(11));
t256 = sin(qJ(6));
t259 = cos(qJ(6));
t231 = t252 * t259 + t254 * t256;
t220 = t231 * qJD(6);
t253 = sin(pkin(10));
t258 = sin(qJ(3));
t255 = cos(pkin(10));
t260 = cos(qJ(3));
t300 = t255 * t260;
t274 = t253 * t258 - t300;
t270 = t274 * qJD(1);
t342 = cos(qJ(4));
t209 = t342 * t270;
t232 = t253 * t260 + t255 * t258;
t218 = t232 * qJD(1);
t257 = sin(qJ(4));
t190 = t218 * t257 + t209;
t390 = t231 * t190;
t406 = t390 + t220;
t275 = t252 * t256 - t254 * t259;
t219 = t275 * qJD(6);
t389 = t275 * t190;
t401 = t389 + t219;
t251 = qJD(3) + qJD(4);
t268 = t257 * t270;
t266 = t218 * t342 - t268;
t172 = t251 * t252 + t254 * t266;
t335 = pkin(7) + qJ(2);
t238 = t335 * t253;
t233 = qJD(1) * t238;
t240 = t335 * t255;
t234 = qJD(1) * t240;
t198 = -t258 * t233 + t260 * t234;
t175 = -pkin(8) * t270 + t198;
t168 = t342 * t175;
t197 = -t260 * t233 - t234 * t258;
t174 = -pkin(8) * t218 + t197;
t169 = qJD(3) * pkin(3) + t174;
t108 = t257 * t169 + t168;
t102 = t251 * qJ(5) + t108;
t290 = -t255 * pkin(2) - pkin(1);
t207 = pkin(3) * t274 + t290;
t199 = qJD(1) * t207 + qJD(2);
t109 = t190 * pkin(4) - qJ(5) * t266 + t199;
t68 = -t102 * t252 + t254 * t109;
t35 = pkin(5) * t190 - pkin(9) * t172 + t68;
t286 = t254 * t251 - t252 * t266;
t69 = t254 * t102 + t252 * t109;
t43 = pkin(9) * t286 + t69;
t12 = -t256 * t43 + t259 * t35;
t13 = t256 * t35 + t259 * t43;
t221 = t274 * qJD(3);
t212 = qJD(1) * t221;
t222 = t232 * qJD(3);
t213 = qJD(1) * t222;
t146 = -t342 * t212 + (-qJD(4) * t218 - t213) * t257 - qJD(4) * t209;
t289 = qJD(4) * t342;
t147 = -qJD(4) * t268 - t257 * t212 + t213 * t342 + t218 * t289;
t243 = qJD(2) * t300;
t295 = qJD(3) * t260;
t296 = qJD(2) * t253;
t165 = -t233 * t295 + qJD(1) * t243 + (-qJD(1) * t296 - qJD(3) * t234) * t258;
t158 = -pkin(8) * t213 + t165;
t271 = t232 * qJD(2);
t166 = -qJD(1) * t271 - qJD(3) * t198;
t265 = pkin(8) * t212 + t166;
t294 = qJD(4) * t257;
t55 = t342 * t158 + t169 * t289 - t175 * t294 + t257 * t265;
t51 = qJD(5) * t251 + t55;
t341 = pkin(3) * t213;
t67 = pkin(4) * t147 - qJ(5) * t146 - qJD(5) * t266 + t341;
t24 = t252 * t67 + t254 * t51;
t310 = t146 * t252;
t14 = -pkin(9) * t310 + t24;
t23 = -t252 * t51 + t254 * t67;
t309 = t146 * t254;
t8 = pkin(5) * t147 - pkin(9) * t309 + t23;
t2 = qJD(6) * t12 + t14 * t259 + t256 * t8;
t3 = -qJD(6) * t13 - t14 * t256 + t259 * t8;
t314 = t24 * t254;
t56 = qJD(4) * t108 + t257 * t158 - t342 * t265;
t32 = pkin(5) * t310 + t56;
t329 = Ifges(6,4) * t254;
t330 = Ifges(6,4) * t252;
t343 = t254 / 0.2e1;
t186 = qJD(6) + t190;
t351 = t186 / 0.2e1;
t355 = t147 / 0.2e1;
t113 = t172 * t259 + t256 * t286;
t356 = t113 / 0.2e1;
t387 = -t172 * t256 + t259 * t286;
t358 = t387 / 0.2e1;
t110 = Ifges(7,4) * t387;
t60 = Ifges(7,1) * t113 + Ifges(7,5) * t186 + t110;
t360 = t60 / 0.2e1;
t328 = Ifges(7,4) * t113;
t59 = Ifges(7,2) * t387 + Ifges(7,6) * t186 + t328;
t362 = t59 / 0.2e1;
t42 = -qJD(6) * t113 - t146 * t231;
t364 = t42 / 0.2e1;
t41 = qJD(6) * t387 - t146 * t275;
t365 = t41 / 0.2e1;
t366 = Ifges(7,1) * t365 + Ifges(7,4) * t364 + Ifges(7,5) * t355;
t367 = Ifges(7,4) * t365 + Ifges(7,2) * t364 + Ifges(7,6) * t355;
t326 = Ifges(6,2) * t252;
t280 = -t326 + t329;
t53 = t147 * Ifges(6,6) + t146 * t280;
t281 = Ifges(6,1) * t254 - t330;
t54 = t147 * Ifges(6,5) + t146 * t281;
t167 = t257 * t175;
t107 = t169 * t342 - t167;
t101 = -t251 * pkin(4) + qJD(5) - t107;
t88 = -pkin(5) * t286 + t101;
t405 = mrSges(6,3) * t314 + (-mrSges(6,1) * t254 + mrSges(6,2) * t252 - mrSges(5,1)) * t56 + (-Ifges(7,5) * t219 - Ifges(7,6) * t220) * t351 + (-Ifges(7,1) * t219 - Ifges(7,4) * t220) * t356 - t389 * t60 / 0.2e1 - t390 * t59 / 0.2e1 + (-Ifges(7,4) * t219 - Ifges(7,2) * t220) * t358 - t219 * t360 + t53 * t343 + (Ifges(7,4) * t231 - Ifges(7,2) * t275) * t364 + (Ifges(7,1) * t231 - Ifges(7,4) * t275) * t365 + t32 * (mrSges(7,1) * t275 + mrSges(7,2) * t231) - t275 * t367 + (Ifges(6,5) * t252 + Ifges(7,5) * t231 + Ifges(6,6) * t254 - Ifges(7,6) * t275) * t355 - t55 * mrSges(5,2) - t220 * t362 + t231 * t366 + Ifges(5,5) * t146 - Ifges(5,6) * t147 + (Ifges(6,1) * t252 + t329) * t309 / 0.2e1 - (Ifges(6,2) * t254 + t330) * t310 / 0.2e1 + (t406 * mrSges(7,1) - mrSges(7,2) * t401) * t88 + t252 * t54 / 0.2e1 + (t12 * t401 - t13 * t406 - t2 * t275 - t3 * t231) * mrSges(7,3);
t404 = t266 / 0.2e1;
t391 = t190 * t252;
t403 = pkin(5) * t391;
t402 = pkin(9) * t391;
t400 = (m(3) * qJ(2) + mrSges(3,3)) * (-t253 ^ 2 - t255 ^ 2);
t398 = Ifges(7,1) * t389 + Ifges(7,4) * t390;
t397 = Ifges(7,4) * t389 + Ifges(7,2) * t390;
t396 = Ifges(7,5) * t389 + Ifges(7,6) * t390;
t185 = Ifges(5,4) * t190;
t394 = -t185 / 0.2e1;
t393 = t251 * Ifges(5,5) / 0.2e1;
t392 = Ifges(5,1) * t404;
t282 = mrSges(6,1) * t252 + mrSges(6,2) * t254;
t272 = t101 * t282;
t344 = -t252 / 0.2e1;
t386 = t394 + t393 + t392;
t96 = t172 * Ifges(6,4) + Ifges(6,2) * t286 + Ifges(6,6) * t190;
t97 = t172 * Ifges(6,1) + Ifges(6,4) * t286 + Ifges(6,5) * t190;
t388 = t272 + t286 * t280 / 0.2e1 + t172 * t281 / 0.2e1 + t199 * mrSges(5,2) + (-t252 * t69 - t254 * t68) * mrSges(6,3) + t386 + t393 + t96 * t344 + t97 * t343;
t148 = pkin(4) * t266 + qJ(5) * t190;
t347 = -t266 / 0.2e1;
t384 = pkin(5) * t266;
t383 = t12 * mrSges(7,1);
t382 = t13 * mrSges(7,2);
t331 = Ifges(5,4) * t266;
t338 = pkin(3) * t257;
t245 = qJ(5) + t338;
t223 = (-pkin(9) - t245) * t252;
t248 = t254 * pkin(9);
t303 = t245 * t254;
t224 = t248 + t303;
t194 = t223 * t256 + t224 * t259;
t285 = pkin(3) * t289;
t242 = t285 + qJD(5);
t305 = t190 * t254;
t116 = t174 * t342 - t167;
t340 = pkin(3) * t218;
t122 = t148 + t340;
t75 = -t116 * t252 + t254 * t122;
t37 = pkin(9) * t305 + t384 + t75;
t76 = t254 * t116 + t252 * t122;
t61 = t76 + t402;
t381 = -qJD(6) * t194 - t231 * t242 + t256 * t61 - t259 * t37;
t193 = t223 * t259 - t224 * t256;
t380 = qJD(6) * t193 - t242 * t275 - t256 * t37 - t259 * t61;
t237 = (-pkin(9) - qJ(5)) * t252;
t311 = qJ(5) * t254;
t239 = t248 + t311;
t202 = t237 * t256 + t239 * t259;
t81 = -t107 * t252 + t254 * t148;
t44 = t190 * t248 + t384 + t81;
t82 = t254 * t107 + t252 * t148;
t64 = t82 + t402;
t379 = -qJD(5) * t231 - qJD(6) * t202 + t256 * t64 - t259 * t44;
t200 = t237 * t259 - t239 * t256;
t378 = -qJD(5) * t275 + qJD(6) * t200 - t256 * t44 - t259 * t64;
t298 = -mrSges(5,1) * t251 - mrSges(6,1) * t286 + mrSges(6,2) * t172 + mrSges(5,3) * t266;
t201 = -t260 * t238 - t240 * t258;
t183 = -pkin(8) * t232 + t201;
t203 = -t258 * t238 + t260 * t240;
t184 = -pkin(8) * t274 + t203;
t373 = t342 * t183 - t257 * t184;
t127 = -mrSges(6,2) * t190 + mrSges(6,3) * t286;
t128 = mrSges(6,1) * t190 - mrSges(6,3) * t172;
t371 = t254 * t127 - t252 * t128;
t370 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t41 + Ifges(7,6) * t42;
t74 = -mrSges(7,1) * t387 + mrSges(7,2) * t113;
t368 = m(6) * t101 + m(7) * t88 + t298 + t74;
t359 = -t387 / 0.2e1;
t357 = -t113 / 0.2e1;
t352 = -t186 / 0.2e1;
t349 = -t190 / 0.2e1;
t348 = t190 / 0.2e1;
t346 = -t221 / 0.2e1;
t345 = -t222 / 0.2e1;
t339 = pkin(3) * t222;
t337 = t254 * pkin(5);
t180 = -t238 * t295 + t243 + (-qJD(3) * t240 - t296) * t258;
t162 = -pkin(8) * t222 + t180;
t181 = -qJD(3) * t203 - t271;
t163 = pkin(8) * t221 + t181;
t72 = qJD(4) * t373 + t342 * t162 + t257 * t163;
t273 = -t257 * t232 - t274 * t342;
t156 = qJD(4) * t273 - t221 * t342 - t257 * t222;
t196 = t232 * t342 - t257 * t274;
t157 = qJD(4) * t196 - t257 * t221 + t222 * t342;
t79 = pkin(4) * t157 - qJ(5) * t156 - qJD(5) * t196 + t339;
t28 = t252 * t79 + t254 * t72;
t333 = mrSges(4,3) * t218;
t332 = Ifges(4,4) * t218;
t327 = Ifges(6,5) * t254;
t325 = Ifges(6,6) * t252;
t324 = t107 * mrSges(5,3);
t323 = t108 * mrSges(5,3);
t322 = t387 * Ifges(7,6);
t321 = t113 * Ifges(7,5);
t320 = t373 * t56;
t319 = t286 * Ifges(6,6);
t318 = t172 * Ifges(6,5);
t317 = t186 * Ifges(7,3);
t315 = t23 * t252;
t312 = t251 * Ifges(5,6);
t308 = t156 * t252;
t304 = t196 * t252;
t129 = -pkin(4) * t273 - qJ(5) * t196 + t207;
t135 = t257 * t183 + t184 * t342;
t84 = t252 * t129 + t254 * t135;
t85 = mrSges(6,1) * t310 + mrSges(6,2) * t309;
t293 = t342 * pkin(3);
t291 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t16 = -t42 * mrSges(7,1) + t41 * mrSges(7,2);
t27 = -t252 * t72 + t254 * t79;
t288 = t213 * mrSges(4,1) - t212 * mrSges(4,2);
t287 = t147 * mrSges(5,1) + t146 * mrSges(5,2);
t83 = t254 * t129 - t135 * t252;
t115 = t174 * t257 + t168;
t247 = -t293 - pkin(4);
t279 = -t325 + t327;
t278 = t23 * t254 + t24 * t252;
t277 = t314 - t315;
t276 = -t252 * t68 + t254 * t69;
t57 = -pkin(5) * t273 - t196 * t248 + t83;
t71 = -pkin(9) * t304 + t84;
t25 = -t256 * t71 + t259 * t57;
t26 = t256 * t57 + t259 * t71;
t269 = mrSges(4,3) * t270;
t73 = qJD(4) * t135 + t257 * t162 - t342 * t163;
t138 = -t190 * Ifges(5,2) + t312 + t331;
t58 = t317 + t321 + t322;
t95 = t190 * Ifges(6,3) + t318 + t319;
t262 = t383 + t199 * mrSges(5,1) + t68 * mrSges(6,1) - t138 / 0.2e1 + t58 / 0.2e1 + t95 / 0.2e1 + t322 / 0.2e1 + t321 / 0.2e1 - t382 + t319 / 0.2e1 + t318 / 0.2e1 + t317 / 0.2e1 - t312 / 0.2e1 - t69 * mrSges(6,2);
t246 = -pkin(4) - t337;
t236 = t247 - t337;
t235 = qJD(1) * t290 + qJD(2);
t214 = Ifges(4,4) * t270;
t206 = qJD(3) * mrSges(4,1) - t333;
t205 = -qJD(3) * mrSges(4,2) - t269;
t188 = t218 * Ifges(4,1) + Ifges(4,5) * qJD(3) - t214;
t187 = -Ifges(4,2) * t270 + Ifges(4,6) * qJD(3) + t332;
t176 = -mrSges(5,2) * t251 - mrSges(5,3) * t190;
t151 = t275 * t196;
t150 = t231 * t196;
t149 = mrSges(5,1) * t190 + mrSges(5,2) * t266;
t141 = Ifges(7,3) * t147;
t98 = pkin(5) * t304 - t373;
t92 = t115 - t403;
t91 = mrSges(7,1) * t186 - mrSges(7,3) * t113;
t90 = -mrSges(7,2) * t186 + mrSges(7,3) * t387;
t89 = t108 - t403;
t87 = mrSges(6,1) * t147 - mrSges(6,3) * t309;
t86 = -mrSges(6,2) * t147 - mrSges(6,3) * t310;
t63 = -t156 * t231 + t196 * t219;
t62 = -t156 * t275 - t196 * t220;
t36 = pkin(5) * t308 + t73;
t31 = -mrSges(7,2) * t147 + mrSges(7,3) * t42;
t30 = mrSges(7,1) * t147 - mrSges(7,3) * t41;
t22 = -pkin(9) * t308 + t28;
t15 = pkin(5) * t157 - t156 * t248 + t27;
t5 = -qJD(6) * t26 + t15 * t259 - t22 * t256;
t4 = qJD(6) * t25 + t15 * t256 + t22 * t259;
t1 = [(t279 * t348 + t388 + t392) * t156 + (-Ifges(7,5) * t151 - Ifges(7,6) * t150) * t355 + (-t12 * t62 + t13 * t63 - t150 * t2 + t151 * t3) * mrSges(7,3) + (-Ifges(7,4) * t151 - Ifges(7,2) * t150) * t364 + (-Ifges(7,1) * t151 - Ifges(7,4) * t150) * t365 + t32 * (mrSges(7,1) * t150 - mrSges(7,2) * t151) + (-t232 * t212 + t218 * t346) * Ifges(4,1) + (-t165 * t274 - t166 * t232 + t197 * t221 - t198 * t222 + t201 * t212 - t203 * t213) * mrSges(4,3) + (t212 * t274 - t232 * t213 + t218 * t345) * Ifges(4,4) + m(5) * (-t107 * t73 + t108 * t72 - t320 + t135 * t55 + (t199 * t222 + t207 * t213) * pkin(3)) + m(6) * (t101 * t73 + t23 * t83 + t24 * t84 + t27 * t68 + t28 * t69 - t320) - (t141 / 0.2e1 + mrSges(5,1) * t341 - t24 * mrSges(6,2) + t23 * mrSges(6,1) + (Ifges(6,3) + Ifges(5,2) + Ifges(7,3) / 0.2e1) * t147 + t370) * t273 + (mrSges(5,2) * t341 - mrSges(6,3) * t278 + t279 * t355 + t282 * t56 + t343 * t54 + t344 * t53) * t196 + qJD(3) * (-Ifges(4,5) * t221 - Ifges(4,6) * t222) / 0.2e1 + t235 * (mrSges(4,1) * t222 - mrSges(4,2) * t221) + t98 * t16 + t188 * t346 + (Ifges(7,5) * t62 + Ifges(7,6) * t63) * t351 + (Ifges(7,1) * t62 + Ifges(7,4) * t63) * t356 + (Ifges(7,4) * t62 + Ifges(7,2) * t63) * t358 + t187 * t345 + t36 * t74 + ((Ifges(5,1) + Ifges(6,1) * t254 ^ 2 / 0.2e1 + (-t329 + t326 / 0.2e1) * t252) * t196 - t273 * t279) * t146 + (t146 * t273 - t196 * t147 + t156 * t349 + t157 * t347) * Ifges(5,4) + t88 * (-mrSges(7,1) * t63 + mrSges(7,2) * t62) + t4 * t90 + t5 * t91 + t84 * t86 + t83 * t87 + t207 * t287 + (-t107 * t156 - t108 * t157 - t135 * t147 - t146 * t373 + t196 * t56 + t273 * t55) * mrSges(5,3) - t373 * t85 + t25 * t30 + t26 * t31 + m(7) * (t12 * t5 + t13 * t4 + t2 * t26 + t25 * t3 + t32 * t98 + t36 * t88) + m(4) * (t165 * t203 + t166 * t201 + t180 * t198 + t181 * t197) + t62 * t360 + t63 * t362 - t151 * t366 - t150 * t367 + (-t274 * (-Ifges(4,4) * t221 - Ifges(4,2) * t222) / 0.2e1 - 0.2e1 * t400 * qJD(2)) * qJD(1) + t290 * t288 + t28 * t127 + t27 * t128 + (t190 * t291 + t262) * t157 + t298 * t73 + t72 * t176 + t180 * t205 + t181 * t206 + t149 * t339 + t274 * Ifges(4,2) * t213; t205 * t270 - m(4) * (-t197 * t218 - t198 * t270) + m(6) * t278 + t287 + t288 + t218 * t206 - t275 * t30 + t231 * t31 + t252 * t86 + t254 * t87 + m(5) * t341 - t406 * t91 - t401 * t90 + (-t12 * t406 - t13 * t401 + t2 * t231 - t275 * t3) * m(7) - (-m(5) * t108 - m(6) * t276 - t176 - t371) * t190 + (m(5) * t107 - t368) * t266 + t400 * qJD(1) ^ 2; -t235 * (t218 * mrSges(4,1) - mrSges(4,2) * t270) - qJD(3) * (-Ifges(4,5) * t270 - Ifges(4,6) * t218) / 0.2e1 + t86 * t303 + (-Ifges(4,2) * t218 + t188 - t214) * t270 / 0.2e1 + t368 * pkin(3) * t294 + (-t269 - t205) * t197 + (-t146 * t293 - t147 * t338) * mrSges(5,3) + (t285 - t116) * t176 + (t333 + t206) * t198 + t138 * t404 + t266 * t382 + ((-t342 * t56 + t257 * t55 + (-t107 * t257 + t108 * t342) * qJD(4)) * pkin(3) + t107 * t115 - t108 * t116 - t199 * t340) * m(5) + (Ifges(7,6) * t266 + t397) * t359 + (Ifges(7,5) * t266 + t398) * t357 + (Ifges(7,3) * t266 + t396) * t352 + (-t101 * t115 + t242 * t276 + t245 * t277 + t247 * t56 - t68 * t75 - t69 * t76) * m(6) - t172 * (Ifges(6,5) * t266 - t190 * t281) / 0.2e1 - t190 * t324 + t190 * t272 + (Ifges(6,3) * t266 - t190 * t279) * t349 - t286 * (Ifges(6,6) * t266 - t190 * t280) / 0.2e1 - t199 * (mrSges(5,1) * t266 - mrSges(5,2) * t190) - t251 * (-Ifges(5,5) * t190 - Ifges(5,6) * t266) / 0.2e1 + t190 * t386 + (-Ifges(5,1) * t190 - t331 + t58 + t95) * t347 - t92 * t74 + (-Ifges(5,2) * t266 - t185) * t348 - t69 * (-mrSges(6,2) * t266 + mrSges(6,3) * t391) - t96 * t391 / 0.2e1 - t218 * (-Ifges(4,1) * t270 - t332) / 0.2e1 - t68 * (mrSges(6,1) * t266 + mrSges(6,3) * t305) + t371 * t242 - t298 * t115 - t76 * t127 - t75 * t128 + t380 * t90 + t381 * t91 + (t12 * t381 + t13 * t380 + t193 * t3 + t194 * t2 + t236 * t32 - t88 * t92) * m(7) - t266 * t383 - t165 * mrSges(4,2) + t166 * mrSges(4,1) + t266 * t323 + t97 * t305 / 0.2e1 + t193 * t30 + t194 * t31 - mrSges(6,3) * t315 - Ifges(4,5) * t212 - Ifges(4,6) * t213 - t252 * t245 * t87 + t405 + t218 * t187 / 0.2e1 + t236 * t16 + t247 * t85 - t149 * t340; (-t324 + t394 + (t327 / 0.2e1 - t325 / 0.2e1) * t190 + (Ifges(5,1) / 0.2e1 - t291) * t266 + t388) * t190 + (-t262 + t331 / 0.2e1 + t323) * t266 + t86 * t311 + (-t23 * mrSges(6,3) - qJ(5) * t87 - qJD(5) * t128) * t252 + t397 * t359 + t398 * t357 + t396 * t352 - t89 * t74 - pkin(4) * t85 + (qJD(5) * t254 - t82) * t127 + (-pkin(4) * t56 + qJ(5) * t277 + qJD(5) * t276 - t101 * t108 - t68 * t81 - t69 * t82) * m(6) - t81 * t128 + t378 * t90 + t379 * t91 + (t12 * t379 + t13 * t378 + t2 * t202 + t200 * t3 + t246 * t32 - t88 * t89) * m(7) - t298 * t108 - t107 * t176 + t200 * t30 + t202 * t31 + t405 + t246 * t16; t113 * t91 - t387 * t90 - t286 * t127 + t172 * t128 + t16 + t85 + (t113 * t12 - t13 * t387 + t32) * m(7) + (t172 * t68 - t286 * t69 + t56) * m(6); t141 - t88 * (mrSges(7,1) * t113 + mrSges(7,2) * t387) + (Ifges(7,1) * t387 - t328) * t357 + t59 * t356 + (Ifges(7,5) * t387 - Ifges(7,6) * t113) * t352 - t12 * t90 + t13 * t91 + (t113 * t13 + t12 * t387) * mrSges(7,3) + (-Ifges(7,2) * t113 + t110 + t60) * t359 + t370;];
tauc  = t1(:);
