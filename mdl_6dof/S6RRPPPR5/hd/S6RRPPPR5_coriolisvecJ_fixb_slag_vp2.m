% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:21
% EndTime: 2019-03-09 08:21:38
% DurationCPUTime: 8.26s
% Computational Cost: add. (3819->591), mult. (9700->777), div. (0->0), fcn. (5915->6), ass. (0->276)
t364 = Ifges(4,1) + Ifges(6,3);
t363 = Ifges(6,1) + Ifges(5,3);
t366 = Ifges(3,5) / 0.2e1;
t365 = -qJD(2) / 0.2e1;
t362 = Ifges(6,4) - Ifges(5,5);
t361 = Ifges(4,5) - Ifges(6,6);
t360 = Ifges(6,5) - Ifges(5,6);
t204 = sin(pkin(9));
t205 = cos(pkin(9));
t206 = sin(qJ(6));
t208 = cos(qJ(6));
t152 = t204 * t208 + t205 * t206;
t207 = sin(qJ(2));
t271 = -pkin(7) * t204 - pkin(3);
t229 = (-qJ(5) + t271) * t207;
t209 = cos(qJ(2));
t294 = t205 * t209;
t277 = pkin(8) * t294;
t215 = t229 + t277;
t284 = qJD(1) * t209;
t268 = t205 * t284;
t237 = pkin(2) * t207 - qJ(3) * t209;
t156 = t237 * qJD(1);
t296 = t205 * t156;
t261 = pkin(4) * t268 - t296;
t43 = qJD(1) * t215 + t261;
t332 = pkin(4) + pkin(8);
t275 = t204 * t332;
t260 = t209 * t275;
t214 = (-pkin(7) * t205 + pkin(5)) * t207 - t260;
t140 = t204 * t156;
t285 = qJD(1) * t207;
t187 = qJ(4) * t285;
t290 = t140 + t187;
t44 = qJD(1) * t214 + t290;
t196 = t204 * qJ(3);
t163 = t204 * pkin(4) + t196;
t147 = pkin(8) * t204 + t163;
t197 = t205 * qJ(3);
t164 = t205 * pkin(4) + t197;
t148 = pkin(8) * t205 + t164;
t75 = -t147 * t206 + t148 * t208;
t359 = qJD(3) * t152 + qJD(6) * t75 - t206 * t44 - t208 * t43;
t153 = t204 * t206 - t205 * t208;
t76 = t147 * t208 + t148 * t206;
t358 = -qJD(3) * t153 - qJD(6) * t76 + t206 * t43 - t208 * t44;
t306 = Ifges(5,6) * t205;
t308 = Ifges(6,5) * t205;
t357 = t204 * t363 - t306 + t308;
t309 = Ifges(6,5) * t204;
t311 = Ifges(4,4) * t204;
t356 = t205 * t364 + t309 - t311;
t270 = t204 * t285;
t279 = t205 * qJD(2);
t150 = t270 - t279;
t355 = -t150 / 0.2e1;
t354 = t150 / 0.2e1;
t269 = t205 * t285;
t280 = t204 * qJD(2);
t151 = t269 + t280;
t353 = -t151 / 0.2e1;
t352 = t151 / 0.2e1;
t278 = qJD(1) * qJD(2);
t351 = -t278 / 0.2e1;
t350 = t278 / 0.2e1;
t267 = Ifges(3,6) * t365;
t348 = qJD(2) * t366;
t315 = pkin(3) + qJ(5);
t347 = t150 * t315;
t346 = t205 * t315;
t262 = -t150 * t206 + t208 * t151;
t295 = t205 * t207;
t345 = pkin(4) * t295 + t209 * qJ(5);
t281 = qJD(4) * t209;
t283 = qJD(2) * t207;
t344 = qJ(4) * t283 - t281;
t255 = Ifges(6,6) / 0.2e1 - Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t256 = Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1 + Ifges(4,6) / 0.2e1;
t343 = -t204 * t256 - t205 * t255;
t137 = t152 * qJD(6);
t185 = qJD(6) - t284;
t198 = t209 * qJD(5);
t264 = t209 * t278;
t257 = t205 * t264;
t135 = qJD(2) * t237 - t207 * qJD(3);
t119 = t135 * qJD(1);
t191 = pkin(7) * t285;
t157 = (qJD(3) - t191) * qJD(2);
t60 = t205 * t119 - t204 * t157;
t228 = pkin(4) * t257 + qJD(1) * t198 - t60;
t17 = (-t207 * t315 + t277) * t278 + t228;
t265 = t207 * t278;
t61 = t204 * t119 + t205 * t157;
t274 = qJ(4) * t265 + t61;
t18 = (-t281 + (pkin(5) * t207 - t260) * qJD(2)) * qJD(1) + t274;
t193 = pkin(3) * t284;
t161 = -pkin(2) * t209 - t207 * qJ(3) - pkin(1);
t144 = t161 * qJD(1);
t192 = pkin(7) * t284;
t171 = qJD(2) * qJ(3) + t192;
t85 = t144 * t205 - t204 * t171;
t62 = qJD(4) + t193 - t85;
t220 = qJ(5) * t284 + t62;
t20 = t151 * t332 + t220;
t314 = -pkin(5) - qJ(4);
t266 = t314 * t209;
t86 = t204 * t144 + t205 * t171;
t21 = qJD(1) * t266 - t150 * t332 + qJD(5) + t86;
t5 = -t20 * t206 + t208 * t21;
t1 = qJD(6) * t5 + t17 * t208 + t18 * t206;
t6 = t20 * t208 + t206 * t21;
t2 = -qJD(6) * t6 - t17 * t206 + t18 * t208;
t222 = t152 * t209;
t218 = qJD(2) * t222;
t35 = qJD(1) * t218 + qJD(6) * t262;
t221 = t153 * t209;
t217 = qJD(2) * t221;
t82 = t150 * t208 + t151 * t206;
t36 = -qJD(1) * t217 - qJD(6) * t82;
t342 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t35 - Ifges(7,6) * t36;
t203 = t209 * pkin(3);
t282 = qJD(2) * t209;
t287 = pkin(7) * t282 + t280 * t203;
t341 = -(qJD(4) * t205 - qJD(5) * t204) * t207 + t287;
t312 = Ifges(3,4) * t207;
t34 = pkin(4) * t151 + t220;
t64 = qJ(4) * t284 - t86;
t39 = -t150 * pkin(4) + qJD(5) - t64;
t340 = -t256 * t150 - t255 * t151 + t39 * mrSges(6,1) + t5 * mrSges(7,1) + t62 * mrSges(5,2) + t85 * mrSges(4,1) + t185 * Ifges(7,3) + t82 * Ifges(7,5) + t262 * Ifges(7,6) + t267 - (t209 * Ifges(3,2) + t312) * qJD(1) / 0.2e1 + Ifges(5,5) * t354 + Ifges(4,5) * t352 - t34 * mrSges(6,3) - t6 * mrSges(7,2) - t64 * mrSges(5,3) - t86 * mrSges(4,2) + (Ifges(4,6) + Ifges(6,4)) * t355 + (Ifges(5,4) + Ifges(6,6)) * t353 - (Ifges(5,1) + Ifges(4,3) + Ifges(6,2)) * t284 / 0.2e1;
t190 = Ifges(3,4) * t284;
t307 = Ifges(5,6) * t204;
t239 = -Ifges(5,2) * t205 + t307;
t310 = Ifges(4,4) * t205;
t241 = -Ifges(4,2) * t204 + t310;
t244 = -mrSges(5,2) * t204 - mrSges(5,3) * t205;
t245 = mrSges(6,1) * t205 - mrSges(6,3) * t204;
t249 = qJD(2) * pkin(2) - qJD(3) - t191;
t320 = t205 / 0.2e1;
t321 = -t205 / 0.2e1;
t322 = t204 / 0.2e1;
t323 = -t204 / 0.2e1;
t219 = qJ(4) * t151 + t249;
t38 = t219 - t347;
t58 = pkin(3) * t150 - t219;
t339 = (t204 * t64 + t205 * t62) * mrSges(5,1) + (t204 * t39 - t205 * t34) * mrSges(6,2) - (t204 * t86 + t205 * t85) * mrSges(4,3) + (-m(4) * t249 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t150 + mrSges(4,2) * t151 + mrSges(3,3) * t285) * pkin(7) - t249 * (mrSges(4,1) * t204 + mrSges(4,2) * t205) + t38 * t245 + t58 * t244 + Ifges(3,1) * t285 / 0.2e1 + t190 / 0.2e1 + t348 + t241 * t355 + t239 * t353 + t357 * t354 + t356 * t352 + (-Ifges(5,4) * t321 - Ifges(4,6) * t323 - t320 * t361 + t322 * t362) * t284 + (Ifges(4,4) * t323 - Ifges(5,2) * t321 + t320 * t364 + t322 * t360) * t151 + (-Ifges(4,2) * t323 + Ifges(5,6) * t321 + t363 * t322 + (-Ifges(4,4) + Ifges(6,5)) * t320) * t150;
t338 = t35 / 0.2e1;
t337 = t36 / 0.2e1;
t336 = -t262 / 0.2e1;
t335 = t262 / 0.2e1;
t334 = -t82 / 0.2e1;
t333 = t82 / 0.2e1;
t331 = pkin(1) * mrSges(3,1);
t330 = pkin(1) * mrSges(3,2);
t127 = t153 * t207;
t329 = -t127 / 0.2e1;
t128 = t152 * t207;
t328 = t128 / 0.2e1;
t325 = -t185 / 0.2e1;
t324 = t185 / 0.2e1;
t319 = Ifges(7,4) * t82;
t22 = -mrSges(7,1) * t262 + mrSges(7,2) * t82;
t90 = mrSges(6,1) * t151 - mrSges(6,3) * t150;
t313 = -t22 + t90;
t302 = qJ(5) * t204;
t301 = qJD(2) * mrSges(3,2);
t299 = t204 * t207;
t298 = t204 * t209;
t297 = t205 * t135;
t112 = t150 * mrSges(5,1) + mrSges(5,3) * t284;
t114 = -mrSges(6,1) * t284 + t150 * mrSges(6,2);
t293 = t114 - t112;
t117 = qJD(1) * t221;
t138 = t153 * qJD(6);
t292 = -t117 + t138;
t118 = qJD(1) * t222;
t291 = t118 - t137;
t289 = -qJ(4) * t257 - t151 * qJD(4);
t258 = t204 * t264;
t130 = mrSges(6,1) * t265 + mrSges(6,2) * t258;
t122 = mrSges(4,1) * t258 + mrSges(4,2) * t257;
t288 = t204 * t193 + t192;
t134 = mrSges(5,1) * t257 + mrSges(5,2) * t265;
t286 = pkin(3) * t299 + t207 * pkin(7);
t110 = pkin(7) * t294 + t204 * t161;
t276 = pkin(7) * t283;
t111 = -t151 * mrSges(6,2) + mrSges(6,3) * t284;
t113 = t151 * mrSges(5,1) - mrSges(5,2) * t284;
t116 = -mrSges(4,1) * t284 - t151 * mrSges(4,3);
t273 = t111 + t113 - t116;
t115 = mrSges(4,2) * t284 - t150 * mrSges(4,3);
t272 = t115 + t293;
t9 = -t36 * mrSges(7,1) + t35 * mrSges(7,2);
t263 = qJ(4) * t204 + pkin(2);
t182 = pkin(7) * t298;
t109 = t161 * t205 - t182;
t253 = t207 * t271;
t252 = t204 * t315 + pkin(7);
t103 = qJ(4) * t209 - t110;
t248 = t209 * pkin(4) * t279 + t198 - t297;
t104 = -t109 + t203;
t123 = t204 * t135;
t247 = t123 + t344;
t246 = t150 * qJD(5) + t289;
t50 = -mrSges(7,2) * t185 + mrSges(7,3) * t262;
t51 = mrSges(7,1) * t185 - mrSges(7,3) * t82;
t233 = -t206 * t51 + t208 * t50;
t232 = -t206 * t50 - t208 * t51;
t52 = t182 + t203 + (pkin(8) * t207 - t161) * t205 + t345;
t53 = -t207 * t275 + t110 + t266;
t13 = t206 * t53 + t208 * t52;
t12 = -t206 * t52 + t208 * t53;
t231 = qJ(4) * t205 - t302;
t106 = -pkin(7) * t269 + t140;
t95 = -t205 * t276 + t123;
t225 = t205 * t314 + t302;
t224 = -pkin(4) * t298 - pkin(7) * t295;
t223 = t231 * t209;
t216 = t209 * t225;
t24 = -t265 * t315 + t228;
t27 = t252 * t264 + t246;
t48 = -pkin(3) * t265 - t60;
t54 = (pkin(3) * t204 + pkin(7)) * t264 + t289;
t213 = (Ifges(5,4) * t207 + t209 * t239) * t351 - t27 * mrSges(6,1) - t54 * mrSges(5,3) - t24 * mrSges(6,2) + t48 * mrSges(5,1) - t60 * mrSges(4,3) + (t207 * t361 + t209 * t356) * t350;
t28 = (-pkin(4) * t280 - qJD(4)) * t284 + t274;
t45 = qJD(1) * t281 - t274;
t212 = (Ifges(4,6) * t207 + t209 * t241) * t351 + t27 * mrSges(6,3) - t54 * mrSges(5,2) + t45 * mrSges(5,1) + t28 * mrSges(6,2) - t61 * mrSges(4,3) + (-t207 * t362 + t209 * t357) * t350;
t180 = Ifges(7,3) * t265;
t172 = mrSges(3,3) * t284 - t301;
t159 = qJD(4) * t204 + qJD(5) * t205;
t158 = -pkin(3) * t205 - t263;
t139 = t263 + t346;
t133 = (mrSges(5,1) * t298 - mrSges(5,3) * t207) * t278;
t132 = (-mrSges(6,2) * t294 - mrSges(6,3) * t207) * t278;
t131 = (mrSges(4,1) * t207 - mrSges(4,3) * t294) * t278;
t129 = (-mrSges(4,2) * t207 - mrSges(4,3) * t298) * t278;
t125 = -qJ(4) * t295 + t286;
t124 = t204 * t314 - pkin(2) - t346;
t121 = t244 * t264;
t120 = t245 * t264;
t108 = -qJ(4) * t268 + t288;
t105 = pkin(7) * t270 + t296;
t102 = t207 * t231 - t286;
t94 = t204 * t276 + t297;
t93 = qJD(1) * t253 - t296;
t92 = -mrSges(5,2) * t150 - mrSges(5,3) * t151;
t89 = -t106 - t187;
t88 = (-qJ(4) * t282 - qJD(4) * t207) * t205 + t287;
t87 = qJD(1) * t223 - t288;
t84 = -pkin(4) * t299 - t103;
t79 = t207 * t225 + t286;
t78 = qJD(2) * t253 - t297;
t77 = Ifges(7,4) * t262;
t65 = t104 + t345;
t63 = qJD(1) * t224 + t290;
t59 = -t95 - t344;
t57 = qJD(1) * t216 + t288;
t56 = -t138 * t207 + t218;
t55 = -t137 * t207 - t217;
t49 = qJD(1) * t229 + t261;
t47 = qJD(2) * t223 - t341;
t46 = qJD(2) * t224 + t247;
t42 = qJD(2) * t229 + t248;
t37 = qJD(2) * t216 + t341;
t30 = qJD(2) * t214 + t247;
t29 = qJD(2) * t215 + t248;
t26 = -mrSges(7,2) * t265 + mrSges(7,3) * t36;
t25 = mrSges(7,1) * t265 - mrSges(7,3) * t35;
t23 = t151 * t314 - t249 + t347;
t19 = (-pkin(5) * t205 + t252) * t264 + t246;
t16 = Ifges(7,1) * t82 + Ifges(7,5) * t185 + t77;
t15 = Ifges(7,2) * t262 + Ifges(7,6) * t185 + t319;
t8 = t35 * Ifges(7,1) + t36 * Ifges(7,4) + Ifges(7,5) * t265;
t7 = t35 * Ifges(7,4) + t36 * Ifges(7,2) + Ifges(7,6) * t265;
t4 = -qJD(6) * t13 - t206 * t29 + t208 * t30;
t3 = qJD(6) * t12 + t206 * t30 + t208 * t29;
t10 = [m(7) * (t1 * t13 + t12 * t2 + t19 * t79 + t23 * t37 + t3 * t6 + t4 * t5) + m(6) * (-t102 * t27 + t24 * t65 + t28 * t84 + t34 * t42 + t38 * t47 + t39 * t46) + m(5) * (t103 * t45 + t104 * t48 + t125 * t54 + t58 * t88 + t59 * t64 + t62 * t78) + (((Ifges(7,5) * t328 + Ifges(7,6) * t329 - 0.3e1 / 0.2e1 * t312 - 0.2e1 * t331 + t343 * t207) * qJD(1) - pkin(7) * t172 + t267 + t340) * t207 + ((0.3e1 / 0.2e1 * Ifges(3,4) * t209 - 0.2e1 * t330 + (0.3e1 / 0.2e1 * Ifges(6,6) + 0.3e1 / 0.2e1 * Ifges(5,4) - 0.3e1 / 0.2e1 * Ifges(4,5)) * t294 + (0.3e1 / 0.2e1 * Ifges(6,4) - 0.3e1 / 0.2e1 * Ifges(5,5) + 0.3e1 / 0.2e1 * Ifges(4,6)) * t298 + (m(4) * pkin(7) ^ 2 - Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + (pkin(7) * mrSges(4,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t205) * t205 + (pkin(7) * mrSges(4,1) + (Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t204 + (-Ifges(4,4) + t360) * t205) * t204) * t207) * qJD(1) + t348 + t339) * t209) * qJD(2) + m(4) * (t60 * t109 + t61 * t110 + t85 * t94 + t86 * t95) + t37 * t22 + t12 * t25 + t13 * t26 + (-t1 * t127 - t128 * t2 - t5 * t56 + t55 * t6) * mrSges(7,3) + t19 * (mrSges(7,1) * t127 + mrSges(7,2) * t128) + (Ifges(7,4) * t128 - Ifges(7,2) * t127) * t337 + (Ifges(7,1) * t128 - Ifges(7,4) * t127) * t338 + (-t180 / 0.2e1 + t45 * mrSges(5,3) - t28 * mrSges(6,1) + t61 * mrSges(4,2) + t24 * mrSges(6,3) - t48 * mrSges(5,2) - t60 * mrSges(4,1) + t342) * t209 + (pkin(7) * t122 + t204 * t212 + t205 * t213) * t207 + t3 * t50 + t4 * t51 + t55 * t15 / 0.2e1 + t56 * t16 / 0.2e1 + t23 * (-mrSges(7,1) * t55 + mrSges(7,2) * t56) + t79 * t9 + t47 * t90 + t88 * t92 + t42 * t111 + t59 * t112 + t78 * t113 + t46 * t114 + t95 * t115 + t94 * t116 + t102 * t120 + t125 * t121 + t110 * t129 + t84 * t130 + t109 * t131 + t65 * t132 + t103 * t133 + t104 * t134 + (Ifges(7,5) * t56 + Ifges(7,6) * t55) * t324 + t8 * t328 + t7 * t329 + (Ifges(7,1) * t56 + Ifges(7,4) * t55) * t333 + (Ifges(7,4) * t56 + Ifges(7,2) * t55) * t335; t358 * t51 + (t1 * t76 + t124 * t19 + t2 * t75 + t359 * t6 + t358 * t5 + (-t159 - t57) * t23) * m(7) + t359 * t50 + (t117 / 0.2e1 - t138 / 0.2e1) * t15 + (Ifges(7,5) * t118 - Ifges(7,6) * t117) * t325 + (Ifges(7,5) * t137 - Ifges(7,6) * t138) * t324 + (Ifges(7,1) * t137 - Ifges(7,4) * t138) * t333 + (Ifges(7,4) * t137 - Ifges(7,2) * t138) * t335 + m(4) * (-t196 * t60 + t197 * t61 + (-t204 * t85 + t205 * t86) * qJD(3)) + (Ifges(7,1) * t118 - Ifges(7,4) * t117) * t334 + (Ifges(7,4) * t118 - Ifges(7,2) * t117) * t336 - m(4) * (t85 * t105 + t86 * t106) + (-t118 / 0.2e1 + t137 / 0.2e1) * t16 + (t158 * t54 + (-qJ(3) * t45 - qJD(3) * t64) * t205 + (qJ(3) * t48 + qJD(3) * t62 - qJD(4) * t58) * t204 - t108 * t58 - t62 * t93 - t64 * t89) * m(5) + (-t139 * t27 + t163 * t24 + t164 * t28 + (t204 * t34 + t205 * t39) * qJD(3) - t34 * t49 - t39 * t63 + (t159 - t87) * t38) * m(6) + (((t331 + t312 / 0.2e1) * qJD(1) + (t172 + t301) * pkin(7) + t267 + (-Ifges(6,4) * t205 + Ifges(6,6) * t204) * t365 - t340) * t207 + (((Ifges(4,2) * t205 + t311) * t323 + (-Ifges(5,2) * t204 - t306) * t321 + t366 + (-m(4) * pkin(2) - mrSges(4,1) * t205 + mrSges(4,2) * t204 - mrSges(3,1)) * pkin(7) + (-t205 * t363 - t307 + t309) * t322 + (t204 * t364 - t308 + t310) * t320) * qJD(2) + (t209 * t343 + t330) * qJD(1) + (-Ifges(3,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t285 - t190 / 0.2e1 - t339) * t209 + (Ifges(7,5) * t153 + Ifges(7,6) * t152 + (-Ifges(5,5) + Ifges(4,6)) * t205 + (-Ifges(5,4) + Ifges(4,5)) * t204) * t283 / 0.2e1) * qJD(1) - t57 * t22 + t75 * t25 + t76 * t26 - t87 * t90 - t108 * t92 - t49 * t111 - t89 * t112 - t93 * t113 - t63 * t114 - t106 * t115 - t105 * t116 - pkin(2) * t122 + t124 * t9 + t139 * t120 + t152 * t7 / 0.2e1 + t153 * t8 / 0.2e1 + t19 * (-mrSges(7,1) * t152 + mrSges(7,2) * t153) + (Ifges(7,4) * t153 + Ifges(7,2) * t152) * t337 + (Ifges(7,1) * t153 + Ifges(7,4) * t152) * t338 + t158 * t121 + t163 * t132 + t164 * t130 + ((t129 - t133) * qJ(3) + t272 * qJD(3) - t212) * t205 + (-qJD(4) * t92 + (-t131 + t134) * qJ(3) + t273 * qJD(3) + t213) * t204 + (mrSges(7,1) * t292 - mrSges(7,2) * t291) * t23 + (t1 * t152 - t153 * t2 + t291 * t5 - t292 * t6) * mrSges(7,3) + t313 * t159; -t262 * t50 + t82 * t51 - t273 * t151 + t272 * t150 + (m(4) * pkin(7) + (-mrSges(6,1) - mrSges(5,3)) * t205 + (-mrSges(5,2) + mrSges(6,3)) * t204) * t264 - m(4) * (-t150 * t86 - t151 * t85) + t9 + t122 + (-t262 * t6 + t5 * t82 + t19) * m(7) + (t150 * t39 - t151 * t34 + t27) * m(6) + (-t150 * t64 - t151 * t62 + t54) * m(5); -t206 * t25 + t208 * t26 + t232 * qJD(6) + (t92 - t313) * t151 + (-mrSges(6,3) * t283 + (-mrSges(6,2) * t279 - t232 + t293) * t209) * qJD(1) + t134 + (t1 * t208 + t23 * t151 - t2 * t206 + t185 * (-t206 * t6 - t208 * t5)) * m(7) + (-t38 * t151 + t284 * t39 + t24) * m(6) + (t58 * t151 - t284 * t64 + t48) * m(5); t206 * t26 + t208 * t25 + t313 * t150 + t233 * qJD(6) + (-t111 - t233) * t284 + t130 + (t1 * t206 - t23 * t150 + t2 * t208 + t185 * (-t206 * t5 + t208 * t6)) * m(7) + (t38 * t150 - t284 * t34 + t28) * m(6); t180 - t23 * (mrSges(7,1) * t82 + mrSges(7,2) * t262) + (Ifges(7,1) * t262 - t319) * t334 + t15 * t333 + (Ifges(7,5) * t262 - Ifges(7,6) * t82) * t325 - t5 * t50 + t6 * t51 + (t262 * t5 + t6 * t82) * mrSges(7,3) + (-Ifges(7,2) * t82 + t16 + t77) * t336 - t342;];
tauc  = t10(:);
