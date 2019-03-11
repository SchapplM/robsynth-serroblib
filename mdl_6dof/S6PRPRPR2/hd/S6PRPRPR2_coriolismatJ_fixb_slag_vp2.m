% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:46
% EndTime: 2019-03-08 19:29:55
% DurationCPUTime: 5.17s
% Computational Cost: add. (10508->461), mult. (25986->683), div. (0->0), fcn. (28506->12), ass. (0->248)
t394 = -m(7) / 0.2e1;
t266 = cos(qJ(4));
t261 = sin(pkin(11));
t249 = pkin(2) * t261 + pkin(8);
t260 = sin(pkin(12));
t304 = pkin(5) * t260 + t249;
t218 = t304 * t266;
t262 = cos(pkin(12));
t263 = sin(qJ(6));
t369 = cos(qJ(6));
t282 = t263 * t260 - t369 * t262;
t211 = t282 * t266;
t358 = t211 * mrSges(7,2);
t283 = t260 * t369 + t263 * t262;
t210 = t283 * t266;
t359 = t210 * mrSges(7,1);
t317 = t359 / 0.2e1 - t358 / 0.2e1;
t405 = t218 * t394 - t317;
t256 = t260 ^ 2;
t257 = t262 ^ 2;
t314 = t256 + t257;
t302 = t314 * mrSges(6,3);
t404 = mrSges(5,2) - t302;
t366 = pkin(9) + qJ(5);
t237 = t366 * t260;
t239 = t366 * t262;
t174 = -t237 * t369 - t263 * t239;
t175 = -t263 * t237 + t239 * t369;
t348 = sin(pkin(6));
t349 = cos(pkin(11));
t297 = t349 * t348;
t265 = sin(qJ(2));
t303 = t265 * t348;
t370 = cos(qJ(2));
t206 = t261 * t303 - t297 * t370;
t299 = t370 * t348;
t207 = t261 * t299 + t265 * t297;
t324 = t262 * t266;
t120 = -t206 * t324 + t207 * t260;
t342 = t120 * t262;
t327 = t260 * t266;
t119 = t206 * t327 + t207 * t262;
t343 = t119 * t260;
t291 = t342 - t343;
t251 = -pkin(5) * t262 - pkin(4);
t264 = sin(qJ(4));
t330 = t251 * t264;
t336 = t206 * t264;
t397 = -m(6) / 0.2e1;
t68 = t119 * t369 - t263 * t120;
t69 = t263 * t119 + t120 * t369;
t403 = (t343 / 0.2e1 - t342 / 0.2e1) * mrSges(6,3) + (pkin(4) * t336 + qJ(5) * t291) * t397 + (t174 * t68 + t175 * t69 - t206 * t330) * t394;
t350 = cos(pkin(6));
t172 = t207 * t266 + t264 * t350;
t101 = -t172 * t260 + t206 * t262;
t102 = t172 * t262 + t206 * t260;
t63 = t101 * t369 - t263 * t102;
t402 = t63 / 0.2e1;
t401 = t257 / 0.2e1;
t301 = t314 * qJ(5);
t367 = t264 * pkin(4);
t240 = -qJ(5) * t266 + t367;
t328 = t260 * t264;
t176 = t262 * t240 + t249 * t328;
t325 = t262 * t264;
t177 = t260 * t240 - t249 * t325;
t400 = -t176 * t260 + t177 * t262;
t398 = 2 * qJD(4);
t396 = m(6) / 0.2e1;
t395 = m(6) / 0.4e1;
t393 = m(7) / 0.2e1;
t392 = m(7) / 0.4e1;
t391 = mrSges(7,1) / 0.2e1;
t390 = -mrSges(7,2) / 0.2e1;
t389 = mrSges(7,3) / 0.2e1;
t64 = t263 * t101 + t102 * t369;
t388 = t64 / 0.2e1;
t171 = t207 * t264 - t266 * t350;
t91 = t283 * t171;
t387 = t91 / 0.2e1;
t386 = t171 / 0.2e1;
t209 = t283 * t264;
t178 = mrSges(7,2) * t266 - t209 * mrSges(7,3);
t385 = t178 / 0.2e1;
t208 = t282 * t264;
t180 = -mrSges(7,1) * t266 + t208 * mrSges(7,3);
t384 = -t180 / 0.2e1;
t383 = -t206 / 0.2e1;
t382 = -t208 / 0.2e1;
t381 = -t209 / 0.2e1;
t380 = -t210 / 0.2e1;
t379 = -t211 / 0.2e1;
t378 = -t282 / 0.2e1;
t377 = -t283 / 0.2e1;
t354 = t266 * mrSges(5,2);
t241 = t264 * mrSges(5,1) + t354;
t376 = t241 / 0.2e1;
t375 = -t260 / 0.2e1;
t374 = t260 / 0.2e1;
t373 = t262 / 0.2e1;
t372 = t264 / 0.2e1;
t371 = -t266 / 0.2e1;
t365 = Ifges(6,4) * t260;
t364 = Ifges(6,4) * t262;
t363 = Ifges(7,4) * t208;
t362 = Ifges(7,4) * t283;
t361 = Ifges(6,5) * t262;
t360 = Ifges(6,6) * t260;
t357 = t260 * mrSges(6,1);
t356 = t260 * Ifges(6,2);
t355 = t262 * mrSges(6,2);
t353 = t68 * t283;
t352 = t69 * t282;
t238 = -mrSges(6,1) * t262 + mrSges(6,2) * t260;
t351 = t238 - mrSges(5,1);
t347 = mrSges(7,3) * qJD(4);
t335 = t206 * t266;
t36 = (-t208 * t69 - t209 * t68) * t393 + 0.2e1 * (t335 * t392 + (t291 + t335) * t395) * t264;
t346 = qJD(1) * t36;
t345 = t101 * t262;
t344 = t102 * t260;
t341 = t172 * t266;
t258 = t264 ^ 2;
t338 = t206 * t258;
t259 = t266 ^ 2;
t337 = t206 * t259;
t334 = t208 * t180;
t333 = t209 * t178;
t141 = -t208 * mrSges(7,1) - mrSges(7,2) * t209;
t22 = t334 / 0.2e1 - t333 / 0.2e1 + t141 * t371 + (-t209 ^ 2 / 0.2e1 - t208 ^ 2 / 0.2e1) * mrSges(7,3);
t332 = t22 * qJD(2);
t331 = t249 * t264;
t233 = t264 * mrSges(6,1) - mrSges(6,3) * t324;
t329 = t260 * t233;
t231 = -t264 * mrSges(6,2) - mrSges(6,3) * t327;
t326 = t262 * t231;
t322 = t264 * t171;
t217 = t304 * t264;
t321 = t264 * t217;
t320 = t266 * t249;
t316 = -Ifges(7,5) * t209 + Ifges(7,6) * t208;
t315 = -Ifges(7,5) * t282 - Ifges(7,6) * t283;
t250 = -pkin(2) * t349 - pkin(3);
t225 = -t266 * pkin(4) - t264 * qJ(5) + t250;
t163 = t260 * t225 + t262 * t320;
t313 = qJD(4) * t266;
t312 = t141 * qJD(6);
t213 = t262 * t225;
t162 = -t260 * t320 + t213;
t179 = -mrSges(7,2) * t264 - mrSges(7,3) * t210;
t181 = mrSges(7,1) * t264 + mrSges(7,3) * t211;
t143 = -t358 + t359;
t296 = t355 + t357;
t220 = t296 * t266;
t230 = t266 * mrSges(6,2) - mrSges(6,3) * t328;
t232 = -t266 * mrSges(6,1) - mrSges(6,3) * t325;
t284 = t232 * t374 - t262 * t230 / 0.2e1;
t277 = t220 / 0.2e1 + t143 / 0.2e1 + t284;
t142 = mrSges(7,1) * t209 - mrSges(7,2) * t208;
t219 = t296 * t264;
t305 = t219 / 0.2e1 + t142 / 0.2e1;
t140 = -pkin(9) * t325 + t213 + (-t249 * t260 - pkin(5)) * t266;
t149 = -pkin(9) * t328 + t163;
t76 = t140 * t369 - t263 * t149;
t77 = t263 * t140 + t149 * t369;
t150 = t264 * pkin(5) - pkin(9) * t324 + t176;
t161 = -pkin(9) * t327 + t177;
t82 = t150 * t369 - t263 * t161;
t83 = t263 * t150 + t161 * t369;
t14 = t180 * t380 + t181 * t381 + t178 * t379 + t179 * t382 - t277 * t266 + (t326 / 0.2e1 - t329 / 0.2e1 + t305) * t264 + (-t208 * t83 - t209 * t82 - t210 * t76 - t211 * t77 - t218 * t266 + t321) * t393 + ((t163 * t266 + t177 * t264) * t262 + (-t162 * t266 - t176 * t264) * t260 + (t258 - t259) * t249) * t396;
t311 = t14 * qJD(4) + t22 * qJD(6);
t310 = -0.1e1 + t314;
t309 = m(6) * t372;
t164 = mrSges(7,1) * t283 - mrSges(7,2) * t282;
t307 = t164 * t386;
t165 = mrSges(7,1) * t282 + mrSges(7,2) * t283;
t300 = qJD(4) * (t165 + t351);
t298 = 0.2e1 * (t392 + t395) * t172;
t285 = m(7) * (t208 * t63 - t209 * t64);
t30 = -t285 / 0.2e1 + (m(7) * t383 + (t383 + t345 / 0.2e1 + t344 / 0.2e1) * m(6)) * t264;
t34 = m(7) * (t208 * t76 - t209 * t77) + t334 - t333 + (-t260 * t230 - t262 * t232 + m(6) * (-t162 * t262 - t163 * t260)) * t264;
t295 = -t30 * qJD(1) + t34 * qJD(2);
t110 = t206 * t322;
t12 = m(7) * (t63 * t68 + t64 * t69 - t110) + m(6) * (t101 * t119 + t102 * t120 - t110) + m(5) * (t207 - t322 - t341) * t206;
t294 = t12 * qJD(1) + t36 * qJD(3);
t292 = t260 * t101 - t262 * t102;
t92 = t282 * t171;
t15 = (-t208 * t92 - t209 * t91 - t210 * t63 - t211 * t64 + t322 - t341) * t393 + ((-t172 - t292) * t266 - t310 * t322) * t396;
t100 = t171 * t172;
t18 = m(7) * (t63 * t91 + t64 * t92 + t100) + m(6) * (t171 * t292 + t100);
t293 = t18 * qJD(1) + t15 * qJD(3);
t290 = -t162 * t260 + t163 * t262;
t289 = qJD(2) * t141 + qJD(4) * t164;
t288 = t390 * t69 + t391 * t68;
t287 = mrSges(7,1) * t387 + t390 * t92;
t286 = Ifges(7,5) * t211 / 0.2e1 + Ifges(7,6) * t210 / 0.2e1;
t267 = t305 * t172 + t277 * t171 + (t172 * t331 + t176 * t101 + t177 * t102 + (-t290 + t320) * t171) * t396 + (t171 * t218 + t172 * t217 + t63 * t82 + t64 * t83 + t76 * t91 + t77 * t92) * t393 + t101 * t233 / 0.2e1 + t102 * t231 / 0.2e1 + t181 * t402 + t179 * t388 + t180 * t387 + t92 * t385;
t2 = t267 + (t376 - t354 / 0.2e1 + (-mrSges(5,1) / 0.2e1 + t238 / 0.2e1 + t165 / 0.2e1) * t264) * t206 + (t352 / 0.2e1 + t353 / 0.2e1) * mrSges(7,3) + t403;
t136 = -Ifges(7,2) * t209 - t266 * Ifges(7,6) - t363;
t137 = -Ifges(7,4) * t211 - Ifges(7,2) * t210 + Ifges(7,6) * t264;
t197 = Ifges(7,4) * t209;
t138 = -Ifges(7,1) * t208 - t266 * Ifges(7,5) - t197;
t139 = -Ifges(7,1) * t211 - Ifges(7,4) * t210 + Ifges(7,5) * t264;
t204 = Ifges(6,6) * t264 + (-t356 + t364) * t266;
t205 = Ifges(6,5) * t264 + (t262 * Ifges(6,1) - t365) * t266;
t3 = t177 * t230 + t163 * t231 + t176 * t232 + t162 * t233 + t138 * t379 + t217 * t143 + t218 * t142 + t137 * t381 + t136 * t380 + t139 * t382 + t83 * t178 + t77 * t179 + t82 * t180 + t76 * t181 + t250 * t241 + m(7) * (t217 * t218 + t76 * t82 + t77 * t83) + m(6) * (t162 * t176 + t163 * t177) + (Ifges(7,5) * t382 + Ifges(7,6) * t381 + t249 * t220 + t204 * t375 + t205 * t373 + (t361 / 0.2e1 - t360 / 0.2e1 - Ifges(5,4)) * t264) * t264 + (t249 * t219 + (-Ifges(6,3) - Ifges(7,3) + Ifges(5,1) - Ifges(5,2) + m(6) * t249 ^ 2 + Ifges(6,1) * t401 + (-t364 + t356 / 0.2e1) * t260) * t264 + t286 + (Ifges(5,4) + t360 - t361) * t266) * t266;
t281 = t2 * qJD(1) + t3 * qJD(2) + t14 * qJD(3);
t271 = (t208 * t388 + t209 * t402) * mrSges(7,3) + t141 * t386 + t63 * t385 + t64 * t384;
t5 = t271 - t288;
t144 = Ifges(7,2) * t208 - t197;
t145 = -Ifges(7,1) * t209 + t363;
t9 = -t77 * t180 + t76 * t178 + t316 * t371 + t217 * t141 - (-t76 * mrSges(7,3) + t138 / 0.2e1 + t144 / 0.2e1) * t209 + (t77 * mrSges(7,3) - t145 / 0.2e1 + t136 / 0.2e1) * t208;
t280 = t5 * qJD(1) + t9 * qJD(2) + t22 * qJD(3);
t279 = (t208 * t282 + t209 * t283) * t393;
t65 = m(7) * (t208 * t211 + t209 * t210) + 0.4e1 * (t310 * t395 - m(7) / 0.4e1) * t266 * t264;
t278 = t15 * qJD(1) + t14 * qJD(2) + t65 * qJD(3);
t10 = t307 - t287;
t224 = Ifges(7,4) * t282;
t166 = -Ifges(7,2) * t283 - t224;
t167 = -Ifges(7,2) * t282 + t362;
t168 = -Ifges(7,1) * t282 - t362;
t169 = Ifges(7,1) * t283 - t224;
t35 = t251 * t164 - (t167 / 0.2e1 - t168 / 0.2e1) * t283 - (t166 / 0.2e1 + t169 / 0.2e1) * t282;
t274 = t164 * t371;
t38 = t274 + t317;
t268 = -(t138 / 0.4e1 + t144 / 0.4e1) * t282 - (-t145 / 0.4e1 + t136 / 0.4e1) * t283 - (-t174 * mrSges(7,3) / 0.2e1 + t169 / 0.4e1 + t166 / 0.4e1) * t209 + (t175 * t389 - t168 / 0.4e1 + t167 / 0.4e1) * t208 + t174 * t385 + t175 * t384 + t217 * t164 / 0.2e1 + t251 * t141 / 0.2e1 - t266 * t315 / 0.4e1;
t273 = Ifges(7,3) * t372 + t390 * t83 + t391 * t82 - t286;
t7 = t268 - t273;
t276 = t10 * qJD(1) + t7 * qJD(2) + t38 * qJD(3) + t35 * qJD(4);
t269 = (t208 * t377 - t209 * t378) * mrSges(7,3) + t290 * t396 + (t174 * t208 - t175 * t209 - t282 * t77 - t283 * t76) * t393 + t178 * t378 + t180 * t377 - t284;
t20 = (-t355 / 0.2e1 - t357 / 0.2e1 + t249 * t397) * t266 + t269 + t405;
t272 = t292 * t396 + (-t282 * t64 - t283 * t63) * t394;
t28 = t298 + t272;
t48 = (t282 ^ 2 + t283 ^ 2) * mrSges(7,3) + t302 + m(7) * (-t174 * t283 - t175 * t282) + m(6) * t301;
t75 = t279 + (t394 + (t256 / 0.2e1 + t401 - 0.1e1 / 0.2e1) * m(6)) * t264;
t275 = -qJD(1) * t28 + qJD(2) * t20 + qJD(3) * t75 + qJD(4) * t48;
t173 = t249 * t338;
t74 = m(7) * t372 + t309 * t314 + t279 + t309;
t39 = t274 - t317;
t31 = t285 / 0.2e1 + (-t344 - t345) * t309 + (t394 + t397) * t336;
t29 = t298 - t272;
t19 = mrSges(6,2) * t324 / 0.2e1 + mrSges(6,1) * t327 / 0.2e1 + t320 * t396 + t269 - t405;
t11 = t307 + t287;
t8 = t268 + t273;
t6 = t271 + t288;
t4 = t36 * qJD(2) + t15 * qJD(4);
t1 = t267 + t206 * t376 + mrSges(5,1) * t336 / 0.2e1 + mrSges(5,2) * t335 / 0.2e1 + (-t352 - t353) * t389 - (t238 + t165) * t336 / 0.2e1 - t403;
t13 = [t12 * qJD(2) + t18 * qJD(4) (-mrSges(3,2) * t299 - mrSges(3,1) * t303 + m(7) * (-t206 * t321 + t68 * t76 + t69 * t77) + t69 * t178 + t68 * t180 + t120 * t230 + t119 * t232 + m(5) * (t207 * t250 - t249 * t337 - t173) + m(6) * (t119 * t162 + t120 * t163 - t173) + t207 * (-t266 * mrSges(5,1) + t264 * mrSges(5,2)) + m(4) * (-t206 * t261 - t207 * t349) * pkin(2) + t206 * mrSges(4,2) - t207 * mrSges(4,1) + (-t142 - t219) * t336 + (-t337 - t338) * mrSges(5,3)) * qJD(2) + t1 * qJD(4) + t31 * qJD(5) + t6 * qJD(6) + t294, t4, t1 * qJD(2) + t29 * qJD(5) + t11 * qJD(6) + (-t282 * t92 - t283 * t91) * t347 + t172 * t300 + t404 * qJD(4) * t171 + ((t172 * t251 + t174 * t91 + t175 * t92) * t393 + (-pkin(4) * t172 - t171 * t301) * t396) * t398 + t293, qJD(2) * t31 + qJD(4) * t29, t6 * qJD(2) + t11 * qJD(4) + (-mrSges(7,1) * t64 - mrSges(7,2) * t63) * qJD(6); qJD(4) * t2 - qJD(5) * t30 + qJD(6) * t5 - t294, qJD(4) * t3 + qJD(5) * t34 + qJD(6) * t9, t311 - t346, t19 * qJD(5) + t8 * qJD(6) + (Ifges(5,5) + (Ifges(6,2) * t262 + t365) * t375 + (Ifges(6,1) * t260 + t364) * t373 + (-m(6) * pkin(4) + t351) * t249) * t313 + t281 + (-Ifges(5,6) * t264 + t205 * t374 + t204 * t373 + t283 * t139 / 0.2e1 - pkin(4) * t220 + t137 * t378 + t169 * t379 + t218 * t165 + t167 * t380 + t175 * t179 + t174 * t181 + t251 * t143 + m(7) * (t174 * t82 + t175 * t83 + t218 * t251) + mrSges(5,2) * t331 + (m(6) * t400 + t326 - t329) * qJ(5) + (Ifges(6,5) * t260 + Ifges(7,5) * t283 + Ifges(6,6) * t262 - Ifges(7,6) * t282) * t372 + (-t282 * t83 - t283 * t82) * mrSges(7,3) + t400 * mrSges(6,3)) * qJD(4), t19 * qJD(4) + t295, t8 * qJD(4) + (-mrSges(7,1) * t77 - mrSges(7,2) * t76 + t316) * qJD(6) + t280; t4, t311 + t346, t65 * qJD(4), t74 * qJD(5) + t39 * qJD(6) + (t210 * t283 + t211 * t282) * t347 - t404 * t313 + t264 * t300 + ((t266 * t301 - t367) * t396 + (-t174 * t210 - t175 * t211 + t330) * t393) * t398 + t278, t74 * qJD(4), t39 * qJD(4) - t312 + t332; -qJD(2) * t2 - qJD(5) * t28 + qJD(6) * t10 - t293, qJD(5) * t20 + qJD(6) * t7 - t281, qJD(5) * t75 + qJD(6) * t38 - t278, qJD(5) * t48 + qJD(6) * t35, t275 (-mrSges(7,1) * t175 - mrSges(7,2) * t174 + t315) * qJD(6) + t276; qJD(2) * t30 + qJD(4) * t28, -t20 * qJD(4) - t295 + t312, -t75 * qJD(4), qJD(6) * t164 - t275, 0, t289; -t5 * qJD(2) - t10 * qJD(4), -qJD(4) * t7 - qJD(5) * t141 - t280, -t38 * qJD(4) - t332, -qJD(5) * t164 - t276, -t289, 0;];
Cq  = t13;
