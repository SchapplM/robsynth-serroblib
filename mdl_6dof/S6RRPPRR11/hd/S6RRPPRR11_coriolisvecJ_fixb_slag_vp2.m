% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:39:04
% EndTime: 2019-03-09 09:39:34
% DurationCPUTime: 16.81s
% Computational Cost: add. (12659->691), mult. (32789->943), div. (0->0), fcn. (24951->10), ass. (0->341)
t265 = sin(qJ(2));
t259 = sin(pkin(6));
t343 = qJD(1) * t259;
t329 = t265 * t343;
t233 = qJD(5) + t329;
t441 = -t233 / 0.2e1;
t261 = cos(pkin(6));
t342 = qJD(1) * t261;
t248 = qJD(2) + t342;
t258 = sin(pkin(11));
t260 = cos(pkin(11));
t267 = cos(qJ(2));
t328 = t267 * t343;
t193 = -t260 * t248 + t258 * t328;
t264 = sin(qJ(5));
t281 = -t258 * t248 - t260 * t328;
t385 = cos(qJ(5));
t436 = -t193 * t385 + t264 * t281;
t442 = -t436 / 0.2e1;
t416 = Ifges(6,4) * t442 + Ifges(6,6) * t441;
t437 = t264 * t193 + t281 * t385;
t443 = -t437 / 0.2e1;
t431 = Ifges(6,2) * t443;
t406 = t431 + t416;
t120 = qJD(6) - t437;
t394 = t120 / 0.2e1;
t263 = sin(qJ(6));
t266 = cos(qJ(6));
t106 = t233 * t263 + t266 * t436;
t396 = t106 / 0.2e1;
t105 = t233 * t266 - t263 * t436;
t398 = t105 / 0.2e1;
t410 = Ifges(7,5) * t396 + Ifges(7,6) * t398 + Ifges(7,3) * t394;
t459 = t406 + 0.2e1 * t410;
t449 = t437 / 0.2e1;
t451 = t233 / 0.2e1;
t444 = Ifges(6,4) * t449 + Ifges(6,5) * t451;
t450 = t436 / 0.2e1;
t447 = Ifges(6,1) * t450;
t405 = t447 + t444;
t351 = t259 * t265;
t402 = pkin(3) + pkin(8);
t287 = (-pkin(4) * t260 - t402) * t351;
t280 = qJD(2) * t287;
t336 = pkin(1) * t342;
t242 = t267 * t336;
t232 = qJD(2) * t242;
t234 = t248 * qJD(3);
t345 = t232 + t234;
t122 = qJD(1) * t280 + t345;
t283 = -t258 * t385 - t264 * t260;
t276 = t283 * t351;
t274 = qJD(2) * t276;
t93 = -qJD(1) * t274 + qJD(5) * t437;
t384 = Ifges(6,1) * t93;
t262 = -pkin(2) - qJ(4);
t134 = t248 * t262 + t329 * t402 + qJD(3) - t242;
t323 = -qJ(3) * t265 - pkin(1);
t174 = (t262 * t267 + t323) * t259;
t161 = qJD(1) * t174;
t85 = t260 * t134 - t161 * t258;
t66 = pkin(4) * t329 + pkin(9) * t193 + t85;
t86 = t258 * t134 + t260 * t161;
t68 = pkin(9) * t281 + t86;
t30 = t264 * t66 + t385 * t68;
t352 = t258 * t265;
t282 = (pkin(4) * t267 - pkin(9) * t352) * t259;
t278 = qJD(2) * t282;
t324 = qJD(2) * t343;
t319 = t265 * t324;
t231 = pkin(2) * t319;
t294 = -qJ(3) * t267 + qJ(4) * t265;
t339 = qJD(3) * t265;
t273 = (qJD(2) * t294 - qJD(4) * t267 - t339) * t259;
t127 = qJD(1) * t273 + t231;
t252 = t261 * t265 * pkin(1);
t350 = t259 * t267;
t275 = (t350 * t402 + t252) * qJD(2);
t147 = qJD(1) * t275 - t248 * qJD(4);
t74 = -t258 * t127 + t260 * t147;
t62 = qJD(1) * t278 + t74;
t292 = t260 * t319;
t75 = t260 * t127 + t258 * t147;
t69 = pkin(9) * t292 + t75;
t8 = -qJD(5) * t30 - t264 * t69 + t385 * t62;
t458 = t384 / 0.2e1 + t122 * mrSges(6,2) - t8 * mrSges(6,3);
t28 = t233 * pkin(10) + t30;
t210 = pkin(8) * t328 + t265 * t336;
t186 = pkin(3) * t328 + t210;
t235 = t248 * qJ(3);
t152 = t235 + qJD(4) + t186;
t113 = -pkin(4) * t281 + t152;
t51 = -pkin(5) * t437 - pkin(10) * t436 + t113;
t12 = -t263 * t28 + t266 * t51;
t330 = t385 * t260;
t296 = t330 * t351;
t286 = qJD(1) * t296;
t293 = t258 * t319;
t94 = -qJD(2) * t286 + qJD(5) * t436 + t264 * t293;
t31 = pkin(5) * t94 - pkin(10) * t93 + t122;
t318 = t267 * t324;
t325 = qJD(5) * t385;
t338 = qJD(5) * t264;
t7 = t264 * t62 + t66 * t325 - t338 * t68 + t385 * t69;
t5 = pkin(10) * t318 + t7;
t1 = qJD(6) * t12 + t263 * t31 + t266 * t5;
t13 = t263 * t51 + t266 * t28;
t2 = -qJD(6) * t13 - t263 * t5 + t266 * t31;
t316 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t383 = Ifges(6,4) * t93;
t46 = -qJD(6) * t106 - t263 * t93 + t266 * t318;
t42 = Ifges(7,6) * t46;
t45 = qJD(6) * t105 + t263 * t318 + t266 * t93;
t43 = Ifges(7,5) * t45;
t9 = Ifges(7,3) * t94 + t42 + t43;
t457 = t316 + t122 * mrSges(6,1) - t7 * mrSges(6,3) + t9 / 0.2e1 - t383 / 0.2e1;
t353 = t258 * t264;
t320 = t351 * t353;
t175 = -qJD(1) * t320 + t286;
t216 = -t258 * t338 + t260 * t325;
t348 = t175 + t216;
t176 = qJD(1) * t276;
t215 = -t258 * t325 - t260 * t338;
t346 = t215 + t176;
t29 = -t264 * t68 + t385 * t66;
t454 = -t113 * mrSges(6,2) + t29 * mrSges(6,3);
t452 = -t113 * mrSges(6,1) - t12 * mrSges(7,1) + t13 * mrSges(7,2) + t30 * mrSges(6,3);
t237 = pkin(2) * t329;
t178 = t294 * t343 + t237;
t112 = t260 * t178 + t258 * t186;
t349 = t260 * t265;
t322 = pkin(9) * t259 * t349;
t101 = qJD(1) * t322 + t112;
t378 = -pkin(9) + t262;
t226 = t378 * t258;
t227 = t378 * t260;
t284 = -t264 * t226 + t227 * t385;
t111 = -t258 * t178 + t260 * t186;
t97 = qJD(1) * t282 + t111;
t425 = qJD(4) * t283 + qJD(5) * t284 - t385 * t101 - t264 * t97;
t148 = qJD(1) * t287 + t242;
t448 = qJD(3) - t148;
t446 = -pkin(10) * t328 + t425;
t445 = pkin(5) * t348 - pkin(10) * t346 + t448;
t432 = -t343 / 0.2e1;
t440 = mrSges(4,1) + mrSges(3,3);
t439 = mrSges(4,2) - mrSges(3,1);
t167 = t226 * t385 + t264 * t227;
t420 = t330 - t353;
t424 = -qJD(4) * t420 - qJD(5) * t167 + t264 * t101 - t385 * t97;
t438 = -t29 * t346 - t30 * t348;
t435 = -t248 / 0.2e1;
t434 = t248 / 0.2e1;
t253 = t258 * pkin(4) + qJ(3);
t154 = -pkin(5) * t283 - pkin(10) * t420 + t253;
t99 = t154 * t266 - t167 * t263;
t430 = qJD(6) * t99 + t263 * t445 + t266 * t446;
t100 = t154 * t263 + t167 * t266;
t429 = -qJD(6) * t100 - t263 * t446 + t266 * t445;
t303 = t12 * t266 + t13 * t263;
t426 = t303 * mrSges(7,3);
t110 = mrSges(6,1) * t233 - mrSges(6,3) * t436;
t54 = -mrSges(7,1) * t105 + mrSges(7,2) * t106;
t354 = t110 - t54;
t423 = pkin(5) * t328 - t424;
t315 = t1 * t266 - t2 * t263;
t209 = pkin(8) * t329 - t242;
t419 = -qJD(3) - t209;
t27 = -t233 * pkin(5) - t29;
t305 = Ifges(7,5) * t266 - Ifges(7,6) * t263;
t370 = Ifges(7,4) * t266;
t307 = -Ifges(7,2) * t263 + t370;
t371 = Ifges(7,4) * t263;
t310 = Ifges(7,1) * t266 - t371;
t312 = mrSges(7,1) * t263 + mrSges(7,2) * t266;
t102 = Ifges(7,4) * t105;
t34 = Ifges(7,1) * t106 + Ifges(7,5) * t120 + t102;
t386 = t266 / 0.2e1;
t372 = Ifges(7,4) * t106;
t33 = Ifges(7,2) * t105 + Ifges(7,6) * t120 + t372;
t409 = -t33 / 0.2e1;
t271 = t263 * t409 + t27 * t312 + t305 * t394 + t307 * t398 + t310 * t396 + t34 * t386;
t418 = t405 + t271 + t444;
t417 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t93 - Ifges(6,6) * t94;
t173 = -pkin(2) * t248 - t419;
t200 = (-pkin(2) * t267 + t323) * t259;
t189 = qJD(1) * t200;
t236 = Ifges(3,4) * t328;
t369 = Ifges(4,6) * t267;
t415 = -(Ifges(4,4) / 0.2e1 - Ifges(3,5) / 0.2e1) * t248 + t173 * mrSges(4,1) + t209 * mrSges(3,3) + t29 * mrSges(6,1) + t85 * mrSges(5,1) - t193 * Ifges(5,5) + t281 * Ifges(5,6) + Ifges(3,5) * t434 + t236 / 0.2e1 + Ifges(4,4) * t435 + (-t265 * Ifges(4,2) - t369) * t432 + t233 * Ifges(6,3) + t436 * Ifges(6,5) + t437 * Ifges(6,6) - t189 * mrSges(4,3) - t30 * mrSges(6,2) - t86 * mrSges(5,2) + (Ifges(5,3) + Ifges(3,1)) * t329 / 0.2e1;
t249 = pkin(8) * t351;
t381 = pkin(1) * t267;
t331 = -pkin(2) - t381;
t162 = pkin(3) * t351 + t249 + (-qJ(4) + t331) * t261;
t107 = t260 * t162 - t174 * t258;
t214 = -t258 * t350 + t261 * t260;
t83 = pkin(4) * t351 - pkin(9) * t214 + t107;
t108 = t258 * t162 + t260 * t174;
t213 = -t261 * t258 - t260 * t350;
t91 = pkin(9) * t213 + t108;
t376 = t264 * t83 + t385 * t91;
t341 = qJD(2) * t265;
t327 = t259 * t341;
t239 = pkin(2) * t327;
t139 = t239 + t273;
t164 = -t261 * qJD(4) + t275;
t95 = -t258 * t139 + t260 * t164;
t73 = t278 + t95;
t96 = t260 * t139 + t258 * t164;
t82 = qJD(2) * t322 + t96;
t17 = -qJD(5) * t376 - t264 * t82 + t385 * t73;
t183 = -t235 - t210;
t299 = t258 * t85 - t260 * t86;
t375 = Ifges(5,4) * t258;
t308 = Ifges(5,2) * t260 + t375;
t374 = Ifges(5,4) * t260;
t311 = Ifges(5,1) * t258 + t374;
t388 = t260 / 0.2e1;
t389 = t258 / 0.2e1;
t412 = Ifges(4,5) / 0.2e1;
t414 = t299 * mrSges(5,3) - (t412 - Ifges(3,6) / 0.2e1) * t248 - (-Ifges(5,4) * t193 + Ifges(5,2) * t281 + Ifges(5,6) * t329) * t388 - (-Ifges(5,1) * t193 + Ifges(5,4) * t281 + Ifges(5,5) * t329) * t389 + t210 * mrSges(3,3) - t183 * mrSges(4,1) - Ifges(3,6) * t435 - Ifges(4,5) * t434 + t193 * t311 / 0.2e1 + t189 * mrSges(4,2) - t152 * (-mrSges(5,1) * t260 + mrSges(5,2) * t258) - t281 * t308 / 0.2e1 + ((-Ifges(3,2) - Ifges(4,3)) * t267 + (-Ifges(3,4) - Ifges(4,6)) * t265) * t432;
t411 = Ifges(6,2) / 0.2e1;
t408 = t45 / 0.2e1;
t407 = t46 / 0.2e1;
t404 = -t94 / 0.2e1;
t403 = t94 / 0.2e1;
t401 = pkin(1) * mrSges(3,1);
t400 = pkin(1) * mrSges(3,2);
t399 = -t105 / 0.2e1;
t397 = -t106 / 0.2e1;
t395 = -t120 / 0.2e1;
t285 = t213 * t385 - t264 * t214;
t393 = t285 / 0.2e1;
t145 = t264 * t213 + t214 * t385;
t392 = t145 / 0.2e1;
t391 = t213 / 0.2e1;
t390 = t214 / 0.2e1;
t387 = t263 / 0.2e1;
t382 = Ifges(6,4) * t94;
t20 = -mrSges(7,1) * t46 + mrSges(7,2) * t45;
t78 = mrSges(6,1) * t318 - t93 * mrSges(6,3);
t377 = t20 - t78;
t356 = t258 * Ifges(5,5);
t355 = t260 * Ifges(5,6);
t347 = t248 * t439 + t329 * t440;
t337 = t261 * t381;
t243 = qJD(2) * t337;
t255 = t261 * qJD(3);
t344 = t243 + t255;
t218 = pkin(8) * t350 + t252;
t340 = qJD(2) * t267;
t333 = 0.3e1 / 0.2e1 * Ifges(4,6) + 0.3e1 / 0.2e1 * Ifges(3,4);
t132 = -mrSges(5,1) * t281 - mrSges(5,2) * t193;
t205 = -mrSges(4,1) * t328 - t248 * mrSges(4,3);
t70 = -mrSges(6,1) * t437 + mrSges(6,2) * t436;
t332 = t205 - t132 - t70;
t199 = -t261 * qJ(3) - t218;
t326 = t259 * t340;
t44 = t94 * mrSges(6,1) + t93 * mrSges(6,2);
t321 = t402 * t351;
t177 = pkin(3) * t350 - t199;
t314 = t1 * t263 + t2 * t266;
t313 = -mrSges(7,1) * t266 + mrSges(7,2) * t263;
t309 = Ifges(7,1) * t263 + t370;
t306 = Ifges(7,2) * t266 + t371;
t304 = Ifges(7,5) * t263 + Ifges(7,6) * t266;
t302 = -t12 * t263 + t13 * t266;
t25 = mrSges(7,1) * t94 - mrSges(7,3) * t45;
t26 = -mrSges(7,2) * t94 + mrSges(7,3) * t46;
t301 = -t263 * t25 + t266 * t26;
t300 = t258 * t75 + t260 * t74;
t36 = pkin(10) * t351 + t376;
t135 = -pkin(4) * t213 + t177;
t57 = -pkin(5) * t285 - pkin(10) * t145 + t135;
t22 = t263 * t57 + t266 * t36;
t21 = -t263 * t36 + t266 * t57;
t58 = -mrSges(7,2) * t120 + mrSges(7,3) * t105;
t59 = mrSges(7,1) * t120 - mrSges(7,3) * t106;
t298 = -t263 * t59 + t266 * t58;
t297 = -t263 * t58 - t266 * t59;
t295 = qJD(2) * t321;
t211 = -pkin(8) * t327 + t243;
t109 = -mrSges(6,2) * t233 + mrSges(6,3) * t437;
t291 = -t109 - t298;
t38 = -t264 * t91 + t385 * t83;
t117 = -t145 * t263 + t266 * t351;
t118 = t145 * t266 + t263 * t351;
t16 = t264 * t73 + t83 * t325 - t338 * t91 + t385 * t82;
t197 = -pkin(8) * t319 + t232;
t184 = -mrSges(5,1) * t292 + mrSges(5,2) * t293;
t279 = (-qJ(3) * t340 - t339) * t259;
t212 = t218 * qJD(2);
t172 = -t197 - t234;
t198 = qJD(1) * t212;
t277 = -t197 * mrSges(3,2) - t172 * mrSges(4,3) + t198 * t439;
t138 = t280 + t344;
t272 = t416 + t459;
t230 = Ifges(3,5) * t318;
t229 = Ifges(4,5) * t319;
t228 = Ifges(6,3) * t318;
t217 = -t249 + t337;
t208 = -qJ(3) * t328 + t237;
t207 = (mrSges(4,2) * t267 - mrSges(4,3) * t265) * t343;
t204 = -t248 * mrSges(3,2) + mrSges(3,3) * t328;
t202 = t261 * t331 + t249;
t192 = -t211 - t255;
t191 = (-mrSges(5,2) * t267 + mrSges(5,3) * t349) * t324;
t190 = (mrSges(5,1) * t267 - mrSges(5,3) * t352) * t324;
t187 = t239 + t279;
t185 = -qJD(1) * t321 + t242;
t169 = qJD(1) * t279 + t231;
t163 = -t295 + t344;
t160 = mrSges(5,1) * t329 + mrSges(5,3) * t193;
t159 = -mrSges(5,2) * t329 + mrSges(5,3) * t281;
t151 = (Ifges(5,5) * t267 + t265 * t311) * t324;
t150 = (Ifges(5,6) * t267 + t265 * t308) * t324;
t146 = -qJD(1) * t295 + t345;
t143 = -t266 * t176 + t263 * t328;
t142 = t263 * t176 + t266 * t328;
t141 = -t175 * t266 + t248 * t263;
t140 = t175 * t263 + t248 * t266;
t104 = qJD(5) * t145 + (-t296 + t320) * qJD(2);
t103 = qJD(5) * t285 - t274;
t79 = -mrSges(6,2) * t318 - t94 * mrSges(6,3);
t71 = pkin(5) * t436 - pkin(10) * t437;
t56 = -qJD(6) * t118 - t263 * t103 + t266 * t326;
t55 = qJD(6) * t117 + t266 * t103 + t263 * t326;
t41 = Ifges(6,5) * t318 - t382 + t384;
t40 = -Ifges(6,2) * t94 + Ifges(6,6) * t318 + t383;
t37 = pkin(5) * t104 - pkin(10) * t103 + t138;
t35 = -pkin(5) * t351 - t38;
t19 = t263 * t71 + t266 * t29;
t18 = -t263 * t29 + t266 * t71;
t15 = -pkin(5) * t326 - t17;
t14 = pkin(10) * t326 + t16;
t11 = t45 * Ifges(7,1) + t46 * Ifges(7,4) + t94 * Ifges(7,5);
t10 = t45 * Ifges(7,4) + t46 * Ifges(7,2) + t94 * Ifges(7,6);
t6 = -pkin(5) * t318 - t8;
t4 = -qJD(6) * t22 - t14 * t263 + t266 * t37;
t3 = qJD(6) * t21 + t14 * t266 + t263 * t37;
t23 = [((Ifges(5,5) * t390 + Ifges(5,6) * t391 + Ifges(6,5) * t392 + Ifges(6,6) * t393 - t200 * mrSges(4,3) - t217 * mrSges(3,3) + t202 * mrSges(4,1) + (-Ifges(4,4) + Ifges(3,5) / 0.2e1) * t261 + (t267 * t333 - 0.2e1 * t400) * t259) * t340 + (-t200 * mrSges(4,2) - t218 * mrSges(3,3) + t199 * mrSges(4,1) + (Ifges(5,1) * t214 + Ifges(5,4) * t213) * t389 + (Ifges(5,4) * t214 + Ifges(5,2) * t213) * t388 + (-Ifges(3,6) + t412) * t261 + (-0.2e1 * t401 + (0.3e1 / 0.2e1 * t356 + 0.3e1 / 0.2e1 * t355 - t333) * t265) * t259 + (0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + Ifges(6,3) / 0.2e1) * t350) * t341) * t343 + m(6) * (t113 * t138 + t122 * t135 + t16 * t30 + t17 * t29 + t376 * t7 + t38 * t8) + t376 * t79 + ((-mrSges(4,1) * t172 + mrSges(4,2) * t169 + mrSges(3,3) * t197) * t267 + (-t75 * mrSges(5,2) + t74 * mrSges(5,1) - t169 * mrSges(4,3) + t228 / 0.2e1 + t440 * t198 + t417) * t265 + (-t265 * t414 + t267 * t415) * qJD(2)) * t259 + m(3) * (t197 * t218 - t198 * t217 + t209 * t212 + t210 * t211) + m(4) * (t169 * t200 + t172 * t199 + t173 * t212 + t183 * t192 + t187 * t189 + t198 * t202) + m(5) * (t107 * t74 + t108 * t75 + t146 * t177 + t152 * t163 + t85 * t95 + t86 * t96) + m(7) * (t1 * t22 + t12 * t4 + t13 * t3 + t15 * t27 + t2 * t21 + t35 * t6) + (t213 * t75 - t214 * t74) * mrSges(5,3) + (t1 * t117 - t118 * t2 - t12 * t55 + t13 * t56) * mrSges(7,3) + (Ifges(6,4) * t404 + t458) * t145 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t394 + (Ifges(7,5) * t118 + Ifges(7,6) * t117) * t403 + t151 * t390 + t150 * t391 + t41 * t392 + t40 * t393 + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t407 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t398 + (0.2e1 * t405 - t454) * t103 + (-Ifges(6,4) * t450 - Ifges(6,2) * t449 - Ifges(6,6) * t451 - t452 + t459) * t104 + t347 * t212 + (Ifges(7,1) * t118 + Ifges(7,4) * t117) * t408 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t396 + (-Ifges(7,5) * t408 + Ifges(6,2) * t404 - Ifges(7,6) * t407 - Ifges(7,3) * t403 - t457) * t285 + t146 * (-mrSges(5,1) * t213 + mrSges(5,2) * t214) + t211 * t204 + t192 * t205 + t187 * t207 + t107 * t190 + t108 * t191 + t177 * t184 + t96 * t159 + t95 * t160 + t163 * t132 + t135 * t44 + t138 * t70 + t117 * t10 / 0.2e1 + t6 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) + t118 * t11 / 0.2e1 + t16 * t109 + t17 * t110 + t38 * t78 + t3 * t58 + t4 * t59 + t15 * t54 + t55 * t34 / 0.2e1 + t27 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t56 * t33 / 0.2e1 + t35 * t20 + t21 * t25 + t22 * t26 + (t229 / 0.2e1 + t230 / 0.2e1 + t277) * t261; -t332 * qJD(3) + (-Ifges(6,1) * t442 - Ifges(6,4) * t443 - Ifges(6,5) * t441 + t405) * t176 + (t447 + t418) * t215 + (-t205 + t204) * t209 + (((t400 - t369 / 0.2e1) * t343 - t236 / 0.2e1 - t415) * t267 + ((t401 + (Ifges(3,4) / 0.2e1 - t356 / 0.2e1 - t355 / 0.2e1 + Ifges(4,6) / 0.2e1) * t265) * t343 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1) * t328 + t414) * t265 + ((Ifges(5,5) * t388 - Ifges(5,6) * t258 / 0.2e1 + Ifges(6,5) * t420 / 0.2e1 + Ifges(6,6) * t283 / 0.2e1 - Ifges(4,4) - pkin(2) * mrSges(4,1)) * t267 + (-Ifges(3,6) - qJ(3) * mrSges(4,1) + (Ifges(5,1) * t260 - t375) * t389 + (-Ifges(5,2) * t258 + t374) * t388) * t265) * qJD(2)) * t343 - (-t6 * t312 + t305 * t404 - t45 * t310 / 0.2e1 - t46 * t307 / 0.2e1 + t382 / 0.2e1 - t41 / 0.2e1 - t266 * t11 / 0.2e1 + t10 * t387 + t314 * mrSges(7,3) + (mrSges(7,3) * t302 + t27 * t313 + t304 * t394 + t306 * t398 + t309 * t396 + t33 * t386 + t34 * t387) * qJD(6) - t458) * t420 + t438 * mrSges(6,3) + t277 + t142 * t409 - t347 * t210 + ((-t266 * t215 + t143) * mrSges(7,3) + t348 * mrSges(7,1)) * t12 + (mrSges(6,1) * t348 + mrSges(6,2) * t346) * t113 + ((-t263 * t215 - t142) * mrSges(7,3) - t348 * mrSges(7,2)) * t13 - (-t40 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1 + (t411 + Ifges(7,3) / 0.2e1) * t94 + t457) * t283 + t253 * t44 - t208 * t207 + qJ(3) * t184 - t185 * t132 - t112 * t159 - t111 * t160 + t167 * t79 - t148 * t70 - t27 * (-mrSges(7,1) * t142 + mrSges(7,2) * t143) - t143 * t34 / 0.2e1 + t99 * t25 + t100 * t26 + (Ifges(7,5) * t143 + Ifges(7,6) * t142) * t395 + (Ifges(7,4) * t143 + Ifges(7,2) * t142) * t399 - t377 * t284 + (t1 * t100 + t12 * t429 + t13 * t430 + t2 * t99 + t27 * t423 - t284 * t6) * m(7) + (t113 * t448 + t122 * t253 + t167 * t7 + t284 * t8 + t424 * t29 + t425 * t30) * m(6) + (-Ifges(7,5) * t397 - Ifges(7,6) * t399 - Ifges(7,3) * t395 + 0.2e1 * t406 + t410) * t175 + (-pkin(2) * t198 - qJ(3) * t172 - t173 * t210 + t183 * t419 - t189 * t208) * m(4) + t423 * t54 + t424 * t110 + t425 * t109 + t229 + t230 + t429 * t59 + t430 * t58 + (t431 + t272) * t216 + (-t111 * t85 - t112 * t86 + qJ(3) * t146 + t300 * t262 + (-t258 * t86 - t260 * t85) * qJD(4) + (-t185 + qJD(3)) * t152) * m(5) + (t146 * mrSges(5,1) - t150 / 0.2e1 + t262 * t191 - qJD(4) * t159 - t75 * mrSges(5,3)) * t258 + (Ifges(7,1) * t143 + Ifges(7,4) * t142) * t397 + (t146 * mrSges(5,2) + t151 / 0.2e1 - t74 * mrSges(5,3) + t262 * t190 - qJD(4) * t160) * t260; t175 * t109 - t140 * t59 - t141 * t58 + t260 * t190 + t258 * t191 - t377 * t420 + t332 * t248 - t291 * t216 + (mrSges(4,1) * t340 + (t159 * t260 - t160 * t258 + t207) * t265) * t343 - (t297 * qJD(6) + t301 + t79) * t283 + t354 * t346 + (-t420 * t6 + t302 * t216 - (-qJD(6) * t303 + t315) * t283 - t12 * t140 - t13 * t141 - t346 * t27) * m(7) + (-t113 * t248 - t283 * t7 + t420 * t8 - t438) * m(6) + (-t152 * t248 - t299 * t329 + t300) * m(5) + (t183 * t248 + t189 * t329 + t198) * m(4); t298 * qJD(6) + t354 * t436 + t291 * t437 - t281 * t159 - t193 * t160 + t266 * t25 + t263 * t26 + t184 + t44 + (t120 * t302 - t436 * t27 + t314) * m(7) + (t29 * t436 - t30 * t437 + t122) * m(6) + (-t193 * t85 - t281 * t86 + t146) * m(5); t354 * t30 + ((-Ifges(6,1) / 0.2e1 + t411) * t436 + t426 - t418 + t454) * t437 + t417 + (-t272 + t452) * t436 + t306 * t407 + t309 * t408 + t304 * t403 + t10 * t386 + t11 * t387 + t6 * t313 + t315 * mrSges(7,3) - t29 * t109 - t19 * t58 - t18 * t59 + (-pkin(5) * t6 - t12 * t18 - t13 * t19 - t27 * t30) * m(7) - pkin(5) * t20 + (t271 - t426) * qJD(6) + (m(7) * t315 + t301 + (-m(7) * t303 + t297) * qJD(6)) * pkin(10) + t228; -t27 * (mrSges(7,1) * t106 + mrSges(7,2) * t105) + (Ifges(7,1) * t105 - t372) * t397 + t33 * t396 + (Ifges(7,5) * t105 - Ifges(7,6) * t106) * t395 - t12 * t58 + t13 * t59 + (t105 * t12 + t106 * t13) * mrSges(7,3) + t316 + t9 + (-Ifges(7,2) * t106 + t102 + t34) * t399;];
tauc  = t23(:);
