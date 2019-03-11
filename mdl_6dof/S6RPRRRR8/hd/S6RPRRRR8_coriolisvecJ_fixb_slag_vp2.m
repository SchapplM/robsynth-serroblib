% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:46
% EndTime: 2019-03-09 07:20:08
% DurationCPUTime: 10.62s
% Computational Cost: add. (12803->610), mult. (27146->852), div. (0->0), fcn. (18417->8), ass. (0->287)
t284 = cos(qJ(5));
t391 = cos(qJ(4));
t324 = qJD(4) * t391;
t316 = pkin(3) * t324;
t285 = -pkin(1) - pkin(7);
t254 = qJD(1) * t285 + qJD(2);
t282 = sin(qJ(3));
t347 = qJD(1) * t282;
t217 = -pkin(8) * t347 + t254 * t282;
t281 = sin(qJ(4));
t211 = t281 * t217;
t392 = cos(qJ(3));
t242 = t392 * t254;
t326 = qJD(1) * t392;
t218 = -pkin(8) * t326 + t242;
t161 = t218 * t391 - t211;
t330 = t391 * t282;
t226 = -qJD(1) * t330 - t281 * t326;
t313 = t391 * t392;
t329 = t281 * t347;
t228 = qJD(1) * t313 - t329;
t190 = t228 * pkin(4) - t226 * pkin(9);
t315 = pkin(3) * t326;
t165 = t315 + t190;
t280 = sin(qJ(5));
t91 = t284 * t161 + t280 * t165;
t460 = t284 * t316 - t91;
t90 = -t161 * t280 + t284 * t165;
t459 = -t280 * t316 - t90;
t390 = pkin(3) * t281;
t268 = pkin(9) + t390;
t383 = -pkin(10) - t268;
t320 = qJD(5) * t383;
t357 = t226 * t280;
t341 = pkin(10) * t357;
t458 = t280 * t320 + t341 + t460;
t356 = t226 * t284;
t312 = t228 * pkin(5) - pkin(10) * t356;
t457 = t284 * t320 - t312 + t459;
t413 = -pkin(10) - pkin(9);
t333 = qJD(5) * t413;
t213 = qJD(3) * pkin(3) + t218;
t152 = t391 * t213 - t211;
t95 = t284 * t152 + t280 * t190;
t456 = t280 * t333 + t341 - t95;
t94 = -t152 * t280 + t284 * t190;
t455 = t284 * t333 - t312 - t94;
t279 = sin(qJ(6));
t283 = cos(qJ(6));
t300 = t279 * t280 - t283 * t284;
t159 = t300 * t226;
t425 = qJD(5) + qJD(6);
t194 = t425 * t300;
t454 = t159 - t194;
t276 = qJD(3) + qJD(4);
t143 = -t276 * pkin(4) - t152;
t202 = -t228 * t280 + t276 * t284;
t107 = -t202 * pkin(5) + t143;
t203 = t228 * t284 + t276 * t280;
t135 = t202 * t279 + t203 * t283;
t212 = t391 * t217;
t153 = t281 * t213 + t212;
t144 = t276 * pkin(9) + t153;
t249 = pkin(3) * t347 + qJD(1) * qJ(2);
t155 = -pkin(4) * t226 - pkin(9) * t228 + t249;
t78 = t144 * t284 + t155 * t280;
t64 = pkin(10) * t202 + t78;
t366 = t279 * t64;
t221 = qJD(5) - t226;
t77 = -t144 * t280 + t284 * t155;
t63 = -pkin(10) * t203 + t77;
t57 = pkin(5) * t221 + t63;
t19 = t283 * t57 - t366;
t363 = t283 * t64;
t20 = t279 * t57 + t363;
t332 = t281 * t392;
t236 = t332 + t330;
t184 = t276 * t236 * qJD(1);
t120 = qJD(5) * t202 - t184 * t284;
t214 = t218 * qJD(3);
t346 = qJD(3) * t282;
t297 = (pkin(8) * qJD(1) - t254) * t346;
t82 = qJD(4) * t152 + t391 * t214 + t281 * t297;
t289 = t276 * t313;
t342 = qJD(1) * qJD(3);
t323 = t282 * t342;
t185 = qJD(1) * t289 - qJD(4) * t329 - t281 * t323;
t277 = qJD(1) * qJD(2);
t325 = qJD(3) * t392;
t311 = qJD(1) * t325;
t243 = pkin(3) * t311 + t277;
t97 = pkin(4) * t185 + pkin(9) * t184 + t243;
t18 = -qJD(5) * t78 - t280 * t82 + t284 * t97;
t12 = pkin(5) * t185 - pkin(10) * t120 + t18;
t121 = -qJD(5) * t203 + t184 * t280;
t343 = qJD(5) * t284;
t344 = qJD(5) * t280;
t17 = -t144 * t344 + t155 * t343 + t280 * t97 + t284 * t82;
t13 = pkin(10) * t121 + t17;
t3 = qJD(6) * t19 + t12 * t279 + t13 * t283;
t318 = t283 * t202 - t203 * t279;
t36 = qJD(6) * t318 + t120 * t283 + t121 * t279;
t37 = -qJD(6) * t135 - t120 * t279 + t121 * t283;
t340 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t185;
t376 = Ifges(7,4) * t135;
t4 = -qJD(6) * t20 + t12 * t283 - t13 * t279;
t216 = qJD(6) + t221;
t400 = -t216 / 0.2e1;
t408 = -t135 / 0.2e1;
t453 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t340 + (Ifges(7,5) * t318 - Ifges(7,6) * t135) * t400 + (t135 * t20 + t19 * t318) * mrSges(7,3) - t107 * (mrSges(7,1) * t135 + mrSges(7,2) * t318) + (Ifges(7,1) * t318 - t376) * t408;
t452 = t276 * Ifges(5,5);
t451 = t276 * Ifges(5,6);
t450 = t203 * Ifges(6,5) + t135 * Ifges(7,5) + t202 * Ifges(6,6) + Ifges(7,6) * t318 + t221 * Ifges(6,3) + t216 * Ifges(7,3);
t237 = t279 * t284 + t280 * t283;
t158 = t237 * t226;
t195 = t425 * t237;
t449 = t158 - t195;
t370 = t228 * mrSges(5,3);
t348 = mrSges(5,1) * t276 + mrSges(6,1) * t202 - mrSges(6,2) * t203 - t370;
t160 = t218 * t281 + t212;
t345 = qJD(4) * t281;
t448 = pkin(3) * t345 - t160;
t201 = Ifges(6,4) * t202;
t111 = Ifges(6,1) * t203 + Ifges(6,5) * t221 + t201;
t219 = Ifges(5,4) * t226;
t447 = t228 * Ifges(5,1) + t284 * t111 + t219 + t452;
t446 = t17 * t284 - t18 * t280;
t445 = qJ(2) * (m(4) + m(3)) + mrSges(3,3);
t131 = Ifges(7,4) * t318;
t444 = -Ifges(7,2) * t135 + t131;
t419 = t36 / 0.2e1;
t418 = t37 / 0.2e1;
t406 = t185 / 0.2e1;
t384 = pkin(8) - t285;
t229 = t383 * t280;
t275 = t284 * pkin(10);
t353 = t268 * t284;
t230 = t275 + t353;
t186 = t229 * t283 - t230 * t279;
t441 = qJD(6) * t186 + t457 * t279 + t458 * t283;
t187 = t229 * t279 + t230 * t283;
t440 = -qJD(6) * t187 - t458 * t279 + t457 * t283;
t250 = t413 * t280;
t251 = pkin(9) * t284 + t275;
t204 = t250 * t283 - t251 * t279;
t435 = qJD(6) * t204 + t455 * t279 + t456 * t283;
t205 = t250 * t279 + t251 * t283;
t434 = -qJD(6) * t205 - t456 * t279 + t455 * t283;
t174 = t300 * t236;
t197 = -t281 * t346 - t282 * t345 + t289;
t433 = t300 * qJD(1) + t174 * t425 - t237 * t197;
t432 = -t237 * qJD(1) - t195 * t236 - t197 * t300;
t382 = mrSges(5,3) * t184;
t62 = -mrSges(6,1) * t121 + mrSges(6,2) * t120;
t431 = t62 - t382;
t309 = mrSges(6,1) * t280 + mrSges(6,2) * t284;
t430 = t143 * t309;
t235 = t281 * t282 - t313;
t173 = t237 * t235;
t266 = t282 * pkin(3) + qJ(2);
t191 = pkin(4) * t236 + pkin(9) * t235 + t266;
t244 = t384 * t282;
t245 = t384 * t392;
t200 = -t244 * t391 - t281 * t245;
t192 = t284 * t200;
t116 = t280 * t191 + t192;
t429 = t281 * t244 - t391 * t245;
t215 = pkin(5) * t357;
t337 = pkin(5) * t344;
t428 = t337 - t215 + t448;
t426 = -t280 * t77 + t284 * t78;
t196 = -qJD(3) * t330 - qJD(4) * t332 - t281 * t325 - t282 * t324;
t83 = qJD(4) * t153 + t281 * t214 - t391 * t297;
t368 = t235 * t83;
t424 = -t152 * t196 - t153 * t197 - t236 * t82 - t368;
t369 = t228 * Ifges(5,4);
t166 = t226 * Ifges(5,2) + t369 + t451;
t306 = Ifges(6,5) * t284 - Ifges(6,6) * t280;
t377 = Ifges(6,4) * t284;
t307 = -Ifges(6,2) * t280 + t377;
t378 = Ifges(6,4) * t280;
t308 = Ifges(6,1) * t284 - t378;
t321 = t343 / 0.2e1;
t379 = Ifges(6,4) * t203;
t110 = Ifges(6,2) * t202 + Ifges(6,6) * t221 + t379;
t352 = t280 * t110;
t327 = -t352 / 0.2e1;
t381 = mrSges(5,3) * t226;
t394 = t228 / 0.2e1;
t396 = t226 / 0.2e1;
t398 = -t221 / 0.2e1;
t399 = t216 / 0.2e1;
t402 = -t203 / 0.2e1;
t403 = -t202 / 0.2e1;
t407 = t135 / 0.2e1;
t409 = t318 / 0.2e1;
t410 = -t318 / 0.2e1;
t411 = t121 / 0.2e1;
t412 = t120 / 0.2e1;
t60 = Ifges(7,1) * t135 + Ifges(7,5) * t216 + t131;
t414 = t60 / 0.2e1;
t415 = -t60 / 0.2e1;
t59 = Ifges(7,2) * t318 + Ifges(7,6) * t216 + t376;
t416 = t59 / 0.2e1;
t417 = -t59 / 0.2e1;
t420 = Ifges(7,1) * t419 + Ifges(7,4) * t418 + Ifges(7,5) * t406;
t421 = Ifges(7,4) * t419 + Ifges(7,2) * t418 + Ifges(7,6) * t406;
t45 = t120 * Ifges(6,4) + t121 * Ifges(6,2) + t185 * Ifges(6,6);
t46 = t120 * Ifges(6,1) + t121 * Ifges(6,4) + t185 * Ifges(6,5);
t51 = -t121 * pkin(5) + t83;
t423 = (-t449 * mrSges(7,1) + t454 * mrSges(7,2)) * t107 + (-t454 * t19 + t449 * t20 - t237 * t4 - t300 * t3) * mrSges(7,3) + t111 * t321 + (t327 + t430) * qJD(5) + (-mrSges(6,1) * t284 + mrSges(6,2) * t280 - mrSges(5,1)) * t83 + (-Ifges(7,1) * t159 - Ifges(7,4) * t158) * t408 + (-Ifges(7,5) * t159 - Ifges(7,6) * t158) * t400 + (-Ifges(7,4) * t159 - Ifges(7,2) * t158) * t410 + (-Ifges(7,5) * t194 - Ifges(7,6) * t195) * t399 + (-Ifges(7,1) * t194 - Ifges(7,4) * t195) * t407 + (-Ifges(7,4) * t194 - Ifges(7,2) * t195) * t409 + (t202 * t307 + t203 * t308 + t221 * t306) * qJD(5) / 0.2e1 + (Ifges(6,5) * t280 + Ifges(7,5) * t237 + Ifges(6,6) * t284 - Ifges(7,6) * t300) * t406 + t51 * (mrSges(7,1) * t300 + mrSges(7,2) * t237) + (Ifges(7,4) * t237 - Ifges(7,2) * t300) * t418 + (Ifges(7,1) * t237 - Ifges(7,4) * t300) * t419 - t300 * t421 + t284 * t45 / 0.2e1 + t280 * t46 / 0.2e1 - Ifges(5,6) * t185 - Ifges(5,5) * t184 - t82 * mrSges(5,2) + (t356 * t77 + t357 * t78 + t446) * mrSges(6,3) - (-Ifges(5,2) * t228 + t219 + t447) * t226 / 0.2e1 - (Ifges(5,1) * t226 - t369 + t450) * t228 / 0.2e1 + (t20 * mrSges(7,2) + Ifges(7,3) * t400 + Ifges(7,5) * t408 + Ifges(7,6) * t410 - t19 * mrSges(7,1) + t78 * mrSges(6,2) - t77 * mrSges(6,1) - t249 * mrSges(5,1) + t451 / 0.2e1 + Ifges(6,3) * t398 + Ifges(6,5) * t402 + Ifges(6,6) * t403) * t228 + (-t430 - t249 * mrSges(5,2) - t452 / 0.2e1 + t306 * t398 + t308 * t402 + t307 * t403) * t226 + t152 * t381 + t166 * t394 + t352 * t396 + (Ifges(6,2) * t284 + t378) * t411 + (Ifges(6,1) * t280 + t377) * t412 - t194 * t414 - t159 * t415 - t195 * t416 - t158 * t417 + t237 * t420;
t422 = qJD(1) ^ 2;
t401 = t203 / 0.2e1;
t388 = t284 * pkin(5);
t387 = t77 * mrSges(6,3);
t386 = t78 * mrSges(6,3);
t380 = Ifges(4,4) * t282;
t373 = t185 * mrSges(5,3);
t371 = t429 * t83;
t73 = mrSges(6,1) * t185 - mrSges(6,3) * t120;
t365 = t280 * t73;
t359 = t196 * t284;
t355 = t235 * t280;
t354 = t235 * t284;
t248 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t326;
t350 = t282 * t248;
t255 = pkin(3) * t325 + qJD(2);
t339 = t391 * pkin(3);
t335 = Ifges(4,4) * t392;
t334 = Ifges(6,5) * t120 + Ifges(6,6) * t121 + Ifges(6,3) * t185;
t114 = pkin(4) * t197 - pkin(9) * t196 + t255;
t232 = t384 * t346;
t233 = t245 * qJD(3);
t124 = qJD(4) * t429 + t281 * t232 - t391 * t233;
t319 = t284 * t114 - t124 * t280;
t115 = t284 * t191 - t200 * t280;
t269 = -t339 - pkin(4);
t79 = pkin(5) * t236 + pkin(10) * t354 + t115;
t96 = pkin(10) * t355 + t116;
t38 = -t279 * t96 + t283 * t79;
t39 = t279 * t79 + t283 * t96;
t74 = -mrSges(6,2) * t185 + mrSges(6,3) * t121;
t305 = t284 * t74 - t365;
t304 = t280 * t78 + t284 * t77;
t147 = -mrSges(6,2) * t221 + mrSges(6,3) * t202;
t148 = mrSges(6,1) * t221 - mrSges(6,3) * t203;
t301 = -t280 * t147 - t284 * t148;
t298 = -Ifges(4,5) * t282 - Ifges(4,6) * t392;
t296 = -t196 * t280 + t235 * t343;
t295 = t235 * t344 + t359;
t32 = t280 * t114 + t284 * t124 + t191 * t343 - t200 * t344;
t294 = qJ(2) * (mrSges(4,1) * t392 - mrSges(4,2) * t282);
t293 = t282 * (-Ifges(4,2) * t392 - t380);
t292 = (Ifges(4,1) * t392 - t380) * qJD(1);
t291 = (-Ifges(4,2) * t282 + t335) * qJD(1);
t238 = (t282 * mrSges(4,1) + mrSges(4,2) * t392) * qJD(1);
t290 = (-Ifges(4,1) * t282 - t335) * t392;
t288 = -qJD(5) * t304 + t446;
t125 = qJD(4) * t200 - t391 * t232 - t281 * t233;
t270 = -pkin(4) - t388;
t247 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t347;
t246 = t269 - t388;
t224 = Ifges(4,5) * qJD(3) + t292;
t223 = Ifges(4,6) * qJD(3) + t291;
t209 = -mrSges(5,2) * t276 + t381;
t189 = -mrSges(5,1) * t226 + mrSges(5,2) * t228;
t175 = t300 * t235;
t172 = t237 * t236;
t154 = -pkin(5) * t355 - t429;
t123 = t215 + t153;
t102 = mrSges(7,1) * t216 - mrSges(7,3) * t135;
t101 = -mrSges(7,2) * t216 + mrSges(7,3) * t318;
t71 = -pkin(5) * t296 + t125;
t70 = -mrSges(7,1) * t318 + mrSges(7,2) * t135;
t56 = -t194 * t235 - t196 * t237;
t54 = t173 * t425 - t300 * t196;
t33 = -qJD(5) * t116 + t319;
t25 = -mrSges(7,2) * t185 + mrSges(7,3) * t37;
t24 = mrSges(7,1) * t185 - mrSges(7,3) * t36;
t23 = t283 * t63 - t366;
t22 = -t279 * t63 - t363;
t21 = pkin(10) * t296 + t32;
t15 = -pkin(10) * t359 + pkin(5) * t197 + (-t192 + (-pkin(10) * t235 - t191) * t280) * qJD(5) + t319;
t11 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t6 = -qJD(6) * t39 + t15 * t283 - t21 * t279;
t5 = qJD(6) * t38 + t15 * t279 + t21 * t283;
t1 = [-t431 * t429 + t196 * t327 + t424 * mrSges(5,3) + t235 * t110 * t321 + m(7) * (t107 * t71 + t154 * t51 + t19 * t6 + t20 * t5 + t3 * t39 + t38 * t4) + (0.2e1 * t294 - t293 + t290) * t342 + m(6) * (t115 * t18 + t116 * t17 + t32 * t78 + t33 * t77 - t371) + m(5) * (t124 * t153 + t200 * t82 + t243 * t266 + t249 * t255 - t371) + (-m(5) * t152 + m(6) * t143 - t348) * t125 + (t184 * t236 + t185 * t235 + t196 * t396 - t197 * t394) * Ifges(5,4) + (t184 * t235 + t196 * t394) * Ifges(5,1) + t266 * (mrSges(5,1) * t185 - mrSges(5,2) * t184) - (t291 + t223) * t325 / 0.2e1 - (t292 + t224) * t346 / 0.2e1 + (qJD(5) * t111 + t45) * t355 / 0.2e1 + (t340 + t334) * t236 / 0.2e1 + (Ifges(7,5) * t175 + Ifges(7,6) * t173 - t235 * t306 + (Ifges(6,3) + Ifges(7,3)) * t236) * t406 - t200 * t373 - t309 * t368 + t276 * (Ifges(5,5) * t196 - Ifges(5,6) * t197) / 0.2e1 + t255 * t189 + t243 * (mrSges(5,1) * t236 - mrSges(5,2) * t235) + t249 * (mrSges(5,1) * t197 + mrSges(5,2) * t196) + t4 * (mrSges(7,1) * t236 - mrSges(7,3) * t175) + t3 * (-mrSges(7,2) * t236 + mrSges(7,3) * t173) + t17 * (-mrSges(6,2) * t236 + mrSges(6,3) * t355) + 0.2e1 * qJD(2) * t238 - t46 * t354 / 0.2e1 + t18 * (mrSges(6,1) * t236 + mrSges(6,3) * t354) + t124 * t209 - t197 * t166 / 0.2e1 + t20 * (-mrSges(7,2) * t197 + mrSges(7,3) * t56) + t19 * (mrSges(7,1) * t197 - mrSges(7,3) * t54) + t51 * (-mrSges(7,1) * t173 + mrSges(7,2) * t175) + t154 * t11 + t32 * t147 + t33 * t148 + t115 * t73 + t116 * t74 + t107 * (-mrSges(7,1) * t56 + mrSges(7,2) * t54) + t5 * t101 + (t247 * t325 - t248 * t346) * t285 + t6 * t102 + t71 * t70 + t38 * t24 + t39 * t25 + 0.2e1 * t445 * t277 + t447 * t196 / 0.2e1 + t450 * t197 / 0.2e1 + t77 * (mrSges(6,1) * t197 - mrSges(6,3) * t295) + t78 * (-mrSges(6,2) * t197 + mrSges(6,3) * t296) + t221 * (Ifges(6,5) * t295 + Ifges(6,6) * t296 + Ifges(6,3) * t197) / 0.2e1 + t202 * (Ifges(6,4) * t295 + Ifges(6,2) * t296 + Ifges(6,6) * t197) / 0.2e1 + t143 * (-mrSges(6,1) * t296 + mrSges(6,2) * t295) + qJD(3) ^ 2 * t298 / 0.2e1 + (Ifges(7,5) * t54 + Ifges(7,6) * t56 + Ifges(7,3) * t197) * t399 + (Ifges(6,1) * t295 + Ifges(6,4) * t296 + Ifges(6,5) * t197) * t401 + (Ifges(7,1) * t54 + Ifges(7,4) * t56 + Ifges(7,5) * t197) * t407 + (Ifges(7,4) * t54 + Ifges(7,2) * t56 + Ifges(7,6) * t197) * t409 + (Ifges(6,6) * t236 - t235 * t307) * t411 + (Ifges(6,5) * t236 - t235 * t308) * t412 + t54 * t414 + t56 * t416 + (Ifges(7,4) * t175 + Ifges(7,2) * t173 + Ifges(7,6) * t236) * t418 + (Ifges(7,1) * t175 + Ifges(7,4) * t173 + Ifges(7,5) * t236) * t419 + t175 * t420 + t173 * t421 + (t185 * t236 - t197 * t396) * Ifges(5,2); -t172 * t24 - t174 * t25 + t433 * t102 + t432 * t101 + (t247 * t392 - t350) * qJD(3) - t445 * t422 + (t11 + t431) * t235 + (t147 * t284 - t148 * t280 + t209) * t197 + (-t70 + t348) * t196 + (qJD(5) * t301 + t305 - t373) * t236 + m(6) * (-t143 * t196 + t197 * t426 + t288 * t236 + t368) - m(5) * t424 + (-m(5) * t249 - m(6) * t304 - t189 - t238 + t301) * qJD(1) + (-t107 * t196 - t172 * t4 - t174 * t3 + t19 * t433 + t20 * t432 + t235 * t51) * m(7); (-t268 * t343 + t459) * t148 + (-t268 * t344 + t460) * t147 + (-mrSges(4,1) * t346 - mrSges(4,2) * t325 + t350) * t254 + t423 + (-t290 / 0.2e1 + t293 / 0.2e1 - t294) * t422 + (t83 * t269 + (t143 * t281 + t391 * t426) * qJD(4) * pkin(3) + t288 * t268 - t143 * t160 - t77 * t90 - t78 * t91) * m(6) + t428 * t70 - t373 * t390 + (t316 - t161) * t209 - t344 * t386 - t343 * t387 + t339 * t382 - t268 * t365 + t269 * t62 + t246 * t11 + t224 * t347 / 0.2e1 - t298 * t342 / 0.2e1 + t186 * t24 + t187 * t25 - t247 * t242 + (t152 * t160 - t153 * t161 - t249 * t315 + (-t391 * t83 + t281 * t82 + (-t152 * t281 + t153 * t391) * qJD(4)) * pkin(3)) * m(5) + t223 * t326 / 0.2e1 - Ifges(4,5) * t323 - t189 * t315 - Ifges(4,6) * t311 + t440 * t102 + t441 * t101 + (t107 * t428 + t186 * t4 + t187 * t3 + t19 * t440 + t20 * t441 + t246 * t51) * m(7) - t348 * t448 + t74 * t353 + t153 * t370; t434 * t102 + t435 * t101 + t270 * t11 + t204 * t24 + t205 * t25 - t152 * t209 - t95 * t147 - t94 * t148 - t123 * t70 + ((-pkin(9) * t148 - t387) * t284 + (pkin(5) * t70 - pkin(9) * t147 - t386) * t280) * qJD(5) - pkin(4) * t62 + (t348 + t370) * t153 + t305 * pkin(9) + (t204 * t4 + t205 * t3 + t270 * t51 + t435 * t20 + t434 * t19 + (-t123 + t337) * t107) * m(7) + (-pkin(4) * t83 + pkin(9) * t288 - t143 * t153 - t77 * t94 - t78 * t95) * m(6) + t423; (t202 * t77 + t203 * t78) * mrSges(6,3) + (-t203 * t70 + t283 * t24 + t279 * t25 + (t101 * t283 - t102 * t279) * qJD(6) + (-t107 * t203 + t279 * t3 + t283 * t4 + (-t19 * t279 + t20 * t283) * qJD(6)) * m(7)) * pkin(5) + t318 * t415 + t334 - t135 * t417 - m(7) * (t19 * t22 + t20 * t23) + (-Ifges(6,2) * t203 + t111 + t201) * t403 - t143 * (mrSges(6,1) * t203 + mrSges(6,2) * t202) - t77 * t147 + t78 * t148 - t23 * t101 - t22 * t102 + t18 * mrSges(6,1) - t17 * mrSges(6,2) + (Ifges(6,5) * t202 - Ifges(6,6) * t203) * t398 + t110 * t401 + (Ifges(6,1) * t202 - t379) * t402 + t444 * t410 + t453; t59 * t407 - t19 * t101 + t20 * t102 + (t444 + t60) * t410 + t453;];
tauc  = t1(:);
