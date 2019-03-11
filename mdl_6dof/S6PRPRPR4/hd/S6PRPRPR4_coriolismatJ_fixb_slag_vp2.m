% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRPR4
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:22
% EndTime: 2019-03-08 19:39:31
% DurationCPUTime: 5.82s
% Computational Cost: add. (15398->479), mult. (35017->701), div. (0->0), fcn. (41029->12), ass. (0->255)
t283 = cos(pkin(11));
t388 = pkin(8) + qJ(3);
t266 = t388 * t283;
t280 = sin(pkin(11));
t328 = t388 * t280;
t390 = sin(qJ(4));
t392 = cos(qJ(4));
t235 = t266 * t392 - t328 * t390;
t259 = t280 * t390 - t283 * t392;
t279 = sin(pkin(12));
t355 = t259 * t279;
t171 = -pkin(5) * t355 + t235;
t419 = -m(7) / 0.2e1;
t421 = -m(6) / 0.2e1;
t282 = cos(pkin(12));
t284 = sin(qJ(6));
t285 = cos(qJ(6));
t316 = t279 * t284 - t282 * t285;
t197 = t316 * t259;
t378 = t197 * mrSges(7,2);
t262 = t279 * t285 + t282 * t284;
t194 = t262 * t259;
t380 = t194 * mrSges(7,1);
t431 = t380 / 0.2e1 - t378 / 0.2e1;
t445 = t171 * t419 + t235 * t421 + t431;
t261 = -t280 * t392 - t283 * t390;
t432 = t316 * t261;
t443 = -t432 / 0.2e1;
t444 = mrSges(7,1) * t443;
t304 = t262 * t261;
t434 = -t304 / 0.2e1;
t442 = mrSges(7,2) * t434;
t340 = t280 ^ 2 + t283 ^ 2;
t441 = t340 * mrSges(4,3);
t397 = -t262 / 0.2e1;
t440 = t397 * t432;
t281 = sin(pkin(6));
t286 = cos(qJ(2));
t349 = t281 * t286;
t240 = t259 * t349;
t391 = sin(qJ(2));
t332 = t281 * t391;
t207 = t279 * t240 + t282 * t332;
t208 = -t282 * t240 + t279 * t332;
t121 = t207 * t285 - t208 * t284;
t122 = t207 * t284 + t208 * t285;
t387 = pkin(9) + qJ(5);
t263 = t387 * t279;
t265 = t387 * t282;
t232 = -t263 * t285 - t265 * t284;
t234 = -t263 * t284 + t265 * t285;
t239 = t261 * t349;
t272 = -pkin(5) * t282 - pkin(4);
t359 = t208 * t282;
t360 = t207 * t279;
t439 = -(t359 / 0.2e1 - t360 / 0.2e1) * mrSges(6,3) + (pkin(4) * t239 + (t359 - t360) * qJ(5)) * t421 + (t121 * t232 + t122 * t234 - t239 * t272) * t419 - t240 * mrSges(5,2) / 0.2e1;
t438 = -(t282 * t207 + t279 * t208) * t421 - (-t121 * t316 + t262 * t122) * t419 + (m(4) + m(5)) * t332 / 0.2e1;
t338 = -m(7) / 0.4e1 - m(6) / 0.4e1;
t437 = 0.2e1 * t338;
t436 = mrSges(5,1) / 0.2e1;
t415 = -mrSges(7,3) / 0.2e1;
t274 = t279 ^ 2;
t277 = t282 ^ 2;
t341 = t274 + t277;
t433 = t341 * mrSges(6,3);
t225 = mrSges(7,1) * t316 + mrSges(7,2) * t262;
t264 = -mrSges(6,1) * t282 + mrSges(6,2) * t279;
t429 = t264 + t225;
t366 = cos(pkin(6));
t249 = t280 * t366 + t283 * t332;
t302 = t280 * t332 - t283 * t366;
t428 = t283 * t249 + t280 * t302;
t389 = pkin(4) * t261;
t223 = qJ(5) * t259 - t389;
t233 = t266 * t390 + t392 * t328;
t142 = t282 * t223 + t233 * t279;
t143 = t279 * t223 - t282 * t233;
t427 = -t142 * t279 + t143 * t282;
t118 = -mrSges(7,1) * t304 + mrSges(7,2) * t432;
t371 = t282 * mrSges(6,2);
t374 = t279 * mrSges(6,1);
t323 = -t371 - t374;
t210 = t323 * t261;
t426 = t261 * mrSges(5,3) - t118 - t210;
t147 = -mrSges(7,2) * t259 + mrSges(7,3) * t304;
t149 = mrSges(7,1) * t259 - mrSges(7,3) * t432;
t400 = -t316 / 0.2e1;
t310 = t147 * t400 + t149 * t397;
t294 = (-t316 * t434 + t440) * mrSges(7,3) + t310;
t17 = t294 - t431;
t424 = t17 * qJD(2);
t423 = 0.2e1 * t261;
t422 = 2 * qJD(4);
t420 = m(6) / 0.2e1;
t418 = m(7) / 0.2e1;
t417 = -mrSges(7,1) / 0.2e1;
t416 = mrSges(7,2) / 0.2e1;
t191 = t249 * t390 + t302 * t392;
t104 = t262 * t191;
t105 = t316 * t191;
t47 = -t104 * t316 + t105 * t262;
t414 = m(7) * t47;
t413 = t147 / 0.2e1;
t412 = t149 / 0.2e1;
t409 = t194 / 0.2e1;
t406 = t197 / 0.2e1;
t222 = t262 * mrSges(7,1) - mrSges(7,2) * t316;
t404 = t222 / 0.2e1;
t255 = Ifges(7,4) * t316;
t229 = Ifges(7,1) * t262 - t255;
t403 = t229 / 0.2e1;
t402 = -t234 / 0.2e1;
t399 = -t261 / 0.2e1;
t398 = t261 / 0.2e1;
t396 = t262 / 0.2e1;
t395 = t279 / 0.2e1;
t394 = -t282 / 0.2e1;
t393 = t282 / 0.2e1;
t386 = qJD(4) * t414 / 0.2e1;
t385 = m(7) * qJD(3);
t384 = Ifges(6,4) * t279;
t383 = Ifges(6,4) * t282;
t382 = Ifges(7,4) * t432;
t381 = Ifges(7,4) * t262;
t379 = t304 * mrSges(7,2);
t377 = t432 * mrSges(7,1);
t376 = t259 * mrSges(5,3);
t373 = t279 * Ifges(6,2);
t372 = t279 * Ifges(6,6);
t370 = t282 * Ifges(6,5);
t192 = t249 * t392 - t302 * t390;
t162 = -t279 * t192 - t282 * t349;
t163 = t282 * t192 - t279 * t349;
t86 = t162 * t285 - t163 * t284;
t369 = t86 * t304;
t87 = t162 * t284 + t163 * t285;
t368 = t87 * t432;
t365 = t121 * t262;
t364 = t122 * t316;
t144 = t191 * t239;
t358 = t233 * t239;
t357 = t233 * t261;
t333 = t281 ^ 2 * t391;
t24 = m(6) * (t162 * t207 + t163 * t208 - t144) + m(7) * (t121 * t86 + t122 * t87 - t144) + m(5) * (-t192 * t240 - t286 * t333 - t144) + m(4) * (t281 * t428 - t333) * t286;
t356 = t24 * qJD(1);
t354 = t259 * t282;
t161 = t261 * t191;
t353 = t261 * t279;
t352 = t261 * t282;
t218 = t259 * mrSges(6,1) + mrSges(6,3) * t352;
t351 = t279 * t218;
t216 = -t259 * mrSges(6,2) + mrSges(6,3) * t353;
t348 = t282 * t216;
t307 = m(7) * (t262 * t304 + t316 * t432);
t329 = m(6) * t341;
t58 = -t307 / 0.2e1 + (-t329 / 0.4e1 + t338) * t423;
t345 = t58 * qJD(2);
t343 = Ifges(7,5) * t304 - Ifges(7,6) * t432;
t273 = -pkin(3) * t283 - pkin(2);
t212 = pkin(4) * t259 + qJ(5) * t261 + t273;
t139 = t279 * t212 + t282 * t235;
t342 = -Ifges(7,5) * t316 - Ifges(7,6) * t262;
t339 = t222 * qJD(6);
t331 = t191 * t404;
t330 = -t225 / 0.2e1 - t264 / 0.2e1;
t327 = qJ(5) * t341;
t326 = t340 * qJ(3);
t138 = t282 * t212 - t235 * t279;
t324 = t239 * t437;
t100 = Ifges(7,1) * t197 + Ifges(7,4) * t194 - Ifges(7,5) * t261;
t190 = Ifges(7,4) * t304;
t101 = Ifges(7,1) * t432 + t259 * Ifges(7,5) + t190;
t117 = t378 - t380;
t146 = mrSges(7,2) * t261 + t194 * mrSges(7,3);
t148 = -mrSges(7,1) * t261 - t197 * mrSges(7,3);
t164 = -Ifges(6,6) * t261 + (t373 - t383) * t259;
t165 = -Ifges(6,5) * t261 + (-t282 * Ifges(6,1) + t384) * t259;
t170 = -pkin(5) * t353 + t233;
t209 = t323 * t259;
t215 = mrSges(6,2) * t261 + mrSges(6,3) * t355;
t217 = -mrSges(6,1) * t261 + mrSges(6,3) * t354;
t251 = t259 * mrSges(5,2);
t224 = -t261 * mrSges(5,1) - t251;
t313 = Ifges(7,5) * t406 + Ifges(7,6) * t409;
t109 = pkin(9) * t353 + t139;
t92 = pkin(5) * t259 + pkin(9) * t352 + t138;
t51 = -t109 * t284 + t285 * t92;
t52 = t109 * t285 + t284 * t92;
t115 = pkin(9) * t355 + t143;
t97 = -pkin(5) * t261 + pkin(9) * t354 + t142;
t56 = -t115 * t284 + t285 * t97;
t57 = t115 * t285 + t284 * t97;
t98 = Ifges(7,4) * t197 + Ifges(7,2) * t194 - Ifges(7,6) * t261;
t99 = Ifges(7,2) * t304 + t259 * Ifges(7,6) + t382;
t1 = t273 * t224 + t235 * t210 + t233 * t209 + t139 * t215 + t143 * t216 + t138 * t217 + t142 * t218 + t101 * t406 + t432 * t100 / 0.2e1 + t99 * t409 + t304 * t98 / 0.2e1 + t170 * t117 + t171 * t118 + t56 * t149 + t52 * t146 + t57 * t147 + t51 * t148 + m(6) * (t138 * t142 + t139 * t143 + t233 * t235) + m(7) * (t170 * t171 + t51 * t56 + t52 * t57) + (Ifges(7,5) * t443 + Ifges(7,6) * t434 + t165 * t394 + t164 * t395 + (t370 / 0.2e1 - t372 / 0.2e1 - Ifges(5,4)) * t261) * t261 + ((-Ifges(6,3) + Ifges(5,1) - Ifges(5,2) - Ifges(7,3) + t277 * Ifges(6,1) / 0.2e1 + (-t383 + t373 / 0.2e1) * t279) * t261 + t313 + (Ifges(5,4) - t370 + t372) * t259) * t259;
t309 = t351 / 0.2e1 - t348 / 0.2e1;
t318 = t138 * t279 - t139 * t282;
t287 = (t210 / 0.2e1 + t118 / 0.2e1) * t192 + (t209 / 0.2e1 + t117 / 0.2e1 + t309) * t191 + (t142 * t162 + t143 * t163 + t192 * t233 + (t235 + t318) * t191) * t420 + (t104 * t51 + t105 * t52 + t170 * t192 + t171 * t191 + t56 * t86 + t57 * t87) * t418 + t104 * t412 + t105 * t413 + t162 * t217 / 0.2e1 + t163 * t215 / 0.2e1 + t86 * t148 / 0.2e1 + t87 * t146 / 0.2e1 - t224 * t349 / 0.2e1;
t3 = t287 - (t436 + t330) * t239 + (t365 / 0.2e1 + t364 / 0.2e1) * mrSges(7,3) + t439;
t322 = t3 * qJD(1) + t1 * qJD(2);
t116 = t377 + t379;
t119 = -Ifges(7,2) * t432 + t190;
t120 = Ifges(7,1) * t304 - t382;
t6 = t51 * t147 - t52 * t149 + t170 * t116 + t259 * t343 / 0.2e1 + (-t52 * mrSges(7,3) + t120 / 0.2e1 - t99 / 0.2e1) * t432 + (-t51 * mrSges(7,3) + t101 / 0.2e1 + t119 / 0.2e1) * t304;
t301 = -t191 * t116 / 0.2e1 - t86 * t147 / 0.2e1 + t87 * t412;
t311 = t121 * mrSges(7,1) / 0.2e1 - t122 * mrSges(7,2) / 0.2e1;
t7 = (t369 / 0.2e1 + t368 / 0.2e1) * mrSges(7,3) + t301 + t311;
t321 = -t7 * qJD(1) + t6 * qJD(2);
t10 = t197 * t147 + t194 * t149 + t441 + t426 * t261 + (-t348 + t351 + t376) * t259 + m(7) * (-t170 * t261 + t194 * t51 + t197 * t52) + m(6) * (t259 * t318 - t357) + m(5) * (-t235 * t259 - t357) + m(4) * t326;
t317 = t162 * t279 - t163 * t282;
t289 = (t259 * t317 - t161) * t420 + (t194 * t86 + t197 * t87 - t161) * t418 + m(5) * (-t192 * t259 - t161) / 0.2e1 + m(4) * t428 / 0.2e1;
t21 = t289 - t438;
t320 = t21 * qJD(1) + t10 * qJD(2);
t22 = m(7) * (t304 * t52 - t432 * t51) + t304 * t147 - t432 * t149 + (t279 * t216 + t282 * t218 + m(6) * (t138 * t282 + t139 * t279)) * t261;
t295 = (t304 * t87 - t432 * t86) * t419 + m(6) * (t162 * t282 + t163 * t279) * t399;
t31 = t324 + t295;
t319 = -qJD(1) * t31 + qJD(2) * t22;
t65 = 0.2e1 * t442 + 0.2e1 * t444;
t315 = qJD(2) * t65 - qJD(4) * t222;
t312 = t104 * t417 + t105 * t416;
t291 = (t194 * t397 + t197 * t400) * mrSges(7,3) + (-t277 / 0.2e1 - t274 / 0.2e1) * t259 * mrSges(6,3) + (-t259 * t327 + t389) * t420 + (t194 * t232 + t197 * t234 - t261 * t272) * t418;
t292 = (t142 * t282 + t143 * t279) * t420 + (t262 * t57 - t316 * t56) * t418 + t148 * t400 + t146 * t396 + t215 * t395 + t217 * t393;
t14 = t251 + (mrSges(5,1) + t330) * t261 + t291 - t292;
t306 = t14 * qJD(2) - qJD(1) * t414 / 0.2e1;
t123 = t191 * t192;
t23 = m(6) * (t191 * t317 + t123) + m(7) * (t104 * t86 + t105 * t87 + t123);
t305 = -t23 * qJD(1) - t385 * t47 / 0.2e1;
t290 = (t304 * t400 - t440) * mrSges(7,3) + t318 * t421 + (-t232 * t432 + t234 * t304 - t262 * t51 - t316 * t52) * t418 - t309 + t310;
t13 = (t371 / 0.2e1 + t374 / 0.2e1) * t259 + t290 + t445;
t298 = t317 * t421 + (-t262 * t86 - t316 * t87) * t418;
t35 = t192 * t437 + t298;
t72 = (t262 ^ 2 + t316 ^ 2) * mrSges(7,3) + t433 + m(7) * (-t232 * t262 - t234 * t316) + qJ(5) * t329;
t303 = t35 * qJD(1) + t13 * qJD(2) + t72 * qJD(4);
t19 = t331 + t312;
t226 = -Ifges(7,2) * t262 - t255;
t227 = -Ifges(7,2) * t316 + t381;
t228 = -Ifges(7,1) * t316 - t381;
t37 = t272 * t222 + (t228 / 0.2e1 - t227 / 0.2e1) * t262 - (t403 + t226 / 0.2e1) * t316;
t288 = -(t101 / 0.4e1 + t119 / 0.4e1) * t316 + (t120 / 0.4e1 - t99 / 0.4e1) * t262 + (t232 * t415 + t229 / 0.4e1 + t226 / 0.4e1) * t304 + (mrSges(7,3) * t402 + t228 / 0.4e1 - t227 / 0.4e1) * t432 + t170 * t404 + t232 * t413 + t149 * t402 + t259 * t342 / 0.4e1 + t272 * t116 / 0.2e1;
t297 = Ifges(7,3) * t398 + t416 * t57 + t417 * t56 - t313;
t5 = t288 + t297;
t299 = -t19 * qJD(1) - t5 * qJD(2) - t37 * qJD(4);
t66 = t379 / 0.2e1 + t377 / 0.2e1 + t444 + t442;
t59 = t307 / 0.2e1 + t329 * t398 + t338 * t423;
t34 = t298 + (m(6) + m(7)) * t192 / 0.2e1;
t32 = t324 - t295;
t20 = t289 + t438;
t18 = t331 - t312;
t16 = t294 + t431;
t15 = t261 * t330 + t291 + t292;
t12 = -mrSges(6,2) * t354 / 0.2e1 - mrSges(6,1) * t355 / 0.2e1 + t290 - t445;
t8 = -t301 + t311 + (t368 + t369) * t415;
t4 = t288 - t297;
t2 = t287 + (t365 + t364) * t415 + (t436 - t429 / 0.2e1) * t239 - t439;
t9 = [t24 * qJD(2) + t23 * qJD(4), t20 * qJD(3) + t2 * qJD(4) + t32 * qJD(5) + t8 * qJD(6) + t356 + (t121 * t149 + m(6) * (t138 * t207 + t139 * t208 - t358) + m(7) * (t121 * t51 + t122 * t52) + t122 * t147 + t207 * t218 + m(5) * (-t235 * t240 - t358) + t240 * t376 + t208 * t216 + m(4) * (-pkin(2) * t391 + t286 * t326) * t281 + (m(5) * t273 - mrSges(4,1) * t283 + mrSges(5,1) * t259 + mrSges(4,2) * t280 - mrSges(5,2) * t261 - mrSges(3,1)) * t332 - (m(7) * t170 - t426) * t239 + (-mrSges(3,2) + t441) * t349) * qJD(2), qJD(2) * t20 + t386, t2 * qJD(2) + t34 * qJD(5) + t18 * qJD(6) + ((-pkin(4) * t192 - t191 * t327) * t420 + (t104 * t232 + t105 * t234 + t192 * t272) * t418) * t422 - t305 + ((-t104 * t262 - t105 * t316) * mrSges(7,3) + (-mrSges(5,1) + t429) * t192 + (mrSges(5,2) - t433) * t191) * qJD(4), qJD(2) * t32 + qJD(4) * t34, t8 * qJD(2) + t18 * qJD(4) + (-mrSges(7,1) * t87 - mrSges(7,2) * t86) * qJD(6); qJD(3) * t21 + qJD(4) * t3 - qJD(5) * t31 - qJD(6) * t7 - t356, qJD(3) * t10 + qJD(4) * t1 + qJD(5) * t22 + qJD(6) * t6 (-t194 * t316 + t197 * t262) * t385 + t15 * qJD(4) + t59 * qJD(5) + t16 * qJD(6) + t320, t15 * qJD(3) + t12 * qJD(5) + t4 * qJD(6) + ((-pkin(4) * t235 + qJ(5) * t427) * t420 + (t171 * t272 + t232 * t56 + t234 * t57) * t418) * t422 + t322 + (t164 * t393 + t165 * t395 + t272 * t117 + t100 * t396 + Ifges(5,6) * t261 + t98 * t400 + t232 * t148 + t233 * mrSges(5,2) + t234 * t146 + t171 * t225 + t227 * t409 + t197 * t403 - pkin(4) * t209 + (-Ifges(5,5) + (Ifges(6,1) * t279 + t383) * t394 + (Ifges(6,2) * t282 + t384) * t395) * t259 + (Ifges(6,5) * t279 + Ifges(7,5) * t262 + Ifges(6,6) * t282 - Ifges(7,6) * t316) * t399 + (t264 - mrSges(5,1)) * t235 + (t282 * t215 - t279 * t217) * qJ(5) + (-t262 * t56 - t316 * t57) * mrSges(7,3) + t427 * mrSges(6,3)) * qJD(4), qJD(3) * t59 + qJD(4) * t12 + qJD(6) * t66 + t319, t16 * qJD(3) + t4 * qJD(4) + t66 * qJD(5) + (-mrSges(7,1) * t52 - mrSges(7,2) * t51 + t343) * qJD(6) + t321; -qJD(2) * t21 + t386, -qJD(4) * t14 - qJD(5) * t58 + qJD(6) * t17 - t320, 0, -t306, -t345, -t339 + t424; -t3 * qJD(2) + t35 * qJD(5) + t19 * qJD(6) + t305, qJD(3) * t14 + qJD(5) * t13 + qJD(6) * t5 - t322, t306, qJD(5) * t72 + qJD(6) * t37, t303 (-mrSges(7,1) * t234 - mrSges(7,2) * t232 + t342) * qJD(6) - t299; qJD(2) * t31 - qJD(4) * t35, qJD(3) * t58 - qJD(4) * t13 - qJD(6) * t65 - t319, t345, -t303 + t339, 0, -t315; t7 * qJD(2) - t19 * qJD(4), -qJD(3) * t17 - qJD(4) * t5 + qJD(5) * t65 - t321, -t424, -t222 * qJD(5) + t299, t315, 0;];
Cq  = t9;
