% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPPR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:21
% EndTime: 2019-03-09 08:17:27
% DurationCPUTime: 3.92s
% Computational Cost: add. (7898->475), mult. (15310->640), div. (0->0), fcn. (14254->6), ass. (0->240)
t409 = pkin(3) + pkin(7);
t426 = Ifges(5,1) + Ifges(6,1);
t425 = Ifges(6,4) + Ifges(5,5);
t278 = cos(qJ(2));
t237 = t409 * t278;
t272 = cos(pkin(9));
t271 = sin(pkin(9));
t366 = qJ(5) * t271;
t408 = pkin(4) + pkin(5);
t299 = t272 * t408 + t366;
t131 = -t278 * t299 - t237;
t424 = Ifges(6,6) / 0.2e1;
t276 = sin(qJ(2));
t423 = -t276 / 0.2e1;
t390 = t276 / 0.2e1;
t422 = -t278 / 0.2e1;
t421 = -Ifges(3,4) - Ifges(4,6);
t277 = cos(qJ(6));
t389 = sin(qJ(6));
t298 = t271 * t389 + t277 * t272;
t181 = t298 * t278;
t172 = Ifges(7,4) * t181;
t330 = t389 * t272;
t355 = t271 * t278;
t182 = t277 * t355 - t278 * t330;
t96 = -Ifges(7,1) * t182 - Ifges(7,5) * t276 + t172;
t420 = Ifges(7,2) * t182 + t172 + t96;
t215 = t271 * t277 - t330;
t135 = -t215 * mrSges(7,1) + mrSges(7,2) * t298;
t353 = t272 * t278;
t226 = -t276 * mrSges(5,2) - mrSges(5,3) * t353;
t227 = -mrSges(6,2) * t353 + t276 * mrSges(6,3);
t418 = t227 + t226;
t317 = t276 * pkin(2) - qJ(3) * t278;
t210 = qJ(4) * t276 + t317;
t129 = t272 * t210 + t271 * t237;
t116 = t278 * qJ(5) + t129;
t354 = t272 * t276;
t100 = -pkin(8) * t354 + t116;
t128 = -t271 * t210 + t237 * t272;
t356 = t271 * t276;
t90 = -pkin(8) * t356 - t278 * t408 - t128;
t44 = -t100 * t389 + t277 * t90;
t45 = t277 * t100 + t389 * t90;
t417 = t215 * t45 - t298 * t44;
t348 = -t272 * qJ(5) + qJ(3);
t201 = t271 * t408 + t348;
t231 = pkin(4) * t271 + t348;
t233 = t271 * mrSges(6,1) - t272 * mrSges(6,3);
t416 = -m(6) * t231 - m(7) * t201 + t135 - t233;
t270 = t272 ^ 2;
t415 = -m(5) / 0.2e1;
t414 = m(5) / 0.2e1;
t413 = -m(6) / 0.2e1;
t412 = m(6) / 0.2e1;
t411 = -m(7) / 0.2e1;
t410 = m(7) / 0.2e1;
t179 = t298 * t276;
t407 = t179 / 0.2e1;
t406 = -t179 / 0.2e1;
t352 = t276 * t277;
t180 = t271 * t352 - t276 * t330;
t405 = -t180 / 0.2e1;
t404 = t180 / 0.2e1;
t403 = -t181 / 0.2e1;
t401 = -t182 / 0.2e1;
t400 = t182 / 0.2e1;
t399 = t298 / 0.2e1;
t398 = -t298 / 0.2e1;
t397 = -t215 / 0.2e1;
t396 = t215 / 0.2e1;
t392 = -t272 / 0.2e1;
t391 = t272 / 0.2e1;
t273 = -pkin(2) - qJ(4);
t388 = pkin(8) + t273;
t387 = Ifges(5,4) * t271;
t386 = Ifges(5,4) * t272;
t385 = Ifges(7,4) * t182;
t384 = Ifges(7,4) * t298;
t383 = Ifges(6,5) * t271;
t382 = Ifges(6,5) * t272;
t379 = t180 * mrSges(7,2);
t380 = t179 * mrSges(7,1);
t101 = t379 + t380;
t376 = t182 * mrSges(7,2);
t378 = t181 * mrSges(7,1);
t102 = -t376 - t378;
t319 = -qJ(3) * t276 - pkin(1);
t209 = t273 * t278 + t319;
t236 = t409 * t276;
t127 = t272 * t209 + t271 * t236;
t113 = t276 * qJ(5) + t127;
t126 = -t271 * t209 + t236 * t272;
t114 = -pkin(4) * t276 - t126;
t118 = -pkin(4) * t278 - t128;
t130 = (t299 + t409) * t276;
t148 = mrSges(7,2) * t278 - t179 * mrSges(7,3);
t149 = -mrSges(7,1) * t278 - t180 * mrSges(7,3);
t370 = t276 * mrSges(7,2);
t377 = t181 * mrSges(7,3);
t150 = t370 + t377;
t375 = t182 * mrSges(7,3);
t151 = -mrSges(7,1) * t276 + t375;
t310 = pkin(4) * t272 + t366;
t154 = (-t310 - t409) * t276;
t155 = t278 * t310 + t237;
t313 = t272 * mrSges(6,1) + t271 * mrSges(6,3);
t202 = t313 * t276;
t203 = (-t272 * mrSges(5,1) + t271 * mrSges(5,2)) * t276;
t204 = t313 * t278;
t221 = t278 * mrSges(5,1) - mrSges(5,3) * t356;
t252 = mrSges(6,2) * t356;
t368 = t278 * mrSges(6,1);
t222 = t252 - t368;
t223 = t276 * mrSges(5,1) + mrSges(5,3) * t355;
t224 = -t276 * mrSges(6,1) - mrSges(6,2) * t355;
t225 = -t278 * mrSges(5,2) + mrSges(5,3) * t354;
t228 = mrSges(6,2) * t354 + t278 * mrSges(6,3);
t232 = -pkin(2) * t278 + t319;
t301 = Ifges(7,5) * t405 + Ifges(7,6) * t407;
t320 = m(4) * t232 + t278 * mrSges(4,2) - t276 * mrSges(4,3);
t325 = (t426 * t271 - t382 + t386) * t423 + t425 * t422;
t326 = t278 * t424 + (-t272 * Ifges(6,3) + t383) * t390 + Ifges(5,6) * t422 + (t272 * Ifges(5,2) + t387) * t423;
t338 = t424 - Ifges(5,6) / 0.2e1;
t339 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t89 = pkin(8) * t355 - t276 * t408 - t126;
t99 = pkin(8) * t353 + t113;
t42 = t277 * t89 - t389 * t99;
t43 = t277 * t99 + t389 * t89;
t93 = Ifges(7,4) * t180 - Ifges(7,2) * t179 - Ifges(7,6) * t278;
t94 = Ifges(7,2) * t181 - Ifges(7,6) * t276 - t385;
t95 = Ifges(7,1) * t180 - Ifges(7,4) * t179 - Ifges(7,5) * t278;
t1 = t320 * t317 + m(5) * (t126 * t128 + t127 * t129 - t236 * t237) + m(6) * (t113 * t116 + t114 * t118 + t154 * t155) + m(7) * (t130 * t131 + t42 * t44 + t43 * t45) + (-pkin(1) * mrSges(3,1) - t232 * mrSges(4,2) + ((Ifges(5,6) - Ifges(6,6)) * t272 + t425 * t271 + t421) * t276 + (Ifges(4,2) - Ifges(4,3) + Ifges(3,1) - Ifges(3,2) + Ifges(6,2) + Ifges(5,3) + Ifges(7,3) + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1) * t270 + ((-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t271 + (-Ifges(5,4) + Ifges(6,5)) * t272) * t271) * t278 + t301) * t276 + t95 * t401 + t96 * t404 + t94 * t406 + (Ifges(7,5) * t400 + Ifges(7,6) * t403 - t232 * mrSges(4,3) - pkin(1) * mrSges(3,2) + (-t236 * mrSges(5,1) + t326) * t272 + (t236 * mrSges(5,2) + t325) * t271 + (-t271 * t339 + t272 * t338 - t421) * t278) * t278 + t237 * t203 + t126 * t221 + t114 * t222 + t128 * t223 + t118 * t224 + t127 * t225 + t129 * t226 + t116 * t227 + t113 * t228 - t155 * t202 + t154 * t204 + t130 * t102 + t131 * t101 + t43 * t148 + t42 * t149 + t45 * t150 + t44 * t151 + t181 * t93 / 0.2e1;
t381 = t1 * qJD(1);
t373 = t298 * mrSges(7,3);
t369 = t277 * mrSges(7,2);
t311 = Ifges(7,5) * t181 + Ifges(7,6) * t182;
t312 = Ifges(7,1) * t181 + t385;
t52 = t182 * mrSges(7,1) - t181 * mrSges(7,2);
t4 = t131 * t52 + t311 * t390 + t312 * t400 + t94 * t401 + (t151 - t375) * t43 + (-t150 + t377) * t42 + t420 * t403;
t367 = t4 * qJD(1);
t333 = t389 * t151;
t351 = t277 * t150;
t16 = (m(6) * t155 - m(7) * t131 - t102 + t204) * t355 + (m(7) * (t277 * t43 - t389 * t42) - t333 + t351 + m(6) * t113 + t227) * t276;
t365 = qJD(1) * t16;
t293 = m(7) * (t215 * t277 + t298 * t389) * t390;
t295 = m(7) * (-t179 * t389 - t277 * t180);
t41 = t293 + m(6) * t356 - t295 / 0.2e1;
t364 = qJD(1) * t41;
t291 = (t181 * t399 + t182 * t396) * mrSges(7,3) + t150 * t398 + t151 * t397;
t300 = mrSges(7,1) * t405 + mrSges(7,2) * t407;
t11 = t291 - t300;
t363 = t11 * qJD(1);
t287 = -t418 * t272 + (t223 - t224) * t271 + m(6) * (-t113 * t272 - t114 * t271) + m(5) * (t126 * t271 - t127 * t272);
t13 = m(7) * (-t179 * t43 - t180 * t42) - t180 * t151 - t179 * t150 + (t287 - t320) * t276;
t362 = t13 * qJD(1);
t361 = t179 * t215;
t360 = t180 * t298;
t359 = t181 * t215;
t358 = t182 * t298;
t342 = -t389 / 0.2e1;
t294 = t333 / 0.2e1 + t342 * t375;
t315 = t276 * t342;
t20 = mrSges(7,1) * t315 + (-t370 / 0.2e1 + t377 / 0.2e1 - t150 / 0.2e1) * t277 + t294;
t357 = t20 * qJD(1);
t350 = t277 * t182;
t349 = t277 * t298;
t347 = Ifges(7,5) * t215 - Ifges(7,6) * t298;
t346 = t271 * mrSges(5,1) + t272 * mrSges(5,2);
t345 = t271 ^ 2 + t270;
t344 = m(6) / 0.4e1 + m(5) / 0.4e1;
t343 = (t358 - t359) * t410;
t341 = -mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1;
t340 = -mrSges(6,3) / 0.2e1 + mrSges(5,2) / 0.2e1;
t337 = t389 * mrSges(7,1);
t336 = -t377 / 0.2e1;
t335 = t373 / 0.2e1;
t334 = mrSges(7,3) * t396;
t332 = t389 * t181;
t331 = t389 * t215;
t329 = -t355 / 0.2e1;
t328 = t355 / 0.2e1;
t324 = -t222 / 0.2e1 + t221 / 0.2e1;
t323 = t228 / 0.2e1 + t225 / 0.2e1;
t322 = -t215 ^ 2 - t298 ^ 2;
t321 = -m(4) * qJ(3) - mrSges(4,3);
t316 = t415 + t413 + t411;
t65 = -mrSges(7,1) * t298 - t215 * mrSges(7,2);
t208 = Ifges(7,4) * t215;
t136 = -Ifges(7,2) * t298 + t208;
t137 = Ifges(7,2) * t215 + t384;
t138 = Ifges(7,1) * t215 - t384;
t139 = Ifges(7,1) * t298 + t208;
t18 = t201 * t65 + (t139 / 0.2e1 + t136 / 0.2e1) * t215 - (-t138 / 0.2e1 + t137 / 0.2e1) * t298;
t229 = t388 * t271;
t230 = t388 * t272;
t132 = -t229 * t389 - t277 * t230;
t133 = t277 * t229 - t230 * t389;
t282 = (t132 * t403 + t133 * t400) * mrSges(7,3) - t131 * t65 / 0.2e1 + t132 * t150 / 0.2e1 - t133 * t151 / 0.2e1 + t201 * t52 / 0.2e1 - t276 * t347 / 0.4e1 + (t136 + t139) * t181 / 0.4e1 + t420 * t215 / 0.4e1 - (t94 / 0.4e1 - t312 / 0.4e1) * t298 + (t137 / 0.4e1 - t138 / 0.4e1) * t182;
t289 = Ifges(7,3) * t422 + t44 * mrSges(7,1) / 0.2e1 - t45 * mrSges(7,2) / 0.2e1 - t301;
t2 = t282 - t289;
t309 = t2 * qJD(1) + t18 * qJD(2);
t28 = t322 * mrSges(7,3) + m(7) * (t132 * t298 - t133 * t215) + (-0.4e1 * t273 * t344 + mrSges(6,2) + mrSges(5,3)) * t345;
t284 = (-t272 * t126 - t271 * t127) * t415 + (-t271 * t113 + t272 * t114) * t413 + (-t132 * t182 - t133 * t181 - t215 * t43 + t298 * t42) * t411 + t151 * t398 + t150 * t396;
t290 = t236 * t415 + t154 * t412 + t130 * t411 - t380 / 0.2e1 - t379 / 0.2e1;
t6 = (t359 / 0.2e1 - t358 / 0.2e1) * mrSges(7,3) + (-t224 / 0.2e1 + t223 / 0.2e1 + t341 * t276) * t272 + (t227 / 0.2e1 + t226 / 0.2e1 + t340 * t276) * t271 + t284 + t290;
t308 = -qJD(1) * t6 + qJD(2) * t28;
t47 = -m(5) * qJ(3) + t321 - t346 + t416;
t285 = t237 * t415 + t155 * t413 + (-t132 * t180 - t133 * t179 - t131) * t411 - t378 / 0.2e1 - t376 / 0.2e1;
t302 = t272 * t128 + t271 * t129;
t303 = t271 * t116 - t272 * t118;
t286 = t148 * t396 + t149 * t398 + t302 * t414 + t303 * t412 + t410 * t417;
t8 = (t361 / 0.2e1 - t360 / 0.2e1) * mrSges(7,3) + (t278 * t341 + t324) * t272 + (t278 * t340 + t323) * t271 + t285 + t286;
t307 = -qJD(1) * t8 - qJD(2) * t47;
t51 = t416 * t272;
t283 = (-t272 * t155 + (t231 * t278 + t273 * t276) * t271) * t412 + (t201 * t355 + t272 * t131 + (-t132 * t389 + t133 * t277) * t276) * t410 + t102 * t391 + t204 * t392 + t135 * t329 + t233 * t328 + t334 * t352 - t315 * t373;
t288 = t118 * t412 + (t277 * t44 + t389 * t45) * t410 + t277 * t149 / 0.2e1 - t368 / 0.2e1 + t389 * t148 / 0.2e1;
t9 = -t252 + t283 - t288;
t306 = -qJD(1) * t9 - qJD(2) * t51;
t305 = qJD(1) * t52 + qJD(2) * t65;
t63 = m(6) * t355 + (t328 + t332 / 0.2e1 + t350 / 0.2e1) * m(7);
t88 = m(6) * t272 + (-t331 / 0.2e1 + t349 / 0.2e1 + t391) * m(7);
t304 = qJD(1) * t63 - qJD(2) * t88;
t14 = m(7) * (-t181 * t43 - t182 * t42) - t182 * t151 - t181 * t150 + t287 * t278;
t297 = t14 * qJD(1) + qJD(3) * t343;
t292 = t322 * t410 - 0.2e1 * t345 * t344;
t49 = t292 + t316;
t296 = qJD(1) * t343 + t49 * qJD(2);
t97 = m(7) * t392 + (-t331 + t349) * t410;
t76 = m(7) * t328 + (-t332 - t350) * t410;
t59 = qJD(4) * t343;
t48 = t292 - t316;
t46 = t295 / 0.2e1 + t293;
t21 = t277 * t336 + t351 / 0.2e1 + (-t369 / 0.2e1 - t337 / 0.2e1) * t276 - t294;
t12 = t291 + t300;
t10 = t283 + t288;
t7 = t215 * t336 + t182 * t335 + t224 * t391 + t223 * t392 + (t271 * t340 + t272 * t341) * t276 - t284 + t290 - t418 * t271 / 0.2e1;
t5 = -t285 + mrSges(5,2) * t329 + mrSges(6,3) * t328 + t323 * t271 + t324 * t272 - t179 * t334 + t180 * t335 + (m(4) * pkin(7) + mrSges(4,1)) * t278 + t286 + (mrSges(5,1) + mrSges(6,1)) * t353 / 0.2e1;
t3 = t282 + t289;
t15 = [qJD(2) * t1 + qJD(3) * t13 + qJD(4) * t14 + qJD(5) * t16 - qJD(6) * t4, t5 * qJD(3) + t7 * qJD(4) + t10 * qJD(5) + t3 * qJD(6) + t381 + (mrSges(7,3) * t417 + qJ(3) * t203 - t201 * t101 + t130 * t135 + t132 * t149 + t133 * t148 + t137 * t406 + t139 * t404 + t154 * t233 - t231 * t202 - t236 * t346 + t93 * t396 + t95 * t399 + (-t128 * mrSges(5,3) + t118 * mrSges(6,2) + (t221 - t222) * t273 - t325) * t272 + (-t129 * mrSges(5,3) - t116 * mrSges(6,2) + (t225 + t228) * t273 + t326) * t271 + 0.2e1 * (t154 * t231 + t273 * t303) * t412 + 0.2e1 * (-qJ(3) * t236 + t273 * t302) * t414 + 0.2e1 * (-t130 * t201 + t132 * t44 + t133 * t45) * t410 + (-pkin(2) * mrSges(4,1) + Ifges(3,5) - Ifges(4,4) + Ifges(7,5) * t398 + Ifges(7,6) * t397 + t339 * t272 + t338 * t271 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t278 + (-qJ(3) * mrSges(4,1) + (-Ifges(5,2) * t271 + t386) * t391 + (Ifges(6,3) * t271 + t382) * t392 - Ifges(3,6) + Ifges(4,5) + (mrSges(3,2) + t321) * pkin(7) + (t426 * t272 + t383 - t387) * t271 / 0.2e1) * t276) * qJD(2), t362 + t5 * qJD(2) + m(7) * (t360 - t361) * qJD(3) + t59 + t46 * qJD(5) + t12 * qJD(6), t7 * qJD(2) + t76 * qJD(5) + t297, qJD(2) * t10 + qJD(3) * t46 + qJD(4) * t76 + qJD(6) * t21 + t365, -t367 + t3 * qJD(2) + t12 * qJD(3) + t21 * qJD(5) + (-mrSges(7,1) * t43 - mrSges(7,2) * t42 + t311) * qJD(6); -qJD(3) * t8 - qJD(4) * t6 + qJD(5) * t9 + qJD(6) * t2 - t381, -qJD(3) * t47 + qJD(4) * t28 + qJD(5) * t51 + qJD(6) * t18, qJD(4) * t48 + t307, qJD(3) * t48 + qJD(5) * t97 + t308, qJD(4) * t97 - t306 (-mrSges(7,1) * t133 - mrSges(7,2) * t132 + t347) * qJD(6) + t309; qJD(2) * t8 + qJD(5) * t41 + qJD(6) * t11 - t362 + t59, qJD(4) * t49 - t307, 0, t296, t364, qJD(6) * t135 + t363; t6 * qJD(2) + t63 * qJD(5) + t52 * qJD(6) - t297, -qJD(3) * t49 - qJD(5) * t88 + qJD(6) * t65 - t308, -t296, 0, t304, t305; -qJD(2) * t9 - qJD(3) * t41 - qJD(4) * t63 - qJD(6) * t20 - t365, qJD(4) * t88 + t306, -t364, -t304, 0, -t357 + (-t337 - t369) * qJD(6); -qJD(2) * t2 - qJD(3) * t11 - qJD(4) * t52 + qJD(5) * t20 + t367, -qJD(4) * t65 - t309, -t363, -t305, t357, 0;];
Cq  = t15;
