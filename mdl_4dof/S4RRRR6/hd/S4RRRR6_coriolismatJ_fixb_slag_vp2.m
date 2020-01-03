% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:25
% EndTime: 2019-12-31 17:29:37
% DurationCPUTime: 5.16s
% Computational Cost: add. (8087->499), mult. (20615->747), div. (0->0), fcn. (20520->8), ass. (0->253)
t290 = sin(qJ(4));
t293 = cos(qJ(4));
t269 = mrSges(5,1) * t290 + mrSges(5,2) * t293;
t407 = Ifges(5,4) * t290;
t332 = Ifges(5,1) * t293 - t407;
t427 = -t293 / 0.2e1;
t431 = -t290 / 0.2e1;
t287 = Ifges(5,4) * t293;
t470 = -Ifges(5,2) * t290 + t287;
t478 = pkin(3) * t269 + t332 * t431 + t427 * t470;
t291 = sin(qJ(3));
t294 = cos(qJ(3));
t336 = mrSges(4,1) * t291 + mrSges(4,2) * t294;
t477 = -t336 / 0.2e1;
t289 = sin(pkin(4));
t292 = sin(qJ(2));
t375 = t289 * t292;
t476 = t375 / 0.2e1;
t295 = cos(qJ(2));
t379 = cos(pkin(4));
t348 = pkin(1) * t379;
t338 = t295 * t348;
t362 = pkin(6) * t375;
t245 = t338 - t362;
t474 = t245 * mrSges(3,2);
t286 = Ifges(5,5) * t293;
t399 = Ifges(5,6) * t290;
t471 = t286 - t399;
t227 = -Ifges(5,3) * t294 + t291 * t471;
t409 = Ifges(4,4) * t291;
t333 = Ifges(4,1) * t294 - t409;
t473 = t227 + t333;
t288 = Ifges(4,4) * t294;
t274 = Ifges(4,1) * t291 + t288;
t342 = -Ifges(4,2) * t291 + t288;
t472 = t274 + t342;
t273 = Ifges(5,1) * t290 + t287;
t277 = pkin(3) * t291 - pkin(8) * t294;
t371 = t290 * t291;
t220 = pkin(7) * t371 + t277 * t293;
t369 = t291 * t293;
t221 = -pkin(7) * t369 + t277 * t290;
t469 = -t220 * t290 + t221 * t293;
t374 = t289 * t295;
t244 = t379 * t291 + t294 * t375;
t389 = t244 * mrSges(4,3);
t325 = -mrSges(4,1) * t374 - t389;
t467 = t389 + t325;
t334 = mrSges(5,1) * t293 - mrSges(5,2) * t290;
t465 = -m(5) * pkin(3) - mrSges(4,1) - t334;
t464 = t289 ^ 2;
t463 = m(5) / 0.2e1;
t462 = -pkin(3) / 0.2e1;
t461 = -Ifges(5,3) / 0.2e1;
t246 = (pkin(2) * t292 - pkin(7) * t295) * t289;
t171 = t294 * t245 + t291 * t246;
t154 = pkin(8) * t375 + t171;
t247 = pkin(6) * t374 + t292 * t348;
t178 = t277 * t374 + t247;
t69 = -t154 * t290 + t178 * t293;
t460 = -t69 / 0.2e1;
t70 = t154 * t293 + t178 * t290;
t459 = t70 / 0.2e1;
t190 = -t244 * t290 - t293 * t374;
t243 = t291 * t375 - t379 * t294;
t191 = t244 * t293 - t290 * t374;
t408 = Ifges(5,4) * t191;
t74 = Ifges(5,2) * t190 + Ifges(5,6) * t243 + t408;
t458 = -t74 / 0.2e1;
t457 = -t74 / 0.4e1;
t456 = pkin(2) * mrSges(4,1);
t455 = pkin(2) * mrSges(4,2);
t454 = pkin(7) * mrSges(4,2);
t453 = pkin(8) * mrSges(5,3);
t106 = mrSges(5,1) * t191 + mrSges(5,2) * t190;
t452 = t106 / 0.2e1;
t189 = Ifges(5,4) * t190;
t108 = -Ifges(5,2) * t191 + t189;
t451 = -t108 / 0.4e1;
t365 = t294 * t295;
t215 = (-t290 * t365 + t292 * t293) * t289;
t216 = (t290 * t292 + t293 * t365) * t289;
t368 = t291 * t295;
t350 = t289 * t368;
t117 = Ifges(5,1) * t216 + Ifges(5,4) * t215 + Ifges(5,5) * t350;
t450 = t117 / 0.4e1;
t223 = t379 * pkin(7) + t247;
t224 = (-pkin(2) * t295 - pkin(7) * t292 - pkin(1)) * t289;
t146 = -t291 * t223 + t294 * t224;
t126 = pkin(3) * t374 - t146;
t449 = -t126 / 0.2e1;
t448 = -t190 / 0.4e1;
t447 = -t191 / 0.4e1;
t267 = -pkin(3) * t294 - pkin(8) * t291 - pkin(2);
t370 = t290 * t294;
t217 = -pkin(7) * t370 + t267 * t293;
t446 = t217 / 0.2e1;
t366 = t293 * t294;
t218 = pkin(7) * t366 + t267 * t290;
t445 = t218 / 0.2e1;
t444 = -t220 / 0.2e1;
t442 = t243 / 0.2e1;
t441 = t243 / 0.4e1;
t260 = mrSges(5,2) * t294 - mrSges(5,3) * t371;
t439 = -t260 / 0.2e1;
t262 = -mrSges(5,1) * t294 - mrSges(5,3) * t369;
t438 = -t262 / 0.2e1;
t398 = Ifges(5,6) * t293;
t404 = Ifges(5,5) * t290;
t270 = t398 + t404;
t437 = t270 / 0.4e1;
t271 = Ifges(5,2) * t293 + t407;
t436 = t271 / 0.4e1;
t272 = Ifges(4,2) * t294 + t409;
t435 = -t272 / 0.2e1;
t434 = -t273 / 0.4e1;
t433 = -t286 / 0.4e1;
t432 = -t289 / 0.2e1;
t430 = -t290 / 0.4e1;
t429 = t290 / 0.2e1;
t426 = t293 / 0.2e1;
t425 = t293 / 0.4e1;
t424 = -t294 / 0.2e1;
t423 = -t294 / 0.4e1;
t422 = t294 / 0.2e1;
t421 = t294 / 0.4e1;
t419 = pkin(7) * t291;
t418 = pkin(7) * t294;
t310 = -t379 * pkin(2) - t338;
t304 = t310 + t362;
t121 = t243 * pkin(3) - t244 * pkin(8) + t304;
t147 = t294 * t223 + t291 * t224;
t127 = -pkin(8) * t374 + t147;
t38 = t121 * t293 - t127 * t290;
t417 = t38 * mrSges(5,3);
t39 = t121 * t290 + t127 * t293;
t416 = t39 * mrSges(5,3);
t415 = mrSges(5,3) * t243;
t413 = Ifges(3,4) * t292;
t412 = Ifges(3,4) * t295;
t411 = Ifges(4,4) * t243;
t410 = Ifges(4,4) * t244;
t406 = Ifges(4,5) * t294;
t405 = Ifges(5,5) * t216;
t403 = Ifges(5,5) * t294;
t400 = Ifges(5,6) * t215;
t397 = Ifges(5,6) * t294;
t396 = Ifges(4,3) * t292;
t395 = Ifges(5,3) * t244;
t394 = Ifges(5,3) * t291;
t393 = t191 * mrSges(5,3);
t392 = t217 * mrSges(5,3);
t391 = t218 * mrSges(5,3);
t390 = t243 * mrSges(4,3);
t387 = t290 * t74;
t75 = Ifges(5,1) * t191 + Ifges(5,5) * t243 + t189;
t385 = t293 * t75;
t116 = Ifges(5,4) * t216 + Ifges(5,2) * t215 + Ifges(5,6) * t350;
t141 = -mrSges(5,2) * t243 + t190 * mrSges(5,3);
t142 = mrSges(5,1) * t243 - t393;
t148 = -mrSges(5,1) * t215 + mrSges(5,2) * t216;
t170 = -t245 * t291 + t246 * t294;
t153 = -pkin(3) * t375 - t170;
t179 = -mrSges(5,2) * t350 + t215 * mrSges(5,3);
t180 = mrSges(5,1) * t350 - t216 * mrSges(5,3);
t302 = t304 * t336;
t360 = Ifges(4,5) * t374;
t305 = t294 * (Ifges(4,1) * t244 - t360 - t411);
t321 = -Ifges(4,6) * t374 + t410;
t306 = t291 * (-Ifges(4,2) * t243 + t321);
t330 = Ifges(5,5) * t191 + Ifges(5,6) * t190;
t308 = t291 * (Ifges(5,3) * t243 + t330);
t314 = t295 * (-Ifges(4,6) * t291 + t406);
t323 = mrSges(4,2) * t374 - t390;
t329 = t400 + t405;
t335 = mrSges(5,1) * t190 - mrSges(5,2) * t191;
t339 = Ifges(5,3) * t350;
t347 = -t374 / 0.2e1;
t3 = -(t295 * (-Ifges(3,2) * t292 + t412) + t292 * (Ifges(3,1) * t295 - t413)) * t464 / 0.2e1 + (t308 + t305) * t347 + t215 * t458 + (pkin(1) * (mrSges(3,1) * t292 + mrSges(3,2) * t295) + t295 * (t314 + t396) / 0.2e1) * t464 - (Ifges(4,5) * t244 - Ifges(4,6) * t243 - Ifges(4,3) * t374) * t375 / 0.2e1 + t306 * t374 / 0.2e1 - t302 * t374 - t243 * (t329 + t339) / 0.2e1 + t153 * t335 - t171 * t323 - t170 * t325 + ((Ifges(3,1) * t292 + t412) * t347 + (Ifges(3,2) * t295 + t413) * t476 - t146 * (mrSges(4,1) * t292 - mrSges(4,3) * t365) - t147 * (-mrSges(4,2) * t292 - mrSges(4,3) * t368) + (Ifges(4,6) * t292 + t342 * t295) * t442) * t289 + (Ifges(3,5) * t347 + Ifges(3,6) * t476 + t474 + t247 * mrSges(3,1) + (Ifges(3,5) * t295 - Ifges(3,6) * t292) * t432) * t379 - m(4) * (t146 * t170 + t147 * t171 + t304 * t247) - m(5) * (t126 * t153 + t38 * t69 + t39 * t70) - t247 * (mrSges(4,1) * t243 + mrSges(4,2) * t244) + t244 * (Ifges(4,5) * t292 + t333 * t295) * t432 - t70 * t141 - t69 * t142 - t126 * t148 - t39 * t179 - t38 * t180 - t190 * t116 / 0.2e1 - t191 * t117 / 0.2e1 - t216 * t75 / 0.2e1;
t384 = t3 * qJD(1);
t313 = t269 * t243;
t349 = -t286 / 0.2e1;
t320 = t349 + t399 / 0.2e1;
t322 = -mrSges(5,2) * t244 + t290 * t415;
t324 = mrSges(5,1) * t244 + t293 * t415;
t337 = mrSges(4,1) * t244 - mrSges(4,2) * t243;
t174 = pkin(3) * t244 + pkin(8) * t243;
t61 = -t146 * t290 + t174 * t293;
t62 = t146 * t293 + t174 * t290;
t4 = t126 * t313 - m(5) * (t38 * t61 + t39 * t62) - t62 * t141 - t39 * t322 - t61 * t142 - t38 * t324 - t146 * t323 - t304 * t337 + (t321 - t330) * t244 + (-t387 / 0.2e1 - t146 * mrSges(4,3) - t360 + t385 / 0.2e1 + t191 * t332 / 0.2e1 + t190 * t470 / 0.2e1 + (-Ifges(4,4) - t320) * t243 + (-Ifges(5,3) + Ifges(4,1) - Ifges(4,2)) * t244) * t243 + (-m(5) * t126 + t335 + t467) * t147;
t383 = t4 * qJD(1);
t382 = t69 * t290;
t107 = Ifges(5,5) * t190 - Ifges(5,6) * t191;
t109 = Ifges(5,1) * t190 - t408;
t7 = t126 * t106 + t38 * t141 - t39 * t142 + t107 * t442 + (t109 / 0.2e1 + t458 - t416) * t191 + (t75 / 0.2e1 + t108 / 0.2e1 - t417) * t190;
t381 = t7 * qJD(1);
t380 = t70 * t293;
t378 = t171 * t294;
t373 = t290 * t190;
t228 = t291 * t470 - t397;
t372 = t290 * t228;
t230 = t332 * t291 - t403;
t367 = t293 * t230;
t363 = -Ifges(4,6) + t270 / 0.2e1;
t361 = -t453 / 0.2e1;
t357 = pkin(7) * t269 / 0.2e1;
t354 = pkin(8) * t179 / 0.2e1;
t353 = pkin(8) * t431;
t351 = -Ifges(4,2) / 0.2e1 + t461;
t346 = -t370 / 0.2e1;
t345 = t366 / 0.2e1;
t254 = t291 * t273;
t344 = -t228 / 0.4e1 - t254 / 0.4e1;
t253 = t271 * t291;
t343 = t230 / 0.4e1 - t253 / 0.4e1;
t229 = t291 * Ifges(5,6) + t294 * t470;
t231 = t291 * Ifges(5,5) + t332 * t294;
t250 = t269 * t291;
t251 = t269 * t294;
t261 = -mrSges(5,2) * t291 - mrSges(5,3) * t370;
t263 = mrSges(5,1) * t291 - mrSges(5,3) * t366;
t307 = t294 * t471 + t394;
t14 = t250 * t418 + t251 * t419 + t230 * t345 + t231 * t369 / 0.2e1 + t228 * t346 - t229 * t371 / 0.2e1 + t307 * t424 + t221 * t260 + t218 * t261 + t220 * t262 + t217 * t263 + m(5) * (t217 * t220 + t218 * t221) - pkin(2) * t336 + t472 * t422 + (m(5) * pkin(7) ^ 2 * t294 + t435 + t473 / 0.2e1) * t291;
t296 = -m(5) * (t217 * t61 + t218 * t62 + t220 * t38 + t221 * t39 + (t126 * t294 + t147 * t291) * pkin(7)) / 0.2e1 + t251 * t449 - t147 * t250 / 0.2e1 + t229 * t448 + t231 * t447 + t142 * t444 - t221 * t141 / 0.2e1 - t38 * t263 / 0.2e1 - t39 * t261 / 0.2e1 + t61 * t438 + t62 * t439 + t335 * t418 / 0.2e1;
t299 = (-pkin(3) * t153 + (t380 - t382) * pkin(8)) * t463 + t148 * t462 - t153 * t334 / 0.2e1 + t170 * mrSges(4,1) / 0.2e1 - t171 * mrSges(4,2) / 0.2e1 + t215 * t436 + t216 * t273 / 0.4e1;
t2 = (t367 / 0.4e1 - t372 / 0.4e1 + t288 / 0.4e1 + t274 / 0.4e1 - t455 / 0.2e1 + (t217 * t427 + t218 * t431) * mrSges(5,3) + (t357 + t332 * t425 + t470 * t430 + Ifges(4,1) / 0.4e1 + t351) * t291 + t320 * t294) * t243 + (-t217 * mrSges(5,1) / 0.2e1 + mrSges(5,2) * t445 + t272 / 0.4e1 + t456 / 0.2e1 - t227 / 0.4e1 + (t433 + t399 / 0.4e1) * t291 + (Ifges(4,2) / 0.4e1 - Ifges(4,1) / 0.2e1 + Ifges(5,3) / 0.4e1) * t294) * t244 + (t75 * t423 + t354 + mrSges(5,3) * t459 + t116 / 0.4e1) * t293 + ((Ifges(4,3) / 0.2e1 + pkin(6) * t477) * t292 + ((Ifges(4,5) - pkin(7) * mrSges(4,1) / 0.2e1) * t294 + (-Ifges(4,6) + t454 / 0.2e1 + t437) * t291) * t295) * t289 + t296 + (t74 * t421 - pkin(8) * t180 / 0.2e1 + mrSges(5,3) * t460 + t450) * t290 + (0.3e1 / 0.4e1 * t291 * t244 + t243 * t422) * Ifges(4,4) + t310 * t477 - t291 * t330 / 0.4e1 + t299;
t327 = -t2 * qJD(1) + t14 * qJD(2);
t249 = t334 * t291;
t252 = t291 * t270;
t17 = -t252 * t424 + t217 * t260 - t218 * t262 + (pkin(7) * t249 + (-t254 / 0.2e1 - t228 / 0.2e1 - t391) * t293 + (-t230 / 0.2e1 + t253 / 0.2e1 + t392) * t290) * t291;
t298 = (-t392 / 0.2e1 + t343) * t190 + (-t391 / 0.2e1 + t344) * t191 + t126 * t249 / 0.2e1 + t141 * t446 - t218 * t142 / 0.2e1 - t252 * t441 + t107 * t423 + t38 * t260 / 0.2e1 + t39 * t438;
t301 = pkin(7) * t452 + (t109 / 0.4e1 + t457 - t416 / 0.2e1) * t293 + (t451 - t75 / 0.4e1 + t417 / 0.2e1) * t290;
t303 = -t405 / 0.2e1 - t400 / 0.2e1 + mrSges(5,1) * t460 + mrSges(5,2) * t459;
t5 = (Ifges(5,3) * t347 + t301) * t291 + t298 + t303;
t326 = t5 * qJD(1) + t17 * qJD(2);
t319 = t249 * t462 + t286 * t423;
t318 = mrSges(5,1) * t444 + t221 * mrSges(5,2) / 0.2e1;
t317 = pkin(8) * t439 + t344;
t316 = pkin(8) * t438 + t343;
t315 = t271 * t429 + t273 * t427;
t312 = Ifges(4,5) - t315;
t300 = t357 + (t434 - t287 / 0.4e1 + (Ifges(5,2) / 0.4e1 + t361) * t290) * t290 + (-t407 / 0.4e1 - t271 / 0.4e1 + (Ifges(5,1) / 0.4e1 + t361) * t293) * t293;
t16 = (-t403 / 0.2e1 + t316) * t293 + (0.3e1 / 0.4e1 * t397 + t317) * t290 + (t461 + t300) * t291 + t318 + t319;
t297 = pkin(3) * t452 + t269 * t449 + t190 * t434 + t470 * t448 + t191 * t436 + t332 * t447 + t109 * t430 + t387 / 0.4e1 + t293 * t451 - t385 / 0.4e1;
t309 = t395 / 0.2e1 + t61 * mrSges(5,1) / 0.2e1 - t62 * mrSges(5,2) / 0.2e1;
t8 = (t141 * t429 + t142 * t426 + (t191 * t426 - t373 / 0.2e1) * mrSges(5,3)) * pkin(8) + t297 + (t349 + 0.3e1 / 0.4e1 * t399 + t433) * t243 + t309;
t96 = t315 + t478;
t311 = t8 * qJD(1) - t16 * qJD(2) + t96 * qJD(3);
t15 = Ifges(5,5) * t345 + Ifges(5,6) * t346 + t394 / 0.2e1 + (t397 / 0.4e1 + t317) * t290 + t316 * t293 + t300 * t291 - t318 + t319;
t9 = t320 * t243 - t297 + t141 * t353 + t373 * t453 / 0.2e1 + t471 * t441 + t309 + (t142 + t393) * pkin(8) * t427;
t6 = t298 + t339 / 0.2e1 + t301 * t291 - t303;
t1 = t308 / 0.4e1 - t306 / 0.4e1 + t305 / 0.4e1 + t299 - t296 + t370 * t457 + t290 * t450 + t322 * t445 + t324 * t446 - (t390 + t323 + t313) * t419 / 0.2e1 + (-t243 * t471 + t395) * t423 - (Ifges(5,6) * t244 - t243 * t470) * t371 / 0.4e1 + t291 * (-Ifges(4,1) * t243 - t410) / 0.4e1 + (t380 / 0.2e1 - t382 / 0.2e1) * mrSges(5,3) + (Ifges(5,5) * t244 - t332 * t243) * t369 / 0.4e1 + t75 * t366 / 0.4e1 + (-Ifges(4,2) * t244 - t411) * t421 + t116 * t425 - pkin(2) * t337 / 0.2e1 - (t367 + t472) * t243 / 0.4e1 + t473 * t244 / 0.4e1 - t467 * t418 / 0.2e1 - t244 * t272 / 0.4e1 + t302 / 0.2e1 + t180 * t353 + t293 * t354 + (t396 / 0.2e1 + (t406 / 0.2e1 + (-Ifges(4,6) / 0.2e1 + t437) * t291) * t295 - t314 / 0.4e1) * t289 + (t372 + t307) * t441;
t10 = [-qJD(2) * t3 - qJD(3) * t4 + qJD(4) * t7, t1 * qJD(3) + t6 * qJD(4) - t384 + (mrSges(4,3) * t378 + t329 * t424 + t70 * t260 + t69 * t262 + t153 * t250 - t474 + t217 * t180 + t218 * t179 + t215 * t228 / 0.2e1 + t216 * t230 / 0.2e1 + (-t170 * mrSges(4,3) + pkin(7) * t148 + t116 * t431 + t117 * t426) * t291 + m(4) * (-t170 * t291 + t378) * pkin(7) + 0.2e1 * (t153 * t419 + t217 * t69 + t218 * t70) * t463 + ((Ifges(3,5) + (t274 / 0.2e1 - t455 + t288 / 0.2e1) * t294 + (t435 + t227 / 0.2e1 - t456 - t409 / 0.2e1 + (Ifges(4,1) / 0.2e1 + t351) * t294) * t291) * t295 + (Ifges(4,5) * t291 + Ifges(4,6) * t294 - t336 * pkin(7) - Ifges(3,6)) * t292) * t289 + (-m(4) * pkin(2) - t294 * mrSges(4,1) + t291 * mrSges(4,2) - mrSges(3,1)) * t247) * qJD(2), t1 * qJD(2) + t9 * qJD(4) - t383 + (-t146 * mrSges(4,2) + (t404 / 0.2e1 + t398 / 0.2e1 - t269 * pkin(8) + t363) * t244 + (-t312 + t478) * t243 + t465 * t147 + (m(5) * pkin(8) + mrSges(5,3)) * (-t290 * t61 + t293 * t62)) * qJD(3), t381 + t6 * qJD(2) + t9 * qJD(3) + (-mrSges(5,1) * t39 - mrSges(5,2) * t38 + t107) * qJD(4); -qJD(3) * t2 + qJD(4) * t5 + t384, qJD(3) * t14 + qJD(4) * t17, t15 * qJD(4) + t327 + (-pkin(3) * t251 + t229 * t426 + t231 * t429 + (m(5) * t469 + t293 * t261 - t290 * t263) * pkin(8) + (t465 * pkin(7) + t312) * t294 + (t363 + t454) * t291 + t469 * mrSges(5,3)) * qJD(3), t15 * qJD(3) + (-mrSges(5,1) * t218 - mrSges(5,2) * t217 - t252) * qJD(4) + t326; qJD(2) * t2 - qJD(4) * t8 + t383, qJD(4) * t16 - t327, -t96 * qJD(4), (-t334 * pkin(8) + t471) * qJD(4) - t311; -qJD(2) * t5 + qJD(3) * t8 - t381, -qJD(3) * t16 - t326, t311, 0;];
Cq = t10;
