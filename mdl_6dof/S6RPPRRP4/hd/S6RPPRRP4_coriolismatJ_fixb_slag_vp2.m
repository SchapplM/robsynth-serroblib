% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:13
% EndTime: 2019-03-09 02:05:23
% DurationCPUTime: 4.76s
% Computational Cost: add. (7210->487), mult. (13666->680), div. (0->0), fcn. (11450->6), ass. (0->263)
t460 = m(7) / 0.2e1;
t482 = Ifges(7,4) + Ifges(6,5);
t470 = m(7) / 0.4e1 + m(6) / 0.4e1;
t287 = sin(qJ(5));
t282 = t287 ^ 2;
t289 = cos(qJ(5));
t283 = t289 ^ 2;
t488 = t282 + t283;
t490 = (-0.1e1 + t488) * t470;
t326 = -pkin(5) * t289 - qJ(6) * t287;
t220 = -pkin(4) + t326;
t286 = cos(pkin(9));
t290 = cos(qJ(4));
t377 = t289 * t290;
t285 = sin(pkin(9));
t388 = t285 * t287;
t185 = t286 * t377 + t388;
t394 = t185 * t289;
t378 = t289 * t285;
t381 = t287 * t290;
t183 = t286 * t381 - t378;
t397 = t183 * t287;
t318 = t394 + t397;
t310 = t318 * pkin(8);
t288 = sin(qJ(4));
t386 = t286 * t288;
t461 = -m(7) / 0.2e1;
t464 = -m(6) / 0.2e1;
t483 = mrSges(7,2) + mrSges(6,3);
t489 = -t483 * (t394 / 0.2e1 + t397 / 0.2e1) + (-pkin(4) * t386 + t310) * t464 + (t220 * t386 + t310) * t461;
t291 = -pkin(1) - pkin(2);
t373 = qJ(2) * t286 + t285 * t291;
t206 = -pkin(7) + t373;
t337 = t206 * t287 - pkin(5);
t333 = -qJ(2) * t285 + t286 * t291;
t205 = pkin(3) - t333;
t125 = t290 * pkin(4) + pkin(8) * t288 + t205;
t399 = t125 * t289;
t64 = t290 * t337 - t399;
t66 = -t206 * t381 + t399;
t480 = t64 + t66;
t463 = m(6) / 0.2e1;
t469 = t460 + t463;
t362 = mrSges(7,3) / 0.2e1 - mrSges(6,2) / 0.2e1;
t487 = t362 * t185 + (-pkin(5) * t183 + qJ(6) * t185) * t460;
t182 = t285 * t381 + t286 * t289;
t184 = t285 * t377 - t286 * t287;
t391 = t206 * t289;
t334 = qJ(6) + t391;
t383 = t287 * t125;
t63 = t290 * t334 + t383;
t67 = t206 * t377 + t383;
t429 = -t63 + t67;
t485 = t182 * t429 + t184 * t480;
t484 = mrSges(6,1) + mrSges(7,1);
t481 = Ifges(7,2) + Ifges(6,3);
t361 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t479 = t290 * t361;
t382 = t287 * t288;
t208 = -mrSges(6,2) * t290 + mrSges(6,3) * t382;
t276 = t290 * mrSges(7,3);
t215 = mrSges(7,2) * t382 + t276;
t477 = t215 + t208;
t414 = t287 * mrSges(7,3);
t222 = -mrSges(7,1) * t289 - t414;
t338 = m(7) * t220 + t222;
t223 = -mrSges(6,1) * t289 + mrSges(6,2) * t287;
t476 = t223 + t222;
t225 = t287 * mrSges(7,1) - t289 * mrSges(7,3);
t403 = qJ(6) * t289;
t433 = pkin(5) * t287;
t224 = -t403 + t433;
t435 = m(7) * t224;
t475 = t225 + t435;
t380 = t288 * t289;
t474 = Ifges(6,6) * t380 + t382 * t482;
t278 = Ifges(7,5) * t287;
t228 = -Ifges(7,3) * t289 + t278;
t233 = Ifges(7,1) * t289 + t278;
t473 = Ifges(7,6) * t287 + t289 * t482;
t281 = Ifges(6,4) * t289;
t234 = Ifges(6,1) * t287 + t281;
t431 = pkin(8) * t290;
t434 = pkin(4) * t288;
t241 = t431 - t434;
t389 = t241 * t289;
t80 = t206 * t382 + t389;
t207 = t287 * t241;
t81 = -t206 * t380 + t207;
t471 = -t287 * t80 + t289 * t81;
t73 = -t288 * t334 + t207;
t74 = -t288 * t337 - t389;
t324 = t287 * t74 + t289 * t73;
t327 = Ifges(6,2) * t287 - t281;
t275 = m(7) * qJ(6) + mrSges(7,3);
t415 = t283 * mrSges(6,3);
t416 = t283 * mrSges(7,2);
t417 = t282 * mrSges(6,3);
t418 = t282 * mrSges(7,2);
t428 = m(7) * qJD(6);
t468 = (-mrSges(5,2) + t415 + t416 + t417 + t418) * qJD(4) + t287 * t428;
t194 = t290 * t225;
t226 = t287 * mrSges(6,1) + t289 * mrSges(6,2);
t195 = t290 * t226;
t343 = t215 / 0.2e1 + t208 / 0.2e1;
t210 = mrSges(6,1) * t290 + mrSges(6,3) * t380;
t211 = -mrSges(7,1) * t290 - mrSges(7,2) * t380;
t346 = t211 / 0.2e1 - t210 / 0.2e1;
t467 = -t346 * t287 - t343 * t289 - t194 / 0.2e1 - t195 / 0.2e1;
t466 = 0.2e1 * m(7);
t465 = t288 ^ 2;
t458 = mrSges(7,1) / 0.2e1;
t455 = Ifges(7,6) / 0.2e1;
t454 = t67 / 0.2e1;
t453 = t74 / 0.2e1;
t452 = -qJ(6) / 0.2e1;
t190 = t288 * t222;
t450 = t190 / 0.2e1;
t192 = t288 * t225;
t449 = -t192 / 0.2e1;
t260 = mrSges(7,2) * t377;
t411 = t288 * mrSges(7,1);
t213 = -t260 + t411;
t448 = t213 / 0.2e1;
t447 = t225 / 0.2e1;
t446 = t226 / 0.2e1;
t408 = t290 * mrSges(5,2);
t227 = -t288 * mrSges(5,1) - t408;
t445 = -t227 / 0.2e1;
t444 = t287 / 0.2e1;
t443 = -t288 / 0.2e1;
t442 = -t289 / 0.2e1;
t441 = t289 / 0.2e1;
t189 = t326 * t288;
t437 = m(7) * t189;
t284 = t290 ^ 2;
t25 = -0.2e1 * t285 * t465 * t490 + t469 * (t182 * t381 + t184 * t377 - t284 * t285);
t430 = t25 * qJD(4);
t427 = mrSges(7,3) * t288;
t424 = Ifges(6,4) * t287;
t423 = Ifges(7,5) * t289;
t421 = Ifges(7,6) * t289;
t419 = pkin(8) * qJD(4);
t407 = t290 * Ifges(6,6);
t406 = t290 * Ifges(7,6);
t405 = t223 - mrSges(5,1);
t82 = (t206 - t224) * t288;
t29 = -t192 * t380 + t290 * t215 + m(7) * (t290 * t63 + t380 * t82);
t402 = qJD(1) * t29;
t384 = t286 * t290;
t31 = 0.2e1 * t470 * (t318 - t384) * t288;
t401 = qJD(1) * t31;
t259 = mrSges(6,2) * t382;
t191 = -mrSges(6,1) * t380 + t259;
t349 = t450 + t191 / 0.2e1;
t315 = t437 / 0.2e1 + t349;
t12 = t315 * t290 + ((t287 * t429 + t289 * t480) * t461 + t211 * t442 + t210 * t441 + (-t416 / 0.2e1 - t418 / 0.2e1 - t415 / 0.2e1 - t417 / 0.2e1) * t288 + t477 * t444) * t288;
t400 = t12 * qJD(1);
t398 = t182 * t287;
t396 = t184 * t289;
t395 = t184 * t290;
t393 = t206 * t284;
t392 = t206 * t288;
t390 = t220 * t288;
t387 = t285 * t288;
t385 = t286 * t465;
t379 = t288 * t290;
t359 = t465 * t378;
t50 = (t183 / 0.4e1 - t359 / 0.4e1 - t395 / 0.4e1) * t466;
t376 = t50 * qJD(1);
t374 = t488 * t431;
t370 = qJD(4) * t288;
t369 = qJD(4) * t290;
t368 = qJD(5) * t288;
t366 = m(7) * t453;
t364 = t435 / 0.2e1;
t363 = mrSges(6,1) / 0.2e1 + t458;
t360 = t455 - Ifges(6,6) / 0.2e1;
t354 = -t381 / 0.2e1;
t353 = t380 / 0.2e1;
t229 = Ifges(7,3) * t287 + t423;
t151 = -t288 * t229 + t406;
t153 = t288 * t327 + t407;
t351 = -t151 / 0.2e1 + t153 / 0.2e1;
t155 = Ifges(7,4) * t290 - t288 * t233;
t264 = Ifges(6,4) * t382;
t157 = -Ifges(6,1) * t380 + Ifges(6,5) * t290 + t264;
t350 = -t155 / 0.2e1 - t157 / 0.2e1;
t193 = t288 * t226;
t348 = t449 - t193 / 0.2e1;
t212 = -mrSges(6,1) * t288 + mrSges(6,3) * t377;
t345 = t448 - t212 / 0.2e1;
t209 = mrSges(6,2) * t288 + mrSges(6,3) * t381;
t214 = mrSges(7,2) * t381 - t427;
t344 = t214 / 0.2e1 + t209 / 0.2e1;
t342 = t447 + t446;
t230 = Ifges(6,2) * t289 + t424;
t341 = -t228 / 0.2e1 + t230 / 0.2e1;
t232 = Ifges(7,1) * t287 - t423;
t340 = -t232 / 0.2e1 - t234 / 0.2e1;
t332 = -m(6) * pkin(4) + t405;
t330 = t362 * t289;
t235 = Ifges(6,1) * t289 - t424;
t325 = -t287 * t66 + t289 * t67;
t299 = t483 * (t396 / 0.2e1 + t398 / 0.2e1);
t300 = -t182 * t343 + t184 * t346;
t9 = t363 * t183 + (t285 * t349 + t299) * t288 + (t189 * t387 + t485) * t460 + t300 - t487;
t323 = t9 * qJD(1);
t314 = m(6) * t471;
t257 = qJ(6) * t377;
t83 = t257 + (t206 - t433) * t290;
t316 = t287 * t64 + t289 * t63 - t83;
t10 = t393 * t464 + (t344 * t289 + t345 * t287 + (t324 + t82) * t460 + t314 / 0.2e1 + t392 * t463 + t348) * t288 + (t316 * t460 + t325 * t463 - t467) * t290;
t152 = -Ifges(7,6) * t288 - t229 * t290;
t154 = -Ifges(6,6) * t288 + t290 * t327;
t156 = -Ifges(7,4) * t288 - t233 * t290;
t158 = -Ifges(6,5) * t288 - t235 * t290;
t3 = -t83 * t192 - t82 * t194 + t205 * t227 + t81 * t208 + t67 * t209 + t80 * t210 + t74 * t211 + t66 * t212 + t64 * t213 + t63 * t214 + t73 * t215 + m(6) * (t66 * t80 + t67 * t81) + m(7) * (t63 * t73 + t64 * t74 + t82 * t83) + (Ifges(5,4) * t290 - t206 * t193 + (t350 - t479) * t289 + (-t290 * t360 + t351) * t287) * t290 + (-Ifges(5,4) * t288 - t206 * t195 + (-t156 / 0.2e1 - t158 / 0.2e1 + t361 * t288) * t289 + (-t152 / 0.2e1 + t154 / 0.2e1 + t360 * t288) * t287 + (m(6) * t206 ^ 2 + Ifges(5,1) - Ifges(5,2) - t481) * t290) * t288;
t322 = qJD(1) * t3 + qJD(3) * t10;
t196 = t288 * t228;
t197 = Ifges(6,2) * t380 + t264;
t198 = t288 * t232;
t199 = t288 * t234;
t5 = -t189 * t192 + t82 * t190 + t66 * t215 + m(7) * (t189 * t82 + t63 * t66 + t64 * t67) + t66 * t208 - t67 * t210 + t67 * t211 + (t206 * t191 + (-t66 * mrSges(6,3) + t64 * mrSges(7,2) + t197 / 0.2e1 - t196 / 0.2e1 - t350) * t287 + (t63 * mrSges(7,2) - t406 / 0.2e1 + t67 * mrSges(6,3) - t198 / 0.2e1 - t199 / 0.2e1 + t351) * t289) * t288 + t474 * t290 / 0.2e1;
t321 = qJD(1) * t5 - qJD(3) * t12;
t163 = t206 * t385;
t11 = m(3) * qJ(2) + mrSges(3,3) + (t290 * mrSges(5,1) - t288 * mrSges(5,2) + mrSges(4,1)) * t285 + t477 * t185 + (-t210 + t211) * t183 + (-t284 * mrSges(5,3) + mrSges(4,2) + (-mrSges(5,3) * t288 - t192 - t193) * t288) * t286 + m(7) * (t183 * t64 + t185 * t63 + t386 * t82) + m(6) * (-t183 * t66 + t185 * t67 + t163) + m(5) * (t205 * t285 + t286 * t393 + t163) + m(4) * (-t285 * t333 + t286 * t373);
t320 = t11 * qJD(1) + t31 * qJD(3);
t301 = (-t287 * t82 + (t390 + t431) * t289) * t461 - t192 * t444;
t26 = -t260 + (t222 * t442 + t458) * t288 + t366 + t301;
t92 = t338 * t287;
t319 = qJD(1) * t26 + qJD(4) * t92;
t33 = -t276 + (t67 / 0.4e1 - t383 / 0.4e1 + (-t391 / 0.4e1 + t452) * t290) * t466;
t317 = qJD(1) * t33 - qJD(5) * t275;
t30 = (m(6) + m(7)) * (t285 ^ 2 * t379 + (-t396 - t398) * t387);
t292 = t345 * t182 + t344 * t184 + ((t460 * t82 + t348) * t290 + ((0.2e1 * t206 * t290 - t325) * t463 + t316 * t461 + t467) * t288) * t285 + (-t182 * t80 + t184 * t81) * t463 + (t182 * t74 + t184 * t73) * t460;
t7 = (t408 / 0.2e1 + t445 + (mrSges(5,1) / 0.2e1 - t223 / 0.2e1 - t222 / 0.2e1) * t288) * t286 + t292 + t489;
t309 = t7 * qJD(1) + t30 * qJD(2) + t25 * qJD(3);
t49 = 0.4e1 * t379 * t490;
t308 = t10 * qJD(1) + t25 * qJD(2) + t49 * qJD(3);
t306 = t330 + t342;
t15 = (-t363 * t287 + (t224 / 0.4e1 - t433 / 0.4e1 + t403 / 0.4e1) * t466 + t306) * t387;
t293 = t206 * t446 + (-t235 / 0.4e1 - t233 / 0.4e1 + t230 / 0.4e1 - t228 / 0.4e1) * t289 + (t234 / 0.4e1 + t232 / 0.4e1 - t327 / 0.4e1 - t229 / 0.4e1) * t287 + t483 * pkin(8) * (t283 / 0.2e1 + t282 / 0.2e1);
t294 = (t189 * t220 + t224 * t82) * t460 - pkin(4) * t191 / 0.2e1 + t189 * t222 / 0.2e1 + t220 * t450 + t224 * t449 + t82 * t447 + t473 * t290 / 0.4e1;
t296 = (-pkin(5) * t74 + qJ(6) * t73) * t461 + pkin(5) * t448 + t214 * t452 - t73 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t453 - t80 * mrSges(6,1) / 0.2e1 + t81 * mrSges(6,2) / 0.2e1;
t297 = (t454 - t63 / 0.2e1) * mrSges(7,2) + (t429 * t460 - t343) * pkin(8) + t151 / 0.4e1 - t153 / 0.4e1 + t198 / 0.4e1 + t199 / 0.4e1;
t298 = t197 / 0.4e1 - t196 / 0.4e1 + t157 / 0.4e1 + t155 / 0.4e1 + (t66 / 0.2e1 + t64 / 0.2e1) * mrSges(7,2) + (t460 * t480 + t346) * pkin(8);
t2 = t294 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + t293) * t288 + (t298 + t479) * t289 + ((-0.3e1 / 0.4e1 * Ifges(6,6) + t455) * t290 + t297) * t287 + t296;
t23 = -pkin(4) * t226 + t224 * t222 + t475 * t220 + (-t229 / 0.2e1 - t327 / 0.2e1 - t340) * t289 + (t235 / 0.2e1 + t233 / 0.2e1 - t341) * t287;
t38 = t257 * t460 + (t364 + (pkin(5) * t461 - t363) * t287 + t306) * t290;
t303 = qJD(1) * t2 + qJD(2) * t15 - qJD(3) * t38 + qJD(4) * t23;
t240 = t285 * t385;
t221 = (m(7) * pkin(8) + mrSges(7,2)) * t289;
t51 = (t359 + t395 + t183) * t460;
t39 = (-pkin(5) * t381 + t257) * t460 + t290 * t330 + t484 * t354 - (t226 + t475) * t290 / 0.2e1;
t32 = (t383 + (0.2e1 * qJ(6) + t391) * t290) * t460 + m(7) * t454 + t215;
t27 = t222 * t353 + t366 + t411 / 0.2e1 - t301;
t14 = -t378 * t427 / 0.2e1 + t285 * mrSges(6,2) * t353 + (t364 + t342) * t387 + (t287 * t484 + t435) * t387 / 0.2e1;
t8 = t485 * t460 + (t285 * t315 + t299) * t288 + t300 - t484 * t183 / 0.2e1 + t487;
t6 = -mrSges(5,2) * t384 / 0.2e1 + t286 * t445 + t292 + (-mrSges(5,1) / 0.2e1 + t476 / 0.2e1) * t386 - t489;
t4 = qJD(2) * t31 + qJD(4) * t10 - qJD(5) * t12;
t1 = -t296 + t294 + t293 * t288 + (-t407 / 0.4e1 + t297) * t287 + t298 * t289 + Ifges(6,6) * t381 / 0.2e1 + Ifges(7,6) * t354 + t481 * t443 - t482 * t377 / 0.2e1;
t13 = [qJD(2) * t11 + qJD(4) * t3 + qJD(5) * t5 + qJD(6) * t29, 0.2e1 * (m(5) * (t240 + (t284 - 0.1e1) * t286 * t285) / 0.2e1 + t469 * (t182 * t183 + t184 * t185 + t240)) * qJD(2) + t6 * qJD(4) + t8 * qJD(5) + t51 * qJD(6) + t320, t4, t6 * qJD(2) + (mrSges(5,2) * t392 + Ifges(5,6) * t288 + pkin(4) * t195 + t152 * t442 + t154 * t441 - t220 * t194 + t338 * t83 + (t156 + t158) * t444 + (Ifges(6,6) * t289 + t287 * t482 - t421) * t443 + t471 * mrSges(6,3) + t324 * mrSges(7,2)) * qJD(4) + t1 * qJD(5) + t27 * qJD(6) + ((t209 + t214) * t289 + (-t212 + t213) * t287 + t314 + m(7) * t324) * t419 + (t206 * t332 + t287 * t341 + t289 * t340 - Ifges(5,5)) * t369 + t322, t8 * qJD(2) + t1 * qJD(4) + ((-m(7) * pkin(5) - t484) * t67 + (-mrSges(6,2) + t275) * t66 + t474) * qJD(5) + t32 * qJD(6) + (-mrSges(7,2) * t224 - t421) * t368 + t321, qJD(2) * t51 + qJD(4) * t27 + qJD(5) * t32 + t402; qJD(4) * t7 + qJD(5) * t9 - qJD(6) * t50 - t320, qJD(4) * t30, -t401 + t430, t14 * qJD(5) + ((t332 + t338) * t369 + (-0.2e1 * t419 * t469 * t488 - t468) * t288) * t285 + t309, t14 * qJD(4) + (-t484 * t184 + (mrSges(6,2) - mrSges(7,3)) * t182) * qJD(5) + ((-pkin(5) * t184 - qJ(6) * t182) * qJD(5) / 0.2e1 + t184 * qJD(6) / 0.2e1) * t466 + t323, -t376 + (qJD(5) * t184 - t370 * t388) * m(7); t4, t401 + t430, t49 * qJD(4), t39 * qJD(5) + (t222 + t405) * t370 + 0.2e1 * ((t374 - t434) * t463 + (t374 + t390) * t460) * qJD(4) + t468 * t290 + t308, -t400 + t39 * qJD(4) + (t259 + t437) * qJD(5) + (-qJD(5) * t414 + (-t484 * qJD(5) + t428) * t289) * t288 (t287 * t369 + t289 * t368) * m(7); -qJD(2) * t7 + qJD(5) * t2 - qJD(6) * t26 - t322, qJD(5) * t15 - t309, -qJD(5) * t38 - t308, qJD(5) * t23 - qJD(6) * t92, t221 * qJD(6) + t303 + (-Ifges(6,6) * t287 + (m(7) * t326 + t476) * pkin(8) + t326 * mrSges(7,2) + t473) * qJD(5), qJD(5) * t221 - t319; -qJD(2) * t9 - qJD(4) * t2 - qJD(6) * t33 - t321, -qJD(4) * t15 - t323, qJD(4) * t38 + t400, -t303, t275 * qJD(6), -t317; qJD(2) * t50 + qJD(4) * t26 + qJD(5) * t33 - t402, t376, 0, t319, t317, 0;];
Cq  = t13;
