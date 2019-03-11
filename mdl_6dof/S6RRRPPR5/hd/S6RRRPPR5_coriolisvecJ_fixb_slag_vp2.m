% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:22
% EndTime: 2019-03-09 15:39:35
% DurationCPUTime: 39.89s
% Computational Cost: add. (22863->913), mult. (60834->1294), div. (0->0), fcn. (48790->12), ass. (0->377)
t358 = cos(qJ(2));
t355 = sin(qJ(2));
t349 = sin(pkin(6));
t403 = qJD(1) * t349;
t389 = t355 * t403;
t352 = cos(pkin(6));
t402 = qJD(1) * t352;
t395 = pkin(1) * t402;
t302 = -pkin(8) * t389 + t358 * t395;
t360 = (pkin(2) * t355 - pkin(9) * t358) * t349;
t303 = qJD(1) * t360;
t354 = sin(qJ(3));
t357 = cos(qJ(3));
t235 = -t354 * t302 + t357 * t303;
t448 = -qJ(4) - pkin(9);
t380 = qJD(3) * t448;
t523 = -(-qJ(4) * t357 * t358 + pkin(3) * t355) * t403 - t235 - qJD(4) * t354 + t357 * t380;
t236 = t357 * t302 + t354 * t303;
t388 = t358 * t403;
t376 = t354 * t388;
t522 = -qJ(4) * t376 - qJD(4) * t357 - t354 * t380 + t236;
t348 = sin(pkin(11));
t351 = cos(pkin(11));
t323 = t348 * t357 + t351 * t354;
t264 = t323 * t388;
t310 = t323 * qJD(3);
t407 = t264 - t310;
t321 = t348 * t354 - t351 * t357;
t265 = t321 * t388;
t311 = t321 * qJD(3);
t406 = -t265 + t311;
t499 = t348 * t523 - t522 * t351;
t305 = pkin(8) * t388 + t355 * t395;
t262 = pkin(3) * t376 + t305;
t398 = qJD(3) * t354;
t521 = pkin(3) * t398 - t262;
t520 = qJ(5) * t389 - t499;
t519 = -pkin(4) * t407 + t406 * qJ(5) - qJD(5) * t323 + t521;
t338 = qJD(2) + t402;
t399 = qJD(2) * t358;
t386 = t357 * t399;
t397 = qJD(3) * t357;
t253 = t338 * t397 + (-t355 * t398 + t386) * t403;
t387 = t354 * t399;
t254 = -t338 * t398 + (-t355 * t397 - t387) * t403;
t184 = t253 * t348 - t351 * t254;
t331 = qJD(3) - t388;
t347 = sin(pkin(12));
t350 = cos(pkin(12));
t284 = t338 * t357 - t354 * t389;
t285 = t338 * t354 + t357 * t389;
t362 = t284 * t348 + t351 * t285;
t194 = t331 * t347 + t350 * t362;
t353 = sin(qJ(6));
t356 = cos(qJ(6));
t379 = t350 * t331 - t347 * t362;
t120 = t194 * t356 + t353 * t379;
t185 = t253 * t351 + t254 * t348;
t401 = qJD(2) * t349;
t382 = qJD(1) * t401;
t375 = t355 * t382;
t155 = -t185 * t347 + t350 * t375;
t156 = t185 * t350 + t347 * t375;
t50 = -qJD(6) * t120 + t155 * t356 - t156 * t353;
t44 = Ifges(7,6) * t50;
t498 = -t194 * t353 + t356 * t379;
t49 = qJD(6) * t498 + t155 * t353 + t156 * t356;
t45 = Ifges(7,5) * t49;
t10 = Ifges(7,3) * t184 + t44 + t45;
t272 = pkin(9) * t338 + t305;
t298 = (-pkin(2) * t358 - pkin(9) * t355 - pkin(1)) * t349;
t279 = qJD(1) * t298;
t304 = qJD(2) * t360;
t292 = qJD(1) * t304;
t415 = t349 * t355;
t339 = pkin(8) * t415;
t453 = pkin(1) * t358;
t318 = t352 * t453 - t339;
t306 = t318 * qJD(2);
t293 = qJD(1) * t306;
t148 = -t272 * t398 + t279 * t397 + t354 * t292 + t357 * t293;
t105 = qJ(4) * t254 + qJD(4) * t284 + t148;
t214 = t272 * t357 + t279 * t354;
t149 = -qJD(3) * t214 + t357 * t292 - t293 * t354;
t96 = pkin(3) * t375 - qJ(4) * t253 - qJD(4) * t285 + t149;
t42 = t351 * t105 + t348 * t96;
t35 = qJ(5) * t375 + qJD(5) * t331 + t42;
t414 = t349 * t358;
t319 = t352 * t355 * pkin(1) + pkin(8) * t414;
t307 = t319 * qJD(2);
t294 = qJD(1) * t307;
t210 = -t254 * pkin(3) + t294;
t70 = t184 * pkin(4) - t185 * qJ(5) - qJD(5) * t362 + t210;
t18 = -t347 * t35 + t350 * t70;
t19 = t347 * t70 + t350 * t35;
t13 = pkin(10) * t155 + t19;
t224 = -t351 * t284 + t285 * t348;
t271 = -t338 * pkin(2) - t302;
t222 = -t284 * pkin(3) + qJD(4) + t271;
t109 = t224 * pkin(4) - qJ(5) * t362 + t222;
t213 = -t272 * t354 + t357 * t279;
t180 = -qJ(4) * t285 + t213;
t160 = pkin(3) * t331 + t180;
t181 = qJ(4) * t284 + t214;
t412 = t351 * t181;
t93 = t348 * t160 + t412;
t90 = qJ(5) * t331 + t93;
t39 = t350 * t109 - t347 * t90;
t27 = pkin(5) * t224 - pkin(10) * t194 + t39;
t40 = t347 * t109 + t350 * t90;
t31 = pkin(10) * t379 + t40;
t5 = t27 * t356 - t31 * t353;
t7 = pkin(5) * t184 - pkin(10) * t156 + t18;
t1 = qJD(6) * t5 + t13 * t356 + t353 * t7;
t6 = t27 * t353 + t31 * t356;
t2 = -qJD(6) * t6 - t13 * t353 + t356 * t7;
t373 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t434 = t185 * Ifges(5,4);
t479 = t184 / 0.2e1;
t480 = t156 / 0.2e1;
t481 = t155 / 0.2e1;
t518 = 0.2e1 * Ifges(6,5) * t480 + 0.2e1 * Ifges(6,6) * t481 + t373 + t18 * mrSges(6,1) + t210 * mrSges(5,1) + t10 / 0.2e1 + Ifges(6,3) * t479 - t19 * mrSges(6,2) - t434 / 0.2e1;
t237 = t265 * t347 + t350 * t389;
t416 = t347 * t311;
t378 = t237 - t416;
t238 = -t265 * t350 + t347 * t389;
t413 = t350 * t311;
t517 = t238 + t413;
t504 = t347 * t520 + t519 * t350;
t503 = t519 * t347 - t350 * t520;
t516 = -pkin(5) * t407 + pkin(10) * t517 + t504;
t515 = -pkin(10) * t378 + t503;
t500 = t522 * t348 + t351 * t523;
t511 = Ifges(4,3) + Ifges(5,3);
t510 = t184 * Ifges(5,2);
t501 = pkin(4) * t389 - t500;
t420 = t331 * Ifges(5,6);
t427 = t362 * Ifges(5,4);
t142 = -t224 * Ifges(5,2) + t420 + t427;
t508 = t40 * mrSges(6,2) + t6 * mrSges(7,2) + t93 * mrSges(5,3) + t142 / 0.2e1 - t222 * mrSges(5,1) - t39 * mrSges(6,1) - t5 * mrSges(7,1);
t346 = -pkin(3) * t357 - pkin(2);
t255 = pkin(4) * t321 - qJ(5) * t323 + t346;
t332 = t448 * t354;
t333 = t448 * t357;
t270 = t332 * t348 - t333 * t351;
t189 = t350 * t255 - t270 * t347;
t451 = pkin(10) * t350;
t153 = pkin(5) * t321 - t323 * t451 + t189;
t190 = t347 * t255 + t350 * t270;
t417 = t323 * t347;
t167 = -pkin(10) * t417 + t190;
t82 = t153 * t356 - t167 * t353;
t507 = qJD(6) * t82 + t353 * t516 + t356 * t515;
t83 = t153 * t353 + t167 * t356;
t506 = -qJD(6) * t83 - t353 * t515 + t356 * t516;
t221 = qJD(6) + t224;
t431 = t221 * Ifges(7,3);
t439 = t120 * Ifges(7,5);
t440 = t498 * Ifges(7,6);
t46 = t431 + t439 + t440;
t432 = t194 * Ifges(6,5);
t433 = t379 * Ifges(6,6);
t86 = t224 * Ifges(6,3) + t432 + t433;
t505 = t86 + t46;
t502 = pkin(5) * t378 + t501;
t489 = t49 / 0.2e1;
t488 = t50 / 0.2e1;
t343 = pkin(3) * t348 + qJ(5);
t449 = pkin(10) + t343;
t316 = t449 * t347;
t317 = t449 * t350;
t248 = -t316 * t356 - t317 * t353;
t173 = t348 * t181;
t104 = t180 * t351 - t173;
t452 = pkin(3) * t285;
t132 = pkin(4) * t362 + qJ(5) * t224 + t452;
t54 = -t104 * t347 + t350 * t132;
t32 = pkin(5) * t362 + t224 * t451 + t54;
t361 = t347 * t353 - t350 * t356;
t418 = t224 * t347;
t55 = t350 * t104 + t347 * t132;
t38 = pkin(10) * t418 + t55;
t497 = -qJD(5) * t361 + qJD(6) * t248 - t32 * t353 - t356 * t38;
t249 = -t316 * t353 + t317 * t356;
t324 = t347 * t356 + t350 * t353;
t496 = -qJD(5) * t324 - qJD(6) * t249 - t32 * t356 + t353 * t38;
t144 = t324 * t224;
t313 = t324 * qJD(6);
t404 = -t144 - t313;
t145 = t361 * t224;
t312 = t361 * qJD(6);
t405 = -t145 - t312;
t297 = pkin(9) * t352 + t319;
t232 = t357 * t297 + t354 * t298;
t493 = Ifges(4,5) * t253 + Ifges(5,5) * t185 + Ifges(4,6) * t254 - Ifges(5,6) * t184 + t375 * t511;
t424 = t285 * Ifges(4,4);
t208 = t284 * Ifges(4,2) + t331 * Ifges(4,6) + t424;
t280 = Ifges(4,4) * t284;
t209 = t285 * Ifges(4,1) + t331 * Ifges(4,5) + t280;
t363 = t213 * t357 + t214 * t354;
t445 = Ifges(4,4) * t357;
t446 = Ifges(4,4) * t354;
t454 = t357 / 0.2e1;
t457 = t331 / 0.2e1;
t461 = t285 / 0.2e1;
t462 = t284 / 0.2e1;
t492 = -t363 * mrSges(4,3) + t271 * (mrSges(4,1) * t354 + mrSges(4,2) * t357) + (-Ifges(4,2) * t354 + t445) * t462 + (Ifges(4,1) * t357 - t446) * t461 + (Ifges(4,5) * t357 - Ifges(4,6) * t354) * t457 - t354 * t208 / 0.2e1 + t209 * t454;
t491 = Ifges(7,4) * t489 + Ifges(7,2) * t488 + Ifges(7,6) * t479;
t490 = Ifges(7,1) * t489 + Ifges(7,4) * t488 + Ifges(7,5) * t479;
t487 = pkin(1) * mrSges(3,1);
t486 = pkin(1) * mrSges(3,2);
t485 = -t498 / 0.2e1;
t484 = t498 / 0.2e1;
t483 = -t120 / 0.2e1;
t482 = t120 / 0.2e1;
t478 = -t379 / 0.2e1;
t477 = t379 / 0.2e1;
t476 = -t194 / 0.2e1;
t475 = t194 / 0.2e1;
t473 = -t221 / 0.2e1;
t472 = t221 / 0.2e1;
t471 = -t224 / 0.2e1;
t470 = t224 / 0.2e1;
t469 = t362 / 0.2e1;
t314 = t352 * t357 - t354 * t415;
t315 = t352 * t354 + t357 * t415;
t243 = -t351 * t314 + t315 * t348;
t468 = -t243 / 0.2e1;
t244 = t314 * t348 + t315 * t351;
t466 = t244 / 0.2e1;
t465 = t253 / 0.2e1;
t464 = t254 / 0.2e1;
t460 = t314 / 0.2e1;
t459 = t315 / 0.2e1;
t458 = -t331 / 0.2e1;
t456 = t347 / 0.2e1;
t455 = t350 / 0.2e1;
t400 = qJD(2) * t355;
t164 = -qJD(3) * t232 + t357 * t304 - t306 * t354;
t261 = qJD(3) * t314 + t349 * t386;
t385 = t349 * t400;
t114 = pkin(3) * t385 - qJ(4) * t261 - qJD(4) * t315 + t164;
t163 = -t297 * t398 + t298 * t397 + t354 * t304 + t357 * t306;
t260 = -qJD(3) * t315 - t349 * t387;
t126 = qJ(4) * t260 + qJD(4) * t314 + t163;
t59 = t348 * t114 + t351 * t126;
t53 = (qJ(5) * t400 - qJD(5) * t358) * t349 + t59;
t198 = -t351 * t260 + t261 * t348;
t199 = t260 * t348 + t261 * t351;
t229 = -t260 * pkin(3) + t307;
t77 = t198 * pkin(4) - t199 * qJ(5) - t244 * qJD(5) + t229;
t24 = t347 * t77 + t350 * t53;
t447 = Ifges(3,4) * t355;
t444 = Ifges(6,4) * t347;
t443 = Ifges(6,4) * t350;
t442 = Ifges(7,4) * t120;
t441 = Ifges(3,5) * t358;
t436 = t184 * Ifges(5,4);
t435 = t185 * Ifges(5,1);
t430 = t224 * Ifges(5,4);
t429 = t224 * Ifges(5,6);
t428 = t362 * Ifges(5,1);
t426 = t362 * Ifges(5,5);
t425 = t284 * Ifges(4,6);
t423 = t285 * Ifges(4,5);
t422 = t293 * mrSges(3,2);
t421 = t331 * Ifges(5,5);
t419 = t338 * Ifges(3,5);
t231 = -t354 * t297 + t357 * t298;
t188 = -pkin(3) * t414 - t315 * qJ(4) + t231;
t202 = qJ(4) * t314 + t232;
t124 = t348 * t188 + t351 * t202;
t112 = -qJ(5) * t414 + t124;
t296 = t339 + (-pkin(2) - t453) * t352;
t252 = -t314 * pkin(3) + t296;
t140 = t243 * pkin(4) - t244 * qJ(5) + t252;
t64 = t350 * t112 + t347 * t140;
t161 = t237 * t356 - t238 * t353;
t169 = t311 * t324 + t312 * t323;
t411 = t161 - t169;
t162 = t237 * t353 + t238 * t356;
t168 = t311 * t361 - t313 * t323;
t410 = t162 - t168;
t125 = -mrSges(6,1) * t379 + mrSges(6,2) * t194;
t201 = mrSges(5,1) * t331 - mrSges(5,3) * t362;
t409 = t201 - t125;
t408 = -mrSges(3,1) * t338 - mrSges(4,1) * t284 + mrSges(4,2) * t285 + mrSges(3,3) * t389;
t393 = -Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t392 = -Ifges(4,3) / 0.2e1 - Ifges(5,3) / 0.2e1;
t345 = -pkin(3) * t351 - pkin(4);
t16 = -t50 * mrSges(7,1) + t49 * mrSges(7,2);
t23 = -t347 * t53 + t350 * t77;
t106 = t184 * mrSges(5,1) + t185 * mrSges(5,2);
t84 = -t155 * mrSges(6,1) + t156 * mrSges(6,2);
t41 = -t348 * t105 + t351 * t96;
t63 = -t112 * t347 + t350 * t140;
t58 = t114 * t351 - t348 * t126;
t92 = t160 * t351 - t173;
t103 = t180 * t348 + t412;
t123 = t188 * t351 - t348 * t202;
t269 = -t351 * t332 - t333 * t348;
t113 = pkin(4) * t414 - t123;
t372 = mrSges(6,1) * t347 + mrSges(6,2) * t350;
t371 = Ifges(6,1) * t350 - t444;
t370 = -Ifges(6,2) * t347 + t443;
t369 = Ifges(6,5) * t350 - Ifges(6,6) * t347;
t368 = t18 * t350 + t19 * t347;
t367 = -t18 * t347 + t19 * t350;
t366 = -t347 * t39 + t350 * t40;
t220 = t350 * t244 - t347 * t414;
t36 = pkin(5) * t243 - pkin(10) * t220 + t63;
t219 = -t347 * t244 - t350 * t414;
t51 = pkin(10) * t219 + t64;
t14 = -t353 * t51 + t356 * t36;
t15 = t353 * t36 + t356 * t51;
t136 = -mrSges(6,2) * t224 + mrSges(6,3) * t379;
t137 = mrSges(6,1) * t224 - mrSges(6,3) * t194;
t365 = t136 * t350 - t137 * t347;
t364 = t148 * t357 - t149 * t354;
t146 = t219 * t356 - t220 * t353;
t147 = t219 * t353 + t220 * t356;
t89 = -pkin(4) * t331 + qJD(5) - t92;
t56 = -pkin(4) * t385 - t58;
t37 = -pkin(4) * t375 - t41;
t334 = Ifges(3,4) * t388;
t330 = t382 * t441;
t327 = -pkin(5) * t350 + t345;
t301 = -t338 * mrSges(3,2) + mrSges(3,3) * t388;
t267 = Ifges(3,1) * t389 + t334 + t419;
t266 = Ifges(3,6) * t338 + (Ifges(3,2) * t358 + t447) * t403;
t257 = mrSges(4,1) * t331 - mrSges(4,3) * t285;
t256 = -mrSges(4,2) * t331 + mrSges(4,3) * t284;
t242 = t361 * t323;
t241 = t324 * t323;
t234 = pkin(5) * t417 + t269;
t228 = -mrSges(4,2) * t375 + mrSges(4,3) * t254;
t227 = mrSges(4,1) * t375 - mrSges(4,3) * t253;
t207 = t331 * Ifges(4,3) + t423 + t425;
t200 = -mrSges(5,2) * t331 - mrSges(5,3) * t224;
t186 = -mrSges(4,1) * t254 + mrSges(4,2) * t253;
t183 = t199 * t350 + t347 * t385;
t182 = -t199 * t347 + t350 * t385;
t172 = t253 * Ifges(4,1) + t254 * Ifges(4,4) + Ifges(4,5) * t375;
t171 = t253 * Ifges(4,4) + t254 * Ifges(4,2) + Ifges(4,6) * t375;
t159 = mrSges(5,1) * t375 - mrSges(5,3) * t185;
t158 = -mrSges(5,2) * t375 - mrSges(5,3) * t184;
t152 = mrSges(5,1) * t224 + mrSges(5,2) * t362;
t143 = t421 + t428 - t430;
t141 = t331 * Ifges(5,3) + t426 - t429;
t115 = Ifges(7,4) * t498;
t99 = Ifges(5,5) * t375 + t435 - t436;
t98 = Ifges(5,6) * t375 + t434 - t510;
t95 = mrSges(6,1) * t184 - mrSges(6,3) * t156;
t94 = -mrSges(6,2) * t184 + mrSges(6,3) * t155;
t88 = t194 * Ifges(6,1) + Ifges(6,4) * t379 + t224 * Ifges(6,5);
t87 = t194 * Ifges(6,4) + Ifges(6,2) * t379 + t224 * Ifges(6,6);
t81 = mrSges(7,1) * t221 - mrSges(7,3) * t120;
t80 = -mrSges(7,2) * t221 + mrSges(7,3) * t498;
t79 = -pkin(5) * t219 + t113;
t78 = -pkin(5) * t418 + t103;
t72 = -pkin(5) * t379 + t89;
t67 = Ifges(6,1) * t156 + Ifges(6,4) * t155 + Ifges(6,5) * t184;
t66 = t156 * Ifges(6,4) + t155 * Ifges(6,2) + t184 * Ifges(6,6);
t62 = -qJD(6) * t147 + t182 * t356 - t183 * t353;
t61 = qJD(6) * t146 + t182 * t353 + t183 * t356;
t60 = -mrSges(7,1) * t498 + mrSges(7,2) * t120;
t48 = Ifges(7,1) * t120 + Ifges(7,5) * t221 + t115;
t47 = Ifges(7,2) * t498 + Ifges(7,6) * t221 + t442;
t33 = -pkin(5) * t182 + t56;
t30 = -mrSges(7,2) * t184 + mrSges(7,3) * t50;
t29 = mrSges(7,1) * t184 - mrSges(7,3) * t49;
t28 = -pkin(5) * t155 + t37;
t20 = pkin(10) * t182 + t24;
t17 = pkin(5) * t198 - pkin(10) * t183 + t23;
t4 = -qJD(6) * t15 + t17 * t356 - t20 * t353;
t3 = qJD(6) * t14 + t17 * t353 + t20 * t356;
t8 = [m(7) * (t1 * t15 + t14 * t2 + t28 * t79 + t3 * t6 + t33 * t72 + t4 * t5) + m(6) * (t113 * t37 + t18 * t63 + t19 * t64 + t23 * t39 + t24 * t40 + t56 * t89) + m(5) * (t123 * t41 + t124 * t42 + t210 * t252 + t222 * t229 + t58 * t92 + t59 * t93) + m(4) * (t148 * t232 + t149 * t231 + t163 * t214 + t164 * t213 + t271 * t307 + t294 * t296) + m(3) * (t293 * t319 - t294 * t318 - t302 * t307 + t305 * t306) + (Ifges(5,1) * t199 + Ifges(5,5) * t385) * t469 + t185 * (Ifges(5,1) * t244 - Ifges(5,5) * t414) / 0.2e1 + (Ifges(5,4) * t199 + Ifges(5,6) * t385) * t471 - t184 * (Ifges(5,4) * t244 - Ifges(5,6) * t414) / 0.2e1 + (t199 * t222 + t210 * t244 - t385 * t93 + t414 * t42) * mrSges(5,2) + (Ifges(7,5) * t482 + Ifges(7,6) * t484 - Ifges(5,4) * t469 + Ifges(6,3) * t470 - Ifges(5,2) * t471 + Ifges(7,3) * t472 + Ifges(6,5) * t475 + Ifges(6,6) * t477 - Ifges(5,6) * t457 + t505 / 0.2e1 - t508) * t198 + (-t18 * t220 + t182 * t40 - t183 * t39 + t19 * t219) * mrSges(6,3) + ((Ifges(6,3) + Ifges(7,3)) * t479 + Ifges(7,6) * t488 + Ifges(7,5) * t489 + t510 / 0.2e1 - t42 * mrSges(5,3) + t518) * t243 + (t1 * t146 - t147 * t2 - t5 * t61 + t6 * t62) * mrSges(7,3) + t199 * t143 / 0.2e1 + t59 * t200 + t182 * t87 / 0.2e1 + t89 * (-mrSges(6,1) * t182 + mrSges(6,2) * t183) + t183 * t88 / 0.2e1 + (t293 * t358 + t294 * t355 + (-t302 * t358 - t305 * t355) * qJD(2)) * t349 * mrSges(3,3) + (Ifges(4,4) * t261 + Ifges(4,2) * t260 + Ifges(4,6) * t385) * t462 + (Ifges(4,4) * t315 + Ifges(4,2) * t314 - Ifges(4,6) * t414) * t464 + (Ifges(4,1) * t315 + Ifges(4,4) * t314 - Ifges(4,5) * t414) * t465 + t99 * t466 + t98 * t468 + (Ifges(6,5) * t220 + Ifges(7,5) * t147 + Ifges(6,6) * t219 + Ifges(7,6) * t146) * t479 + t172 * t459 + t171 * t460 + (Ifges(4,1) * t261 + Ifges(4,4) * t260 + Ifges(4,5) * t385) * t461 + t219 * t66 / 0.2e1 + t58 * t201 + (Ifges(7,5) * t61 + Ifges(7,6) * t62) * t472 + t37 * (-mrSges(6,1) * t219 + mrSges(6,2) * t220) + t220 * t67 / 0.2e1 + (Ifges(6,4) * t220 + Ifges(6,2) * t219) * t481 + (Ifges(6,4) * t183 + Ifges(6,2) * t182) * t477 + (Ifges(7,4) * t61 + Ifges(7,2) * t62) * t484 + (Ifges(7,4) * t147 + Ifges(7,2) * t146) * t488 + (Ifges(6,5) * t183 + Ifges(6,6) * t182) * t470 + (t358 * t267 + t338 * (-Ifges(3,6) * t355 + t441) + (t207 + t141) * t355) * t401 / 0.2e1 + t124 * t158 + t123 * t159 + t28 * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) + t24 * t136 + t23 * t137 + t56 * t125 + t113 * t84 + t64 * t94 + t63 * t95 + t3 * t80 + t4 * t81 + (Ifges(7,1) * t61 + Ifges(7,4) * t62) * t482 + (Ifges(7,1) * t147 + Ifges(7,4) * t146) * t489 + t79 * t16 + t72 * (-mrSges(7,1) * t62 + mrSges(7,2) * t61) + t33 * t60 + t61 * t48 / 0.2e1 + t62 * t47 / 0.2e1 + t14 * t29 + t15 * t30 + (Ifges(6,1) * t220 + Ifges(6,4) * t219) * t480 + (Ifges(6,1) * t183 + Ifges(6,4) * t182) * t475 + ((Ifges(3,5) * t352 / 0.2e1 - t318 * mrSges(3,3) + (-0.2e1 * t486 + 0.3e1 / 0.2e1 * Ifges(3,4) * t358) * t349) * t358 + (Ifges(5,5) * t466 + Ifges(5,6) * t468 + Ifges(4,5) * t459 + Ifges(4,6) * t460 - Ifges(3,6) * t352 - t319 * mrSges(3,3) + (-0.2e1 * t487 - 0.3e1 / 0.2e1 * t447) * t349 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + t392) * t414) * t355) * t382 + (t330 / 0.2e1 - t294 * mrSges(3,1) - t422) * t352 + t41 * (-mrSges(5,1) * t414 - t244 * mrSges(5,3)) + t149 * (-mrSges(4,1) * t414 - t315 * mrSges(4,3)) + t148 * (mrSges(4,2) * t414 + t314 * mrSges(4,3)) + t408 * t307 - t266 * t385 / 0.2e1 + t214 * (-mrSges(4,2) * t385 + mrSges(4,3) * t260) + t92 * (mrSges(5,1) * t385 - mrSges(5,3) * t199) + t213 * (mrSges(4,1) * t385 - mrSges(4,3) * t261) + t147 * t490 + t146 * t491 + t229 * t152 + t231 * t227 + t232 * t228 - t493 * t414 / 0.2e1 + t252 * t106 + t163 * t256 + t164 * t257 + t260 * t208 / 0.2e1 + t261 * t209 / 0.2e1 + t271 * (-mrSges(4,1) * t260 + mrSges(4,2) * t261) + (Ifges(4,5) * t261 + Ifges(5,5) * t199 + Ifges(4,6) * t260 + t385 * t511) * t457 + t296 * t186 + t306 * t301 + t294 * (-mrSges(4,1) * t314 + mrSges(4,2) * t315); -t422 + (t241 * t28 - t407 * t5 + t411 * t72) * mrSges(7,1) + (-Ifges(7,5) * t242 - Ifges(7,6) * t241) * t479 + (-Ifges(7,4) * t242 - Ifges(7,2) * t241) * t488 + (t84 - t159) * t269 + (-mrSges(6,1) * t407 + mrSges(6,3) * t517) * t39 + (mrSges(6,1) * t378 - mrSges(6,2) * t517) * t89 + (-mrSges(4,1) * t357 + mrSges(4,2) * t354 - mrSges(3,1)) * t294 + (t45 / 0.2e1 + t44 / 0.2e1 - t98 / 0.2e1 + (Ifges(7,3) / 0.2e1 - t393) * t184 + t518) * t321 + (-Ifges(5,4) * t311 + Ifges(6,5) * t238 - Ifges(5,2) * t310 + Ifges(6,6) * t237 + Ifges(6,3) * t264) * t471 + (-t227 * t354 + t228 * t357) * pkin(9) - pkin(2) * t186 + t189 * t95 + t190 * t94 + (-Ifges(7,1) * t242 - Ifges(7,4) * t241) * t489 + (-t1 * t241 + t2 * t242 + t410 * t5 - t411 * t6) * mrSges(7,3) + (-t242 * t28 + t407 * t6 - t410 * t72) * mrSges(7,2) - t362 * (-Ifges(5,1) * t265 - Ifges(5,4) * t264) / 0.2e1 + (-Ifges(5,5) * t265 - Ifges(5,6) * t264) * t458 + (t265 / 0.2e1 - t311 / 0.2e1) * t143 + (Ifges(6,6) * t310 - t311 * t370) * t477 + (-Ifges(5,1) * t311 - Ifges(5,4) * t310) * t469 + (Ifges(6,5) * t310 - t311 * t371) * t475 + (t168 / 0.2e1 - t162 / 0.2e1) * t48 + (Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t264) * t478 + (Ifges(7,1) * t168 + Ifges(7,4) * t169 + Ifges(7,5) * t310) * t482 + (Ifges(7,1) * t162 + Ifges(7,4) * t161 + Ifges(7,5) * t264) * t483 + (Ifges(7,4) * t168 + Ifges(7,2) * t169 + Ifges(7,6) * t310) * t484 + (Ifges(7,4) * t162 + Ifges(7,2) * t161 + Ifges(7,6) * t264) * t485 + (Ifges(4,2) * t357 + t446) * t464 + (Ifges(4,1) * t354 + t445) * t465 + (Ifges(7,5) * t168 + Ifges(7,6) * t169 + Ifges(7,3) * t310) * t472 + (Ifges(7,5) * t162 + Ifges(7,6) * t161 + Ifges(7,3) * t264) * t473 + (Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t264) * t476 + t171 * t454 + (t210 * t346 + t222 * t521 - t269 * t41 + t270 * t42 + t499 * t93 + t500 * t92) * m(5) + (-Ifges(5,5) * t311 - Ifges(5,6) * t310) * t457 + (-pkin(2) * t294 - t213 * t235 - t214 * t236 - t271 * t305 + (-qJD(3) * t363 + t364) * pkin(9)) * m(4) + (-Ifges(5,4) * t265 - Ifges(5,2) * t264 + Ifges(6,3) * t310 - t311 * t369) * t470 + (-t238 / 0.2e1 - t413 / 0.2e1) * t88 + (-t237 / 0.2e1 + t416 / 0.2e1) * t87 + t82 * t29 + t83 * t30 + (t99 / 0.2e1 + t435 / 0.2e1 - t436 / 0.2e1 + t210 * mrSges(5,2) + t37 * t372 + t369 * t479 + t371 * t480 + t370 * t481 - t347 * t66 / 0.2e1 + t67 * t455 - t368 * mrSges(6,3)) * t323 - t408 * t305 + (-t321 * t42 - t323 * t41 + t406 * t92 + t407 * t93) * mrSges(5,3) + (mrSges(6,2) * t407 - mrSges(6,3) * t378) * t40 + (-mrSges(5,1) * t407 - mrSges(5,2) * t406) * t222 + t330 - t242 * t490 - t241 * t491 + (t354 * pkin(3) * t152 + (-t256 * t354 - t257 * t357) * pkin(9) + t492) * qJD(3) + ((-t419 / 0.2e1 - t334 / 0.2e1 - t267 / 0.2e1 + t302 * mrSges(3,3) + t403 * t486 - t492) * t358 + (t214 * mrSges(4,2) - t213 * mrSges(4,1) - t425 / 0.2e1 - t423 / 0.2e1 + t266 / 0.2e1 - t207 / 0.2e1 - t141 / 0.2e1 - t92 * mrSges(5,1) + t93 * mrSges(5,2) + t429 / 0.2e1 - t426 / 0.2e1 + t305 * mrSges(3,3) + t392 * t331 + (t487 + t447 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t358) * t403 + (t338 / 0.2e1 - qJD(2)) * Ifges(3,6) + (Ifges(4,5) * t354 + Ifges(5,5) * t323 + Ifges(4,6) * t357 - Ifges(5,6) * t321) * qJD(2) / 0.2e1) * t355) * t403 + t234 * t16 + t499 * t200 + t500 * t201 + t501 * t125 + (t169 / 0.2e1 - t161 / 0.2e1) * t47 - t236 * t256 - t235 * t257 + t502 * t60 + t503 * t136 + t504 * t137 + (t18 * t189 + t19 * t190 + t269 * t37 + t39 * t504 + t40 * t503 + t501 * t89) * m(6) - t262 * t152 + t270 * t158 - t302 * t301 + (t142 - t505) * (t264 / 0.2e1 - t310 / 0.2e1) + t506 * t81 + t507 * t80 + (t1 * t83 + t2 * t82 + t234 * t28 + t506 * t5 + t502 * t72 + t507 * t6) * m(7) + t346 * t106 + t354 * t172 / 0.2e1 + t364 * mrSges(4,3); (t420 / 0.2e1 - t86 / 0.2e1 - t440 / 0.2e1 - t439 / 0.2e1 - t432 / 0.2e1 - t433 / 0.2e1 - t431 / 0.2e1 + t508 - t46 / 0.2e1 + t393 * t224 + t427 / 0.2e1) * t362 - (-Ifges(4,2) * t285 + t209 + t280) * t284 / 0.2e1 + (-t152 * t285 + t158 * t348 + t159 * t351) * pkin(3) + t493 + (qJD(5) * t366 - t103 * t89 + t343 * t367 + t345 * t37 - t39 * t54 - t40 * t55) * m(6) - (-t143 / 0.2e1 - t222 * mrSges(5,2) + t430 / 0.2e1 - t428 / 0.2e1 + t92 * mrSges(5,3) - t421 / 0.2e1 + t370 * t478 + t371 * t476 - t89 * t372 + t369 * t471 + t87 * t456 - t350 * t88 / 0.2e1 + (t347 * t40 + t350 * t39) * mrSges(6,3)) * t224 + (Ifges(7,4) * t324 - Ifges(7,2) * t361) * t488 + (Ifges(7,1) * t324 - Ifges(7,4) * t361) * t489 + (Ifges(6,5) * t347 + Ifges(7,5) * t324 + Ifges(6,6) * t350 - Ifges(7,6) * t361) * t479 + (-t1 * t361 - t2 * t324 + t404 * t6 - t405 * t5) * mrSges(7,3) + t28 * (mrSges(7,1) * t361 + mrSges(7,2) * t324) - t361 * t491 + (-t144 / 0.2e1 - t313 / 0.2e1) * t47 + (-t145 / 0.2e1 - t312 / 0.2e1) * t48 - t104 * t200 + (t103 * t92 - t104 * t93 - t222 * t452 + (t348 * t42 + t351 * t41) * pkin(3)) * m(5) + (Ifges(6,1) * t347 + t443) * t480 + (Ifges(6,2) * t350 + t444) * t481 + (Ifges(7,1) * t145 + Ifges(7,4) * t144) * t483 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t485 + t324 * t490 + t208 * t461 + (Ifges(7,5) * t145 + Ifges(7,6) * t144) * t473 + (t213 * t284 + t214 * t285) * mrSges(4,3) + t66 * t455 + t67 * t456 + (Ifges(4,5) * t284 - Ifges(4,6) * t285) * t458 + (-Ifges(7,1) * t312 - Ifges(7,4) * t313) * t482 + (-Ifges(7,4) * t312 - Ifges(7,2) * t313) * t484 + (-Ifges(7,5) * t312 - Ifges(7,6) * t313) * t472 + (-t347 * t95 + t350 * t94) * t343 + t149 * mrSges(4,1) - t148 * mrSges(4,2) - t55 * t136 - t54 * t137 - t78 * t60 + t41 * mrSges(5,1) - t42 * mrSges(5,2) - t285 * (Ifges(4,1) * t284 - t424) / 0.2e1 + t409 * t103 + (-mrSges(7,1) * t404 + mrSges(7,2) * t405) * t72 + t248 * t29 + t249 * t30 - t213 * t256 + t214 * t257 - t271 * (mrSges(4,1) * t285 + mrSges(4,2) * t284) + t496 * t81 + t497 * t80 + (t1 * t249 + t2 * t248 + t28 * t327 + t496 * t5 + t497 * t6 - t72 * t78) * m(7) + t327 * t16 + t345 * t84 + t37 * (-mrSges(6,1) * t350 + mrSges(6,2) * t347) + t365 * qJD(5) + t367 * mrSges(6,3); -t361 * t29 + t324 * t30 + t347 * t94 + t350 * t95 + t404 * t81 + t405 * t80 - (-t200 - t365) * t224 + (-t60 + t409) * t362 + t106 + (t1 * t324 - t2 * t361 - t362 * t72 + t404 * t5 + t405 * t6) * m(7) + (t224 * t366 - t362 * t89 + t368) * m(6) + (t224 * t93 + t362 * t92 + t210) * m(5); t120 * t81 - t498 * t80 - t379 * t136 + t194 * t137 + t16 + t84 + (t120 * t5 - t498 * t6 + t28) * m(7) + (t194 * t39 - t379 * t40 + t37) * m(6); -t72 * (mrSges(7,1) * t120 + mrSges(7,2) * t498) + (Ifges(7,1) * t498 - t442) * t483 + t47 * t482 + (Ifges(7,5) * t498 - Ifges(7,6) * t120) * t473 - t5 * t80 + t6 * t81 + (t120 * t6 + t498 * t5) * mrSges(7,3) + t373 + t10 + (-Ifges(7,2) * t120 + t115 + t48) * t485;];
tauc  = t8(:);
