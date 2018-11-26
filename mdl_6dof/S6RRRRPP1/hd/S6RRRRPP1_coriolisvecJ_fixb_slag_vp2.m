% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2018-11-23 18:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:04:18
% EndTime: 2018-11-23 18:04:29
% DurationCPUTime: 11.60s
% Computational Cost: add. (14535->637), mult. (35901->811), div. (0->0), fcn. (25619->8), ass. (0->318)
t343 = sin(qJ(3));
t344 = sin(qJ(2));
t346 = cos(qJ(3));
t347 = cos(qJ(2));
t312 = t343 * t344 - t346 * t347;
t340 = qJD(2) + qJD(3);
t269 = t340 * t312;
t254 = t269 * qJD(1);
t313 = t343 * t347 + t346 * t344;
t301 = t313 * qJD(1);
t342 = sin(qJ(4));
t345 = cos(qJ(4));
t277 = -t301 * t342 + t340 * t345;
t183 = qJD(4) * t277 - t254 * t345;
t278 = t301 * t345 + t340 * t342;
t184 = -qJD(4) * t278 + t254 * t342;
t341 = sin(pkin(10));
t420 = cos(pkin(10));
t96 = t183 * t420 + t341 * t184;
t484 = t96 / 0.2e1;
t563 = 0.2e1 * t484;
t300 = t312 * qJD(1);
t411 = t300 * t342;
t562 = qJ(5) * t411 - t345 * qJD(5);
t482 = -pkin(8) - pkin(7);
t326 = t482 * t347;
t316 = qJD(1) * t326;
t302 = t343 * t316;
t325 = t482 * t344;
t315 = qJD(1) * t325;
t306 = qJD(2) * pkin(2) + t315;
t262 = t306 * t346 + t302;
t232 = -pkin(3) * t340 - t262;
t180 = -pkin(4) * t277 + qJD(5) + t232;
t294 = qJD(4) + t300;
t334 = -pkin(2) * t347 - pkin(1);
t324 = qJD(1) * t334;
t220 = t300 * pkin(3) - t301 * pkin(9) + t324;
t402 = t346 * t316;
t263 = t343 * t306 - t402;
t233 = pkin(9) * t340 + t263;
t142 = t345 * t220 - t233 * t342;
t116 = -qJ(5) * t278 + t142;
t107 = pkin(4) * t294 + t116;
t143 = t220 * t342 + t233 * t345;
t117 = qJ(5) * t277 + t143;
t113 = t420 * t117;
t37 = t341 * t107 + t113;
t32 = qJ(6) * t294 + t37;
t459 = t294 / 0.2e1;
t460 = -t294 / 0.2e1;
t356 = t341 * t277 + t278 * t420;
t467 = t356 / 0.2e1;
t468 = -t356 / 0.2e1;
t194 = -t420 * t277 + t278 * t341;
t471 = -t194 / 0.2e1;
t56 = pkin(5) * t194 - qJ(6) * t356 + t180;
t561 = mrSges(6,3) * t37 - mrSges(6,1) * t180 - mrSges(7,1) * t56 + mrSges(7,2) * t32 + Ifges(6,4) * t467 + Ifges(7,5) * t468 + Ifges(6,6) * t459 + Ifges(7,6) * t460 + (Ifges(6,2) + Ifges(7,3)) * t471;
t95 = t183 * t341 - t184 * t420;
t486 = -t95 / 0.2e1;
t270 = t340 * t313;
t255 = t270 * qJD(1);
t464 = t255 / 0.2e1;
t517 = Ifges(6,1) + Ifges(7,1);
t516 = Ifges(6,4) - Ifges(7,5);
t515 = Ifges(7,4) + Ifges(6,5);
t541 = Ifges(7,2) + Ifges(6,3);
t560 = Ifges(5,3) + t541;
t258 = pkin(3) * t301 + pkin(9) * t300;
t399 = qJD(1) * t344;
t227 = pkin(2) * t399 + t258;
t266 = t315 * t346 + t302;
t160 = t345 * t227 - t266 * t342;
t339 = t345 * qJ(5);
t372 = t301 * pkin(4) + t300 * t339;
t331 = pkin(2) * t343 + pkin(9);
t401 = -qJ(5) - t331;
t374 = qJD(4) * t401;
t432 = pkin(2) * qJD(3);
t390 = t346 * t432;
t559 = -t160 - t372 + (-qJD(5) - t390) * t342 + t345 * t374;
t164 = t345 * t258 - t262 * t342;
t441 = -qJ(5) - pkin(9);
t379 = qJD(4) * t441;
t558 = -qJD(5) * t342 + t345 * t379 - t164 - t372;
t161 = t342 * t227 + t345 * t266;
t557 = -t342 * t374 - t345 * t390 + t161 + t562;
t165 = t342 * t258 + t345 * t262;
t556 = -t342 * t379 + t165 + t562;
t470 = t194 / 0.2e1;
t514 = Ifges(6,6) - Ifges(7,6);
t555 = -Ifges(6,2) * t471 + Ifges(7,3) * t470 - t459 * t514 - t467 * t516 - t561;
t554 = t340 * Ifges(4,6) / 0.2e1;
t285 = pkin(4) * t411;
t395 = qJD(4) * t342;
t336 = pkin(4) * t395;
t553 = t285 + t336;
t101 = Ifges(7,1) * t356 + Ifges(7,4) * t294 + Ifges(7,5) * t194;
t102 = Ifges(6,1) * t356 - Ifges(6,4) * t194 + Ifges(6,5) * t294;
t407 = t341 * t117;
t36 = t107 * t420 - t407;
t524 = t180 * mrSges(6,2) - t36 * mrSges(6,3) - t56 * mrSges(7,3);
t552 = t524 + t102 / 0.2e1 + t101 / 0.2e1;
t424 = t340 * Ifges(4,5);
t551 = t324 * mrSges(4,2) + t424 / 0.2e1;
t31 = -t294 * pkin(5) + qJD(6) - t36;
t449 = mrSges(7,2) * t31;
t550 = Ifges(6,4) * t470 + Ifges(7,5) * t471 + t517 * t468 - t449;
t425 = t301 * Ifges(4,4);
t549 = t554 + t425 / 0.2e1 - t300 * Ifges(4,2) / 0.2e1;
t548 = -Ifges(6,2) * t470 + Ifges(7,3) * t471 - t468 * t516 + t561;
t546 = Ifges(6,4) * t471 + Ifges(7,5) * t470 + t515 * t459 + t517 * t467 + t449;
t397 = qJD(2) * t344;
t391 = pkin(2) * t397;
t531 = qJD(1) * t391;
t151 = pkin(3) * t255 + pkin(9) * t254 + t531;
t386 = qJD(2) * t482;
t373 = qJD(1) * t386;
t307 = t344 * t373;
t361 = t347 * t373;
t170 = qJD(3) * t262 + t346 * t307 + t343 * t361;
t40 = -qJD(4) * t143 + t345 * t151 - t170 * t342;
t14 = pkin(4) * t255 - qJ(5) * t183 - qJD(5) * t278 + t40;
t394 = qJD(4) * t345;
t39 = t342 * t151 + t345 * t170 + t220 * t394 - t233 * t395;
t17 = qJ(5) * t184 + qJD(5) * t277 + t39;
t6 = t341 * t14 + t420 * t17;
t1 = qJ(6) * t255 + qJD(6) * t294 + t6;
t171 = qJD(3) * t263 + t307 * t343 - t346 * t361;
t103 = -pkin(4) * t184 + t171;
t12 = pkin(5) * t95 - qJ(6) * t96 - qJD(6) * t356 + t103;
t485 = t95 / 0.2e1;
t520 = -Ifges(6,6) / 0.2e1;
t545 = t516 * t563 - mrSges(6,1) * t103 - mrSges(7,1) * t12 + mrSges(7,2) * t1 + mrSges(6,3) * t6 - 0.2e1 * Ifges(7,3) * t485 - t255 * t520 + (-t485 + t486) * Ifges(6,2) + (t514 - Ifges(7,6)) * t464;
t5 = t14 * t420 - t341 * t17;
t3 = -t255 * pkin(5) - t5;
t544 = mrSges(6,2) * t103 + mrSges(7,2) * t3 - mrSges(6,3) * t5 - mrSges(7,3) * t12 + Ifges(7,5) * t485 + 0.2e1 * t464 * t515 + t517 * t563 + (Ifges(6,4) + t516) * t486;
t453 = -t342 / 0.2e1;
t508 = t341 * t557 + t420 * t559;
t507 = t341 * t559 - t420 * t557;
t506 = t341 * t556 + t420 * t558;
t505 = t341 * t558 - t420 * t556;
t537 = t142 * t345;
t378 = t420 * t342;
t310 = t341 * t345 + t378;
t216 = t310 * t300;
t377 = t420 * t345;
t406 = t341 * t342;
t355 = t377 - t406;
t217 = t355 * t300;
t296 = t310 * qJD(4);
t297 = t355 * qJD(4);
t536 = -qJD(6) * t310 + t553 + (-t217 - t297) * qJ(6) + (t216 + t296) * pkin(5);
t265 = t315 * t343 - t402;
t535 = t343 * t432 - t265;
t435 = Ifges(5,4) * t278;
t167 = t277 * Ifges(5,2) + t294 * Ifges(5,6) + t435;
t370 = mrSges(5,1) * t342 + mrSges(5,2) * t345;
t534 = t167 * t453 + t232 * t370;
t533 = -t40 * t342 + t345 * t39;
t365 = Ifges(5,5) * t345 - Ifges(5,6) * t342;
t433 = Ifges(5,4) * t345;
t367 = -Ifges(5,2) * t342 + t433;
t434 = Ifges(5,4) * t342;
t369 = Ifges(5,1) * t345 - t434;
t274 = Ifges(5,4) * t277;
t168 = Ifges(5,1) * t278 + Ifges(5,5) * t294 + t274;
t403 = t345 * t168;
t461 = t278 / 0.2e1;
t526 = t365 * t459 + t369 * t461 + t277 * t367 / 0.2e1 + t403 / 0.2e1 + t534;
t521 = -t324 * mrSges(4,1) - t142 * mrSges(5,1) - t36 * mrSges(6,1) + t31 * mrSges(7,1) + t143 * mrSges(5,2) + t37 * mrSges(6,2) - t32 * mrSges(7,3) + t549;
t513 = mrSges(4,2) * t301;
t290 = t301 * qJ(6);
t512 = -t290 + t507;
t444 = t301 * pkin(5);
t511 = t444 - t508;
t510 = -t290 + t505;
t509 = t444 - t506;
t504 = -t263 + t536;
t42 = t116 * t420 - t407;
t503 = qJD(6) - t42;
t501 = t535 + t536;
t154 = mrSges(6,1) * t294 - mrSges(6,3) * t356;
t155 = -mrSges(7,1) * t294 + mrSges(7,2) * t356;
t500 = t154 - t155;
t499 = Ifges(5,5) * t183 + Ifges(5,6) * t184;
t498 = t535 + t553;
t261 = t312 * pkin(3) - t313 * pkin(9) + t334;
t280 = t325 * t343 - t326 * t346;
t273 = t345 * t280;
t186 = t342 * t261 + t273;
t497 = t346 * t325 + t326 * t343;
t496 = -t142 * t342 + t143 * t345;
t105 = -mrSges(5,1) * t184 + mrSges(5,2) * t183;
t495 = m(5) * t171 + t105;
t494 = t255 * t560 - t514 * t95 + t515 * t96 + t499;
t493 = t40 * mrSges(5,1) + t5 * mrSges(6,1) - t3 * mrSges(7,1) - t39 * mrSges(5,2) - t6 * mrSges(6,2) + t1 * mrSges(7,3);
t490 = m(4) / 0.2e1;
t489 = Ifges(5,3) / 0.2e1;
t481 = pkin(1) * mrSges(3,1);
t480 = pkin(1) * mrSges(3,2);
t473 = t183 / 0.2e1;
t472 = t184 / 0.2e1;
t463 = -t277 / 0.2e1;
t462 = -t278 / 0.2e1;
t452 = t345 / 0.2e1;
t451 = m(4) * t324;
t447 = pkin(2) * t346;
t446 = pkin(4) * t278;
t445 = pkin(4) * t341;
t360 = qJ(5) * t269 - qJD(5) * t313;
t177 = pkin(3) * t270 + pkin(9) * t269 + t391;
t318 = t344 * t386;
t319 = t347 * t386;
t197 = qJD(3) * t497 + t318 * t346 + t319 * t343;
t375 = t345 * t177 - t197 * t342;
t27 = pkin(4) * t270 + t360 * t345 + (-t273 + (qJ(5) * t313 - t261) * t342) * qJD(4) + t375;
t384 = t313 * t394;
t388 = t342 * t177 + t345 * t197 + t261 * t394;
t35 = -qJ(5) * t384 + (-qJD(4) * t280 + t360) * t342 + t388;
t10 = t341 * t27 + t420 * t35;
t63 = -mrSges(7,2) * t95 + mrSges(7,3) * t255;
t64 = -mrSges(6,2) * t255 - mrSges(6,3) * t95;
t440 = t63 + t64;
t65 = mrSges(6,1) * t255 - mrSges(6,3) * t96;
t66 = -t255 * mrSges(7,1) + t96 * mrSges(7,2);
t439 = t66 - t65;
t438 = mrSges(4,3) * t300;
t437 = mrSges(5,3) * t277;
t436 = Ifges(3,4) * t344;
t291 = Ifges(4,4) * t300;
t431 = t143 * mrSges(5,3);
t430 = t277 * Ifges(5,6);
t429 = t278 * Ifges(5,5);
t427 = t301 * mrSges(4,3);
t426 = t301 * Ifges(4,1);
t419 = Ifges(3,5) * qJD(2);
t418 = Ifges(3,6) * qJD(2);
t417 = qJD(2) * mrSges(3,1);
t416 = qJD(2) * mrSges(3,2);
t413 = t171 * t497;
t409 = t313 * t342;
t185 = t345 * t261 - t280 * t342;
t135 = pkin(4) * t312 - t313 * t339 + t185;
t149 = -qJ(5) * t409 + t186;
t58 = t341 * t135 + t420 * t149;
t400 = -mrSges(4,1) * t340 - mrSges(5,1) * t277 + mrSges(5,2) * t278 + t427;
t398 = qJD(1) * t347;
t396 = qJD(4) * t313;
t333 = -pkin(4) * t345 - pkin(3);
t385 = t420 * pkin(4);
t383 = t419 / 0.2e1;
t382 = -t418 / 0.2e1;
t30 = t95 * mrSges(6,1) + t96 * mrSges(6,2);
t29 = t95 * mrSges(7,1) - t96 * mrSges(7,3);
t376 = t401 * t342;
t221 = pkin(4) * t409 - t497;
t371 = mrSges(5,1) * t345 - mrSges(5,2) * t342;
t368 = Ifges(5,1) * t342 + t433;
t366 = Ifges(5,2) * t345 + t434;
t364 = Ifges(5,5) * t342 + Ifges(5,6) * t345;
t124 = mrSges(5,1) * t255 - mrSges(5,3) * t183;
t125 = -mrSges(5,2) * t255 + mrSges(5,3) * t184;
t363 = -t124 * t342 + t125 * t345;
t362 = -t143 * t342 - t537;
t9 = t27 * t420 - t341 * t35;
t357 = t362 * mrSges(5,3);
t57 = t135 * t420 - t341 * t149;
t253 = -pkin(5) * t355 - qJ(6) * t310 + t333;
t198 = qJD(3) * t280 + t318 * t343 - t346 * t319;
t123 = t198 + (-t269 * t342 + t384) * pkin(4);
t349 = m(5) * (qJD(4) * t362 + t533);
t166 = t294 * Ifges(5,3) + t429 + t430;
t229 = -t291 + t424 + t426;
t69 = t183 * Ifges(5,4) + t184 * Ifges(5,2) + t255 * Ifges(5,6);
t70 = t183 * Ifges(5,1) + t184 * Ifges(5,4) + t255 * Ifges(5,5);
t98 = Ifges(6,5) * t356 - t194 * Ifges(6,6) + t294 * Ifges(6,3);
t99 = Ifges(7,4) * t356 + t294 * Ifges(7,2) + t194 * Ifges(7,6);
t348 = (-t143 * t411 - t300 * t537 + t533) * mrSges(5,3) - (-Ifges(4,1) * t300 + t166 - t425 + t98 + t99) * t301 / 0.2e1 + (-Ifges(4,2) * t301 + t229 - t291 + t403) * t300 / 0.2e1 + t364 * t464 + t366 * t472 + t368 * t473 + (-t367 * t463 - t369 * t462 + t534 + t551) * t300 - t548 * t216 - (-t524 + t550) * t217 + (-t371 - mrSges(4,1)) * t171 - t262 * t438 + t69 * t452 + (t514 * t216 - t515 * t217 - t300 * t365 + t301 * t560) * t460 + (Ifges(5,5) * t462 + Ifges(5,6) * t463 + Ifges(6,6) * t470 + Ifges(7,6) * t471 + t515 * t468 + t521 + t554) * t301 + t545 * t355 + t555 * t296 + (t524 + t546) * t297 + t526 * qJD(4) - t170 * mrSges(4,2) - Ifges(4,5) * t254 - Ifges(4,6) * t255 + t342 * t70 / 0.2e1 + (t102 + t101) * (t217 / 0.2e1 + t297 / 0.2e1) + t544 * t310;
t335 = Ifges(3,4) * t398;
t329 = -t385 - pkin(5);
t327 = qJ(6) + t445;
t323 = pkin(9) * t345 + t339;
t322 = mrSges(3,3) * t398 - t416;
t321 = -mrSges(3,3) * t399 + t417;
t320 = t333 - t447;
t308 = t331 * t345 + t339;
t299 = Ifges(3,1) * t399 + t335 + t419;
t298 = t418 + (t347 * Ifges(3,2) + t436) * qJD(1);
t281 = -mrSges(4,2) * t340 - t438;
t276 = t323 * t420 + t406 * t441;
t275 = t323 * t341 - t378 * t441;
t257 = mrSges(4,1) * t300 + t513;
t250 = t308 * t420 + t341 * t376;
t249 = t308 * t341 - t376 * t420;
t236 = t355 * t313;
t235 = t310 * t313;
t231 = t253 - t447;
t213 = mrSges(5,1) * t294 - mrSges(5,3) * t278;
t212 = -mrSges(5,2) * t294 + t437;
t205 = -t285 + t263;
t153 = -mrSges(6,2) * t294 - mrSges(6,3) * t194;
t152 = -mrSges(7,2) * t194 + mrSges(7,3) * t294;
t121 = -t269 * t355 - t310 * t396;
t120 = t313 * t341 * t395 + t269 * t310 - t377 * t396;
t112 = pkin(5) * t235 - qJ(6) * t236 + t221;
t111 = mrSges(6,1) * t194 + mrSges(6,2) * t356;
t110 = mrSges(7,1) * t194 - mrSges(7,3) * t356;
t86 = pkin(5) * t356 + qJ(6) * t194 + t446;
t55 = -t312 * pkin(5) - t57;
t54 = -qJD(4) * t186 + t375;
t53 = -t280 * t395 + t388;
t52 = qJ(6) * t312 + t58;
t41 = t116 * t341 + t113;
t15 = -pkin(5) * t120 - qJ(6) * t121 - qJD(6) * t236 + t123;
t8 = -t270 * pkin(5) - t9;
t7 = qJ(6) * t270 + qJD(6) * t312 + t10;
t2 = [-(mrSges(4,2) * t334 - mrSges(4,3) * t497) * t254 - t497 * t105 + (mrSges(4,1) * t531 - mrSges(4,3) * t170 + Ifges(4,4) * t254 + Ifges(6,6) * t486 + Ifges(7,6) * t485 + t541 * t464 + t515 * t484 + t493) * t312 + (t70 * t452 + t69 * t453 - Ifges(4,1) * t254 + t369 * t473 + t367 * t472 + t365 * t464 + (mrSges(4,3) + t370) * t171 + (-t342 * t39 - t345 * t40) * mrSges(5,3) + (t168 * t453 - t345 * t167 / 0.2e1 + t232 * t371 + t366 * t463 + t368 * t462 + t364 * t460 - t496 * mrSges(5,3)) * qJD(4)) * t313 + (t262 * t269 - t263 * t270) * mrSges(4,3) + m(4) * (t170 * t280 + t197 * t263 - t198 * t262 - t413) + m(5) * (t142 * t54 + t143 * t53 + t185 * t40 + t186 * t39 + t198 * t232 - t413) + (t166 / 0.2e1 + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t356 + t429 / 0.2e1 + t430 / 0.2e1 + t99 / 0.2e1 - t521 + t98 / 0.2e1 + (Ifges(7,6) / 0.2e1 + t520) * t194 + (t489 + Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t294 - t549) * t270 + m(7) * (t1 * t52 + t112 * t12 + t15 * t56 + t3 * t55 + t31 * t8 + t32 * t7) + m(6) * (t10 * t37 + t103 * t221 + t123 * t180 + t36 * t9 + t5 * t57 + t58 * t6) + t544 * t236 - t545 * t235 - (t229 / 0.2e1 - t291 / 0.2e1 + t426 / 0.2e1 + t357 + t526 + t551) * t269 + (t546 + t552) * t121 - t555 * t120 + (t299 / 0.2e1 - pkin(7) * t321 + t383 + (-0.2e1 * t480 + 0.3e1 / 0.2e1 * Ifges(3,4) * t347) * qJD(1)) * t347 * qJD(2) + (t494 + t499) * t312 / 0.2e1 + t52 * t63 + t58 * t64 + t57 * t65 + t55 * t66 + (-pkin(7) * t322 - t298 / 0.2e1 + t382 + (-0.2e1 * t481 - 0.3e1 / 0.2e1 * t436 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t347) * qJD(1) + (t257 + 0.2e1 * t451 + t513) * pkin(2)) * t397 + t15 * t110 + t112 * t29 + t123 * t111 + t7 * t152 + t10 * t153 + t9 * t154 + t8 * t155 + t185 * t124 + t186 * t125 + t53 * t212 + t54 * t213 + (-t280 * mrSges(4,3) + (Ifges(4,2) + t489) * t312 + t334 * mrSges(4,1) - Ifges(4,4) * t313) * t255 + t400 * t198 + t221 * t30 + t197 * t281; ((t254 * t346 - t255 * t343) * mrSges(4,3) + (t400 * t343 + (t212 * t345 - t213 * t342 + t281) * t346) * qJD(3)) * pkin(2) - m(5) * (t142 * t160 + t143 * t161 + t232 * t265) + t263 * t427 - m(4) * (-t262 * t265 + t263 * t266) + t348 + t357 * qJD(4) + (t363 + t349 + (-t212 * t342 - t213 * t345) * qJD(4)) * t331 + t495 * (-pkin(3) - t447) + 0.2e1 * ((t170 * t343 - t171 * t346) * t490 + (m(5) * (t232 * t343 + t346 * t496) / 0.2e1 + (-t262 * t343 + t263 * t346) * t490) * qJD(3)) * pkin(2) + t498 * t111 + t501 * t110 + t507 * t153 + t508 * t154 + (t103 * t320 + t180 * t498 - t249 * t5 + t250 * t6 + t36 * t508 + t37 * t507) * m(6) + t511 * t155 + t512 * t152 + (t1 * t250 + t12 * t231 + t249 * t3 + t31 * t511 + t32 * t512 + t501 * t56) * m(7) - t161 * t212 - t160 * t213 - t400 * t265 + t231 * t29 - t266 * t281 + t439 * t249 + t440 * t250 + t320 * t30 + ((-t335 / 0.2e1 - t299 / 0.2e1 + t383 + qJD(1) * t480 + (t321 - t417) * pkin(7)) * t347 + (t382 + t298 / 0.2e1 + (t481 + t436 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t347) * qJD(1) + (t322 + t416) * pkin(7) + (-t257 - t451) * pkin(2)) * t344) * qJD(1); t363 * pkin(9) + t506 * t154 + t509 * t155 + t510 * t152 + t504 * t110 + pkin(9) * t349 + t348 + t505 * t153 - m(5) * (t142 * t164 + t143 * t165 + t232 * t263) - t205 * t111 - t165 * t212 - t164 * t213 + t253 * t29 + (-t400 + t427) * t263 + ((-mrSges(5,3) * t142 - pkin(9) * t213) * t345 + (pkin(4) * t111 - pkin(9) * t212 - t431) * t342) * qJD(4) - t262 * t281 + t439 * t275 + t440 * t276 + t333 * t30 - t495 * pkin(3) + (t1 * t276 + t12 * t253 + t275 * t3 + t31 * t509 + t32 * t510 + t504 * t56) * m(7) + (t103 * t333 - t275 * t5 + t276 * t6 + t505 * t37 + t506 * t36 + (-t205 + t336) * t180) * m(6); ((t341 * t6 + t420 * t5) * pkin(4) - t180 * t446 + t36 * t41 - t37 * t42) * m(6) + t500 * t41 + (Ifges(5,5) * t277 - Ifges(5,6) * t278) * t460 + t167 * t461 + (Ifges(5,1) * t277 - t435) * t462 + (-t460 * t514 + t548) * t356 + t64 * t445 + t278 * t431 + t65 * t385 + (t437 - t212) * t142 + (t1 * t327 + t3 * t329 - t31 * t41 + t32 * t503 - t56 * t86) * m(7) + (-Ifges(5,2) * t278 + t168 + t274) * t463 + t494 + t493 + t503 * t152 - t86 * t110 - t42 * t153 + t143 * t213 - t232 * (mrSges(5,1) * t278 + mrSges(5,2) * t277) - t111 * t446 + t327 * t63 + t329 * t66 + (-t515 * t460 - t550 + t552) * t194; -(-t152 - t153) * t194 + t500 * t356 + t29 + t30 + (t194 * t32 - t31 * t356 + t12) * m(7) + (t194 * t37 + t356 * t36 + t103) * m(6); t356 * t110 - t294 * t152 + 0.2e1 * (t3 / 0.2e1 + t56 * t467 + t32 * t460) * m(7) + t66;];
tauc  = t2(:);
