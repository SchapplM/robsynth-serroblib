% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:48:04
% EndTime: 2019-03-09 03:48:18
% DurationCPUTime: 7.05s
% Computational Cost: add. (19297->497), mult. (37159->661), div. (0->0), fcn. (42896->8), ass. (0->289)
t316 = sin(pkin(10));
t317 = cos(pkin(10));
t489 = sin(qJ(3));
t491 = cos(qJ(3));
t285 = t316 * t489 - t317 * t491;
t287 = t316 * t491 + t317 * t489;
t480 = -t317 * pkin(2) - pkin(1);
t228 = t285 * pkin(3) - t287 * qJ(4) + t480;
t170 = -pkin(4) * t285 - t228;
t320 = sin(qJ(5));
t490 = cos(qJ(5));
t238 = t285 * t320 + t490 * t287;
t355 = -t285 * t490 + t320 * t287;
t575 = m(6) * t170 + mrSges(6,1) * t355 + mrSges(6,2) * t238;
t319 = sin(qJ(6));
t321 = cos(qJ(6));
t474 = Ifges(7,4) * t321;
t300 = -Ifges(7,2) * t319 + t474;
t89 = Ifges(7,6) * t238 - t300 * t355;
t574 = -t89 / 0.4e1;
t475 = Ifges(7,4) * t319;
t302 = Ifges(7,1) * t321 - t475;
t92 = Ifges(7,5) * t238 - t302 * t355;
t573 = -t92 / 0.4e1;
t479 = pkin(7) + qJ(2);
t382 = t479 * t317;
t383 = t479 * t316;
t246 = t382 * t491 - t383 * t489;
t196 = t285 * pkin(8) + t246;
t245 = t382 * t489 + t491 * t383;
t348 = t287 * pkin(8) - t245;
t524 = t196 * t320 + t348 * t490;
t428 = t320 * t524;
t105 = -t490 * t196 + t320 * t348;
t560 = t490 * t105;
t571 = -t560 + t428;
t553 = t524 * mrSges(6,2);
t456 = t321 * mrSges(7,1);
t478 = mrSges(7,2) * t319;
t371 = t456 - t478;
t563 = t105 * t371;
t565 = t105 * mrSges(6,1);
t570 = t553 + t563 + t565;
t569 = -t563 / 0.2e1 - t565 / 0.2e1;
t492 = t321 / 0.2e1;
t496 = -t319 / 0.2e1;
t310 = Ifges(7,6) * t319;
t473 = Ifges(7,5) * t321;
t530 = -t310 + t473;
t557 = -t238 / 0.2e1;
t540 = Ifges(6,4) * t557;
t541 = t355 / 0.2e1;
t556 = t238 / 0.2e1;
t568 = t92 * t492 + t89 * t496 + t530 * t556 + t540 + (Ifges(6,2) + Ifges(7,3)) * t541;
t511 = mrSges(7,3) / 0.2e1;
t567 = pkin(5) * t105;
t566 = Ifges(7,3) * t556;
t322 = -pkin(3) - pkin(4);
t293 = -t320 * qJ(4) + t322 * t490;
t291 = pkin(5) - t293;
t564 = t105 * t291;
t562 = t105 * t524;
t297 = Ifges(7,5) * t319 + Ifges(7,6) * t321;
t561 = t238 * t297;
t314 = t319 ^ 2;
t315 = t321 ^ 2;
t416 = t314 + t315;
t381 = t416 * t355;
t288 = t320 * t371;
t353 = t416 * t490;
t338 = -t320 * mrSges(6,1) - t490 * mrSges(6,2) + mrSges(7,3) * t353;
t559 = t338 - t288;
t484 = pkin(5) * t238;
t164 = pkin(9) * t355 + t484;
t455 = t321 * mrSges(7,2);
t460 = t319 * mrSges(7,1);
t296 = t455 + t460;
t151 = t296 * t238;
t534 = t296 * t355;
t558 = -t105 * t151 - t524 * t534;
t485 = pkin(5) * t534;
t555 = mrSges(7,1) * t238;
t554 = mrSges(7,2) * t238;
t294 = qJ(4) * t490 + t320 * t322;
t551 = t294 * t524;
t550 = t319 * t524;
t549 = t321 * t524;
t548 = t371 * t238;
t546 = m(5) * t228 + mrSges(5,1) * t285 - mrSges(5,3) * t287;
t545 = t548 / 0.2e1 + t381 * t511;
t504 = -t355 / 0.2e1;
t542 = -t355 / 0.4e1;
t537 = -Ifges(5,5) + Ifges(4,4);
t372 = mrSges(7,3) * (t315 / 0.2e1 + t314 / 0.2e1);
t239 = pkin(3) * t287 + t285 * qJ(4);
t519 = m(5) / 0.2e1;
t535 = t239 * t519;
t408 = -t473 / 0.2e1;
t360 = t408 + t310 / 0.2e1;
t531 = t360 * t355;
t459 = t319 * mrSges(7,3);
t406 = -t459 / 0.2e1;
t74 = pkin(5) * t355 - pkin(9) * t238 + t170;
t49 = -t105 * t321 + t319 * t74;
t458 = t319 * t49;
t529 = t49 * t406 + t458 * t511;
t66 = t164 * t321 + t550;
t67 = t164 * t319 - t549;
t368 = -t66 * t319 + t67 * t321;
t403 = t319 * t490;
t374 = -t403 / 0.2e1;
t402 = t321 * t490;
t342 = mrSges(7,1) * t374 - mrSges(7,2) * t402 / 0.2e1;
t396 = t490 * t296;
t335 = t396 / 0.2e1 + t342;
t528 = qJD(6) * t335;
t527 = t335 * qJD(4);
t336 = -t396 / 0.2e1 + t342;
t526 = t336 * qJD(6);
t467 = t355 * Ifges(7,5);
t94 = t238 * t302 + t467;
t453 = t321 * t94;
t466 = t355 * Ifges(7,6);
t91 = t238 * t300 + t466;
t457 = t319 * t91;
t525 = t170 * mrSges(6,2) + Ifges(6,1) * t556 - Ifges(6,4) * t355 - t457 / 0.2e1 + t453 / 0.2e1;
t454 = t321 * mrSges(7,3);
t160 = mrSges(7,1) * t355 - t238 * t454;
t494 = t320 / 0.2e1;
t412 = t238 * t459;
t157 = -mrSges(7,2) * t355 - t412;
t508 = t157 / 0.2e1;
t337 = t151 * t494 + t160 * t374 + t402 * t508;
t48 = t105 * t319 + t321 * t74;
t344 = t402 * t49 - t403 * t48;
t415 = -t490 / 0.2e1;
t427 = t320 * t238;
t515 = m(7) / 0.2e1;
t517 = m(6) / 0.2e1;
t523 = t571 * t517 + (t344 + t428) * t515 + t337 + (t355 * t415 + t427 / 0.2e1) * mrSges(6,3);
t301 = Ifges(7,1) * t319 + t474;
t497 = t301 / 0.4e1;
t299 = Ifges(7,2) * t321 + t475;
t498 = -t299 / 0.4e1;
t522 = -t319 * (t497 + t300 / 0.4e1) + t321 * (t302 / 0.4e1 + t498);
t152 = t238 * t299;
t153 = t238 * t301;
t521 = -t321 * t152 / 0.4e1 - t319 * t153 / 0.4e1 + t524 * t296 / 0.2e1 + t453 / 0.4e1 - t457 / 0.4e1;
t520 = 0.2e1 * m(7);
t518 = -m(6) / 0.2e1;
t516 = -m(7) / 0.2e1;
t514 = -pkin(5) / 0.2e1;
t513 = -mrSges(7,1) / 0.2e1;
t512 = mrSges(7,2) / 0.2e1;
t510 = t66 / 0.2e1;
t509 = -t67 / 0.2e1;
t421 = t321 * t355;
t159 = -mrSges(7,3) * t421 - t555;
t507 = t159 / 0.2e1;
t502 = t355 / 0.4e1;
t500 = t291 / 0.2e1;
t499 = -t293 / 0.2e1;
t495 = t319 / 0.2e1;
t493 = -t321 / 0.2e1;
t488 = m(5) * t287;
t487 = m(7) * (-t490 + t353) * t320;
t486 = m(7) * (-pkin(5) * t320 + pkin(9) * t353);
t483 = pkin(5) * t296;
t482 = pkin(9) * t319;
t431 = t319 * t355;
t156 = -mrSges(7,3) * t431 + t554;
t189 = -pkin(4) * t287 - t239;
t231 = t238 * mrSges(6,1);
t265 = t285 * mrSges(4,2);
t266 = t287 * mrSges(5,1);
t411 = Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1;
t84 = t189 - t164;
t52 = t321 * t84 - t550;
t53 = t319 * t84 + t549;
t1 = t228 * t266 + t49 * t156 + t53 * t157 + t48 * t159 + t52 * t160 - t170 * t231 - t480 * t265 + m(7) * (t48 * t52 + t49 * t53 + t562) + (t480 * mrSges(4,1) - t537 * t287) * t287 + (t228 * mrSges(5,3) + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,3)) * t287 + t537 * t285) * t285 + (t530 * t541 + t525) * t355 + t546 * t239 + (Ifges(6,1) * t541 + Ifges(6,4) * t556 - t411 * t355 - t568) * t238 - t558 + t575 * t189;
t472 = t1 * qJD(1);
t465 = t355 * mrSges(6,2);
t292 = -pkin(9) + t294;
t464 = t292 * t66;
t463 = t292 * t67;
t462 = t293 * t48;
t461 = t293 * t49;
t407 = t466 / 0.2e1;
t5 = t49 * t160 - t524 * t548 + ((t94 / 0.2e1 - t152 / 0.2e1 + t467 / 0.2e1) * t319 + (t153 / 0.2e1 + t91 / 0.2e1 + t407 + t49 * mrSges(7,3)) * t321) * t238 + (-t157 - t412) * t48;
t452 = t5 * qJD(1);
t451 = t53 * t321;
t423 = t321 * t157;
t432 = t319 * t160;
t444 = t524 * t238;
t14 = -(t238 * mrSges(6,3) + t151) * t238 + (-mrSges(6,3) * t355 + t423 - t432) * t355 + m(7) * (-t444 + (-t319 * t48 + t321 * t49) * t355) + m(6) * (-t105 * t355 - t444) + (t285 ^ 2 + t287 ^ 2) * (mrSges(5,2) + mrSges(4,3)) + (m(3) * qJ(2) + mrSges(3,3)) * (t316 ^ 2 + t317 ^ 2) + (m(5) + m(4)) * (t245 * t287 - t246 * t285);
t448 = qJD(1) * t14;
t422 = t321 * t160;
t433 = t319 * t157;
t19 = (t433 + t422 + m(7) * (t321 * t48 + t458) - t546 + t575) * t287;
t447 = qJD(1) * t19;
t326 = -t535 + (t238 * t293 + t294 * t355) * t517 + (-t238 * t291 + t292 * t381) * t515 - t545;
t331 = t535 + t189 * t518 + (-t319 * t53 - t321 * t52) * t515 + t156 * t496 + t159 * t493;
t15 = -t287 * mrSges(4,1) - t285 * mrSges(5,3) - t231 + t265 - t266 + t326 - t331 + t465;
t441 = t15 * qJD(1);
t155 = t355 * t459 - t554;
t158 = t355 * t454 + t555;
t332 = (-t319 * t67 - t321 * t66) * t516 + mrSges(6,1) * t556 + t155 * t495 + t158 * t492;
t334 = t231 / 0.2e1 + (pkin(9) * t381 + t484) * t515 + t545;
t20 = 0.2e1 * mrSges(6,2) * t504 + t332 + t334;
t439 = t20 * qJD(1);
t437 = t355 * t530;
t359 = -t455 / 0.2e1 - t460 / 0.2e1;
t351 = t359 * t355;
t357 = t432 / 0.2e1 - t423 / 0.2e1;
t25 = t351 - t357;
t436 = t25 * qJD(1);
t435 = t291 * t296;
t430 = t319 * t299;
t429 = t319 * t320;
t425 = t321 * t155;
t424 = t321 * t156;
t420 = t321 * t301;
t397 = t490 * t238;
t329 = (t320 * t355 + t397) * t517 + (t320 * t381 + t397) * t515;
t384 = m(7) * t416;
t40 = -t488 - t329 + (-t384 / 0.2e1 + t518) * t287;
t419 = t40 * qJD(1);
t414 = t490 / 0.2e1;
t413 = qJD(4) * t487;
t404 = t294 * t490;
t399 = t490 * t548;
t393 = -t429 / 0.2e1;
t391 = -t422 / 0.2e1;
t390 = t502 + t542;
t389 = t301 / 0.2e1 + t300 / 0.2e1;
t387 = t302 / 0.2e1 - t299 / 0.2e1;
t380 = t416 * t293;
t369 = -t52 * t319 + t451;
t328 = (t551 - t564) * t515 - t534 * t500 + t294 * t151 / 0.2e1 + t569;
t349 = t573 - t292 * t158 / 0.2e1 + t160 * t499;
t350 = t574 + t292 * t155 / 0.2e1 + t293 * t508;
t3 = t485 / 0.2e1 + (t556 + t557) * Ifges(6,6) + (t541 + t504) * Ifges(6,5) + (pkin(5) * t515 + mrSges(6,1) / 0.2e1 + t371 / 0.2e1) * t105 + (-pkin(9) * t156 / 0.2e1 + t89 / 0.4e1 - t390 * t301 + (t509 - t53 / 0.2e1) * mrSges(7,3) + (-pkin(9) * t53 / 0.4e1 + t463 / 0.4e1 + t461 / 0.4e1) * t520 + t350) * t321 + (pkin(9) * t507 + t92 / 0.4e1 + t390 * t299 + (t510 + t52 / 0.2e1) * mrSges(7,3) + (pkin(9) * t52 / 0.4e1 - t464 / 0.4e1 - t462 / 0.4e1) * t520 + t349) * t319 + t328;
t339 = mrSges(7,3) * t380 - t293 * mrSges(6,2) + (-mrSges(6,1) - t371) * t294;
t55 = -m(7) * (t291 * t294 + t292 * t380) + t339;
t367 = t3 * qJD(1) - t55 * qJD(3);
t13 = (t560 + (t524 + t368) * t320 + t344) * t515 + t425 * t494 + t158 * t393 - t534 * t415 + t337;
t4 = t67 * t157 + t49 * t155 + t66 * t160 + t48 * t158 + m(7) * (t48 * t66 + t49 * t67 - t562) + (t170 * mrSges(6,1) + t540 + t568) * t238 + (t531 + (-Ifges(6,1) / 0.2e1 + t411) * t238 - t525) * t355 + t558;
t366 = t4 * qJD(1) + t13 * qJD(4);
t323 = t571 * t518 + (t320 * t369 - t560) * t516 + t534 * t414 + t429 * t507 - t320 * t424 / 0.2e1 + (-t238 * t494 + t355 * t414) * mrSges(6,3);
t11 = t323 + t523;
t345 = t353 * t292;
t72 = -m(7) * (t291 * t320 + t345) - m(6) * (-t293 * t320 + t404) - mrSges(5,3) - m(5) * qJ(4) + t559;
t365 = -t11 * qJD(1) + t72 * qJD(3);
t340 = t387 * t319 + t389 * t321;
t154 = t340 - t435;
t358 = -t433 / 0.2e1 + t391;
t324 = t358 * t292 + (-t292 * t372 - t522) * t238 - t437 / 0.4e1 + t548 * t500 - t521 + t529;
t346 = t512 * t53 + t513 * t52 + t566;
t7 = t324 + t346 + t531;
t364 = t7 * qJD(1) + t154 * qJD(3);
t352 = (-t478 / 0.2e1 + t456 / 0.2e1) * t287;
t23 = t399 / 0.2e1 + t352 + (t238 * t372 - t358) * t320;
t363 = -t23 * qJD(1) - qJD(3) * t336;
t362 = t512 * t67 + t513 * t66;
t361 = t486 / 0.2e1 - t288 / 0.2e1;
t356 = t430 / 0.2e1 - t420 / 0.2e1;
t354 = -Ifges(6,5) + t356;
t327 = t288 / 0.2e1 + (-t404 + t345 + (t291 + t380) * t320) * t515;
t54 = -t327 + t338 + t361;
t347 = t13 * qJD(1) - t54 * qJD(3) + t413;
t165 = t340 - t483;
t62 = (t514 - t291 / 0.2e1) * t296 + (mrSges(7,2) * t499 + t389) * t321 + (mrSges(7,1) * t499 + t387) * t319;
t330 = t358 * pkin(9) + t514 * t548 + t521 + t529;
t333 = -pkin(9) * t372 + t522;
t9 = (t530 / 0.4e1 - t360) * t355 + (-Ifges(7,3) / 0.2e1 + t333) * t238 + t330 + t362;
t341 = t9 * qJD(1) - t62 * qJD(3) + t165 * qJD(5) - t527;
t63 = t483 / 0.2e1 + t302 * t496 + t300 * t493 + t435 / 0.2e1 + t359 * t293 + t356;
t59 = t327 + t361;
t39 = t488 / 0.2e1 + 0.2e1 * (-t384 / 0.4e1 - m(6) / 0.4e1 - m(5) / 0.4e1) * t287 + t329;
t26 = t351 + t357;
t24 = t157 * t393 + t320 * t391 - t399 / 0.2e1 + t352 - t416 * mrSges(7,3) * t427 / 0.2e1;
t21 = mrSges(6,2) * t541 - t465 / 0.2e1 - t332 + t334;
t16 = t326 + t331;
t12 = t13 * qJD(5);
t10 = m(5) * t246 - t285 * mrSges(5,2) - t323 + t523;
t8 = t333 * t238 + t330 + t355 * t408 + t319 * t407 + t566 + t437 / 0.4e1 - t362;
t6 = t324 + Ifges(7,5) * t421 / 0.2e1 - Ifges(7,6) * t431 / 0.2e1 - t346;
t2 = (t369 * pkin(9) - t567) * t515 + t328 + t451 * t511 + t420 * t502 + pkin(9) * t424 / 0.2e1 + t430 * t542 - t159 * t482 / 0.2e1 + t52 * t406 - t561 / 0.2e1 - t553 - t485 / 0.2e1 + (mrSges(7,3) * t509 + (t461 + t463) * t515 + t355 * t497 + t350 + t574) * t321 + ((-t462 - t464) * t515 + t355 * t498 + mrSges(7,3) * t510 + t349 + t573) * t319 + 0.2e1 * t541 * Ifges(6,5) + 0.2e1 * t556 * Ifges(6,6) + t569;
t17 = [qJD(2) * t14 + qJD(3) * t1 + qJD(4) * t19 + qJD(5) * t4 - qJD(6) * t5, qJD(3) * t16 + qJD(4) * t39 + qJD(5) * t21 + qJD(6) * t26 + t448, t16 * qJD(2) + t10 * qJD(4) + t2 * qJD(5) + t6 * qJD(6) + t472 + (-t246 * mrSges(4,1) - t246 * mrSges(5,1) + t245 * mrSges(4,2) - t245 * mrSges(5,3) + t291 * t534 + (t292 * t156 - t53 * mrSges(7,3) + t89 / 0.2e1) * t321 + (-t292 * t159 + t52 * mrSges(7,3) + t92 / 0.2e1) * t319 + (-mrSges(5,2) * qJ(4) - Ifges(4,6) + Ifges(5,6)) * t287 + (mrSges(5,2) * pkin(3) - Ifges(5,4) - Ifges(4,5)) * t285 - (Ifges(6,6) - t294 * mrSges(6,3) - t297 / 0.2e1) * t238 + 0.2e1 * (t292 * t369 + t564) * t515 + 0.2e1 * (-pkin(3) * t246 - qJ(4) * t245) * t519 + 0.2e1 * (-t105 * t293 + t551) * t517 + (-t293 * mrSges(6,3) + t354) * t355 + t570) * qJD(3), qJD(2) * t39 + qJD(3) * t10 + qJD(6) * t24 + t12 + t447, t21 * qJD(2) + t2 * qJD(3) + t8 * qJD(6) + t366 + (t297 * t556 + m(7) * (pkin(9) * t368 + t567) + pkin(9) * t425 - t158 * t482 + t92 * t495 + t89 * t492 + t485 - Ifges(6,6) * t238 + t354 * t355 + t368 * mrSges(7,3) + t570) * qJD(5), -t452 + t26 * qJD(2) + t6 * qJD(3) + t24 * qJD(4) + t8 * qJD(5) + (-mrSges(7,1) * t49 - mrSges(7,2) * t48 - t561) * qJD(6); -qJD(3) * t15 + qJD(4) * t40 - qJD(5) * t20 - qJD(6) * t25 - t448, 0, -t441, t419, -t439, qJD(6) * t296 - t436; qJD(2) * t15 + qJD(4) * t11 + qJD(5) * t3 + qJD(6) * t7 - t472, t441, -qJD(4) * t72 - qJD(5) * t55 + qJD(6) * t154, t59 * qJD(5) - t365 + t413 + t528, t59 * qJD(4) + (m(7) * (-pkin(5) * t294 + pkin(9) * t380) + t339) * qJD(5) + t63 * qJD(6) + t367, t527 + t63 * qJD(5) + (-t292 * t371 - t530) * qJD(6) + t364; -qJD(2) * t40 - qJD(3) * t11 - qJD(6) * t23 + t12 - t447, -t419, -qJD(5) * t54 + t365 - t526, qJD(5) * t487 (t486 + t559) * qJD(5) + t526 + t347, qJD(5) * t336 - qJD(6) * t288 + t363; qJD(2) * t20 - qJD(3) * t3 + qJD(6) * t9 - t366, t439, qJD(4) * t54 - qJD(6) * t62 - t367, -t347 - t528, t165 * qJD(6) (-pkin(9) * t371 + t530) * qJD(6) + t341; qJD(2) * t25 - qJD(3) * t7 + qJD(4) * t23 - qJD(5) * t9 + t452, t436, qJD(4) * t336 + qJD(5) * t62 - t364, qJD(5) * t335 - t363, -t341, 0;];
Cq  = t17;
