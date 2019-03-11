% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:11
% EndTime: 2019-03-09 10:12:29
% DurationCPUTime: 9.84s
% Computational Cost: add. (22241->524), mult. (42701->672), div. (0->0), fcn. (50324->8), ass. (0->292)
t308 = sin(qJ(4));
t309 = sin(qJ(2));
t479 = -qJ(3) - pkin(7);
t290 = t479 * t309;
t306 = sin(pkin(10));
t311 = cos(qJ(2));
t440 = cos(pkin(10));
t373 = t440 * t311;
t254 = t306 * t290 - t479 * t373;
t281 = -t306 * t309 + t373;
t367 = -pkin(8) * t281 - t254;
t487 = cos(qJ(4));
t420 = t306 * t311;
t282 = -t440 * t309 - t420;
t515 = t440 * t290 + t479 * t420;
t545 = t282 * pkin(8) + t515;
t120 = t308 * t367 + t487 * t545;
t340 = t308 * t281 - t487 * t282;
t71 = -pkin(5) * t340 + t120;
t581 = qJ(5) * t71;
t246 = t487 * t281 + t282 * t308;
t566 = t308 * t545 - t487 * t367;
t570 = pkin(5) * t246 + t566;
t580 = t570 * t71;
t310 = cos(qJ(6));
t452 = t310 * mrSges(7,2);
t307 = sin(qJ(6));
t459 = t307 * mrSges(7,1);
t291 = t452 + t459;
t579 = t71 * t291;
t453 = t310 * mrSges(7,1);
t458 = t307 * mrSges(7,2);
t289 = -t453 + t458;
t148 = t289 * t246;
t300 = -pkin(2) * t311 - pkin(1);
t261 = -pkin(3) * t281 + t300;
t438 = qJ(5) * t340;
t331 = t261 - t438;
t500 = pkin(4) + pkin(9);
t66 = -t246 * t500 + t331;
t44 = -t307 * t66 - t310 * t71;
t525 = t289 * t340;
t417 = t307 * t340;
t543 = mrSges(7,1) * t246 - mrSges(7,3) * t417;
t578 = -t71 * t148 + t44 * t543 + t525 * t570;
t116 = -pkin(4) * t246 + t331;
t163 = mrSges(6,2) * t246 - mrSges(6,3) * t340;
t577 = m(6) * t116 + t163;
t304 = t307 ^ 2;
t305 = t310 ^ 2;
t408 = t304 + t305;
t575 = t408 * mrSges(7,3);
t574 = t307 * t570;
t573 = t310 * t570;
t572 = -pkin(4) * t566 + qJ(5) * t120;
t486 = m(6) * t566;
t569 = -t282 * mrSges(4,1) + t281 * mrSges(4,2);
t568 = -m(6) / 0.2e1;
t467 = Ifges(7,6) * t310;
t471 = Ifges(7,5) * t307;
t360 = t467 + t471;
t489 = t310 / 0.2e1;
t491 = t307 / 0.2e1;
t538 = -t246 / 0.2e1;
t563 = Ifges(5,4) + Ifges(6,6);
t517 = t563 * t538;
t535 = t340 / 0.2e1;
t470 = Ifges(7,2) * t310;
t476 = Ifges(7,4) * t307;
t361 = t470 + t476;
t541 = Ifges(7,6) * t246 + t340 * t361;
t475 = Ifges(7,4) * t310;
t478 = Ifges(7,1) * t307;
t362 = t475 + t478;
t542 = Ifges(7,5) * t246 + t340 * t362;
t567 = t489 * t541 + t491 * t542 - t261 * mrSges(5,2) + t116 * mrSges(6,3) - t360 * t538 + t517 - (Ifges(5,1) + Ifges(6,2) + Ifges(7,3)) * t535;
t376 = m(7) * t408;
t565 = t568 - t376 / 0.2e1;
t303 = t309 * pkin(2);
t262 = -pkin(3) * t282 + t303;
t547 = t246 * qJ(5);
t347 = t262 - t547;
t522 = t340 * t500;
t67 = t347 + t522;
t46 = -t307 * t67 + t573;
t564 = t46 / 0.2e1;
t559 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t558 = -t575 / 0.2e1;
t483 = pkin(4) * t340;
t557 = t547 - t483;
t554 = t340 * t565;
t553 = t246 / 0.2e1;
t390 = t440 * pkin(2);
t365 = t390 + pkin(3);
t484 = pkin(2) * t306;
t269 = t308 * t484 - t487 * t365;
t552 = t269 / 0.2e1;
t530 = mrSges(6,2) - mrSges(5,1);
t550 = mrSges(7,2) * t246;
t548 = qJ(5) * t525;
t546 = t246 * t291;
t47 = t310 * t67 + t574;
t503 = t47 / 0.2e1;
t504 = -Ifges(7,3) / 0.2e1;
t346 = t471 / 0.2e1 + t467 / 0.2e1;
t521 = t346 * t340;
t544 = mrSges(7,1) * t564 - mrSges(7,2) * t503 - t246 * t504 + t521;
t294 = -Ifges(7,2) * t307 + t475;
t295 = Ifges(7,1) * t310 - t476;
t341 = t294 * t489 + t295 * t491;
t427 = t246 * t310;
t160 = -mrSges(7,2) * t340 - mrSges(7,3) * t427;
t414 = t310 * t160;
t428 = t246 * t307;
t157 = mrSges(7,1) * t340 + mrSges(7,3) * t428;
t419 = t307 * t157;
t344 = t419 / 0.2e1 - t414 / 0.2e1;
t363 = mrSges(7,3) * (t304 / 0.2e1 + t305 / 0.2e1);
t540 = t246 * t363 - t344;
t270 = t308 * t365 + t487 * t484;
t372 = t408 * t270;
t531 = mrSges(5,2) - mrSges(6,3);
t539 = -mrSges(7,3) * t372 + t270 * t530 + (-t291 + t531) * t269;
t537 = -t295 / 0.4e1;
t536 = -t340 / 0.2e1;
t528 = -Ifges(5,6) + Ifges(6,5);
t240 = t340 * mrSges(6,2);
t519 = t558 * t246;
t516 = t340 * t360 / 0.4e1 + t570 * t289 / 0.2e1;
t514 = 0.2e1 * m(7);
t513 = 2 * qJD(4);
t512 = m(5) / 0.2e1;
t510 = m(6) / 0.2e1;
t509 = -m(7) / 0.2e1;
t508 = m(7) / 0.2e1;
t507 = m(4) * pkin(2);
t506 = -mrSges(6,1) / 0.2e1;
t95 = -t547 + t522;
t55 = t310 * t95 + t574;
t502 = -t55 / 0.2e1;
t501 = t570 / 0.4e1;
t268 = -pkin(4) + t269;
t263 = -pkin(9) + t268;
t495 = t263 / 0.2e1;
t267 = qJ(5) + t270;
t494 = t267 / 0.2e1;
t493 = t270 / 0.2e1;
t492 = -t307 / 0.2e1;
t490 = -t310 / 0.2e1;
t488 = t500 / 0.2e1;
t485 = m(6) * t557;
t54 = -t307 * t95 + t573;
t481 = t54 * mrSges(7,3);
t480 = t55 * mrSges(7,3);
t477 = Ifges(6,4) * t246;
t474 = Ifges(5,5) * t246;
t473 = Ifges(6,5) * t340;
t472 = Ifges(7,5) * t340;
t469 = Ifges(5,6) * t340;
t468 = Ifges(7,6) * t340;
t117 = t347 + t483;
t413 = t310 * t340;
t159 = mrSges(7,3) * t413 - t550;
t241 = t340 * mrSges(5,1);
t401 = -Ifges(5,4) / 0.2e1 - Ifges(6,6) / 0.2e1;
t85 = -t246 * t361 + t468;
t447 = t310 * t85;
t218 = Ifges(7,4) * t427;
t88 = -Ifges(7,1) * t428 - t218 + t472;
t454 = t307 * t88;
t320 = t447 / 0.2e1 + t454 / 0.2e1 + (t346 + t401) * t340 + (Ifges(5,2) + Ifges(6,3)) * t538 + t563 * t536;
t45 = -t307 * t71 + t310 * t66;
t1 = m(7) * (t44 * t46 + t45 * t47 + t580) + t46 * t157 + t45 * t159 + t47 * t160 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t309) * t309 + (-mrSges(4,1) * t303 + Ifges(4,4) * t281 + (-Ifges(4,1) + Ifges(4,2)) * t282) * t281 + (-mrSges(4,2) * t303 - Ifges(4,4) * t282) * t282 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t311 + (Ifges(3,1) - Ifges(3,2)) * t309) * t311 + (t262 * mrSges(5,2) + t320) * t340 - t116 * t240 + (-t262 * mrSges(5,1) - t246 * t401 + t559 * t340 - t567) * t246 + (m(4) * t303 + t569) * t300 + (m(5) * t262 + t241) * t261 + t577 * t117 + t578;
t466 = t1 * qJD(1);
t465 = t566 * mrSges(5,1);
t464 = t566 * mrSges(6,2);
t463 = t120 * mrSges(5,2);
t462 = t120 * mrSges(6,3);
t461 = t246 * mrSges(6,1);
t457 = t307 * t47;
t456 = t307 * t55;
t455 = t307 * t541;
t451 = t310 * mrSges(7,3);
t450 = t310 * t45;
t449 = t310 * t46;
t448 = t310 * t54;
t446 = t310 * t542;
t161 = t340 * t451 - t550;
t4 = t54 * t157 + t55 * t160 + t45 * t161 - t557 * t163 - t116 * t485 + m(7) * (t44 * t54 + t45 * t55 + t580) - (t517 + t567) * t246 + (t261 * mrSges(5,1) - t116 * mrSges(6,2) + t559 * t246 + t320) * t340 + t578;
t445 = t4 * qJD(1);
t151 = Ifges(7,2) * t428 - t218;
t152 = t295 * t246;
t217 = Ifges(7,5) * t427;
t5 = -t570 * t546 - t217 * t535 + t44 * t160 - t45 * t157 - ((t88 / 0.2e1 + t151 / 0.2e1 - t44 * mrSges(7,3)) * t310 + (-t152 / 0.2e1 - t85 / 0.2e1 - t468 / 0.2e1 - t45 * mrSges(7,3)) * t307) * t246;
t444 = t5 * qJD(1);
t439 = qJ(5) * t546;
t436 = qJ(5) * t269;
t358 = t307 * t45 + t310 * t44;
t10 = -t160 * t417 - t157 * t413 + t246 * t148 - m(7) * (t246 * t570 + t340 * t358) - m(4) * (t254 * t281 + t282 * t515) + (-m(6) - m(5)) * (-t120 * t340 + t246 * t566) + (-t281 ^ 2 - t282 ^ 2) * mrSges(4,3) + (-t246 ^ 2 - t340 ^ 2) * (mrSges(5,3) + mrSges(6,1));
t435 = qJD(1) * t10;
t18 = (-t414 + m(7) * (t307 * t44 - t450) + t419 - t577) * t340;
t434 = qJD(1) * t18;
t407 = t507 / 0.2e1;
t422 = t267 * t246;
t313 = (t270 * t246 + t269 * t340) * t512 + (t268 * t340 + t422) * t510 + (t263 * t340 * t408 + t422) * t508 + (t281 * t306 + t440 * t282) * t407;
t317 = t262 * t512 + t117 * t510 + (-t307 * t46 + t310 * t47) * t508 + t543 * t492 + t159 * t489 + t309 * t407;
t11 = -t546 / 0.2e1 + t241 - t240 - t313 + t317 + t535 * t575 + t531 * t246 + t569;
t433 = t11 * qJD(1);
t371 = t408 * t500;
t392 = t546 / 0.2e1 + t558 * t340;
t316 = t240 / 0.2e1 - t241 / 0.2e1 + t557 * t510 + (-t340 * t371 + t547) * t508 + t392;
t322 = t485 / 0.2e1 + (-t307 * t54 + t310 * t55) * t509 + t543 * t491 + t161 * t490;
t403 = mrSges(6,2) / 0.2e1 - mrSges(5,1) / 0.2e1;
t14 = t403 * t340 - 0.2e1 * t531 * t553 + t316 + t322;
t432 = t14 * qJD(1);
t336 = (t452 / 0.2e1 + t459 / 0.2e1) * t340;
t21 = t336 - t540;
t431 = t21 * qJD(1);
t345 = t458 / 0.2e1 - t453 / 0.2e1;
t338 = t345 * t340;
t343 = t157 * t490 + t160 * t492;
t23 = -t338 - t343;
t430 = t23 * qJD(1);
t293 = Ifges(7,5) * t310 - Ifges(7,6) * t307;
t424 = t246 * t293;
t421 = t267 * t269;
t418 = t307 * t159;
t415 = t310 * t543;
t411 = t500 * t157;
t410 = t500 * t159;
t57 = 0.2e1 * t554;
t409 = t57 * qJD(1);
t406 = m(7) / 0.4e1 + m(6) / 0.4e1;
t405 = mrSges(7,3) * t457;
t404 = t570 * t508;
t402 = -mrSges(6,3) / 0.2e1 + mrSges(5,2) / 0.2e1;
t394 = -t461 / 0.2e1;
t393 = -t451 / 0.2e1;
t385 = t500 * t489;
t379 = -t268 / 0.2e1 + t552;
t378 = t493 - t267 / 0.2e1;
t366 = t376 / 0.2e1;
t359 = t267 * t71 - t269 * t570;
t357 = t449 + t457;
t356 = t448 + t456;
t318 = t572 * t568 + (-t357 * t500 + t581) * t509 - t548 / 0.2e1;
t321 = ((t268 - t269) * t566 - (-t267 + t270) * t120) * t510 + t525 * t494 + t148 * t552;
t334 = t542 / 0.4e1 + t543 * t495 + t157 * t493;
t335 = -t541 / 0.4e1 + t161 * t495 + t160 * t493;
t3 = (-t542 / 0.4e1 + t543 * t488 + (-t54 / 0.2e1 + t564) * mrSges(7,3) + t334) * t310 + (t541 / 0.4e1 + t410 / 0.2e1 + (t502 + t503) * mrSges(7,3) + t335) * t307 + (t356 * t263 + t358 * t270 + t359) * t508 + (pkin(4) * t553 + t438 / 0.2e1 + t378 * t340 - t379 * t246) * mrSges(6,1) + t318 + t321 + t528 * (t535 + t536);
t37 = -m(7) * (t263 * t372 - t421) - m(6) * (t268 * t270 - t421) - t539;
t355 = t3 * qJD(1) - t37 * qJD(2);
t249 = t267 * t289;
t319 = t361 * t492 + t362 * t489 + t341;
t123 = t249 + t319;
t323 = t45 * t393 - t307 * t151 / 0.4e1 - t310 * t152 / 0.4e1 + mrSges(7,3) * t450 / 0.2e1 - t454 / 0.4e1 - t447 / 0.4e1 - t516 + (t294 + t362) * t428 / 0.4e1 + (t361 / 0.4e1 + t537) * t427;
t314 = t263 * t540 - t494 * t546 + t323;
t7 = t314 - t544;
t354 = t7 * qJD(1) - t123 * qJD(2);
t328 = t394 - t418 / 0.2e1 - t415 / 0.2e1;
t329 = (t506 + t345) * t246;
t16 = -t329 + (t501 - t457 / 0.4e1 - t449 / 0.4e1) * t514 + t328;
t348 = mrSges(6,3) + t291;
t173 = 0.4e1 * t406 * t267 + t348;
t351 = qJD(1) * t16 + qJD(2) * t173;
t349 = t54 * mrSges(7,1) / 0.2e1 + mrSges(7,2) * t502;
t342 = t161 * t492 + t490 * t543;
t337 = t345 * t270;
t288 = qJ(5) * t289;
t170 = t288 + t319;
t315 = -t249 / 0.2e1 - t288 / 0.2e1 - t319;
t60 = -t337 - t315;
t8 = t439 / 0.2e1 + (t468 / 0.2e1 + t160 * t488 + t152 / 0.4e1 + t85 / 0.4e1) * t310 + (t472 / 0.2e1 - t411 / 0.2e1 + t88 / 0.4e1 + t151 / 0.4e1) * t307 - (t504 + (t537 + t470 / 0.4e1) * t310 - t500 * t363 + (t475 / 0.2e1 + t294 / 0.4e1 + t478 / 0.4e1) * t307) * t246 + t349 + t516;
t333 = t8 * qJD(1) + t60 * qJD(2) + t170 * qJD(4);
t20 = -t345 * t246 + (t501 - t456 / 0.4e1 - t448 / 0.4e1) * t514 + t342;
t273 = (m(6) + m(7)) * qJ(5) + t348;
t325 = 0.2e1 * t406 * (0.2e1 * qJ(5) + t270) + t348;
t74 = t270 * t565 + t325;
t332 = qJD(1) * t20 + qJD(2) * t74 + qJD(4) * t273;
t73 = m(6) * t493 + t270 * t366 + t325;
t61 = -t337 + t315;
t56 = m(6) * t535 + t340 * t366 + t554;
t24 = -t338 + t343;
t22 = t336 - t344 - t519;
t19 = t356 * t508 + t404 + t486 - (-mrSges(6,1) + t345) * t246 - t342;
t17 = t357 * t508 + t566 * t510 + t404 + t486 / 0.2e1 - t329 - t328;
t15 = mrSges(6,3) * t538 - t240 / 0.2e1 + mrSges(5,2) * t553 + mrSges(5,1) * t535 - t402 * t246 + t316 - t322;
t13 = t313 + t317 + t392;
t9 = t323 + t521 + Ifges(7,3) * t553 - t439 / 0.2e1 - t411 * t492 - t160 * t385 + t349 + t519 * t500;
t6 = t314 + t544;
t2 = t321 - t318 - (mrSges(6,1) * t379 - t293 / 0.4e1 + Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t246 - t543 * t385 - t477 / 0.2e1 + t473 / 0.2e1 + t474 / 0.2e1 - t469 / 0.2e1 + t446 / 0.4e1 - t405 / 0.2e1 + t438 * t506 + t359 * t508 - t463 / 0.2e1 + t464 / 0.2e1 - t465 / 0.2e1 + t462 / 0.2e1 + t579 - t410 * t491 - t402 * t120 + (Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1 + t378 * mrSges(6,1) + t341) * t340 - t455 / 0.4e1 + t46 * t393 + t403 * t566 + ((t263 * t55 + t270 * t45) * t508 - t480 / 0.2e1 + t335) * t307 + ((t263 * t54 + t270 * t44) * t508 - t481 / 0.2e1 + t334) * t310 + pkin(4) * t394 + t424 / 0.4e1;
t12 = [qJD(2) * t1 - qJD(3) * t10 + qJD(4) * t4 + qJD(5) * t18 + qJD(6) * t5, t13 * qJD(3) + t2 * qJD(4) + t17 * qJD(5) + t6 * qJD(6) + t466 + (-t477 + t473 + t474 - t469 + t446 / 0.2e1 - t405 - t463 + t464 - t465 + t462 + (-t270 * mrSges(5,3) + t341) * t340 + Ifges(3,5) * t311 - Ifges(3,6) * t309 + Ifges(4,6) * t282 + Ifges(4,5) * t281 + t579 - t515 * mrSges(4,2) + (m(6) * t120 + m(7) * t71 - mrSges(6,1) * t340 + t525) * t267 - t455 / 0.2e1 + t269 * t246 * mrSges(5,3) + m(5) * (t120 * t270 + t269 * t566) - mrSges(7,3) * t449 + (-t311 * mrSges(3,1) + t309 * mrSges(3,2)) * pkin(7) + (-t281 * t390 + t282 * t484) * mrSges(4,3) + (m(7) * t357 + t415 + t418) * t263 - t254 * mrSges(4,1) + (-t254 * t440 + t306 * t515) * t507 + (t461 + t486) * t268 + t424 / 0.2e1) * qJD(2), qJD(2) * t13 + qJD(4) * t15 + qJD(5) * t56 + qJD(6) * t24 - t435, t445 + t2 * qJD(2) + t15 * qJD(3) + t19 * qJD(5) + t9 * qJD(6) + ((-t356 * t500 + t581) * t508 + t572 * t510) * t513 + (t548 + t579 + (-t481 - t500 * t543 + t542 / 0.2e1) * t310 + (-t480 - t500 * t161 - t541 / 0.2e1) * t307 - (Ifges(6,4) - Ifges(5,5) + pkin(4) * mrSges(6,1) - t293 / 0.2e1) * t246 + (-qJ(5) * mrSges(6,1) + t341 + t528) * t340 + t530 * t566 - t531 * t120) * qJD(4), qJD(2) * t17 + qJD(3) * t56 + qJD(4) * t19 + qJD(6) * t22 + t434, t444 + t6 * qJD(2) + t24 * qJD(3) + t9 * qJD(4) + t22 * qJD(5) + (-mrSges(7,1) * t45 - mrSges(7,2) * t44 + Ifges(7,6) * t428 - t217) * qJD(6); -qJD(3) * t11 + qJD(4) * t3 + qJD(5) * t16 + qJD(6) * t7 - t466, -qJD(4) * t37 + qJD(5) * t173 - qJD(6) * t123, -t433, t539 * qJD(4) + t73 * qJD(5) + t61 * qJD(6) + ((-t270 * t371 - t436) * t508 + (-pkin(4) * t270 - t436) * t510) * t513 + t355, qJD(4) * t73 + t351, t61 * qJD(4) + (-t263 * t291 - t360) * qJD(6) + t354; qJD(2) * t11 - qJD(4) * t14 + qJD(5) * t57 - qJD(6) * t23 + t435, t433, 0, -t432, t409, qJD(6) * t289 - t430; -qJD(2) * t3 + qJD(3) * t14 + qJD(5) * t20 - qJD(6) * t8 - t445, qJD(5) * t74 - qJD(6) * t60 - t355, t432, qJD(5) * t273 - qJD(6) * t170, t332 ((mrSges(7,2) * t500 - Ifges(7,6)) * t310 + (mrSges(7,1) * t500 - Ifges(7,5)) * t307) * qJD(6) - t333; -qJD(2) * t16 - qJD(3) * t57 - qJD(4) * t20 - qJD(6) * t21 - t434, -qJD(4) * t74 - t351, -t409, -t332, 0, -qJD(6) * t291 - t431; -qJD(2) * t7 + qJD(3) * t23 + qJD(4) * t8 + qJD(5) * t21 - t444, qJD(4) * t60 - t354, t430, t333, t431, 0;];
Cq  = t12;
