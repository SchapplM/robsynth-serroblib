% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:25:03
% EndTime: 2019-03-09 15:25:21
% DurationCPUTime: 10.17s
% Computational Cost: add. (22741->511), mult. (43808->673), div. (0->0), fcn. (51146->8), ass. (0->306)
t309 = sin(pkin(10));
t310 = cos(pkin(10));
t312 = sin(qJ(3));
t313 = sin(qJ(2));
t315 = cos(qJ(3));
t316 = cos(qJ(2));
t377 = t312 * t313 - t315 * t316;
t444 = t312 * t316 + t313 * t315;
t248 = t309 * t444 + t310 * t377;
t311 = sin(qJ(6));
t314 = cos(qJ(6));
t525 = Ifges(7,4) * t314;
t583 = Ifges(7,1) * t311;
t389 = t525 + t583;
t566 = -t309 * t377 + t310 * t444;
t341 = -Ifges(7,5) * t248 + t389 * t566;
t388 = -Ifges(7,2) * t311 + t525;
t365 = t314 * t388;
t526 = Ifges(7,4) * t311;
t390 = Ifges(7,1) * t314 - t526;
t366 = t311 * t390;
t538 = t314 / 0.2e1;
t544 = -t248 / 0.2e1;
t581 = Ifges(7,2) * t314;
t387 = t526 + t581;
t543 = -t566 / 0.2e1;
t605 = t248 / 0.2e1;
t587 = Ifges(7,6) * t605 + t387 * t543;
t597 = t566 / 0.2e1;
t330 = t311 * t587 + t341 * t538 + (Ifges(7,5) * t314 - Ifges(7,6) * t311) * t544 - Ifges(4,6) * t444 - Ifges(4,5) * t377 + (t366 + t365) * t597 - (Ifges(5,5) - Ifges(6,4)) * t248 + (Ifges(6,5) - Ifges(5,6)) * t566;
t552 = -pkin(8) - pkin(7);
t296 = t552 * t313;
t297 = t552 * t316;
t578 = t296 * t315 + t297 * t312;
t594 = t578 * mrSges(4,2);
t378 = t312 * t296 - t315 * t297;
t595 = t378 * mrSges(4,1);
t219 = -qJ(4) * t377 + t378;
t591 = -qJ(4) * t444 + t578;
t600 = t219 * t310 + t309 * t591;
t615 = t600 * mrSges(6,2);
t616 = t600 * mrSges(5,1);
t599 = -t219 * t309 + t310 * t591;
t617 = t599 * mrSges(6,3);
t618 = t599 * mrSges(5,2);
t305 = t311 * mrSges(7,1);
t306 = t314 * mrSges(7,2);
t577 = t306 + t305;
t612 = -pkin(5) * t566 + t599;
t623 = t612 * t577;
t634 = t330 - t594 - t595 + t615 + t617 - t616 - t618 + t623;
t419 = -t316 * pkin(2) - pkin(1);
t337 = pkin(3) * t377 + t419;
t329 = -qJ(5) * t566 + t337;
t533 = t248 * pkin(4);
t328 = t329 + t533;
t507 = t248 * mrSges(6,2);
t633 = -m(6) * t328 + mrSges(6,3) * t566 + t507;
t632 = -m(5) * t337 - mrSges(5,1) * t248 - mrSges(5,2) * t566;
t631 = -t594 / 0.2e1 - t595 / 0.2e1 + t623 / 0.2e1;
t536 = pkin(2) * t315;
t304 = pkin(3) + t536;
t459 = t310 * t312;
t441 = pkin(2) * t459;
t274 = t304 * t309 + t441;
t270 = qJ(5) + t274;
t628 = t270 * t612;
t300 = pkin(3) * t309 + qJ(5);
t627 = t300 * t612;
t611 = -pkin(5) * t248 + t600;
t626 = t311 * t611;
t625 = t314 * t611;
t624 = t611 * t612;
t499 = t314 * mrSges(7,1);
t500 = t311 * mrSges(7,2);
t393 = t499 - t500;
t172 = t393 * t566;
t174 = -mrSges(7,3) * t311 * t566 - mrSges(7,1) * t248;
t498 = t314 * mrSges(7,3);
t176 = mrSges(7,2) * t248 + t498 * t566;
t364 = t393 * t248;
t508 = t566 * mrSges(5,3);
t509 = t566 * mrSges(6,1);
t91 = (pkin(4) + pkin(9)) * t248 + t329;
t57 = -t311 * t91 - t314 * t612;
t58 = -t311 * t612 + t314 * t91;
t580 = t444 * mrSges(4,3);
t593 = t378 * t580;
t622 = t612 * t364 + t611 * t172 - t57 * t174 - t58 * t176 - t337 * (mrSges(5,1) * t566 - mrSges(5,2) * t248) + t593 - t419 * (mrSges(4,1) * t444 - mrSges(4,2) * t377) - t328 * (-mrSges(6,2) * t566 + mrSges(6,3) * t248) + (t508 + t509) * t600;
t620 = t615 / 0.2e1 - t616 / 0.2e1 + t617 / 0.2e1 - t618 / 0.2e1;
t619 = Ifges(6,2) - Ifges(5,2);
t585 = mrSges(6,1) + mrSges(5,3);
t614 = t585 * t600;
t560 = m(6) / 0.2e1;
t613 = t600 * t560;
t537 = pkin(2) * t312;
t298 = t309 * t537;
t273 = t304 * t310 - t298;
t271 = -pkin(4) - t273;
t609 = t270 * t599 + t271 * t600;
t608 = -t273 * t600 + t274 * t599;
t607 = t309 * t599 - t310 * t600;
t301 = -pkin(3) * t310 - pkin(4);
t606 = (t300 * t599 + t301 * t600) * t560;
t308 = t314 ^ 2;
t501 = t308 * mrSges(7,3);
t307 = t311 ^ 2;
t502 = t307 * mrSges(7,3);
t603 = -t501 / 0.2e1 - t502 / 0.2e1;
t529 = mrSges(6,2) - mrSges(5,1);
t596 = Ifges(6,6) * t566;
t521 = Ifges(7,6) * t314;
t523 = Ifges(7,5) * t311;
t386 = t521 + t523;
t579 = -Ifges(5,4) + t386;
t475 = t248 ^ 2;
t406 = t444 * pkin(3);
t342 = pkin(4) * t566 + qJ(5) * t248 + t406;
t532 = t313 * pkin(2);
t331 = t342 + t532;
t372 = t406 + t532;
t446 = -t311 * t174 / 0.2e1 + t176 * t538;
t557 = m(7) / 0.2e1;
t562 = m(5) / 0.2e1;
t245 = t566 * pkin(9);
t92 = t245 + t331;
t59 = -t311 * t92 + t625;
t60 = t314 * t92 + t626;
t592 = t372 * t562 + t331 * t560 + (-t311 * t59 + t314 * t60) * t557 + t446;
t280 = (t309 * t315 + t459) * pkin(2);
t281 = t310 * t536 - t298;
t530 = mrSges(5,2) - mrSges(6,3);
t590 = -(mrSges(4,1) * t312 + mrSges(4,2) * t315) * pkin(2) + t529 * t280 + (-t530 + t577) * t281;
t556 = m(7) / 0.4e1;
t559 = m(6) / 0.4e1;
t589 = 0.4e1 * t556 + 0.4e1 * t559;
t588 = m(6) + m(5);
t474 = t248 * t311;
t586 = -t474 / 0.2e1;
t584 = Ifges(5,1) + Ifges(7,3);
t582 = Ifges(4,4) * t444;
t539 = t311 / 0.2e1;
t166 = t176 * t539;
t369 = t174 * t538 + t166;
t576 = t603 * t248;
t175 = mrSges(7,1) * t566 - mrSges(7,3) * t474;
t452 = t314 * t175;
t473 = t248 * t314;
t177 = -mrSges(7,2) * t566 + mrSges(7,3) * t473;
t456 = t311 * t177;
t368 = -t452 / 0.2e1 - t456 / 0.2e1;
t553 = -mrSges(7,2) / 0.2e1;
t554 = mrSges(7,1) / 0.2e1;
t104 = t245 + t342;
t63 = -t104 * t311 + t625;
t64 = t104 * t314 + t626;
t574 = t553 * t64 + t554 * t63;
t573 = t553 * t60 + t554 * t59;
t494 = t60 * t311;
t495 = t59 * t314;
t384 = t494 + t495;
t567 = -t280 * t599 + t281 * t600;
t299 = -pkin(9) + t301;
t460 = t310 * t248;
t462 = t309 * t566;
t463 = t301 * t248;
t464 = t300 * t566;
t465 = t300 * t172;
t565 = t606 + (t299 * t384 + t627) * t557 - t465 / 0.2e1 + (t607 * t562 + (-t462 / 0.2e1 + t460 / 0.2e1) * mrSges(5,3)) * pkin(3) + (-t464 / 0.2e1 - t463 / 0.2e1) * mrSges(6,1) + (-t495 / 0.2e1 - t494 / 0.2e1) * mrSges(7,3) + t631;
t564 = 2 * qJD(3);
t563 = -m(5) / 0.2e1;
t561 = -m(6) / 0.2e1;
t558 = -m(7) / 0.2e1;
t555 = m(5) * pkin(3);
t551 = pkin(5) * mrSges(7,1);
t550 = t175 / 0.2e1;
t549 = -t177 / 0.2e1;
t541 = -t270 / 0.2e1;
t540 = -t300 / 0.2e1;
t524 = Ifges(7,5) * t566;
t522 = Ifges(7,6) * t566;
t520 = Ifges(7,3) * t248;
t359 = t377 * mrSges(4,3);
t361 = Ifges(4,2) * t377;
t363 = Ifges(4,1) * t377;
t476 = t566 ^ 2;
t1 = -m(4) * t419 * t532 - m(7) * (t57 * t59 + t58 * t60 + t624) + (t313 ^ 2 - t316 ^ 2) * Ifges(3,4) + t578 * t359 - t59 * t175 - t60 * t177 + pkin(1) * (mrSges(3,1) * t313 + mrSges(3,2) * t316) + (Ifges(3,2) - Ifges(3,1)) * t313 * t316 + (-mrSges(4,1) * t532 - mrSges(4,3) * t578 - Ifges(4,4) * t377) * t377 + (-mrSges(4,2) * t532 - mrSges(4,3) * t378 - t361 + t363 + t582) * t444 + (Ifges(6,6) - t579) * (t476 - t475) + (-t614 + (-0.2e1 * t526 - t581) * t473 - t474 * t583 + (t584 - Ifges(6,3) + t619) * t248) * t566 + t622 + t632 * t372 + t633 * t331;
t519 = t1 * qJD(1);
t332 = t314 * (t248 * t387 + t522);
t333 = t311 * (t248 * t389 + t524);
t2 = t473 * t587 + t341 * t586 - t593 + t596 * t597 - m(7) * (t57 * t63 + t58 * t64 + t624) - t377 ^ 2 * Ifges(4,4) - t63 * t175 - t64 * t177 - Ifges(6,6) * t475 / 0.2e1 + (-t361 / 0.2e1 + t363 / 0.2e1 + t582 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t377) * t444 + (t332 + t333 - t520 - t596) * t543 + (Ifges(5,4) * t597 - Ifges(6,2) * t544 + t579 * t543 - t614 + (-Ifges(5,2) + t584) * t605) * t566 + (Ifges(6,6) * t544 - Ifges(6,3) * t566 - Ifges(5,1) * t543 + t619 * t597 + (-Ifges(5,4) + t579) * t605) * t248 + t622 + t632 * t406 + t633 * t342;
t510 = t2 * qJD(1);
t506 = t248 * mrSges(6,1);
t497 = t314 * t58;
t236 = Ifges(7,5) * t473;
t371 = -t526 + (Ifges(7,1) - Ifges(7,2)) * t314;
t404 = t600 * mrSges(7,1);
t5 = t236 * t597 + t57 * t177 - t58 * t175 + ((t524 / 0.2e1 - t57 * mrSges(7,3) + Ifges(7,4) * t473 + t611 * mrSges(7,2)) * t314 + (-t522 + t404 - t58 * mrSges(7,3) + (t371 - t551) * t248) * t311) * t248;
t496 = t5 * qJD(1);
t493 = t63 * t314;
t492 = t64 * t311;
t385 = t311 * t58 + t314 * t57;
t10 = -(m(7) * t611 + t588 * t600 - t364) * t248 + t585 * (t475 + t476) + (m(7) * t385 - t588 * t599 + t452 + t456) * t566;
t491 = qJD(1) * t10;
t450 = t314 * t177;
t457 = t311 * t175;
t22 = (m(7) * (t311 * t57 - t497) - t450 + t457 - m(6) * (t337 + t533) + t507 + (m(6) * qJ(5) + mrSges(6,3)) * t566) * t566;
t490 = qJD(1) * t22;
t269 = -pkin(9) + t271;
t443 = t307 + t308;
t402 = t443 * t566;
t421 = t544 * t577 + t566 * t603;
t472 = t248 * t270;
t321 = (-t248 * t274 - t273 * t566) * t562 + (t271 * t566 - t472) * t560 + (t269 * t402 - t472) * t557 + t421;
t11 = t248 * t530 + t529 * t566 + t321 - t592;
t487 = t11 * qJD(1);
t400 = t443 * t299;
t471 = t248 * t300;
t322 = (t301 * t566 - t471) * t561 + (t400 * t566 - t471) * t558 - (-t248 * t309 - t310 * t566) * t555 / 0.2e1;
t323 = t342 * t560 + (-t311 * t63 + t314 * t64) * t557 - t406 * t563 + t446;
t395 = (t307 / 0.2e1 + t308 / 0.2e1) * mrSges(7,3);
t13 = -(-t577 / 0.2e1 + t530) * t248 + (t395 - t529) * t566 + t322 + t323;
t486 = t13 * qJD(1);
t354 = (t306 / 0.2e1 + t305 / 0.2e1) * t566;
t370 = t457 / 0.2e1 - t450 / 0.2e1;
t375 = t248 * t395;
t23 = t354 + t375 + t370;
t477 = t23 * qJD(1);
t373 = -t500 / 0.2e1 + t499 / 0.2e1;
t355 = t373 * t566;
t25 = t355 - t368;
t470 = t25 * qJD(1);
t469 = t270 * t172;
t468 = t270 * t281;
t466 = t281 * t300;
t455 = t311 * t269;
t454 = t311 * t299;
t449 = t314 * t269;
t448 = t314 * t299;
t356 = (-m(7) * t443 / 0.4e1 - m(6) / 0.4e1) * t566;
t374 = t443 * t556 + t559;
t66 = -0.2e1 * t374 * t566 + 0.2e1 * t356;
t447 = t66 * qJD(1);
t440 = mrSges(7,3) * t492;
t439 = t270 * t509;
t438 = t273 * t248 * mrSges(5,3);
t437 = t274 * t508;
t432 = Ifges(7,2) / 0.2e1 - Ifges(7,1) / 0.2e1;
t431 = mrSges(6,3) + t577;
t429 = t174 * t449;
t428 = -t520 / 0.2e1;
t427 = t506 / 0.2e1;
t424 = t498 / 0.2e1;
t411 = -t455 / 0.2e1;
t409 = t448 / 0.2e1;
t408 = 0.2e1 * t597;
t407 = t280 / 0.2e1 + t541;
t401 = t443 * t280;
t397 = t387 * t539 - t366 / 0.2e1 - t365 / 0.2e1 - t314 * t389 / 0.2e1;
t383 = t492 + t493;
t317 = t585 * (t280 * t543 + t281 * t605) - t429 / 0.2e1 - t620 + t469 / 0.2e1 + (t269 * t383 + t281 * t611 + t628) * t558 + t439 / 0.2e1 + (t385 * t558 + t368) * t280 + t440 / 0.2e1 + (t567 + t608) * t563 + (t567 + t609) * t561 + t437 / 0.2e1 + t271 * t427 + t176 * t411 + t63 * t424 - t438 / 0.2e1 + t281 * t364 / 0.2e1 - t631;
t4 = t299 * t166 + t174 * t409 + t317 + t565 + t620;
t40 = mrSges(7,3) * t401 - m(7) * (t269 * t401 + t468) - m(6) * (t271 * t280 + t468) - m(5) * (-t273 * t280 + t274 * t281) - t590;
t382 = -qJD(1) * t4 - qJD(2) * t40;
t367 = t270 * t393;
t158 = t367 + t397;
t334 = pkin(5) * t553 - t311 * t432 + 0.2e1 * t525;
t338 = t408 * Ifges(7,6) - t404 / 0.2e1;
t339 = t408 * Ifges(7,5) + t600 * mrSges(7,2) / 0.2e1;
t351 = t432 * t314 + t551 / 0.2e1;
t6 = t428 + t269 * t375 + (t269 * t549 + (mrSges(7,2) * t541 + t351) * t248 + t338) * t314 + (t269 * t550 + (mrSges(7,1) * t541 + t334) * t248 + t339) * t311 + t573;
t380 = -qJD(1) * t6 + qJD(2) * t158;
t327 = (-mrSges(6,1) / 0.2e1 - t373) * t248 + t613 + t611 * t557;
t319 = t427 + t327 - t369;
t335 = t384 * t557 + t613;
t18 = t319 - t335;
t183 = t270 * t589 + t431;
t376 = -qJD(1) * t18 - qJD(2) * t183;
t352 = t577 * t605;
t277 = t300 * t393;
t179 = Ifges(7,4) * t308 + t311 * t371 - t277;
t8 = t428 + t299 * t375 + (t299 * t549 + (mrSges(7,2) * t540 + t351) * t248 + t338) * t314 + (t299 * t550 + (mrSges(7,1) * t540 + t334) * t248 + t339) * t311 + t574;
t82 = -t277 / 0.2e1 + (mrSges(7,1) * t407 + t525) * t314 + (-mrSges(7,2) * t407 + t371) * t311;
t349 = qJD(1) * t8 + qJD(2) * t82 + qJD(3) * t179;
t348 = 0.2e1 * t374 * t280;
t336 = t383 * t557 + t613;
t20 = t319 - t336;
t251 = t300 * t589 + t431;
t324 = -t431 + (t561 + t558) * (t441 + 0.2e1 * qJ(5) + (pkin(3) + t304) * t309);
t89 = t348 + t324;
t347 = qJD(1) * t20 - qJD(2) * t89 + qJD(3) * t251;
t326 = -t333 / 0.4e1 - t332 / 0.4e1 + t388 * t586 - t389 * t474 / 0.4e1 - mrSges(7,3) * t497 / 0.2e1 + t58 * t424 + t611 * t393 / 0.2e1 + t428 + (t390 / 0.2e1 - t387 / 0.4e1) * t473 + (-t386 / 0.4e1 + t523 / 0.2e1 + t521 / 0.2e1) * t566;
t325 = -t506 / 0.2e1 + t327 + t369;
t90 = t348 - t324;
t83 = t367 / 0.2e1 + t277 / 0.2e1 + t373 * t280 + t397;
t65 = m(6) * t597 + t402 * t557 + 0.2e1 * t356;
t26 = t355 + t368;
t24 = t354 - t370 + t576;
t21 = t325 + t336;
t19 = t325 + t335;
t17 = -t322 + t323 + t421;
t14 = t321 + t592;
t9 = t300 * t352 + t177 * t409 - t175 * t454 / 0.2e1 + t326 + t576 * t299 + t574;
t7 = t270 * t352 + t177 * t449 / 0.2e1 + t175 * t411 + t326 + t576 * t269 + t573;
t3 = t369 * t299 + (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t599 + (mrSges(6,2) / 0.2e1 - mrSges(5,1) / 0.2e1) * t600 - t317 + t330 + t565;
t12 = [-qJD(2) * t1 - qJD(3) * t2 + qJD(4) * t10 + qJD(5) * t22 + qJD(6) * t5, -t519 + (-t580 * t537 + (-mrSges(3,1) * t316 + mrSges(3,2) * t313) * pkin(7) + t429 - t271 * t506 - t469 - t439 + t176 * t455 + t359 * t536 - t384 * mrSges(7,3) - t437 + t438 - Ifges(3,6) * t313 + Ifges(3,5) * t316 + t634) * qJD(2) + t3 * qJD(3) + t14 * qJD(4) + t19 * qJD(5) + t7 * qJD(6) + 0.2e1 * ((t269 * t384 + t628) * t557 + t609 * t560 + t608 * t562 + m(4) * (t312 * t578 - t315 * t378) * pkin(2) / 0.2e1) * qJD(2), -t510 + t3 * qJD(2) + t17 * qJD(4) + t21 * qJD(5) + t9 * qJD(6) + ((t299 * t383 + t627) * t557 + t606) * t564 + (-mrSges(7,3) * t493 + t174 * t448 + t176 * t454 - t440 - t465 + (m(5) * t607 + (t460 - t462) * mrSges(5,3)) * pkin(3) + (-t463 - t464) * mrSges(6,1) + t634) * qJD(3), qJD(2) * t14 + qJD(3) * t17 + qJD(5) * t65 + qJD(6) * t26 + t491, qJD(2) * t19 + qJD(3) * t21 + qJD(4) * t65 + qJD(6) * t24 + t490, t496 + t7 * qJD(2) + t9 * qJD(3) + t26 * qJD(4) + t24 * qJD(5) + (-mrSges(7,1) * t58 - mrSges(7,2) * t57 - Ifges(7,6) * t474 + t236) * qJD(6); -qJD(3) * t4 + qJD(4) * t11 + qJD(5) * t18 - qJD(6) * t6 + t519, -qJD(3) * t40 + qJD(5) * t183 + qJD(6) * t158 ((-t501 - t502) * t280 + t590) * qJD(3) + t90 * qJD(5) + t83 * qJD(6) + ((t280 * t400 + t466) * t557 + (t280 * t301 + t466) * t560 + (-t280 * t310 + t281 * t309) * t555 / 0.2e1) * t564 + t382, t487, qJD(3) * t90 - t376, t83 * qJD(3) + (-t269 * t577 - t386) * qJD(6) + t380; qJD(2) * t4 - qJD(4) * t13 + qJD(5) * t20 - qJD(6) * t8 + t510, -qJD(5) * t89 - qJD(6) * t82 - t382, qJD(5) * t251 - qJD(6) * t179, -t486, t347 (-t299 * t577 - t386) * qJD(6) - t349; -qJD(2) * t11 + qJD(3) * t13 + qJD(5) * t66 - qJD(6) * t25 - t491, -t487, t486, 0, t447, -qJD(6) * t393 - t470; -qJD(2) * t18 - qJD(3) * t20 - qJD(4) * t66 - qJD(6) * t23 - t490, qJD(3) * t89 + t376, -t347, -t447, 0, -qJD(6) * t577 - t477; qJD(2) * t6 + qJD(3) * t8 + qJD(4) * t25 + qJD(5) * t23 - t496, qJD(3) * t82 - t380, t349, t470, t477, 0;];
Cq  = t12;
