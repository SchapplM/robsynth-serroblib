% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:44
% EndTime: 2019-12-31 22:00:07
% DurationCPUTime: 10.79s
% Computational Cost: add. (13184->610), mult. (28577->805), div. (0->0), fcn. (27948->6), ass. (0->321)
t331 = sin(qJ(4));
t333 = sin(qJ(2));
t334 = cos(qJ(4));
t335 = cos(qJ(3));
t450 = t334 * t335;
t332 = sin(qJ(3));
t453 = t332 * t333;
t270 = t331 * t453 - t333 * t450;
t513 = Ifges(5,4) + Ifges(6,4);
t628 = t513 * t270;
t298 = -t331 * t335 - t332 * t334;
t271 = t298 * t333;
t636 = t513 * t271;
t297 = -t331 * t332 + t450;
t620 = t513 * t297;
t621 = t513 * t298;
t572 = -pkin(8) - pkin(7);
t311 = t572 * t332;
t312 = t572 * t335;
t220 = t311 * t331 - t312 * t334;
t164 = qJ(5) * t297 + t220;
t287 = Ifges(6,6) * t298;
t288 = Ifges(5,6) * t298;
t289 = Ifges(6,5) * t297;
t290 = Ifges(5,5) * t297;
t411 = t290 + t288 + t289 + t287;
t598 = t311 * t334 + t312 * t331;
t605 = qJ(5) * t298 + t598;
t635 = -t220 * mrSges(5,1) - t164 * mrSges(6,1) - mrSges(5,2) * t598 - t605 * mrSges(6,2) + t411;
t632 = Ifges(5,1) + Ifges(6,1);
t630 = Ifges(5,2) + Ifges(6,2);
t633 = -t605 / 0.2e1;
t631 = Ifges(5,5) + Ifges(6,5);
t600 = Ifges(5,6) + Ifges(6,6);
t508 = mrSges(6,3) * t297;
t507 = mrSges(6,3) * t298;
t477 = t297 * mrSges(5,3);
t624 = t271 * t630 - t628;
t623 = -t270 * t632 + t636;
t596 = -t298 * t632 + t620;
t619 = t298 * t630 + t596 + t620;
t336 = cos(qJ(2));
t480 = t271 * mrSges(6,3);
t221 = mrSges(6,2) * t336 + t480;
t481 = t271 * mrSges(5,3);
t222 = mrSges(5,2) * t336 + t481;
t618 = t221 * t633 - t222 * t598 / 0.2e1;
t616 = t605 / 0.2e1;
t528 = pkin(4) * t298;
t615 = m(6) * t528;
t530 = pkin(3) * t334;
t253 = t270 * qJ(5);
t524 = t333 * pkin(7);
t308 = -pkin(2) * t336 - pkin(1) - t524;
t296 = t335 * t308;
t451 = t333 * t335;
t408 = -pkin(8) * t451 + t296;
t186 = (-pkin(6) * t332 - pkin(3)) * t336 + t408;
t449 = t335 * t336;
t249 = pkin(6) * t449 + t308 * t332;
t209 = -pkin(8) * t453 + t249;
t191 = t331 * t209;
t74 = t186 * t334 - t191;
t53 = t253 + t74;
t45 = -pkin(4) * t336 + t53;
t512 = -t45 + t53;
t585 = m(6) / 0.2e1;
t613 = t512 * t585;
t272 = t298 * t336;
t273 = t297 * t336;
t611 = t272 * t630 + t273 * t513 + t333 * t600;
t610 = t272 * t513 + t273 * t632 + t333 * t631;
t609 = -t336 * t600 + t624;
t597 = t297 * t630 - t621;
t607 = t271 * t632 + t628;
t606 = t297 * t632 + t621;
t476 = t298 * mrSges(5,3);
t604 = -t164 * t507 - t220 * t476;
t603 = t270 * t630 - t336 * t631 + t623 + t636;
t532 = pkin(3) * t331;
t525 = t333 * pkin(2);
t313 = -pkin(7) * t336 + t525;
t252 = -pkin(6) * t451 + t313 * t332;
t484 = t252 * mrSges(4,2);
t251 = pkin(6) * t453 + t313 * t335;
t485 = t251 * mrSges(4,1);
t595 = -t484 / 0.2e1 + t485 / 0.2e1;
t255 = Ifges(6,6) * t270;
t256 = Ifges(5,6) * t270;
t257 = Ifges(6,5) * t271;
t258 = Ifges(5,5) * t271;
t412 = t258 + t256 + t257 + t255;
t594 = -t251 * t332 + t252 * t335;
t472 = t335 * mrSges(4,1);
t474 = t332 * mrSges(4,2);
t404 = t472 - t474;
t452 = t332 * t336;
t442 = pkin(6) * t452;
t208 = t408 - t442;
t91 = t208 * t334 - t191;
t432 = t91 / 0.2e1 - t74 / 0.2e1;
t324 = -pkin(3) * t335 - pkin(2);
t250 = -pkin(4) * t297 + t324;
t401 = -mrSges(5,1) * t270 + mrSges(5,2) * t271;
t254 = t271 * mrSges(6,2);
t419 = -mrSges(6,1) * t270 + t254;
t593 = -t324 * t401 / 0.2e1 - t250 * t419 / 0.2e1;
t592 = t411 * t336 / 0.4e1;
t506 = Ifges(4,4) * t332;
t591 = (Ifges(4,1) - Ifges(4,2)) * t335 - t506;
t589 = t335 ^ 2;
t588 = -m(5) / 0.2e1;
t587 = m(5) / 0.2e1;
t586 = -m(6) / 0.2e1;
t584 = -mrSges(5,2) / 0.2e1;
t583 = -mrSges(6,2) / 0.2e1;
t582 = mrSges(5,3) / 0.2e1;
t581 = mrSges(6,3) / 0.2e1;
t580 = Ifges(4,1) / 0.2e1;
t579 = t45 / 0.2e1;
t190 = pkin(3) * t333 - pkin(8) * t449 + t251;
t211 = -pkin(8) * t452 + t252;
t82 = t190 * t334 - t211 * t331;
t49 = pkin(4) * t333 - qJ(5) * t273 + t82;
t578 = t49 / 0.2e1;
t577 = -t53 / 0.2e1;
t468 = qJ(5) * t271;
t193 = t334 * t209;
t75 = t186 * t331 + t193;
t54 = t75 + t468;
t576 = -t54 / 0.2e1;
t574 = -t75 / 0.2e1;
t571 = m(6) * t54;
t570 = pkin(2) * mrSges(4,1);
t569 = pkin(2) * mrSges(4,2);
t90 = -t208 * t331 - t193;
t64 = t90 - t468;
t568 = pkin(4) * t64;
t567 = pkin(6) * mrSges(4,1);
t566 = -t164 / 0.2e1;
t565 = t164 / 0.2e1;
t165 = -mrSges(6,1) * t271 - mrSges(6,2) * t270;
t564 = t165 / 0.2e1;
t329 = t333 * pkin(6);
t306 = pkin(3) * t453 + t329;
t188 = -pkin(4) * t271 + t306;
t563 = t188 / 0.2e1;
t562 = t220 / 0.2e1;
t482 = t270 * mrSges(6,3);
t225 = -mrSges(6,1) * t336 + t482;
t559 = -t225 / 0.2e1;
t483 = t270 * mrSges(5,3);
t226 = -mrSges(5,1) * t336 + t483;
t558 = -t226 / 0.2e1;
t227 = mrSges(6,1) * t333 - mrSges(6,3) * t273;
t557 = t227 / 0.2e1;
t556 = -t270 / 0.2e1;
t552 = t271 / 0.2e1;
t550 = t272 / 0.2e1;
t549 = t273 / 0.2e1;
t548 = t297 / 0.2e1;
t546 = -t298 / 0.2e1;
t544 = t298 / 0.2e1;
t542 = t306 / 0.2e1;
t323 = pkin(4) + t530;
t541 = -t323 / 0.2e1;
t540 = -t332 / 0.2e1;
t539 = -t332 / 0.4e1;
t537 = t335 / 0.2e1;
t536 = t335 / 0.4e1;
t534 = m(6) * t323;
t533 = m(6) * t331;
t531 = pkin(3) * t332;
t529 = pkin(4) * t270;
t527 = pkin(7) * t332;
t526 = pkin(7) * t335;
t330 = t336 * pkin(6);
t523 = t53 * mrSges(6,2);
t522 = t54 * mrSges(6,1);
t521 = t64 * mrSges(6,1);
t65 = t253 + t91;
t520 = t65 * mrSges(6,2);
t519 = t74 * mrSges(5,2);
t518 = t75 * mrSges(5,1);
t517 = t90 * mrSges(5,1);
t516 = t91 * mrSges(5,2);
t515 = mrSges(5,2) + mrSges(6,2);
t511 = mrSges(5,1) * t297;
t510 = mrSges(4,2) * t335;
t509 = mrSges(5,2) * t298;
t505 = Ifges(4,4) * t335;
t496 = Ifges(4,5) * t335;
t495 = Ifges(4,6) * t332;
t494 = pkin(3) * qJD(3);
t479 = t272 * mrSges(6,1);
t478 = t273 * mrSges(6,2);
t83 = t190 * t331 + t211 * t334;
t59 = qJ(5) * t272 + t83;
t475 = t331 * t59;
t473 = t332 * Ifges(4,2);
t166 = t478 - t479;
t167 = -mrSges(5,1) * t272 + mrSges(5,2) * t273;
t307 = pkin(3) * t452 + t330;
t189 = -pkin(4) * t272 + t307;
t223 = -mrSges(6,2) * t333 + mrSges(6,3) * t272;
t224 = -mrSges(5,2) * t333 + mrSges(5,3) * t272;
t228 = mrSges(5,1) * t333 - mrSges(5,3) * t273;
t248 = t296 - t442;
t388 = -t473 + t505;
t267 = Ifges(4,6) * t333 + t336 * t388;
t398 = Ifges(4,1) * t335 - t506;
t268 = Ifges(4,5) * t333 + t336 * t398;
t403 = mrSges(4,1) * t332 + t510;
t284 = t403 * t336;
t304 = -mrSges(4,2) * t333 - mrSges(4,3) * t452;
t305 = mrSges(4,1) * t333 - mrSges(4,3) * t449;
t366 = t496 / 0.2e1 - t495 / 0.2e1;
t378 = -t495 + t496;
t402 = -mrSges(5,1) * t271 - mrSges(5,2) * t270;
t437 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t439 = -Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t5 = m(6) * (t188 * t189 + t45 * t49 + t54 * t59) + m(5) * (t306 * t307 + t74 * t82 + t75 * t83) + (-pkin(1) * mrSges(3,1) + pkin(6) * t284 + t267 * t540 + t268 * t537 + (-Ifges(3,4) + t366) * t333 + t437 * t271 + t439 * t270 + (-t251 * t335 - t252 * t332) * mrSges(4,3)) * t333 + m(4) * (t248 * t251 + t249 * t252) + t249 * t304 + t248 * t305 + t306 * t167 + t74 * t228 + t59 * t221 + t83 * t222 + t54 * t223 + t75 * t224 + t49 * t225 + t82 * t226 + t45 * t227 + t188 * t166 + t189 * t165 + t307 * t402 + (-pkin(1) * mrSges(3,2) - t485 + t484 - t631 * t273 - t600 * t272 + (t589 * t580 + Ifges(3,1) - Ifges(3,2) - Ifges(4,3) - Ifges(5,3) - Ifges(6,3) + (m(4) * pkin(6) + t510) * pkin(6) + (t567 - t505 + t473 / 0.2e1) * t332) * t333 + (Ifges(3,4) - t378) * t336) * t336 + t610 * t556 + t611 * t552 + t624 * t550 + t623 * t549;
t471 = t5 * qJD(1);
t440 = mrSges(4,3) * t453;
t367 = mrSges(4,2) * t336 - t440;
t368 = -mrSges(4,1) * t336 - mrSges(4,3) * t451;
t443 = pkin(3) * t451;
t370 = t443 - t529;
t377 = Ifges(4,5) * t332 + Ifges(4,6) * t335;
t406 = -t74 * mrSges(5,3) - t45 * mrSges(6,3);
t415 = -t188 * t419 - t306 * t401 - t482 * t54 - t483 * t75;
t420 = -t630 + t632;
t6 = -m(5) * (t306 * t443 + t74 * t90 + t75 * t91) - t402 * t443 + (mrSges(4,3) * t249 + Ifges(4,4) * t451) * t451 + t415 - t65 * t221 - t91 * t222 - t64 * t225 - t90 * t226 + t513 * t270 ^ 2 + (t270 * t420 - t336 * t439 - t406 - t636) * t271 + (t258 / 0.2e1 + t256 / 0.2e1 + t257 / 0.2e1 + t255 / 0.2e1 - t333 * t377 + t437 * t270) * t336 - t370 * t165 - m(6) * (t188 * t370 + t45 * t64 + t54 * t65) + t249 * t368 + (-pkin(6) * t404 + t332 * t591) * t333 ^ 2 + (-t367 - t440) * t248;
t470 = t6 * qJD(1);
t7 = -t415 + t53 * t221 + t74 * t222 - t54 * t225 - t75 * t226 + t406 * t271 - t165 * t529 + m(6) * (-t188 * t529 + t512 * t54) + t607 * t556 + t609 * t270 / 0.2e1 - t412 * t336 / 0.2e1 + t603 * t552;
t469 = t7 * qJD(1);
t23 = m(6) * (t270 * t45 + t271 * t54) + t271 * t221 + t270 * t225;
t466 = qJD(1) * t23;
t465 = t164 * t270;
t463 = t164 * t323;
t462 = t220 * t270;
t461 = t250 * t270;
t202 = -mrSges(6,1) * t297 - mrSges(6,2) * t298;
t458 = t270 * t202;
t457 = t270 * t331;
t456 = t298 * t331;
t455 = t306 * t332;
t441 = pkin(4) * t461;
t438 = -Ifges(4,2) / 0.2e1 + t580;
t436 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t435 = pkin(4) / 0.2e1 + t541;
t434 = -t64 / 0.2e1 + t576;
t433 = -t90 / 0.2e1 + t574;
t431 = m(6) * pkin(4) + mrSges(6,1);
t430 = -t508 / 0.2e1;
t429 = t508 / 0.2e1;
t428 = t507 / 0.2e1;
t425 = -t481 / 0.2e1;
t424 = -t480 / 0.2e1;
t423 = t476 / 0.2e1;
t422 = t458 / 0.2e1;
t286 = t298 * mrSges(6,1);
t418 = mrSges(6,2) * t297 - t286;
t400 = -mrSges(5,1) * t298 + mrSges(5,2) * t297;
t413 = -t250 * t418 - t324 * t400 + t604;
t410 = m(6) * t578 + t557;
t409 = Ifges(4,3) / 0.2e1 + t436;
t407 = t436 * t333;
t405 = -t528 + t531;
t399 = -t509 - t511;
t397 = Ifges(4,1) * t332 + t505;
t387 = Ifges(4,2) * t335 + t506;
t343 = -t439 * t273 + t437 * t272 + mrSges(6,1) * t578 + t59 * t583 + t82 * mrSges(5,1) / 0.2e1 + t83 * t584;
t337 = (-t465 / 0.2e1 + t605 * t552) * mrSges(6,3) + (-t462 / 0.2e1 + t598 * t552) * mrSges(5,3) + t286 * t563 + t343 + t593;
t369 = -Ifges(5,2) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1;
t341 = mrSges(5,1) * t542 + (Ifges(5,6) / 0.4e1 + Ifges(6,6) / 0.4e1) * t336 + t628 + (m(6) * t563 + t564) * pkin(4) + t369 * t271;
t342 = -t636 + (Ifges(5,5) / 0.4e1 + Ifges(6,5) / 0.4e1) * t336 + t369 * t270 + t188 * t583 + t306 * t584;
t345 = -t164 * t559 - t220 * t558 + t618;
t346 = (t75 + t90) * t598 + (-t74 + t91) * t220;
t347 = (t54 + t64) * t605 + (-t45 + t65) * t164;
t349 = t334 * t228 / 0.2e1 + (t224 / 0.2e1 + t223 / 0.2e1) * t331;
t357 = (t331 * t83 + t334 * t82) * t587;
t363 = t332 * t402;
t1 = (pkin(3) * t455 + t346) * t588 + (t188 * t531 + t347 - t441) * t586 + t337 + (pkin(3) * t475 + t323 * t49) * t585 + (t165 * t540 - t363 / 0.2e1 + t349) * pkin(3) + ((-t65 / 0.2e1 + t579) * mrSges(6,3) - t432 * mrSges(5,3) + t342) * t297 + t323 * t557 + pkin(4) * t422 + ((-t569 / 0.2e1 - t567 / 0.2e1 + t438 * t332) * t332 + (0.2e1 * t506 + t570 / 0.2e1 - pkin(6) * mrSges(4,2) / 0.2e1 - t438 * t335 + (-t202 / 0.2e1 + t511 / 0.2e1 + t509 / 0.2e1 + t324 * t588 + t250 * t586) * pkin(3)) * t335 + t409) * t333 + (mrSges(5,3) * t433 + mrSges(6,3) * t434 + t341) * t298 + pkin(3) * t357 + (t290 / 0.4e1 + t288 / 0.4e1 + t289 / 0.4e1 + t287 / 0.4e1 + (-t472 / 0.2e1 + t474 / 0.2e1) * pkin(7) + t378) * t336 + t345 + t595;
t365 = pkin(3) * t399;
t10 = (mrSges(5,3) * t220 + mrSges(6,3) * t164 + t621) * t298 - m(5) * t324 * t531 + (-t505 + t569) * t335 + (t298 * t420 - t620) * t297 + t413 + (-t365 + t570 - t591) * t332 + (-m(6) * t250 - t202) * t405;
t375 = -qJD(1) * t1 - qJD(2) * t10;
t11 = -t202 * t528 - t250 * t615 + t597 * t544 + t546 * t606 + t548 * t619 - t413 + t604;
t344 = t225 * t565 + t226 * t562 + t618;
t3 = ((t577 + t579) * mrSges(6,3) + t342) * t297 + t337 + (t557 + t422) * pkin(4) + t407 + (pkin(4) * t578 + t441 / 0.2e1 + t605 * t576 + t54 * t616 + t45 * t565 + t164 * t577) * m(6) + t341 * t298 + t344 + t592;
t374 = -t3 * qJD(1) + t11 * qJD(2);
t339 = (t270 * t544 + t271 * t548) * mrSges(6,3) + (t164 * t271 + t270 * t605 + t297 * t54 + t298 * t45) * t585 + t221 * t548 + t225 * t544;
t348 = t189 * t585 - t479 / 0.2e1 + t478 / 0.2e1;
t17 = t339 - t348;
t28 = m(6) * (t164 * t297 + t298 * t605) + (t297 ^ 2 + t298 ^ 2) * mrSges(6,3);
t373 = qJD(1) * t17 + qJD(2) * t28;
t350 = (t270 * t323 + t271 * t532) * t585;
t360 = m(6) * t370;
t31 = t350 - t360 / 0.2e1 - t419;
t358 = m(6) * (t297 * t532 + t298 * t323);
t359 = t405 * t585;
t42 = t359 - t358 / 0.2e1 + t418;
t372 = qJD(1) * t31 - qJD(2) * t42;
t139 = t418 - t615;
t97 = t270 * t431 - t254;
t371 = qJD(1) * t97 - qJD(2) * t139;
t351 = t164 * t530;
t13 = (t616 + t633) * mrSges(6,2) + (t565 + t566) * mrSges(6,1) + (-t220 / 0.2e1 + t562) * mrSges(5,1) + (t463 / 0.2e1 - t351 / 0.2e1 + pkin(4) * t566) * m(6) + (-t530 / 0.2e1 - t435) * t508;
t153 = t515 * t530 + (mrSges(5,1) + mrSges(6,1) + (t323 - t530) * m(6)) * t532;
t340 = ((t221 / 0.2e1 + t222 / 0.2e1 + t571 / 0.2e1 + t425) * t334 + (t558 + t559 + t613 + (t581 + t582) * t270) * t331) * pkin(3);
t9 = t435 * t480 + (t65 / 0.2e1 + t577) * mrSges(6,2) + t432 * mrSges(5,2) + t434 * mrSges(6,1) + t433 * mrSges(5,1) + (-t568 / 0.2e1 + t54 * t541) * m(6) + t340;
t352 = t9 * qJD(1) - t13 * qJD(2) - t153 * qJD(3);
t338 = t605 * t424 + t598 * t425 + t400 * t542 + t418 * t563 + t75 * t423 + t54 * t428 + t45 * t430 + t462 * t582 + t465 * t581 + t343 - t592 - t593 + t619 * t271 / 0.4e1 + t603 * t297 / 0.4e1 + (-t607 / 0.4e1 + t609 / 0.4e1) * t298 + (-t606 / 0.4e1 + t597 / 0.4e1) * t270;
t85 = t358 / 0.2e1 + t359;
t60 = t350 + t360 / 0.2e1;
t16 = t339 + t348;
t12 = (t351 - t463) * t585 + t323 * t430 + t429 * t530 + (m(6) * t566 + t430) * pkin(4) + t635;
t8 = t568 * t585 + pkin(4) * t424 - t520 / 0.2e1 + t521 / 0.2e1 - t516 / 0.2e1 + t517 / 0.2e1 - t523 / 0.2e1 - t522 / 0.2e1 - t519 / 0.2e1 - t518 / 0.2e1 + (t424 - t571 / 0.2e1) * t323 + t340 + t412;
t4 = t338 + t407 + t53 * t429 - t165 * t528 / 0.2e1 + t476 * t574 + t507 * t576 + t164 * t613 - t344 + (t410 - t458 / 0.2e1 + (-t188 * t298 - t461) * t585) * pkin(4);
t2 = (-t387 / 0.2e1 + t365 / 0.2e1) * t451 + (-Ifges(4,5) * t536 - Ifges(4,6) * t539 - t378 / 0.4e1 + t366) * t336 + t346 * t587 + ((t324 * t451 + t455) * t587 + t475 * t585 + t349 + t357 + t363 / 0.2e1) * pkin(3) + t595 + t64 * t428 + t65 * t429 + t90 * t423 + (t188 * t405 + t250 * t370 + t347) * t585 + t405 * t564 - (t332 ^ 2 + t589) * mrSges(4,3) * t524 / 0.2e1 + t338 + t403 * t329 / 0.2e1 - t404 * t525 / 0.2e1 - t368 * t526 / 0.2e1 - t367 * t527 / 0.2e1 - t397 * t453 / 0.2e1 + t409 * t333 + t410 * t323 + t370 * t202 / 0.2e1 - t345 + t432 * t477 + 0.2e1 * (t388 * t539 + t398 * t536) * t333;
t14 = [qJD(2) * t5 - qJD(3) * t6 + qJD(4) * t7 + qJD(5) * t23, t2 * qJD(3) + t4 * qJD(4) + t16 * qJD(5) + t471 + (0.2e1 * (t220 * t83 + t307 * t324 + t598 * t82) * t587 + t598 * t228 + (t600 * t297 - t298 * t631 + t377) * t333 / 0.2e1 + mrSges(3,2) * t329 + t304 * t526 + t267 * t537 + t59 * t508 + t49 * t507 + t82 * t476 + t83 * t477 + (Ifges(3,5) + t397 * t537 + t387 * t540 + (-mrSges(3,1) - t404) * pkin(6)) * t336 - t305 * t527 + t610 * t546 + t611 * t548 + 0.2e1 * (t164 * t59 + t189 * t250 + t49 * t605) * t585 + t605 * t227 - Ifges(3,6) * t333 + t332 * t268 / 0.2e1 + t324 * t167 - pkin(2) * t284 + t250 * t166 + t164 * t223 + t220 * t224 + t189 * t202 + t307 * t399 + t594 * mrSges(4,3) + m(4) * (-pkin(2) * t330 + pkin(7) * t594) + t596 * t549 + t597 * t550) * qJD(2), -t470 + t2 * qJD(2) + (-t249 * mrSges(4,1) - t248 * mrSges(4,2) - Ifges(4,5) * t453 - Ifges(4,6) * t451 - t323 * t480 + t534 * t64 + t412 - t516 + t517 - t520 + t521) * qJD(3) + t8 * qJD(4) + t60 * qJD(5) + (mrSges(6,3) * t457 + (-t271 * t334 + t457) * mrSges(5,3) + t65 * t533 + m(5) * (t331 * t91 + t334 * t90)) * t494, t469 + t4 * qJD(2) + t8 * qJD(3) + (-t518 - t522 - t519 - t523 + (-t480 - t571) * pkin(4) + t412) * qJD(4), qJD(2) * t16 + qJD(3) * t60 + t466; -qJD(3) * t1 - qJD(4) * t3 + qJD(5) * t17 - t471, -qJD(3) * t10 + qJD(4) * t11 + qJD(5) * t28, (-pkin(7) * t404 - t164 * t534 - t323 * t508 + t378 + t635) * qJD(3) + t12 * qJD(4) + t85 * qJD(5) + (mrSges(6,3) * t456 + (-t297 * t334 + t456) * mrSges(5,3) + t605 * t533 + m(5) * (-t220 * t334 + t331 * t598)) * t494 + t375, t12 * qJD(3) + ((-m(6) * t164 - t508) * pkin(4) + t635) * qJD(4) + t374, qJD(3) * t85 + t373; qJD(2) * t1 + qJD(4) * t9 + qJD(5) * t31 + t470, -qJD(4) * t13 - qJD(5) * t42 - t375, -t153 * qJD(4), (-t515 * t334 + (-mrSges(5,1) - t431) * t331) * qJD(4) * pkin(3) + t352, t372; qJD(2) * t3 - qJD(3) * t9 + qJD(5) * t97 - t469, qJD(3) * t13 - qJD(5) * t139 - t374, -t352, 0, t371; -qJD(2) * t17 - qJD(3) * t31 - qJD(4) * t97 - t466, qJD(3) * t42 + qJD(4) * t139 - t373, -t372, -t371, 0;];
Cq = t14;
