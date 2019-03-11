% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:18
% EndTime: 2019-03-08 20:45:34
% DurationCPUTime: 8.99s
% Computational Cost: add. (14446->642), mult. (33331->903), div. (0->0), fcn. (34340->10), ass. (0->364)
t390 = cos(qJ(4));
t586 = t390 / 0.2e1;
t384 = sin(qJ(6));
t388 = cos(qJ(6));
t385 = sin(qJ(5));
t389 = cos(qJ(5));
t500 = t388 * t389;
t431 = t384 * t385 - t500;
t432 = t384 * t389 + t388 * t385;
t433 = t384 * t432 - t388 * t431;
t619 = m(7) * pkin(5);
t489 = t619 / 0.2e1;
t545 = t389 * mrSges(6,1);
t549 = t385 * mrSges(6,2);
t551 = t432 * mrSges(7,2);
t552 = t431 * mrSges(7,1);
t635 = -t552 / 0.2e1 - t551 / 0.2e1;
t653 = t549 / 0.2e1 - t545 / 0.2e1 - t433 * t489 - t635;
t652 = -mrSges(7,1) / 0.2e1;
t651 = mrSges(7,2) / 0.2e1;
t383 = cos(pkin(6));
t386 = sin(qJ(4));
t382 = sin(pkin(6));
t391 = cos(qJ(2));
t515 = t382 * t391;
t312 = t383 * t390 - t386 * t515;
t387 = sin(qJ(2));
t503 = t387 * t389;
t232 = -t312 * t385 + t382 * t503;
t507 = t385 * t387;
t233 = t312 * t389 + t382 * t507;
t130 = t232 * t384 + t233 * t388;
t454 = t388 * t232 - t233 * t384;
t40 = -t130 * mrSges(7,1) - t454 * mrSges(7,2);
t650 = t40 * qJD(6);
t615 = -pkin(10) - pkin(9);
t354 = t615 * t385;
t356 = t615 * t389;
t247 = t354 * t384 - t356 * t388;
t453 = t388 * t354 + t356 * t384;
t320 = Ifges(7,6) * t432;
t321 = Ifges(7,5) * t431;
t491 = -t321 - t320;
t59 = -t247 * mrSges(7,1) - t453 * mrSges(7,2) + t491;
t649 = t59 * qJD(6);
t581 = pkin(9) * t390;
t584 = pkin(4) * t386;
t344 = qJ(3) - t581 + t584;
t325 = t389 * t344;
t392 = -pkin(2) - pkin(8);
t458 = -t385 * t392 + pkin(5);
t497 = t389 * t390;
t486 = pkin(10) * t497;
t209 = t386 * t458 + t325 - t486;
t504 = t386 * t392;
t262 = t385 * t344 + t389 * t504;
t506 = t385 * t390;
t225 = -pkin(10) * t506 + t262;
t536 = t225 * t384;
t110 = t209 * t388 - t536;
t261 = -t385 * t504 + t325;
t224 = t261 - t486;
t123 = t224 * t388 - t536;
t648 = -t110 + t123;
t299 = t432 * t386;
t508 = t385 * t386;
t301 = -t384 * t508 + t386 * t500;
t456 = -t301 * mrSges(7,1) + t299 * mrSges(7,2);
t647 = qJD(6) * t456;
t355 = pkin(4) * t390 + pkin(9) * t386;
t333 = t389 * t355;
t505 = t386 * t389;
t217 = pkin(10) * t505 + t390 * t458 + t333;
t496 = t390 * t392;
t269 = t385 * t355 + t389 * t496;
t240 = pkin(10) * t508 + t269;
t116 = t217 * t388 - t240 * t384;
t117 = t217 * t384 + t240 * t388;
t599 = -t301 / 0.2e1;
t601 = t299 / 0.2e1;
t428 = Ifges(7,5) * t599 + Ifges(7,6) * t601;
t631 = t116 * t652 + t117 * t651 - t428;
t627 = Ifges(7,3) * t586 - t631;
t646 = -t454 / 0.2e1;
t535 = t225 * t388;
t111 = t209 * t384 + t535;
t122 = -t224 * t384 - t535;
t643 = t111 + t122;
t378 = t385 ^ 2;
t380 = t389 ^ 2;
t490 = t378 + t380;
t641 = qJD(4) * (mrSges(6,3) * t490 - mrSges(5,2));
t270 = (-t386 * t507 + t389 * t391) * t382;
t271 = (t385 * t391 + t386 * t503) * t382;
t151 = t270 * t388 - t271 * t384;
t152 = t270 * t384 + t271 * t388;
t582 = pkin(5) * t389;
t371 = -pkin(4) - t582;
t527 = t271 * t389;
t528 = t270 * t385;
t435 = t527 - t528;
t516 = t382 * t387;
t481 = t390 * t516;
t622 = -m(7) / 0.2e1;
t624 = -m(6) / 0.2e1;
t640 = -(t527 / 0.2e1 - t528 / 0.2e1) * mrSges(6,3) + (pkin(4) * t481 + pkin(9) * t435) * t624 + (t151 * t453 + t152 * t247 - t371 * t481) * t622;
t550 = t432 * mrSges(7,3);
t637 = t386 * t390;
t300 = t432 * t390;
t249 = -mrSges(7,2) * t386 - mrSges(7,3) * t300;
t636 = t249 * t646;
t377 = Ifges(6,4) * t389;
t634 = -Ifges(6,2) * t385 + t377;
t350 = Ifges(6,1) * t385 + t377;
t268 = -t385 * t496 + t333;
t436 = -t268 * t385 + t269 * t389;
t447 = t545 - t549;
t298 = t384 * t506 - t388 * t497;
t193 = -mrSges(7,1) * t298 - mrSges(7,2) * t300;
t311 = t383 * t386 + t390 * t515;
t521 = t311 * t193;
t524 = t300 * t454;
t526 = t298 * t130;
t616 = mrSges(7,3) / 0.2e1;
t617 = -mrSges(7,3) / 0.2e1;
t632 = -t524 * t617 + t526 * t616 + t521 / 0.2e1;
t553 = t300 * mrSges(7,1);
t554 = t298 * mrSges(7,2);
t493 = -t553 / 0.2e1 + t554 / 0.2e1;
t178 = t432 * t311;
t179 = t431 * t311;
t426 = t178 * t652 + t179 * t651;
t576 = Ifges(7,4) * t298;
t181 = -Ifges(7,2) * t300 + t386 * Ifges(7,6) - t576;
t280 = Ifges(7,4) * t300;
t183 = -Ifges(7,1) * t298 + t386 * Ifges(7,5) - t280;
t196 = Ifges(7,2) * t298 - t280;
t197 = -Ifges(7,1) * t300 + t576;
t583 = pkin(5) * t385;
t457 = -t392 + t583;
t327 = t457 * t390;
t592 = t371 / 0.2e1;
t606 = t249 / 0.2e1;
t226 = mrSges(7,1) * t432 - mrSges(7,2) * t431;
t613 = t226 / 0.2e1;
t628 = t453 * t606 - (t183 / 0.4e1 + t196 / 0.4e1) * t431 - (-t197 / 0.4e1 + t181 / 0.4e1) * t432 + t327 * t613 + t193 * t592;
t626 = 0.2e1 * qJD(4);
t625 = m(5) / 0.2e1;
t623 = m(6) / 0.2e1;
t621 = m(7) / 0.2e1;
t618 = mrSges(6,1) / 0.2e1;
t614 = t454 / 0.2e1;
t227 = t551 + t552;
t612 = t227 / 0.2e1;
t611 = t232 / 0.2e1;
t610 = t233 / 0.2e1;
t609 = -t453 / 0.2e1;
t248 = -mrSges(7,2) * t390 + mrSges(7,3) * t299;
t608 = t248 / 0.2e1;
t251 = mrSges(7,1) * t386 + mrSges(7,3) * t298;
t605 = -t251 / 0.2e1;
t604 = t251 / 0.2e1;
t603 = -t298 / 0.2e1;
t602 = -t299 / 0.2e1;
t600 = -t300 / 0.2e1;
t598 = t311 / 0.2e1;
t315 = t447 * t390;
t597 = -t315 / 0.2e1;
t319 = t390 * t350;
t596 = -t319 / 0.4e1;
t595 = t432 / 0.2e1;
t485 = mrSges(6,3) * t506;
t335 = -mrSges(6,2) * t386 - t485;
t594 = -t335 / 0.2e1;
t579 = mrSges(6,2) * t389;
t580 = mrSges(6,1) * t385;
t347 = t579 + t580;
t593 = t347 / 0.2e1;
t591 = t385 / 0.2e1;
t590 = t386 / 0.4e1;
t589 = -t389 / 0.2e1;
t588 = t389 / 0.2e1;
t587 = -t390 / 0.2e1;
t577 = Ifges(6,4) * t385;
t575 = Ifges(7,4) * t432;
t376 = Ifges(6,5) * t389;
t573 = Ifges(6,6) * t385;
t571 = pkin(5) * qJD(5);
t570 = t110 * mrSges(7,2);
t569 = t111 * mrSges(7,1);
t566 = t122 * mrSges(7,1);
t565 = t123 * mrSges(7,2);
t548 = t386 * mrSges(5,2);
t547 = t386 * Ifges(6,5);
t546 = t386 * Ifges(6,6);
t544 = t390 * mrSges(5,1);
t543 = -t447 - mrSges(5,1);
t542 = mrSges(7,3) * qJD(4);
t541 = t110 * t300;
t540 = t116 * t388;
t539 = t117 * t384;
t538 = t152 * t431;
t523 = t301 * t298;
t525 = t299 * t300;
t409 = (t523 / 0.2e1 - t525 / 0.2e1) * mrSges(7,3) + t249 * t602 + t251 * t599;
t400 = t193 * t587 + t409;
t22 = t400 - t635;
t537 = t22 * qJD(2);
t534 = t232 * t385;
t533 = t233 * t389;
t532 = t453 * t300;
t531 = t247 * t298;
t522 = t301 * t432;
t520 = t311 * t385;
t519 = t311 * t386;
t518 = t327 * t385;
t452 = t311 * t481;
t35 = m(6) * (t232 * t270 + t233 * t271 - t452) + m(7) * (t130 * t152 + t151 * t454 - t452) + m(5) * (-t311 * t390 + t312 * t386 + t515) * t516;
t517 = t35 * qJD(1);
t514 = t384 * t248;
t513 = t384 * t298;
t295 = t390 * t634 + t546;
t512 = t385 * t295;
t511 = t385 * t335;
t336 = mrSges(6,1) * t390 + mrSges(6,3) * t505;
t510 = t385 * t336;
t348 = Ifges(6,2) * t389 + t577;
t509 = t385 * t348;
t250 = mrSges(7,1) * t390 + mrSges(7,3) * t301;
t502 = t388 * t250;
t501 = t388 * t300;
t334 = -mrSges(6,2) * t390 + mrSges(6,3) * t508;
t499 = t389 * t334;
t484 = mrSges(6,3) * t497;
t337 = mrSges(6,1) * t386 - t484;
t498 = t389 * t337;
t492 = -Ifges(7,5) * t300 + Ifges(7,6) * t298;
t438 = -t533 + t534;
t487 = -0.1e1 + t490;
t31 = ((-t312 - t438) * t390 - t487 * t519) * t623 + (-t178 * t299 + t179 * t301 - t312 * t390 + t519 - t524 - t526) * t621;
t488 = t31 * qJD(4);
t480 = t392 * t516;
t479 = t573 / 0.2e1;
t474 = -t550 / 0.2e1;
t473 = -t547 / 0.2e1;
t472 = t226 * t587 + t301 * t474 - t522 * t617;
t471 = t226 * t598;
t470 = -t516 / 0.2e1;
t195 = t553 - t554;
t468 = t195 * t591;
t467 = -t511 / 0.2e1;
t466 = -t110 / 0.2e1 + t123 / 0.2e1;
t465 = t111 / 0.2e1 + t122 / 0.2e1;
t463 = t646 + t614;
t462 = -t193 / 0.2e1 + t597;
t317 = t390 * t347;
t461 = t195 / 0.2e1 + t317 / 0.2e1;
t229 = -Ifges(7,2) * t431 + t575;
t230 = -Ifges(7,1) * t431 - t575;
t460 = -t230 / 0.4e1 + t229 / 0.4e1;
t322 = Ifges(7,4) * t431;
t228 = -Ifges(7,2) * t432 - t322;
t231 = Ifges(7,1) * t432 - t322;
t459 = t231 / 0.4e1 + t228 / 0.4e1;
t455 = t376 - t573;
t450 = mrSges(6,3) * (-t378 / 0.2e1 - t380 / 0.2e1);
t449 = qJD(4) * (t227 + t543);
t448 = t491 * t590;
t351 = Ifges(6,1) * t389 - t577;
t346 = t544 - t548;
t326 = t457 * t386;
t194 = -mrSges(7,1) * t299 - mrSges(7,2) * t301;
t316 = t347 * t386;
t417 = -t316 / 0.2e1 + t194 / 0.2e1 + t335 * t589 + t337 * t591;
t437 = -t261 * t385 + t262 * t389;
t393 = t461 * t312 + t417 * t311 + (-t312 * t496 + t268 * t232 + t269 * t233 + (-t437 + t504) * t311) * t623 + (t110 * t178 + t111 * t179 + t116 * t454 + t117 * t130 - t311 * t326 + t312 * t327) * t621 + t250 * t614 + t130 * t608 + t178 * t604 + t179 * t606 + t336 * t611 + t334 * t610;
t4 = t393 + (t548 / 0.2e1 + t346 / 0.2e1 + (-mrSges(5,1) / 0.2e1 - t447 / 0.2e1 + t612) * t390) * t516 + (t538 / 0.2e1 + t151 * t595) * mrSges(7,3) + t640;
t180 = -Ifges(7,4) * t301 + Ifges(7,2) * t299 + Ifges(7,6) * t390;
t182 = -Ifges(7,1) * t301 + Ifges(7,4) * t299 + Ifges(7,5) * t390;
t294 = Ifges(6,6) * t390 - t386 * t634;
t296 = Ifges(6,5) * t390 - t351 * t386;
t297 = t390 * t351 + t547;
t422 = Ifges(5,4) - t376 / 0.2e1 + t479;
t5 = t262 * t334 + t269 * t335 + t261 * t336 + t268 * t337 + qJ(3) * t346 - t326 * t195 + t327 * t194 + t183 * t599 + t180 * t600 + t182 * t603 + t181 * t601 + t111 * t248 + t117 * t249 + t110 * t250 + t116 * t251 + m(6) * (t261 * t268 + t262 * t269) + m(7) * (t110 * t116 + t111 * t117 - t326 * t327) + (Ifges(7,5) * t603 + Ifges(7,6) * t600 + t392 * t316 + t296 * t588 - t385 * t294 / 0.2e1 - t422 * t390) * t390 + (t392 * t317 + t297 * t589 + t512 / 0.2e1 + t422 * t386 + (-m(6) * t392 ^ 2 - Ifges(5,1) + Ifges(5,2) + Ifges(6,3) + Ifges(7,3)) * t390 + t428) * t386;
t445 = t4 * qJD(1) + t5 * qJD(2);
t318 = t390 * t348;
t404 = t327 * t193 - (t183 / 0.2e1 + t196 / 0.2e1) * t300 + (t111 * mrSges(7,3) - t197 / 0.2e1 + t181 / 0.2e1) * t298 + mrSges(7,3) * t541 + t386 * t492 / 0.2e1;
t10 = m(7) * (t110 * t122 + t111 * t123) + t123 * t249 + t122 * t251 - t262 * t337 + t261 * t335 + (-t392 * t315 + (t473 + t261 * mrSges(6,3) - t297 / 0.2e1 + t318 / 0.2e1) * t385 + (-t546 / 0.2e1 - t262 * mrSges(6,3) - t319 / 0.2e1 - t295 / 0.2e1 + (m(7) * t327 + t195) * pkin(5)) * t389) * t390 + t404;
t398 = (pkin(5) * t311 * t497 + t643 * t454) * t622 + t636 + t232 * t594 + t337 * t610 + (t648 * t622 - t605) * t130;
t427 = t151 * mrSges(7,1) / 0.2e1 - t152 * mrSges(7,2) / 0.2e1;
t405 = t270 * t618 - t271 * mrSges(6,2) / 0.2e1 + (t151 * t388 + t152 * t384) * t489 + t427;
t420 = (-t524 / 0.2e1 - t526 / 0.2e1) * mrSges(7,3);
t6 = t462 * t311 + t420 + (t533 / 0.2e1 - t534 / 0.2e1) * t390 * mrSges(6,3) + t398 + t405;
t444 = -t6 * qJD(1) + t10 * qJD(2);
t20 = t471 - t426;
t424 = t130 * t604 + t636;
t13 = -t521 / 0.2e1 + t420 + t424 + t427;
t15 = t110 * t249 - t111 * t251 + t404;
t443 = -t13 * qJD(1) + t15 * qJD(2);
t381 = t390 ^ 2;
t396 = t462 * t390 + (-t498 / 0.2e1 + t467 + t390 * t450) * t386 + (-t643 * t299 + t648 * t301 - t381 * t582) * t621 + t409;
t19 = t396 + t653;
t442 = t19 * qJD(2);
t358 = t381 * t516;
t379 = t386 ^ 2;
t401 = (t386 * t435 + t358) * t623 + (-t151 * t299 + t152 * t301 + t358) * t621 + (t379 * t516 + t358) * t625;
t406 = (t232 * t389 + t233 * t385) * t624 + (t130 * t432 - t431 * t454) * t622 + m(5) * t470;
t33 = t401 + t406;
t430 = t386 * mrSges(5,1) + t390 * mrSges(5,2) + mrSges(4,3);
t51 = t432 * t249 - t431 * t251 + t511 + t498 + (m(5) + m(4)) * qJ(3) + m(7) * (-t110 * t431 + t111 * t432) + m(6) * (t261 * t389 + t262 * t385) + t430;
t441 = -t33 * qJD(1) + t51 * qJD(2);
t218 = t311 * t312;
t34 = m(6) * (t311 * t438 + t218) + m(7) * (t130 * t179 + t178 * t454 + t218);
t440 = t34 * qJD(1) + t31 * qJD(3);
t434 = t501 + t513;
t429 = -t579 / 0.2e1 - t580 / 0.2e1;
t425 = -t268 * mrSges(6,1) / 0.2e1 + t269 * mrSges(6,2) / 0.2e1;
t423 = m(7) * (t247 * t432 - t431 * t453);
t421 = (t178 * t388 + t179 * t384) * t619;
t394 = -t417 * t390 + (-t510 / 0.2e1 + t499 / 0.2e1 + t461) * t386 + (t437 * t390 + (t436 - 0.2e1 * t496) * t386) * t623 + (-t111 * t298 - t116 * t299 + t117 * t301 + t326 * t390 + t327 * t386 - t541) * t621 + t249 * t603 + t250 * t602 + t251 * t600 + t301 * t608;
t17 = -t423 / 0.2e1 + t394;
t74 = m(6) * t487 * t637 + (-t523 + t525 - t637) * m(7);
t419 = t31 * qJD(1) + t17 * qJD(2) + t74 * qJD(3);
t418 = t593 + t613 + t429;
t395 = t460 * t298 - t459 * t300 + (-t318 / 0.4e1 + t297 / 0.4e1 - pkin(9) * t337 / 0.2e1) * t389 + (t531 / 0.2e1 + t532 / 0.2e1 - t465 * t432 - t466 * t431) * mrSges(7,3) + pkin(4) * t597 - t247 * t604 + t628;
t402 = -t392 * t347 / 0.2e1 + (-t350 / 0.4e1 - t634 / 0.4e1) * t385 + pkin(9) * t450 + (t351 / 0.4e1 - t348 / 0.4e1 + (m(7) * t592 + t612) * pkin(5)) * t389;
t413 = t648 * t247 + t643 * t453;
t1 = t395 + t413 * t621 + (-t502 / 0.2e1 + t468 - t514 / 0.2e1 + 0.2e1 * (-t540 / 0.4e1 - t539 / 0.4e1 + t518 / 0.4e1) * m(7)) * pkin(5) + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t402) * t390 + (t596 - t295 / 0.4e1 + pkin(9) * t594) * t385 + (-0.3e1 / 0.4e1 * t573 + 0.3e1 / 0.4e1 * t376 - t321 / 0.4e1 - t320 / 0.4e1) * t386 + t425 + t631;
t399 = mrSges(7,3) * t431 * t463 + pkin(5) * t520 * t621;
t12 = t418 * t311 - t421 / 0.2e1 + t399 + t426;
t43 = t371 * t226 - (-t230 / 0.2e1 + t229 / 0.2e1) * t432 - (t231 / 0.2e1 + t228 / 0.2e1) * t431;
t24 = -pkin(4) * t347 + t351 * t591 - t509 / 0.2e1 + (t350 / 0.2e1 + t634 / 0.2e1) * t389 + t43 + (m(7) * t371 + t227) * t583;
t412 = t506 * t619;
t414 = -t434 * t619 / 0.2e1 + t493;
t41 = t418 * t390 + t412 / 0.2e1 + t414;
t416 = t12 * qJD(1) + t1 * qJD(2) - t41 * qJD(3) + t24 * qJD(4);
t21 = t471 + t426;
t56 = t472 - t493;
t397 = -(mrSges(7,3) * t609 + t459) * t300 + (t247 * t616 + t460) * t298 + t247 * t605 + t448 + t628;
t9 = t397 - t627;
t415 = t21 * qJD(1) + t9 * qJD(2) + t56 * qJD(3) + t43 * qJD(4);
t408 = (t388 * t606 + t384 * t605 + (t513 / 0.2e1 + t501 / 0.2e1) * mrSges(7,3)) * pkin(5);
t29 = -mrSges(7,1) * t465 + mrSges(7,2) * t466 + t408;
t341 = (mrSges(7,1) * t384 + mrSges(7,2) * t388) * pkin(5);
t39 = mrSges(7,2) * t463;
t62 = (t609 + t453 / 0.2e1) * mrSges(7,2);
t410 = t39 * qJD(1) - t29 * qJD(2) - t62 * qJD(4) + t341 * qJD(5);
t367 = qJ(3) * t515;
t340 = t381 * t480;
t328 = t341 * qJD(6);
t57 = t472 + t493;
t42 = t347 * t587 - t412 / 0.2e1 + t429 * t390 + t414 + t472;
t32 = m(4) * t516 + t401 - t406;
t25 = -t570 / 0.2e1 - t569 / 0.2e1 - t565 / 0.2e1 + t566 / 0.2e1 + t408 + t492;
t23 = t400 + t635;
t18 = t396 - t653;
t16 = t423 / 0.2e1 + t394;
t14 = -t424 + t427 + t632;
t11 = t311 * t593 + t579 * t598 + t520 * t618 + t421 / 0.2e1 + t399 + t20;
t8 = t397 + t627;
t7 = -t233 * t484 / 0.2e1 + t485 * t611 + t315 * t598 - t398 + t405 + t632;
t3 = t151 * t474 - t538 * t616 + t393 + (t544 + t346) * t516 / 0.2e1 + (t548 + (-t447 + t227) * t390) * t470 - t640;
t2 = t395 + t627 + pkin(9) * t467 + pkin(5) * t468 + t389 * t473 + t386 * t479 - t512 / 0.4e1 + t448 + (t502 + t514) * pkin(5) / 0.2e1 + (pkin(5) * t518 + t413) * t621 + t385 * t596 + Ifges(6,3) * t586 + t455 * t590 - t425 + (t539 + t540) * t489 + t402 * t390;
t26 = [t35 * qJD(2) + t34 * qJD(4), t32 * qJD(3) + t3 * qJD(4) + t7 * qJD(5) + t14 * qJD(6) + t517 + (t151 * t251 + t152 * t249 + t270 * t337 + t271 * t335 + ((-mrSges(3,2) + t430) * t391 + (-mrSges(3,1) + mrSges(4,2) + (-t195 - t317) * t390 + (-t379 - t381) * mrSges(5,3)) * t387) * t382 + 0.2e1 * (t261 * t270 + t262 * t271 + t340) * t623 + 0.2e1 * (t110 * t151 + t111 * t152 - t327 * t481) * t621 + 0.2e1 * (t379 * t480 + t340 + t367) * t625 + m(4) * (-pkin(2) * t516 + t367)) * qJD(2), qJD(2) * t32 + t488, t3 * qJD(2) + t11 * qJD(5) + t20 * qJD(6) + (-t178 * t432 - t179 * t431) * t542 + t312 * t449 - t311 * t641 + ((-pkin(9) * t311 * t490 - pkin(4) * t312) * t623 + (t178 * t453 + t179 * t247 + t312 * t371) * t621) * t626 + t440, t7 * qJD(2) + t11 * qJD(4) + (-t233 * mrSges(6,1) - t232 * mrSges(6,2) + (-t130 * t388 + t384 * t454) * t619 + t40) * qJD(5) + t650, t14 * qJD(2) + t20 * qJD(4) + t40 * qJD(5) + t650; -qJD(3) * t33 + qJD(4) * t4 - qJD(5) * t6 - qJD(6) * t13 - t517, qJD(3) * t51 + qJD(4) * t5 + qJD(5) * t10 + qJD(6) * t15, m(7) * (t299 * t431 + t522) * qJD(3) + t16 * qJD(4) + t18 * qJD(5) + t23 * qJD(6) + t441, t16 * qJD(3) + t2 * qJD(5) + t8 * qJD(6) + t445 + (-Ifges(5,6) * t390 + t294 * t588 + t296 * t591 + t371 * t194 + t182 * t595 - t326 * t227 + pkin(4) * t316 + t231 * t599 + t229 * t601 + t247 * t248 + t453 * t250 - mrSges(5,2) * t496 + m(7) * (t116 * t453 + t117 * t247 - t326 * t371) - t116 * t550 + (m(6) * t436 + t499 - t510) * pkin(9) + (-Ifges(5,5) + t350 * t589 + t509 / 0.2e1 + (-m(6) * pkin(4) + t543) * t392) * t386 + (Ifges(6,5) * t385 + Ifges(7,5) * t432 + Ifges(6,6) * t389) * t586 - (Ifges(7,6) * t586 + t180 / 0.2e1 + t117 * mrSges(7,3)) * t431 + t436 * mrSges(6,3)) * qJD(4), t18 * qJD(3) + t2 * qJD(4) + (-t262 * mrSges(6,1) - t261 * mrSges(6,2) - Ifges(6,5) * t506 - Ifges(6,6) * t497 + t492 - t565 + t566) * qJD(5) + t25 * qJD(6) + (m(7) * (t122 * t388 + t123 * t384) + t434 * mrSges(7,3)) * t571 + t444, t23 * qJD(3) + t8 * qJD(4) + t25 * qJD(5) + (t492 - t569 - t570) * qJD(6) + t443; qJD(2) * t33 + t488, qJD(4) * t17 + qJD(5) * t19 + qJD(6) * t22 - t441, t74 * qJD(4), t42 * qJD(5) + t57 * qJD(6) + (t298 * t431 + t300 * t432) * t542 + t390 * t641 + t386 * t449 + ((t490 * t581 - t584) * t623 + (t371 * t386 - t531 - t532) * t621) * t626 + t419, t42 * qJD(4) + (-t447 * t386 + t456 + (-t299 * t384 - t301 * t388) * t619) * qJD(5) + t647 + t442, t57 * qJD(4) + qJD(5) * t456 + t537 + t647; -qJD(2) * t4 + qJD(5) * t12 + qJD(6) * t21 - t440, -qJD(3) * t17 + qJD(5) * t1 + qJD(6) * t9 - t445, -qJD(5) * t41 + qJD(6) * t56 - t419, qJD(5) * t24 + qJD(6) * t43 (-pkin(9) * t447 + t455 + t59) * qJD(5) + t649 + (m(7) * (-t247 * t388 + t384 * t453) - t433 * mrSges(7,3)) * t571 + t416, t59 * qJD(5) + t415 + t649; qJD(2) * t6 - qJD(4) * t12 - qJD(6) * t39, -qJD(3) * t19 - qJD(4) * t1 + qJD(6) * t29 - t444, t41 * qJD(4) - t442, qJD(6) * t62 - t416, -t328, -t328 - t410; t13 * qJD(2) - t21 * qJD(4) + t39 * qJD(5), -qJD(3) * t22 - qJD(4) * t9 - qJD(5) * t29 - t443, -t56 * qJD(4) - t537, -qJD(5) * t62 - t415, t410, 0;];
Cq  = t26;
