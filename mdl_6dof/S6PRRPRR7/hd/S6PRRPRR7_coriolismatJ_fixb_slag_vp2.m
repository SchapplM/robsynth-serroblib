% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:25
% EndTime: 2019-03-08 22:30:40
% DurationCPUTime: 9.38s
% Computational Cost: add. (15494->644), mult. (33428->885), div. (0->0), fcn. (34256->10), ass. (0->344)
t399 = sin(qJ(6));
t400 = sin(qJ(5));
t403 = cos(qJ(6));
t404 = cos(qJ(5));
t509 = t403 * t404;
t338 = t399 * t400 - t509;
t451 = t399 * t404 + t403 * t400;
t239 = -mrSges(7,1) * t338 - mrSges(7,2) * t451;
t685 = t239 / 0.2e1;
t401 = sin(qJ(3));
t405 = cos(qJ(3));
t651 = -t401 * pkin(3) + qJ(4) * t405;
t684 = m(5) * t651;
t398 = sin(pkin(6));
t402 = sin(qJ(2));
t406 = cos(qJ(2));
t504 = t404 * t406;
t280 = (-t400 * t402 + t401 * t504) * t398;
t513 = t400 * t406;
t281 = (t401 * t513 + t402 * t404) * t398;
t154 = t280 * t403 - t281 * t399;
t155 = t280 * t399 + t281 * t403;
t638 = m(7) * pkin(5);
t673 = mrSges(7,1) / 0.2e1;
t682 = -mrSges(7,2) / 0.2e1;
t656 = t154 * t673 + t155 * t682;
t683 = -t280 * mrSges(6,1) / 0.2e1 + t281 * mrSges(6,2) / 0.2e1 - (t154 * t403 + t155 * t399) * t638 / 0.2e1 - t656;
t308 = t451 * t405;
t620 = t308 / 0.2e1;
t407 = -pkin(3) - pkin(9);
t512 = t400 * t407;
t354 = -t400 * pkin(10) + t512;
t355 = (-pkin(10) + t407) * t404;
t471 = -t354 * t399 + t403 * t355;
t681 = -t471 / 0.2e1;
t394 = t400 ^ 2;
t396 = t404 ^ 2;
t498 = t394 + t396;
t680 = mrSges(6,3) * t498;
t554 = t404 * mrSges(6,2);
t561 = t400 * mrSges(6,1);
t361 = t554 + t561;
t474 = -mrSges(7,1) * t451 + t338 * mrSges(7,2);
t427 = -t361 + t474;
t650 = mrSges(5,3) - t427;
t549 = qJ(4) * t401;
t328 = t407 * t405 - pkin(2) - t549;
t633 = pkin(4) + pkin(8);
t367 = t633 * t401;
t231 = -t328 * t400 + t404 * t367;
t514 = t400 * t405;
t202 = pkin(10) * t514 + t231;
t232 = t328 * t404 + t367 * t400;
t505 = t404 * t405;
t203 = -pkin(10) * t505 + t232;
t543 = t203 * t399;
t103 = t202 * t403 - t543;
t188 = pkin(5) * t401 + t202;
t94 = t188 * t403 - t543;
t679 = t103 - t94;
t521 = t398 * t402;
t550 = cos(pkin(6));
t324 = t401 * t521 - t550 * t405;
t246 = t324 * t404 + t398 * t513;
t247 = -t324 * t400 + t398 * t504;
t129 = t246 * t399 - t247 * t403;
t472 = t403 * t246 + t247 * t399;
t39 = -t129 * mrSges(7,1) - t472 * mrSges(7,2);
t678 = t39 * qJD(6);
t241 = -Ifges(7,5) * t451 + Ifges(7,6) * t338;
t250 = t354 * t403 + t355 * t399;
t53 = -t250 * mrSges(7,1) - t471 * mrSges(7,2) + t241;
t677 = t53 * qJD(6);
t519 = t398 * t406;
t479 = -t519 / 0.2e1;
t676 = t479 * t684;
t675 = qJD(6) * t474;
t307 = t338 * t405;
t591 = Ifges(7,4) * t308;
t174 = Ifges(7,2) * t307 + t401 * Ifges(7,6) - t591;
t190 = -mrSges(7,1) * t308 + mrSges(7,2) * t307;
t194 = Ifges(7,1) * t307 + t591;
t590 = Ifges(7,4) * t338;
t243 = -Ifges(7,2) * t451 - t590;
t244 = -Ifges(7,1) * t451 + t590;
t368 = t633 * t405;
t323 = pkin(5) * t505 + t368;
t382 = pkin(5) * t400 + qJ(4);
t674 = (t243 / 0.4e1 - t244 / 0.4e1) * t308 + (t174 / 0.4e1 - t194 / 0.4e1) * t338 + t382 * t190 / 0.2e1 + t323 * t685;
t566 = t307 * mrSges(7,3);
t262 = -mrSges(7,2) * t401 + t566;
t628 = t262 / 0.2e1;
t564 = t308 * mrSges(7,3);
t264 = mrSges(7,1) * t401 + t564;
t625 = t264 / 0.2e1;
t602 = t405 / 0.2e1;
t600 = m(7) * t382;
t562 = t451 * mrSges(7,3);
t542 = t203 * t403;
t102 = -t202 * t399 - t542;
t95 = t188 * t399 + t542;
t670 = t102 + t95;
t626 = -t264 / 0.2e1;
t669 = t250 * t626;
t294 = Ifges(7,4) * t307;
t176 = -Ifges(7,1) * t308 + Ifges(7,5) * t401 + t294;
t193 = Ifges(7,2) * t308 + t294;
t668 = t176 + t193;
t666 = t401 ^ 2 + t405 ^ 2;
t665 = t338 * t399 + t451 * t403;
t336 = pkin(9) * t401 - t651;
t342 = t404 * t368;
t189 = pkin(5) * t405 + t342 + (-pkin(10) * t401 - t336) * t400;
t238 = t404 * t336 + t400 * t368;
t511 = t401 * t404;
t207 = pkin(10) * t511 + t238;
t100 = t189 * t403 - t207 * t399;
t101 = t189 * t399 + t207 * t403;
t309 = t451 * t401;
t619 = t309 / 0.2e1;
t515 = t400 * t401;
t306 = -t399 * t515 + t401 * t509;
t624 = t306 / 0.2e1;
t450 = Ifges(7,5) * t619 + Ifges(7,6) * t624;
t465 = Ifges(7,3) * t602 + t100 * t673 + t101 * t682 + t450;
t662 = t250 * t620 + t307 * t681;
t661 = -0.2e1 * t405;
t660 = -Ifges(4,4) - Ifges(5,6);
t563 = t338 * mrSges(7,3);
t360 = t405 * mrSges(5,2) - t401 * mrSges(5,3);
t655 = -t405 * mrSges(4,1) + t401 * mrSges(4,2) + t360;
t654 = t602 * t680;
t507 = t404 * t280;
t517 = t400 * t281;
t452 = t507 + t517;
t652 = -t100 * t338 + t101 * t451;
t364 = t401 * mrSges(4,1) + t405 * mrSges(4,2);
t648 = t338 * t628 + t451 * t625;
t520 = t398 * t405;
t325 = t550 * t401 + t402 * t520;
t171 = t338 * t325;
t172 = t451 * t325;
t448 = t171 * t673 + t172 * mrSges(7,2) / 0.2e1;
t618 = t325 / 0.2e1;
t631 = -t472 / 0.2e1;
t647 = t472 * t628 + (t129 * t620 + t307 * t631) * mrSges(7,3) + t190 * t618;
t645 = 0.2e1 * t325;
t644 = 2 * qJD(3);
t643 = m(5) / 0.2e1;
t642 = -m(6) / 0.2e1;
t641 = m(6) / 0.2e1;
t640 = -m(7) / 0.2e1;
t639 = m(7) / 0.2e1;
t637 = mrSges(6,1) / 0.2e1;
t635 = -Ifges(6,2) / 0.2e1;
t634 = -t95 / 0.2e1;
t632 = qJ(4) / 0.2e1;
t630 = t474 / 0.2e1;
t261 = -mrSges(7,2) * t405 + mrSges(7,3) * t306;
t629 = t261 / 0.2e1;
t263 = mrSges(7,1) * t405 - mrSges(7,3) * t309;
t627 = t263 / 0.2e1;
t623 = t307 / 0.2e1;
t621 = -t308 / 0.2e1;
t617 = -t338 / 0.2e1;
t616 = -t451 / 0.2e1;
t493 = mrSges(6,3) * t514;
t558 = t401 * mrSges(6,1);
t344 = t493 + t558;
t614 = t344 / 0.2e1;
t551 = t405 * mrSges(6,2);
t345 = mrSges(6,3) * t511 - t551;
t613 = -t345 / 0.2e1;
t492 = mrSges(6,3) * t505;
t556 = t401 * mrSges(6,2);
t346 = -t492 - t556;
t612 = t346 / 0.2e1;
t555 = t404 * mrSges(6,1);
t560 = t400 * mrSges(6,2);
t359 = t555 - t560;
t611 = t359 / 0.2e1;
t610 = -t361 / 0.2e1;
t609 = -t400 / 0.2e1;
t607 = t400 / 0.2e1;
t606 = t403 / 0.2e1;
t605 = -t404 / 0.2e1;
t603 = t404 / 0.2e1;
t599 = pkin(5) * t404;
t598 = t94 * mrSges(7,2);
t597 = t95 * mrSges(7,1);
t594 = Ifges(6,1) * t400;
t593 = Ifges(6,4) * t400;
t592 = Ifges(6,4) * t404;
t589 = Ifges(6,5) * t400;
t588 = Ifges(6,5) * t401;
t587 = Ifges(6,6) * t401;
t586 = Ifges(6,6) * t404;
t584 = pkin(5) * qJD(5);
t581 = t102 * mrSges(7,1);
t580 = t103 * mrSges(7,2);
t567 = t307 * mrSges(7,1);
t565 = t308 * mrSges(7,2);
t557 = t401 * mrSges(5,2);
t553 = t405 * mrSges(6,1);
t546 = t100 * t403;
t544 = t101 * t399;
t541 = t246 * t404;
t540 = t247 * t400;
t447 = mrSges(7,1) * t619 + mrSges(7,2) * t624;
t28 = (t307 * t617 + t308 * t616) * mrSges(7,3) + t447 + t648;
t537 = t28 * qJD(2);
t225 = t324 * t325;
t30 = m(7) * (t129 * t172 - t171 * t472 - t225) + m(6) * (-t225 + (-t540 + t541) * t325);
t536 = t30 * qJD(1);
t489 = t405 * t519;
t260 = t325 * t489;
t445 = (t324 * t401 - t521) * t398;
t31 = m(7) * (t129 * t155 + t154 * t472 + t260) + m(6) * (t246 * t280 - t247 * t281 + t260) + m(5) * (t325 * t520 + t445) * t406 + (t445 * t406 + t260) * m(4);
t535 = t31 * qJD(1);
t533 = t323 * t404;
t532 = t324 * qJ(4);
t531 = t325 * t404;
t530 = t338 * t154;
t528 = t338 * t309;
t526 = t451 * t155;
t524 = t451 * t306;
t518 = t399 * t308;
t516 = t400 * t344;
t510 = t403 * t307;
t192 = -t565 - t567;
t508 = t404 * t192;
t506 = t404 * t346;
t503 = t404 * t407;
t500 = Ifges(7,5) * t307 + Ifges(7,6) * t308;
t499 = t666 * pkin(8) * t519;
t496 = m(6) / 0.4e1 + m(7) / 0.4e1;
t495 = pkin(5) * t514;
t485 = t563 / 0.2e1;
t484 = t562 / 0.2e1;
t483 = -t562 / 0.2e1;
t482 = -t94 / 0.2e1 + t103 / 0.2e1;
t481 = t634 - t102 / 0.2e1;
t480 = t239 * t618;
t477 = t631 + t472 / 0.2e1;
t363 = -t405 * mrSges(5,3) - t557;
t475 = -t364 / 0.2e1 - t363 / 0.2e1;
t473 = t498 * t325;
t466 = mrSges(6,3) * (-t394 / 0.2e1 - t396 / 0.2e1);
t366 = Ifges(6,1) * t404 - t593;
t463 = t592 + t594;
t365 = -t400 * Ifges(6,2) + t592;
t462 = Ifges(6,2) * t404 + t593;
t461 = t586 + t589;
t460 = -pkin(3) * t405 - t549;
t191 = -mrSges(7,1) * t306 + mrSges(7,2) * t309;
t237 = -t336 * t400 + t342;
t322 = (-t599 - t633) * t401;
t326 = t359 * t401;
t327 = t359 * t405;
t343 = -mrSges(6,3) * t515 + t553;
t408 = (-t327 / 0.2e1 - t192 / 0.2e1) * t324 + (-t326 / 0.2e1 + t191 / 0.2e1 + t344 * t603 + t346 * t607) * t325 + (t237 * t246 - t238 * t247 - t324 * t368 + (t231 * t404 + t232 * t400 - t367) * t325) * t641 + (t100 * t472 + t101 * t129 - t171 * t94 + t172 * t95 + t322 * t325 - t323 * t324) * t639 + t472 * t627 + t129 * t629 - t171 * t625 + t172 * t628 + t246 * t343 / 0.2e1 + t247 * t613 - t676;
t410 = (qJ(4) * t489 + t452 * t407) * t642 + (t154 * t471 + t155 * t250 + t382 * t489) * t640 + t676;
t4 = ((mrSges(4,1) / 0.2e1 - mrSges(5,2) / 0.2e1) * t401 + (mrSges(4,2) / 0.2e1 - mrSges(5,3) / 0.2e1 + t610 + t630) * t405 + t475) * t519 + t408 + (-t530 / 0.2e1 + t526 / 0.2e1) * mrSges(7,3) + (t507 / 0.2e1 + t517 / 0.2e1) * mrSges(6,3) + t410;
t173 = Ifges(7,4) * t309 + Ifges(7,2) * t306 + Ifges(7,6) * t405;
t175 = Ifges(7,1) * t309 + Ifges(7,4) * t306 + Ifges(7,5) * t405;
t303 = Ifges(6,6) * t405 + t462 * t401;
t304 = Ifges(6,5) * t405 + t463 * t401;
t356 = -pkin(2) + t460;
t449 = t586 / 0.2e1 + t589 / 0.2e1;
t5 = -t368 * t326 - t651 * t360 - pkin(2) * t364 - t367 * t327 + t231 * t343 + t237 * t344 + t232 * t345 + t238 * t346 + t322 * t192 + t323 * t191 + t176 * t619 + t173 * t623 + t175 * t621 + t174 * t624 + t95 * t261 + t101 * t262 + t94 * t263 + t100 * t264 + (t363 - t684) * t356 + m(7) * (t100 * t94 + t101 * t95 + t322 * t323) + m(6) * (t231 * t237 + t232 * t238 - t367 * t368) + (Ifges(7,5) * t621 + Ifges(7,6) * t623 + t303 * t605 + t304 * t609 + (-t449 - t660) * t405) * t405 + ((t461 + t660) * t401 + (Ifges(6,3) - Ifges(5,3) + Ifges(4,1) + Ifges(5,2) - Ifges(4,2) + Ifges(7,3) + t396 * t635 + (-t592 - t594 / 0.2e1) * t400) * t405 + t450) * t401;
t459 = t4 * qJD(1) + t5 * qJD(2);
t437 = t323 * t190 + t95 * t564 + t401 * t500 / 0.2e1;
t440 = t405 * t361;
t10 = -t102 * t264 - m(7) * (t102 * t94 + t103 * t95) - t103 * t262 + t368 * t440 + (-t174 / 0.2e1 + t194 / 0.2e1) * t308 + (-t193 / 0.2e1 - t176 / 0.2e1 + t94 * mrSges(7,3)) * t307 + ((-Ifges(6,4) * t505 + t588) * t404 + (Ifges(6,4) * t514 - t587 + (-Ifges(6,1) + Ifges(6,2)) * t505 + (m(7) * t323 + t192) * pkin(5)) * t400) * t405 - t437 + (-t493 + t344) * t232 + (-t346 - t492) * t231;
t409 = (-t325 * t495 + t670 * t472) * t639 + t246 * t612 + t247 * t614 + t647 + (t639 * t679 - t625) * t129;
t431 = (-t540 / 0.2e1 + t541 / 0.2e1) * mrSges(6,3);
t7 = (t325 * t610 + t431) * t405 + t409 + t683;
t458 = t7 * qJD(1) - t10 * qJD(2);
t19 = t480 - t448;
t13 = -t264 * t95 + t174 * t620 + t194 * t621 + (t262 - t566) * t94 + t437 + t668 * t623;
t415 = t129 * t626 + t647;
t14 = t415 - t656;
t457 = t14 * qJD(1) + t13 * qJD(2);
t34 = t309 * t264 + m(7) * (-t306 * t95 + t309 * t94) - t306 * t262 + (-t506 + t516 + m(6) * (t231 * t400 - t232 * t404) - m(5) * t356 - t360) * t401;
t417 = t452 * t642 + (t526 - t530) * t640;
t434 = (-t129 * t306 + t309 * t472) * t639;
t435 = (t246 * t400 + t247 * t404) * t641;
t36 = t401 * t435 + t417 + t434;
t456 = t36 * qJD(1) + t34 * qJD(2);
t454 = -t171 * t338 - t172 * t451;
t453 = t237 * t404 + t238 * t400;
t446 = t560 / 0.2e1 - t555 / 0.2e1;
t444 = t263 * t606 + t399 * t629;
t443 = t365 * t603 + t366 * t607;
t442 = -t506 / 0.2e1 + t516 / 0.2e1;
t441 = m(7) * t454;
t436 = (-t171 * t403 + t172 * t399) * t638;
t433 = -t440 / 0.2e1;
t18 = t480 + t448;
t332 = Ifges(7,4) * t451;
t242 = Ifges(7,2) * t338 - t332;
t245 = -Ifges(7,1) * t338 - t332;
t37 = t382 * t239 - (t245 / 0.2e1 + t242 / 0.2e1) * t451 + (-t244 / 0.2e1 + t243 / 0.2e1) * t338;
t420 = t95 * t485 + t401 * t241 / 0.4e1 + (t242 + t245) * t307 / 0.4e1 - t668 * t451 / 0.4e1 + t662 * mrSges(7,3) + t674;
t413 = t471 * t628 + t563 * t634 + t420 + t669;
t9 = t413 - t465;
t432 = t18 * qJD(1) + t9 * qJD(2) + t37 * qJD(3);
t426 = m(7) * (-t670 * t338 + t451 * t679);
t428 = (-t306 * t399 + t309 * t403) * t638 / 0.2e1;
t16 = (t556 / 0.2e1 - t346 / 0.2e1) * t404 + (t558 / 0.2e1 + t614) * t400 + t405 * t466 - t426 / 0.2e1 + t428 + t28;
t430 = t16 * qJD(2);
t147 = t600 + (m(6) + m(5)) * qJ(4) + t650;
t412 = (t528 / 0.2e1 + t524 / 0.2e1) * mrSges(7,3) + t368 * t641 + (-t250 * t306 + t309 * t471 + t323) * t639 - t567 / 0.2e1 - t565 / 0.2e1;
t414 = t261 * t616 + t338 * t627 + t453 * t642 + t652 * t640;
t24 = (t553 / 0.2e1 - t343 / 0.2e1) * t404 + (-t551 / 0.2e1 + t613) * t400 + t412 + t414;
t61 = t441 / 0.2e1 + (-m(6) * t498 / 0.4e1 + t496) * t645;
t429 = -qJD(1) * t61 - qJD(2) * t24 - qJD(3) * t147;
t425 = -t669 + t262 * t681 - t368 * t359 / 0.2e1;
t29 = t307 * t485 + t308 * t484 + t447 - t648;
t424 = t250 * t679 + t670 * t471;
t422 = t237 * t637 - t238 * mrSges(6,2) / 0.2e1 + t465;
t1 = t422 + t424 * t640 + (-t508 / 0.2e1 + 0.2e1 * (-t533 / 0.4e1 + t546 / 0.4e1 + t544 / 0.4e1) * m(7) + t444) * pkin(5) + (Ifges(6,3) / 0.2e1 + t407 * t466 + (mrSges(6,2) * t632 + t366 / 0.4e1 + (Ifges(6,1) / 0.4e1 + t635) * t404) * t404 + (mrSges(6,1) * t632 - 0.3e1 / 0.2e1 * t592 - t365 / 0.4e1 + (Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.2e1) * t400 + (-t474 / 0.2e1 + t600 / 0.2e1) * pkin(5)) * t400) * t405 + (t481 * t338 + t482 * t451 - t662) * mrSges(7,3) - (-t193 / 0.4e1 - t176 / 0.4e1) * t451 + (-t245 / 0.4e1 - t242 / 0.4e1) * t307 + t442 * t407 + (-t241 / 0.4e1 + t461) * t401 + t425 - t674;
t411 = pkin(5) * t531 * t639 - t477 * t562;
t12 = (t685 + t611 + t446) * t325 - t436 / 0.2e1 + t411 + t448;
t27 = qJ(4) * t359 + t462 * t607 + t463 * t605 + t37 - t443 + (-t474 + t600) * t599;
t423 = t12 * qJD(1) - t1 * qJD(2) + t27 * qJD(3);
t416 = (t399 * t626 + t262 * t606 + (-t510 / 0.2e1 + t518 / 0.2e1) * mrSges(7,3)) * pkin(5);
t21 = t481 * mrSges(7,1) + t482 * mrSges(7,2) + t416;
t351 = (mrSges(7,1) * t399 + mrSges(7,2) * t403) * pkin(5);
t38 = t477 * mrSges(7,2);
t418 = -t38 * qJD(1) - t21 * qJD(2) + t351 * qJD(5);
t337 = t351 * qJD(6);
t60 = -t441 / 0.2e1 + t473 * t641 + (t643 + t496) * t645;
t35 = t434 + (m(5) * t519 + t435) * t401 - t417;
t22 = t345 * t607 + t343 * t603 + (m(5) * pkin(8) + mrSges(5,1) - t446) * t405 + t412 - t414;
t20 = -t597 / 0.2e1 - t598 / 0.2e1 + t581 / 0.2e1 - t580 / 0.2e1 + t416 + t500;
t17 = t426 / 0.2e1 + (t554 / 0.2e1 + t561 / 0.2e1) * t401 + t428 + t29 - t442 + t654;
t15 = t415 + t656;
t11 = t436 / 0.2e1 + t531 * t637 + t411 + t19 + (t611 - t560 / 0.2e1) * t325;
t8 = t413 + t465;
t6 = t325 * t433 + t405 * t431 + t409 - t683;
t3 = t154 * t485 + t155 * t483 + t475 * t519 + t408 - t410 - t452 * mrSges(6,3) / 0.2e1 + t364 * t479 + (t650 * t405 + t557) * t519 / 0.2e1;
t2 = qJ(4) * t433 + t424 * t639 + t422 + t420 + t495 * t630 + Ifges(6,3) * t602 + t503 * t612 - t344 * t512 / 0.2e1 - t366 * t505 / 0.2e1 + t365 * t514 / 0.2e1 - t425 + t102 * t485 + t103 * t483 + t94 * t484 - (t463 * t661 + t588) * t400 / 0.4e1 - (t462 * t661 + t587) * t404 / 0.4e1 + t654 * t407 + (t449 - t461 / 0.4e1) * t401 + (t508 / 0.2e1 + t444 + (-t382 * t514 + t533 + t544 + t546) * t639) * pkin(5);
t23 = [t31 * qJD(2) + t30 * qJD(3), t3 * qJD(3) + t35 * qJD(4) + t6 * qJD(5) + t15 * qJD(6) + t535 + (t154 * t264 + t155 * t262 + t280 * t344 + t281 * t346 + ((-mrSges(3,1) + t655) * t402 + (-mrSges(3,2) + (t192 + t327) * t405 + (mrSges(4,3) + mrSges(5,1)) * t666) * t406) * t398 + 0.2e1 * (t154 * t94 + t155 * t95 + t323 * t489) * t639 + 0.2e1 * (t231 * t280 + t232 * t281 + t368 * t489) * t641 + 0.2e1 * (t356 * t521 + t499) * t643 + m(4) * (-pkin(2) * t521 + t499)) * qJD(2), t536 + t3 * qJD(2) + t60 * qJD(4) + t11 * qJD(5) + t19 * qJD(6) + ((-t171 * t471 + t172 * t250 - t324 * t382) * t639 + (t407 * t473 - t532) * t641 + (-pkin(3) * t325 - t532) * t643) * t644 + (t454 * mrSges(7,3) + (-mrSges(4,1) + mrSges(5,2) - t680) * t325 + (mrSges(4,2) - t650) * t324) * qJD(3), qJD(2) * t35 + qJD(3) * t60, t6 * qJD(2) + t11 * qJD(3) + (t247 * mrSges(6,1) - t246 * mrSges(6,2) + (-t129 * t403 + t399 * t472) * t638 + t39) * qJD(5) + t678, t15 * qJD(2) + t19 * qJD(3) + t39 * qJD(5) + t678; qJD(3) * t4 + qJD(4) * t36 + qJD(5) * t7 + qJD(6) * t14 - t535, qJD(3) * t5 + qJD(4) * t34 - qJD(5) * t10 + qJD(6) * t13, t22 * qJD(4) + t2 * qJD(5) + t8 * qJD(6) + ((-qJ(4) * t367 + t407 * t453) * t641 + (t100 * t471 + t101 * t250 + t322 * t382) * t639) * t644 + t459 + (-qJ(4) * t326 + t173 * t616 + t175 * t617 + t382 * t191 - t322 * t474 + t243 * t624 + t245 * t619 + t471 * t263 + t250 * t261 - t367 * t361 + (t304 / 0.2e1 + t407 * t343 - t237 * mrSges(6,3)) * t404 + (-t303 / 0.2e1 + t407 * t345 - t238 * mrSges(6,3)) * t400 + (-pkin(3) * mrSges(5,1) + Ifges(6,5) * t603 + Ifges(7,5) * t617 + Ifges(6,6) * t609 + Ifges(7,6) * t616 - Ifges(5,4) + Ifges(4,5)) * t405 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + t443) * t401 + (m(5) * t460 + t655) * pkin(8) - t652 * mrSges(7,3)) * qJD(3), t22 * qJD(3) + m(7) * (-t524 - t528) * qJD(4) + t17 * qJD(5) + t29 * qJD(6) + t456, t2 * qJD(3) + t17 * qJD(4) + (-t232 * mrSges(6,1) - t231 * mrSges(6,2) - Ifges(6,5) * t505 + Ifges(6,6) * t514 + t500 - t580 + t581) * qJD(5) + t20 * qJD(6) + (m(7) * (t102 * t403 + t103 * t399) + (-t510 + t518) * mrSges(7,3)) * t584 + t458, t8 * qJD(3) + t29 * qJD(4) + t20 * qJD(5) + (t500 - t597 - t598) * qJD(6) + t457; -qJD(2) * t4 + qJD(4) * t61 + qJD(5) * t12 + qJD(6) * t18 - t536, qJD(4) * t24 - qJD(5) * t1 + qJD(6) * t9 - t459, qJD(4) * t147 + qJD(5) * t27 + qJD(6) * t37, -t429 (-mrSges(6,1) * t512 - mrSges(6,2) * t503 - t461 + t53) * qJD(5) + t677 + (m(7) * (-t250 * t403 + t399 * t471) + t665 * mrSges(7,3)) * t584 + t423, t53 * qJD(5) + t432 + t677; -qJD(2) * t36 - qJD(3) * t61, -qJD(3) * t24 - qJD(5) * t16 - qJD(6) * t28 - t456, t429, 0 (-t665 * t638 + t427) * qJD(5) + t675 - t430, qJD(5) * t474 - t537 + t675; -qJD(2) * t7 - qJD(3) * t12 + qJD(6) * t38, qJD(3) * t1 + qJD(4) * t16 + qJD(6) * t21 - t458, -t423, t430, -t337, -t337 - t418; -t14 * qJD(2) - t18 * qJD(3) - t38 * qJD(5), -qJD(3) * t9 + qJD(4) * t28 - qJD(5) * t21 - t457, -t432, t537, t418, 0;];
Cq  = t23;
