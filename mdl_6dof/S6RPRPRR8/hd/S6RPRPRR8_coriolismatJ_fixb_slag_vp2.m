% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:50
% EndTime: 2019-03-09 03:58:05
% DurationCPUTime: 8.58s
% Computational Cost: add. (23925->549), mult. (44092->734), div. (0->0), fcn. (49429->8), ass. (0->299)
t376 = sin(pkin(10));
t379 = sin(qJ(3));
t533 = cos(pkin(10));
t589 = cos(qJ(3));
t352 = t376 * t379 - t533 * t589;
t595 = -t352 / 0.2e1;
t353 = t376 * t589 + t379 * t533;
t377 = sin(qJ(6));
t378 = sin(qJ(5));
t380 = cos(qJ(5));
t588 = cos(qJ(6));
t618 = -t377 * t378 + t588 * t380;
t643 = t618 * t353;
t649 = -t643 / 0.2e1;
t410 = t377 * t380 + t378 * t588;
t620 = t410 * t353;
t648 = t410 * t643 - t618 * t620;
t647 = t643 / 0.2e1;
t592 = t410 / 0.2e1;
t646 = -t618 / 0.2e1;
t635 = -t620 / 0.2e1;
t483 = t379 * pkin(3) + qJ(2);
t279 = pkin(4) * t353 + pkin(8) * t352 + t483;
t381 = -pkin(1) - pkin(7);
t361 = (-qJ(4) + t381) * t379;
t457 = t589 * t381;
t362 = -qJ(4) * t589 + t457;
t619 = t533 * t361 + t376 * t362;
t162 = t380 * t279 - t378 * t619;
t508 = t352 * t380;
t142 = pkin(9) * t508 + t162;
t116 = pkin(5) * t353 + t142;
t163 = t279 * t378 + t380 * t619;
t509 = t352 * t378;
t143 = pkin(9) * t509 + t163;
t493 = t377 * t143;
t73 = t116 * t588 - t493;
t83 = t142 * t588 - t493;
t574 = t73 - t83;
t374 = t378 ^ 2;
t375 = t380 ^ 2;
t478 = t374 + t375;
t645 = t478 * mrSges(6,3) * t352;
t367 = pkin(3) * t376 + pkin(8);
t578 = pkin(9) + t367;
t338 = t578 * t378;
t339 = t578 * t380;
t278 = -t377 * t338 + t339 * t588;
t300 = Ifges(7,5) * t618 - Ifges(7,6) * t410;
t435 = -t588 * t338 - t339 * t377;
t48 = -t278 * mrSges(7,1) - t435 * mrSges(7,2) + t300;
t644 = t48 * qJD(6);
t343 = Ifges(7,4) * t618;
t302 = Ifges(7,1) * t410 + t343;
t642 = -Ifges(7,2) * t410 + t302 + t343;
t414 = Ifges(7,5) * t647 + Ifges(7,6) * t635;
t472 = t589 * pkin(3);
t281 = -t352 * pkin(4) + t353 * pkin(8) + t472;
t292 = t361 * t376 - t533 * t362;
t164 = t380 * t281 + t292 * t378;
t504 = t353 * t380;
t118 = -pkin(5) * t352 + pkin(9) * t504 + t164;
t165 = t378 * t281 - t292 * t380;
t505 = t353 * t378;
t144 = pkin(9) * t505 + t165;
t75 = t118 * t588 - t377 * t144;
t76 = t377 * t118 + t144 * t588;
t613 = Ifges(7,3) * t595 - t76 * mrSges(7,2) / 0.2e1 + t75 * mrSges(7,1) / 0.2e1 - t414;
t255 = t618 * t352;
t258 = t410 * t352;
t148 = -mrSges(7,1) * t255 + mrSges(7,2) * t258;
t641 = t148 / 0.2e1;
t191 = mrSges(7,2) * t352 + mrSges(7,3) * t620;
t640 = t191 / 0.2e1;
t638 = t353 / 0.4e1;
t636 = -t478 / 0.2e1;
t634 = m(5) * t483;
t456 = t588 * t143;
t74 = t377 * t116 + t456;
t82 = -t377 * t142 - t456;
t573 = t74 + t82;
t552 = t258 * mrSges(7,3);
t192 = -mrSges(7,2) * t353 + t552;
t598 = t435 / 0.2e1;
t631 = t192 * t598;
t517 = t278 * t410;
t439 = -mrSges(7,1) * t643 + mrSges(7,2) * t620;
t628 = qJD(6) * t439;
t298 = mrSges(7,1) * t410 + mrSges(7,2) * t618;
t627 = t298 * qJD(6);
t556 = t255 * mrSges(7,3);
t194 = mrSges(7,1) * t353 + t556;
t626 = t192 * t646 + t194 * t592;
t538 = t380 * mrSges(6,1);
t540 = t378 * mrSges(6,2);
t364 = -t538 + t540;
t624 = t364 / 0.2e1;
t617 = t352 * t376 + t533 * t353;
t416 = -t164 * t378 + t165 * t380;
t542 = t410 * mrSges(7,2);
t544 = t618 * mrSges(7,1);
t479 = t544 / 0.2e1 - t542 / 0.2e1;
t553 = t258 * mrSges(7,1);
t557 = t255 * mrSges(7,2);
t481 = t553 / 0.2e1 + t557 / 0.2e1;
t541 = t410 * mrSges(7,3);
t460 = t541 / 0.2e1;
t543 = t618 * mrSges(7,3);
t462 = -t543 / 0.2e1;
t616 = t255 * t460 + t258 * t462 - t626;
t287 = t353 * mrSges(6,1) + mrSges(6,3) * t508;
t485 = t380 * t287;
t285 = -mrSges(6,2) * t353 + mrSges(6,3) * t509;
t490 = t378 * t285;
t615 = -t485 / 0.2e1 - t490 / 0.2e1 + t645 / 0.2e1;
t614 = t192 * t635 + t194 * t649 + t352 * t641;
t496 = t410 * t255;
t500 = t618 * t258;
t612 = (-t496 / 0.2e1 + t500 / 0.2e1) * mrSges(7,3) + t626;
t555 = t255 * Ifges(7,4);
t130 = t258 * Ifges(7,2) + t353 * Ifges(7,6) - t555;
t250 = Ifges(7,4) * t258;
t131 = -t255 * Ifges(7,1) + t353 * Ifges(7,5) + t250;
t151 = Ifges(7,2) * t255 + t250;
t152 = Ifges(7,1) * t258 + t555;
t475 = pkin(5) * t509;
t218 = t292 - t475;
t564 = Ifges(7,4) * t410;
t301 = Ifges(7,2) * t618 + t564;
t368 = -pkin(3) * t533 - pkin(4);
t585 = t380 * pkin(5);
t363 = t368 - t585;
t424 = Ifges(7,1) * t618 - t564;
t611 = t363 * t641 + t300 * t638 + t218 * t298 / 0.2e1 + t642 * t258 / 0.4e1 + (t151 + t131) * t618 / 0.4e1 + (t152 / 0.4e1 - t130 / 0.4e1) * t410 + (-t424 / 0.4e1 + t301 / 0.4e1) * t255;
t308 = t352 ^ 2;
t609 = t353 ^ 2;
t608 = m(5) / 0.2e1;
t607 = m(6) / 0.2e1;
t606 = -m(7) / 0.2e1;
t605 = m(7) / 0.2e1;
t604 = m(5) * pkin(3);
t603 = m(7) * pkin(5);
t150 = -t553 - t557;
t601 = t150 / 0.2e1;
t600 = -t192 / 0.2e1;
t599 = t194 / 0.2e1;
t597 = -t278 / 0.2e1;
t593 = t618 / 0.2e1;
t591 = t378 / 0.2e1;
t590 = t380 / 0.2e1;
t587 = pkin(5) * t377;
t586 = pkin(5) * t378;
t584 = t73 * mrSges(7,2);
t583 = t74 * mrSges(7,1);
t580 = t82 * mrSges(7,1);
t579 = t83 * mrSges(7,2);
t572 = m(7) * qJD(2);
t571 = m(7) * qJD(4);
t566 = Ifges(6,4) * t378;
t565 = Ifges(6,4) * t380;
t563 = Ifges(6,5) * t353;
t561 = Ifges(6,2) * t378;
t560 = Ifges(6,6) * t353;
t554 = t620 * mrSges(7,1);
t551 = t258 * t73;
t550 = t643 * mrSges(7,2);
t149 = -t550 - t554;
t193 = -mrSges(7,1) * t352 + mrSges(7,3) * t643;
t219 = -pkin(5) * t505 + t619;
t428 = t378 * mrSges(6,1) + t380 * mrSges(6,2);
t273 = t428 * t353;
t274 = t352 * t428;
t284 = mrSges(6,2) * t352 + mrSges(6,3) * t505;
t286 = -mrSges(6,1) * t352 + mrSges(6,3) * t504;
t420 = Ifges(6,5) * t380 - Ifges(6,6) * t378;
t415 = -Ifges(5,4) + t420;
t421 = -Ifges(7,4) * t643 + Ifges(7,2) * t620;
t425 = -Ifges(7,1) * t643 + Ifges(7,4) * t620;
t440 = t353 * mrSges(5,1) - t352 * mrSges(5,2);
t441 = -t352 * mrSges(5,1) - t353 * mrSges(5,2);
t3 = (t353 * t415 + t414) * t353 - m(7) * (t218 * t219 + t73 * t75 + t74 * t76) - m(6) * (t162 * t164 + t163 * t165 + t292 * t619) - t483 * t441 + t130 * t635 + (-Ifges(7,5) * t255 + Ifges(7,6) * t258 - t415 * t352 + (-t375 * Ifges(6,1) - Ifges(5,1) + Ifges(5,2) + Ifges(6,3) + Ifges(7,3) + (-t561 + 0.2e1 * t565) * t378) * t353) * t352 + (mrSges(4,2) * qJ(2) - Ifges(4,4) * t379) * t379 - t162 * t286 - t164 * t287 - t163 * t284 - t165 * t285 + t131 * t647 - t218 * t149 - t219 * t150 - t74 * t191 - t76 * t192 - t73 * t193 - t75 * t194 - t258 * t421 / 0.2e1 + t255 * t425 / 0.2e1 + ((Ifges(4,1) - Ifges(4,2)) * t379 - qJ(2) * mrSges(4,1) + Ifges(4,4) * t589) * t589 + t619 * t274 + (-t634 - t440) * t472 + t292 * t273;
t545 = t3 * qJD(1);
t539 = t379 * mrSges(4,1);
t480 = Ifges(7,5) * t258 + Ifges(7,6) * t255;
t394 = t218 * t148 + (t151 / 0.2e1 + t131 / 0.2e1) * t258 + (t130 / 0.2e1 + t74 * mrSges(7,3) - t152 / 0.2e1) * t255 - mrSges(7,3) * t551 + t353 * t480 / 0.2e1;
t412 = -t565 + (-Ifges(6,1) + Ifges(6,2)) * t378;
t4 = t82 * t194 + t83 * t192 + m(7) * (t73 * t82 + t74 * t83) + t162 * t285 - t163 * t287 + ((t292 * mrSges(6,2) - mrSges(6,3) * t162 + Ifges(6,4) * t509 + t563) * t378 + (-t292 * mrSges(6,1) + t163 * mrSges(6,3) + t560 + t412 * t352 + (-m(7) * t218 - t150) * pkin(5)) * t380) * t352 + t394;
t537 = t4 * qJD(1);
t9 = t73 * t192 - t74 * t194 + t394;
t536 = t9 * qJD(1);
t393 = (-t620 ^ 2 - t643 ^ 2 - t308) * t605 + (-t478 * t609 - t308) * t607 + (-t308 - t609) * t608;
t401 = -m(5) / 0.2e1 + m(6) * t636 + (t410 ^ 2 + t618 ^ 2) * t606;
t43 = t393 + t401;
t532 = qJD(1) * t43;
t387 = (-t308 * t585 - t573 * t620 - t574 * t643) * t605 + t308 * t624 + t615 * t353;
t476 = t603 / 0.2e1;
t396 = -t540 / 0.2e1 + t538 / 0.2e1 + (t377 * t410 + t588 * t618) * t476 + t479;
t464 = t552 / 0.2e1;
t466 = t556 / 0.2e1;
t10 = -t464 * t620 - t466 * t643 - t387 + t396 - t614;
t531 = t10 * qJD(1);
t413 = t554 / 0.2e1 + t550 / 0.2e1;
t450 = t505 / 0.2e1;
t390 = (-t377 * t643 + t588 * t620) * t476 + mrSges(6,1) * t450 + mrSges(6,2) * t504 / 0.2e1 + t413;
t486 = t380 * t285;
t488 = t378 * t287;
t411 = t488 / 0.2e1 - t486 / 0.2e1;
t397 = (-t410 * t574 + t573 * t618) * t605 - t411;
t12 = t390 - t397 + t612;
t530 = t12 * qJD(1);
t465 = -t552 / 0.2e1;
t467 = -t556 / 0.2e1;
t418 = -t465 * t620 - t467 * t643 + t614;
t15 = t418 - t479;
t529 = t15 * qJD(1);
t417 = t162 * t378 - t163 * t380;
t515 = t292 * t352;
t18 = -t643 * t192 + t620 * t194 + (mrSges(5,3) * t353 - t486 + t488) * t353 + (mrSges(5,3) * t352 - t150 + t274) * t352 + m(7) * (-t218 * t352 + t620 * t73 - t643 * t74) + m(6) * (t353 * t417 - t515) + m(5) * (-t353 * t619 - t515);
t526 = t18 * qJD(1);
t21 = t413 + t612;
t525 = t21 * qJD(1);
t521 = t255 * t377;
t518 = t435 * t618;
t402 = mrSges(4,2) * t589 + t440 + t539;
t30 = t410 * t192 + t618 * t194 + t490 + t485 + mrSges(3,3) + (m(3) + m(4)) * qJ(2) + m(7) * (t410 * t74 + t618 * t73) + m(6) * (t162 * t380 + t163 * t378) + t634 + t402;
t514 = t30 * qJD(1);
t512 = t352 * t298;
t307 = t352 * t353;
t511 = t352 * t368;
t506 = t353 * t364;
t492 = t377 * t194;
t489 = t378 * t286;
t487 = t380 * t284;
t477 = mrSges(7,3) * t587;
t117 = -t496 + t500;
t474 = t117 * qJD(3) * t605;
t473 = pkin(5) * t588;
t214 = mrSges(7,3) * t517;
t461 = t543 / 0.2e1;
t459 = -t541 / 0.2e1;
t458 = t512 / 0.2e1 - (t459 + t460) * t643;
t455 = t588 * t258;
t299 = t542 - t544;
t453 = t299 * t595;
t452 = t509 / 0.2e1;
t451 = t508 / 0.2e1;
t449 = -t504 / 0.2e1;
t437 = t478 * t367;
t436 = t255 * t477;
t434 = mrSges(7,3) * t473;
t433 = t74 * t460;
t432 = t473 / 0.2e1;
t429 = t258 * t434;
t427 = Ifges(6,1) * t380 - t566;
t426 = Ifges(6,1) * t378 + t565;
t423 = -t561 + t565;
t422 = Ifges(6,2) * t380 + t566;
t408 = m(7) * (t517 + t518);
t385 = (-t353 * t437 - t511) * t607 + (-t278 * t643 - t352 * t363 + t435 * t620) * t605 + t453 + t364 * t595 + (t352 * t533 - t353 * t376) * t604 / 0.2e1 + t620 * t459 - t643 * t461 + t353 * mrSges(6,3) * t636;
t388 = (t164 * t380 + t165 * t378) * t607 + (t410 * t76 + t618 * t75) * t605 + t193 * t593 + t191 * t592 + t284 * t591 + t286 * t590 + t472 * t608;
t19 = t385 - t388 - t441;
t405 = -t19 * qJD(1) + t117 * t572 / 0.2e1;
t17 = (-t83 / 0.2e1 + t73 / 0.2e1) * mrSges(7,2) + (t82 / 0.2e1 + t74 / 0.2e1) * mrSges(7,1) + (t588 * t600 + t492 / 0.2e1 + (-t521 / 0.2e1 + t455 / 0.2e1) * mrSges(7,3)) * pkin(5);
t360 = (mrSges(7,1) * t377 + mrSges(7,2) * t588) * pkin(5);
t49 = (-t435 / 0.2e1 + t598) * mrSges(7,2) + (t597 + t278 / 0.2e1) * mrSges(7,1);
t404 = t17 * qJD(1) - t49 * qJD(3) + t360 * qJD(5);
t45 = m(7) * (-t255 * t643 - t258 * t620 + t307) + m(6) * (-t478 + 0.1e1) * t307;
t384 = (t149 / 0.2e1 - t273 / 0.2e1 + t411) * t352 + (t601 - t274 / 0.2e1 + t487 / 0.2e1 - t489 / 0.2e1) * t353 + ((t292 + t416) * t353 + (t619 + t417) * t352) * t607 + (t218 * t353 + t219 * t352 - t255 * t74 - t620 * t75 + t643 * t76 + t551) * t605 + t643 * t640 + t255 * t600 + t193 * t635 + t258 * t599;
t8 = -t408 / 0.2e1 + t384;
t403 = -t8 * qJD(1) - t45 * qJD(2) - t117 * t571 / 0.2e1;
t383 = ((t218 * t378 - t363 * t508) * pkin(5) + t573 * t435 - t574 * t278) * t605 - t278 * t599 + t631 + t292 * t428 / 0.2e1 + t420 * t638 - t378 * (-t352 * t423 + t560) / 0.4e1 + t380 * (-t352 * t427 + t563) / 0.4e1 + t511 * t624 + t586 * t601 + t426 * t452 + t423 * t509 / 0.4e1 - t427 * t508 / 0.4e1 + t422 * t451 + t73 * t462 + t82 * t459 + t83 * t461 + t453 * t585 + t615 * t367;
t386 = Ifges(6,3) * t595 + t164 * mrSges(6,1) / 0.2e1 - t165 * mrSges(6,2) / 0.2e1 + (t377 * t76 + t588 * t75) * t476 + Ifges(6,5) * t449 + Ifges(6,6) * t450 + t587 * t640 + t193 * t432 + t613;
t2 = t278 * t467 + t435 * t464 - t383 + t386 + t433 - t611;
t392 = t428 * t595 + t475 * t606;
t395 = (t455 - t521) * t476 + mrSges(6,1) * t452 + mrSges(6,2) * t451 + t481;
t32 = -t512 / 0.2e1 + t392 + t395;
t391 = t214 - t363 * t298 + t301 * t592 - t410 * t424 / 0.2e1 + t642 * t646;
t34 = t374 * Ifges(6,4) - t368 * t428 + (t518 - t517) * mrSges(7,3) + t412 * t380 + t391 + (-m(7) * t363 - t299) * t586 - t435 * t543;
t400 = t2 * qJD(1) + t32 * qJD(2) + t34 * qJD(3);
t35 = t458 - t481;
t37 = t391 - t214;
t398 = t278 * t466 + t435 * t465 + t74 * t459 + t611;
t389 = t194 * t597 + t398 + t433 + t631;
t6 = t389 - t613;
t399 = -t6 * qJD(1) - t35 * qJD(2) + t37 * qJD(3);
t355 = t360 * qJD(6);
t42 = t393 - t401;
t36 = t458 + t481;
t33 = -t392 + t395 + t458;
t22 = t413 + t616;
t20 = t385 + t388;
t16 = t418 + t479;
t14 = -t583 / 0.2e1 - t584 / 0.2e1 + t192 * t432 + t436 / 0.2e1 - pkin(5) * t492 / 0.2e1 - t429 / 0.2e1 + t580 / 0.2e1 - t579 / 0.2e1 + t480;
t13 = t390 + t397 + t616;
t11 = t387 + t396 + t418;
t7 = t408 / 0.2e1 + t384;
t5 = t389 + t613;
t1 = t398 + t383 + t386;
t23 = [qJD(2) * t30 - qJD(3) * t3 + qJD(4) * t18 + qJD(5) * t4 + qJD(6) * t9, t7 * qJD(3) + t42 * qJD(4) + t11 * qJD(5) + t16 * qJD(6) + t648 * t572 + t514, -t545 + t7 * qJD(2) + (m(7) * (t219 * t363 + t278 * t76 + t435 * t75) + t435 * t193 + t617 * mrSges(5,3) * pkin(3) + t426 * t449 + t422 * t450 - Ifges(4,6) * t589 - t381 * t539 - mrSges(4,2) * t457 + (-Ifges(6,6) * t352 - t353 * t423) * t590 + (-Ifges(6,5) * t352 - t353 * t427) * t591 + (-Ifges(7,5) * t352 + t425) * t592 + (-Ifges(7,6) * t352 + t421) * t593 - t75 * t541 + t76 * t543 + t416 * mrSges(6,3) + (m(6) * t416 + t487 - t489) * t367 + (Ifges(6,5) * t378 + Ifges(7,5) * t410 + Ifges(6,6) * t380 + Ifges(7,6) * t618) * t595 + t620 * t301 / 0.2e1 + (m(6) * t368 - t533 * t604 - mrSges(5,1) + t364) * t619 - (t376 * t604 - mrSges(5,2)) * t292 - Ifges(4,5) * t379 - t368 * t273 + t363 * t149 - Ifges(5,5) * t353 + Ifges(5,6) * t352 + t219 * t299 + t278 * t191 + t302 * t649) * qJD(3) + t20 * qJD(4) + t1 * qJD(5) + t5 * qJD(6), t42 * qJD(2) + t20 * qJD(3) + t13 * qJD(5) + t22 * qJD(6) - t648 * t571 + t526, t537 + t11 * qJD(2) + t1 * qJD(3) + t13 * qJD(4) + (Ifges(6,5) * t509 + Ifges(6,6) * t508 + t580 - t579 + t436 - t429 + (t377 * t83 + t588 * t82) * t603 - t162 * mrSges(6,2) - t163 * mrSges(6,1) + t480) * qJD(5) + t14 * qJD(6), t536 + t16 * qJD(2) + t5 * qJD(3) + t22 * qJD(4) + t14 * qJD(5) + (t480 - t583 - t584) * qJD(6); qJD(3) * t8 + qJD(4) * t43 - qJD(5) * t10 + qJD(6) * t15 - t514, t45 * qJD(3) (m(7) * (-t255 * t278 + t258 * t435 + t353 * t363) + t353 * t299 - t255 * t543 - t258 * t541 + m(6) * (-t352 * t437 + t353 * t368) + t506 - t617 * t604 - t402 - t645) * qJD(3) + t33 * qJD(5) + t36 * qJD(6) - t403, t474 + t532, -t531 + t33 * qJD(3) + (m(7) * (-t473 * t643 - t587 * t620) + t506 + t439) * qJD(5) + t628, t36 * qJD(3) + qJD(5) * t439 + t529 + t628; -qJD(2) * t8 + qJD(4) * t19 - qJD(5) * t2 + qJD(6) * t6 + t545, -t32 * qJD(5) + t35 * qJD(6) + t403, -qJD(5) * t34 - qJD(6) * t37, -t405 (-t618 * t434 + (-t278 * t588 + t377 * t435) * t603 - t410 * t477 + t420 + t364 * t367 + t48) * qJD(5) + t644 - t400, t48 * qJD(5) - t399 + t644; -qJD(2) * t43 - qJD(3) * t19 - qJD(5) * t12 - qJD(6) * t21 - t526, t474 - t532, t405, 0, -t530 - t627 + (-t298 - t428 + (-t410 * t473 + t587 * t618) * m(7)) * qJD(5), -qJD(5) * t298 - t525 - t627; qJD(2) * t10 + qJD(3) * t2 + qJD(4) * t12 - qJD(6) * t17 - t537, t32 * qJD(3) + t531, t49 * qJD(6) + t400, t530, -t355, -t355 - t404; -qJD(2) * t15 - qJD(3) * t6 + qJD(4) * t21 + qJD(5) * t17 - t536, -t35 * qJD(3) - t529, -t49 * qJD(5) + t399, t525, t404, 0;];
Cq  = t23;
