% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:25
% EndTime: 2019-03-09 10:45:46
% DurationCPUTime: 11.14s
% Computational Cost: add. (22212->515), mult. (41393->682), div. (0->0), fcn. (46936->8), ass. (0->296)
t378 = sin(qJ(2));
t380 = cos(qJ(2));
t599 = -t380 * pkin(2) - t378 * qJ(3);
t340 = -pkin(1) + t599;
t313 = t380 * pkin(3) - t340;
t557 = cos(qJ(4));
t459 = t380 * t557;
t556 = sin(qJ(4));
t333 = -t378 * t556 - t459;
t273 = -pkin(4) * t333 + t313;
t458 = t380 * t556;
t334 = t378 * t557 - t458;
t375 = sin(pkin(10));
t516 = cos(pkin(10));
t601 = t333 * t375 + t516 * t334;
t602 = -t516 * t333 + t375 * t334;
t647 = m(6) * t273 + mrSges(6,1) * t602 + mrSges(6,2) * t601;
t331 = -t375 * t557 - t516 * t556;
t577 = pkin(7) - pkin(8);
t443 = t577 * t378;
t283 = -t557 * t443 + t577 * t458;
t230 = t334 * qJ(5) + t283;
t589 = t443 * t556 + t577 * t459;
t612 = t333 * qJ(5) + t589;
t652 = t230 * t516 + t375 * t612;
t506 = t652 * t331;
t137 = t230 * t375 - t516 * t612;
t330 = t375 * t556 - t516 * t557;
t510 = t137 * t330;
t670 = t510 - t506;
t430 = Ifges(5,5) * t333 - Ifges(5,6) * t334;
t622 = t589 * mrSges(5,1);
t635 = t283 * mrSges(5,2);
t658 = t652 * mrSges(6,2);
t379 = cos(qJ(6));
t523 = t379 * mrSges(7,1);
t377 = sin(qJ(6));
t528 = t377 * mrSges(7,2);
t431 = t523 - t528;
t665 = t137 * t431;
t667 = t137 * mrSges(6,1);
t669 = t430 - t622 + t635 + t658 + t665 + t667;
t345 = Ifges(7,5) * t377 + Ifges(7,6) * t379;
t500 = t601 * t345;
t668 = -t500 / 0.4e1 + t622 / 0.2e1 - t635 / 0.2e1 - t658 / 0.2e1 - t665 / 0.2e1 - t667 / 0.2e1;
t381 = -pkin(2) - pkin(3);
t338 = -qJ(3) * t556 + t557 * t381;
t337 = -pkin(4) + t338;
t339 = qJ(3) * t557 + t381 * t556;
t490 = t375 * t339;
t266 = t337 * t516 - t490;
t263 = pkin(5) - t266;
t666 = t137 * t263;
t664 = t137 * t652;
t522 = t379 * mrSges(7,2);
t529 = t377 * mrSges(7,1);
t343 = t522 + t529;
t163 = t343 * t601;
t619 = t377 * t602;
t167 = -mrSges(7,2) * t601 + mrSges(7,3) * t619;
t618 = t379 * t602;
t170 = mrSges(7,1) * t601 + mrSges(7,3) * t618;
t565 = t333 / 0.2e1;
t569 = t602 / 0.2e1;
t366 = Ifges(7,6) * t377;
t545 = Ifges(7,5) * t379;
t598 = -t366 + t545;
t606 = t598 / 0.2e1;
t116 = pkin(5) * t602 - pkin(9) * t601 + t273;
t62 = t116 * t379 + t137 * t377;
t620 = t343 * t602;
t63 = t116 * t377 - t137 * t379;
t648 = Ifges(5,1) - Ifges(5,2);
t649 = 2 * Ifges(5,4);
t661 = -t137 * t163 + t63 * t167 + t62 * t170 - t652 * t620 + (t273 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t602 + (Ifges(7,3) + Ifges(6,2)) * t569 + (t606 - Ifges(6,4) / 0.2e1) * t601) * t601 + t313 * (mrSges(5,1) * t334 + mrSges(5,2) * t333) - (-t648 * t333 + t334 * t649) * t334 / 0.2e1 + (t333 * t649 + t648 * t334) * t565;
t311 = t516 * t339;
t267 = t375 * t337 + t311;
t660 = -t137 * t266 + t267 * t652;
t659 = t137 * t516 - t375 * t652;
t655 = t377 * t652;
t654 = t379 * t652;
t558 = t379 / 0.2e1;
t560 = t377 / 0.2e1;
t651 = t167 * t560 + t170 * t558;
t373 = t377 ^ 2;
t374 = t379 ^ 2;
t471 = t373 + t374;
t636 = mrSges(7,3) * t471;
t646 = m(5) * t313 - mrSges(5,1) * t333 + mrSges(5,2) * t334;
t561 = -t377 / 0.2e1;
t548 = Ifges(7,4) * t379;
t348 = -Ifges(7,2) * t377 + t548;
t96 = Ifges(7,6) * t601 - t348 * t602;
t549 = Ifges(7,4) * t377;
t350 = Ifges(7,1) * t379 - t549;
t99 = Ifges(7,5) * t601 - t350 * t602;
t643 = -t558 * t99 - t561 * t96;
t641 = -t377 / 0.4e1;
t640 = -t379 / 0.4e1;
t616 = t431 * t601;
t639 = t616 / 0.2e1;
t527 = t377 * mrSges(7,3);
t468 = t601 * t527;
t168 = -mrSges(7,2) * t602 - t468;
t478 = t379 * t168;
t521 = t379 * mrSges(7,3);
t171 = mrSges(7,1) * t602 - t521 * t601;
t486 = t377 * t171;
t413 = t486 / 0.2e1 - t478 / 0.2e1;
t537 = t602 * mrSges(6,3);
t631 = t537 / 0.2e1 + t413;
t628 = t283 * t556 + t557 * t589;
t627 = -t379 / 0.2e1;
t626 = -t601 / 0.2e1;
t625 = t601 / 0.2e1;
t570 = -t602 / 0.2e1;
t289 = t331 * t431;
t617 = (-t545 / 0.2e1 + t366 / 0.2e1) * t602;
t535 = t601 * mrSges(6,3);
t613 = -t163 / 0.2e1 - t535 / 0.2e1;
t611 = t602 * t636;
t610 = -pkin(5) * t601 - pkin(9) * t602;
t605 = Ifges(3,4) - Ifges(4,5);
t264 = -pkin(9) + t267;
t604 = t264 * t471;
t269 = t338 * t516 - t490;
t440 = t471 * t269;
t77 = m(7) * (-0.1e1 + t471) * t331 * t330;
t603 = t77 * qJD(3);
t398 = -mrSges(5,1) * t556 + t331 * mrSges(6,1) - mrSges(5,2) * t557 + (mrSges(6,2) - t636) * t330;
t597 = t398 + t289;
t461 = -t527 / 0.2e1;
t526 = t377 * t63;
t596 = mrSges(7,3) * t526 / 0.2e1 + t63 * t461;
t349 = Ifges(7,1) * t377 + t548;
t473 = t379 * t349;
t347 = Ifges(7,2) * t379 + t549;
t482 = t377 * t347;
t595 = t473 / 0.4e1 - t482 / 0.4e1;
t479 = t379 * t167;
t487 = t377 * t170;
t594 = -t479 / 0.2e1 + t487 / 0.2e1;
t476 = t379 * t171;
t488 = t377 * t168;
t593 = t488 + t476;
t555 = pkin(4) * t334;
t143 = t555 - t610;
t72 = t143 * t379 + t655;
t73 = t143 * t377 - t654;
t426 = -t72 * t377 + t73 * t379;
t412 = -t167 * t627 - t170 * t560;
t588 = t377 * (t349 / 0.4e1 + t348 / 0.4e1) - t379 * (t350 / 0.4e1 - t347 / 0.4e1);
t447 = t348 / 0.2e1 + t349 / 0.2e1;
t448 = t347 / 0.2e1 - t350 / 0.2e1;
t587 = t448 * t377 - t447 * t379;
t165 = t347 * t601;
t166 = t349 * t601;
t546 = Ifges(7,5) * t602;
t100 = t350 * t601 + t546;
t481 = t379 * t100;
t542 = Ifges(7,6) * t602;
t97 = t348 * t601 + t542;
t525 = t377 * t97;
t563 = t343 / 0.2e1;
t586 = t652 * t563 + t602 * t598 / 0.4e1 + t165 * t640 + t481 / 0.4e1 + t166 * t641 - t525 / 0.4e1;
t585 = m(5) / 0.2e1;
t584 = -m(6) / 0.2e1;
t583 = m(6) / 0.2e1;
t582 = -m(7) / 0.2e1;
t581 = m(7) / 0.2e1;
t580 = m(6) * pkin(4);
t579 = -mrSges(7,1) / 0.2e1;
t578 = mrSges(7,2) / 0.2e1;
t575 = t620 / 0.2e1;
t568 = -t263 / 0.2e1;
t567 = -t269 / 0.2e1;
t457 = t516 * pkin(4);
t357 = -t457 - pkin(5);
t562 = t357 / 0.2e1;
t554 = pkin(4) * t375;
t547 = Ifges(6,5) * t602;
t543 = Ifges(6,6) * t601;
t364 = t380 * qJ(3);
t329 = t378 * t381 + t364;
t278 = t329 - t555;
t393 = Ifges(6,1) * t625 - t525 / 0.2e1 + t481 / 0.2e1 + (t606 - Ifges(6,4)) * t602;
t342 = -t380 * mrSges(4,1) - t378 * mrSges(4,3);
t444 = m(4) * t340 + t342;
t117 = t278 + t610;
t64 = t117 * t379 - t655;
t65 = t117 * t377 + t654;
t3 = (-mrSges(3,2) * pkin(1) - mrSges(4,3) * t340 + t605 * t380) * t380 + (-mrSges(3,1) * pkin(1) + mrSges(4,1) * t340 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t380 - t605 * t378) * t378 + t65 * t168 + t64 * t171 + m(7) * (t62 * t64 + t63 * t65 + t664) + t444 * (pkin(2) * t378 - t364) + (Ifges(6,1) * t569 + Ifges(6,4) * t625 + t643) * t601 + t646 * t329 + (t273 * mrSges(6,2) + t393) * t602 - t661 + t647 * t278;
t530 = t3 * qJD(1);
t524 = t377 * t99;
t520 = t379 * t96;
t5 = -t652 * t616 + t63 * t171 + ((t100 / 0.2e1 - t165 / 0.2e1 + t546 / 0.2e1) * t377 + (t166 / 0.2e1 + t97 / 0.2e1 + t63 * mrSges(7,3) + t542 / 0.2e1) * t379) * t601 + (-t168 - t468) * t62;
t519 = t5 * qJD(1);
t428 = t377 * t62 - t379 * t63;
t508 = t652 * t601;
t14 = (t163 + t535) * t601 - (t478 - t486 - t537) * t602 + m(7) * (t428 * t602 + t508) + m(6) * (t137 * t602 + t508);
t515 = qJD(1) * t14;
t24 = (m(7) * (t379 * t62 + t526) - t444 + t593 + t646 + t647) * t378;
t514 = qJD(1) * t24;
t268 = t338 * t375 + t311;
t507 = t652 * t268;
t253 = t602 * mrSges(6,2);
t391 = -t253 / 0.2e1 + (-t266 * t601 - t267 * t602) * t583 + (t263 * t601 - t602 * t604) * t581 + t611 / 0.2e1;
t392 = t278 * t583 + (t377 * t65 + t379 * t64) * t581 + mrSges(6,2) * t569 - t651;
t16 = 0.2e1 * t625 * mrSges(6,1) + t391 - t392 + t639;
t505 = t16 * qJD(1);
t356 = pkin(9) + t554;
t439 = t471 * t356;
t469 = t580 / 0.2e1;
t389 = (t357 * t601 - t439 * t602) * t581 - t616 / 0.2e1 + (-t375 * t602 - t516 * t601) * t469 - t611 / 0.2e1;
t395 = (t377 * t73 + t379 * t72) * t581 + t334 * t469 + t651;
t19 = -mrSges(6,1) * t601 + t253 + t389 - t395;
t504 = t19 * qJD(1);
t502 = t601 * t330;
t499 = t602 * t331;
t498 = t263 * t343;
t417 = t522 / 0.2e1 + t529 / 0.2e1;
t409 = t417 * t602;
t27 = t409 + t413;
t497 = t27 * qJD(1);
t496 = t330 * t616;
t495 = t330 * t343;
t492 = t357 * t343;
t394 = (t471 * t499 + t502) * t581 + (t499 + t502) * t583;
t406 = -(m(7) * t471 + m(6)) * t378 / 0.2e1;
t44 = t394 + t406;
t472 = t44 * qJD(1);
t470 = mrSges(6,3) * t554;
t465 = mrSges(6,3) * t625;
t463 = mrSges(6,3) * t570;
t460 = t521 / 0.2e1;
t454 = t619 / 0.2e1;
t450 = -t618 / 0.2e1;
t441 = t471 * t330;
t436 = mrSges(6,3) * t457;
t433 = mrSges(7,3) * (-t374 / 0.2e1 - t373 / 0.2e1);
t427 = t377 * t64 - t379 * t65;
t425 = t330 * t631 + t613 * t331;
t382 = t524 / 0.4e1 - t547 / 0.2e1 + t520 / 0.4e1 + (t507 - t660) * t584 - t543 / 0.2e1 - t620 * t568 + (t507 - t666) * t582 + t73 * t460 + t72 * t461 + t266 * t463 + t267 * t465 + t613 * t268 + (t426 * t582 + t594) * t264 - t595 * t602 + (-t137 * t584 - t428 * t582 + t631) * t269 - t668;
t383 = t137 * t357 * t581 + Ifges(6,5) * t569 + t620 * t562 + t99 * t641 + t96 * t640 - t659 * t469 + t64 * t461 + t65 * t460 + (Ifges(6,6) + t470) * t625 + (-t427 * t581 - t412) * t356 + (-t436 / 0.2e1 + t595) * t602 + t668;
t2 = t383 + t382;
t401 = -t339 * mrSges(5,1) - t269 * mrSges(6,2) - t338 * mrSges(5,2) + (-mrSges(6,1) - t431) * t268;
t29 = mrSges(7,3) * t440 - m(7) * (t263 * t268 + t264 * t440) - m(6) * (-t266 * t268 + t267 * t269) + t401;
t424 = -t2 * qJD(1) - t29 * qJD(2);
t40 = -mrSges(4,3) - m(4) * qJ(3) - m(7) * (-t263 * t331 - t264 * t441) - m(6) * (t266 * t331 - t267 * t330) - m(5) * (-t338 * t556 + t339 * t557) + t597;
t386 = (t331 * t427 + t510) * t581 + t670 * t583 + t628 * t585;
t387 = -m(5) * t628 / 0.2e1 + t670 * t584 + (t330 * t428 - t506) * t582;
t9 = (t163 / 0.2e1 + (t626 + t625) * mrSges(6,3) + t412) * t331 + (t575 + (t569 + t570) * mrSges(6,3) - t413) * t330 + t386 + t387;
t423 = t9 * qJD(1) + t40 * qJD(2);
t13 = t425 + (t463 - t620 / 0.2e1 + (-t137 + t428) * t581) * t330 + (t465 + t594 + (-t652 - t426) * t581) * t331;
t4 = -t273 * t253 + t73 * t168 + t72 * t171 - t393 * t602 + (Ifges(6,1) * t570 + Ifges(6,4) * t626 - t643) * t601 + m(7) * (t62 * t72 + t63 * t73 - t664) + t647 * t555 + t661;
t422 = t4 * qJD(1) + t13 * qJD(3);
t118 = -t498 - t587;
t414 = -t488 / 0.2e1 - t476 / 0.2e1;
t384 = t414 * t264 + (t264 * t433 + t588) * t601 + t263 * t639 - t586 + t596;
t403 = Ifges(7,3) * t625 + t578 * t65 + t579 * t64;
t7 = t384 + t403 + t617;
t420 = t7 * qJD(1) + t118 * qJD(2);
t188 = (t563 + t417) * t330;
t410 = (-t528 / 0.2e1 + t523 / 0.2e1) * t378;
t21 = -t496 / 0.2e1 + t410 + (t433 * t601 + t414) * t331;
t419 = -t21 * qJD(1) - t188 * qJD(2);
t411 = t482 / 0.2e1 - t473 / 0.2e1;
t408 = t417 * t330;
t388 = -t289 / 0.2e1 + ((t266 - t269) * t331 + (-t267 + t268) * t330) * t583 + ((-t263 - t440) * t331 + (t268 - t604) * t330) * t581;
t402 = m(7) * (-t357 * t331 - t356 * t441);
t405 = (-t330 * t375 + t331 * t516) * t580;
t390 = t402 / 0.2e1 + t289 / 0.2e1 + t405 / 0.2e1;
t25 = t388 - t390 - t398;
t407 = t13 * qJD(1) + t25 * qJD(2) + t603;
t404 = Ifges(7,3) * t626 + t578 * t73 + t579 * t72;
t385 = t414 * t356 + (t356 * t433 - t588) * t601 + t616 * t562 + t586 + t596;
t11 = t385 + t404 - t617;
t174 = -t492 + t587;
t187 = (-t343 / 0.2e1 + t417) * t330;
t47 = (t568 + t562) * t343 + (mrSges(7,2) * t567 + t447) * t379 + (mrSges(7,1) * t567 - t448) * t377;
t399 = t11 * qJD(1) - t47 * qJD(2) - t187 * qJD(3) - t174 * qJD(4);
t190 = -t495 / 0.2e1 + t408;
t189 = t495 / 0.2e1 + t408;
t48 = t348 * t627 + t350 * t561 + t498 / 0.2e1 - t492 / 0.2e1 - t417 * t269 + t411;
t43 = t394 - t406;
t28 = t409 - t413;
t26 = t388 + t390;
t22 = t496 / 0.2e1 + t410 + (t601 * t636 + t593) * t331 / 0.2e1;
t20 = t389 + t395;
t15 = mrSges(6,1) * t626 + (t431 / 0.2e1 + mrSges(6,1) / 0.2e1) * t601 + t391 + t392;
t12 = t13 * qJD(4);
t10 = Ifges(7,5) * t450 + Ifges(7,6) * t454 + t385 - t404;
t8 = (m(4) * pkin(7) + mrSges(4,2)) * t380 + (mrSges(6,3) * t569 + t575) * t330 + (mrSges(6,3) * t626 + t412) * t331 + t386 - t387 + t425 + (t334 * t556 + 0.2e1 * t557 * t565) * mrSges(5,3);
t6 = t384 - Ifges(7,6) * t619 / 0.2e1 + Ifges(7,5) * t618 / 0.2e1 - t403;
t1 = -t382 + t383 - t430;
t17 = [qJD(2) * t3 + qJD(3) * t24 + qJD(4) * t4 + qJD(5) * t14 - qJD(6) * t5, t8 * qJD(3) + t1 * qJD(4) + t15 * qJD(5) + t6 * qJD(6) + t530 + (t263 * t620 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t380 + (t96 / 0.2e1 - t65 * mrSges(7,3) - t264 * t167) * t379 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t378 + (t99 / 0.2e1 + t64 * mrSges(7,3) + t264 * t170) * t377 - (-t345 / 0.2e1 + Ifges(6,6) - t267 * mrSges(6,3)) * t601 + 0.2e1 * t660 * t583 + 0.2e1 * (t283 * t339 + t338 * t589) * t585 + 0.2e1 * (-t264 * t427 + t666) * t581 + (-t266 * mrSges(6,3) - Ifges(6,5) + t411) * t602 + (m(4) * t599 - t380 * mrSges(3,1) + t378 * mrSges(3,2) + t342) * pkin(7) + (t333 * t338 + t334 * t339) * mrSges(5,3) + t669) * qJD(2), qJD(2) * t8 + qJD(5) * t43 + qJD(6) * t22 + t12 + t514, t1 * qJD(2) + (t659 * t580 + t520 / 0.2e1 + t524 / 0.2e1 + t500 / 0.2e1 - t547 - t543 + t602 * t436 + t347 * t454 + t349 * t450 - t601 * t470 + (-m(7) * t137 - t620) * t357 + (m(7) * t426 + t479 - t487) * t356 + t426 * mrSges(7,3) + t669) * qJD(4) + t20 * qJD(5) + t10 * qJD(6) + t422, qJD(2) * t15 + qJD(3) * t43 + qJD(4) * t20 + qJD(6) * t28 + t515, -t519 + t6 * qJD(2) + t22 * qJD(3) + t10 * qJD(4) + t28 * qJD(5) + (-mrSges(7,1) * t63 - mrSges(7,2) * t62 - t500) * qJD(6); -qJD(3) * t9 - qJD(4) * t2 + qJD(5) * t16 + qJD(6) * t7 - t530, -qJD(3) * t40 - qJD(4) * t29 + qJD(6) * t118, t26 * qJD(4) + t190 * qJD(6) - t423 + t603, t26 * qJD(3) + (t401 + (m(7) * t357 - t516 * t580) * t268 + (m(7) * t439 + t375 * t580 + t636) * t269) * qJD(4) + t48 * qJD(6) + t424, t505, t190 * qJD(3) + t48 * qJD(4) + (-t264 * t431 - t598) * qJD(6) + t420; qJD(2) * t9 + qJD(5) * t44 - qJD(6) * t21 + t12 - t514, qJD(4) * t25 - qJD(6) * t188 + t423, t77 * qJD(4) (t402 + t405 + t597) * qJD(4) + t189 * qJD(6) + t407, t472, t189 * qJD(4) + qJD(6) * t289 + t419; qJD(2) * t2 + qJD(5) * t19 + qJD(6) * t11 - t422, -qJD(3) * t25 - qJD(6) * t47 - t424, -qJD(6) * t187 - t407, -t174 * qJD(6), t504 (-t356 * t431 + t598) * qJD(6) + t399; -qJD(2) * t16 - qJD(3) * t44 - qJD(4) * t19 - qJD(6) * t27 - t515, -t505, -t472, -t504, 0, -qJD(6) * t343 - t497; -qJD(2) * t7 + qJD(3) * t21 - qJD(4) * t11 + qJD(5) * t27 + t519, qJD(3) * t188 + qJD(4) * t47 - t420, qJD(4) * t187 - t419, -t399, t497, 0;];
Cq  = t17;
