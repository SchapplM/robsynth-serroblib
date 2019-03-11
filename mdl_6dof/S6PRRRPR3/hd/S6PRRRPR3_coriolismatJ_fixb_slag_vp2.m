% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:11:12
% EndTime: 2019-03-08 23:11:28
% DurationCPUTime: 8.71s
% Computational Cost: add. (14040->557), mult. (31511->751), div. (0->0), fcn. (33838->10), ass. (0->328)
t348 = sin(qJ(3));
t349 = sin(qJ(2));
t345 = sin(pkin(6));
t351 = cos(qJ(3));
t492 = t345 * t351;
t520 = cos(pkin(6));
t295 = t348 * t520 + t349 * t492;
t347 = sin(qJ(4));
t493 = t345 * t349;
t384 = t348 * t493 - t351 * t520;
t565 = cos(qJ(4));
t370 = t295 * t565 - t347 * t384;
t562 = m(7) * t370;
t311 = t347 * t351 + t348 * t565;
t352 = cos(qJ(2));
t491 = t345 * t352;
t265 = t311 * t491;
t662 = (mrSges(6,2) / 0.2e1 - mrSges(5,1) / 0.2e1) * t265;
t350 = cos(qJ(6));
t525 = t350 * mrSges(7,1);
t346 = sin(qJ(6));
t530 = t346 * mrSges(7,2);
t421 = t525 - t530;
t221 = t421 * t311;
t453 = t565 * t351;
t310 = t347 * t348 - t453;
t500 = t310 * t346;
t225 = mrSges(7,1) * t311 - mrSges(7,3) * t500;
t499 = t310 * t350;
t227 = -t311 * mrSges(7,2) + mrSges(7,3) * t499;
t566 = t350 / 0.2e1;
t567 = t346 / 0.2e1;
t395 = t225 * t566 + t227 * t567;
t382 = -t221 / 0.2e1 + t395;
t581 = -pkin(9) - pkin(8);
t323 = t581 * t348;
t324 = t581 * t351;
t258 = -t565 * t323 - t324 * t347;
t195 = pkin(5) * t311 + t258;
t332 = -pkin(3) * t351 - pkin(2);
t516 = qJ(5) * t311;
t405 = t332 - t516;
t582 = pkin(4) + pkin(10);
t146 = t310 * t582 + t405;
t75 = -t146 * t346 + t195 * t350;
t76 = t146 * t350 + t195 * t346;
t414 = t346 * t76 + t350 * t75;
t402 = -t195 + t414;
t591 = m(7) / 0.2e1;
t661 = (t402 * t591 + t382) * t370;
t561 = pkin(3) * t347;
t329 = qJ(5) + t561;
t392 = t329 * t421;
t394 = qJ(5) * t421;
t554 = Ifges(7,4) * t346;
t418 = Ifges(7,2) * t350 + t554;
t553 = Ifges(7,4) * t350;
t555 = Ifges(7,1) * t346;
t419 = t553 + t555;
t319 = -Ifges(7,2) * t346 + t553;
t320 = Ifges(7,1) * t350 - t554;
t624 = t319 * t566 + t320 * t567;
t427 = -t350 * t419 / 0.2e1 + t418 * t567 - t624;
t660 = t427 + t394 / 0.2e1 + t392 / 0.2e1;
t202 = t295 * t347 + t565 * t384;
t135 = -t346 * t202 + t350 * t491;
t399 = t530 / 0.2e1 - t525 / 0.2e1;
t387 = t399 * t370;
t393 = t370 * t421;
t524 = t350 * mrSges(7,3);
t458 = t524 / 0.2e1;
t459 = -t524 / 0.2e1;
t623 = t458 + t459;
t625 = -t387 + t393 / 0.2e1 + t623 * t135;
t514 = t135 * t346;
t134 = t202 * t350 + t346 * t491;
t515 = t134 * t350;
t409 = -t514 + t515;
t646 = (-t202 + t409) * t562;
t656 = t646 * qJD(1);
t659 = t625 * qJD(6) + t656;
t658 = qJD(3) + qJD(4);
t336 = t346 * mrSges(7,1);
t337 = t350 * mrSges(7,2);
t607 = t337 + t336;
t468 = mrSges(6,3) + t607;
t608 = t347 * t323 - t565 * t324;
t622 = -t310 * pkin(5) + t608;
t653 = t346 * t622;
t652 = t350 * t622;
t223 = pkin(4) * t310 + t405;
t236 = -mrSges(6,2) * t310 - mrSges(6,3) * t311;
t651 = m(6) * t223 + t236;
t343 = t350 ^ 2;
t556 = mrSges(7,3) * t343;
t341 = t346 ^ 2;
t557 = mrSges(7,3) * t341;
t650 = -t556 / 0.2e1 - t557 / 0.2e1;
t475 = t341 + t343;
t432 = t475 * t582;
t634 = qJ(5) * t202;
t649 = -t370 * t432 - t634;
t648 = -pkin(4) * t370 - t634;
t144 = -Ifges(7,6) * t310 + t311 * t418;
t145 = -Ifges(7,5) * t310 + t311 * t419;
t570 = -t310 / 0.2e1;
t377 = -t346 * t144 / 0.2e1 + t145 * t566 + (Ifges(7,5) * t350 - Ifges(7,6) * t346) * t570 + (Ifges(6,4) - Ifges(5,5)) * t310 + (Ifges(6,5) - Ifges(5,6) + t624) * t311;
t632 = t195 * t607;
t637 = t608 * mrSges(6,2);
t638 = t608 * mrSges(5,1);
t639 = t258 * mrSges(6,3);
t640 = t258 * mrSges(5,2);
t645 = t377 + t637 + t640 - t632 - t638 - t639;
t644 = t632 / 0.2e1 - t637 / 0.2e1 + t638 / 0.2e1 + t639 / 0.2e1 - t640 / 0.2e1;
t614 = mrSges(5,3) + mrSges(6,1);
t641 = -Ifges(5,4) - Ifges(6,6);
t603 = mrSges(7,3) * t475;
t636 = mrSges(5,1) + t603;
t635 = qJ(5) * t195;
t633 = t195 * t622;
t631 = t622 * t202;
t630 = t329 * t195;
t594 = m(6) / 0.2e1;
t627 = t491 * t594;
t233 = pkin(4) * t311 + qJ(5) * t310;
t339 = t348 * pkin(3);
t224 = t233 + t339;
t529 = t346 * mrSges(7,3);
t226 = -mrSges(7,1) * t310 - t311 * t529;
t228 = t310 * mrSges(7,2) + t311 * t524;
t580 = t134 / 0.2e1;
t397 = t226 * t580 - t135 * t228 / 0.2e1;
t466 = t348 * t491;
t220 = t421 * t310;
t573 = -t220 / 0.2e1;
t596 = m(5) / 0.2e1;
t307 = t311 * pkin(10);
t147 = t224 + t307;
t77 = -t147 * t346 + t652;
t78 = t147 * t350 + t653;
t626 = -pkin(3) * t466 * t596 - t224 * t627 + (t134 * t77 - t135 * t78 - t631) * t591 - t202 * t573 + t397 + t661;
t621 = -pkin(4) * t608 - qJ(5) * t258;
t470 = t565 * pkin(3);
t331 = -t470 - pkin(4);
t620 = -t329 * t258 + t331 * t608;
t619 = t258 * t347 + t565 * t608;
t266 = -t347 * t466 + t453 * t491;
t618 = t662 + (-mrSges(5,2) / 0.2e1 + t468 / 0.2e1) * t266;
t232 = -mrSges(6,2) * t311 + mrSges(6,3) * t310;
t234 = mrSges(5,1) * t311 - mrSges(5,2) * t310;
t549 = Ifges(7,6) * t350;
t551 = Ifges(7,5) * t346;
t401 = t551 / 0.2e1 + t549 / 0.2e1;
t417 = -t549 - t551;
t585 = Ifges(7,2) / 0.2e1;
t617 = t195 * t220 + (t144 * t566 + t145 * t567 + (-t401 - t641) * t310) * t310 + ((-t417 + t641) * t311 + (-Ifges(7,3) - Ifges(5,1) + Ifges(6,3) + Ifges(5,2) - Ifges(6,2) + t343 * t585 + (t553 + t555 / 0.2e1) * t346) * t310) * t311 - t622 * t221 + t223 * t232 + t75 * t226 + t76 * t228 + t332 * t234;
t583 = -m(7) - m(6);
t616 = t370 / 0.2e1;
t606 = t650 * t310;
t318 = t348 * mrSges(4,1) + t351 * mrSges(4,2);
t169 = t233 + t307;
t82 = -t169 * t346 + t652;
t522 = t350 * t82;
t83 = t169 * t350 + t653;
t527 = t346 * t83;
t412 = t522 + t527;
t605 = -mrSges(4,1) * t351 + mrSges(4,2) * t348;
t593 = m(6) / 0.4e1;
t472 = m(7) / 0.4e1 + t593;
t213 = t265 * t350 - t346 * t493;
t483 = t350 * t213;
t214 = t265 * t346 + t350 * t493;
t490 = t346 * t214;
t602 = (t490 / 0.2e1 + t483 / 0.2e1) * mrSges(7,3);
t601 = (t627 + (t134 * t346 + t135 * t350) * t591) * t311;
t294 = -t583 * qJ(5) + t468;
t460 = -t529 / 0.2e1;
t599 = t213 * t459 + t214 * t460 + t618;
t598 = 0.2e1 * m(7);
t597 = 2 * qJD(4);
t595 = -m(6) / 0.2e1;
t592 = -m(7) / 0.2e1;
t589 = m(5) * pkin(3);
t588 = mrSges(7,1) / 0.2e1;
t587 = mrSges(5,2) / 0.2e1;
t586 = -mrSges(7,2) / 0.2e1;
t584 = -Ifges(7,3) / 0.2e1;
t579 = t622 / 0.2e1;
t578 = t622 / 0.4e1;
t576 = -t202 / 0.2e1;
t575 = -t370 / 0.2e1;
t571 = t225 / 0.2e1;
t569 = t311 / 0.2e1;
t568 = -t607 / 0.2e1;
t564 = m(6) * t233;
t563 = m(6) * t608;
t552 = Ifges(7,5) * t311;
t550 = Ifges(7,6) * t311;
t548 = t370 * mrSges(5,1);
t547 = t370 * mrSges(6,2);
t546 = t202 * mrSges(5,2);
t545 = t202 * mrSges(6,3);
t532 = t310 * mrSges(6,1);
t531 = t311 * mrSges(6,1);
t528 = t346 * t78;
t523 = t350 * t77;
t518 = qJ(5) * t221;
t222 = t607 * t310;
t517 = qJ(5) * t222;
t508 = t202 * t607;
t114 = t370 * t266;
t494 = t345 ^ 2 * t349;
t23 = m(7) * (t134 * t213 - t135 * t214 + t114) + m(4) * (t345 * t348 * t384 + t295 * t492 - t494) * t352 + (m(6) + m(5)) * (t202 * t265 - t352 * t494 + t114);
t505 = t23 * qJD(1);
t503 = t266 * qJ(5);
t501 = t266 * t329;
t498 = t329 * t202;
t497 = t329 * t221;
t496 = t329 * t222;
t386 = (t337 / 0.2e1 + t336 / 0.2e1) * t311;
t481 = t350 * t227;
t489 = t346 * t225;
t396 = t489 / 0.2e1 - t481 / 0.2e1;
t436 = t341 / 0.2e1 + t343 / 0.2e1;
t423 = mrSges(7,3) * t436;
t33 = t310 * t423 + t386 + t396;
t495 = t33 * qJD(2);
t488 = t346 * t228;
t327 = -pkin(10) + t331;
t486 = t346 * t327;
t485 = t346 * t582;
t482 = t350 * t226;
t479 = t350 * t327;
t478 = t582 * t227;
t474 = t348 ^ 2 + t351 ^ 2;
t471 = mrSges(7,3) * t528;
t469 = t329 * t531;
t467 = t565 * mrSges(5,2);
t462 = -t532 / 0.2e1;
t461 = t531 / 0.2e1;
t456 = mrSges(6,2) * t561 + t468 * t470;
t455 = t565 * t370;
t454 = t565 * t329;
t447 = -t491 / 0.2e1;
t216 = t488 / 0.2e1;
t443 = t485 / 0.2e1;
t440 = -t482 / 0.2e1;
t438 = t479 / 0.2e1;
t437 = -t234 / 0.2e1 - t232 / 0.2e1;
t435 = m(7) * t475;
t433 = t475 * t327;
t431 = t310 * t470;
t424 = t435 / 0.4e1;
t408 = t483 + t490;
t359 = (t265 * t331 + t501) * t594 + (t327 * t408 + t501) * t591 + (-t265 * t565 + t266 * t347) * t589 / 0.2e1 + t318 * t447;
t2 = -t359 + t602 + (t318 + t234 + t232) * t447 + t614 * (t311 * t575 + t370 * t569) - t618 + t626;
t235 = mrSges(5,1) * t310 + mrSges(5,2) * t311;
t3 = -pkin(2) * t318 + t77 * t225 + t78 * t227 + (-Ifges(4,4) * t348 + pkin(3) * t235) * t348 + m(7) * (t75 * t77 + t76 * t78 - t633) + m(5) * t332 * t339 + (Ifges(4,4) * t351 + (Ifges(4,1) - Ifges(4,2)) * t348) * t351 + t617 + t651 * t224;
t416 = t2 * qJD(1) + t3 * qJD(2);
t362 = (-pkin(4) * t265 + t503) * t595 + (-t408 * t582 + t503) * t592;
t376 = -t220 * t576 + t397;
t388 = t134 * t82 - t135 * t83 - t631;
t6 = t437 * t491 - t662 + t602 + (t568 + t587 - mrSges(6,3) / 0.2e1) * t266 + t382 * t370 + t447 * t564 + (t370 * t402 + t388) * t591 + t362 + t376;
t7 = t82 * t225 + t83 * t227 + t233 * t236 + m(7) * (t75 * t82 + t76 * t83 - t633) + t223 * t564 + t617;
t415 = t6 * qJD(1) + t7 * qJD(2);
t413 = t523 + t528;
t292 = Ifges(7,5) * t499;
t11 = t75 * t227 - t76 * t225 + t292 * t569 + t622 * t222 + ((-t75 * mrSges(7,3) + t552 / 0.2e1 + Ifges(7,4) * t499) * t350 + (-t76 * mrSges(7,3) - t550 + (-t554 + (Ifges(7,1) - Ifges(7,2)) * t350) * t310) * t346) * t310;
t363 = (t514 / 0.2e1 - t515 / 0.2e1) * t310 * mrSges(7,3) + t227 * t580 + t135 * t571 + t222 * t616;
t400 = t213 * t588 + t214 * t586;
t16 = t363 - t400;
t411 = t16 * qJD(1) + t11 * qJD(2);
t26 = (-t481 + t489 + m(7) * (t346 * t75 - t350 * t76) - t651) * t311;
t374 = t265 * t594 + t408 * t591;
t41 = -t374 + t601;
t410 = qJD(1) * t41 + qJD(2) * t26;
t406 = t258 * t265 + t266 * t608;
t404 = t586 * t78 + t588 * t77;
t403 = t586 * t83 + t588 * t82;
t389 = (t593 + t424) * t370;
t355 = (-t498 + (t202 * t347 + t455) * pkin(3)) * t594 + (pkin(3) * t455 + t409 * t561 - t498) * t591 + t202 * t587 + mrSges(6,3) * t576 + t202 * t568 + mrSges(5,1) * t575 + mrSges(6,2) * t616 + (t331 * t594 + t433 * t591 + t650) * t370;
t358 = t648 * t595 + t649 * t592 + t548 / 0.2e1 - t547 / 0.2e1 - t546 / 0.2e1 + t545 / 0.2e1 + t508 / 0.2e1 + t616 * t603;
t10 = t355 + t358;
t61 = -t456 + (t467 - m(7) * (t347 * t433 + t454) - m(6) * (t331 * t347 + t454)) * pkin(3) + t636 * t561;
t354 = (pkin(3) * t619 + t620) * t594 + (-t630 + t412 * t327 + (t347 * t414 + t565 * t622) * pkin(3)) * t591 - t497 / 0.2e1 - t469 / 0.2e1 + t331 * t462 + t327 * t216 + t226 * t438 + t82 * t459 + t83 * t460 + t470 * t573 - mrSges(6,1) * t431 / 0.2e1 + (t461 + t395) * t561 - t644;
t357 = t621 * t595 + (-t413 * t582 - t635) * t592 + t518 / 0.2e1 + pkin(4) * t462 + qJ(5) * t461 + t228 * t443 - t582 * t440 + t77 * t458 + t471 / 0.2e1 + t644;
t8 = t354 + t357;
t383 = t10 * qJD(1) + t8 * qJD(2) - t61 * qJD(3);
t102 = t392 + t427;
t364 = t552 + (0.3e1 / 0.2e1 * t553 + t319 / 0.4e1 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.4e1) * t346) * t310 + mrSges(7,2) * t579;
t368 = t550 + (-t320 / 0.4e1 + (-Ifges(7,1) / 0.4e1 + t585) * t350) * t310 - t622 * mrSges(7,1) / 0.2e1;
t12 = -t496 / 0.2e1 + (t327 * t423 + t584) * t310 + (-t327 * t227 / 0.2e1 + t368) * t350 + (t327 * t571 + t364) * t346 + t404;
t29 = -t393 / 0.2e1 - t387;
t381 = -t29 * qJD(1) - t12 * qJD(2) + t102 * qJD(3);
t229 = 0.4e1 * t329 * t472 + t468;
t367 = t399 * t310 - t488 / 0.2e1 + t440;
t27 = (t578 - t528 / 0.4e1 - t523 / 0.4e1) * t598 + t367;
t44 = -0.2e1 * t370 * t472 + 0.2e1 * t389;
t380 = qJD(1) * t44 - qJD(2) * t27 - qJD(3) * t229;
t110 = t394 + t427;
t14 = -t517 / 0.2e1 + (-t423 * t582 + t584) * t310 + (t478 / 0.2e1 + t368) * t350 + (-t571 * t582 + t364) * t346 + t403;
t378 = t399 * t561;
t67 = -t378 - t660;
t379 = t14 * qJD(2) + t67 * qJD(3) - t110 * qJD(4);
t136 = (-0.1e1 / 0.2e1 + t436) * m(7) * t561 - t294;
t35 = (t578 - t527 / 0.4e1 - t522 / 0.4e1) * t598 + t367;
t56 = 0.2e1 * (t341 / 0.4e1 + t343 / 0.4e1 - 0.1e1 / 0.4e1) * t562;
t372 = qJD(1) * t56 - qJD(2) * t35 + qJD(3) * t136 - qJD(4) * t294;
t369 = t216 + t482 / 0.2e1 + (-mrSges(6,1) + t399) * t310 + m(7) * t579;
t361 = Ifges(7,3) * t570 - t346 * (t310 * t419 + t552) / 0.4e1 - t350 * (t310 * t418 + t550) / 0.4e1 + t421 * t579 + t623 * t76 + (t320 / 0.2e1 - t418 / 0.4e1) * t499 - (0.2e1 * t319 + t419) * t500 / 0.4e1 + (t401 + t417 / 0.4e1) * t311;
t333 = qJ(5) * t470;
t158 = t562 / 0.2e1;
t137 = (t435 / 0.2e1 + t594) * t561 + t468 + t472 * (0.4e1 * qJ(5) + 0.2e1 * t561);
t68 = -t378 + t660;
t53 = t158 + 0.2e1 * (t594 + t424) * t370;
t45 = m(6) * t616 + t158 + 0.2e1 * t389;
t40 = t601 + t374;
t34 = t386 - t396 + t606;
t28 = t412 * t591 + t369 + t563;
t25 = t413 * t591 + t608 * t594 + t563 / 0.2e1 + t369;
t17 = t363 + t400;
t15 = t517 / 0.2e1 - t478 * t566 + t225 * t443 + t361 + t403 - t606 * t582;
t13 = t496 / 0.2e1 - t225 * t486 / 0.2e1 + t227 * t438 + t361 + t404 + t606 * t327;
t9 = t355 - t358;
t5 = t388 * t591 + t661 + (-t564 / 0.2e1 + t437) * t491 - t362 + t376 + t599;
t4 = t354 - t357 + t377;
t1 = (-t318 / 0.2e1 + t437) * t491 + t359 + t614 * ((t575 + t616) * t311 + (t576 + t202 / 0.2e1) * t310) + t599 + t626;
t18 = [t23 * qJD(2) + t646 * t658, t1 * qJD(3) + t5 * qJD(4) + t40 * qJD(5) + t17 * qJD(6) + t505 + (t213 * t225 + t214 * t227 + (-t614 * t310 - t220) * t266 + 0.2e1 * (t223 * t493 + t406) * t594 + 0.2e1 * (t213 * t75 + t214 * t76 + t266 * t622) * t591 + 0.2e1 * (t332 * t493 + t406) * t596 + t614 * t265 * t311 + ((mrSges(4,3) * t474 - mrSges(3,2)) * t352 + (-mrSges(3,1) + t235 + t236 + t605) * t349 + m(4) * (pkin(8) * t352 * t474 - pkin(2) * t349)) * t345) * qJD(2), t1 * qJD(2) + (-t295 * mrSges(4,1) + t384 * mrSges(4,2) - t508 - t545 + t546 + t547 - t548 + (m(6) * t331 + m(7) * t433 - t565 * t589 - t556 - t557) * t370 + (t583 * t329 - t347 * t589) * t202) * qJD(3) + t9 * qJD(4) + t45 * qJD(5) + t659, t5 * qJD(2) + t9 * qJD(3) + t53 * qJD(5) + (t591 * t649 + t594 * t648) * t597 + ((mrSges(5,2) - t468) * t202 + (mrSges(6,2) - t636) * t370) * qJD(4) + t659, qJD(2) * t40 + qJD(3) * t45 + qJD(4) * t53, t17 * qJD(2) + (mrSges(7,1) * t135 - mrSges(7,2) * t134) * qJD(6) + t658 * t625; qJD(3) * t2 + qJD(4) * t6 + qJD(5) * t41 + qJD(6) * t16 - t505, qJD(3) * t3 + qJD(4) * t7 + qJD(5) * t26 + qJD(6) * t11 (-t469 + t645 + (-t311 * t561 + t431) * mrSges(5,3) - t471 + Ifges(4,5) * t351 - Ifges(4,6) * t348 + t605 * pkin(8) - t331 * t532 - mrSges(7,3) * t523 + m(7) * (t327 * t413 - t630) - t619 * t589 + m(6) * t620 + t226 * t479 + t228 * t486 - t497) * qJD(3) + t4 * qJD(4) + t25 * qJD(5) + t13 * qJD(6) + t416, t4 * qJD(3) + (-mrSges(6,1) * t516 - t412 * mrSges(7,3) + pkin(4) * t532 - t228 * t485 - t482 * t582 - t518 + t645) * qJD(4) + t28 * qJD(5) + t15 * qJD(6) + ((-t412 * t582 - t635) * t591 + t621 * t594) * t597 + t415, qJD(3) * t25 + qJD(4) * t28 + qJD(6) * t34 + t410, t13 * qJD(3) + t15 * qJD(4) + t34 * qJD(5) + (-mrSges(7,1) * t76 - mrSges(7,2) * t75 - Ifges(7,6) * t500 + t292) * qJD(6) + t411; -qJD(2) * t2 + qJD(4) * t10 - qJD(5) * t44 - qJD(6) * t29 - t656, qJD(4) * t8 + qJD(5) * t27 - qJD(6) * t12 - t416, -qJD(4) * t61 + qJD(5) * t229 + qJD(6) * t102, t137 * qJD(5) + t68 * qJD(6) + ((-t432 * t561 + t333) * t591 + (-pkin(4) * t561 + t333) * t594) * t597 + t383 + (t456 + (-t347 * t636 - t467) * pkin(3)) * qJD(4), qJD(4) * t137 - t380, t68 * qJD(4) + (-t327 * t607 + t417) * qJD(6) + t381; -qJD(2) * t6 - qJD(3) * t10 - qJD(5) * t56 - t656, -qJD(3) * t8 + qJD(5) * t35 - qJD(6) * t14 - t415, -qJD(5) * t136 - qJD(6) * t67 - t383, qJD(5) * t294 + qJD(6) * t110, -t372 ((mrSges(7,2) * t582 - Ifges(7,6)) * t350 + (mrSges(7,1) * t582 - Ifges(7,5)) * t346) * qJD(6) - t379; -qJD(2) * t41 + qJD(3) * t44 + qJD(4) * t56, -qJD(3) * t27 - qJD(4) * t35 - qJD(6) * t33 - t410, qJD(4) * t136 + t380, t372, 0, -qJD(6) * t607 - t495; -t16 * qJD(2) + t29 * qJD(3), qJD(3) * t12 + qJD(4) * t14 + qJD(5) * t33 - t411, qJD(4) * t67 - t381, t379, t495, 0;];
Cq  = t18;
