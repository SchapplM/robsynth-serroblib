% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:21
% EndTime: 2019-03-09 03:51:43
% DurationCPUTime: 12.94s
% Computational Cost: add. (41287->631), mult. (82982->872), div. (0->0), fcn. (99101->10), ass. (0->332)
t371 = sin(pkin(11));
t373 = cos(pkin(11));
t376 = sin(qJ(5));
t379 = cos(qJ(5));
t358 = t371 * t379 + t376 * t373;
t375 = sin(qJ(6));
t378 = cos(qJ(6));
t416 = t376 * t371 - t373 * t379;
t318 = t358 * t378 - t375 * t416;
t443 = -t358 * t375 - t378 * t416;
t662 = t318 * mrSges(7,1) + t443 * mrSges(7,2);
t659 = qJD(6) * t662;
t466 = Ifges(7,5) * t443 - Ifges(7,6) * t318;
t546 = pkin(8) + qJ(4);
t362 = t546 * t371;
t363 = t546 * t373;
t331 = -t376 * t362 + t363 * t379;
t277 = -pkin(9) * t416 + t331;
t329 = -t379 * t362 - t363 * t376;
t414 = -pkin(9) * t358 + t329;
t183 = t277 * t378 + t375 * t414;
t608 = -t277 * t375 + t378 * t414;
t655 = -t183 * mrSges(7,1) - t608 * mrSges(7,2);
t32 = t466 + t655;
t663 = t32 * qJD(6);
t661 = t608 / 0.2e1;
t372 = sin(pkin(10));
t374 = cos(pkin(10));
t377 = sin(qJ(3));
t557 = cos(qJ(3));
t356 = t372 * t377 - t374 * t557;
t359 = t372 * t557 + t377 * t374;
t283 = t416 * t359;
t403 = t358 * t359;
t418 = t378 * t283 + t375 * t403;
t526 = t418 * Ifges(7,4);
t643 = t283 * t375 - t378 * t403;
t102 = Ifges(7,2) * t643 + t356 * Ifges(7,6) - t526;
t161 = -mrSges(7,2) * t356 + mrSges(7,3) * t643;
t547 = pkin(7) + qJ(2);
t364 = t547 * t374;
t449 = t547 * t372;
t330 = t364 * t377 + t557 * t449;
t480 = t359 * t371;
t263 = pkin(4) * t480 + t330;
t218 = pkin(5) * t403 + t263;
t366 = -pkin(4) * t373 - pkin(3);
t334 = pkin(5) * t416 + t366;
t618 = t643 * mrSges(7,2);
t636 = t418 * mrSges(7,1);
t644 = -t636 + t618;
t647 = Ifges(7,1) * t643 + t526;
t658 = (t647 / 0.4e1 - t102 / 0.4e1) * t318 + t161 * t661 + t218 * t662 / 0.2e1 + t334 * t644 / 0.2e1;
t600 = m(7) * pkin(5);
t463 = t600 / 0.2e1;
t657 = t218 * t644;
t284 = t356 * t358;
t287 = t416 * t356;
t214 = t284 * t375 + t287 * t378;
t528 = t214 * mrSges(7,2);
t211 = t284 * t378 - t287 * t375;
t531 = t211 * mrSges(7,1);
t473 = t531 / 0.2e1 - t528 / 0.2e1;
t523 = t287 * mrSges(6,2);
t525 = t284 * mrSges(6,1);
t656 = t525 / 0.2e1 - t523 / 0.2e1 + t473;
t516 = t318 * Ifges(7,4);
t237 = Ifges(7,2) * t443 + t516;
t648 = Ifges(7,1) * t443 - t516;
t654 = t237 / 0.4e1 - t648 / 0.4e1;
t640 = -mrSges(7,2) / 0.2e1;
t641 = mrSges(7,1) / 0.2e1;
t474 = t418 * t641 + t640 * t643;
t598 = -mrSges(6,2) / 0.2e1;
t599 = mrSges(6,1) / 0.2e1;
t650 = -t403 * t598 + (t599 + (t375 ^ 2 + t378 ^ 2) * t463) * t283 + t474;
t649 = (t211 * t378 + t214 * t375) * t463 + t656;
t620 = Ifges(7,5) * t643;
t638 = Ifges(7,6) * t418;
t472 = t620 + t638;
t205 = Ifges(7,4) * t643;
t104 = -Ifges(7,1) * t418 + t356 * Ifges(7,5) + t205;
t629 = Ifges(7,2) * t418 + t205;
t646 = t629 + t104;
t409 = t618 / 0.2e1 - t636 / 0.2e1;
t452 = -t638 / 0.2e1 - t620 / 0.2e1;
t163 = mrSges(7,1) * t356 + mrSges(7,3) * t418;
t573 = -t318 / 0.2e1;
t574 = t443 / 0.2e1;
t642 = t161 * t574 + t163 * t573;
t566 = -t358 / 0.2e1;
t562 = t359 / 0.2e1;
t569 = -t416 / 0.2e1;
t453 = -pkin(2) * t374 - pkin(1);
t308 = pkin(3) * t356 - qJ(4) * t359 + t453;
t332 = t364 * t557 - t377 * t449;
t240 = t373 * t308 - t332 * t371;
t241 = t371 * t308 + t373 * t332;
t419 = t240 * t371 - t241 * t373;
t479 = t359 * t373;
t159 = pkin(4) * t356 - pkin(8) * t479 + t240;
t200 = -pkin(8) * t480 + t241;
t89 = t376 * t159 + t200 * t379;
t71 = -pkin(9) * t403 + t89;
t508 = t375 * t71;
t88 = t379 * t159 - t200 * t376;
t70 = pkin(9) * t283 + t88;
t67 = pkin(5) * t356 + t70;
t51 = t378 * t67 - t508;
t507 = t378 * t71;
t52 = t375 * t67 + t507;
t425 = -t318 * t51 + t443 * t52;
t248 = -mrSges(6,2) * t356 - mrSges(6,3) * t403;
t250 = mrSges(6,1) * t356 + mrSges(6,3) * t283;
t441 = t248 * t569 + t250 * t566 + t642;
t320 = -t356 * mrSges(5,2) - mrSges(5,3) * t480;
t477 = t373 * t320;
t322 = t356 * mrSges(5,1) - mrSges(5,3) * t479;
t478 = t371 * t322;
t602 = -m(7) / 0.2e1;
t604 = -m(6) / 0.2e1;
t605 = m(5) / 0.2e1;
t631 = t419 * t605 + (t283 * t329 - t331 * t403 - t358 * t88 - t416 * t89) * t604 + (t183 * t643 + t418 * t608 + t425) * t602 + t478 / 0.2e1 - t477 / 0.2e1 - t441;
t502 = qJ(4) * t356;
t556 = pkin(3) * t359;
t326 = t502 + t556;
t244 = t373 * t326 + t330 * t371;
t245 = t371 * t326 - t373 * t330;
t630 = -t244 * t371 + t245 * t373;
t587 = t214 / 0.2e1;
t590 = t211 / 0.2e1;
t412 = Ifges(7,5) * t587 + Ifges(7,6) * t590;
t482 = t356 * t373;
t168 = pkin(4) * t359 + pkin(8) * t482 + t244;
t483 = t356 * t371;
t220 = pkin(8) * t483 + t245;
t105 = t379 * t168 - t220 * t376;
t69 = pkin(5) * t359 - pkin(9) * t287 + t105;
t106 = t376 * t168 + t379 * t220;
t72 = pkin(9) * t284 + t106;
t59 = -t375 * t72 + t378 * t69;
t60 = t375 * t69 + t378 * t72;
t439 = Ifges(7,3) * t562 + t59 * t641 + t60 * t640 + t412;
t312 = Ifges(7,4) * t443;
t628 = -Ifges(7,2) * t318 + t312;
t627 = -t102 / 0.2e1;
t624 = t643 / 0.2e1;
t615 = t378 * t643;
t571 = t416 / 0.2e1;
t614 = t403 * t571;
t575 = -t443 / 0.2e1;
t613 = t643 * t575;
t280 = Ifges(6,4) * t403;
t172 = -Ifges(6,1) * t283 + t356 * Ifges(6,5) - t280;
t612 = Ifges(6,2) * t283 + t172 - t280;
t367 = t371 ^ 2;
t369 = t373 ^ 2;
t465 = t369 + t367;
t611 = -Ifges(6,5) * t416 - Ifges(6,6) * t358 + t466;
t351 = Ifges(6,4) * t416;
t328 = t358 * Ifges(6,1) - t351;
t610 = -Ifges(6,2) * t358 + t328 - t351;
t609 = -Ifges(6,5) * t403 + Ifges(6,6) * t283 + t472;
t606 = 0.2e1 * t359;
t603 = m(6) / 0.2e1;
t601 = m(7) / 0.2e1;
t596 = -mrSges(7,3) / 0.2e1;
t595 = mrSges(7,3) / 0.2e1;
t593 = t183 / 0.2e1;
t592 = -t608 / 0.2e1;
t591 = -t183 / 0.2e1;
t586 = -t418 / 0.2e1;
t584 = t237 / 0.2e1;
t239 = t318 * Ifges(7,1) + t312;
t583 = t239 / 0.2e1;
t582 = t284 / 0.2e1;
t581 = t403 / 0.4e1;
t580 = -t403 / 0.2e1;
t579 = t287 / 0.2e1;
t578 = t283 / 0.2e1;
t577 = -t283 / 0.2e1;
t576 = -t283 / 0.4e1;
t572 = t318 / 0.2e1;
t570 = t416 / 0.4e1;
t567 = t356 / 0.2e1;
t565 = t358 / 0.2e1;
t564 = t358 / 0.4e1;
t561 = t371 / 0.2e1;
t560 = t373 / 0.2e1;
t559 = -t375 / 0.2e1;
t558 = t378 / 0.2e1;
t555 = pkin(5) * t283;
t554 = t358 * pkin(5);
t553 = t51 * mrSges(7,2);
t552 = t52 * mrSges(7,1);
t55 = -t375 * t70 - t507;
t551 = t55 * mrSges(7,1);
t56 = t378 * t70 - t508;
t550 = t56 * mrSges(7,2);
t545 = mrSges(7,3) * t318;
t544 = Ifges(5,4) * t371;
t543 = Ifges(5,4) * t373;
t542 = Ifges(6,4) * t283;
t541 = Ifges(6,4) * t358;
t537 = pkin(5) * qJD(5);
t524 = t403 * mrSges(6,2);
t522 = t283 * mrSges(6,1);
t101 = Ifges(7,4) * t214 + Ifges(7,2) * t211 + Ifges(7,6) * t359;
t103 = Ifges(7,1) * t214 + Ifges(7,4) * t211 + Ifges(7,5) * t359;
t114 = t528 - t531;
t115 = -mrSges(7,1) * t643 - mrSges(7,2) * t418;
t160 = -mrSges(7,2) * t359 + t211 * mrSges(7,3);
t162 = mrSges(7,1) * t359 - t214 * mrSges(7,3);
t169 = Ifges(6,4) * t287 + Ifges(6,2) * t284 + Ifges(6,6) * t359;
t170 = -Ifges(6,2) * t403 + t356 * Ifges(6,6) - t542;
t171 = Ifges(6,1) * t287 + Ifges(6,4) * t284 + Ifges(6,5) * t359;
t264 = -pkin(4) * t483 + t332;
t219 = -t284 * pkin(5) + t264;
t221 = t523 - t525;
t222 = mrSges(6,1) * t403 - mrSges(6,2) * t283;
t247 = -t359 * mrSges(6,2) + t284 * mrSges(6,3);
t249 = t359 * mrSges(6,1) - t287 * mrSges(6,3);
t512 = t371 * Ifges(5,2);
t258 = Ifges(5,6) * t359 + (t512 - t543) * t356;
t259 = Ifges(5,5) * t359 + (-Ifges(5,1) * t373 + t544) * t356;
t510 = t373 * mrSges(5,2);
t513 = t371 * mrSges(5,1);
t437 = t510 + t513;
t306 = t437 * t356;
t307 = t437 * t359;
t319 = -t359 * mrSges(5,2) + mrSges(5,3) * t483;
t321 = t359 * mrSges(5,1) + mrSges(5,3) * t482;
t350 = t359 * mrSges(4,1);
t413 = Ifges(6,5) * t579 + Ifges(6,6) * t582;
t509 = t373 * Ifges(5,5);
t511 = t371 * Ifges(5,6);
t3 = t102 * t590 + (-t453 * mrSges(4,2) + (-Ifges(5,1) * t369 / 0.2e1 + Ifges(5,3) + Ifges(6,3) + Ifges(7,3) - Ifges(4,1) + Ifges(4,2) + (t543 - t512 / 0.2e1) * t371) * t359 + t412 + t413 + (Ifges(4,4) - t509 + t511) * t356) * t356 + t453 * t350 + m(5) * (t240 * t244 + t241 * t245 + t330 * t332) + m(6) * (t105 * t88 + t106 * t89 + t263 * t264) + m(7) * (t218 * t219 + t51 * t59 + t52 * t60) + (-t371 * t258 / 0.2e1 + t259 * t560 + Ifges(7,5) * t586 + Ifges(7,6) * t624 + Ifges(6,5) * t577 + Ifges(6,6) * t580 + (t509 / 0.2e1 - t511 / 0.2e1 - Ifges(4,4)) * t359) * t359 + t101 * t624 + t332 * t307 - t330 * t306 + t241 * t319 + t245 * t320 + t240 * t321 + t244 * t322 + t264 * t222 + t263 * t221 + t89 * t247 + t106 * t248 + t88 * t249 + t105 * t250 + t218 * t114 + t219 * t115 + t59 * t163 + t52 * t160 + t60 * t161 + t51 * t162 + t171 * t577 + t172 * t579 + t169 * t580 + t170 * t582 + t103 * t586 + t104 * t587;
t521 = t3 * qJD(1);
t519 = t443 * mrSges(7,3);
t515 = t416 * mrSges(6,3);
t514 = t358 * mrSges(6,3);
t433 = -Ifges(6,1) * t403 + t542;
t436 = -t522 - t524;
t4 = m(7) * (-t218 * t555 + t51 * t55 + t52 * t56) + t56 * t161 + t55 * t163 - t115 * t555 + t657 + t647 * t586 - t418 * t627 + t263 * t436 + t433 * t577 + t88 * t248 - t89 * t250 + t170 * t578 + (t418 * t52 - t51 * t643) * mrSges(7,3) + (t283 * t89 + t403 * t88) * mrSges(6,3) + t612 * t580 + t609 * t567 + t646 * t624;
t506 = t4 * qJD(1);
t505 = t51 * t443;
t504 = t52 * t318;
t7 = t51 * t161 - t52 * t163 + t657 + t472 * t567 - (-t52 * mrSges(7,3) + t647 / 0.2e1 + t627) * t418 + (-t51 * mrSges(7,3) + t104 / 0.2e1 + t629 / 0.2e1) * t643;
t503 = t7 * qJD(1);
t235 = -mrSges(7,1) * t443 + mrSges(7,2) * t318;
t435 = mrSges(6,1) * t416 + t358 * mrSges(6,2);
t438 = -t373 * mrSges(5,1) + t371 * mrSges(5,2);
t381 = -m(5) * (-t465 * t502 - t556) / 0.2e1 + (t284 * t329 + t287 * t331 + t359 * t366) * t604 + (t183 * t214 + t211 * t608 + t334 * t359) * t602 - (t235 + t438 + t435) * t359 / 0.2e1;
t384 = (t244 * t373 + t245 * t371) * t605 + (-t105 * t416 + t106 * t358) * t603 + (t318 * t60 + t443 * t59) * t601 + t162 * t574 + t160 * t572 + t249 * t569 + t247 * t565 + t319 * t561 + t321 * t560;
t14 = (-mrSges(4,2) + (t369 / 0.2e1 + t367 / 0.2e1) * mrSges(5,3)) * t356 + (t211 * t572 + t214 * t575) * mrSges(7,3) + (t284 * t565 + t287 * t571) * mrSges(6,3) + t350 + t381 + t384;
t501 = qJD(1) * t14;
t19 = t643 * t161 + t418 * t163 - t403 * t248 + t283 * t250 + m(7) * (t418 * t51 + t52 * t643) + m(6) * (t283 * t88 - t403 * t89) + (-t373 * t322 - t371 * t320 + m(5) * (-t240 * t373 - t241 * t371)) * t359;
t500 = qJD(1) * t19;
t486 = t330 * t359;
t16 = t214 * t161 + t211 * t163 + t287 * t248 + t284 * t250 + (mrSges(4,3) * t356 - t477 + t478) * t356 + (mrSges(4,3) * t359 + t115 + t222 + t307) * t359 + m(7) * (t211 * t51 + t214 * t52 + t218 * t359) + m(6) * (t263 * t359 + t284 * t88 + t287 * t89) + m(5) * (t356 * t419 + t486) + m(4) * (-t332 * t356 + t486) + (m(3) * qJ(2) + mrSges(3,3)) * (t372 ^ 2 + t374 ^ 2);
t499 = t16 * qJD(1);
t391 = (-t418 * t573 + t613) * mrSges(7,3) + t642;
t17 = t391 - t473;
t498 = t17 * qJD(1);
t497 = t183 * t418;
t492 = t283 * t358;
t490 = t318 * t378;
t488 = t443 * t375;
t487 = t329 * t403;
t485 = t331 * t283;
t461 = -t555 / 0.2e1;
t460 = t554 / 0.2e1;
t459 = -t545 / 0.2e1;
t458 = mrSges(7,3) * t592;
t457 = t519 / 0.2e1;
t456 = -t515 / 0.2e1;
t454 = -t514 / 0.2e1;
t450 = m(5) * t465;
t448 = t358 * mrSges(6,1) - mrSges(6,2) * t416;
t442 = m(7) / 0.4e1 + m(6) / 0.4e1 + m(5) / 0.4e1;
t432 = -Ifges(6,1) * t416 - t541;
t26 = t334 * t662 + (t648 / 0.2e1 - t237 / 0.2e1) * t318 + (t583 + t628 / 0.2e1) * t443;
t383 = (t104 / 0.4e1 + t629 / 0.4e1) * t443 + (t458 + t239 / 0.4e1 + t628 / 0.4e1) * t643 - (mrSges(7,3) * t591 - t654) * t418 + t163 * t591 + t356 * t466 / 0.4e1 + t658;
t6 = t383 - t439;
t424 = t6 * qJD(1) + t26 * qJD(3);
t397 = -t524 / 0.2e1 + t409;
t35 = m(7) * t461 - t522 / 0.2e1 + t397 - t650;
t402 = -t662 - t448;
t92 = (t565 + t490 / 0.2e1 - t488 / 0.2e1) * t600 - t402;
t423 = qJD(1) * t35 + qJD(3) * t92;
t387 = (-t283 * t416 - t358 * t403) * t604 + (t318 * t643 + t418 * t443) * t602;
t61 = (t450 / 0.4e1 + t442) * t606 + t387;
t422 = qJD(1) * t61;
t420 = t218 * t358 - t283 * t334;
t65 = -t474 + t409;
t415 = qJD(1) * t65 + qJD(3) * t662;
t405 = t418 * t318;
t404 = t643 * t443;
t401 = (t488 - t490) * t600;
t327 = -Ifges(6,2) * t416 + t541;
t380 = t163 * t593 - t263 * t448 / 0.2e1 + t327 * t576 - t329 * t248 / 0.2e1 + t331 * t250 / 0.2e1 + t170 * t564 - t366 * t436 / 0.2e1 + t610 * t581 + t612 * t570 - t611 * t356 / 0.4e1 - t646 * t443 / 0.4e1 - (t628 + t239) * t643 / 0.4e1 - t654 * t418 - t658;
t388 = Ifges(6,3) * t562 + t105 * t599 + t106 * t598 + t413 + t439;
t393 = t162 * t558 + (t375 * t60 + t378 * t59) * t601 + t375 * t160 / 0.2e1;
t396 = (t52 + t55) * t608 + (-t51 + t56) * t183;
t2 = (-t283 * t570 + t358 * t581) * Ifges(6,1) + (t505 / 0.2e1 + t55 * t572 + t504 / 0.2e1 + t56 * t575 + t608 * t624 - t497 / 0.2e1) * mrSges(7,3) + (-t487 / 0.2e1 - t485 / 0.2e1) * mrSges(6,3) + t380 + (t115 * t566 + t235 * t578 + t420 * t602 + t393) * pkin(5) + t396 * t602 + t388 - Ifges(6,4) * t492 / 0.2e1;
t24 = t235 * t554 - t318 * t584 + t327 * t566 + t366 * t448 + t648 * t572 + t432 * t565 + t443 * t583 + t610 * t569 + t574 * t628 + (m(7) * t554 + t662) * t334;
t400 = -t2 * qJD(1) + t24 * qJD(3);
t385 = (t418 * t572 + t613) * mrSges(7,3) + (t492 / 0.2e1 - t614) * mrSges(6,3) + (t318 * t56 + t443 * t55 + t425) * t601 + t441;
t8 = t385 - t649;
t399 = t8 * qJD(1);
t386 = (-t510 / 0.2e1 - t513 / 0.2e1) * t356 + t332 * t605 + t264 * t603 + t219 * t601 - t656;
t12 = (-t404 / 0.2e1 + t405 / 0.2e1) * mrSges(7,3) + (t283 * t565 - t614) * mrSges(6,3) + t386 + t631;
t39 = (t318 ^ 2 + t443 ^ 2) * mrSges(7,3) + (t358 ^ 2 + t416 ^ 2) * mrSges(6,3) + m(7) * (t183 * t443 - t318 * t608) + m(6) * (-t329 * t358 - t331 * t416) + (m(5) * qJ(4) + mrSges(5,3)) * t465;
t398 = qJD(1) * t12 - qJD(3) * t39;
t389 = (t161 * t558 + t163 * t559 + (-t418 * t559 - t615 / 0.2e1) * mrSges(7,3)) * pkin(5) - t452;
t11 = (-t51 / 0.2e1 + t56 / 0.2e1) * mrSges(7,2) + (-t52 / 0.2e1 - t55 / 0.2e1) * mrSges(7,1) + t389 + t452;
t33 = (t592 + t661) * mrSges(7,2) + (t591 + t593) * mrSges(7,1);
t361 = (t375 * mrSges(7,1) + t378 * mrSges(7,2)) * pkin(5);
t395 = -qJD(1) * t11 - qJD(3) * t33 + qJD(5) * t361;
t360 = t361 * qJD(6);
t173 = m(7) * t460 + t401 / 0.2e1;
t66 = t409 + t474;
t62 = t442 * t606 - t450 * t562 - t387;
t34 = -(t463 + t599) * t283 + t397 + t650;
t18 = t391 + t473;
t15 = t211 * t459 + t214 * t457 + t284 * t454 + t287 * t456 - t381 + t384 - t465 * t356 * mrSges(5,3) / 0.2e1;
t13 = t283 * t454 - t403 * t456 + t404 * t595 + t405 * t596 + t386 - t631;
t10 = -t553 / 0.2e1 - t552 / 0.2e1 - t550 / 0.2e1 + t551 / 0.2e1 + t389 - t452;
t9 = t385 + t649;
t5 = t383 + t439;
t1 = -t380 + t393 * pkin(5) + (pkin(5) * t420 + t396) * t601 + t433 * t564 + t432 * t576 + t115 * t460 + t388 + t505 * t596 + t55 * t459 + t56 * t457 + t235 * t461 + t643 * t458 - (-t487 - t485) * mrSges(6,3) / 0.2e1 + (-t504 + t497) * t595;
t20 = [qJD(2) * t16 + qJD(3) * t3 + qJD(4) * t19 + qJD(5) * t4 + qJD(6) * t7, t499 + 0.2e1 * ((t211 * t443 + t214 * t318) * t601 + (-t284 * t416 + t287 * t358) * t603) * qJD(2) + t15 * qJD(3) + t62 * qJD(4) + t9 * qJD(5) + t18 * qJD(6), t15 * qJD(2) + t13 * qJD(4) + t1 * qJD(5) + t5 * qJD(6) + t521 + (0.2e1 * (t105 * t329 + t106 * t331 + t264 * t366) * t603 + (-t373 * (Ifges(5,1) * t371 + t543) / 0.2e1 + (Ifges(5,2) * t373 + t544) * t561 - Ifges(4,5)) * t356 + t171 * t565 + t169 * t569 + t103 * t572 + t101 * t574 + t258 * t560 + t259 * t561 - t59 * t545 - t105 * t514 - t106 * t515 + t60 * t519 + t264 * t435 + t332 * t438 - t371 * qJ(4) * t321 + t373 * qJ(4) * t319 + (Ifges(5,5) * t371 + Ifges(6,5) * t358 + Ifges(7,5) * t318 + Ifges(5,6) * t373 - Ifges(6,6) * t416 + Ifges(7,6) * t443) * t562 + 0.2e1 * (t183 * t60 + t219 * t334 + t59 * t608) * t601 + t608 * t162 + 0.2e1 * (-pkin(3) * t332 + qJ(4) * t630) * t605 + t630 * mrSges(5,3) + t366 * t221 - Ifges(4,6) * t359 - t332 * mrSges(4,1) + t334 * t114 + t329 * t249 + t330 * mrSges(4,2) + t331 * t247 + pkin(3) * t306 + t219 * t235 + t183 * t160 + t328 * t579 + t327 * t582 + t214 * t583 + t211 * t584) * qJD(3), qJD(2) * t62 + qJD(3) * t13 + qJD(5) * t34 + qJD(6) * t66 + t500, t506 + t9 * qJD(2) + t1 * qJD(3) + t34 * qJD(4) + (-t89 * mrSges(6,1) - t88 * mrSges(6,2) - t550 + t551 + t609) * qJD(5) + t10 * qJD(6) + (m(7) * (t375 * t56 + t378 * t55) + (t375 * t418 - t615) * mrSges(7,3)) * t537, t503 + t18 * qJD(2) + t5 * qJD(3) + t66 * qJD(4) + t10 * qJD(5) + (t472 - t552 - t553) * qJD(6); qJD(3) * t14 - qJD(4) * t61 + qJD(5) * t8 + qJD(6) * t17 - t499, 0, t501, -t422 (t401 + t402) * qJD(5) - t659 + t399, -qJD(5) * t662 + t498 - t659; -qJD(2) * t14 - qJD(4) * t12 - qJD(5) * t2 + qJD(6) * t6 - t521, -t501, qJD(4) * t39 + qJD(5) * t24 + qJD(6) * t26, qJD(5) * t173 - t398, t173 * qJD(4) + (-t331 * mrSges(6,1) - t329 * mrSges(6,2) + t611 + t655) * qJD(5) + t663 + (m(7) * (-t183 * t378 + t375 * t608) + (-t318 * t375 - t378 * t443) * mrSges(7,3)) * t537 + t400, t32 * qJD(5) + t424 + t663; qJD(2) * t61 + qJD(3) * t12 + qJD(5) * t35 + qJD(6) * t65 - t500, t422, qJD(5) * t92 + t398 + t659, 0, t423, t415; -qJD(2) * t8 + qJD(3) * t2 - qJD(4) * t35 + qJD(6) * t11 - t506, -t399, -qJD(4) * t92 + qJD(6) * t33 - t400, -t423, -t360, -t360 - t395; -qJD(2) * t17 - qJD(3) * t6 - qJD(4) * t65 - qJD(5) * t11 - t503, -t498, -qJD(4) * t662 - t33 * qJD(5) - t424, -t415, t395, 0;];
Cq  = t20;
