% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPRRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:15
% EndTime: 2019-03-08 18:56:30
% DurationCPUTime: 8.34s
% Computational Cost: add. (10726->512), mult. (30555->729), div. (0->0), fcn. (33658->12), ass. (0->269)
t363 = cos(qJ(5));
t356 = t363 ^ 2;
t360 = sin(qJ(5));
t632 = t360 ^ 2 + t356;
t361 = sin(qJ(4));
t364 = cos(qJ(4));
t359 = sin(pkin(7));
t362 = sin(qJ(3));
t490 = t359 * t362;
t505 = cos(pkin(7));
t289 = t361 * t505 + t364 * t490;
t543 = cos(qJ(3));
t465 = t359 * t543;
t195 = t360 * t289 + t363 * t465;
t433 = t360 * t465;
t196 = t289 * t363 - t433;
t288 = t361 * t490 - t364 * t505;
t418 = pkin(5) * t363 + qJ(6) * t360;
t293 = t418 * t361;
t486 = t361 * t363;
t318 = -mrSges(6,1) * t364 - mrSges(6,3) * t486;
t319 = mrSges(7,1) * t364 + mrSges(7,2) * t486;
t448 = t318 / 0.2e1 - t319 / 0.2e1;
t424 = t363 * mrSges(6,1) - t360 * mrSges(6,2);
t295 = t424 * t361;
t423 = t363 * mrSges(7,1) + t360 * mrSges(7,3);
t294 = t423 * t361;
t567 = t294 / 0.2e1;
t452 = t295 / 0.2e1 + t567;
t466 = pkin(9) * t360 + pkin(5);
t532 = t361 * pkin(10);
t328 = -pkin(4) * t364 - pkin(3) - t532;
t485 = t363 * t328;
t224 = t364 * t466 - t485;
t487 = t360 * t364;
t253 = -pkin(9) * t487 + t485;
t479 = t224 + t253;
t538 = pkin(9) * t363;
t443 = -qJ(6) + t538;
t489 = t360 * t328;
t223 = t364 * t443 + t489;
t484 = t363 * t364;
t254 = pkin(9) * t484 + t489;
t480 = -t223 + t254;
t579 = m(7) / 0.2e1;
t633 = -t448 * t196 + t452 * t288 + (t195 * t480 + t196 * t479 + t288 * t293) * t579;
t504 = sin(pkin(6));
t442 = sin(pkin(12)) * t504;
t425 = cos(pkin(12)) * t504;
t506 = cos(pkin(6));
t599 = t506 * t359 + t505 * t425;
t188 = t362 * t442 - t543 * t599;
t189 = t362 * t599 + t543 * t442;
t101 = -t188 * t484 + t189 * t360;
t502 = t101 * t363;
t100 = -t188 * t487 - t189 * t363;
t503 = t100 * t360;
t393 = (t502 + t503) * pkin(10);
t326 = -pkin(4) - t418;
t492 = t326 * t361;
t500 = t188 * t361;
t613 = mrSges(7,2) + mrSges(6,3);
t631 = t613 * (-t503 / 0.2e1 - t502 / 0.2e1) - m(6) * (pkin(4) * t500 + t393) / 0.2e1 - m(7) * (-t188 * t492 + t393) / 0.2e1;
t417 = pkin(5) * t360 - qJ(6) * t363;
t407 = pkin(9) + t417;
t270 = t407 * t364;
t522 = mrSges(7,3) * t363;
t523 = mrSges(7,1) * t360;
t422 = -t522 + t523;
t298 = t422 * t364;
t512 = t363 * mrSges(6,2);
t517 = t360 * mrSges(6,1);
t333 = t512 + t517;
t299 = t333 * t364;
t537 = pkin(9) * t364;
t581 = m(6) / 0.2e1;
t630 = (t253 * t360 - t254 * t363 + t537) * t581 + (-t223 * t363 - t224 * t360 + t270) * t579 + t299 / 0.2e1 + t298 / 0.2e1;
t615 = mrSges(6,1) + mrSges(7,1);
t629 = -t615 / 0.2e1;
t628 = Ifges(6,1) + Ifges(7,1);
t525 = Ifges(6,6) - Ifges(7,6);
t627 = t360 * t525;
t528 = Ifges(7,4) + Ifges(6,5);
t626 = t363 * t528;
t625 = t528 * t364;
t377 = -t359 * t425 + t505 * t506;
t133 = t189 * t364 + t361 * t377;
t85 = t133 * t363 + t188 * t360;
t511 = t363 * t85;
t624 = t133 - t511;
t132 = t189 * t361 - t364 * t377;
t535 = pkin(10) * t364;
t340 = t361 * pkin(4) - t535;
t314 = t360 * t340;
t226 = -t361 * t443 + t314;
t491 = t340 * t363;
t231 = -t361 * t466 - t491;
t488 = t360 * t361;
t259 = pkin(9) * t488 + t491;
t260 = -pkin(9) * t486 + t314;
t269 = t407 * t361;
t316 = mrSges(6,2) * t364 - mrSges(6,3) * t488;
t509 = t364 * mrSges(7,3);
t323 = -mrSges(7,2) * t488 - t509;
t450 = t316 / 0.2e1 + t323 / 0.2e1;
t369 = t360 * t448 - t363 * t450 + t630;
t320 = t361 * mrSges(6,1) - mrSges(6,3) * t484;
t475 = mrSges(7,2) * t484;
t513 = t361 * mrSges(7,1);
t321 = t475 - t513;
t447 = -t320 / 0.2e1 + t321 / 0.2e1;
t317 = -t361 * mrSges(6,2) - mrSges(6,3) * t487;
t322 = -mrSges(7,2) * t487 + t361 * mrSges(7,3);
t554 = t322 / 0.2e1;
t449 = t317 / 0.2e1 + t554;
t297 = t333 * t361;
t296 = t422 * t361;
t565 = t296 / 0.2e1;
t451 = t565 + t297 / 0.2e1;
t539 = pkin(9) * t361;
t84 = t133 * t360 - t188 * t363;
t623 = t451 * t133 + t449 * t85 + t447 * t84 + (t133 * t539 - t259 * t84 + t260 * t85) * t581 + (t133 * t269 + t226 * t85 + t231 * t84) * t579 + t369 * t132;
t601 = -t423 - t424;
t622 = t632 * pkin(10);
t510 = t364 * mrSges(5,2);
t334 = t361 * mrSges(5,1) + t510;
t410 = -t334 * t465 / 0.2e1;
t619 = (-t259 * t195 + t260 * t196 + t289 * t539) * t581 + (t195 * t231 + t196 * t226 + t269 * t289) * t579 + t410 + t289 * t451 + t196 * t449 + t195 * t447;
t618 = 0.2e1 * t361;
t617 = m(6) + m(7);
t550 = -t360 / 0.2e1;
t548 = t360 / 0.4e1;
t546 = t363 / 0.2e1;
t616 = t422 / 0.2e1;
t614 = mrSges(6,2) - mrSges(7,3);
t612 = Ifges(7,2) + Ifges(6,3);
t611 = -t364 * mrSges(5,1) + t361 * mrSges(5,2) - mrSges(4,1);
t527 = -Ifges(7,5) + Ifges(6,4);
t605 = t527 * t356;
t445 = t527 * t360;
t604 = -t296 - t297;
t603 = t316 + t323;
t602 = -t318 + t319;
t600 = m(7) * t326 - t423;
t597 = -t632 * t532 / 0.2e1;
t592 = -Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t574 = mrSges(7,3) / 0.2e1;
t575 = -mrSges(6,2) / 0.2e1;
t591 = t575 + t574;
t401 = m(7) * t417;
t469 = -t522 / 0.2e1;
t552 = t333 / 0.2e1;
t588 = t401 + t523 / 0.2e1 + t469 + t552;
t353 = m(7) * qJ(6) + mrSges(7,3);
t587 = t361 * t613;
t586 = -Ifges(5,4) + t626;
t584 = 0.2e1 * m(7);
t355 = t361 ^ 2;
t583 = 2 * qJD(4);
t582 = m(5) / 0.2e1;
t577 = mrSges(6,1) / 0.2e1;
t576 = -mrSges(7,1) / 0.2e1;
t570 = -t195 / 0.2e1;
t568 = t254 / 0.2e1;
t553 = -t423 / 0.2e1;
t551 = t334 / 0.2e1;
t542 = m(7) * t231;
t357 = t364 ^ 2;
t540 = pkin(9) * t357;
t531 = qJD(5) / 0.2e1;
t530 = qJD(6) / 0.2e1;
t526 = -Ifges(6,2) - Ifges(7,3);
t524 = m(7) * qJD(5);
t521 = Ifges(6,4) * t360;
t520 = Ifges(6,4) * t363;
t519 = Ifges(7,5) * t360;
t515 = t360 * t84;
t508 = t364 * t85;
t507 = -t424 - mrSges(5,1);
t501 = t132 * t360;
t499 = t195 * t360;
t498 = t196 * t363;
t251 = -t363 * t490 + t364 * t433;
t496 = t251 * t360;
t463 = t364 * t543;
t252 = (t360 * t362 + t363 * t463) * t359;
t495 = t252 * t363;
t494 = t288 * t360;
t493 = t288 * t363;
t483 = t364 * t196;
t481 = t622 * t132;
t478 = t622 * t288;
t474 = t542 / 0.2e1;
t471 = t577 + mrSges(7,1) / 0.2e1;
t470 = mrSges(6,2) / 0.2e1 - mrSges(7,3) / 0.2e1;
t468 = t512 / 0.2e1;
t464 = t361 * t543;
t455 = -t486 / 0.2e1;
t454 = t486 / 0.4e1;
t444 = -mrSges(7,2) * qJ(6) - Ifges(6,6);
t440 = t526 + t628;
t436 = mrSges(7,2) * pkin(5) - t528;
t435 = t355 * t465;
t434 = t359 * t464;
t432 = t288 * t464;
t427 = qJD(4) * (-t423 + t507);
t12 = m(5) * (-t132 * t361 - t133 * t364 + t189) * t188 + t617 * (t100 * t84 + t85 * t101 - t132 * t500);
t6 = ((-t288 * t361 - t289 * t364) * t188 + (t132 * t464 + t133 * t463 + t188 * t362 - t189 * t543) * t359) * t582 + (t581 + t579) * (t195 * t100 + t196 * t101 + t251 * t84 + t252 * t85 + (t132 * t465 - t188 * t288) * t361);
t416 = t12 * qJD(1) + t6 * qJD(2);
t15 = t617 * (t132 * t624 - t501 * t84);
t7 = 0.2e1 * (t288 * (-t515 + t624) + t132 * (t289 - t498 - t499)) * (m(7) / 0.4e1 + m(6) / 0.4e1);
t415 = t15 * qJD(1) + t7 * qJD(2);
t34 = t617 * (-t195 * t494 - t196 * t493 + t288 * t289);
t414 = t7 * qJD(1) + t34 * qJD(2);
t35 = m(5) * (t289 * t463 - t362 * t465 + t432) * t359 + t617 * (t195 * t251 + t196 * t252 + t359 * t432);
t413 = t6 * qJD(1) + t35 * qJD(2);
t197 = t600 * t360;
t373 = (-t360 * t269 + (-t492 - t535) * t363) * t579 + t296 * t550;
t57 = t475 + (-t423 * t546 + t576) * t361 + t474 - t373;
t409 = qJD(3) * t57 + qJD(4) * t197;
t127 = t509 + (t254 / 0.4e1 - t489 / 0.4e1 + (-t538 / 0.4e1 + qJ(6) / 0.2e1) * t364) * t584;
t408 = qJD(3) * t127 - qJD(5) * t353;
t404 = t517 / 0.2e1 + t468;
t402 = m(7) * (-pkin(5) * t100 + qJ(6) * t101);
t392 = (t495 + t496) * pkin(10);
t390 = (-pkin(5) * t251 + qJ(6) * t252) * t579;
t386 = -t333 / 0.2e1 + t404;
t366 = (-pkin(4) * t434 + t392) * t581 + (t326 * t434 + t392) * t579 + t410 + t601 * t434 / 0.2e1 + t613 * (t496 / 0.2e1 + t495 / 0.2e1);
t11 = -t366 + t448 * t494 - t603 * t493 / 0.2e1 + t630 * t288 + t619;
t348 = Ifges(7,5) * t486;
t17 = pkin(3) * t334 - t223 * t322 - t224 * t321 - t226 * t323 - t231 * t319 - t269 * t298 - t270 * t296 - t253 * t320 - t254 * t317 - t259 * t318 - t260 * t316 - t299 * t539 - t297 * t537 - m(7) * (t223 * t226 + t224 * t231 + t269 * t270) - m(6) * (t253 * t259 + t254 * t260) + (-t586 + t627) * t355 + (t586 * t364 + (-t348 / 0.2e1 - t525 * t364) * t360 + (-Ifges(5,1) + Ifges(5,2) - m(6) * pkin(9) ^ 2 - t628 * t356 + (t526 * t360 + (0.2e1 * Ifges(6,4) - 0.3e1 / 0.2e1 * Ifges(7,5)) * t363) * t360 + t612) * t361) * t364;
t3 = (t551 - t510 / 0.2e1 + (-mrSges(5,1) / 0.2e1 - t424 / 0.2e1 + t553) * t361) * t188 + t623 + t631;
t385 = t3 * qJD(1) + t11 * qJD(2) - t17 * qJD(3);
t18 = -t470 * t252 - t471 * t251 + t450 * t195 + t390 + (t498 / 0.2e1 + t499 / 0.2e1) * t587 - t633;
t347 = Ifges(7,6) * t486;
t20 = t293 * t296 + t269 * t294 + t253 * t323 + m(7) * (t223 * t253 + t224 * t254 + t269 * t293) - t364 * t347 / 0.2e1 + t253 * t316 - t254 * t318 + t254 * t319 + (pkin(9) * t295 + (-t224 * mrSges(7,2) + t253 * mrSges(6,3) + t361 * t445 + t625) * t360 + (-Ifges(6,4) * t486 - t223 * mrSges(7,2) - t254 * mrSges(6,3) + t348 + (Ifges(6,6) - Ifges(7,6) / 0.2e1) * t364 - t440 * t488) * t363) * t361;
t367 = t452 * t132 - t450 * t84 + (t132 * t293 + t480 * t84) * t579 + (-t515 / 0.2e1 - t511 / 0.2e1) * t587 + (t479 * t579 - t448) * t85;
t5 = t470 * t101 + t471 * t100 - t402 / 0.2e1 + t367;
t384 = t5 * qJD(1) - t18 * qJD(2) + t20 * qJD(3);
t42 = (t100 / 0.4e1 + t132 * t454 + t508 / 0.4e1) * t584;
t83 = t296 * t486 + t364 * t323 - m(7) * (-t223 * t364 - t269 * t486);
t87 = (t251 / 0.4e1 + t288 * t454 + t483 / 0.4e1) * t584;
t383 = -qJD(1) * t42 - qJD(2) * t87 - qJD(3) * t83;
t365 = (t417 * t269 + t326 * t293) * t579 + t417 * t565 + t326 * t567 + t539 * t552 + t293 * t553 + (-Ifges(7,3) * t363 + t519) * t486 / 0.2e1 + (Ifges(6,2) * t363 + t521) * t455 + t269 * t616 - pkin(4) * t295 / 0.2e1 - (-Ifges(6,6) * t364 + (-Ifges(6,2) * t360 + t520) * t618) * t360 / 0.4e1 + (-Ifges(7,6) * t364 + 0.2e1 * t348) * t548 - (t626 - t627) * t364 / 0.4e1 + t597 * mrSges(6,3) + (Ifges(7,5) * t546 + Ifges(6,1) * t550 - t520 / 0.2e1 + 0.2e1 * (Ifges(7,3) - Ifges(7,1)) * t548) * t488 + (-t625 + (t363 * t628 + t519 - t521) * t618) * t363 / 0.4e1 + ((t360 * t480 + t363 * t479) * t579 - t448 * t363 + t603 * t550) * pkin(10) + ((t568 - t223 / 0.2e1) * t360 + t479 * t546 + t597) * mrSges(7,2);
t371 = (-pkin(5) * t231 + qJ(6) * t226) * t579 - pkin(5) * t321 / 0.2e1 + qJ(6) * t554 + t226 * t574 + t231 * t576 + t259 * t577 + t260 * t575;
t13 = -t365 + t371 + t612 * t361 / 0.2e1 + t592 * t487 + t528 * t484 / 0.2e1;
t31 = t386 * t288;
t50 = pkin(4) * t333 + t417 * t423 - t605 + (-t363 * t440 + t445) * t360 + (-t401 - t422) * t326;
t9 = t386 * t132;
t379 = -t9 * qJD(1) - t31 * qJD(2) - t13 * qJD(3) - t50 * qJD(4);
t378 = qJD(4) * (-t613 * t632 + mrSges(5,2));
t329 = pkin(9) * t435;
t327 = (m(7) * pkin(10) + mrSges(7,2)) * t363;
t164 = t355 * pkin(9) * t188;
t152 = m(7) * t494;
t115 = (t489 + (-0.2e1 * qJ(6) + t538) * t364) * t579 + m(7) * t568 + t323;
t88 = (-t288 * t486 + t251 - t483) * t579;
t62 = -t423 * t455 + t474 - t513 / 0.2e1 + t373;
t60 = m(7) * t501;
t43 = (-t132 * t486 + t100 - t508) * t579;
t32 = (t616 + t404 + t588) * t288;
t19 = t390 + t603 * t570 + t251 * t629 + t591 * t252 + t613 * (t196 * t455 + t488 * t570) + t633;
t14 = ((Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t363 + t592 * t360) * t364 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t361 + t365 + t371;
t10 = t288 * t369 + t366 + t619;
t8 = t615 * t501 / 0.2e1 + (t468 + t469 + t588) * t132;
t4 = t402 / 0.2e1 + t367 + t100 * t629 + t591 * t101;
t2 = (t510 / 0.2e1 + t551) * t188 + (mrSges(5,1) / 0.2e1 - t601 / 0.2e1) * t500 + t623 - t631;
t1 = qJD(3) * t6 + qJD(4) * t7;
t16 = [qJD(3) * t12 + qJD(4) * t15, t1 (t602 * t100 + t603 * t101 + t611 * t189) * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t43 * qJD(6) + (mrSges(4,2) + t604 * t361 + (-t355 - t357) * mrSges(5,3)) * qJD(3) * t188 + 0.2e1 * ((t100 * t224 + t101 * t223 - t269 * t500) * t579 + (-t100 * t253 + t101 * t254 - t164) * t581 + (-pkin(3) * t189 - t188 * t540 - t164) * t582) * qJD(3) + t416, t2 * qJD(3) + t8 * qJD(5) - t60 * qJD(6) + t133 * t427 + ((t133 * t326 - t481) * t579 + (-pkin(4) * t133 - t481) * t581) * t583 + t132 * t378 + t415, t4 * qJD(3) + t8 * qJD(4) + (t614 * t84 - t615 * t85) * qJD(5) + ((-pkin(5) * t85 - qJ(6) * t84) * t531 + t85 * t530) * t584, t43 * qJD(3) - t60 * qJD(4) + t524 * t85; t1, qJD(3) * t35 + qJD(4) * t34 (-mrSges(4,2) * t465 + m(6) * t329 + m(5) * (t329 + (-pkin(3) * t362 + t540 * t543) * t359) + t611 * t490 + (m(7) * t269 - t604) * t434 + (m(6) * t254 + m(7) * t223 + t603) * t252 + (-m(6) * t253 + m(7) * t224 + t602) * t251 + (t357 * t465 + t435) * mrSges(5,3)) * qJD(3) + t10 * qJD(4) + t19 * qJD(5) + t88 * qJD(6) + t413, t10 * qJD(3) + t32 * qJD(5) - t152 * qJD(6) + t289 * t427 + ((-pkin(4) * t289 - t478) * t581 + (t289 * t326 - t478) * t579) * t583 + t288 * t378 + t414, t19 * qJD(3) + t32 * qJD(4) + (t195 * t614 - t196 * t615) * qJD(5) + ((-pkin(5) * t196 - qJ(6) * t195) * t531 + t196 * t530) * t584, t88 * qJD(3) - t152 * qJD(4) + t196 * t524; qJD(4) * t3 + qJD(5) * t5 - qJD(6) * t42 - t416, qJD(4) * t11 - qJD(5) * t18 - qJD(6) * t87 - t413, -qJD(4) * t17 + qJD(5) * t20 - qJD(6) * t83, t14 * qJD(5) + t62 * qJD(6) + t385 + ((t226 * mrSges(7,2) + t260 * mrSges(6,3) + (m(6) * t260 + m(7) * t226 + t317 + t322) * pkin(10) + t440 * t487) * t363 - pkin(4) * t299 + t326 * t298 + (Ifges(5,5) + t605 + (-m(6) * pkin(4) + t507) * pkin(9)) * t364 + (pkin(9) * mrSges(5,2) + t363 * t525 - Ifges(5,6)) * t361 + (t231 * mrSges(7,2) - t259 * mrSges(6,3) + t528 * t361 + (-m(6) * t259 - t320 + t321 + t542) * pkin(10) - t364 * t445) * t360 + t600 * t270) * qJD(4), t14 * qJD(4) + t115 * qJD(6) + t384 + (t347 + (t360 * t436 + t363 * t444) * t361 + (-m(7) * pkin(5) - t615) * t254 + (-mrSges(6,2) + t353) * t253) * qJD(5), qJD(4) * t62 + qJD(5) * t115 + t383; -qJD(3) * t3 - qJD(5) * t9 - t415, -qJD(3) * t11 - qJD(5) * t31 - t414, -qJD(5) * t13 - qJD(6) * t57 - t385, -qJD(5) * t50 - qJD(6) * t197, t327 * qJD(6) + t379 + (-t436 * t363 + (Ifges(7,6) + t444) * t360 + (-m(7) * t418 + t601) * pkin(10)) * qJD(5), qJD(5) * t327 - t409; -qJD(3) * t5 + qJD(4) * t9, qJD(3) * t18 + qJD(4) * t31, qJD(4) * t13 - qJD(6) * t127 - t384, -t379, t353 * qJD(6), -t408; qJD(3) * t42, qJD(3) * t87, qJD(4) * t57 + qJD(5) * t127 - t383, t409, t408, 0;];
Cq  = t16;
