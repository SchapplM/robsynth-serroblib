% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:58
% EndTime: 2019-03-09 08:46:17
% DurationCPUTime: 10.00s
% Computational Cost: add. (20256->506), mult. (38360->665), div. (0->0), fcn. (44010->8), ass. (0->296)
t484 = sin(pkin(10));
t485 = cos(pkin(10));
t525 = sin(qJ(2));
t527 = cos(qJ(2));
t298 = t484 * t525 - t485 * t527;
t300 = -t484 * t527 - t485 * t525;
t328 = -t527 * pkin(2) - pkin(1);
t234 = t298 * pkin(3) + t300 * qJ(4) + t328;
t176 = -pkin(4) * t298 - t234;
t336 = sin(qJ(5));
t526 = cos(qJ(5));
t241 = t298 * t336 - t526 * t300;
t373 = -t298 * t526 - t336 * t300;
t621 = -m(6) * t176 - mrSges(6,1) * t373 - mrSges(6,2) * t241;
t337 = cos(qJ(6));
t528 = t337 / 0.2e1;
t445 = t525 * pkin(7);
t371 = -qJ(3) * t525 - t445;
t447 = t527 * pkin(7);
t372 = qJ(3) * t527 + t447;
t253 = -t485 * t371 + t484 * t372;
t353 = t300 * pkin(8) + t253;
t579 = t484 * t371 + t485 * t372;
t580 = t298 * pkin(8) + t579;
t601 = t336 * t580 - t353 * t526;
t335 = sin(qJ(6));
t608 = t335 / 0.2e1;
t516 = Ifges(7,4) * t337;
t314 = -Ifges(7,2) * t335 + t516;
t92 = Ifges(7,6) * t241 - t314 * t373;
t517 = Ifges(7,4) * t335;
t316 = Ifges(7,1) * t337 - t517;
t95 = Ifges(7,5) * t241 - t316 * t373;
t619 = t601 * mrSges(6,2) + t92 * t528 + t608 * t95;
t463 = t336 * t601;
t118 = t336 * t353 + t526 * t580;
t612 = t526 * t118;
t618 = t612 + t463;
t513 = Ifges(6,6) * t241;
t494 = t337 * mrSges(7,1);
t502 = t335 * mrSges(7,2);
t390 = t494 - t502;
t614 = t118 * t390;
t615 = t118 * mrSges(6,1);
t617 = -t513 + t619 - t614 - t615;
t616 = -t614 / 0.2e1 - t615 / 0.2e1;
t613 = t118 * t601;
t493 = t337 * mrSges(7,2);
t503 = t335 * mrSges(7,1);
t310 = t493 + t503;
t568 = t310 * t373;
t610 = -m(7) * t118 + t568;
t154 = t310 * t241;
t609 = t118 * t154 - t601 * t568;
t422 = t485 * pkin(2);
t327 = -t422 - pkin(3);
t319 = -pkin(4) + t327;
t421 = t484 * pkin(2);
t322 = t421 + qJ(4);
t264 = t336 * t319 + t322 * t526;
t605 = t264 * t601;
t604 = t335 * t601;
t603 = t337 * t601;
t333 = t335 ^ 2;
t334 = t337 ^ 2;
t450 = t333 + t334;
t597 = m(5) * t234 + mrSges(5,1) * t298 + mrSges(5,3) * t300;
t596 = -t300 * mrSges(4,1) - t298 * mrSges(4,2);
t531 = -t335 / 0.2e1;
t330 = Ifges(7,6) * t335;
t514 = Ifges(7,5) * t337;
t562 = -t330 + t514;
t587 = -t241 / 0.2e1;
t573 = Ifges(6,4) * t587;
t574 = t373 / 0.2e1;
t586 = t241 / 0.2e1;
t595 = -t95 * t528 - t92 * t531 - t562 * t586 - t573 - (Ifges(6,2) + Ifges(7,3)) * t574;
t582 = t390 * t241;
t594 = t582 / 0.2e1;
t446 = t525 * pkin(2);
t593 = m(4) * t446;
t591 = Ifges(7,3) * t586;
t403 = t450 * t373;
t523 = pkin(5) * t241;
t166 = pkin(9) * t373 + t523;
t588 = t601 / 0.2e1;
t585 = mrSges(7,2) * t241;
t584 = t241 * mrSges(7,1);
t506 = t241 * mrSges(6,3);
t306 = t336 * t390;
t581 = -t568 * t526 / 0.2e1;
t370 = t450 * t526;
t237 = -t300 * pkin(3) + t298 * qJ(4) + t446;
t177 = -t300 * pkin(4) + t237;
t77 = -t177 - t166;
t56 = t335 * t77 + t603;
t491 = t337 * t56;
t55 = t337 * t77 - t604;
t499 = t335 * t55;
t387 = t491 - t499;
t578 = -m(5) * t322 - mrSges(5,3);
t577 = t594 + mrSges(7,3) * t403 / 0.2e1;
t575 = -t373 / 0.2e1;
t570 = -Ifges(5,5) + Ifges(4,4);
t391 = mrSges(7,3) * (t334 / 0.2e1 + t333 / 0.2e1);
t315 = Ifges(7,1) * t335 + t516;
t454 = t337 * t315;
t313 = Ifges(7,2) * t337 + t517;
t464 = t335 * t313;
t374 = t464 / 0.2e1 - t454 / 0.2e1;
t569 = -Ifges(6,5) + t374;
t378 = -t514 / 0.2e1 + t330 / 0.2e1;
t565 = t378 * t373;
t367 = t336 * mrSges(6,1) + mrSges(6,2) * t526 - mrSges(7,3) * t370;
t563 = -t306 - t367;
t501 = t335 * mrSges(7,3);
t434 = t501 / 0.2e1;
t76 = pkin(5) * t373 - pkin(9) * t241 + t176;
t54 = t118 * t337 + t335 * t76;
t500 = t335 * t54;
t561 = t54 * t434 - mrSges(7,3) * t500 / 0.2e1;
t430 = t337 * t526;
t392 = -t430 / 0.2e1;
t431 = t335 * t526;
t393 = -t431 / 0.2e1;
t359 = mrSges(7,1) * t393 + mrSges(7,2) * t392;
t425 = t526 * t310;
t350 = t425 / 0.2e1 + t359;
t559 = qJD(6) * t350;
t558 = t350 * qJD(4);
t351 = -t425 / 0.2e1 + t359;
t557 = t351 * qJD(6);
t515 = Ifges(7,5) * t373;
t97 = t241 * t316 + t515;
t487 = t337 * t97;
t512 = Ifges(7,6) * t373;
t94 = t241 * t314 + t512;
t497 = t335 * t94;
t556 = t176 * mrSges(6,2) + Ifges(6,1) * t586 - Ifges(6,4) * t373 - t497 / 0.2e1 + t487 / 0.2e1;
t555 = -(t315 / 0.4e1 + t314 / 0.4e1) * t335 + (t316 / 0.4e1 - t313 / 0.4e1) * t337;
t406 = t316 / 0.2e1 - t313 / 0.2e1;
t408 = t315 / 0.2e1 + t314 / 0.2e1;
t554 = -t406 * t335 - t408 * t337;
t155 = t241 * t313;
t156 = t241 * t315;
t553 = -t337 * t155 / 0.4e1 - t335 * t156 / 0.4e1 + t310 * t588 + t487 / 0.4e1 - t497 / 0.4e1;
t552 = 0.2e1 * m(7);
t551 = m(5) / 0.2e1;
t550 = m(6) / 0.2e1;
t549 = -m(7) / 0.2e1;
t548 = m(7) / 0.2e1;
t547 = -pkin(5) / 0.2e1;
t546 = m(4) * pkin(2);
t545 = -t55 / 0.2e1;
t544 = t56 / 0.2e1;
t68 = t166 * t337 + t604;
t543 = -t68 / 0.2e1;
t69 = t335 * t166 - t603;
t542 = t69 / 0.2e1;
t263 = t319 * t526 - t336 * t322;
t261 = pkin(5) - t263;
t535 = -t261 / 0.2e1;
t262 = -pkin(9) + t264;
t534 = -t262 / 0.2e1;
t533 = -t263 / 0.2e1;
t530 = -t336 / 0.2e1;
t529 = -t337 / 0.2e1;
t524 = m(7) * (-t526 + t370) * t336;
t522 = pkin(5) * t310;
t465 = t335 * t373;
t158 = -mrSges(7,3) * t465 + t585;
t443 = t241 * t501;
t159 = -mrSges(7,2) * t373 - t443;
t455 = t337 * t373;
t161 = -mrSges(7,3) * t455 - t584;
t492 = t337 * mrSges(7,3);
t162 = mrSges(7,1) * t373 - t241 * t492;
t235 = t241 * mrSges(6,1);
t277 = t300 * mrSges(5,1);
t440 = Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1;
t53 = -t335 * t118 + t337 * t76;
t1 = (mrSges(4,1) * t446 + t234 * mrSges(5,3) + (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t300 + t570 * t298) * t298 + (-mrSges(4,2) * t446 - t570 * t300) * t300 - t234 * t277 - pkin(1) * (mrSges(3,1) * t525 + mrSges(3,2) * t527) + m(7) * (t53 * t55 + t54 * t56 - t613) - t176 * t235 + t54 * t158 + t56 * t159 + t53 * t161 + t55 * t162 + (Ifges(3,1) - Ifges(3,2)) * t527 * t525 + (t562 * t574 + t556) * t373 + (-t525 ^ 2 + t527 ^ 2) * Ifges(3,4) + (Ifges(6,1) * t574 + Ifges(6,4) * t586 - t440 * t373 + t595) * t241 + (t593 + t596) * t328 + t597 * t237 - t609 + t621 * t177;
t511 = t1 * qJD(1);
t505 = t373 * mrSges(6,2);
t504 = t298 * mrSges(5,2);
t498 = t335 * t68;
t490 = t337 * t69;
t5 = t54 * t162 - t601 * t582 + ((t97 / 0.2e1 - t155 / 0.2e1 + t515 / 0.2e1) * t335 + (t156 / 0.2e1 + t94 / 0.2e1 + t54 * mrSges(7,3) + t512 / 0.2e1) * t337) * t241 + (-t159 - t443) * t53;
t486 = t5 * qJD(1);
t388 = -t335 * t53 + t337 * t54;
t457 = t337 * t159;
t466 = t335 * t162;
t479 = t601 * t241;
t14 = -(t154 + t506) * t241 + (-mrSges(6,3) * t373 + t457 - t466) * t373 + m(7) * (t373 * t388 - t479) + m(6) * (t118 * t373 - t479) + (t298 ^ 2 + t300 ^ 2) * (mrSges(4,3) + mrSges(5,2)) + (m(5) + m(4)) * (-t253 * t300 - t298 * t579);
t483 = qJD(1) * t14;
t456 = t337 * t162;
t468 = t335 * t159;
t19 = (-t468 - t456 + m(7) * (-t337 * t53 - t500) + t597 + t621) * t300;
t482 = qJD(1) * t19;
t340 = (-t298 * t322 - t300 * t327) * t551 + (t241 * t263 + t264 * t373) * t550 + (-t241 * t261 + t262 * t403) * t548 + (-t298 * t484 + t300 * t485) * t546 / 0.2e1 - t577;
t343 = t237 * t551 + t177 * t550 + (-t335 * t56 - t337 * t55) * t548 + t158 * t531 + t161 * t529 + t593 / 0.2e1;
t15 = -t298 * mrSges(5,3) - t235 + t277 + t340 - t343 + t505 - t596;
t476 = t15 * qJD(1);
t157 = t373 * t501 - t585;
t160 = t373 * t492 + t584;
t347 = (-t335 * t69 - t337 * t68) * t549 + mrSges(6,1) * t586 + t157 * t608 + t160 * t528;
t349 = t235 / 0.2e1 + (pkin(9) * t403 + t523) * t548 + t577;
t20 = 0.2e1 * mrSges(6,2) * t575 + t347 + t349;
t475 = t20 * qJD(1);
t473 = t373 * t562;
t311 = Ifges(7,5) * t335 + Ifges(7,6) * t337;
t472 = t241 * t311;
t377 = -t493 / 0.2e1 - t503 / 0.2e1;
t368 = t377 * t373;
t420 = t466 / 0.2e1;
t375 = t420 - t457 / 0.2e1;
t25 = t368 - t375;
t471 = t25 * qJD(1);
t470 = t261 * t310;
t467 = t335 * t160;
t462 = t336 * t154;
t461 = t336 * t241;
t459 = t337 * t157;
t458 = t337 * t158;
t426 = t526 * t241;
t344 = (t336 * t373 + t426) * t550 + (t336 * t403 + t426) * t548;
t366 = (m(7) * t450 + m(6)) * t300 / 0.2e1;
t40 = -m(5) * t300 + t344 - t366;
t453 = t40 * qJD(1);
t449 = mrSges(7,3) * t490;
t444 = qJD(4) * t524;
t442 = t568 * t547;
t441 = -pkin(9) * t161 / 0.2e1;
t439 = mrSges(6,3) * t526;
t432 = t264 * t526;
t427 = t526 * t582;
t417 = t335 * t530;
t416 = -t461 / 0.2e1;
t415 = t459 / 0.2e1;
t414 = t458 / 0.2e1;
t413 = t457 / 0.2e1;
t412 = -t456 / 0.2e1;
t402 = t450 * t263;
t396 = mrSges(6,3) * t416;
t395 = t162 * t393 + t526 * t413 + t462 / 0.2e1 + mrSges(6,3) * t461 / 0.2e1;
t386 = t490 - t498;
t356 = -t568 * t535 - t264 * t154 / 0.2e1 + t616;
t381 = t261 * t118 + t605;
t2 = t381 * t549 + t442 + (t586 + t587) * Ifges(6,6) + (t575 + t574) * Ifges(6,5) + (t588 - t601 / 0.2e1) * mrSges(6,2) - (-mrSges(6,1) / 0.2e1 + m(7) * t547 - t390 / 0.2e1) * t118 + (pkin(9) * t158 / 0.2e1 + t159 * t533 + t157 * t534 + (t542 + t544) * mrSges(7,3) + (pkin(9) * t56 / 0.4e1 - t262 * t69 / 0.4e1 - t263 * t54 / 0.4e1) * t552) * t337 + (t441 + t263 * t162 / 0.2e1 + t262 * t160 / 0.2e1 + (t543 + t545) * mrSges(7,3) + (-pkin(9) * t55 / 0.4e1 + t262 * t68 / 0.4e1 + t263 * t53 / 0.4e1) * t552) * t335 + t356;
t357 = mrSges(7,3) * t402 - t263 * mrSges(6,2) + (-mrSges(6,1) - t390) * t264;
t36 = -m(7) * (t261 * t264 + t262 * t402) + t357;
t385 = -t2 * qJD(1) - t36 * qJD(2);
t362 = t430 * t54 - t431 * t53;
t13 = t396 + (-t612 + (t601 + t386) * t336 + t362) * t548 - t581 + t336 * t415 + t160 * t417 + t395;
t4 = t69 * t159 + t54 * t157 + t68 * t162 + t53 * t160 + m(7) * (t53 * t68 + t54 * t69 + t613) + (t176 * mrSges(6,1) + t573 - t595) * t241 + (t565 + (-Ifges(6,1) / 0.2e1 + t440) * t241 - t556) * t373 + t609;
t384 = t4 * qJD(1) + t13 * qJD(4);
t338 = t579 * t551 + t618 * t550 + (t336 * t387 + t612) * t548 + t581 + t161 * t417 - t506 * t530 + t336 * t414 + t439 * t575;
t341 = -m(5) * t579 / 0.2e1 - m(6) * t618 / 0.2e1 + (t362 + t463) * t549 + t439 * t574;
t10 = t338 + t526 * t420 + t396 + t159 * t392 - t462 / 0.2e1 + t341;
t363 = t370 * t262;
t62 = -m(7) * (t261 * t336 + t363) - m(6) * (-t263 * t336 + t432) + t563 + t578;
t383 = t10 * qJD(1) + t62 * qJD(2);
t120 = t470 + t554;
t376 = -t468 / 0.2e1 + t412;
t339 = t376 * t262 + (-t262 * t391 - t555) * t241 - t473 / 0.4e1 + t261 * t594 - t553 + t561;
t364 = mrSges(7,1) * t545 + mrSges(7,2) * t544 + t591;
t7 = t339 + t364 + t565;
t382 = t7 * qJD(1) - t120 * qJD(2);
t369 = (t502 / 0.2e1 - t494 / 0.2e1) * t300;
t23 = t427 / 0.2e1 + t369 + (t241 * t391 - t376) * t336;
t380 = -t23 * qJD(1) - qJD(2) * t351;
t379 = mrSges(7,1) * t543 + mrSges(7,2) * t542;
t342 = t306 / 0.2e1 + (-t432 + t363 + (t261 + t402) * t336) * t548;
t355 = m(7) * (-pkin(5) * t336 + pkin(9) * t370);
t346 = t355 / 0.2e1 - t306 / 0.2e1;
t42 = t342 - t346 + t367;
t365 = t13 * qJD(1) + t42 * qJD(2) + t444;
t167 = -t522 - t554;
t57 = (t535 + t547) * t310 + (mrSges(7,2) * t533 + t408) * t337 + (mrSges(7,1) * t533 + t406) * t335;
t345 = t376 * pkin(9) + t547 * t582 + t553 + t561;
t348 = -pkin(9) * t391 + t555;
t9 = (t562 / 0.4e1 - t378) * t373 + (-Ifges(7,3) / 0.2e1 + t348) * t241 + t345 + t379;
t358 = t9 * qJD(1) - t57 * qJD(2) + t167 * qJD(5) - t558;
t58 = t316 * t531 + t314 * t529 + t470 / 0.2e1 + t522 / 0.2e1 + t377 * t263 + t374;
t43 = t342 + t346;
t41 = t344 + t366;
t26 = t368 + t375;
t24 = -t427 / 0.2e1 + t159 * t417 + t336 * t412 + t369 + t450 * mrSges(7,3) * t416;
t21 = mrSges(6,2) * t574 - t505 / 0.2e1 - t347 + t349;
t16 = t340 + t343;
t12 = t13 * qJD(5);
t11 = t338 - t341 + t395 - t504;
t8 = t345 + t348 * t241 + t591 + t473 / 0.4e1 - t379 + t565;
t6 = t339 - Ifges(7,6) * t465 / 0.2e1 + Ifges(7,5) * t455 / 0.2e1 - t364;
t3 = (pkin(5) * t118 + pkin(9) * t387 + t262 * t386 + t263 * t388 + t381) * t548 + t442 + t513 / 0.2e1 - t449 / 0.2e1 - t472 / 0.4e1 + 0.2e1 * (-t464 / 0.4e1 + t454 / 0.4e1) * t373 + 0.2e1 * t574 * Ifges(6,5) - (-Ifges(6,6) / 0.2e1 + t311 / 0.4e1) * t241 + t335 * t441 + pkin(9) * t414 + t262 * t415 + t68 * t434 + (-t499 / 0.2e1 + t491 / 0.2e1) * mrSges(7,3) + t466 * t533 + t467 * t534 - t356 + t263 * t413 - t616 - t619;
t17 = [qJD(2) * t1 + qJD(3) * t14 + qJD(4) * t19 + qJD(5) * t4 - qJD(6) * t5, t16 * qJD(3) + t11 * qJD(4) + t3 * qJD(5) + t6 * qJD(6) + t511 + ((mrSges(5,2) * t322 + mrSges(4,3) * t421 + Ifges(4,6) - Ifges(5,6)) * t300 + (mrSges(4,3) * t422 - Ifges(5,4) - Ifges(4,5)) * t298 + (m(7) * t387 - t335 * t161 + t458) * t262 + t610 * t261 - t387 * mrSges(7,3) + m(6) * (t118 * t263 + t605) + Ifges(3,5) * t527 - Ifges(3,6) * t525 + t264 * t506 - t327 * t504 - mrSges(3,1) * t447 + mrSges(3,2) * t445 + (-t263 * mrSges(6,3) + t569) * t373 + (m(5) * t327 - t485 * t546 - mrSges(4,1) - mrSges(5,1)) * t579 - (t484 * t546 - mrSges(4,2) - t578) * t253 + t311 * t586 + t617) * qJD(2), qJD(2) * t16 + qJD(4) * t41 + qJD(5) * t21 + qJD(6) * t26 + t483, qJD(2) * t11 + qJD(3) * t41 + qJD(6) * t24 + t12 + t482, t3 * qJD(2) + t21 * qJD(3) + t8 * qJD(6) + t384 + (t449 - mrSges(7,3) * t498 + t472 / 0.2e1 + t569 * t373 + t610 * pkin(5) + (m(7) * t386 + t459 - t467) * pkin(9) + t617) * qJD(5), -t486 + t6 * qJD(2) + t26 * qJD(3) + t24 * qJD(4) + t8 * qJD(5) + (-t54 * mrSges(7,1) - t53 * mrSges(7,2) - t472) * qJD(6); qJD(3) * t15 - qJD(4) * t10 - qJD(5) * t2 + qJD(6) * t7 - t511, -qJD(4) * t62 - qJD(5) * t36 - qJD(6) * t120, t476, t43 * qJD(5) - t383 + t444 + t559, t43 * qJD(4) + (m(7) * (-pkin(5) * t264 + pkin(9) * t402) + t357) * qJD(5) + t58 * qJD(6) + t385, t558 + t58 * qJD(5) + (-t262 * t390 - t562) * qJD(6) + t382; -qJD(2) * t15 - qJD(4) * t40 - qJD(5) * t20 - qJD(6) * t25 - t483, -t476, 0, -t453, -t475, qJD(6) * t310 - t471; qJD(2) * t10 + qJD(3) * t40 - qJD(6) * t23 + t12 - t482, qJD(5) * t42 + t383 - t557, t453, qJD(5) * t524 (t355 + t563) * qJD(5) + t557 + t365, qJD(5) * t351 - qJD(6) * t306 + t380; qJD(2) * t2 + qJD(3) * t20 + qJD(6) * t9 - t384, -qJD(4) * t42 - qJD(6) * t57 - t385, t475, -t365 - t559, t167 * qJD(6) (-pkin(9) * t390 + t562) * qJD(6) + t358; -qJD(2) * t7 + qJD(3) * t25 + qJD(4) * t23 - qJD(5) * t9 + t486, qJD(4) * t351 + qJD(5) * t57 - t382, t471, qJD(5) * t350 - t380, -t358, 0;];
Cq  = t17;
