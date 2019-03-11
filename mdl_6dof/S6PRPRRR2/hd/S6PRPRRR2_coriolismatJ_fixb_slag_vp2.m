% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:59
% EndTime: 2019-03-08 20:27:15
% DurationCPUTime: 9.20s
% Computational Cost: add. (17363->602), mult. (41980->854), div. (0->0), fcn. (45580->12), ass. (0->325)
t376 = sin(qJ(4));
t561 = t376 / 0.2e1;
t373 = sin(pkin(12));
t505 = sin(pkin(6));
t506 = cos(pkin(12));
t425 = t506 * t505;
t377 = sin(qJ(2));
t434 = t377 * t505;
t556 = cos(qJ(2));
t291 = t373 * t434 - t425 * t556;
t427 = t556 * t505;
t292 = t373 * t427 + t377 * t425;
t379 = cos(qJ(5));
t375 = sin(qJ(5));
t380 = cos(qJ(4));
t476 = t375 * t380;
t176 = t291 * t476 + t292 * t379;
t468 = t379 * t380;
t177 = -t291 * t468 + t292 * t375;
t374 = sin(qJ(6));
t378 = cos(qJ(6));
t103 = t176 * t378 - t177 * t374;
t104 = t176 * t374 + t177 * t378;
t587 = mrSges(6,2) / 0.2e1;
t588 = -mrSges(6,1) / 0.2e1;
t589 = m(7) * pkin(5);
t605 = t103 * mrSges(7,1) / 0.2e1 - t104 * mrSges(7,2) / 0.2e1;
t622 = t176 * t588 + t177 * t587 - (t103 * t378 + t104 * t374) * t589 / 0.2e1 - t605;
t552 = t376 * pkin(4);
t553 = pkin(9) * t380;
t349 = t552 - t553;
t359 = pkin(2) * t373 + pkin(8);
t477 = t375 * t376;
t257 = t379 * t349 + t359 * t477;
t222 = t376 * pkin(5) - pkin(10) * t468 + t257;
t473 = t376 * t379;
t258 = t375 * t349 - t359 * t473;
t234 = -pkin(10) * t476 + t258;
t127 = t222 * t378 - t234 * t374;
t128 = t222 * t374 + t234 * t378;
t415 = t374 * t379 + t378 * t375;
t301 = t415 * t380;
t470 = t378 * t379;
t414 = t374 * t375 - t470;
t302 = t414 * t380;
t412 = Ifges(7,5) * t302 / 0.2e1 + Ifges(7,6) * t301 / 0.2e1;
t620 = mrSges(7,2) / 0.2e1;
t621 = -mrSges(7,1) / 0.2e1;
t600 = t127 * t621 + t128 * t620 + t412;
t595 = Ifges(7,3) * t561 - t600;
t507 = cos(pkin(6));
t243 = t292 * t380 + t376 * t507;
t157 = -t243 * t375 + t291 * t379;
t158 = t243 * t379 + t291 * t375;
t431 = t378 * t157 - t158 * t374;
t96 = t157 * t374 + t158 * t378;
t30 = -t96 * mrSges(7,1) - t431 * mrSges(7,2);
t619 = t30 * qJD(6);
t585 = -pkin(10) - pkin(9);
t347 = t585 * t375;
t348 = t585 * t379;
t256 = t347 * t374 - t348 * t378;
t430 = t378 * t347 + t348 * t374;
t319 = Ifges(7,6) * t415;
t320 = Ifges(7,5) * t414;
t464 = -t320 - t319;
t61 = -t256 * mrSges(7,1) - t430 * mrSges(7,2) + t464;
t618 = t61 * qJD(6);
t360 = -pkin(2) * t506 - pkin(3);
t318 = -t380 * pkin(4) - t376 * pkin(9) + t360;
t303 = t379 * t318;
t457 = pkin(10) * t473;
t207 = -t457 + t303 + (-t359 * t375 - pkin(5)) * t380;
t237 = t375 * t318 + t359 * t468;
t219 = -pkin(10) * t477 + t237;
t495 = t219 * t374;
t115 = t207 * t378 - t495;
t236 = -t359 * t476 + t303;
t218 = t236 - t457;
t121 = t218 * t378 - t495;
t617 = t115 - t121;
t299 = t374 * t477 - t376 * t470;
t300 = t415 * t376;
t209 = -t299 * mrSges(7,1) - t300 * mrSges(7,2);
t616 = t209 * qJD(6);
t494 = t219 * t378;
t116 = t207 * t374 + t494;
t120 = -t218 * t374 - t494;
t613 = t116 + t120;
t496 = t177 * t379;
t497 = t176 * t375;
t416 = t496 - t497;
t363 = -pkin(5) * t379 - pkin(4);
t483 = t363 * t376;
t487 = t291 * t376;
t610 = -(t496 / 0.2e1 - t497 / 0.2e1) * mrSges(6,3) - m(6) * (pkin(4) * t487 + pkin(9) * t416) / 0.2e1 - m(7) * (t103 * t430 + t104 * t256 - t291 * t483) / 0.2e1;
t539 = Ifges(7,4) * t415;
t609 = t256 * mrSges(7,3) - t539;
t514 = t415 * mrSges(7,3);
t366 = Ifges(6,5) * t379;
t451 = -t366 / 0.2e1;
t536 = Ifges(6,6) * t375;
t606 = Ifges(5,4) + t451 + t536 / 0.2e1;
t601 = -mrSges(6,1) * t379 + mrSges(6,2) * t375;
t313 = t601 * t376;
t367 = Ifges(6,4) * t379;
t604 = -Ifges(6,2) * t375 + t367;
t343 = Ifges(6,1) * t375 + t367;
t602 = -t257 * t375 + t258 * t379;
t516 = t302 * mrSges(7,2);
t517 = t301 * mrSges(7,1);
t466 = -t517 / 0.2e1 + t516 / 0.2e1;
t242 = t292 * t376 - t380 * t507;
t133 = t415 * t242;
t134 = t414 * t242;
t411 = t133 * t621 + t134 * t620;
t369 = t375 ^ 2;
t371 = t379 ^ 2;
t463 = t369 + t371;
t598 = mrSges(6,3) * t463 - mrSges(5,2);
t541 = Ifges(7,4) * t299;
t199 = -Ifges(7,2) * t300 - t380 * Ifges(7,6) - t541;
t279 = Ifges(7,4) * t300;
t201 = -Ifges(7,1) * t299 - t380 * Ifges(7,5) - t279;
t212 = t299 * Ifges(7,2) - t279;
t213 = -t300 * Ifges(7,1) + t541;
t554 = pkin(5) * t375;
t435 = t359 + t554;
t307 = t435 * t376;
t454 = Ifges(7,2) / 0.4e1 - Ifges(7,1) / 0.4e1;
t565 = t363 / 0.2e1;
t518 = t300 * mrSges(7,3);
t261 = mrSges(7,2) * t380 - t518;
t580 = t261 / 0.2e1;
t244 = mrSges(7,1) * t415 - mrSges(7,2) * t414;
t582 = t244 / 0.2e1;
t597 = t430 * t580 - (t201 / 0.4e1 - t279 / 0.2e1 + t212 / 0.4e1 + t454 * t299) * t414 - (-t213 / 0.4e1 - t541 / 0.2e1 + t199 / 0.4e1 - t454 * t300) * t415 + t307 * t582 + t209 * t565;
t575 = t300 / 0.2e1;
t576 = t299 / 0.2e1;
t596 = t431 * t580 + (t431 * t575 + t576 * t96) * mrSges(7,3);
t594 = 2 * qJD(4);
t593 = m(6) / 0.2e1;
t592 = m(6) / 0.4e1;
t591 = m(7) / 0.2e1;
t586 = t96 / 0.2e1;
t584 = t157 / 0.2e1;
t583 = t242 / 0.2e1;
t245 = mrSges(7,1) * t414 + mrSges(7,2) * t415;
t581 = t245 / 0.2e1;
t263 = -mrSges(7,1) * t380 + t299 * mrSges(7,3);
t579 = -t263 / 0.2e1;
t578 = t263 / 0.2e1;
t577 = -t299 / 0.2e1;
t574 = -t300 / 0.2e1;
t573 = -t301 / 0.2e1;
t572 = -t302 / 0.2e1;
t571 = -t313 / 0.2e1;
t317 = t376 * t343;
t570 = -t317 / 0.4e1;
t569 = t415 / 0.2e1;
t335 = -mrSges(6,1) * t380 - mrSges(6,3) * t473;
t568 = -t335 / 0.2e1;
t512 = t379 * mrSges(6,2);
t513 = t375 * mrSges(6,1);
t340 = t512 + t513;
t567 = t340 / 0.2e1;
t511 = t380 * mrSges(5,2);
t341 = t376 * mrSges(5,1) + t511;
t566 = t341 / 0.2e1;
t564 = -t375 / 0.2e1;
t562 = t375 / 0.2e1;
t560 = -t379 / 0.2e1;
t559 = t379 / 0.2e1;
t558 = -t380 / 0.2e1;
t557 = -t380 / 0.4e1;
t547 = Ifges(7,1) - Ifges(7,2);
t543 = Ifges(6,1) * t379;
t542 = Ifges(6,4) * t375;
t540 = Ifges(7,4) * t414;
t533 = pkin(5) * qJD(5);
t530 = t115 * mrSges(7,2);
t529 = t116 * mrSges(7,1);
t528 = t120 * mrSges(7,1);
t527 = t121 * mrSges(7,2);
t510 = t380 * Ifges(6,5);
t509 = t380 * Ifges(6,6);
t508 = t601 - mrSges(5,1);
t504 = mrSges(7,3) * qJD(4);
t486 = t291 * t380;
t49 = (-t300 * t103 - t299 * t104) * t591 + 0.2e1 * (m(7) * t486 / 0.4e1 + (t416 + t486) * t592) * t376;
t503 = qJD(1) * t49;
t502 = t104 * t414;
t501 = t127 * t378;
t500 = t128 * t374;
t499 = t157 * t375;
t498 = t158 * t379;
t493 = t242 * t375;
t492 = t243 * t380;
t370 = t376 ^ 2;
t489 = t291 * t370;
t372 = t380 ^ 2;
t488 = t291 * t372;
t485 = t307 * t375;
t484 = t359 * t376;
t37 = t261 * t574 + t263 * t576 + t209 * t558 + (-t299 ^ 2 / 0.2e1 - t300 ^ 2 / 0.2e1) * mrSges(7,3);
t482 = t37 * qJD(2);
t262 = -mrSges(7,2) * t376 - mrSges(7,3) * t301;
t481 = t374 * t262;
t480 = t374 * t299;
t295 = t376 * t604 - t509;
t479 = t375 * t295;
t336 = t376 * mrSges(6,1) - mrSges(6,3) * t468;
t478 = t375 * t336;
t475 = t376 * t242;
t474 = t376 * t307;
t264 = mrSges(7,1) * t376 + mrSges(7,3) * t302;
t472 = t378 * t264;
t471 = t378 * t300;
t334 = -t376 * mrSges(6,2) - mrSges(6,3) * t476;
t469 = t379 * t334;
t467 = -t363 * t244 + t256 * t514;
t465 = -Ifges(7,5) * t300 + Ifges(7,6) * t299;
t462 = qJD(4) * t380;
t461 = t589 / 0.2e1;
t308 = t435 * t380;
t211 = -t516 + t517;
t315 = t380 * t340;
t333 = mrSges(6,2) * t380 - mrSges(6,3) * t477;
t399 = t315 / 0.2e1 + t211 / 0.2e1 + t333 * t560 + t335 * t562;
t210 = mrSges(7,1) * t300 - mrSges(7,2) * t299;
t314 = t340 * t376;
t436 = t314 / 0.2e1 + t210 / 0.2e1;
t21 = t261 * t572 + t262 * t577 + t263 * t573 + t264 * t574 - t399 * t380 + (t469 / 0.2e1 - t478 / 0.2e1 + t436) * t376 + (-t301 * t115 - t302 * t116 - t300 * t127 - t299 * t128 - t308 * t380 + t474) * t591 + ((t237 * t380 + t258 * t376) * t379 + (-t236 * t380 - t257 * t376) * t375 + (t370 - t372) * t359) * t593;
t429 = mrSges(6,3) * (-t371 / 0.2e1 - t369 / 0.2e1);
t439 = t333 * t564;
t458 = pkin(5) * t473;
t28 = (t299 * t617 - t613 * t300 - t380 * t458) * t591 + t376 * t439 + t473 * t568 - t313 * t558 + t370 * t429 + t37;
t460 = t21 * qJD(4) + t28 * qJD(5) + t37 * qJD(6);
t459 = -0.1e1 + t463;
t452 = t586 - t96 / 0.2e1;
t446 = -t514 / 0.2e1;
t445 = t510 / 0.2e1;
t444 = t244 * t558 + t299 * t446 - t514 * t577;
t443 = t209 * t583;
t442 = t242 * t582;
t440 = t210 * t562;
t438 = -t115 / 0.2e1 + t121 / 0.2e1;
t437 = t116 / 0.2e1 + t120 / 0.2e1;
t432 = t366 - t536;
t428 = qJD(4) * (t245 + t508);
t426 = t464 * t557;
t423 = -t542 + t543;
t342 = Ifges(6,2) * t379 + t542;
t18 = t442 - t411;
t421 = t28 * qJD(2);
t168 = t291 * t475;
t19 = m(7) * (t103 * t431 + t104 * t96 - t168) + m(6) * (t157 * t176 + t158 * t177 - t168) + m(5) * (t292 - t475 - t492) * t291;
t420 = t19 * qJD(1) + t49 * qJD(3);
t417 = -t498 + t499;
t22 = (-t300 * t133 - t299 * t134 - t301 * t431 - t302 * t96 + t475 - t492) * t591 + ((-t243 - t417) * t380 - t459 * t475) * t593;
t153 = t242 * t243;
t25 = m(7) * (t133 * t431 + t134 * t96 + t153) + m(6) * (t242 * t417 + t153);
t419 = t25 * qJD(1) + t22 * qJD(3);
t410 = t257 * t588 + t258 * t587;
t409 = -t512 / 0.2e1 - t513 / 0.2e1;
t408 = t256 * t576 + t430 * t575;
t407 = t342 * t564 + t343 * t559;
t406 = (t133 * t378 + t134 * t374) * t589;
t200 = -Ifges(7,4) * t302 - Ifges(7,2) * t301 + Ifges(7,6) * t376;
t202 = -Ifges(7,1) * t302 - Ifges(7,4) * t301 + Ifges(7,5) * t376;
t296 = Ifges(6,6) * t376 + t380 * t604;
t297 = t376 * t423 - t510;
t298 = Ifges(6,5) * t376 + t380 * t423;
t10 = t360 * t341 + t258 * t333 + t237 * t334 + t257 * t335 + t236 * t336 + t307 * t211 + t308 * t210 + t200 * t574 + t199 * t573 + t201 * t572 + t202 * t577 + t128 * t261 + t116 * t262 + t127 * t263 + t115 * t264 + m(6) * (t236 * t257 + t237 * t258) + m(7) * (t115 * t127 + t116 * t128 + t307 * t308) + (t297 * t559 - t479 / 0.2e1 + t359 * t314 + t412 + t606 * t380) * t380 + (Ifges(7,5) * t577 + Ifges(7,6) * t574 + t298 * t559 + t296 * t564 + t359 * t315 - t606 * t376 + (m(6) * t359 ^ 2 + Ifges(5,1) - Ifges(5,2) - Ifges(6,3) - Ifges(7,3)) * t380) * t376;
t381 = t436 * t243 + t399 * t242 + (t243 * t484 + t257 * t157 + t258 * t158 + (t236 * t375 - t237 * t379 + t359 * t380) * t242) * t593 + (t115 * t133 + t116 * t134 + t127 * t431 + t128 * t96 + t242 * t308 + t243 * t307) * t591 + t133 * t578 + t134 * t580 + t336 * t584 + t158 * t334 / 0.2e1 + t431 * t264 / 0.2e1 + t262 * t586;
t2 = t381 + (t566 - t511 / 0.2e1 + (-mrSges(5,1) / 0.2e1 + t601 / 0.2e1 + t581) * t376) * t291 + (t502 / 0.2e1 + t103 * t569) * mrSges(7,3) + t610;
t404 = t2 * qJD(1) + t10 * qJD(2) + t21 * qJD(3);
t316 = t376 * t342;
t389 = t307 * t209 - (t201 / 0.2e1 + t212 / 0.2e1) * t300 + (-t213 / 0.2e1 + t199 / 0.2e1 + t116 * mrSges(7,3)) * t299 + t115 * t518 + t465 * t558;
t13 = t121 * t261 + m(7) * (t115 * t120 + t116 * t121) + t120 * t263 + t236 * t333 - t237 * t335 + (-t359 * t313 + (t236 * mrSges(6,3) + t445 - t297 / 0.2e1 + t316 / 0.2e1) * t375 + (-t237 * mrSges(6,3) + t509 / 0.2e1 - t317 / 0.2e1 - t295 / 0.2e1 + (m(7) * t307 + t210) * pkin(5)) * t379) * t376 + t389;
t384 = (-t498 / 0.2e1 + t499 / 0.2e1) * t376 * mrSges(6,3) + (t242 * t458 + t613 * t431) * t591 + t333 * t584 + t158 * t568 + t596 + (-t591 * t617 - t578) * t96;
t4 = (t571 + t209 / 0.2e1) * t242 + t384 + t622;
t403 = t4 * qJD(1) + t13 * qJD(2) + t28 * qJD(3);
t390 = t579 * t96 + t443 + t596;
t11 = t390 - t605;
t16 = t115 * t261 - t116 * t263 + t389;
t402 = t11 * qJD(1) + t16 * qJD(2) + t37 * qJD(3);
t98 = m(7) * (t299 * t302 + t300 * t301) + 0.4e1 * (t459 * t592 - m(7) / 0.4e1) * t380 * t376;
t401 = t22 * qJD(1) + t21 * qJD(2) + t98 * qJD(3);
t400 = t567 + t582 + t409;
t39 = pkin(4) * t340 + t604 * t560 + t423 * t564 - t609 * t415 - (-t415 * t547 + t540) * t414 - t407 + t467 + (-m(7) * t363 - t245) * t554;
t382 = (-t316 / 0.4e1 + t297 / 0.4e1 + pkin(9) * t568) * t379 + (-t414 * t438 - t415 * t437 + t408) * mrSges(7,3) + pkin(4) * t313 / 0.2e1 - t256 * t578 + t597;
t385 = t359 * t567 + pkin(9) * t429 + (-t542 / 0.4e1 - t342 / 0.4e1 + t543 / 0.4e1 + (m(7) * t565 + t581) * pkin(5)) * t379 - (t343 + t604) * t375 / 0.4e1;
t395 = -t256 * t617 + t613 * t430;
t5 = t382 + t395 * t591 + (-t472 / 0.2e1 + t440 - t481 / 0.2e1 + 0.2e1 * (t485 / 0.4e1 - t501 / 0.4e1 - t500 / 0.4e1) * m(7)) * pkin(5) + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t385) * t376 + (t570 - t295 / 0.4e1 - pkin(9) * t333 / 0.2e1) * t375 + (0.3e1 / 0.4e1 * t536 - t366 / 0.4e1 + t320 / 0.4e1 + t319 / 0.4e1 + t451) * t380 + t410 + t600;
t394 = t476 * t589;
t397 = (-t301 * t378 - t302 * t374) * t461 + t466;
t50 = t400 * t380 + t394 / 0.2e1 + t397;
t387 = pkin(5) * t493 * t591 - t452 * t514;
t9 = t400 * t242 - t406 / 0.2e1 + t387 + t411;
t398 = t9 * qJD(1) + t5 * qJD(2) - t50 * qJD(3) - t39 * qJD(4);
t383 = mrSges(7,3) * t408 + t256 * t579 + t426 + t597;
t15 = t383 - t595;
t17 = t442 + t411;
t52 = -Ifges(7,4) * t414 ^ 2 - (-t414 * t547 + t609) * t415 + t467;
t55 = t444 - t466;
t396 = t17 * qJD(1) + t15 * qJD(2) + t55 * qJD(3) - t52 * qJD(4);
t29 = t452 * mrSges(7,1);
t338 = (mrSges(7,1) * t374 + mrSges(7,2) * t378) * pkin(5);
t391 = (t378 * t580 + t374 * t579 + (t480 / 0.2e1 + t471 / 0.2e1) * mrSges(7,3)) * pkin(5);
t34 = -mrSges(7,1) * t437 + mrSges(7,2) * t438 + t391;
t392 = t29 * qJD(1) - t34 * qJD(2) + t338 * qJD(5);
t327 = t338 * qJD(6);
t250 = t359 * t489;
t56 = t444 + t466;
t51 = t340 * t558 - t394 / 0.2e1 + t409 * t380 + t397 + t444;
t33 = -t530 / 0.2e1 - t529 / 0.2e1 - t527 / 0.2e1 + t528 / 0.2e1 + t391 + t465;
t14 = t383 + t595;
t12 = t390 + t605;
t8 = t242 * t567 + t406 / 0.2e1 + t512 * t583 + mrSges(6,1) * t493 / 0.2e1 + t387 + t18;
t7 = t49 * qJD(2) + t22 * qJD(4);
t6 = t385 * t376 + t382 - t479 / 0.4e1 + t595 - t410 + (t472 + t481) * pkin(5) / 0.2e1 + t426 + (pkin(5) * t485 + t395) * t591 + pkin(9) * t439 + pkin(5) * t440 + t379 * t445 + t375 * t570 + Ifges(6,3) * t561 + t432 * t557 + (t500 + t501) * t461 - Ifges(6,6) * t476 / 0.2e1;
t3 = t242 * t571 + t384 + t443 - t622;
t1 = t381 + t291 * t566 + mrSges(5,2) * t486 / 0.2e1 + mrSges(5,1) * t487 / 0.2e1 - mrSges(7,3) * t502 / 0.2e1 + t103 * t446 - (t601 + t245) * t487 / 0.2e1 - t610;
t20 = [t19 * qJD(2) + t25 * qJD(4) (-mrSges(3,2) * t427 - mrSges(3,1) * t434 + m(7) * (t103 * t115 + t104 * t116 - t291 * t474) + m(6) * (t176 * t236 + t177 * t237 - t250) + t177 * t333 + t104 * t261 + t103 * t263 + t176 * t335 + m(5) * (t292 * t360 - t359 * t488 - t250) + t292 * (-t380 * mrSges(5,1) + t376 * mrSges(5,2)) + m(4) * (-t291 * t373 - t292 * t506) * pkin(2) - t292 * mrSges(4,1) + t291 * mrSges(4,2) + (-t314 - t210) * t487 + (-t488 - t489) * mrSges(5,3)) * qJD(2) + t1 * qJD(4) + t3 * qJD(5) + t12 * qJD(6) + t420, t7, t1 * qJD(2) + t8 * qJD(5) + t18 * qJD(6) + (-t133 * t415 - t134 * t414) * t504 + t243 * t428 - t598 * qJD(4) * t242 + ((t133 * t430 + t134 * t256 + t243 * t363) * t591 + (-pkin(9) * t242 * t463 - pkin(4) * t243) * t593) * t594 + t419, t3 * qJD(2) + t8 * qJD(4) + (-t158 * mrSges(6,1) - t157 * mrSges(6,2) + (t374 * t431 - t378 * t96) * t589 + t30) * qJD(5) + t619, t12 * qJD(2) + t18 * qJD(4) + t30 * qJD(5) + t619; qJD(4) * t2 + qJD(5) * t4 + qJD(6) * t11 - t420, qJD(4) * t10 + qJD(5) * t13 + qJD(6) * t16, t460 - t503, t6 * qJD(5) + t14 * qJD(6) + (Ifges(5,5) + (-m(6) * pkin(4) + t508) * t359 + t407) * t462 + t404 + (t296 * t559 - Ifges(5,6) * t376 + t298 * t562 + t363 * t211 + t539 * t573 + (Ifges(7,1) * t415 - t540) * t572 + t202 * t569 + t308 * t245 - pkin(4) * t315 + t256 * t262 + t430 * t264 + mrSges(5,2) * t484 + m(7) * (t127 * t430 + t128 * t256 + t308 * t363) - t127 * t514 + (m(6) * t602 + t469 - t478) * pkin(9) + (Ifges(6,5) * t375 + Ifges(7,5) * t415 + Ifges(6,6) * t379) * t561 + t602 * mrSges(6,3) + (-t200 / 0.2e1 - Ifges(7,2) * t573 - t128 * mrSges(7,3) - Ifges(7,6) * t561) * t414) * qJD(4), t6 * qJD(4) + (-t237 * mrSges(6,1) - t236 * mrSges(6,2) - Ifges(6,5) * t477 - Ifges(6,6) * t473 + t465 - t527 + t528) * qJD(5) + t33 * qJD(6) + (m(7) * (t120 * t378 + t121 * t374) + (t471 + t480) * mrSges(7,3)) * t533 + t403, t14 * qJD(4) + t33 * qJD(5) + (t465 - t529 - t530) * qJD(6) + t402; t7, t460 + t503, t98 * qJD(4), t51 * qJD(5) + t56 * qJD(6) + (t301 * t415 + t302 * t414) * t504 + t598 * t462 + t376 * t428 + ((t463 * t553 - t552) * t593 + (-t256 * t302 - t301 * t430 + t483) * t591) * t594 + t401, t51 * qJD(4) + (t313 - t209 + (t299 * t378 - t300 * t374) * t589) * qJD(5) - t616 + t421, t56 * qJD(4) - qJD(5) * t209 + t482 - t616; -qJD(2) * t2 + qJD(5) * t9 + qJD(6) * t17 - t419, qJD(5) * t5 + qJD(6) * t15 - t404, -qJD(5) * t50 + qJD(6) * t55 - t401, -qJD(5) * t39 - qJD(6) * t52 (pkin(9) * t601 + t432 + t61) * qJD(5) + t618 + (m(7) * (-t256 * t378 + t374 * t430) + (-t374 * t415 + t378 * t414) * mrSges(7,3)) * t533 + t398, t61 * qJD(5) + t396 + t618; -qJD(2) * t4 - qJD(4) * t9 - qJD(6) * t29, -qJD(4) * t5 + qJD(6) * t34 - t403, t50 * qJD(4) - t421, -t398, -t327, -t327 - t392; -t11 * qJD(2) - t17 * qJD(4) + t29 * qJD(5), -qJD(4) * t15 - qJD(5) * t34 - t402, -t55 * qJD(4) - t482, -t396, t392, 0;];
Cq  = t20;
