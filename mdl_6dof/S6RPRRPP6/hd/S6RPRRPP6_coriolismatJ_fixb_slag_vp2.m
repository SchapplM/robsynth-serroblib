% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:26
% EndTime: 2019-03-09 04:46:46
% DurationCPUTime: 10.61s
% Computational Cost: add. (12732->655), mult. (26892->865), div. (0->0), fcn. (25971->6), ass. (0->343)
t654 = Ifges(6,1) + Ifges(7,1);
t377 = cos(qJ(4));
t375 = sin(qJ(4));
t529 = sin(pkin(9));
t460 = t529 * t375;
t530 = cos(pkin(9));
t320 = -t530 * t377 + t460;
t570 = pkin(4) * t377;
t359 = -pkin(3) - t570;
t462 = t530 * t375;
t398 = t377 * t529 + t462;
t187 = pkin(5) * t320 - qJ(6) * t398 + t359;
t543 = t398 * mrSges(7,3);
t546 = t320 * mrSges(7,1);
t216 = -t543 + t546;
t653 = m(7) * t187 + t216;
t643 = Ifges(7,4) + Ifges(6,5);
t650 = Ifges(6,6) - Ifges(7,6);
t378 = cos(qJ(3));
t569 = pkin(8) * t378;
t376 = sin(qJ(3));
t572 = pkin(3) * t376;
t337 = qJ(2) - t569 + t572;
t379 = -pkin(1) - pkin(7);
t366 = t376 * t379;
t486 = t377 * t366;
t395 = (qJ(5) * t378 - t337) * t375 - t486;
t649 = t395 * t529;
t648 = t395 * t530;
t373 = t375 ^ 2;
t374 = t377 ^ 2;
t493 = t373 + t374;
t644 = mrSges(6,3) + mrSges(7,2);
t314 = Ifges(6,4) * t320;
t561 = Ifges(7,5) * t320;
t647 = t398 * t654 - t314 + t561;
t290 = t398 * t376;
t231 = mrSges(7,2) * t290 + mrSges(7,3) * t378;
t232 = -mrSges(6,2) * t378 + mrSges(6,3) * t290;
t646 = t231 + t232;
t343 = pkin(3) * t378 + pkin(8) * t376;
t330 = t377 * t343;
t502 = t378 * t379;
t262 = -t375 * t502 + t330;
t263 = t375 * t343 + t377 * t502;
t416 = -t262 * t375 + t263 * t377;
t645 = m(6) + m(7);
t291 = t398 * t378;
t609 = mrSges(6,3) / 0.2e1;
t610 = mrSges(7,2) / 0.2e1;
t489 = t609 + t610;
t642 = t291 * t489;
t289 = t320 * t378;
t163 = -t289 * mrSges(7,1) + t291 * mrSges(7,3);
t164 = -t289 * mrSges(6,1) - t291 * mrSges(6,2);
t641 = t163 + t164;
t540 = t376 * mrSges(7,3);
t552 = t291 * mrSges(7,2);
t230 = t540 - t552;
t233 = -mrSges(6,2) * t376 - t291 * mrSges(6,3);
t498 = t230 + t233;
t371 = Ifges(5,4) * t377;
t640 = -Ifges(5,2) * t375 + t371;
t458 = mrSges(6,1) * t398 - t320 * mrSges(6,2);
t459 = mrSges(7,1) * t398 + t320 * mrSges(7,3);
t639 = t459 + t458;
t637 = t289 * t650 - t291 * t643;
t480 = t529 * pkin(4);
t353 = t480 + qJ(6);
t481 = t530 * pkin(4);
t358 = -t481 - pkin(5);
t635 = -t320 * t358 - t353 * t398;
t634 = -t289 * t353 + t291 * t358;
t319 = t377 * t337;
t463 = -t375 * t379 + pkin(4);
t504 = t377 * t378;
t485 = qJ(5) * t504;
t203 = t376 * t463 + t319 - t485;
t77 = t203 * t530 + t649;
t78 = t529 * t203 - t648;
t633 = t77 * t320 - t398 * t78;
t434 = -t375 * Ifges(5,1) - t371;
t539 = t377 * mrSges(5,1);
t542 = t375 * mrSges(5,2);
t632 = t542 - t539;
t465 = t230 / 0.2e1 + t233 / 0.2e1;
t613 = m(7) / 0.2e1;
t615 = m(6) / 0.2e1;
t631 = t615 + t613;
t630 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t559 = Ifges(5,6) * t375;
t563 = Ifges(5,5) * t377;
t629 = -t320 * t643 - t398 * t650 - t559 + t563;
t313 = Ifges(7,5) * t398;
t218 = t320 * Ifges(7,3) + t313;
t564 = Ifges(6,4) * t398;
t628 = -t320 * t654 + t218 + t313 - t564;
t268 = Ifges(7,5) * t289;
t154 = t376 * Ifges(7,6) + t291 * Ifges(7,3) - t268;
t565 = Ifges(6,4) * t289;
t627 = -t291 * t654 + t154 - t268 + t565;
t334 = mrSges(5,1) * t376 - mrSges(5,3) * t504;
t505 = t377 * t334;
t468 = -t505 / 0.2e1;
t537 = t378 * mrSges(5,3);
t626 = t468 - t493 * t537 / 0.2e1;
t625 = -Ifges(6,2) * t398 - t314 + t647;
t562 = Ifges(7,5) * t291;
t158 = -t289 * Ifges(7,1) + t376 * Ifges(7,4) + t562;
t271 = Ifges(6,4) * t291;
t160 = -t289 * Ifges(6,1) + t376 * Ifges(6,5) - t271;
t624 = Ifges(6,2) * t289 + t158 + t160 - t271;
t328 = m(7) * t353 + mrSges(7,3);
t370 = t377 * qJ(5);
t494 = t377 * pkin(8) + t370;
t567 = -qJ(5) - pkin(8);
t227 = -t567 * t462 + t494 * t529;
t623 = t567 * t460 + t494 * t530;
t134 = t623 * t289;
t525 = t227 * t291;
t585 = t398 / 0.2e1;
t589 = -t320 / 0.2e1;
t245 = -t366 * t375 + t319;
t222 = t245 - t485;
t83 = t222 * t529 - t648;
t84 = t222 * t530 + t649;
t622 = t83 * t585 + t134 / 0.2e1 - t525 / 0.2e1 + t84 * t589;
t612 = m(6) * pkin(4);
t621 = m(7) * t358 - t530 * t612 - mrSges(6,1) - mrSges(7,1);
t620 = t529 * t612 - mrSges(6,2) + t328;
t618 = 2 * qJD(3);
t617 = m(5) / 0.2e1;
t616 = -m(6) / 0.2e1;
t614 = -m(7) / 0.2e1;
t611 = mrSges(6,1) / 0.2e1;
t208 = t370 * t376 + t378 * t463 + t330;
t511 = t375 * t376;
t223 = qJ(5) * t511 + t263;
t81 = t208 * t530 - t223 * t529;
t73 = -t378 * pkin(5) - t81;
t608 = -t73 / 0.2e1;
t607 = m(7) * t83;
t126 = pkin(4) * t504 - pkin(5) * t289 + qJ(6) * t291;
t606 = t126 / 0.2e1;
t571 = pkin(4) * t375;
t194 = pkin(5) * t398 + qJ(6) * t320 + t571;
t605 = t194 / 0.2e1;
t604 = t623 / 0.2e1;
t292 = t376 * t320;
t234 = mrSges(6,1) * t378 - mrSges(6,3) * t292;
t601 = t234 / 0.2e1;
t600 = -t289 / 0.2e1;
t598 = t289 / 0.2e1;
t597 = -t290 / 0.2e1;
t596 = t290 / 0.2e1;
t595 = -t291 / 0.2e1;
t593 = t291 / 0.2e1;
t592 = t292 / 0.2e1;
t591 = -t292 / 0.2e1;
t315 = t632 * t378;
t590 = t315 / 0.2e1;
t587 = t320 / 0.2e1;
t586 = t320 / 0.4e1;
t583 = -t398 / 0.2e1;
t582 = -t375 / 0.2e1;
t581 = t375 / 0.2e1;
t580 = -t376 / 0.2e1;
t579 = t376 / 0.2e1;
t577 = -t377 / 0.2e1;
t576 = t377 / 0.2e1;
t575 = -t378 / 0.2e1;
t574 = t378 / 0.2e1;
t573 = m(7) * t292;
t566 = Ifges(5,4) * t375;
t558 = t289 * mrSges(6,2);
t557 = t289 * mrSges(7,3);
t556 = t290 * mrSges(6,1);
t555 = t290 * mrSges(7,1);
t554 = t291 * mrSges(6,1);
t553 = t291 * mrSges(7,1);
t551 = t292 * mrSges(6,2);
t550 = t292 * mrSges(7,2);
t549 = t292 * mrSges(7,3);
t323 = -pkin(4) * t511 + t366;
t106 = -pkin(5) * t290 - qJ(6) * t292 + t323;
t510 = t375 * t378;
t352 = pkin(4) * t510;
t325 = t352 - t502;
t107 = t291 * pkin(5) + t289 * qJ(6) + t325;
t153 = Ifges(7,5) * t292 + t378 * Ifges(7,6) - Ifges(7,3) * t290;
t155 = Ifges(6,4) * t292 + Ifges(6,2) * t290 + t378 * Ifges(6,6);
t156 = -t291 * Ifges(6,2) + t376 * Ifges(6,6) - t565;
t157 = Ifges(7,1) * t292 + t378 * Ifges(7,4) - Ifges(7,5) * t290;
t159 = Ifges(6,1) * t292 + Ifges(6,4) * t290 + t378 * Ifges(6,5);
t165 = -t549 - t555;
t166 = t551 - t556;
t167 = t553 + t557;
t168 = t554 - t558;
t538 = t378 * mrSges(7,1);
t235 = -t538 + t550;
t236 = mrSges(6,1) * t376 + t289 * mrSges(6,3);
t237 = -mrSges(7,1) * t376 - t289 * mrSges(7,2);
t246 = t375 * t337 + t486;
t283 = Ifges(5,6) * t378 - t376 * t640;
t435 = Ifges(5,1) * t377 - t566;
t285 = t378 * Ifges(5,5) - t376 * t435;
t436 = mrSges(5,1) * t375 + mrSges(5,2) * t377;
t316 = t376 * t436;
t331 = -mrSges(5,2) * t378 + mrSges(5,3) * t511;
t332 = -mrSges(5,2) * t376 - mrSges(5,3) * t510;
t509 = t376 * t377;
t333 = mrSges(5,1) * t378 + mrSges(5,3) * t509;
t405 = t563 / 0.2e1 - t559 / 0.2e1 - Ifges(4,4);
t487 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t488 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t407 = t435 * t378;
t286 = t376 * Ifges(5,5) + t407;
t507 = t377 * t286;
t284 = t376 * Ifges(5,6) + t378 * t640;
t514 = t375 * t284;
t70 = t376 * qJ(6) + t78;
t71 = -t376 * pkin(5) - t77;
t82 = t529 * t208 + t530 * t223;
t72 = qJ(6) * t378 + t82;
t3 = (qJ(2) * mrSges(4,1) + t283 * t582 + t285 * t576 + t379 * t316 + t405 * t378 - t487 * t291 + t488 * t289 + (-Ifges(4,1) + Ifges(4,2) + (-m(5) * t379 + t436) * t379 + t630) * t376) * t378 + (t514 / 0.2e1 - qJ(2) * mrSges(4,2) - t507 / 0.2e1 - t405 * t376 - t488 * t292 + t487 * t290) * t376 - (-t158 / 0.2e1 - t160 / 0.2e1) * t292 + (-t154 / 0.2e1 + t156 / 0.2e1) * t290 + m(5) * (t245 * t262 + t246 * t263) + (-t155 / 0.2e1 + t153 / 0.2e1) * t291 + (-t157 / 0.2e1 - t159 / 0.2e1) * t289 + m(6) * (t323 * t325 + t77 * t81 + t78 * t82) + m(7) * (t106 * t107 + t70 * t72 + t71 * t73) + t246 * t331 + t263 * t332 + t245 * t333 + t262 * t334 + t323 * t168 + t325 * t166 + t81 * t236 + t73 * t237 + t72 * t230 + t70 * t231 + t78 * t232 + t82 * t233 + t77 * t234 + t71 * t235 + t107 * t165 + t106 * t167;
t548 = t3 * qJD(1);
t547 = t320 * mrSges(6,1);
t545 = t320 * mrSges(7,2);
t544 = t398 * mrSges(6,2);
t417 = t245 * t375 - t246 * t377;
t422 = t289 * t78 + t291 * t77;
t423 = t289 * t70 - t291 * t71;
t426 = -Ifges(7,3) * t289 - t562;
t428 = Ifges(5,5) * t375 + Ifges(5,6) * t377;
t339 = t377 * Ifges(5,2) + t566;
t473 = t339 * t581;
t497 = t236 - t237;
t4 = t245 * t332 - t246 * t334 + t325 * t164 + t426 * t593 + t156 * t598 + t107 * t163 + t126 * t167 + t498 * t84 - t497 * t83 + t422 * mrSges(6,3) + t423 * mrSges(7,2) + m(7) * (t107 * t126 + t70 * t84 + t71 * t83) + m(6) * (-t77 * t83 + t78 * t84) + (t428 * t580 + t286 * t582 + t284 * t577 + t379 * t315 + (t434 * t576 + t473) * t378 + (m(6) * t325 + t168) * t570 + t417 * mrSges(5,3)) * t378 + t637 * t579 + t627 * t600 + t624 * t595;
t536 = t4 * qJD(1);
t535 = t70 * t398;
t534 = t71 * t320;
t531 = t632 - mrSges(4,1);
t28 = t376 * t230 + m(7) * (t107 * t289 + t376 * t70) + t289 * t167;
t528 = qJD(1) * t28;
t39 = t631 * (-t289 * t290 + t291 * t292);
t527 = qJD(1) * t39;
t513 = t375 * t332;
t14 = t376 * mrSges(4,1) + t378 * mrSges(4,2) + t513 + t505 + mrSges(3,3) + t498 * t398 - t497 * t320 + (m(4) + m(3)) * qJ(2) + m(7) * (t534 + t535) - m(6) * t633 + m(5) * (t245 * t377 + t246 * t375);
t526 = t14 * qJD(1);
t520 = t290 * t291;
t185 = t292 * t289;
t518 = t292 * t376;
t512 = t375 * t333;
t508 = t376 * t378;
t506 = t377 * t331;
t503 = t378 * t289;
t85 = 0.2e1 * (t586 + t503 / 0.4e1 + t518 / 0.4e1) * m(7);
t501 = t85 * qJD(1);
t492 = qJD(3) * t376;
t491 = t612 / 0.2e1;
t490 = t571 / 0.2e1;
t475 = t398 * t575;
t474 = -t513 / 0.2e1;
t472 = -t511 / 0.2e1;
t471 = -t510 / 0.2e1;
t469 = t320 * t580;
t467 = -t504 / 0.2e1;
t466 = t504 / 0.2e1;
t464 = -t237 / 0.2e1 + t236 / 0.2e1;
t450 = -t227 * t289 - t291 * t623;
t449 = -t134 + t525;
t443 = pkin(4) * t466;
t442 = mrSges(6,3) * t481;
t441 = mrSges(6,3) * t480;
t438 = t489 * t289;
t425 = Ifges(7,3) * t398 - t561;
t424 = t227 * t83 + t623 * t84;
t32 = m(5) * (-0.1e1 + t493) * t508 + t645 * (t185 - t508 + t520);
t382 = -m(5) * (-t417 * t378 + (t416 - 0.2e1 * t502) * t376) / 0.2e1 + (-t290 * t81 - t292 * t82 - t323 * t378 + t325 * t376 - t422) * t616 + (-t106 * t378 + t107 * t376 + t290 * t73 - t292 * t72 - t423) * t614;
t392 = t631 * (t227 * t320 + t398 * t623);
t5 = -(-t231 / 0.2e1 - t232 / 0.2e1) * t292 + t464 * t291 + (-t235 / 0.2e1 + t601) * t290 + t465 * t289 + (-t167 / 0.2e1 - t168 / 0.2e1 + t512 / 0.2e1 - t506 / 0.2e1) * t376 + (t165 / 0.2e1 + t436 * t580 - t316 / 0.2e1 + t166 / 0.2e1 + t334 * t581 + t332 * t577) * t378 + t382 + t392;
t421 = -t5 * qJD(1) + t32 * qJD(2);
t386 = -t635 * t613 - t547 / 0.2e1 - t546 / 0.2e1 - t544 / 0.2e1 + t543 / 0.2e1 - t542 / 0.2e1 + t539 / 0.2e1 + (-t320 * t530 + t398 * t529) * t491;
t56 = t292 * t84;
t387 = (-t378 ^ 2 * t570 + t292 * t77 - t56 + (-t78 + t83) * t290) * t615 + (-t126 * t378 - t292 * t71 - t56 + (-t70 + t83) * t290) * t613;
t8 = t236 * t592 + t237 * t591 + t332 * t472 - t386 + t387 + t498 * t597 + (-t315 + t641) * t575 + t626 * t376 + t644 * (-t185 / 0.2e1 - t520 / 0.2e1);
t420 = t8 * qJD(1);
t16 = -t498 * t291 + t497 * t289 + m(7) * (-t289 * t71 - t291 * t70) + m(6) * (t289 * t77 - t291 * t78);
t419 = qJD(1) * t16 + qJD(2) * t39;
t389 = (-t289 * t358 - t291 * t353) * t613 + (t289 * t530 - t291 * t529) * t491;
t404 = m(6) * t443 + m(7) * t606;
t27 = -t389 + t404 + t641;
t388 = (-t320 * t353 + t358 * t398) * t613 + (-t320 * t529 - t398 * t530) * t491;
t406 = m(6) * t490 + m(7) * t605;
t34 = -t388 + t406 + t639;
t418 = qJD(1) * t27 + qJD(3) * t34;
t129 = m(7) * t289;
t189 = m(7) * t398;
t414 = qJD(1) * t129 - qJD(3) * t189;
t412 = t106 * t614 + t323 * t616;
t411 = m(7) * t608 + t538 / 0.2e1;
t409 = -t434 * t577 + t473;
t408 = t378 * t436;
t217 = t544 + t547;
t219 = -t320 * Ifges(6,2) + t564;
t402 = -t408 / 0.2e1;
t380 = -t514 / 0.4e1 + t507 / 0.4e1 + t633 * t609 + (t627 / 0.4e1 - t156 / 0.4e1) * t398 + (-t623 * t77 + (t325 * t375 + t359 * t504) * pkin(4) + t424) * t615 + (t107 * t194 + t126 * t187 + t623 * t71 + t424) * t613 - t623 * t236 / 0.2e1 + t107 * t459 / 0.2e1 + t379 * t402 + t325 * t458 / 0.2e1 + (-t625 / 0.4e1 + t425 / 0.4e1) * t291 + (-t628 / 0.4e1 + t219 / 0.4e1) * t289 - (t613 * t70 + t615 * t78 + t465) * t227 + t168 * t490 + t237 * t604 + t167 * t605 + t216 * t606 - t535 * t610 + pkin(3) * t590 + t426 * t586 - t624 * t320 / 0.4e1 + (t474 + t626) * pkin(8) + t629 * t376 / 0.4e1 + (-t534 / 0.2e1 + t622) * mrSges(7,2) + t622 * mrSges(6,3) + t359 * t164 / 0.2e1 + t217 * t443 + t187 * t163 / 0.2e1 + t377 * t407 / 0.4e1 + t339 * t467 + (-t640 / 0.4e1 + t434 / 0.2e1) * t510;
t381 = (t353 * t72 + t358 * t73) * t613 + Ifges(6,6) * t596 + Ifges(7,6) * t597 + t262 * mrSges(5,1) / 0.2e1 - t263 * mrSges(5,2) / 0.2e1 + t353 * t231 / 0.2e1 + t358 * t235 / 0.2e1 + t72 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t608 + t81 * t611 - t82 * mrSges(6,2) / 0.2e1 + (t529 * t82 + t530 * t81) * t491 - Ifges(5,5) * t509 / 0.2e1 + Ifges(5,6) * t511 / 0.2e1 + t481 * t601 + t232 * t480 / 0.2e1 + t643 * t592 + t630 * t574;
t1 = -t380 + t381;
t383 = -t194 * t378 * t613 - t352 * t615 + t639 * t575 + t402;
t385 = t634 * t613 + t558 / 0.2e1 - t557 / 0.2e1 - t554 / 0.2e1 - t553 / 0.2e1 + (-t289 * t529 - t291 * t530) * t491 + mrSges(5,1) * t471 + mrSges(5,2) * t467;
t12 = -t383 + t385;
t9 = -pkin(3) * t436 + t187 * t459 + t217 * t571 + t219 * t583 + t425 * t587 + t435 * t581 + t640 * t576 + t628 * t585 + t625 * t589 - t409 + (m(6) * t571 + t458) * t359 + t653 * t194;
t403 = -t1 * qJD(1) - t12 * qJD(2) + t9 * qJD(3);
t384 = (-t465 + t642) * t320 - (t438 + t464) * t398 + (-t320 * t78 - t398 * t77 + t450) * t615 + (-t320 * t70 + t398 * t71 + t450) * t613;
t11 = -(mrSges(6,2) / 0.2e1 - mrSges(7,3) / 0.2e1) * t292 + (mrSges(7,1) / 0.2e1 + t611) * t290 + t384 + t412;
t26 = (t320 ^ 2 + t398 ^ 2) * t644 + t645 * (t227 * t398 - t320 * t623);
t391 = t631 * (t290 * t398 + t320 * t292);
t43 = (t616 + t614) * t376 + t391;
t401 = qJD(1) * t11 + qJD(2) * t43 + qJD(3) * t26;
t162 = (-t475 + t595) * m(7);
t390 = (-t107 * t398 + t187 * t289 + t376 * t623) * t613 + t216 * t598 + t167 * t583;
t22 = (t469 + t591) * mrSges(7,2) + t390 + t411;
t46 = t653 * t398;
t400 = qJD(1) * t22 + qJD(2) * t162 - qJD(3) * t46;
t393 = t540 + (t353 * t376 + t70) * t613;
t31 = -t607 / 0.2e1 + t393;
t399 = qJD(1) * t31 + qJD(4) * t328;
t161 = (-t475 + t593) * m(7);
t86 = (-t503 - t518) * t613 + m(7) * t587;
t47 = m(7) * t604 + t613 * t623 - t545;
t44 = t388 + t406;
t42 = t579 * t645 + t391;
t37 = t39 * qJD(5);
t35 = t389 + t404;
t30 = -t552 + t607 / 0.2e1 + t393;
t21 = mrSges(7,2) * t469 + t550 / 0.2e1 + t390 - t411;
t13 = t383 + t385;
t10 = -t549 / 0.2e1 - t555 / 0.2e1 + t551 / 0.2e1 - t556 / 0.2e1 + t384 - t412;
t7 = (-t164 / 0.2e1 - t163 / 0.2e1 + t590) * t378 - (t438 - t464) * t292 + (-t465 - t642) * t290 + (t474 + t468 + (-t374 / 0.2e1 - t373 / 0.2e1) * t537) * t376 + t386 + t387;
t6 = t234 * t597 + t235 * t596 + t236 * t595 + t237 * t593 + t332 * t466 + t333 * t472 + t334 * t471 - t382 + t392 + t498 * t600 + t646 * t591 + (t165 + t166 - t316) * t575 + (t167 + t168 + t408 + t506) * t579;
t2 = t380 + t381;
t15 = [qJD(2) * t14 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t16 + qJD(6) * t28, t6 * qJD(3) + t7 * qJD(4) + t86 * qJD(6) + t37 + t526 + 0.2e1 * t631 * qJD(2) * (t290 * t320 - t292 * t398) t548 + t6 * qJD(2) + t2 * qJD(4) + t10 * qJD(5) + t21 * qJD(6) + (t379 * t531 - Ifges(4,5) + t409) * t492 + ((-t227 * t81 + t323 * t359 + t623 * t82) * t615 + (t106 * t187 + t227 * t73 + t623 * t72) * t613 + (-pkin(3) * t366 + pkin(8) * t416) * t617) * t618 + (t73 * t398 * mrSges(7,2) - t72 * t545 + t155 * t589 + t219 * t596 + t218 * t597 + t285 * t581 + t153 * t587 - Ifges(4,6) * t378 + t359 * t166 + t323 * t217 + pkin(3) * t316 + t106 * t216 + t187 * t165 - mrSges(4,2) * t502 + t283 * t576 + (-t320 * t650 + t643 * t398 + t428) * t574 + t647 * t592 + (t159 + t157) * t585 + t646 * t623 + (-t234 + t235) * t227 + (-t512 + t506) * pkin(8) + (-t320 * t82 - t398 * t81) * mrSges(6,3) + t416 * mrSges(5,3)) * qJD(3), t536 + t7 * qJD(2) + t2 * qJD(3) + (-t246 * mrSges(5,1) - t245 * mrSges(5,2) - mrSges(7,2) * t634 - Ifges(5,5) * t510 - Ifges(5,6) * t504 + t289 * t441 + t291 * t442 + t620 * t84 + t621 * t83 + t637) * qJD(4) + t35 * qJD(5) + t30 * qJD(6), qJD(3) * t10 + qJD(4) * t35 + t419, qJD(2) * t86 + qJD(3) * t21 + qJD(4) * t30 + t528; -qJD(3) * t5 + qJD(4) * t8 - qJD(6) * t85 + t37 - t526, qJD(3) * t32, t13 * qJD(4) + t42 * qJD(5) + t161 * qJD(6) + (t216 + t217 + t531) * t492 + ((t493 * t569 - t572) * t617 + (t359 * t376 + t449) * t615 + (t187 * t376 + t449) * t613) * t618 + t421 + ((mrSges(5,3) * t493 - mrSges(4,2)) * t378 + t644 * (t289 * t320 + t291 * t398)) * qJD(3), t13 * qJD(3) + (-mrSges(5,1) * t509 + mrSges(5,2) * t511 - t290 * t620 - t292 * t621) * qJD(4) - qJD(6) * t573 + t420, qJD(3) * t42 + t527, t161 * qJD(3) - qJD(4) * t573 - t501; qJD(2) * t5 - qJD(4) * t1 + qJD(5) * t11 + qJD(6) * t22 - t548, -qJD(4) * t12 + qJD(5) * t43 + qJD(6) * t162 - t421, qJD(4) * t9 + qJD(5) * t26 - qJD(6) * t46 (mrSges(7,2) * t635 + pkin(8) * t632 - t227 * t620 + t320 * t442 - t398 * t441 + t621 * t623 + t629) * qJD(4) + t44 * qJD(5) + t47 * qJD(6) + t403, qJD(4) * t44 + t401, qJD(4) * t47 + t400; -qJD(2) * t8 + qJD(3) * t1 - qJD(5) * t27 + qJD(6) * t31 - t536, qJD(3) * t12 - t420, -qJD(5) * t34 - t403, t328 * qJD(6), -t418, t399; -qJD(3) * t11 + qJD(4) * t27 + qJD(6) * t129 - t419, -qJD(3) * t43 - t527, qJD(4) * t34 - qJD(6) * t189 - t401, t418, 0, t414; qJD(2) * t85 - qJD(3) * t22 - qJD(4) * t31 - qJD(5) * t129 - t528, -qJD(3) * t162 + t501, qJD(5) * t189 - t400, -t399, -t414, 0;];
Cq  = t15;
