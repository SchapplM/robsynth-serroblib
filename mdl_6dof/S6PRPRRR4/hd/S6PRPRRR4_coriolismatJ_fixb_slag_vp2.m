% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRR4
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:36:14
% EndTime: 2019-03-08 20:36:29
% DurationCPUTime: 9.93s
% Computational Cost: add. (24981->614), mult. (55617->870), div. (0->0), fcn. (64865->12), ass. (0->334)
t381 = sin(pkin(12));
t383 = cos(pkin(12));
t576 = sin(qJ(4));
t579 = cos(qJ(4));
t347 = t381 * t576 - t383 * t579;
t349 = -t381 * t579 - t383 * t576;
t575 = pkin(4) * t349;
t296 = pkin(9) * t347 - t575;
t563 = pkin(8) + qJ(3);
t355 = t563 * t383;
t445 = t563 * t381;
t307 = t355 * t576 + t579 * t445;
t385 = sin(qJ(5));
t386 = cos(qJ(5));
t199 = t386 * t296 + t307 * t385;
t200 = t385 * t296 - t386 * t307;
t384 = sin(qJ(6));
t578 = cos(qJ(6));
t418 = t384 * t386 + t385 * t578;
t264 = t418 * t347;
t213 = mrSges(7,2) * t349 + t264 * mrSges(7,3);
t641 = -t384 * t385 + t578 * t386;
t266 = t641 * t347;
t215 = -mrSges(7,1) * t349 + t266 * mrSges(7,3);
t491 = pkin(5) * t578;
t437 = t491 / 0.2e1;
t515 = t347 * t385;
t456 = t515 / 0.2e1;
t374 = Ifges(6,5) * t386;
t482 = -t374 / 0.2e1;
t624 = m(7) * pkin(5);
t494 = t624 / 0.2e1;
t574 = pkin(5) * t384;
t591 = -t349 / 0.2e1;
t622 = -mrSges(6,2) / 0.2e1;
t623 = mrSges(6,1) / 0.2e1;
t606 = -t266 / 0.2e1;
t609 = t264 / 0.2e1;
t421 = Ifges(7,5) * t606 + Ifges(7,6) * t609;
t675 = mrSges(7,2) / 0.2e1;
t514 = t347 * t386;
t136 = -pkin(5) * t349 + pkin(10) * t514 + t199;
t163 = pkin(10) * t515 + t200;
t81 = t136 * t578 - t384 * t163;
t82 = t384 * t136 + t163 * t578;
t633 = Ifges(7,3) * t349 / 0.2e1 + t82 * t675 - t81 * mrSges(7,1) / 0.2e1 - t421;
t680 = Ifges(6,3) * t591 + t199 * t623 + t200 * t622 + (t384 * t82 + t578 * t81) * t494 + t347 * t482 + Ifges(6,6) * t456 + t213 * t574 / 0.2e1 + t215 * t437 - t633;
t529 = t386 * mrSges(6,2);
t470 = t529 / 0.2e1;
t540 = t266 * mrSges(7,2);
t541 = t264 * mrSges(7,1);
t647 = t541 / 0.2e1 + t540 / 0.2e1;
t679 = -(t264 * t578 - t266 * t384) * t494 - mrSges(6,1) * t456 - t347 * t470 - t647;
t382 = sin(pkin(6));
t577 = sin(qJ(2));
t468 = t382 * t577;
t528 = cos(pkin(6));
t337 = t381 * t528 + t383 * t468;
t412 = t381 * t468 - t383 * t528;
t259 = t337 * t576 + t412 * t579;
t145 = t418 * t259;
t146 = t641 * t259;
t522 = t259 * t385;
t459 = t522 / 0.2e1;
t676 = mrSges(7,1) / 0.2e1;
t649 = t145 * t676 + t146 * t675;
t678 = -(t145 * t578 - t146 * t384) * t494 - mrSges(6,1) * t459 - t259 * t470 - t649;
t387 = cos(qJ(2));
t511 = t382 * t387;
t322 = t347 * t511;
t279 = t385 * t322 + t386 * t468;
t280 = -t386 * t322 + t385 * t468;
t169 = t279 * t578 - t384 * t280;
t170 = t384 * t279 + t280 * t578;
t648 = t169 * t676 - t170 * mrSges(7,2) / 0.2e1;
t677 = -t279 * t623 - t280 * t622 - (t169 * t578 + t170 * t384) * t494 - t648;
t379 = t385 ^ 2;
t380 = t386 ^ 2;
t645 = t380 + t379;
t446 = pkin(9) * t645;
t370 = -pkin(3) * t383 - pkin(2);
t287 = pkin(4) * t347 + pkin(9) * t349 + t370;
t308 = t355 * t579 - t445 * t576;
t187 = t386 * t287 - t308 * t385;
t512 = t349 * t386;
t153 = pkin(10) * t512 + t187;
t132 = pkin(5) * t347 + t153;
t188 = t287 * t385 + t308 * t386;
t513 = t349 * t385;
t154 = pkin(10) * t513 + t188;
t509 = t384 * t154;
t76 = t132 * t578 - t509;
t86 = t153 * t578 - t509;
t674 = -t76 + t86;
t673 = mrSges(6,3) * t645;
t260 = t337 * t579 - t412 * t576;
t226 = -t385 * t260 - t386 * t511;
t227 = t260 * t386 - t385 * t511;
t129 = t384 * t226 + t227 * t578;
t441 = t578 * t226 - t227 * t384;
t39 = -t129 * mrSges(7,1) - t441 * mrSges(7,2);
t672 = t39 * qJD(6);
t621 = -pkin(10) - pkin(9);
t362 = t621 * t385;
t363 = t621 * t386;
t318 = t384 * t362 - t363 * t578;
t440 = t578 * t362 + t363 * t384;
t499 = Ifges(7,5) * t641 - Ifges(7,6) * t418;
t59 = -t318 * mrSges(7,1) - t440 * mrSges(7,2) + t499;
t671 = t59 * qJD(6);
t265 = t418 * t349;
t557 = mrSges(7,3) * t265;
t214 = -mrSges(7,2) * t347 + t557;
t263 = t641 * t349;
t216 = mrSges(7,1) * t347 + mrSges(7,3) * t263;
t586 = -t418 / 0.2e1;
t587 = t418 / 0.2e1;
t589 = t641 / 0.2e1;
t590 = -t641 / 0.2e1;
t396 = (t263 * t587 + t265 * t590) * mrSges(7,3) + t214 * t589 + t216 * t586;
t164 = -mrSges(7,1) * t263 + mrSges(7,2) * t265;
t608 = -t265 / 0.2e1;
t610 = t263 / 0.2e1;
t613 = t259 / 0.2e1;
t616 = t214 / 0.2e1;
t670 = (t129 * t610 + t441 * t608) * mrSges(7,3) + t164 * t613 + t441 * t616;
t593 = t318 / 0.2e1;
t467 = t578 * t154;
t77 = t384 * t132 + t467;
t85 = -t384 * t153 - t467;
t669 = t85 + t77;
t531 = t641 * mrSges(7,3);
t498 = t381 ^ 2 + t383 ^ 2;
t666 = t498 * mrSges(4,3);
t292 = -t347 * mrSges(6,2) + mrSges(6,3) * t513;
t502 = t386 * t292;
t294 = mrSges(6,1) * t347 + mrSges(6,3) * t512;
t505 = t385 * t294;
t419 = -t505 / 0.2e1 + t502 / 0.2e1;
t625 = m(7) / 0.2e1;
t665 = (t418 * t674 + t669 * t641) * t625 + t419 + t396;
t356 = -mrSges(6,1) * t386 + mrSges(6,2) * t385;
t281 = t349 * t356;
t524 = t227 * t386;
t525 = t226 * t385;
t558 = mrSges(6,3) * t349;
t604 = -t294 / 0.2e1;
t614 = t216 / 0.2e1;
t664 = (-pkin(5) * t259 * t512 + t669 * t441) * t625 + t226 * t292 / 0.2e1 + t227 * t604 + t281 * t613 + (t524 / 0.2e1 - t525 / 0.2e1) * t558 + (t625 * t674 - t614) * t129 + t670;
t297 = mrSges(7,1) * t418 + mrSges(7,2) * t641;
t662 = t297 * qJD(6);
t166 = -mrSges(7,1) * t265 - mrSges(7,2) * t263;
t375 = Ifges(6,4) * t386;
t360 = Ifges(6,1) * t385 + t375;
t285 = t349 * t360;
t240 = -pkin(5) * t513 + t307;
t572 = pkin(5) * t386;
t373 = -pkin(4) - t572;
t592 = t347 / 0.4e1;
t603 = t297 / 0.2e1;
t409 = t240 * t603 + t499 * t592 + t373 * t164 / 0.2e1;
t551 = Ifges(6,6) * t385;
t444 = t374 - t551;
t553 = Ifges(7,4) * t418;
t300 = Ifges(7,2) * t641 + t553;
t301 = Ifges(7,1) * t641 - t553;
t447 = t300 / 0.4e1 - t301 / 0.4e1;
t345 = Ifges(7,4) * t641;
t299 = -Ifges(7,2) * t418 + t345;
t302 = Ifges(7,1) * t418 + t345;
t448 = t299 / 0.4e1 + t302 / 0.4e1;
t431 = Ifges(6,2) * t385 - t375;
t533 = t347 * Ifges(6,6);
t231 = t349 * t431 + t533;
t506 = t385 * t231;
t571 = pkin(9) * t385;
t573 = pkin(5) * t385;
t357 = t385 * mrSges(6,1) + t529;
t585 = t357 / 0.2e1;
t595 = t440 / 0.2e1;
t661 = -pkin(4) * t281 / 0.2e1 + t307 * t585 - t318 * t614 + t214 * t595 + t444 * t592 - t506 / 0.4e1 + t385 * t285 / 0.4e1 + t166 * t573 / 0.2e1 - t292 * t571 / 0.2e1 + t409 + t448 * t265 + t447 * t263;
t626 = -m(7) / 0.2e1;
t628 = -m(6) / 0.2e1;
t660 = -(t279 * t386 + t280 * t385) * t628 - (t169 * t641 + t170 * t418) * t626 + (m(4) + m(5)) * t468 / 0.2e1;
t530 = t418 * mrSges(7,3);
t655 = Ifges(5,4) + t482 + t551 / 0.2e1;
t654 = t356 - mrSges(5,1);
t283 = t357 * t349;
t646 = -t283 + t166;
t644 = t383 * t337 + t381 * t412;
t256 = Ifges(7,4) * t265;
t140 = -Ifges(7,1) * t263 + t347 * Ifges(7,5) + t256;
t167 = Ifges(7,2) * t263 + t256;
t554 = Ifges(7,4) * t263;
t138 = Ifges(7,2) * t265 + t347 * Ifges(7,6) - t554;
t168 = Ifges(7,1) * t265 + t554;
t433 = (-t168 / 0.4e1 + t138 / 0.4e1) * t418;
t643 = -t433 + (t140 / 0.4e1 + t167 / 0.4e1) * t641;
t596 = -t440 / 0.2e1;
t642 = t263 * t593 + t265 * t596;
t555 = Ifges(6,4) * t385;
t358 = Ifges(6,2) * t386 + t555;
t284 = t349 * t358;
t637 = pkin(9) * t604 + t284 / 0.4e1;
t636 = -t199 * t385 + t200 * t386;
t635 = -m(7) * t240 - t166;
t629 = 2 * qJD(4);
t627 = m(6) / 0.2e1;
t620 = -t441 / 0.2e1;
t617 = -t214 / 0.2e1;
t615 = -t216 / 0.2e1;
t611 = -t263 / 0.2e1;
t607 = t265 / 0.2e1;
t601 = t300 / 0.2e1;
t598 = t302 / 0.2e1;
t594 = -t318 / 0.2e1;
t584 = t358 / 0.4e1;
t361 = Ifges(6,1) * t386 - t555;
t583 = -t361 / 0.4e1;
t582 = t385 / 0.2e1;
t581 = -t386 / 0.2e1;
t580 = t386 / 0.2e1;
t570 = t76 * mrSges(7,2);
t569 = t77 * mrSges(7,1);
t566 = t85 * mrSges(7,1);
t565 = t86 * mrSges(7,2);
t561 = m(7) * qJD(3);
t535 = t347 * mrSges(5,3);
t534 = t347 * Ifges(6,5);
t532 = t349 * mrSges(5,3);
t523 = t240 * t385;
t321 = t349 * t511;
t198 = t259 * t321;
t27 = t396 - t647;
t521 = t27 * qJD(2);
t520 = t279 * t385;
t519 = t280 * t386;
t518 = t307 * t321;
t517 = t307 * t349;
t469 = t382 ^ 2 * t577;
t34 = m(7) * (t129 * t170 + t169 * t441 - t198) + m(6) * (t226 * t279 + t227 * t280 - t198) + m(5) * (-t260 * t322 - t387 * t469 - t198) + m(4) * (t382 * t644 - t469) * t387;
t516 = t34 * qJD(1);
t225 = t349 * t259;
t508 = t384 * t216;
t504 = t385 * t358;
t233 = -t349 * t361 + t534;
t503 = t386 * t233;
t501 = Ifges(7,5) * t265 + Ifges(7,6) * t263;
t495 = mrSges(7,3) * t574;
t63 = t145 * t641 - t146 * t418;
t493 = t63 * qJD(4) * t625;
t486 = t85 / 0.2e1 + t77 / 0.2e1;
t485 = t86 / 0.2e1 - t76 / 0.2e1;
t473 = -t531 / 0.2e1;
t472 = t531 / 0.2e1;
t471 = -t530 / 0.2e1;
t460 = t259 * t603;
t298 = -mrSges(7,1) * t641 + mrSges(7,2) * t418;
t455 = t298 * t591;
t451 = t620 + t441 / 0.2e1;
t449 = t298 / 0.2e1 + t356 / 0.2e1;
t443 = t498 * qJ(3);
t442 = t263 * t495;
t439 = mrSges(7,3) * t491;
t435 = (t380 / 0.2e1 + t379 / 0.2e1) * mrSges(6,3);
t432 = t265 * t439;
t137 = -Ifges(7,4) * t266 + Ifges(7,2) * t264 - Ifges(7,6) * t349;
t139 = -Ifges(7,1) * t266 + Ifges(7,4) * t264 - Ifges(7,5) * t349;
t165 = -t540 - t541;
t230 = -Ifges(6,6) * t349 + t347 * t431;
t232 = -Ifges(6,5) * t349 - t347 * t361;
t241 = -pkin(5) * t515 + t308;
t282 = t357 * t347;
t291 = mrSges(6,2) * t349 + mrSges(6,3) * t515;
t293 = -mrSges(6,1) * t349 + mrSges(6,3) * t514;
t340 = t349 * mrSges(5,1);
t3 = -t370 * t340 - t307 * t282 - t308 * t283 + t199 * t294 + t188 * t291 + t200 * t292 + t187 * t293 + t138 * t609 + t137 * t607 + t140 * t606 + t139 * t611 + t240 * t165 + t241 * t166 + t81 * t216 + t77 * t213 + t82 * t214 + t76 * t215 + m(7) * (t240 * t241 + t76 * t81 + t77 * t82) + m(6) * (t187 * t199 + t188 * t200 + t307 * t308) + (-t370 * mrSges(5,2) + t506 / 0.2e1 - t503 / 0.2e1 + t421 + t655 * t347) * t347 + (Ifges(7,5) * t610 + Ifges(7,6) * t608 + t230 * t582 + t232 * t581 - t655 * t349 + (-Ifges(6,3) - Ifges(7,3) + Ifges(5,1) - Ifges(5,2)) * t347) * t349;
t424 = t187 * t385 - t188 * t386;
t388 = (t199 * t226 + t200 * t227 + t260 * t307 + (t308 + t424) * t259) * t628 + (t129 * t82 + t145 * t76 - t146 * t77 + t240 * t260 + t241 * t259 + t441 * t81) * t626 + t215 * t620 - t129 * t213 / 0.2e1 + t145 * t615 - t146 * t617 - t226 * t293 / 0.2e1 - t227 * t291 / 0.2e1 + (-t347 * mrSges(5,2) - t340) * t511 / 0.2e1;
t389 = (t169 * t586 + t170 * t589) * mrSges(7,3) + (-t520 / 0.2e1 + t519 / 0.2e1) * mrSges(6,3) - (-mrSges(5,1) / 0.2e1 + t449) * t321 + (pkin(4) * t321 + (t519 - t520) * pkin(9)) * t627 + (t169 * t440 + t170 * t318 - t321 * t373) * t625 + t322 * mrSges(5,2) / 0.2e1;
t4 = t388 + (t283 / 0.2e1 - t166 / 0.2e1) * t260 + (t282 / 0.2e1 - t165 / 0.2e1 + t419) * t259 + t389;
t430 = -t4 * qJD(1) + t3 * qJD(2);
t399 = t240 * t164 + (t140 / 0.2e1 + t167 / 0.2e1) * t265 + (-t168 / 0.2e1 + t138 / 0.2e1 + t77 * mrSges(7,3)) * t263 - t76 * t557 + t347 * t501 / 0.2e1;
t6 = t85 * t216 + m(7) * (t76 * t85 + t77 * t86) + t86 * t214 + t187 * t292 - t188 * t294 + t307 * t281 + ((t534 / 0.2e1 - t187 * mrSges(6,3) + t233 / 0.2e1 + t284 / 0.2e1) * t385 + (t533 / 0.2e1 + t188 * mrSges(6,3) - t285 / 0.2e1 + t231 / 0.2e1 + t635 * pkin(5)) * t386) * t349 + t399;
t7 = t664 + t677;
t429 = t7 * qJD(1) + t6 * qJD(2);
t11 = t76 * t214 - t77 * t216 + t399;
t401 = t129 * t615 + t670;
t12 = t401 - t648;
t428 = t12 * qJD(1) + t11 * qJD(2);
t17 = t665 + t679;
t427 = t17 * qJD(2);
t21 = -t266 * t214 + t264 * t216 + t666 + (t532 - t646) * t349 + (-t502 + t505 + t535) * t347 + m(7) * (-t240 * t349 + t264 * t76 - t266 * t77) + m(6) * (t347 * t424 - t517) + m(5) * (-t308 * t347 - t517) + m(4) * t443;
t423 = -t524 + t525;
t392 = (-t129 * t266 + t264 * t441 - t225) * t625 + (t347 * t423 - t225) * t627 + m(5) * (-t260 * t347 - t225) / 0.2e1 + m(4) * t644 / 0.2e1;
t31 = t392 - t660;
t426 = t31 * qJD(1) + t21 * qJD(2);
t395 = (t199 * t386 + t200 * t385) * t627 + (t418 * t82 + t641 * t81) * t625 + t215 * t589 + t213 * t587 + t291 * t582 + t293 * t580;
t400 = (-t347 * t446 + t575) * t628 + (t264 * t440 - t266 * t318 - t349 * t373) * t626;
t24 = -t340 + t449 * t349 + (t264 * t587 - t266 * t590) * mrSges(7,3) + (-mrSges(5,2) + t435) * t347 + t395 + t400;
t416 = qJD(1) * t626 * t63 - t24 * qJD(2);
t171 = t259 * t260;
t32 = m(7) * (-t129 * t146 + t145 * t441 + t171) + m(6) * (t259 * t423 + t171);
t415 = -t32 * qJD(1) - t561 * t63 / 0.2e1;
t411 = t318 * t674 + t669 * t440;
t390 = (mrSges(7,3) * t593 + t447) * t263 + (mrSges(7,3) * t596 + t448) * t265 + t440 * t616 + t216 * t594 + t409 + t643;
t10 = t390 + t633;
t28 = t460 - t649;
t54 = t373 * t297 - (-t301 / 0.2e1 + t601) * t418 + (t598 + t299 / 0.2e1) * t641;
t408 = -t28 * qJD(1) - t10 * qJD(2) - t54 * qJD(4);
t397 = pkin(5) * t522 * t625 + t259 * t585 + t460;
t14 = t397 + (t472 + t473) * t441 + t678;
t2 = t503 / 0.4e1 + t76 * t473 + t86 * t472 - t433 + t558 * t446 / 0.2e1 + (t167 + t140) * t641 / 0.4e1 + (t360 - t431) * t513 / 0.4e1 + t669 * t471 + t637 * t386 + t642 * mrSges(7,3) + t455 * t572 + t512 * t583 + t512 * t584 + ((-t373 * t512 + t523) * pkin(5) + t411) * t625 + t661 - t680;
t45 = -pkin(4) * t357 - t504 / 0.2e1 + t361 * t582 + (-t431 / 0.2e1 + t360 / 0.2e1) * t386 + t54 + (m(7) * t373 + t298) * t573;
t407 = -t14 * qJD(1) - t2 * qJD(2) - t45 * qJD(4);
t19 = -t485 * mrSges(7,2) + t486 * mrSges(7,1) + (t578 * t617 + t508 / 0.2e1 + (t384 * t611 + t578 * t607) * mrSges(7,3)) * pkin(5);
t354 = (mrSges(7,1) * t384 + mrSges(7,2) * t578) * pkin(5);
t38 = mrSges(7,2) * t451;
t67 = (t596 + t595) * mrSges(7,2) + (t594 + t593) * mrSges(7,1);
t406 = -t38 * qJD(1) + t19 * qJD(2) - t67 * qJD(4) + t354 * qJD(5);
t350 = t354 * qJD(6);
t30 = t392 + t660;
t29 = t460 + t649;
t26 = t396 + t647;
t25 = t264 * t471 - t266 * t472 + t356 * t591 + t395 - t400 + t455 - t347 * t673 / 0.2e1;
t18 = -t569 / 0.2e1 + t214 * t437 + t442 / 0.2e1 - pkin(5) * t508 / 0.2e1 - t432 / 0.2e1 - t570 / 0.2e1 + t566 / 0.2e1 - t565 / 0.2e1 + t501;
t16 = t665 - t679;
t15 = t451 * t531 + t397 - t678;
t13 = t401 + t648;
t9 = t390 - t633;
t8 = t664 - t677;
t5 = -t388 + t294 * t459 - t259 * t502 / 0.2e1 + t389 + (-t282 + t165) * t613 + t646 * t260 / 0.2e1;
t1 = ((t360 / 0.4e1 - t431 / 0.4e1) * t385 + pkin(9) * t435 + (t583 + t584 + (-t298 / 0.2e1 + t373 * t626) * pkin(5)) * t386) * t349 + (pkin(5) * t523 + t411) * t625 + (-t418 * t486 + t485 * t641 + t642) * mrSges(7,3) + (t233 / 0.4e1 + t637) * t386 + t643 + t661 + t680;
t20 = [t34 * qJD(2) + t32 * qJD(4), t30 * qJD(3) + t5 * qJD(4) + t8 * qJD(5) + t13 * qJD(6) + t516 + (t170 * t214 + t169 * t216 + m(7) * (t169 * t76 + t170 * t77) + t280 * t292 + t279 * t294 + m(6) * (t187 * t279 + t188 * t280 - t518) + t322 * t535 + m(5) * (-t308 * t322 - t518) + m(4) * (-pkin(2) * t577 + t387 * t443) * t382 + (m(5) * t370 - mrSges(4,1) * t383 + mrSges(5,1) * t347 + mrSges(4,2) * t381 - mrSges(5,2) * t349 - mrSges(3,1)) * t468 - (-t283 - t532 - t635) * t321 + (-mrSges(3,2) + t666) * t511) * qJD(2), qJD(2) * t30 + t493, t5 * qJD(2) + t15 * qJD(5) + t29 * qJD(6) + ((t145 * t440 - t318 * t146 + t373 * t260) * t625 + (-pkin(4) * t260 - t259 * t446) * t627) * t629 - t415 + ((-t145 * t418 - t146 * t641) * mrSges(7,3) + (t298 + t654) * t260 + (mrSges(5,2) - t673) * t259) * qJD(4), t8 * qJD(2) + t15 * qJD(4) + ((-t129 * t578 + t384 * t441) * t624 - t226 * mrSges(6,2) - t227 * mrSges(6,1) + t39) * qJD(5) + t672, t13 * qJD(2) + t29 * qJD(4) + t39 * qJD(5) + t672; qJD(3) * t31 - qJD(4) * t4 + qJD(5) * t7 + qJD(6) * t12 - t516, qJD(3) * t21 + qJD(4) * t3 + qJD(5) * t6 + qJD(6) * t11 (t264 * t641 - t266 * t418) * t561 + t25 * qJD(4) + t16 * qJD(5) + t26 * qJD(6) + t426, t25 * qJD(3) + t1 * qJD(5) + t9 * qJD(6) + ((t241 * t373 + t318 * t82 + t440 * t81) * t625 + (-pkin(4) * t308 + pkin(9) * t636) * t627) * t629 + t430 + (t230 * t580 + t232 * t582 + t373 * t165 + t139 * t587 + t137 * t589 + Ifges(5,6) * t349 + t440 * t215 + t318 * t213 + t307 * mrSges(5,2) + t241 * t298 + t264 * t601 - t266 * t598 + pkin(4) * t282 + t82 * t531 - t81 * t530 - t293 * t571 + t386 * pkin(9) * t291 + (-Ifges(5,5) + t504 / 0.2e1 + t360 * t581) * t347 + (Ifges(6,5) * t385 + Ifges(7,5) * t418 + Ifges(6,6) * t386 + Ifges(7,6) * t641) * t591 + t654 * t308 + t636 * mrSges(6,3)) * qJD(4), t16 * qJD(3) + t1 * qJD(4) + (t566 + (t384 * t86 + t578 * t85) * t624 + t442 - t432 - t565 + Ifges(6,5) * t513 + Ifges(6,6) * t512 - t187 * mrSges(6,2) - t188 * mrSges(6,1) + t501) * qJD(5) + t18 * qJD(6) + t429, t26 * qJD(3) + t9 * qJD(4) + t18 * qJD(5) + (t501 - t569 - t570) * qJD(6) + t428; -qJD(2) * t31 + t493, qJD(4) * t24 + qJD(5) * t17 + qJD(6) * t27 - t426, 0, -t416, -t662 + t427 + (-t297 - t357 + (-t418 * t491 + t574 * t641) * m(7)) * qJD(5), -qJD(5) * t297 + t521 - t662; t4 * qJD(2) + t14 * qJD(5) + t28 * qJD(6) + t415, -qJD(3) * t24 + qJD(5) * t2 + qJD(6) * t10 - t430, t416, qJD(5) * t45 + qJD(6) * t54 ((-t318 * t578 + t384 * t440) * t624 - t641 * t439 - t418 * t495 + t444 + t356 * pkin(9) + t59) * qJD(5) + t671 - t407, t59 * qJD(5) - t408 + t671; -qJD(2) * t7 - qJD(4) * t14 + qJD(6) * t38, -qJD(3) * t17 - qJD(4) * t2 - qJD(6) * t19 - t429, -t427, t67 * qJD(6) + t407, -t350, -t350 - t406; -t12 * qJD(2) - t28 * qJD(4) - t38 * qJD(5), -qJD(3) * t27 - qJD(4) * t10 + qJD(5) * t19 - t428, -t521, -t67 * qJD(5) + t408, t406, 0;];
Cq  = t20;
