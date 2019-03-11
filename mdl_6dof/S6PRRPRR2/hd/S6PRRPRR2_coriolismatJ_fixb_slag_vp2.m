% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:56:13
% EndTime: 2019-03-08 21:56:32
% DurationCPUTime: 11.23s
% Computational Cost: add. (26217->630), mult. (57352->886), div. (0->0), fcn. (66578->12), ass. (0->339)
t633 = m(7) * pkin(5);
t683 = -t633 / 0.2e1;
t413 = cos(qJ(3));
t409 = sin(qJ(3));
t541 = cos(pkin(12));
t469 = t541 * t409;
t540 = sin(pkin(12));
t366 = -t413 * t540 - t469;
t601 = -t366 / 0.2e1;
t569 = pkin(5) * qJD(5);
t468 = t540 * t409;
t364 = -t413 * t541 + t468;
t405 = sin(pkin(6));
t414 = cos(qJ(2));
t516 = t405 * t414;
t334 = t364 * t516;
t408 = sin(qJ(5));
t412 = cos(qJ(5));
t410 = sin(qJ(2));
t518 = t405 * t410;
t291 = t334 * t408 + t412 * t518;
t292 = -t334 * t412 + t408 * t518;
t407 = sin(qJ(6));
t411 = cos(qJ(6));
t170 = t291 * t411 - t292 * t407;
t171 = t291 * t407 + t292 * t411;
t631 = mrSges(6,2) / 0.2e1;
t632 = -mrSges(6,1) / 0.2e1;
t679 = -mrSges(7,2) / 0.2e1;
t657 = t170 * mrSges(7,1) / 0.2e1 + t171 * t679;
t682 = t291 * t632 + t292 * t631 + (t170 * t411 + t171 * t407) * t683 - t657;
t452 = t407 * t412 + t411 * t408;
t275 = t452 * t364;
t651 = -t407 * t408 + t411 * t412;
t277 = t651 * t364;
t558 = t277 * mrSges(7,2);
t560 = t275 * mrSges(7,1);
t656 = t560 / 0.2e1 + t558 / 0.2e1;
t681 = (t275 * t411 - t277 * t407) * t683 - t656;
t680 = -mrSges(7,1) / 0.2e1;
t394 = -pkin(3) * t413 - pkin(2);
t303 = pkin(4) * t364 + pkin(9) * t366 + t394;
t582 = -qJ(4) - pkin(8);
t378 = t582 * t413;
t655 = -t541 * t378 + t582 * t468;
t197 = t412 * t303 - t408 * t655;
t522 = t366 * t412;
t158 = pkin(10) * t522 + t197;
t132 = pkin(5) * t364 + t158;
t198 = t303 * t408 + t412 * t655;
t523 = t366 * t408;
t159 = pkin(10) * t523 + t198;
t539 = t159 * t407;
t82 = t132 * t411 - t539;
t90 = t158 * t411 - t539;
t678 = -t82 + t90;
t406 = cos(pkin(6));
t354 = t406 * t413 - t409 * t518;
t517 = t405 * t413;
t355 = t406 * t409 + t410 * t517;
t438 = t354 * t540 + t355 * t541;
t229 = -t408 * t438 - t412 * t516;
t230 = -t408 * t516 + t412 * t438;
t129 = t229 * t407 + t230 * t411;
t466 = t411 * t229 - t230 * t407;
t39 = -t129 * mrSges(7,1) - t466 * mrSges(7,2);
t677 = t39 * qJD(6);
t483 = t540 * pkin(3);
t390 = t483 + pkin(9);
t583 = pkin(10) + t390;
t358 = t583 * t408;
t359 = t583 * t412;
t301 = -t358 * t407 + t359 * t411;
t465 = -t411 * t358 - t359 * t407;
t361 = Ifges(7,6) * t452;
t362 = Ifges(7,5) * t651;
t503 = t362 - t361;
t57 = -t301 * mrSges(7,1) - t465 * mrSges(7,2) + t503;
t676 = t57 * qJD(6);
t543 = t412 * mrSges(6,2);
t546 = t408 * mrSges(6,1);
t377 = t543 + t546;
t294 = t377 * t364;
t307 = -mrSges(6,2) * t364 + mrSges(6,3) * t523;
t507 = t412 * t307;
t309 = t364 * mrSges(6,1) + mrSges(6,3) * t522;
t511 = t408 * t309;
t444 = t507 / 0.2e1 - t511 / 0.2e1;
t675 = -t294 / 0.2e1 - t444;
t314 = -mrSges(7,1) * t651 + mrSges(7,2) * t452;
t484 = t541 * pkin(3);
t391 = -t484 - pkin(4);
t375 = -t412 * pkin(5) + t391;
t674 = m(7) * t375 + t314;
t615 = -t277 / 0.2e1;
t617 = t275 / 0.2e1;
t448 = Ifges(7,5) * t615 + Ifges(7,6) * t617;
t593 = pkin(3) * t409;
t304 = -pkin(4) * t366 + pkin(9) * t364 + t593;
t326 = -t378 * t540 - t582 * t469;
t199 = t412 * t304 + t326 * t408;
t524 = t364 * t412;
t133 = -pkin(5) * t366 + pkin(10) * t524 + t199;
t200 = t408 * t304 - t326 * t412;
t525 = t364 * t408;
t161 = pkin(10) * t525 + t200;
t85 = t133 * t411 - t161 * t407;
t86 = t133 * t407 + t161 * t411;
t646 = t86 * mrSges(7,2) / 0.2e1 + t85 * t680 - t448;
t643 = Ifges(7,3) * t601 - t646;
t271 = -t541 * t354 + t355 * t540;
t673 = t271 / 0.2e1;
t607 = t301 / 0.2e1;
t672 = m(5) * t394;
t538 = t159 * t411;
t83 = t132 * t407 + t538;
t89 = -t158 * t407 - t538;
t670 = t83 + t89;
t553 = t364 * mrSges(5,3);
t172 = t271 * t438;
t401 = t408 ^ 2;
t403 = t412 ^ 2;
t667 = -t403 - t401;
t313 = mrSges(7,1) * t452 + mrSges(7,2) * t651;
t665 = t313 * qJD(6);
t663 = m(4) * t414;
t548 = t452 * mrSges(7,3);
t396 = Ifges(6,5) * t412;
t571 = Ifges(6,6) * t408;
t660 = Ifges(5,4) - t396 / 0.2e1 + t571 / 0.2e1;
t658 = t390 * t667;
t397 = Ifges(6,4) * t412;
t382 = Ifges(6,1) * t408 + t397;
t654 = t409 ^ 2 + t413 ^ 2;
t379 = t409 * mrSges(4,1) + t413 * mrSges(4,2);
t462 = Ifges(6,2) * t408 - t397;
t649 = -mrSges(4,1) * t413 + mrSges(4,2) * t409;
t297 = t366 * t382;
t597 = -t390 / 0.2e1;
t648 = t307 * t597 + t297 / 0.4e1;
t150 = t452 * t271;
t151 = t651 * t271;
t447 = t150 * t680 + t151 * t679;
t634 = m(5) * pkin(3);
t645 = t540 * t634 - mrSges(5,2);
t276 = t452 * t366;
t616 = -t276 / 0.2e1;
t274 = t651 * t366;
t618 = t274 / 0.2e1;
t559 = t276 * mrSges(7,3);
t213 = -mrSges(7,2) * t364 + t559;
t626 = t213 / 0.2e1;
t644 = t466 * t626 + (t129 * t618 + t466 * t616) * mrSges(7,3);
t167 = -mrSges(7,1) * t276 - mrSges(7,2) * t274;
t212 = mrSges(7,2) * t366 + t275 * mrSges(7,3);
t214 = -mrSges(7,1) * t366 + t277 * mrSges(7,3);
t247 = -pkin(5) * t525 + t655;
t248 = -pkin(5) * t523 + t326;
t295 = t377 * t366;
t306 = mrSges(6,2) * t366 + mrSges(6,3) * t525;
t308 = -mrSges(6,1) * t366 + mrSges(6,3) * t524;
t455 = t197 * t408 - t198 * t412;
t623 = t229 / 0.2e1;
t215 = mrSges(7,1) * t364 + t274 * mrSges(7,3);
t624 = t215 / 0.2e1;
t629 = t129 / 0.2e1;
t636 = m(7) / 0.2e1;
t638 = m(6) / 0.2e1;
t639 = m(5) / 0.2e1;
t642 = -t516 * t593 * t639 + (t199 * t229 + t200 * t230 + t326 * t438 + (t455 + t655) * t271) * t638 + (t129 * t86 + t150 * t82 - t151 * t83 + t247 * t271 + t248 * t438 + t466 * t85) * t636 + t466 * t214 / 0.2e1 + t212 * t629 + t150 * t624 - t151 * t626 + t308 * t623 + t230 * t306 / 0.2e1 + (-t295 / 0.2e1 + t167 / 0.2e1) * t438;
t376 = -mrSges(6,1) * t412 + mrSges(6,2) * t408;
t641 = m(6) * t391 - t541 * t634 - mrSges(5,1) + t376;
t640 = 0.2e1 * m(7);
t637 = -m(7) / 0.2e1;
t630 = -mrSges(6,3) / 0.2e1;
t165 = -mrSges(7,1) * t274 + mrSges(7,2) * t276;
t628 = t165 / 0.2e1;
t625 = -t215 / 0.2e1;
t293 = t366 * t376;
t614 = t293 / 0.2e1;
t610 = -t465 / 0.2e1;
t609 = t465 / 0.2e1;
t608 = -t301 / 0.2e1;
t606 = t313 / 0.2e1;
t573 = Ifges(7,4) * t452;
t316 = Ifges(7,2) * t651 + t573;
t605 = t316 / 0.2e1;
t363 = Ifges(7,4) * t651;
t318 = Ifges(7,1) * t452 + t363;
t604 = t318 / 0.2e1;
t602 = t364 / 0.4e1;
t600 = t651 / 0.2e1;
t599 = t452 / 0.2e1;
t598 = t377 / 0.2e1;
t596 = t408 / 0.2e1;
t595 = -t412 / 0.2e1;
t594 = t412 / 0.2e1;
t592 = pkin(5) * t408;
t590 = t82 * mrSges(7,2);
t589 = t83 * mrSges(7,1);
t586 = t89 * mrSges(7,1);
t585 = t90 * mrSges(7,2);
t580 = m(7) * qJD(4);
t575 = Ifges(6,4) * t408;
t574 = Ifges(7,4) * t274;
t552 = t364 * Ifges(6,5);
t551 = t364 * Ifges(6,6);
t550 = t366 * mrSges(5,3);
t549 = t651 * mrSges(7,3);
t547 = t407 * t86;
t544 = t411 * t85;
t535 = t229 * t408;
t534 = t230 * t412;
t533 = t248 * t408;
t432 = (t274 * t599 + t616 * t651) * mrSges(7,3) + t213 * t600 - t452 * t624;
t26 = t432 - t656;
t532 = t26 * qJD(2);
t333 = t366 * t516;
t202 = t271 * t333;
t531 = t291 * t408;
t530 = t292 * t412;
t528 = t326 * t333;
t527 = t326 * t366;
t519 = t405 ^ 2 * t410;
t35 = m(7) * (t129 * t171 + t170 * t466 - t202) + m(6) * (t229 * t291 + t230 * t292 - t202) + (-t354 * t405 * t409 + t355 * t517 - t519) * t663 + (-t334 * t438 - t414 * t519 - t202) * m(5);
t526 = t35 * qJD(1);
t228 = t366 * t271;
t521 = t390 * t408;
t520 = t390 * t412;
t515 = t407 * t212;
t514 = t407 * t274;
t234 = t366 * t462 + t551;
t512 = t408 * t234;
t510 = t411 * t214;
t509 = t411 * t276;
t504 = Ifges(7,5) * t276 + Ifges(7,6) * t274;
t501 = t634 / 0.2e1;
t67 = t150 * t651 - t151 * t452;
t500 = t67 * qJD(3) * t636;
t497 = -t82 / 0.2e1 + t90 / 0.2e1;
t496 = t83 / 0.2e1 + t89 / 0.2e1;
t487 = t549 / 0.2e1;
t486 = -t548 / 0.2e1;
t485 = t543 / 0.2e1;
t482 = t165 * t673;
t481 = t271 * t606;
t480 = t525 / 0.2e1;
t479 = -t524 / 0.2e1;
t477 = -t516 / 0.2e1;
t476 = t167 * t596;
t472 = t629 - t129 / 0.2e1;
t315 = -Ifges(7,2) * t452 + t363;
t471 = t315 / 0.4e1 + t318 / 0.4e1;
t317 = Ifges(7,1) * t651 - t573;
t470 = t316 / 0.4e1 - t317 / 0.4e1;
t310 = -t366 * mrSges(5,1) - t364 * mrSges(5,2);
t467 = t396 - t571;
t464 = -t167 + t295 + t550;
t463 = t503 * t602;
t383 = Ifges(6,1) * t412 - t575;
t380 = Ifges(6,2) * t412 + t575;
t166 = -t558 - t560;
t416 = (-t333 * t391 + (t530 - t531) * t390) * t638 + (t170 * t465 + t171 * t301 - t333 * t375) * t636 + t333 * mrSges(5,1) / 0.2e1 + t334 * mrSges(5,2) / 0.2e1 + (t333 * t541 - t334 * t540) * t501 + t170 * t486 + t171 * t487 + t531 * t630 + mrSges(6,3) * t530 / 0.2e1 - (t314 + t376) * t333 / 0.2e1 + t379 * t477;
t4 = t166 * t673 - t416 + (t310 + t379) * t477 + t642 + t675 * t271;
t137 = -Ifges(7,4) * t277 + Ifges(7,2) * t275 - Ifges(7,6) * t366;
t138 = Ifges(7,2) * t276 + t364 * Ifges(7,6) - t574;
t139 = -Ifges(7,1) * t277 + Ifges(7,4) * t275 - Ifges(7,5) * t366;
t265 = Ifges(7,4) * t276;
t140 = -Ifges(7,1) * t274 + t364 * Ifges(7,5) + t265;
t233 = -Ifges(6,6) * t366 + t364 * t462;
t235 = -t366 * Ifges(6,5) - t364 * t383;
t236 = -t366 * t383 + t552;
t311 = mrSges(5,1) * t364 - mrSges(5,2) * t366;
t5 = (-Ifges(4,4) * t409 + pkin(3) * t311) * t409 + (Ifges(4,4) * t413 + (Ifges(4,1) - Ifges(4,2)) * t409) * t413 + t394 * t310 - pkin(2) * t379 - t326 * t294 + t198 * t306 + t200 * t307 + t197 * t308 + t199 * t309 - t274 * t139 / 0.2e1 + t276 * t137 / 0.2e1 + t247 * t167 + t248 * t166 + t86 * t213 + t82 * t214 + t85 * t215 + t83 * t212 + m(7) * (t247 * t248 + t82 * t85 + t83 * t86) - t655 * t295 + m(6) * (t197 * t199 + t198 * t200 + t655 * t326) + (t236 * t595 + t512 / 0.2e1 + t448 + t660 * t364) * t364 + (Ifges(7,5) * t618 + Ifges(7,6) * t616 + t235 * t595 + t233 * t596 - t660 * t366 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2) - Ifges(7,3)) * t364) * t366 + t593 * t672 + t140 * t615 + t138 * t617;
t461 = t4 * qJD(1) + t5 * qJD(2);
t296 = t366 * t380;
t168 = Ifges(7,2) * t274 + t265;
t169 = Ifges(7,1) * t276 + t574;
t427 = t248 * t165 + (t168 / 0.2e1 + t140 / 0.2e1) * t276 + (t138 / 0.2e1 - t169 / 0.2e1 + t83 * mrSges(7,3)) * t274 - t82 * t559 + t364 * t504 / 0.2e1;
t6 = m(7) * (t82 * t89 + t83 * t90) + t90 * t213 + t89 * t215 + t326 * t293 + t197 * t307 - t198 * t309 + ((t552 / 0.2e1 - t197 * mrSges(6,3) + t236 / 0.2e1 + t296 / 0.2e1) * t408 + (t551 / 0.2e1 + t198 * mrSges(6,3) - t297 / 0.2e1 + t234 / 0.2e1 + (-m(7) * t248 - t167) * pkin(5)) * t412) * t366 + t427;
t419 = (t534 / 0.2e1 - t535 / 0.2e1) * t366 * mrSges(6,3) + (-pkin(5) * t271 * t522 + t670 * t466) * t636 + t307 * t623 - t230 * t309 / 0.2e1 + t644 + (t636 * t678 - t624) * t129;
t7 = (t628 + t614) * t271 + t419 + t682;
t460 = t7 * qJD(1) + t6 * qJD(2);
t29 = t481 - t447;
t11 = t82 * t213 - t83 * t215 + t427;
t429 = t129 * t625 + t482 + t644;
t12 = t429 - t657;
t459 = t12 * qJD(1) + t11 * qJD(2);
t422 = (t452 * t678 + t670 * t651) * t636 + t432 + t444;
t445 = t485 + t546 / 0.2e1;
t17 = -t445 * t364 + t422 + t681;
t458 = t17 * qJD(2);
t22 = -t277 * t213 + t275 * t215 + t464 * t366 + (-t507 + t511 + t553) * t364 + m(7) * (-t248 * t366 + t275 * t82 - t277 * t83) + m(6) * (t364 * t455 - t527) + m(5) * (-t364 * t655 - t527);
t454 = -t534 + t535;
t424 = (-t129 * t277 + t275 * t466 - t228) * t636 + (t364 * t454 - t228) * t638 + (-t364 * t438 - t228) * t639;
t428 = (t291 * t412 + t292 * t408) * t638 + (t170 * t651 + t171 * t452) * t636 + t518 * t639;
t33 = t424 - t428;
t457 = t33 * qJD(1) + t22 * qJD(2);
t446 = t199 * t632 + t200 * t631;
t418 = (t364 * t658 - t366 * t391) * t638 + (t275 * t465 - t277 * t301 - t366 * t375) * t636 + t314 * t601 - t293 / 0.2e1 + (-t364 * t540 + t366 * t541) * t501 + t275 * t486 - t277 * t487 - t667 * t364 * t630;
t421 = (t199 * t412 + t200 * t408) * t638 + (t452 * t86 + t651 * t85) * t636 + t214 * t600 + t212 * t599 + t306 * t596 + t308 * t594 + t409 * t501;
t21 = t310 - t418 + t421;
t443 = qJD(1) * t637 * t67 - t21 * qJD(2);
t30 = m(7) * (-t129 * t151 + t150 * t466 + t172) + m(6) * (t271 * t454 + t172);
t442 = -t30 * qJD(1) - t580 * t67 / 0.2e1;
t441 = (t150 * t411 - t151 * t407) * t633;
t437 = t301 * t678 + t670 * t465;
t430 = (t140 / 0.4e1 + t168 / 0.4e1) * t651 - (-t169 / 0.4e1 + t138 / 0.4e1) * t452 + t248 * t606 + t375 * t628;
t417 = t471 * t276 + t470 * t274 + (t296 / 0.4e1 + t236 / 0.4e1 + t309 * t597) * t412 + (t274 * t607 + t276 * t610 - t452 * t496 + t497 * t651) * mrSges(7,3) - t301 * t624 + t213 * t609 + t326 * t598 + t391 * t614 + t430;
t426 = (t382 / 0.4e1 - t462 / 0.4e1) * t408 + (t401 / 0.2e1 + t403 / 0.2e1) * t390 * mrSges(6,3) + (-t383 / 0.4e1 + t380 / 0.4e1 + (t375 * t637 - t314 / 0.2e1) * pkin(5)) * t412;
t1 = t417 + (-0.3e1 / 0.4e1 * t571 + 0.3e1 / 0.4e1 * t396 + t362 / 0.4e1 - t361 / 0.4e1) * t364 + (-t234 / 0.4e1 + t648) * t408 + t437 * t636 + (-t510 / 0.2e1 + t476 - t515 / 0.2e1 + (-t547 / 0.4e1 - t544 / 0.4e1 + t533 / 0.4e1) * t640) * pkin(5) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t426) * t366 + t446 + t646;
t423 = t271 * t592 * t636 - t472 * t548;
t14 = -t441 / 0.2e1 + t423 + t447 + (-t445 + t606 + t598) * t271;
t51 = t375 * t313 - (t605 - t317 / 0.2e1) * t452 + (t315 / 0.2e1 + t604) * t651;
t40 = -t408 * t380 / 0.2e1 + t391 * t377 + t383 * t596 + (-t462 / 0.2e1 + t382 / 0.2e1) * t412 + t51 + t674 * t592;
t436 = -t14 * qJD(1) - t1 * qJD(2) - t40 * qJD(3);
t420 = (mrSges(7,3) * t610 + t471) * t276 + (mrSges(7,3) * t607 + t470) * t274 + t465 * t626 + t215 * t608 + t463 + t430;
t10 = t420 - t643;
t28 = t481 + t447;
t435 = -t28 * qJD(1) - t10 * qJD(2) - t51 * qJD(3);
t431 = (t411 * t626 + t407 * t625 + (t514 / 0.2e1 - t509 / 0.2e1) * mrSges(7,3)) * pkin(5);
t19 = -mrSges(7,1) * t496 + mrSges(7,2) * t497 + t431;
t374 = (mrSges(7,1) * t407 + mrSges(7,2) * t411) * pkin(5);
t38 = t472 * mrSges(7,1);
t60 = (t610 + t609) * mrSges(7,2) + (t608 + t607) * mrSges(7,1);
t433 = t38 * qJD(1) - t19 * qJD(2) - t60 * qJD(3) + t374 * qJD(5);
t367 = t374 * qJD(6);
t32 = t424 + t428;
t27 = t432 + t656;
t23 = t418 + t421;
t18 = -t590 / 0.2e1 - t589 / 0.2e1 - t585 / 0.2e1 + t586 / 0.2e1 + t431 + t504;
t16 = mrSges(6,1) * t480 + t364 * t485 + t422 - t681;
t15 = t271 * t598 + t441 / 0.2e1 + t423 + t29 + t377 * t673;
t13 = t429 + t657;
t9 = t420 + t643;
t8 = t271 * t614 + t419 + t482 - t682;
t3 = t416 + t673 * t553 + (-t379 / 0.2e1 - t310 / 0.2e1) * t516 + t642 + (-t553 / 0.2e1 + t166 / 0.2e1 + t675) * t271;
t2 = t417 + t426 * t366 + t648 * t408 - t512 / 0.4e1 + (t544 + t547) * t633 / 0.2e1 - t446 + (pkin(5) * t533 + t437) * t636 + t463 + t643 + (t510 + t515) * pkin(5) / 0.2e1 + pkin(5) * t476 + Ifges(6,5) * t479 + Ifges(6,6) * t480 + Ifges(6,3) * t601 + t467 * t602;
t20 = [t35 * qJD(2) + t30 * qJD(3), t3 * qJD(3) + t32 * qJD(4) + t8 * qJD(5) + t13 * qJD(6) + t526 + (t170 * t215 + t171 * t213 + t291 * t309 + t292 * t307 + t334 * t553 + t464 * t333 + 0.2e1 * (t170 * t82 + t171 * t83 - t248 * t333) * t636 + 0.2e1 * (t197 * t291 + t198 * t292 - t528) * t638 + 0.2e1 * (-t334 * t655 - t528) * t639 + ((mrSges(4,3) * t654 - mrSges(3,2)) * t414 + t654 * pkin(8) * t663 + (-m(4) * pkin(2) - mrSges(3,1) + t311 + t649 + t672) * t410) * t405) * qJD(2), t3 * qJD(2) + t15 * qJD(5) + t29 * qJD(6) - t442 + (-t151 * t549 - t150 * t548 + m(7) * (t150 * t465 - t151 * t301) - t354 * mrSges(4,2) - t355 * mrSges(4,1) + (t641 + t674) * t438 + (m(6) * t658 + mrSges(6,3) * t667 - t645) * t271) * qJD(3), qJD(2) * t32 + t500, t8 * qJD(2) + t15 * qJD(3) + (-t230 * mrSges(6,1) - t229 * mrSges(6,2) + (-t129 * t411 + t407 * t466) * t633 + t39) * qJD(5) + t677, t13 * qJD(2) + t29 * qJD(3) + t39 * qJD(5) + t677; qJD(3) * t4 + qJD(4) * t33 + qJD(5) * t7 + qJD(6) * t12 - t526, qJD(3) * t5 + qJD(4) * t22 + qJD(5) * t6 + qJD(6) * t11 (t649 * pkin(8) + Ifges(4,5) * t413 - Ifges(4,6) * t409 - t391 * t294 + t375 * t166 - Ifges(5,5) * t364 + Ifges(5,6) * t366 + t247 * t314 + t301 * t212 - t308 * t521 + t641 * t655 - t645 * t326 + (m(6) * t390 + mrSges(6,3)) * (-t199 * t408 + t200 * t412) + (Ifges(6,5) * t408 + Ifges(7,5) * t452 + Ifges(6,6) * t412 + Ifges(7,6) * t651) * t601 + t465 * t214 + m(7) * (t247 * t375 + t301 * t86 + t465 * t85) + t382 * t479 + t380 * t480 + t306 * t520 - t85 * t548 + t86 * t549 + t483 * t550 + t484 * t553 + t233 * t594 + t235 * t596 + t139 * t599 + t137 * t600 - t277 * t604 + t275 * t605) * qJD(3) + t23 * qJD(4) + t2 * qJD(5) + t9 * qJD(6) + t461, t23 * qJD(3) + (t275 * t651 - t277 * t452) * t580 + t16 * qJD(5) + t27 * qJD(6) + t457, t2 * qJD(3) + t16 * qJD(4) + (-t198 * mrSges(6,1) - t197 * mrSges(6,2) + Ifges(6,5) * t523 + Ifges(6,6) * t522 + t504 - t585 + t586) * qJD(5) + t18 * qJD(6) + (m(7) * (t407 * t90 + t411 * t89) + (-t509 + t514) * mrSges(7,3)) * t569 + t460, t9 * qJD(3) + t27 * qJD(4) + t18 * qJD(5) + (t504 - t589 - t590) * qJD(6) + t459; -t4 * qJD(2) + t14 * qJD(5) + t28 * qJD(6) + t442, -qJD(4) * t21 + qJD(5) * t1 + qJD(6) * t10 - t461, qJD(5) * t40 + qJD(6) * t51, t443 (-mrSges(6,1) * t520 + mrSges(6,2) * t521 + t467 + t57) * qJD(5) + t676 + (m(7) * (-t301 * t411 + t407 * t465) + (-t407 * t452 - t411 * t651) * mrSges(7,3)) * t569 - t436, t57 * qJD(5) - t435 + t676; -qJD(2) * t33 + t500, qJD(3) * t21 + qJD(5) * t17 + qJD(6) * t26 - t457, -t443, 0 (-t313 - t377) * qJD(5) - t665 + (t407 * t651 - t411 * t452) * t640 * t569 / 0.2e1 + t458, -qJD(5) * t313 + t532 - t665; -qJD(2) * t7 - qJD(3) * t14 - qJD(6) * t38, -qJD(3) * t1 - qJD(4) * t17 + qJD(6) * t19 - t460, t60 * qJD(6) + t436, -t458, -t367, -t367 - t433; -t12 * qJD(2) - t28 * qJD(3) + t38 * qJD(5), -qJD(3) * t10 - qJD(4) * t26 - qJD(5) * t19 - t459, -t60 * qJD(5) + t435, -t532, t433, 0;];
Cq  = t20;
