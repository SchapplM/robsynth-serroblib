% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:29
% EndTime: 2019-12-31 22:24:47
% DurationCPUTime: 9.23s
% Computational Cost: add. (25156->575), mult. (50479->761), div. (0->0), fcn. (54772->8), ass. (0->355)
t632 = -m(6) / 0.2e1;
t683 = mrSges(6,1) / 0.2e1;
t682 = -mrSges(6,2) / 0.2e1;
t589 = sin(qJ(3));
t590 = sin(qJ(2));
t591 = cos(qJ(3));
t592 = cos(qJ(2));
t347 = t589 * t590 - t591 * t592;
t348 = -t589 * t592 - t591 * t590;
t383 = sin(qJ(5));
t384 = sin(qJ(4));
t385 = cos(qJ(5));
t386 = cos(qJ(4));
t640 = -t383 * t384 + t385 * t386;
t242 = t640 * t348;
t563 = t242 * mrSges(6,3);
t181 = mrSges(6,1) * t347 + t563;
t620 = -t181 / 0.2e1;
t494 = t589 * pkin(2);
t372 = t494 + pkin(8);
t576 = pkin(9) + t372;
t342 = t576 * t384;
t343 = t576 * t386;
t262 = -t342 * t383 + t343 * t385;
t454 = -t385 * t342 - t343 * t383;
t681 = -t262 * mrSges(6,1) - t454 * mrSges(6,2);
t622 = -pkin(9) - pkin(8);
t364 = t622 * t384;
t365 = t622 * t386;
t303 = t364 * t383 - t365 * t385;
t453 = t385 * t364 + t365 * t383;
t680 = -t303 * mrSges(6,1) - t453 * mrSges(6,2);
t434 = t383 * t386 + t385 * t384;
t243 = t434 * t347;
t245 = t640 * t347;
t100 = -Ifges(6,1) * t245 + Ifges(6,4) * t243 - t348 * Ifges(6,5);
t377 = Ifges(5,4) * t386;
t442 = Ifges(5,2) * t384 - t377;
t195 = -Ifges(5,6) * t348 + t442 * t347;
t573 = Ifges(5,4) * t384;
t361 = t386 * Ifges(5,1) - t573;
t197 = -Ifges(5,5) * t348 - t361 * t347;
t568 = Ifges(5,2) * t386;
t358 = t568 + t573;
t574 = Ifges(5,1) * t384;
t360 = t377 + t574;
t441 = Ifges(5,5) * t384 + Ifges(5,6) * t386;
t521 = t347 * t386;
t469 = -t521 / 0.2e1;
t522 = t347 * t384;
t470 = t522 / 0.2e1;
t593 = t386 / 0.2e1;
t595 = t384 / 0.2e1;
t601 = -t348 / 0.2e1;
t334 = Ifges(6,4) * t640;
t280 = Ifges(6,1) * t434 + t334;
t607 = t280 / 0.2e1;
t571 = Ifges(6,4) * t434;
t278 = Ifges(6,2) * t640 + t571;
t608 = t278 / 0.2e1;
t98 = -Ifges(6,4) * t245 + Ifges(6,2) * t243 - Ifges(6,6) * t348;
t405 = t243 * t608 - t245 * t607 + t197 * t595 + t195 * t593 + t358 * t470 + t360 * t469 + Ifges(4,6) * t348 + t640 * t98 / 0.2e1 + t434 * t100 / 0.2e1 - Ifges(4,5) * t347 + (Ifges(6,5) * t434 + Ifges(6,6) * t640 + t441) * t601;
t356 = -mrSges(5,1) * t386 + t384 * mrSges(5,2);
t498 = t592 * pkin(6);
t366 = t592 * pkin(7) + t498;
t495 = t590 * pkin(6);
t419 = -t590 * pkin(7) - t495;
t643 = t591 * t366 + t589 * t419;
t652 = t643 * t356;
t658 = t643 * mrSges(4,1);
t301 = t589 * t366 - t591 * t419;
t660 = t301 * mrSges(4,2);
t276 = -mrSges(6,1) * t640 + mrSges(6,2) * t434;
t651 = -pkin(4) * t522 + t643;
t670 = t651 * t276;
t679 = t405 + t652 + t660 - t658 + t670;
t474 = t386 * t591;
t475 = t384 * t591;
t312 = (-t383 * t474 - t385 * t475) * pkin(2);
t313 = (-t383 * t475 + t385 * t474) * pkin(2);
t676 = -t313 / 0.2e1;
t645 = mrSges(6,2) * t676 + t312 * t683;
t678 = (t312 * t385 + t313 * t383) * pkin(4) * t632 - t645;
t677 = t660 / 0.2e1 + t670 / 0.2e1;
t520 = t348 * t384;
t220 = -pkin(4) * t520 + t301;
t675 = t220 * t651;
t497 = t591 * pkin(2);
t373 = -t497 - pkin(3);
t585 = t386 * pkin(4);
t354 = t373 - t585;
t674 = t354 * t651;
t502 = Ifges(6,5) * t640 - Ifges(6,6) * t434;
t37 = t502 + t681;
t673 = t37 * qJD(5);
t374 = -pkin(3) - t585;
t672 = t374 * t651;
t42 = t502 + t680;
t671 = t42 * qJD(5);
t281 = -t348 * pkin(3) + t347 * pkin(8);
t655 = t384 * t301;
t171 = t386 * t281 + t655;
t446 = -t348 * pkin(4) + pkin(9) * t521;
t102 = t171 + t446;
t654 = t386 * t301;
t172 = t384 * t281 - t654;
t501 = pkin(9) * t522;
t131 = t501 + t172;
t66 = t102 * t385 - t131 * t383;
t67 = t102 * t383 + t131 * t385;
t662 = t66 * t683 + t67 * t682;
t496 = t590 * pkin(2);
t267 = t496 + t281;
t164 = t384 * t267 - t654;
t119 = t501 + t164;
t163 = t386 * t267 + t655;
t93 = t163 + t446;
t62 = -t119 * t383 + t385 * t93;
t63 = t119 * t385 + t383 * t93;
t664 = t62 * t683 + t63 * t682;
t669 = t652 / 0.2e1 - t658 / 0.2e1;
t244 = t434 * t348;
t235 = Ifges(6,4) * t244;
t101 = -Ifges(6,1) * t242 + t347 * Ifges(6,5) + t235;
t127 = -mrSges(6,1) * t243 - mrSges(6,2) * t245;
t128 = -mrSges(6,1) * t244 - mrSges(6,2) * t242;
t375 = -t592 * pkin(2) - pkin(1);
t266 = t347 * pkin(3) + t348 * pkin(8) + t375;
t153 = t386 * t266 - t384 * t643;
t154 = t384 * t266 + t386 * t643;
t178 = mrSges(6,2) * t348 + mrSges(6,3) * t243;
t180 = -mrSges(6,1) * t348 + mrSges(6,3) * t245;
t545 = t384 * mrSges(5,1);
t357 = mrSges(5,2) * t386 + t545;
t264 = t357 * t347;
t265 = t357 * t348;
t270 = mrSges(5,2) * t348 + mrSges(5,3) * t522;
t272 = -t348 * mrSges(5,1) + mrSges(5,3) * t521;
t118 = pkin(9) * t520 + t154;
t538 = t118 * t383;
t519 = t348 * t386;
t117 = pkin(9) * t519 + t153;
t92 = t347 * pkin(4) + t117;
t60 = t385 * t92 - t538;
t537 = t118 * t385;
t61 = t383 * t92 + t537;
t614 = -t245 / 0.2e1;
t615 = t244 / 0.2e1;
t616 = t243 / 0.2e1;
t618 = -t242 / 0.2e1;
t572 = Ifges(6,4) * t242;
t99 = Ifges(6,2) * t244 + Ifges(6,6) * t347 - t572;
t668 = t651 * t128 - t643 * t265 + t100 * t618 + t101 * t614 + t220 * t127 + t153 * t272 + t154 * t270 + t61 * t178 + t60 * t180 - t301 * t264 + t375 * (-mrSges(4,1) * t348 - mrSges(4,2) * t347) + t98 * t615 + t99 * t616;
t631 = m(6) / 0.2e1;
t71 = t117 * t385 - t538;
t666 = (-t60 + t71) * t631 + t620;
t610 = t262 / 0.2e1;
t603 = t303 / 0.2e1;
t665 = pkin(3) * t643;
t70 = -t117 * t383 - t537;
t663 = t61 + t70;
t376 = Ifges(5,5) * t386;
t566 = Ifges(5,6) * t384;
t657 = Ifges(4,4) - t376 / 0.2e1 + t566 / 0.2e1;
t525 = t301 * t643;
t656 = t373 * t643;
t653 = t589 * t301;
t546 = t434 * mrSges(6,3);
t430 = Ifges(6,5) * t614 + Ifges(6,6) * t616;
t381 = t384 ^ 2;
t382 = t386 ^ 2;
t644 = t381 + t382;
t432 = pkin(4) * t276 + t361 / 0.2e1 - t358 / 0.2e1;
t456 = t360 / 0.2e1 - t442 / 0.2e1;
t642 = t432 * t384 + t456 * t386;
t534 = t164 * t386;
t535 = t163 * t384;
t436 = t534 - t535;
t639 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t565 = Ifges(6,3) * t348;
t638 = -t565 / 0.2e1 + t430;
t637 = -t430 - t662;
t636 = -t430 - t664;
t455 = t376 - t566;
t587 = pkin(4) * t128;
t489 = t587 / 0.2e1;
t570 = Ifges(5,5) * t347;
t198 = -t348 * t361 + t570;
t505 = t386 * t198;
t508 = t385 * t180;
t567 = Ifges(5,6) * t347;
t196 = t442 * t348 + t567;
t510 = t384 * t196;
t513 = t383 * t178;
t602 = t347 / 0.4e1;
t629 = pkin(4) / 0.2e1;
t635 = Ifges(5,5) * t469 + Ifges(5,6) * t470 + Ifges(5,3) * t601 + t384 * t489 + t455 * t602 - t510 / 0.4e1 + t505 / 0.4e1 + t638 + (t508 + t513) * t629;
t504 = t386 * t270;
t509 = t384 * t272;
t515 = t374 * t127;
t524 = t303 * t178;
t527 = t453 * t180;
t541 = t63 * t640;
t542 = t62 * t434;
t588 = pkin(3) * t264;
t633 = m(5) / 0.2e1;
t634 = (t436 * pkin(8) - t665) * t633 + (t303 * t63 + t453 * t62 + t672) * t631 + t588 / 0.2e1 + t527 / 0.2e1 + t524 / 0.2e1 + t515 / 0.2e1 + (t504 / 0.2e1 - t509 / 0.2e1) * pkin(8) + (t534 / 0.2e1 - t535 / 0.2e1) * mrSges(5,3) + (-t542 / 0.2e1 + t541 / 0.2e1) * mrSges(6,3) + t677;
t630 = -pkin(3) / 0.2e1;
t628 = -pkin(8) / 0.2e1;
t627 = m(5) * pkin(2);
t626 = -mrSges(5,1) / 0.2e1;
t625 = mrSges(5,2) / 0.2e1;
t562 = t244 * mrSges(6,3);
t179 = -mrSges(6,2) * t347 + t562;
t621 = t179 / 0.2e1;
t617 = t242 / 0.2e1;
t613 = -t454 / 0.2e1;
t612 = t454 / 0.2e1;
t611 = -t262 / 0.2e1;
t609 = -t276 / 0.2e1;
t606 = -t453 / 0.2e1;
t605 = t453 / 0.2e1;
t604 = -t303 / 0.2e1;
t600 = t354 / 0.2e1;
t599 = -t372 / 0.2e1;
t598 = t373 / 0.2e1;
t597 = t374 / 0.2e1;
t596 = -t384 / 0.2e1;
t594 = -t386 / 0.2e1;
t586 = pkin(4) * t384;
t584 = t60 * mrSges(6,2);
t583 = t61 * mrSges(6,1);
t578 = t70 * mrSges(6,1);
t577 = t71 * mrSges(6,2);
t564 = pkin(4) * qJD(4);
t547 = t640 * mrSges(6,3);
t493 = mrSges(5,3) * t520;
t271 = -mrSges(5,2) * t347 + t493;
t273 = t347 * mrSges(5,1) + mrSges(5,3) * t519;
t399 = t197 * t594 + t195 * t595 + Ifges(6,5) * t617 - Ifges(6,6) * t244 / 0.2e1 - t657 * t348 + (-Ifges(5,3) + Ifges(4,1) - Ifges(4,2) - Ifges(6,3)) * t347;
t401 = -t505 / 0.2e1 + t510 / 0.2e1 + t430 + t657 * t347;
t6 = m(5) * (t153 * t163 + t154 * t164 + t525) + m(6) * (t60 * t62 + t61 * t63 + t675) + m(4) * t375 * t496 + (mrSges(4,1) * t496 + t401) * t347 - pkin(1) * (t590 * mrSges(3,1) + t592 * mrSges(3,2)) + (-mrSges(4,2) * t496 + t399) * t348 + t164 * t271 + t163 * t273 + t62 * t181 + t63 * t179 + (Ifges(3,1) - Ifges(3,2)) * t592 * t590 + (-t590 ^ 2 + t592 ^ 2) * Ifges(3,4) + t668;
t543 = t6 * qJD(1);
t8 = t401 * t347 + m(5) * (t153 * t171 + t154 * t172 + t525) + m(6) * (t60 * t66 + t61 * t67 + t675) + t172 * t271 + t171 * t273 + t66 * t181 + t67 * t179 + t399 * t348 + t668;
t540 = t8 * qJD(1);
t129 = Ifges(6,2) * t242 + t235;
t130 = Ifges(6,1) * t244 + t572;
t263 = t356 * t348;
t126 = -mrSges(6,1) * t242 + mrSges(6,2) * t244;
t503 = Ifges(6,5) * t244 + Ifges(6,6) * t242;
t421 = t61 * t563 + t220 * t126 + t347 * t503 / 0.2e1;
t9 = -t71 * t179 - t70 * t181 + t519 * t587 - m(6) * (t60 * t70 + t61 * t71) + t154 * t273 - t301 * t263 + (-t99 / 0.2e1 + t130 / 0.2e1) * t242 + (-t129 / 0.2e1 + t60 * mrSges(6,3) - t101 / 0.2e1) * t244 + (t196 * t594 + t198 * t596 + m(6) * t220 * t585 - t347 * t441 / 0.2e1 - t154 * t386 * mrSges(5,3) + (t358 * t596 + t360 * t593) * t348) * t348 - t421 + (t493 - t271) * t153;
t539 = t9 * qJD(1);
t14 = -t181 * t61 + t99 * t617 + t130 * t618 + (t179 - t562) * t60 + t421 + (t101 + t129) * t615;
t536 = t14 * qJD(1);
t533 = t171 * t384;
t532 = t172 * t386;
t529 = t454 * t180;
t528 = t262 * t178;
t518 = t354 * t127;
t274 = mrSges(6,1) * t434 + mrSges(6,2) * t640;
t517 = t354 * t274;
t516 = t373 * t264;
t514 = t374 * t274;
t512 = t383 * t242;
t507 = t385 * t244;
t500 = t66 * t546;
t499 = t67 * t547;
t492 = mrSges(5,3) * t533;
t491 = mrSges(5,3) * t532;
t490 = t385 * t547;
t484 = -t60 / 0.2e1 + t71 / 0.2e1;
t483 = t70 / 0.2e1 + t61 / 0.2e1;
t482 = t372 * t509;
t481 = t372 * t504;
t480 = t565 / 0.2e1;
t471 = t376 * t602;
t468 = t126 * t600;
t467 = t126 * t597;
t462 = t613 + t612;
t461 = t611 + t610;
t279 = Ifges(6,1) * t640 - t571;
t460 = t278 / 0.4e1 - t279 / 0.4e1;
t277 = -Ifges(6,2) * t434 + t334;
t459 = t280 / 0.4e1 + t277 / 0.4e1;
t458 = t606 + t605;
t457 = t604 + t603;
t450 = -t497 / 0.2e1;
t449 = t497 / 0.2e1;
t445 = mrSges(5,3) * (t382 / 0.2e1 + t381 / 0.2e1);
t444 = t198 / 0.4e1 + t570 / 0.2e1;
t440 = t386 * t450;
t200 = t220 * t586;
t402 = (t101 / 0.4e1 + t129 / 0.4e1) * t640 - (-t130 / 0.4e1 + t99 / 0.4e1) * t434 + t220 * t274 / 0.2e1 + t502 * t602;
t395 = t460 * t242 + t459 * t244 + t301 * t357 / 0.2e1 + t402;
t408 = -t434 * t483 + t484 * t640;
t388 = (t242 * t610 + t244 * t613 + t408) * mrSges(6,3) + (t663 * t454 + t200) * t631 + t179 * t612 + t468 + t263 * t598 + t395 + t666 * t262;
t410 = t358 / 0.4e1 - t361 / 0.4e1 + t573 / 0.2e1 + t568 / 0.4e1;
t413 = (t360 / 0.4e1 - t442 / 0.4e1 + t574 / 0.4e1) * t384;
t396 = t372 * t445 + t413 + ((t354 * t632 + t609) * pkin(4) + t410) * t386;
t415 = -t196 / 0.4e1 - 0.3e1 / 0.4e1 * t567 + t489;
t423 = -t513 / 0.2e1 - t508 / 0.2e1;
t427 = m(6) * (t383 * t63 + t385 * t62);
t429 = t163 * t626 + t164 * t625;
t2 = t388 + (-t427 / 0.2e1 + t423) * pkin(4) + (t396 + t639) * t348 + t471 + (t271 * t599 + t415) * t384 + (t273 * t599 + t444) * t386 + t429 + t636;
t340 = t354 * t586;
t406 = -(-t279 / 0.2e1 + t608) * t434 + (t607 + t277 / 0.2e1) * t640;
t400 = t406 + t642;
t27 = m(6) * t340 + t373 * t357 + t400 + t517;
t439 = t2 * qJD(1) + t27 * qJD(2);
t394 = -t312 * t546 + t313 * t547 + (-mrSges(4,1) + t276 + t356) * t494 + (t644 * mrSges(5,3) - mrSges(4,2)) * t497;
t416 = t644 * t591;
t36 = m(6) * (t262 * t313 + t312 * t454 + t354 * t494) + (t416 * t372 + t589 * t373) * t627 + t394;
t435 = t532 - t533;
t387 = -m(5) * (t656 + t435 * t372 + (-t153 * t475 + t154 * t474 + t653) * pkin(2)) / 0.2e1 + (t220 * t494 + t262 * t67 + t312 * t60 + t313 * t61 + t454 * t66 + t674) * t632 - t529 / 0.2e1 - t528 / 0.2e1 + t312 * t620 + t179 * t676 - t518 / 0.2e1 + t516 / 0.2e1 + t492 / 0.2e1 - t491 / 0.2e1 + t482 / 0.2e1 - t481 / 0.2e1 + t500 / 0.2e1 - t499 / 0.2e1 + t384 * t273 * t449 + t271 * t440 - (t128 - t265) * t494 / 0.2e1 - t669 - t677;
t7 = t387 + t634 + t669;
t438 = -t7 * qJD(1) + t36 * qJD(2);
t392 = (mrSges(6,3) * t613 + t459) * t244 + (mrSges(6,3) * t610 + t460) * t242 + t454 * t621 + t181 * t611 + t468 + t402;
t11 = t480 + t392 + t636;
t30 = t406 + t517;
t437 = t11 * qJD(1) + t30 * qJD(2);
t428 = t171 * t626 + t172 * t625;
t426 = m(6) * (t383 * t67 + t385 * t66);
t422 = t271 * t596 + t273 * t594;
t420 = -t383 * pkin(4) * t546 + t455 + t502;
t355 = t374 * t586;
t403 = (t600 + t597) * t274 + t406;
t390 = (t598 + t630) * t357 + (-(-t457 - t461) * t434 + (t458 + t462) * t640) * mrSges(6,3) + (t340 + t355) * t631 + t403;
t18 = (mrSges(5,2) * t449 + t456) * t386 + (mrSges(5,1) * t449 + t432) * t384 + t390 + t678;
t28 = m(6) * t355 - pkin(3) * t357 + t400 + t514;
t389 = (t242 * t603 + t244 * t606 + t408) * mrSges(6,3) + (t663 * t453 + t200) * t631 + t263 * t630 + t179 * t605 + t467 + t395 + t666 * t303;
t397 = pkin(8) * t445 + t413 + ((t374 * t632 + t609) * pkin(4) + t410) * t386;
t4 = t389 + (-t426 / 0.2e1 + t423) * pkin(4) + (t397 + t639) * t348 + (t273 * t628 + t444) * t386 + (t271 * t628 + t415) * t384 + t471 + t428 + t637;
t414 = t4 * qJD(1) + t18 * qJD(2) + t28 * qJD(3);
t391 = (mrSges(6,3) * t606 + t459) * t244 + (mrSges(6,3) * t603 + t460) * t242 + t453 * t621 + t181 * t604 + t467 + t402;
t13 = t480 + t391 + t637;
t26 = t403 - t645;
t31 = t406 + t514;
t412 = t13 * qJD(1) + t26 * qJD(2) + t31 * qJD(3);
t404 = (t385 * t621 + t383 * t620 + (t512 / 0.2e1 - t507 / 0.2e1) * mrSges(6,3)) * pkin(4);
t16 = -t483 * mrSges(6,1) + t484 * mrSges(6,2) + t404;
t353 = (mrSges(6,1) * t383 + mrSges(6,2) * t385) * pkin(4);
t40 = t461 * mrSges(6,1) + t462 * mrSges(6,2);
t46 = t457 * mrSges(6,1) + t458 * mrSges(6,2);
t407 = -qJD(1) * t16 - qJD(2) * t40 - qJD(3) * t46 + qJD(4) * t353;
t344 = t353 * qJD(5);
t25 = t403 + t645;
t17 = mrSges(5,2) * t440 + t450 * t545 + t390 + t642 - t678;
t15 = -t584 / 0.2e1 - t583 / 0.2e1 - t577 / 0.2e1 + t578 / 0.2e1 + t404 + t503;
t12 = t391 + t638 + t662;
t10 = t392 + t638 + t664;
t5 = t422 * pkin(8) + t397 * t348 + t426 * t629 + t389 - t428 + t635 + t662;
t3 = t396 * t348 + t422 * t372 + t427 * t629 + t388 - t429 + t635 + t664;
t1 = t405 - t387 + (t356 / 0.2e1 - mrSges(4,1) / 0.2e1) * t643 + t634;
t19 = [qJD(2) * t6 + qJD(3) * t8 - qJD(4) * t9 + qJD(5) * t14, t543 + (m(6) * (t262 * t63 + t454 * t62 + t674) + t518 + t481 + m(5) * (t436 * t372 + t656) + m(4) * (-t591 * t643 - t653) * pkin(2) + (t347 * t497 + t348 * t494) * mrSges(4,3) + mrSges(3,2) * t495 - Ifges(3,6) * t590 + Ifges(3,5) * t592 - mrSges(3,1) * t498 + t436 * mrSges(5,3) - t482 + (-t542 + t541) * mrSges(6,3) - t516 + t528 + t529 + t679) * qJD(2) + t1 * qJD(3) + t3 * qJD(4) + t10 * qJD(5), t540 + t1 * qJD(2) + (t491 - t492 + t499 - t500 + t515 + t524 + t527 + t588 + (t504 - t509) * pkin(8) + t679) * qJD(3) + t5 * qJD(4) + t12 * qJD(5) + 0.2e1 * ((pkin(8) * t435 - t665) * t633 + (t303 * t67 + t453 * t66 + t672) * t631) * qJD(3), -t539 + t3 * qJD(2) + t5 * qJD(3) + (-t154 * mrSges(5,1) - t153 * mrSges(5,2) + Ifges(5,5) * t520 + Ifges(5,6) * t519 + t503 - t577 + t578) * qJD(4) + t15 * qJD(5) + (m(6) * (t383 * t71 + t385 * t70) + (-t507 + t512) * mrSges(6,3)) * t564, t536 + t10 * qJD(2) + t12 * qJD(3) + t15 * qJD(4) + (t503 - t583 - t584) * qJD(5); -qJD(3) * t7 + qJD(4) * t2 + qJD(5) * t11 - t543, qJD(3) * t36 + qJD(4) * t27 + qJD(5) * t30, (m(6) * (t303 * t313 + t312 * t453 + t374 * t494) + (-t589 * pkin(3) + t416 * pkin(8)) * t627 + t394) * qJD(3) + t17 * qJD(4) + t25 * qJD(5) + t438, t17 * qJD(3) + (t356 * t372 + t420 + t681) * qJD(4) + t673 + (-t490 + m(6) * (-t262 * t385 + t383 * t454)) * t564 + t439, t25 * qJD(3) + t37 * qJD(4) + t437 + t673; qJD(2) * t7 + qJD(4) * t4 + qJD(5) * t13 - t540, qJD(4) * t18 + qJD(5) * t26 - t438, qJD(4) * t28 + qJD(5) * t31, (t356 * pkin(8) + t420 + t680) * qJD(4) + t671 + (-t490 + m(6) * (-t303 * t385 + t383 * t453)) * t564 + t414, t42 * qJD(4) + t412 + t671; -qJD(2) * t2 - qJD(3) * t4 + qJD(5) * t16 + t539, -qJD(3) * t18 + qJD(5) * t40 - t439, qJD(5) * t46 - t414, -t344, -t344 - t407; -qJD(2) * t11 - qJD(3) * t13 - qJD(4) * t16 - t536, -qJD(3) * t26 - qJD(4) * t40 - t437, -qJD(4) * t46 - t412, t407, 0;];
Cq = t19;
