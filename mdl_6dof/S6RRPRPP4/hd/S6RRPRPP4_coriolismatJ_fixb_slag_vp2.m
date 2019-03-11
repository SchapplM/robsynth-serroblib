% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:56
% EndTime: 2019-03-09 09:59:13
% DurationCPUTime: 9.78s
% Computational Cost: add. (13361->602), mult. (26275->774), div. (0->0), fcn. (25287->6), ass. (0->319)
t679 = Ifges(6,1) + Ifges(7,1);
t678 = Ifges(6,5) + Ifges(7,4);
t386 = sin(pkin(9));
t387 = sin(qJ(4));
t547 = cos(pkin(9));
t594 = cos(qJ(4));
t470 = t547 * t594;
t340 = t386 * t387 - t470;
t491 = t547 * t387;
t341 = t386 * t594 + t491;
t222 = t341 * mrSges(7,1) + t340 * mrSges(7,3);
t376 = t387 * pkin(4) + qJ(3);
t193 = pkin(5) * t341 + qJ(6) * t340 + t376;
t593 = m(7) * t193;
t676 = t222 + t593;
t223 = t341 * mrSges(6,1) - t340 * mrSges(6,2);
t675 = m(6) * t376;
t680 = t223 + t675;
t667 = Ifges(6,6) - Ifges(7,6);
t677 = m(4) * pkin(7) + mrSges(4,1);
t600 = t387 / 0.2e1;
t522 = t594 / 0.2e1;
t674 = mrSges(6,3) + mrSges(7,2);
t325 = Ifges(6,4) * t341;
t579 = Ifges(7,5) * t341;
t673 = -t679 * t340 - t325 + t579;
t388 = sin(qJ(2));
t534 = t387 * t388;
t298 = t386 * t534 - t388 * t470;
t512 = t388 * t594;
t300 = -t386 * t512 - t388 * t491;
t389 = cos(qJ(2));
t672 = t678 * t389 - t679 * t300 + (-Ifges(6,4) + Ifges(7,5)) * t298;
t545 = qJ(3) * t389;
t626 = pkin(2) + pkin(8);
t339 = t388 * t626 - t545;
t384 = t389 * pkin(7);
t355 = t389 * pkin(3) + t384;
t345 = t594 * t355;
t230 = -t339 * t387 + t345;
t233 = t594 * t339 + t387 * t355;
t432 = t594 * t230 + t387 * t233;
t546 = qJ(3) * t388;
t457 = -pkin(2) * t389 - t546;
t170 = pkin(4) * t389 + t345 + (-qJ(5) * t388 - t339) * t387;
t510 = t594 * qJ(5);
t190 = t388 * t510 + t233;
t88 = t386 * t170 + t547 * t190;
t79 = qJ(6) * t389 + t88;
t87 = t170 * t547 - t386 * t190;
t80 = -t389 * pkin(5) - t87;
t671 = t340 * t80 + t341 * t79;
t670 = t340 * t87 - t341 * t88;
t669 = mrSges(6,1) + mrSges(7,1);
t436 = Ifges(5,5) * t387 + Ifges(5,6) * t594;
t666 = t436 * t388;
t552 = t388 * mrSges(7,3);
t533 = t387 * t389;
t299 = t386 * t533 - t389 * t470;
t567 = t299 * mrSges(7,2);
t243 = t552 + t567;
t566 = t299 * mrSges(6,3);
t246 = -mrSges(6,2) * t388 + t566;
t665 = t243 + t246;
t301 = t341 * t389;
t559 = t301 * mrSges(6,3);
t249 = mrSges(6,1) * t388 + t559;
t560 = t301 * mrSges(7,2);
t250 = -mrSges(7,1) * t388 - t560;
t664 = t249 - t250;
t290 = Ifges(6,4) * t299;
t580 = Ifges(7,5) * t299;
t663 = -t679 * t301 + t678 * t388 + t290 - t580;
t662 = t667 * t340 - t678 * t341;
t443 = pkin(1) - t457;
t464 = t389 * mrSges(4,2) - mrSges(4,3) * t388;
t661 = m(4) * t443 - t464;
t590 = pkin(4) * t386;
t372 = qJ(6) + t590;
t505 = t547 * pkin(4);
t375 = -t505 - pkin(5);
t658 = -t340 * t372 + t341 * t375;
t657 = -pkin(5) * t299 + qJ(6) * t301;
t656 = (Ifges(3,4) + Ifges(4,6)) * t388;
t655 = -t250 / 0.2e1 + t249 / 0.2e1;
t654 = -t246 / 0.2e1 - t243 / 0.2e1;
t631 = m(7) / 0.2e1;
t633 = m(6) / 0.2e1;
t653 = t633 + t631;
t652 = -t340 ^ 2 - t341 ^ 2;
t324 = Ifges(7,5) * t340;
t224 = Ifges(7,3) * t341 - t324;
t583 = Ifges(6,4) * t340;
t651 = -t679 * t341 + t224 - t324 + t583;
t650 = Ifges(6,2) * t340 - t325 + t673;
t649 = Ifges(6,2) * t301 + t290 + t663;
t287 = Ifges(7,5) * t301;
t439 = Ifges(7,6) * t388 - Ifges(7,3) * t299 - t287;
t584 = Ifges(6,4) * t301;
t648 = t679 * t299 - t287 + t439 + t584;
t511 = t389 * t594;
t647 = -Ifges(5,5) * t511 + Ifges(5,6) * t533 + t678 * t299 + t667 * t301;
t516 = Ifges(5,4) * t594;
t352 = -Ifges(5,2) * t387 + t516;
t585 = Ifges(5,4) * t387;
t353 = Ifges(5,1) * t594 - t585;
t646 = t352 * t522 + t353 * t600;
t342 = m(7) * t372 + mrSges(7,3);
t506 = t594 * t626;
t423 = -t506 - t510;
t381 = t387 * t626;
t482 = qJ(5) * t387 + t381;
t232 = t386 * t482 + t547 * t423;
t645 = t386 * t423 - t482 * t547;
t520 = mrSges(5,3) * t533;
t347 = t388 * mrSges(5,1) + t520;
t476 = mrSges(5,3) * t511;
t349 = -t388 * mrSges(5,2) - t476;
t508 = t594 * t349;
t596 = -t389 / 0.2e1;
t644 = -t508 / 0.2e1 + t347 * t600 + (t387 ^ 2 + t594 ^ 2) * mrSges(5,3) * t596;
t602 = t341 / 0.2e1;
t605 = t340 / 0.2e1;
t609 = -t301 / 0.2e1;
t614 = -t299 / 0.2e1;
t328 = -t389 * t626 - pkin(1) - t546;
t625 = pkin(3) + pkin(7);
t354 = t625 * t388;
t221 = t328 * t594 + t387 * t354;
t189 = -t389 * t510 + t221;
t166 = t547 * t189;
t220 = -t328 * t387 + t594 * t354;
t188 = qJ(5) * t533 + t220;
t89 = t188 * t386 + t166;
t537 = t386 * t189;
t90 = t188 * t547 - t537;
t643 = -t232 * t614 + t90 * t602 + t89 * t605 + t609 * t645;
t165 = pkin(4) * t388 + t188;
t86 = t386 * t165 + t166;
t77 = t388 * qJ(6) + t86;
t641 = m(6) * t86 + m(7) * t77 + t665;
t85 = t165 * t547 - t537;
t78 = -t388 * pkin(5) - t85;
t640 = -m(6) * t85 + m(7) * t78 - t664;
t629 = m(6) * pkin(4);
t639 = m(7) * t375 - t547 * t629 - t669;
t638 = t386 * t629 - mrSges(6,2) + t342;
t634 = -m(6) / 0.2e1;
t632 = -m(7) / 0.2e1;
t628 = mrSges(6,1) / 0.2e1;
t627 = -t80 / 0.2e1;
t624 = m(7) * t89;
t244 = -t298 * mrSges(7,2) + mrSges(7,3) * t389;
t622 = t244 / 0.2e1;
t245 = -mrSges(6,2) * t389 - t298 * mrSges(6,3);
t621 = t245 / 0.2e1;
t551 = t389 * mrSges(7,1);
t563 = t300 * mrSges(7,2);
t248 = -t551 - t563;
t619 = t248 / 0.2e1;
t616 = -t298 / 0.2e1;
t615 = t298 / 0.2e1;
t612 = t299 / 0.2e1;
t611 = -t300 / 0.2e1;
t610 = t300 / 0.2e1;
t608 = t301 / 0.2e1;
t606 = -t340 / 0.2e1;
t603 = -t341 / 0.2e1;
t597 = t388 / 0.2e1;
t595 = t389 / 0.2e1;
t592 = m(7) * t341;
t582 = Ifges(7,4) * t300;
t581 = Ifges(6,5) * t300;
t578 = Ifges(7,2) * t389;
t575 = Ifges(6,6) * t298;
t574 = Ifges(7,6) * t298;
t573 = Ifges(5,3) * t389;
t572 = Ifges(6,3) * t389;
t571 = t298 * mrSges(6,1);
t570 = t298 * mrSges(7,1);
t569 = t299 * mrSges(6,1);
t568 = t299 * mrSges(7,1);
t521 = t594 * pkin(4);
t316 = (-t521 - t625) * t388;
t128 = t298 * pkin(5) + t300 * qJ(6) + t316;
t317 = pkin(4) * t511 + t355;
t129 = t317 + t657;
t558 = t301 * mrSges(7,3);
t171 = t558 - t568;
t247 = mrSges(6,1) * t389 + t300 * mrSges(6,3);
t346 = mrSges(5,1) * t389 - mrSges(5,3) * t534;
t549 = t389 * mrSges(5,2);
t348 = mrSges(5,3) * t512 - t549;
t437 = Ifges(5,2) * t594 + t585;
t418 = t437 * t389;
t407 = Ifges(5,6) * t388 - t418;
t403 = t594 * t407;
t438 = Ifges(5,1) * t387 + t516;
t419 = t438 * t389;
t409 = Ifges(5,5) * t388 - t419;
t405 = t387 * t409;
t406 = Ifges(5,6) * t389 + t388 * t437;
t408 = Ifges(5,5) * t389 + t388 * t438;
t435 = -mrSges(5,1) * t594 + t387 * mrSges(5,2);
t421 = t388 * t435;
t425 = -Ifges(7,5) * t300 + Ifges(7,6) * t389 + Ifges(7,3) * t298;
t426 = Ifges(6,2) * t299 + Ifges(6,6) * t388 - t584;
t427 = -Ifges(6,4) * t300 - Ifges(6,2) * t298 + Ifges(6,6) * t389;
t562 = t300 * mrSges(7,3);
t465 = t562 + t570;
t561 = t301 * mrSges(6,2);
t466 = t561 + t569;
t564 = t300 * mrSges(6,2);
t467 = -t564 + t571;
t472 = t511 / 0.2e1;
t494 = t533 / 0.2e1;
t3 = -(t666 + t403 + t405 + t572 + t573 + t574 - t575 + t578 - t581 - t582 + (Ifges(3,1) - Ifges(4,3)) * t389 - t656) * t388 / 0.2e1 + t661 * (pkin(2) * t388 - t545) + t663 * t610 + ((Ifges(3,2) - Ifges(4,2)) * t389 + t656) * t597 + (-0.2e1 * Ifges(4,6) * t389 + (-Ifges(4,2) + Ifges(4,3)) * t388) * t595 - t129 * t465 + t316 * t466 - t317 * t467 - m(5) * (t220 * t230 + t221 * t233 - t354 * t355) - m(6) * (t316 * t317 + t85 * t87 + t86 * t88) - m(7) * (t128 * t129 + t77 * t79 + t78 * t80) + t408 * t494 + t443 * (-mrSges(4,2) * t388 - mrSges(4,3) * t389) - t355 * t421 + t406 * t472 + t426 * t615 + t439 * t616 + t425 * t612 + t427 * t614 + pkin(1) * (mrSges(3,1) * t388 + mrSges(3,2) * t389) - t230 * t347 - t221 * t348 - t233 * t349 - t220 * t346 - t86 * t245 - t88 * t246 - t85 * t247 - t78 * t248 - t87 * t249 - t80 * t250 - t79 * t243 - t77 * t244 - t128 * t171 + (-t678 * t301 + t667 * t299 + (Ifges(3,1) - Ifges(3,2) + Ifges(7,2) + Ifges(5,3) + Ifges(6,3)) * t388 + (-t436 + 0.2e1 * Ifges(3,4)) * t389) * t596 + t672 * t608 - t354 * t435 * t389;
t565 = t3 * qJD(1);
t555 = t341 * mrSges(7,2);
t524 = pkin(4) * t533;
t154 = -pkin(5) * t301 - qJ(6) * t299 - t524;
t434 = t387 * mrSges(5,1) + mrSges(5,2) * t594;
t420 = t389 * t434;
t433 = pkin(4) * t466;
t459 = -Ifges(7,3) * t301 + t580;
t473 = -t511 / 0.2e1;
t488 = -t301 * mrSges(6,1) + t299 * mrSges(6,2);
t490 = -t301 * mrSges(7,1) - t299 * mrSges(7,3);
t4 = t154 * t171 - t355 * t420 + t407 * t494 + t409 * t473 + t426 * t608 + t433 * t533 + t459 * t614 + t86 * t559 + t77 * t560 - t85 * t566 + t78 * t567 + (-m(6) * t524 + t488) * t317 + (m(7) * t154 + t490) * t129 + t641 * t90 + t640 * t89 + t646 * t389 ^ 2 + (t520 - t347) * t221 + (t476 + t349) * t220 + t649 * t612 + t648 * t609 + t647 * t597;
t548 = t4 * qJD(1);
t27 = t388 * t243 + m(7) * (t129 * t301 + t388 * t77) + t301 * t171;
t543 = qJD(1) * t27;
t12 = t347 * t534 + t640 * t300 + t641 * t298 + (-t508 + m(5) * (t220 * t387 - t221 * t594) + t661) * t388;
t542 = t12 * qJD(1);
t541 = t300 * t340;
t217 = t341 * t298;
t538 = t341 * t388;
t498 = t538 / 0.2e1;
t468 = t498 + t611;
t169 = t468 * m(7);
t528 = t169 * qJD(1);
t527 = mrSges(6,3) * t590;
t526 = t629 / 0.2e1;
t525 = t634 + t632;
t523 = -t594 / 0.2e1;
t519 = mrSges(7,1) / 0.2e1 + t628;
t518 = mrSges(6,3) / 0.2e1 + mrSges(7,2) / 0.2e1;
t517 = -mrSges(7,3) / 0.2e1 + mrSges(6,2) / 0.2e1;
t496 = t534 / 0.2e1;
t495 = -t533 / 0.2e1;
t493 = t621 + t622;
t492 = -t247 / 0.2e1 + t619;
t489 = -t340 * mrSges(7,1) + t341 * mrSges(7,3);
t487 = -t340 * mrSges(6,1) - t341 * mrSges(6,2);
t480 = t232 * t301 + t299 * t645;
t477 = mrSges(6,3) * t505;
t474 = t512 / 0.2e1;
t469 = 0.2e1 * (m(6) / 0.4e1 + m(7) / 0.4e1) * t652;
t458 = -Ifges(7,3) * t340 - t579;
t455 = -t232 * t89 + t645 * t90;
t404 = -t222 - t223 - t434;
t45 = mrSges(4,3) + (m(5) + m(4)) * qJ(3) + t675 + t593 - t404;
t414 = -t232 * t300 + t298 * t645 + t317;
t398 = -m(5) * t355 / 0.2e1 + t414 * t634 + (t414 + t657) * t632 + mrSges(5,1) * t473;
t399 = m(5) * t432 / 0.2e1 - t670 * t633 + t671 * t631 + t346 * t522;
t9 = (t348 / 0.2e1 + t549 / 0.2e1) * t387 + t517 * t301 + t519 * t299 + (t298 * t518 + t493) * t341 + (t300 * t518 + t492) * t340 + t398 + t399;
t454 = -qJD(1) * t9 + qJD(2) * t45;
t396 = (-t299 * t518 + t654) * t341 + (t301 * t518 + t655) * t340 + (t340 * t85 - t341 * t86 + t480) * t633 + (-t340 * t78 - t341 * t77 + t480) * t631;
t442 = t128 * t632 + t316 * t634;
t11 = -t298 * t519 + t300 * t517 + t396 + t442;
t26 = (m(6) + m(7)) * (-t232 * t340 + t341 * t645) + t674 * t652;
t452 = qJD(1) * t11 - qJD(2) * t26;
t402 = (t129 * t340 + t193 * t301 + t388 * t645) * t631 + t222 * t608 + t171 * t605;
t441 = m(7) * t627 + t551 / 0.2e1;
t21 = -mrSges(7,2) * t468 + t402 + t441;
t71 = t676 * t340;
t451 = qJD(1) * t21 + qJD(2) * t71;
t401 = (t299 * t372 - t301 * t375) * t631 + (t299 * t386 + t301 * t547) * t526;
t424 = t154 * t631 + t495 * t629;
t28 = -t401 + t424 + t488 + t490;
t400 = (-t340 * t375 - t341 * t372) * t631 + (t340 * t547 - t341 * t386) * t526;
t198 = -t340 * pkin(5) + t341 * qJ(6) + t521;
t417 = t198 * t631 + t521 * t633;
t34 = -t400 + t417 + t487 + t489;
t450 = qJD(1) * t28 + qJD(2) * t34;
t43 = t653 * (t299 * t341 - t301 * t340);
t47 = t469 + t525;
t449 = qJD(1) * t43 + qJD(2) * t47;
t14 = t664 * t301 + t665 * t299 + m(7) * (t299 * t77 - t301 * t78) + m(6) * (t299 * t86 + t301 * t85);
t448 = qJD(1) * t14 + qJD(3) * t43;
t157 = m(7) * t301;
t195 = m(7) * t340;
t445 = qJD(1) * t157 + qJD(2) * t195;
t225 = -Ifges(6,2) * t341 - t583;
t393 = -t662 * t388 / 0.4e1 + (t632 * t77 + t634 * t86 + t654) * t232 + (t603 * t85 + t606 * t86 + t643) * mrSges(6,3) + (t606 * t77 + t643) * mrSges(7,2) + (t129 * t198 + t154 * t193 + t455) * t632 + ((t317 * t594 - t376 * t533) * pkin(4) + t455) * t634 + (t632 * t78 - t634 * t85 + t655) * t645 - t644 * t626 + t403 / 0.4e1 - t317 * t487 / 0.2e1 - t376 * t488 / 0.2e1 - t129 * t489 / 0.2e1 - t193 * t490 / 0.2e1 + pkin(4) * t223 * t494 + t405 / 0.4e1 + t352 * t495 + t353 * t472 + t355 * t435 / 0.2e1 + qJ(3) * t420 / 0.2e1 - t387 * t419 / 0.4e1 + (t649 / 0.4e1 - t459 / 0.4e1) * t341 + (t648 / 0.4e1 - t426 / 0.4e1) * t340 + (t651 / 0.4e1 - t225 / 0.4e1) * t301 + (-t650 / 0.4e1 + t458 / 0.4e1) * t299 + t666 / 0.4e1 - t154 * t222 / 0.2e1 - t198 * t171 / 0.2e1 - t433 * t523 + t78 * t555 / 0.2e1 - t594 * t418 / 0.4e1;
t394 = (t372 * t79 + t375 * t80) * t631 - t582 / 0.2e1 - t581 / 0.2e1 + t578 / 0.2e1 - t575 / 0.2e1 + t574 / 0.2e1 + t573 / 0.2e1 + t572 / 0.2e1 + t230 * mrSges(5,1) / 0.2e1 - t233 * mrSges(5,2) / 0.2e1 + t372 * t622 + t375 * t619 + t79 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t627 + t87 * t628 - t88 * mrSges(6,2) / 0.2e1 + (t386 * t88 + t547 * t87) * t526 + Ifges(5,5) * t496 + Ifges(5,6) * t474 + t590 * t621 + t247 * t505 / 0.2e1;
t1 = t394 + t393;
t7 = -qJ(3) * t435 + t193 * t489 + t676 * t198 + t225 * t605 + t376 * t487 + t437 * t600 + t438 * t523 + t458 * t602 + t680 * t521 + t650 * t603 + t651 * t606 - t646;
t422 = -t1 * qJD(1) + t7 * qJD(2);
t395 = ((-t85 + t90) * t341 + (-t86 + t89) * t340) * t634 + ((t78 + t90) * t341 + (-t77 + t89) * t340) * t632 + t249 * t602 + t250 * t603 + t665 * t605 + t644 + t674 * (t299 * t606 + t301 * t603);
t397 = (t298 * t372 + t300 * t375) * t631 + mrSges(6,2) * t616 + mrSges(7,3) * t615 + (t298 * t386 - t300 * t547) * t526 + mrSges(5,1) * t496 + mrSges(5,2) * t474 + t669 * t611;
t5 = t395 + t397;
t416 = t5 * qJD(1);
t410 = t552 + (t372 * t388 + t77) * t631;
t30 = -t624 / 0.2e1 + t410;
t415 = qJD(1) * t30 + qJD(4) * t342;
t168 = (t498 + t610) * m(7);
t64 = 0.2e1 * t645 * t631 - t555;
t46 = t469 - t525;
t44 = t400 + t417;
t42 = t43 * qJD(5);
t38 = t401 + t424;
t29 = t567 + t624 / 0.2e1 + t410;
t20 = -mrSges(7,2) * t538 / 0.2e1 - t563 / 0.2e1 + t402 - t441;
t10 = t562 / 0.2e1 + t570 / 0.2e1 - t564 / 0.2e1 + t571 / 0.2e1 + t396 - t442;
t8 = -t398 + t677 * t389 + mrSges(5,2) * t495 + t348 * t600 + t492 * t340 - t561 / 0.2e1 + t558 / 0.2e1 - t568 / 0.2e1 - t569 / 0.2e1 + t493 * t341 + t399 + t674 * (-t217 / 0.2e1 - t541 / 0.2e1);
t6 = -t395 + t397;
t2 = t394 - t393;
t13 = [-qJD(2) * t3 + qJD(3) * t12 + qJD(4) * t4 + qJD(5) * t14 + qJD(6) * t27, t8 * qJD(3) + t2 * qJD(4) + t10 * qJD(5) + t20 * qJD(6) - t565 + (t670 * mrSges(6,3) - t671 * mrSges(7,2) + t680 * t316 + (m(6) * t88 + m(7) * t79 + t244 + t245) * t645 + m(5) * (-qJ(3) * t354 - t432 * t626) + t464 * pkin(7) + (pkin(7) * mrSges(3,2) + Ifges(4,5) - Ifges(3,6)) * t388 + t193 * t465 + t376 * t467 + t353 * t496 + t352 * t474 - t354 * t434 - (-m(6) * t87 + m(7) * t80 - t247 + t248) * t232 + t224 * t615 + t225 * t616 - t387 * t406 / 0.2e1 + (-Ifges(4,4) + Ifges(3,5)) * t389 + t425 * t602 + t427 * t603 + (Ifges(5,5) * t594 - Ifges(5,6) * t387 - t340 * t678 - t667 * t341) * t595 - t432 * mrSges(5,3) + t672 * t606 + t673 * t611 + t676 * t128 + t677 * t457 + qJ(3) * t421 - t346 * t506 + t408 * t522 - t348 * t381 - mrSges(3,1) * t384) * qJD(2), t8 * qJD(2) + t6 * qJD(4) + t168 * qJD(6) + t42 + t542 + 0.2e1 * t653 * qJD(3) * (t217 + t541) t548 + t2 * qJD(2) + t6 * qJD(3) + (-t221 * mrSges(5,1) - t220 * mrSges(5,2) - t299 * t477 + t301 * t527 + t372 * t560 + t375 * t567 + t638 * t90 + t639 * t89 + t647) * qJD(4) + t38 * qJD(5) + t29 * qJD(6), qJD(2) * t10 + qJD(4) * t38 + t448, qJD(2) * t20 + qJD(3) * t168 + qJD(4) * t29 + t543; -qJD(3) * t9 - qJD(4) * t1 + qJD(5) * t11 + qJD(6) * t21 + t565, qJD(3) * t45 + qJD(4) * t7 - qJD(5) * t26 + qJD(6) * t71, qJD(5) * t46 + t454 (mrSges(5,1) * t381 + mrSges(5,2) * t506 - mrSges(7,2) * t658 + t232 * t638 + t340 * t527 + t341 * t477 + t639 * t645 - t436 + t662) * qJD(4) + t44 * qJD(5) + t64 * qJD(6) + t422, qJD(3) * t46 + qJD(4) * t44 + t452, qJD(4) * t64 + t451; qJD(2) * t9 - qJD(4) * t5 + qJD(6) * t169 + t42 - t542, qJD(5) * t47 - t454, 0 ((-t340 * t386 - t341 * t547) * t629 + m(7) * t658 + t404) * qJD(4) + qJD(6) * t592 - t416, t449, qJD(4) * t592 + t528; qJD(2) * t1 + qJD(3) * t5 - qJD(5) * t28 + qJD(6) * t30 - t548, -qJD(5) * t34 - t422, t416, t342 * qJD(6), -t450, t415; -qJD(2) * t11 + qJD(4) * t28 + qJD(6) * t157 - t448, -qJD(3) * t47 + qJD(4) * t34 + qJD(6) * t195 - t452, -t449, t450, 0, t445; -qJD(2) * t21 - qJD(3) * t169 - qJD(4) * t30 - qJD(5) * t157 - t543, -qJD(5) * t195 - t451, -t528, -t415, -t445, 0;];
Cq  = t13;
