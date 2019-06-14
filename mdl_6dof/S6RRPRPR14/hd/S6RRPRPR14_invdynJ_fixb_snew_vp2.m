% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 17:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR14_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:55:13
% EndTime: 2019-05-06 16:55:21
% DurationCPUTime: 4.14s
% Computational Cost: add. (27328->330), mult. (61367->391), div. (0->0), fcn. (43441->10), ass. (0->150)
t620 = Ifges(5,4) + Ifges(6,6);
t640 = -Ifges(5,2) - Ifges(6,3);
t634 = Ifges(5,6) - Ifges(6,5);
t639 = -2 * qJD(3);
t638 = Ifges(3,1) + Ifges(4,2);
t637 = Ifges(5,1) + Ifges(6,2);
t621 = Ifges(3,4) + Ifges(4,6);
t619 = Ifges(3,5) - Ifges(4,4);
t636 = Ifges(5,5) - Ifges(6,4);
t635 = Ifges(3,2) + Ifges(4,3);
t618 = Ifges(3,6) - Ifges(4,5);
t633 = Ifges(3,3) + Ifges(4,1);
t632 = Ifges(5,3) + Ifges(6,1);
t568 = cos(pkin(6));
t562 = qJD(1) * t568 + qJD(2);
t570 = sin(qJ(4));
t574 = cos(qJ(4));
t567 = sin(pkin(6));
t575 = cos(qJ(2));
t602 = qJD(1) * t575;
t596 = t567 * t602;
t528 = t562 * t570 + t574 * t596;
t529 = t562 * t574 - t570 * t596;
t571 = sin(qJ(2));
t603 = qJD(1) * t567;
t597 = t571 * t603;
t555 = qJD(4) + t597;
t631 = t640 * t528 + t620 * t529 + t634 * t555;
t630 = (pkin(2) * t562 + t639) * t597;
t577 = qJD(1) ^ 2;
t572 = sin(qJ(1));
t576 = cos(qJ(1));
t592 = -g(1) * t576 - g(2) * t572;
t600 = qJDD(1) * t567;
t544 = -pkin(1) * t577 + pkin(8) * t600 + t592;
t612 = t567 * t571;
t595 = t572 * g(1) - g(2) * t576;
t625 = pkin(8) * t567;
t543 = qJDD(1) * pkin(1) + t577 * t625 + t595;
t614 = t543 * t568;
t498 = -g(3) * t612 + t575 * t544 + t571 * t614;
t545 = (-pkin(2) * t575 - qJ(3) * t571) * t603;
t560 = t562 ^ 2;
t561 = qJDD(1) * t568 + qJDD(2);
t464 = pkin(2) * t560 - t561 * qJ(3) - t545 * t596 + t562 * t639 - t498;
t550 = -qJD(2) * t597 + t575 * t600;
t495 = qJD(4) * t529 + t574 * t550 + t561 * t570;
t505 = mrSges(6,1) * t528 - mrSges(6,3) * t555;
t629 = -t495 * mrSges(6,2) - t528 * t505;
t507 = -mrSges(5,2) * t555 - mrSges(5,3) * t528;
t608 = t505 - t507;
t622 = mrSges(5,1) - mrSges(6,2);
t628 = t622 * t495 - t608 * t528;
t627 = -2 * qJD(5);
t626 = -pkin(2) - pkin(9);
t624 = g(3) * t568;
t623 = mrSges(3,1) - mrSges(4,2);
t615 = t528 * t555;
t613 = t567 ^ 2 * t577;
t611 = t567 * t575;
t548 = pkin(3) * t597 - pkin(9) * t562;
t549 = (qJD(2) * t602 + qJDD(1) * t571) * t567;
t598 = t575 ^ 2 * t613;
t454 = -pkin(3) * t598 - t624 - qJ(3) * t549 + t626 * t550 + (-t543 + (-qJ(3) * t562 * t575 - t548 * t571) * qJD(1)) * t567 + t630;
t604 = g(3) * t611 + t571 * t544;
t589 = -qJ(3) * t560 + t545 * t597 + qJDD(3) + t604;
t456 = pkin(3) * t549 + t626 * t561 + (-pkin(3) * t562 * t603 - pkin(9) * t571 * t613 - t614) * t575 + t589;
t449 = t574 * t454 + t570 * t456;
t610 = t634 * t528 - t636 * t529 - t632 * t555;
t609 = -t620 * t528 + t637 * t529 + t636 * t555;
t607 = (t571 * t619 + t575 * t618) * t603 + t633 * t562;
t606 = (t571 * t621 + t575 * t635) * t603 + t618 * t562;
t605 = (t571 * t638 + t575 * t621) * t603 + t619 * t562;
t599 = t575 * t614;
t448 = -t570 * t454 + t456 * t574;
t496 = -qJD(4) * t528 - t550 * t570 + t561 * t574;
t500 = mrSges(5,1) * t528 + mrSges(5,2) * t529;
t538 = qJDD(4) + t549;
t499 = pkin(4) * t528 - qJ(5) * t529;
t552 = t555 ^ 2;
t444 = -pkin(4) * t538 - qJ(5) * t552 + t529 * t499 + qJDD(5) - t448;
t438 = (t528 * t529 - t538) * pkin(10) + (t496 + t615) * pkin(5) + t444;
t509 = pkin(5) * t529 - pkin(10) * t555;
t527 = t528 ^ 2;
t453 = pkin(3) * t550 - pkin(9) * t598 + t562 * t548 - t464;
t579 = (-t496 + t615) * qJ(5) + t453 + (t555 * pkin(4) + t627) * t529;
t441 = t579 + (pkin(4) + pkin(10)) * t495 - pkin(5) * t527 - t509 * t529;
t569 = sin(qJ(6));
t573 = cos(qJ(6));
t436 = t438 * t573 - t441 * t569;
t503 = t528 * t573 - t555 * t569;
t462 = qJD(6) * t503 + t495 * t569 + t538 * t573;
t504 = t528 * t569 + t555 * t573;
t471 = -mrSges(7,1) * t503 + mrSges(7,2) * t504;
t526 = qJD(6) + t529;
t476 = -mrSges(7,2) * t526 + mrSges(7,3) * t503;
t492 = qJDD(6) + t496;
t433 = m(7) * t436 + mrSges(7,1) * t492 - mrSges(7,3) * t462 - t471 * t504 + t476 * t526;
t437 = t438 * t569 + t441 * t573;
t461 = -qJD(6) * t504 + t495 * t573 - t538 * t569;
t477 = mrSges(7,1) * t526 - mrSges(7,3) * t504;
t434 = m(7) * t437 - mrSges(7,2) * t492 + mrSges(7,3) * t461 + t471 * t503 - t477 * t526;
t425 = t433 * t573 + t434 * t569;
t501 = -mrSges(6,2) * t528 - mrSges(6,3) * t529;
t586 = -m(6) * t444 - t496 * mrSges(6,1) - t529 * t501 - t425;
t421 = m(5) * t448 - mrSges(5,3) * t496 - t500 * t529 + t538 * t622 - t555 * t608 + t586;
t508 = mrSges(5,1) * t555 - mrSges(5,3) * t529;
t585 = -pkin(4) * t552 + qJ(5) * t538 - t499 * t528 + t449;
t442 = t555 * t627 - t585;
t506 = mrSges(6,1) * t529 + mrSges(6,2) * t555;
t440 = -pkin(5) * t495 - pkin(10) * t527 + ((2 * qJD(5)) + t509) * t555 + t585;
t588 = -m(7) * t440 + mrSges(7,1) * t461 - t462 * mrSges(7,2) + t476 * t503 - t504 * t477;
t581 = -m(6) * t442 + t538 * mrSges(6,3) + t555 * t506 - t588;
t430 = m(5) * t449 - mrSges(5,2) * t538 - t508 * t555 + (-t500 - t501) * t528 + (-mrSges(5,3) - mrSges(6,1)) * t495 + t581;
t594 = -t570 * t421 + t574 * t430;
t593 = -t569 * t433 + t573 * t434;
t517 = -t543 * t567 - t624;
t419 = t421 * t574 + t430 * t570;
t465 = -pkin(2) * t550 + (-t562 * t596 - t549) * qJ(3) + t517 + t630;
t541 = -mrSges(4,1) * t596 - mrSges(4,3) * t562;
t591 = -m(4) * t465 + t549 * mrSges(4,3) - t541 * t596 - t594;
t446 = pkin(4) * t495 + t579;
t590 = -m(6) * t446 + t496 * mrSges(6,3) + t529 * t506 - t593;
t470 = -pkin(2) * t561 + t589 - t599;
t587 = -m(4) * t470 - t549 * mrSges(4,1) - t419;
t583 = -m(5) * t453 - t496 * mrSges(5,2) - t529 * t508 + t590;
t467 = Ifges(7,4) * t504 + Ifges(7,2) * t503 + Ifges(7,6) * t526;
t468 = Ifges(7,1) * t504 + Ifges(7,4) * t503 + Ifges(7,5) * t526;
t582 = mrSges(7,1) * t436 - mrSges(7,2) * t437 + Ifges(7,5) * t462 + Ifges(7,6) * t461 + Ifges(7,3) * t492 + t504 * t467 - t503 * t468;
t542 = mrSges(4,1) * t597 + mrSges(4,2) * t562;
t546 = (mrSges(4,2) * t575 - mrSges(4,3) * t571) * t603;
t580 = -m(4) * t464 + t561 * mrSges(4,3) + t562 * t542 + t546 * t596 - t583;
t423 = mrSges(6,2) * t538 + t505 * t555 - t586;
t466 = Ifges(7,5) * t504 + Ifges(7,6) * t503 + Ifges(7,3) * t526;
t427 = -mrSges(7,1) * t440 + mrSges(7,3) * t437 + Ifges(7,4) * t462 + Ifges(7,2) * t461 + Ifges(7,6) * t492 - t466 * t504 + t468 * t526;
t428 = mrSges(7,2) * t440 - mrSges(7,3) * t436 + Ifges(7,1) * t462 + Ifges(7,4) * t461 + Ifges(7,5) * t492 + t466 * t503 - t467 * t526;
t578 = -mrSges(5,2) * t449 - mrSges(6,3) * t442 - pkin(4) * t423 - pkin(10) * t425 - t569 * t427 + t573 * t428 + t609 * t528 + qJ(5) * (-t501 * t528 + t581) + mrSges(6,2) * t444 + mrSges(5,1) * t448 + t632 * t538 + t631 * t529 + t636 * t496 + (-qJ(5) * mrSges(6,1) - t634) * t495;
t547 = (-mrSges(3,1) * t575 + mrSges(3,2) * t571) * t603;
t540 = -mrSges(3,2) * t562 + mrSges(3,3) * t596;
t539 = mrSges(3,1) * t562 - mrSges(3,3) * t597;
t497 = t599 - t604;
t424 = -t590 + t629;
t420 = t580 - t562 * t539 - t561 * mrSges(3,2) + m(3) * t498 + (mrSges(3,3) + mrSges(4,1)) * t550 + t547 * t596 + t628;
t418 = mrSges(4,2) * t561 + t541 * t562 + t546 * t597 - t587;
t417 = t550 * mrSges(4,2) - t542 * t597 - t591;
t416 = m(3) * t497 - mrSges(3,3) * t549 + (t540 - t541) * t562 + t623 * t561 + (-t546 - t547) * t597 + t587;
t415 = mrSges(6,1) * t444 + mrSges(5,2) * t453 - mrSges(5,3) * t448 - mrSges(6,3) * t446 + pkin(5) * t425 - qJ(5) * t424 - t620 * t495 + t637 * t496 + t610 * t528 + t636 * t538 - t631 * t555 + t582;
t414 = -mrSges(5,1) * t453 - mrSges(6,1) * t442 + mrSges(6,2) * t446 + mrSges(5,3) * t449 - pkin(4) * t424 - pkin(5) * t588 - pkin(10) * t593 - t573 * t427 - t569 * t428 + t640 * t495 + t620 * t496 + t610 * t529 + t634 * t538 + t609 * t555;
t413 = mrSges(3,1) * t497 - mrSges(3,2) * t498 + mrSges(4,2) * t470 - mrSges(4,3) * t464 + t574 * t415 - t570 * t414 - pkin(9) * t419 - pkin(2) * t418 + qJ(3) * (t495 * mrSges(5,1) + t528 * t507 + t580 + t629) + t633 * t561 + (mrSges(4,1) * qJ(3) + t618) * t550 + t619 * t549 + (t571 * t606 - t575 * t605) * t603;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t595 - mrSges(2,2) * t592 + (mrSges(4,1) * t470 + mrSges(3,2) * t517 - mrSges(3,3) * t497 - mrSges(4,3) * t465 + pkin(3) * t419 - qJ(3) * t417 + t638 * t549 + t621 * t550 + t619 * t561 - t606 * t562 + t607 * t596 + t578) * t612 + (-mrSges(3,1) * t517 + mrSges(3,3) * t498 - mrSges(4,1) * t464 + mrSges(4,2) * t465 - t570 * t415 - t574 * t414 - pkin(3) * (t583 - t628) - pkin(9) * t594 - pkin(2) * t417 + t605 * t562 + t618 * t561 + t635 * t550 + t621 * t549 - t607 * t597) * t611 + t568 * t413 + pkin(1) * ((t416 * t575 + t420 * t571) * t568 + (-m(3) * t517 - t549 * mrSges(3,2) + t623 * t550 + (t540 * t575 + (-t539 + t542) * t571) * t603 + t591) * t567) + (-t416 * t571 + t420 * t575) * t625; t413; t418; t578; t423; t582;];
tauJ  = t1;
