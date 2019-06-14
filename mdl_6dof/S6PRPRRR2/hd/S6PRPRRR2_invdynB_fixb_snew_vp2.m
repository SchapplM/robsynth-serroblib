% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRR2
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 00:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:22:13
% EndTime: 2019-05-05 00:22:24
% DurationCPUTime: 11.31s
% Computational Cost: add. (212496->298), mult. (395135->382), div. (0->0), fcn. (281147->14), ass. (0->130)
t590 = sin(pkin(6));
t597 = sin(qJ(2));
t622 = t590 * t597;
t601 = cos(qJ(2));
t621 = t590 * t601;
t593 = cos(pkin(6));
t620 = t593 * t597;
t619 = t593 * t601;
t589 = sin(pkin(11));
t592 = cos(pkin(11));
t577 = g(1) * t589 - g(2) * t592;
t578 = -g(1) * t592 - g(2) * t589;
t587 = -g(3) + qJDD(1);
t540 = t577 * t619 - t578 * t597 + t587 * t621;
t538 = qJDD(2) * pkin(2) + t540;
t541 = t577 * t620 + t601 * t578 + t587 * t622;
t603 = qJD(2) ^ 2;
t539 = -pkin(2) * t603 + t541;
t588 = sin(pkin(12));
t591 = cos(pkin(12));
t527 = t588 * t538 + t591 * t539;
t525 = -pkin(3) * t603 + qJDD(2) * pkin(8) + t527;
t558 = -t577 * t590 + t593 * t587;
t557 = qJDD(3) + t558;
t596 = sin(qJ(4));
t600 = cos(qJ(4));
t518 = t600 * t525 + t596 * t557;
t573 = (-mrSges(5,1) * t600 + mrSges(5,2) * t596) * qJD(2);
t615 = qJD(2) * qJD(4);
t586 = t596 * t615;
t576 = qJDD(2) * t600 - t586;
t617 = qJD(2) * t596;
t579 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t617;
t574 = (-pkin(4) * t600 - pkin(9) * t596) * qJD(2);
t602 = qJD(4) ^ 2;
t616 = qJD(2) * t600;
t513 = -pkin(4) * t602 + qJDD(4) * pkin(9) + t574 * t616 + t518;
t526 = t591 * t538 - t588 * t539;
t524 = -qJDD(2) * pkin(3) - t603 * pkin(8) - t526;
t613 = t600 * t615;
t575 = qJDD(2) * t596 + t613;
t516 = (-t575 - t613) * pkin(9) + (-t576 + t586) * pkin(4) + t524;
t595 = sin(qJ(5));
t599 = cos(qJ(5));
t508 = -t595 * t513 + t599 * t516;
t571 = qJD(4) * t599 - t595 * t617;
t548 = qJD(5) * t571 + qJDD(4) * t595 + t575 * t599;
t569 = qJDD(5) - t576;
t572 = qJD(4) * t595 + t599 * t617;
t584 = qJD(5) - t616;
t506 = (t571 * t584 - t548) * pkin(10) + (t571 * t572 + t569) * pkin(5) + t508;
t509 = t599 * t513 + t595 * t516;
t547 = -qJD(5) * t572 + qJDD(4) * t599 - t575 * t595;
t556 = pkin(5) * t584 - pkin(10) * t572;
t568 = t571 ^ 2;
t507 = -pkin(5) * t568 + pkin(10) * t547 - t556 * t584 + t509;
t594 = sin(qJ(6));
t598 = cos(qJ(6));
t504 = t506 * t598 - t507 * t594;
t549 = t571 * t598 - t572 * t594;
t521 = qJD(6) * t549 + t547 * t594 + t548 * t598;
t550 = t571 * t594 + t572 * t598;
t532 = -mrSges(7,1) * t549 + mrSges(7,2) * t550;
t583 = qJD(6) + t584;
t536 = -mrSges(7,2) * t583 + mrSges(7,3) * t549;
t565 = qJDD(6) + t569;
t502 = m(7) * t504 + mrSges(7,1) * t565 - mrSges(7,3) * t521 - t532 * t550 + t536 * t583;
t505 = t506 * t594 + t507 * t598;
t520 = -qJD(6) * t550 + t547 * t598 - t548 * t594;
t537 = mrSges(7,1) * t583 - mrSges(7,3) * t550;
t503 = m(7) * t505 - mrSges(7,2) * t565 + mrSges(7,3) * t520 + t532 * t549 - t537 * t583;
t494 = t598 * t502 + t594 * t503;
t551 = -mrSges(6,1) * t571 + mrSges(6,2) * t572;
t554 = -mrSges(6,2) * t584 + mrSges(6,3) * t571;
t492 = m(6) * t508 + mrSges(6,1) * t569 - mrSges(6,3) * t548 - t551 * t572 + t554 * t584 + t494;
t555 = mrSges(6,1) * t584 - mrSges(6,3) * t572;
t608 = -t502 * t594 + t598 * t503;
t493 = m(6) * t509 - mrSges(6,2) * t569 + mrSges(6,3) * t547 + t551 * t571 - t555 * t584 + t608;
t609 = -t492 * t595 + t599 * t493;
t489 = m(5) * t518 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t576 - qJD(4) * t579 + t573 * t616 + t609;
t517 = -t596 * t525 + t557 * t600;
t580 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t616;
t512 = -qJDD(4) * pkin(4) - pkin(9) * t602 + t574 * t617 - t517;
t510 = -pkin(5) * t547 - pkin(10) * t568 + t556 * t572 + t512;
t606 = m(7) * t510 - t520 * mrSges(7,1) + mrSges(7,2) * t521 - t549 * t536 + t537 * t550;
t604 = -m(6) * t512 + t547 * mrSges(6,1) - mrSges(6,2) * t548 + t571 * t554 - t555 * t572 - t606;
t498 = m(5) * t517 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t575 + qJD(4) * t580 - t573 * t617 + t604;
t610 = t600 * t489 - t498 * t596;
t479 = m(4) * t527 - mrSges(4,1) * t603 - qJDD(2) * mrSges(4,2) + t610;
t490 = t492 * t599 + t493 * t595;
t605 = -m(5) * t524 + t576 * mrSges(5,1) - mrSges(5,2) * t575 - t579 * t617 + t580 * t616 - t490;
t486 = m(4) * t526 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t603 + t605;
t475 = t588 * t479 + t591 * t486;
t473 = m(3) * t540 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t603 + t475;
t611 = t591 * t479 - t486 * t588;
t474 = m(3) * t541 - mrSges(3,1) * t603 - qJDD(2) * mrSges(3,2) + t611;
t482 = t596 * t489 + t600 * t498;
t614 = m(4) * t557 + t482;
t481 = m(3) * t558 + t614;
t461 = t473 * t619 + t474 * t620 - t481 * t590;
t459 = m(2) * t577 + t461;
t466 = -t473 * t597 + t601 * t474;
t465 = m(2) * t578 + t466;
t618 = t592 * t459 + t589 * t465;
t460 = t473 * t621 + t474 * t622 + t593 * t481;
t612 = -t459 * t589 + t592 * t465;
t528 = Ifges(7,5) * t550 + Ifges(7,6) * t549 + Ifges(7,3) * t583;
t530 = Ifges(7,1) * t550 + Ifges(7,4) * t549 + Ifges(7,5) * t583;
t495 = -mrSges(7,1) * t510 + mrSges(7,3) * t505 + Ifges(7,4) * t521 + Ifges(7,2) * t520 + Ifges(7,6) * t565 - t528 * t550 + t530 * t583;
t529 = Ifges(7,4) * t550 + Ifges(7,2) * t549 + Ifges(7,6) * t583;
t496 = mrSges(7,2) * t510 - mrSges(7,3) * t504 + Ifges(7,1) * t521 + Ifges(7,4) * t520 + Ifges(7,5) * t565 + t528 * t549 - t529 * t583;
t542 = Ifges(6,5) * t572 + Ifges(6,6) * t571 + Ifges(6,3) * t584;
t544 = Ifges(6,1) * t572 + Ifges(6,4) * t571 + Ifges(6,5) * t584;
t483 = -mrSges(6,1) * t512 + mrSges(6,3) * t509 + Ifges(6,4) * t548 + Ifges(6,2) * t547 + Ifges(6,6) * t569 - pkin(5) * t606 + pkin(10) * t608 + t598 * t495 + t594 * t496 - t572 * t542 + t584 * t544;
t543 = Ifges(6,4) * t572 + Ifges(6,2) * t571 + Ifges(6,6) * t584;
t484 = mrSges(6,2) * t512 - mrSges(6,3) * t508 + Ifges(6,1) * t548 + Ifges(6,4) * t547 + Ifges(6,5) * t569 - pkin(10) * t494 - t495 * t594 + t496 * t598 + t542 * t571 - t543 * t584;
t562 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t596 + Ifges(5,6) * t600) * qJD(2);
t563 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t596 + Ifges(5,2) * t600) * qJD(2);
t467 = mrSges(5,2) * t524 - mrSges(5,3) * t517 + Ifges(5,1) * t575 + Ifges(5,4) * t576 + Ifges(5,5) * qJDD(4) - pkin(9) * t490 - qJD(4) * t563 - t483 * t595 + t484 * t599 + t562 * t616;
t564 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t596 + Ifges(5,4) * t600) * qJD(2);
t476 = Ifges(5,4) * t575 + Ifges(5,2) * t576 + Ifges(5,6) * qJDD(4) - t562 * t617 + qJD(4) * t564 - mrSges(5,1) * t524 + mrSges(5,3) * t518 - Ifges(6,5) * t548 - Ifges(6,6) * t547 - Ifges(6,3) * t569 - t572 * t543 + t571 * t544 - mrSges(6,1) * t508 + mrSges(6,2) * t509 - Ifges(7,5) * t521 - Ifges(7,6) * t520 - Ifges(7,3) * t565 - t550 * t529 + t549 * t530 - mrSges(7,1) * t504 + mrSges(7,2) * t505 - pkin(5) * t494 - pkin(4) * t490;
t457 = mrSges(4,2) * t557 - mrSges(4,3) * t526 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t603 - pkin(8) * t482 + t467 * t600 - t476 * t596;
t462 = Ifges(4,6) * qJDD(2) + t603 * Ifges(4,5) - mrSges(4,1) * t557 + mrSges(4,3) * t527 - Ifges(5,5) * t575 - Ifges(5,6) * t576 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t517 + mrSges(5,2) * t518 - t595 * t484 - t599 * t483 - pkin(4) * t604 - pkin(9) * t609 - pkin(3) * t482 + (-t563 * t596 + t564 * t600) * qJD(2);
t454 = -mrSges(3,1) * t558 + mrSges(3,3) * t541 + t603 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t614 + qJ(3) * t611 + t588 * t457 + t591 * t462;
t455 = mrSges(3,2) * t558 - mrSges(3,3) * t540 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t603 - qJ(3) * t475 + t457 * t591 - t462 * t588;
t607 = pkin(7) * t466 + t454 * t601 + t455 * t597;
t456 = mrSges(3,1) * t540 - mrSges(3,2) * t541 + mrSges(4,1) * t526 - mrSges(4,2) * t527 + t596 * t467 + t600 * t476 + pkin(3) * t605 + pkin(8) * t610 + pkin(2) * t475 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t453 = mrSges(2,2) * t587 - mrSges(2,3) * t577 - t597 * t454 + t601 * t455 + (-t460 * t590 - t461 * t593) * pkin(7);
t452 = -mrSges(2,1) * t587 + mrSges(2,3) * t578 - pkin(1) * t460 - t590 * t456 + t593 * t607;
t1 = [-m(1) * g(1) + t612; -m(1) * g(2) + t618; -m(1) * g(3) + m(2) * t587 + t460; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t618 - t589 * t452 + t592 * t453; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t612 + t592 * t452 + t589 * t453; -mrSges(1,1) * g(2) + mrSges(2,1) * t577 + mrSges(1,2) * g(1) - mrSges(2,2) * t578 + pkin(1) * t461 + t593 * t456 + t590 * t607;];
tauB  = t1;
