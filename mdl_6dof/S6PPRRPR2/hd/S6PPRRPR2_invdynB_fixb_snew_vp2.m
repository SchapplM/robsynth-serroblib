% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:16:06
% EndTime: 2019-05-04 20:16:15
% DurationCPUTime: 7.81s
% Computational Cost: add. (132595->280), mult. (241910->351), div. (0->0), fcn. (174617->14), ass. (0->134)
t658 = Ifges(5,1) + Ifges(6,2);
t650 = Ifges(5,4) + Ifges(6,6);
t649 = Ifges(5,5) - Ifges(6,4);
t657 = Ifges(5,2) + Ifges(6,3);
t648 = Ifges(5,6) - Ifges(6,5);
t656 = Ifges(5,3) + Ifges(6,1);
t602 = sin(pkin(11));
t606 = cos(pkin(11));
t588 = -g(1) * t606 - g(2) * t602;
t601 = sin(pkin(12));
t605 = cos(pkin(12));
t587 = g(1) * t602 - g(2) * t606;
t600 = -g(3) + qJDD(1);
t604 = sin(pkin(6));
t608 = cos(pkin(6));
t628 = t587 * t608 + t600 * t604;
t551 = -t601 * t588 + t628 * t605;
t552 = t605 * t588 + t628 * t601;
t564 = -t587 * t604 + t600 * t608 + qJDD(2);
t614 = cos(qJ(3));
t607 = cos(pkin(7));
t611 = sin(qJ(3));
t643 = t607 * t611;
t603 = sin(pkin(7));
t644 = t603 * t611;
t545 = t551 * t643 + t614 * t552 + t564 * t644;
t616 = qJD(3) ^ 2;
t543 = -pkin(3) * t616 + qJDD(3) * pkin(9) + t545;
t610 = sin(qJ(4));
t540 = t610 * t543;
t547 = -t551 * t603 + t564 * t607;
t613 = cos(qJ(4));
t642 = t613 * t547;
t537 = -t540 + t642;
t582 = (mrSges(6,2) * t613 - mrSges(6,3) * t610) * qJD(3);
t583 = (-mrSges(5,1) * t613 + mrSges(5,2) * t610) * qJD(3);
t634 = qJD(3) * qJD(4);
t633 = t613 * t634;
t584 = qJDD(3) * t610 + t633;
t635 = qJD(3) * t613;
t590 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t635;
t591 = -mrSges(6,1) * t635 - qJD(4) * mrSges(6,3);
t581 = (-pkin(4) * t613 - qJ(5) * t610) * qJD(3);
t615 = qJD(4) ^ 2;
t636 = qJD(3) * t610;
t627 = -t615 * qJ(5) + t581 * t636 + qJDD(5) + t540;
t652 = pkin(10) * t616;
t653 = -pkin(4) - pkin(10);
t533 = t584 * pkin(5) + t653 * qJDD(4) + (-pkin(5) * t634 - t610 * t652 - t547) * t613 + t627;
t632 = t610 * t634;
t585 = qJDD(3) * t613 - t632;
t593 = pkin(5) * t636 - qJD(4) * pkin(10);
t599 = t613 ^ 2;
t544 = -t611 * t552 + (t551 * t607 + t564 * t603) * t614;
t621 = -qJDD(3) * pkin(3) - t544;
t654 = -2 * qJD(5);
t618 = pkin(4) * t632 + t636 * t654 + (-t584 - t633) * qJ(5) + t621;
t536 = -t593 * t636 + (-pkin(5) * t599 - pkin(9)) * t616 + t653 * t585 + t618;
t609 = sin(qJ(6));
t612 = cos(qJ(6));
t529 = t533 * t612 - t536 * t609;
t579 = -qJD(4) * t609 - t612 * t635;
t559 = qJD(6) * t579 + qJDD(4) * t612 - t585 * t609;
t580 = qJD(4) * t612 - t609 * t635;
t560 = -mrSges(7,1) * t579 + mrSges(7,2) * t580;
t595 = qJD(6) + t636;
t562 = -mrSges(7,2) * t595 + mrSges(7,3) * t579;
t578 = qJDD(6) + t584;
t527 = m(7) * t529 + mrSges(7,1) * t578 - mrSges(7,3) * t559 - t560 * t580 + t562 * t595;
t530 = t533 * t609 + t536 * t612;
t558 = -qJD(6) * t580 - qJDD(4) * t609 - t585 * t612;
t563 = mrSges(7,1) * t595 - mrSges(7,3) * t580;
t528 = m(7) * t530 - mrSges(7,2) * t578 + mrSges(7,3) * t558 + t560 * t579 - t563 * t595;
t519 = t612 * t527 + t609 * t528;
t535 = -qJDD(4) * pkin(4) + t627 - t642;
t622 = -m(6) * t535 - t584 * mrSges(6,1) - t519;
t517 = m(5) * t537 - t584 * mrSges(5,3) + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t590 - t591) * qJD(4) + (-t582 - t583) * t636 + t622;
t538 = t613 * t543 + t610 * t547;
t589 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t636;
t620 = -t615 * pkin(4) + qJDD(4) * qJ(5) + t581 * t635 + t538;
t534 = qJD(4) * t654 - t620;
t592 = mrSges(6,1) * t636 + qJD(4) * mrSges(6,2);
t532 = -t599 * t652 + t585 * pkin(5) + ((2 * qJD(5)) + t593) * qJD(4) + t620;
t623 = -m(7) * t532 + t558 * mrSges(7,1) - t559 * mrSges(7,2) + t579 * t562 - t580 * t563;
t619 = -m(6) * t534 + qJDD(4) * mrSges(6,3) + qJD(4) * t592 + t582 * t635 - t623;
t522 = t583 * t635 + m(5) * t538 - qJDD(4) * mrSges(5,2) - qJD(4) * t589 + (mrSges(5,3) + mrSges(6,1)) * t585 + t619;
t630 = -t517 * t610 + t613 * t522;
t509 = m(4) * t545 - mrSges(4,1) * t616 - qJDD(3) * mrSges(4,2) + t630;
t512 = t613 * t517 + t610 * t522;
t511 = m(4) * t547 + t512;
t651 = t616 * pkin(9);
t542 = t621 - t651;
t539 = -t585 * pkin(4) + t618 - t651;
t640 = -t609 * t527 + t612 * t528;
t626 = -m(6) * t539 - t585 * mrSges(6,2) + t592 * t636 - t640;
t617 = -m(5) * t542 + t590 * t635 + t585 * mrSges(5,1) + (-mrSges(5,2) + mrSges(6,3)) * t584 + (-t589 * t610 - t591 * t613) * qJD(3) + t626;
t515 = m(4) * t544 + qJDD(3) * mrSges(4,1) - t616 * mrSges(4,2) + t617;
t645 = t515 * t614;
t498 = t509 * t643 - t511 * t603 + t607 * t645;
t494 = m(3) * t551 + t498;
t505 = t614 * t509 - t515 * t611;
t504 = m(3) * t552 + t505;
t655 = t494 * t605 + t504 * t601;
t497 = t509 * t644 + t607 * t511 + t603 * t645;
t496 = m(3) * t564 + t497;
t484 = -t496 * t604 + t608 * t655;
t482 = m(2) * t587 + t484;
t490 = -t494 * t601 + t605 * t504;
t489 = m(2) * t588 + t490;
t641 = t606 * t482 + t602 * t489;
t639 = t656 * qJD(4) + (t610 * t649 + t613 * t648) * qJD(3);
t638 = -t648 * qJD(4) + (-t610 * t650 - t613 * t657) * qJD(3);
t637 = t649 * qJD(4) + (t610 * t658 + t650 * t613) * qJD(3);
t483 = t608 * t496 + t604 * t655;
t631 = -t482 * t602 + t606 * t489;
t518 = -t584 * mrSges(6,3) + t591 * t635 - t626;
t553 = Ifges(7,5) * t580 + Ifges(7,6) * t579 + Ifges(7,3) * t595;
t555 = Ifges(7,1) * t580 + Ifges(7,4) * t579 + Ifges(7,5) * t595;
t523 = -mrSges(7,1) * t532 + mrSges(7,3) * t530 + Ifges(7,4) * t559 + Ifges(7,2) * t558 + Ifges(7,6) * t578 - t553 * t580 + t555 * t595;
t554 = Ifges(7,4) * t580 + Ifges(7,2) * t579 + Ifges(7,6) * t595;
t524 = mrSges(7,2) * t532 - mrSges(7,3) * t529 + Ifges(7,1) * t559 + Ifges(7,4) * t558 + Ifges(7,5) * t578 + t553 * t579 - t554 * t595;
t499 = -mrSges(5,1) * t542 - mrSges(6,1) * t534 + mrSges(6,2) * t539 + mrSges(5,3) * t538 - pkin(4) * t518 - pkin(5) * t623 - pkin(10) * t640 + t637 * qJD(4) + t648 * qJDD(4) - t612 * t523 - t609 * t524 + t650 * t584 + t585 * t657 - t639 * t636;
t500 = mrSges(6,1) * t535 + mrSges(7,1) * t529 + mrSges(5,2) * t542 - mrSges(7,2) * t530 - mrSges(5,3) * t537 - mrSges(6,3) * t539 + Ifges(7,5) * t559 + Ifges(7,6) * t558 + Ifges(7,3) * t578 + pkin(5) * t519 - qJ(5) * t518 + t580 * t554 - t579 * t555 + t650 * t585 + t658 * t584 + t649 * qJDD(4) + t638 * qJD(4) + t639 * t635;
t486 = mrSges(4,2) * t547 - mrSges(4,3) * t544 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t616 - pkin(9) * t512 - t499 * t610 + t500 * t613;
t491 = t616 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t512 + mrSges(4,3) * t545 - mrSges(4,1) * t547 - pkin(4) * (-qJD(4) * t591 + t622) - qJ(5) * t619 - mrSges(6,2) * t535 + mrSges(6,3) * t534 - t612 * t524 + t609 * t523 + pkin(10) * t519 - mrSges(5,1) * t537 + mrSges(5,2) * t538 + (-mrSges(6,1) * qJ(5) - t648) * t585 - t649 * t584 + (mrSges(6,2) * pkin(4) - t656) * qJDD(4) + (t637 * t613 + (pkin(4) * t582 + t638) * t610) * qJD(3);
t625 = pkin(8) * t505 + t486 * t611 + t491 * t614;
t485 = mrSges(4,1) * t544 - mrSges(4,2) * t545 + Ifges(4,3) * qJDD(3) + pkin(3) * t617 + pkin(9) * t630 + t613 * t499 + t610 * t500;
t479 = -mrSges(3,1) * t564 + mrSges(3,3) * t552 - pkin(2) * t497 - t603 * t485 + t625 * t607;
t480 = mrSges(3,2) * t564 - mrSges(3,3) * t551 + t614 * t486 - t611 * t491 + (-t497 * t603 - t498 * t607) * pkin(8);
t624 = qJ(2) * t490 + t479 * t605 + t480 * t601;
t478 = mrSges(3,1) * t551 - mrSges(3,2) * t552 + pkin(2) * t498 + t607 * t485 + t625 * t603;
t477 = mrSges(2,2) * t600 - mrSges(2,3) * t587 - t601 * t479 + t605 * t480 + (-t483 * t604 - t484 * t608) * qJ(2);
t476 = -mrSges(2,1) * t600 + mrSges(2,3) * t588 - pkin(1) * t483 - t604 * t478 + t624 * t608;
t1 = [-m(1) * g(1) + t631; -m(1) * g(2) + t641; -m(1) * g(3) + m(2) * t600 + t483; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t641 - t602 * t476 + t606 * t477; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t631 + t606 * t476 + t602 * t477; -mrSges(1,1) * g(2) + mrSges(2,1) * t587 + mrSges(1,2) * g(1) - mrSges(2,2) * t588 + pkin(1) * t484 + t608 * t478 + t624 * t604;];
tauB  = t1;
