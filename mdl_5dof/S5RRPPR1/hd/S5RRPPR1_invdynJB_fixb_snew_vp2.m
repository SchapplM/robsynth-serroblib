% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m [6x1]
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:34
% EndTime: 2022-01-20 09:51:37
% DurationCPUTime: 3.20s
% Computational Cost: add. (45658->207), mult. (62159->259), div. (0->0), fcn. (35962->10), ass. (0->98)
t573 = qJD(1) + qJD(2);
t569 = t573 ^ 2;
t578 = cos(pkin(9));
t611 = pkin(4) * t578;
t576 = sin(pkin(9));
t610 = mrSges(5,2) * t576;
t572 = t578 ^ 2;
t609 = t569 * t572;
t570 = qJDD(1) + qJDD(2);
t608 = t570 * t578;
t596 = Ifges(5,5) * t576 + Ifges(5,6) * t578;
t607 = t569 * t596;
t582 = sin(qJ(1));
t585 = cos(qJ(1));
t558 = t582 * g(1) - t585 * g(2);
t553 = qJDD(1) * pkin(1) + t558;
t559 = -t585 * g(1) - t582 * g(2);
t586 = qJD(1) ^ 2;
t554 = -t586 * pkin(1) + t559;
t581 = sin(qJ(2));
t584 = cos(qJ(2));
t542 = t584 * t553 - t581 * t554;
t539 = t570 * pkin(2) + t542;
t543 = t581 * t553 + t584 * t554;
t540 = -t569 * pkin(2) + t543;
t577 = sin(pkin(8));
t579 = cos(pkin(8));
t527 = t577 * t539 + t579 * t540;
t524 = -t569 * pkin(3) + t570 * qJ(4) + t527;
t575 = -g(3) + qJDD(3);
t604 = qJD(4) * t573;
t605 = t578 * t575 - 0.2e1 * t576 * t604;
t517 = (-pkin(7) * t570 + t569 * t611 - t524) * t576 + t605;
t521 = t576 * t575 + (t524 + 0.2e1 * t604) * t578;
t518 = -pkin(4) * t609 + pkin(7) * t608 + t521;
t580 = sin(qJ(5));
t583 = cos(qJ(5));
t515 = t583 * t517 - t580 * t518;
t592 = -t576 * t580 + t578 * t583;
t546 = t592 * t573;
t593 = t576 * t583 + t578 * t580;
t547 = t593 * t573;
t533 = -t546 * mrSges(6,1) + t547 * mrSges(6,2);
t535 = t546 * qJD(5) + t593 * t570;
t544 = -qJD(5) * mrSges(6,2) + t546 * mrSges(6,3);
t513 = m(6) * t515 + qJDD(5) * mrSges(6,1) - t535 * mrSges(6,3) + qJD(5) * t544 - t547 * t533;
t516 = t580 * t517 + t583 * t518;
t534 = -t547 * qJD(5) + t592 * t570;
t545 = qJD(5) * mrSges(6,1) - t547 * mrSges(6,3);
t514 = m(6) * t516 - qJDD(5) * mrSges(6,2) + t534 * mrSges(6,3) - qJD(5) * t545 + t546 * t533;
t503 = t583 * t513 + t580 * t514;
t520 = -t576 * t524 + t605;
t594 = mrSges(5,3) * t570 + (-mrSges(5,1) * t578 + t610) * t569;
t501 = m(5) * t520 - t594 * t576 + t503;
t599 = -t580 * t513 + t583 * t514;
t502 = m(5) * t521 + t594 * t578 + t599;
t600 = -t576 * t501 + t578 * t502;
t494 = m(4) * t527 - t569 * mrSges(4,1) - t570 * mrSges(4,2) + t600;
t526 = t579 * t539 - t577 * t540;
t595 = qJDD(4) - t526;
t523 = -t570 * pkin(3) - t569 * qJ(4) + t595;
t571 = t576 ^ 2;
t519 = (-pkin(3) - t611) * t570 + (-qJ(4) + (-t571 - t572) * pkin(7)) * t569 + t595;
t591 = m(6) * t519 - t534 * mrSges(6,1) + t535 * mrSges(6,2) - t546 * t544 + t547 * t545;
t589 = -m(5) * t523 + mrSges(5,1) * t608 - t591 + (t569 * t571 + t609) * mrSges(5,3);
t507 = m(4) * t526 - t569 * mrSges(4,2) + (mrSges(4,1) - t610) * t570 + t589;
t489 = t577 * t494 + t579 * t507;
t484 = m(3) * t542 + t570 * mrSges(3,1) - t569 * mrSges(3,2) + t489;
t601 = t579 * t494 - t577 * t507;
t485 = m(3) * t543 - t569 * mrSges(3,1) - t570 * mrSges(3,2) + t601;
t479 = t584 * t484 + t581 * t485;
t476 = m(2) * t558 + qJDD(1) * mrSges(2,1) - t586 * mrSges(2,2) + t479;
t602 = -t581 * t484 + t584 * t485;
t477 = m(2) * t559 - t586 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t602;
t606 = t585 * t476 + t582 * t477;
t497 = t578 * t501 + t576 * t502;
t495 = m(4) * t575 + t497;
t603 = -t582 * t476 + t585 * t477;
t598 = Ifges(5,1) * t576 + Ifges(5,4) * t578;
t597 = Ifges(5,4) * t576 + Ifges(5,2) * t578;
t528 = Ifges(6,5) * t547 + Ifges(6,6) * t546 + Ifges(6,3) * qJD(5);
t530 = Ifges(6,1) * t547 + Ifges(6,4) * t546 + Ifges(6,5) * qJD(5);
t504 = -mrSges(6,1) * t519 + mrSges(6,3) * t516 + Ifges(6,4) * t535 + Ifges(6,2) * t534 + Ifges(6,6) * qJDD(5) + qJD(5) * t530 - t547 * t528;
t529 = Ifges(6,4) * t547 + Ifges(6,2) * t546 + Ifges(6,6) * qJD(5);
t505 = mrSges(6,2) * t519 - mrSges(6,3) * t515 + Ifges(6,1) * t535 + Ifges(6,4) * t534 + Ifges(6,5) * qJDD(5) - qJD(5) * t529 + t546 * t528;
t487 = -mrSges(5,1) * t523 + mrSges(5,3) * t521 - pkin(4) * t591 + pkin(7) * t599 + t583 * t504 + t580 * t505 + t597 * t570 - t576 * t607;
t491 = mrSges(5,2) * t523 - mrSges(5,3) * t520 - pkin(7) * t503 - t580 * t504 + t583 * t505 + t598 * t570 + t578 * t607;
t509 = t570 * t610 - t589;
t590 = mrSges(3,1) * t542 + mrSges(4,1) * t526 - mrSges(3,2) * t543 - mrSges(4,2) * t527 + pkin(2) * t489 - pkin(3) * t509 + qJ(4) * t600 + t578 * t487 + t576 * t491 + (Ifges(3,3) + Ifges(4,3)) * t570;
t588 = mrSges(6,1) * t515 - mrSges(6,2) * t516 + Ifges(6,5) * t535 + Ifges(6,6) * t534 + Ifges(6,3) * qJDD(5) + t547 * t529 - t546 * t530;
t587 = mrSges(2,1) * t558 - mrSges(2,2) * t559 + Ifges(2,3) * qJDD(1) + pkin(1) * t479 + t590;
t480 = -mrSges(4,1) * t575 - mrSges(5,1) * t520 + mrSges(5,2) * t521 + mrSges(4,3) * t527 - pkin(3) * t497 - pkin(4) * t503 + (Ifges(4,6) - t596) * t570 - t588 + (-t576 * t597 + t578 * t598 + Ifges(4,5)) * t569;
t472 = mrSges(4,2) * t575 - mrSges(4,3) * t526 + Ifges(4,5) * t570 - t569 * Ifges(4,6) - qJ(4) * t497 - t576 * t487 + t578 * t491;
t471 = -mrSges(3,2) * g(3) - mrSges(3,3) * t542 + Ifges(3,5) * t570 - t569 * Ifges(3,6) - qJ(3) * t489 + t579 * t472 - t577 * t480;
t470 = mrSges(3,1) * g(3) + mrSges(3,3) * t543 + t569 * Ifges(3,5) + Ifges(3,6) * t570 - pkin(2) * t495 + qJ(3) * t601 + t577 * t472 + t579 * t480;
t469 = -mrSges(2,2) * g(3) - mrSges(2,3) * t558 + Ifges(2,5) * qJDD(1) - t586 * Ifges(2,6) - pkin(6) * t479 - t581 * t470 + t584 * t471;
t468 = Ifges(2,6) * qJDD(1) + t586 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t559 + t581 * t471 + t584 * t470 - pkin(1) * (-m(3) * g(3) + t495) + pkin(6) * t602;
t1 = [-m(1) * g(1) + t603; -m(1) * g(2) + t606; (-m(1) - m(2) - m(3)) * g(3) + t495; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t606 - t582 * t468 + t585 * t469; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t603 + t585 * t468 + t582 * t469; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t587; t587; t590; t495; t509; t588;];
tauJB = t1;
