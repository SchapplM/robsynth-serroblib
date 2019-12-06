% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:09
% EndTime: 2019-12-05 17:38:10
% DurationCPUTime: 1.13s
% Computational Cost: add. (8786->216), mult. (16715->254), div. (0->0), fcn. (8129->6), ass. (0->86)
t619 = 2 * qJD(1);
t587 = sin(qJ(1));
t590 = cos(qJ(1));
t561 = -t590 * g(1) - t587 * g(2);
t618 = qJDD(1) * qJ(2) + (qJD(2) * t619) + t561;
t560 = t587 * g(1) - t590 * g(2);
t591 = qJD(1) ^ 2;
t541 = -qJDD(1) * pkin(1) - t591 * qJ(2) + qJDD(2) - t560;
t598 = qJDD(1) * qJ(3) + (qJD(3) * t619) - t541;
t617 = -m(3) - m(4);
t616 = mrSges(2,1) - mrSges(3,2);
t615 = t591 * mrSges(4,3);
t539 = t591 * pkin(1) - t618;
t534 = qJDD(3) + (-pkin(1) - qJ(3)) * t591 + t618;
t530 = -qJDD(1) * pkin(6) + t534;
t586 = sin(qJ(4));
t589 = cos(qJ(4));
t523 = t586 * g(3) + t589 * t530;
t611 = qJD(1) * qJD(4);
t606 = t586 * t611;
t555 = t589 * qJDD(1) - t606;
t506 = (-t555 - t606) * pkin(7) + (-t586 * t589 * t591 + qJDD(4)) * pkin(4) + t523;
t524 = -t589 * g(3) + t586 * t530;
t554 = -t586 * qJDD(1) - t589 * t611;
t612 = qJD(1) * t589;
t559 = qJD(4) * pkin(4) - pkin(7) * t612;
t580 = t586 ^ 2;
t507 = -t580 * t591 * pkin(4) + t554 * pkin(7) - qJD(4) * t559 + t524;
t585 = sin(qJ(5));
t588 = cos(qJ(5));
t504 = t588 * t506 - t585 * t507;
t545 = (-t585 * t589 - t586 * t588) * qJD(1);
t516 = t545 * qJD(5) + t585 * t554 + t588 * t555;
t546 = (-t585 * t586 + t588 * t589) * qJD(1);
t525 = -t545 * mrSges(6,1) + t546 * mrSges(6,2);
t569 = qJD(4) + qJD(5);
t535 = -t569 * mrSges(6,2) + t545 * mrSges(6,3);
t568 = qJDD(4) + qJDD(5);
t501 = m(6) * t504 + t568 * mrSges(6,1) - t516 * mrSges(6,3) - t546 * t525 + t569 * t535;
t505 = t585 * t506 + t588 * t507;
t515 = -t546 * qJD(5) + t588 * t554 - t585 * t555;
t536 = t569 * mrSges(6,1) - t546 * mrSges(6,3);
t502 = m(6) * t505 - t568 * mrSges(6,2) + t515 * mrSges(6,3) + t545 * t525 - t569 * t536;
t490 = t588 * t501 + t585 * t502;
t553 = (mrSges(5,1) * t586 + mrSges(5,2) * t589) * qJD(1);
t613 = qJD(1) * t586;
t557 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t613;
t487 = m(5) * t523 + qJDD(4) * mrSges(5,1) - t555 * mrSges(5,3) + qJD(4) * t557 - t553 * t612 + t490;
t558 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t612;
t603 = -t585 * t501 + t588 * t502;
t488 = m(5) * t524 - qJDD(4) * mrSges(5,2) + t554 * mrSges(5,3) - qJD(4) * t558 - t553 * t613 + t603;
t483 = t589 * t487 + t586 * t488;
t602 = m(4) * t534 + qJDD(1) * mrSges(4,2) + t483;
t597 = -m(3) * t539 + (t591 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t602;
t479 = m(2) * t561 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t591) + t597;
t529 = -t591 * pkin(6) + t598;
t509 = t559 * t612 - t554 * pkin(4) + (-pkin(7) * t580 - pkin(6)) * t591 + t598;
t601 = m(6) * t509 - t515 * mrSges(6,1) + t516 * mrSges(6,2) - t545 * t535 + t546 * t536;
t596 = -m(5) * t529 + t554 * mrSges(5,1) - t555 * mrSges(5,2) - t557 * t613 - t558 * t612 - t601;
t497 = -m(4) * t598 - t591 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t596;
t594 = -m(3) * t541 + t591 * mrSges(3,3) - t497;
t494 = m(2) * t560 - t591 * mrSges(2,2) + t616 * qJDD(1) + t594;
t614 = t587 * t479 + t590 * t494;
t608 = (Ifges(2,5) - Ifges(3,4) + Ifges(4,5));
t607 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t605 = t590 * t479 - t587 * t494;
t604 = -t586 * t487 + t589 * t488;
t518 = Ifges(6,4) * t546 + Ifges(6,2) * t545 + Ifges(6,6) * t569;
t519 = Ifges(6,1) * t546 + Ifges(6,4) * t545 + Ifges(6,5) * t569;
t595 = mrSges(6,1) * t504 - mrSges(6,2) * t505 + Ifges(6,5) * t516 + Ifges(6,6) * t515 + Ifges(6,3) * t568 + t546 * t518 - t545 * t519;
t543 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t589 - Ifges(5,2) * t586) * qJD(1);
t544 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t589 - Ifges(5,4) * t586) * qJD(1);
t593 = mrSges(5,1) * t523 - mrSges(5,2) * t524 + Ifges(5,5) * t555 + Ifges(5,6) * t554 + Ifges(5,3) * qJDD(4) + pkin(4) * t490 + t543 * t612 + t544 * t613 + t595;
t517 = Ifges(6,5) * t546 + Ifges(6,6) * t545 + Ifges(6,3) * t569;
t491 = -mrSges(6,1) * t509 + mrSges(6,3) * t505 + Ifges(6,4) * t516 + Ifges(6,2) * t515 + Ifges(6,6) * t568 - t546 * t517 + t569 * t519;
t492 = mrSges(6,2) * t509 - mrSges(6,3) * t504 + Ifges(6,1) * t516 + Ifges(6,4) * t515 + Ifges(6,5) * t568 + t545 * t517 - t569 * t518;
t542 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t589 - Ifges(5,6) * t586) * qJD(1);
t474 = -mrSges(5,1) * t529 + mrSges(5,3) * t524 + Ifges(5,4) * t555 + Ifges(5,2) * t554 + Ifges(5,6) * qJDD(4) - pkin(4) * t601 + pkin(7) * t603 + qJD(4) * t544 + t588 * t491 + t585 * t492 - t542 * t612;
t476 = mrSges(5,2) * t529 - mrSges(5,3) * t523 + Ifges(5,1) * t555 + Ifges(5,4) * t554 + Ifges(5,5) * qJDD(4) - pkin(7) * t490 - qJD(4) * t543 - t585 * t491 + t588 * t492 - t542 * t613;
t496 = qJDD(1) * mrSges(3,2) - t594;
t592 = -mrSges(2,2) * t561 - mrSges(3,3) * t539 + mrSges(4,3) * t598 - pkin(6) * t483 - qJ(3) * t497 - t586 * t474 + t589 * t476 + qJ(2) * (t597 - t615) - pkin(1) * t496 + mrSges(4,2) * t534 + mrSges(3,2) * t541 + mrSges(2,1) * t560 + (Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);
t482 = t617 * g(3) + t604;
t481 = t602 - t615;
t473 = mrSges(2,3) * t561 - qJ(3) * t604 + pkin(3) * t483 + pkin(2) * t481 - pkin(1) * t482 + mrSges(4,1) * t534 - mrSges(3,1) * t539 + (t608 * t591) + t607 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t616) * g(3) + t593;
t472 = -qJ(2) * t482 - mrSges(2,3) * t560 + pkin(2) * t497 + mrSges(3,1) * t541 + t589 * t474 + pkin(3) * t596 + pkin(6) * t604 + t586 * t476 - mrSges(4,1) * t598 - t607 * t591 + t608 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t605; -m(1) * g(2) + t614; (-m(1) - m(2) + t617) * g(3) + t604; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t614 + t590 * t472 - t587 * t473; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t605 + t587 * t472 + t590 * t473; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t592; t592; t496; t481; t593; t595;];
tauJB = t1;
