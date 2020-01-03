% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:19
% EndTime: 2019-12-31 17:15:21
% DurationCPUTime: 1.67s
% Computational Cost: add. (11861->218), mult. (25020->267), div. (0->0), fcn. (14883->6), ass. (0->88)
t595 = Ifges(4,4) + Ifges(5,4);
t605 = Ifges(4,2) + Ifges(5,2);
t601 = Ifges(4,6) + Ifges(5,6);
t568 = sin(qJ(3));
t569 = sin(qJ(2));
t571 = cos(qJ(3));
t572 = cos(qJ(2));
t544 = (-t568 * t569 + t571 * t572) * qJD(1);
t588 = qJD(1) * qJD(2);
t552 = t569 * qJDD(1) + t572 * t588;
t553 = t572 * qJDD(1) - t569 * t588;
t518 = t544 * qJD(3) + t571 * t552 + t568 * t553;
t545 = (t568 * t572 + t569 * t571) * qJD(1);
t529 = -t544 * mrSges(5,1) + t545 * mrSges(5,2);
t570 = sin(qJ(1));
t573 = cos(qJ(1));
t559 = -t573 * g(1) - t570 * g(2);
t574 = qJD(1) ^ 2;
t547 = -t574 * pkin(1) + qJDD(1) * pkin(5) + t559;
t594 = t569 * t547;
t596 = pkin(2) * t574;
t506 = qJDD(2) * pkin(2) - t552 * pkin(6) - t594 + (pkin(6) * t588 + t569 * t596 - g(3)) * t572;
t533 = -t569 * g(3) + t572 * t547;
t590 = qJD(1) * t569;
t557 = qJD(2) * pkin(2) - pkin(6) * t590;
t567 = t572 ^ 2;
t507 = t553 * pkin(6) - qJD(2) * t557 - t567 * t596 + t533;
t499 = t571 * t506 - t568 * t507;
t564 = qJDD(2) + qJDD(3);
t565 = qJD(2) + qJD(3);
t493 = -0.2e1 * qJD(4) * t545 + (t544 * t565 - t518) * qJ(4) + (t544 * t545 + t564) * pkin(3) + t499;
t534 = -t565 * mrSges(5,2) + t544 * mrSges(5,3);
t587 = m(5) * t493 + t564 * mrSges(5,1) + t565 * t534;
t489 = -t518 * mrSges(5,3) - t545 * t529 + t587;
t500 = t568 * t506 + t571 * t507;
t517 = -t545 * qJD(3) - t568 * t552 + t571 * t553;
t536 = t565 * pkin(3) - t545 * qJ(4);
t540 = t544 ^ 2;
t495 = -t540 * pkin(3) + t517 * qJ(4) + 0.2e1 * qJD(4) * t544 - t565 * t536 + t500;
t602 = Ifges(4,5) + Ifges(5,5);
t603 = Ifges(4,1) + Ifges(5,1);
t591 = -t595 * t544 - t603 * t545 - t602 * t565;
t599 = t605 * t544 + t595 * t545 + t601 * t565;
t600 = Ifges(4,3) + Ifges(5,3);
t604 = mrSges(4,1) * t499 + mrSges(5,1) * t493 - mrSges(4,2) * t500 - mrSges(5,2) * t495 + pkin(3) * t489 + t601 * t517 + t602 * t518 + t591 * t544 + t599 * t545 + t600 * t564;
t530 = -t544 * mrSges(4,1) + t545 * mrSges(4,2);
t535 = -t565 * mrSges(4,2) + t544 * mrSges(4,3);
t484 = m(4) * t499 + t564 * mrSges(4,1) + t565 * t535 + (-t529 - t530) * t545 + (-mrSges(4,3) - mrSges(5,3)) * t518 + t587;
t537 = t565 * mrSges(5,1) - t545 * mrSges(5,3);
t538 = t565 * mrSges(4,1) - t545 * mrSges(4,3);
t586 = m(5) * t495 + t517 * mrSges(5,3) + t544 * t529;
t487 = m(4) * t500 + t517 * mrSges(4,3) + t544 * t530 + (-t537 - t538) * t565 + (-mrSges(4,2) - mrSges(5,2)) * t564 + t586;
t479 = t571 * t484 + t568 * t487;
t532 = -t572 * g(3) - t594;
t542 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t569 + Ifges(3,2) * t572) * qJD(1);
t543 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t569 + Ifges(3,4) * t572) * qJD(1);
t598 = mrSges(3,1) * t532 - mrSges(3,2) * t533 + Ifges(3,5) * t552 + Ifges(3,6) * t553 + Ifges(3,3) * qJDD(2) + pkin(2) * t479 + (t569 * t542 - t572 * t543) * qJD(1) + t604;
t551 = (-mrSges(3,1) * t572 + mrSges(3,2) * t569) * qJD(1);
t589 = qJD(1) * t572;
t556 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t589;
t477 = m(3) * t532 + qJDD(2) * mrSges(3,1) - t552 * mrSges(3,3) + qJD(2) * t556 - t551 * t590 + t479;
t555 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t590;
t582 = -t568 * t484 + t571 * t487;
t478 = m(3) * t533 - qJDD(2) * mrSges(3,2) + t553 * mrSges(3,3) - qJD(2) * t555 + t551 * t589 + t582;
t583 = -t569 * t477 + t572 * t478;
t469 = m(2) * t559 - t574 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t583;
t558 = t570 * g(1) - t573 * g(2);
t580 = -qJDD(1) * pkin(1) - t558;
t546 = -t574 * pkin(5) + t580;
t519 = -t553 * pkin(2) + t557 * t590 + (-pkin(6) * t567 - pkin(5)) * t574 + t580;
t497 = -t517 * pkin(3) - t540 * qJ(4) + t545 * t536 + qJDD(4) + t519;
t490 = m(5) * t497 - t517 * mrSges(5,1) + t518 * mrSges(5,2) - t544 * t534 + t545 * t537;
t578 = m(4) * t519 - t517 * mrSges(4,1) + t518 * mrSges(4,2) - t544 * t535 + t545 * t538 + t490;
t576 = -m(3) * t546 + t553 * mrSges(3,1) - t552 * mrSges(3,2) - t555 * t590 + t556 * t589 - t578;
t481 = m(2) * t558 + qJDD(1) * mrSges(2,1) - t574 * mrSges(2,2) + t576;
t593 = t570 * t469 + t573 * t481;
t471 = t572 * t477 + t569 * t478;
t592 = -t601 * t544 - t602 * t545 - t600 * t565;
t584 = t573 * t469 - t570 * t481;
t472 = -mrSges(4,1) * t519 + mrSges(4,3) * t500 - mrSges(5,1) * t497 + mrSges(5,3) * t495 - pkin(3) * t490 + qJ(4) * t586 + (-qJ(4) * t537 - t591) * t565 + (-qJ(4) * mrSges(5,2) + t601) * t564 + t592 * t545 + t595 * t518 + t605 * t517;
t473 = mrSges(4,2) * t519 + mrSges(5,2) * t497 - mrSges(4,3) * t499 - mrSges(5,3) * t493 - qJ(4) * t489 + t595 * t517 + t603 * t518 - t592 * t544 + t602 * t564 - t599 * t565;
t541 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t569 + Ifges(3,6) * t572) * qJD(1);
t464 = -mrSges(3,1) * t546 + mrSges(3,3) * t533 + Ifges(3,4) * t552 + Ifges(3,2) * t553 + Ifges(3,6) * qJDD(2) - pkin(2) * t578 + pkin(6) * t582 + qJD(2) * t543 + t571 * t472 + t568 * t473 - t541 * t590;
t466 = mrSges(3,2) * t546 - mrSges(3,3) * t532 + Ifges(3,1) * t552 + Ifges(3,4) * t553 + Ifges(3,5) * qJDD(2) - pkin(6) * t479 - qJD(2) * t542 - t568 * t472 + t571 * t473 + t541 * t589;
t579 = mrSges(2,1) * t558 - mrSges(2,2) * t559 + Ifges(2,3) * qJDD(1) + pkin(1) * t576 + pkin(5) * t583 + t572 * t464 + t569 * t466;
t462 = mrSges(2,1) * g(3) + mrSges(2,3) * t559 + t574 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t471 - t598;
t461 = -mrSges(2,2) * g(3) - mrSges(2,3) * t558 + Ifges(2,5) * qJDD(1) - t574 * Ifges(2,6) - pkin(5) * t471 - t569 * t464 + t572 * t466;
t1 = [-m(1) * g(1) + t584; -m(1) * g(2) + t593; (-m(1) - m(2)) * g(3) + t471; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t593 + t573 * t461 - t570 * t462; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t584 + t570 * t461 + t573 * t462; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t579; t579; t598; t604; t490;];
tauJB = t1;
