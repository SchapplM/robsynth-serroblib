% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR8
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:56
% EndTime: 2019-12-31 17:07:58
% DurationCPUTime: 1.24s
% Computational Cost: add. (7789->221), mult. (16127->273), div. (0->0), fcn. (8131->6), ass. (0->92)
t612 = Ifges(3,1) + Ifges(4,1);
t603 = Ifges(3,4) - Ifges(4,5);
t602 = Ifges(3,5) + Ifges(4,4);
t611 = Ifges(3,2) + Ifges(4,3);
t601 = Ifges(3,6) - Ifges(4,6);
t610 = Ifges(3,3) + Ifges(4,2);
t572 = sin(qJ(2));
t575 = cos(qJ(2));
t544 = (-mrSges(4,1) * t575 - mrSges(4,3) * t572) * qJD(1);
t592 = qJD(1) * qJD(2);
t591 = t575 * t592;
t546 = t572 * qJDD(1) + t591;
t573 = sin(qJ(1));
t576 = cos(qJ(1));
t557 = -t576 * g(1) - t573 * g(2);
t578 = qJD(1) ^ 2;
t537 = -t578 * pkin(1) + qJDD(1) * pkin(5) + t557;
t520 = -t572 * g(3) + t575 * t537;
t543 = (-pkin(2) * t575 - qJ(3) * t572) * qJD(1);
t577 = qJD(2) ^ 2;
t593 = qJD(1) * t575;
t606 = 2 * qJD(3);
t505 = -t577 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t606 + t543 * t593 + t520;
t547 = t575 * qJDD(1) - t572 * t592;
t594 = qJD(1) * t572;
t555 = -qJD(2) * pkin(3) - pkin(6) * t594;
t570 = t575 ^ 2;
t500 = -t570 * t578 * pkin(3) - t547 * pkin(6) + qJD(2) * t555 + t505;
t519 = -t575 * g(3) - t572 * t537;
t510 = -qJDD(2) * pkin(2) - t577 * qJ(3) + t543 * t594 + qJDD(3) - t519;
t501 = (-t546 + t591) * pkin(6) + (-t572 * t575 * t578 - qJDD(2)) * pkin(3) + t510;
t571 = sin(qJ(4));
t574 = cos(qJ(4));
t496 = -t571 * t500 + t574 * t501;
t534 = (-t571 * t572 - t574 * t575) * qJD(1);
t512 = t534 * qJD(4) + t574 * t546 - t571 * t547;
t535 = (-t571 * t575 + t572 * t574) * qJD(1);
t518 = -t534 * mrSges(5,1) + t535 * mrSges(5,2);
t565 = -qJD(2) + qJD(4);
t521 = -t565 * mrSges(5,2) + t534 * mrSges(5,3);
t564 = -qJDD(2) + qJDD(4);
t493 = m(5) * t496 + t564 * mrSges(5,1) - t512 * mrSges(5,3) - t535 * t518 + t565 * t521;
t497 = t574 * t500 + t571 * t501;
t511 = -t535 * qJD(4) - t571 * t546 - t574 * t547;
t522 = t565 * mrSges(5,1) - t535 * mrSges(5,3);
t494 = m(5) * t497 - t564 * mrSges(5,2) + t511 * mrSges(5,3) + t534 * t518 - t565 * t522;
t485 = t574 * t493 + t571 * t494;
t554 = mrSges(4,2) * t593 + qJD(2) * mrSges(4,3);
t582 = -m(4) * t510 + qJDD(2) * mrSges(4,1) + qJD(2) * t554 - t485;
t484 = t546 * mrSges(4,2) + t544 * t594 - t582;
t514 = Ifges(5,4) * t535 + Ifges(5,2) * t534 + Ifges(5,6) * t565;
t515 = Ifges(5,1) * t535 + Ifges(5,4) * t534 + Ifges(5,5) * t565;
t581 = -mrSges(5,1) * t496 + mrSges(5,2) * t497 - Ifges(5,5) * t512 - Ifges(5,6) * t511 - Ifges(5,3) * t564 - t535 * t514 + t534 * t515;
t552 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t594;
t588 = -t571 * t493 + t574 * t494;
t584 = m(4) * t505 + qJDD(2) * mrSges(4,3) + qJD(2) * t552 + t544 * t593 + t588;
t595 = t602 * qJD(2) + (t612 * t572 + t603 * t575) * qJD(1);
t596 = -t601 * qJD(2) + (-t603 * t572 - t611 * t575) * qJD(1);
t609 = -(t596 * t572 + t595 * t575) * qJD(1) + t610 * qJDD(2) + t602 * t546 + t601 * t547 + mrSges(3,1) * t519 - mrSges(4,1) * t510 - mrSges(3,2) * t520 + mrSges(4,3) * t505 - pkin(2) * t484 - pkin(3) * t485 + qJ(3) * (t547 * mrSges(4,2) + t584) + t581;
t605 = t578 * pkin(5);
t604 = mrSges(3,3) + mrSges(4,2);
t599 = qJ(3) * t575;
t545 = (-mrSges(3,1) * t575 + mrSges(3,2) * t572) * qJD(1);
t551 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t594;
t481 = m(3) * t520 - qJDD(2) * mrSges(3,2) - qJD(2) * t551 + t545 * t593 + t604 * t547 + t584;
t553 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t593;
t482 = m(3) * t519 + qJDD(2) * mrSges(3,1) + qJD(2) * t553 - t604 * t546 + (-t544 - t545) * t594 + t582;
t589 = t575 * t481 - t572 * t482;
t475 = m(2) * t557 - t578 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t589;
t556 = t573 * g(1) - t576 * g(2);
t587 = qJDD(1) * pkin(1) + t556;
t585 = -t546 * qJ(3) - t587;
t502 = -t547 * pkin(2) - t605 + (-0.2e1 * qJD(3) * t572 + (pkin(2) * t572 - t599) * qJD(2)) * qJD(1) + t585;
t499 = (-pkin(6) * t570 + pkin(5)) * t578 + (pkin(2) + pkin(3)) * t547 + (qJD(2) * t599 + (-pkin(2) * qJD(2) + t555 + t606) * t572) * qJD(1) - t585;
t586 = -m(5) * t499 + t511 * mrSges(5,1) - t512 * mrSges(5,2) + t534 * t521 - t535 * t522;
t491 = m(4) * t502 - t547 * mrSges(4,1) - t546 * mrSges(4,3) - t552 * t594 - t554 * t593 + t586;
t536 = -t587 - t605;
t580 = -m(3) * t536 + t547 * mrSges(3,1) - t546 * mrSges(3,2) - t551 * t594 + t553 * t593 - t491;
t487 = m(2) * t556 + qJDD(1) * mrSges(2,1) - t578 * mrSges(2,2) + t580;
t598 = t573 * t475 + t576 * t487;
t477 = t572 * t481 + t575 * t482;
t597 = t610 * qJD(2) + (t602 * t572 + t601 * t575) * qJD(1);
t590 = t576 * t475 - t573 * t487;
t513 = Ifges(5,5) * t535 + Ifges(5,6) * t534 + Ifges(5,3) * t565;
t488 = -mrSges(5,1) * t499 + mrSges(5,3) * t497 + Ifges(5,4) * t512 + Ifges(5,2) * t511 + Ifges(5,6) * t564 - t535 * t513 + t565 * t515;
t489 = mrSges(5,2) * t499 - mrSges(5,3) * t496 + Ifges(5,1) * t512 + Ifges(5,4) * t511 + Ifges(5,5) * t564 + t534 * t513 - t565 * t514;
t470 = -mrSges(3,1) * t536 - mrSges(4,1) * t502 + mrSges(4,2) * t505 + mrSges(3,3) * t520 - pkin(2) * t491 - pkin(3) * t586 - pkin(6) * t588 + t595 * qJD(2) + t601 * qJDD(2) - t574 * t488 - t571 * t489 + t603 * t546 + t611 * t547 - t597 * t594;
t472 = mrSges(3,2) * t536 + mrSges(4,2) * t510 - mrSges(3,3) * t519 - mrSges(4,3) * t502 - pkin(6) * t485 - qJ(3) * t491 + t596 * qJD(2) + t602 * qJDD(2) - t571 * t488 + t574 * t489 + t612 * t546 + t603 * t547 + t597 * t593;
t583 = mrSges(2,1) * t556 - mrSges(2,2) * t557 + Ifges(2,3) * qJDD(1) + pkin(1) * t580 + pkin(5) * t589 + t575 * t470 + t572 * t472;
t468 = mrSges(2,1) * g(3) + mrSges(2,3) * t557 + t578 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t477 - t609;
t467 = -mrSges(2,2) * g(3) - mrSges(2,3) * t556 + Ifges(2,5) * qJDD(1) - t578 * Ifges(2,6) - pkin(5) * t477 - t572 * t470 + t575 * t472;
t1 = [-m(1) * g(1) + t590; -m(1) * g(2) + t598; (-m(1) - m(2)) * g(3) + t477; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t598 + t576 * t467 - t573 * t468; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t590 + t573 * t467 + t576 * t468; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t583; t583; t609; t484; -t581;];
tauJB = t1;
