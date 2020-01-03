% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:50
% EndTime: 2019-12-31 17:50:51
% DurationCPUTime: 1.19s
% Computational Cost: add. (7139->201), mult. (12733->231), div. (0->0), fcn. (5413->6), ass. (0->86)
t611 = Ifges(5,1) + Ifges(6,1);
t599 = Ifges(5,4) + Ifges(6,4);
t609 = Ifges(5,5) + Ifges(6,5);
t610 = Ifges(5,2) + Ifges(6,2);
t608 = Ifges(5,6) + Ifges(6,6);
t607 = Ifges(5,3) + Ifges(6,3);
t570 = sin(qJ(4));
t572 = cos(qJ(4));
t606 = t608 * qJD(4) + (-t610 * t570 + t599 * t572) * qJD(1);
t605 = t609 * qJD(4) + (-t599 * t570 + t611 * t572) * qJD(1);
t571 = sin(qJ(1));
t573 = cos(qJ(1));
t549 = t571 * g(1) - t573 * g(2);
t536 = qJDD(1) * pkin(1) + t549;
t550 = -t573 * g(1) - t571 * g(2);
t574 = qJD(1) ^ 2;
t539 = -t574 * pkin(1) + t550;
t568 = sin(pkin(7));
t569 = cos(pkin(7));
t511 = t568 * t536 + t569 * t539;
t581 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t511;
t603 = -pkin(2) - pkin(6);
t504 = t603 * t574 + t581;
t591 = qJD(1) * qJD(4);
t540 = -t570 * qJDD(1) - t572 * t591;
t541 = t572 * qJDD(1) - t570 * t591;
t593 = qJD(1) * t570;
t545 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t593;
t592 = qJD(1) * t572;
t548 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t592;
t546 = qJD(4) * pkin(4) - qJ(5) * t592;
t564 = t570 ^ 2;
t497 = t546 * t592 - t540 * pkin(4) + qJDD(5) + (-qJ(5) * t564 + t603) * t574 + t581;
t544 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t593;
t547 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t592;
t582 = m(6) * t497 + t541 * mrSges(6,2) + t544 * t593 + t547 * t592;
t604 = -m(5) * t504 - t541 * mrSges(5,2) + (mrSges(5,1) + mrSges(6,1)) * t540 - t545 * t593 - t548 * t592 - t582;
t602 = pkin(4) * t574;
t601 = mrSges(3,1) - mrSges(4,2);
t598 = Ifges(3,5) - Ifges(4,4);
t597 = Ifges(3,6) - Ifges(4,5);
t510 = t569 * t536 - t568 * t539;
t580 = -t574 * qJ(3) + qJDD(3) - t510;
t505 = t603 * qJDD(1) + t580;
t502 = t572 * t505;
t565 = -g(3) + qJDD(2);
t499 = -t570 * t565 + t502;
t537 = (mrSges(6,1) * t570 + mrSges(6,2) * t572) * qJD(1);
t583 = qJD(1) * (-t537 - (mrSges(5,1) * t570 + mrSges(5,2) * t572) * qJD(1));
t589 = -0.2e1 * qJD(1) * qJD(5);
t494 = t572 * t589 + qJDD(4) * pkin(4) - t541 * qJ(5) + t502 + (-qJ(5) * t591 - t572 * t602 - t565) * t570;
t588 = m(6) * t494 + qJDD(4) * mrSges(6,1) + qJD(4) * t544;
t486 = m(5) * t499 + qJDD(4) * mrSges(5,1) + qJD(4) * t545 + (-mrSges(5,3) - mrSges(6,3)) * t541 + t572 * t583 + t588;
t500 = t570 * t505 + t572 * t565;
t495 = t540 * qJ(5) - qJD(4) * t546 - t564 * t602 + t570 * t589 + t500;
t595 = m(6) * t495 + t540 * mrSges(6,3);
t487 = m(5) * t500 + t540 * mrSges(5,3) + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t547 - t548) * qJD(4) + t570 * t583 + t595;
t480 = t572 * t486 + t570 * t487;
t508 = -qJDD(1) * pkin(2) + t580;
t578 = -m(4) * t508 + t574 * mrSges(4,3) - t480;
t473 = m(3) * t510 - t574 * mrSges(3,2) + t601 * qJDD(1) + t578;
t506 = t574 * pkin(2) - t581;
t577 = -m(4) * t506 + t574 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t604;
t483 = m(3) * t511 - t574 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t577;
t470 = t569 * t473 + t568 * t483;
t467 = m(2) * t549 + qJDD(1) * mrSges(2,1) - t574 * mrSges(2,2) + t470;
t585 = -t568 * t473 + t569 * t483;
t468 = m(2) * t550 - t574 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t585;
t596 = t573 * t467 + t571 * t468;
t594 = -t607 * qJD(4) + (t608 * t570 - t609 * t572) * qJD(1);
t586 = -t571 * t467 + t573 * t468;
t584 = -t570 * t486 + t572 * t487;
t479 = m(4) * t565 + t584;
t478 = m(3) * t565 + t479;
t489 = -t541 * mrSges(6,3) - t537 * t592 + t588;
t576 = mrSges(5,1) * t499 + mrSges(6,1) * t494 - mrSges(5,2) * t500 - mrSges(6,2) * t495 + pkin(4) * t489 + t607 * qJDD(4) + t608 * t540 + t609 * t541 + t606 * t592 + t605 * t593;
t490 = -t540 * mrSges(6,1) + t582;
t471 = -mrSges(5,1) * t504 + mrSges(5,3) * t500 - mrSges(6,1) * t497 + mrSges(6,3) * t495 - pkin(4) * t490 + qJ(5) * t595 + t599 * t541 + t610 * t540 + (-qJ(5) * mrSges(6,2) + t608) * qJDD(4) + (-qJ(5) * t547 + t605) * qJD(4) + (-qJ(5) * t537 * t570 + t594 * t572) * qJD(1);
t475 = mrSges(5,2) * t504 + mrSges(6,2) * t497 - mrSges(5,3) * t499 - mrSges(6,3) * t494 - qJ(5) * t489 - t606 * qJD(4) + t609 * qJDD(4) + t599 * t540 + t611 * t541 + t594 * t593;
t477 = qJDD(1) * mrSges(4,2) - t578;
t575 = mrSges(2,1) * t549 + mrSges(3,1) * t510 - mrSges(2,2) * t550 - mrSges(3,2) * t511 + mrSges(4,2) * t508 - mrSges(4,3) * t506 + pkin(1) * t470 - pkin(2) * t477 - pkin(6) * t480 + qJ(3) * t577 - t570 * t471 + t572 * t475 + (Ifges(2,3) + Ifges(4,1) + Ifges(3,3)) * qJDD(1);
t463 = t576 + t598 * qJDD(1) - t597 * t574 + (mrSges(3,2) - mrSges(4,3)) * t565 - mrSges(3,3) * t510 + mrSges(4,1) * t508 - qJ(3) * t479 + pkin(3) * t480;
t462 = -mrSges(4,1) * t506 + mrSges(3,3) * t511 - pkin(2) * t479 - pkin(3) * t604 - pkin(6) * t584 + t597 * qJDD(1) - t572 * t471 - t570 * t475 - t601 * t565 + t598 * t574;
t461 = -mrSges(2,2) * g(3) - mrSges(2,3) * t549 + Ifges(2,5) * qJDD(1) - t574 * Ifges(2,6) - qJ(2) * t470 - t568 * t462 + t569 * t463;
t460 = mrSges(2,1) * g(3) + mrSges(2,3) * t550 + t574 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t478 + qJ(2) * t585 + t569 * t462 + t568 * t463;
t1 = [-m(1) * g(1) + t586; -m(1) * g(2) + t596; (-m(1) - m(2)) * g(3) + t478; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t596 - t571 * t460 + t573 * t461; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t586 + t573 * t460 + t571 * t461; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t575; t575; t478; t477; t576; t490;];
tauJB = t1;
