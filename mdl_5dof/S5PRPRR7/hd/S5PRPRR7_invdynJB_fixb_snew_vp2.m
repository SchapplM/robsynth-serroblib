% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:53
% EndTime: 2019-12-05 15:59:56
% DurationCPUTime: 1.50s
% Computational Cost: add. (13788->213), mult. (25344->263), div. (0->0), fcn. (14707->8), ass. (0->94)
t543 = sin(pkin(8));
t570 = cos(pkin(8));
t525 = -t570 * g(1) - t543 * g(2);
t540 = -g(3) + qJDD(1);
t546 = sin(qJ(2));
t549 = cos(qJ(2));
t505 = t549 * t525 + t546 * t540;
t558 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t505;
t504 = -t546 * t525 + t549 * t540;
t550 = qJD(2) ^ 2;
t557 = -t550 * qJ(3) + qJDD(3) - t504;
t576 = -pkin(2) - pkin(6);
t499 = t576 * qJDD(2) + t557;
t548 = cos(qJ(4));
t496 = t548 * t499;
t545 = sin(qJ(4));
t565 = qJD(2) * qJD(4);
t523 = t548 * qJDD(2) - t545 * t565;
t524 = t543 * g(1) - t570 * g(2);
t575 = pkin(4) * t550;
t477 = (qJDD(4) * pkin(4)) - t523 * pkin(7) + t496 + (-pkin(7) * t565 - t548 * t575 + t524) * t545;
t487 = t545 * t499 - t548 * t524;
t522 = -t545 * qJDD(2) - t548 * t565;
t566 = qJD(2) * t548;
t528 = (qJD(4) * pkin(4)) - pkin(7) * t566;
t539 = t545 ^ 2;
t478 = t522 * pkin(7) - qJD(4) * t528 - t539 * t575 + t487;
t544 = sin(qJ(5));
t547 = cos(qJ(5));
t476 = t544 * t477 + t547 * t478;
t480 = t528 * t566 - t522 * pkin(4) + (-pkin(7) * t539 + t576) * t550 + t558;
t512 = (-t544 * t545 + t547 * t548) * qJD(2);
t488 = -t512 * qJD(5) + t547 * t522 - t544 * t523;
t511 = (-t544 * t548 - t545 * t547) * qJD(2);
t489 = t511 * qJD(5) + t544 * t522 + t547 * t523;
t535 = qJD(4) + qJD(5);
t490 = Ifges(6,5) * t512 + Ifges(6,6) * t511 + Ifges(6,3) * t535;
t492 = Ifges(6,1) * t512 + Ifges(6,4) * t511 + Ifges(6,5) * t535;
t534 = qJDD(4) + qJDD(5);
t463 = -mrSges(6,1) * t480 + mrSges(6,3) * t476 + Ifges(6,4) * t489 + Ifges(6,2) * t488 + Ifges(6,6) * t534 - t512 * t490 + t535 * t492;
t475 = t547 * t477 - t544 * t478;
t491 = Ifges(6,4) * t512 + Ifges(6,2) * t511 + Ifges(6,6) * t535;
t464 = mrSges(6,2) * t480 - mrSges(6,3) * t475 + Ifges(6,1) * t489 + Ifges(6,4) * t488 + Ifges(6,5) * t534 + t511 * t490 - t535 * t491;
t498 = t576 * t550 + t558;
t508 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t548 - Ifges(5,6) * t545) * qJD(2);
t510 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t548 - Ifges(5,4) * t545) * qJD(2);
t502 = -t535 * mrSges(6,2) + t511 * mrSges(6,3);
t503 = t535 * mrSges(6,1) - t512 * mrSges(6,3);
t556 = m(6) * t480 - t488 * mrSges(6,1) + t489 * mrSges(6,2) - t511 * t502 + t512 * t503;
t494 = -t511 * mrSges(6,1) + t512 * mrSges(6,2);
t472 = m(6) * t475 + t534 * mrSges(6,1) - t489 * mrSges(6,3) - t512 * t494 + t535 * t502;
t473 = m(6) * t476 - t534 * mrSges(6,2) + t488 * mrSges(6,3) + t511 * t494 - t535 * t503;
t559 = -t544 * t472 + t547 * t473;
t447 = -mrSges(5,1) * t498 + mrSges(5,3) * t487 + Ifges(5,4) * t523 + Ifges(5,2) * t522 + Ifges(5,6) * qJDD(4) - pkin(4) * t556 + pkin(7) * t559 + qJD(4) * t510 + t547 * t463 + t544 * t464 - t508 * t566;
t462 = t547 * t472 + t544 * t473;
t486 = t545 * t524 + t496;
t509 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t548 - Ifges(5,2) * t545) * qJD(2);
t567 = qJD(2) * t545;
t449 = mrSges(5,2) * t498 - mrSges(5,3) * t486 + Ifges(5,1) * t523 + Ifges(5,4) * t522 + Ifges(5,5) * qJDD(4) - pkin(7) * t462 - qJD(4) * t509 - t544 * t463 + t547 * t464 - t508 * t567;
t521 = (mrSges(5,1) * t545 + mrSges(5,2) * t548) * qJD(2);
t526 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t567;
t459 = m(5) * t486 + qJDD(4) * mrSges(5,1) - t523 * mrSges(5,3) + qJD(4) * t526 - t521 * t566 + t462;
t527 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t566;
t460 = m(5) * t487 - qJDD(4) * mrSges(5,2) + t522 * mrSges(5,3) - qJD(4) * t527 - t521 * t567 + t559;
t456 = t548 * t459 + t545 * t460;
t501 = -qJDD(2) * pkin(2) + t557;
t555 = -m(4) * t501 + (t550 * mrSges(4,3)) - t456;
t453 = qJDD(2) * mrSges(4,2) - t555;
t500 = t550 * pkin(2) - t558;
t553 = -m(5) * t498 + t522 * mrSges(5,1) - t523 * mrSges(5,2) - t526 * t567 - t527 * t566 - t556;
t468 = -m(4) * t500 + t550 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t553;
t578 = mrSges(3,1) * t504 - mrSges(3,2) * t505 + mrSges(4,2) * t501 - mrSges(4,3) * t500 - pkin(2) * t453 - pkin(6) * t456 + qJ(3) * t468 - t545 * t447 + t548 * t449 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t577 = m(3) + m(4);
t574 = mrSges(3,1) - mrSges(4,2);
t573 = Ifges(3,5) - Ifges(4,4);
t572 = (-Ifges(3,6) + Ifges(4,5));
t451 = m(3) * t504 - t550 * mrSges(3,2) + t574 * qJDD(2) + t555;
t467 = m(3) * t505 - t550 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t468;
t560 = -t546 * t451 + t549 * t467;
t446 = m(2) * t525 + t560;
t568 = -t545 * t459 + t548 * t460;
t454 = (m(2) + t577) * t524 - t568;
t569 = t543 * t446 + t570 * t454;
t448 = t549 * t451 + t546 * t467;
t562 = m(2) * t540 + t448;
t561 = t570 * t446 - t543 * t454;
t554 = mrSges(6,1) * t475 - mrSges(6,2) * t476 + Ifges(6,5) * t489 + Ifges(6,6) * t488 + Ifges(6,3) * t534 + t512 * t491 - t511 * t492;
t552 = mrSges(5,1) * t486 - mrSges(5,2) * t487 + Ifges(5,5) * t523 + Ifges(5,6) * t522 + (Ifges(5,3) * qJDD(4)) + pkin(4) * t462 + t509 * t566 + t510 * t567 + t554;
t455 = -m(4) * t524 + t568;
t443 = mrSges(4,1) * t501 - mrSges(3,3) * t504 + t552 - qJ(3) * t455 + pkin(3) * t456 + t573 * qJDD(2) + (t572 * t550) + (-mrSges(3,2) + mrSges(4,3)) * t524;
t442 = -mrSges(4,1) * t500 + mrSges(3,3) * t505 - pkin(2) * t455 - pkin(3) * t553 - pkin(6) * t568 - t572 * qJDD(2) - t548 * t447 - t545 * t449 + t574 * t524 + t573 * t550;
t441 = -mrSges(2,1) * t540 + mrSges(2,3) * t525 - pkin(1) * t448 - t578;
t440 = mrSges(2,2) * t540 - mrSges(2,3) * t524 - pkin(5) * t448 - t546 * t442 + t549 * t443;
t1 = [-m(1) * g(1) + t561; -m(1) * g(2) + t569; -m(1) * g(3) + t562; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t569 + t570 * t440 - t543 * t441; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t561 + t543 * t440 + t570 * t441; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t525 + t546 * t443 + t549 * t442 - pkin(1) * t568 + pkin(5) * t560 + (pkin(1) * t577 + mrSges(2,1)) * t524; t562; t578; t453; t552; t554;];
tauJB = t1;
