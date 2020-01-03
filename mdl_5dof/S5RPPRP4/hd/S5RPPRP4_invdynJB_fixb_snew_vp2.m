% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:10
% EndTime: 2019-12-31 17:52:11
% DurationCPUTime: 1.06s
% Computational Cost: add. (8142->204), mult. (14231->233), div. (0->0), fcn. (5150->6), ass. (0->88)
t572 = Ifges(5,1) + Ifges(6,1);
t561 = Ifges(5,4) + Ifges(6,4);
t560 = Ifges(5,5) + Ifges(6,5);
t571 = Ifges(5,2) + Ifges(6,2);
t558 = Ifges(5,6) + Ifges(6,6);
t570 = Ifges(5,3) + Ifges(6,3);
t532 = sin(qJ(4));
t534 = cos(qJ(4));
t549 = qJD(1) * qJD(4);
t547 = t534 * t549;
t504 = -t532 * qJDD(1) - t547;
t536 = qJD(1) ^ 2;
t533 = sin(qJ(1));
t535 = cos(qJ(1));
t514 = -t535 * g(1) - t533 * g(2);
t540 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t514;
t566 = -pkin(1) - pkin(2);
t484 = t566 * t536 + t540;
t513 = t533 * g(1) - t535 * g(2);
t539 = -t536 * qJ(2) + qJDD(2) - t513;
t487 = t566 * qJDD(1) + t539;
t530 = sin(pkin(7));
t531 = cos(pkin(7));
t480 = t531 * t484 + t530 * t487;
t477 = -t536 * pkin(3) - qJDD(1) * pkin(6) + t480;
t527 = g(3) + qJDD(3);
t516 = t534 * t527;
t548 = qJD(1) * qJD(5);
t565 = pkin(4) * t536;
t470 = qJDD(4) * pkin(4) + t516 + (-t504 - t547) * qJ(5) + (t534 * t565 - t477 + (2 * t548)) * t532;
t502 = (mrSges(6,1) * t534 - mrSges(6,2) * t532) * qJD(1);
t550 = qJD(1) * t534;
t511 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t550;
t551 = qJD(1) * t532;
t544 = m(6) * t470 + qJDD(4) * mrSges(6,1) + qJD(4) * t511 + t502 * t551;
t466 = -t504 * mrSges(6,3) + t544;
t474 = t534 * t477 + t532 * t527;
t505 = -t534 * qJDD(1) + t532 * t549;
t508 = qJD(4) * pkin(4) + qJ(5) * t551;
t526 = t534 ^ 2;
t471 = t505 * qJ(5) - qJD(4) * t508 - t526 * t565 - 0.2e1 * t534 * t548 + t474;
t473 = -t532 * t477 + t516;
t552 = t560 * qJD(4) + (-t572 * t532 - t561 * t534) * qJD(1);
t553 = t558 * qJD(4) + (-t561 * t532 - t571 * t534) * qJD(1);
t569 = mrSges(5,1) * t473 + mrSges(6,1) * t470 - mrSges(5,2) * t474 - mrSges(6,2) * t471 + pkin(4) * t466 - (t553 * t532 - t552 * t534) * qJD(1) + t570 * qJDD(4) + t560 * t504 + t558 * t505;
t564 = mrSges(2,1) + mrSges(3,1);
t563 = -mrSges(5,2) - mrSges(6,2);
t562 = Ifges(3,4) + Ifges(2,5);
t559 = Ifges(2,6) - Ifges(3,6);
t488 = -t536 * pkin(1) + t540;
t503 = (mrSges(5,1) * t534 - mrSges(5,2) * t532) * qJD(1);
t512 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t550;
t463 = t503 * t551 + m(5) * t473 + qJDD(4) * mrSges(5,1) + qJD(4) * t512 + (-mrSges(5,3) - mrSges(6,3)) * t504 + t544;
t509 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t551;
t510 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t551;
t555 = m(6) * t471 + t505 * mrSges(6,3);
t464 = m(5) * t474 + t505 * mrSges(5,3) + t563 * qJDD(4) + (-t509 - t510) * qJD(4) + (-t502 - t503) * t550 + t555;
t460 = -t532 * t463 + t534 * t464;
t456 = m(4) * t480 - t536 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t460;
t479 = -t530 * t484 + t531 * t487;
t542 = qJDD(1) * pkin(3) - t479;
t476 = -t536 * pkin(6) + t542;
t472 = -t508 * t551 - t505 * pkin(4) + qJDD(5) + (-qJ(5) * t526 - pkin(6)) * t536 + t542;
t543 = -m(6) * t472 + t505 * mrSges(6,1) + t509 * t551;
t465 = -m(5) * t476 + t510 * t551 + t505 * mrSges(5,1) + t563 * t504 + (-t511 - t512) * t550 + t543;
t461 = m(4) * t479 - qJDD(1) * mrSges(4,1) - t536 * mrSges(4,2) + t465;
t545 = t531 * t456 - t530 * t461;
t541 = m(3) * t488 + qJDD(1) * mrSges(3,3) + t545;
t447 = m(2) * t514 - qJDD(1) * mrSges(2,2) - t564 * t536 + t541;
t452 = t530 * t456 + t531 * t461;
t489 = -qJDD(1) * pkin(1) + t539;
t451 = m(3) * t489 - qJDD(1) * mrSges(3,1) - t536 * mrSges(3,3) + t452;
t448 = m(2) * t513 + qJDD(1) * mrSges(2,1) - t536 * mrSges(2,2) - t451;
t556 = t533 * t447 + t535 * t448;
t554 = t570 * qJD(4) + (-t560 * t532 - t558 * t534) * qJD(1);
t546 = t535 * t447 - t533 * t448;
t459 = t534 * t463 + t532 * t464;
t458 = m(4) * t527 + t459;
t467 = t504 * mrSges(6,2) + t511 * t550 - t543;
t453 = -mrSges(5,1) * t476 + mrSges(5,3) * t474 - mrSges(6,1) * t472 + mrSges(6,3) * t471 - pkin(4) * t467 + qJ(5) * t555 + t571 * t505 + t561 * t504 + (-qJ(5) * mrSges(6,2) + t558) * qJDD(4) + (-qJ(5) * t509 + t552) * qJD(4) + (-qJ(5) * t502 * t534 + t554 * t532) * qJD(1);
t454 = mrSges(5,2) * t476 + mrSges(6,2) * t472 - mrSges(5,3) * t473 - mrSges(6,3) * t470 - qJ(5) * t466 - t553 * qJD(4) + t560 * qJDD(4) + t572 * t504 + t561 * t505 - t554 * t550;
t537 = -mrSges(3,1) * t489 - mrSges(4,1) * t479 - mrSges(2,2) * t514 - pkin(2) * t452 - pkin(3) * t465 - pkin(6) * t460 - t534 * t453 - t532 * t454 + qJ(2) * (-t536 * mrSges(3,1) + t541) - pkin(1) * t451 + mrSges(4,2) * t480 + mrSges(3,3) * t488 + mrSges(2,1) * t513 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t457 = -m(3) * g(3) - t458;
t443 = -mrSges(4,1) * t527 + mrSges(4,3) * t480 + t536 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t459 - t569;
t442 = mrSges(4,2) * t527 - mrSges(4,3) * t479 - Ifges(4,5) * qJDD(1) - t536 * Ifges(4,6) - pkin(6) * t459 - t532 * t453 + t534 * t454;
t441 = mrSges(3,2) * t489 - mrSges(2,3) * t513 - qJ(2) * t457 - qJ(3) * t452 + t531 * t442 - t530 * t443 - t559 * t536 + t562 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t440 = mrSges(3,2) * t488 + mrSges(2,3) * t514 - pkin(1) * t457 + pkin(2) * t458 + t564 * g(3) - qJ(3) * t545 + t559 * qJDD(1) - t530 * t442 - t531 * t443 + t562 * t536;
t1 = [-m(1) * g(1) + t546; -m(1) * g(2) + t556; (-m(1) - m(2) - m(3)) * g(3) - t458; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t556 - t533 * t440 + t535 * t441; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t546 + t535 * t440 + t533 * t441; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t537; t537; t451; t458; t569; t467;];
tauJB = t1;
