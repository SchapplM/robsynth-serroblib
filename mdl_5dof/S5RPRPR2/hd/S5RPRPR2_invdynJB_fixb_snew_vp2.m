% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:56
% EndTime: 2020-01-03 11:34:00
% DurationCPUTime: 3.22s
% Computational Cost: add. (42687->206), mult. (61438->259), div. (0->0), fcn. (35545->10), ass. (0->99)
t574 = qJD(1) + qJD(3);
t570 = t574 ^ 2;
t580 = cos(pkin(9));
t614 = pkin(4) * t580;
t578 = sin(pkin(9));
t613 = mrSges(5,2) * t578;
t573 = t580 ^ 2;
t612 = t570 * t573;
t571 = qJDD(1) + qJDD(3);
t611 = t571 * t580;
t598 = Ifges(5,5) * t578 + Ifges(5,6) * t580;
t610 = t570 * t598;
t584 = sin(qJ(1));
t587 = cos(qJ(1));
t559 = -t584 * g(2) + t587 * g(3);
t588 = qJD(1) ^ 2;
t560 = -t587 * g(2) - t584 * g(3);
t554 = qJDD(1) * pkin(1) + t560;
t555 = -t588 * pkin(1) + t559;
t579 = sin(pkin(8));
t581 = cos(pkin(8));
t543 = t581 * t554 - t579 * t555;
t540 = qJDD(1) * pkin(2) + t543;
t544 = t579 * t554 + t581 * t555;
t541 = -t588 * pkin(2) + t544;
t583 = sin(qJ(3));
t586 = cos(qJ(3));
t528 = t583 * t540 + t586 * t541;
t525 = -t570 * pkin(3) + t571 * qJ(4) + t528;
t577 = -g(1) + qJDD(2);
t607 = qJD(4) * t574;
t608 = t580 * t577 - 0.2e1 * t578 * t607;
t518 = (-pkin(7) * t571 + t570 * t614 - t525) * t578 + t608;
t522 = t578 * t577 + (t525 + 0.2e1 * t607) * t580;
t519 = -pkin(4) * t612 + pkin(7) * t611 + t522;
t582 = sin(qJ(5));
t585 = cos(qJ(5));
t516 = t585 * t518 - t582 * t519;
t594 = -t578 * t582 + t580 * t585;
t547 = t594 * t574;
t595 = t578 * t585 + t580 * t582;
t548 = t595 * t574;
t534 = -t547 * mrSges(6,1) + t548 * mrSges(6,2);
t536 = t547 * qJD(5) + t595 * t571;
t545 = -qJD(5) * mrSges(6,2) + t547 * mrSges(6,3);
t514 = m(6) * t516 + qJDD(5) * mrSges(6,1) - t536 * mrSges(6,3) + qJD(5) * t545 - t548 * t534;
t517 = t582 * t518 + t585 * t519;
t535 = -t548 * qJD(5) + t594 * t571;
t546 = qJD(5) * mrSges(6,1) - t548 * mrSges(6,3);
t515 = m(6) * t517 - qJDD(5) * mrSges(6,2) + t535 * mrSges(6,3) - qJD(5) * t546 + t547 * t534;
t504 = t585 * t514 + t582 * t515;
t521 = -t578 * t525 + t608;
t596 = mrSges(5,3) * t571 + (-mrSges(5,1) * t580 + t613) * t570;
t502 = m(5) * t521 - t596 * t578 + t504;
t601 = -t582 * t514 + t585 * t515;
t503 = m(5) * t522 + t596 * t580 + t601;
t602 = -t578 * t502 + t580 * t503;
t495 = m(4) * t528 - t570 * mrSges(4,1) - t571 * mrSges(4,2) + t602;
t527 = t586 * t540 - t583 * t541;
t597 = qJDD(4) - t527;
t524 = -t571 * pkin(3) - t570 * qJ(4) + t597;
t572 = t578 ^ 2;
t520 = (-pkin(3) - t614) * t571 + (-qJ(4) + (-t572 - t573) * pkin(7)) * t570 + t597;
t592 = m(6) * t520 - t535 * mrSges(6,1) + t536 * mrSges(6,2) - t547 * t545 + t548 * t546;
t591 = -m(5) * t524 + mrSges(5,1) * t611 - t592 + (t570 * t572 + t612) * mrSges(5,3);
t508 = m(4) * t527 - t570 * mrSges(4,2) + (mrSges(4,1) - t613) * t571 + t591;
t490 = t583 * t495 + t586 * t508;
t485 = m(3) * t543 + qJDD(1) * mrSges(3,1) - t588 * mrSges(3,2) + t490;
t603 = t586 * t495 - t583 * t508;
t486 = m(3) * t544 - t588 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t603;
t604 = -t579 * t485 + t581 * t486;
t477 = m(2) * t559 - t588 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t604;
t480 = t581 * t485 + t579 * t486;
t478 = m(2) * t560 + qJDD(1) * mrSges(2,1) - t588 * mrSges(2,2) + t480;
t609 = t584 * t477 + t587 * t478;
t498 = t580 * t502 + t578 * t503;
t606 = m(4) * t577 + t498;
t605 = -t587 * t477 + t584 * t478;
t496 = m(3) * t577 + t606;
t600 = Ifges(5,1) * t578 + Ifges(5,4) * t580;
t599 = Ifges(5,4) * t578 + Ifges(5,2) * t580;
t529 = Ifges(6,5) * t548 + Ifges(6,6) * t547 + Ifges(6,3) * qJD(5);
t531 = Ifges(6,1) * t548 + Ifges(6,4) * t547 + Ifges(6,5) * qJD(5);
t505 = -mrSges(6,1) * t520 + mrSges(6,3) * t517 + Ifges(6,4) * t536 + Ifges(6,2) * t535 + Ifges(6,6) * qJDD(5) + qJD(5) * t531 - t548 * t529;
t530 = Ifges(6,4) * t548 + Ifges(6,2) * t547 + Ifges(6,6) * qJD(5);
t506 = mrSges(6,2) * t520 - mrSges(6,3) * t516 + Ifges(6,1) * t536 + Ifges(6,4) * t535 + Ifges(6,5) * qJDD(5) - qJD(5) * t530 + t547 * t529;
t488 = -mrSges(5,1) * t524 + mrSges(5,3) * t522 - pkin(4) * t592 + pkin(7) * t601 + t585 * t505 + t582 * t506 + t599 * t571 - t578 * t610;
t492 = mrSges(5,2) * t524 - mrSges(5,3) * t521 - pkin(7) * t504 - t582 * t505 + t585 * t506 + t600 * t571 + t580 * t610;
t510 = t571 * t613 - t591;
t593 = mrSges(4,1) * t527 - mrSges(4,2) * t528 + Ifges(4,3) * t571 - pkin(3) * t510 + qJ(4) * t602 + t580 * t488 + t578 * t492;
t590 = mrSges(6,1) * t516 - mrSges(6,2) * t517 + Ifges(6,5) * t536 + Ifges(6,6) * t535 + Ifges(6,3) * qJDD(5) + t548 * t530 - t547 * t531;
t589 = mrSges(2,1) * t560 + mrSges(3,1) * t543 - mrSges(2,2) * t559 - mrSges(3,2) * t544 + pkin(1) * t480 + pkin(2) * t490 + t593 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t481 = -mrSges(4,1) * t577 - mrSges(5,1) * t521 + mrSges(5,2) * t522 + mrSges(4,3) * t528 - pkin(3) * t498 - pkin(4) * t504 + (Ifges(4,6) - t598) * t571 - t590 + (-t578 * t599 + t580 * t600 + Ifges(4,5)) * t570;
t473 = mrSges(4,2) * t577 - mrSges(4,3) * t527 + Ifges(4,5) * t571 - t570 * Ifges(4,6) - qJ(4) * t498 - t578 * t488 + t580 * t492;
t472 = mrSges(3,2) * t577 - mrSges(3,3) * t543 + Ifges(3,5) * qJDD(1) - t588 * Ifges(3,6) - pkin(6) * t490 + t586 * t473 - t583 * t481;
t471 = -mrSges(3,1) * t577 + mrSges(3,3) * t544 + t588 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t606 + pkin(6) * t603 + t583 * t473 + t586 * t481;
t470 = -mrSges(2,2) * g(1) - mrSges(2,3) * t560 + Ifges(2,5) * qJDD(1) - t588 * Ifges(2,6) - qJ(2) * t480 - t579 * t471 + t581 * t472;
t469 = mrSges(2,1) * g(1) + mrSges(2,3) * t559 + t588 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t496 + qJ(2) * t604 + t581 * t471 + t579 * t472;
t1 = [(-m(1) - m(2)) * g(1) + t496; -m(1) * g(2) + t609; -m(1) * g(3) + t605; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t589; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t605 + t587 * t469 + t584 * t470; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t609 + t584 * t469 - t587 * t470; t589; t496; t593; t510; t590;];
tauJB = t1;
