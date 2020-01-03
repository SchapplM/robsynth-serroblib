% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:15
% EndTime: 2019-12-31 18:47:17
% DurationCPUTime: 1.16s
% Computational Cost: add. (11272->202), mult. (14383->235), div. (0->0), fcn. (5151->6), ass. (0->86)
t558 = Ifges(5,1) + Ifges(6,1);
t546 = Ifges(5,4) - Ifges(6,5);
t545 = Ifges(5,5) + Ifges(6,4);
t557 = Ifges(5,2) + Ifges(6,3);
t543 = Ifges(5,6) - Ifges(6,6);
t556 = Ifges(5,3) + Ifges(6,2);
t508 = -qJD(1) + qJD(3);
t517 = sin(qJ(4));
t520 = cos(qJ(4));
t489 = (-mrSges(6,1) * t520 - mrSges(6,3) * t517) * t508;
t507 = -qJDD(1) + qJDD(3);
t535 = qJD(4) * t508;
t491 = t517 * t507 + t520 * t535;
t524 = qJD(1) ^ 2;
t519 = sin(qJ(1));
t522 = cos(qJ(1));
t501 = -t522 * g(1) - t519 * g(2);
t529 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t501;
t551 = -pkin(1) - pkin(2);
t471 = t551 * t524 + t529;
t500 = t519 * g(1) - t522 * g(2);
t528 = -t524 * qJ(2) + qJDD(2) - t500;
t474 = t551 * qJDD(1) + t528;
t518 = sin(qJ(3));
t521 = cos(qJ(3));
t466 = t521 * t471 + t518 * t474;
t506 = t508 ^ 2;
t463 = -t506 * pkin(3) + t507 * pkin(7) + t466;
t488 = (-pkin(4) * t520 - qJ(5) * t517) * t508;
t523 = qJD(4) ^ 2;
t550 = t520 * g(3);
t458 = -qJDD(4) * pkin(4) - t550 - t523 * qJ(5) + qJDD(5) + (t488 * t508 + t463) * t517;
t540 = t508 * t520;
t498 = mrSges(6,2) * t540 + qJD(4) * mrSges(6,3);
t531 = -m(6) * t458 + qJDD(4) * mrSges(6,1) + qJD(4) * t498;
t541 = t508 * t517;
t454 = t491 * mrSges(6,2) + t489 * t541 - t531;
t460 = t517 * g(3) + t520 * t463;
t457 = -t523 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t488 * t540 + t460;
t459 = -t517 * t463 + t550;
t492 = t520 * t507 - t517 * t535;
t496 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t541;
t532 = m(6) * t457 + qJDD(4) * mrSges(6,3) + qJD(4) * t496 + t489 * t540;
t536 = (t558 * t517 + t546 * t520) * t508 + t545 * qJD(4);
t538 = (-t546 * t517 - t557 * t520) * t508 - t543 * qJD(4);
t555 = -(t538 * t517 + t536 * t520) * t508 + t556 * qJDD(4) + t545 * t491 + t543 * t492 + mrSges(5,1) * t459 - mrSges(6,1) * t458 - mrSges(5,2) * t460 + mrSges(6,3) * t457 - pkin(4) * t454 + qJ(5) * (t492 * mrSges(6,2) + t532);
t552 = -m(3) - m(4);
t549 = -mrSges(2,1) - mrSges(3,1);
t548 = mrSges(5,3) + mrSges(6,2);
t547 = Ifges(3,4) + Ifges(2,5);
t544 = Ifges(2,6) - Ifges(3,6);
t481 = -t524 * pkin(1) + t529;
t490 = (-mrSges(5,1) * t520 + mrSges(5,2) * t517) * t508;
t495 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t541;
t450 = m(5) * t460 - qJDD(4) * mrSges(5,2) - qJD(4) * t495 + t490 * t540 + t548 * t492 + t532;
t497 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t540;
t451 = m(5) * t459 + qJDD(4) * mrSges(5,1) + qJD(4) * t497 + (-t489 - t490) * t541 - t548 * t491 + t531;
t446 = t520 * t450 - t517 * t451;
t443 = m(4) * t466 - t506 * mrSges(4,1) - t507 * mrSges(4,2) + t446;
t465 = -t518 * t471 + t521 * t474;
t462 = -t507 * pkin(3) - t506 * pkin(7) - t465;
t455 = -t492 * pkin(4) - t491 * qJ(5) + (-0.2e1 * qJD(5) * t517 + (pkin(4) * t517 - qJ(5) * t520) * qJD(4)) * t508 + t462;
t452 = m(6) * t455 - t492 * mrSges(6,1) - t491 * mrSges(6,3) - t496 * t541 - t498 * t540;
t448 = -m(5) * t462 + t492 * mrSges(5,1) - t491 * mrSges(5,2) - t495 * t541 + t497 * t540 - t452;
t447 = m(4) * t465 + t507 * mrSges(4,1) - t506 * mrSges(4,2) + t448;
t533 = t521 * t443 - t518 * t447;
t530 = m(3) * t481 + qJDD(1) * mrSges(3,3) + t533;
t434 = m(2) * t501 - qJDD(1) * mrSges(2,2) + t549 * t524 + t530;
t439 = t518 * t443 + t521 * t447;
t487 = -qJDD(1) * pkin(1) + t528;
t438 = m(3) * t487 - qJDD(1) * mrSges(3,1) - t524 * mrSges(3,3) + t439;
t435 = m(2) * t500 + qJDD(1) * mrSges(2,1) - t524 * mrSges(2,2) - t438;
t539 = t519 * t434 + t522 * t435;
t537 = (t545 * t517 + t543 * t520) * t508 + t556 * qJD(4);
t534 = t522 * t434 - t519 * t435;
t445 = t517 * t450 + t520 * t451;
t440 = -mrSges(5,1) * t462 - mrSges(6,1) * t455 + mrSges(6,2) * t457 + mrSges(5,3) * t460 - pkin(4) * t452 + t536 * qJD(4) + t543 * qJDD(4) + t546 * t491 + t557 * t492 - t537 * t541;
t441 = mrSges(5,2) * t462 + mrSges(6,2) * t458 - mrSges(5,3) * t459 - mrSges(6,3) * t455 - qJ(5) * t452 + t538 * qJD(4) + t545 * qJDD(4) + t558 * t491 + t546 * t492 + t537 * t540;
t526 = mrSges(4,1) * t465 - mrSges(4,2) * t466 + Ifges(4,3) * t507 + pkin(3) * t448 + pkin(7) * t446 + t520 * t440 + t517 * t441;
t525 = -mrSges(3,1) * t487 - mrSges(2,2) * t501 - pkin(2) * t439 + qJ(2) * (-t524 * mrSges(3,1) + t530) - pkin(1) * t438 + mrSges(3,3) * t481 + mrSges(2,1) * t500 - t526 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);
t444 = t552 * g(3) - t445;
t430 = -mrSges(4,1) * g(3) + mrSges(4,3) * t466 + t506 * Ifges(4,5) + Ifges(4,6) * t507 - pkin(3) * t445 - t555;
t429 = mrSges(4,2) * g(3) - mrSges(4,3) * t465 + Ifges(4,5) * t507 - t506 * Ifges(4,6) - pkin(7) * t445 - t517 * t440 + t520 * t441;
t428 = mrSges(3,2) * t487 - mrSges(2,3) * t500 - pkin(6) * t439 - qJ(2) * t444 + t521 * t429 - t518 * t430 - t544 * t524 + t547 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t427 = mrSges(2,3) * t501 + mrSges(3,2) * t481 - t518 * t429 - t521 * t430 + pkin(2) * t445 - pkin(6) * t533 - pkin(1) * t444 + t547 * t524 + t544 * qJDD(1) + (pkin(2) * m(4) - t549) * g(3);
t1 = [-m(1) * g(1) + t534; -m(1) * g(2) + t539; (-m(1) - m(2) + t552) * g(3) - t445; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t539 - t519 * t427 + t522 * t428; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t534 + t522 * t427 + t519 * t428; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t525; t525; t438; t526; t555; t454;];
tauJB = t1;
