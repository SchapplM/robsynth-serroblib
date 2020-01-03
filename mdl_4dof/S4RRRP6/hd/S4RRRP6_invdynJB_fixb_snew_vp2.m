% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:19
% EndTime: 2019-12-31 17:18:21
% DurationCPUTime: 1.36s
% Computational Cost: add. (10660->216), mult. (20953->261), div. (0->0), fcn. (12024->6), ass. (0->89)
t568 = Ifges(4,1) + Ifges(5,1);
t562 = Ifges(4,4) + Ifges(5,4);
t561 = Ifges(4,5) + Ifges(5,5);
t567 = Ifges(4,2) + Ifges(5,2);
t560 = Ifges(4,6) + Ifges(5,6);
t566 = Ifges(4,3) + Ifges(5,3);
t531 = sin(qJ(3));
t534 = cos(qJ(3));
t532 = sin(qJ(2));
t553 = qJD(1) * t532;
t516 = t534 * qJD(2) - t531 * t553;
t535 = cos(qJ(2));
t551 = qJD(1) * qJD(2);
t547 = t535 * t551;
t520 = t532 * qJDD(1) + t547;
t492 = t516 * qJD(3) + t531 * qJDD(2) + t534 * t520;
t517 = t531 * qJD(2) + t534 * t553;
t552 = t535 * qJD(1);
t527 = qJD(3) - t552;
t499 = t527 * mrSges(5,1) - t517 * mrSges(5,3);
t533 = sin(qJ(1));
t536 = cos(qJ(1));
t526 = -t536 * g(1) - t533 * g(2);
t538 = qJD(1) ^ 2;
t510 = -t538 * pkin(1) + qJDD(1) * pkin(5) + t526;
t501 = -t535 * g(3) - t532 * t510;
t519 = (-pkin(2) * t535 - pkin(6) * t532) * qJD(1);
t537 = qJD(2) ^ 2;
t477 = -qJDD(2) * pkin(2) - t537 * pkin(6) + t519 * t553 - t501;
t491 = -t517 * qJD(3) + t534 * qJDD(2) - t531 * t520;
t498 = t527 * pkin(3) - t517 * qJ(4);
t514 = t516 ^ 2;
t470 = -t491 * pkin(3) - t514 * qJ(4) + t517 * t498 + qJDD(4) + t477;
t496 = -t527 * mrSges(5,2) + t516 * mrSges(5,3);
t544 = -m(5) * t470 + t491 * mrSges(5,1) + t516 * t496;
t465 = t492 * mrSges(5,2) + t517 * t499 - t544;
t525 = t533 * g(1) - t536 * g(2);
t509 = -qJDD(1) * pkin(1) - t538 * pkin(5) - t525;
t548 = t532 * t551;
t521 = t535 * qJDD(1) - t548;
t475 = (-t520 - t547) * pkin(6) + (-t521 + t548) * pkin(2) + t509;
t502 = -t532 * g(3) + t535 * t510;
t478 = -t537 * pkin(2) + qJDD(2) * pkin(6) + t519 * t552 + t502;
t472 = t531 * t475 + t534 * t478;
t469 = -t514 * pkin(3) + t491 * qJ(4) + 0.2e1 * qJD(4) * t516 - t527 * t498 + t472;
t515 = qJDD(3) - t521;
t494 = -t516 * mrSges(5,1) + t517 * mrSges(5,2);
t549 = m(5) * t469 + t491 * mrSges(5,3) + t516 * t494;
t555 = t562 * t516 + t568 * t517 + t561 * t527;
t557 = -t560 * t516 - t561 * t517 - t566 * t527;
t448 = -mrSges(4,1) * t477 + mrSges(4,3) * t472 - mrSges(5,1) * t470 + mrSges(5,3) * t469 - pkin(3) * t465 + qJ(4) * t549 + (-qJ(4) * t499 + t555) * t527 + t557 * t517 + (-qJ(4) * mrSges(5,2) + t560) * t515 + t562 * t492 + t567 * t491;
t471 = t534 * t475 - t531 * t478;
t467 = -0.2e1 * qJD(4) * t517 + (t516 * t527 - t492) * qJ(4) + (t516 * t517 + t515) * pkin(3) + t471;
t550 = m(5) * t467 + t515 * mrSges(5,1) + t527 * t496;
t464 = -t492 * mrSges(5,3) - t517 * t494 + t550;
t556 = -t567 * t516 - t562 * t517 - t560 * t527;
t451 = mrSges(4,2) * t477 + mrSges(5,2) * t470 - mrSges(4,3) * t471 - mrSges(5,3) * t467 - qJ(4) * t464 + t562 * t491 + t568 * t492 + t561 * t515 - t557 * t516 + t556 * t527;
t495 = -t516 * mrSges(4,1) + t517 * mrSges(4,2);
t497 = -t527 * mrSges(4,2) + t516 * mrSges(4,3);
t458 = m(4) * t471 + t515 * mrSges(4,1) + t527 * t497 + (-t494 - t495) * t517 + (-mrSges(4,3) - mrSges(5,3)) * t492 + t550;
t554 = -t527 * mrSges(4,1) + t517 * mrSges(4,3) - t499;
t563 = -mrSges(4,2) - mrSges(5,2);
t460 = m(4) * t472 + t491 * mrSges(4,3) + t516 * t495 + t563 * t515 + t554 * t527 + t549;
t457 = -t531 * t458 + t534 * t460;
t463 = -m(4) * t477 + t491 * mrSges(4,1) + t563 * t492 + t516 * t497 + t554 * t517 + t544;
t507 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t532 + Ifges(3,2) * t535) * qJD(1);
t508 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t532 + Ifges(3,4) * t535) * qJD(1);
t565 = mrSges(3,1) * t501 - mrSges(3,2) * t502 + Ifges(3,5) * t520 + Ifges(3,6) * t521 + Ifges(3,3) * qJDD(2) + pkin(2) * t463 + pkin(6) * t457 + t534 * t448 + t531 * t451 + (t532 * t507 - t535 * t508) * qJD(1);
t564 = mrSges(4,1) * t471 + mrSges(5,1) * t467 - mrSges(4,2) * t472 - mrSges(5,2) * t469 + pkin(3) * t464 + t560 * t491 + t561 * t492 + t566 * t515 - t555 * t516 - t556 * t517;
t518 = (-mrSges(3,1) * t535 + mrSges(3,2) * t532) * qJD(1);
t523 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t553;
t455 = m(3) * t502 - qJDD(2) * mrSges(3,2) + t521 * mrSges(3,3) - qJD(2) * t523 + t518 * t552 + t457;
t524 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t552;
t462 = m(3) * t501 + qJDD(2) * mrSges(3,1) - t520 * mrSges(3,3) + qJD(2) * t524 - t518 * t553 + t463;
t545 = t535 * t455 - t532 * t462;
t445 = m(2) * t526 - t538 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t545;
t456 = t534 * t458 + t531 * t460;
t540 = -m(3) * t509 + t521 * mrSges(3,1) - t520 * mrSges(3,2) - t523 * t553 + t524 * t552 - t456;
t450 = m(2) * t525 + qJDD(1) * mrSges(2,1) - t538 * mrSges(2,2) + t540;
t558 = t533 * t445 + t536 * t450;
t447 = t532 * t455 + t535 * t462;
t546 = t536 * t445 - t533 * t450;
t506 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t532 + Ifges(3,6) * t535) * qJD(1);
t440 = mrSges(3,2) * t509 - mrSges(3,3) * t501 + Ifges(3,1) * t520 + Ifges(3,4) * t521 + Ifges(3,5) * qJDD(2) - pkin(6) * t456 - qJD(2) * t507 - t531 * t448 + t534 * t451 + t506 * t552;
t442 = -mrSges(3,1) * t509 + mrSges(3,3) * t502 + Ifges(3,4) * t520 + Ifges(3,2) * t521 + Ifges(3,6) * qJDD(2) - pkin(2) * t456 + qJD(2) * t508 - t506 * t553 - t564;
t542 = mrSges(2,1) * t525 - mrSges(2,2) * t526 + Ifges(2,3) * qJDD(1) + pkin(1) * t540 + pkin(5) * t545 + t532 * t440 + t535 * t442;
t438 = mrSges(2,1) * g(3) + mrSges(2,3) * t526 + t538 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t447 - t565;
t437 = -mrSges(2,2) * g(3) - mrSges(2,3) * t525 + Ifges(2,5) * qJDD(1) - t538 * Ifges(2,6) - pkin(5) * t447 + t535 * t440 - t532 * t442;
t1 = [-m(1) * g(1) + t546; -m(1) * g(2) + t558; (-m(1) - m(2)) * g(3) + t447; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t558 + t536 * t437 - t533 * t438; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t546 + t533 * t437 + t536 * t438; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t542; t542; t565; t564; t465;];
tauJB = t1;
