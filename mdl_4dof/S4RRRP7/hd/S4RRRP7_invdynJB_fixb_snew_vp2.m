% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:19
% EndTime: 2019-12-31 17:20:21
% DurationCPUTime: 1.35s
% Computational Cost: add. (10382->214), mult. (20203->261), div. (0->0), fcn. (11499->6), ass. (0->88)
t563 = Ifges(4,1) + Ifges(5,1);
t556 = Ifges(4,4) - Ifges(5,5);
t555 = -Ifges(4,5) - Ifges(5,4);
t562 = Ifges(4,2) + Ifges(5,3);
t554 = Ifges(4,6) - Ifges(5,6);
t561 = -Ifges(4,3) - Ifges(5,2);
t529 = sin(qJ(1));
t531 = cos(qJ(1));
t519 = -t531 * g(1) - t529 * g(2);
t533 = qJD(1) ^ 2;
t503 = -t533 * pkin(1) + qJDD(1) * pkin(5) + t519;
t528 = sin(qJ(2));
t530 = cos(qJ(2));
t494 = -t530 * g(3) - t528 * t503;
t512 = (-pkin(2) * t530 - pkin(6) * t528) * qJD(1);
t532 = qJD(2) ^ 2;
t547 = qJD(1) * t528;
t471 = -qJDD(2) * pkin(2) - t532 * pkin(6) + t512 * t547 - t494;
t527 = sin(qJ(3));
t558 = cos(qJ(3));
t510 = t527 * qJD(2) + t558 * t547;
t545 = qJD(1) * qJD(2);
t542 = t530 * t545;
t513 = t528 * qJDD(1) + t542;
t483 = t510 * qJD(3) - t558 * qJDD(2) + t527 * t513;
t509 = -t558 * qJD(2) + t527 * t547;
t484 = -t509 * qJD(3) + t527 * qJDD(2) + t558 * t513;
t546 = t530 * qJD(1);
t521 = qJD(3) - t546;
t464 = -0.2e1 * qJD(4) * t510 + (t509 * t521 - t484) * qJ(4) + (t510 * t521 + t483) * pkin(3) + t471;
t492 = -t521 * mrSges(5,1) + t510 * mrSges(5,2);
t493 = -t509 * mrSges(5,2) + t521 * mrSges(5,3);
t459 = m(5) * t464 + t483 * mrSges(5,1) - t484 * mrSges(5,3) - t510 * t492 + t509 * t493;
t518 = t529 * g(1) - t531 * g(2);
t502 = -qJDD(1) * pkin(1) - t533 * pkin(5) - t518;
t543 = t528 * t545;
t514 = t530 * qJDD(1) - t543;
t469 = (-t513 - t542) * pkin(6) + (-t514 + t543) * pkin(2) + t502;
t495 = -t528 * g(3) + t530 * t503;
t472 = -t532 * pkin(2) + qJDD(2) * pkin(6) + t512 * t546 + t495;
t467 = t527 * t469 + t558 * t472;
t487 = t509 * pkin(3) - t510 * qJ(4);
t508 = qJDD(3) - t514;
t520 = t521 ^ 2;
t463 = -t520 * pkin(3) + t508 * qJ(4) + 0.2e1 * qJD(4) * t521 - t509 * t487 + t467;
t549 = t556 * t509 - t563 * t510 + t555 * t521;
t551 = t554 * t509 + t555 * t510 + t561 * t521;
t443 = -mrSges(4,1) * t471 - mrSges(5,1) * t464 + mrSges(5,2) * t463 + mrSges(4,3) * t467 - pkin(3) * t459 - t562 * t483 + t556 * t484 + t554 * t508 + t551 * t510 - t549 * t521;
t466 = t558 * t469 - t527 * t472;
t465 = -t508 * pkin(3) - t520 * qJ(4) + t510 * t487 + qJDD(4) - t466;
t550 = t562 * t509 - t556 * t510 - t554 * t521;
t444 = mrSges(4,2) * t471 + mrSges(5,2) * t465 - mrSges(4,3) * t466 - mrSges(5,3) * t464 - qJ(4) * t459 - t556 * t483 + t563 * t484 - t555 * t508 + t551 * t509 + t550 * t521;
t491 = t521 * mrSges(4,1) - t510 * mrSges(4,3);
t544 = m(5) * t463 + t508 * mrSges(5,3) + t521 * t492;
t488 = t509 * mrSges(5,1) - t510 * mrSges(5,3);
t548 = -t509 * mrSges(4,1) - t510 * mrSges(4,2) - t488;
t557 = -mrSges(4,3) - mrSges(5,2);
t456 = m(4) * t467 - t508 * mrSges(4,2) + t557 * t483 - t521 * t491 + t548 * t509 + t544;
t490 = -t521 * mrSges(4,2) - t509 * mrSges(4,3);
t539 = -m(5) * t465 + t508 * mrSges(5,1) + t521 * t493;
t457 = m(4) * t466 + t508 * mrSges(4,1) + t557 * t484 + t521 * t490 + t548 * t510 + t539;
t452 = t558 * t456 - t527 * t457;
t458 = -m(4) * t471 - t483 * mrSges(4,1) - t484 * mrSges(4,2) - t509 * t490 - t510 * t491 - t459;
t500 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t528 + Ifges(3,2) * t530) * qJD(1);
t501 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t528 + Ifges(3,4) * t530) * qJD(1);
t560 = mrSges(3,1) * t494 - mrSges(3,2) * t495 + Ifges(3,5) * t513 + Ifges(3,6) * t514 + Ifges(3,3) * qJDD(2) + pkin(2) * t458 + pkin(6) * t452 + (t500 * t528 - t501 * t530) * qJD(1) + t558 * t443 + t527 * t444;
t461 = t484 * mrSges(5,2) + t510 * t488 - t539;
t559 = -t554 * t483 - t555 * t484 - t561 * t508 - t549 * t509 - t550 * t510 + mrSges(4,1) * t466 - mrSges(5,1) * t465 - mrSges(4,2) * t467 + mrSges(5,3) * t463 - pkin(3) * t461 + qJ(4) * (-t483 * mrSges(5,2) - t509 * t488 + t544);
t511 = (-mrSges(3,1) * t530 + mrSges(3,2) * t528) * qJD(1);
t516 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t547;
t450 = m(3) * t495 - qJDD(2) * mrSges(3,2) + t514 * mrSges(3,3) - qJD(2) * t516 + t511 * t546 + t452;
t517 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t546;
t454 = m(3) * t494 + qJDD(2) * mrSges(3,1) - t513 * mrSges(3,3) + qJD(2) * t517 - t511 * t547 + t458;
t540 = t530 * t450 - t528 * t454;
t440 = m(2) * t519 - t533 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t540;
t451 = t527 * t456 + t558 * t457;
t535 = -m(3) * t502 + t514 * mrSges(3,1) - t513 * mrSges(3,2) - t516 * t547 + t517 * t546 - t451;
t446 = m(2) * t518 + qJDD(1) * mrSges(2,1) - t533 * mrSges(2,2) + t535;
t552 = t529 * t440 + t531 * t446;
t442 = t528 * t450 + t530 * t454;
t541 = t531 * t440 - t529 * t446;
t499 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t528 + Ifges(3,6) * t530) * qJD(1);
t435 = mrSges(3,2) * t502 - mrSges(3,3) * t494 + Ifges(3,1) * t513 + Ifges(3,4) * t514 + Ifges(3,5) * qJDD(2) - pkin(6) * t451 - qJD(2) * t500 - t527 * t443 + t558 * t444 + t499 * t546;
t437 = -mrSges(3,1) * t502 + mrSges(3,3) * t495 + Ifges(3,4) * t513 + Ifges(3,2) * t514 + Ifges(3,6) * qJDD(2) - pkin(2) * t451 + qJD(2) * t501 - t499 * t547 - t559;
t537 = mrSges(2,1) * t518 - mrSges(2,2) * t519 + Ifges(2,3) * qJDD(1) + pkin(1) * t535 + pkin(5) * t540 + t528 * t435 + t530 * t437;
t433 = mrSges(2,1) * g(3) + mrSges(2,3) * t519 + t533 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t442 - t560;
t432 = -mrSges(2,2) * g(3) - mrSges(2,3) * t518 + Ifges(2,5) * qJDD(1) - t533 * Ifges(2,6) - pkin(5) * t442 + t530 * t435 - t528 * t437;
t1 = [-m(1) * g(1) + t541; -m(1) * g(2) + t552; (-m(1) - m(2)) * g(3) + t442; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t552 + t531 * t432 - t529 * t433; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t541 + t529 * t432 + t531 * t433; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t537; t537; t560; t559; t461;];
tauJB = t1;
