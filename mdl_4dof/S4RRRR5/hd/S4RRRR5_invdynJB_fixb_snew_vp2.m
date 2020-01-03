% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:36
% EndTime: 2019-12-31 17:27:39
% DurationCPUTime: 2.44s
% Computational Cost: add. (26103->239), mult. (51965->303), div. (0->0), fcn. (33070->8), ass. (0->98)
t541 = sin(qJ(1));
t545 = cos(qJ(1));
t530 = t541 * g(1) - t545 * g(2);
t547 = qJD(1) ^ 2;
t514 = -qJDD(1) * pkin(1) - t547 * pkin(5) - t530;
t540 = sin(qJ(2));
t544 = cos(qJ(2));
t559 = qJD(1) * qJD(2);
t558 = t544 * t559;
t525 = t540 * qJDD(1) + t558;
t534 = t540 * t559;
t526 = t544 * qJDD(1) - t534;
t486 = (-t525 - t558) * pkin(6) + (-t526 + t534) * pkin(2) + t514;
t531 = -t545 * g(1) - t541 * g(2);
t515 = -t547 * pkin(1) + qJDD(1) * pkin(5) + t531;
t505 = -g(3) * t540 + t544 * t515;
t524 = (-pkin(2) * t544 - pkin(6) * t540) * qJD(1);
t546 = qJD(2) ^ 2;
t560 = t544 * qJD(1);
t489 = -pkin(2) * t546 + qJDD(2) * pkin(6) + t524 * t560 + t505;
t539 = sin(qJ(3));
t543 = cos(qJ(3));
t477 = t543 * t486 - t539 * t489;
t561 = qJD(1) * t540;
t521 = t543 * qJD(2) - t539 * t561;
t498 = qJD(3) * t521 + qJDD(2) * t539 + t525 * t543;
t520 = qJDD(3) - t526;
t522 = t539 * qJD(2) + t543 * t561;
t533 = qJD(3) - t560;
t468 = (t521 * t533 - t498) * pkin(7) + (t521 * t522 + t520) * pkin(3) + t477;
t478 = t539 * t486 + t543 * t489;
t497 = -qJD(3) * t522 + qJDD(2) * t543 - t525 * t539;
t506 = pkin(3) * t533 - pkin(7) * t522;
t519 = t521 ^ 2;
t469 = -pkin(3) * t519 + pkin(7) * t497 - t506 * t533 + t478;
t538 = sin(qJ(4));
t542 = cos(qJ(4));
t467 = t468 * t538 + t469 * t542;
t504 = -t544 * g(3) - t540 * t515;
t488 = -qJDD(2) * pkin(2) - pkin(6) * t546 + t524 * t561 - t504;
t470 = -pkin(3) * t497 - pkin(7) * t519 + t506 * t522 + t488;
t500 = t521 * t538 + t522 * t542;
t475 = -qJD(4) * t500 + t497 * t542 - t498 * t538;
t499 = t521 * t542 - t522 * t538;
t476 = qJD(4) * t499 + t497 * t538 + t498 * t542;
t532 = qJD(4) + t533;
t479 = Ifges(5,5) * t500 + Ifges(5,6) * t499 + Ifges(5,3) * t532;
t481 = Ifges(5,1) * t500 + Ifges(5,4) * t499 + Ifges(5,5) * t532;
t516 = qJDD(4) + t520;
t455 = -mrSges(5,1) * t470 + mrSges(5,3) * t467 + Ifges(5,4) * t476 + Ifges(5,2) * t475 + Ifges(5,6) * t516 - t479 * t500 + t481 * t532;
t466 = t468 * t542 - t469 * t538;
t480 = Ifges(5,4) * t500 + Ifges(5,2) * t499 + Ifges(5,6) * t532;
t456 = mrSges(5,2) * t470 - mrSges(5,3) * t466 + Ifges(5,1) * t476 + Ifges(5,4) * t475 + Ifges(5,5) * t516 + t479 * t499 - t480 * t532;
t492 = Ifges(4,5) * t522 + Ifges(4,6) * t521 + Ifges(4,3) * t533;
t494 = Ifges(4,1) * t522 + Ifges(4,4) * t521 + Ifges(4,5) * t533;
t490 = -mrSges(5,2) * t532 + mrSges(5,3) * t499;
t491 = mrSges(5,1) * t532 - mrSges(5,3) * t500;
t552 = m(5) * t470 - t475 * mrSges(5,1) + t476 * mrSges(5,2) - t499 * t490 + t491 * t500;
t483 = -mrSges(5,1) * t499 + mrSges(5,2) * t500;
t463 = m(5) * t466 + mrSges(5,1) * t516 - t476 * mrSges(5,3) - t483 * t500 + t490 * t532;
t464 = m(5) * t467 - mrSges(5,2) * t516 + t475 * mrSges(5,3) + t483 * t499 - t491 * t532;
t555 = -t463 * t538 + t542 * t464;
t436 = -mrSges(4,1) * t488 + mrSges(4,3) * t478 + Ifges(4,4) * t498 + Ifges(4,2) * t497 + Ifges(4,6) * t520 - pkin(3) * t552 + pkin(7) * t555 + t542 * t455 + t538 * t456 - t522 * t492 + t533 * t494;
t454 = t542 * t463 + t538 * t464;
t493 = Ifges(4,4) * t522 + Ifges(4,2) * t521 + Ifges(4,6) * t533;
t442 = mrSges(4,2) * t488 - mrSges(4,3) * t477 + Ifges(4,1) * t498 + Ifges(4,4) * t497 + Ifges(4,5) * t520 - pkin(7) * t454 - t455 * t538 + t456 * t542 + t492 * t521 - t493 * t533;
t501 = -mrSges(4,1) * t521 + mrSges(4,2) * t522;
t502 = -mrSges(4,2) * t533 + mrSges(4,3) * t521;
t452 = m(4) * t477 + mrSges(4,1) * t520 - mrSges(4,3) * t498 - t501 * t522 + t502 * t533 + t454;
t503 = mrSges(4,1) * t533 - mrSges(4,3) * t522;
t453 = m(4) * t478 - mrSges(4,2) * t520 + mrSges(4,3) * t497 + t501 * t521 - t503 * t533 + t555;
t450 = -t452 * t539 + t543 * t453;
t459 = -m(4) * t488 + t497 * mrSges(4,1) - mrSges(4,2) * t498 + t521 * t502 - t503 * t522 - t552;
t512 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t540 + Ifges(3,2) * t544) * qJD(1);
t513 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t540 + Ifges(3,4) * t544) * qJD(1);
t563 = mrSges(3,1) * t504 - mrSges(3,2) * t505 + Ifges(3,5) * t525 + Ifges(3,6) * t526 + Ifges(3,3) * qJDD(2) + pkin(2) * t459 + pkin(6) * t450 + t543 * t436 + t539 * t442 + (t512 * t540 - t513 * t544) * qJD(1);
t523 = (-mrSges(3,1) * t544 + mrSges(3,2) * t540) * qJD(1);
t528 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t561;
t448 = m(3) * t505 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t526 - qJD(2) * t528 + t523 * t560 + t450;
t529 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t560;
t458 = m(3) * t504 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t525 + qJD(2) * t529 - t523 * t561 + t459;
t556 = t544 * t448 - t458 * t540;
t439 = m(2) * t531 - mrSges(2,1) * t547 - qJDD(1) * mrSges(2,2) + t556;
t449 = t452 * t543 + t453 * t539;
t550 = -m(3) * t514 + t526 * mrSges(3,1) - mrSges(3,2) * t525 - t528 * t561 + t529 * t560 - t449;
t444 = m(2) * t530 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t547 + t550;
t562 = t541 * t439 + t545 * t444;
t441 = t540 * t448 + t544 * t458;
t557 = t545 * t439 - t444 * t541;
t511 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t540 + Ifges(3,6) * t544) * qJD(1);
t433 = mrSges(3,2) * t514 - mrSges(3,3) * t504 + Ifges(3,1) * t525 + Ifges(3,4) * t526 + Ifges(3,5) * qJDD(2) - pkin(6) * t449 - qJD(2) * t512 - t436 * t539 + t442 * t543 + t511 * t560;
t551 = -mrSges(5,1) * t466 + mrSges(5,2) * t467 - Ifges(5,5) * t476 - Ifges(5,6) * t475 - Ifges(5,3) * t516 - t500 * t480 + t499 * t481;
t548 = mrSges(4,1) * t477 - mrSges(4,2) * t478 + Ifges(4,5) * t498 + Ifges(4,6) * t497 + Ifges(4,3) * t520 + pkin(3) * t454 + t522 * t493 - t521 * t494 - t551;
t435 = -mrSges(3,1) * t514 + mrSges(3,3) * t505 + Ifges(3,4) * t525 + Ifges(3,2) * t526 + Ifges(3,6) * qJDD(2) - pkin(2) * t449 + qJD(2) * t513 - t511 * t561 - t548;
t553 = mrSges(2,1) * t530 - mrSges(2,2) * t531 + Ifges(2,3) * qJDD(1) + pkin(1) * t550 + pkin(5) * t556 + t540 * t433 + t544 * t435;
t431 = mrSges(2,1) * g(3) + mrSges(2,3) * t531 + t547 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t441 - t563;
t430 = -mrSges(2,2) * g(3) - mrSges(2,3) * t530 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t547 - pkin(5) * t441 + t433 * t544 - t435 * t540;
t1 = [-m(1) * g(1) + t557; -m(1) * g(2) + t562; (-m(1) - m(2)) * g(3) + t441; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t562 + t545 * t430 - t541 * t431; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t557 + t541 * t430 + t545 * t431; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t553; t553; t563; t548; -t551;];
tauJB = t1;
