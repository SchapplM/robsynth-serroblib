% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:47
% EndTime: 2019-12-31 18:10:49
% DurationCPUTime: 1.20s
% Computational Cost: add. (7857->223), mult. (15096->268), div. (0->0), fcn. (6511->6), ass. (0->90)
t557 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t543 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t542 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t556 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t541 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t555 = -Ifges(4,3) - Ifges(5,2) - Ifges(6,3);
t554 = 2 * qJD(4);
t553 = pkin(3) + pkin(4);
t524 = qJD(1) ^ 2;
t552 = pkin(4) * t524;
t551 = t524 * pkin(6);
t550 = mrSges(4,3) + mrSges(5,2);
t521 = cos(qJ(3));
t549 = qJ(4) * t521;
t516 = -g(3) + qJDD(2);
t548 = t521 * t516;
t520 = sin(qJ(1));
t522 = cos(qJ(1));
t503 = t520 * g(1) - t522 * g(2);
t484 = qJDD(1) * pkin(1) + t503;
t504 = -t522 * g(1) - t520 * g(2);
t489 = -t524 * pkin(1) + t504;
t517 = sin(pkin(7));
t518 = cos(pkin(7));
t455 = t517 * t484 + t518 * t489;
t453 = -t524 * pkin(2) + qJDD(1) * pkin(6) + t455;
t519 = sin(qJ(3));
t449 = t521 * t453 + t519 * t516;
t487 = (mrSges(6,1) * t521 + mrSges(6,2) * t519) * qJD(1);
t488 = (-mrSges(4,1) * t521 + mrSges(4,2) * t519) * qJD(1);
t544 = qJD(1) * qJD(3);
t491 = t521 * qJDD(1) - t519 * t544;
t546 = qJD(1) * t519;
t498 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t546;
t485 = (-pkin(3) * t521 - qJ(4) * t519) * qJD(1);
t523 = qJD(3) ^ 2;
t545 = qJD(1) * t521;
t446 = -t523 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t554 + t485 * t545 + t449;
t486 = (-mrSges(5,1) * t521 - mrSges(5,3) * t519) * qJD(1);
t499 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t546;
t496 = -qJD(3) * pkin(4) - qJ(5) * t546;
t515 = t521 ^ 2;
t540 = -0.2e1 * qJD(1) * qJD(5);
t442 = -t491 * qJ(5) + qJD(3) * t496 - t515 * t552 + t521 * t540 + t446;
t497 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t546;
t532 = m(6) * t442 + qJDD(3) * mrSges(6,2) - t491 * mrSges(6,3) + qJD(3) * t497;
t527 = m(5) * t446 + qJDD(3) * mrSges(5,3) + qJD(3) * t499 + t486 * t545 + t532;
t435 = m(4) * t449 - qJDD(3) * mrSges(4,2) - qJD(3) * t498 + t550 * t491 + (-t487 + t488) * t545 + t527;
t450 = t519 * t453;
t448 = -t450 + t548;
t490 = t519 * qJDD(1) + t521 * t544;
t501 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t545;
t529 = -t523 * qJ(4) + t485 * t546 + qJDD(4) + t450;
t443 = t519 * t540 - t490 * qJ(5) - t553 * qJDD(3) + (qJ(5) * t544 - t519 * t552 - t516) * t521 + t529;
t500 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t545;
t438 = m(6) * t443 - qJDD(3) * mrSges(6,1) - t490 * mrSges(6,3) - qJD(3) * t500 - t487 * t546;
t447 = -qJDD(3) * pkin(3) + t529 - t548;
t502 = mrSges(5,2) * t545 + qJD(3) * mrSges(5,3);
t526 = -m(5) * t447 + qJDD(3) * mrSges(5,1) + qJD(3) * t502 - t438;
t436 = m(4) * t448 + qJDD(3) * mrSges(4,1) + qJD(3) * t501 - t550 * t490 + (-t486 - t488) * t546 + t526;
t533 = t521 * t435 - t519 * t436;
t428 = m(3) * t455 - t524 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t533;
t454 = t518 * t484 - t517 * t489;
t531 = qJDD(1) * pkin(2) + t454;
t528 = -t490 * qJ(4) - t531;
t444 = -t491 * pkin(3) - t551 + (-0.2e1 * qJD(4) * t519 + (pkin(3) * t519 - t549) * qJD(3)) * qJD(1) + t528;
t440 = qJDD(5) + (-qJ(5) * t515 + pkin(6)) * t524 + t553 * t491 + (qJD(3) * t549 + (-pkin(3) * qJD(3) + t496 + t554) * t519) * qJD(1) - t528;
t530 = -m(6) * t440 - t491 * mrSges(6,1) - t490 * mrSges(6,2) - t497 * t546 - t500 * t545;
t437 = m(5) * t444 - t491 * mrSges(5,1) - t490 * mrSges(5,3) - t499 * t546 - t502 * t545 + t530;
t452 = -t531 - t551;
t525 = -m(4) * t452 + t491 * mrSges(4,1) - t490 * mrSges(4,2) - t498 * t546 + t501 * t545 - t437;
t431 = m(3) * t454 + qJDD(1) * mrSges(3,1) - t524 * mrSges(3,2) + t525;
t424 = t517 * t428 + t518 * t431;
t422 = m(2) * t503 + qJDD(1) * mrSges(2,1) - t524 * mrSges(2,2) + t424;
t534 = t518 * t428 - t517 * t431;
t423 = m(2) * t504 - t524 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t534;
t547 = t522 * t422 + t520 * t423;
t429 = t519 * t435 + t521 * t436;
t539 = t555 * qJD(3) + (-t542 * t519 - t541 * t521) * qJD(1);
t538 = -t541 * qJD(3) + (-t543 * t519 - t556 * t521) * qJD(1);
t537 = t542 * qJD(3) + (t557 * t519 + t543 * t521) * qJD(1);
t536 = m(3) * t516 + t429;
t535 = -t520 * t422 + t522 * t423;
t425 = mrSges(4,2) * t452 + mrSges(5,2) * t447 + mrSges(6,2) * t440 - mrSges(4,3) * t448 - mrSges(5,3) * t444 - mrSges(6,3) * t443 - qJ(4) * t437 - qJ(5) * t438 + t538 * qJD(3) + t542 * qJDD(3) + t557 * t490 + t543 * t491 - t539 * t545;
t418 = -mrSges(4,1) * t452 + mrSges(4,3) * t449 - mrSges(5,1) * t444 + mrSges(5,2) * t446 + mrSges(6,1) * t440 - mrSges(6,3) * t442 - pkin(4) * t530 - qJ(5) * t532 - pkin(3) * t437 + t556 * t491 + t543 * t490 + t541 * qJDD(3) + t537 * qJD(3) + (qJ(5) * t487 * t521 + t539 * t519) * qJD(1);
t417 = -pkin(3) * t526 - qJ(4) * t527 + pkin(4) * t438 - mrSges(6,2) * t442 + mrSges(6,1) * t443 - mrSges(5,3) * t446 + mrSges(5,1) * t447 - mrSges(4,1) * t448 + mrSges(4,2) * t449 + mrSges(3,3) * t455 + t524 * Ifges(3,5) - pkin(2) * t429 - mrSges(3,1) * t516 + Ifges(3,6) * qJDD(1) + (-qJ(4) * mrSges(5,2) - t541) * t491 + (pkin(3) * mrSges(5,2) - t542) * t490 + t555 * qJDD(3) + ((qJ(4) * t487 + t537) * t521 + (pkin(3) * t486 + t538) * t519) * qJD(1);
t416 = mrSges(3,2) * t516 - mrSges(3,3) * t454 + Ifges(3,5) * qJDD(1) - t524 * Ifges(3,6) - pkin(6) * t429 - t519 * t418 + t521 * t425;
t415 = -mrSges(2,2) * g(3) - mrSges(2,3) * t503 + Ifges(2,5) * qJDD(1) - t524 * Ifges(2,6) - qJ(2) * t424 + t518 * t416 - t517 * t417;
t414 = mrSges(2,1) * g(3) + mrSges(2,3) * t504 + t524 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t536 + qJ(2) * t534 + t517 * t416 + t518 * t417;
t1 = [-m(1) * g(1) + t535; -m(1) * g(2) + t547; (-m(1) - m(2)) * g(3) + t536; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t547 - t520 * t414 + t522 * t415; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t535 + t522 * t414 + t520 * t415; pkin(1) * t424 + mrSges(2,1) * t503 - mrSges(2,2) * t504 + t519 * t425 + t521 * t418 + pkin(2) * t525 + pkin(6) * t533 + mrSges(3,1) * t454 - mrSges(3,2) * t455 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
