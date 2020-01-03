% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR16_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:52
% EndTime: 2019-12-31 18:38:54
% DurationCPUTime: 1.34s
% Computational Cost: add. (7945->248), mult. (15474->288), div. (0->0), fcn. (7213->6), ass. (0->105)
t547 = Ifges(4,1) + Ifges(5,2);
t535 = Ifges(4,4) + Ifges(5,6);
t533 = (Ifges(4,5) - Ifges(5,4));
t546 = Ifges(4,2) + Ifges(5,3);
t531 = (Ifges(4,6) - Ifges(5,5));
t545 = (-Ifges(4,3) - Ifges(5,1));
t496 = sin(qJ(1));
t499 = cos(qJ(1));
t480 = -g(1) * t499 - g(2) * t496;
t512 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t480;
t501 = qJD(1) ^ 2;
t541 = (-pkin(1) - pkin(6));
t522 = t541 * t501;
t449 = t522 + t512;
t495 = sin(qJ(3));
t498 = cos(qJ(3));
t524 = qJD(1) * qJD(3);
t521 = t498 * t524;
t471 = qJDD(1) * t495 + t521;
t485 = t495 * t524;
t472 = qJDD(1) * t498 - t485;
t526 = qJD(1) * t495;
t474 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t526;
t525 = qJD(1) * t498;
t475 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t525;
t542 = -2 * qJD(4);
t505 = pkin(3) * t521 + t525 * t542 + t512 + (-t472 + t485) * qJ(4);
t432 = pkin(3) * t471 + t505 + t522;
t476 = mrSges(5,1) * t526 - (qJD(3) * mrSges(5,3));
t477 = mrSges(5,1) * t525 + (qJD(3) * mrSges(5,2));
t478 = pkin(4) * t525 - (qJD(3) * pkin(7));
t491 = t495 ^ 2;
t540 = pkin(3) + pkin(7);
t427 = -t478 * t525 + t540 * t471 + (-pkin(4) * t491 + t541) * t501 + t505;
t468 = (pkin(3) * t495 - qJ(4) * t498) * qJD(1);
t500 = qJD(3) ^ 2;
t479 = g(1) * t496 - t499 * g(2);
t511 = -qJ(2) * t501 + qJDD(2) - t479;
t450 = t541 * qJDD(1) + t511;
t530 = t450 * t498;
t510 = -qJ(4) * t500 + t468 * t525 + qJDD(4) - t530;
t539 = pkin(7) * t501;
t430 = pkin(4) * t472 - t540 * qJDD(3) + (pkin(4) * t524 + t498 * t539 - g(3)) * t495 + t510;
t494 = sin(qJ(5));
t497 = cos(qJ(5));
t425 = -t427 * t494 + t430 * t497;
t466 = -qJD(3) * t494 + t497 * t526;
t441 = qJD(5) * t466 + qJDD(3) * t497 + t471 * t494;
t467 = qJD(3) * t497 + t494 * t526;
t442 = -mrSges(6,1) * t466 + mrSges(6,2) * t467;
t483 = qJD(5) + t525;
t445 = -mrSges(6,2) * t483 + mrSges(6,3) * t466;
t465 = qJDD(5) + t472;
t423 = m(6) * t425 + mrSges(6,1) * t465 - mrSges(6,3) * t441 - t442 * t467 + t445 * t483;
t426 = t427 * t497 + t430 * t494;
t440 = -qJD(5) * t467 - qJDD(3) * t494 + t471 * t497;
t446 = mrSges(6,1) * t483 - mrSges(6,3) * t467;
t424 = m(6) * t426 - mrSges(6,2) * t465 + mrSges(6,3) * t440 + t442 * t466 - t446 * t483;
t518 = -t423 * t494 + t497 * t424;
t503 = m(5) * t432 - t472 * mrSges(5,3) - (t476 * t495 + t477 * t498) * qJD(1) + t518;
t536 = mrSges(4,1) - mrSges(5,2);
t544 = -m(4) * t449 - t472 * mrSges(4,2) - t536 * t471 - t474 * t526 - t475 * t525 - t503;
t538 = g(3) * t495;
t537 = mrSges(2,1) - mrSges(3,2);
t534 = Ifges(2,5) - Ifges(3,4);
t532 = (-Ifges(2,6) + Ifges(3,5));
t443 = t530 + t538;
t417 = t423 * t497 + t424 * t494;
t434 = -qJDD(3) * pkin(3) + t510 - t538;
t507 = -m(5) * t434 - t472 * mrSges(5,1) - t417;
t469 = (-mrSges(5,2) * t495 - mrSges(5,3) * t498) * qJD(1);
t516 = qJD(1) * (-t469 - (mrSges(4,1) * t495 + mrSges(4,2) * t498) * qJD(1));
t415 = m(4) * t443 - mrSges(4,3) * t472 + t536 * qJDD(3) + (t474 - t476) * qJD(3) + t498 * t516 + t507;
t444 = -g(3) * t498 + t495 * t450;
t506 = -pkin(3) * t500 + qJDD(3) * qJ(4) - t468 * t526 + t444;
t433 = (qJD(3) * t542) - t506;
t429 = -t491 * t539 - pkin(4) * t471 + ((2 * qJD(4)) + t478) * qJD(3) + t506;
t509 = -m(6) * t429 + mrSges(6,1) * t440 - t441 * mrSges(6,2) + t445 * t466 - t467 * t446;
t504 = -m(5) * t433 + qJDD(3) * mrSges(5,3) + qJD(3) * t477 - t509;
t421 = m(4) * t444 - qJDD(3) * mrSges(4,2) - qJD(3) * t475 + (-mrSges(4,3) - mrSges(5,1)) * t471 + t495 * t516 + t504;
t411 = t415 * t498 + t421 * t495;
t452 = -qJDD(1) * pkin(1) + t511;
t508 = -m(3) * t452 + (t501 * mrSges(3,3)) - t411;
t409 = m(2) * t479 - (mrSges(2,2) * t501) + t537 * qJDD(1) + t508;
t451 = pkin(1) * t501 - t512;
t502 = -m(3) * t451 + (t501 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t544;
t414 = m(2) * t480 - (mrSges(2,1) * t501) - qJDD(1) * mrSges(2,2) + t502;
t529 = t499 * t409 + t496 * t414;
t528 = -(t531 * qJD(3)) + (t546 * t495 - t535 * t498) * qJD(1);
t527 = (t533 * qJD(3)) + (-t535 * t495 + t547 * t498) * qJD(1);
t520 = -t409 * t496 + t499 * t414;
t519 = -t495 * t415 + t498 * t421;
t517 = qJD(1) * ((t545 * qJD(3)) + (t531 * t495 - t533 * t498) * qJD(1));
t437 = Ifges(6,1) * t467 + Ifges(6,4) * t466 + Ifges(6,5) * t483;
t436 = Ifges(6,4) * t467 + Ifges(6,2) * t466 + Ifges(6,6) * t483;
t435 = Ifges(6,5) * t467 + Ifges(6,6) * t466 + Ifges(6,3) * t483;
t419 = mrSges(6,2) * t429 - mrSges(6,3) * t425 + Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * t465 + t435 * t466 - t436 * t483;
t418 = -mrSges(6,1) * t429 + mrSges(6,3) * t426 + Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * t465 - t435 * t467 + t437 * t483;
t416 = -mrSges(5,2) * t471 + t503;
t410 = -m(3) * g(3) + t519;
t407 = mrSges(5,1) * t434 + mrSges(6,1) * t425 + mrSges(4,2) * t449 - mrSges(6,2) * t426 - mrSges(4,3) * t443 - mrSges(5,3) * t432 + Ifges(6,5) * t441 + Ifges(6,6) * t440 + Ifges(6,3) * t465 + pkin(4) * t417 - qJ(4) * t416 + t467 * t436 - t466 * t437 + t547 * t472 - t535 * t471 + t533 * qJDD(3) + t528 * qJD(3) + t495 * t517;
t406 = -mrSges(4,1) * t449 - mrSges(5,1) * t433 + mrSges(5,2) * t432 + mrSges(4,3) * t444 - pkin(3) * t416 - pkin(4) * t509 - pkin(7) * t518 + t527 * qJD(3) + t531 * qJDD(3) - t497 * t418 - t494 * t419 - t546 * t471 + t535 * t472 + t498 * t517;
t405 = -t494 * t418 + t497 * t419 - mrSges(2,3) * t479 + mrSges(4,1) * t443 - mrSges(4,2) * t444 + mrSges(3,1) * t452 - mrSges(5,3) * t433 + mrSges(5,2) * t434 - pkin(7) * t417 - qJ(2) * t410 + pkin(2) * t411 + qJ(4) * t504 + pkin(3) * (-qJD(3) * t476 + t507) + (t532 * t501) + t533 * t472 + (-mrSges(5,1) * qJ(4) - t531) * t471 + (-mrSges(5,2) * pkin(3) - t545) * qJDD(3) + t534 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + ((-pkin(3) * t469 - t528) * t498 + (-qJ(4) * t469 + t527) * t495) * qJD(1);
t404 = -mrSges(3,1) * t451 + mrSges(2,3) * t480 - pkin(1) * t410 - pkin(2) * t544 - pkin(6) * t519 + t537 * g(3) - t532 * qJDD(1) - t498 * t406 - t495 * t407 + t534 * t501;
t1 = [-m(1) * g(1) + t520; -m(1) * g(2) + t529; (-m(1) - m(2) - m(3)) * g(3) + t519; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t529 - t496 * t404 + t499 * t405; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t520 + t499 * t404 + t496 * t405; pkin(1) * t508 + qJ(2) * t502 - pkin(6) * t411 + mrSges(2,1) * t479 - mrSges(2,2) * t480 + t498 * t407 - t495 * t406 + mrSges(3,2) * t452 - mrSges(3,3) * t451 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
