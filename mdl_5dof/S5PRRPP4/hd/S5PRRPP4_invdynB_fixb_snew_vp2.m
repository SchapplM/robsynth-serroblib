% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRPP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRPP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:52
% EndTime: 2019-12-31 17:40:54
% DurationCPUTime: 1.12s
% Computational Cost: add. (7103->212), mult. (13958->259), div. (0->0), fcn. (6511->6), ass. (0->88)
t537 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t523 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t522 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t536 = Ifges(4,2) + Ifges(5,3) + Ifges(6,2);
t521 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t535 = -Ifges(4,3) - Ifges(5,2) - Ifges(6,3);
t534 = 2 * qJD(4);
t533 = pkin(3) + pkin(4);
t504 = qJD(2) ^ 2;
t532 = pkin(4) * t504;
t531 = t504 * pkin(6);
t530 = mrSges(4,3) + mrSges(5,2);
t501 = cos(qJ(3));
t529 = qJ(4) * t501;
t496 = -g(3) + qJDD(1);
t528 = t501 * t496;
t497 = sin(pkin(7));
t498 = cos(pkin(7));
t477 = t497 * g(1) - t498 * g(2);
t478 = -t498 * g(1) - t497 * g(2);
t500 = sin(qJ(2));
t502 = cos(qJ(2));
t438 = t500 * t477 + t502 * t478;
t436 = -t504 * pkin(2) + qJDD(2) * pkin(6) + t438;
t499 = sin(qJ(3));
t432 = t501 * t436 + t499 * t496;
t469 = (mrSges(6,1) * t501 + mrSges(6,2) * t499) * qJD(2);
t470 = (-mrSges(4,1) * t501 + mrSges(4,2) * t499) * qJD(2);
t524 = qJD(2) * qJD(3);
t472 = t501 * qJDD(2) - t499 * t524;
t526 = qJD(2) * t499;
t481 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t526;
t467 = (-pkin(3) * t501 - qJ(4) * t499) * qJD(2);
t503 = qJD(3) ^ 2;
t525 = qJD(2) * t501;
t429 = -t503 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t534 + t467 * t525 + t432;
t468 = (-mrSges(5,1) * t501 - mrSges(5,3) * t499) * qJD(2);
t482 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t526;
t479 = -qJD(3) * pkin(4) - qJ(5) * t526;
t495 = t501 ^ 2;
t520 = -0.2e1 * qJD(2) * qJD(5);
t425 = -t472 * qJ(5) + qJD(3) * t479 - t495 * t532 + t501 * t520 + t429;
t480 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t526;
t512 = m(6) * t425 + qJDD(3) * mrSges(6,2) - t472 * mrSges(6,3) + qJD(3) * t480;
t507 = m(5) * t429 + qJDD(3) * mrSges(5,3) + qJD(3) * t482 + t468 * t525 + t512;
t418 = m(4) * t432 - qJDD(3) * mrSges(4,2) - qJD(3) * t481 + t530 * t472 + (-t469 + t470) * t525 + t507;
t433 = t499 * t436;
t431 = -t433 + t528;
t471 = t499 * qJDD(2) + t501 * t524;
t484 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t525;
t509 = -t503 * qJ(4) + t467 * t526 + qJDD(4) + t433;
t426 = t499 * t520 - t471 * qJ(5) - t533 * qJDD(3) + (qJ(5) * t524 - t499 * t532 - t496) * t501 + t509;
t483 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t525;
t421 = m(6) * t426 - qJDD(3) * mrSges(6,1) - t471 * mrSges(6,3) - qJD(3) * t483 - t469 * t526;
t430 = -qJDD(3) * pkin(3) + t509 - t528;
t485 = mrSges(5,2) * t525 + qJD(3) * mrSges(5,3);
t506 = -m(5) * t430 + qJDD(3) * mrSges(5,1) + qJD(3) * t485 - t421;
t419 = m(4) * t431 + qJDD(3) * mrSges(4,1) + qJD(3) * t484 - t530 * t471 + (-t468 - t470) * t526 + t506;
t513 = t501 * t418 - t499 * t419;
t411 = m(3) * t438 - t504 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t513;
t437 = t502 * t477 - t500 * t478;
t511 = qJDD(2) * pkin(2) + t437;
t508 = -t471 * qJ(4) - t511;
t427 = -t472 * pkin(3) - t531 + (-0.2e1 * qJD(4) * t499 + (pkin(3) * t499 - t529) * qJD(3)) * qJD(2) + t508;
t423 = qJDD(5) + (-qJ(5) * t495 + pkin(6)) * t504 + t533 * t472 + (qJD(3) * t529 + (-pkin(3) * qJD(3) + t479 + t534) * t499) * qJD(2) - t508;
t510 = -m(6) * t423 - t472 * mrSges(6,1) - t471 * mrSges(6,2) - t480 * t526 - t483 * t525;
t420 = m(5) * t427 - t472 * mrSges(5,1) - t471 * mrSges(5,3) - t482 * t526 - t485 * t525 + t510;
t435 = -t511 - t531;
t505 = -m(4) * t435 + t472 * mrSges(4,1) - t471 * mrSges(4,2) - t481 * t526 + t484 * t525 - t420;
t414 = m(3) * t437 + qJDD(2) * mrSges(3,1) - t504 * mrSges(3,2) + t505;
t407 = t500 * t411 + t502 * t414;
t405 = m(2) * t477 + t407;
t514 = t502 * t411 - t500 * t414;
t406 = m(2) * t478 + t514;
t527 = t498 * t405 + t497 * t406;
t412 = t499 * t418 + t501 * t419;
t519 = t535 * qJD(3) + (-t522 * t499 - t521 * t501) * qJD(2);
t518 = -t521 * qJD(3) + (-t523 * t499 - t536 * t501) * qJD(2);
t517 = t522 * qJD(3) + (t537 * t499 + t523 * t501) * qJD(2);
t516 = m(3) * t496 + t412;
t515 = -t497 * t405 + t498 * t406;
t408 = mrSges(4,2) * t435 + mrSges(5,2) * t430 + mrSges(6,2) * t423 - mrSges(4,3) * t431 - mrSges(5,3) * t427 - mrSges(6,3) * t426 - qJ(4) * t420 - qJ(5) * t421 + t518 * qJD(3) + t522 * qJDD(3) + t537 * t471 + t523 * t472 - t519 * t525;
t401 = -mrSges(4,1) * t435 + mrSges(4,3) * t432 - mrSges(5,1) * t427 + mrSges(5,2) * t429 + mrSges(6,1) * t423 - mrSges(6,3) * t425 - pkin(4) * t510 - qJ(5) * t512 - pkin(3) * t420 + t536 * t472 + t523 * t471 + t521 * qJDD(3) + t517 * qJD(3) + (qJ(5) * t469 * t501 + t519 * t499) * qJD(2);
t400 = mrSges(3,3) * t438 - mrSges(6,2) * t425 + mrSges(6,1) * t426 - mrSges(5,3) * t429 + mrSges(5,1) * t430 - mrSges(4,1) * t431 + mrSges(4,2) * t432 + pkin(4) * t421 - pkin(3) * t506 - qJ(4) * t507 + t504 * Ifges(3,5) - mrSges(3,1) * t496 - pkin(2) * t412 + Ifges(3,6) * qJDD(2) + (-qJ(4) * mrSges(5,2) - t521) * t472 + (pkin(3) * mrSges(5,2) - t522) * t471 + t535 * qJDD(3) + ((qJ(4) * t469 + t517) * t501 + (pkin(3) * t468 + t518) * t499) * qJD(2);
t399 = mrSges(3,2) * t496 - mrSges(3,3) * t437 + Ifges(3,5) * qJDD(2) - t504 * Ifges(3,6) - pkin(6) * t412 - t499 * t401 + t501 * t408;
t398 = mrSges(2,2) * t496 - mrSges(2,3) * t477 - pkin(5) * t407 + t502 * t399 - t500 * t400;
t397 = -mrSges(2,1) * t496 + mrSges(2,3) * t478 - pkin(1) * t516 + pkin(5) * t514 + t500 * t399 + t502 * t400;
t1 = [-m(1) * g(1) + t515; -m(1) * g(2) + t527; -m(1) * g(3) + m(2) * t496 + t516; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t527 - t497 * t397 + t498 * t398; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t515 + t498 * t397 + t497 * t398; -mrSges(1,1) * g(2) + mrSges(2,1) * t477 + mrSges(3,1) * t437 + mrSges(1,2) * g(1) - mrSges(2,2) * t478 - mrSges(3,2) * t438 + Ifges(3,3) * qJDD(2) + pkin(1) * t407 + pkin(2) * t505 + pkin(6) * t513 + t501 * t401 + t499 * t408;];
tauB = t1;
