% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:39
% EndTime: 2019-12-31 16:41:40
% DurationCPUTime: 0.77s
% Computational Cost: add. (5260->166), mult. (11081->203), div. (0->0), fcn. (6151->6), ass. (0->78)
t468 = sin(qJ(1));
t470 = cos(qJ(1));
t447 = t468 * g(1) - t470 * g(2);
t471 = qJD(1) ^ 2;
t478 = -t471 * qJ(2) + qJDD(2) - t447;
t497 = -pkin(1) - qJ(3);
t505 = -(2 * qJD(1) * qJD(3)) + t497 * qJDD(1) + t478;
t448 = -t470 * g(1) - t468 * g(2);
t504 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t448;
t439 = t471 * pkin(1) - t504;
t503 = -m(3) * t439 + t471 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t501 = pkin(3) * t471;
t500 = mrSges(2,1) - mrSges(3,2);
t499 = -Ifges(3,4) + Ifges(2,5);
t498 = -Ifges(2,6) + Ifges(3,5);
t466 = cos(pkin(6));
t496 = mrSges(4,2) * t466;
t465 = sin(pkin(6));
t428 = t465 * g(3) + t505 * t466;
t416 = (-pkin(5) * qJDD(1) - t465 * t501) * t466 + t428;
t429 = -t466 * g(3) + t505 * t465;
t456 = t465 ^ 2;
t491 = qJDD(1) * t465;
t417 = -pkin(5) * t491 - t456 * t501 + t429;
t467 = sin(qJ(4));
t469 = cos(qJ(4));
t413 = t469 * t416 - t467 * t417;
t482 = -t465 * t469 - t466 * t467;
t442 = t482 * qJD(1);
t481 = -t465 * t467 + t466 * t469;
t443 = t481 * qJD(1);
t424 = -t442 * mrSges(5,1) + t443 * mrSges(5,2);
t431 = t442 * qJD(4) + t481 * qJDD(1);
t436 = -qJD(4) * mrSges(5,2) + t442 * mrSges(5,3);
t410 = m(5) * t413 + qJDD(4) * mrSges(5,1) - t431 * mrSges(5,3) + qJD(4) * t436 - t443 * t424;
t414 = t467 * t416 + t469 * t417;
t430 = -t443 * qJD(4) + t482 * qJDD(1);
t437 = qJD(4) * mrSges(5,1) - t443 * mrSges(5,3);
t411 = m(5) * t414 - qJDD(4) * mrSges(5,2) + t430 * mrSges(5,3) - qJD(4) * t437 + t442 * t424;
t399 = t469 * t410 + t467 * t411;
t480 = -qJDD(1) * mrSges(4,3) - t471 * (mrSges(4,1) * t465 + t496);
t397 = m(4) * t428 + t480 * t466 + t399;
t486 = -t467 * t410 + t469 * t411;
t398 = m(4) * t429 + t480 * t465 + t486;
t395 = t466 * t397 + t465 * t398;
t441 = -qJDD(1) * pkin(1) + t478;
t475 = -m(3) * t441 + t471 * mrSges(3,3) - t395;
t391 = m(2) * t447 - t471 * mrSges(2,2) + t500 * qJDD(1) + t475;
t477 = qJDD(3) + t504;
t435 = t497 * t471 + t477;
t494 = -t466 ^ 2 - t456;
t419 = pkin(3) * t491 + (t494 * pkin(5) + t497) * t471 + t477;
t476 = m(5) * t419 - t430 * mrSges(5,1) + t431 * mrSges(5,2) - t442 * t436 + t443 * t437;
t473 = m(4) * t435 + mrSges(4,1) * t491 + qJDD(1) * t496 + t476;
t489 = t494 * mrSges(4,3);
t404 = m(2) * t448 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) + t489) * t471 + t473 + t503;
t495 = t470 * t391 + t468 * t404;
t483 = Ifges(4,5) * t466 - Ifges(4,6) * t465;
t493 = t471 * t483;
t488 = -t468 * t391 + t470 * t404;
t487 = -t465 * t397 + t466 * t398;
t485 = Ifges(4,1) * t466 - Ifges(4,4) * t465;
t484 = Ifges(4,4) * t466 - Ifges(4,2) * t465;
t421 = Ifges(5,4) * t443 + Ifges(5,2) * t442 + Ifges(5,6) * qJD(4);
t422 = Ifges(5,1) * t443 + Ifges(5,4) * t442 + Ifges(5,5) * qJD(4);
t474 = mrSges(5,1) * t413 - mrSges(5,2) * t414 + Ifges(5,5) * t431 + Ifges(5,6) * t430 + Ifges(5,3) * qJDD(4) + t443 * t421 - t442 * t422;
t406 = t471 * t489 + t473;
t420 = Ifges(5,5) * t443 + Ifges(5,6) * t442 + Ifges(5,3) * qJD(4);
t400 = -mrSges(5,1) * t419 + mrSges(5,3) * t414 + Ifges(5,4) * t431 + Ifges(5,2) * t430 + Ifges(5,6) * qJDD(4) + qJD(4) * t422 - t443 * t420;
t401 = mrSges(5,2) * t419 - mrSges(5,3) * t413 + Ifges(5,1) * t431 + Ifges(5,4) * t430 + Ifges(5,5) * qJDD(4) - qJD(4) * t421 + t442 * t420;
t387 = -mrSges(4,1) * t435 + mrSges(4,3) * t429 - pkin(3) * t476 + pkin(5) * t486 + t484 * qJDD(1) + t469 * t400 + t467 * t401 - t466 * t493;
t389 = mrSges(4,2) * t435 - mrSges(4,3) * t428 - pkin(5) * t399 + t485 * qJDD(1) - t467 * t400 + t469 * t401 - t465 * t493;
t393 = qJDD(1) * mrSges(3,2) - t475;
t472 = -mrSges(2,2) * t448 - mrSges(3,3) * t439 - pkin(1) * t393 - qJ(3) * t395 - t465 * t387 + t466 * t389 + qJ(2) * (t406 + t503) + mrSges(3,2) * t441 + mrSges(2,1) * t447 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t394 = -m(3) * g(3) + t487;
t386 = -qJ(2) * t394 + mrSges(3,1) * t441 - mrSges(2,3) * t447 + t474 + pkin(2) * t395 + pkin(3) * t399 + mrSges(4,1) * t428 - mrSges(4,2) * t429 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t483 + t499) * qJDD(1) + (t465 * t485 + t466 * t484 + t498) * t471;
t385 = -mrSges(3,1) * t439 + mrSges(2,3) * t448 - pkin(1) * t394 + pkin(2) * t406 + t500 * g(3) - qJ(3) * t487 - t498 * qJDD(1) - t466 * t387 - t465 * t389 + t499 * t471;
t1 = [-m(1) * g(1) + t488; -m(1) * g(2) + t495; (-m(1) - m(2) - m(3)) * g(3) + t487; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t495 - t468 * t385 + t470 * t386; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t488 + t470 * t385 + t468 * t386; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t472; t472; t393; t406; t474;];
tauJB = t1;
