% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRP6
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:03
% EndTime: 2019-12-31 16:46:04
% DurationCPUTime: 0.73s
% Computational Cost: add. (2694->170), mult. (5067->195), div. (0->0), fcn. (1960->4), ass. (0->71)
t507 = Ifges(4,1) + Ifges(5,1);
t496 = Ifges(4,4) + Ifges(5,4);
t505 = Ifges(4,5) + Ifges(5,5);
t506 = Ifges(4,2) + Ifges(5,2);
t504 = Ifges(4,6) + Ifges(5,6);
t503 = Ifges(4,3) + Ifges(5,3);
t467 = sin(qJ(3));
t469 = cos(qJ(3));
t502 = t504 * qJD(3) + (-t506 * t467 + t496 * t469) * qJD(1);
t501 = t505 * qJD(3) + (-t496 * t467 + t507 * t469) * qJD(1);
t468 = sin(qJ(1));
t470 = cos(qJ(1));
t453 = -t470 * g(1) - t468 * g(2);
t478 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t453;
t471 = qJD(1) ^ 2;
t499 = -pkin(1) - pkin(5);
t415 = t499 * t471 + t478;
t488 = qJD(1) * qJD(3);
t443 = -t467 * qJDD(1) - t469 * t488;
t484 = t467 * t488;
t444 = t469 * qJDD(1) - t484;
t490 = qJD(1) * t467;
t448 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t490;
t489 = qJD(1) * t469;
t451 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t489;
t449 = qJD(3) * pkin(3) - qJ(4) * t489;
t464 = t467 ^ 2;
t408 = t449 * t489 - t443 * pkin(3) + qJDD(4) + (-qJ(4) * t464 + t499) * t471 + t478;
t447 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t490;
t450 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t489;
t480 = m(5) * t408 + t444 * mrSges(5,2) + t447 * t490 + t450 * t489;
t500 = -m(4) * t415 - t444 * mrSges(4,2) + (mrSges(4,1) + mrSges(5,1)) * t443 - t448 * t490 - t451 * t489 - t480;
t498 = mrSges(2,1) - mrSges(3,2);
t495 = Ifges(2,5) - Ifges(3,4);
t494 = -Ifges(2,6) + Ifges(3,5);
t452 = t468 * g(1) - t470 * g(2);
t477 = -t471 * qJ(2) + qJDD(2) - t452;
t416 = t499 * qJDD(1) + t477;
t410 = t467 * g(3) + t469 * t416;
t441 = (mrSges(5,1) * t467 + mrSges(5,2) * t469) * qJD(1);
t481 = qJD(1) * (-t441 - (mrSges(4,1) * t467 + mrSges(4,2) * t469) * qJD(1));
t486 = -0.2e1 * qJD(1) * qJD(4);
t405 = t469 * t486 + (-t444 - t484) * qJ(4) + (-t467 * t469 * t471 + qJDD(3)) * pkin(3) + t410;
t485 = m(5) * t405 + qJDD(3) * mrSges(5,1) + qJD(3) * t447;
t397 = m(4) * t410 + qJDD(3) * mrSges(4,1) + qJD(3) * t448 + (-mrSges(4,3) - mrSges(5,3)) * t444 + t469 * t481 + t485;
t411 = -t469 * g(3) + t467 * t416;
t406 = -t464 * t471 * pkin(3) + t443 * qJ(4) - qJD(3) * t449 + t467 * t486 + t411;
t492 = m(5) * t406 + t443 * mrSges(5,3);
t398 = m(4) * t411 + t443 * mrSges(4,3) + (-mrSges(4,2) - mrSges(5,2)) * qJDD(3) + (-t450 - t451) * qJD(3) + t467 * t481 + t492;
t391 = t469 * t397 + t467 * t398;
t423 = -qJDD(1) * pkin(1) + t477;
t475 = -m(3) * t423 + t471 * mrSges(3,3) - t391;
t387 = m(2) * t452 - t471 * mrSges(2,2) + t498 * qJDD(1) + t475;
t421 = t471 * pkin(1) - t478;
t474 = -m(3) * t421 + t471 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t500;
t394 = m(2) * t453 - t471 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t474;
t493 = t470 * t387 + t468 * t394;
t491 = -t503 * qJD(3) + (t504 * t467 - t505 * t469) * qJD(1);
t483 = -t468 * t387 + t470 * t394;
t482 = -t467 * t397 + t469 * t398;
t400 = -t443 * mrSges(5,1) + t480;
t383 = -mrSges(4,1) * t415 + mrSges(4,3) * t411 - mrSges(5,1) * t408 + mrSges(5,3) * t406 - pkin(3) * t400 + qJ(4) * t492 + t496 * t444 + t506 * t443 + (-qJ(4) * mrSges(5,2) + t504) * qJDD(3) + (-qJ(4) * t450 + t501) * qJD(3) + (-qJ(4) * t441 * t467 + t491 * t469) * qJD(1);
t401 = -t444 * mrSges(5,3) - t441 * t489 + t485;
t385 = mrSges(4,2) * t415 + mrSges(5,2) * t408 - mrSges(4,3) * t410 - mrSges(5,3) * t405 - qJ(4) * t401 - t502 * qJD(3) + t505 * qJDD(3) + t496 * t443 + t507 * t444 + t491 * t490;
t389 = qJDD(1) * mrSges(3,2) - t475;
t473 = mrSges(2,1) * t452 - mrSges(2,2) * t453 + mrSges(3,2) * t423 - mrSges(3,3) * t421 - pkin(1) * t389 - pkin(5) * t391 + qJ(2) * t474 - t467 * t383 + t469 * t385 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t472 = mrSges(4,1) * t410 + mrSges(5,1) * t405 - mrSges(4,2) * t411 - mrSges(5,2) * t406 + pkin(3) * t401 + t503 * qJDD(3) + t504 * t443 + t505 * t444 + t502 * t489 + t501 * t490;
t390 = -m(3) * g(3) + t482;
t382 = t472 + t494 * t471 + mrSges(3,1) * t423 - mrSges(2,3) * t452 - qJ(2) * t390 + pkin(2) * t391 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t495 * qJDD(1);
t381 = -mrSges(3,1) * t421 + mrSges(2,3) * t453 - pkin(1) * t390 - pkin(2) * t500 - pkin(5) * t482 + t498 * g(3) - t494 * qJDD(1) - t469 * t383 - t467 * t385 + t495 * t471;
t1 = [-m(1) * g(1) + t483; -m(1) * g(2) + t493; (-m(1) - m(2) - m(3)) * g(3) + t482; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t493 - t468 * t381 + t470 * t382; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t483 + t470 * t381 + t468 * t382; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t473; t473; t389; t472; t400;];
tauJB = t1;
