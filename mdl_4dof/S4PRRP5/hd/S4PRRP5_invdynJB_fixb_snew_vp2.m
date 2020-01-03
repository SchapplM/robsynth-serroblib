% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRP5
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:56
% EndTime: 2019-12-31 16:28:57
% DurationCPUTime: 0.69s
% Computational Cost: add. (4016->165), mult. (7636->201), div. (0->0), fcn. (3684->6), ass. (0->75)
t472 = Ifges(4,1) + Ifges(5,1);
t464 = Ifges(4,4) + Ifges(5,4);
t463 = Ifges(4,5) + Ifges(5,5);
t471 = Ifges(4,2) + Ifges(5,2);
t462 = Ifges(4,6) + Ifges(5,6);
t470 = Ifges(4,3) + Ifges(5,3);
t436 = sin(qJ(3));
t438 = cos(qJ(3));
t417 = (-mrSges(5,1) * t438 + mrSges(5,2) * t436) * qJD(2);
t452 = qJD(2) * qJD(3);
t447 = t438 * t452;
t419 = t436 * qJDD(2) + t447;
t435 = sin(pkin(6));
t460 = cos(pkin(6));
t423 = -t460 * g(1) - t435 * g(2);
t434 = -g(3) + qJDD(1);
t437 = sin(qJ(2));
t439 = cos(qJ(2));
t400 = t439 * t423 + t437 * t434;
t440 = qJD(2) ^ 2;
t398 = -t440 * pkin(2) + qJDD(2) * pkin(5) + t400;
t422 = t435 * g(1) - t460 * g(2);
t415 = t438 * t422;
t451 = qJD(2) * qJD(4);
t466 = pkin(3) * t440;
t391 = qJDD(3) * pkin(3) - t415 + (-t419 + t447) * qJ(4) + (t438 * t466 - t398 - 0.2e1 * t451) * t436;
t453 = qJD(2) * t438;
t427 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t453;
t450 = m(5) * t391 + qJDD(3) * mrSges(5,1) + qJD(3) * t427;
t454 = qJD(2) * t436;
t387 = -t419 * mrSges(5,3) - t417 * t454 + t450;
t395 = t438 * t398 - t436 * t422;
t420 = t438 * qJDD(2) - t436 * t452;
t424 = qJD(3) * pkin(3) - qJ(4) * t454;
t433 = t438 ^ 2;
t392 = t420 * qJ(4) - qJD(3) * t424 - t433 * t466 + 0.2e1 * t438 * t451 + t395;
t394 = -t436 * t398 - t415;
t456 = t463 * qJD(3) + (t472 * t436 + t464 * t438) * qJD(2);
t457 = t462 * qJD(3) + (t464 * t436 + t471 * t438) * qJD(2);
t469 = mrSges(4,1) * t394 + mrSges(5,1) * t391 - mrSges(4,2) * t395 - mrSges(5,2) * t392 + pkin(3) * t387 + (t457 * t436 - t456 * t438) * qJD(2) + t470 * qJDD(3) + t463 * t419 + t462 * t420;
t465 = -mrSges(4,2) - mrSges(5,2);
t418 = (-mrSges(4,1) * t438 + mrSges(4,2) * t436) * qJD(2);
t428 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t453;
t384 = m(4) * t394 + qJDD(3) * mrSges(4,1) + qJD(3) * t428 + (-mrSges(4,3) - mrSges(5,3)) * t419 + (-t417 - t418) * t454 + t450;
t449 = m(5) * t392 + t420 * mrSges(5,3) + t417 * t453;
t425 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t454;
t455 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t454 - t425;
t385 = m(4) * t395 + t420 * mrSges(4,3) + t455 * qJD(3) + t465 * qJDD(3) + t418 * t453 + t449;
t380 = -t436 * t384 + t438 * t385;
t376 = m(3) * t400 - t440 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t380;
t399 = -t437 * t423 + t439 * t434;
t443 = -qJDD(2) * pkin(2) - t399;
t397 = -t440 * pkin(5) + t443;
t393 = t424 * t454 - t420 * pkin(3) + qJDD(4) + (-qJ(4) * t433 - pkin(5)) * t440 + t443;
t444 = -m(5) * t393 + t420 * mrSges(5,1) + t427 * t453;
t386 = -m(4) * t397 + t420 * mrSges(4,1) + t465 * t419 + t428 * t453 + t455 * t454 + t444;
t382 = m(3) * t399 + qJDD(2) * mrSges(3,1) - t440 * mrSges(3,2) + t386;
t445 = t439 * t376 - t437 * t382;
t370 = m(2) * t423 + t445;
t379 = t438 * t384 + t436 * t385;
t378 = (m(2) + m(3)) * t422 - t379;
t459 = t435 * t370 + t460 * t378;
t371 = t437 * t376 + t439 * t382;
t458 = t470 * qJD(3) + (t463 * t436 + t462 * t438) * qJD(2);
t448 = m(2) * t434 + t371;
t446 = t460 * t370 - t435 * t378;
t388 = t419 * mrSges(5,2) + t425 * t454 - t444;
t372 = -mrSges(4,1) * t397 + mrSges(4,3) * t395 - mrSges(5,1) * t393 + mrSges(5,3) * t392 - pkin(3) * t388 + qJ(4) * t449 + t471 * t420 + t464 * t419 + (-qJ(4) * mrSges(5,2) + t462) * qJDD(3) + (-qJ(4) * t425 + t456) * qJD(3) - t458 * t454;
t373 = mrSges(4,2) * t397 + mrSges(5,2) * t393 - mrSges(4,3) * t394 - mrSges(5,3) * t391 - qJ(4) * t387 - t457 * qJD(3) + t463 * qJDD(3) + t472 * t419 + t464 * t420 + t458 * t453;
t441 = mrSges(3,1) * t399 - mrSges(3,2) * t400 + Ifges(3,3) * qJDD(2) + pkin(2) * t386 + pkin(5) * t380 + t438 * t372 + t436 * t373;
t367 = mrSges(3,1) * t422 + mrSges(3,3) * t400 + t440 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t379 - t469;
t366 = -mrSges(3,2) * t422 - mrSges(3,3) * t399 + Ifges(3,5) * qJDD(2) - t440 * Ifges(3,6) - pkin(5) * t379 - t436 * t372 + t438 * t373;
t365 = -mrSges(2,1) * t434 + mrSges(2,3) * t423 - pkin(1) * t371 - t441;
t364 = mrSges(2,2) * t434 - mrSges(2,3) * t422 - pkin(4) * t371 + t439 * t366 - t437 * t367;
t1 = [-m(1) * g(1) + t446; -m(1) * g(2) + t459; -m(1) * g(3) + t448; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t459 + t460 * t364 - t435 * t365; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t446 + t435 * t364 + t460 * t365; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t422 - mrSges(2,2) * t423 + t437 * t366 + t439 * t367 + pkin(1) * (m(3) * t422 - t379) + pkin(4) * t445; t448; t441; t469; t388;];
tauJB = t1;
