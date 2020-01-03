% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPP4
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:57
% EndTime: 2019-12-31 16:58:59
% DurationCPUTime: 0.82s
% Computational Cost: add. (2830->190), mult. (5787->232), div. (0->0), fcn. (2323->4), ass. (0->72)
t462 = Ifges(3,1) + Ifges(4,1) + Ifges(5,1);
t451 = Ifges(3,4) - Ifges(4,5) - Ifges(5,4);
t450 = Ifges(3,5) + Ifges(4,4) - Ifges(5,5);
t461 = Ifges(3,2) + Ifges(4,3) + Ifges(5,2);
t449 = Ifges(3,6) - Ifges(4,6) + Ifges(5,6);
t460 = -Ifges(3,3) - Ifges(4,2) - Ifges(5,3);
t459 = 2 * qJD(3);
t434 = qJD(1) ^ 2;
t458 = t434 * pkin(5);
t457 = mrSges(3,3) + mrSges(4,2);
t431 = cos(qJ(2));
t456 = qJ(3) * t431;
t430 = sin(qJ(1));
t432 = cos(qJ(1));
t417 = -t432 * g(1) - t430 * g(2);
t387 = -t434 * pkin(1) + qJDD(1) * pkin(5) + t417;
t429 = sin(qJ(2));
t370 = -t429 * g(3) + t431 * t387;
t401 = (mrSges(5,1) * t431 + mrSges(5,2) * t429) * qJD(1);
t402 = (-mrSges(3,1) * t431 + mrSges(3,2) * t429) * qJD(1);
t452 = qJD(1) * qJD(2);
t404 = t431 * qJDD(1) - t429 * t452;
t454 = qJD(1) * t429;
t411 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t454;
t399 = (-pkin(2) * t431 - qJ(3) * t429) * qJD(1);
t433 = qJD(2) ^ 2;
t453 = qJD(1) * t431;
t367 = -t433 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t459 + t399 * t453 + t370;
t400 = (-mrSges(4,1) * t431 - mrSges(4,3) * t429) * qJD(1);
t412 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t454;
t409 = -qJD(2) * pkin(3) - qJ(4) * t454;
t428 = t431 ^ 2;
t448 = -0.2e1 * qJD(1) * qJD(4);
t363 = -t428 * t434 * pkin(3) - t404 * qJ(4) + qJD(2) * t409 + t431 * t448 + t367;
t410 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t454;
t441 = m(5) * t363 + qJDD(2) * mrSges(5,2) - t404 * mrSges(5,3) + qJD(2) * t410;
t437 = m(4) * t367 + qJDD(2) * mrSges(4,3) + qJD(2) * t412 + t400 * t453 + t441;
t356 = m(3) * t370 - qJDD(2) * mrSges(3,2) - qJD(2) * t411 + t457 * t404 + (-t401 + t402) * t453 + t437;
t369 = -t431 * g(3) - t429 * t387;
t444 = t431 * t452;
t403 = t429 * qJDD(1) + t444;
t414 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t453;
t368 = -qJDD(2) * pkin(2) - t433 * qJ(3) + t399 * t454 + qJDD(3) - t369;
t364 = t429 * t448 + (-t403 + t444) * qJ(4) + (-t429 * t431 * t434 - qJDD(2)) * pkin(3) + t368;
t413 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t453;
t359 = m(5) * t364 - qJDD(2) * mrSges(5,1) - t403 * mrSges(5,3) - qJD(2) * t413 - t401 * t454;
t415 = mrSges(4,2) * t453 + qJD(2) * mrSges(4,3);
t436 = -m(4) * t368 + qJDD(2) * mrSges(4,1) + qJD(2) * t415 - t359;
t357 = m(3) * t369 + qJDD(2) * mrSges(3,1) + qJD(2) * t414 - t457 * t403 + (-t400 - t402) * t454 + t436;
t442 = t431 * t356 - t429 * t357;
t349 = m(2) * t417 - t434 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t442;
t416 = t430 * g(1) - t432 * g(2);
t440 = qJDD(1) * pkin(1) + t416;
t438 = -t403 * qJ(3) - t440;
t365 = -t404 * pkin(2) - t458 + (-0.2e1 * qJD(3) * t429 + (pkin(2) * t429 - t456) * qJD(2)) * qJD(1) + t438;
t361 = qJDD(4) + (-qJ(4) * t428 + pkin(5)) * t434 + (pkin(2) + pkin(3)) * t404 + (qJD(2) * t456 + (-pkin(2) * qJD(2) + t409 + t459) * t429) * qJD(1) - t438;
t439 = -m(5) * t361 - t404 * mrSges(5,1) - t403 * mrSges(5,2) - t410 * t454 - t413 * t453;
t358 = m(4) * t365 - t404 * mrSges(4,1) - t403 * mrSges(4,3) - t412 * t454 - t415 * t453 + t439;
t386 = -t440 - t458;
t435 = -m(3) * t386 + t404 * mrSges(3,1) - t403 * mrSges(3,2) - t411 * t454 + t414 * t453 - t358;
t352 = m(2) * t416 + qJDD(1) * mrSges(2,1) - t434 * mrSges(2,2) + t435;
t455 = t430 * t349 + t432 * t352;
t350 = t429 * t356 + t431 * t357;
t447 = t460 * qJD(2) + (-t450 * t429 - t449 * t431) * qJD(1);
t446 = -t449 * qJD(2) + (-t451 * t429 - t461 * t431) * qJD(1);
t445 = t450 * qJD(2) + (t462 * t429 + t451 * t431) * qJD(1);
t443 = t432 * t349 - t430 * t352;
t346 = mrSges(3,2) * t386 + mrSges(4,2) * t368 + mrSges(5,2) * t361 - mrSges(3,3) * t369 - mrSges(4,3) * t365 - mrSges(5,3) * t364 - qJ(3) * t358 - qJ(4) * t359 + t446 * qJD(2) + t450 * qJDD(2) + t462 * t403 + t451 * t404 - t447 * t453;
t345 = -mrSges(3,1) * t386 + mrSges(3,3) * t370 - mrSges(4,1) * t365 + mrSges(4,2) * t367 + mrSges(5,1) * t361 - mrSges(5,3) * t363 - pkin(3) * t439 - qJ(4) * t441 - pkin(2) * t358 + t461 * t404 + t451 * t403 + t449 * qJDD(2) + t445 * qJD(2) + (qJ(4) * t401 * t431 + t447 * t429) * qJD(1);
t344 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) - pkin(1) * t350 + pkin(3) * t359 - mrSges(5,2) * t363 + mrSges(5,1) * t364 - mrSges(4,3) * t367 + mrSges(4,1) * t368 - mrSges(3,1) * t369 + mrSges(3,2) * t370 + mrSges(2,3) * t417 - pkin(2) * t436 - qJ(3) * t437 + t434 * Ifges(2,5) + (-qJ(3) * mrSges(4,2) - t449) * t404 + (pkin(2) * mrSges(4,2) - t450) * t403 + t460 * qJDD(2) + ((qJ(3) * t401 + t445) * t431 + (pkin(2) * t400 + t446) * t429) * qJD(1);
t343 = -mrSges(2,2) * g(3) - mrSges(2,3) * t416 + Ifges(2,5) * qJDD(1) - t434 * Ifges(2,6) - pkin(5) * t350 - t429 * t345 + t431 * t346;
t1 = [-m(1) * g(1) + t443; -m(1) * g(2) + t455; (-m(1) - m(2)) * g(3) + t350; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t455 + t432 * t343 - t430 * t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t443 + t430 * t343 + t432 * t344; -mrSges(1,1) * g(2) + mrSges(2,1) * t416 + mrSges(1,2) * g(1) - mrSges(2,2) * t417 + Ifges(2,3) * qJDD(1) + pkin(1) * t435 + pkin(5) * t442 + t431 * t345 + t429 * t346;];
tauB = t1;
