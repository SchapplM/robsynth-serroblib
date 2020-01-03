% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRP4
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:17
% EndTime: 2019-12-31 17:15:19
% DurationCPUTime: 1.26s
% Computational Cost: add. (9331->218), mult. (19658->267), div. (0->0), fcn. (11678->6), ass. (0->85)
t460 = Ifges(4,1) + Ifges(5,1);
t456 = Ifges(4,4) + Ifges(5,4);
t455 = Ifges(4,5) + Ifges(5,5);
t459 = Ifges(4,2) + Ifges(5,2);
t454 = -Ifges(4,6) - Ifges(5,6);
t458 = -Ifges(4,3) - Ifges(5,3);
t436 = qJD(1) ^ 2;
t457 = pkin(2) * t436;
t432 = sin(qJ(1));
t435 = cos(qJ(1));
t424 = -t435 * g(1) - t432 * g(2);
t413 = -t436 * pkin(1) + qJDD(1) * pkin(5) + t424;
t431 = sin(qJ(2));
t453 = t431 * t413;
t434 = cos(qJ(2));
t446 = qJD(1) * qJD(2);
t418 = t431 * qJDD(1) + t434 * t446;
t376 = qJDD(2) * pkin(2) - t418 * pkin(6) - t453 + (pkin(6) * t446 + t431 * t457 - g(3)) * t434;
t399 = -t431 * g(3) + t434 * t413;
t419 = t434 * qJDD(1) - t431 * t446;
t448 = qJD(1) * t431;
t422 = qJD(2) * pkin(2) - pkin(6) * t448;
t429 = t434 ^ 2;
t377 = t419 * pkin(6) - qJD(2) * t422 - t429 * t457 + t399;
t430 = sin(qJ(3));
t433 = cos(qJ(3));
t371 = t433 * t376 - t430 * t377;
t410 = (-t430 * t431 + t433 * t434) * qJD(1);
t384 = t410 * qJD(3) + t433 * t418 + t430 * t419;
t411 = (t430 * t434 + t431 * t433) * qJD(1);
t395 = -t410 * mrSges(5,1) + t411 * mrSges(5,2);
t396 = -t410 * mrSges(4,1) + t411 * mrSges(4,2);
t428 = qJD(2) + qJD(3);
t401 = -t428 * mrSges(4,2) + t410 * mrSges(4,3);
t427 = qJDD(2) + qJDD(3);
t366 = -0.2e1 * qJD(4) * t411 + (t410 * t428 - t384) * qJ(4) + (t410 * t411 + t427) * pkin(3) + t371;
t400 = -t428 * mrSges(5,2) + t410 * mrSges(5,3);
t445 = m(5) * t366 + t427 * mrSges(5,1) + t428 * t400;
t360 = m(4) * t371 + t427 * mrSges(4,1) + t428 * t401 + (-t395 - t396) * t411 + (-mrSges(4,3) - mrSges(5,3)) * t384 + t445;
t372 = t430 * t376 + t433 * t377;
t383 = -t411 * qJD(3) - t430 * t418 + t433 * t419;
t403 = t428 * mrSges(5,1) - t411 * mrSges(5,3);
t404 = t428 * mrSges(4,1) - t411 * mrSges(4,3);
t402 = t428 * pkin(3) - t411 * qJ(4);
t406 = t410 ^ 2;
t368 = -t406 * pkin(3) + t383 * qJ(4) + 0.2e1 * qJD(4) * t410 - t428 * t402 + t372;
t444 = m(5) * t368 + t383 * mrSges(5,3) + t410 * t395;
t363 = m(4) * t372 + t383 * mrSges(4,3) + t410 * t396 + (-t403 - t404) * t428 + (-mrSges(4,2) - mrSges(5,2)) * t427 + t444;
t356 = t433 * t360 + t430 * t363;
t398 = -t434 * g(3) - t453;
t417 = (-mrSges(3,1) * t434 + mrSges(3,2) * t431) * qJD(1);
t447 = qJD(1) * t434;
t421 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t447;
t354 = m(3) * t398 + qJDD(2) * mrSges(3,1) - t418 * mrSges(3,3) + qJD(2) * t421 - t417 * t448 + t356;
t420 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t448;
t441 = -t430 * t360 + t433 * t363;
t355 = m(3) * t399 - qJDD(2) * mrSges(3,2) + t419 * mrSges(3,3) - qJD(2) * t420 + t417 * t447 + t441;
t442 = -t431 * t354 + t434 * t355;
t347 = m(2) * t424 - t436 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t442;
t423 = t432 * g(1) - t435 * g(2);
t439 = -qJDD(1) * pkin(1) - t423;
t412 = -t436 * pkin(5) + t439;
t385 = -t419 * pkin(2) + t422 * t448 + (-pkin(6) * t429 - pkin(5)) * t436 + t439;
t370 = -t383 * pkin(3) - t406 * qJ(4) + t411 * t402 + qJDD(4) + t385;
t440 = m(5) * t370 - t383 * mrSges(5,1) + t384 * mrSges(5,2) - t410 * t400 + t411 * t403;
t438 = m(4) * t385 - t383 * mrSges(4,1) + t384 * mrSges(4,2) - t410 * t401 + t411 * t404 + t440;
t437 = -m(3) * t412 + t419 * mrSges(3,1) - t418 * mrSges(3,2) - t420 * t448 + t421 * t447 - t438;
t358 = m(2) * t423 + qJDD(1) * mrSges(2,1) - t436 * mrSges(2,2) + t437;
t452 = t432 * t347 + t435 * t358;
t348 = t434 * t354 + t431 * t355;
t451 = t454 * t410 - t455 * t411 + t458 * t428;
t450 = -t459 * t410 - t456 * t411 + t454 * t428;
t449 = t456 * t410 + t460 * t411 + t455 * t428;
t443 = t435 * t347 - t432 * t358;
t409 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t431 + Ifges(3,4) * t434) * qJD(1);
t408 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t431 + Ifges(3,2) * t434) * qJD(1);
t407 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t431 + Ifges(3,6) * t434) * qJD(1);
t364 = -t384 * mrSges(5,3) - t411 * t395 + t445;
t350 = mrSges(4,2) * t385 + mrSges(5,2) * t370 - mrSges(4,3) * t371 - mrSges(5,3) * t366 - qJ(4) * t364 + t456 * t383 + t460 * t384 - t451 * t410 + t455 * t427 + t450 * t428;
t349 = -mrSges(4,1) * t385 + mrSges(4,3) * t372 - mrSges(5,1) * t370 + mrSges(5,3) * t368 - pkin(3) * t440 + qJ(4) * t444 + (-qJ(4) * t403 + t449) * t428 + (-qJ(4) * mrSges(5,2) - t454) * t427 + t451 * t411 + t456 * t384 + t459 * t383;
t344 = mrSges(3,2) * t412 - mrSges(3,3) * t398 + Ifges(3,1) * t418 + Ifges(3,4) * t419 + Ifges(3,5) * qJDD(2) - pkin(6) * t356 - qJD(2) * t408 - t430 * t349 + t433 * t350 + t407 * t447;
t343 = -mrSges(3,1) * t412 + mrSges(3,3) * t399 + Ifges(3,4) * t418 + Ifges(3,2) * t419 + Ifges(3,6) * qJDD(2) - pkin(2) * t438 + pkin(6) * t441 + qJD(2) * t409 + t433 * t349 + t430 * t350 - t407 * t448;
t342 = -pkin(2) * t356 - pkin(1) * t348 + t436 * Ifges(2,5) + mrSges(2,3) * t424 - Ifges(3,5) * t418 - Ifges(3,6) * t419 - mrSges(3,1) * t398 + mrSges(3,2) * t399 + mrSges(5,2) * t368 - mrSges(4,1) * t371 + mrSges(4,2) * t372 - pkin(3) * t364 - mrSges(5,1) * t366 + Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) + t458 * t427 + t450 * t411 + t449 * t410 - t455 * t384 + t454 * t383 + (-t431 * t408 + t434 * t409) * qJD(1);
t341 = -mrSges(2,2) * g(3) - mrSges(2,3) * t423 + Ifges(2,5) * qJDD(1) - t436 * Ifges(2,6) - pkin(5) * t348 - t431 * t343 + t434 * t344;
t1 = [-m(1) * g(1) + t443; -m(1) * g(2) + t452; (-m(1) - m(2)) * g(3) + t348; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t452 + t435 * t341 - t432 * t342; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t443 + t432 * t341 + t435 * t342; -mrSges(1,1) * g(2) + mrSges(2,1) * t423 + mrSges(1,2) * g(1) - mrSges(2,2) * t424 + Ifges(2,3) * qJDD(1) + pkin(1) * t437 + pkin(5) * t442 + t434 * t343 + t431 * t344;];
tauB = t1;
