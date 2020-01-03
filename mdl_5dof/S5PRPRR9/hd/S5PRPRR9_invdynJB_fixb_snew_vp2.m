% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:42
% EndTime: 2019-12-31 17:39:43
% DurationCPUTime: 1.02s
% Computational Cost: add. (11684->170), mult. (16087->208), div. (0->0), fcn. (7836->8), ass. (0->78)
t423 = qJD(2) ^ 2;
t415 = sin(pkin(8));
t416 = cos(pkin(8));
t399 = t415 * g(1) - t416 * g(2);
t400 = -t416 * g(1) - t415 * g(2);
t419 = sin(qJ(2));
t422 = cos(qJ(2));
t385 = t419 * t399 + t422 * t400;
t432 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t385;
t443 = -pkin(2) - pkin(3);
t377 = t443 * t423 + t432;
t384 = t422 * t399 - t419 * t400;
t429 = -t423 * qJ(3) + qJDD(3) - t384;
t380 = t443 * qJDD(2) + t429;
t418 = sin(qJ(4));
t421 = cos(qJ(4));
t374 = t421 * t377 + t418 * t380;
t407 = -qJD(2) + qJD(4);
t405 = t407 ^ 2;
t406 = -qJDD(2) + qJDD(4);
t371 = -(t405 * pkin(4)) + t406 * pkin(7) + t374;
t413 = g(3) - qJDD(1);
t417 = sin(qJ(5));
t420 = cos(qJ(5));
t368 = -t417 * t371 + t420 * t413;
t369 = t420 * t371 + t417 * t413;
t387 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t417 + Ifges(6,2) * t420) * t407;
t388 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t417 + Ifges(6,4) * t420) * t407;
t436 = qJD(5) * t407;
t392 = t417 * t406 + t420 * t436;
t393 = t420 * t406 - t417 * t436;
t444 = mrSges(6,1) * t368 - mrSges(6,2) * t369 + Ifges(6,5) * t392 + Ifges(6,6) * t393 + Ifges(6,3) * qJDD(5) + (t387 * t417 - t388 * t420) * t407;
t442 = -mrSges(3,1) - mrSges(4,1);
t441 = Ifges(4,4) + Ifges(3,5);
t440 = Ifges(3,6) - Ifges(4,6);
t439 = t407 * t417;
t438 = t407 * t420;
t381 = -t423 * pkin(2) + t432;
t391 = (-mrSges(6,1) * t420 + mrSges(6,2) * t417) * t407;
t398 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t438;
t366 = m(6) * t368 + qJDD(5) * mrSges(6,1) - t392 * mrSges(6,3) + qJD(5) * t398 - t391 * t439;
t397 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t439;
t367 = m(6) * t369 - qJDD(5) * mrSges(6,2) + t393 * mrSges(6,3) - qJD(5) * t397 + t391 * t438;
t360 = -t417 * t366 + t420 * t367;
t357 = m(5) * t374 - (t405 * mrSges(5,1)) - t406 * mrSges(5,2) + t360;
t373 = -t418 * t377 + t421 * t380;
t370 = -t406 * pkin(4) - t405 * pkin(7) - t373;
t364 = -m(6) * t370 + t393 * mrSges(6,1) - t392 * mrSges(6,2) - t397 * t439 + t398 * t438;
t363 = m(5) * t373 + t406 * mrSges(5,1) - t405 * mrSges(5,2) + t364;
t433 = t421 * t357 - t418 * t363;
t430 = m(4) * t381 + qJDD(2) * mrSges(4,3) + t433;
t349 = m(3) * t385 - qJDD(2) * mrSges(3,2) + t442 * t423 + t430;
t355 = t418 * t357 + t421 * t363;
t382 = -qJDD(2) * pkin(2) + t429;
t353 = m(4) * t382 - qJDD(2) * mrSges(4,1) - t423 * mrSges(4,3) + t355;
t350 = m(3) * t384 + qJDD(2) * mrSges(3,1) - t423 * mrSges(3,2) - t353;
t344 = t419 * t349 + t422 * t350;
t342 = m(2) * t399 + t344;
t434 = t422 * t349 - t419 * t350;
t343 = m(2) * t400 + t434;
t437 = t416 * t342 + t415 * t343;
t435 = -t415 * t342 + t416 * t343;
t359 = t420 * t366 + t417 * t367;
t358 = -t359 + (-m(4) - m(5)) * t413;
t428 = -m(3) * t413 + t358;
t427 = -m(2) * t413 + t428;
t386 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t417 + Ifges(6,6) * t420) * t407;
t361 = -mrSges(6,1) * t370 + mrSges(6,3) * t369 + Ifges(6,4) * t392 + Ifges(6,2) * t393 + Ifges(6,6) * qJDD(5) + qJD(5) * t388 - t386 * t439;
t362 = mrSges(6,2) * t370 - mrSges(6,3) * t368 + Ifges(6,1) * t392 + Ifges(6,4) * t393 + Ifges(6,5) * qJDD(5) - qJD(5) * t387 + t386 * t438;
t425 = mrSges(5,1) * t373 - mrSges(5,2) * t374 + Ifges(5,3) * t406 + pkin(4) * t364 + pkin(7) * t360 + t420 * t361 + t417 * t362;
t424 = -mrSges(4,1) * t382 - mrSges(3,2) * t385 - pkin(3) * t355 + qJ(3) * (-t423 * mrSges(4,1) + t430) - pkin(2) * t353 + mrSges(4,3) * t381 + mrSges(3,1) * t384 - t425 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2);
t354 = -mrSges(5,1) * t413 + mrSges(5,3) * t374 + (t405 * Ifges(5,5)) + Ifges(5,6) * t406 - pkin(4) * t359 - t444;
t345 = mrSges(5,2) * t413 - mrSges(5,3) * t373 + Ifges(5,5) * t406 - t405 * Ifges(5,6) - pkin(7) * t359 - t417 * t361 + t420 * t362;
t338 = mrSges(4,2) * t382 - mrSges(3,3) * t384 - pkin(6) * t355 - qJ(3) * t358 + t421 * t345 - t418 * t354 - t440 * t423 + (-mrSges(3,2) + mrSges(4,3)) * t413 + t441 * qJDD(2);
t337 = mrSges(3,3) * t385 + mrSges(4,2) * t381 - t418 * t345 - t421 * t354 + pkin(3) * t359 - pkin(6) * t433 - pkin(2) * t358 + t441 * t423 + (pkin(3) * m(5) - t442) * t413 + t440 * qJDD(2);
t336 = -mrSges(2,2) * t413 - mrSges(2,3) * t399 - pkin(5) * t344 - t419 * t337 + t422 * t338;
t335 = mrSges(2,1) * t413 + mrSges(2,3) * t400 - pkin(1) * t428 + pkin(5) * t434 + t422 * t337 + t419 * t338;
t1 = [-m(1) * g(1) + t435; -m(1) * g(2) + t437; -m(1) * g(3) + t427; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t437 - t415 * t335 + t416 * t336; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t435 + t416 * t335 + t415 * t336; -mrSges(1,1) * g(2) + mrSges(2,1) * t399 + mrSges(1,2) * g(1) - mrSges(2,2) * t400 + pkin(1) * t344 + t424; t427; t424; t353; t425; t444;];
tauJB = t1;
