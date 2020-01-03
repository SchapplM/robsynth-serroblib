% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR6
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:36
% EndTime: 2019-12-31 16:52:38
% DurationCPUTime: 1.99s
% Computational Cost: add. (18486->217), mult. (44608->276), div. (0->0), fcn. (30738->8), ass. (0->94)
t428 = qJD(1) ^ 2;
t421 = cos(pkin(7));
t451 = pkin(2) * t421;
t420 = sin(pkin(7));
t450 = mrSges(3,2) * t420;
t418 = t421 ^ 2;
t449 = t418 * t428;
t424 = sin(qJ(1));
t427 = cos(qJ(1));
t409 = -t427 * g(1) - t424 * g(2);
t405 = -t428 * pkin(1) + qJDD(1) * qJ(2) + t409;
t445 = qJD(1) * qJD(2);
t443 = -t421 * g(3) - 0.2e1 * t420 * t445;
t382 = (-pkin(5) * qJDD(1) + t428 * t451 - t405) * t420 + t443;
t396 = -t420 * g(3) + (t405 + 0.2e1 * t445) * t421;
t444 = qJDD(1) * t421;
t383 = -pkin(2) * t449 + pkin(5) * t444 + t396;
t423 = sin(qJ(3));
t426 = cos(qJ(3));
t369 = t426 * t382 - t423 * t383;
t433 = t420 * t426 + t421 * t423;
t432 = -t420 * t423 + t421 * t426;
t403 = t432 * qJD(1);
t446 = t403 * qJD(3);
t394 = t433 * qJDD(1) + t446;
t404 = t433 * qJD(1);
t361 = (-t394 + t446) * pkin(6) + (t403 * t404 + qJDD(3)) * pkin(3) + t369;
t370 = t423 * t382 + t426 * t383;
t393 = -t404 * qJD(3) + t432 * qJDD(1);
t399 = qJD(3) * pkin(3) - t404 * pkin(6);
t402 = t403 ^ 2;
t362 = -t402 * pkin(3) + t393 * pkin(6) - qJD(3) * t399 + t370;
t422 = sin(qJ(4));
t425 = cos(qJ(4));
t359 = t425 * t361 - t422 * t362;
t387 = t425 * t403 - t422 * t404;
t368 = t387 * qJD(4) + t422 * t393 + t425 * t394;
t388 = t422 * t403 + t425 * t404;
t376 = -t387 * mrSges(5,1) + t388 * mrSges(5,2);
t419 = qJD(3) + qJD(4);
t379 = -t419 * mrSges(5,2) + t387 * mrSges(5,3);
t416 = qJDD(3) + qJDD(4);
t357 = m(5) * t359 + t416 * mrSges(5,1) - t368 * mrSges(5,3) - t388 * t376 + t419 * t379;
t360 = t422 * t361 + t425 * t362;
t367 = -t388 * qJD(4) + t425 * t393 - t422 * t394;
t380 = t419 * mrSges(5,1) - t388 * mrSges(5,3);
t358 = m(5) * t360 - t416 * mrSges(5,2) + t367 * mrSges(5,3) + t387 * t376 - t419 * t380;
t349 = t425 * t357 + t422 * t358;
t390 = -t403 * mrSges(4,1) + t404 * mrSges(4,2);
t397 = -qJD(3) * mrSges(4,2) + t403 * mrSges(4,3);
t347 = m(4) * t369 + qJDD(3) * mrSges(4,1) - t394 * mrSges(4,3) + qJD(3) * t397 - t404 * t390 + t349;
t398 = qJD(3) * mrSges(4,1) - t404 * mrSges(4,3);
t439 = -t422 * t357 + t425 * t358;
t348 = m(4) * t370 - qJDD(3) * mrSges(4,2) + t393 * mrSges(4,3) - qJD(3) * t398 + t403 * t390 + t439;
t343 = t426 * t347 + t423 * t348;
t395 = -t420 * t405 + t443;
t431 = mrSges(3,3) * qJDD(1) + t428 * (-mrSges(3,1) * t421 + t450);
t341 = m(3) * t395 - t431 * t420 + t343;
t440 = -t423 * t347 + t426 * t348;
t342 = m(3) * t396 + t431 * t421 + t440;
t441 = -t420 * t341 + t421 * t342;
t334 = m(2) * t409 - t428 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t441;
t408 = t424 * g(1) - t427 * g(2);
t438 = qJDD(2) - t408;
t401 = -qJDD(1) * pkin(1) - t428 * qJ(2) + t438;
t417 = t420 ^ 2;
t392 = (-pkin(1) - t451) * qJDD(1) + (-qJ(2) + (-t417 - t418) * pkin(5)) * t428 + t438;
t364 = -t393 * pkin(3) - t402 * pkin(6) + t404 * t399 + t392;
t434 = m(5) * t364 - t367 * mrSges(5,1) + t368 * mrSges(5,2) - t387 * t379 + t388 * t380;
t430 = m(4) * t392 - t393 * mrSges(4,1) + t394 * mrSges(4,2) - t403 * t397 + t404 * t398 + t434;
t429 = -m(3) * t401 + mrSges(3,1) * t444 - t430 + (t417 * t428 + t449) * mrSges(3,3);
t353 = m(2) * t408 - t428 * mrSges(2,2) + t429 + (mrSges(2,1) - t450) * qJDD(1);
t448 = t424 * t334 + t427 * t353;
t335 = t421 * t341 + t420 * t342;
t435 = Ifges(3,5) * t420 + Ifges(3,6) * t421;
t447 = t428 * t435;
t442 = t427 * t334 - t424 * t353;
t437 = Ifges(3,1) * t420 + Ifges(3,4) * t421;
t436 = Ifges(3,4) * t420 + Ifges(3,2) * t421;
t386 = Ifges(4,1) * t404 + Ifges(4,4) * t403 + Ifges(4,5) * qJD(3);
t385 = Ifges(4,4) * t404 + Ifges(4,2) * t403 + Ifges(4,6) * qJD(3);
t384 = Ifges(4,5) * t404 + Ifges(4,6) * t403 + Ifges(4,3) * qJD(3);
t373 = Ifges(5,1) * t388 + Ifges(5,4) * t387 + Ifges(5,5) * t419;
t372 = Ifges(5,4) * t388 + Ifges(5,2) * t387 + Ifges(5,6) * t419;
t371 = Ifges(5,5) * t388 + Ifges(5,6) * t387 + Ifges(5,3) * t419;
t351 = mrSges(5,2) * t364 - mrSges(5,3) * t359 + Ifges(5,1) * t368 + Ifges(5,4) * t367 + Ifges(5,5) * t416 + t387 * t371 - t419 * t372;
t350 = -mrSges(5,1) * t364 + mrSges(5,3) * t360 + Ifges(5,4) * t368 + Ifges(5,2) * t367 + Ifges(5,6) * t416 - t388 * t371 + t419 * t373;
t337 = mrSges(4,2) * t392 - mrSges(4,3) * t369 + Ifges(4,1) * t394 + Ifges(4,4) * t393 + Ifges(4,5) * qJDD(3) - pkin(6) * t349 - qJD(3) * t385 - t422 * t350 + t425 * t351 + t403 * t384;
t336 = -mrSges(4,1) * t392 + mrSges(4,3) * t370 + Ifges(4,4) * t394 + Ifges(4,2) * t393 + Ifges(4,6) * qJDD(3) - pkin(3) * t434 + pkin(6) * t439 + qJD(3) * t386 + t425 * t350 + t422 * t351 - t404 * t384;
t331 = mrSges(3,2) * t401 - mrSges(3,3) * t395 - pkin(5) * t343 + t437 * qJDD(1) - t423 * t336 + t426 * t337 + t421 * t447;
t330 = -mrSges(3,1) * t401 + mrSges(3,3) * t396 - pkin(2) * t430 + pkin(5) * t440 + t436 * qJDD(1) + t426 * t336 + t423 * t337 - t420 * t447;
t329 = -pkin(1) * t335 + t403 * t386 - t404 * t385 + mrSges(2,3) * t409 - Ifges(5,3) * t416 - pkin(2) * t343 - pkin(3) * t349 + t387 * t373 - t388 * t372 - Ifges(4,6) * t393 - Ifges(4,5) * t394 - mrSges(3,1) * t395 + mrSges(3,2) * t396 - Ifges(4,3) * qJDD(3) + mrSges(2,1) * g(3) - mrSges(5,1) * t359 + mrSges(5,2) * t360 - Ifges(5,6) * t367 - Ifges(5,5) * t368 - mrSges(4,1) * t369 + mrSges(4,2) * t370 + (Ifges(2,6) - t435) * qJDD(1) + (-t420 * t436 + t421 * t437 + Ifges(2,5)) * t428;
t328 = -mrSges(2,2) * g(3) - mrSges(2,3) * t408 + Ifges(2,5) * qJDD(1) - t428 * Ifges(2,6) - qJ(2) * t335 - t420 * t330 + t421 * t331;
t1 = [-m(1) * g(1) + t442; -m(1) * g(2) + t448; (-m(1) - m(2)) * g(3) + t335; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t448 + t427 * t328 - t424 * t329; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t442 + t424 * t328 + t427 * t329; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t408 - mrSges(2,2) * t409 + t420 * t331 + t421 * t330 + pkin(1) * (-qJDD(1) * t450 + t429) + qJ(2) * t441;];
tauB = t1;
