% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:07
% EndTime: 2019-12-05 17:38:09
% DurationCPUTime: 1.03s
% Computational Cost: add. (7053->216), mult. (13414->256), div. (0->0), fcn. (6523->6), ass. (0->80)
t471 = 2 * qJD(1);
t443 = sin(qJ(1));
t446 = cos(qJ(1));
t422 = -t446 * g(1) - t443 * g(2);
t470 = qJDD(1) * qJ(2) + (qJD(2) * t471) + t422;
t421 = t443 * g(1) - t446 * g(2);
t447 = qJD(1) ^ 2;
t405 = -qJDD(1) * pkin(1) - t447 * qJ(2) + qJDD(2) - t421;
t451 = qJDD(1) * qJ(3) + (qJD(3) * t471) - t405;
t469 = -m(3) - m(4);
t468 = mrSges(2,1) - mrSges(3,2);
t404 = t447 * pkin(1) - t470;
t401 = qJDD(3) + (-pkin(1) - qJ(3)) * t447 + t470;
t398 = -qJDD(1) * pkin(6) + t401;
t442 = sin(qJ(4));
t445 = cos(qJ(4));
t391 = t442 * g(3) + t445 * t398;
t464 = qJD(1) * qJD(4);
t459 = t442 * t464;
t417 = t445 * qJDD(1) - t459;
t378 = (-t417 - t459) * pkin(7) + (-t442 * t445 * t447 + qJDD(4)) * pkin(4) + t391;
t392 = -t445 * g(3) + t442 * t398;
t416 = -t442 * qJDD(1) - t445 * t464;
t465 = qJD(1) * t445;
t420 = qJD(4) * pkin(4) - pkin(7) * t465;
t436 = t442 ^ 2;
t379 = -t436 * t447 * pkin(4) + t416 * pkin(7) - qJD(4) * t420 + t392;
t441 = sin(qJ(5));
t444 = cos(qJ(5));
t376 = t444 * t378 - t441 * t379;
t409 = (-t441 * t445 - t442 * t444) * qJD(1);
t385 = t409 * qJD(5) + t441 * t416 + t444 * t417;
t410 = (-t441 * t442 + t444 * t445) * qJD(1);
t393 = -t409 * mrSges(6,1) + t410 * mrSges(6,2);
t429 = qJD(4) + qJD(5);
t402 = -t429 * mrSges(6,2) + t409 * mrSges(6,3);
t428 = qJDD(4) + qJDD(5);
t374 = m(6) * t376 + t428 * mrSges(6,1) - t385 * mrSges(6,3) - t410 * t393 + t429 * t402;
t377 = t441 * t378 + t444 * t379;
t384 = -t410 * qJD(5) + t444 * t416 - t441 * t417;
t403 = t429 * mrSges(6,1) - t410 * mrSges(6,3);
t375 = m(6) * t377 - t428 * mrSges(6,2) + t384 * mrSges(6,3) + t409 * t393 - t429 * t403;
t365 = t444 * t374 + t441 * t375;
t415 = (mrSges(5,1) * t442 + mrSges(5,2) * t445) * qJD(1);
t466 = qJD(1) * t442;
t418 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t466;
t363 = m(5) * t391 + qJDD(4) * mrSges(5,1) - t417 * mrSges(5,3) + qJD(4) * t418 - t415 * t465 + t365;
t419 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t465;
t456 = -t441 * t374 + t444 * t375;
t364 = m(5) * t392 - qJDD(4) * mrSges(5,2) + t416 * mrSges(5,3) - qJD(4) * t419 - t415 * t466 + t456;
t359 = t445 * t363 + t442 * t364;
t455 = -m(4) * t401 - qJDD(1) * mrSges(4,2) - t359;
t450 = -m(3) * t404 + (t447 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t455;
t357 = m(2) * t422 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t447) + t450;
t397 = -t447 * pkin(6) + t451;
t381 = t420 * t465 - t416 * pkin(4) + (-pkin(7) * t436 - pkin(6)) * t447 + t451;
t454 = m(6) * t381 - t384 * mrSges(6,1) + t385 * mrSges(6,2) - t409 * t402 + t410 * t403;
t449 = -m(5) * t397 + t416 * mrSges(5,1) - t417 * mrSges(5,2) - t418 * t466 - t419 * t465 - t454;
t370 = -m(4) * t451 - t447 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t449;
t448 = -m(3) * t405 + t447 * mrSges(3,3) - t370;
t369 = m(2) * t421 - t447 * mrSges(2,2) + t468 * qJDD(1) + t448;
t467 = t443 * t357 + t446 * t369;
t461 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t460 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t458 = t446 * t357 - t443 * t369;
t457 = -t442 * t363 + t445 * t364;
t408 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t445 - Ifges(5,4) * t442) * qJD(1);
t407 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t445 - Ifges(5,2) * t442) * qJD(1);
t406 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t445 - Ifges(5,6) * t442) * qJD(1);
t388 = Ifges(6,1) * t410 + Ifges(6,4) * t409 + Ifges(6,5) * t429;
t387 = Ifges(6,4) * t410 + Ifges(6,2) * t409 + Ifges(6,6) * t429;
t386 = Ifges(6,5) * t410 + Ifges(6,6) * t409 + Ifges(6,3) * t429;
t367 = mrSges(6,2) * t381 - mrSges(6,3) * t376 + Ifges(6,1) * t385 + Ifges(6,4) * t384 + Ifges(6,5) * t428 + t409 * t386 - t429 * t387;
t366 = -mrSges(6,1) * t381 + mrSges(6,3) * t377 + Ifges(6,4) * t385 + Ifges(6,2) * t384 + Ifges(6,6) * t428 - t410 * t386 + t429 * t388;
t358 = t469 * g(3) + t457;
t354 = mrSges(5,2) * t397 - mrSges(5,3) * t391 + Ifges(5,1) * t417 + Ifges(5,4) * t416 + Ifges(5,5) * qJDD(4) - pkin(7) * t365 - qJD(4) * t407 - t441 * t366 + t444 * t367 - t406 * t466;
t353 = -mrSges(5,1) * t397 + mrSges(5,3) * t392 + Ifges(5,4) * t417 + Ifges(5,2) * t416 + Ifges(5,6) * qJDD(4) - pkin(4) * t454 + pkin(7) * t456 + qJD(4) * t408 + t444 * t366 + t441 * t367 - t406 * t465;
t352 = Ifges(5,3) * qJDD(4) - qJ(3) * t457 - pkin(2) * t455 + mrSges(2,3) * t422 + Ifges(6,3) * t428 + Ifges(5,6) * t416 + Ifges(5,5) * t417 - mrSges(3,1) * t404 - t409 * t388 + t410 * t387 + mrSges(4,1) * t401 + Ifges(6,6) * t384 + Ifges(6,5) * t385 + mrSges(5,1) * t391 - mrSges(5,2) * t392 + mrSges(6,1) * t376 - mrSges(6,2) * t377 + pkin(4) * t365 + pkin(3) * t359 - pkin(1) * t358 + (t445 * t407 + t442 * t408) * qJD(1) + (-pkin(2) * mrSges(4,3) + t461) * t447 + t460 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t468) * g(3);
t351 = -qJ(2) * t358 - mrSges(2,3) * t421 + pkin(2) * t370 + mrSges(3,1) * t405 + t445 * t353 + pkin(3) * t449 + pkin(6) * t457 + t442 * t354 - mrSges(4,1) * t451 - t460 * t447 + t461 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t458; -m(1) * g(2) + t467; (-m(1) - m(2) + t469) * g(3) + t457; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t467 + t446 * t351 - t443 * t352; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t458 + t443 * t351 + t446 * t352; qJ(2) * (-(t447 * mrSges(4,3)) + t450) + pkin(1) * t448 + mrSges(2,1) * t421 - mrSges(2,2) * t422 - qJ(3) * t370 + mrSges(3,2) * t405 - mrSges(3,3) * t404 + t445 * t354 - t442 * t353 - pkin(6) * t359 + mrSges(4,3) * t451 + mrSges(4,2) * t401 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
