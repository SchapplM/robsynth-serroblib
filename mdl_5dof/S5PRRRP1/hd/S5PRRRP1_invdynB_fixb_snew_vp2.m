% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:59
% EndTime: 2019-12-05 16:40:01
% DurationCPUTime: 1.61s
% Computational Cost: add. (16875->198), mult. (22715->241), div. (0->0), fcn. (12664->8), ass. (0->89)
t458 = Ifges(5,1) + Ifges(6,1);
t451 = Ifges(5,4) + Ifges(6,4);
t450 = Ifges(5,5) + Ifges(6,5);
t457 = Ifges(5,2) + Ifges(6,2);
t456 = Ifges(5,6) + Ifges(6,6);
t455 = Ifges(5,3) + Ifges(6,3);
t415 = qJD(2) + qJD(3);
t413 = t415 ^ 2;
t419 = sin(pkin(8));
t420 = cos(pkin(8));
t407 = g(1) * t419 - g(2) * t420;
t408 = -g(1) * t420 - g(2) * t419;
t423 = sin(qJ(2));
t426 = cos(qJ(2));
t380 = t426 * t407 - t408 * t423;
t378 = qJDD(2) * pkin(2) + t380;
t381 = t423 * t407 + t426 * t408;
t427 = qJD(2) ^ 2;
t379 = -pkin(2) * t427 + t381;
t422 = sin(qJ(3));
t425 = cos(qJ(3));
t373 = t378 * t425 - t422 * t379;
t414 = qJDD(2) + qJDD(3);
t429 = -pkin(3) * t414 - t373;
t371 = -pkin(7) * t413 + t429;
t421 = sin(qJ(4));
t424 = cos(qJ(4));
t441 = qJD(4) * t415;
t436 = t424 * t441;
t396 = t414 * t421 + t436;
t397 = t414 * t424 - t421 * t441;
t447 = t415 * t424;
t406 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t447;
t448 = t415 * t421;
t402 = qJD(4) * pkin(4) - qJ(5) * t448;
t417 = t424 ^ 2;
t367 = t402 * t448 - pkin(4) * t397 + qJDD(5) + (-qJ(5) * t417 - pkin(7)) * t413 + t429;
t405 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t447;
t430 = m(6) * t367 - t397 * mrSges(6,1) - t405 * t447;
t452 = -mrSges(5,2) - mrSges(6,2);
t454 = -m(5) * t371 + t397 * mrSges(5,1) + t452 * t396 + t406 * t447 - t430;
t453 = pkin(4) * t413;
t374 = t422 * t378 + t425 * t379;
t372 = -pkin(3) * t413 + pkin(7) * t414 + t374;
t418 = -g(3) + qJDD(1);
t369 = t424 * t372 + t421 * t418;
t395 = (-mrSges(5,1) * t424 + mrSges(5,2) * t421) * t415;
t440 = qJD(5) * t415;
t366 = qJ(5) * t397 - qJD(4) * t402 - t417 * t453 + 0.2e1 * t424 * t440 + t369;
t394 = (-mrSges(6,1) * t424 + mrSges(6,2) * t421) * t415;
t438 = m(6) * t366 + t397 * mrSges(6,3) + t394 * t447;
t403 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t448;
t442 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t448 - t403;
t361 = m(5) * t369 + mrSges(5,3) * t397 + t442 * qJD(4) + t452 * qJDD(4) + t395 * t447 + t438;
t359 = t424 * t361;
t410 = t424 * t418;
t368 = -t372 * t421 + t410;
t365 = qJDD(4) * pkin(4) + t410 + (-t396 + t436) * qJ(5) + (t424 * t453 - t372 - 0.2e1 * t440) * t421;
t439 = m(6) * t365 + qJDD(4) * mrSges(6,1) + qJD(4) * t405;
t360 = m(5) * t368 + qJDD(4) * mrSges(5,1) + qJD(4) * t406 + (-t394 - t395) * t448 + (-mrSges(5,3) - mrSges(6,3)) * t396 + t439;
t353 = m(4) * t374 - mrSges(4,1) * t413 - mrSges(4,2) * t414 - t360 * t421 + t359;
t435 = t415 * t442;
t356 = m(4) * t373 + mrSges(4,1) * t414 - mrSges(4,2) * t413 + t421 * t435 + t454;
t348 = t422 * t353 + t425 * t356;
t346 = m(3) * t380 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t427 + t348;
t432 = t425 * t353 - t356 * t422;
t347 = m(3) * t381 - mrSges(3,1) * t427 - qJDD(2) * mrSges(3,2) + t432;
t341 = t426 * t346 + t423 * t347;
t339 = m(2) * t407 + t341;
t433 = -t423 * t346 + t426 * t347;
t340 = m(2) * t408 + t433;
t446 = t420 * t339 + t419 * t340;
t354 = t424 * t360 + t421 * t361;
t445 = (t450 * t421 + t456 * t424) * t415 + t455 * qJD(4);
t444 = (-t451 * t421 - t457 * t424) * t415 - t456 * qJD(4);
t443 = (t458 * t421 + t451 * t424) * t415 + t450 * qJD(4);
t437 = m(4) * t418 + t354;
t434 = -t339 * t419 + t420 * t340;
t431 = m(3) * t418 + t437;
t362 = -mrSges(6,3) * t396 - t394 * t448 + t439;
t350 = mrSges(5,2) * t371 + mrSges(6,2) * t367 - mrSges(5,3) * t368 - mrSges(6,3) * t365 - qJ(5) * t362 + t444 * qJD(4) + t450 * qJDD(4) + t458 * t396 + t451 * t397 + t445 * t447;
t349 = -mrSges(5,1) * t371 + mrSges(5,3) * t369 - mrSges(6,1) * t367 + mrSges(6,3) * t366 - pkin(4) * t430 + qJ(5) * t438 + (-pkin(4) * t403 - t445) * t448 + t457 * t397 + (-mrSges(6,2) * pkin(4) + t451) * t396 + (-mrSges(6,2) * qJ(5) + t456) * qJDD(4) + (-qJ(5) * t403 + t443) * qJD(4);
t342 = -mrSges(4,1) * t418 - mrSges(5,1) * t368 - mrSges(6,1) * t365 + mrSges(5,2) * t369 + mrSges(6,2) * t366 + mrSges(4,3) * t374 + t413 * Ifges(4,5) + Ifges(4,6) * t414 - pkin(3) * t354 - pkin(4) * t362 - t456 * t397 - t450 * t396 - t455 * qJDD(4) + (t444 * t421 + t443 * t424) * t415;
t335 = mrSges(4,2) * t418 - mrSges(4,3) * t373 + Ifges(4,5) * t414 - Ifges(4,6) * t413 - pkin(7) * t354 - t349 * t421 + t350 * t424;
t334 = mrSges(3,2) * t418 - mrSges(3,3) * t380 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t427 - pkin(6) * t348 + t335 * t425 - t342 * t422;
t333 = -mrSges(3,1) * t418 + mrSges(3,3) * t381 + t427 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t437 + pkin(6) * t432 + t422 * t335 + t425 * t342;
t332 = mrSges(2,2) * t418 - mrSges(2,3) * t407 - pkin(5) * t341 - t333 * t423 + t334 * t426;
t331 = -mrSges(2,1) * t418 + mrSges(2,3) * t408 - pkin(1) * t431 + pkin(5) * t433 + t426 * t333 + t423 * t334;
t1 = [-m(1) * g(1) + t434; -m(1) * g(2) + t446; -m(1) * g(3) + m(2) * t418 + t431; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t446 - t419 * t331 + t420 * t332; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t434 + t420 * t331 + t419 * t332; pkin(1) * t341 + mrSges(2,1) * t407 - mrSges(2,2) * t408 + pkin(2) * t348 - mrSges(3,2) * t381 + mrSges(3,1) * t380 + t424 * t349 + pkin(3) * t454 + pkin(7) * t359 + mrSges(4,1) * t373 - mrSges(4,2) * t374 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t414 + Ifges(3,3) * qJDD(2) + (pkin(3) * t435 - pkin(7) * t360 + t350) * t421;];
tauB = t1;
