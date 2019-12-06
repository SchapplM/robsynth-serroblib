% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:27
% EndTime: 2019-12-05 15:40:29
% DurationCPUTime: 0.85s
% Computational Cost: add. (4852->188), mult. (8432->224), div. (0->0), fcn. (3823->6), ass. (0->78)
t429 = Ifges(5,1) + Ifges(6,1);
t421 = Ifges(5,4) - Ifges(6,5);
t420 = (Ifges(6,4) + Ifges(5,5));
t428 = Ifges(5,2) + Ifges(6,3);
t417 = (Ifges(5,6) - Ifges(6,6));
t427 = (Ifges(5,3) + Ifges(6,2));
t389 = sin(pkin(7));
t416 = cos(pkin(7));
t374 = -t416 * g(1) - t389 * g(2);
t386 = -g(3) + qJDD(1);
t391 = sin(qJ(2));
t393 = cos(qJ(2));
t351 = t393 * t374 + t391 * t386;
t426 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t351;
t425 = m(3) + m(4);
t424 = -pkin(2) - pkin(6);
t423 = mrSges(3,1) - mrSges(4,2);
t422 = -mrSges(5,3) - mrSges(6,2);
t419 = Ifges(3,5) - Ifges(4,4);
t418 = (-Ifges(3,6) + Ifges(4,5));
t373 = t389 * g(1) - t416 * g(2);
t390 = sin(qJ(4));
t415 = t390 * t373;
t350 = -t391 * t374 + t393 * t386;
t395 = qJD(2) ^ 2;
t399 = -t395 * qJ(3) + qJDD(3) - t350;
t347 = t424 * qJDD(2) + t399;
t392 = cos(qJ(4));
t343 = t390 * t347 - t392 * t373;
t408 = qJD(2) * qJD(4);
t369 = t390 * qJDD(2) + t392 * t408;
t409 = qJD(2) * t392;
t376 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t409;
t367 = (mrSges(6,1) * t390 - mrSges(6,3) * t392) * qJD(2);
t401 = qJD(2) * (-t367 - (mrSges(5,1) * t390 + mrSges(5,2) * t392) * qJD(2));
t366 = (pkin(4) * t390 - qJ(5) * t392) * qJD(2);
t394 = qJD(4) ^ 2;
t410 = qJD(2) * t390;
t340 = -t394 * pkin(4) + qJDD(4) * qJ(5) + (2 * qJD(5) * qJD(4)) - t366 * t410 + t343;
t377 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t409;
t406 = m(6) * t340 + qJDD(4) * mrSges(6,3) + qJD(4) * t377;
t334 = m(5) * t343 - qJDD(4) * mrSges(5,2) - qJD(4) * t376 + t422 * t369 + t390 * t401 + t406;
t342 = t392 * t347 + t415;
t370 = t392 * qJDD(2) - t390 * t408;
t375 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t410;
t341 = -qJDD(4) * pkin(4) - t394 * qJ(5) - t415 + qJDD(5) + (qJD(2) * t366 - t347) * t392;
t378 = -mrSges(6,2) * t410 + qJD(4) * mrSges(6,3);
t400 = -m(6) * t341 + qJDD(4) * mrSges(6,1) + qJD(4) * t378;
t335 = m(5) * t342 + qJDD(4) * mrSges(5,1) + qJD(4) * t375 + t422 * t370 + t392 * t401 + t400;
t328 = t390 * t334 + t392 * t335;
t349 = -qJDD(2) * pkin(2) + t399;
t398 = -m(4) * t349 + (t395 * mrSges(4,3)) - t328;
t324 = m(3) * t350 - (t395 * mrSges(3,2)) + t423 * qJDD(2) + t398;
t348 = t395 * pkin(2) + t426;
t346 = t424 * t395 - t426;
t338 = t369 * pkin(4) - t370 * qJ(5) + (-0.2e1 * qJD(5) * t392 + (pkin(4) * t392 + qJ(5) * t390) * qJD(4)) * qJD(2) + t346;
t336 = m(6) * t338 + t369 * mrSges(6,1) - t370 * mrSges(6,3) - t377 * t409 + t378 * t410;
t397 = -m(5) * t346 - t369 * mrSges(5,1) - t370 * mrSges(5,2) - t375 * t410 - t376 * t409 - t336;
t396 = -m(4) * t348 + (t395 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t397;
t331 = m(3) * t351 - (t395 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t396;
t403 = -t391 * t324 + t393 * t331;
t319 = m(2) * t374 + t403;
t413 = t392 * t334 - t390 * t335;
t326 = (m(2) + t425) * t373 - t413;
t414 = t389 * t319 + t416 * t326;
t320 = t393 * t324 + t391 * t331;
t412 = -(t417 * qJD(4)) + (t428 * t390 - t421 * t392) * qJD(2);
t411 = (t420 * qJD(4)) + (-t421 * t390 + t429 * t392) * qJD(2);
t404 = t416 * t319 - t389 * t326;
t402 = qJD(2) * (-(t427 * qJD(4)) + (t417 * t390 - t420 * t392) * qJD(2));
t327 = -m(4) * t373 + t413;
t322 = mrSges(5,2) * t346 + mrSges(6,2) * t341 - mrSges(5,3) * t342 - mrSges(6,3) * t338 - qJ(5) * t336 + t412 * qJD(4) + t420 * qJDD(4) - t421 * t369 + t429 * t370 + t390 * t402;
t321 = -mrSges(5,1) * t346 - mrSges(6,1) * t338 + mrSges(6,2) * t340 + mrSges(5,3) * t343 - pkin(4) * t336 + t411 * qJD(4) + t417 * qJDD(4) - t428 * t369 + t421 * t370 + t392 * t402;
t316 = -qJ(3) * t327 + pkin(3) * t328 + mrSges(6,3) * t340 - mrSges(6,1) * t341 + mrSges(5,1) * t342 - mrSges(5,2) * t343 + mrSges(4,1) * t349 - mrSges(3,3) * t350 + qJ(5) * t406 + pkin(4) * t400 + (t418 * t395) + (-mrSges(3,2) + mrSges(4,3)) * t373 + (-pkin(4) * mrSges(6,2) + t420) * t370 + (-qJ(5) * mrSges(6,2) - t417) * t369 + t427 * qJDD(4) + t419 * qJDD(2) + ((-pkin(4) * t367 - t412) * t392 + (-qJ(5) * t367 + t411) * t390) * qJD(2);
t315 = -mrSges(4,1) * t348 + mrSges(3,3) * t351 - pkin(2) * t327 - pkin(3) * t397 - pkin(6) * t413 - t418 * qJDD(2) - t392 * t321 - t390 * t322 + t423 * t373 + t419 * t395;
t314 = -pkin(1) * t320 + mrSges(2,3) * t374 - pkin(2) * t398 - qJ(3) * t396 - t392 * t322 + t390 * t321 + pkin(6) * t328 - mrSges(3,1) * t350 + mrSges(3,2) * t351 - mrSges(4,2) * t349 + mrSges(4,3) * t348 - mrSges(2,1) * t386 + (pkin(2) * mrSges(4,2) - Ifges(4,1) - Ifges(3,3)) * qJDD(2);
t313 = mrSges(2,2) * t386 - mrSges(2,3) * t373 - pkin(5) * t320 - t391 * t315 + t393 * t316;
t1 = [-m(1) * g(1) + t404; -m(1) * g(2) + t414; -m(1) * g(3) + m(2) * t386 + t320; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t414 + t416 * t313 - t389 * t314; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t404 + t389 * t313 + t416 * t314; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t374 + t391 * t316 + t393 * t315 - pkin(1) * t413 + pkin(5) * t403 + (pkin(1) * t425 + mrSges(2,1)) * t373;];
tauB = t1;
