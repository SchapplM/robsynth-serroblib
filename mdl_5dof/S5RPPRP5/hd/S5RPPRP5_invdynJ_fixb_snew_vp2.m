% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:27
% EndTime: 2019-12-31 17:53:28
% DurationCPUTime: 0.87s
% Computational Cost: add. (2226->179), mult. (5271->217), div. (0->0), fcn. (3125->6), ass. (0->81)
t399 = qJD(1) ^ 2;
t395 = sin(qJ(1));
t397 = cos(qJ(1));
t407 = -t397 * g(1) - t395 * g(2);
t441 = -t399 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t407;
t440 = Ifges(5,1) + Ifges(6,1);
t433 = Ifges(5,4) - Ifges(6,5);
t432 = Ifges(5,5) + Ifges(6,4);
t439 = -Ifges(5,2) - Ifges(6,3);
t431 = Ifges(5,6) - Ifges(6,6);
t438 = Ifges(5,3) + Ifges(6,2);
t437 = mrSges(4,2) + mrSges(3,3);
t393 = cos(pkin(7));
t387 = t393 ^ 2;
t392 = sin(pkin(7));
t411 = (-t392 ^ 2 - t387) * t399;
t422 = t395 * g(1) - t397 * g(2);
t368 = -qJDD(1) * pkin(1) - t399 * qJ(2) + qJDD(2) - t422;
t415 = qJDD(1) * t393;
t420 = qJD(1) * t392;
t430 = t392 * qJ(3);
t353 = -pkin(2) * t415 - 0.2e1 * qJD(3) * t420 - qJDD(1) * t430 + t368;
t356 = -t393 * g(3) - t441 * t392;
t436 = pkin(3) * t399;
t435 = -mrSges(5,3) - mrSges(6,2);
t434 = Ifges(3,4) - Ifges(4,5);
t372 = (-t393 * pkin(2) - t430) * qJD(1);
t336 = t372 * t420 + qJDD(3) - t356;
t332 = (-pkin(6) * qJDD(1) - t393 * t436) * t392 + t336;
t357 = -t392 * g(3) + t441 * t393;
t419 = qJD(1) * t393;
t337 = t372 * t419 + t357;
t334 = -pkin(6) * t415 - t387 * t436 + t337;
t394 = sin(qJ(4));
t396 = cos(qJ(4));
t330 = t394 * t332 + t396 * t334;
t404 = t392 * t394 + t393 * t396;
t405 = t392 * t396 - t393 * t394;
t370 = t405 * qJD(1);
t417 = t370 * qJD(4);
t354 = t404 * qJDD(1) + t417;
t361 = qJD(4) * mrSges(5,1) - t370 * mrSges(5,3);
t369 = t404 * qJD(1);
t347 = t369 * pkin(4) - t370 * qJ(5);
t398 = qJD(4) ^ 2;
t325 = -t398 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t369 * t347 + t330;
t362 = -qJD(4) * mrSges(6,1) + t370 * mrSges(6,2);
t414 = m(6) * t325 + qJDD(4) * mrSges(6,3) + qJD(4) * t362;
t348 = t369 * mrSges(6,1) - t370 * mrSges(6,3);
t425 = -t369 * mrSges(5,1) - t370 * mrSges(5,2) - t348;
t320 = m(5) * t330 - qJDD(4) * mrSges(5,2) - qJD(4) * t361 + t435 * t354 + t425 * t369 + t414;
t329 = t396 * t332 - t394 * t334;
t418 = t369 * qJD(4);
t355 = t405 * qJDD(1) - t418;
t360 = -qJD(4) * mrSges(5,2) - t369 * mrSges(5,3);
t326 = -qJDD(4) * pkin(4) - t398 * qJ(5) + t370 * t347 + qJDD(5) - t329;
t363 = -t369 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t408 = -m(6) * t326 + qJDD(4) * mrSges(6,1) + qJD(4) * t363;
t321 = m(5) * t329 + qJDD(4) * mrSges(5,1) + qJD(4) * t360 + t435 * t355 + t425 * t370 + t408;
t429 = t394 * t320 + t396 * t321;
t428 = -t438 * qJD(4) + t431 * t369 - t432 * t370;
t427 = t431 * qJD(4) + t439 * t369 + t433 * t370;
t426 = t432 * qJD(4) - t433 * t369 + t440 * t370;
t423 = ((Ifges(3,6) - Ifges(4,6)) * t393 + (Ifges(4,4) + Ifges(3,5)) * t392) * qJD(1);
t410 = t396 * t320 - t394 * t321;
t409 = m(4) * t336 + t429;
t406 = -t393 * mrSges(4,1) - t392 * mrSges(4,3);
t335 = pkin(3) * t415 + pkin(6) * t411 - t353;
t328 = -0.2e1 * qJD(5) * t370 + (-t355 + t418) * qJ(5) + (t354 + t417) * pkin(4) + t335;
t322 = m(6) * t328 + t354 * mrSges(6,1) - t355 * mrSges(6,3) - t370 * t362 + t369 * t363;
t373 = t406 * qJD(1);
t403 = qJDD(1) * mrSges(4,2) + qJD(1) * t373;
t401 = -m(5) * t335 - t354 * mrSges(5,1) - t355 * mrSges(5,2) - t369 * t360 - t370 * t361 - t322;
t400 = m(4) * t353 + t401;
t374 = (-t393 * mrSges(3,1) + t392 * mrSges(3,2)) * qJD(1);
t323 = t355 * mrSges(6,2) + t370 * t348 - t408;
t316 = mrSges(4,2) * t411 + t406 * qJDD(1) + t400;
t315 = t400 + ((-mrSges(3,1) - mrSges(4,1)) * t393 + (mrSges(3,2) - mrSges(4,3)) * t392) * qJDD(1) + m(3) * t368 + t437 * t411;
t314 = mrSges(5,2) * t335 + mrSges(6,2) * t326 - mrSges(5,3) * t329 - mrSges(6,3) * t328 - qJ(5) * t322 - t427 * qJD(4) + t432 * qJDD(4) - t433 * t354 + t440 * t355 + t428 * t369;
t313 = -mrSges(5,1) * t335 - mrSges(6,1) * t328 + mrSges(6,2) * t325 + mrSges(5,3) * t330 - pkin(4) * t322 + t426 * qJD(4) + t431 * qJDD(4) + t439 * t354 + t433 * t355 + t428 * t370;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t422 - mrSges(2,2) * t407 + t392 * (mrSges(3,2) * t368 - mrSges(3,3) * t356 + mrSges(4,2) * t336 - mrSges(4,3) * t353 + t396 * t314 - t394 * t313 - pkin(6) * t429 - qJ(3) * t316 + t423 * t419 + (t393 * t434 + (Ifges(3,1) + Ifges(4,1)) * t392) * qJDD(1)) + t393 * (-mrSges(3,1) * t368 + mrSges(3,3) * t357 - mrSges(4,1) * t353 + mrSges(4,2) * t337 - t394 * t314 - t396 * t313 - pkin(3) * t401 - pkin(6) * t410 - pkin(2) * t316 - t423 * t420 + ((Ifges(3,2) + Ifges(4,3)) * t393 + t392 * t434) * qJDD(1)) - pkin(1) * t315 + qJ(2) * (t393 * (m(3) * t357 + m(4) * t337 + (t437 * qJDD(1) + (t373 + t374) * qJD(1)) * t393 + t410) + (-m(3) * t356 + (qJDD(1) * mrSges(3,3) + qJD(1) * t374 + t403) * t392 + t409) * t392); t315; t403 * t392 + t409; mrSges(5,1) * t329 - mrSges(5,2) * t330 - mrSges(6,1) * t326 + mrSges(6,3) * t325 - pkin(4) * t323 + qJ(5) * t414 + t427 * t370 + (-qJ(5) * t348 + t426) * t369 + t432 * t355 + (-qJ(5) * mrSges(6,2) - t431) * t354 + t438 * qJDD(4); t323;];
tauJ = t1;
