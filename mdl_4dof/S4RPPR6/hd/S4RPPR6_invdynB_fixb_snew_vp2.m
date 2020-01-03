% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:41
% EndTime: 2019-12-31 16:40:42
% DurationCPUTime: 0.86s
% Computational Cost: add. (4293->183), mult. (9911->228), div. (0->0), fcn. (5547->6), ass. (0->80)
t390 = cos(pkin(6));
t429 = (Ifges(3,6) - Ifges(4,6)) * t390;
t392 = sin(qJ(1));
t394 = cos(qJ(1));
t373 = -t394 * g(1) - t392 * g(2);
t395 = qJD(1) ^ 2;
t428 = -t395 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t373;
t427 = Ifges(4,4) + Ifges(3,5);
t389 = sin(pkin(6));
t387 = t389 ^ 2;
t388 = t390 ^ 2;
t418 = t388 * t395;
t425 = t387 * t395 + t418;
t354 = -t390 * g(3) - t428 * t389;
t424 = Ifges(3,1) + Ifges(4,1);
t423 = Ifges(3,4) - Ifges(4,5);
t422 = Ifges(3,2) + Ifges(4,3);
t421 = mrSges(3,2) * t389;
t420 = qJ(3) * t389;
t417 = t395 * qJ(2);
t368 = (-mrSges(4,1) * t390 - mrSges(4,3) * t389) * qJD(1);
t369 = (-mrSges(3,1) * t390 + t421) * qJD(1);
t403 = -pkin(2) * t390 - t420;
t367 = t403 * qJD(1);
t413 = t389 * qJD(1);
t340 = t367 * t413 + qJDD(3) - t354;
t336 = (-pkin(3) * t390 * t395 - pkin(5) * qJDD(1)) * t389 + t340;
t355 = -t389 * g(3) + t428 * t390;
t412 = t390 * qJD(1);
t342 = t367 * t412 + t355;
t409 = qJDD(1) * t390;
t337 = -pkin(3) * t418 - pkin(5) * t409 + t342;
t391 = sin(qJ(4));
t393 = cos(qJ(4));
t334 = t393 * t336 - t391 * t337;
t400 = -t389 * t391 - t390 * t393;
t364 = t400 * qJD(1);
t401 = t389 * t393 - t390 * t391;
t365 = t401 * qJD(1);
t348 = -t364 * mrSges(5,1) + t365 * mrSges(5,2);
t353 = t364 * qJD(4) + qJDD(1) * t401;
t356 = -qJD(4) * mrSges(5,2) + t364 * mrSges(5,3);
t332 = m(5) * t334 + qJDD(4) * mrSges(5,1) - t353 * mrSges(5,3) + qJD(4) * t356 - t365 * t348;
t335 = t391 * t336 + t393 * t337;
t352 = -t365 * qJD(4) + qJDD(1) * t400;
t357 = qJD(4) * mrSges(5,1) - t365 * mrSges(5,3);
t333 = m(5) * t335 - qJDD(4) * mrSges(5,2) + t352 * mrSges(5,3) - qJD(4) * t357 + t364 * t348;
t325 = t393 * t332 + t391 * t333;
t397 = -m(4) * t340 - t325;
t323 = m(3) * t354 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t368 - t369) * qJD(1)) * t389 + t397;
t404 = -t391 * t332 + t393 * t333;
t399 = m(4) * t342 + mrSges(4,2) * t409 + t368 * t412 + t404;
t324 = m(3) * t355 + (qJDD(1) * mrSges(3,3) + qJD(1) * t369) * t390 + t399;
t405 = -t389 * t323 + t390 * t324;
t318 = m(2) * t373 - t395 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t405;
t372 = t392 * g(1) - t394 * g(2);
t407 = -qJDD(2) + t372;
t398 = -0.2e1 * qJD(3) * t413 - t407;
t351 = -t417 + (-pkin(1) + t403) * qJDD(1) + t398;
t339 = (qJ(2) + (-t387 - t388) * pkin(5)) * t395 + (t420 + pkin(1) + (pkin(2) + pkin(3)) * t390) * qJDD(1) - t398;
t402 = -m(5) * t339 + t352 * mrSges(5,1) - t353 * mrSges(5,2) + t364 * t356 - t365 * t357;
t410 = qJDD(1) * t389;
t330 = m(4) * t351 - mrSges(4,1) * t409 - t425 * mrSges(4,2) - mrSges(4,3) * t410 + t402;
t363 = -qJDD(1) * pkin(1) - t407 - t417;
t396 = -m(3) * t363 + mrSges(3,1) * t409 + t425 * mrSges(3,3) - t330;
t329 = (mrSges(2,1) - t421) * qJDD(1) + t396 - t395 * mrSges(2,2) + m(2) * t372;
t416 = t392 * t318 + t394 * t329;
t319 = t390 * t323 + t389 * t324;
t414 = (t427 * t389 + t429) * qJD(1);
t406 = t394 * t318 - t392 * t329;
t345 = Ifges(5,1) * t365 + Ifges(5,4) * t364 + Ifges(5,5) * qJD(4);
t344 = Ifges(5,4) * t365 + Ifges(5,2) * t364 + Ifges(5,6) * qJD(4);
t343 = Ifges(5,5) * t365 + Ifges(5,6) * t364 + Ifges(5,3) * qJD(4);
t327 = mrSges(5,2) * t339 - mrSges(5,3) * t334 + Ifges(5,1) * t353 + Ifges(5,4) * t352 + Ifges(5,5) * qJDD(4) - qJD(4) * t344 + t364 * t343;
t326 = -mrSges(5,1) * t339 + mrSges(5,3) * t335 + Ifges(5,4) * t353 + Ifges(5,2) * t352 + Ifges(5,6) * qJDD(4) + qJD(4) * t345 - t365 * t343;
t315 = mrSges(3,2) * t363 + mrSges(4,2) * t340 - mrSges(3,3) * t354 - mrSges(4,3) * t351 - pkin(5) * t325 - qJ(3) * t330 - t391 * t326 + t393 * t327 + t414 * t412 + (t424 * t389 + t423 * t390) * qJDD(1);
t314 = -mrSges(3,1) * t363 + mrSges(3,3) * t355 - mrSges(4,1) * t351 + mrSges(4,2) * t342 - t391 * t327 - t393 * t326 - pkin(3) * t402 - pkin(5) * t404 - pkin(2) * t330 - t414 * t413 + (t423 * t389 + t422 * t390) * qJDD(1);
t313 = Ifges(5,3) * qJDD(4) + mrSges(2,1) * g(3) - pkin(2) * t397 + t395 * Ifges(2,5) - qJ(3) * t399 + mrSges(2,3) * t373 - t364 * t345 + t365 * t344 + Ifges(5,6) * t352 + Ifges(5,5) * t353 - mrSges(3,1) * t354 + mrSges(3,2) * t355 + mrSges(4,1) * t340 - mrSges(4,3) * t342 + mrSges(5,1) * t334 - mrSges(5,2) * t335 + pkin(3) * t325 - pkin(1) * t319 + (Ifges(2,6) - t429 + (mrSges(4,2) * pkin(2) - t427) * t389) * qJDD(1) + (t423 * t388 * qJD(1) + (pkin(2) * t368 - t423 * t413 + (-t422 + t424) * t412) * t389) * qJD(1);
t312 = -mrSges(2,2) * g(3) - mrSges(2,3) * t372 + Ifges(2,5) * qJDD(1) - t395 * Ifges(2,6) - qJ(2) * t319 - t389 * t314 + t390 * t315;
t1 = [-m(1) * g(1) + t406; -m(1) * g(2) + t416; (-m(1) - m(2)) * g(3) + t319; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t416 + t394 * t312 - t392 * t313; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t406 + t392 * t312 + t394 * t313; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t372 - mrSges(2,2) * t373 + t389 * t315 + t390 * t314 + pkin(1) * (-mrSges(3,2) * t410 + t396) + qJ(2) * t405;];
tauB = t1;
