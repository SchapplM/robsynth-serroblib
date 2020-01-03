% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:12
% EndTime: 2019-12-31 17:22:13
% DurationCPUTime: 0.84s
% Computational Cost: add. (12697->156), mult. (14090->198), div. (0->0), fcn. (7174->8), ass. (0->71)
t391 = sin(qJ(1));
t395 = cos(qJ(1));
t376 = t391 * g(1) - t395 * g(2);
t373 = qJDD(1) * pkin(1) + t376;
t377 = -t395 * g(1) - t391 * g(2);
t396 = qJD(1) ^ 2;
t374 = -t396 * pkin(1) + t377;
t390 = sin(qJ(2));
t394 = cos(qJ(2));
t358 = t394 * t373 - t390 * t374;
t385 = qJDD(1) + qJDD(2);
t355 = t385 * pkin(2) + t358;
t359 = t390 * t373 + t394 * t374;
t386 = qJD(1) + qJD(2);
t384 = t386 ^ 2;
t356 = -t384 * pkin(2) + t359;
t389 = sin(qJ(3));
t393 = cos(qJ(3));
t352 = t389 * t355 + t393 * t356;
t382 = qJD(3) + t386;
t380 = t382 ^ 2;
t381 = qJDD(3) + t385;
t349 = -t380 * pkin(3) + t381 * pkin(7) + t352;
t388 = sin(qJ(4));
t392 = cos(qJ(4));
t346 = -t392 * g(3) - t388 * t349;
t347 = -t388 * g(3) + t392 * t349;
t361 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t388 + Ifges(5,2) * t392) * t382;
t362 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t388 + Ifges(5,4) * t392) * t382;
t407 = qJD(4) * t382;
t366 = t388 * t381 + t392 * t407;
t367 = t392 * t381 - t388 * t407;
t412 = mrSges(5,1) * t346 - mrSges(5,2) * t347 + Ifges(5,5) * t366 + Ifges(5,6) * t367 + Ifges(5,3) * qJDD(4) + (t361 * t388 - t362 * t392) * t382;
t411 = -m(3) - m(4);
t410 = t382 * t388;
t409 = t382 * t392;
t365 = (-mrSges(5,1) * t392 + mrSges(5,2) * t388) * t382;
t372 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t409;
t344 = m(5) * t346 + qJDD(4) * mrSges(5,1) - t366 * mrSges(5,3) + qJD(4) * t372 - t365 * t410;
t371 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t410;
t345 = m(5) * t347 - qJDD(4) * mrSges(5,2) + t367 * mrSges(5,3) - qJD(4) * t371 + t365 * t409;
t403 = -t388 * t344 + t392 * t345;
t331 = m(4) * t352 - t380 * mrSges(4,1) - t381 * mrSges(4,2) + t403;
t351 = t393 * t355 - t389 * t356;
t348 = -t381 * pkin(3) - t380 * pkin(7) - t351;
t400 = -m(5) * t348 + t367 * mrSges(5,1) - t366 * mrSges(5,2) - t371 * t410 + t372 * t409;
t339 = m(4) * t351 + t381 * mrSges(4,1) - t380 * mrSges(4,2) + t400;
t328 = t389 * t331 + t393 * t339;
t324 = m(3) * t358 + t385 * mrSges(3,1) - t384 * mrSges(3,2) + t328;
t404 = t393 * t331 - t389 * t339;
t325 = m(3) * t359 - t384 * mrSges(3,1) - t385 * mrSges(3,2) + t404;
t319 = t394 * t324 + t390 * t325;
t316 = m(2) * t376 + qJDD(1) * mrSges(2,1) - t396 * mrSges(2,2) + t319;
t405 = -t390 * t324 + t394 * t325;
t317 = m(2) * t377 - t396 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t405;
t408 = t395 * t316 + t391 * t317;
t333 = t392 * t344 + t388 * t345;
t406 = -t391 * t316 + t395 * t317;
t360 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t388 + Ifges(5,6) * t392) * t382;
t336 = -mrSges(5,1) * t348 + mrSges(5,3) * t347 + Ifges(5,4) * t366 + Ifges(5,2) * t367 + Ifges(5,6) * qJDD(4) + qJD(4) * t362 - t360 * t410;
t337 = mrSges(5,2) * t348 - mrSges(5,3) * t346 + Ifges(5,1) * t366 + Ifges(5,4) * t367 + Ifges(5,5) * qJDD(4) - qJD(4) * t361 + t360 * t409;
t401 = mrSges(4,1) * t351 - mrSges(4,2) * t352 + Ifges(4,3) * t381 + pkin(3) * t400 + pkin(7) * t403 + t392 * t336 + t388 * t337;
t398 = mrSges(3,1) * t358 - mrSges(3,2) * t359 + Ifges(3,3) * t385 + pkin(2) * t328 + t401;
t397 = mrSges(2,1) * t376 - mrSges(2,2) * t377 + Ifges(2,3) * qJDD(1) + pkin(1) * t319 + t398;
t326 = mrSges(4,1) * g(3) + mrSges(4,3) * t352 + t380 * Ifges(4,5) + Ifges(4,6) * t381 - pkin(3) * t333 - t412;
t320 = -mrSges(4,2) * g(3) - mrSges(4,3) * t351 + Ifges(4,5) * t381 - t380 * Ifges(4,6) - pkin(7) * t333 - t388 * t336 + t392 * t337;
t312 = -mrSges(3,2) * g(3) - mrSges(3,3) * t358 + Ifges(3,5) * t385 - t384 * Ifges(3,6) - pkin(6) * t328 + t393 * t320 - t389 * t326;
t311 = Ifges(3,6) * t385 + t384 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t359 + t389 * t320 + t393 * t326 - pkin(2) * (-m(4) * g(3) + t333) + pkin(6) * t404;
t310 = -mrSges(2,2) * g(3) - mrSges(2,3) * t376 + Ifges(2,5) * qJDD(1) - t396 * Ifges(2,6) - pkin(5) * t319 - t390 * t311 + t394 * t312;
t309 = Ifges(2,6) * qJDD(1) + t396 * Ifges(2,5) + mrSges(2,3) * t377 + t390 * t312 + t394 * t311 - pkin(1) * t333 + pkin(5) * t405 + (-pkin(1) * t411 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t406; -m(1) * g(2) + t408; (-m(1) - m(2) + t411) * g(3) + t333; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t408 - t391 * t309 + t395 * t310; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t406 + t395 * t309 + t391 * t310; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t397; t397; t398; t401; t412;];
tauJB = t1;
