% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR2
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
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:11
% EndTime: 2019-12-31 16:48:12
% DurationCPUTime: 0.82s
% Computational Cost: add. (9703->153), mult. (13625->196), div. (0->0), fcn. (6934->8), ass. (0->69)
t389 = sin(qJ(1));
t392 = cos(qJ(1));
t373 = t389 * g(1) - g(2) * t392;
t368 = qJDD(1) * pkin(1) + t373;
t374 = -g(1) * t392 - g(2) * t389;
t393 = qJD(1) ^ 2;
t369 = -pkin(1) * t393 + t374;
t385 = sin(pkin(7));
t386 = cos(pkin(7));
t355 = t386 * t368 - t369 * t385;
t352 = qJDD(1) * pkin(2) + t355;
t356 = t385 * t368 + t386 * t369;
t353 = -pkin(2) * t393 + t356;
t388 = sin(qJ(3));
t391 = cos(qJ(3));
t349 = t388 * t352 + t391 * t353;
t381 = qJD(1) + qJD(3);
t379 = t381 ^ 2;
t380 = qJDD(1) + qJDD(3);
t346 = -pkin(3) * t379 + pkin(6) * t380 + t349;
t384 = -g(3) + qJDD(2);
t387 = sin(qJ(4));
t390 = cos(qJ(4));
t343 = -t346 * t387 + t384 * t390;
t344 = t346 * t390 + t384 * t387;
t358 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t387 + Ifges(5,2) * t390) * t381;
t359 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t387 + Ifges(5,4) * t390) * t381;
t404 = qJD(4) * t381;
t363 = t380 * t387 + t390 * t404;
t364 = t380 * t390 - t387 * t404;
t408 = mrSges(5,1) * t343 - mrSges(5,2) * t344 + Ifges(5,5) * t363 + Ifges(5,6) * t364 + Ifges(5,3) * qJDD(4) + (t358 * t387 - t359 * t390) * t381;
t407 = t381 * t387;
t406 = t381 * t390;
t362 = (-mrSges(5,1) * t390 + mrSges(5,2) * t387) * t381;
t371 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t406;
t341 = m(5) * t343 + qJDD(4) * mrSges(5,1) - t363 * mrSges(5,3) + qJD(4) * t371 - t362 * t407;
t370 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t407;
t342 = m(5) * t344 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t364 - qJD(4) * t370 + t362 * t406;
t399 = -t341 * t387 + t390 * t342;
t327 = m(4) * t349 - mrSges(4,1) * t379 - mrSges(4,2) * t380 + t399;
t348 = t352 * t391 - t353 * t388;
t345 = -pkin(3) * t380 - pkin(6) * t379 - t348;
t396 = -m(5) * t345 + t364 * mrSges(5,1) - t363 * mrSges(5,2) - t370 * t407 + t371 * t406;
t336 = m(4) * t348 + mrSges(4,1) * t380 - mrSges(4,2) * t379 + t396;
t324 = t327 * t388 + t391 * t336;
t320 = m(3) * t355 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t393 + t324;
t400 = t327 * t391 - t336 * t388;
t321 = m(3) * t356 - mrSges(3,1) * t393 - qJDD(1) * mrSges(3,2) + t400;
t315 = t320 * t386 + t321 * t385;
t312 = m(2) * t373 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t393 + t315;
t401 = -t320 * t385 + t321 * t386;
t313 = m(2) * t374 - mrSges(2,1) * t393 - qJDD(1) * mrSges(2,2) + t401;
t405 = t312 * t392 + t313 * t389;
t330 = t390 * t341 + t387 * t342;
t403 = m(4) * t384 + t330;
t402 = -t312 * t389 + t313 * t392;
t328 = m(3) * t384 + t403;
t357 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t387 + Ifges(5,6) * t390) * t381;
t333 = -mrSges(5,1) * t345 + mrSges(5,3) * t344 + Ifges(5,4) * t363 + Ifges(5,2) * t364 + Ifges(5,6) * qJDD(4) + qJD(4) * t359 - t357 * t407;
t334 = mrSges(5,2) * t345 - mrSges(5,3) * t343 + Ifges(5,1) * t363 + Ifges(5,4) * t364 + Ifges(5,5) * qJDD(4) - qJD(4) * t358 + t357 * t406;
t397 = mrSges(4,1) * t348 - mrSges(4,2) * t349 + Ifges(4,3) * t380 + pkin(3) * t396 + pkin(6) * t399 + t390 * t333 + t387 * t334;
t394 = mrSges(2,1) * t373 + mrSges(3,1) * t355 - mrSges(2,2) * t374 - mrSges(3,2) * t356 + pkin(1) * t315 + pkin(2) * t324 + t397 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t322 = -mrSges(4,1) * t384 + mrSges(4,3) * t349 + t379 * Ifges(4,5) + Ifges(4,6) * t380 - pkin(3) * t330 - t408;
t316 = mrSges(4,2) * t384 - mrSges(4,3) * t348 + Ifges(4,5) * t380 - Ifges(4,6) * t379 - pkin(6) * t330 - t333 * t387 + t334 * t390;
t308 = mrSges(3,2) * t384 - mrSges(3,3) * t355 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t393 - pkin(5) * t324 + t316 * t391 - t322 * t388;
t307 = -mrSges(3,1) * t384 + mrSges(3,3) * t356 + t393 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t403 + pkin(5) * t400 + t388 * t316 + t391 * t322;
t306 = -mrSges(2,2) * g(3) - mrSges(2,3) * t373 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t393 - qJ(2) * t315 - t307 * t385 + t308 * t386;
t305 = mrSges(2,1) * g(3) + mrSges(2,3) * t374 + t393 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t328 + qJ(2) * t401 + t386 * t307 + t385 * t308;
t1 = [-m(1) * g(1) + t402; -m(1) * g(2) + t405; (-m(1) - m(2)) * g(3) + t328; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t405 - t305 * t389 + t306 * t392; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t402 + t392 * t305 + t389 * t306; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t394; t394; t328; t397; t408;];
tauJB = t1;
