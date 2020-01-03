% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR3
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:32
% EndTime: 2019-12-31 17:01:33
% DurationCPUTime: 0.81s
% Computational Cost: add. (10586->154), mult. (13837->196), div. (0->0), fcn. (7044->8), ass. (0->68)
t387 = sin(qJ(1));
t390 = cos(qJ(1));
t372 = t387 * g(1) - t390 * g(2);
t367 = qJDD(1) * pkin(1) + t372;
t373 = -t390 * g(1) - t387 * g(2);
t391 = qJD(1) ^ 2;
t368 = -t391 * pkin(1) + t373;
t386 = sin(qJ(2));
t389 = cos(qJ(2));
t354 = t389 * t367 - t386 * t368;
t379 = qJDD(1) + qJDD(2);
t351 = t379 * pkin(2) + t354;
t355 = t386 * t367 + t389 * t368;
t380 = qJD(1) + qJD(2);
t378 = t380 ^ 2;
t352 = -t378 * pkin(2) + t355;
t383 = sin(pkin(7));
t384 = cos(pkin(7));
t348 = t383 * t351 + t384 * t352;
t345 = -t378 * pkin(3) + t379 * pkin(6) + t348;
t382 = -g(3) + qJDD(3);
t385 = sin(qJ(4));
t388 = cos(qJ(4));
t342 = -t385 * t345 + t388 * t382;
t343 = t388 * t345 + t385 * t382;
t357 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t385 + Ifges(5,2) * t388) * t380;
t358 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t385 + Ifges(5,4) * t388) * t380;
t401 = qJD(4) * t380;
t362 = t385 * t379 + t388 * t401;
t363 = t388 * t379 - t385 * t401;
t405 = mrSges(5,1) * t342 - mrSges(5,2) * t343 + Ifges(5,5) * t362 + Ifges(5,6) * t363 + Ifges(5,3) * qJDD(4) + (t357 * t385 - t358 * t388) * t380;
t404 = t380 * t385;
t403 = t380 * t388;
t361 = (-mrSges(5,1) * t388 + mrSges(5,2) * t385) * t380;
t370 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t403;
t340 = m(5) * t342 + qJDD(4) * mrSges(5,1) - t362 * mrSges(5,3) + qJD(4) * t370 - t361 * t404;
t369 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t404;
t341 = m(5) * t343 - qJDD(4) * mrSges(5,2) + t363 * mrSges(5,3) - qJD(4) * t369 + t361 * t403;
t397 = -t385 * t340 + t388 * t341;
t326 = m(4) * t348 - t378 * mrSges(4,1) - t379 * mrSges(4,2) + t397;
t347 = t384 * t351 - t383 * t352;
t344 = -t379 * pkin(3) - t378 * pkin(6) - t347;
t395 = -m(5) * t344 + t363 * mrSges(5,1) - t362 * mrSges(5,2) - t369 * t404 + t370 * t403;
t335 = m(4) * t347 + t379 * mrSges(4,1) - t378 * mrSges(4,2) + t395;
t323 = t383 * t326 + t384 * t335;
t319 = m(3) * t354 + t379 * mrSges(3,1) - t378 * mrSges(3,2) + t323;
t398 = t384 * t326 - t383 * t335;
t320 = m(3) * t355 - t378 * mrSges(3,1) - t379 * mrSges(3,2) + t398;
t314 = t389 * t319 + t386 * t320;
t311 = m(2) * t372 + qJDD(1) * mrSges(2,1) - t391 * mrSges(2,2) + t314;
t399 = -t386 * t319 + t389 * t320;
t312 = m(2) * t373 - t391 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t399;
t402 = t390 * t311 + t387 * t312;
t329 = t388 * t340 + t385 * t341;
t327 = m(4) * t382 + t329;
t400 = -t387 * t311 + t390 * t312;
t356 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t385 + Ifges(5,6) * t388) * t380;
t332 = -mrSges(5,1) * t344 + mrSges(5,3) * t343 + Ifges(5,4) * t362 + Ifges(5,2) * t363 + Ifges(5,6) * qJDD(4) + qJD(4) * t358 - t356 * t404;
t333 = mrSges(5,2) * t344 - mrSges(5,3) * t342 + Ifges(5,1) * t362 + Ifges(5,4) * t363 + Ifges(5,5) * qJDD(4) - qJD(4) * t357 + t356 * t403;
t393 = mrSges(3,1) * t354 + mrSges(4,1) * t347 - mrSges(3,2) * t355 - mrSges(4,2) * t348 + pkin(2) * t323 + pkin(3) * t395 + pkin(6) * t397 + t388 * t332 + t385 * t333 + (Ifges(4,3) + Ifges(3,3)) * t379;
t392 = mrSges(2,1) * t372 - mrSges(2,2) * t373 + Ifges(2,3) * qJDD(1) + pkin(1) * t314 + t393;
t321 = -mrSges(4,1) * t382 + mrSges(4,3) * t348 + t378 * Ifges(4,5) + Ifges(4,6) * t379 - pkin(3) * t329 - t405;
t315 = mrSges(4,2) * t382 - mrSges(4,3) * t347 + Ifges(4,5) * t379 - t378 * Ifges(4,6) - pkin(6) * t329 - t385 * t332 + t388 * t333;
t307 = -mrSges(3,2) * g(3) - mrSges(3,3) * t354 + Ifges(3,5) * t379 - t378 * Ifges(3,6) - qJ(3) * t323 + t384 * t315 - t383 * t321;
t306 = mrSges(3,1) * g(3) + mrSges(3,3) * t355 + t378 * Ifges(3,5) + Ifges(3,6) * t379 - pkin(2) * t327 + qJ(3) * t398 + t383 * t315 + t384 * t321;
t305 = -mrSges(2,2) * g(3) - mrSges(2,3) * t372 + Ifges(2,5) * qJDD(1) - t391 * Ifges(2,6) - pkin(5) * t314 - t386 * t306 + t389 * t307;
t304 = Ifges(2,6) * qJDD(1) + t391 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t373 + t386 * t307 + t389 * t306 - pkin(1) * (-m(3) * g(3) + t327) + pkin(5) * t399;
t1 = [-m(1) * g(1) + t400; -m(1) * g(2) + t402; (-m(1) - m(2) - m(3)) * g(3) + t327; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t402 - t387 * t304 + t390 * t305; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t400 + t390 * t304 + t387 * t305; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t392; t392; t393; t327; t405;];
tauJB = t1;
