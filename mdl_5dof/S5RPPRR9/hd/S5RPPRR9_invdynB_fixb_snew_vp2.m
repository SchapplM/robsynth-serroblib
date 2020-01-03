% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:26
% EndTime: 2019-12-31 18:02:28
% DurationCPUTime: 1.20s
% Computational Cost: add. (13703->224), mult. (23932->270), div. (0->0), fcn. (10899->8), ass. (0->92)
t409 = -pkin(1) - pkin(2);
t408 = mrSges(2,1) + mrSges(3,1);
t407 = Ifges(3,4) + Ifges(2,5);
t406 = Ifges(2,6) - Ifges(3,6);
t375 = g(3) + qJDD(3);
t384 = cos(qJ(4));
t405 = t384 * t375;
t382 = sin(qJ(1));
t385 = cos(qJ(1));
t367 = -t385 * g(1) - t382 * g(2);
t387 = qJD(1) ^ 2;
t393 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t367;
t350 = -t387 * pkin(1) + t393;
t347 = t409 * t387 + t393;
t366 = t382 * g(1) - t385 * g(2);
t392 = -t387 * qJ(2) + qJDD(2) - t366;
t349 = t409 * qJDD(1) + t392;
t378 = sin(pkin(8));
t379 = cos(pkin(8));
t333 = t379 * t347 + t378 * t349;
t331 = -t387 * pkin(3) - qJDD(1) * pkin(6) + t333;
t381 = sin(qJ(4));
t328 = t384 * t331 + t381 * t375;
t360 = (mrSges(5,1) * t384 - mrSges(5,2) * t381) * qJD(1);
t401 = qJD(1) * qJD(4);
t400 = t381 * t401;
t363 = -t384 * qJDD(1) + t400;
t403 = qJD(1) * t381;
t364 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t403;
t332 = -t378 * t347 + t379 * t349;
t330 = qJDD(1) * pkin(3) - t387 * pkin(6) - t332;
t399 = t384 * t401;
t362 = -t381 * qJDD(1) - t399;
t324 = (-t362 + t399) * pkin(7) + (-t363 - t400) * pkin(4) + t330;
t361 = (pkin(4) * t384 + pkin(7) * t381) * qJD(1);
t386 = qJD(4) ^ 2;
t402 = t384 * qJD(1);
t326 = -t386 * pkin(4) + qJDD(4) * pkin(7) - t361 * t402 + t328;
t380 = sin(qJ(5));
t383 = cos(qJ(5));
t322 = t383 * t324 - t380 * t326;
t358 = t383 * qJD(4) + t380 * t403;
t340 = t358 * qJD(5) + t380 * qJDD(4) + t383 * t362;
t359 = t380 * qJD(4) - t383 * t403;
t341 = -t358 * mrSges(6,1) + t359 * mrSges(6,2);
t368 = qJD(5) + t402;
t345 = -t368 * mrSges(6,2) + t358 * mrSges(6,3);
t357 = qJDD(5) - t363;
t320 = m(6) * t322 + t357 * mrSges(6,1) - t340 * mrSges(6,3) - t359 * t341 + t368 * t345;
t323 = t380 * t324 + t383 * t326;
t339 = -t359 * qJD(5) + t383 * qJDD(4) - t380 * t362;
t346 = t368 * mrSges(6,1) - t359 * mrSges(6,3);
t321 = m(6) * t323 - t357 * mrSges(6,2) + t339 * mrSges(6,3) + t358 * t341 - t368 * t346;
t395 = -t380 * t320 + t383 * t321;
t314 = m(5) * t328 - qJDD(4) * mrSges(5,2) + t363 * mrSges(5,3) - qJD(4) * t364 - t360 * t402 + t395;
t327 = -t381 * t331 + t405;
t365 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t402;
t325 = -qJDD(4) * pkin(4) - t386 * pkin(7) - t405 + (-qJD(1) * t361 + t331) * t381;
t390 = -m(6) * t325 + t339 * mrSges(6,1) - t340 * mrSges(6,2) + t358 * t345 - t359 * t346;
t318 = m(5) * t327 + qJDD(4) * mrSges(5,1) - t362 * mrSges(5,3) + qJD(4) * t365 + t360 * t403 + t390;
t396 = t384 * t314 - t381 * t318;
t309 = m(4) * t333 - t387 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t396;
t315 = t383 * t320 + t380 * t321;
t388 = -m(5) * t330 + t363 * mrSges(5,1) - t362 * mrSges(5,2) + t364 * t403 - t365 * t402 - t315;
t312 = m(4) * t332 - qJDD(1) * mrSges(4,1) - t387 * mrSges(4,2) + t388;
t397 = t379 * t309 - t378 * t312;
t394 = m(3) * t350 + qJDD(1) * mrSges(3,3) + t397;
t303 = m(2) * t367 - qJDD(1) * mrSges(2,2) - t408 * t387 + t394;
t305 = t378 * t309 + t379 * t312;
t351 = -qJDD(1) * pkin(1) + t392;
t389 = -m(3) * t351 + qJDD(1) * mrSges(3,1) + t387 * mrSges(3,3) - t305;
t304 = m(2) * t366 + qJDD(1) * mrSges(2,1) - t387 * mrSges(2,2) + t389;
t404 = t382 * t303 + t385 * t304;
t398 = t385 * t303 - t382 * t304;
t311 = t381 * t314 + t384 * t318;
t391 = -m(4) * t375 - t311;
t354 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t381 - Ifges(5,4) * t384) * qJD(1);
t353 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t381 - Ifges(5,2) * t384) * qJD(1);
t352 = (Ifges(5,3) * qJD(4)) + (-Ifges(5,5) * t381 - Ifges(5,6) * t384) * qJD(1);
t336 = Ifges(6,1) * t359 + Ifges(6,4) * t358 + Ifges(6,5) * t368;
t335 = Ifges(6,4) * t359 + Ifges(6,2) * t358 + Ifges(6,6) * t368;
t334 = Ifges(6,5) * t359 + Ifges(6,6) * t358 + Ifges(6,3) * t368;
t317 = mrSges(6,2) * t325 - mrSges(6,3) * t322 + Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t357 + t358 * t334 - t368 * t335;
t316 = -mrSges(6,1) * t325 + mrSges(6,3) * t323 + Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t357 - t359 * t334 + t368 * t336;
t310 = -m(3) * g(3) + t391;
t307 = -mrSges(5,1) * t330 - mrSges(6,1) * t322 + mrSges(6,2) * t323 + mrSges(5,3) * t328 + Ifges(5,4) * t362 - Ifges(6,5) * t340 + Ifges(5,2) * t363 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t339 - Ifges(6,3) * t357 - pkin(4) * t315 + qJD(4) * t354 - t359 * t335 + t358 * t336 + t352 * t403;
t306 = mrSges(5,2) * t330 - mrSges(5,3) * t327 + Ifges(5,1) * t362 + Ifges(5,4) * t363 + Ifges(5,5) * qJDD(4) - pkin(7) * t315 - qJD(4) * t353 - t380 * t316 + t383 * t317 - t352 * t402;
t299 = -Ifges(4,6) * qJDD(1) + t387 * Ifges(4,5) - mrSges(4,1) * t375 + mrSges(4,3) * t333 - Ifges(5,5) * t362 - Ifges(5,6) * t363 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t327 + mrSges(5,2) * t328 - t380 * t317 - t383 * t316 - pkin(4) * t390 - pkin(7) * t395 - pkin(3) * t311 + (t381 * t353 - t384 * t354) * qJD(1);
t298 = mrSges(4,2) * t375 - mrSges(4,3) * t332 - Ifges(4,5) * qJDD(1) - t387 * Ifges(4,6) - pkin(6) * t311 + t384 * t306 - t381 * t307;
t297 = mrSges(3,2) * t351 - mrSges(2,3) * t366 - qJ(2) * t310 - qJ(3) * t305 + t379 * t298 - t378 * t299 - t406 * t387 + t407 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t296 = mrSges(3,2) * t350 + mrSges(2,3) * t367 - pkin(1) * t310 - pkin(2) * t391 + t408 * g(3) - qJ(3) * t397 + t406 * qJDD(1) - t378 * t298 - t379 * t299 + t407 * t387;
t1 = [-m(1) * g(1) + t398; -m(1) * g(2) + t404; (-m(1) - m(2) - m(3)) * g(3) + t391; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t404 - t382 * t296 + t385 * t297; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t398 + t385 * t296 + t382 * t297; pkin(1) * t389 + qJ(2) * (-t387 * mrSges(3,1) + t394) - mrSges(2,2) * t367 + mrSges(2,1) * t366 - pkin(2) * t305 - mrSges(3,1) * t351 + mrSges(3,3) * t350 - t381 * t306 - t384 * t307 - pkin(3) * t388 - pkin(6) * t396 - mrSges(4,1) * t332 + mrSges(4,2) * t333 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB = t1;
