% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:45
% EndTime: 2019-12-31 18:10:47
% DurationCPUTime: 1.06s
% Computational Cost: add. (1336->175), mult. (2586->208), div. (0->0), fcn. (1114->6), ass. (0->75)
t373 = sin(qJ(3));
t375 = cos(qJ(3));
t396 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t413 = t373 * (Ifges(4,1) + Ifges(5,1) + Ifges(6,1)) + t375 * t396;
t412 = -t375 * (Ifges(4,2) + Ifges(5,3) + Ifges(6,2)) - t373 * t396;
t395 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t394 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t407 = t373 * (t412 * qJD(1) - t394 * qJD(3));
t406 = 2 * qJD(4);
t405 = pkin(3) + pkin(4);
t378 = qJD(1) ^ 2;
t404 = pkin(4) * t378;
t403 = pkin(6) * t378;
t402 = mrSges(4,3) + mrSges(5,2);
t401 = qJ(4) * t375;
t370 = -g(3) + qJDD(2);
t400 = t370 * t375;
t374 = sin(qJ(1));
t376 = cos(qJ(1));
t389 = t374 * g(1) - g(2) * t376;
t341 = qJDD(1) * pkin(1) + t389;
t386 = -g(1) * t376 - g(2) * t374;
t346 = -pkin(1) * t378 + t386;
t371 = sin(pkin(7));
t372 = cos(pkin(7));
t312 = t371 * t341 + t372 * t346;
t310 = -pkin(2) * t378 + qJDD(1) * pkin(6) + t312;
t306 = t375 * t310 + t373 * t370;
t311 = t372 * t341 - t371 * t346;
t399 = qJD(1) * t373;
t398 = t375 * qJD(1);
t397 = qJD(1) * qJD(3);
t393 = -0.2e1 * qJD(1) * qJD(5);
t391 = t413 * qJD(1) + t395 * qJD(3);
t390 = pkin(1) * t372 + pkin(2);
t344 = (t375 * mrSges(6,1) + t373 * mrSges(6,2)) * qJD(1);
t345 = (-t375 * mrSges(4,1) + t373 * mrSges(4,2)) * qJD(1);
t348 = qJDD(1) * t375 - t373 * t397;
t355 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t399;
t342 = (-t375 * pkin(3) - t373 * qJ(4)) * qJD(1);
t377 = qJD(3) ^ 2;
t303 = -pkin(3) * t377 + qJDD(3) * qJ(4) + qJD(3) * t406 + t342 * t398 + t306;
t343 = (-t375 * mrSges(5,1) - t373 * mrSges(5,3)) * qJD(1);
t356 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t399;
t353 = -qJD(3) * pkin(4) - qJ(5) * t399;
t369 = t375 ^ 2;
t299 = -qJ(5) * t348 + qJD(3) * t353 - t369 * t404 + t375 * t393 + t303;
t354 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t399;
t387 = m(6) * t299 + qJDD(3) * mrSges(6,2) - t348 * mrSges(6,3) + qJD(3) * t354;
t382 = m(5) * t303 + qJDD(3) * mrSges(5,3) + qJD(3) * t356 + t343 * t398 + t387;
t290 = m(4) * t306 - qJDD(3) * mrSges(4,2) - qJD(3) * t355 + t402 * t348 + (-t344 + t345) * t398 + t382;
t307 = t373 * t310;
t305 = -t307 + t400;
t347 = qJDD(1) * t373 + t375 * t397;
t358 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t398;
t384 = -qJ(4) * t377 + t342 * t399 + qJDD(4) + t307;
t300 = t373 * t393 - qJ(5) * t347 - t405 * qJDD(3) + (qJ(5) * t397 - t373 * t404 - t370) * t375 + t384;
t357 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t398;
t295 = m(6) * t300 - qJDD(3) * mrSges(6,1) - t347 * mrSges(6,3) - qJD(3) * t357 - t344 * t399;
t304 = -qJDD(3) * pkin(3) + t384 - t400;
t359 = mrSges(5,2) * t398 + qJD(3) * mrSges(5,3);
t380 = -m(5) * t304 + qJDD(3) * mrSges(5,1) + qJD(3) * t359 - t295;
t291 = m(4) * t305 + qJDD(3) * mrSges(4,1) + qJD(3) * t358 - t402 * t347 + (-t343 - t345) * t399 + t380;
t388 = t375 * t290 - t373 * t291;
t385 = qJDD(1) * pkin(2) + t311;
t383 = -qJ(4) * t347 - t385;
t297 = qJDD(5) + (-qJ(5) * t369 + pkin(6)) * t378 + t405 * t348 + (qJD(3) * t401 + (-pkin(3) * qJD(3) + t353 + t406) * t373) * qJD(1) - t383;
t294 = m(6) * t297 + t348 * mrSges(6,1) + t347 * mrSges(6,2) + t354 * t399 + t357 * t398;
t301 = -pkin(3) * t348 - t403 + (-0.2e1 * qJD(4) * t373 + (pkin(3) * t373 - t401) * qJD(3)) * qJD(1) + t383;
t381 = m(5) * t301 - t347 * mrSges(5,3) - t356 * t399 - t359 * t398 - t294;
t309 = -t385 - t403;
t379 = -m(4) * t309 + t348 * mrSges(4,1) - t355 * t399 + t358 * t398 - t381;
t293 = mrSges(5,2) * t347 + t343 * t399 - t380;
t292 = -mrSges(5,1) * t348 + t381;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t389 - mrSges(2,2) * t386 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t311 - mrSges(3,2) * t312 + t373 * (mrSges(4,2) * t309 + mrSges(5,2) * t304 + mrSges(6,2) * t297 - mrSges(4,3) * t305 - mrSges(5,3) * t301 - mrSges(6,3) * t300 - qJ(4) * t292 - qJ(5) * t295) + t375 * (-mrSges(4,1) * t309 + mrSges(4,3) * t306 - mrSges(5,1) * t301 + mrSges(5,2) * t303 + mrSges(6,1) * t297 - mrSges(6,3) * t299 + pkin(4) * t294 - qJ(5) * (-t344 * t398 + t387) - pkin(3) * t292) + pkin(2) * t379 + pkin(6) * t388 + pkin(1) * (t371 * (m(3) * t312 - mrSges(3,1) * t378 - qJDD(1) * mrSges(3,2) + t388) + t372 * (m(3) * t311 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t378 + t379)) + (t390 * mrSges(5,1) - t412) * t348 + (-t390 * mrSges(4,2) + t413) * t347 + (t373 * t395 + t375 * t394) * qJDD(3) + (t375 * t391 + t407) * qJD(3); m(3) * t370 + t290 * t373 + t291 * t375; mrSges(4,1) * t305 - mrSges(4,2) * t306 - mrSges(5,1) * t304 + mrSges(5,3) * t303 - mrSges(6,1) * t300 + mrSges(6,2) * t299 - pkin(4) * t295 - pkin(3) * t293 + qJ(4) * t382 + (qJ(4) * mrSges(5,2) + t394) * t348 + t395 * t347 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * qJDD(3) + (-t407 + (-qJ(4) * t344 - t391) * t375) * qJD(1); t293; t294;];
tauJ = t1;
