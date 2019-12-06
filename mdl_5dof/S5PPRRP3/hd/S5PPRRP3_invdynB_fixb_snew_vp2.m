% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:40
% EndTime: 2019-12-05 15:10:42
% DurationCPUTime: 1.09s
% Computational Cost: add. (8519->182), mult. (14610->228), div. (0->0), fcn. (8262->8), ass. (0->78)
t388 = Ifges(5,1) + Ifges(6,1);
t384 = Ifges(5,4) - Ifges(6,5);
t383 = Ifges(6,4) + Ifges(5,5);
t387 = Ifges(5,2) + Ifges(6,3);
t382 = Ifges(5,6) - Ifges(6,6);
t386 = Ifges(5,3) + Ifges(6,2);
t385 = mrSges(5,3) + mrSges(6,2);
t357 = sin(pkin(7));
t359 = cos(pkin(7));
t346 = -g(1) * t359 - g(2) * t357;
t355 = -g(3) + qJDD(1);
t356 = sin(pkin(8));
t358 = cos(pkin(8));
t322 = t346 * t356 - t358 * t355;
t362 = cos(qJ(4));
t381 = t322 * t362;
t323 = t346 * t358 + t355 * t356;
t345 = g(1) * t357 - g(2) * t359;
t344 = qJDD(2) - t345;
t361 = sin(qJ(3));
t363 = cos(qJ(3));
t318 = t363 * t323 + t361 * t344;
t365 = qJD(3) ^ 2;
t316 = -pkin(3) * t365 + qJDD(3) * pkin(6) + t318;
t360 = sin(qJ(4));
t313 = t362 * t316 + t360 * t322;
t339 = (-mrSges(5,1) * t362 + mrSges(5,2) * t360) * qJD(3);
t374 = qJD(3) * qJD(4);
t341 = qJDD(3) * t362 - t360 * t374;
t376 = qJD(3) * t360;
t347 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t376;
t337 = (-pkin(4) * t362 - qJ(5) * t360) * qJD(3);
t364 = qJD(4) ^ 2;
t375 = qJD(3) * t362;
t309 = -pkin(4) * t364 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t337 * t375 + t313;
t338 = (-mrSges(6,1) * t362 - mrSges(6,3) * t360) * qJD(3);
t348 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t376;
t369 = m(6) * t309 + qJDD(4) * mrSges(6,3) + qJD(4) * t348 + t338 * t375;
t305 = m(5) * t313 - qJDD(4) * mrSges(5,2) - qJD(4) * t347 + t339 * t375 + t385 * t341 + t369;
t312 = -t316 * t360 + t381;
t340 = qJDD(3) * t360 + t362 * t374;
t349 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t375;
t310 = -qJDD(4) * pkin(4) - qJ(5) * t364 - t381 + qJDD(5) + (qJD(3) * t337 + t316) * t360;
t350 = mrSges(6,2) * t375 + qJD(4) * mrSges(6,3);
t368 = -m(6) * t310 + qJDD(4) * mrSges(6,1) + qJD(4) * t350;
t306 = m(5) * t312 + qJDD(4) * mrSges(5,1) + qJD(4) * t349 - t385 * t340 + (-t338 - t339) * t376 + t368;
t370 = t362 * t305 - t306 * t360;
t300 = m(4) * t318 - mrSges(4,1) * t365 - qJDD(3) * mrSges(4,2) + t370;
t317 = -t361 * t323 + t344 * t363;
t315 = -qJDD(3) * pkin(3) - pkin(6) * t365 - t317;
t311 = -pkin(4) * t341 - qJ(5) * t340 + (-0.2e1 * qJD(5) * t360 + (pkin(4) * t360 - qJ(5) * t362) * qJD(4)) * qJD(3) + t315;
t307 = m(6) * t311 - mrSges(6,1) * t341 - t340 * mrSges(6,3) - t348 * t376 - t350 * t375;
t366 = -m(5) * t315 + t341 * mrSges(5,1) - mrSges(5,2) * t340 - t347 * t376 + t349 * t375 - t307;
t303 = m(4) * t317 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t365 + t366;
t371 = t363 * t300 - t361 * t303;
t294 = m(3) * t323 + t371;
t302 = t360 * t305 + t306 * t362;
t301 = (-m(3) - m(4)) * t322 - t302;
t372 = t358 * t294 - t301 * t356;
t287 = m(2) * t346 + t372;
t295 = t300 * t361 + t303 * t363;
t367 = -m(3) * t344 - t295;
t293 = m(2) * t345 + t367;
t380 = t357 * t287 + t359 * t293;
t288 = t356 * t294 + t358 * t301;
t379 = -t382 * qJD(4) + (-t384 * t360 - t387 * t362) * qJD(3);
t378 = t386 * qJD(4) + (t383 * t360 + t382 * t362) * qJD(3);
t377 = t383 * qJD(4) + (t388 * t360 + t384 * t362) * qJD(3);
t373 = t359 * t287 - t293 * t357;
t297 = mrSges(5,2) * t315 + mrSges(6,2) * t310 - mrSges(5,3) * t312 - mrSges(6,3) * t311 - qJ(5) * t307 + t379 * qJD(4) + t383 * qJDD(4) + t388 * t340 + t384 * t341 + t378 * t375;
t296 = -mrSges(5,1) * t315 - mrSges(6,1) * t311 + mrSges(6,2) * t309 + mrSges(5,3) * t313 - pkin(4) * t307 + t377 * qJD(4) + t382 * qJDD(4) + t384 * t340 + t387 * t341 - t378 * t376;
t289 = Ifges(4,6) * qJDD(3) + t365 * Ifges(4,5) - mrSges(4,1) * t322 + mrSges(4,3) * t318 - mrSges(5,1) * t312 + mrSges(5,2) * t313 + mrSges(6,1) * t310 - mrSges(6,3) * t309 - pkin(4) * t368 - qJ(5) * t369 - pkin(3) * t302 + (-qJ(5) * mrSges(6,2) - t382) * t341 + (pkin(4) * mrSges(6,2) - t383) * t340 - t386 * qJDD(4) + (t377 * t362 + (pkin(4) * t338 + t379) * t360) * qJD(3);
t284 = mrSges(4,2) * t322 - mrSges(4,3) * t317 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t365 - pkin(6) * t302 - t296 * t360 + t297 * t362;
t283 = -mrSges(3,1) * t344 - mrSges(4,1) * t317 + mrSges(4,2) * t318 + mrSges(3,3) * t323 - Ifges(4,3) * qJDD(3) - pkin(2) * t295 - pkin(3) * t366 - pkin(6) * t370 - t362 * t296 - t360 * t297;
t282 = mrSges(3,2) * t344 + mrSges(3,3) * t322 - pkin(5) * t295 + t284 * t363 - t289 * t361;
t281 = -mrSges(2,1) * t355 + mrSges(2,3) * t346 + mrSges(3,1) * t322 + mrSges(3,2) * t323 - t361 * t284 - t363 * t289 - pkin(2) * (-m(4) * t322 - t302) - pkin(5) * t371 - pkin(1) * t288;
t280 = mrSges(2,2) * t355 - mrSges(2,3) * t345 - qJ(2) * t288 + t282 * t358 - t283 * t356;
t1 = [-m(1) * g(1) + t373; -m(1) * g(2) + t380; -m(1) * g(3) + m(2) * t355 + t288; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t380 + t359 * t280 - t357 * t281; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t373 + t357 * t280 + t359 * t281; -mrSges(1,1) * g(2) + mrSges(2,1) * t345 + mrSges(1,2) * g(1) - mrSges(2,2) * t346 + pkin(1) * t367 + qJ(2) * t372 + t356 * t282 + t358 * t283;];
tauB = t1;
