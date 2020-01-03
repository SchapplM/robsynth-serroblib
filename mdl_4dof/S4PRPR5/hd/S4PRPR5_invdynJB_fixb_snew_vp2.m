% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRPR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:03
% EndTime: 2019-12-31 16:23:04
% DurationCPUTime: 0.63s
% Computational Cost: add. (5341->142), mult. (8720->185), div. (0->0), fcn. (4942->8), ass. (0->64)
t343 = sin(pkin(6));
t345 = cos(pkin(6));
t334 = -t345 * g(1) - t343 * g(2);
t341 = -g(3) + qJDD(1);
t347 = sin(qJ(2));
t349 = cos(qJ(2));
t321 = -t347 * t334 + t349 * t341;
t319 = qJDD(2) * pkin(2) + t321;
t322 = t349 * t334 + t347 * t341;
t350 = qJD(2) ^ 2;
t320 = -t350 * pkin(2) + t322;
t342 = sin(pkin(7));
t344 = cos(pkin(7));
t316 = t342 * t319 + t344 * t320;
t314 = -t350 * pkin(3) + qJDD(2) * pkin(5) + t316;
t333 = t343 * g(1) - t345 * g(2);
t332 = qJDD(3) - t333;
t346 = sin(qJ(4));
t348 = cos(qJ(4));
t311 = -t346 * t314 + t348 * t332;
t312 = t348 * t314 + t346 * t332;
t324 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t346 + Ifges(5,2) * t348) * qJD(2);
t325 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t346 + Ifges(5,4) * t348) * qJD(2);
t358 = qJD(2) * qJD(4);
t330 = t346 * qJDD(2) + t348 * t358;
t331 = t348 * qJDD(2) - t346 * t358;
t364 = mrSges(5,1) * t311 - mrSges(5,2) * t312 + Ifges(5,5) * t330 + Ifges(5,6) * t331 + Ifges(5,3) * qJDD(4) + (t324 * t346 - t325 * t348) * qJD(2);
t329 = (-mrSges(5,1) * t348 + mrSges(5,2) * t346) * qJD(2);
t359 = qJD(2) * t348;
t336 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t359;
t360 = qJD(2) * t346;
t308 = m(5) * t311 + qJDD(4) * mrSges(5,1) - t330 * mrSges(5,3) + qJD(4) * t336 - t329 * t360;
t335 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t360;
t309 = m(5) * t312 - qJDD(4) * mrSges(5,2) + t331 * mrSges(5,3) - qJD(4) * t335 + t329 * t359;
t300 = -t346 * t308 + t348 * t309;
t295 = m(4) * t316 - t350 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t300;
t315 = t344 * t319 - t342 * t320;
t313 = -qJDD(2) * pkin(3) - t350 * pkin(5) - t315;
t310 = -m(5) * t313 + t331 * mrSges(5,1) - t330 * mrSges(5,2) - t335 * t360 + t336 * t359;
t304 = m(4) * t315 + qJDD(2) * mrSges(4,1) - t350 * mrSges(4,2) + t310;
t292 = t342 * t295 + t344 * t304;
t323 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t346 + Ifges(5,6) * t348) * qJD(2);
t301 = -mrSges(5,1) * t313 + mrSges(5,3) * t312 + Ifges(5,4) * t330 + Ifges(5,2) * t331 + Ifges(5,6) * qJDD(4) + qJD(4) * t325 - t323 * t360;
t302 = mrSges(5,2) * t313 - mrSges(5,3) * t311 + Ifges(5,1) * t330 + Ifges(5,4) * t331 + Ifges(5,5) * qJDD(4) - qJD(4) * t324 + t323 * t359;
t363 = mrSges(3,1) * t321 + mrSges(4,1) * t315 - mrSges(3,2) * t322 - mrSges(4,2) * t316 + pkin(2) * t292 + pkin(3) * t310 + pkin(5) * t300 + t348 * t301 + t346 * t302 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2);
t290 = m(3) * t321 + qJDD(2) * mrSges(3,1) - t350 * mrSges(3,2) + t292;
t354 = t344 * t295 - t342 * t304;
t291 = m(3) * t322 - t350 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t354;
t355 = -t347 * t290 + t349 * t291;
t283 = m(2) * t334 + t355;
t299 = t348 * t308 + t346 * t309;
t298 = m(4) * t332 + t299;
t297 = (m(2) + m(3)) * t333 - t298;
t361 = t343 * t283 + t345 * t297;
t284 = t349 * t290 + t347 * t291;
t357 = m(2) * t341 + t284;
t356 = t345 * t283 - t343 * t297;
t286 = -mrSges(4,1) * t332 + mrSges(4,3) * t316 + t350 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t299 - t364;
t285 = mrSges(4,2) * t332 - mrSges(4,3) * t315 + Ifges(4,5) * qJDD(2) - t350 * Ifges(4,6) - pkin(5) * t299 - t346 * t301 + t348 * t302;
t280 = -mrSges(3,2) * t333 - mrSges(3,3) * t321 + Ifges(3,5) * qJDD(2) - t350 * Ifges(3,6) - qJ(3) * t292 + t344 * t285 - t342 * t286;
t279 = mrSges(3,1) * t333 + mrSges(3,3) * t322 + t350 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t298 + qJ(3) * t354 + t342 * t285 + t344 * t286;
t278 = -mrSges(2,1) * t341 + mrSges(2,3) * t334 - pkin(1) * t284 - t363;
t277 = mrSges(2,2) * t341 - mrSges(2,3) * t333 - pkin(4) * t284 - t347 * t279 + t349 * t280;
t1 = [-m(1) * g(1) + t356; -m(1) * g(2) + t361; -m(1) * g(3) + t357; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t361 + t345 * t277 - t343 * t278; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t356 + t343 * t277 + t345 * t278; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t333 - mrSges(2,2) * t334 + t347 * t280 + t349 * t279 + pkin(1) * (m(3) * t333 - t298) + pkin(4) * t355; t357; t363; t298; t364;];
tauJB = t1;
