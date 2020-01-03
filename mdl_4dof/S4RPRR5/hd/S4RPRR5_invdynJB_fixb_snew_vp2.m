% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:35
% EndTime: 2019-12-31 16:51:36
% DurationCPUTime: 0.55s
% Computational Cost: add. (4982->151), mult. (6517->182), div. (0->0), fcn. (2328->6), ass. (0->66)
t347 = qJD(1) ^ 2;
t343 = sin(qJ(1));
t346 = cos(qJ(1));
t329 = -t346 * g(1) - t343 * g(2);
t352 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t329;
t364 = -pkin(1) - pkin(2);
t311 = t364 * t347 + t352;
t328 = t343 * g(1) - t346 * g(2);
t351 = -t347 * qJ(2) + qJDD(2) - t328;
t314 = t364 * qJDD(1) + t351;
t342 = sin(qJ(3));
t345 = cos(qJ(3));
t308 = t345 * t311 + t342 * t314;
t334 = -qJD(1) + qJD(3);
t332 = t334 ^ 2;
t333 = -qJDD(1) + qJDD(3);
t305 = -(t332 * pkin(3)) + t333 * pkin(6) + t308;
t341 = sin(qJ(4));
t344 = cos(qJ(4));
t302 = t344 * g(3) - t341 * t305;
t303 = t341 * g(3) + t344 * t305;
t316 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t341 + Ifges(5,2) * t344) * t334;
t317 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t341 + Ifges(5,4) * t344) * t334;
t357 = qJD(4) * t334;
t323 = t341 * t333 + t344 * t357;
t324 = t344 * t333 - t341 * t357;
t366 = mrSges(5,1) * t302 - mrSges(5,2) * t303 + Ifges(5,5) * t323 + Ifges(5,6) * t324 + Ifges(5,3) * qJDD(4) + (t316 * t341 - t317 * t344) * t334;
t365 = -m(3) - m(4);
t363 = -mrSges(2,1) - mrSges(3,1);
t362 = Ifges(3,4) + Ifges(2,5);
t361 = Ifges(2,6) - Ifges(3,6);
t360 = t334 * t341;
t359 = t334 * t344;
t318 = -t347 * pkin(1) + t352;
t322 = (-mrSges(5,1) * t344 + mrSges(5,2) * t341) * t334;
t326 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t359;
t300 = m(5) * t302 + qJDD(4) * mrSges(5,1) - t323 * mrSges(5,3) + qJD(4) * t326 - t322 * t360;
t325 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t360;
t301 = m(5) * t303 - qJDD(4) * mrSges(5,2) + t324 * mrSges(5,3) - qJD(4) * t325 + t322 * t359;
t294 = -t341 * t300 + t344 * t301;
t291 = m(4) * t308 - (t332 * mrSges(4,1)) - t333 * mrSges(4,2) + t294;
t307 = -t342 * t311 + t345 * t314;
t304 = -t333 * pkin(3) - t332 * pkin(6) - t307;
t298 = -m(5) * t304 + t324 * mrSges(5,1) - t323 * mrSges(5,2) - t325 * t360 + t326 * t359;
t297 = m(4) * t307 + t333 * mrSges(4,1) - t332 * mrSges(4,2) + t298;
t355 = t345 * t291 - t342 * t297;
t353 = m(3) * t318 + qJDD(1) * mrSges(3,3) + t355;
t283 = m(2) * t329 - qJDD(1) * mrSges(2,2) + t363 * t347 + t353;
t289 = t342 * t291 + t345 * t297;
t321 = -qJDD(1) * pkin(1) + t351;
t288 = m(3) * t321 - qJDD(1) * mrSges(3,1) - t347 * mrSges(3,3) + t289;
t284 = m(2) * t328 + qJDD(1) * mrSges(2,1) - t347 * mrSges(2,2) - t288;
t358 = t343 * t283 + t346 * t284;
t356 = t346 * t283 - t343 * t284;
t293 = t344 * t300 + t341 * t301;
t315 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t341 + Ifges(5,6) * t344) * t334;
t295 = -mrSges(5,1) * t304 + mrSges(5,3) * t303 + Ifges(5,4) * t323 + Ifges(5,2) * t324 + Ifges(5,6) * qJDD(4) + qJD(4) * t317 - t315 * t360;
t296 = mrSges(5,2) * t304 - mrSges(5,3) * t302 + Ifges(5,1) * t323 + Ifges(5,4) * t324 + Ifges(5,5) * qJDD(4) - qJD(4) * t316 + t315 * t359;
t349 = mrSges(4,1) * t307 - mrSges(4,2) * t308 + Ifges(4,3) * t333 + pkin(3) * t298 + pkin(6) * t294 + t344 * t295 + t341 * t296;
t348 = -mrSges(3,1) * t321 - mrSges(2,2) * t329 - pkin(2) * t289 + qJ(2) * (-t347 * mrSges(3,1) + t353) - pkin(1) * t288 + mrSges(3,3) * t318 + mrSges(2,1) * t328 - t349 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);
t292 = t365 * g(3) - t293;
t285 = -mrSges(4,1) * g(3) + mrSges(4,3) * t308 + (t332 * Ifges(4,5)) + Ifges(4,6) * t333 - pkin(3) * t293 - t366;
t279 = mrSges(4,2) * g(3) - mrSges(4,3) * t307 + Ifges(4,5) * t333 - t332 * Ifges(4,6) - pkin(6) * t293 - t341 * t295 + t344 * t296;
t278 = mrSges(3,2) * t321 - mrSges(2,3) * t328 - pkin(5) * t289 - qJ(2) * t292 + t345 * t279 - t342 * t285 - t361 * t347 + t362 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t277 = mrSges(2,3) * t329 + mrSges(3,2) * t318 - t342 * t279 - t345 * t285 + pkin(2) * t293 - pkin(5) * t355 - pkin(1) * t292 + t362 * t347 + t361 * qJDD(1) + (pkin(2) * m(4) - t363) * g(3);
t1 = [-m(1) * g(1) + t356; -m(1) * g(2) + t358; (-m(1) - m(2) + t365) * g(3) - t293; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t358 - t343 * t277 + t346 * t278; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t356 + t346 * t277 + t343 * t278; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t348; t348; t288; t349; t366;];
tauJB = t1;
