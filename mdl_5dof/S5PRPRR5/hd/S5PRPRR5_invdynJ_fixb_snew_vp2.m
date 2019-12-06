% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:48
% EndTime: 2019-12-05 15:53:50
% DurationCPUTime: 0.92s
% Computational Cost: add. (5958->184), mult. (13583->239), div. (0->0), fcn. (9694->10), ass. (0->87)
t359 = qJD(2) ^ 2;
t351 = cos(pkin(9));
t346 = t351 ^ 2;
t379 = 0.2e1 * t351;
t378 = pkin(3) * t351;
t349 = sin(pkin(9));
t377 = mrSges(4,2) * t349;
t376 = t346 * t359;
t350 = sin(pkin(8));
t352 = cos(pkin(8));
t336 = -t352 * g(1) - t350 * g(2);
t348 = -g(3) + qJDD(1);
t355 = sin(qJ(2));
t358 = cos(qJ(2));
t325 = t358 * t336 + t355 * t348;
t320 = -t359 * pkin(2) + qJDD(2) * qJ(3) + t325;
t335 = -t350 * g(1) + t352 * g(2);
t372 = qJD(2) * qJD(3);
t374 = t351 * t335 - 0.2e1 * t349 * t372;
t301 = (-pkin(6) * qJDD(2) + t359 * t378 - t320) * t349 + t374;
t304 = t351 * t320 + t349 * t335 + t372 * t379;
t371 = qJDD(2) * t351;
t302 = -pkin(3) * t376 + pkin(6) * t371 + t304;
t354 = sin(qJ(4));
t357 = cos(qJ(4));
t283 = t357 * t301 - t354 * t302;
t365 = t349 * t357 + t351 * t354;
t364 = -t349 * t354 + t351 * t357;
t327 = t364 * qJD(2);
t373 = t327 * qJD(4);
t317 = t365 * qJDD(2) + t373;
t328 = t365 * qJD(2);
t279 = (-t317 + t373) * pkin(7) + (t327 * t328 + qJDD(4)) * pkin(4) + t283;
t284 = t354 * t301 + t357 * t302;
t316 = -t328 * qJD(4) + t364 * qJDD(2);
t323 = qJD(4) * pkin(4) - t328 * pkin(7);
t326 = t327 ^ 2;
t280 = -t326 * pkin(4) + t316 * pkin(7) - qJD(4) * t323 + t284;
t353 = sin(qJ(5));
t356 = cos(qJ(5));
t277 = t356 * t279 - t353 * t280;
t311 = t356 * t327 - t353 * t328;
t291 = t311 * qJD(5) + t353 * t316 + t356 * t317;
t312 = t353 * t327 + t356 * t328;
t297 = -t311 * mrSges(6,1) + t312 * mrSges(6,2);
t347 = qJD(4) + qJD(5);
t305 = -t347 * mrSges(6,2) + t311 * mrSges(6,3);
t344 = qJDD(4) + qJDD(5);
t274 = m(6) * t277 + t344 * mrSges(6,1) - t291 * mrSges(6,3) - t312 * t297 + t347 * t305;
t278 = t353 * t279 + t356 * t280;
t290 = -t312 * qJD(5) + t356 * t316 - t353 * t317;
t306 = t347 * mrSges(6,1) - t312 * mrSges(6,3);
t275 = m(6) * t278 - t344 * mrSges(6,2) + t290 * mrSges(6,3) + t311 * t297 - t347 * t306;
t267 = t356 * t274 + t353 * t275;
t314 = -t327 * mrSges(5,1) + t328 * mrSges(5,2);
t321 = -qJD(4) * mrSges(5,2) + t327 * mrSges(5,3);
t265 = m(5) * t283 + qJDD(4) * mrSges(5,1) - t317 * mrSges(5,3) + qJD(4) * t321 - t328 * t314 + t267;
t322 = qJD(4) * mrSges(5,1) - t328 * mrSges(5,3);
t368 = -t353 * t274 + t356 * t275;
t266 = m(5) * t284 - qJDD(4) * mrSges(5,2) + t316 * mrSges(5,3) - qJD(4) * t322 + t327 * t314 + t368;
t375 = t357 * t265 + t354 * t266;
t303 = -t349 * t320 + t374;
t363 = mrSges(4,3) * qJDD(2) + t359 * (-t351 * mrSges(4,1) + t377);
t369 = -t354 * t265 + t357 * t266;
t370 = -t349 * (m(4) * t303 - t363 * t349 + t375) + t351 * (m(4) * t304 + t363 * t351 + t369);
t324 = -t355 * t336 + t358 * t348;
t367 = qJDD(3) - t324;
t345 = t349 ^ 2;
t307 = (-pkin(2) - t378) * qJDD(2) + (-qJ(3) + (-t345 - t346) * pkin(6)) * t359 + t367;
t282 = -t316 * pkin(4) - t326 * pkin(7) + t328 * t323 + t307;
t366 = m(6) * t282 - t290 * mrSges(6,1) + t291 * mrSges(6,2) - t311 * t305 + t312 * t306;
t293 = Ifges(6,4) * t312 + Ifges(6,2) * t311 + Ifges(6,6) * t347;
t294 = Ifges(6,1) * t312 + Ifges(6,4) * t311 + Ifges(6,5) * t347;
t362 = mrSges(6,1) * t277 - mrSges(6,2) * t278 + Ifges(6,5) * t291 + Ifges(6,6) * t290 + Ifges(6,3) * t344 + t312 * t293 - t294 * t311;
t361 = m(5) * t307 - t316 * mrSges(5,1) + t317 * mrSges(5,2) - t327 * t321 + t328 * t322 + t366;
t319 = -qJDD(2) * pkin(2) - t359 * qJ(3) + t367;
t360 = -m(4) * t319 + mrSges(4,1) * t371 - t361 + (t345 * t359 + t376) * mrSges(4,3);
t310 = Ifges(5,1) * t328 + Ifges(5,4) * t327 + Ifges(5,5) * qJD(4);
t309 = Ifges(5,4) * t328 + Ifges(5,2) * t327 + Ifges(5,6) * qJD(4);
t308 = Ifges(5,5) * t328 + Ifges(5,6) * t327 + Ifges(5,3) * qJD(4);
t292 = Ifges(6,5) * t312 + Ifges(6,6) * t311 + Ifges(6,3) * t347;
t270 = qJDD(2) * t377 - t360;
t269 = mrSges(6,2) * t282 - mrSges(6,3) * t277 + Ifges(6,1) * t291 + Ifges(6,4) * t290 + Ifges(6,5) * t344 + t311 * t292 - t347 * t293;
t268 = -mrSges(6,1) * t282 + mrSges(6,3) * t278 + Ifges(6,4) * t291 + Ifges(6,2) * t290 + Ifges(6,6) * t344 - t312 * t292 + t347 * t294;
t259 = mrSges(5,2) * t307 - mrSges(5,3) * t283 + Ifges(5,1) * t317 + Ifges(5,4) * t316 + Ifges(5,5) * qJDD(4) - pkin(7) * t267 - qJD(4) * t309 - t353 * t268 + t356 * t269 + t327 * t308;
t258 = -mrSges(5,1) * t307 + mrSges(5,3) * t284 + Ifges(5,4) * t317 + Ifges(5,2) * t316 + Ifges(5,6) * qJDD(4) - pkin(4) * t366 + pkin(7) * t368 + qJD(4) * t310 + t356 * t268 + t353 * t269 - t328 * t308;
t1 = [m(2) * t348 + t355 * (m(3) * t325 - t359 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t370) + t358 * ((mrSges(3,1) - t377) * qJDD(2) - t359 * mrSges(3,2) + m(3) * t324 + t360); mrSges(3,1) * t324 - mrSges(3,2) * t325 + t349 * (mrSges(4,2) * t319 - mrSges(4,3) * t303 - pkin(6) * t375 - t354 * t258 + t357 * t259) + t351 * (-mrSges(4,1) * t319 + mrSges(4,3) * t304 - pkin(3) * t361 + pkin(6) * t369 + t357 * t258 + t354 * t259) - pkin(2) * t270 + qJ(3) * t370 + (Ifges(4,2) * t346 + Ifges(3,3) + (Ifges(4,1) * t349 + Ifges(4,4) * t379) * t349) * qJDD(2); t270; mrSges(5,1) * t283 - mrSges(5,2) * t284 + Ifges(5,5) * t317 + Ifges(5,6) * t316 + Ifges(5,3) * qJDD(4) + pkin(4) * t267 + t309 * t328 - t310 * t327 + t362; t362;];
tauJ = t1;
