% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:11
% EndTime: 2019-12-31 16:36:14
% DurationCPUTime: 1.53s
% Computational Cost: add. (15957->194), mult. (29954->258), div. (0->0), fcn. (19424->10), ass. (0->88)
t354 = sin(pkin(8));
t356 = cos(pkin(8));
t347 = t354 * g(1) - t356 * g(2);
t348 = -t356 * g(1) - t354 * g(2);
t353 = -g(3) + qJDD(1);
t355 = sin(pkin(4));
t357 = cos(pkin(4));
t360 = sin(qJ(2));
t363 = cos(qJ(2));
t318 = -t360 * t348 + (t347 * t357 + t353 * t355) * t363;
t365 = qJD(2) ^ 2;
t380 = t357 * t360;
t381 = t355 * t360;
t319 = t347 * t380 + t363 * t348 + t353 * t381;
t317 = -t365 * pkin(2) + qJDD(2) * pkin(6) + t319;
t331 = -t355 * t347 + t357 * t353;
t359 = sin(qJ(3));
t362 = cos(qJ(3));
t314 = t362 * t317 + t359 * t331;
t344 = (-pkin(3) * t362 - pkin(7) * t359) * qJD(2);
t364 = qJD(3) ^ 2;
t376 = t362 * qJD(2);
t311 = -t364 * pkin(3) + qJDD(3) * pkin(7) + t344 * t376 + t314;
t316 = -qJDD(2) * pkin(2) - t365 * pkin(6) - t318;
t375 = qJD(2) * qJD(3);
t373 = t362 * t375;
t345 = t359 * qJDD(2) + t373;
t374 = t359 * t375;
t346 = t362 * qJDD(2) - t374;
t312 = (-t345 - t373) * pkin(7) + (-t346 + t374) * pkin(3) + t316;
t358 = sin(qJ(4));
t361 = cos(qJ(4));
t308 = -t358 * t311 + t361 * t312;
t377 = qJD(2) * t359;
t341 = t361 * qJD(3) - t358 * t377;
t326 = t341 * qJD(4) + t358 * qJDD(3) + t361 * t345;
t342 = t358 * qJD(3) + t361 * t377;
t327 = -t341 * mrSges(5,1) + t342 * mrSges(5,2);
t352 = qJD(4) - t376;
t329 = -t352 * mrSges(5,2) + t341 * mrSges(5,3);
t338 = qJDD(4) - t346;
t306 = m(5) * t308 + t338 * mrSges(5,1) - t326 * mrSges(5,3) - t342 * t327 + t352 * t329;
t309 = t361 * t311 + t358 * t312;
t325 = -t342 * qJD(4) + t361 * qJDD(3) - t358 * t345;
t330 = t352 * mrSges(5,1) - t342 * mrSges(5,3);
t307 = m(5) * t309 - t338 * mrSges(5,2) + t325 * mrSges(5,3) + t341 * t327 - t352 * t330;
t300 = t361 * t306 + t358 * t307;
t349 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t377;
t350 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t376;
t366 = -m(4) * t316 + t346 * mrSges(4,1) - t345 * mrSges(4,2) - t349 * t377 + t350 * t376 - t300;
t296 = m(3) * t318 + qJDD(2) * mrSges(3,1) - t365 * mrSges(3,2) + t366;
t382 = t296 * t363;
t379 = t362 * t331;
t343 = (-mrSges(4,1) * t362 + mrSges(4,2) * t359) * qJD(2);
t370 = -t358 * t306 + t361 * t307;
t299 = m(4) * t314 - qJDD(3) * mrSges(4,2) + t346 * mrSges(4,3) - qJD(3) * t349 + t343 * t376 + t370;
t313 = -t359 * t317 + t379;
t310 = -qJDD(3) * pkin(3) - t364 * pkin(7) - t379 + (qJD(2) * t344 + t317) * t359;
t367 = -m(5) * t310 + t325 * mrSges(5,1) - t326 * mrSges(5,2) + t341 * t329 - t342 * t330;
t304 = m(4) * t313 + qJDD(3) * mrSges(4,1) - t345 * mrSges(4,3) + qJD(3) * t350 - t343 * t377 + t367;
t371 = t362 * t299 - t359 * t304;
t290 = m(3) * t319 - t365 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t371;
t293 = t359 * t299 + t362 * t304;
t292 = m(3) * t331 + t293;
t280 = t290 * t380 - t355 * t292 + t357 * t382;
t278 = m(2) * t347 + t280;
t284 = t363 * t290 - t360 * t296;
t283 = m(2) * t348 + t284;
t378 = t356 * t278 + t354 * t283;
t279 = t290 * t381 + t357 * t292 + t355 * t382;
t372 = -t354 * t278 + t356 * t283;
t320 = Ifges(5,5) * t342 + Ifges(5,6) * t341 + Ifges(5,3) * t352;
t322 = Ifges(5,1) * t342 + Ifges(5,4) * t341 + Ifges(5,5) * t352;
t301 = -mrSges(5,1) * t310 + mrSges(5,3) * t309 + Ifges(5,4) * t326 + Ifges(5,2) * t325 + Ifges(5,6) * t338 - t342 * t320 + t352 * t322;
t321 = Ifges(5,4) * t342 + Ifges(5,2) * t341 + Ifges(5,6) * t352;
t302 = mrSges(5,2) * t310 - mrSges(5,3) * t308 + Ifges(5,1) * t326 + Ifges(5,4) * t325 + Ifges(5,5) * t338 + t341 * t320 - t352 * t321;
t333 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t359 + Ifges(4,6) * t362) * qJD(2);
t334 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t359 + Ifges(4,2) * t362) * qJD(2);
t285 = mrSges(4,2) * t316 - mrSges(4,3) * t313 + Ifges(4,1) * t345 + Ifges(4,4) * t346 + Ifges(4,5) * qJDD(3) - pkin(7) * t300 - qJD(3) * t334 - t358 * t301 + t361 * t302 + t333 * t376;
t335 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t359 + Ifges(4,4) * t362) * qJD(2);
t286 = -mrSges(4,1) * t316 - mrSges(5,1) * t308 + mrSges(5,2) * t309 + mrSges(4,3) * t314 + Ifges(4,4) * t345 - Ifges(5,5) * t326 + Ifges(4,2) * t346 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t325 - Ifges(5,3) * t338 - pkin(3) * t300 + qJD(3) * t335 - t342 * t321 + t341 * t322 - t333 * t377;
t275 = mrSges(3,2) * t331 - mrSges(3,3) * t318 + Ifges(3,5) * qJDD(2) - t365 * Ifges(3,6) - pkin(6) * t293 + t362 * t285 - t359 * t286;
t276 = Ifges(3,6) * qJDD(2) + t365 * Ifges(3,5) - mrSges(3,1) * t331 + mrSges(3,3) * t319 - Ifges(4,5) * t345 - Ifges(4,6) * t346 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t313 + mrSges(4,2) * t314 - t358 * t302 - t361 * t301 - pkin(3) * t367 - pkin(7) * t370 - pkin(2) * t293 + (-t359 * t334 + t362 * t335) * qJD(2);
t368 = pkin(5) * t284 + t275 * t360 + t276 * t363;
t274 = mrSges(3,1) * t318 - mrSges(3,2) * t319 + Ifges(3,3) * qJDD(2) + pkin(2) * t366 + pkin(6) * t371 + t359 * t285 + t362 * t286;
t273 = mrSges(2,2) * t353 - mrSges(2,3) * t347 + t363 * t275 - t360 * t276 + (-t279 * t355 - t280 * t357) * pkin(5);
t272 = -mrSges(2,1) * t353 + mrSges(2,3) * t348 - pkin(1) * t279 - t355 * t274 + t368 * t357;
t1 = [-m(1) * g(1) + t372; -m(1) * g(2) + t378; -m(1) * g(3) + m(2) * t353 + t279; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t378 - t354 * t272 + t356 * t273; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t372 + t356 * t272 + t354 * t273; -mrSges(1,1) * g(2) + mrSges(2,1) * t347 + mrSges(1,2) * g(1) - mrSges(2,2) * t348 + pkin(1) * t280 + t357 * t274 + t368 * t355;];
tauB = t1;
