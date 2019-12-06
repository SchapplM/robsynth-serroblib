% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:56
% EndTime: 2019-12-05 16:25:59
% DurationCPUTime: 1.24s
% Computational Cost: add. (7399->217), mult. (15661->286), div. (0->0), fcn. (10785->12), ass. (0->94)
t365 = sin(pkin(9));
t368 = cos(pkin(9));
t356 = g(1) * t365 - g(2) * t368;
t363 = -g(3) + qJDD(1);
t366 = sin(pkin(5));
t369 = cos(pkin(5));
t395 = t356 * t369 + t363 * t366;
t364 = sin(pkin(10));
t367 = cos(pkin(10));
t371 = sin(qJ(3));
t374 = cos(qJ(3));
t342 = (t364 * t371 - t367 * t374) * qJD(2);
t357 = -g(1) * t368 - g(2) * t365;
t372 = sin(qJ(2));
t375 = cos(qJ(2));
t324 = -t372 * t357 + t375 * t395;
t394 = 2 * qJD(4);
t325 = t375 * t357 + t372 * t395;
t377 = qJD(2) ^ 2;
t320 = -pkin(2) * t377 + qJDD(2) * pkin(7) + t325;
t339 = -t356 * t366 + t363 * t369;
t306 = -t371 * t320 + t374 * t339;
t389 = qJD(2) * qJD(3);
t388 = t374 * t389;
t354 = qJDD(2) * t371 + t388;
t303 = (-t354 + t388) * qJ(4) + (t371 * t374 * t377 + qJDD(3)) * pkin(3) + t306;
t307 = t374 * t320 + t371 * t339;
t355 = qJDD(2) * t374 - t371 * t389;
t390 = t371 * qJD(2);
t358 = qJD(3) * pkin(3) - qJ(4) * t390;
t362 = t374 ^ 2;
t304 = -pkin(3) * t362 * t377 + qJ(4) * t355 - qJD(3) * t358 + t307;
t299 = t303 * t364 + t304 * t367 - t342 * t394;
t343 = (t364 * t374 + t367 * t371) * qJD(2);
t327 = mrSges(5,1) * t342 + mrSges(5,2) * t343;
t331 = -t354 * t364 + t355 * t367;
t338 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t343;
t328 = pkin(4) * t342 - pkin(8) * t343;
t376 = qJD(3) ^ 2;
t297 = -pkin(4) * t376 + qJDD(3) * pkin(8) - t328 * t342 + t299;
t380 = -qJDD(2) * pkin(2) - t324;
t305 = -t355 * pkin(3) + qJDD(4) + t358 * t390 + (-qJ(4) * t362 - pkin(7)) * t377 + t380;
t332 = t354 * t367 + t355 * t364;
t300 = (qJD(3) * t342 - t332) * pkin(8) + (qJD(3) * t343 - t331) * pkin(4) + t305;
t370 = sin(qJ(5));
t373 = cos(qJ(5));
t294 = -t297 * t370 + t300 * t373;
t333 = qJD(3) * t373 - t343 * t370;
t314 = qJD(5) * t333 + qJDD(3) * t370 + t332 * t373;
t334 = qJD(3) * t370 + t343 * t373;
t315 = -mrSges(6,1) * t333 + mrSges(6,2) * t334;
t341 = qJD(5) + t342;
t317 = -mrSges(6,2) * t341 + mrSges(6,3) * t333;
t330 = qJDD(5) - t331;
t292 = m(6) * t294 + mrSges(6,1) * t330 - mrSges(6,3) * t314 - t315 * t334 + t317 * t341;
t295 = t297 * t373 + t300 * t370;
t313 = -qJD(5) * t334 + qJDD(3) * t373 - t332 * t370;
t318 = mrSges(6,1) * t341 - mrSges(6,3) * t334;
t293 = m(6) * t295 - mrSges(6,2) * t330 + mrSges(6,3) * t313 + t315 * t333 - t318 * t341;
t385 = -t292 * t370 + t293 * t373;
t282 = m(5) * t299 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t331 - qJD(3) * t338 - t327 * t342 + t385;
t384 = -t367 * t303 + t364 * t304;
t298 = -0.2e1 * qJD(4) * t343 - t384;
t337 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t342;
t296 = -qJDD(3) * pkin(4) - t376 * pkin(8) + (t394 + t328) * t343 + t384;
t381 = -m(6) * t296 + mrSges(6,1) * t313 - mrSges(6,2) * t314 + t317 * t333 - t318 * t334;
t288 = m(5) * t298 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t332 + qJD(3) * t337 - t327 * t343 + t381;
t279 = t282 * t364 + t288 * t367;
t284 = t292 * t373 + t293 * t370;
t391 = qJD(2) * t374;
t353 = (-mrSges(4,1) * t374 + mrSges(4,2) * t371) * qJD(2);
t360 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t391;
t277 = m(4) * t306 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t354 + qJD(3) * t360 - t353 * t390 + t279;
t359 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t390;
t386 = t282 * t367 - t288 * t364;
t278 = m(4) * t307 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t355 - qJD(3) * t359 + t353 * t391 + t386;
t387 = -t277 * t371 + t278 * t374;
t283 = m(5) * t305 - t331 * mrSges(5,1) + mrSges(5,2) * t332 + t342 * t337 + t338 * t343 + t284;
t309 = Ifges(6,4) * t334 + Ifges(6,2) * t333 + Ifges(6,6) * t341;
t310 = Ifges(6,1) * t334 + Ifges(6,4) * t333 + Ifges(6,5) * t341;
t379 = mrSges(6,1) * t294 - mrSges(6,2) * t295 + Ifges(6,5) * t314 + Ifges(6,6) * t313 + Ifges(6,3) * t330 + t309 * t334 - t310 * t333;
t319 = -t377 * pkin(7) + t380;
t378 = -m(4) * t319 + t355 * mrSges(4,1) - mrSges(4,2) * t354 - t359 * t390 + t360 * t391 - t283;
t347 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t371 + Ifges(4,4) * t374) * qJD(2);
t346 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t371 + Ifges(4,2) * t374) * qJD(2);
t323 = Ifges(5,1) * t343 - Ifges(5,4) * t342 + Ifges(5,5) * qJD(3);
t322 = Ifges(5,4) * t343 - Ifges(5,2) * t342 + Ifges(5,6) * qJD(3);
t321 = Ifges(5,5) * t343 - Ifges(5,6) * t342 + Ifges(5,3) * qJD(3);
t308 = Ifges(6,5) * t334 + Ifges(6,6) * t333 + Ifges(6,3) * t341;
t286 = mrSges(6,2) * t296 - mrSges(6,3) * t294 + Ifges(6,1) * t314 + Ifges(6,4) * t313 + Ifges(6,5) * t330 + t308 * t333 - t309 * t341;
t285 = -mrSges(6,1) * t296 + mrSges(6,3) * t295 + Ifges(6,4) * t314 + Ifges(6,2) * t313 + Ifges(6,6) * t330 - t308 * t334 + t310 * t341;
t275 = -mrSges(5,1) * t305 + mrSges(5,3) * t299 + Ifges(5,4) * t332 + Ifges(5,2) * t331 + Ifges(5,6) * qJDD(3) - pkin(4) * t284 + qJD(3) * t323 - t321 * t343 - t379;
t274 = mrSges(5,2) * t305 - mrSges(5,3) * t298 + Ifges(5,1) * t332 + Ifges(5,4) * t331 + Ifges(5,5) * qJDD(3) - pkin(8) * t284 - qJD(3) * t322 - t285 * t370 + t286 * t373 - t321 * t342;
t1 = [m(2) * t363 + t369 * (m(3) * t339 + t277 * t374 + t278 * t371) + (t372 * (m(3) * t325 - mrSges(3,1) * t377 - qJDD(2) * mrSges(3,2) + t387) + t375 * (m(3) * t324 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t377 + t378)) * t366; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t324 - mrSges(3,2) * t325 + t371 * (mrSges(4,2) * t319 - mrSges(4,3) * t306 + Ifges(4,1) * t354 + Ifges(4,4) * t355 + Ifges(4,5) * qJDD(3) - qJ(4) * t279 - qJD(3) * t346 + t274 * t367 - t275 * t364) + t374 * (-mrSges(4,1) * t319 + mrSges(4,3) * t307 + Ifges(4,4) * t354 + Ifges(4,2) * t355 + Ifges(4,6) * qJDD(3) - pkin(3) * t283 + qJ(4) * t386 + qJD(3) * t347 + t364 * t274 + t367 * t275) + pkin(2) * t378 + pkin(7) * t387; Ifges(4,5) * t354 + Ifges(4,6) * t355 + mrSges(4,1) * t306 - mrSges(4,2) * t307 + Ifges(5,5) * t332 + Ifges(5,6) * t331 + t343 * t322 + t342 * t323 + mrSges(5,1) * t298 - mrSges(5,2) * t299 + t370 * t286 + t373 * t285 + pkin(4) * t381 + pkin(8) * t385 + pkin(3) * t279 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t346 * t371 - t347 * t374) * qJD(2); t283; t379;];
tauJ = t1;
