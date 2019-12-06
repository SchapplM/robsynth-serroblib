% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:58
% EndTime: 2019-12-05 16:51:00
% DurationCPUTime: 0.92s
% Computational Cost: add. (3308->188), mult. (6577->233), div. (0->0), fcn. (4065->8), ass. (0->79)
t378 = Ifges(5,1) + Ifges(6,1);
t369 = Ifges(5,4) - Ifges(6,5);
t376 = Ifges(6,4) + Ifges(5,5);
t377 = Ifges(5,2) + Ifges(6,3);
t375 = Ifges(5,6) - Ifges(6,6);
t374 = Ifges(5,3) + Ifges(6,2);
t349 = sin(qJ(4));
t350 = sin(qJ(3));
t365 = t350 * qJD(2);
t352 = cos(qJ(3));
t366 = qJD(2) * t352;
t371 = cos(qJ(4));
t318 = t349 * t365 - t371 * t366;
t319 = (t349 * t352 + t371 * t350) * qJD(2);
t344 = qJD(3) + qJD(4);
t373 = t377 * t318 - t369 * t319 - t375 * t344;
t372 = -t369 * t318 + t378 * t319 + t376 * t344;
t370 = -mrSges(5,3) - mrSges(6,2);
t347 = sin(pkin(8));
t348 = cos(pkin(8));
t332 = -t348 * g(1) - t347 * g(2);
t346 = -g(3) + qJDD(1);
t351 = sin(qJ(2));
t353 = cos(qJ(2));
t314 = t353 * t332 + t351 * t346;
t354 = qJD(2) ^ 2;
t308 = -t354 * pkin(2) + qJDD(2) * pkin(6) + t314;
t331 = -t347 * g(1) + t348 * g(2);
t290 = -t350 * t308 + t352 * t331;
t364 = qJD(2) * qJD(3);
t362 = t352 * t364;
t329 = t350 * qJDD(2) + t362;
t276 = (-t329 + t362) * pkin(7) + (t350 * t352 * t354 + qJDD(3)) * pkin(3) + t290;
t291 = t352 * t308 + t350 * t331;
t330 = t352 * qJDD(2) - t350 * t364;
t335 = qJD(3) * pkin(3) - pkin(7) * t365;
t345 = t352 ^ 2;
t277 = -t345 * t354 * pkin(3) + t330 * pkin(7) - qJD(3) * t335 + t291;
t273 = t349 * t276 + t371 * t277;
t288 = t319 * qJD(4) + t349 * t329 - t371 * t330;
t310 = t344 * mrSges(5,1) - t319 * mrSges(5,3);
t343 = qJDD(3) + qJDD(4);
t301 = t318 * pkin(4) - t319 * qJ(5);
t342 = t344 ^ 2;
t267 = -t342 * pkin(4) + t343 * qJ(5) + 0.2e1 * qJD(5) * t344 - t318 * t301 + t273;
t311 = -t344 * mrSges(6,1) + t319 * mrSges(6,2);
t363 = m(6) * t267 + t343 * mrSges(6,3) + t344 * t311;
t302 = t318 * mrSges(6,1) - t319 * mrSges(6,3);
t367 = -t318 * mrSges(5,1) - t319 * mrSges(5,2) - t302;
t258 = m(5) * t273 - t343 * mrSges(5,2) + t370 * t288 - t344 * t310 + t367 * t318 + t363;
t272 = t371 * t276 - t349 * t277;
t289 = -t318 * qJD(4) + t371 * t329 + t349 * t330;
t309 = -t344 * mrSges(5,2) - t318 * mrSges(5,3);
t268 = -t343 * pkin(4) - t342 * qJ(5) + t319 * t301 + qJDD(5) - t272;
t312 = -t318 * mrSges(6,2) + t344 * mrSges(6,3);
t359 = -m(6) * t268 + t343 * mrSges(6,1) + t344 * t312;
t260 = m(5) * t272 + t343 * mrSges(5,1) + t370 * t289 + t344 * t309 + t367 * t319 + t359;
t255 = t349 * t258 + t371 * t260;
t368 = t375 * t318 - t376 * t319 - t374 * t344;
t328 = (-t352 * mrSges(4,1) + t350 * mrSges(4,2)) * qJD(2);
t333 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t365;
t334 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t366;
t360 = t371 * t258 - t349 * t260;
t361 = -t350 * (m(4) * t290 + qJDD(3) * mrSges(4,1) - t329 * mrSges(4,3) + qJD(3) * t334 - t328 * t365 + t255) + t352 * (m(4) * t291 - qJDD(3) * mrSges(4,2) + t330 * mrSges(4,3) - qJD(3) * t333 + t328 * t366 + t360);
t313 = -t351 * t332 + t353 * t346;
t358 = -qJDD(2) * pkin(2) - t313;
t278 = -t330 * pkin(3) + t335 * t365 + (-pkin(7) * t345 - pkin(6)) * t354 + t358;
t270 = -0.2e1 * qJD(5) * t319 + (t318 * t344 - t289) * qJ(5) + (t319 * t344 + t288) * pkin(4) + t278;
t261 = m(6) * t270 + t288 * mrSges(6,1) - t289 * mrSges(6,3) - t319 * t311 + t318 * t312;
t357 = m(5) * t278 + t288 * mrSges(5,1) + t289 * mrSges(5,2) + t318 * t309 + t319 * t310 + t261;
t264 = t289 * mrSges(6,2) + t319 * t302 - t359;
t356 = mrSges(5,1) * t272 - mrSges(6,1) * t268 - mrSges(5,2) * t273 + mrSges(6,3) * t267 - pkin(4) * t264 + qJ(5) * t363 + t374 * t343 - t373 * t319 + (-qJ(5) * t302 + t372) * t318 + t376 * t289 + (-qJ(5) * mrSges(6,2) - t375) * t288;
t307 = -t354 * pkin(6) + t358;
t355 = -m(4) * t307 + t330 * mrSges(4,1) - t329 * mrSges(4,2) - t333 * t365 + t334 * t366 - t357;
t317 = Ifges(4,5) * qJD(3) + (t350 * Ifges(4,1) + t352 * Ifges(4,4)) * qJD(2);
t316 = Ifges(4,6) * qJD(3) + (t350 * Ifges(4,4) + t352 * Ifges(4,2)) * qJD(2);
t252 = mrSges(5,2) * t278 + mrSges(6,2) * t268 - mrSges(5,3) * t272 - mrSges(6,3) * t270 - qJ(5) * t261 - t369 * t288 + t378 * t289 + t368 * t318 + t376 * t343 + t373 * t344;
t251 = -mrSges(5,1) * t278 - mrSges(6,1) * t270 + mrSges(6,2) * t267 + mrSges(5,3) * t273 - pkin(4) * t261 - t377 * t288 + t369 * t289 + t368 * t319 + t375 * t343 + t372 * t344;
t1 = [m(2) * t346 + t351 * (m(3) * t314 - t354 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t361) + t353 * (m(3) * t313 + qJDD(2) * mrSges(3,1) - t354 * mrSges(3,2) + t355); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t313 - mrSges(3,2) * t314 + t350 * (mrSges(4,2) * t307 - mrSges(4,3) * t290 + Ifges(4,1) * t329 + Ifges(4,4) * t330 + Ifges(4,5) * qJDD(3) - pkin(7) * t255 - qJD(3) * t316 - t349 * t251 + t371 * t252) + t352 * (-mrSges(4,1) * t307 + mrSges(4,3) * t291 + Ifges(4,4) * t329 + Ifges(4,2) * t330 + Ifges(4,6) * qJDD(3) - pkin(3) * t357 + pkin(7) * t360 + qJD(3) * t317 + t371 * t251 + t349 * t252) + pkin(2) * t355 + pkin(6) * t361; mrSges(4,1) * t290 - mrSges(4,2) * t291 + Ifges(4,3) * qJDD(3) + (t350 * t316 - t352 * t317) * qJD(2) + Ifges(4,5) * t329 + Ifges(4,6) * t330 + t356 + pkin(3) * t255; t356; t264;];
tauJ = t1;
