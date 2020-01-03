% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:24
% EndTime: 2019-12-31 18:56:25
% DurationCPUTime: 0.99s
% Computational Cost: add. (3283->193), mult. (6217->230), div. (0->0), fcn. (3398->6), ass. (0->80)
t380 = Ifges(5,1) + Ifges(6,1);
t370 = Ifges(5,4) + Ifges(6,4);
t369 = Ifges(5,5) + Ifges(6,5);
t379 = Ifges(5,2) + Ifges(6,2);
t368 = Ifges(5,6) + Ifges(6,6);
t339 = sin(qJ(4));
t342 = cos(qJ(4));
t343 = cos(qJ(3));
t362 = t343 * qJD(1);
t327 = qJD(3) * t342 - t339 * t362;
t340 = sin(qJ(3));
t361 = qJD(1) * qJD(3);
t357 = t340 * t361;
t332 = qJDD(1) * t343 - t357;
t304 = qJD(4) * t327 + qJDD(3) * t339 + t332 * t342;
t328 = qJD(3) * t339 + t342 * t362;
t306 = -mrSges(6,1) * t327 + mrSges(6,2) * t328;
t346 = qJD(1) ^ 2;
t372 = -pkin(1) - pkin(6);
t341 = sin(qJ(1));
t344 = cos(qJ(1));
t352 = -t344 * g(1) - t341 * g(2);
t373 = -qJDD(1) * qJ(2) - 0.2e1 * qJD(2) * qJD(1) - t352;
t316 = t346 * t372 - t373;
t356 = t343 * t361;
t331 = -qJDD(1) * t340 - t356;
t287 = (-t332 + t357) * pkin(7) + (-t331 + t356) * pkin(3) + t316;
t355 = t341 * g(1) - t344 * g(2);
t349 = -t346 * qJ(2) + qJDD(2) - t355;
t317 = qJDD(1) * t372 + t349;
t309 = -g(3) * t343 + t340 * t317;
t330 = (pkin(3) * t340 - pkin(7) * t343) * qJD(1);
t345 = qJD(3) ^ 2;
t363 = t340 * qJD(1);
t290 = -pkin(3) * t345 + qJDD(3) * pkin(7) - t330 * t363 + t309;
t283 = t287 * t342 - t339 * t290;
t326 = qJDD(4) - t331;
t335 = qJD(4) + t363;
t279 = -0.2e1 * qJD(5) * t328 + (t327 * t335 - t304) * qJ(5) + (t327 * t328 + t326) * pkin(4) + t283;
t310 = -mrSges(6,2) * t335 + mrSges(6,3) * t327;
t359 = m(6) * t279 + t326 * mrSges(6,1) + t310 * t335;
t276 = -t304 * mrSges(6,3) - t328 * t306 + t359;
t284 = t287 * t339 + t290 * t342;
t303 = -qJD(4) * t328 + qJDD(3) * t342 - t332 * t339;
t312 = pkin(4) * t335 - qJ(5) * t328;
t325 = t327 ^ 2;
t281 = -pkin(4) * t325 + qJ(5) * t303 + 0.2e1 * qJD(5) * t327 - t312 * t335 + t284;
t366 = -t327 * t379 - t328 * t370 - t335 * t368;
t374 = t327 * t370 + t328 * t380 + t335 * t369;
t376 = Ifges(5,3) + Ifges(6,3);
t378 = mrSges(5,1) * t283 + mrSges(6,1) * t279 - mrSges(5,2) * t284 - mrSges(6,2) * t281 + pkin(4) * t276 + t368 * t303 + t369 * t304 + t326 * t376 - t327 * t374 - t366 * t328;
t371 = -mrSges(5,2) - mrSges(6,2);
t307 = -mrSges(5,1) * t327 + mrSges(5,2) * t328;
t311 = -mrSges(5,2) * t335 + mrSges(5,3) * t327;
t271 = m(5) * t283 + t326 * mrSges(5,1) + t335 * t311 + (-t306 - t307) * t328 + (-mrSges(5,3) - mrSges(6,3)) * t304 + t359;
t358 = m(6) * t281 + mrSges(6,3) * t303 + t306 * t327;
t313 = mrSges(6,1) * t335 - mrSges(6,3) * t328;
t364 = -mrSges(5,1) * t335 + mrSges(5,3) * t328 - t313;
t274 = m(5) * t284 + t303 * mrSges(5,3) + t327 * t307 + t326 * t371 + t335 * t364 + t358;
t269 = t271 * t342 + t274 * t339;
t367 = -t327 * t368 - t328 * t369 - t335 * t376;
t354 = -t271 * t339 + t274 * t342;
t308 = g(3) * t340 + t317 * t343;
t289 = -qJDD(3) * pkin(3) - pkin(7) * t345 + t330 * t362 - t308;
t282 = -pkin(4) * t303 - qJ(5) * t325 + t312 * t328 + qJDD(5) + t289;
t353 = -m(6) * t282 + mrSges(6,1) * t303 + t310 * t327;
t329 = (mrSges(4,1) * t340 + mrSges(4,2) * t343) * qJD(1);
t333 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t363;
t334 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t362;
t347 = -m(5) * t289 + mrSges(5,1) * t303 + t304 * t371 + t311 * t327 + t328 * t364 + t353;
t351 = t340 * (m(4) * t309 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t331 - qJD(3) * t334 - t329 * t363 + t354) + t343 * (m(4) * t308 + qJDD(3) * mrSges(4,1) - t332 * mrSges(4,3) + qJD(3) * t333 - t329 * t362 + t347);
t323 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t343 - Ifges(4,4) * t340) * qJD(1);
t322 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t343 - Ifges(4,2) * t340) * qJD(1);
t319 = -qJDD(1) * pkin(1) + t349;
t318 = t346 * pkin(1) + t373;
t277 = t304 * mrSges(6,2) + t328 * t313 - t353;
t267 = mrSges(5,2) * t289 + mrSges(6,2) * t282 - mrSges(5,3) * t283 - mrSges(6,3) * t279 - qJ(5) * t276 + t370 * t303 + t304 * t380 + t369 * t326 - t367 * t327 + t366 * t335;
t266 = -mrSges(5,1) * t289 + mrSges(5,3) * t284 - mrSges(6,1) * t282 + mrSges(6,3) * t281 - pkin(4) * t277 + qJ(5) * t358 + (-qJ(5) * t313 + t374) * t335 + t367 * t328 + (-mrSges(6,2) * qJ(5) + t368) * t326 + t370 * t304 + t379 * t303;
t265 = m(3) * t319 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t346 + t351;
t1 = [mrSges(2,1) * t355 - mrSges(2,2) * t352 + mrSges(3,2) * t319 - mrSges(3,3) * t318 + t343 * (mrSges(4,2) * t316 - mrSges(4,3) * t308 + Ifges(4,1) * t332 + Ifges(4,4) * t331 + Ifges(4,5) * qJDD(3) - pkin(7) * t269 - qJD(3) * t322 - t266 * t339 + t267 * t342) - t340 * (-mrSges(4,1) * t316 + mrSges(4,3) * t309 + Ifges(4,4) * t332 + Ifges(4,2) * t331 + Ifges(4,6) * qJDD(3) - pkin(3) * t269 + qJD(3) * t323 - t378) - pkin(6) * t351 - pkin(1) * t265 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t318 + m(4) * t316 - t331 * mrSges(4,1) + t346 * mrSges(3,2) + t332 * mrSges(4,2) + t269 + qJDD(1) * mrSges(3,3) + (t333 * t340 + t334 * t343) * qJD(1)) * qJ(2); t265; Ifges(4,5) * t332 + Ifges(4,6) * t331 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t308 - mrSges(4,2) * t309 + t339 * t267 + t342 * t266 + pkin(3) * t347 + pkin(7) * t354 + (t322 * t343 + t323 * t340) * qJD(1); t378; t277;];
tauJ = t1;
