% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PRPR6
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:24
% EndTime: 2019-12-31 16:24:25
% DurationCPUTime: 0.89s
% Computational Cost: add. (7000->164), mult. (14958->213), div. (0->0), fcn. (9429->8), ass. (0->76)
t336 = qJD(2) ^ 2;
t331 = cos(pkin(7));
t358 = pkin(3) * t331;
t329 = sin(pkin(7));
t357 = mrSges(4,2) * t329;
t356 = cos(pkin(6));
t327 = t331 ^ 2;
t355 = t327 * t336;
t330 = sin(pkin(6));
t318 = -t356 * g(1) - t330 * g(2);
t328 = -g(3) + qJDD(1);
t333 = sin(qJ(2));
t335 = cos(qJ(2));
t308 = t335 * t318 + t333 * t328;
t304 = -t336 * pkin(2) + qJDD(2) * qJ(3) + t308;
t317 = t330 * g(1) - t356 * g(2);
t351 = qJD(2) * qJD(3);
t353 = -t331 * t317 - 0.2e1 * t329 * t351;
t289 = (-pkin(5) * qJDD(2) + t336 * t358 - t304) * t329 + t353;
t292 = -t329 * t317 + (t304 + 0.2e1 * t351) * t331;
t350 = qJDD(2) * t331;
t290 = -pkin(3) * t355 + pkin(5) * t350 + t292;
t332 = sin(qJ(4));
t334 = cos(qJ(4));
t287 = t334 * t289 - t332 * t290;
t340 = -t329 * t332 + t331 * t334;
t309 = t340 * qJD(2);
t341 = t329 * t334 + t331 * t332;
t310 = t341 * qJD(2);
t298 = -t309 * mrSges(5,1) + t310 * mrSges(5,2);
t301 = t309 * qJD(4) + t341 * qJDD(2);
t305 = -qJD(4) * mrSges(5,2) + t309 * mrSges(5,3);
t285 = m(5) * t287 + qJDD(4) * mrSges(5,1) - t301 * mrSges(5,3) + qJD(4) * t305 - t310 * t298;
t288 = t332 * t289 + t334 * t290;
t300 = -t310 * qJD(4) + t340 * qJDD(2);
t306 = qJD(4) * mrSges(5,1) - t310 * mrSges(5,3);
t286 = m(5) * t288 - qJDD(4) * mrSges(5,2) + t300 * mrSges(5,3) - qJD(4) * t306 + t309 * t298;
t277 = t334 * t285 + t332 * t286;
t291 = -t329 * t304 + t353;
t339 = mrSges(4,3) * qJDD(2) + t336 * (-mrSges(4,1) * t331 + t357);
t275 = m(4) * t291 - t339 * t329 + t277;
t346 = -t332 * t285 + t334 * t286;
t276 = m(4) * t292 + t339 * t331 + t346;
t347 = -t329 * t275 + t331 * t276;
t270 = m(3) * t308 - t336 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t347;
t307 = -t333 * t318 + t335 * t328;
t342 = qJDD(3) - t307;
t303 = -qJDD(2) * pkin(2) - t336 * qJ(3) + t342;
t326 = t329 ^ 2;
t293 = (-pkin(2) - t358) * qJDD(2) + (-qJ(3) + (-t326 - t327) * pkin(5)) * t336 + t342;
t338 = m(5) * t293 - t300 * mrSges(5,1) + t301 * mrSges(5,2) - t309 * t305 + t310 * t306;
t337 = -m(4) * t303 + mrSges(4,1) * t350 - t338 + (t326 * t336 + t355) * mrSges(4,3);
t281 = m(3) * t307 - t336 * mrSges(3,2) + (mrSges(3,1) - t357) * qJDD(2) + t337;
t348 = t335 * t270 - t333 * t281;
t265 = m(2) * t318 + t348;
t273 = t331 * t275 + t329 * t276;
t272 = (m(2) + m(3)) * t317 - t273;
t354 = t330 * t265 + t356 * t272;
t266 = t333 * t270 + t335 * t281;
t343 = Ifges(4,5) * t329 + Ifges(4,6) * t331;
t352 = t336 * t343;
t349 = t356 * t265 - t330 * t272;
t345 = Ifges(4,1) * t329 + Ifges(4,4) * t331;
t344 = Ifges(4,4) * t329 + Ifges(4,2) * t331;
t296 = Ifges(5,1) * t310 + Ifges(5,4) * t309 + Ifges(5,5) * qJD(4);
t295 = Ifges(5,4) * t310 + Ifges(5,2) * t309 + Ifges(5,6) * qJD(4);
t294 = Ifges(5,5) * t310 + Ifges(5,6) * t309 + Ifges(5,3) * qJD(4);
t279 = mrSges(5,2) * t293 - mrSges(5,3) * t287 + Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * qJDD(4) - qJD(4) * t295 + t309 * t294;
t278 = -mrSges(5,1) * t293 + mrSges(5,3) * t288 + Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * qJDD(4) + qJD(4) * t296 - t310 * t294;
t267 = mrSges(4,2) * t303 - mrSges(4,3) * t291 - pkin(5) * t277 + t345 * qJDD(2) - t332 * t278 + t334 * t279 + t331 * t352;
t262 = -mrSges(4,1) * t303 + mrSges(4,3) * t292 - pkin(3) * t338 + pkin(5) * t346 + t344 * qJDD(2) + t334 * t278 + t332 * t279 - t329 * t352;
t261 = mrSges(3,1) * t317 - mrSges(4,1) * t291 - mrSges(5,1) * t287 + mrSges(4,2) * t292 + mrSges(5,2) * t288 + mrSges(3,3) * t308 - Ifges(5,5) * t301 - Ifges(5,6) * t300 - Ifges(5,3) * qJDD(4) - pkin(2) * t273 - pkin(3) * t277 - t310 * t295 + t309 * t296 + (Ifges(3,6) - t343) * qJDD(2) + (-t329 * t344 + t331 * t345 + Ifges(3,5)) * t336;
t260 = -mrSges(3,2) * t317 - mrSges(3,3) * t307 + Ifges(3,5) * qJDD(2) - t336 * Ifges(3,6) - qJ(3) * t273 - t329 * t262 + t331 * t267;
t259 = -mrSges(2,1) * t328 + mrSges(2,3) * t318 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t307 + mrSges(3,2) * t308 - t329 * t267 - t331 * t262 - pkin(2) * (-qJDD(2) * t357 + t337) - qJ(3) * t347 - pkin(1) * t266;
t258 = mrSges(2,2) * t328 - mrSges(2,3) * t317 - pkin(4) * t266 + t335 * t260 - t333 * t261;
t1 = [-m(1) * g(1) + t349; -m(1) * g(2) + t354; -m(1) * g(3) + m(2) * t328 + t266; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t354 + t356 * t258 - t330 * t259; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t349 + t330 * t258 + t356 * t259; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t317 - mrSges(2,2) * t318 + t333 * t260 + t335 * t261 + pkin(1) * (m(3) * t317 - t273) + pkin(4) * t348;];
tauB = t1;
