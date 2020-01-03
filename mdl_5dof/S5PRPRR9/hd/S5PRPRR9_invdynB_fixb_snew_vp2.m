% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:41
% EndTime: 2019-12-31 17:39:42
% DurationCPUTime: 0.86s
% Computational Cost: add. (10383->172), mult. (14294->209), div. (0->0), fcn. (6962->8), ass. (0->75)
t334 = -pkin(2) - pkin(3);
t333 = -mrSges(3,1) - mrSges(4,1);
t332 = Ifges(4,4) + Ifges(3,5);
t331 = Ifges(3,6) - Ifges(4,6);
t302 = -qJD(2) + qJD(4);
t310 = sin(qJ(5));
t330 = t302 * t310;
t313 = cos(qJ(5));
t329 = t302 * t313;
t308 = sin(pkin(8));
t309 = cos(pkin(8));
t295 = t308 * g(1) - t309 * g(2);
t296 = -t309 * g(1) - t308 * g(2);
t312 = sin(qJ(2));
t315 = cos(qJ(2));
t281 = t312 * t295 + t315 * t296;
t316 = qJD(2) ^ 2;
t322 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t281;
t278 = -t316 * pkin(2) + t322;
t275 = t334 * t316 + t322;
t280 = t315 * t295 - t312 * t296;
t319 = -t316 * qJ(3) + qJDD(3) - t280;
t277 = t334 * qJDD(2) + t319;
t311 = sin(qJ(4));
t314 = cos(qJ(4));
t272 = t314 * t275 + t311 * t277;
t300 = t302 ^ 2;
t301 = -qJDD(2) + qJDD(4);
t270 = -t300 * pkin(4) + t301 * pkin(7) + t272;
t306 = g(3) - qJDD(1);
t267 = -t310 * t270 + t313 * t306;
t287 = (-mrSges(6,1) * t313 + mrSges(6,2) * t310) * t302;
t327 = qJD(5) * t302;
t288 = t310 * t301 + t313 * t327;
t294 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t329;
t265 = m(6) * t267 + qJDD(5) * mrSges(6,1) - t288 * mrSges(6,3) + qJD(5) * t294 - t287 * t330;
t268 = t313 * t270 + t310 * t306;
t289 = t313 * t301 - t310 * t327;
t293 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t330;
t266 = m(6) * t268 - qJDD(5) * mrSges(6,2) + t289 * mrSges(6,3) - qJD(5) * t293 + t287 * t329;
t323 = -t310 * t265 + t313 * t266;
t258 = m(5) * t272 - t300 * mrSges(5,1) - t301 * mrSges(5,2) + t323;
t271 = -t311 * t275 + t314 * t277;
t269 = -t301 * pkin(4) - t300 * pkin(7) - t271;
t317 = -m(6) * t269 + t289 * mrSges(6,1) - t288 * mrSges(6,2) - t293 * t330 + t294 * t329;
t263 = m(5) * t271 + t301 * mrSges(5,1) - t300 * mrSges(5,2) + t317;
t324 = t314 * t258 - t311 * t263;
t321 = m(4) * t278 + qJDD(2) * mrSges(4,3) + t324;
t253 = m(3) * t281 - qJDD(2) * mrSges(3,2) + t333 * t316 + t321;
t256 = t311 * t258 + t314 * t263;
t279 = -qJDD(2) * pkin(2) + t319;
t318 = -m(4) * t279 + qJDD(2) * mrSges(4,1) + t316 * mrSges(4,3) - t256;
t254 = m(3) * t280 + qJDD(2) * mrSges(3,1) - t316 * mrSges(3,2) + t318;
t248 = t312 * t253 + t315 * t254;
t246 = m(2) * t295 + t248;
t325 = t315 * t253 - t312 * t254;
t247 = m(2) * t296 + t325;
t328 = t309 * t246 + t308 * t247;
t326 = -t308 * t246 + t309 * t247;
t260 = t313 * t265 + t310 * t266;
t320 = -m(4) * t306 - t260;
t259 = -m(5) * t306 + t320;
t298 = m(3) * t306;
t284 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t310 + Ifges(6,4) * t313) * t302;
t283 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t310 + Ifges(6,2) * t313) * t302;
t282 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t310 + Ifges(6,6) * t313) * t302;
t262 = mrSges(6,2) * t269 - mrSges(6,3) * t267 + Ifges(6,1) * t288 + Ifges(6,4) * t289 + Ifges(6,5) * qJDD(5) - qJD(5) * t283 + t282 * t329;
t261 = -mrSges(6,1) * t269 + mrSges(6,3) * t268 + Ifges(6,4) * t288 + Ifges(6,2) * t289 + Ifges(6,6) * qJDD(5) + qJD(5) * t284 - t282 * t330;
t255 = -mrSges(5,1) * t306 - mrSges(6,1) * t267 + mrSges(6,2) * t268 + mrSges(5,3) * t272 + t300 * Ifges(5,5) - Ifges(6,5) * t288 + Ifges(5,6) * t301 - Ifges(6,6) * t289 - Ifges(6,3) * qJDD(5) - pkin(4) * t260 + (-t283 * t310 + t284 * t313) * t302;
t249 = mrSges(5,2) * t306 - mrSges(5,3) * t271 + Ifges(5,5) * t301 - t300 * Ifges(5,6) - pkin(7) * t260 - t310 * t261 + t313 * t262;
t242 = mrSges(4,2) * t279 - mrSges(3,3) * t280 - pkin(6) * t256 - qJ(3) * t259 + t314 * t249 - t311 * t255 - t331 * t316 + (-mrSges(3,2) + mrSges(4,3)) * t306 + t332 * qJDD(2);
t241 = mrSges(3,3) * t281 + mrSges(4,2) * t278 - t311 * t249 - t314 * t255 + pkin(3) * t260 - pkin(6) * t324 - pkin(2) * t259 + t332 * t316 + (pkin(3) * m(5) - t333) * t306 + t331 * qJDD(2);
t240 = -mrSges(2,2) * t306 - mrSges(2,3) * t295 - pkin(5) * t248 - t312 * t241 + t315 * t242;
t239 = mrSges(2,1) * t306 + mrSges(2,3) * t296 + t312 * t242 + t315 * t241 - pkin(1) * (t259 - t298) + pkin(5) * t325;
t1 = [-m(1) * g(1) + t326; -m(1) * g(2) + t328; -m(1) * g(3) - t298 + (-m(2) - m(5)) * t306 + t320; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t328 - t308 * t239 + t309 * t240; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t326 + t309 * t239 + t308 * t240; pkin(1) * t248 + mrSges(2,1) * t295 - mrSges(2,2) * t296 + pkin(2) * t318 + qJ(3) * (-t316 * mrSges(4,1) + t321) + mrSges(3,1) * t280 - mrSges(3,2) * t281 - pkin(3) * t256 - mrSges(4,1) * t279 + mrSges(4,3) * t278 - mrSges(5,1) * t271 + mrSges(5,2) * t272 - t310 * t262 - t313 * t261 - pkin(4) * t317 - pkin(7) * t323 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(5,3) * t301 + (Ifges(4,2) + Ifges(3,3)) * qJDD(2);];
tauB = t1;
