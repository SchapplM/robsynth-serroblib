% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:30
% EndTime: 2019-12-05 17:49:31
% DurationCPUTime: 0.56s
% Computational Cost: add. (4896->137), mult. (7045->180), div. (0->0), fcn. (4082->10), ass. (0->73)
t312 = qJD(1) + qJD(3);
t308 = t312 ^ 2;
t316 = cos(pkin(9));
t344 = pkin(4) * t316;
t314 = sin(pkin(9));
t343 = mrSges(5,2) * t314;
t311 = t316 ^ 2;
t342 = t308 * t311;
t309 = qJDD(1) + qJDD(3);
t341 = t309 * t316;
t320 = sin(qJ(1));
t323 = cos(qJ(1));
t336 = t323 * g(2) + t320 * g(3);
t295 = qJDD(1) * pkin(1) + t336;
t324 = qJD(1) ^ 2;
t334 = t320 * g(2) - t323 * g(3);
t296 = -t324 * pkin(1) + t334;
t315 = sin(pkin(8));
t317 = cos(pkin(8));
t284 = t317 * t295 - t315 * t296;
t282 = qJDD(1) * pkin(2) + t284;
t285 = t315 * t295 + t317 * t296;
t283 = -t324 * pkin(2) + t285;
t319 = sin(qJ(3));
t322 = cos(qJ(3));
t270 = t319 * t282 + t322 * t283;
t267 = -t308 * pkin(3) + t309 * qJ(4) + t270;
t313 = -g(1) + qJDD(2);
t335 = qJD(4) * t312;
t337 = t316 * t313 - 0.2e1 * t314 * t335;
t263 = -t314 * t267 + t337;
t330 = mrSges(5,3) * t309 + (-mrSges(5,1) * t316 + t343) * t308;
t260 = (-pkin(7) * t309 + t308 * t344 - t267) * t314 + t337;
t264 = t314 * t313 + (t267 + 0.2e1 * t335) * t316;
t261 = -pkin(4) * t342 + pkin(7) * t341 + t264;
t318 = sin(qJ(5));
t321 = cos(qJ(5));
t258 = t321 * t260 - t318 * t261;
t328 = -t314 * t318 + t316 * t321;
t288 = t328 * t312;
t329 = t314 * t321 + t316 * t318;
t289 = t329 * t312;
t276 = -t288 * mrSges(6,1) + t289 * mrSges(6,2);
t278 = t288 * qJD(5) + t329 * t309;
t286 = -qJD(5) * mrSges(6,2) + t288 * mrSges(6,3);
t256 = m(6) * t258 + qJDD(5) * mrSges(6,1) - t278 * mrSges(6,3) + qJD(5) * t286 - t289 * t276;
t259 = t318 * t260 + t321 * t261;
t277 = -t289 * qJD(5) + t328 * t309;
t287 = qJD(5) * mrSges(6,1) - t289 * mrSges(6,3);
t257 = m(6) * t259 - qJDD(5) * mrSges(6,2) + t277 * mrSges(6,3) - qJD(5) * t287 + t288 * t276;
t338 = t321 * t256 + t318 * t257;
t245 = m(5) * t263 - t330 * t314 + t338;
t332 = -t318 * t256 + t321 * t257;
t246 = m(5) * t264 + t330 * t316 + t332;
t333 = -t314 * t245 + t316 * t246;
t242 = m(4) * t270 - t308 * mrSges(4,1) - t309 * mrSges(4,2) + t333;
t269 = t322 * t282 - t319 * t283;
t331 = qJDD(4) - t269;
t266 = -t309 * pkin(3) - t308 * qJ(4) + t331;
t310 = t314 ^ 2;
t262 = (-pkin(3) - t344) * t309 + (-qJ(4) + (-t310 - t311) * pkin(7)) * t308 + t331;
t326 = m(6) * t262 - t277 * mrSges(6,1) + t278 * mrSges(6,2) - t288 * t286 + t289 * t287;
t325 = -m(5) * t266 + mrSges(5,1) * t341 - t326 + (t308 * t310 + t342) * mrSges(5,3);
t250 = m(4) * t269 - t308 * mrSges(4,2) + (mrSges(4,1) - t343) * t309 + t325;
t339 = t319 * t242 + t322 * t250;
t271 = Ifges(6,5) * t289 + Ifges(6,6) * t288 + Ifges(6,3) * qJD(5);
t273 = Ifges(6,1) * t289 + Ifges(6,4) * t288 + Ifges(6,5) * qJD(5);
t247 = -mrSges(6,1) * t262 + mrSges(6,3) * t259 + Ifges(6,4) * t278 + Ifges(6,2) * t277 + Ifges(6,6) * qJDD(5) + qJD(5) * t273 - t289 * t271;
t272 = Ifges(6,4) * t289 + Ifges(6,2) * t288 + Ifges(6,6) * qJD(5);
t248 = mrSges(6,2) * t262 - mrSges(6,3) * t258 + Ifges(6,1) * t278 + Ifges(6,4) * t277 + Ifges(6,5) * qJDD(5) - qJD(5) * t272 + t288 * t271;
t252 = t309 * t343 - t325;
t327 = -mrSges(4,2) * t270 + t316 * (-mrSges(5,1) * t266 + mrSges(5,3) * t264 + t318 * t248 + t321 * t247 - pkin(4) * t326 + pkin(7) * t332 + (Ifges(5,4) * t314 + Ifges(5,2) * t316) * t309) + t314 * (mrSges(5,2) * t266 - mrSges(5,3) * t263 + t321 * t248 - t318 * t247 - pkin(7) * t338 + (Ifges(5,1) * t314 + Ifges(5,4) * t316) * t309) + qJ(4) * t333 - pkin(3) * t252 + mrSges(4,1) * t269 + Ifges(4,3) * t309;
t1 = [pkin(1) * (t315 * (m(3) * t285 - mrSges(3,1) * t324 - qJDD(1) * mrSges(3,2) + t242 * t322 - t250 * t319) + t317 * (m(3) * t284 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t324 + t339)) + mrSges(2,1) * t336 - mrSges(2,2) * t334 + pkin(2) * t339 + mrSges(3,1) * t284 - mrSges(3,2) * t285 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t327; t316 * t245 + t314 * t246 + (m(3) + m(4)) * t313; t327; t252; mrSges(6,1) * t258 - mrSges(6,2) * t259 + Ifges(6,5) * t278 + Ifges(6,6) * t277 + Ifges(6,3) * qJDD(5) + t272 * t289 - t273 * t288;];
tauJ = t1;
