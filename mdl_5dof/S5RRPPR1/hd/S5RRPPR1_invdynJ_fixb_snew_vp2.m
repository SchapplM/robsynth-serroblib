% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:20
% EndTime: 2019-12-05 18:18:20
% DurationCPUTime: 0.63s
% Computational Cost: add. (5699->135), mult. (7766->178), div. (0->0), fcn. (4499->10), ass. (0->72)
t318 = cos(pkin(9));
t313 = t318 ^ 2;
t314 = qJD(1) + qJD(2);
t310 = t314 ^ 2;
t345 = pkin(4) * t318;
t316 = sin(pkin(9));
t344 = mrSges(5,2) * t316;
t343 = t310 * t313;
t311 = qJDD(1) + qJDD(2);
t342 = t311 * t318;
t322 = sin(qJ(1));
t325 = cos(qJ(1));
t337 = t325 * g(2) + t322 * g(3);
t296 = qJDD(1) * pkin(1) + t337;
t335 = t322 * g(2) - t325 * g(3);
t297 = -qJD(1) ^ 2 * pkin(1) + t335;
t321 = sin(qJ(2));
t324 = cos(qJ(2));
t285 = t324 * t296 - t321 * t297;
t282 = t311 * pkin(2) + t285;
t286 = t321 * t296 + t324 * t297;
t283 = -t310 * pkin(2) + t286;
t317 = sin(pkin(8));
t319 = cos(pkin(8));
t270 = t317 * t282 + t319 * t283;
t267 = -t310 * pkin(3) + t311 * qJ(4) + t270;
t315 = -g(1) + qJDD(3);
t336 = qJD(4) * t314;
t338 = t318 * t315 - 0.2e1 * t316 * t336;
t263 = -t316 * t267 + t338;
t331 = mrSges(5,3) * t311 + (-mrSges(5,1) * t318 + t344) * t310;
t260 = (-pkin(7) * t311 + t310 * t345 - t267) * t316 + t338;
t264 = t316 * t315 + (t267 + 0.2e1 * t336) * t318;
t261 = -pkin(4) * t343 + pkin(7) * t342 + t264;
t320 = sin(qJ(5));
t323 = cos(qJ(5));
t258 = t323 * t260 - t320 * t261;
t329 = -t316 * t320 + t318 * t323;
t289 = t329 * t314;
t330 = t316 * t323 + t318 * t320;
t290 = t330 * t314;
t276 = -t289 * mrSges(6,1) + t290 * mrSges(6,2);
t278 = t289 * qJD(5) + t330 * t311;
t287 = -qJD(5) * mrSges(6,2) + t289 * mrSges(6,3);
t256 = m(6) * t258 + qJDD(5) * mrSges(6,1) - t278 * mrSges(6,3) + qJD(5) * t287 - t290 * t276;
t259 = t320 * t260 + t323 * t261;
t277 = -t290 * qJD(5) + t329 * t311;
t288 = qJD(5) * mrSges(6,1) - t290 * mrSges(6,3);
t257 = m(6) * t259 - qJDD(5) * mrSges(6,2) + t277 * mrSges(6,3) - qJD(5) * t288 + t289 * t276;
t339 = t323 * t256 + t320 * t257;
t245 = m(5) * t263 - t331 * t316 + t339;
t333 = -t320 * t256 + t323 * t257;
t246 = m(5) * t264 + t331 * t318 + t333;
t334 = -t316 * t245 + t318 * t246;
t242 = m(4) * t270 - t310 * mrSges(4,1) - t311 * mrSges(4,2) + t334;
t269 = t319 * t282 - t317 * t283;
t332 = qJDD(4) - t269;
t266 = -t311 * pkin(3) - t310 * qJ(4) + t332;
t312 = t316 ^ 2;
t262 = (-pkin(3) - t345) * t311 + (-qJ(4) + (-t312 - t313) * pkin(7)) * t310 + t332;
t328 = m(6) * t262 - t277 * mrSges(6,1) + t278 * mrSges(6,2) - t289 * t287 + t290 * t288;
t326 = -m(5) * t266 + mrSges(5,1) * t342 - t328 + (t310 * t312 + t343) * mrSges(5,3);
t250 = m(4) * t269 - t310 * mrSges(4,2) + (mrSges(4,1) - t344) * t311 + t326;
t340 = t317 * t242 + t319 * t250;
t271 = Ifges(6,5) * t290 + Ifges(6,6) * t289 + Ifges(6,3) * qJD(5);
t273 = Ifges(6,1) * t290 + Ifges(6,4) * t289 + Ifges(6,5) * qJD(5);
t247 = -mrSges(6,1) * t262 + mrSges(6,3) * t259 + Ifges(6,4) * t278 + Ifges(6,2) * t277 + Ifges(6,6) * qJDD(5) + qJD(5) * t273 - t290 * t271;
t272 = Ifges(6,4) * t290 + Ifges(6,2) * t289 + Ifges(6,6) * qJD(5);
t248 = mrSges(6,2) * t262 - mrSges(6,3) * t258 + Ifges(6,1) * t278 + Ifges(6,4) * t277 + Ifges(6,5) * qJDD(5) - qJD(5) * t272 + t289 * t271;
t252 = t311 * t344 - t326;
t327 = -mrSges(3,2) * t286 - mrSges(4,2) * t270 + t318 * (-mrSges(5,1) * t266 + mrSges(5,3) * t264 - pkin(4) * t328 + pkin(7) * t333 + t323 * t247 + t320 * t248) + pkin(2) * t340 + t316 * (mrSges(5,2) * t266 - mrSges(5,3) * t263 - pkin(7) * t339 - t320 * t247 + t323 * t248) + qJ(4) * t334 - pkin(3) * t252 + mrSges(4,1) * t269 + mrSges(3,1) * t285 + (Ifges(5,2) * t313 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t316 + 0.2e1 * Ifges(5,4) * t318) * t316) * t311;
t1 = [t327 + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t337 + pkin(1) * (t321 * (m(3) * t286 - t310 * mrSges(3,1) - t311 * mrSges(3,2) + t319 * t242 - t317 * t250) + t324 * (m(3) * t285 + t311 * mrSges(3,1) - t310 * mrSges(3,2) + t340)) - mrSges(2,2) * t335; t327; m(4) * t315 + t318 * t245 + t316 * t246; t252; mrSges(6,1) * t258 - mrSges(6,2) * t259 + Ifges(6,5) * t278 + Ifges(6,6) * t277 + Ifges(6,3) * qJDD(5) + t290 * t272 - t289 * t273;];
tauJ = t1;
