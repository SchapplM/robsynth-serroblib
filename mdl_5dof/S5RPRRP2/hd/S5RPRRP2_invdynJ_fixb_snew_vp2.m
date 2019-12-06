% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:32
% EndTime: 2019-12-05 18:01:33
% DurationCPUTime: 0.54s
% Computational Cost: add. (2571->147), mult. (3501->178), div. (0->0), fcn. (1699->8), ass. (0->72)
t350 = Ifges(5,1) + Ifges(6,1);
t345 = Ifges(5,4) + Ifges(6,4);
t344 = Ifges(5,5) + Ifges(6,5);
t349 = Ifges(5,2) + Ifges(6,2);
t343 = Ifges(5,6) + Ifges(6,6);
t348 = Ifges(5,3) + Ifges(6,3);
t311 = qJD(1) + qJD(3);
t309 = t311 ^ 2;
t347 = pkin(4) * t309;
t346 = -mrSges(5,2) - mrSges(6,2);
t317 = sin(qJ(4));
t342 = t311 * t317;
t320 = cos(qJ(4));
t341 = t311 * t320;
t319 = sin(qJ(1));
t322 = cos(qJ(1));
t335 = t322 * g(2) + t319 * g(3);
t296 = qJDD(1) * pkin(1) + t335;
t323 = qJD(1) ^ 2;
t329 = t319 * g(2) - t322 * g(3);
t297 = -t323 * pkin(1) + t329;
t315 = sin(pkin(8));
t316 = cos(pkin(8));
t274 = t316 * t296 - t297 * t315;
t272 = qJDD(1) * pkin(2) + t274;
t275 = t315 * t296 + t316 * t297;
t273 = -t323 * pkin(2) + t275;
t318 = sin(qJ(3));
t321 = cos(qJ(3));
t268 = t318 * t272 + t321 * t273;
t310 = qJDD(1) + qJDD(3);
t265 = -pkin(3) * t309 + pkin(7) * t310 + t268;
t314 = -g(1) + qJDD(2);
t304 = t320 * t314;
t261 = -t317 * t265 + t304;
t288 = (-mrSges(6,1) * t320 + mrSges(6,2) * t317) * t311;
t289 = (-mrSges(5,1) * t320 + mrSges(5,2) * t317) * t311;
t334 = qJD(4) * t311;
t330 = t320 * t334;
t290 = t317 * t310 + t330;
t302 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t341;
t333 = qJD(5) * t311;
t258 = qJDD(4) * pkin(4) + t304 + (-t290 + t330) * qJ(5) + (t320 * t347 - t265 - 0.2e1 * t333) * t317;
t301 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t341;
t332 = m(6) * t258 + qJDD(4) * mrSges(6,1) + qJD(4) * t301;
t252 = m(5) * t261 + qJDD(4) * mrSges(5,1) + qJD(4) * t302 + (-t288 - t289) * t342 + (-mrSges(5,3) - mrSges(6,3)) * t290 + t332;
t262 = t320 * t265 + t317 * t314;
t291 = t320 * t310 - t317 * t334;
t298 = qJD(4) * pkin(4) - qJ(5) * t342;
t313 = t320 ^ 2;
t259 = t291 * qJ(5) - qJD(4) * t298 - t313 * t347 + 0.2e1 * t320 * t333 + t262;
t331 = m(6) * t259 + t291 * mrSges(6,3) + t288 * t341;
t299 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t342;
t336 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t342 - t299;
t253 = m(5) * t262 + t291 * mrSges(5,3) + t336 * qJD(4) + t346 * qJDD(4) + t289 * t341 + t331;
t328 = -t317 * t252 + t320 * t253;
t246 = m(4) * t268 - t309 * mrSges(4,1) - t310 * mrSges(4,2) + t328;
t267 = t321 * t272 - t318 * t273;
t326 = -t310 * pkin(3) - t267;
t264 = -t309 * pkin(7) + t326;
t260 = t298 * t342 - t291 * pkin(4) + qJDD(5) + (-qJ(5) * t313 - pkin(7)) * t309 + t326;
t327 = -m(6) * t260 + t291 * mrSges(6,1) + t301 * t341;
t324 = -m(5) * t264 + t291 * mrSges(5,1) + t346 * t290 + t302 * t341 + t336 * t342 + t327;
t249 = m(4) * t267 + t310 * mrSges(4,1) - t309 * mrSges(4,2) + t324;
t340 = t318 * t246 + t321 * t249;
t339 = (t317 * t344 + t343 * t320) * t311 + t348 * qJD(4);
t338 = (t317 * t345 + t349 * t320) * t311 + t343 * qJD(4);
t337 = (-t350 * t317 - t320 * t345) * t311 - t344 * qJD(4);
t254 = t290 * mrSges(6,2) + t299 * t342 - t327;
t255 = -t290 * mrSges(6,3) - t288 * t342 + t332;
t325 = -mrSges(4,2) * t268 + t320 * (-mrSges(5,1) * t264 + mrSges(5,3) * t262 - mrSges(6,1) * t260 + mrSges(6,3) * t259 - pkin(4) * t254 + qJ(5) * t331 - t339 * t342 + t349 * t291 + t345 * t290 + (-qJ(5) * mrSges(6,2) + t343) * qJDD(4) + (-qJ(5) * t299 - t337) * qJD(4)) + t317 * (mrSges(5,2) * t264 + mrSges(6,2) * t260 - mrSges(5,3) * t261 - mrSges(6,3) * t258 - qJ(5) * t255 - t338 * qJD(4) + t344 * qJDD(4) + t350 * t290 + t345 * t291 + t339 * t341) + pkin(7) * t328 + pkin(3) * t324 + mrSges(4,1) * t267 + Ifges(4,3) * t310;
t1 = [mrSges(2,1) * t335 - mrSges(2,2) * t329 + pkin(1) * (t315 * (m(3) * t275 - t323 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t321 * t246 - t318 * t249) + t316 * (m(3) * t274 + qJDD(1) * mrSges(3,1) - t323 * mrSges(3,2) + t340)) + pkin(2) * t340 + mrSges(3,1) * t274 - mrSges(3,2) * t275 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t325; t320 * t252 + t317 * t253 + (m(3) + m(4)) * t314; t325; mrSges(5,1) * t261 + mrSges(6,1) * t258 - mrSges(5,2) * t262 - mrSges(6,2) * t259 + pkin(4) * t255 + t343 * t291 + t344 * t290 + t348 * qJDD(4) + (t338 * t317 + t337 * t320) * t311; t254;];
tauJ = t1;
