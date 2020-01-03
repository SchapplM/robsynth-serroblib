% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:43
% EndTime: 2019-12-31 18:44:44
% DurationCPUTime: 1.04s
% Computational Cost: add. (3797->194), mult. (7114->237), div. (0->0), fcn. (4052->8), ass. (0->83)
t371 = Ifges(5,1) + Ifges(6,1);
t361 = Ifges(5,4) - Ifges(6,5);
t367 = -Ifges(5,5) - Ifges(6,4);
t370 = Ifges(5,2) + Ifges(6,3);
t359 = Ifges(5,6) - Ifges(6,6);
t334 = sin(qJ(4));
t335 = sin(qJ(3));
t354 = t335 * qJD(1);
t363 = cos(qJ(4));
t315 = -qJD(3) * t363 + t334 * t354;
t337 = cos(qJ(3));
t352 = qJD(1) * qJD(3);
t349 = t337 * t352;
t321 = qJDD(1) * t335 + t349;
t291 = -qJD(4) * t315 + qJDD(3) * t334 + t321 * t363;
t316 = qJD(3) * t334 + t354 * t363;
t297 = mrSges(6,1) * t315 - mrSges(6,3) * t316;
t336 = sin(qJ(1));
t338 = cos(qJ(1));
t348 = t336 * g(1) - g(2) * t338;
t317 = qJDD(1) * pkin(1) + t348;
t340 = qJD(1) ^ 2;
t344 = -g(1) * t338 - g(2) * t336;
t319 = -pkin(1) * t340 + t344;
t332 = sin(pkin(8));
t333 = cos(pkin(8));
t292 = t333 * t317 - t319 * t332;
t278 = -qJDD(1) * pkin(2) - t340 * pkin(6) - t292;
t350 = t335 * t352;
t322 = qJDD(1) * t337 - t350;
t270 = (-t321 - t349) * pkin(7) + (-t322 + t350) * pkin(3) + t278;
t293 = t317 * t332 + t319 * t333;
t279 = -pkin(2) * t340 + qJDD(1) * pkin(6) + t293;
t331 = -g(3) + qJDD(2);
t275 = t279 * t337 + t335 * t331;
t320 = (-pkin(3) * t337 - pkin(7) * t335) * qJD(1);
t339 = qJD(3) ^ 2;
t353 = t337 * qJD(1);
t273 = -pkin(3) * t339 + qJDD(3) * pkin(7) + t320 * t353 + t275;
t267 = t270 * t363 - t273 * t334;
t296 = pkin(4) * t315 - qJ(5) * t316;
t314 = qJDD(4) - t322;
t326 = qJD(4) - t353;
t325 = t326 ^ 2;
t265 = -pkin(4) * t314 - qJ(5) * t325 + t296 * t316 + qJDD(5) - t267;
t302 = -mrSges(6,2) * t315 + mrSges(6,3) * t326;
t345 = -m(6) * t265 + t314 * mrSges(6,1) + t302 * t326;
t261 = t291 * mrSges(6,2) + t316 * t297 - t345;
t268 = t270 * t334 + t273 * t363;
t264 = -pkin(4) * t325 + qJ(5) * t314 + 0.2e1 * qJD(5) * t326 - t296 * t315 + t268;
t290 = qJD(4) * t316 - qJDD(3) * t363 + t321 * t334;
t301 = -mrSges(6,1) * t326 + mrSges(6,2) * t316;
t351 = m(6) * t264 + t314 * mrSges(6,3) + t301 * t326;
t357 = t315 * t370 - t316 * t361 - t326 * t359;
t364 = t315 * t361 - t316 * t371 + t326 * t367;
t366 = -Ifges(5,3) - Ifges(6,2);
t369 = -t367 * t291 - t364 * t315 - t359 * t290 - t366 * t314 + mrSges(5,1) * t267 - mrSges(6,1) * t265 - mrSges(5,2) * t268 + mrSges(6,3) * t264 - pkin(4) * t261 + qJ(5) * (-mrSges(6,2) * t290 - t315 * t297 + t351) - t357 * t316;
t362 = -mrSges(5,3) - mrSges(6,2);
t358 = t315 * t359 + t316 * t367 + t326 * t366;
t355 = -mrSges(5,1) * t315 - mrSges(5,2) * t316 - t297;
t318 = (-mrSges(4,1) * t337 + mrSges(4,2) * t335) * qJD(1);
t323 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t354;
t300 = mrSges(5,1) * t326 - mrSges(5,3) * t316;
t258 = m(5) * t268 - t314 * mrSges(5,2) + t290 * t362 - t326 * t300 + t315 * t355 + t351;
t299 = -mrSges(5,2) * t326 - mrSges(5,3) * t315;
t259 = m(5) * t267 + t314 * mrSges(5,1) + t291 * t362 + t326 * t299 + t316 * t355 + t345;
t346 = t258 * t363 - t259 * t334;
t254 = m(4) * t275 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t322 - qJD(3) * t323 + t318 * t353 + t346;
t274 = -t279 * t335 + t337 * t331;
t324 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t353;
t272 = -qJDD(3) * pkin(3) - t339 * pkin(7) + t320 * t354 - t274;
t266 = -0.2e1 * qJD(5) * t316 + (t315 * t326 - t291) * qJ(5) + (t316 * t326 + t290) * pkin(4) + t272;
t262 = m(6) * t266 + mrSges(6,1) * t290 - mrSges(6,3) * t291 - t301 * t316 + t302 * t315;
t341 = -m(5) * t272 - mrSges(5,1) * t290 - mrSges(5,2) * t291 - t299 * t315 - t300 * t316 - t262;
t256 = m(4) * t274 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t321 + qJD(3) * t324 - t318 * t354 + t341;
t347 = t254 * t337 - t256 * t335;
t255 = t258 * t334 + t259 * t363;
t342 = -m(4) * t278 + t322 * mrSges(4,1) - mrSges(4,2) * t321 - t323 * t354 + t324 * t353 - t255;
t309 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t335 + Ifges(4,4) * t337) * qJD(1);
t308 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t335 + Ifges(4,2) * t337) * qJD(1);
t252 = mrSges(5,2) * t272 + mrSges(6,2) * t265 - mrSges(5,3) * t267 - mrSges(6,3) * t266 - qJ(5) * t262 - t361 * t290 + t291 * t371 - t367 * t314 + t358 * t315 + t357 * t326;
t251 = -mrSges(5,1) * t272 - mrSges(6,1) * t266 + mrSges(6,2) * t264 + mrSges(5,3) * t268 - pkin(4) * t262 - t290 * t370 + t361 * t291 + t359 * t314 + t358 * t316 - t364 * t326;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t348 - mrSges(2,2) * t344 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t292 - mrSges(3,2) * t293 + t335 * (mrSges(4,2) * t278 - mrSges(4,3) * t274 + Ifges(4,1) * t321 + Ifges(4,4) * t322 + Ifges(4,5) * qJDD(3) - pkin(7) * t255 - qJD(3) * t308 - t334 * t251 + t252 * t363) + t337 * (-mrSges(4,1) * t278 + mrSges(4,3) * t275 + Ifges(4,4) * t321 + Ifges(4,2) * t322 + Ifges(4,6) * qJDD(3) - pkin(3) * t255 + qJD(3) * t309 - t369) + pkin(2) * t342 + pkin(6) * t347 + pkin(1) * (t332 * (m(3) * t293 - mrSges(3,1) * t340 - qJDD(1) * mrSges(3,2) + t347) + t333 * (m(3) * t292 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t340 + t342)); m(3) * t331 + t254 * t335 + t256 * t337; Ifges(4,5) * t321 + Ifges(4,6) * t322 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t274 - mrSges(4,2) * t275 + t334 * t252 + t363 * t251 + pkin(3) * t341 + pkin(7) * t346 + (t308 * t335 - t309 * t337) * qJD(1); t369; t261;];
tauJ = t1;
