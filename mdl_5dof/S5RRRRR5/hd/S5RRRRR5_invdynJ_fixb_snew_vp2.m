% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:22
% EndTime: 2019-12-05 18:58:23
% DurationCPUTime: 0.70s
% Computational Cost: add. (10318->167), mult. (10575->220), div. (0->0), fcn. (5990->10), ass. (0->79)
t317 = qJD(1) + qJD(2);
t309 = qJD(3) + t317;
t307 = t309 ^ 2;
t344 = pkin(4) * t307;
t320 = sin(qJ(4));
t343 = t309 * t320;
t325 = cos(qJ(4));
t342 = t309 * t325;
t323 = sin(qJ(1));
t328 = cos(qJ(1));
t339 = t328 * g(2) + t323 * g(3);
t302 = qJDD(1) * pkin(1) + t339;
t337 = t323 * g(2) - t328 * g(3);
t303 = -qJD(1) ^ 2 * pkin(1) + t337;
t322 = sin(qJ(2));
t327 = cos(qJ(2));
t282 = t327 * t302 - t322 * t303;
t315 = qJDD(1) + qJDD(2);
t279 = t315 * pkin(2) + t282;
t283 = t322 * t302 + t327 * t303;
t313 = t317 ^ 2;
t280 = -t313 * pkin(2) + t283;
t321 = sin(qJ(3));
t326 = cos(qJ(3));
t264 = t321 * t279 + t326 * t280;
t308 = qJDD(3) + t315;
t261 = -t307 * pkin(3) + t308 * pkin(8) + t264;
t341 = t320 * t261;
t338 = qJD(4) * t309;
t294 = t320 * t308 + t325 * t338;
t254 = qJDD(4) * pkin(4) - t294 * pkin(9) - t341 + (pkin(9) * t338 + t320 * t344 - g(1)) * t325;
t258 = -t320 * g(1) + t325 * t261;
t295 = t325 * t308 - t320 * t338;
t301 = qJD(4) * pkin(4) - pkin(9) * t343;
t318 = t325 ^ 2;
t255 = t295 * pkin(9) - qJD(4) * t301 - t318 * t344 + t258;
t319 = sin(qJ(5));
t324 = cos(qJ(5));
t252 = t324 * t254 - t319 * t255;
t289 = (-t319 * t320 + t324 * t325) * t309;
t270 = t289 * qJD(5) + t324 * t294 + t319 * t295;
t290 = (t319 * t325 + t320 * t324) * t309;
t275 = -t289 * mrSges(6,1) + t290 * mrSges(6,2);
t316 = qJD(4) + qJD(5);
t284 = -t316 * mrSges(6,2) + t289 * mrSges(6,3);
t314 = qJDD(4) + qJDD(5);
t249 = m(6) * t252 + t314 * mrSges(6,1) - t270 * mrSges(6,3) - t290 * t275 + t316 * t284;
t253 = t319 * t254 + t324 * t255;
t269 = -t290 * qJD(5) - t319 * t294 + t324 * t295;
t285 = t316 * mrSges(6,1) - t290 * mrSges(6,3);
t250 = m(6) * t253 - t314 * mrSges(6,2) + t269 * mrSges(6,3) + t289 * t275 - t316 * t285;
t240 = t324 * t249 + t319 * t250;
t257 = -t325 * g(1) - t341;
t293 = (-mrSges(5,1) * t325 + mrSges(5,2) * t320) * t309;
t299 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t343;
t300 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t342;
t335 = -t319 * t249 + t324 * t250;
t336 = -t320 * (m(5) * t257 + qJDD(4) * mrSges(5,1) - t294 * mrSges(5,3) + qJD(4) * t300 - t293 * t343 + t240) + t325 * (m(5) * t258 - qJDD(4) * mrSges(5,2) + t295 * mrSges(5,3) - qJD(4) * t299 + t293 * t342 + t335);
t236 = m(4) * t264 - t307 * mrSges(4,1) - t308 * mrSges(4,2) + t336;
t263 = t326 * t279 - t321 * t280;
t334 = -t308 * pkin(3) - t263;
t260 = -t307 * pkin(8) + t334;
t256 = t301 * t343 - t295 * pkin(4) + (-pkin(9) * t318 - pkin(8)) * t307 + t334;
t332 = m(6) * t256 - t269 * mrSges(6,1) + t270 * mrSges(6,2) - t289 * t284 + t290 * t285;
t329 = -m(5) * t260 + t295 * mrSges(5,1) - t294 * mrSges(5,2) - t299 * t343 + t300 * t342 - t332;
t244 = m(4) * t263 + t308 * mrSges(4,1) - t307 * mrSges(4,2) + t329;
t340 = t321 * t236 + t326 * t244;
t271 = Ifges(6,5) * t290 + Ifges(6,6) * t289 + Ifges(6,3) * t316;
t273 = Ifges(6,1) * t290 + Ifges(6,4) * t289 + Ifges(6,5) * t316;
t241 = -mrSges(6,1) * t256 + mrSges(6,3) * t253 + Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t314 - t290 * t271 + t316 * t273;
t272 = Ifges(6,4) * t290 + Ifges(6,2) * t289 + Ifges(6,6) * t316;
t242 = mrSges(6,2) * t256 - mrSges(6,3) * t252 + Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t314 + t289 * t271 - t316 * t272;
t286 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t320 + Ifges(5,6) * t325) * t309;
t287 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t320 + Ifges(5,2) * t325) * t309;
t288 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t320 + Ifges(5,4) * t325) * t309;
t333 = -mrSges(4,2) * t264 + t325 * (-mrSges(5,1) * t260 + mrSges(5,3) * t258 + Ifges(5,4) * t294 + Ifges(5,2) * t295 + Ifges(5,6) * qJDD(4) - pkin(4) * t332 + pkin(9) * t335 + qJD(4) * t288 + t324 * t241 + t319 * t242 - t286 * t343) + t320 * (mrSges(5,2) * t260 - mrSges(5,3) * t257 + Ifges(5,1) * t294 + Ifges(5,4) * t295 + Ifges(5,5) * qJDD(4) - pkin(9) * t240 - qJD(4) * t287 - t319 * t241 + t324 * t242 + t286 * t342) + pkin(8) * t336 + pkin(3) * t329 + mrSges(4,1) * t263 + Ifges(4,3) * t308;
t331 = mrSges(6,1) * t252 - mrSges(6,2) * t253 + Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t314 + t290 * t272 - t289 * t273;
t330 = mrSges(3,1) * t282 - mrSges(3,2) * t283 + Ifges(3,3) * t315 + pkin(2) * t340 + t333;
t1 = [t330 + Ifges(2,3) * qJDD(1) + pkin(1) * (t322 * (m(3) * t283 - t313 * mrSges(3,1) - t315 * mrSges(3,2) + t326 * t236 - t321 * t244) + t327 * (m(3) * t282 + t315 * mrSges(3,1) - t313 * mrSges(3,2) + t340)) - mrSges(2,2) * t337 + mrSges(2,1) * t339; t330; t333; mrSges(5,1) * t257 - mrSges(5,2) * t258 + Ifges(5,5) * t294 + Ifges(5,6) * t295 + Ifges(5,3) * qJDD(4) + pkin(4) * t240 + (t287 * t320 - t288 * t325) * t309 + t331; t331;];
tauJ = t1;
