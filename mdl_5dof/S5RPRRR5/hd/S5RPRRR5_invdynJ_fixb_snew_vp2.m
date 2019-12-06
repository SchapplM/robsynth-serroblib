% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:20
% EndTime: 2019-12-05 18:16:21
% DurationCPUTime: 0.60s
% Computational Cost: add. (6034->168), mult. (8145->221), div. (0->0), fcn. (4613->10), ass. (0->78)
t316 = qJD(1) + qJD(3);
t322 = sin(qJ(4));
t343 = t316 * t322;
t326 = cos(qJ(4));
t342 = t316 * t326;
t324 = sin(qJ(1));
t328 = cos(qJ(1));
t340 = t328 * g(2) + t324 * g(3);
t300 = qJDD(1) * pkin(1) + t340;
t329 = qJD(1) ^ 2;
t337 = t324 * g(2) - t328 * g(3);
t301 = -t329 * pkin(1) + t337;
t319 = sin(pkin(9));
t320 = cos(pkin(9));
t283 = t320 * t300 - t319 * t301;
t281 = qJDD(1) * pkin(2) + t283;
t284 = t319 * t300 + t320 * t301;
t282 = -t329 * pkin(2) + t284;
t323 = sin(qJ(3));
t327 = cos(qJ(3));
t266 = t323 * t281 + t327 * t282;
t312 = t316 ^ 2;
t314 = qJDD(1) + qJDD(3);
t263 = -t312 * pkin(3) + t314 * pkin(7) + t266;
t318 = -g(1) + qJDD(2);
t259 = -t322 * t263 + t326 * t318;
t339 = qJD(4) * t316;
t338 = t326 * t339;
t295 = t322 * t314 + t338;
t256 = (-t295 + t338) * pkin(8) + (t312 * t322 * t326 + qJDD(4)) * pkin(4) + t259;
t260 = t326 * t263 + t322 * t318;
t296 = t326 * t314 - t322 * t339;
t304 = qJD(4) * pkin(4) - pkin(8) * t343;
t317 = t326 ^ 2;
t257 = -t317 * t312 * pkin(4) + t296 * pkin(8) - qJD(4) * t304 + t260;
t321 = sin(qJ(5));
t325 = cos(qJ(5));
t254 = t325 * t256 - t321 * t257;
t290 = (-t321 * t322 + t325 * t326) * t316;
t272 = t290 * qJD(5) + t325 * t295 + t321 * t296;
t291 = (t321 * t326 + t322 * t325) * t316;
t277 = -t290 * mrSges(6,1) + t291 * mrSges(6,2);
t315 = qJD(4) + qJD(5);
t285 = -t315 * mrSges(6,2) + t290 * mrSges(6,3);
t313 = qJDD(4) + qJDD(5);
t251 = m(6) * t254 + t313 * mrSges(6,1) - t272 * mrSges(6,3) - t291 * t277 + t315 * t285;
t255 = t321 * t256 + t325 * t257;
t271 = -t291 * qJD(5) - t321 * t295 + t325 * t296;
t286 = t315 * mrSges(6,1) - t291 * mrSges(6,3);
t252 = m(6) * t255 - t313 * mrSges(6,2) + t271 * mrSges(6,3) + t290 * t277 - t315 * t286;
t242 = t325 * t251 + t321 * t252;
t294 = (-mrSges(5,1) * t326 + mrSges(5,2) * t322) * t316;
t303 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t342;
t240 = m(5) * t259 + qJDD(4) * mrSges(5,1) - t295 * mrSges(5,3) + qJD(4) * t303 - t294 * t343 + t242;
t302 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t343;
t335 = -t321 * t251 + t325 * t252;
t241 = m(5) * t260 - qJDD(4) * mrSges(5,2) + t296 * mrSges(5,3) - qJD(4) * t302 + t294 * t342 + t335;
t336 = -t322 * t240 + t326 * t241;
t237 = m(4) * t266 - t312 * mrSges(4,1) - t314 * mrSges(4,2) + t336;
t265 = t327 * t281 - t323 * t282;
t334 = -t314 * pkin(3) - t265;
t262 = -t312 * pkin(7) + t334;
t258 = t304 * t343 - t296 * pkin(4) + (-pkin(8) * t317 - pkin(7)) * t312 + t334;
t332 = m(6) * t258 - t271 * mrSges(6,1) + t272 * mrSges(6,2) - t290 * t285 + t291 * t286;
t330 = -m(5) * t262 + t296 * mrSges(5,1) - t295 * mrSges(5,2) - t302 * t343 + t303 * t342 - t332;
t246 = m(4) * t265 + t314 * mrSges(4,1) - t312 * mrSges(4,2) + t330;
t341 = t323 * t237 + t327 * t246;
t273 = Ifges(6,5) * t291 + Ifges(6,6) * t290 + Ifges(6,3) * t315;
t275 = Ifges(6,1) * t291 + Ifges(6,4) * t290 + Ifges(6,5) * t315;
t243 = -mrSges(6,1) * t258 + mrSges(6,3) * t255 + Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t313 - t291 * t273 + t315 * t275;
t274 = Ifges(6,4) * t291 + Ifges(6,2) * t290 + Ifges(6,6) * t315;
t244 = mrSges(6,2) * t258 - mrSges(6,3) * t254 + Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t313 + t290 * t273 - t315 * t274;
t287 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t322 + Ifges(5,6) * t326) * t316;
t288 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t322 + Ifges(5,2) * t326) * t316;
t289 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t322 + Ifges(5,4) * t326) * t316;
t333 = -mrSges(4,2) * t266 + t326 * (-mrSges(5,1) * t262 + mrSges(5,3) * t260 + Ifges(5,4) * t295 + Ifges(5,2) * t296 + Ifges(5,6) * qJDD(4) - pkin(4) * t332 + pkin(8) * t335 + qJD(4) * t289 + t325 * t243 + t321 * t244 - t287 * t343) + t322 * (mrSges(5,2) * t262 - mrSges(5,3) * t259 + Ifges(5,1) * t295 + Ifges(5,4) * t296 + Ifges(5,5) * qJDD(4) - pkin(8) * t242 - qJD(4) * t288 - t321 * t243 + t325 * t244 + t287 * t342) + pkin(7) * t336 + pkin(3) * t330 + mrSges(4,1) * t265 + Ifges(4,3) * t314;
t331 = mrSges(6,1) * t254 - mrSges(6,2) * t255 + Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t313 + t291 * t274 - t290 * t275;
t1 = [pkin(1) * (t319 * (m(3) * t284 - t329 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t327 * t237 - t323 * t246) + t320 * (m(3) * t283 + qJDD(1) * mrSges(3,1) - t329 * mrSges(3,2) + t341)) + mrSges(2,1) * t340 - mrSges(2,2) * t337 + pkin(2) * t341 + mrSges(3,1) * t283 - mrSges(3,2) * t284 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t333; t326 * t240 + t322 * t241 + (m(3) + m(4)) * t318; t333; mrSges(5,1) * t259 - mrSges(5,2) * t260 + Ifges(5,5) * t295 + Ifges(5,6) * t296 + Ifges(5,3) * qJDD(4) + pkin(4) * t242 + (t288 * t322 - t289 * t326) * t316 + t331; t331;];
tauJ = t1;
