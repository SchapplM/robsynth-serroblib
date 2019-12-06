% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:08
% EndTime: 2019-12-05 15:19:09
% DurationCPUTime: 0.78s
% Computational Cost: add. (4877->155), mult. (8653->210), div. (0->0), fcn. (6737->14), ass. (0->80)
t296 = sin(pkin(10));
t300 = cos(pkin(10));
t290 = -t300 * g(1) - t296 * g(2);
t295 = sin(pkin(11));
t299 = cos(pkin(11));
t289 = t296 * g(1) - t300 * g(2);
t294 = -g(3) + qJDD(1);
t298 = sin(pkin(5));
t302 = cos(pkin(5));
t314 = t289 * t302 + t294 * t298;
t263 = -t295 * t290 + t299 * t314;
t276 = -t298 * t289 + t302 * t294 + qJDD(2);
t297 = sin(pkin(6));
t301 = cos(pkin(6));
t327 = t263 * t301 + t276 * t297;
t264 = t299 * t290 + t295 * t314;
t305 = sin(qJ(3));
t308 = cos(qJ(3));
t256 = -t305 * t264 + t308 * t327;
t259 = -t297 * t263 + t301 * t276;
t307 = cos(qJ(4));
t324 = t307 * t259;
t257 = t308 * t264 + t327 * t305;
t310 = qJD(3) ^ 2;
t255 = -t310 * pkin(3) + qJDD(3) * pkin(8) + t257;
t304 = sin(qJ(4));
t251 = t307 * t255 + t304 * t259;
t323 = qJD(3) * t304;
t322 = t307 * qJD(3);
t321 = qJD(3) * qJD(4);
t320 = t304 * t321;
t319 = t307 * t321;
t285 = (-t307 * mrSges(5,1) + t304 * mrSges(5,2)) * qJD(3);
t288 = t307 * qJDD(3) - t320;
t291 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t323;
t286 = (-t307 * pkin(4) - t304 * pkin(9)) * qJD(3);
t309 = qJD(4) ^ 2;
t249 = -t309 * pkin(4) + qJDD(4) * pkin(9) + t286 * t322 + t251;
t254 = -qJDD(3) * pkin(3) - t310 * pkin(8) - t256;
t287 = t304 * qJDD(3) + t319;
t252 = (-t287 - t319) * pkin(9) + (-t288 + t320) * pkin(4) + t254;
t303 = sin(qJ(5));
t306 = cos(qJ(5));
t246 = -t303 * t249 + t306 * t252;
t283 = t306 * qJD(4) - t303 * t323;
t271 = t283 * qJD(5) + t303 * qJDD(4) + t306 * t287;
t284 = t303 * qJD(4) + t306 * t323;
t272 = -t283 * mrSges(6,1) + t284 * mrSges(6,2);
t293 = qJD(5) - t322;
t274 = -t293 * mrSges(6,2) + t283 * mrSges(6,3);
t282 = qJDD(5) - t288;
t244 = m(6) * t246 + t282 * mrSges(6,1) - t271 * mrSges(6,3) - t284 * t272 + t293 * t274;
t247 = t306 * t249 + t303 * t252;
t270 = -t284 * qJD(5) + t306 * qJDD(4) - t303 * t287;
t275 = t293 * mrSges(6,1) - t284 * mrSges(6,3);
t245 = m(6) * t247 - t282 * mrSges(6,2) + t270 * mrSges(6,3) + t283 * t272 - t293 * t275;
t317 = -t303 * t244 + t306 * t245;
t238 = m(5) * t251 - qJDD(4) * mrSges(5,2) + t288 * mrSges(5,3) - qJD(4) * t291 + t285 * t322 + t317;
t250 = -t304 * t255 + t324;
t292 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t322;
t248 = -qJDD(4) * pkin(4) - t309 * pkin(9) - t324 + (qJD(3) * t286 + t255) * t304;
t313 = -m(6) * t248 + t270 * mrSges(6,1) - t271 * mrSges(6,2) + t283 * t274 - t284 * t275;
t242 = m(5) * t250 + qJDD(4) * mrSges(5,1) - t287 * mrSges(5,3) + qJD(4) * t292 - t285 * t323 + t313;
t318 = t307 * t238 - t304 * t242;
t234 = m(4) * t257 - t310 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t318;
t239 = t306 * t244 + t303 * t245;
t312 = -m(5) * t254 + t288 * mrSges(5,1) - t287 * mrSges(5,2) - t291 * t323 + t292 * t322 - t239;
t236 = m(4) * t256 + qJDD(3) * mrSges(4,1) - t310 * mrSges(4,2) + t312;
t316 = t234 * t305 + t236 * t308;
t266 = Ifges(6,4) * t284 + Ifges(6,2) * t283 + Ifges(6,6) * t293;
t267 = Ifges(6,1) * t284 + Ifges(6,4) * t283 + Ifges(6,5) * t293;
t311 = mrSges(6,1) * t246 - mrSges(6,2) * t247 + Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t282 + t284 * t266 - t283 * t267;
t279 = Ifges(5,5) * qJD(4) + (t304 * Ifges(5,1) + t307 * Ifges(5,4)) * qJD(3);
t278 = Ifges(5,6) * qJD(4) + (t304 * Ifges(5,4) + t307 * Ifges(5,2)) * qJD(3);
t265 = Ifges(6,5) * t284 + Ifges(6,6) * t283 + Ifges(6,3) * t293;
t241 = mrSges(6,2) * t248 - mrSges(6,3) * t246 + Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t282 + t283 * t265 - t293 * t266;
t240 = -mrSges(6,1) * t248 + mrSges(6,3) * t247 + Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t282 - t284 * t265 + t293 * t267;
t235 = m(4) * t259 + t304 * t238 + t307 * t242;
t233 = m(3) * t276 + t301 * t235 + t297 * t316;
t1 = [m(2) * t294 + t302 * t233 + (t295 * (m(3) * t264 + t308 * t234 - t305 * t236) + t299 * (m(3) * t263 - t297 * t235 + t301 * t316)) * t298; t233; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t256 - mrSges(4,2) * t257 + t304 * (mrSges(5,2) * t254 - mrSges(5,3) * t250 + Ifges(5,1) * t287 + Ifges(5,4) * t288 + Ifges(5,5) * qJDD(4) - pkin(9) * t239 - qJD(4) * t278 - t303 * t240 + t306 * t241) + t307 * (-mrSges(5,1) * t254 + mrSges(5,3) * t251 + Ifges(5,4) * t287 + Ifges(5,2) * t288 + Ifges(5,6) * qJDD(4) - pkin(4) * t239 + qJD(4) * t279 - t311) + pkin(3) * t312 + pkin(8) * t318; Ifges(5,5) * t287 + Ifges(5,6) * t288 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t250 - mrSges(5,2) * t251 + t303 * t241 + t306 * t240 + pkin(4) * t313 + pkin(9) * t317 + (t304 * t278 - t307 * t279) * qJD(3); t311;];
tauJ = t1;
