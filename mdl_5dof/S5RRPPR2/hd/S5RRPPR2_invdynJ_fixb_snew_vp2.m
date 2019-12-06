% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPPR2
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:15
% EndTime: 2019-12-05 18:20:15
% DurationCPUTime: 0.69s
% Computational Cost: add. (4951->131), mult. (6524->179), div. (0->0), fcn. (3631->10), ass. (0->73)
t333 = 2 * qJD(4);
t307 = sin(qJ(1));
t310 = cos(qJ(1));
t323 = t310 * g(2) + t307 * g(3);
t283 = qJDD(1) * pkin(1) + t323;
t318 = t307 * g(2) - g(3) * t310;
t284 = -qJD(1) ^ 2 * pkin(1) + t318;
t306 = sin(qJ(2));
t309 = cos(qJ(2));
t270 = t309 * t283 - t284 * t306;
t298 = qJDD(1) + qJDD(2);
t264 = pkin(2) * t298 + t270;
t271 = t306 * t283 + t309 * t284;
t299 = qJD(1) + qJD(2);
t297 = t299 ^ 2;
t265 = -pkin(2) * t297 + t271;
t302 = sin(pkin(8));
t304 = cos(pkin(8));
t260 = t302 * t264 + t304 * t265;
t257 = -pkin(3) * t297 + qJ(4) * t298 + t260;
t332 = t299 * t333 + t257;
t301 = sin(pkin(9));
t331 = mrSges(5,2) * t301;
t329 = mrSges(5,3) * t298;
t303 = cos(pkin(9));
t328 = t298 * t303;
t327 = t299 * t301;
t326 = t299 * t303;
t300 = -g(1) + qJDD(3);
t325 = t300 * t303;
t253 = t301 * t300 + t332 * t303;
t316 = -pkin(4) * t303 - pkin(7) * t301;
t279 = t316 * t299;
t251 = t279 * t326 + t253;
t259 = t264 * t304 - t302 * t265;
t314 = -qJ(4) * t297 + qJDD(4) - t259;
t254 = (-pkin(3) + t316) * t298 + t314;
t305 = sin(qJ(5));
t308 = cos(qJ(5));
t248 = -t251 * t305 + t254 * t308;
t286 = qJD(5) - t326;
t320 = t305 * t327;
t272 = -mrSges(6,2) * t286 - mrSges(6,3) * t320;
t274 = (t305 * mrSges(6,1) + t308 * mrSges(6,2)) * t327;
t321 = qJD(5) * t299;
t276 = (t298 * t308 - t305 * t321) * t301;
t285 = qJDD(5) - t328;
t319 = t308 * t327;
t246 = m(6) * t248 + mrSges(6,1) * t285 - mrSges(6,3) * t276 + t272 * t286 - t274 * t319;
t249 = t251 * t308 + t254 * t305;
t273 = mrSges(6,1) * t286 - mrSges(6,3) * t319;
t275 = (-t298 * t305 - t308 * t321) * t301;
t247 = m(6) * t249 - mrSges(6,2) * t285 + mrSges(6,3) * t275 - t273 * t286 - t274 * t320;
t277 = (-mrSges(5,1) * t303 + t331) * t299;
t241 = m(5) * t253 - t246 * t305 + t247 * t308 + (t277 * t299 + t329) * t303;
t250 = -t325 + (t257 + (t333 + t279) * t299) * t301;
t252 = -t332 * t301 + t325;
t245 = m(5) * t252 - m(6) * t250 + mrSges(6,1) * t275 - mrSges(6,2) * t276 + (-t329 + (-t272 * t305 - t273 * t308 - t277) * t299) * t301;
t317 = t303 * t241 - t245 * t301;
t236 = m(4) * t260 - mrSges(4,1) * t297 - mrSges(4,2) * t298 + t317;
t244 = t246 * t308 + t247 * t305;
t256 = -pkin(3) * t298 + t314;
t313 = -m(5) * t256 + mrSges(5,1) * t328 - t244 + (t301 ^ 2 + t303 ^ 2) * mrSges(5,3) * t297;
t239 = m(4) * t259 - mrSges(4,2) * t297 + (mrSges(4,1) - t331) * t298 + t313;
t324 = t302 * t236 + t304 * t239;
t267 = Ifges(6,6) * t286 + (t308 * Ifges(6,4) - t305 * Ifges(6,2)) * t327;
t268 = Ifges(6,5) * t286 + (Ifges(6,1) * t308 - Ifges(6,4) * t305) * t327;
t315 = t308 * t267 + t305 * t268;
t312 = mrSges(6,1) * t248 - mrSges(6,2) * t249 + Ifges(6,5) * t276 + Ifges(6,6) * t275 + Ifges(6,3) * t285;
t243 = t298 * t331 - t313;
t278 = (Ifges(5,5) * t301 + Ifges(5,6) * t303) * t299;
t311 = -mrSges(3,2) * t271 - mrSges(4,2) * t260 + pkin(2) * t324 + t301 * (t278 * t326 + mrSges(5,2) * t256 - mrSges(5,3) * t252 + t308 * (mrSges(6,2) * t250 - mrSges(6,3) * t248 + Ifges(6,1) * t276 + Ifges(6,4) * t275 + Ifges(6,5) * t285 - t267 * t286) - t305 * (-mrSges(6,1) * t250 + mrSges(6,3) * t249 + Ifges(6,4) * t276 + Ifges(6,2) * t275 + Ifges(6,6) * t285 + t268 * t286) - pkin(7) * t244 + (Ifges(5,1) * t301 + Ifges(5,4) * t303) * t298) + t303 * (Ifges(5,2) * t328 - mrSges(5,1) * t256 + mrSges(5,3) * t253 - pkin(4) * t244 + (Ifges(5,4) * t298 + (-t278 - t315) * t299) * t301 - t312) + qJ(4) * t317 - pkin(3) * t243 + mrSges(4,1) * t259 + mrSges(3,1) * t270 + Ifges(4,3) * t298 + Ifges(3,3) * t298;
t1 = [t311 + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t323 + pkin(1) * (t306 * (m(3) * t271 - mrSges(3,1) * t297 - mrSges(3,2) * t298 + t236 * t304 - t239 * t302) + t309 * (m(3) * t270 + mrSges(3,1) * t298 - mrSges(3,2) * t297 + t324)) - mrSges(2,2) * t318; t311; m(4) * t300 + t241 * t301 + t245 * t303; t243; t315 * t327 + t312;];
tauJ = t1;
