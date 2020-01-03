% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:43
% EndTime: 2019-12-31 17:16:44
% DurationCPUTime: 0.92s
% Computational Cost: add. (2465->180), mult. (5126->222), div. (0->0), fcn. (3021->6), ass. (0->74)
t334 = Ifges(4,1) + Ifges(5,1);
t322 = Ifges(4,4) - Ifges(5,5);
t332 = Ifges(5,4) + Ifges(4,5);
t333 = Ifges(4,2) + Ifges(5,3);
t331 = Ifges(4,6) - Ifges(5,6);
t330 = Ifges(4,3) + Ifges(5,2);
t300 = sin(qJ(3));
t303 = cos(qJ(2));
t316 = qJD(1) * t303;
t301 = sin(qJ(2));
t317 = qJD(1) * t301;
t326 = cos(qJ(3));
t278 = t300 * t317 - t326 * t316;
t279 = (t300 * t303 + t326 * t301) * qJD(1);
t298 = qJD(2) + qJD(3);
t329 = t333 * t278 - t322 * t279 - t331 * t298;
t328 = -t322 * t278 + t334 * t279 + t332 * t298;
t315 = qJD(1) * qJD(2);
t285 = qJDD(1) * t301 + t303 * t315;
t286 = qJDD(1) * t303 - t301 * t315;
t253 = qJD(3) * t279 + t285 * t300 - t326 * t286;
t254 = -t278 * qJD(3) + t326 * t285 + t300 * t286;
t289 = qJD(2) * pkin(2) - pkin(6) * t317;
t299 = t303 ^ 2;
t305 = qJD(1) ^ 2;
t302 = sin(qJ(1));
t304 = cos(qJ(1));
t312 = g(1) * t302 - t304 * g(2);
t308 = -qJDD(1) * pkin(1) - t312;
t255 = -pkin(2) * t286 + t289 * t317 + (-pkin(6) * t299 - pkin(5)) * t305 + t308;
t270 = -mrSges(4,2) * t298 - mrSges(4,3) * t278;
t271 = mrSges(4,1) * t298 - mrSges(4,3) * t279;
t272 = -mrSges(5,1) * t298 + mrSges(5,2) * t279;
t234 = -0.2e1 * qJD(4) * t279 + (t278 * t298 - t254) * qJ(4) + (t279 * t298 + t253) * pkin(3) + t255;
t273 = -mrSges(5,2) * t278 + mrSges(5,3) * t298;
t314 = m(5) * t234 + t253 * mrSges(5,1) + t278 * t273;
t327 = m(4) * t255 + mrSges(4,1) * t253 + (mrSges(4,2) - mrSges(5,3)) * t254 + t270 * t278 + (t271 - t272) * t279 + t314;
t325 = pkin(2) * t305;
t323 = -mrSges(4,3) - mrSges(5,2);
t309 = -g(1) * t304 - g(2) * t302;
t281 = -pkin(1) * t305 + qJDD(1) * pkin(5) + t309;
t321 = t281 * t301;
t246 = qJDD(2) * pkin(2) - pkin(6) * t285 - t321 + (pkin(6) * t315 + t301 * t325 - g(3)) * t303;
t269 = -g(3) * t301 + t303 * t281;
t247 = pkin(6) * t286 - qJD(2) * t289 - t299 * t325 + t269;
t241 = t300 * t246 + t326 * t247;
t297 = qJDD(2) + qJDD(3);
t263 = pkin(3) * t278 - qJ(4) * t279;
t296 = t298 ^ 2;
t237 = -pkin(3) * t296 + qJ(4) * t297 + 0.2e1 * qJD(4) * t298 - t263 * t278 + t241;
t313 = m(5) * t237 + t297 * mrSges(5,3) + t298 * t272;
t264 = mrSges(5,1) * t278 - mrSges(5,3) * t279;
t319 = -mrSges(4,1) * t278 - mrSges(4,2) * t279 - t264;
t226 = m(4) * t241 - mrSges(4,2) * t297 + t323 * t253 - t271 * t298 + t319 * t278 + t313;
t240 = t326 * t246 - t300 * t247;
t238 = -t297 * pkin(3) - t296 * qJ(4) + t279 * t263 + qJDD(4) - t240;
t310 = -m(5) * t238 + t297 * mrSges(5,1) + t298 * t273;
t228 = m(4) * t240 + mrSges(4,1) * t297 + t323 * t254 + t270 * t298 + t319 * t279 + t310;
t223 = t300 * t226 + t326 * t228;
t320 = t331 * t278 - t332 * t279 - t330 * t298;
t311 = t326 * t226 - t228 * t300;
t232 = mrSges(5,2) * t254 + t264 * t279 - t310;
t306 = mrSges(4,1) * t240 - mrSges(5,1) * t238 - mrSges(4,2) * t241 + mrSges(5,3) * t237 - pkin(3) * t232 + qJ(4) * t313 + t330 * t297 - t329 * t279 + (-qJ(4) * t264 + t328) * t278 + t332 * t254 + (-qJ(4) * mrSges(5,2) - t331) * t253;
t288 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t316;
t287 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t317;
t284 = (-t303 * mrSges(3,1) + t301 * mrSges(3,2)) * qJD(1);
t280 = -pkin(5) * t305 + t308;
t277 = Ifges(3,5) * qJD(2) + (t301 * Ifges(3,1) + t303 * Ifges(3,4)) * qJD(1);
t276 = Ifges(3,6) * qJD(2) + (t301 * Ifges(3,4) + t303 * Ifges(3,2)) * qJD(1);
t268 = -g(3) * t303 - t321;
t229 = -t254 * mrSges(5,3) - t272 * t279 + t314;
t222 = mrSges(4,2) * t255 + mrSges(5,2) * t238 - mrSges(4,3) * t240 - mrSges(5,3) * t234 - qJ(4) * t229 - t322 * t253 + t334 * t254 + t320 * t278 + t332 * t297 + t329 * t298;
t221 = -mrSges(4,1) * t255 - mrSges(5,1) * t234 + mrSges(5,2) * t237 + mrSges(4,3) * t241 - pkin(3) * t229 - t333 * t253 + t322 * t254 + t320 * t279 + t331 * t297 + t328 * t298;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t312 - mrSges(2,2) * t309 + t301 * (mrSges(3,2) * t280 - mrSges(3,3) * t268 + Ifges(3,1) * t285 + Ifges(3,4) * t286 + Ifges(3,5) * qJDD(2) - pkin(6) * t223 - qJD(2) * t276 - t300 * t221 + t326 * t222) + t303 * (-mrSges(3,1) * t280 + mrSges(3,3) * t269 + Ifges(3,4) * t285 + Ifges(3,2) * t286 + Ifges(3,6) * qJDD(2) - pkin(2) * t327 + pkin(6) * t311 + qJD(2) * t277 + t326 * t221 + t300 * t222) + pkin(1) * (-m(3) * t280 + mrSges(3,1) * t286 - mrSges(3,2) * t285 + (-t287 * t301 + t288 * t303) * qJD(1) - t327) + pkin(5) * (t303 * (m(3) * t269 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t286 - qJD(2) * t287 + t284 * t316 + t311) - t301 * (m(3) * t268 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t285 + qJD(2) * t288 - t284 * t317 + t223)); t306 + pkin(2) * t223 + Ifges(3,5) * t285 + Ifges(3,6) * t286 + Ifges(3,3) * qJDD(2) + (t301 * t276 - t303 * t277) * qJD(1) + mrSges(3,1) * t268 - mrSges(3,2) * t269; t306; t232;];
tauJ = t1;
