% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:45
% EndTime: 2019-12-05 15:06:45
% DurationCPUTime: 0.52s
% Computational Cost: add. (934->130), mult. (1651->160), div. (0->0), fcn. (938->8), ass. (0->63)
t279 = cos(qJ(4));
t296 = Ifges(5,4) + Ifges(6,4);
t303 = t279 * t296;
t302 = Ifges(5,1) + Ifges(6,1);
t295 = Ifges(5,5) + Ifges(6,5);
t301 = Ifges(5,2) + Ifges(6,2);
t300 = Ifges(5,6) + Ifges(6,6);
t277 = sin(qJ(4));
t299 = -t295 * qJD(4) + (-t302 * t277 - t303) * qJD(3);
t281 = qJD(3) ^ 2;
t298 = pkin(4) * t281;
t297 = -mrSges(5,2) - mrSges(6,2);
t274 = sin(pkin(7));
t276 = cos(pkin(7));
t264 = -t276 * g(1) - t274 * g(2);
t272 = -g(3) + qJDD(1);
t273 = sin(pkin(8));
t275 = cos(pkin(8));
t242 = -t273 * t264 + t275 * t272;
t243 = t275 * t264 + t273 * t272;
t278 = sin(qJ(3));
t280 = cos(qJ(3));
t238 = t278 * t242 + t280 * t243;
t236 = -t281 * pkin(3) + qJDD(3) * pkin(6) + t238;
t263 = -t274 * g(1) + t276 * g(2) + qJDD(2);
t233 = t279 * t236 + t277 * t263;
t294 = t300 * qJD(4) + (t296 * t277 + t279 * t301) * qJD(3);
t291 = t277 * qJD(3);
t266 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t291;
t293 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t291 - t266;
t292 = qJD(3) * t279;
t290 = qJD(3) * qJD(4);
t289 = qJD(3) * qJD(5);
t252 = t279 * t263;
t286 = t279 * t290;
t260 = t277 * qJDD(3) + t286;
t229 = qJDD(4) * pkin(4) + t252 + (-t260 + t286) * qJ(5) + (t279 * t298 - t236 - 0.2e1 * t289) * t277;
t268 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t292;
t288 = m(6) * t229 + qJDD(4) * mrSges(6,1) + qJD(4) * t268;
t261 = t279 * qJDD(3) - t277 * t290;
t265 = qJD(4) * pkin(4) - qJ(5) * t291;
t271 = t279 ^ 2;
t230 = t261 * qJ(5) - qJD(4) * t265 - t271 * t298 + 0.2e1 * t279 * t289 + t233;
t258 = (-t279 * mrSges(6,1) + t277 * mrSges(6,2)) * qJD(3);
t287 = m(6) * t230 + t261 * mrSges(6,3) + t258 * t292;
t237 = t280 * t242 - t278 * t243;
t285 = qJD(3) * t293;
t283 = -qJDD(3) * pkin(3) - t237;
t231 = t265 * t291 - t261 * pkin(4) + qJDD(5) + (-qJ(5) * t271 - pkin(6)) * t281 + t283;
t284 = m(6) * t231 - t261 * mrSges(6,1) - t268 * t292;
t235 = -t281 * pkin(6) + t283;
t269 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t292;
t282 = -m(5) * t235 + t261 * mrSges(5,1) + t269 * t292 - t284;
t259 = (-t279 * mrSges(5,1) + t277 * mrSges(5,2)) * qJD(3);
t232 = -t277 * t236 + t252;
t226 = t260 * mrSges(6,2) + t266 * t291 + t284;
t225 = -t260 * mrSges(6,3) - t258 * t291 + t288;
t224 = m(5) * t233 + t261 * mrSges(5,3) + t293 * qJD(4) + t297 * qJDD(4) + t259 * t292 + t287;
t223 = m(5) * t232 + qJDD(4) * mrSges(5,1) + qJD(4) * t269 + (-mrSges(5,3) - mrSges(6,3)) * t260 + (-t258 - t259) * t291 + t288;
t222 = t279 * t224;
t221 = m(4) * t237 + qJDD(3) * mrSges(4,1) - t281 * mrSges(4,2) + t297 * t260 + t277 * t285 + t282;
t220 = m(4) * t238 - t281 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t277 * t223 + t222;
t1 = [m(2) * t272 + t273 * (m(3) * t243 + t280 * t220 - t278 * t221) + t275 * (m(3) * t242 + t278 * t220 + t280 * t221); t279 * t223 + t277 * t224 + (m(3) + m(4)) * t263; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t237 - mrSges(4,2) * t238 + t279 * (-mrSges(5,1) * t235 - mrSges(6,1) * t231 + mrSges(5,3) * t233 + mrSges(6,3) * t230 - pkin(4) * t226 + qJ(5) * t287 + t301 * t261 + (-qJ(5) * mrSges(6,2) + t300) * qJDD(4) + (-qJ(5) * t266 - t299) * qJD(4)) + pkin(3) * t282 + pkin(6) * t222 + (pkin(3) * t297 + t303) * t260 + (mrSges(5,2) * t235 + mrSges(6,2) * t231 - mrSges(5,3) * t232 - mrSges(6,3) * t229 + pkin(3) * t285 - pkin(6) * t223 - qJ(5) * t225 - t294 * qJD(4) + t295 * qJDD(4) + t302 * t260 + t296 * t261) * t277; mrSges(5,1) * t232 + mrSges(6,1) * t229 - mrSges(5,2) * t233 - mrSges(6,2) * t230 + pkin(4) * t225 + t300 * t261 + t295 * t260 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t294 * t277 + t279 * t299) * qJD(3); t226;];
tauJ = t1;
