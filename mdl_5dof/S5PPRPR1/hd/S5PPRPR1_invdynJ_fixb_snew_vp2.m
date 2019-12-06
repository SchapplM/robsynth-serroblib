% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:01
% EndTime: 2019-12-05 15:01:02
% DurationCPUTime: 0.40s
% Computational Cost: add. (1721->120), mult. (3332->163), div. (0->0), fcn. (2283->10), ass. (0->65)
t282 = qJD(3) ^ 2;
t275 = cos(pkin(9));
t270 = t275 ^ 2;
t298 = 0.2e1 * t275;
t297 = pkin(4) * t275;
t272 = sin(pkin(9));
t296 = t272 * mrSges(5,2);
t295 = t270 * t282;
t274 = sin(pkin(7));
t277 = cos(pkin(7));
t263 = -t277 * g(1) - t274 * g(2);
t271 = -g(3) + qJDD(1);
t273 = sin(pkin(8));
t276 = cos(pkin(8));
t254 = -t273 * t263 + t276 * t271;
t255 = t276 * t263 + t273 * t271;
t279 = sin(qJ(3));
t281 = cos(qJ(3));
t240 = t279 * t254 + t281 * t255;
t238 = -t282 * pkin(3) + qJDD(3) * qJ(4) + t240;
t262 = -t274 * g(1) + t277 * g(2) + qJDD(2);
t292 = qJD(3) * qJD(4);
t293 = t275 * t262 - 0.2e1 * t272 * t292;
t231 = (-pkin(6) * qJDD(3) + t282 * t297 - t238) * t272 + t293;
t234 = t275 * t238 + t272 * t262 + t292 * t298;
t291 = t275 * qJDD(3);
t232 = -pkin(4) * t295 + pkin(6) * t291 + t234;
t278 = sin(qJ(5));
t280 = cos(qJ(5));
t229 = t280 * t231 - t278 * t232;
t286 = -t272 * t278 + t275 * t280;
t256 = t286 * qJD(3);
t287 = t272 * t280 + t275 * t278;
t257 = t287 * qJD(3);
t245 = -t256 * mrSges(6,1) + t257 * mrSges(6,2);
t248 = t256 * qJD(5) + t287 * qJDD(3);
t252 = -qJD(5) * mrSges(6,2) + t256 * mrSges(6,3);
t226 = m(6) * t229 + qJDD(5) * mrSges(6,1) - t248 * mrSges(6,3) + qJD(5) * t252 - t257 * t245;
t230 = t278 * t231 + t280 * t232;
t247 = -t257 * qJD(5) + t286 * qJDD(3);
t253 = qJD(5) * mrSges(6,1) - t257 * mrSges(6,3);
t227 = m(6) * t230 - qJDD(5) * mrSges(6,2) + t247 * mrSges(6,3) - qJD(5) * t253 + t256 * t245;
t294 = t280 * t226 + t278 * t227;
t233 = -t272 * t238 + t293;
t285 = mrSges(5,3) * qJDD(3) + t282 * (-t275 * mrSges(5,1) + t296);
t218 = m(5) * t233 - t285 * t272 + t294;
t289 = -t278 * t226 + t280 * t227;
t219 = m(5) * t234 + t285 * t275 + t289;
t290 = -t272 * t218 + t275 * t219;
t239 = t281 * t254 - t279 * t255;
t288 = qJDD(4) - t239;
t269 = t272 ^ 2;
t235 = (-pkin(3) - t297) * qJDD(3) + (-qJ(4) + (-t269 - t270) * pkin(6)) * t282 + t288;
t284 = m(6) * t235 - t247 * mrSges(6,1) + t248 * mrSges(6,2) - t256 * t252 + t257 * t253;
t237 = -qJDD(3) * pkin(3) - t282 * qJ(4) + t288;
t283 = -m(5) * t237 + mrSges(5,1) * t291 - t284 + (t269 * t282 + t295) * mrSges(5,3);
t243 = Ifges(6,1) * t257 + Ifges(6,4) * t256 + Ifges(6,5) * qJD(5);
t242 = Ifges(6,4) * t257 + Ifges(6,2) * t256 + Ifges(6,6) * qJD(5);
t241 = Ifges(6,5) * t257 + Ifges(6,6) * t256 + Ifges(6,3) * qJD(5);
t228 = qJDD(3) * t296 - t283;
t222 = m(4) * t239 - t282 * mrSges(4,2) + (mrSges(4,1) - t296) * qJDD(3) + t283;
t221 = mrSges(6,2) * t235 - mrSges(6,3) * t229 + Ifges(6,1) * t248 + Ifges(6,4) * t247 + Ifges(6,5) * qJDD(5) - qJD(5) * t242 + t256 * t241;
t220 = -mrSges(6,1) * t235 + mrSges(6,3) * t230 + Ifges(6,4) * t248 + Ifges(6,2) * t247 + Ifges(6,6) * qJDD(5) + qJD(5) * t243 - t257 * t241;
t216 = m(4) * t240 - t282 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t290;
t1 = [m(2) * t271 + t273 * (m(3) * t255 + t281 * t216 - t279 * t222) + t276 * (m(3) * t254 + t279 * t216 + t281 * t222); t275 * t218 + t272 * t219 + (m(3) + m(4)) * t262; mrSges(4,1) * t239 - mrSges(4,2) * t240 + t272 * (mrSges(5,2) * t237 - mrSges(5,3) * t233 - pkin(6) * t294 - t278 * t220 + t280 * t221) + t275 * (-mrSges(5,1) * t237 + mrSges(5,3) * t234 - pkin(4) * t284 + pkin(6) * t289 + t280 * t220 + t278 * t221) - pkin(3) * t228 + qJ(4) * t290 + (Ifges(5,2) * t270 + Ifges(4,3) + (Ifges(5,1) * t272 + Ifges(5,4) * t298) * t272) * qJDD(3); t228; mrSges(6,1) * t229 - mrSges(6,2) * t230 + Ifges(6,5) * t248 + Ifges(6,6) * t247 + Ifges(6,3) * qJDD(5) + t257 * t242 - t256 * t243;];
tauJ = t1;
