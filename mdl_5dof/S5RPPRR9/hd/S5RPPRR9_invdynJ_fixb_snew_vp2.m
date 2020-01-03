% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:25
% EndTime: 2019-12-31 18:02:26
% DurationCPUTime: 0.52s
% Computational Cost: add. (3148->164), mult. (5508->207), div. (0->0), fcn. (2507->8), ass. (0->74)
t295 = -pkin(1) - pkin(2);
t267 = g(3) + qJDD(3);
t275 = cos(qJ(4));
t294 = t275 * t267;
t278 = qJD(1) ^ 2;
t273 = sin(qJ(1));
t276 = cos(qJ(1));
t285 = -t276 * g(1) - t273 * g(2);
t283 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t285;
t244 = t295 * t278 + t283;
t288 = t273 * g(1) - t276 * g(2);
t282 = -t278 * qJ(2) + qJDD(2) - t288;
t245 = t295 * qJDD(1) + t282;
t269 = sin(pkin(8));
t270 = cos(pkin(8));
t230 = t270 * t244 + t269 * t245;
t228 = -t278 * pkin(3) - qJDD(1) * pkin(6) + t230;
t272 = sin(qJ(4));
t225 = t275 * t228 + t272 * t267;
t293 = qJD(1) * t272;
t292 = t275 * qJD(1);
t291 = qJD(1) * qJD(4);
t290 = t272 * t291;
t289 = t275 * t291;
t256 = (t275 * mrSges(5,1) - t272 * mrSges(5,2)) * qJD(1);
t259 = -t275 * qJDD(1) + t290;
t260 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t293;
t229 = -t269 * t244 + t270 * t245;
t227 = qJDD(1) * pkin(3) - t278 * pkin(6) - t229;
t258 = -t272 * qJDD(1) - t289;
t221 = (-t258 + t289) * pkin(7) + (-t259 - t290) * pkin(4) + t227;
t257 = (t275 * pkin(4) + t272 * pkin(7)) * qJD(1);
t277 = qJD(4) ^ 2;
t223 = -t277 * pkin(4) + qJDD(4) * pkin(7) - t257 * t292 + t225;
t271 = sin(qJ(5));
t274 = cos(qJ(5));
t219 = t274 * t221 - t271 * t223;
t254 = t274 * qJD(4) + t271 * t293;
t237 = t254 * qJD(5) + t271 * qJDD(4) + t274 * t258;
t255 = t271 * qJD(4) - t274 * t293;
t238 = -t254 * mrSges(6,1) + t255 * mrSges(6,2);
t262 = qJD(5) + t292;
t242 = -t262 * mrSges(6,2) + t254 * mrSges(6,3);
t253 = qJDD(5) - t259;
t217 = m(6) * t219 + t253 * mrSges(6,1) - t237 * mrSges(6,3) - t255 * t238 + t262 * t242;
t220 = t271 * t221 + t274 * t223;
t236 = -t255 * qJD(5) + t274 * qJDD(4) - t271 * t258;
t243 = t262 * mrSges(6,1) - t255 * mrSges(6,3);
t218 = m(6) * t220 - t253 * mrSges(6,2) + t236 * mrSges(6,3) + t254 * t238 - t262 * t243;
t286 = -t271 * t217 + t274 * t218;
t211 = m(5) * t225 - qJDD(4) * mrSges(5,2) + t259 * mrSges(5,3) - qJD(4) * t260 - t256 * t292 + t286;
t224 = -t272 * t228 + t294;
t261 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t292;
t222 = -qJDD(4) * pkin(4) - t277 * pkin(7) - t294 + (-qJD(1) * t257 + t228) * t272;
t281 = -m(6) * t222 + t236 * mrSges(6,1) - t237 * mrSges(6,2) + t254 * t242 - t255 * t243;
t215 = m(5) * t224 + qJDD(4) * mrSges(5,1) - t258 * mrSges(5,3) + qJD(4) * t261 + t256 * t293 + t281;
t287 = t275 * t211 - t272 * t215;
t208 = m(4) * t230 - t278 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t287;
t212 = t274 * t217 + t271 * t218;
t280 = -m(5) * t227 + t259 * mrSges(5,1) - t258 * mrSges(5,2) + t260 * t293 - t261 * t292 - t212;
t209 = m(4) * t229 - qJDD(1) * mrSges(4,1) - t278 * mrSges(4,2) + t280;
t284 = t269 * t208 + t270 * t209;
t232 = Ifges(6,4) * t255 + Ifges(6,2) * t254 + Ifges(6,6) * t262;
t233 = Ifges(6,1) * t255 + Ifges(6,4) * t254 + Ifges(6,5) * t262;
t279 = mrSges(6,1) * t219 - mrSges(6,2) * t220 + Ifges(6,5) * t237 + Ifges(6,6) * t236 + Ifges(6,3) * t253 + t255 * t232 - t254 * t233;
t250 = (Ifges(5,5) * qJD(4)) + (-t272 * Ifges(5,1) - t275 * Ifges(5,4)) * qJD(1);
t249 = (Ifges(5,6) * qJD(4)) + (-t272 * Ifges(5,4) - t275 * Ifges(5,2)) * qJD(1);
t247 = -qJDD(1) * pkin(1) + t282;
t246 = -t278 * pkin(1) + t283;
t231 = Ifges(6,5) * t255 + Ifges(6,6) * t254 + Ifges(6,3) * t262;
t214 = mrSges(6,2) * t222 - mrSges(6,3) * t219 + Ifges(6,1) * t237 + Ifges(6,4) * t236 + Ifges(6,5) * t253 + t254 * t231 - t262 * t232;
t213 = -mrSges(6,1) * t222 + mrSges(6,3) * t220 + Ifges(6,4) * t237 + Ifges(6,2) * t236 + Ifges(6,6) * t253 - t255 * t231 + t262 * t233;
t207 = m(3) * t247 - qJDD(1) * mrSges(3,1) - t278 * mrSges(3,3) + t284;
t1 = [-pkin(1) * t207 + qJ(2) * (m(3) * t246 - t278 * mrSges(3,1) + t270 * t208 - t269 * t209) - mrSges(2,2) * t285 + mrSges(2,1) * t288 - pkin(2) * t284 - mrSges(3,1) * t247 + mrSges(3,3) * t246 - t272 * (mrSges(5,2) * t227 - mrSges(5,3) * t224 + Ifges(5,1) * t258 + Ifges(5,4) * t259 + Ifges(5,5) * qJDD(4) - pkin(7) * t212 - qJD(4) * t249 - t271 * t213 + t274 * t214) - t275 * (-mrSges(5,1) * t227 + mrSges(5,3) * t225 + Ifges(5,4) * t258 + Ifges(5,2) * t259 + Ifges(5,6) * qJDD(4) - pkin(4) * t212 + qJD(4) * t250 - t279) - pkin(3) * t280 - pkin(6) * t287 - mrSges(4,1) * t229 + mrSges(4,2) * t230 + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); t207; m(4) * t267 + t272 * t211 + t275 * t215; Ifges(5,5) * t258 + Ifges(5,6) * t259 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t224 - mrSges(5,2) * t225 + t271 * t214 + t274 * t213 + pkin(4) * t281 + pkin(7) * t286 + (-t272 * t249 + t275 * t250) * qJD(1); t279;];
tauJ = t1;
