% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:36
% EndTime: 2019-12-31 17:59:37
% DurationCPUTime: 0.53s
% Computational Cost: add. (2132->162), mult. (3830->203), div. (0->0), fcn. (1994->8), ass. (0->72)
t299 = sin(qJ(1));
t302 = cos(qJ(1));
t313 = t299 * g(1) - g(2) * t302;
t278 = qJDD(1) * pkin(1) + t313;
t304 = qJD(1) ^ 2;
t311 = -g(1) * t302 - g(2) * t299;
t280 = -pkin(1) * t304 + t311;
t295 = sin(pkin(8));
t296 = cos(pkin(8));
t262 = t295 * t278 + t296 * t280;
t323 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t262;
t322 = -pkin(2) - pkin(6);
t292 = -g(3) + qJDD(2);
t298 = sin(qJ(4));
t321 = t292 * t298;
t250 = t322 * t304 - t323;
t301 = cos(qJ(4));
t318 = qJD(1) * qJD(4);
t314 = t301 * t318;
t282 = -qJDD(1) * t298 - t314;
t315 = t298 * t318;
t283 = qJDD(1) * t301 - t315;
t243 = (-t283 + t315) * pkin(7) + (-t282 + t314) * pkin(4) + t250;
t261 = t278 * t296 - t295 * t280;
t309 = -qJ(3) * t304 + qJDD(3) - t261;
t251 = t322 * qJDD(1) + t309;
t247 = t298 * t251 + t301 * t292;
t281 = (pkin(4) * t298 - pkin(7) * t301) * qJD(1);
t303 = qJD(4) ^ 2;
t320 = qJD(1) * t298;
t245 = -pkin(4) * t303 + qJDD(4) * pkin(7) - t281 * t320 + t247;
t297 = sin(qJ(5));
t300 = cos(qJ(5));
t241 = t243 * t300 - t245 * t297;
t319 = qJD(1) * t301;
t276 = qJD(4) * t300 - t297 * t319;
t260 = qJD(5) * t276 + qJDD(4) * t297 + t283 * t300;
t277 = qJD(4) * t297 + t300 * t319;
t263 = -mrSges(6,1) * t276 + mrSges(6,2) * t277;
t286 = qJD(5) + t320;
t264 = -mrSges(6,2) * t286 + mrSges(6,3) * t276;
t275 = qJDD(5) - t282;
t239 = m(6) * t241 + mrSges(6,1) * t275 - mrSges(6,3) * t260 - t263 * t277 + t264 * t286;
t242 = t243 * t297 + t245 * t300;
t259 = -qJD(5) * t277 + qJDD(4) * t300 - t283 * t297;
t265 = mrSges(6,1) * t286 - mrSges(6,3) * t277;
t240 = m(6) * t242 - mrSges(6,2) * t275 + mrSges(6,3) * t259 + t263 * t276 - t265 * t286;
t232 = t300 * t239 + t297 * t240;
t312 = -t239 * t297 + t300 * t240;
t279 = (t298 * mrSges(5,1) + t301 * mrSges(5,2)) * qJD(1);
t285 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t319;
t231 = m(5) * t247 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t282 - qJD(4) * t285 - t279 * t320 + t312;
t246 = t251 * t301 - t321;
t284 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t320;
t244 = -qJDD(4) * pkin(4) - pkin(7) * t303 + t321 + (qJD(1) * t281 - t251) * t301;
t307 = -m(6) * t244 + t259 * mrSges(6,1) - mrSges(6,2) * t260 + t276 * t264 - t265 * t277;
t235 = m(5) * t246 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t283 + qJD(4) * t284 - t279 * t319 + t307;
t310 = t231 * t298 + t235 * t301;
t253 = -qJDD(1) * pkin(2) + t309;
t308 = m(4) * t253 - t304 * mrSges(4,3) + t310;
t252 = pkin(2) * t304 + t323;
t306 = -m(4) * t252 + m(5) * t250 - mrSges(5,1) * t282 + t304 * mrSges(4,2) + t283 * mrSges(5,2) + qJDD(1) * mrSges(4,3) + t284 * t320 + t285 * t319 + t232;
t255 = Ifges(6,4) * t277 + Ifges(6,2) * t276 + Ifges(6,6) * t286;
t256 = Ifges(6,1) * t277 + Ifges(6,4) * t276 + Ifges(6,5) * t286;
t305 = mrSges(6,1) * t241 - mrSges(6,2) * t242 + Ifges(6,5) * t260 + Ifges(6,6) * t259 + Ifges(6,3) * t275 + t277 * t255 - t276 * t256;
t271 = (Ifges(5,5) * qJD(4)) + (t301 * Ifges(5,1) - t298 * Ifges(5,4)) * qJD(1);
t270 = (Ifges(5,6) * qJD(4)) + (t301 * Ifges(5,4) - t298 * Ifges(5,2)) * qJD(1);
t254 = Ifges(6,5) * t277 + Ifges(6,6) * t276 + Ifges(6,3) * t286;
t234 = mrSges(6,2) * t244 - mrSges(6,3) * t241 + Ifges(6,1) * t260 + Ifges(6,4) * t259 + Ifges(6,5) * t275 + t254 * t276 - t255 * t286;
t233 = -mrSges(6,1) * t244 + mrSges(6,3) * t242 + Ifges(6,4) * t260 + Ifges(6,2) * t259 + Ifges(6,6) * t275 - t254 * t277 + t256 * t286;
t230 = qJDD(1) * mrSges(4,2) + t308;
t1 = [pkin(1) * (t295 * (m(3) * t262 - mrSges(3,1) * t304 + t306) + t296 * (m(3) * t261 - mrSges(3,2) * t304 - t308)) + mrSges(2,1) * t313 - mrSges(2,2) * t311 - pkin(2) * t230 + qJ(3) * t306 + t301 * (mrSges(5,2) * t250 - mrSges(5,3) * t246 + Ifges(5,1) * t283 + Ifges(5,4) * t282 + Ifges(5,5) * qJDD(4) - pkin(7) * t232 - qJD(4) * t270 - t233 * t297 + t300 * t234) - t298 * (-mrSges(5,1) * t250 + mrSges(5,3) * t247 + Ifges(5,4) * t283 + Ifges(5,2) * t282 + Ifges(5,6) * qJDD(4) - pkin(4) * t232 + qJD(4) * t271 - t305) - pkin(6) * t310 + mrSges(3,1) * t261 - mrSges(3,2) * t262 + mrSges(4,2) * t253 - mrSges(4,3) * t252 + (pkin(1) * (-t295 * mrSges(3,2) + t296 * (mrSges(3,1) - mrSges(4,2))) + Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1); t231 * t301 - t235 * t298 + (m(3) + m(4)) * t292; t230; Ifges(5,5) * t283 + Ifges(5,6) * t282 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t246 - mrSges(5,2) * t247 + t297 * t234 + t300 * t233 + pkin(4) * t307 + pkin(7) * t312 + (t270 * t301 + t271 * t298) * qJD(1); t305;];
tauJ = t1;
