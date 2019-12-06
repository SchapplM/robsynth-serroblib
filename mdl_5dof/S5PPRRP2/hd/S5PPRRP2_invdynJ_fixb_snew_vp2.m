% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:40
% EndTime: 2019-12-05 15:08:41
% DurationCPUTime: 0.51s
% Computational Cost: add. (926->124), mult. (1603->158), div. (0->0), fcn. (914->8), ass. (0->57)
t266 = sin(qJ(4));
t268 = cos(qJ(4));
t286 = Ifges(5,4) - Ifges(6,5);
t292 = t266 * (Ifges(5,1) + Ifges(6,1)) + t268 * t286;
t291 = -t266 * t286 - t268 * (Ifges(5,2) + Ifges(6,3));
t285 = Ifges(6,4) + Ifges(5,5);
t284 = Ifges(5,6) - Ifges(6,6);
t288 = t266 * (t291 * qJD(3) - t284 * qJD(4)) + t268 * (t292 * qJD(3) + t285 * qJD(4));
t287 = mrSges(5,3) + mrSges(6,2);
t277 = qJD(3) * qJD(4);
t249 = t268 * qJDD(3) - t266 * t277;
t283 = t249 * mrSges(6,1);
t263 = sin(pkin(7));
t265 = cos(pkin(7));
t252 = -t263 * g(1) + t265 * g(2) + qJDD(2);
t282 = t268 * t252;
t253 = -t265 * g(1) - t263 * g(2);
t261 = -g(3) + qJDD(1);
t262 = sin(pkin(8));
t264 = cos(pkin(8));
t230 = -t262 * t253 + t264 * t261;
t231 = t264 * t253 + t262 * t261;
t267 = sin(qJ(3));
t269 = cos(qJ(3));
t226 = t267 * t230 + t269 * t231;
t271 = qJD(3) ^ 2;
t224 = -t271 * pkin(3) + qJDD(3) * pkin(6) + t226;
t221 = t268 * t224 + t266 * t252;
t279 = qJD(3) * t266;
t278 = qJD(3) * t268;
t247 = (-t268 * mrSges(5,1) + t266 * mrSges(5,2)) * qJD(3);
t254 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t279;
t245 = (-t268 * pkin(4) - t266 * qJ(5)) * qJD(3);
t270 = qJD(4) ^ 2;
t218 = -t270 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t245 * t278 + t221;
t246 = (-t268 * mrSges(6,1) - t266 * mrSges(6,3)) * qJD(3);
t255 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t279;
t275 = m(6) * t218 + qJDD(4) * mrSges(6,3) + qJD(4) * t255 + t246 * t278;
t212 = m(5) * t221 - qJDD(4) * mrSges(5,2) - qJD(4) * t254 + t247 * t278 + t287 * t249 + t275;
t220 = -t266 * t224 + t282;
t248 = t266 * qJDD(3) + t268 * t277;
t256 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t278;
t219 = -qJDD(4) * pkin(4) - t270 * qJ(5) - t282 + qJDD(5) + (qJD(3) * t245 + t224) * t266;
t257 = mrSges(6,2) * t278 + qJD(4) * mrSges(6,3);
t274 = -m(6) * t219 + qJDD(4) * mrSges(6,1) + qJD(4) * t257;
t213 = m(5) * t220 + qJDD(4) * mrSges(5,1) + qJD(4) * t256 - t287 * t248 + (-t246 - t247) * t279 + t274;
t276 = t268 * t212 - t266 * t213;
t225 = t269 * t230 - t267 * t231;
t223 = -qJDD(3) * pkin(3) - t271 * pkin(6) - t225;
t216 = -t249 * pkin(4) - t248 * qJ(5) + (-0.2e1 * qJD(5) * t266 + (pkin(4) * t266 - qJ(5) * t268) * qJD(4)) * qJD(3) + t223;
t273 = m(6) * t216 - t248 * mrSges(6,3) - t255 * t279 - t257 * t278;
t272 = -m(5) * t223 + t249 * mrSges(5,1) - t254 * t279 + t256 * t278 - t273;
t215 = t248 * mrSges(6,2) + t246 * t279 - t274;
t214 = t273 - t283;
t210 = m(4) * t225 + qJDD(3) * mrSges(4,1) - t271 * mrSges(4,2) - t248 * mrSges(5,2) + t272 + t283;
t209 = m(4) * t226 - t271 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t276;
t1 = [m(2) * t261 + t262 * (m(3) * t231 + t269 * t209 - t267 * t210) + t264 * (m(3) * t230 + t267 * t209 + t269 * t210); t266 * t212 + t268 * t213 + (m(3) + m(4)) * t252; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t225 - mrSges(4,2) * t226 + t266 * (mrSges(5,2) * t223 + mrSges(6,2) * t219 - mrSges(5,3) * t220 - mrSges(6,3) * t216 - qJ(5) * t214) + t268 * (-mrSges(5,1) * t223 - mrSges(6,1) * t216 + mrSges(6,2) * t218 + mrSges(5,3) * t221 - pkin(4) * t214) + pkin(3) * t272 + pkin(6) * t276 + (pkin(3) * mrSges(6,1) - t291) * t249 + (-pkin(3) * mrSges(5,2) + t292) * t248 + (t266 * t285 + t268 * t284) * qJDD(4) + t288 * qJD(4); mrSges(5,1) * t220 - mrSges(5,2) * t221 - mrSges(6,1) * t219 + mrSges(6,3) * t218 - pkin(4) * t215 + qJ(5) * t275 + (qJ(5) * mrSges(6,2) + t284) * t249 + t285 * t248 + (Ifges(5,3) + Ifges(6,2)) * qJDD(4) - t288 * qJD(3); t215;];
tauJ = t1;
