% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPPR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:40
% EndTime: 2019-12-31 16:40:41
% DurationCPUTime: 0.52s
% Computational Cost: add. (937->135), mult. (2163->172), div. (0->0), fcn. (1208->6), ass. (0->64)
t283 = qJD(1) ^ 2;
t280 = sin(qJ(1));
t282 = cos(qJ(1));
t292 = -g(1) * t282 - g(2) * t280;
t311 = -pkin(1) * t283 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t292;
t277 = sin(pkin(6));
t278 = cos(pkin(6));
t253 = -t278 * g(3) - t277 * t311;
t310 = pkin(3) * t283;
t309 = Ifges(3,4) - Ifges(4,5);
t308 = pkin(5) * qJDD(1);
t307 = t277 * qJ(3);
t306 = t283 * qJ(2);
t290 = -pkin(2) * t278 - t307;
t265 = t290 * qJD(1);
t300 = qJD(1) * t277;
t242 = t265 * t300 + qJDD(3) - t253;
t239 = (-t278 * t310 - t308) * t277 + t242;
t254 = -t277 * g(3) + t278 * t311;
t299 = qJD(1) * t278;
t243 = t265 * t299 + t254;
t276 = t278 ^ 2;
t240 = -t276 * t310 - t278 * t308 + t243;
t279 = sin(qJ(4));
t281 = cos(qJ(4));
t237 = t239 * t281 - t240 * t279;
t288 = -t277 * t279 - t278 * t281;
t262 = t288 * qJD(1);
t289 = t277 * t281 - t278 * t279;
t263 = t289 * qJD(1);
t248 = -mrSges(5,1) * t262 + mrSges(5,2) * t263;
t252 = t262 * qJD(4) + qJDD(1) * t289;
t255 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t262;
t235 = m(5) * t237 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t252 + qJD(4) * t255 - t248 * t263;
t238 = t239 * t279 + t240 * t281;
t251 = -t263 * qJD(4) + qJDD(1) * t288;
t256 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t263;
t236 = m(5) * t238 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t251 - qJD(4) * t256 + t248 * t262;
t305 = t281 * t235 + t279 * t236;
t303 = ((Ifges(3,6) - Ifges(4,6)) * t278 + (Ifges(4,4) + Ifges(3,5)) * t277) * qJD(1);
t302 = t280 * g(1) - t282 * g(2);
t301 = -t277 ^ 2 - t276;
t296 = -qJDD(2) + t302;
t295 = t301 * mrSges(4,2);
t294 = -t279 * t235 + t281 * t236;
t293 = m(4) * t242 + t305;
t291 = -t278 * mrSges(4,1) - t277 * mrSges(4,3);
t266 = t291 * qJD(1);
t287 = qJDD(1) * mrSges(4,2) + qJD(1) * t266;
t286 = -0.2e1 * qJD(3) * t300 - t296;
t241 = (pkin(5) * t301 + qJ(2)) * t283 + (t307 + pkin(1) + (pkin(2) + pkin(3)) * t278) * qJDD(1) - t286;
t285 = -m(5) * t241 + t251 * mrSges(5,1) - t252 * mrSges(5,2) + t262 * t255 - t263 * t256;
t250 = -t306 + (-pkin(1) + t290) * qJDD(1) + t286;
t284 = m(4) * t250 + t285;
t267 = (-t278 * mrSges(3,1) + t277 * mrSges(3,2)) * qJD(1);
t261 = -qJDD(1) * pkin(1) - t296 - t306;
t246 = Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * qJD(4);
t245 = Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * qJD(4);
t244 = Ifges(5,5) * t263 + Ifges(5,6) * t262 + Ifges(5,3) * qJD(4);
t231 = qJDD(1) * t291 + t283 * t295 + t284;
t230 = m(3) * t261 + (mrSges(3,3) * t301 + t295) * t283 + ((-mrSges(3,1) - mrSges(4,1)) * t278 + (mrSges(3,2) - mrSges(4,3)) * t277) * qJDD(1) + t284;
t229 = mrSges(5,2) * t241 - mrSges(5,3) * t237 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * qJDD(4) - qJD(4) * t245 + t244 * t262;
t228 = -mrSges(5,1) * t241 + mrSges(5,3) * t238 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * qJDD(4) + qJD(4) * t246 - t244 * t263;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t302 - mrSges(2,2) * t292 + t277 * (mrSges(3,2) * t261 - mrSges(3,3) * t253 + mrSges(4,2) * t242 - mrSges(4,3) * t250 + t281 * t229 - t279 * t228 - pkin(5) * t305 - qJ(3) * t231 + t303 * t299 + (t278 * t309 + (Ifges(3,1) + Ifges(4,1)) * t277) * qJDD(1)) + t278 * (-mrSges(3,1) * t261 + mrSges(3,3) * t254 - mrSges(4,1) * t250 + mrSges(4,2) * t243 - t279 * t229 - t281 * t228 - pkin(3) * t285 - pkin(5) * t294 - pkin(2) * t231 - t303 * t300 + ((Ifges(3,2) + Ifges(4,3)) * t278 + t277 * t309) * qJDD(1)) - pkin(1) * t230 + qJ(2) * (t278 * (m(3) * t254 + m(4) * t243 + ((mrSges(4,2) + mrSges(3,3)) * qJDD(1) + (t266 + t267) * qJD(1)) * t278 + t294) + (-m(3) * t253 + (qJDD(1) * mrSges(3,3) + qJD(1) * t267 + t287) * t277 + t293) * t277); t230; t277 * t287 + t293; mrSges(5,1) * t237 - mrSges(5,2) * t238 + Ifges(5,5) * t252 + Ifges(5,6) * t251 + Ifges(5,3) * qJDD(4) + t245 * t263 - t246 * t262;];
tauJ = t1;
