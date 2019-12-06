% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:46
% EndTime: 2019-12-05 16:41:47
% DurationCPUTime: 0.48s
% Computational Cost: add. (1729->131), mult. (2293->164), div. (0->0), fcn. (1272->8), ass. (0->61)
t315 = Ifges(5,1) + Ifges(6,1);
t311 = Ifges(5,4) - Ifges(6,5);
t310 = Ifges(5,5) + Ifges(6,4);
t314 = Ifges(5,2) + Ifges(6,3);
t309 = Ifges(5,6) - Ifges(6,6);
t313 = Ifges(5,3) + Ifges(6,2);
t312 = mrSges(5,3) + mrSges(6,2);
t282 = qJD(2) + qJD(3);
t288 = sin(qJ(4));
t308 = t282 * t288;
t291 = cos(qJ(4));
t307 = t282 * t291;
t285 = -g(3) + qJDD(1);
t306 = t291 * t285;
t286 = sin(pkin(8));
t287 = cos(pkin(8));
t275 = t286 * g(1) - t287 * g(2);
t276 = -t287 * g(1) - t286 * g(2);
t290 = sin(qJ(2));
t293 = cos(qJ(2));
t299 = t293 * t275 - t290 * t276;
t247 = qJDD(2) * pkin(2) + t299;
t302 = t290 * t275 + t293 * t276;
t248 = -qJD(2) ^ 2 * pkin(2) + t302;
t289 = sin(qJ(3));
t292 = cos(qJ(3));
t243 = t289 * t247 + t292 * t248;
t280 = t282 ^ 2;
t281 = qJDD(2) + qJDD(3);
t240 = -t280 * pkin(3) + t281 * pkin(7) + t243;
t237 = t291 * t240 + t288 * t285;
t305 = (-t288 * t311 - t314 * t291) * t282 - t309 * qJD(4);
t304 = (t288 * t310 + t291 * t309) * t282 + t313 * qJD(4);
t303 = (t315 * t288 + t291 * t311) * t282 + t310 * qJD(4);
t301 = qJD(4) * t282;
t263 = (-mrSges(5,1) * t291 + mrSges(5,2) * t288) * t282;
t265 = t291 * t281 - t288 * t301;
t271 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t308;
t261 = (-pkin(4) * t291 - qJ(5) * t288) * t282;
t294 = qJD(4) ^ 2;
t234 = -t294 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t261 * t307 + t237;
t262 = (-mrSges(6,1) * t291 - mrSges(6,3) * t288) * t282;
t272 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t308;
t298 = m(6) * t234 + qJDD(4) * mrSges(6,3) + qJD(4) * t272 + t262 * t307;
t228 = m(5) * t237 - qJDD(4) * mrSges(5,2) - qJD(4) * t271 + t263 * t307 + t312 * t265 + t298;
t236 = -t288 * t240 + t306;
t264 = t288 * t281 + t291 * t301;
t273 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t307;
t235 = -qJDD(4) * pkin(4) - t294 * qJ(5) - t306 + qJDD(5) + (t261 * t282 + t240) * t288;
t274 = mrSges(6,2) * t307 + qJD(4) * mrSges(6,3);
t297 = -m(6) * t235 + qJDD(4) * mrSges(6,1) + qJD(4) * t274;
t229 = m(5) * t236 + qJDD(4) * mrSges(5,1) + qJD(4) * t273 + (-t262 - t263) * t308 - t312 * t264 + t297;
t300 = t291 * t228 - t288 * t229;
t242 = t292 * t247 - t289 * t248;
t239 = -t281 * pkin(3) - t280 * pkin(7) - t242;
t232 = -t265 * pkin(4) - t264 * qJ(5) + (-0.2e1 * qJD(5) * t288 + (pkin(4) * t288 - qJ(5) * t291) * qJD(4)) * t282 + t239;
t230 = m(6) * t232 - t265 * mrSges(6,1) - t264 * mrSges(6,3) - t272 * t308 - t274 * t307;
t295 = -m(5) * t239 + t265 * mrSges(5,1) - t264 * mrSges(5,2) - t271 * t308 + t273 * t307 - t230;
t296 = -mrSges(4,2) * t243 + t291 * (-mrSges(5,1) * t239 - mrSges(6,1) * t232 + mrSges(6,2) * t234 + mrSges(5,3) * t237 - pkin(4) * t230 + t303 * qJD(4) + t309 * qJDD(4) + t311 * t264 + t314 * t265 - t304 * t308) + t288 * (mrSges(5,2) * t239 + mrSges(6,2) * t235 - mrSges(5,3) * t236 - mrSges(6,3) * t232 - qJ(5) * t230 + t305 * qJD(4) + t310 * qJDD(4) + t315 * t264 + t311 * t265 + t304 * t307) + pkin(7) * t300 + pkin(3) * t295 + mrSges(4,1) * t242 + Ifges(4,3) * t281;
t231 = t264 * mrSges(6,2) + t262 * t308 - t297;
t1 = [t288 * t228 + t291 * t229 + (m(2) + m(3) + m(4)) * t285; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t299 - mrSges(3,2) * t302 + pkin(2) * (t289 * (m(4) * t243 - t280 * mrSges(4,1) - t281 * mrSges(4,2) + t300) + t292 * (m(4) * t242 + t281 * mrSges(4,1) - t280 * mrSges(4,2) + t295)) + t296; t296; mrSges(5,1) * t236 - mrSges(5,2) * t237 - mrSges(6,1) * t235 + mrSges(6,3) * t234 - pkin(4) * t231 + qJ(5) * t298 + (qJ(5) * mrSges(6,2) + t309) * t265 + t310 * t264 + t313 * qJDD(4) + (-t305 * t288 - t303 * t291) * t282; t231;];
tauJ = t1;
