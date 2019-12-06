% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:26
% EndTime: 2019-12-05 15:40:27
% DurationCPUTime: 0.54s
% Computational Cost: add. (792->131), mult. (1374->160), div. (0->0), fcn. (617->6), ass. (0->57)
t279 = sin(qJ(4));
t281 = cos(qJ(4));
t304 = Ifges(5,4) - Ifges(6,5);
t312 = t281 * (Ifges(5,1) + Ifges(6,1)) - t279 * t304;
t311 = (Ifges(5,2) + Ifges(6,3)) * t279 - t281 * t304;
t303 = Ifges(6,4) + Ifges(5,5);
t302 = Ifges(5,6) - Ifges(6,6);
t308 = t281 * (t311 * qJD(2) - t302 * qJD(4));
t277 = sin(pkin(7));
t278 = cos(pkin(7));
t263 = -t278 * g(1) - t277 * g(2);
t274 = -g(3) + qJDD(1);
t280 = sin(qJ(2));
t282 = cos(qJ(2));
t240 = t282 * t263 + t280 * t274;
t307 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t240;
t306 = -pkin(2) - pkin(6);
t305 = -mrSges(5,3) - mrSges(6,2);
t295 = qJD(2) * qJD(4);
t259 = t281 * qJDD(2) - t279 * t295;
t301 = t259 * mrSges(6,3);
t262 = -t277 * g(1) + t278 * g(2);
t300 = t279 * t262;
t298 = t312 * qJD(2) + t303 * qJD(4);
t239 = -t280 * t263 + t282 * t274;
t284 = qJD(2) ^ 2;
t287 = -t284 * qJ(3) + qJDD(3) - t239;
t236 = t306 * qJDD(2) + t287;
t232 = t279 * t236 + t281 * t262;
t297 = qJD(2) * t279;
t296 = qJD(2) * t281;
t255 = (t279 * pkin(4) - t281 * qJ(5)) * qJD(2);
t283 = qJD(4) ^ 2;
t229 = -t283 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t255 * t297 + t232;
t266 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t296;
t293 = m(6) * t229 + qJDD(4) * mrSges(6,3) + qJD(4) * t266;
t256 = (t279 * mrSges(6,1) - t281 * mrSges(6,3)) * qJD(2);
t291 = qJD(2) * (-t256 - (t279 * mrSges(5,1) + t281 * mrSges(5,2)) * qJD(2));
t230 = -qJDD(4) * pkin(4) - t283 * qJ(5) + t300 + qJDD(5) + (qJD(2) * t255 - t236) * t281;
t267 = -mrSges(6,2) * t297 + qJD(4) * mrSges(6,3);
t290 = -m(6) * t230 + qJDD(4) * mrSges(6,1) + qJD(4) * t267;
t231 = t281 * t236 - t300;
t258 = t279 * qJDD(2) + t281 * t295;
t264 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t297;
t265 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t296;
t289 = t279 * (m(5) * t232 - qJDD(4) * mrSges(5,2) - qJD(4) * t265 + t305 * t258 + t279 * t291 + t293) + t281 * (m(5) * t231 + qJDD(4) * mrSges(5,1) + qJD(4) * t264 + t305 * t259 + t281 * t291 + t290);
t235 = t306 * t284 - t307;
t227 = t258 * pkin(4) - t259 * qJ(5) + (-0.2e1 * qJD(5) * t281 + (pkin(4) * t281 + qJ(5) * t279) * qJD(4)) * qJD(2) + t235;
t288 = m(6) * t227 + t258 * mrSges(6,1) - t266 * t296 + t267 * t297;
t238 = -qJDD(2) * pkin(2) + t287;
t286 = -m(4) * t238 + t284 * mrSges(4,3) - t289;
t237 = t284 * pkin(2) + t307;
t285 = -m(4) * t237 + m(5) * t235 + t284 * mrSges(4,2) + t259 * mrSges(5,2) + qJDD(2) * mrSges(4,3) + t264 * t297 + t265 * t296 + t288;
t225 = t259 * mrSges(6,2) + t256 * t296 - t290;
t224 = t288 - t301;
t221 = qJDD(2) * mrSges(4,2) - t286;
t1 = [m(2) * t274 + t280 * (m(3) * t240 - t284 * mrSges(3,1) + t258 * mrSges(5,1) - qJDD(2) * mrSges(3,2) + t285 - t301) + t282 * (m(3) * t239 - t284 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + t286); mrSges(3,1) * t239 - mrSges(3,2) * t240 + mrSges(4,2) * t238 - mrSges(4,3) * t237 + t281 * (mrSges(5,2) * t235 + mrSges(6,2) * t230 - mrSges(5,3) * t231 - mrSges(6,3) * t227 - qJ(5) * t224) - t279 * (-mrSges(5,1) * t235 - mrSges(6,1) * t227 + mrSges(6,2) * t229 + mrSges(5,3) * t232 - pkin(4) * t224) - pkin(6) * t289 - pkin(2) * t221 + qJ(3) * t285 + (-qJ(3) * mrSges(6,3) + t312) * t259 + (qJ(3) * mrSges(5,1) + t311) * t258 + (-t279 * t302 + t281 * t303) * qJDD(4) + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + (-t279 * t298 + t308) * qJD(4); t221; mrSges(5,1) * t231 - mrSges(5,2) * t232 - mrSges(6,1) * t230 + mrSges(6,3) * t229 - pkin(4) * t225 + qJ(5) * t293 + t303 * t259 + (-qJ(5) * mrSges(6,2) - t302) * t258 + (Ifges(5,3) + Ifges(6,2)) * qJDD(4) + (-t308 + (-qJ(5) * t256 + t298) * t279) * qJD(2); t225;];
tauJ = t1;
