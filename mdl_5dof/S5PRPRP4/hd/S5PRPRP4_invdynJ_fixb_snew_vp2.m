% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:11
% EndTime: 2019-12-05 15:35:14
% DurationCPUTime: 0.63s
% Computational Cost: add. (1195->133), mult. (2042->167), div. (0->0), fcn. (1073->8), ass. (0->60)
t284 = sin(qJ(4));
t286 = cos(qJ(4));
t305 = Ifges(5,4) - Ifges(6,5);
t312 = t284 * (Ifges(5,1) + Ifges(6,1)) + t286 * t305;
t311 = -t284 * t305 - t286 * (Ifges(5,2) + Ifges(6,3));
t304 = Ifges(6,4) + Ifges(5,5);
t309 = Ifges(6,6) - Ifges(5,6);
t307 = t284 * (t311 * qJD(2) + t309 * qJD(4)) + t286 * (t312 * qJD(2) + t304 * qJD(4));
t306 = mrSges(5,3) + mrSges(6,2);
t295 = qJD(2) * qJD(4);
t265 = t286 * qJDD(2) - t284 * t295;
t302 = t265 * mrSges(6,1);
t281 = sin(pkin(7));
t283 = cos(pkin(7));
t268 = -t281 * g(1) + t283 * g(2) + qJDD(3);
t301 = t286 * t268;
t269 = -t283 * g(1) - t281 * g(2);
t279 = -g(3) + qJDD(1);
t285 = sin(qJ(2));
t287 = cos(qJ(2));
t245 = -t285 * t269 + t287 * t279;
t243 = qJDD(2) * pkin(2) + t245;
t246 = t287 * t269 + t285 * t279;
t289 = qJD(2) ^ 2;
t244 = -t289 * pkin(2) + t246;
t280 = sin(pkin(8));
t282 = cos(pkin(8));
t239 = t280 * t243 + t282 * t244;
t237 = -t289 * pkin(3) + qJDD(2) * pkin(6) + t239;
t234 = t286 * t237 + t284 * t268;
t263 = (-t286 * mrSges(5,1) + t284 * mrSges(5,2)) * qJD(2);
t297 = qJD(2) * t284;
t270 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t297;
t261 = (-t286 * pkin(4) - t284 * qJ(5)) * qJD(2);
t288 = qJD(4) ^ 2;
t296 = qJD(2) * t286;
t231 = -t288 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t261 * t296 + t234;
t262 = (-t286 * mrSges(6,1) - t284 * mrSges(6,3)) * qJD(2);
t271 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t297;
t293 = m(6) * t231 + qJDD(4) * mrSges(6,3) + qJD(4) * t271 + t262 * t296;
t225 = m(5) * t234 - qJDD(4) * mrSges(5,2) - qJD(4) * t270 + t263 * t296 + t306 * t265 + t293;
t233 = -t284 * t237 + t301;
t264 = t284 * qJDD(2) + t286 * t295;
t272 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t296;
t232 = -qJDD(4) * pkin(4) - t288 * qJ(5) - t301 + qJDD(5) + (qJD(2) * t261 + t237) * t284;
t273 = mrSges(6,2) * t296 + qJD(4) * mrSges(6,3);
t292 = -m(6) * t232 + qJDD(4) * mrSges(6,1) + qJD(4) * t273;
t226 = m(5) * t233 + qJDD(4) * mrSges(5,1) + qJD(4) * t272 - t306 * t264 + (-t262 - t263) * t297 + t292;
t294 = t286 * t225 - t284 * t226;
t221 = m(4) * t239 - t289 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t294;
t238 = t282 * t243 - t280 * t244;
t236 = -qJDD(2) * pkin(3) - t289 * pkin(6) - t238;
t229 = -t265 * pkin(4) - t264 * qJ(5) + (-0.2e1 * qJD(5) * t284 + (pkin(4) * t284 - qJ(5) * t286) * qJD(4)) * qJD(2) + t236;
t291 = m(6) * t229 - t264 * mrSges(6,3) - t271 * t297 - t273 * t296;
t290 = -m(5) * t236 + t265 * mrSges(5,1) - t270 * t297 + t272 * t296 - t291;
t223 = m(4) * t238 + qJDD(2) * mrSges(4,1) - t289 * mrSges(4,2) - t264 * mrSges(5,2) + t290 + t302;
t300 = t280 * t221 + t282 * t223;
t228 = t264 * mrSges(6,2) + t262 * t297 - t292;
t227 = t291 - t302;
t1 = [m(2) * t279 + t285 * (m(3) * t246 - t289 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t282 * t221 - t280 * t223) + t287 * (m(3) * t245 + qJDD(2) * mrSges(3,1) - t289 * mrSges(3,2) + t300); mrSges(3,1) * t245 - mrSges(3,2) * t246 + mrSges(4,1) * t238 - mrSges(4,2) * t239 + t284 * (mrSges(5,2) * t236 + mrSges(6,2) * t232 - mrSges(5,3) * t233 - mrSges(6,3) * t229 - qJ(5) * t227) + t286 * (-mrSges(5,1) * t236 - mrSges(6,1) * t229 + mrSges(6,2) * t231 + mrSges(5,3) * t234 - pkin(4) * t227) + pkin(3) * t290 + pkin(6) * t294 + pkin(2) * t300 + (pkin(3) * mrSges(6,1) - t311) * t265 + (-pkin(3) * mrSges(5,2) + t312) * t264 + (t284 * t304 - t286 * t309) * qJDD(4) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t307 * qJD(4); m(4) * t268 + t284 * t225 + t286 * t226; mrSges(5,1) * t233 - mrSges(5,2) * t234 - mrSges(6,1) * t232 + mrSges(6,3) * t231 - pkin(4) * t228 + qJ(5) * t293 + (qJ(5) * mrSges(6,2) - t309) * t265 + t304 * t264 + (Ifges(5,3) + Ifges(6,2)) * qJDD(4) - t307 * qJD(2); t228;];
tauJ = t1;
