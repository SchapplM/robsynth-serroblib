% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:26
% EndTime: 2019-12-31 17:46:27
% DurationCPUTime: 0.51s
% Computational Cost: add. (2664->135), mult. (5057->175), div. (0->0), fcn. (2413->8), ass. (0->69)
t287 = cos(pkin(8));
t281 = t287 ^ 2;
t293 = qJD(1) ^ 2;
t315 = -pkin(1) - pkin(2);
t285 = sin(pkin(8));
t314 = t285 * mrSges(5,2);
t313 = t287 * mrSges(5,1);
t312 = t281 * t293;
t290 = sin(qJ(1));
t292 = cos(qJ(1));
t303 = -t292 * g(1) - t290 * g(2);
t297 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t303;
t261 = t293 * t315 + t297;
t306 = t290 * g(1) - t292 * g(2);
t296 = -t293 * qJ(2) + qJDD(2) - t306;
t264 = qJDD(1) * t315 + t296;
t286 = sin(pkin(7));
t288 = cos(pkin(7));
t249 = t288 * t261 + t286 * t264;
t247 = -pkin(3) * t293 - qJDD(1) * qJ(4) + t249;
t283 = g(3) + qJDD(3);
t308 = qJD(1) * qJD(4);
t310 = t287 * t283 + 0.2e1 * t285 * t308;
t240 = (pkin(4) * t287 * t293 + pkin(6) * qJDD(1) - t247) * t285 + t310;
t244 = t285 * t283 + (t247 - (2 * t308)) * t287;
t307 = t287 * qJDD(1);
t241 = -pkin(4) * t312 - pkin(6) * t307 + t244;
t289 = sin(qJ(5));
t291 = cos(qJ(5));
t238 = t240 * t291 - t241 * t289;
t300 = t285 * t289 - t287 * t291;
t267 = t300 * qJD(1);
t301 = -t285 * t291 - t287 * t289;
t268 = t301 * qJD(1);
t254 = -mrSges(6,1) * t267 + mrSges(6,2) * t268;
t257 = t267 * qJD(5) + qJDD(1) * t301;
t262 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t267;
t236 = m(6) * t238 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t257 + qJD(5) * t262 - t254 * t268;
t239 = t240 * t289 + t241 * t291;
t256 = -t268 * qJD(5) + qJDD(1) * t300;
t263 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t268;
t237 = m(6) * t239 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t256 - qJD(5) * t263 + t254 * t267;
t311 = t291 * t236 + t289 * t237;
t243 = -t247 * t285 + t310;
t299 = mrSges(5,3) * qJDD(1) + t293 * (t313 - t314);
t227 = m(5) * t243 + t285 * t299 + t311;
t304 = -t289 * t236 + t291 * t237;
t228 = m(5) * t244 - t287 * t299 + t304;
t305 = -t227 * t285 + t287 * t228;
t248 = -t286 * t261 + t288 * t264;
t225 = m(4) * t249 - mrSges(4,1) * t293 + qJDD(1) * mrSges(4,2) + t305;
t298 = qJDD(1) * pkin(3) + qJDD(4) - t248;
t246 = -t293 * qJ(4) + t298;
t280 = t285 ^ 2;
t242 = pkin(4) * t307 + (-qJ(4) + (-t280 - t281) * pkin(6)) * t293 + t298;
t295 = m(6) * t242 - t256 * mrSges(6,1) + t257 * mrSges(6,2) - t267 * t262 + t268 * t263;
t294 = -m(5) * t246 + qJDD(1) * t314 - t295 + (t280 * t293 + t312) * mrSges(5,3);
t231 = m(4) * t248 - t293 * mrSges(4,2) + (-mrSges(4,1) - t313) * qJDD(1) + t294;
t302 = t225 * t286 + t231 * t288;
t266 = -qJDD(1) * pkin(1) + t296;
t265 = -pkin(1) * t293 + t297;
t252 = Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * qJD(5);
t251 = Ifges(6,4) * t268 + Ifges(6,2) * t267 + Ifges(6,6) * qJD(5);
t250 = Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * qJD(5);
t232 = mrSges(5,1) * t307 - t294;
t230 = mrSges(6,2) * t242 - mrSges(6,3) * t238 + Ifges(6,1) * t257 + Ifges(6,4) * t256 + Ifges(6,5) * qJDD(5) - qJD(5) * t251 + t250 * t267;
t229 = -mrSges(6,1) * t242 + mrSges(6,3) * t239 + Ifges(6,4) * t257 + Ifges(6,2) * t256 + Ifges(6,6) * qJDD(5) + qJD(5) * t252 - t250 * t268;
t224 = m(3) * t266 - qJDD(1) * mrSges(3,1) - mrSges(3,3) * t293 + t302;
t1 = [-pkin(1) * t224 + qJ(2) * (m(3) * t265 - mrSges(3,1) * t293 + t225 * t288 - t231 * t286) + mrSges(2,1) * t306 - mrSges(2,2) * t303 - pkin(2) * t302 - mrSges(3,1) * t266 + mrSges(3,3) * t265 - t285 * (mrSges(5,2) * t246 - mrSges(5,3) * t243 - pkin(6) * t311 - t289 * t229 + t291 * t230) - t287 * (-mrSges(5,1) * t246 + mrSges(5,3) * t244 - pkin(4) * t295 + pkin(6) * t304 + t291 * t229 + t289 * t230) + pkin(3) * t232 - qJ(4) * t305 - mrSges(4,1) * t248 + mrSges(4,2) * t249 + (t281 * Ifges(5,2) + qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (Ifges(5,1) * t285 + 0.2e1 * Ifges(5,4) * t287) * t285) * qJDD(1); t224; m(4) * t283 + t227 * t287 + t228 * t285; t232; mrSges(6,1) * t238 - mrSges(6,2) * t239 + Ifges(6,5) * t257 + Ifges(6,6) * t256 + Ifges(6,3) * qJDD(5) + t251 * t268 - t252 * t267;];
tauJ = t1;
