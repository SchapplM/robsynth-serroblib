% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR10_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:15
% EndTime: 2019-12-31 17:11:16
% DurationCPUTime: 0.93s
% Computational Cost: add. (1792->188), mult. (3665->228), div. (0->0), fcn. (1809->6), ass. (0->79)
t329 = Ifges(3,1) + Ifges(4,2);
t323 = Ifges(3,4) + Ifges(4,6);
t322 = Ifges(3,5) - Ifges(4,4);
t328 = Ifges(3,2) + Ifges(4,3);
t321 = Ifges(3,6) - Ifges(4,5);
t327 = Ifges(3,3) + Ifges(4,1);
t326 = -2 * qJD(3);
t300 = qJD(1) ^ 2;
t325 = pkin(5) * t300;
t324 = mrSges(3,1) - mrSges(4,2);
t295 = sin(qJ(1));
t298 = cos(qJ(1));
t309 = -g(1) * t298 - g(2) * t295;
t268 = -pkin(1) * t300 + qJDD(1) * pkin(5) + t309;
t294 = sin(qJ(2));
t297 = cos(qJ(2));
t255 = -t297 * g(3) - t294 * t268;
t320 = t327 * qJD(2) + (t294 * t322 + t297 * t321) * qJD(1);
t319 = t321 * qJD(2) + (t294 * t323 + t297 * t328) * qJD(1);
t318 = t322 * qJD(2) + (t294 * t329 + t297 * t323) * qJD(1);
t315 = qJD(1) * t297;
t282 = -mrSges(4,1) * t315 - qJD(2) * mrSges(4,3);
t317 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t315 - t282;
t316 = qJD(1) * t294;
t314 = qJD(1) * qJD(2);
t313 = t297 * t314;
t312 = t294 * t314;
t311 = g(1) * t295 - t298 * g(2);
t256 = -g(3) * t294 + t297 * t268;
t278 = qJDD(1) * t297 - t312;
t284 = pkin(3) * t316 - qJD(2) * pkin(6);
t292 = t297 ^ 2;
t277 = qJDD(1) * t294 + t313;
t307 = -qJDD(1) * pkin(1) - t311;
t301 = pkin(2) * t312 + t316 * t326 + (-t277 - t313) * qJ(3) + t307;
t235 = -t284 * t316 + (-pkin(3) * t292 - pkin(5)) * t300 + (-pkin(2) - pkin(6)) * t278 + t301;
t274 = (-t297 * pkin(2) - t294 * qJ(3)) * qJD(1);
t299 = qJD(2) ^ 2;
t242 = -qJDD(2) * pkin(2) - qJ(3) * t299 + t274 * t316 + qJDD(3) - t255;
t238 = (-t294 * t297 * t300 - qJDD(2)) * pkin(6) + (t277 - t313) * pkin(3) + t242;
t293 = sin(qJ(4));
t296 = cos(qJ(4));
t233 = -t235 * t293 + t238 * t296;
t272 = -qJD(2) * t293 - t296 * t315;
t251 = qJD(4) * t272 + qJDD(2) * t296 - t278 * t293;
t273 = qJD(2) * t296 - t293 * t315;
t252 = -mrSges(5,1) * t272 + mrSges(5,2) * t273;
t286 = qJD(4) + t316;
t253 = -mrSges(5,2) * t286 + mrSges(5,3) * t272;
t271 = qJDD(4) + t277;
t230 = m(5) * t233 + mrSges(5,1) * t271 - mrSges(5,3) * t251 - t252 * t273 + t253 * t286;
t234 = t235 * t296 + t238 * t293;
t250 = -qJD(4) * t273 - qJDD(2) * t293 - t278 * t296;
t254 = mrSges(5,1) * t286 - mrSges(5,3) * t273;
t231 = m(5) * t234 - mrSges(5,2) * t271 + mrSges(5,3) * t250 + t252 * t272 - t254 * t286;
t310 = -t230 * t293 + t296 * t231;
t226 = t296 * t230 + t293 * t231;
t239 = -pkin(2) * t278 + t301 - t325;
t308 = m(4) * t239 + t310;
t303 = -pkin(2) * t299 + qJDD(2) * qJ(3) + t274 * t315 + t256;
t237 = -pkin(6) * t292 * t300 + pkin(3) * t278 + ((2 * qJD(3)) + t284) * qJD(2) + t303;
t306 = -m(5) * t237 + t250 * mrSges(5,1) - t251 * mrSges(5,2) + t272 * t253 - t273 * t254;
t305 = m(4) * t242 + t277 * mrSges(4,1) + t226;
t244 = Ifges(5,4) * t273 + Ifges(5,2) * t272 + Ifges(5,6) * t286;
t245 = Ifges(5,1) * t273 + Ifges(5,4) * t272 + Ifges(5,5) * t286;
t304 = mrSges(5,1) * t233 - mrSges(5,2) * t234 + Ifges(5,5) * t251 + Ifges(5,6) * t250 + Ifges(5,3) * t271 + t273 * t244 - t272 * t245;
t241 = qJD(2) * t326 - t303;
t275 = (t297 * mrSges(4,2) - t294 * mrSges(4,3)) * qJD(1);
t283 = mrSges(4,1) * t316 + qJD(2) * mrSges(4,2);
t302 = -m(4) * t241 + qJDD(2) * mrSges(4,3) + qJD(2) * t283 + t275 * t315 - t306;
t280 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t316;
t276 = (-t297 * mrSges(3,1) + t294 * mrSges(3,2)) * qJD(1);
t267 = t307 - t325;
t243 = Ifges(5,5) * t273 + Ifges(5,6) * t272 + Ifges(5,3) * t286;
t228 = mrSges(5,2) * t237 - mrSges(5,3) * t233 + Ifges(5,1) * t251 + Ifges(5,4) * t250 + Ifges(5,5) * t271 + t243 * t272 - t244 * t286;
t227 = -mrSges(5,1) * t237 + mrSges(5,3) * t234 + Ifges(5,4) * t251 + Ifges(5,2) * t250 + Ifges(5,6) * t271 - t243 * t273 + t245 * t286;
t225 = qJDD(2) * mrSges(4,2) + qJD(2) * t282 + t275 * t316 + t305;
t224 = mrSges(4,2) * t278 - mrSges(4,3) * t277 + (t282 * t297 - t283 * t294) * qJD(1) + t308;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t311 - mrSges(2,2) * t309 + t294 * (mrSges(4,1) * t242 + mrSges(3,2) * t267 - mrSges(3,3) * t255 - mrSges(4,3) * t239 + pkin(3) * t226 - qJ(3) * t224 - t319 * qJD(2) + t322 * qJDD(2) + t329 * t277 + t323 * t278 + t320 * t315 + t304) + t297 * (-mrSges(3,1) * t267 - mrSges(4,1) * t241 + mrSges(4,2) * t239 + mrSges(3,3) * t256 - pkin(2) * t224 - pkin(3) * t306 - pkin(6) * t310 + t318 * qJD(2) + t321 * qJDD(2) - t296 * t227 - t293 * t228 + t323 * t277 + t328 * t278 - t320 * t316) + pkin(1) * (-m(3) * t267 + t324 * t278 + (-mrSges(3,2) + mrSges(4,3)) * t277 + (t317 * t297 + (-t280 + t283) * t294) * qJD(1) - t308) + pkin(5) * (t297 * (t276 * t315 + m(3) * t256 - qJDD(2) * mrSges(3,2) - qJD(2) * t280 + (mrSges(3,3) + mrSges(4,1)) * t278 + t302) + (-m(3) * t255 + t277 * mrSges(3,3) - t324 * qJDD(2) - t317 * qJD(2) + (t275 + t276) * t316 + t305) * t294); mrSges(3,1) * t255 - mrSges(3,2) * t256 + mrSges(4,2) * t242 - mrSges(4,3) * t241 + t296 * t228 - t293 * t227 - pkin(6) * t226 - pkin(2) * t225 + qJ(3) * t302 + (qJ(3) * mrSges(4,1) + t321) * t278 + t322 * t277 + t327 * qJDD(2) + (t319 * t294 - t318 * t297) * qJD(1); t225; t304;];
tauJ = t1;
