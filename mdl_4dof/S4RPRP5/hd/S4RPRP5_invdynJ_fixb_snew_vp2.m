% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRP5
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:50
% EndTime: 2019-12-31 16:44:51
% DurationCPUTime: 0.55s
% Computational Cost: add. (1606->153), mult. (3820->190), div. (0->0), fcn. (2347->6), ass. (0->72)
t329 = Ifges(4,1) + Ifges(5,1);
t323 = Ifges(4,4) - Ifges(5,5);
t322 = Ifges(4,5) + Ifges(5,4);
t328 = -Ifges(4,2) - Ifges(5,3);
t321 = Ifges(4,6) - Ifges(5,6);
t327 = Ifges(4,3) + Ifges(5,2);
t296 = qJD(1) ^ 2;
t326 = cos(qJ(3));
t325 = pkin(2) * t296;
t324 = -mrSges(4,3) - mrSges(5,2);
t320 = pkin(5) * qJDD(1);
t293 = sin(qJ(1));
t294 = cos(qJ(1));
t302 = -t294 * g(1) - t293 * g(2);
t278 = -t296 * pkin(1) + qJDD(1) * qJ(2) + t302;
t290 = sin(pkin(6));
t291 = cos(pkin(6));
t310 = qJD(1) * qJD(2);
t305 = -t291 * g(3) - 0.2e1 * t290 * t310;
t250 = (t291 * t325 - t278 - t320) * t290 + t305;
t267 = -t290 * g(3) + (t278 + 0.2e1 * t310) * t291;
t287 = t291 ^ 2;
t251 = -t287 * t325 + t291 * t320 + t267;
t292 = sin(qJ(3));
t247 = t292 * t250 + t326 * t251;
t307 = t291 * t326;
t298 = t290 * t326 + t291 * t292;
t277 = t298 * qJD(1);
t311 = t277 * qJD(3);
t264 = t311 + (t290 * t292 - t307) * qJDD(1);
t271 = qJD(3) * mrSges(4,1) - t277 * mrSges(4,3);
t313 = qJD(1) * t290;
t276 = -qJD(1) * t307 + t292 * t313;
t259 = t276 * pkin(3) - t277 * qJ(4);
t295 = qJD(3) ^ 2;
t244 = -t295 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) - t276 * t259 + t247;
t272 = -qJD(3) * mrSges(5,1) + t277 * mrSges(5,2);
t308 = m(5) * t244 + qJDD(3) * mrSges(5,3) + qJD(3) * t272;
t260 = t276 * mrSges(5,1) - t277 * mrSges(5,3);
t315 = -t276 * mrSges(4,1) - t277 * mrSges(4,2) - t260;
t237 = m(4) * t247 - qJDD(3) * mrSges(4,2) - qJD(3) * t271 + t264 * t324 + t276 * t315 + t308;
t246 = t250 * t326 - t292 * t251;
t312 = t276 * qJD(3);
t265 = qJDD(1) * t298 - t312;
t270 = -qJD(3) * mrSges(4,2) - t276 * mrSges(4,3);
t245 = -qJDD(3) * pkin(3) - t295 * qJ(4) + t277 * t259 + qJDD(4) - t246;
t273 = -t276 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t303 = -m(5) * t245 + qJDD(3) * mrSges(5,1) + qJD(3) * t273;
t238 = m(4) * t246 + qJDD(3) * mrSges(4,1) + qJD(3) * t270 + t265 * t324 + t277 * t315 + t303;
t319 = t292 * t237 + t326 * t238;
t318 = -t327 * qJD(3) + t321 * t276 - t322 * t277;
t317 = t321 * qJD(3) + t328 * t276 + t323 * t277;
t316 = t322 * qJD(3) - t323 * t276 + t329 * t277;
t314 = -t290 ^ 2 - t287;
t306 = t293 * g(1) - t294 * g(2);
t301 = qJDD(2) - t306;
t263 = (-pkin(2) * t291 - pkin(1)) * qJDD(1) + (pkin(5) * t314 - qJ(2)) * t296 + t301;
t242 = -0.2e1 * qJD(4) * t277 + (-t265 + t312) * qJ(4) + (t264 + t311) * pkin(3) + t263;
t309 = m(5) * t242 + t264 * mrSges(5,1) + t276 * t273;
t304 = t326 * t237 - t292 * t238;
t300 = -t291 * mrSges(3,1) + t290 * mrSges(3,2);
t299 = mrSges(3,3) * qJDD(1) + t296 * t300;
t297 = m(4) * t263 + t264 * mrSges(4,1) + t276 * t270 + (t271 - t272) * t277 + (mrSges(4,2) - mrSges(5,3)) * t265 + t309;
t280 = (Ifges(3,5) * t290 + Ifges(3,6) * t291) * qJD(1);
t275 = -qJDD(1) * pkin(1) - t296 * qJ(2) + t301;
t266 = -t290 * t278 + t305;
t240 = t265 * mrSges(5,2) + t277 * t260 - t303;
t239 = -t265 * mrSges(5,3) - t277 * t272 + t309;
t233 = mrSges(3,3) * t296 * t314 + m(3) * t275 + qJDD(1) * t300 + t297;
t232 = mrSges(4,2) * t263 + mrSges(5,2) * t245 - mrSges(4,3) * t246 - mrSges(5,3) * t242 - qJ(4) * t239 - t317 * qJD(3) + t322 * qJDD(3) - t323 * t264 + t329 * t265 + t318 * t276;
t231 = -mrSges(4,1) * t263 - mrSges(5,1) * t242 + mrSges(5,2) * t244 + mrSges(4,3) * t247 - pkin(3) * t239 + t316 * qJD(3) + t321 * qJDD(3) + t328 * t264 + t323 * t265 + t318 * t277;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t306 - mrSges(2,2) * t302 + t290 * (t291 * qJD(1) * t280 + mrSges(3,2) * t275 - mrSges(3,3) * t266 + t326 * t232 - t292 * t231 - pkin(5) * t319 + (Ifges(3,1) * t290 + Ifges(3,4) * t291) * qJDD(1)) + t291 * (-t280 * t313 - mrSges(3,1) * t275 + mrSges(3,3) * t267 + t292 * t232 + t326 * t231 - pkin(2) * t297 + pkin(5) * t304 + (Ifges(3,4) * t290 + Ifges(3,2) * t291) * qJDD(1)) - pkin(1) * t233 + qJ(2) * ((m(3) * t267 + t291 * t299 + t304) * t291 + (-m(3) * t266 + t290 * t299 - t319) * t290); t233; mrSges(4,1) * t246 - mrSges(4,2) * t247 - mrSges(5,1) * t245 + mrSges(5,3) * t244 - pkin(3) * t240 + qJ(4) * t308 + t317 * t277 + (-qJ(4) * t260 + t316) * t276 + t322 * t265 + (-qJ(4) * mrSges(5,2) - t321) * t264 + t327 * qJDD(3); t240;];
tauJ = t1;
