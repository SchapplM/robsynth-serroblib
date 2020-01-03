% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:32
% EndTime: 2019-12-31 19:49:33
% DurationCPUTime: 0.54s
% Computational Cost: add. (2945->143), mult. (3768->179), div. (0->0), fcn. (1833->8), ass. (0->66)
t330 = Ifges(5,1) + Ifges(6,1);
t326 = Ifges(5,4) - Ifges(6,5);
t325 = Ifges(5,5) + Ifges(6,4);
t329 = Ifges(5,2) + Ifges(6,3);
t324 = Ifges(5,6) - Ifges(6,6);
t328 = Ifges(5,3) + Ifges(6,2);
t327 = mrSges(5,3) + mrSges(6,2);
t296 = qJD(1) + qJD(2);
t302 = sin(qJ(4));
t323 = t296 * t302;
t305 = cos(qJ(4));
t322 = t296 * t305;
t299 = -g(3) + qJDD(3);
t321 = t305 * t299;
t304 = sin(qJ(1));
t307 = cos(qJ(1));
t315 = t304 * g(1) - t307 * g(2);
t283 = qJDD(1) * pkin(1) + t315;
t311 = -t307 * g(1) - t304 * g(2);
t284 = -qJD(1) ^ 2 * pkin(1) + t311;
t303 = sin(qJ(2));
t306 = cos(qJ(2));
t259 = t306 * t283 - t303 * t284;
t295 = qJDD(1) + qJDD(2);
t256 = t295 * pkin(2) + t259;
t260 = t303 * t283 + t306 * t284;
t294 = t296 ^ 2;
t257 = -t294 * pkin(2) + t260;
t300 = sin(pkin(8));
t301 = cos(pkin(8));
t252 = t300 * t256 + t301 * t257;
t249 = -t294 * pkin(3) + t295 * pkin(7) + t252;
t246 = t305 * t249 + t302 * t299;
t275 = (-mrSges(5,1) * t305 + mrSges(5,2) * t302) * t296;
t316 = qJD(4) * t296;
t277 = t305 * t295 - t302 * t316;
t285 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t323;
t273 = (-pkin(4) * t305 - qJ(5) * t302) * t296;
t308 = qJD(4) ^ 2;
t243 = -t308 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t273 * t322 + t246;
t274 = (-mrSges(6,1) * t305 - mrSges(6,3) * t302) * t296;
t286 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t323;
t313 = m(6) * t243 + qJDD(4) * mrSges(6,3) + qJD(4) * t286 + t274 * t322;
t237 = m(5) * t246 - qJDD(4) * mrSges(5,2) - qJD(4) * t285 + t275 * t322 + t327 * t277 + t313;
t245 = -t302 * t249 + t321;
t276 = t302 * t295 + t305 * t316;
t287 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t322;
t244 = -qJDD(4) * pkin(4) - t308 * qJ(5) - t321 + qJDD(5) + (t273 * t296 + t249) * t302;
t288 = mrSges(6,2) * t322 + qJD(4) * mrSges(6,3);
t312 = -m(6) * t244 + qJDD(4) * mrSges(6,1) + qJD(4) * t288;
t238 = m(5) * t245 + qJDD(4) * mrSges(5,1) + qJD(4) * t287 + (-t274 - t275) * t323 - t327 * t276 + t312;
t314 = t305 * t237 - t302 * t238;
t231 = m(4) * t252 - t294 * mrSges(4,1) - t295 * mrSges(4,2) + t314;
t251 = t301 * t256 - t300 * t257;
t248 = -t295 * pkin(3) - t294 * pkin(7) - t251;
t241 = -t277 * pkin(4) - t276 * qJ(5) + (-0.2e1 * qJD(5) * t302 + (pkin(4) * t302 - qJ(5) * t305) * qJD(4)) * t296 + t248;
t239 = m(6) * t241 - t277 * mrSges(6,1) - t276 * mrSges(6,3) - t286 * t323 - t288 * t322;
t309 = -m(5) * t248 + t277 * mrSges(5,1) - t276 * mrSges(5,2) - t285 * t323 + t287 * t322 - t239;
t234 = m(4) * t251 + t295 * mrSges(4,1) - t294 * mrSges(4,2) + t309;
t320 = t300 * t231 + t301 * t234;
t319 = (-t302 * t326 - t329 * t305) * t296 - t324 * qJD(4);
t318 = (t302 * t325 + t305 * t324) * t296 + t328 * qJD(4);
t317 = (t330 * t302 + t305 * t326) * t296 + t325 * qJD(4);
t310 = -mrSges(3,2) * t260 - mrSges(4,2) * t252 + pkin(2) * t320 + t305 * (-mrSges(5,1) * t248 - mrSges(6,1) * t241 + mrSges(6,2) * t243 + mrSges(5,3) * t246 - pkin(4) * t239 + t317 * qJD(4) + t324 * qJDD(4) + t326 * t276 + t329 * t277 - t318 * t323) + t302 * (mrSges(5,2) * t248 + mrSges(6,2) * t244 - mrSges(5,3) * t245 - mrSges(6,3) * t241 - qJ(5) * t239 + t319 * qJD(4) + t325 * qJDD(4) + t330 * t276 + t326 * t277 + t318 * t322) + pkin(7) * t314 + pkin(3) * t309 + mrSges(4,1) * t251 + mrSges(3,1) * t259 + (Ifges(4,3) + Ifges(3,3)) * t295;
t240 = t276 * mrSges(6,2) + t274 * t323 - t312;
t1 = [pkin(1) * (t303 * (m(3) * t260 - t294 * mrSges(3,1) - t295 * mrSges(3,2) + t301 * t231 - t300 * t234) + t306 * (m(3) * t259 + t295 * mrSges(3,1) - t294 * mrSges(3,2) + t320)) + mrSges(2,1) * t315 - mrSges(2,2) * t311 + Ifges(2,3) * qJDD(1) + t310; t310; m(4) * t299 + t302 * t237 + t305 * t238; mrSges(5,1) * t245 - mrSges(5,2) * t246 - mrSges(6,1) * t244 + mrSges(6,3) * t243 - pkin(4) * t240 + qJ(5) * t313 + (qJ(5) * mrSges(6,2) + t324) * t277 + t325 * t276 + t328 * qJDD(4) + (-t319 * t302 - t317 * t305) * t296; t240;];
tauJ = t1;
