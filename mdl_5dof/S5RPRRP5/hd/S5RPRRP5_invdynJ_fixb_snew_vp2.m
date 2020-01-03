% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:48
% EndTime: 2019-12-31 18:40:49
% DurationCPUTime: 0.54s
% Computational Cost: add. (2504->144), mult. (3386->180), div. (0->0), fcn. (1644->8), ass. (0->67)
t329 = Ifges(5,1) + Ifges(6,1);
t325 = Ifges(5,4) - Ifges(6,5);
t324 = Ifges(5,5) + Ifges(6,4);
t328 = Ifges(5,2) + Ifges(6,3);
t323 = Ifges(5,6) - Ifges(6,6);
t327 = Ifges(5,3) + Ifges(6,2);
t326 = mrSges(5,3) + mrSges(6,2);
t294 = qJD(1) + qJD(3);
t300 = sin(qJ(4));
t322 = t294 * t300;
t303 = cos(qJ(4));
t321 = t294 * t303;
t297 = -g(3) + qJDD(2);
t320 = t303 * t297;
t302 = sin(qJ(1));
t305 = cos(qJ(1));
t314 = t302 * g(1) - t305 * g(2);
t282 = qJDD(1) * pkin(1) + t314;
t307 = qJD(1) ^ 2;
t310 = -t305 * g(1) - t302 * g(2);
t283 = -t307 * pkin(1) + t310;
t298 = sin(pkin(8));
t299 = cos(pkin(8));
t258 = t299 * t282 - t298 * t283;
t256 = qJDD(1) * pkin(2) + t258;
t259 = t298 * t282 + t299 * t283;
t257 = -t307 * pkin(2) + t259;
t301 = sin(qJ(3));
t304 = cos(qJ(3));
t252 = t301 * t256 + t304 * t257;
t292 = t294 ^ 2;
t293 = qJDD(1) + qJDD(3);
t249 = -t292 * pkin(3) + t293 * pkin(7) + t252;
t246 = t303 * t249 + t300 * t297;
t274 = (-mrSges(5,1) * t303 + mrSges(5,2) * t300) * t294;
t315 = qJD(4) * t294;
t276 = t303 * t293 - t300 * t315;
t284 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t322;
t272 = (-pkin(4) * t303 - qJ(5) * t300) * t294;
t306 = qJD(4) ^ 2;
t243 = -t306 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t272 * t321 + t246;
t273 = (-mrSges(6,1) * t303 - mrSges(6,3) * t300) * t294;
t285 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t322;
t312 = m(6) * t243 + qJDD(4) * mrSges(6,3) + qJD(4) * t285 + t273 * t321;
t237 = m(5) * t246 - qJDD(4) * mrSges(5,2) - qJD(4) * t284 + t274 * t321 + t276 * t326 + t312;
t245 = -t300 * t249 + t320;
t275 = t300 * t293 + t303 * t315;
t286 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t321;
t244 = -qJDD(4) * pkin(4) - t306 * qJ(5) - t320 + qJDD(5) + (t272 * t294 + t249) * t300;
t287 = mrSges(6,2) * t321 + qJD(4) * mrSges(6,3);
t311 = -m(6) * t244 + qJDD(4) * mrSges(6,1) + qJD(4) * t287;
t238 = m(5) * t245 + qJDD(4) * mrSges(5,1) + qJD(4) * t286 + (-t273 - t274) * t322 - t326 * t275 + t311;
t313 = t303 * t237 - t300 * t238;
t231 = m(4) * t252 - t292 * mrSges(4,1) - t293 * mrSges(4,2) + t313;
t251 = t304 * t256 - t301 * t257;
t248 = -t293 * pkin(3) - t292 * pkin(7) - t251;
t241 = -t276 * pkin(4) - t275 * qJ(5) + (-0.2e1 * qJD(5) * t300 + (pkin(4) * t300 - qJ(5) * t303) * qJD(4)) * t294 + t248;
t239 = m(6) * t241 - t276 * mrSges(6,1) - t275 * mrSges(6,3) - t285 * t322 - t287 * t321;
t308 = -m(5) * t248 + t276 * mrSges(5,1) - t275 * mrSges(5,2) - t284 * t322 + t286 * t321 - t239;
t234 = m(4) * t251 + t293 * mrSges(4,1) - t292 * mrSges(4,2) + t308;
t319 = t301 * t231 + t304 * t234;
t318 = (-t300 * t325 - t328 * t303) * t294 - t323 * qJD(4);
t317 = (t300 * t324 + t303 * t323) * t294 + t327 * qJD(4);
t316 = (t329 * t300 + t303 * t325) * t294 + t324 * qJD(4);
t309 = -mrSges(4,2) * t252 + t303 * (-mrSges(5,1) * t248 - mrSges(6,1) * t241 + mrSges(6,2) * t243 + mrSges(5,3) * t246 - pkin(4) * t239 + t316 * qJD(4) + t323 * qJDD(4) + t325 * t275 + t328 * t276 - t317 * t322) + t300 * (mrSges(5,2) * t248 + mrSges(6,2) * t244 - mrSges(5,3) * t245 - mrSges(6,3) * t241 - qJ(5) * t239 + t318 * qJD(4) + t324 * qJDD(4) + t329 * t275 + t325 * t276 + t317 * t321) + pkin(7) * t313 + pkin(3) * t308 + mrSges(4,1) * t251 + Ifges(4,3) * t293;
t240 = t275 * mrSges(6,2) + t273 * t322 - t311;
t1 = [pkin(1) * (t298 * (m(3) * t259 - t307 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t304 * t231 - t301 * t234) + t299 * (m(3) * t258 + qJDD(1) * mrSges(3,1) - t307 * mrSges(3,2) + t319)) + mrSges(2,1) * t314 - mrSges(2,2) * t310 + pkin(2) * t319 + mrSges(3,1) * t258 - mrSges(3,2) * t259 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t309; t300 * t237 + t303 * t238 + (m(3) + m(4)) * t297; t309; mrSges(5,1) * t245 - mrSges(5,2) * t246 - mrSges(6,1) * t244 + mrSges(6,3) * t243 - pkin(4) * t240 + qJ(5) * t312 + (qJ(5) * mrSges(6,2) + t323) * t276 + t324 * t275 + t327 * qJDD(4) + (-t300 * t318 - t303 * t316) * t294; t240;];
tauJ = t1;
