% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:14
% EndTime: 2019-12-31 21:49:15
% DurationCPUTime: 0.57s
% Computational Cost: add. (4267->142), mult. (4419->178), div. (0->0), fcn. (2148->8), ass. (0->67)
t331 = Ifges(5,1) + Ifges(6,1);
t326 = Ifges(5,4) - Ifges(6,5);
t325 = Ifges(5,5) + Ifges(6,4);
t330 = Ifges(5,2) + Ifges(6,3);
t324 = Ifges(5,6) - Ifges(6,6);
t329 = Ifges(5,3) + Ifges(6,2);
t304 = cos(qJ(4));
t328 = t304 * g(3);
t327 = mrSges(5,3) + mrSges(6,2);
t297 = qJD(1) + qJD(2);
t293 = qJD(3) + t297;
t300 = sin(qJ(4));
t323 = t293 * t300;
t322 = t293 * t304;
t303 = sin(qJ(1));
t307 = cos(qJ(1));
t316 = t303 * g(1) - g(2) * t307;
t286 = qJDD(1) * pkin(1) + t316;
t312 = -g(1) * t307 - g(2) * t303;
t287 = -qJD(1) ^ 2 * pkin(1) + t312;
t302 = sin(qJ(2));
t306 = cos(qJ(2));
t258 = t306 * t286 - t287 * t302;
t296 = qJDD(1) + qJDD(2);
t255 = pkin(2) * t296 + t258;
t259 = t302 * t286 + t306 * t287;
t295 = t297 ^ 2;
t256 = -pkin(2) * t295 + t259;
t301 = sin(qJ(3));
t305 = cos(qJ(3));
t251 = t301 * t255 + t305 * t256;
t291 = t293 ^ 2;
t292 = qJDD(3) + t296;
t248 = -pkin(3) * t291 + pkin(8) * t292 + t251;
t244 = -t248 * t300 - t328;
t245 = -g(3) * t300 + t304 * t248;
t273 = (-mrSges(6,1) * t304 - mrSges(6,3) * t300) * t293;
t274 = (-mrSges(5,1) * t304 + mrSges(5,2) * t300) * t293;
t317 = qJD(4) * t293;
t275 = t292 * t300 + t304 * t317;
t276 = t292 * t304 - t300 * t317;
t282 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t323;
t284 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t322;
t272 = (-pkin(4) * t304 - qJ(5) * t300) * t293;
t308 = qJD(4) ^ 2;
t243 = -qJDD(4) * pkin(4) + t328 - t308 * qJ(5) + qJDD(5) + (t272 * t293 + t248) * t300;
t285 = mrSges(6,2) * t322 + qJD(4) * mrSges(6,3);
t313 = -m(6) * t243 + qJDD(4) * mrSges(6,1) + qJD(4) * t285;
t242 = -pkin(4) * t308 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t272 * t322 + t245;
t283 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t323;
t314 = m(6) * t242 + qJDD(4) * mrSges(6,3) + qJD(4) * t283 + t273 * t322;
t315 = -t300 * (m(5) * t244 + qJDD(4) * mrSges(5,1) + qJD(4) * t284 + (-t273 - t274) * t323 - t327 * t275 + t313) + t304 * (m(5) * t245 - qJDD(4) * mrSges(5,2) - qJD(4) * t282 + t274 * t322 + t276 * t327 + t314);
t231 = m(4) * t251 - mrSges(4,1) * t291 - mrSges(4,2) * t292 + t315;
t250 = t305 * t255 - t301 * t256;
t247 = -t292 * pkin(3) - t291 * pkin(8) - t250;
t240 = -t276 * pkin(4) - t275 * qJ(5) + (-0.2e1 * qJD(5) * t300 + (pkin(4) * t300 - qJ(5) * t304) * qJD(4)) * t293 + t247;
t238 = m(6) * t240 - mrSges(6,1) * t276 - t275 * mrSges(6,3) - t283 * t323 - t285 * t322;
t309 = -m(5) * t247 + t276 * mrSges(5,1) - mrSges(5,2) * t275 - t282 * t323 + t284 * t322 - t238;
t234 = m(4) * t250 + mrSges(4,1) * t292 - mrSges(4,2) * t291 + t309;
t321 = t301 * t231 + t305 * t234;
t320 = (-t300 * t326 - t304 * t330) * t293 - t324 * qJD(4);
t319 = (t300 * t325 + t304 * t324) * t293 + t329 * qJD(4);
t318 = (t300 * t331 + t304 * t326) * t293 + t325 * qJD(4);
t311 = -mrSges(4,2) * t251 + t304 * (-mrSges(5,1) * t247 - mrSges(6,1) * t240 + mrSges(6,2) * t242 + mrSges(5,3) * t245 - pkin(4) * t238 + t318 * qJD(4) + t324 * qJDD(4) + t326 * t275 + t276 * t330 - t319 * t323) + t300 * (mrSges(5,2) * t247 + mrSges(6,2) * t243 - mrSges(5,3) * t244 - mrSges(6,3) * t240 - qJ(5) * t238 + t320 * qJD(4) + t325 * qJDD(4) + t275 * t331 + t326 * t276 + t319 * t322) + pkin(8) * t315 + pkin(3) * t309 + mrSges(4,1) * t250 + Ifges(4,3) * t292;
t310 = mrSges(3,1) * t258 - mrSges(3,2) * t259 + Ifges(3,3) * t296 + pkin(2) * t321 + t311;
t239 = t275 * mrSges(6,2) + t273 * t323 - t313;
t1 = [pkin(1) * (t302 * (m(3) * t259 - mrSges(3,1) * t295 - mrSges(3,2) * t296 + t231 * t305 - t234 * t301) + t306 * (m(3) * t258 + mrSges(3,1) * t296 - mrSges(3,2) * t295 + t321)) + mrSges(2,1) * t316 - mrSges(2,2) * t312 + Ifges(2,3) * qJDD(1) + t310; t310; t311; mrSges(5,1) * t244 - mrSges(5,2) * t245 - mrSges(6,1) * t243 + mrSges(6,3) * t242 - pkin(4) * t239 + qJ(5) * t314 + (mrSges(6,2) * qJ(5) + t324) * t276 + t325 * t275 + t329 * qJDD(4) + (-t300 * t320 - t304 * t318) * t293; t239;];
tauJ = t1;
