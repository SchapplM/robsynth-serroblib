% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PPRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:32
% EndTime: 2019-12-31 16:18:33
% DurationCPUTime: 0.66s
% Computational Cost: add. (4662->131), mult. (7673->175), div. (0->0), fcn. (4854->8), ass. (0->62)
t312 = sin(pkin(6));
t314 = cos(pkin(6));
t306 = -t314 * g(1) - t312 * g(2);
t310 = -g(3) + qJDD(1);
t311 = sin(pkin(7));
t313 = cos(pkin(7));
t294 = -t311 * t306 + t313 * t310;
t295 = t313 * t306 + t311 * t310;
t316 = sin(qJ(3));
t318 = cos(qJ(3));
t291 = t316 * t294 + t318 * t295;
t319 = qJD(3) ^ 2;
t289 = -t319 * pkin(3) + qJDD(3) * pkin(5) + t291;
t305 = t312 * g(1) - t314 * g(2);
t304 = qJDD(2) - t305;
t315 = sin(qJ(4));
t317 = cos(qJ(4));
t286 = -t315 * t289 + t317 * t304;
t301 = (-mrSges(5,1) * t317 + mrSges(5,2) * t315) * qJD(3);
t327 = qJD(3) * qJD(4);
t302 = t315 * qJDD(3) + t317 * t327;
t328 = qJD(3) * t317;
t308 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t328;
t329 = qJD(3) * t315;
t283 = m(5) * t286 + qJDD(4) * mrSges(5,1) - t302 * mrSges(5,3) + qJD(4) * t308 - t301 * t329;
t287 = t317 * t289 + t315 * t304;
t303 = t317 * qJDD(3) - t315 * t327;
t307 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t329;
t284 = m(5) * t287 - qJDD(4) * mrSges(5,2) + t303 * mrSges(5,3) - qJD(4) * t307 + t301 * t328;
t274 = t317 * t283 + t315 * t284;
t273 = (m(3) + m(4)) * t304 + t274;
t297 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t315 + Ifges(5,2) * t317) * qJD(3);
t298 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t315 + Ifges(5,4) * t317) * qJD(3);
t332 = mrSges(5,1) * t286 - mrSges(5,2) * t287 + Ifges(5,5) * t302 + Ifges(5,6) * t303 + Ifges(5,3) * qJDD(4) + (t297 * t315 - t298 * t317) * qJD(3);
t275 = -t315 * t283 + t317 * t284;
t270 = m(4) * t291 - t319 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t275;
t290 = t318 * t294 - t316 * t295;
t288 = -qJDD(3) * pkin(3) - t319 * pkin(5) - t290;
t285 = -m(5) * t288 + t303 * mrSges(5,1) - t302 * mrSges(5,2) - t307 * t329 + t308 * t328;
t279 = m(4) * t290 + qJDD(3) * mrSges(4,1) - t319 * mrSges(4,2) + t285;
t267 = t316 * t270 + t318 * t279;
t265 = m(3) * t294 + t267;
t323 = t318 * t270 - t316 * t279;
t266 = m(3) * t295 + t323;
t324 = -t311 * t265 + t313 * t266;
t258 = m(2) * t306 + t324;
t272 = m(2) * t305 - t273;
t330 = t312 * t258 + t314 * t272;
t259 = t313 * t265 + t311 * t266;
t326 = m(2) * t310 + t259;
t325 = t314 * t258 - t272 * t312;
t296 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t315 + Ifges(5,6) * t317) * qJD(3);
t276 = -mrSges(5,1) * t288 + mrSges(5,3) * t287 + Ifges(5,4) * t302 + Ifges(5,2) * t303 + Ifges(5,6) * qJDD(4) + qJD(4) * t298 - t296 * t329;
t277 = mrSges(5,2) * t288 - mrSges(5,3) * t286 + Ifges(5,1) * t302 + Ifges(5,4) * t303 + Ifges(5,5) * qJDD(4) - qJD(4) * t297 + t296 * t328;
t320 = mrSges(4,1) * t290 - mrSges(4,2) * t291 + Ifges(4,3) * qJDD(3) + pkin(3) * t285 + pkin(5) * t275 + t317 * t276 + t315 * t277;
t261 = -mrSges(4,1) * t304 + mrSges(4,3) * t291 + t319 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t274 - t332;
t260 = mrSges(4,2) * t304 - mrSges(4,3) * t290 + Ifges(4,5) * qJDD(3) - t319 * Ifges(4,6) - pkin(5) * t274 - t315 * t276 + t317 * t277;
t255 = mrSges(3,2) * t304 - mrSges(3,3) * t294 - pkin(4) * t267 + t318 * t260 - t316 * t261;
t254 = -mrSges(3,1) * t304 + mrSges(3,3) * t295 + t316 * t260 + t318 * t261 - pkin(2) * (m(4) * t304 + t274) + pkin(4) * t323;
t253 = -mrSges(2,1) * t310 - mrSges(3,1) * t294 + mrSges(3,2) * t295 + mrSges(2,3) * t306 - pkin(1) * t259 - pkin(2) * t267 - t320;
t252 = mrSges(2,2) * t310 - mrSges(2,3) * t305 - qJ(2) * t259 - t254 * t311 + t255 * t313;
t1 = [-m(1) * g(1) + t325; -m(1) * g(2) + t330; -m(1) * g(3) + t326; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t330 + t314 * t252 - t312 * t253; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t325 + t312 * t252 + t314 * t253; -mrSges(1,1) * g(2) + mrSges(2,1) * t305 + mrSges(1,2) * g(1) - mrSges(2,2) * t306 - pkin(1) * t273 + qJ(2) * t324 + t313 * t254 + t311 * t255; t326; t273; t320; t332;];
tauJB = t1;
