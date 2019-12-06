% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:05
% EndTime: 2019-12-05 15:03:06
% DurationCPUTime: 0.71s
% Computational Cost: add. (6436->156), mult. (10014->194), div. (0->0), fcn. (5968->8), ass. (0->70)
t318 = -pkin(3) - pkin(6);
t317 = mrSges(4,1) - mrSges(5,2);
t316 = -Ifges(5,4) + Ifges(4,5);
t315 = Ifges(5,5) - Ifges(4,6);
t293 = sin(pkin(7));
t295 = cos(pkin(7));
t285 = -t295 * g(1) - t293 * g(2);
t289 = -g(3) + qJDD(1);
t292 = sin(pkin(8));
t294 = cos(pkin(8));
t271 = -t292 * t285 + t294 * t289;
t272 = t294 * t285 + t292 * t289;
t297 = sin(qJ(3));
t299 = cos(qJ(3));
t266 = t299 * t271 - t297 * t272;
t300 = qJD(3) ^ 2;
t304 = -t300 * qJ(4) + qJDD(4) - t266;
t263 = t318 * qJDD(3) + t304;
t284 = t293 * g(1) - t295 * g(2);
t283 = qJDD(2) - t284;
t296 = sin(qJ(5));
t298 = cos(qJ(5));
t259 = t298 * t263 - t296 * t283;
t280 = (mrSges(6,1) * t296 + mrSges(6,2) * t298) * qJD(3);
t310 = qJD(3) * qJD(5);
t282 = t298 * qJDD(3) - t296 * t310;
t312 = qJD(3) * t296;
t286 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t312;
t311 = qJD(3) * t298;
t257 = m(6) * t259 + qJDD(5) * mrSges(6,1) - t282 * mrSges(6,3) + qJD(5) * t286 - t280 * t311;
t260 = t296 * t263 + t298 * t283;
t281 = -t296 * qJDD(3) - t298 * t310;
t287 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t311;
t258 = m(6) * t260 - qJDD(5) * mrSges(6,2) + t281 * mrSges(6,3) - qJD(5) * t287 - t280 * t312;
t249 = t298 * t257 + t296 * t258;
t265 = -qJDD(3) * pkin(3) + t304;
t302 = -m(5) * t265 + t300 * mrSges(5,3) - t249;
t245 = m(4) * t266 - t300 * mrSges(4,2) + t317 * qJDD(3) + t302;
t267 = t297 * t271 + t299 * t272;
t303 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t267;
t264 = t300 * pkin(3) - t303;
t262 = t318 * t300 + t303;
t306 = -m(6) * t262 + t281 * mrSges(6,1) - t282 * mrSges(6,2) - t286 * t312 - t287 * t311;
t301 = -m(5) * t264 + t300 * mrSges(5,2) + qJDD(3) * mrSges(5,3) - t306;
t252 = m(4) * t267 - t300 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t301;
t243 = t299 * t245 + t297 * t252;
t241 = m(3) * t271 + t243;
t307 = -t297 * t245 + t299 * t252;
t242 = m(3) * t272 + t307;
t308 = -t292 * t241 + t294 * t242;
t234 = m(2) * t285 + t308;
t313 = -t296 * t257 + t298 * t258;
t248 = m(5) * t283 + t313;
t305 = (-m(3) - m(4)) * t283 - t248;
t247 = m(2) * t284 + t305;
t314 = t293 * t234 + t295 * t247;
t235 = t294 * t241 + t292 * t242;
t309 = t295 * t234 - t293 * t247;
t275 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t298 - Ifges(6,4) * t296) * qJD(3);
t274 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t298 - Ifges(6,2) * t296) * qJD(3);
t273 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t298 - Ifges(6,6) * t296) * qJD(3);
t254 = mrSges(6,2) * t262 - mrSges(6,3) * t259 + Ifges(6,1) * t282 + Ifges(6,4) * t281 + Ifges(6,5) * qJDD(5) - qJD(5) * t274 - t273 * t312;
t253 = -mrSges(6,1) * t262 + mrSges(6,3) * t260 + Ifges(6,4) * t282 + Ifges(6,2) * t281 + Ifges(6,6) * qJDD(5) + qJD(5) * t275 - t273 * t311;
t237 = mrSges(5,1) * t265 + mrSges(6,1) * t259 - mrSges(6,2) * t260 - mrSges(4,3) * t266 + Ifges(6,5) * t282 + Ifges(6,6) * t281 + Ifges(6,3) * qJDD(5) + pkin(4) * t249 - qJ(4) * t248 + t315 * t300 + (mrSges(4,2) - mrSges(5,3)) * t283 + t316 * qJDD(3) + (t298 * t274 + t296 * t275) * qJD(3);
t236 = -mrSges(5,1) * t264 + mrSges(4,3) * t267 - pkin(3) * t248 - pkin(4) * t306 - pkin(6) * t313 - t315 * qJDD(3) - t298 * t253 - t296 * t254 - t317 * t283 + t316 * t300;
t231 = mrSges(3,2) * t283 - mrSges(3,3) * t271 - pkin(5) * t243 - t297 * t236 + t299 * t237;
t230 = -mrSges(3,1) * t283 + mrSges(3,3) * t272 + t297 * t237 + t299 * t236 - pkin(2) * (m(4) * t283 + t248) + pkin(5) * t307;
t229 = -pkin(1) * t235 - pkin(2) * t243 - mrSges(3,1) * t271 + mrSges(3,2) * t272 - pkin(3) * t302 - qJ(4) * t301 - t298 * t254 + t296 * t253 + pkin(6) * t249 - mrSges(4,1) * t266 + mrSges(4,2) * t267 - mrSges(5,2) * t265 + mrSges(5,3) * t264 + mrSges(2,3) * t285 - mrSges(2,1) * t289 + (pkin(3) * mrSges(5,2) - Ifges(5,1) - Ifges(4,3)) * qJDD(3);
t228 = mrSges(2,2) * t289 - mrSges(2,3) * t284 - qJ(2) * t235 - t292 * t230 + t294 * t231;
t1 = [-m(1) * g(1) + t309; -m(1) * g(2) + t314; -m(1) * g(3) + m(2) * t289 + t235; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t314 + t295 * t228 - t293 * t229; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t309 + t293 * t228 + t295 * t229; -mrSges(1,1) * g(2) + mrSges(2,1) * t284 + mrSges(1,2) * g(1) - mrSges(2,2) * t285 + pkin(1) * t305 + qJ(2) * t308 + t294 * t230 + t292 * t231;];
tauB = t1;
