% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:43
% EndTime: 2019-12-05 16:15:44
% DurationCPUTime: 0.46s
% Computational Cost: add. (3469->124), mult. (4956->164), div. (0->0), fcn. (3254->10), ass. (0->67)
t294 = qJD(2) + qJD(3);
t290 = t294 ^ 2;
t298 = cos(pkin(9));
t324 = pkin(4) * t298;
t296 = sin(pkin(9));
t323 = mrSges(5,2) * t296;
t293 = t298 ^ 2;
t321 = t290 * t293;
t291 = qJDD(2) + qJDD(3);
t320 = t291 * t298;
t297 = sin(pkin(8));
t299 = cos(pkin(8));
t280 = t297 * g(1) - t299 * g(2);
t281 = -t299 * g(1) - t297 * g(2);
t302 = sin(qJ(2));
t305 = cos(qJ(2));
t313 = t305 * t280 - t302 * t281;
t269 = qJDD(2) * pkin(2) + t313;
t318 = t302 * t280 + t305 * t281;
t270 = -qJD(2) ^ 2 * pkin(2) + t318;
t301 = sin(qJ(3));
t304 = cos(qJ(3));
t257 = t301 * t269 + t304 * t270;
t254 = -t290 * pkin(3) + t291 * qJ(4) + t257;
t295 = -g(3) + qJDD(1);
t316 = qJD(4) * t294;
t317 = t298 * t295 - 0.2e1 * t296 * t316;
t247 = (-pkin(7) * t291 + t290 * t324 - t254) * t296 + t317;
t251 = t296 * t295 + (t254 + 0.2e1 * t316) * t298;
t248 = -pkin(4) * t321 + pkin(7) * t320 + t251;
t300 = sin(qJ(5));
t303 = cos(qJ(5));
t245 = t303 * t247 - t300 * t248;
t309 = -t296 * t300 + t298 * t303;
t273 = t309 * t294;
t310 = t296 * t303 + t298 * t300;
t274 = t310 * t294;
t263 = -t273 * mrSges(6,1) + t274 * mrSges(6,2);
t265 = t273 * qJD(5) + t310 * t291;
t271 = -qJD(5) * mrSges(6,2) + t273 * mrSges(6,3);
t243 = m(6) * t245 + qJDD(5) * mrSges(6,1) - t265 * mrSges(6,3) + qJD(5) * t271 - t274 * t263;
t246 = t300 * t247 + t303 * t248;
t264 = -t274 * qJD(5) + t309 * t291;
t272 = qJD(5) * mrSges(6,1) - t274 * mrSges(6,3);
t244 = m(6) * t246 - qJDD(5) * mrSges(6,2) + t264 * mrSges(6,3) - qJD(5) * t272 + t273 * t263;
t319 = t303 * t243 + t300 * t244;
t250 = -t296 * t254 + t317;
t311 = mrSges(5,3) * t291 + (-mrSges(5,1) * t298 + t323) * t290;
t234 = m(5) * t250 - t311 * t296 + t319;
t314 = -t300 * t243 + t303 * t244;
t235 = m(5) * t251 + t311 * t298 + t314;
t315 = -t234 * t296 + t298 * t235;
t256 = t304 * t269 - t301 * t270;
t312 = qJDD(4) - t256;
t292 = t296 ^ 2;
t249 = (-pkin(3) - t324) * t291 + (-qJ(4) + (-t292 - t293) * pkin(7)) * t290 + t312;
t258 = Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * qJD(5);
t260 = Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * qJD(5);
t236 = -mrSges(6,1) * t249 + mrSges(6,3) * t246 + Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * qJDD(5) + qJD(5) * t260 - t258 * t274;
t259 = Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * qJD(5);
t237 = mrSges(6,2) * t249 - mrSges(6,3) * t245 + Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * qJDD(5) - qJD(5) * t259 + t258 * t273;
t253 = -t291 * pkin(3) - t290 * qJ(4) + t312;
t307 = m(6) * t249 - t264 * mrSges(6,1) + t265 * mrSges(6,2) - t273 * t271 + t274 * t272;
t306 = -m(5) * t253 + mrSges(5,1) * t320 - t307 + (t290 * t292 + t321) * mrSges(5,3);
t239 = t291 * t323 - t306;
t308 = -mrSges(4,2) * t257 + t298 * (-mrSges(5,1) * t253 + mrSges(5,3) * t251 + t300 * t237 + t303 * t236 - pkin(4) * t307 + pkin(7) * t314 + (Ifges(5,4) * t296 + Ifges(5,2) * t298) * t291) + t296 * (mrSges(5,2) * t253 - mrSges(5,3) * t250 + t303 * t237 - t300 * t236 - pkin(7) * t319 + (Ifges(5,1) * t296 + Ifges(5,4) * t298) * t291) + qJ(4) * t315 - pkin(3) * t239 + mrSges(4,1) * t256 + Ifges(4,3) * t291;
t1 = [t298 * t234 + t296 * t235 + (m(2) + m(3) + m(4)) * t295; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t313 - mrSges(3,2) * t318 + pkin(2) * (t301 * (m(4) * t257 - mrSges(4,1) * t290 - mrSges(4,2) * t291 + t315) + t304 * (m(4) * t256 - t290 * mrSges(4,2) + (mrSges(4,1) - t323) * t291 + t306)) + t308; t308; t239; mrSges(6,1) * t245 - mrSges(6,2) * t246 + Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * qJDD(5) + t259 * t274 - t260 * t273;];
tauJ = t1;
