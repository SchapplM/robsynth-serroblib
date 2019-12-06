% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR2
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:28
% EndTime: 2019-12-05 16:17:29
% DurationCPUTime: 0.55s
% Computational Cost: add. (2883->119), mult. (3981->163), div. (0->0), fcn. (2563->10), ass. (0->68)
t312 = 2 * qJD(4);
t282 = sin(pkin(8));
t284 = cos(pkin(8));
t267 = t282 * g(1) - t284 * g(2);
t269 = -t284 * g(1) - t282 * g(2);
t287 = sin(qJ(2));
t290 = cos(qJ(2));
t297 = t290 * t267 - t287 * t269;
t254 = qJDD(2) * pkin(2) + t297;
t303 = t287 * t267 + t290 * t269;
t255 = -qJD(2) ^ 2 * pkin(2) + t303;
t286 = sin(qJ(3));
t289 = cos(qJ(3));
t247 = t286 * t254 + t289 * t255;
t279 = qJD(2) + qJD(3);
t277 = t279 ^ 2;
t278 = qJDD(2) + qJDD(3);
t244 = -t277 * pkin(3) + t278 * qJ(4) + t247;
t311 = t279 * t312 + t244;
t281 = sin(pkin(9));
t310 = mrSges(5,2) * t281;
t308 = mrSges(5,3) * t278;
t307 = t279 * t281;
t283 = cos(pkin(9));
t306 = t283 * t278;
t305 = t283 * t279;
t280 = -g(3) + qJDD(1);
t304 = t283 * t280;
t301 = qJD(5) * t279;
t285 = sin(qJ(5));
t300 = t285 * t307;
t288 = cos(qJ(5));
t299 = t288 * t307;
t240 = t281 * t280 + t311 * t283;
t296 = -pkin(4) * t283 - pkin(7) * t281;
t263 = t296 * t279;
t238 = t263 * t305 + t240;
t246 = t289 * t254 - t286 * t255;
t294 = -t277 * qJ(4) + qJDD(4) - t246;
t241 = (-pkin(3) + t296) * t278 + t294;
t235 = -t285 * t238 + t288 * t241;
t270 = qJD(5) - t305;
t256 = -t270 * mrSges(6,2) - mrSges(6,3) * t300;
t258 = (t285 * mrSges(6,1) + t288 * mrSges(6,2)) * t307;
t260 = (t278 * t288 - t285 * t301) * t281;
t268 = qJDD(5) - t306;
t233 = m(6) * t235 + t268 * mrSges(6,1) - t260 * mrSges(6,3) + t270 * t256 - t258 * t299;
t236 = t288 * t238 + t285 * t241;
t257 = t270 * mrSges(6,1) - mrSges(6,3) * t299;
t259 = (-t278 * t285 - t288 * t301) * t281;
t234 = m(6) * t236 - t268 * mrSges(6,2) + t259 * mrSges(6,3) - t270 * t257 - t258 * t300;
t261 = (-mrSges(5,1) * t283 + t310) * t279;
t228 = m(5) * t240 - t285 * t233 + t288 * t234 + (t261 * t279 + t308) * t283;
t237 = -t304 + (t244 + (t312 + t263) * t279) * t281;
t239 = -t311 * t281 + t304;
t232 = m(5) * t239 - m(6) * t237 + t259 * mrSges(6,1) - t260 * mrSges(6,2) + (-t308 + (-t256 * t285 - t257 * t288 - t261) * t279) * t281;
t298 = t283 * t228 - t281 * t232;
t231 = t288 * t233 + t285 * t234;
t249 = Ifges(6,6) * t270 + (t288 * Ifges(6,4) - t285 * Ifges(6,2)) * t307;
t250 = Ifges(6,5) * t270 + (t288 * Ifges(6,1) - t285 * Ifges(6,4)) * t307;
t295 = t288 * t249 + t285 * t250;
t243 = -t278 * pkin(3) + t294;
t292 = -m(5) * t243 + mrSges(5,1) * t306 - t231 + (t281 ^ 2 + t283 ^ 2) * mrSges(5,3) * t277;
t230 = t278 * t310 - t292;
t262 = (Ifges(5,5) * t281 + Ifges(5,6) * t283) * t279;
t291 = mrSges(6,1) * t235 - mrSges(6,2) * t236 + Ifges(6,5) * t260 + Ifges(6,6) * t259 + Ifges(6,3) * t268;
t293 = -mrSges(4,2) * t247 + t281 * (t262 * t305 + mrSges(5,2) * t243 - mrSges(5,3) * t239 + t288 * (mrSges(6,2) * t237 - mrSges(6,3) * t235 + Ifges(6,1) * t260 + Ifges(6,4) * t259 + Ifges(6,5) * t268 - t270 * t249) - t285 * (-mrSges(6,1) * t237 + mrSges(6,3) * t236 + Ifges(6,4) * t260 + Ifges(6,2) * t259 + Ifges(6,6) * t268 + t270 * t250) - pkin(7) * t231 + (Ifges(5,1) * t281 + Ifges(5,4) * t283) * t278) + t283 * (Ifges(5,2) * t306 - mrSges(5,1) * t243 + mrSges(5,3) * t240 - pkin(4) * t231 + (Ifges(5,4) * t278 + (-t262 - t295) * t279) * t281 - t291) + qJ(4) * t298 - pkin(3) * t230 + mrSges(4,1) * t246 + Ifges(4,3) * t278;
t1 = [t281 * t228 + t283 * t232 + (m(2) + m(3) + m(4)) * t280; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t297 - mrSges(3,2) * t303 + pkin(2) * (t286 * (m(4) * t247 - t277 * mrSges(4,1) - t278 * mrSges(4,2) + t298) + t289 * (m(4) * t246 - t277 * mrSges(4,2) + (mrSges(4,1) - t310) * t278 + t292)) + t293; t293; t230; t295 * t307 + t291;];
tauJ = t1;
