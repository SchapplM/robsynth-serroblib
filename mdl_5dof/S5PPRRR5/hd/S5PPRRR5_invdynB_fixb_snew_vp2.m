% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:42
% EndTime: 2019-12-31 17:35:43
% DurationCPUTime: 0.74s
% Computational Cost: add. (8852->158), mult. (11295->197), div. (0->0), fcn. (6806->8), ass. (0->69)
t313 = -m(4) - m(5);
t312 = -mrSges(2,2) + mrSges(3,3);
t289 = qJD(3) + qJD(4);
t293 = sin(qJ(5));
t311 = t289 * t293;
t296 = cos(qJ(5));
t310 = t289 * t296;
t291 = sin(pkin(8));
t292 = cos(pkin(8));
t284 = g(1) * t291 - g(2) * t292;
t282 = qJDD(2) - t284;
t285 = -g(1) * t292 - g(2) * t291;
t295 = sin(qJ(3));
t298 = cos(qJ(3));
t267 = t298 * t282 - t285 * t295;
t265 = qJDD(3) * pkin(3) + t267;
t268 = t295 * t282 + t298 * t285;
t299 = qJD(3) ^ 2;
t266 = -pkin(3) * t299 + t268;
t294 = sin(qJ(4));
t297 = cos(qJ(4));
t262 = t294 * t265 + t297 * t266;
t287 = t289 ^ 2;
t288 = qJDD(3) + qJDD(4);
t260 = -pkin(4) * t287 + pkin(7) * t288 + t262;
t290 = g(3) - qJDD(1);
t257 = -t260 * t293 + t290 * t296;
t274 = (-mrSges(6,1) * t296 + mrSges(6,2) * t293) * t289;
t308 = qJD(5) * t289;
t275 = t288 * t293 + t296 * t308;
t281 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t310;
t255 = m(6) * t257 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t275 + qJD(5) * t281 - t274 * t311;
t258 = t260 * t296 + t290 * t293;
t276 = t288 * t296 - t293 * t308;
t280 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t311;
t256 = m(6) * t258 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t276 - qJD(5) * t280 + t274 * t310;
t303 = -t255 * t293 + t296 * t256;
t245 = m(5) * t262 - mrSges(5,1) * t287 - mrSges(5,2) * t288 + t303;
t261 = t265 * t297 - t266 * t294;
t259 = -pkin(4) * t288 - pkin(7) * t287 - t261;
t300 = -m(6) * t259 + t276 * mrSges(6,1) - mrSges(6,2) * t275 - t280 * t311 + t281 * t310;
t251 = m(5) * t261 + mrSges(5,1) * t288 - mrSges(5,2) * t287 + t300;
t242 = t294 * t245 + t297 * t251;
t240 = m(4) * t267 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t299 + t242;
t304 = t297 * t245 - t251 * t294;
t241 = m(4) * t268 - mrSges(4,1) * t299 - qJDD(3) * mrSges(4,2) + t304;
t236 = t240 * t298 + t241 * t295;
t301 = -m(3) * t282 - t236;
t234 = m(2) * t284 + t301;
t305 = -t240 * t295 + t298 * t241;
t302 = m(3) * t285 + t305;
t235 = m(2) * t285 + t302;
t309 = t292 * t234 + t291 * t235;
t247 = t296 * t255 + t293 * t256;
t307 = -m(3) * t290 - t247;
t306 = -t234 * t291 + t292 * t235;
t271 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t293 + Ifges(6,4) * t296) * t289;
t270 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t293 + Ifges(6,2) * t296) * t289;
t269 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t293 + Ifges(6,6) * t296) * t289;
t249 = mrSges(6,2) * t259 - mrSges(6,3) * t257 + Ifges(6,1) * t275 + Ifges(6,4) * t276 + Ifges(6,5) * qJDD(5) - qJD(5) * t270 + t269 * t310;
t248 = -mrSges(6,1) * t259 + mrSges(6,3) * t258 + Ifges(6,4) * t275 + Ifges(6,2) * t276 + Ifges(6,6) * qJDD(5) + qJD(5) * t271 - t269 * t311;
t246 = t290 * t313 + t307;
t239 = -mrSges(5,1) * t290 - mrSges(6,1) * t257 + mrSges(6,2) * t258 + mrSges(5,3) * t262 + t287 * Ifges(5,5) - Ifges(6,5) * t275 + Ifges(5,6) * t288 - Ifges(6,6) * t276 - Ifges(6,3) * qJDD(5) - pkin(4) * t247 + (-t270 * t293 + t271 * t296) * t289;
t237 = mrSges(5,2) * t290 - mrSges(5,3) * t261 + Ifges(5,5) * t288 - Ifges(5,6) * t287 - pkin(7) * t247 - t248 * t293 + t249 * t296;
t230 = mrSges(4,2) * t290 - mrSges(4,3) * t267 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t299 - pkin(6) * t242 + t237 * t297 - t239 * t294;
t229 = Ifges(4,6) * qJDD(3) + t299 * Ifges(4,5) - mrSges(4,1) * t290 + mrSges(4,3) * t268 + t294 * t237 + t297 * t239 - pkin(3) * (m(5) * t290 + t247) + pkin(6) * t304;
t228 = mrSges(3,2) * t282 - mrSges(2,3) * t284 - pkin(5) * t236 - qJ(2) * t246 - t295 * t229 + t298 * t230 + t290 * t312;
t227 = -t295 * t230 - t298 * t229 + pkin(2) * t247 - pkin(5) * t305 - pkin(1) * t246 + (mrSges(3,2) + mrSges(2,3)) * t285 + (-pkin(2) * t313 + mrSges(2,1) + mrSges(3,1)) * t290;
t1 = [-m(1) * g(1) + t306; -m(1) * g(2) + t309; -m(1) * g(3) + (-m(2) + t313) * t290 + t307; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t309 - t291 * t227 + t292 * t228; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t306 + t292 * t227 + t291 * t228; pkin(1) * t301 + qJ(2) * t302 + mrSges(2,1) * t284 - pkin(2) * t236 - mrSges(3,1) * t282 - pkin(3) * t242 - mrSges(4,1) * t267 + mrSges(4,2) * t268 - t293 * t249 - t296 * t248 - pkin(4) * t300 - pkin(7) * t303 - mrSges(5,1) * t261 + mrSges(5,2) * t262 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(5,3) * t288 - Ifges(4,3) * qJDD(3) + t312 * t285;];
tauB = t1;
