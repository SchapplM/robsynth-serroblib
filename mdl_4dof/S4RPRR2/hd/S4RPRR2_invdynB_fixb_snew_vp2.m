% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:10
% EndTime: 2019-12-31 16:48:11
% DurationCPUTime: 0.71s
% Computational Cost: add. (8642->153), mult. (12133->196), div. (0->0), fcn. (6170->8), ass. (0->66)
t283 = qJD(1) + qJD(3);
t287 = sin(qJ(4));
t304 = t283 * t287;
t290 = cos(qJ(4));
t303 = t283 * t290;
t289 = sin(qJ(1));
t292 = cos(qJ(1));
t276 = t289 * g(1) - t292 * g(2);
t272 = qJDD(1) * pkin(1) + t276;
t277 = -t292 * g(1) - t289 * g(2);
t293 = qJD(1) ^ 2;
t273 = -t293 * pkin(1) + t277;
t285 = sin(pkin(7));
t286 = cos(pkin(7));
t259 = t286 * t272 - t285 * t273;
t257 = qJDD(1) * pkin(2) + t259;
t260 = t285 * t272 + t286 * t273;
t258 = -t293 * pkin(2) + t260;
t288 = sin(qJ(3));
t291 = cos(qJ(3));
t254 = t288 * t257 + t291 * t258;
t281 = t283 ^ 2;
t282 = qJDD(1) + qJDD(3);
t252 = -t281 * pkin(3) + t282 * pkin(6) + t254;
t284 = -g(3) + qJDD(2);
t249 = -t287 * t252 + t290 * t284;
t266 = (-mrSges(5,1) * t290 + mrSges(5,2) * t287) * t283;
t301 = qJD(4) * t283;
t267 = t287 * t282 + t290 * t301;
t275 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t303;
t247 = m(5) * t249 + qJDD(4) * mrSges(5,1) - t267 * mrSges(5,3) + qJD(4) * t275 - t266 * t304;
t250 = t290 * t252 + t287 * t284;
t268 = t290 * t282 - t287 * t301;
t274 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t304;
t248 = m(5) * t250 - qJDD(4) * mrSges(5,2) + t268 * mrSges(5,3) - qJD(4) * t274 + t266 * t303;
t296 = -t287 * t247 + t290 * t248;
t238 = m(4) * t254 - t281 * mrSges(4,1) - t282 * mrSges(4,2) + t296;
t253 = t291 * t257 - t288 * t258;
t251 = -t282 * pkin(3) - t281 * pkin(6) - t253;
t294 = -m(5) * t251 + t268 * mrSges(5,1) - t267 * mrSges(5,2) - t274 * t304 + t275 * t303;
t243 = m(4) * t253 + t282 * mrSges(4,1) - t281 * mrSges(4,2) + t294;
t235 = t288 * t238 + t291 * t243;
t232 = m(3) * t259 + qJDD(1) * mrSges(3,1) - t293 * mrSges(3,2) + t235;
t297 = t291 * t238 - t288 * t243;
t233 = m(3) * t260 - t293 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t297;
t227 = t286 * t232 + t285 * t233;
t225 = m(2) * t276 + qJDD(1) * mrSges(2,1) - t293 * mrSges(2,2) + t227;
t298 = -t285 * t232 + t286 * t233;
t226 = m(2) * t277 - t293 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t298;
t302 = t292 * t225 + t289 * t226;
t239 = t290 * t247 + t287 * t248;
t300 = m(4) * t284 + t239;
t299 = -t289 * t225 + t292 * t226;
t295 = m(3) * t284 + t300;
t263 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t287 + Ifges(5,4) * t290) * t283;
t262 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t287 + Ifges(5,2) * t290) * t283;
t261 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t287 + Ifges(5,6) * t290) * t283;
t241 = mrSges(5,2) * t251 - mrSges(5,3) * t249 + Ifges(5,1) * t267 + Ifges(5,4) * t268 + Ifges(5,5) * qJDD(4) - qJD(4) * t262 + t261 * t303;
t240 = -mrSges(5,1) * t251 + mrSges(5,3) * t250 + Ifges(5,4) * t267 + Ifges(5,2) * t268 + Ifges(5,6) * qJDD(4) + qJD(4) * t263 - t261 * t304;
t234 = -mrSges(4,1) * t284 - mrSges(5,1) * t249 + mrSges(5,2) * t250 + mrSges(4,3) * t254 + t281 * Ifges(4,5) - Ifges(5,5) * t267 + Ifges(4,6) * t282 - Ifges(5,6) * t268 - Ifges(5,3) * qJDD(4) - pkin(3) * t239 + (-t262 * t287 + t263 * t290) * t283;
t228 = mrSges(4,2) * t284 - mrSges(4,3) * t253 + Ifges(4,5) * t282 - t281 * Ifges(4,6) - pkin(6) * t239 - t287 * t240 + t290 * t241;
t221 = mrSges(3,2) * t284 - mrSges(3,3) * t259 + Ifges(3,5) * qJDD(1) - t293 * Ifges(3,6) - pkin(5) * t235 + t291 * t228 - t288 * t234;
t220 = -mrSges(3,1) * t284 + mrSges(3,3) * t260 + t293 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t300 + pkin(5) * t297 + t288 * t228 + t291 * t234;
t219 = -mrSges(2,2) * g(3) - mrSges(2,3) * t276 + Ifges(2,5) * qJDD(1) - t293 * Ifges(2,6) - qJ(2) * t227 - t285 * t220 + t286 * t221;
t218 = mrSges(2,1) * g(3) + mrSges(2,3) * t277 + t293 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t295 + qJ(2) * t298 + t286 * t220 + t285 * t221;
t1 = [-m(1) * g(1) + t299; -m(1) * g(2) + t302; (-m(1) - m(2)) * g(3) + t295; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t302 - t289 * t218 + t292 * t219; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t299 + t292 * t218 + t289 * t219; pkin(1) * t227 - mrSges(2,2) * t277 + mrSges(2,1) * t276 + pkin(2) * t235 + mrSges(3,1) * t259 - mrSges(3,2) * t260 + pkin(3) * t294 + pkin(6) * t296 + mrSges(4,1) * t253 - mrSges(4,2) * t254 + t287 * t241 + t290 * t240 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t282 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
