% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR3
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:31
% EndTime: 2019-12-31 17:01:32
% DurationCPUTime: 0.72s
% Computational Cost: add. (9288->154), mult. (12133->196), div. (0->0), fcn. (6170->8), ass. (0->65)
t283 = qJD(1) + qJD(2);
t287 = sin(qJ(4));
t303 = t283 * t287;
t290 = cos(qJ(4));
t302 = t283 * t290;
t289 = sin(qJ(1));
t292 = cos(qJ(1));
t277 = t289 * g(1) - t292 * g(2);
t273 = qJDD(1) * pkin(1) + t277;
t278 = -t292 * g(1) - t289 * g(2);
t293 = qJD(1) ^ 2;
t274 = -t293 * pkin(1) + t278;
t288 = sin(qJ(2));
t291 = cos(qJ(2));
t260 = t291 * t273 - t288 * t274;
t282 = qJDD(1) + qJDD(2);
t258 = t282 * pkin(2) + t260;
t261 = t288 * t273 + t291 * t274;
t281 = t283 ^ 2;
t259 = -t281 * pkin(2) + t261;
t285 = sin(pkin(7));
t286 = cos(pkin(7));
t255 = t285 * t258 + t286 * t259;
t253 = -t281 * pkin(3) + t282 * pkin(6) + t255;
t284 = -g(3) + qJDD(3);
t250 = -t287 * t253 + t290 * t284;
t267 = (-mrSges(5,1) * t290 + mrSges(5,2) * t287) * t283;
t300 = qJD(4) * t283;
t268 = t287 * t282 + t290 * t300;
t276 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t302;
t248 = m(5) * t250 + qJDD(4) * mrSges(5,1) - t268 * mrSges(5,3) + qJD(4) * t276 - t267 * t303;
t251 = t290 * t253 + t287 * t284;
t269 = t290 * t282 - t287 * t300;
t275 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t303;
t249 = m(5) * t251 - qJDD(4) * mrSges(5,2) + t269 * mrSges(5,3) - qJD(4) * t275 + t267 * t302;
t295 = -t287 * t248 + t290 * t249;
t239 = m(4) * t255 - t281 * mrSges(4,1) - t282 * mrSges(4,2) + t295;
t254 = t286 * t258 - t285 * t259;
t252 = -t282 * pkin(3) - t281 * pkin(6) - t254;
t294 = -m(5) * t252 + t269 * mrSges(5,1) - t268 * mrSges(5,2) - t275 * t303 + t276 * t302;
t244 = m(4) * t254 + t282 * mrSges(4,1) - t281 * mrSges(4,2) + t294;
t236 = t285 * t239 + t286 * t244;
t233 = m(3) * t260 + t282 * mrSges(3,1) - t281 * mrSges(3,2) + t236;
t296 = t286 * t239 - t285 * t244;
t234 = m(3) * t261 - t281 * mrSges(3,1) - t282 * mrSges(3,2) + t296;
t228 = t291 * t233 + t288 * t234;
t226 = m(2) * t277 + qJDD(1) * mrSges(2,1) - t293 * mrSges(2,2) + t228;
t297 = -t288 * t233 + t291 * t234;
t227 = m(2) * t278 - t293 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t297;
t301 = t292 * t226 + t289 * t227;
t240 = t290 * t248 + t287 * t249;
t299 = m(4) * t284 + t240;
t298 = -t289 * t226 + t292 * t227;
t264 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t287 + Ifges(5,4) * t290) * t283;
t263 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t287 + Ifges(5,2) * t290) * t283;
t262 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t287 + Ifges(5,6) * t290) * t283;
t242 = mrSges(5,2) * t252 - mrSges(5,3) * t250 + Ifges(5,1) * t268 + Ifges(5,4) * t269 + Ifges(5,5) * qJDD(4) - qJD(4) * t263 + t262 * t302;
t241 = -mrSges(5,1) * t252 + mrSges(5,3) * t251 + Ifges(5,4) * t268 + Ifges(5,2) * t269 + Ifges(5,6) * qJDD(4) + qJD(4) * t264 - t262 * t303;
t235 = -mrSges(4,1) * t284 - mrSges(5,1) * t250 + mrSges(5,2) * t251 + mrSges(4,3) * t255 + t281 * Ifges(4,5) - Ifges(5,5) * t268 + Ifges(4,6) * t282 - Ifges(5,6) * t269 - Ifges(5,3) * qJDD(4) - pkin(3) * t240 + (-t263 * t287 + t264 * t290) * t283;
t229 = mrSges(4,2) * t284 - mrSges(4,3) * t254 + Ifges(4,5) * t282 - t281 * Ifges(4,6) - pkin(6) * t240 - t287 * t241 + t290 * t242;
t222 = -mrSges(3,2) * g(3) - mrSges(3,3) * t260 + Ifges(3,5) * t282 - t281 * Ifges(3,6) - qJ(3) * t236 + t286 * t229 - t285 * t235;
t221 = mrSges(3,1) * g(3) + mrSges(3,3) * t261 + t281 * Ifges(3,5) + Ifges(3,6) * t282 - pkin(2) * t299 + qJ(3) * t296 + t285 * t229 + t286 * t235;
t220 = -mrSges(2,2) * g(3) - mrSges(2,3) * t277 + Ifges(2,5) * qJDD(1) - t293 * Ifges(2,6) - pkin(5) * t228 - t288 * t221 + t291 * t222;
t219 = Ifges(2,6) * qJDD(1) + t293 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t278 + t288 * t222 + t291 * t221 - pkin(1) * (-m(3) * g(3) + t299) + pkin(5) * t297;
t1 = [-m(1) * g(1) + t298; -m(1) * g(2) + t301; (-m(1) - m(2) - m(3)) * g(3) + t299; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t301 - t289 * t219 + t292 * t220; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t298 + t292 * t219 + t289 * t220; pkin(1) * t228 - mrSges(2,2) * t278 + mrSges(2,1) * t277 + pkin(2) * t236 + mrSges(3,1) * t260 - mrSges(3,2) * t261 + t287 * t242 + t290 * t241 + pkin(3) * t294 + pkin(6) * t295 + mrSges(4,1) * t254 - mrSges(4,2) * t255 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (Ifges(4,3) + Ifges(3,3)) * t282;];
tauB = t1;
