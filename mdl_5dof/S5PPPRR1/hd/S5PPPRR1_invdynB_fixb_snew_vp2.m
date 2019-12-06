% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:59
% EndTime: 2019-12-05 14:58:01
% DurationCPUTime: 1.08s
% Computational Cost: add. (12028->151), mult. (18419->200), div. (0->0), fcn. (13408->10), ass. (0->71)
t276 = sin(pkin(7));
t279 = cos(pkin(7));
t270 = -t279 * g(1) - t276 * g(2);
t273 = -g(3) + qJDD(1);
t275 = sin(pkin(8));
t278 = cos(pkin(8));
t259 = t278 * t270 + t275 * t273;
t269 = t276 * g(1) - t279 * g(2);
t268 = qJDD(2) - t269;
t274 = sin(pkin(9));
t277 = cos(pkin(9));
t255 = -t274 * t259 + t277 * t268;
t256 = t277 * t259 + t274 * t268;
t281 = sin(qJ(4));
t283 = cos(qJ(4));
t252 = t281 * t255 + t283 * t256;
t284 = qJD(4) ^ 2;
t250 = -t284 * pkin(4) + qJDD(4) * pkin(6) + t252;
t258 = -t275 * t270 + t278 * t273;
t257 = qJDD(3) - t258;
t280 = sin(qJ(5));
t282 = cos(qJ(5));
t247 = -t280 * t250 + t282 * t257;
t265 = (-mrSges(6,1) * t282 + mrSges(6,2) * t280) * qJD(4);
t293 = qJD(4) * qJD(5);
t266 = t280 * qJDD(4) + t282 * t293;
t294 = qJD(4) * t282;
t272 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t294;
t295 = qJD(4) * t280;
t245 = m(6) * t247 + qJDD(5) * mrSges(6,1) - t266 * mrSges(6,3) + qJD(5) * t272 - t265 * t295;
t248 = t282 * t250 + t280 * t257;
t267 = t282 * qJDD(4) - t280 * t293;
t271 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t295;
t246 = m(6) * t248 - qJDD(5) * mrSges(6,2) + t267 * mrSges(6,3) - qJD(5) * t271 + t265 * t294;
t288 = -t280 * t245 + t282 * t246;
t234 = m(5) * t252 - t284 * mrSges(5,1) - qJDD(4) * mrSges(5,2) + t288;
t251 = t283 * t255 - t281 * t256;
t249 = -qJDD(4) * pkin(4) - t284 * pkin(6) - t251;
t285 = -m(6) * t249 + t267 * mrSges(6,1) - t266 * mrSges(6,2) - t271 * t295 + t272 * t294;
t241 = m(5) * t251 + qJDD(4) * mrSges(5,1) - t284 * mrSges(5,2) + t285;
t231 = t281 * t234 + t283 * t241;
t229 = m(4) * t255 + t231;
t289 = t283 * t234 - t281 * t241;
t230 = m(4) * t256 + t289;
t290 = -t274 * t229 + t277 * t230;
t224 = m(3) * t259 + t290;
t237 = t282 * t245 + t280 * t246;
t287 = (-m(4) - m(5)) * t257 - t237;
t236 = m(3) * t258 + t287;
t291 = t278 * t224 - t275 * t236;
t218 = m(2) * t270 + t291;
t225 = t277 * t229 + t274 * t230;
t286 = -m(3) * t268 - t225;
t223 = m(2) * t269 + t286;
t296 = t276 * t218 + t279 * t223;
t219 = t275 * t224 + t278 * t236;
t292 = t279 * t218 - t276 * t223;
t262 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t280 + Ifges(6,4) * t282) * qJD(4);
t261 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t280 + Ifges(6,2) * t282) * qJD(4);
t260 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t280 + Ifges(6,6) * t282) * qJD(4);
t239 = mrSges(6,2) * t249 - mrSges(6,3) * t247 + Ifges(6,1) * t266 + Ifges(6,4) * t267 + Ifges(6,5) * qJDD(5) - qJD(5) * t261 + t260 * t294;
t238 = -mrSges(6,1) * t249 + mrSges(6,3) * t248 + Ifges(6,4) * t266 + Ifges(6,2) * t267 + Ifges(6,6) * qJDD(5) + qJD(5) * t262 - t260 * t295;
t227 = -mrSges(5,1) * t257 - mrSges(6,1) * t247 + mrSges(6,2) * t248 + mrSges(5,3) * t252 + t284 * Ifges(5,5) - Ifges(6,5) * t266 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t267 - Ifges(6,3) * qJDD(5) - pkin(4) * t237 + (-t261 * t280 + t262 * t282) * qJD(4);
t226 = mrSges(5,2) * t257 - mrSges(5,3) * t251 + Ifges(5,5) * qJDD(4) - t284 * Ifges(5,6) - pkin(6) * t237 - t280 * t238 + t282 * t239;
t215 = mrSges(4,2) * t257 - mrSges(4,3) * t255 - pkin(5) * t231 + t283 * t226 - t281 * t227;
t214 = -mrSges(4,1) * t257 + mrSges(4,3) * t256 + t281 * t226 + t283 * t227 - pkin(3) * (m(5) * t257 + t237) + pkin(5) * t289;
t213 = -mrSges(3,1) * t268 - mrSges(4,1) * t255 - mrSges(5,1) * t251 + mrSges(4,2) * t256 + mrSges(5,2) * t252 + mrSges(3,3) * t259 - Ifges(5,3) * qJDD(4) - pkin(2) * t225 - pkin(3) * t231 - pkin(4) * t285 - pkin(6) * t288 - t282 * t238 - t280 * t239;
t212 = mrSges(3,2) * t268 - mrSges(3,3) * t258 - qJ(3) * t225 - t274 * t214 + t277 * t215;
t211 = -mrSges(2,1) * t273 - mrSges(3,1) * t258 + mrSges(3,2) * t259 + mrSges(2,3) * t270 - pkin(1) * t219 - pkin(2) * t287 - qJ(3) * t290 - t277 * t214 - t274 * t215;
t210 = mrSges(2,2) * t273 - mrSges(2,3) * t269 - qJ(2) * t219 + t278 * t212 - t275 * t213;
t1 = [-m(1) * g(1) + t292; -m(1) * g(2) + t296; -m(1) * g(3) + m(2) * t273 + t219; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t296 + t279 * t210 - t276 * t211; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t292 + t276 * t210 + t279 * t211; -mrSges(1,1) * g(2) + mrSges(2,1) * t269 + mrSges(1,2) * g(1) - mrSges(2,2) * t270 + pkin(1) * t286 + qJ(2) * t291 + t275 * t212 + t278 * t213;];
tauB = t1;
