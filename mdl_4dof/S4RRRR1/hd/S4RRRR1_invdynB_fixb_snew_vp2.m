% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:11
% EndTime: 2019-12-31 17:22:12
% DurationCPUTime: 0.72s
% Computational Cost: add. (10936->156), mult. (12133->198), div. (0->0), fcn. (6170->8), ass. (0->67)
t308 = -m(3) - m(4);
t289 = qJD(1) + qJD(2);
t285 = qJD(3) + t289;
t290 = sin(qJ(4));
t307 = t285 * t290;
t294 = cos(qJ(4));
t306 = t285 * t294;
t293 = sin(qJ(1));
t297 = cos(qJ(1));
t281 = t293 * g(1) - t297 * g(2);
t279 = qJDD(1) * pkin(1) + t281;
t282 = -t297 * g(1) - t293 * g(2);
t298 = qJD(1) ^ 2;
t280 = -t298 * pkin(1) + t282;
t292 = sin(qJ(2));
t296 = cos(qJ(2));
t264 = t296 * t279 - t292 * t280;
t288 = qJDD(1) + qJDD(2);
t262 = t288 * pkin(2) + t264;
t265 = t292 * t279 + t296 * t280;
t287 = t289 ^ 2;
t263 = -t287 * pkin(2) + t265;
t291 = sin(qJ(3));
t295 = cos(qJ(3));
t259 = t291 * t262 + t295 * t263;
t283 = t285 ^ 2;
t284 = qJDD(3) + t288;
t257 = -t283 * pkin(3) + t284 * pkin(7) + t259;
t254 = -t294 * g(3) - t290 * t257;
t271 = (-mrSges(5,1) * t294 + mrSges(5,2) * t290) * t285;
t304 = qJD(4) * t285;
t272 = t290 * t284 + t294 * t304;
t278 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t306;
t252 = m(5) * t254 + qJDD(4) * mrSges(5,1) - t272 * mrSges(5,3) + qJD(4) * t278 - t271 * t307;
t255 = -t290 * g(3) + t294 * t257;
t273 = t294 * t284 - t290 * t304;
t277 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t307;
t253 = m(5) * t255 - qJDD(4) * mrSges(5,2) + t273 * mrSges(5,3) - qJD(4) * t277 + t271 * t306;
t300 = -t290 * t252 + t294 * t253;
t243 = m(4) * t259 - t283 * mrSges(4,1) - t284 * mrSges(4,2) + t300;
t258 = t295 * t262 - t291 * t263;
t256 = -t284 * pkin(3) - t283 * pkin(7) - t258;
t299 = -m(5) * t256 + t273 * mrSges(5,1) - t272 * mrSges(5,2) - t277 * t307 + t278 * t306;
t248 = m(4) * t258 + t284 * mrSges(4,1) - t283 * mrSges(4,2) + t299;
t240 = t291 * t243 + t295 * t248;
t237 = m(3) * t264 + t288 * mrSges(3,1) - t287 * mrSges(3,2) + t240;
t301 = t295 * t243 - t291 * t248;
t238 = m(3) * t265 - t287 * mrSges(3,1) - t288 * mrSges(3,2) + t301;
t232 = t296 * t237 + t292 * t238;
t230 = m(2) * t281 + qJDD(1) * mrSges(2,1) - t298 * mrSges(2,2) + t232;
t302 = -t292 * t237 + t296 * t238;
t231 = m(2) * t282 - t298 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t302;
t305 = t297 * t230 + t293 * t231;
t244 = t294 * t252 + t290 * t253;
t303 = -t293 * t230 + t297 * t231;
t268 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t290 + Ifges(5,4) * t294) * t285;
t267 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t290 + Ifges(5,2) * t294) * t285;
t266 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t290 + Ifges(5,6) * t294) * t285;
t246 = mrSges(5,2) * t256 - mrSges(5,3) * t254 + Ifges(5,1) * t272 + Ifges(5,4) * t273 + Ifges(5,5) * qJDD(4) - qJD(4) * t267 + t266 * t306;
t245 = -mrSges(5,1) * t256 + mrSges(5,3) * t255 + Ifges(5,4) * t272 + Ifges(5,2) * t273 + Ifges(5,6) * qJDD(4) + qJD(4) * t268 - t266 * t307;
t239 = mrSges(4,1) * g(3) - mrSges(5,1) * t254 + mrSges(5,2) * t255 + mrSges(4,3) * t259 + t283 * Ifges(4,5) - Ifges(5,5) * t272 + Ifges(4,6) * t284 - Ifges(5,6) * t273 - Ifges(5,3) * qJDD(4) - pkin(3) * t244 + (-t267 * t290 + t268 * t294) * t285;
t233 = -mrSges(4,2) * g(3) - mrSges(4,3) * t258 + Ifges(4,5) * t284 - t283 * Ifges(4,6) - pkin(7) * t244 - t290 * t245 + t294 * t246;
t226 = -mrSges(3,2) * g(3) - mrSges(3,3) * t264 + Ifges(3,5) * t288 - t287 * Ifges(3,6) - pkin(6) * t240 + t295 * t233 - t291 * t239;
t225 = Ifges(3,6) * t288 + t287 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t265 + t291 * t233 + t295 * t239 - pkin(2) * (-m(4) * g(3) + t244) + pkin(6) * t301;
t224 = -mrSges(2,2) * g(3) - mrSges(2,3) * t281 + Ifges(2,5) * qJDD(1) - t298 * Ifges(2,6) - pkin(5) * t232 - t292 * t225 + t296 * t226;
t223 = Ifges(2,6) * qJDD(1) + t298 * Ifges(2,5) + mrSges(2,3) * t282 + t292 * t226 + t296 * t225 - pkin(1) * t244 + pkin(5) * t302 + (-pkin(1) * t308 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t303; -m(1) * g(2) + t305; (-m(1) - m(2) + t308) * g(3) + t244; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t305 - t293 * t223 + t297 * t224; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t303 + t297 * t223 + t293 * t224; pkin(1) * t232 - mrSges(2,2) * t282 + mrSges(2,1) * t281 + pkin(2) * t240 + mrSges(3,1) * t264 - mrSges(3,2) * t265 + pkin(7) * t300 + mrSges(4,1) * t258 - mrSges(4,2) * t259 + t290 * t246 + t294 * t245 + pkin(3) * t299 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t284 + Ifges(3,3) * t288 + Ifges(2,3) * qJDD(1);];
tauB = t1;
