% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR9
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:44
% EndTime: 2019-12-31 19:40:49
% DurationCPUTime: 2.32s
% Computational Cost: add. (18116->319), mult. (37571->371), div. (0->0), fcn. (17431->6), ass. (0->121)
t301 = sin(qJ(1));
t304 = cos(qJ(1));
t276 = -t304 * g(1) - t301 * g(2);
t306 = qJD(1) ^ 2;
t231 = -t306 * pkin(1) + qJDD(1) * pkin(6) + t276;
t300 = sin(qJ(2));
t303 = cos(qJ(2));
t207 = -t303 * g(3) - t300 * t231;
t208 = -t300 * g(3) + t303 * t231;
t222 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t300 - Ifges(4,3) * t303) * qJD(1);
t226 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t300 + Ifges(3,2) * t303) * qJD(1);
t258 = (-mrSges(4,1) * t303 - mrSges(4,3) * t300) * qJD(1);
t332 = qJD(1) * qJD(2);
t327 = t303 * t332;
t262 = t300 * qJDD(1) + t327;
t328 = t300 * t332;
t263 = t303 * qJDD(1) - t328;
t334 = qJD(1) * t303;
t272 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t334;
t305 = qJD(2) ^ 2;
t257 = (-pkin(2) * t303 - qJ(3) * t300) * qJD(1);
t333 = t300 * qJD(1);
t325 = t257 * t333 + qJDD(3) - t207;
t194 = -qJDD(2) * pkin(2) - t305 * qJ(3) + t325;
t339 = t303 * t306;
t331 = qJD(1) * qJD(4);
t354 = -0.2e1 * t300 * t331 + (-t262 + t327) * qJ(4);
t187 = (-t300 * t339 - qJDD(2)) * pkin(3) + t194 + t354;
t260 = (mrSges(5,1) * t300 - mrSges(5,2) * t303) * qJD(1);
t268 = -qJD(2) * pkin(3) - qJ(4) * t333;
t275 = t301 * g(1) - t304 * g(2);
t230 = -qJDD(1) * pkin(1) - t306 * pkin(6) - t275;
t323 = -t263 * pkin(2) + t230 + (-t262 - t327) * qJ(3);
t340 = t303 ^ 2 * t306;
t350 = 2 * qJD(3);
t314 = -qJ(4) * t340 + qJDD(4) - t323 + (t268 + t350) * t333;
t348 = pkin(3) + pkin(7);
t349 = -pkin(2) - pkin(7);
t176 = (pkin(4) * t303 + t300 * t349) * t332 + t348 * t263 + t314 + t262 * pkin(4);
t261 = (pkin(4) * t300 + pkin(7) * t303) * qJD(1);
t179 = (-pkin(4) - qJ(3)) * t305 + (-pkin(3) * t339 - qJD(1) * t261) * t300 + (-pkin(2) - t348) * qJDD(2) + t325 + t354;
t299 = sin(qJ(5));
t302 = cos(qJ(5));
t174 = t302 * t176 - t299 * t179;
t255 = -t302 * qJD(2) + t299 * t334;
t203 = t255 * qJD(5) - t299 * qJDD(2) - t302 * t263;
t256 = -t299 * qJD(2) - t302 * t334;
t204 = -t255 * mrSges(6,1) + t256 * mrSges(6,2);
t279 = qJD(5) + t333;
t205 = -t279 * mrSges(6,2) + t255 * mrSges(6,3);
t253 = qJDD(5) + t262;
t169 = m(6) * t174 + t253 * mrSges(6,1) - t203 * mrSges(6,3) - t256 * t204 + t279 * t205;
t175 = t299 * t176 + t302 * t179;
t202 = -t256 * qJD(5) - t302 * qJDD(2) + t299 * t263;
t206 = t279 * mrSges(6,1) - t256 * mrSges(6,3);
t170 = m(6) * t175 - t253 * mrSges(6,2) + t202 * mrSges(6,3) + t255 * t204 - t279 * t206;
t338 = -t299 * t169 + t302 * t170;
t324 = -m(5) * t187 + t260 * t333 - t338;
t319 = -qJDD(2) * mrSges(5,2) - qJD(2) * t272 + t324;
t151 = -t262 * mrSges(5,3) - t319;
t283 = qJD(2) * t350;
t346 = t305 * pkin(2);
t356 = qJDD(2) * qJ(3) + t257 * t334 + t208;
t192 = t283 - t346 + t356;
t317 = pkin(3) * t340 + t263 * qJ(4) - t356;
t178 = qJDD(2) * pkin(4) + qJD(2) * t268 + t283 + t349 * t305 + (-0.2e1 * qJD(4) - t261) * t334 - t317;
t195 = Ifges(6,5) * t256 + Ifges(6,6) * t255 + Ifges(6,3) * t279;
t197 = Ifges(6,1) * t256 + Ifges(6,4) * t255 + Ifges(6,5) * t279;
t161 = -mrSges(6,1) * t178 + mrSges(6,3) * t175 + Ifges(6,4) * t203 + Ifges(6,2) * t202 + Ifges(6,6) * t253 - t256 * t195 + t279 * t197;
t196 = Ifges(6,4) * t256 + Ifges(6,2) * t255 + Ifges(6,6) * t279;
t162 = mrSges(6,2) * t178 - mrSges(6,3) * t174 + Ifges(6,1) * t203 + Ifges(6,4) * t202 + Ifges(6,5) * t253 + t255 * t195 - t279 * t196;
t351 = -2 * qJD(3);
t186 = 0.2e1 * t303 * t331 + t346 + (t351 - t268) * qJD(2) + t317;
t224 = -Ifges(5,6) * qJD(2) + (-Ifges(5,4) * t303 - Ifges(5,2) * t300) * qJD(1);
t227 = -Ifges(5,5) * qJD(2) + (-Ifges(5,1) * t303 - Ifges(5,4) * t300) * qJD(1);
t322 = -m(6) * t178 + t202 * mrSges(6,1) - t203 * mrSges(6,2) + t255 * t205 - t256 * t206;
t318 = mrSges(5,1) * t186 - mrSges(5,2) * t187 - Ifges(5,5) * t263 - Ifges(5,6) * t262 - Ifges(5,3) * qJDD(2) + pkin(4) * t322 + pkin(7) * t338 + t302 * t161 + t299 * t162 - t224 * t334 + t227 * t333;
t311 = -mrSges(4,1) * t194 + mrSges(4,3) * t192 + Ifges(4,4) * t262 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t263 - pkin(3) * t151 - t318;
t271 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t333;
t269 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t333;
t315 = -m(5) * t186 + qJDD(2) * mrSges(5,1) - t263 * mrSges(5,3) + qJD(2) * t269 - t322;
t312 = m(4) * t192 + qJDD(2) * mrSges(4,3) + qJD(2) * t271 + t258 * t334 + t315;
t329 = t260 * t334;
t228 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t300 - Ifges(4,5) * t303) * qJD(1);
t335 = t228 + Ifges(3,5) * qJD(2) + (Ifges(3,1) * t300 + Ifges(3,4) * t303) * qJD(1);
t274 = mrSges(4,2) * t334 + qJD(2) * mrSges(4,3);
t355 = -m(4) * t194 + qJDD(2) * mrSges(4,1) + qJD(2) * t274;
t357 = -((t222 - t226) * t300 + t303 * t335) * qJD(1) + mrSges(3,1) * t207 - mrSges(3,2) * t208 + Ifges(3,5) * t262 + Ifges(3,6) * t263 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t258 * t333 + (-mrSges(4,2) + mrSges(5,3)) * t262 + t319 + t355) + qJ(3) * (t263 * mrSges(4,2) + t312 - t329) + t311;
t221 = -Ifges(5,3) * qJD(2) + (-Ifges(5,5) * t303 - Ifges(5,6) * t300) * qJD(1);
t353 = Ifges(5,4) * t263 + Ifges(5,2) * t262 - t221 * t334;
t345 = mrSges(3,3) + mrSges(4,2);
t344 = Ifges(4,6) - Ifges(5,5);
t154 = t302 * t169 + t299 * t170;
t225 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t300 - Ifges(4,6) * t303) * qJD(1);
t337 = -t221 + t225;
t259 = (-mrSges(3,1) * t303 + mrSges(3,2) * t300) * qJD(1);
t273 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t334;
t147 = m(3) * t207 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t272 + t273) * qJD(2) + (-t258 - t259) * t333 + (mrSges(5,3) - t345) * t262 + t324 + t355;
t270 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t333;
t157 = (t259 - t260) * t334 + t345 * t263 - qJDD(2) * mrSges(3,2) + t312 - qJD(2) * t270 + m(3) * t208;
t326 = -t300 * t147 + t303 * t157;
t183 = -pkin(2) * t328 + t263 * pkin(3) + t314;
t150 = -m(5) * t183 - t262 * mrSges(5,1) + t263 * mrSges(5,2) - t269 * t333 + t272 * t334 - t154;
t188 = (pkin(2) * qJD(2) + t351) * t333 + t323;
t148 = m(4) * t188 - t263 * mrSges(4,1) - t262 * mrSges(4,3) - t271 * t333 - t274 * t334 + t150;
t223 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t300 + Ifges(3,6) * t303) * qJD(1);
t316 = -mrSges(5,2) * t183 + mrSges(5,3) * t186 + Ifges(5,1) * t263 + Ifges(5,4) * t262 + pkin(7) * t154 - qJD(2) * t224 + t299 * t161 - t302 * t162;
t310 = mrSges(4,1) * t188 - mrSges(4,2) * t192 + pkin(3) * t150 + qJ(4) * (t315 - t329) - t316;
t139 = (-t223 - t337) * t333 - t310 + (Ifges(4,3) + Ifges(3,2)) * t263 + (Ifges(3,4) - Ifges(4,5)) * t262 + (Ifges(3,6) - t344) * qJDD(2) + t335 * qJD(2) + mrSges(3,3) * t208 - mrSges(3,1) * t230 - pkin(2) * t148;
t320 = mrSges(6,1) * t174 - mrSges(6,2) * t175 + Ifges(6,5) * t203 + Ifges(6,6) * t202 + Ifges(6,3) * t253 + t256 * t196 - t255 * t197;
t313 = mrSges(5,1) * t183 - mrSges(5,3) * t187 + Ifges(5,6) * qJDD(2) + pkin(4) * t154 + qJD(2) * t227 + t320;
t308 = mrSges(4,2) * t194 - mrSges(4,3) * t188 + Ifges(4,1) * t262 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t263 - qJ(4) * t151 + qJD(2) * t222 + t225 * t334 + t313;
t141 = (-t221 + t223) * t334 + Ifges(3,5) * qJDD(2) + t308 + (Ifges(5,4) + Ifges(3,4)) * t263 + (Ifges(3,1) + Ifges(5,2)) * t262 - qJD(2) * t226 + mrSges(3,2) * t230 - mrSges(3,3) * t207 - qJ(3) * t148;
t309 = -m(3) * t230 + t263 * mrSges(3,1) - t262 * mrSges(3,2) - t270 * t333 + t273 * t334 - t148;
t321 = mrSges(2,1) * t275 - mrSges(2,2) * t276 + Ifges(2,3) * qJDD(1) + pkin(1) * t309 + pkin(6) * t326 + t303 * t139 + t300 * t141;
t145 = m(2) * t275 + qJDD(1) * mrSges(2,1) - t306 * mrSges(2,2) + t309;
t144 = t303 * t147 + t300 * t157;
t142 = m(2) * t276 - t306 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t326;
t137 = mrSges(2,1) * g(3) + mrSges(2,3) * t276 + t306 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t144 - t357;
t136 = -mrSges(2,2) * g(3) - mrSges(2,3) * t275 + Ifges(2,5) * qJDD(1) - t306 * Ifges(2,6) - pkin(6) * t144 - t300 * t139 + t303 * t141;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t136 - t301 * t137 - pkin(5) * (t301 * t142 + t304 * t145), t136, t141, t308 + t353, -Ifges(5,5) * qJDD(2) - t221 * t333 - t316, t162; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t136 + t304 * t137 + pkin(5) * (t304 * t142 - t301 * t145), t137, t139, (-t300 * t222 - t303 * t228) * qJD(1) + t311, -t313 - t353, t161; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t321, t321, t357, Ifges(4,5) * t262 - Ifges(4,3) * t263 - qJD(2) * t228 + qJDD(2) * t344 + t333 * t337 + t310, t318, t320;];
m_new = t1;
