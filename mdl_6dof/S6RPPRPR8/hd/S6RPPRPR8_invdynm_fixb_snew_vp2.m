% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:37:37
% EndTime: 2019-05-05 14:37:48
% DurationCPUTime: 6.13s
% Computational Cost: add. (68777->322), mult. (157263->374), div. (0->0), fcn. (103211->8), ass. (0->134)
t290 = sin(pkin(9));
t291 = cos(pkin(9));
t293 = sin(qJ(1));
t296 = cos(qJ(1));
t266 = t293 * g(1) - t296 * g(2);
t298 = qJD(1) ^ 2;
t318 = -t298 * qJ(2) + qJDD(2) - t266;
t342 = -pkin(1) - qJ(3);
t354 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t342 + t318;
t233 = t290 * g(3) + t354 * t291;
t347 = pkin(3) * t298;
t206 = (-pkin(7) * qJDD(1) - t290 * t347) * t291 + t233;
t234 = -g(3) * t291 + t354 * t290;
t279 = t290 ^ 2;
t332 = qJDD(1) * t290;
t207 = -pkin(7) * t332 - t279 * t347 + t234;
t295 = cos(qJ(4));
t348 = sin(qJ(4));
t186 = t295 * t206 - t348 * t207;
t187 = t348 * t206 + t295 * t207;
t320 = t290 * t295 + t348 * t291;
t259 = t320 * qJD(1);
t329 = t290 * t348;
t336 = qJD(1) * t291;
t260 = -qJD(1) * t329 + t295 * t336;
t215 = Ifges(5,4) * t260 - Ifges(5,2) * t259 + Ifges(5,6) * qJD(4);
t224 = -mrSges(6,2) * t259 - mrSges(6,3) * t260;
t334 = t260 * qJD(4);
t235 = qJDD(1) * t320 + t334;
t331 = qJDD(1) * t291;
t335 = qJD(4) * t259;
t236 = -qJDD(1) * t329 + t295 * t331 - t335;
t246 = mrSges(6,1) * t259 - qJD(4) * mrSges(6,3);
t248 = pkin(5) * t260 - qJD(4) * pkin(8);
t258 = t259 ^ 2;
t267 = -t296 * g(1) - t293 * g(2);
t351 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t267;
t317 = qJDD(3) + t351;
t337 = t291 ^ 2 + t279;
t211 = pkin(3) * t332 + (-pkin(7) * t337 + t342) * t298 + t317;
t349 = -2 * qJD(5);
t302 = pkin(4) * t334 + t260 * t349 + (-t236 + t335) * qJ(5) + t211;
t176 = -t258 * pkin(5) - t260 * t248 + (pkin(4) + pkin(8)) * t235 + t302;
t222 = pkin(4) * t259 - qJ(5) * t260;
t297 = qJD(4) ^ 2;
t184 = -qJDD(4) * pkin(4) - t297 * qJ(5) + t260 * t222 + qJDD(5) - t186;
t177 = (t259 * t260 - qJDD(4)) * pkin(8) + (t236 + t335) * pkin(5) + t184;
t292 = sin(qJ(6));
t294 = cos(qJ(6));
t174 = -t176 * t292 + t177 * t294;
t238 = -qJD(4) * t292 + t259 * t294;
t197 = qJD(6) * t238 + qJDD(4) * t294 + t235 * t292;
t239 = qJD(4) * t294 + t259 * t292;
t201 = -mrSges(7,1) * t238 + mrSges(7,2) * t239;
t256 = qJD(6) + t260;
t208 = -mrSges(7,2) * t256 + mrSges(7,3) * t238;
t232 = qJDD(6) + t236;
t170 = m(7) * t174 + mrSges(7,1) * t232 - mrSges(7,3) * t197 - t201 * t239 + t208 * t256;
t175 = t176 * t294 + t177 * t292;
t196 = -qJD(6) * t239 - qJDD(4) * t292 + t235 * t294;
t209 = mrSges(7,1) * t256 - mrSges(7,3) * t239;
t171 = m(7) * t175 - mrSges(7,2) * t232 + mrSges(7,3) * t196 + t201 * t238 - t209 * t256;
t159 = t294 * t170 + t292 * t171;
t313 = -t297 * pkin(4) + qJDD(4) * qJ(5) - t259 * t222 + t187;
t179 = -t235 * pkin(5) - t258 * pkin(8) + ((2 * qJD(5)) + t248) * qJD(4) + t313;
t189 = Ifges(7,5) * t239 + Ifges(7,6) * t238 + Ifges(7,3) * t256;
t191 = Ifges(7,1) * t239 + Ifges(7,4) * t238 + Ifges(7,5) * t256;
t162 = -mrSges(7,1) * t179 + mrSges(7,3) * t175 + Ifges(7,4) * t197 + Ifges(7,2) * t196 + Ifges(7,6) * t232 - t189 * t239 + t191 * t256;
t190 = Ifges(7,4) * t239 + Ifges(7,2) * t238 + Ifges(7,6) * t256;
t163 = mrSges(7,2) * t179 - mrSges(7,3) * t174 + Ifges(7,1) * t197 + Ifges(7,4) * t196 + Ifges(7,5) * t232 + t189 * t238 - t190 * t256;
t182 = qJD(4) * t349 - t313;
t212 = Ifges(6,5) * qJD(4) - Ifges(6,6) * t260 + Ifges(6,3) * t259;
t308 = -mrSges(6,2) * t184 + mrSges(6,3) * t182 - Ifges(6,1) * qJDD(4) + Ifges(6,4) * t236 - Ifges(6,5) * t235 + pkin(8) * t159 + t292 * t162 - t294 * t163 + t260 * t212;
t172 = -m(7) * t179 + t196 * mrSges(7,1) - t197 * mrSges(7,2) + t238 * t208 - t239 * t209;
t247 = mrSges(6,1) * t260 + qJD(4) * mrSges(6,2);
t310 = -m(6) * t182 + qJDD(4) * mrSges(6,3) + qJD(4) * t247 - t172;
t315 = -m(6) * t184 - t236 * mrSges(6,1) - t260 * t224 - t159;
t214 = Ifges(6,4) * qJD(4) - Ifges(6,2) * t260 + Ifges(6,6) * t259;
t339 = Ifges(5,1) * t260 - Ifges(5,4) * t259 + Ifges(5,5) * qJD(4) - t214;
t357 = t339 * t259 - mrSges(5,2) * t187 + pkin(4) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t246 + t315) + qJ(5) * (-t235 * mrSges(6,1) - t259 * t224 + t310) + mrSges(5,1) * t186 + t260 * t215 - Ifges(5,6) * t235 + Ifges(5,5) * t236 + Ifges(5,3) * qJDD(4) - t308;
t223 = mrSges(5,1) * t259 + mrSges(5,2) * t260;
t338 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t259 - t246;
t345 = mrSges(5,1) - mrSges(6,2);
t155 = m(5) * t186 - t236 * mrSges(5,3) + qJD(4) * t338 + qJDD(4) * t345 - t260 * t223 + t315;
t245 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t260;
t166 = m(5) * t187 - qJDD(4) * mrSges(5,2) - qJD(4) * t245 + (-t223 - t224) * t259 + (-mrSges(5,3) - mrSges(6,1)) * t235 + t310;
t150 = t295 * t155 + t348 * t166;
t323 = Ifges(4,4) * t291 - Ifges(4,2) * t290;
t324 = Ifges(4,1) * t291 - Ifges(4,4) * t290;
t356 = mrSges(4,1) * t233 - mrSges(4,2) * t234 + Ifges(4,5) * t331 + pkin(3) * t150 + (t290 * t324 + t291 * t323) * t298 + t357;
t327 = t337 * mrSges(4,3);
t355 = t298 * t327;
t242 = t298 * t342 + t317;
t160 = -t292 * t170 + t294 * t171;
t181 = t235 * pkin(4) + t302;
t319 = -m(6) * t181 + t236 * mrSges(6,3) + t260 * t247 - t160;
t304 = m(5) * t211 + t236 * mrSges(5,2) + t345 * t235 + t260 * t245 + t338 * t259 - t319;
t353 = -m(4) * t242 - mrSges(4,1) * t332 - mrSges(4,2) * t331 - t304;
t322 = -qJDD(1) * mrSges(4,3) - t298 * (mrSges(4,1) * t290 + mrSges(4,2) * t291);
t147 = m(4) * t233 + t291 * t322 + t150;
t325 = -t348 * t155 + t295 * t166;
t148 = m(4) * t234 + t290 * t322 + t325;
t143 = t291 * t147 + t290 * t148;
t257 = -qJDD(1) * pkin(1) + t318;
t352 = mrSges(3,1) * t257 + pkin(2) * t143 + t356;
t346 = mrSges(2,1) - mrSges(3,2);
t344 = Ifges(5,4) + Ifges(6,6);
t343 = -Ifges(2,6) + Ifges(3,5);
t341 = Ifges(4,6) * t290;
t216 = Ifges(6,1) * qJD(4) - Ifges(6,4) * t260 + Ifges(6,5) * t259;
t340 = -Ifges(5,5) * t260 + Ifges(5,6) * t259 - Ifges(5,3) * qJD(4) - t216;
t328 = Ifges(3,4) + t341;
t144 = -t147 * t290 + t291 * t148;
t316 = -m(3) * t257 + t298 * mrSges(3,3) - t143;
t312 = mrSges(7,1) * t174 - mrSges(7,2) * t175 + Ifges(7,5) * t197 + Ifges(7,6) * t196 + Ifges(7,3) * t232 + t239 * t190 - t238 * t191;
t156 = -t235 * mrSges(6,2) - t259 * t246 - t319;
t307 = -mrSges(6,1) * t182 + mrSges(6,2) * t181 - pkin(5) * t172 - pkin(8) * t160 - t294 * t162 - t292 * t163;
t139 = -mrSges(5,1) * t211 + mrSges(5,3) * t187 - pkin(4) * t156 + t340 * t260 + t344 * t236 + (-Ifges(5,2) - Ifges(6,3)) * t235 + (Ifges(5,6) - Ifges(6,5)) * qJDD(4) + t339 * qJD(4) + t307;
t305 = mrSges(6,1) * t184 - mrSges(6,3) * t181 + pkin(5) * t159 + t312;
t145 = t340 * t259 + (Ifges(5,1) + Ifges(6,2)) * t236 - t344 * t235 + (Ifges(5,5) - Ifges(6,4)) * qJDD(4) + (-t215 + t212) * qJD(4) + mrSges(5,2) * t211 - mrSges(5,3) * t186 - qJ(5) * t156 + t305;
t262 = (Ifges(4,5) * t291 - t341) * qJD(1);
t136 = -mrSges(4,1) * t242 + mrSges(4,3) * t234 - pkin(3) * t304 + pkin(7) * t325 + qJDD(1) * t323 + t295 * t139 + t145 * t348 - t262 * t336;
t138 = -t290 * qJD(1) * t262 + mrSges(4,2) * t242 - mrSges(4,3) * t233 - pkin(7) * t150 + qJDD(1) * t324 - t139 * t348 + t295 * t145;
t250 = t298 * pkin(1) - t351;
t311 = mrSges(3,2) * t257 - mrSges(3,3) * t250 + Ifges(3,1) * qJDD(1) - qJ(3) * t143 - t136 * t290 + t291 * t138;
t309 = -mrSges(3,1) * t250 - pkin(2) * (t353 + t355) - qJ(3) * t144 - t291 * t136 - t290 * t138;
t301 = -m(3) * t250 + t298 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t353;
t306 = -mrSges(2,2) * t267 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t316) + qJ(2) * (t301 - t355) + mrSges(2,1) * t266 + Ifges(2,3) * qJDD(1) + t311;
t151 = -qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t327) * t298 + m(2) * t267 + t301;
t142 = -m(3) * g(3) + t144;
t140 = m(2) * t266 - t298 * mrSges(2,2) + qJDD(1) * t346 + t316;
t135 = t343 * t298 + (Ifges(2,5) - t328) * qJDD(1) - mrSges(2,3) * t266 - qJ(2) * t142 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t352;
t134 = mrSges(2,3) * t267 - pkin(1) * t142 + (-Ifges(3,4) + Ifges(2,5)) * t298 - t343 * qJDD(1) + t346 * g(3) + t309;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t296 * t135 - t293 * t134 - pkin(6) * (t140 * t296 + t151 * t293), t135, t311, t138, t145, -t259 * t214 - t308, t163; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t293 * t135 + t296 * t134 + pkin(6) * (-t140 * t293 + t151 * t296), t134, -mrSges(3,3) * g(3) - t298 * Ifges(3,5) + qJDD(1) * t328 - t352, t136, t139, Ifges(6,4) * qJDD(4) - Ifges(6,2) * t236 + Ifges(6,6) * t235 - qJD(4) * t212 + t259 * t216 - t305, t162; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t306, t306, mrSges(3,2) * g(3) + t298 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t309, -Ifges(4,6) * t332 + t356, t357, Ifges(6,5) * qJDD(4) - Ifges(6,6) * t236 + Ifges(6,3) * t235 + qJD(4) * t214 + t260 * t216 - t307, t312;];
m_new  = t1;
