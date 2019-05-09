% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 08:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:22:15
% EndTime: 2019-05-05 08:22:39
% DurationCPUTime: 11.38s
% Computational Cost: add. (204536->338), mult. (390305->416), div. (0->0), fcn. (265670->12), ass. (0->138)
t309 = sin(pkin(11));
t311 = cos(pkin(11));
t296 = t309 * g(1) - t311 * g(2);
t297 = -t311 * g(1) - t309 * g(2);
t308 = -g(3) + qJDD(1);
t319 = cos(qJ(2));
t312 = cos(pkin(6));
t316 = sin(qJ(2));
t345 = t312 * t316;
t310 = sin(pkin(6));
t346 = t310 * t316;
t239 = t296 * t345 + t319 * t297 + t308 * t346;
t321 = qJD(2) ^ 2;
t233 = -t321 * pkin(2) + qJDD(2) * pkin(8) + t239;
t270 = -t310 * t296 + t312 * t308;
t315 = sin(qJ(3));
t318 = cos(qJ(3));
t226 = t318 * t233 + t315 * t270;
t292 = (-pkin(3) * t318 - pkin(9) * t315) * qJD(2);
t320 = qJD(3) ^ 2;
t340 = t318 * qJD(2);
t217 = -t320 * pkin(3) + qJDD(3) * pkin(9) + t292 * t340 + t226;
t238 = -t316 * t297 + (t296 * t312 + t308 * t310) * t319;
t232 = -qJDD(2) * pkin(2) - t321 * pkin(8) - t238;
t339 = qJD(2) * qJD(3);
t337 = t318 * t339;
t293 = t315 * qJDD(2) + t337;
t338 = t315 * t339;
t294 = t318 * qJDD(2) - t338;
t219 = (-t293 - t337) * pkin(9) + (-t294 + t338) * pkin(3) + t232;
t314 = sin(qJ(4));
t352 = cos(qJ(4));
t205 = t352 * t217 + t314 * t219;
t341 = qJD(2) * t315;
t290 = t314 * qJD(3) + t341 * t352;
t254 = t290 * qJD(4) - qJDD(3) * t352 + t314 * t293;
t304 = qJD(4) - t340;
t266 = t304 * mrSges(5,1) - t290 * mrSges(5,3);
t286 = qJDD(4) - t294;
t289 = -qJD(3) * t352 + t314 * t341;
t204 = -t314 * t217 + t219 * t352;
t260 = t289 * pkin(4) - t290 * qJ(5);
t303 = t304 ^ 2;
t202 = -t286 * pkin(4) - t303 * qJ(5) + t290 * t260 + qJDD(5) - t204;
t255 = -t289 * qJD(4) + t314 * qJDD(3) + t293 * t352;
t347 = t289 * t304;
t194 = (-t255 - t347) * pkin(10) + (t289 * t290 - t286) * pkin(5) + t202;
t353 = 2 * qJD(5);
t200 = -t303 * pkin(4) + t286 * qJ(5) - t289 * t260 + t304 * t353 + t205;
t269 = -t304 * pkin(5) - t290 * pkin(10);
t285 = t289 ^ 2;
t195 = -t285 * pkin(5) + t254 * pkin(10) + t304 * t269 + t200;
t313 = sin(qJ(6));
t317 = cos(qJ(6));
t191 = t317 * t194 - t313 * t195;
t256 = t317 * t289 - t313 * t290;
t213 = t256 * qJD(6) + t313 * t254 + t317 * t255;
t257 = t313 * t289 + t317 * t290;
t227 = -t256 * mrSges(7,1) + t257 * mrSges(7,2);
t302 = qJD(6) - t304;
t234 = -t302 * mrSges(7,2) + t256 * mrSges(7,3);
t282 = qJDD(6) - t286;
t187 = m(7) * t191 + t282 * mrSges(7,1) - t213 * mrSges(7,3) - t257 * t227 + t302 * t234;
t192 = t313 * t194 + t317 * t195;
t212 = -t257 * qJD(6) + t317 * t254 - t313 * t255;
t235 = t302 * mrSges(7,1) - t257 * mrSges(7,3);
t188 = m(7) * t192 - t282 * mrSges(7,2) + t212 * mrSges(7,3) + t256 * t227 - t302 * t235;
t178 = -t313 * t187 + t317 * t188;
t267 = -t304 * mrSges(6,1) + t290 * mrSges(6,2);
t332 = m(6) * t200 + t286 * mrSges(6,3) + t304 * t267 + t178;
t261 = t289 * mrSges(6,1) - t290 * mrSges(6,3);
t342 = -t289 * mrSges(5,1) - t290 * mrSges(5,2) - t261;
t350 = -mrSges(5,3) - mrSges(6,2);
t173 = m(5) * t205 - t286 * mrSges(5,2) + t254 * t350 - t304 * t266 + t289 * t342 + t332;
t265 = -t304 * mrSges(5,2) - t289 * mrSges(5,3);
t177 = t317 * t187 + t313 * t188;
t268 = -t289 * mrSges(6,2) + t304 * mrSges(6,3);
t328 = -m(6) * t202 + t286 * mrSges(6,1) + t304 * t268 - t177;
t174 = m(5) * t204 + t286 * mrSges(5,1) + t255 * t350 + t304 * t265 + t290 * t342 + t328;
t171 = t352 * t173 - t314 * t174;
t291 = (-mrSges(4,1) * t318 + mrSges(4,2) * t315) * qJD(2);
t298 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t341;
t169 = m(4) * t226 - qJDD(3) * mrSges(4,2) + t294 * mrSges(4,3) - qJD(3) * t298 + t291 * t340 + t171;
t225 = -t315 * t233 + t318 * t270;
t330 = qJDD(3) * pkin(3) + t320 * pkin(9) - t292 * t341 + t225;
t356 = (-t255 + t347) * qJ(5) - t330;
t197 = -t285 * pkin(10) + (-pkin(4) - pkin(5)) * t254 + (-pkin(4) * t304 + t269 + t353) * t290 - t356;
t193 = -m(7) * t197 + t212 * mrSges(7,1) - t213 * mrSges(7,2) + t256 * t234 - t257 * t235;
t203 = -0.2e1 * qJD(5) * t290 + (t290 * t304 + t254) * pkin(4) + t356;
t185 = m(6) * t203 + t254 * mrSges(6,1) - t255 * mrSges(6,3) - t290 * t267 + t289 * t268 + t193;
t184 = m(5) * t330 - t254 * mrSges(5,1) - t255 * mrSges(5,2) - t289 * t265 - t290 * t266 - t185;
t299 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t340;
t183 = m(4) * t225 + qJDD(3) * mrSges(4,1) - t293 * mrSges(4,3) + qJD(3) * t299 - t291 * t341 + t184;
t164 = t315 * t169 + t318 * t183;
t244 = Ifges(6,1) * t290 + Ifges(6,4) * t304 + Ifges(6,5) * t289;
t245 = Ifges(5,1) * t290 - Ifges(5,4) * t289 + Ifges(5,5) * t304;
t220 = Ifges(7,5) * t257 + Ifges(7,6) * t256 + Ifges(7,3) * t302;
t222 = Ifges(7,1) * t257 + Ifges(7,4) * t256 + Ifges(7,5) * t302;
t180 = -mrSges(7,1) * t197 + mrSges(7,3) * t192 + Ifges(7,4) * t213 + Ifges(7,2) * t212 + Ifges(7,6) * t282 - t257 * t220 + t302 * t222;
t221 = Ifges(7,4) * t257 + Ifges(7,2) * t256 + Ifges(7,6) * t302;
t181 = mrSges(7,2) * t197 - mrSges(7,3) * t191 + Ifges(7,1) * t213 + Ifges(7,4) * t212 + Ifges(7,5) * t282 + t256 * t220 - t302 * t221;
t326 = -mrSges(6,1) * t203 + mrSges(6,2) * t200 - pkin(5) * t193 - pkin(10) * t178 - t317 * t180 - t313 * t181;
t242 = Ifges(6,4) * t290 + Ifges(6,2) * t304 + Ifges(6,6) * t289;
t344 = -Ifges(5,5) * t290 + Ifges(5,6) * t289 - Ifges(5,3) * t304 - t242;
t158 = mrSges(5,1) * t330 + mrSges(5,3) * t205 - pkin(4) * t185 + (t245 + t244) * t304 + t344 * t290 + (Ifges(5,6) - Ifges(6,6)) * t286 + (Ifges(5,4) - Ifges(6,5)) * t255 + (-Ifges(5,2) - Ifges(6,3)) * t254 + t326;
t243 = Ifges(5,4) * t290 - Ifges(5,2) * t289 + Ifges(5,6) * t304;
t240 = Ifges(6,5) * t290 + Ifges(6,6) * t304 + Ifges(6,3) * t289;
t327 = mrSges(6,2) * t202 - mrSges(6,3) * t203 + Ifges(6,1) * t255 + Ifges(6,4) * t286 + Ifges(6,5) * t254 - pkin(10) * t177 - t313 * t180 + t317 * t181 + t304 * t240;
t159 = -mrSges(5,2) * t330 - mrSges(5,3) * t204 + Ifges(5,1) * t255 - Ifges(5,4) * t254 + Ifges(5,5) * t286 - qJ(5) * t185 - t304 * t243 + t289 * t344 + t327;
t275 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t315 + Ifges(4,2) * t318) * qJD(2);
t276 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t315 + Ifges(4,4) * t318) * qJD(2);
t354 = mrSges(4,1) * t225 - mrSges(4,2) * t226 + Ifges(4,5) * t293 + Ifges(4,6) * t294 + Ifges(4,3) * qJDD(3) + pkin(3) * t184 + pkin(9) * t171 + (t275 * t315 - t276 * t318) * qJD(2) + t158 * t352 + t314 * t159;
t148 = -mrSges(3,1) * t270 + mrSges(3,3) * t239 + t321 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t164 - t354;
t336 = t318 * t169 - t315 * t183;
t162 = m(3) * t239 - t321 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t336;
t170 = t314 * t173 + t174 * t352;
t325 = -m(4) * t232 + t294 * mrSges(4,1) - t293 * mrSges(4,2) - t298 * t341 + t299 * t340 - t170;
t166 = m(3) * t238 + qJDD(2) * mrSges(3,1) - t321 * mrSges(3,2) + t325;
t157 = t319 * t162 - t316 * t166;
t357 = pkin(7) * t157 + t148 * t319;
t331 = mrSges(7,1) * t191 - mrSges(7,2) * t192 + Ifges(7,5) * t213 + Ifges(7,6) * t212 + Ifges(7,3) * t282 + t257 * t221 - t256 * t222;
t324 = mrSges(6,1) * t202 - mrSges(6,3) * t200 - Ifges(6,4) * t255 - Ifges(6,2) * t286 - Ifges(6,6) * t254 + pkin(5) * t177 - t289 * t244 + t331;
t355 = (t243 - t240) * t290 + mrSges(5,1) * t204 - mrSges(5,2) * t205 + Ifges(5,5) * t255 - Ifges(5,6) * t254 + Ifges(5,3) * t286 + pkin(4) * (-t255 * mrSges(6,2) - t290 * t261 + t328) + qJ(5) * (-t254 * mrSges(6,2) - t289 * t261 + t332) + t289 * t245 - t324;
t348 = t166 * t319;
t163 = m(3) * t270 + t164;
t153 = t162 * t345 - t310 * t163 + t312 * t348;
t274 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t315 + Ifges(4,6) * t318) * qJD(2);
t149 = mrSges(4,2) * t232 - mrSges(4,3) * t225 + Ifges(4,1) * t293 + Ifges(4,4) * t294 + Ifges(4,5) * qJDD(3) - pkin(9) * t170 - qJD(3) * t275 - t314 * t158 + t159 * t352 + t274 * t340;
t154 = -mrSges(4,1) * t232 + mrSges(4,3) * t226 + Ifges(4,4) * t293 + Ifges(4,2) * t294 + Ifges(4,6) * qJDD(3) - pkin(3) * t170 + qJD(3) * t276 - t274 * t341 - t355;
t144 = mrSges(3,1) * t238 - mrSges(3,2) * t239 + Ifges(3,3) * qJDD(2) + pkin(2) * t325 + pkin(8) * t336 + t315 * t149 + t318 * t154;
t146 = mrSges(3,2) * t270 - mrSges(3,3) * t238 + Ifges(3,5) * qJDD(2) - t321 * Ifges(3,6) - pkin(8) * t164 + t318 * t149 - t315 * t154;
t329 = mrSges(2,1) * t296 - mrSges(2,2) * t297 + pkin(1) * t153 + t312 * t144 + t146 * t346 + t357 * t310;
t155 = m(2) * t297 + t157;
t152 = t312 * t163 + (t162 * t316 + t348) * t310;
t150 = m(2) * t296 + t153;
t142 = mrSges(2,2) * t308 - mrSges(2,3) * t296 + t319 * t146 - t316 * t148 + (-t152 * t310 - t153 * t312) * pkin(7);
t141 = -mrSges(2,1) * t308 + mrSges(2,3) * t297 - pkin(1) * t152 - t310 * t144 + (t146 * t316 + t357) * t312;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t311 * t142 - t309 * t141 - qJ(1) * (t311 * t150 + t309 * t155), t142, t146, t149, t159, -t289 * t242 + t327, t181; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t309 * t142 + t311 * t141 + qJ(1) * (-t309 * t150 + t311 * t155), t141, t148, t154, t158, -t290 * t240 - t324, t180; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t329, t329, t144, t354, t355, Ifges(6,5) * t255 + Ifges(6,6) * t286 + Ifges(6,3) * t254 + t290 * t242 - t304 * t244 - t326, t331;];
m_new  = t1;
