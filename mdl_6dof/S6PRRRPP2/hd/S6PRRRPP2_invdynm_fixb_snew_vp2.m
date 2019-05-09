% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRPP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:42:07
% EndTime: 2019-05-05 06:42:24
% DurationCPUTime: 7.62s
% Computational Cost: add. (102492->337), mult. (193801->398), div. (0->0), fcn. (125828->10), ass. (0->129)
t308 = sin(pkin(10));
t310 = cos(pkin(10));
t295 = t308 * g(1) - t310 * g(2);
t296 = -t310 * g(1) - t308 * g(2);
t307 = -g(3) + qJDD(1);
t316 = cos(qJ(2));
t311 = cos(pkin(6));
t314 = sin(qJ(2));
t343 = t311 * t314;
t309 = sin(pkin(6));
t344 = t309 * t314;
t221 = t295 * t343 + t316 * t296 + t307 * t344;
t318 = qJD(2) ^ 2;
t215 = -t318 * pkin(2) + qJDD(2) * pkin(8) + t221;
t268 = -t309 * t295 + t311 * t307;
t313 = sin(qJ(3));
t315 = cos(qJ(3));
t205 = t315 * t215 + t313 * t268;
t291 = (-pkin(3) * t315 - pkin(9) * t313) * qJD(2);
t317 = qJD(3) ^ 2;
t338 = t315 * qJD(2);
t201 = -t317 * pkin(3) + qJDD(3) * pkin(9) + t291 * t338 + t205;
t220 = -t314 * t296 + (t295 * t311 + t307 * t309) * t316;
t214 = -qJDD(2) * pkin(2) - t318 * pkin(8) - t220;
t337 = qJD(2) * qJD(3);
t335 = t315 * t337;
t292 = t313 * qJDD(2) + t335;
t304 = t313 * t337;
t293 = t315 * qJDD(2) - t304;
t203 = (-t292 - t335) * pkin(9) + (-t293 + t304) * pkin(3) + t214;
t312 = sin(qJ(4));
t352 = cos(qJ(4));
t196 = -t312 * t201 + t203 * t352;
t339 = qJD(2) * t313;
t288 = -qJD(3) * t352 + t312 * t339;
t252 = -t288 * qJD(4) + t312 * qJDD(3) + t292 * t352;
t302 = -qJD(4) + t338;
t261 = -t302 * mrSges(7,2) + t288 * mrSges(7,3);
t262 = t302 * mrSges(5,2) - t288 * mrSges(5,3);
t285 = -qJDD(4) + t293;
t289 = t312 * qJD(3) + t339 * t352;
t255 = t288 * pkin(4) - t289 * qJ(5);
t301 = t302 ^ 2;
t194 = t285 * pkin(4) - t301 * qJ(5) + t289 * t255 + qJDD(5) - t196;
t346 = t288 * t302;
t182 = -0.2e1 * qJD(6) * t289 + (-t252 + t346) * qJ(6) + (t288 * t289 + t285) * pkin(5) + t194;
t257 = -t288 * mrSges(7,1) + t289 * mrSges(7,2);
t332 = -m(7) * t182 + t252 * mrSges(7,3) + t289 * t257;
t256 = t288 * mrSges(6,1) - t289 * mrSges(6,3);
t340 = -t288 * mrSges(5,1) - t289 * mrSges(5,2) - t256;
t349 = -mrSges(5,3) - mrSges(6,2);
t267 = -t288 * mrSges(6,2) - t302 * mrSges(6,3);
t357 = -m(6) * t194 - t285 * mrSges(6,1) - t302 * t267;
t172 = m(5) * t196 + (-t261 - t262) * t302 + t340 * t289 + (-mrSges(5,1) - mrSges(7,1)) * t285 + t349 * t252 + t332 + t357;
t197 = t352 * t201 + t312 * t203;
t251 = t289 * qJD(4) - qJDD(3) * t352 + t312 * t292;
t264 = t302 * mrSges(7,1) - t289 * mrSges(7,3);
t265 = -t302 * mrSges(5,1) - t289 * mrSges(5,3);
t353 = -2 * qJD(5);
t192 = -t301 * pkin(4) - t285 * qJ(5) - t288 * t255 + t302 * t353 + t197;
t263 = t302 * pkin(5) - t289 * qJ(6);
t284 = t288 ^ 2;
t186 = -t284 * pkin(5) + t251 * qJ(6) + 0.2e1 * qJD(6) * t288 - t302 * t263 + t192;
t336 = m(7) * t186 + t251 * mrSges(7,3) + t288 * t257;
t266 = t302 * mrSges(6,1) + t289 * mrSges(6,2);
t358 = m(6) * t192 - t285 * mrSges(6,3) - t302 * t266;
t173 = m(5) * t197 + (-t264 + t265) * t302 + t340 * t288 + (mrSges(5,2) - mrSges(7,2)) * t285 + t349 * t251 + t336 + t358;
t168 = -t312 * t172 + t352 * t173;
t290 = (-mrSges(4,1) * t315 + mrSges(4,2) * t313) * qJD(2);
t297 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t339;
t166 = m(4) * t205 - qJDD(3) * mrSges(4,2) + t293 * mrSges(4,3) - qJD(3) * t297 + t290 * t338 + t168;
t204 = -t313 * t215 + t315 * t268;
t326 = qJDD(3) * pkin(3) + t317 * pkin(9) - t291 * t339 + t204;
t356 = -(t252 + t346) * qJ(5) - t326;
t189 = -t284 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t251 + (pkin(4) * t302 + (2 * qJD(5)) + t263) * t289 - t356;
t180 = -m(7) * t189 + t251 * mrSges(7,1) - t252 * mrSges(7,2) + t288 * t261 - t289 * t264;
t195 = t289 * t353 + (-t289 * t302 + t251) * pkin(4) + t356;
t177 = m(6) * t195 + t251 * mrSges(6,1) - t252 * mrSges(6,3) - t289 * t266 + t288 * t267 + t180;
t174 = m(5) * t326 - t251 * mrSges(5,1) - t252 * mrSges(5,2) - t288 * t262 - t289 * t265 - t177;
t298 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t338;
t170 = m(4) * t204 + qJDD(3) * mrSges(4,1) - t292 * mrSges(4,3) + qJD(3) * t298 - t290 * t339 + t174;
t160 = t313 * t166 + t315 * t170;
t222 = Ifges(7,5) * t289 + Ifges(7,6) * t288 + Ifges(7,3) * t302;
t229 = Ifges(6,1) * t289 - Ifges(6,4) * t302 + Ifges(6,5) * t288;
t230 = Ifges(5,1) * t289 - Ifges(5,4) * t288 - Ifges(5,5) * t302;
t179 = -t285 * mrSges(7,2) - t302 * t264 + t336;
t228 = Ifges(7,1) * t289 + Ifges(7,4) * t288 + Ifges(7,5) * t302;
t328 = mrSges(7,1) * t189 - mrSges(7,3) * t186 - Ifges(7,4) * t252 - Ifges(7,2) * t251 - Ifges(7,6) * t285 - t302 * t228;
t322 = mrSges(6,1) * t195 - mrSges(6,2) * t192 + pkin(5) * t180 + qJ(6) * t179 - t328;
t226 = Ifges(6,4) * t289 - Ifges(6,2) * t302 + Ifges(6,6) * t288;
t342 = -Ifges(5,5) * t289 + Ifges(5,6) * t288 + Ifges(5,3) * t302 - t226;
t155 = (-t230 - t229) * t302 + (t222 + t342) * t289 + (-Ifges(5,6) + Ifges(6,6)) * t285 + (Ifges(5,4) - Ifges(6,5)) * t252 + (-Ifges(5,2) - Ifges(6,3)) * t251 + mrSges(5,3) * t197 + mrSges(5,1) * t326 - pkin(4) * t177 - t322;
t225 = Ifges(7,4) * t289 + Ifges(7,2) * t288 + Ifges(7,6) * t302;
t227 = Ifges(5,4) * t289 - Ifges(5,2) * t288 - Ifges(5,6) * t302;
t178 = t285 * mrSges(7,1) + t302 * t261 - t332;
t223 = Ifges(6,5) * t289 - Ifges(6,6) * t302 + Ifges(6,3) * t288;
t327 = mrSges(7,2) * t189 - mrSges(7,3) * t182 + Ifges(7,1) * t252 + Ifges(7,4) * t251 + Ifges(7,5) * t285 + t288 * t222;
t321 = mrSges(6,2) * t194 - mrSges(6,3) * t195 + Ifges(6,1) * t252 - Ifges(6,4) * t285 + Ifges(6,5) * t251 - qJ(6) * t178 - t302 * t223 + t327;
t161 = (t227 - t225) * t302 + t342 * t288 - Ifges(5,5) * t285 - Ifges(5,4) * t251 + Ifges(5,1) * t252 - mrSges(5,3) * t196 - mrSges(5,2) * t326 - qJ(5) * t177 + t321;
t272 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t313 + Ifges(4,2) * t315) * qJD(2);
t273 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t313 + Ifges(4,4) * t315) * qJD(2);
t354 = mrSges(4,1) * t204 - mrSges(4,2) * t205 + Ifges(4,5) * t292 + Ifges(4,6) * t293 + Ifges(4,3) * qJDD(3) + pkin(3) * t174 + pkin(9) * t168 + (t272 * t313 - t273 * t315) * qJD(2) + t155 * t352 + t312 * t161;
t145 = -mrSges(3,1) * t268 + mrSges(3,3) * t221 + t318 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t160 - t354;
t334 = t315 * t166 - t313 * t170;
t158 = m(3) * t221 - t318 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t334;
t167 = t172 * t352 + t312 * t173;
t324 = -m(4) * t214 + t293 * mrSges(4,1) - t292 * mrSges(4,2) - t297 * t339 + t298 * t338 - t167;
t163 = m(3) * t220 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,2) + t324;
t154 = t316 * t158 - t314 * t163;
t359 = pkin(7) * t154 + t145 * t316;
t329 = mrSges(7,1) * t182 - mrSges(7,2) * t186 + Ifges(7,5) * t252 + Ifges(7,6) * t251 + Ifges(7,3) * t285 + t289 * t225 - t288 * t228;
t323 = mrSges(6,1) * t194 - mrSges(6,3) * t192 - Ifges(6,4) * t252 + Ifges(6,2) * t285 - Ifges(6,6) * t251 + pkin(5) * t178 - t288 * t229 + t329;
t355 = (t227 - t223) * t289 + mrSges(5,1) * t196 - mrSges(5,2) * t197 + Ifges(5,5) * t252 - Ifges(5,6) * t251 - Ifges(5,3) * t285 + pkin(4) * (-t252 * mrSges(6,2) - t289 * t256 - t178 + t357) + qJ(5) * (-t251 * mrSges(6,2) - t288 * t256 + t179 + t358) + t288 * t230 - t323;
t347 = t163 * t316;
t345 = t302 * t225;
t159 = m(3) * t268 + t160;
t149 = t158 * t343 - t309 * t159 + t311 * t347;
t271 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t313 + Ifges(4,6) * t315) * qJD(2);
t150 = mrSges(4,2) * t214 - mrSges(4,3) * t204 + Ifges(4,1) * t292 + Ifges(4,4) * t293 + Ifges(4,5) * qJDD(3) - pkin(9) * t167 - qJD(3) * t272 - t312 * t155 + t161 * t352 + t271 * t338;
t151 = -mrSges(4,1) * t214 + mrSges(4,3) * t205 + Ifges(4,4) * t292 + Ifges(4,2) * t293 + Ifges(4,6) * qJDD(3) - pkin(3) * t167 + qJD(3) * t273 - t271 * t339 - t355;
t141 = mrSges(3,1) * t220 - mrSges(3,2) * t221 + Ifges(3,3) * qJDD(2) + pkin(2) * t324 + pkin(8) * t334 + t313 * t150 + t315 * t151;
t143 = mrSges(3,2) * t268 - mrSges(3,3) * t220 + Ifges(3,5) * qJDD(2) - t318 * Ifges(3,6) - pkin(8) * t160 + t315 * t150 - t313 * t151;
t325 = mrSges(2,1) * t295 - mrSges(2,2) * t296 + pkin(1) * t149 + t311 * t141 + t143 * t344 + t309 * t359;
t152 = m(2) * t296 + t154;
t148 = t311 * t159 + (t158 * t314 + t347) * t309;
t146 = m(2) * t295 + t149;
t139 = mrSges(2,2) * t307 - mrSges(2,3) * t295 + t316 * t143 - t314 * t145 + (-t148 * t309 - t149 * t311) * pkin(7);
t138 = -mrSges(2,1) * t307 + mrSges(2,3) * t296 - pkin(1) * t148 - t309 * t141 + (t143 * t314 + t359) * t311;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t310 * t139 - t308 * t138 - qJ(1) * (t310 * t146 + t308 * t152), t139, t143, t150, t161, -t288 * t226 + t321 - t345, t327 - t345; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t308 * t139 + t310 * t138 + qJ(1) * (-t308 * t146 + t310 * t152), t138, t145, t151, t155, -t289 * t223 - t323, -t289 * t222 - t328; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t325, t325, t141, t354, t355, t302 * t229 - Ifges(6,6) * t285 + Ifges(6,3) * t251 + Ifges(6,5) * t252 + t322 + (-t222 + t226) * t289, t329;];
m_new  = t1;
