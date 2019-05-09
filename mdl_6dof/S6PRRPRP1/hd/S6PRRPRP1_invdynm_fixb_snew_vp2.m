% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:36:32
% EndTime: 2019-05-05 03:36:55
% DurationCPUTime: 13.83s
% Computational Cost: add. (242094->335), mult. (518892->422), div. (0->0), fcn. (365456->12), ass. (0->135)
t356 = -2 * qJD(4);
t306 = sin(pkin(10));
t309 = cos(pkin(10));
t296 = t306 * g(1) - t309 * g(2);
t297 = -t309 * g(1) - t306 * g(2);
t304 = -g(3) + qJDD(1);
t316 = cos(qJ(2));
t310 = cos(pkin(6));
t313 = sin(qJ(2));
t347 = t310 * t313;
t307 = sin(pkin(6));
t348 = t307 * t313;
t256 = t296 * t347 + t316 * t297 + t304 * t348;
t318 = qJD(2) ^ 2;
t251 = -t318 * pkin(2) + qJDD(2) * pkin(8) + t256;
t277 = -t307 * t296 + t310 * t304;
t312 = sin(qJ(3));
t315 = cos(qJ(3));
t218 = -t312 * t251 + t315 * t277;
t341 = qJD(2) * qJD(3);
t338 = t315 * t341;
t293 = t312 * qJDD(2) + t338;
t211 = (-t293 + t338) * qJ(4) + (t312 * t315 * t318 + qJDD(3)) * pkin(3) + t218;
t219 = t315 * t251 + t312 * t277;
t294 = t315 * qJDD(2) - t312 * t341;
t344 = qJD(2) * t312;
t298 = qJD(3) * pkin(3) - qJ(4) * t344;
t303 = t315 ^ 2;
t212 = -t303 * t318 * pkin(3) + t294 * qJ(4) - qJD(3) * t298 + t219;
t305 = sin(pkin(11));
t308 = cos(pkin(11));
t281 = (t305 * t312 - t308 * t315) * qJD(2);
t204 = t305 * t211 + t308 * t212 + t281 * t356;
t282 = (t305 * t315 + t308 * t312) * qJD(2);
t258 = t281 * mrSges(5,1) + t282 * mrSges(5,2);
t268 = -t305 * t293 + t308 * t294;
t276 = qJD(3) * mrSges(5,1) - t282 * mrSges(5,3);
t259 = t281 * pkin(4) - t282 * pkin(9);
t317 = qJD(3) ^ 2;
t201 = -t317 * pkin(4) + qJDD(3) * pkin(9) - t281 * t259 + t204;
t255 = -t313 * t297 + (t296 * t310 + t304 * t307) * t316;
t325 = -qJDD(2) * pkin(2) - t255;
t216 = -t294 * pkin(3) + qJDD(4) + t298 * t344 + (-qJ(4) * t303 - pkin(8)) * t318 + t325;
t269 = t308 * t293 + t305 * t294;
t207 = (qJD(3) * t281 - t269) * pkin(9) + (qJD(3) * t282 - t268) * pkin(4) + t216;
t311 = sin(qJ(5));
t314 = cos(qJ(5));
t195 = -t311 * t201 + t314 * t207;
t271 = t314 * qJD(3) - t311 * t282;
t238 = t271 * qJD(5) + t311 * qJDD(3) + t314 * t269;
t272 = t311 * qJD(3) + t314 * t282;
t241 = -t271 * mrSges(7,1) + t272 * mrSges(7,2);
t242 = -t271 * mrSges(6,1) + t272 * mrSges(6,2);
t280 = qJD(5) + t281;
t246 = -t280 * mrSges(6,2) + t271 * mrSges(6,3);
t267 = qJDD(5) - t268;
t191 = -0.2e1 * qJD(6) * t272 + (t271 * t280 - t238) * qJ(6) + (t271 * t272 + t267) * pkin(5) + t195;
t245 = -t280 * mrSges(7,2) + t271 * mrSges(7,3);
t340 = m(7) * t191 + t267 * mrSges(7,1) + t280 * t245;
t180 = m(6) * t195 + t267 * mrSges(6,1) + t280 * t246 + (-t241 - t242) * t272 + (-mrSges(6,3) - mrSges(7,3)) * t238 + t340;
t196 = t314 * t201 + t311 * t207;
t237 = -t272 * qJD(5) + t314 * qJDD(3) - t311 * t269;
t247 = t280 * pkin(5) - t272 * qJ(6);
t270 = t271 ^ 2;
t194 = -t270 * pkin(5) + t237 * qJ(6) + 0.2e1 * qJD(6) * t271 - t280 * t247 + t196;
t339 = m(7) * t194 + t237 * mrSges(7,3) + t271 * t241;
t248 = t280 * mrSges(7,1) - t272 * mrSges(7,3);
t345 = -t280 * mrSges(6,1) + t272 * mrSges(6,3) - t248;
t351 = -mrSges(6,2) - mrSges(7,2);
t185 = m(6) * t196 + t237 * mrSges(6,3) + t271 * t242 + t267 * t351 + t280 * t345 + t339;
t335 = -t311 * t180 + t314 * t185;
t173 = m(5) * t204 - qJDD(3) * mrSges(5,2) + t268 * mrSges(5,3) - qJD(3) * t276 - t281 * t258 + t335;
t203 = t308 * t211 - t305 * t212 + t282 * t356;
t275 = -qJD(3) * mrSges(5,2) - t281 * mrSges(5,3);
t200 = -qJDD(3) * pkin(4) - t317 * pkin(9) + t282 * t259 - t203;
t198 = -t237 * pkin(5) - t270 * qJ(6) + t272 * t247 + qJDD(6) + t200;
t333 = -m(7) * t198 + t237 * mrSges(7,1) + t271 * t245;
t322 = -m(6) * t200 + t237 * mrSges(6,1) + t351 * t238 + t271 * t246 + t345 * t272 + t333;
t182 = m(5) * t203 + qJDD(3) * mrSges(5,1) - t269 * mrSges(5,3) + qJD(3) * t275 - t282 * t258 + t322;
t166 = t305 * t173 + t308 * t182;
t292 = (-mrSges(4,1) * t315 + mrSges(4,2) * t312) * qJD(2);
t343 = qJD(2) * t315;
t300 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t343;
t164 = m(4) * t218 + qJDD(3) * mrSges(4,1) - t293 * mrSges(4,3) + qJD(3) * t300 - t292 * t344 + t166;
t299 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t344;
t336 = t308 * t173 - t305 * t182;
t165 = m(4) * t219 - qJDD(3) * mrSges(4,2) + t294 * mrSges(4,3) - qJD(3) * t299 + t292 * t343 + t336;
t159 = t315 * t164 + t312 * t165;
t285 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t312 + Ifges(4,2) * t315) * qJD(2);
t286 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t312 + Ifges(4,4) * t315) * qJD(2);
t220 = Ifges(7,5) * t272 + Ifges(7,6) * t271 + Ifges(7,3) * t280;
t221 = Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t280;
t225 = Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t280;
t224 = Ifges(7,1) * t272 + Ifges(7,4) * t271 + Ifges(7,5) * t280;
t329 = -mrSges(7,1) * t198 + mrSges(7,3) * t194 + Ifges(7,4) * t238 + Ifges(7,2) * t237 + Ifges(7,6) * t267 + t280 * t224;
t168 = Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t267 + t280 * t225 - mrSges(6,1) * t200 + mrSges(6,3) * t196 - pkin(5) * (t238 * mrSges(7,2) - t333) + qJ(6) * (-t267 * mrSges(7,2) - t280 * t248 + t339) + (-pkin(5) * t248 - t220 - t221) * t272 + t329;
t188 = -t238 * mrSges(7,3) - t272 * t241 + t340;
t222 = Ifges(7,4) * t272 + Ifges(7,2) * t271 + Ifges(7,6) * t280;
t223 = Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t280;
t327 = mrSges(7,2) * t198 - mrSges(7,3) * t191 + Ifges(7,1) * t238 + Ifges(7,4) * t237 + Ifges(7,5) * t267 + t271 * t220;
t175 = mrSges(6,2) * t200 - mrSges(6,3) * t195 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t267 - qJ(6) * t188 + t271 * t221 + (-t222 - t223) * t280 + t327;
t253 = Ifges(5,4) * t282 - Ifges(5,2) * t281 + Ifges(5,6) * qJD(3);
t254 = Ifges(5,1) * t282 - Ifges(5,4) * t281 + Ifges(5,5) * qJD(3);
t323 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t269 - Ifges(5,6) * t268 - Ifges(5,3) * qJDD(3) - pkin(4) * t322 - pkin(9) * t335 - t314 * t168 - t311 * t175 - t282 * t253 - t281 * t254;
t353 = mrSges(4,1) * t218 - mrSges(4,2) * t219 + Ifges(4,5) * t293 + Ifges(4,6) * t294 + Ifges(4,3) * qJDD(3) + pkin(3) * t166 + (t312 * t285 - t315 * t286) * qJD(2) - t323;
t144 = -mrSges(3,1) * t277 + mrSges(3,3) * t256 + t318 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t159 - t353;
t337 = -t312 * t164 + t315 * t165;
t157 = m(3) * t256 - t318 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t337;
t250 = -t318 * pkin(8) + t325;
t177 = t314 * t180 + t311 * t185;
t324 = m(5) * t216 - t268 * mrSges(5,1) + t269 * mrSges(5,2) + t281 * t275 + t282 * t276 + t177;
t321 = -m(4) * t250 + t294 * mrSges(4,1) - t293 * mrSges(4,2) - t299 * t344 + t300 * t343 - t324;
t170 = m(3) * t255 + qJDD(2) * mrSges(3,1) - t318 * mrSges(3,2) + t321;
t153 = t316 * t157 - t313 * t170;
t355 = pkin(7) * t153 + t144 * t316;
t328 = -mrSges(7,1) * t191 + mrSges(7,2) * t194 - Ifges(7,5) * t238 - Ifges(7,6) * t237 - Ifges(7,3) * t267 - t272 * t222;
t354 = mrSges(6,1) * t195 - mrSges(6,2) * t196 + Ifges(6,5) * t238 + Ifges(6,6) * t237 + Ifges(6,3) * t267 + pkin(5) * t188 + t272 * t223 - (t225 + t224) * t271 - t328;
t349 = t170 * t316;
t158 = m(3) * t277 + t159;
t149 = t157 * t347 - t307 * t158 + t310 * t349;
t252 = Ifges(5,5) * t282 - Ifges(5,6) * t281 + Ifges(5,3) * qJD(3);
t154 = mrSges(5,2) * t216 - mrSges(5,3) * t203 + Ifges(5,1) * t269 + Ifges(5,4) * t268 + Ifges(5,5) * qJDD(3) - pkin(9) * t177 - qJD(3) * t253 - t311 * t168 + t314 * t175 - t281 * t252;
t160 = -mrSges(5,1) * t216 + mrSges(5,3) * t204 + Ifges(5,4) * t269 + Ifges(5,2) * t268 + Ifges(5,6) * qJDD(3) - pkin(4) * t177 + qJD(3) * t254 - t282 * t252 - t354;
t284 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t312 + Ifges(4,6) * t315) * qJD(2);
t145 = -mrSges(4,1) * t250 + mrSges(4,3) * t219 + Ifges(4,4) * t293 + Ifges(4,2) * t294 + Ifges(4,6) * qJDD(3) - pkin(3) * t324 + qJ(4) * t336 + qJD(3) * t286 + t305 * t154 + t308 * t160 - t284 * t344;
t150 = mrSges(4,2) * t250 - mrSges(4,3) * t218 + Ifges(4,1) * t293 + Ifges(4,4) * t294 + Ifges(4,5) * qJDD(3) - qJ(4) * t166 - qJD(3) * t285 + t308 * t154 - t305 * t160 + t284 * t343;
t140 = mrSges(3,1) * t255 - mrSges(3,2) * t256 + Ifges(3,3) * qJDD(2) + pkin(2) * t321 + pkin(8) * t337 + t315 * t145 + t312 * t150;
t142 = mrSges(3,2) * t277 - mrSges(3,3) * t255 + Ifges(3,5) * qJDD(2) - t318 * Ifges(3,6) - pkin(8) * t159 - t312 * t145 + t315 * t150;
t326 = mrSges(2,1) * t296 - mrSges(2,2) * t297 + pkin(1) * t149 + t310 * t140 + t142 * t348 + t307 * t355;
t151 = m(2) * t297 + t153;
t148 = t310 * t158 + (t157 * t313 + t349) * t307;
t146 = m(2) * t296 + t149;
t138 = mrSges(2,2) * t304 - mrSges(2,3) * t296 + t316 * t142 - t313 * t144 + (-t148 * t307 - t149 * t310) * pkin(7);
t137 = -mrSges(2,1) * t304 + mrSges(2,3) * t297 - pkin(1) * t148 - t307 * t140 + (t142 * t313 + t355) * t310;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t309 * t138 - t306 * t137 - qJ(1) * (t309 * t146 + t306 * t151), t138, t142, t150, t154, t175, -t280 * t222 + t327; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t306 * t138 + t309 * t137 + qJ(1) * (-t306 * t146 + t309 * t151), t137, t144, t145, t160, t168, -t272 * t220 + t329; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t326, t326, t140, t353, -t323, t354, -t271 * t224 - t328;];
m_new  = t1;
