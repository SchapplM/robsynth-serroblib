% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:41:44
% EndTime: 2019-05-04 23:42:05
% DurationCPUTime: 12.69s
% Computational Cost: add. (223150->315), mult. (492138->390), div. (0->0), fcn. (365152->12), ass. (0->137)
t310 = qJD(2) ^ 2;
t298 = sin(pkin(10));
t301 = cos(pkin(10));
t285 = t298 * g(1) - t301 * g(2);
t286 = -t301 * g(1) - t298 * g(2);
t296 = -g(3) + qJDD(1);
t308 = cos(qJ(2));
t302 = cos(pkin(6));
t305 = sin(qJ(2));
t344 = t302 * t305;
t299 = sin(pkin(6));
t345 = t299 * t305;
t252 = t285 * t344 + t308 * t286 + t296 * t345;
t247 = -t310 * pkin(2) + qJDD(2) * qJ(3) + t252;
t297 = sin(pkin(11));
t273 = -t299 * t285 + t302 * t296;
t300 = cos(pkin(11));
t337 = qJD(2) * qJD(3);
t341 = t300 * t273 - 0.2e1 * t297 * t337;
t352 = pkin(3) * t300;
t209 = (-pkin(8) * qJDD(2) + t310 * t352 - t247) * t297 + t341;
t213 = t297 * t273 + (t247 + 0.2e1 * t337) * t300;
t336 = qJDD(2) * t300;
t294 = t300 ^ 2;
t346 = t294 * t310;
t210 = -pkin(3) * t346 + pkin(8) * t336 + t213;
t304 = sin(qJ(4));
t307 = cos(qJ(4));
t200 = t304 * t209 + t307 * t210;
t323 = t297 * t304 - t300 * t307;
t275 = t323 * qJD(2);
t324 = t297 * t307 + t300 * t304;
t276 = t324 * qJD(2);
t258 = t275 * mrSges(5,1) + t276 * mrSges(5,2);
t338 = t276 * qJD(4);
t264 = -qJDD(2) * t323 - t338;
t272 = qJD(4) * mrSges(5,1) - t276 * mrSges(5,3);
t263 = t275 * pkin(4) - t276 * pkin(9);
t309 = qJD(4) ^ 2;
t197 = -t309 * pkin(4) + qJDD(4) * pkin(9) - t275 * t263 + t200;
t293 = t297 ^ 2;
t251 = -t305 * t286 + (t285 * t302 + t296 * t299) * t308;
t317 = qJDD(3) - t251;
t234 = (-pkin(2) - t352) * qJDD(2) + (-qJ(3) + (-t293 - t294) * pkin(8)) * t310 + t317;
t339 = t275 * qJD(4);
t265 = qJDD(2) * t324 - t339;
t203 = (-t265 + t339) * pkin(9) + (-t264 + t338) * pkin(4) + t234;
t303 = sin(qJ(5));
t306 = cos(qJ(5));
t191 = -t303 * t197 + t306 * t203;
t267 = t306 * qJD(4) - t303 * t276;
t233 = t267 * qJD(5) + t303 * qJDD(4) + t306 * t265;
t268 = t303 * qJD(4) + t306 * t276;
t237 = -t267 * mrSges(7,1) + t268 * mrSges(7,2);
t238 = -t267 * mrSges(6,1) + t268 * mrSges(6,2);
t274 = qJD(5) + t275;
t243 = -t274 * mrSges(6,2) + t267 * mrSges(6,3);
t262 = qJDD(5) - t264;
t187 = -0.2e1 * qJD(6) * t268 + (t267 * t274 - t233) * qJ(6) + (t267 * t268 + t262) * pkin(5) + t191;
t242 = -t274 * mrSges(7,2) + t267 * mrSges(7,3);
t335 = m(7) * t187 + t262 * mrSges(7,1) + t274 * t242;
t176 = m(6) * t191 + t262 * mrSges(6,1) + t274 * t243 + (-t237 - t238) * t268 + (-mrSges(6,3) - mrSges(7,3)) * t233 + t335;
t192 = t306 * t197 + t303 * t203;
t232 = -t268 * qJD(5) + t306 * qJDD(4) - t303 * t265;
t244 = t274 * pkin(5) - t268 * qJ(6);
t266 = t267 ^ 2;
t190 = -t266 * pkin(5) + t232 * qJ(6) + 0.2e1 * qJD(6) * t267 - t274 * t244 + t192;
t334 = m(7) * t190 + t232 * mrSges(7,3) + t267 * t237;
t245 = t274 * mrSges(7,1) - t268 * mrSges(7,3);
t342 = -t274 * mrSges(6,1) + t268 * mrSges(6,3) - t245;
t350 = -mrSges(6,2) - mrSges(7,2);
t180 = m(6) * t192 + t232 * mrSges(6,3) + t267 * t238 + t262 * t350 + t274 * t342 + t334;
t331 = -t303 * t176 + t306 * t180;
t169 = m(5) * t200 - qJDD(4) * mrSges(5,2) + t264 * mrSges(5,3) - qJD(4) * t272 - t275 * t258 + t331;
t199 = t307 * t209 - t304 * t210;
t271 = -qJD(4) * mrSges(5,2) - t275 * mrSges(5,3);
t196 = -qJDD(4) * pkin(4) - t309 * pkin(9) + t276 * t263 - t199;
t194 = -t232 * pkin(5) - t266 * qJ(6) + t268 * t244 + qJDD(6) + t196;
t330 = -m(7) * t194 + t232 * mrSges(7,1) + t267 * t242;
t314 = -m(6) * t196 + t232 * mrSges(6,1) + t233 * t350 + t267 * t243 + t268 * t342 + t330;
t181 = m(5) * t199 + qJDD(4) * mrSges(5,1) - t265 * mrSges(5,3) + qJD(4) * t271 - t276 * t258 + t314;
t162 = t304 * t169 + t307 * t181;
t212 = -t297 * t247 + t341;
t349 = mrSges(4,2) * t297;
t322 = mrSges(4,3) * qJDD(2) + t310 * (-mrSges(4,1) * t300 + t349);
t160 = m(4) * t212 - t297 * t322 + t162;
t332 = t307 * t169 - t304 * t181;
t161 = m(4) * t213 + t300 * t322 + t332;
t155 = t300 * t160 + t297 * t161;
t327 = Ifges(4,5) * t297 + Ifges(4,6) * t300;
t215 = Ifges(7,5) * t268 + Ifges(7,6) * t267 + Ifges(7,3) * t274;
t216 = Ifges(6,5) * t268 + Ifges(6,6) * t267 + Ifges(6,3) * t274;
t220 = Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * t274;
t219 = Ifges(7,1) * t268 + Ifges(7,4) * t267 + Ifges(7,5) * t274;
t321 = -mrSges(7,1) * t194 + mrSges(7,3) * t190 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t262 + t274 * t219;
t164 = Ifges(6,4) * t233 + Ifges(6,2) * t232 + Ifges(6,6) * t262 + t274 * t220 - mrSges(6,1) * t196 + mrSges(6,3) * t192 - pkin(5) * (t233 * mrSges(7,2) - t330) + qJ(6) * (-t262 * mrSges(7,2) - t274 * t245 + t334) + (-pkin(5) * t245 - t215 - t216) * t268 + t321;
t184 = -t233 * mrSges(7,3) - t268 * t237 + t335;
t217 = Ifges(7,4) * t268 + Ifges(7,2) * t267 + Ifges(7,6) * t274;
t218 = Ifges(6,4) * t268 + Ifges(6,2) * t267 + Ifges(6,6) * t274;
t319 = mrSges(7,2) * t194 - mrSges(7,3) * t187 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t262 + t267 * t215;
t171 = mrSges(6,2) * t196 - mrSges(6,3) * t191 + Ifges(6,1) * t233 + Ifges(6,4) * t232 + Ifges(6,5) * t262 - qJ(6) * t184 + t267 * t216 + (-t217 - t218) * t274 + t319;
t249 = Ifges(5,4) * t276 - Ifges(5,2) * t275 + Ifges(5,6) * qJD(4);
t250 = Ifges(5,1) * t276 - Ifges(5,4) * t275 + Ifges(5,5) * qJD(4);
t315 = -mrSges(5,1) * t199 + mrSges(5,2) * t200 - Ifges(5,5) * t265 - Ifges(5,6) * t264 - Ifges(5,3) * qJDD(4) - pkin(4) * t314 - pkin(9) * t331 - t306 * t164 - t303 * t171 - t276 * t249 - t275 * t250;
t328 = Ifges(4,4) * t297 + Ifges(4,2) * t300;
t329 = Ifges(4,1) * t297 + Ifges(4,4) * t300;
t353 = -mrSges(4,1) * t212 + mrSges(4,2) * t213 - pkin(3) * t162 - (t297 * t328 - t300 * t329) * t310 + t315;
t140 = (Ifges(3,6) - t327) * qJDD(2) + t310 * Ifges(3,5) - mrSges(3,1) * t273 + mrSges(3,3) * t252 - pkin(2) * t155 + t353;
t333 = -t297 * t160 + t300 * t161;
t153 = m(3) * t252 - t310 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t333;
t241 = -qJDD(2) * pkin(2) - t310 * qJ(3) + t317;
t173 = t306 * t176 + t303 * t180;
t316 = m(5) * t234 - t264 * mrSges(5,1) + t265 * mrSges(5,2) + t275 * t271 + t276 * t272 + t173;
t313 = -m(4) * t241 + mrSges(4,1) * t336 - t316 + (t293 * t310 + t346) * mrSges(4,3);
t166 = (mrSges(3,1) - t349) * qJDD(2) + t313 - t310 * mrSges(3,2) + m(3) * t251;
t149 = t308 * t153 - t305 * t166;
t355 = pkin(7) * t149 + t140 * t308;
t320 = -mrSges(7,1) * t187 + mrSges(7,2) * t190 - Ifges(7,5) * t233 - Ifges(7,6) * t232 - Ifges(7,3) * t262 - t268 * t217;
t354 = mrSges(6,1) * t191 - mrSges(6,2) * t192 + Ifges(6,5) * t233 + Ifges(6,6) * t232 + Ifges(6,3) * t262 + pkin(5) * t184 + t268 * t218 - (t220 + t219) * t267 - t320;
t347 = t166 * t308;
t340 = t310 * t327;
t154 = m(3) * t273 + t155;
t145 = t153 * t344 - t299 * t154 + t302 * t347;
t248 = Ifges(5,5) * t276 - Ifges(5,6) * t275 + Ifges(5,3) * qJD(4);
t150 = mrSges(5,2) * t234 - mrSges(5,3) * t199 + Ifges(5,1) * t265 + Ifges(5,4) * t264 + Ifges(5,5) * qJDD(4) - pkin(9) * t173 - qJD(4) * t249 - t303 * t164 + t306 * t171 - t275 * t248;
t156 = -mrSges(5,1) * t234 + mrSges(5,3) * t200 + Ifges(5,4) * t265 + Ifges(5,2) * t264 + Ifges(5,6) * qJDD(4) - pkin(4) * t173 + qJD(4) * t250 - t276 * t248 - t354;
t141 = -mrSges(4,1) * t241 + mrSges(4,3) * t213 - pkin(3) * t316 + pkin(8) * t332 + qJDD(2) * t328 + t304 * t150 + t307 * t156 - t297 * t340;
t146 = mrSges(4,2) * t241 - mrSges(4,3) * t212 - pkin(8) * t162 + qJDD(2) * t329 + t307 * t150 - t304 * t156 + t300 * t340;
t136 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t251 - mrSges(3,2) * t252 + t297 * t146 + t300 * t141 + pkin(2) * (-qJDD(2) * t349 + t313) + qJ(3) * t333;
t138 = mrSges(3,2) * t273 - mrSges(3,3) * t251 + Ifges(3,5) * qJDD(2) - t310 * Ifges(3,6) - qJ(3) * t155 - t297 * t141 + t300 * t146;
t318 = mrSges(2,1) * t285 - mrSges(2,2) * t286 + pkin(1) * t145 + t302 * t136 + t138 * t345 + t299 * t355;
t147 = m(2) * t286 + t149;
t144 = t302 * t154 + (t153 * t305 + t347) * t299;
t142 = m(2) * t285 + t145;
t134 = mrSges(2,2) * t296 - mrSges(2,3) * t285 + t308 * t138 - t305 * t140 + (-t144 * t299 - t145 * t302) * pkin(7);
t133 = -mrSges(2,1) * t296 + mrSges(2,3) * t286 - pkin(1) * t144 - t299 * t136 + (t138 * t305 + t355) * t302;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t301 * t134 - t298 * t133 - qJ(1) * (t301 * t142 + t298 * t147), t134, t138, t146, t150, t171, -t274 * t217 + t319; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t298 * t134 + t301 * t133 + qJ(1) * (-t298 * t142 + t301 * t147), t133, t140, t141, t156, t164, -t268 * t215 + t321; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t318, t318, t136, qJDD(2) * t327 - t353, -t315, t354, -t267 * t219 - t320;];
m_new  = t1;
