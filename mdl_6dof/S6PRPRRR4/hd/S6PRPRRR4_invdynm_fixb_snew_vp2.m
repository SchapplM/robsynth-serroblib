% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 01:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:57:24
% EndTime: 2019-05-05 00:57:59
% DurationCPUTime: 28.18s
% Computational Cost: add. (500839->321), mult. (1116051->411), div. (0->0), fcn. (854620->14), ass. (0->149)
t303 = sin(pkin(11));
t306 = cos(pkin(11));
t288 = t303 * g(1) - t306 * g(2);
t289 = -t306 * g(1) - t303 * g(2);
t301 = -g(3) + qJDD(1);
t315 = cos(qJ(2));
t307 = cos(pkin(6));
t311 = sin(qJ(2));
t346 = t307 * t311;
t304 = sin(pkin(6));
t347 = t304 * t311;
t256 = t288 * t346 + t315 * t289 + t301 * t347;
t317 = qJD(2) ^ 2;
t250 = -t317 * pkin(2) + qJDD(2) * qJ(3) + t256;
t302 = sin(pkin(12));
t274 = -t304 * t288 + t307 * t301;
t305 = cos(pkin(12));
t341 = qJD(2) * qJD(3);
t345 = t305 * t274 - 0.2e1 * t302 * t341;
t353 = pkin(3) * t305;
t226 = (-pkin(8) * qJDD(2) + t317 * t353 - t250) * t302 + t345;
t229 = t302 * t274 + (t250 + 0.2e1 * t341) * t305;
t339 = qJDD(2) * t305;
t299 = t305 ^ 2;
t348 = t299 * t317;
t227 = -pkin(3) * t348 + pkin(8) * t339 + t229;
t310 = sin(qJ(4));
t314 = cos(qJ(4));
t213 = t310 * t226 + t314 * t227;
t343 = qJD(2) * t305;
t344 = qJD(2) * t302;
t278 = -t310 * t344 + t314 * t343;
t329 = t302 * t314 + t305 * t310;
t279 = t329 * qJD(2);
t259 = -t278 * mrSges(5,1) + t279 * mrSges(5,2);
t276 = t279 * qJD(4);
t340 = qJDD(2) * t302;
t265 = -t310 * t340 + t314 * t339 - t276;
t273 = qJD(4) * mrSges(5,1) - t279 * mrSges(5,3);
t264 = -t278 * pkin(4) - t279 * pkin(9);
t316 = qJD(4) ^ 2;
t205 = -t316 * pkin(4) + qJDD(4) * pkin(9) + t278 * t264 + t213;
t298 = t302 ^ 2;
t255 = -t311 * t289 + (t288 * t307 + t301 * t304) * t315;
t325 = qJDD(3) - t255;
t239 = (-pkin(2) - t353) * qJDD(2) + (-qJ(3) + (-t298 - t299) * pkin(8)) * t317 + t325;
t342 = t278 * qJD(4);
t266 = qJDD(2) * t329 + t342;
t217 = (-t266 - t342) * pkin(9) + (-t265 + t276) * pkin(4) + t239;
t309 = sin(qJ(5));
t313 = cos(qJ(5));
t200 = -t309 * t205 + t313 * t217;
t268 = t313 * qJD(4) - t309 * t279;
t238 = t268 * qJD(5) + t309 * qJDD(4) + t313 * t266;
t263 = qJDD(5) - t265;
t269 = t309 * qJD(4) + t313 * t279;
t277 = qJD(5) - t278;
t198 = (t268 * t277 - t238) * pkin(10) + (t268 * t269 + t263) * pkin(5) + t200;
t201 = t313 * t205 + t309 * t217;
t237 = -t269 * qJD(5) + t313 * qJDD(4) - t309 * t266;
t249 = t277 * pkin(5) - t269 * pkin(10);
t267 = t268 ^ 2;
t199 = -t267 * pkin(5) + t237 * pkin(10) - t277 * t249 + t201;
t308 = sin(qJ(6));
t312 = cos(qJ(6));
t196 = t312 * t198 - t308 * t199;
t240 = t312 * t268 - t308 * t269;
t210 = t240 * qJD(6) + t308 * t237 + t312 * t238;
t241 = t308 * t268 + t312 * t269;
t222 = -t240 * mrSges(7,1) + t241 * mrSges(7,2);
t275 = qJD(6) + t277;
t230 = -t275 * mrSges(7,2) + t240 * mrSges(7,3);
t258 = qJDD(6) + t263;
t191 = m(7) * t196 + t258 * mrSges(7,1) - t210 * mrSges(7,3) - t241 * t222 + t275 * t230;
t197 = t308 * t198 + t312 * t199;
t209 = -t241 * qJD(6) + t312 * t237 - t308 * t238;
t231 = t275 * mrSges(7,1) - t241 * mrSges(7,3);
t192 = m(7) * t197 - t258 * mrSges(7,2) + t209 * mrSges(7,3) + t240 * t222 - t275 * t231;
t183 = t312 * t191 + t308 * t192;
t243 = -t268 * mrSges(6,1) + t269 * mrSges(6,2);
t247 = -t277 * mrSges(6,2) + t268 * mrSges(6,3);
t181 = m(6) * t200 + t263 * mrSges(6,1) - t238 * mrSges(6,3) - t269 * t243 + t277 * t247 + t183;
t248 = t277 * mrSges(6,1) - t269 * mrSges(6,3);
t335 = -t308 * t191 + t312 * t192;
t182 = m(6) * t201 - t263 * mrSges(6,2) + t237 * mrSges(6,3) + t268 * t243 - t277 * t248 + t335;
t336 = -t309 * t181 + t313 * t182;
t174 = m(5) * t213 - qJDD(4) * mrSges(5,2) + t265 * mrSges(5,3) - qJD(4) * t273 + t278 * t259 + t336;
t212 = t314 * t226 - t310 * t227;
t272 = -qJD(4) * mrSges(5,2) + t278 * mrSges(5,3);
t204 = -qJDD(4) * pkin(4) - t316 * pkin(9) + t279 * t264 - t212;
t202 = -t237 * pkin(5) - t267 * pkin(10) + t269 * t249 + t204;
t326 = m(7) * t202 - t209 * mrSges(7,1) + t210 * mrSges(7,2) - t240 * t230 + t241 * t231;
t320 = -m(6) * t204 + t237 * mrSges(6,1) - t238 * mrSges(6,2) + t268 * t247 - t269 * t248 - t326;
t187 = m(5) * t212 + qJDD(4) * mrSges(5,1) - t266 * mrSges(5,3) + qJD(4) * t272 - t279 * t259 + t320;
t165 = t310 * t174 + t314 * t187;
t228 = -t302 * t250 + t345;
t351 = mrSges(4,2) * t302;
t328 = mrSges(4,3) * qJDD(2) + t317 * (-mrSges(4,1) * t305 + t351);
t163 = m(4) * t228 - t302 * t328 + t165;
t337 = t314 * t174 - t310 * t187;
t164 = m(4) * t229 + t305 * t328 + t337;
t158 = t305 * t163 + t302 * t164;
t332 = Ifges(4,5) * t302 + Ifges(4,6) * t305;
t218 = Ifges(7,5) * t241 + Ifges(7,6) * t240 + Ifges(7,3) * t275;
t220 = Ifges(7,1) * t241 + Ifges(7,4) * t240 + Ifges(7,5) * t275;
t184 = -mrSges(7,1) * t202 + mrSges(7,3) * t197 + Ifges(7,4) * t210 + Ifges(7,2) * t209 + Ifges(7,6) * t258 - t241 * t218 + t275 * t220;
t219 = Ifges(7,4) * t241 + Ifges(7,2) * t240 + Ifges(7,6) * t275;
t185 = mrSges(7,2) * t202 - mrSges(7,3) * t196 + Ifges(7,1) * t210 + Ifges(7,4) * t209 + Ifges(7,5) * t258 + t240 * t218 - t275 * t219;
t232 = Ifges(6,5) * t269 + Ifges(6,6) * t268 + Ifges(6,3) * t277;
t234 = Ifges(6,1) * t269 + Ifges(6,4) * t268 + Ifges(6,5) * t277;
t167 = -mrSges(6,1) * t204 + mrSges(6,3) * t201 + Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t263 - pkin(5) * t326 + pkin(10) * t335 + t312 * t184 + t308 * t185 - t269 * t232 + t277 * t234;
t233 = Ifges(6,4) * t269 + Ifges(6,2) * t268 + Ifges(6,6) * t277;
t169 = mrSges(6,2) * t204 - mrSges(6,3) * t200 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t263 - pkin(10) * t183 - t308 * t184 + t312 * t185 + t268 * t232 - t277 * t233;
t252 = Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * qJD(4);
t253 = Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * qJD(4);
t322 = -mrSges(5,1) * t212 + mrSges(5,2) * t213 - Ifges(5,5) * t266 - Ifges(5,6) * t265 - Ifges(5,3) * qJDD(4) - pkin(4) * t320 - pkin(9) * t336 - t313 * t167 - t309 * t169 - t279 * t252 + t278 * t253;
t333 = Ifges(4,4) * t302 + Ifges(4,2) * t305;
t334 = Ifges(4,1) * t302 + Ifges(4,4) * t305;
t354 = -mrSges(4,1) * t228 + mrSges(4,2) * t229 - pkin(3) * t165 - (t333 * t344 - t334 * t343) * qJD(2) + t322;
t144 = (Ifges(3,6) - t332) * qJDD(2) + t317 * Ifges(3,5) - mrSges(3,1) * t274 + mrSges(3,3) * t256 - pkin(2) * t158 + t354;
t338 = -t302 * t163 + t305 * t164;
t156 = m(3) * t256 - t317 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t338;
t246 = -qJDD(2) * pkin(2) - t317 * qJ(3) + t325;
t176 = t313 * t181 + t309 * t182;
t323 = m(5) * t239 - t265 * mrSges(5,1) + t266 * mrSges(5,2) - t278 * t272 + t279 * t273 + t176;
t321 = -m(4) * t246 + mrSges(4,1) * t339 - t323 + (t298 * t317 + t348) * mrSges(4,3);
t171 = -t317 * mrSges(3,2) + m(3) * t255 + t321 + (mrSges(3,1) - t351) * qJDD(2);
t152 = t315 * t156 - t311 * t171;
t355 = pkin(7) * t152 + t144 * t315;
t349 = t171 * t315;
t157 = m(3) * t274 + t158;
t149 = t156 * t346 - t304 * t157 + t307 * t349;
t251 = Ifges(5,5) * t279 + Ifges(5,6) * t278 + Ifges(5,3) * qJD(4);
t153 = mrSges(5,2) * t239 - mrSges(5,3) * t212 + Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * qJDD(4) - pkin(9) * t176 - qJD(4) * t252 - t309 * t167 + t313 * t169 + t278 * t251;
t324 = -mrSges(7,1) * t196 + mrSges(7,2) * t197 - Ifges(7,5) * t210 - Ifges(7,6) * t209 - Ifges(7,3) * t258 - t241 * t219 + t240 * t220;
t318 = mrSges(6,1) * t200 - mrSges(6,2) * t201 + Ifges(6,5) * t238 + Ifges(6,6) * t237 + Ifges(6,3) * t263 + pkin(5) * t183 + t269 * t233 - t268 * t234 - t324;
t159 = -mrSges(5,1) * t239 + mrSges(5,3) * t213 + Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * qJDD(4) - pkin(4) * t176 + qJD(4) * t253 - t279 * t251 - t318;
t284 = t332 * qJD(2);
t142 = -mrSges(4,1) * t246 + mrSges(4,3) * t229 - pkin(3) * t323 + pkin(8) * t337 + qJDD(2) * t333 + t310 * t153 + t314 * t159 - t284 * t344;
t145 = mrSges(4,2) * t246 - mrSges(4,3) * t228 - pkin(8) * t165 + qJDD(2) * t334 + t314 * t153 - t310 * t159 + t284 * t343;
t139 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t255 - mrSges(3,2) * t256 + t302 * t145 + t305 * t142 + pkin(2) * (-mrSges(4,2) * t340 + t321) + qJ(3) * t338;
t141 = mrSges(3,2) * t274 - mrSges(3,3) * t255 + Ifges(3,5) * qJDD(2) - t317 * Ifges(3,6) - qJ(3) * t158 - t302 * t142 + t305 * t145;
t327 = mrSges(2,1) * t288 - mrSges(2,2) * t289 + pkin(1) * t149 + t307 * t139 + t141 * t347 + t304 * t355;
t150 = m(2) * t289 + t152;
t148 = t307 * t157 + (t156 * t311 + t349) * t304;
t146 = m(2) * t288 + t149;
t137 = mrSges(2,2) * t301 - mrSges(2,3) * t288 + t315 * t141 - t311 * t144 + (-t148 * t304 - t149 * t307) * pkin(7);
t136 = -mrSges(2,1) * t301 + mrSges(2,3) * t289 - pkin(1) * t148 - t304 * t139 + (t141 * t311 + t355) * t307;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t306 * t137 - t303 * t136 - qJ(1) * (t306 * t146 + t303 * t150), t137, t141, t145, t153, t169, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t303 * t137 + t306 * t136 + qJ(1) * (-t303 * t146 + t306 * t150), t136, t144, t142, t159, t167, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t327, t327, t139, qJDD(2) * t332 - t354, -t322, t318, -t324;];
m_new  = t1;
