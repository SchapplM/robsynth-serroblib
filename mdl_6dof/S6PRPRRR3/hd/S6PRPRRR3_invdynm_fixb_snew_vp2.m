% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRR3
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
% Datum: 2019-05-05 00:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:40:07
% EndTime: 2019-05-05 00:40:43
% DurationCPUTime: 30.54s
% Computational Cost: add. (547754->320), mult. (1237788->409), div. (0->0), fcn. (972711->14), ass. (0->146)
t314 = qJD(2) ^ 2;
t301 = sin(pkin(11));
t304 = cos(pkin(11));
t284 = t301 * g(1) - t304 * g(2);
t285 = -t304 * g(1) - t301 * g(2);
t299 = -g(3) + qJDD(1);
t313 = cos(qJ(2));
t305 = cos(pkin(6));
t309 = sin(qJ(2));
t342 = t305 * t309;
t302 = sin(pkin(6));
t343 = t302 * t309;
t258 = t284 * t342 + t313 * t285 + t299 * t343;
t250 = -t314 * pkin(2) + qJDD(2) * qJ(3) + t258;
t300 = sin(pkin(12));
t272 = -t302 * t284 + t305 * t299;
t303 = cos(pkin(12));
t338 = qJD(2) * qJD(3);
t341 = t303 * t272 - 0.2e1 * t300 * t338;
t349 = pkin(3) * t303;
t229 = (-pkin(8) * qJDD(2) + t314 * t349 - t250) * t300 + t341;
t241 = t300 * t272 + (t250 + 0.2e1 * t338) * t303;
t337 = qJDD(2) * t303;
t296 = t303 ^ 2;
t344 = t296 * t314;
t232 = -pkin(3) * t344 + pkin(8) * t337 + t241;
t308 = sin(qJ(4));
t312 = cos(qJ(4));
t206 = t312 * t229 - t308 * t232;
t327 = t300 * t312 + t303 * t308;
t326 = -t300 * t308 + t303 * t312;
t274 = t326 * qJD(2);
t339 = t274 * qJD(4);
t266 = t327 * qJDD(2) + t339;
t275 = t327 * qJD(2);
t202 = (-t266 + t339) * pkin(9) + (t274 * t275 + qJDD(4)) * pkin(4) + t206;
t207 = t308 * t229 + t312 * t232;
t265 = -t275 * qJD(4) + t326 * qJDD(2);
t271 = qJD(4) * pkin(4) - t275 * pkin(9);
t273 = t274 ^ 2;
t204 = -t273 * pkin(4) + t265 * pkin(9) - qJD(4) * t271 + t207;
t307 = sin(qJ(5));
t311 = cos(qJ(5));
t199 = t307 * t202 + t311 * t204;
t256 = t307 * t274 + t311 * t275;
t223 = -t256 * qJD(5) + t311 * t265 - t307 * t266;
t255 = t311 * t274 - t307 * t275;
t238 = -mrSges(6,1) * t255 + mrSges(6,2) * t256;
t297 = qJD(4) + qJD(5);
t249 = t297 * mrSges(6,1) - t256 * mrSges(6,3);
t294 = qJDD(4) + qJDD(5);
t239 = -pkin(5) * t255 - pkin(10) * t256;
t293 = t297 ^ 2;
t196 = -t293 * pkin(5) + t294 * pkin(10) + t255 * t239 + t199;
t295 = t300 ^ 2;
t257 = -t309 * t285 + (t284 * t305 + t299 * t302) * t313;
t321 = qJDD(3) - t257;
t242 = (-pkin(2) - t349) * qJDD(2) + (-qJ(3) + (-t295 - t296) * pkin(8)) * t314 + t321;
t212 = -t265 * pkin(4) - t273 * pkin(9) + t275 * t271 + t242;
t224 = t255 * qJD(5) + t307 * t265 + t311 * t266;
t200 = (-t255 * t297 - t224) * pkin(10) + (t256 * t297 - t223) * pkin(5) + t212;
t306 = sin(qJ(6));
t310 = cos(qJ(6));
t193 = -t306 * t196 + t310 * t200;
t244 = -t306 * t256 + t310 * t297;
t210 = t244 * qJD(6) + t310 * t224 + t306 * t294;
t222 = qJDD(6) - t223;
t245 = t310 * t256 + t306 * t297;
t225 = -mrSges(7,1) * t244 + mrSges(7,2) * t245;
t251 = qJD(6) - t255;
t230 = -mrSges(7,2) * t251 + mrSges(7,3) * t244;
t189 = m(7) * t193 + t222 * mrSges(7,1) - t210 * mrSges(7,3) - t225 * t245 + t230 * t251;
t194 = t310 * t196 + t306 * t200;
t209 = -t245 * qJD(6) - t306 * t224 + t310 * t294;
t231 = mrSges(7,1) * t251 - mrSges(7,3) * t245;
t190 = m(7) * t194 - t222 * mrSges(7,2) + t209 * mrSges(7,3) + t225 * t244 - t231 * t251;
t333 = -t306 * t189 + t310 * t190;
t176 = m(6) * t199 - t294 * mrSges(6,2) + t223 * mrSges(6,3) + t255 * t238 - t297 * t249 + t333;
t198 = t311 * t202 - t307 * t204;
t248 = -t297 * mrSges(6,2) + t255 * mrSges(6,3);
t195 = -t294 * pkin(5) - t293 * pkin(10) + t256 * t239 - t198;
t322 = -m(7) * t195 + t209 * mrSges(7,1) - t210 * mrSges(7,2) + t244 * t230 - t245 * t231;
t185 = m(6) * t198 + t294 * mrSges(6,1) - t224 * mrSges(6,3) - t256 * t238 + t297 * t248 + t322;
t171 = t307 * t176 + t311 * t185;
t261 = -t274 * mrSges(5,1) + t275 * mrSges(5,2);
t269 = -qJD(4) * mrSges(5,2) + t274 * mrSges(5,3);
t168 = m(5) * t206 + qJDD(4) * mrSges(5,1) - t266 * mrSges(5,3) + qJD(4) * t269 - t275 * t261 + t171;
t270 = qJD(4) * mrSges(5,1) - t275 * mrSges(5,3);
t334 = t311 * t176 - t307 * t185;
t169 = m(5) * t207 - qJDD(4) * mrSges(5,2) + t265 * mrSges(5,3) - qJD(4) * t270 + t274 * t261 + t334;
t162 = t312 * t168 + t308 * t169;
t240 = -t300 * t250 + t341;
t347 = mrSges(4,2) * t300;
t325 = mrSges(4,3) * qJDD(2) + t314 * (-mrSges(4,1) * t303 + t347);
t160 = m(4) * t240 - t325 * t300 + t162;
t335 = -t308 * t168 + t312 * t169;
t161 = m(4) * t241 + t325 * t303 + t335;
t155 = t303 * t160 + t300 * t161;
t330 = Ifges(4,5) * t300 + Ifges(4,6) * t303;
t253 = Ifges(5,4) * t275 + Ifges(5,2) * t274 + Ifges(5,6) * qJD(4);
t254 = Ifges(5,1) * t275 + Ifges(5,4) * t274 + Ifges(5,5) * qJD(4);
t213 = Ifges(7,5) * t245 + Ifges(7,6) * t244 + Ifges(7,3) * t251;
t215 = Ifges(7,1) * t245 + Ifges(7,4) * t244 + Ifges(7,5) * t251;
t182 = -mrSges(7,1) * t195 + mrSges(7,3) * t194 + Ifges(7,4) * t210 + Ifges(7,2) * t209 + Ifges(7,6) * t222 - t213 * t245 + t215 * t251;
t214 = Ifges(7,4) * t245 + Ifges(7,2) * t244 + Ifges(7,6) * t251;
t183 = mrSges(7,2) * t195 - mrSges(7,3) * t193 + Ifges(7,1) * t210 + Ifges(7,4) * t209 + Ifges(7,5) * t222 + t213 * t244 - t214 * t251;
t234 = Ifges(6,4) * t256 + Ifges(6,2) * t255 + Ifges(6,6) * t297;
t235 = Ifges(6,1) * t256 + Ifges(6,4) * t255 + Ifges(6,5) * t297;
t320 = -mrSges(6,1) * t198 + mrSges(6,2) * t199 - Ifges(6,5) * t224 - Ifges(6,6) * t223 - Ifges(6,3) * t294 - pkin(5) * t322 - pkin(10) * t333 - t310 * t182 - t306 * t183 - t256 * t234 + t255 * t235;
t316 = -mrSges(5,1) * t206 + mrSges(5,2) * t207 - Ifges(5,5) * t266 - Ifges(5,6) * t265 - Ifges(5,3) * qJDD(4) - pkin(4) * t171 - t275 * t253 + t274 * t254 + t320;
t331 = Ifges(4,4) * t300 + Ifges(4,2) * t303;
t332 = Ifges(4,1) * t300 + Ifges(4,4) * t303;
t350 = -mrSges(4,1) * t240 + mrSges(4,2) * t241 - pkin(3) * t162 - (t300 * t331 - t303 * t332) * t314 + t316;
t142 = -pkin(2) * t155 + (Ifges(3,6) - t330) * qJDD(2) + t314 * Ifges(3,5) - mrSges(3,1) * t272 + mrSges(3,3) * t258 + t350;
t336 = -t300 * t160 + t303 * t161;
t153 = m(3) * t258 - t314 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t336;
t247 = -qJDD(2) * pkin(2) - t314 * qJ(3) + t321;
t178 = t310 * t189 + t306 * t190;
t324 = m(6) * t212 - t223 * mrSges(6,1) + t224 * mrSges(6,2) - t255 * t248 + t256 * t249 + t178;
t319 = m(5) * t242 - t265 * mrSges(5,1) + t266 * mrSges(5,2) - t274 * t269 + t275 * t270 + t324;
t317 = -m(4) * t247 + mrSges(4,1) * t337 - t319 + (t295 * t314 + t344) * mrSges(4,3);
t173 = (mrSges(3,1) - t347) * qJDD(2) + t317 - t314 * mrSges(3,2) + m(3) * t257;
t149 = t313 * t153 - t309 * t173;
t351 = pkin(7) * t149 + t142 * t313;
t345 = t173 * t313;
t340 = t314 * t330;
t154 = m(3) * t272 + t155;
t146 = t153 * t342 - t302 * t154 + t305 * t345;
t233 = Ifges(6,5) * t256 + Ifges(6,6) * t255 + Ifges(6,3) * t297;
t163 = mrSges(6,2) * t212 - mrSges(6,3) * t198 + Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t294 - pkin(10) * t178 - t306 * t182 + t310 * t183 + t255 * t233 - t297 * t234;
t318 = mrSges(7,1) * t193 - mrSges(7,2) * t194 + Ifges(7,5) * t210 + Ifges(7,6) * t209 + Ifges(7,3) * t222 + t245 * t214 - t244 * t215;
t164 = -mrSges(6,1) * t212 + mrSges(6,3) * t199 + Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t294 - pkin(5) * t178 - t256 * t233 + t297 * t235 - t318;
t252 = Ifges(5,5) * t275 + Ifges(5,6) * t274 + Ifges(5,3) * qJD(4);
t150 = -mrSges(5,1) * t242 + mrSges(5,3) * t207 + Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * qJDD(4) - pkin(4) * t324 + pkin(9) * t334 + qJD(4) * t254 + t307 * t163 + t311 * t164 - t275 * t252;
t156 = mrSges(5,2) * t242 - mrSges(5,3) * t206 + Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * qJDD(4) - pkin(9) * t171 - qJD(4) * t253 + t311 * t163 - t307 * t164 + t274 * t252;
t139 = -mrSges(4,1) * t247 + mrSges(4,3) * t241 - pkin(3) * t319 + pkin(8) * t335 + t331 * qJDD(2) + t312 * t150 + t308 * t156 - t300 * t340;
t140 = mrSges(4,2) * t247 - mrSges(4,3) * t240 - pkin(8) * t162 + t332 * qJDD(2) - t308 * t150 + t312 * t156 + t303 * t340;
t136 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t257 - mrSges(3,2) * t258 + t300 * t140 + t303 * t139 + pkin(2) * (-qJDD(2) * t347 + t317) + qJ(3) * t336;
t138 = mrSges(3,2) * t272 - mrSges(3,3) * t257 + Ifges(3,5) * qJDD(2) - t314 * Ifges(3,6) - qJ(3) * t155 - t300 * t139 + t303 * t140;
t323 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + pkin(1) * t146 + t305 * t136 + t138 * t343 + t351 * t302;
t147 = m(2) * t285 + t149;
t145 = t305 * t154 + (t153 * t309 + t345) * t302;
t143 = m(2) * t284 + t146;
t134 = mrSges(2,2) * t299 - mrSges(2,3) * t284 + t313 * t138 - t309 * t142 + (-t145 * t302 - t146 * t305) * pkin(7);
t133 = -mrSges(2,1) * t299 + mrSges(2,3) * t285 - pkin(1) * t145 - t302 * t136 + (t138 * t309 + t351) * t305;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t134 - t301 * t133 - qJ(1) * (t304 * t143 + t301 * t147), t134, t138, t140, t156, t163, t183; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t134 + t304 * t133 + qJ(1) * (-t301 * t143 + t304 * t147), t133, t142, t139, t150, t164, t182; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t323, t323, t136, t330 * qJDD(2) - t350, -t316, -t320, t318;];
m_new  = t1;
