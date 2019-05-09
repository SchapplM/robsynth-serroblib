% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:23:46
% EndTime: 2019-05-05 10:24:32
% DurationCPUTime: 35.06s
% Computational Cost: add. (663199->341), mult. (1381381->439), div. (0->0), fcn. (1038857->14), ass. (0->146)
t305 = sin(pkin(12));
t307 = cos(pkin(12));
t290 = g(1) * t305 - g(2) * t307;
t291 = -g(1) * t307 - g(2) * t305;
t304 = -g(3) + qJDD(1);
t318 = cos(qJ(2));
t308 = cos(pkin(6));
t313 = sin(qJ(2));
t340 = t308 * t313;
t306 = sin(pkin(6));
t341 = t306 * t313;
t262 = t290 * t340 + t318 * t291 + t304 * t341;
t319 = qJD(2) ^ 2;
t256 = -pkin(2) * t319 + qJDD(2) * pkin(8) + t262;
t273 = -t290 * t306 + t304 * t308;
t312 = sin(qJ(3));
t317 = cos(qJ(3));
t242 = -t312 * t256 + t317 * t273;
t337 = qJD(2) * qJD(3);
t336 = t317 * t337;
t287 = qJDD(2) * t312 + t336;
t229 = (-t287 + t336) * pkin(9) + (t312 * t317 * t319 + qJDD(3)) * pkin(3) + t242;
t243 = t317 * t256 + t312 * t273;
t288 = qJDD(2) * t317 - t312 * t337;
t339 = qJD(2) * t312;
t295 = qJD(3) * pkin(3) - pkin(9) * t339;
t303 = t317 ^ 2;
t231 = -pkin(3) * t303 * t319 + pkin(9) * t288 - qJD(3) * t295 + t243;
t311 = sin(qJ(4));
t316 = cos(qJ(4));
t207 = t316 * t229 - t311 * t231;
t279 = (-t311 * t312 + t316 * t317) * qJD(2);
t252 = qJD(4) * t279 + t287 * t316 + t288 * t311;
t280 = (t311 * t317 + t312 * t316) * qJD(2);
t301 = qJDD(3) + qJDD(4);
t302 = qJD(3) + qJD(4);
t203 = (t279 * t302 - t252) * pkin(10) + (t279 * t280 + t301) * pkin(4) + t207;
t208 = t311 * t229 + t316 * t231;
t251 = -qJD(4) * t280 - t287 * t311 + t288 * t316;
t272 = pkin(4) * t302 - pkin(10) * t280;
t275 = t279 ^ 2;
t205 = -pkin(4) * t275 + pkin(10) * t251 - t272 * t302 + t208;
t310 = sin(qJ(5));
t315 = cos(qJ(5));
t200 = t310 * t203 + t315 * t205;
t266 = t279 * t310 + t280 * t315;
t223 = -qJD(5) * t266 + t251 * t315 - t252 * t310;
t265 = t279 * t315 - t280 * t310;
t239 = -mrSges(6,1) * t265 + mrSges(6,2) * t266;
t300 = qJD(5) + t302;
t254 = mrSges(6,1) * t300 - mrSges(6,3) * t266;
t299 = qJDD(5) + t301;
t240 = -pkin(5) * t265 - pkin(11) * t266;
t298 = t300 ^ 2;
t197 = -pkin(5) * t298 + pkin(11) * t299 + t240 * t265 + t200;
t261 = -t313 * t291 + (t290 * t308 + t304 * t306) * t318;
t326 = -qJDD(2) * pkin(2) - t261;
t241 = -t288 * pkin(3) + t295 * t339 + (-pkin(9) * t303 - pkin(8)) * t319 + t326;
t213 = -t251 * pkin(4) - t275 * pkin(10) + t280 * t272 + t241;
t224 = qJD(5) * t265 + t251 * t310 + t252 * t315;
t201 = (-t265 * t300 - t224) * pkin(11) + (t266 * t300 - t223) * pkin(5) + t213;
t309 = sin(qJ(6));
t314 = cos(qJ(6));
t194 = -t197 * t309 + t201 * t314;
t245 = -t266 * t309 + t300 * t314;
t211 = qJD(6) * t245 + t224 * t314 + t299 * t309;
t221 = qJDD(6) - t223;
t246 = t266 * t314 + t300 * t309;
t230 = -mrSges(7,1) * t245 + mrSges(7,2) * t246;
t260 = qJD(6) - t265;
t232 = -mrSges(7,2) * t260 + mrSges(7,3) * t245;
t190 = m(7) * t194 + mrSges(7,1) * t221 - mrSges(7,3) * t211 - t230 * t246 + t232 * t260;
t195 = t197 * t314 + t201 * t309;
t210 = -qJD(6) * t246 - t224 * t309 + t299 * t314;
t233 = mrSges(7,1) * t260 - mrSges(7,3) * t246;
t191 = m(7) * t195 - mrSges(7,2) * t221 + mrSges(7,3) * t210 + t230 * t245 - t233 * t260;
t332 = -t190 * t309 + t314 * t191;
t177 = m(6) * t200 - mrSges(6,2) * t299 + mrSges(6,3) * t223 + t239 * t265 - t254 * t300 + t332;
t199 = t203 * t315 - t205 * t310;
t253 = -mrSges(6,2) * t300 + mrSges(6,3) * t265;
t196 = -pkin(5) * t299 - pkin(11) * t298 + t240 * t266 - t199;
t327 = -m(7) * t196 + t210 * mrSges(7,1) - mrSges(7,2) * t211 + t245 * t232 - t233 * t246;
t186 = m(6) * t199 + mrSges(6,1) * t299 - mrSges(6,3) * t224 - t239 * t266 + t253 * t300 + t327;
t172 = t310 * t177 + t315 * t186;
t267 = -mrSges(5,1) * t279 + mrSges(5,2) * t280;
t270 = -mrSges(5,2) * t302 + mrSges(5,3) * t279;
t169 = m(5) * t207 + mrSges(5,1) * t301 - mrSges(5,3) * t252 - t267 * t280 + t270 * t302 + t172;
t271 = mrSges(5,1) * t302 - mrSges(5,3) * t280;
t333 = t315 * t177 - t186 * t310;
t170 = m(5) * t208 - mrSges(5,2) * t301 + mrSges(5,3) * t251 + t267 * t279 - t271 * t302 + t333;
t163 = t316 * t169 + t311 * t170;
t286 = (-mrSges(4,1) * t317 + mrSges(4,2) * t312) * qJD(2);
t338 = qJD(2) * t317;
t293 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t338;
t161 = m(4) * t242 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t287 + qJD(3) * t293 - t286 * t339 + t163;
t292 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t339;
t334 = -t169 * t311 + t316 * t170;
t162 = m(4) * t243 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t288 - qJD(3) * t292 + t286 * t338 + t334;
t156 = t317 * t161 + t312 * t162;
t277 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t312 + Ifges(4,2) * t317) * qJD(2);
t278 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t312 + Ifges(4,4) * t317) * qJD(2);
t258 = Ifges(5,4) * t280 + Ifges(5,2) * t279 + Ifges(5,6) * t302;
t259 = Ifges(5,1) * t280 + Ifges(5,4) * t279 + Ifges(5,5) * t302;
t214 = Ifges(7,5) * t246 + Ifges(7,6) * t245 + Ifges(7,3) * t260;
t216 = Ifges(7,1) * t246 + Ifges(7,4) * t245 + Ifges(7,5) * t260;
t183 = -mrSges(7,1) * t196 + mrSges(7,3) * t195 + Ifges(7,4) * t211 + Ifges(7,2) * t210 + Ifges(7,6) * t221 - t214 * t246 + t216 * t260;
t215 = Ifges(7,4) * t246 + Ifges(7,2) * t245 + Ifges(7,6) * t260;
t184 = mrSges(7,2) * t196 - mrSges(7,3) * t194 + Ifges(7,1) * t211 + Ifges(7,4) * t210 + Ifges(7,5) * t221 + t214 * t245 - t215 * t260;
t235 = Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t300;
t236 = Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t300;
t325 = -mrSges(6,1) * t199 + mrSges(6,2) * t200 - Ifges(6,5) * t224 - Ifges(6,6) * t223 - Ifges(6,3) * t299 - pkin(5) * t327 - pkin(11) * t332 - t314 * t183 - t309 * t184 - t266 * t235 + t265 * t236;
t322 = -mrSges(5,1) * t207 + mrSges(5,2) * t208 - Ifges(5,5) * t252 - Ifges(5,6) * t251 - Ifges(5,3) * t301 - pkin(4) * t172 - t280 * t258 + t279 * t259 + t325;
t345 = mrSges(4,1) * t242 - mrSges(4,2) * t243 + Ifges(4,5) * t287 + Ifges(4,6) * t288 + Ifges(4,3) * qJDD(3) + pkin(3) * t163 + (t277 * t312 - t278 * t317) * qJD(2) - t322;
t143 = -mrSges(3,1) * t273 + mrSges(3,3) * t262 + t319 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t156 - t345;
t335 = -t161 * t312 + t317 * t162;
t154 = m(3) * t262 - mrSges(3,1) * t319 - qJDD(2) * mrSges(3,2) + t335;
t255 = -t319 * pkin(8) + t326;
t179 = t314 * t190 + t309 * t191;
t329 = m(6) * t213 - t223 * mrSges(6,1) + t224 * mrSges(6,2) - t265 * t253 + t266 * t254 + t179;
t324 = m(5) * t241 - t251 * mrSges(5,1) + mrSges(5,2) * t252 - t279 * t270 + t271 * t280 + t329;
t321 = -m(4) * t255 + t288 * mrSges(4,1) - mrSges(4,2) * t287 - t292 * t339 + t293 * t338 - t324;
t174 = m(3) * t261 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t319 + t321;
t150 = t318 * t154 - t174 * t313;
t346 = pkin(7) * t150 + t143 * t318;
t342 = t174 * t318;
t155 = m(3) * t273 + t156;
t147 = t154 * t340 - t155 * t306 + t308 * t342;
t234 = Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t300;
t164 = mrSges(6,2) * t213 - mrSges(6,3) * t199 + Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t299 - pkin(11) * t179 - t183 * t309 + t184 * t314 + t234 * t265 - t235 * t300;
t323 = mrSges(7,1) * t194 - mrSges(7,2) * t195 + Ifges(7,5) * t211 + Ifges(7,6) * t210 + Ifges(7,3) * t221 + t215 * t246 - t216 * t245;
t165 = -mrSges(6,1) * t213 + mrSges(6,3) * t200 + Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t299 - pkin(5) * t179 - t234 * t266 + t236 * t300 - t323;
t257 = Ifges(5,5) * t280 + Ifges(5,6) * t279 + Ifges(5,3) * t302;
t151 = -mrSges(5,1) * t241 + mrSges(5,3) * t208 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * t301 - pkin(4) * t329 + pkin(10) * t333 + t310 * t164 + t315 * t165 - t280 * t257 + t302 * t259;
t157 = mrSges(5,2) * t241 - mrSges(5,3) * t207 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * t301 - pkin(10) * t172 + t164 * t315 - t165 * t310 + t257 * t279 - t258 * t302;
t276 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t312 + Ifges(4,6) * t317) * qJD(2);
t140 = -mrSges(4,1) * t255 + mrSges(4,3) * t243 + Ifges(4,4) * t287 + Ifges(4,2) * t288 + Ifges(4,6) * qJDD(3) - pkin(3) * t324 + pkin(9) * t334 + qJD(3) * t278 + t316 * t151 + t311 * t157 - t276 * t339;
t141 = mrSges(4,2) * t255 - mrSges(4,3) * t242 + Ifges(4,1) * t287 + Ifges(4,4) * t288 + Ifges(4,5) * qJDD(3) - pkin(9) * t163 - qJD(3) * t277 - t151 * t311 + t157 * t316 + t276 * t338;
t137 = mrSges(3,1) * t261 - mrSges(3,2) * t262 + Ifges(3,3) * qJDD(2) + pkin(2) * t321 + pkin(8) * t335 + t317 * t140 + t312 * t141;
t139 = mrSges(3,2) * t273 - mrSges(3,3) * t261 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t319 - pkin(8) * t156 - t140 * t312 + t141 * t317;
t328 = mrSges(2,1) * t290 - mrSges(2,2) * t291 + pkin(1) * t147 + t308 * t137 + t139 * t341 + t346 * t306;
t148 = m(2) * t291 + t150;
t146 = t308 * t155 + (t154 * t313 + t342) * t306;
t144 = m(2) * t290 + t147;
t135 = mrSges(2,2) * t304 - mrSges(2,3) * t290 + t318 * t139 - t313 * t143 + (-t146 * t306 - t147 * t308) * pkin(7);
t134 = -mrSges(2,1) * t304 + mrSges(2,3) * t291 - pkin(1) * t146 - t306 * t137 + (t139 * t313 + t346) * t308;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t307 * t135 - t305 * t134 - qJ(1) * (t144 * t307 + t148 * t305), t135, t139, t141, t157, t164, t184; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t305 * t135 + t307 * t134 + qJ(1) * (-t144 * t305 + t148 * t307), t134, t143, t140, t151, t165, t183; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t328, t328, t137, t345, -t322, -t325, t323;];
m_new  = t1;
