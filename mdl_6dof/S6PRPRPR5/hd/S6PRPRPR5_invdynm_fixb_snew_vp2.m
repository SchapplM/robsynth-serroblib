% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 23:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:57:06
% EndTime: 2019-05-04 22:57:21
% DurationCPUTime: 10.68s
% Computational Cost: add. (171876->322), mult. (384749->394), div. (0->0), fcn. (277200->12), ass. (0->142)
t301 = sin(pkin(10));
t304 = cos(pkin(10));
t284 = g(1) * t301 - g(2) * t304;
t285 = -g(1) * t304 - g(2) * t301;
t299 = -g(3) + qJDD(1);
t310 = cos(qJ(2));
t305 = cos(pkin(6));
t308 = sin(qJ(2));
t347 = t305 * t308;
t302 = sin(pkin(6));
t348 = t302 * t308;
t240 = t284 * t347 + t310 * t285 + t299 * t348;
t312 = qJD(2) ^ 2;
t232 = -pkin(2) * t312 + qJDD(2) * qJ(3) + t240;
t300 = sin(pkin(11));
t269 = -t284 * t302 + t299 * t305;
t303 = cos(pkin(11));
t340 = qJD(2) * qJD(3);
t344 = t303 * t269 - 0.2e1 * t300 * t340;
t355 = pkin(3) * t303;
t209 = (-pkin(8) * qJDD(2) + t312 * t355 - t232) * t300 + t344;
t212 = t300 * t269 + (t232 + 0.2e1 * t340) * t303;
t338 = qJDD(2) * t303;
t295 = t303 ^ 2;
t349 = t295 * t312;
t210 = -pkin(3) * t349 + pkin(8) * t338 + t212;
t307 = sin(qJ(4));
t356 = cos(qJ(4));
t201 = t356 * t209 - t307 * t210;
t202 = t307 * t209 + t356 * t210;
t337 = t303 * t356;
t343 = qJD(2) * t300;
t274 = -qJD(2) * t337 + t307 * t343;
t327 = t356 * t300 + t303 * t307;
t275 = t327 * qJD(2);
t236 = Ifges(5,4) * t275 - Ifges(5,2) * t274 + Ifges(5,6) * qJD(4);
t247 = -mrSges(6,2) * t274 - mrSges(6,3) * t275;
t339 = qJDD(2) * t300;
t342 = qJD(4) * t275;
t257 = -qJDD(2) * t337 + t307 * t339 + t342;
t341 = t274 * qJD(4);
t258 = t327 * qJDD(2) - t341;
t266 = mrSges(6,1) * t274 - qJD(4) * mrSges(6,3);
t245 = pkin(4) * t274 - qJ(5) * t275;
t311 = qJD(4) ^ 2;
t199 = -qJDD(4) * pkin(4) - t311 * qJ(5) + t275 * t245 + qJDD(5) - t201;
t193 = (t274 * t275 - qJDD(4)) * pkin(9) + (t258 + t341) * pkin(5) + t199;
t268 = pkin(5) * t275 - qJD(4) * pkin(9);
t273 = t274 ^ 2;
t294 = t300 ^ 2;
t239 = -t308 * t285 + (t284 * t305 + t299 * t302) * t310;
t324 = qJDD(3) - t239;
t222 = (-pkin(2) - t355) * qJDD(2) + (-qJ(3) + (-t294 - t295) * pkin(8)) * t312 + t324;
t357 = -2 * qJD(5);
t314 = pkin(4) * t342 + t275 * t357 + (-t258 + t341) * qJ(5) + t222;
t196 = -pkin(5) * t273 - t268 * t275 + (pkin(4) + pkin(9)) * t257 + t314;
t306 = sin(qJ(6));
t309 = cos(qJ(6));
t190 = t193 * t309 - t196 * t306;
t259 = -qJD(4) * t306 + t274 * t309;
t221 = qJD(6) * t259 + qJDD(4) * t309 + t257 * t306;
t260 = qJD(4) * t309 + t274 * t306;
t225 = -mrSges(7,1) * t259 + mrSges(7,2) * t260;
t272 = qJD(6) + t275;
t230 = -mrSges(7,2) * t272 + mrSges(7,3) * t259;
t256 = qJDD(6) + t258;
t187 = m(7) * t190 + mrSges(7,1) * t256 - mrSges(7,3) * t221 - t225 * t260 + t230 * t272;
t191 = t193 * t306 + t196 * t309;
t220 = -qJD(6) * t260 - qJDD(4) * t306 + t257 * t309;
t231 = mrSges(7,1) * t272 - mrSges(7,3) * t260;
t188 = m(7) * t191 - mrSges(7,2) * t256 + mrSges(7,3) * t220 + t225 * t259 - t231 * t272;
t175 = t187 * t309 + t188 * t306;
t323 = -pkin(4) * t311 + qJDD(4) * qJ(5) - t245 * t274 + t202;
t195 = -pkin(5) * t257 - pkin(9) * t273 + ((2 * qJD(5)) + t268) * qJD(4) + t323;
t213 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t272;
t215 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t272;
t178 = -mrSges(7,1) * t195 + mrSges(7,3) * t191 + Ifges(7,4) * t221 + Ifges(7,2) * t220 + Ifges(7,6) * t256 - t213 * t260 + t215 * t272;
t214 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t272;
t179 = mrSges(7,2) * t195 - mrSges(7,3) * t190 + Ifges(7,1) * t221 + Ifges(7,4) * t220 + Ifges(7,5) * t256 + t213 * t259 - t214 * t272;
t197 = qJD(4) * t357 - t323;
t233 = Ifges(6,5) * qJD(4) - Ifges(6,6) * t275 + Ifges(6,3) * t274;
t320 = -mrSges(6,2) * t199 + mrSges(6,3) * t197 - Ifges(6,1) * qJDD(4) + Ifges(6,4) * t258 - Ifges(6,5) * t257 + pkin(9) * t175 + t306 * t178 - t309 * t179 + t275 * t233;
t192 = -m(7) * t195 + mrSges(7,1) * t220 - t221 * mrSges(7,2) + t230 * t259 - t260 * t231;
t267 = mrSges(6,1) * t275 + qJD(4) * mrSges(6,2);
t321 = -m(6) * t197 + qJDD(4) * mrSges(6,3) + qJD(4) * t267 - t192;
t325 = -m(6) * t199 - t258 * mrSges(6,1) - t275 * t247 - t175;
t235 = Ifges(6,4) * qJD(4) - Ifges(6,2) * t275 + Ifges(6,6) * t274;
t345 = Ifges(5,1) * t275 - Ifges(5,4) * t274 + Ifges(5,5) * qJD(4) - t235;
t362 = t345 * t274 - mrSges(5,2) * t202 + pkin(4) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t266 + t325) + qJ(5) * (-mrSges(6,1) * t257 - t247 * t274 + t321) + mrSges(5,1) * t201 + t275 * t236 - Ifges(5,6) * t257 + Ifges(5,5) * t258 + Ifges(5,3) * qJDD(4) - t320;
t246 = mrSges(5,1) * t274 + mrSges(5,2) * t275;
t264 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t274;
t169 = m(5) * t201 - mrSges(5,3) * t258 - t246 * t275 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t264 - t266) * qJD(4) + t325;
t265 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t275;
t182 = m(5) * t202 - qJDD(4) * mrSges(5,2) - qJD(4) * t265 + (-t246 - t247) * t274 + (-mrSges(5,3) - mrSges(6,1)) * t257 + t321;
t167 = t356 * t169 + t307 * t182;
t211 = -t232 * t300 + t344;
t352 = mrSges(4,2) * t300;
t328 = mrSges(4,3) * qJDD(2) + t312 * (-mrSges(4,1) * t303 + t352);
t165 = m(4) * t211 - t328 * t300 + t167;
t334 = -t169 * t307 + t356 * t182;
t166 = m(4) * t212 + t328 * t303 + t334;
t160 = t303 * t165 + t300 * t166;
t331 = Ifges(4,5) * t300 + Ifges(4,6) * t303;
t332 = Ifges(4,4) * t300 + Ifges(4,2) * t303;
t333 = Ifges(4,1) * t300 + Ifges(4,4) * t303;
t359 = qJD(2) * t303;
t360 = -mrSges(4,1) * t211 + mrSges(4,2) * t212 - pkin(3) * t167 - (t332 * t343 - t333 * t359) * qJD(2) - t362;
t146 = (Ifges(3,6) - t331) * qJDD(2) + t312 * Ifges(3,5) - mrSges(3,1) * t269 + mrSges(3,3) * t240 - pkin(2) * t160 + t360;
t335 = -t165 * t300 + t303 * t166;
t158 = m(3) * t240 - mrSges(3,1) * t312 - qJDD(2) * mrSges(3,2) + t335;
t229 = -qJDD(2) * pkin(2) - qJ(3) * t312 + t324;
t176 = -t306 * t187 + t309 * t188;
t204 = pkin(4) * t257 + t314;
t174 = m(6) * t204 - t257 * mrSges(6,2) - t258 * mrSges(6,3) - t274 * t266 - t275 * t267 + t176;
t318 = m(5) * t222 + t257 * mrSges(5,1) + mrSges(5,2) * t258 + t274 * t264 + t265 * t275 + t174;
t316 = -m(4) * t229 + mrSges(4,1) * t338 - t318 + (t294 * t312 + t349) * mrSges(4,3);
t171 = (mrSges(3,1) - t352) * qJDD(2) + m(3) * t239 - mrSges(3,2) * t312 + t316;
t154 = t310 * t158 - t171 * t308;
t361 = pkin(7) * t154 + t146 * t310;
t353 = Ifges(5,4) + Ifges(6,6);
t350 = t171 * t310;
t237 = Ifges(6,1) * qJD(4) - Ifges(6,4) * t275 + Ifges(6,5) * t274;
t346 = -Ifges(5,5) * t275 + Ifges(5,6) * t274 - Ifges(5,3) * qJD(4) - t237;
t159 = m(3) * t269 + t160;
t151 = t158 * t347 - t159 * t302 + t305 * t350;
t319 = -mrSges(6,1) * t197 + mrSges(6,2) * t204 - pkin(5) * t192 - pkin(9) * t176 - t309 * t178 - t306 * t179;
t155 = -mrSges(5,1) * t222 + mrSges(5,3) * t202 - pkin(4) * t174 + t346 * t275 + t353 * t258 + (-Ifges(5,2) - Ifges(6,3)) * t257 + (Ifges(5,6) - Ifges(6,5)) * qJDD(4) + t345 * qJD(4) + t319;
t322 = mrSges(7,1) * t190 - mrSges(7,2) * t191 + Ifges(7,5) * t221 + Ifges(7,6) * t220 + Ifges(7,3) * t256 + t260 * t214 - t259 * t215;
t317 = mrSges(6,1) * t199 - mrSges(6,3) * t204 + pkin(5) * t175 + t322;
t161 = t346 * t274 + (Ifges(5,1) + Ifges(6,2)) * t258 - t353 * t257 + (Ifges(5,5) - Ifges(6,4)) * qJDD(4) + (-t236 + t233) * qJD(4) + mrSges(5,2) * t222 - mrSges(5,3) * t201 - qJ(5) * t174 + t317;
t280 = t331 * qJD(2);
t144 = -mrSges(4,1) * t229 + mrSges(4,3) * t212 - pkin(3) * t318 + pkin(8) * t334 + t332 * qJDD(2) + t356 * t155 + t307 * t161 - t280 * t343;
t147 = mrSges(4,2) * t229 - mrSges(4,3) * t211 - pkin(8) * t167 + t333 * qJDD(2) - t307 * t155 + t356 * t161 + t280 * t359;
t141 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t239 - mrSges(3,2) * t240 + t300 * t147 + t303 * t144 + pkin(2) * (-mrSges(4,2) * t339 + t316) + qJ(3) * t335;
t143 = mrSges(3,2) * t269 - mrSges(3,3) * t239 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t312 - qJ(3) * t160 - t144 * t300 + t147 * t303;
t326 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + pkin(1) * t151 + t305 * t141 + t143 * t348 + t361 * t302;
t152 = m(2) * t285 + t154;
t150 = t159 * t305 + (t158 * t308 + t350) * t302;
t148 = m(2) * t284 + t151;
t139 = mrSges(2,2) * t299 - mrSges(2,3) * t284 + t143 * t310 - t146 * t308 + (-t150 * t302 - t151 * t305) * pkin(7);
t138 = -mrSges(2,1) * t299 + mrSges(2,3) * t285 - pkin(1) * t150 - t141 * t302 + (t143 * t308 + t361) * t305;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t304 * t139 - t301 * t138 - qJ(1) * (t148 * t304 + t152 * t301), t139, t143, t147, t161, -t274 * t235 - t320, t179; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t301 * t139 + t304 * t138 + qJ(1) * (-t148 * t301 + t152 * t304), t138, t146, t144, t155, Ifges(6,4) * qJDD(4) - Ifges(6,2) * t258 + Ifges(6,6) * t257 - qJD(4) * t233 + t274 * t237 - t317, t178; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t326, t326, t141, t331 * qJDD(2) - t360, t362, Ifges(6,5) * qJDD(4) - Ifges(6,6) * t258 + Ifges(6,3) * t257 + qJD(4) * t235 + t275 * t237 - t319, t322;];
m_new  = t1;
