% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 00:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:17:20
% EndTime: 2019-05-07 00:17:56
% DurationCPUTime: 18.01s
% Computational Cost: add. (324034->391), mult. (672943->470), div. (0->0), fcn. (436020->10), ass. (0->148)
t380 = -2 * qJD(3);
t340 = sin(qJ(1));
t345 = cos(qJ(1));
t321 = -g(1) * t345 - g(2) * t340;
t347 = qJD(1) ^ 2;
t292 = -pkin(1) * t347 + qJDD(1) * pkin(7) + t321;
t339 = sin(qJ(2));
t344 = cos(qJ(2));
t273 = -g(3) * t344 - t292 * t339;
t274 = -g(3) * t339 + t292 * t344;
t287 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t339 + Ifges(3,4) * t344) * qJD(1);
t307 = (mrSges(4,2) * t344 - mrSges(4,3) * t339) * qJD(1);
t371 = qJD(1) * qJD(2);
t369 = t344 * t371;
t309 = qJDD(1) * t339 + t369;
t368 = t339 * t371;
t310 = qJDD(1) * t344 - t368;
t372 = qJD(1) * t344;
t317 = -mrSges(4,1) * t372 - qJD(2) * mrSges(4,3);
t327 = t339 * qJD(1);
t306 = (-pkin(2) * t344 - qJ(3) * t339) * qJD(1);
t346 = qJD(2) ^ 2;
t251 = pkin(2) * t346 - qJDD(2) * qJ(3) + qJD(2) * t380 - t306 * t372 - t274;
t319 = pkin(3) * t327 - qJD(2) * pkin(8);
t335 = t344 ^ 2;
t240 = -pkin(8) * t335 * t347 + pkin(3) * t310 + qJD(2) * t319 - t251;
t338 = sin(qJ(4));
t343 = cos(qJ(4));
t305 = qJD(2) * t343 - t338 * t372;
t265 = -qJD(4) * t305 - qJDD(2) * t338 - t310 * t343;
t304 = -qJD(2) * t338 - t343 * t372;
t266 = qJD(4) * t304 + qJDD(2) * t343 - t310 * t338;
t324 = t327 + qJD(4);
t271 = -mrSges(5,2) * t324 + mrSges(5,3) * t304;
t272 = mrSges(5,1) * t324 - mrSges(5,3) * t305;
t275 = pkin(4) * t324 - pkin(9) * t305;
t302 = t304 ^ 2;
t221 = -pkin(4) * t265 - pkin(9) * t302 + t275 * t305 + t240;
t337 = sin(qJ(5));
t342 = cos(qJ(5));
t269 = t304 * t337 + t305 * t342;
t231 = -qJD(5) * t269 + t265 * t342 - t266 * t337;
t268 = t304 * t342 - t305 * t337;
t232 = qJD(5) * t268 + t265 * t337 + t266 * t342;
t322 = qJD(5) + t324;
t254 = -mrSges(6,2) * t322 + mrSges(6,3) * t268;
t255 = mrSges(6,1) * t322 - mrSges(6,3) * t269;
t256 = pkin(5) * t322 - pkin(10) * t269;
t267 = t268 ^ 2;
t205 = -pkin(5) * t231 - pkin(10) * t267 + t256 * t269 + t221;
t336 = sin(qJ(6));
t341 = cos(qJ(6));
t247 = t268 * t336 + t269 * t341;
t209 = -qJD(6) * t247 + t231 * t341 - t232 * t336;
t246 = t268 * t341 - t269 * t336;
t210 = qJD(6) * t246 + t231 * t336 + t232 * t341;
t314 = qJD(6) + t322;
t235 = -mrSges(7,2) * t314 + mrSges(7,3) * t246;
t236 = mrSges(7,1) * t314 - mrSges(7,3) * t247;
t361 = m(7) * t205 - t209 * mrSges(7,1) + mrSges(7,2) * t210 - t246 * t235 + t236 * t247;
t354 = m(6) * t221 - t231 * mrSges(6,1) + mrSges(6,2) * t232 - t268 * t254 + t255 * t269 + t361;
t195 = -m(5) * t240 + t265 * mrSges(5,1) - mrSges(5,2) * t266 + t304 * t271 - t272 * t305 - t354;
t318 = mrSges(4,1) * t327 + qJD(2) * mrSges(4,2);
t351 = -m(4) * t251 + qJDD(2) * mrSges(4,3) + qJD(2) * t318 + t307 * t372 - t195;
t320 = t340 * g(1) - g(2) * t345;
t364 = -qJDD(1) * pkin(1) - t320;
t355 = pkin(2) * t368 + t327 * t380 + (-t309 - t369) * qJ(3) + t364;
t234 = -t319 * t327 + (-pkin(3) * t335 - pkin(7)) * t347 + (-pkin(2) - pkin(8)) * t310 + t355;
t253 = -qJDD(2) * pkin(2) - t346 * qJ(3) + t306 * t327 + qJDD(3) - t273;
t241 = (-t339 * t344 * t347 - qJDD(2)) * pkin(8) + (t309 - t369) * pkin(3) + t253;
t218 = -t338 * t234 + t241 * t343;
t303 = qJDD(4) + t309;
t213 = (t304 * t324 - t266) * pkin(9) + (t304 * t305 + t303) * pkin(4) + t218;
t219 = t234 * t343 + t241 * t338;
t215 = -pkin(4) * t302 + pkin(9) * t265 - t275 * t324 + t219;
t202 = t213 * t342 - t337 * t215;
t294 = qJDD(5) + t303;
t199 = (t268 * t322 - t232) * pkin(10) + (t268 * t269 + t294) * pkin(5) + t202;
t203 = t213 * t337 + t215 * t342;
t200 = -pkin(5) * t267 + pkin(10) * t231 - t256 * t322 + t203;
t198 = t199 * t336 + t200 * t341;
t222 = Ifges(7,5) * t247 + Ifges(7,6) * t246 + Ifges(7,3) * t314;
t224 = Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t314;
t284 = qJDD(6) + t294;
t185 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t210 + Ifges(7,2) * t209 + Ifges(7,6) * t284 - t222 * t247 + t224 * t314;
t197 = t199 * t341 - t200 * t336;
t223 = Ifges(7,4) * t247 + Ifges(7,2) * t246 + Ifges(7,6) * t314;
t186 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t210 + Ifges(7,4) * t209 + Ifges(7,5) * t284 + t222 * t246 - t223 * t314;
t242 = Ifges(6,5) * t269 + Ifges(6,6) * t268 + Ifges(6,3) * t322;
t244 = Ifges(6,1) * t269 + Ifges(6,4) * t268 + Ifges(6,5) * t322;
t226 = -mrSges(7,1) * t246 + mrSges(7,2) * t247;
t193 = m(7) * t197 + mrSges(7,1) * t284 - mrSges(7,3) * t210 - t226 * t247 + t235 * t314;
t194 = m(7) * t198 - mrSges(7,2) * t284 + mrSges(7,3) * t209 + t226 * t246 - t236 * t314;
t365 = -t193 * t336 + t194 * t341;
t170 = -mrSges(6,1) * t221 + mrSges(6,3) * t203 + Ifges(6,4) * t232 + Ifges(6,2) * t231 + Ifges(6,6) * t294 - pkin(5) * t361 + pkin(10) * t365 + t341 * t185 + t336 * t186 - t269 * t242 + t322 * t244;
t184 = t193 * t341 + t194 * t336;
t243 = Ifges(6,4) * t269 + Ifges(6,2) * t268 + Ifges(6,6) * t322;
t171 = mrSges(6,2) * t221 - mrSges(6,3) * t202 + Ifges(6,1) * t232 + Ifges(6,4) * t231 + Ifges(6,5) * t294 - pkin(10) * t184 - t185 * t336 + t186 * t341 + t242 * t268 - t243 * t322;
t257 = Ifges(5,5) * t305 + Ifges(5,6) * t304 + Ifges(5,3) * t324;
t259 = Ifges(5,1) * t305 + Ifges(5,4) * t304 + Ifges(5,5) * t324;
t248 = -mrSges(6,1) * t268 + mrSges(6,2) * t269;
t181 = m(6) * t202 + mrSges(6,1) * t294 - mrSges(6,3) * t232 - t248 * t269 + t254 * t322 + t184;
t182 = m(6) * t203 - mrSges(6,2) * t294 + mrSges(6,3) * t231 + t248 * t268 - t255 * t322 + t365;
t366 = -t181 * t337 + t182 * t342;
t157 = -mrSges(5,1) * t240 + mrSges(5,3) * t219 + Ifges(5,4) * t266 + Ifges(5,2) * t265 + Ifges(5,6) * t303 - pkin(4) * t354 + pkin(9) * t366 + t342 * t170 + t337 * t171 - t305 * t257 + t324 * t259;
t177 = t181 * t342 + t182 * t337;
t258 = Ifges(5,4) * t305 + Ifges(5,2) * t304 + Ifges(5,6) * t324;
t159 = mrSges(5,2) * t240 - mrSges(5,3) * t218 + Ifges(5,1) * t266 + Ifges(5,4) * t265 + Ifges(5,5) * t303 - pkin(9) * t177 - t170 * t337 + t171 * t342 + t257 * t304 - t258 * t324;
t270 = -mrSges(5,1) * t304 + mrSges(5,2) * t305;
t174 = m(5) * t218 + mrSges(5,1) * t303 - mrSges(5,3) * t266 - t270 * t305 + t271 * t324 + t177;
t175 = m(5) * t219 - mrSges(5,2) * t303 + mrSges(5,3) * t265 + t270 * t304 - t272 * t324 + t366;
t168 = t343 * t174 + t338 * t175;
t289 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t339 - Ifges(4,6) * t344) * qJD(1);
t357 = -mrSges(4,2) * t253 + mrSges(4,3) * t251 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t309 + Ifges(4,5) * t310 + pkin(8) * t168 + t338 * t157 - t159 * t343 - t289 * t372;
t360 = -m(4) * t253 - mrSges(4,1) * t309 - t168;
t288 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t339 - Ifges(4,3) * t344) * qJD(1);
t373 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t339 + Ifges(3,2) * t344) * qJD(1) - t288;
t379 = qJD(1) * (-t344 * t287 + t339 * t373) + mrSges(3,1) * t273 - mrSges(3,2) * t274 + Ifges(3,5) * t309 + Ifges(3,6) * t310 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t317 - t307 * t327 + t360) + qJ(3) * (t310 * mrSges(4,1) + t351) - t357;
t377 = t347 * pkin(7);
t376 = Ifges(3,4) + Ifges(4,6);
t169 = -t174 * t338 + t175 * t343;
t290 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t339 - Ifges(4,5) * t344) * qJD(1);
t374 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t339 + Ifges(3,6) * t344) * qJD(1) + t290;
t308 = (-mrSges(3,1) * t344 + mrSges(3,2) * t339) * qJD(1);
t316 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t372;
t165 = m(3) * t273 - t309 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t316 - t317) * qJD(2) + (-t307 - t308) * t327 + t360;
t315 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t327;
t188 = t351 + t308 * t372 + (mrSges(3,3) + mrSges(4,1)) * t310 - qJDD(2) * mrSges(3,2) - qJD(2) * t315 + m(3) * t274;
t367 = -t165 * t339 + t188 * t344;
t249 = -t310 * pkin(2) + t355 - t377;
t363 = -m(4) * t249 - mrSges(4,2) * t310 + t318 * t327 - t169;
t166 = -t309 * mrSges(4,3) + t317 * t372 - t363;
t291 = t364 - t377;
t356 = -mrSges(4,1) * t251 + mrSges(4,2) * t249 - pkin(3) * t195 - pkin(8) * t169 - t343 * t157 - t338 * t159;
t154 = -mrSges(3,1) * t291 + mrSges(3,3) * t274 - pkin(2) * t166 + (Ifges(3,2) + Ifges(4,3)) * t310 + t376 * t309 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t287 - t289) * qJD(2) - t374 * t327 + t356;
t358 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t210 - Ifges(7,6) * t209 - Ifges(7,3) * t284 - t223 * t247 + t246 * t224;
t353 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t232 - Ifges(6,6) * t231 - Ifges(6,3) * t294 - pkin(5) * t184 - t243 * t269 + t268 * t244 + t358;
t349 = -mrSges(5,1) * t218 + mrSges(5,2) * t219 - Ifges(5,5) * t266 - Ifges(5,6) * t265 - Ifges(5,3) * t303 - pkin(4) * t177 - t258 * t305 + t304 * t259 + t353;
t348 = -mrSges(4,1) * t253 + mrSges(4,3) * t249 - pkin(3) * t168 + t349;
t156 = t374 * t372 - t348 - qJ(3) * t166 + t376 * t310 + (Ifges(3,1) + Ifges(4,2)) * t309 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t373 * qJD(2) + mrSges(3,2) * t291 - mrSges(3,3) * t273;
t352 = -m(3) * t291 + t316 * t372 + t310 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t309 + (-t315 * t339 - t317 * t344) * qJD(1) + t363;
t359 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t367 + t154 * t344 + t156 * t339;
t163 = m(2) * t320 + qJDD(1) * mrSges(2,1) - t347 * mrSges(2,2) + t352;
t162 = t165 * t344 + t188 * t339;
t160 = m(2) * t321 - mrSges(2,1) * t347 - qJDD(1) * mrSges(2,2) + t367;
t152 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t379;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t347 - pkin(7) * t162 - t154 * t339 + t156 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t151 - t340 * t152 - pkin(6) * (t160 * t340 + t163 * t345), t151, t156, -t288 * t327 - t357, t159, t171, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t340 * t151 + t345 * t152 + pkin(6) * (t160 * t345 - t163 * t340), t152, t154, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t309 - Ifges(4,6) * t310 - qJD(2) * t288 - t290 * t372 + t348, t157, t170, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, t379, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t309 - Ifges(4,3) * t310 + qJD(2) * t289 + t290 * t327 - t356, -t349, -t353, -t358;];
m_new  = t1;
