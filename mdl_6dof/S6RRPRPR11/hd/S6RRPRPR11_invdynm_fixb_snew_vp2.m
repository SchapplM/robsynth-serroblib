% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:59:45
% EndTime: 2019-05-06 16:00:20
% DurationCPUTime: 17.79s
% Computational Cost: add. (308029->390), mult. (653922->472), div. (0->0), fcn. (417420->10), ass. (0->146)
t380 = -2 * qJD(3);
t341 = sin(qJ(1));
t345 = cos(qJ(1));
t321 = -t345 * g(1) - t341 * g(2);
t347 = qJD(1) ^ 2;
t292 = -t347 * pkin(1) + qJDD(1) * pkin(7) + t321;
t340 = sin(qJ(2));
t344 = cos(qJ(2));
t276 = -t344 * g(3) - t340 * t292;
t277 = -t340 * g(3) + t344 * t292;
t287 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t340 + Ifges(3,4) * t344) * qJD(1);
t308 = (mrSges(4,2) * t344 - mrSges(4,3) * t340) * qJD(1);
t371 = qJD(1) * qJD(2);
t368 = t344 * t371;
t310 = t340 * qJDD(1) + t368;
t369 = t340 * t371;
t311 = t344 * qJDD(1) - t369;
t372 = qJD(1) * t344;
t317 = -mrSges(4,1) * t372 - qJD(2) * mrSges(4,3);
t327 = t340 * qJD(1);
t307 = (-pkin(2) * t344 - qJ(3) * t340) * qJD(1);
t346 = qJD(2) ^ 2;
t251 = t346 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t380 - t307 * t372 - t277;
t319 = pkin(3) * t327 - qJD(2) * pkin(8);
t335 = t344 ^ 2;
t238 = -t335 * t347 * pkin(8) + t311 * pkin(3) + qJD(2) * t319 - t251;
t339 = sin(qJ(4));
t343 = cos(qJ(4));
t306 = t343 * qJD(2) - t339 * t372;
t267 = -t306 * qJD(4) - t339 * qJDD(2) - t343 * t311;
t305 = -t339 * qJD(2) - t343 * t372;
t268 = t305 * qJD(4) + t343 * qJDD(2) - t339 * t311;
t324 = t327 + qJD(4);
t273 = -t324 * mrSges(5,2) + t305 * mrSges(5,3);
t275 = t324 * mrSges(5,1) - t306 * mrSges(5,3);
t274 = t324 * pkin(4) - t306 * qJ(5);
t303 = t305 ^ 2;
t221 = -t267 * pkin(4) - t303 * qJ(5) + t306 * t274 + qJDD(5) + t238;
t336 = sin(pkin(10));
t337 = cos(pkin(10));
t243 = t337 * t267 - t336 * t268;
t244 = t336 * t267 + t337 * t268;
t270 = t337 * t305 - t336 * t306;
t254 = -t324 * mrSges(6,2) + t270 * mrSges(6,3);
t271 = t336 * t305 + t337 * t306;
t255 = t324 * mrSges(6,1) - t271 * mrSges(6,3);
t256 = t324 * pkin(5) - t271 * pkin(9);
t269 = t270 ^ 2;
t205 = -t243 * pkin(5) - t269 * pkin(9) + t271 * t256 + t221;
t338 = sin(qJ(6));
t342 = cos(qJ(6));
t247 = t338 * t270 + t342 * t271;
t214 = -t247 * qJD(6) + t342 * t243 - t338 * t244;
t246 = t342 * t270 - t338 * t271;
t215 = t246 * qJD(6) + t338 * t243 + t342 * t244;
t322 = qJD(6) + t324;
t230 = -t322 * mrSges(7,2) + t246 * mrSges(7,3);
t231 = t322 * mrSges(7,1) - t247 * mrSges(7,3);
t361 = m(7) * t205 - t214 * mrSges(7,1) + t215 * mrSges(7,2) - t246 * t230 + t247 * t231;
t354 = m(6) * t221 - t243 * mrSges(6,1) + t244 * mrSges(6,2) - t270 * t254 + t271 * t255 + t361;
t195 = -m(5) * t238 + t267 * mrSges(5,1) - t268 * mrSges(5,2) + t305 * t273 - t306 * t275 - t354;
t318 = mrSges(4,1) * t327 + qJD(2) * mrSges(4,2);
t351 = -m(4) * t251 + qJDD(2) * mrSges(4,3) + qJD(2) * t318 + t308 * t372 - t195;
t320 = t341 * g(1) - t345 * g(2);
t364 = -qJDD(1) * pkin(1) - t320;
t355 = pkin(2) * t369 + t327 * t380 + (-t310 - t368) * qJ(3) + t364;
t229 = -t319 * t327 + (-pkin(3) * t335 - pkin(7)) * t347 + (-pkin(2) - pkin(8)) * t311 + t355;
t253 = -qJDD(2) * pkin(2) - t346 * qJ(3) + t307 * t327 + qJDD(3) - t276;
t239 = (-t340 * t344 * t347 - qJDD(2)) * pkin(8) + (t310 - t368) * pkin(3) + t253;
t218 = -t339 * t229 + t343 * t239;
t304 = qJDD(4) + t310;
t208 = (t305 * t324 - t268) * qJ(5) + (t305 * t306 + t304) * pkin(4) + t218;
t219 = t343 * t229 + t339 * t239;
t210 = -t303 * pkin(4) + t267 * qJ(5) - t324 * t274 + t219;
t202 = -0.2e1 * qJD(5) * t271 + t337 * t208 - t336 * t210;
t199 = (t270 * t324 - t244) * pkin(9) + (t270 * t271 + t304) * pkin(5) + t202;
t203 = 0.2e1 * qJD(5) * t270 + t336 * t208 + t337 * t210;
t200 = -t269 * pkin(5) + t243 * pkin(9) - t324 * t256 + t203;
t198 = t338 * t199 + t342 * t200;
t222 = Ifges(7,5) * t247 + Ifges(7,6) * t246 + Ifges(7,3) * t322;
t224 = Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t322;
t295 = qJDD(6) + t304;
t185 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t215 + Ifges(7,2) * t214 + Ifges(7,6) * t295 - t247 * t222 + t322 * t224;
t197 = t342 * t199 - t338 * t200;
t223 = Ifges(7,4) * t247 + Ifges(7,2) * t246 + Ifges(7,6) * t322;
t186 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t215 + Ifges(7,4) * t214 + Ifges(7,5) * t295 + t246 * t222 - t322 * t223;
t240 = Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t324;
t242 = Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t324;
t226 = -t246 * mrSges(7,1) + t247 * mrSges(7,2);
t190 = m(7) * t197 + t295 * mrSges(7,1) - t215 * mrSges(7,3) - t247 * t226 + t322 * t230;
t191 = m(7) * t198 - t295 * mrSges(7,2) + t214 * mrSges(7,3) + t246 * t226 - t322 * t231;
t365 = -t338 * t190 + t342 * t191;
t170 = -mrSges(6,1) * t221 + mrSges(6,3) * t203 + Ifges(6,4) * t244 + Ifges(6,2) * t243 + Ifges(6,6) * t304 - pkin(5) * t361 + pkin(9) * t365 + t342 * t185 + t338 * t186 - t271 * t240 + t324 * t242;
t184 = t342 * t190 + t338 * t191;
t241 = Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t324;
t171 = mrSges(6,2) * t221 - mrSges(6,3) * t202 + Ifges(6,1) * t244 + Ifges(6,4) * t243 + Ifges(6,5) * t304 - pkin(9) * t184 - t338 * t185 + t342 * t186 + t270 * t240 - t324 * t241;
t257 = Ifges(5,5) * t306 + Ifges(5,6) * t305 + Ifges(5,3) * t324;
t259 = Ifges(5,1) * t306 + Ifges(5,4) * t305 + Ifges(5,5) * t324;
t248 = -t270 * mrSges(6,1) + t271 * mrSges(6,2);
t181 = m(6) * t202 + t304 * mrSges(6,1) - t244 * mrSges(6,3) - t271 * t248 + t324 * t254 + t184;
t182 = m(6) * t203 - t304 * mrSges(6,2) + t243 * mrSges(6,3) + t270 * t248 - t324 * t255 + t365;
t366 = -t336 * t181 + t337 * t182;
t157 = -mrSges(5,1) * t238 + mrSges(5,3) * t219 + Ifges(5,4) * t268 + Ifges(5,2) * t267 + Ifges(5,6) * t304 - pkin(4) * t354 + qJ(5) * t366 + t337 * t170 + t336 * t171 - t306 * t257 + t324 * t259;
t177 = t337 * t181 + t336 * t182;
t258 = Ifges(5,4) * t306 + Ifges(5,2) * t305 + Ifges(5,6) * t324;
t159 = mrSges(5,2) * t238 - mrSges(5,3) * t218 + Ifges(5,1) * t268 + Ifges(5,4) * t267 + Ifges(5,5) * t304 - qJ(5) * t177 - t336 * t170 + t337 * t171 + t305 * t257 - t324 * t258;
t272 = -t305 * mrSges(5,1) + t306 * mrSges(5,2);
t174 = m(5) * t218 + t304 * mrSges(5,1) - t268 * mrSges(5,3) - t306 * t272 + t324 * t273 + t177;
t175 = m(5) * t219 - t304 * mrSges(5,2) + t267 * mrSges(5,3) + t305 * t272 - t324 * t275 + t366;
t168 = t343 * t174 + t339 * t175;
t289 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t340 - Ifges(4,6) * t344) * qJD(1);
t357 = -mrSges(4,2) * t253 + mrSges(4,3) * t251 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t310 + Ifges(4,5) * t311 + pkin(8) * t168 + t339 * t157 - t343 * t159 - t289 * t372;
t360 = -m(4) * t253 - t310 * mrSges(4,1) - t168;
t288 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t340 - Ifges(4,3) * t344) * qJD(1);
t373 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t340 + Ifges(3,2) * t344) * qJD(1) - t288;
t379 = (-t344 * t287 + t340 * t373) * qJD(1) + mrSges(3,1) * t276 - mrSges(3,2) * t277 + Ifges(3,5) * t310 + Ifges(3,6) * t311 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t317 - t308 * t327 + t360) + qJ(3) * (t311 * mrSges(4,1) + t351) - t357;
t377 = t347 * pkin(7);
t376 = Ifges(3,4) + Ifges(4,6);
t169 = -t339 * t174 + t343 * t175;
t290 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t340 - Ifges(4,5) * t344) * qJD(1);
t374 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t340 + Ifges(3,6) * t344) * qJD(1) + t290;
t309 = (-mrSges(3,1) * t344 + mrSges(3,2) * t340) * qJD(1);
t316 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t372;
t165 = m(3) * t276 - t310 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t316 - t317) * qJD(2) + (-t308 - t309) * t327 + t360;
t315 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t327;
t193 = -qJDD(2) * mrSges(3,2) + t309 * t372 + (mrSges(3,3) + mrSges(4,1)) * t311 - qJD(2) * t315 + m(3) * t277 + t351;
t367 = -t340 * t165 + t344 * t193;
t249 = -t311 * pkin(2) + t355 - t377;
t363 = -m(4) * t249 - t311 * mrSges(4,2) + t318 * t327 - t169;
t166 = -t310 * mrSges(4,3) + t317 * t372 - t363;
t291 = t364 - t377;
t356 = -mrSges(4,1) * t251 + mrSges(4,2) * t249 - pkin(3) * t195 - pkin(8) * t169 - t343 * t157 - t339 * t159;
t154 = -mrSges(3,1) * t291 + mrSges(3,3) * t277 - pkin(2) * t166 + (Ifges(3,2) + Ifges(4,3)) * t311 + t376 * t310 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t287 - t289) * qJD(2) - t374 * t327 + t356;
t358 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t215 - Ifges(7,6) * t214 - Ifges(7,3) * t295 - t247 * t223 + t246 * t224;
t353 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t244 - Ifges(6,6) * t243 - Ifges(6,3) * t304 - pkin(5) * t184 - t271 * t241 + t270 * t242 + t358;
t349 = -mrSges(5,1) * t218 + mrSges(5,2) * t219 - Ifges(5,5) * t268 - Ifges(5,6) * t267 - Ifges(5,3) * t304 - pkin(4) * t177 - t306 * t258 + t305 * t259 + t353;
t348 = -mrSges(4,1) * t253 + mrSges(4,3) * t249 - pkin(3) * t168 + t349;
t156 = t376 * t311 + (Ifges(3,1) + Ifges(4,2)) * t310 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t373 * qJD(2) - t348 + mrSges(3,2) * t291 - mrSges(3,3) * t276 + t374 * t372 - qJ(3) * t166;
t352 = -m(3) * t291 + t316 * t372 + t311 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t310 + (-t315 * t340 - t317 * t344) * qJD(1) + t363;
t359 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t367 + t344 * t154 + t340 * t156;
t163 = m(2) * t320 + qJDD(1) * mrSges(2,1) - t347 * mrSges(2,2) + t352;
t162 = t344 * t165 + t340 * t193;
t160 = m(2) * t321 - t347 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t367;
t152 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t379;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t347 - pkin(7) * t162 - t154 * t340 + t156 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t151 - t341 * t152 - pkin(6) * (t160 * t341 + t163 * t345), t151, t156, -t288 * t327 - t357, t159, t171, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t341 * t151 + t345 * t152 + pkin(6) * (t160 * t345 - t163 * t341), t152, t154, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t310 - Ifges(4,6) * t311 - qJD(2) * t288 - t290 * t372 + t348, t157, t170, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, t379, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t310 - Ifges(4,3) * t311 + qJD(2) * t289 + t290 * t327 - t356, -t349, -t353, -t358;];
m_new  = t1;
