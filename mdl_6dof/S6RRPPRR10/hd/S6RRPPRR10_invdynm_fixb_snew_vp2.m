% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 11:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:46:44
% EndTime: 2019-05-06 11:47:16
% DurationCPUTime: 17.29s
% Computational Cost: add. (291467->390), mult. (641751->472), div. (0->0), fcn. (409665->10), ass. (0->146)
t380 = -2 * qJD(3);
t341 = sin(qJ(1));
t345 = cos(qJ(1));
t321 = -g(1) * t345 - g(2) * t341;
t347 = qJD(1) ^ 2;
t292 = -pkin(1) * t347 + qJDD(1) * pkin(7) + t321;
t340 = sin(qJ(2));
t344 = cos(qJ(2));
t266 = -t344 * g(3) - t340 * t292;
t267 = -g(3) * t340 + t344 * t292;
t287 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t340 + Ifges(3,4) * t344) * qJD(1);
t308 = (mrSges(4,2) * t344 - mrSges(4,3) * t340) * qJD(1);
t371 = qJD(1) * qJD(2);
t369 = t344 * t371;
t310 = qJDD(1) * t340 + t369;
t368 = t340 * t371;
t311 = qJDD(1) * t344 - t368;
t372 = qJD(1) * t344;
t318 = -mrSges(4,1) * t372 - qJD(2) * mrSges(4,3);
t327 = t340 * qJD(1);
t307 = (-pkin(2) * t344 - qJ(3) * t340) * qJD(1);
t346 = qJD(2) ^ 2;
t251 = pkin(2) * t346 - qJDD(2) * qJ(3) + qJD(2) * t380 - t307 * t372 - t267;
t317 = pkin(3) * t327 - qJD(2) * qJ(4);
t335 = t344 ^ 2;
t240 = -qJ(4) * t335 * t347 + pkin(3) * t311 + qJD(2) * t317 + qJDD(4) - t251;
t336 = sin(pkin(10));
t337 = cos(pkin(10));
t298 = -qJD(2) * t336 - t337 * t372;
t271 = -mrSges(5,2) * t327 + mrSges(5,3) * t298;
t299 = qJD(2) * t337 - t336 * t372;
t272 = mrSges(5,1) * t327 - mrSges(5,3) * t299;
t273 = -qJDD(2) * t336 - t311 * t337;
t274 = qJDD(2) * t337 - t311 * t336;
t275 = pkin(4) * t327 - pkin(8) * t299;
t297 = t298 ^ 2;
t225 = -pkin(4) * t273 - pkin(8) * t297 + t299 * t275 + t240;
t339 = sin(qJ(5));
t343 = cos(qJ(5));
t264 = t298 * t339 + t299 * t343;
t235 = -qJD(5) * t264 + t273 * t343 - t274 * t339;
t263 = t298 * t343 - t299 * t339;
t236 = qJD(5) * t263 + t273 * t339 + t274 * t343;
t324 = t327 + qJD(5);
t254 = -mrSges(6,2) * t324 + mrSges(6,3) * t263;
t255 = mrSges(6,1) * t324 - mrSges(6,3) * t264;
t256 = pkin(5) * t324 - pkin(9) * t264;
t262 = t263 ^ 2;
t205 = -pkin(5) * t235 - pkin(9) * t262 + t256 * t264 + t225;
t338 = sin(qJ(6));
t342 = cos(qJ(6));
t247 = t263 * t338 + t264 * t342;
t214 = -qJD(6) * t247 + t235 * t342 - t236 * t338;
t246 = t263 * t342 - t264 * t338;
t215 = qJD(6) * t246 + t235 * t338 + t236 * t342;
t322 = qJD(6) + t324;
t230 = -mrSges(7,2) * t322 + mrSges(7,3) * t246;
t231 = mrSges(7,1) * t322 - mrSges(7,3) * t247;
t361 = m(7) * t205 - t214 * mrSges(7,1) + t215 * mrSges(7,2) - t246 * t230 + t247 * t231;
t354 = m(6) * t225 - t235 * mrSges(6,1) + t236 * mrSges(6,2) - t263 * t254 + t264 * t255 + t361;
t195 = -m(5) * t240 + t273 * mrSges(5,1) - t274 * mrSges(5,2) + t298 * t271 - t299 * t272 - t354;
t319 = mrSges(4,1) * t327 + qJD(2) * mrSges(4,2);
t351 = -m(4) * t251 + qJDD(2) * mrSges(4,3) + qJD(2) * t319 + t308 * t372 - t195;
t320 = t341 * g(1) - t345 * g(2);
t364 = -qJDD(1) * pkin(1) - t320;
t356 = pkin(2) * t368 + t327 * t380 + (-t310 - t369) * qJ(3) + t364;
t229 = -t317 * t327 + (-pkin(3) * t335 - pkin(7)) * t347 + (-pkin(2) - qJ(4)) * t311 + t356;
t253 = -qJDD(2) * pkin(2) - t346 * qJ(3) + t307 * t327 + qJDD(3) - t266;
t244 = (-t340 * t344 * t347 - qJDD(2)) * qJ(4) + (t310 - t369) * pkin(3) + t253;
t218 = -0.2e1 * qJD(4) * t299 - t336 * t229 + t337 * t244;
t208 = (t298 * t327 - t274) * pkin(8) + (t298 * t299 + t310) * pkin(4) + t218;
t219 = 0.2e1 * qJD(4) * t298 + t337 * t229 + t336 * t244;
t210 = -pkin(4) * t297 + pkin(8) * t273 - t275 * t327 + t219;
t202 = t343 * t208 - t339 * t210;
t306 = qJDD(5) + t310;
t199 = (t263 * t324 - t236) * pkin(9) + (t263 * t264 + t306) * pkin(5) + t202;
t203 = t339 * t208 + t343 * t210;
t200 = -pkin(5) * t262 + pkin(9) * t235 - t256 * t324 + t203;
t198 = t199 * t338 + t200 * t342;
t220 = Ifges(7,5) * t247 + Ifges(7,6) * t246 + Ifges(7,3) * t322;
t222 = Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t322;
t294 = qJDD(6) + t306;
t185 = -mrSges(7,1) * t205 + mrSges(7,3) * t198 + Ifges(7,4) * t215 + Ifges(7,2) * t214 + Ifges(7,6) * t294 - t220 * t247 + t222 * t322;
t197 = t199 * t342 - t200 * t338;
t221 = Ifges(7,4) * t247 + Ifges(7,2) * t246 + Ifges(7,6) * t322;
t186 = mrSges(7,2) * t205 - mrSges(7,3) * t197 + Ifges(7,1) * t215 + Ifges(7,4) * t214 + Ifges(7,5) * t294 + t220 * t246 - t221 * t322;
t241 = Ifges(6,5) * t264 + Ifges(6,6) * t263 + Ifges(6,3) * t324;
t243 = Ifges(6,1) * t264 + Ifges(6,4) * t263 + Ifges(6,5) * t324;
t226 = -mrSges(7,1) * t246 + mrSges(7,2) * t247;
t190 = m(7) * t197 + mrSges(7,1) * t294 - mrSges(7,3) * t215 - t226 * t247 + t230 * t322;
t191 = m(7) * t198 - mrSges(7,2) * t294 + mrSges(7,3) * t214 + t226 * t246 - t231 * t322;
t365 = -t190 * t338 + t342 * t191;
t170 = -mrSges(6,1) * t225 + mrSges(6,3) * t203 + Ifges(6,4) * t236 + Ifges(6,2) * t235 + Ifges(6,6) * t306 - pkin(5) * t361 + pkin(9) * t365 + t342 * t185 + t338 * t186 - t264 * t241 + t324 * t243;
t184 = t342 * t190 + t338 * t191;
t242 = Ifges(6,4) * t264 + Ifges(6,2) * t263 + Ifges(6,6) * t324;
t171 = mrSges(6,2) * t225 - mrSges(6,3) * t202 + Ifges(6,1) * t236 + Ifges(6,4) * t235 + Ifges(6,5) * t306 - pkin(9) * t184 - t185 * t338 + t186 * t342 + t241 * t263 - t242 * t324;
t257 = Ifges(5,5) * t299 + Ifges(5,6) * t298 + Ifges(5,3) * t327;
t259 = Ifges(5,1) * t299 + Ifges(5,4) * t298 + Ifges(5,5) * t327;
t248 = -mrSges(6,1) * t263 + mrSges(6,2) * t264;
t181 = m(6) * t202 + mrSges(6,1) * t306 - mrSges(6,3) * t236 - t248 * t264 + t254 * t324 + t184;
t182 = m(6) * t203 - mrSges(6,2) * t306 + mrSges(6,3) * t235 + t248 * t263 - t255 * t324 + t365;
t366 = -t181 * t339 + t343 * t182;
t157 = -mrSges(5,1) * t240 + mrSges(5,3) * t219 + Ifges(5,4) * t274 + Ifges(5,2) * t273 + Ifges(5,6) * t310 - pkin(4) * t354 + pkin(8) * t366 + t343 * t170 + t339 * t171 - t299 * t257 + t259 * t327;
t177 = t343 * t181 + t339 * t182;
t258 = Ifges(5,4) * t299 + Ifges(5,2) * t298 + Ifges(5,6) * t327;
t159 = mrSges(5,2) * t240 - mrSges(5,3) * t218 + Ifges(5,1) * t274 + Ifges(5,4) * t273 + Ifges(5,5) * t310 - pkin(8) * t177 - t170 * t339 + t171 * t343 + t257 * t298 - t258 * t327;
t265 = -mrSges(5,1) * t298 + mrSges(5,2) * t299;
t174 = m(5) * t218 + mrSges(5,1) * t310 - mrSges(5,3) * t274 - t265 * t299 + t271 * t327 + t177;
t175 = m(5) * t219 - mrSges(5,2) * t310 + mrSges(5,3) * t273 + t265 * t298 - t272 * t327 + t366;
t168 = t174 * t337 + t175 * t336;
t289 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t340 - Ifges(4,6) * t344) * qJD(1);
t357 = -mrSges(4,2) * t253 + mrSges(4,3) * t251 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t310 + Ifges(4,5) * t311 + qJ(4) * t168 + t336 * t157 - t337 * t159 - t289 * t372;
t360 = -m(4) * t253 - t310 * mrSges(4,1) - t168;
t288 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t340 - Ifges(4,3) * t344) * qJD(1);
t373 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t340 + Ifges(3,2) * t344) * qJD(1) - t288;
t379 = (-t344 * t287 + t340 * t373) * qJD(1) + mrSges(3,1) * t266 - mrSges(3,2) * t267 + Ifges(3,5) * t310 + Ifges(3,6) * t311 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t318 - t308 * t327 + t360) + qJ(3) * (t311 * mrSges(4,1) + t351) - t357;
t377 = t347 * pkin(7);
t376 = Ifges(3,4) + Ifges(4,6);
t169 = -t336 * t174 + t337 * t175;
t290 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t340 - Ifges(4,5) * t344) * qJD(1);
t374 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t340 + Ifges(3,6) * t344) * qJD(1) + t290;
t309 = (-mrSges(3,1) * t344 + mrSges(3,2) * t340) * qJD(1);
t316 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t372;
t165 = m(3) * t266 - mrSges(3,3) * t310 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t316 - t318) * qJD(2) + (-t308 - t309) * t327 + t360;
t315 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t327;
t193 = -qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t311 + t351 - qJD(2) * t315 + m(3) * t267 + t309 * t372;
t367 = -t165 * t340 + t344 * t193;
t249 = -t311 * pkin(2) + t356 - t377;
t363 = -m(4) * t249 - t311 * mrSges(4,2) + t319 * t327 - t169;
t166 = -t310 * mrSges(4,3) + t318 * t372 - t363;
t291 = t364 - t377;
t355 = -mrSges(4,1) * t251 + mrSges(4,2) * t249 - pkin(3) * t195 - qJ(4) * t169 - t337 * t157 - t336 * t159;
t154 = -mrSges(3,1) * t291 + mrSges(3,3) * t267 - pkin(2) * t166 + (Ifges(3,2) + Ifges(4,3)) * t311 + t376 * t310 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t287 - t289) * qJD(2) - t374 * t327 + t355;
t358 = -mrSges(7,1) * t197 + mrSges(7,2) * t198 - Ifges(7,5) * t215 - Ifges(7,6) * t214 - Ifges(7,3) * t294 - t247 * t221 + t246 * t222;
t353 = -mrSges(6,1) * t202 + mrSges(6,2) * t203 - Ifges(6,5) * t236 - Ifges(6,6) * t235 - Ifges(6,3) * t306 - pkin(5) * t184 - t264 * t242 + t263 * t243 + t358;
t349 = -mrSges(5,1) * t218 + mrSges(5,2) * t219 - Ifges(5,5) * t274 - Ifges(5,6) * t273 - Ifges(5,3) * t310 - pkin(4) * t177 - t299 * t258 + t298 * t259 + t353;
t348 = -mrSges(4,1) * t253 + mrSges(4,3) * t249 - pkin(3) * t168 + t349;
t156 = -qJ(3) * t166 + t376 * t311 + (Ifges(3,1) + Ifges(4,2)) * t310 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t348 - t373 * qJD(2) + mrSges(3,2) * t291 - mrSges(3,3) * t266 + t374 * t372;
t352 = -m(3) * t291 + t316 * t372 + t311 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t310 + (-t315 * t340 - t318 * t344) * qJD(1) + t363;
t359 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t367 + t344 * t154 + t340 * t156;
t163 = m(2) * t320 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t347 + t352;
t162 = t165 * t344 + t193 * t340;
t160 = m(2) * t321 - mrSges(2,1) * t347 - qJDD(1) * mrSges(2,2) + t367;
t152 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t347 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t162 - t379;
t151 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t347 - pkin(7) * t162 - t154 * t340 + t156 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t151 - t341 * t152 - pkin(6) * (t160 * t341 + t163 * t345), t151, t156, -t288 * t327 - t357, t159, t171, t186; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t341 * t151 + t345 * t152 + pkin(6) * (t160 * t345 - t163 * t341), t152, t154, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t310 - Ifges(4,6) * t311 - qJD(2) * t288 - t290 * t372 + t348, t157, t170, t185; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t359, t359, t379, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t310 - Ifges(4,3) * t311 + qJD(2) * t289 + t290 * t327 - t355, -t349, -t353, -t358;];
m_new  = t1;
