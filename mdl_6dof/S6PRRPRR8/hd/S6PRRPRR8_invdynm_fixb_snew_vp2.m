% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:15:35
% EndTime: 2019-05-05 06:16:11
% DurationCPUTime: 22.51s
% Computational Cost: add. (404043->361), mult. (867798->459), div. (0->0), fcn. (641556->14), ass. (0->158)
t379 = -2 * qJD(4);
t324 = sin(pkin(12));
t327 = cos(pkin(12));
t309 = g(1) * t324 - g(2) * t327;
t323 = -g(3) + qJDD(1);
t326 = sin(pkin(6));
t329 = cos(pkin(6));
t284 = -t309 * t326 + t323 * t329;
t325 = sin(pkin(7));
t336 = cos(qJ(3));
t310 = -g(1) * t327 - g(2) * t324;
t333 = sin(qJ(2));
t337 = cos(qJ(2));
t364 = t329 * t337;
t367 = t326 * t337;
t260 = t309 * t364 - t310 * t333 + t323 * t367;
t338 = qJD(2) ^ 2;
t256 = pkin(9) * t325 * t338 + qJDD(2) * pkin(2) + t260;
t328 = cos(pkin(7));
t371 = t256 * t328;
t378 = t336 * (t284 * t325 + t371);
t320 = qJD(2) * t328 + qJD(3);
t332 = sin(qJ(3));
t361 = qJD(2) * t325;
t356 = t332 * t361;
t377 = (pkin(3) * t320 + t379) * t356;
t365 = t329 * t333;
t368 = t326 * t333;
t261 = t309 * t365 + t337 * t310 + t323 * t368;
t358 = qJDD(2) * t325;
t257 = -pkin(2) * t338 + pkin(9) * t358 + t261;
t366 = t328 * t332;
t369 = t325 * t332;
t230 = t256 * t366 + t336 * t257 + t284 * t369;
t295 = (-pkin(3) * t336 - qJ(4) * t332) * t361;
t318 = t320 ^ 2;
t319 = qJDD(2) * t328 + qJDD(3);
t360 = qJD(2) * t336;
t355 = t325 * t360;
t224 = pkin(3) * t318 - t319 * qJ(4) - t295 * t355 + t320 * t379 - t230;
t376 = -pkin(3) - pkin(10);
t280 = t328 * t284;
t240 = -t256 * t325 + t280;
t291 = mrSges(4,1) * t320 - mrSges(4,3) * t356;
t292 = -mrSges(4,2) * t320 + mrSges(4,3) * t355;
t294 = mrSges(5,1) * t356 + mrSges(5,2) * t320;
t299 = (qJD(3) * t360 + qJDD(2) * t332) * t325;
t300 = -qJD(3) * t356 + t336 * t358;
t254 = t332 * t257;
t351 = -qJ(4) * t318 + t295 * t356 + qJDD(4) + t254;
t370 = t325 ^ 2 * t338;
t221 = pkin(4) * t299 + t376 * t319 + (-pkin(10) * t332 * t370 - t371 + (-pkin(4) * qJD(2) * t320 - t284) * t325) * t336 + t351;
t298 = pkin(4) * t356 - pkin(10) * t320;
t357 = t336 ^ 2 * t370;
t223 = -pkin(4) * t357 - qJ(4) * t299 + t280 + t376 * t300 + (-t256 + (-qJ(4) * t320 * t336 - t298 * t332) * qJD(2)) * t325 + t377;
t331 = sin(qJ(5));
t335 = cos(qJ(5));
t216 = t331 * t221 + t335 * t223;
t283 = t320 * t335 - t331 * t355;
t251 = -qJD(5) * t283 - t300 * t335 - t319 * t331;
t282 = -t320 * t331 - t335 * t355;
t258 = -mrSges(6,1) * t282 + mrSges(6,2) * t283;
t308 = qJD(5) + t356;
t265 = mrSges(6,1) * t308 - mrSges(6,3) * t283;
t290 = qJDD(5) + t299;
t259 = -pkin(5) * t282 - pkin(11) * t283;
t306 = t308 ^ 2;
t213 = -pkin(5) * t306 + pkin(11) * t290 + t259 * t282 + t216;
t219 = pkin(4) * t300 - pkin(10) * t357 + t320 * t298 - t224;
t252 = qJD(5) * t282 - t300 * t331 + t319 * t335;
t217 = (-t282 * t308 - t252) * pkin(11) + (t283 * t308 - t251) * pkin(5) + t219;
t330 = sin(qJ(6));
t334 = cos(qJ(6));
t210 = -t213 * t330 + t217 * t334;
t262 = -t283 * t330 + t308 * t334;
t233 = qJD(6) * t262 + t252 * t334 + t290 * t330;
t263 = t283 * t334 + t308 * t330;
t238 = -mrSges(7,1) * t262 + mrSges(7,2) * t263;
t281 = qJD(6) - t282;
t241 = -mrSges(7,2) * t281 + mrSges(7,3) * t262;
t249 = qJDD(6) - t251;
t206 = m(7) * t210 + mrSges(7,1) * t249 - mrSges(7,3) * t233 - t238 * t263 + t241 * t281;
t211 = t213 * t334 + t217 * t330;
t232 = -qJD(6) * t263 - t252 * t330 + t290 * t334;
t242 = mrSges(7,1) * t281 - mrSges(7,3) * t263;
t207 = m(7) * t211 - mrSges(7,2) * t249 + mrSges(7,3) * t232 + t238 * t262 - t242 * t281;
t353 = -t206 * t330 + t334 * t207;
t193 = m(6) * t216 - mrSges(6,2) * t290 + mrSges(6,3) * t251 + t258 * t282 - t265 * t308 + t353;
t215 = t221 * t335 - t223 * t331;
t264 = -mrSges(6,2) * t308 + mrSges(6,3) * t282;
t212 = -pkin(5) * t290 - pkin(11) * t306 + t259 * t283 - t215;
t345 = -m(7) * t212 + t232 * mrSges(7,1) - mrSges(7,2) * t233 + t262 * t241 - t242 * t263;
t202 = m(6) * t215 + mrSges(6,1) * t290 - mrSges(6,3) * t252 - t258 * t283 + t264 * t308 + t345;
t187 = t335 * t193 - t202 * t331;
t228 = -pkin(3) * t300 + (-t320 * t355 - t299) * qJ(4) + t240 + t377;
t293 = -mrSges(5,1) * t355 - mrSges(5,3) * t320;
t350 = m(5) * t228 - t299 * mrSges(5,3) + t293 * t355 + t187;
t374 = mrSges(4,1) - mrSges(5,2);
t184 = m(4) * t240 + mrSges(4,2) * t299 - t374 * t300 + (-t292 * t336 + (t291 - t294) * t332) * t361 + t350;
t297 = (-mrSges(4,1) * t336 + mrSges(4,2) * t332) * t361;
t196 = t334 * t206 + t330 * t207;
t194 = -m(6) * t219 + mrSges(6,1) * t251 - t252 * mrSges(6,2) + t264 * t282 - t283 * t265 - t196;
t296 = (mrSges(5,2) * t336 - mrSges(5,3) * t332) * t361;
t341 = -m(5) * t224 + t319 * mrSges(5,3) + t320 * t294 + t296 * t355 - t194;
t191 = (mrSges(4,3) + mrSges(5,1)) * t300 + m(4) * t230 - mrSges(4,2) * t319 - t291 * t320 + t297 * t355 + t341;
t229 = -t254 + t378;
t186 = t193 * t331 + t202 * t335;
t226 = -pkin(3) * t319 + t351 - t378;
t347 = -m(5) * t226 - t299 * mrSges(5,1) - t186;
t182 = m(4) * t229 - mrSges(4,3) * t299 + (t292 - t293) * t320 + t374 * t319 + (-t296 - t297) * t356 + t347;
t372 = t182 * t336;
t172 = -t184 * t325 + t191 * t366 + t328 * t372;
t169 = m(3) * t260 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t338 + t172;
t176 = -t182 * t332 + t336 * t191;
t175 = m(3) * t261 - mrSges(3,1) * t338 - qJDD(2) * mrSges(3,2) + t176;
t166 = -t169 * t333 + t337 * t175;
t375 = pkin(8) * t166;
t373 = Ifges(4,4) + Ifges(5,6);
t273 = Ifges(5,1) * t320 + (-Ifges(5,4) * t332 - Ifges(5,5) * t336) * t361;
t363 = Ifges(4,3) * t320 + (Ifges(4,5) * t332 + Ifges(4,6) * t336) * t361 + t273;
t271 = Ifges(5,5) * t320 + (-Ifges(5,6) * t332 - Ifges(5,3) * t336) * t361;
t362 = -Ifges(4,6) * t320 - (Ifges(4,4) * t332 + Ifges(4,2) * t336) * t361 + t271;
t171 = t328 * t184 + t191 * t369 + t325 * t372;
t354 = t273 * t361;
t170 = m(3) * t284 + t171;
t161 = t169 * t364 - t170 * t326 + t175 * t365;
t185 = mrSges(5,2) * t300 - t294 * t356 + t350;
t270 = Ifges(4,5) * t320 + (Ifges(4,1) * t332 + Ifges(4,4) * t336) * t361;
t272 = Ifges(5,4) * t320 + (-Ifges(5,2) * t332 - Ifges(5,6) * t336) * t361;
t234 = Ifges(7,5) * t263 + Ifges(7,6) * t262 + Ifges(7,3) * t281;
t236 = Ifges(7,1) * t263 + Ifges(7,4) * t262 + Ifges(7,5) * t281;
t200 = -mrSges(7,1) * t212 + mrSges(7,3) * t211 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t249 - t234 * t263 + t236 * t281;
t235 = Ifges(7,4) * t263 + Ifges(7,2) * t262 + Ifges(7,6) * t281;
t201 = mrSges(7,2) * t212 - mrSges(7,3) * t210 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t249 + t234 * t262 - t235 * t281;
t243 = Ifges(6,5) * t283 + Ifges(6,6) * t282 + Ifges(6,3) * t308;
t244 = Ifges(6,4) * t283 + Ifges(6,2) * t282 + Ifges(6,6) * t308;
t178 = mrSges(6,2) * t219 - mrSges(6,3) * t215 + Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t290 - pkin(11) * t196 - t200 * t330 + t201 * t334 + t243 * t282 - t244 * t308;
t245 = Ifges(6,1) * t283 + Ifges(6,4) * t282 + Ifges(6,5) * t308;
t340 = mrSges(7,1) * t210 - mrSges(7,2) * t211 + Ifges(7,5) * t233 + Ifges(7,6) * t232 + Ifges(7,3) * t249 + t235 * t263 - t236 * t262;
t179 = -mrSges(6,1) * t219 + mrSges(6,3) * t216 + Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t290 - pkin(5) * t196 - t243 * t283 + t245 * t308 - t340;
t342 = -mrSges(5,1) * t224 + mrSges(5,2) * t228 - pkin(4) * t194 - pkin(10) * t187 - t331 * t178 - t335 * t179;
t163 = -mrSges(4,1) * t240 + mrSges(4,3) * t230 - pkin(3) * t185 + (t270 - t272) * t320 + (Ifges(4,6) - Ifges(5,5)) * t319 + (Ifges(4,2) + Ifges(5,3)) * t300 + t373 * t299 - t363 * t356 + t342;
t343 = -mrSges(6,1) * t215 + mrSges(6,2) * t216 - Ifges(6,5) * t252 - Ifges(6,6) * t251 - Ifges(6,3) * t290 - pkin(5) * t345 - pkin(11) * t353 - t334 * t200 - t330 * t201 - t283 * t244 + t282 * t245;
t339 = -mrSges(5,1) * t226 + mrSges(5,3) * t228 - pkin(4) * t186 + t343;
t167 = t362 * t320 + (Ifges(4,5) - Ifges(5,4)) * t319 + t373 * t300 + (Ifges(4,1) + Ifges(5,2)) * t299 + t363 * t355 + mrSges(4,2) * t240 - mrSges(4,3) * t229 - qJ(4) * t185 - t339;
t348 = pkin(9) * t176 + t163 * t336 + t167 * t332;
t344 = mrSges(5,2) * t226 - mrSges(5,3) * t224 + Ifges(5,1) * t319 - Ifges(5,4) * t299 - Ifges(5,5) * t300 - pkin(10) * t186 + t335 * t178 - t331 * t179 + t272 * t355;
t162 = pkin(3) * (-mrSges(5,2) * t319 - t293 * t320 + t347) + t344 + qJ(4) * (mrSges(5,1) * t300 + t341) + Ifges(4,3) * t319 + Ifges(4,5) * t299 + Ifges(4,6) * t300 + mrSges(4,1) * t229 - mrSges(4,2) * t230 + (-t270 * t336 + (-pkin(3) * t296 - t362) * t332) * t361;
t153 = mrSges(3,1) * t260 - mrSges(3,2) * t261 + Ifges(3,3) * qJDD(2) + pkin(2) * t172 + t162 * t328 + t325 * t348;
t155 = -mrSges(3,1) * t284 + mrSges(3,3) * t261 + Ifges(3,5) * t338 + Ifges(3,6) * qJDD(2) - pkin(2) * t171 - t162 * t325 + t328 * t348;
t157 = mrSges(3,2) * t284 - mrSges(3,3) * t260 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t338 - t163 * t332 + t167 * t336 + (-t171 * t325 - t172 * t328) * pkin(9);
t346 = mrSges(2,1) * t309 - mrSges(2,2) * t310 + pkin(1) * t161 + t329 * t153 + t155 * t367 + t157 * t368 + t326 * t375;
t164 = m(2) * t310 + t166;
t160 = t170 * t329 + (t169 * t337 + t175 * t333) * t326;
t158 = m(2) * t309 + t161;
t151 = mrSges(2,2) * t323 - mrSges(2,3) * t309 - t155 * t333 + t157 * t337 + (-t160 * t326 - t161 * t329) * pkin(8);
t150 = -mrSges(2,1) * t323 + mrSges(2,3) * t310 - pkin(1) * t160 - t153 * t326 + (t155 * t337 + t157 * t333 + t375) * t329;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t327 * t151 - t324 * t150 - qJ(1) * (t158 * t327 + t164 * t324), t151, t157, t167, -t271 * t356 + t344, t178, t201; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t324 * t151 + t327 * t150 + qJ(1) * (-t158 * t324 + t164 * t327), t150, t155, t163, Ifges(5,4) * t319 - Ifges(5,2) * t299 - Ifges(5,6) * t300 - t320 * t271 - t336 * t354 + t339, t179, t200; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t346, t346, t153, t162, Ifges(5,5) * t319 - Ifges(5,6) * t299 - Ifges(5,3) * t300 + t320 * t272 + t332 * t354 - t342, -t343, t340;];
m_new  = t1;
