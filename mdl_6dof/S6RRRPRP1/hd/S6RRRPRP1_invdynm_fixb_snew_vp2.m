% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:24:38
% EndTime: 2019-05-07 07:25:15
% DurationCPUTime: 20.98s
% Computational Cost: add. (359618->381), mult. (808035->470), div. (0->0), fcn. (592743->10), ass. (0->143)
t334 = sin(qJ(2));
t338 = cos(qJ(2));
t363 = qJD(1) * qJD(2);
t314 = t334 * qJDD(1) + t338 * t363;
t335 = sin(qJ(1));
t339 = cos(qJ(1));
t321 = -t339 * g(1) - t335 * g(2);
t340 = qJD(1) ^ 2;
t309 = -t340 * pkin(1) + qJDD(1) * pkin(7) + t321;
t369 = t334 * t309;
t371 = pkin(2) * t340;
t272 = qJDD(2) * pkin(2) - t314 * pkin(8) - t369 + (pkin(8) * t363 + t334 * t371 - g(3)) * t338;
t297 = -t334 * g(3) + t338 * t309;
t315 = t338 * qJDD(1) - t334 * t363;
t366 = qJD(1) * t334;
t319 = qJD(2) * pkin(2) - pkin(8) * t366;
t329 = t338 ^ 2;
t273 = t315 * pkin(8) - qJD(2) * t319 - t329 * t371 + t297;
t333 = sin(qJ(3));
t337 = cos(qJ(3));
t241 = t337 * t272 - t333 * t273;
t306 = (-t333 * t334 + t337 * t338) * qJD(1);
t281 = t306 * qJD(3) + t337 * t314 + t333 * t315;
t307 = (t333 * t338 + t334 * t337) * qJD(1);
t326 = qJDD(2) + qJDD(3);
t327 = qJD(2) + qJD(3);
t208 = (t306 * t327 - t281) * qJ(4) + (t306 * t307 + t326) * pkin(3) + t241;
t242 = t333 * t272 + t337 * t273;
t280 = -t307 * qJD(3) - t333 * t314 + t337 * t315;
t299 = t327 * pkin(3) - t307 * qJ(4);
t302 = t306 ^ 2;
t211 = -t302 * pkin(3) + t280 * qJ(4) - t327 * t299 + t242;
t330 = sin(pkin(10));
t331 = cos(pkin(10));
t294 = t330 * t306 + t331 * t307;
t202 = -0.2e1 * qJD(4) * t294 + t331 * t208 - t330 * t211;
t253 = t330 * t280 + t331 * t281;
t332 = sin(qJ(5));
t336 = cos(qJ(5));
t278 = -t332 * t294 + t336 * t327;
t225 = t278 * qJD(5) + t336 * t253 + t332 * t326;
t279 = t336 * t294 + t332 * t327;
t254 = -t278 * mrSges(7,1) + t279 * mrSges(7,2);
t293 = t331 * t306 - t330 * t307;
t203 = 0.2e1 * qJD(4) * t293 + t330 * t208 + t331 * t211;
t267 = -t293 * pkin(4) - t294 * pkin(9);
t325 = t327 ^ 2;
t200 = -t325 * pkin(4) + t326 * pkin(9) + t293 * t267 + t203;
t320 = t335 * g(1) - t339 * g(2);
t353 = -qJDD(1) * pkin(1) - t320;
t282 = -t315 * pkin(2) + t319 * t366 + (-pkin(8) * t329 - pkin(7)) * t340 + t353;
t227 = -t280 * pkin(3) - t302 * qJ(4) + t307 * t299 + qJDD(4) + t282;
t252 = t331 * t280 - t330 * t281;
t206 = (-t293 * t327 - t253) * pkin(9) + (t294 * t327 - t252) * pkin(4) + t227;
t194 = -t332 * t200 + t336 * t206;
t251 = qJDD(5) - t252;
t287 = qJD(5) - t293;
t190 = -0.2e1 * qJD(6) * t279 + (t278 * t287 - t225) * qJ(6) + (t278 * t279 + t251) * pkin(5) + t194;
t256 = -t287 * mrSges(7,2) + t278 * mrSges(7,3);
t362 = m(7) * t190 + t251 * mrSges(7,1) + t287 * t256;
t187 = -t225 * mrSges(7,3) - t279 * t254 + t362;
t195 = t336 * t200 + t332 * t206;
t224 = -t279 * qJD(5) - t332 * t253 + t336 * t326;
t232 = Ifges(6,4) * t279 + Ifges(6,2) * t278 + Ifges(6,6) * t287;
t233 = Ifges(7,1) * t279 + Ifges(7,4) * t278 + Ifges(7,5) * t287;
t234 = Ifges(6,1) * t279 + Ifges(6,4) * t278 + Ifges(6,5) * t287;
t258 = t287 * pkin(5) - t279 * qJ(6);
t274 = t278 ^ 2;
t193 = -t274 * pkin(5) + t224 * qJ(6) + 0.2e1 * qJD(6) * t278 - t287 * t258 + t195;
t231 = Ifges(7,4) * t279 + Ifges(7,2) * t278 + Ifges(7,6) * t287;
t350 = -mrSges(7,1) * t190 + mrSges(7,2) * t193 - Ifges(7,5) * t225 - Ifges(7,6) * t224 - Ifges(7,3) * t251 - t279 * t231;
t373 = mrSges(6,1) * t194 - mrSges(6,2) * t195 + Ifges(6,5) * t225 + Ifges(6,6) * t224 + Ifges(6,3) * t251 + pkin(5) * t187 + t279 * t232 - (t234 + t233) * t278 - t350;
t266 = -t293 * mrSges(5,1) + t294 * mrSges(5,2);
t284 = t327 * mrSges(5,1) - t294 * mrSges(5,3);
t255 = -t278 * mrSges(6,1) + t279 * mrSges(6,2);
t257 = -t287 * mrSges(6,2) + t278 * mrSges(6,3);
t179 = m(6) * t194 + t251 * mrSges(6,1) + t287 * t257 + (-t254 - t255) * t279 + (-mrSges(6,3) - mrSges(7,3)) * t225 + t362;
t361 = m(7) * t193 + t224 * mrSges(7,3) + t278 * t254;
t259 = t287 * mrSges(7,1) - t279 * mrSges(7,3);
t367 = -t287 * mrSges(6,1) + t279 * mrSges(6,3) - t259;
t370 = -mrSges(6,2) - mrSges(7,2);
t182 = m(6) * t195 + t224 * mrSges(6,3) + t251 * t370 + t278 * t255 + t287 * t367 + t361;
t357 = -t332 * t179 + t336 * t182;
t172 = m(5) * t203 - t326 * mrSges(5,2) + t252 * mrSges(5,3) + t293 * t266 - t327 * t284 + t357;
t283 = -t327 * mrSges(5,2) + t293 * mrSges(5,3);
t199 = -t326 * pkin(4) - t325 * pkin(9) + t294 * t267 - t202;
t197 = -t224 * pkin(5) - t274 * qJ(6) + t279 * t258 + qJDD(6) + t199;
t355 = -m(7) * t197 + t224 * mrSges(7,1) + t278 * t256;
t345 = -m(6) * t199 + t224 * mrSges(6,1) + t225 * t370 + t278 * t257 + t279 * t367 + t355;
t184 = m(5) * t202 + t326 * mrSges(5,1) - t253 * mrSges(5,3) - t294 * t266 + t327 * t283 + t345;
t165 = t330 * t172 + t331 * t184;
t295 = -t306 * mrSges(4,1) + t307 * mrSges(4,2);
t298 = -t327 * mrSges(4,2) + t306 * mrSges(4,3);
t162 = m(4) * t241 + t326 * mrSges(4,1) - t281 * mrSges(4,3) - t307 * t295 + t327 * t298 + t165;
t300 = t327 * mrSges(4,1) - t307 * mrSges(4,3);
t358 = t331 * t172 - t330 * t184;
t163 = m(4) * t242 - t326 * mrSges(4,2) + t280 * mrSges(4,3) + t306 * t295 - t327 * t300 + t358;
t157 = t337 * t162 + t333 * t163;
t296 = -t338 * g(3) - t369;
t304 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t334 + Ifges(3,2) * t338) * qJD(1);
t305 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t334 + Ifges(3,4) * t338) * qJD(1);
t289 = Ifges(4,4) * t307 + Ifges(4,2) * t306 + Ifges(4,6) * t327;
t290 = Ifges(4,1) * t307 + Ifges(4,4) * t306 + Ifges(4,5) * t327;
t229 = Ifges(7,5) * t279 + Ifges(7,6) * t278 + Ifges(7,3) * t287;
t230 = Ifges(6,5) * t279 + Ifges(6,6) * t278 + Ifges(6,3) * t287;
t351 = -mrSges(7,1) * t197 + mrSges(7,3) * t193 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t251 + t287 * t233;
t167 = Ifges(6,4) * t225 + Ifges(6,2) * t224 + Ifges(6,6) * t251 + t287 * t234 - mrSges(6,1) * t199 + mrSges(6,3) * t195 - pkin(5) * (t225 * mrSges(7,2) - t355) + qJ(6) * (-t251 * mrSges(7,2) - t287 * t259 + t361) + (-pkin(5) * t259 - t229 - t230) * t279 + t351;
t349 = mrSges(7,2) * t197 - mrSges(7,3) * t190 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t251 + t278 * t229;
t174 = mrSges(6,2) * t199 - mrSges(6,3) * t194 + Ifges(6,1) * t225 + Ifges(6,4) * t224 + Ifges(6,5) * t251 - qJ(6) * t187 + t278 * t230 + (-t231 - t232) * t287 + t349;
t262 = Ifges(5,4) * t294 + Ifges(5,2) * t293 + Ifges(5,6) * t327;
t263 = Ifges(5,1) * t294 + Ifges(5,4) * t293 + Ifges(5,5) * t327;
t347 = -mrSges(5,1) * t202 + mrSges(5,2) * t203 - Ifges(5,5) * t253 - Ifges(5,6) * t252 - Ifges(5,3) * t326 - pkin(4) * t345 - pkin(9) * t357 - t336 * t167 - t332 * t174 - t294 * t262 + t293 * t263;
t343 = -mrSges(4,1) * t241 + mrSges(4,2) * t242 - Ifges(4,5) * t281 - Ifges(4,6) * t280 - Ifges(4,3) * t326 - pkin(3) * t165 - t307 * t289 + t306 * t290 + t347;
t372 = mrSges(3,1) * t296 - mrSges(3,2) * t297 + Ifges(3,5) * t314 + Ifges(3,6) * t315 + Ifges(3,3) * qJDD(2) + pkin(2) * t157 + (t334 * t304 - t338 * t305) * qJD(1) - t343;
t176 = t336 * t179 + t332 * t182;
t365 = qJD(1) * t338;
t313 = (-mrSges(3,1) * t338 + mrSges(3,2) * t334) * qJD(1);
t318 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t365;
t155 = m(3) * t296 + qJDD(2) * mrSges(3,1) - t314 * mrSges(3,3) + qJD(2) * t318 - t313 * t366 + t157;
t317 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t366;
t359 = -t333 * t162 + t337 * t163;
t156 = m(3) * t297 - qJDD(2) * mrSges(3,2) + t315 * mrSges(3,3) - qJD(2) * t317 + t313 * t365 + t359;
t360 = -t334 * t155 + t338 * t156;
t352 = m(5) * t227 - t252 * mrSges(5,1) + t253 * mrSges(5,2) - t293 * t283 + t294 * t284 + t176;
t261 = Ifges(5,5) * t294 + Ifges(5,6) * t293 + Ifges(5,3) * t327;
t153 = mrSges(5,2) * t227 - mrSges(5,3) * t202 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * t326 - pkin(9) * t176 - t332 * t167 + t336 * t174 + t293 * t261 - t327 * t262;
t158 = -mrSges(5,1) * t227 + mrSges(5,3) * t203 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * t326 - pkin(4) * t176 - t294 * t261 + t327 * t263 - t373;
t288 = Ifges(4,5) * t307 + Ifges(4,6) * t306 + Ifges(4,3) * t327;
t148 = -mrSges(4,1) * t282 + mrSges(4,3) * t242 + Ifges(4,4) * t281 + Ifges(4,2) * t280 + Ifges(4,6) * t326 - pkin(3) * t352 + qJ(4) * t358 + t330 * t153 + t331 * t158 - t307 * t288 + t327 * t290;
t149 = mrSges(4,2) * t282 - mrSges(4,3) * t241 + Ifges(4,1) * t281 + Ifges(4,4) * t280 + Ifges(4,5) * t326 - qJ(4) * t165 + t331 * t153 - t330 * t158 + t306 * t288 - t327 * t289;
t303 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t334 + Ifges(3,6) * t338) * qJD(1);
t308 = -t340 * pkin(7) + t353;
t346 = m(4) * t282 - t280 * mrSges(4,1) + t281 * mrSges(4,2) - t306 * t298 + t307 * t300 + t352;
t144 = -mrSges(3,1) * t308 + mrSges(3,3) * t297 + Ifges(3,4) * t314 + Ifges(3,2) * t315 + Ifges(3,6) * qJDD(2) - pkin(2) * t346 + pkin(8) * t359 + qJD(2) * t305 + t337 * t148 + t333 * t149 - t303 * t366;
t146 = mrSges(3,2) * t308 - mrSges(3,3) * t296 + Ifges(3,1) * t314 + Ifges(3,4) * t315 + Ifges(3,5) * qJDD(2) - pkin(8) * t157 - qJD(2) * t304 - t333 * t148 + t337 * t149 + t303 * t365;
t342 = -m(3) * t308 + t315 * mrSges(3,1) - t314 * mrSges(3,2) - t317 * t366 + t318 * t365 - t346;
t348 = mrSges(2,1) * t320 - mrSges(2,2) * t321 + Ifges(2,3) * qJDD(1) + pkin(1) * t342 + pkin(7) * t360 + t338 * t144 + t334 * t146;
t168 = m(2) * t320 + qJDD(1) * mrSges(2,1) - t340 * mrSges(2,2) + t342;
t152 = t338 * t155 + t334 * t156;
t150 = m(2) * t321 - t340 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t360;
t147 = mrSges(2,1) * g(3) + mrSges(2,3) * t321 + t340 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t152 - t372;
t142 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - t340 * Ifges(2,6) - pkin(7) * t152 - t334 * t144 + t338 * t146;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t339 * t142 - t335 * t147 - pkin(6) * (t335 * t150 + t339 * t168), t142, t146, t149, t153, t174, -t287 * t231 + t349; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t142 + t339 * t147 + pkin(6) * (t339 * t150 - t335 * t168), t147, t144, t148, t158, t167, -t279 * t229 + t351; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t348, t348, t372, -t343, -t347, t373, -t278 * t233 - t350;];
m_new  = t1;
