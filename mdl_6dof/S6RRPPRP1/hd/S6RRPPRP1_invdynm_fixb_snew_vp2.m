% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-05-06 09:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:02:26
% EndTime: 2019-05-06 09:02:58
% DurationCPUTime: 17.71s
% Computational Cost: add. (286531->382), mult. (675962->470), div. (0->0), fcn. (478159->10), ass. (0->141)
t373 = -2 * qJD(3);
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t361 = qJD(1) * qJD(2);
t317 = qJDD(1) * t334 + t336 * t361;
t335 = sin(qJ(1));
t337 = cos(qJ(1));
t324 = -g(1) * t337 - g(2) * t335;
t339 = qJD(1) ^ 2;
t312 = -pkin(1) * t339 + qJDD(1) * pkin(7) + t324;
t367 = t312 * t334;
t370 = pkin(2) * t339;
t261 = qJDD(2) * pkin(2) - qJ(3) * t317 - t367 + (qJ(3) * t361 + t334 * t370 - g(3)) * t336;
t297 = -g(3) * t334 + t336 * t312;
t318 = qJDD(1) * t336 - t334 * t361;
t364 = qJD(1) * t334;
t320 = qJD(2) * pkin(2) - qJ(3) * t364;
t329 = t336 ^ 2;
t263 = qJ(3) * t318 - qJD(2) * t320 - t329 * t370 + t297;
t331 = sin(pkin(9));
t368 = cos(pkin(9));
t306 = (t331 * t336 + t368 * t334) * qJD(1);
t239 = t368 * t261 - t331 * t263 + t306 * t373;
t363 = qJD(1) * t336;
t305 = t331 * t364 - t368 * t363;
t240 = t331 * t261 + t368 * t263 + t305 * t373;
t278 = mrSges(4,1) * t305 + mrSges(4,2) * t306;
t288 = t317 * t331 - t368 * t318;
t299 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t306;
t277 = pkin(3) * t305 - qJ(4) * t306;
t338 = qJD(2) ^ 2;
t213 = -pkin(3) * t338 + qJDD(2) * qJ(4) - t277 * t305 + t240;
t323 = g(1) * t335 - t337 * g(2);
t351 = -qJDD(1) * pkin(1) - t323;
t265 = -pkin(2) * t318 + qJDD(3) + t320 * t364 + (-qJ(3) * t329 - pkin(7)) * t339 + t351;
t289 = t368 * t317 + t331 * t318;
t217 = (qJD(2) * t305 - t289) * qJ(4) + (qJD(2) * t306 + t288) * pkin(3) + t265;
t330 = sin(pkin(10));
t332 = cos(pkin(10));
t295 = qJD(2) * t330 + t306 * t332;
t206 = -0.2e1 * qJD(4) * t295 - t213 * t330 + t332 * t217;
t272 = qJDD(2) * t330 + t289 * t332;
t294 = qJD(2) * t332 - t306 * t330;
t203 = (t294 * t305 - t272) * pkin(8) + (t294 * t295 + t288) * pkin(4) + t206;
t207 = 0.2e1 * qJD(4) * t294 + t332 * t213 + t330 * t217;
t268 = pkin(4) * t305 - pkin(8) * t295;
t271 = qJDD(2) * t332 - t289 * t330;
t293 = t294 ^ 2;
t205 = -pkin(4) * t293 + pkin(8) * t271 - t268 * t305 + t207;
t333 = sin(qJ(5));
t371 = cos(qJ(5));
t199 = t333 * t203 + t371 * t205;
t257 = t333 * t294 + t371 * t295;
t227 = qJD(5) * t257 - t371 * t271 + t272 * t333;
t304 = qJD(5) + t305;
t248 = mrSges(6,1) * t304 - mrSges(6,3) * t257;
t256 = -t371 * t294 + t295 * t333;
t287 = qJDD(5) + t288;
t241 = pkin(5) * t256 - qJ(6) * t257;
t303 = t304 ^ 2;
t194 = -pkin(5) * t303 + qJ(6) * t287 + 0.2e1 * qJD(6) * t304 - t241 * t256 + t199;
t249 = -mrSges(7,1) * t304 + mrSges(7,2) * t257;
t360 = m(7) * t194 + t287 * mrSges(7,3) + t304 * t249;
t242 = mrSges(7,1) * t256 - mrSges(7,3) * t257;
t365 = -mrSges(6,1) * t256 - mrSges(6,2) * t257 - t242;
t369 = -mrSges(6,3) - mrSges(7,2);
t180 = m(6) * t199 - mrSges(6,2) * t287 + t227 * t369 - t248 * t304 + t256 * t365 + t360;
t198 = t371 * t203 - t333 * t205;
t228 = -t256 * qJD(5) + t333 * t271 + t371 * t272;
t247 = -mrSges(6,2) * t304 - mrSges(6,3) * t256;
t196 = -t287 * pkin(5) - t303 * qJ(6) + t257 * t241 + qJDD(6) - t198;
t246 = -mrSges(7,2) * t256 + mrSges(7,3) * t304;
t355 = -m(7) * t196 + t287 * mrSges(7,1) + t304 * t246;
t182 = m(6) * t198 + mrSges(6,1) * t287 + t228 * t369 + t247 * t304 + t257 * t365 + t355;
t175 = t333 * t180 + t371 * t182;
t262 = -mrSges(5,1) * t294 + mrSges(5,2) * t295;
t266 = -mrSges(5,2) * t305 + mrSges(5,3) * t294;
t173 = m(5) * t206 + mrSges(5,1) * t288 - mrSges(5,3) * t272 - t262 * t295 + t266 * t305 + t175;
t267 = mrSges(5,1) * t305 - mrSges(5,3) * t295;
t356 = t371 * t180 - t333 * t182;
t174 = m(5) * t207 - mrSges(5,2) * t288 + mrSges(5,3) * t271 + t262 * t294 - t267 * t305 + t356;
t357 = -t173 * t330 + t332 * t174;
t166 = m(4) * t240 - qJDD(2) * mrSges(4,2) - mrSges(4,3) * t288 - qJD(2) * t299 - t278 * t305 + t357;
t298 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t305;
t212 = -qJDD(2) * pkin(3) - t338 * qJ(4) + t306 * t277 + qJDD(4) - t239;
t208 = -t271 * pkin(4) - t293 * pkin(8) + t295 * t268 + t212;
t201 = -0.2e1 * qJD(6) * t257 + (t256 * t304 - t228) * qJ(6) + (t257 * t304 + t227) * pkin(5) + t208;
t191 = m(7) * t201 + t227 * mrSges(7,1) - t228 * mrSges(7,3) + t256 * t246 - t257 * t249;
t345 = m(6) * t208 + t227 * mrSges(6,1) + t228 * mrSges(6,2) + t256 * t247 + t257 * t248 + t191;
t342 = -m(5) * t212 + t271 * mrSges(5,1) - mrSges(5,2) * t272 + t294 * t266 - t267 * t295 - t345;
t184 = m(4) * t239 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t289 + qJD(2) * t298 - t278 * t306 + t342;
t161 = t331 * t166 + t368 * t184;
t296 = -g(3) * t336 - t367;
t308 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t334 + Ifges(3,2) * t336) * qJD(1);
t309 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t334 + Ifges(3,4) * t336) * qJD(1);
t233 = Ifges(7,1) * t257 + Ifges(7,4) * t304 + Ifges(7,5) * t256;
t234 = Ifges(6,1) * t257 - Ifges(6,4) * t256 + Ifges(6,5) * t304;
t354 = -mrSges(7,1) * t201 + mrSges(7,2) * t194;
t231 = Ifges(7,4) * t257 + Ifges(7,2) * t304 + Ifges(7,6) * t256;
t366 = -Ifges(6,5) * t257 + Ifges(6,6) * t256 - Ifges(6,3) * t304 - t231;
t176 = -mrSges(6,1) * t208 + mrSges(6,3) * t199 - pkin(5) * t191 + (t233 + t234) * t304 + (Ifges(6,6) - Ifges(7,6)) * t287 + t366 * t257 + (Ifges(6,4) - Ifges(7,5)) * t228 + (-Ifges(6,2) - Ifges(7,3)) * t227 + t354;
t232 = Ifges(6,4) * t257 - Ifges(6,2) * t256 + Ifges(6,6) * t304;
t229 = Ifges(7,5) * t257 + Ifges(7,6) * t304 + Ifges(7,3) * t256;
t350 = mrSges(7,2) * t196 - mrSges(7,3) * t201 + Ifges(7,1) * t228 + Ifges(7,4) * t287 + Ifges(7,5) * t227 + t304 * t229;
t177 = mrSges(6,2) * t208 - mrSges(6,3) * t198 + Ifges(6,1) * t228 - Ifges(6,4) * t227 + Ifges(6,5) * t287 - qJ(6) * t191 - t232 * t304 + t256 * t366 + t350;
t250 = Ifges(5,5) * t295 + Ifges(5,6) * t294 + Ifges(5,3) * t305;
t252 = Ifges(5,1) * t295 + Ifges(5,4) * t294 + Ifges(5,5) * t305;
t155 = -mrSges(5,1) * t212 + mrSges(5,3) * t207 + Ifges(5,4) * t272 + Ifges(5,2) * t271 + Ifges(5,6) * t288 - pkin(4) * t345 + pkin(8) * t356 + t176 * t371 + t333 * t177 - t295 * t250 + t305 * t252;
t251 = Ifges(5,4) * t295 + Ifges(5,2) * t294 + Ifges(5,6) * t305;
t157 = mrSges(5,2) * t212 - mrSges(5,3) * t206 + Ifges(5,1) * t272 + Ifges(5,4) * t271 + Ifges(5,5) * t288 - pkin(8) * t175 - t333 * t176 + t177 * t371 + t294 * t250 - t305 * t251;
t274 = Ifges(4,4) * t306 - Ifges(4,2) * t305 + Ifges(4,6) * qJD(2);
t275 = Ifges(4,1) * t306 - Ifges(4,4) * t305 + Ifges(4,5) * qJD(2);
t346 = -mrSges(4,1) * t239 + mrSges(4,2) * t240 - Ifges(4,5) * t289 + Ifges(4,6) * t288 - Ifges(4,3) * qJDD(2) - pkin(3) * t342 - qJ(4) * t357 - t332 * t155 - t330 * t157 - t306 * t274 - t305 * t275;
t372 = mrSges(3,1) * t296 - mrSges(3,2) * t297 + Ifges(3,5) * t317 + Ifges(3,6) * t318 + Ifges(3,3) * qJDD(2) + pkin(2) * t161 + (t334 * t308 - t336 * t309) * qJD(1) - t346;
t168 = t332 * t173 + t330 * t174;
t316 = (-mrSges(3,1) * t336 + mrSges(3,2) * t334) * qJD(1);
t322 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t159 = m(3) * t296 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t317 + qJD(2) * t322 - t316 * t364 + t161;
t321 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t358 = t368 * t166 - t184 * t331;
t160 = m(3) * t297 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t318 - qJD(2) * t321 + t316 * t363 + t358;
t359 = -t159 * t334 + t336 * t160;
t273 = Ifges(4,5) * t306 - Ifges(4,6) * t305 + Ifges(4,3) * qJD(2);
t149 = mrSges(4,2) * t265 - mrSges(4,3) * t239 + Ifges(4,1) * t289 - Ifges(4,4) * t288 + Ifges(4,5) * qJDD(2) - qJ(4) * t168 - qJD(2) * t274 - t155 * t330 + t157 * t332 - t273 * t305;
t348 = mrSges(7,1) * t196 - mrSges(7,3) * t194 - Ifges(7,4) * t228 - Ifges(7,2) * t287 - Ifges(7,6) * t227 + t257 * t229 - t256 * t233;
t343 = mrSges(6,2) * t199 - t256 * t234 - qJ(6) * (-mrSges(7,2) * t227 - t242 * t256 + t360) - pkin(5) * (-mrSges(7,2) * t228 - t242 * t257 + t355) - mrSges(6,1) * t198 - t257 * t232 + Ifges(6,6) * t227 - Ifges(6,5) * t228 - Ifges(6,3) * t287 + t348;
t340 = mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t272 + Ifges(5,6) * t271 + pkin(4) * t175 + t295 * t251 - t294 * t252 - t343;
t153 = Ifges(4,6) * qJDD(2) + (-Ifges(5,3) - Ifges(4,2)) * t288 - t306 * t273 + Ifges(4,4) * t289 + qJD(2) * t275 - mrSges(4,1) * t265 + mrSges(4,3) * t240 - t340 - pkin(3) * t168;
t307 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t334 + Ifges(3,6) * t336) * qJD(1);
t311 = -pkin(7) * t339 + t351;
t347 = m(4) * t265 + t288 * mrSges(4,1) + mrSges(4,2) * t289 + t305 * t298 + t299 * t306 + t168;
t145 = -mrSges(3,1) * t311 + mrSges(3,3) * t297 + Ifges(3,4) * t317 + Ifges(3,2) * t318 + Ifges(3,6) * qJDD(2) - pkin(2) * t347 + qJ(3) * t358 + qJD(2) * t309 + t331 * t149 + t153 * t368 - t307 * t364;
t148 = mrSges(3,2) * t311 - mrSges(3,3) * t296 + Ifges(3,1) * t317 + Ifges(3,4) * t318 + Ifges(3,5) * qJDD(2) - qJ(3) * t161 - qJD(2) * t308 + t149 * t368 - t331 * t153 + t307 * t363;
t344 = -m(3) * t311 + t318 * mrSges(3,1) - mrSges(3,2) * t317 - t321 * t364 + t322 * t363 - t347;
t349 = mrSges(2,1) * t323 - mrSges(2,2) * t324 + Ifges(2,3) * qJDD(1) + pkin(1) * t344 + pkin(7) * t359 + t336 * t145 + t334 * t148;
t162 = m(2) * t323 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t339 + t344;
t152 = t159 * t336 + t160 * t334;
t150 = m(2) * t324 - mrSges(2,1) * t339 - qJDD(1) * mrSges(2,2) + t359;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t324 + t339 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t152 - t372;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t323 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t339 - pkin(7) * t152 - t145 * t334 + t148 * t336;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t143 - t335 * t146 - pkin(6) * (t150 * t335 + t162 * t337), t143, t148, t149, t157, t177, -t231 * t256 + t350; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t335 * t143 + t337 * t146 + pkin(6) * (t150 * t337 - t335 * t162), t146, t145, t153, t155, t176, -t348; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t349, t349, t372, -t346, Ifges(5,3) * t288 + t340, -t343, Ifges(7,5) * t228 + Ifges(7,6) * t287 + Ifges(7,3) * t227 + t231 * t257 - t233 * t304 - t354;];
m_new  = t1;
