% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 14:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:56:00
% EndTime: 2019-05-06 13:57:38
% DurationCPUTime: 71.76s
% Computational Cost: add. (1214793->396), mult. (3218897->519), div. (0->0), fcn. (2566244->14), ass. (0->163)
t385 = -2 * qJD(3);
t345 = sin(pkin(11));
t348 = cos(pkin(11));
t352 = sin(qJ(2));
t355 = cos(qJ(2));
t346 = sin(pkin(6));
t377 = qJD(1) * t346;
t320 = (t345 * t352 - t348 * t355) * t377;
t375 = qJD(1) * qJD(2);
t329 = (qJDD(1) * t352 + t355 * t375) * t346;
t349 = cos(pkin(6));
t338 = t349 * qJDD(1) + qJDD(2);
t339 = t349 * qJD(1) + qJD(2);
t353 = sin(qJ(1));
t356 = cos(qJ(1));
t335 = t353 * g(1) - t356 * g(2);
t357 = qJD(1) ^ 2;
t383 = pkin(8) * t346;
t326 = qJDD(1) * pkin(1) + t357 * t383 + t335;
t336 = -t356 * g(1) - t353 * g(2);
t327 = -t357 * pkin(1) + qJDD(1) * t383 + t336;
t378 = t349 * t355;
t367 = t326 * t378 - t352 * t327;
t382 = t346 ^ 2 * t357;
t263 = t338 * pkin(2) - t329 * qJ(3) + (pkin(2) * t352 * t382 + (qJ(3) * qJD(1) * t339 - g(3)) * t346) * t355 + t367;
t379 = t349 * t352;
t381 = t346 * t352;
t295 = -g(3) * t381 + t326 * t379 + t355 * t327;
t373 = t352 * t377;
t323 = t339 * pkin(2) - qJ(3) * t373;
t330 = (qJDD(1) * t355 - t352 * t375) * t346;
t374 = t355 ^ 2 * t382;
t267 = -pkin(2) * t374 + t330 * qJ(3) - t339 * t323 + t295;
t321 = (t345 * t355 + t348 * t352) * t377;
t243 = t348 * t263 - t345 * t267 + t321 * t385;
t384 = cos(qJ(4));
t380 = t346 * t355;
t244 = t345 * t263 + t348 * t267 + t320 * t385;
t296 = t320 * mrSges(4,1) + t321 * mrSges(4,2);
t302 = -t345 * t329 + t348 * t330;
t308 = t339 * mrSges(4,1) - t321 * mrSges(4,3);
t297 = t320 * pkin(3) - t321 * pkin(9);
t337 = t339 ^ 2;
t236 = -t337 * pkin(3) + t338 * pkin(9) - t320 * t297 + t244;
t312 = -t349 * g(3) - t346 * t326;
t278 = -t330 * pkin(2) - qJ(3) * t374 + t323 * t373 + qJDD(3) + t312;
t303 = t348 * t329 + t345 * t330;
t246 = (t320 * t339 - t303) * pkin(9) + (t321 * t339 - t302) * pkin(3) + t278;
t351 = sin(qJ(4));
t226 = t384 * t236 + t351 * t246;
t305 = t351 * t321 - t384 * t339;
t306 = t384 * t321 + t351 * t339;
t280 = t305 * pkin(4) - t306 * qJ(5);
t301 = qJDD(4) - t302;
t319 = qJD(4) + t320;
t318 = t319 ^ 2;
t221 = -t318 * pkin(4) + t301 * qJ(5) - t305 * t280 + t226;
t235 = -t338 * pkin(3) - t337 * pkin(9) + t321 * t297 - t243;
t275 = t306 * qJD(4) + t351 * t303 - t384 * t338;
t276 = -t305 * qJD(4) + t384 * t303 + t351 * t338;
t224 = (t305 * t319 - t276) * qJ(5) + (t306 * t319 + t275) * pkin(4) + t235;
t344 = sin(pkin(12));
t347 = cos(pkin(12));
t287 = t347 * t306 + t344 * t319;
t216 = -0.2e1 * qJD(5) * t287 - t344 * t221 + t347 * t224;
t255 = t347 * t276 + t344 * t301;
t286 = -t344 * t306 + t347 * t319;
t214 = (t286 * t305 - t255) * pkin(10) + (t286 * t287 + t275) * pkin(5) + t216;
t217 = 0.2e1 * qJD(5) * t286 + t347 * t221 + t344 * t224;
t254 = -t344 * t276 + t347 * t301;
t266 = t305 * pkin(5) - t287 * pkin(10);
t285 = t286 ^ 2;
t215 = -t285 * pkin(5) + t254 * pkin(10) - t305 * t266 + t217;
t350 = sin(qJ(6));
t354 = cos(qJ(6));
t212 = t354 * t214 - t350 * t215;
t256 = t354 * t286 - t350 * t287;
t232 = t256 * qJD(6) + t350 * t254 + t354 * t255;
t257 = t350 * t286 + t354 * t287;
t241 = -t256 * mrSges(7,1) + t257 * mrSges(7,2);
t304 = qJD(6) + t305;
t247 = -t304 * mrSges(7,2) + t256 * mrSges(7,3);
t274 = qJDD(6) + t275;
t208 = m(7) * t212 + t274 * mrSges(7,1) - t232 * mrSges(7,3) - t257 * t241 + t304 * t247;
t213 = t350 * t214 + t354 * t215;
t231 = -t257 * qJD(6) + t354 * t254 - t350 * t255;
t248 = t304 * mrSges(7,1) - t257 * mrSges(7,3);
t209 = m(7) * t213 - t274 * mrSges(7,2) + t231 * mrSges(7,3) + t256 * t241 - t304 * t248;
t200 = t354 * t208 + t350 * t209;
t258 = -t286 * mrSges(6,1) + t287 * mrSges(6,2);
t264 = -t305 * mrSges(6,2) + t286 * mrSges(6,3);
t198 = m(6) * t216 + t275 * mrSges(6,1) - t255 * mrSges(6,3) - t287 * t258 + t305 * t264 + t200;
t265 = t305 * mrSges(6,1) - t287 * mrSges(6,3);
t369 = -t350 * t208 + t354 * t209;
t199 = m(6) * t217 - t275 * mrSges(6,2) + t254 * mrSges(6,3) + t286 * t258 - t305 * t265 + t369;
t196 = -t344 * t198 + t347 * t199;
t281 = t305 * mrSges(5,1) + t306 * mrSges(5,2);
t289 = t319 * mrSges(5,1) - t306 * mrSges(5,3);
t193 = m(5) * t226 - t301 * mrSges(5,2) - t275 * mrSges(5,3) - t305 * t281 - t319 * t289 + t196;
t225 = -t351 * t236 + t384 * t246;
t220 = -t301 * pkin(4) - t318 * qJ(5) + t306 * t280 + qJDD(5) - t225;
t218 = -t254 * pkin(5) - t285 * pkin(10) + t287 * t266 + t220;
t365 = m(7) * t218 - t231 * mrSges(7,1) + t232 * mrSges(7,2) - t256 * t247 + t257 * t248;
t210 = -m(6) * t220 + t254 * mrSges(6,1) - t255 * mrSges(6,2) + t286 * t264 - t287 * t265 - t365;
t288 = -t319 * mrSges(5,2) - t305 * mrSges(5,3);
t207 = m(5) * t225 + t301 * mrSges(5,1) - t276 * mrSges(5,3) - t306 * t281 + t319 * t288 + t210;
t370 = t384 * t193 - t351 * t207;
t183 = m(4) * t244 - t338 * mrSges(4,2) + t302 * mrSges(4,3) - t320 * t296 - t339 * t308 + t370;
t307 = -t339 * mrSges(4,2) - t320 * mrSges(4,3);
t195 = t347 * t198 + t344 * t199;
t360 = -m(5) * t235 - t275 * mrSges(5,1) - t276 * mrSges(5,2) - t305 * t288 - t306 * t289 - t195;
t190 = m(4) * t243 + t338 * mrSges(4,1) - t303 * mrSges(4,3) - t321 * t296 + t339 * t307 + t360;
t178 = t345 * t183 + t348 * t190;
t186 = t351 * t193 + t384 * t207;
t372 = t355 * t377;
t294 = -g(3) * t380 + t367;
t325 = -t339 * mrSges(3,2) + mrSges(3,3) * t372;
t328 = (-mrSges(3,1) * t355 + mrSges(3,2) * t352) * t377;
t176 = m(3) * t294 + t338 * mrSges(3,1) - t329 * mrSges(3,3) + t339 * t325 - t328 * t373 + t178;
t324 = t339 * mrSges(3,1) - mrSges(3,3) * t373;
t371 = t348 * t183 - t345 * t190;
t177 = m(3) * t295 - t338 * mrSges(3,2) + t330 * mrSges(3,3) - t339 * t324 + t328 * t372 + t371;
t170 = -t352 * t176 + t355 * t177;
t363 = m(4) * t278 - t302 * mrSges(4,1) + t303 * mrSges(4,2) + t320 * t307 + t321 * t308 + t186;
t184 = m(3) * t312 - t330 * mrSges(3,1) + t329 * mrSges(3,2) + (t324 * t352 - t325 * t355) * t377 + t363;
t166 = t176 * t378 + t177 * t379 - t346 * t184;
t237 = Ifges(7,5) * t257 + Ifges(7,6) * t256 + Ifges(7,3) * t304;
t239 = Ifges(7,1) * t257 + Ifges(7,4) * t256 + Ifges(7,5) * t304;
t201 = -mrSges(7,1) * t218 + mrSges(7,3) * t213 + Ifges(7,4) * t232 + Ifges(7,2) * t231 + Ifges(7,6) * t274 - t257 * t237 + t304 * t239;
t238 = Ifges(7,4) * t257 + Ifges(7,2) * t256 + Ifges(7,6) * t304;
t202 = mrSges(7,2) * t218 - mrSges(7,3) * t212 + Ifges(7,1) * t232 + Ifges(7,4) * t231 + Ifges(7,5) * t274 + t256 * t237 - t304 * t238;
t249 = Ifges(6,5) * t287 + Ifges(6,6) * t286 + Ifges(6,3) * t305;
t251 = Ifges(6,1) * t287 + Ifges(6,4) * t286 + Ifges(6,5) * t305;
t187 = -mrSges(6,1) * t220 + mrSges(6,3) * t217 + Ifges(6,4) * t255 + Ifges(6,2) * t254 + Ifges(6,6) * t275 - pkin(5) * t365 + pkin(10) * t369 + t354 * t201 + t350 * t202 - t287 * t249 + t305 * t251;
t250 = Ifges(6,4) * t287 + Ifges(6,2) * t286 + Ifges(6,6) * t305;
t188 = mrSges(6,2) * t220 - mrSges(6,3) * t216 + Ifges(6,1) * t255 + Ifges(6,4) * t254 + Ifges(6,5) * t275 - pkin(10) * t200 - t350 * t201 + t354 * t202 + t286 * t249 - t305 * t250;
t268 = Ifges(5,5) * t306 - Ifges(5,6) * t305 + Ifges(5,3) * t319;
t269 = Ifges(5,4) * t306 - Ifges(5,2) * t305 + Ifges(5,6) * t319;
t172 = mrSges(5,2) * t235 - mrSges(5,3) * t225 + Ifges(5,1) * t276 - Ifges(5,4) * t275 + Ifges(5,5) * t301 - qJ(5) * t195 - t344 * t187 + t347 * t188 - t305 * t268 - t319 * t269;
t270 = Ifges(5,1) * t306 - Ifges(5,4) * t305 + Ifges(5,5) * t319;
t362 = -mrSges(7,1) * t212 + mrSges(7,2) * t213 - Ifges(7,5) * t232 - Ifges(7,6) * t231 - Ifges(7,3) * t274 - t257 * t238 + t256 * t239;
t359 = -mrSges(6,1) * t216 + mrSges(6,2) * t217 - Ifges(6,5) * t255 - Ifges(6,6) * t254 - pkin(5) * t200 - t287 * t250 + t286 * t251 + t362;
t180 = (-Ifges(5,2) - Ifges(6,3)) * t275 + t319 * t270 - t306 * t268 + Ifges(5,6) * t301 + Ifges(5,4) * t276 - mrSges(5,1) * t235 + mrSges(5,3) * t226 + t359 - pkin(4) * t195;
t290 = Ifges(4,5) * t321 - Ifges(4,6) * t320 + Ifges(4,3) * t339;
t291 = Ifges(4,4) * t321 - Ifges(4,2) * t320 + Ifges(4,6) * t339;
t162 = mrSges(4,2) * t278 - mrSges(4,3) * t243 + Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * t338 - pkin(9) * t186 + t384 * t172 - t351 * t180 - t320 * t290 - t339 * t291;
t292 = Ifges(4,1) * t321 - Ifges(4,4) * t320 + Ifges(4,5) * t339;
t358 = mrSges(5,1) * t225 - mrSges(5,2) * t226 + Ifges(5,5) * t276 - Ifges(5,6) * t275 + Ifges(5,3) * t301 + pkin(4) * t210 + qJ(5) * t196 + t347 * t187 + t344 * t188 + t306 * t269 + t305 * t270;
t167 = -mrSges(4,1) * t278 + mrSges(4,3) * t244 + Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * t338 - pkin(3) * t186 - t321 * t290 + t339 * t292 - t358;
t309 = Ifges(3,3) * t339 + (Ifges(3,5) * t352 + Ifges(3,6) * t355) * t377;
t311 = Ifges(3,5) * t339 + (Ifges(3,1) * t352 + Ifges(3,4) * t355) * t377;
t157 = -mrSges(3,1) * t312 + mrSges(3,3) * t295 + Ifges(3,4) * t329 + Ifges(3,2) * t330 + Ifges(3,6) * t338 - pkin(2) * t363 + qJ(3) * t371 + t345 * t162 + t348 * t167 - t309 * t373 + t339 * t311;
t310 = Ifges(3,6) * t339 + (Ifges(3,4) * t352 + Ifges(3,2) * t355) * t377;
t159 = mrSges(3,2) * t312 - mrSges(3,3) * t294 + Ifges(3,1) * t329 + Ifges(3,4) * t330 + Ifges(3,5) * t338 - qJ(3) * t178 + t348 * t162 - t345 * t167 + t309 * t372 - t339 * t310;
t361 = mrSges(4,1) * t243 - mrSges(4,2) * t244 + Ifges(4,5) * t303 + Ifges(4,6) * t302 + Ifges(4,3) * t338 + pkin(3) * t360 + pkin(9) * t370 + t351 * t172 + t384 * t180 + t321 * t291 + t320 * t292;
t161 = Ifges(3,3) * t338 + Ifges(3,5) * t329 + Ifges(3,6) * t330 + mrSges(3,1) * t294 - mrSges(3,2) * t295 + pkin(2) * t178 + t361 + (t310 * t352 - t311 * t355) * t377;
t364 = mrSges(2,1) * t335 - mrSges(2,2) * t336 + Ifges(2,3) * qJDD(1) + pkin(1) * t166 + t157 * t380 + t159 * t381 + t349 * t161 + t170 * t383;
t168 = m(2) * t336 - t357 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t170;
t165 = t349 * t184 + (t176 * t355 + t177 * t352) * t346;
t163 = m(2) * t335 + qJDD(1) * mrSges(2,1) - t357 * mrSges(2,2) + t166;
t155 = -mrSges(2,2) * g(3) - mrSges(2,3) * t335 + Ifges(2,5) * qJDD(1) - t357 * Ifges(2,6) - t352 * t157 + t355 * t159 + (-t165 * t346 - t166 * t349) * pkin(8);
t154 = mrSges(2,1) * g(3) + mrSges(2,3) * t336 + t357 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t165 - t346 * t161 + (pkin(8) * t170 + t157 * t355 + t159 * t352) * t349;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t356 * t155 - t353 * t154 - pkin(7) * (t356 * t163 + t353 * t168), t155, t159, t162, t172, t188, t202; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t353 * t155 + t356 * t154 + pkin(7) * (-t353 * t163 + t356 * t168), t154, t157, t167, t180, t187, t201; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t364, t364, t161, t361, t358, Ifges(6,3) * t275 - t359, -t362;];
m_new  = t1;
