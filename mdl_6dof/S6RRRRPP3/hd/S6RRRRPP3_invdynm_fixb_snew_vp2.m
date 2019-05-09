% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:07:54
% EndTime: 2019-05-07 18:08:13
% DurationCPUTime: 8.02s
% Computational Cost: add. (137915->390), mult. (276153->448), div. (0->0), fcn. (190419->8), ass. (0->142)
t325 = sin(qJ(3));
t326 = sin(qJ(2));
t328 = cos(qJ(3));
t329 = cos(qJ(2));
t299 = (t325 * t329 + t326 * t328) * qJD(1);
t356 = qJD(1) * qJD(2);
t306 = qJDD(1) * t326 + t329 * t356;
t307 = qJDD(1) * t329 - t326 * t356;
t267 = -qJD(3) * t299 - t306 * t325 + t307 * t328;
t357 = qJD(1) * t329;
t358 = qJD(1) * t326;
t298 = -t325 * t358 + t328 * t357;
t268 = qJD(3) * t298 + t306 * t328 + t307 * t325;
t311 = qJD(2) * pkin(2) - pkin(8) * t358;
t323 = t329 ^ 2;
t331 = qJD(1) ^ 2;
t327 = sin(qJ(1));
t330 = cos(qJ(1));
t312 = t327 * g(1) - t330 * g(2);
t349 = -qJDD(1) * pkin(1) - t312;
t269 = -t307 * pkin(2) + t311 * t358 + (-pkin(8) * t323 - pkin(7)) * t331 + t349;
t321 = qJD(2) + qJD(3);
t200 = (-t298 * t321 - t268) * pkin(9) + (t299 * t321 - t267) * pkin(3) + t269;
t313 = -g(1) * t330 - g(2) * t327;
t301 = -pkin(1) * t331 + qJDD(1) * pkin(7) + t313;
t363 = t326 * t301;
t368 = pkin(2) * t331;
t253 = qJDD(2) * pkin(2) - t306 * pkin(8) - t363 + (pkin(8) * t356 + t326 * t368 - g(3)) * t329;
t289 = -g(3) * t326 + t329 * t301;
t254 = pkin(8) * t307 - qJD(2) * t311 - t323 * t368 + t289;
t207 = t325 * t253 + t328 * t254;
t283 = -pkin(3) * t298 - pkin(9) * t299;
t319 = t321 ^ 2;
t320 = qJDD(2) + qJDD(3);
t204 = -pkin(3) * t319 + pkin(9) * t320 + t283 * t298 + t207;
t324 = sin(qJ(4));
t370 = cos(qJ(4));
t197 = t200 * t370 - t324 * t204;
t286 = t299 * t324 - t321 * t370;
t287 = t299 * t370 + t324 * t321;
t249 = pkin(4) * t286 - qJ(5) * t287;
t266 = qJDD(4) - t267;
t294 = qJD(4) - t298;
t293 = t294 ^ 2;
t195 = -t266 * pkin(4) - t293 * qJ(5) + t287 * t249 + qJDD(5) - t197;
t226 = -t286 * qJD(4) + t268 * t370 + t324 * t320;
t251 = -mrSges(6,2) * t286 - mrSges(6,3) * t287;
t374 = -m(6) * t195 - t226 * mrSges(6,1) - t287 * t251;
t282 = -mrSges(4,1) * t298 + mrSges(4,2) * t299;
t291 = mrSges(4,1) * t321 - mrSges(4,3) * t299;
t248 = -mrSges(7,2) * t287 + mrSges(7,3) * t286;
t250 = mrSges(5,1) * t286 + mrSges(5,2) * t287;
t272 = mrSges(6,1) * t286 - mrSges(6,3) * t294;
t275 = -mrSges(5,2) * t294 - mrSges(5,3) * t286;
t365 = t286 * t294;
t186 = -0.2e1 * qJD(6) * t294 + (t286 * t287 - t266) * qJ(6) + (t226 + t365) * pkin(5) + t195;
t273 = -mrSges(7,1) * t286 + mrSges(7,2) * t294;
t351 = -m(7) * t186 + t266 * mrSges(7,3) + t294 * t273;
t367 = -mrSges(7,1) - mrSges(5,3);
t173 = m(5) * t197 + (-t272 + t275) * t294 + (-t248 - t250) * t287 + (mrSges(5,1) - mrSges(6,2)) * t266 + t367 * t226 + t351 + t374;
t198 = t324 * t200 + t370 * t204;
t225 = qJD(4) * t287 + t268 * t324 - t320 * t370;
t276 = mrSges(5,1) * t294 - mrSges(5,3) * t287;
t342 = -t293 * pkin(4) + t266 * qJ(5) - t286 * t249 + t198;
t371 = -2 * qJD(5);
t193 = t294 * t371 - t342;
t274 = mrSges(6,1) * t287 + mrSges(6,2) * t294;
t270 = pkin(5) * t287 - qJ(6) * t294;
t285 = t286 ^ 2;
t189 = -t225 * pkin(5) - t285 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t270) * t294 + t342;
t271 = mrSges(7,1) * t287 - mrSges(7,3) * t294;
t355 = m(7) * t189 + t266 * mrSges(7,2) + t294 * t271;
t348 = -m(6) * t193 + t266 * mrSges(6,3) + t294 * t274 + t355;
t359 = -t248 - t251;
t176 = m(5) * t198 - t266 * mrSges(5,2) - t294 * t276 + (-t250 + t359) * t286 + (-mrSges(6,1) + t367) * t225 + t348;
t352 = -t173 * t324 + t370 * t176;
t166 = m(4) * t207 - mrSges(4,2) * t320 + mrSges(4,3) * t267 + t282 * t298 - t291 * t321 + t352;
t206 = t328 * t253 - t325 * t254;
t290 = -mrSges(4,2) * t321 + mrSges(4,3) * t298;
t203 = -t320 * pkin(3) - t319 * pkin(9) + t299 * t283 - t206;
t338 = (-t226 + t365) * qJ(5) + t203 + (t294 * pkin(4) + t371) * t287;
t192 = -t285 * pkin(5) + 0.2e1 * qJD(6) * t286 - t287 * t270 + (pkin(4) + qJ(6)) * t225 + t338;
t183 = m(7) * t192 - t226 * mrSges(7,2) + t225 * mrSges(7,3) - t287 * t271 + t286 * t273;
t196 = t225 * pkin(4) + t338;
t343 = -m(6) * t196 + t225 * mrSges(6,2) + t286 * t272 - t183;
t335 = -m(5) * t203 - t225 * mrSges(5,1) - t286 * t275 + (t274 - t276) * t287 + (-mrSges(5,2) + mrSges(6,3)) * t226 + t343;
t171 = m(4) * t206 + t320 * mrSges(4,1) - t268 * mrSges(4,3) - t299 * t282 + t321 * t290 + t335;
t159 = t325 * t166 + t328 * t171;
t288 = -t329 * g(3) - t363;
t296 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t326 + Ifges(3,2) * t329) * qJD(1);
t297 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t326 + Ifges(3,4) * t329) * qJD(1);
t181 = -t226 * mrSges(6,3) - t287 * t274 - t343;
t229 = Ifges(7,5) * t294 + Ifges(7,6) * t286 + Ifges(7,3) * t287;
t231 = Ifges(5,5) * t287 - Ifges(5,6) * t286 + Ifges(5,3) * t294;
t236 = Ifges(6,1) * t294 - Ifges(6,4) * t287 + Ifges(6,5) * t286;
t235 = Ifges(7,1) * t294 + Ifges(7,4) * t286 + Ifges(7,5) * t287;
t347 = mrSges(7,1) * t189 - mrSges(7,3) * t192 - Ifges(7,4) * t266 - Ifges(7,2) * t225 - Ifges(7,6) * t226 - t287 * t235;
t337 = mrSges(6,1) * t193 - mrSges(6,2) * t196 + pkin(5) * (t225 * mrSges(7,1) + t286 * t248 - t355) + qJ(6) * t183 - t347;
t233 = Ifges(6,4) * t294 - Ifges(6,2) * t287 + Ifges(6,6) * t286;
t360 = Ifges(5,1) * t287 - Ifges(5,4) * t286 + Ifges(5,5) * t294 - t233;
t366 = Ifges(5,4) + Ifges(6,6);
t155 = -t337 + (t229 + t360) * t294 + (-t231 - t236) * t287 + (Ifges(5,6) - Ifges(6,5)) * t266 + t366 * t226 + (-Ifges(5,2) - Ifges(6,3)) * t225 - pkin(4) * t181 + mrSges(5,3) * t198 - mrSges(5,1) * t203;
t230 = Ifges(6,5) * t294 - Ifges(6,6) * t287 + Ifges(6,3) * t286;
t234 = Ifges(5,4) * t287 - Ifges(5,2) * t286 + Ifges(5,6) * t294;
t182 = t226 * mrSges(7,1) + t287 * t248 - t351;
t232 = Ifges(7,4) * t294 + Ifges(7,2) * t286 + Ifges(7,6) * t287;
t346 = -mrSges(7,1) * t186 + mrSges(7,2) * t192 - Ifges(7,5) * t266 - Ifges(7,6) * t225 - Ifges(7,3) * t226 - t294 * t232;
t339 = -mrSges(6,1) * t195 + mrSges(6,3) * t196 - pkin(5) * t182 + t346;
t361 = t235 + t236;
t161 = -t339 + (-t234 + t230) * t294 + (-t231 - t361) * t286 + (Ifges(5,5) - Ifges(6,4)) * t266 + (Ifges(5,1) + Ifges(6,2)) * t226 - t366 * t225 - qJ(5) * t181 - mrSges(5,3) * t197 + mrSges(5,2) * t203;
t279 = Ifges(4,4) * t299 + Ifges(4,2) * t298 + Ifges(4,6) * t321;
t280 = Ifges(4,1) * t299 + Ifges(4,4) * t298 + Ifges(4,5) * t321;
t340 = -mrSges(4,1) * t206 + mrSges(4,2) * t207 - Ifges(4,5) * t268 - Ifges(4,6) * t267 - Ifges(4,3) * t320 - pkin(3) * t335 - pkin(9) * t352 - t370 * t155 - t324 * t161 - t299 * t279 + t298 * t280;
t373 = mrSges(3,1) * t288 - mrSges(3,2) * t289 + Ifges(3,5) * t306 + Ifges(3,6) * t307 + Ifges(3,3) * qJDD(2) + pkin(2) * t159 + (t296 * t326 - t297 * t329) * qJD(1) - t340;
t345 = -mrSges(7,2) * t189 + mrSges(7,3) * t186 - Ifges(7,1) * t266 - Ifges(7,4) * t225 - Ifges(7,5) * t226 - t286 * t229;
t336 = -mrSges(6,2) * t195 + mrSges(6,3) * t193 - Ifges(6,1) * t266 + Ifges(6,4) * t226 - Ifges(6,5) * t225 + qJ(6) * t182 + t287 * t230 + t345;
t372 = t286 * t360 + (t234 - t232) * t287 + mrSges(5,1) * t197 - mrSges(5,2) * t198 + Ifges(5,5) * t226 - Ifges(5,6) * t225 + Ifges(5,3) * t266 + pkin(4) * (-t266 * mrSges(6,2) - t294 * t272 - t182 + t374) + qJ(5) * (t359 * t286 + (-mrSges(6,1) - mrSges(7,1)) * t225 + t348) - t336;
t364 = t287 * t232;
t168 = t370 * t173 + t324 * t176;
t305 = (-mrSges(3,1) * t329 + mrSges(3,2) * t326) * qJD(1);
t310 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t357;
t157 = m(3) * t288 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t306 + qJD(2) * t310 - t305 * t358 + t159;
t309 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t358;
t353 = t328 * t166 - t171 * t325;
t158 = m(3) * t289 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t307 - qJD(2) * t309 + t305 * t357 + t353;
t354 = -t157 * t326 + t329 * t158;
t278 = Ifges(4,5) * t299 + Ifges(4,6) * t298 + Ifges(4,3) * t321;
t149 = mrSges(4,2) * t269 - mrSges(4,3) * t206 + Ifges(4,1) * t268 + Ifges(4,4) * t267 + Ifges(4,5) * t320 - pkin(9) * t168 - t324 * t155 + t161 * t370 + t298 * t278 - t321 * t279;
t153 = -mrSges(4,1) * t269 + mrSges(4,3) * t207 + Ifges(4,4) * t268 + Ifges(4,2) * t267 + Ifges(4,6) * t320 - pkin(3) * t168 - t299 * t278 + t321 * t280 - t372;
t295 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t326 + Ifges(3,6) * t329) * qJD(1);
t300 = -t331 * pkin(7) + t349;
t341 = m(4) * t269 - t267 * mrSges(4,1) + mrSges(4,2) * t268 - t298 * t290 + t291 * t299 + t168;
t145 = -mrSges(3,1) * t300 + mrSges(3,3) * t289 + Ifges(3,4) * t306 + Ifges(3,2) * t307 + Ifges(3,6) * qJDD(2) - pkin(2) * t341 + pkin(8) * t353 + qJD(2) * t297 + t325 * t149 + t328 * t153 - t295 * t358;
t148 = mrSges(3,2) * t300 - mrSges(3,3) * t288 + Ifges(3,1) * t306 + Ifges(3,4) * t307 + Ifges(3,5) * qJDD(2) - pkin(8) * t159 - qJD(2) * t296 + t149 * t328 - t153 * t325 + t295 * t357;
t334 = -m(3) * t300 + t307 * mrSges(3,1) - mrSges(3,2) * t306 - t309 * t358 + t310 * t357 - t341;
t344 = mrSges(2,1) * t312 - mrSges(2,2) * t313 + Ifges(2,3) * qJDD(1) + pkin(1) * t334 + pkin(7) * t354 + t329 * t145 + t326 * t148;
t162 = m(2) * t312 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t331 + t334;
t152 = t157 * t329 + t158 * t326;
t150 = m(2) * t313 - mrSges(2,1) * t331 - qJDD(1) * mrSges(2,2) + t354;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t313 + t331 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t152 - t373;
t143 = -mrSges(2,2) * g(3) - mrSges(2,3) * t312 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t331 - pkin(7) * t152 - t145 * t326 + t148 * t329;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t330 * t143 - t327 * t146 - pkin(6) * (t150 * t327 + t162 * t330), t143, t148, t149, t161, -t286 * t233 - t336 - t364, -t345 - t364; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t327 * t143 + t330 * t146 + pkin(6) * (t150 * t330 - t162 * t327), t146, t145, t153, t155, Ifges(6,4) * t266 - Ifges(6,2) * t226 + Ifges(6,6) * t225 - t294 * t230 + t286 * t361 + t339, -t294 * t229 - t347; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t344, t344, t373, -t340, t372, t337 + (-t229 + t233) * t294 + Ifges(6,3) * t225 - Ifges(6,6) * t226 + Ifges(6,5) * t266 + t287 * t236, -t286 * t235 - t346;];
m_new  = t1;
