% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:39:56
% EndTime: 2019-05-06 01:40:30
% DurationCPUTime: 17.05s
% Computational Cost: add. (288384->362), mult. (678201->439), div. (0->0), fcn. (512915->10), ass. (0->146)
t334 = qJD(1) ^ 2;
t324 = sin(pkin(10));
t366 = qJD(1) * t324;
t325 = cos(pkin(10));
t365 = qJD(1) * t325;
t329 = sin(qJ(1));
t332 = cos(qJ(1));
t311 = -g(1) * t332 - g(2) * t329;
t304 = -pkin(1) * t334 + qJDD(1) * qJ(2) + t311;
t363 = qJD(1) * qJD(2);
t359 = -t325 * g(3) - 0.2e1 * t324 * t363;
t372 = pkin(2) * t325;
t267 = (-pkin(7) * qJDD(1) + t334 * t372 - t304) * t324 + t359;
t289 = -g(3) * t324 + (t304 + 0.2e1 * t363) * t325;
t361 = qJDD(1) * t325;
t321 = t325 ^ 2;
t369 = t321 * t334;
t268 = -pkin(2) * t369 + pkin(7) * t361 + t289;
t328 = sin(qJ(3));
t331 = cos(qJ(3));
t242 = t328 * t267 + t331 * t268;
t302 = -t328 * t366 + t331 * t365;
t347 = t324 * t331 + t325 * t328;
t303 = t347 * qJD(1);
t279 = -mrSges(4,1) * t302 + mrSges(4,2) * t303;
t299 = t303 * qJD(3);
t362 = qJDD(1) * t324;
t286 = -t328 * t362 + t331 * t361 - t299;
t295 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t303;
t284 = -pkin(3) * t302 - pkin(8) * t303;
t333 = qJD(3) ^ 2;
t222 = -pkin(3) * t333 + qJDD(3) * pkin(8) + t284 * t302 + t242;
t320 = t324 ^ 2;
t310 = t329 * g(1) - t332 * g(2);
t353 = qJDD(2) - t310;
t285 = (-pkin(1) - t372) * qJDD(1) + (-qJ(2) + (-t320 - t321) * pkin(7)) * t334 + t353;
t364 = t302 * qJD(3);
t287 = t347 * qJDD(1) + t364;
t225 = (-t287 - t364) * pkin(8) + (-t286 + t299) * pkin(3) + t285;
t327 = sin(qJ(4));
t330 = cos(qJ(4));
t204 = -t327 * t222 + t330 * t225;
t292 = qJD(3) * t330 - t303 * t327;
t255 = qJD(4) * t292 + qJDD(3) * t327 + t287 * t330;
t283 = qJDD(4) - t286;
t293 = qJD(3) * t327 + t303 * t330;
t300 = qJD(4) - t302;
t200 = (t292 * t300 - t255) * pkin(9) + (t292 * t293 + t283) * pkin(4) + t204;
t205 = t330 * t222 + t327 * t225;
t254 = -qJD(4) * t293 + qJDD(3) * t330 - t287 * t327;
t266 = pkin(4) * t300 - pkin(9) * t293;
t291 = t292 ^ 2;
t202 = -pkin(4) * t291 + pkin(9) * t254 - t266 * t300 + t205;
t326 = sin(qJ(5));
t373 = cos(qJ(5));
t196 = t326 * t200 + t373 * t202;
t257 = t326 * t292 + t373 * t293;
t215 = qJD(5) * t257 - t373 * t254 + t255 * t326;
t298 = qJD(5) + t300;
t245 = mrSges(6,1) * t298 - mrSges(6,3) * t257;
t256 = -t373 * t292 + t293 * t326;
t278 = qJDD(5) + t283;
t235 = pkin(5) * t256 - qJ(6) * t257;
t296 = t298 ^ 2;
t191 = -pkin(5) * t296 + qJ(6) * t278 + 0.2e1 * qJD(6) * t298 - t235 * t256 + t196;
t246 = -mrSges(7,1) * t298 + mrSges(7,2) * t257;
t360 = m(7) * t191 + t278 * mrSges(7,3) + t298 * t246;
t236 = mrSges(7,1) * t256 - mrSges(7,3) * t257;
t367 = -mrSges(6,1) * t256 - mrSges(6,2) * t257 - t236;
t371 = -mrSges(6,3) - mrSges(7,2);
t179 = m(6) * t196 - t278 * mrSges(6,2) + t371 * t215 - t298 * t245 + t367 * t256 + t360;
t195 = t373 * t200 - t326 * t202;
t216 = -t256 * qJD(5) + t326 * t254 + t373 * t255;
t244 = -mrSges(6,2) * t298 - mrSges(6,3) * t256;
t193 = -t278 * pkin(5) - t296 * qJ(6) + t257 * t235 + qJDD(6) - t195;
t243 = -mrSges(7,2) * t256 + mrSges(7,3) * t298;
t354 = -m(7) * t193 + t278 * mrSges(7,1) + t298 * t243;
t181 = m(6) * t195 + t278 * mrSges(6,1) + t371 * t216 + t298 * t244 + t367 * t257 + t354;
t174 = t326 * t179 + t373 * t181;
t259 = -mrSges(5,1) * t292 + mrSges(5,2) * t293;
t262 = -mrSges(5,2) * t300 + mrSges(5,3) * t292;
t170 = m(5) * t204 + mrSges(5,1) * t283 - mrSges(5,3) * t255 - t259 * t293 + t262 * t300 + t174;
t263 = mrSges(5,1) * t300 - mrSges(5,3) * t293;
t355 = t373 * t179 - t181 * t326;
t171 = m(5) * t205 - mrSges(5,2) * t283 + mrSges(5,3) * t254 + t259 * t292 - t263 * t300 + t355;
t356 = -t170 * t327 + t330 * t171;
t163 = m(4) * t242 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t286 - qJD(3) * t295 + t279 * t302 + t356;
t241 = t331 * t267 - t328 * t268;
t294 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t302;
t221 = -qJDD(3) * pkin(3) - t333 * pkin(8) + t303 * t284 - t241;
t203 = -t254 * pkin(4) - t291 * pkin(9) + t293 * t266 + t221;
t198 = -0.2e1 * qJD(6) * t257 + (t256 * t298 - t216) * qJ(6) + (t257 * t298 + t215) * pkin(5) + t203;
t188 = m(7) * t198 + t215 * mrSges(7,1) - t216 * mrSges(7,3) + t256 * t243 - t257 * t246;
t340 = m(6) * t203 + t215 * mrSges(6,1) + mrSges(6,2) * t216 + t256 * t244 + t245 * t257 + t188;
t336 = -m(5) * t221 + t254 * mrSges(5,1) - mrSges(5,2) * t255 + t292 * t262 - t263 * t293 - t340;
t176 = m(4) * t241 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t287 + qJD(3) * t294 - t279 * t303 + t336;
t158 = t328 * t163 + t331 * t176;
t288 = -t324 * t304 + t359;
t230 = Ifges(7,1) * t257 + Ifges(7,4) * t298 + Ifges(7,5) * t256;
t231 = Ifges(6,1) * t257 - Ifges(6,4) * t256 + Ifges(6,5) * t298;
t352 = -mrSges(7,1) * t198 + mrSges(7,2) * t191;
t228 = Ifges(7,4) * t257 + Ifges(7,2) * t298 + Ifges(7,6) * t256;
t368 = -Ifges(6,5) * t257 + Ifges(6,6) * t256 - Ifges(6,3) * t298 - t228;
t172 = -mrSges(6,1) * t203 + mrSges(6,3) * t196 - pkin(5) * t188 + (t230 + t231) * t298 + (Ifges(6,6) - Ifges(7,6)) * t278 + t368 * t257 + (Ifges(6,4) - Ifges(7,5)) * t216 + (-Ifges(6,2) - Ifges(7,3)) * t215 + t352;
t229 = Ifges(6,4) * t257 - Ifges(6,2) * t256 + Ifges(6,6) * t298;
t226 = Ifges(7,5) * t257 + Ifges(7,6) * t298 + Ifges(7,3) * t256;
t345 = mrSges(7,2) * t193 - mrSges(7,3) * t198 + Ifges(7,1) * t216 + Ifges(7,4) * t278 + Ifges(7,5) * t215 + t298 * t226;
t173 = mrSges(6,2) * t203 - mrSges(6,3) * t195 + Ifges(6,1) * t216 - Ifges(6,4) * t215 + Ifges(6,5) * t278 - qJ(6) * t188 - t298 * t229 + t368 * t256 + t345;
t248 = Ifges(5,5) * t293 + Ifges(5,6) * t292 + Ifges(5,3) * t300;
t250 = Ifges(5,1) * t293 + Ifges(5,4) * t292 + Ifges(5,5) * t300;
t152 = -mrSges(5,1) * t221 + mrSges(5,3) * t205 + Ifges(5,4) * t255 + Ifges(5,2) * t254 + Ifges(5,6) * t283 - pkin(4) * t340 + pkin(9) * t355 + t373 * t172 + t326 * t173 - t293 * t248 + t300 * t250;
t249 = Ifges(5,4) * t293 + Ifges(5,2) * t292 + Ifges(5,6) * t300;
t154 = mrSges(5,2) * t221 - mrSges(5,3) * t204 + Ifges(5,1) * t255 + Ifges(5,4) * t254 + Ifges(5,5) * t283 - pkin(9) * t174 - t326 * t172 + t373 * t173 + t292 * t248 - t300 * t249;
t270 = Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * qJD(3);
t271 = Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * qJD(3);
t341 = -mrSges(4,1) * t241 + mrSges(4,2) * t242 - Ifges(4,5) * t287 - Ifges(4,6) * t286 - Ifges(4,3) * qJDD(3) - pkin(3) * t336 - pkin(8) * t356 - t330 * t152 - t327 * t154 - t303 * t270 + t302 * t271;
t350 = Ifges(3,4) * t324 + Ifges(3,2) * t325;
t351 = Ifges(3,1) * t324 + Ifges(3,4) * t325;
t374 = -mrSges(3,1) * t288 + mrSges(3,2) * t289 - pkin(2) * t158 - (t350 * t366 - t351 * t365) * qJD(1) + t341;
t370 = mrSges(3,2) * t324;
t165 = t330 * t170 + t327 * t171;
t346 = mrSges(3,3) * qJDD(1) + t334 * (-mrSges(3,1) * t325 + t370);
t156 = m(3) * t288 - t346 * t324 + t158;
t357 = t331 * t163 - t328 * t176;
t157 = m(3) * t289 + t346 * t325 + t357;
t358 = -t156 * t324 + t325 * t157;
t349 = Ifges(3,5) * t324 + Ifges(3,6) * t325;
t269 = Ifges(4,5) * t303 + Ifges(4,6) * t302 + Ifges(4,3) * qJD(3);
t146 = mrSges(4,2) * t285 - mrSges(4,3) * t241 + Ifges(4,1) * t287 + Ifges(4,4) * t286 + Ifges(4,5) * qJDD(3) - pkin(8) * t165 - qJD(3) * t270 - t152 * t327 + t154 * t330 + t269 * t302;
t343 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t216 - Ifges(7,2) * t278 - Ifges(7,6) * t215 + t257 * t226 - t256 * t230;
t337 = mrSges(6,2) * t196 - t256 * t231 - qJ(6) * (-t215 * mrSges(7,2) - t256 * t236 + t360) - pkin(5) * (-t216 * mrSges(7,2) - t257 * t236 + t354) - mrSges(6,1) * t195 + Ifges(6,6) * t215 - Ifges(6,5) * t216 - t257 * t229 - Ifges(6,3) * t278 + t343;
t335 = mrSges(5,1) * t204 - mrSges(5,2) * t205 + Ifges(5,5) * t255 + Ifges(5,6) * t254 + Ifges(5,3) * t283 + pkin(4) * t174 + t293 * t249 - t292 * t250 - t337;
t150 = -mrSges(4,1) * t285 + mrSges(4,3) * t242 + Ifges(4,4) * t287 + Ifges(4,2) * t286 + Ifges(4,6) * qJDD(3) - pkin(3) * t165 + qJD(3) * t271 - t303 * t269 - t335;
t301 = -qJDD(1) * pkin(1) - t334 * qJ(2) + t353;
t306 = t349 * qJD(1);
t342 = m(4) * t285 - t286 * mrSges(4,1) + t287 * mrSges(4,2) - t302 * t294 + t303 * t295 + t165;
t142 = -mrSges(3,1) * t301 + mrSges(3,3) * t289 - pkin(2) * t342 + pkin(7) * t357 + t350 * qJDD(1) + t328 * t146 + t331 * t150 - t306 * t366;
t145 = mrSges(3,2) * t301 - mrSges(3,3) * t288 - pkin(7) * t158 + t351 * qJDD(1) + t331 * t146 - t328 * t150 + t306 * t365;
t339 = -m(3) * t301 + mrSges(3,1) * t361 - t342 + (t320 * t334 + t369) * mrSges(3,3);
t344 = -mrSges(2,2) * t311 + qJ(2) * t358 + t325 * t142 + t324 * t145 + pkin(1) * (-mrSges(3,2) * t362 + t339) + mrSges(2,1) * t310 + Ifges(2,3) * qJDD(1);
t159 = (mrSges(2,1) - t370) * qJDD(1) + t339 - t334 * mrSges(2,2) + m(2) * t310;
t149 = t156 * t325 + t157 * t324;
t147 = m(2) * t311 - mrSges(2,1) * t334 - qJDD(1) * mrSges(2,2) + t358;
t143 = -pkin(1) * t149 + mrSges(2,1) * g(3) + (Ifges(2,6) - t349) * qJDD(1) + t334 * Ifges(2,5) + mrSges(2,3) * t311 + t374;
t140 = -mrSges(2,2) * g(3) - mrSges(2,3) * t310 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t334 - qJ(2) * t149 - t142 * t324 + t145 * t325;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t332 * t140 - t329 * t143 - pkin(6) * (t147 * t329 + t159 * t332), t140, t145, t146, t154, t173, -t228 * t256 + t345; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t329 * t140 + t332 * t143 + pkin(6) * (t147 * t332 - t159 * t329), t143, t142, t150, t152, t172, -t343; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t344, t344, t349 * qJDD(1) - t374, -t341, t335, -t337, Ifges(7,5) * t216 + Ifges(7,6) * t278 + Ifges(7,3) * t215 + t257 * t228 - t298 * t230 - t352;];
m_new  = t1;
