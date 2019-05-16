% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:52:09
% EndTime: 2019-05-05 18:52:30
% DurationCPUTime: 12.58s
% Computational Cost: add. (200784->359), mult. (486441->439), div. (0->0), fcn. (360580->10), ass. (0->146)
t350 = sin(pkin(10));
t360 = qJD(1) ^ 2;
t351 = cos(pkin(10));
t394 = t351 ^ 2 * t360;
t402 = t350 ^ 2 * t360 + t394;
t391 = qJD(1) * t350;
t401 = qJD(1) * t351;
t355 = sin(qJ(1));
t358 = cos(qJ(1));
t321 = -t358 * g(1) - t355 * g(2);
t314 = -t360 * pkin(1) + qJDD(1) * qJ(2) + t321;
t388 = qJD(1) * qJD(2);
t384 = -t351 * g(3) - 0.2e1 * t350 * t388;
t265 = (pkin(2) * t351 * t360 - pkin(7) * qJDD(1) - t314) * t350 + t384;
t299 = -t350 * g(3) + (t314 + 0.2e1 * t388) * t351;
t386 = qJDD(1) * t351;
t266 = -pkin(2) * t394 + pkin(7) * t386 + t299;
t354 = sin(qJ(3));
t398 = cos(qJ(3));
t244 = t354 * t265 + t398 * t266;
t385 = t351 * t398;
t387 = qJDD(1) * t350;
t375 = t398 * t350 + t351 * t354;
t313 = t375 * qJD(1);
t389 = t313 * qJD(3);
t296 = -qJDD(1) * t385 + t354 * t387 + t389;
t303 = qJD(3) * mrSges(4,1) - t313 * mrSges(4,3);
t312 = -qJD(1) * t385 + t354 * t391;
t243 = t398 * t265 - t354 * t266;
t282 = t312 * pkin(3) - t313 * qJ(4);
t359 = qJD(3) ^ 2;
t227 = -qJDD(3) * pkin(3) - t359 * qJ(4) + t313 * t282 + qJDD(4) - t243;
t390 = t312 * qJD(3);
t297 = t375 * qJDD(1) - t390;
t215 = (-t297 - t390) * pkin(8) + (t312 * t313 - qJDD(3)) * pkin(4) + t227;
t399 = 2 * qJD(4);
t225 = -t359 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t399 - t312 * t282 + t244;
t306 = -qJD(3) * pkin(4) - t313 * pkin(8);
t311 = t312 ^ 2;
t217 = -t311 * pkin(4) + t296 * pkin(8) + qJD(3) * t306 + t225;
t353 = sin(qJ(5));
t357 = cos(qJ(5));
t211 = t353 * t215 + t357 * t217;
t276 = t353 * t312 + t357 * t313;
t239 = -t276 * qJD(5) + t357 * t296 - t353 * t297;
t275 = t357 * t312 - t353 * t313;
t252 = -t275 * mrSges(6,1) + t276 * mrSges(6,2);
t340 = -qJD(3) + qJD(5);
t262 = t340 * mrSges(6,1) - t276 * mrSges(6,3);
t337 = -qJDD(3) + qJDD(5);
t253 = -t275 * pkin(5) - t276 * pkin(9);
t336 = t340 ^ 2;
t206 = -t336 * pkin(5) + t337 * pkin(9) + t275 * t253 + t211;
t320 = t355 * g(1) - t358 * g(2);
t310 = -qJDD(1) * pkin(1) - t360 * qJ(2) + qJDD(2) - t320;
t295 = -pkin(2) * t386 - pkin(7) * t402 + t310;
t369 = t296 * pkin(3) + t295 + (-t297 + t390) * qJ(4);
t213 = -pkin(3) * t389 - t296 * pkin(4) - t311 * pkin(8) - t369 + (t306 + t399) * t313;
t240 = t275 * qJD(5) + t353 * t296 + t357 * t297;
t207 = t213 + (t276 * t340 - t239) * pkin(5) + (-t275 * t340 - t240) * pkin(9);
t352 = sin(qJ(6));
t356 = cos(qJ(6));
t203 = -t352 * t206 + t356 * t207;
t258 = -t352 * t276 + t356 * t340;
t220 = t258 * qJD(6) + t356 * t240 + t352 * t337;
t238 = qJDD(6) - t239;
t259 = t356 * t276 + t352 * t340;
t241 = -t258 * mrSges(7,1) + t259 * mrSges(7,2);
t268 = qJD(6) - t275;
t245 = -t268 * mrSges(7,2) + t258 * mrSges(7,3);
t199 = m(7) * t203 + t238 * mrSges(7,1) - t220 * mrSges(7,3) - t259 * t241 + t268 * t245;
t204 = t356 * t206 + t352 * t207;
t219 = -t259 * qJD(6) - t352 * t240 + t356 * t337;
t246 = t268 * mrSges(7,1) - t259 * mrSges(7,3);
t200 = m(7) * t204 - t238 * mrSges(7,2) + t219 * mrSges(7,3) + t258 * t241 - t268 * t246;
t381 = -t352 * t199 + t356 * t200;
t186 = m(6) * t211 - t337 * mrSges(6,2) + t239 * mrSges(6,3) + t275 * t252 - t340 * t262 + t381;
t210 = t357 * t215 - t353 * t217;
t261 = -t340 * mrSges(6,2) + t275 * mrSges(6,3);
t205 = -t337 * pkin(5) - t336 * pkin(9) + t276 * t253 - t210;
t372 = -m(7) * t205 + t219 * mrSges(7,1) - t220 * mrSges(7,2) + t258 * t245 - t259 * t246;
t195 = m(6) * t210 + t337 * mrSges(6,1) - t240 * mrSges(6,3) - t276 * t252 + t340 * t261 + t372;
t181 = t357 * t186 - t353 * t195;
t304 = -qJD(3) * mrSges(5,1) + t313 * mrSges(5,2);
t374 = m(5) * t225 + qJDD(3) * mrSges(5,3) + qJD(3) * t304 + t181;
t283 = t312 * mrSges(5,1) - t313 * mrSges(5,3);
t392 = -t312 * mrSges(4,1) - t313 * mrSges(4,2) - t283;
t397 = -mrSges(4,3) - mrSges(5,2);
t174 = m(4) * t244 - qJDD(3) * mrSges(4,2) - qJD(3) * t303 + t397 * t296 + t392 * t312 + t374;
t302 = -qJD(3) * mrSges(4,2) - t312 * mrSges(4,3);
t180 = t353 * t186 + t357 * t195;
t305 = -t312 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t371 = -m(5) * t227 + qJDD(3) * mrSges(5,1) + qJD(3) * t305 - t180;
t175 = m(4) * t243 + qJDD(3) * mrSges(4,1) + qJD(3) * t302 + t397 * t297 + t392 * t313 + t371;
t167 = t354 * t174 + t398 * t175;
t298 = -t350 * t314 + t384;
t272 = Ifges(4,4) * t313 - Ifges(4,2) * t312 + Ifges(4,6) * qJD(3);
t274 = Ifges(4,1) * t313 - Ifges(4,4) * t312 + Ifges(4,5) * qJD(3);
t269 = Ifges(5,5) * t313 + Ifges(5,6) * qJD(3) + Ifges(5,3) * t312;
t273 = Ifges(5,1) * t313 + Ifges(5,4) * qJD(3) + Ifges(5,5) * t312;
t228 = Ifges(7,5) * t259 + Ifges(7,6) * t258 + Ifges(7,3) * t268;
t230 = Ifges(7,1) * t259 + Ifges(7,4) * t258 + Ifges(7,5) * t268;
t193 = -mrSges(7,1) * t205 + mrSges(7,3) * t204 + Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t238 - t259 * t228 + t268 * t230;
t229 = Ifges(7,4) * t259 + Ifges(7,2) * t258 + Ifges(7,6) * t268;
t194 = mrSges(7,2) * t205 - mrSges(7,3) * t203 + Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t238 + t258 * t228 - t268 * t229;
t248 = Ifges(6,4) * t276 + Ifges(6,2) * t275 + Ifges(6,6) * t340;
t249 = Ifges(6,1) * t276 + Ifges(6,4) * t275 + Ifges(6,5) * t340;
t370 = mrSges(6,1) * t210 - mrSges(6,2) * t211 + Ifges(6,5) * t240 + Ifges(6,6) * t239 + Ifges(6,3) * t337 + pkin(5) * t372 + pkin(9) * t381 + t356 * t193 + t352 * t194 + t276 * t248 - t275 * t249;
t364 = mrSges(5,1) * t227 - mrSges(5,3) * t225 - Ifges(5,4) * t297 - Ifges(5,2) * qJDD(3) - Ifges(5,6) * t296 + pkin(4) * t180 + t313 * t269 - t312 * t273 + t370;
t362 = -mrSges(4,2) * t244 + t312 * t274 + qJ(4) * (-t296 * mrSges(5,2) - t312 * t283 + t374) + pkin(3) * (-t297 * mrSges(5,2) - t313 * t283 + t371) + mrSges(4,1) * t243 + t313 * t272 - Ifges(4,6) * t296 + Ifges(4,5) * t297 + Ifges(4,3) * qJDD(3) - t364;
t379 = Ifges(3,4) * t350 + Ifges(3,2) * t351;
t380 = Ifges(3,1) * t350 + Ifges(3,4) * t351;
t400 = -mrSges(3,1) * t298 + mrSges(3,2) * t299 - pkin(2) * t167 - (t379 * t391 - t380 * t401) * qJD(1) - t362;
t396 = t350 * mrSges(3,2);
t189 = t356 * t199 + t352 * t200;
t271 = Ifges(5,4) * t313 + Ifges(5,2) * qJD(3) + Ifges(5,6) * t312;
t393 = -Ifges(4,5) * t313 + Ifges(4,6) * t312 - Ifges(4,3) * qJD(3) - t271;
t376 = mrSges(3,3) * qJDD(1) + t360 * (-mrSges(3,1) * t351 + t396);
t165 = m(3) * t298 - t376 * t350 + t167;
t382 = t398 * t174 - t354 * t175;
t166 = m(3) * t299 + t376 * t351 + t382;
t383 = -t350 * t165 + t351 * t166;
t378 = Ifges(3,5) * t350 + Ifges(3,6) * t351;
t187 = -m(6) * t213 + t239 * mrSges(6,1) - t240 * mrSges(6,2) + t275 * t261 - t276 * t262 - t189;
t222 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t313 + t369;
t184 = m(5) * t222 + t296 * mrSges(5,1) - t297 * mrSges(5,3) - t313 * t304 + t312 * t305 + t187;
t247 = Ifges(6,5) * t276 + Ifges(6,6) * t275 + Ifges(6,3) * t340;
t169 = mrSges(6,2) * t213 - mrSges(6,3) * t210 + Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t337 - pkin(9) * t189 - t352 * t193 + t356 * t194 + t275 * t247 - t340 * t248;
t366 = mrSges(7,1) * t203 - mrSges(7,2) * t204 + Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t238 + t259 * t229 - t258 * t230;
t170 = -mrSges(6,1) * t213 + mrSges(6,3) * t211 + Ifges(6,4) * t240 + Ifges(6,2) * t239 + Ifges(6,6) * t337 - pkin(5) * t189 - t276 * t247 + t340 * t249 - t366;
t367 = -mrSges(5,1) * t222 + mrSges(5,2) * t225 - pkin(4) * t187 - pkin(8) * t181 - t353 * t169 - t357 * t170;
t159 = -mrSges(4,1) * t295 + mrSges(4,3) * t244 - pkin(3) * t184 + t393 * t313 + (Ifges(4,4) - Ifges(5,5)) * t297 + (-Ifges(4,2) - Ifges(5,3)) * t296 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + (t274 + t273) * qJD(3) + t367;
t368 = mrSges(5,2) * t227 - mrSges(5,3) * t222 + Ifges(5,1) * t297 + Ifges(5,4) * qJDD(3) + Ifges(5,5) * t296 - pkin(8) * t180 + qJD(3) * t269 + t357 * t169 - t353 * t170;
t160 = mrSges(4,2) * t295 - mrSges(4,3) * t243 + Ifges(4,1) * t297 - Ifges(4,4) * t296 + Ifges(4,5) * qJDD(3) - qJ(4) * t184 - qJD(3) * t272 + t393 * t312 + t368;
t316 = t378 * qJD(1);
t365 = m(4) * t295 + t296 * mrSges(4,1) + t297 * mrSges(4,2) + t312 * t302 + t313 * t303 + t184;
t155 = -mrSges(3,1) * t310 + mrSges(3,3) * t299 - pkin(2) * t365 + pkin(7) * t382 + t379 * qJDD(1) + t398 * t159 + t354 * t160 - t316 * t391;
t157 = mrSges(3,2) * t310 - mrSges(3,3) * t298 - pkin(7) * t167 + t380 * qJDD(1) - t354 * t159 + t398 * t160 + t316 * t401;
t363 = -m(3) * t310 + mrSges(3,1) * t386 + mrSges(3,3) * t402 - t365;
t373 = -mrSges(2,2) * t321 + qJ(2) * t383 + t351 * t155 + t350 * t157 + pkin(1) * (-mrSges(3,2) * t387 + t363) + mrSges(2,1) * t320 + Ifges(2,3) * qJDD(1);
t182 = (mrSges(2,1) - t396) * qJDD(1) - t360 * mrSges(2,2) + m(2) * t320 + t363;
t163 = t351 * t165 + t350 * t166;
t161 = m(2) * t321 - t360 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t383;
t158 = t360 * Ifges(2,5) + mrSges(2,3) * t321 + mrSges(2,1) * g(3) - pkin(1) * t163 + (Ifges(2,6) - t378) * qJDD(1) + t400;
t153 = -mrSges(2,2) * g(3) - mrSges(2,3) * t320 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t360 - qJ(2) * t163 - t155 * t350 + t157 * t351;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t358 * t153 - t355 * t158 - pkin(6) * (t161 * t355 + t182 * t358), t153, t157, t160, -t271 * t312 + t368, t169, t194; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t355 * t153 + t358 * t158 + pkin(6) * (t161 * t358 - t182 * t355), t158, t155, t159, -t364, t170, t193; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t373, t373, t378 * qJDD(1) - t400, t362, Ifges(5,5) * t297 + Ifges(5,6) * qJDD(3) + Ifges(5,3) * t296 - qJD(3) * t273 + t313 * t271 - t367, t370, t366;];
m_new  = t1;
