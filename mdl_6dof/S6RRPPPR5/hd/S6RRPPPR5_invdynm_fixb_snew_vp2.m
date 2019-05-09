% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-05-06 08:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:51:46
% EndTime: 2019-05-06 08:52:00
% DurationCPUTime: 6.77s
% Computational Cost: add. (76365->394), mult. (168385->447), div. (0->0), fcn. (103800->8), ass. (0->149)
t349 = sin(pkin(9));
t351 = sin(qJ(2));
t390 = qJD(1) * t351;
t400 = cos(pkin(9));
t322 = t349 * qJD(2) + t390 * t400;
t411 = -0.2e1 * t322;
t410 = 2 * qJD(4);
t352 = sin(qJ(1));
t355 = cos(qJ(1));
t338 = -g(1) * t355 - g(2) * t352;
t357 = qJD(1) ^ 2;
t315 = -pkin(1) * t357 + qJDD(1) * pkin(7) + t338;
t354 = cos(qJ(2));
t277 = -t354 * g(3) - t351 * t315;
t330 = (-pkin(2) * t354 - qJ(3) * t351) * qJD(1);
t356 = qJD(2) ^ 2;
t370 = qJDD(2) * pkin(2) + t356 * qJ(3) - t330 * t390 - qJDD(3) + t277;
t321 = -qJD(2) * t400 + t349 * t390;
t389 = qJD(1) * t354;
t385 = t321 * t389;
t409 = -qJ(4) * t385 - t370;
t273 = pkin(3) * t321 - qJ(4) * t322;
t386 = qJD(1) * qJD(2);
t342 = t351 * t386;
t333 = qJDD(1) * t354 - t342;
t337 = t352 * g(1) - t355 * g(2);
t314 = -qJDD(1) * pkin(1) - t357 * pkin(7) - t337;
t382 = t354 * t386;
t332 = qJDD(1) * t351 + t382;
t240 = (-t332 - t382) * qJ(3) + (-t333 + t342) * pkin(2) + t314;
t278 = -g(3) * t351 + t354 * t315;
t248 = -pkin(2) * t356 + qJDD(2) * qJ(3) + t330 * t389 + t278;
t395 = t349 * t240 + t400 * t248;
t397 = t354 ^ 2 * t357;
t408 = pkin(3) * t397 + t333 * qJ(4) + t321 * t273 + t389 * t410 - t395;
t223 = qJD(3) * t411 + t240 * t400 - t349 * t248;
t388 = qJD(3) * t321;
t309 = -0.2e1 * t388;
t224 = t309 + t395;
t263 = Ifges(4,4) * t322 - Ifges(4,2) * t321 - Ifges(4,6) * t389;
t275 = -mrSges(5,2) * t321 - mrSges(5,3) * t322;
t295 = -mrSges(6,2) * t322 + mrSges(6,3) * t389;
t297 = mrSges(5,1) * t322 - mrSges(5,2) * t389;
t301 = -qJDD(2) * t400 + t332 * t349;
t302 = t349 * qJDD(2) + t332 * t400;
t221 = t333 * pkin(3) - qJ(4) * t397 + t322 * t273 + qJDD(4) - t223;
t398 = t321 * t322;
t405 = 2 * qJD(5);
t211 = t389 * t405 + (t333 + t398) * qJ(5) + (t302 - t385) * pkin(4) + t221;
t272 = mrSges(6,1) * t322 - mrSges(6,3) * t321;
t303 = -pkin(5) * t389 - pkin(8) * t321;
t319 = t322 ^ 2;
t206 = -t319 * pkin(5) + t302 * pkin(8) + t303 * t389 + t211;
t294 = pkin(4) * t322 + qJ(5) * t389;
t318 = t321 ^ 2;
t215 = -t301 * pkin(4) - t318 * qJ(5) - t294 * t389 + qJDD(5) + t309 - t408;
t384 = t322 * t389;
t207 = t215 + (-t333 + t398) * pkin(5) + (-t301 - t384) * pkin(8);
t350 = sin(qJ(6));
t353 = cos(qJ(6));
t202 = -t206 * t350 + t207 * t353;
t270 = -t321 * t350 + t322 * t353;
t232 = qJD(6) * t270 + t301 * t353 + t302 * t350;
t271 = t321 * t353 + t322 * t350;
t238 = -mrSges(7,1) * t270 + mrSges(7,2) * t271;
t339 = qJD(6) - t389;
t249 = -mrSges(7,2) * t339 + mrSges(7,3) * t270;
t329 = qJDD(6) - t333;
t198 = m(7) * t202 + mrSges(7,1) * t329 - mrSges(7,3) * t232 - t238 * t271 + t249 * t339;
t203 = t206 * t353 + t207 * t350;
t231 = -qJD(6) * t271 - t301 * t350 + t302 * t353;
t250 = mrSges(7,1) * t339 - mrSges(7,3) * t271;
t199 = m(7) * t203 - mrSges(7,2) * t329 + mrSges(7,3) * t231 + t238 * t270 - t250 * t339;
t396 = -t350 * t198 + t353 * t199;
t375 = -m(6) * t211 + t302 * mrSges(6,2) + t322 * t272 - t396;
t298 = -mrSges(6,1) * t389 + mrSges(6,2) * t321;
t383 = t298 * t389;
t179 = mrSges(6,3) * t333 - t375 + t383;
t219 = 0.2e1 * t388 + t408;
t257 = -Ifges(5,5) * t389 - Ifges(5,6) * t322 + Ifges(5,3) * t321;
t183 = t353 * t198 + t350 * t199;
t260 = Ifges(6,5) * t321 + Ifges(6,6) * t389 + Ifges(6,3) * t322;
t264 = Ifges(6,1) * t321 + Ifges(6,4) * t389 + Ifges(6,5) * t322;
t234 = Ifges(7,4) * t271 + Ifges(7,2) * t270 + Ifges(7,6) * t339;
t235 = Ifges(7,1) * t271 + Ifges(7,4) * t270 + Ifges(7,5) * t339;
t374 = mrSges(7,1) * t202 - mrSges(7,2) * t203 + Ifges(7,5) * t232 + Ifges(7,6) * t231 + Ifges(7,3) * t329 + t271 * t234 - t270 * t235;
t365 = -mrSges(6,1) * t215 + mrSges(6,3) * t211 + Ifges(6,4) * t301 + Ifges(6,2) * t333 + Ifges(6,6) * t302 - pkin(5) * t183 - t321 * t260 + t322 * t264 - t374;
t359 = -mrSges(5,2) * t221 + mrSges(5,3) * t219 + Ifges(5,1) * t333 + Ifges(5,4) * t302 - Ifges(5,5) * t301 + qJ(5) * t179 + t322 * t257 + t365;
t296 = mrSges(5,1) * t321 + mrSges(5,3) * t389;
t368 = -m(5) * t221 - t302 * mrSges(5,1) - t322 * t275 + t296 * t389 + t375;
t378 = m(6) * t215 + t301 * mrSges(6,2) + t321 * t272 + t183;
t372 = -m(5) * t219 - t333 * mrSges(5,3) + t378;
t258 = -Ifges(5,4) * t389 - Ifges(5,2) * t322 + Ifges(5,6) * t321;
t392 = Ifges(4,1) * t322 - Ifges(4,4) * t321 - Ifges(4,5) * t389 - t258;
t401 = t333 * mrSges(6,1);
t403 = mrSges(5,2) - mrSges(6,3);
t407 = -t392 * t321 - mrSges(4,1) * t223 + mrSges(4,2) * t224 - Ifges(4,5) * t302 + Ifges(4,6) * t301 - pkin(3) * (t333 * t403 + t368 - t383) - qJ(4) * (-t301 * mrSges(5,1) - t401 - t321 * t275 + (-t295 - t297) * t389 + t372) - t322 * t263 + t359;
t399 = t302 * qJ(4);
t404 = pkin(4) * t318;
t218 = t404 + t399 - 0.2e1 * qJD(5) * t321 + (-pkin(3) - qJ(5)) * t301 + (pkin(3) * t389 + t294 + t410) * t322 - t409;
t363 = qJD(4) * t411 + (t301 - t384) * pkin(3) + t409;
t209 = (-pkin(5) - qJ(4)) * t302 + (t405 + t303) * t321 - pkin(8) * t319 + t363 - t404 + qJ(5) * t301 - t294 * t322;
t377 = -m(7) * t209 + t231 * mrSges(7,1) - t232 * mrSges(7,2) + t270 * t249 - t271 * t250;
t193 = -m(6) * t218 - t302 * mrSges(6,1) + t301 * mrSges(6,3) - t322 * t295 + t321 * t298 - t377;
t222 = t363 - t399;
t366 = -m(5) * t222 + t301 * mrSges(5,2) + t321 * t296 - t193;
t192 = -t302 * mrSges(5,3) - t322 * t297 - t366;
t262 = Ifges(6,4) * t321 + Ifges(6,2) * t389 + Ifges(6,6) * t322;
t233 = Ifges(7,5) * t271 + Ifges(7,6) * t270 + Ifges(7,3) * t339;
t190 = -mrSges(7,1) * t209 + mrSges(7,3) * t203 + Ifges(7,4) * t232 + Ifges(7,2) * t231 + Ifges(7,6) * t329 - t233 * t271 + t235 * t339;
t191 = mrSges(7,2) * t209 - mrSges(7,3) * t202 + Ifges(7,1) * t232 + Ifges(7,4) * t231 + Ifges(7,5) * t329 + t233 * t270 - t234 * t339;
t367 = -mrSges(6,2) * t215 + mrSges(6,3) * t218 - Ifges(6,1) * t301 - Ifges(6,4) * t333 - Ifges(6,5) * t302 + pkin(8) * t183 + t350 * t190 - t353 * t191 - t260 * t389;
t361 = mrSges(5,1) * t219 - mrSges(5,2) * t222 + pkin(4) * (t295 * t389 - t378 + t401) + qJ(5) * t193 - t367;
t259 = -Ifges(5,1) * t389 - Ifges(5,4) * t322 + Ifges(5,5) * t321;
t393 = -Ifges(4,5) * t322 + Ifges(4,6) * t321 + Ifges(4,3) * t389 - t259;
t402 = -Ifges(5,6) - Ifges(4,4);
t163 = -t361 + (-Ifges(4,6) + Ifges(5,5)) * t333 + (t262 + t393) * t322 - t402 * t302 + (-Ifges(5,3) - Ifges(4,2)) * t301 - t392 * t389 + mrSges(4,1) * t370 + mrSges(4,3) * t224 - pkin(3) * t192;
t369 = -mrSges(6,1) * t218 + mrSges(6,2) * t211 - Ifges(6,5) * t301 - Ifges(6,6) * t333 - Ifges(6,3) * t302 - pkin(5) * t377 - pkin(8) * t396 - t353 * t190 - t350 * t191 - t321 * t262;
t362 = -mrSges(5,1) * t221 + mrSges(5,3) * t222 - pkin(4) * t179 + t369;
t394 = t257 + t264;
t164 = -t362 + (-Ifges(4,5) + Ifges(5,4)) * t333 + t393 * t321 + (Ifges(4,1) + Ifges(5,2)) * t302 + t402 * t301 - mrSges(4,2) * t370 - mrSges(4,3) * t223 - qJ(4) * t192 + (t263 - t394) * t389;
t274 = mrSges(4,1) * t321 + mrSges(4,2) * t322;
t299 = mrSges(4,2) * t389 - mrSges(4,3) * t321;
t174 = m(4) * t223 - mrSges(4,3) * t302 - t274 * t322 + (-t298 - t299) * t389 + (-mrSges(4,1) + t403) * t333 + t368;
t391 = mrSges(4,1) * t389 + mrSges(4,3) * t322 + t297;
t176 = m(4) * t224 + (-mrSges(6,1) + mrSges(4,2)) * t333 + (-t274 - t275) * t321 + (-mrSges(4,3) - mrSges(5,1)) * t301 + (-t295 - t391) * t389 + t372;
t173 = -t174 * t349 + t400 * t176;
t186 = -t321 * t299 + t366 + t391 * t322 + (-mrSges(4,2) + mrSges(5,3)) * t302 + m(4) * t370 - t301 * mrSges(4,1);
t312 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t351 + Ifges(3,2) * t354) * qJD(1);
t313 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t351 + Ifges(3,4) * t354) * qJD(1);
t406 = mrSges(3,1) * t277 - mrSges(3,2) * t278 + Ifges(3,5) * t332 + Ifges(3,6) * t333 + Ifges(3,3) * qJDD(2) + pkin(2) * t186 + qJ(3) * t173 + (t312 * t351 - t313 * t354) * qJD(1) + t163 * t400 + t349 * t164;
t331 = (-mrSges(3,1) * t354 + mrSges(3,2) * t351) * qJD(1);
t335 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t390;
t171 = m(3) * t278 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t333 - qJD(2) * t335 + t331 * t389 + t173;
t336 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t389;
t185 = m(3) * t277 + qJDD(2) * mrSges(3,1) - t332 * mrSges(3,3) + qJD(2) * t336 - t331 * t390 + t186;
t380 = t354 * t171 - t185 * t351;
t172 = t174 * t400 + t349 * t176;
t311 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t351 + Ifges(3,6) * t354) * qJD(1);
t160 = mrSges(3,2) * t314 - mrSges(3,3) * t277 + Ifges(3,1) * t332 + Ifges(3,4) * t333 + Ifges(3,5) * qJDD(2) - qJ(3) * t172 - qJD(2) * t312 - t349 * t163 + t164 * t400 + t311 * t389;
t162 = (Ifges(4,3) + Ifges(3,2)) * t333 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t332 - mrSges(3,1) * t314 + qJD(2) * t313 + mrSges(3,3) * t278 - pkin(2) * t172 - t311 * t390 + t407;
t364 = -m(3) * t314 + t333 * mrSges(3,1) - t332 * mrSges(3,2) - t335 * t390 + t336 * t389 - t172;
t371 = mrSges(2,1) * t337 - mrSges(2,2) * t338 + Ifges(2,3) * qJDD(1) + pkin(1) * t364 + pkin(7) * t380 + t351 * t160 + t354 * t162;
t168 = m(2) * t337 + qJDD(1) * mrSges(2,1) - t357 * mrSges(2,2) + t364;
t167 = t171 * t351 + t185 * t354;
t165 = m(2) * t338 - mrSges(2,1) * t357 - qJDD(1) * mrSges(2,2) + t380;
t158 = mrSges(2,1) * g(3) + mrSges(2,3) * t338 + t357 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t167 - t406;
t157 = -mrSges(2,2) * g(3) - mrSges(2,3) * t337 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t357 - pkin(7) * t167 + t160 * t354 - t162 * t351;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t355 * t157 - t352 * t158 - pkin(6) * (t165 * t352 + t168 * t355), t157, t160, t164, -t321 * t258 - t359, -t322 * t262 - t367, t191; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t352 * t157 + t355 * t158 + pkin(6) * (t165 * t355 - t168 * t352), t158, t162, t163, -Ifges(5,4) * t333 - Ifges(5,2) * t302 + Ifges(5,6) * t301 + t321 * t259 + t389 * t394 + t362, t365, t190; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t371, t371, t406, -Ifges(4,3) * t333 - t407, t361 + (t259 - t262) * t322 - Ifges(5,5) * t333 - Ifges(5,6) * t302 + Ifges(5,3) * t301 - t258 * t389, -t264 * t389 - t369, t374;];
m_new  = t1;
