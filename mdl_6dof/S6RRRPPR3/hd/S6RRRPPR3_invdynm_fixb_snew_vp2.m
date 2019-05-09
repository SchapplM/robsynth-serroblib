% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:30:52
% EndTime: 2019-05-07 04:31:10
% DurationCPUTime: 7.79s
% Computational Cost: add. (96550->392), mult. (202250->449), div. (0->0), fcn. (132220->8), ass. (0->146)
t363 = (qJD(2) + qJD(3));
t421 = 2 * qJD(4);
t347 = t363 * t421;
t371 = sin(qJ(2));
t374 = cos(qJ(2));
t402 = qJD(1) * qJD(2);
t338 = qJDD(1) * t371 + t374 * t402;
t376 = qJD(1) ^ 2;
t372 = sin(qJ(1));
t375 = cos(qJ(1));
t345 = -g(1) * t375 - g(2) * t372;
t331 = -pkin(1) * t376 + qJDD(1) * pkin(7) + t345;
t409 = t371 * t331;
t253 = qJDD(2) * pkin(2) - t338 * pkin(8) - t409 + (pkin(2) * t371 * t376 + pkin(8) * t402 - g(3)) * t374;
t311 = -g(3) * t371 + t374 * t331;
t339 = qJDD(1) * t374 - t371 * t402;
t404 = qJD(1) * t371;
t343 = qJD(2) * pkin(2) - pkin(8) * t404;
t410 = t374 ^ 2 * t376;
t254 = -pkin(2) * t410 + pkin(8) * t339 - qJD(2) * t343 + t311;
t370 = sin(qJ(3));
t418 = cos(qJ(3));
t232 = t370 * t253 + t418 * t254;
t403 = qJD(1) * t374;
t328 = t370 * t404 - t418 * t403;
t329 = (t370 * t374 + t418 * t371) * qJD(1);
t300 = pkin(3) * t328 - qJ(4) * t329;
t362 = qJDD(2) + qJDD(3);
t423 = t363 ^ 2;
t397 = pkin(3) * t423 - t362 * qJ(4) + t328 * t300 - t232;
t227 = t347 - t397;
t316 = -(mrSges(5,1) * t363) + mrSges(5,2) * t329;
t426 = m(5) * t227 + t362 * mrSges(5,3) + t363 * t316;
t231 = t418 * t253 - t370 * t254;
t279 = -t328 * qJD(3) + t418 * t338 + t370 * t339;
t299 = mrSges(6,1) * t329 + mrSges(6,2) * t328;
t313 = -(mrSges(4,2) * t363) - mrSges(4,3) * t328;
t394 = -qJ(4) * t423 + t329 * t300 + qJDD(4) - t231;
t229 = -t362 * pkin(3) + t394;
t318 = -mrSges(5,2) * t328 + mrSges(5,3) * t363;
t412 = t328 * t363;
t420 = -2 * qJD(5);
t384 = (-t279 - t412) * qJ(5) + t394 + (t328 * pkin(4) + t420) * t329;
t220 = (-pkin(3) - pkin(4)) * t362 + t384;
t312 = -mrSges(6,1) * t363 - mrSges(6,3) * t328;
t278 = qJD(3) * t329 + t338 * t370 - t418 * t339;
t314 = -pkin(4) * t363 - qJ(5) * t329;
t324 = t328 ^ 2;
t344 = t372 * g(1) - t375 * g(2);
t330 = -qJDD(1) * pkin(1) - t376 * pkin(7) - t344;
t280 = -t339 * pkin(2) - pkin(8) * t410 + t343 * t404 + t330;
t393 = t278 * pkin(3) + t280 + (-t279 + t412) * qJ(4);
t385 = -t324 * qJ(5) + qJDD(5) - t393 + (t314 + t421) * t329;
t419 = -pkin(4) - pkin(9);
t211 = t385 + (-pkin(5) * t328 + (-pkin(3) - pkin(9)) * t329) * t363 + t419 * t278 + t279 * pkin(5);
t303 = pkin(5) * t329 - pkin(9) * t328;
t212 = -t423 * pkin(5) - t329 * t303 + (-pkin(3) + t419) * t362 + t384;
t369 = sin(qJ(6));
t373 = cos(qJ(6));
t207 = t211 * t373 - t212 * t369;
t308 = -t328 * t369 - t363 * t373;
t238 = qJD(6) * t308 + t278 * t373 - t362 * t369;
t309 = t328 * t373 - t363 * t369;
t246 = -mrSges(7,1) * t308 + mrSges(7,2) * t309;
t276 = qJDD(6) + t279;
t323 = qJD(6) + t329;
t281 = -mrSges(7,2) * t323 + mrSges(7,3) * t308;
t203 = m(7) * t207 + mrSges(7,1) * t276 - mrSges(7,3) * t238 - t246 * t309 + t281 * t323;
t208 = t211 * t369 + t212 * t373;
t237 = -qJD(6) * t309 - t278 * t369 - t362 * t373;
t282 = mrSges(7,1) * t323 - mrSges(7,3) * t309;
t204 = m(7) * t208 - mrSges(7,2) * t276 + mrSges(7,3) * t237 + t246 * t308 - t282 * t323;
t408 = -t369 * t203 + t373 * t204;
t398 = -m(6) * t220 - t362 * mrSges(6,2) - t363 * t312 - t408;
t391 = -m(5) * t229 + t362 * mrSges(5,1) + t363 * t318 + t398;
t301 = mrSges(5,1) * t328 - mrSges(5,3) * t329;
t405 = -mrSges(4,1) * t328 - mrSges(4,2) * t329 - t301;
t416 = -mrSges(4,3) - mrSges(5,2);
t179 = m(4) * t231 + t362 * mrSges(4,1) + t363 * t313 + (t299 + t405) * t329 + (mrSges(6,3) + t416) * t279 + t391;
t315 = mrSges(4,1) * t363 - mrSges(4,3) * t329;
t317 = mrSges(6,2) * t363 - mrSges(6,3) * t329;
t389 = t324 * pkin(4) - t278 * qJ(5) + t397;
t422 = -2 * qJD(4);
t222 = t328 * t420 + (t422 - t314) * t363 + t389;
t214 = t362 * pkin(5) - t423 * pkin(9) + t363 * t314 + t347 + ((2 * qJD(5)) + t303) * t328 - t389;
t396 = -m(7) * t214 + t237 * mrSges(7,1) - t238 * mrSges(7,2) + t308 * t281 - t309 * t282;
t388 = -m(6) * t222 + t278 * mrSges(6,3) + t328 * t299 - t396;
t191 = t388 + m(4) * t232 + (-t315 + t317) * t363 + (-mrSges(4,2) + mrSges(6,1)) * t362 + t405 * t328 + t416 * t278 + t426;
t175 = t418 * t179 + t370 * t191;
t310 = -t374 * g(3) - t409;
t326 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t371 + Ifges(3,2) * t374) * qJD(1);
t327 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t371 + Ifges(3,4) * t374) * qJD(1);
t205 = t362 * mrSges(6,1) + t363 * t317 + t388;
t289 = Ifges(4,4) * t329 - Ifges(4,2) * t328 + (Ifges(4,6) * t363);
t292 = Ifges(4,1) * t329 - Ifges(4,4) * t328 + (Ifges(4,5) * t363);
t184 = -t279 * mrSges(6,3) - t329 * t299 - t398;
t285 = Ifges(5,5) * t329 + (Ifges(5,6) * t363) + Ifges(5,3) * t328;
t291 = Ifges(5,1) * t329 + (Ifges(5,4) * t363) + Ifges(5,5) * t328;
t239 = Ifges(7,5) * t309 + Ifges(7,6) * t308 + Ifges(7,3) * t323;
t241 = Ifges(7,1) * t309 + Ifges(7,4) * t308 + Ifges(7,5) * t323;
t195 = -mrSges(7,1) * t214 + mrSges(7,3) * t208 + Ifges(7,4) * t238 + Ifges(7,2) * t237 + Ifges(7,6) * t276 - t239 * t309 + t241 * t323;
t240 = Ifges(7,4) * t309 + Ifges(7,2) * t308 + Ifges(7,6) * t323;
t196 = mrSges(7,2) * t214 - mrSges(7,3) * t207 + Ifges(7,1) * t238 + Ifges(7,4) * t237 + Ifges(7,5) * t276 + t239 * t308 - t240 * t323;
t287 = Ifges(6,4) * t328 - Ifges(6,2) * t329 - (Ifges(6,6) * t363);
t290 = Ifges(6,1) * t328 - Ifges(6,4) * t329 - (Ifges(6,5) * t363);
t390 = mrSges(6,1) * t222 - mrSges(6,2) * t220 + Ifges(6,5) * t278 - Ifges(6,6) * t279 - Ifges(6,3) * t362 + pkin(5) * t396 + pkin(9) * t408 + t373 * t195 + t369 * t196 + t328 * t287 + t329 * t290;
t381 = mrSges(5,1) * t229 - mrSges(5,3) * t227 - Ifges(5,4) * t279 - Ifges(5,2) * t362 - Ifges(5,6) * t278 + pkin(4) * t184 + t329 * t285 - t328 * t291 + t390;
t378 = -mrSges(4,2) * t232 + t328 * t292 + pkin(3) * ((t299 - t301) * t329 + (-mrSges(5,2) + mrSges(6,3)) * t279 + t391) + qJ(4) * (-t278 * mrSges(5,2) - t328 * t301 + t205 + t426) + mrSges(4,1) * t231 + t329 * t289 - Ifges(4,6) * t278 + Ifges(4,5) * t279 + Ifges(4,3) * t362 - t381;
t425 = mrSges(3,1) * t310 - mrSges(3,2) * t311 + Ifges(3,5) * t338 + Ifges(3,6) * t339 + Ifges(3,3) * qJDD(2) + pkin(2) * t175 + (t326 * t371 - t327 * t374) * qJD(1) + t378;
t424 = Ifges(6,4) * t278 - Ifges(6,2) * t279 - Ifges(6,6) * t362 - t363 * t290;
t417 = pkin(3) * t363;
t187 = t373 * t203 + t369 * t204;
t288 = Ifges(5,4) * t329 + Ifges(5,2) * t363 + Ifges(5,6) * t328;
t407 = -Ifges(4,5) * t329 + Ifges(4,6) * t328 - Ifges(4,3) * t363 - t288;
t406 = t287 - t291;
t337 = (-mrSges(3,1) * t374 + mrSges(3,2) * t371) * qJD(1);
t342 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t403;
t173 = m(3) * t310 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t338 + qJD(2) * t342 - t337 * t404 + t175;
t341 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t404;
t400 = -t179 * t370 + t418 * t191;
t174 = m(3) * t311 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t339 - qJD(2) * t341 + t337 * t403 + t400;
t401 = -t173 * t371 + t374 * t174;
t218 = -t278 * pkin(4) - t329 * t417 + t385;
t183 = -m(6) * t218 - t279 * mrSges(6,1) - t278 * mrSges(6,2) - t328 * t312 - t329 * t317 - t187;
t224 = (t422 + t417) * t329 + t393;
t180 = m(5) * t224 + t278 * mrSges(5,1) - t279 * mrSges(5,3) - t329 * t316 + t328 * t318 + t183;
t284 = Ifges(6,5) * t328 - Ifges(6,6) * t329 - Ifges(6,3) * t363;
t387 = -mrSges(6,2) * t218 + mrSges(6,3) * t222 - Ifges(6,1) * t278 + Ifges(6,4) * t279 + Ifges(6,5) * t362 + pkin(9) * t187 + t369 * t195 - t373 * t196 + t329 * t284;
t382 = mrSges(5,1) * t224 - mrSges(5,2) * t227 + pkin(4) * t183 + qJ(5) * t205 - t387;
t167 = -t382 + mrSges(4,3) * t232 - pkin(3) * t180 + (t292 - t406) * t363 + (Ifges(4,6) - Ifges(5,6)) * t362 + t407 * t329 + (-Ifges(5,5) + Ifges(4,4)) * t279 + (-Ifges(5,3) - Ifges(4,2)) * t278 - mrSges(4,1) * t280;
t392 = mrSges(7,1) * t207 - mrSges(7,2) * t208 + Ifges(7,5) * t238 + Ifges(7,6) * t237 + Ifges(7,3) * t276 + t309 * t240 - t308 * t241;
t386 = mrSges(6,1) * t218 - mrSges(6,3) * t220 + pkin(5) * t187 + t328 * t284 + t392;
t380 = mrSges(5,2) * t229 - mrSges(5,3) * t224 + Ifges(5,1) * t279 + Ifges(5,4) * t362 + Ifges(5,5) * t278 - qJ(5) * t184 + t363 * t285 + t386;
t168 = t380 - mrSges(4,3) * t231 - qJ(4) * t180 + (-t289 + t290) * t363 + (Ifges(4,5) + Ifges(6,6)) * t362 + t407 * t328 + (Ifges(4,1) + Ifges(6,2)) * t279 + (-Ifges(6,4) - Ifges(4,4)) * t278 + mrSges(4,2) * t280;
t325 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t371 + Ifges(3,6) * t374) * qJD(1);
t383 = m(4) * t280 + t278 * mrSges(4,1) + t279 * mrSges(4,2) + t328 * t313 + t329 * t315 + t180;
t163 = -mrSges(3,1) * t330 + mrSges(3,3) * t311 + Ifges(3,4) * t338 + Ifges(3,2) * t339 + Ifges(3,6) * qJDD(2) - pkin(2) * t383 + pkin(8) * t400 + qJD(2) * t327 + t418 * t167 + t370 * t168 - t325 * t404;
t165 = mrSges(3,2) * t330 - mrSges(3,3) * t310 + Ifges(3,1) * t338 + Ifges(3,4) * t339 + Ifges(3,5) * qJDD(2) - pkin(8) * t175 - qJD(2) * t326 - t370 * t167 + t418 * t168 + t325 * t403;
t379 = -m(3) * t330 + t339 * mrSges(3,1) - t338 * mrSges(3,2) - t341 * t404 + t342 * t403 - t383;
t395 = mrSges(2,1) * t344 - mrSges(2,2) * t345 + Ifges(2,3) * qJDD(1) + pkin(1) * t379 + pkin(7) * t401 + t374 * t163 + t371 * t165;
t176 = m(2) * t344 + qJDD(1) * mrSges(2,1) - t376 * mrSges(2,2) + t379;
t171 = t173 * t374 + t174 * t371;
t169 = m(2) * t345 - mrSges(2,1) * t376 - qJDD(1) * mrSges(2,2) + t401;
t166 = mrSges(2,1) * g(3) + mrSges(2,3) * t345 + t376 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t171 - t425;
t161 = -mrSges(2,2) * g(3) - mrSges(2,3) * t344 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t376 - pkin(7) * t171 - t163 * t371 + t165 * t374;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t375 * t161 - t372 * t166 - pkin(6) * (t169 * t372 + t176 * t375), t161, t165, t168, -t328 * t288 + t380 - t424, t363 * t287 - t387, t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t372 * t161 + t375 * t166 + pkin(6) * (t169 * t375 - t176 * t372), t166, t163, t167, -t381, -t386 + t424, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t395, t395, t425, t378, Ifges(5,5) * t279 + Ifges(5,6) * t362 + Ifges(5,3) * t278 + t329 * t288 + t406 * t363 + t382, t390, t392;];
m_new  = t1;
