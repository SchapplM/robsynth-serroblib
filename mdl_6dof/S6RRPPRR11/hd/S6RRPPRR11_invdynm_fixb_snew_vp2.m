% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR11
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 12:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:59:07
% EndTime: 2019-05-06 11:59:53
% DurationCPUTime: 27.12s
% Computational Cost: add. (453554->405), mult. (1070524->507), div. (0->0), fcn. (794202->12), ass. (0->162)
t401 = -2 * qJD(3);
t353 = cos(pkin(6));
t344 = t353 * qJD(1) + qJD(2);
t356 = sin(qJ(2));
t351 = sin(pkin(6));
t386 = qJD(1) * t351;
t379 = t356 * t386;
t400 = (pkin(2) * t344 + t401) * t379;
t357 = sin(qJ(1));
t361 = cos(qJ(1));
t339 = t357 * g(1) - t361 * g(2);
t362 = qJD(1) ^ 2;
t399 = pkin(8) * t351;
t323 = qJDD(1) * pkin(1) + t362 * t399 + t339;
t340 = -t361 * g(1) - t357 * g(2);
t383 = qJDD(1) * t351;
t324 = -t362 * pkin(1) + pkin(8) * t383 + t340;
t360 = cos(qJ(2));
t391 = t353 * t356;
t393 = t351 * t356;
t276 = -g(3) * t393 + t323 * t391 + t360 * t324;
t325 = (-pkin(2) * t360 - qJ(3) * t356) * t386;
t342 = t344 ^ 2;
t343 = t353 * qJDD(1) + qJDD(2);
t385 = qJD(1) * t360;
t380 = t351 * t385;
t261 = t342 * pkin(2) - t343 * qJ(3) - t325 * t380 + t344 * t401 - t276;
t398 = t353 * g(3);
t397 = mrSges(3,1) - mrSges(4,2);
t396 = Ifges(3,4) + Ifges(4,6);
t395 = -pkin(2) - qJ(4);
t394 = t351 ^ 2 * t362;
t392 = t351 * t360;
t390 = t353 * t360;
t320 = pkin(3) * t379 - t344 * qJ(4);
t328 = (qJD(2) * t385 + qJDD(1) * t356) * t351;
t329 = -qJD(2) * t379 + t360 * t383;
t382 = t360 ^ 2 * t394;
t241 = -pkin(3) * t382 - t398 - t328 * qJ(3) + t395 * t329 + (-t323 + (-qJ(3) * t344 * t360 - t320 * t356) * qJD(1)) * t351 + t400;
t387 = g(3) * t392 + t356 * t324;
t375 = -t342 * qJ(3) + t325 * t379 + qJDD(3) + t387;
t244 = t328 * pkin(3) + t395 * t343 + (-pkin(3) * t344 * t386 - qJ(4) * t356 * t394 - t323 * t353) * t360 + t375;
t350 = sin(pkin(11));
t352 = cos(pkin(11));
t308 = t352 * t344 - t350 * t380;
t226 = -0.2e1 * qJD(4) * t308 - t350 * t241 + t352 * t244;
t288 = -t350 * t329 + t352 * t343;
t307 = -t350 * t344 - t352 * t380;
t222 = (t307 * t379 - t288) * pkin(9) + (t307 * t308 + t328) * pkin(4) + t226;
t227 = 0.2e1 * qJD(4) * t307 + t352 * t241 + t350 * t244;
t287 = -t352 * t329 - t350 * t343;
t289 = pkin(4) * t379 - t308 * pkin(9);
t306 = t307 ^ 2;
t224 = -t306 * pkin(4) + t287 * pkin(9) - t289 * t379 + t227;
t355 = sin(qJ(5));
t359 = cos(qJ(5));
t219 = t355 * t222 + t359 * t224;
t279 = t355 * t307 + t359 * t308;
t250 = -t279 * qJD(5) + t359 * t287 - t355 * t288;
t278 = t359 * t307 - t355 * t308;
t259 = -t278 * mrSges(6,1) + t279 * mrSges(6,2);
t334 = qJD(5) + t379;
t269 = t334 * mrSges(6,1) - t279 * mrSges(6,3);
t317 = qJDD(5) + t328;
t260 = -t278 * pkin(5) - t279 * pkin(10);
t332 = t334 ^ 2;
t216 = -t332 * pkin(5) + t317 * pkin(10) + t278 * t260 + t219;
t240 = t329 * pkin(3) - qJ(4) * t382 + t344 * t320 + qJDD(4) - t261;
t229 = -t287 * pkin(4) - t306 * pkin(9) + t308 * t289 + t240;
t251 = t278 * qJD(5) + t355 * t287 + t359 * t288;
t220 = t229 + (-t278 * t334 - t251) * pkin(10) + (t279 * t334 - t250) * pkin(5);
t354 = sin(qJ(6));
t358 = cos(qJ(6));
t213 = -t354 * t216 + t358 * t220;
t266 = -t354 * t279 + t358 * t334;
t232 = t266 * qJD(6) + t358 * t251 + t354 * t317;
t267 = t358 * t279 + t354 * t334;
t245 = -t266 * mrSges(7,1) + t267 * mrSges(7,2);
t249 = qJDD(6) - t250;
t277 = qJD(6) - t278;
t252 = -t277 * mrSges(7,2) + t266 * mrSges(7,3);
t209 = m(7) * t213 + t249 * mrSges(7,1) - t232 * mrSges(7,3) - t267 * t245 + t277 * t252;
t214 = t358 * t216 + t354 * t220;
t231 = -t267 * qJD(6) - t354 * t251 + t358 * t317;
t253 = t277 * mrSges(7,1) - t267 * mrSges(7,3);
t210 = m(7) * t214 - t249 * mrSges(7,2) + t231 * mrSges(7,3) + t266 * t245 - t277 * t253;
t377 = -t354 * t209 + t358 * t210;
t195 = m(6) * t219 - t317 * mrSges(6,2) + t250 * mrSges(6,3) + t278 * t259 - t334 * t269 + t377;
t218 = t359 * t222 - t355 * t224;
t268 = -t334 * mrSges(6,2) + t278 * mrSges(6,3);
t215 = -t317 * pkin(5) - t332 * pkin(10) + t279 * t260 - t218;
t373 = -m(7) * t215 + t231 * mrSges(7,1) - t232 * mrSges(7,2) + t266 * t252 - t267 * t253;
t205 = m(6) * t218 + t317 * mrSges(6,1) - t251 * mrSges(6,3) - t279 * t259 + t334 * t268 + t373;
t189 = t355 * t195 + t359 * t205;
t198 = t358 * t209 + t354 * t210;
t295 = Ifges(4,1) * t344 + (-Ifges(4,4) * t356 - Ifges(4,5) * t360) * t386;
t389 = Ifges(3,3) * t344 + (Ifges(3,5) * t356 + Ifges(3,6) * t360) * t386 + t295;
t293 = Ifges(4,5) * t344 + (-Ifges(4,6) * t356 - Ifges(4,3) * t360) * t386;
t388 = -Ifges(3,6) * t344 - (Ifges(3,4) * t356 + Ifges(3,2) * t360) * t386 + t293;
t381 = t323 * t390;
t275 = t381 - t387;
t319 = -t344 * mrSges(3,2) + mrSges(3,3) * t380;
t321 = -mrSges(4,1) * t380 - t344 * mrSges(4,3);
t326 = (mrSges(4,2) * t360 - mrSges(4,3) * t356) * t386;
t327 = (-mrSges(3,1) * t360 + mrSges(3,2) * t356) * t386;
t280 = -t307 * mrSges(5,1) + t308 * mrSges(5,2);
t285 = -mrSges(5,2) * t379 + t307 * mrSges(5,3);
t186 = m(5) * t226 + t328 * mrSges(5,1) - t288 * mrSges(5,3) - t308 * t280 + t285 * t379 + t189;
t286 = mrSges(5,1) * t379 - t308 * mrSges(5,3);
t378 = t359 * t195 - t355 * t205;
t187 = m(5) * t227 - t328 * mrSges(5,2) + t287 * mrSges(5,3) + t307 * t280 - t286 * t379 + t378;
t181 = t352 * t186 + t350 * t187;
t264 = -t343 * pkin(2) + t375 - t381;
t374 = -m(4) * t264 - t328 * mrSges(4,1) - t181;
t179 = m(3) * t275 - t328 * mrSges(3,3) + (t319 - t321) * t344 + t397 * t343 + (-t326 - t327) * t379 + t374;
t318 = t344 * mrSges(3,1) - mrSges(3,3) * t379;
t370 = m(6) * t229 - t250 * mrSges(6,1) + t251 * mrSges(6,2) - t278 * t268 + t279 * t269 + t198;
t196 = -m(5) * t240 + t287 * mrSges(5,1) - t288 * mrSges(5,2) + t307 * t285 - t308 * t286 - t370;
t322 = mrSges(4,1) * t379 + t344 * mrSges(4,2);
t365 = -m(4) * t261 + t343 * mrSges(4,3) + t344 * t322 + t326 * t380 - t196;
t192 = (mrSges(3,3) + mrSges(4,1)) * t329 + t327 * t380 - t343 * mrSges(3,2) - t344 * t318 + m(3) * t276 + t365;
t176 = -t356 * t179 + t360 * t192;
t182 = -t350 * t186 + t352 * t187;
t296 = -t351 * t323 - t398;
t262 = -t329 * pkin(2) + (-t344 * t380 - t328) * qJ(3) + t296 + t400;
t376 = m(4) * t262 - t328 * mrSges(4,3) + t321 * t380 + t182;
t178 = m(3) * t296 + t328 * mrSges(3,2) - t397 * t329 + (-t319 * t360 + (t318 - t322) * t356) * t386 + t376;
t170 = -t351 * t178 + t179 * t390 + t192 * t391;
t292 = Ifges(3,5) * t344 + (Ifges(3,1) * t356 + Ifges(3,4) * t360) * t386;
t233 = Ifges(7,5) * t267 + Ifges(7,6) * t266 + Ifges(7,3) * t277;
t235 = Ifges(7,1) * t267 + Ifges(7,4) * t266 + Ifges(7,5) * t277;
t202 = -mrSges(7,1) * t215 + mrSges(7,3) * t214 + Ifges(7,4) * t232 + Ifges(7,2) * t231 + Ifges(7,6) * t249 - t267 * t233 + t277 * t235;
t234 = Ifges(7,4) * t267 + Ifges(7,2) * t266 + Ifges(7,6) * t277;
t203 = mrSges(7,2) * t215 - mrSges(7,3) * t213 + Ifges(7,1) * t232 + Ifges(7,4) * t231 + Ifges(7,5) * t249 + t266 * t233 - t277 * t234;
t254 = Ifges(6,5) * t279 + Ifges(6,6) * t278 + Ifges(6,3) * t334;
t255 = Ifges(6,4) * t279 + Ifges(6,2) * t278 + Ifges(6,6) * t334;
t183 = mrSges(6,2) * t229 - mrSges(6,3) * t218 + Ifges(6,1) * t251 + Ifges(6,4) * t250 + Ifges(6,5) * t317 - pkin(10) * t198 - t354 * t202 + t358 * t203 + t278 * t254 - t334 * t255;
t256 = Ifges(6,1) * t279 + Ifges(6,4) * t278 + Ifges(6,5) * t334;
t366 = mrSges(7,1) * t213 - mrSges(7,2) * t214 + Ifges(7,5) * t232 + Ifges(7,6) * t231 + Ifges(7,3) * t249 + t267 * t234 - t266 * t235;
t184 = -mrSges(6,1) * t229 + mrSges(6,3) * t219 + Ifges(6,4) * t251 + Ifges(6,2) * t250 + Ifges(6,6) * t317 - pkin(5) * t198 - t279 * t254 + t334 * t256 - t366;
t270 = Ifges(5,5) * t308 + Ifges(5,6) * t307 + Ifges(5,3) * t379;
t272 = Ifges(5,1) * t308 + Ifges(5,4) * t307 + Ifges(5,5) * t379;
t171 = -mrSges(5,1) * t240 + mrSges(5,3) * t227 + Ifges(5,4) * t288 + Ifges(5,2) * t287 + Ifges(5,6) * t328 - pkin(4) * t370 + pkin(9) * t378 + t355 * t183 + t359 * t184 - t308 * t270 + t272 * t379;
t271 = Ifges(5,4) * t308 + Ifges(5,2) * t307 + Ifges(5,6) * t379;
t173 = mrSges(5,2) * t240 - mrSges(5,3) * t226 + Ifges(5,1) * t288 + Ifges(5,4) * t287 + Ifges(5,5) * t328 - pkin(9) * t189 + t359 * t183 - t355 * t184 + t307 * t270 - t271 * t379;
t294 = Ifges(4,4) * t344 + (-Ifges(4,2) * t356 - Ifges(4,6) * t360) * t386;
t369 = mrSges(4,2) * t264 - mrSges(4,3) * t261 + Ifges(4,1) * t343 - Ifges(4,4) * t328 - Ifges(4,5) * t329 - qJ(4) * t181 - t350 * t171 + t352 * t173 + t294 * t380;
t162 = (-t292 * t360 + (-pkin(2) * t326 - t388) * t356) * t386 + Ifges(3,3) * t343 + Ifges(3,5) * t328 + Ifges(3,6) * t329 + mrSges(3,1) * t275 - mrSges(3,2) * t276 + t369 + qJ(3) * (t329 * mrSges(4,1) + t365) + pkin(2) * (-t343 * mrSges(4,2) - t344 * t321 + t374);
t180 = t329 * mrSges(4,2) - t322 * t379 + t376;
t367 = -mrSges(4,1) * t261 + mrSges(4,2) * t262 - pkin(3) * t196 - qJ(4) * t182 - t352 * t171 - t350 * t173;
t164 = -mrSges(3,1) * t296 + mrSges(3,3) * t276 - pkin(2) * t180 + (t292 - t294) * t344 + (Ifges(3,6) - Ifges(4,5)) * t343 + (Ifges(3,2) + Ifges(4,3)) * t329 + t396 * t328 - t389 * t379 + t367;
t368 = -mrSges(6,1) * t218 + mrSges(6,2) * t219 - Ifges(6,5) * t251 - Ifges(6,6) * t250 - Ifges(6,3) * t317 - pkin(5) * t373 - pkin(10) * t377 - t358 * t202 - t354 * t203 - t279 * t255 + t278 * t256;
t364 = -mrSges(5,1) * t226 + mrSges(5,2) * t227 - Ifges(5,5) * t288 - Ifges(5,6) * t287 - Ifges(5,3) * t328 - pkin(4) * t189 - t308 * t271 + t307 * t272 + t368;
t363 = -mrSges(4,1) * t264 + mrSges(4,3) * t262 - pkin(3) * t181 + t364;
t166 = t388 * t344 + (Ifges(3,5) - Ifges(4,4)) * t343 + t396 * t329 + (Ifges(3,1) + Ifges(4,2)) * t328 + t389 * t380 + mrSges(3,2) * t296 - mrSges(3,3) * t275 - qJ(3) * t180 - t363;
t372 = mrSges(2,1) * t339 - mrSges(2,2) * t340 + Ifges(2,3) * qJDD(1) + pkin(1) * t170 + t353 * t162 + t164 * t392 + t166 * t393 + t176 * t399;
t174 = m(2) * t340 - t362 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t176;
t169 = t353 * t178 + (t179 * t360 + t192 * t356) * t351;
t167 = m(2) * t339 + qJDD(1) * mrSges(2,1) - t362 * mrSges(2,2) + t170;
t160 = -mrSges(2,2) * g(3) - mrSges(2,3) * t339 + Ifges(2,5) * qJDD(1) - t362 * Ifges(2,6) - t356 * t164 + t360 * t166 + (-t169 * t351 - t170 * t353) * pkin(8);
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t340 + t362 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t169 - t351 * t162 + (pkin(8) * t176 + t164 * t360 + t166 * t356) * t353;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t361 * t160 - t357 * t159 - pkin(7) * (t361 * t167 + t357 * t174), t160, t166, -t293 * t379 + t369, t173, t183, t203; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t357 * t160 + t361 * t159 + pkin(7) * (-t357 * t167 + t361 * t174), t159, t164, Ifges(4,4) * t343 - Ifges(4,2) * t328 - Ifges(4,6) * t329 - t344 * t293 - t295 * t380 + t363, t171, t184, t202; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t372, t372, t162, Ifges(4,5) * t343 - Ifges(4,6) * t328 - Ifges(4,3) * t329 + t344 * t294 + t295 * t379 - t367, -t364, -t368, t366;];
m_new  = t1;
