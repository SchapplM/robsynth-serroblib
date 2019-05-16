% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 00:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:37:22
% EndTime: 2019-05-07 00:38:18
% DurationCPUTime: 28.62s
% Computational Cost: add. (506469->406), mult. (1130788->505), div. (0->0), fcn. (847512->12), ass. (0->164)
t401 = -2 * qJD(3);
t356 = sin(qJ(2));
t351 = sin(pkin(6));
t386 = qJD(1) * t351;
t343 = t356 * t386;
t352 = cos(pkin(6));
t345 = qJD(1) * t352 + qJD(2);
t400 = (pkin(2) * t345 + t401) * t343;
t357 = sin(qJ(1));
t362 = cos(qJ(1));
t339 = t357 * g(1) - g(2) * t362;
t363 = qJD(1) ^ 2;
t398 = pkin(8) * t351;
t321 = qJDD(1) * pkin(1) + t363 * t398 + t339;
t340 = -g(1) * t362 - g(2) * t357;
t383 = qJDD(1) * t351;
t322 = -pkin(1) * t363 + pkin(8) * t383 + t340;
t361 = cos(qJ(2));
t391 = t352 * t356;
t393 = t351 * t356;
t281 = -g(3) * t393 + t321 * t391 + t361 * t322;
t323 = (-pkin(2) * t361 - qJ(3) * t356) * t386;
t342 = t345 ^ 2;
t344 = qJDD(1) * t352 + qJDD(2);
t385 = qJD(1) * t361;
t380 = t351 * t385;
t259 = pkin(2) * t342 - t344 * qJ(3) - t323 * t380 + t345 * t401 - t281;
t399 = -pkin(2) - pkin(9);
t397 = g(3) * t352;
t396 = mrSges(3,1) - mrSges(4,2);
t395 = Ifges(3,4) + Ifges(4,6);
t394 = t351 ^ 2 * t363;
t392 = t351 * t361;
t390 = t352 * t361;
t326 = pkin(3) * t343 - pkin(9) * t345;
t327 = (qJD(2) * t385 + qJDD(1) * t356) * t351;
t328 = -qJD(2) * t343 + t361 * t383;
t381 = t361 ^ 2 * t394;
t247 = -pkin(3) * t381 - t397 - qJ(3) * t327 + t399 * t328 + (-t321 + (-qJ(3) * t345 * t361 - t326 * t356) * qJD(1)) * t351 + t400;
t387 = g(3) * t392 + t356 * t322;
t376 = -qJ(3) * t342 + t323 * t343 + qJDD(3) + t387;
t250 = pkin(3) * t327 + t399 * t344 + (-pkin(3) * t345 * t386 - pkin(9) * t356 * t394 - t321 * t352) * t361 + t376;
t355 = sin(qJ(4));
t360 = cos(qJ(4));
t226 = -t247 * t355 + t360 * t250;
t306 = -t345 * t355 - t360 * t380;
t279 = qJD(4) * t306 - t328 * t355 + t344 * t360;
t307 = t345 * t360 - t355 * t380;
t316 = qJDD(4) + t327;
t334 = t343 + qJD(4);
t222 = (t306 * t334 - t279) * pkin(10) + (t306 * t307 + t316) * pkin(4) + t226;
t227 = t360 * t247 + t355 * t250;
t278 = -qJD(4) * t307 - t328 * t360 - t344 * t355;
t288 = pkin(4) * t334 - pkin(10) * t307;
t305 = t306 ^ 2;
t224 = -pkin(4) * t305 + pkin(10) * t278 - t288 * t334 + t227;
t354 = sin(qJ(5));
t359 = cos(qJ(5));
t219 = t354 * t222 + t359 * t224;
t284 = t306 * t354 + t307 * t359;
t241 = -qJD(5) * t284 + t278 * t359 - t279 * t354;
t283 = t306 * t359 - t307 * t354;
t261 = -mrSges(6,1) * t283 + mrSges(6,2) * t284;
t331 = qJD(5) + t334;
t269 = mrSges(6,1) * t331 - mrSges(6,3) * t284;
t311 = qJDD(5) + t316;
t262 = -pkin(5) * t283 - pkin(11) * t284;
t330 = t331 ^ 2;
t216 = -pkin(5) * t330 + pkin(11) * t311 + t262 * t283 + t219;
t246 = pkin(3) * t328 - pkin(9) * t381 + t345 * t326 - t259;
t229 = -pkin(4) * t278 - pkin(10) * t305 + t307 * t288 + t246;
t242 = qJD(5) * t283 + t278 * t354 + t279 * t359;
t220 = (t284 * t331 - t241) * pkin(5) + (-t283 * t331 - t242) * pkin(11) + t229;
t353 = sin(qJ(6));
t358 = cos(qJ(6));
t213 = -t216 * t353 + t220 * t358;
t266 = -t284 * t353 + t331 * t358;
t232 = qJD(6) * t266 + t242 * t358 + t311 * t353;
t240 = qJDD(6) - t241;
t267 = t284 * t358 + t331 * t353;
t251 = -mrSges(7,1) * t266 + mrSges(7,2) * t267;
t282 = qJD(6) - t283;
t252 = -mrSges(7,2) * t282 + mrSges(7,3) * t266;
t209 = m(7) * t213 + mrSges(7,1) * t240 - mrSges(7,3) * t232 - t251 * t267 + t252 * t282;
t214 = t216 * t358 + t220 * t353;
t231 = -qJD(6) * t267 - t242 * t353 + t311 * t358;
t253 = mrSges(7,1) * t282 - mrSges(7,3) * t267;
t210 = m(7) * t214 - mrSges(7,2) * t240 + mrSges(7,3) * t231 + t251 * t266 - t253 * t282;
t378 = -t209 * t353 + t358 * t210;
t195 = m(6) * t219 - mrSges(6,2) * t311 + mrSges(6,3) * t241 + t261 * t283 - t269 * t331 + t378;
t218 = t222 * t359 - t224 * t354;
t268 = -mrSges(6,2) * t331 + mrSges(6,3) * t283;
t215 = -pkin(5) * t311 - pkin(11) * t330 + t262 * t284 - t218;
t374 = -m(7) * t215 + t231 * mrSges(7,1) - mrSges(7,2) * t232 + t266 * t252 - t253 * t267;
t205 = m(6) * t218 + mrSges(6,1) * t311 - mrSges(6,3) * t242 - t261 * t284 + t268 * t331 + t374;
t189 = t354 * t195 + t359 * t205;
t198 = t358 * t209 + t353 * t210;
t295 = Ifges(4,1) * t345 + (-Ifges(4,4) * t356 - Ifges(4,5) * t361) * t386;
t389 = Ifges(3,3) * t345 + (Ifges(3,5) * t356 + Ifges(3,6) * t361) * t386 + t295;
t293 = Ifges(4,5) * t345 + (-Ifges(4,6) * t356 - Ifges(4,3) * t361) * t386;
t388 = -Ifges(3,6) * t345 - (Ifges(3,4) * t356 + Ifges(3,2) * t361) * t386 + t293;
t382 = t321 * t390;
t280 = t382 - t387;
t318 = -mrSges(3,2) * t345 + mrSges(3,3) * t380;
t319 = -mrSges(4,1) * t380 - mrSges(4,3) * t345;
t324 = (mrSges(4,2) * t361 - mrSges(4,3) * t356) * t386;
t325 = (-mrSges(3,1) * t361 + mrSges(3,2) * t356) * t386;
t285 = -mrSges(5,1) * t306 + mrSges(5,2) * t307;
t286 = -mrSges(5,2) * t334 + mrSges(5,3) * t306;
t186 = m(5) * t226 + mrSges(5,1) * t316 - mrSges(5,3) * t279 - t285 * t307 + t286 * t334 + t189;
t287 = mrSges(5,1) * t334 - mrSges(5,3) * t307;
t379 = t359 * t195 - t205 * t354;
t187 = m(5) * t227 - mrSges(5,2) * t316 + mrSges(5,3) * t278 + t285 * t306 - t287 * t334 + t379;
t181 = t186 * t360 + t187 * t355;
t264 = -pkin(2) * t344 + t376 - t382;
t375 = -m(4) * t264 - t327 * mrSges(4,1) - t181;
t179 = m(3) * t280 - mrSges(3,3) * t327 + (t318 - t319) * t345 + t396 * t344 + (-t324 - t325) * t343 + t375;
t317 = mrSges(3,1) * t345 - mrSges(3,3) * t343;
t371 = m(6) * t229 - t241 * mrSges(6,1) + t242 * mrSges(6,2) - t283 * t268 + t284 * t269 + t198;
t196 = -m(5) * t246 + t278 * mrSges(5,1) - t279 * mrSges(5,2) + t306 * t286 - t307 * t287 - t371;
t320 = mrSges(4,1) * t343 + mrSges(4,2) * t345;
t366 = -m(4) * t259 + t344 * mrSges(4,3) + t345 * t320 + t324 * t380 - t196;
t192 = t366 + t325 * t380 - t344 * mrSges(3,2) - t345 * t317 + m(3) * t281 + (mrSges(3,3) + mrSges(4,1)) * t328;
t176 = -t179 * t356 + t361 * t192;
t182 = -t186 * t355 + t360 * t187;
t296 = -t321 * t351 - t397;
t260 = -pkin(2) * t328 + (-t345 * t380 - t327) * qJ(3) + t296 + t400;
t377 = m(4) * t260 - t327 * mrSges(4,3) + t319 * t380 + t182;
t178 = m(3) * t296 + mrSges(3,2) * t327 - t396 * t328 + (-t318 * t361 + (t317 - t320) * t356) * t386 + t377;
t170 = -t178 * t351 + t179 * t390 + t192 * t391;
t292 = Ifges(3,5) * t345 + (Ifges(3,1) * t356 + Ifges(3,4) * t361) * t386;
t233 = Ifges(7,5) * t267 + Ifges(7,6) * t266 + Ifges(7,3) * t282;
t235 = Ifges(7,1) * t267 + Ifges(7,4) * t266 + Ifges(7,5) * t282;
t202 = -mrSges(7,1) * t215 + mrSges(7,3) * t214 + Ifges(7,4) * t232 + Ifges(7,2) * t231 + Ifges(7,6) * t240 - t233 * t267 + t235 * t282;
t234 = Ifges(7,4) * t267 + Ifges(7,2) * t266 + Ifges(7,6) * t282;
t203 = mrSges(7,2) * t215 - mrSges(7,3) * t213 + Ifges(7,1) * t232 + Ifges(7,4) * t231 + Ifges(7,5) * t240 + t233 * t266 - t234 * t282;
t254 = Ifges(6,5) * t284 + Ifges(6,6) * t283 + Ifges(6,3) * t331;
t255 = Ifges(6,4) * t284 + Ifges(6,2) * t283 + Ifges(6,6) * t331;
t183 = mrSges(6,2) * t229 - mrSges(6,3) * t218 + Ifges(6,1) * t242 + Ifges(6,4) * t241 + Ifges(6,5) * t311 - pkin(11) * t198 - t202 * t353 + t203 * t358 + t254 * t283 - t255 * t331;
t256 = Ifges(6,1) * t284 + Ifges(6,4) * t283 + Ifges(6,5) * t331;
t367 = mrSges(7,1) * t213 - mrSges(7,2) * t214 + Ifges(7,5) * t232 + Ifges(7,6) * t231 + Ifges(7,3) * t240 + t234 * t267 - t235 * t266;
t184 = -mrSges(6,1) * t229 + mrSges(6,3) * t219 + Ifges(6,4) * t242 + Ifges(6,2) * t241 + Ifges(6,6) * t311 - pkin(5) * t198 - t254 * t284 + t256 * t331 - t367;
t270 = Ifges(5,5) * t307 + Ifges(5,6) * t306 + Ifges(5,3) * t334;
t272 = Ifges(5,1) * t307 + Ifges(5,4) * t306 + Ifges(5,5) * t334;
t171 = -mrSges(5,1) * t246 + mrSges(5,3) * t227 + Ifges(5,4) * t279 + Ifges(5,2) * t278 + Ifges(5,6) * t316 - pkin(4) * t371 + pkin(10) * t379 + t354 * t183 + t359 * t184 - t307 * t270 + t334 * t272;
t271 = Ifges(5,4) * t307 + Ifges(5,2) * t306 + Ifges(5,6) * t334;
t173 = mrSges(5,2) * t246 - mrSges(5,3) * t226 + Ifges(5,1) * t279 + Ifges(5,4) * t278 + Ifges(5,5) * t316 - pkin(10) * t189 + t183 * t359 - t184 * t354 + t270 * t306 - t271 * t334;
t294 = Ifges(4,4) * t345 + (-Ifges(4,2) * t356 - Ifges(4,6) * t361) * t386;
t370 = mrSges(4,2) * t264 - mrSges(4,3) * t259 + Ifges(4,1) * t344 - Ifges(4,4) * t327 - Ifges(4,5) * t328 - pkin(9) * t181 - t355 * t171 + t360 * t173 + t294 * t380;
t162 = qJ(3) * (t328 * mrSges(4,1) + t366) + pkin(2) * (-mrSges(4,2) * t344 - t319 * t345 + t375) + (-t292 * t361 + (-pkin(2) * t324 - t388) * t356) * t386 + Ifges(3,3) * t344 + Ifges(3,5) * t327 + Ifges(3,6) * t328 + mrSges(3,1) * t280 - mrSges(3,2) * t281 + t370;
t180 = mrSges(4,2) * t328 - t320 * t343 + t377;
t368 = -mrSges(4,1) * t259 + mrSges(4,2) * t260 - pkin(3) * t196 - pkin(9) * t182 - t360 * t171 - t355 * t173;
t164 = -mrSges(3,1) * t296 + mrSges(3,3) * t281 - pkin(2) * t180 + (t292 - t294) * t345 + (Ifges(3,6) - Ifges(4,5)) * t344 + (Ifges(3,2) + Ifges(4,3)) * t328 + t395 * t327 - t389 * t343 + t368;
t369 = -mrSges(6,1) * t218 + mrSges(6,2) * t219 - Ifges(6,5) * t242 - Ifges(6,6) * t241 - Ifges(6,3) * t311 - pkin(5) * t374 - pkin(11) * t378 - t358 * t202 - t353 * t203 - t284 * t255 + t283 * t256;
t365 = -mrSges(5,1) * t226 + mrSges(5,2) * t227 - Ifges(5,5) * t279 - Ifges(5,6) * t278 - Ifges(5,3) * t316 - pkin(4) * t189 - t307 * t271 + t306 * t272 + t369;
t364 = -mrSges(4,1) * t264 + mrSges(4,3) * t260 - pkin(3) * t181 + t365;
t166 = t389 * t380 - t364 - qJ(3) * t180 + mrSges(3,2) * t296 - mrSges(3,3) * t280 + (Ifges(3,1) + Ifges(4,2)) * t327 + t395 * t328 + (Ifges(3,5) - Ifges(4,4)) * t344 + t388 * t345;
t373 = mrSges(2,1) * t339 - mrSges(2,2) * t340 + Ifges(2,3) * qJDD(1) + pkin(1) * t170 + t352 * t162 + t164 * t392 + t166 * t393 + t176 * t398;
t174 = m(2) * t340 - mrSges(2,1) * t363 - qJDD(1) * mrSges(2,2) + t176;
t169 = t178 * t352 + (t179 * t361 + t192 * t356) * t351;
t167 = m(2) * t339 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t363 + t170;
t160 = -mrSges(2,2) * g(3) - mrSges(2,3) * t339 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t363 - t164 * t356 + t166 * t361 + (-t169 * t351 - t170 * t352) * pkin(8);
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t340 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t169 - t351 * t162 + (pkin(8) * t176 + t164 * t361 + t166 * t356) * t352;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t362 * t160 - t357 * t159 - pkin(7) * (t167 * t362 + t174 * t357), t160, t166, -t293 * t343 + t370, t173, t183, t203; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t357 * t160 + t362 * t159 + pkin(7) * (-t357 * t167 + t174 * t362), t159, t164, Ifges(4,4) * t344 - Ifges(4,2) * t327 - Ifges(4,6) * t328 - t345 * t293 - t295 * t380 + t364, t171, t184, t202; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t373, t373, t162, Ifges(4,5) * t344 - Ifges(4,6) * t327 - Ifges(4,3) * t328 + t345 * t294 + t295 * t343 - t368, -t365, -t369, t367;];
m_new  = t1;
