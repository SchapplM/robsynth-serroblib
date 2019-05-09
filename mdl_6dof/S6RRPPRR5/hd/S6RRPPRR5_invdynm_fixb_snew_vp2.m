% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 10:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:38:19
% EndTime: 2019-05-06 10:38:48
% DurationCPUTime: 9.94s
% Computational Cost: add. (151878->417), mult. (350112->501), div. (0->0), fcn. (237043->10), ass. (0->161)
t377 = cos(pkin(6));
t368 = qJD(1) * t377 + qJD(2);
t380 = sin(qJ(2));
t384 = cos(qJ(2));
t376 = sin(pkin(6));
t418 = qJD(1) * t376;
t282 = -Ifges(5,3) * t368 + (Ifges(5,5) * t380 - Ifges(5,6) * t384) * t418;
t285 = -Ifges(5,6) * t368 + (Ifges(5,4) * t380 - Ifges(5,2) * t384) * t418;
t417 = qJD(1) * t384;
t412 = t376 * t417;
t443 = -t282 * t412 + t368 * t285;
t341 = (qJD(2) * t417 + qJDD(1) * t380) * t376;
t413 = t380 * t418;
t416 = qJDD(1) * t376;
t342 = -qJD(2) * t413 + t384 * t416;
t367 = qJDD(1) * t377 + qJDD(2);
t442 = -Ifges(5,5) * t341 + Ifges(5,6) * t342 + Ifges(5,3) * t367;
t441 = t368 ^ 2;
t440 = -2 * qJD(4);
t439 = pkin(2) * t368;
t438 = pkin(8) * t376;
t437 = -mrSges(5,2) - mrSges(4,3);
t436 = mrSges(3,3) + mrSges(4,2);
t328 = -mrSges(5,1) * t368 - mrSges(5,3) * t413;
t431 = t368 * t328;
t331 = mrSges(5,2) * t368 - mrSges(5,3) * t412;
t430 = t368 * t331;
t386 = qJD(1) ^ 2;
t429 = t376 ^ 2 * t386;
t428 = t376 * t380;
t427 = t376 * t384;
t426 = t377 * t380;
t425 = t377 * t384;
t327 = -pkin(3) * t368 - qJ(4) * t413;
t381 = sin(qJ(1));
t385 = cos(qJ(1));
t359 = t381 * g(1) - g(2) * t385;
t334 = qJDD(1) * pkin(1) + t386 * t438 + t359;
t291 = -t377 * g(3) - t376 * t334;
t409 = t368 * t412;
t407 = -t342 * pkin(2) + t291 + (-t341 - t409) * qJ(3);
t410 = 0.2e1 * t413;
t414 = t384 ^ 2 * t429;
t395 = -qJ(4) * t414 + qJD(3) * t410 + t327 * t413 + qJDD(4) - t407;
t228 = -t341 * pkin(9) + (pkin(3) + pkin(4)) * t342 + (-pkin(9) * t384 + (-pkin(2) - pkin(4)) * t380) * t368 * t418 + t395;
t340 = (pkin(4) * t384 - pkin(9) * t380) * t418;
t336 = (-pkin(2) * t384 - qJ(3) * t380) * t418;
t360 = -g(1) * t385 - g(2) * t381;
t335 = -pkin(1) * t386 + pkin(8) * t416 + t360;
t421 = t334 * t426 + t384 * t335;
t405 = -pkin(2) * t441 + t367 * qJ(3) + 0.2e1 * qJD(3) * t368 + t336 * t412 + t421;
t393 = -pkin(3) * t414 - t342 * qJ(4) + t368 * t327 + t412 * t440 + t405;
t230 = -t441 * pkin(4) - t367 * pkin(9) + (-g(3) * t380 - t340 * t417) * t376 + t393;
t379 = sin(qJ(5));
t383 = cos(qJ(5));
t225 = t379 * t228 + t383 * t230;
t305 = -t368 * t383 - t379 * t413;
t306 = -t368 * t379 + t383 * t413;
t274 = -pkin(5) * t305 - pkin(10) * t306;
t325 = qJDD(5) + t342;
t348 = qJD(5) + t412;
t346 = t348 ^ 2;
t221 = -pkin(5) * t346 + pkin(10) * t325 + t274 * t305 + t225;
t271 = -g(3) * t427 + t334 * t425 - t380 * t335;
t253 = -t367 * pkin(2) - qJ(3) * t441 + t336 * t413 + qJDD(3) - t271;
t400 = -t341 * qJ(4) + t253 + (-t380 * t384 * t429 - t367) * pkin(3);
t232 = t367 * pkin(4) - pkin(9) * t441 - qJ(4) * t409 + qJD(4) * t410 + t340 * t413 - t400;
t269 = -qJD(5) * t306 - t341 * t379 - t367 * t383;
t270 = qJD(5) * t305 + t341 * t383 - t367 * t379;
t226 = t232 + (t306 * t348 - t269) * pkin(5) + (-t305 * t348 - t270) * pkin(10);
t378 = sin(qJ(6));
t382 = cos(qJ(6));
t218 = -t221 * t378 + t226 * t382;
t275 = -t306 * t378 + t348 * t382;
t242 = qJD(6) * t275 + t270 * t382 + t325 * t378;
t276 = t306 * t382 + t348 * t378;
t254 = -mrSges(7,1) * t275 + mrSges(7,2) * t276;
t303 = qJD(6) - t305;
t257 = -mrSges(7,2) * t303 + mrSges(7,3) * t275;
t266 = qJDD(6) - t269;
t214 = m(7) * t218 + mrSges(7,1) * t266 - t242 * mrSges(7,3) - t254 * t276 + t257 * t303;
t219 = t221 * t382 + t226 * t378;
t241 = -qJD(6) * t276 - t270 * t378 + t325 * t382;
t258 = mrSges(7,1) * t303 - mrSges(7,3) * t276;
t215 = m(7) * t219 - mrSges(7,2) * t266 + t241 * mrSges(7,3) + t254 * t275 - t258 * t303;
t204 = t382 * t214 + t378 * t215;
t286 = Ifges(4,2) * t368 + (Ifges(4,4) * t380 - Ifges(4,6) * t384) * t418;
t424 = -t282 + t286;
t423 = t285 - Ifges(3,6) * t368 - (Ifges(3,4) * t380 + Ifges(3,2) * t384) * t418;
t289 = Ifges(4,4) * t368 + (Ifges(4,1) * t380 - Ifges(4,5) * t384) * t418;
t422 = t289 + Ifges(3,5) * t368 + (Ifges(3,1) * t380 + Ifges(3,4) * t384) * t418;
t420 = -mrSges(3,1) * t368 + mrSges(3,3) * t413 + t328;
t419 = -mrSges(3,2) * t368 + mrSges(3,3) * t412 + t331;
t415 = g(3) * t428;
t272 = -t415 + t421;
t338 = (mrSges(5,1) * t384 + mrSges(5,2) * t380) * t418;
t339 = (-mrSges(3,1) * t384 + mrSges(3,2) * t380) * t418;
t246 = t405 - t415;
t330 = -mrSges(4,1) * t368 + mrSges(4,2) * t413;
t337 = (-mrSges(4,1) * t384 - mrSges(4,3) * t380) * t418;
t273 = -mrSges(6,1) * t305 + mrSges(6,2) * t306;
t278 = mrSges(6,1) * t348 - mrSges(6,3) * t306;
t411 = -t214 * t378 + t382 * t215;
t201 = m(6) * t225 - mrSges(6,2) * t325 + mrSges(6,3) * t269 + t273 * t305 - t278 * t348 + t411;
t224 = t228 * t383 - t230 * t379;
t277 = -mrSges(6,2) * t348 + mrSges(6,3) * t305;
t220 = -pkin(5) * t325 - pkin(10) * t346 + t274 * t306 - t224;
t403 = -m(7) * t220 + t241 * mrSges(7,1) - t242 * mrSges(7,2) + t275 * t257 - t258 * t276;
t210 = m(6) * t224 + mrSges(6,1) * t325 - mrSges(6,3) * t270 - t273 * t306 + t277 * t348 + t403;
t196 = t383 * t201 - t379 * t210;
t236 = t393 - t415;
t406 = m(5) * t236 - t342 * mrSges(5,3) + t196;
t401 = m(4) * t246 + t367 * mrSges(4,3) + t368 * t330 + t337 * t412 + t406;
t190 = m(3) * t272 + t420 * t368 + (-mrSges(3,2) + mrSges(5,2)) * t367 + t436 * t342 + (-t338 + t339) * t412 + t401;
t333 = mrSges(4,2) * t412 + mrSges(4,3) * t368;
t202 = -m(6) * t232 + t269 * mrSges(6,1) - t270 * mrSges(6,2) + t305 * t277 - t306 * t278 - t204;
t237 = (qJ(4) * t368 * t384 + t380 * t440) * t418 + t400;
t397 = -m(5) * t237 + t341 * mrSges(5,3) + t338 * t413 - t202;
t392 = -m(4) * t253 + t367 * mrSges(4,1) + t368 * t333 + t397;
t198 = t392 + (mrSges(3,1) + mrSges(5,1)) * t367 + t419 * t368 - t436 * t341 + m(3) * t271 + (-t337 - t339) * t413;
t184 = t384 * t190 - t198 * t380;
t247 = (-0.2e1 * qJD(3) + t439) * t413 + t407;
t195 = t379 * t201 + t383 * t210;
t239 = t342 * pkin(3) - t413 * t439 + t395;
t404 = -m(5) * t239 - t342 * mrSges(5,1) - t195;
t399 = m(4) * t247 - t342 * mrSges(4,1) + t404;
t191 = m(3) * t291 - t342 * mrSges(3,1) + (mrSges(3,2) + t437) * t341 + ((-t333 - t419) * t384 + (-t330 - t420) * t380) * t418 + t399;
t181 = t190 * t426 - t191 * t376 + t198 * t425;
t283 = Ifges(4,6) * t368 + (Ifges(4,5) * t380 - Ifges(4,3) * t384) * t418;
t288 = -Ifges(5,5) * t368 + (Ifges(5,1) * t380 - Ifges(5,4) * t384) * t418;
t199 = -t367 * mrSges(5,1) - t397 - t430;
t248 = Ifges(7,5) * t276 + Ifges(7,6) * t275 + Ifges(7,3) * t303;
t250 = Ifges(7,1) * t276 + Ifges(7,4) * t275 + Ifges(7,5) * t303;
t208 = -mrSges(7,1) * t220 + mrSges(7,3) * t219 + Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t266 - t248 * t276 + t250 * t303;
t249 = Ifges(7,4) * t276 + Ifges(7,2) * t275 + Ifges(7,6) * t303;
t209 = mrSges(7,2) * t220 - mrSges(7,3) * t218 + Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t266 + t248 * t275 - t249 * t303;
t259 = Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t348;
t260 = Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t348;
t186 = mrSges(6,2) * t232 - mrSges(6,3) * t224 + Ifges(6,1) * t270 + Ifges(6,4) * t269 + Ifges(6,5) * t325 - pkin(10) * t204 - t208 * t378 + t209 * t382 + t259 * t305 - t260 * t348;
t261 = Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t348;
t391 = mrSges(7,1) * t218 - mrSges(7,2) * t219 + Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t266 + t249 * t276 - t250 * t275;
t187 = -mrSges(6,1) * t232 + mrSges(6,3) * t225 + Ifges(6,4) * t270 + Ifges(6,2) * t269 + Ifges(6,6) * t325 - pkin(5) * t204 - t259 * t306 + t261 * t348 - t391;
t394 = mrSges(5,1) * t237 - mrSges(5,2) * t236 + pkin(4) * t202 + pkin(9) * t196 + t379 * t186 + t383 * t187;
t389 = -mrSges(4,1) * t253 + mrSges(4,3) * t246 + Ifges(4,4) * t341 + Ifges(4,2) * t367 - Ifges(4,6) * t342 - pkin(3) * t199 - t394;
t173 = t389 + pkin(2) * (t392 + t430) + ((-qJ(3) * t338 - t288 - t422) * t384 + (-pkin(2) * t337 - t283 - t423) * t380) * t418 + (mrSges(5,1) * pkin(2) + mrSges(5,2) * qJ(3) + Ifges(3,3) + Ifges(5,3)) * t367 + (mrSges(4,2) * qJ(3) + Ifges(3,6) + Ifges(5,6)) * t342 + (-mrSges(4,2) * pkin(2) + Ifges(3,5) - Ifges(5,5)) * t341 + qJ(3) * (t401 + t431) + mrSges(3,1) * t271 - mrSges(3,2) * t272;
t192 = t437 * t341 + ((-t331 - t333) * t384 + (-t328 - t330) * t380) * t418 + t399;
t284 = Ifges(3,3) * t368 + (Ifges(3,5) * t380 + Ifges(3,6) * t384) * t418;
t398 = mrSges(6,1) * t224 - mrSges(6,2) * t225 + Ifges(6,5) * t270 + Ifges(6,6) * t269 + Ifges(6,3) * t325 + pkin(5) * t403 + pkin(10) * t411 + t382 * t208 + t378 * t209 + t306 * t260 - t305 * t261;
t390 = mrSges(5,1) * t239 - mrSges(5,3) * t236 - Ifges(5,4) * t341 + Ifges(5,2) * t342 + Ifges(5,6) * t367 + pkin(4) * t195 + t368 * t288 + t398;
t387 = mrSges(4,1) * t247 - mrSges(4,2) * t246 + pkin(3) * (-t341 * mrSges(5,2) + (-t328 * t380 - t331 * t384) * t418 + t404) + qJ(4) * (t367 * mrSges(5,2) - t338 * t412 + t406 + t431) - t390;
t175 = -t387 + (Ifges(3,6) - Ifges(4,6)) * t367 + t422 * t368 + (-t284 - t424) * t413 + (Ifges(4,3) + Ifges(3,2)) * t342 + (Ifges(3,4) - Ifges(4,5)) * t341 - mrSges(3,1) * t291 + mrSges(3,3) * t272 - pkin(2) * t192;
t396 = mrSges(5,2) * t239 - mrSges(5,3) * t237 + Ifges(5,1) * t341 - Ifges(5,4) * t342 - Ifges(5,5) * t367 - pkin(9) * t195 + t383 * t186 - t379 * t187;
t388 = mrSges(4,2) * t253 - mrSges(4,3) * t247 + Ifges(4,1) * t341 + Ifges(4,4) * t367 - Ifges(4,5) * t342 - qJ(4) * t199 + t368 * t283 + t286 * t412 + t396;
t177 = t388 + t423 * t368 + Ifges(3,5) * t367 + Ifges(3,1) * t341 + Ifges(3,4) * t342 + mrSges(3,2) * t291 - mrSges(3,3) * t271 - qJ(3) * t192 + (-t282 + t284) * t412;
t402 = mrSges(2,1) * t359 - mrSges(2,2) * t360 + Ifges(2,3) * qJDD(1) + pkin(1) * t181 + t377 * t173 + t175 * t427 + t177 * t428 + t184 * t438;
t182 = m(2) * t360 - mrSges(2,1) * t386 - qJDD(1) * mrSges(2,2) + t184;
t180 = t377 * t191 + (t190 * t380 + t198 * t384) * t376;
t178 = m(2) * t359 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t386 + t181;
t171 = -mrSges(2,2) * g(3) - mrSges(2,3) * t359 + Ifges(2,5) * qJDD(1) - t386 * Ifges(2,6) - t380 * t175 + t384 * t177 + (-t180 * t376 - t181 * t377) * pkin(8);
t170 = mrSges(2,1) * g(3) + mrSges(2,3) * t360 + t386 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t180 - t376 * t173 + (pkin(8) * t184 + t175 * t384 + t177 * t380) * t377;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t385 * t171 - t381 * t170 - pkin(7) * (t178 * t385 + t182 * t381), t171, t177, t388 + t443, t396 + t443, t186, t209; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t381 * t171 + t385 * t170 + pkin(7) * (-t178 * t381 + t182 * t385), t170, t175, t389 + ((-t288 - t289) * t384 + (-t283 - t285) * t380) * t418 + t442, -t282 * t413 - t390, t187, t208; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t402, t402, t173, Ifges(4,5) * t341 + Ifges(4,6) * t367 - Ifges(4,3) * t342 - t368 * t289 + t424 * t413 + t387, (t285 * t380 + t288 * t384) * t418 + t394 - t442, t398, t391;];
m_new  = t1;
