% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:04:36
% EndTime: 2019-03-09 03:04:45
% DurationCPUTime: 4.95s
% Computational Cost: add. (10312->412), mult. (19792->543), div. (0->0), fcn. (20522->8), ass. (0->220)
t393 = sin(qJ(3));
t333 = t393 * pkin(3);
t437 = m(5) * t333;
t429 = Ifges(7,4) + Ifges(6,5);
t436 = Ifges(6,6) - Ifges(7,6);
t359 = sin(pkin(10));
t360 = cos(pkin(10));
t394 = cos(qJ(3));
t218 = -t359 * t394 - t360 * t393;
t242 = sin(qJ(5));
t350 = t218 * t242;
t329 = mrSges(7,2) * t350;
t216 = t359 * t393 - t360 * t394;
t375 = t216 * mrSges(7,3);
t160 = t329 + t375;
t332 = mrSges(6,3) * t350;
t376 = t216 * mrSges(6,2);
t268 = t332 - t376;
t435 = t160 + t268;
t239 = t242 ^ 2;
t243 = cos(qJ(5));
t240 = t243 ^ 2;
t434 = t239 + t240;
t379 = Ifges(6,4) * t242;
t290 = Ifges(6,1) * t243 - t379;
t235 = Ifges(7,5) * t242;
t301 = Ifges(7,1) * t243 + t235;
t433 = -t301 - t290;
t320 = t359 * pkin(3);
t322 = t360 * pkin(3);
t432 = mrSges(7,2) + mrSges(6,3);
t220 = -t243 * mrSges(7,1) - t242 * mrSges(7,3);
t221 = -t243 * mrSges(6,1) + t242 * mrSges(6,2);
t252 = -t220 - t221;
t282 = t243 * pkin(5) + t242 * qJ(6);
t431 = (m(7) * t282 + t252) * qJD(5);
t430 = -t220 / 0.2e1;
t385 = mrSges(6,1) + mrSges(7,1);
t428 = Ifges(7,2) + Ifges(6,3);
t390 = t216 * pkin(5);
t324 = -cos(pkin(9)) * pkin(1) - pkin(2);
t219 = -pkin(3) * t394 + t324;
t133 = t216 * pkin(4) + t218 * pkin(8) + t219;
t232 = sin(pkin(9)) * pkin(1) + pkin(7);
t321 = t393 * t232;
t214 = -qJ(4) * t393 - t321;
t323 = t394 * t232;
t215 = qJ(4) * t394 + t323;
t419 = t359 * t214 + t360 * t215;
t62 = t133 * t243 - t242 * t419;
t48 = -t62 - t390;
t383 = t48 + t62;
t352 = t216 * t242;
t157 = mrSges(6,2) * t218 + mrSges(6,3) * t352;
t161 = mrSges(7,2) * t352 - t218 * mrSges(7,3);
t427 = t157 + t161;
t236 = Ifges(6,4) * t243;
t287 = Ifges(6,2) * t242 - t236;
t424 = t242 * Ifges(6,1) + t236;
t426 = t424 - t287;
t425 = -t243 * Ifges(7,3) + t235;
t423 = t429 * t216 + t433 * t218;
t400 = -t218 / 0.2e1;
t422 = (t429 * t242 + t436 * t243) * t400;
t421 = -t436 * t242 + t429 * t243;
t231 = t320 + pkin(8);
t420 = t434 * t218 * t231 / 0.2e1;
t224 = t243 * Ifges(6,2) + t379;
t378 = Ifges(7,5) * t243;
t225 = Ifges(7,1) * t242 - t378;
t396 = t243 / 0.2e1;
t399 = -t242 / 0.2e1;
t418 = t224 * t399 + t225 * t396;
t144 = -t360 * t214 + t215 * t359;
t156 = -t218 * pkin(4) + t216 * pkin(8) + t333;
t68 = t144 * t242 + t156 * t243;
t69 = -t144 * t243 + t242 * t156;
t275 = -t242 * t68 + t243 * t69;
t52 = -qJ(6) * t218 + t69;
t53 = t218 * pkin(5) - t68;
t277 = t242 * t53 + t243 * t52;
t416 = t425 - t433;
t234 = m(7) * qJ(6) + mrSges(7,3);
t366 = t243 * mrSges(7,3);
t372 = t242 * mrSges(7,1);
t291 = -t366 + t372;
t154 = t291 * t218;
t358 = qJ(6) * t243;
t389 = t242 * pkin(5);
t281 = -t358 + t389;
t73 = -t218 * t281 + t144;
t415 = -m(7) * t73 + t154;
t414 = t216 * t432;
t349 = t218 * t243;
t328 = mrSges(7,2) * t349;
t269 = -t216 * mrSges(7,1) - t328;
t331 = mrSges(6,3) * t349;
t270 = t216 * mrSges(6,1) + t331;
t397 = -t243 / 0.2e1;
t412 = t269 * t396 + t270 * t397 + t435 * t399;
t411 = t218 ^ 2;
t409 = m(6) / 0.2e1;
t408 = -m(7) / 0.2e1;
t407 = m(7) / 0.2e1;
t406 = mrSges(7,1) / 0.2e1;
t123 = t242 * t133;
t132 = t243 * t419;
t63 = t132 + t123;
t405 = t63 / 0.2e1;
t404 = t154 / 0.2e1;
t403 = -t216 / 0.2e1;
t398 = t242 / 0.2e1;
t392 = m(7) * t218;
t391 = m(7) * t281;
t388 = t53 * mrSges(7,1);
t387 = t68 * mrSges(6,1);
t386 = t69 * mrSges(6,2);
t353 = t216 * qJ(6);
t47 = t63 + t353;
t384 = t47 - t63;
t382 = m(7) * qJD(6);
t381 = mrSges(5,3) * t218;
t374 = t218 * mrSges(7,1);
t373 = t242 * mrSges(6,1);
t367 = t243 * mrSges(6,2);
t266 = t218 * t220;
t267 = t218 * t221;
t335 = t411 / 0.2e1;
t7 = -((-t383 - t390) * t243 + (-t353 + t384) * t242) * t392 / 0.2e1 + (t266 + t267) * t403 + t412 * t218 + t432 * t434 * t335;
t362 = t7 * qJD(1);
t211 = t216 * mrSges(5,2);
t233 = -t322 - pkin(4);
t348 = t231 * t216;
t341 = t434 * t348;
t213 = -t282 + t233;
t354 = t213 * t218;
t260 = (-t216 * t359 + t218 * t360) * pkin(3) * m(5) / 0.2e1 + (-t341 - t354) * t407 + (-t233 * t218 - t341) * t409;
t247 = t260 + (-t240 / 0.2e1 - t239 / 0.2e1) * t414;
t248 = -m(6) * (t242 * t69 + t243 * t68) / 0.2e1 + (t242 * t52 - t243 * t53) * t408 - t437 / 0.2e1;
t304 = t430 - t221 / 0.2e1;
t351 = t216 * t243;
t158 = -t218 * mrSges(6,1) + mrSges(6,3) * t351;
t330 = mrSges(7,2) * t351;
t159 = -t330 + t374;
t305 = t159 / 0.2e1 - t158 / 0.2e1;
t306 = -t157 / 0.2e1 - t161 / 0.2e1;
t9 = t211 + t305 * t243 + t306 * t242 + (mrSges(5,1) + t304) * t218 + t247 + t248;
t361 = t9 * qJD(1);
t22 = m(7) * (t47 * t216 + t349 * t73) + t216 * t160 - t154 * t349;
t357 = qJD(1) * t22;
t344 = t243 * t160;
t250 = (t242 * t383 + t243 * t384) * t408 - t344 / 0.2e1;
t251 = t281 * t407 - t366 / 0.2e1;
t292 = t367 + t373;
t312 = t349 / 0.2e1;
t297 = mrSges(7,2) * t312;
t11 = t242 * t297 + (t251 + t292 + t372) * t216 + t250;
t356 = t11 * qJD(1);
t298 = (m(7) / 0.4e1 + m(6) / 0.4e1) * (-0.1e1 + t434) * t218 * t216;
t24 = 0.2e1 * t298;
t345 = t24 * qJD(1);
t339 = qJD(5) * t242;
t338 = qJD(5) * t243;
t148 = m(7) * t352;
t337 = t148 * qJD(1);
t336 = t148 * qJD(6);
t334 = t53 * t407;
t327 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t326 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t319 = -t352 / 0.2e1;
t318 = t352 / 0.2e1;
t317 = -t351 / 0.2e1;
t299 = -pkin(5) * mrSges(7,2) + t429;
t294 = qJ(6) * mrSges(7,2) + t436;
t284 = Ifges(7,3) * t242 + t378;
t111 = -Ifges(7,6) * t218 - t216 * t284;
t112 = -Ifges(6,6) * t218 + t216 * t287;
t113 = -Ifges(7,4) * t218 - t216 * t301;
t114 = -Ifges(6,5) * t218 - t216 * t290;
t152 = t216 * t291;
t153 = t216 * t292;
t155 = t292 * t218;
t264 = mrSges(4,1) * t393 + mrSges(4,2) * t394;
t72 = -t216 * t281 + t419;
t1 = t324 * t264 + t63 * t157 + t62 * t158 + t48 * t159 + t52 * t160 + t47 * t161 - t73 * t152 - t144 * t153 - t72 * t154 - t419 * t155 + m(6) * (t144 * t419 + t62 * t68 + t63 * t69) + m(7) * (t47 * t52 + t48 * t53 + t72 * t73) + (-mrSges(5,2) * t333 - Ifges(5,4) * t218 + (t69 * mrSges(6,3) - t111 / 0.2e1 + t112 / 0.2e1 + t326 * t218) * t242 + (-t53 * mrSges(7,2) + t68 * mrSges(6,3) - t113 / 0.2e1 - t114 / 0.2e1 + t327 * t218) * t243) * t218 + (mrSges(5,1) * t333 - t388 + t387 - t386 + (Ifges(5,1) - Ifges(5,2) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t240 + ((Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t242 + (Ifges(7,5) - Ifges(6,4)) * t243) * t242 - t428) * t218 + (Ifges(5,4) - t421) * t216) * t216 + (-Ifges(4,2) + Ifges(4,1)) * t394 * t393 + (-t393 ^ 2 + t394 ^ 2) * Ifges(4,4) + (-t218 * mrSges(5,1) - t211 + t437) * t219;
t276 = t242 * t62 - t243 * t63;
t278 = -t242 * t48 - t243 * t47;
t6 = (t404 + t155 / 0.2e1 + t306 * t243 - t305 * t242 + (-t277 - t73) * t407 + (-t144 - t275) * t409) * t218 + (-t152 / 0.2e1 - t153 / 0.2e1 + (-t160 / 0.2e1 + t376 / 0.2e1) * t243 + (t297 + (mrSges(6,1) / 0.2e1 + t406) * t216) * t242 + (t278 + t72) * t407 + (t419 + t276) * t409) * t216;
t280 = t1 * qJD(1) + t6 * qJD(2);
t253 = Ifges(7,6) * t216 - t218 * t284;
t254 = Ifges(6,6) * t216 + t218 * t287;
t265 = t282 * t218;
t5 = t253 * t312 - t254 * t349 / 0.2e1 - t48 * t329 - t47 * t328 - t73 * t266 - t144 * t267 + t418 * t411 + (t242 * t425 + t243 * t424) * t335 - t423 * t350 / 0.2e1 - t415 * t265 + t422 * t216 + (-m(7) * t48 - t269 + t270 - t331) * t63 + (-m(7) * t47 + t332 - t435) * t62;
t279 = -t5 * qJD(1) - t7 * qJD(2);
t23 = 0.4e1 * t298;
t274 = t6 * qJD(1) + t23 * qJD(2);
t8 = (-m(6) * t144 + t155 + t381 + t415) * t218 + (-t344 + (mrSges(5,3) + t367) * t216 + (t216 * t385 + t328) * t242 + m(7) * t278 + m(6) * t276) * t216 + (-t144 * t218 - t216 * t419) * m(5);
t273 = qJD(1) * t8 + qJD(2) * t24;
t162 = (m(7) * t213 + t220) * t242;
t249 = (-t242 * t73 + (t348 + t354) * t243) * t407 - t154 * t399;
t19 = -t330 + (t220 * t397 + t406) * t218 + t334 - t249;
t272 = qJD(1) * t19 + qJD(3) * t162;
t27 = t375 + 0.2e1 * (t353 / 0.2e1 + t123 / 0.4e1 + t132 / 0.4e1 - t63 / 0.4e1) * m(7);
t271 = qJD(1) * t27 + qJD(5) * t234;
t244 = t224 * t312 + (-t213 * t265 + t281 * t73) * t407 - t281 * t404 + t144 * t292 / 0.2e1 + t73 * t291 / 0.2e1 + t233 * t267 / 0.2e1 + t213 * t266 / 0.2e1 + t265 * t430 + t421 * t216 / 0.4e1 + t423 * t243 / 0.4e1 + (-t284 / 0.4e1 + t225 / 0.2e1) * t350 + (t253 / 0.4e1 - t254 / 0.4e1) * t242 + t420 * mrSges(6,3) + (t424 + t426) * t350 / 0.4e1 - (t425 + t416) * t349 / 0.4e1 + ((-t242 * t384 + t243 * t383) * t407 + t412) * t231 + ((-t47 / 0.2e1 + t405) * t242 + t383 * t396 + t420) * mrSges(7,2);
t245 = (-pkin(5) * t53 + qJ(6) * t52) * t407 - pkin(5) * t159 / 0.2e1 + qJ(6) * t161 / 0.2e1 + t52 * mrSges(7,3) / 0.2e1 - t388 / 0.2e1 + t387 / 0.2e1 - t386 / 0.2e1;
t3 = Ifges(6,6) * t318 + Ifges(7,6) * t319 + t317 * t429 + t400 * t428 - t244 + t245;
t31 = t284 * t397 + t233 * t292 + t281 * t220 + (t291 + t391) * t213 + t426 * t396 + t416 * t398 + t418;
t33 = (t389 / 0.2e1 - t358 / 0.2e1 - t281 / 0.2e1) * m(7) * t216;
t261 = t3 * qJD(1) + t33 * qJD(2) - t31 * qJD(3);
t246 = (t367 / 0.2e1 + t373 / 0.2e1 + t372 / 0.2e1 + t251) * t216;
t210 = (m(7) * t231 + mrSges(7,2)) * t243;
t34 = -t391 * t403 + t246 + (t291 + t292) * t216 / 0.2e1;
t25 = (t63 + 0.2e1 * t353) * t407 + m(7) * t405 + t160;
t20 = t220 * t312 + t374 / 0.2e1 + t334 + t249;
t12 = t268 * t396 + t269 * t398 + t270 * t399 + t246 - t250;
t10 = t158 * t396 + t159 * t397 + t304 * t218 + t398 * t427 + t247 - t248;
t4 = t244 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t218 + (-t242 * t326 - t243 * t327) * t216 + t245;
t2 = qJD(3) * t6 + qJD(4) * t24 - qJD(5) * t7;
t13 = [qJD(3) * t1 + qJD(4) * t8 - qJD(5) * t5 + qJD(6) * t22, t2 (Ifges(4,5) * t394 - Ifges(4,6) * t393 - mrSges(4,1) * t323 + t277 * mrSges(7,2) + t275 * mrSges(6,3) + t112 * t396 + t111 * t397 + (t114 + t113) * t398 + (t424 + t225) * t317 + t320 * t381 + mrSges(4,2) * t321 + (mrSges(5,3) * t322 - Ifges(5,5)) * t216 + t422 - t233 * t153 + t72 * t220 + Ifges(5,6) * t218 + t224 * t318 + t425 * t319 + (m(7) * t72 - t152) * t213 + (-m(5) * t320 + mrSges(5,2)) * t144 + (-m(5) * t322 + m(6) * t233 - mrSges(5,1) + t221) * t419 + ((-t158 + t159) * t242 + t427 * t243 + m(6) * t275 + m(7) * t277) * t231) * qJD(3) + t10 * qJD(4) + t4 * qJD(5) + t20 * qJD(6) + t280, qJD(3) * t10 + qJD(5) * t12 + t273, t4 * qJD(3) + t12 * qJD(4) + t25 * qJD(6) + t279 + ((-m(7) * pkin(5) - t385) * t63 + (-mrSges(6,2) + t234) * t62 + (t242 * t299 + t243 * t294) * t218) * qJD(5), qJD(3) * t20 + qJD(5) * t25 + t357; t2, t23 * qJD(3), t34 * qJD(5) + t274 - t336 + (t211 - t264 + (mrSges(5,1) + t252) * t218 + 0.2e1 * t260 - t434 * t414) * qJD(3), t345, -t362 + t34 * qJD(3) + (-t243 * t382 + t431) * t218, -t148 * qJD(3) - t338 * t392; qJD(4) * t9 - qJD(5) * t3 - qJD(6) * t19 - t280, -qJD(5) * t33 - t274, qJD(5) * t31 - qJD(6) * t162, t361, t210 * qJD(6) - t231 * t431 - t294 * t339 + t299 * t338 - t261, qJD(5) * t210 - t272; -qJD(3) * t9 - qJD(5) * t11 - t273 + t336, -t345, -t361, 0, -t356 + (t366 - t367 - t391) * qJD(5) + (-qJD(5) * t385 + t382) * t242, m(7) * t339 + t337; qJD(3) * t3 + qJD(4) * t11 + qJD(6) * t27 - t279, qJD(3) * t33 + t362, t261, t356, t234 * qJD(6), t271; qJD(3) * t19 - qJD(4) * t148 - qJD(5) * t27 - t357, 0, t272, -t337, -t271, 0;];
Cq  = t13;
