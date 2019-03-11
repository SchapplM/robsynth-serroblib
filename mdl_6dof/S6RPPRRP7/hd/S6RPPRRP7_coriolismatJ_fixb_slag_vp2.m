% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:45
% EndTime: 2019-03-09 02:12:53
% DurationCPUTime: 5.44s
% Computational Cost: add. (10086->431), mult. (18954->584), div. (0->0), fcn. (19696->6), ass. (0->228)
t266 = sin(qJ(5));
t395 = Ifges(7,4) + Ifges(6,4);
t461 = t395 * t266;
t268 = cos(qJ(5));
t460 = t395 * t268;
t458 = Ifges(6,1) + Ifges(7,1);
t457 = Ifges(6,2) + Ifges(7,2);
t459 = -m(7) / 0.2e1;
t447 = Ifges(6,5) + Ifges(7,5);
t393 = -Ifges(7,6) - Ifges(6,6);
t264 = sin(pkin(9));
t265 = cos(pkin(9));
t267 = sin(qJ(4));
t409 = cos(qJ(4));
t236 = t264 * t267 - t265 * t409;
t451 = t236 / 0.2e1;
t415 = -t266 / 0.2e1;
t411 = t268 / 0.2e1;
t441 = -t457 * t266 + t460;
t439 = t458 * t268 - t461;
t391 = -qJ(6) - pkin(8);
t244 = t391 * t266;
t247 = t391 * t268;
t292 = t244 * t268 - t247 * t266;
t237 = t264 * t409 + t267 * t265;
t396 = t237 * pkin(5);
t361 = t268 * t236;
t248 = t264 * pkin(3) + qJ(2);
t397 = pkin(8) * t236;
t402 = pkin(4) * t237;
t146 = t248 + t397 + t402;
t392 = pkin(1) + qJ(3);
t351 = -pkin(7) - t392;
t243 = t351 * t264;
t314 = t265 * t351;
t169 = t243 * t409 + t267 * t314;
t65 = t268 * t146 - t169 * t266;
t53 = qJ(6) * t361 + t65;
t44 = t53 + t396;
t363 = t266 * t236;
t66 = t146 * t266 + t169 * t268;
t54 = qJ(6) * t363 + t66;
t296 = t266 * t44 - t268 * t54;
t347 = mrSges(7,3) * t363;
t156 = -mrSges(7,2) * t237 + t347;
t160 = mrSges(7,1) * t237 + mrSges(7,3) * t361;
t359 = t156 * t411 + t160 * t415;
t453 = (t236 * t292 - t296) * t459 - t359;
t430 = m(7) * pkin(5);
t422 = t160 / 0.2e1;
t448 = mrSges(6,3) + mrSges(7,3);
t446 = -mrSges(7,1) - t430;
t445 = t393 * t237;
t444 = -t244 * t266 - t247 * t268;
t258 = t266 * mrSges(7,1);
t259 = t268 * mrSges(7,2);
t355 = t259 + t258;
t354 = t264 ^ 2 + t265 ^ 2;
t262 = t266 ^ 2;
t263 = t268 ^ 2;
t353 = t262 + t263;
t443 = t393 * t266 + t268 * t447;
t442 = t268 * t457 + t461;
t440 = t266 * t458 + t460;
t338 = -t53 / 0.2e1 + t44 / 0.2e1;
t342 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t437 = t236 * t342 + t447 * t451 + t439 * t237 / 0.2e1;
t184 = t236 ^ 2;
t436 = t237 ^ 2;
t435 = 2 * qJD(4);
t434 = -m(4) / 0.2e1;
t433 = -m(6) / 0.2e1;
t432 = m(6) / 0.2e1;
t431 = m(7) / 0.2e1;
t429 = mrSges(6,1) / 0.2e1;
t428 = -mrSges(7,1) / 0.2e1;
t427 = -mrSges(7,2) / 0.2e1;
t367 = t237 * t268;
t168 = t243 * t267 - t409 * t314;
t403 = pkin(4) * t236;
t170 = pkin(8) * t237 - t403;
t73 = t168 * t266 + t268 * t170;
t47 = -pkin(5) * t236 + qJ(6) * t367 + t73;
t426 = t47 / 0.2e1;
t425 = pkin(4) * mrSges(6,1);
t424 = pkin(4) * mrSges(6,2);
t149 = t355 * t236;
t423 = -t149 / 0.2e1;
t421 = -t168 / 0.2e1;
t218 = mrSges(7,1) * t361;
t420 = -t218 / 0.2e1;
t417 = -t244 / 0.2e1;
t400 = pkin(5) * t268;
t251 = -pkin(4) - t400;
t416 = t251 / 0.2e1;
t414 = -t266 / 0.4e1;
t413 = t266 / 0.2e1;
t412 = -t268 / 0.2e1;
t410 = t268 / 0.4e1;
t119 = -pkin(5) * t363 + t168;
t408 = m(7) * t119;
t406 = m(7) * t236;
t405 = m(7) * t237;
t404 = m(7) * t251;
t401 = pkin(5) * t266;
t348 = mrSges(6,3) * t363;
t382 = t237 * mrSges(6,2);
t157 = t348 - t382;
t399 = pkin(8) * t157;
t383 = t237 * mrSges(6,1);
t161 = mrSges(6,3) * t361 + t383;
t398 = pkin(8) * t161;
t390 = -t44 + t53;
t376 = t268 * mrSges(6,2);
t308 = t266 * mrSges(6,1) + t376;
t150 = t236 * t308;
t295 = t266 * t65 - t268 * t66;
t357 = t160 + t161;
t358 = t156 + t157;
t369 = t168 * t236;
t9 = (mrSges(5,3) * t236 + t149 + t150) * t236 + (mrSges(5,3) * t237 + t266 * t357 - t268 * t358) * t237 + m(6) * (t237 * t295 - t369) + m(7) * (-t119 * t236 + t237 * t296) + m(5) * (-t169 * t237 - t369) + (m(4) * t392 + mrSges(4,3)) * t354;
t385 = qJD(1) * t9;
t368 = t237 * t266;
t120 = -pkin(5) * t368 + t169;
t147 = t355 * t237;
t148 = t308 * t237;
t154 = t236 * mrSges(7,2) + mrSges(7,3) * t368;
t155 = t236 * mrSges(6,2) + mrSges(6,3) * t368;
t158 = -t236 * mrSges(7,1) + mrSges(7,3) * t367;
t159 = -t236 * mrSges(6,1) + mrSges(6,3) * t367;
t223 = t237 * mrSges(5,2);
t340 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t281 = t340 * t236 + t393 * t451 - t441 * t237 / 0.2e1;
t341 = Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1;
t343 = Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1;
t74 = -t268 * t168 + t266 * t170;
t59 = qJ(6) * t368 + t74;
t1 = -t119 * t147 - t120 * t149 - t168 * t148 - t169 * t150 + t54 * t154 + t66 * t155 + t59 * t156 + t74 * t157 + t44 * t158 + t65 * t159 + t47 * t160 + t73 * t161 - t248 * t223 + m(6) * (t168 * t169 + t65 * t73 + t66 * t74) + m(7) * (t119 * t120 + t44 * t47 + t54 * t59) + (-t248 * mrSges(5,1) - Ifges(5,4) * t236 + t281 * t266 + t437 * t268) * t236 + ((-Ifges(7,3) - Ifges(6,3) + Ifges(5,1) - Ifges(5,2) + t343 * t263 + (t266 * t341 - t460) * t266) * t236 + (Ifges(5,4) - t443) * t237) * t237;
t384 = t1 * qJD(1);
t381 = t247 * mrSges(7,3);
t380 = t266 * mrSges(6,2);
t379 = t266 * mrSges(7,2);
t378 = t268 * mrSges(6,1);
t377 = t268 * mrSges(7,1);
t375 = t268 * t44;
t246 = -t378 + t380;
t285 = t236 * t246;
t311 = mrSges(7,2) * t363 - t218;
t317 = t457 - t458;
t318 = t447 * t237;
t337 = m(7) * t390;
t349 = pkin(5) * t361;
t4 = t66 * t161 - t168 * t285 - t149 * t349 - t119 * t311 - t53 * t156 + t54 * t160 + t44 * t347 - t54 * t337 + ((-t236 * t461 - t318) * t266 + (-t66 * mrSges(6,3) - t54 * mrSges(7,3) + pkin(5) * t408 - t317 * t363 + t361 * t395 + t445) * t268) * t236 + (-t157 + t348) * t65;
t374 = t4 * qJD(1);
t323 = t263 / 0.2e1 + t262 / 0.2e1;
t280 = t448 * t323;
t326 = t161 / 0.2e1 + t422;
t328 = -t157 / 0.2e1 - t156 / 0.2e1;
t344 = mrSges(6,2) / 0.2e1 + mrSges(7,2) / 0.2e1;
t7 = (t237 * t328 + t344) * t266 + (t237 * t280 + t344 * t363 + t420) * t236 + (t428 + (-m(7) * t338 - t326) * t237 + (mrSges(6,1) + t430) * (-t184 / 0.2e1 - 0.1e1 / 0.2e1)) * t268;
t373 = t7 * qJD(1);
t310 = t323 * mrSges(6,3);
t321 = pkin(8) * t353;
t272 = (-mrSges(7,3) * t323 - t310) * t237 + (-t237 * t321 + t403) * t432 + (-t236 * t251 - t237 * t444) * t431;
t274 = (t266 * t74 + t268 * t73) * t433 + (t266 * t59 + t268 * t47) * t459;
t245 = -t377 + t379;
t325 = -t246 / 0.2e1 - t245 / 0.2e1;
t327 = -t159 / 0.2e1 - t158 / 0.2e1;
t329 = t155 / 0.2e1 + t154 / 0.2e1;
t12 = t223 + t327 * t268 - t329 * t266 + (mrSges(5,1) + t325) * t236 + t272 + t274;
t372 = qJD(1) * t12;
t224 = t237 * mrSges(5,1);
t287 = m(7) * (t266 * t54 + t375);
t15 = t264 * mrSges(4,1) + t265 * mrSges(4,2) - t236 * mrSges(5,2) + mrSges(3,3) + t224 + t357 * t268 + t358 * t266 + (m(4) + m(3)) * qJ(2) + m(6) * (t266 * t66 + t268 * t65) + t287 + m(5) * t248;
t371 = qJD(1) * t15;
t23 = (t266 * t156 + t268 * t160 + t287) * t236;
t370 = qJD(1) * t23;
t350 = m(7) / 0.4e1 + m(6) / 0.4e1;
t275 = m(5) * (-t184 - t436) / 0.2e1 + t354 * t434 + 0.2e1 * t350 * (-t353 * t436 - t184);
t322 = m(7) * t353;
t279 = t434 - m(5) / 0.2e1 + t353 * t433 - t322 / 0.2e1;
t26 = t275 + t279;
t365 = t26 * qJD(1);
t324 = -t262 / 0.4e1 - t263 / 0.4e1;
t79 = 0.2e1 * (0.1e1 / 0.4e1 - t324) * t406;
t360 = t79 * qJD(1);
t356 = (t259 / 0.2e1 + t258 / 0.2e1) * t237;
t345 = -t400 / 0.2e1;
t339 = -Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1;
t226 = -m(7) * t401 - t355;
t334 = t363 / 0.2e1;
t315 = -mrSges(7,3) * pkin(5) + t447;
t313 = t429 + t430 / 0.2e1;
t312 = t322 / 0.2e1;
t309 = t120 * t431 + t367 * t427 + t368 * t428;
t294 = -t266 * t73 + t268 * t74;
t28 = 0.4e1 * t350 * (0.1e1 - t353) * t237 * t236;
t271 = (-t148 / 0.2e1 - t147 / 0.2e1 + t328 * t268 + t326 * t266 + (t169 + t295) * t432 + (t120 + t296) * t431) * t236 + (-t150 / 0.2e1 + t423 + t329 * t268 + t327 * t266 + (t168 + t294) * t432 + (-t266 * t47 + t268 * t59 + t119) * t431) * t237;
t286 = m(7) * t292;
t6 = -t286 / 0.2e1 + t271;
t293 = t6 * qJD(1) + t28 * qJD(2);
t88 = m(7) * t349 - t311;
t290 = qJD(1) * t88 + qJD(4) * t226;
t289 = -t337 / 0.2e1 + t422;
t288 = -t341 + t343;
t10 = (t382 / 0.2e1 + t328) * t268 + (t383 / 0.2e1 + (t396 / 0.2e1 + t338) * m(7) + t326) * t266 + t356;
t284 = t10 * qJD(1);
t283 = -t119 * t355 / 0.2e1 + t156 * t417;
t130 = m(7) * t444 + t353 * mrSges(7,3);
t20 = t309 + t453;
t80 = 0.2e1 * (0.1e1 / 0.4e1 + t324) * t405;
t282 = -qJD(1) * t20 - qJD(2) * t80 + qJD(4) * t130;
t273 = (m(7) * t426 + t158 / 0.2e1) * pkin(5) + mrSges(7,1) * t426 + t59 * t427 + t73 * t429 - t74 * mrSges(6,2) / 0.2e1;
t2 = t218 * t416 - t289 * t247 + (-pkin(8) * t310 + t339) * t236 + (t399 / 0.2e1 + mrSges(6,1) * t421 - t445 + (t149 / 0.2e1 - t408 / 0.2e1) * pkin(5) + (t244 * mrSges(7,3) / 0.2e1 + t251 * t427 + t424 / 0.2e1 - t288 * t266) * t236) * t266 + (t398 / 0.2e1 + mrSges(6,2) * t421 - t318 + t338 * mrSges(7,3) + (t381 / 0.2e1 - t425 / 0.2e1 + t288 * t268 - 0.2e1 * t461 + (t404 / 0.2e1 + t245 / 0.2e1) * pkin(5)) * t236) * t268 + t273 + t283;
t29 = -t251 * t355 - t245 * t401 + (-t460 + t424) * t268 + (-pkin(5) * t404 + t268 * t317 + t425 + t461) * t266;
t278 = t2 * qJD(1) + t29 * qJD(4);
t277 = -t378 / 0.2e1 + t380 / 0.2e1 + t379 / 0.2e1 + m(7) * t345;
t81 = t237 * t312 + t405 / 0.2e1;
t78 = -t406 / 0.2e1 + t236 * t312;
t31 = t334 * t430 + (t344 * t268 + (mrSges(7,1) / 0.2e1 + t313) * t266) * t236 + (t308 + t355) * t451;
t25 = t275 - t279;
t21 = t309 - t453;
t13 = t325 * t236 + t272 - t274 + (t154 + t155) * t413 + (t158 + t159) * t411;
t11 = t157 * t411 + t161 * t415 + t337 * t413 + (t376 / 0.2e1 + t313 * t266) * t237 + t356 + t359;
t8 = t377 / 0.2e1 + (t236 * t277 + t420) * t236 + (t328 * t266 + (-t161 / 0.2e1 - t289) * t268 + t280 * t236) * t237 - t277;
t5 = t286 / 0.2e1 + t271;
t3 = -pkin(4) * t285 / 0.2e1 + (-t390 * t247 + (t119 * t266 - t251 * t361) * pkin(5)) * t431 + t398 * t412 + t399 * t415 + t311 * t416 + t347 * t417 + t401 * t423 - t283 + t273 + t168 * t308 / 0.2e1 + t247 * t422 + t440 * t334 + t441 * t363 / 0.4e1 - t439 * t361 / 0.4e1 + t353 * mrSges(6,3) * t397 / 0.2e1 + (t245 * t345 - t410 * t439 - t414 * t441 + t339) * t236 + (t53 * t411 - t375 / 0.2e1) * mrSges(7,3) + (-t381 + t442) * t361 / 0.2e1 + (t443 / 0.4e1 - t340 * t266 - t342 * t268 - t393 * t414 + t447 * t410) * t237;
t14 = [qJD(2) * t15 + qJD(3) * t9 + qJD(4) * t1 - qJD(5) * t4 + qJD(6) * t23, qJD(3) * t25 + qJD(4) * t5 + qJD(5) * t8 + t371, qJD(2) * t25 + qJD(4) * t13 + qJD(5) * t11 + qJD(6) * t78 + t385, t384 + t5 * qJD(2) + t13 * qJD(3) + t3 * qJD(5) + t21 * qJD(6) + ((-pkin(4) * t169 + pkin(8) * t294) * t432 + (t120 * t251 + t244 * t47 - t247 * t59) * t431) * t435 + (t168 * mrSges(5,2) + Ifges(5,6) * t236 + pkin(4) * t148 + t120 * t245 - t251 * t147 - t247 * t154 + t244 * t158 + (t74 * mrSges(6,3) + t59 * mrSges(7,3) + pkin(8) * t155 + t281) * t268 + (-t73 * mrSges(6,3) - t47 * mrSges(7,3) - pkin(8) * t159 - t437) * t266 + (t412 * t440 + t413 * t442 - Ifges(5,5)) * t237 + (t246 - mrSges(5,1)) * t169) * qJD(4), t8 * qJD(2) + t11 * qJD(3) + t3 * qJD(4) - t374 + (-mrSges(6,1) * t66 - mrSges(6,2) * t65 - mrSges(7,2) * t53 + (t266 * t315 - t268 * t393) * t236 + t446 * t54) * qJD(5), qJD(3) * t78 + qJD(4) * t21 + t370; qJD(3) * t26 + qJD(4) * t6 + qJD(5) * t7 - t371, t28 * qJD(4), t365, t31 * qJD(5) + t81 * qJD(6) + ((-t236 * t321 - t402) * t432 + (-t236 * t444 + t237 * t251) * t431) * t435 + t293 + (-t224 + (-t353 * t448 + mrSges(5,2)) * t236 + (t245 + t246) * t237) * qJD(4), t373 + t31 * qJD(4) + ((mrSges(6,2) + mrSges(7,2)) * t266 + (-mrSges(6,1) + t446) * t268) * qJD(5) * t237, t81 * qJD(4); -qJD(2) * t26 - qJD(4) * t12 - qJD(5) * t10 + qJD(6) * t79 - t385, -t365, 0, -t372 (-t308 + t226) * qJD(5) - t284, t360; -qJD(2) * t6 + qJD(3) * t12 - qJD(5) * t2 - qJD(6) * t20 - t384, -qJD(6) * t80 - t293, t372, -qJD(5) * t29 + qJD(6) * t130, -t278 + (-mrSges(7,2) * t244 + (mrSges(6,2) * pkin(8) + t393) * t266 + (-mrSges(6,1) * pkin(8) + t315) * t268 - t446 * t247) * qJD(5), t282; -qJD(2) * t7 + qJD(3) * t10 + qJD(4) * t2 + qJD(6) * t88 + t374, -t373, t284, t226 * qJD(6) + t278, 0, t290; -qJD(3) * t79 + qJD(4) * t20 - qJD(5) * t88 - t370, t80 * qJD(4), -t360, -qJD(5) * t226 - t282, -t290, 0;];
Cq  = t14;
