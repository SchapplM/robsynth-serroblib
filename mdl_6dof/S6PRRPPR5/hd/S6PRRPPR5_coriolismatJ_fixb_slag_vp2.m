% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:14
% EndTime: 2019-03-08 21:18:23
% DurationCPUTime: 5.57s
% Computational Cost: add. (9261->495), mult. (20644->700), div. (0->0), fcn. (21270->10), ass. (0->262)
t283 = sin(qJ(3));
t285 = cos(qJ(3));
t433 = -t283 * pkin(3) + qJ(4) * t285;
t442 = m(5) * t433;
t279 = sin(pkin(6));
t286 = cos(qJ(2));
t358 = t279 * t286;
t340 = -t358 / 0.2e1;
t441 = t340 * t442;
t417 = pkin(4) + pkin(8);
t284 = sin(qJ(2));
t360 = t279 * t284;
t440 = t283 ^ 2 + t285 ^ 2;
t377 = cos(pkin(6));
t215 = t283 * t360 - t285 * t377;
t278 = sin(pkin(11));
t280 = cos(pkin(11));
t165 = t215 * t280 + t278 * t358;
t166 = t215 * t278 - t280 * t358;
t282 = sin(qJ(6));
t401 = cos(qJ(6));
t81 = t165 * t401 - t282 * t166;
t439 = t81 / 0.2e1;
t438 = -t285 / 0.2e1;
t437 = -Ifges(4,4) - Ifges(5,6);
t275 = t280 ^ 2;
t349 = t278 ^ 2 + t275;
t436 = t349 * mrSges(6,3);
t251 = t285 * mrSges(5,2) - t283 * mrSges(5,3);
t435 = -t285 * mrSges(4,1) + t283 * mrSges(4,2) + t251;
t311 = -t282 * t278 + t401 * t280;
t353 = t283 * t286;
t175 = (-t278 * t284 + t280 * t353) * t279;
t357 = t280 * t175;
t176 = (t278 * t353 + t280 * t284) * t279;
t362 = t278 * t176;
t321 = t357 + t362;
t107 = t282 * t175 + t176 * t401;
t312 = t278 * t401 + t282 * t280;
t363 = t312 * t107;
t106 = t175 * t401 - t282 * t176;
t364 = t311 * t106;
t434 = t363 + t364;
t225 = qJ(5) * t283 - t433;
t256 = t417 * t285;
t234 = t280 * t256;
t124 = pkin(5) * t285 + t234 + (-pkin(9) * t283 - t225) * t278;
t154 = t280 * t225 + t278 * t256;
t356 = t280 * t283;
t138 = pkin(9) * t356 + t154;
t62 = t124 * t401 - t282 * t138;
t63 = t282 * t124 + t138 * t401;
t432 = t311 * t62 + t312 * t63;
t254 = t283 * mrSges(4,1) + t285 * mrSges(4,2);
t158 = mrSges(7,1) * t312 + mrSges(7,2) * t311;
t248 = t278 * mrSges(6,1) + t280 * mrSges(6,2);
t431 = mrSges(5,3) + t158 + t248;
t359 = t279 * t285;
t216 = t283 * t377 + t284 * t359;
t430 = 0.2e1 * t216;
t429 = 2 * qJD(3);
t428 = m(5) / 0.2e1;
t427 = -m(6) / 0.2e1;
t426 = -m(6) / 0.4e1;
t425 = m(6) / 0.2e1;
t424 = -m(7) / 0.2e1;
t423 = m(7) / 0.2e1;
t422 = -mrSges(7,1) / 0.2e1;
t421 = -mrSges(7,2) / 0.2e1;
t419 = mrSges(7,3) / 0.2e1;
t82 = t282 * t165 + t166 * t401;
t418 = t82 / 0.2e1;
t118 = t311 * t216;
t416 = t118 / 0.2e1;
t202 = t311 * t285;
t171 = -mrSges(7,2) * t283 - t202 * mrSges(7,3);
t415 = t171 / 0.2e1;
t203 = t312 * t285;
t173 = mrSges(7,1) * t283 + t203 * mrSges(7,3);
t414 = -t173 / 0.2e1;
t398 = Ifges(6,4) * t278;
t413 = Ifges(6,6) * t438 - (t280 * Ifges(6,2) + t398) * t283 / 0.2e1;
t201 = t311 * t283;
t412 = t201 / 0.2e1;
t411 = -t202 / 0.2e1;
t410 = -t203 / 0.2e1;
t204 = t312 * t283;
t409 = t204 / 0.2e1;
t408 = t216 / 0.2e1;
t407 = t311 / 0.2e1;
t406 = -t311 / 0.2e1;
t405 = -t312 / 0.2e1;
t404 = -t278 / 0.2e1;
t403 = t278 / 0.2e1;
t402 = t280 / 0.2e1;
t281 = -pkin(3) - qJ(5);
t400 = -pkin(9) + t281;
t399 = Ifges(6,1) * t278;
t397 = Ifges(6,4) * t280;
t396 = Ifges(7,4) * t203;
t395 = Ifges(7,4) * t311;
t394 = t201 * mrSges(7,1);
t393 = t202 * mrSges(7,1);
t392 = t203 * mrSges(7,2);
t391 = t204 * mrSges(7,2);
t390 = t278 * mrSges(6,2);
t389 = t278 * Ifges(6,5);
t388 = t280 * mrSges(6,1);
t387 = t280 * Ifges(6,6);
t385 = t283 * mrSges(6,1);
t384 = t283 * mrSges(5,2);
t383 = t283 * mrSges(6,2);
t382 = t285 * mrSges(6,1);
t380 = t285 * mrSges(6,2);
t376 = qJ(4) * t215;
t375 = qJ(4) * t283;
t373 = t165 * t278;
t372 = t166 * t280;
t299 = (-t202 * t406 - t312 * t410) * mrSges(7,3) + t171 * t407 + t173 * t405;
t316 = t201 * t421 + t204 * t422;
t17 = t299 + t316;
t371 = t17 * qJD(2);
t119 = t312 * t216;
t150 = t215 * t216;
t322 = t165 * t280 + t166 * t278;
t19 = m(7) * (t118 * t81 + t119 * t82 - t150) + m(6) * (t216 * t322 - t150);
t370 = t19 * qJD(1);
t343 = t285 * t358;
t174 = t216 * t343;
t314 = (t215 * t283 - t360) * t279;
t20 = m(7) * (t106 * t81 + t107 * t82 + t174) + m(6) * (t165 * t175 + t166 * t176 + t174) + m(5) * (t216 * t359 + t314) * t286 + (t286 * t314 + t174) * m(4);
t369 = t20 * qJD(1);
t368 = t201 * t312;
t367 = t202 * t312;
t366 = t203 * t311;
t365 = t204 * t311;
t361 = t278 * t283;
t355 = t280 * t285;
t352 = -Ifges(7,5) * t202 + Ifges(7,6) * t203;
t224 = t281 * t285 - pkin(2) - t375;
t255 = t417 * t283;
t152 = t280 * t224 + t278 * t255;
t351 = -Ifges(7,5) * t312 - Ifges(7,6) * t311;
t350 = t440 * pkin(8) * t358;
t347 = -m(7) / 0.4e1 + t426;
t346 = t423 + t425;
t344 = (t366 - t367) * t423;
t214 = pkin(5) * t355 + t256;
t157 = mrSges(7,1) * t311 - mrSges(7,2) * t312;
t341 = t157 * t408;
t339 = t358 / 0.2e1;
t253 = -t285 * mrSges(5,3) - t384;
t338 = -t254 / 0.2e1 - t253 / 0.2e1;
t337 = t311 ^ 2 + t312 ^ 2;
t336 = m(6) * t349;
t335 = t349 * t216;
t332 = 0.2e1 * t347 * t215;
t330 = t388 - t390;
t329 = -pkin(3) * t285 - t375;
t126 = t391 - t394;
t127 = -t392 + t393;
t153 = -t225 * t278 + t234;
t170 = -mrSges(7,2) * t285 + t201 * mrSges(7,3);
t172 = mrSges(7,1) * t285 - t204 * mrSges(7,3);
t213 = (-pkin(5) * t280 - t417) * t283;
t217 = t330 * t283;
t218 = t330 * t285;
t235 = -mrSges(6,3) * t361 + t382;
t237 = mrSges(6,3) * t356 - t380;
t236 = mrSges(6,3) * t278 * t285 + t385;
t238 = -mrSges(6,3) * t355 - t383;
t313 = t236 * t402 + t238 * t403;
t233 = t280 * t255;
t151 = -t224 * t278 + t233;
t324 = -t151 * t280 - t152 * t278;
t123 = pkin(5) * t283 + t233 + (pkin(9) * t285 - t224) * t278;
t137 = -pkin(9) * t355 + t152;
t60 = t123 * t401 - t282 * t137;
t61 = t282 * t123 + t137 * t401;
t289 = (-t127 / 0.2e1 - t218 / 0.2e1) * t215 + (t126 / 0.2e1 - t217 / 0.2e1 + t313) * t216 + (t153 * t165 + t154 * t166 - t215 * t256 + (-t255 - t324) * t216) * t425 + (t118 * t60 + t119 * t61 + t213 * t216 - t214 * t215 + t62 * t81 + t63 * t82) * t423 + t173 * t416 + t119 * t415 + t165 * t235 / 0.2e1 + t166 * t237 / 0.2e1 + t172 * t439 + t170 * t418 - t441;
t239 = t400 * t278;
t240 = t400 * t280;
t155 = -t282 * t239 + t240 * t401;
t156 = t239 * t401 + t282 * t240;
t263 = pkin(5) * t278 + qJ(4);
t291 = (qJ(4) * t343 + t281 * t321) * t427 + (t106 * t155 + t107 * t156 + t263 * t343) * t424 + t441;
t2 = t289 + (t364 / 0.2e1 + t363 / 0.2e1) * mrSges(7,3) + (t357 / 0.2e1 + t362 / 0.2e1) * mrSges(6,3) + ((mrSges(4,1) / 0.2e1 - mrSges(5,2) / 0.2e1) * t283 + (mrSges(4,2) / 0.2e1 - mrSges(5,3) / 0.2e1 - t248 / 0.2e1 - t158 / 0.2e1) * t285 + t338) * t358 + t291;
t114 = Ifges(7,4) * t204 + Ifges(7,2) * t201 + t285 * Ifges(7,6);
t115 = -Ifges(7,2) * t202 + t283 * Ifges(7,6) - t396;
t116 = Ifges(7,1) * t204 + Ifges(7,4) * t201 + t285 * Ifges(7,5);
t198 = Ifges(7,4) * t202;
t117 = -Ifges(7,1) * t203 + t283 * Ifges(7,5) - t198;
t200 = t285 * Ifges(6,5) + (t397 + t399) * t283;
t247 = -pkin(2) + t329;
t319 = Ifges(7,5) * t409 + Ifges(7,6) * t412;
t3 = -t433 * t251 - pkin(2) * t254 - t255 * t218 - t256 * t217 + t151 * t235 + t153 * t236 + t152 * t237 + t154 * t238 + t214 * t126 + t213 * t127 + t114 * t411 + t116 * t410 + t117 * t409 + t115 * t412 + t63 * t171 + t60 * t172 + t62 * t173 + t61 * t170 + (t253 - t442) * t247 + m(6) * (t151 * t153 + t152 * t154 - t255 * t256) + m(7) * (t213 * t214 + t60 * t62 + t61 * t63) + (Ifges(7,5) * t410 + Ifges(7,6) * t411 + t200 * t404 + t280 * t413 + (-t389 / 0.2e1 - t387 / 0.2e1 - t437) * t285) * t285 + ((t387 + t389 + t437) * t283 + (-t275 * Ifges(6,2) / 0.2e1 - Ifges(4,2) + Ifges(7,3) - Ifges(5,3) + Ifges(4,1) + Ifges(5,2) + Ifges(6,3) + (-t397 - t399 / 0.2e1) * t278) * t285 + t319) * t283;
t328 = t2 * qJD(1) + t3 * qJD(2);
t125 = -t203 * mrSges(7,1) - mrSges(7,2) * t202;
t128 = Ifges(7,2) * t203 - t198;
t129 = -Ifges(7,1) * t202 + t396;
t6 = t283 * t352 / 0.2e1 + t60 * t171 - t61 * t173 + t214 * t125 + (t61 * mrSges(7,3) + t115 / 0.2e1 - t129 / 0.2e1) * t203 - (-t60 * mrSges(7,3) + t128 / 0.2e1 + t117 / 0.2e1) * t202;
t295 = (t202 * t439 + t203 * t418) * mrSges(7,3) + t125 * t408 + t81 * t415 + t82 * t414;
t318 = t106 * mrSges(7,1) / 0.2e1 + t107 * t421;
t7 = t295 - t318;
t327 = t7 * qJD(1) + t6 * qJD(2);
t300 = -t280 * t238 + t278 * t236 + m(6) * (t151 * t278 - t152 * t280);
t22 = m(7) * (-t201 * t61 + t204 * t60) - t201 * t171 + t204 * t173 + (-m(5) * t247 - t251 + t300) * t283;
t296 = t321 * t427 + t424 * t434;
t306 = (-t372 + t373) * t425;
t308 = (-t201 * t82 + t204 * t81) * t423;
t24 = t283 * t306 + t296 + t308;
t326 = t24 * qJD(1) + t22 * qJD(2);
t325 = -t118 * t311 - t119 * t312;
t323 = t153 * t280 + t154 * t278;
t320 = qJD(2) * t125 + qJD(3) * t157;
t317 = mrSges(7,1) * t416 + t119 * t421;
t315 = m(7) * (-t202 * t82 + t203 * t81);
t310 = m(7) * t325;
t302 = t337 * t424 - t336 / 0.2e1;
t72 = t302 - t346;
t309 = qJD(2) * t344 + t72 * qJD(3);
t223 = Ifges(7,4) * t312;
t159 = -Ifges(7,2) * t311 - t223;
t160 = -Ifges(7,2) * t312 + t395;
t161 = -Ifges(7,1) * t312 - t395;
t162 = Ifges(7,1) * t311 - t223;
t25 = t263 * t157 - (t162 / 0.2e1 + t159 / 0.2e1) * t312 - (-t161 / 0.2e1 + t160 / 0.2e1) * t311;
t290 = -(t117 / 0.4e1 + t128 / 0.4e1) * t312 - (-t129 / 0.4e1 + t115 / 0.4e1) * t311 - (-t155 * mrSges(7,3) / 0.2e1 + t159 / 0.4e1 + t162 / 0.4e1) * t202 + (t156 * t419 + t160 / 0.4e1 - t161 / 0.4e1) * t203 + t155 * t415 + t156 * t414 + t214 * t157 / 0.2e1 + t263 * t125 / 0.2e1 + t283 * t351 / 0.4e1;
t298 = Ifges(7,3) * t438 + t62 * t422 + t63 * mrSges(7,2) / 0.2e1 - t319;
t5 = t290 + t298;
t9 = t341 - t317;
t307 = t9 * qJD(1) + t5 * qJD(2) + t25 * qJD(3);
t292 = (t367 / 0.2e1 - t366 / 0.2e1) * mrSges(7,3) + t324 * t425 + (t155 * t203 - t156 * t202 - t311 * t60 - t312 * t61) * t423 + t173 * t406 + t171 * t405;
t301 = t255 * t425 + t213 * t424 + t394 / 0.2e1 - t391 / 0.2e1;
t13 = (-t236 / 0.2e1 + t385 / 0.2e1) * t280 + (-t238 / 0.2e1 - t383 / 0.2e1) * t278 + t292 + t301;
t297 = t322 * t425 + (-t311 * t81 - t312 * t82) * t424;
t33 = t332 + t297;
t48 = t337 * mrSges(7,3) + t436 + m(7) * (-t155 * t311 - t156 * t312) - t281 * t336;
t305 = -qJD(1) * t33 + qJD(2) * t13 + qJD(3) * t48;
t293 = (t368 / 0.2e1 - t365 / 0.2e1) * mrSges(7,3) + t256 * t425 + (t155 * t204 - t156 * t201 + t214) * t423 + t393 / 0.2e1 - t392 / 0.2e1;
t294 = t170 * t405 + t172 * t406 + t323 * t427 + t424 * t432;
t14 = (t382 / 0.2e1 - t235 / 0.2e1) * t280 + (-t380 / 0.2e1 - t237 / 0.2e1) * t278 + t293 + t294;
t42 = t310 / 0.2e1 + (t349 * t426 - t347) * t430;
t93 = m(7) * t263 + (m(6) + m(5)) * qJ(4) + t431;
t304 = qJD(1) * t42 + qJD(2) * t14 + qJD(3) * t93;
t26 = m(7) * (-t202 * t61 + t203 * t60) - t202 * t171 + t203 * t173 + t300 * t285;
t35 = -t315 / 0.2e1 + (m(7) * t339 + (t339 - t373 / 0.2e1 + t372 / 0.2e1) * m(6)) * t285;
t303 = -t35 * qJD(1) + t26 * qJD(2) + qJD(4) * t344;
t77 = qJD(5) * t344;
t71 = t302 + t346;
t40 = -t310 / 0.2e1 + t335 * t425 + (t428 - t347) * t430;
t36 = t315 / 0.2e1 + t285 * t306 + t346 * t343;
t34 = t332 - t297;
t23 = t308 + (m(5) * t358 + t306) * t283 - t296;
t18 = t299 - t316;
t12 = mrSges(6,2) * t361 / 0.2e1 - mrSges(6,1) * t356 / 0.2e1 + t292 - t301 - t313;
t11 = t237 * t403 + t235 * t402 + (-t390 / 0.2e1 + t388 / 0.2e1 + m(5) * pkin(8) + mrSges(5,1)) * t285 + t293 - t294;
t10 = t341 + t317;
t8 = t295 + t318;
t4 = t290 - t298;
t1 = t338 * t358 + t289 - t291 - t321 * mrSges(6,3) / 0.2e1 - t434 * t419 + t254 * t340 + (t285 * t431 + t384) * t339;
t15 = [t20 * qJD(2) + t19 * qJD(3), t1 * qJD(3) + t23 * qJD(4) + t36 * qJD(5) + t8 * qJD(6) + t369 + (t106 * t173 + t107 * t171 + t175 * t236 + t176 * t238 + ((-mrSges(3,1) + t435) * t284 + (-mrSges(3,2) + (t127 + t218) * t285 + (mrSges(5,1) + mrSges(4,3)) * t440) * t286) * t279 + 0.2e1 * (t106 * t60 + t107 * t61 + t214 * t343) * t423 + 0.2e1 * (t247 * t360 + t350) * t428 + 0.2e1 * (t151 * t175 + t152 * t176 + t256 * t343) * t425 + m(4) * (-pkin(2) * t360 + t350)) * qJD(2), t370 + t1 * qJD(2) + t40 * qJD(4) + t34 * qJD(5) + t10 * qJD(6) + ((t118 * t155 + t119 * t156 - t215 * t263) * t423 + (-pkin(3) * t216 - t376) * t428 + (t281 * t335 - t376) * t425) * t429 + (t325 * mrSges(7,3) + (-mrSges(4,1) + mrSges(5,2) - t436) * t216 + (mrSges(4,2) - t431) * t215) * qJD(3), qJD(2) * t23 + qJD(3) * t40, qJD(2) * t36 + qJD(3) * t34, t8 * qJD(2) + t10 * qJD(3) + (-mrSges(7,1) * t82 - mrSges(7,2) * t81) * qJD(6); qJD(3) * t2 + qJD(4) * t24 - qJD(5) * t35 + qJD(6) * t7 - t369, qJD(3) * t3 + qJD(4) * t22 + qJD(5) * t26 + qJD(6) * t6, t11 * qJD(4) + t12 * qJD(5) + t4 * qJD(6) + ((-qJ(4) * t255 + t281 * t323) * t425 + (t155 * t62 + t156 * t63 + t213 * t263) * t423) * t429 + t328 + (-qJ(4) * t217 + t114 * t405 + t116 * t407 + t263 * t126 + t155 * t172 + t156 * t170 + t213 * t158 + t160 * t412 + t162 * t409 - t255 * t248 + (-t153 * mrSges(6,3) + t281 * t235 + t200 / 0.2e1) * t280 + (-t154 * mrSges(6,3) + t281 * t237 + t413) * t278 + (-pkin(3) * mrSges(5,1) + Ifges(6,5) * t402 + Ifges(7,5) * t407 + Ifges(6,6) * t404 + Ifges(7,6) * t405 - Ifges(5,4) + Ifges(4,5)) * t285 + (-qJ(4) * mrSges(5,1) + (Ifges(6,1) * t280 - t398) * t403 + (-Ifges(6,2) * t278 + t397) * t402 - Ifges(4,6) + Ifges(5,5)) * t283 + (m(5) * t329 + t435) * pkin(8) - t432 * mrSges(7,3)) * qJD(3), t11 * qJD(3) + m(7) * (t365 - t368) * qJD(4) + t77 + t18 * qJD(6) + t326, t12 * qJD(3) + t303, t4 * qJD(3) + t18 * qJD(4) + (-mrSges(7,1) * t61 - mrSges(7,2) * t60 + t352) * qJD(6) + t327; -qJD(2) * t2 + qJD(4) * t42 - qJD(5) * t33 + qJD(6) * t9 - t370, qJD(4) * t14 + qJD(5) * t13 + qJD(6) * t5 - t328, qJD(4) * t93 + qJD(5) * t48 + qJD(6) * t25, qJD(5) * t71 + t304, qJD(4) * t71 + t305 (-mrSges(7,1) * t156 - mrSges(7,2) * t155 + t351) * qJD(6) + t307; -qJD(2) * t24 - qJD(3) * t42, -qJD(3) * t14 + qJD(6) * t17 - t326 + t77, qJD(5) * t72 - t304, 0, t309, -t158 * qJD(6) + t371; qJD(2) * t35 + qJD(3) * t33, -t13 * qJD(3) + t125 * qJD(6) - t303, -qJD(4) * t72 + qJD(6) * t157 - t305, -t309, 0, t320; -t7 * qJD(2) - t9 * qJD(3), -qJD(3) * t5 - qJD(4) * t17 - qJD(5) * t125 - t327, -qJD(5) * t157 - t307, -t371, -t320, 0;];
Cq  = t15;
