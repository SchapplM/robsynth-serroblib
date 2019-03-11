% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:48
% EndTime: 2019-03-08 19:46:57
% DurationCPUTime: 5.22s
% Computational Cost: add. (8809->502), mult. (20744->731), div. (0->0), fcn. (21611->10), ass. (0->262)
t266 = sin(pkin(11));
t262 = t266 ^ 2;
t268 = cos(pkin(11));
t263 = t268 ^ 2;
t335 = t262 + t263;
t320 = t335 * mrSges(6,3);
t425 = (mrSges(5,2) - t320) * qJD(4);
t267 = sin(pkin(6));
t273 = cos(qJ(2));
t270 = sin(qJ(4));
t271 = sin(qJ(2));
t340 = t270 * t271;
t178 = (-t266 * t340 + t268 * t273) * t267;
t179 = (t266 * t273 + t268 * t340) * t267;
t387 = sin(qJ(6));
t388 = cos(qJ(6));
t104 = t178 * t388 - t179 * t387;
t105 = t178 * t387 + t179 * t388;
t385 = pkin(9) + qJ(5);
t240 = t385 * t266;
t242 = t385 * t268;
t167 = -t240 * t388 - t242 * t387;
t168 = -t240 * t387 + t242 * t388;
t258 = -pkin(5) * t268 - pkin(4);
t354 = t179 * t268;
t355 = t178 * t266;
t305 = t354 - t355;
t272 = cos(qJ(4));
t345 = t267 * t271;
t331 = t272 * t345;
t414 = -m(7) / 0.2e1;
t416 = -m(6) / 0.2e1;
t423 = -(t354 / 0.2e1 - t355 / 0.2e1) * mrSges(6,3) + (pkin(4) * t331 + qJ(5) * t305) * t416 + (t104 * t167 + t105 * t168 - t258 * t331) * t414;
t319 = t335 * qJ(5);
t294 = t266 * t388 + t268 * t387;
t205 = t294 * t272;
t377 = t205 * mrSges(7,1);
t421 = t387 * t266 - t388 * t268;
t203 = t421 * t272;
t379 = t203 * mrSges(7,2);
t134 = t377 - t379;
t370 = t268 * mrSges(6,2);
t373 = t266 * mrSges(6,1);
t315 = t370 + t373;
t217 = t315 * t272;
t422 = t134 + t217;
t243 = pkin(4) * t272 + qJ(5) * t270;
t229 = t268 * t243;
t274 = -pkin(2) - pkin(8);
t338 = t272 * t274;
t176 = -t266 * t338 + t229;
t177 = t266 * t243 + t268 * t338;
t306 = -t176 * t266 + t177 * t268;
t386 = pkin(4) * t270;
t239 = -qJ(5) * t272 + qJ(3) + t386;
t223 = t268 * t239;
t339 = t270 * t274;
t173 = -t266 * t339 + t223;
t174 = t266 * t239 + t268 * t339;
t346 = t266 * t272;
t231 = -t270 * mrSges(6,2) - mrSges(6,3) * t346;
t341 = t268 * t272;
t233 = t270 * mrSges(6,1) - mrSges(6,3) * t341;
t419 = -m(6) * (t173 * t268 + t174 * t266) - t266 * t231 - t268 * t233;
t418 = 0.2e1 * qJD(4);
t417 = m(5) / 0.2e1;
t415 = m(6) / 0.2e1;
t413 = -m(7) / 0.4e1;
t412 = m(7) / 0.2e1;
t411 = mrSges(7,1) / 0.2e1;
t410 = -mrSges(7,2) / 0.2e1;
t409 = -mrSges(7,3) / 0.2e1;
t408 = mrSges(7,3) / 0.2e1;
t204 = t294 * t270;
t206 = t421 * t270;
t407 = m(7) * (-t203 * t204 + t205 * t206);
t170 = -mrSges(7,2) * t270 - t205 * mrSges(7,3);
t406 = -t170 / 0.2e1;
t405 = t170 / 0.2e1;
t172 = mrSges(7,1) * t270 + t203 * mrSges(7,3);
t404 = -t172 / 0.2e1;
t403 = t172 / 0.2e1;
t402 = -t203 / 0.2e1;
t401 = t204 / 0.2e1;
t400 = -t205 / 0.2e1;
t399 = t206 / 0.2e1;
t397 = -t421 / 0.2e1;
t396 = -t294 / 0.2e1;
t347 = t266 * t270;
t230 = -mrSges(6,2) * t272 + mrSges(6,3) * t347;
t395 = t230 / 0.2e1;
t394 = t266 / 0.2e1;
t393 = -t268 / 0.2e1;
t392 = t268 / 0.2e1;
t391 = t270 / 0.2e1;
t390 = -t272 / 0.2e1;
t389 = t272 / 0.2e1;
t269 = cos(pkin(6));
t344 = t267 * t273;
t214 = t269 * t270 + t272 * t344;
t127 = t294 * t214;
t128 = t421 * t214;
t215 = t269 * t272 - t270 * t344;
t163 = -t215 * t266 + t268 * t345;
t164 = t215 * t268 + t266 * t345;
t309 = t163 * t266 - t164 * t268;
t334 = -0.1e1 + t335;
t350 = t214 * t270;
t87 = t163 * t387 + t164 * t388;
t365 = t87 * t203;
t86 = t163 * t388 - t164 * t387;
t366 = t86 * t205;
t20 = ((-t215 - t309) * t272 - t334 * t350) * t415 + (-t127 * t204 - t128 * t206 - t215 * t272 + t350 - t365 - t366) * t412;
t384 = t20 * qJD(4);
t383 = Ifges(6,4) * t266;
t382 = Ifges(6,4) * t268;
t381 = Ifges(7,4) * t203;
t380 = Ifges(7,4) * t294;
t378 = t204 * mrSges(7,1);
t376 = t206 * mrSges(7,2);
t375 = t421 * mrSges(7,1);
t374 = t294 * mrSges(7,2);
t372 = t266 * Ifges(6,2);
t371 = t266 * Ifges(6,6);
t369 = t268 * Ifges(6,5);
t368 = t270 * mrSges(5,2);
t367 = t272 * mrSges(5,1);
t241 = -mrSges(6,1) * t268 + mrSges(6,2) * t266;
t364 = t241 - mrSges(5,1);
t363 = mrSges(7,3) * qJD(4);
t362 = t104 * t294;
t361 = t105 * t421;
t132 = mrSges(7,1) * t203 + t205 * mrSges(7,2);
t352 = t206 * t203;
t353 = t204 * t205;
t282 = (-t352 / 0.2e1 - t353 / 0.2e1) * mrSges(7,3) + t204 * t406 + t172 * t399 + t132 * t389;
t297 = t375 / 0.2e1 + t374 / 0.2e1;
t14 = t282 + t297;
t360 = t14 * qJD(2);
t359 = t163 * t268;
t358 = t164 * t266;
t156 = mrSges(7,1) * t294 - mrSges(7,2) * t421;
t351 = t214 * t156;
t318 = t214 * t331;
t24 = m(6) * (t163 * t178 + t164 * t179 - t318) + m(7) * (t104 * t86 + t105 * t87 - t318) + m(5) * (-t214 * t272 + t215 * t270 + t344) * t345;
t349 = t24 * qJD(1);
t342 = t268 * t270;
t232 = mrSges(6,1) * t272 + mrSges(6,3) * t342;
t348 = t266 * t232;
t343 = t268 * t230;
t337 = -Ifges(7,5) * t205 + Ifges(7,6) * t203;
t336 = -Ifges(7,5) * t421 - Ifges(7,6) * t294;
t333 = t407 / 0.2e1;
t332 = m(6) * t391;
t330 = t274 * t345;
t327 = -t347 / 0.2e1;
t326 = -t345 / 0.2e1;
t324 = t156 * t390;
t323 = t134 / 0.2e1 + t217 / 0.2e1;
t322 = -t266 * t274 + pkin(5);
t321 = pkin(5) * t266 - t274;
t157 = t374 + t375;
t316 = qJD(4) * (t157 + t364);
t244 = t367 - t368;
t169 = -mrSges(7,2) * t272 + mrSges(7,3) * t204;
t171 = mrSges(7,1) * t272 - mrSges(7,3) * t206;
t224 = t321 * t270;
t225 = t321 * t272;
t133 = t376 - t378;
t216 = t315 * t270;
t295 = t231 * t393 + t233 * t394;
t287 = t133 / 0.2e1 - t216 / 0.2e1 + t295;
t307 = -t173 * t266 + t174 * t268;
t144 = -pkin(9) * t341 + t270 * t322 + t223;
t155 = -pkin(9) * t346 + t174;
t73 = t144 * t388 - t155 * t387;
t74 = t144 * t387 + t155 * t388;
t149 = pkin(9) * t342 + t272 * t322 + t229;
t162 = pkin(9) * t347 + t177;
t79 = t149 * t388 - t162 * t387;
t80 = t149 * t387 + t162 * t388;
t275 = t323 * t215 + t287 * t214 + (-t215 * t338 + t176 * t163 + t177 * t164 + (-t307 + t339) * t214) * t415 + (t127 * t73 + t128 * t74 - t214 * t224 + t215 * t225 + t79 * t86 + t80 * t87) * t412 + t127 * t403 + t128 * t405 + t163 * t232 / 0.2e1 + t164 * t395 + t86 * t171 / 0.2e1 + t87 * t169 / 0.2e1;
t2 = t275 + (t361 / 0.2e1 + t362 / 0.2e1) * mrSges(7,3) + (t244 / 0.2e1 + t368 / 0.2e1 + (-mrSges(5,1) / 0.2e1 + t241 / 0.2e1 + t157 / 0.2e1) * t272) * t345 + t423;
t123 = Ifges(7,4) * t206 + Ifges(7,2) * t204 + Ifges(7,6) * t272;
t124 = -Ifges(7,2) * t205 + t270 * Ifges(7,6) - t381;
t125 = Ifges(7,1) * t206 + Ifges(7,4) * t204 + Ifges(7,5) * t272;
t194 = Ifges(7,4) * t205;
t126 = -Ifges(7,1) * t203 + t270 * Ifges(7,5) - t194;
t201 = Ifges(6,6) * t272 + (t372 - t382) * t270;
t202 = Ifges(6,5) * t272 + (-t268 * Ifges(6,1) + t383) * t270;
t301 = Ifges(7,5) * t399 + Ifges(7,6) * t401;
t3 = qJ(3) * t244 + t173 * t232 + t176 * t233 + t174 * t230 + t177 * t231 - t224 * t134 + t225 * t133 + t126 * t399 + t123 * t400 + t125 * t402 + t124 * t401 + t73 * t171 + t79 * t172 + t74 * t169 + t80 * t170 + m(6) * (t173 * t176 + t174 * t177) + m(7) * (-t224 * t225 + t73 * t79 + t74 * t80) + (t274 * t216 + t202 * t392 - t266 * t201 / 0.2e1 + Ifges(7,5) * t402 + Ifges(7,6) * t400 + (-Ifges(5,4) + t369 / 0.2e1 - t371 / 0.2e1) * t272) * t272 + (t274 * t217 + (Ifges(5,4) - t369 + t371) * t270 + (-m(6) * t274 ^ 2 - t263 * Ifges(6,1) / 0.2e1 + Ifges(6,3) + Ifges(7,3) - Ifges(5,1) + Ifges(5,2) + (t382 - t372 / 0.2e1) * t266) * t272 + t301) * t270;
t314 = t2 * qJD(1) + t3 * qJD(2);
t288 = t214 * t132 / 0.2e1 + t86 * t406 + t87 * t403;
t300 = t104 * t411 + t105 * t410;
t6 = (-t365 / 0.2e1 - t366 / 0.2e1) * mrSges(7,3) + t288 + t300;
t135 = Ifges(7,2) * t203 - t194;
t136 = -Ifges(7,1) * t205 + t381;
t8 = -t74 * t172 + t73 * t170 + t337 * t391 - t225 * t132 - (-t73 * mrSges(7,3) + t126 / 0.2e1 + t135 / 0.2e1) * t205 + (t74 * mrSges(7,3) - t136 / 0.2e1 + t124 / 0.2e1) * t203;
t313 = -t6 * qJD(1) + t8 * qJD(2);
t265 = t272 ^ 2;
t248 = t265 * t345;
t264 = t270 ^ 2;
t279 = (t270 * t305 + t248) * t415 + (-t104 * t204 - t105 * t206 + t248) * t412 + (t264 * t345 + t248) * t417;
t310 = t358 + t359;
t281 = t310 * t416 + (t294 * t87 - t421 * t86) * t414 + m(5) * t326;
t22 = t279 + t281;
t302 = t270 * mrSges(5,1) + t272 * mrSges(5,2) + mrSges(4,3);
t31 = t294 * t170 - t421 * t172 + (m(5) + m(4)) * qJ(3) + m(7) * (t294 * t74 - t421 * t73) + t302 - t419;
t312 = -t22 * qJD(1) + t31 * qJD(2);
t154 = t214 * t215;
t23 = m(6) * (t214 * t309 + t154) + m(7) * (t127 * t86 + t128 * t87 + t154);
t311 = t23 * qJD(1) + t20 * qJD(3);
t303 = qJD(2) * t132 - qJD(4) * t156;
t299 = t127 * t411 + t128 * t410;
t298 = -t379 / 0.2e1 + t377 / 0.2e1;
t296 = m(7) * (t203 * t86 - t205 * t87);
t293 = (-t167 * t421 + t168 * t294) * t412;
t292 = (t204 * t294 + t206 * t421) * t412;
t53 = m(7) * (t352 + t353) + 0.4e1 * (m(6) * t334 / 0.4e1 + t413) * t272 * t270;
t276 = (t307 * t272 + (t306 - 0.2e1 * t338) * t270) * t416 + (-t74 * t203 - t204 * t79 - t73 * t205 - t206 * t80 + t224 * t272 + t225 * t270) * t414 + t203 * t405 + t171 * t401 - t205 * t404 + t169 * t399;
t9 = t287 * t272 + (-t343 / 0.2e1 + t348 / 0.2e1 - t323) * t270 + t293 + t276;
t291 = t20 * qJD(1) - t9 * qJD(2) + t53 * qJD(3);
t290 = -t224 * t414 + t378 / 0.2e1 - t376 / 0.2e1;
t29 = t203 * t172 - t205 * t170 + m(7) * (t203 * t73 - t205 * t74) + t419 * t272;
t39 = -t296 / 0.2e1 + (m(7) * t326 + (t326 + t359 / 0.2e1 + t358 / 0.2e1) * m(6)) * t272;
t289 = -t39 * qJD(1) + t29 * qJD(2) + qJD(3) * t333;
t11 = -t351 / 0.2e1 + t299;
t221 = Ifges(7,4) * t421;
t158 = -Ifges(7,2) * t294 - t221;
t159 = -Ifges(7,2) * t421 + t380;
t160 = -Ifges(7,1) * t421 - t380;
t161 = Ifges(7,1) * t294 - t221;
t25 = t258 * t156 - (-t160 / 0.2e1 + t159 / 0.2e1) * t294 - (t161 / 0.2e1 + t158 / 0.2e1) * t421;
t35 = t324 + t298;
t277 = -(t124 / 0.4e1 - t136 / 0.4e1) * t294 - (t135 / 0.4e1 + t126 / 0.4e1) * t421 + (-t160 / 0.4e1 + t159 / 0.4e1 + t168 * t408) * t203 - (t161 / 0.4e1 + t158 / 0.4e1 + t167 * t409) * t205 + t167 * t405 + t168 * t404 + t225 * t156 / 0.2e1 - t258 * t132 / 0.2e1 + t270 * t336 / 0.4e1;
t284 = Ifges(7,3) * t390 - t79 * mrSges(7,1) / 0.2e1 + t80 * mrSges(7,2) / 0.2e1 - t301;
t5 = t277 + t284;
t286 = -t11 * qJD(1) + t5 * qJD(2) + t35 * qJD(3) + t25 * qJD(4);
t278 = (t203 * t396 - t205 * t397) * mrSges(7,3) + t307 * t415 + (t167 * t203 - t168 * t205 - t294 * t73 - t421 * t74) * t412 + t170 * t397 + t172 * t396 - t295;
t16 = (t274 * t416 + t370 / 0.2e1 + t373 / 0.2e1) * t270 + t278 + t290;
t283 = t309 * t416 + (-t294 * t86 - t421 * t87) * t412;
t38 = 0.2e1 * (-m(6) / 0.4e1 + t413) * t215 + t283;
t50 = (t294 ^ 2 + t421 ^ 2) * mrSges(7,3) + t320 + m(7) * (-t167 * t294 - t168 * t421) + m(6) * t319;
t55 = t292 + (t414 + (t262 / 0.2e1 + t263 / 0.2e1 - 0.1e1 / 0.2e1) * m(6)) * t270;
t285 = qJD(1) * t38 + qJD(2) * t16 + qJD(3) * t55 + qJD(4) * t50;
t255 = qJ(3) * t344;
t236 = t265 * t330;
t67 = qJD(5) * t333;
t54 = m(7) * t391 + t332 * t335 + t292 + t332;
t40 = m(6) * t310 * t390 + t296 / 0.2e1 + (t416 + t414) * t331;
t37 = t283 + (m(6) + m(7)) * t215 / 0.2e1;
t36 = t324 - t298;
t21 = m(4) * t345 + t279 - t281;
t15 = t274 * t332 - mrSges(6,2) * t342 / 0.2e1 + mrSges(6,1) * t327 + t278 - t290;
t13 = t282 - t297;
t12 = t351 / 0.2e1 + t299;
t10 = t231 * t341 / 0.2e1 + t342 * t395 - t233 * t346 / 0.2e1 + t232 * t327 + t293 - t276 + t422 * t391 + (t133 - t216) * t390;
t7 = t365 * t408 - t366 * t409 - t288 + t300;
t4 = t277 - t284;
t1 = t275 + (-t361 - t362) * t408 + (t244 + t367) * t345 / 0.2e1 + (t368 + (t241 + t157) * t272) * t326 - t423;
t17 = [t24 * qJD(2) + t23 * qJD(4), t21 * qJD(3) + t1 * qJD(4) + t40 * qJD(5) + t7 * qJD(6) + t349 + (t104 * t172 + t105 * t170 + t178 * t233 + t179 * t231 + ((-mrSges(3,2) + t302) * t273 + (-mrSges(3,1) + mrSges(4,2) - t422 * t272 + (-t264 - t265) * mrSges(5,3)) * t271) * t267 + 0.2e1 * (t173 * t178 + t174 * t179 + t236) * t415 + 0.2e1 * (t104 * t73 + t105 * t74 - t225 * t331) * t412 + 0.2e1 * (t264 * t330 + t236 + t255) * t417 + m(4) * (-pkin(2) * t345 + t255)) * qJD(2), qJD(2) * t21 + t384, t1 * qJD(2) + t37 * qJD(5) + t12 * qJD(6) + (-t127 * t294 - t128 * t421) * t363 + t215 * t316 + t214 * t425 + ((-pkin(4) * t215 - t214 * t319) * t415 + (t127 * t167 + t128 * t168 + t215 * t258) * t412) * t418 + t311, qJD(2) * t40 + qJD(4) * t37, t7 * qJD(2) + t12 * qJD(4) + (-mrSges(7,1) * t87 - mrSges(7,2) * t86) * qJD(6); -qJD(3) * t22 + qJD(4) * t2 - qJD(5) * t39 - qJD(6) * t6 - t349, qJD(3) * t31 + qJD(4) * t3 + qJD(5) * t29 + qJD(6) * t8, m(7) * (t204 * t421 - t206 * t294) * qJD(3) + t10 * qJD(4) + t67 + t13 * qJD(6) + t312, t10 * qJD(3) + t15 * qJD(5) + t4 * qJD(6) + t314 + (-mrSges(5,2) * t338 + m(7) * (t167 * t79 + t168 * t80 - t224 * t258) - Ifges(5,6) * t272 + t201 * t392 + t202 * t394 + t258 * t133 + t294 * t125 / 0.2e1 + t123 * t397 - t224 * t157 + pkin(4) * t216 + t161 * t399 + t159 * t401 + t167 * t171 + t168 * t169 + (m(6) * t306 + t343 - t348) * qJ(5) + ((Ifges(6,2) * t268 + t383) * t394 + (Ifges(6,1) * t266 + t382) * t393 - Ifges(5,5) + (-m(6) * pkin(4) + t364) * t274) * t270 + (Ifges(6,5) * t266 + Ifges(7,5) * t294 + Ifges(6,6) * t268 - Ifges(7,6) * t421) * t389 + (-t294 * t79 - t421 * t80) * mrSges(7,3) + t306 * mrSges(6,3)) * qJD(4), t15 * qJD(4) + t289, t13 * qJD(3) + t4 * qJD(4) + (-mrSges(7,1) * t74 - mrSges(7,2) * t73 + t337) * qJD(6) + t313; qJD(2) * t22 + t384, -qJD(4) * t9 + qJD(6) * t14 - t312 + t67, t53 * qJD(4), t54 * qJD(5) + t36 * qJD(6) + (t203 * t421 + t205 * t294) * t363 - t272 * t425 + t270 * t316 + ((t272 * t319 - t386) * t415 + (-t167 * t205 - t168 * t203 + t258 * t270) * t412) * t418 + t291, qJD(2) * t333 + t54 * qJD(4), t360 + t36 * qJD(4) + (mrSges(7,1) * t206 + mrSges(7,2) * t204) * qJD(6); -qJD(2) * t2 + qJD(5) * t38 - qJD(6) * t11 - t311, qJD(3) * t9 + qJD(5) * t16 + qJD(6) * t5 - t314, qJD(5) * t55 + qJD(6) * t35 - t291, qJD(5) * t50 + qJD(6) * t25, t285 (-mrSges(7,1) * t168 - mrSges(7,2) * t167 + t336) * qJD(6) + t286; qJD(2) * t39 - qJD(4) * t38, -t16 * qJD(4) - t132 * qJD(6) - t289, -qJD(2) * t407 / 0.2e1 - t55 * qJD(4), qJD(6) * t156 - t285, 0, -t303; t6 * qJD(2) + t11 * qJD(4), -qJD(3) * t14 - qJD(4) * t5 + qJD(5) * t132 - t313, -t35 * qJD(4) - t360, -qJD(5) * t156 - t286, t303, 0;];
Cq  = t17;
