% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR13_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:28
% EndTime: 2024-09-27 17:30:35
% DurationCPUTime: 5.88s
% Computational Cost: add. (19514->453), mult. (45026->589), div. (0->0), fcn. (42333->10), ass. (0->286)
t441 = mrSges(6,1) / 0.2e1;
t439 = -mrSges(6,2) / 0.2e1;
t253 = sin(pkin(5));
t254 = sin(qJ(5));
t257 = cos(qJ(5));
t423 = sin(qJ(4));
t424 = cos(qJ(4));
t327 = t424 * t254 + t423 * t257;
t463 = t253 * t327;
t314 = mrSges(6,3) * t463;
t328 = t423 * t254 - t424 * t257;
t169 = (t328 * mrSges(6,1) + t327 * mrSges(6,2)) * t253;
t462 = t169 + (-mrSges(5,1) * t424 + mrSges(5,2) * t423) * t253;
t280 = t327 * mrSges(6,1) - t328 * mrSges(6,2);
t259 = cos(qJ(2));
t422 = pkin(1) * t259;
t249 = pkin(2) + t422;
t258 = cos(qJ(3));
t255 = sin(qJ(3));
t256 = sin(qJ(2));
t382 = t255 * t256;
t231 = -pkin(1) * t382 + t258 * t249;
t225 = pkin(3) + t231;
t252 = t253 ^ 2;
t373 = t424 * pkin(4);
t337 = t252 * (-t373 - t225);
t117 = t280 * t337;
t381 = t256 * t258;
t232 = pkin(1) * t381 + t249 * t255;
t250 = t253 * pkin(9);
t216 = t232 + t250;
t403 = cos(pkin(5));
t359 = t403 * t424;
t357 = -t423 * t216 + t225 * t359;
t369 = t253 * t423;
t364 = pkin(10) * t369;
t145 = -t364 + t357;
t251 = t403 * pkin(4);
t131 = t251 + t145;
t358 = t403 * t423;
t310 = t424 * t216 + t225 * t358;
t370 = t253 * t424;
t365 = pkin(10) * t370;
t146 = t365 + t310;
t396 = t146 * t257;
t68 = t131 * t254 + t396;
t62 = t68 * t314;
t472 = -t117 / 0.2e1 + t62 / 0.2e1;
t248 = t258 * pkin(2) + pkin(3);
t230 = (-t373 - t248) * t253;
t277 = t253 * t280;
t136 = t230 * t277;
t421 = pkin(2) * t255;
t243 = t250 + t421;
t200 = -t423 * t243 + t248 * t359;
t182 = -t364 + t200;
t177 = t251 + t182;
t201 = t424 * t243 + t248 * t358;
t183 = t365 + t201;
t390 = t183 * t257;
t111 = t177 * t254 + t390;
t93 = t111 * t314;
t471 = -t136 / 0.2e1 + t93 / 0.2e1;
t245 = pkin(3) * t359;
t465 = pkin(9) + pkin(10);
t213 = -t465 * t369 + t245;
t196 = t251 + t213;
t340 = pkin(3) * t358;
t214 = t465 * t370 + t340;
t387 = t214 * t257;
t144 = t196 * t254 + t387;
t112 = t144 * t314;
t240 = (-t373 - pkin(3)) * t253;
t150 = t240 * t277;
t470 = t112 / 0.2e1 - t150 / 0.2e1;
t199 = t403 * mrSges(6,1) - t314;
t401 = t111 * t199;
t391 = t183 * t254;
t110 = t177 * t257 - t391;
t318 = t253 * t328;
t313 = mrSges(6,3) * t318;
t198 = -t403 * mrSges(6,2) - t313;
t402 = t110 * t198;
t464 = t402 / 0.2e1 - t401 / 0.2e1;
t404 = t68 * t199;
t397 = t146 * t254;
t67 = t131 * t257 - t397;
t405 = t67 * t198;
t469 = -t404 / 0.2e1 + t405 / 0.2e1;
t362 = mrSges(5,3) * t370;
t322 = -t403 * mrSges(5,2) + t362;
t174 = t200 * t322;
t361 = mrSges(5,3) * t369;
t321 = t403 * mrSges(5,1) - t361;
t175 = t201 * t321;
t343 = -t182 * t254 - t390;
t86 = t343 * t199;
t344 = t182 * t257 - t391;
t87 = t344 * t198;
t468 = t174 / 0.2e1 - t175 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1;
t341 = -t213 * t254 - t387;
t106 = t341 * t199;
t388 = t214 * t254;
t342 = t213 * t257 - t388;
t107 = t342 * t198;
t235 = -pkin(9) * t369 + t245;
t192 = t235 * t322;
t236 = pkin(9) * t370 + t340;
t193 = t236 * t321;
t467 = -t106 / 0.2e1 - t107 / 0.2e1 - t192 / 0.2e1 + t193 / 0.2e1;
t366 = pkin(4) * t369;
t371 = t423 * Ifges(5,4);
t372 = t424 * Ifges(5,4);
t451 = t169 * t366 - (Ifges(5,6) * t403 + (t424 * Ifges(5,2) + t371) * t253) * t369 / 0.2e1 + (Ifges(5,5) * t403 + (t423 * Ifges(5,1) + t372) * t253) * t370 / 0.2e1 + t403 * (t424 * Ifges(5,5) - t423 * Ifges(5,6)) * t253 / 0.2e1 + (t423 * (t424 * Ifges(5,1) - t371) + t424 * (-t423 * Ifges(5,2) + t372)) * t252 / 0.2e1;
t127 = t357 * t322;
t128 = t310 * t321;
t287 = t310 * t361;
t315 = t357 * t362;
t75 = -t145 * t254 - t396;
t58 = t75 * t199;
t76 = t145 * t257 - t397;
t59 = t76 * t198;
t466 = t127 / 0.2e1 - t128 / 0.2e1 - t287 / 0.2e1 - t315 / 0.2e1 + t58 / 0.2e1 + t59 / 0.2e1;
t229 = (-t255 * t358 + t424 * t258) * pkin(2);
t308 = t229 * t322;
t228 = (-t255 * t359 - t423 * t258) * pkin(2);
t309 = t228 * t321;
t168 = t228 * t254 + t229 * t257;
t394 = t168 * t198;
t167 = t228 * t257 - t229 * t254;
t395 = t167 * t199;
t406 = t258 * mrSges(4,2);
t461 = t308 + t309 + t394 + t395 + (-t406 + (t462 * t253 - mrSges(4,1)) * t255) * pkin(2);
t298 = t313 / 0.2e1;
t459 = t67 * t298 - t472;
t458 = t110 * t298 - t471;
t143 = t196 * t257 - t388;
t297 = -t313 / 0.2e1;
t457 = t143 * t297 + t470;
t456 = t110 * t297 + t471;
t455 = t67 * t297 + t472;
t398 = t144 * t199;
t399 = t143 * t198;
t454 = -t398 / 0.2e1 + t399 / 0.2e1;
t380 = t167 * t441 + t168 * t439;
t237 = (-t255 * t259 - t381) * pkin(1);
t383 = t253 * t237;
t453 = t462 * t383;
t452 = t459 + t469;
t420 = pkin(4) * t254;
t209 = t314 * t420;
t419 = pkin(4) * t257;
t295 = t313 * t419;
t360 = -t199 * t420 / 0.2e1 + t198 * t419 / 0.2e1 - t209 / 0.2e1 + t295 / 0.2e1;
t320 = t252 * (mrSges(5,1) * t423 + mrSges(5,2) * t424);
t312 = -t320 / 0.2e1;
t346 = -t361 / 0.2e1;
t348 = -t362 / 0.2e1;
t450 = t200 * t348 + t201 * t346 + t248 * t312 + t468;
t379 = -Ifges(6,5) * t318 - Ifges(6,6) * t463;
t263 = -t403 * t379 / 0.2e1 + (Ifges(6,1) * t463 - Ifges(6,4) * t318 + Ifges(6,5) * t403) * t318 / 0.2e1 + (Ifges(6,4) * t463 - Ifges(6,2) * t318 + Ifges(6,6) * t403) * t463 / 0.2e1 + (t328 * (-Ifges(6,4) * t328 - Ifges(6,2) * t327) / 0.2e1 - t327 * (-Ifges(6,1) * t328 - Ifges(6,4) * t327) / 0.2e1) * t252;
t449 = t225 * t312 + t459 + t466;
t311 = t320 / 0.2e1;
t345 = t361 / 0.2e1;
t347 = t362 / 0.2e1;
t448 = pkin(3) * t311 + t235 * t347 + t236 * t345 + t457 + t467;
t447 = t200 * t347 + t201 * t345 + t248 * t311 + t456 - t468;
t446 = 2 * qJD(3);
t445 = m(5) / 0.2e1;
t444 = m(6) / 0.2e1;
t443 = m(6) * pkin(4);
t442 = mrSges(5,1) / 0.2e1;
t440 = -mrSges(5,2) / 0.2e1;
t418 = t67 * mrSges(6,2);
t417 = t68 * mrSges(6,1);
t416 = t75 * mrSges(6,1);
t415 = t76 * mrSges(6,2);
t414 = t110 * mrSges(6,2);
t413 = t111 * mrSges(6,1);
t412 = t143 * mrSges(6,2);
t411 = t144 * mrSges(6,1);
t408 = t231 * mrSges(4,2);
t238 = (t258 * t259 - t382) * pkin(1);
t407 = t238 * mrSges(4,2);
t261 = t263 - t451;
t291 = t67 * t313 + t117 - t62;
t319 = t423 * pkin(4) * t337;
t13 = -m(6) * (t67 * t75 + t68 * t76 + t319) - t58 + t225 * t320 + t128 - t127 + t287 + t315 - t59 + t261 - t291;
t400 = t13 * qJD(1);
t17 = t291 - t263 - t404 + t405;
t393 = t17 * qJD(1);
t163 = -t423 * t231 - t232 * t359;
t164 = t424 * t231 - t232 * t358;
t141 = t163 * t321;
t142 = t164 * t322;
t224 = t232 * mrSges(4,1);
t384 = t253 * t232;
t94 = t163 * t257 - t164 * t254;
t73 = t94 * t199;
t95 = t163 * t254 + t164 * t257;
t74 = t95 * t198;
t288 = t462 * t384 + t141 + t142 - t224 - t408 + t73 + t74;
t386 = t232 * t252;
t18 = -m(6) * (t232 * t337 + t67 * t94 + t68 * t95) - m(5) * (t357 * t163 + t310 * t164 - t225 * t386) - t288;
t392 = t18 * qJD(1);
t172 = t237 * t359 - t423 * t238;
t173 = t237 * t358 + t424 * t238;
t153 = t172 * t321;
t154 = t173 * t322;
t233 = t237 * mrSges(4,1);
t98 = t172 * t257 - t173 * t254;
t79 = t98 * t199;
t99 = t172 * t254 + t173 * t257;
t80 = t99 * t198;
t276 = t256 * pkin(1) * mrSges(3,1) + mrSges(3,2) * t422 - t153 - t154 - t233 + t407 + t453 - t79 - t80;
t385 = t237 * t252;
t20 = -m(6) * (-t237 * t337 + t67 * t98 + t68 * t99) - m(5) * (t357 * t172 + t310 * t173 + t225 * t385) - m(4) * (t231 * t237 + t232 * t238) + t276;
t389 = t20 * qJD(1);
t378 = t443 / 0.2e1;
t377 = t252 * t421;
t376 = t253 * t421;
t356 = t230 * t366;
t355 = t240 * t366;
t290 = t110 * t313 + t136 - t93;
t16 = -m(6) * (t110 * t343 + t111 * t344 + t356) - t174 + t175 + t248 * t320 + t201 * t361 + t200 * t362 - t87 - t86 + t261 - t290;
t260 = t225 * t311 + t261 + t455 - t466;
t268 = m(6) * (t110 * t75 + t111 * t76 + t343 * t67 + t344 * t68 + t319 + t356);
t338 = t99 * t439 + t98 * t441;
t278 = t172 * t442 + t173 * t440 + (t254 * t99 + t257 * t98) * t378 + t338;
t2 = -t268 / 0.2e1 + t260 + t278 + t447;
t352 = -t2 * qJD(1) - t16 * qJD(2);
t21 = t290 - t263 - t401 + t402;
t262 = t263 + t455 - t469;
t9 = t262 + t338 + t456 - t464;
t351 = -t9 * qJD(1) + t21 * qJD(2);
t30 = m(6) * (t110 * t167 + t111 * t168 + t230 * t376) + m(5) * (t200 * t228 + t201 * t229 - t248 * t377) + t461;
t264 = t141 / 0.2e1 + t142 / 0.2e1 - t224 / 0.2e1 + t73 / 0.2e1 + t74 / 0.2e1 + (t229 * t310 + t201 * t164 + t228 * t357 + t200 * t163 + (-t225 * t421 - t232 * t248) * t252) * t445 + (t110 * t94 + t111 * t95 + t167 * t67 + t168 * t68 + t230 * t384 + t337 * t421) * t444 + t395 / 0.2e1 + t394 / 0.2e1 + t309 / 0.2e1 + t308 / 0.2e1 - t408 / 0.2e1 - mrSges(4,1) * t421 / 0.2e1 - pkin(2) * t406 / 0.2e1 + t462 * (t384 / 0.2e1 + t376 / 0.2e1);
t266 = -t153 / 0.2e1 - t154 / 0.2e1 - t233 / 0.2e1 - t79 / 0.2e1 - t80 / 0.2e1 - m(5) * (pkin(3) * t385 + t172 * t235 + t173 * t236) / 0.2e1 - m(6) * (t143 * t98 + t144 * t99 - t240 * t383) / 0.2e1 + t407 / 0.2e1 + t453 / 0.2e1;
t6 = t264 + t266;
t350 = t6 * qJD(1) + t30 * qJD(2);
t349 = (-t254 ^ 2 - t257 ^ 2) * t443;
t339 = t95 * t439 + t94 * t441;
t332 = t344 * mrSges(6,2);
t331 = t343 * mrSges(6,1);
t330 = t342 * mrSges(6,2);
t329 = t341 * mrSges(6,1);
t289 = t143 * t313 - t112 + t150;
t19 = -m(6) * (t143 * t341 + t144 * t342 + t355) + pkin(3) * t320 + t236 * t361 + t235 * t362 - t192 + t193 - t107 - t106 + t261 - t289;
t267 = m(6) * (t143 * t75 + t144 * t76 + t341 * t67 + t342 * t68 + t319 + t355);
t279 = t163 * t442 + t164 * t440 + (t254 * t95 + t257 * t94) * t378 + t339;
t4 = -t267 / 0.2e1 + t260 + t279 + t448;
t265 = m(6) * (t341 * t110 + t342 * t111 + t143 * t343 + t144 * t344 + (t230 + t240) * t366);
t283 = t228 * t442 + t229 * t440 + (t167 * t257 + t168 * t254) * t378 + t380;
t8 = t447 + t448 + t283 - t265 / 0.2e1 + t261;
t326 = -t4 * qJD(1) - t8 * qJD(2) - t19 * qJD(3);
t325 = -t263 + t458;
t11 = t262 + t339 - t454 + t457;
t323 = t143 * t298 - t263 - t470;
t305 = t323 + t454;
t286 = t305 + t458 + t464;
t15 = t286 - t380;
t23 = t289 - t263 - t398 + t399;
t324 = -t11 * qJD(1) + t15 * qJD(2) + t23 * qJD(3);
t307 = -t414 / 0.2e1 - t413 / 0.2e1 + t360;
t306 = -t412 / 0.2e1 - t411 / 0.2e1 + t360;
t24 = (-t76 / 0.2e1 + t67 / 0.2e1) * mrSges(6,2) + (t75 / 0.2e1 + t68 / 0.2e1) * mrSges(6,1) - t360;
t242 = (mrSges(6,1) * t254 + mrSges(6,2) * t257) * pkin(4);
t285 = t332 / 0.2e1 - t331 / 0.2e1;
t32 = t285 + t307;
t284 = t330 / 0.2e1 - t329 / 0.2e1;
t39 = t284 + t306;
t304 = qJD(1) * t24 - qJD(2) * t32 - qJD(3) * t39 + qJD(4) * t242;
t275 = pkin(3) * t312 + t235 * t348 + t236 * t346 + t323 + t451 - t467;
t274 = Ifges(5,5) * t370 - Ifges(5,6) * t369 - t209 + t295 + t379;
t239 = t242 * qJD(5);
t36 = -t284 + t306 + t379;
t31 = -t285 + t307 + t379;
t22 = -t417 / 0.2e1 - t418 / 0.2e1 + t416 / 0.2e1 - t415 / 0.2e1 + t360 + t379;
t14 = t286 + t380;
t12 = t305 + t339 + t452;
t10 = t325 + t338 + t452 + t464;
t7 = t275 + t265 / 0.2e1 + t283 + t450 + t458;
t5 = t264 - t266;
t3 = t275 + t267 / 0.2e1 + t279 + t449;
t1 = t450 + t449 + t325 + t278 + t268 / 0.2e1 + t451;
t25 = [-qJD(2) * t20 - qJD(3) * t18 - qJD(4) * t13 + qJD(5) * t17, t5 * qJD(3) + t1 * qJD(4) + t10 * qJD(5) - t389 + (-t276 + 0.2e1 * (t110 * t98 + t111 * t99 - t230 * t383) * t444 + 0.2e1 * (t172 * t200 + t173 * t201 + t248 * t385) * t445 + m(4) * (t237 * t258 + t238 * t255) * pkin(2)) * qJD(2), -t392 + t5 * qJD(2) + t288 * qJD(3) + t3 * qJD(4) + t12 * qJD(5) + ((t143 * t94 + t144 * t95 + t240 * t384) * t444 + (-pkin(3) * t386 + t163 * t235 + t164 * t236) * t445) * t446, -t400 + t1 * qJD(2) + t3 * qJD(3) + (t416 - t415 + (t254 * t76 + t257 * t75) * t443 - t357 * mrSges(5,2) - t310 * mrSges(5,1) + t274) * qJD(4) + t22 * qJD(5), t393 + t10 * qJD(2) + t12 * qJD(3) + t22 * qJD(4) + (t379 - t417 - t418) * qJD(5); qJD(3) * t6 - qJD(4) * t2 - qJD(5) * t9 + t389, qJD(3) * t30 - qJD(4) * t16 + qJD(5) * t21, t7 * qJD(4) + t14 * qJD(5) + ((t143 * t167 + t144 * t168 + t240 * t376) * t444 + (-pkin(3) * t377 + t228 * t235 + t229 * t236) * t445) * t446 + t350 + t461 * qJD(3), t7 * qJD(3) + (-t201 * mrSges(5,1) - t200 * mrSges(5,2) + t183 * t349 + t274 + t331 - t332) * qJD(4) + t31 * qJD(5) + t352, t14 * qJD(3) + t31 * qJD(4) + (t379 - t413 - t414) * qJD(5) + t351; -qJD(2) * t6 - qJD(4) * t4 - qJD(5) * t11 + t392, -qJD(4) * t8 + qJD(5) * t15 - t350, -qJD(4) * t19 + qJD(5) * t23, (-t236 * mrSges(5,1) - t235 * mrSges(5,2) + t214 * t349 + t274 + t329 - t330) * qJD(4) + t36 * qJD(5) + t326, t36 * qJD(4) + (t379 - t411 - t412) * qJD(5) + t324; qJD(2) * t2 + qJD(3) * t4 - qJD(5) * t24 + t400, qJD(3) * t8 + qJD(5) * t32 - t352, qJD(5) * t39 - t326, -t239, -t239 - t304; qJD(2) * t9 + qJD(3) * t11 + qJD(4) * t24 - t393, -qJD(3) * t15 - qJD(4) * t32 - t351, -qJD(4) * t39 - t324, t304, 0;];
Cq = t25;
