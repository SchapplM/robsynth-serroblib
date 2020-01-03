% Calculate time derivative of joint inertia matrix for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:23
% EndTime: 2019-12-31 17:27:41
% DurationCPUTime: 9.95s
% Computational Cost: add. (21605->729), mult. (35520->1029), div. (0->0), fcn. (34578->8), ass. (0->377)
t286 = sin(qJ(2));
t473 = -qJD(1) * t286 / 0.2e1;
t290 = cos(qJ(1));
t289 = cos(qJ(2));
t396 = qJD(2) * t289;
t365 = t396 / 0.2e1;
t287 = sin(qJ(1));
t401 = qJD(1) * t287;
t377 = t286 * t401;
t472 = t290 * t365 - t377 / 0.2e1;
t399 = qJD(1) * t290;
t366 = t399 / 0.2e1;
t471 = t286 * t366 + t287 * t365;
t370 = t287 * t396;
t300 = t286 * t399 + t370;
t443 = Icges(3,4) * t286;
t335 = Icges(3,1) * t289 - t443;
t229 = Icges(3,5) * t287 + t290 * t335;
t424 = t229 * t289;
t442 = Icges(3,4) * t289;
t331 = -Icges(3,2) * t286 + t442;
t226 = Icges(3,6) * t287 + t290 * t331;
t430 = t226 * t286;
t315 = -t424 + t430;
t470 = t287 * t315;
t291 = -pkin(7) - pkin(6);
t448 = pkin(6) + t291;
t469 = t289 * t448;
t228 = -Icges(3,5) * t290 + t287 * t335;
t426 = t228 * t289;
t225 = -Icges(3,6) * t290 + t287 * t331;
t432 = t225 * t286;
t316 = -t426 + t432;
t468 = t290 * t316;
t418 = t286 * t290;
t467 = -rSges(3,2) * t418 + t287 * rSges(3,3);
t327 = Icges(3,5) * t289 - Icges(3,6) * t286;
t222 = -Icges(3,3) * t290 + t287 * t327;
t400 = qJD(1) * t289;
t356 = -qJD(3) + t400;
t395 = qJD(2) * t290;
t372 = t286 * t395;
t466 = t287 * t356 + t372;
t283 = qJD(3) + qJD(4);
t361 = -t283 + t400;
t465 = t287 * t361 + t372;
t397 = qJD(2) * t287;
t373 = t286 * t397;
t464 = t290 * t361 - t373;
t284 = qJ(3) + qJ(4);
t277 = sin(t284);
t278 = cos(t284);
t415 = t289 * t290;
t235 = -t277 * t415 + t278 * t287;
t236 = t277 * t287 + t278 * t415;
t171 = t236 * rSges(5,1) + t235 * rSges(5,2) + rSges(5,3) * t418;
t273 = pkin(2) * t415;
t247 = pkin(6) * t418 + t273;
t285 = sin(qJ(3));
t422 = t285 * t287;
t271 = pkin(3) * t422;
t288 = cos(qJ(3));
t274 = pkin(3) * t288 + pkin(2);
t309 = t274 * t415 - t291 * t418 + t271;
t188 = t309 - t247;
t410 = t171 + t188;
t417 = t287 * t289;
t233 = -t277 * t417 - t278 * t290;
t234 = -t277 * t290 + t278 * t417;
t344 = -t234 * rSges(5,1) - t233 * rSges(5,2);
t419 = t286 * t287;
t170 = rSges(5,3) * t419 - t344;
t449 = pkin(2) - t274;
t298 = -t286 * t448 - t289 * t449;
t421 = t285 * t290;
t392 = pkin(3) * t421;
t187 = t287 * t298 - t392;
t411 = t170 + t187;
t463 = -t287 * t411 - t290 * t410;
t462 = 2 * m(3);
t461 = 2 * m(4);
t460 = 2 * m(5);
t459 = t287 ^ 2;
t458 = t290 ^ 2;
t457 = t287 / 0.2e1;
t456 = -t289 / 0.2e1;
t455 = t290 / 0.2e1;
t454 = -rSges(4,3) - pkin(6);
t259 = rSges(3,1) * t286 + rSges(3,2) * t289;
t453 = m(3) * t259;
t452 = pkin(2) * t289;
t451 = pkin(3) * t285;
t280 = t287 * pkin(5);
t450 = qJD(1) / 0.2e1;
t447 = rSges(3,3) * t290;
t446 = rSges(5,3) * t286;
t445 = -rSges(5,3) + t291;
t362 = -t283 * t289 + qJD(1);
t314 = t287 * t362;
t151 = -t277 * t464 + t278 * t314;
t152 = t277 * t314 + t278 * t464;
t345 = t152 * rSges(5,1) + t151 * rSges(5,2);
t104 = rSges(5,3) * t300 + t345;
t369 = t289 * t395;
t444 = t104 * t418 + t170 * t369;
t441 = Icges(4,4) * t285;
t440 = Icges(4,4) * t288;
t439 = Icges(5,4) * t277;
t438 = Icges(5,4) * t278;
t416 = t288 * t290;
t242 = -t285 * t417 - t416;
t243 = t288 * t417 - t421;
t347 = -rSges(4,1) * t243 - rSges(4,2) * t242;
t185 = rSges(4,3) * t419 - t347;
t434 = t185 * t290;
t332 = Icges(5,1) * t278 - t439;
t217 = -Icges(5,5) * t289 + t286 * t332;
t433 = t217 * t278;
t431 = t225 * t289;
t429 = t226 * t289;
t333 = Icges(4,1) * t288 - t441;
t227 = -Icges(4,5) * t289 + t286 * t333;
t428 = t227 * t288;
t427 = t228 * t286;
t425 = t229 * t286;
t423 = t283 * t286;
t420 = t286 * t274;
t313 = t290 * t362;
t149 = t277 * t465 + t278 * t313;
t150 = t277 * t313 - t278 * t465;
t384 = t150 * rSges(5,1) + t149 * rSges(5,2) + rSges(5,3) * t369;
t103 = -rSges(5,3) * t377 + t384;
t266 = pkin(6) * t369;
t393 = qJD(3) * t288;
t389 = pkin(3) * t393;
t379 = qJD(1) * t392 + t287 * t389 + t291 * t377;
t390 = qJD(3) * t451;
t414 = t103 - t266 + (pkin(6) * t401 + t395 * t449) * t286 + ((-qJD(2) * t291 - t390) * t290 + t449 * t401) * t289 + t379;
t343 = rSges(5,1) * t278 - rSges(5,2) * t277;
t218 = -rSges(5,3) * t289 + t286 * t343;
t398 = qJD(2) * t286;
t413 = t171 * t398 + t218 * t377;
t161 = (-rSges(5,1) * t277 - rSges(5,2) * t278) * t423 + (t289 * t343 + t446) * qJD(2);
t394 = qJD(3) * t286;
t371 = t285 * t394;
t199 = -pkin(3) * t371 + qJD(2) * t298;
t412 = -t161 - t199;
t244 = -t285 * t415 + t287 * t288;
t245 = t288 * t415 + t422;
t186 = t245 * rSges(4,1) + t244 * rSges(4,2) + rSges(4,3) * t418;
t409 = -t186 - t247;
t346 = rSges(4,1) * t288 - rSges(4,2) * t285;
t198 = (-rSges(4,1) * t285 - rSges(4,2) * t288) * t394 + (rSges(4,3) * t286 + t289 * t346) * qJD(2);
t351 = pkin(6) * t286 + t452;
t253 = t351 * qJD(2);
t408 = -t198 - t253;
t132 = t289 * t170 + t218 * t419;
t214 = -t286 * t449 + t469;
t407 = t214 + t218;
t230 = -rSges(4,3) * t289 + t286 * t346;
t260 = t286 * pkin(2) - t289 * pkin(6);
t406 = -t230 - t260;
t246 = t351 * t287;
t405 = t287 * t246 + t290 * t247;
t404 = rSges(3,2) * t377 + rSges(3,3) * t399;
t403 = t290 * pkin(1) + t280;
t223 = Icges(3,3) * t287 + t290 * t327;
t402 = qJD(1) * t223;
t357 = -qJD(3) * t289 + qJD(1);
t311 = t357 * t288;
t177 = t287 * t311 + (-t290 * t356 + t373) * t285;
t312 = t357 * t285;
t178 = t356 * t416 + (-t288 * t398 + t312) * t287;
t113 = Icges(4,5) * t178 + Icges(4,6) * t177 + Icges(4,3) * t300;
t115 = Icges(4,4) * t178 + Icges(4,2) * t177 + Icges(4,6) * t300;
t117 = Icges(4,1) * t178 + Icges(4,4) * t177 + Icges(4,5) * t300;
t179 = Icges(4,5) * t243 + Icges(4,6) * t242 + Icges(4,3) * t419;
t181 = Icges(4,4) * t243 + Icges(4,2) * t242 + Icges(4,6) * t419;
t183 = Icges(4,1) * t243 + Icges(4,4) * t242 + Icges(4,5) * t419;
t322 = -t181 * t285 + t183 * t288;
t31 = (qJD(2) * t322 - t113) * t289 + (qJD(2) * t179 - t115 * t285 + t117 * t288 + (-t181 * t288 - t183 * t285) * qJD(3)) * t286;
t326 = Icges(4,5) * t288 - Icges(4,6) * t285;
t189 = (-Icges(4,5) * t285 - Icges(4,6) * t288) * t394 + (Icges(4,3) * t286 + t289 * t326) * qJD(2);
t329 = -Icges(4,2) * t285 + t440;
t192 = (-Icges(4,2) * t288 - t441) * t394 + (Icges(4,6) * t286 + t289 * t329) * qJD(2);
t195 = (-Icges(4,1) * t285 - t440) * t394 + (Icges(4,5) * t286 + t289 * t333) * qJD(2);
t221 = -Icges(4,3) * t289 + t286 * t326;
t224 = -Icges(4,6) * t289 + t286 * t329;
t60 = t177 * t224 + t178 * t227 + t189 * t419 + t192 * t242 + t195 * t243 + t221 * t300;
t388 = t31 / 0.2e1 + t60 / 0.2e1;
t175 = t285 * t466 + t290 * t311;
t176 = -t288 * t466 + t290 * t312;
t299 = t369 - t377;
t112 = Icges(4,5) * t176 + Icges(4,6) * t175 + Icges(4,3) * t299;
t114 = Icges(4,4) * t176 + Icges(4,2) * t175 + Icges(4,6) * t299;
t116 = Icges(4,1) * t176 + Icges(4,4) * t175 + Icges(4,5) * t299;
t180 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t418;
t182 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t418;
t184 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t418;
t321 = -t182 * t285 + t184 * t288;
t32 = (qJD(2) * t321 - t112) * t289 + (qJD(2) * t180 - t114 * t285 + t116 * t288 + (-t182 * t288 - t184 * t285) * qJD(3)) * t286;
t59 = t175 * t224 + t176 * t227 + t189 * t418 + t192 * t244 + t195 * t245 + t221 * t299;
t387 = t32 / 0.2e1 + t59 / 0.2e1;
t123 = t221 * t419 + t224 * t242 + t227 * t243;
t95 = -t179 * t289 + t286 * t322;
t386 = t123 / 0.2e1 + t95 / 0.2e1;
t124 = t221 * t418 + t224 * t244 + t227 * t245;
t96 = -t180 * t289 + t286 * t321;
t385 = t96 / 0.2e1 + t124 / 0.2e1;
t383 = -t253 + t412;
t382 = t176 * rSges(4,1) + t175 * rSges(4,2) + rSges(4,3) * t369;
t265 = pkin(2) * t373;
t301 = -t287 * t400 - t372;
t381 = t287 * (pkin(6) * t300 + qJD(1) * t273 - t265) + t290 * (pkin(2) * t301 - pkin(6) * t377 + t266) + t246 * t399;
t380 = -t260 - t407;
t378 = t230 * t401;
t328 = -Icges(5,2) * t277 + t438;
t216 = -Icges(5,6) * t289 + t286 * t328;
t375 = t216 * t396;
t374 = t224 * t396;
t368 = t419 / 0.2e1;
t367 = t418 / 0.2e1;
t364 = -t274 * t289 - pkin(1);
t363 = t407 * t290;
t203 = t406 * t290;
t360 = t289 * t390;
t359 = t290 * t389;
t358 = t289 * t104 + t161 * t419 + t300 * t218;
t143 = t380 * t290;
t325 = Icges(5,5) * t278 - Icges(5,6) * t277;
t215 = -Icges(5,3) * t289 + t286 * t325;
t129 = -t215 * t289 + (-t216 * t277 + t433) * t286;
t164 = Icges(5,5) * t234 + Icges(5,6) * t233 + Icges(5,3) * t419;
t166 = Icges(5,4) * t234 + Icges(5,2) * t233 + Icges(5,6) * t419;
t168 = Icges(5,1) * t234 + Icges(5,4) * t233 + Icges(5,5) * t419;
t324 = -t166 * t277 + t168 * t278;
t83 = -t164 * t289 + t286 * t324;
t165 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t418;
t167 = Icges(5,4) * t236 + Icges(5,2) * t235 + Icges(5,6) * t418;
t169 = Icges(5,1) * t236 + Icges(5,4) * t235 + Icges(5,5) * t418;
t323 = -t167 * t277 + t169 * t278;
t84 = -t165 * t289 + t286 * t323;
t339 = t83 * t287 + t84 * t290;
t110 = t215 * t419 + t216 * t233 + t217 * t234;
t73 = t164 * t419 + t166 * t233 + t168 * t234;
t74 = t165 * t419 + t167 * t233 + t169 * t234;
t342 = t287 * t73 + t290 * t74;
t40 = -t110 * t289 + t286 * t342;
t111 = t215 * t418 + t216 * t235 + t217 * t236;
t75 = t164 * t418 + t166 * t235 + t168 * t236;
t76 = t165 * t418 + t167 * t235 + t169 * t236;
t341 = t287 * t75 + t290 * t76;
t41 = -t111 * t289 + t286 * t341;
t100 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t300;
t102 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t300;
t98 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t300;
t19 = t100 * t235 + t102 * t236 + t149 * t166 + t150 * t168 + t164 * t299 + t418 * t98;
t101 = Icges(5,1) * t150 + Icges(5,4) * t149 + Icges(5,5) * t299;
t97 = Icges(5,5) * t150 + Icges(5,6) * t149 + Icges(5,3) * t299;
t99 = Icges(5,4) * t150 + Icges(5,2) * t149 + Icges(5,6) * t299;
t20 = t101 * t236 + t149 * t167 + t150 * t169 + t165 * t299 + t235 * t99 + t418 * t97;
t157 = (-Icges(5,5) * t277 - Icges(5,6) * t278) * t423 + (Icges(5,3) * t286 + t289 * t325) * qJD(2);
t158 = (-Icges(5,2) * t278 - t439) * t423 + (Icges(5,6) * t286 + t289 * t328) * qJD(2);
t159 = (-Icges(5,1) * t277 - t438) * t423 + (Icges(5,5) * t286 + t289 * t332) * qJD(2);
t49 = t149 * t216 + t150 * t217 + t157 * t418 + t158 * t235 + t159 * t236 + t215 * t299;
t57 = t287 * t76 - t290 * t75;
t5 = (qJD(2) * t341 - t49) * t289 + (-qJD(1) * t57 + qJD(2) * t111 + t19 * t287 + t20 * t290) * t286;
t21 = t100 * t233 + t102 * t234 + t151 * t166 + t152 * t168 + t164 * t300 + t419 * t98;
t22 = t101 * t234 + t151 * t167 + t152 * t169 + t165 * t300 + t233 * t99 + t419 * t97;
t50 = t151 * t216 + t152 * t217 + t157 * t419 + t158 * t233 + t159 * t234 + t215 * t300;
t56 = t287 * t74 - t290 * t73;
t6 = (qJD(2) * t342 - t50) * t289 + (-qJD(1) * t56 + qJD(2) * t110 + t21 * t287 + t22 * t290) * t286;
t350 = t41 * t369 + t5 * t418 + t6 * t419 + (-t129 * t289 + t286 * t339) * t398 + t300 * t40;
t349 = rSges(3,1) * t289 - rSges(3,2) * t286;
t348 = rSges(4,1) * t178 + rSges(4,2) * t177;
t340 = t287 * t84 - t290 * t83;
t85 = t179 * t419 + t181 * t242 + t183 * t243;
t86 = t180 * t419 + t182 * t242 + t184 * t243;
t62 = t287 * t86 - t290 * t85;
t338 = t287 * t85 + t290 * t86;
t87 = t179 * t418 + t181 * t244 + t183 * t245;
t88 = t180 * t418 + t182 * t244 + t184 * t245;
t63 = t287 * t88 - t290 * t87;
t337 = t287 * t87 + t290 * t88;
t336 = t287 * t95 + t290 * t96;
t330 = Icges(3,2) * t289 + t443;
t320 = -t186 * t287 + t434;
t319 = -t185 * t287 - t186 * t290;
t232 = rSges(3,1) * t415 + t467;
t310 = -pkin(1) - t349;
t308 = qJD(2) * t259;
t307 = t187 * t290 - t287 * t410;
t306 = t286 * t454 - pkin(1) - t452;
t304 = qJD(2) * t330;
t303 = qJD(2) * (-Icges(3,5) * t286 - Icges(3,6) * t289);
t302 = t286 * t445 + t364;
t297 = t306 * t287;
t12 = qJD(1) * t341 - t19 * t290 + t20 * t287;
t13 = qJD(1) * t342 - t21 * t290 + t22 * t287;
t25 = (qJD(2) * t324 - t98) * t289 + (qJD(2) * t164 + (-t166 * t283 + t102) * t278 + (-t168 * t283 - t100) * t277) * t286;
t26 = (qJD(2) * t323 - t97) * t289 + (qJD(2) * t165 + (-t167 * t283 + t101) * t278 + (-t169 * t283 - t99) * t277) * t286;
t296 = -t290 * t6 / 0.2e1 + t5 * t457 + t13 * t368 + t12 * t367 + (qJD(1) * t339 - t25 * t290 + t26 * t287) * t456 + t40 * t401 / 0.2e1 + t41 * t366 + t340 * t398 / 0.2e1 + t472 * t57 + t471 * t56;
t295 = -t289 * t157 + t215 * t398 + t396 * t433 + (t159 * t286 - t216 * t423) * t278;
t294 = -t289 * t189 + t221 * t398 + t396 * t428 + (t195 * t288 - t224 * t393) * t286;
t128 = t129 * t398;
t61 = (-t375 + (-t217 * t283 - t158) * t286) * t277 + t295;
t7 = t128 + (qJD(2) * t339 - t61) * t289 + (-qJD(1) * t340 + t25 * t287 + t26 * t290) * t286;
t293 = -t289 * t7 - t377 * t41 + t350;
t292 = t128 + (t25 + t50) * t368 + (t26 + t49) * t367 + (t111 + t84) * t472 + (t110 + t83) * t471;
t281 = t290 * pkin(5);
t276 = pkin(5) * t399;
t252 = t349 * qJD(2);
t248 = t260 * t401;
t231 = t287 * t349 - t447;
t211 = t232 + t403;
t210 = t287 * t310 + t281 + t447;
t202 = t406 * t287;
t191 = t287 * t303 + t402;
t190 = -qJD(1) * t222 + t290 * t303;
t163 = t259 * t397 + ((-rSges(3,3) - pkin(5)) * t287 + t310 * t290) * qJD(1);
t162 = rSges(3,1) * t301 - rSges(3,2) * t369 - pkin(1) * t401 + t276 + t404;
t155 = t170 * t418;
t145 = t403 - t409;
t144 = t281 + t297 + t347;
t142 = t380 * t287;
t141 = -t186 * t289 - t230 * t418;
t140 = t185 * t289 + t230 * t419;
t139 = t223 * t287 - t290 * t315;
t138 = t222 * t287 - t468;
t137 = -t223 * t290 - t470;
t136 = -t222 * t290 - t287 * t316;
t135 = -t221 * t289 + (-t224 * t285 + t428) * t286;
t134 = t135 * t398;
t133 = -t171 * t289 - t218 * t418;
t131 = t309 + t171 + t403;
t130 = t287 * t302 + t281 + t344 + t392;
t127 = t320 * t286;
t126 = qJD(1) * t203 + t287 * t408;
t125 = t290 * t408 + t248 + t378;
t122 = -t359 + t265 + (-t360 + (-t420 - t469) * qJD(2)) * t287 + (t290 * t298 + t271) * qJD(1);
t120 = -t171 * t419 + t155;
t119 = rSges(4,3) * t300 + t348;
t118 = -rSges(4,3) * t377 + t382;
t105 = -t319 + t405;
t92 = t265 + t454 * t370 + (t290 * t306 - t280) * qJD(1) - t348;
t91 = -pkin(2) * t372 + qJD(1) * t297 + t266 + t276 + t382;
t90 = -t286 * t363 - t289 * t410;
t89 = t187 * t289 + t214 * t419 + t132;
t78 = qJD(1) * t143 + t287 * t383;
t77 = t290 * t383 + t401 * t407 + t248;
t72 = t286 * t307 + t155;
t71 = t359 + (t360 + (t289 * t445 + t420) * qJD(2)) * t287 + ((-pkin(5) - t451) * t287 + t302 * t290) * qJD(1) - t345;
t70 = t276 + (-t360 + (-t289 * t291 - t420) * qJD(2)) * t290 + (t364 - t446) * t401 + t379 + t384;
t69 = t405 - t463;
t68 = (t230 * t397 + t119) * t289 + (-qJD(2) * t185 + t198 * t287 + t230 * t399) * t286;
t67 = (-t230 * t395 - t118) * t289 + (qJD(2) * t186 - t198 * t290 + t378) * t286;
t66 = (-t374 + (-qJD(3) * t227 - t192) * t286) * t285 + t294;
t65 = -t170 * t398 + t358;
t64 = -t161 * t418 + (-t218 * t395 - t103) * t289 + t413;
t51 = t320 * t396 + (qJD(1) * t319 - t118 * t287 + t119 * t290) * t286;
t48 = t118 * t290 + t119 * t287 + (t287 * t409 + t434) * qJD(1) + t381;
t45 = -t124 * t289 + t286 * t337;
t44 = -t123 * t289 + t286 * t338;
t42 = -t171 * t370 + (-t103 * t287 + (-t170 * t287 - t171 * t290) * qJD(1)) * t286 + t444;
t34 = (t214 * t397 + t122) * t289 + (-qJD(2) * t411 + t199 * t287 + t214 * t399) * t286 + t358;
t33 = (-qJD(2) * t363 - t414) * t289 + (qJD(2) * t188 + t214 * t401 + t290 * t412) * t286 + t413;
t30 = t112 * t419 + t114 * t242 + t116 * t243 + t177 * t182 + t178 * t184 + t180 * t300;
t29 = t113 * t419 + t115 * t242 + t117 * t243 + t177 * t181 + t178 * t183 + t179 * t300;
t28 = t112 * t418 + t114 * t244 + t116 * t245 + t175 * t182 + t176 * t184 + t180 * t299;
t27 = t113 * t418 + t115 * t244 + t117 * t245 + t175 * t181 + t176 * t183 + t179 * t299;
t18 = t414 * t290 + (t104 + t122) * t287 + (t411 * t290 + (-t247 - t410) * t287) * qJD(1) + t381;
t17 = t307 * t396 + (qJD(1) * t463 + t122 * t290 - t414 * t287) * t286 + t444;
t16 = qJD(1) * t338 + t287 * t30 - t29 * t290;
t15 = qJD(1) * t337 - t27 * t290 + t28 * t287;
t9 = (qJD(2) * t338 - t60) * t289 + (-qJD(1) * t62 + qJD(2) * t123 + t287 * t29 + t290 * t30) * t286;
t8 = (qJD(2) * t337 - t59) * t289 + (-t63 * qJD(1) + qJD(2) * t124 + t27 * t287 + t28 * t290) * t286;
t1 = [t294 + t295 - t227 * t371 + (t130 * t71 + t131 * t70) * t460 + (t144 * t92 + t145 * t91) * t461 + (t162 * t211 + t163 * t210) * t462 + (-t330 + t335) * t398 + (Icges(3,1) * t286 + t331 + t442) * t396 + (-t192 * t286 - t374) * t285 + (-t158 * t286 - t217 * t423 - t375) * t277; m(4) * (t125 * t144 + t126 * t145 + t202 * t91 + t203 * t92) + m(5) * (t130 * t77 + t131 * t78 + t142 * t70 + t143 * t71) + (t459 / 0.2e1 + t458 / 0.2e1) * t327 * qJD(2) + (-t25 / 0.2e1 + (qJD(1) * t226 - t287 * t304) * t456 + t229 * t473 - t50 / 0.2e1 + (t432 / 0.2e1 - t426 / 0.2e1) * qJD(2) - t388 + m(3) * (-t163 * t259 - t210 * t252) + (t84 / 0.2e1 + t429 / 0.2e1 + t425 / 0.2e1 - t211 * t453 + t111 / 0.2e1 + t385) * qJD(1)) * t290 + ((-qJD(1) * t225 - t290 * t304) * t289 / 0.2e1 + t228 * t473 + t49 / 0.2e1 + t26 / 0.2e1 + (-t430 / 0.2e1 + t424 / 0.2e1) * qJD(2) + t387 + m(3) * (-t162 * t259 - t211 * t252) + (t431 / 0.2e1 + t427 / 0.2e1 + t110 / 0.2e1 + t83 / 0.2e1 + t210 * t453 + t386) * qJD(1)) * t287; (t142 * t78 + t143 * t77 + t18 * t69) * t460 - t290 * t13 + t287 * t12 - t290 * t16 + t287 * t15 + (t105 * t48 + t125 * t203 + t126 * t202) * t461 + ((t231 * t287 + t232 * t290) * ((qJD(1) * t231 - t290 * t308 + t404) * t290 + (-t287 * t308 + (-t232 + t467) * qJD(1)) * t287) + (t458 + t459) * t259 * t252) * t462 + t287 * ((t190 * t287 + (t138 + t470) * qJD(1)) * t287 + (t139 * qJD(1) + (t225 * t396 + t228 * t398) * t290 + (-t191 + (-t425 - t429) * qJD(2) + (t223 - t316) * qJD(1)) * t287) * t290) - t290 * ((t191 * t290 + (t137 + t468) * qJD(1)) * t290 + (t136 * qJD(1) + (-t226 * t396 - t229 * t398 + t402) * t287 + (-t190 + (t427 + t431) * qJD(2) - t315 * qJD(1)) * t290) * t287) + (-t136 * t290 + t137 * t287 + t56 + t62) * t401 + (-t138 * t290 + t139 * t287 + t57 + t63) * t399; (t387 * t290 + t388 * t287 + (-t287 * t385 + t290 * t386) * qJD(1)) * t286 + (-t61 - t66 + (t287 * t386 + t290 * t385) * qJD(2)) * t289 + t134 + m(4) * (t140 * t92 + t141 * t91 + t144 * t68 + t145 * t67) + m(5) * (t130 * t34 + t131 * t33 + t70 * t90 + t71 * t89) + t292; (qJD(2) * (t287 * t96 - t290 * t95) / 0.2e1 + t15 * t455 + t16 * t457 + (t62 * t455 - t287 * t63 / 0.2e1) * qJD(1)) * t286 + (t45 * t450 + t63 * t365 + (qJD(1) * t96 - t31) * t456 - t9 / 0.2e1) * t290 + (t44 * t450 + t62 * t365 + (qJD(1) * t95 + t32) * t456 + t8 / 0.2e1) * t287 + m(4) * (t105 * t51 + t125 * t140 + t126 * t141 + t127 * t48 + t202 * t67 + t203 * t68) + m(5) * (t142 * t33 + t143 * t34 + t17 * t69 + t18 * t72 + t77 * t89 + t78 * t90) + t296; (t17 * t72 + t33 * t90 + t34 * t89) * t460 + (t127 * t51 + t140 * t68 + t141 * t67) * t461 + (t66 * t289 - t134 - t7 + (t287 * t44 - t289 * t336 + t290 * t45) * qJD(2)) * t289 + (t290 * t8 + t287 * t9 - t289 * (t287 * t31 + t290 * t32) + (-t135 * t289 + t286 * t336) * qJD(2) + ((-t289 * t95 + t44) * t290 + (t289 * t96 - t41 - t45) * t287) * qJD(1)) * t286 + t350; t292 - t61 * t289 + m(5) * (t130 * t65 + t131 * t64 + t132 * t71 + t133 * t70); m(5) * (t120 * t18 + t132 * t77 + t133 * t78 + t142 * t64 + t143 * t65 + t42 * t69) + t296; m(5) * (t120 * t17 + t132 * t34 + t133 * t33 + t42 * t72 + t64 * t90 + t65 * t89) + t293; (t120 * t42 + t132 * t65 + t133 * t64) * t460 + t293;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
