% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:27
% EndTime: 2019-03-09 18:40:13
% DurationCPUTime: 21.91s
% Computational Cost: add. (38115->735), mult. (98935->1010), div. (0->0), fcn. (80419->12), ass. (0->321)
t320 = sin(qJ(6));
t324 = cos(qJ(6));
t318 = sin(pkin(6));
t325 = cos(qJ(2));
t409 = qJD(1) * t325;
t387 = t318 * t409;
t298 = -qJD(3) + t387;
t346 = -qJD(5) + t298;
t338 = t324 * t346;
t322 = sin(qJ(3));
t319 = cos(pkin(6));
t411 = qJD(1) * t319;
t365 = qJD(2) + t411;
t323 = sin(qJ(2));
t410 = qJD(1) * t323;
t388 = t318 * t410;
t458 = cos(qJ(3));
t238 = t322 * t388 - t365 * t458;
t391 = t318 * t458;
t359 = qJD(1) * t391;
t240 = t322 * t365 + t323 * t359;
t317 = sin(pkin(12));
t439 = cos(pkin(12));
t192 = -t238 * t439 - t240 * t317;
t321 = sin(qJ(5));
t341 = -t238 * t317 + t240 * t439;
t457 = cos(qJ(5));
t476 = t192 * t321 + t341 * t457;
t120 = t320 * t476 + t338;
t403 = qJD(6) * t324;
t408 = qJD(2) * t323;
t386 = t318 * t408;
t302 = qJD(1) * t386;
t122 = -t320 * t346 + t324 * t476;
t405 = qJD(6) * t122;
t352 = t325 * t359;
t477 = qJD(3) * t238;
t210 = -qJD(2) * t352 + t477;
t401 = qJD(1) * qJD(2);
t377 = t325 * t401;
t358 = t318 * t377;
t211 = qJD(3) * t240 + t322 * t358;
t156 = -t210 * t317 + t211 * t439;
t157 = -t210 * t439 - t211 * t317;
t189 = t192 * t457;
t406 = qJD(5) * t321;
t72 = -qJD(5) * t189 + t156 * t321 - t157 * t457 + t341 * t406;
t48 = -t302 * t324 - t320 * t72 + t405;
t443 = -t120 * t403 - t320 * t48;
t131 = -t321 * t341 + t189;
t492 = qJD(6) - t131;
t502 = t492 * t320;
t469 = t122 * t502;
t404 = qJD(6) * t320;
t47 = qJD(6) * t338 - t302 * t320 + t324 * t72 + t404 * t476;
t494 = t131 * t324;
t507 = t120 * t494 - t324 * t47 + t443 - t469;
t397 = pkin(1) * t411;
t305 = t323 * t397;
t250 = pkin(8) * t358 + qJD(2) * t305;
t181 = pkin(3) * t211 + t250;
t117 = pkin(4) * t156 + t181;
t370 = t156 * t457 + t321 * t157;
t483 = t476 * qJD(5);
t73 = t370 + t483;
t22 = pkin(5) * t73 + pkin(11) * t72 + t117;
t472 = pkin(10) * t341;
t268 = pkin(8) * t387 + t305;
t226 = pkin(9) * t365 + t268;
t257 = (-pkin(2) * t325 - pkin(9) * t323 - pkin(1)) * t318;
t234 = qJD(1) * t257;
t184 = -t226 * t322 + t234 * t458;
t154 = -qJ(4) * t240 + t184;
t140 = -pkin(3) * t298 + t154;
t185 = t226 * t458 + t234 * t322;
t155 = -qJ(4) * t238 + t185;
t148 = t317 * t155;
t98 = t140 * t439 - t148;
t76 = -pkin(4) * t298 - t472 + t98;
t479 = t192 * pkin(10);
t150 = t439 * t155;
t99 = t140 * t317 + t150;
t78 = t99 + t479;
t35 = t321 * t76 + t457 * t78;
t33 = -pkin(11) * t346 + t35;
t265 = -pkin(8) * t388 + t325 * t397;
t225 = -pkin(2) * t365 - t265;
t190 = t238 * pkin(3) + qJD(4) + t225;
t135 = -pkin(4) * t192 + t190;
t63 = -pkin(5) * t131 - pkin(11) * t476 + t135;
t354 = t320 * t33 - t324 * t63;
t380 = qJD(5) * t457;
t345 = t318 * (pkin(2) * t323 - pkin(9) * t325);
t267 = qJD(2) * t345;
t248 = qJD(1) * t267;
t419 = t318 * t323;
t306 = pkin(8) * t419;
t417 = t319 * t325;
t279 = pkin(1) * t417 - t306;
t269 = t279 * qJD(2);
t249 = qJD(1) * t269;
t134 = -qJD(3) * t185 + t248 * t458 - t322 * t249;
t100 = pkin(3) * t302 + t210 * qJ(4) - t240 * qJD(4) + t134;
t381 = qJD(3) * t458;
t407 = qJD(3) * t322;
t133 = -t226 * t407 + t234 * t381 + t248 * t322 + t249 * t458;
t105 = -qJ(4) * t211 - qJD(4) * t238 + t133;
t57 = t100 * t439 - t105 * t317;
t40 = pkin(4) * t302 - pkin(10) * t157 + t57;
t58 = t100 * t317 + t105 * t439;
t42 = -pkin(10) * t156 + t58;
t374 = -t321 * t40 - t380 * t76 + t406 * t78 - t42 * t457;
t7 = pkin(11) * t302 - t374;
t2 = -qJD(6) * t354 + t22 * t320 + t324 * t7;
t448 = t492 * t354;
t506 = t2 + t448;
t266 = qJD(1) * t345;
t200 = -t265 * t322 + t266 * t458;
t390 = t458 * qJ(4);
t300 = pkin(9) * t458 + t390;
t412 = qJD(1) * t318;
t505 = (pkin(3) * t323 - t325 * t390) * t412 + t200 + qJD(3) * t300 + t322 * qJD(4);
t201 = t265 * t458 + t266 * t322;
t299 = (-qJ(4) - pkin(9)) * t322;
t361 = t322 * t387;
t504 = qJ(4) * t361 + qJD(3) * t299 + qJD(4) * t458 - t201;
t46 = t48 * t324;
t503 = t120 * t502 - t46;
t282 = t317 * t458 + t322 * t439;
t461 = t298 * t282;
t336 = t317 * t322 - t439 * t458;
t329 = t336 * t412;
t221 = t325 * t329;
t493 = -qJD(3) * t336 + t221;
t433 = t131 ^ 2;
t434 = t476 ^ 2;
t501 = -t433 + t434;
t16 = t320 * t63 + t324 * t33;
t3 = -qJD(6) * t16 + t22 * t324 - t320 * t7;
t500 = t16 * t492 + t3;
t44 = t47 * t320;
t498 = -t44 + (t403 - t494) * t122;
t436 = t122 * t476;
t69 = t320 * t73;
t442 = t403 * t492 + t69;
t497 = -t492 * t494 - t436 + t442;
t34 = -t321 * t78 + t457 * t76;
t32 = pkin(5) * t346 - t34;
t496 = t131 * t32;
t468 = t131 * t346;
t495 = -t72 + t468;
t430 = t476 * t131;
t464 = t317 * t504 + t439 * t505;
t463 = -t317 * t505 + t439 * t504;
t85 = pkin(5) * t476 - pkin(11) * t131;
t491 = -t135 * t131 + t374;
t311 = pkin(3) * t439 + pkin(4);
t455 = pkin(3) * t317;
t272 = t311 * t321 + t455 * t457;
t103 = -t154 * t317 - t150;
t335 = t103 - t479;
t104 = t154 * t439 - t148;
t84 = t104 - t472;
t440 = qJD(5) * t272 - t321 * t84 + t335 * t457;
t490 = t35 + t440;
t489 = pkin(4) * t388 + pkin(10) * t493 + t464;
t488 = -pkin(10) * t461 - t463;
t344 = t476 * t298;
t487 = -t344 - t370;
t438 = t120 * t476;
t484 = t492 * t476;
t330 = t457 * t336;
t400 = qJD(3) + qJD(5);
t416 = -t221 * t457 + t282 * t406 - t321 * t461 + t330 * t400;
t334 = t321 * t336;
t415 = -qJD(5) * t334 + t282 * t380 + t321 * t493 - t457 * t461;
t462 = -t268 + (-t361 + t407) * pkin(3);
t376 = t321 * t42 - t40 * t457;
t10 = -qJD(5) * t35 - t376;
t8 = -pkin(5) * t302 - t10;
t482 = t16 * t476 + t32 * t403 + t320 * t8;
t481 = -t135 * t476 - t376;
t480 = t32 * t404 - t324 * t8 + t354 * t476;
t223 = t299 * t439 - t300 * t317;
t204 = -pkin(10) * t282 + t223;
t224 = t299 * t317 + t300 * t439;
t205 = -pkin(10) * t336 + t224;
t451 = -t204 * t380 + t205 * t406 + t321 * t489 + t457 * t488;
t478 = t341 * t192;
t414 = -pkin(4) * t461 + t462;
t475 = pkin(11) * t388 + t451;
t474 = pkin(5) * t415 + pkin(11) * t416 + t414;
t314 = t318 ^ 2;
t473 = -0.2e1 * t314 * t401;
t145 = t204 * t321 + t205 * t457;
t450 = -qJD(5) * t145 + t321 * t488 - t457 * t489;
t271 = t311 * t457 - t321 * t455;
t253 = t271 * qJD(5);
t37 = t321 * t335 + t457 * t84;
t441 = t253 - t37;
t71 = t324 * t73;
t470 = t404 * t492 - t71;
t465 = t185 * t298 - t134;
t418 = t318 * t325;
t280 = pkin(1) * t319 * t323 + pkin(8) * t418;
t256 = pkin(9) * t319 + t280;
t196 = t256 * t458 + t257 * t322;
t195 = -t256 * t322 + t257 * t458;
t277 = t319 * t322 + t323 * t391;
t166 = -pkin(3) * t418 - qJ(4) * t277 + t195;
t276 = -t319 * t458 + t322 * t419;
t175 = -qJ(4) * t276 + t196;
t114 = t166 * t439 - t175 * t317;
t207 = -t276 * t317 + t277 * t439;
t90 = -pkin(4) * t418 - pkin(10) * t207 + t114;
t115 = t166 * t317 + t175 * t439;
t206 = t276 * t439 + t277 * t317;
t94 = -pkin(10) * t206 + t115;
t449 = t321 * t90 + t457 * t94;
t385 = qJD(2) * t418;
t216 = qJD(3) * t277 + t322 * t385;
t379 = t458 * qJD(2);
t217 = -qJD(3) * t276 + t379 * t418;
t174 = -t216 * t317 + t217 * t439;
t143 = -qJD(3) * t196 + t267 * t458 - t269 * t322;
t112 = pkin(3) * t386 - qJ(4) * t217 - qJD(4) * t277 + t143;
t142 = -t256 * t407 + t257 * t381 + t267 * t322 + t269 * t458;
t116 = -qJ(4) * t216 - qJD(4) * t276 + t142;
t67 = t112 * t439 - t116 * t317;
t54 = pkin(4) * t386 - pkin(10) * t174 + t67;
t173 = t216 * t439 + t217 * t317;
t68 = t112 * t317 + t116 * t439;
t56 = -pkin(10) * t173 + t68;
t14 = -qJD(5) * t449 - t321 * t56 + t457 * t54;
t326 = qJD(1) ^ 2;
t456 = pkin(1) * t325;
t1 = t2 * t324;
t212 = t282 * t321 + t330;
t213 = t282 * t457 - t334;
t312 = -pkin(3) * t458 - pkin(2);
t247 = pkin(4) * t336 + t312;
t141 = pkin(5) * t212 - pkin(11) * t213 + t247;
t91 = t141 * t324 - t145 * t320;
t454 = qJD(6) * t91 + t320 * t474 - t324 * t475;
t92 = t141 * t320 + t145 * t324;
t453 = -qJD(6) * t92 + t320 * t475 + t324 * t474;
t452 = pkin(5) * t388 - t450;
t146 = t206 * t457 + t207 * t321;
t446 = t146 * t73;
t445 = t212 * t73;
t437 = t122 * t120;
t428 = t341 * t298;
t427 = t192 * t298;
t426 = t341 ^ 2;
t425 = t213 * t320;
t424 = t213 * t324;
t423 = t240 * t238;
t422 = t240 * t298;
t421 = t298 * t322;
t420 = t314 * t326;
t413 = t323 ^ 2 - t325 ^ 2;
t396 = t323 * t420;
t395 = t320 * t418;
t394 = t133 * t458;
t393 = t211 * t458;
t392 = t298 * t458;
t389 = t458 * t225;
t383 = t314 * t409;
t159 = pkin(3) * t240 + pkin(4) * t341;
t66 = t159 + t85;
t17 = -t320 * t37 + t324 * t66;
t373 = -t253 * t320 - t17;
t18 = t320 * t66 + t324 * t37;
t372 = t253 * t324 - t18;
t369 = t320 * t416 - t324 * t388;
t368 = t320 * t388 + t324 * t416;
t364 = qJD(2) + 0.2e1 * t411;
t363 = t325 * t396;
t360 = t298 * t388;
t357 = pkin(1) * t473;
t263 = pkin(11) + t272;
t355 = -t263 * t73 - t496;
t50 = -pkin(11) * t418 + t449;
t147 = -t206 * t321 + t207 * t457;
t255 = t306 + (-pkin(2) - t456) * t319;
t209 = pkin(3) * t276 + t255;
t158 = pkin(4) * t206 + t209;
t74 = pkin(5) * t146 - pkin(11) * t147 + t158;
t24 = t320 * t74 + t324 * t50;
t23 = -t320 * t50 + t324 * t74;
t351 = t314 * t323 * t377;
t349 = t131 * t502 - t470;
t51 = -t321 * t94 + t457 * t90;
t136 = t147 * t320 + t324 * t418;
t13 = t321 * t54 + t380 * t90 - t406 * t94 + t457 * t56;
t343 = t192 ^ 2;
t340 = t213 * t403 - t369;
t339 = -t213 * t404 - t368;
t337 = t346 * t318;
t270 = t280 * qJD(2);
t331 = t352 - t381;
t194 = pkin(3) * t216 + t270;
t328 = t1 - t3 * t320 + (-t16 * t320 + t324 * t354) * qJD(6);
t125 = pkin(4) * t173 + t194;
t262 = -pkin(5) - t271;
t215 = (-t298 * t318 - t383) * t408;
t144 = -t204 * t457 + t205 * t321;
t137 = t147 * t324 - t395;
t81 = qJD(5) * t147 + t173 * t457 + t321 * t174;
t80 = t173 * t321 - t174 * t457 + t206 * t380 + t207 * t406;
t60 = -qJD(6) * t395 + t147 * t403 - t320 * t80 - t324 * t386;
t59 = qJD(6) * t136 - t320 * t386 + t324 * t80;
t49 = pkin(5) * t418 - t51;
t27 = pkin(5) * t81 + pkin(11) * t80 + t125;
t20 = t320 * t85 + t324 * t34;
t19 = -t320 * t34 + t324 * t85;
t12 = -pkin(5) * t386 - t14;
t11 = pkin(11) * t386 + t13;
t5 = -qJD(6) * t24 - t11 * t320 + t27 * t324;
t4 = qJD(6) * t23 + t11 * t324 + t27 * t320;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t351, t413 * t473, t364 * t385, -0.2e1 * t351, -t364 * t386, 0, -t250 * t319 - t270 * t365 + t323 * t357, -t249 * t319 - t269 * t365 + t325 * t357 (t249 * t325 + t250 * t323 + (-t265 * t325 - t268 * t323) * qJD(2) + (t269 * t325 + t270 * t323 + (-t279 * t325 - t280 * t323) * qJD(2)) * qJD(1)) * t318, t249 * t280 - t250 * t279 - t265 * t270 + t268 * t269, -t210 * t277 + t217 * t240, t210 * t276 - t211 * t277 - t216 * t240 - t217 * t238, -t217 * t298 + (t210 * t325 + (qJD(1) * t277 + t240) * t408) * t318, t211 * t276 + t216 * t238, t216 * t298 + (t211 * t325 + (-qJD(1) * t276 - t238) * t408) * t318, t215, -t143 * t298 + t211 * t255 + t216 * t225 + t238 * t270 + t250 * t276 + (-t134 * t325 + (qJD(1) * t195 + t184) * t408) * t318, t142 * t298 - t210 * t255 + t217 * t225 + t240 * t270 + t250 * t277 + (t133 * t325 + (-qJD(1) * t196 - t185) * t408) * t318, -t133 * t276 - t134 * t277 - t142 * t238 - t143 * t240 - t184 * t217 - t185 * t216 + t195 * t210 - t196 * t211, t133 * t196 + t134 * t195 + t142 * t185 + t143 * t184 + t225 * t270 + t250 * t255, t157 * t207 + t174 * t341, -t156 * t207 - t157 * t206 - t173 * t341 + t174 * t192, -t174 * t298 + (-t157 * t325 + (qJD(1) * t207 + t341) * t408) * t318, t156 * t206 - t173 * t192, t173 * t298 + (t156 * t325 + (-qJD(1) * t206 + t192) * t408) * t318, t215, -t67 * t298 - t194 * t192 + t209 * t156 + t181 * t206 + t190 * t173 + (-t57 * t325 + (qJD(1) * t114 + t98) * t408) * t318, t157 * t209 + t174 * t190 + t181 * t207 + t341 * t194 + t298 * t68 + (t325 * t58 + (-qJD(1) * t115 - t99) * t408) * t318, -t114 * t157 - t115 * t156 - t173 * t99 - t174 * t98 + t192 * t68 - t206 * t58 - t207 * t57 - t341 * t67, t114 * t57 + t115 * t58 + t181 * t209 + t190 * t194 + t67 * t98 + t68 * t99, -t147 * t72 - t476 * t80, -t131 * t80 + t146 * t72 - t147 * t73 - t476 * t81, -t80 * t400 + (t476 * t408 + t72 * t325 + (t147 * t408 + t325 * t80) * qJD(1)) * t318, -t131 * t81 + t446, -t81 * t400 + (t131 * t408 + t73 * t325 + (-t146 * t408 + t325 * t81) * qJD(1)) * t318 (-t337 - t383) * t408, t14 * t400 - t125 * t131 + t158 * t73 + t117 * t146 + t135 * t81 + (t34 * t408 - t10 * t325 + (-t14 * t325 + t408 * t51) * qJD(1)) * t318, -t13 * t400 + t125 * t476 - t158 * t72 + t117 * t147 - t135 * t80 + (-t35 * t408 - t374 * t325 + (t13 * t325 - t408 * t449) * qJD(1)) * t318, -t10 * t147 + t13 * t131 - t14 * t476 + t146 * t374 + t34 * t80 - t35 * t81 - t449 * t73 + t51 * t72, t10 * t51 + t117 * t158 + t125 * t135 + t13 * t35 + t14 * t34 - t374 * t449, -t122 * t59 - t137 * t47, t120 * t59 - t122 * t60 + t136 * t47 - t137 * t48, t122 * t81 + t137 * t73 - t146 * t47 - t492 * t59, t120 * t60 + t136 * t48, -t120 * t81 - t136 * t73 - t146 * t48 - t492 * t60, t492 * t81 + t446, t12 * t120 + t136 * t8 + t146 * t3 + t23 * t73 + t32 * t60 - t354 * t81 + t48 * t49 + t492 * t5, t12 * t122 + t137 * t8 - t146 * t2 - t16 * t81 - t24 * t73 - t32 * t59 - t4 * t492 - t47 * t49, -t120 * t4 - t122 * t5 - t136 * t2 - t137 * t3 - t16 * t60 + t23 * t47 - t24 * t48 - t354 * t59, t12 * t32 + t16 * t4 + t2 * t24 + t23 * t3 - t354 * t5 + t49 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t363, t413 * t420, -t318 * t326 * t417, t363, t365 * t388 - t302, 0, pkin(1) * t396 + t268 * t365 - t250, pkin(8) * t302 + t265 * t365 + (-t319 * t401 + t420) * t456, 0, 0, -t210 * t322 - t240 * t331, -t210 * t458 + t331 * t238 + (-t211 + t422) * t322, -t298 * t381 + (t325 * t392 + (qJD(2) * t322 - t240) * t323) * t412, -t238 * t421 - t393, t298 * t407 + (-t325 * t421 + (t379 + t238) * t323) * t412, t360, -t250 * t458 - pkin(2) * t211 + t200 * t298 - t268 * t238 + (pkin(9) * t392 + t225 * t322) * qJD(3) + (-t184 * t323 + (-pkin(9) * t408 - t225 * t325) * t322) * t412, pkin(2) * t210 - t201 * t298 - t268 * t240 + t250 * t322 + (-pkin(9) * t421 + t389) * qJD(3) + (-t325 * t389 + (-pkin(9) * t379 + t185) * t323) * t412, t394 + t200 * t240 + t201 * t238 + t331 * t184 + (t240 * t381 - t393) * pkin(9) + ((-t210 + t477) * pkin(9) + t465) * t322, -t250 * pkin(2) - t184 * t200 - t185 * t201 - t225 * t268 + (t394 - t134 * t322 + (-t184 * t458 - t185 * t322) * qJD(3)) * pkin(9), t157 * t282 + t341 * t493, -t282 * t156 - t157 * t336 + t192 * t493 + t341 * t461, t282 * t302 - t298 * t493 - t341 * t388, t156 * t336 + t192 * t461, -t192 * t388 - t298 * t461 - t329 * t408, t360, t312 * t156 + t181 * t336 - t190 * t461 - t192 * t462 + t223 * t302 + t298 * t464 - t388 * t98, t312 * t157 + t181 * t282 + t190 * t493 - t224 * t302 + t298 * t463 + t341 * t462 + t388 * t99, -t224 * t156 - t223 * t157 - t98 * t221 - t57 * t282 + t461 * t99 + t463 * t192 + (qJD(3) * t98 - t58) * t336 + t464 * t341, t181 * t312 + t190 * t462 + t223 * t57 + t224 * t58 + t463 * t99 - t464 * t98, -t213 * t72 - t416 * t476, -t131 * t416 + t212 * t72 - t213 * t73 - t415 * t476 (t416 * t325 + (qJD(2) * t213 - t476) * t323) * t412 - t416 * t400, -t131 * t415 + t445 (t415 * t325 + (-qJD(2) * t212 - t131) * t323) * t412 - t415 * t400, t337 * t410, t247 * t73 + t117 * t212 + t415 * t135 - t414 * t131 + (-t450 * t325 + (-qJD(2) * t144 - t34) * t323) * t412 + t450 * t400, -t247 * t72 + t117 * t213 - t416 * t135 + t414 * t476 + (-t451 * t325 + (-qJD(2) * t145 + t35) * t323) * t412 + t451 * t400, -t10 * t213 - t131 * t451 - t144 * t72 - t145 * t73 + t212 * t374 + t34 * t416 - t35 * t415 - t450 * t476, -t10 * t144 + t117 * t247 + t135 * t414 - t145 * t374 + t34 * t450 - t35 * t451, t122 * t339 - t424 * t47, t369 * t122 + t368 * t120 + (t44 - t46 + (t120 * t320 - t122 * t324) * qJD(6)) * t213, t122 * t415 - t212 * t47 + t339 * t492 + t424 * t73, t120 * t340 + t425 * t48, -t120 * t415 - t212 * t48 - t340 * t492 - t425 * t73, t415 * t492 + t445, t120 * t452 + t144 * t48 + t212 * t3 + t32 * t340 - t354 * t415 + t425 * t8 + t453 * t492 + t73 * t91, t122 * t452 - t144 * t47 - t16 * t415 - t2 * t212 + t32 * t339 + t424 * t8 - t454 * t492 - t73 * t92, t47 * t91 - t48 * t92 + t369 * t16 - t368 * t354 - t453 * t122 - t454 * t120 + (-t2 * t320 - t3 * t324 + (-t16 * t324 - t320 * t354) * qJD(6)) * t213, t144 * t8 + t16 * t454 + t2 * t92 + t3 * t91 + t32 * t452 - t354 * t453; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t423, -t238 ^ 2 + t240 ^ 2, -t238 * t298 - t210, -t423, -t211 - t422, t302, -t225 * t240 - t465, -t184 * t298 + t225 * t238 - t133, 0, 0, -t478, -t343 + t426, t157 + t427, t478, -t156 - t428, t302, t103 * t298 - t190 * t341 + (t192 * t240 + t302 * t439) * pkin(3) + t57, -t104 * t298 - t190 * t192 + (-t240 * t341 - t302 * t317) * pkin(3) - t58 (-t156 * t317 - t157 * t439) * pkin(3) + (t99 + t103) * t341 + (-t104 + t98) * t192, -t98 * t103 - t99 * t104 + (-t190 * t240 + t317 * t58 + t439 * t57) * pkin(3), -t430, t501, t495, t430, t487, t302, -qJD(5) * t490 + t131 * t159 + t271 * t302 + t440 * t298 + t481, -t159 * t476 - t272 * t302 + t346 * t441 + t491, t271 * t72 - t272 * t73 + t490 * t476 + (t34 + t441) * t131, t10 * t271 - t135 * t159 - t272 * t374 - t34 * t440 + t35 * t441, t498, t507, t497, t503, t349 + t438, -t484, t262 * t48 + t355 * t320 + t440 * t120 + (-t263 * t403 + t373) * t492 + t480, -t262 * t47 + t355 * t324 + t440 * t122 + (t263 * t404 - t372) * t492 + t482, t120 * t18 + t122 * t17 + t1 + (-t120 * t253 - t131 * t354 - t263 * t48 + (t122 * t263 + t354) * qJD(6)) * t324 + (t122 * t253 + t131 * t16 - t263 * t47 - t3 + (t120 * t263 - t16) * qJD(6)) * t320, t16 * t372 + t262 * t8 + t263 * t328 + t32 * t440 - t354 * t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156 - t428, t157 - t427, -t343 - t426, -t192 * t99 + t341 * t98 + t181, 0, 0, 0, 0, 0, 0, -t344 + t370 + 0.2e1 * t483, -t72 - t468, -t433 - t434, -t131 * t35 + t34 * t476 + t117, 0, 0, 0, 0, 0, 0, t349 - t438, -t324 * t492 ^ 2 - t436 - t69 (t120 * t131 + t47) * t324 + t469 + t443, -t476 * t32 + t320 * t506 + t324 * t500; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t430, t501, t495, t430, t487, t302, -t298 * t35 + t481, -t34 * t346 + t491, 0, 0, t498, t507, t497, t503, -t492 * t502 + t438 + t71, -t484, -pkin(5) * t48 - pkin(11) * t442 - t120 * t35 - t19 * t492 - t320 * t496 + t480, pkin(5) * t47 + pkin(11) * t470 - t122 * t35 + t20 * t492 - t32 * t494 + t482, t120 * t20 + t122 * t19 + t1 + (t448 + (-t48 + t405) * pkin(11)) * t324 + ((qJD(6) * t120 - t47) * pkin(11) - t500) * t320, -pkin(5) * t8 + pkin(11) * t328 - t16 * t20 + t19 * t354 - t32 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t437, -t120 ^ 2 + t122 ^ 2, t120 * t492 - t47, -t437, t122 * t492 - t48, t73, -t122 * t32 + t500, t120 * t32 - t506, 0, 0;];
tauc_reg  = t6;
