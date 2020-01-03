% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:50
% DurationCPUTime: 15.75s
% Computational Cost: add. (13737->691), mult. (18251->967), div. (0->0), fcn. (16979->8), ass. (0->340)
t248 = qJ(1) + pkin(8);
t245 = sin(t248);
t246 = cos(t248);
t253 = cos(qJ(4));
t252 = cos(qJ(5));
t249 = sin(qJ(5));
t250 = sin(qJ(4));
t393 = t249 * t250;
t159 = t245 * t252 + t246 * t393;
t392 = t250 * t252;
t160 = -t245 * t249 + t246 * t392;
t154 = Icges(6,4) * t160;
t394 = t246 * t253;
t88 = Icges(6,2) * t159 + Icges(6,6) * t394 - t154;
t153 = Icges(6,4) * t159;
t90 = Icges(6,1) * t160 - Icges(6,5) * t394 - t153;
t313 = -t159 * t88 - t160 * t90;
t157 = -t245 * t393 + t246 * t252;
t158 = t245 * t392 + t246 * t249;
t396 = t245 * t253;
t405 = Icges(6,4) * t158;
t86 = Icges(6,2) * t157 - Icges(6,6) * t396 + t405;
t152 = Icges(6,4) * t157;
t89 = Icges(6,1) * t158 - Icges(6,5) * t396 + t152;
t427 = t157 * t86 + t158 * t89;
t83 = Icges(6,5) * t158 + Icges(6,6) * t157 - Icges(6,3) * t396;
t85 = -Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t394;
t476 = t313 + t427 + (-t245 * t83 - t246 * t85) * t253;
t364 = qJD(5) * t253;
t368 = qJD(4) * t246;
t182 = -t245 * t364 + t368;
t365 = qJD(5) * t250;
t233 = qJD(1) + t365;
t27 = t159 * t86 - t160 * t89 + t83 * t394;
t297 = Icges(6,5) * t252 - Icges(6,6) * t249;
t162 = Icges(6,3) * t250 + t253 * t297;
t403 = Icges(6,4) * t252;
t299 = -Icges(6,2) * t249 + t403;
t164 = Icges(6,6) * t250 + t253 * t299;
t404 = Icges(6,4) * t249;
t301 = Icges(6,1) * t252 - t404;
t166 = Icges(6,5) * t250 + t253 * t301;
t53 = t159 * t164 - t160 * t166 + t162 * t394;
t475 = -t182 * t27 - t233 * t53;
t329 = -rSges(4,2) * t246 + t245 * rSges(4,3);
t254 = cos(qJ(1));
t247 = t254 * pkin(1);
t461 = t246 * pkin(2) + t245 * qJ(3);
t354 = t247 + t461;
t470 = t329 + t354;
t303 = t249 * t88 + t252 * t90;
t34 = t250 * t85 - t253 * t303;
t26 = t157 * t88 - t158 * t90 - t396 * t85;
t255 = qJD(1) ^ 2;
t366 = qJD(4) * t253;
t371 = qJD(1) * t250;
t469 = t245 * t366 + t246 * t371;
t406 = Icges(5,4) * t253;
t211 = -Icges(5,2) * t250 + t406;
t302 = Icges(5,1) * t250 + t406;
t380 = t211 + t302;
t407 = Icges(5,4) * t250;
t213 = Icges(5,1) * t253 - t407;
t300 = Icges(5,2) * t253 + t407;
t381 = -t300 + t213;
t468 = (t250 * t380 - t253 * t381) * qJD(1);
t316 = rSges(6,1) * t252 - rSges(6,2) * t249;
t176 = rSges(6,3) * t250 + t253 * t316;
t369 = qJD(4) * t245;
t181 = t246 * t364 + t369;
t432 = t253 * pkin(4);
t222 = pkin(7) * t250 + t432;
t317 = rSges(6,1) * t160 - rSges(6,2) * t159;
t94 = rSges(6,3) * t394 - t317;
t467 = t176 * t181 + t222 * t369 - t233 * t94;
t465 = 0.2e1 * qJD(4);
t52 = t157 * t164 + t158 * t166 - t162 * t396;
t464 = t181 * t26 + t52 * t233;
t161 = Icges(6,3) * t253 - t250 * t297;
t292 = t164 * t249 - t166 * t252;
t304 = t249 * t86 - t252 * t89;
t256 = t181 * (-t162 * t246 + t303) + t182 * (t162 * t245 + t304) + t233 * (t161 + t292);
t463 = t256 * t253;
t143 = -Icges(5,6) * t245 + t246 * t300;
t145 = -Icges(5,5) * t245 + t246 * t302;
t294 = t143 * t253 + t145 * t250;
t462 = t294 * t246;
t235 = t246 * qJ(3);
t435 = pkin(2) * t245;
t191 = -t235 + t435;
t251 = sin(qJ(1));
t436 = pkin(1) * t251;
t323 = -pkin(6) * t245 - t436;
t287 = -t191 + t323;
t330 = t246 * rSges(3,1) - rSges(3,2) * t245;
t460 = t247 + t330;
t172 = t213 * t246;
t402 = Icges(5,5) * t246;
t101 = -qJD(4) * t172 + (t245 * t302 + t402) * qJD(1);
t298 = Icges(5,5) * t250 + Icges(5,6) * t253;
t141 = -Icges(5,3) * t245 + t246 * t298;
t375 = qJD(1) * t141;
t73 = t143 * t250 - t145 * t253;
t142 = Icges(5,6) * t246 + t245 * t300;
t170 = t211 * t246;
t99 = qJD(1) * t142 - qJD(4) * t170;
t458 = qJD(4) * t73 + t101 * t250 + t253 * t99 + t375;
t456 = t233 * t249 - t252 * t366;
t455 = t233 * t252 + t249 * t366;
t196 = t300 * qJD(4);
t197 = t302 * qJD(4);
t209 = Icges(5,5) * t253 - Icges(5,6) * t250;
t454 = qJD(1) * t209 + qJD(4) * (t211 * t250 - t213 * t253) + t196 * t253 + t197 * t250;
t100 = qJD(1) * t143 + t211 * t369;
t171 = t213 * t245;
t102 = qJD(1) * t145 + qJD(4) * t171;
t218 = Icges(5,4) * t396;
t397 = t245 * t250;
t144 = Icges(5,1) * t397 + t218 + t402;
t295 = t142 * t250 - t144 * t253;
t140 = Icges(5,3) * t246 + t245 * t298;
t376 = qJD(1) * t140;
t453 = qJD(4) * t295 - t100 * t253 - t102 * t250 + t376;
t387 = t145 + t170;
t389 = t143 - t172;
t451 = t250 * t389 - t253 * t387;
t388 = -Icges(5,2) * t397 + t144 + t218;
t390 = t142 - t171;
t450 = t250 * t390 - t253 * t388;
t188 = (-Icges(6,2) * t252 - t404) * t253;
t262 = t181 * (Icges(6,2) * t160 + t153 - t90) + t182 * (-Icges(6,2) * t158 + t152 + t89) + t233 * (t166 + t188);
t189 = (-Icges(6,1) * t249 - t403) * t253;
t261 = t181 * (Icges(6,1) * t159 + t154 - t88) + t182 * (Icges(6,1) * t157 - t405 - t86) - t233 * (t164 - t189);
t449 = -pkin(2) - pkin(6);
t361 = qJD(4) * qJD(5);
t343 = t250 * t361;
t138 = -qJD(1) * t181 + t245 * t343;
t448 = t138 / 0.2e1;
t139 = qJD(1) * t182 - t246 * t343;
t447 = t139 / 0.2e1;
t446 = -t181 / 0.2e1;
t445 = t181 / 0.2e1;
t444 = -t182 / 0.2e1;
t443 = t182 / 0.2e1;
t442 = -t233 / 0.2e1;
t441 = t233 / 0.2e1;
t440 = t245 / 0.2e1;
t439 = -t246 / 0.2e1;
t437 = rSges(6,3) + pkin(7);
t434 = pkin(4) * t250;
t433 = pkin(7) * t253;
t243 = t246 * pkin(6);
t367 = qJD(4) * t250;
t347 = t245 * t367;
t370 = qJD(1) * t253;
t350 = t246 * t370;
t274 = t347 - t350;
t328 = qJD(5) + t371;
t289 = t328 * t249;
t80 = -t245 * t455 - t246 * t289;
t288 = t328 * t252;
t81 = -t245 * t456 + t246 * t288;
t42 = Icges(6,5) * t81 + Icges(6,6) * t80 + Icges(6,3) * t274;
t44 = Icges(6,4) * t81 + Icges(6,2) * t80 + Icges(6,6) * t274;
t46 = Icges(6,1) * t81 + Icges(6,4) * t80 + Icges(6,5) * t274;
t7 = (qJD(4) * t304 + t42) * t250 + (qJD(4) * t83 - t249 * t44 + t252 * t46 + (-t249 * t89 - t252 * t86) * qJD(5)) * t253;
t431 = t7 * t182;
t352 = t245 * t370;
t273 = -t246 * t367 - t352;
t78 = -t245 * t289 + t246 * t455;
t79 = t245 * t288 + t246 * t456;
t41 = Icges(6,5) * t79 + Icges(6,6) * t78 + Icges(6,3) * t273;
t43 = Icges(6,4) * t79 + Icges(6,2) * t78 + Icges(6,6) * t273;
t45 = Icges(6,1) * t79 + Icges(6,4) * t78 + Icges(6,5) * t273;
t8 = (qJD(4) * t303 + t41) * t250 + (qJD(4) * t85 - t249 * t43 + t252 * t45 + (t249 * t90 - t252 * t88) * qJD(5)) * t253;
t430 = t8 * t181;
t429 = -qJD(1) / 0.2e1;
t187 = (-Icges(6,5) * t249 - Icges(6,6) * t252) * t253;
t119 = qJD(4) * t161 + qJD(5) * t187;
t163 = Icges(6,6) * t253 - t250 * t299;
t120 = qJD(4) * t163 + qJD(5) * t188;
t165 = Icges(6,5) * t253 - t250 * t301;
t121 = qJD(4) * t165 + qJD(5) * t189;
t24 = (qJD(4) * t292 + t119) * t250 + (qJD(4) * t162 - t120 * t249 + t121 * t252 + (-t164 * t252 - t166 * t249) * qJD(5)) * t253;
t342 = t253 * t361;
t63 = t162 * t250 - t253 * t292;
t428 = t24 * t233 + t63 * t342;
t424 = rSges(4,3) * t246;
t422 = rSges(6,3) * t253;
t216 = rSges(5,1) * t253 - rSges(5,2) * t250;
t174 = t216 * t246;
t241 = t246 * rSges(5,3);
t318 = rSges(5,1) * t250 + rSges(5,2) * t253;
t105 = -qJD(4) * t174 + (t245 * t318 + t241) * qJD(1);
t232 = qJD(3) * t246;
t149 = qJD(1) * t461 - t232;
t198 = t318 * qJD(4);
t362 = qJD(1) * qJD(3);
t226 = t246 * t362;
t322 = -t247 - t243;
t271 = t255 * t322 + t226;
t349 = t216 * t368;
t48 = -t198 * t369 + (-t105 - t149 + t349) * qJD(1) + t271;
t419 = t245 * t48;
t356 = rSges(5,1) * t469 + rSges(5,2) * t350;
t106 = (-rSges(5,2) * t367 - rSges(5,3) * qJD(1)) * t245 + t356;
t183 = t216 * t369;
t373 = qJD(1) * t245;
t231 = qJD(3) * t245;
t372 = qJD(1) * t246;
t379 = qJ(3) * t372 + t231;
t386 = qJD(1) * (-pkin(2) * t373 + t379) + t245 * t362;
t263 = t255 * t323 + t386;
t47 = t198 * t368 + (t106 + t183) * qJD(1) + t263;
t417 = t246 * t47;
t239 = t245 * rSges(5,3);
t147 = t318 * t246 - t239;
t60 = t183 + t231 + (t147 + t287) * qJD(1);
t416 = t246 * t60;
t33 = t250 * t83 - t253 * t304;
t415 = t33 * t138;
t414 = t34 * t139;
t223 = pkin(4) * t397;
t177 = -pkin(7) * t396 + t223;
t385 = t158 * rSges(6,1) + t157 * rSges(6,2);
t92 = -rSges(6,3) * t396 + t385;
t409 = -t177 - t92;
t224 = pkin(7) * t394;
t395 = t246 * t250;
t179 = pkin(4) * t395 - t224;
t408 = t179 - t94;
t398 = t209 * t246;
t167 = t245 * t209;
t175 = -t250 * t316 + t422;
t190 = (-rSges(6,1) * t249 - rSges(6,2) * t252) * t253;
t122 = qJD(4) * t175 + qJD(5) * t190;
t221 = t433 - t434;
t200 = qJD(4) * t221;
t391 = t122 + t200;
t382 = t176 + t222;
t378 = rSges(4,2) * t373 + rSges(4,3) * t372;
t377 = -qJD(1) * t191 + t231;
t374 = qJD(1) * t298;
t291 = t211 * t253 + t213 * t250;
t104 = t246 * t291 - t167;
t363 = t104 * qJD(1);
t360 = -rSges(5,3) + t449;
t359 = t255 * t436;
t358 = t255 * t247;
t357 = t81 * rSges(6,1) + t80 * rSges(6,2) + rSges(6,3) * t347;
t54 = t246 * t140 + t142 * t396 + t144 * t397;
t55 = -t246 * t141 - t143 * t396 - t145 * t397;
t355 = pkin(4) * t469 + pkin(7) * t347;
t146 = rSges(5,1) * t397 + rSges(5,2) * t396 + t241;
t353 = t437 * t253;
t348 = t222 * t368;
t339 = -t370 / 0.2e1;
t338 = -t369 / 0.2e1;
t337 = t369 / 0.2e1;
t336 = -t368 / 0.2e1;
t334 = t366 / 0.2e1;
t332 = t235 - t436;
t327 = t243 + t354;
t325 = qJD(5) * t334;
t324 = rSges(6,1) * t79 + rSges(6,2) * t78;
t321 = -t422 + t434;
t193 = rSges(3,1) * t245 + rSges(3,2) * t246;
t25 = -t396 * t83 + t427;
t312 = t245 * t26 + t246 * t25;
t311 = t245 * t25 - t246 * t26;
t28 = t394 * t85 - t313;
t310 = t245 * t28 + t246 * t27;
t309 = t245 * t27 - t246 * t28;
t308 = t245 * t34 + t246 * t33;
t307 = t245 * t33 - t246 * t34;
t286 = t461 - t322;
t61 = -t349 - t232 + (t146 + t286) * qJD(1);
t306 = t245 * t60 - t246 * t61;
t305 = t245 * t94 + t246 * t92;
t296 = t142 * t253 + t144 * t250;
t293 = -t146 * t245 - t147 * t246;
t283 = -t253 * t41 + t367 * t85;
t282 = -t253 * t42 + t367 * t83;
t277 = (t245 * t55 + t246 * t54) * qJD(4);
t135 = t245 * t140;
t56 = -t296 * t246 + t135;
t57 = -t141 * t245 + t462;
t276 = (t245 * t57 + t246 * t56) * qJD(4);
t272 = t162 * t233 + t181 * t85 + t182 * t83;
t269 = (Icges(6,5) * t157 - Icges(6,6) * t158) * t182 + (Icges(6,5) * t159 + Icges(6,6) * t160) * t181 + t187 * t233;
t266 = -qJD(1) * t294 - qJD(4) * t398 + t376;
t265 = qJD(1) * t296 + qJD(4) * t167 + t375;
t264 = t291 * qJD(1) - t298 * qJD(4);
t260 = t105 * t246 - t106 * t245 + (-t146 * t246 + t147 * t245) * qJD(1);
t29 = -t181 * t92 + t182 * t94 + qJD(2) + (-t177 * t245 - t179 * t246) * qJD(4);
t37 = t231 + (t179 + t287) * qJD(1) + t467;
t38 = -t348 - t176 * t182 + t233 * t92 - t232 + (t177 + t286) * qJD(1);
t257 = t29 * t305 + (-t245 * t38 - t246 * t37) * t176;
t180 = t222 * t246;
t178 = t222 * t245;
t173 = t216 * t245;
t134 = t176 * t246;
t133 = t176 * t245;
t132 = t166 * t246;
t131 = t166 * t245;
t130 = t164 * t246;
t129 = t164 * t245;
t118 = -pkin(7) * t350 + t355;
t117 = t273 * pkin(7) + (t245 * t371 - t246 * t366) * pkin(4);
t116 = rSges(6,1) * t159 + rSges(6,2) * t160;
t115 = rSges(6,1) * t157 - rSges(6,2) * t158;
t103 = t245 * t291 + t398;
t82 = t103 * qJD(1);
t75 = -t358 + t226 + (-qJD(1) * t329 - t149) * qJD(1);
t74 = qJD(1) * t378 - t359 + t386;
t70 = qJD(4) * t293 + qJD(2);
t50 = -rSges(6,3) * t350 + t357;
t49 = rSges(6,3) * t273 + t324;
t40 = -t245 * t454 + t264 * t246;
t39 = t264 * t245 + t246 * t454;
t36 = t294 * qJD(4) + t101 * t253 - t250 * t99;
t35 = -qJD(4) * t296 - t100 * t250 + t102 * t253;
t30 = t260 * qJD(4);
t23 = t276 - t363;
t22 = t82 + t277;
t16 = -t119 * t396 + t120 * t157 + t121 * t158 + t162 * t274 + t164 * t80 + t166 * t81;
t15 = t119 * t394 + t120 * t159 - t121 * t160 + t162 * t273 + t164 * t78 + t166 * t79;
t14 = t122 * t181 + t139 * t176 - t233 * t49 + (t200 * t245 - t364 * t94) * qJD(4) + (-t117 - t149 + t348) * qJD(1) + t271;
t13 = qJD(1) * t118 - t122 * t182 - t138 * t176 + t233 * t50 + (-t200 * t246 + t222 * t373 + t364 * t92) * qJD(4) + t263;
t12 = t181 * t34 + t182 * t33 + t233 * t63;
t11 = t181 * t28 - t475;
t10 = t182 * t25 + t464;
t9 = t138 * t94 - t139 * t92 - t181 * t50 + t182 * t49 + (t117 * t246 - t118 * t245 + (-t177 * t246 + t179 * t245) * qJD(1)) * qJD(4);
t6 = t157 * t43 + t158 * t45 + t245 * t283 - t350 * t85 + t80 * t88 - t81 * t90;
t5 = t157 * t44 + t158 * t46 + t245 * t282 - t350 * t83 + t80 * t86 + t81 * t89;
t4 = t159 * t43 - t160 * t45 - t246 * t283 - t352 * t85 + t78 * t88 - t79 * t90;
t3 = t159 * t44 - t160 * t46 - t246 * t282 - t352 * t83 + t78 * t86 + t79 * t89;
t2 = t138 * t25 + t139 * t26 + t16 * t233 + t181 * t6 + t182 * t5 + t342 * t52;
t1 = t138 * t27 + t139 * t28 + t15 * t233 + t181 * t4 + t182 * t3 + t342 * t53;
t17 = [m(3) * ((-t193 * t255 - t359) * t460 + (-t358 + (-0.2e1 * t330 - t247 + t460) * t255) * (-t193 - t436)) + (-qJD(4) * t291 + t196 * t250 - t197 * t253) * qJD(1) + (t82 + ((-t56 + t135 + t55) * t245 + (t57 - t462 + (t141 - t296) * t245 + t54) * t246) * qJD(4)) * t338 + t415 / 0.2e1 + t16 * t443 + t53 * t447 + t52 * t448 + t414 / 0.2e1 + t428 + t430 / 0.2e1 + ((t28 + t476) * t182 + t464) * t446 + t431 / 0.2e1 + (t15 + t10) * t445 + (t11 + (-t25 + t476) * t181 + t475) * t444 + (t363 + (t141 * t245 ^ 2 + (-t135 + t55 + (t141 + t296) * t246) * t246) * qJD(4) + t23) * t336 + (-(-t37 + t377 + (t179 + t323) * qJD(1) + t467) * t38 + t14 * (-t224 + t317 + t332) + t37 * (t232 - t324) + t13 * (t223 + t327 + t385) + t38 * (t355 + t357 + t379) + (-t13 * t353 + t14 * t449) * t245 + (t14 * t321 + t37 * (t250 * t437 + t432) * qJD(4)) * t246 + ((-t251 * t38 - t254 * t37) * pkin(1) + (-t353 * t38 + t37 * t449) * t246 + (t37 * (-qJ(3) - t321 + t433) + t38 * t449) * t245) * qJD(1)) * m(6) + (t48 * (-t239 + t287) + t60 * t232 + t47 * (t327 + t146) + t61 * (-rSges(5,2) * t347 + t356 + t379) + (qJD(4) * t216 * t60 + t318 * t48) * t246 + ((-t251 * t61 - t254 * t60) * pkin(1) + t360 * t416 + (t60 * (-qJ(3) - t318) + t61 * t360) * t245) * qJD(1) - (t183 - t60 + (t147 + t323) * qJD(1) + t377) * t61) * m(5) + (t75 * (t424 + (rSges(4,2) - pkin(2)) * t245 + t332) + t74 * t470 + (-t377 + t378 + t379 + (-rSges(4,2) * t245 - t424 - t435) * qJD(1)) * (qJD(1) * t470 - t232)) * m(4) + (t36 + t39 + t22) * t337 + (qJD(1) * t73 + t35 + t40) * t368 / 0.2e1 + (t246 * t104 + (-t295 + t103) * t245) * qJD(4) * t429; m(5) * t30 + m(6) * t9; 0.2e1 * (t13 * t439 + t14 * t440) * m(6) + 0.2e1 * (t419 / 0.2e1 - t417 / 0.2e1) * m(5) + 0.2e1 * (t439 * t74 + t440 * t75) * m(4); -t12 * t364 / 0.2e1 + ((t129 * t157 + t131 * t158) * t182 + (-t130 * t157 - t132 * t158) * t181 + (t157 * t163 + t158 * t165) * t233 + (t253 * t52 - t26 * t395) * qJD(5) + ((qJD(5) * t25 + t272) * t250 - t463) * t245) * t444 + (-qJD(1) * t309 + t245 * t4 + t246 * t3) * t445 + ((t129 * t159 - t131 * t160) * t182 + (-t130 * t159 + t132 * t160) * t181 + (t159 * t163 - t160 * t165) * t233 + (t253 * t53 + t27 * t397) * qJD(5) + ((-qJD(5) * t28 - t272) * t250 + t463) * t246) * t446 + t310 * t447 + t312 * t448 + (-qJD(1) * t307 + t245 * t8 + t246 * t7) * t441 + (((-t129 * t249 + t131 * t252 + t83) * t182 + (t130 * t249 - t132 * t252 + t85) * t181 + (-t163 * t249 + t165 * t252 + t162) * t233 + t63 * qJD(5)) * t253 + (qJD(5) * t307 + t256) * t250) * t442 + (-qJD(1) * t311 + t245 * t6 + t246 * t5) * t443 + ((-t250 * t381 - t253 * t380) * qJD(1) + ((t245 * t389 - t246 * t390) * t253 + (t245 * t387 - t246 * t388) * t250) * qJD(4)) * t429 - t245 * t10 * t365 / 0.2e1 + t308 * t325 + ((-t369 * t398 - t374) * t245 + (t468 + (t450 * t246 + (t167 - t451) * t245) * qJD(4)) * t246) * t338 + ((t167 * t368 - t374) * t246 + (-t468 + (t451 * t245 + (-t398 - t450) * t246) * qJD(4)) * t245) * t336 + qJD(1) * (t245 * t36 + t246 * t35 + (t245 * t295 + t246 * t73) * qJD(1)) / 0.2e1 + (-t37 * (qJD(1) * t180 + t134 * t233 + t175 * t181 + t221 * t369) - t38 * (qJD(1) * t178 + t133 * t233 - t175 * t182 - t221 * t368) - t29 * (-t133 * t181 - t134 * t182 - t178 * t369 - t180 * t368) - ((-t37 * t94 + t38 * t92) * t253 + t257 * t250) * qJD(5) + (-t13 * t382 - t38 * t391 - t9 * t408 + t29 * (t117 + t49) + (t29 * t409 + t37 * t382) * qJD(1)) * t246 + (t14 * t382 + t37 * t391 + t9 * t409 + t29 * (-t118 - t50) + (t29 * t408 + t38 * t382) * qJD(1)) * t245) * m(6) + (t30 * t293 + t70 * t260 - t306 * t198 + (t419 - t417 + (t245 * t61 + t416) * qJD(1)) * t216 - (t173 * t61 + t174 * t60) * qJD(1) - (t70 * (-t173 * t245 - t174 * t246) - t306 * t318) * qJD(4)) * m(5) + (qJD(1) * t39 + t1 + ((t265 * t245 + t246 * t453) * t246 + (t266 * t245 - t246 * t458) * t245 + (-t56 * t245 + t57 * t246) * qJD(1)) * t465) * t440 - (t277 + t22 + t10) * t373 / 0.2e1 + (t276 + t23 + t11) * t372 / 0.2e1 + (qJD(1) * t40 + t11 * t365 + t2 + ((-t245 * t453 + t265 * t246) * t246 + (t245 * t458 + t266 * t246) * t245 + (-t54 * t245 + t55 * t246) * qJD(1)) * t465) * t246 / 0.2e1; -t2 * t396 / 0.2e1 + (t250 * t52 - t253 * t311) * t448 + ((qJD(4) * t311 + t16) * t250 + (-qJD(1) * t312 + qJD(4) * t52 - t245 * t5 + t246 * t6) * t253) * t443 + t1 * t394 / 0.2e1 + (t250 * t53 - t253 * t309) * t447 + ((qJD(4) * t309 + t15) * t250 + (-qJD(1) * t310 + qJD(4) * t53 - t245 * t3 + t246 * t4) * t253) * t445 + t12 * t334 + t250 * (t414 + t415 + t428 + t430 + t431) / 0.2e1 + (t250 * t63 - t253 * t307) * t325 + ((qJD(4) * t307 + t24) * t250 + (-qJD(1) * t308 + qJD(4) * t63 - t245 * t7 + t246 * t8) * t253) * t441 + (t157 * t262 + t158 * t261 - t269 * t396) * t444 + (t262 * t159 - t160 * t261 + t269 * t394) * t446 + (t269 * t250 + (-t249 * t262 + t261 * t252) * t253) * t442 + (t245 * t339 + t250 * t336) * t11 + (t246 * t339 + t250 * t337) * t10 + ((qJD(4) * t257 + t13 * t92 - t14 * t94 - t37 * t49 + t38 * t50) * t250 + (t37 * (-qJD(4) * t94 + t122 * t246) + t38 * (qJD(4) * t92 + t122 * t245) - t9 * t305 + t29 * (-t245 * t49 - t246 * t50 - t372 * t94 + t373 * t92) + (t13 * t245 + t14 * t246 + (-t245 * t37 + t246 * t38) * qJD(1)) * t176) * t253 - t37 * (-t116 * t233 + t181 * t190) - t38 * (t115 * t233 - t182 * t190) - t29 * (-t115 * t181 + t116 * t182)) * m(6);];
tauc = t17(:);
