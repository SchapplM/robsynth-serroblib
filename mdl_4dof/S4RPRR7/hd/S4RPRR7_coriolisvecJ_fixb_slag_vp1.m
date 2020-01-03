% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:54:00
% DurationCPUTime: 14.20s
% Computational Cost: add. (12511->674), mult. (18081->946), div. (0->0), fcn. (16955->8), ass. (0->336)
t253 = pkin(7) + qJ(3);
t239 = sin(t253);
t240 = cos(t253);
t257 = sin(qJ(4));
t260 = cos(qJ(1));
t391 = t260 * t257;
t258 = sin(qJ(1));
t259 = cos(qJ(4));
t394 = t258 * t259;
t184 = t240 * t394 - t391;
t171 = Icges(5,4) * t184;
t393 = t259 * t260;
t395 = t257 * t258;
t183 = t240 * t395 + t393;
t399 = t239 * t258;
t101 = -Icges(5,2) * t183 + Icges(5,6) * t399 + t171;
t170 = Icges(5,4) * t183;
t105 = -Icges(5,1) * t184 - Icges(5,5) * t399 + t170;
t468 = t101 * t257 + t105 * t259;
t98 = Icges(5,5) * t184 - Icges(5,6) * t183 + Icges(5,3) * t399;
t38 = -t239 * t468 - t240 * t98;
t437 = t240 * pkin(3);
t204 = pkin(6) * t239 + t437;
t180 = t204 * t258;
t244 = t260 * qJ(2);
t211 = pkin(1) * t258 - t244;
t255 = cos(pkin(7));
t231 = pkin(2) * t255 + pkin(1);
t256 = -pkin(5) - qJ(2);
t234 = t260 * t256;
t372 = -t258 * t231 - t234;
t150 = t211 + t372;
t380 = t150 - t211;
t320 = rSges(5,1) * t184 - rSges(5,2) * t183;
t107 = rSges(5,3) * t399 + t320;
t319 = rSges(5,1) * t259 - rSges(5,2) * t257;
t148 = -rSges(5,3) * t240 + t239 * t319;
t361 = qJD(4) * t239;
t362 = qJD(3) * t260;
t194 = -t258 * t361 + t362;
t438 = t239 * pkin(3);
t203 = -pkin(6) * t240 + t438;
t360 = qJD(4) * t240;
t223 = qJD(1) - t360;
t241 = qJD(2) * t258;
t469 = -t107 * t223 - t148 * t194 - t203 * t362 + t241;
t33 = (-t180 + t380) * qJD(1) + t469;
t470 = t260 * t33;
t363 = qJD(3) * t258;
t193 = t260 * t361 + t363;
t25 = -t101 * t183 - t105 * t184 + t399 * t98;
t185 = -t240 * t391 + t394;
t186 = t240 * t393 + t395;
t398 = t239 * t260;
t100 = Icges(5,5) * t186 + Icges(5,6) * t185 + Icges(5,3) * t398;
t413 = Icges(5,4) * t186;
t103 = Icges(5,2) * t185 + Icges(5,6) * t398 + t413;
t172 = Icges(5,4) * t185;
t106 = Icges(5,1) * t186 + Icges(5,5) * t398 + t172;
t26 = t100 * t399 - t183 * t103 + t184 * t106;
t307 = Icges(5,5) * t259 - Icges(5,6) * t257;
t142 = -Icges(5,3) * t240 + t239 * t307;
t411 = Icges(5,4) * t259;
t308 = -Icges(5,2) * t257 + t411;
t144 = -Icges(5,6) * t240 + t239 * t308;
t412 = Icges(5,4) * t257;
t310 = Icges(5,1) * t259 - t412;
t146 = -Icges(5,5) * t240 + t239 * t310;
t51 = t142 * t399 - t144 * t183 + t146 * t184;
t10 = t193 * t26 - t194 * t25 + t223 * t51;
t27 = t185 * t101 - t105 * t186 + t98 * t398;
t28 = t100 * t398 + t185 * t103 + t186 * t106;
t52 = t142 * t398 + t144 * t185 + t146 * t186;
t11 = t193 * t28 - t194 * t27 + t52 * t223;
t462 = 0.2e1 * qJD(3);
t206 = qJD(1) * t211;
t460 = qJD(1) * t150 - t206;
t248 = t258 * rSges(4,3);
t396 = t240 * t260;
t159 = rSges(4,1) * t396 - rSges(4,2) * t398 + t248;
t219 = t260 * t231;
t330 = -t256 * t258 + t219;
t459 = t159 + t330;
t243 = t258 * qJ(2);
t213 = t260 * pkin(1) + t243;
t429 = rSges(3,2) * sin(pkin(7));
t431 = rSges(3,1) * t255;
t297 = t258 * rSges(3,3) + (-t429 + t431) * t260;
t458 = t213 + t297;
t230 = Icges(4,4) * t240;
t309 = -Icges(4,2) * t239 + t230;
t199 = Icges(4,1) * t239 + t230;
t196 = Icges(4,5) * t240 - Icges(4,6) * t239;
t195 = Icges(4,5) * t239 + Icges(4,6) * t240;
t286 = qJD(3) * t195;
t414 = Icges(4,4) * t239;
t200 = Icges(4,1) * t240 - t414;
t157 = Icges(4,5) * t258 + t200 * t260;
t155 = Icges(4,6) * t258 + t260 * t309;
t402 = t155 * t239;
t301 = -t157 * t240 + t402;
t405 = Icges(4,3) * t260;
t456 = -t260 * t286 + (-t196 * t258 + t301 + t405) * qJD(1);
t217 = Icges(4,4) * t399;
t397 = t240 * t258;
t410 = Icges(4,5) * t260;
t156 = Icges(4,1) * t397 - t217 - t410;
t407 = Icges(4,6) * t260;
t154 = Icges(4,4) * t397 - Icges(4,2) * t399 - t407;
t403 = t154 * t239;
t302 = -t156 * t240 + t403;
t153 = Icges(4,3) * t258 + t196 * t260;
t367 = qJD(1) * t153;
t455 = qJD(1) * t302 - t258 * t286 + t367;
t152 = Icges(4,5) * t397 - Icges(4,6) * t399 - t405;
t57 = -t260 * t152 - t258 * t302;
t329 = qJD(1) * t240 - qJD(4);
t345 = t239 * t362;
t454 = t258 * t329 + t345;
t197 = Icges(4,2) * t240 + t414;
t300 = t197 * t239 - t240 * t199;
t453 = t300 * qJD(1) + t196 * qJD(3);
t452 = t258 * (-t197 * t260 + t157) - t260 * (-Icges(4,2) * t397 + t156 - t217);
t143 = Icges(5,3) * t239 + t240 * t307;
t303 = -t144 * t257 + t146 * t259;
t305 = -t103 * t257 + t106 * t259;
t451 = t193 * (-t142 * t260 - t305) - t194 * (-t142 * t258 + t468) + t223 * (t143 - t303);
t164 = (-Icges(5,2) * t259 - t412) * t239;
t450 = t193 * (-Icges(5,2) * t186 + t106 + t172) - t194 * (-Icges(5,2) * t184 - t105 - t170) + t223 * (t146 + t164);
t357 = qJD(3) * qJD(4);
t342 = t240 * t357;
t137 = qJD(1) * t193 + t258 * t342;
t449 = t137 / 0.2e1;
t138 = qJD(1) * t194 + t260 * t342;
t448 = t138 / 0.2e1;
t447 = -t193 / 0.2e1;
t446 = t193 / 0.2e1;
t445 = -t194 / 0.2e1;
t444 = t194 / 0.2e1;
t443 = -t223 / 0.2e1;
t442 = t223 / 0.2e1;
t441 = t258 / 0.2e1;
t440 = -t260 / 0.2e1;
t439 = -rSges(5,3) - pkin(6);
t348 = t240 * t363;
t365 = qJD(1) * t260;
t283 = t239 * t365 + t348;
t298 = t223 * t259;
t346 = t239 * t363;
t85 = t258 * t298 + (-t260 * t329 + t346) * t257;
t299 = t223 * t257;
t364 = qJD(3) * t239;
t86 = t329 * t393 + (-t259 * t364 + t299) * t258;
t44 = Icges(5,5) * t86 + Icges(5,6) * t85 + Icges(5,3) * t283;
t46 = Icges(5,4) * t86 + Icges(5,2) * t85 + Icges(5,6) * t283;
t48 = Icges(5,1) * t86 + Icges(5,4) * t85 + Icges(5,5) * t283;
t7 = (-qJD(3) * t468 - t44) * t240 + (qJD(3) * t98 - t257 * t46 + t259 * t48 + (-t101 * t259 + t105 * t257) * qJD(4)) * t239;
t436 = t7 * t194;
t347 = t240 * t362;
t366 = qJD(1) * t258;
t352 = t239 * t366;
t282 = t347 - t352;
t83 = t454 * t257 + t260 * t298;
t84 = -t454 * t259 + t260 * t299;
t43 = Icges(5,5) * t84 + Icges(5,6) * t83 + Icges(5,3) * t282;
t45 = Icges(5,4) * t84 + Icges(5,2) * t83 + Icges(5,6) * t282;
t47 = Icges(5,1) * t84 + Icges(5,4) * t83 + Icges(5,5) * t282;
t8 = (qJD(3) * t305 - t43) * t240 + (qJD(3) * t100 - t257 * t45 + t259 * t47 + (-t103 * t259 - t106 * t257) * qJD(4)) * t239;
t435 = t8 * t193;
t434 = qJD(1) / 0.2e1;
t433 = pkin(1) - t231;
t161 = (-Icges(5,5) * t257 - Icges(5,6) * t259) * t239;
t80 = qJD(3) * t143 + qJD(4) * t161;
t145 = Icges(5,6) * t239 + t240 * t308;
t81 = qJD(3) * t145 + qJD(4) * t164;
t147 = Icges(5,5) * t239 + t240 * t310;
t167 = (-Icges(5,1) * t257 - t411) * t239;
t82 = qJD(3) * t147 + qJD(4) * t167;
t18 = (qJD(3) * t303 - t80) * t240 + (qJD(3) * t142 - t257 * t81 + t259 * t82 + (-t144 * t259 - t146 * t257) * qJD(4)) * t239;
t343 = t239 * t357;
t54 = -t142 * t240 + t239 * t303;
t432 = t18 * t223 + t54 * t343;
t430 = rSges(4,1) * t240;
t428 = rSges(5,3) * t239;
t109 = t186 * rSges(5,1) + t185 * rSges(5,2) + rSges(5,3) * t398;
t209 = pkin(6) * t347;
t284 = -t240 * t366 - t345;
t112 = pkin(3) * t284 - pkin(6) * t352 + t209;
t192 = qJD(3) * t204;
t235 = qJ(2) * t365;
t358 = qJD(1) * qJD(2);
t368 = t235 + t241;
t376 = qJD(1) * (-pkin(1) * t366 + t368) + t258 * t358;
t353 = qJD(1) * (-t235 + (t258 * t433 - t234) * qJD(1)) + t376;
t356 = t84 * rSges(5,1) + t83 * rSges(5,2) + rSges(5,3) * t347;
t49 = -rSges(5,3) * t352 + t356;
t149 = t240 * t319 + t428;
t173 = (-rSges(5,1) * t257 - rSges(5,2) * t259) * t239;
t89 = qJD(3) * t149 + qJD(4) * t173;
t13 = qJD(1) * t112 - t138 * t148 - t193 * t89 + t223 * t49 + (t109 * t361 - t192 * t258 - t203 * t365) * qJD(3) + t353;
t426 = t13 * t260;
t113 = t283 * pkin(6) + (t240 * t365 - t346) * pkin(3);
t233 = t260 * t358;
t349 = t203 * t363;
t242 = qJD(2) * t260;
t190 = t213 * qJD(1) - t242;
t227 = t256 * t366;
t383 = t227 - (-t260 * t433 - t243) * qJD(1) - t190;
t325 = rSges(5,1) * t86 + rSges(5,2) * t85;
t50 = rSges(5,3) * t283 + t325;
t14 = t137 * t148 - t194 * t89 - t223 * t50 + t233 + (-t107 * t361 - t192 * t260) * qJD(3) + (-t113 + t349 + t383) * qJD(1);
t425 = t14 * t258;
t201 = rSges(4,1) * t239 + rSges(4,2) * t240;
t175 = t201 * t260;
t351 = t201 * t363;
t62 = t459 * qJD(1) - t242 - t351;
t424 = t175 * t62;
t182 = pkin(3) * t396 + pkin(6) * t398;
t34 = -t349 + t109 * t223 - t148 * t193 - t242 + (t182 + t330) * qJD(1);
t420 = t258 * t34;
t158 = rSges(4,1) * t397 - rSges(4,2) * t399 - t260 * rSges(4,3);
t350 = t201 * t362;
t324 = t241 - t350;
t61 = (-t158 + t380) * qJD(1) + t324;
t419 = t258 * t61;
t418 = t38 * t137;
t39 = -t100 * t240 + t239 * t305;
t417 = t39 * t138;
t416 = -t192 - t89;
t401 = t195 * t258;
t400 = t195 * t260;
t70 = -t258 * t300 - t400;
t390 = t70 * qJD(1);
t387 = t107 + t180;
t386 = t109 + t182;
t385 = -t258 * t152 - t156 * t396;
t384 = t258 * t153 + t157 * t396;
t381 = t148 + t203;
t375 = -t197 + t200;
t374 = t199 + t309;
t373 = rSges(4,2) * t352 + rSges(4,3) * t365;
t228 = t258 * t429;
t371 = rSges(3,3) * t365 + qJD(1) * t228;
t370 = t227 + t242;
t369 = t260 * rSges(3,3) + t228;
t359 = t196 * qJD(1);
t355 = t258 * t431;
t344 = -pkin(1) - t431;
t340 = t365 / 0.2e1;
t339 = t364 / 0.2e1;
t338 = -t363 / 0.2e1;
t337 = t363 / 0.2e1;
t335 = t362 / 0.2e1;
t133 = t157 * t397;
t332 = t260 * t153 - t133;
t331 = -t152 + t402;
t326 = qJD(4) * t339;
t321 = -rSges(4,2) * t239 + t430;
t318 = t25 * t260 - t258 * t26;
t317 = t25 * t258 + t26 * t260;
t316 = t258 * t28 - t260 * t27;
t315 = t258 * t27 + t260 * t28;
t314 = t258 * t39 - t260 * t38;
t313 = t258 * t38 + t260 * t39;
t312 = -t258 * t62 - t260 * t61;
t304 = t107 * t260 - t109 * t258;
t72 = t154 * t240 + t156 * t239;
t73 = t155 * t240 + t157 * t239;
t174 = t201 * t258;
t58 = -t155 * t399 - t332;
t290 = (t258 * t58 - t260 * t57) * qJD(3);
t59 = -t154 * t398 - t385;
t60 = -t155 * t398 + t384;
t289 = (t258 * t60 - t260 * t59) * qJD(3);
t288 = qJD(3) * t199;
t287 = qJD(3) * t197;
t74 = (t158 * t258 + t159 * t260) * qJD(3);
t285 = -t204 - t428;
t281 = t100 * t193 + t142 * t223 - t194 * t98;
t280 = (-Icges(5,5) * t183 - Icges(5,6) * t184) * t194 - (Icges(5,5) * t185 - Icges(5,6) * t186) * t193 - t161 * t223;
t279 = t154 * t260 - t155 * t258;
t278 = t239 * t280;
t273 = (-t239 * t374 + t240 * t375) * qJD(1);
t271 = (Icges(5,1) * t185 - t103 - t413) * t193 - (-Icges(5,1) * t183 - t101 - t171) * t194 + (-t144 + t167) * t223;
t35 = t107 * t193 + t109 * t194 + (t180 * t258 + t182 * t260) * qJD(3);
t267 = t35 * t304 + (t258 * t33 - t260 * t34) * t148;
t93 = qJD(1) * t155 - t258 * t287;
t95 = qJD(1) * t157 - t258 * t288;
t266 = qJD(1) * t152 - qJD(3) * t72 - t239 * t93 + t240 * t95;
t92 = -t260 * t287 + (-t258 * t309 + t407) * qJD(1);
t94 = -t260 * t288 + (-t200 * t258 + t410) * qJD(1);
t265 = -qJD(3) * t73 - t239 * t92 + t240 * t94 + t367;
t188 = t309 * qJD(3);
t189 = t200 * qJD(3);
t264 = qJD(1) * t195 - t188 * t239 + t189 * t240 + (-t197 * t240 - t199 * t239) * qJD(3);
t263 = -t452 * t239 + t279 * t240;
t262 = t451 * t239;
t191 = t321 * qJD(3);
t181 = t203 * t260;
t179 = t203 * t258;
t160 = t355 - t369;
t131 = t148 * t260;
t130 = t148 * t258;
t129 = t146 * t260;
t128 = t146 * t258;
t127 = t144 * t260;
t126 = t144 * t258;
t123 = t458 * qJD(1) - t242;
t122 = t241 + (-t160 - t211) * qJD(1);
t121 = rSges(5,1) * t185 - rSges(5,2) * t186;
t120 = -rSges(5,1) * t183 - rSges(5,2) * t184;
t97 = -qJD(3) * t174 + (t260 * t321 + t248) * qJD(1);
t96 = rSges(4,1) * t284 - rSges(4,2) * t347 + t373;
t88 = t233 + (-t297 * qJD(1) - t190) * qJD(1);
t87 = qJD(1) * (-qJD(1) * t355 + t371) + t376;
t71 = -t260 * t300 + t401;
t69 = t71 * qJD(1);
t41 = -t191 * t362 + t233 + (-t97 + t351 + t383) * qJD(1);
t40 = -t191 * t363 + (t96 - t350) * qJD(1) + t353;
t37 = -qJD(3) * t301 + t239 * t94 + t240 * t92;
t36 = -t302 * qJD(3) + t239 * t95 + t240 * t93;
t32 = t264 * t258 - t453 * t260;
t31 = t453 * t258 + t264 * t260;
t24 = t69 + t289;
t23 = t290 + t390;
t16 = t142 * t283 + t144 * t85 + t146 * t86 - t183 * t81 + t184 * t82 + t399 * t80;
t15 = t142 * t282 + t144 * t83 + t146 * t84 + t185 * t81 + t186 * t82 + t398 * t80;
t12 = t193 * t39 - t194 * t38 + t223 * t54;
t9 = t107 * t138 - t109 * t137 + t193 * t50 + t194 * t49 + (t112 * t260 + t113 * t258 + (t180 * t260 - t182 * t258) * qJD(1)) * qJD(3);
t6 = t100 * t283 + t103 * t85 + t106 * t86 - t183 * t45 + t184 * t47 + t399 * t43;
t5 = t98 * t348 + t101 * t85 - t105 * t86 - t183 * t46 + t184 * t48 + (t258 * t44 + t365 * t98) * t239;
t4 = t100 * t282 + t103 * t83 + t106 * t84 + t185 * t45 + t186 * t47 + t398 * t43;
t3 = t98 * t347 + t101 * t83 - t105 * t84 + t185 * t46 + t186 * t48 + (t260 * t44 - t366 * t98) * t239;
t2 = t137 * t25 + t138 * t26 + t16 * t223 + t193 * t6 - t194 * t5 + t343 * t51;
t1 = t137 * t27 + t138 * t28 + t15 * t223 + t193 * t4 - t194 * t3 + t343 * t52;
t17 = [-t436 / 0.2e1 + (t69 + ((t58 - t133 + (t153 + t403) * t260 + t385) * t260 + t384 * t258) * qJD(3)) * t335 + t432 + t417 / 0.2e1 + t418 / 0.2e1 + t435 / 0.2e1 + t11 * t444 + t15 * t446 + t52 * t448 + t51 * t449 + (-qJD(3) * t300 + t188 * t240 + t189 * t239) * qJD(1) + (t16 + t11) * t445 + (-t390 + ((t260 * t331 - t384 + t60) * t260 + (t258 * t331 + t332 + t59) * t258) * qJD(3) + t23) * t338 + (t37 + t31) * t337 + (t14 * (-t320 + t372) + t33 * (-t325 + t370) + t13 * (t219 + t386) + (t14 * t285 - t13 * t256 + t33 * (t240 * t439 + t438) * qJD(3)) * t258 + ((t239 * t439 - t231 - t437) * t420 + (-t231 + t285) * t470) * qJD(1) + (-pkin(3) * t345 + t209 + t241 + t33 + t356 - t460 + (t180 - t234) * qJD(1) - t469) * t34) * m(5) + (t41 * (-t158 + t372) + t61 * t370 + t40 * t459 + t62 * (t241 + t373) + (t201 * t419 - t424) * qJD(3) + ((-t61 * rSges(4,3) + t62 * (-t231 - t430)) * t258 + (t61 * (-t231 - t321) - t62 * t256) * t260) * qJD(1) - (-qJD(1) * t158 + t324 + t460 - t61) * t62) * m(4) + (t88 * (t258 * t344 + t244 + t369) + t122 * t242 + t87 * t458 + t123 * (t368 + t371) + (t122 * (t344 + t429) * t260 + (t122 * (-rSges(3,3) - qJ(2)) + t123 * t344) * t258) * qJD(1) - (-qJD(1) * t160 - t122 - t206 + t241) * t123) * m(3) - (t36 + t32 + t24) * t362 / 0.2e1 + ((t72 + t70) * t258 + (t73 + t71) * t260) * qJD(3) * t434; 0.2e1 * (-t426 / 0.2e1 + t425 / 0.2e1) * m(5) + 0.2e1 * (t40 * t440 + t41 * t441) * m(4) + 0.2e1 * (t440 * t87 + t441 * t88) * m(3); t314 * t326 - t137 * t318 / 0.2e1 - qJD(1) * ((t239 * t375 + t240 * t374) * qJD(1) + (t279 * t239 + t452 * t240) * qJD(3)) / 0.2e1 + ((-t363 * t400 + t359) * t258 + (t273 + (t258 * t401 + t263) * qJD(3)) * t260) * t338 - t12 * t361 / 0.2e1 + ((-t362 * t401 - t359) * t260 + (t273 + (t260 * t400 + t263) * qJD(3)) * t258) * t335 + (t258 * t37 - t260 * t36 + (t72 * t258 + t260 * t73) * qJD(1)) * t434 + (qJD(1) * t313 + t258 * t8 - t260 * t7) * t442 + (((t127 * t257 - t129 * t259 + t100) * t193 - (t126 * t257 - t128 * t259 + t98) * t194 + (-t145 * t257 + t147 * t259 + t142) * t223 + t54 * qJD(4)) * t239 + (qJD(4) * t313 - t451) * t240) * t443 + ((t127 * t183 - t129 * t184) * t193 - (t126 * t183 - t128 * t184) * t194 + (-t145 * t183 + t147 * t184) * t223 + (t239 * t51 + t26 * t396) * qJD(4) + ((qJD(4) * t25 + t281) * t240 + t262) * t258) * t444 + (qJD(1) * t317 + t258 * t6 - t260 * t5) * t445 + (qJD(1) * t315 + t258 * t4 - t260 * t3) * t446 + ((-t127 * t185 - t129 * t186) * t193 - (-t126 * t185 - t128 * t186) * t194 + (t145 * t185 + t147 * t186) * t223 + (t239 * t52 + t27 * t397) * qJD(4) + ((qJD(4) * t28 + t281) * t240 + t262) * t260) * t447 + t316 * t448 - (t10 * t258 + t11 * t260) * t360 / 0.2e1 + ((-t14 * t381 + t33 * t416 + t9 * t386 + t35 * (t112 + t49) + (-t34 * t381 + t35 * t387) * qJD(1)) * t260 + (-t13 * t381 + t34 * t416 + t9 * t387 + t35 * (t113 + t50) + (t33 * t381 - t35 * t386) * qJD(1)) * t258 - t33 * (qJD(1) * t179 + t130 * t223 - t149 * t194 - t204 * t362) - t34 * (-qJD(1) * t181 - t131 * t223 - t149 * t193 - t204 * t363) - t35 * (-t130 * t193 - t131 * t194 - t179 * t363 - t181 * t362) - ((-t107 * t33 + t109 * t34) * t239 + t267 * t240) * qJD(4)) * m(5) + (-(t174 * t61 - t424) * qJD(1) - (t74 * (-t174 * t258 - t175 * t260) + t312 * t321) * qJD(3) + 0.2e1 * t74 * (t258 * t97 + t260 * t96 + (t158 * t260 - t159 * t258) * qJD(1)) + t312 * t191 + (-t40 * t258 - t41 * t260 + (-t260 * t62 + t419) * qJD(1)) * t201) * m(4) + (qJD(1) * t31 + t1 + (-(t455 * t258 + t266 * t260) * t260 + (t456 * t258 + t265 * t260) * t258 + (t59 * t258 + t60 * t260) * qJD(1)) * t462) * t441 + (qJD(1) * t32 + t2 + (-(t266 * t258 - t455 * t260) * t260 + (t265 * t258 - t456 * t260) * t258 + (t57 * t258 + t58 * t260) * qJD(1)) * t462) * t440 + (t290 + t23 + t10) * t366 / 0.2e1 + (t289 + t24 + t11) * t340; t1 * t398 / 0.2e1 + (t239 * t315 - t240 * t52) * t448 + ((qJD(3) * t315 - t15) * t240 + (-qJD(1) * t316 + qJD(3) * t52 + t258 * t3 + t260 * t4) * t239) * t446 + t2 * t399 / 0.2e1 + (t239 * t317 - t240 * t51) * t449 + ((qJD(3) * t317 - t16) * t240 + (qJD(1) * t318 + qJD(3) * t51 + t258 * t5 + t260 * t6) * t239) * t445 + t12 * t339 - t240 * (t417 + t418 + t432 + t435 - t436) / 0.2e1 + (t239 * t313 - t240 * t54) * t326 + ((qJD(3) * t313 - t18) * t240 + (-qJD(1) * t314 + qJD(3) * t54 + t258 * t7 + t260 * t8) * t239) * t442 + (t450 * t185 + t271 * t186 - t260 * t278) * t447 + (-t183 * t450 + t184 * t271 - t258 * t278) * t444 + (t280 * t240 + (-t257 * t450 + t259 * t271) * t239) * t443 + (-t352 / 0.2e1 + t240 * t335) * t11 + (t239 * t340 + t240 * t337) * t10 + ((qJD(3) * t267 + t14 * t107 - t13 * t109 + t33 * t50 - t34 * t49) * t240 + (t33 * (-qJD(3) * t107 + t258 * t89) + t34 * (qJD(3) * t109 - t260 * t89) + t9 * t304 + t35 * (-t107 * t366 - t109 * t365 - t258 * t49 + t260 * t50) + (-t426 + t425 + (t420 + t470) * qJD(1)) * t148) * t239 - t33 * (-t120 * t223 - t173 * t194) - t34 * (t121 * t223 - t173 * t193) - t35 * (t120 * t193 + t121 * t194)) * m(5);];
tauc = t17(:);
