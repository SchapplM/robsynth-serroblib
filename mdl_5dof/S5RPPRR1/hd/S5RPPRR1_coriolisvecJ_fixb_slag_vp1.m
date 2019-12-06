% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:12
% DurationCPUTime: 7.60s
% Computational Cost: add. (6097->530), mult. (9982->716), div. (0->0), fcn. (7534->6), ass. (0->313)
t229 = sin(qJ(1));
t231 = cos(qJ(1));
t227 = qJ(4) + qJ(5);
t218 = cos(t227);
t217 = sin(t227);
t385 = Icges(6,4) * t217;
t277 = Icges(6,2) * t218 + t385;
t117 = Icges(6,6) * t231 + t229 * t277;
t163 = Icges(6,1) * t218 - t385;
t441 = t163 * t229 - t117;
t118 = -Icges(6,6) * t229 + t231 * t277;
t440 = t163 * t231 - t118;
t384 = Icges(6,4) * t218;
t279 = Icges(6,1) * t217 + t384;
t119 = Icges(6,5) * t231 + t229 * t279;
t161 = -Icges(6,2) * t217 + t384;
t439 = t161 * t229 + t119;
t438 = -t161 - t279;
t437 = -t277 + t163;
t367 = t218 * t231;
t191 = Icges(6,4) * t367;
t368 = t217 * t231;
t381 = Icges(6,5) * t229;
t120 = Icges(6,1) * t368 + t191 - t381;
t436 = t118 * t218 + t120 * t217;
t335 = t229 * rSges(4,2) + t231 * rSges(4,3);
t188 = t231 * pkin(1) + t229 * qJ(2);
t431 = t231 * qJ(3) + t188;
t435 = t431 + t335;
t434 = 2 * qJD(4);
t433 = t117 * t367 + t119 * t368;
t303 = -rSges(3,2) * t231 + t229 * rSges(3,3);
t432 = t188 + t303;
t275 = Icges(6,5) * t217 + Icges(6,6) * t218;
t115 = Icges(6,3) * t231 + t229 * t275;
t430 = qJD(1) * t115;
t286 = rSges(6,1) * t217 + rSges(6,2) * t218;
t121 = rSges(6,3) * t231 + t229 * t286;
t429 = qJD(1) * t121;
t228 = sin(qJ(4));
t230 = cos(qJ(4));
t276 = Icges(5,5) * t228 + Icges(5,6) * t230;
t129 = Icges(5,3) * t231 + t229 * t276;
t428 = qJD(1) * t129;
t397 = rSges(6,2) * t217;
t400 = rSges(6,1) * t218;
t165 = -t397 + t400;
t144 = t165 * t229;
t145 = t165 * t231;
t226 = qJD(4) + qJD(5);
t175 = t226 * t229;
t176 = t226 * t231;
t345 = rSges(6,1) * t368 + rSges(6,2) * t367;
t393 = rSges(6,3) * t229;
t122 = t345 - t393;
t232 = -pkin(7) - pkin(6);
t365 = t228 * t231;
t341 = pkin(4) * t365 + t229 * t232;
t146 = pkin(6) * t229 + t341;
t366 = t228 * t229;
t203 = pkin(4) * t366;
t147 = -t203 + (pkin(6) + t232) * t231;
t41 = -t121 * t175 - t122 * t176 + (-t146 * t231 + t147 * t229) * qJD(4);
t221 = t231 * qJ(2);
t184 = pkin(1) * t229 - t221;
t375 = qJ(3) * t229;
t404 = pkin(6) * t231;
t285 = -t375 - t404;
t353 = -t121 + t147;
t254 = t285 + t353;
t326 = qJD(4) * t231;
t312 = t230 * t326;
t195 = pkin(4) * t312;
t214 = qJD(3) * t231;
t215 = qJD(2) * t229;
t337 = t214 + t215;
t266 = t165 * t176 + t195 + t337;
t42 = (-t184 + t254) * qJD(1) + t266;
t327 = qJD(4) * t230;
t289 = pkin(4) * t327 + qJD(3);
t332 = qJD(1) * t229;
t213 = pkin(6) * t332;
t216 = qJD(2) * t231;
t338 = t213 + t216;
t352 = -t122 - t146;
t43 = t165 * t175 + t289 * t229 + (t431 - t352) * qJD(1) - t338;
t427 = -t43 * (qJD(1) * t145 - t175 * t286) - t41 * (-t175 * t144 - t145 * t176) - t42 * (-qJD(1) * t144 - t176 * t286);
t178 = Icges(5,5) * t230 - Icges(5,6) * t228;
t255 = qJD(4) * t178;
t362 = t230 * t231;
t198 = Icges(5,4) * t362;
t383 = Icges(5,5) * t229;
t134 = Icges(5,1) * t365 + t198 - t383;
t387 = Icges(5,4) * t228;
t278 = Icges(5,2) * t230 + t387;
t132 = -Icges(5,6) * t229 + t231 * t278;
t369 = t132 * t230;
t271 = t134 * t228 + t369;
t426 = qJD(1) * t271 + t231 * t255 - t428;
t131 = Icges(5,6) * t231 + t229 * t278;
t386 = Icges(5,4) * t230;
t280 = Icges(5,1) * t228 + t386;
t133 = Icges(5,5) * t231 + t229 * t280;
t272 = t131 * t230 + t133 * t228;
t130 = -Icges(5,3) * t229 + t231 * t276;
t333 = qJD(1) * t130;
t425 = qJD(1) * t272 + t229 * t255 + t333;
t159 = Icges(6,5) * t218 - Icges(6,6) * t217;
t139 = t159 * t231;
t424 = qJD(1) * t436 + t226 * t139 - t430;
t138 = t159 * t229;
t274 = t117 * t218 + t119 * t217;
t116 = -Icges(6,3) * t229 + t231 * t275;
t334 = qJD(1) * t116;
t423 = qJD(1) * t274 + t138 * t226 + t334;
t287 = rSges(5,1) * t228 + rSges(5,2) * t230;
t135 = rSges(5,3) * t231 + t229 * t287;
t270 = t161 * t218 + t163 * t217;
t422 = qJD(1) * t270 - t275 * t226;
t180 = -Icges(5,2) * t228 + t386;
t182 = Icges(5,1) * t230 - t387;
t269 = t230 * t180 + t182 * t228;
t421 = t269 * qJD(1) - t276 * qJD(4);
t420 = t229 * (-Icges(5,2) * t365 + t134 + t198) - t231 * (t180 * t229 + t133);
t419 = qJD(1) * t437 - t175 * (-Icges(6,2) * t368 + t120 + t191) + t176 * t439;
t293 = qJD(1) * t226;
t166 = t229 * t293;
t418 = -t166 / 0.2e1;
t167 = t231 * t293;
t417 = -t167 / 0.2e1;
t416 = t175 / 0.2e1;
t415 = -t175 / 0.2e1;
t414 = -t176 / 0.2e1;
t413 = t176 / 0.2e1;
t412 = -t229 / 0.2e1;
t411 = t229 / 0.2e1;
t410 = -t231 / 0.2e1;
t409 = t231 / 0.2e1;
t408 = rSges(3,2) - pkin(1);
t407 = -rSges(5,3) - pkin(6);
t406 = pkin(4) * t228;
t403 = -qJD(1) / 0.2e1;
t402 = qJD(1) / 0.2e1;
t401 = -pkin(1) - qJ(3);
t398 = rSges(4,2) * t231;
t396 = rSges(3,3) * t231;
t395 = rSges(5,3) * t229;
t391 = -rSges(6,3) + t232;
t390 = t118 * t367 + t120 * t368;
t313 = t229 * t327;
t330 = qJD(1) * t232;
t331 = qJD(1) * t231;
t101 = t229 * t330 + t213 + (t228 * t331 + t313) * pkin(4);
t71 = t226 * t144 + (t231 * t286 - t393) * qJD(1);
t389 = -t101 - t71;
t344 = t231 * t330 + t195;
t102 = (-t203 + t404) * qJD(1) + t344;
t321 = t226 * t400;
t171 = t231 * t321;
t320 = t226 * t397;
t70 = -t231 * t320 + t171 - t429;
t388 = -t102 - t70;
t234 = qJD(1) ^ 2;
t374 = qJ(3) * t234;
t371 = t129 * t229;
t370 = t129 * t231;
t149 = t178 * t229;
t150 = t178 * t231;
t364 = t229 * t115;
t363 = t229 * t234;
t361 = t231 * t115;
t123 = t231 * t130;
t360 = t231 * t234;
t58 = t229 * t270 + t139;
t359 = t58 * qJD(1);
t76 = t229 * t269 + t150;
t358 = t76 * qJD(1);
t357 = t131 * t362 + t133 * t365;
t356 = t132 * t362 + t134 * t365;
t324 = qJD(1) * qJD(2);
t207 = qJ(2) * t331;
t340 = t207 + t215;
t349 = qJD(1) * (-pkin(1) * t332 + t340) + t229 * t324;
t347 = -t278 + t182;
t346 = -t180 - t280;
t343 = rSges(5,1) * t365 + rSges(5,2) * t362;
t206 = t231 * t324;
t323 = qJD(1) * qJD(3);
t342 = -0.2e1 * t229 * t323 + t206;
t339 = rSges(3,2) * t332 + rSges(3,3) * t331;
t174 = qJD(1) * t184;
t336 = t215 - t174;
t329 = qJD(4) * t228;
t328 = qJD(4) * t229;
t325 = t276 * qJD(1);
t322 = -rSges(4,3) + t401;
t45 = t231 * t116 + t436 * t229;
t319 = rSges(5,2) * t329;
t49 = t134 * t366 + t229 * t369 + t123;
t318 = 0.2e1 * t231 * t323 + t349;
t317 = pkin(6) * t363 + t342;
t316 = t207 + t337;
t315 = t214 + t336;
t187 = rSges(5,1) * t230 - rSges(5,2) * t228;
t157 = t187 * t326;
t311 = -t332 / 0.2e1;
t310 = -t331 / 0.2e1;
t308 = t328 / 0.2e1;
t307 = -t326 / 0.2e1;
t305 = pkin(4) * t230 + t165;
t302 = (t231 * t279 - t381) * qJD(1) + t441 * t226;
t301 = -qJD(1) * t119 + t440 * t226;
t300 = qJD(1) * t118 + t439 * t226;
t299 = -qJD(1) * t117 + t120 * t226 + t161 * t176;
t297 = -rSges(4,3) * t229 - t375 + t398;
t296 = t437 * t226;
t295 = t438 * t226;
t128 = t286 * t226;
t290 = -pkin(4) * t329 - t128;
t263 = -pkin(6) * t360 + t318;
t265 = -(qJD(4) ^ 2) * t406 - t374;
t24 = -t128 * t175 + t165 * t167 + t265 * t229 + (t195 - t388) * qJD(1) + t263;
t148 = qJD(1) * t188 - t216;
t25 = -t128 * t176 - t165 * t166 + t265 * t231 + (-pkin(4) * t313 - t148 + t389) * qJD(1) + t317;
t284 = t24 * t229 + t25 * t231;
t172 = t287 * qJD(4);
t268 = -qJD(4) * t172 - t374;
t194 = rSges(5,1) * t312;
t89 = -qJD(1) * t135 - t231 * t319 + t194;
t36 = t268 * t229 + (t89 + t157) * qJD(1) + t263;
t155 = t187 * t229;
t90 = qJD(4) * t155 + (t231 * t287 - t395) * qJD(1);
t37 = t268 * t231 + (-t187 * t328 - t148 - t90) * qJD(1) + t317;
t283 = t36 * t229 + t37 * t231;
t282 = t229 * t43 + t231 * t42;
t264 = -t135 + t285;
t55 = t157 + (-t184 + t264) * qJD(1) + t337;
t136 = t343 - t395;
t56 = (qJD(4) * t187 + qJD(3)) * t229 + (t136 + t431) * qJD(1) - t338;
t281 = t229 * t56 + t231 * t55;
t60 = -t117 * t217 + t119 * t218;
t74 = -t131 * t228 + t133 * t230;
t75 = -t132 * t228 + t134 * t230;
t267 = -t45 + t364;
t48 = t229 * t272 + t370;
t261 = (-t229 * t49 + t231 * t48) * qJD(4);
t50 = t357 - t371;
t51 = -t130 * t229 + t356;
t260 = (-t229 * t51 + t231 * t50) * qJD(4);
t259 = -t287 + t401;
t258 = -t286 + t401;
t257 = qJD(4) * t182;
t256 = qJD(4) * t180;
t72 = (-t135 * t229 - t136 * t231) * qJD(4);
t253 = -qJD(1) * t275 + t138 * t176 - t139 * t175;
t252 = -t131 * t231 + t132 * t229;
t251 = t258 - t406;
t242 = t217 * t302 + t218 * t300 - t430;
t10 = t242 * t229 + t231 * t423;
t241 = t217 * t301 + t218 * t299 - t334;
t11 = t241 * t229 + t231 * t424;
t44 = t229 * t274 + t361;
t20 = -t175 * t45 + t176 * t44 + t359;
t46 = -t364 + t433;
t47 = -t116 * t229 + t390;
t59 = t231 * t270 - t138;
t57 = t59 * qJD(1);
t21 = -t175 * t47 + t176 * t46 + t57;
t246 = t438 * qJD(1) - t440 * t175 + t441 * t176;
t235 = t246 * t217 + t218 * t419;
t240 = -qJD(1) * t159 + t217 * t295 + t218 * t296;
t28 = -t229 * t422 + t240 * t231;
t29 = t240 * t229 + t231 * t422;
t30 = -t217 * t300 + t218 * t302;
t31 = -t217 * t299 + t218 * t301;
t61 = -t118 * t217 + t120 * t218;
t8 = -t229 * t423 + t242 * t231;
t9 = -t229 * t424 + t241 * t231;
t250 = (qJD(1) * t28 - t166 * t46 - t167 * t47 - t175 * t9 + t176 * t8) * t412 + (-t217 * t419 + t246 * t218) * t403 + t20 * t311 + t21 * t310 + (qJD(1) * t29 + t10 * t176 - t11 * t175 - t166 * t44 - t167 * t45) * t409 + (-t229 * t45 + t231 * t44) * t418 + (-t229 * t47 + t231 * t46) * t417 + (-t229 * t9 + t231 * t8 + (-t229 * t46 - t231 * t47) * qJD(1)) * t415 + (t10 * t231 - t11 * t229 + (-t229 * t44 - t231 * t45) * qJD(1)) * t413 + (-t229 * t31 + t231 * t30 + (-t229 * t60 - t231 * t61) * qJD(1)) * t402 + (-t229 * t253 + t231 * t235) * t416 + (t229 * t235 + t231 * t253) * t414;
t249 = (t228 * t346 + t230 * t347) * qJD(1);
t245 = -t361 + (-t116 - t274) * t229 + t390;
t84 = qJD(1) * t132 + t229 * t256;
t86 = t229 * t257 + (t231 * t280 - t383) * qJD(1);
t239 = qJD(4) * t74 + t228 * t86 + t230 * t84 - t428;
t83 = -qJD(1) * t131 + t231 * t256;
t85 = -qJD(1) * t133 + t231 * t257;
t238 = qJD(4) * t75 + t228 * t85 + t230 * t83 - t333;
t169 = t278 * qJD(4);
t170 = t280 * qJD(4);
t237 = -qJD(1) * t178 - t169 * t230 - t170 * t228 + (-t180 * t228 + t182 * t230) * qJD(4);
t236 = t252 * t228 - t230 * t420;
t210 = rSges(4,2) * t331;
t185 = rSges(3,2) * t229 + t396;
t156 = t187 * t231;
t106 = t122 * t332;
t105 = qJD(1) * t432 - t216;
t104 = t215 + (-t184 + t185) * qJD(1);
t92 = t435 * qJD(1) + qJD(3) * t229 - t216;
t91 = (-t184 + t297) * qJD(1) + t337;
t88 = t206 + (-qJD(1) * t303 - t148) * qJD(1);
t87 = qJD(1) * t339 + t349;
t77 = t231 * t269 - t149;
t73 = t77 * qJD(1);
t69 = -qJ(3) * t360 + (-qJD(1) * t335 - t148) * qJD(1) + t342;
t68 = -qJ(3) * t363 + qJD(1) * (-rSges(4,3) * t332 + t210) + t318;
t35 = t237 * t229 + t231 * t421;
t34 = -t229 * t421 + t237 * t231;
t33 = -qJD(4) * t271 - t228 * t83 + t230 * t85;
t32 = -qJD(4) * t272 - t228 * t84 + t230 * t86;
t23 = t73 + t260;
t22 = t261 + t358;
t13 = -t121 * t167 + t122 * t166 - t175 * t71 - t176 * t70 + (-t101 * t229 - t102 * t231 + (t146 * t229 + t147 * t231) * qJD(1)) * qJD(4);
t1 = [(t57 + (-t267 - t45 + t433) * t176 - (t44 + t245) * t175) * t414 + (t73 + (t357 * t231 + (-t48 + (t130 + t272) * t229 - t356) * t229) * qJD(4)) * t307 + (t60 + t58) * t418 + (t61 + t59) * t417 + (t20 - t359 + (-t47 + t245) * t176 - (-t231 * t274 + t267 + t46) * t175) * t416 + (t31 + t28) * t415 - (t33 + t34) * t328 / 0.2e1 + (t22 - t358 + ((t356 - t51 - t370) * t231 + (-t123 + t49 - t50 - t371) * t229) * qJD(4)) * t308 + (-t269 * qJD(4) + t169 * t228 - t170 * t230 - t296 * t217 + t295 * t218) * qJD(1) + (t25 * (-t203 + t221) + t42 * t216 + t24 * (t431 + t341 + t345) + t43 * (t171 + t316 + t344) + (t25 * t391 - t43 * t320 + (-t43 * rSges(6,3) + t251 * t42) * qJD(1)) * t231 + (t25 * t258 + t42 * (-t289 + t320 - t321) - t24 * rSges(6,3) + (t42 * (-qJ(2) - t391) + t43 * t251) * qJD(1)) * t229 - (t254 * qJD(1) - t174 + t266 - t42) * t43) * m(6) + (t37 * t221 + t55 * t338 + t36 * (t431 + t343) + t56 * (t194 + t316) + (t37 * t407 - t56 * t319 + (t259 * t55 + t407 * t56) * qJD(1)) * t231 + (t37 * t259 + t55 * (-rSges(5,1) * t327 - qJD(3) + t319) + t36 * t407 + (t55 * (rSges(5,3) - qJ(2)) + t56 * t259) * qJD(1)) * t229 - (qJD(1) * t264 + t157 + t315 - t55) * t56) * m(5) + (t69 * (t221 + t398) + t91 * t216 + t68 * t435 + t92 * (t210 + t316) + (-t91 * qJD(3) + t69 * t322) * t229 + (t91 * t322 * t231 + (t91 * (-rSges(4,2) - qJ(2)) + t92 * t322) * t229) * qJD(1) - (qJD(1) * t297 + t315 - t91) * t92) * m(4) + (t88 * (t229 * t408 + t221 + t396) + t104 * t216 + t87 * t432 + t105 * (t339 + t340) + (t104 * t408 * t231 + (t104 * (-rSges(3,3) - qJ(2)) - t105 * pkin(1)) * t229) * qJD(1) - (qJD(1) * t185 - t104 + t336) * t105) * m(3) + (t30 + t29 + t21) * t413 + (t32 + t35 + t23) * t326 / 0.2e1 + ((t74 + t76) * t229 + (t75 + t77) * t231) * qJD(4) * t403; 0.2e1 * (t24 * t410 + t25 * t411) * m(6) + 0.2e1 * (t36 * t410 + t37 * t411) * m(5) + 0.2e1 * (t410 * t68 + t411 * t69) * m(4) + 0.2e1 * (t410 * t87 + t411 * t88) * m(3); m(4) * (t229 * t68 + t231 * t69) + m(5) * t283 + m(6) * t284; (-t229 * t33 + t231 * t32 + (-t74 * t229 - t231 * t75) * qJD(1)) * t402 + ((-t228 * t347 + t230 * t346) * qJD(1) + (t228 * t420 + t252 * t230) * qJD(4)) * t403 + ((t149 * t326 - t325) * t231 + (t249 + (-t231 * t150 + t236) * qJD(4)) * t229) * t307 + ((t150 * t328 + t325) * t229 + (t249 + (-t229 * t149 + t236) * qJD(4)) * t231) * t308 + t250 + (qJD(1) * t34 + ((-t229 * t425 + t239 * t231) * t231 - (-t229 * t426 + t238 * t231) * t229 + (-t50 * t229 - t51 * t231) * qJD(1)) * t434) * t412 + (qJD(1) * t35 + ((t239 * t229 + t231 * t425) * t231 - (t238 * t229 + t231 * t426) * t229 + (-t48 * t229 - t49 * t231) * qJD(1)) * t434) * t409 + (t261 + t22) * t311 + (t260 + t23) * t310 + (-(-t42 * t230 * t332 + (t41 * (-t229 ^ 2 - t231 ^ 2) * t230 - t282 * t228) * qJD(4)) * pkin(4) + t41 * t106 + (t25 * t305 + t42 * t290 + t13 * t352 + t41 * t388 + (t43 * t165 + t353 * t41) * qJD(1)) * t231 + (t24 * t305 + t43 * t290 + t13 * t353 + t41 * t389 + (t41 * t146 - t305 * t42) * qJD(1)) * t229 + t427) * m(6) + (0.2e1 * t72 * (-t229 * t90 - t231 * t89 + (-t135 * t231 + t136 * t229) * qJD(1)) - t281 * t172 + ((-t229 * t55 + t231 * t56) * qJD(1) + t283) * t187 - (-t155 * t55 + t156 * t56) * qJD(1) - (t72 * (-t155 * t229 - t156 * t231) - t281 * t287) * qJD(4)) * m(5); t250 + (t13 * (-t121 * t229 - t122 * t231) + t41 * (-t229 * t71 + t106 + (-t70 - t429) * t231) - t282 * t128 + ((-t229 * t42 + t231 * t43) * qJD(1) + t284) * t165 + t427) * m(6);];
tauc = t1(:);
