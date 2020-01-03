% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR6
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:38
% DurationCPUTime: 6.28s
% Computational Cost: add. (9102->511), mult. (9454->686), div. (0->0), fcn. (7388->8), ass. (0->303)
t240 = pkin(7) + qJ(3);
t220 = qJ(4) + t240;
t210 = cos(t220);
t245 = sin(qJ(1));
t371 = t210 * t245;
t209 = sin(t220);
t373 = t209 * t245;
t246 = cos(qJ(1));
t383 = Icges(5,6) * t246;
t114 = Icges(5,4) * t371 - Icges(5,2) * t373 - t383;
t203 = Icges(5,4) * t210;
t162 = Icges(5,1) * t209 + t203;
t437 = -t162 * t245 - t114;
t283 = -Icges(5,2) * t209 + t203;
t115 = Icges(5,6) * t245 + t246 * t283;
t436 = -t162 * t246 - t115;
t389 = Icges(5,4) * t209;
t163 = Icges(5,1) * t210 - t389;
t117 = Icges(5,5) * t245 + t163 * t246;
t160 = Icges(5,2) * t210 + t389;
t435 = -t160 * t246 + t117;
t434 = -t160 + t163;
t433 = t162 + t283;
t432 = 2 * qJD(3);
t218 = sin(t240);
t431 = rSges(4,2) * t218;
t225 = t246 * qJ(2);
t244 = -pkin(5) - qJ(2);
t364 = t244 * t246;
t243 = cos(pkin(7));
t211 = t243 * pkin(2) + pkin(1);
t406 = pkin(1) - t211;
t274 = t245 * t406 - t364;
t120 = -t225 + t274;
t192 = pkin(1) * t245 - t225;
t181 = qJD(1) * t192;
t430 = qJD(1) * t120 - t181;
t233 = t245 * rSges(4,3);
t219 = cos(t240);
t366 = t219 * t246;
t368 = t218 * t246;
t133 = rSges(4,1) * t366 - rSges(4,2) * t368 + t233;
t199 = t246 * t211;
t298 = -t244 * t245 + t199;
t429 = t133 + t298;
t223 = t245 * qJ(2);
t194 = t246 * pkin(1) + t223;
t402 = rSges(3,2) * sin(pkin(7));
t405 = rSges(3,1) * t243;
t275 = t245 * rSges(3,3) + (-t402 + t405) * t246;
t428 = t194 + t275;
t208 = Icges(4,4) * t219;
t284 = -Icges(4,2) * t218 + t208;
t172 = Icges(4,1) * t218 + t208;
t169 = Icges(4,5) * t219 - Icges(4,6) * t218;
t168 = Icges(4,5) * t218 + Icges(4,6) * t219;
t267 = qJD(3) * t168;
t390 = Icges(4,4) * t218;
t173 = Icges(4,1) * t219 - t390;
t127 = Icges(4,5) * t245 + t173 * t246;
t125 = Icges(4,6) * t245 + t246 * t284;
t378 = t125 * t218;
t279 = -t127 * t219 + t378;
t382 = Icges(4,3) * t246;
t427 = -t246 * t267 + (-t169 * t245 + t279 + t382) * qJD(1);
t369 = t218 * t245;
t198 = Icges(4,4) * t369;
t367 = t219 * t245;
t388 = Icges(4,5) * t246;
t126 = Icges(4,1) * t367 - t198 - t388;
t384 = Icges(4,6) * t246;
t124 = Icges(4,4) * t367 - Icges(4,2) * t369 - t384;
t379 = t124 * t218;
t280 = -t126 * t219 + t379;
t123 = Icges(4,3) * t245 + t169 * t246;
t335 = qJD(1) * t123;
t426 = qJD(1) * t280 - t245 * t267 + t335;
t159 = Icges(5,5) * t210 - Icges(5,6) * t209;
t241 = qJD(3) + qJD(4);
t158 = Icges(5,5) * t209 + Icges(5,6) * t210;
t376 = t158 * t246;
t380 = t115 * t209;
t381 = Icges(5,3) * t246;
t425 = -t241 * t376 + (-t117 * t210 - t159 * t245 + t380 + t381) * qJD(1);
t185 = Icges(5,4) * t373;
t387 = Icges(5,5) * t246;
t116 = Icges(5,1) * t371 - t185 - t387;
t282 = t114 * t209 - t116 * t210;
t113 = Icges(5,3) * t245 + t159 * t246;
t336 = qJD(1) * t113;
t377 = t158 * t245;
t424 = qJD(1) * t282 - t241 * t377 + t336;
t122 = Icges(4,5) * t367 - Icges(4,6) * t369 - t382;
t49 = -t246 * t122 - t245 * t280;
t278 = t160 * t209 - t162 * t210;
t423 = qJD(1) * t278 + t159 * t241;
t170 = Icges(4,2) * t219 + t390;
t277 = t218 * t170 - t172 * t219;
t422 = t277 * qJD(1) + t169 * qJD(3);
t188 = rSges(5,2) * t373;
t118 = rSges(5,1) * t371 - t246 * rSges(5,3) - t188;
t232 = t245 * rSges(5,3);
t370 = t210 * t246;
t372 = t209 * t246;
t119 = rSges(5,1) * t370 - rSges(5,2) * t372 + t232;
t400 = rSges(5,2) * t210;
t164 = rSges(5,1) * t209 + t400;
t141 = t164 * t245;
t142 = t164 * t246;
t403 = rSges(5,1) * t210;
t165 = -rSges(5,2) * t209 + t403;
t190 = t241 * t245;
t191 = t241 * t246;
t409 = pkin(3) * t219;
t178 = t211 + t409;
t239 = -pkin(6) + t244;
t348 = -t245 * t178 - t246 * t239;
t90 = t211 * t245 + t348 + t364;
t167 = t246 * t178;
t337 = -t239 + t244;
t91 = t245 * t337 + t167 - t199;
t36 = t118 * t190 + t119 * t191 + (-t245 * t90 + t246 * t91) * qJD(3);
t221 = qJD(2) * t245;
t331 = qJD(3) * t246;
t320 = t218 * t331;
t296 = pkin(3) * t320;
t276 = t221 - t296;
t265 = -t191 * t164 + t276;
t354 = t120 - t192;
t394 = -t118 + t90;
t41 = (t354 + t394) * qJD(1) + t265;
t332 = qJD(3) * t245;
t410 = pkin(3) * t218;
t186 = t332 * t410;
t222 = qJD(2) * t246;
t342 = t186 + t222;
t393 = -t119 - t91;
t42 = -t164 * t190 + (t298 - t393) * qJD(1) - t342;
t295 = qJD(1) * t241;
t176 = t245 * t295;
t177 = t246 * t295;
t334 = qJD(1) * t245;
t333 = qJD(1) * t246;
t344 = rSges(5,3) * t333 + qJD(1) * t188;
t73 = -t191 * t400 + (-t191 * t209 - t210 * t334) * rSges(5,1) + t344;
t74 = -t241 * t141 + (t165 * t246 + t232) * qJD(1);
t345 = t178 - t211;
t76 = -t296 + (-t245 * t345 + t246 * t337) * qJD(1);
t204 = t244 * t334;
t365 = t239 * t245;
t77 = -t186 + t204 + (t246 * t345 - t365) * qJD(1);
t8 = t118 * t177 - t119 * t176 + t190 * t74 + t191 * t73 + (t245 * t77 + t246 * t76 + (-t245 * t91 - t246 * t90) * qJD(1)) * qJD(3);
t421 = -(qJD(1) * t141 - t191 * t165) * t41 - t36 * (-t190 * t141 - t142 * t191) - t42 * (-qJD(1) * t142 - t165 * t190) + t8 * (t245 * t118 + t246 * t119);
t420 = t245 * (-t170 * t246 + t127) - t246 * (-Icges(4,2) * t367 + t126 - t198);
t419 = qJD(1) * t433 + t190 * t435 - t191 * (-Icges(5,2) * t371 + t116 - t185);
t418 = t176 / 0.2e1;
t417 = t177 / 0.2e1;
t416 = -t190 / 0.2e1;
t415 = t190 / 0.2e1;
t414 = -t191 / 0.2e1;
t413 = t191 / 0.2e1;
t412 = t245 / 0.2e1;
t411 = -t246 / 0.2e1;
t408 = -qJD(1) / 0.2e1;
t407 = qJD(1) / 0.2e1;
t404 = rSges(4,1) * t219;
t401 = rSges(4,2) * t219;
t174 = rSges(4,1) * t218 + t401;
t151 = t174 * t246;
t322 = t174 * t332;
t54 = qJD(1) * t429 - t222 - t322;
t399 = t151 * t54;
t341 = rSges(4,2) * t369 + t246 * rSges(4,3);
t132 = rSges(4,1) * t367 - t341;
t321 = t174 * t331;
t291 = t221 - t321;
t53 = (-t132 + t354) * qJD(1) + t291;
t398 = t245 * t53;
t397 = t41 * t164;
t112 = Icges(5,5) * t371 - Icges(5,6) * t373 - t381;
t396 = -t245 * t112 - t116 * t370;
t395 = t245 * t113 + t117 * t370;
t375 = t168 * t245;
t374 = t168 * t246;
t363 = t246 * t112;
t56 = -t245 * t278 - t376;
t361 = t56 * qJD(1);
t61 = -t245 * t277 - t374;
t360 = t61 * qJD(1);
t359 = -t245 * t122 - t126 * t366;
t358 = t245 * t123 + t127 * t366;
t156 = qJD(1) * t194 - t222;
t357 = t204 - (-t246 * t406 - t223) * qJD(1) - t156;
t329 = qJD(1) * qJD(2);
t214 = qJ(2) * t333;
t338 = t214 + t221;
t350 = qJD(1) * (-pkin(1) * t334 + t338) + t245 * t329;
t347 = -t170 + t173;
t346 = t172 + t284;
t343 = rSges(4,3) * t333 + t334 * t431;
t205 = t245 * t402;
t340 = rSges(3,3) * t333 + qJD(1) * t205;
t339 = t246 * rSges(3,3) + t205;
t330 = t169 * qJD(1);
t328 = (qJD(3) ^ 2) * t409;
t327 = t118 * t333 + t245 * t74 + t246 * t73;
t326 = t245 * t405;
t324 = qJD(1) * (qJD(1) * t274 - t214) + t350;
t319 = -pkin(1) - t405;
t317 = t334 / 0.2e1;
t316 = t333 / 0.2e1;
t315 = -t332 / 0.2e1;
t312 = t331 / 0.2e1;
t311 = -t164 - t410;
t309 = -t211 - t404;
t94 = t117 * t371;
t308 = t246 * t113 - t94;
t307 = qJD(1) * t117 + t241 * t437;
t306 = (-t163 * t245 + t387) * qJD(1) + t436 * t241;
t305 = qJD(1) * t115 + t116 * t241 - t160 * t190;
t304 = (-t245 * t283 + t383) * qJD(1) + t435 * t241;
t99 = t127 * t367;
t303 = t246 * t123 - t99;
t302 = -t112 + t380;
t301 = -t122 + t378;
t300 = t433 * t241;
t299 = t434 * t241;
t134 = t165 * t241;
t292 = -qJD(3) * t409 - t134;
t289 = t404 - t431;
t288 = -t245 * t42 - t246 * t41;
t287 = -t245 * t54 - t246 * t53;
t58 = t114 * t210 + t116 * t209;
t65 = t124 * t219 + t126 * t218;
t66 = t125 * t219 + t127 * t218;
t150 = t174 * t245;
t50 = -t125 * t369 - t303;
t272 = (t245 * t50 - t246 * t49) * qJD(3);
t51 = -t124 * t368 - t359;
t52 = -t125 * t368 + t358;
t271 = (t245 * t52 - t246 * t51) * qJD(3);
t270 = t282 * t245;
t269 = qJD(3) * t172;
t268 = qJD(3) * t170;
t75 = (t132 * t245 + t133 * t246) * qJD(3);
t266 = qJD(1) * t159 - t190 * t376 + t191 * t377;
t264 = t124 * t246 - t125 * t245;
t255 = -t209 * t304 + t210 * t306 + t336;
t10 = t245 * t425 + t255 * t246;
t256 = qJD(1) * t112 - t209 * t305 + t210 * t307;
t11 = t256 * t245 - t246 * t424;
t12 = t255 * t245 - t246 * t425;
t43 = -t270 - t363;
t44 = -t115 * t373 - t308;
t20 = t190 * t44 - t191 * t43 + t361;
t45 = -t114 * t372 - t396;
t46 = -t115 * t372 + t395;
t57 = -t246 * t278 + t377;
t55 = t57 * qJD(1);
t21 = t190 * t46 - t191 * t45 + t55;
t259 = t434 * qJD(1) + t436 * t190 - t191 * t437;
t249 = -t209 * t419 + t259 * t210;
t254 = qJD(1) * t158 - t209 * t300 + t210 * t299;
t26 = t245 * t423 + t254 * t246;
t27 = t254 * t245 - t246 * t423;
t30 = t209 * t307 + t210 * t305;
t31 = t209 * t306 + t210 * t304;
t59 = t115 * t210 + t117 * t209;
t9 = t245 * t424 + t256 * t246;
t263 = (qJD(1) * t26 + t10 * t190 + t176 * t45 + t177 * t46 - t191 * t9) * t412 + (t259 * t209 + t210 * t419) * t408 + t20 * t317 + t21 * t316 + (qJD(1) * t27 - t11 * t191 + t12 * t190 + t176 * t43 + t177 * t44) * t411 + (t245 * t44 - t246 * t43) * t418 + (t245 * t46 - t246 * t45) * t417 + (t10 * t245 - t246 * t9 + (t245 * t45 + t246 * t46) * qJD(1)) * t415 + (-t11 * t246 + t12 * t245 + (t245 * t43 + t246 * t44) * qJD(1)) * t414 + (t245 * t31 - t246 * t30 + (t245 * t58 + t246 * t59) * qJD(1)) * t407 + (t245 * t266 + t246 * t249) * t416 + (t245 * t249 - t246 * t266) * t413;
t262 = (-t218 * t346 + t219 * t347) * qJD(1);
t83 = qJD(1) * t125 - t245 * t268;
t85 = qJD(1) * t127 - t245 * t269;
t253 = qJD(1) * t122 - qJD(3) * t65 - t218 * t83 + t219 * t85;
t82 = -t246 * t268 + (-t245 * t284 + t384) * qJD(1);
t84 = -t246 * t269 + (-t173 * t245 + t388) * qJD(1);
t252 = -qJD(3) * t66 - t218 * t82 + t219 * t84 + t335;
t154 = t284 * qJD(3);
t155 = t173 * qJD(3);
t251 = qJD(1) * t168 - t154 * t218 + t155 * t219 + (-t170 * t219 - t172 * t218) * qJD(3);
t250 = -t218 * t420 + t264 * t219;
t213 = t246 * t329;
t157 = t289 * qJD(3);
t143 = t326 - t339;
t93 = qJD(1) * t428 - t222;
t92 = t221 + (-t143 - t192) * qJD(1);
t87 = -qJD(3) * t150 + (t246 * t289 + t233) * qJD(1);
t86 = -t331 * t401 + (-t219 * t334 - t320) * rSges(4,1) + t343;
t79 = t213 + (-qJD(1) * t275 - t156) * qJD(1);
t78 = qJD(1) * (-qJD(1) * t326 + t340) + t350;
t62 = -t246 * t277 + t375;
t60 = t62 * qJD(1);
t40 = -t157 * t331 + t213 + (-t87 + t322 + t357) * qJD(1);
t39 = -t157 * t332 + (t86 - t321) * qJD(1) + t324;
t35 = -qJD(3) * t279 + t218 * t84 + t219 * t82;
t34 = -t280 * qJD(3) + t218 * t85 + t219 * t83;
t33 = t251 * t245 - t246 * t422;
t32 = t245 * t422 + t251 * t246;
t25 = t60 + t271;
t24 = t272 + t360;
t23 = -t246 * t328 - t134 * t191 + t164 * t176 + t213 + (-t74 - t77 + t186 + t357) * qJD(1);
t22 = -t245 * t328 - t134 * t190 - t164 * t177 + (t73 + t76 - t296) * qJD(1) + t324;
t1 = [(t60 + ((t50 - t99 + (t123 + t379) * t246 + t359) * t246 + t358 * t245) * qJD(3)) * t312 + (t55 + (t44 + (t114 * t246 + t115 * t245) * t209 + t308 + t396) * t191 + (-t116 * t371 + t363 + t43 + (t114 * t245 - t115 * t246) * t209 + t395) * t190) * t413 + (t58 + t56) * t418 + (t59 + t57) * t417 + (-t361 + (t46 - t270 - t395) * t191 + (t245 * t302 + t45 - t94) * t190 + ((t113 + t282) * t190 + t302 * t191) * t246 + t20) * t416 + (t31 + t26) * t415 + (t24 - t360 + ((t246 * t301 - t358 + t52) * t246 + (t245 * t301 + t303 + t51) * t245) * qJD(3)) * t315 + (t35 + t32) * t332 / 0.2e1 + (-qJD(3) * t277 + t154 * t219 + t155 * t218 + t209 * t299 + t210 * t300) * qJD(1) + (-(qJD(1) * t394 + t265 - t41 + t430) * t42 + t23 * (-t118 + t348) + t41 * t342 + t22 * (t119 + t167 - t365) + t42 * (t276 + t344) + (-t142 * t42 + t245 * t397) * t241 + ((t41 * (-t165 - t178) - t42 * t239) * t246 + (t41 * (-rSges(5,3) + t239) + t42 * (-t178 - t403)) * t245) * qJD(1)) * m(5) + (-(-qJD(1) * t132 + t291 + t430 - t53) * t54 + t40 * (t245 * t309 + t341 - t364) + t53 * (t204 + t222) + t39 * t429 + t54 * (t221 + t343) + (t174 * t398 - t399) * qJD(3) + ((-t53 * rSges(4,3) + t309 * t54) * t245 + (t53 * (-t211 - t289) - t54 * t244) * t246) * qJD(1)) * m(4) + (-(-qJD(1) * t143 - t181 + t221 - t92) * t93 + t79 * (t245 * t319 + t225 + t339) + t92 * t222 + t78 * t428 + t93 * (t338 + t340) + (t92 * (t319 + t402) * t246 + (t92 * (-rSges(3,3) - qJ(2)) + t93 * t319) * t245) * qJD(1)) * m(3) + (t30 + t27 + t21) * t414 - (t34 + t33 + t25) * t331 / 0.2e1 + ((t65 + t61) * t245 + (t66 + t62) * t246) * qJD(3) * t407; 0.2e1 * (t22 * t411 + t23 * t412) * m(5) + 0.2e1 * (t39 * t411 + t40 * t412) * m(4) + 0.2e1 * (t411 * t78 + t412 * t79) * m(3); t263 + ((-t331 * t375 - t330) * t246 + (t262 + (t246 * t374 + t250) * qJD(3)) * t245) * t312 + ((-t332 * t374 + t330) * t245 + (t262 + (t245 * t375 + t250) * qJD(3)) * t246) * t315 + ((t218 * t347 + t219 * t346) * qJD(1) + (t264 * t218 + t219 * t420) * qJD(3)) * t408 + (t245 * t35 - t246 * t34 + (t65 * t245 + t246 * t66) * qJD(1)) * t407 + (qJD(1) * t32 + (-(t245 * t426 + t253 * t246) * t246 + (t245 * t427 + t252 * t246) * t245 + (t51 * t245 + t52 * t246) * qJD(1)) * t432) * t412 + (qJD(1) * t33 + (-(t253 * t245 - t246 * t426) * t246 + (t252 * t245 - t246 * t427) * t245 + (t49 * t245 + t50 * t246) * qJD(1)) * t432) * t411 + (t272 + t24) * t317 + (t271 + t25) * t316 + (t36 * t327 + (t23 * t311 + t41 * t292 + t8 * t91 + t36 * t76 + (t311 * t42 - t36 * t90) * qJD(1)) * t246 + (t22 * t311 + t42 * t292 - t8 * t90 + t36 * t77 + (t36 * t393 + t397) * qJD(1)) * t245 - (-t42 * t218 * t333 + (t288 * t219 + t36 * (-t245 ^ 2 - t246 ^ 2) * t218) * qJD(3)) * pkin(3) + t421) * m(5) + (0.2e1 * t75 * (t245 * t87 + t246 * t86 + (t132 * t246 - t133 * t245) * qJD(1)) + t287 * t157 + (-t39 * t245 - t40 * t246 + (-t246 * t54 + t398) * qJD(1)) * t174 - (t150 * t53 - t399) * qJD(1) - (t75 * (-t150 * t245 - t151 * t246) + t287 * t289) * qJD(3)) * m(4); t263 + (t36 * (-t119 * t334 + t327) + t288 * t134 + (-t22 * t245 - t23 * t246 + (t245 * t41 - t246 * t42) * qJD(1)) * t164 + t421) * m(5);];
tauc = t1(:);
