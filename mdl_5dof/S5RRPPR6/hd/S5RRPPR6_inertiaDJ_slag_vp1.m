% Calculate time derivative of joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:37
% EndTime: 2019-12-31 19:32:00
% DurationCPUTime: 12.92s
% Computational Cost: add. (15866->770), mult. (18158->1092), div. (0->0), fcn. (16827->10), ass. (0->379)
t260 = sin(qJ(1));
t445 = t260 / 0.2e1;
t262 = cos(qJ(1));
t466 = -t262 / 0.2e1;
t465 = -qJD(1) / 0.2e1;
t256 = cos(pkin(9));
t238 = pkin(4) * t256 + pkin(3);
t252 = qJ(2) + pkin(8);
t245 = cos(t252);
t401 = t245 * t262;
t255 = sin(pkin(9));
t395 = t260 * t255;
t236 = pkin(4) * t395;
t258 = -pkin(7) - qJ(4);
t243 = sin(t252);
t405 = t243 * t262;
t460 = -t258 * t405 + t236;
t251 = pkin(9) + qJ(5);
t242 = sin(t251);
t244 = cos(t251);
t396 = t260 * t244;
t181 = -t242 * t401 + t396;
t397 = t260 * t242;
t182 = t244 * t401 + t397;
t98 = t182 * rSges(6,1) + t181 * rSges(6,2) + rSges(6,3) * t405;
t464 = t238 * t401 + t460 + t98;
t435 = pkin(3) - t238;
t463 = t243 * t435;
t393 = qJ(4) + t258;
t462 = t245 * t393;
t261 = cos(qJ(2));
t239 = pkin(2) * t261 + pkin(1);
t436 = pkin(1) - t239;
t461 = t260 * t436;
t259 = sin(qJ(2));
t423 = Icges(3,4) * t261;
t313 = -Icges(3,2) * t259 + t423;
t188 = Icges(3,6) * t260 + t262 * t313;
t424 = Icges(3,4) * t259;
t319 = Icges(3,1) * t261 - t424;
t190 = Icges(3,5) * t260 + t262 * t319;
t296 = t188 * t259 - t190 * t261;
t278 = t296 * t260;
t187 = -Icges(3,6) * t262 + t260 * t313;
t189 = -Icges(3,5) * t262 + t260 * t319;
t297 = t187 * t259 - t189 * t261;
t279 = t297 * t262;
t421 = Icges(4,4) * t245;
t311 = -Icges(4,2) * t243 + t421;
t168 = Icges(4,6) * t260 + t262 * t311;
t422 = Icges(4,4) * t243;
t317 = Icges(4,1) * t245 - t422;
t170 = Icges(4,5) * t260 + t262 * t317;
t298 = t168 * t243 - t170 * t245;
t280 = t298 * t260;
t167 = -Icges(4,6) * t262 + t260 * t311;
t169 = -Icges(4,5) * t262 + t260 * t317;
t299 = t167 * t243 - t169 * t245;
t281 = t299 * t262;
t248 = t260 * rSges(4,3);
t459 = -rSges(4,2) * t405 + t248;
t306 = Icges(4,5) * t245 - Icges(4,6) * t243;
t165 = -Icges(4,3) * t262 + t260 * t306;
t307 = Icges(3,5) * t261 - Icges(3,6) * t259;
t185 = -Icges(3,3) * t262 + t260 * t307;
t342 = qJD(1) * t245 - qJD(5);
t374 = qJD(2) * t262;
t351 = t243 * t374;
t458 = t260 * t342 + t351;
t376 = qJD(2) * t260;
t352 = t243 * t376;
t457 = t262 * t342 - t352;
t456 = 2 * m(3);
t455 = 2 * m(4);
t454 = 2 * m(5);
t453 = 2 * m(6);
t253 = t260 ^ 2;
t254 = t262 ^ 2;
t452 = m(5) / 0.2e1;
t451 = m(6) / 0.2e1;
t309 = Icges(5,4) * t256 - Icges(5,2) * t255;
t158 = -Icges(5,6) * t245 + t243 * t309;
t450 = t158 / 0.2e1;
t315 = Icges(5,1) * t256 - Icges(5,4) * t255;
t159 = -Icges(5,5) * t245 + t243 * t315;
t449 = t159 / 0.2e1;
t448 = -t245 / 0.2e1;
t447 = -t255 / 0.2e1;
t446 = t256 / 0.2e1;
t443 = t262 / 0.2e1;
t224 = rSges(3,1) * t259 + rSges(3,2) * t261;
t442 = m(3) * t224;
t441 = pkin(2) * t259;
t440 = pkin(3) * t243;
t439 = pkin(3) * t245;
t438 = t260 * pkin(6);
t250 = t262 * pkin(6);
t437 = qJD(1) / 0.2e1;
t257 = -qJ(3) - pkin(6);
t434 = -pkin(6) - t257;
t433 = rSges(3,1) * t261;
t432 = rSges(4,1) * t245;
t431 = rSges(3,2) * t259;
t430 = rSges(3,3) * t262;
t429 = rSges(6,3) * t243;
t249 = t260 * rSges(3,3);
t428 = -rSges(5,3) - qJ(4);
t427 = -rSges(6,3) + t258;
t266 = -t243 * t393 - t245 * t435;
t400 = t255 * t262;
t371 = pkin(4) * t400;
t179 = -t244 * t262 - t245 * t397;
t180 = -t242 * t262 + t245 * t396;
t332 = -t180 * rSges(6,1) - t179 * rSges(6,2);
t406 = t243 * t260;
t97 = rSges(6,3) * t406 - t332;
t426 = t260 * t266 - t371 + t97;
t229 = pkin(3) * t401;
t192 = qJ(4) * t405 + t229;
t425 = -t192 + t464;
t420 = Icges(6,4) * t242;
t419 = Icges(6,4) * t244;
t314 = Icges(6,1) * t244 - t420;
t149 = -Icges(6,5) * t245 + t243 * t314;
t412 = t149 * t244;
t411 = t187 * t261;
t410 = t188 * t261;
t409 = t189 * t259;
t408 = t190 * t259;
t407 = t243 * t238;
t404 = t245 * t255;
t403 = t245 * t256;
t402 = t245 * t258;
t399 = t256 * t262;
t398 = t257 * t262;
t394 = t260 * t256;
t331 = rSges(6,1) * t244 - rSges(6,2) * t242;
t151 = -rSges(6,3) * t245 + t243 * t331;
t392 = t151 + t462 - t463;
t163 = t250 + t398 - t461;
t230 = t262 * t239;
t164 = -pkin(1) * t262 + t260 * t434 + t230;
t391 = t260 * t163 + t262 * t164;
t390 = -t164 - t192;
t208 = -qJ(4) * t245 + t440;
t381 = qJD(1) * t260;
t354 = t259 * t381;
t234 = pkin(2) * t354;
t389 = t208 * t381 + t234;
t355 = t243 * t381;
t388 = qJD(1) * t371 + t258 * t355;
t380 = qJD(1) * t262;
t387 = rSges(4,2) * t355 + rSges(4,3) * t380;
t373 = qJD(4) * t243;
t223 = t262 * t373;
t246 = qJD(3) * t260;
t386 = t223 + t246;
t385 = t262 * t433 + t249;
t384 = t253 + t254;
t166 = Icges(4,3) * t260 + t262 * t306;
t383 = qJD(1) * t166;
t186 = Icges(3,3) * t260 + t262 * t307;
t382 = qJD(1) * t186;
t379 = qJD(2) * t243;
t378 = qJD(2) * t245;
t377 = qJD(2) * t259;
t375 = qJD(2) * t261;
t372 = qJD(5) * t243;
t349 = t245 * t374;
t343 = -qJD(5) * t245 + qJD(1);
t293 = t262 * t343;
t83 = t242 * t458 + t244 * t293;
t84 = t242 * t293 - t244 * t458;
t370 = t84 * rSges(6,1) + t83 * rSges(6,2) + rSges(6,3) * t349;
t369 = t262 * t431;
t367 = pkin(2) * t377;
t366 = pkin(2) * t375;
t93 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t406;
t95 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t406;
t329 = -t242 * t93 + t244 * t95;
t91 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t406;
t33 = t243 * t329 - t245 * t91;
t304 = Icges(6,5) * t244 - Icges(6,6) * t242;
t147 = -Icges(6,3) * t245 + t243 * t304;
t308 = -Icges(6,2) * t242 + t419;
t148 = -Icges(6,6) * t245 + t243 * t308;
t49 = t147 * t406 + t148 * t179 + t149 * t180;
t365 = t49 / 0.2e1 + t33 / 0.2e1;
t94 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t405;
t96 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t405;
t328 = -t242 * t94 + t244 * t96;
t92 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t405;
t34 = t243 * t328 - t245 * t92;
t50 = t147 * t405 + t181 * t148 + t182 * t149;
t364 = t50 / 0.2e1 + t34 / 0.2e1;
t197 = t245 * t394 - t400;
t285 = t245 * t395 + t399;
t109 = Icges(5,5) * t197 - Icges(5,6) * t285 + Icges(5,3) * t406;
t363 = t109 * t406;
t362 = t109 * t405;
t198 = -t245 * t400 + t394;
t199 = t245 * t399 + t395;
t110 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t405;
t361 = t110 * t406;
t360 = t110 * t405;
t356 = qJD(3) * t262 + t257 * t381 + t260 * t367;
t358 = t260 * ((-t262 * t436 - t438) * qJD(1) - t356) + t262 * (-t262 * t367 + t246 + (t262 * t434 + t461) * qJD(1)) + t163 * t380;
t140 = qJD(1) * t285 + t255 * t351;
t141 = -qJD(1) * t197 - t256 * t351;
t357 = t141 * rSges(5,1) + t140 * rSges(5,2) + rSges(5,3) * t349;
t122 = t199 * rSges(5,1) + t198 * rSges(5,2) + rSges(5,3) * t405;
t353 = t148 * t378;
t350 = t245 * t376;
t348 = t378 / 0.2e1;
t347 = -t208 - t441;
t209 = rSges(4,1) * t243 + rSges(4,2) * t245;
t346 = -t209 - t441;
t345 = -t260 * t257 + t230;
t344 = -t238 * t245 - t239;
t330 = qJ(4) * t243 + t439;
t191 = t330 * t260;
t341 = t260 * t191 + t262 * t192 + t391;
t333 = rSges(5,1) * t256 - rSges(5,2) * t255;
t160 = -rSges(5,3) * t245 + t243 * t333;
t340 = -t160 + t347;
t339 = -qJD(2) * t330 + qJD(4) * t245 - t366;
t294 = t260 * t343;
t85 = -t242 * t457 + t244 * t294;
t86 = t242 * t294 + t244 * t457;
t338 = t86 * rSges(6,1) + t85 * rSges(6,2);
t337 = -t431 + t433;
t336 = -rSges(4,2) * t243 + t432;
t142 = qJD(1) * t198 + t255 * t352;
t143 = qJD(1) * t199 - t256 * t352;
t335 = -t143 * rSges(5,1) - t142 * rSges(5,2);
t334 = -t197 * rSges(5,1) + rSges(5,2) * t285;
t28 = t179 * t93 + t180 * t95 + t406 * t91;
t29 = t179 * t94 + t180 * t96 + t406 * t92;
t18 = t29 * t260 - t262 * t28;
t327 = t260 * t28 + t262 * t29;
t30 = t181 * t93 + t182 * t95 + t405 * t91;
t31 = t181 * t94 + t182 * t96 + t405 * t92;
t19 = t31 * t260 - t262 * t30;
t326 = t260 * t30 + t262 * t31;
t325 = t34 * t260 - t33 * t262;
t324 = t260 * t33 + t262 * t34;
t268 = t243 * t427 + t344;
t56 = (pkin(4) * t255 - t257) * t262 + t268 * t260 + t332;
t57 = t345 + t464;
t323 = t260 * t57 + t262 * t56;
t58 = t151 * t406 + t245 * t97;
t59 = -t151 * t405 - t245 * t98;
t322 = t260 * t59 + t262 * t58;
t271 = t243 * t428 - t239 - t439;
t263 = t260 * t271 - t398;
t72 = t263 + t334;
t73 = t345 + t122 + t192;
t321 = t260 * t73 + t262 * t72;
t320 = -t260 * t98 + t262 * t97;
t318 = Icges(3,1) * t259 + t423;
t316 = Icges(4,1) * t243 + t421;
t312 = Icges(3,2) * t261 + t424;
t310 = Icges(4,2) * t245 + t422;
t305 = Icges(5,5) * t256 - Icges(5,6) * t255;
t295 = t347 - t392;
t176 = rSges(4,1) * t401 + t459;
t292 = -(rSges(5,3) * t243 + t245 * t333) * qJD(2) + t339;
t291 = -pkin(1) - t337;
t104 = t340 * t262;
t215 = qJ(4) * t349;
t219 = pkin(3) * t352;
t270 = t243 * t380 + t350;
t289 = t262 * (-qJ(4) * t355 + t215 + t223 + (-t245 * t381 - t351) * pkin(3)) + t191 * t380 + t260 * (qJ(4) * t270 + qJD(1) * t229 + t260 * t373 - t219) + t358;
t288 = -t239 - t336;
t90 = (-rSges(6,1) * t242 - rSges(6,2) * t244) * t372 + (t245 * t331 + t429) * qJD(2);
t284 = -t266 * qJD(2) + t339 - t90;
t283 = qJD(2) * t224;
t282 = qJD(2) * t209;
t277 = qJD(2) * t318;
t276 = qJD(2) * t316;
t275 = qJD(2) * t312;
t274 = qJD(2) * t310;
t273 = qJD(2) * (-Icges(3,5) * t259 - Icges(3,6) * t261);
t272 = qJD(2) * (-Icges(4,5) * t243 - Icges(4,6) * t245);
t61 = t295 * t262;
t269 = t349 - t355;
t87 = (-Icges(6,5) * t242 - Icges(6,6) * t244) * t372 + (Icges(6,3) * t243 + t245 * t304) * qJD(2);
t89 = (-Icges(6,1) * t242 - t419) * t372 + (Icges(6,5) * t243 + t245 * t314) * qJD(2);
t265 = t147 * t379 - t245 * t87 + t378 * t412 + (-t148 * t372 + t243 * t89) * t244;
t264 = rSges(3,2) * t354 + rSges(3,3) * t380 - t262 * t283;
t214 = t337 * qJD(2);
t203 = t336 * qJD(2);
t194 = -t369 + t385;
t193 = t260 * t337 - t430;
t175 = -rSges(4,3) * t262 + t260 * t336;
t162 = t346 * t262;
t161 = t346 * t260;
t156 = t438 + (pkin(1) - t431) * t262 + t385;
t155 = t260 * t291 + t250 + t430;
t146 = (Icges(5,5) * t243 + t245 * t315) * qJD(2);
t145 = (Icges(5,6) * t243 + t245 * t309) * qJD(2);
t135 = t176 + t345;
t134 = (rSges(4,3) - t257) * t262 + t288 * t260;
t128 = t260 * t273 + t382;
t127 = -qJD(1) * t185 + t262 * t273;
t121 = rSges(5,3) * t406 - t334;
t116 = t260 * t272 + t383;
t115 = -qJD(1) * t165 + t262 * t272;
t114 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t405;
t113 = Icges(5,1) * t197 - Icges(5,4) * t285 + Icges(5,5) * t406;
t112 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t405;
t111 = Icges(5,4) * t197 - Icges(5,2) * t285 + Icges(5,6) * t406;
t108 = t224 * t376 + ((-rSges(3,3) - pkin(6)) * t260 + t291 * t262) * qJD(1);
t107 = (t250 + (-pkin(1) - t433) * t260) * qJD(1) + t264;
t103 = t340 * t260;
t102 = -t209 * t380 - t260 * t203 + (-t259 * t380 - t260 * t375) * pkin(2);
t101 = t209 * t381 + t234 + (-t203 - t366) * t262;
t88 = (-Icges(6,2) * t244 - t420) * t372 + (Icges(6,6) * t243 + t245 * t308) * qJD(2);
t79 = t260 * t186 - t296 * t262;
t78 = t260 * t185 - t279;
t77 = -t186 * t262 - t278;
t76 = -t185 * t262 - t260 * t297;
t75 = t209 * t376 + (t262 * t288 - t248) * qJD(1) + t356;
t74 = t246 + (-t239 - t432) * t381 + (-qJD(1) * t257 + qJD(2) * t346) * t262 + t387;
t71 = t260 * t166 - t298 * t262;
t70 = t260 * t165 - t281;
t69 = -t166 * t262 - t280;
t68 = -t165 * t262 - t260 * t299;
t67 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t270;
t66 = Icges(5,1) * t141 + Icges(5,4) * t140 + Icges(5,5) * t269;
t65 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t270;
t64 = Icges(5,4) * t141 + Icges(5,2) * t140 + Icges(5,6) * t269;
t60 = t295 * t260;
t55 = -t147 * t245 + (-t148 * t242 + t412) * t243;
t54 = t55 * t379;
t53 = qJD(1) * t104 + t260 * t292;
t52 = t160 * t381 + t262 * t292 + t389;
t51 = t320 * t243;
t48 = rSges(6,3) * t270 + t338;
t47 = -rSges(6,3) * t355 + t370;
t46 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t270;
t45 = Icges(6,1) * t84 + Icges(6,4) * t83 + Icges(6,5) * t269;
t44 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t270;
t43 = Icges(6,4) * t84 + Icges(6,2) * t83 + Icges(6,6) * t269;
t42 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t270;
t41 = Icges(6,5) * t84 + Icges(6,6) * t83 + Icges(6,3) * t269;
t40 = t219 + (t378 * t428 - t373) * t260 + t271 * t380 + t335 + t356;
t39 = t215 + (-t440 - t441) * t374 + t263 * qJD(1) + t357 + t386;
t38 = t198 * t112 + t199 * t114 + t360;
t37 = t198 * t111 + t199 * t113 + t362;
t36 = -t112 * t285 + t114 * t197 + t361;
t35 = -t111 * t285 + t113 * t197 + t363;
t32 = t260 * t121 + t122 * t262 + t341;
t27 = qJD(1) * t61 + t260 * t284;
t26 = t262 * t284 + t381 * t392 + t389;
t25 = (-t373 + (t245 * t427 + t407) * qJD(2)) * t260 + (t262 * t268 - t236) * qJD(1) - t338 + t356;
t24 = (-t402 - t407 - t441) * t374 + (-t398 + (t344 - t429) * t260) * qJD(1) + t370 + t386 + t388;
t23 = t260 * t426 + t262 * t425 + t341;
t22 = (t151 * t376 + t48) * t245 + (-qJD(2) * t97 + t151 * t380 + t260 * t90) * t243;
t21 = (-t151 * t374 - t47) * t245 + (qJD(2) * t98 + t151 * t381 - t262 * t90) * t243;
t20 = (-t353 + (-qJD(5) * t149 - t88) * t243) * t242 + t265;
t17 = t147 * t270 + t85 * t148 + t86 * t149 + t179 * t88 + t180 * t89 + t406 * t87;
t16 = t147 * t269 + t83 * t148 + t84 * t149 + t181 * t88 + t182 * t89 + t405 * t87;
t15 = t320 * t378 + (-t260 * t47 + t262 * t48 + (-t260 * t97 - t262 * t98) * qJD(1)) * t243;
t14 = t260 * (rSges(5,3) * t350 - t335) + t262 * t357 + (t262 * t121 + (-t122 + t390) * t260) * qJD(1) + t289;
t13 = t243 * t326 - t50 * t245;
t12 = t243 * t327 - t49 * t245;
t11 = (qJD(2) * t328 - t41) * t245 + (qJD(2) * t92 - t242 * t43 + t244 * t45 + (-t242 * t96 - t244 * t94) * qJD(5)) * t243;
t10 = (qJD(2) * t329 - t42) * t245 + (qJD(2) * t91 - t242 * t44 + t244 * t46 + (-t242 * t95 - t244 * t93) * qJD(5)) * t243;
t9 = t92 * t350 + t179 * t43 + t180 * t45 + t85 * t94 + t86 * t96 + (t260 * t41 + t380 * t92) * t243;
t8 = t91 * t350 + t179 * t44 + t180 * t46 + t85 * t93 + t86 * t95 + (t260 * t42 + t380 * t91) * t243;
t7 = t92 * t349 + t181 * t43 + t182 * t45 + t83 * t94 + t84 * t96 + (t262 * t41 - t381 * t92) * t243;
t6 = t91 * t349 + t181 * t44 + t182 * t46 + t83 * t93 + t84 * t95 + (t262 * t42 - t381 * t91) * t243;
t5 = (-t215 + t47 + t388) * t262 + (t219 + t48) * t260 + (t254 * (-t402 + t463) + (-t407 - t462) * t253) * qJD(2) + (t426 * t262 + (t390 - t425 + t460) * t260) * qJD(1) + t289;
t4 = qJD(1) * t327 + t9 * t260 - t262 * t8;
t3 = qJD(1) * t326 + t7 * t260 - t262 * t6;
t2 = (qJD(2) * t327 - t17) * t245 + (-qJD(1) * t18 + qJD(2) * t49 + t260 * t8 + t262 * t9) * t243;
t1 = (qJD(2) * t326 - t16) * t245 + (-qJD(1) * t19 + qJD(2) * t50 + t260 * t6 + t262 * t7) * t243;
t62 = [(t107 * t156 + t108 * t155) * t456 + (t134 * t75 + t135 * t74) * t455 + (t39 * t73 + t40 * t72) * t454 + (t24 * t57 + t25 * t56) * t453 + t265 + (-t312 + t319) * t377 + (t318 + t313) * t375 + (-t255 * t145 + t256 * t146) * t243 + (-Icges(5,3) * t245 + t243 * t305 - t310 + t317) * t379 + (-t149 * t372 - t243 * t88 - t353) * t242 + (-Icges(5,3) * t243 - t255 * t158 + t256 * t159 - t245 * t305 + t311 + t316) * t378; m(4) * (t101 * t134 + t102 * t135 + t161 * t74 + t162 * t75) + m(5) * (t103 * t39 + t104 * t40 + t52 * t72 + t53 * t73) + m(6) * (t24 * t60 + t25 * t61 + t26 * t56 + t27 * t57) + m(3) * ((-t107 * t260 - t108 * t262) * t224 + (-t155 * t262 - t156 * t260) * t214) + ((Icges(5,5) * t143 / 0.2e1 + Icges(5,6) * t142 / 0.2e1 + Icges(5,3) * t270 / 0.2e1 + t168 * t465 + t274 * t445) * t262 + (-Icges(5,5) * t141 / 0.2e1 - Icges(5,6) * t140 / 0.2e1 - Icges(5,3) * t269 / 0.2e1 + t167 * t465 + t274 * t466) * t260) * t245 + ((t198 * t450 + t199 * t449 - t156 * t442 + t410 / 0.2e1 + t408 / 0.2e1 + (t168 / 0.2e1 - t110 / 0.2e1) * t245 + (t170 / 0.2e1 + t112 * t447 + t114 * t446) * t243 + t364) * t262 + (-t285 * t450 + t197 * t449 + t411 / 0.2e1 + t409 / 0.2e1 + t155 * t442 + (t167 / 0.2e1 - t109 / 0.2e1) * t245 + (t169 / 0.2e1 + t111 * t447 + t113 * t446) * t243 + t365) * t260) * qJD(1) + (t140 * t158 + t141 * t159 + t198 * t145 + t199 * t146 + t16 + (-qJD(1) * t187 - t262 * t275) * t261 + (-qJD(1) * t189 - t262 * t277) * t259 + t11 + (-qJD(1) * t169 - t255 * t64 + t256 * t66 - t262 * t276) * t243) * t445 + (t142 * t158 + t143 * t159 - t285 * t145 + t197 * t146 + (qJD(1) * t188 - t260 * t275) * t261 + (qJD(1) * t190 - t260 * t277) * t259 + t17 + t10 + (qJD(1) * t170 - t255 * t65 + t256 * t67 - t260 * t276) * t243) * t466 + (t279 / 0.2e1 + t281 / 0.2e1 - t280 / 0.2e1 - t278 / 0.2e1 + (t306 + t307) * (t254 / 0.2e1 + t253 / 0.2e1) + (t110 * t243 - t112 * t404 + t114 * t403) * t445 + (t109 * t243 - t111 * t404 + t113 * t403) * t466) * qJD(2); (t23 * t5 + t26 * t61 + t27 * t60) * t453 - t262 * t4 + t260 * t3 + (t103 * t53 + t104 * t52 + t14 * t32) * t454 + (t162 * t101 + t161 * t102 + (t260 * t175 + t176 * t262 + t391) * ((qJD(1) * t175 - t262 * t282 + t387) * t262 + (-t260 * t282 + (-t164 - t176 + t459) * qJD(1)) * t260 + t358)) * t455 + ((t260 * t193 + t194 * t262) * ((qJD(1) * t193 + t264) * t262 + (-t260 * t283 + (-t194 - t369 + t249) * qJD(1)) * t260) + t384 * t224 * t214) * t456 + t260 * ((t140 * t112 + t141 * t114 + t198 * t64 + t199 * t66 + (t37 - t361) * qJD(1)) * t260 + (-t140 * t111 - t141 * t113 - t198 * t65 - t199 * t67 + (t38 + t363) * qJD(1)) * t262) + t260 * ((t260 * t115 + (t70 + t280) * qJD(1)) * t260 + (t71 * qJD(1) + (t167 * t378 + t169 * t379) * t262 + (-t116 + (-t168 * t245 - t170 * t243) * qJD(2) + (t166 - t299) * qJD(1)) * t260) * t262) - t262 * ((t262 * t116 + (t69 + t281) * qJD(1)) * t262 + (t68 * qJD(1) + (-t168 * t378 - t170 * t379 + t383) * t260 + (-t115 + (t167 * t245 + t169 * t243) * qJD(2) - t298 * qJD(1)) * t262) * t260) - t262 * ((t262 * t128 + (t77 + t279) * qJD(1)) * t262 + (t76 * qJD(1) + (-t188 * t375 - t190 * t377 + t382) * t260 + (-t127 + (t409 + t411) * qJD(2) - t296 * qJD(1)) * t262) * t260) + t260 * ((t260 * t127 + (t78 + t278) * qJD(1)) * t260 + (t79 * qJD(1) + (t187 * t375 + t189 * t377) * t262 + (-t128 + (-t408 - t410) * qJD(2) + (t186 - t297) * qJD(1)) * t260) * t262) - t262 * ((-t142 * t111 - t143 * t113 + t285 * t65 - t197 * t67 + (t36 - t362) * qJD(1)) * t262 + (t142 * t112 + t143 * t114 - t285 * t64 + t197 * t66 + (t35 + t360) * qJD(1)) * t260) + (t18 + (-t35 - t68 - t76) * t262 + (t36 + t69 + t77) * t260) * t381 + (t19 + (-t37 - t70 - t78) * t262 + (t38 + t71 + t79) * t260) * t380; m(6) * (qJD(1) * t323 - t24 * t262 + t260 * t25) + m(5) * (qJD(1) * t321 + t260 * t40 - t262 * t39) + m(4) * (t260 * t75 - t262 * t74 + (t134 * t262 + t135 * t260) * qJD(1)); m(6) * (t260 * t26 - t262 * t27 + (t260 * t60 + t262 * t61) * qJD(1)) + m(5) * (t260 * t52 - t262 * t53 + (t103 * t260 + t104 * t262) * qJD(1)) + m(4) * (t260 * t101 - t102 * t262 + (t161 * t260 + t162 * t262) * qJD(1)); 0; 0.2e1 * (t321 * t452 + t323 * t451) * t378 + 0.2e1 * ((t24 * t260 + t25 * t262 + t380 * t57 - t381 * t56) * t451 + (t260 * t39 + t262 * t40 + t380 * t73 - t381 * t72) * t452) * t243; 0.2e1 * ((t374 * t61 + t376 * t60 - t5) * t451 + (t103 * t376 + t104 * t374 - t14) * t452) * t245 + 0.2e1 * ((qJD(2) * t23 + t26 * t262 + t260 * t27 + t380 * t60 - t381 * t61) * t451 + (qJD(2) * t32 + t103 * t380 - t104 * t381 + t260 * t53 + t262 * t52) * t452) * t243; 0; 0.4e1 * (t452 + t451) * (-0.1e1 + t384) * t243 * t378; m(6) * (t21 * t57 + t22 * t56 + t24 * t59 + t25 * t58) + t54 + (-t20 + (t260 * t365 + t262 * t364) * qJD(2)) * t245 + ((t16 / 0.2e1 + t11 / 0.2e1) * t262 + (t17 / 0.2e1 + t10 / 0.2e1) * t260 + (-t260 * t364 + t262 * t365) * qJD(1)) * t243; m(6) * (t15 * t23 + t21 * t60 + t22 * t61 + t26 * t58 + t27 * t59 + t5 * t51) + (-t2 / 0.2e1 + t19 * t348 + (t34 * qJD(1) - t10) * t448 + t13 * t437) * t262 + (t12 * t437 + (qJD(1) * t33 + t11) * t448 + t18 * t348 + t1 / 0.2e1) * t260 + (t3 * t443 + qJD(2) * t325 / 0.2e1 + t4 * t445 + (-t260 * t19 / 0.2e1 + t18 * t443) * qJD(1)) * t243; m(6) * (qJD(1) * t322 - t21 * t262 + t22 * t260); m(6) * ((qJD(2) * t322 - t15) * t245 + (qJD(2) * t51 + t21 * t260 + t22 * t262 + (-t260 * t58 + t262 * t59) * qJD(1)) * t243); (t15 * t51 + t21 * t59 + t22 * t58) * t453 + (t20 * t245 - t54 + (t260 * t12 + t262 * t13 - t245 * t324) * qJD(2)) * t245 + (t262 * t1 + t260 * t2 - t245 * (t10 * t260 + t11 * t262) + (t243 * t324 - t55 * t245) * qJD(2) + (t262 * t12 - t260 * t13 + t245 * t325) * qJD(1)) * t243;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t62(1), t62(2), t62(4), t62(7), t62(11); t62(2), t62(3), t62(5), t62(8), t62(12); t62(4), t62(5), t62(6), t62(9), t62(13); t62(7), t62(8), t62(9), t62(10), t62(14); t62(11), t62(12), t62(13), t62(14), t62(15);];
Mq = res;
