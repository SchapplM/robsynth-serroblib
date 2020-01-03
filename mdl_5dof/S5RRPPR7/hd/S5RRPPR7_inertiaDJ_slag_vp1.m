% Calculate time derivative of joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:58
% EndTime: 2019-12-31 19:35:31
% DurationCPUTime: 17.67s
% Computational Cost: add. (11140->689), mult. (16788->976), div. (0->0), fcn. (15327->8), ass. (0->352)
t238 = qJ(2) + pkin(8);
t230 = cos(t238);
t229 = sin(t238);
t404 = Icges(5,6) * t229;
t414 = Icges(4,4) * t229;
t468 = -t404 - t414 + (-Icges(4,2) - Icges(5,3)) * t230;
t403 = Icges(5,6) * t230;
t413 = Icges(4,4) * t230;
t467 = -t403 - t413 + (-Icges(4,1) - Icges(5,2)) * t229;
t466 = t468 * qJD(2);
t465 = t467 * qJD(2);
t244 = sin(qJ(1));
t247 = cos(qJ(1));
t297 = -Icges(4,2) * t229 + t413;
t143 = -Icges(4,6) * t247 + t244 * t297;
t288 = -Icges(5,3) * t229 + t403;
t449 = Icges(5,5) * t247 + t244 * t288;
t464 = -t143 - t449;
t144 = Icges(4,6) * t244 + t247 * t297;
t147 = Icges(5,5) * t244 - t247 * t288;
t463 = t144 - t147;
t302 = Icges(4,1) * t230 - t414;
t146 = Icges(4,5) * t244 + t247 * t302;
t290 = Icges(5,2) * t230 - t404;
t149 = Icges(5,4) * t244 - t247 * t290;
t462 = t146 - t149;
t145 = -Icges(4,5) * t247 + t244 * t302;
t448 = Icges(5,4) * t247 + t244 * t290;
t461 = -t448 - t145;
t361 = qJD(2) * t247;
t340 = t229 * t361;
t369 = qJD(1) * t244;
t460 = t230 * t369 + t340;
t459 = qJD(2) / 0.2e1;
t246 = cos(qJ(2));
t243 = sin(qJ(2));
t416 = Icges(3,4) * t243;
t304 = Icges(3,1) * t246 - t416;
t166 = Icges(3,5) * t244 + t247 * t304;
t390 = t166 * t246;
t415 = Icges(3,4) * t246;
t299 = -Icges(3,2) * t243 + t415;
t164 = Icges(3,6) * t244 + t247 * t299;
t395 = t164 * t243;
t276 = -t390 + t395;
t458 = t244 * t276;
t279 = t147 * t229 - t149 * t230;
t457 = t244 * t279;
t280 = t144 * t229 - t146 * t230;
t456 = t244 * t280;
t433 = pkin(2) * t246;
t224 = pkin(1) + t433;
t428 = pkin(1) - t224;
t455 = t244 * t428;
t165 = -Icges(3,5) * t247 + t244 * t304;
t392 = t165 * t246;
t163 = -Icges(3,6) * t247 + t244 * t299;
t397 = t163 * t243;
t277 = -t392 + t397;
t454 = t247 * t277;
t278 = -t229 * t449 + t230 * t448;
t453 = t247 * t278;
t281 = t143 * t229 - t145 * t230;
t452 = t247 * t281;
t233 = t244 * rSges(4,3);
t389 = t229 * t247;
t451 = -rSges(4,2) * t389 + t233;
t235 = t244 * rSges(5,1);
t387 = t230 * t247;
t450 = -rSges(5,2) * t387 + t235;
t292 = Icges(4,5) * t230 - Icges(4,6) * t229;
t141 = -Icges(4,3) * t247 + t244 * t292;
t293 = Icges(3,5) * t246 - Icges(3,6) * t243;
t161 = -Icges(3,3) * t247 + t244 * t293;
t295 = Icges(5,4) * t230 - Icges(5,5) * t229;
t447 = Icges(5,1) * t247 + t244 * t295;
t446 = 2 * m(3);
t445 = 2 * m(4);
t444 = 2 * m(5);
t443 = 2 * m(6);
t442 = m(5) / 0.2e1;
t441 = m(6) / 0.2e1;
t440 = t229 / 0.2e1;
t439 = t244 / 0.2e1;
t438 = -t247 / 0.2e1;
t437 = t247 / 0.2e1;
t436 = rSges(6,3) + pkin(7);
t211 = rSges(3,1) * t243 + rSges(3,2) * t246;
t435 = m(3) * t211;
t434 = pkin(2) * t243;
t432 = pkin(3) * t229;
t431 = pkin(3) * t230;
t236 = t244 * pkin(4);
t430 = t244 * pkin(6);
t237 = t247 * pkin(6);
t429 = qJD(1) / 0.2e1;
t241 = -qJ(3) - pkin(6);
t427 = -pkin(6) - t241;
t242 = sin(qJ(5));
t245 = cos(qJ(5));
t410 = Icges(6,4) * t245;
t300 = Icges(6,1) * t242 + t410;
t135 = Icges(6,5) * t229 - t230 * t300;
t291 = Icges(6,5) * t242 + Icges(6,6) * t245;
t133 = Icges(6,3) * t229 - t230 * t291;
t411 = Icges(6,4) * t242;
t294 = Icges(6,2) * t245 + t411;
t134 = Icges(6,6) * t229 - t230 * t294;
t359 = qJD(5) * t230;
t363 = qJD(2) * t245;
t366 = qJD(2) * t230;
t367 = qJD(2) * t229;
t398 = t135 * t242;
t81 = (-Icges(6,5) * t245 + Icges(6,6) * t242) * t359 + (Icges(6,3) * t230 + t229 * t291) * qJD(2);
t324 = t133 * t366 + t229 * t81 + t367 * t398 + (t229 * t363 + t242 * t359) * t134;
t83 = (-Icges(6,1) * t245 + t411) * t359 + (Icges(6,5) * t230 + t229 * t300) * qJD(2);
t419 = t242 * t83;
t49 = t133 * t229 + (-t134 * t245 - t398) * t230;
t82 = (Icges(6,2) * t242 - t410) * t359 + (Icges(6,6) * t230 + t229 * t294) * qJD(2);
t426 = ((-t419 + (-qJD(5) * t135 - t82) * t245) * t230 + t324) * t229 + t49 * t366;
t331 = qJD(1) * t229 + qJD(5);
t338 = t230 * t361;
t249 = -t244 * t331 + t338;
t332 = qJD(5) * t229 + qJD(1);
t274 = t332 * t242;
t86 = t245 * t249 - t247 * t274;
t275 = t245 * t332;
t87 = t242 * t249 + t247 * t275;
t425 = t87 * rSges(6,1) + t86 * rSges(6,2);
t424 = rSges(3,1) * t246;
t423 = rSges(4,1) * t230;
t422 = rSges(3,2) * t243;
t421 = rSges(5,2) * t229;
t420 = rSges(3,3) * t247;
t234 = t244 * rSges(3,3);
t418 = -rSges(5,3) - qJ(4);
t400 = qJ(4) * t229;
t399 = qJ(4) * t230;
t388 = t230 * t244;
t386 = t241 * t247;
t385 = t242 * t247;
t384 = t244 * t242;
t383 = t244 * t245;
t382 = t245 * t247;
t174 = t229 * t382 - t384;
t175 = t229 * t385 + t383;
t109 = t175 * rSges(6,1) + t174 * rSges(6,2) + rSges(6,3) * t387;
t381 = pkin(7) * t387 + t109 + t236;
t176 = t229 * t383 + t385;
t177 = t229 * t384 - t382;
t321 = -t177 * rSges(6,1) - t176 * rSges(6,2);
t110 = rSges(6,3) * t388 - t321;
t380 = -pkin(4) * t247 + pkin(7) * t388 + t110;
t139 = t237 + t386 - t455;
t218 = t247 * t224;
t140 = -pkin(1) * t247 + t244 * t427 + t218;
t379 = t244 * t139 + t247 * t140;
t168 = pkin(3) * t387 + qJ(4) * t389;
t378 = -t140 - t168;
t192 = -t399 + t432;
t342 = t243 * t369;
t221 = pkin(2) * t342;
t377 = t192 * t369 + t221;
t360 = qJD(4) * t229;
t376 = qJ(4) * t338 + t247 * t360;
t345 = t229 * t369;
t368 = qJD(1) * t247;
t375 = rSges(4,2) * t345 + rSges(4,3) * t368;
t374 = t247 * t424 + t234;
t373 = t244 ^ 2 + t247 ^ 2;
t142 = Icges(4,3) * t244 + t247 * t292;
t372 = qJD(1) * t142;
t151 = Icges(5,1) * t244 - t247 * t295;
t371 = qJD(1) * t151;
t162 = Icges(3,3) * t244 + t247 * t293;
t370 = qJD(1) * t162;
t365 = qJD(2) * t243;
t364 = qJD(2) * t244;
t362 = qJD(2) * t246;
t358 = -pkin(3) - t436;
t357 = t247 * t422;
t354 = pkin(2) * t365;
t353 = pkin(2) * t362;
t103 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t387;
t105 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t387;
t107 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t387;
t286 = t105 * t245 + t107 * t242;
t38 = Icges(6,5) * t87 + Icges(6,6) * t86 - Icges(6,3) * t460;
t40 = Icges(6,4) * t87 + Icges(6,2) * t86 - Icges(6,6) * t460;
t42 = Icges(6,1) * t87 + Icges(6,4) * t86 - Icges(6,5) * t460;
t10 = (qJD(2) * t286 + t38) * t229 + (qJD(2) * t103 - t242 * t42 - t245 * t40 + (t105 * t242 - t107 * t245) * qJD(5)) * t230;
t17 = -t133 * t460 + t86 * t134 + t87 * t135 + t174 * t82 + t175 * t83 + t387 * t81;
t352 = t10 / 0.2e1 + t17 / 0.2e1;
t104 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t388;
t106 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t388;
t108 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t388;
t285 = t106 * t245 + t108 * t242;
t341 = t229 * t364;
t343 = t230 * t368;
t255 = -t341 + t343;
t273 = t331 * t247;
t84 = t245 * t273 + (t230 * t363 - t274) * t244;
t339 = t230 * t364;
t85 = t244 * t275 + (t273 + t339) * t242;
t37 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t255;
t39 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t255;
t41 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t255;
t11 = (qJD(2) * t285 + t37) * t229 + (qJD(2) * t104 - t242 * t41 - t245 * t39 + (t106 * t242 - t108 * t245) * qJD(5)) * t230;
t16 = t133 * t255 + t84 * t134 + t85 * t135 + t176 * t82 + t177 * t83 + t388 * t81;
t351 = t11 / 0.2e1 + t16 / 0.2e1;
t30 = t103 * t229 - t230 * t286;
t45 = t133 * t387 + t174 * t134 + t175 * t135;
t350 = -t30 / 0.2e1 - t45 / 0.2e1;
t31 = t104 * t229 - t230 * t285;
t46 = t133 * t388 + t134 * t176 + t135 * t177;
t349 = t31 / 0.2e1 + t46 / 0.2e1;
t231 = qJD(3) * t244;
t346 = qJD(3) * t247 + t241 * t369 + t244 * t354;
t348 = t244 * ((-t247 * t428 - t430) * qJD(1) - t346) + t247 * (-t247 * t354 + t231 + (t247 * t427 + t455) * qJD(1)) + t139 * t368;
t347 = t231 + t376;
t337 = -t367 / 0.2e1;
t336 = -t192 - t434;
t194 = rSges(4,1) * t229 + rSges(4,2) * t230;
t335 = -t194 - t434;
t320 = rSges(6,1) * t242 + rSges(6,2) * t245;
t136 = rSges(6,3) * t229 - t230 * t320;
t334 = pkin(7) * t229 + t136;
t333 = -t244 * t241 + t218;
t317 = t400 + t431;
t167 = t317 * t244;
t330 = t244 * t167 + t247 * t168 + t379;
t329 = rSges(5,1) * t368 + rSges(5,2) * t460 + rSges(5,3) * t338;
t206 = pkin(3) * t341;
t328 = t206 + t346;
t318 = rSges(5,3) * t230 + t421;
t327 = t318 + t336;
t326 = t85 * rSges(6,1) + t84 * rSges(6,2);
t325 = -t295 * qJD(2) / 0.2e1 + (t292 + t293) * t459;
t323 = -t422 + t424;
t322 = -rSges(4,2) * t229 + t423;
t319 = -rSges(5,2) * t230 + rSges(5,3) * t229;
t26 = t103 * t387 + t174 * t105 + t175 * t107;
t27 = t104 * t387 + t174 * t106 + t175 * t108;
t313 = t244 * t27 + t247 * t26;
t18 = t26 * t244 - t247 * t27;
t28 = t103 * t388 + t105 * t176 + t107 * t177;
t29 = t104 * t388 + t106 * t176 + t108 * t177;
t312 = t244 * t29 + t247 * t28;
t19 = t28 * t244 - t247 * t29;
t311 = t244 * t31 + t247 * t30;
t310 = t30 * t244 - t247 * t31;
t253 = t230 * t358 - t224 - t400;
t250 = t253 * t244;
t52 = (pkin(4) - t241) * t247 + t250 + t321;
t270 = t333 + t168;
t53 = t270 + t381;
t309 = t244 * t53 + t247 * t52;
t62 = -t110 * t229 + t136 * t388;
t63 = t229 * t109 - t136 * t387;
t308 = t244 * t63 + t247 * t62;
t306 = t229 * t418 - t224;
t251 = (rSges(5,2) - pkin(3)) * t230 + t306;
t73 = (rSges(5,1) - t241) * t247 + t251 * t244;
t157 = rSges(5,3) * t389 + t450;
t74 = t157 + t270;
t307 = t244 * t74 + t247 * t73;
t284 = t109 * t244 - t110 * t247;
t156 = rSges(4,1) * t387 + t451;
t154 = qJD(2) * t317 - qJD(4) * t230;
t272 = -t319 * qJD(2) - t154 - t353;
t271 = -pkin(1) - t323;
t269 = t167 * t368 + t244 * (pkin(3) * t343 + t244 * t360 - t206 + (t229 * t368 + t339) * qJ(4)) + t247 * (-pkin(3) * t460 - qJ(4) * t345 + t376) + t348;
t114 = t327 * t247;
t268 = -t224 - t322;
t267 = qJD(2) * t211;
t266 = qJD(2) * t194;
t265 = -t334 + t336;
t260 = qJD(2) * (Icges(5,4) * t229 + Icges(5,5) * t230);
t259 = qJD(2) * (-Icges(3,5) * t243 - Icges(3,6) * t246);
t258 = qJD(2) * (-Icges(4,5) * t229 - Icges(4,6) * t230);
t71 = t265 * t247;
t88 = (-rSges(6,1) * t245 + rSges(6,2) * t242) * t359 + (rSges(6,3) * t230 + t229 * t320) * qJD(2);
t252 = -t154 - t88 + (-pkin(7) * t230 - t433) * qJD(2);
t248 = rSges(3,2) * t342 + rSges(3,3) * t368 - t247 * t267;
t228 = pkin(4) * t368;
t200 = t323 * qJD(2);
t185 = t322 * qJD(2);
t172 = -t357 + t374;
t171 = t244 * t323 - t420;
t158 = -rSges(5,1) * t247 + t244 * t319;
t155 = -rSges(4,3) * t247 + t244 * t322;
t138 = t335 * t247;
t137 = t335 * t244;
t132 = t430 + (pkin(1) - t422) * t247 + t374;
t131 = t244 * t271 + t237 + t420;
t123 = t156 + t333;
t122 = (rSges(4,3) - t241) * t247 + t268 * t244;
t117 = t244 * t259 + t370;
t116 = -qJD(1) * t161 + t247 * t259;
t113 = t327 * t244;
t102 = qJD(1) * t447 + t247 * t260;
t101 = t244 * t260 + t371;
t92 = t244 * t258 + t372;
t91 = -qJD(1) * t141 + t247 * t258;
t90 = t211 * t364 + ((-rSges(3,3) - pkin(6)) * t244 + t271 * t247) * qJD(1);
t89 = (t237 + (-pkin(1) - t424) * t244) * qJD(1) + t248;
t78 = -t194 * t368 - t244 * t185 + (-t243 * t368 - t244 * t362) * pkin(2);
t77 = t194 * t369 + t221 + (-t185 - t353) * t247;
t70 = t265 * t244;
t69 = t244 * t162 - t276 * t247;
t68 = t244 * t161 - t454;
t67 = -t162 * t247 - t458;
t66 = -t161 * t247 - t244 * t277;
t65 = t194 * t364 + (t247 * t268 - t233) * qJD(1) + t346;
t64 = t231 + (-t224 - t423) * t369 + (-qJD(1) * t241 + qJD(2) * t335) * t247 + t375;
t61 = t244 * t278 + t247 * t447;
t60 = -t151 * t247 + t457;
t59 = -t244 * t447 + t453;
t58 = t244 * t151 + t279 * t247;
t57 = t244 * t142 - t280 * t247;
t56 = t244 * t141 - t452;
t55 = -t142 * t247 - t456;
t54 = -t141 * t247 - t244 * t281;
t51 = qJD(1) * t114 + t244 * t272;
t50 = t247 * t272 - t318 * t369 + t377;
t47 = t284 * t230;
t44 = -rSges(6,3) * t460 + t425;
t43 = rSges(6,3) * t255 + t326;
t36 = (-t360 + (t230 * t418 - t421) * qJD(2)) * t244 + (t247 * t251 - t235) * qJD(1) + t328;
t35 = (-t432 - t434) * t361 + (-t386 + (t306 - t431) * t244) * qJD(1) + t329 + t347;
t34 = t157 * t247 + t244 * t158 + t330;
t33 = qJD(1) * t71 + t244 * t252;
t32 = t247 * t252 + t334 * t369 + t377;
t25 = t244 * t380 + t247 * t381 + t330;
t24 = (-t360 + (t229 * t436 - t399) * qJD(2)) * t244 + (t247 * t253 - t236) * qJD(1) - t326 + t328;
t23 = t228 + (t229 * t358 - t434) * t361 + (t250 - t386) * qJD(1) + t347 + t425;
t22 = (-t136 * t364 - t43) * t229 + (-qJD(2) * t110 + t136 * t368 + t244 * t88) * t230;
t21 = (t136 * t361 + t44) * t229 + (qJD(2) * t109 + t136 * t369 - t247 * t88) * t230;
t15 = (qJD(1) * t158 + t329) * t247 + (t318 * t364 + (-t157 + t378 + t450) * qJD(1)) * t244 + t269;
t14 = t284 * t367 + (-t244 * t44 + t247 * t43 + (-t109 * t247 - t244 * t110) * qJD(1)) * t230;
t13 = t46 * t229 + t230 * t312;
t12 = t45 * t229 + t230 * t313;
t9 = -t104 * t460 + t86 * t106 + t87 * t108 + t174 * t39 + t175 * t41 + t37 * t387;
t8 = -t103 * t460 + t86 * t105 + t87 * t107 + t174 * t40 + t175 * t42 + t38 * t387;
t7 = t104 * t255 + t84 * t106 + t85 * t108 + t176 * t39 + t177 * t41 + t37 * t388;
t6 = t103 * t255 + t84 * t105 + t85 * t107 + t176 * t40 + t177 * t42 + t38 * t388;
t5 = (-pkin(7) * t340 + qJD(1) * t380 + t228 + t44) * t247 + (-pkin(7) * t341 + t43 + (t378 - t381 + t236) * qJD(1)) * t244 + t269;
t4 = qJD(1) * t313 + t8 * t244 - t247 * t9;
t3 = qJD(1) * t312 + t6 * t244 - t247 * t7;
t2 = (-qJD(2) * t313 + t17) * t229 + (-qJD(1) * t18 + qJD(2) * t45 + t244 * t9 + t247 * t8) * t230;
t1 = (-qJD(2) * t312 + t16) * t229 + (-qJD(1) * t19 + qJD(2) * t46 + t244 * t7 + t247 * t6) * t230;
t20 = [(t131 * t90 + t132 * t89) * t446 + (t122 * t65 + t123 * t64) * t445 + (t35 * t74 + t36 * t73) * t444 + (t23 * t53 + t24 * t52) * t443 + t324 - t245 * t135 * t359 + (-Icges(3,2) * t246 + t304 - t416) * t365 + (Icges(3,1) * t243 + t299 + t415) * t362 + (-t245 * t82 - t419) * t230 + (t302 + t290 + t468) * t367 + (t297 + t288 - t467) * t366; m(4) * (t122 * t77 + t123 * t78 + t137 * t64 + t138 * t65) + m(5) * (t113 * t35 + t114 * t36 + t50 * t73 + t51 * t74) + m(6) * (t23 * t70 + t24 * t71 + t32 * t52 + t33 * t53) + (m(3) * (-t131 * t200 - t211 * t90) + t325 * t247 + (t397 / 0.2e1 - t392 / 0.2e1) * qJD(2) - t351) * t247 + (m(3) * (-t132 * t200 - t211 * t89) + t325 * t244 + (-t395 / 0.2e1 + t390 / 0.2e1) * qJD(2) + t352) * t244 + ((-qJD(2) * t463 + t465 * t247) * t439 + (qJD(2) * t464 + t465 * t244) * t438) * t229 + ((qJD(2) * t462 + t466 * t247) * t439 + (-qJD(2) * t461 + t466 * t244) * t438) * t230 + ((t438 * t462 + t439 * t461) * t229 + (t463 * t438 + t439 * t464) * t230 + (-t132 * t435 + (t144 / 0.2e1 - t147 / 0.2e1) * t230 + (t146 / 0.2e1 - t149 / 0.2e1) * t229 - t350) * t247 + (t131 * t435 + (t449 / 0.2e1 + t143 / 0.2e1) * t230 + (t448 / 0.2e1 + t145 / 0.2e1) * t229 + t349) * t244) * qJD(1); t244 * t4 - t247 * t3 + (t25 * t5 + t32 * t71 + t33 * t70) * t443 + (t113 * t51 + t114 * t50 + t15 * t34) * t444 + (t138 * t77 + t137 * t78 + (t244 * t155 + t156 * t247 + t379) * ((qJD(1) * t155 - t247 * t266 + t375) * t247 + (-t244 * t266 + (-t140 - t156 + t451) * qJD(1)) * t244 + t348)) * t445 - t247 * ((t247 * t117 + (t67 + t454) * qJD(1)) * t247 + (t66 * qJD(1) + (-t164 * t362 - t166 * t365 + t370) * t244 + (-t116 + (t163 * t246 + t165 * t243) * qJD(2) - t276 * qJD(1)) * t247) * t244) + t244 * ((t244 * t91 + (t56 + t456) * qJD(1)) * t244 + (t57 * qJD(1) + (t143 * t366 + t145 * t367) * t247 + (-t92 + (-t144 * t230 - t146 * t229) * qJD(2) + (t142 - t281) * qJD(1)) * t244) * t247) + t244 * ((t244 * t102 + (t59 - t457) * qJD(1)) * t244 + (t58 * qJD(1) + (t366 * t449 + t367 * t448) * t247 + (-t101 + (t147 * t230 + t149 * t229) * qJD(2) + (t151 + t278) * qJD(1)) * t244) * t247) + t244 * ((t244 * t116 + (t68 + t458) * qJD(1)) * t244 + (t69 * qJD(1) + (t163 * t362 + t165 * t365) * t247 + (-t117 + (-t164 * t246 - t166 * t243) * qJD(2) + (t162 - t277) * qJD(1)) * t244) * t247) + ((t244 * t171 + t172 * t247) * ((qJD(1) * t171 + t248) * t247 + (-t244 * t267 + (-t172 - t357 + t234) * qJD(1)) * t244) + t373 * t211 * t200) * t446 - t247 * ((t247 * t92 + (t55 + t452) * qJD(1)) * t247 + (t54 * qJD(1) + (-t144 * t366 - t146 * t367 + t372) * t244 + (-t91 + (t143 * t230 + t145 * t229) * qJD(2) - t280 * qJD(1)) * t247) * t244) - t247 * ((t247 * t101 + (t60 - t453) * qJD(1)) * t247 + (t61 * qJD(1) + (t147 * t366 + t149 * t367 + t371) * t244 + (-t102 + (t229 * t448 + t230 * t449) * qJD(2) + t279 * qJD(1)) * t247) * t244) + (t19 + (-t54 - t61 - t66) * t247 + (t55 + t60 + t67) * t244) * t369 + (t18 + (-t56 - t59 - t68) * t247 + (t57 + t58 + t69) * t244) * t368; m(6) * (qJD(1) * t309 - t23 * t247 + t244 * t24) + m(4) * (t244 * t65 - t247 * t64 + (t122 * t247 + t123 * t244) * qJD(1)) + m(5) * (qJD(1) * t307 + t244 * t36 - t247 * t35); m(6) * (t244 * t32 - t247 * t33 + (t244 * t70 + t247 * t71) * qJD(1)) + m(5) * (t244 * t50 - t247 * t51 + (t113 * t244 + t114 * t247) * qJD(1)) + m(4) * (t244 * t77 - t247 * t78 + (t137 * t244 + t138 * t247) * qJD(1)); 0; 0.2e1 * (t307 * t442 + t309 * t441) * t366 + 0.2e1 * ((t23 * t244 + t24 * t247 + t368 * t53 - t369 * t52) * t441 + (t244 * t35 + t247 * t36 + t368 * t74 - t369 * t73) * t442) * t229; 0.2e1 * ((t361 * t71 + t364 * t70 - t5) * t441 + (t113 * t364 + t114 * t361 - t15) * t442) * t230 + 0.2e1 * ((qJD(2) * t25 + t244 * t33 + t247 * t32 + t368 * t70 - t369 * t71) * t441 + (qJD(2) * t34 + t113 * t368 - t114 * t369 + t244 * t51 + t247 * t50) * t442) * t229; 0; 0.4e1 * (t442 + t441) * (-0.1e1 + t373) * t229 * t366; m(6) * (t21 * t53 + t22 * t52 + t23 * t63 + t24 * t62) + (-t244 * t349 + t247 * t350) * t367 + (t352 * t247 + t351 * t244 + (t244 * t350 + t247 * t349) * qJD(1)) * t230 + t426; m(6) * (t14 * t25 + t21 * t70 + t22 * t71 + t32 * t62 + t33 * t63 - t47 * t5) + (t12 * t429 - t1 / 0.2e1 + t18 * t337 + (qJD(1) * t30 - t11) * t440) * t247 + (t2 / 0.2e1 + t19 * t337 + t13 * t429 + (qJD(1) * t31 + t10) * t440) * t244 + (t3 * t439 + t4 * t437 + t310 * t459 + (t19 * t437 - t244 * t18 / 0.2e1) * qJD(1)) * t230; m(6) * (qJD(1) * t308 - t21 * t247 + t22 * t244); m(6) * ((qJD(2) * t308 - t14) * t230 + (-qJD(2) * t47 + t21 * t244 + t22 * t247 + (-t244 * t62 + t247 * t63) * qJD(1)) * t229); (-t14 * t47 + t21 * t63 + t22 * t62) * t443 + ((-t247 * t12 - t244 * t13 - t229 * t311) * qJD(2) + t426) * t229 + (t247 * t2 + t244 * t1 + t229 * (t10 * t247 + t11 * t244) + (t49 * t229 + t230 * t311) * qJD(2) + (-t244 * t12 + t247 * t13 - t229 * t310) * qJD(1)) * t230;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t20(1), t20(2), t20(4), t20(7), t20(11); t20(2), t20(3), t20(5), t20(8), t20(12); t20(4), t20(5), t20(6), t20(9), t20(13); t20(7), t20(8), t20(9), t20(10), t20(14); t20(11), t20(12), t20(13), t20(14), t20(15);];
Mq = res;
