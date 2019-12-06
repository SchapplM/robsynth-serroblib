% Calculate time derivative of joint inertia matrix for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:05
% EndTime: 2019-12-05 18:27:36
% DurationCPUTime: 11.26s
% Computational Cost: add. (14448->620), mult. (12086->853), div. (0->0), fcn. (9434->10), ass. (0->356)
t257 = sin(qJ(1));
t259 = cos(qJ(1));
t252 = qJ(2) + pkin(9);
t240 = qJ(4) + t252;
t231 = qJ(5) + t240;
t225 = sin(t231);
t431 = rSges(6,2) * t225;
t226 = cos(t231);
t436 = rSges(6,1) * t226;
t328 = -t431 + t436;
t142 = t257 * rSges(6,3) + t259 * t328;
t245 = t257 * rSges(5,3);
t229 = sin(t240);
t230 = cos(t240);
t437 = rSges(5,1) * t230;
t329 = -rSges(5,2) * t229 + t437;
t154 = t259 * t329 + t245;
t430 = rSges(6,2) * t226;
t186 = rSges(6,1) * t225 + t430;
t251 = qJD(2) + qJD(4);
t239 = qJD(5) + t251;
t284 = t186 * t239;
t432 = rSges(5,2) * t230;
t195 = rSges(5,1) * t229 + t432;
t283 = t195 * t251;
t256 = sin(qJ(2));
t258 = cos(qJ(2));
t425 = Icges(3,4) * t258;
t318 = -Icges(3,2) * t256 + t425;
t178 = Icges(3,6) * t257 + t259 * t318;
t426 = Icges(3,4) * t256;
t324 = Icges(3,1) * t258 - t426;
t180 = Icges(3,5) * t257 + t259 * t324;
t296 = t178 * t256 - t180 * t258;
t467 = t257 * t296;
t237 = sin(t252);
t238 = cos(t252);
t423 = Icges(4,4) * t238;
t316 = -Icges(4,2) * t237 + t423;
t163 = Icges(4,6) * t257 + t259 * t316;
t424 = Icges(4,4) * t237;
t322 = Icges(4,1) * t238 - t424;
t165 = Icges(4,5) * t257 + t259 * t322;
t298 = t163 * t237 - t165 * t238;
t466 = t257 * t298;
t421 = Icges(5,4) * t230;
t314 = -Icges(5,2) * t229 + t421;
t150 = Icges(5,6) * t257 + t259 * t314;
t422 = Icges(5,4) * t229;
t320 = Icges(5,1) * t230 - t422;
t152 = Icges(5,5) * t257 + t259 * t320;
t300 = t150 * t229 - t152 * t230;
t465 = t257 * t300;
t419 = Icges(6,4) * t226;
t313 = -Icges(6,2) * t225 + t419;
t132 = Icges(6,6) * t257 + t259 * t313;
t420 = Icges(6,4) * t225;
t319 = Icges(6,1) * t226 - t420;
t134 = Icges(6,5) * t257 + t259 * t319;
t302 = t132 * t225 - t134 * t226;
t464 = t257 * t302;
t248 = t258 * pkin(2);
t232 = t248 + pkin(1);
t441 = pkin(1) - t232;
t463 = t257 * t441;
t177 = -Icges(3,6) * t259 + t257 * t318;
t179 = -Icges(3,5) * t259 + t257 * t324;
t297 = t177 * t256 - t179 * t258;
t462 = t259 * t297;
t162 = -Icges(4,6) * t259 + t257 * t316;
t164 = -Icges(4,5) * t259 + t257 * t322;
t299 = t162 * t237 - t164 * t238;
t461 = t259 * t299;
t149 = -Icges(5,6) * t259 + t257 * t314;
t151 = -Icges(5,5) * t259 + t257 * t320;
t301 = t149 * t229 - t151 * t230;
t460 = t259 * t301;
t131 = -Icges(6,6) * t259 + t257 * t313;
t133 = -Icges(6,5) * t259 + t257 * t319;
t303 = t131 * t225 - t133 * t226;
t459 = t259 * t303;
t246 = t257 * rSges(4,3);
t434 = rSges(4,2) * t237;
t458 = -t259 * t434 + t246;
t444 = pkin(2) * t256;
t334 = -pkin(3) * t237 - t444;
t291 = -t195 + t334;
t119 = t291 * t257;
t120 = t291 * t259;
t457 = t119 * t257 + t120 * t259;
t309 = Icges(6,5) * t226 - Icges(6,6) * t225;
t129 = -Icges(6,3) * t259 + t257 * t309;
t456 = qJD(1) * t129;
t310 = Icges(5,5) * t230 - Icges(5,6) * t229;
t147 = -Icges(5,3) * t259 + t257 * t310;
t455 = qJD(1) * t147;
t311 = Icges(4,5) * t238 - Icges(4,6) * t237;
t160 = -Icges(4,3) * t259 + t257 * t311;
t312 = Icges(3,5) * t258 - Icges(3,6) * t256;
t175 = -Icges(3,3) * t259 + t257 * t312;
t182 = Icges(6,2) * t226 + t420;
t183 = Icges(6,1) * t225 + t419;
t295 = t182 * t225 - t183 * t226;
t454 = qJD(1) * t295 + t309 * t239;
t193 = Icges(5,2) * t230 + t422;
t194 = Icges(5,1) * t229 + t421;
t294 = t193 * t229 - t194 * t230;
t453 = qJD(1) * t294 + t310 * t251;
t255 = -qJ(3) - pkin(6);
t250 = -pkin(7) + t255;
t373 = t250 - t255;
t228 = pkin(3) * t238;
t206 = t228 + t232;
t380 = t206 - t232;
t117 = t257 * t380 + t259 * t373;
t243 = -pkin(8) + t250;
t374 = t243 - t250;
t174 = pkin(4) * t230 + t206;
t382 = t174 - t206;
t90 = t257 * t382 + t259 * t374;
t452 = 2 * m(3);
t451 = 2 * m(4);
t450 = 2 * m(5);
t449 = 2 * m(6);
t253 = t257 ^ 2;
t254 = t259 ^ 2;
t448 = t257 / 0.2e1;
t447 = -t259 / 0.2e1;
t217 = rSges(3,1) * t256 + rSges(3,2) * t258;
t446 = m(3) * t217;
t445 = m(5) * t195;
t443 = t257 * pkin(6);
t249 = t259 * pkin(6);
t29 = -t129 * t259 - t257 * t303;
t130 = Icges(6,3) * t257 + t259 * t309;
t30 = -t130 * t259 - t464;
t396 = t182 * t239;
t70 = qJD(1) * t132 - t257 * t396;
t343 = t133 * t239 + t70;
t395 = t183 * t239;
t72 = qJD(1) * t134 - t257 * t395;
t345 = t131 * t239 - t72;
t371 = qJD(1) * t130;
t391 = t226 * t239;
t392 = t225 * t239;
t181 = Icges(6,5) * t225 + Icges(6,6) * t226;
t278 = t239 * t181;
t67 = -t259 * t278 - t456;
t68 = -t257 * t278 + t371;
t69 = -qJD(1) * t131 - t259 * t396;
t71 = -qJD(1) * t133 - t259 * t395;
t2 = (t259 * t68 + (t30 + t459) * qJD(1)) * t259 + (t29 * qJD(1) + (-t132 * t391 - t134 * t392 - t225 * t69 + t226 * t71 + t371) * t257 + (-t67 + t345 * t226 + t343 * t225 + (-t129 - t302) * qJD(1)) * t259) * t257;
t442 = t259 * t2;
t440 = -pkin(6) - t255;
t439 = rSges(3,1) * t258;
t438 = rSges(4,1) * t238;
t435 = rSges(3,2) * t256;
t429 = rSges(3,3) * t259;
t247 = t257 * rSges(3,3);
t428 = rSges(6,3) - t243;
t173 = t259 * t174;
t196 = t259 * t206;
t91 = -t257 * t374 + t173 - t196;
t427 = -t142 - t91;
t404 = t162 * t238;
t403 = t163 * t238;
t402 = t164 * t237;
t401 = t165 * t237;
t400 = t177 * t258;
t399 = t178 * t258;
t398 = t179 * t256;
t397 = t180 * t256;
t394 = t193 * t251;
t393 = t194 * t251;
t390 = t229 * t251;
t389 = t230 * t251;
t388 = t239 * t259;
t387 = t251 * t259;
t386 = t257 * t243;
t219 = t259 * t232;
t118 = -t257 * t373 + t196 - t219;
t159 = -pkin(1) * t259 + t257 * t440 + t219;
t385 = -t118 - t159;
t141 = -rSges(6,3) * t259 + t257 * t328;
t73 = t257 * t141 + t259 * t142;
t153 = -rSges(5,3) * t259 + t257 * t329;
t78 = t257 * t153 + t259 * t154;
t158 = t255 * t259 + t249 - t463;
t384 = t257 * t158 + t259 * t159;
t367 = qJD(1) * t257;
t352 = t229 * t367;
t383 = pkin(4) * t352 + t186 * t367;
t366 = qJD(1) * t259;
t381 = rSges(6,3) * t366 + t367 * t431;
t379 = rSges(5,2) * t352 + rSges(5,3) * t366;
t351 = t237 * t367;
t378 = rSges(4,2) * t351 + rSges(4,3) * t366;
t350 = t256 * t367;
t223 = pkin(2) * t350;
t377 = pkin(3) * t351 + t223;
t363 = qJD(2) * t256;
t355 = pkin(2) * t363;
t376 = t255 * t367 + t257 * t355;
t375 = t259 * t439 + t247;
t372 = t253 + t254;
t148 = Icges(5,3) * t257 + t259 * t310;
t370 = qJD(1) * t148;
t161 = Icges(4,3) * t257 + t259 * t311;
t369 = qJD(1) * t161;
t176 = Icges(3,3) * t257 + t259 * t312;
t368 = qJD(1) * t176;
t365 = qJD(2) * t237;
t364 = qJD(2) * t238;
t362 = qJD(2) * t257;
t361 = qJD(2) * t258;
t31 = t257 * t129 - t459;
t32 = t257 * t130 - t259 * t302;
t342 = t134 * t239 + t69;
t344 = -t132 * t239 + t71;
t360 = t257 * ((t257 * t67 + (t31 + t464) * qJD(1)) * t257 + (t32 * qJD(1) + (t131 * t391 + t133 * t392 + t225 * t70 - t226 * t72 - t456) * t259 + (-t68 + t344 * t226 - t342 * t225 + (t130 - t303) * qJD(1)) * t257) * t259) + (t30 * t257 - t259 * t29) * t367 + (t32 * t257 - t259 * t31) * t366;
t359 = t141 * t366 + t257 * (qJD(1) * t142 - t257 * t284) + t259 * (-t388 * t430 + (-t225 * t388 - t226 * t367) * rSges(6,1) + t381);
t358 = t153 * t366 + t257 * (qJD(1) * t154 - t257 * t283) + t259 * (-t387 * t432 + (-t229 * t387 - t230 * t367) * rSges(5,1) + t379);
t357 = t259 * t435;
t241 = qJD(3) * t257;
t337 = t259 * t355;
t242 = qJD(3) * t259;
t353 = -t242 - t376;
t354 = t257 * ((-t259 * t441 - t443) * qJD(1) + t353) + t259 * (-t337 + t241 + (t259 * t440 + t463) * qJD(1)) + t158 * t366;
t349 = t367 / 0.2e1;
t348 = t366 / 0.2e1;
t201 = rSges(4,1) * t237 + rSges(4,2) * t238;
t347 = -t201 - t444;
t346 = -pkin(4) * t229 - t186;
t84 = qJD(1) * t152 - t257 * t393;
t341 = t149 * t251 - t84;
t83 = -qJD(1) * t151 - t259 * t393;
t340 = -t150 * t251 + t83;
t82 = qJD(1) * t150 - t257 * t394;
t339 = t151 * t251 + t82;
t81 = -qJD(1) * t149 - t259 * t394;
t338 = t152 * t251 + t81;
t24 = t257 * t90 + t259 * t91 + t73;
t336 = t257 * t117 + t259 * t118 + t384;
t139 = t328 * t239;
t335 = -pkin(4) * t389 - t139;
t333 = t360 - t442;
t37 = -t147 * t259 - t257 * t301;
t38 = -t148 * t259 - t465;
t39 = t257 * t147 - t460;
t40 = t257 * t148 - t259 * t300;
t192 = Icges(5,5) * t229 + Icges(5,6) * t230;
t275 = t251 * t192;
t79 = -t259 * t275 - t455;
t80 = -t257 * t275 + t370;
t332 = (t38 * t257 - t259 * t37) * t367 + (t40 * t257 - t259 * t39) * t366 + t257 * ((t257 * t79 + (t39 + t465) * qJD(1)) * t257 + (t40 * qJD(1) + (t149 * t389 + t151 * t390 + t229 * t82 - t230 * t84 - t455) * t259 + (-t80 + t340 * t230 - t338 * t229 + (t148 - t301) * qJD(1)) * t257) * t259) + t360;
t331 = -t435 + t439;
t330 = -t434 + t438;
t287 = -t174 - t328;
t76 = t257 * t287 + t259 * t428;
t77 = t142 + t173 - t386;
t325 = t257 * t77 + t259 * t76;
t323 = Icges(3,1) * t256 + t425;
t321 = Icges(4,1) * t237 + t423;
t317 = Icges(3,2) * t258 + t426;
t315 = Icges(4,2) * t238 + t424;
t266 = t334 + t346;
t100 = t266 * t257;
t101 = t266 * t259;
t308 = t100 * t257 + t101 * t259;
t288 = -t206 - t329;
t102 = (rSges(5,3) - t250) * t259 + t288 * t257;
t103 = -t257 * t250 + t154 + t196;
t307 = t102 * t259 + t103 * t257;
t121 = t346 * t257;
t122 = t346 * t259;
t304 = t121 * t257 + t122 * t259;
t200 = t334 * qJD(2);
t155 = -pkin(4) * t390 + t200;
t140 = t259 * t155;
t187 = t259 * t200;
t221 = t250 * t367;
t293 = t257 * (t221 + (t155 - t200) * t257 + (t259 * t382 - t386) * qJD(1)) + t259 * (-qJD(1) * t90 + t140 - t187) + t90 * t366 + t359;
t170 = t259 * t438 + t458;
t292 = -pkin(1) - t331;
t290 = t117 * t366 + t257 * (t257 * t200 + t366 * t380 - t221 + t376) + t259 * (-qJD(1) * t117 + t187 + t337) + t354;
t289 = -t232 - t330;
t285 = (-t248 - t228) * qJD(2);
t282 = qJD(2) * t217;
t281 = qJD(2) * t201;
t136 = t313 * t239;
t137 = t319 * t239;
t261 = qJD(1) * t181 + (t137 - t396) * t226 + (-t136 - t395) * t225;
t274 = (t225 * t344 + t226 * t342 + t257 * t454 + t261 * t259) * t448 + (-t225 * t345 + t226 * t343 + t261 * t257 - t259 * t454) * t447 + (t131 * t226 + t133 * t225 - t181 * t259 - t257 * t295) * t349 + (t132 * t226 + t134 * t225 + t257 * t181 - t259 * t295) * t348;
t273 = qJD(2) * t323;
t272 = qJD(2) * t321;
t271 = qJD(2) * t317;
t270 = qJD(2) * t315;
t269 = qJD(2) * (-Icges(3,5) * t256 - Icges(3,6) * t258);
t268 = qJD(2) * (-Icges(4,5) * t237 - Icges(4,6) * t238);
t172 = t329 * t251;
t267 = -t172 + t285;
t4 = (t259 * t80 + (t38 + t460) * qJD(1)) * t259 + (t37 * qJD(1) + (-t150 * t389 - t152 * t390 - t229 * t81 + t230 * t83 + t370) * t257 + (-t79 + t341 * t230 + t339 * t229 + (-t147 - t300) * qJD(1)) * t259) * t257;
t265 = (-t4 - t2) * t259 + t332;
t264 = rSges(3,2) * t350 + rSges(3,3) * t366 - t259 * t282;
t263 = t285 + t335;
t167 = t314 * t251;
t168 = t320 * t251;
t260 = qJD(1) * t192 + (t168 - t394) * t230 + (-t167 - t393) * t229;
t262 = t274 + (t229 * t340 + t230 * t338 + t257 * t453 + t260 * t259) * t448 + (-t229 * t341 + t230 * t339 + t260 * t257 - t259 * t453) * t447 + (t149 * t230 + t151 * t229 - t192 * t259 - t257 * t294) * t349 + (t150 * t230 + t152 * t229 + t257 * t192 - t259 * t294) * t348;
t207 = t331 * qJD(2);
t191 = t330 * qJD(2);
t185 = -t357 + t375;
t184 = t257 * t331 - t429;
t169 = -rSges(4,3) * t259 + t257 * t330;
t157 = t347 * t259;
t156 = t347 * t257;
t146 = t443 + (pkin(1) - t435) * t259 + t375;
t145 = t257 * t292 + t249 + t429;
t116 = -t257 * t255 + t170 + t219;
t115 = (rSges(4,3) - t255) * t259 + t289 * t257;
t110 = t257 * t269 + t368;
t109 = -qJD(1) * t175 + t259 * t269;
t95 = t257 * t268 + t369;
t94 = -qJD(1) * t160 + t259 * t268;
t93 = t217 * t362 + ((-rSges(3,3) - pkin(6)) * t257 + t292 * t259) * qJD(1);
t92 = (t249 + (-pkin(1) - t439) * t257) * qJD(1) + t264;
t89 = -t201 * t366 - t257 * t191 + (-t256 * t366 - t257 * t361) * pkin(2);
t88 = t201 * t367 + t223 + (-pkin(2) * t361 - t191) * t259;
t60 = t257 * t176 - t296 * t259;
t59 = t257 * t175 - t462;
t58 = -t176 * t259 - t467;
t57 = -t175 * t259 - t257 * t297;
t54 = t201 * t362 + (t259 * t289 - t246) * qJD(1) - t353;
t53 = t241 + (-t232 - t438) * t367 + (-qJD(1) * t255 + qJD(2) * t347) * t259 + t378;
t50 = -t186 * t366 - t257 * t139 + (-t229 * t366 - t257 * t389) * pkin(4);
t49 = t259 * t335 + t383;
t48 = qJD(1) * t120 + t257 * t267;
t47 = t195 * t367 + t259 * t267 + t377;
t44 = t257 * t161 - t298 * t259;
t43 = t257 * t160 - t461;
t42 = -t161 * t259 - t466;
t41 = -t160 * t259 - t257 * t299;
t36 = t221 + t242 + (-t200 + t283) * t257 + (t259 * t288 - t245) * qJD(1);
t35 = t187 + t241 - t259 * t283 + (-t250 * t259 + (-t206 - t437) * t257) * qJD(1) + t379;
t28 = qJD(1) * t101 + t257 * t263;
t27 = t259 * t263 + t377 + t383;
t26 = t242 + (-t155 + t284) * t257 + (-t257 * t428 + t259 * t287) * qJD(1);
t25 = t140 + t241 - t259 * t284 + (-t243 * t259 + (-t174 - t436) * t257) * qJD(1) + t381;
t23 = -t154 * t367 + t358;
t22 = t336 + t78;
t21 = -t142 * t367 + t359;
t8 = t336 + t24;
t7 = t367 * t427 + t293;
t6 = (-t154 + t385) * t367 + t290 + t358;
t5 = (t385 + t427) * t367 + t290 + t293;
t1 = [(t25 * t77 + t26 * t76) * t449 + t183 * t391 + t225 * t137 - t182 * t392 + t226 * t136 + (t102 * t36 + t103 * t35) * t450 + t194 * t389 + t229 * t168 - t193 * t390 + t230 * t167 + (t115 * t54 + t116 * t53) * t451 + (t145 * t93 + t146 * t92) * t452 + (t322 - t315) * t365 + (t321 + t316) * t364 + (t324 - t317) * t363 + (t323 + t318) * t361; m(3) * ((-t257 * t92 - t259 * t93) * t217 + (-t145 * t259 - t146 * t257) * t207) + ((t399 / 0.2e1 + t397 / 0.2e1 - t146 * t446 + t403 / 0.2e1 + t401 / 0.2e1) * t259 + (t145 * t446 + t404 / 0.2e1 + t402 / 0.2e1 + t400 / 0.2e1 + t398 / 0.2e1) * t257) * qJD(1) + t262 + m(6) * (t100 * t25 + t101 * t26 + t27 * t76 + t28 * t77) + m(5) * (t102 * t47 + t103 * t48 + t119 * t35 + t120 * t36) + m(4) * (t115 * t88 + t116 * t89 + t156 * t53 + t157 * t54) + ((-qJD(1) * t177 - t259 * t271) * t258 + (-qJD(1) * t179 - t259 * t273) * t256 + t237 * (-qJD(1) * t164 - t259 * t272) + t238 * (-qJD(1) * t162 - t259 * t270) + (-t296 - t298) * qJD(2)) * t448 + ((qJD(1) * t178 - t257 * t271) * t258 + (qJD(1) * t180 - t257 * t273) * t256 + t237 * (qJD(1) * t165 - t257 * t272) + t238 * (qJD(1) * t163 - t257 * t270) + (-t297 - t299) * qJD(2)) * t447 + (t312 + t311) * qJD(2) * (t254 / 0.2e1 + t253 / 0.2e1); t332 + (t100 * t28 + t101 * t27 + t5 * t8) * t449 + (t119 * t48 + t120 * t47 + t22 * t6) * t450 + t257 * ((t257 * t94 + (t43 + t466) * qJD(1)) * t257 + (t44 * qJD(1) + (t162 * t364 + t164 * t365) * t259 + (-t95 + (-t401 - t403) * qJD(2) + (t161 - t299) * qJD(1)) * t257) * t259) - t259 * ((t259 * t95 + (t42 + t461) * qJD(1)) * t259 + (t41 * qJD(1) + (-t163 * t364 - t165 * t365 + t369) * t257 + (-t94 + (t402 + t404) * qJD(2) - t298 * qJD(1)) * t259) * t257) + t257 * ((t257 * t109 + (t59 + t467) * qJD(1)) * t257 + (t60 * qJD(1) + (t177 * t361 + t179 * t363) * t259 + (-t110 + (-t397 - t399) * qJD(2) + (t176 - t297) * qJD(1)) * t257) * t259) - t259 * ((t259 * t110 + (t58 + t462) * qJD(1)) * t259 + (t57 * qJD(1) + (-t178 * t361 - t180 * t363 + t368) * t257 + (-t109 + (t398 + t400) * qJD(2) - t296 * qJD(1)) * t259) * t257) + (t157 * t88 + t156 * t89 + (t257 * t169 + t170 * t259 + t384) * ((qJD(1) * t169 - t259 * t281 + t378) * t259 + (-t257 * t281 + (-t159 - t170 + t458) * qJD(1)) * t257 + t354)) * t451 + ((t257 * t184 + t185 * t259) * ((qJD(1) * t184 + t264) * t259 + (-t257 * t282 + (-t185 - t357 + t247) * qJD(1)) * t257) + t372 * t217 * t207) * t452 - t442 - t259 * t4 + ((-t41 - t57) * t259 + (t42 + t58) * t257) * t367 + ((-t43 - t59) * t259 + (t44 + t60) * t257) * t366; m(6) * (qJD(1) * t325 - t25 * t259 + t257 * t26) + m(5) * (qJD(1) * t307 + t257 * t36 - t259 * t35) + m(4) * (t257 * t54 - t259 * t53 + (t115 * t259 + t116 * t257) * qJD(1)); m(6) * (qJD(1) * t308 + t257 * t27 - t259 * t28) + m(5) * (qJD(1) * t457 + t257 * t47 - t259 * t48) + m(4) * (t257 * t88 - t259 * t89 + (t156 * t257 + t157 * t259) * qJD(1)); 0; m(6) * (t121 * t25 + t122 * t26 + t49 * t76 + t50 * t77) + (-t257 * t35 - t259 * t36 + (t102 * t257 - t103 * t259) * qJD(1)) * t445 - m(5) * t307 * t172 + t262; m(6) * (t100 * t50 + t101 * t49 + t121 * t28 + t122 * t27 + t24 * t5 + t7 * t8) + m(5) * (-t172 * t457 + t23 * t22 + t78 * t6) + (-t257 * t48 - t259 * t47 + (-t119 * t259 + t120 * t257) * qJD(1)) * t445 + t265; m(6) * (qJD(1) * t304 + t49 * t257 - t259 * t50); (t172 * t195 * t372 + t23 * t78) * t450 + (t121 * t50 + t122 * t49 + t24 * t7) * t449 + t265; m(6) * (-t325 * t139 + (-t25 * t257 - t259 * t26 + (t257 * t76 - t259 * t77) * qJD(1)) * t186) + t274; m(6) * (t21 * t8 + t73 * t5 - t308 * t139 + (-t257 * t28 - t259 * t27 + (-t100 * t259 + t101 * t257) * qJD(1)) * t186) + t333; 0; m(6) * (t21 * t24 + t73 * t7 - t304 * t139 + (-t257 * t50 - t259 * t49 + (-t121 * t259 + t122 * t257) * qJD(1)) * t186) + t333; (t139 * t186 * t372 + t21 * t73) * t449 + t333;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
