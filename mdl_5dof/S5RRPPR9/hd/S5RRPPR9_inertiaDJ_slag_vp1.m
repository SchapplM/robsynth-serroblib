% Calculate time derivative of joint inertia matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:42
% EndTime: 2019-12-31 19:41:05
% DurationCPUTime: 15.05s
% Computational Cost: add. (6398->699), mult. (17240->986), div. (0->0), fcn. (15675->6), ass. (0->339)
t246 = cos(qJ(2));
t243 = sin(qJ(2));
t408 = Icges(5,4) * t243;
t411 = Icges(3,4) * t243;
t451 = t408 + t411 + (Icges(5,1) + Icges(3,2)) * t246;
t402 = Icges(4,5) * t246;
t407 = Icges(5,4) * t246;
t410 = Icges(3,4) * t246;
t450 = -t402 + t407 + t410 + (Icges(3,1) + Icges(4,1) + Icges(5,2)) * t243;
t449 = t451 * qJD(2);
t244 = sin(qJ(1));
t447 = -t244 / 0.2e1;
t446 = t244 / 0.2e1;
t247 = cos(qJ(1));
t445 = -t247 / 0.2e1;
t429 = t247 / 0.2e1;
t444 = -qJD(1) / 0.2e1;
t425 = qJD(1) / 0.2e1;
t363 = qJD(2) * t247;
t345 = t243 * t363;
t368 = qJD(1) * t246;
t443 = t244 * t368 + t345;
t242 = sin(qJ(5));
t245 = cos(qJ(5));
t405 = Icges(6,4) * t245;
t298 = -Icges(6,2) * t242 + t405;
t142 = Icges(6,6) * t243 - t246 * t298;
t406 = Icges(6,4) * t242;
t304 = Icges(6,1) * t245 - t406;
t149 = Icges(6,5) * t243 - t246 * t304;
t442 = -t142 * t242 + t149 * t245;
t441 = qJD(2) / 0.2e1;
t303 = -Icges(3,2) * t243 + t410;
t148 = Icges(3,6) * t244 + t247 * t303;
t310 = Icges(3,1) * t246 - t411;
t155 = Icges(3,5) * t244 + t247 * t310;
t281 = t148 * t243 - t155 * t246;
t265 = t281 * t244;
t147 = -Icges(3,6) * t247 + t244 * t303;
t154 = -Icges(3,5) * t247 + t244 * t310;
t282 = t147 * t243 - t154 * t246;
t266 = t282 * t247;
t299 = -Icges(5,2) * t246 + t408;
t144 = -Icges(5,6) * t244 + t247 * t299;
t305 = Icges(5,1) * t243 - t407;
t151 = -Icges(5,5) * t244 + t247 * t305;
t283 = t144 * t246 - t151 * t243;
t267 = t283 * t244;
t143 = Icges(5,6) * t247 + t244 * t299;
t150 = Icges(5,5) * t247 + t244 * t305;
t284 = t143 * t246 - t150 * t243;
t268 = t284 * t247;
t296 = Icges(4,3) * t243 + t402;
t139 = Icges(4,6) * t244 + t247 * t296;
t403 = Icges(4,5) * t243;
t308 = Icges(4,1) * t246 + t403;
t153 = Icges(4,4) * t244 + t247 * t308;
t285 = t139 * t243 + t153 * t246;
t269 = t285 * t244;
t138 = -Icges(4,6) * t247 + t244 * t296;
t152 = -Icges(4,4) * t247 + t244 * t308;
t286 = t138 * t243 + t152 * t246;
t270 = t286 * t247;
t391 = t243 * t247;
t440 = -rSges(3,2) * t391 + rSges(3,3) * t244;
t294 = Icges(5,5) * t243 - Icges(5,6) * t246;
t136 = Icges(5,3) * t247 + t244 * t294;
t297 = Icges(3,5) * t246 - Icges(3,6) * t243;
t140 = -Icges(3,3) * t247 + t244 * t297;
t301 = Icges(4,4) * t246 + Icges(4,6) * t243;
t145 = -Icges(4,2) * t247 + t244 * t301;
t338 = qJD(1) * t243 + qJD(5);
t343 = t246 * t363;
t439 = t244 * t338 - t343;
t438 = 2 * m(3);
t437 = 2 * m(4);
t436 = 2 * m(5);
t435 = 2 * m(6);
t240 = t244 ^ 2;
t241 = t247 ^ 2;
t434 = m(4) / 0.2e1;
t433 = m(5) / 0.2e1;
t432 = m(6) / 0.2e1;
t431 = -pkin(2) - pkin(3);
t430 = t243 / 0.2e1;
t428 = -rSges(4,1) - pkin(2);
t427 = rSges(6,3) + pkin(7);
t207 = rSges(3,1) * t243 + rSges(3,2) * t246;
t426 = m(3) * t207;
t237 = t244 * pkin(6);
t424 = -pkin(4) - qJ(3);
t293 = Icges(6,5) * t245 - Icges(6,6) * t242;
t360 = qJD(5) * t246;
t100 = (Icges(6,5) * t242 + Icges(6,6) * t245) * t360 + (Icges(6,3) * t246 + t243 * t293) * qJD(2);
t107 = (Icges(6,2) * t245 + t406) * t360 + (Icges(6,6) * t246 + t243 * t298) * qJD(2);
t114 = (Icges(6,1) * t242 + t405) * t360 + (Icges(6,5) * t246 + t243 * t304) * qJD(2);
t135 = Icges(6,3) * t243 - t246 * t293;
t364 = qJD(2) * t246;
t366 = qJD(2) * t243;
t249 = t243 * t100 + t135 * t364 + t442 * t366 + (-t114 * t246 + t142 * t360) * t245 + (t107 * t246 + t149 * t360) * t242;
t63 = t135 * t243 - t246 * t442;
t423 = t243 * t249 + t364 * t63;
t339 = -qJD(5) * t243 - qJD(1);
t279 = t339 * t247;
t83 = t242 * t439 + t245 * t279;
t84 = t242 * t279 - t245 * t439;
t422 = rSges(6,1) * t84 + rSges(6,2) * t83;
t421 = rSges(4,1) * t243;
t420 = rSges(4,2) * t247;
t419 = rSges(5,2) * t243;
t418 = rSges(3,3) * t247;
t417 = rSges(5,3) * t247;
t236 = t244 * rSges(4,2);
t416 = -rSges(5,1) - qJ(3);
t415 = -rSges(4,3) - qJ(3);
t414 = -rSges(5,3) - qJ(4);
t331 = pkin(4) * t243 + pkin(7) * t246;
t386 = t245 * t247;
t389 = t244 * t242;
t172 = -t243 * t389 + t386;
t388 = t244 * t245;
t173 = t242 * t247 + t243 * t388;
t324 = -rSges(6,1) * t173 - rSges(6,2) * t172;
t387 = t244 * t246;
t96 = rSges(6,3) * t387 - t324;
t413 = t244 * t331 + t96;
t385 = t246 * t247;
t174 = -t242 * t391 - t388;
t175 = t243 * t386 - t389;
t97 = rSges(6,1) * t175 + rSges(6,2) * t174 + rSges(6,3) * t385;
t412 = pkin(4) * t391 + pkin(7) * t385 + t97;
t394 = qJ(4) * t247;
t390 = t244 * qJ(4);
t323 = rSges(6,1) * t245 - rSges(6,2) * t242;
t160 = rSges(6,3) * t243 - t246 * t323;
t332 = pkin(4) * t246 - pkin(7) * t243;
t384 = t160 - t332;
t322 = pkin(2) * t246 + qJ(3) * t243;
t176 = t322 * t244;
t177 = pkin(2) * t385 + qJ(3) * t391;
t383 = t176 * t244 + t177 * t247;
t166 = qJD(2) * t322 - qJD(3) * t246;
t327 = rSges(4,1) * t246 + rSges(4,3) * t243;
t382 = -qJD(2) * t327 - t166;
t229 = pkin(3) * t385;
t194 = t229 - t390;
t381 = -t177 - t194;
t205 = pkin(2) * t243 - qJ(3) * t246;
t369 = qJD(1) * t244;
t180 = t205 * t369;
t349 = t243 * t369;
t380 = pkin(3) * t349 + t180;
t206 = -rSges(4,3) * t246 + t421;
t379 = -t205 - t206;
t362 = qJD(3) * t243;
t378 = qJ(3) * t343 + t247 * t362;
t367 = qJD(1) * t247;
t377 = rSges(4,2) * t367 + rSges(4,3) * t343;
t376 = rSges(3,2) * t349 + rSges(3,3) * t367;
t365 = qJD(2) * t244;
t346 = t243 * t365;
t375 = -pkin(3) * t346 - qJ(4) * t369;
t374 = pkin(1) * t247 + t237;
t373 = t240 + t241;
t137 = -Icges(5,3) * t244 + t247 * t294;
t372 = qJD(1) * t137;
t141 = Icges(3,3) * t244 + t247 * t297;
t371 = qJD(1) * t141;
t146 = Icges(4,2) * t244 + t247 * t301;
t370 = qJD(1) * t146;
t361 = qJD(4) * t244;
t217 = pkin(2) * t346;
t344 = t244 * t364;
t347 = t246 * t367;
t359 = t176 * t367 + t244 * (pkin(2) * t347 + t244 * t362 - t217 + (t243 * t367 + t344) * qJ(3)) + t247 * (-pkin(2) * t443 - qJ(3) * t349 + t378);
t357 = rSges(5,2) * t385;
t93 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t385;
t95 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t385;
t320 = t242 * t93 - t245 * t95;
t38 = Icges(6,5) * t84 + Icges(6,6) * t83 - Icges(6,3) * t443;
t40 = Icges(6,4) * t84 + Icges(6,2) * t83 - Icges(6,6) * t443;
t42 = Icges(6,1) * t84 + Icges(6,4) * t83 - Icges(6,5) * t443;
t91 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t385;
t11 = (-qJD(2) * t320 + t38) * t243 + (qJD(2) * t91 + t242 * t40 - t245 * t42 + (t242 * t95 + t245 * t93) * qJD(5)) * t246;
t16 = t100 * t385 + t107 * t174 + t114 * t175 - t135 * t443 + t142 * t83 + t149 * t84;
t356 = t11 / 0.2e1 + t16 / 0.2e1;
t92 = Icges(6,4) * t173 + Icges(6,2) * t172 + Icges(6,6) * t387;
t94 = Icges(6,1) * t173 + Icges(6,4) * t172 + Icges(6,5) * t387;
t321 = t242 * t92 - t245 * t94;
t255 = -t346 + t347;
t85 = t339 * t388 + (-t247 * t338 - t344) * t242;
t86 = t338 * t386 + (t242 * t339 + t245 * t364) * t244;
t39 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t255;
t41 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t255;
t43 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t255;
t90 = Icges(6,5) * t173 + Icges(6,6) * t172 + Icges(6,3) * t387;
t10 = (-qJD(2) * t321 + t39) * t243 + (qJD(2) * t90 + t242 * t41 - t245 * t43 + (t242 * t94 + t245 * t92) * qJD(5)) * t246;
t17 = t100 * t387 + t107 * t172 + t114 * t173 + t135 * t255 + t142 * t85 + t149 * t86;
t355 = t17 / 0.2e1 + t10 / 0.2e1;
t34 = t243 * t91 + t246 * t320;
t47 = t135 * t385 + t142 * t174 + t149 * t175;
t354 = -t34 / 0.2e1 - t47 / 0.2e1;
t33 = t243 * t90 + t246 * t321;
t46 = t135 * t387 + t142 * t172 + t149 * t173;
t353 = t46 / 0.2e1 + t33 / 0.2e1;
t234 = pkin(6) * t367;
t352 = t234 + t378;
t351 = rSges(5,1) * t343 + rSges(5,2) * t443;
t162 = rSges(4,1) * t385 + rSges(4,3) * t391 + t236;
t350 = -t427 + t431;
t342 = -t366 / 0.2e1;
t341 = -pkin(3) * t243 - t205;
t340 = t414 * t247;
t133 = t379 * t247;
t193 = pkin(3) * t387 + t394;
t337 = t193 * t244 + t194 * t247 + t383;
t336 = t374 + t177;
t326 = rSges(5,1) * t246 + t419;
t335 = t326 + t341;
t334 = -pkin(3) * t364 - t166;
t333 = rSges(6,1) * t86 + rSges(6,2) * t85;
t330 = t243 * t416 - pkin(1);
t329 = -t294 * qJD(2) / 0.2e1 + (t301 + t297) * t441;
t328 = rSges(3,1) * t246 - rSges(3,2) * t243;
t325 = rSges(5,1) * t243 - rSges(5,2) * t246;
t27 = t172 * t92 + t173 * t94 + t387 * t90;
t28 = t172 * t93 + t173 * t95 + t387 * t91;
t18 = t244 * t28 - t247 * t27;
t319 = t244 * t27 + t247 * t28;
t29 = t174 * t92 + t175 * t94 + t385 * t90;
t30 = t174 * t93 + t175 * t95 + t385 * t91;
t19 = t244 * t30 - t247 * t29;
t318 = t244 * t29 + t247 * t30;
t317 = t244 * t34 - t247 * t33;
t316 = t244 * t33 + t247 * t34;
t238 = t247 * pkin(6);
t251 = t243 * t424 + t246 * t350 - pkin(1);
t248 = t244 * t251 - t394;
t51 = t238 + t248 + t324;
t311 = t229 + t336;
t52 = t311 - t390 + t412;
t315 = t244 * t52 + t247 * t51;
t70 = t160 * t387 - t243 * t96;
t71 = -t160 * t385 + t243 * t97;
t314 = t244 * t71 + t247 * t70;
t252 = (rSges(5,2) + t431) * t246 + t330;
t76 = t244 * t252 + t238 + t340;
t223 = rSges(5,1) * t391;
t77 = t244 * t414 + t223 + t311 - t357;
t313 = t244 * t77 + t247 * t76;
t312 = t244 * t97 - t247 * t96;
t295 = -Icges(4,3) * t246 + t403;
t280 = t341 - t384;
t163 = rSges(3,1) * t385 + t440;
t278 = -qJD(2) * t325 + t334;
t277 = -pkin(1) - t328;
t276 = t244 * ((pkin(3) * t368 + qJD(4)) * t247 + t375) + t247 * (-pkin(3) * t443 - qJ(4) * t367 - t361) + t193 * t367 + t359;
t275 = t352 - t361;
t274 = -qJD(4) * t247 + t217 - t375;
t123 = t335 * t247;
t273 = -rSges(5,3) * t244 - t357;
t121 = (rSges(6,1) * t242 + rSges(6,2) * t245) * t360 + (rSges(6,3) * t246 + t243 * t323) * qJD(2);
t272 = -qJD(2) * t331 - t121 + t334;
t271 = qJD(2) * t207;
t260 = qJD(2) * (-Icges(4,4) * t243 + Icges(4,6) * t246);
t258 = qJD(2) * (-Icges(3,5) * t243 - Icges(3,6) * t246);
t257 = qJD(2) * t295;
t256 = qJD(2) * (Icges(5,5) * t246 + Icges(5,6) * t243);
t75 = t280 * t247;
t253 = t243 * t415 + t246 * t428 - pkin(1);
t250 = t253 * t244;
t218 = pkin(4) * t343;
t192 = t328 * qJD(2);
t161 = t223 + t273;
t159 = t244 * t328 - t418;
t158 = t244 * t327 - t420;
t157 = t244 * t325 + t417;
t132 = t379 * t244;
t128 = t163 + t374;
t127 = t244 * t277 + t238 + t418;
t122 = t335 * t244;
t111 = t244 * t260 + t370;
t110 = -qJD(1) * t145 + t247 * t260;
t106 = t244 * t258 + t371;
t105 = -qJD(1) * t140 + t247 * t258;
t102 = t244 * t256 + t372;
t101 = -qJD(1) * t136 + t247 * t256;
t99 = t336 + t162;
t98 = t238 + t250 + t420;
t79 = t207 * t365 + ((-rSges(3,3) - pkin(6)) * t244 + t277 * t247) * qJD(1);
t78 = -rSges(3,1) * t443 - rSges(3,2) * t343 - pkin(1) * t369 + t234 + t376;
t74 = t280 * t244;
t73 = qJD(1) * t133 + t244 * t382;
t72 = t206 * t369 + t247 * t382 + t180;
t69 = t141 * t244 - t247 * t281;
t68 = t140 * t244 - t266;
t67 = t146 * t244 + t247 * t285;
t66 = t145 * t244 + t270;
t65 = -t137 * t244 - t247 * t283;
t64 = -t136 * t244 - t268;
t62 = -t141 * t247 - t265;
t61 = -t140 * t247 - t244 * t282;
t60 = -t146 * t247 + t269;
t59 = -t145 * t247 + t244 * t286;
t58 = t137 * t247 - t267;
t57 = t136 * t247 - t244 * t284;
t55 = qJD(1) * t123 + t244 * t278;
t54 = t247 * t278 - t326 * t369 + t380;
t53 = t158 * t244 + t162 * t247 + t383;
t50 = t217 + (-t362 + (t246 * t415 + t421) * qJD(2)) * t244 + ((-rSges(4,2) - pkin(6)) * t244 + t253 * t247) * qJD(1);
t49 = qJD(1) * t250 + t345 * t428 + t352 + t377;
t48 = t312 * t246;
t45 = rSges(6,3) * t255 + t333;
t44 = -rSges(6,3) * t443 + t422;
t37 = t157 * t244 + t161 * t247 + t337;
t36 = (-t362 + (t246 * t416 - t419) * qJD(2)) * t244 + ((rSges(5,3) - pkin(6)) * t244 + t252 * t247) * qJD(1) + t274;
t35 = t431 * t345 + (t340 + (t246 * t431 + t330) * t244) * qJD(1) + t275 + t351;
t32 = qJD(1) * t75 + t244 * t272;
t31 = t247 * t272 + t369 * t384 + t380;
t26 = t244 * t413 + t247 * t412 + t337;
t25 = (-t362 + (t243 * t427 + t246 * t424) * qJD(2)) * t244 + (t247 * t251 - t237) * qJD(1) + t274 - t333;
t24 = qJD(1) * t248 + t345 * t350 + t218 + t275 + t422;
t23 = (-t160 * t365 - t45) * t243 + (-qJD(2) * t96 + t121 * t244 + t160 * t367) * t246;
t22 = (t160 * t363 + t44) * t243 + (qJD(2) * t97 - t121 * t247 + t160 * t369) * t246;
t21 = t247 * t377 + (-t206 * t240 - t241 * t421) * qJD(2) + (t247 * t158 + (-t162 - t177 + t236) * t244) * qJD(1) + t359;
t15 = ((t157 - t417) * qJD(1) + t351) * t247 + (t326 * t365 + (-t161 + t273 + t381) * qJD(1)) * t244 + t276;
t14 = t312 * t366 + (-t244 * t44 + t247 * t45 + (-t244 * t96 - t247 * t97) * qJD(1)) * t246;
t13 = t243 * t47 + t246 * t318;
t12 = t243 * t46 + t246 * t319;
t9 = t91 * t347 + t172 * t40 + t173 * t42 + t85 * t93 + t86 * t95 + (t246 * t38 - t366 * t91) * t244;
t8 = t90 * t347 + t172 * t41 + t173 * t43 + t85 * t92 + t86 * t94 + (t246 * t39 - t366 * t90) * t244;
t7 = -t91 * t345 + t174 * t40 + t175 * t42 + t83 * t93 + t84 * t95 + (t247 * t38 - t369 * t91) * t246;
t6 = -t90 * t345 + t174 * t41 + t175 * t43 + t83 * t92 + t84 * t94 + (t247 * t39 - t369 * t90) * t246;
t5 = (-pkin(7) * t345 + t218 + t44) * t247 + (t332 * t365 + t45) * t244 + (t413 * t247 + (t381 - t412) * t244) * qJD(1) + t276;
t4 = qJD(1) * t319 + t244 * t9 - t247 * t8;
t3 = qJD(1) * t318 + t244 * t7 - t247 * t6;
t2 = (-qJD(2) * t319 + t17) * t243 + (-qJD(1) * t18 + qJD(2) * t46 + t244 * t8 + t247 * t9) * t246;
t1 = (-qJD(2) * t318 + t16) * t243 + (-qJD(1) * t19 + qJD(2) * t47 + t244 * t6 + t247 * t7) * t246;
t20 = [(t127 * t79 + t128 * t78) * t438 + (t49 * t99 + t50 * t98) * t437 + (t35 * t77 + t36 * t76) * t436 + (t24 * t52 + t25 * t51) * t435 + t249 + (t295 + t308 + t310 - t299 - t451) * t366 + (t303 - t305 - t296 + t450) * t364; m(4) * (t132 * t49 + t133 * t50 + t72 * t98 + t73 * t99) + m(5) * (t122 * t35 + t123 * t36 + t54 * t76 + t55 * t77) + m(6) * (t24 * t74 + t25 * t75 + t31 * t51 + t32 * t52) + (m(3) * (-t127 * t192 - t207 * t79) + t329 * t247 + (t148 * t444 + t257 * t447 + t449 * t446 + (t139 + t151) * t425) * t246 - t355) * t247 + (m(3) * (-t128 * t192 - t207 * t78) + t329 * t244 + (t147 * t444 + t257 * t429 + t449 * t445 + (t138 + t150) * t425) * t246 + t356) * t244 + ((t143 * t244 + t144 * t247) * t425 + ((t153 + t155) * t247 + (t152 + t154) * t244) * t444 + (t244 * t445 + t247 * t446) * t450 * qJD(2)) * t243 + (-t267 / 0.2e1 + t269 / 0.2e1 - t265 / 0.2e1 + t268 / 0.2e1 - t270 / 0.2e1 + t266 / 0.2e1) * qJD(2) + ((-t128 * t426 + (-t151 / 0.2e1 - t139 / 0.2e1 + t148 / 0.2e1) * t246 + (-t144 / 0.2e1 + t153 / 0.2e1 + t155 / 0.2e1) * t243 - t354) * t247 + (t127 * t426 + (-t138 / 0.2e1 + t147 / 0.2e1 - t150 / 0.2e1) * t246 + (t152 / 0.2e1 + t154 / 0.2e1 - t143 / 0.2e1) * t243 + t353) * t244) * qJD(1); t244 * t3 - t247 * t4 + (t26 * t5 + t31 * t75 + t32 * t74) * t435 + (t122 * t55 + t123 * t54 + t15 * t37) * t436 + (t132 * t73 + t133 * t72 + t21 * t53) * t437 - t247 * ((-t247 * t102 + (t58 + t268) * qJD(1)) * t247 + (t57 * qJD(1) + (t144 * t366 + t151 * t364 - t372) * t244 + (t101 + (-t143 * t243 - t150 * t246) * qJD(2) - t283 * qJD(1)) * t247) * t244) - t247 * ((t247 * t106 + (t62 + t266) * qJD(1)) * t247 + (t61 * qJD(1) + (-t148 * t364 - t155 * t366 + t371) * t244 + (-t105 + (t147 * t246 + t154 * t243) * qJD(2) - t281 * qJD(1)) * t247) * t244) - t247 * ((t247 * t111 + (t60 - t270) * qJD(1)) * t247 + (t59 * qJD(1) + (t139 * t364 - t153 * t366 + t370) * t244 + (-t110 + (-t138 * t246 + t152 * t243) * qJD(2) + t285 * qJD(1)) * t247) * t244) + t244 * ((-t244 * t101 + (t64 + t267) * qJD(1)) * t244 + (t65 * qJD(1) + (-t143 * t366 - t150 * t364) * t247 + (t102 + (t144 * t243 + t151 * t246) * qJD(2) + (-t137 - t284) * qJD(1)) * t244) * t247) + t244 * ((t244 * t105 + (t68 + t265) * qJD(1)) * t244 + (t69 * qJD(1) + (t147 * t364 + t154 * t366) * t247 + (-t106 + (-t148 * t246 - t155 * t243) * qJD(2) + (t141 - t282) * qJD(1)) * t244) * t247) + t244 * ((t244 * t110 + (t66 - t269) * qJD(1)) * t244 + (t67 * qJD(1) + (-t138 * t364 + t152 * t366) * t247 + (-t111 + (t139 * t246 - t153 * t243) * qJD(2) + (t146 + t286) * qJD(1)) * t244) * t247) + ((t159 * t244 + t163 * t247) * ((qJD(1) * t159 - t247 * t271 + t376) * t247 + (-t244 * t271 + (-t163 + t440) * qJD(1)) * t244) + t373 * t207 * t192) * t438 + (t18 + (-t57 - t59 - t61) * t247 + (t58 + t60 + t62) * t244) * t369 + (t19 + (-t64 - t66 - t68) * t247 + (t65 + t67 + t69) * t244) * t367; 0.2e1 * (t315 * t432 + t313 * t433 + (t244 * t99 + t247 * t98) * t434) * t364 + 0.2e1 * ((t24 * t244 + t247 * t25 + t367 * t52 - t369 * t51) * t432 + (t244 * t35 + t247 * t36 + t367 * t77 - t369 * t76) * t433 + (t244 * t49 + t247 * t50 + t367 * t99 - t369 * t98) * t434) * t243; 0.2e1 * ((t363 * t75 + t365 * t74 - t5) * t432 + (t122 * t365 + t123 * t363 - t15) * t433 + (t132 * t365 + t133 * t363 - t21) * t434) * t246 + 0.2e1 * ((qJD(2) * t26 + t244 * t32 + t247 * t31 + t367 * t74 - t369 * t75) * t432 + (qJD(2) * t37 + t122 * t367 - t123 * t369 + t244 * t55 + t247 * t54) * t433 + (qJD(2) * t53 + t132 * t367 - t133 * t369 + t244 * t73 + t247 * t72) * t434) * t243; 0.4e1 * (t434 + t433 + t432) * (-0.1e1 + t373) * t243 * t364; m(6) * (-qJD(1) * t315 + t24 * t247 - t244 * t25) + m(5) * (-qJD(1) * t313 - t244 * t36 + t247 * t35); m(6) * (-t244 * t31 + t247 * t32 + (-t244 * t74 - t247 * t75) * qJD(1)) + m(5) * (-t244 * t54 + t247 * t55 + (-t122 * t244 - t123 * t247) * qJD(1)); 0; 0; m(6) * (t22 * t52 + t23 * t51 + t24 * t71 + t25 * t70) + (-t244 * t353 + t247 * t354) * t366 + (t356 * t247 + t355 * t244 + (t244 * t354 + t247 * t353) * qJD(1)) * t246 + t423; m(6) * (t14 * t26 + t22 * t74 + t23 * t75 + t31 * t70 + t32 * t71 - t48 * t5) + (t19 * t342 + (qJD(1) * t34 - t10) * t430 - t2 / 0.2e1 + t13 * t425) * t247 + ((qJD(1) * t33 + t11) * t430 + t18 * t342 + t12 * t425 + t1 / 0.2e1) * t244 + (t3 * t429 + t317 * t441 + t4 * t446 + (t18 * t429 + t19 * t447) * qJD(1)) * t246; m(6) * ((qJD(2) * t314 - t14) * t246 + (-qJD(2) * t48 + t22 * t244 + t23 * t247 + (-t244 * t70 + t247 * t71) * qJD(1)) * t243); m(6) * (-qJD(1) * t314 + t22 * t247 - t23 * t244); (-t14 * t48 + t22 * t71 + t23 * t70) * t435 + ((-t12 * t244 - t13 * t247 - t243 * t316) * qJD(2) + t423) * t243 + (t247 * t1 + t244 * t2 + t243 * (t10 * t244 + t11 * t247) + (t243 * t63 + t246 * t316) * qJD(2) + (t12 * t247 - t13 * t244 - t243 * t317) * qJD(1)) * t246;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t20(1), t20(2), t20(4), t20(7), t20(11); t20(2), t20(3), t20(5), t20(8), t20(12); t20(4), t20(5), t20(6), t20(9), t20(13); t20(7), t20(8), t20(9), t20(10), t20(14); t20(11), t20(12), t20(13), t20(14), t20(15);];
Mq = res;
