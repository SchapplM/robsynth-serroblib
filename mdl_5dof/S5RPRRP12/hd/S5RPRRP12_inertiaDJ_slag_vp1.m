% Calculate time derivative of joint inertia matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP12_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP12_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP12_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:17
% EndTime: 2019-12-31 18:56:36
% DurationCPUTime: 11.47s
% Computational Cost: add. (10216->705), mult. (26647->988), div. (0->0), fcn. (25314->6), ass. (0->342)
t244 = sin(qJ(4));
t245 = sin(qJ(3));
t248 = cos(qJ(3));
t247 = cos(qJ(4));
t389 = Icges(6,4) * t247;
t284 = -Icges(6,2) * t244 + t389;
t171 = Icges(6,6) * t245 + t248 * t284;
t391 = Icges(5,4) * t247;
t285 = -Icges(5,2) * t244 + t391;
t172 = Icges(5,6) * t245 + t248 * t285;
t453 = -t172 - t171;
t455 = t453 * t244;
t280 = Icges(6,5) * t247 - Icges(6,6) * t244;
t167 = Icges(6,3) * t245 + t248 * t280;
t281 = Icges(5,5) * t247 - Icges(5,6) * t244;
t168 = Icges(5,3) * t245 + t248 * t281;
t454 = t167 + t168;
t390 = Icges(6,4) * t244;
t288 = Icges(6,1) * t247 - t390;
t175 = Icges(6,5) * t245 + t248 * t288;
t392 = Icges(5,4) * t244;
t289 = Icges(5,1) * t247 - t392;
t176 = Icges(5,5) * t245 + t248 * t289;
t449 = -t176 - t175;
t243 = -qJ(5) - pkin(7);
t395 = rSges(6,3) - t243;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t364 = t247 * t249;
t368 = t246 * t244;
t190 = -t245 * t368 + t364;
t367 = t246 * t247;
t372 = t244 * t249;
t191 = t245 * t367 + t372;
t366 = t246 * t248;
t124 = Icges(5,5) * t191 + Icges(5,6) * t190 - Icges(5,3) * t366;
t128 = Icges(5,4) * t191 + Icges(5,2) * t190 - Icges(5,6) * t366;
t132 = Icges(5,1) * t191 + Icges(5,4) * t190 - Icges(5,5) * t366;
t53 = -t124 * t366 + t128 * t190 + t132 * t191;
t370 = t245 * t249;
t192 = t244 * t370 + t367;
t193 = -t245 * t364 + t368;
t363 = t248 * t249;
t125 = Icges(5,5) * t193 + Icges(5,6) * t192 + Icges(5,3) * t363;
t129 = Icges(5,4) * t193 + Icges(5,2) * t192 + Icges(5,6) * t363;
t133 = Icges(5,1) * t193 + Icges(5,4) * t192 + Icges(5,5) * t363;
t54 = -t125 * t366 + t129 * t190 + t133 * t191;
t296 = t246 * t53 - t249 * t54;
t122 = Icges(6,5) * t191 + Icges(6,6) * t190 - Icges(6,3) * t366;
t126 = Icges(6,4) * t191 + Icges(6,2) * t190 - Icges(6,6) * t366;
t130 = Icges(6,1) * t191 + Icges(6,4) * t190 - Icges(6,5) * t366;
t51 = -t122 * t366 + t126 * t190 + t130 * t191;
t123 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t363;
t127 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t363;
t131 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t363;
t52 = -t123 * t366 + t127 * t190 + t131 * t191;
t297 = t246 * t51 - t249 * t52;
t84 = -t167 * t366 + t171 * t190 + t175 * t191;
t85 = -t168 * t366 + t172 * t190 + t176 * t191;
t452 = (-t296 - t297) * t248 + (t84 + t85) * t245;
t57 = t124 * t363 + t192 * t128 + t193 * t132;
t58 = t125 * t363 + t192 * t129 + t193 * t133;
t294 = t246 * t57 - t249 * t58;
t55 = t122 * t363 + t192 * t126 + t193 * t130;
t56 = t123 * t363 + t192 * t127 + t193 * t131;
t295 = t246 * t55 - t249 * t56;
t86 = t167 * t363 + t192 * t171 + t193 * t175;
t87 = t168 * t363 + t192 * t172 + t193 * t176;
t451 = (-t294 - t295) * t248 + (t86 + t87) * t245;
t450 = (-t449 * t247 + t455) * t248 + t454 * t245;
t345 = qJD(4) * t248;
t140 = (-Icges(6,5) * t244 - Icges(6,6) * t247) * t345 + (Icges(6,3) * t248 - t245 * t280) * qJD(3);
t141 = (-Icges(5,5) * t244 - Icges(5,6) * t247) * t345 + (Icges(5,3) * t248 - t245 * t281) * qJD(3);
t148 = (-Icges(6,1) * t244 - t389) * t345 + (Icges(6,5) * t248 - t245 * t288) * qJD(3);
t149 = (-Icges(5,1) * t244 - t391) * t345 + (Icges(5,5) * t248 - t245 * t289) * qJD(3);
t347 = qJD(3) * t248;
t349 = qJD(3) * t245;
t448 = t454 * t347 - t349 * t455 + (t140 + t141) * t245 + ((t148 + t149) * t248 + t449 * t349) * t247;
t322 = t246 * t347;
t350 = qJD(1) * t249;
t447 = t245 * t350 + t322;
t346 = qJD(3) * t249;
t323 = t245 * t346;
t351 = qJD(1) * t246;
t446 = t248 * t351 + t323;
t348 = qJD(3) * t246;
t325 = t245 * t348;
t327 = t248 * t350;
t257 = t325 - t327;
t344 = qJD(5) * t248;
t310 = qJD(4) * t245 + qJD(1);
t434 = t310 * t244;
t445 = pkin(4) * t434 + t344;
t144 = (-Icges(6,2) * t247 - t390) * t345 + (Icges(6,6) * t248 - t245 * t284) * qJD(3);
t145 = (-Icges(5,2) * t247 - t392) * t345 + (Icges(5,6) * t248 - t245 * t285) * qJD(3);
t444 = (-t145 - t144) * t244;
t443 = t453 * t247;
t277 = t127 * t244 - t131 * t247;
t62 = t123 * t245 - t248 * t277;
t275 = t129 * t244 - t133 * t247;
t64 = t125 * t245 - t248 * t275;
t402 = t62 + t64;
t278 = t126 * t244 - t130 * t247;
t61 = t122 * t245 - t248 * t278;
t276 = t128 * t244 - t132 * t247;
t63 = t124 * t245 - t248 * t276;
t403 = t61 + t63;
t442 = -t246 * t403 + t249 * t402;
t441 = t246 * t402 + t249 * t403;
t309 = qJD(1) * t245 + qJD(4);
t117 = -t310 * t367 + (-t249 * t309 - t322) * t244;
t118 = t309 * t364 + (t247 * t347 - t434) * t246;
t69 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t257;
t73 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t257;
t77 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t257;
t19 = (qJD(3) * t278 + t69) * t245 + (qJD(3) * t122 - t244 * t73 + t247 * t77 + (-t126 * t247 - t130 * t244) * qJD(4)) * t248;
t71 = Icges(5,5) * t118 + Icges(5,6) * t117 + Icges(5,3) * t257;
t75 = Icges(5,4) * t118 + Icges(5,2) * t117 + Icges(5,6) * t257;
t79 = Icges(5,1) * t118 + Icges(5,4) * t117 + Icges(5,5) * t257;
t21 = (qJD(3) * t276 + t71) * t245 + (qJD(3) * t124 - t244 * t75 + t247 * t79 + (-t128 * t247 - t132 * t244) * qJD(4)) * t248;
t440 = t19 + t21;
t321 = t248 * t346;
t431 = t246 * t309 - t321;
t115 = -t431 * t244 + t310 * t364;
t116 = t431 * t247 + t249 * t434;
t68 = Icges(6,5) * t116 + Icges(6,6) * t115 - Icges(6,3) * t446;
t72 = Icges(6,4) * t116 + Icges(6,2) * t115 - Icges(6,6) * t446;
t76 = Icges(6,1) * t116 + Icges(6,4) * t115 - Icges(6,5) * t446;
t20 = (qJD(3) * t277 + t68) * t245 + (qJD(3) * t123 - t244 * t72 + t247 * t76 + (-t127 * t247 - t131 * t244) * qJD(4)) * t248;
t70 = Icges(5,5) * t116 + Icges(5,6) * t115 - Icges(5,3) * t446;
t74 = Icges(5,4) * t116 + Icges(5,2) * t115 - Icges(5,6) * t446;
t78 = Icges(5,1) * t116 + Icges(5,4) * t115 - Icges(5,5) * t446;
t22 = (qJD(3) * t275 + t70) * t245 + (qJD(3) * t125 - t244 * t74 + t247 * t78 + (-t129 * t247 - t133 * t244) * qJD(4)) * t248;
t439 = t20 + t22;
t233 = pkin(4) * t247 + pkin(3);
t397 = pkin(4) * qJD(4);
t340 = t247 * t397;
t438 = t118 * rSges(6,1) + t117 * rSges(6,2) + t257 * rSges(6,3) + t447 * t233 + t243 * t327 + t249 * t340;
t407 = pkin(3) - t233;
t437 = t245 * t407;
t393 = Icges(4,4) * t248;
t290 = Icges(4,1) * t245 + t393;
t177 = Icges(4,5) * t249 + t246 * t290;
t379 = t177 * t245;
t394 = Icges(4,4) * t245;
t286 = Icges(4,2) * t248 + t394;
t173 = Icges(4,6) * t249 + t246 * t286;
t382 = t173 * t248;
t271 = t379 + t382;
t436 = t249 * t271;
t307 = rSges(4,1) * t245 + rSges(4,2) * t248;
t261 = t249 * t307;
t231 = pkin(3) * t370;
t302 = -t193 * rSges(6,1) - t192 * rSges(6,2);
t406 = -pkin(7) - t243;
t317 = t406 * t248;
t375 = t233 * t245;
t360 = pkin(4) * t368 + t231 + (t317 - t375) * t249 + rSges(6,3) * t363 - t302;
t435 = t249 * t360;
t371 = t245 * t246;
t433 = t191 * rSges(6,1) + t190 * rSges(6,2) + pkin(4) * t372 + t233 * t371 - t395 * t366;
t334 = t447 * rSges(4,1) + rSges(4,2) * t327;
t419 = -pkin(1) - pkin(6);
t343 = -rSges(4,3) + t419;
t355 = qJ(2) * t350 + qJD(2) * t246;
t107 = (-rSges(4,2) * t349 + qJD(1) * t343) * t246 + t334 + t355;
t399 = rSges(4,2) * t245;
t212 = rSges(4,1) * t248 - t399;
t235 = qJD(2) * t249;
t108 = t235 + t212 * t346 + (t343 * t249 + (-qJ(2) - t307) * t246) * qJD(1);
t432 = -t107 * t249 + t246 * t108;
t430 = -t248 * t395 + t375;
t282 = Icges(4,5) * t245 + Icges(4,6) * t248;
t429 = -Icges(4,3) * t246 + t249 * t282;
t428 = -Icges(4,6) * t246 + t249 * t286;
t427 = -Icges(4,5) * t246 + t249 * t290;
t213 = pkin(3) * t248 + pkin(7) * t245;
t204 = t246 * t213;
t301 = rSges(6,1) * t247 - rSges(6,2) * t244;
t358 = (t301 - t407) * t248 + (rSges(6,3) + t406) * t245;
t314 = t358 * t246;
t105 = t204 + t314;
t106 = (-t213 - t358) * t249;
t197 = t213 * t351;
t411 = pkin(3) * t245;
t206 = (pkin(7) * t248 - t411) * qJD(3);
t313 = qJD(1) * t358;
t320 = t244 * t345;
t362 = -pkin(4) * t320 + qJD(5) * t245 + (-rSges(6,1) * t244 - rSges(6,2) * t247) * t345 + (rSges(6,3) * t248 - t245 * t301 + t317 + t437) * qJD(3);
t47 = t197 + t246 * t313 + (-t206 - t362) * t249;
t250 = t246 * t362 + t249 * t313;
t356 = t246 * t206 + t213 * t350;
t48 = t250 + t356;
t426 = qJD(1) * (t105 * t249 + t106 * t246) + t48 * t246 - t249 * t47;
t311 = -pkin(4) * t244 + t419;
t253 = qJD(1) * t311 - t344;
t341 = t244 * t397;
t44 = ((-qJD(3) * t243 - t341) * t245 + t253) * t246 + t355 + t438;
t303 = t116 * rSges(6,1) + t115 * rSges(6,2);
t374 = t233 * t248;
t45 = t235 + (-t340 + (-qJ(2) - t430) * qJD(1)) * t246 + (-t245 * t341 + (t245 * t395 + t374) * qJD(3) + t253) * t249 - t303;
t237 = t249 * qJ(2);
t91 = t311 * t246 + t430 * t249 + t237 + t302;
t354 = t249 * pkin(1) + t246 * qJ(2);
t331 = t249 * pkin(6) + t354;
t92 = t331 + t433;
t425 = qJD(1) * (t246 * t92 + t249 * t91) + t246 * t45 - t249 * t44;
t312 = qJD(3) * t358;
t332 = pkin(3) * t321 + t446 * pkin(7);
t401 = (t340 + (t243 * t248 - t437) * qJD(1)) * t246 + ((t243 * t245 - t374) * qJD(3) + t445) * t249 + t332 - rSges(6,3) * t446 + t303;
t23 = (-t249 * t312 - t401) * t245 + (-qJD(3) * t360 + t249 * t362 - t351 * t358) * t248;
t229 = pkin(3) * t371;
t195 = -pkin(7) * t366 + t229;
t361 = -t195 + t433;
t333 = t447 * pkin(3) + pkin(7) * t325;
t255 = pkin(7) * t327 - t333;
t400 = -(-t243 * t349 - t445) * t246 - t255 - t438;
t24 = (-t246 * t312 - t400) * t245 + (qJD(3) * t361 + t250) * t248;
t59 = t245 * t361 + t248 * t314;
t60 = -t245 * t360 + t358 * t363;
t424 = qJD(1) * (t246 * t59 + t249 * t60) + t23 * t246 - t24 * t249;
t423 = t246 * t360 + t249 * t361;
t422 = 2 * m(4);
t421 = 2 * m(5);
t420 = 2 * m(6);
t241 = t246 ^ 2;
t242 = t249 ^ 2;
t418 = -t245 / 0.2e1;
t415 = t248 / 0.2e1;
t413 = rSges(3,2) - pkin(1);
t412 = m(4) * t212;
t398 = rSges(5,3) * t248;
t396 = t246 * rSges(4,3);
t238 = t249 * rSges(4,3);
t305 = -t193 * rSges(5,1) - t192 * rSges(5,2);
t139 = rSges(5,3) * t363 - t305;
t384 = t139 * t249;
t383 = t173 * t245;
t381 = t428 * t245;
t380 = t428 * t248;
t378 = t177 * t248;
t377 = t427 * t245;
t376 = t427 * t248;
t357 = t191 * rSges(5,1) + t190 * rSges(5,2);
t137 = -rSges(5,3) * t366 + t357;
t359 = -t137 - t195;
t353 = t241 + t242;
t169 = Icges(4,3) * t249 + t246 * t282;
t352 = qJD(1) * t169;
t342 = m(6) * t349;
t35 = t52 * t246 + t249 * t51;
t36 = t54 * t246 + t249 * t53;
t339 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t56 * t246 + t249 * t55;
t38 = t58 * t246 + t249 * t57;
t338 = -t37 / 0.2e1 - t38 / 0.2e1;
t337 = t118 * rSges(5,1) + t117 * rSges(5,2) + rSges(5,3) * t325;
t336 = -t195 - t361;
t179 = rSges(4,1) * t371 + rSges(4,2) * t366 + t238;
t330 = t248 * (-rSges(5,3) - pkin(7));
t304 = rSges(5,1) * t247 - rSges(5,2) * t244;
t181 = rSges(5,3) * t245 + t248 * t304;
t329 = t181 * t351;
t316 = t450 * t347 + ((t444 + (t449 * t244 + t443) * qJD(4)) * t248 + t448) * t245;
t205 = t307 * qJD(3);
t315 = t205 * t353;
t306 = t116 * rSges(5,1) + t115 * rSges(5,2);
t291 = Icges(4,1) * t248 - t394;
t287 = -Icges(4,2) * t245 + t393;
t283 = Icges(4,5) * t248 - Icges(4,6) * t245;
t274 = t137 * t249 + t139 * t246;
t270 = -t377 - t380;
t33 = t117 * t171 + t118 * t175 - t140 * t366 + t190 * t144 + t191 * t148 + t167 * t257;
t34 = t117 * t172 + t118 * t176 - t141 * t366 + t190 * t145 + t191 * t149 + t168 * t257;
t267 = -t19 / 0.2e1 - t21 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1;
t31 = t115 * t171 + t116 * t175 + t140 * t363 + t192 * t144 + t193 * t148 - t167 * t446;
t32 = t115 * t172 + t116 * t176 + t141 * t363 + t192 * t145 + t193 * t149 - t168 * t446;
t266 = t20 / 0.2e1 + t22 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t265 = t61 / 0.2e1 + t63 / 0.2e1 + t84 / 0.2e1 + t85 / 0.2e1;
t264 = -t62 / 0.2e1 - t64 / 0.2e1 - t86 / 0.2e1 - t87 / 0.2e1;
t263 = rSges(3,3) * t249 + t246 * t413;
t153 = (-rSges(5,1) * t244 - rSges(5,2) * t247) * t345 + (-t245 * t304 + t398) * qJD(3);
t262 = t246 * t153 + t181 * t350;
t260 = t270 * t246;
t259 = qJD(3) * t291;
t258 = qJD(3) * t287;
t254 = t246 * t419 + t249 * t330;
t196 = pkin(7) * t363 - t231;
t185 = t249 * t196;
t184 = -rSges(3,2) * t249 + t246 * rSges(3,3) + t354;
t183 = t237 + t263;
t182 = t396 - t261;
t161 = t235 + (t413 * t249 + (-rSges(3,3) - qJ(2)) * t246) * qJD(1);
t160 = qJD(1) * t263 + t355;
t159 = t331 + t179;
t158 = t246 * t343 + t237 + t261;
t156 = (-t181 - t213) * t249;
t155 = t181 * t246 + t204;
t154 = t249 * (qJD(1) * t229 - t332);
t143 = t429 * qJD(1) + t283 * t348;
t142 = -t283 * t346 + t352;
t104 = t246 * t330 + t229 + t331 + t357;
t103 = t231 + t237 + t254 + t305;
t102 = -t245 * t139 + t181 * t363;
t101 = t137 * t245 + t181 * t366;
t100 = -t246 * t429 - t249 * t270;
t99 = t246 * t169 - t436;
t96 = -t249 * t429 + t260;
t95 = t169 * t249 + t246 * t271;
t90 = t274 * t248;
t89 = t262 + t356;
t88 = t329 + t197 + (-t153 - t206) * t249;
t83 = -rSges(5,3) * t327 + t337;
t81 = -rSges(5,3) * t446 + t306;
t65 = t246 * t359 + t185 + t384;
t50 = rSges(5,3) * t323 + t235 + (t419 * t249 + (-qJ(2) + t398 - t411) * t246) * qJD(1) - t306 + t332;
t49 = qJD(1) * t254 + t333 + t337 + t355;
t46 = t423 * t248;
t43 = t246 * t336 + t185 + t435;
t42 = (-t181 * t348 + t83) * t245 + (qJD(3) * t137 + t262) * t248;
t41 = (-t181 * t346 - t81) * t245 + (-qJD(3) * t139 + t153 * t249 - t329) * t248;
t30 = t274 * t349 + (-t246 * t81 - t249 * t83 + (t246 * t137 - t384) * qJD(1)) * t248;
t29 = t249 * t81 + t154 + (t255 - t83) * t246 + (t359 * t249 + (-t139 - t196) * t246) * qJD(1);
t18 = t117 * t129 + t118 * t133 + t125 * t257 + t190 * t74 + t191 * t78 - t366 * t70;
t17 = t117 * t128 + t118 * t132 + t124 * t257 + t190 * t75 + t191 * t79 - t366 * t71;
t16 = t117 * t127 + t118 * t131 + t123 * t257 + t190 * t72 + t191 * t76 - t366 * t68;
t15 = t117 * t126 + t118 * t130 + t122 * t257 + t190 * t73 + t191 * t77 - t366 * t69;
t14 = t115 * t129 + t116 * t133 - t125 * t446 + t192 * t74 + t193 * t78 + t363 * t70;
t13 = t115 * t128 + t116 * t132 - t124 * t446 + t192 * t75 + t193 * t79 + t363 * t71;
t12 = t115 * t127 + t116 * t131 - t123 * t446 + t192 * t72 + t193 * t76 + t363 * t68;
t11 = t115 * t126 + t116 * t130 - t122 * t446 + t192 * t73 + t193 * t77 + t363 * t69;
t10 = t154 + t401 * t249 + (t255 + t400) * t246 + (t336 * t249 + (-t196 - t360) * t246) * qJD(1);
t9 = t423 * t349 + (t400 * t249 - t401 * t246 + (t246 * t361 - t435) * qJD(1)) * t248;
t8 = -qJD(1) * t296 + t17 * t249 + t18 * t246;
t7 = -qJD(1) * t297 + t15 * t249 + t16 * t246;
t6 = -qJD(1) * t294 + t13 * t249 + t14 * t246;
t5 = -qJD(1) * t295 + t11 * t249 + t12 * t246;
t4 = (qJD(3) * t296 + t34) * t245 + (-qJD(1) * t36 + qJD(3) * t85 - t17 * t246 + t18 * t249) * t248;
t3 = (qJD(3) * t297 + t33) * t245 + (-qJD(1) * t35 + qJD(3) * t84 - t15 * t246 + t16 * t249) * t248;
t2 = (qJD(3) * t294 + t32) * t245 + (-qJD(1) * t38 + qJD(3) * t87 - t13 * t246 + t14 * t249) * t248;
t1 = (qJD(3) * t295 + t31) * t245 + (-qJD(1) * t37 + qJD(3) * t86 - t11 * t246 + t12 * t249) * t248;
t25 = [0.2e1 * m(3) * (t160 * t184 + t161 * t183) + (t107 * t159 + t108 * t158) * t422 + (t103 * t50 + t104 * t49) * t421 + (t44 * t92 + t45 * t91) * t420 + t286 * t349 - t290 * t347 - t245 * t259 + t449 * t320 + t345 * t443 + (-t258 + t444) * t248 + t448; m(6) * t425 + m(5) * (t246 * t50 - t249 * t49 + (t103 * t249 + t104 * t246) * qJD(1)) + m(4) * ((t158 * t249 + t159 * t246) * qJD(1) + t432) + m(3) * (-t160 * t249 + t246 * t161 + (t183 * t249 + t184 * t246) * qJD(1)); 0; ((t428 * qJD(1) + t246 * t258) * t418 + (t427 * qJD(1) + t246 * t259) * t415 + (-t382 / 0.2e1 - t379 / 0.2e1) * qJD(3) - t267) * t249 + ((qJD(1) * t173 - t287 * t346) * t418 + (qJD(1) * t177 - t291 * t346) * t415 + (t380 / 0.2e1 + t377 / 0.2e1) * qJD(3) + t266) * t246 + m(5) * (t103 * t89 + t104 * t88 + t155 * t50 + t156 * t49) + m(6) * (t105 * t45 + t106 * t44 + t47 * t92 + t48 * t91) + m(4) * (t432 * t212 - (t158 * t246 - t159 * t249) * t205) - (t241 / 0.2e1 + t242 / 0.2e1) * t282 * qJD(3) + ((t383 / 0.2e1 - t378 / 0.2e1 + t159 * t412 - t265) * t246 + (t381 / 0.2e1 - t376 / 0.2e1 + t158 * t412 - t264) * t249) * qJD(1); m(5) * (t89 * t246 - t249 * t88 + (t155 * t249 + t156 * t246) * qJD(1)) + m(6) * t426 - m(4) * t315; (t43 * t10 + t105 * t48 + t106 * t47) * t420 + t246 * t6 + t246 * t5 + t249 * t7 + t249 * t8 + (t155 * t89 + t156 * t88 + t65 * t29) * t421 + ((-t246 * t179 + t182 * t249) * (-t246 * t334 + (-t212 * t242 + t241 * t399) * qJD(3) + ((-t179 + t238) * t249 + (-t182 + t261 + t396) * t246) * qJD(1)) - t212 * t315) * t422 + t249 * ((t249 * t143 + (t96 + t436) * qJD(1)) * t249 + (-t95 * qJD(1) + (-t347 * t427 + t349 * t428) * t246 + (t142 + (t378 - t383) * qJD(3) + (-t169 + t270) * qJD(1)) * t249) * t246) + t246 * ((t246 * t142 + (-t99 + t260) * qJD(1)) * t246 + (t100 * qJD(1) + (t173 * t349 - t177 * t347 + t352) * t249 + (t143 + (t376 - t381) * qJD(3) + t271 * qJD(1)) * t246) * t249) + (-t96 * t246 - t249 * t95 - t35 - t36) * t351 + (t100 * t246 + t249 * t99 + t37 + t38) * t350; m(5) * (t101 * t49 + t102 * t50 + t103 * t41 + t104 * t42) + m(6) * (t23 * t91 + t24 * t92 + t59 * t44 + t60 * t45) + (t246 * t265 + t249 * t264) * t349 + (t266 * t249 + t267 * t246 + (t246 * t264 - t249 * t265) * qJD(1)) * t248 + t316; m(5) * (t41 * t246 - t249 * t42 + (t101 * t246 + t102 * t249) * qJD(1)) + m(6) * t424; m(5) * (t101 * t88 + t102 * t89 + t155 * t41 + t156 * t42 - t29 * t90 + t30 * t65) + m(6) * (-t46 * t10 + t105 * t23 + t106 * t24 + t9 * t43 + t59 * t47 + t60 * t48) + (t4 / 0.2e1 + t3 / 0.2e1 + t338 * t349 + t451 * qJD(1) / 0.2e1) * t249 + (t2 / 0.2e1 + t1 / 0.2e1 + t339 * t349 - t452 * qJD(1) / 0.2e1) * t246 + ((t246 * t338 - t249 * t339) * qJD(1) - (t7 + t8) * t246 / 0.2e1 + (t5 + t6) * t249 / 0.2e1 + t441 * qJD(3) / 0.2e1) * t248 + (t442 * qJD(1) + t439 * t246 + t440 * t249) * t245 / 0.2e1; (t23 * t60 + t24 * t59 - t46 * t9) * t420 + (t101 * t42 + t102 * t41 - t30 * t90) * t421 + (((-t245 * t402 - t451) * t249 + (t245 * t403 + t452) * t246) * qJD(3) + t316) * t245 + ((t1 + t2) * t249 + (-t3 - t4) * t246 + (-t440 * t246 + t439 * t249) * t245 + (t450 * t245 + t442 * t248) * qJD(3) + (-t441 * t245 - t246 * t451 - t249 * t452) * qJD(1)) * t248; m(6) * ((t246 * t91 - t249 * t92) * t349 - t425 * t248); t353 * t342; m(6) * ((t10 + (t105 * t246 - t106 * t249) * qJD(3)) * t245 + (qJD(3) * t43 - t426) * t248); m(6) * ((t9 + (t246 * t60 - t249 * t59) * qJD(3)) * t245 + (-qJD(3) * t46 - t424) * t248); 0.2e1 * (0.1e1 - t353) * t248 * t342;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t25(1), t25(2), t25(4), t25(7), t25(11); t25(2), t25(3), t25(5), t25(8), t25(12); t25(4), t25(5), t25(6), t25(9), t25(13); t25(7), t25(8), t25(9), t25(10), t25(14); t25(11), t25(12), t25(13), t25(14), t25(15);];
Mq = res;
