% Calculate time derivative of joint inertia matrix for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:24
% EndTime: 2019-12-31 17:29:41
% DurationCPUTime: 9.11s
% Computational Cost: add. (37138->798), mult. (103397->1135), div. (0->0), fcn. (116156->10), ass. (0->324)
t306 = cos(pkin(4));
t312 = cos(qJ(2));
t313 = cos(qJ(1));
t363 = t312 * t313;
t309 = sin(qJ(2));
t310 = sin(qJ(1));
t366 = t309 * t310;
t289 = -t306 * t366 + t363;
t344 = t306 * t363;
t320 = t344 - t366;
t249 = qJD(1) * t289 + qJD(2) * t320;
t308 = sin(qJ(3));
t364 = t310 * t312;
t365 = t309 * t313;
t287 = t306 * t365 + t364;
t305 = sin(pkin(4));
t379 = cos(qJ(3));
t339 = t305 * t379;
t317 = -t287 * t308 - t313 * t339;
t354 = qJD(1) * t310;
t338 = t305 * t354;
t187 = qJD(3) * t317 + t249 * t379 + t308 * t338;
t367 = t305 * t313;
t262 = t287 * t379 - t308 * t367;
t307 = sin(qJ(4));
t311 = cos(qJ(4));
t222 = t262 * t311 - t307 * t320;
t288 = t306 * t364 + t365;
t248 = qJD(1) * t288 + qJD(2) * t287;
t135 = -qJD(4) * t222 - t187 * t307 + t248 * t311;
t221 = -t262 * t307 - t311 * t320;
t136 = qJD(4) * t221 + t187 * t311 + t248 * t307;
t152 = Icges(5,5) * t222 + Icges(5,6) * t221 - Icges(5,3) * t317;
t154 = Icges(5,4) * t222 + Icges(5,2) * t221 - Icges(5,6) * t317;
t156 = Icges(5,1) * t222 + Icges(5,4) * t221 - Icges(5,5) * t317;
t332 = qJD(1) * t339;
t186 = qJD(3) * t262 + t249 * t308 - t310 * t332;
t76 = Icges(5,5) * t136 + Icges(5,6) * t135 + Icges(5,3) * t186;
t78 = Icges(5,4) * t136 + Icges(5,2) * t135 + Icges(5,6) * t186;
t80 = Icges(5,1) * t136 + Icges(5,4) * t135 + Icges(5,5) * t186;
t18 = t135 * t154 + t136 * t156 + t152 * t186 + t221 * t78 + t222 * t80 - t317 * t76;
t247 = -qJD(1) * t287 - qJD(2) * t288;
t369 = t305 * t310;
t264 = t289 * t379 + t308 * t369;
t184 = qJD(3) * t264 + t247 * t308 - t313 * t332;
t223 = -t264 * t307 + t288 * t311;
t224 = t264 * t311 + t288 * t307;
t318 = -t289 * t308 + t310 * t339;
t153 = Icges(5,5) * t224 + Icges(5,6) * t223 - Icges(5,3) * t318;
t155 = Icges(5,4) * t224 + Icges(5,2) * t223 - Icges(5,6) * t318;
t157 = Icges(5,1) * t224 + Icges(5,4) * t223 - Icges(5,5) * t318;
t353 = qJD(1) * t313;
t337 = t305 * t353;
t185 = qJD(3) * t318 + t247 * t379 + t308 * t337;
t350 = qJD(2) * t312;
t246 = -qJD(1) * t344 - t313 * t350 + (qJD(2) * t306 + qJD(1)) * t366;
t133 = -qJD(4) * t224 - t185 * t307 - t246 * t311;
t134 = qJD(4) * t223 + t185 * t311 - t246 * t307;
t75 = Icges(5,5) * t134 + Icges(5,6) * t133 + Icges(5,3) * t184;
t77 = Icges(5,4) * t134 + Icges(5,2) * t133 + Icges(5,6) * t184;
t79 = Icges(5,1) * t134 + Icges(5,4) * t133 + Icges(5,5) * t184;
t19 = t135 * t155 + t136 * t157 + t153 * t186 + t221 * t77 + t222 * t79 - t317 * t75;
t285 = t306 * t308 + t309 * t339;
t335 = t305 * t350;
t257 = qJD(3) * t285 + t308 * t335;
t370 = t305 * t309;
t284 = -t306 * t379 + t308 * t370;
t258 = -qJD(3) * t284 + t335 * t379;
t368 = t305 * t312;
t321 = -t285 * t311 + t307 * t368;
t351 = qJD(2) * t309;
t336 = t305 * t351;
t203 = qJD(4) * t321 - t258 * t307 + t311 * t336;
t259 = -t285 * t307 - t311 * t368;
t204 = qJD(4) * t259 + t258 * t311 + t307 * t336;
t126 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t257;
t127 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t257;
t128 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t257;
t188 = -Icges(5,5) * t321 + Icges(5,6) * t259 + Icges(5,3) * t284;
t189 = -Icges(5,4) * t321 + Icges(5,2) * t259 + Icges(5,6) * t284;
t190 = -Icges(5,1) * t321 + Icges(5,4) * t259 + Icges(5,5) * t284;
t33 = -t126 * t317 + t127 * t221 + t128 * t222 + t135 * t189 + t136 * t190 + t186 * t188;
t61 = -t152 * t317 + t154 * t221 + t156 * t222;
t62 = -t153 * t317 + t155 * t221 + t157 * t222;
t87 = -t188 * t317 + t189 * t221 + t190 * t222;
t2 = -t18 * t317 + t184 * t62 + t186 * t61 - t19 * t318 + t257 * t87 + t284 * t33;
t392 = -t2 / 0.2e1;
t304 = t313 * pkin(1);
t355 = pkin(6) * t369 + t304;
t391 = 2 * m(3);
t390 = 2 * m(4);
t389 = 2 * m(5);
t388 = t184 / 0.2e1;
t387 = t186 / 0.2e1;
t386 = t257 / 0.2e1;
t385 = -t317 / 0.2e1;
t384 = -t318 / 0.2e1;
t383 = t284 / 0.2e1;
t382 = t306 / 0.2e1;
t381 = t310 / 0.2e1;
t380 = -rSges(5,3) - pkin(8);
t378 = pkin(3) * t262;
t377 = t187 * pkin(3);
t46 = t284 * t126 + t259 * t127 - t128 * t321 + t257 * t188 + t203 * t189 + t204 * t190;
t94 = t188 * t284 + t189 * t259 - t190 * t321;
t376 = t94 * t257 + t46 * t284;
t209 = Icges(4,5) * t258 - Icges(4,6) * t257 + Icges(4,3) * t336;
t210 = Icges(4,4) * t258 - Icges(4,2) * t257 + Icges(4,6) * t336;
t211 = Icges(4,1) * t258 - Icges(4,4) * t257 + Icges(4,5) * t336;
t225 = Icges(4,5) * t285 - Icges(4,6) * t284 - Icges(4,3) * t368;
t226 = Icges(4,4) * t285 - Icges(4,2) * t284 - Icges(4,6) * t368;
t227 = Icges(4,1) * t285 - Icges(4,4) * t284 - Icges(4,5) * t368;
t74 = -t209 * t368 - t284 * t210 + t285 * t211 + t225 * t336 - t257 * t226 + t258 * t227;
t375 = -t46 - t74;
t81 = t134 * rSges(5,1) + t133 * rSges(5,2) + t184 * rSges(5,3);
t374 = t185 * pkin(3) + pkin(8) * t184 + t81;
t328 = -t136 * rSges(5,1) - t135 * rSges(5,2);
t82 = t186 * rSges(5,3) - t328;
t373 = t186 * pkin(8) + t377 + t82;
t372 = Icges(3,4) * t309;
t371 = Icges(3,4) * t312;
t129 = rSges(5,1) * t204 + rSges(5,2) * t203 + rSges(5,3) * t257;
t362 = pkin(3) * t258 + pkin(8) * t257 + t129;
t327 = -rSges(5,1) * t222 - rSges(5,2) * t221;
t158 = -rSges(5,3) * t317 - t327;
t361 = -pkin(8) * t317 + t158 + t378;
t159 = t224 * rSges(5,1) + t223 * rSges(5,2) - rSges(5,3) * t318;
t360 = t264 * pkin(3) - pkin(8) * t318 + t159;
t191 = -rSges(5,1) * t321 + rSges(5,2) * t259 + rSges(5,3) * t284;
t359 = pkin(3) * t285 + pkin(8) * t284 + t191;
t212 = rSges(4,1) * t258 - rSges(4,2) * t257 + rSges(4,3) * t336;
t352 = qJD(2) * t305;
t278 = (pkin(2) * t312 + pkin(7) * t309) * t352;
t358 = -t212 - t278;
t228 = rSges(4,1) * t285 - rSges(4,2) * t284 - rSges(4,3) * t368;
t290 = (pkin(2) * t309 - pkin(7) * t312) * t305;
t357 = -t228 - t290;
t251 = pkin(2) * t287 - pkin(7) * t320;
t252 = t289 * pkin(2) + pkin(7) * t288;
t356 = t251 * t369 + t252 * t367;
t20 = t152 * t257 + t154 * t203 + t156 * t204 + t259 * t78 + t284 * t76 - t321 * t80;
t349 = t20 / 0.2e1 + t33 / 0.2e1;
t21 = t153 * t257 + t155 * t203 + t157 * t204 + t259 * t77 + t284 * t75 - t321 * t79;
t32 = -t126 * t318 + t127 * t223 + t128 * t224 + t133 * t189 + t134 * t190 + t184 * t188;
t348 = t32 / 0.2e1 + t21 / 0.2e1;
t66 = t152 * t284 + t154 * t259 - t156 * t321;
t347 = t66 / 0.2e1 + t87 / 0.2e1;
t67 = t153 * t284 + t155 * t259 - t157 * t321;
t88 = -t188 * t318 + t189 * t223 + t190 * t224;
t346 = t88 / 0.2e1 + t67 / 0.2e1;
t343 = -t278 - t362;
t120 = t185 * rSges(4,1) - t184 * rSges(4,2) - t246 * rSges(4,3);
t201 = t247 * pkin(2) - pkin(7) * t246;
t202 = t249 * pkin(2) + t248 * pkin(7);
t342 = t201 * t367 + t202 * t369 + t251 * t337;
t341 = -t290 - t359;
t176 = t247 * rSges(3,1) + t246 * rSges(3,2) + rSges(3,3) * t337;
t200 = t264 * rSges(4,1) + rSges(4,2) * t318 + t288 * rSges(4,3);
t236 = t289 * rSges(3,1) - t288 * rSges(3,2) + rSges(3,3) * t369;
t334 = -t310 * pkin(1) + pkin(6) * t367;
t333 = t357 * t313;
t331 = -pkin(1) * t354 + pkin(6) * t337;
t330 = t341 * t313;
t329 = -t249 * rSges(3,1) + t248 * rSges(3,2);
t326 = t252 + t355;
t115 = Icges(4,5) * t187 - Icges(4,6) * t186 + Icges(4,3) * t248;
t117 = Icges(4,4) * t187 - Icges(4,2) * t186 + Icges(4,6) * t248;
t119 = Icges(4,1) * t187 - Icges(4,4) * t186 + Icges(4,5) * t248;
t193 = Icges(4,5) * t262 + Icges(4,6) * t317 - Icges(4,3) * t320;
t195 = Icges(4,4) * t262 + Icges(4,2) * t317 - Icges(4,6) * t320;
t197 = Icges(4,1) * t262 + Icges(4,4) * t317 - Icges(4,5) * t320;
t49 = -t117 * t284 + t119 * t285 - t195 * t257 + t197 * t258 + (-t115 * t312 + t193 * t351) * t305;
t60 = -t186 * t226 + t187 * t227 - t209 * t320 + t210 * t317 + t211 * t262 + t225 * t248;
t325 = t60 / 0.2e1 + t49 / 0.2e1 + t349;
t114 = Icges(4,5) * t185 - Icges(4,6) * t184 - Icges(4,3) * t246;
t116 = Icges(4,4) * t185 - Icges(4,2) * t184 - Icges(4,6) * t246;
t118 = Icges(4,1) * t185 - Icges(4,4) * t184 - Icges(4,5) * t246;
t194 = Icges(4,5) * t264 + Icges(4,6) * t318 + Icges(4,3) * t288;
t196 = Icges(4,4) * t264 + Icges(4,2) * t318 + Icges(4,6) * t288;
t198 = Icges(4,1) * t264 + Icges(4,4) * t318 + Icges(4,5) * t288;
t50 = -t116 * t284 + t118 * t285 - t196 * t257 + t198 * t258 + (-t114 * t312 + t194 * t351) * t305;
t59 = -t184 * t226 + t185 * t227 + t209 * t288 + t210 * t318 + t211 * t264 - t225 * t246;
t324 = t50 / 0.2e1 + t59 / 0.2e1 + t348;
t101 = -t193 * t368 - t195 * t284 + t197 * t285;
t112 = -t225 * t320 + t226 * t317 + t227 * t262;
t323 = t112 / 0.2e1 + t101 / 0.2e1 + t347;
t102 = -t194 * t368 - t196 * t284 + t198 * t285;
t113 = t225 * t288 + t226 * t318 + t227 * t264;
t322 = -t102 / 0.2e1 - t113 / 0.2e1 - t346;
t319 = -t251 + t334;
t121 = t187 * rSges(4,1) - t186 * rSges(4,2) + t248 * rSges(4,3);
t199 = rSges(4,1) * t262 + rSges(4,2) * t317 - rSges(4,3) * t320;
t316 = t201 + t331;
t235 = t287 * rSges(3,1) + rSges(3,2) * t320 - rSges(3,3) * t367;
t271 = Icges(3,6) * t306 + (Icges(3,2) * t312 + t372) * t305;
t272 = Icges(3,5) * t306 + (Icges(3,1) * t309 + t371) * t305;
t274 = (Icges(3,5) * t312 - Icges(3,6) * t309) * t352;
t275 = (-Icges(3,2) * t309 + t371) * t352;
t276 = (Icges(3,1) * t312 - t372) * t352;
t315 = -t271 * t336 + t272 * t335 + t306 * t274 + t275 * t368 + t276 * t370;
t314 = -qJD(1) * t355 - t202;
t277 = (rSges(3,1) * t312 - rSges(3,2) * t309) * t352;
t273 = rSges(3,3) * t306 + (rSges(3,1) * t309 + rSges(3,2) * t312) * t305;
t270 = Icges(3,3) * t306 + (Icges(3,5) * t309 + Icges(3,6) * t312) * t305;
t268 = t290 * t338;
t240 = t306 * t252;
t234 = Icges(3,1) * t289 - Icges(3,4) * t288 + Icges(3,5) * t369;
t233 = Icges(3,1) * t287 + Icges(3,4) * t320 - Icges(3,5) * t367;
t232 = Icges(3,4) * t289 - Icges(3,2) * t288 + Icges(3,6) * t369;
t231 = Icges(3,4) * t287 + Icges(3,2) * t320 - Icges(3,6) * t367;
t230 = Icges(3,5) * t289 - Icges(3,6) * t288 + Icges(3,3) * t369;
t229 = Icges(3,5) * t287 + Icges(3,6) * t320 - Icges(3,3) * t367;
t217 = t236 + t355;
t216 = -t235 + t334;
t208 = -t306 * t235 - t273 * t367;
t207 = t236 * t306 - t273 * t369;
t192 = t306 * t201;
t177 = rSges(3,3) * t338 - t329;
t175 = Icges(3,1) * t249 - Icges(3,4) * t248 + Icges(3,5) * t338;
t174 = Icges(3,1) * t247 + Icges(3,4) * t246 + Icges(3,5) * t337;
t173 = Icges(3,4) * t249 - Icges(3,2) * t248 + Icges(3,6) * t338;
t172 = Icges(3,4) * t247 + Icges(3,2) * t246 + Icges(3,6) * t337;
t171 = Icges(3,5) * t249 - Icges(3,6) * t248 + Icges(3,3) * t338;
t170 = Icges(3,5) * t247 + Icges(3,6) * t246 + Icges(3,3) * t337;
t167 = (-t304 + (-rSges(3,3) - pkin(6)) * t369) * qJD(1) + t329;
t166 = t331 + t176;
t165 = t270 * t369 - t271 * t288 + t272 * t289;
t164 = -t270 * t367 + t271 * t320 + t287 * t272;
t162 = t326 + t200;
t161 = -t199 + t319;
t160 = t315 * t306;
t151 = t306 * t176 + (-t273 * t353 - t277 * t310) * t305;
t150 = -t306 * t177 + (t273 * t354 - t277 * t313) * t305;
t149 = -t200 * t368 - t228 * t288;
t148 = t199 * t368 - t228 * t320;
t147 = t230 * t306 + (t232 * t312 + t234 * t309) * t305;
t146 = t229 * t306 + (t231 * t312 + t233 * t309) * t305;
t141 = t230 * t369 - t232 * t288 + t234 * t289;
t140 = t229 * t369 - t231 * t288 + t233 * t289;
t139 = -t230 * t367 + t232 * t320 + t287 * t234;
t138 = -t229 * t367 + t231 * t320 + t287 * t233;
t137 = -t225 * t368 - t226 * t284 + t227 * t285;
t130 = t137 * t336;
t125 = t199 * t288 + t200 * t320;
t124 = (-t199 - t251) * t306 + t305 * t333;
t123 = t200 * t306 + t357 * t369 + t240;
t109 = -t248 * t271 + t249 * t272 + t320 * t275 + t287 * t276 + (t270 * t354 - t274 * t313) * t305;
t108 = t246 * t271 + t247 * t272 - t288 * t275 + t289 * t276 + (t270 * t353 + t274 * t310) * t305;
t107 = t326 + t360;
t106 = -t317 * t380 + t319 + t327 - t378;
t105 = (t199 * t310 + t200 * t313) * t305 + t356;
t104 = t159 * t284 + t191 * t318;
t103 = -t158 * t284 - t191 * t317;
t100 = -t121 + t314;
t99 = t316 + t120;
t98 = t194 * t288 + t196 * t318 + t198 * t264;
t97 = t193 * t288 + t195 * t318 + t197 * t264;
t96 = -t194 * t320 + t196 * t317 + t198 * t262;
t95 = -t193 * t320 + t195 * t317 + t197 * t262;
t93 = t94 * t336;
t92 = -t158 * t318 + t159 * t317;
t90 = -t288 * t359 - t360 * t368;
t89 = -t320 * t359 + t361 * t368;
t86 = t170 * t306 + (t172 * t312 + t174 * t309 + (-t232 * t309 + t234 * t312) * qJD(2)) * t305;
t85 = t171 * t306 + (t173 * t312 + t175 * t309 + (-t231 * t309 + t233 * t312) * qJD(2)) * t305;
t84 = (-t251 - t361) * t306 + t305 * t330;
t83 = t306 * t360 + t341 * t369 + t240;
t73 = t74 * t306;
t72 = t288 * t361 + t320 * t360;
t71 = t306 * t120 + t192 + (qJD(1) * t333 + t310 * t358) * t305;
t70 = t268 + (-t121 - t202) * t306 + (t228 * t354 + t313 * t358) * t305;
t69 = -t212 * t320 + t228 * t248 + (t121 * t312 - t199 * t351) * t305;
t68 = -t212 * t288 + t228 * t246 + (-t120 * t312 + t200 * t351) * t305;
t65 = (t310 * t361 + t313 * t360) * t305 + t356;
t64 = -t153 * t318 + t155 * t223 + t157 * t224;
t63 = -t152 * t318 + t154 * t223 + t156 * t224;
t58 = t120 * t320 + t121 * t288 - t199 * t246 - t200 * t248;
t57 = t186 * t380 + t314 + t328 - t377;
t56 = t316 + t374;
t55 = (t120 * t313 + t121 * t310 + (t199 * t313 + (-t200 - t252) * t310) * qJD(1)) * t305 + t342;
t54 = t113 * t306 + (t310 * t98 - t313 * t97) * t305;
t53 = t112 * t306 + (t310 * t96 - t313 * t95) * t305;
t52 = -t113 * t368 + t288 * t98 - t320 * t97;
t51 = -t112 * t368 + t288 * t96 - t320 * t95;
t48 = -t129 * t317 - t158 * t257 + t186 * t191 - t284 * t82;
t47 = t129 * t318 + t159 * t257 - t184 * t191 + t284 * t81;
t45 = t46 * t306;
t43 = -t114 * t320 + t116 * t317 + t118 * t262 - t186 * t196 + t187 * t198 + t194 * t248;
t42 = -t115 * t320 + t117 * t317 + t119 * t262 - t186 * t195 + t187 * t197 + t193 * t248;
t41 = t114 * t288 + t116 * t318 + t118 * t264 - t184 * t196 + t185 * t198 - t194 * t246;
t40 = t115 * t288 + t117 * t318 + t119 * t264 - t184 * t195 + t185 * t197 - t193 * t246;
t39 = t192 + t374 * t306 + (qJD(1) * t330 + t310 * t343) * t305;
t38 = t268 + (-t202 - t373) * t306 + (t313 * t343 + t354 * t359) * t305;
t37 = -t362 * t320 + t359 * t248 + (t312 * t373 - t351 * t361) * t305;
t36 = -t362 * t288 + t359 * t246 + (-t312 * t374 + t351 * t360) * t305;
t35 = t94 * t306 + (t310 * t67 - t313 * t66) * t305;
t34 = t288 * t67 - t320 * t66 - t368 * t94;
t31 = t284 * t94 - t317 * t66 - t318 * t67;
t30 = t158 * t184 - t159 * t186 + t317 * t81 - t318 * t82;
t29 = t88 * t306 + (t310 * t64 - t313 * t63) * t305;
t28 = t87 * t306 + (t310 * t62 - t313 * t61) * t305;
t27 = t288 * t64 - t320 * t63 - t368 * t88;
t26 = t288 * t62 - t320 * t61 - t368 * t87;
t25 = t284 * t88 - t317 * t63 - t318 * t64;
t24 = t284 * t87 - t317 * t61 - t318 * t62;
t23 = -t246 * t361 - t248 * t360 + t288 * t373 + t320 * t374;
t22 = (t374 * t313 + t373 * t310 + (t361 * t313 + (-t252 - t360) * t310) * qJD(1)) * t305 + t342;
t17 = t133 * t155 + t134 * t157 + t153 * t184 + t223 * t77 + t224 * t79 - t318 * t75;
t16 = t133 * t154 + t134 * t156 + t152 * t184 + t223 * t78 + t224 * t80 - t318 * t76;
t15 = t73 + (t50 * t310 - t49 * t313 + (t101 * t310 + t102 * t313) * qJD(1)) * t305;
t14 = t101 * t248 - t102 * t246 + t50 * t288 - t320 * t49 - t368 * t74 + t130;
t13 = t60 * t306 + (t310 * t43 - t313 * t42 + (t310 * t95 + t313 * t96) * qJD(1)) * t305;
t12 = t59 * t306 + (t310 * t41 - t313 * t40 + (t310 * t97 + t313 * t98) * qJD(1)) * t305;
t11 = -t246 * t96 + t248 * t95 - t320 * t42 + t288 * t43 + (t112 * t351 - t312 * t60) * t305;
t10 = -t246 * t98 + t248 * t97 - t320 * t40 + t288 * t41 + (t113 * t351 - t312 * t59) * t305;
t9 = t45 + (-t20 * t313 + t21 * t310 + (t66 * t310 + t67 * t313) * qJD(1)) * t305;
t8 = -t20 * t320 + t21 * t288 - t67 * t246 + t66 * t248 - t368 * t46 + t93;
t7 = t67 * t184 + t66 * t186 - t20 * t317 - t21 * t318 + t376;
t6 = t33 * t306 + (-t18 * t313 + t19 * t310 + (t310 * t61 + t313 * t62) * qJD(1)) * t305;
t5 = t32 * t306 + (-t16 * t313 + t17 * t310 + (t310 * t63 + t313 * t64) * qJD(1)) * t305;
t4 = -t18 * t320 + t19 * t288 - t246 * t62 + t248 * t61 + (-t312 * t33 + t351 * t87) * t305;
t3 = -t16 * t320 + t17 * t288 - t246 * t64 + t248 * t63 + (-t312 * t32 + t351 * t88) * t305;
t1 = -t16 * t317 - t17 * t318 + t184 * t64 + t186 * t63 + t257 * t88 + t284 * t32;
t44 = [(t106 * t57 + t107 * t56) * t389 + (t100 * t161 + t162 * t99) * t390 + (t166 * t217 + t167 * t216) * t391 + t315 - t375; t73 + t160 + t45 + m(4) * (t100 * t124 + t123 * t99 + t161 * t70 + t162 * t71) + m(5) * (t106 * t38 + t107 * t39 + t56 * t83 + t57 * t84) + m(3) * (t150 * t216 + t151 * t217 + t166 * t207 + t167 * t208) + ((-t85 / 0.2e1 - t109 / 0.2e1 - t325) * t313 + (t86 / 0.2e1 + t108 / 0.2e1 + t324) * t310 + ((t147 / 0.2e1 + t165 / 0.2e1 - t322) * t313 + (t146 / 0.2e1 + t164 / 0.2e1 + t323) * t310) * qJD(1)) * t305; (t22 * t65 + t38 * t84 + t39 * t83) * t389 + (t105 * t55 + t123 * t71 + t124 * t70) * t390 + (t208 * t150 + t207 * t151) * t391 + (t5 + t12) * t369 + (-t6 - t13) * t367 + (t28 + t53) * t338 + (t29 + t54) * t337 + ((-t140 * t313 + t141 * t310) * t337 + ((-t288 * t172 + t289 * t174 + t246 * t232 + t247 * t234) * t310 + t141 * t353 - (-t288 * t173 + t289 * t175 + t246 * t231 + t247 * t233) * t313 + t140 * t354) * t369 + (-t138 * t313 + t139 * t310) * t338 - ((t172 * t320 + t287 * t174 - t248 * t232 + t249 * t234) * t310 + t139 * t353 - (t173 * t320 + t287 * t175 - t248 * t231 + t249 * t233) * t313 + t138 * t354) * t367 + (((t170 * t310 + t230 * t353) * t310 - (t171 * t310 + t229 * t353) * t313) * t369 - ((-t170 * t313 + t230 * t354) * t310 - (-t171 * t313 + t229 * t354) * t313) * t367 + (t235 * t310 + t236 * t313) * (t176 * t313 + t177 * t310 + (t235 * t313 - t236 * t310) * qJD(1)) * t391) * t305) * t305 + (t165 * t337 + t108 * t369 + t164 * t338 - t109 * t367 + t9 + t15 + t160 + (t310 * t86 - t313 * t85 + (t146 * t310 + t147 * t313) * qJD(1)) * t305) * t306; t130 + t93 + t375 * t368 + m(5) * (t106 * t37 + t107 * t36 + t56 * t90 + t57 * t89) + m(4) * (t100 * t148 + t149 * t99 + t161 * t69 + t162 * t68) + t324 * t288 - t325 * t320 + t323 * t248 + t322 * t246; (t8 / 0.2e1 + t14 / 0.2e1) * t306 + (t5 / 0.2e1 + t12 / 0.2e1) * t288 - (t6 / 0.2e1 + t13 / 0.2e1) * t320 + (t28 / 0.2e1 + t53 / 0.2e1) * t248 + (-t29 / 0.2e1 - t54 / 0.2e1) * t246 + m(5) * (t22 * t72 + t23 * t65 + t36 * t83 + t37 * t84 + t38 * t89 + t39 * t90) + m(4) * (t105 * t58 + t123 * t68 + t124 * t69 + t125 * t55 + t148 * t70 + t149 * t71) + ((-t4 / 0.2e1 - t11 / 0.2e1) * t313 + (-t9 / 0.2e1 - t15 / 0.2e1) * t312 + (t3 / 0.2e1 + t10 / 0.2e1) * t310 + (t35 / 0.2e1 + t137 * t382 + (-t101 * t313 + t102 * t310) * t305 / 0.2e1) * t351 + ((t27 / 0.2e1 + t52 / 0.2e1) * t313 + (t26 / 0.2e1 + t51 / 0.2e1) * t310) * qJD(1)) * t305; (-t14 - t8) * t368 + (t10 + t3) * t288 - (t11 + t4) * t320 + (t26 + t51) * t248 + (-t27 - t52) * t246 + (-t101 * t320 + t102 * t288 - t137 * t368 + t34) * t336 + (t23 * t72 + t36 * t90 + t37 * t89) * t389 + (t125 * t58 + t148 * t69 + t149 * t68) * t390; m(5) * (t103 * t57 + t104 * t56 + t106 * t48 + t107 * t47) - t348 * t318 - t349 * t317 + t347 * t186 + t346 * t184 + t376; t29 * t388 + t5 * t384 + t28 * t387 + t6 * t385 + t7 * t382 + t35 * t386 + t9 * t383 + m(5) * (t103 * t38 + t104 * t39 + t22 * t92 + t30 * t65 + t47 * t83 + t48 * t84) + (t1 * t381 + t313 * t392 + (t313 * t25 / 0.2e1 + t24 * t381) * qJD(1)) * t305; t248 * t24 / 0.2e1 + t320 * t392 + t27 * t388 + t3 * t384 + t26 * t387 + t4 * t385 + t34 * t386 + t8 * t383 - t246 * t25 / 0.2e1 + t288 * t1 / 0.2e1 + m(5) * (t103 * t37 + t104 * t36 + t23 * t92 + t30 * t72 + t47 * t90 + t48 * t89) + (t31 * t351 / 0.2e1 - t312 * t7 / 0.2e1) * t305; (t103 * t48 + t104 * t47 + t30 * t92) * t389 + t184 * t25 - t318 * t1 + t186 * t24 - t317 * t2 + t257 * t31 + t284 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t44(1), t44(2), t44(4), t44(7); t44(2), t44(3), t44(5), t44(8); t44(4), t44(5), t44(6), t44(9); t44(7), t44(8), t44(9), t44(10);];
Mq = res;
