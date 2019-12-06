% Calculate time derivative of joint inertia matrix for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:00
% EndTime: 2019-12-05 16:12:27
% DurationCPUTime: 11.64s
% Computational Cost: add. (10323->593), mult. (29999->829), div. (0->0), fcn. (31126->8), ass. (0->265)
t396 = Icges(5,1) + Icges(6,1);
t390 = -Icges(5,4) + Icges(6,5);
t389 = Icges(6,4) + Icges(5,5);
t395 = Icges(5,2) + Icges(6,3);
t394 = Icges(6,2) + Icges(5,3);
t388 = -Icges(5,6) + Icges(6,6);
t242 = sin(pkin(7));
t246 = cos(qJ(3));
t247 = cos(qJ(2));
t320 = t246 * t247;
t243 = cos(pkin(7));
t244 = sin(qJ(3));
t324 = t243 * t244;
t221 = t242 * t320 - t324;
t241 = sin(pkin(8));
t245 = sin(qJ(2));
t335 = cos(pkin(8));
t289 = t245 * t335;
t190 = t221 * t241 - t242 * t289;
t326 = t242 * t245;
t191 = t221 * t335 + t241 * t326;
t321 = t244 * t247;
t220 = t242 * t321 + t243 * t246;
t108 = Icges(6,5) * t191 + Icges(6,6) * t220 + Icges(6,3) * t190;
t114 = Icges(5,4) * t191 - Icges(5,2) * t190 + Icges(5,6) * t220;
t384 = t108 - t114;
t225 = t242 * t244 + t243 * t320;
t192 = t225 * t241 - t243 * t289;
t323 = t243 * t245;
t193 = t225 * t335 + t241 * t323;
t224 = -t242 * t246 + t243 * t321;
t109 = Icges(6,5) * t193 + Icges(6,6) * t224 + Icges(6,3) * t192;
t115 = Icges(5,4) * t193 - Icges(5,2) * t192 + Icges(5,6) * t224;
t383 = t109 - t115;
t110 = Icges(5,5) * t191 - Icges(5,6) * t190 + Icges(5,3) * t220;
t112 = Icges(6,4) * t191 + Icges(6,2) * t220 + Icges(6,6) * t190;
t382 = t110 + t112;
t111 = Icges(5,5) * t193 - Icges(5,6) * t192 + Icges(5,3) * t224;
t113 = Icges(6,4) * t193 + Icges(6,2) * t224 + Icges(6,6) * t192;
t381 = t111 + t113;
t116 = Icges(6,1) * t191 + Icges(6,4) * t220 + Icges(6,5) * t190;
t118 = Icges(5,1) * t191 - Icges(5,4) * t190 + Icges(5,5) * t220;
t380 = t116 + t118;
t117 = Icges(6,1) * t193 + Icges(6,4) * t224 + Icges(6,5) * t192;
t119 = Icges(5,1) * t193 - Icges(5,4) * t192 + Icges(5,5) * t224;
t379 = t117 + t119;
t331 = Icges(4,4) * t246;
t264 = -Icges(4,2) * t244 + t331;
t209 = -Icges(4,6) * t247 + t245 * t264;
t288 = t247 * t335;
t327 = t241 * t246;
t222 = t245 * t327 + t288;
t223 = -t247 * t241 + t246 * t289;
t322 = t244 * t245;
t357 = t388 * t222 + t389 * t223 + t394 * t322;
t393 = -t209 + t357;
t392 = t384 * t190 + t380 * t191 + t382 * t220;
t391 = t383 * t192 + t379 * t193 + t381 * t224;
t370 = t395 * t222 + t390 * t223 + t388 * t322;
t369 = t390 * t222 + t396 * t223 + t389 * t322;
t159 = Icges(4,5) * t225 - Icges(4,6) * t224 + Icges(4,3) * t323;
t161 = Icges(4,4) * t225 - Icges(4,2) * t224 + Icges(4,6) * t323;
t163 = Icges(4,1) * t225 - Icges(4,4) * t224 + Icges(4,5) * t323;
t387 = t159 * t326 + t163 * t221 + t379 * t191 + t383 * t190 + (-t161 + t381) * t220;
t158 = Icges(4,5) * t221 - Icges(4,6) * t220 + Icges(4,3) * t326;
t160 = Icges(4,4) * t221 - Icges(4,2) * t220 + Icges(4,6) * t326;
t162 = Icges(4,1) * t221 - Icges(4,4) * t220 + Icges(4,5) * t326;
t386 = t158 * t323 + t162 * t225 + t380 * t193 + t384 * t192 + (-t160 + t382) * t224;
t287 = qJD(2) * t335;
t303 = qJD(3) * t245;
t305 = qJD(2) * t247;
t198 = -t241 * t244 * t303 - t245 * t287 + t305 * t327;
t304 = qJD(3) * t244;
t199 = -t289 * t304 + (t241 * t245 + t246 * t288) * qJD(2);
t302 = qJD(3) * t246;
t291 = t245 * t302;
t248 = t244 * t305 + t291;
t378 = t395 * t198 + t390 * t199 + t388 * t248;
t377 = t390 * t198 + t396 * t199 + t389 * t248;
t263 = Icges(4,5) * t246 - Icges(4,6) * t244;
t208 = -Icges(4,3) * t247 + t245 * t263;
t332 = Icges(4,4) * t244;
t266 = Icges(4,1) * t246 - t332;
t210 = -Icges(4,5) * t247 + t245 * t266;
t376 = t370 * t190 + t369 * t191 + t208 * t326 + t210 * t221 + t220 * t393;
t375 = -t370 * t192 - t369 * t193 - t208 * t323 - t210 * t225 - t224 * t393;
t374 = -(-Icges(4,2) * t246 - t332) * t303 - (Icges(4,6) * t245 + t247 * t264) * qJD(2) + t394 * t248 + t389 * t199 + t388 * t198;
t372 = t386 * t242 + t391 * t243;
t371 = t392 * t242 + t387 * t243;
t262 = -t160 * t244 + t162 * t246;
t368 = -t158 * t247 + t384 * t222 + t380 * t223 + t245 * t262 + t382 * t322;
t261 = -t161 * t244 + t163 * t246;
t367 = -t159 * t247 + t383 * t222 + t379 * t223 + t245 * t261 + t381 * t322;
t366 = rSges(6,1) + pkin(4);
t365 = rSges(6,3) + qJ(5);
t364 = t243 * t247;
t346 = t243 ^ 2;
t347 = t242 ^ 2;
t356 = t346 + t347;
t306 = qJD(2) * t245;
t292 = t246 * t306;
t195 = -qJD(3) * t220 - t242 * t292;
t278 = t247 * t287;
t154 = t195 * t241 - t242 * t278;
t294 = t242 * t305;
t155 = t195 * t335 + t241 * t294;
t290 = t247 * t302;
t194 = -t242 * t290 + (qJD(3) * t243 + t242 * t306) * t244;
t82 = Icges(6,5) * t155 - Icges(6,6) * t194 + Icges(6,3) * t154;
t86 = Icges(6,4) * t155 - Icges(6,2) * t194 + Icges(6,6) * t154;
t90 = Icges(6,1) * t155 - Icges(6,4) * t194 + Icges(6,5) * t154;
t16 = t108 * t154 - t112 * t194 + t116 * t155 + t190 * t82 + t191 * t90 + t220 * t86;
t197 = -qJD(3) * t224 - t243 * t292;
t156 = t197 * t241 - t243 * t278;
t293 = t243 * t305;
t157 = t197 * t335 + t241 * t293;
t196 = -t242 * t304 - t243 * t290 + t306 * t324;
t83 = Icges(6,5) * t157 - Icges(6,6) * t196 + Icges(6,3) * t156;
t87 = Icges(6,4) * t157 - Icges(6,2) * t196 + Icges(6,6) * t156;
t91 = Icges(6,1) * t157 - Icges(6,4) * t196 + Icges(6,5) * t156;
t17 = t109 * t154 - t113 * t194 + t117 * t155 + t190 * t83 + t191 * t91 + t220 * t87;
t179 = (-Icges(4,1) * t244 - t331) * t303 + (Icges(4,5) * t245 + t247 * t266) * qJD(2);
t84 = Icges(5,5) * t155 - Icges(5,6) * t154 - Icges(5,3) * t194;
t88 = Icges(5,4) * t155 - Icges(5,2) * t154 - Icges(5,6) * t194;
t92 = Icges(5,1) * t155 - Icges(5,4) * t154 - Icges(5,5) * t194;
t18 = -t110 * t194 - t114 * t154 + t118 * t155 - t190 * t88 + t191 * t92 + t220 * t84;
t85 = Icges(5,5) * t157 - Icges(5,6) * t156 - Icges(5,3) * t196;
t89 = Icges(5,4) * t157 - Icges(5,2) * t156 - Icges(5,6) * t196;
t93 = Icges(5,1) * t157 - Icges(5,4) * t156 - Icges(5,5) * t196;
t19 = -t111 * t194 - t115 * t154 + t119 * t155 - t190 * t89 + t191 * t93 + t220 * t85;
t329 = t208 * t247;
t330 = ((-Icges(4,5) * t244 - Icges(4,6) * t246) * t303 + (Icges(4,3) * t245 + t247 * t263) * qJD(2)) * t247;
t142 = Icges(4,4) * t195 + Icges(4,2) * t194 + Icges(4,6) * t294;
t144 = Icges(4,1) * t195 + Icges(4,4) * t194 + Icges(4,5) * t294;
t140 = Icges(4,5) * t195 + Icges(4,6) * t194 + Icges(4,3) * t294;
t254 = t140 * t245 + t158 * t305;
t35 = -t220 * t142 + t221 * t144 + t194 * t160 + t195 * t162 + t242 * t254;
t143 = Icges(4,4) * t197 + Icges(4,2) * t196 + Icges(4,6) * t293;
t145 = Icges(4,1) * t197 + Icges(4,4) * t196 + Icges(4,5) * t293;
t141 = Icges(4,5) * t197 + Icges(4,6) * t196 + Icges(4,3) * t293;
t253 = t141 * t245 + t159 * t305;
t36 = -t220 * t143 + t221 * t145 + t194 * t161 + t195 * t163 + t242 * t253;
t69 = t158 * t326 - t160 * t220 + t162 * t221;
t363 = (-t370 * t154 - t369 * t155 - t221 * t179 - t378 * t190 - t377 * t191 + t194 * t393 - t195 * t210 - t374 * t220) * t247 + ((t17 + t19 + t36) * t243 + (t16 + t18 + t35 - t330) * t242) * t245 + (((t69 - t329) * t242 + t371) * t247 + t376 * t245) * qJD(2);
t20 = t108 * t156 - t112 * t196 + t116 * t157 + t192 * t82 + t193 * t90 + t224 * t86;
t21 = t109 * t156 - t113 * t196 + t117 * t157 + t192 * t83 + t193 * t91 + t224 * t87;
t22 = -t110 * t196 - t114 * t156 + t118 * t157 - t192 * t88 + t193 * t92 + t224 * t84;
t23 = -t111 * t196 - t115 * t156 + t119 * t157 - t192 * t89 + t193 * t93 + t224 * t85;
t37 = -t224 * t142 + t225 * t144 + t196 * t160 + t197 * t162 + t243 * t254;
t38 = -t224 * t143 + t225 * t145 + t196 * t161 + t197 * t163 + t243 * t253;
t72 = t159 * t323 - t161 * t224 + t163 * t225;
t362 = (-t370 * t156 - t369 * t157 - t225 * t179 - t378 * t192 - t377 * t193 + t196 * t393 - t197 * t210 - t374 * t224) * t247 + ((t21 + t23 + t38 - t330) * t243 + (t20 + t22 + t37) * t242) * t245 + (((t72 - t329) * t243 + t372) * t247 - t375 * t245) * qJD(2);
t361 = -(qJD(2) * t262 - t140) * t247 - (qJD(2) * t158 - t142 * t244 + t144 * t246 + (-t160 * t246 - t162 * t244) * qJD(3)) * t245 + (-t84 - t86) * t322 - t382 * t248 + (-t90 - t92) * t223 + (-t82 + t88) * t222 - t380 * t199 - t384 * t198;
t360 = (qJD(2) * t261 - t141) * t247 + (qJD(2) * t159 - t143 * t244 + t145 * t246 + (-t161 * t246 - t163 * t244) * qJD(3)) * t245 + (t85 + t87) * t322 + t381 * t248 + (t91 + t93) * t223 + (t83 - t89) * t222 + t379 * t199 + t383 * t198;
t355 = qJD(2) * (t245 * rSges(3,1) + rSges(3,2) * t247);
t258 = t209 * t244 - t210 * t246;
t354 = -t370 * t222 - t369 * t223 + t245 * t258 - t357 * t322 + t329;
t353 = t368 * t242 + t367 * t243;
t350 = 2 * m(4);
t349 = 2 * m(5);
t348 = 2 * m(6);
t340 = -rSges(6,2) * t194 + qJD(5) * t190 + t365 * t154 + t366 * t155;
t339 = -rSges(6,2) * t196 + qJD(5) * t192 + t365 * t156 + t366 * t157;
t125 = pkin(3) * t197 - qJ(4) * t196 + qJD(4) * t224;
t97 = rSges(5,1) * t157 - rSges(5,2) * t156 - rSges(5,3) * t196;
t336 = -t125 - t97;
t325 = t242 * t247;
t124 = pkin(3) * t195 - qJ(4) * t194 + qJD(4) * t220;
t187 = pkin(3) * t221 + qJ(4) * t220;
t319 = t124 * t323 + t187 * t293;
t318 = rSges(6,2) * t220 + t365 * t190 + t366 * t191;
t317 = rSges(6,2) * t224 + t365 * t192 + t366 * t193;
t123 = rSges(5,1) * t193 - rSges(5,2) * t192 + rSges(5,3) * t224;
t189 = pkin(3) * t225 + qJ(4) * t224;
t316 = -t123 - t189;
t315 = rSges(6,2) * t248 + qJD(5) * t222 + t365 * t198 + t366 * t199;
t134 = t199 * rSges(5,1) - t198 * rSges(5,2) + rSges(5,3) * t248;
t176 = (pkin(3) * t305 + qJ(4) * t303) * t246 + (qJ(4) * t305 + (-pkin(3) * qJD(3) + qJD(4)) * t245) * t244;
t314 = -t134 - t176;
t313 = rSges(6,2) * t322 + t365 * t222 + t366 * t223;
t174 = rSges(5,1) * t223 - rSges(5,2) * t222 + rSges(5,3) * t322;
t226 = (pkin(3) * t246 + qJ(4) * t244) * t245;
t312 = -t174 - t226;
t275 = rSges(4,1) * t246 - rSges(4,2) * t244;
t181 = (-rSges(4,1) * t244 - rSges(4,2) * t246) * t303 + (rSges(4,3) * t245 + t247 * t275) * qJD(2);
t277 = pkin(2) * t247 + pkin(6) * t245;
t228 = t277 * qJD(2);
t311 = -t181 - t228;
t310 = t247 * t187 + t226 * t326;
t237 = t245 * pkin(2) - pkin(6) * t247;
t309 = t356 * qJD(2) * t237;
t211 = -rSges(4,3) * t247 + t245 * t275;
t308 = -t211 - t237;
t307 = t356 * t277;
t301 = -t125 - t339;
t300 = t247 * t124 + t176 * t326 + t226 * t294;
t299 = -t189 - t317;
t298 = -t176 - t315;
t297 = -t228 + t314;
t296 = -t226 - t313;
t295 = -t237 + t312;
t286 = t242 * t313;
t285 = t318 * t247;
t284 = t316 * t247;
t283 = t242 * t124 + t243 * t125 - t309;
t282 = -t228 + t298;
t281 = -t237 + t296;
t280 = t242 * t187 + t243 * t189 + t307;
t279 = t299 * t247;
t171 = rSges(4,1) * t221 - rSges(4,2) * t220 + rSges(4,3) * t326;
t172 = rSges(4,1) * t225 - rSges(4,2) * t224 + rSges(4,3) * t323;
t260 = t171 * t243 - t172 * t242;
t249 = qJD(2) * (-Icges(3,5) * t245 - Icges(3,6) * t247);
t215 = t243 * t249;
t214 = t242 * t249;
t186 = t308 * t243;
t185 = t308 * t242;
t180 = t189 * t306;
t175 = t187 * t323;
t153 = t356 * t355;
t151 = t311 * t243;
t150 = t311 * t242;
t147 = t197 * rSges(4,1) + t196 * rSges(4,2) + rSges(4,3) * t293;
t146 = t195 * rSges(4,1) + t194 * rSges(4,2) + rSges(4,3) * t294;
t139 = -t172 * t247 - t211 * t323;
t138 = t171 * t247 + t211 * t326;
t136 = t295 * t243;
t135 = t295 * t242;
t121 = rSges(5,1) * t191 - rSges(5,2) * t190 + rSges(5,3) * t220;
t101 = t260 * t245;
t99 = t281 * t243;
t98 = t281 * t242;
t95 = rSges(5,1) * t155 - rSges(5,2) * t154 - rSges(5,3) * t194;
t79 = t171 * t242 + t172 * t243 + t307;
t78 = t297 * t243;
t77 = t297 * t242;
t68 = t146 * t242 + t147 * t243 - t309;
t67 = t312 * t323 + t284;
t66 = t121 * t247 + t174 * t326 + t310;
t65 = -t181 * t323 - t147 * t247 + (t172 * t245 - t211 * t364) * qJD(2);
t64 = t181 * t326 + t146 * t247 + (-t171 * t245 + t211 * t325) * qJD(2);
t59 = t282 * t243;
t58 = t282 * t242;
t57 = (t146 * t243 - t147 * t242) * t245 + t260 * t305;
t56 = t175 + (t121 * t243 + t242 * t316) * t245;
t51 = t121 * t242 + t123 * t243 + t280;
t50 = t296 * t323 + t279;
t49 = t245 * t286 + t285 + t310;
t40 = t175 + (t242 * t299 + t243 * t318) * t245;
t39 = t242 * t318 + t243 * t317 + t280;
t34 = t242 * t95 + t243 * t97 + t283;
t31 = t123 * t306 + t180 + t336 * t247 + (t245 * t314 + t305 * t312) * t243;
t30 = t134 * t326 + t247 * t95 + (t174 * t325 + (-t121 - t187) * t245) * qJD(2) + t300;
t29 = (t121 * t305 + t245 * t95) * t243 + (qJD(2) * t284 + t245 * t336) * t242 + t319;
t28 = t242 * t340 + t243 * t339 + t283;
t15 = t180 + t317 * t306 + t301 * t247 + (t245 * t298 + t296 * t305) * t243;
t14 = t340 * t247 + t315 * t326 + (t247 * t286 + (-t187 - t318) * t245) * qJD(2) + t300;
t13 = t242 * t38 - t243 * t37;
t12 = t242 * t36 - t243 * t35;
t11 = (qJD(2) * t285 + t245 * t340) * t243 + (qJD(2) * t279 + t245 * t301) * t242 + t319;
t10 = -t22 * t243 + t23 * t242;
t9 = -t20 * t243 + t21 * t242;
t8 = -t18 * t243 + t19 * t242;
t7 = -t16 * t243 + t17 * t242;
t1 = [0; -m(3) * t153 + m(4) * t68 + m(5) * t34 + m(6) * t28; (t39 * t28 + t58 * t98 + t59 * t99) * t348 + (t135 * t77 + t136 * t78 + t51 * t34) * t349 + (t150 * t185 + t151 * t186 + t68 * t79) * t350 + 0.2e1 * m(3) * (-t153 + t355) * t356 * (rSges(3,1) * t247 - rSges(3,2) * t245) + (-t346 * t214 - t12 - t7 - t8) * t243 + (t347 * t215 + t10 + t13 + t9 + (-t242 * t214 + t243 * t215) * t243) * t242; m(4) * t57 + m(5) * t29 + m(6) * t11; m(6) * (t11 * t39 + t14 * t99 + t15 * t98 + t40 * t28 + t49 * t59 + t50 * t58) + m(5) * (t135 * t31 + t136 * t30 + t29 * t51 + t56 * t34 + t66 * t78 + t67 * t77) + m(4) * (t101 * t68 + t138 * t151 + t139 * t150 + t185 * t65 + t186 * t64 + t57 * t79) + ((t9 / 0.2e1 + t10 / 0.2e1 + t13 / 0.2e1) * t243 + (t7 / 0.2e1 + t8 / 0.2e1 + t12 / 0.2e1) * t242) * t245 + ((-t386 * t243 + (t72 + t391) * t242) * t364 / 0.2e1 + (t367 * t242 - t368 * t243) * t245 / 0.2e1) * qJD(2) + (((-t69 - t392) * t243 + t387 * t242) * t305 + t362) * t242 / 0.2e1 - t363 * t243 / 0.2e1 - (t360 * t242 + t361 * t243) * t247 / 0.2e1; (t11 * t40 + t14 * t49 + t15 * t50) * t348 + (t56 * t29 + t30 * t66 + t31 * t67) * t349 + (t101 * t57 + t138 * t64 + t139 * t65) * t350 + t363 * t326 + t362 * t323 + (t353 * t306 + (t242 * t69 + t371) * t294 + (t243 * t72 + t372) * t293) * t245 + ((t370 * t198 + t369 * t199 + t378 * t222 + t377 * t223 - t330) * t247 + t354 * t306 - t376 * t294 + t375 * t293 + (-t247 * t258 + t321 * t357 - t353) * t305 + ((t179 * t246 - t210 * t304 + t374 * t244 + t302 * t393) * t247 - t360 * t243 + t361 * t242 + (t354 + t329) * qJD(2)) * t245) * t247; (m(5) + m(6)) * t248; m(6) * (t39 * t291 - t194 * t98 - t196 * t99 + t220 * t58 + t224 * t59 + (t245 * t28 + t305 * t39) * t244) + m(5) * (t51 * t291 - t194 * t135 - t196 * t136 + t220 * t77 + t224 * t78 + (t245 * t34 + t305 * t51) * t244); m(6) * (t40 * t291 + t224 * t14 + t220 * t15 - t194 * t50 - t196 * t49 + (t11 * t245 + t305 * t40) * t244) + m(5) * (t56 * t291 - t194 * t67 - t196 * t66 + t220 * t31 + t224 * t30 + (t245 * t29 + t305 * t56) * t244); 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t220 * t194 - t224 * t196 + t248 * t322); m(6) * t198; m(6) * (t154 * t98 + t156 * t99 + t190 * t58 + t192 * t59 + t198 * t39 + t222 * t28); m(6) * (t11 * t222 + t14 * t192 + t15 * t190 + t154 * t50 + t156 * t49 + t198 * t40); m(6) * (t154 * t220 + t156 * t224 - t190 * t194 - t192 * t196 + t198 * t322 + t222 * t248); (t154 * t190 + t156 * t192 + t198 * t222) * t348;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
