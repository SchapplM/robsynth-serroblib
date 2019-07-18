% Calculate time derivative of joint inertia matrix for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:24:55
% EndTime: 2019-07-18 13:25:15
% DurationCPUTime: 8.94s
% Computational Cost: add. (20761->716), mult. (57358->1064), div. (0->0), fcn. (61082->8), ass. (0->326)
t241 = sin(qJ(3));
t244 = cos(qJ(4));
t313 = qJD(4) * t244;
t297 = t241 * t313;
t240 = sin(qJ(4));
t245 = cos(qJ(3));
t317 = qJD(3) * t245;
t301 = t240 * t317;
t250 = t297 + t301;
t246 = cos(qJ(1));
t316 = qJD(3) * t246;
t300 = t245 * t316;
t242 = sin(qJ(1));
t323 = qJD(1) * t242;
t302 = t241 * t323;
t251 = t300 - t302;
t270 = Icges(5,5) * t244 - Icges(5,6) * t240;
t189 = -Icges(5,3) * t245 + t241 * t270;
t357 = Icges(5,4) * t244;
t272 = -Icges(5,2) * t240 + t357;
t192 = -Icges(5,6) * t245 + t241 * t272;
t358 = Icges(5,4) * t240;
t275 = Icges(5,1) * t244 - t358;
t195 = -Icges(5,5) * t245 + t241 * t275;
t331 = t244 * t246;
t333 = t242 * t245;
t206 = t240 * t333 + t331;
t207 = -t240 * t246 + t244 * t333;
t337 = t241 * t242;
t115 = t189 * t337 - t192 * t206 + t195 * t207;
t239 = sin(qJ(5));
t243 = cos(qJ(5));
t336 = t241 * t243;
t179 = -t207 * t239 + t242 * t336;
t180 = t207 * t243 + t239 * t337;
t124 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t206;
t126 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t206;
t128 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t206;
t54 = t124 * t206 + t126 * t179 + t128 * t180;
t330 = t245 * t246;
t303 = t244 * t330;
t209 = t242 * t240 + t303;
t334 = t241 * t246;
t181 = -t209 * t239 + t243 * t334;
t182 = t209 * t243 + t239 * t334;
t208 = t240 * t330 - t242 * t244;
t125 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t208;
t127 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t208;
t129 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t208;
t55 = t125 * t206 + t127 * t179 + t129 * t180;
t282 = t242 * t54 + t246 * t55;
t332 = t243 * t245;
t335 = t241 * t244;
t204 = -t239 * t335 - t332;
t205 = -t239 * t245 + t243 * t335;
t338 = t240 * t241;
t155 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t338;
t158 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t338;
t161 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t338;
t73 = t155 * t206 + t158 * t179 + t161 * t180;
t23 = t241 * t282 - t73 * t245;
t156 = Icges(5,5) * t207 - Icges(5,6) * t206 + Icges(5,3) * t337;
t159 = Icges(5,4) * t207 - Icges(5,2) * t206 + Icges(5,6) * t337;
t162 = Icges(5,1) * t207 - Icges(5,4) * t206 + Icges(5,5) * t337;
t81 = t156 * t337 - t159 * t206 + t162 * t207;
t157 = Icges(5,5) * t209 - Icges(5,6) * t208 + Icges(5,3) * t334;
t160 = Icges(5,4) * t209 - Icges(5,2) * t208 + Icges(5,6) * t334;
t163 = Icges(5,1) * t209 - Icges(5,4) * t208 + Icges(5,5) * t334;
t82 = t157 * t337 - t160 * t206 + t163 * t207;
t279 = t242 * t81 + t246 * t82;
t386 = -t115 * t245 + t241 * t279 + t23;
t116 = t189 * t334 - t208 * t192 + t209 * t195;
t56 = t124 * t208 + t126 * t181 + t128 * t182;
t57 = t125 * t208 + t127 * t181 + t129 * t182;
t281 = t242 * t56 + t246 * t57;
t74 = t155 * t208 + t158 * t181 + t161 * t182;
t24 = t241 * t281 - t74 * t245;
t83 = t156 * t334 - t208 * t159 + t209 * t162;
t84 = t157 * t334 - t208 * t160 + t209 * t163;
t278 = t242 * t83 + t246 * t84;
t385 = -t116 * t245 + t241 * t278 + t24;
t360 = Icges(4,4) * t241;
t277 = Icges(4,1) * t245 - t360;
t197 = Icges(4,5) * t242 + t246 * t277;
t339 = t197 * t245;
t359 = Icges(4,4) * t245;
t274 = -Icges(4,2) * t241 + t359;
t194 = Icges(4,6) * t242 + t246 * t274;
t345 = t194 * t241;
t262 = -t339 + t345;
t384 = t242 * t262;
t196 = -Icges(4,5) * t246 + t242 * t277;
t341 = t196 * t245;
t193 = -Icges(4,6) * t246 + t242 * t274;
t347 = t193 * t241;
t263 = -t341 + t347;
t383 = t246 * t263;
t382 = -rSges(4,2) * t334 + t242 * rSges(4,3);
t200 = rSges(4,1) * t330 + t382;
t234 = t242 * qJ(2);
t185 = t200 + t234;
t285 = rSges(4,1) * t245 - rSges(4,2) * t241;
t199 = -rSges(4,3) * t246 + t242 * t285;
t235 = t246 * qJ(2);
t186 = -t199 + t235;
t381 = t185 * t242 + t186 * t246;
t131 = t182 * rSges(6,1) + t181 * rSges(6,2) + t208 * rSges(6,3);
t121 = t234 + t131;
t130 = rSges(6,1) * t180 + rSges(6,2) * t179 + rSges(6,3) * t206;
t122 = -t130 + t235;
t380 = t121 * t242 + t122 * t246;
t315 = qJD(4) * t240;
t319 = qJD(3) * t242;
t321 = qJD(1) * t246;
t153 = -t319 * t338 - t246 * t315 - t244 * t323 + (t240 * t321 + t242 * t313) * t245;
t261 = (-qJD(4) * t245 + qJD(1)) * t240;
t322 = qJD(1) * t245;
t294 = -qJD(4) + t322;
t318 = qJD(3) * t244;
t154 = t294 * t331 + (-t241 * t318 + t261) * t242;
t252 = t241 * t321 + t242 * t317;
t108 = t154 * rSges(5,1) - t153 * rSges(5,2) + rSges(5,3) * t252;
t271 = Icges(4,5) * t245 - Icges(4,6) * t241;
t190 = -Icges(4,3) * t246 + t242 * t271;
t379 = 2 * m(4);
t378 = 2 * m(5);
t377 = 2 * m(6);
t237 = t242 ^ 2;
t238 = t246 ^ 2;
t299 = t241 * t316;
t151 = qJD(1) * t206 - qJD(4) * t303 + t240 * t299 - t242 * t315;
t376 = -t151 / 0.2e1;
t375 = t153 / 0.2e1;
t374 = t206 / 0.2e1;
t373 = t208 / 0.2e1;
t372 = t240 / 0.2e1;
t371 = t241 / 0.2e1;
t370 = t242 / 0.2e1;
t369 = -t245 / 0.2e1;
t368 = t246 / 0.2e1;
t221 = rSges(4,1) * t241 + rSges(4,2) * t245;
t367 = m(4) * t221;
t283 = rSges(5,1) * t244 - rSges(5,2) * t240;
t198 = -rSges(5,3) * t245 + t241 * t283;
t366 = m(5) * t198;
t164 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t338;
t365 = m(6) * t164;
t364 = qJD(1) / 0.2e1;
t363 = qJD(3) / 0.2e1;
t362 = rSges(3,1) * t246;
t314 = qJD(4) * t241;
t176 = (-rSges(5,1) * t240 - rSges(5,2) * t244) * t314 + (rSges(5,3) * t241 + t245 * t283) * qJD(3);
t351 = t176 * t242;
t350 = t176 * t246;
t346 = t193 * t245;
t344 = t194 * t245;
t343 = t195 * t244;
t342 = t196 * t241;
t340 = t197 * t241;
t329 = rSges(4,2) * t302 + rSges(4,3) * t321;
t328 = qJ(2) * t321 + qJD(2) * t242;
t327 = t237 + t238;
t326 = qJD(1) * t164;
t191 = Icges(4,3) * t242 + t246 * t271;
t325 = qJD(1) * t191;
t324 = qJD(1) * t198;
t320 = qJD(3) * t241;
t312 = qJD(5) * t241;
t292 = -qJD(5) + t318;
t293 = -qJD(5) * t244 + qJD(3);
t298 = t240 * t314;
t149 = t293 * t336 + (-t245 * t292 + t298) * t239;
t150 = t292 * t332 + (t239 * t293 - t243 * t315) * t241;
t100 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t250;
t103 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t250;
t97 = Icges(6,5) * t150 + Icges(6,6) * t149 + Icges(6,3) * t250;
t32 = t204 * t100 + t205 * t103 + t149 * t158 + t150 * t161 + t155 * t250 + t97 * t338;
t80 = t155 * t338 + t158 * t204 + t161 * t205;
t311 = t250 * t80 + t32 * t338;
t249 = -qJD(5) * t207 + t252;
t287 = t242 * t312 + t154;
t113 = -t239 * t287 + t243 * t249;
t114 = t239 * t249 + t243 * t287;
t66 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t153;
t68 = Icges(6,4) * t114 + Icges(6,2) * t113 + Icges(6,6) * t153;
t70 = Icges(6,1) * t114 + Icges(6,4) * t113 + Icges(6,5) * t153;
t18 = t124 * t250 + t126 * t149 + t128 * t150 + t204 * t68 + t205 * t70 + t338 * t66;
t26 = t100 * t179 + t103 * t180 + t113 * t158 + t114 * t161 + t153 * t155 + t206 * t97;
t309 = t18 / 0.2e1 + t26 / 0.2e1;
t248 = -qJD(5) * t209 + t251;
t152 = t246 * t261 + (-t242 * t294 - t299) * t244;
t288 = t246 * t312 + t152;
t111 = -t239 * t288 + t243 * t248;
t112 = t239 * t248 + t243 * t288;
t65 = Icges(6,5) * t112 + Icges(6,6) * t111 - Icges(6,3) * t151;
t67 = Icges(6,4) * t112 + Icges(6,2) * t111 - Icges(6,6) * t151;
t69 = Icges(6,1) * t112 + Icges(6,4) * t111 - Icges(6,5) * t151;
t19 = t125 * t250 + t127 * t149 + t129 * t150 + t204 * t67 + t205 * t69 + t338 * t65;
t25 = t100 * t181 + t103 * t182 + t111 * t158 + t112 * t161 - t151 * t155 + t208 * t97;
t308 = t25 / 0.2e1 + t19 / 0.2e1;
t39 = t55 * t242 - t246 * t54;
t51 = t82 * t242 - t246 * t81;
t307 = t39 / 0.2e1 + t51 / 0.2e1;
t40 = t57 * t242 - t246 * t56;
t52 = t84 * t242 - t246 * t83;
t306 = t40 / 0.2e1 + t52 / 0.2e1;
t59 = t124 * t338 + t126 * t204 + t128 * t205;
t305 = t59 / 0.2e1 + t73 / 0.2e1;
t60 = t125 * t338 + t127 * t204 + t129 * t205;
t304 = -t74 / 0.2e1 - t60 / 0.2e1;
t71 = t112 * rSges(6,1) + t111 * rSges(6,2) - t151 * rSges(6,3);
t166 = t209 * rSges(5,1) - t208 * rSges(5,2) + rSges(5,3) * t334;
t45 = t60 * t242 - t59 * t246;
t296 = t45 * t363;
t295 = t313 / 0.2e1;
t268 = -t159 * t240 + t162 * t244;
t86 = -t156 * t245 + t241 * t268;
t291 = -t245 * t86 + t386;
t267 = -t160 * t240 + t163 * t244;
t87 = -t157 * t245 + t241 * t267;
t290 = t245 * t87 - t385;
t233 = qJD(2) * t246;
t289 = -qJ(2) * t323 + t233;
t286 = -t242 * rSges(3,1) + rSges(3,3) * t246;
t280 = t59 * t242 + t60 * t246;
t276 = Icges(4,1) * t241 + t359;
t273 = Icges(4,2) * t245 + t360;
t269 = t130 * t246 - t131 * t242;
t85 = t242 * t130 + t131 * t246;
t165 = rSges(5,1) * t207 - rSges(5,2) * t206 + rSges(5,3) * t337;
t266 = t165 * t246 - t166 * t242;
t118 = t242 * t165 + t166 * t246;
t102 = Icges(5,4) * t154 - Icges(5,2) * t153 + Icges(5,6) * t252;
t105 = Icges(5,1) * t154 - Icges(5,4) * t153 + Icges(5,5) * t252;
t99 = Icges(5,5) * t154 - Icges(5,6) * t153 + Icges(5,3) * t252;
t41 = (qJD(3) * t268 - t99) * t245 + (qJD(3) * t156 - t102 * t240 + t105 * t244 + (-t159 * t244 - t162 * t240) * qJD(4)) * t241;
t167 = (-Icges(5,5) * t240 - Icges(5,6) * t244) * t314 + (Icges(5,3) * t241 + t245 * t270) * qJD(3);
t170 = (-Icges(5,2) * t244 - t358) * t314 + (Icges(5,6) * t241 + t245 * t272) * qJD(3);
t173 = (-Icges(5,1) * t240 - t357) * t314 + (Icges(5,5) * t241 + t245 * t275) * qJD(3);
t50 = -t153 * t192 + t154 * t195 + t167 * t337 - t206 * t170 + t207 * t173 + t189 * t252;
t260 = t41 / 0.2e1 + t50 / 0.2e1 + t309;
t101 = Icges(5,4) * t152 + Icges(5,2) * t151 + Icges(5,6) * t251;
t104 = Icges(5,1) * t152 + Icges(5,4) * t151 + Icges(5,5) * t251;
t98 = Icges(5,5) * t152 + Icges(5,6) * t151 + Icges(5,3) * t251;
t42 = (qJD(3) * t267 - t98) * t245 + (qJD(3) * t157 - t101 * t240 + t104 * t244 + (-t160 * t244 - t163 * t240) * qJD(4)) * t241;
t49 = t151 * t192 + t152 * t195 + t167 * t334 - t208 * t170 + t209 * t173 + t189 * t251;
t259 = t42 / 0.2e1 + t49 / 0.2e1 + t308;
t258 = t86 / 0.2e1 + t115 / 0.2e1 + t305;
t257 = -t116 / 0.2e1 - t87 / 0.2e1 + t304;
t256 = qJD(3) * t221;
t255 = qJD(3) * t276;
t254 = qJD(3) * t273;
t253 = qJD(3) * (-Icges(4,5) * t241 - Icges(4,6) * t245);
t72 = rSges(6,1) * t114 + rSges(6,2) * t113 + rSges(6,3) * t153;
t107 = t152 * rSges(5,1) + t151 * rSges(5,2) + rSges(5,3) * t251;
t247 = -t245 * t167 + t173 * t335 + t189 * t320 - t192 * t297 + t317 * t343;
t215 = t285 * qJD(3);
t211 = t235 + t286;
t210 = t242 * rSges(3,3) + t234 + t362;
t188 = qJD(1) * t286 + t328;
t187 = t233 + (-t362 + (-rSges(3,3) - qJ(2)) * t242) * qJD(1);
t169 = t242 * t253 + t325;
t168 = -qJD(1) * t190 + t246 * t253;
t148 = -t165 + t235;
t147 = t234 + t166;
t140 = -rSges(4,2) * t300 + (-t242 * t322 - t299) * rSges(4,1) + t328 + t329;
t139 = t233 + t221 * t319 + (-t285 * t246 + (-rSges(4,3) - qJ(2)) * t242) * qJD(1);
t138 = -t245 * t166 - t198 * t334;
t137 = t165 * t245 + t198 * t337;
t136 = t242 * t191 - t246 * t262;
t135 = t242 * t190 - t383;
t134 = -t191 * t246 - t384;
t133 = -t190 * t246 - t242 * t263;
t132 = -t189 * t245 + (-t192 * t240 + t343) * t241;
t123 = t132 * t320;
t117 = t266 * t241;
t106 = rSges(6,1) * t150 + rSges(6,2) * t149 + rSges(6,3) * t250;
t95 = t107 + t328;
t94 = -t108 + t289;
t91 = -t245 * t131 - t164 * t334;
t90 = t130 * t245 + t164 * t337;
t89 = t131 * t338 - t164 * t208;
t88 = -t130 * t338 + t164 * t206;
t79 = t269 * t241;
t78 = t80 * t320;
t75 = t130 * t208 - t131 * t206;
t64 = t71 + t328;
t63 = t289 - t72;
t62 = (t198 * t319 + t108) * t245 + (-qJD(3) * t165 + t198 * t321 + t351) * t241;
t61 = (-t198 * t316 - t107) * t245 + (qJD(3) * t166 + t198 * t323 - t350) * t241;
t58 = (-t192 * t317 + (-qJD(4) * t195 - t170) * t241) * t240 + t247;
t53 = qJD(1) * t266 + t107 * t246 + t242 * t108;
t48 = t266 * t317 + (-qJD(1) * t118 - t107 * t242 + t108 * t246) * t241;
t44 = (t164 * t319 + t72) * t245 + (-qJD(3) * t130 + t106 * t242 + t164 * t321) * t241;
t43 = (-t164 * t316 - t71) * t245 + (qJD(3) * t131 - t106 * t246 + t164 * t323) * t241;
t38 = t106 * t206 - t130 * t250 + t153 * t164 - t338 * t72;
t37 = -t106 * t208 + t131 * t250 + t151 * t164 + t338 * t71;
t36 = -t206 * t101 + t207 * t104 - t153 * t160 + t154 * t163 + t157 * t252 + t337 * t98;
t35 = -t206 * t102 + t207 * t105 - t153 * t159 + t154 * t162 + t156 * t252 + t337 * t99;
t34 = -t208 * t101 + t209 * t104 + t151 * t160 + t152 * t163 + t157 * t251 + t334 * t98;
t33 = -t208 * t102 + t209 * t105 + t151 * t159 + t152 * t162 + t156 * t251 + t334 * t99;
t30 = qJD(1) * t269 + t242 * t72 + t246 * t71;
t29 = t241 * t280 - t80 * t245;
t28 = t206 * t59 + t208 * t60 + t338 * t80;
t27 = -t130 * t151 - t131 * t153 - t206 * t71 + t208 * t72;
t22 = t269 * t317 + (-qJD(1) * t85 - t242 * t71 + t246 * t72) * t241;
t21 = t206 * t56 + t208 * t57 + t338 * t74;
t20 = t206 * t54 + t208 * t55 + t338 * t73;
t17 = t113 * t127 + t114 * t129 + t125 * t153 + t179 * t67 + t180 * t69 + t206 * t65;
t16 = t113 * t126 + t114 * t128 + t124 * t153 + t179 * t68 + t180 * t70 + t206 * t66;
t15 = t111 * t127 + t112 * t129 - t125 * t151 + t181 * t67 + t182 * t69 + t208 * t65;
t14 = t111 * t126 + t112 * t128 - t124 * t151 + t181 * t68 + t182 * t70 + t208 * t66;
t13 = qJD(1) * t279 + t36 * t242 - t246 * t35;
t12 = qJD(1) * t278 + t34 * t242 - t246 * t33;
t11 = (qJD(3) * t279 - t50) * t245 + (-qJD(1) * t51 + qJD(3) * t115 + t242 * t35 + t246 * t36) * t241;
t10 = (qJD(3) * t278 - t49) * t245 + (-qJD(1) * t52 + qJD(3) * t116 + t242 * t33 + t246 * t34) * t241;
t9 = qJD(1) * t280 - t18 * t246 + t19 * t242;
t8 = qJD(1) * t282 - t16 * t246 + t17 * t242;
t7 = qJD(1) * t281 - t14 * t246 + t15 * t242;
t6 = -t60 * t151 + t59 * t153 + t18 * t206 + t19 * t208 + t311;
t5 = t78 + (qJD(3) * t280 - t32) * t245 + (-qJD(1) * t45 + t18 * t242 + t19 * t246) * t241;
t4 = t73 * t297 - t151 * t55 + t153 * t54 + t16 * t206 + t17 * t208 + (t241 * t26 + t317 * t73) * t240;
t3 = t74 * t297 + t14 * t206 + t15 * t208 - t151 * t57 + t153 * t56 + (t241 * t25 + t317 * t74) * t240;
t2 = (qJD(3) * t282 - t26) * t245 + (-qJD(1) * t39 + qJD(3) * t73 + t16 * t242 + t17 * t246) * t241;
t1 = (qJD(3) * t281 - t25) * t245 + (-qJD(1) * t40 + qJD(3) * t74 + t14 * t242 + t15 * t246) * t241;
t31 = [t32 + t247 + (t121 * t64 + t122 * t63) * t377 + (t147 * t95 + t148 * t94) * t378 + (t139 * t186 + t140 * t185) * t379 + 0.2e1 * m(3) * (t187 * t211 + t188 * t210) - t170 * t338 - t192 * t301 - t195 * t298 + (t277 - t273) * t320 + (t274 + t276) * t317; m(6) * (qJD(1) * t380 + t242 * t63 - t246 * t64) + m(5) * (t242 * t94 - t246 * t95 + (t147 * t242 + t148 * t246) * qJD(1)) + m(4) * (qJD(1) * t381 + t242 * t139 - t140 * t246) + m(3) * (t242 * t187 - t188 * t246 + (t210 * t242 + t211 * t246) * qJD(1)); 0; ((qJD(1) * t194 - t242 * t254) * t369 - (qJD(1) * t197 - t242 * t255) * t241 / 0.2e1 + (t347 / 0.2e1 - t341 / 0.2e1) * qJD(3) - t260) * t246 + ((-qJD(1) * t193 - t246 * t254) * t245 / 0.2e1 + (-qJD(1) * t196 - t246 * t255) * t371 + (-t345 / 0.2e1 + t339 / 0.2e1) * qJD(3) + t259) * t242 + m(5) * (-t147 * t351 - t148 * t350 + (-t242 * t95 - t246 * t94) * t198) + m(6) * ((-t242 * t64 - t246 * t63) * t164 - t380 * t106) + m(4) * ((-t139 * t246 - t140 * t242) * t221 - t381 * t215) + (t238 / 0.2e1 + t237 / 0.2e1) * t271 * qJD(3) + ((t344 / 0.2e1 + t340 / 0.2e1 - t185 * t367 - t147 * t366 - t121 * t365 - t257) * t246 + (t346 / 0.2e1 + t342 / 0.2e1 + t186 * t367 + t148 * t366 + t122 * t365 + t258) * t242) * qJD(1); 0; (t106 * t164 * t327 + t30 * t85) * t377 - t246 * t8 + t242 * t7 + (t176 * t198 * t327 + t118 * t53) * t378 - t246 * t13 + t242 * t12 + t242 * ((t242 * t168 + (t135 + t384) * qJD(1)) * t242 + (t136 * qJD(1) + (t193 * t317 + t196 * t320) * t246 + (-t169 + (-t340 - t344) * qJD(3) + (t191 - t263) * qJD(1)) * t242) * t246) - t246 * ((t246 * t169 + (t134 + t383) * qJD(1)) * t246 + (t133 * qJD(1) + (-t194 * t317 - t197 * t320 + t325) * t242 + (-t168 + (t342 + t346) * qJD(3) - t262 * qJD(1)) * t246) * t242) + ((t242 * t199 + t200 * t246) * ((qJD(1) * t199 - t246 * t256 + t329) * t246 + (-t242 * t256 + (-t200 + t382) * qJD(1)) * t242) + t327 * t221 * t215) * t379 + (-t133 * t246 + t134 * t242 + t39 + t51) * t323 + (-t135 * t246 + t136 * t242 + t40 + t52) * t321; t78 + t123 + m(5) * (t137 * t94 + t138 * t95 + t147 * t61 + t148 * t62) + m(6) * (t121 * t43 + t122 * t44 + t63 * t90 + t64 * t91) + (-t32 - t58 + (t242 * t258 - t246 * t257) * qJD(3)) * t245 + (t259 * t246 + t260 * t242 + (t242 * t257 + t246 * t258) * qJD(1)) * t241; m(5) * (t62 * t242 - t246 * t61 + (t137 * t246 + t138 * t242) * qJD(1)) + m(6) * (t44 * t242 - t246 * t43 + (t242 * t91 + t246 * t90) * qJD(1)); t9 * t369 + m(5) * (t117 * t53 + t48 * t118) + m(6) * (t22 * t85 + t79 * t30) + ((qJD(1) * t87 - t41) * t369 - t11 / 0.2e1 - t2 / 0.2e1 + t306 * t317 + m(5) * (-t137 * t176 - t138 * t324 - t198 * t62) + m(6) * (-t106 * t90 - t164 * t44 - t326 * t91)) * t246 + ((qJD(1) * t86 + t42) * t369 + t10 / 0.2e1 + t1 / 0.2e1 + t307 * t317 + m(5) * (t137 * t324 - t138 * t176 - t198 * t61) + m(6) * (-t106 * t91 - t164 * t43 + t326 * t90)) * t242 + ((t87 * t242 - t246 * t86) * t363 + t296 + (-t242 * t306 + t246 * t307) * qJD(1) + (t13 + t8) * t370 + (t12 + t7) * t368) * t241 + (t242 * t386 + t385 * t246) * t364; (t22 * t79 + t43 * t91 + t44 * t90) * t377 + (t117 * t48 + t137 * t62 + t138 * t61) * t378 + (t58 * t245 - t123 - t5 + (t242 * t291 - t246 * t290) * qJD(3)) * t245 + ((-t245 * t42 + t1 + t10) * t246 + (-t245 * t41 + t11 + t2) * t242 + (t29 - t132 * t245 + (t242 * t86 + t246 * t87) * t241) * qJD(3) + (t242 * t290 + t246 * t291) * qJD(1)) * t241; m(6) * (t121 * t37 + t122 * t38 + t63 * t88 + t64 * t89) + t308 * t208 + t309 * t206 + t305 * t153 + t304 * t151 + t311; m(6) * (t38 * t242 - t246 * t37 + (t242 * t89 + t246 * t88) * qJD(1)); m(6) * (t27 * t85 + t75 * t30) + t39 * t375 + t8 * t374 + t40 * t376 + t7 * t373 + t241 * t45 * t295 + (m(6) * (-t106 * t88 - t164 * t38 - t326 * t89) - t4 / 0.2e1 + t21 * t364) * t246 + (m(6) * (-t106 * t89 - t164 * t37 + t326 * t88) + t20 * t364 + t3 / 0.2e1) * t242 + (t245 * t296 + t371 * t9) * t240; m(6) * (t22 * t75 + t27 * t79 + t37 * t91 + t38 * t90 + t43 * t89 + t44 * t88) + t23 * t375 + t2 * t374 + t24 * t376 + t1 * t373 + (-t6 / 0.2e1 + (t20 * t370 + t21 * t368 + t29 * t372) * qJD(3)) * t245 + (t28 * t363 + t29 * t295 + t5 * t372 + t3 * t368 + t4 * t370 + (-t242 * t21 / 0.2e1 + t20 * t368) * qJD(1)) * t241; (t27 * t75 + t37 * t89 + t38 * t88) * t377 - t151 * t21 + t208 * t3 + t153 * t20 + t206 * t4 + t28 * t297 + (t241 * t6 + t28 * t317) * t240;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t31(1), t31(2), t31(4), t31(7), t31(11); t31(2), t31(3), t31(5), t31(8), t31(12); t31(4), t31(5), t31(6), t31(9), t31(13); t31(7), t31(8), t31(9), t31(10), t31(14); t31(11), t31(12), t31(13), t31(14), t31(15);];
Mq  = res;
