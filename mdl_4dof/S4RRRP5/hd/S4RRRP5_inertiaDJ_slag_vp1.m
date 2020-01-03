% Calculate time derivative of joint inertia matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:50
% DurationCPUTime: 6.81s
% Computational Cost: add. (5797->434), mult. (7811->600), div. (0->0), fcn. (6098->6), ass. (0->252)
t191 = qJ(2) + qJ(3);
t183 = cos(t191);
t182 = sin(t191);
t308 = Icges(5,5) * t182;
t146 = -Icges(5,3) * t183 + t308;
t313 = Icges(4,4) * t182;
t149 = Icges(4,2) * t183 + t313;
t367 = -t149 + t146;
t307 = Icges(5,5) * t183;
t150 = Icges(5,1) * t182 - t307;
t312 = Icges(4,4) * t183;
t151 = Icges(4,1) * t182 + t312;
t366 = t151 + t150;
t329 = rSges(5,1) + pkin(3);
t365 = rSges(5,3) + qJ(4);
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t228 = Icges(5,3) * t182 + t307;
t108 = Icges(5,6) * t193 + t195 * t228;
t232 = -Icges(4,2) * t182 + t312;
t114 = Icges(4,6) * t193 + t195 * t232;
t364 = t108 - t114;
t229 = Icges(4,5) * t183 - Icges(4,6) * t182;
t110 = Icges(4,3) * t193 + t195 * t229;
t231 = Icges(5,4) * t183 + Icges(5,6) * t182;
t112 = Icges(5,2) * t193 + t195 * t231;
t363 = t110 + t112;
t235 = Icges(5,1) * t183 + t308;
t116 = Icges(5,4) * t193 + t195 * t235;
t236 = Icges(4,1) * t183 - t313;
t118 = Icges(4,5) * t193 + t195 * t236;
t362 = t116 + t118;
t188 = qJD(2) + qJD(3);
t361 = (-t228 + t232) * t188;
t360 = (t235 + t236) * t188;
t147 = Icges(4,5) * t182 + Icges(4,6) * t183;
t148 = Icges(5,4) * t182 - Icges(5,6) * t183;
t354 = t147 + t148;
t351 = t367 * t182 + t366 * t183;
t107 = -Icges(5,6) * t195 + t193 * t228;
t113 = -Icges(4,6) * t195 + t193 * t232;
t359 = -t107 + t113;
t115 = -Icges(5,4) * t195 + t193 * t235;
t117 = -Icges(4,5) * t195 + t193 * t236;
t358 = t115 + t117;
t109 = -Icges(4,3) * t195 + t193 * t229;
t111 = -Icges(5,2) * t195 + t193 * t231;
t227 = t107 * t182 + t115 * t183;
t339 = t195 * t227;
t225 = t113 * t182 - t117 * t183;
t340 = t195 * t225;
t357 = -t339 + t340 + (-t109 - t111) * t193;
t224 = t114 * t182 - t118 * t183;
t226 = t108 * t182 + t116 * t183;
t356 = (-t224 + t226) * t195 + t363 * t193;
t186 = t193 * rSges(5,2);
t288 = t183 * t195;
t290 = t182 * t195;
t355 = t329 * t288 + t365 * t290 + t186;
t153 = rSges(5,1) * t182 - rSges(5,3) * t183;
t353 = pkin(3) * t182 - qJ(4) * t183 + t153;
t292 = t151 * t188;
t293 = t150 * t188;
t294 = t149 * t188;
t295 = t146 * t188;
t352 = (t295 - t294 + t360) * t183 + (-t293 - t292 - t361) * t182 + t354 * qJD(1);
t63 = -t107 * qJD(1) - t195 * t295;
t69 = -t113 * qJD(1) - t195 * t294;
t350 = -t362 * t188 + t63 - t69;
t71 = -t115 * qJD(1) - t195 * t293;
t73 = -t117 * qJD(1) - t195 * t292;
t349 = t364 * t188 + t71 + t73;
t196 = -pkin(6) - pkin(5);
t192 = sin(qJ(2));
t269 = qJD(2) * t192;
t263 = pkin(2) * t269;
t348 = qJD(1) * t196 + t263;
t347 = (-t229 - t231) * t188 + t351 * qJD(1);
t194 = cos(qJ(2));
t165 = rSges(3,1) * t192 + rSges(3,2) * t194;
t214 = qJD(2) * t165;
t346 = t193 * t214;
t314 = Icges(3,4) * t194;
t234 = -Icges(3,2) * t192 + t314;
t137 = Icges(3,6) * t193 + t195 * t234;
t315 = Icges(3,4) * t192;
t238 = Icges(3,1) * t194 - t315;
t139 = Icges(3,5) * t193 + t195 * t238;
t222 = t137 * t192 - t139 * t194;
t345 = t193 * t222;
t344 = t193 * t224;
t343 = t193 * t226;
t178 = pkin(2) * t194 + pkin(1);
t324 = pkin(1) - t178;
t342 = t193 * t324;
t136 = -Icges(3,6) * t195 + t193 * t234;
t138 = -Icges(3,5) * t195 + t193 * t238;
t223 = t136 * t192 - t138 * t194;
t341 = t195 * t223;
t338 = qJD(1) * t109;
t337 = qJD(1) * t111;
t230 = Icges(3,5) * t194 - Icges(3,6) * t192;
t134 = -Icges(3,3) * t195 + t193 * t230;
t334 = 2 * m(3);
t333 = 2 * m(4);
t332 = 2 * m(5);
t189 = t193 ^ 2;
t190 = t195 ^ 2;
t331 = t193 / 0.2e1;
t330 = -t195 / 0.2e1;
t328 = m(3) * t165;
t154 = rSges(4,1) * t182 + rSges(4,2) * t183;
t327 = m(4) * t154;
t326 = pkin(2) * t192;
t325 = t193 * pkin(5);
t187 = t195 * pkin(5);
t323 = -pkin(5) - t196;
t285 = t195 * t196;
t105 = t187 + t285 - t342;
t168 = t195 * t178;
t106 = -t195 * pkin(1) + t193 * t323 + t168;
t322 = t193 * t105 + t195 * t106;
t321 = rSges(3,1) * t194;
t320 = rSges(4,1) * t183;
t319 = rSges(3,2) * t192;
t318 = rSges(3,3) * t195;
t185 = t193 * rSges(3,3);
t184 = t193 * rSges(4,3);
t241 = rSges(5,1) * t183 + rSges(5,3) * t182;
t255 = pkin(3) * t188 - qJD(4);
t291 = t182 * t188;
t316 = -qJ(4) * t291 - t183 * t255 - t241 * t188;
t242 = -rSges(4,2) * t182 + t320;
t133 = t242 * t188;
t300 = t133 * t193;
t299 = t136 * t194;
t298 = t137 * t194;
t297 = t138 * t192;
t296 = t139 * t192;
t289 = t183 * t188;
t287 = t188 * t193;
t286 = t188 * t195;
t121 = -t195 * rSges(4,3) + t193 * t242;
t123 = rSges(4,1) * t288 - rSges(4,2) * t290 + t184;
t62 = t193 * t121 + t195 * t123;
t272 = qJD(1) * t193;
t283 = t353 * t272;
t261 = t183 * t286;
t267 = qJD(4) * t182;
t281 = qJ(4) * t261 + t195 * t267;
t271 = qJD(1) * t195;
t280 = rSges(5,2) * t271 + rSges(5,3) * t261;
t260 = t182 * t272;
t279 = rSges(4,2) * t260 + rSges(4,3) * t271;
t278 = t348 * t193;
t277 = t195 * t321 + t185;
t276 = t189 + t190;
t275 = qJD(1) * t110;
t274 = qJD(1) * t112;
t135 = Icges(3,3) * t193 + t195 * t230;
t273 = qJD(1) * t135;
t268 = qJD(2) * t194;
t204 = -t182 * t286 - t183 * t272;
t215 = t154 * t188;
t266 = t193 * (-t193 * t215 + (t195 * t242 + t184) * qJD(1)) + t195 * (rSges(4,1) * t204 - rSges(4,2) * t261 + t279) + t121 * t271;
t265 = t193 * ((-t195 * t324 - t325) * qJD(1) - t278) + t195 * (-t195 * t263 + (t195 * t323 + t342) * qJD(1)) + t105 * t271;
t264 = t195 * t319;
t262 = pkin(2) * t268;
t259 = t192 * t272;
t256 = -t154 - t326;
t72 = qJD(1) * t116 - t193 * t293;
t254 = t107 * t188 + t72;
t74 = qJD(1) * t118 - t193 * t292;
t252 = t113 * t188 - t74;
t64 = qJD(1) * t108 - t193 * t295;
t250 = t115 * t188 - t64;
t70 = qJD(1) * t114 - t193 * t294;
t248 = t117 * t188 + t70;
t90 = t353 * t195;
t246 = -t193 * t196 + t168;
t120 = -t195 * rSges(5,2) + t193 * t241;
t142 = (pkin(3) * t183 + qJ(4) * t182) * t193;
t27 = t355 * t195 + (t120 + t142) * t193;
t245 = -t353 - t326;
t30 = -t111 * t195 + t193 * t227;
t31 = -t112 * t195 + t343;
t32 = -t109 * t195 - t193 * t225;
t33 = -t110 * t195 - t344;
t209 = t188 * t147;
t65 = -t195 * t209 - t338;
t66 = -t193 * t209 + t275;
t210 = t188 * t148;
t67 = -t195 * t210 - t337;
t68 = -t193 * t210 + t274;
t244 = ((-t30 - t32) * t272 + t357 * t271) * t195 + (((t67 + t65) * t193 + (-t343 + t344 - t357) * qJD(1)) * t193 + (t31 + t33) * t272 + t356 * t271 + ((-t66 - t68) * t193 + (t359 * t289 + t358 * t291 - t337 - t338) * t195 + (t349 * t193 + (-t72 - t74) * t195) * t183 + (t350 * t193 + (-t64 + t70) * t195) * t182 + ((t227 - t225 + t363) * t193 + t356) * qJD(1)) * t195) * t193;
t243 = -t319 + t321;
t237 = Icges(3,1) * t192 + t314;
t233 = Icges(3,2) * t194 + t315;
t219 = t142 * t271 + t193 * ((pkin(3) * t271 + qJ(4) * t287) * t183 + (qJ(4) * t271 - t193 * t255) * t182) + t195 * (pkin(3) * t204 - qJ(4) * t260 + t281) + t193 * (-t153 * t287 + (t195 * t241 + t186) * qJD(1)) + t195 * (rSges(5,1) * t204 - rSges(5,3) * t260 + t280) + t120 * t271;
t218 = -t262 + t316;
t217 = -pkin(1) - t243;
t86 = t245 * t195;
t216 = -t178 - t242;
t207 = qJD(2) * t237;
t206 = qJD(2) * t233;
t205 = qJD(2) * (-Icges(3,5) * t192 - Icges(3,6) * t194);
t3 = (t195 * t68 + (t31 - t339) * qJD(1)) * t195 + (t30 * qJD(1) + (t108 * t289 - t116 * t291 + t182 * t63 + t183 * t71 + t274) * t193 + (-t67 - t254 * t183 + t250 * t182 + (-t111 + t226) * qJD(1)) * t195) * t193;
t4 = (t195 * t66 + (t33 + t340) * qJD(1)) * t195 + (t32 * qJD(1) + (-t114 * t289 - t118 * t291 - t182 * t69 + t183 * t73 + t275) * t193 + (-t65 + t252 * t183 + t248 * t182 + (-t109 - t224) * qJD(1)) * t195) * t193;
t203 = (-t4 - t3) * t195 + t244;
t202 = -t182 * t365 - t183 * t329 - t178;
t201 = rSges(3,2) * t259 + rSges(3,3) * t271 - t195 * t214;
t200 = t202 * t193;
t199 = (t349 * t182 - t350 * t183 - t347 * t193 + t352 * t195) * t331 + (t347 * t195 + t352 * t193 + (t248 + t250) * t183 + (-t252 + t254) * t182) * t330 + (t358 * t182 + t359 * t183 + t351 * t193 - t354 * t195) * t272 / 0.2e1 + (t362 * t182 - t183 * t364 + t354 * t193 + t351 * t195) * t271 / 0.2e1;
t175 = pkin(2) * t259;
t160 = t243 * qJD(2);
t141 = -t264 + t277;
t140 = t193 * t243 - t318;
t104 = t256 * t195;
t103 = t256 * t193;
t95 = t325 + (pkin(1) - t319) * t195 + t277;
t94 = t193 * t217 + t187 + t318;
t89 = t353 * t193;
t88 = t123 + t246;
t87 = (rSges(4,3) - t196) * t195 + t216 * t193;
t85 = t245 * t193;
t80 = t193 * t205 + t273;
t79 = -qJD(1) * t134 + t195 * t205;
t76 = t346 + ((-rSges(3,3) - pkin(5)) * t193 + t217 * t195) * qJD(1);
t75 = (t187 + (-pkin(1) - t321) * t193) * qJD(1) + t201;
t57 = t246 + t355;
t56 = (rSges(5,2) - t196) * t195 + t200;
t53 = -t154 * t271 - t300 + (-t192 * t271 - t193 * t268) * pkin(2);
t52 = t154 * t272 + t175 + (-t133 - t262) * t195;
t43 = t135 * t193 - t222 * t195;
t42 = t134 * t193 - t341;
t41 = -t135 * t195 - t345;
t40 = -t134 * t195 - t193 * t223;
t39 = t154 * t287 + (t195 * t216 - t184) * qJD(1) + t278;
t38 = (-t178 - t320) * t272 + (-t215 - t348) * t195 + t279;
t29 = -qJD(1) * t90 + t193 * t316;
t28 = t195 * t316 + t283;
t26 = qJD(1) * t86 + t193 * t218;
t25 = t195 * t218 + t175 + t283;
t24 = t62 + t322;
t23 = (-t267 + (t182 * t329 - t183 * t365) * t188) * t193 + (t195 * t202 - t186) * qJD(1) + t278;
t22 = (-t291 * t329 - t263) * t195 + (t200 - t285) * qJD(1) + t280 + t281;
t21 = t27 + t322;
t20 = -t123 * t272 + t266;
t7 = (-t106 - t123) * t272 + t265 + t266;
t6 = -t272 * t355 + t219;
t5 = (-t106 - t355) * t272 + t219 + t265;
t1 = [(t75 * t95 + t76 * t94) * t334 + (t38 * t88 + t39 * t87) * t333 + (t22 * t57 + t23 * t56) * t332 + t367 * t291 + t366 * t289 + (t238 - t233) * t269 + (t237 + t234) * t268 + t361 * t183 + t360 * t182; (t189 / 0.2e1 + t190 / 0.2e1) * t230 * qJD(2) + m(3) * ((-t193 * t75 - t195 * t76) * t165 + (-t193 * t95 - t195 * t94) * t160) + m(4) * (t103 * t38 + t104 * t39 + t52 * t87 + t53 * t88) + m(5) * (t22 * t85 + t23 * t86 + t25 * t56 + t26 * t57) + (-qJD(2) * t222 + t192 * (-t138 * qJD(1) - t195 * t207) + t194 * (-t136 * qJD(1) - t195 * t206)) * t331 + (-qJD(2) * t223 + t192 * (qJD(1) * t139 - t193 * t207) + t194 * (qJD(1) * t137 - t193 * t206)) * t330 + ((-t95 * t328 + t298 / 0.2e1 + t296 / 0.2e1) * t195 + (t299 / 0.2e1 + t297 / 0.2e1 + t94 * t328) * t193) * qJD(1) + t199; (t21 * t5 + t25 * t86 + t26 * t85) * t332 - t195 * t3 - t195 * t4 + (t103 * t53 + t104 * t52 + t24 * t7) * t333 + (t193 * t43 - t195 * t42) * t271 + t193 * ((t193 * t79 + (t42 + t345) * qJD(1)) * t193 + (t43 * qJD(1) + (t136 * t268 + t138 * t269) * t195 + (-t80 + (-t296 - t298) * qJD(2) + (t135 - t223) * qJD(1)) * t193) * t195) + (t193 * t41 - t195 * t40) * t272 - t195 * ((t195 * t80 + (t41 + t341) * qJD(1)) * t195 + (t40 * qJD(1) + (-t137 * t268 - t139 * t269 + t273) * t193 + (-t79 + (t297 + t299) * qJD(2) - t222 * qJD(1)) * t195) * t193) + ((t140 * t193 + t141 * t195) * ((qJD(1) * t140 + t201) * t195 + (-t346 + (-t141 - t264 + t185) * qJD(1)) * t193) + t276 * t165 * t160) * t334 + t244; m(4) * (-t193 * t88 - t195 * t87) * t133 + m(5) * (-t22 * t89 - t23 * t90 + t28 * t56 + t29 * t57) + (-t193 * t38 - t195 * t39 + (t193 * t87 - t195 * t88) * qJD(1)) * t327 + t199; m(5) * (t6 * t21 - t25 * t90 - t26 * t89 + t27 * t5 + t28 * t86 + t29 * t85) + m(4) * (-t104 * t133 * t195 - t103 * t300 + t20 * t24 + t62 * t7) + (-t193 * t53 - t195 * t52 + (-t103 * t195 + t104 * t193) * qJD(1)) * t327 + t203; (t133 * t154 * t276 + t20 * t62) * t333 + (t27 * t6 - t28 * t90 - t29 * t89) * t332 + t203; m(5) * ((t193 * t57 + t195 * t56) * t289 + (t193 * t22 + t195 * t23 + (-t193 * t56 + t195 * t57) * qJD(1)) * t182); m(5) * ((-t5 + (t193 * t85 + t195 * t86) * t188) * t183 + (t188 * t21 + t193 * t26 + t195 * t25 + (-t193 * t86 + t195 * t85) * qJD(1)) * t182); m(5) * ((-t6 + (-t193 * t89 - t195 * t90) * t188) * t183 + (t188 * t27 + t193 * t29 + t195 * t28 + (t193 * t90 - t195 * t89) * qJD(1)) * t182); (-0.1e1 + t276) * t182 * t289 * t332;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
