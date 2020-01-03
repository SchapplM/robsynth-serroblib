% Calculate time derivative of joint inertia matrix for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:22
% DurationCPUTime: 6.72s
% Computational Cost: add. (5454->425), mult. (7342->573), div. (0->0), fcn. (5726->6), ass. (0->257)
t191 = qJ(2) + qJ(3);
t181 = cos(t191);
t180 = sin(t191);
t312 = Icges(5,4) * t180;
t149 = Icges(5,2) * t181 + t312;
t314 = Icges(4,4) * t180;
t150 = Icges(4,2) * t181 + t314;
t372 = -t150 - t149;
t311 = Icges(5,4) * t181;
t151 = Icges(5,1) * t180 + t311;
t313 = Icges(4,4) * t181;
t152 = Icges(4,1) * t180 + t313;
t371 = t152 + t151;
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t228 = Icges(5,5) * t181 - Icges(5,6) * t180;
t112 = Icges(5,3) * t193 + t195 * t228;
t229 = Icges(4,5) * t181 - Icges(4,6) * t180;
t114 = Icges(4,3) * t193 + t195 * t229;
t370 = t112 + t114;
t231 = -Icges(5,2) * t180 + t311;
t116 = Icges(5,6) * t193 + t195 * t231;
t232 = -Icges(4,2) * t180 + t313;
t118 = Icges(4,6) * t193 + t195 * t232;
t369 = t116 + t118;
t235 = Icges(5,1) * t181 - t312;
t120 = Icges(5,5) * t193 + t195 * t235;
t236 = Icges(4,1) * t181 - t314;
t122 = Icges(4,5) * t193 + t195 * t236;
t368 = t120 + t122;
t188 = qJD(2) + qJD(3);
t367 = (t231 + t232) * t188;
t366 = (t235 + t236) * t188;
t147 = Icges(5,5) * t180 + Icges(5,6) * t181;
t148 = Icges(4,5) * t180 + Icges(4,6) * t181;
t355 = t147 + t148;
t365 = -t372 * t180 - t371 * t181;
t115 = -Icges(5,6) * t195 + t193 * t231;
t117 = -Icges(4,6) * t195 + t193 * t232;
t364 = t115 + t117;
t119 = -Icges(5,5) * t195 + t193 * t235;
t121 = -Icges(4,5) * t195 + t193 * t236;
t363 = t119 + t121;
t194 = cos(qJ(2));
t175 = t194 * pkin(2) + pkin(1);
t159 = pkin(3) * t181 + t175;
t182 = t193 * rSges(5,3);
t289 = t181 * t195;
t291 = t180 * t195;
t362 = rSges(5,1) * t289 - rSges(5,2) * t291 + t195 * t159 + t182;
t192 = sin(qJ(2));
t270 = qJD(2) * t192;
t264 = pkin(2) * t270;
t292 = t180 * t188;
t145 = -pkin(3) * t292 - t264;
t153 = rSges(5,1) * t180 + rSges(5,2) * t181;
t215 = t153 * t188;
t361 = qJD(4) * t195 - (t145 - t215) * t193;
t111 = -Icges(5,3) * t195 + t193 * t228;
t113 = -Icges(4,3) * t195 + t193 * t229;
t227 = t115 * t180 - t119 * t181;
t340 = t195 * t227;
t225 = t117 * t180 - t121 * t181;
t341 = t195 * t225;
t360 = t340 + t341 + (-t111 - t113) * t193;
t224 = t118 * t180 - t122 * t181;
t226 = t116 * t180 - t120 * t181;
t359 = (-t224 - t226) * t195 + t370 * t193;
t167 = t195 * t175;
t196 = -pkin(6) - pkin(5);
t187 = -qJ(4) + t196;
t278 = t187 - t196;
t358 = -t193 * t278 - t167 + t362;
t322 = rSges(5,1) * t181;
t242 = -rSges(5,2) * t180 + t322;
t283 = t159 - t175;
t90 = t193 * t283 + t195 * t278;
t357 = -rSges(5,3) * t195 + t193 * t242 + t90;
t293 = t152 * t188;
t294 = t151 * t188;
t295 = t150 * t188;
t296 = t149 * t188;
t354 = (-t295 - t296 + t366) * t181 + (-t293 - t294 - t367) * t180 + t355 * qJD(1);
t65 = -t115 * qJD(1) - t195 * t296;
t67 = -t117 * qJD(1) - t195 * t295;
t352 = -t368 * t188 - t65 - t67;
t69 = -t119 * qJD(1) - t195 * t294;
t71 = -t121 * qJD(1) - t195 * t293;
t351 = -t369 * t188 + t69 + t71;
t350 = qJD(1) * t196 + t264;
t273 = qJD(1) * t193;
t261 = t180 * t273;
t272 = qJD(1) * t195;
t349 = rSges(5,2) * t261 + rSges(5,3) * t272 + qJD(4) * t193 + t195 * t145;
t348 = (t228 + t229) * t188 + t365 * qJD(1);
t166 = rSges(3,1) * t192 + rSges(3,2) * t194;
t213 = qJD(2) * t166;
t347 = t193 * t213;
t315 = Icges(3,4) * t194;
t234 = -Icges(3,2) * t192 + t315;
t138 = Icges(3,6) * t193 + t195 * t234;
t316 = Icges(3,4) * t192;
t238 = Icges(3,1) * t194 - t316;
t140 = Icges(3,5) * t193 + t195 * t238;
t222 = t138 * t192 - t140 * t194;
t346 = t193 * t222;
t345 = t193 * t224;
t344 = t193 * t226;
t326 = pkin(1) - t175;
t343 = t193 * t326;
t137 = -Icges(3,6) * t195 + t193 * t234;
t139 = -Icges(3,5) * t195 + t193 * t238;
t223 = t137 * t192 - t139 * t194;
t342 = t195 * t223;
t339 = qJD(1) * t111;
t338 = qJD(1) * t113;
t230 = Icges(3,5) * t194 - Icges(3,6) * t192;
t135 = -Icges(3,3) * t195 + t193 * t230;
t335 = 2 * m(3);
t334 = 2 * m(4);
t333 = 2 * m(5);
t189 = t193 ^ 2;
t190 = t195 ^ 2;
t332 = t193 / 0.2e1;
t331 = -t195 / 0.2e1;
t330 = m(3) * t166;
t154 = rSges(4,1) * t180 + rSges(4,2) * t181;
t329 = m(4) * t154;
t328 = pkin(2) * t192;
t327 = t193 * pkin(5);
t186 = t195 * pkin(5);
t325 = -pkin(5) - t196;
t324 = rSges(3,1) * t194;
t323 = rSges(4,1) * t181;
t321 = rSges(3,2) * t192;
t320 = rSges(3,3) * t195;
t184 = t193 * rSges(3,3);
t183 = t193 * rSges(4,3);
t319 = rSges(5,3) - t187;
t109 = t195 * t196 + t186 - t343;
t110 = -t195 * pkin(1) + t193 * t325 + t167;
t318 = t193 * t109 + t195 * t110;
t243 = -rSges(4,2) * t180 + t323;
t134 = t243 * t188;
t301 = t134 * t193;
t300 = t137 * t194;
t299 = t138 * t194;
t298 = t139 * t192;
t297 = t140 * t192;
t290 = t181 * t188;
t288 = t188 * t193;
t287 = t188 * t195;
t286 = t193 * t187;
t124 = -t195 * rSges(4,3) + t193 * t243;
t126 = rSges(4,1) * t289 - rSges(4,2) * t291 + t183;
t60 = t193 * t124 + t195 * t126;
t284 = pkin(3) * t261 + t153 * t273;
t281 = rSges(4,2) * t261 + rSges(4,3) * t272;
t280 = t350 * t193;
t279 = t195 * t324 + t184;
t277 = t189 + t190;
t276 = qJD(1) * t112;
t275 = qJD(1) * t114;
t136 = Icges(3,3) * t193 + t195 * t230;
t274 = qJD(1) * t136;
t269 = qJD(2) * t194;
t203 = -t180 * t287 - t181 * t273;
t214 = t154 * t188;
t262 = t181 * t287;
t267 = t193 * (-t193 * t214 + (t195 * t243 + t183) * qJD(1)) + t195 * (rSges(4,1) * t203 - rSges(4,2) * t262 + t281) + t124 * t272;
t247 = t195 * t264;
t266 = t193 * ((-t195 * t326 - t327) * qJD(1) - t280) + t195 * (-t247 + (t195 * t325 + t343) * qJD(1)) + t109 * t272;
t265 = t195 * t321;
t263 = pkin(2) * t269;
t260 = t192 * t273;
t257 = -t154 - t328;
t256 = -pkin(3) * t180 - t153;
t70 = qJD(1) * t120 - t193 * t294;
t255 = t115 * t188 - t70;
t72 = qJD(1) * t122 - t193 * t293;
t253 = t117 * t188 - t72;
t66 = qJD(1) * t116 - t193 * t296;
t251 = t119 * t188 + t66;
t68 = qJD(1) * t118 - t193 * t295;
t249 = t121 * t188 + t68;
t22 = t193 * t357 + t195 * t358;
t133 = t242 * t188;
t246 = -pkin(3) * t290 - t133;
t26 = -t111 * t195 - t193 * t227;
t27 = -t112 * t195 - t344;
t28 = -t113 * t195 - t193 * t225;
t29 = -t114 * t195 - t345;
t207 = t188 * t147;
t61 = -t195 * t207 - t339;
t62 = -t193 * t207 + t276;
t208 = t188 * t148;
t63 = -t195 * t208 - t338;
t64 = -t193 * t208 + t275;
t245 = ((-t26 - t28) * t273 + t360 * t272) * t195 + (((t61 + t63) * t193 + (t344 + t345 - t360) * qJD(1)) * t193 + (t27 + t29) * t273 + t359 * t272 + ((-t62 - t64) * t193 + (t364 * t290 + t363 * t292 - t338 - t339) * t195 + (t351 * t193 + (-t70 - t72) * t195) * t181 + (t352 * t193 + (t66 + t68) * t195) * t180 + ((-t227 - t225 + t370) * t193 + t359) * qJD(1)) * t195) * t193;
t244 = -t321 + t324;
t239 = t357 * t272 + (rSges(5,1) * t203 - rSges(5,2) * t262 - t90 * qJD(1) + t247 + t349) * t195 + (t280 + (t182 - t286 + (t242 + t283) * t195) * qJD(1) - t361) * t193;
t237 = Icges(3,1) * t192 + t315;
t233 = Icges(3,2) * t194 + t316;
t219 = -pkin(1) - t244;
t218 = t256 - t328;
t217 = -t175 - t243;
t216 = -t159 - t242;
t206 = qJD(2) * t237;
t205 = qJD(2) * t233;
t204 = qJD(2) * (-Icges(3,5) * t192 - Icges(3,6) * t194);
t93 = t218 * t195;
t202 = t246 - t263;
t3 = (t195 * t62 + (t27 + t340) * qJD(1)) * t195 + (t26 * qJD(1) + (-t116 * t290 - t120 * t292 - t180 * t65 + t181 * t69 + t276) * t193 + (-t61 + t255 * t181 + t251 * t180 + (-t111 - t226) * qJD(1)) * t195) * t193;
t4 = (t195 * t64 + (t29 + t341) * qJD(1)) * t195 + (t28 * qJD(1) + (-t118 * t290 - t122 * t292 - t180 * t67 + t181 * t71 + t275) * t193 + (-t63 + t253 * t181 + t249 * t180 + (-t113 - t224) * qJD(1)) * t195) * t193;
t201 = (-t4 - t3) * t195 + t245;
t200 = rSges(3,2) * t260 + rSges(3,3) * t272 - t195 * t213;
t199 = (t180 * t351 - t352 * t181 + t193 * t348 + t354 * t195) * t332 + (-t348 * t195 + t354 * t193 + (t249 + t251) * t181 + (-t253 - t255) * t180) * t331 + (t363 * t180 + t364 * t181 - t193 * t365 - t355 * t195) * t273 / 0.2e1 + (t368 * t180 + t369 * t181 + t355 * t193 - t195 * t365) * t272 / 0.2e1;
t172 = pkin(2) * t260;
t158 = t244 * qJD(2);
t143 = -t265 + t279;
t142 = t193 * t244 - t320;
t108 = t257 * t195;
t107 = t257 * t193;
t106 = t256 * t195;
t105 = t256 * t193;
t98 = t327 + (pkin(1) - t321) * t195 + t279;
t97 = t193 * t219 + t186 + t320;
t92 = t218 * t193;
t89 = -t193 * t196 + t126 + t167;
t88 = (rSges(4,3) - t196) * t195 + t217 * t193;
t83 = t193 * t204 + t274;
t82 = -qJD(1) * t135 + t195 * t204;
t76 = -t286 + t362;
t75 = t193 * t216 + t195 * t319;
t74 = t347 + ((-rSges(3,3) - pkin(5)) * t193 + t219 * t195) * qJD(1);
t73 = (t186 + (-pkin(1) - t324) * t193) * qJD(1) + t200;
t55 = -t154 * t272 - t301 + (-t192 * t272 - t193 * t269) * pkin(2);
t54 = t154 * t273 + t172 + (-t134 - t263) * t195;
t49 = -t153 * t272 - t133 * t193 + (-t180 * t272 - t181 * t288) * pkin(3);
t48 = t195 * t246 + t284;
t43 = t136 * t193 - t222 * t195;
t42 = t135 * t193 - t342;
t41 = -t136 * t195 - t346;
t40 = -t135 * t195 - t193 * t223;
t39 = t154 * t288 + (t195 * t217 - t183) * qJD(1) + t280;
t38 = (-t175 - t323) * t273 + (-t214 - t350) * t195 + t281;
t35 = qJD(1) * t93 + t193 * t202;
t34 = t195 * t202 + t172 + t284;
t25 = (-t193 * t319 + t195 * t216) * qJD(1) + t361;
t24 = -t195 * t215 + (-t195 * t187 + (-t159 - t322) * t193) * qJD(1) + t349;
t23 = t60 + t318;
t21 = -t126 * t273 + t267;
t16 = t22 + t318;
t7 = (-t110 - t126) * t273 + t266 + t267;
t6 = -t273 * t358 + t239;
t5 = (-t110 - t358) * t273 + t239 + t266;
t1 = [(t73 * t98 + t74 * t97) * t335 + (t38 * t89 + t39 * t88) * t334 + (t24 * t76 + t25 * t75) * t333 + t372 * t292 + t371 * t290 + (t238 - t233) * t270 + (t237 + t234) * t269 + t367 * t181 + t366 * t180; m(3) * ((-t193 * t73 - t195 * t74) * t166 + (-t193 * t98 - t195 * t97) * t158) + ((-t98 * t330 + t299 / 0.2e1 + t297 / 0.2e1) * t195 + (t300 / 0.2e1 + t298 / 0.2e1 + t97 * t330) * t193) * qJD(1) + (-qJD(2) * t222 + t192 * (-t139 * qJD(1) - t195 * t206) + t194 * (-t137 * qJD(1) - t195 * t205)) * t332 + (-qJD(2) * t223 + t192 * (qJD(1) * t140 - t193 * t206) + t194 * (qJD(1) * t138 - t193 * t205)) * t331 + m(4) * (t107 * t38 + t108 * t39 + t54 * t88 + t55 * t89) + m(5) * (t24 * t92 + t25 * t93 + t34 * t75 + t35 * t76) + (t189 / 0.2e1 + t190 / 0.2e1) * t230 * qJD(2) + t199; (t16 * t5 + t34 * t93 + t35 * t92) * t333 - t195 * t4 - t195 * t3 + (t107 * t55 + t108 * t54 + t23 * t7) * t334 + (t193 * t43 - t195 * t42) * t272 + t193 * ((t193 * t82 + (t42 + t346) * qJD(1)) * t193 + (t43 * qJD(1) + (t137 * t269 + t139 * t270) * t195 + (-t83 + (-t297 - t299) * qJD(2) + (t136 - t223) * qJD(1)) * t193) * t195) + ((t142 * t193 + t143 * t195) * ((qJD(1) * t142 + t200) * t195 + (-t347 + (-t143 - t265 + t184) * qJD(1)) * t193) + t277 * t166 * t158) * t335 + (t193 * t41 - t195 * t40) * t273 - t195 * ((t195 * t83 + (t41 + t342) * qJD(1)) * t195 + (t40 * qJD(1) + (-t138 * t269 - t140 * t270 + t274) * t193 + (-t82 + (t298 + t300) * qJD(2) - t222 * qJD(1)) * t195) * t193) + t245; (-t193 * t38 - t195 * t39 + (t193 * t88 - t195 * t89) * qJD(1)) * t329 + m(4) * (-t193 * t89 - t195 * t88) * t134 + m(5) * (t105 * t24 + t106 * t25 + t48 * t75 + t49 * t76) + t199; m(5) * (t105 * t35 + t106 * t34 + t16 * t6 + t22 * t5 + t48 * t93 + t49 * t92) + m(4) * (-t108 * t134 * t195 - t107 * t301 + t21 * t23 + t60 * t7) + (-t193 * t55 - t195 * t54 + (-t107 * t195 + t108 * t193) * qJD(1)) * t329 + t201; (t134 * t154 * t277 + t21 * t60) * t334 + (t105 * t49 + t106 * t48 + t22 * t6) * t333 + t201; m(5) * (t193 * t25 - t195 * t24 + (t193 * t76 + t195 * t75) * qJD(1)); m(5) * (t193 * t34 - t195 * t35 + (t193 * t92 + t195 * t93) * qJD(1)); m(5) * (t193 * t48 - t195 * t49 + (t105 * t193 + t106 * t195) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
