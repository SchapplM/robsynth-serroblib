% Calculate time derivative of joint inertia matrix for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:19
% EndTime: 2019-12-31 18:21:33
% DurationCPUTime: 7.54s
% Computational Cost: add. (15652->612), mult. (16292->880), div. (0->0), fcn. (15361->10), ass. (0->302)
t207 = qJ(1) + pkin(8);
t204 = cos(t207);
t202 = sin(t207);
t211 = sin(qJ(3));
t213 = cos(qJ(3));
t336 = Icges(4,4) * t213;
t249 = -Icges(4,2) * t211 + t336;
t130 = -Icges(4,6) * t204 + t202 * t249;
t337 = Icges(4,4) * t211;
t253 = Icges(4,1) * t213 - t337;
t132 = -Icges(4,5) * t204 + t202 * t253;
t240 = t130 * t211 - t132 * t213;
t229 = t240 * t204;
t323 = t204 * t211;
t209 = cos(pkin(9));
t320 = t209 * t213;
t208 = sin(pkin(9));
t324 = t204 * t208;
t157 = t202 * t320 - t324;
t321 = t208 * t213;
t231 = t202 * t321 + t204 * t209;
t327 = t202 * t211;
t94 = Icges(5,5) * t157 - Icges(5,6) * t231 + Icges(5,3) * t327;
t380 = t94 * t323 - t229;
t131 = Icges(4,6) * t202 + t204 * t249;
t133 = Icges(4,5) * t202 + t204 * t253;
t239 = t131 * t211 - t133 * t213;
t228 = t239 * t202;
t158 = t202 * t209 - t204 * t321;
t328 = t202 * t208;
t159 = t204 * t320 + t328;
t95 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t323;
t379 = -t95 * t327 + t228;
t355 = t202 / 0.2e1;
t378 = -t204 / 0.2e1;
t377 = -qJD(1) / 0.2e1;
t245 = Icges(4,5) * t213 - Icges(4,6) * t211;
t129 = Icges(4,3) * t202 + t204 * t245;
t97 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t323;
t99 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t323;
t376 = -t129 * t204 + t157 * t99 - t231 * t97 - t379;
t128 = -Icges(4,3) * t204 + t202 * t245;
t96 = Icges(5,4) * t157 - Icges(5,2) * t231 + Icges(5,6) * t327;
t98 = Icges(5,1) * t157 - Icges(5,4) * t231 + Icges(5,5) * t327;
t375 = -t128 * t202 - t158 * t96 - t159 * t98 - t380;
t374 = t239 * t204 - t95 * t323;
t294 = t94 * t327;
t373 = t128 * t204 - t157 * t98 + t202 * t240 + t231 * t96 - t294;
t372 = t129 * t202 + t158 * t97 + t159 * t99 - t374;
t362 = 2 * m(4);
t344 = rSges(4,1) * t213;
t267 = -rSges(4,2) * t211 + t344;
t343 = rSges(4,3) * t204;
t137 = t202 * t267 - t343;
t322 = t204 * t213;
t365 = -rSges(4,2) * t323 + t202 * rSges(4,3);
t138 = rSges(4,1) * t322 + t365;
t185 = rSges(4,1) * t211 + rSges(4,2) * t213;
t230 = qJD(3) * t185;
t306 = qJD(1) * t211;
t283 = t202 * t306;
t307 = qJD(1) * t204;
t217 = rSges(4,2) * t283 + rSges(4,3) * t307 - t204 * t230;
t38 = (qJD(1) * t137 + t217) * t204 + (-t202 * t230 + (-t138 + t365) * qJD(1)) * t202;
t371 = t362 * t38;
t198 = pkin(4) * t209 + pkin(3);
t210 = -pkin(7) - qJ(4);
t366 = pkin(4) * t328 - t210 * t323;
t206 = pkin(9) + qJ(5);
t201 = sin(t206);
t203 = cos(t206);
t147 = -t201 * t322 + t202 * t203;
t148 = t201 * t202 + t203 * t322;
t90 = t148 * rSges(6,1) + t147 * rSges(6,2) + rSges(6,3) * t323;
t370 = t198 * t322 + t366 + t90;
t345 = pkin(3) - t198;
t369 = t211 * t345;
t318 = qJ(4) + t210;
t368 = t213 * t318;
t205 = cos(qJ(1)) * pkin(1);
t364 = t202 * pkin(6) + t205;
t305 = qJD(1) * t213;
t270 = -qJD(5) + t305;
t302 = qJD(3) * t211;
t278 = t204 * t302;
t363 = t202 * t270 + t278;
t361 = 2 * m(5);
t360 = 2 * m(6);
t199 = t202 ^ 2;
t200 = t204 ^ 2;
t359 = m(5) / 0.2e1;
t358 = m(6) / 0.2e1;
t247 = Icges(5,4) * t209 - Icges(5,2) * t208;
t161 = -Icges(5,6) * t213 + t211 * t247;
t357 = t161 / 0.2e1;
t251 = Icges(5,1) * t209 - Icges(5,4) * t208;
t162 = -Icges(5,5) * t213 + t211 * t251;
t356 = t162 / 0.2e1;
t353 = t204 / 0.2e1;
t352 = -t208 / 0.2e1;
t351 = t209 / 0.2e1;
t350 = -t213 / 0.2e1;
t349 = m(4) * t185;
t348 = sin(qJ(1)) * pkin(1);
t347 = pkin(3) * t213;
t346 = qJD(1) / 0.2e1;
t342 = rSges(6,3) * t211;
t341 = -rSges(5,3) - qJ(4);
t340 = -rSges(6,3) + t210;
t218 = -t211 * t318 - t213 * t345;
t298 = pkin(4) * t324;
t325 = t203 * t204;
t326 = t202 * t213;
t145 = -t201 * t326 - t325;
t146 = -t201 * t204 + t203 * t326;
t263 = -rSges(6,1) * t146 - rSges(6,2) * t145;
t89 = rSges(6,3) * t327 - t263;
t339 = t202 * t218 - t298 + t89;
t191 = pkin(3) * t322;
t165 = qJ(4) * t323 + t191;
t338 = -t165 + t370;
t335 = Icges(6,4) * t201;
t334 = Icges(6,4) * t203;
t250 = Icges(6,1) * t203 - t335;
t151 = -Icges(6,5) * t213 + t211 * t250;
t330 = t151 * t203;
t329 = t198 * t211;
t319 = t210 * t213;
t101 = t159 * rSges(5,1) + t158 * rSges(5,2) + rSges(5,3) * t323;
t317 = -t101 - t165;
t262 = rSges(6,1) * t203 - rSges(6,2) * t201;
t153 = -rSges(6,3) * t213 + t211 * t262;
t316 = t153 + t368 - t369;
t261 = qJ(4) * t211 + t347;
t164 = t261 * t202;
t315 = t202 * t164 + t204 * t165;
t166 = qJD(3) * t261 - qJD(4) * t213;
t264 = rSges(5,1) * t209 - rSges(5,2) * t208;
t314 = -(rSges(5,3) * t211 + t213 * t264) * qJD(3) - t166;
t163 = -rSges(5,3) * t213 + t211 * t264;
t184 = pkin(3) * t211 - qJ(4) * t213;
t313 = -t163 - t184;
t312 = qJD(1) * t298 + t210 * t283;
t300 = qJD(4) * t211;
t183 = t204 * t300;
t193 = pkin(6) * t307;
t311 = t183 + t193;
t310 = t199 + t200;
t309 = qJD(1) * t129;
t308 = qJD(1) * t202;
t304 = qJD(3) * t202;
t303 = qJD(3) * t204;
t301 = qJD(3) * t213;
t299 = qJD(5) * t211;
t277 = t204 * t301;
t174 = qJ(4) * t277;
t280 = t202 * t302;
t179 = pkin(3) * t280;
t279 = t202 * t301;
t282 = t204 * t306;
t221 = t279 + t282;
t297 = t164 * t307 + t202 * (qJ(4) * t221 + qJD(1) * t191 + t202 * t300 - t179) + t204 * (-qJ(4) * t283 + t174 + t183 + (-t202 * t305 - t278) * pkin(3));
t271 = -qJD(5) * t213 + qJD(1);
t237 = t271 * t203;
t79 = t201 * t363 + t204 * t237;
t238 = t271 * t201;
t80 = -t203 * t363 + t204 * t238;
t296 = t80 * rSges(6,1) + t79 * rSges(6,2) + rSges(6,3) * t277;
t86 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t323;
t88 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t323;
t259 = -t201 * t86 + t203 * t88;
t84 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t323;
t35 = t211 * t259 - t213 * t84;
t243 = Icges(6,5) * t203 - Icges(6,6) * t201;
t149 = -Icges(6,3) * t213 + t211 * t243;
t246 = -Icges(6,2) * t201 + t334;
t150 = -Icges(6,6) * t213 + t211 * t246;
t52 = t147 * t150 + t148 * t151 + t149 * t323;
t290 = t35 / 0.2e1 + t52 / 0.2e1;
t85 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t327;
t87 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t327;
t260 = -t201 * t85 + t203 * t87;
t83 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t327;
t34 = t211 * t260 - t213 * t83;
t51 = t145 * t150 + t146 * t151 + t149 * t327;
t289 = t51 / 0.2e1 + t34 / 0.2e1;
t113 = (-rSges(6,1) * t201 - rSges(6,2) * t203) * t299 + (t213 * t262 + t342) * qJD(3);
t287 = -t218 * qJD(3) - t113 - t166;
t276 = t208 * t302;
t120 = qJD(1) * t231 + t204 * t276;
t275 = t209 * t302;
t121 = -qJD(1) * t157 - t204 * t275;
t286 = t121 * rSges(5,1) + t120 * rSges(5,2) + rSges(5,3) * t277;
t285 = -t184 - t316;
t284 = t204 * pkin(2) + t364;
t281 = t150 * t301;
t274 = t301 / 0.2e1;
t196 = t204 * pkin(6);
t273 = t196 - t348;
t272 = -t198 * t213 - pkin(2);
t117 = t313 * t204;
t74 = t285 * t204;
t81 = t202 * t237 + (-t204 * t270 + t280) * t201;
t82 = t270 * t325 + (-t203 * t302 + t238) * t202;
t269 = t82 * rSges(6,1) + t81 * rSges(6,2);
t122 = qJD(1) * t158 + t202 * t276;
t123 = qJD(1) * t159 - t202 * t275;
t266 = -t123 * rSges(5,1) - t122 * rSges(5,2);
t265 = -rSges(5,1) * t157 + rSges(5,2) * t231;
t26 = t145 * t85 + t146 * t87 + t327 * t83;
t27 = t145 * t86 + t146 * t88 + t327 * t84;
t15 = t202 * t27 - t204 * t26;
t258 = t202 * t26 + t204 * t27;
t28 = t147 * t85 + t148 * t87 + t323 * t83;
t29 = t147 * t86 + t148 * t88 + t323 * t84;
t16 = t202 * t29 - t204 * t28;
t257 = t202 * t28 + t204 * t29;
t256 = t202 * t35 - t204 * t34;
t255 = t202 * t34 + t204 * t35;
t254 = -t202 * t90 + t204 * t89;
t252 = Icges(4,1) * t211 + t336;
t248 = Icges(4,2) * t213 + t337;
t244 = Icges(5,5) * t209 - Icges(5,6) * t208;
t236 = -pkin(2) - t267;
t220 = t277 - t283;
t42 = Icges(6,5) * t80 + Icges(6,6) * t79 + Icges(6,3) * t220;
t234 = t211 * t42 + t301 * t84;
t43 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t221;
t233 = t211 * t43 + t301 * t83;
t227 = qJD(3) * t252;
t226 = qJD(3) * t248;
t225 = qJD(3) * (-Icges(4,5) * t211 - Icges(4,6) * t213);
t224 = t211 * t341 - pkin(2) - t347;
t222 = t211 * t340 + t272;
t110 = (-Icges(6,5) * t201 - Icges(6,6) * t203) * t299 + (Icges(6,3) * t211 + t213 * t243) * qJD(3);
t112 = (-Icges(6,1) * t201 - t334) * t299 + (Icges(6,5) * t211 + t213 * t250) * qJD(3);
t216 = -t213 * t110 + t149 * t302 + t301 * t330 + (t112 * t211 - t150 * t299) * t203;
t215 = t202 * t224 - t348;
t173 = t267 * qJD(3);
t167 = t184 * t308;
t144 = (Icges(5,5) * t211 + t213 * t251) * qJD(3);
t143 = (Icges(5,6) * t211 + t213 * t247) * qJD(3);
t116 = t313 * t202;
t115 = t138 + t284;
t114 = t202 * t236 + t273 + t343;
t111 = (-Icges(6,2) * t203 - t335) * t299 + (Icges(6,6) * t211 + t213 * t246) * qJD(3);
t105 = t202 * t225 + t309;
t104 = -qJD(1) * t128 + t204 * t225;
t100 = rSges(5,3) * t327 - t265;
t76 = t185 * t304 + (-t205 + (-rSges(4,3) - pkin(6)) * t202 + t236 * t204) * qJD(1);
t75 = t193 + (-t348 + (-pkin(2) - t344) * t202) * qJD(1) + t217;
t73 = t285 * t202;
t72 = t284 - t317;
t71 = t196 + t215 + t265;
t70 = -t149 * t213 + (-t150 * t201 + t330) * t211;
t69 = qJD(1) * t117 + t202 * t314;
t68 = t163 * t308 + t204 * t314 + t167;
t67 = t70 * t302;
t66 = -t153 * t323 - t213 * t90;
t65 = t153 * t327 + t213 * t89;
t64 = Icges(5,1) * t123 + Icges(5,4) * t122 + Icges(5,5) * t221;
t63 = Icges(5,1) * t121 + Icges(5,4) * t120 + Icges(5,5) * t220;
t62 = Icges(5,4) * t123 + Icges(5,2) * t122 + Icges(5,6) * t221;
t61 = Icges(5,4) * t121 + Icges(5,2) * t120 + Icges(5,6) * t220;
t54 = t284 + t370;
t53 = t202 * t222 + t263 + t273 + t298;
t50 = t254 * t211;
t49 = rSges(6,3) * t221 + t269;
t48 = -rSges(6,3) * t283 + t296;
t47 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t221;
t46 = Icges(6,1) * t80 + Icges(6,4) * t79 + Icges(6,5) * t220;
t45 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t221;
t44 = Icges(6,4) * t80 + Icges(6,2) * t79 + Icges(6,6) * t220;
t41 = t100 * t202 + t101 * t204 + t315;
t40 = t179 + (t301 * t341 - t300) * t202 + (t204 * t224 - t364) * qJD(1) + t266;
t39 = -pkin(3) * t278 + qJD(1) * t215 + t174 + t286 + t311;
t37 = qJD(1) * t74 + t202 * t287;
t36 = t204 * t287 + t308 * t316 + t167;
t25 = (-t300 + (t213 * t340 + t329) * qJD(3)) * t202 + (-t205 + (-pkin(4) * t208 - pkin(6)) * t202 + t222 * t204) * qJD(1) - t269;
t24 = (-t319 - t329) * t303 + (-t348 + (t272 - t342) * t202) * qJD(1) + t296 + t311 + t312;
t23 = t202 * t339 + t204 * t338 + t315;
t22 = (-t281 + (-qJD(5) * t151 - t111) * t211) * t201 + t216;
t21 = (t153 * t304 + t49) * t213 + (-qJD(3) * t89 + t113 * t202 + t153 * t307) * t211;
t20 = (-t153 * t303 - t48) * t213 + (qJD(3) * t90 - t113 * t204 + t153 * t308) * t211;
t19 = t110 * t327 + t111 * t145 + t112 * t146 + t149 * t221 + t150 * t81 + t151 * t82;
t18 = t110 * t323 + t111 * t147 + t112 * t148 + t149 * t220 + t150 * t79 + t151 * t80;
t17 = t202 * (rSges(5,3) * t279 - t266) + t204 * t286 + (t204 * t100 + t202 * t317) * qJD(1) + t297;
t14 = t254 * t301 + (-t202 * t48 + t204 * t49 + (-t202 * t89 - t204 * t90) * qJD(1)) * t211;
t13 = t211 * t257 - t213 * t52;
t12 = t211 * t258 - t213 * t51;
t11 = (qJD(3) * t259 - t42) * t213 + (qJD(3) * t84 - t201 * t44 + t203 * t46 + (-t201 * t88 - t203 * t86) * qJD(5)) * t211;
t10 = (qJD(3) * t260 - t43) * t213 + (qJD(3) * t83 - t201 * t45 + t203 * t47 + (-t201 * t87 - t203 * t85) * qJD(5)) * t211;
t9 = t145 * t44 + t146 * t46 + t202 * t234 + t282 * t84 + t81 * t86 + t82 * t88;
t8 = t145 * t45 + t146 * t47 + t202 * t233 + t282 * t83 + t81 * t85 + t82 * t87;
t7 = t147 * t44 + t148 * t46 + t204 * t234 - t283 * t84 + t79 * t86 + t80 * t88;
t6 = t147 * t45 + t148 * t47 + t204 * t233 - t283 * t83 + t79 * t85 + t80 * t87;
t5 = (-t174 + t48 + t312) * t204 + (t179 + t49) * t202 + (t200 * (-t319 + t369) + (-t329 - t368) * t199) * qJD(3) + (t339 * t204 + (-t165 - t338 + t366) * t202) * qJD(1) + t297;
t4 = qJD(1) * t258 + t202 * t9 - t204 * t8;
t3 = qJD(1) * t257 + t202 * t7 - t204 * t6;
t2 = (qJD(3) * t258 - t19) * t213 + (-qJD(1) * t15 + qJD(3) * t51 + t202 * t8 + t204 * t9) * t211;
t1 = (qJD(3) * t257 - t18) * t213 + (-qJD(1) * t16 + qJD(3) * t52 + t202 * t6 + t204 * t7) * t211;
t30 = [(t24 * t54 + t25 * t53) * t360 + (t39 * t72 + t40 * t71) * t361 + (t114 * t76 + t115 * t75) * t362 + t216 + (-t143 * t208 + t144 * t209) * t211 + (-Icges(5,3) * t213 + t211 * t244 - t248 + t253) * t302 + (-t111 * t211 - t151 * t299 - t281) * t201 + (-Icges(5,3) * t211 - t161 * t208 + t162 * t209 - t213 * t244 + t249 + t252) * t301; 0; 0; m(4) * ((-t202 * t75 - t204 * t76) * t185 + (-t114 * t204 - t115 * t202) * t173) + m(5) * (t116 * t39 + t117 * t40 + t68 * t71 + t69 * t72) + m(6) * (t24 * t73 + t25 * t74 + t36 * t53 + t37 * t54) + ((t131 * t377 + t226 * t355 + Icges(5,5) * t123 / 0.2e1 + Icges(5,6) * t122 / 0.2e1 + Icges(5,3) * t221 / 0.2e1) * t204 + (t130 * t377 + t226 * t378 - Icges(5,5) * t121 / 0.2e1 - Icges(5,6) * t120 / 0.2e1 - Icges(5,3) * t220 / 0.2e1) * t202) * t213 + ((t158 * t357 + t159 * t356 - t115 * t349 + (t131 / 0.2e1 - t95 / 0.2e1) * t213 + (t133 / 0.2e1 + t97 * t352 + t99 * t351) * t211 + t290) * t204 + (-t231 * t357 + t157 * t356 + t114 * t349 + (t130 / 0.2e1 - t94 / 0.2e1) * t213 + (t132 / 0.2e1 + t96 * t352 + t98 * t351) * t211 + t289) * t202) * qJD(1) + ((t200 / 0.2e1 + t199 / 0.2e1) * t245 - t228 / 0.2e1 + t229 / 0.2e1) * qJD(3) + (t11 + t18 + t120 * t161 + t121 * t162 + t143 * t158 + t144 * t159 + (-qJD(1) * t132 - t204 * t227 - t208 * t61 + t209 * t63) * t211 + (t211 * t95 + t320 * t99 - t321 * t97) * qJD(3)) * t355 + (t10 + t19 + t122 * t161 + t123 * t162 - t143 * t231 + t144 * t157 + (qJD(1) * t133 - t202 * t227 - t208 * t62 + t209 * t64) * t211 + (t211 * t94 + t320 * t98 - t321 * t96) * qJD(3)) * t378; m(4) * t38 + m(5) * t17 + m(6) * t5; (t23 * t5 + t36 * t74 + t37 * t73) * t360 + t15 * t308 + t16 * t307 + (t116 * t69 + t117 * t68 + t17 * t41) * t361 + t310 * t185 * t173 * t362 + (t137 * t371 + t3 + (t202 * t104 + t120 * t97 + t121 * t99 + t158 * t61 + t159 * t63 + (-t375 + t379) * qJD(1)) * t202 + t376 * t308 + t372 * t307) * t202 + (-t4 + t138 * t371 + (-t204 * t105 + t122 * t96 + t123 * t98 - t231 * t62 + t157 * t64 + (-t376 + t380) * qJD(1)) * t204 + t373 * t308 + t375 * t307 + (-t120 * t96 - t121 * t98 - t158 * t62 - t159 * t64 - t122 * t97 - t123 * t99 + t231 * t61 - t157 * t63 + (t130 * t301 + t132 * t302 + t104 - (t130 * t213 + t132 * t211) * qJD(3)) * t204 + (-t105 + (-t131 * t213 - t133 * t211) * qJD(3) + t131 * t301 + t133 * t302 - t309) * t202 + (t294 + (t129 - t240) * t202 + t372 + t373 + t374) * qJD(1)) * t202) * t204; 0.2e1 * ((t202 * t54 + t204 * t53) * t358 + (t202 * t72 + t204 * t71) * t359) * t301 + 0.2e1 * ((t202 * t24 + t204 * t25 + t307 * t54 - t308 * t53) * t358 + (t202 * t39 + t204 * t40 + t307 * t72 - t308 * t71) * t359) * t211; (m(5) + m(6)) * t302; 0.2e1 * ((t303 * t74 + t304 * t73 - t5) * t358 + (t116 * t304 + t117 * t303 - t17) * t359) * t213 + 0.2e1 * ((qJD(3) * t23 + t202 * t37 + t204 * t36 + t307 * t73 - t308 * t74) * t358 + (qJD(3) * t41 + t116 * t307 - t117 * t308 + t202 * t69 + t204 * t68) * t359) * t211; 0.4e1 * (t359 + t358) * (-0.1e1 + t310) * t211 * t301; m(6) * (t20 * t54 + t21 * t53 + t24 * t66 + t25 * t65) + t67 + (-t22 + (t202 * t289 + t204 * t290) * qJD(3)) * t213 + ((t11 / 0.2e1 + t18 / 0.2e1) * t204 + (t19 / 0.2e1 + t10 / 0.2e1) * t202 + (-t202 * t290 + t204 * t289) * qJD(1)) * t211; m(6) * t14; m(6) * (t14 * t23 + t20 * t73 + t21 * t74 + t36 * t65 + t37 * t66 + t5 * t50) + (-t2 / 0.2e1 + t16 * t274 + (qJD(1) * t35 - t10) * t350 + t13 * t346) * t204 + (t12 * t346 + t15 * t274 + (qJD(1) * t34 + t11) * t350 + t1 / 0.2e1) * t202 + (t3 * t353 + t4 * t355 + qJD(3) * t256 / 0.2e1 + (-t202 * t16 / 0.2e1 + t15 * t353) * qJD(1)) * t211; m(6) * ((-t14 + (t202 * t66 + t204 * t65) * qJD(3)) * t213 + (qJD(3) * t50 + t20 * t202 + t204 * t21 + (-t202 * t65 + t204 * t66) * qJD(1)) * t211); (t14 * t50 + t20 * t66 + t21 * t65) * t360 + (t22 * t213 - t67 + (t202 * t12 + t204 * t13 - t213 * t255) * qJD(3)) * t213 + (t204 * t1 + t202 * t2 - t213 * (t10 * t202 + t11 * t204) + (t211 * t255 - t70 * t213) * qJD(3) + (t204 * t12 - t202 * t13 + t213 * t256) * qJD(1)) * t211;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t30(1), t30(2), t30(4), t30(7), t30(11); t30(2), t30(3), t30(5), t30(8), t30(12); t30(4), t30(5), t30(6), t30(9), t30(13); t30(7), t30(8), t30(9), t30(10), t30(14); t30(11), t30(12), t30(13), t30(14), t30(15);];
Mq = res;
