% Calculate time derivative of joint inertia matrix for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:42
% EndTime: 2019-12-05 15:50:12
% DurationCPUTime: 11.03s
% Computational Cost: add. (45896->745), mult. (130581->1108), div. (0->0), fcn. (155509->12), ass. (0->317)
t333 = sin(pkin(10));
t334 = cos(pkin(10));
t341 = sin(qJ(2));
t343 = cos(qJ(2));
t282 = t341 * t333 - t343 * t334;
t262 = t282 * qJD(2);
t274 = sin(pkin(9));
t276 = cos(pkin(9));
t335 = cos(pkin(5));
t293 = t335 * t333;
t294 = t335 * t334;
t248 = t293 * t343 + t294 * t341;
t280 = qJD(2) * t248;
t206 = t274 * t262 - t276 * t280;
t208 = t262 * t276 + t274 * t280;
t281 = -t293 * t341 + t294 * t343;
t283 = t333 * t343 + t334 * t341;
t227 = -t274 * t283 + t276 * t281;
t229 = -t274 * t281 - t276 * t283;
t275 = sin(pkin(5));
t247 = t283 * t275;
t240 = qJD(2) * t247;
t246 = t282 * t275;
t228 = t248 * t276 - t274 * t282;
t278 = sin(qJ(4));
t329 = t275 * t278;
t342 = cos(qJ(4));
t192 = t228 * t342 - t276 * t329;
t277 = sin(qJ(5));
t279 = cos(qJ(5));
t155 = -t192 * t277 - t227 * t279;
t156 = t192 * t279 - t227 * t277;
t314 = t275 * t342;
t289 = -t228 * t278 - t276 * t314;
t112 = Icges(6,5) * t156 + Icges(6,6) * t155 - Icges(6,3) * t289;
t114 = Icges(6,4) * t156 + Icges(6,2) * t155 - Icges(6,6) * t289;
t116 = Icges(6,1) * t156 + Icges(6,4) * t155 - Icges(6,5) * t289;
t242 = t281 * qJD(2);
t263 = t283 * qJD(2);
t207 = t242 * t276 - t263 * t274;
t150 = qJD(4) * t289 + t207 * t342;
t118 = -qJD(5) * t156 - t150 * t277 - t206 * t279;
t119 = qJD(5) * t155 + t150 * t279 - t206 * t277;
t149 = qJD(4) * t192 + t207 * t278;
t70 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t149;
t72 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t149;
t74 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t149;
t16 = t112 * t149 + t114 * t118 + t116 * t119 + t155 * t72 + t156 * t74 - t289 * t70;
t230 = -t248 * t274 - t276 * t282;
t194 = t230 * t342 + t274 * t329;
t157 = -t194 * t277 - t229 * t279;
t158 = t194 * t279 - t229 * t277;
t290 = -t230 * t278 + t274 * t314;
t113 = Icges(6,5) * t158 + Icges(6,6) * t157 - Icges(6,3) * t290;
t115 = Icges(6,4) * t158 + Icges(6,2) * t157 - Icges(6,6) * t290;
t117 = Icges(6,1) * t158 + Icges(6,4) * t157 - Icges(6,5) * t290;
t209 = -t242 * t274 - t263 * t276;
t152 = qJD(4) * t290 + t209 * t342;
t120 = -qJD(5) * t158 - t152 * t277 - t208 * t279;
t121 = qJD(5) * t157 + t152 * t279 - t208 * t277;
t151 = qJD(4) * t194 + t209 * t278;
t71 = Icges(6,5) * t121 + Icges(6,6) * t120 + Icges(6,3) * t151;
t73 = Icges(6,4) * t121 + Icges(6,2) * t120 + Icges(6,6) * t151;
t75 = Icges(6,1) * t121 + Icges(6,4) * t120 + Icges(6,5) * t151;
t17 = t113 * t149 + t115 * t118 + t117 * t119 + t155 * t73 + t156 * t75 - t289 * t71;
t235 = t247 * t342 + t278 * t335;
t189 = -t235 * t277 + t246 * t279;
t190 = t235 * t279 + t246 * t277;
t288 = -t247 * t278 + t335 * t342;
t136 = Icges(6,5) * t190 + Icges(6,6) * t189 - Icges(6,3) * t288;
t137 = Icges(6,4) * t190 + Icges(6,2) * t189 - Icges(6,6) * t288;
t138 = Icges(6,1) * t190 + Icges(6,4) * t189 - Icges(6,5) * t288;
t318 = qJD(2) * t275;
t241 = t282 * t318;
t188 = qJD(4) * t288 - t241 * t342;
t140 = -qJD(5) * t190 - t188 * t277 + t240 * t279;
t141 = qJD(5) * t189 + t188 * t279 + t240 * t277;
t187 = qJD(4) * t235 - t241 * t278;
t96 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t187;
t97 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t187;
t98 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t187;
t32 = t118 * t137 + t119 * t138 + t136 * t149 + t155 * t97 + t156 * t98 - t289 * t96;
t52 = -t112 * t289 + t114 * t155 + t116 * t156;
t53 = -t113 * t289 + t115 * t155 + t117 * t156;
t68 = -t136 * t289 + t137 * t155 + t138 * t156;
t3 = -t16 * t227 - t17 * t229 - t206 * t52 - t208 * t53 + t240 * t68 + t246 * t32;
t104 = Icges(5,5) * t150 - Icges(5,6) * t149 - Icges(5,3) * t206;
t106 = Icges(5,4) * t150 - Icges(5,2) * t149 - Icges(5,6) * t206;
t108 = Icges(5,1) * t150 - Icges(5,4) * t149 - Icges(5,5) * t206;
t128 = Icges(5,5) * t192 + Icges(5,6) * t289 - Icges(5,3) * t227;
t130 = Icges(5,4) * t192 + Icges(5,2) * t289 - Icges(5,6) * t227;
t132 = Icges(5,1) * t192 + Icges(5,4) * t289 - Icges(5,5) * t227;
t39 = -t104 * t227 + t106 * t289 + t108 * t192 - t128 * t206 - t130 * t149 + t132 * t150;
t105 = Icges(5,5) * t152 - Icges(5,6) * t151 - Icges(5,3) * t208;
t107 = Icges(5,4) * t152 - Icges(5,2) * t151 - Icges(5,6) * t208;
t109 = Icges(5,1) * t152 - Icges(5,4) * t151 - Icges(5,5) * t208;
t129 = Icges(5,5) * t194 + Icges(5,6) * t290 - Icges(5,3) * t229;
t131 = Icges(5,4) * t194 + Icges(5,2) * t290 - Icges(5,6) * t229;
t133 = Icges(5,1) * t194 + Icges(5,4) * t290 - Icges(5,5) * t229;
t40 = -t105 * t227 + t107 * t289 + t109 * t192 - t129 * t206 - t131 * t149 + t133 * t150;
t142 = Icges(5,5) * t188 - Icges(5,6) * t187 + Icges(5,3) * t240;
t143 = Icges(5,4) * t188 - Icges(5,2) * t187 + Icges(5,6) * t240;
t144 = Icges(5,1) * t188 - Icges(5,4) * t187 + Icges(5,5) * t240;
t182 = Icges(5,5) * t235 + Icges(5,6) * t288 + Icges(5,3) * t246;
t183 = Icges(5,4) * t235 + Icges(5,2) * t288 + Icges(5,6) * t246;
t184 = Icges(5,1) * t235 + Icges(5,4) * t288 + Icges(5,5) * t246;
t49 = -t142 * t227 + t143 * t289 + t144 * t192 - t149 * t183 + t150 * t184 - t182 * t206;
t81 = -t128 * t227 + t130 * t289 + t132 * t192;
t82 = -t129 * t227 + t131 * t289 + t133 * t192;
t93 = -t182 * t227 + t183 * t289 + t184 * t192;
t360 = -t206 * t81 - t208 * t82 - t227 * t39 - t229 * t40 + t240 * t93 + t246 * t49 + t3;
t18 = t112 * t151 + t114 * t120 + t116 * t121 + t157 * t72 + t158 * t74 - t290 * t70;
t19 = t113 * t151 + t115 * t120 + t117 * t121 + t157 * t73 + t158 * t75 - t290 * t71;
t33 = t120 * t137 + t121 * t138 + t136 * t151 + t157 * t97 + t158 * t98 - t290 * t96;
t54 = -t112 * t290 + t114 * t157 + t116 * t158;
t55 = -t113 * t290 + t115 * t157 + t117 * t158;
t69 = -t136 * t290 + t137 * t157 + t138 * t158;
t4 = -t18 * t227 - t19 * t229 - t206 * t54 - t208 * t55 + t240 * t69 + t246 * t33;
t41 = -t104 * t229 + t106 * t290 + t108 * t194 - t128 * t208 - t130 * t151 + t132 * t152;
t42 = -t105 * t229 + t107 * t290 + t109 * t194 - t129 * t208 - t131 * t151 + t133 * t152;
t50 = -t142 * t229 + t143 * t290 + t144 * t194 - t151 * t183 + t152 * t184 - t182 * t208;
t83 = -t128 * t229 + t130 * t290 + t132 * t194;
t84 = -t129 * t229 + t131 * t290 + t133 * t194;
t94 = -t182 * t229 + t183 * t290 + t184 * t194;
t359 = -t206 * t83 - t208 * t84 - t227 * t41 - t229 * t42 + t240 * t94 + t246 * t50 + t4;
t100 = t182 * t246 + t183 * t288 + t184 * t235;
t332 = t100 * t240;
t47 = t104 * t246 + t106 * t288 + t108 * t235 + t128 * t240 - t130 * t187 + t132 * t188;
t48 = t105 * t246 + t107 * t288 + t109 * t235 + t129 * t240 - t131 * t187 + t133 * t188;
t20 = t112 * t187 + t114 * t140 + t116 * t141 + t189 * t72 + t190 * t74 - t288 * t70;
t21 = t113 * t187 + t115 * t140 + t117 * t141 + t189 * t73 + t190 * t75 - t288 * t71;
t38 = t136 * t187 + t137 * t140 + t138 * t141 + t189 * t97 + t190 * t98 - t288 * t96;
t57 = -t112 * t288 + t114 * t189 + t116 * t190;
t58 = -t113 * t288 + t115 * t189 + t117 * t190;
t86 = -t136 * t288 + t137 * t189 + t138 * t190;
t6 = -t20 * t227 - t206 * t57 - t208 * t58 - t21 * t229 + t240 * t86 + t246 * t38;
t60 = t142 * t246 + t143 * t288 + t144 * t235 + t182 * t240 - t183 * t187 + t184 * t188;
t87 = t128 * t246 + t130 * t288 + t132 * t235;
t88 = t129 * t246 + t131 * t288 + t133 * t235;
t358 = -t206 * t87 - t208 * t88 - t227 * t47 - t229 * t48 + t246 * t60 + t332 + t6;
t76 = rSges(6,1) * t119 + rSges(6,2) * t118 + rSges(6,3) * t149;
t338 = pkin(4) * t150 + pkin(8) * t149 + t76;
t77 = rSges(6,1) * t121 + rSges(6,2) * t120 + rSges(6,3) * t151;
t337 = pkin(4) * t152 + pkin(8) * t151 + t77;
t122 = rSges(6,1) * t156 + rSges(6,2) * t155 - rSges(6,3) * t289;
t328 = pkin(4) * t192 - pkin(8) * t289 + t122;
t123 = rSges(6,1) * t158 + rSges(6,2) * t157 - rSges(6,3) * t290;
t327 = pkin(4) * t194 - pkin(8) * t290 + t123;
t357 = 2 * m(5);
t356 = 2 * m(6);
t355 = t149 / 0.2e1;
t354 = t151 / 0.2e1;
t353 = t187 / 0.2e1;
t352 = -t289 / 0.2e1;
t351 = -t290 / 0.2e1;
t350 = -t206 / 0.2e1;
t349 = -t208 / 0.2e1;
t348 = -t227 / 0.2e1;
t347 = -t229 / 0.2e1;
t346 = -t288 / 0.2e1;
t345 = t240 / 0.2e1;
t344 = t246 / 0.2e1;
t340 = pkin(2) * t275;
t339 = pkin(2) * t343;
t99 = rSges(6,1) * t141 + rSges(6,2) * t140 + rSges(6,3) * t187;
t336 = pkin(4) * t188 + pkin(8) * t187 + t99;
t331 = t274 * t275;
t330 = t275 * t276;
t139 = rSges(6,1) * t190 + rSges(6,2) * t189 - rSges(6,3) * t288;
t326 = pkin(4) * t235 - pkin(8) * t288 + t139;
t170 = pkin(3) * t209 - pkin(7) * t208;
t298 = t335 * t343;
t265 = pkin(2) * qJD(2) * t298 - t275 * qJD(3);
t308 = qJD(2) * t341;
t303 = pkin(2) * t308;
t237 = -t274 * t265 - t276 * t303;
t233 = t335 * t237;
t325 = t335 * t170 + t233;
t181 = pkin(3) * t230 - pkin(7) * t229;
t297 = t335 * t341;
t307 = pkin(2) * t297 - qJ(3) * t275;
t217 = -t274 * t307 + t276 * t339;
t205 = t335 * t217;
t324 = t335 * t181 + t205;
t216 = t274 * t339 + t276 * t307;
t323 = t216 * t331 + t217 * t330;
t309 = qJD(2) * t343;
t264 = qJD(3) * t335 + t309 * t340;
t322 = pkin(3) * t241 - pkin(7) * t240 - t264;
t268 = qJ(3) * t335 + t340 * t341;
t321 = -t247 * rSges(4,1) + t246 * rSges(4,2) - rSges(4,3) * t335 - t268;
t320 = -pkin(3) * t247 - pkin(7) * t246 - t268;
t236 = t276 * t265 - t274 * t303;
t319 = t236 * t331 + t237 * t330;
t316 = t335 / 0.2e1;
t313 = t343 * Icges(3,4);
t312 = t341 * Icges(3,4);
t306 = t335 * t216;
t305 = t335 * t236;
t304 = t275 * (rSges(4,1) * t241 + rSges(4,2) * t240 - t264);
t169 = pkin(3) * t207 - pkin(7) * t206;
t302 = t169 * t331 + t170 * t330 + t319;
t180 = pkin(3) * t228 - pkin(7) * t227;
t301 = t180 * t331 + t181 * t330 + t323;
t145 = rSges(5,1) * t188 - rSges(5,2) * t187 + rSges(5,3) * t240;
t300 = (-t145 + t322) * t275;
t185 = rSges(5,1) * t235 + rSges(5,2) * t288 + rSges(5,3) * t246;
t299 = (-t185 + t320) * t275;
t292 = t275 * (t322 - t336);
t291 = (t320 - t326) * t275;
t287 = -t169 * t335 - t305;
t286 = -t180 * t335 - t306;
t258 = -t274 * t341 + t276 * t298;
t260 = -t274 * t298 - t276 * t341;
t285 = t274 * t297 - t276 * t343;
t284 = -t274 * t343 - t276 * t297;
t256 = (rSges(3,1) * t343 - rSges(3,2) * t341) * t318;
t255 = (Icges(3,1) * t343 - t312) * t318;
t254 = (-Icges(3,2) * t341 + t313) * t318;
t253 = (Icges(3,5) * t343 - Icges(3,6) * t341) * t318;
t252 = t285 * qJD(2);
t251 = t260 * qJD(2);
t250 = t284 * qJD(2);
t249 = t258 * qJD(2);
t245 = t335 * rSges(3,3) + (rSges(3,1) * t341 + rSges(3,2) * t343) * t275;
t244 = Icges(3,5) * t335 + (Icges(3,1) * t341 + t313) * t275;
t243 = Icges(3,6) * t335 + (Icges(3,2) * t343 + t312) * t275;
t226 = rSges(3,1) * t251 + rSges(3,2) * t252;
t224 = rSges(3,1) * t249 + rSges(3,2) * t250;
t223 = Icges(3,1) * t251 + Icges(3,4) * t252;
t222 = Icges(3,1) * t249 + Icges(3,4) * t250;
t221 = Icges(3,4) * t251 + Icges(3,2) * t252;
t220 = Icges(3,4) * t249 + Icges(3,2) * t250;
t219 = Icges(3,5) * t251 + Icges(3,6) * t252;
t218 = Icges(3,5) * t249 + Icges(3,6) * t250;
t215 = -rSges(3,1) * t285 + rSges(3,2) * t260 + rSges(3,3) * t331;
t214 = -rSges(3,1) * t284 + rSges(3,2) * t258 - rSges(3,3) * t330;
t213 = -Icges(3,1) * t285 + Icges(3,4) * t260 + Icges(3,5) * t331;
t212 = -Icges(3,1) * t284 + Icges(3,4) * t258 - Icges(3,5) * t330;
t211 = -Icges(3,4) * t285 + Icges(3,2) * t260 + Icges(3,6) * t331;
t210 = -Icges(3,4) * t284 + Icges(3,2) * t258 - Icges(3,6) * t330;
t202 = Icges(4,1) * t247 - Icges(4,4) * t246 + Icges(4,5) * t335;
t201 = Icges(4,4) * t247 - Icges(4,2) * t246 + Icges(4,6) * t335;
t199 = -Icges(4,1) * t241 - Icges(4,4) * t240;
t198 = -Icges(4,4) * t241 - Icges(4,2) * t240;
t197 = -Icges(4,5) * t241 - Icges(4,6) * t240;
t176 = rSges(4,1) * t230 + rSges(4,2) * t229 + rSges(4,3) * t331;
t175 = rSges(4,1) * t228 + rSges(4,2) * t227 - rSges(4,3) * t330;
t174 = Icges(4,1) * t230 + Icges(4,4) * t229 + Icges(4,5) * t331;
t173 = Icges(4,1) * t228 + Icges(4,4) * t227 - Icges(4,5) * t330;
t172 = Icges(4,4) * t230 + Icges(4,2) * t229 + Icges(4,6) * t331;
t171 = Icges(4,4) * t228 + Icges(4,2) * t227 - Icges(4,6) * t330;
t168 = rSges(4,1) * t209 + rSges(4,2) * t208;
t167 = rSges(4,1) * t207 + rSges(4,2) * t206;
t166 = Icges(4,1) * t209 + Icges(4,4) * t208;
t165 = Icges(4,1) * t207 + Icges(4,4) * t206;
t164 = Icges(4,4) * t209 + Icges(4,2) * t208;
t163 = Icges(4,4) * t207 + Icges(4,2) * t206;
t162 = Icges(4,5) * t209 + Icges(4,6) * t208;
t161 = Icges(4,5) * t207 + Icges(4,6) * t206;
t160 = (t224 * t274 + t226 * t276) * t275;
t135 = rSges(5,1) * t194 + rSges(5,2) * t290 - rSges(5,3) * t229;
t134 = rSges(5,1) * t192 + rSges(5,2) * t289 - rSges(5,3) * t227;
t127 = -t167 * t335 + t276 * t304 - t305;
t126 = t168 * t335 + t274 * t304 + t233;
t111 = rSges(5,1) * t152 - rSges(5,2) * t151 - rSges(5,3) * t208;
t110 = rSges(5,1) * t150 - rSges(5,2) * t149 - rSges(5,3) * t206;
t103 = (t167 * t274 + t168 * t276) * t275 + t319;
t102 = t135 * t246 + t185 * t229;
t101 = -t134 * t246 - t185 * t227;
t95 = -t134 * t229 + t135 * t227;
t92 = -t134 * t335 + t276 * t299 + t286;
t91 = t135 * t335 + t274 * t299 + t324;
t90 = -t123 * t288 + t139 * t290;
t89 = t122 * t288 - t139 * t289;
t85 = (t134 * t274 + t135 * t276) * t275 + t301;
t80 = -t110 * t335 + t276 * t300 + t287;
t79 = t111 * t335 + t274 * t300 + t325;
t78 = -t122 * t290 + t123 * t289;
t67 = t229 * t326 + t246 * t327;
t66 = -t227 * t326 - t246 * t328;
t65 = t111 * t246 + t135 * t240 + t145 * t229 + t185 * t208;
t64 = -t110 * t246 - t134 * t240 - t145 * t227 - t185 * t206;
t63 = t276 * t291 - t328 * t335 + t286;
t62 = t274 * t291 + t327 * t335 + t324;
t61 = (t110 * t274 + t111 * t276) * t275 + t302;
t59 = t227 * t327 - t229 * t328;
t56 = -t110 * t229 + t111 * t227 - t134 * t208 + t135 * t206;
t51 = (t274 * t328 + t276 * t327) * t275 + t301;
t46 = t276 * t292 - t338 * t335 + t287;
t45 = t274 * t292 + t337 * t335 + t325;
t44 = t123 * t187 - t139 * t151 - t288 * t77 + t290 * t99;
t43 = -t122 * t187 + t139 * t149 + t288 * t76 - t289 * t99;
t37 = (t274 * t338 + t276 * t337) * t275 + t302;
t36 = t122 * t151 - t123 * t149 + t289 * t77 - t290 * t76;
t35 = t208 * t326 + t229 * t336 + t240 * t327 + t246 * t337;
t34 = -t206 * t326 - t227 * t336 - t240 * t328 - t246 * t338;
t31 = t86 * t335 + (t274 * t58 - t276 * t57) * t275;
t30 = -t227 * t57 - t229 * t58 + t246 * t86;
t29 = -t288 * t86 - t289 * t57 - t290 * t58;
t28 = t206 * t327 - t208 * t328 + t227 * t337 - t229 * t338;
t27 = t69 * t335 + (t274 * t55 - t276 * t54) * t275;
t26 = t68 * t335 + (t274 * t53 - t276 * t52) * t275;
t25 = -t227 * t54 - t229 * t55 + t246 * t69;
t24 = -t227 * t52 - t229 * t53 + t246 * t68;
t23 = -t288 * t69 - t289 * t54 - t290 * t55;
t22 = -t288 * t68 - t289 * t52 - t290 * t53;
t15 = t60 * t335 + (t274 * t48 - t276 * t47) * t275;
t14 = t50 * t335 + (t274 * t42 - t276 * t41) * t275;
t13 = t49 * t335 + (t274 * t40 - t276 * t39) * t275;
t9 = t38 * t335 + (-t20 * t276 + t21 * t274) * t275;
t8 = t33 * t335 + (-t18 * t276 + t19 * t274) * t275;
t7 = t32 * t335 + (-t16 * t276 + t17 * t274) * t275;
t5 = t149 * t57 + t151 * t58 + t187 * t86 - t20 * t289 - t21 * t290 - t288 * t38;
t2 = t149 * t54 + t151 * t55 - t18 * t289 + t187 * t69 - t19 * t290 - t288 * t33;
t1 = t149 * t52 + t151 * t53 - t16 * t289 - t17 * t290 + t187 * t68 - t288 * t32;
t10 = [0; m(3) * t160 + m(4) * t103 + m(5) * t61 + m(6) * t37; t8 * t331 - t7 * t330 + t14 * t331 - t13 * t330 + ((t211 * t252 + t213 * t251 + t219 * t331 + t221 * t260 - t223 * t285) * t331 - (t210 * t252 + t212 * t251 + t218 * t331 + t220 * t260 - t222 * t285) * t330 + (t243 * t252 + t244 * t251 + t253 * t331 + t254 * t260 - t255 * t285) * t335) * t331 - ((t211 * t250 + t213 * t249 - t219 * t330 + t221 * t258 - t223 * t284) * t331 - (t210 * t250 + t212 * t249 - t218 * t330 + t220 * t258 - t222 * t284) * t330 + (t243 * t250 + t244 * t249 - t253 * t330 + t254 * t258 - t255 * t284) * t335) * t330 + ((t162 * t331 + t164 * t229 + t166 * t230 + t172 * t208 + t174 * t209) * t331 - (t161 * t331 + t163 * t229 + t165 * t230 + t171 * t208 + t173 * t209) * t330 + (t197 * t331 + t198 * t229 + t199 * t230 + t201 * t208 + t202 * t209) * t335) * t331 - ((-t162 * t330 + t164 * t227 + t166 * t228 + t172 * t206 + t174 * t207) * t331 - (-t161 * t330 + t163 * t227 + t165 * t228 + t171 * t206 + t173 * t207) * t330 + (-t197 * t330 + t198 * t227 + t199 * t228 + t201 * t206 + t202 * t207) * t335) * t330 + t335 * t9 + (t37 * t51 + t45 * t62 + t46 * t63) * t356 + (t61 * t85 + t79 * t91 + t80 * t92) * t357 + t335 * t15 + 0.2e1 * m(4) * ((-t175 * t335 - t306) * t127 + (t176 * t335 + t205) * t126 + t323 * t103 + ((t176 * t103 + t127 * t321) * t276 + (t175 * t103 + t126 * t321) * t274) * t275) + 0.2e1 * m(3) * ((-t214 * t335 - t245 * t330) * (-t224 * t335 - t256 * t330) + (t215 * t335 - t245 * t331) * (t226 * t335 - t256 * t331) + (t214 * t274 + t215 * t276) * t275 * t160) + t335 * (((-t211 * t308 + t213 * t309 + t221 * t343 + t223 * t341) * t274 - (-t210 * t308 + t212 * t309 + t220 * t343 + t222 * t341) * t276) * t275 ^ 2 + ((-t218 * t276 + t219 * t274 - t243 * t308 + t244 * t309 + t254 * t343 + t255 * t341) * t275 + t253 * t335) * t335) + t335 * ((t335 * t197 - t246 * t198 + t247 * t199 - t240 * t201 - t241 * t202) * t335 + ((t162 * t335 - t246 * t164 + t247 * t166 - t240 * t172 - t241 * t174) * t274 - (t161 * t335 - t246 * t163 + t247 * t165 - t240 * t171 - t241 * t173) * t276) * t275); 0; m(6) * (t335 * t37 + (t274 * t46 - t276 * t45) * t275) + m(5) * (t335 * t61 + (t274 * t80 - t276 * t79) * t275) + m(4) * (t335 * t103 + (-t126 * t276 + t127 * t274) * t275); 0; m(5) * t56 + m(6) * t28; m(6) * (t28 * t51 + t34 * t63 + t35 * t62 + t37 * t59 + t45 * t67 + t46 * t66) + m(5) * (t101 * t80 + t102 * t79 + t56 * t85 + t61 * t95 + t64 * t92 + t65 * t91) + (t26 + t93 * t335 + (t274 * t82 - t276 * t81) * t275) * t350 + (t27 + t94 * t335 + (t274 * t84 - t276 * t83) * t275) * t349 + (t7 + t13) * t348 + (t8 + t14) * t347 + (t31 + t100 * t335 + (t274 * t88 - t276 * t87) * t275) * t345 + (t9 + t15) * t344 + t358 * t316 + t359 * t331 / 0.2e1 - t360 * t330 / 0.2e1; m(5) * (t56 * t335 + (t274 * t64 - t276 * t65) * t275) + m(6) * (t28 * t335 + (t274 * t34 - t276 * t35) * t275); t240 * t30 + (t332 + t358) * t246 + (-t240 * t88 - t359) * t229 + (-t240 * t87 - t360) * t227 + (t227 * t83 + t229 * t84 - t246 * t94 - t25) * t208 + (t227 * t81 + t229 * t82 - t246 * t93 - t24) * t206 + (t28 * t59 + t34 * t66 + t35 * t67) * t356 + (t101 * t64 + t102 * t65 + t56 * t95) * t357; m(6) * t36; t27 * t354 + t8 * t351 + t26 * t355 + t7 * t352 + t31 * t353 + t9 * t346 + t5 * t316 + m(6) * (t36 * t51 + t37 * t78 + t43 * t63 + t44 * t62 + t45 * t90 + t46 * t89) + (t274 * t2 / 0.2e1 - t276 * t1 / 0.2e1) * t275; m(6) * (t36 * t335 + (t274 * t43 - t276 * t44) * t275); m(6) * (t28 * t78 + t34 * t89 + t35 * t90 + t36 * t59 + t43 * t66 + t44 * t67) + t29 * t345 + t5 * t344 + t22 * t350 + t1 * t348 + t24 * t355 + t3 * t352 + t30 * t353 + t6 * t346 + t25 * t354 + t4 * t351 + t23 * t349 + t2 * t347; (t36 * t78 + t43 * t89 + t44 * t90) * t356 + t151 * t23 - t290 * t2 + t149 * t22 - t289 * t1 + t187 * t29 - t288 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
