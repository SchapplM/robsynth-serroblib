% Calculate time derivative of joint inertia matrix for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:32
% EndTime: 2019-03-09 01:48:47
% DurationCPUTime: 9.99s
% Computational Cost: add. (10628->699), mult. (17714->991), div. (0->0), fcn. (16449->8), ass. (0->330)
t224 = sin(qJ(1));
t385 = -t224 / 0.2e1;
t226 = cos(qJ(1));
t409 = t226 / 0.2e1;
t375 = -qJD(1) / 0.2e1;
t225 = cos(qJ(4));
t223 = sin(qJ(4));
t321 = qJD(4) * t226;
t296 = t223 * t321;
t326 = qJD(1) * t224;
t231 = t225 * t326 + t296;
t217 = pkin(9) + qJ(6);
t208 = sin(t217);
t209 = cos(t217);
t359 = Icges(7,4) * t209;
t256 = -Icges(7,2) * t208 + t359;
t137 = Icges(7,6) * t223 + t225 * t256;
t360 = Icges(7,4) * t208;
t259 = Icges(7,1) * t209 - t360;
t138 = Icges(7,5) * t223 + t225 * t259;
t408 = -t137 * t208 + t138 * t209;
t221 = cos(pkin(9));
t205 = pkin(5) * t221 + pkin(4);
t297 = t225 * t321;
t292 = -qJD(6) * t223 - qJD(1);
t245 = t292 * t226;
t291 = qJD(1) * t223 + qJD(6);
t399 = t224 * t291 - t297;
t79 = t208 * t399 + t209 * t245;
t80 = t208 * t245 - t209 * t399;
t46 = t80 * rSges(7,1) + t79 * rSges(7,2) + t231 * rSges(7,3);
t407 = t205 * t297 + t46;
t222 = -pkin(8) - qJ(5);
t339 = t225 * t226;
t345 = t223 * t226;
t343 = t224 * t209;
t153 = -t208 * t345 - t343;
t344 = t224 * t208;
t154 = t209 * t345 - t344;
t91 = t154 * rSges(7,1) + t153 * rSges(7,2) - rSges(7,3) * t339;
t406 = -t205 * t345 - t222 * t339 - t91;
t373 = pkin(4) - t205;
t405 = t223 * t373;
t362 = Icges(5,4) * t223;
t258 = Icges(5,2) * t225 + t362;
t158 = -Icges(5,6) * t224 + t226 * t258;
t361 = Icges(5,4) * t225;
t261 = Icges(5,1) * t223 + t361;
t160 = -Icges(5,5) * t224 + t226 * t261;
t249 = t158 * t225 + t160 * t223;
t238 = t249 * t224;
t218 = t224 ^ 2;
t219 = t226 ^ 2;
t328 = t218 + t219;
t404 = -t137 * t209 - t138 * t208;
t215 = t226 * qJ(2);
t285 = rSges(5,1) * t223 + rSges(5,2) * t225;
t374 = -pkin(1) - qJ(3);
t240 = -t285 + t374;
t381 = -rSges(5,3) - pkin(7);
t227 = t224 * t240 + t226 * t381;
t116 = t215 + t227;
t213 = t224 * qJ(2);
t329 = t226 * pkin(1) + t213;
t303 = t226 * qJ(3) + t329;
t332 = rSges(5,1) * t345 + rSges(5,2) * t339;
t117 = t224 * t381 + t303 + t332;
t403 = t116 * t226 + t117 * t224;
t196 = rSges(5,1) * t297;
t325 = qJD(1) * t226;
t331 = qJ(2) * t325 + qJD(2) * t224;
t304 = qJD(3) * t226 + t331;
t72 = -rSges(5,2) * t296 + qJD(1) * t227 + t196 + t304;
t371 = rSges(5,2) * t223;
t189 = rSges(5,1) * t225 - t371;
t212 = qJD(2) * t226;
t330 = pkin(7) * t326 + t212;
t73 = (-t189 * qJD(4) - qJD(3)) * t224 + ((rSges(5,3) - qJ(2)) * t224 + t240 * t226) * qJD(1) + t330;
t402 = t224 * t72 + t226 * t73;
t255 = Icges(5,5) * t223 + Icges(5,6) * t225;
t155 = Icges(5,3) * t226 + t224 * t255;
t401 = qJD(1) * t155;
t157 = Icges(5,6) * t226 + t224 * t258;
t159 = Icges(5,5) * t226 + t224 * t261;
t322 = qJD(4) * t225;
t298 = t224 * t322;
t398 = t226 * t291 + t298;
t280 = rSges(7,1) * t209 - rSges(7,2) * t208;
t140 = rSges(7,3) * t223 + t225 * t280;
t319 = qJD(6) * t225;
t95 = (-rSges(7,1) * t208 - rSges(7,2) * t209) * t319 + (rSges(7,3) * t225 - t223 * t280) * qJD(4);
t21 = (-t140 * t321 + t46) * t223 + (qJD(4) * t91 - t140 * t326 + t226 * t95) * t225;
t323 = qJD(4) * t224;
t299 = t223 * t323;
t300 = t225 * t325;
t232 = t299 - t300;
t246 = t292 * t224;
t81 = -t208 * t398 + t209 * t246;
t82 = t208 * t246 + t209 * t398;
t286 = t82 * rSges(7,1) + t81 * rSges(7,2);
t47 = rSges(7,3) * t232 + t286;
t151 = t209 * t226 - t223 * t344;
t152 = t208 * t226 + t223 * t343;
t281 = -t152 * rSges(7,1) - t151 * rSges(7,2);
t340 = t224 * t225;
t90 = -rSges(7,3) * t340 - t281;
t22 = (t140 * t323 - t47) * t223 + (-qJD(4) * t90 - t140 * t325 - t224 * t95) * t225;
t58 = -t140 * t340 - t223 * t90;
t59 = t140 * t339 + t223 * t91;
t397 = -qJD(1) * (t224 * t58 - t226 * t59) + t21 * t224 + t22 * t226;
t338 = qJ(5) + t222;
t396 = t223 * t338 + t225 * t373;
t395 = 2 * m(5);
t394 = 2 * m(6);
t393 = 2 * m(7);
t392 = m(6) / 0.2e1;
t391 = m(7) / 0.2e1;
t220 = sin(pkin(9));
t257 = Icges(6,4) * t221 - Icges(6,2) * t220;
t146 = Icges(6,6) * t223 + t225 * t257;
t390 = -t146 / 0.2e1;
t260 = Icges(6,1) * t221 - Icges(6,4) * t220;
t147 = Icges(6,5) * t223 + t225 * t260;
t389 = -t147 / 0.2e1;
t388 = t220 / 0.2e1;
t387 = -t221 / 0.2e1;
t386 = t223 / 0.2e1;
t384 = -t226 / 0.2e1;
t382 = rSges(3,2) - pkin(1);
t380 = m(5) * t189;
t379 = pkin(4) * t223;
t378 = pkin(4) * t225;
t377 = pkin(5) * t220;
t376 = pkin(7) * t226;
t253 = Icges(7,5) * t209 - Icges(7,6) * t208;
t136 = Icges(7,3) * t223 + t225 * t253;
t324 = qJD(4) * t223;
t92 = (-Icges(7,5) * t208 - Icges(7,6) * t209) * t319 + (Icges(7,3) * t225 - t223 * t253) * qJD(4);
t94 = (-Icges(7,1) * t208 - t359) * t319 + (Icges(7,5) * t225 - t223 * t259) * qJD(4);
t228 = t225 * t209 * t94 + t136 * t322 + t223 * t92 - t408 * t324;
t93 = (-Icges(7,2) * t209 - t360) * t319 + (Icges(7,6) * t225 - t223 * t256) * qJD(4);
t369 = t208 * t93;
t55 = t136 * t223 + t408 * t225;
t372 = ((qJD(6) * t404 - t369) * t225 + t228) * t223 + t55 * t322;
t370 = rSges(6,3) * t225;
t366 = -rSges(6,3) - qJ(5);
t365 = rSges(7,3) - t222;
t204 = pkin(4) * t345;
t175 = -qJ(5) * t339 + t204;
t342 = t224 * t220;
t316 = pkin(5) * t342;
t364 = t175 + t316 + t406;
t363 = (-t225 * t338 + t405) * qJD(4) + t95;
t349 = t220 * t223;
t348 = t220 * t226;
t347 = t221 * t223;
t346 = t222 * t225;
t341 = t224 * t221;
t201 = qJ(5) * t340;
t315 = t224 * t379;
t174 = -t201 + t315;
t170 = t223 * t341 + t348;
t242 = -t221 * t226 + t223 * t342;
t283 = -t170 * rSges(6,1) + rSges(6,2) * t242;
t314 = rSges(6,3) * t340;
t337 = t283 + t314 - t174;
t336 = t140 - t396;
t168 = qJD(5) * t223 + (qJ(5) * t225 - t379) * qJD(4);
t188 = qJ(5) * t223 + t378;
t335 = t224 * t168 + t188 * t325;
t171 = -t220 * t345 - t341;
t172 = t221 * t345 - t342;
t334 = t172 * rSges(6,1) + t171 * rSges(6,2);
t156 = -Icges(5,3) * t224 + t226 * t255;
t327 = qJD(1) * t156;
t320 = qJD(5) * t225;
t318 = -rSges(4,3) + t374;
t317 = pkin(5) * t348;
t96 = Icges(6,5) * t170 - Icges(6,6) * t242 - Icges(6,3) * t340;
t313 = t96 * t340;
t97 = Icges(6,5) * t172 + Icges(6,6) * t171 - Icges(6,3) * t339;
t312 = t97 * t340;
t311 = t96 * t339;
t310 = t97 * t339;
t86 = Icges(7,4) * t152 + Icges(7,2) * t151 - Icges(7,6) * t340;
t88 = Icges(7,1) * t152 + Icges(7,4) * t151 - Icges(7,5) * t340;
t279 = t208 * t86 - t209 * t88;
t84 = Icges(7,5) * t152 + Icges(7,6) * t151 - Icges(7,3) * t340;
t32 = t223 * t84 - t225 * t279;
t49 = -t136 * t340 + t137 * t151 + t138 * t152;
t309 = t32 / 0.2e1 + t49 / 0.2e1;
t87 = Icges(7,4) * t154 + Icges(7,2) * t153 - Icges(7,6) * t339;
t89 = Icges(7,1) * t154 + Icges(7,4) * t153 - Icges(7,5) * t339;
t278 = t208 * t87 - t209 * t89;
t85 = Icges(7,5) * t154 + Icges(7,6) * t153 - Icges(7,3) * t339;
t33 = t223 * t85 - t225 * t278;
t50 = -t136 * t339 + t153 * t137 + t154 * t138;
t308 = t33 / 0.2e1 + t50 / 0.2e1;
t307 = -t317 - t201 - t224 * (t346 - t405) - t174 - t90;
t306 = pkin(4) * t297 + t231 * qJ(5);
t200 = t224 * t320;
t305 = t200 + t330;
t302 = -pkin(7) - t377;
t295 = t324 / 0.2e1;
t294 = t336 * t226;
t183 = t285 * qJD(4);
t293 = t328 * t183;
t124 = qJD(1) * t242 - t220 * t297;
t125 = -qJD(1) * t170 + t221 * t297;
t290 = -t125 * rSges(6,1) - t124 * rSges(6,2) - t231 * rSges(6,3);
t289 = t374 - t379;
t288 = -t205 * t223 + t374;
t287 = t302 * t226;
t126 = qJD(1) * t171 - t220 * t298;
t127 = qJD(1) * t172 + t221 * t298;
t284 = -t127 * rSges(6,1) - t126 * rSges(6,2);
t282 = rSges(6,1) * t221 - rSges(6,2) * t220;
t24 = (-t222 * t324 - t320) * t226 + (t287 + (t288 - t346) * t224) * qJD(1) + t304 + t407;
t229 = t225 * t365 + t288;
t25 = (-qJD(3) + (-t205 * t225 - t223 * t365) * qJD(4)) * t224 + ((-qJ(2) + t377) * t224 + t229 * t226) * qJD(1) - t286 + t305;
t276 = t224 * t24 + t226 * t25;
t26 = t151 * t86 + t152 * t88 - t340 * t84;
t27 = t151 * t87 + t152 * t89 - t340 * t85;
t275 = t27 * t224 - t226 * t26;
t274 = t224 * t26 + t226 * t27;
t28 = t153 * t86 + t154 * t88 - t339 * t84;
t29 = t153 * t87 + t154 * t89 - t339 * t85;
t273 = t29 * t224 - t226 * t28;
t272 = t224 * t28 + t226 * t29;
t149 = t226 * t168;
t30 = t149 + t363 * t226 + (-t188 - t336) * t326;
t31 = qJD(1) * t294 + t224 * t363 + t335;
t271 = t31 * t224 + t226 * t30;
t270 = t33 * t224 - t226 * t32;
t269 = t32 * t224 + t226 * t33;
t237 = -t226 * t320 + t306;
t38 = (t224 * t289 - t376) * qJD(1) + t237 - t290 + t304;
t193 = qJ(5) * t300;
t241 = t289 + t370;
t39 = t193 + (-qJD(3) + (t223 * t366 - t378) * qJD(4)) * t224 + (t226 * t241 - t213) * qJD(1) + t284 + t305;
t268 = t224 * t38 + t226 * t39;
t53 = t224 * t229 + t215 + t281 + t287;
t54 = t224 * t302 + t303 - t406;
t267 = t224 * t54 + t226 * t53;
t139 = (-t223 * t282 + t370) * qJD(4);
t150 = rSges(6,3) * t223 + t225 * t282;
t56 = t139 * t226 + t149 + (-t150 - t188) * t326;
t57 = t224 * t139 + t150 * t325 + t335;
t266 = t57 * t224 + t226 * t56;
t265 = t224 * t59 + t226 * t58;
t70 = t224 * t241 + t201 + t215 + t283 - t376;
t71 = -t224 * pkin(7) + t339 * t366 + t204 + t303 + t334;
t263 = t224 * t71 + t226 * t70;
t262 = t224 * t91 - t226 * t90;
t254 = Icges(6,5) * t221 - Icges(6,6) * t220;
t250 = t157 * t225 + t159 * t223;
t247 = (t392 + t391) * t324;
t243 = rSges(3,3) * t226 + t224 * t382;
t239 = t250 * t226;
t236 = qJD(4) * (Icges(5,1) * t225 - t362);
t235 = qJD(4) * (-Icges(5,2) * t223 + t361);
t234 = qJD(4) * (Icges(5,5) * t225 - Icges(5,6) * t223);
t233 = rSges(4,2) * t226 + t224 * t318;
t179 = t226 * t188;
t178 = t224 * t188;
t165 = -rSges(3,2) * t226 + t224 * rSges(3,3) + t329;
t164 = t215 + t243;
t163 = -rSges(5,3) * t224 + t332;
t162 = rSges(5,3) * t226 + t224 * t285;
t161 = t175 * t326;
t142 = t224 * rSges(4,2) + rSges(4,3) * t226 + t303;
t141 = t215 + t233;
t135 = (Icges(6,5) * t225 - t223 * t260) * qJD(4);
t134 = (Icges(6,6) * t225 - t223 * t257) * qJD(4);
t130 = t212 + (t382 * t226 + (-rSges(3,3) - qJ(2)) * t224) * qJD(1);
t129 = qJD(1) * t243 + t331;
t119 = t150 * t226 + t179;
t118 = t150 * t224 + t178;
t115 = -qJD(3) * t224 + t212 + ((-rSges(4,2) - qJ(2)) * t224 + t318 * t226) * qJD(1);
t114 = qJD(1) * t233 + t304;
t109 = t224 * t234 + t327;
t108 = t226 * t234 - t401;
t107 = qJ(5) * t299 - t193 - t200 + (t223 * t325 + t298) * pkin(4);
t106 = -qJD(1) * t315 + t237;
t103 = -rSges(6,3) * t339 + t334;
t101 = Icges(6,1) * t172 + Icges(6,4) * t171 - Icges(6,5) * t339;
t100 = Icges(6,1) * t170 - Icges(6,4) * t242 - Icges(6,5) * t340;
t99 = Icges(6,4) * t172 + Icges(6,2) * t171 - Icges(6,6) * t339;
t98 = Icges(6,4) * t170 - Icges(6,2) * t242 - Icges(6,6) * t340;
t75 = t179 + t294;
t74 = t224 * t336 + t178;
t69 = -t224 * t156 + t249 * t226;
t68 = -t224 * t155 + t239;
t67 = t156 * t226 + t238;
t66 = t155 * t226 + t224 * t250;
t65 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t232;
t64 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t231;
t63 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t232;
t62 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t231;
t51 = t262 * t225;
t48 = (-t103 - t175) * t226 + t337 * t224;
t45 = Icges(7,1) * t82 + Icges(7,4) * t81 + Icges(7,5) * t232;
t44 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t231;
t43 = Icges(7,4) * t82 + Icges(7,2) * t81 + Icges(7,6) * t232;
t42 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t231;
t41 = Icges(7,5) * t82 + Icges(7,6) * t81 + Icges(7,3) * t232;
t40 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t231;
t37 = t172 * t101 + t171 * t99 - t310;
t36 = t172 * t100 + t171 * t98 - t311;
t35 = t101 * t170 - t242 * t99 - t312;
t34 = t100 * t170 - t242 * t98 - t313;
t23 = (-t175 + t364) * t226 + t307 * t224;
t19 = t161 + (-rSges(6,3) * t299 + qJD(1) * t103 - t107 + t284) * t224 + (-t106 + (t314 + t337) * qJD(1) + t290) * t226;
t16 = t136 * t232 + t81 * t137 + t82 * t138 + t151 * t93 + t152 * t94 - t340 * t92;
t15 = t136 * t231 + t79 * t137 + t80 * t138 + t153 * t93 + t154 * t94 - t339 * t92;
t14 = -t262 * t324 + (t224 * t46 - t226 * t47 + (t224 * t90 + t226 * t91) * qJD(1)) * t225;
t13 = t50 * t223 - t225 * t272;
t12 = t49 * t223 - t225 * t274;
t11 = (qJD(4) * t278 + t40) * t223 + (qJD(4) * t85 - t208 * t42 + t209 * t44 + (-t208 * t89 - t209 * t87) * qJD(6)) * t225;
t10 = (qJD(4) * t279 + t41) * t223 + (qJD(4) * t84 - t208 * t43 + t209 * t45 + (-t208 * t88 - t209 * t86) * qJD(6)) * t225;
t9 = -t85 * t300 + t151 * t42 + t152 * t44 + t81 * t87 + t82 * t89 + (-t225 * t40 + t324 * t85) * t224;
t8 = -t84 * t300 + t151 * t43 + t152 * t45 + t81 * t86 + t82 * t88 + (-t225 * t41 + t324 * t84) * t224;
t7 = t85 * t296 + t153 * t42 + t154 * t44 + t79 * t87 + t80 * t89 + (-t226 * t40 + t326 * t85) * t225;
t6 = t84 * t296 + t153 * t43 + t154 * t45 + t79 * t86 + t80 * t88 + (-t226 * t41 + t326 * t84) * t225;
t5 = t161 + (t222 * t296 - t106 + t306 - t407) * t226 + (t323 * t396 - t107 - t193 - t47) * t224 + ((t316 - t364) * t224 + (t307 + t317) * t226) * qJD(1);
t4 = -qJD(1) * t274 - t9 * t224 + t226 * t8;
t3 = -qJD(1) * t272 - t7 * t224 + t226 * t6;
t2 = (qJD(4) * t274 + t16) * t223 + (qJD(1) * t275 + qJD(4) * t49 - t224 * t8 - t226 * t9) * t225;
t1 = (qJD(4) * t272 + t15) * t223 + (qJD(1) * t273 + qJD(4) * t50 - t224 * t6 - t226 * t7) * t225;
t17 = [-t261 * t322 + t228 + (t24 * t54 + t25 * t53) * t393 + (t38 * t71 + t39 * t70) * t394 + (t116 * t73 + t117 * t72) * t395 + 0.2e1 * m(4) * (t114 * t142 + t115 * t141) + 0.2e1 * m(3) * (t129 * t165 + t130 * t164) + t404 * t319 + (Icges(6,3) * t322 - t236) * t223 + (t146 * t220 - t147 * t221 - t223 * t254 + t258) * t324 + (Icges(6,3) * t324 - t134 * t220 + t135 * t221 + t254 * t322 - t235 - t369) * t225; m(7) * (qJD(1) * t267 + t224 * t25 - t226 * t24) + m(6) * (qJD(1) * t263 + t224 * t39 - t226 * t38) + m(5) * (qJD(1) * t403 + t224 * t73 - t226 * t72) + m(4) * (-t114 * t226 + t224 * t115 + (t141 * t226 + t142 * t224) * qJD(1)) + m(3) * (-t129 * t226 + t224 * t130 + (t164 * t226 + t165 * t224) * qJD(1)); 0; m(7) * ((-t224 * t53 + t226 * t54) * qJD(1) + t276) + m(6) * ((-t224 * t70 + t226 * t71) * qJD(1) + t268) + m(5) * ((-t116 * t224 + t117 * t226) * qJD(1) + t402) + m(4) * (t224 * t114 + t115 * t226 + (-t141 * t224 + t142 * t226) * qJD(1)); 0; 0; m(5) * (-t183 * t403 + t189 * t402) + m(6) * (t118 * t38 + t119 * t39 + t56 * t70 + t57 * t71) + m(7) * (t24 * t74 + t25 * t75 + t30 * t53 + t31 * t54) + ((Icges(6,5) * t127 / 0.2e1 + Icges(6,6) * t126 / 0.2e1 + Icges(6,3) * t232 / 0.2e1 + t158 * t375 + t235 * t385) * t226 + (-Icges(6,5) * t125 / 0.2e1 - Icges(6,6) * t124 / 0.2e1 - Icges(6,3) * t231 / 0.2e1 + t157 * t375 + t235 * t409) * t224) * t223 + ((t171 * t390 + t172 * t389 + t117 * t380 + (-t97 / 0.2e1 + t158 / 0.2e1) * t223 + (t101 * t387 + t99 * t388 - t160 / 0.2e1) * t225 - t308) * t226 + (-t116 * t380 - t242 * t390 + t170 * t389 + (t157 / 0.2e1 - t96 / 0.2e1) * t223 + (-t159 / 0.2e1 + t100 * t387 + t98 * t388) * t225 - t309) * t224) * qJD(1) + (-(t218 / 0.2e1 + t219 / 0.2e1) * t255 + t250 * t384 + t238 / 0.2e1) * qJD(4) + (t124 * t146 + t125 * t147 + t171 * t134 + t172 * t135 + t15 + t11 + (-qJD(1) * t159 - t220 * t62 + t221 * t64 + t226 * t236) * t225 + (-t101 * t347 + t225 * t97 + t349 * t99) * qJD(4)) * t385 + (t126 * t146 + t127 * t147 - t242 * t134 + t170 * t135 + t16 + t10 + (qJD(1) * t160 - t220 * t63 + t221 * t65 + t224 * t236) * t225 + (-t100 * t347 + t225 * t96 + t349 * t98) * qJD(4)) * t409; m(6) * (t56 * t224 - t226 * t57 + (t118 * t224 + t119 * t226) * qJD(1)) + m(7) * (t30 * t224 - t226 * t31 + (t224 * t74 + t226 * t75) * qJD(1)); m(6) * ((t118 * t226 - t119 * t224) * qJD(1) + t266) + m(7) * ((-t224 * t75 + t226 * t74) * qJD(1) + t271) - m(5) * t293; -t224 * t3 + t226 * t4 + (t23 * t5 + t30 * t75 + t31 * t74) * t393 + (t118 * t57 + t119 * t56 + t19 * t48) * t394 + t226 * ((t127 * t100 + t126 * t98 - t242 * t63 + t170 * t65 + (-t35 - t311) * qJD(1)) * t226 + (-t127 * t101 - t126 * t99 + t242 * t62 - t170 * t64 + (-t34 + t310) * qJD(1)) * t224) + (-t189 * t293 + (-t226 * t196 + (-t189 * t218 + t219 * t371) * qJD(4) + (rSges(5,3) * t328 - t226 * t162 + t224 * t163) * qJD(1)) * (-t224 * t162 - t163 * t226)) * t395 + t226 * ((t226 * t109 + (-t67 + t239) * qJD(1)) * t226 + (-t66 * qJD(1) + (t158 * t324 - t160 * t322 + t327) * t224 + (-t108 + (-t157 * t223 + t159 * t225) * qJD(4) + (-t155 - t249) * qJD(1)) * t226) * t224) - t224 * ((-t125 * t101 - t124 * t99 - t171 * t62 - t172 * t64 + (-t36 - t312) * qJD(1)) * t224 + (t125 * t100 + t124 * t98 + t171 * t63 + t172 * t65 + (-t37 + t313) * qJD(1)) * t226) - t224 * ((t224 * t108 + (-t68 + t238) * qJD(1)) * t224 + (-t69 * qJD(1) + (-t157 * t324 + t159 * t322 - t401) * t226 + (-t109 + (t158 * t223 - t160 * t225) * qJD(4) + (t156 - t250) * qJD(1)) * t224) * t226) + (t275 + (-t34 - t66) * t226 + (t35 + t67) * t224) * t326 + (t273 + (-t36 - t68) * t226 + (t37 + t69) * t224) * t325; 0.2e1 * (t263 * t392 + t267 * t391) * t324 + 0.2e1 * ((-t325 * t54 + t326 * t53 - t276) * t391 + (-t325 * t71 + t326 * t70 - t268) * t392) * t225; 0; 0.2e1 * t328 * t247; 0.2e1 * ((t321 * t75 + t323 * t74 + t5) * t391 + (t118 * t323 + t119 * t321 + t19) * t392) * t223 + 0.2e1 * ((qJD(4) * t23 - t325 * t74 + t326 * t75 - t271) * t391 + (qJD(4) * t48 - t118 * t325 + t119 * t326 - t266) * t392) * t225; 0.4e1 * (0.1e1 - t328) * t225 * t247; m(7) * (t21 * t54 + t22 * t53 + t24 * t59 + t25 * t58) + (t224 * t309 + t226 * t308) * t324 + ((-t11 / 0.2e1 - t15 / 0.2e1) * t226 + (-t10 / 0.2e1 - t16 / 0.2e1) * t224 + (t224 * t308 - t226 * t309) * qJD(1)) * t225 + t372; m(7) * (qJD(1) * t265 - t21 * t226 + t22 * t224); m(7) * t397; m(7) * (t14 * t23 + t21 * t74 + t22 * t75 + t30 * t58 + t31 * t59 + t5 * t51) + (t13 * t375 + t2 / 0.2e1 - t273 * t295 + (-qJD(1) * t33 + t10) * t386) * t226 + (-t1 / 0.2e1 + t12 * t375 - t275 * t295 + (-qJD(1) * t32 - t11) * t386) * t224 + (t4 * t385 + t3 * t384 - qJD(4) * t270 / 0.2e1 + (t273 * t385 - t275 * t384) * qJD(1)) * t225; m(7) * ((qJD(4) * t265 + t14) * t223 + (qJD(4) * t51 - t397) * t225); (t14 * t51 + t21 * t59 + t22 * t58) * t393 + ((t224 * t12 + t226 * t13 + t223 * t269) * qJD(4) + t372) * t223 + (-t226 * t1 - t224 * t2 + t223 * (-t10 * t224 - t11 * t226) + (t55 * t223 - t225 * t269) * qJD(4) + (-t226 * t12 + t224 * t13 + t223 * t270) * qJD(1)) * t225;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
