% Calculate time derivative of joint inertia matrix for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:08
% EndTime: 2019-03-09 01:39:22
% DurationCPUTime: 8.23s
% Computational Cost: add. (21806->653), mult. (17280->935), div. (0->0), fcn. (16151->12), ass. (0->316)
t215 = qJ(1) + pkin(9);
t211 = cos(t215);
t208 = sin(t215);
t214 = pkin(10) + qJ(4);
t207 = sin(t214);
t210 = cos(t214);
t351 = Icges(5,4) * t210;
t256 = -Icges(5,2) * t207 + t351;
t140 = -Icges(5,6) * t211 + t208 * t256;
t352 = Icges(5,4) * t207;
t260 = Icges(5,1) * t210 - t352;
t143 = -Icges(5,5) * t211 + t208 * t260;
t247 = t140 * t207 - t143 * t210;
t238 = t247 * t211;
t340 = t207 * t211;
t216 = sin(pkin(11));
t331 = t211 * t216;
t218 = cos(pkin(11));
t333 = t210 * t218;
t168 = t208 * t333 - t331;
t330 = t211 * t218;
t334 = t210 * t216;
t240 = t208 * t334 + t330;
t341 = t207 * t208;
t96 = Icges(6,5) * t168 - Icges(6,6) * t240 + Icges(6,3) * t341;
t395 = t96 * t340 - t238;
t141 = Icges(5,6) * t208 + t211 * t256;
t144 = Icges(5,5) * t208 + t211 * t260;
t246 = t141 * t207 - t144 * t210;
t237 = t246 * t208;
t169 = t208 * t218 - t210 * t331;
t337 = t208 * t216;
t170 = t210 * t330 + t337;
t97 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t340;
t394 = -t97 * t341 + t237;
t370 = t208 / 0.2e1;
t393 = -t211 / 0.2e1;
t392 = -qJD(1) / 0.2e1;
t101 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t340;
t252 = Icges(5,5) * t210 - Icges(5,6) * t207;
t138 = Icges(5,3) * t208 + t211 * t252;
t99 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t340;
t391 = t101 * t168 - t138 * t211 - t240 * t99 - t394;
t100 = Icges(6,1) * t168 - Icges(6,4) * t240 + Icges(6,5) * t341;
t137 = -Icges(5,3) * t211 + t208 * t252;
t98 = Icges(6,4) * t168 - Icges(6,2) * t240 + Icges(6,6) * t341;
t390 = -t100 * t170 - t137 * t208 - t169 * t98 - t395;
t389 = t246 * t211 - t97 * t340;
t388 = t101 * t170 + t138 * t208 + t169 * t99 - t389;
t303 = t96 * t341;
t387 = -t100 * t168 + t137 * t211 + t208 * t247 + t240 * t98 - t303;
t377 = 2 * m(5);
t357 = rSges(5,1) * t210;
t278 = -rSges(5,2) * t207 + t357;
t148 = -rSges(5,3) * t211 + t208 * t278;
t335 = t210 * t211;
t201 = t208 * rSges(5,3);
t380 = -rSges(5,2) * t340 + t201;
t149 = rSges(5,1) * t335 + t380;
t182 = rSges(5,1) * t207 + rSges(5,2) * t210;
t239 = qJD(4) * t182;
t318 = qJD(1) * t208;
t293 = t207 * t318;
t316 = qJD(1) * t211;
t226 = rSges(5,2) * t293 + rSges(5,3) * t316 - t211 * t239;
t40 = (qJD(1) * t148 + t226) * t211 + (-t208 * t239 + (-t149 + t380) * qJD(1)) * t208;
t386 = t377 * t40;
t202 = pkin(5) * t218 + pkin(4);
t197 = pkin(5) * t337;
t220 = -pkin(8) - qJ(5);
t381 = -t220 * t340 + t197;
t213 = pkin(11) + qJ(6);
t206 = sin(t213);
t209 = cos(t213);
t339 = t208 * t209;
t160 = -t206 * t335 + t339;
t161 = t206 * t208 + t209 * t335;
t93 = t161 * rSges(7,1) + t160 * rSges(7,2) + rSges(7,3) * t340;
t385 = t202 * t335 + t381 + t93;
t360 = pkin(4) - t202;
t384 = t207 * t360;
t329 = qJ(5) + t220;
t383 = t210 * t329;
t221 = -pkin(7) - qJ(3);
t219 = cos(pkin(10));
t203 = pkin(3) * t219 + pkin(2);
t242 = -t203 - t278;
t363 = sin(qJ(1)) * pkin(1);
t116 = -t363 + (rSges(5,3) - t221) * t211 + t242 * t208;
t212 = cos(qJ(1)) * pkin(1);
t280 = t211 * t203 - t208 * t221 + t212;
t117 = t149 + t280;
t379 = t116 * t211 + t117 * t208;
t317 = qJD(1) * t210;
t282 = -qJD(6) + t317;
t312 = qJD(4) * t211;
t290 = t207 * t312;
t378 = t208 * t282 + t290;
t244 = rSges(4,1) * t219 - rSges(4,2) * sin(pkin(10)) + pkin(2);
t355 = rSges(4,3) + qJ(3);
t122 = t208 * t355 + t211 * t244 + t212;
t376 = 2 * m(6);
t375 = 2 * m(7);
t204 = t208 ^ 2;
t205 = t211 ^ 2;
t374 = m(6) / 0.2e1;
t373 = m(7) / 0.2e1;
t254 = Icges(6,4) * t218 - Icges(6,2) * t216;
t155 = -Icges(6,6) * t210 + t207 * t254;
t372 = t155 / 0.2e1;
t258 = Icges(6,1) * t218 - Icges(6,4) * t216;
t156 = -Icges(6,5) * t210 + t207 * t258;
t371 = t156 / 0.2e1;
t369 = -t210 / 0.2e1;
t367 = t211 / 0.2e1;
t366 = -t216 / 0.2e1;
t365 = t218 / 0.2e1;
t364 = m(5) * t182;
t362 = pkin(4) * t210;
t361 = qJD(1) / 0.2e1;
t227 = -t207 * t329 - t210 * t360;
t307 = pkin(5) * t331;
t336 = t209 * t211;
t338 = t208 * t210;
t158 = -t206 * t338 - t336;
t159 = -t206 * t211 + t209 * t338;
t274 = -rSges(7,1) * t159 - rSges(7,2) * t158;
t92 = rSges(7,3) * t341 - t274;
t359 = t208 * t227 - t307 + t92;
t194 = pkin(4) * t335;
t166 = qJ(5) * t340 + t194;
t358 = -t166 + t385;
t356 = rSges(7,3) * t207;
t354 = -rSges(6,3) - qJ(5);
t353 = -rSges(7,3) + t220;
t350 = Icges(7,4) * t206;
t349 = Icges(7,4) * t209;
t257 = Icges(7,1) * t209 - t350;
t142 = -Icges(7,5) * t210 + t207 * t257;
t343 = t142 * t209;
t342 = t207 * t202;
t332 = t210 * t220;
t112 = t170 * rSges(6,1) + t169 * rSges(6,2) + rSges(6,3) * t340;
t328 = -t112 - t166;
t273 = rSges(7,1) * t209 - rSges(7,2) * t206;
t147 = -rSges(7,3) * t210 + t207 * t273;
t327 = t147 + t383 - t384;
t271 = qJ(5) * t207 + t362;
t164 = qJD(4) * t271 - qJD(5) * t210;
t275 = rSges(6,1) * t218 - rSges(6,2) * t216;
t326 = -(rSges(6,3) * t207 + t210 * t275) * qJD(4) - t164;
t165 = t271 * t208;
t325 = t208 * t165 + t211 * t166;
t157 = -rSges(6,3) * t210 + t207 * t275;
t181 = pkin(4) * t207 - qJ(5) * t210;
t324 = -t157 - t181;
t323 = qJD(1) * t307 + t220 * t293;
t309 = qJD(5) * t207;
t188 = t211 * t309;
t199 = qJD(3) * t208;
t322 = t188 + t199;
t200 = qJD(3) * t211;
t321 = t221 * t318 + t200;
t320 = t204 + t205;
t319 = qJD(1) * t138;
t315 = qJD(4) * t207;
t314 = qJD(4) * t208;
t313 = qJD(4) * t210;
t311 = qJD(4) * t216;
t310 = qJD(4) * t218;
t308 = qJD(6) * t207;
t286 = t210 * t312;
t183 = qJ(5) * t286;
t291 = t207 * t314;
t187 = pkin(4) * t291;
t287 = t208 * t313;
t231 = t207 * t316 + t287;
t306 = t165 * t316 + t208 * (qJ(5) * t231 + qJD(1) * t194 + t208 * t309 - t187) + t211 * (-qJ(5) * t293 + t183 + t188 + (-t208 * t317 - t290) * pkin(4));
t283 = -qJD(6) * t210 + qJD(1);
t245 = t283 * t211;
t79 = t206 * t378 + t209 * t245;
t80 = t206 * t245 - t209 * t378;
t305 = t80 * rSges(7,1) + t79 * rSges(7,2) + rSges(7,3) * t286;
t89 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t340;
t91 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t340;
t269 = -t206 * t89 + t209 * t91;
t87 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t340;
t37 = t207 * t269 - t210 * t87;
t250 = Icges(7,5) * t209 - Icges(7,6) * t206;
t136 = -Icges(7,3) * t210 + t207 * t250;
t253 = -Icges(7,2) * t206 + t349;
t139 = -Icges(7,6) * t210 + t207 * t253;
t52 = t136 * t340 + t139 * t160 + t142 * t161;
t299 = t37 / 0.2e1 + t52 / 0.2e1;
t88 = Icges(7,4) * t159 + Icges(7,2) * t158 + Icges(7,6) * t341;
t90 = Icges(7,1) * t159 + Icges(7,4) * t158 + Icges(7,5) * t341;
t270 = -t206 * t88 + t209 * t90;
t86 = Icges(7,5) * t159 + Icges(7,6) * t158 + Icges(7,3) * t341;
t36 = t207 * t270 - t210 * t86;
t51 = t136 * t341 + t139 * t158 + t142 * t159;
t298 = t51 / 0.2e1 + t36 / 0.2e1;
t113 = (-rSges(7,1) * t206 - rSges(7,2) * t209) * t308 + (t210 * t273 + t356) * qJD(4);
t296 = -t227 * qJD(4) - t113 - t164;
t289 = t207 * t311;
t126 = qJD(1) * t240 + t211 * t289;
t288 = t207 * t310;
t127 = -qJD(1) * t168 - t211 * t288;
t295 = t127 * rSges(6,1) + t126 * rSges(6,2) + rSges(6,3) * t286;
t294 = -t181 - t327;
t292 = t139 * t313;
t285 = t313 / 0.2e1;
t119 = t324 * t211;
t284 = -t202 * t210 - t203;
t76 = t294 * t211;
t81 = t283 * t339 + (-t211 * t282 + t291) * t206;
t82 = t282 * t336 + (t206 * t283 - t209 * t315) * t208;
t281 = t82 * rSges(7,1) + t81 * rSges(7,2);
t128 = qJD(1) * t169 + t208 * t289;
t129 = qJD(1) * t170 - t208 * t288;
t277 = -t129 * rSges(6,1) - t128 * rSges(6,2);
t276 = -rSges(6,1) * t168 + rSges(6,2) * t240;
t272 = -t211 * t221 - t363;
t26 = t158 * t88 + t159 * t90 + t341 * t86;
t27 = t158 * t89 + t159 * t91 + t341 * t87;
t16 = t208 * t27 - t211 * t26;
t268 = t208 * t26 + t211 * t27;
t28 = t160 * t88 + t161 * t90 + t340 * t86;
t29 = t160 * t89 + t161 * t91 + t340 * t87;
t17 = t208 * t29 - t211 * t28;
t267 = t208 * t28 + t211 * t29;
t266 = t208 * t37 - t211 * t36;
t265 = t208 * t36 + t211 * t37;
t228 = t207 * t353 + t284;
t53 = -t363 + (pkin(5) * t216 - t221) * t211 + t228 * t208 + t274;
t54 = t280 + t385;
t264 = t208 * t54 + t211 * t53;
t63 = t147 * t341 + t210 * t92;
t64 = -t147 * t340 - t210 * t93;
t263 = t208 * t64 + t211 * t63;
t232 = t207 * t354 - t203 - t362;
t224 = t208 * t232 + t272;
t71 = t224 + t276;
t72 = t280 - t328;
t262 = t208 * t72 + t211 * t71;
t261 = -t208 * t93 + t211 * t92;
t259 = Icges(5,1) * t207 + t351;
t255 = Icges(5,2) * t210 + t352;
t251 = Icges(6,5) * t218 - Icges(6,6) * t216;
t236 = qJD(4) * t259;
t235 = qJD(4) * t255;
t234 = qJD(4) * (-Icges(5,5) * t207 - Icges(5,6) * t210);
t230 = t286 - t293;
t102 = (-Icges(7,5) * t206 - Icges(7,6) * t209) * t308 + (Icges(7,3) * t207 + t210 * t250) * qJD(4);
t108 = (-Icges(7,1) * t206 - t349) * t308 + (Icges(7,5) * t207 + t210 * t257) * qJD(4);
t225 = -t210 * t102 + t136 * t315 + t313 * t343 + (t108 * t207 - t139 * t308) * t209;
t121 = -t208 * t244 + t211 * t355 - t363;
t175 = t278 * qJD(4);
t171 = t181 * t318;
t135 = (Icges(6,5) * t207 + t210 * t258) * qJD(4);
t134 = (Icges(6,6) * t207 + t210 * t254) * qJD(4);
t118 = t324 * t208;
t115 = -qJD(1) * t122 + t200;
t114 = qJD(1) * t121 + t199;
t111 = rSges(6,3) * t341 - t276;
t105 = (-Icges(7,2) * t209 - t350) * t308 + (Icges(7,6) * t207 + t210 * t253) * qJD(4);
t104 = t208 * t234 + t319;
t103 = -qJD(1) * t137 + t211 * t234;
t75 = t294 * t208;
t74 = t182 * t314 + (t211 * t242 - t201 - t212) * qJD(1) + t321;
t73 = t199 + ((-t203 - t357) * t208 + t272) * qJD(1) + t226;
t70 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t231;
t69 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t230;
t68 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t231;
t67 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t230;
t58 = -t136 * t210 + (-t139 * t206 + t343) * t207;
t57 = qJD(1) * t119 + t208 * t326;
t56 = t157 * t318 + t211 * t326 + t171;
t55 = t58 * t315;
t50 = t261 * t207;
t49 = rSges(7,3) * t231 + t281;
t48 = -rSges(7,3) * t293 + t305;
t47 = Icges(7,1) * t82 + Icges(7,4) * t81 + Icges(7,5) * t231;
t46 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t230;
t45 = Icges(7,4) * t82 + Icges(7,2) * t81 + Icges(7,6) * t231;
t44 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t230;
t43 = Icges(7,5) * t82 + Icges(7,6) * t81 + Icges(7,3) * t231;
t42 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t230;
t41 = t111 * t208 + t112 * t211 + t325;
t39 = t187 + (t313 * t354 - t309) * t208 + (t211 * t232 - t212) * qJD(1) + t277 + t321;
t38 = -pkin(4) * t290 + qJD(1) * t224 + t183 + t295 + t322;
t35 = qJD(1) * t76 + t208 * t296;
t34 = t211 * t296 + t318 * t327 + t171;
t25 = (-t309 + (t210 * t353 + t342) * qJD(4)) * t208 + (t211 * t228 - t197 - t212) * qJD(1) - t281 + t321;
t24 = (-t332 - t342) * t312 + ((t284 - t356) * t208 + t272) * qJD(1) + t305 + t322 + t323;
t23 = t208 * t359 + t211 * t358 + t325;
t22 = (t147 * t314 + t49) * t210 + (-qJD(4) * t92 + t113 * t208 + t147 * t316) * t207;
t21 = (-t147 * t312 - t48) * t210 + (qJD(4) * t93 - t113 * t211 + t147 * t318) * t207;
t20 = (-t292 + (-qJD(6) * t142 - t105) * t207) * t206 + t225;
t19 = t102 * t341 + t105 * t158 + t108 * t159 + t136 * t231 + t139 * t81 + t142 * t82;
t18 = t102 * t340 + t105 * t160 + t108 * t161 + t136 * t230 + t139 * t79 + t142 * t80;
t15 = t208 * (rSges(6,3) * t287 - t277) + t211 * t295 + (t211 * t111 + t208 * t328) * qJD(1) + t306;
t14 = t261 * t313 + (-t208 * t48 + t211 * t49 + (-t208 * t92 - t211 * t93) * qJD(1)) * t207;
t13 = t207 * t267 - t210 * t52;
t12 = t207 * t268 - t210 * t51;
t11 = (qJD(4) * t269 - t42) * t210 + (qJD(4) * t87 - t206 * t44 + t209 * t46 + (-t206 * t91 - t209 * t89) * qJD(6)) * t207;
t10 = (qJD(4) * t270 - t43) * t210 + (qJD(4) * t86 - t206 * t45 + t209 * t47 + (-t206 * t90 - t209 * t88) * qJD(6)) * t207;
t9 = t87 * t287 + t158 * t44 + t159 * t46 + t81 * t89 + t82 * t91 + (t208 * t42 + t316 * t87) * t207;
t8 = t86 * t287 + t158 * t45 + t159 * t47 + t81 * t88 + t82 * t90 + (t208 * t43 + t316 * t86) * t207;
t7 = t87 * t286 + t160 * t44 + t161 * t46 + t79 * t89 + t80 * t91 + (t211 * t42 - t318 * t87) * t207;
t6 = t86 * t286 + t160 * t45 + t161 * t47 + t79 * t88 + t80 * t90 + (t211 * t43 - t318 * t86) * t207;
t5 = (-t183 + t48 + t323) * t211 + (t187 + t49) * t208 + (t205 * (-t332 + t384) + (-t342 - t383) * t204) * qJD(4) + (t359 * t211 + (-t166 - t358 + t381) * t208) * qJD(1) + t306;
t4 = qJD(1) * t268 + t208 * t9 - t211 * t8;
t3 = qJD(1) * t267 + t208 * t7 - t211 * t6;
t2 = (qJD(4) * t268 - t19) * t210 + (-qJD(1) * t16 + qJD(4) * t51 + t208 * t8 + t211 * t9) * t207;
t1 = (qJD(4) * t267 - t18) * t210 + (-qJD(1) * t17 + qJD(4) * t52 + t208 * t6 + t211 * t7) * t207;
t30 = [t225 + (t24 * t54 + t25 * t53) * t375 + (t38 * t72 + t39 * t71) * t376 + (t116 * t74 + t117 * t73) * t377 + 0.2e1 * m(4) * (t114 * t122 + t115 * t121) + (-t155 * t311 + t156 * t310) * t210 + (-t134 * t216 + t135 * t218) * t207 + (-Icges(6,3) * t210 + t207 * t251 - t255 + t260) * t315 + (-Icges(6,3) * t207 - t210 * t251 + t256 + t259) * t313 + (-t105 * t207 - t142 * t308 - t292) * t206; 0; 0; m(7) * (qJD(1) * t264 + t208 * t25 - t211 * t24) + m(6) * (qJD(1) * t262 + t208 * t39 - t211 * t38) + m(5) * (qJD(1) * t379 + t208 * t74 - t211 * t73) + m(4) * (-t114 * t211 + t115 * t208 + (t121 * t211 + t122 * t208) * qJD(1)); 0; 0; m(7) * (t24 * t75 + t25 * t76 + t34 * t53 + t35 * t54) + m(6) * (t118 * t38 + t119 * t39 + t56 * t71 + t57 * t72) + m(5) * ((-t208 * t73 - t211 * t74) * t182 - t379 * t175) + ((t141 * t392 + t235 * t370 + Icges(6,5) * t129 / 0.2e1 + Icges(6,6) * t128 / 0.2e1 + Icges(6,3) * t231 / 0.2e1) * t211 + (-Icges(6,5) * t127 / 0.2e1 - Icges(6,6) * t126 / 0.2e1 - Icges(6,3) * t230 / 0.2e1 + t140 * t392 + t235 * t393) * t208) * t210 + ((t169 * t372 + t170 * t371 - t117 * t364 + (t141 / 0.2e1 - t97 / 0.2e1) * t210 + (t144 / 0.2e1 + t101 * t365 + t99 * t366) * t207 + t299) * t211 + (-t240 * t372 + t168 * t371 + t116 * t364 + (t140 / 0.2e1 - t96 / 0.2e1) * t210 + (t143 / 0.2e1 + t100 * t365 + t98 * t366) * t207 + t298) * t208) * qJD(1) + ((t204 / 0.2e1 + t205 / 0.2e1) * t252 + t238 / 0.2e1 - t237 / 0.2e1) * qJD(4) + (t126 * t155 + t127 * t156 + t134 * t169 + t135 * t170 + t11 + t18 + (-qJD(1) * t143 - t211 * t236 - t216 * t67 + t218 * t69) * t207 + (t101 * t333 + t207 * t97 - t334 * t99) * qJD(4)) * t370 + (t128 * t155 + t129 * t156 - t134 * t240 + t135 * t168 + t19 + t10 + (qJD(1) * t144 - t208 * t236 - t216 * t68 + t218 * t70) * t207 + (t100 * t333 + t207 * t96 - t334 * t98) * qJD(4)) * t393; m(5) * t40 + m(6) * t15 + m(7) * t5; m(6) * (t208 * t56 - t211 * t57 + (t118 * t208 + t119 * t211) * qJD(1)) + m(7) * (t208 * t34 - t211 * t35 + (t208 * t75 + t211 * t76) * qJD(1)); t17 * t316 + t16 * t318 + (t23 * t5 + t34 * t76 + t35 * t75) * t375 + (t118 * t57 + t119 * t56 + t41 * t15) * t376 + t320 * t182 * t175 * t377 + (t148 * t386 + t3 + (t127 * t101 + t208 * t103 + t126 * t99 + t169 * t67 + t170 * t69 + (-t390 + t394) * qJD(1)) * t208 + t391 * t318 + t388 * t316) * t208 + (-t4 + t149 * t386 + (t129 * t100 - t211 * t104 + t128 * t98 - t240 * t68 + t168 * t70 + (-t391 + t395) * qJD(1)) * t211 + t387 * t318 + t390 * t316 + (-t129 * t101 - t128 * t99 + t240 * t67 - t168 * t69 - t127 * t100 - t126 * t98 - t169 * t68 - t170 * t70 + (t140 * t313 + t143 * t315 + t103 - (t140 * t210 + t143 * t207) * qJD(4)) * t211 + (-t104 + (-t141 * t210 - t144 * t207) * qJD(4) + t141 * t313 + t144 * t315 - t319) * t208 + (t303 + (t138 - t247) * t208 + t387 + t388 + t389) * qJD(1)) * t208) * t211; 0.2e1 * (t262 * t374 + t264 * t373) * t313 + 0.2e1 * ((t208 * t24 + t211 * t25 + t316 * t54 - t318 * t53) * t373 + (t208 * t38 + t211 * t39 + t316 * t72 - t318 * t71) * t374) * t207; (m(6) + m(7)) * t315; 0; 0.2e1 * ((t312 * t76 + t314 * t75 - t5) * t373 + (t118 * t314 + t119 * t312 - t15) * t374) * t210 + 0.2e1 * ((qJD(4) * t23 + t208 * t35 + t211 * t34 + t316 * t75 - t318 * t76) * t373 + (qJD(4) * t41 + t118 * t316 - t119 * t318 + t208 * t57 + t211 * t56) * t374) * t207; 0.4e1 * (t374 + t373) * (-0.1e1 + t320) * t207 * t313; m(7) * (t21 * t54 + t22 * t53 + t24 * t64 + t25 * t63) + t55 + (-t20 + (t208 * t298 + t211 * t299) * qJD(4)) * t210 + ((t11 / 0.2e1 + t18 / 0.2e1) * t211 + (t19 / 0.2e1 + t10 / 0.2e1) * t208 + (-t208 * t299 + t211 * t298) * qJD(1)) * t207; m(7) * t14; m(7) * (qJD(1) * t263 + t208 * t22 - t21 * t211); m(7) * (t14 * t23 + t21 * t75 + t22 * t76 + t63 * t34 + t64 * t35 + t50 * t5) + (t13 * t361 - t2 / 0.2e1 + (qJD(1) * t37 - t10) * t369 + t17 * t285) * t211 + (t1 / 0.2e1 + t12 * t361 + (qJD(1) * t36 + t11) * t369 + t16 * t285) * t208 + (qJD(4) * t266 / 0.2e1 + t4 * t370 + t3 * t367 + (t16 * t367 - t208 * t17 / 0.2e1) * qJD(1)) * t207; m(7) * ((qJD(4) * t263 - t14) * t210 + (qJD(4) * t50 + t208 * t21 + t211 * t22 + (-t208 * t63 + t211 * t64) * qJD(1)) * t207); (t14 * t50 + t21 * t64 + t22 * t63) * t375 + (t20 * t210 - t55 + (t208 * t12 + t211 * t13 - t210 * t265) * qJD(4)) * t210 + (t211 * t1 + t208 * t2 - t210 * (t10 * t208 + t11 * t211) + (t207 * t265 - t58 * t210) * qJD(4) + (t211 * t12 - t208 * t13 + t210 * t266) * qJD(1)) * t207;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t30(1) t30(2) t30(4) t30(7) t30(11) t30(16); t30(2) t30(3) t30(5) t30(8) t30(12) t30(17); t30(4) t30(5) t30(6) t30(9) t30(13) t30(18); t30(7) t30(8) t30(9) t30(10) t30(14) t30(19); t30(11) t30(12) t30(13) t30(14) t30(15) t30(20); t30(16) t30(17) t30(18) t30(19) t30(20) t30(21);];
Mq  = res;
