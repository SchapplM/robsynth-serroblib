% Calculate time derivative of joint inertia matrix for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:09
% EndTime: 2019-12-05 16:30:49
% DurationCPUTime: 16.58s
% Computational Cost: add. (41699->884), mult. (98998->1254), div. (0->0), fcn. (112330->12), ass. (0->362)
t315 = sin(pkin(9));
t318 = cos(pkin(9));
t323 = cos(qJ(2));
t319 = cos(pkin(5));
t322 = sin(qJ(2));
t374 = t319 * t322;
t299 = t315 * t323 + t318 * t374;
t316 = sin(pkin(5));
t321 = sin(qJ(3));
t376 = t316 * t321;
t386 = cos(qJ(3));
t281 = t299 * t386 - t318 * t376;
t314 = sin(pkin(10));
t317 = cos(pkin(10));
t373 = t319 * t323;
t326 = -t315 * t322 + t318 * t373;
t240 = -t281 * t314 - t317 * t326;
t380 = t326 * t314;
t241 = t281 * t317 - t380;
t339 = t316 * t386;
t324 = -t299 * t321 - t318 * t339;
t171 = Icges(5,5) * t241 + Icges(5,6) * t240 - Icges(5,3) * t324;
t204 = Icges(4,4) * t281 + Icges(4,2) * t324 - Icges(4,6) * t326;
t411 = -t171 + t204;
t348 = t315 * t374;
t301 = t318 * t323 - t348;
t283 = t301 * t386 + t315 * t376;
t300 = t315 * t373 + t318 * t322;
t242 = -t283 * t314 + t300 * t317;
t379 = t300 * t314;
t243 = t283 * t317 + t379;
t325 = -t301 * t321 + t315 * t339;
t172 = Icges(5,5) * t243 + Icges(5,6) * t242 - Icges(5,3) * t325;
t205 = Icges(4,4) * t283 + Icges(4,2) * t325 + Icges(4,6) * t300;
t410 = -t172 + t205;
t173 = Icges(5,4) * t241 + Icges(5,2) * t240 - Icges(5,6) * t324;
t175 = Icges(5,1) * t241 + Icges(5,4) * t240 - Icges(5,5) * t324;
t202 = Icges(4,5) * t281 + Icges(4,6) * t324 - Icges(4,3) * t326;
t206 = Icges(4,1) * t281 + Icges(4,4) * t324 - Icges(4,5) * t326;
t407 = t173 * t240 + t175 * t241 - t202 * t326 + t206 * t281 + t324 * t411;
t174 = Icges(5,4) * t243 + Icges(5,2) * t242 - Icges(5,6) * t325;
t176 = Icges(5,1) * t243 + Icges(5,4) * t242 - Icges(5,5) * t325;
t203 = Icges(4,5) * t283 + Icges(4,6) * t325 + Icges(4,3) * t300;
t207 = Icges(4,1) * t283 + Icges(4,4) * t325 + Icges(4,5) * t300;
t406 = t174 * t240 + t176 * t241 - t203 * t326 + t207 * t281 + t324 * t410;
t405 = t173 * t242 + t175 * t243 + t202 * t300 + t206 * t283 + t325 * t411;
t404 = t174 * t242 + t176 * t243 + t203 * t300 + t207 * t283 + t325 * t410;
t302 = -t319 * t386 + t322 * t376;
t303 = t319 * t321 + t322 * t339;
t375 = t316 * t323;
t256 = Icges(4,5) * t303 - Icges(4,6) * t302 - Icges(4,3) * t375;
t257 = Icges(4,4) * t303 - Icges(4,2) * t302 - Icges(4,6) * t375;
t258 = Icges(4,1) * t303 - Icges(4,4) * t302 - Icges(4,5) * t375;
t138 = -t256 * t326 + t257 * t324 + t258 * t281;
t278 = -t303 * t314 - t317 * t375;
t349 = t314 * t375;
t279 = t303 * t317 - t349;
t212 = Icges(5,5) * t279 + Icges(5,6) * t278 + Icges(5,3) * t302;
t213 = Icges(5,4) * t279 + Icges(5,2) * t278 + Icges(5,6) * t302;
t214 = Icges(5,1) * t279 + Icges(5,4) * t278 + Icges(5,5) * t302;
t98 = -t212 * t324 + t213 * t240 + t214 * t241;
t409 = t138 + t98;
t139 = t256 * t300 + t257 * t325 + t258 * t283;
t99 = -t212 * t325 + t213 * t242 + t214 * t243;
t408 = t139 + t99;
t124 = -t202 * t375 - t302 * t204 + t303 * t206;
t89 = t171 * t302 + t173 * t278 + t175 * t279;
t403 = t124 + t89;
t125 = -t203 * t375 - t302 * t205 + t303 * t207;
t90 = t172 * t302 + t174 * t278 + t176 * t279;
t402 = t125 + t90;
t352 = qJD(2) * t322;
t338 = t316 * t352;
t290 = t299 * qJD(2);
t351 = qJD(2) * t323;
t292 = -qJD(2) * t348 + t318 * t351;
t313 = pkin(10) + qJ(5);
t311 = sin(t313);
t312 = cos(t313);
t232 = t281 * t312 - t311 * t326;
t289 = t326 * qJD(2);
t237 = qJD(3) * t324 + t289 * t386;
t165 = -qJD(5) * t232 - t237 * t311 + t290 * t312;
t231 = -t281 * t311 - t312 * t326;
t166 = qJD(5) * t231 + t237 * t312 + t290 * t311;
t236 = t281 * qJD(3) + t289 * t321;
t107 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t236;
t109 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t236;
t111 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t236;
t154 = Icges(6,5) * t232 + Icges(6,6) * t231 - Icges(6,3) * t324;
t156 = Icges(6,4) * t232 + Icges(6,2) * t231 - Icges(6,6) * t324;
t158 = Icges(6,1) * t232 + Icges(6,4) * t231 - Icges(6,5) * t324;
t22 = -t107 * t324 + t109 * t231 + t111 * t232 + t154 * t236 + t156 * t165 + t158 * t166;
t234 = t283 * t312 + t300 * t311;
t291 = t300 * qJD(2);
t239 = qJD(3) * t325 - t291 * t386;
t167 = -qJD(5) * t234 - t239 * t311 + t292 * t312;
t233 = -t283 * t311 + t300 * t312;
t168 = qJD(5) * t233 + t239 * t312 + t292 * t311;
t238 = qJD(3) * t283 - t291 * t321;
t108 = Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t238;
t110 = Icges(6,4) * t168 + Icges(6,2) * t167 + Icges(6,6) * t238;
t112 = Icges(6,1) * t168 + Icges(6,4) * t167 + Icges(6,5) * t238;
t155 = Icges(6,5) * t234 + Icges(6,6) * t233 - Icges(6,3) * t325;
t157 = Icges(6,4) * t234 + Icges(6,2) * t233 - Icges(6,6) * t325;
t159 = Icges(6,1) * t234 + Icges(6,4) * t233 - Icges(6,5) * t325;
t23 = -t108 * t324 + t110 * t231 + t112 * t232 + t155 * t236 + t157 * t165 + t159 * t166;
t337 = t316 * t351;
t285 = -qJD(3) * t302 + t337 * t386;
t327 = -t303 * t312 + t311 * t375;
t199 = qJD(5) * t327 - t285 * t311 + t312 * t338;
t276 = -t303 * t311 - t312 * t375;
t200 = qJD(5) * t276 + t285 * t312 + t311 * t338;
t284 = qJD(3) * t303 + t321 * t337;
t140 = Icges(6,5) * t200 + Icges(6,6) * t199 + Icges(6,3) * t284;
t141 = Icges(6,4) * t200 + Icges(6,2) * t199 + Icges(6,6) * t284;
t142 = Icges(6,1) * t200 + Icges(6,4) * t199 + Icges(6,5) * t284;
t194 = -Icges(6,5) * t327 + Icges(6,6) * t276 + Icges(6,3) * t302;
t195 = -Icges(6,4) * t327 + Icges(6,2) * t276 + Icges(6,6) * t302;
t196 = -Icges(6,1) * t327 + Icges(6,4) * t276 + Icges(6,5) * t302;
t45 = -t140 * t324 + t141 * t231 + t142 * t232 + t165 * t195 + t166 * t196 + t194 * t236;
t76 = -t154 * t324 + t156 * t231 + t158 * t232;
t77 = -t155 * t324 + t157 * t231 + t159 * t232;
t94 = -t194 * t324 + t195 * t231 + t196 * t232;
t3 = -t22 * t326 + t23 * t300 + t76 * t290 + t77 * t292 + (-t323 * t45 + t352 * t94) * t316;
t208 = -t237 * t314 + t290 * t317;
t382 = t290 * t314;
t209 = t237 * t317 + t382;
t129 = Icges(5,5) * t209 + Icges(5,6) * t208 + Icges(5,3) * t236;
t131 = Icges(5,4) * t209 + Icges(5,2) * t208 + Icges(5,6) * t236;
t133 = Icges(5,1) * t209 + Icges(5,4) * t208 + Icges(5,5) * t236;
t40 = -t129 * t324 + t131 * t240 + t133 * t241 + t171 * t236 + t173 * t208 + t175 * t209;
t210 = -t239 * t314 + t292 * t317;
t381 = t292 * t314;
t211 = t239 * t317 + t381;
t130 = Icges(5,5) * t211 + Icges(5,6) * t210 + Icges(5,3) * t238;
t132 = Icges(5,4) * t211 + Icges(5,2) * t210 + Icges(5,6) * t238;
t134 = Icges(5,1) * t211 + Icges(5,4) * t210 + Icges(5,5) * t238;
t41 = -t130 * t324 + t132 * t240 + t134 * t241 + t172 * t236 + t174 * t208 + t176 * t209;
t254 = -t285 * t314 + t317 * t338;
t332 = t314 * t338;
t255 = t285 * t317 + t332;
t188 = Icges(5,5) * t255 + Icges(5,6) * t254 + Icges(5,3) * t284;
t189 = Icges(5,4) * t255 + Icges(5,2) * t254 + Icges(5,6) * t284;
t190 = Icges(5,1) * t255 + Icges(5,4) * t254 + Icges(5,5) * t284;
t56 = -t188 * t324 + t189 * t240 + t190 * t241 + t208 * t213 + t209 * t214 + t212 * t236;
t180 = Icges(4,5) * t237 - Icges(4,6) * t236 + Icges(4,3) * t290;
t182 = Icges(4,4) * t237 - Icges(4,2) * t236 + Icges(4,6) * t290;
t184 = Icges(4,1) * t237 - Icges(4,4) * t236 + Icges(4,5) * t290;
t62 = -t180 * t326 + t182 * t324 + t184 * t281 + t202 * t290 - t204 * t236 + t206 * t237;
t181 = Icges(4,5) * t239 - Icges(4,6) * t238 + Icges(4,3) * t292;
t183 = Icges(4,4) * t239 - Icges(4,2) * t238 + Icges(4,6) * t292;
t185 = Icges(4,1) * t239 - Icges(4,4) * t238 + Icges(4,5) * t292;
t63 = -t181 * t326 + t183 * t324 + t185 * t281 + t203 * t290 - t205 * t236 + t207 * t237;
t220 = Icges(4,5) * t285 - Icges(4,6) * t284 + Icges(4,3) * t338;
t221 = Icges(4,4) * t285 - Icges(4,2) * t284 + Icges(4,6) * t338;
t222 = Icges(4,1) * t285 - Icges(4,4) * t284 + Icges(4,5) * t338;
t84 = -t220 * t326 + t221 * t324 + t222 * t281 - t236 * t257 + t237 * t258 + t256 * t290;
t401 = t3 + (-t62 - t40) * t326 + (t409 * t352 + (-t56 - t84) * t323) * t316 + (t63 + t41) * t300 + t406 * t292 + t407 * t290;
t24 = -t107 * t325 + t109 * t233 + t111 * t234 + t154 * t238 + t156 * t167 + t158 * t168;
t25 = -t108 * t325 + t110 * t233 + t112 * t234 + t155 * t238 + t157 * t167 + t159 * t168;
t46 = -t140 * t325 + t141 * t233 + t142 * t234 + t167 * t195 + t168 * t196 + t194 * t238;
t78 = -t154 * t325 + t156 * t233 + t158 * t234;
t79 = -t155 * t325 + t157 * t233 + t159 * t234;
t95 = -t194 * t325 + t195 * t233 + t196 * t234;
t4 = -t24 * t326 + t25 * t300 + t78 * t290 + t79 * t292 + (-t323 * t46 + t352 * t95) * t316;
t42 = -t129 * t325 + t131 * t242 + t133 * t243 + t171 * t238 + t173 * t210 + t175 * t211;
t43 = -t130 * t325 + t132 * t242 + t134 * t243 + t172 * t238 + t174 * t210 + t176 * t211;
t57 = -t188 * t325 + t189 * t242 + t190 * t243 + t210 * t213 + t211 * t214 + t212 * t238;
t64 = t180 * t300 + t182 * t325 + t184 * t283 + t202 * t292 - t204 * t238 + t206 * t239;
t65 = t181 * t300 + t183 * t325 + t185 * t283 + t203 * t292 - t205 * t238 + t207 * t239;
t85 = t220 * t300 + t221 * t325 + t222 * t283 - t238 * t257 + t239 * t258 + t256 * t292;
t400 = t4 + (-t64 - t42) * t326 + (t408 * t352 + (-t57 - t85) * t323) * t316 + (t65 + t43) * t300 + t404 * t292 + t405 * t290;
t399 = 2 * m(4);
t398 = 2 * m(5);
t397 = 2 * m(6);
t396 = t236 / 0.2e1;
t395 = t238 / 0.2e1;
t394 = -t324 / 0.2e1;
t393 = -t325 / 0.2e1;
t392 = t284 / 0.2e1;
t391 = t290 / 0.2e1;
t390 = t292 / 0.2e1;
t389 = t302 / 0.2e1;
t388 = t315 / 0.2e1;
t387 = -t318 / 0.2e1;
t385 = pkin(4) * t317;
t384 = Icges(3,4) * t322;
t383 = Icges(3,4) * t323;
t378 = t315 * t316;
t377 = t316 * t318;
t113 = rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t236;
t371 = pkin(4) * t382 + pkin(8) * t236 + t237 * t385 + t113;
t114 = rSges(6,1) * t168 + rSges(6,2) * t167 + rSges(6,3) * t238;
t370 = pkin(4) * t381 + pkin(8) * t238 + t239 * t385 + t114;
t136 = rSges(5,1) * t211 + rSges(5,2) * t210 + rSges(5,3) * t238;
t170 = pkin(3) * t239 + qJ(4) * t238 - qJD(4) * t325;
t369 = -t136 - t170;
t143 = rSges(6,1) * t200 + rSges(6,2) * t199 + rSges(6,3) * t284;
t368 = pkin(4) * t332 + pkin(8) * t284 + t285 * t385 + t143;
t169 = pkin(3) * t237 + qJ(4) * t236 - qJD(4) * t324;
t229 = pkin(3) * t281 - qJ(4) * t324;
t367 = t300 * t169 + t292 * t229;
t160 = rSges(6,1) * t232 + rSges(6,2) * t231 - rSges(6,3) * t324;
t366 = -pkin(4) * t380 - pkin(8) * t324 + t281 * t385 + t160;
t161 = rSges(6,1) * t234 + rSges(6,2) * t233 - rSges(6,3) * t325;
t365 = pkin(4) * t379 - pkin(8) * t325 + t283 * t385 + t161;
t271 = -pkin(2) * t291 + pkin(7) * t292;
t247 = t319 * t271;
t364 = t319 * t170 + t247;
t270 = pkin(2) * t289 + pkin(7) * t290;
t363 = -t169 - t270;
t177 = rSges(5,1) * t241 + rSges(5,2) * t240 - rSges(5,3) * t324;
t362 = -t177 - t229;
t178 = rSges(5,1) * t243 + rSges(5,2) * t242 - rSges(5,3) * t325;
t230 = pkin(3) * t283 - qJ(4) * t325;
t361 = -t178 - t230;
t191 = rSges(5,1) * t255 + rSges(5,2) * t254 + rSges(5,3) * t284;
t218 = pkin(3) * t285 + qJ(4) * t284 + qJD(4) * t302;
t360 = -t191 - t218;
t197 = -rSges(6,1) * t327 + rSges(6,2) * t276 + rSges(6,3) * t302;
t359 = -pkin(4) * t349 + pkin(8) * t302 + t303 * t385 + t197;
t217 = rSges(5,1) * t279 + rSges(5,2) * t278 + rSges(5,3) * t302;
t275 = t303 * pkin(3) + t302 * qJ(4);
t358 = -t217 - t275;
t357 = t229 * t375 - t275 * t326;
t274 = pkin(2) * t301 + pkin(7) * t300;
t272 = t319 * t274;
t356 = t319 * t230 + t272;
t355 = t270 * t378 + t271 * t377;
t273 = pkin(2) * t299 - pkin(7) * t326;
t354 = t273 * t378 + t274 * t377;
t353 = qJD(2) * t316;
t350 = m(5) / 0.2e1 + m(6) / 0.2e1;
t346 = -t170 - t370;
t345 = -t218 - t368;
t344 = -t229 - t366;
t343 = -t230 - t365;
t342 = t169 * t375 - t218 * t326 + t290 * t275;
t341 = -t275 - t359;
t223 = rSges(4,1) * t285 - rSges(4,2) * t284 + rSges(4,3) * t338;
t297 = (pkin(2) * t323 + pkin(7) * t322) * t353;
t336 = (-t223 - t297) * t316;
t259 = t303 * rSges(4,1) - t302 * rSges(4,2) - rSges(4,3) * t375;
t304 = (pkin(2) * t322 - pkin(7) * t323) * t316;
t335 = (-t259 - t304) * t316;
t334 = t169 * t378 + t170 * t377 + t355;
t333 = t229 * t378 + t230 * t377 + t354;
t331 = (-t297 + t360) * t316;
t330 = (-t304 + t358) * t316;
t329 = (-t297 + t345) * t316;
t328 = (-t304 + t341) * t316;
t296 = (rSges(3,1) * t323 - rSges(3,2) * t322) * t353;
t295 = (Icges(3,1) * t323 - t384) * t353;
t294 = (-Icges(3,2) * t322 + t383) * t353;
t293 = (Icges(3,5) * t323 - Icges(3,6) * t322) * t353;
t288 = t319 * rSges(3,3) + (rSges(3,1) * t322 + rSges(3,2) * t323) * t316;
t287 = Icges(3,5) * t319 + (Icges(3,1) * t322 + t383) * t316;
t286 = Icges(3,6) * t319 + (Icges(3,2) * t323 + t384) * t316;
t269 = -rSges(3,1) * t291 - rSges(3,2) * t292;
t268 = rSges(3,1) * t289 - rSges(3,2) * t290;
t267 = -Icges(3,1) * t291 - Icges(3,4) * t292;
t266 = Icges(3,1) * t289 - Icges(3,4) * t290;
t265 = -Icges(3,4) * t291 - Icges(3,2) * t292;
t264 = Icges(3,4) * t289 - Icges(3,2) * t290;
t263 = -Icges(3,5) * t291 - Icges(3,6) * t292;
t262 = Icges(3,5) * t289 - Icges(3,6) * t290;
t253 = rSges(3,1) * t301 - rSges(3,2) * t300 + rSges(3,3) * t378;
t252 = rSges(3,1) * t299 + rSges(3,2) * t326 - rSges(3,3) * t377;
t251 = Icges(3,1) * t301 - Icges(3,4) * t300 + Icges(3,5) * t378;
t250 = Icges(3,1) * t299 + Icges(3,4) * t326 - Icges(3,5) * t377;
t249 = Icges(3,4) * t301 - Icges(3,2) * t300 + Icges(3,6) * t378;
t248 = Icges(3,4) * t299 + Icges(3,2) * t326 - Icges(3,6) * t377;
t224 = t230 * t338;
t219 = t300 * t229;
t216 = rSges(4,1) * t283 + rSges(4,2) * t325 + rSges(4,3) * t300;
t215 = rSges(4,1) * t281 + rSges(4,2) * t324 - rSges(4,3) * t326;
t192 = (t268 * t315 + t269 * t318) * t316;
t187 = rSges(4,1) * t239 - rSges(4,2) * t238 + rSges(4,3) * t292;
t186 = rSges(4,1) * t237 - rSges(4,2) * t236 + rSges(4,3) * t290;
t164 = -t216 * t375 - t300 * t259;
t163 = t215 * t375 - t259 * t326;
t147 = -t256 * t375 - t302 * t257 + t303 * t258;
t146 = t215 * t300 + t216 * t326;
t145 = (-t215 - t273) * t319 + t318 * t335;
t144 = t216 * t319 + t315 * t335 + t272;
t135 = rSges(5,1) * t209 + rSges(5,2) * t208 + rSges(5,3) * t236;
t126 = (t215 * t315 + t216 * t318) * t316 + t354;
t123 = (-t186 - t270) * t319 + t318 * t336;
t122 = t187 * t319 + t315 * t336 + t247;
t121 = t161 * t302 + t197 * t325;
t120 = -t160 * t302 - t197 * t324;
t119 = t212 * t302 + t213 * t278 + t214 * t279;
t106 = t194 * t302 + t195 * t276 - t196 * t327;
t105 = -t160 * t325 + t161 * t324;
t104 = t300 * t358 + t361 * t375;
t103 = t177 * t375 - t217 * t326 + t357;
t102 = (t186 * t315 + t187 * t318) * t316 + t355;
t101 = -t300 * t223 - t292 * t259 + (-t187 * t323 + t216 * t352) * t316;
t100 = -t326 * t223 + t290 * t259 + (t186 * t323 - t215 * t352) * t316;
t97 = (-t273 + t362) * t319 + t318 * t330;
t96 = t178 * t319 + t315 * t330 + t356;
t93 = -t302 * t221 + t303 * t222 - t284 * t257 + t285 * t258 + (-t220 * t323 + t256 * t352) * t316;
t92 = t177 * t300 - t326 * t361 + t219;
t91 = t186 * t300 + t187 * t326 + t215 * t292 - t216 * t290;
t88 = (t177 * t315 + t178 * t318) * t316 + t333;
t87 = t155 * t302 + t157 * t276 - t159 * t327;
t86 = t154 * t302 + t156 * t276 - t158 * t327;
t75 = t300 * t341 + t343 * t375;
t74 = -t326 * t359 + t366 * t375 + t357;
t73 = (-t135 + t363) * t319 + t318 * t331;
t72 = t136 * t319 + t315 * t331 + t364;
t71 = (-t273 + t344) * t319 + t318 * t328;
t70 = t315 * t328 + t319 * t365 + t356;
t69 = -t302 * t183 + t303 * t185 - t284 * t205 + t285 * t207 + (-t181 * t323 + t203 * t352) * t316;
t68 = -t302 * t182 + t303 * t184 - t284 * t204 + t285 * t206 + (-t180 * t323 + t202 * t352) * t316;
t67 = t188 * t302 + t189 * t278 + t190 * t279 + t212 * t284 + t213 * t254 + t214 * t255;
t66 = t300 * t366 - t326 * t343 + t219;
t61 = (t315 * t366 + t318 * t365) * t316 + t333;
t60 = (t135 * t315 + t136 * t318) * t316 + t334;
t59 = t114 * t302 + t143 * t325 + t161 * t284 - t197 * t238;
t58 = -t113 * t302 - t143 * t324 - t160 * t284 + t197 * t236;
t55 = t224 + t360 * t300 + t358 * t292 + (t178 * t352 + t323 * t369) * t316;
t54 = -t326 * t191 + t290 * t217 + (t135 * t323 + t352 * t362) * t316 + t342;
t53 = t140 * t302 + t141 * t276 - t142 * t327 + t194 * t284 + t195 * t199 + t196 * t200;
t52 = -t113 * t325 + t114 * t324 + t160 * t238 - t161 * t236;
t51 = (t363 - t371) * t319 + t318 * t329;
t50 = t315 * t329 + t319 * t370 + t364;
t49 = t135 * t300 + t177 * t292 + t290 * t361 - t326 * t369 + t367;
t48 = t130 * t302 + t132 * t278 + t134 * t279 + t172 * t284 + t174 * t254 + t176 * t255;
t47 = t129 * t302 + t131 * t278 + t133 * t279 + t171 * t284 + t173 * t254 + t175 * t255;
t44 = t106 * t319 + (t315 * t87 - t318 * t86) * t316;
t39 = -t106 * t375 + t87 * t300 - t326 * t86;
t38 = t106 * t302 - t324 * t86 - t325 * t87;
t37 = t319 * t95 + (t315 * t79 - t318 * t78) * t316;
t36 = t319 * t94 + (t315 * t77 - t318 * t76) * t316;
t35 = (t315 * t371 + t318 * t370) * t316 + t334;
t34 = t79 * t300 - t326 * t78 - t375 * t95;
t33 = t77 * t300 - t326 * t76 - t375 * t94;
t32 = t302 * t95 - t324 * t78 - t325 * t79;
t31 = t302 * t94 - t324 * t76 - t325 * t77;
t30 = t108 * t302 + t110 * t276 - t112 * t327 + t155 * t284 + t157 * t199 + t159 * t200;
t29 = t107 * t302 + t109 * t276 - t111 * t327 + t154 * t284 + t156 * t199 + t158 * t200;
t28 = t224 + t345 * t300 + t341 * t292 + (t323 * t346 + t352 * t365) * t316;
t27 = -t368 * t326 + t359 * t290 + (t323 * t371 + t344 * t352) * t316 + t342;
t26 = t319 * t93 + (t315 * t69 - t318 * t68) * t316;
t21 = t319 * t85 + (t315 * t65 - t318 * t64) * t316;
t20 = t319 * t84 + (t315 * t63 - t318 * t62) * t316;
t19 = t290 * t343 + t292 * t366 + t300 * t371 - t326 * t346 + t367;
t18 = t124 * t290 + t125 * t292 - t68 * t326 + t69 * t300 + (t147 * t352 - t323 * t93) * t316;
t15 = t319 * t67 + (t315 * t48 - t318 * t47) * t316;
t14 = t319 * t57 + (t315 * t43 - t318 * t42) * t316;
t13 = t319 * t56 + (t315 * t41 - t318 * t40) * t316;
t12 = t319 * t53 + (-t29 * t318 + t30 * t315) * t316;
t11 = t89 * t290 + t90 * t292 - t47 * t326 + t48 * t300 + (t119 * t352 - t323 * t67) * t316;
t10 = t319 * t46 + (-t24 * t318 + t25 * t315) * t316;
t9 = t319 * t45 + (-t22 * t318 + t23 * t315) * t316;
t6 = -t29 * t326 + t86 * t290 + t87 * t292 + t30 * t300 + (t106 * t352 - t323 * t53) * t316;
t5 = t106 * t284 + t236 * t86 + t238 * t87 - t29 * t324 - t30 * t325 + t302 * t53;
t2 = t236 * t78 + t238 * t79 - t24 * t324 - t25 * t325 + t284 * t95 + t302 * t46;
t1 = -t22 * t324 - t23 * t325 + t236 * t76 + t238 * t77 + t284 * t94 + t302 * t45;
t7 = [0; m(3) * t192 + m(4) * t102 + m(5) * t60 + m(6) * t35; t10 * t378 - t9 * t377 + t14 * t378 - t20 * t377 - t13 * t377 + t21 * t378 - ((-t249 * t290 + t251 * t289 - t263 * t377 + t265 * t326 + t267 * t299) * t378 - (-t248 * t290 + t250 * t289 - t262 * t377 + t264 * t326 + t266 * t299) * t377 + (-t286 * t290 + t287 * t289 - t293 * t377 + t294 * t326 + t295 * t299) * t319) * t377 + ((-t249 * t292 - t251 * t291 + t263 * t378 - t265 * t300 + t267 * t301) * t378 - (-t248 * t292 - t250 * t291 + t262 * t378 - t264 * t300 + t266 * t301) * t377 + (-t286 * t292 - t287 * t291 + t293 * t378 - t294 * t300 + t295 * t301) * t319) * t378 + (t35 * t61 + t50 * t70 + t51 * t71) * t397 + t319 * t12 + (t60 * t88 + t72 * t96 + t73 * t97) * t398 + t319 * t26 + t319 * t15 + (t102 * t126 + t122 * t144 + t123 * t145) * t399 + t319 * (t319 ^ 2 * t293 + (((t265 * t323 + t267 * t322) * t315 - (t264 * t323 + t266 * t322) * t318 + ((-t249 * t322 + t251 * t323) * t315 - (-t248 * t322 + t250 * t323) * t318) * qJD(2)) * t316 + (-t262 * t318 + t263 * t315 + t294 * t323 + t295 * t322 + (-t286 * t322 + t287 * t323) * qJD(2)) * t319) * t316) + 0.2e1 * m(3) * ((-t252 * t319 - t288 * t377) * (-t268 * t319 - t296 * t377) + (t253 * t319 - t288 * t378) * (t269 * t319 - t296 * t378) + (t252 * t315 + t253 * t318) * t316 * t192); m(4) * t91 + m(5) * t49 + m(6) * t19; (t21 / 0.2e1 + t14 / 0.2e1 + t10 / 0.2e1) * t300 - (t20 / 0.2e1 + t13 / 0.2e1 + t9 / 0.2e1) * t326 + m(6) * (t19 * t61 + t27 * t71 + t28 * t70 + t35 * t66 + t50 * t75 + t51 * t74) + m(5) * (t103 * t73 + t104 * t72 + t49 * t88 + t54 * t97 + t55 * t96 + t60 * t92) + m(4) * (t100 * t145 + t101 * t144 + t102 * t146 + t122 * t164 + t123 * t163 + t126 * t91) + (t18 / 0.2e1 + t11 / 0.2e1 + t6 / 0.2e1 + (t139 / 0.2e1 + t99 / 0.2e1) * t292 + (t138 / 0.2e1 + t98 / 0.2e1) * t290) * t319 + t36 * t391 + t37 * t390 + ((-t26 / 0.2e1 - t15 / 0.2e1 - t12 / 0.2e1) * t323 + (t44 / 0.2e1 + (t119 / 0.2e1 + t147 / 0.2e1) * t319) * t352 + (t402 * t338 + t400) * t388 + (t403 * t338 + t401) * t387 + (-t405 * t390 - t407 * t391) * t318 + (t404 * t390 + t406 * t391) * t315) * t316; (-t11 - t18 - t6) * t375 + t400 * t300 - t401 * t326 + (t404 * t300 - t405 * t326 - t408 * t375 + t34) * t292 + (t406 * t300 - t407 * t326 - t409 * t375 + t33) * t290 + (t39 + (-t119 - t147) * t375 + t402 * t300 - t403 * t326) * t338 + (t19 * t66 + t27 * t74 + t28 * t75) * t397 + (t103 * t54 + t104 * t55 + t49 * t92) * t398 + (t100 * t163 + t101 * t164 + t146 * t91) * t399; 0.2e1 * t350 * t284; m(6) * (t236 * t70 + t238 * t71 + t284 * t61 + t302 * t35 - t324 * t50 - t325 * t51) + m(5) * (t236 * t96 + t238 * t97 + t284 * t88 + t302 * t60 - t324 * t72 - t325 * t73); m(6) * (t19 * t302 + t236 * t75 + t238 * t74 - t27 * t325 - t28 * t324 + t284 * t66) + m(5) * (t103 * t238 + t104 * t236 + t284 * t92 + t302 * t49 - t324 * t55 - t325 * t54); 0.4e1 * t350 * (-t236 * t324 - t238 * t325 + t284 * t302); m(6) * t52; m(6) * (t105 * t35 + t120 * t51 + t121 * t50 + t52 * t61 + t58 * t71 + t59 * t70) + t44 * t392 + t12 * t389 + t319 * t5 / 0.2e1 + t37 * t395 + t10 * t393 + t36 * t396 + t9 * t394 + (t1 * t387 + t2 * t388) * t316; t31 * t391 - t326 * t1 / 0.2e1 + t33 * t396 + t3 * t394 + t39 * t392 + t6 * t389 + t34 * t395 + t4 * t393 + t32 * t390 + t300 * t2 / 0.2e1 + m(6) * (t105 * t19 + t120 * t27 + t121 * t28 + t52 * t66 + t58 * t74 + t59 * t75) + (t38 * t352 / 0.2e1 - t323 * t5 / 0.2e1) * t316; m(6) * (t105 * t284 + t120 * t238 + t121 * t236 + t302 * t52 - t324 * t59 - t325 * t58); (t105 * t52 + t120 * t58 + t121 * t59) * t397 + t238 * t32 - t325 * t2 + t236 * t31 - t324 * t1 + t284 * t38 + t302 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
