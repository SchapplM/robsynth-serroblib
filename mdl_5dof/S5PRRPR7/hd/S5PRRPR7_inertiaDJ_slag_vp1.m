% Calculate time derivative of joint inertia matrix for
% S5PRRPR7
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:04
% EndTime: 2019-12-05 16:35:48
% DurationCPUTime: 17.40s
% Computational Cost: add. (46825->877), mult. (133558->1248), div. (0->0), fcn. (155950->12), ass. (0->355)
t316 = sin(pkin(9));
t318 = cos(pkin(9));
t324 = cos(qJ(2));
t319 = cos(pkin(5));
t322 = sin(qJ(2));
t370 = t319 * t322;
t303 = t316 * t324 + t318 * t370;
t317 = sin(pkin(5));
t321 = sin(qJ(3));
t372 = t317 * t321;
t380 = cos(qJ(3));
t281 = t303 * t380 - t318 * t372;
t315 = sin(pkin(10));
t369 = t319 * t324;
t327 = -t316 * t322 + t318 * t369;
t377 = cos(pkin(10));
t240 = t281 * t315 + t327 * t377;
t241 = t281 * t377 - t315 * t327;
t339 = t317 * t380;
t325 = -t303 * t321 - t318 * t339;
t169 = Icges(5,5) * t241 - Icges(5,6) * t240 - Icges(5,3) * t325;
t207 = Icges(4,4) * t281 + Icges(4,2) * t325 - Icges(4,6) * t327;
t405 = -t169 + t207;
t347 = t316 * t370;
t305 = t318 * t324 - t347;
t283 = t305 * t380 + t316 * t372;
t304 = t316 * t369 + t318 * t322;
t242 = t283 * t315 - t304 * t377;
t243 = t283 * t377 + t304 * t315;
t326 = -t305 * t321 + t316 * t339;
t170 = Icges(5,5) * t243 - Icges(5,6) * t242 - Icges(5,3) * t326;
t208 = Icges(4,4) * t283 + Icges(4,2) * t326 + Icges(4,6) * t304;
t404 = -t170 + t208;
t171 = Icges(5,4) * t241 - Icges(5,2) * t240 - Icges(5,6) * t325;
t173 = Icges(5,1) * t241 - Icges(5,4) * t240 - Icges(5,5) * t325;
t205 = Icges(4,5) * t281 + Icges(4,6) * t325 - Icges(4,3) * t327;
t209 = Icges(4,1) * t281 + Icges(4,4) * t325 - Icges(4,5) * t327;
t401 = -t171 * t240 + t173 * t241 - t205 * t327 + t209 * t281 + t405 * t325;
t172 = Icges(5,4) * t243 - Icges(5,2) * t242 - Icges(5,6) * t326;
t174 = Icges(5,1) * t243 - Icges(5,4) * t242 - Icges(5,5) * t326;
t206 = Icges(4,5) * t283 + Icges(4,6) * t326 + Icges(4,3) * t304;
t210 = Icges(4,1) * t283 + Icges(4,4) * t326 + Icges(4,5) * t304;
t400 = -t172 * t240 + t174 * t241 - t206 * t327 + t210 * t281 + t404 * t325;
t399 = -t171 * t242 + t173 * t243 + t205 * t304 + t209 * t283 + t405 * t326;
t398 = -t172 * t242 + t174 * t243 + t206 * t304 + t210 * t283 + t404 * t326;
t307 = t319 * t321 + t322 * t339;
t336 = t317 * t377;
t278 = t307 * t315 + t324 * t336;
t371 = t317 * t324;
t279 = t307 * t377 - t315 * t371;
t306 = -t319 * t380 + t322 * t372;
t215 = Icges(5,5) * t279 - Icges(5,6) * t278 + Icges(5,3) * t306;
t216 = Icges(5,4) * t279 - Icges(5,2) * t278 + Icges(5,6) * t306;
t217 = Icges(5,1) * t279 - Icges(5,4) * t278 + Icges(5,5) * t306;
t110 = -t215 * t325 - t216 * t240 + t217 * t241;
t258 = Icges(4,5) * t307 - Icges(4,6) * t306 - Icges(4,3) * t371;
t259 = Icges(4,4) * t307 - Icges(4,2) * t306 - Icges(4,6) * t371;
t260 = Icges(4,1) * t307 - Icges(4,4) * t306 - Icges(4,5) * t371;
t152 = -t258 * t327 + t259 * t325 + t260 * t281;
t403 = t110 + t152;
t111 = -t215 * t326 - t216 * t242 + t217 * t243;
t153 = t258 * t304 + t259 * t326 + t260 * t283;
t402 = t111 + t153;
t351 = qJD(2) * t322;
t338 = t317 * t351;
t101 = t169 * t306 - t171 * t278 + t173 * t279;
t128 = -t205 * t371 - t306 * t207 + t307 * t209;
t397 = t101 + t128;
t102 = t170 * t306 - t172 * t278 + t174 * t279;
t129 = -t206 * t371 - t306 * t208 + t307 * t210;
t396 = t102 + t129;
t294 = t303 * qJD(2);
t350 = qJD(2) * t324;
t296 = -qJD(2) * t347 + t318 * t350;
t320 = sin(qJ(5));
t323 = cos(qJ(5));
t200 = -t241 * t320 - t323 * t325;
t201 = t241 * t323 - t320 * t325;
t131 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t240;
t133 = Icges(6,4) * t201 + Icges(6,2) * t200 + Icges(6,6) * t240;
t135 = Icges(6,1) * t201 + Icges(6,4) * t200 + Icges(6,5) * t240;
t293 = t327 * qJD(2);
t237 = qJD(3) * t325 + t293 * t380;
t212 = t237 * t377 + t294 * t315;
t236 = t281 * qJD(3) + t293 * t321;
t139 = -qJD(5) * t201 - t212 * t320 + t236 * t323;
t140 = qJD(5) * t200 + t212 * t323 + t236 * t320;
t211 = t237 * t315 - t294 * t377;
t90 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t211;
t92 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t211;
t94 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t211;
t18 = t131 * t211 + t133 * t139 + t135 * t140 + t200 * t92 + t201 * t94 + t240 * t90;
t202 = -t243 * t320 - t323 * t326;
t203 = t243 * t323 - t320 * t326;
t132 = Icges(6,5) * t203 + Icges(6,6) * t202 + Icges(6,3) * t242;
t134 = Icges(6,4) * t203 + Icges(6,2) * t202 + Icges(6,6) * t242;
t136 = Icges(6,1) * t203 + Icges(6,4) * t202 + Icges(6,5) * t242;
t295 = t304 * qJD(2);
t239 = qJD(3) * t326 - t295 * t380;
t214 = t239 * t377 + t296 * t315;
t238 = qJD(3) * t283 - t295 * t321;
t141 = -qJD(5) * t203 - t214 * t320 + t238 * t323;
t142 = qJD(5) * t202 + t214 * t323 + t238 * t320;
t213 = t239 * t315 - t296 * t377;
t91 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t213;
t93 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t213;
t95 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t213;
t19 = t132 * t211 + t134 * t139 + t136 * t140 + t200 * t93 + t201 * t95 + t240 * t91;
t246 = t279 * t323 + t306 * t320;
t337 = t317 * t350;
t285 = -qJD(3) * t306 + t337 * t380;
t257 = t285 * t377 + t315 * t338;
t284 = qJD(3) * t307 + t321 * t337;
t181 = -qJD(5) * t246 - t257 * t320 + t284 * t323;
t245 = -t279 * t320 + t306 * t323;
t182 = qJD(5) * t245 + t257 * t323 + t284 * t320;
t256 = t285 * t315 - t336 * t351;
t122 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t256;
t123 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t256;
t124 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t256;
t177 = Icges(6,5) * t246 + Icges(6,6) * t245 + Icges(6,3) * t278;
t178 = Icges(6,4) * t246 + Icges(6,2) * t245 + Icges(6,6) * t278;
t179 = Icges(6,1) * t246 + Icges(6,4) * t245 + Icges(6,5) * t278;
t37 = t122 * t240 + t123 * t200 + t124 * t201 + t139 * t178 + t140 * t179 + t177 * t211;
t60 = t131 * t240 + t133 * t200 + t135 * t201;
t61 = t132 * t240 + t134 * t200 + t136 * t201;
t82 = t177 * t240 + t178 * t200 + t179 * t201;
t3 = -t18 * t327 + t19 * t304 + t60 * t294 + t61 * t296 + (-t324 * t37 + t351 * t82) * t317;
t143 = Icges(5,5) * t212 - Icges(5,6) * t211 + Icges(5,3) * t236;
t145 = Icges(5,4) * t212 - Icges(5,2) * t211 + Icges(5,6) * t236;
t147 = Icges(5,1) * t212 - Icges(5,4) * t211 + Icges(5,5) * t236;
t49 = -t143 * t325 - t145 * t240 + t147 * t241 + t169 * t236 - t171 * t211 + t173 * t212;
t144 = Icges(5,5) * t214 - Icges(5,6) * t213 + Icges(5,3) * t238;
t146 = Icges(5,4) * t214 - Icges(5,2) * t213 + Icges(5,6) * t238;
t148 = Icges(5,1) * t214 - Icges(5,4) * t213 + Icges(5,5) * t238;
t50 = -t144 * t325 - t146 * t240 + t148 * t241 + t170 * t236 - t172 * t211 + t174 * t212;
t191 = Icges(5,5) * t257 - Icges(5,6) * t256 + Icges(5,3) * t284;
t192 = Icges(5,4) * t257 - Icges(5,2) * t256 + Icges(5,6) * t284;
t193 = Icges(5,1) * t257 - Icges(5,4) * t256 + Icges(5,5) * t284;
t58 = -t191 * t325 - t192 * t240 + t193 * t241 - t211 * t216 + t212 * t217 + t215 * t236;
t183 = Icges(4,5) * t237 - Icges(4,6) * t236 + Icges(4,3) * t294;
t185 = Icges(4,4) * t237 - Icges(4,2) * t236 + Icges(4,6) * t294;
t187 = Icges(4,1) * t237 - Icges(4,4) * t236 + Icges(4,5) * t294;
t69 = -t183 * t327 + t185 * t325 + t187 * t281 + t205 * t294 - t207 * t236 + t209 * t237;
t184 = Icges(4,5) * t239 - Icges(4,6) * t238 + Icges(4,3) * t296;
t186 = Icges(4,4) * t239 - Icges(4,2) * t238 + Icges(4,6) * t296;
t188 = Icges(4,1) * t239 - Icges(4,4) * t238 + Icges(4,5) * t296;
t70 = -t184 * t327 + t186 * t325 + t188 * t281 + t206 * t294 - t208 * t236 + t210 * t237;
t223 = Icges(4,5) * t285 - Icges(4,6) * t284 + Icges(4,3) * t338;
t224 = Icges(4,4) * t285 - Icges(4,2) * t284 + Icges(4,6) * t338;
t225 = Icges(4,1) * t285 - Icges(4,4) * t284 + Icges(4,5) * t338;
t98 = -t223 * t327 + t224 * t325 + t225 * t281 - t236 * t259 + t237 * t260 + t258 * t294;
t395 = t3 + (-t49 - t69) * t327 + (t403 * t351 + (-t58 - t98) * t324) * t317 + (t50 + t70) * t304 + t400 * t296 + t401 * t294;
t20 = t131 * t213 + t133 * t141 + t135 * t142 + t202 * t92 + t203 * t94 + t242 * t90;
t21 = t132 * t213 + t134 * t141 + t136 * t142 + t202 * t93 + t203 * t95 + t242 * t91;
t38 = t122 * t242 + t123 * t202 + t124 * t203 + t141 * t178 + t142 * t179 + t177 * t213;
t62 = t131 * t242 + t133 * t202 + t135 * t203;
t63 = t132 * t242 + t134 * t202 + t136 * t203;
t83 = t177 * t242 + t178 * t202 + t179 * t203;
t4 = -t20 * t327 + t21 * t304 + t62 * t294 + t63 * t296 + (-t324 * t38 + t351 * t83) * t317;
t51 = -t143 * t326 - t145 * t242 + t147 * t243 + t169 * t238 - t171 * t213 + t173 * t214;
t52 = -t144 * t326 - t146 * t242 + t148 * t243 + t170 * t238 - t172 * t213 + t174 * t214;
t59 = -t191 * t326 - t192 * t242 + t193 * t243 - t213 * t216 + t214 * t217 + t215 * t238;
t71 = t183 * t304 + t185 * t326 + t187 * t283 + t205 * t296 - t207 * t238 + t209 * t239;
t72 = t184 * t304 + t186 * t326 + t188 * t283 + t206 * t296 - t208 * t238 + t210 * t239;
t99 = t223 * t304 + t224 * t326 + t225 * t283 - t238 * t259 + t239 * t260 + t258 * t296;
t394 = t4 + (-t51 - t71) * t327 + (t402 * t351 + (-t59 - t99) * t324) * t317 + (t52 + t72) * t304 + t398 * t296 + t399 * t294;
t393 = 2 * m(4);
t392 = 2 * m(5);
t391 = 2 * m(6);
t390 = t211 / 0.2e1;
t389 = t213 / 0.2e1;
t388 = t240 / 0.2e1;
t387 = t242 / 0.2e1;
t386 = t256 / 0.2e1;
t385 = t278 / 0.2e1;
t384 = t294 / 0.2e1;
t383 = t296 / 0.2e1;
t382 = t316 / 0.2e1;
t381 = -t318 / 0.2e1;
t96 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t211;
t379 = pkin(4) * t212 + pkin(8) * t211 + t96;
t97 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t213;
t378 = pkin(4) * t214 + pkin(8) * t213 + t97;
t376 = Icges(3,4) * t322;
t375 = Icges(3,4) * t324;
t374 = t316 * t317;
t373 = t317 * t318;
t125 = rSges(6,1) * t182 + rSges(6,2) * t181 + rSges(6,3) * t256;
t368 = pkin(4) * t257 + pkin(8) * t256 + t125;
t137 = rSges(6,1) * t201 + rSges(6,2) * t200 + rSges(6,3) * t240;
t367 = pkin(4) * t241 + pkin(8) * t240 + t137;
t138 = rSges(6,1) * t203 + rSges(6,2) * t202 + rSges(6,3) * t242;
t366 = pkin(4) * t243 + pkin(8) * t242 + t138;
t150 = rSges(5,1) * t214 - rSges(5,2) * t213 + rSges(5,3) * t238;
t168 = pkin(3) * t239 + qJ(4) * t238 - qJD(4) * t326;
t365 = -t150 - t168;
t167 = pkin(3) * t237 + qJ(4) * t236 - qJD(4) * t325;
t232 = pkin(3) * t281 - qJ(4) * t325;
t364 = t304 * t167 + t296 * t232;
t273 = -pkin(2) * t295 + pkin(7) * t296;
t249 = t319 * t273;
t363 = t319 * t168 + t249;
t272 = pkin(2) * t293 + pkin(7) * t294;
t362 = -t167 - t272;
t175 = rSges(5,1) * t241 - rSges(5,2) * t240 - rSges(5,3) * t325;
t361 = -t175 - t232;
t176 = rSges(5,1) * t243 - rSges(5,2) * t242 - rSges(5,3) * t326;
t233 = pkin(3) * t283 - qJ(4) * t326;
t360 = -t176 - t233;
t180 = rSges(6,1) * t246 + rSges(6,2) * t245 + rSges(6,3) * t278;
t359 = pkin(4) * t279 + pkin(8) * t278 + t180;
t194 = rSges(5,1) * t257 - rSges(5,2) * t256 + rSges(5,3) * t284;
t221 = pkin(3) * t285 + qJ(4) * t284 + qJD(4) * t306;
t358 = -t194 - t221;
t220 = rSges(5,1) * t279 - rSges(5,2) * t278 + rSges(5,3) * t306;
t277 = pkin(3) * t307 + qJ(4) * t306;
t357 = -t220 - t277;
t356 = t232 * t371 - t277 * t327;
t276 = pkin(2) * t305 + pkin(7) * t304;
t274 = t319 * t276;
t355 = t319 * t233 + t274;
t354 = t272 * t374 + t273 * t373;
t275 = pkin(2) * t303 - pkin(7) * t327;
t353 = t275 * t374 + t276 * t373;
t352 = qJD(2) * t317;
t349 = m(5) / 0.2e1 + m(6) / 0.2e1;
t348 = -t168 - t378;
t345 = -t221 - t368;
t344 = -t232 - t367;
t343 = -t233 - t366;
t342 = t167 * t371 - t221 * t327 + t294 * t277;
t341 = -t277 - t359;
t226 = rSges(4,1) * t285 - rSges(4,2) * t284 + rSges(4,3) * t338;
t301 = (pkin(2) * t324 + pkin(7) * t322) * t352;
t335 = (-t226 - t301) * t317;
t261 = t307 * rSges(4,1) - t306 * rSges(4,2) - rSges(4,3) * t371;
t308 = (pkin(2) * t322 - pkin(7) * t324) * t317;
t334 = (-t261 - t308) * t317;
t333 = t167 * t374 + t168 * t373 + t354;
t332 = t232 * t374 + t233 * t373 + t353;
t331 = (-t301 + t358) * t317;
t330 = (-t308 + t357) * t317;
t329 = (-t301 + t345) * t317;
t328 = (-t308 + t341) * t317;
t300 = (rSges(3,1) * t324 - rSges(3,2) * t322) * t352;
t299 = (Icges(3,1) * t324 - t376) * t352;
t298 = (-Icges(3,2) * t322 + t375) * t352;
t297 = (Icges(3,5) * t324 - Icges(3,6) * t322) * t352;
t290 = t319 * rSges(3,3) + (rSges(3,1) * t322 + rSges(3,2) * t324) * t317;
t289 = Icges(3,5) * t319 + (Icges(3,1) * t322 + t375) * t317;
t288 = Icges(3,6) * t319 + (Icges(3,2) * t324 + t376) * t317;
t271 = -rSges(3,1) * t295 - rSges(3,2) * t296;
t270 = rSges(3,1) * t293 - rSges(3,2) * t294;
t269 = -Icges(3,1) * t295 - Icges(3,4) * t296;
t268 = Icges(3,1) * t293 - Icges(3,4) * t294;
t267 = -Icges(3,4) * t295 - Icges(3,2) * t296;
t266 = Icges(3,4) * t293 - Icges(3,2) * t294;
t265 = -Icges(3,5) * t295 - Icges(3,6) * t296;
t264 = Icges(3,5) * t293 - Icges(3,6) * t294;
t255 = rSges(3,1) * t305 - rSges(3,2) * t304 + rSges(3,3) * t374;
t254 = rSges(3,1) * t303 + rSges(3,2) * t327 - rSges(3,3) * t373;
t253 = Icges(3,1) * t305 - Icges(3,4) * t304 + Icges(3,5) * t374;
t252 = Icges(3,1) * t303 + Icges(3,4) * t327 - Icges(3,5) * t373;
t251 = Icges(3,4) * t305 - Icges(3,2) * t304 + Icges(3,6) * t374;
t250 = Icges(3,4) * t303 + Icges(3,2) * t327 - Icges(3,6) * t373;
t227 = t233 * t338;
t222 = t304 * t232;
t219 = rSges(4,1) * t283 + rSges(4,2) * t326 + rSges(4,3) * t304;
t218 = rSges(4,1) * t281 + rSges(4,2) * t325 - rSges(4,3) * t327;
t197 = (t270 * t316 + t271 * t318) * t317;
t190 = rSges(4,1) * t239 - rSges(4,2) * t238 + rSges(4,3) * t296;
t189 = rSges(4,1) * t237 - rSges(4,2) * t236 + rSges(4,3) * t294;
t166 = -t219 * t371 - t304 * t261;
t165 = t218 * t371 - t261 * t327;
t157 = -t258 * t371 - t306 * t259 + t307 * t260;
t156 = t218 * t304 + t219 * t327;
t155 = (-t218 - t275) * t319 + t318 * t334;
t154 = t219 * t319 + t316 * t334 + t274;
t149 = rSges(5,1) * t212 - rSges(5,2) * t211 + rSges(5,3) * t236;
t130 = (t218 * t316 + t219 * t318) * t317 + t353;
t127 = (-t189 - t272) * t319 + t318 * t335;
t126 = t190 * t319 + t316 * t335 + t249;
t121 = t215 * t306 - t216 * t278 + t217 * t279;
t116 = t304 * t357 + t360 * t371;
t115 = t175 * t371 - t220 * t327 + t356;
t114 = (t189 * t316 + t190 * t318) * t317 + t354;
t113 = -t304 * t226 - t296 * t261 + (-t190 * t324 + t219 * t351) * t317;
t112 = -t327 * t226 + t294 * t261 + (t189 * t324 - t218 * t351) * t317;
t109 = (-t275 + t361) * t319 + t318 * t330;
t108 = t176 * t319 + t316 * t330 + t355;
t107 = -t306 * t224 + t307 * t225 - t284 * t259 + t285 * t260 + (-t223 * t324 + t258 * t351) * t317;
t106 = t138 * t278 - t180 * t242;
t105 = -t137 * t278 + t180 * t240;
t104 = t175 * t304 - t327 * t360 + t222;
t103 = t189 * t304 + t190 * t327 + t218 * t296 - t219 * t294;
t100 = (t175 * t316 + t176 * t318) * t317 + t332;
t89 = t177 * t278 + t178 * t245 + t179 * t246;
t84 = t137 * t242 - t138 * t240;
t81 = (-t149 + t362) * t319 + t318 * t331;
t80 = t150 * t319 + t316 * t331 + t363;
t79 = t304 * t341 + t343 * t371;
t78 = -t327 * t359 + t367 * t371 + t356;
t77 = (-t275 + t344) * t319 + t318 * t328;
t76 = t316 * t328 + t319 * t366 + t355;
t75 = -t306 * t186 + t307 * t188 - t284 * t208 + t285 * t210 + (-t184 * t324 + t206 * t351) * t317;
t74 = -t306 * t185 + t307 * t187 - t284 * t207 + t285 * t209 + (-t183 * t324 + t205 * t351) * t317;
t73 = t191 * t306 - t192 * t278 + t193 * t279 + t215 * t284 - t216 * t256 + t217 * t257;
t68 = t132 * t278 + t134 * t245 + t136 * t246;
t67 = t131 * t278 + t133 * t245 + t135 * t246;
t66 = (t149 * t316 + t150 * t318) * t317 + t333;
t65 = t304 * t367 - t327 * t343 + t222;
t64 = (t316 * t367 + t318 * t366) * t317 + t332;
t57 = t227 + t358 * t304 + t357 * t296 + (t176 * t351 + t324 * t365) * t317;
t56 = -t327 * t194 + t294 * t220 + (t149 * t324 + t351 * t361) * t317 + t342;
t55 = t149 * t304 + t175 * t296 + t294 * t360 - t327 * t365 + t364;
t54 = t144 * t306 - t146 * t278 + t148 * t279 + t170 * t284 - t172 * t256 + t174 * t257;
t53 = t143 * t306 - t145 * t278 + t147 * t279 + t169 * t284 - t171 * t256 + t173 * t257;
t48 = (t362 - t379) * t319 + t318 * t329;
t47 = t316 * t329 + t319 * t378 + t363;
t46 = -t125 * t242 + t138 * t256 - t180 * t213 + t278 * t97;
t45 = t125 * t240 - t137 * t256 + t180 * t211 - t278 * t96;
t44 = t122 * t278 + t123 * t245 + t124 * t246 + t177 * t256 + t178 * t181 + t179 * t182;
t43 = t137 * t213 - t138 * t211 - t240 * t97 + t242 * t96;
t42 = (t316 * t379 + t318 * t378) * t317 + t333;
t41 = t107 * t319 + (t316 * t75 - t318 * t74) * t317;
t40 = t227 + t345 * t304 + t341 * t296 + (t324 * t348 + t351 * t366) * t317;
t39 = -t368 * t327 + t359 * t294 + (t324 * t379 + t344 * t351) * t317 + t342;
t36 = t319 * t99 + (t316 * t72 - t318 * t71) * t317;
t35 = t319 * t98 + (t316 * t70 - t318 * t69) * t317;
t34 = t319 * t89 + (t316 * t68 - t318 * t67) * t317;
t33 = t68 * t304 - t327 * t67 - t371 * t89;
t32 = t240 * t67 + t242 * t68 + t278 * t89;
t31 = t319 * t83 + (t316 * t63 - t318 * t62) * t317;
t30 = t319 * t82 + (t316 * t61 - t318 * t60) * t317;
t29 = t63 * t304 - t327 * t62 - t371 * t83;
t28 = t61 * t304 - t327 * t60 - t371 * t82;
t27 = t240 * t62 + t242 * t63 + t278 * t83;
t26 = t240 * t60 + t242 * t61 + t278 * t82;
t25 = t128 * t294 + t129 * t296 - t74 * t327 + t75 * t304 + (-t107 * t324 + t157 * t351) * t317;
t24 = t294 * t343 + t296 * t367 + t304 * t379 - t327 * t348 + t364;
t23 = t132 * t256 + t134 * t181 + t136 * t182 + t245 * t93 + t246 * t95 + t278 * t91;
t22 = t131 * t256 + t133 * t181 + t135 * t182 + t245 * t92 + t246 * t94 + t278 * t90;
t15 = t319 * t73 + (t316 * t54 - t318 * t53) * t317;
t14 = t319 * t59 + (t316 * t52 - t318 * t51) * t317;
t13 = t319 * t58 + (t316 * t50 - t318 * t49) * t317;
t12 = t101 * t294 + t102 * t296 - t53 * t327 + t54 * t304 + (t121 * t351 - t324 * t73) * t317;
t9 = t319 * t44 + (-t22 * t318 + t23 * t316) * t317;
t8 = t319 * t38 + (-t20 * t318 + t21 * t316) * t317;
t7 = t319 * t37 + (-t18 * t318 + t19 * t316) * t317;
t6 = -t22 * t327 + t23 * t304 + t67 * t294 + t68 * t296 + (-t324 * t44 + t351 * t89) * t317;
t5 = t211 * t67 + t213 * t68 + t22 * t240 + t23 * t242 + t256 * t89 + t278 * t44;
t2 = t20 * t240 + t21 * t242 + t211 * t62 + t213 * t63 + t256 * t83 + t278 * t38;
t1 = t18 * t240 + t19 * t242 + t211 * t60 + t213 * t61 + t256 * t82 + t278 * t37;
t10 = [0; m(3) * t197 + m(4) * t114 + m(5) * t66 + m(6) * t42; t8 * t374 - t7 * t373 + t14 * t374 - t13 * t373 + t36 * t374 - t35 * t373 - ((-t251 * t294 + t253 * t293 - t265 * t373 + t267 * t327 + t269 * t303) * t374 - (-t250 * t294 + t252 * t293 - t264 * t373 + t266 * t327 + t268 * t303) * t373 + (-t294 * t288 + t293 * t289 - t297 * t373 + t298 * t327 + t303 * t299) * t319) * t373 + ((-t251 * t296 - t253 * t295 + t265 * t374 - t267 * t304 + t269 * t305) * t374 - (-t250 * t296 - t252 * t295 + t264 * t374 - t266 * t304 + t268 * t305) * t373 + (-t296 * t288 - t295 * t289 + t297 * t374 - t304 * t298 + t305 * t299) * t319) * t374 + (t42 * t64 + t47 * t76 + t48 * t77) * t391 + t319 * t9 + (t100 * t66 + t108 * t80 + t109 * t81) * t392 + (t114 * t130 + t126 * t154 + t127 * t155) * t393 + t319 * t41 + t319 * t15 + t319 * (t319 ^ 2 * t297 + (((t267 * t324 + t269 * t322) * t316 - (t266 * t324 + t268 * t322) * t318 + ((-t251 * t322 + t253 * t324) * t316 - (-t250 * t322 + t252 * t324) * t318) * qJD(2)) * t317 + (-t264 * t318 + t265 * t316 + t298 * t324 + t299 * t322 + (-t288 * t322 + t289 * t324) * qJD(2)) * t319) * t317) + 0.2e1 * m(3) * ((-t254 * t319 - t290 * t373) * (-t270 * t319 - t300 * t373) + (t255 * t319 - t290 * t374) * (t271 * t319 - t300 * t374) + (t254 * t316 + t255 * t318) * t317 * t197); m(4) * t103 + m(5) * t55 + m(6) * t24; (t36 / 0.2e1 + t8 / 0.2e1 + t14 / 0.2e1) * t304 - (t7 / 0.2e1 + t13 / 0.2e1 + t35 / 0.2e1) * t327 + m(6) * (t24 * t64 + t39 * t77 + t40 * t76 + t42 * t65 + t47 * t79 + t48 * t78) + m(5) * (t100 * t55 + t104 * t66 + t108 * t57 + t109 * t56 + t115 * t81 + t116 * t80) + m(4) * (t103 * t130 + t112 * t155 + t113 * t154 + t114 * t156 + t126 * t166 + t127 * t165) + (t6 / 0.2e1 + t12 / 0.2e1 + t25 / 0.2e1 + (t153 / 0.2e1 + t111 / 0.2e1) * t296 + (t152 / 0.2e1 + t110 / 0.2e1) * t294) * t319 + t30 * t384 + t31 * t383 + ((-t9 / 0.2e1 - t15 / 0.2e1 - t41 / 0.2e1) * t324 + (t34 / 0.2e1 + (t121 / 0.2e1 + t157 / 0.2e1) * t319) * t351 + (t396 * t338 + t394) * t382 + (t397 * t338 + t395) * t381 + (-t399 * t383 - t401 * t384) * t318 + (t398 * t383 + t400 * t384) * t316) * t317; (-t12 - t25 - t6) * t371 + t394 * t304 - t395 * t327 + (t398 * t304 - t399 * t327 - t402 * t371 + t29) * t296 + (t400 * t304 - t401 * t327 - t403 * t371 + t28) * t294 + (t33 + (-t121 - t157) * t371 + t396 * t304 - t397 * t327) * t338 + (t24 * t65 + t39 * t78 + t40 * t79) * t391 + (t104 * t55 + t115 * t56 + t116 * t57) * t392 + (t103 * t156 + t112 * t165 + t113 * t166) * t393; 0.2e1 * t349 * t284; m(6) * (t236 * t76 + t238 * t77 + t284 * t64 + t306 * t42 - t325 * t47 - t326 * t48) + m(5) * (t100 * t284 + t108 * t236 + t109 * t238 + t306 * t66 - t325 * t80 - t326 * t81); m(6) * (t236 * t79 + t238 * t78 + t24 * t306 + t284 * t65 - t325 * t40 - t326 * t39) + m(5) * (t104 * t284 + t115 * t238 + t116 * t236 + t306 * t55 - t325 * t57 - t326 * t56); 0.4e1 * t349 * (-t236 * t325 - t238 * t326 + t284 * t306); m(6) * t43; m(6) * (t105 * t48 + t106 * t47 + t42 * t84 + t43 * t64 + t45 * t77 + t46 * t76) + t31 * t389 + t8 * t387 + t30 * t390 + t7 * t388 + t319 * t5 / 0.2e1 + t34 * t386 + t9 * t385 + (t1 * t381 + t2 * t382) * t317; m(6) * (t105 * t39 + t106 * t40 + t24 * t84 + t43 * t65 + t45 * t78 + t46 * t79) + t33 * t386 + t6 * t385 + t26 * t384 - t327 * t1 / 0.2e1 + t28 * t390 + t3 * t388 + t29 * t389 + t4 * t387 + t27 * t383 + t304 * t2 / 0.2e1 + (t32 * t351 / 0.2e1 - t324 * t5 / 0.2e1) * t317; m(6) * (t105 * t238 + t106 * t236 + t284 * t84 + t306 * t43 - t325 * t46 - t326 * t45); (t105 * t45 + t106 * t46 + t43 * t84) * t391 + t213 * t27 + t242 * t2 + t211 * t26 + t240 * t1 + t256 * t32 + t278 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
