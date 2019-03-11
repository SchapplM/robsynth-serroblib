% Calculate joint inertia matrix for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR13_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR13_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:52:03
% EndTime: 2019-03-09 19:52:24
% DurationCPUTime: 8.44s
% Computational Cost: add. (54575->713), mult. (129834->979), div. (0->0), fcn. (170166->16), ass. (0->336)
t316 = cos(qJ(1));
t308 = cos(pkin(6));
t315 = cos(qJ(2));
t381 = t315 * t316;
t312 = sin(qJ(2));
t313 = sin(qJ(1));
t383 = t313 * t312;
t326 = t308 * t381 - t383;
t306 = sin(pkin(6));
t395 = cos(pkin(7));
t345 = t306 * t395;
t394 = sin(pkin(7));
t271 = -t316 * t345 - t326 * t394;
t382 = t313 * t315;
t384 = t312 * t316;
t290 = -t308 * t382 - t384;
t272 = -t290 * t394 + t313 * t345;
t388 = t306 * t313;
t289 = t308 * t384 + t382;
t311 = sin(qJ(3));
t322 = t326 * t395;
t344 = t306 * t394;
t398 = cos(qJ(3));
t251 = t289 * t398 + (-t316 * t344 + t322) * t311;
t365 = pkin(13) + qJ(5);
t303 = sin(t365);
t343 = cos(t365);
t216 = t251 * t303 - t271 * t343;
t291 = -t308 * t383 + t381;
t253 = t291 * t398 + (t290 * t395 + t313 * t344) * t311;
t218 = t253 * t303 - t272 * t343;
t270 = t308 * t394 * t311 + (t311 * t315 * t395 + t312 * t398) * t306;
t288 = t308 * t395 - t315 * t344;
t235 = t270 * t303 - t288 * t343;
t217 = t251 * t343 + t271 * t303;
t334 = t398 * t394;
t330 = t306 * t334;
t250 = t289 * t311 + t316 * t330 - t322 * t398;
t310 = sin(qJ(6));
t314 = cos(qJ(6));
t176 = -t217 * t310 + t250 * t314;
t177 = t217 * t314 + t250 * t310;
t107 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t216;
t109 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t216;
t111 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t216;
t34 = t107 * t216 + t109 * t176 + t111 * t177;
t219 = t253 * t343 + t272 * t303;
t335 = t395 * t398;
t252 = -t290 * t335 + t291 * t311 - t313 * t330;
t178 = -t219 * t310 + t252 * t314;
t179 = t219 * t314 + t252 * t310;
t108 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t218;
t110 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t218;
t112 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t218;
t35 = t108 * t216 + t110 * t176 + t112 * t177;
t236 = t270 * t343 + t288 * t303;
t387 = t306 * t315;
t389 = t306 * t312;
t269 = -t308 * t334 + t311 * t389 - t335 * t387;
t209 = -t236 * t310 + t269 * t314;
t210 = t236 * t314 + t269 * t310;
t138 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t235;
t139 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t235;
t140 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t235;
t49 = t138 * t216 + t139 * t176 + t140 * t177;
t1 = t216 * t34 + t218 * t35 + t235 * t49;
t404 = t1 / 0.2e1;
t36 = t107 * t218 + t109 * t178 + t111 * t179;
t37 = t108 * t218 + t110 * t178 + t112 * t179;
t50 = t138 * t218 + t139 * t178 + t140 * t179;
t2 = t216 * t36 + t218 * t37 + t235 * t50;
t403 = t2 / 0.2e1;
t42 = t107 * t235 + t109 * t209 + t111 * t210;
t43 = t108 * t235 + t110 * t209 + t112 * t210;
t61 = t235 * t138 + t209 * t139 + t210 * t140;
t57 = t61 * t235;
t9 = t42 * t216 + t43 * t218 + t57;
t402 = t9 / 0.2e1;
t401 = t216 / 0.2e1;
t400 = t218 / 0.2e1;
t399 = t235 / 0.2e1;
t397 = pkin(5) * t217;
t307 = cos(pkin(13));
t302 = pkin(4) * t307 + pkin(3);
t396 = -pkin(3) + t302;
t309 = -pkin(11) - qJ(4);
t393 = t250 * t309;
t305 = sin(pkin(13));
t392 = t271 * t305;
t391 = t272 * t305;
t390 = t288 * t305;
t386 = t306 * t316;
t385 = t308 * t315;
t332 = -rSges(7,1) * t177 - rSges(7,2) * t176;
t113 = rSges(7,3) * t216 - t332;
t380 = pkin(12) * t216 + t113 + t397;
t114 = t179 * rSges(7,1) + t178 * rSges(7,2) + t218 * rSges(7,3);
t379 = t219 * pkin(5) + t218 * pkin(12) + t114;
t240 = t250 * qJ(4);
t364 = pkin(4) * t392;
t136 = t251 * t396 - t240 + t364 - t393;
t203 = pkin(3) * t251 + t240;
t194 = t272 * t203;
t378 = t272 * t136 + t194;
t204 = t253 * pkin(3) + t252 * qJ(4);
t349 = pkin(4) * t391 - t252 * t309 + t253 * t302;
t137 = -t204 + t349;
t196 = t288 * t204;
t377 = t288 * t137 + t196;
t376 = -t136 - t203;
t375 = -t137 - t204;
t141 = rSges(7,1) * t210 + rSges(7,2) * t209 + rSges(7,3) * t235;
t374 = pkin(5) * t236 + pkin(12) * t235 + t141;
t222 = -t251 * t305 + t271 * t307;
t223 = t251 * t307 + t392;
t156 = rSges(5,1) * t223 + rSges(5,2) * t222 + rSges(5,3) * t250;
t373 = -t156 - t203;
t224 = -t253 * t305 + t272 * t307;
t225 = t253 * t307 + t391;
t157 = t225 * rSges(5,1) + t224 * rSges(5,2) + t252 * rSges(5,3);
t372 = -t157 - t204;
t169 = pkin(4) * t390 + t396 * t270 + (-qJ(4) - t309) * t269;
t233 = pkin(3) * t270 + qJ(4) * t269;
t215 = t271 * t233;
t371 = t271 * t169 + t215;
t370 = -t169 - t233;
t238 = -t270 * t305 + t288 * t307;
t239 = t270 * t307 + t390;
t183 = rSges(5,1) * t239 + rSges(5,2) * t238 + rSges(5,3) * t269;
t369 = -t183 - t233;
t256 = t291 * pkin(2) + pkin(10) * t272;
t254 = t308 * t256;
t368 = t308 * t204 + t254;
t255 = pkin(2) * t289 + pkin(10) * t271;
t367 = t255 * t388 + t256 * t386;
t366 = t316 * pkin(1) + pkin(9) * t388;
t142 = Icges(6,5) * t217 - Icges(6,6) * t216 + Icges(6,3) * t250;
t144 = Icges(6,4) * t217 - Icges(6,2) * t216 + Icges(6,6) * t250;
t146 = Icges(6,1) * t217 - Icges(6,4) * t216 + Icges(6,5) * t250;
t63 = t142 * t250 - t144 * t216 + t146 * t217;
t143 = Icges(6,5) * t219 - Icges(6,6) * t218 + Icges(6,3) * t252;
t145 = Icges(6,4) * t219 - Icges(6,2) * t218 + Icges(6,6) * t252;
t147 = Icges(6,1) * t219 - Icges(6,4) * t218 + Icges(6,5) * t252;
t64 = t143 * t250 - t145 * t216 + t147 * t217;
t172 = Icges(6,5) * t236 - Icges(6,6) * t235 + Icges(6,3) * t269;
t173 = Icges(6,4) * t236 - Icges(6,2) * t235 + Icges(6,6) * t269;
t174 = Icges(6,1) * t236 - Icges(6,4) * t235 + Icges(6,5) * t269;
t81 = t172 * t250 - t173 * t216 + t174 * t217;
t13 = t250 * t63 + t252 * t64 + t269 * t81;
t3 = t250 * t34 + t252 * t35 + t269 * t49;
t363 = t3 / 0.2e1 + t13 / 0.2e1;
t65 = t142 * t252 - t144 * t218 + t146 * t219;
t66 = t143 * t252 - t145 * t218 + t147 * t219;
t82 = t172 * t252 - t173 * t218 + t174 * t219;
t14 = t250 * t65 + t252 * t66 + t269 * t82;
t4 = t250 * t36 + t252 * t37 + t269 * t50;
t362 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t271 * t63 + t272 * t64 + t288 * t81;
t5 = t271 * t34 + t272 * t35 + t288 * t49;
t361 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t271 * t65 + t272 * t66 + t288 * t82;
t6 = t271 * t36 + t272 * t37 + t288 * t50;
t360 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t81 * t308 + (t313 * t64 - t316 * t63) * t306;
t7 = t49 * t308 + (t313 * t35 - t316 * t34) * t306;
t359 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t82 * t308 + (t313 * t66 - t316 * t65) * t306;
t8 = t50 * t308 + (t313 * t37 - t316 * t36) * t306;
t358 = t8 / 0.2e1 + t18 / 0.2e1;
t58 = t61 * t269;
t10 = t42 * t250 + t43 * t252 + t58;
t71 = t142 * t269 - t144 * t235 + t146 * t236;
t72 = t143 * t269 - t145 * t235 + t147 * t236;
t92 = t269 * t172 - t235 * t173 + t236 * t174;
t89 = t92 * t269;
t23 = t71 * t250 + t72 * t252 + t89;
t357 = t10 / 0.2e1 + t23 / 0.2e1;
t59 = t61 * t288;
t11 = t42 * t271 + t43 * t272 + t59;
t90 = t92 * t288;
t24 = t71 * t271 + t72 * t272 + t90;
t356 = t11 / 0.2e1 + t24 / 0.2e1;
t60 = t61 * t308;
t12 = t60 + (t43 * t313 - t42 * t316) * t306;
t91 = t92 * t308;
t25 = t91 + (t72 * t313 - t71 * t316) * t306;
t355 = t12 / 0.2e1 + t25 / 0.2e1;
t354 = t42 / 0.2e1 + t49 / 0.2e1;
t353 = t43 / 0.2e1 + t50 / 0.2e1;
t352 = t308 * t137 + t368;
t333 = -rSges(6,1) * t217 + rSges(6,2) * t216;
t148 = rSges(6,3) * t250 - t333;
t351 = -t148 + t376;
t180 = Icges(5,5) * t239 + Icges(5,6) * t238 + Icges(5,3) * t269;
t181 = Icges(5,4) * t239 + Icges(5,2) * t238 + Icges(5,6) * t269;
t182 = Icges(5,1) * t239 + Icges(5,4) * t238 + Icges(5,5) * t269;
t96 = t269 * t180 + t238 * t181 + t239 * t182;
t175 = rSges(6,1) * t236 - rSges(6,2) * t235 + rSges(6,3) * t269;
t350 = -t175 + t370;
t226 = Icges(4,5) * t270 - Icges(4,6) * t269 + Icges(4,3) * t288;
t227 = Icges(4,4) * t270 - Icges(4,2) * t269 + Icges(4,6) * t288;
t228 = Icges(4,1) * t270 - Icges(4,4) * t269 + Icges(4,5) * t288;
t125 = t288 * t226 - t269 * t227 + t270 * t228;
t149 = t219 * rSges(6,1) - t218 * rSges(6,2) + t252 * rSges(6,3);
t193 = t253 * rSges(4,1) - t252 * rSges(4,2) + t272 * rSges(4,3);
t280 = Icges(3,3) * t308 + (Icges(3,5) * t312 + Icges(3,6) * t315) * t306;
t281 = Icges(3,6) * t308 + (Icges(3,4) * t312 + Icges(3,2) * t315) * t306;
t282 = Icges(3,5) * t308 + (Icges(3,1) * t312 + Icges(3,4) * t315) * t306;
t348 = t308 * t280 + t281 * t387 + t282 * t389;
t264 = t291 * rSges(3,1) + t290 * rSges(3,2) + rSges(3,3) * t388;
t347 = -t313 * pkin(1) + pkin(9) * t386;
t229 = rSges(4,1) * t270 - rSges(4,2) * t269 + rSges(4,3) * t288;
t275 = pkin(2) * t389 + pkin(10) * t288;
t342 = t306 * (-t229 - t275);
t341 = t376 - t380;
t340 = t370 - t374;
t339 = t203 * t388 + t204 * t386 + t367;
t336 = t306 * (-t275 + t369);
t331 = t306 * (-t275 + t350);
t329 = t71 / 0.2e1 + t81 / 0.2e1 + t354;
t328 = t72 / 0.2e1 + t82 / 0.2e1 + t353;
t327 = t136 * t388 + t137 * t386 + t339;
t325 = t306 * (-t275 + t340);
t192 = rSges(4,1) * t251 - rSges(4,2) * t250 + rSges(4,3) * t271;
t324 = -t255 + t347;
t323 = t256 + t366;
t320 = t323 + t349;
t186 = Icges(4,5) * t251 - Icges(4,6) * t250 + Icges(4,3) * t271;
t188 = Icges(4,4) * t251 - Icges(4,2) * t250 + Icges(4,6) * t271;
t190 = Icges(4,1) * t251 - Icges(4,4) * t250 + Icges(4,5) * t271;
t101 = t186 * t288 - t188 * t269 + t190 * t270;
t118 = t226 * t271 - t227 * t250 + t228 * t251;
t150 = Icges(5,5) * t223 + Icges(5,6) * t222 + Icges(5,3) * t250;
t152 = Icges(5,4) * t223 + Icges(5,2) * t222 + Icges(5,6) * t250;
t154 = Icges(5,1) * t223 + Icges(5,4) * t222 + Icges(5,5) * t250;
t73 = t150 * t269 + t152 * t238 + t154 * t239;
t83 = t180 * t250 + t181 * t222 + t182 * t223;
t319 = t101 / 0.2e1 + t118 / 0.2e1 + t83 / 0.2e1 + t73 / 0.2e1 + t329;
t187 = Icges(4,5) * t253 - Icges(4,6) * t252 + Icges(4,3) * t272;
t189 = Icges(4,4) * t253 - Icges(4,2) * t252 + Icges(4,6) * t272;
t191 = Icges(4,1) * t253 - Icges(4,4) * t252 + Icges(4,5) * t272;
t102 = t187 * t288 - t189 * t269 + t191 * t270;
t119 = t226 * t272 - t227 * t252 + t228 * t253;
t151 = Icges(5,5) * t225 + Icges(5,6) * t224 + Icges(5,3) * t252;
t153 = Icges(5,4) * t225 + Icges(5,2) * t224 + Icges(5,6) * t252;
t155 = Icges(5,1) * t225 + Icges(5,4) * t224 + Icges(5,5) * t252;
t74 = t151 * t269 + t153 * t238 + t155 * t239;
t84 = t180 * t252 + t181 * t224 + t182 * t225;
t318 = t102 / 0.2e1 + t119 / 0.2e1 + t84 / 0.2e1 + t74 / 0.2e1 + t328;
t317 = -t251 * t302 + t324 - t364;
t263 = t289 * rSges(3,1) + rSges(3,2) * t326 - rSges(3,3) * t386;
t297 = rSges(2,1) * t316 - t313 * rSges(2,2);
t296 = -t313 * rSges(2,1) - rSges(2,2) * t316;
t283 = rSges(3,3) * t308 + (rSges(3,1) * t312 + rSges(3,2) * t315) * t306;
t262 = Icges(3,1) * t291 + Icges(3,4) * t290 + Icges(3,5) * t388;
t261 = Icges(3,1) * t289 + Icges(3,4) * t326 - Icges(3,5) * t386;
t260 = Icges(3,4) * t291 + Icges(3,2) * t290 + Icges(3,6) * t388;
t259 = Icges(3,4) * t289 + Icges(3,2) * t326 - Icges(3,6) * t386;
t258 = Icges(3,5) * t291 + Icges(3,6) * t290 + Icges(3,3) * t388;
t257 = Icges(3,5) * t289 + Icges(3,6) * t326 - Icges(3,3) * t386;
t249 = t264 + t366;
t248 = -t263 + t347;
t232 = -t308 * t263 - t283 * t386;
t231 = t264 * t308 - t283 * t388;
t230 = t348 * t308;
t208 = (t263 * t313 + t264 * t316) * t306;
t207 = t280 * t388 + t281 * t290 + t282 * t291;
t206 = -t280 * t386 + t281 * t326 + t289 * t282;
t185 = t258 * t308 + (t260 * t315 + t262 * t312) * t306;
t184 = t257 * t308 + (t259 * t315 + t261 * t312) * t306;
t163 = t323 + t193;
t162 = -t192 + t324;
t132 = t193 * t288 - t229 * t272;
t131 = -t192 * t288 + t229 * t271;
t128 = (-t192 - t255) * t308 + t316 * t342;
t127 = t193 * t308 + t313 * t342 + t254;
t124 = t125 * t308;
t123 = t192 * t272 - t193 * t271;
t122 = t125 * t288;
t117 = (t192 * t313 + t193 * t316) * t306 + t367;
t116 = t323 - t372;
t115 = t324 + t373;
t106 = t320 + t149;
t105 = (-rSges(6,3) + t309) * t250 + t317 + t333;
t104 = t149 * t269 - t175 * t252;
t103 = -t148 * t269 + t175 * t250;
t100 = t187 * t272 - t189 * t252 + t191 * t253;
t99 = t186 * t272 - t188 * t252 + t190 * t253;
t98 = t187 * t271 - t189 * t250 + t191 * t251;
t97 = t186 * t271 - t188 * t250 + t190 * t251;
t95 = t96 * t308;
t94 = t148 * t252 - t149 * t250;
t93 = t96 * t288;
t88 = t157 * t288 + t272 * t369 + t196;
t87 = t183 * t271 + t288 * t373 + t215;
t86 = (-t255 + t373) * t308 + t316 * t336;
t85 = t157 * t308 + t313 * t336 + t368;
t80 = t320 + t379;
t79 = -t397 + t393 + (-rSges(7,3) - pkin(12)) * t216 + t317 + t332;
t78 = t156 * t272 + t271 * t372 + t194;
t77 = t114 * t235 - t141 * t218;
t76 = -t113 * t235 + t141 * t216;
t75 = (t156 * t313 + t157 * t316) * t306 + t339;
t70 = t151 * t252 + t153 * t224 + t155 * t225;
t69 = t150 * t252 + t152 * t224 + t154 * t225;
t68 = t151 * t250 + t153 * t222 + t155 * t223;
t67 = t150 * t250 + t152 * t222 + t154 * t223;
t62 = t113 * t218 - t114 * t216;
t56 = (-t255 + t351) * t308 + t316 * t331;
t55 = t149 * t308 + t313 * t331 + t352;
t54 = t149 * t288 + t272 * t350 + t377;
t53 = t175 * t271 + t288 * t351 + t371;
t52 = -t252 * t374 + t269 * t379;
t51 = t250 * t374 - t269 * t380;
t48 = t124 + (-t101 * t316 + t102 * t313) * t306;
t47 = (t148 * t313 + t149 * t316) * t306 + t327;
t46 = t148 * t272 + (-t149 + t375) * t271 + t378;
t45 = -t250 * t379 + t252 * t380;
t44 = t101 * t271 + t102 * t272 + t122;
t41 = t119 * t308 + (t100 * t313 - t316 * t99) * t306;
t40 = t118 * t308 + (t313 * t98 - t316 * t97) * t306;
t39 = t100 * t272 + t119 * t288 + t271 * t99;
t38 = t118 * t288 + t271 * t97 + t272 * t98;
t33 = (-t255 + t341) * t308 + t316 * t325;
t32 = t308 * t379 + t313 * t325 + t352;
t31 = t272 * t340 + t288 * t379 + t377;
t30 = t271 * t374 + t288 * t341 + t371;
t29 = (t313 * t380 + t316 * t379) * t306 + t327;
t28 = t380 * t272 + (t375 - t379) * t271 + t378;
t27 = t95 + (t74 * t313 - t73 * t316) * t306;
t26 = t73 * t271 + t74 * t272 + t93;
t22 = t84 * t308 + (t313 * t70 - t316 * t69) * t306;
t21 = t83 * t308 + (t313 * t68 - t316 * t67) * t306;
t20 = t271 * t69 + t272 * t70 + t288 * t84;
t19 = t271 * t67 + t272 * t68 + t288 * t83;
t120 = [t348 + Icges(2,3) + m(7) * (t79 ^ 2 + t80 ^ 2) + m(6) * (t105 ^ 2 + t106 ^ 2) + m(5) * (t115 ^ 2 + t116 ^ 2) + m(4) * (t162 ^ 2 + t163 ^ 2) + m(3) * (t248 ^ 2 + t249 ^ 2) + m(2) * (t296 ^ 2 + t297 ^ 2) + t125 + t96 + t92 + t61; t91 + t230 + t124 + t60 + t95 + m(4) * (t127 * t163 + t128 * t162) + m(5) * (t115 * t86 + t116 * t85) + m(3) * (t231 * t249 + t232 * t248) + m(6) * (t105 * t56 + t106 * t55) + m(7) * (t32 * t80 + t33 * t79) + ((-t184 / 0.2e1 - t206 / 0.2e1 - t319) * t316 + (t185 / 0.2e1 + t207 / 0.2e1 + t318) * t313) * t306; (t12 + t25 + t48 + t27 + t230) * t308 + m(7) * (t29 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t47 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t75 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t117 ^ 2 + t127 ^ 2 + t128 ^ 2) + m(3) * (t208 ^ 2 + t231 ^ 2 + t232 ^ 2) + ((t18 + t22 + t41 + t8 + (t258 * t388 + t260 * t290 + t262 * t291) * t388 + (t185 + t207) * t308) * t313 + (-t7 - t17 - t21 - t40 + (t289 * t261 * t306 + (-t257 * t306 + t259 * t385) * t386) * t316 + (-t259 * t290 - t261 * t291 - t289 * t262 + (-t259 * t312 - t260 * t385) * t316 + t260 * t383 + (-t257 * t313 + t258 * t316) * t306) * t388 + (-t206 - t184) * t308) * t316) * t306; t59 + t90 + t122 + t93 + m(7) * (t30 * t79 + t31 * t80) + m(6) * (t105 * t53 + t106 * t54) + m(5) * (t115 * t87 + t116 * t88) + m(4) * (t131 * t162 + t132 * t163) + t318 * t272 + t319 * t271; (t26 / 0.2e1 + t44 / 0.2e1 + t356) * t308 + (t27 / 0.2e1 + t48 / 0.2e1 + t355) * t288 + (t22 / 0.2e1 + t41 / 0.2e1 + t358) * t272 + (t21 / 0.2e1 + t40 / 0.2e1 + t359) * t271 + m(5) * (t75 * t78 + t85 * t88 + t86 * t87) + m(6) * (t46 * t47 + t53 * t56 + t54 * t55) + m(4) * (t117 * t123 + t127 * t132 + t128 * t131) + m(7) * (t28 * t29 + t30 * t33 + t31 * t32) + ((-t38 / 0.2e1 - t19 / 0.2e1 - t361) * t316 + (t39 / 0.2e1 + t20 / 0.2e1 + t360) * t313) * t306; (t24 + t26 + t44 + t11) * t288 + (t6 + t16 + t39 + t20) * t272 + (t5 + t15 + t38 + t19) * t271 + m(6) * (t46 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t78 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t123 ^ 2 + t131 ^ 2 + t132 ^ 2) + m(7) * (t28 ^ 2 + t30 ^ 2 + t31 ^ 2); m(7) * (t250 * t80 + t252 * t79) + m(6) * (t105 * t252 + t106 * t250) + m(5) * (t115 * t252 + t116 * t250); m(7) * (t250 * t32 + t252 * t33 + t269 * t29) + m(6) * (t250 * t55 + t252 * t56 + t269 * t47) + m(5) * (t250 * t85 + t252 * t86 + t269 * t75); m(7) * (t250 * t31 + t252 * t30 + t269 * t28) + m(6) * (t250 * t54 + t252 * t53 + t269 * t46) + m(5) * (t250 * t88 + t252 * t87 + t269 * t78); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t250 ^ 2 + t252 ^ 2 + t269 ^ 2); t58 + t89 + m(7) * (t51 * t79 + t52 * t80) + m(6) * (t103 * t105 + t104 * t106) + t328 * t252 + t329 * t250; t357 * t308 + t355 * t269 + t358 * t252 + t359 * t250 + m(7) * (t29 * t45 + t32 * t52 + t33 * t51) + m(6) * (t103 * t56 + t104 * t55 + t47 * t94) + (t313 * t362 - t316 * t363) * t306; t357 * t288 + t362 * t272 + t363 * t271 + t356 * t269 + t360 * t252 + t361 * t250 + m(7) * (t28 * t45 + t30 * t51 + t31 * t52) + m(6) * (t103 * t53 + t104 * t54 + t46 * t94); m(6) * (t103 * t252 + t104 * t250 + t269 * t94) + m(7) * (t250 * t52 + t252 * t51 + t269 * t45); (t10 + t23) * t269 + (t4 + t14) * t252 + (t3 + t13) * t250 + m(7) * (t45 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(6) * (t103 ^ 2 + t104 ^ 2 + t94 ^ 2); m(7) * (t76 * t79 + t77 * t80) + t57 + t353 * t218 + t354 * t216; m(7) * (t29 * t62 + t32 * t77 + t33 * t76) + t308 * t402 + t12 * t399 + t7 * t401 + t8 * t400 + (-t316 * t1 / 0.2e1 + t313 * t403) * t306; t11 * t399 + m(7) * (t28 * t62 + t30 * t76 + t31 * t77) + t5 * t401 + t288 * t402 + t271 * t404 + t272 * t403 + t6 * t400; m(7) * (t250 * t77 + t252 * t76 + t269 * t62); m(7) * (t45 * t62 + t51 * t76 + t52 * t77) + t4 * t400 + t3 * t401 + t269 * t402 + t10 * t399 + t252 * t403 + t250 * t404; t218 * t2 + t216 * t1 + t235 * t9 + m(7) * (t62 ^ 2 + t76 ^ 2 + t77 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t120(1) t120(2) t120(4) t120(7) t120(11) t120(16); t120(2) t120(3) t120(5) t120(8) t120(12) t120(17); t120(4) t120(5) t120(6) t120(9) t120(13) t120(18); t120(7) t120(8) t120(9) t120(10) t120(14) t120(19); t120(11) t120(12) t120(13) t120(14) t120(15) t120(20); t120(16) t120(17) t120(18) t120(19) t120(20) t120(21);];
Mq  = res;
