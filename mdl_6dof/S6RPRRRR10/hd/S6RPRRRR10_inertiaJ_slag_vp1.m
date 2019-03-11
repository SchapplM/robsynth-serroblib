% Calculate joint inertia matrix for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:50
% EndTime: 2019-03-09 07:28:07
% DurationCPUTime: 7.12s
% Computational Cost: add. (50430->587), mult. (119499->804), div. (0->0), fcn. (157184->16), ass. (0->291)
t289 = cos(pkin(6));
t296 = cos(qJ(1));
t360 = sin(pkin(13));
t319 = t296 * t360;
t293 = sin(qJ(1));
t362 = cos(pkin(13));
t322 = t293 * t362;
t273 = t289 * t319 + t322;
t292 = sin(qJ(3));
t320 = t296 * t362;
t321 = t293 * t360;
t307 = -t289 * t320 + t321;
t363 = cos(pkin(7));
t304 = t307 * t363;
t288 = sin(pkin(6));
t361 = sin(pkin(7));
t370 = cos(qJ(3));
t316 = t370 * t361;
t310 = t288 * t316;
t254 = t273 * t292 + t296 * t310 + t304 * t370;
t274 = -t289 * t321 + t320;
t306 = t289 * t322 + t319;
t302 = t306 * t363;
t256 = t274 * t292 - t293 * t310 + t302 * t370;
t314 = t363 * t362;
t323 = t288 * t360;
t262 = -t288 * t314 * t370 - t289 * t316 + t292 * t323;
t324 = t288 * t361;
t255 = t273 * t370 + (-t296 * t324 - t304) * t292;
t325 = t288 * t363;
t264 = -t296 * t325 + t307 * t361;
t353 = qJ(4) + qJ(5);
t285 = sin(t353);
t327 = cos(t353);
t226 = t255 * t327 + t264 * t285;
t290 = sin(qJ(6));
t294 = cos(qJ(6));
t190 = -t226 * t290 + t254 * t294;
t191 = t226 * t294 + t254 * t290;
t225 = t255 * t285 - t264 * t327;
t123 = Icges(7,5) * t191 + Icges(7,6) * t190 + Icges(7,3) * t225;
t125 = Icges(7,4) * t191 + Icges(7,2) * t190 + Icges(7,6) * t225;
t127 = Icges(7,1) * t191 + Icges(7,4) * t190 + Icges(7,5) * t225;
t52 = t123 * t225 + t125 * t190 + t127 * t191;
t257 = t274 * t370 + (t293 * t324 - t302) * t292;
t265 = t293 * t325 + t306 * t361;
t228 = t257 * t327 + t265 * t285;
t192 = -t228 * t290 + t256 * t294;
t193 = t228 * t294 + t256 * t290;
t227 = t257 * t285 - t265 * t327;
t124 = Icges(7,5) * t193 + Icges(7,6) * t192 + Icges(7,3) * t227;
t126 = Icges(7,4) * t193 + Icges(7,2) * t192 + Icges(7,6) * t227;
t128 = Icges(7,1) * t193 + Icges(7,4) * t192 + Icges(7,5) * t227;
t53 = t124 * t225 + t126 * t190 + t128 * t191;
t263 = t289 * t361 * t292 + (t292 * t314 + t360 * t370) * t288;
t272 = t289 * t363 - t324 * t362;
t242 = t263 * t327 + t272 * t285;
t218 = -t242 * t290 + t262 * t294;
t219 = t242 * t294 + t262 * t290;
t241 = t263 * t285 - t272 * t327;
t151 = Icges(7,5) * t219 + Icges(7,6) * t218 + Icges(7,3) * t241;
t152 = Icges(7,4) * t219 + Icges(7,2) * t218 + Icges(7,6) * t241;
t153 = Icges(7,1) * t219 + Icges(7,4) * t218 + Icges(7,5) * t241;
t65 = t151 * t225 + t152 * t190 + t153 * t191;
t11 = t254 * t52 + t256 * t53 + t262 * t65;
t155 = Icges(6,5) * t226 - Icges(6,6) * t225 + Icges(6,3) * t254;
t157 = Icges(6,4) * t226 - Icges(6,2) * t225 + Icges(6,6) * t254;
t159 = Icges(6,1) * t226 - Icges(6,4) * t225 + Icges(6,5) * t254;
t78 = t155 * t254 - t157 * t225 + t159 * t226;
t156 = Icges(6,5) * t228 - Icges(6,6) * t227 + Icges(6,3) * t256;
t158 = Icges(6,4) * t228 - Icges(6,2) * t227 + Icges(6,6) * t256;
t160 = Icges(6,1) * t228 - Icges(6,4) * t227 + Icges(6,5) * t256;
t79 = t156 * t254 - t158 * t225 + t160 * t226;
t194 = Icges(6,5) * t242 - Icges(6,6) * t241 + Icges(6,3) * t262;
t195 = Icges(6,4) * t242 - Icges(6,2) * t241 + Icges(6,6) * t262;
t196 = Icges(6,1) * t242 - Icges(6,4) * t241 + Icges(6,5) * t262;
t99 = t194 * t254 - t195 * t225 + t196 * t226;
t390 = t254 * t78 + t256 * t79 + t262 * t99 + t11;
t100 = t194 * t256 - t195 * t227 + t196 * t228;
t54 = t123 * t227 + t125 * t192 + t127 * t193;
t55 = t124 * t227 + t126 * t192 + t128 * t193;
t66 = t151 * t227 + t152 * t192 + t153 * t193;
t12 = t254 * t54 + t256 * t55 + t262 * t66;
t80 = t155 * t256 - t157 * t227 + t159 * t228;
t81 = t156 * t256 - t158 * t227 + t160 * t228;
t389 = t100 * t262 + t254 * t80 + t256 * t81 + t12;
t15 = t264 * t52 + t265 * t53 + t272 * t65;
t388 = t264 * t78 + t265 * t79 + t272 * t99 + t15;
t16 = t264 * t54 + t265 * t55 + t272 * t66;
t387 = t100 * t272 + t264 * t80 + t265 * t81 + t16;
t108 = t262 * t194 - t241 * t195 + t242 * t196;
t105 = t108 * t262;
t59 = t124 * t241 + t126 * t218 + t128 * t219;
t366 = t59 * t256;
t58 = t123 * t241 + t125 * t218 + t127 * t219;
t367 = t58 * t254;
t74 = t241 * t151 + t218 * t152 + t219 * t153;
t72 = t74 * t262;
t22 = t366 + t72 + t367;
t89 = t156 * t262 - t158 * t241 + t160 * t242;
t364 = t89 * t256;
t88 = t155 * t262 - t157 * t241 + t159 * t242;
t365 = t88 * t254;
t386 = t22 + t105 + t364 + t365;
t106 = t108 * t272;
t73 = t74 * t272;
t24 = t58 * t264 + t59 * t265 + t73;
t385 = t88 * t264 + t89 * t265 + t106 + t24;
t311 = -t191 * rSges(7,1) - t190 * rSges(7,2);
t131 = t225 * rSges(7,3) - t311;
t369 = t226 * pkin(5);
t350 = t225 * pkin(12) + t131 + t369;
t132 = t193 * rSges(7,1) + t192 * rSges(7,2) + t227 * rSges(7,3);
t349 = t228 * pkin(5) + t227 * pkin(12) + t132;
t154 = rSges(7,1) * t219 + rSges(7,2) * t218 + rSges(7,3) * t241;
t343 = pkin(5) * t242 + pkin(12) * t241 + t154;
t384 = m(3) / 0.2e1;
t383 = m(4) / 0.2e1;
t382 = m(5) / 0.2e1;
t381 = m(6) / 0.2e1;
t380 = m(7) / 0.2e1;
t379 = t225 / 0.2e1;
t378 = t227 / 0.2e1;
t377 = t241 / 0.2e1;
t376 = t254 / 0.2e1;
t375 = t256 / 0.2e1;
t374 = t262 / 0.2e1;
t373 = t264 / 0.2e1;
t372 = t265 / 0.2e1;
t371 = t272 / 0.2e1;
t295 = cos(qJ(4));
t284 = pkin(4) * t295 + pkin(3);
t368 = -pkin(3) + t284;
t297 = -pkin(11) - pkin(10);
t359 = t254 * t297;
t291 = sin(qJ(4));
t358 = t264 * t291;
t357 = t265 * t291;
t356 = t272 * t291;
t355 = t288 * t293;
t354 = t288 * t296;
t352 = t350 * t256;
t351 = t349 * t262;
t348 = t343 * t254;
t248 = t254 * pkin(10);
t337 = pkin(4) * t358;
t149 = t255 * t368 - t248 + t337 - t359;
t215 = t255 * pkin(3) + t248;
t210 = t265 * t215;
t347 = t265 * t149 + t210;
t216 = t257 * pkin(3) + t256 * pkin(10);
t328 = pkin(4) * t357 - t256 * t297 + t257 * t284;
t150 = -t216 + t328;
t211 = t272 * t216;
t346 = t272 * t150 + t211;
t312 = -t226 * rSges(6,1) + t225 * rSges(6,2);
t161 = t254 * rSges(6,3) - t312;
t345 = -t149 - t161;
t162 = t228 * rSges(6,1) - t227 * rSges(6,2) + t256 * rSges(6,3);
t344 = -t150 - t162;
t231 = -t255 * t291 + t264 * t295;
t232 = t255 * t295 + t358;
t170 = rSges(5,1) * t232 + rSges(5,2) * t231 + rSges(5,3) * t254;
t342 = -t170 - t215;
t233 = -t257 * t291 + t265 * t295;
t234 = t257 * t295 + t357;
t171 = t234 * rSges(5,1) + t233 * rSges(5,2) + t256 * rSges(5,3);
t341 = -t171 - t216;
t187 = pkin(4) * t356 + t368 * t263 + (-pkin(10) - t297) * t262;
t239 = t263 * pkin(3) + t262 * pkin(10);
t224 = t264 * t239;
t340 = t264 * t187 + t224;
t197 = rSges(6,1) * t242 - rSges(6,2) * t241 + rSges(6,3) * t262;
t339 = -t187 - t197;
t338 = t296 * pkin(1) + qJ(2) * t355;
t336 = t58 / 0.2e1 + t65 / 0.2e1;
t335 = t59 / 0.2e1 + t66 / 0.2e1;
t252 = -t263 * t291 + t272 * t295;
t253 = t263 * t295 + t356;
t198 = Icges(5,5) * t253 + Icges(5,6) * t252 + Icges(5,3) * t262;
t199 = Icges(5,4) * t253 + Icges(5,2) * t252 + Icges(5,6) * t262;
t200 = Icges(5,1) * t253 + Icges(5,4) * t252 + Icges(5,5) * t262;
t101 = t198 * t254 + t199 * t231 + t200 * t232;
t164 = Icges(5,5) * t232 + Icges(5,6) * t231 + Icges(5,3) * t254;
t166 = Icges(5,4) * t232 + Icges(5,2) * t231 + Icges(5,6) * t254;
t168 = Icges(5,1) * t232 + Icges(5,4) * t231 + Icges(5,5) * t254;
t90 = t164 * t262 + t166 * t252 + t168 * t253;
t334 = t101 / 0.2e1 + t90 / 0.2e1;
t102 = t198 * t256 + t199 * t233 + t200 * t234;
t165 = Icges(5,5) * t234 + Icges(5,6) * t233 + Icges(5,3) * t256;
t167 = Icges(5,4) * t234 + Icges(5,2) * t233 + Icges(5,6) * t256;
t169 = Icges(5,1) * t234 + Icges(5,4) * t233 + Icges(5,5) * t256;
t91 = t165 * t262 + t167 * t252 + t169 * t253;
t333 = t102 / 0.2e1 + t91 / 0.2e1;
t332 = -t149 - t350;
t331 = -t150 - t349;
t330 = -t187 - t343;
t112 = t262 * t198 + t252 * t199 + t253 * t200;
t235 = Icges(4,5) * t263 - Icges(4,6) * t262 + Icges(4,3) * t272;
t236 = Icges(4,4) * t263 - Icges(4,2) * t262 + Icges(4,6) * t272;
t237 = Icges(4,1) * t263 - Icges(4,4) * t262 + Icges(4,5) * t272;
t329 = t272 * t235 - t262 * t236 + t263 * t237;
t209 = t257 * rSges(4,1) - t256 * rSges(4,2) + t265 * rSges(4,3);
t326 = -t293 * pkin(1) + qJ(2) * t354;
t71 = t74 * t241;
t18 = t58 * t225 + t59 * t227 + t71;
t3 = t225 * t52 + t227 * t53 + t241 * t65;
t4 = t225 * t54 + t227 * t55 + t241 * t66;
t315 = t11 * t379 + t12 * t378 + t18 * t374 + t22 * t377 + t3 * t376 + t4 * t375;
t313 = t254 * t390 + t389 * t256 + t386 * t262;
t208 = rSges(4,1) * t255 - rSges(4,2) * t254 + rSges(4,3) * t264;
t309 = -t273 * pkin(2) - pkin(9) * t264 + t326;
t308 = t105 + t367 / 0.2e1 + t366 / 0.2e1 + t72 + t365 / 0.2e1 + t364 / 0.2e1 + (t65 + t99) * t376 + (t66 + t100) * t375;
t305 = t386 * t371 + t389 * t372 + t373 * t390 + t385 * t374 + t387 * t375 + t388 * t376;
t300 = -t255 * t284 + t309 - t337;
t299 = t274 * pkin(2) + pkin(9) * t265 + t338;
t298 = t299 + t328;
t280 = rSges(2,1) * t296 - rSges(2,2) * t293;
t279 = -rSges(2,1) * t293 - rSges(2,2) * t296;
t251 = t274 * rSges(3,1) - rSges(3,2) * t306 + rSges(3,3) * t355 + t338;
t250 = -t273 * rSges(3,1) + rSges(3,2) * t307 + rSges(3,3) * t354 + t326;
t238 = rSges(4,1) * t263 - rSges(4,2) * t262 + rSges(4,3) * t272;
t207 = Icges(4,1) * t257 - Icges(4,4) * t256 + Icges(4,5) * t265;
t206 = Icges(4,1) * t255 - Icges(4,4) * t254 + Icges(4,5) * t264;
t205 = Icges(4,4) * t257 - Icges(4,2) * t256 + Icges(4,6) * t265;
t204 = Icges(4,4) * t255 - Icges(4,2) * t254 + Icges(4,6) * t264;
t203 = Icges(4,5) * t257 - Icges(4,6) * t256 + Icges(4,3) * t265;
t202 = Icges(4,5) * t255 - Icges(4,6) * t254 + Icges(4,3) * t264;
t201 = rSges(5,1) * t253 + rSges(5,2) * t252 + rSges(5,3) * t262;
t180 = t299 + t209;
t179 = -t208 + t309;
t176 = t254 * t197;
t175 = t254 * t187;
t148 = t209 * t272 - t238 * t265;
t147 = -t208 * t272 + t238 * t264;
t145 = t262 * t162;
t143 = t262 * t150;
t142 = t256 * t161;
t140 = t256 * t149;
t138 = t208 * t265 - t209 * t264;
t137 = t329 * t272;
t134 = t235 * t265 - t236 * t256 + t237 * t257;
t133 = t235 * t264 - t236 * t254 + t237 * t255;
t130 = t299 - t341;
t129 = t309 + t342;
t121 = t298 + t162;
t120 = (-rSges(6,3) + t297) * t254 + t300 + t312;
t118 = t171 * t262 - t201 * t256;
t117 = -t170 * t262 + t201 * t254;
t116 = -t197 * t256 + t145;
t115 = -t161 * t262 + t176;
t114 = t203 * t272 - t205 * t262 + t207 * t263;
t113 = t202 * t272 - t204 * t262 + t206 * t263;
t111 = t170 * t256 - t171 * t254;
t110 = t112 * t272;
t109 = -t162 * t254 + t142;
t107 = t112 * t262;
t104 = t171 * t272 + t211 + (-t201 - t239) * t265;
t103 = t201 * t264 + t272 * t342 + t224;
t96 = t298 + t349;
t95 = -t369 + t359 + (-rSges(7,3) - pkin(12)) * t225 + t300 + t311;
t94 = t132 * t241 - t154 * t227;
t93 = -t131 * t241 + t154 * t225;
t92 = t170 * t265 + t264 * t341 + t210;
t85 = t165 * t256 + t167 * t233 + t169 * t234;
t84 = t164 * t256 + t166 * t233 + t168 * t234;
t83 = t165 * t254 + t167 * t231 + t169 * t232;
t82 = t164 * t254 + t166 * t231 + t168 * t232;
t77 = t256 * t339 + t143 + t145;
t76 = t262 * t345 + t175 + t176;
t75 = t131 * t227 - t132 * t225;
t70 = t162 * t272 + (-t239 + t339) * t265 + t346;
t69 = t197 * t264 + (-t215 + t345) * t272 + t340;
t68 = -t256 * t343 + t351;
t67 = -t262 * t350 + t348;
t64 = t254 * t344 + t140 + t142;
t61 = -t254 * t349 + t352;
t60 = t161 * t265 + (-t216 + t344) * t264 + t347;
t51 = t256 * t330 + t143 + t351;
t50 = t262 * t332 + t175 + t348;
t49 = t349 * t272 + (-t239 + t330) * t265 + t346;
t48 = t343 * t264 + (-t215 + t332) * t272 + t340;
t47 = t254 * t331 + t140 + t352;
t46 = t350 * t265 + (-t216 + t331) * t264 + t347;
t45 = t90 * t264 + t91 * t265 + t110;
t44 = t90 * t254 + t91 * t256 + t107;
t38 = t102 * t272 + t264 * t84 + t265 * t85;
t37 = t101 * t272 + t264 * t82 + t265 * t83;
t36 = t102 * t262 + t254 * t84 + t256 * t85;
t35 = t101 * t262 + t254 * t82 + t256 * t83;
t1 = [(Icges(3,5) * t289 + (Icges(3,1) * t360 + Icges(3,4) * t362) * t288) * t323 + t329 + m(7) * (t95 ^ 2 + t96 ^ 2) + m(6) * (t120 ^ 2 + t121 ^ 2) + m(5) * (t129 ^ 2 + t130 ^ 2) + m(4) * (t179 ^ 2 + t180 ^ 2) + m(3) * (t250 ^ 2 + t251 ^ 2) + m(2) * (t279 ^ 2 + t280 ^ 2) + t289 * (Icges(3,3) * t289 + (Icges(3,5) * t360 + Icges(3,6) * t362) * t288) + t288 * t362 * (Icges(3,6) * t289 + (Icges(3,4) * t360 + Icges(3,2) * t362) * t288) + Icges(2,3) + t112 + t108 + t74; 0.2e1 * ((t293 * t95 - t296 * t96) * t380 + (t120 * t293 - t121 * t296) * t381 + (t129 * t293 - t130 * t296) * t382 + (t179 * t293 - t180 * t296) * t383 + (t250 * t293 - t251 * t296) * t384) * t288; 0.2e1 * (t384 + t383 + t382 + t381 + t380) * (t289 ^ 2 + (t293 ^ 2 + t296 ^ 2) * t288 ^ 2); t73 + t106 + t110 + t137 + m(7) * (t48 * t95 + t49 * t96) + m(6) * (t120 * t69 + t121 * t70) + m(5) * (t103 * t129 + t104 * t130) + m(4) * (t147 * t179 + t148 * t180) + (t89 / 0.2e1 + t100 / 0.2e1 + t134 / 0.2e1 + t114 / 0.2e1 + t333 + t335) * t265 + (t88 / 0.2e1 + t99 / 0.2e1 + t133 / 0.2e1 + t113 / 0.2e1 + t334 + t336) * t264; m(4) * (t138 * t289 + (t147 * t293 - t148 * t296) * t288) + m(5) * (t289 * t92 + (t103 * t293 - t104 * t296) * t288) + m(6) * (t289 * t60 + (t293 * t69 - t296 * t70) * t288) + m(7) * (t289 * t46 + (t293 * t48 - t296 * t49) * t288); (t137 + t45 + t385) * t272 + (t38 + (t203 * t265 - t205 * t256 + t207 * t257) * t265 + (t114 + t134) * t272 + t387) * t265 + (t37 + (t202 * t264 - t204 * t254 + t206 * t255) * t264 + (t113 + t133) * t272 + (t202 * t265 + t203 * t264 - t204 * t256 - t205 * t254 + t206 * t257 + t207 * t255) * t265 + t388) * t264 + m(7) * (t46 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t60 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t103 ^ 2 + t104 ^ 2 + t92 ^ 2) + m(4) * (t138 ^ 2 + t147 ^ 2 + t148 ^ 2); t107 + t308 + m(6) * (t120 * t76 + t121 * t77) + m(5) * (t117 * t129 + t118 * t130) + t333 * t256 + m(7) * (t50 * t95 + t51 * t96) + t334 * t254; m(5) * (t111 * t289 + (t117 * t293 - t118 * t296) * t288) + m(6) * (t289 * t64 + (t293 * t76 - t296 * t77) * t288) + m(7) * (t289 * t47 + (t293 * t50 - t296 * t51) * t288); t44 * t371 + t36 * t372 + t35 * t373 + t45 * t374 + t38 * t375 + t305 + m(7) * (t46 * t47 + t48 * t50 + t49 * t51) + m(6) * (t64 * t60 + t69 * t76 + t70 * t77) + m(5) * (t103 * t117 + t104 * t118 + t111 * t92) + t37 * t376; t254 * t35 + t256 * t36 + t262 * t44 + m(7) * (t47 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t64 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t111 ^ 2 + t117 ^ 2 + t118 ^ 2) + t313; m(7) * (t67 * t95 + t68 * t96) + m(6) * (t115 * t120 + t116 * t121) + t308; m(6) * (t109 * t289 + (t115 * t293 - t116 * t296) * t288) + m(7) * (t289 * t61 + (t293 * t67 - t296 * t68) * t288); m(7) * (t61 * t46 + t48 * t67 + t49 * t68) + m(6) * (t109 * t60 + t115 * t69 + t116 * t70) + t305; m(7) * (t61 * t47 + t50 * t67 + t51 * t68) + m(6) * (t109 * t64 + t115 * t76 + t116 * t77) + t313; m(7) * (t61 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t109 ^ 2 + t115 ^ 2 + t116 ^ 2) + t313; m(7) * (t93 * t95 + t94 * t96) + t71 + t335 * t227 + t336 * t225; m(7) * (t289 * t75 + (t293 * t93 - t296 * t94) * t288); m(7) * (t46 * t75 + t48 * t93 + t49 * t94) + t24 * t377 + t4 * t372 + t16 * t378 + t15 * t379 + t18 * t371 + t3 * t373; m(7) * (t47 * t75 + t50 * t93 + t51 * t94) + t315; m(7) * (t61 * t75 + t67 * t93 + t68 * t94) + t315; t227 * t4 + t225 * t3 + t241 * t18 + m(7) * (t75 ^ 2 + t93 ^ 2 + t94 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
