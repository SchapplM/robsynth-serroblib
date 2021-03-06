% Calculate joint inertia matrix for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:50
% EndTime: 2019-03-10 05:17:13
% DurationCPUTime: 10.54s
% Computational Cost: add. (84348->783), mult. (220218->1052), div. (0->0), fcn. (290820->16), ass. (0->371)
t357 = cos(pkin(6));
t362 = sin(qJ(1));
t364 = cos(qJ(2));
t440 = t362 * t364;
t361 = sin(qJ(2));
t365 = cos(qJ(1));
t441 = t361 * t365;
t339 = -t357 * t440 - t441;
t356 = sin(pkin(6));
t452 = cos(pkin(7));
t408 = t356 * t452;
t451 = sin(pkin(7));
t320 = -t339 * t451 + t362 * t408;
t439 = t364 * t365;
t442 = t361 * t362;
t337 = t357 * t439 - t442;
t319 = -t337 * t451 - t365 * t408;
t443 = t356 * t365;
t338 = t357 * t441 + t440;
t360 = sin(qJ(3));
t407 = t356 * t451;
t457 = cos(qJ(3));
t300 = t338 * t457 + (t337 * t452 - t365 * t407) * t360;
t359 = sin(qJ(4));
t456 = cos(qJ(4));
t275 = t300 * t359 - t319 * t456;
t467 = t275 / 0.2e1;
t340 = -t357 * t442 + t439;
t302 = t340 * t457 + (t339 * t452 + t362 * t407) * t360;
t277 = t302 * t359 - t320 * t456;
t466 = t277 / 0.2e1;
t318 = t357 * t451 * t360 + (t360 * t364 * t452 + t361 * t457) * t356;
t336 = t357 * t452 - t364 * t407;
t293 = t318 * t359 - t336 * t456;
t465 = t293 / 0.2e1;
t383 = t457 * t451;
t377 = t356 * t383;
t384 = t452 * t457;
t299 = -t337 * t384 + t338 * t360 + t365 * t377;
t464 = t299 / 0.2e1;
t301 = -t339 * t384 + t340 * t360 - t362 * t377;
t463 = t301 / 0.2e1;
t444 = t356 * t364;
t446 = t356 * t361;
t317 = -t357 * t383 + t360 * t446 - t384 * t444;
t462 = t317 / 0.2e1;
t461 = t319 / 0.2e1;
t460 = t320 / 0.2e1;
t459 = t336 / 0.2e1;
t458 = t357 / 0.2e1;
t363 = cos(qJ(5));
t351 = pkin(5) * t363 + pkin(4);
t455 = -pkin(4) + t351;
t276 = t300 * t456 + t319 * t359;
t355 = qJ(5) + qJ(6);
t352 = sin(t355);
t353 = cos(t355);
t223 = -t276 * t352 + t299 * t353;
t224 = t276 * t353 + t299 * t352;
t156 = Icges(7,5) * t224 + Icges(7,6) * t223 + Icges(7,3) * t275;
t158 = Icges(7,4) * t224 + Icges(7,2) * t223 + Icges(7,6) * t275;
t160 = Icges(7,1) * t224 + Icges(7,4) * t223 + Icges(7,5) * t275;
t294 = t318 * t456 + t336 * t359;
t262 = -t294 * t352 + t317 * t353;
t263 = t294 * t353 + t317 * t352;
t83 = t156 * t293 + t158 * t262 + t160 * t263;
t454 = t83 * t275;
t278 = t302 * t456 + t320 * t359;
t225 = -t278 * t352 + t301 * t353;
t226 = t278 * t353 + t301 * t352;
t157 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t277;
t159 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t277;
t161 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t277;
t84 = t157 * t293 + t159 * t262 + t161 * t263;
t453 = t84 * t277;
t358 = sin(qJ(5));
t450 = t299 * t358;
t449 = t301 * t358;
t306 = Icges(3,5) * t338 + Icges(3,6) * t337 - Icges(3,3) * t443;
t448 = t306 * t365;
t447 = t317 * t358;
t445 = t356 * t362;
t273 = t275 * pkin(12);
t366 = -pkin(13) - pkin(12);
t422 = pkin(5) * t450;
t151 = -t275 * t366 + t276 * t455 - t273 + t422;
t379 = -t224 * rSges(7,1) - t223 * rSges(7,2);
t162 = t275 * rSges(7,3) - t379;
t438 = t151 + t162;
t218 = t278 * pkin(4) + t277 * pkin(12);
t413 = pkin(5) * t449 - t277 * t366 + t278 * t351;
t152 = -t218 + t413;
t163 = t226 * rSges(7,1) + t225 * rSges(7,2) + t277 * rSges(7,3);
t437 = t152 + t163;
t229 = -t276 * t358 + t299 * t363;
t230 = t276 * t363 + t450;
t170 = rSges(6,1) * t230 + rSges(6,2) * t229 + rSges(6,3) * t275;
t217 = t276 * pkin(4) + t273;
t436 = -t170 - t217;
t231 = -t278 * t358 + t301 * t363;
t232 = t278 * t363 + t449;
t171 = t232 * rSges(6,1) + t231 * rSges(6,2) + t277 * rSges(6,3);
t435 = -t171 - t218;
t189 = pkin(5) * t447 + t455 * t294 + (-pkin(12) - t366) * t293;
t193 = rSges(7,1) * t263 + rSges(7,2) * t262 + rSges(7,3) * t293;
t434 = t189 + t193;
t267 = -t294 * t358 + t317 * t363;
t268 = t294 * t363 + t447;
t197 = rSges(6,1) * t268 + rSges(6,2) * t267 + rSges(6,3) * t293;
t255 = t294 * pkin(4) + t293 * pkin(12);
t433 = -t197 - t255;
t204 = rSges(5,1) * t276 - rSges(5,2) * t275 + rSges(5,3) * t299;
t256 = t300 * pkin(3) + t299 * pkin(11);
t432 = -t204 - t256;
t248 = t320 * t256;
t431 = t320 * t217 + t248;
t257 = t302 * pkin(3) + t301 * pkin(11);
t249 = t336 * t257;
t430 = t336 * t218 + t249;
t238 = rSges(5,1) * t294 - rSges(5,2) * t293 + rSges(5,3) * t317;
t286 = pkin(3) * t318 + pkin(11) * t317;
t429 = -t238 - t286;
t265 = t319 * t286;
t428 = t319 * t255 + t265;
t305 = t340 * pkin(2) + t320 * pkin(10);
t303 = t357 * t305;
t427 = t357 * t257 + t303;
t304 = t338 * pkin(2) + t319 * pkin(10);
t426 = t304 * t445 + t305 * t443;
t425 = t365 * pkin(1) + pkin(9) * t445;
t69 = t156 * t277 + t158 * t225 + t160 * t226;
t70 = t157 * t277 + t159 * t225 + t161 * t226;
t190 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t293;
t191 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t293;
t192 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t293;
t97 = t190 * t277 + t191 * t225 + t192 * t226;
t10 = t275 * t69 + t277 * t70 + t293 * t97;
t108 = t293 * t190 + t262 * t191 + t263 * t192;
t102 = t108 * t293;
t35 = t102 + t453 + t454;
t67 = t156 * t275 + t158 * t223 + t160 * t224;
t68 = t157 * t275 + t159 * t223 + t161 * t224;
t96 = t190 * t275 + t191 * t223 + t192 * t224;
t9 = t275 * t67 + t277 * t68 + t293 * t96;
t423 = t277 * t10 + t275 * t9 + t293 * t35;
t165 = Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t277;
t167 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t277;
t169 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t277;
t86 = t165 * t293 + t167 * t267 + t169 * t268;
t194 = Icges(6,5) * t268 + Icges(6,6) * t267 + Icges(6,3) * t293;
t195 = Icges(6,4) * t268 + Icges(6,2) * t267 + Icges(6,6) * t293;
t196 = Icges(6,1) * t268 + Icges(6,4) * t267 + Icges(6,5) * t293;
t99 = t194 * t277 + t195 * t231 + t196 * t232;
t421 = t86 / 0.2e1 + t99 / 0.2e1;
t164 = Icges(6,5) * t230 + Icges(6,6) * t229 + Icges(6,3) * t275;
t166 = Icges(6,4) * t230 + Icges(6,2) * t229 + Icges(6,6) * t275;
t168 = Icges(6,1) * t230 + Icges(6,4) * t229 + Icges(6,5) * t275;
t85 = t164 * t293 + t166 * t267 + t168 * t268;
t98 = t194 * t275 + t195 * t229 + t196 * t230;
t420 = t98 / 0.2e1 + t85 / 0.2e1;
t419 = -t217 - t438;
t418 = -t218 - t437;
t417 = -t256 + t436;
t112 = t293 * t194 + t267 * t195 + t268 * t196;
t416 = -t255 - t434;
t415 = -t286 + t433;
t233 = Icges(5,5) * t294 - Icges(5,6) * t293 + Icges(5,3) * t317;
t234 = Icges(5,4) * t294 - Icges(5,2) * t293 + Icges(5,6) * t317;
t235 = Icges(5,1) * t294 - Icges(5,4) * t293 + Icges(5,5) * t317;
t140 = t317 * t233 - t293 * t234 + t294 * t235;
t414 = t357 * t218 + t427;
t279 = Icges(4,5) * t318 - Icges(4,6) * t317 + Icges(4,3) * t336;
t280 = Icges(4,4) * t318 - Icges(4,2) * t317 + Icges(4,6) * t336;
t281 = Icges(4,1) * t318 - Icges(4,4) * t317 + Icges(4,5) * t336;
t182 = t336 * t279 - t317 * t280 + t318 * t281;
t205 = t278 * rSges(5,1) - t277 * rSges(5,2) + t301 * rSges(5,3);
t247 = t302 * rSges(4,1) - t301 * rSges(4,2) + t320 * rSges(4,3);
t327 = Icges(3,3) * t357 + (Icges(3,5) * t361 + Icges(3,6) * t364) * t356;
t328 = Icges(3,6) * t357 + (Icges(3,4) * t361 + Icges(3,2) * t364) * t356;
t329 = Icges(3,5) * t357 + (Icges(3,1) * t361 + Icges(3,4) * t364) * t356;
t412 = t357 * t327 + t328 * t444 + t329 * t446;
t313 = t340 * rSges(3,1) + t339 * rSges(3,2) + rSges(3,3) * t445;
t411 = -t362 * pkin(1) + pkin(9) * t443;
t282 = rSges(4,1) * t318 - rSges(4,2) * t317 + rSges(4,3) * t336;
t323 = pkin(2) * t446 + pkin(10) * t336;
t401 = t356 * (-t282 - t323);
t400 = -t256 + t419;
t399 = -t286 + t416;
t398 = t256 * t445 + t257 * t443 + t426;
t17 = t299 * t67 + t301 * t68 + t317 * t96;
t73 = t164 * t275 + t166 * t229 + t168 * t230;
t74 = t165 * t275 + t167 * t229 + t169 * t230;
t23 = t299 * t73 + t301 * t74 + t317 * t98;
t198 = Icges(5,5) * t276 - Icges(5,6) * t275 + Icges(5,3) * t299;
t200 = Icges(5,4) * t276 - Icges(5,2) * t275 + Icges(5,6) * t299;
t202 = Icges(5,1) * t276 - Icges(5,4) * t275 + Icges(5,5) * t299;
t114 = t198 * t299 - t200 * t275 + t202 * t276;
t199 = Icges(5,5) * t278 - Icges(5,6) * t277 + Icges(5,3) * t301;
t201 = Icges(5,4) * t278 - Icges(5,2) * t277 + Icges(5,6) * t301;
t203 = Icges(5,1) * t278 - Icges(5,4) * t277 + Icges(5,5) * t301;
t115 = t199 * t299 - t201 * t275 + t203 * t276;
t130 = t233 * t299 - t234 * t275 + t235 * t276;
t46 = t114 * t299 + t115 * t301 + t130 * t317;
t397 = t17 / 0.2e1 + t23 / 0.2e1 + t46 / 0.2e1;
t18 = t299 * t69 + t301 * t70 + t317 * t97;
t75 = t164 * t277 + t166 * t231 + t168 * t232;
t76 = t165 * t277 + t167 * t231 + t169 * t232;
t24 = t299 * t75 + t301 * t76 + t317 * t99;
t116 = t198 * t301 - t200 * t277 + t202 * t278;
t117 = t199 * t301 - t201 * t277 + t203 * t278;
t131 = t233 * t301 - t234 * t277 + t235 * t278;
t47 = t116 * t299 + t117 * t301 + t131 * t317;
t396 = t18 / 0.2e1 + t24 / 0.2e1 + t47 / 0.2e1;
t21 = t319 * t67 + t320 * t68 + t336 * t96;
t27 = t319 * t73 + t320 * t74 + t336 * t98;
t48 = t114 * t319 + t115 * t320 + t130 * t336;
t395 = t21 / 0.2e1 + t27 / 0.2e1 + t48 / 0.2e1;
t22 = t319 * t69 + t320 * t70 + t336 * t97;
t28 = t319 * t75 + t320 * t76 + t336 * t99;
t49 = t116 * t319 + t117 * t320 + t131 * t336;
t394 = t22 / 0.2e1 + t28 / 0.2e1 + t49 / 0.2e1;
t25 = t357 * t96 + (t362 * t68 - t365 * t67) * t356;
t29 = t357 * t98 + (t362 * t74 - t365 * t73) * t356;
t50 = t130 * t357 + (-t114 * t365 + t115 * t362) * t356;
t393 = t25 / 0.2e1 + t29 / 0.2e1 + t50 / 0.2e1;
t26 = t357 * t97 + (t362 * t70 - t365 * t69) * t356;
t30 = t357 * t99 + (t362 * t76 - t365 * t75) * t356;
t51 = t131 * t357 + (-t116 * t365 + t117 * t362) * t356;
t392 = t26 / 0.2e1 + t30 / 0.2e1 + t51 / 0.2e1;
t103 = t108 * t317;
t38 = t83 * t299 + t84 * t301 + t103;
t107 = t112 * t317;
t42 = t85 * t299 + t86 * t301 + t107;
t118 = t198 * t317 - t200 * t293 + t202 * t294;
t119 = t199 * t317 - t201 * t293 + t203 * t294;
t136 = t140 * t317;
t53 = t118 * t299 + t119 * t301 + t136;
t391 = t38 / 0.2e1 + t42 / 0.2e1 + t53 / 0.2e1;
t105 = t108 * t336;
t40 = t83 * t319 + t84 * t320 + t105;
t109 = t112 * t336;
t44 = t85 * t319 + t86 * t320 + t109;
t137 = t140 * t336;
t55 = t118 * t319 + t119 * t320 + t137;
t390 = t40 / 0.2e1 + t44 / 0.2e1 + t55 / 0.2e1;
t106 = t108 * t357;
t43 = t106 + (t84 * t362 - t83 * t365) * t356;
t111 = t112 * t357;
t45 = t111 + (t86 * t362 - t85 * t365) * t356;
t139 = t140 * t357;
t57 = t139 + (-t118 * t365 + t119 * t362) * t356;
t389 = t43 / 0.2e1 + t45 / 0.2e1 + t57 / 0.2e1;
t386 = t102 + t454 / 0.2e1 + t453 / 0.2e1 + t96 * t467 + t97 * t466;
t385 = t356 * (-t323 + t429);
t382 = t10 * t463 + t17 * t467 + t18 * t466 + t35 * t462 + t38 * t465 + t9 * t464;
t381 = t10 * t460 + t21 * t467 + t22 * t466 + t35 * t459 + t40 * t465 + t9 * t461;
t380 = t25 * t467 + t26 * t466 + t35 * t458 + t43 * t465 + t10 * t445 / 0.2e1 - t9 * t443 / 0.2e1;
t378 = t356 * (-t323 + t415);
t376 = t217 * t445 + t218 * t443 + t398;
t375 = t356 * (-t323 + t399);
t246 = rSges(4,1) * t300 - rSges(4,2) * t299 + rSges(4,3) * t319;
t374 = -t304 + t411;
t312 = rSges(3,1) * t338 + rSges(3,2) * t337 - rSges(3,3) * t443;
t373 = t305 + t425;
t372 = t83 / 0.2e1 + t96 / 0.2e1 + t130 / 0.2e1 + t118 / 0.2e1 + t420;
t371 = t84 / 0.2e1 + t97 / 0.2e1 + t119 / 0.2e1 + t131 / 0.2e1 + t421;
t370 = -t256 + t374;
t369 = t257 + t373;
t240 = Icges(4,5) * t300 - Icges(4,6) * t299 + Icges(4,3) * t319;
t242 = Icges(4,4) * t300 - Icges(4,2) * t299 + Icges(4,6) * t319;
t244 = Icges(4,1) * t300 - Icges(4,4) * t299 + Icges(4,5) * t319;
t145 = t240 * t336 - t242 * t317 + t244 * t318;
t172 = t279 * t319 - t280 * t299 + t281 * t300;
t368 = t172 / 0.2e1 + t145 / 0.2e1 + t372;
t241 = Icges(4,5) * t302 - Icges(4,6) * t301 + Icges(4,3) * t320;
t243 = Icges(4,4) * t302 - Icges(4,2) * t301 + Icges(4,6) * t320;
t245 = Icges(4,1) * t302 - Icges(4,4) * t301 + Icges(4,5) * t320;
t146 = t241 * t336 - t243 * t317 + t245 * t318;
t173 = t279 * t320 - t280 * t301 + t281 * t302;
t367 = t173 / 0.2e1 + t146 / 0.2e1 + t371;
t346 = rSges(2,1) * t365 - rSges(2,2) * t362;
t345 = -rSges(2,1) * t362 - rSges(2,2) * t365;
t330 = rSges(3,3) * t357 + (rSges(3,1) * t361 + rSges(3,2) * t364) * t356;
t311 = Icges(3,1) * t340 + Icges(3,4) * t339 + Icges(3,5) * t445;
t310 = Icges(3,1) * t338 + Icges(3,4) * t337 - Icges(3,5) * t443;
t309 = Icges(3,4) * t340 + Icges(3,2) * t339 + Icges(3,6) * t445;
t308 = Icges(3,4) * t338 + Icges(3,2) * t337 - Icges(3,6) * t443;
t307 = Icges(3,5) * t340 + Icges(3,6) * t339 + Icges(3,3) * t445;
t298 = t313 + t425;
t297 = -t312 + t411;
t285 = -t312 * t357 - t330 * t443;
t284 = t313 * t357 - t330 * t445;
t283 = t412 * t357;
t261 = (t312 * t362 + t313 * t365) * t356;
t260 = t327 * t445 + t328 * t339 + t329 * t340;
t259 = -t327 * t443 + t328 * t337 + t329 * t338;
t237 = t307 * t357 + (t309 * t364 + t311 * t361) * t356;
t236 = t306 * t357 + (t308 * t364 + t310 * t361) * t356;
t220 = t299 * t255;
t213 = t373 + t247;
t212 = -t246 + t374;
t207 = t317 * t218;
t206 = t301 * t217;
t188 = t247 * t336 - t282 * t320;
t187 = -t246 * t336 + t282 * t319;
t184 = (-t246 - t304) * t357 + t365 * t401;
t183 = t247 * t357 + t362 * t401 + t303;
t181 = t182 * t357;
t178 = t246 * t320 - t247 * t319;
t177 = t275 * t193;
t174 = t182 * t336;
t155 = (t246 * t362 + t247 * t365) * t356 + t426;
t154 = t369 + t205;
t153 = -t204 + t370;
t150 = t293 * t163;
t149 = t205 * t317 - t238 * t301;
t148 = -t204 * t317 + t238 * t299;
t147 = t277 * t162;
t144 = t241 * t320 - t243 * t301 + t245 * t302;
t143 = t240 * t320 - t242 * t301 + t244 * t302;
t142 = t241 * t319 - t243 * t299 + t245 * t300;
t141 = t240 * t319 - t242 * t299 + t244 * t300;
t138 = t204 * t301 - t205 * t299;
t135 = t205 * t336 + t320 * t429 + t249;
t134 = t238 * t319 + t336 * t432 + t265;
t133 = (-t304 + t432) * t357 + t365 * t385;
t132 = t205 * t357 + t362 * t385 + t427;
t129 = t369 - t435;
t128 = t370 + t436;
t127 = t171 * t293 - t197 * t277;
t126 = -t170 * t293 + t197 * t275;
t125 = -t193 * t277 + t150;
t124 = -t162 * t293 + t177;
t123 = t204 * t320 + t248 + (-t205 - t257) * t319;
t122 = t369 + t413 + t163;
t121 = -t422 - t276 * t351 + (-rSges(7,3) + t366) * t275 + t370 + t379;
t120 = (t204 * t362 + t205 * t365) * t356 + t398;
t113 = t170 * t277 - t171 * t275;
t110 = -t163 * t275 + t147;
t104 = t112 * t293;
t101 = t171 * t317 + t301 * t433 + t207;
t100 = t197 * t299 + t317 * t436 + t220;
t95 = (-t304 + t417) * t357 + t365 * t378;
t94 = t171 * t357 + t362 * t378 + t414;
t91 = t171 * t336 + t320 * t415 + t430;
t90 = t197 * t319 + t336 * t417 + t428;
t89 = t181 + (-t145 * t365 + t146 * t362) * t356;
t88 = t170 * t301 + t299 * t435 + t206;
t87 = t145 * t319 + t146 * t320 + t174;
t82 = (t170 * t362 + t171 * t365) * t356 + t376;
t81 = t173 * t357 + (-t143 * t365 + t144 * t362) * t356;
t80 = t172 * t357 + (-t141 * t365 + t142 * t362) * t356;
t77 = t170 * t320 + (-t257 + t435) * t319 + t431;
t72 = t143 * t319 + t144 * t320 + t173 * t336;
t71 = t141 * t319 + t142 * t320 + t172 * t336;
t66 = t152 * t293 - t277 * t434 + t150;
t65 = t189 * t275 - t293 * t438 + t177;
t64 = t301 * t416 + t317 * t437 + t207;
t63 = t299 * t434 + t317 * t419 + t220;
t62 = (-t304 + t400) * t357 + t365 * t375;
t61 = t357 * t437 + t362 * t375 + t414;
t60 = t320 * t399 + t336 * t437 + t430;
t59 = t319 * t434 + t336 * t400 + t428;
t58 = t151 * t277 - t275 * t437 + t147;
t56 = t299 * t418 + t301 * t438 + t206;
t54 = (t362 * t438 + t365 * t437) * t356 + t376;
t52 = t438 * t320 + (-t257 + t418) * t319 + t431;
t37 = t85 * t275 + t86 * t277 + t104;
t16 = t275 * t75 + t277 * t76 + t293 * t99;
t15 = t275 * t73 + t277 * t74 + t293 * t98;
t1 = [(t121 ^ 2 + t122 ^ 2) * m(7) + m(2) * (t345 ^ 2 + t346 ^ 2) + (t153 ^ 2 + t154 ^ 2) * m(5) + (t212 ^ 2 + t213 ^ 2) * m(4) + (t128 ^ 2 + t129 ^ 2) * m(6) + t412 + (t297 ^ 2 + t298 ^ 2) * m(3) + Icges(2,3) + t182 + t140 + t112 + t108; t283 + t139 + t111 + t106 + t181 + (t284 * t298 + t285 * t297) * m(3) + (t121 * t62 + t122 * t61) * m(7) + (t128 * t95 + t129 * t94) * m(6) + (t132 * t154 + t133 * t153) * m(5) + (t183 * t213 + t184 * t212) * m(4) + ((-t236 / 0.2e1 - t259 / 0.2e1 - t368) * t365 + (t260 / 0.2e1 + t237 / 0.2e1 + t367) * t362) * t356; (t54 ^ 2 + t61 ^ 2 + t62 ^ 2) * m(7) + (t82 ^ 2 + t94 ^ 2 + t95 ^ 2) * m(6) + (t120 ^ 2 + t132 ^ 2 + t133 ^ 2) * m(5) + (t155 ^ 2 + t183 ^ 2 + t184 ^ 2) * m(4) + (t261 ^ 2 + t284 ^ 2 + t285 ^ 2) * m(3) + (t43 + t45 + t57 + t89 + t283) * t357 + ((-t25 - t29 - t50 - t80 + (t308 * t337 + t310 * t338 - t356 * t448) * t443) * t365 + (t26 + t30 + t51 + t81 + ((t309 * t339 + t311 * t340 + (t307 * t362 - t448) * t356) * t362 + (t307 * t443 - t308 * t339 - t309 * t337 - t310 * t340 - t311 * t338) * t365) * t356) * t362 + ((-t236 - t259) * t365 + (t237 + t260) * t362) * t357) * t356; t105 + (t121 * t59 + t122 * t60) * m(7) + (t128 * t90 + t129 * t91) * m(6) + t109 + (t134 * t153 + t135 * t154) * m(5) + t137 + t174 + (t187 * t212 + t188 * t213) * m(4) + t367 * t320 + t368 * t319; (t77 * t82 + t90 * t95 + t91 * t94) * m(6) + (t52 * t54 + t59 * t62 + t60 * t61) * m(7) + (t120 * t123 + t132 * t135 + t133 * t134) * m(5) + (t155 * t178 + t183 * t188 + t184 * t187) * m(4) + (t87 / 0.2e1 + t390) * t357 + (t89 / 0.2e1 + t389) * t336 + (t81 / 0.2e1 + t392) * t320 + (t80 / 0.2e1 + t393) * t319 + ((-t71 / 0.2e1 - t395) * t365 + (t72 / 0.2e1 + t394) * t362) * t356; (t52 ^ 2 + t59 ^ 2 + t60 ^ 2) * m(7) + (t77 ^ 2 + t90 ^ 2 + t91 ^ 2) * m(6) + (t123 ^ 2 + t134 ^ 2 + t135 ^ 2) * m(5) + (t178 ^ 2 + t187 ^ 2 + t188 ^ 2) * m(4) + (t40 + t44 + t55 + t87) * t336 + (t22 + t28 + t49 + t72) * t320 + (t21 + t27 + t48 + t71) * t319; t103 + (t121 * t63 + t122 * t64) * m(7) + (t100 * t128 + t101 * t129) * m(6) + t107 + (t148 * t153 + t149 * t154) * m(5) + t136 + t371 * t301 + t372 * t299; (t54 * t56 + t61 * t64 + t62 * t63) * m(7) + (t100 * t95 + t101 * t94 + t82 * t88) * m(6) + (t120 * t138 + t132 * t149 + t133 * t148) * m(5) + t391 * t357 + t389 * t317 + t392 * t301 + t393 * t299 + (t362 * t396 - t365 * t397) * t356; (t52 * t56 + t59 * t63 + t60 * t64) * m(7) + (t100 * t90 + t101 * t91 + t77 * t88) * m(6) + (t123 * t138 + t134 * t148 + t135 * t149) * m(5) + t391 * t336 + t396 * t320 + t397 * t319 + t390 * t317 + t394 * t301 + t395 * t299; (t56 ^ 2 + t63 ^ 2 + t64 ^ 2) * m(7) + (t100 ^ 2 + t101 ^ 2 + t88 ^ 2) * m(6) + (t138 ^ 2 + t148 ^ 2 + t149 ^ 2) * m(5) + (t38 + t42 + t53) * t317 + (t18 + t24 + t47) * t301 + (t17 + t23 + t46) * t299; (t121 * t65 + t122 * t66) * m(7) + (t126 * t128 + t127 * t129) * m(6) + t104 + t421 * t277 + t420 * t275 + t386; (t54 * t58 + t61 * t66 + t62 * t65) * m(7) + (t113 * t82 + t126 * t95 + t127 * t94) * m(6) + t29 * t467 + t45 * t465 + t30 * t466 + t37 * t458 + (-t365 * t15 / 0.2e1 + t362 * t16 / 0.2e1) * t356 + t380; (t52 * t58 + t59 * t65 + t60 * t66) * m(7) + t16 * t460 + t15 * t461 + (t113 * t77 + t126 * t90 + t127 * t91) * m(6) + t37 * t459 + t28 * t466 + t27 * t467 + t44 * t465 + t381; (t56 * t58 + t63 * t65 + t64 * t66) * m(7) + t23 * t467 + t24 * t466 + t37 * t462 + t16 * t463 + (t100 * t126 + t101 * t127 + t113 * t88) * m(6) + t15 * t464 + t42 * t465 + t382; (t58 ^ 2 + t65 ^ 2 + t66 ^ 2) * m(7) + t293 * t37 + t277 * t16 + (t113 ^ 2 + t126 ^ 2 + t127 ^ 2) * m(6) + t275 * t15 + t423; (t121 * t124 + t122 * t125) * m(7) + t386; (t110 * t54 + t124 * t62 + t125 * t61) * m(7) + t380; (t110 * t52 + t124 * t59 + t125 * t60) * m(7) + t381; (t110 * t56 + t124 * t63 + t125 * t64) * m(7) + t382; (t110 * t58 + t124 * t65 + t125 * t66) * m(7) + t423; (t110 ^ 2 + t124 ^ 2 + t125 ^ 2) * m(7) + t423;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
