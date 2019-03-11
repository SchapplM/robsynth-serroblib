% Calculate joint inertia matrix for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:35
% EndTime: 2019-03-09 22:26:59
% DurationCPUTime: 10.36s
% Computational Cost: add. (37908->670), mult. (51392->899), div. (0->0), fcn. (63694->14), ass. (0->329)
t347 = qJ(3) + qJ(4);
t380 = pkin(12) + t347;
t340 = sin(t380);
t349 = cos(pkin(6));
t374 = cos(t380);
t348 = sin(pkin(6));
t352 = sin(qJ(2));
t430 = t348 * t352;
t301 = t340 * t430 - t349 * t374;
t365 = t348 * t374;
t302 = t349 * t340 + t352 * t365;
t356 = cos(qJ(2));
t427 = t348 * t356;
t245 = Icges(6,5) * t302 - Icges(6,6) * t301 - Icges(6,3) * t427;
t342 = sin(t347);
t343 = cos(t347);
t309 = -t342 * t430 + t343 * t349;
t310 = t342 * t349 + t343 * t430;
t251 = Icges(5,5) * t310 + Icges(5,6) * t309 - Icges(5,3) * t427;
t469 = -t245 - t251;
t353 = sin(qJ(1));
t422 = t353 * t356;
t357 = cos(qJ(1));
t423 = t352 * t357;
t324 = t349 * t423 + t422;
t279 = t324 * t340 + t357 * t365;
t426 = t348 * t357;
t280 = t324 * t374 - t340 * t426;
t421 = t356 * t357;
t424 = t352 * t353;
t323 = -t349 * t421 + t424;
t189 = Icges(6,5) * t280 - Icges(6,6) * t279 + Icges(6,3) * t323;
t287 = -t324 * t342 - t343 * t426;
t288 = t324 * t343 - t342 * t426;
t204 = Icges(5,5) * t288 + Icges(5,6) * t287 + Icges(5,3) * t323;
t468 = t189 + t204;
t326 = -t349 * t424 + t421;
t281 = t326 * t340 - t353 * t365;
t429 = t348 * t353;
t282 = t326 * t374 + t340 * t429;
t325 = t349 * t422 + t423;
t190 = Icges(6,5) * t282 - Icges(6,6) * t281 + Icges(6,3) * t325;
t289 = -t326 * t342 + t343 * t429;
t290 = t326 * t343 + t342 * t429;
t205 = Icges(5,5) * t290 + Icges(5,6) * t289 + Icges(5,3) * t325;
t467 = t190 + t205;
t246 = Icges(6,4) * t302 - Icges(6,2) * t301 - Icges(6,6) * t427;
t247 = Icges(6,1) * t302 - Icges(6,4) * t301 - Icges(6,5) * t427;
t252 = Icges(5,4) * t310 + Icges(5,2) * t309 - Icges(5,6) * t427;
t253 = Icges(5,1) * t310 + Icges(5,4) * t309 - Icges(5,5) * t427;
t466 = -t301 * t246 + t302 * t247 + t309 * t252 + t310 * t253;
t465 = t427 * t469 + t466;
t191 = Icges(6,4) * t280 - Icges(6,2) * t279 + Icges(6,6) * t323;
t193 = Icges(6,1) * t280 - Icges(6,4) * t279 + Icges(6,5) * t323;
t206 = Icges(5,4) * t288 + Icges(5,2) * t287 + Icges(5,6) * t323;
t208 = Icges(5,1) * t288 + Icges(5,4) * t287 + Icges(5,5) * t323;
t464 = -t191 * t279 + t193 * t280 + t206 * t287 + t208 * t288 + t323 * t468;
t192 = Icges(6,4) * t282 - Icges(6,2) * t281 + Icges(6,6) * t325;
t194 = Icges(6,1) * t282 - Icges(6,4) * t281 + Icges(6,5) * t325;
t207 = Icges(5,4) * t290 + Icges(5,2) * t289 + Icges(5,6) * t325;
t209 = Icges(5,1) * t290 + Icges(5,4) * t289 + Icges(5,5) * t325;
t463 = -t192 * t279 + t194 * t280 + t207 * t287 + t209 * t288 + t323 * t467;
t462 = -t191 * t281 + t193 * t282 + t206 * t289 + t208 * t290 + t325 * t468;
t461 = -t192 * t281 + t194 * t282 + t207 * t289 + t209 * t290 + t325 * t467;
t126 = t245 * t323 - t246 * t279 + t247 * t280;
t131 = t251 * t323 + t252 * t287 + t253 * t288;
t460 = t126 + t131;
t127 = t245 * t325 - t246 * t281 + t247 * t282;
t132 = t251 * t325 + t252 * t289 + t253 * t290;
t459 = t127 + t132;
t457 = t465 * t349;
t350 = sin(qJ(6));
t354 = cos(qJ(6));
t241 = -t280 * t350 + t323 * t354;
t242 = t280 * t354 + t323 * t350;
t368 = -rSges(7,1) * t242 - rSges(7,2) * t241;
t157 = rSges(7,3) * t279 - t368;
t439 = pkin(5) * t280;
t420 = pkin(11) * t279 + t157 + t439;
t277 = -t302 * t350 - t354 * t427;
t278 = t302 * t354 - t350 * t427;
t186 = rSges(7,1) * t278 + rSges(7,2) * t277 + rSges(7,3) * t301;
t456 = pkin(5) * t302 + pkin(11) * t301 + t186;
t151 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t279;
t153 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t279;
t155 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t279;
t69 = t151 * t279 + t153 * t241 + t155 * t242;
t243 = -t282 * t350 + t325 * t354;
t244 = t282 * t354 + t325 * t350;
t152 = Icges(7,5) * t244 + Icges(7,6) * t243 + Icges(7,3) * t281;
t154 = Icges(7,4) * t244 + Icges(7,2) * t243 + Icges(7,6) * t281;
t156 = Icges(7,1) * t244 + Icges(7,4) * t243 + Icges(7,5) * t281;
t70 = t152 * t279 + t154 * t241 + t156 * t242;
t183 = Icges(7,5) * t278 + Icges(7,6) * t277 + Icges(7,3) * t301;
t184 = Icges(7,4) * t278 + Icges(7,2) * t277 + Icges(7,6) * t301;
t185 = Icges(7,1) * t278 + Icges(7,4) * t277 + Icges(7,5) * t301;
t84 = t183 * t279 + t184 * t241 + t185 * t242;
t11 = t323 * t69 + t325 * t70 - t427 * t84;
t455 = t323 * t464 + t463 * t325 - t460 * t427 + t11;
t71 = t151 * t281 + t153 * t243 + t155 * t244;
t72 = t152 * t281 + t154 * t243 + t156 * t244;
t85 = t183 * t281 + t184 * t243 + t185 * t244;
t12 = t323 * t71 + t325 * t72 - t427 * t85;
t454 = t323 * t462 + t325 * t461 - t427 * t459 + t12;
t15 = t349 * t84 + (t353 * t70 - t357 * t69) * t348;
t453 = t15 + t460 * t349 + (t463 * t353 - t357 * t464) * t348;
t16 = t349 * t85 + (t353 * t72 - t357 * t71) * t348;
t452 = t16 + t459 * t349 + (t353 * t461 - t357 * t462) * t348;
t76 = t152 * t301 + t154 * t277 + t156 * t278;
t436 = t76 * t325;
t75 = t151 * t301 + t153 * t277 + t155 * t278;
t437 = t75 * t323;
t90 = t301 * t183 + t277 * t184 + t278 * t185;
t21 = -t427 * t90 + t436 + t437;
t121 = -t205 * t427 + t207 * t309 + t209 * t310;
t432 = t121 * t325;
t120 = -t204 * t427 + t206 * t309 + t208 * t310;
t433 = t120 * t323;
t107 = -t190 * t427 - t192 * t301 + t194 * t302;
t434 = t107 * t325;
t106 = -t189 * t427 - t191 * t301 + t193 * t302;
t435 = t106 * t323;
t451 = -t427 * t465 + t21 + t432 + t433 + t434 + t435;
t89 = t90 * t349;
t23 = t89 + (t76 * t353 - t75 * t357) * t348;
t450 = t23 + t457 + ((-t106 - t120) * t357 + (t107 + t121) * t353) * t348;
t449 = t348 ^ 2;
t358 = -pkin(10) - pkin(9);
t448 = t279 / 0.2e1;
t447 = t281 / 0.2e1;
t446 = t301 / 0.2e1;
t445 = t323 / 0.2e1;
t444 = t325 / 0.2e1;
t443 = t349 / 0.2e1;
t442 = t353 / 0.2e1;
t441 = -t357 / 0.2e1;
t351 = sin(qJ(3));
t440 = pkin(3) * t351;
t355 = cos(qJ(3));
t341 = t355 * pkin(3) + pkin(2);
t438 = -pkin(2) + t341;
t260 = Icges(3,5) * t324 - Icges(3,6) * t323 - Icges(3,3) * t426;
t431 = t260 * t357;
t428 = t348 * t355;
t425 = t349 * t351;
t158 = t244 * rSges(7,1) + t243 * rSges(7,2) + t281 * rSges(7,3);
t419 = t282 * pkin(5) + pkin(11) * t281 + t158;
t330 = pkin(4) * t342 + t440;
t313 = t330 * t426;
t396 = t351 * t426;
t332 = pkin(3) * t396;
t346 = -qJ(5) + t358;
t402 = t346 - t358;
t329 = pkin(4) * t343 + t341;
t404 = t329 - t341;
t176 = -t323 * t402 + t324 * t404 - t313 + t332;
t167 = t325 * t176;
t369 = -rSges(6,1) * t280 + rSges(6,2) * t279;
t195 = rSges(6,3) * t323 - t369;
t418 = t325 * t195 + t167;
t233 = (t330 - t440) * t349 + (t352 * t404 + t356 * t402) * t348;
t417 = t176 * t427 + t323 * t233;
t397 = t351 * t429;
t386 = -pkin(3) * t397 + t325 * t358 - t326 * t341;
t388 = -t325 * t346 + t326 * t329 + t330 * t429;
t177 = t386 + t388;
t196 = t282 * rSges(6,1) - t281 * rSges(6,2) + t325 * rSges(6,3);
t416 = -t177 - t196;
t319 = t323 * pkin(9);
t212 = -t323 * t358 + t324 * t438 - t319 - t332;
t266 = pkin(3) * t425 + ((pkin(9) + t358) * t356 + t438 * t352) * t348;
t415 = t212 * t427 + t323 * t266;
t284 = t326 * pkin(2) + pkin(9) * t325;
t213 = -t284 - t386;
t272 = t349 * t284;
t414 = t349 * t213 + t272;
t211 = t290 * rSges(5,1) + t289 * rSges(5,2) + t325 * rSges(5,3);
t413 = -t211 - t213;
t283 = t324 * pkin(2) + t319;
t412 = -t212 - t283;
t293 = -t324 * t351 - t355 * t426;
t294 = t324 * t355 - t396;
t224 = rSges(4,1) * t294 + rSges(4,2) * t293 + rSges(4,3) * t323;
t410 = -t224 - t283;
t248 = rSges(6,1) * t302 - rSges(6,2) * t301 - rSges(6,3) * t427;
t408 = -t233 - t248;
t370 = -t288 * rSges(5,1) - t287 * rSges(5,2);
t210 = t323 * rSges(5,3) - t370;
t254 = rSges(5,1) * t310 + rSges(5,2) * t309 - rSges(5,3) * t427;
t149 = t210 * t427 + t323 * t254;
t321 = t349 * t355 - t351 * t430;
t322 = t352 * t428 + t425;
t257 = Icges(4,4) * t322 + Icges(4,2) * t321 - Icges(4,6) * t427;
t258 = Icges(4,1) * t322 + Icges(4,4) * t321 - Icges(4,5) * t427;
t407 = t321 * t257 + t322 * t258;
t406 = -t254 - t266;
t405 = t283 * t429 + t284 * t426;
t403 = t357 * pkin(1) + pkin(8) * t429;
t400 = t76 / 0.2e1 + t85 / 0.2e1;
t399 = t84 / 0.2e1 + t75 / 0.2e1;
t398 = -t90 - t465;
t395 = t325 * t420 + t167;
t394 = -t177 - t419;
t393 = t349 * t177 + t414;
t392 = -t176 + t412;
t391 = -t213 + t416;
t390 = -t233 - t456;
t389 = -t266 + t408;
t295 = -t326 * t351 + t353 * t428;
t296 = t326 * t355 + t397;
t225 = t296 * rSges(4,1) + t295 * rSges(4,2) + t325 * rSges(4,3);
t304 = Icges(3,3) * t349 + (Icges(3,5) * t352 + Icges(3,6) * t356) * t348;
t305 = Icges(3,6) * t349 + (Icges(3,4) * t352 + Icges(3,2) * t356) * t348;
t306 = Icges(3,5) * t349 + (Icges(3,1) * t352 + Icges(3,4) * t356) * t348;
t387 = t349 * t304 + t305 * t427 + t306 * t430;
t268 = t326 * rSges(3,1) - t325 * rSges(3,2) + rSges(3,3) * t429;
t384 = -t427 / 0.2e1;
t219 = Icges(4,5) * t296 + Icges(4,6) * t295 + Icges(4,3) * t325;
t221 = Icges(4,4) * t296 + Icges(4,2) * t295 + Icges(4,6) * t325;
t223 = Icges(4,1) * t296 + Icges(4,4) * t295 + Icges(4,5) * t325;
t125 = -t219 * t427 + t221 * t321 + t223 * t322;
t256 = Icges(4,5) * t322 + Icges(4,6) * t321 - Icges(4,3) * t427;
t136 = t256 * t325 + t257 * t295 + t258 * t296;
t382 = t125 / 0.2e1 + t136 / 0.2e1;
t218 = Icges(4,5) * t294 + Icges(4,6) * t293 + Icges(4,3) * t323;
t220 = Icges(4,4) * t294 + Icges(4,2) * t293 + Icges(4,6) * t323;
t222 = Icges(4,1) * t294 + Icges(4,4) * t293 + Icges(4,5) * t323;
t124 = -t218 * t427 + t220 * t321 + t222 * t322;
t135 = t256 * t323 + t257 * t293 + t258 * t294;
t381 = t135 / 0.2e1 + t124 / 0.2e1;
t379 = -t353 * pkin(1) + pkin(8) * t426;
t259 = rSges(4,1) * t322 + rSges(4,2) * t321 - rSges(4,3) * t427;
t327 = (pkin(2) * t352 - pkin(9) * t356) * t348;
t378 = t348 * (-t259 - t327);
t377 = -t213 + t394;
t376 = -t266 + t390;
t375 = t212 * t429 + t213 * t426 + t405;
t92 = t195 * t427 + t323 * t248 + t417;
t373 = t348 * (-t327 + t406);
t87 = t90 * t301;
t18 = t75 * t279 + t76 * t281 + t87;
t3 = t279 * t69 + t281 * t70 + t301 * t84;
t4 = t279 * t71 + t281 * t72 + t301 * t85;
t372 = t11 * t448 + t12 * t447 + t18 * t384 + t21 * t446 + t3 * t445 + t4 * t444;
t371 = t323 * t455 + t325 * t454;
t367 = t388 + t403;
t366 = t348 * (-t327 + t389);
t364 = t176 * t429 + t177 * t426 + t375;
t67 = t456 * t323 + t420 * t427 + t417;
t363 = t348 * (-t327 + t376);
t362 = -t324 * t329 + t313 + t379;
t267 = rSges(3,1) * t324 - rSges(3,2) * t323 - rSges(3,3) * t426;
t361 = -t427 * t451 + t371;
t360 = t433 / 0.2e1 + t432 / 0.2e1 + t437 / 0.2e1 + t436 / 0.2e1 + t435 / 0.2e1 + t434 / 0.2e1 + (t84 + t460) * t445 + (t85 + t459) * t444;
t359 = t453 * t445 + t452 * t444 + t451 * t443 + t454 * t429 / 0.2e1 + t450 * t384 - t455 * t426 / 0.2e1;
t335 = rSges(2,1) * t357 - rSges(2,2) * t353;
t334 = -rSges(2,1) * t353 - rSges(2,2) * t357;
t307 = rSges(3,3) * t349 + (rSges(3,1) * t352 + rSges(3,2) * t356) * t348;
t265 = Icges(3,1) * t326 - Icges(3,4) * t325 + Icges(3,5) * t429;
t264 = Icges(3,1) * t324 - Icges(3,4) * t323 - Icges(3,5) * t426;
t263 = Icges(3,4) * t326 - Icges(3,2) * t325 + Icges(3,6) * t429;
t262 = Icges(3,4) * t324 - Icges(3,2) * t323 - Icges(3,6) * t426;
t261 = Icges(3,5) * t326 - Icges(3,6) * t325 + Icges(3,3) * t429;
t250 = t268 + t403;
t249 = -t267 + t379;
t230 = -t267 * t349 - t307 * t426;
t229 = t268 * t349 - t307 * t429;
t214 = t387 * t349;
t188 = (t267 * t353 + t268 * t357) * t348;
t182 = t304 * t429 - t305 * t325 + t306 * t326;
t181 = -t304 * t426 - t305 * t323 + t306 * t324;
t180 = t325 * t212;
t179 = t325 * t210;
t170 = t284 + t225 + t403;
t169 = t379 + t410;
t164 = -t225 * t427 - t259 * t325;
t163 = t224 * t427 + t259 * t323;
t162 = t261 * t349 + (t263 * t356 + t265 * t352) * t348;
t161 = t260 * t349 + (t262 * t356 + t264 * t352) * t348;
t160 = -t386 + t211 + t403;
t159 = -t324 * t341 + t332 + (-rSges(5,3) + t358) * t323 + t370 + t379;
t150 = -t211 * t427 - t254 * t325;
t147 = t367 + t196;
t146 = (-rSges(6,3) + t346) * t323 + t362 + t369;
t145 = -t256 * t427 + t407;
t143 = t145 * t349;
t142 = t224 * t325 - t225 * t323;
t141 = t349 * t410 + t357 * t378;
t140 = t225 * t349 + t353 * t378 + t272;
t138 = -t211 * t323 + t179;
t130 = (t224 * t353 + t225 * t357) * t348 + t405;
t119 = t367 + t419;
t118 = -t439 + t323 * t346 + (-rSges(7,3) - pkin(11)) * t279 + t362 + t368;
t117 = t325 * t406 + t413 * t427;
t116 = t149 + t415;
t115 = t219 * t325 + t221 * t295 + t223 * t296;
t114 = t218 * t325 + t220 * t295 + t222 * t296;
t113 = t219 * t323 + t221 * t293 + t223 * t294;
t112 = t218 * t323 + t220 * t293 + t222 * t294;
t111 = t158 * t301 - t186 * t281;
t110 = -t157 * t301 + t186 * t279;
t101 = (-t210 + t412) * t349 + t357 * t373;
t100 = t211 * t349 + t353 * t373 + t414;
t93 = t325 * t408 + t416 * t427;
t91 = t157 * t281 - t158 * t279;
t88 = t323 * t413 + t179 + t180;
t86 = (t210 * t353 + t211 * t357) * t348 + t375;
t81 = t323 * t416 + t418;
t80 = t325 * t389 + t391 * t427;
t79 = t92 + t415;
t78 = (-t195 + t392) * t349 + t357 * t366;
t77 = t196 * t349 + t353 * t366 + t393;
t68 = t325 * t390 + t394 * t427;
t66 = t323 * t391 + t180 + t418;
t65 = (t195 * t353 + t196 * t357) * t348 + t364;
t64 = t143 + (-t124 * t357 + t125 * t353) * t348;
t63 = t124 * t323 + t125 * t325 - t145 * t427;
t62 = t325 * t376 + t377 * t427;
t61 = t67 + t415;
t60 = t323 * t394 + t395;
t59 = (t392 - t420) * t349 + t357 * t363;
t58 = t349 * t419 + t353 * t363 + t393;
t55 = t136 * t349 + (-t114 * t357 + t115 * t353) * t348;
t54 = t135 * t349 + (-t112 * t357 + t113 * t353) * t348;
t51 = t114 * t323 + t115 * t325 - t136 * t427;
t50 = t112 * t323 + t113 * t325 - t135 * t427;
t35 = t323 * t377 + t180 + t395;
t34 = (t353 * t420 + t357 * t419) * t348 + t364;
t1 = [Icges(2,3) + (-t256 + t469) * t427 + m(7) * (t118 ^ 2 + t119 ^ 2) + m(6) * (t146 ^ 2 + t147 ^ 2) + m(5) * (t159 ^ 2 + t160 ^ 2) + m(4) * (t169 ^ 2 + t170 ^ 2) + m(3) * (t249 ^ 2 + t250 ^ 2) + m(2) * (t334 ^ 2 + t335 ^ 2) + t387 + t90 + t407 + t466; t89 + t143 + t214 + m(4) * (t140 * t170 + t141 * t169) + m(5) * (t100 * t160 + t101 * t159) + m(3) * (t229 * t250 + t230 * t249) + m(6) * (t146 * t78 + t147 * t77) + m(7) * (t118 * t59 + t119 * t58) + ((-t120 / 0.2e1 - t106 / 0.2e1 - t161 / 0.2e1 - t181 / 0.2e1 - t131 / 0.2e1 - t126 / 0.2e1 - t381 - t399) * t357 + (t121 / 0.2e1 + t107 / 0.2e1 + t162 / 0.2e1 + t182 / 0.2e1 + t127 / 0.2e1 + t132 / 0.2e1 + t382 + t400) * t353) * t348 + t457; (t64 + t214 + t450) * t349 + m(7) * (t34 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(6) * (t65 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2 + t86 ^ 2) + m(4) * (t130 ^ 2 + t140 ^ 2 + t141 ^ 2) + m(3) * (t188 ^ 2 + t229 ^ 2 + t230 ^ 2) + ((-t54 + ((-t262 * t323 + t264 * t324) * t348 - t449 * t431) * t357 - t453) * t357 + (t55 + ((-t263 * t325 + t265 * t326 + (t261 * t353 - t431) * t348) * t353 + (t261 * t426 + t262 * t325 + t263 * t323 - t264 * t326 - t265 * t324) * t357) * t348 + t452) * t353 + ((-t161 - t181) * t357 + (t162 + t182) * t353) * t349) * t348; t382 * t325 + t381 * t323 + t360 + m(7) * (t118 * t61 + t119 * t62) + m(6) * (t146 * t79 + t147 * t80) + m(5) * (t116 * t159 + t117 * t160) + m(4) * (t163 * t169 + t164 * t170) + (-t145 + t398) * t427; t359 + t63 * t443 + (t51 * t442 - t356 * t64 / 0.2e1 + t50 * t441) * t348 + m(4) * (t130 * t142 + t140 * t164 + t141 * t163) + m(5) * (t100 * t117 + t101 * t116 + t86 * t88) + m(7) * (t34 * t35 + t58 * t62 + t59 * t61) + m(6) * (t65 * t66 + t77 * t80 + t78 * t79) + t54 * t445 + t55 * t444; t323 * t50 + t325 * t51 + (-t63 - t451) * t427 + m(7) * (t35 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t66 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t116 ^ 2 + t117 ^ 2 + t88 ^ 2) + m(4) * (t142 ^ 2 + t163 ^ 2 + t164 ^ 2) + t371; m(7) * (t118 * t67 + t119 * t68) + m(6) * (t146 * t92 + t147 * t93) + m(5) * (t149 * t159 + t150 * t160) + t398 * t427 + t360; m(7) * (t34 * t60 + t58 * t68 + t59 * t67) + m(6) * (t65 * t81 + t77 * t93 + t78 * t92) + m(5) * (t100 * t150 + t101 * t149 + t138 * t86) + t359; m(7) * (t35 * t60 + t61 * t67 + t62 * t68) + m(6) * (t66 * t81 + t79 * t92 + t80 * t93) + m(5) * (t116 * t149 + t117 * t150 + t138 * t88) + t361; m(7) * (t60 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t81 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(5) * (t138 ^ 2 + t149 ^ 2 + t150 ^ 2) + t361; m(7) * (t118 * t325 + t119 * t323) + m(6) * (t146 * t325 + t147 * t323); m(7) * (t323 * t58 + t325 * t59 - t34 * t427) + m(6) * (t323 * t77 + t325 * t78 - t427 * t65); m(7) * (t323 * t62 + t325 * t61 - t35 * t427) + m(6) * (t323 * t80 + t325 * t79 - t427 * t66); m(7) * (t323 * t68 + t325 * t67 - t427 * t60) + m(6) * (t323 * t93 + t325 * t92 - t427 * t81); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t356 ^ 2 * t449 + t323 ^ 2 + t325 ^ 2); m(7) * (t110 * t118 + t111 * t119) + t87 + t400 * t281 + t399 * t279; m(7) * (t110 * t59 + t111 * t58 + t34 * t91) + t15 * t448 + t23 * t446 + t18 * t443 + t16 * t447 + (t3 * t441 + t4 * t442) * t348; m(7) * (t110 * t61 + t111 * t62 + t35 * t91) + t372; m(7) * (t110 * t67 + t111 * t68 + t60 * t91) + t372; m(7) * (t110 * t325 + t111 * t323 - t427 * t91); t281 * t4 + t279 * t3 + t301 * t18 + m(7) * (t110 ^ 2 + t111 ^ 2 + t91 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
