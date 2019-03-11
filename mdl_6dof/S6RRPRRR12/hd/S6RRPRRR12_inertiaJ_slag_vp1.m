% Calculate joint inertia matrix for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR12_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR12_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:18
% EndTime: 2019-03-09 14:35:30
% DurationCPUTime: 6.05s
% Computational Cost: add. (23188->625), mult. (40469->854), div. (0->0), fcn. (50594->12), ass. (0->296)
t366 = qJ(4) + qJ(5);
t303 = sin(t366);
t306 = cos(pkin(6));
t313 = cos(qJ(2));
t305 = sin(pkin(6));
t338 = cos(t366);
t330 = t305 * t338;
t269 = t306 * t303 + t313 * t330;
t373 = t305 * t313;
t270 = -t303 * t373 + t306 * t338;
t309 = sin(qJ(2));
t375 = t305 * t309;
t202 = Icges(6,5) * t270 - Icges(6,6) * t269 + Icges(6,3) * t375;
t203 = Icges(6,4) * t270 - Icges(6,2) * t269 + Icges(6,6) * t375;
t204 = Icges(6,1) * t270 - Icges(6,4) * t269 + Icges(6,5) * t375;
t310 = sin(qJ(1));
t368 = t310 * t313;
t314 = cos(qJ(1));
t369 = t309 * t314;
t286 = t306 * t368 + t369;
t374 = t305 * t310;
t241 = -t286 * t338 + t303 * t374;
t242 = t286 * t303 + t310 * t330;
t367 = t313 * t314;
t370 = t309 * t310;
t287 = -t306 * t370 + t367;
t101 = t202 * t287 - t203 * t241 + t204 * t242;
t285 = t306 * t369 + t368;
t307 = sin(qJ(6));
t311 = cos(qJ(6));
t196 = -t242 * t307 + t287 * t311;
t197 = t242 * t311 + t287 * t307;
t120 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t241;
t122 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t241;
t124 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t241;
t50 = t120 * t241 + t122 * t196 + t124 * t197;
t284 = -t306 * t367 + t370;
t244 = t284 * t303 - t314 * t330;
t198 = -t244 * t307 + t285 * t311;
t199 = t244 * t311 + t285 * t307;
t372 = t305 * t314;
t243 = t284 * t338 + t303 * t372;
t121 = Icges(7,5) * t199 + Icges(7,6) * t198 - Icges(7,3) * t243;
t123 = Icges(7,4) * t199 + Icges(7,2) * t198 - Icges(7,6) * t243;
t125 = Icges(7,1) * t199 + Icges(7,4) * t198 - Icges(7,5) * t243;
t51 = t121 * t241 + t123 * t196 + t125 * t197;
t239 = -t270 * t307 + t311 * t375;
t240 = t270 * t311 + t307 * t375;
t150 = Icges(7,5) * t240 + Icges(7,6) * t239 + Icges(7,3) * t269;
t151 = Icges(7,4) * t240 + Icges(7,2) * t239 + Icges(7,6) * t269;
t152 = Icges(7,1) * t240 + Icges(7,4) * t239 + Icges(7,5) * t269;
t64 = t150 * t241 + t151 * t196 + t152 * t197;
t11 = t285 * t51 + t287 * t50 + t375 * t64;
t156 = Icges(6,5) * t242 - Icges(6,6) * t241 + Icges(6,3) * t287;
t158 = Icges(6,4) * t242 - Icges(6,2) * t241 + Icges(6,6) * t287;
t160 = Icges(6,1) * t242 - Icges(6,4) * t241 + Icges(6,5) * t287;
t76 = t156 * t287 - t158 * t241 + t160 * t242;
t157 = Icges(6,5) * t244 + Icges(6,6) * t243 + Icges(6,3) * t285;
t159 = Icges(6,4) * t244 + Icges(6,2) * t243 + Icges(6,6) * t285;
t161 = Icges(6,1) * t244 + Icges(6,4) * t243 + Icges(6,5) * t285;
t77 = t157 * t287 - t159 * t241 + t161 * t242;
t404 = t101 * t375 + t285 * t77 + t287 * t76 + t11;
t102 = t202 * t285 + t203 * t243 + t204 * t244;
t52 = -t120 * t243 + t122 * t198 + t124 * t199;
t53 = -t121 * t243 + t123 * t198 + t125 * t199;
t65 = -t150 * t243 + t151 * t198 + t152 * t199;
t12 = t285 * t53 + t287 * t52 + t375 * t65;
t78 = t156 * t285 + t158 * t243 + t160 * t244;
t79 = t157 * t285 + t159 * t243 + t161 * t244;
t403 = t102 * t375 + t285 * t79 + t287 * t78 + t12;
t15 = t306 * t64 + (t310 * t50 - t314 * t51) * t305;
t402 = t15 + t101 * t306 + (t310 * t76 - t314 * t77) * t305;
t16 = t306 * t65 + (t310 * t52 - t314 * t53) * t305;
t401 = t16 + t102 * t306 + (t310 * t78 - t314 * t79) * t305;
t108 = t202 * t375 - t269 * t203 + t270 * t204;
t103 = t108 * t375;
t59 = t121 * t269 + t123 * t239 + t125 * t240;
t380 = t59 * t285;
t58 = t120 * t269 + t122 * t239 + t124 * t240;
t381 = t58 * t287;
t75 = t269 * t150 + t239 * t151 + t240 * t152;
t72 = t75 * t375;
t22 = t380 + t72 + t381;
t90 = t157 * t375 - t159 * t269 + t161 * t270;
t378 = t90 * t285;
t89 = t156 * t375 - t158 * t269 + t160 * t270;
t379 = t89 * t287;
t400 = t22 + t103 + t378 + t379;
t106 = t108 * t306;
t73 = t75 * t306;
t24 = t73 + (t58 * t310 - t59 * t314) * t305;
t399 = t24 + t106 + (t89 * t310 - t90 * t314) * t305;
t128 = t197 * rSges(7,1) + t196 * rSges(7,2) + t241 * rSges(7,3);
t363 = t242 * pkin(5) + pkin(11) * t241 + t128;
t326 = -t199 * rSges(7,1) - t198 * rSges(7,2);
t129 = -t243 * rSges(7,3) - t326;
t383 = t244 * pkin(5);
t362 = -t243 * pkin(11) + t129 + t383;
t153 = rSges(7,1) * t240 + rSges(7,2) * t239 + rSges(7,3) * t269;
t398 = pkin(5) * t270 + pkin(11) * t269 + t153;
t261 = Icges(3,3) * t306 + (Icges(3,5) * t309 + Icges(3,6) * t313) * t305;
t262 = Icges(3,6) * t306 + (Icges(3,4) * t309 + Icges(3,2) * t313) * t305;
t263 = Icges(3,5) * t306 + (Icges(3,1) * t309 + Icges(3,4) * t313) * t305;
t264 = Icges(4,5) * t306 + (-Icges(4,6) * t309 - Icges(4,3) * t313) * t305;
t265 = Icges(4,4) * t306 + (-Icges(4,2) * t309 - Icges(4,6) * t313) * t305;
t266 = Icges(4,1) * t306 + (-Icges(4,4) * t309 - Icges(4,5) * t313) * t305;
t397 = (-t264 * t313 - t265 * t309) * t305 + t262 * t373 + t263 * t375 + (t266 + t261) * t306;
t213 = -Icges(4,5) * t372 - Icges(4,6) * t285 + Icges(4,3) * t284;
t220 = Icges(3,4) * t285 - Icges(3,2) * t284 - Icges(3,6) * t372;
t396 = t213 - t220;
t215 = -Icges(4,4) * t372 - Icges(4,2) * t285 + Icges(4,6) * t284;
t222 = Icges(3,1) * t285 - Icges(3,4) * t284 - Icges(3,5) * t372;
t395 = t215 - t222;
t212 = Icges(4,5) * t374 - Icges(4,6) * t287 + Icges(4,3) * t286;
t221 = Icges(3,4) * t287 - Icges(3,2) * t286 + Icges(3,6) * t374;
t394 = -t221 + t212;
t214 = Icges(4,4) * t374 - Icges(4,2) * t287 + Icges(4,6) * t286;
t223 = Icges(3,1) * t287 - Icges(3,4) * t286 + Icges(3,5) * t374;
t393 = t223 - t214;
t392 = t305 ^ 2;
t391 = t241 / 0.2e1;
t390 = -t243 / 0.2e1;
t389 = t269 / 0.2e1;
t388 = t285 / 0.2e1;
t387 = t287 / 0.2e1;
t386 = t306 / 0.2e1;
t385 = t310 / 0.2e1;
t384 = -t314 / 0.2e1;
t315 = -pkin(10) - pkin(9);
t382 = -pkin(2) + t315;
t308 = sin(qJ(4));
t377 = t284 * t308;
t376 = t286 * t308;
t371 = t308 * t313;
t365 = t362 * t287;
t364 = t363 * t375;
t361 = t398 * t285;
t359 = t397 * t306;
t217 = -Icges(4,1) * t372 - Icges(4,4) * t285 + Icges(4,5) * t284;
t218 = Icges(3,5) * t285 - Icges(3,6) * t284 - Icges(3,3) * t372;
t358 = -t218 - t217;
t216 = Icges(4,1) * t374 - Icges(4,4) * t287 + Icges(4,5) * t286;
t219 = Icges(3,5) * t287 - Icges(3,6) * t286 + Icges(3,3) * t374;
t357 = t219 + t216;
t273 = t284 * qJ(3);
t233 = t285 * pkin(2) + t273;
t234 = t287 * pkin(2) + qJ(3) * t286;
t356 = t233 * t374 + t234 * t372;
t232 = t306 * t234;
t256 = pkin(3) * t374 + pkin(9) * t287;
t355 = t306 * t256 + t232;
t351 = pkin(3) * t372 - t285 * pkin(9);
t354 = -t233 + t351;
t288 = (pkin(2) * t309 - qJ(3) * t313) * t305;
t353 = -t306 * pkin(3) - pkin(9) * t375 - t288;
t352 = t314 * pkin(1) + pkin(8) * t374;
t350 = t58 / 0.2e1 + t64 / 0.2e1;
t349 = -t59 / 0.2e1 - t65 / 0.2e1;
t312 = cos(qJ(4));
t282 = -t306 * t308 - t312 * t373;
t283 = -t305 * t371 + t306 * t312;
t207 = Icges(5,5) * t283 + Icges(5,6) * t282 + Icges(5,3) * t375;
t208 = Icges(5,4) * t283 + Icges(5,2) * t282 + Icges(5,6) * t375;
t209 = Icges(5,1) * t283 + Icges(5,4) * t282 + Icges(5,5) * t375;
t250 = t286 * t312 - t308 * t374;
t251 = t312 * t374 + t376;
t104 = t207 * t287 + t208 * t250 + t209 * t251;
t171 = Icges(5,5) * t251 + Icges(5,6) * t250 + Icges(5,3) * t287;
t173 = Icges(5,4) * t251 + Icges(5,2) * t250 + Icges(5,6) * t287;
t175 = Icges(5,1) * t251 + Icges(5,4) * t250 + Icges(5,5) * t287;
t95 = t171 * t375 + t173 * t282 + t175 * t283;
t348 = t104 / 0.2e1 + t95 / 0.2e1;
t252 = t284 * t312 + t308 * t372;
t253 = -t312 * t372 + t377;
t105 = t207 * t285 + t208 * t252 + t209 * t253;
t172 = Icges(5,5) * t253 + Icges(5,6) * t252 + Icges(5,3) * t285;
t174 = Icges(5,4) * t253 + Icges(5,2) * t252 + Icges(5,6) * t285;
t176 = Icges(5,1) * t253 + Icges(5,4) * t252 + Icges(5,5) * t285;
t96 = t172 * t375 + t174 * t282 + t176 * t283;
t347 = t105 / 0.2e1 + t96 / 0.2e1;
t302 = pkin(4) * t312 + pkin(3);
t342 = pkin(4) * t376 - t287 * t315 + t302 * t374;
t179 = -t256 + t342;
t346 = t306 * t179 + t355;
t334 = -pkin(4) * t377 + t302 * t372;
t180 = -t285 * t315 - t334 + t351;
t345 = -t180 + t354;
t115 = t207 * t375 + t282 * t208 + t283 * t209;
t230 = (-pkin(3) + t302) * t306 + (-pkin(4) * t371 + (-pkin(9) - t315) * t309) * t305;
t344 = -t230 + t353;
t162 = t242 * rSges(6,1) - t241 * rSges(6,2) + t287 * rSges(6,3);
t177 = t251 * rSges(5,1) + t250 * rSges(5,2) + t287 * rSges(5,3);
t227 = t287 * rSges(3,1) - t286 * rSges(3,2) + rSges(3,3) * t374;
t224 = rSges(4,1) * t374 - t287 * rSges(4,2) + t286 * rSges(4,3);
t341 = t375 / 0.2e1;
t337 = -t310 * pkin(1) + pkin(8) * t372;
t336 = t305 * (-rSges(4,1) * t306 - (-rSges(4,2) * t309 - rSges(4,3) * t313) * t305 - t288);
t335 = t256 * t372 - t351 * t374 + t356;
t333 = -t273 + t337;
t210 = rSges(5,1) * t283 + rSges(5,2) * t282 + rSges(5,3) * t375;
t332 = t305 * (-t210 + t353);
t69 = t75 * t269;
t18 = t58 * t241 - t59 * t243 + t69;
t3 = t241 * t50 - t243 * t51 + t269 * t64;
t4 = t241 * t52 - t243 * t53 + t269 * t65;
t331 = t11 * t391 + t12 * t390 + t18 * t341 + t22 * t389 + t3 * t387 + t4 * t388;
t329 = t403 * t285 + t287 * t404 + t400 * t375;
t328 = -rSges(5,1) * t253 - rSges(5,2) * t252;
t327 = -t244 * rSges(6,1) - t243 * rSges(6,2);
t205 = rSges(6,1) * t270 - rSges(6,2) * t269 + rSges(6,3) * t375;
t325 = t305 * (-t205 + t344);
t324 = t234 + t352;
t323 = rSges(4,1) * t372 - rSges(4,3) * t284;
t322 = t179 * t372 + t180 * t374 + t335;
t321 = t305 * (t344 - t398);
t226 = rSges(3,1) * t285 - rSges(3,2) * t284 - rSges(3,3) * t372;
t319 = t333 + t334;
t318 = t324 + t342;
t317 = t103 + t381 / 0.2e1 + t380 / 0.2e1 + t72 + t379 / 0.2e1 + t378 / 0.2e1 + (t65 + t102) * t388 + (t64 + t101) * t387;
t316 = t401 * t388 + t402 * t387 + t400 * t386 + t399 * t341 + t404 * t374 / 0.2e1 - t403 * t372 / 0.2e1;
t293 = rSges(2,1) * t314 - rSges(2,2) * t310;
t292 = -rSges(2,1) * t310 - rSges(2,2) * t314;
t267 = rSges(3,3) * t306 + (rSges(3,1) * t309 + rSges(3,2) * t313) * t305;
t225 = -rSges(4,2) * t285 - t323;
t201 = t227 + t352;
t200 = -t226 + t337;
t195 = t285 * t230;
t188 = t285 * t205;
t183 = -t226 * t306 - t267 * t372;
t182 = t227 * t306 - t267 * t374;
t178 = rSges(5,3) * t285 - t328;
t168 = t179 * t375;
t163 = t285 * rSges(6,3) - t327;
t155 = t324 + t224;
t154 = (rSges(4,2) - pkin(2)) * t285 + t323 + t333;
t149 = t287 * t180;
t148 = t162 * t375;
t147 = (t226 * t310 + t227 * t314) * t305;
t146 = t261 * t374 - t262 * t286 + t263 * t287;
t145 = -t261 * t372 - t262 * t284 + t263 * t285;
t144 = t264 * t284 - t265 * t285 - t266 * t372;
t143 = t264 * t286 - t265 * t287 + t266 * t374;
t142 = t287 * t163;
t137 = (-t225 - t233) * t306 + t314 * t336;
t136 = t224 * t306 + t310 * t336 + t232;
t135 = t177 * t375 - t210 * t287;
t134 = -t178 * t375 + t210 * t285;
t133 = t217 * t306 + (-t213 * t313 - t215 * t309) * t305;
t132 = t216 * t306 + (-t212 * t313 - t214 * t309) * t305;
t131 = t219 * t306 + (t221 * t313 + t223 * t309) * t305;
t130 = t218 * t306 + (t220 * t313 + t222 * t309) * t305;
t127 = t256 + t324 + t177;
t126 = (-rSges(5,3) - pkin(2)) * t285 + t328 + t333 + t351;
t118 = -t205 * t287 + t148;
t117 = -t163 * t375 + t188;
t114 = t115 * t306;
t113 = t115 * t375;
t112 = (t224 * t314 + t225 * t310) * t305 + t356;
t111 = t318 + t162;
t110 = (-rSges(6,3) + t382) * t285 + t319 + t327;
t109 = -t177 * t285 + t178 * t287;
t107 = -t162 * t285 + t142;
t100 = (-t178 + t354) * t306 + t314 * t332;
t99 = t177 * t306 + t310 * t332 + t355;
t94 = -t129 * t269 - t153 * t243;
t93 = t128 * t269 - t153 * t241;
t92 = t148 + t168 + (-t205 - t230) * t287;
t91 = t188 + t195 + (-t163 - t180) * t375;
t88 = t172 * t285 + t174 * t252 + t176 * t253;
t87 = t171 * t285 + t173 * t252 + t175 * t253;
t86 = t172 * t287 + t174 * t250 + t176 * t251;
t85 = t171 * t287 + t173 * t250 + t175 * t251;
t84 = t318 + t363;
t83 = -t383 + t382 * t285 + (rSges(7,3) + pkin(11)) * t243 + t319 + t326;
t82 = (t177 * t314 + t178 * t310) * t305 + t335;
t74 = t128 * t243 + t129 * t241;
t71 = (-t163 + t345) * t306 + t314 * t325;
t70 = t162 * t306 + t310 * t325 + t346;
t68 = t142 + t149 + (-t162 - t179) * t285;
t67 = -t287 * t398 + t364;
t66 = -t362 * t375 + t361;
t61 = (t162 * t314 + t163 * t310) * t305 + t322;
t60 = -t285 * t363 + t365;
t57 = t168 + (-t230 - t398) * t287 + t364;
t56 = t195 + (-t180 - t362) * t375 + t361;
t49 = (t345 - t362) * t306 + t314 * t321;
t48 = t306 * t363 + t310 * t321 + t346;
t47 = t149 + (-t179 - t363) * t285 + t365;
t46 = t114 + (t95 * t310 - t96 * t314) * t305;
t45 = t96 * t285 + t95 * t287 + t113;
t44 = (t310 * t362 + t314 * t363) * t305 + t322;
t41 = t105 * t306 + (t310 * t87 - t314 * t88) * t305;
t40 = t104 * t306 + (t310 * t85 - t314 * t86) * t305;
t36 = t105 * t375 + t285 * t88 + t287 * t87;
t35 = t104 * t375 + t285 * t86 + t287 * t85;
t1 = [Icges(2,3) + m(7) * (t83 ^ 2 + t84 ^ 2) + m(6) * (t110 ^ 2 + t111 ^ 2) + m(5) * (t126 ^ 2 + t127 ^ 2) + m(4) * (t154 ^ 2 + t155 ^ 2) + m(3) * (t200 ^ 2 + t201 ^ 2) + m(2) * (t292 ^ 2 + t293 ^ 2) + t115 + t108 + t75 + t397; t106 + t73 + t114 + m(3) * (t182 * t201 + t183 * t200) + m(4) * (t136 * t155 + t137 * t154) + m(7) * (t48 * t84 + t49 * t83) + m(6) * (t110 * t71 + t111 * t70) + m(5) * (t100 * t126 + t127 * t99) + ((-t90 / 0.2e1 - t133 / 0.2e1 - t130 / 0.2e1 - t145 / 0.2e1 - t102 / 0.2e1 - t144 / 0.2e1 - t347 + t349) * t314 + (t89 / 0.2e1 + t132 / 0.2e1 + t131 / 0.2e1 + t146 / 0.2e1 + t101 / 0.2e1 + t143 / 0.2e1 + t348 + t350) * t310) * t305 + t359; (t46 + t359 + t399) * t306 + m(7) * (t44 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t61 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t100 ^ 2 + t82 ^ 2 + t99 ^ 2) + m(4) * (t112 ^ 2 + t136 ^ 2 + t137 ^ 2) + m(3) * (t147 ^ 2 + t182 ^ 2 + t183 ^ 2) + ((-t41 + ((t284 * t396 - t395 * t285) * t305 + t358 * t392 * t314) * t314 + (-t130 - t133 - t144 - t145) * t306 - t401) * t314 + (t40 + ((t286 * t394 + t287 * t393) * t305 + t357 * t392 * t310) * t310 + (t146 + t143 + t132 + t131) * t306 + ((t310 * t358 + t314 * t357) * t305 + t395 * t287 - t396 * t286 - t393 * t285 - t394 * t284) * t372 + t402) * t310) * t305; m(7) * (t284 * t84 + t286 * t83) + m(6) * (t110 * t286 + t111 * t284) + m(5) * (t126 * t286 + t127 * t284) + m(4) * (t154 * t286 + t155 * t284); m(7) * (t284 * t48 + t286 * t49 - t373 * t44) + m(6) * (t284 * t70 + t286 * t71 - t373 * t61) + m(5) * (t100 * t286 + t284 * t99 - t373 * t82) + m(4) * (-t112 * t373 + t136 * t284 + t137 * t286); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t313 ^ 2 * t392 + t284 ^ 2 + t286 ^ 2); t317 + t113 + t347 * t285 + t348 * t287 + m(7) * (t56 * t83 + t57 * t84) + m(6) * (t110 * t91 + t111 * t92) + m(5) * (t126 * t134 + t127 * t135); t316 + m(7) * (t44 * t47 + t48 * t57 + t49 * t56) + m(6) * (t61 * t68 + t70 * t92 + t71 * t91) + m(5) * (t100 * t134 + t109 * t82 + t135 * t99) + (t309 * t46 / 0.2e1 + t35 * t385 + t36 * t384) * t305 + t45 * t386 + t40 * t387 + t41 * t388; m(5) * (-t109 * t373 + t134 * t286 + t135 * t284) + m(6) * (t284 * t92 + t286 * t91 - t373 * t68) + m(7) * (t284 * t57 + t286 * t56 - t373 * t47); t45 * t375 + t285 * t36 + t287 * t35 + m(7) * (t47 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t68 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t109 ^ 2 + t134 ^ 2 + t135 ^ 2) + t329; m(7) * (t66 * t83 + t67 * t84) + m(6) * (t110 * t117 + t111 * t118) + t317; t316 + m(7) * (t44 * t60 + t48 * t67 + t49 * t66) + m(6) * (t107 * t61 + t117 * t71 + t118 * t70); m(6) * (-t107 * t373 + t117 * t286 + t118 * t284) + m(7) * (t284 * t67 + t286 * t66 - t373 * t60); m(7) * (t47 * t60 + t56 * t66 + t57 * t67) + m(6) * (t107 * t68 + t117 * t91 + t118 * t92) + t329; m(7) * (t60 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t107 ^ 2 + t117 ^ 2 + t118 ^ 2) + t329; m(7) * (t83 * t94 + t84 * t93) + t69 + t349 * t243 + t350 * t241; t18 * t386 + t24 * t389 + t15 * t391 + t16 * t390 + m(7) * (t44 * t74 + t48 * t93 + t49 * t94) + (t3 * t385 + t384 * t4) * t305; m(7) * (t284 * t93 + t286 * t94 - t373 * t74); m(7) * (t47 * t74 + t56 * t94 + t57 * t93) + t331; m(7) * (t60 * t74 + t66 * t94 + t67 * t93) + t331; t241 * t3 - t243 * t4 + t269 * t18 + m(7) * (t74 ^ 2 + t93 ^ 2 + t94 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
