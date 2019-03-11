% Calculate joint inertia matrix for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR14_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR14_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:36
% EndTime: 2019-03-09 20:12:48
% DurationCPUTime: 5.37s
% Computational Cost: add. (23913->686), mult. (52791->931), div. (0->0), fcn. (67122->12), ass. (0->313)
t307 = cos(pkin(6));
t309 = sin(qJ(2));
t306 = sin(pkin(6));
t379 = sin(qJ(3));
t338 = t306 * t379;
t380 = cos(qJ(3));
t285 = -t307 * t380 + t309 * t338;
t305 = qJ(5) + qJ(6);
t302 = sin(t305);
t303 = cos(t305);
t312 = cos(qJ(2));
t370 = t306 * t312;
t250 = t285 * t303 + t302 * t370;
t251 = t285 * t302 - t303 * t370;
t339 = t306 * t380;
t286 = t307 * t379 + t309 * t339;
t165 = Icges(7,5) * t251 + Icges(7,6) * t250 + Icges(7,3) * t286;
t166 = Icges(7,4) * t251 + Icges(7,2) * t250 + Icges(7,6) * t286;
t167 = Icges(7,1) * t251 + Icges(7,4) * t250 + Icges(7,5) * t286;
t85 = t286 * t165 + t250 * t166 + t251 * t167;
t308 = sin(qJ(5));
t311 = cos(qJ(5));
t263 = t285 * t311 + t308 * t370;
t372 = t285 * t308;
t264 = -t311 * t370 + t372;
t170 = Icges(6,5) * t264 + Icges(6,6) * t263 + Icges(6,3) * t286;
t171 = Icges(6,4) * t264 + Icges(6,2) * t263 + Icges(6,6) * t286;
t172 = Icges(6,1) * t264 + Icges(6,4) * t263 + Icges(6,5) * t286;
t92 = t286 * t170 + t263 * t171 + t264 * t172;
t387 = -t85 - t92;
t313 = cos(qJ(1));
t369 = t306 * t313;
t310 = sin(qJ(1));
t366 = t310 * t312;
t367 = t309 * t313;
t288 = t307 * t367 + t366;
t266 = t288 * t380 - t313 * t338;
t386 = t266 / 0.2e1;
t365 = t312 * t313;
t368 = t309 * t310;
t290 = -t307 * t368 + t365;
t268 = t290 * t380 + t310 * t338;
t385 = t268 / 0.2e1;
t384 = t286 / 0.2e1;
t287 = -t307 * t365 + t368;
t383 = t287 / 0.2e1;
t289 = t307 * t366 + t367;
t382 = t289 / 0.2e1;
t381 = t307 / 0.2e1;
t301 = pkin(5) * t311 + pkin(4);
t378 = pkin(4) - t301;
t265 = t288 * t379 + t313 * t339;
t208 = t265 * t303 - t287 * t302;
t209 = t265 * t302 + t287 * t303;
t133 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t266;
t135 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t266;
t137 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t266;
t67 = t133 * t286 + t135 * t250 + t137 * t251;
t377 = t67 * t266;
t267 = t290 * t379 - t310 * t339;
t210 = t267 * t303 - t289 * t302;
t211 = t267 * t302 + t289 * t303;
t134 = Icges(7,5) * t211 + Icges(7,6) * t210 + Icges(7,3) * t268;
t136 = Icges(7,4) * t211 + Icges(7,2) * t210 + Icges(7,6) * t268;
t138 = Icges(7,1) * t211 + Icges(7,4) * t210 + Icges(7,5) * t268;
t68 = t134 * t286 + t136 * t250 + t138 * t251;
t376 = t68 * t268;
t235 = Icges(3,5) * t288 - Icges(3,6) * t287 - Icges(3,3) * t369;
t375 = t235 * t313;
t374 = t265 * t308;
t373 = t267 * t308;
t371 = t306 * t310;
t324 = -t209 * rSges(7,1) - t208 * rSges(7,2);
t141 = t266 * rSges(7,3) - t324;
t261 = t266 * pkin(10);
t314 = -pkin(11) - pkin(10);
t348 = pkin(5) * t374;
t151 = -t266 * t314 - t287 * t378 - t261 + t348;
t364 = t141 + t151;
t142 = t211 * rSges(7,1) + t210 * rSges(7,2) + t268 * rSges(7,3);
t226 = t289 * pkin(4) + pkin(10) * t268;
t341 = pkin(5) * t373 - t268 * t314 + t289 * t301;
t152 = -t226 + t341;
t363 = t142 + t152;
t168 = rSges(7,1) * t251 + rSges(7,2) * t250 + rSges(7,3) * t286;
t191 = pkin(5) * t372 + t378 * t370 + (-pkin(10) - t314) * t286;
t362 = t168 + t191;
t187 = t289 * rSges(5,1) - t268 * rSges(5,2) + t267 * rSges(5,3);
t201 = t268 * pkin(3) + qJ(4) * t267;
t361 = -t187 - t201;
t254 = t265 * qJ(4);
t200 = t266 * pkin(3) + t254;
t190 = t289 * t200;
t225 = t287 * pkin(4) + t261;
t360 = t289 * t225 + t190;
t247 = pkin(3) * t286 + qJ(4) * t285;
t359 = t200 * t370 + t287 * t247;
t249 = t290 * pkin(2) + pkin(9) * t289;
t246 = t307 * t249;
t358 = t307 * t201 + t246;
t248 = t288 * pkin(2) + t287 * pkin(9);
t357 = -t200 - t248;
t356 = -t201 - t226;
t227 = -Icges(5,5) * t370 - Icges(5,6) * t286 + Icges(5,3) * t285;
t228 = -Icges(5,4) * t370 - Icges(5,2) * t286 + Icges(5,6) * t285;
t355 = t285 * t227 - t286 * t228;
t231 = Icges(4,4) * t286 - Icges(4,2) * t285 - Icges(4,6) * t370;
t232 = Icges(4,1) * t286 - Icges(4,4) * t285 - Icges(4,5) * t370;
t354 = -t285 * t231 + t286 * t232;
t233 = -rSges(5,1) * t370 - rSges(5,2) * t286 + rSges(5,3) * t285;
t353 = -t233 - t247;
t352 = t248 * t371 + t249 * t369;
t271 = -pkin(4) * t370 + t286 * pkin(10);
t351 = -t247 - t271;
t350 = t313 * pkin(1) + pkin(8) * t371;
t82 = t85 * t286;
t26 = t376 + t82 + t377;
t55 = t133 * t266 + t135 * t208 + t137 * t209;
t56 = t134 * t266 + t136 * t208 + t138 * t209;
t77 = t165 * t266 + t166 * t208 + t167 * t209;
t7 = t266 * t55 + t268 * t56 + t286 * t77;
t57 = t133 * t268 + t135 * t210 + t137 * t211;
t58 = t134 * t268 + t136 * t210 + t138 * t211;
t78 = t165 * t268 + t166 * t210 + t167 * t211;
t8 = t266 * t57 + t268 * t58 + t286 * t78;
t349 = t286 * t26 + t266 * t7 + t268 * t8;
t221 = t265 * t311 - t287 * t308;
t222 = t287 * t311 + t374;
t143 = Icges(6,5) * t222 + Icges(6,6) * t221 + Icges(6,3) * t266;
t145 = Icges(6,4) * t222 + Icges(6,2) * t221 + Icges(6,6) * t266;
t147 = Icges(6,1) * t222 + Icges(6,4) * t221 + Icges(6,5) * t266;
t69 = t143 * t286 + t145 * t263 + t147 * t264;
t80 = t170 * t266 + t171 * t221 + t172 * t222;
t347 = t80 / 0.2e1 + t69 / 0.2e1;
t223 = t267 * t311 - t289 * t308;
t224 = t289 * t311 + t373;
t144 = Icges(6,5) * t224 + Icges(6,6) * t223 + Icges(6,3) * t268;
t146 = Icges(6,4) * t224 + Icges(6,2) * t223 + Icges(6,6) * t268;
t148 = Icges(6,1) * t224 + Icges(6,4) * t223 + Icges(6,5) * t268;
t70 = t144 * t286 + t146 * t263 + t148 * t264;
t81 = t170 * t268 + t171 * t223 + t172 * t224;
t346 = t81 / 0.2e1 + t70 / 0.2e1;
t150 = t224 * rSges(6,1) + t223 * rSges(6,2) + t268 * rSges(6,3);
t345 = -t150 + t356;
t173 = rSges(6,1) * t264 + rSges(6,2) * t263 + rSges(6,3) * t286;
t344 = -t173 + t351;
t343 = t307 * t226 + t358;
t342 = -t225 + t357;
t189 = t268 * rSges(4,1) - t267 * rSges(4,2) + t289 * rSges(4,3);
t274 = Icges(3,3) * t307 + (Icges(3,5) * t309 + Icges(3,6) * t312) * t306;
t275 = Icges(3,6) * t307 + (Icges(3,4) * t309 + Icges(3,2) * t312) * t306;
t276 = Icges(3,5) * t307 + (Icges(3,1) * t309 + Icges(3,4) * t312) * t306;
t340 = t306 * t309 * t276 + t307 * t274 + t275 * t370;
t242 = t290 * rSges(3,1) - t289 * rSges(3,2) + rSges(3,3) * t371;
t337 = -t370 / 0.2e1;
t336 = -t310 * pkin(1) + pkin(8) * t369;
t234 = rSges(4,1) * t286 - rSges(4,2) * t285 - rSges(4,3) * t370;
t291 = (pkin(2) * t309 - pkin(9) * t312) * t306;
t335 = t306 * (-t234 - t291);
t334 = t356 - t363;
t333 = t351 - t362;
t332 = t200 * t371 + t201 * t369 + t352;
t331 = t225 * t370 + t287 * t271 + t359;
t330 = t377 / 0.2e1 + t376 / 0.2e1 + t77 * t386 + t78 * t385 + t82;
t329 = t306 * (-t291 + t353);
t13 = t287 * t55 + t289 * t56 - t370 * t77;
t14 = t287 * t57 + t289 * t58 - t370 * t78;
t28 = t67 * t287 + t68 * t289 - t370 * t85;
t328 = t13 * t386 + t14 * t385 + t26 * t337 + t28 * t384 + t8 * t382 + t7 * t383;
t15 = t77 * t307 + (t310 * t56 - t313 * t55) * t306;
t16 = t78 * t307 + (t310 * t58 - t313 * t57) * t306;
t83 = t85 * t307;
t30 = t83 + (t68 * t310 - t67 * t313) * t306;
t327 = t15 * t386 + t16 * t385 + t26 * t381 + t30 * t384 + t8 * t371 / 0.2e1 - t7 * t369 / 0.2e1;
t326 = -t287 * rSges(5,1) - t265 * rSges(5,3);
t325 = -t222 * rSges(6,1) - t221 * rSges(6,2);
t323 = t249 + t350;
t322 = t306 * (-t291 + t344);
t321 = t225 * t371 + t226 * t369 + t332;
t320 = t306 * (-t291 + t333);
t319 = -t248 + t336;
t188 = rSges(4,1) * t266 - rSges(4,2) * t265 + rSges(4,3) * t287;
t318 = -t254 + t319;
t241 = rSges(3,1) * t288 - rSges(3,2) * t287 - rSges(3,3) * t369;
t317 = t201 + t323;
t174 = Icges(5,5) * t287 - Icges(5,6) * t266 + Icges(5,3) * t265;
t178 = Icges(5,4) * t287 - Icges(5,2) * t266 + Icges(5,6) * t265;
t182 = Icges(5,1) * t287 - Icges(5,4) * t266 + Icges(5,5) * t265;
t107 = t174 * t285 - t178 * t286 - t182 * t370;
t176 = Icges(4,5) * t266 - Icges(4,6) * t265 + Icges(4,3) * t287;
t180 = Icges(4,4) * t266 - Icges(4,2) * t265 + Icges(4,6) * t287;
t184 = Icges(4,1) * t266 - Icges(4,4) * t265 + Icges(4,5) * t287;
t109 = -t176 * t370 - t180 * t285 + t184 * t286;
t229 = -Icges(5,1) * t370 - Icges(5,4) * t286 + Icges(5,5) * t285;
t116 = t227 * t265 - t228 * t266 + t229 * t287;
t230 = Icges(4,5) * t286 - Icges(4,6) * t285 - Icges(4,3) * t370;
t118 = t230 * t287 - t231 * t265 + t232 * t266;
t316 = t67 / 0.2e1 + t77 / 0.2e1 + t118 / 0.2e1 + t116 / 0.2e1 + t109 / 0.2e1 + t107 / 0.2e1 + t347;
t175 = Icges(5,5) * t289 - Icges(5,6) * t268 + Icges(5,3) * t267;
t179 = Icges(5,4) * t289 - Icges(5,2) * t268 + Icges(5,6) * t267;
t183 = Icges(5,1) * t289 - Icges(5,4) * t268 + Icges(5,5) * t267;
t108 = t175 * t285 - t179 * t286 - t183 * t370;
t177 = Icges(4,5) * t268 - Icges(4,6) * t267 + Icges(4,3) * t289;
t181 = Icges(4,4) * t268 - Icges(4,2) * t267 + Icges(4,6) * t289;
t185 = Icges(4,1) * t268 - Icges(4,4) * t267 + Icges(4,5) * t289;
t110 = -t177 * t370 - t181 * t285 + t185 * t286;
t117 = t227 * t267 - t228 * t268 + t229 * t289;
t119 = t230 * t289 - t231 * t267 + t232 * t268;
t315 = t68 / 0.2e1 + t78 / 0.2e1 + t108 / 0.2e1 + t119 / 0.2e1 + t117 / 0.2e1 + t110 / 0.2e1 + t346;
t293 = rSges(2,1) * t313 - rSges(2,2) * t310;
t292 = -rSges(2,1) * t310 - rSges(2,2) * t313;
t277 = t307 * rSges(3,3) + (rSges(3,1) * t309 + rSges(3,2) * t312) * t306;
t240 = Icges(3,1) * t290 - Icges(3,4) * t289 + Icges(3,5) * t371;
t239 = Icges(3,1) * t288 - Icges(3,4) * t287 - Icges(3,5) * t369;
t238 = Icges(3,4) * t290 - Icges(3,2) * t289 + Icges(3,6) * t371;
t237 = Icges(3,4) * t288 - Icges(3,2) * t287 - Icges(3,6) * t369;
t236 = Icges(3,5) * t290 - Icges(3,6) * t289 + Icges(3,3) * t371;
t216 = t242 + t350;
t215 = -t241 + t336;
t194 = -t241 * t307 - t277 * t369;
t193 = t242 * t307 - t277 * t371;
t186 = -rSges(5,2) * t266 - t326;
t169 = t340 * t307;
t163 = (t241 * t310 + t242 * t313) * t306;
t162 = t274 * t371 - t275 * t289 + t276 * t290;
t161 = -t274 * t369 - t275 * t287 + t276 * t288;
t157 = t266 * t168;
t154 = t323 + t189;
t153 = -t188 + t319;
t149 = rSges(6,3) * t266 - t325;
t140 = -t189 * t370 - t234 * t289;
t139 = t188 * t370 + t234 * t287;
t132 = t307 * t236 + (t238 * t312 + t240 * t309) * t306;
t131 = t307 * t235 + (t237 * t312 + t239 * t309) * t306;
t130 = t286 * t142;
t129 = -t230 * t370 + t354;
t128 = -t229 * t370 + t355;
t127 = t268 * t141;
t126 = t129 * t307;
t125 = t128 * t307;
t124 = t317 + t187;
t123 = (rSges(5,2) - pkin(3)) * t266 + t318 + t326;
t122 = t188 * t289 - t189 * t287;
t121 = (-t188 - t248) * t307 + t313 * t335;
t120 = t307 * t189 + t310 * t335 + t246;
t115 = (t188 * t310 + t189 * t313) * t306 + t352;
t114 = t150 * t286 - t173 * t268;
t113 = -t149 * t286 + t173 * t266;
t112 = t289 * t353 + t361 * t370;
t111 = t186 * t370 + t233 * t287 + t359;
t106 = -t168 * t268 + t130;
t105 = -t141 * t286 + t157;
t104 = t226 + t317 + t150;
t103 = (-rSges(6,3) - pkin(3)) * t266 + t318 + t325 - t225;
t102 = (-t186 + t357) * t307 + t313 * t329;
t101 = t307 * t187 + t310 * t329 + t358;
t100 = t177 * t289 - t181 * t267 + t185 * t268;
t99 = t176 * t289 - t180 * t267 + t184 * t268;
t98 = t177 * t287 - t181 * t265 + t185 * t266;
t97 = t176 * t287 - t180 * t265 + t184 * t266;
t96 = t175 * t267 - t179 * t268 + t183 * t289;
t95 = t174 * t267 - t178 * t268 + t182 * t289;
t94 = t175 * t265 - t179 * t266 + t183 * t287;
t93 = t174 * t265 - t178 * t266 + t182 * t287;
t91 = t92 * t307;
t90 = t92 * t286;
t89 = t149 * t268 - t150 * t266;
t88 = t317 + t341 + t142;
t87 = -t348 - t287 * t301 + (-rSges(7,3) - pkin(3) + t314) * t266 + t318 + t324;
t86 = t289 * t186 + t287 * t361 + t190;
t84 = -t142 * t266 + t127;
t79 = (t186 * t310 + t187 * t313) * t306 + t332;
t74 = t289 * t344 + t345 * t370;
t73 = t149 * t370 + t173 * t287 + t331;
t72 = (-t149 + t342) * t307 + t313 * t322;
t71 = t307 * t150 + t310 * t322 + t343;
t66 = t286 * t152 - t268 * t362 + t130;
t65 = t266 * t191 - t286 * t364 + t157;
t64 = t144 * t268 + t146 * t223 + t148 * t224;
t63 = t143 * t268 + t145 * t223 + t147 * t224;
t62 = t144 * t266 + t146 * t221 + t148 * t222;
t61 = t143 * t266 + t145 * t221 + t147 * t222;
t54 = t289 * t149 + t287 * t345 + t360;
t53 = (t149 * t310 + t150 * t313) * t306 + t321;
t52 = t268 * t151 - t266 * t363 + t127;
t51 = t289 * t333 + t334 * t370;
t50 = t287 * t362 + t364 * t370 + t331;
t49 = (t342 - t364) * t307 + t313 * t320;
t48 = t307 * t363 + t310 * t320 + t343;
t47 = t126 + (-t109 * t313 + t110 * t310) * t306;
t46 = t125 + (-t107 * t313 + t108 * t310) * t306;
t45 = t109 * t287 + t110 * t289 - t129 * t370;
t44 = t107 * t287 + t108 * t289 - t128 * t370;
t43 = t119 * t307 + (t100 * t310 - t313 * t99) * t306;
t42 = t118 * t307 + (t310 * t98 - t313 * t97) * t306;
t41 = t117 * t307 + (t310 * t96 - t313 * t95) * t306;
t40 = t116 * t307 + (t310 * t94 - t313 * t93) * t306;
t39 = t100 * t289 - t119 * t370 + t287 * t99;
t38 = -t118 * t370 + t287 * t97 + t289 * t98;
t37 = -t117 * t370 + t287 * t95 + t289 * t96;
t36 = -t116 * t370 + t287 * t93 + t289 * t94;
t35 = t287 * t334 + t289 * t364 + t360;
t34 = (t310 * t364 + t313 * t363) * t306 + t321;
t33 = t91 + (t70 * t310 - t69 * t313) * t306;
t32 = t69 * t287 + t70 * t289 - t370 * t92;
t31 = t69 * t266 + t70 * t268 + t90;
t22 = t81 * t307 + (t310 * t64 - t313 * t63) * t306;
t21 = t80 * t307 + (t310 * t62 - t313 * t61) * t306;
t20 = t287 * t63 + t289 * t64 - t370 * t81;
t19 = t287 * t61 + t289 * t62 - t370 * t80;
t18 = t266 * t63 + t268 * t64 + t286 * t81;
t17 = t266 * t61 + t268 * t62 + t286 * t80;
t1 = [t340 + m(7) * (t87 ^ 2 + t88 ^ 2) + m(6) * (t103 ^ 2 + t104 ^ 2) + m(5) * (t123 ^ 2 + t124 ^ 2) + m(4) * (t153 ^ 2 + t154 ^ 2) + m(3) * (t215 ^ 2 + t216 ^ 2) + m(2) * (t292 ^ 2 + t293 ^ 2) + (-t229 - t230) * t370 + Icges(2,3) + t354 + t355 - t387; t91 + t83 + t126 + t125 + t169 + m(7) * (t48 * t88 + t49 * t87) + m(6) * (t103 * t72 + t104 * t71) + m(5) * (t101 * t124 + t102 * t123) + m(4) * (t120 * t154 + t121 * t153) + m(3) * (t193 * t216 + t194 * t215) + ((-t131 / 0.2e1 - t161 / 0.2e1 - t316) * t313 + (t132 / 0.2e1 + t162 / 0.2e1 + t315) * t310) * t306; (t30 + t33 + t46 + t47 + t169) * t307 + m(7) * (t34 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t53 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t101 ^ 2 + t102 ^ 2 + t79 ^ 2) + m(4) * (t115 ^ 2 + t120 ^ 2 + t121 ^ 2) + m(3) * (t163 ^ 2 + t193 ^ 2 + t194 ^ 2) + ((-t15 - t21 - t42 - t40 + (-t237 * t287 + t239 * t288 - t306 * t375) * t369) * t313 + (t16 + t22 + t43 + t41 + ((-t238 * t289 + t240 * t290 + (t236 * t310 - t375) * t306) * t310 + (t236 * t369 + t237 * t289 + t238 * t287 - t239 * t290 - t240 * t288) * t313) * t306) * t310 + ((-t131 - t161) * t313 + (t132 + t162) * t310) * t307) * t306; (-t128 - t129 + t387) * t370 + m(7) * (t50 * t87 + t51 * t88) + m(6) * (t103 * t73 + t104 * t74) + m(5) * (t111 * t123 + t112 * t124) + m(4) * (t139 * t153 + t140 * t154) + t315 * t289 + t316 * t287; (t28 / 0.2e1 + t32 / 0.2e1 + t44 / 0.2e1 + t45 / 0.2e1) * t307 + (t16 / 0.2e1 + t22 / 0.2e1 + t41 / 0.2e1 + t43 / 0.2e1) * t289 + (t15 / 0.2e1 + t21 / 0.2e1 + t40 / 0.2e1 + t42 / 0.2e1) * t287 + m(7) * (t34 * t35 + t48 * t51 + t49 * t50) + m(6) * (t54 * t53 + t71 * t74 + t72 * t73) + m(5) * (t101 * t112 + t102 * t111 + t79 * t86) + m(4) * (t115 * t122 + t120 * t140 + t121 * t139) + ((-t13 / 0.2e1 - t19 / 0.2e1 - t36 / 0.2e1 - t38 / 0.2e1) * t313 + (-t30 / 0.2e1 - t33 / 0.2e1 - t46 / 0.2e1 - t47 / 0.2e1) * t312 + (t14 / 0.2e1 + t20 / 0.2e1 + t37 / 0.2e1 + t39 / 0.2e1) * t310) * t306; (-t28 - t32 - t44 - t45) * t370 + (t14 + t20 + t37 + t39) * t289 + (t13 + t19 + t36 + t38) * t287 + m(7) * (t35 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t54 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t111 ^ 2 + t112 ^ 2 + t86 ^ 2) + m(4) * (t122 ^ 2 + t139 ^ 2 + t140 ^ 2); m(7) * (t265 * t88 + t267 * t87) + m(6) * (t103 * t267 + t104 * t265) + m(5) * (t123 * t267 + t124 * t265); m(7) * (t265 * t48 + t267 * t49 + t285 * t34) + m(6) * (t265 * t71 + t267 * t72 + t285 * t53) + m(5) * (t101 * t265 + t102 * t267 + t285 * t79); m(7) * (t265 * t51 + t267 * t50 + t285 * t35) + m(6) * (t265 * t74 + t267 * t73 + t285 * t54) + m(5) * (t111 * t267 + t112 * t265 + t285 * t86); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1) * (t265 ^ 2 + t267 ^ 2 + t285 ^ 2); t90 + t346 * t268 + t347 * t266 + m(7) * (t65 * t87 + t66 * t88) + m(6) * (t103 * t113 + t104 * t114) + t330; t31 * t381 + t22 * t385 + t21 * t386 + t33 * t384 + (t310 * t18 / 0.2e1 - t313 * t17 / 0.2e1) * t306 + m(7) * (t34 * t52 + t48 * t66 + t49 * t65) + m(6) * (t113 * t72 + t114 * t71 + t53 * t89) + t327; t31 * t337 + t17 * t383 + t32 * t384 + t20 * t385 + t19 * t386 + t18 * t382 + m(7) * (t35 * t52 + t50 * t65 + t51 * t66) + m(6) * (t113 * t73 + t114 * t74 + t54 * t89) + t328; m(6) * (t113 * t267 + t114 * t265 + t285 * t89) + m(7) * (t265 * t66 + t267 * t65 + t285 * t52); t266 * t17 + t268 * t18 + t286 * t31 + m(7) * (t52 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t113 ^ 2 + t114 ^ 2 + t89 ^ 2) + t349; m(7) * (t105 * t87 + t106 * t88) + t330; m(7) * (t105 * t49 + t106 * t48 + t34 * t84) + t327; m(7) * (t105 * t50 + t106 * t51 + t35 * t84) + t328; m(7) * (t105 * t267 + t106 * t265 + t285 * t84); m(7) * (t105 * t65 + t106 * t66 + t52 * t84) + t349; m(7) * (t105 ^ 2 + t106 ^ 2 + t84 ^ 2) + t349;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
