% Calculate joint inertia matrix for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:32
% EndTime: 2019-03-09 13:50:43
% DurationCPUTime: 4.64s
% Computational Cost: add. (10763->441), mult. (15713->637), div. (0->0), fcn. (18442->10), ass. (0->220)
t351 = Icges(3,1) + Icges(4,1);
t348 = Icges(4,4) + Icges(3,5);
t211 = sin(qJ(2));
t350 = (Icges(3,4) - Icges(4,5)) * t211;
t212 = sin(qJ(1));
t216 = cos(qJ(1));
t285 = qJ(4) + qJ(5);
t200 = sin(t285);
t215 = cos(qJ(2));
t267 = cos(t285);
t168 = t211 * t200 + t215 * t267;
t146 = t168 * t216;
t209 = sin(qJ(6));
t213 = cos(qJ(6));
t126 = -t146 * t209 - t212 * t213;
t127 = t146 * t213 - t209 * t212;
t258 = t211 * t267;
t286 = t215 * t216;
t145 = t200 * t286 - t216 * t258;
t144 = t168 * t212;
t244 = -t144 * t209 + t213 * t216;
t231 = Icges(7,6) * t244;
t245 = t144 * t213 + t209 * t216;
t232 = Icges(7,5) * t245;
t169 = -t215 * t200 + t258;
t143 = t169 * t212;
t292 = Icges(7,3) * t143;
t220 = t231 + t232 - t292;
t221 = Icges(7,4) * t245 + Icges(7,2) * t244 - Icges(7,6) * t143;
t233 = Icges(7,4) * t244;
t338 = Icges(7,1) * t245;
t222 = -Icges(7,5) * t143 + t233 + t338;
t218 = t126 * t221 + t127 * t222 + t145 * t220;
t75 = Icges(7,5) * t127 + Icges(7,6) * t126 + Icges(7,3) * t145;
t76 = Icges(7,4) * t127 + Icges(7,2) * t126 + Icges(7,6) * t145;
t77 = Icges(7,1) * t127 + Icges(7,4) * t126 + Icges(7,5) * t145;
t28 = t126 * t76 + t127 * t77 + t145 * t75;
t13 = t212 * t28 - t218 * t216;
t293 = Icges(6,6) * t216;
t295 = Icges(6,2) * t143;
t301 = Icges(6,4) * t144;
t225 = t293 + t295 + t301;
t297 = Icges(6,5) * t216;
t305 = Icges(6,1) * t144;
t227 = Icges(6,4) * t143 + t297 + t305;
t349 = t212 * (Icges(6,5) * t146 - Icges(6,6) * t145 - Icges(6,3) * t212);
t98 = Icges(6,4) * t146 - Icges(6,2) * t145 - Icges(6,6) * t212;
t99 = Icges(6,1) * t146 - Icges(6,4) * t145 - Icges(6,5) * t212;
t321 = t13 + (-t145 * t98 + t146 * t99 - t349) * t212 - (t146 * t227 - t145 * t225 - t212 * (Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t216)) * t216;
t219 = Icges(7,2) * t244 ^ 2 - (0.2e1 * t231 + 0.2e1 * t232 - t292) * t143 + (0.2e1 * t233 + t338) * t245;
t27 = -t143 * t75 + t244 * t76 + t245 * t77;
t12 = t212 * t27 - t219 * t216;
t208 = t216 ^ 2;
t322 = (-t12 - (t143 * t98 + t144 * t99) * t212 + (-t349 + Icges(6,3) * t208 + (0.2e1 * t297 + t305) * t144 - (-0.2e1 * t293 - t295 - 0.2e1 * t301) * t143) * t216) * t216;
t237 = t321 * t212 + t322;
t347 = Icges(4,2) + Icges(3,3);
t214 = cos(qJ(4));
t210 = sin(qJ(4));
t290 = t210 * t215;
t174 = t211 * t214 - t290;
t166 = t174 * t216;
t291 = t210 * t211;
t239 = t214 * t215 + t291;
t167 = t239 * t216;
t105 = Icges(5,5) * t167 + Icges(5,6) * t166 - Icges(5,3) * t212;
t106 = Icges(5,4) * t167 + Icges(5,2) * t166 - Icges(5,6) * t212;
t107 = Icges(5,1) * t167 + Icges(5,4) * t166 - Icges(5,5) * t212;
t164 = t174 * t212;
t165 = t239 * t212;
t294 = Icges(5,6) * t216;
t296 = Icges(5,2) * t164;
t302 = Icges(5,4) * t165;
t226 = t294 + t296 + t302;
t298 = Icges(5,5) * t216;
t306 = Icges(5,1) * t165;
t228 = Icges(5,4) * t164 + t298 + t306;
t30 = (-t105 * t212 + t106 * t166 + t107 * t167) * t212 - (t167 * t228 + t166 * t226 - t212 * (Icges(5,5) * t165 + Icges(5,6) * t164 + Icges(5,3) * t216)) * t216;
t311 = t216 * ((t105 * t216 + t106 * t164 + t107 * t165) * t212 - (Icges(5,3) * t208 + (0.2e1 * t298 + t306) * t165 + (0.2e1 * t294 + t296 + 0.2e1 * t302) * t164) * t216);
t346 = -t212 * t30 - t237 + t311;
t345 = t348 * t215 + (-Icges(3,6) + Icges(4,6)) * t211;
t344 = t215 * t351 - t350;
t343 = -t212 / 0.2e1;
t342 = t212 / 0.2e1;
t341 = -t216 / 0.2e1;
t329 = t216 / 0.2e1;
t94 = rSges(7,3) * t168 + (rSges(7,1) * t213 - rSges(7,2) * t209) * t169;
t340 = pkin(5) * t169 + pkin(10) * t168 + t94;
t217 = -pkin(9) - pkin(8);
t199 = t216 * t217;
t195 = pkin(3) * t286;
t265 = -pkin(8) * t212 + t195;
t197 = pkin(4) * t214 + pkin(3);
t277 = pkin(4) * t291;
t270 = t197 * t286 + t212 * t217 + t216 * t277;
t323 = -pkin(3) + t197;
t327 = pkin(8) * t216;
t284 = t212 * (-t327 - t199 + (t323 * t215 + t277) * t212) + t216 * (-t265 + t270);
t79 = t127 * rSges(7,1) + t126 * rSges(7,2) + t145 * rSges(7,3);
t307 = t146 * pkin(5) + pkin(10) * t145 + t79;
t326 = t144 * pkin(5);
t223 = t245 * rSges(7,1) + t244 * rSges(7,2);
t78 = -t143 * rSges(7,3) + t223;
t35 = -(-t143 * pkin(10) + t326 + t78) * t212 - t307 * t216;
t24 = t35 - t284;
t337 = -t212 * t345 + t216 * t347;
t336 = t212 * t347 + t216 * t345;
t91 = Icges(7,3) * t168 + (Icges(7,5) * t213 - Icges(7,6) * t209) * t169;
t92 = Icges(7,6) * t168 + (Icges(7,4) * t213 - Icges(7,2) * t209) * t169;
t93 = Icges(7,5) * t168 + (Icges(7,1) * t213 - Icges(7,4) * t209) * t169;
t38 = -t143 * t91 + t244 * t92 + t245 * t93;
t3 = -t219 * t143 + t27 * t145 + t38 * t168;
t34 = t168 * t75 + (-t209 * t76 + t213 * t77) * t169;
t309 = t34 * t212;
t33 = t168 * t220 + (-t209 * t221 + t213 * t222) * t169;
t310 = t33 * t216;
t39 = t126 * t92 + t127 * t93 + t145 * t91;
t4 = -t218 * t143 + t28 * t145 + t39 * t168;
t263 = t4 * t343 - t168 * (t309 - t310) / 0.2e1 + t3 * t329 + t143 * t12 / 0.2e1 - t145 * t13 / 0.2e1;
t207 = t212 ^ 2;
t334 = m(4) / 0.2e1;
t333 = m(7) / 0.2e1;
t328 = -rSges(5,3) - pkin(8);
t320 = t169 * t213 * t93 + t168 * t91;
t316 = t209 * t92;
t313 = t216 * rSges(4,2);
t312 = t216 * rSges(3,3);
t70 = t340 * t212;
t71 = t340 * t216;
t303 = Icges(3,4) * t215;
t299 = Icges(4,5) * t215;
t289 = t211 * qJ(3);
t288 = t211 * t216;
t283 = t167 * rSges(5,1) + t166 * rSges(5,2);
t280 = pkin(2) * t286 + qJ(3) * t288;
t282 = t207 * (t215 * pkin(2) + t289) + t216 * t280;
t185 = pkin(2) * t211 - qJ(3) * t215;
t281 = -rSges(4,1) * t211 + rSges(4,3) * t215 - t185;
t279 = t216 * pkin(1) + t212 * pkin(7);
t278 = t207 + t208;
t276 = t33 / 0.2e1 + t38 / 0.2e1;
t275 = t39 / 0.2e1 + t34 / 0.2e1;
t129 = Icges(5,5) * t174 - Icges(5,6) * t239;
t130 = Icges(5,4) * t174 - Icges(5,2) * t239;
t131 = Icges(5,1) * t174 - Icges(5,4) * t239;
t274 = t129 * t329 + t130 * t164 / 0.2e1 + t131 * t165 / 0.2e1 - t239 * t226 / 0.2e1 + t174 * t228 / 0.2e1;
t273 = t129 * t342 - t130 * t166 / 0.2e1 - t131 * t167 / 0.2e1 + t106 * t239 / 0.2e1 - t107 * t174 / 0.2e1;
t147 = -pkin(4) * t290 + t323 * t211;
t137 = t212 * t147;
t65 = t137 + t70;
t138 = t216 * t147;
t66 = t138 + t71;
t269 = rSges(4,1) * t286 + t212 * rSges(4,2) + rSges(4,3) * t288;
t268 = t348 * t211 / 0.2e1 + (-Icges(4,6) / 0.2e1 + Icges(3,6) / 0.2e1) * t215;
t266 = -pkin(3) * t211 - t185;
t264 = t278 * t211;
t262 = t212 * (pkin(3) * t212 * t215 + t327) + t216 * t265 + t282;
t261 = t279 + t280;
t122 = rSges(6,1) * t169 - rSges(6,2) * t168;
t260 = -t122 + t266;
t132 = rSges(5,1) * t174 - rSges(5,2) * t239;
t259 = -t132 + t266;
t101 = t146 * rSges(6,1) - t145 * rSges(6,2) - rSges(6,3) * t212;
t257 = rSges(3,1) * t215 - rSges(3,2) * t211;
t256 = -rSges(5,1) * t165 - rSges(5,2) * t164;
t80 = t260 * t212 - t137;
t81 = t260 * t216 - t138;
t254 = t212 * t80 + t216 * t81;
t87 = t122 * t212 + t137;
t88 = t122 * t216 + t138;
t253 = t212 * t87 + t216 * t88;
t95 = t259 * t212;
t96 = t259 * t216;
t252 = t212 * t95 + t216 * t96;
t249 = -Icges(3,2) * t211 + t303;
t246 = Icges(4,3) * t211 + t299;
t100 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t216;
t69 = -t100 * t212 - t101 * t216;
t74 = t256 * t212 - t283 * t216;
t238 = rSges(3,1) * t286 - rSges(3,2) * t288 + t212 * rSges(3,3);
t204 = t216 * pkin(7);
t85 = t204 + t328 * t216 + (-t289 - pkin(1) + (-pkin(2) - pkin(3)) * t215) * t212 + t256;
t86 = t328 * t212 + t195 + t261 + t283;
t236 = m(5) * (t212 * t86 + t216 * t85);
t224 = t199 + t204 + (-pkin(1) + (-pkin(2) - t197) * t215 + (-pkin(4) * t210 - qJ(3)) * t211) * t212;
t72 = -t100 + t224;
t229 = t261 + t270;
t73 = t101 + t229;
t235 = m(6) * (t212 * t73 + t216 * t72);
t119 = Icges(6,5) * t169 - Icges(6,6) * t168;
t120 = Icges(6,4) * t169 - Icges(6,2) * t168;
t121 = Icges(6,1) * t169 - Icges(6,4) * t168;
t55 = t119 * t216 + t120 * t143 + t121 * t144;
t56 = -t119 * t212 - t120 * t145 + t121 * t146;
t61 = -t168 * t225 + t169 * t227;
t62 = -t168 * t98 + t169 * t99;
t230 = t310 / 0.2e1 - t309 / 0.2e1 + (t39 + t56 + t62) * t343 + (t38 + t55 + t61) * t329;
t45 = t69 - t284;
t189 = rSges(2,1) * t216 - rSges(2,2) * t212;
t188 = -rSges(2,1) * t212 - rSges(2,2) * t216;
t187 = rSges(3,1) * t211 + rSges(3,2) * t215;
t136 = t281 * t216;
t135 = t281 * t212;
t134 = t238 + t279;
t133 = t312 + t204 + (-pkin(1) - t257) * t212;
t113 = t261 + t269;
t112 = t313 + t204 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t215 + (-rSges(4,3) - qJ(3)) * t211) * t212;
t104 = t216 * t238 + (t257 * t212 - t312) * t212;
t82 = t216 * t269 + (-t313 + (rSges(4,1) * t215 + rSges(4,3) * t211) * t212) * t212 + t282;
t58 = t266 * t216 - t66;
t57 = t266 * t212 - t65;
t54 = -t74 + t262;
t49 = t229 + t307;
t48 = -t326 - (-pkin(10) - rSges(7,3)) * t143 - t223 + t224;
t47 = -t145 * t94 + t168 * t79;
t46 = -t143 * t94 - t168 * t78;
t42 = t143 * t79 + t145 * t78;
t41 = (-t169 * t316 + t320) * t168;
t40 = -t45 + t262;
t16 = t262 - t24;
t1 = [-t168 * t120 - t239 * t130 + t174 * t131 + Icges(2,3) + (t121 - t316) * t169 + m(7) * (t48 ^ 2 + t49 ^ 2) + m(6) * (t72 ^ 2 + t73 ^ 2) + m(5) * (t85 ^ 2 + t86 ^ 2) + m(4) * (t112 ^ 2 + t113 ^ 2) + m(3) * (t133 ^ 2 + t134 ^ 2) + m(2) * (t188 ^ 2 + t189 ^ 2) + t320 + ((Icges(3,2) + Icges(4,3)) * t215 + t350) * t215 + (t211 * t351 - t299 + t303) * t211; (-t61 / 0.2e1 - t55 / 0.2e1 + t268 * t216 + (Icges(3,6) * t329 + Icges(4,6) * t341 + t246 * t342 + t249 * t343) * t215 + (t329 * t348 + t343 * t344) * t211 - t274 - t276) * t216 + (t62 / 0.2e1 + t56 / 0.2e1 + t268 * t212 + (Icges(3,6) * t342 + Icges(4,6) * t343 + t246 * t341 + t249 * t329) * t215 + (t329 * t344 + t342 * t348) * t211 - t273 + t275) * t212 + m(7) * (t48 * t58 + t49 * t57) + m(6) * (t72 * t81 + t73 * t80) + m(5) * (t85 * t96 + t86 * t95) + m(4) * (t112 * t136 + t113 * t135) + m(3) * (-t133 * t216 - t134 * t212) * t187; m(7) * (t16 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t40 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t54 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(4) * (t135 ^ 2 + t136 ^ 2 + t82 ^ 2) + m(3) * (t278 * t187 ^ 2 + t104 ^ 2) + t336 * t212 * t207 + (t337 * t208 + (t337 * t212 + t336 * t216) * t212) * t216 - t346; 0.2e1 * ((t212 * t49 + t216 * t48) * t333 + t235 / 0.2e1 + t236 / 0.2e1 + (t112 * t216 + t113 * t212) * t334) * t211; m(7) * (-t16 * t215 + (t212 * t57 + t216 * t58) * t211) + m(6) * (t254 * t211 - t215 * t40) + m(5) * (t252 * t211 - t215 * t54) + m(4) * (-t215 * t82 + (t135 * t212 + t136 * t216) * t211); 0.2e1 * (t334 + m(5) / 0.2e1 + m(6) / 0.2e1 + t333) * (t278 * t211 ^ 2 + t215 ^ 2); t274 * t216 + t273 * t212 + m(7) * (t48 * t66 + t49 * t65) + m(6) * (t72 * t88 + t73 * t87) + t132 * t236 + t230; m(7) * (t16 * t24 + t57 * t65 + t58 * t66) + m(6) * (t40 * t45 + t80 * t87 + t81 * t88) + m(5) * (t252 * t132 + t54 * t74) + t346; m(5) * (t132 * t264 - t215 * t74) + m(6) * (t253 * t211 - t215 * t45) + m(7) * (-t215 * t24 + (t212 * t65 + t216 * t66) * t211); -t311 + (t30 + t321) * t212 + m(7) * (t24 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t45 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t278 * t132 ^ 2 + t74 ^ 2) + t322; m(7) * (t48 * t71 + t49 * t70) + t122 * t235 + t230; m(7) * (t16 * t35 + t57 * t70 + t58 * t71) + m(6) * (t254 * t122 + t40 * t69) - t237; m(6) * (t122 * t264 - t215 * t69) + m(7) * (-t215 * t35 + (t212 * t70 + t216 * t71) * t211); m(7) * (t24 * t35 + t65 * t70 + t66 * t71) + m(6) * (t253 * t122 + t45 * t69) + t237; m(6) * (t278 * t122 ^ 2 + t69 ^ 2) + m(7) * (t35 ^ 2 + t70 ^ 2 + t71 ^ 2) + t237; m(7) * (t46 * t48 + t47 * t49) + t41 + t275 * t145 - t276 * t143; m(7) * (t16 * t42 + t46 * t58 + t47 * t57) - t263; m(7) * (-t215 * t42 + (t212 * t47 + t216 * t46) * t211); m(7) * (t24 * t42 + t46 * t66 + t47 * t65) + t263; m(7) * (t35 * t42 + t46 * t71 + t47 * t70) + t263; t145 * t4 - t143 * t3 + t168 * (-t33 * t143 + t34 * t145 + t41) + m(7) * (t42 ^ 2 + t46 ^ 2 + t47 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
