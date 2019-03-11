% Calculate joint inertia matrix for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:37
% EndTime: 2019-03-09 16:43:51
% DurationCPUTime: 6.12s
% Computational Cost: add. (7197->415), mult. (9428->600), div. (0->0), fcn. (10025->8), ass. (0->205)
t353 = Icges(4,4) + Icges(5,6);
t352 = Icges(4,1) + Icges(5,2);
t351 = -Icges(4,2) - Icges(5,3);
t213 = qJ(2) + qJ(3);
t204 = cos(t213);
t350 = t353 * t204;
t203 = sin(t213);
t349 = t353 * t203;
t348 = Icges(5,4) - Icges(4,5);
t347 = Icges(5,5) - Icges(4,6);
t346 = t351 * t203 + t350;
t345 = t352 * t204 - t349;
t217 = cos(qJ(5));
t219 = cos(qJ(1));
t275 = t217 * t219;
t214 = sin(qJ(5));
t216 = sin(qJ(1));
t279 = t214 * t216;
t166 = -t203 * t275 + t279;
t277 = t216 * t217;
t278 = t214 * t219;
t167 = t203 * t278 + t277;
t280 = t204 * t219;
t91 = Icges(7,5) * t167 + Icges(7,6) * t280 + Icges(7,3) * t166;
t97 = Icges(6,4) * t167 - Icges(6,2) * t166 + Icges(6,6) * t280;
t344 = t91 - t97;
t168 = t203 * t277 + t278;
t170 = t203 * t279 - t275;
t282 = t204 * t216;
t92 = Icges(7,5) * t170 + Icges(7,6) * t282 - Icges(7,3) * t168;
t98 = Icges(6,4) * t170 + Icges(6,2) * t168 + Icges(6,6) * t282;
t343 = t92 - t98;
t93 = Icges(6,5) * t167 - Icges(6,6) * t166 + Icges(6,3) * t280;
t95 = Icges(7,4) * t167 + Icges(7,2) * t280 + Icges(7,6) * t166;
t342 = t93 + t95;
t94 = Icges(6,5) * t170 + Icges(6,6) * t168 + Icges(6,3) * t282;
t96 = Icges(7,4) * t170 + Icges(7,2) * t282 - Icges(7,6) * t168;
t341 = t94 + t96;
t101 = Icges(6,1) * t167 - Icges(6,4) * t166 + Icges(6,5) * t280;
t99 = Icges(7,1) * t167 + Icges(7,4) * t280 + Icges(7,5) * t166;
t340 = t101 + t99;
t100 = Icges(7,1) * t170 + Icges(7,4) * t282 - Icges(7,5) * t168;
t102 = Icges(6,1) * t170 + Icges(6,4) * t168 + Icges(6,5) * t282;
t339 = t100 + t102;
t312 = rSges(7,3) + qJ(6);
t330 = rSges(7,1) + pkin(5);
t338 = t312 * t168 - t170 * t330;
t336 = Icges(5,1) + Icges(4,3);
t335 = t216 * t346 + t219 * t347;
t334 = -t216 * t347 + t219 * t346;
t333 = t345 * t216 + t219 * t348;
t332 = -t216 * t348 + t345 * t219;
t331 = t347 * t203 - t204 * t348;
t329 = t344 * t166 + t340 * t167 + t342 * t280;
t328 = t343 * t166 + t339 * t167 + t341 * t280;
t327 = -t344 * t168 + t340 * t170 + t342 * t282;
t326 = -t343 * t168 + t339 * t170 + t341 * t282;
t129 = Icges(7,6) * t203 + (-Icges(7,5) * t214 + Icges(7,3) * t217) * t204;
t131 = Icges(7,2) * t203 + (-Icges(7,4) * t214 + Icges(7,6) * t217) * t204;
t133 = Icges(7,4) * t203 + (-Icges(7,1) * t214 + Icges(7,5) * t217) * t204;
t55 = t129 * t166 + t131 * t280 + t133 * t167;
t130 = Icges(6,3) * t203 + (-Icges(6,5) * t214 - Icges(6,6) * t217) * t204;
t132 = Icges(6,6) * t203 + (-Icges(6,4) * t214 - Icges(6,2) * t217) * t204;
t134 = Icges(6,5) * t203 + (-Icges(6,1) * t214 - Icges(6,4) * t217) * t204;
t56 = t130 * t280 - t132 * t166 + t134 * t167;
t325 = t55 + t56;
t57 = -t129 * t168 + t131 * t282 + t133 * t170;
t58 = t130 * t282 + t132 * t168 + t134 * t170;
t324 = t57 + t58;
t323 = t351 * t204 - t349;
t322 = t352 * t203 + t350;
t321 = -t331 * t216 + t336 * t219;
t320 = t336 * t216 + t331 * t219;
t319 = t335 * t203 - t333 * t204;
t318 = t334 * t203 - t332 * t204;
t317 = (t216 * t328 + t219 * t329) * t204 + t325 * t203;
t316 = (t216 * t326 + t219 * t327) * t204 + t324 * t203;
t314 = t216 * t329 - t219 * t328;
t313 = t216 * t327 - t219 * t326;
t273 = rSges(7,2) * t280 + t312 * t166 + t167 * t330;
t272 = rSges(7,2) * t282 - t338;
t311 = -t203 * t348 - t347 * t204;
t310 = t203 * t323 + t204 * t322;
t281 = t204 * t217;
t309 = t129 * t281 + (t130 + t131) * t203;
t308 = -t217 * t132 + (-t133 - t134) * t214;
t212 = t219 ^ 2;
t307 = -t313 + t321 * t212 + (t318 * t216 + (-t319 + t320) * t219) * t216;
t211 = t216 ^ 2;
t306 = m(5) / 0.2e1;
t305 = m(6) / 0.2e1;
t304 = m(7) / 0.2e1;
t303 = -pkin(3) - pkin(9);
t301 = t216 / 0.2e1;
t300 = -t219 / 0.2e1;
t215 = sin(qJ(2));
t299 = pkin(2) * t215;
t298 = (t204 * t308 + t309) * t203;
t218 = cos(qJ(2));
t297 = rSges(3,1) * t218;
t296 = rSges(3,2) * t215;
t295 = t219 * rSges(3,3);
t42 = t203 * t95 + (-t214 * t99 + t217 * t91) * t204;
t294 = t42 * t216;
t43 = t203 * t96 + (-t100 * t214 + t217 * t92) * t204;
t293 = t43 * t219;
t44 = t203 * t93 + (-t101 * t214 - t217 * t97) * t204;
t292 = t44 * t216;
t45 = t203 * t94 + (-t102 * t214 - t217 * t98) * t204;
t291 = t45 * t219;
t290 = Icges(3,4) * t215;
t289 = Icges(3,4) * t218;
t284 = qJ(4) * t203;
t283 = t203 * t219;
t220 = -pkin(8) - pkin(7);
t274 = t219 * t220;
t200 = pkin(2) * t218 + pkin(1);
t192 = t219 * t200;
t209 = t219 * pkin(7);
t270 = t216 * (t274 + t209 + (-pkin(1) + t200) * t216) + t219 * (-t219 * pkin(1) + t192 + (-pkin(7) - t220) * t216);
t230 = rSges(4,1) * t280 - rSges(4,2) * t283 + t216 * rSges(4,3);
t250 = rSges(4,1) * t204 - rSges(4,2) * t203;
t86 = t216 * (-t219 * rSges(4,3) + t216 * t250) + t219 * t230;
t269 = t203 * rSges(7,2) + (-t214 * t330 + t217 * t312) * t204;
t266 = pkin(3) * t280 + qJ(4) * t283;
t268 = t211 * (pkin(3) * t204 + t284) + t219 * t266;
t180 = pkin(3) * t203 - qJ(4) * t204;
t267 = rSges(5,2) * t203 + rSges(5,3) * t204 - t180;
t265 = t216 * rSges(3,3) + t219 * t297;
t264 = t216 * pkin(4) + pkin(9) * t280;
t263 = t211 + t212;
t104 = t167 * rSges(6,1) - t166 * rSges(6,2) + rSges(6,3) * t280;
t182 = rSges(4,1) * t203 + rSges(4,2) * t204;
t260 = -t182 - t299;
t259 = -pkin(9) * t203 - t180;
t258 = (t320 * t211 + ((-t318 + t321) * t216 + t319 * t219) * t219 + t314) * t216;
t257 = -t200 - t284;
t256 = -t216 * t220 + t192;
t210 = t219 * pkin(4);
t255 = t210 - t274;
t229 = t216 * rSges(5,1) - rSges(5,2) * t280 + rSges(5,3) * t283;
t63 = t216 * (-t219 * rSges(5,1) + (-rSges(5,2) * t204 + rSges(5,3) * t203) * t216) + t219 * t229 + t268;
t254 = t216 * (pkin(9) * t282 - t210) + t219 * t264 + t268;
t253 = t267 - t299;
t136 = t203 * rSges(6,3) + (-rSges(6,1) * t214 - rSges(6,2) * t217) * t204;
t252 = -t136 + t259;
t251 = -t296 + t297;
t249 = -t170 * rSges(6,1) - t168 * rSges(6,2);
t248 = Icges(3,1) * t218 - t290;
t246 = -Icges(3,2) * t215 + t289;
t243 = Icges(3,5) * t218 - Icges(3,6) * t215;
t231 = t259 - t269;
t228 = t259 - t299;
t227 = t256 + t266;
t106 = rSges(6,3) * t282 - t249;
t37 = t219 * t104 + t216 * t106 + t254;
t226 = -t136 + t228;
t225 = t228 - t269;
t224 = t227 + t264;
t22 = t216 * t272 + t219 * t273 + t254;
t223 = (-t293 + t294 - t291 + t292) * t203 / 0.2e1 + t317 * t301 + t316 * t300 + t313 * t282 / 0.2e1 + t314 * t280 / 0.2e1;
t222 = t307 * t219 + t258;
t221 = t294 / 0.2e1 - t293 / 0.2e1 + t292 / 0.2e1 - t291 / 0.2e1 + (t332 * t203 + t334 * t204 + t311 * t216 + t310 * t219 + t325) * t301 + (t333 * t203 + t335 * t204 + t310 * t216 - t311 * t219 + t324) * t300;
t202 = t204 ^ 2;
t190 = rSges(2,1) * t219 - rSges(2,2) * t216;
t189 = -rSges(2,1) * t216 - rSges(2,2) * t219;
t188 = rSges(3,1) * t215 + rSges(3,2) * t218;
t155 = Icges(3,3) * t216 + t219 * t243;
t154 = -Icges(3,3) * t219 + t216 * t243;
t138 = t260 * t219;
t137 = t260 * t216;
t122 = t216 * pkin(7) + (pkin(1) - t296) * t219 + t265;
t121 = t295 + t209 + (-pkin(1) - t251) * t216;
t117 = t267 * t219;
t116 = t267 * t216;
t115 = t230 + t256;
t114 = (rSges(4,3) - t220) * t219 + (-t200 - t250) * t216;
t111 = t253 * t219;
t110 = t253 * t216;
t107 = t219 * (-t219 * t296 + t265) + (t216 * t251 - t295) * t216;
t85 = t252 * t219;
t84 = t252 * t216;
t83 = t227 + t229;
t82 = (rSges(5,1) - t220) * t219 + (-t200 + (rSges(5,2) - pkin(3)) * t204 + (-rSges(5,3) - qJ(4)) * t203) * t216;
t81 = t226 * t219;
t80 = t226 * t216;
t71 = t231 * t219;
t70 = t231 * t216;
t69 = t225 * t219;
t68 = t225 * t216;
t67 = t104 * t203 - t136 * t280;
t66 = -t106 * t203 + t136 * t282;
t65 = t224 + t104;
t64 = ((-rSges(6,3) + t303) * t204 + t257) * t216 + t249 + t255;
t60 = t86 + t270;
t59 = (-t104 * t216 + t106 * t219) * t204;
t50 = t224 + t273;
t49 = ((-rSges(7,2) + t303) * t204 + t257) * t216 + t255 + t338;
t48 = t63 + t270;
t47 = t203 * t273 - t269 * t280;
t46 = -t203 * t272 + t269 * t282;
t24 = (-t216 * t273 + t219 * t272) * t204;
t23 = t37 + t270;
t21 = t22 + t270;
t1 = [t218 * (Icges(3,2) * t218 + t290) + t215 * (Icges(3,1) * t215 + t289) + Icges(2,3) + t322 * t203 + (t308 - t323) * t204 + m(7) * (t49 ^ 2 + t50 ^ 2) + m(6) * (t64 ^ 2 + t65 ^ 2) + m(5) * (t82 ^ 2 + t83 ^ 2) + m(4) * (t114 ^ 2 + t115 ^ 2) + m(3) * (t121 ^ 2 + t122 ^ 2) + m(2) * (t189 ^ 2 + t190 ^ 2) + t309; m(3) * (-t121 * t219 - t122 * t216) * t188 + m(7) * (t49 * t69 + t50 * t68) + m(6) * (t64 * t81 + t65 * t80) + m(5) * (t110 * t83 + t111 * t82) + m(4) * (t114 * t138 + t115 * t137) + (t212 / 0.2e1 + t211 / 0.2e1) * (Icges(3,5) * t215 + Icges(3,6) * t218) + (t218 * (Icges(3,6) * t216 + t219 * t246) + t215 * (Icges(3,5) * t216 + t219 * t248)) * t301 + (t218 * (-Icges(3,6) * t219 + t216 * t246) + t215 * (-Icges(3,5) * t219 + t216 * t248)) * t300 + t221; m(7) * (t21 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t23 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t110 ^ 2 + t111 ^ 2 + t48 ^ 2) + m(4) * (t137 ^ 2 + t138 ^ 2 + t60 ^ 2) + m(3) * (t188 ^ 2 * t263 + t107 ^ 2) + t216 * t211 * t155 + t258 + (-t212 * t154 + (-t216 * t154 + t219 * t155) * t216 + t307) * t219; m(7) * (t49 * t71 + t50 * t70) + m(6) * (t64 * t85 + t65 * t84) + m(5) * (t116 * t83 + t117 * t82) + m(4) * (-t114 * t219 - t115 * t216) * t182 + t221; m(7) * (t21 * t22 + t68 * t70 + t69 * t71) + m(6) * (t23 * t37 + t80 * t84 + t81 * t85) + m(5) * (t110 * t116 + t111 * t117 + t48 * t63) + m(4) * (t86 * t60 + (-t137 * t216 - t138 * t219) * t182) + t222; m(7) * (t22 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t37 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(5) * (t116 ^ 2 + t117 ^ 2 + t63 ^ 2) + m(4) * (t182 ^ 2 * t263 + t86 ^ 2) + t222; 0.2e1 * ((t216 * t50 + t219 * t49) * t304 + (t216 * t65 + t219 * t64) * t305 + (t216 * t83 + t219 * t82) * t306) * t203; m(7) * (-t204 * t21 + (t216 * t68 + t219 * t69) * t203) + m(6) * (-t204 * t23 + (t216 * t80 + t219 * t81) * t203) + m(5) * (-t204 * t48 + (t110 * t216 + t111 * t219) * t203); m(7) * (-t204 * t22 + (t216 * t70 + t219 * t71) * t203) + m(6) * (-t204 * t37 + (t216 * t84 + t219 * t85) * t203) + m(5) * (-t204 * t63 + (t116 * t216 + t117 * t219) * t203); 0.2e1 * (t306 + t305 + t304) * (t203 ^ 2 * t263 + t202); m(7) * (t46 * t49 + t47 * t50) + m(6) * (t64 * t66 + t65 * t67) + ((t44 / 0.2e1 + t42 / 0.2e1 + t56 / 0.2e1 + t55 / 0.2e1) * t219 + (t57 / 0.2e1 + t58 / 0.2e1 + t45 / 0.2e1 + t43 / 0.2e1) * t216) * t204 + t298; m(7) * (t21 * t24 + t46 * t69 + t47 * t68) + m(6) * (t23 * t59 + t66 * t81 + t67 * t80) + t223; m(7) * (t22 * t24 + t46 * t71 + t47 * t70) + m(6) * (t37 * t59 + t66 * t85 + t67 * t84) + t223; m(6) * (-t59 * t204 + (t216 * t67 + t219 * t66) * t203) + m(7) * (-t24 * t204 + (t216 * t47 + t219 * t46) * t203); t298 * t203 + m(7) * (t24 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t59 ^ 2 + t66 ^ 2 + t67 ^ 2) + (t317 * t219 + t316 * t216 + ((t42 + t44) * t219 + (t43 + t45) * t216) * t203) * t204; m(7) * (t166 * t49 - t168 * t50); m(7) * (t166 * t69 - t168 * t68 + t21 * t281); m(7) * (t166 * t71 - t168 * t70 + t22 * t281); m(7) * (-t202 * t217 + (t166 * t219 - t168 * t216) * t203); m(7) * (t166 * t46 - t168 * t47 + t24 * t281); m(7) * (t202 * t217 ^ 2 + t166 ^ 2 + t168 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
