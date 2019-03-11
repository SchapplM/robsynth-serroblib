% Calculate joint inertia matrix for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:47
% EndTime: 2019-03-10 01:00:57
% DurationCPUTime: 4.76s
% Computational Cost: add. (13677->452), mult. (12686->648), div. (0->0), fcn. (13664->10), ass. (0->232)
t221 = qJ(2) + qJ(3);
t212 = qJ(4) + t221;
t208 = cos(t212);
t225 = cos(qJ(5));
t227 = cos(qJ(1));
t288 = t227 * t225;
t222 = sin(qJ(5));
t224 = sin(qJ(1));
t292 = t224 * t222;
t174 = t208 * t292 + t288;
t289 = t227 * t222;
t291 = t224 * t225;
t175 = t208 * t291 - t289;
t207 = sin(t212);
t297 = t207 * t224;
t96 = Icges(6,5) * t175 - Icges(6,6) * t174 + Icges(6,3) * t297;
t98 = Icges(7,4) * t175 + Icges(7,2) * t297 + Icges(7,6) * t174;
t344 = t96 + t98;
t176 = t208 * t289 - t291;
t177 = t208 * t288 + t292;
t295 = t207 * t227;
t97 = Icges(6,5) * t177 - Icges(6,6) * t176 + Icges(6,3) * t295;
t99 = Icges(7,4) * t177 + Icges(7,2) * t295 + Icges(7,6) * t176;
t343 = t97 + t99;
t100 = Icges(6,4) * t175 - Icges(6,2) * t174 + Icges(6,6) * t297;
t94 = Icges(7,5) * t175 + Icges(7,6) * t297 + Icges(7,3) * t174;
t342 = -t100 + t94;
t101 = Icges(6,4) * t177 - Icges(6,2) * t176 + Icges(6,6) * t295;
t95 = Icges(7,5) * t177 + Icges(7,6) * t295 + Icges(7,3) * t176;
t341 = -t101 + t95;
t102 = Icges(7,1) * t175 + Icges(7,4) * t297 + Icges(7,5) * t174;
t104 = Icges(6,1) * t175 - Icges(6,4) * t174 + Icges(6,5) * t297;
t340 = t102 + t104;
t103 = Icges(7,1) * t177 + Icges(7,4) * t295 + Icges(7,5) * t176;
t105 = Icges(6,1) * t177 - Icges(6,4) * t176 + Icges(6,5) * t295;
t339 = t103 + t105;
t325 = rSges(7,3) + qJ(6);
t329 = rSges(7,1) + pkin(5);
t338 = -t325 * t174 - t329 * t175;
t337 = t342 * t174 + t340 * t175 + t297 * t344;
t336 = t174 * t341 + t175 * t339 + t297 * t343;
t335 = t342 * t176 + t340 * t177 + t295 * t344;
t334 = t176 * t341 + t177 * t339 + t295 * t343;
t128 = -Icges(7,6) * t208 + (Icges(7,5) * t225 + Icges(7,3) * t222) * t207;
t130 = -Icges(7,2) * t208 + (Icges(7,4) * t225 + Icges(7,6) * t222) * t207;
t132 = -Icges(7,4) * t208 + (Icges(7,1) * t225 + Icges(7,5) * t222) * t207;
t56 = t174 * t128 + t130 * t297 + t175 * t132;
t129 = -Icges(6,3) * t208 + (Icges(6,5) * t225 - Icges(6,6) * t222) * t207;
t131 = -Icges(6,6) * t208 + (Icges(6,4) * t225 - Icges(6,2) * t222) * t207;
t133 = -Icges(6,5) * t208 + (Icges(6,1) * t225 - Icges(6,4) * t222) * t207;
t57 = t129 * t297 - t174 * t131 + t175 * t133;
t333 = -t56 - t57;
t58 = t176 * t128 + t130 * t295 + t177 * t132;
t59 = t129 * t295 - t176 * t131 + t177 * t133;
t332 = -t58 - t59;
t219 = t224 ^ 2;
t331 = t333 * t208 + (t224 * t337 + t227 * t336) * t207;
t330 = t332 * t208 + (t224 * t335 + t227 * t334) * t207;
t328 = t224 * pkin(7);
t327 = t224 * t336 - t227 * t337;
t326 = t224 * t334 - t227 * t335;
t287 = rSges(7,2) * t297 - t338;
t324 = rSges(7,2) * t295 + t325 * t176 + t177 * t329;
t323 = -t129 - t130;
t210 = sin(t221);
t211 = cos(t221);
t261 = rSges(4,1) * t211 - rSges(4,2) * t210;
t251 = Icges(4,5) * t211 - Icges(4,6) * t210;
t156 = -Icges(4,3) * t227 + t251 * t224;
t157 = Icges(4,3) * t224 + t251 * t227;
t220 = t227 ^ 2;
t301 = Icges(4,4) * t211;
t254 = -Icges(4,2) * t210 + t301;
t159 = Icges(4,6) * t224 + t254 * t227;
t302 = Icges(4,4) * t210;
t257 = Icges(4,1) * t211 - t302;
t161 = Icges(4,5) * t224 + t257 * t227;
t246 = -t159 * t210 + t161 * t211;
t158 = -Icges(4,6) * t227 + t254 * t224;
t160 = -Icges(4,5) * t227 + t257 * t224;
t247 = t158 * t210 - t160 * t211;
t250 = Icges(5,5) * t208 - Icges(5,6) * t207;
t146 = -Icges(5,3) * t227 + t250 * t224;
t147 = Icges(5,3) * t224 + t250 * t227;
t299 = Icges(5,4) * t208;
t253 = -Icges(5,2) * t207 + t299;
t149 = Icges(5,6) * t224 + t253 * t227;
t300 = Icges(5,4) * t207;
t256 = Icges(5,1) * t208 - t300;
t151 = Icges(5,5) * t224 + t256 * t227;
t248 = -t149 * t207 + t151 * t208;
t148 = -Icges(5,6) * t227 + t253 * t224;
t150 = -Icges(5,5) * t227 + t256 * t224;
t249 = t148 * t207 - t150 * t208;
t274 = -t220 * t146 - (t248 * t224 + (-t147 + t249) * t227) * t224 - t327;
t322 = -t220 * t156 - (t246 * t224 + (-t157 + t247) * t227) * t224 + t274;
t298 = t207 * t222;
t321 = t128 * t298 + (t132 + t133) * t207 * t225;
t228 = -pkin(8) - pkin(7);
t319 = t224 / 0.2e1;
t318 = -t227 / 0.2e1;
t223 = sin(qJ(2));
t317 = pkin(2) * t223;
t316 = pkin(3) * t210;
t315 = pkin(4) * t208;
t226 = cos(qJ(2));
t209 = t226 * pkin(2) + pkin(1);
t293 = t222 * t131;
t314 = -t207 * t293 + t208 * t323 + t321;
t313 = rSges(3,1) * t226;
t311 = rSges(3,2) * t223;
t309 = t227 * rSges(3,3);
t44 = -t208 * t98 + (t102 * t225 + t222 * t94) * t207;
t308 = t44 * t227;
t45 = -t208 * t99 + (t103 * t225 + t222 * t95) * t207;
t307 = t45 * t224;
t46 = -t208 * t96 + (-t100 * t222 + t104 * t225) * t207;
t306 = t46 * t227;
t47 = -t208 * t97 + (-t101 * t222 + t105 * t225) * t207;
t305 = t47 * t224;
t304 = Icges(3,4) * t223;
t303 = Icges(3,4) * t226;
t294 = t208 * t227;
t218 = -pkin(9) + t228;
t290 = t227 * t218;
t190 = pkin(3) * t211 + t209;
t184 = t227 * t190;
t203 = t227 * t209;
t285 = t227 * (t184 - t203) + (t190 - t209) * t219;
t239 = rSges(5,1) * t294 - rSges(5,2) * t295 + t224 * rSges(5,3);
t260 = rSges(5,1) * t208 - rSges(5,2) * t207;
t88 = t224 * (-t227 * rSges(5,3) + t260 * t224) + t227 * t239;
t283 = -t208 * rSges(7,2) + (t222 * t325 + t225 * t329) * t207;
t137 = -t208 * rSges(6,3) + (rSges(6,1) * t225 - rSges(6,2) * t222) * t207;
t183 = t207 * pkin(4) - t208 * pkin(10);
t282 = -t137 - t183;
t217 = t227 * pkin(7);
t281 = t224 * (t217 + (-pkin(1) + t209) * t224) + t227 * (-t227 * pkin(1) + t203 - t328);
t240 = t224 * rSges(4,3) + t227 * t261;
t93 = t224 * (-t227 * rSges(4,3) + t261 * t224) + t227 * t240;
t279 = pkin(4) * t294 + pkin(10) * t295;
t280 = t219 * (pkin(10) * t207 + t315) + t227 * t279;
t278 = t224 * rSges(3,3) + t227 * t313;
t276 = t219 + t220;
t275 = (t219 * t147 + (t249 * t227 + (-t146 + t248) * t224) * t227 + t326) * t224;
t273 = -t183 - t283;
t111 = t177 * rSges(6,1) - t176 * rSges(6,2) + rSges(6,3) * t295;
t189 = t210 * rSges(4,1) + t211 * rSges(4,2);
t270 = -t189 - t317;
t182 = t207 * rSges(5,1) + t208 * rSges(5,2);
t269 = -t182 - t316;
t268 = -t183 - t316;
t267 = -t190 - t315;
t266 = t224 * (t219 * t157 + (t247 * t227 + (-t156 + t246) * t224) * t227) + t275;
t265 = -t224 * t218 + t184;
t259 = -t175 * rSges(6,1) + t174 * rSges(6,2);
t109 = rSges(6,3) * t297 - t259;
t50 = t224 * t109 + t227 * t111 + t280;
t55 = t88 + t285;
t264 = -t137 + t268;
t263 = -t316 - t317;
t262 = -t311 + t313;
t258 = Icges(3,1) * t226 - t304;
t255 = -Icges(3,2) * t223 + t303;
t252 = Icges(3,5) * t226 - Icges(3,6) * t223;
t180 = Icges(5,2) * t208 + t300;
t181 = Icges(5,1) * t207 + t299;
t243 = -t180 * t207 + t181 * t208;
t187 = Icges(4,2) * t211 + t302;
t188 = Icges(4,1) * t210 + t301;
t242 = -t187 * t210 + t188 * t211;
t241 = t268 - t283;
t238 = -t182 + t263;
t237 = -t183 + t263;
t236 = t265 + t279;
t25 = t50 + t285;
t24 = t224 * t287 + t227 * t324 + t280;
t235 = -t137 + t237;
t234 = t274 * t227 + t275;
t233 = t237 - t283;
t22 = t24 + t285;
t232 = -(t307 - t308 + t305 - t306) * t208 / 0.2e1 + t330 * t319 + t331 * t318 + t327 * t297 / 0.2e1 + t326 * t295 / 0.2e1;
t231 = t322 * t227 + t266;
t179 = Icges(5,5) * t207 + Icges(5,6) * t208;
t230 = -t308 / 0.2e1 + t307 / 0.2e1 - t306 / 0.2e1 + t305 / 0.2e1 + (t208 * t149 + t207 * t151 + t224 * t179 + t243 * t227 - t332) * t319 + (t208 * t148 + t207 * t150 - t227 * t179 + t243 * t224 - t333) * t318;
t186 = Icges(4,5) * t210 + Icges(4,6) * t211;
t229 = t230 + (t211 * t159 + t210 * t161 + t224 * t186 + t242 * t227) * t319 + (t211 * t158 + t210 * t160 - t227 * t186 + t242 * t224) * t318;
t202 = t227 * rSges(2,1) - t224 * rSges(2,2);
t201 = -t224 * rSges(2,1) - t227 * rSges(2,2);
t200 = t223 * rSges(3,1) + t226 * rSges(3,2);
t169 = Icges(3,3) * t224 + t252 * t227;
t168 = -Icges(3,3) * t227 + t252 * t224;
t153 = t270 * t227;
t152 = t270 * t224;
t141 = t328 + (pkin(1) - t311) * t227 + t278;
t140 = t309 + t217 + (-pkin(1) - t262) * t224;
t139 = t269 * t227;
t138 = t269 * t224;
t127 = t238 * t227;
t126 = t238 * t224;
t125 = -t224 * t228 + t203 + t240;
t124 = (rSges(4,3) - t228) * t227 + (-t209 - t261) * t224;
t114 = t227 * (-t227 * t311 + t278) + (t262 * t224 - t309) * t224;
t113 = t239 + t265;
t112 = (rSges(5,3) - t218) * t227 + (-t190 - t260) * t224;
t107 = t282 * t227;
t106 = t282 * t224;
t87 = t264 * t227;
t86 = t264 * t224;
t81 = t235 * t227;
t80 = t235 * t224;
t75 = t273 * t227;
t74 = t273 * t224;
t73 = t241 * t227;
t72 = t241 * t224;
t71 = t233 * t227;
t70 = t233 * t224;
t69 = t236 + t111;
t68 = -t290 + ((-rSges(6,3) - pkin(10)) * t207 + t267) * t224 + t259;
t67 = -t208 * t111 - t137 * t295;
t66 = t208 * t109 + t137 * t297;
t65 = t93 + t281;
t62 = (t109 * t227 - t111 * t224) * t207;
t61 = t236 + t324;
t60 = -t290 + ((-rSges(7,2) - pkin(10)) * t207 + t267) * t224 + t338;
t49 = -t208 * t324 - t283 * t295;
t48 = t287 * t208 + t283 * t297;
t37 = t55 + t281;
t26 = (-t224 * t324 + t287 * t227) * t207;
t23 = t25 + t281;
t19 = t22 + t281;
t1 = [t211 * t187 + t210 * t188 + t226 * (Icges(3,2) * t226 + t304) + t223 * (Icges(3,1) * t223 + t303) + Icges(2,3) + (t181 - t293) * t207 + (t180 + t323) * t208 + m(7) * (t60 ^ 2 + t61 ^ 2) + m(6) * (t68 ^ 2 + t69 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2) + m(4) * (t124 ^ 2 + t125 ^ 2) + m(3) * (t140 ^ 2 + t141 ^ 2) + m(2) * (t201 ^ 2 + t202 ^ 2) + t321; t229 + (t220 / 0.2e1 + t219 / 0.2e1) * (Icges(3,5) * t223 + Icges(3,6) * t226) + m(4) * (t153 * t124 + t152 * t125) + m(5) * (t127 * t112 + t126 * t113) + m(6) * (t81 * t68 + t80 * t69) + m(7) * (t71 * t60 + t70 * t61) + m(3) * (-t140 * t227 - t141 * t224) * t200 + (t226 * (Icges(3,6) * t224 + t255 * t227) + t223 * (Icges(3,5) * t224 + t258 * t227)) * t319 + (t226 * (-Icges(3,6) * t227 + t255 * t224) + t223 * (-Icges(3,5) * t227 + t258 * t224)) * t318; m(7) * (t19 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t23 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t126 ^ 2 + t127 ^ 2 + t37 ^ 2) + m(4) * (t152 ^ 2 + t153 ^ 2 + t65 ^ 2) + t224 * t219 * t169 + m(3) * (t276 * t200 ^ 2 + t114 ^ 2) + t266 + (-t220 * t168 + (-t224 * t168 + t227 * t169) * t224 + t322) * t227; t229 + m(4) * (-t124 * t227 - t125 * t224) * t189 + m(7) * (t73 * t60 + t72 * t61) + m(6) * (t87 * t68 + t86 * t69) + m(5) * (t139 * t112 + t138 * t113); m(7) * (t22 * t19 + t72 * t70 + t73 * t71) + m(6) * (t25 * t23 + t86 * t80 + t87 * t81) + m(5) * (t138 * t126 + t139 * t127 + t55 * t37) + m(4) * (t93 * t65 + (-t152 * t224 - t153 * t227) * t189) + t231; m(7) * (t22 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(6) * (t25 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(5) * (t138 ^ 2 + t139 ^ 2 + t55 ^ 2) + m(4) * (t276 * t189 ^ 2 + t93 ^ 2) + t231; t230 + m(7) * (t75 * t60 + t74 * t61) + m(6) * (t106 * t69 + t107 * t68) + m(5) * (-t112 * t227 - t113 * t224) * t182; m(7) * (t24 * t19 + t74 * t70 + t75 * t71) + m(6) * (t106 * t80 + t107 * t81 + t50 * t23) + m(5) * (t88 * t37 + (-t126 * t224 - t127 * t227) * t182) + t234; m(7) * (t24 * t22 + t74 * t72 + t75 * t73) + m(6) * (t106 * t86 + t107 * t87 + t50 * t25) + m(5) * (t88 * t55 + (-t138 * t224 - t139 * t227) * t182) + t234; m(7) * (t24 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(6) * (t106 ^ 2 + t107 ^ 2 + t50 ^ 2) + m(5) * (t182 ^ 2 * t276 + t88 ^ 2) + t234; -t314 * t208 + m(7) * (t48 * t60 + t49 * t61) + m(6) * (t66 * t68 + t67 * t69) + ((t45 / 0.2e1 + t47 / 0.2e1 + t59 / 0.2e1 + t58 / 0.2e1) * t227 + (t46 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1 + t44 / 0.2e1) * t224) * t207; m(7) * (t26 * t19 + t48 * t71 + t49 * t70) + m(6) * (t62 * t23 + t66 * t81 + t67 * t80) + t232; m(7) * (t26 * t22 + t48 * t73 + t49 * t72) + m(6) * (t62 * t25 + t66 * t87 + t67 * t86) + t232; m(7) * (t26 * t24 + t48 * t75 + t49 * t74) + m(6) * (t67 * t106 + t66 * t107 + t62 * t50) + t232; m(7) * (t26 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t62 ^ 2 + t66 ^ 2 + t67 ^ 2) + t314 * t208 ^ 2 + (t330 * t227 + t331 * t224 + ((-t45 - t47) * t227 + (-t44 - t46) * t224) * t208) * t207; m(7) * (t174 * t61 + t176 * t60); m(7) * (t174 * t70 + t176 * t71 + t19 * t298); m(7) * (t174 * t72 + t176 * t73 + t22 * t298); m(7) * (t174 * t74 + t176 * t75 + t24 * t298); m(7) * (t174 * t49 + t176 * t48 + t26 * t298); m(7) * (t207 ^ 2 * t222 ^ 2 + t174 ^ 2 + t176 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
