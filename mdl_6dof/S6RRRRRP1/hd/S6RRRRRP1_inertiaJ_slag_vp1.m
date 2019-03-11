% Calculate joint inertia matrix for
% S6RRRRRP1
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:36
% EndTime: 2019-03-10 00:55:46
% DurationCPUTime: 5.02s
% Computational Cost: add. (14353->462), mult. (12802->658), div. (0->0), fcn. (13639->10), ass. (0->233)
t222 = qJ(2) + qJ(3);
t213 = qJ(4) + t222;
t208 = cos(t213);
t227 = cos(qJ(5));
t229 = cos(qJ(1));
t287 = t229 * t227;
t224 = sin(qJ(5));
t226 = sin(qJ(1));
t290 = t226 * t224;
t172 = -t208 * t290 - t287;
t288 = t229 * t224;
t289 = t226 * t227;
t173 = t208 * t289 - t288;
t207 = sin(t213);
t294 = t207 * t226;
t100 = Icges(6,5) * t173 + Icges(6,6) * t172 + Icges(6,3) * t294;
t98 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t294;
t346 = t100 + t98;
t174 = -t208 * t288 + t289;
t175 = t208 * t287 + t290;
t292 = t207 * t229;
t101 = Icges(6,5) * t175 + Icges(6,6) * t174 + Icges(6,3) * t292;
t99 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t292;
t345 = t101 + t99;
t102 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t294;
t104 = Icges(6,4) * t173 + Icges(6,2) * t172 + Icges(6,6) * t294;
t344 = t102 + t104;
t103 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t292;
t105 = Icges(6,4) * t175 + Icges(6,2) * t174 + Icges(6,6) * t292;
t343 = t103 + t105;
t106 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t294;
t108 = Icges(6,1) * t173 + Icges(6,4) * t172 + Icges(6,5) * t294;
t342 = t106 + t108;
t107 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t292;
t109 = Icges(6,1) * t175 + Icges(6,4) * t174 + Icges(6,5) * t292;
t341 = t107 + t109;
t223 = -qJ(6) - pkin(10);
t340 = rSges(7,3) - t223;
t339 = t344 * t172 + t342 * t173 + t346 * t294;
t338 = t343 * t172 + t341 * t173 + t345 * t294;
t337 = t344 * t174 + t342 * t175 + t346 * t292;
t336 = t343 * t174 + t341 * t175 + t345 * t292;
t128 = -Icges(7,3) * t208 + (Icges(7,5) * t227 - Icges(7,6) * t224) * t207;
t130 = -Icges(7,6) * t208 + (Icges(7,4) * t227 - Icges(7,2) * t224) * t207;
t132 = -Icges(7,5) * t208 + (Icges(7,1) * t227 - Icges(7,4) * t224) * t207;
t56 = t128 * t294 + t172 * t130 + t173 * t132;
t129 = -Icges(6,3) * t208 + (Icges(6,5) * t227 - Icges(6,6) * t224) * t207;
t131 = -Icges(6,6) * t208 + (Icges(6,4) * t227 - Icges(6,2) * t224) * t207;
t133 = -Icges(6,5) * t208 + (Icges(6,1) * t227 - Icges(6,4) * t224) * t207;
t57 = t129 * t294 + t172 * t131 + t173 * t133;
t335 = -t56 - t57;
t58 = t128 * t292 + t174 * t130 + t175 * t132;
t59 = t129 * t292 + t174 * t131 + t175 * t133;
t334 = -t58 - t59;
t209 = t227 * pkin(5) + pkin(4);
t291 = t208 * t229;
t333 = t175 * rSges(7,1) + t174 * rSges(7,2) + pkin(5) * t290 + t209 * t291 + t340 * t292;
t220 = t226 ^ 2;
t332 = t335 * t208 + (t226 * t339 + t338 * t229) * t207;
t331 = t334 * t208 + (t226 * t337 + t229 * t336) * t207;
t330 = t226 * pkin(7);
t329 = t338 * t226 - t229 * t339;
t328 = t226 * t336 - t229 * t337;
t261 = -t173 * rSges(7,1) - t172 * rSges(7,2);
t314 = pkin(10) + t223;
t315 = -pkin(4) + t209;
t303 = rSges(7,3) * t294 - t261 - pkin(5) * t288 + (-t314 * t207 + t315 * t208) * t226;
t281 = pkin(4) * t291 + pkin(10) * t292;
t327 = -t281 + t333;
t326 = -t128 - t129;
t211 = sin(t222);
t212 = cos(t222);
t264 = rSges(4,1) * t212 - rSges(4,2) * t211;
t325 = (-t130 - t131) * t224;
t253 = Icges(4,5) * t212 - Icges(4,6) * t211;
t156 = -Icges(4,3) * t229 + t253 * t226;
t157 = Icges(4,3) * t226 + t253 * t229;
t221 = t229 ^ 2;
t298 = Icges(4,4) * t212;
t256 = -Icges(4,2) * t211 + t298;
t159 = Icges(4,6) * t226 + t256 * t229;
t299 = Icges(4,4) * t211;
t259 = Icges(4,1) * t212 - t299;
t161 = Icges(4,5) * t226 + t259 * t229;
t248 = -t159 * t211 + t161 * t212;
t158 = -Icges(4,6) * t229 + t256 * t226;
t160 = -Icges(4,5) * t229 + t259 * t226;
t249 = t158 * t211 - t160 * t212;
t252 = Icges(5,5) * t208 - Icges(5,6) * t207;
t146 = -Icges(5,3) * t229 + t252 * t226;
t147 = Icges(5,3) * t226 + t252 * t229;
t296 = Icges(5,4) * t208;
t255 = -Icges(5,2) * t207 + t296;
t149 = Icges(5,6) * t226 + t255 * t229;
t297 = Icges(5,4) * t207;
t258 = Icges(5,1) * t208 - t297;
t151 = Icges(5,5) * t226 + t258 * t229;
t250 = -t149 * t207 + t151 * t208;
t148 = -Icges(5,6) * t229 + t255 * t226;
t150 = -Icges(5,5) * t229 + t258 * t226;
t251 = t148 * t207 - t150 * t208;
t276 = -t221 * t146 - (t250 * t226 + (-t147 + t251) * t229) * t226 - t329;
t324 = -t221 * t156 - (t248 * t226 + (-t157 + t249) * t229) * t226 + t276;
t323 = (t132 + t133) * t207 * t227;
t322 = t208 ^ 2;
t230 = -pkin(8) - pkin(7);
t320 = t226 / 0.2e1;
t319 = -t229 / 0.2e1;
t225 = sin(qJ(2));
t318 = pkin(2) * t225;
t317 = pkin(3) * t211;
t316 = pkin(4) * t208;
t228 = cos(qJ(2));
t210 = t228 * pkin(2) + pkin(1);
t313 = t207 * t325 + t208 * t326 + t323;
t312 = rSges(3,1) * t228;
t310 = rSges(3,2) * t225;
t308 = t229 * rSges(3,3);
t46 = -t208 * t98 + (-t102 * t224 + t106 * t227) * t207;
t307 = t46 * t229;
t47 = -t208 * t99 + (-t103 * t224 + t107 * t227) * t207;
t306 = t47 * t226;
t48 = -t208 * t100 + (-t104 * t224 + t108 * t227) * t207;
t305 = t48 * t229;
t49 = -t208 * t101 + (-t105 * t224 + t109 * t227) * t207;
t304 = t49 * t226;
t301 = Icges(3,4) * t225;
t300 = Icges(3,4) * t228;
t189 = pkin(3) * t212 + t210;
t182 = t229 * t189;
t202 = t229 * t210;
t286 = t229 * (t182 - t202) + (t189 - t210) * t220;
t285 = (t314 - rSges(7,3)) * t208 + (rSges(7,1) * t227 - rSges(7,2) * t224 + t315) * t207;
t241 = rSges(5,1) * t291 - rSges(5,2) * t292 + t226 * rSges(5,3);
t263 = rSges(5,1) * t208 - rSges(5,2) * t207;
t88 = t226 * (-t229 * rSges(5,3) + t263 * t226) + t229 * t241;
t137 = -t208 * rSges(6,3) + (rSges(6,1) * t227 - rSges(6,2) * t224) * t207;
t181 = t207 * pkin(4) - t208 * pkin(10);
t284 = -t137 - t181;
t218 = t229 * pkin(7);
t283 = t226 * (t218 + (-pkin(1) + t210) * t226) + t229 * (-t229 * pkin(1) + t202 - t330);
t242 = t226 * rSges(4,3) + t229 * t264;
t97 = t226 * (-t229 * rSges(4,3) + t264 * t226) + t229 * t242;
t282 = t220 * (pkin(10) * t207 + t316) + t229 * t281;
t280 = t226 * rSges(3,3) + t229 * t312;
t278 = t220 + t221;
t277 = (t220 * t147 + (t251 * t229 + (-t146 + t250) * t226) * t229 + t328) * t226;
t275 = -t181 - t285;
t115 = t175 * rSges(6,1) + t174 * rSges(6,2) + rSges(6,3) * t292;
t188 = t211 * rSges(4,1) + t212 * rSges(4,2);
t272 = -t188 - t318;
t180 = t207 * rSges(5,1) + t208 * rSges(5,2);
t271 = -t180 - t317;
t270 = -t181 - t317;
t269 = t226 * (t220 * t157 + (t249 * t229 + (-t156 + t248) * t226) * t229) + t277;
t219 = -pkin(9) + t230;
t268 = -t226 * t219 + t182;
t262 = -t173 * rSges(6,1) - t172 * rSges(6,2);
t113 = rSges(6,3) * t294 - t262;
t50 = t226 * t113 + t229 * t115 + t282;
t55 = t88 + t286;
t267 = -t137 + t270;
t266 = -t317 - t318;
t265 = -t310 + t312;
t260 = Icges(3,1) * t228 - t301;
t257 = -Icges(3,2) * t225 + t300;
t254 = Icges(3,5) * t228 - Icges(3,6) * t225;
t178 = Icges(5,2) * t208 + t297;
t179 = Icges(5,1) * t207 + t296;
t245 = -t178 * t207 + t179 * t208;
t186 = Icges(4,2) * t212 + t299;
t187 = Icges(4,1) * t211 + t298;
t244 = -t186 * t211 + t187 * t212;
t243 = t270 - t285;
t24 = t226 * t303 + t229 * t327 + t282;
t240 = -t180 + t266;
t239 = -t181 + t266;
t26 = t50 + t286;
t237 = -t137 + t239;
t236 = t276 * t229 + t277;
t22 = t24 + t286;
t235 = t239 - t285;
t234 = -(t306 - t307 + t304 - t305) * t208 / 0.2e1 + t331 * t320 + t332 * t319 + t329 * t294 / 0.2e1 + t328 * t292 / 0.2e1;
t233 = t324 * t229 + t269;
t177 = Icges(5,5) * t207 + Icges(5,6) * t208;
t232 = -t307 / 0.2e1 + t306 / 0.2e1 - t305 / 0.2e1 + t304 / 0.2e1 + (t208 * t149 + t207 * t151 + t226 * t177 + t245 * t229 - t334) * t320 + (t208 * t148 + t207 * t150 - t229 * t177 + t245 * t226 - t335) * t319;
t185 = Icges(4,5) * t211 + Icges(4,6) * t212;
t231 = t232 + (t212 * t159 + t211 * t161 + t226 * t185 + t244 * t229) * t320 + (t212 * t158 + t211 * t160 - t229 * t185 + t244 * t226) * t319;
t201 = t229 * rSges(2,1) - t226 * rSges(2,2);
t200 = -t226 * rSges(2,1) - t229 * rSges(2,2);
t199 = t225 * rSges(3,1) + t228 * rSges(3,2);
t167 = Icges(3,3) * t226 + t254 * t229;
t166 = -Icges(3,3) * t229 + t254 * t226;
t153 = t272 * t229;
t152 = t272 * t226;
t141 = t330 + (pkin(1) - t310) * t229 + t280;
t140 = t308 + t218 + (-pkin(1) - t265) * t226;
t139 = t271 * t229;
t138 = t271 * t226;
t127 = t240 * t229;
t126 = t240 * t226;
t124 = -t226 * t230 + t202 + t242;
t123 = (rSges(4,3) - t230) * t229 + (-t210 - t264) * t226;
t118 = t229 * (-t229 * t310 + t280) + (t265 * t226 - t308) * t226;
t117 = t241 + t268;
t116 = (rSges(5,3) - t219) * t229 + (-t189 - t263) * t226;
t111 = t284 * t229;
t110 = t284 * t226;
t87 = t267 * t229;
t86 = t267 * t226;
t81 = t237 * t229;
t80 = t237 * t226;
t75 = t275 * t229;
t74 = t275 * t226;
t73 = t268 + t115 + t281;
t72 = -t229 * t219 + (-t316 - t189 + (-rSges(6,3) - pkin(10)) * t207) * t226 + t262;
t71 = -t208 * t115 - t137 * t292;
t70 = t208 * t113 + t137 * t294;
t69 = t243 * t229;
t68 = t243 * t226;
t67 = t268 + t333;
t66 = (pkin(5) * t224 - t219) * t229 + (-t340 * t207 - t208 * t209 - t189) * t226 + t261;
t65 = t235 * t229;
t64 = t235 * t226;
t63 = t97 + t283;
t60 = (t113 * t229 - t115 * t226) * t207;
t39 = t55 + t283;
t38 = -t208 * t327 - t285 * t292;
t37 = t303 * t208 + t285 * t294;
t25 = (-t226 * t327 + t303 * t229) * t207;
t23 = t26 + t283;
t19 = t22 + t283;
t1 = [t212 * t186 + t211 * t187 + t228 * (Icges(3,2) * t228 + t301) + t225 * (Icges(3,1) * t225 + t300) + Icges(2,3) + (t178 + t326) * t208 + (t179 + t325) * t207 + m(7) * (t66 ^ 2 + t67 ^ 2) + m(6) * (t72 ^ 2 + t73 ^ 2) + m(5) * (t116 ^ 2 + t117 ^ 2) + m(4) * (t123 ^ 2 + t124 ^ 2) + m(3) * (t140 ^ 2 + t141 ^ 2) + m(2) * (t200 ^ 2 + t201 ^ 2) + t323; t231 + (t221 / 0.2e1 + t220 / 0.2e1) * (Icges(3,5) * t225 + Icges(3,6) * t228) + m(3) * (-t140 * t229 - t141 * t226) * t199 + m(4) * (t153 * t123 + t152 * t124) + m(5) * (t127 * t116 + t126 * t117) + m(6) * (t81 * t72 + t80 * t73) + m(7) * (t64 * t67 + t65 * t66) + (t228 * (Icges(3,6) * t226 + t257 * t229) + t225 * (Icges(3,5) * t226 + t260 * t229)) * t320 + (t228 * (-Icges(3,6) * t229 + t257 * t226) + t225 * (-Icges(3,5) * t229 + t260 * t226)) * t319; m(7) * (t19 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(6) * (t23 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t126 ^ 2 + t127 ^ 2 + t39 ^ 2) + m(4) * (t152 ^ 2 + t153 ^ 2 + t63 ^ 2) + m(3) * (t199 ^ 2 * t278 + t118 ^ 2) + t226 * t220 * t167 + t269 + (-t221 * t166 + (-t226 * t166 + t229 * t167) * t226 + t324) * t229; t231 + m(4) * (-t123 * t229 - t124 * t226) * t188 + m(7) * (t69 * t66 + t68 * t67) + m(6) * (t87 * t72 + t86 * t73) + m(5) * (t139 * t116 + t138 * t117); m(7) * (t22 * t19 + t68 * t64 + t69 * t65) + m(6) * (t26 * t23 + t86 * t80 + t87 * t81) + m(5) * (t138 * t126 + t139 * t127 + t55 * t39) + m(4) * (t97 * t63 + (-t152 * t226 - t153 * t229) * t188) + t233; m(7) * (t22 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t26 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(5) * (t138 ^ 2 + t139 ^ 2 + t55 ^ 2) + m(4) * (t188 ^ 2 * t278 + t97 ^ 2) + t233; t232 + m(7) * (t75 * t66 + t74 * t67) + m(6) * (t110 * t73 + t111 * t72) + m(5) * (-t116 * t229 - t117 * t226) * t180; m(7) * (t24 * t19 + t74 * t64 + t75 * t65) + m(6) * (t110 * t80 + t111 * t81 + t50 * t23) + m(5) * (t88 * t39 + (-t126 * t226 - t127 * t229) * t180) + t236; m(7) * (t24 * t22 + t74 * t68 + t75 * t69) + m(6) * (t110 * t86 + t111 * t87 + t50 * t26) + m(5) * (t88 * t55 + (-t138 * t226 - t139 * t229) * t180) + t236; m(7) * (t24 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(6) * (t110 ^ 2 + t111 ^ 2 + t50 ^ 2) + m(5) * (t180 ^ 2 * t278 + t88 ^ 2) + t236; -t313 * t208 + m(7) * (t37 * t66 + t38 * t67) + m(6) * (t70 * t72 + t71 * t73) + ((t49 / 0.2e1 + t47 / 0.2e1 + t59 / 0.2e1 + t58 / 0.2e1) * t229 + (t48 / 0.2e1 + t46 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1) * t226) * t207; m(7) * (t25 * t19 + t37 * t65 + t38 * t64) + m(6) * (t60 * t23 + t70 * t81 + t71 * t80) + t234; m(7) * (t25 * t22 + t37 * t69 + t38 * t68) + m(6) * (t60 * t26 + t70 * t87 + t71 * t86) + t234; m(7) * (t25 * t24 + t37 * t75 + t38 * t74) + m(6) * (t71 * t110 + t70 * t111 + t60 * t50) + t234; m(7) * (t25 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(6) * (t60 ^ 2 + t70 ^ 2 + t71 ^ 2) + t313 * t322 + (t331 * t229 + t332 * t226 + ((-t47 - t49) * t229 + (-t46 - t48) * t226) * t208) * t207; m(7) * (t226 * t67 + t229 * t66) * t207; m(7) * (-t208 * t19 + (t226 * t64 + t229 * t65) * t207); m(7) * (-t208 * t22 + (t226 * t68 + t229 * t69) * t207); m(7) * (-t208 * t24 + (t226 * t74 + t229 * t75) * t207); m(7) * (-t208 * t25 + (t226 * t38 + t229 * t37) * t207); m(7) * (t207 ^ 2 * t278 + t322);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
