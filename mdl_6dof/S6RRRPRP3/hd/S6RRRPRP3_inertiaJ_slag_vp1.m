% Calculate joint inertia matrix for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:09
% EndTime: 2019-03-09 16:39:22
% DurationCPUTime: 5.22s
% Computational Cost: add. (10670->466), mult. (10442->670), div. (0->0), fcn. (11208->10), ass. (0->224)
t219 = pkin(10) + qJ(5);
t213 = cos(t219);
t222 = qJ(2) + qJ(3);
t215 = cos(t222);
t229 = cos(qJ(1));
t212 = sin(t219);
t227 = sin(qJ(1));
t290 = t212 * t227;
t167 = t213 * t229 + t215 * t290;
t281 = t227 * t213;
t168 = -t212 * t229 + t215 * t281;
t214 = sin(t222);
t288 = t214 * t227;
t90 = Icges(7,5) * t168 + Icges(7,6) * t288 + Icges(7,3) * t167;
t96 = Icges(6,4) * t168 - Icges(6,2) * t167 + Icges(6,6) * t288;
t337 = t90 - t96;
t286 = t215 * t229;
t169 = t212 * t286 - t281;
t170 = t213 * t286 + t290;
t287 = t214 * t229;
t91 = Icges(7,5) * t170 + Icges(7,6) * t287 + Icges(7,3) * t169;
t97 = Icges(6,4) * t170 - Icges(6,2) * t169 + Icges(6,6) * t287;
t336 = t91 - t97;
t92 = Icges(6,5) * t168 - Icges(6,6) * t167 + Icges(6,3) * t288;
t94 = Icges(7,4) * t168 + Icges(7,2) * t288 + Icges(7,6) * t167;
t335 = t92 + t94;
t93 = Icges(6,5) * t170 - Icges(6,6) * t169 + Icges(6,3) * t287;
t95 = Icges(7,4) * t170 + Icges(7,2) * t287 + Icges(7,6) * t169;
t334 = t93 + t95;
t100 = Icges(6,1) * t168 - Icges(6,4) * t167 + Icges(6,5) * t288;
t98 = Icges(7,1) * t168 + Icges(7,4) * t288 + Icges(7,5) * t167;
t333 = t100 + t98;
t101 = Icges(6,1) * t170 - Icges(6,4) * t169 + Icges(6,5) * t287;
t99 = Icges(7,1) * t170 + Icges(7,4) * t287 + Icges(7,5) * t169;
t332 = t101 + t99;
t319 = rSges(7,3) + qJ(6);
t322 = rSges(7,1) + pkin(5);
t331 = -t319 * t167 - t322 * t168;
t330 = t337 * t167 + t333 * t168 + t335 * t288;
t329 = t336 * t167 + t332 * t168 + t334 * t288;
t328 = t337 * t169 + t333 * t170 + t335 * t287;
t327 = t336 * t169 + t332 * t170 + t334 * t287;
t132 = -Icges(7,6) * t215 + (Icges(7,5) * t213 + Icges(7,3) * t212) * t214;
t134 = -Icges(7,2) * t215 + (Icges(7,4) * t213 + Icges(7,6) * t212) * t214;
t136 = -Icges(7,4) * t215 + (Icges(7,1) * t213 + Icges(7,5) * t212) * t214;
t58 = t132 * t167 + t134 * t288 + t136 * t168;
t133 = -Icges(6,3) * t215 + (Icges(6,5) * t213 - Icges(6,6) * t212) * t214;
t135 = -Icges(6,6) * t215 + (Icges(6,4) * t213 - Icges(6,2) * t212) * t214;
t137 = -Icges(6,5) * t215 + (Icges(6,1) * t213 - Icges(6,4) * t212) * t214;
t59 = t133 * t288 - t135 * t167 + t137 * t168;
t326 = -t58 - t59;
t60 = t132 * t169 + t134 * t287 + t136 * t170;
t61 = t133 * t287 - t135 * t169 + t137 * t170;
t325 = -t60 - t61;
t324 = t326 * t215 + (t330 * t227 + t329 * t229) * t214;
t323 = t325 * t215 + (t328 * t227 + t327 * t229) * t214;
t321 = t329 * t227 - t330 * t229;
t320 = t327 * t227 - t328 * t229;
t279 = rSges(7,2) * t288 - t331;
t318 = rSges(7,2) * t287 + t319 * t169 + t322 * t170;
t317 = -t133 - t134;
t291 = t212 * t214;
t316 = t132 * t291 + (t136 + t137) * t213 * t214;
t224 = cos(pkin(10));
t282 = t224 * t229;
t223 = sin(pkin(10));
t285 = t223 * t227;
t179 = -t215 * t285 - t282;
t283 = t224 * t227;
t284 = t223 * t229;
t180 = t215 * t283 - t284;
t115 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t288;
t181 = -t215 * t284 + t283;
t182 = t215 * t282 + t285;
t116 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t287;
t117 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t288;
t118 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t287;
t119 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t288;
t120 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t287;
t244 = Icges(4,5) * t215 - Icges(4,6) * t214;
t152 = -Icges(4,3) * t229 + t244 * t227;
t153 = Icges(4,3) * t227 + t244 * t229;
t221 = t229 ^ 2;
t293 = Icges(4,4) * t215;
t246 = -Icges(4,2) * t214 + t293;
t155 = Icges(4,6) * t227 + t246 * t229;
t294 = Icges(4,4) * t214;
t248 = Icges(4,1) * t215 - t294;
t157 = Icges(4,5) * t227 + t248 * t229;
t242 = -t155 * t214 + t157 * t215;
t154 = -Icges(4,6) * t229 + t246 * t227;
t156 = -Icges(4,5) * t229 + t248 * t227;
t243 = t154 * t214 - t156 * t215;
t315 = (t115 * t288 + t117 * t179 + t119 * t180) * t229 - t221 * t152 - t321 + (-t116 * t288 - t118 * t179 - t120 * t180 - (-t153 + t243) * t229 - t242 * t227) * t227;
t314 = t215 ^ 2;
t220 = t227 ^ 2;
t313 = m(5) / 0.2e1;
t312 = m(6) / 0.2e1;
t311 = m(7) / 0.2e1;
t309 = t227 / 0.2e1;
t308 = -t229 / 0.2e1;
t226 = sin(qJ(2));
t307 = pkin(2) * t226;
t306 = pkin(3) * t215;
t209 = pkin(4) * t224 + pkin(3);
t305 = -pkin(3) + t209;
t304 = -t135 * t291 + t215 * t317 + t316;
t228 = cos(qJ(2));
t303 = rSges(3,1) * t228;
t302 = rSges(3,2) * t226;
t301 = t229 * rSges(3,3);
t41 = -t215 * t94 + (t212 * t90 + t213 * t98) * t214;
t300 = t41 * t229;
t42 = -t215 * t95 + (t212 * t91 + t213 * t99) * t214;
t299 = t42 * t227;
t43 = -t215 * t92 + (t100 * t213 - t212 * t96) * t214;
t298 = t43 * t229;
t44 = -t215 * t93 + (t101 * t213 - t212 * t97) * t214;
t297 = t44 * t227;
t296 = Icges(3,4) * t226;
t295 = Icges(3,4) * t228;
t292 = qJ(4) * t214;
t230 = -pkin(8) - pkin(7);
t280 = t229 * t230;
t225 = -pkin(9) - qJ(4);
t131 = (qJ(4) + t225) * t215 + t305 * t214;
t188 = pkin(3) * t214 - qJ(4) * t215;
t276 = -t131 - t188;
t275 = -rSges(7,2) * t215 + (t212 * t319 + t213 * t322) * t214;
t210 = pkin(2) * t228 + pkin(1);
t199 = t229 * t210;
t218 = t229 * pkin(7);
t274 = t227 * (t280 + t218 + (-pkin(1) + t210) * t227) + t229 * (-t229 * pkin(1) + t199 + (-pkin(7) - t230) * t227);
t237 = rSges(4,1) * t286 - rSges(4,2) * t287 + t227 * rSges(4,3);
t252 = rSges(4,1) * t215 - rSges(4,2) * t214;
t108 = t227 * (-t229 * rSges(4,3) + t252 * t227) + t229 * t237;
t149 = -rSges(5,3) * t215 + (rSges(5,1) * t224 - rSges(5,2) * t223) * t214;
t273 = -t149 - t188;
t270 = pkin(3) * t286 + qJ(4) * t287;
t272 = t220 * (t292 + t306) + t229 * t270;
t271 = -pkin(4) * t284 - t225 * t288;
t269 = t227 * rSges(3,3) + t229 * t303;
t268 = t220 + t221;
t139 = -rSges(6,3) * t215 + (rSges(6,1) * t213 - rSges(6,2) * t212) * t214;
t267 = -t139 + t276;
t105 = t170 * rSges(6,1) - t169 * rSges(6,2) + rSges(6,3) * t287;
t266 = t182 * rSges(5,1) + t181 * rSges(5,2) + rSges(5,3) * t287;
t263 = -t188 - t307;
t189 = rSges(4,1) * t214 + rSges(4,2) * t215;
t262 = -t189 - t307;
t261 = (t220 * t153 + (t116 * t287 + t118 * t181 + t120 * t182) * t227 + (-t115 * t287 - t117 * t181 - t119 * t182 + (-t152 + t242) * t227 + t243 * t229) * t229 + t320) * t227;
t260 = -t227 * t230 + t199;
t259 = -t209 * t215 - t210;
t236 = pkin(4) * t285 + t209 * t286 - t225 * t287;
t258 = t227 * ((t305 * t215 - t292) * t227 + t271) + t229 * (t236 - t270) + t272;
t257 = -t275 + t276;
t251 = -t180 * rSges(5,1) - t179 * rSges(5,2);
t53 = t227 * (rSges(5,3) * t288 - t251) + t229 * t266 + t272;
t256 = -t131 + t263;
t255 = -t149 + t263;
t254 = -t271 - t280;
t253 = -t302 + t303;
t250 = -t168 * rSges(6,1) + t167 * rSges(6,2);
t249 = Icges(3,1) * t228 - t296;
t247 = -Icges(3,2) * t226 + t295;
t245 = Icges(3,5) * t228 - Icges(3,6) * t226;
t186 = Icges(4,2) * t215 + t294;
t187 = Icges(4,1) * t214 + t293;
t239 = -t186 * t214 + t187 * t215;
t238 = -t139 + t256;
t103 = rSges(6,3) * t288 - t250;
t26 = t227 * t103 + t229 * t105 + t258;
t235 = t256 - t275;
t24 = t227 * t279 + t229 * t318 + t258;
t234 = -(t299 - t300 + t297 - t298) * t215 / 0.2e1 + t323 * t309 + t324 * t308 + t321 * t288 / 0.2e1 + t320 * t287 / 0.2e1;
t233 = t315 * t229 + t261;
t232 = t236 + t260;
t146 = -Icges(5,3) * t215 + (Icges(5,5) * t224 - Icges(5,6) * t223) * t214;
t147 = -Icges(5,6) * t215 + (Icges(5,4) * t224 - Icges(5,2) * t223) * t214;
t148 = -Icges(5,5) * t215 + (Icges(5,1) * t224 - Icges(5,4) * t223) * t214;
t185 = Icges(4,5) * t214 + Icges(4,6) * t215;
t231 = -t298 / 0.2e1 - t300 / 0.2e1 + t297 / 0.2e1 + t299 / 0.2e1 + (t146 * t287 + t147 * t181 + t148 * t182 + t185 * t227 + t239 * t229 + (-t116 + t155) * t215 + (-t118 * t223 + t120 * t224 + t157) * t214 - t325) * t309 + (t146 * t288 + t147 * t179 + t148 * t180 - t185 * t229 + t239 * t227 + (t154 - t115) * t215 + (-t117 * t223 + t119 * t224 + t156) * t214 - t326) * t308;
t211 = t214 ^ 2;
t197 = rSges(2,1) * t229 - rSges(2,2) * t227;
t196 = -rSges(2,1) * t227 - rSges(2,2) * t229;
t195 = rSges(3,1) * t226 + rSges(3,2) * t228;
t172 = Icges(3,3) * t227 + t245 * t229;
t171 = -Icges(3,3) * t229 + t245 * t227;
t151 = t262 * t229;
t150 = t262 * t227;
t141 = t227 * pkin(7) + (pkin(1) - t302) * t229 + t269;
t140 = t301 + t218 + (-pkin(1) - t253) * t227;
t130 = t237 + t260;
t129 = (rSges(4,3) - t230) * t229 + (-t210 - t252) * t227;
t125 = t229 * (-t229 * t302 + t269) + (t253 * t227 - t301) * t227;
t124 = t273 * t229;
t123 = t273 * t227;
t110 = t255 * t229;
t109 = t255 * t227;
t81 = t260 + t266 + t270;
t80 = -t280 + (-t306 - t210 + (-rSges(5,3) - qJ(4)) * t214) * t227 + t251;
t79 = t267 * t229;
t78 = t267 * t227;
t77 = t238 * t229;
t76 = t238 * t227;
t75 = -t105 * t215 - t139 * t287;
t74 = t103 * t215 + t139 * t288;
t73 = t232 + t105;
t72 = (-rSges(6,3) * t214 + t259) * t227 + t250 + t254;
t71 = t108 + t274;
t70 = t257 * t229;
t69 = t257 * t227;
t66 = t235 * t229;
t65 = t235 * t227;
t62 = (t103 * t229 - t105 * t227) * t214;
t52 = t232 + t318;
t51 = (-rSges(7,2) * t214 + t259) * t227 + t254 + t331;
t48 = -t215 * t318 - t275 * t287;
t47 = t279 * t215 + t275 * t288;
t36 = t53 + t274;
t27 = (-t227 * t318 + t279 * t229) * t214;
t25 = t26 + t274;
t19 = t24 + t274;
t1 = [t228 * (Icges(3,2) * t228 + t296) + t226 * (Icges(3,1) * t226 + t295) + Icges(2,3) + (t186 - t146 + t317) * t215 + (-t212 * t135 - t223 * t147 + t224 * t148 + t187) * t214 + m(7) * (t51 ^ 2 + t52 ^ 2) + m(6) * (t72 ^ 2 + t73 ^ 2) + m(5) * (t80 ^ 2 + t81 ^ 2) + m(4) * (t129 ^ 2 + t130 ^ 2) + m(3) * (t140 ^ 2 + t141 ^ 2) + m(2) * (t196 ^ 2 + t197 ^ 2) + t316; m(7) * (t51 * t66 + t52 * t65) + m(6) * (t72 * t77 + t73 * t76) + m(5) * (t109 * t81 + t110 * t80) + m(4) * (t129 * t151 + t130 * t150) + (t220 / 0.2e1 + t221 / 0.2e1) * (Icges(3,5) * t226 + Icges(3,6) * t228) + t231 + m(3) * (-t140 * t229 - t141 * t227) * t195 + ((-Icges(3,6) * t229 + t247 * t227) * t228 + (-Icges(3,5) * t229 + t249 * t227) * t226) * t308 + ((Icges(3,6) * t227 + t247 * t229) * t228 + (Icges(3,5) * t227 + t249 * t229) * t226) * t309; m(7) * (t19 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t25 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2 + t36 ^ 2) + m(4) * (t150 ^ 2 + t151 ^ 2 + t71 ^ 2) + t227 * t220 * t172 + m(3) * (t268 * t195 ^ 2 + t125 ^ 2) + t261 + (-t221 * t171 + (-t227 * t171 + t229 * t172) * t227 + t315) * t229; m(7) * (t51 * t70 + t52 * t69) + m(6) * (t72 * t79 + t73 * t78) + m(5) * (t123 * t81 + t124 * t80) + t231 + m(4) * (-t129 * t229 - t130 * t227) * t189; m(7) * (t24 * t19 + t65 * t69 + t66 * t70) + m(6) * (t26 * t25 + t76 * t78 + t77 * t79) + m(5) * (t109 * t123 + t110 * t124 + t53 * t36) + m(4) * (t108 * t71 + (-t150 * t227 - t151 * t229) * t189) + t233; m(7) * (t24 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t26 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t123 ^ 2 + t124 ^ 2 + t53 ^ 2) + m(4) * (t268 * t189 ^ 2 + t108 ^ 2) + t233; 0.2e1 * ((t227 * t52 + t229 * t51) * t311 + (t227 * t73 + t229 * t72) * t312 + (t227 * t81 + t229 * t80) * t313) * t214; m(7) * (-t19 * t215 + (t227 * t65 + t229 * t66) * t214) + m(6) * (-t215 * t25 + (t227 * t76 + t229 * t77) * t214) + m(5) * (-t215 * t36 + (t109 * t227 + t110 * t229) * t214); m(7) * (-t215 * t24 + (t227 * t69 + t229 * t70) * t214) + m(6) * (-t215 * t26 + (t227 * t78 + t229 * t79) * t214) + m(5) * (-t215 * t53 + (t123 * t227 + t124 * t229) * t214); 0.2e1 * (t313 + t312 + t311) * (t268 * t211 + t314); -t304 * t215 + m(7) * (t47 * t51 + t48 * t52) + m(6) * (t72 * t74 + t73 * t75) + ((t61 / 0.2e1 + t60 / 0.2e1 + t44 / 0.2e1 + t42 / 0.2e1) * t229 + (t59 / 0.2e1 + t58 / 0.2e1 + t43 / 0.2e1 + t41 / 0.2e1) * t227) * t214; m(7) * (t19 * t27 + t47 * t66 + t48 * t65) + m(6) * (t62 * t25 + t74 * t77 + t75 * t76) + t234; m(7) * (t27 * t24 + t47 * t70 + t48 * t69) + m(6) * (t62 * t26 + t74 * t79 + t75 * t78) + t234; m(6) * (-t215 * t62 + (t227 * t75 + t229 * t74) * t214) + m(7) * (-t215 * t27 + (t227 * t48 + t229 * t47) * t214); m(7) * (t27 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t62 ^ 2 + t74 ^ 2 + t75 ^ 2) + t304 * t314 + (t323 * t229 + t324 * t227 + ((-t42 - t44) * t229 + (-t41 - t43) * t227) * t215) * t214; m(7) * (t167 * t52 + t169 * t51); m(7) * (t167 * t65 + t169 * t66 + t19 * t291); m(7) * (t167 * t69 + t169 * t70 + t24 * t291); m(7) * (t167 * t227 + t169 * t229 - t212 * t215) * t214; m(7) * (t167 * t48 + t169 * t47 + t27 * t291); m(7) * (t211 * t212 ^ 2 + t167 ^ 2 + t169 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
