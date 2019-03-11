% Calculate joint inertia matrix for
% S6RRRPRP2
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:08
% EndTime: 2019-03-09 16:35:17
% DurationCPUTime: 4.37s
% Computational Cost: add. (9751->421), mult. (9075->605), div. (0->0), fcn. (9654->10), ass. (0->212)
t214 = qJ(2) + qJ(3);
t203 = pkin(10) + t214;
t200 = cos(t203);
t218 = cos(qJ(5));
t220 = cos(qJ(1));
t276 = t218 * t220;
t215 = sin(qJ(5));
t217 = sin(qJ(1));
t279 = t215 * t217;
t167 = t200 * t279 + t276;
t277 = t217 * t218;
t278 = t215 * t220;
t168 = t200 * t277 - t278;
t199 = sin(t203);
t283 = t199 * t217;
t89 = Icges(7,5) * t168 + Icges(7,6) * t283 + Icges(7,3) * t167;
t95 = Icges(6,4) * t168 - Icges(6,2) * t167 + Icges(6,6) * t283;
t338 = t89 - t95;
t169 = t200 * t278 - t277;
t170 = t200 * t276 + t279;
t281 = t199 * t220;
t90 = Icges(7,5) * t170 + Icges(7,6) * t281 + Icges(7,3) * t169;
t96 = Icges(6,4) * t170 - Icges(6,2) * t169 + Icges(6,6) * t281;
t337 = t90 - t96;
t91 = Icges(6,5) * t168 - Icges(6,6) * t167 + Icges(6,3) * t283;
t93 = Icges(7,4) * t168 + Icges(7,2) * t283 + Icges(7,6) * t167;
t336 = t91 + t93;
t92 = Icges(6,5) * t170 - Icges(6,6) * t169 + Icges(6,3) * t281;
t94 = Icges(7,4) * t170 + Icges(7,2) * t281 + Icges(7,6) * t169;
t335 = t92 + t94;
t97 = Icges(7,1) * t168 + Icges(7,4) * t283 + Icges(7,5) * t167;
t99 = Icges(6,1) * t168 - Icges(6,4) * t167 + Icges(6,5) * t283;
t334 = t97 + t99;
t100 = Icges(6,1) * t170 - Icges(6,4) * t169 + Icges(6,5) * t281;
t98 = Icges(7,1) * t170 + Icges(7,4) * t281 + Icges(7,5) * t169;
t333 = t100 + t98;
t332 = Icges(4,3) + Icges(5,3);
t204 = sin(t214);
t205 = cos(t214);
t331 = Icges(4,5) * t205 + Icges(5,5) * t200 - Icges(4,6) * t204 - Icges(5,6) * t199;
t313 = rSges(7,3) + qJ(6);
t317 = rSges(7,1) + pkin(5);
t330 = -t313 * t167 - t317 * t168;
t329 = t338 * t167 + t334 * t168 + t336 * t283;
t328 = t337 * t167 + t333 * t168 + t335 * t283;
t327 = t338 * t169 + t334 * t170 + t336 * t281;
t326 = t337 * t169 + t333 * t170 + t335 * t281;
t121 = -Icges(7,6) * t200 + (Icges(7,5) * t218 + Icges(7,3) * t215) * t199;
t123 = -Icges(7,2) * t200 + (Icges(7,4) * t218 + Icges(7,6) * t215) * t199;
t125 = -Icges(7,4) * t200 + (Icges(7,1) * t218 + Icges(7,5) * t215) * t199;
t54 = t121 * t167 + t123 * t283 + t125 * t168;
t122 = -Icges(6,3) * t200 + (Icges(6,5) * t218 - Icges(6,6) * t215) * t199;
t124 = -Icges(6,6) * t200 + (Icges(6,4) * t218 - Icges(6,2) * t215) * t199;
t126 = -Icges(6,5) * t200 + (Icges(6,1) * t218 - Icges(6,4) * t215) * t199;
t55 = t122 * t283 - t124 * t167 + t126 * t168;
t325 = -t54 - t55;
t56 = t121 * t169 + t123 * t281 + t125 * t170;
t57 = t122 * t281 - t124 * t169 + t126 * t170;
t324 = -t56 - t57;
t323 = t217 * t332 + t331 * t220;
t322 = -t331 * t217 + t220 * t332;
t285 = Icges(5,4) * t200;
t244 = -Icges(5,2) * t199 + t285;
t140 = Icges(5,6) * t217 + t220 * t244;
t286 = Icges(5,4) * t199;
t247 = Icges(5,1) * t200 - t286;
t142 = Icges(5,5) * t217 + t220 * t247;
t287 = Icges(4,4) * t205;
t245 = -Icges(4,2) * t204 + t287;
t152 = Icges(4,6) * t217 + t220 * t245;
t288 = Icges(4,4) * t204;
t248 = Icges(4,1) * t205 - t288;
t154 = Icges(4,5) * t217 + t220 * t248;
t321 = t140 * t199 - t142 * t200 + t152 * t204 - t154 * t205;
t139 = -Icges(5,6) * t220 + t217 * t244;
t141 = -Icges(5,5) * t220 + t217 * t247;
t151 = -Icges(4,6) * t220 + t217 * t245;
t153 = -Icges(4,5) * t220 + t217 * t248;
t320 = t139 * t199 - t141 * t200 + t151 * t204 - t153 * t205;
t212 = t217 ^ 2;
t319 = t325 * t200 + (t217 * t329 + t220 * t328) * t199;
t318 = t324 * t200 + (t217 * t327 + t220 * t326) * t199;
t316 = t217 * pkin(7);
t315 = t217 * t328 - t220 * t329;
t314 = t217 * t326 - t220 * t327;
t274 = rSges(7,2) * t283 - t330;
t312 = rSges(7,2) * t281 + t313 * t169 + t170 * t317;
t311 = -t122 - t123;
t310 = Icges(4,5) * t204 + Icges(5,5) * t199 + Icges(4,6) * t205 + Icges(5,6) * t200;
t172 = Icges(5,2) * t200 + t286;
t173 = Icges(5,1) * t199 + t285;
t179 = Icges(4,2) * t205 + t288;
t180 = Icges(4,1) * t204 + t287;
t309 = -t172 * t199 + t173 * t200 - t179 * t204 + t180 * t205;
t252 = rSges(4,1) * t205 - rSges(4,2) * t204;
t284 = t199 * t215;
t308 = t121 * t284 + (t125 + t126) * t199 * t218;
t213 = t220 ^ 2;
t307 = -t315 + t322 * t213 + (t321 * t217 + (-t320 + t323) * t220) * t217;
t221 = -pkin(8) - pkin(7);
t305 = t217 / 0.2e1;
t304 = -t220 / 0.2e1;
t216 = sin(qJ(2));
t303 = pkin(2) * t216;
t302 = pkin(3) * t204;
t301 = pkin(4) * t200;
t219 = cos(qJ(2));
t202 = t219 * pkin(2) + pkin(1);
t300 = -t124 * t284 + t200 * t311 + t308;
t299 = rSges(3,1) * t219;
t297 = rSges(3,2) * t216;
t295 = t220 * rSges(3,3);
t43 = -t200 * t93 + (t215 * t89 + t218 * t97) * t199;
t294 = t43 * t220;
t44 = -t200 * t94 + (t215 * t90 + t218 * t98) * t199;
t293 = t44 * t217;
t45 = -t200 * t91 + (-t215 * t95 + t218 * t99) * t199;
t292 = t45 * t220;
t46 = -t200 * t92 + (t100 * t218 - t215 * t96) * t199;
t291 = t46 * t217;
t290 = Icges(3,4) * t216;
t289 = Icges(3,4) * t219;
t280 = t200 * t220;
t211 = -qJ(4) + t221;
t275 = t220 * t211;
t182 = pkin(3) * t205 + t202;
t176 = t220 * t182;
t195 = t220 * t202;
t272 = t220 * (t176 - t195) + (t182 - t202) * t212;
t270 = -t200 * rSges(7,2) + (t215 * t313 + t218 * t317) * t199;
t210 = t220 * pkin(7);
t269 = t217 * (t210 + (-pkin(1) + t202) * t217) + t220 * (-t220 * pkin(1) + t195 - t316);
t231 = t217 * rSges(4,3) + t220 * t252;
t88 = t217 * (-t220 * rSges(4,3) + t217 * t252) + t220 * t231;
t268 = pkin(4) * t280 + pkin(9) * t281;
t267 = t217 * rSges(3,3) + t220 * t299;
t265 = t212 + t213;
t104 = t170 * rSges(6,1) - t169 * rSges(6,2) + rSges(6,3) * t281;
t181 = rSges(4,1) * t204 + rSges(4,2) * t205;
t262 = -t181 - t303;
t174 = rSges(5,1) * t199 + rSges(5,2) * t200;
t261 = -t174 - t302;
t175 = pkin(4) * t199 - pkin(9) * t200;
t260 = -t175 - t302;
t259 = -t182 - t301;
t258 = (t323 * t212 + ((-t321 + t322) * t217 + t320 * t220) * t220 + t314) * t217;
t257 = -t217 * t211 + t176;
t230 = rSges(5,1) * t280 - rSges(5,2) * t281 + t217 * rSges(5,3);
t251 = rSges(5,1) * t200 - rSges(5,2) * t199;
t53 = t217 * (-t220 * rSges(5,3) + t217 * t251) + t220 * t230 + t272;
t256 = t212 * (pkin(9) * t199 + t301) + t220 * t268 + t272;
t128 = -t200 * rSges(6,3) + (rSges(6,1) * t218 - rSges(6,2) * t215) * t199;
t255 = -t128 + t260;
t254 = -t302 - t303;
t253 = -t297 + t299;
t250 = -rSges(6,1) * t168 + rSges(6,2) * t167;
t249 = Icges(3,1) * t219 - t290;
t246 = -Icges(3,2) * t216 + t289;
t243 = Icges(3,5) * t219 - Icges(3,6) * t216;
t232 = t260 - t270;
t229 = -t174 + t254;
t228 = -t175 + t254;
t227 = t257 + t268;
t102 = rSges(6,3) * t283 - t250;
t24 = t217 * t102 + t220 * t104 + t256;
t226 = -t128 + t228;
t225 = t228 - t270;
t22 = t217 * t274 + t220 * t312 + t256;
t224 = -(t293 - t294 + t291 - t292) * t200 / 0.2e1 + t318 * t305 + t319 * t304 + t315 * t283 / 0.2e1 + t314 * t281 / 0.2e1;
t223 = t307 * t220 + t258;
t222 = -t294 / 0.2e1 + t293 / 0.2e1 - t292 / 0.2e1 + t291 / 0.2e1 + (t140 * t200 + t142 * t199 + t152 * t205 + t154 * t204 + t217 * t310 + t220 * t309 - t324) * t305 + (t139 * t200 + t141 * t199 + t151 * t205 + t153 * t204 + t217 * t309 - t220 * t310 - t325) * t304;
t194 = rSges(2,1) * t220 - rSges(2,2) * t217;
t193 = -rSges(2,1) * t217 - rSges(2,2) * t220;
t192 = rSges(3,1) * t216 + rSges(3,2) * t219;
t162 = Icges(3,3) * t217 + t220 * t243;
t161 = -Icges(3,3) * t220 + t217 * t243;
t148 = t262 * t220;
t147 = t262 * t217;
t134 = t316 + (pkin(1) - t297) * t220 + t267;
t133 = t295 + t210 + (-pkin(1) - t253) * t217;
t132 = t261 * t220;
t131 = t261 * t217;
t120 = t229 * t220;
t119 = t229 * t217;
t118 = -t217 * t221 + t195 + t231;
t117 = (rSges(4,3) - t221) * t220 + (-t202 - t252) * t217;
t107 = t220 * (-t220 * t297 + t267) + (t217 * t253 - t295) * t217;
t106 = t230 + t257;
t105 = (rSges(5,3) - t211) * t220 + (-t182 - t251) * t217;
t83 = t255 * t220;
t82 = t255 * t217;
t77 = t226 * t220;
t76 = t226 * t217;
t71 = t232 * t220;
t70 = t232 * t217;
t69 = t225 * t220;
t68 = t225 * t217;
t67 = t227 + t104;
t66 = -t275 + ((-rSges(6,3) - pkin(9)) * t199 + t259) * t217 + t250;
t65 = -t104 * t200 - t128 * t281;
t64 = t102 * t200 + t128 * t283;
t63 = t88 + t269;
t60 = (t102 * t220 - t104 * t217) * t199;
t59 = t227 + t312;
t58 = -t275 + ((-rSges(7,2) - pkin(9)) * t199 + t259) * t217 + t330;
t48 = -t200 * t312 - t270 * t281;
t47 = t200 * t274 + t270 * t283;
t36 = t53 + t269;
t25 = (-t217 * t312 + t220 * t274) * t199;
t23 = t24 + t269;
t19 = t22 + t269;
t1 = [t205 * t179 + t204 * t180 + t219 * (Icges(3,2) * t219 + t290) + t216 * (Icges(3,1) * t216 + t289) + Icges(2,3) + (-t124 * t215 + t173) * t199 + (t172 + t311) * t200 + m(7) * (t58 ^ 2 + t59 ^ 2) + m(6) * (t66 ^ 2 + t67 ^ 2) + m(5) * (t105 ^ 2 + t106 ^ 2) + m(4) * (t117 ^ 2 + t118 ^ 2) + m(3) * (t133 ^ 2 + t134 ^ 2) + m(2) * (t193 ^ 2 + t194 ^ 2) + t308; m(7) * (t58 * t69 + t59 * t68) + m(6) * (t66 * t77 + t67 * t76) + m(5) * (t105 * t120 + t106 * t119) + m(4) * (t117 * t148 + t118 * t147) + t222 + m(3) * (-t133 * t220 - t134 * t217) * t192 + (t212 / 0.2e1 + t213 / 0.2e1) * (Icges(3,5) * t216 + Icges(3,6) * t219) + (t219 * (Icges(3,6) * t217 + t220 * t246) + t216 * (Icges(3,5) * t217 + t220 * t249)) * t305 + (t219 * (-Icges(3,6) * t220 + t217 * t246) + t216 * (-Icges(3,5) * t220 + t217 * t249)) * t304; m(7) * (t19 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t23 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t119 ^ 2 + t120 ^ 2 + t36 ^ 2) + m(4) * (t147 ^ 2 + t148 ^ 2 + t63 ^ 2) + t217 * t212 * t162 + m(3) * (t192 ^ 2 * t265 + t107 ^ 2) + t258 + (-t213 * t161 + (-t217 * t161 + t220 * t162) * t217 + t307) * t220; m(7) * (t58 * t71 + t59 * t70) + m(6) * (t66 * t83 + t67 * t82) + m(5) * (t105 * t132 + t106 * t131) + t222 + m(4) * (-t117 * t220 - t118 * t217) * t181; m(7) * (t19 * t22 + t68 * t70 + t69 * t71) + m(6) * (t23 * t24 + t76 * t82 + t77 * t83) + m(5) * (t119 * t131 + t120 * t132 + t36 * t53) + m(4) * (t88 * t63 + (-t147 * t217 - t148 * t220) * t181) + t223; m(7) * (t22 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t24 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(5) * (t131 ^ 2 + t132 ^ 2 + t53 ^ 2) + m(4) * (t181 ^ 2 * t265 + t88 ^ 2) + t223; m(7) * (t217 * t58 - t220 * t59) + m(6) * (t217 * t66 - t220 * t67) + m(5) * (t105 * t217 - t106 * t220); m(7) * (t217 * t69 - t220 * t68) + m(6) * (t217 * t77 - t220 * t76) + m(5) * (-t119 * t220 + t120 * t217); m(7) * (t217 * t71 - t220 * t70) + m(6) * (t217 * t83 - t220 * t82) + m(5) * (-t131 * t220 + t132 * t217); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t265; -t300 * t200 + m(7) * (t47 * t58 + t48 * t59) + m(6) * (t64 * t66 + t65 * t67) + ((t46 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1 + t44 / 0.2e1) * t220 + (t54 / 0.2e1 + t43 / 0.2e1 + t45 / 0.2e1 + t55 / 0.2e1) * t217) * t199; m(7) * (t19 * t25 + t47 * t69 + t48 * t68) + m(6) * (t23 * t60 + t64 * t77 + t65 * t76) + t224; m(7) * (t22 * t25 + t47 * t71 + t48 * t70) + m(6) * (t24 * t60 + t64 * t83 + t65 * t82) + t224; m(6) * (t217 * t64 - t220 * t65) + m(7) * (t217 * t47 - t220 * t48); m(7) * (t25 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t60 ^ 2 + t64 ^ 2 + t65 ^ 2) + t300 * t200 ^ 2 + (t318 * t220 + t319 * t217 + ((-t44 - t46) * t220 + (-t43 - t45) * t217) * t200) * t199; m(7) * (t167 * t59 + t169 * t58); m(7) * (t167 * t68 + t169 * t69 + t19 * t284); m(7) * (t167 * t70 + t169 * t71 + t22 * t284); m(7) * (-t167 * t220 + t169 * t217); m(7) * (t167 * t48 + t169 * t47 + t25 * t284); m(7) * (t199 ^ 2 * t215 ^ 2 + t167 ^ 2 + t169 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
