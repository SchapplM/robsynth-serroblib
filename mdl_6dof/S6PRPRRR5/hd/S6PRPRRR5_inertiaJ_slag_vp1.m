% Calculate joint inertia matrix for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:57
% EndTime: 2019-03-08 20:41:06
% DurationCPUTime: 5.33s
% Computational Cost: add. (19600->487), mult. (34336->712), div. (0->0), fcn. (43090->12), ass. (0->233)
t310 = Icges(3,1) + Icges(4,2);
t308 = Icges(3,4) + Icges(4,6);
t307 = Icges(3,5) - Icges(4,4);
t309 = Icges(3,2) + Icges(4,3);
t306 = Icges(3,6) - Icges(4,5);
t305 = Icges(3,3) + Icges(4,1);
t227 = sin(pkin(11));
t229 = cos(pkin(11));
t236 = cos(qJ(2));
t230 = cos(pkin(6));
t233 = sin(qJ(2));
t268 = t230 * t233;
t216 = t227 * t236 + t229 * t268;
t218 = -t227 * t268 + t229 * t236;
t228 = sin(pkin(6));
t271 = t228 * t233;
t267 = t230 * t236;
t217 = t227 * t267 + t229 * t233;
t265 = qJ(4) + qJ(5);
t226 = sin(t265);
t248 = cos(t265);
t243 = t228 * t248;
t190 = t217 * t226 + t227 * t243;
t231 = sin(qJ(6));
t234 = cos(qJ(6));
t153 = -t190 * t231 + t218 * t234;
t154 = t190 * t234 + t218 * t231;
t274 = t227 * t228;
t189 = -t217 * t248 + t226 * t274;
t101 = Icges(7,5) * t154 + Icges(7,6) * t153 + Icges(7,3) * t189;
t103 = Icges(7,4) * t154 + Icges(7,2) * t153 + Icges(7,6) * t189;
t105 = Icges(7,1) * t154 + Icges(7,4) * t153 + Icges(7,5) * t189;
t50 = t101 * t189 + t103 * t153 + t105 * t154;
t215 = t227 * t233 - t229 * t267;
t192 = t215 * t226 - t229 * t243;
t155 = -t192 * t231 + t216 * t234;
t156 = t192 * t234 + t216 * t231;
t273 = t228 * t229;
t191 = t215 * t248 + t226 * t273;
t102 = Icges(7,5) * t156 + Icges(7,6) * t155 - Icges(7,3) * t191;
t104 = Icges(7,4) * t156 + Icges(7,2) * t155 - Icges(7,6) * t191;
t106 = Icges(7,1) * t156 + Icges(7,4) * t155 - Icges(7,5) * t191;
t51 = t102 * t189 + t104 * t153 + t106 * t154;
t269 = t228 * t236;
t214 = -t226 * t269 + t230 * t248;
t193 = -t214 * t231 + t234 * t271;
t194 = t214 * t234 + t231 * t271;
t213 = t226 * t230 + t236 * t243;
t118 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t213;
t119 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t213;
t120 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t213;
t60 = t118 * t189 + t119 * t153 + t120 * t154;
t11 = t216 * t51 + t218 * t50 + t271 * t60;
t122 = Icges(6,5) * t190 - Icges(6,6) * t189 + Icges(6,3) * t218;
t124 = Icges(6,4) * t190 - Icges(6,2) * t189 + Icges(6,6) * t218;
t126 = Icges(6,1) * t190 - Icges(6,4) * t189 + Icges(6,5) * t218;
t69 = t122 * t218 - t124 * t189 + t126 * t190;
t123 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t216;
t125 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t216;
t127 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t216;
t70 = t123 * t218 - t125 * t189 + t127 * t190;
t157 = Icges(6,5) * t214 - Icges(6,6) * t213 + Icges(6,3) * t271;
t158 = Icges(6,4) * t214 - Icges(6,2) * t213 + Icges(6,6) * t271;
t159 = Icges(6,1) * t214 - Icges(6,4) * t213 + Icges(6,5) * t271;
t88 = t157 * t218 - t158 * t189 + t159 * t190;
t304 = t216 * t70 + t218 * t69 + t271 * t88 + t11;
t52 = -t101 * t191 + t103 * t155 + t105 * t156;
t53 = -t102 * t191 + t104 * t155 + t106 * t156;
t61 = -t118 * t191 + t119 * t155 + t120 * t156;
t12 = t216 * t53 + t218 * t52 + t271 * t61;
t71 = t122 * t216 + t124 * t191 + t126 * t192;
t72 = t123 * t216 + t125 * t191 + t127 * t192;
t89 = t157 * t216 + t158 * t191 + t159 * t192;
t303 = t216 * t72 + t218 * t71 + t271 * t89 + t12;
t15 = t60 * t230 + (t227 * t50 - t229 * t51) * t228;
t302 = t15 + t88 * t230 + (t227 * t69 - t229 * t70) * t228;
t16 = t61 * t230 + (t227 * t52 - t229 * t53) * t228;
t301 = t16 + t89 * t230 + (t227 * t71 - t229 * t72) * t228;
t54 = t101 * t213 + t103 * t193 + t105 * t194;
t55 = t102 * t213 + t104 * t193 + t106 * t194;
t68 = t118 * t213 + t119 * t193 + t120 * t194;
t22 = t216 * t55 + t218 * t54 + t271 * t68;
t78 = t122 * t271 - t124 * t213 + t126 * t214;
t79 = t123 * t271 - t125 * t213 + t127 * t214;
t93 = t157 * t271 - t158 * t213 + t159 * t214;
t300 = t216 * t79 + t218 * t78 + t271 * t93 + t22;
t24 = t68 * t230 + (t227 * t54 - t229 * t55) * t228;
t299 = t24 + t93 * t230 + (t227 * t78 - t229 * t79) * t228;
t107 = rSges(7,1) * t154 + rSges(7,2) * t153 + rSges(7,3) * t189;
t263 = pkin(5) * t190 + pkin(10) * t189 + t107;
t108 = rSges(7,1) * t156 + rSges(7,2) * t155 - rSges(7,3) * t191;
t262 = pkin(5) * t192 - pkin(10) * t191 + t108;
t121 = rSges(7,1) * t194 + rSges(7,2) * t193 + rSges(7,3) * t213;
t298 = pkin(5) * t214 + pkin(10) * t213 + t121;
t297 = t215 * t309 - t216 * t308 + t273 * t306;
t296 = t215 * t308 - t216 * t310 + t273 * t307;
t295 = t217 * t309 - t218 * t308 - t274 * t306;
t294 = -t217 * t308 + t218 * t310 + t274 * t307;
t293 = t215 * t306 - t216 * t307 + t273 * t305;
t292 = -t217 * t306 + t218 * t307 + t274 * t305;
t291 = t305 * t230 + (t233 * t307 + t236 * t306) * t228;
t290 = t306 * t230 + (t233 * t308 + t236 * t309) * t228;
t289 = t307 * t230 + (t233 * t310 + t236 * t308) * t228;
t287 = t189 / 0.2e1;
t286 = -t191 / 0.2e1;
t285 = t213 / 0.2e1;
t284 = t216 / 0.2e1;
t283 = t218 / 0.2e1;
t282 = t227 / 0.2e1;
t281 = -t229 / 0.2e1;
t280 = t230 / 0.2e1;
t235 = cos(qJ(4));
t279 = pkin(4) * t235;
t277 = t262 * t218;
t232 = sin(qJ(4));
t276 = t215 * t232;
t275 = t217 * t232;
t272 = t228 * t232;
t270 = t228 * t235;
t266 = t232 * t236;
t264 = t263 * t271;
t261 = t298 * t216;
t187 = pkin(2) * t216 + qJ(3) * t215;
t188 = pkin(2) * t218 + qJ(3) * t217;
t259 = t187 * t274 + t188 * t273;
t186 = t230 * t188;
t202 = pkin(3) * t274 + t218 * pkin(8);
t258 = t202 * t230 + t186;
t203 = -pkin(3) * t273 + t216 * pkin(8);
t257 = -t187 - t203;
t221 = (pkin(2) * t233 - qJ(3) * t236) * t228;
t256 = -t230 * pkin(3) - pkin(8) * t271 - t221;
t255 = -m(4) - m(5) - m(6) - m(7);
t143 = pkin(4) * t275 + pkin(9) * t218 + t274 * t279;
t254 = t143 * t230 + t258;
t144 = pkin(4) * t276 + pkin(9) * t216 - t273 * t279;
t253 = -t144 + t257;
t185 = t279 * t230 + (-pkin(4) * t266 + pkin(9) * t233) * t228;
t252 = -t185 + t256;
t249 = t271 / 0.2e1;
t247 = (-t230 * rSges(4,1) - (-rSges(4,2) * t233 - rSges(4,3) * t236) * t228 - t221) * t228;
t246 = t202 * t273 + t203 * t274 + t259;
t219 = -t230 * t232 - t235 * t269;
t220 = -t228 * t266 + t230 * t235;
t180 = rSges(5,1) * t220 + rSges(5,2) * t219 + rSges(5,3) * t271;
t245 = (-t180 + t256) * t228;
t18 = t189 * t54 - t191 * t55 + t213 * t68;
t3 = t189 * t50 - t191 * t51 + t213 * t60;
t4 = t189 * t52 - t191 * t53 + t213 * t61;
t244 = t11 * t287 + t12 * t286 + t18 * t249 + t22 * t285 + t283 * t3 + t284 * t4;
t242 = t216 * t303 + t218 * t304 + t271 * t300;
t160 = rSges(6,1) * t214 - rSges(6,2) * t213 + rSges(6,3) * t271;
t241 = (-t160 + t252) * t228;
t240 = t143 * t273 + t144 * t274 + t246;
t239 = (t252 - t298) * t228;
t238 = t301 * t284 + t302 * t283 + t300 * t280 + t304 * t274 / 0.2e1 - t303 * t273 / 0.2e1 + t299 * t249;
t211 = t230 * rSges(3,3) + (rSges(3,1) * t233 + rSges(3,2) * t236) * t228;
t201 = -t229 * t270 + t276;
t200 = t215 * t235 + t229 * t272;
t199 = t227 * t270 + t275;
t198 = t217 * t235 - t227 * t272;
t179 = Icges(5,1) * t220 + Icges(5,4) * t219 + Icges(5,5) * t271;
t178 = Icges(5,4) * t220 + Icges(5,2) * t219 + Icges(5,6) * t271;
t177 = Icges(5,5) * t220 + Icges(5,6) * t219 + Icges(5,3) * t271;
t176 = rSges(3,1) * t218 - rSges(3,2) * t217 + rSges(3,3) * t274;
t175 = rSges(3,1) * t216 - rSges(3,2) * t215 - rSges(3,3) * t273;
t174 = -rSges(4,1) * t273 - rSges(4,2) * t216 + rSges(4,3) * t215;
t173 = rSges(4,1) * t274 - rSges(4,2) * t218 + rSges(4,3) * t217;
t152 = t216 * t185;
t150 = t216 * t160;
t147 = -t175 * t230 - t211 * t273;
t146 = t176 * t230 - t211 * t274;
t142 = rSges(5,1) * t201 + rSges(5,2) * t200 + rSges(5,3) * t216;
t141 = rSges(5,1) * t199 + rSges(5,2) * t198 + rSges(5,3) * t218;
t140 = Icges(5,1) * t201 + Icges(5,4) * t200 + Icges(5,5) * t216;
t139 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t218;
t138 = Icges(5,4) * t201 + Icges(5,2) * t200 + Icges(5,6) * t216;
t137 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t218;
t136 = Icges(5,5) * t201 + Icges(5,6) * t200 + Icges(5,3) * t216;
t135 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t218;
t133 = t143 * t271;
t129 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t216;
t128 = rSges(6,1) * t190 - rSges(6,2) * t189 + rSges(6,3) * t218;
t117 = t218 * t144;
t116 = t128 * t271;
t115 = (t175 * t227 + t176 * t229) * t228;
t114 = t218 * t129;
t112 = (-t174 - t187) * t230 + t229 * t247;
t111 = t230 * t173 + t227 * t247 + t186;
t110 = t141 * t271 - t180 * t218;
t109 = -t142 * t271 + t180 * t216;
t99 = -t160 * t218 + t116;
t98 = -t129 * t271 + t150;
t96 = t177 * t271 + t178 * t219 + t179 * t220;
t95 = (t173 * t229 + t174 * t227) * t228 + t259;
t94 = -t141 * t216 + t142 * t218;
t92 = t177 * t216 + t178 * t200 + t179 * t201;
t91 = t177 * t218 + t178 * t198 + t179 * t199;
t90 = -t128 * t216 + t114;
t87 = (-t142 + t257) * t230 + t229 * t245;
t86 = t230 * t141 + t227 * t245 + t258;
t85 = t136 * t271 + t138 * t219 + t140 * t220;
t84 = t135 * t271 + t137 * t219 + t139 * t220;
t83 = -t108 * t213 - t121 * t191;
t82 = t107 * t213 - t121 * t189;
t81 = t116 + t133 + (-t160 - t185) * t218;
t80 = t150 + t152 + (-t129 - t144) * t271;
t77 = t136 * t216 + t138 * t200 + t140 * t201;
t76 = t135 * t216 + t137 * t200 + t139 * t201;
t75 = t136 * t218 + t138 * t198 + t140 * t199;
t74 = t135 * t218 + t137 * t198 + t139 * t199;
t73 = (t141 * t229 + t142 * t227) * t228 + t246;
t67 = t107 * t191 + t108 * t189;
t66 = (-t129 + t253) * t230 + t229 * t241;
t65 = t230 * t128 + t227 * t241 + t254;
t64 = t114 + t117 + (-t128 - t143) * t216;
t63 = -t218 * t298 + t264;
t62 = -t262 * t271 + t261;
t59 = (t128 * t229 + t129 * t227) * t228 + t240;
t58 = -t216 * t263 + t277;
t57 = t133 + (-t185 - t298) * t218 + t264;
t56 = t152 + (-t144 - t262) * t271 + t261;
t49 = (t253 - t262) * t230 + t229 * t239;
t48 = t227 * t239 + t230 * t263 + t254;
t47 = t117 + (-t143 - t263) * t216 + t277;
t46 = t96 * t230 + (t227 * t84 - t229 * t85) * t228;
t45 = t216 * t85 + t218 * t84 + t271 * t96;
t44 = (t227 * t262 + t229 * t263) * t228 + t240;
t41 = t92 * t230 + (t227 * t76 - t229 * t77) * t228;
t40 = t91 * t230 + (t227 * t74 - t229 * t75) * t228;
t36 = t216 * t77 + t218 * t76 + t271 * t92;
t35 = t216 * t75 + t218 * t74 + t271 * t91;
t1 = [m(2) + m(3) - t255; m(3) * t115 + m(4) * t95 + m(5) * t73 + m(6) * t59 + m(7) * t44; m(7) * (t44 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t59 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t73 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2 + t95 ^ 2) + m(3) * (t115 ^ 2 + t146 ^ 2 + t147 ^ 2) + (t46 + t291 * t230 ^ 2 + ((t233 * t289 + t236 * t290) * t230 + (t293 * t230 + (t233 * t296 + t236 * t297) * t228) * t229 + (t292 * t230 + (t233 * t294 - t236 * t295) * t228) * t227) * t228 + t299) * t230 + (t40 + (t217 * t295 + t218 * t294 + t292 * t274) * t274 + (-t217 * t290 + t218 * t289 + t274 * t291) * t230 + t302) * t274 + (-t41 + (t215 * t297 - t216 * t296 + t293 * t273) * t273 + (t215 * t290 - t216 * t289 + t273 * t291) * t230 + (-t215 * t295 - t216 * t294 - t217 * t297 + t218 * t296 + t292 * t273 + t293 * t274) * t274 - t301) * t273; t255 * t269; m(7) * (t215 * t48 + t217 * t49 - t269 * t44) + m(6) * (t215 * t65 + t217 * t66 - t269 * t59) + m(5) * (t215 * t86 + t217 * t87 - t269 * t73) + m(4) * (t111 * t215 + t112 * t217 - t269 * t95); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t228 ^ 2 * t236 ^ 2 + t215 ^ 2 + t217 ^ 2); m(5) * t94 + m(6) * t64 + m(7) * t47; (t233 * t46 / 0.2e1 + t35 * t282 + t36 * t281) * t228 + t238 + m(7) * (t44 * t47 + t48 * t57 + t49 * t56) + m(6) * (t59 * t64 + t65 * t81 + t66 * t80) + m(5) * (t109 * t87 + t110 * t86 + t73 * t94) + t41 * t284 + t40 * t283 + t45 * t280; m(5) * (t109 * t217 + t110 * t215 - t269 * t94) + m(6) * (t215 * t81 + t217 * t80 - t269 * t64) + m(7) * (t215 * t57 + t217 * t56 - t269 * t47); t45 * t271 + t216 * t36 + t218 * t35 + m(7) * (t47 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t64 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2 + t94 ^ 2) + t242; m(6) * t90 + m(7) * t58; m(7) * (t44 * t58 + t48 * t63 + t49 * t62) + m(6) * (t59 * t90 + t65 * t99 + t66 * t98) + t238; m(6) * (t215 * t99 + t217 * t98 - t269 * t90) + m(7) * (t215 * t63 + t217 * t62 - t269 * t58); m(7) * (t47 * t58 + t56 * t62 + t57 * t63) + m(6) * (t64 * t90 + t80 * t98 + t81 * t99) + t242; m(7) * (t58 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t90 ^ 2 + t98 ^ 2 + t99 ^ 2) + t242; m(7) * t67; m(7) * (t44 * t67 + t48 * t82 + t49 * t83) + t24 * t285 + t16 * t286 + t15 * t287 + t18 * t280 + (t281 * t4 + t282 * t3) * t228; m(7) * (t215 * t82 + t217 * t83 - t269 * t67); m(7) * (t47 * t67 + t56 * t83 + t57 * t82) + t244; m(7) * (t58 * t67 + t62 * t83 + t63 * t82) + t244; t189 * t3 - t191 * t4 + t213 * t18 + m(7) * (t67 ^ 2 + t82 ^ 2 + t83 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
