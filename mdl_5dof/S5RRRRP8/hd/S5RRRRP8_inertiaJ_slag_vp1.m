% Calculate joint inertia matrix for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:59:41
% EndTime: 2019-12-31 21:59:51
% DurationCPUTime: 4.06s
% Computational Cost: add. (7557->399), mult. (10679->559), div. (0->0), fcn. (11570->8), ass. (0->208)
t212 = qJ(3) + qJ(4);
t201 = sin(t212);
t202 = cos(t212);
t218 = cos(qJ(1));
t215 = sin(qJ(1));
t217 = cos(qJ(2));
t259 = t215 * t217;
t166 = -t201 * t259 - t202 * t218;
t167 = -t201 * t218 + t202 * t259;
t214 = sin(qJ(2));
t262 = t214 * t215;
t100 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t262;
t102 = Icges(5,5) * t167 + Icges(5,6) * t166 + Icges(5,3) * t262;
t312 = t100 + t102;
t258 = t217 * t218;
t168 = -t201 * t258 + t202 * t215;
t169 = t201 * t215 + t202 * t258;
t261 = t214 * t218;
t101 = Icges(6,5) * t169 + Icges(6,6) * t168 + Icges(6,3) * t261;
t103 = Icges(5,5) * t169 + Icges(5,6) * t168 + Icges(5,3) * t261;
t311 = t101 + t103;
t104 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t262;
t106 = Icges(5,4) * t167 + Icges(5,2) * t166 + Icges(5,6) * t262;
t310 = t104 + t106;
t105 = Icges(6,4) * t169 + Icges(6,2) * t168 + Icges(6,6) * t261;
t107 = Icges(5,4) * t169 + Icges(5,2) * t168 + Icges(5,6) * t261;
t309 = t105 + t107;
t108 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t262;
t110 = Icges(5,1) * t167 + Icges(5,4) * t166 + Icges(5,5) * t262;
t308 = t108 + t110;
t109 = Icges(6,1) * t169 + Icges(6,4) * t168 + Icges(6,5) * t261;
t111 = Icges(5,1) * t169 + Icges(5,4) * t168 + Icges(5,5) * t261;
t307 = t109 + t111;
t306 = t310 * t166 + t308 * t167 + t312 * t262;
t305 = t309 * t166 + t307 * t167 + t311 * t262;
t304 = t310 * t168 + t308 * t169 + t312 * t261;
t303 = t309 * t168 + t307 * t169 + t311 * t261;
t144 = -Icges(6,3) * t217 + (Icges(6,5) * t202 - Icges(6,6) * t201) * t214;
t146 = -Icges(6,6) * t217 + (Icges(6,4) * t202 - Icges(6,2) * t201) * t214;
t148 = -Icges(6,5) * t217 + (Icges(6,1) * t202 - Icges(6,4) * t201) * t214;
t67 = t144 * t262 + t146 * t166 + t148 * t167;
t145 = -Icges(5,3) * t217 + (Icges(5,5) * t202 - Icges(5,6) * t201) * t214;
t147 = -Icges(5,6) * t217 + (Icges(5,4) * t202 - Icges(5,2) * t201) * t214;
t149 = -Icges(5,5) * t217 + (Icges(5,1) * t202 - Icges(5,4) * t201) * t214;
t68 = t145 * t262 + t147 * t166 + t149 * t167;
t302 = -t68 - t67;
t69 = t144 * t261 + t146 * t168 + t148 * t169;
t70 = t145 * t261 + t147 * t168 + t149 * t169;
t301 = -t69 - t70;
t296 = (t148 + t149) * t202 * t214;
t298 = (-t146 - t147) * t201;
t299 = -t144 - t145;
t287 = t214 * t298 + t299 * t217 + t296;
t300 = t287 * t217;
t297 = Icges(3,5) * t214;
t213 = sin(qJ(3));
t182 = pkin(3) * t213 + pkin(4) * t201;
t295 = -rSges(6,1) * t167 - rSges(6,2) * t166 + t182 * t218;
t294 = t297 / 0.2e1;
t293 = t302 * t217 + (t306 * t215 + t305 * t218) * t214;
t292 = t301 * t217 + (t304 * t215 + t303 * t218) * t214;
t291 = t305 * t215 - t306 * t218;
t290 = t303 * t215 - t304 * t218;
t50 = -t217 * t100 + (-t104 * t201 + t108 * t202) * t214;
t52 = -t217 * t102 + (-t106 * t201 + t110 * t202) * t214;
t289 = -t50 - t52;
t51 = -t217 * t101 + (-t105 * t201 + t109 * t202) * t214;
t53 = -t217 * t103 + (-t107 * t201 + t111 * t202) * t214;
t288 = t51 + t53;
t219 = -pkin(8) - pkin(7);
t208 = -qJ(5) + t219;
t260 = t214 * t219;
t263 = t213 * t218;
t248 = pkin(3) * t263 + t215 * t260;
t216 = cos(qJ(3));
t200 = t216 * pkin(3) + pkin(2);
t180 = pkin(4) * t202 + t200;
t250 = t180 - t200;
t286 = rSges(6,3) * t262 + (-t208 * t214 + t250 * t217) * t215 + t248 - t295;
t245 = t208 - t219;
t285 = (t245 - rSges(6,3)) * t217 + (rSges(6,1) * t202 - rSges(6,2) * t201 + t250) * t214;
t284 = t169 * rSges(6,1) + t168 * rSges(6,2) + rSges(6,3) * t261 + t180 * t258 + t215 * t182;
t210 = t215 ^ 2;
t211 = t218 ^ 2;
t283 = t215 / 0.2e1;
t282 = -t217 / 0.2e1;
t281 = -t218 / 0.2e1;
t280 = t218 / 0.2e1;
t187 = rSges(3,1) * t214 + rSges(3,2) * t217;
t279 = m(3) * t187;
t278 = pkin(2) * t217;
t277 = pkin(7) * t214;
t276 = -pkin(2) + t200;
t275 = t300 + (t289 * t215 - t288 * t218) * t214;
t273 = t286 * t261;
t272 = t218 * rSges(3,3);
t264 = t213 * t215;
t249 = -pkin(3) * t264 - t200 * t258;
t271 = -t245 * t261 + t249 + t284;
t231 = -rSges(5,1) * t167 - rSges(5,2) * t166;
t113 = rSges(5,3) * t262 - t231;
t151 = -t217 * rSges(5,3) + (rSges(5,1) * t202 - rSges(5,2) * t201) * t214;
t83 = t217 * t113 + t151 * t262;
t269 = Icges(3,4) * t217;
t159 = -Icges(4,6) * t217 + (Icges(4,4) * t216 - Icges(4,2) * t213) * t214;
t265 = t213 * t159;
t115 = t169 * rSges(5,1) + t168 * rSges(5,2) + rSges(5,3) * t261;
t223 = -t218 * t260 - t249;
t247 = pkin(2) * t258 + pkin(7) * t261;
t131 = t223 - t247;
t257 = -t115 - t131;
t130 = (t276 * t217 - t277) * t215 - t248;
t143 = (pkin(7) + t219) * t217 + t276 * t214;
t256 = t217 * t130 + t143 * t262;
t254 = -t143 - t151;
t165 = -t217 * rSges(4,3) + (rSges(4,1) * t216 - rSges(4,2) * t213) * t214;
t190 = t214 * pkin(2) - t217 * pkin(7);
t253 = -t165 - t190;
t251 = t210 * (t277 + t278) + t218 * t247;
t246 = t218 * pkin(1) + t215 * pkin(6);
t244 = t210 + t211;
t175 = -t213 * t259 - t216 * t218;
t176 = t216 * t259 - t263;
t122 = Icges(4,5) * t176 + Icges(4,6) * t175 + Icges(4,3) * t262;
t124 = Icges(4,4) * t176 + Icges(4,2) * t175 + Icges(4,6) * t262;
t126 = Icges(4,1) * t176 + Icges(4,4) * t175 + Icges(4,5) * t262;
t60 = -t217 * t122 + (-t124 * t213 + t126 * t216) * t214;
t156 = -Icges(4,3) * t217 + (Icges(4,5) * t216 - Icges(4,6) * t213) * t214;
t162 = -Icges(4,5) * t217 + (Icges(4,1) * t216 - Icges(4,4) * t213) * t214;
t76 = t156 * t262 + t159 * t175 + t162 * t176;
t243 = t60 / 0.2e1 + t76 / 0.2e1;
t177 = -t213 * t258 + t215 * t216;
t178 = t216 * t258 + t264;
t123 = Icges(4,5) * t178 + Icges(4,6) * t177 + Icges(4,3) * t261;
t125 = Icges(4,4) * t178 + Icges(4,2) * t177 + Icges(4,6) * t261;
t127 = Icges(4,1) * t178 + Icges(4,4) * t177 + Icges(4,5) * t261;
t61 = -t217 * t123 + (-t125 * t213 + t127 * t216) * t214;
t77 = t156 * t261 + t159 * t177 + t162 * t178;
t242 = t61 / 0.2e1 + t77 / 0.2e1;
t241 = -t131 - t271;
t240 = t292 * t261 + t293 * t262;
t239 = -t143 - t285;
t238 = -t190 + t254;
t129 = t178 * rSges(4,1) + t177 * rSges(4,2) + rSges(4,3) * t261;
t237 = t262 / 0.2e1;
t236 = t261 / 0.2e1;
t36 = t286 * t217 + t285 * t262;
t235 = t215 * t130 + t218 * t131 + t251;
t234 = -t190 + t239;
t233 = rSges(3,1) * t217 - rSges(3,2) * t214;
t232 = -rSges(4,1) * t176 - rSges(4,2) * t175;
t228 = -Icges(3,2) * t214 + t269;
t227 = Icges(3,5) * t217 - Icges(3,6) * t214;
t224 = rSges(3,1) * t258 - rSges(3,2) * t261 + t215 * rSges(3,3);
t222 = t275 * t217 + t240;
t221 = (-t289 - t302) * t237 + (t288 - t301) * t236;
t220 = t292 * t283 + (t288 * t215 + t289 * t218) * t282 + t293 * t281 + t291 * t237 + t290 * t236;
t206 = t218 * pkin(6);
t189 = rSges(2,1) * t218 - rSges(2,2) * t215;
t188 = -rSges(2,1) * t215 - rSges(2,2) * t218;
t184 = Icges(3,6) * t217 + t297;
t158 = Icges(3,3) * t215 + t227 * t218;
t157 = -Icges(3,3) * t218 + t227 * t215;
t142 = t214 * t216 * t162;
t141 = t224 + t246;
t140 = t272 + t206 + (-pkin(1) - t233) * t215;
t134 = t253 * t218;
t133 = t253 * t215;
t128 = rSges(4,3) * t262 - t232;
t117 = t218 * t224 + (t233 * t215 - t272) * t215;
t116 = t130 * t261;
t97 = t113 * t261;
t95 = t129 + t246 + t247;
t94 = t206 + (-t278 - pkin(1) + (-rSges(4,3) - pkin(7)) * t214) * t215 + t232;
t93 = t238 * t218;
t92 = t238 * t215;
t89 = -t129 * t217 - t165 * t261;
t88 = t128 * t217 + t165 * t262;
t86 = -t217 * t156 - t214 * t265 + t142;
t84 = -t115 * t217 - t151 * t261;
t82 = t223 + t115 + t246;
t81 = t206 + (-rSges(5,3) * t214 - t200 * t217 - pkin(1)) * t215 + t231 + t248;
t78 = (t128 * t218 - t129 * t215) * t214;
t75 = -t208 * t261 + t246 + t284;
t74 = t206 + (-t180 * t217 - pkin(1) + (-rSges(6,3) + t208) * t214) * t215 + t295;
t73 = t234 * t218;
t72 = t234 * t215;
t71 = -t115 * t262 + t97;
t62 = t128 * t215 + t129 * t218 + t251;
t59 = t257 * t217 + t254 * t261;
t58 = t256 + t83;
t57 = t123 * t261 + t125 * t177 + t127 * t178;
t56 = t122 * t261 + t124 * t177 + t126 * t178;
t55 = t123 * t262 + t125 * t175 + t127 * t176;
t54 = t122 * t262 + t124 * t175 + t126 * t176;
t37 = -t271 * t217 - t261 * t285;
t35 = t257 * t262 + t116 + t97;
t34 = t113 * t215 + t115 * t218 + t235;
t33 = -t271 * t262 + t273;
t32 = t241 * t217 + t239 * t261;
t31 = t36 + t256;
t30 = t215 * t57 - t218 * t56;
t29 = t215 * t55 - t218 * t54;
t26 = t241 * t262 + t116 + t273;
t25 = t286 * t215 + t271 * t218 + t235;
t16 = -t77 * t217 + (t215 * t56 + t218 * t57) * t214;
t15 = -t76 * t217 + (t215 * t54 + t218 * t55) * t214;
t1 = [Icges(2,3) + t142 + (Icges(3,4) * t214 + Icges(3,2) * t217 - t156 + t299) * t217 + (Icges(3,1) * t214 - t265 + t269 + t298) * t214 + m(6) * (t74 ^ 2 + t75 ^ 2) + m(5) * (t81 ^ 2 + t82 ^ 2) + m(4) * (t94 ^ 2 + t95 ^ 2) + m(3) * (t140 ^ 2 + t141 ^ 2) + m(2) * (t188 ^ 2 + t189 ^ 2) + t296; m(6) * (t72 * t75 + t73 * t74) + m(5) * (t81 * t93 + t82 * t92) + m(4) * (t133 * t95 + t134 * t94) + (-t52 / 0.2e1 - t50 / 0.2e1 - t68 / 0.2e1 - t67 / 0.2e1 + (-Icges(3,6) * t218 + t228 * t215) * t282 + t218 * t294 - t140 * t279 + t184 * t280 - t243) * t218 + (t70 / 0.2e1 + t69 / 0.2e1 + t53 / 0.2e1 + t51 / 0.2e1 + t217 * (Icges(3,6) * t215 + t228 * t218) / 0.2e1 + t215 * t294 - t141 * t279 + t184 * t283 + t242) * t215; m(6) * (t25 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(5) * (t34 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(4) * (t133 ^ 2 + t134 ^ 2 + t62 ^ 2) + m(3) * (t244 * t187 ^ 2 + t117 ^ 2) + (-t211 * t157 - t29 - t291) * t218 + (t210 * t158 + t30 + (-t215 * t157 + t218 * t158) * t218 + t290) * t215; (-t86 - t287) * t217 + m(6) * (t31 * t74 + t32 * t75) + m(5) * (t58 * t81 + t59 * t82) + m(4) * (t88 * t94 + t89 * t95) + (t243 * t215 + t242 * t218) * t214 + t221; m(6) * (t25 * t26 + t31 * t73 + t32 * t72) + m(5) * (t35 * t34 + t58 * t93 + t59 * t92) + m(4) * (t133 * t89 + t134 * t88 + t62 * t78) + (t30 * t280 + t29 * t283) * t214 + t220 + t16 * t283 + (t61 * t215 - t60 * t218) * t282 + t15 * t281; (t86 * t217 + t275) * t217 + (t218 * t16 + t215 * t15 - t217 * (t215 * t60 + t218 * t61)) * t214 + m(6) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t35 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(4) * (t78 ^ 2 + t88 ^ 2 + t89 ^ 2) + t240; -t300 + m(6) * (t36 * t74 + t37 * t75) + m(5) * (t81 * t83 + t82 * t84) + t221; m(6) * (t25 * t33 + t36 * t73 + t37 * t72) + m(5) * (t71 * t34 + t83 * t93 + t84 * t92) + t220; m(6) * (t26 * t33 + t31 * t36 + t32 * t37) + m(5) * (t71 * t35 + t58 * t83 + t59 * t84) + t222; m(6) * (t33 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t71 ^ 2 + t83 ^ 2 + t84 ^ 2) + t222; m(6) * (t215 * t75 + t218 * t74) * t214; m(6) * (-t217 * t25 + (t215 * t72 + t218 * t73) * t214); m(6) * (-t217 * t26 + (t215 * t32 + t218 * t31) * t214); m(6) * (-t217 * t33 + (t215 * t37 + t218 * t36) * t214); m(6) * (t244 * t214 ^ 2 + t217 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
