% Calculate joint inertia matrix for
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP9_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP9_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:49
% EndTime: 2019-12-31 22:03:59
% DurationCPUTime: 3.87s
% Computational Cost: add. (7103->385), mult. (10320->547), div. (0->0), fcn. (11330->8), ass. (0->208)
t206 = qJ(3) + qJ(4);
t200 = sin(t206);
t201 = cos(t206);
t212 = cos(qJ(1));
t253 = t212 * t201;
t209 = sin(qJ(1));
t211 = cos(qJ(2));
t256 = t209 * t211;
t167 = t200 * t256 + t253;
t254 = t212 * t200;
t168 = t201 * t256 - t254;
t208 = sin(qJ(2));
t260 = t208 * t209;
t100 = Icges(6,4) * t168 + Icges(6,2) * t260 + Icges(6,6) * t167;
t98 = Icges(5,5) * t168 - Icges(5,6) * t167 + Icges(5,3) * t260;
t307 = t100 + t98;
t169 = -t201 * t209 + t211 * t254;
t170 = t200 * t209 + t211 * t253;
t259 = t208 * t212;
t101 = Icges(6,4) * t170 + Icges(6,2) * t259 + Icges(6,6) * t169;
t99 = Icges(5,5) * t170 - Icges(5,6) * t169 + Icges(5,3) * t259;
t306 = t101 + t99;
t102 = Icges(5,4) * t168 - Icges(5,2) * t167 + Icges(5,6) * t260;
t96 = Icges(6,5) * t168 + Icges(6,6) * t260 + Icges(6,3) * t167;
t305 = -t102 + t96;
t103 = Icges(5,4) * t170 - Icges(5,2) * t169 + Icges(5,6) * t259;
t97 = Icges(6,5) * t170 + Icges(6,6) * t259 + Icges(6,3) * t169;
t304 = -t103 + t97;
t104 = Icges(6,1) * t168 + Icges(6,4) * t260 + Icges(6,5) * t167;
t106 = Icges(5,1) * t168 - Icges(5,4) * t167 + Icges(5,5) * t260;
t303 = t104 + t106;
t105 = Icges(6,1) * t170 + Icges(6,4) * t259 + Icges(6,5) * t169;
t107 = Icges(5,1) * t170 - Icges(5,4) * t169 + Icges(5,5) * t259;
t302 = t105 + t107;
t282 = rSges(6,3) + qJ(5);
t289 = rSges(6,1) + pkin(4);
t301 = -t282 * t167 - t289 * t168;
t300 = t167 * t305 + t168 * t303 + t260 * t307;
t299 = t167 * t304 + t168 * t302 + t260 * t306;
t298 = t169 * t305 + t170 * t303 + t259 * t307;
t297 = t169 * t304 + t170 * t302 + t259 * t306;
t144 = -Icges(6,6) * t211 + (Icges(6,5) * t201 + Icges(6,3) * t200) * t208;
t146 = -Icges(6,2) * t211 + (Icges(6,4) * t201 + Icges(6,6) * t200) * t208;
t148 = -Icges(6,4) * t211 + (Icges(6,1) * t201 + Icges(6,5) * t200) * t208;
t69 = t144 * t167 + t146 * t260 + t148 * t168;
t145 = -Icges(5,3) * t211 + (Icges(5,5) * t201 - Icges(5,6) * t200) * t208;
t147 = -Icges(5,6) * t211 + (Icges(5,4) * t201 - Icges(5,2) * t200) * t208;
t149 = -Icges(5,5) * t211 + (Icges(5,1) * t201 - Icges(5,4) * t200) * t208;
t70 = t145 * t260 - t147 * t167 + t149 * t168;
t296 = -t70 - t69;
t71 = t144 * t169 + t146 * t259 + t148 * t170;
t72 = t145 * t259 - t147 * t169 + t149 * t170;
t295 = -t71 - t72;
t294 = t145 + t146;
t293 = Icges(3,5) * t208;
t263 = t200 * t208;
t291 = t144 * t263 + (t148 + t149) * t201 * t208;
t268 = t147 * t263 + t211 * t294 - t291;
t292 = t268 * t211;
t290 = t293 / 0.2e1;
t288 = t296 * t211 + (t209 * t300 + t212 * t299) * t208;
t287 = t295 * t211 + (t209 * t298 + t212 * t297) * t208;
t286 = t209 * t299 - t212 * t300;
t285 = t209 * t297 - t212 * t298;
t48 = -t211 * t100 + (t104 * t201 + t200 * t96) * t208;
t50 = -t211 * t98 + (-t102 * t200 + t106 * t201) * t208;
t284 = -t48 - t50;
t49 = -t211 * t101 + (t105 * t201 + t200 * t97) * t208;
t51 = -t211 * t99 + (-t103 * t200 + t107 * t201) * t208;
t283 = t49 + t51;
t281 = rSges(6,2) * t260 - t301;
t280 = -t211 * rSges(6,2) + (t200 * t282 + t201 * t289) * t208;
t279 = t209 ^ 2;
t278 = t212 ^ 2;
t277 = t209 / 0.2e1;
t276 = -t211 / 0.2e1;
t275 = -t212 / 0.2e1;
t274 = t212 / 0.2e1;
t185 = rSges(3,1) * t208 + rSges(3,2) * t211;
t273 = m(3) * t185;
t272 = pkin(2) * t211;
t271 = pkin(7) * t208;
t210 = cos(qJ(3));
t199 = pkin(3) * t210 + pkin(2);
t270 = -pkin(2) + t199;
t269 = -t292 + (t209 * t284 - t212 * t283) * t208;
t267 = t212 * rSges(3,3);
t266 = t281 * t259;
t225 = -t168 * rSges(5,1) + t167 * rSges(5,2);
t109 = rSges(5,3) * t260 - t225;
t151 = -t211 * rSges(5,3) + (rSges(5,1) * t201 - rSges(5,2) * t200) * t208;
t83 = t109 * t211 + t151 * t260;
t264 = Icges(3,4) * t211;
t207 = sin(qJ(3));
t160 = -Icges(4,6) * t211 + (Icges(4,4) * t210 - Icges(4,2) * t207) * t208;
t261 = t207 * t160;
t213 = -pkin(8) - pkin(7);
t258 = t208 * t213;
t257 = t209 * t207;
t255 = t211 * t212;
t252 = t212 * t207;
t251 = t212 * t210;
t250 = rSges(6,2) * t259 + t169 * t282 + t170 * t289;
t111 = rSges(5,1) * t170 - rSges(5,2) * t169 + rSges(5,3) * t259;
t218 = pkin(3) * t257 + t199 * t255 - t212 * t258;
t241 = pkin(2) * t255 + pkin(7) * t259;
t130 = t218 - t241;
t249 = -t111 - t130;
t242 = -pkin(3) * t252 - t209 * t258;
t129 = (t211 * t270 - t271) * t209 + t242;
t143 = (pkin(7) + t213) * t211 + t270 * t208;
t248 = t129 * t211 + t143 * t260;
t246 = -t143 - t151;
t166 = -t211 * rSges(4,3) + (rSges(4,1) * t210 - rSges(4,2) * t207) * t208;
t188 = t208 * pkin(2) - t211 * pkin(7);
t244 = -t166 - t188;
t243 = t279 * (t271 + t272) + t212 * t241;
t240 = pkin(1) * t212 + pkin(6) * t209;
t176 = -t207 * t256 - t251;
t177 = t210 * t256 - t252;
t121 = Icges(4,5) * t177 + Icges(4,6) * t176 + Icges(4,3) * t260;
t123 = Icges(4,4) * t177 + Icges(4,2) * t176 + Icges(4,6) * t260;
t125 = Icges(4,1) * t177 + Icges(4,4) * t176 + Icges(4,5) * t260;
t60 = -t211 * t121 + (-t123 * t207 + t125 * t210) * t208;
t157 = -Icges(4,3) * t211 + (Icges(4,5) * t210 - Icges(4,6) * t207) * t208;
t163 = -Icges(4,5) * t211 + (Icges(4,1) * t210 - Icges(4,4) * t207) * t208;
t74 = t157 * t260 + t160 * t176 + t163 * t177;
t239 = t60 / 0.2e1 + t74 / 0.2e1;
t178 = t209 * t210 - t211 * t252;
t179 = t211 * t251 + t257;
t122 = Icges(4,5) * t179 + Icges(4,6) * t178 + Icges(4,3) * t259;
t124 = Icges(4,4) * t179 + Icges(4,2) * t178 + Icges(4,6) * t259;
t126 = Icges(4,1) * t179 + Icges(4,4) * t178 + Icges(4,5) * t259;
t61 = -t211 * t122 + (-t124 * t207 + t126 * t210) * t208;
t75 = t157 * t259 + t160 * t178 + t163 * t179;
t238 = t75 / 0.2e1 + t61 / 0.2e1;
t237 = t259 * t287 + t260 * t288;
t236 = -t130 - t250;
t235 = -t143 - t280;
t234 = -t188 + t246;
t128 = rSges(4,1) * t179 + rSges(4,2) * t178 + rSges(4,3) * t259;
t204 = t212 * pkin(6);
t233 = t204 - t242;
t232 = t260 / 0.2e1;
t231 = t259 / 0.2e1;
t230 = -t199 * t211 - pkin(1);
t58 = t211 * t281 + t260 * t280;
t229 = t129 * t209 + t130 * t212 + t243;
t228 = -t188 + t235;
t227 = rSges(3,1) * t211 - rSges(3,2) * t208;
t226 = -t177 * rSges(4,1) - t176 * rSges(4,2);
t223 = -Icges(3,2) * t208 + t264;
t222 = Icges(3,5) * t211 - Icges(3,6) * t208;
t219 = rSges(3,1) * t255 - rSges(3,2) * t259 + rSges(3,3) * t209;
t217 = t211 * t269 + t237;
t216 = (-t284 - t296) * t232 + (t283 - t295) * t231;
t215 = t218 + t240;
t214 = t287 * t277 + (t209 * t283 + t212 * t284) * t276 + t288 * t275 + t286 * t232 + t285 * t231;
t187 = rSges(2,1) * t212 - rSges(2,2) * t209;
t186 = -rSges(2,1) * t209 - rSges(2,2) * t212;
t182 = Icges(3,6) * t211 + t293;
t159 = Icges(3,3) * t209 + t212 * t222;
t158 = -Icges(3,3) * t212 + t209 * t222;
t141 = t208 * t210 * t163;
t140 = t219 + t240;
t139 = t267 + t204 + (-pkin(1) - t227) * t209;
t132 = t244 * t212;
t131 = t244 * t209;
t127 = rSges(4,3) * t260 - t226;
t115 = t212 * t219 + (t209 * t227 - t267) * t209;
t113 = t129 * t259;
t93 = t109 * t259;
t91 = t128 + t240 + t241;
t90 = t204 + (-t272 - pkin(1) + (-rSges(4,3) - pkin(7)) * t208) * t209 + t226;
t89 = t234 * t212;
t88 = t234 * t209;
t87 = -t128 * t211 - t166 * t259;
t86 = t127 * t211 + t166 * t260;
t85 = -t211 * t157 - t208 * t261 + t141;
t84 = -t211 * t111 - t151 * t259;
t82 = t215 + t111;
t81 = (-rSges(5,3) * t208 + t230) * t209 + t225 + t233;
t80 = t228 * t212;
t79 = t228 * t209;
t76 = (t127 * t212 - t128 * t209) * t208;
t73 = -t111 * t260 + t93;
t64 = t127 * t209 + t128 * t212 + t243;
t63 = t215 + t250;
t62 = (-rSges(6,2) * t208 + t230) * t209 + t233 + t301;
t59 = -t211 * t250 - t259 * t280;
t57 = t211 * t249 + t246 * t259;
t56 = t248 + t83;
t55 = t122 * t259 + t124 * t178 + t126 * t179;
t54 = t121 * t259 + t123 * t178 + t125 * t179;
t53 = t122 * t260 + t124 * t176 + t126 * t177;
t52 = t121 * t260 + t123 * t176 + t125 * t177;
t35 = t249 * t260 + t113 + t93;
t34 = -t250 * t260 + t266;
t33 = t211 * t236 + t235 * t259;
t32 = t58 + t248;
t31 = t109 * t209 + t111 * t212 + t229;
t30 = t236 * t260 + t113 + t266;
t29 = t209 * t281 + t212 * t250 + t229;
t28 = t209 * t55 - t212 * t54;
t27 = t209 * t53 - t212 * t52;
t16 = -t75 * t211 + (t209 * t54 + t212 * t55) * t208;
t15 = -t74 * t211 + (t209 * t52 + t212 * t53) * t208;
t1 = [Icges(2,3) + t141 + (Icges(3,1) * t208 - t147 * t200 - t261 + t264) * t208 + (Icges(3,4) * t208 + Icges(3,2) * t211 - t157 - t294) * t211 + m(6) * (t62 ^ 2 + t63 ^ 2) + m(5) * (t81 ^ 2 + t82 ^ 2) + m(4) * (t90 ^ 2 + t91 ^ 2) + m(3) * (t139 ^ 2 + t140 ^ 2) + m(2) * (t186 ^ 2 + t187 ^ 2) + t291; m(6) * (t62 * t80 + t63 * t79) + m(5) * (t81 * t89 + t82 * t88) + m(4) * (t131 * t91 + t132 * t90) + (-t70 / 0.2e1 - t69 / 0.2e1 - t50 / 0.2e1 - t48 / 0.2e1 + (-Icges(3,6) * t212 + t209 * t223) * t276 + t212 * t290 - t139 * t273 + t182 * t274 - t239) * t212 + (t72 / 0.2e1 + t71 / 0.2e1 + t51 / 0.2e1 + t49 / 0.2e1 + t211 * (Icges(3,6) * t209 + t212 * t223) / 0.2e1 + t209 * t290 - t140 * t273 + t182 * t277 + t238) * t209; m(6) * (t29 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t31 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t131 ^ 2 + t132 ^ 2 + t64 ^ 2) + m(3) * (t115 ^ 2 + (t278 + t279) * t185 ^ 2) + (-t278 * t158 - t27 - t286) * t212 + (t279 * t159 + t28 + (-t209 * t158 + t212 * t159) * t212 + t285) * t209; (-t85 + t268) * t211 + m(6) * (t32 * t62 + t33 * t63) + m(5) * (t56 * t81 + t57 * t82) + m(4) * (t86 * t90 + t87 * t91) + (t209 * t239 + t212 * t238) * t208 + t216; t214 + m(6) * (t29 * t30 + t32 * t80 + t33 * t79) + m(5) * (t31 * t35 + t56 * t89 + t57 * t88) + m(4) * (t131 * t87 + t132 * t86 + t64 * t76) + (t27 * t277 + t274 * t28) * t208 + t16 * t277 + t15 * t275 + (t61 * t209 - t60 * t212) * t276; (t85 * t211 + t269) * t211 + (t212 * t16 + t209 * t15 - t211 * (t209 * t60 + t212 * t61)) * t208 + m(6) * (t30 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(5) * (t35 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(4) * (t76 ^ 2 + t86 ^ 2 + t87 ^ 2) + t237; t292 + m(6) * (t58 * t62 + t59 * t63) + m(5) * (t81 * t83 + t82 * t84) + t216; m(6) * (t29 * t34 + t58 * t80 + t59 * t79) + m(5) * (t31 * t73 + t83 * t89 + t84 * t88) + t214; m(6) * (t30 * t34 + t32 * t58 + t33 * t59) + m(5) * (t35 * t73 + t56 * t83 + t57 * t84) + t217; m(6) * (t34 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t73 ^ 2 + t83 ^ 2 + t84 ^ 2) + t217; m(6) * (t167 * t63 + t169 * t62); m(6) * (t167 * t79 + t169 * t80 + t263 * t29); m(6) * (t167 * t33 + t169 * t32 + t263 * t30); m(6) * (t167 * t59 + t169 * t58 + t263 * t34); m(6) * (t200 ^ 2 * t208 ^ 2 + t167 ^ 2 + t169 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
