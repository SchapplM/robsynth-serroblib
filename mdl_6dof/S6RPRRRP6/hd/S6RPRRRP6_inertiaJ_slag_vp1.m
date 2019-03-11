% Calculate joint inertia matrix for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:16
% EndTime: 2019-03-09 06:14:26
% DurationCPUTime: 4.83s
% Computational Cost: add. (12254->439), mult. (11443->622), div. (0->0), fcn. (12288->10), ass. (0->221)
t214 = pkin(10) + qJ(3);
t208 = cos(t214);
t217 = qJ(4) + qJ(5);
t210 = cos(t217);
t224 = cos(qJ(1));
t268 = t224 * t210;
t209 = sin(t217);
t222 = sin(qJ(1));
t273 = t222 * t209;
t170 = -t208 * t273 - t268;
t269 = t224 * t209;
t272 = t222 * t210;
t171 = t208 * t272 - t269;
t207 = sin(t214);
t278 = t207 * t222;
t101 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t278;
t103 = Icges(6,5) * t171 + Icges(6,6) * t170 + Icges(6,3) * t278;
t323 = t101 + t103;
t172 = -t208 * t269 + t272;
t173 = t208 * t268 + t273;
t277 = t207 * t224;
t102 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t277;
t104 = Icges(6,5) * t173 + Icges(6,6) * t172 + Icges(6,3) * t277;
t322 = t102 + t104;
t105 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t278;
t107 = Icges(6,4) * t171 + Icges(6,2) * t170 + Icges(6,6) * t278;
t321 = t105 + t107;
t106 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t277;
t108 = Icges(6,4) * t173 + Icges(6,2) * t172 + Icges(6,6) * t277;
t320 = t106 + t108;
t109 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t278;
t111 = Icges(6,1) * t171 + Icges(6,4) * t170 + Icges(6,5) * t278;
t319 = t109 + t111;
t110 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t277;
t112 = Icges(6,1) * t173 + Icges(6,4) * t172 + Icges(6,5) * t277;
t318 = t110 + t112;
t317 = t321 * t170 + t319 * t171 + t323 * t278;
t316 = t320 * t170 + t318 * t171 + t322 * t278;
t315 = t321 * t172 + t319 * t173 + t323 * t277;
t314 = t320 * t172 + t318 * t173 + t322 * t277;
t144 = -Icges(7,3) * t208 + (Icges(7,5) * t210 - Icges(7,6) * t209) * t207;
t146 = -Icges(7,6) * t208 + (Icges(7,4) * t210 - Icges(7,2) * t209) * t207;
t148 = -Icges(7,5) * t208 + (Icges(7,1) * t210 - Icges(7,4) * t209) * t207;
t67 = t144 * t278 + t170 * t146 + t171 * t148;
t145 = -Icges(6,3) * t208 + (Icges(6,5) * t210 - Icges(6,6) * t209) * t207;
t147 = -Icges(6,6) * t208 + (Icges(6,4) * t210 - Icges(6,2) * t209) * t207;
t149 = -Icges(6,5) * t208 + (Icges(6,1) * t210 - Icges(6,4) * t209) * t207;
t68 = t145 * t278 + t170 * t147 + t171 * t149;
t313 = -t68 - t67;
t69 = t144 * t277 + t172 * t146 + t173 * t148;
t70 = t145 * t277 + t172 * t147 + t173 * t149;
t312 = -t69 - t70;
t307 = (t148 + t149) * t207 * t210;
t309 = (-t146 - t147) * t209;
t310 = -t144 - t145;
t299 = t207 * t309 + t208 * t310 + t307;
t311 = t299 * t208;
t308 = Icges(4,5) * t207;
t306 = t308 / 0.2e1;
t305 = t313 * t208 + (t222 * t317 + t316 * t224) * t207;
t304 = t312 * t208 + (t222 * t315 + t224 * t314) * t207;
t303 = t316 * t222 - t224 * t317;
t302 = t222 * t314 - t224 * t315;
t52 = -t208 * t101 + (-t105 * t209 + t109 * t210) * t207;
t54 = -t208 * t103 + (-t107 * t209 + t111 * t210) * t207;
t301 = -t52 - t54;
t53 = -t208 * t102 + (-t106 * t209 + t110 * t210) * t207;
t55 = -t208 * t104 + (-t108 * t209 + t112 * t210) * t207;
t300 = t53 + t55;
t221 = sin(qJ(4));
t190 = t221 * pkin(4) + pkin(5) * t209;
t225 = -pkin(9) - pkin(8);
t213 = -qJ(6) + t225;
t237 = -t171 * rSges(7,1) - t170 * rSges(7,2);
t266 = t224 * t221;
t276 = t207 * t225;
t255 = pkin(4) * t266 + t222 * t276;
t223 = cos(qJ(4));
t205 = t223 * pkin(4) + pkin(3);
t189 = pkin(5) * t210 + t205;
t256 = t189 - t205;
t298 = rSges(7,3) * t278 - t237 - t224 * t190 + (-t207 * t213 + t208 * t256) * t222 + t255;
t253 = t213 - t225;
t297 = (t253 - rSges(7,3)) * t208 + (rSges(7,1) * t210 - rSges(7,2) * t209 + t256) * t207;
t275 = t208 * t224;
t296 = t173 * rSges(7,1) + t172 * rSges(7,2) + rSges(7,3) * t277 + t189 * t275 + t222 * t190;
t215 = t222 ^ 2;
t216 = t224 ^ 2;
t295 = -t208 / 0.2e1;
t294 = t222 / 0.2e1;
t293 = -t224 / 0.2e1;
t292 = t224 / 0.2e1;
t186 = t207 * rSges(4,1) + t208 * rSges(4,2);
t291 = m(4) * t186;
t290 = pkin(3) * t208;
t289 = pkin(8) * t207;
t288 = -pkin(3) + t205;
t287 = t311 + (t222 * t301 - t224 * t300) * t207;
t285 = t298 * t277;
t284 = rSges(3,3) + qJ(2);
t271 = t222 * t221;
t257 = -pkin(4) * t271 - t205 * t275;
t283 = -t253 * t277 + t257 + t296;
t238 = -t171 * rSges(6,1) - t170 * rSges(6,2);
t114 = rSges(6,3) * t278 - t238;
t151 = -t208 * rSges(6,3) + (rSges(6,1) * t210 - rSges(6,2) * t209) * t207;
t84 = t208 * t114 + t151 * t278;
t281 = Icges(4,4) * t208;
t155 = -Icges(5,6) * t208 + (Icges(5,4) * t223 - Icges(5,2) * t221) * t207;
t274 = t221 * t155;
t270 = t222 * t223;
t220 = -pkin(7) - qJ(2);
t267 = t224 * t220;
t265 = t224 * t223;
t116 = t173 * rSges(6,1) + t172 * rSges(6,2) + rSges(6,3) * t277;
t229 = -t224 * t276 - t257;
t254 = pkin(3) * t275 + pkin(8) * t277;
t123 = t229 - t254;
t264 = -t116 - t123;
t122 = (t208 * t288 - t289) * t222 - t255;
t143 = (pkin(8) + t225) * t208 + t288 * t207;
t263 = t208 * t122 + t143 * t278;
t261 = -t143 - t151;
t157 = -t208 * rSges(5,3) + (rSges(5,1) * t223 - rSges(5,2) * t221) * t207;
t188 = t207 * pkin(3) - t208 * pkin(8);
t260 = -t157 - t188;
t259 = t215 * (t289 + t290) + t224 * t254;
t252 = t215 + t216;
t177 = -t208 * t271 - t265;
t178 = t208 * t270 - t266;
t124 = Icges(5,5) * t178 + Icges(5,6) * t177 + Icges(5,3) * t278;
t126 = Icges(5,4) * t178 + Icges(5,2) * t177 + Icges(5,6) * t278;
t128 = Icges(5,1) * t178 + Icges(5,4) * t177 + Icges(5,5) * t278;
t60 = -t208 * t124 + (-t126 * t221 + t128 * t223) * t207;
t154 = -Icges(5,3) * t208 + (Icges(5,5) * t223 - Icges(5,6) * t221) * t207;
t156 = -Icges(5,5) * t208 + (Icges(5,1) * t223 - Icges(5,4) * t221) * t207;
t74 = t154 * t278 + t177 * t155 + t178 * t156;
t251 = t74 / 0.2e1 + t60 / 0.2e1;
t179 = -t208 * t266 + t270;
t180 = t208 * t265 + t271;
t125 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t277;
t127 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t277;
t129 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t277;
t61 = -t208 * t125 + (-t127 * t221 + t129 * t223) * t207;
t75 = t154 * t277 + t179 * t155 + t180 * t156;
t250 = t75 / 0.2e1 + t61 / 0.2e1;
t249 = -t123 - t283;
t248 = t304 * t277 + t278 * t305;
t247 = -t143 - t297;
t246 = -t188 + t261;
t131 = t180 * rSges(5,1) + t179 * rSges(5,2) + rSges(5,3) * t277;
t245 = t278 / 0.2e1;
t244 = t277 / 0.2e1;
t219 = cos(pkin(10));
t204 = t219 * pkin(2) + pkin(1);
t243 = t224 * t204 - t222 * t220;
t36 = t208 * t298 + t278 * t297;
t242 = t222 * t122 + t224 * t123 + t259;
t241 = -t188 + t247;
t240 = rSges(4,1) * t208 - rSges(4,2) * t207;
t239 = -t178 * rSges(5,1) - t177 * rSges(5,2);
t235 = -Icges(4,2) * t207 + t281;
t234 = Icges(4,5) * t208 - Icges(4,6) * t207;
t231 = rSges(4,1) * t275 - rSges(4,2) * t277 + t222 * rSges(4,3);
t218 = sin(pkin(10));
t230 = rSges(3,1) * t219 - rSges(3,2) * t218 + pkin(1);
t228 = t208 * t287 + t248;
t227 = (-t301 - t313) * t245 + (t300 - t312) * t244;
t226 = (t222 * t300 + t224 * t301) * t295 + t304 * t294 + t305 * t293 + t303 * t245 + t302 * t244;
t193 = t224 * rSges(2,1) - t222 * rSges(2,2);
t192 = -t222 * rSges(2,1) - t224 * rSges(2,2);
t182 = Icges(4,6) * t208 + t308;
t159 = Icges(4,3) * t222 + t224 * t234;
t158 = -Icges(4,3) * t224 + t222 * t234;
t153 = t222 * t284 + t224 * t230;
t152 = -t222 * t230 + t224 * t284;
t142 = t207 * t223 * t156;
t139 = t231 + t243;
t138 = (rSges(4,3) - t220) * t224 + (-t204 - t240) * t222;
t134 = t260 * t224;
t133 = t260 * t222;
t130 = rSges(5,3) * t278 - t239;
t119 = t224 * t231 + (-t224 * rSges(4,3) + t222 * t240) * t222;
t100 = t122 * t277;
t97 = t114 * t277;
t93 = t243 + t131 + t254;
t92 = -t267 + (-t290 - t204 + (-rSges(5,3) - pkin(8)) * t207) * t222 + t239;
t89 = t246 * t224;
t88 = t246 * t222;
t87 = -t208 * t131 - t157 * t277;
t86 = t208 * t130 + t157 * t278;
t85 = -t208 * t116 - t151 * t277;
t83 = -t208 * t154 - t207 * t274 + t142;
t82 = t229 + t243 + t116;
t81 = -t267 + (-rSges(6,3) * t207 - t205 * t208 - t204) * t222 + t238 + t255;
t80 = (t130 * t224 - t131 * t222) * t207;
t77 = -t213 * t277 + t243 + t296;
t76 = (t190 - t220) * t224 + (-t189 * t208 - t204 + (-rSges(7,3) + t213) * t207) * t222 + t237;
t73 = -t116 * t278 + t97;
t72 = t241 * t224;
t71 = t241 * t222;
t66 = t222 * t130 + t224 * t131 + t259;
t59 = t125 * t277 + t179 * t127 + t180 * t129;
t58 = t124 * t277 + t179 * t126 + t180 * t128;
t57 = t125 * t278 + t177 * t127 + t178 * t129;
t56 = t124 * t278 + t177 * t126 + t178 * t128;
t51 = t208 * t264 + t261 * t277;
t50 = t263 + t84;
t37 = -t208 * t283 - t277 * t297;
t35 = t264 * t278 + t100 + t97;
t34 = t222 * t114 + t224 * t116 + t242;
t33 = -t278 * t283 + t285;
t32 = t208 * t249 + t247 * t277;
t31 = t36 + t263;
t30 = t59 * t222 - t58 * t224;
t29 = t57 * t222 - t56 * t224;
t26 = t249 * t278 + t100 + t285;
t25 = t222 * t298 + t283 * t224 + t242;
t16 = -t75 * t208 + (t222 * t58 + t224 * t59) * t207;
t15 = -t74 * t208 + (t222 * t56 + t224 * t57) * t207;
t1 = [Icges(3,2) * t219 ^ 2 + Icges(2,3) + t142 + (Icges(3,1) * t218 + 0.2e1 * Icges(3,4) * t219) * t218 + (Icges(4,4) * t207 + Icges(4,2) * t208 - t154 + t310) * t208 + (Icges(4,1) * t207 - t274 + t281 + t309) * t207 + m(7) * (t76 ^ 2 + t77 ^ 2) + m(6) * (t81 ^ 2 + t82 ^ 2) + m(5) * (t92 ^ 2 + t93 ^ 2) + m(4) * (t138 ^ 2 + t139 ^ 2) + m(3) * (t152 ^ 2 + t153 ^ 2) + m(2) * (t192 ^ 2 + t193 ^ 2) + t307; m(7) * (t222 * t76 - t224 * t77) + m(6) * (t222 * t81 - t224 * t82) + m(5) * (t222 * t92 - t224 * t93) + m(4) * (t222 * t138 - t224 * t139) + m(3) * (t222 * t152 - t224 * t153); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t252; m(7) * (t71 * t77 + t72 * t76) + m(6) * (t81 * t89 + t82 * t88) + m(5) * (t133 * t93 + t134 * t92) + (-t68 / 0.2e1 - t67 / 0.2e1 - t54 / 0.2e1 - t52 / 0.2e1 + (-Icges(4,6) * t224 + t222 * t235) * t295 + t224 * t306 - t138 * t291 + t182 * t292 - t251) * t224 + (t70 / 0.2e1 + t69 / 0.2e1 + t55 / 0.2e1 + t53 / 0.2e1 + t208 * (Icges(4,6) * t222 + t224 * t235) / 0.2e1 + t222 * t306 - t139 * t291 + t182 * t294 + t250) * t222; m(5) * (-t133 * t224 + t134 * t222) + m(6) * (t89 * t222 - t88 * t224) + m(7) * (t72 * t222 - t71 * t224); m(7) * (t25 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t34 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(5) * (t133 ^ 2 + t134 ^ 2 + t66 ^ 2) + m(4) * (t186 ^ 2 * t252 + t119 ^ 2) + (-t216 * t158 - t29 - t303) * t224 + (t215 * t159 + t30 + (-t222 * t158 + t224 * t159) * t224 + t302) * t222; (-t83 - t299) * t208 + m(7) * (t31 * t76 + t32 * t77) + m(6) * (t50 * t81 + t51 * t82) + m(5) * (t86 * t92 + t87 * t93) + (t222 * t251 + t224 * t250) * t207 + t227; m(5) * (t86 * t222 - t87 * t224) + m(6) * (t50 * t222 - t51 * t224) + m(7) * (t31 * t222 - t32 * t224); t226 + t15 * t293 + (t61 * t222 - t60 * t224) * t295 + t16 * t294 + m(6) * (t35 * t34 + t50 * t89 + t51 * t88) + m(5) * (t133 * t87 + t134 * t86 + t80 * t66) + m(7) * (t26 * t25 + t31 * t72 + t32 * t71) + (t29 * t294 + t292 * t30) * t207; (t83 * t208 + t287) * t208 + (t224 * t16 + t222 * t15 - t208 * (t222 * t60 + t224 * t61)) * t207 + m(7) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t35 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t80 ^ 2 + t86 ^ 2 + t87 ^ 2) + t248; -t311 + m(7) * (t36 * t76 + t37 * t77) + m(6) * (t81 * t84 + t82 * t85) + t227; m(6) * (t84 * t222 - t85 * t224) + m(7) * (t36 * t222 - t37 * t224); m(7) * (t33 * t25 + t36 * t72 + t37 * t71) + m(6) * (t34 * t73 + t84 * t89 + t85 * t88) + t226; m(7) * (t26 * t33 + t31 * t36 + t32 * t37) + m(6) * (t35 * t73 + t50 * t84 + t51 * t85) + t228; m(7) * (t33 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t73 ^ 2 + t84 ^ 2 + t85 ^ 2) + t228; m(7) * (t222 * t77 + t224 * t76) * t207; 0; m(7) * (-t208 * t25 + (t222 * t71 + t224 * t72) * t207); m(7) * (-t208 * t26 + (t222 * t32 + t224 * t31) * t207); m(7) * (-t208 * t33 + (t222 * t37 + t224 * t36) * t207); m(7) * (t207 ^ 2 * t252 + t208 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
