% Calculate joint inertia matrix for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:46
% EndTime: 2019-03-09 11:11:55
% DurationCPUTime: 4.30s
% Computational Cost: add. (8103->498), mult. (10310->692), div. (0->0), fcn. (10939->10), ass. (0->247)
t329 = Icges(4,1) + Icges(3,3);
t227 = sin(qJ(2));
t230 = cos(qJ(2));
t328 = (-Icges(4,4) + Icges(3,5)) * t230 + (Icges(4,5) - Icges(3,6)) * t227;
t228 = sin(qJ(1));
t327 = -t228 / 0.2e1;
t310 = t228 / 0.2e1;
t231 = cos(qJ(1));
t309 = -t231 / 0.2e1;
t326 = t231 / 0.2e1;
t220 = qJ(4) + pkin(10);
t208 = sin(t220);
t209 = cos(t220);
t129 = Icges(6,3) * t227 + (-Icges(6,5) * t208 - Icges(6,6) * t209) * t230;
t226 = sin(qJ(4));
t229 = cos(qJ(4));
t148 = Icges(5,3) * t227 + (-Icges(5,5) * t226 - Icges(5,6) * t229) * t230;
t325 = (t129 + t148) * t227;
t130 = Icges(6,6) * t227 + (-Icges(6,4) * t208 - Icges(6,2) * t209) * t230;
t131 = Icges(6,5) * t227 + (-Icges(6,1) * t208 - Icges(6,4) * t209) * t230;
t151 = Icges(5,6) * t227 + (-Icges(5,4) * t226 - Icges(5,2) * t229) * t230;
t154 = Icges(5,5) * t227 + (-Icges(5,1) * t226 - Icges(5,4) * t229) * t230;
t324 = -t130 * t209 - t131 * t208 - t151 * t229 - t154 * t226;
t311 = t227 / 0.2e1;
t293 = t228 * t208;
t296 = t227 * t231;
t144 = t209 * t296 - t293;
t292 = t228 * t209;
t145 = t208 * t296 + t292;
t287 = t230 * t231;
t90 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t287;
t92 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t287;
t94 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t287;
t31 = t144 * t92 + t145 * t94 + t90 * t287;
t146 = t208 * t231 + t227 * t292;
t147 = -t209 * t231 + t227 * t293;
t289 = t228 * t230;
t91 = Icges(6,5) * t147 + Icges(6,6) * t146 + Icges(6,3) * t289;
t93 = Icges(6,4) * t147 + Icges(6,2) * t146 + Icges(6,6) * t289;
t95 = Icges(6,1) * t147 + Icges(6,4) * t146 + Icges(6,5) * t289;
t32 = t144 * t93 + t145 * t95 + t91 * t287;
t288 = t229 * t231;
t291 = t228 * t226;
t173 = t227 * t288 - t291;
t274 = t226 * t296;
t290 = t228 * t229;
t174 = t274 + t290;
t105 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t287;
t107 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t287;
t109 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t287;
t41 = t105 * t287 + t173 * t107 + t174 * t109;
t175 = t226 * t231 + t227 * t290;
t176 = t227 * t291 - t288;
t106 = Icges(5,5) * t176 + Icges(5,6) * t175 + Icges(5,3) * t289;
t108 = Icges(5,4) * t176 + Icges(5,2) * t175 + Icges(5,6) * t289;
t110 = Icges(5,1) * t176 + Icges(5,4) * t175 + Icges(5,5) * t289;
t42 = t106 * t287 + t173 * t108 + t174 * t110;
t57 = t129 * t287 + t144 * t130 + t145 * t131;
t64 = t148 * t287 + t173 * t151 + t174 * t154;
t323 = ((t31 + t41) * t231 + (t32 + t42) * t228) * t230 + (t64 + t57) * t227;
t33 = t146 * t92 + t147 * t94 + t90 * t289;
t34 = t146 * t93 + t147 * t95 + t91 * t289;
t43 = t105 * t289 + t107 * t175 + t109 * t176;
t44 = t106 * t289 + t108 * t175 + t110 * t176;
t58 = t129 * t289 + t130 * t146 + t131 * t147;
t65 = t148 * t289 + t151 * t175 + t154 * t176;
t322 = ((t33 + t43) * t231 + (t34 + t44) * t228) * t230 + (t65 + t58) * t227;
t39 = t227 * t90 + (-t208 * t94 - t209 * t92) * t230;
t48 = t105 * t227 + (-t107 * t229 - t109 * t226) * t230;
t321 = t39 + t48;
t40 = t227 * t91 + (-t208 * t95 - t209 * t93) * t230;
t49 = t106 * t227 + (-t108 * t229 - t110 * t226) * t230;
t320 = t40 + t49;
t308 = pkin(4) * t226;
t180 = pkin(5) * t208 + t308;
t270 = -t180 + t308;
t225 = -qJ(5) - pkin(8);
t219 = -pkin(9) + t225;
t278 = -t219 + t225;
t117 = t278 * t227 + t270 * t230;
t210 = qJ(6) + t220;
t205 = sin(t210);
t206 = cos(t210);
t128 = rSges(7,3) * t227 + (-rSges(7,1) * t205 - rSges(7,2) * t206) * t230;
t319 = -t117 - t128;
t318 = -t328 * t228 + t329 * t231;
t317 = t329 * t228 + t328 * t231;
t222 = t228 ^ 2;
t224 = t231 ^ 2;
t316 = 0.2e1 * t230;
t315 = m(4) / 0.2e1;
t314 = m(5) / 0.2e1;
t313 = m(6) / 0.2e1;
t312 = m(7) / 0.2e1;
t207 = t229 * pkin(4) + pkin(3);
t307 = (t230 * t324 + t325) * t227;
t283 = -pkin(4) * t274 - t228 * t207;
t179 = pkin(5) * t209 + t207;
t285 = t228 * t179 + t180 * t296;
t76 = t278 * t287 + t283 + t285;
t295 = t228 * t205;
t138 = t206 * t296 - t295;
t294 = t228 * t206;
t139 = t205 * t296 + t294;
t88 = t139 * rSges(7,1) + t138 * rSges(7,2) + rSges(7,3) * t287;
t306 = t76 + t88;
t282 = t231 * t207 + t225 * t289;
t297 = t179 * t231;
t77 = -t297 + (-t219 * t230 - t270 * t227) * t228 + t282;
t140 = t205 * t231 + t227 * t294;
t141 = -t206 * t231 + t227 * t295;
t259 = -t141 * rSges(7,1) - t140 * rSges(7,2);
t89 = rSges(7,3) * t289 - t259;
t305 = t77 + t89;
t304 = t231 * rSges(4,1);
t303 = t231 * rSges(3,3);
t302 = Icges(3,4) * t227;
t301 = Icges(3,4) * t230;
t300 = Icges(4,6) * t227;
t299 = Icges(4,6) * t230;
t298 = qJ(3) * t227;
t281 = pkin(2) * t287 + qJ(3) * t296;
t286 = t222 * (pkin(2) * t230 + t298) + t231 * t281;
t188 = pkin(2) * t227 - qJ(3) * t230;
t284 = rSges(4,2) * t227 + rSges(4,3) * t230 - t188;
t280 = t231 * pkin(1) + t228 * pkin(7);
t279 = t228 * pkin(3) + pkin(8) * t287;
t277 = t222 + t224;
t82 = Icges(7,5) * t139 + Icges(7,6) * t138 + Icges(7,3) * t287;
t84 = Icges(7,4) * t139 + Icges(7,2) * t138 + Icges(7,6) * t287;
t86 = Icges(7,1) * t139 + Icges(7,4) * t138 + Icges(7,5) * t287;
t37 = t227 * t82 + (-t205 * t86 - t206 * t84) * t230;
t83 = Icges(7,5) * t141 + Icges(7,6) * t140 + Icges(7,3) * t289;
t85 = Icges(7,4) * t141 + Icges(7,2) * t140 + Icges(7,6) * t289;
t87 = Icges(7,1) * t141 + Icges(7,4) * t140 + Icges(7,5) * t289;
t38 = t227 * t83 + (-t205 * t87 - t206 * t85) * t230;
t26 = t138 * t84 + t139 * t86 + t82 * t287;
t27 = t138 * t85 + t139 * t87 + t83 * t287;
t123 = Icges(7,3) * t227 + (-Icges(7,5) * t205 - Icges(7,6) * t206) * t230;
t124 = Icges(7,6) * t227 + (-Icges(7,4) * t205 - Icges(7,2) * t206) * t230;
t125 = Icges(7,5) * t227 + (-Icges(7,1) * t205 - Icges(7,4) * t206) * t230;
t52 = t123 * t287 + t138 * t124 + t139 * t125;
t5 = t52 * t227 + (t228 * t27 + t231 * t26) * t230;
t28 = t140 * t84 + t141 * t86 + t82 * t289;
t29 = t140 * t85 + t141 * t87 + t83 * t289;
t53 = t123 * t289 + t124 * t140 + t125 * t141;
t6 = t53 * t227 + (t228 * t29 + t231 * t28) * t230;
t119 = t227 * t123;
t247 = -t124 * t206 - t125 * t205;
t61 = (t230 * t247 + t119) * t227;
t276 = t5 * t287 + t6 * t289 + t227 * (t61 + (t228 * t38 + t231 * t37) * t230);
t275 = t313 + t312;
t96 = t145 * rSges(6,1) + t144 * rSges(6,2) + rSges(6,3) * t287;
t111 = t174 * rSges(5,1) + t173 * rSges(5,2) + rSges(5,3) * t287;
t273 = t289 / 0.2e1;
t272 = t287 / 0.2e1;
t271 = -Icges(4,4) * t227 / 0.2e1 + Icges(3,5) * t311 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t230;
t269 = -pkin(8) * t227 - t188;
t217 = t231 * pkin(3);
t268 = t228 * (pkin(8) * t289 - t217) + t231 * t279 + t286;
t267 = t280 + t281;
t12 = t26 * t228 - t231 * t27;
t13 = t28 * t228 - t231 * t29;
t266 = t12 * t272 + t13 * t273 + t6 * t309 + t5 * t310 + (t37 * t228 - t38 * t231) * t311;
t265 = t61 + (t38 + t53) * t273 + (t37 + t52) * t272;
t163 = rSges(5,3) * t227 + (-rSges(5,1) * t226 - rSges(5,2) * t229) * t230;
t264 = -t163 + t269;
t166 = -t230 * t308 + (-pkin(8) - t225) * t227;
t263 = -t166 + t269;
t262 = rSges(3,1) * t230 - rSges(3,2) * t227;
t261 = -rSges(5,1) * t176 - rSges(5,2) * t175;
t260 = -rSges(6,1) * t147 - rSges(6,2) * t146;
t116 = t217 + (-pkin(8) * t230 + t227 * t308) * t228 - t282;
t118 = t128 * t289;
t135 = t166 * t289;
t23 = t117 * t289 + t118 + t135 + (-t116 - t305) * t227;
t237 = -t225 * t287 - t283;
t115 = t237 - t279;
t102 = t227 * t115;
t81 = t227 * t88;
t24 = t227 * t76 + t102 + t81 + (-t166 + t319) * t287;
t258 = t228 * t24 + t23 * t231;
t136 = rSges(6,3) * t227 + (-rSges(6,1) * t208 - rSges(6,2) * t209) * t230;
t97 = rSges(6,3) * t289 - t260;
t46 = t136 * t289 + t135 + (-t116 - t97) * t227;
t47 = t227 * t96 + t102 + (-t136 - t166) * t287;
t257 = t228 * t47 + t231 * t46;
t233 = t263 + t319;
t59 = t233 * t228;
t60 = t233 * t231;
t256 = t228 * t59 + t231 * t60;
t68 = -t227 * t89 + t118;
t69 = -t128 * t287 + t81;
t255 = t228 * t69 + t231 * t68;
t240 = -t136 + t263;
t78 = t240 * t228;
t79 = t240 * t231;
t254 = t228 * t78 + t231 * t79;
t253 = Icges(3,1) * t230 - t302;
t252 = -Icges(3,2) * t227 + t301;
t249 = -Icges(4,2) * t230 + t300;
t248 = Icges(4,3) * t227 - t299;
t239 = rSges(3,1) * t287 - rSges(3,2) * t296 + t228 * rSges(3,3);
t238 = t228 * rSges(4,1) - rSges(4,2) * t287 + rSges(4,3) * t296;
t236 = t48 / 0.2e1 + t39 / 0.2e1 + t64 / 0.2e1 + t57 / 0.2e1;
t235 = t65 / 0.2e1 + t58 / 0.2e1 + t49 / 0.2e1 + t40 / 0.2e1;
t234 = t231 * t115 + t228 * t116 + t268;
t216 = t231 * pkin(7);
t55 = t297 + t216 + (-pkin(1) + (-qJ(3) - t180) * t227 + (-rSges(7,3) - pkin(2) + t219) * t230) * t228 + t259;
t56 = -t219 * t287 + t267 + t285 + t88;
t62 = t216 + (-pkin(1) + (-rSges(6,3) - pkin(2)) * t230 + (-qJ(3) - t308) * t227) * t228 + t260 + t282;
t63 = t237 + t267 + t96;
t232 = (t228 * t56 + t231 * t55) * t312 + (t228 * t63 + t231 * t62) * t313;
t223 = t230 ^ 2;
t221 = t227 ^ 2;
t192 = rSges(2,1) * t231 - t228 * rSges(2,2);
t191 = -t228 * rSges(2,1) - rSges(2,2) * t231;
t190 = rSges(3,1) * t227 + rSges(3,2) * t230;
t127 = t284 * t231;
t126 = t284 * t228;
t121 = t239 + t280;
t120 = t303 + t216 + (-pkin(1) - t262) * t228;
t114 = t238 + t267;
t113 = t304 + t216 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t230 + (-rSges(4,3) - qJ(3)) * t227) * t228;
t112 = rSges(5,3) * t289 - t261;
t101 = t264 * t231;
t100 = t264 * t228;
t99 = t116 * t287;
t98 = t231 * t239 + (t262 * t228 - t303) * t228;
t80 = t89 * t287;
t75 = t227 * t111 - t163 * t287;
t74 = -t112 * t227 + t163 * t289;
t73 = t267 + t111 + t279;
t72 = t216 + t217 + (-t298 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(8)) * t230) * t228 + t261;
t70 = t231 * t238 + (-t304 + (-rSges(4,2) * t230 + rSges(4,3) * t227) * t228) * t228 + t286;
t67 = (-t111 * t228 + t112 * t231) * t230;
t54 = -t88 * t289 + t80;
t45 = t111 * t231 + t228 * t112 + t268;
t30 = t99 + (t231 * t97 + (-t115 - t96) * t228) * t230;
t25 = t228 * t97 + t231 * t96 + t234;
t22 = t43 * t228 - t231 * t44;
t21 = t41 * t228 - t231 * t42;
t20 = t80 + t99 + (t231 * t77 + (-t115 - t306) * t228) * t230;
t18 = t33 * t228 - t231 * t34;
t17 = t31 * t228 - t231 * t32;
t16 = t305 * t228 + t306 * t231 + t234;
t1 = [Icges(2,3) + t119 + m(7) * (t55 ^ 2 + t56 ^ 2) + m(6) * (t62 ^ 2 + t63 ^ 2) + m(5) * (t72 ^ 2 + t73 ^ 2) + m(4) * (t113 ^ 2 + t114 ^ 2) + m(3) * (t120 ^ 2 + t121 ^ 2) + m(2) * (t191 ^ 2 + t192 ^ 2) + (t247 + t300 + t302 + (Icges(3,2) + Icges(4,3)) * t230 + t324) * t230 + (t299 + t301 + (Icges(3,1) + Icges(4,2)) * t227) * t227 + t325; (-t38 / 0.2e1 - t53 / 0.2e1 + t271 * t231 + (Icges(4,5) * t309 + Icges(3,6) * t326 + t248 * t310 + t252 * t327) * t230 + (Icges(4,4) * t309 + Icges(3,5) * t326 + t249 * t310 + t253 * t327) * t227 - t235) * t231 + (t37 / 0.2e1 + t52 / 0.2e1 + t271 * t228 + (Icges(4,5) * t327 + Icges(3,6) * t310 + t248 * t309 + t252 * t326) * t230 + (Icges(4,4) * t327 + Icges(3,5) * t310 + t249 * t309 + t253 * t326) * t227 + t236) * t228 + m(7) * (t55 * t60 + t56 * t59) + m(6) * (t62 * t79 + t63 * t78) + m(5) * (t100 * t73 + t101 * t72) + m(4) * (t113 * t127 + t114 * t126) + m(3) * (-t120 * t231 - t121 * t228) * t190; m(7) * (t16 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(6) * (t25 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2 + t45 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2 + t70 ^ 2) + m(3) * (t277 * t190 ^ 2 + t98 ^ 2) + (t318 * t224 - t13 - t18 - t22) * t231 + (t12 + t17 + t21 + t317 * t222 + (t318 * t228 + t317 * t231) * t231) * t228; 0.2e1 * ((t228 * t73 + t231 * t72) * t314 + (t113 * t231 + t114 * t228) * t315 + t232) * t227; m(7) * (-t230 * t16 + t256 * t227) + m(6) * (t254 * t227 - t230 * t25) + m(5) * (-t230 * t45 + (t100 * t228 + t101 * t231) * t227) + m(4) * (-t230 * t70 + (t126 * t228 + t127 * t231) * t227); 0.2e1 * (t315 + t314 + t275) * (t277 * t221 + t223); m(7) * (t23 * t55 + t24 * t56) + m(6) * (t46 * t62 + t47 * t63) + m(5) * (t72 * t74 + t73 * t75) + (t228 * t235 + t231 * t236) * t230 + t265 + t307; m(7) * (t16 * t20 + t23 * t60 + t24 * t59) + m(6) * (t30 * t25 + t46 * t79 + t47 * t78) + m(5) * (t100 * t75 + t101 * t74 + t67 * t45) + ((t21 / 0.2e1 + t17 / 0.2e1) * t231 + (t18 / 0.2e1 + t22 / 0.2e1) * t228) * t230 + t266 + (t321 * t228 - t320 * t231) * t311 + t323 * t310 + t322 * t309; m(5) * (-t67 * t230 + (t228 * t75 + t231 * t74) * t227) + m(6) * (t257 * t227 - t30 * t230) + m(7) * (-t20 * t230 + t258 * t227); t307 * t227 + m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t30 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t67 ^ 2 + t74 ^ 2 + t75 ^ 2) + (t323 * t231 + t322 * t228 + (t320 * t228 + t321 * t231) * t227) * t230 + t276; t232 * t316; m(7) * (t227 * t16 + t256 * t230) + m(6) * (t227 * t25 + t254 * t230); t275 * (-0.1e1 + t277) * t227 * t316; m(7) * (t227 * t20 + t258 * t230) + m(6) * (t227 * t30 + t257 * t230); 0.2e1 * t275 * (t223 * t277 + t221); m(7) * (t55 * t68 + t56 * t69) + t265; m(7) * (t54 * t16 + t59 * t69 + t60 * t68) + t266; m(7) * (t227 * t255 - t54 * t230); m(7) * (t54 * t20 + t23 * t68 + t24 * t69) + t276; m(7) * (t54 * t227 + t230 * t255); m(7) * (t54 ^ 2 + t68 ^ 2 + t69 ^ 2) + t276;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
