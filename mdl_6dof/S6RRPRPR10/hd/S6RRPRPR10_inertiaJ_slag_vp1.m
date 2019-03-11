% Calculate joint inertia matrix for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:27
% EndTime: 2019-03-09 11:04:41
% DurationCPUTime: 5.92s
% Computational Cost: add. (18094->615), mult. (30478->852), div. (0->0), fcn. (37983->12), ass. (0->267)
t255 = sin(pkin(11));
t257 = cos(pkin(11));
t258 = cos(pkin(6));
t256 = sin(pkin(6));
t261 = sin(qJ(2));
t316 = t256 * t261;
t235 = -t255 * t316 + t257 * t258;
t317 = t255 * t258;
t236 = t257 * t316 + t317;
t264 = cos(qJ(2));
t314 = t256 * t264;
t177 = Icges(4,4) * t236 + Icges(4,2) * t235 - Icges(4,6) * t314;
t178 = Icges(4,1) * t236 + Icges(4,4) * t235 - Icges(4,5) * t314;
t220 = Icges(3,3) * t258 + (Icges(3,5) * t261 + Icges(3,6) * t264) * t256;
t221 = Icges(3,6) * t258 + (Icges(3,4) * t261 + Icges(3,2) * t264) * t256;
t222 = Icges(3,5) * t258 + (Icges(3,1) * t261 + Icges(3,4) * t264) * t256;
t328 = t177 * t235 + t178 * t236 + t220 * t258 + t221 * t314 + t222 * t316;
t295 = m(6) / 0.2e1 + m(7) / 0.2e1;
t327 = 0.2e1 * t295;
t176 = Icges(4,5) * t236 + Icges(4,6) * t235 - Icges(4,3) * t314;
t326 = (-t176 * t314 + t328) * t258;
t262 = sin(qJ(1));
t310 = t262 * t264;
t265 = cos(qJ(1));
t311 = t261 * t265;
t238 = t258 * t311 + t310;
t296 = pkin(11) + qJ(4);
t282 = sin(t296);
t275 = t256 * t282;
t283 = cos(t296);
t206 = t238 * t283 - t265 * t275;
t309 = t264 * t265;
t312 = t261 * t262;
t240 = -t258 * t312 + t309;
t208 = t240 * t283 + t262 * t275;
t276 = t256 * t283;
t224 = t258 * t282 + t261 * t276;
t207 = t240 * t282 - t262 * t276;
t239 = t258 * t310 + t311;
t260 = sin(qJ(6));
t263 = cos(qJ(6));
t162 = t207 * t263 - t239 * t260;
t163 = t207 * t260 + t239 * t263;
t205 = t238 * t282 + t265 * t276;
t237 = -t258 * t309 + t312;
t160 = t205 * t263 - t237 * t260;
t161 = t205 * t260 + t237 * t263;
t86 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t206;
t88 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t206;
t90 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t206;
t29 = t162 * t88 + t163 * t90 + t208 * t86;
t87 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t208;
t89 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t208;
t91 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t208;
t30 = t162 * t89 + t163 * t91 + t208 * t87;
t223 = -t258 * t283 + t261 * t275;
t203 = t223 * t263 + t260 * t314;
t204 = t223 * t260 - t263 * t314;
t107 = Icges(7,5) * t204 + Icges(7,6) * t203 + Icges(7,3) * t224;
t108 = Icges(7,4) * t204 + Icges(7,2) * t203 + Icges(7,6) * t224;
t109 = Icges(7,1) * t204 + Icges(7,4) * t203 + Icges(7,5) * t224;
t39 = t107 * t208 + t108 * t162 + t109 * t163;
t2 = t206 * t29 + t208 * t30 + t224 * t39;
t325 = t2 / 0.2e1;
t324 = t206 / 0.2e1;
t323 = t208 / 0.2e1;
t322 = t224 / 0.2e1;
t253 = pkin(3) * t257 + pkin(2);
t321 = -pkin(2) + t253;
t320 = rSges(6,3) * t205;
t277 = -rSges(7,1) * t161 - rSges(7,2) * t160;
t92 = rSges(7,3) * t206 - t277;
t319 = pkin(5) * t237 + pkin(10) * t206 + t92;
t93 = rSges(7,1) * t163 + rSges(7,2) * t162 + rSges(7,3) * t208;
t318 = pkin(5) * t239 + pkin(10) * t208 + t93;
t315 = t256 * t262;
t313 = t256 * t265;
t110 = rSges(7,1) * t204 + rSges(7,2) * t203 + rSges(7,3) * t224;
t308 = -pkin(5) * t314 + pkin(10) * t224 + t110;
t195 = pkin(2) * t240 + qJ(3) * t239;
t259 = -pkin(9) - qJ(3);
t291 = t255 * t315;
t285 = pkin(3) * t291 - t239 * t259 + t240 * t253;
t132 = -t195 + t285;
t193 = t258 * t195;
t307 = t132 * t258 + t193;
t128 = rSges(6,1) * t239 - rSges(6,2) * t208 + rSges(6,3) * t207;
t154 = pkin(4) * t208 + qJ(5) * t207;
t306 = -t128 - t154;
t227 = t237 * qJ(3);
t290 = t255 * t313;
t245 = pkin(3) * t290;
t131 = -t237 * t259 + t238 * t321 - t227 - t245;
t194 = pkin(2) * t238 + t227;
t305 = -t131 - t194;
t212 = -t238 * t255 - t257 * t313;
t213 = t238 * t257 - t290;
t141 = rSges(4,1) * t213 + rSges(4,2) * t212 + rSges(4,3) * t237;
t304 = -t141 - t194;
t196 = t205 * qJ(5);
t153 = pkin(4) * t206 + t196;
t179 = pkin(4) * t224 + qJ(5) * t223;
t303 = t153 * t314 + t179 * t237;
t166 = -Icges(6,5) * t314 - Icges(6,6) * t224 + Icges(6,3) * t223;
t167 = -Icges(6,4) * t314 - Icges(6,2) * t224 + Icges(6,6) * t223;
t302 = t166 * t223 - t167 * t224;
t170 = Icges(5,4) * t224 - Icges(5,2) * t223 - Icges(5,6) * t314;
t171 = Icges(5,1) * t224 - Icges(5,4) * t223 - Icges(5,5) * t314;
t301 = -t170 * t223 + t171 * t224;
t241 = (pkin(2) * t261 - qJ(3) * t264) * t256;
t299 = -pkin(3) * t317 - ((qJ(3) + t259) * t264 + t321 * t261) * t256 - t241;
t298 = t194 * t315 + t195 * t313;
t297 = pkin(1) * t265 + pkin(8) * t315;
t43 = t107 * t224 + t108 * t203 + t109 * t204;
t31 = t203 * t88 + t204 * t90 + t224 * t86;
t38 = t107 * t206 + t108 * t160 + t109 * t161;
t294 = t38 / 0.2e1 + t31 / 0.2e1;
t32 = t203 * t89 + t204 * t91 + t224 * t87;
t293 = t39 / 0.2e1 + t32 / 0.2e1;
t292 = -t154 - t318;
t289 = t154 * t258 + t307;
t288 = -t153 + t305;
t287 = -t179 + t299;
t130 = rSges(5,1) * t208 - rSges(5,2) * t207 + rSges(5,3) * t239;
t214 = -t240 * t255 + t257 * t315;
t215 = t240 * t257 + t291;
t142 = rSges(4,1) * t215 + rSges(4,2) * t214 + rSges(4,3) * t239;
t189 = rSges(3,1) * t240 - rSges(3,2) * t239 + rSges(3,3) * t315;
t284 = -t262 * pkin(1) + pkin(8) * t313;
t281 = t256 * (-rSges(4,1) * t236 - rSges(4,2) * t235 + rSges(4,3) * t314 - t241);
t280 = t131 * t315 + t132 * t313 + t298;
t173 = rSges(5,1) * t224 - rSges(5,2) * t223 - rSges(5,3) * t314;
t279 = t256 * (-t173 + t299);
t278 = -rSges(5,1) * t206 + rSges(5,2) * t205;
t274 = t285 + t297;
t172 = -rSges(6,1) * t314 - rSges(6,2) * t224 + rSges(6,3) * t223;
t273 = t256 * (-t172 + t287);
t272 = t153 * t315 + t154 * t313 + t280;
t271 = t256 * (t287 - t308);
t270 = -t238 * t253 + t245 + t284;
t269 = -t196 + t270;
t188 = rSges(3,1) * t238 - rSges(3,2) * t237 - rSges(3,3) * t313;
t115 = Icges(6,5) * t237 - Icges(6,6) * t206 + Icges(6,3) * t205;
t119 = Icges(6,4) * t237 - Icges(6,2) * t206 + Icges(6,6) * t205;
t123 = Icges(6,1) * t237 - Icges(6,4) * t206 + Icges(6,5) * t205;
t57 = t115 * t223 - t119 * t224 - t123 * t314;
t117 = Icges(5,5) * t206 - Icges(5,6) * t205 + Icges(5,3) * t237;
t121 = Icges(5,4) * t206 - Icges(5,2) * t205 + Icges(5,6) * t237;
t125 = Icges(5,1) * t206 - Icges(5,4) * t205 + Icges(5,5) * t237;
t59 = -t117 * t314 - t121 * t223 + t125 * t224;
t168 = -Icges(6,1) * t314 - Icges(6,4) * t224 + Icges(6,5) * t223;
t68 = t166 * t205 - t167 * t206 + t168 * t237;
t169 = Icges(5,5) * t224 - Icges(5,6) * t223 - Icges(5,3) * t314;
t70 = t169 * t237 - t170 * t205 + t171 * t206;
t268 = t70 / 0.2e1 + t68 / 0.2e1 + t59 / 0.2e1 + t57 / 0.2e1 + t294;
t116 = Icges(6,5) * t239 - Icges(6,6) * t208 + Icges(6,3) * t207;
t120 = Icges(6,4) * t239 - Icges(6,2) * t208 + Icges(6,6) * t207;
t124 = Icges(6,1) * t239 - Icges(6,4) * t208 + Icges(6,5) * t207;
t58 = t116 * t223 - t120 * t224 - t124 * t314;
t118 = Icges(5,5) * t208 - Icges(5,6) * t207 + Icges(5,3) * t239;
t122 = Icges(5,4) * t208 - Icges(5,2) * t207 + Icges(5,6) * t239;
t126 = Icges(5,1) * t208 - Icges(5,4) * t207 + Icges(5,5) * t239;
t60 = -t118 * t314 - t122 * t223 + t126 * t224;
t69 = t166 * t207 - t167 * t208 + t168 * t239;
t71 = t169 * t239 - t170 * t207 + t171 * t208;
t267 = t69 / 0.2e1 + t60 / 0.2e1 + t58 / 0.2e1 + t71 / 0.2e1 + t293;
t266 = t154 + t274;
t247 = rSges(2,1) * t265 - rSges(2,2) * t262;
t246 = -rSges(2,1) * t262 - rSges(2,2) * t265;
t225 = rSges(3,3) * t258 + (rSges(3,1) * t261 + rSges(3,2) * t264) * t256;
t187 = Icges(3,1) * t240 - Icges(3,4) * t239 + Icges(3,5) * t315;
t186 = Icges(3,1) * t238 - Icges(3,4) * t237 - Icges(3,5) * t313;
t185 = Icges(3,4) * t240 - Icges(3,2) * t239 + Icges(3,6) * t315;
t184 = Icges(3,4) * t238 - Icges(3,2) * t237 - Icges(3,6) * t313;
t183 = Icges(3,5) * t240 - Icges(3,6) * t239 + Icges(3,3) * t315;
t182 = Icges(3,5) * t238 - Icges(3,6) * t237 - Icges(3,3) * t313;
t165 = t189 + t297;
t164 = -t188 + t284;
t152 = -t188 * t258 - t225 * t313;
t151 = t189 * t258 - t225 * t315;
t140 = Icges(4,1) * t215 + Icges(4,4) * t214 + Icges(4,5) * t239;
t139 = Icges(4,1) * t213 + Icges(4,4) * t212 + Icges(4,5) * t237;
t138 = Icges(4,4) * t215 + Icges(4,2) * t214 + Icges(4,6) * t239;
t137 = Icges(4,4) * t213 + Icges(4,2) * t212 + Icges(4,6) * t237;
t136 = Icges(4,5) * t215 + Icges(4,6) * t214 + Icges(4,3) * t239;
t135 = Icges(4,5) * t213 + Icges(4,6) * t212 + Icges(4,3) * t237;
t133 = t239 * t153;
t129 = rSges(5,3) * t237 - t278;
t127 = rSges(6,1) * t237 - rSges(6,2) * t206 + t320;
t106 = (t188 * t262 + t189 * t265) * t256;
t105 = t220 * t315 - t221 * t239 + t222 * t240;
t104 = -t220 * t313 - t221 * t237 + t222 * t238;
t102 = t195 + t142 + t297;
t101 = t284 + t304;
t97 = t183 * t258 + (t185 * t264 + t187 * t261) * t256;
t96 = t182 * t258 + (t184 * t264 + t186 * t261) * t256;
t95 = t274 + t130;
t94 = (-rSges(5,3) + t259) * t237 + t270 + t278;
t85 = -t130 * t314 - t173 * t239;
t84 = t129 * t314 + t173 * t237;
t82 = t258 * t304 + t265 * t281;
t81 = t142 * t258 + t262 * t281 + t193;
t80 = -t169 * t314 + t301;
t79 = -t168 * t314 + t302;
t78 = t129 * t239 - t130 * t237;
t77 = t80 * t258;
t76 = t79 * t258;
t75 = t176 * t239 + t177 * t214 + t178 * t215;
t74 = t176 * t237 + t177 * t212 + t178 * t213;
t73 = t266 + t128;
t72 = -t320 + (-rSges(6,1) + t259) * t237 + (rSges(6,2) - pkin(4)) * t206 + t269;
t67 = (t141 * t262 + t142 * t265) * t256 + t298;
t66 = -t136 * t314 + t138 * t235 + t140 * t236;
t65 = -t135 * t314 + t137 * t235 + t139 * t236;
t64 = -t110 * t208 + t224 * t93;
t63 = t110 * t206 - t224 * t92;
t62 = t306 * t314 + (-t172 - t179) * t239;
t61 = t127 * t314 + t172 * t237 + t303;
t56 = t266 + t318;
t55 = (-pkin(5) + t259) * t237 + (-rSges(7,3) - pkin(4) - pkin(10)) * t206 + t269 + t277;
t54 = t118 * t239 - t122 * t207 + t126 * t208;
t53 = t117 * t239 - t121 * t207 + t125 * t208;
t52 = t118 * t237 - t122 * t205 + t126 * t206;
t51 = t117 * t237 - t121 * t205 + t125 * t206;
t50 = t116 * t207 - t120 * t208 + t124 * t239;
t49 = t115 * t207 - t119 * t208 + t123 * t239;
t48 = t116 * t205 - t120 * t206 + t124 * t237;
t47 = t115 * t205 - t119 * t206 + t123 * t237;
t46 = (-t129 + t305) * t258 + t265 * t279;
t45 = t130 * t258 + t262 * t279 + t307;
t44 = -t206 * t93 + t208 * t92;
t42 = t43 * t258;
t41 = t127 * t239 + t237 * t306 + t133;
t40 = t43 * t224;
t37 = (t129 * t262 + t130 * t265) * t256 + t280;
t36 = (-t127 + t288) * t258 + t265 * t273;
t35 = t128 * t258 + t262 * t273 + t289;
t34 = t292 * t314 + (-t179 - t308) * t239;
t33 = t237 * t308 + t314 * t319 + t303;
t28 = t160 * t89 + t161 * t91 + t206 * t87;
t27 = t160 * t88 + t161 * t90 + t206 * t86;
t26 = (t127 * t262 + t128 * t265) * t256 + t272;
t25 = t237 * t292 + t239 * t319 + t133;
t24 = (t288 - t319) * t258 + t265 * t271;
t23 = t258 * t318 + t262 * t271 + t289;
t22 = t77 + (t60 * t262 - t59 * t265) * t256;
t21 = t76 + (t58 * t262 - t57 * t265) * t256;
t20 = (t262 * t319 + t265 * t318) * t256 + t272;
t19 = t59 * t237 + t60 * t239 - t314 * t80;
t18 = t57 * t237 + t58 * t239 - t314 * t79;
t17 = t71 * t258 + (t262 * t54 - t265 * t53) * t256;
t16 = t70 * t258 + (t262 * t52 - t265 * t51) * t256;
t15 = t69 * t258 + (t262 * t50 - t265 * t49) * t256;
t14 = t68 * t258 + (t262 * t48 - t265 * t47) * t256;
t13 = t237 * t53 + t239 * t54 - t314 * t71;
t12 = t237 * t51 + t239 * t52 - t314 * t70;
t11 = t237 * t49 + t239 * t50 - t314 * t69;
t10 = t237 * t47 + t239 * t48 - t314 * t68;
t9 = t42 + (t32 * t262 - t31 * t265) * t256;
t8 = t31 * t237 + t32 * t239 - t314 * t43;
t7 = t31 * t206 + t32 * t208 + t40;
t6 = t39 * t258 + (t262 * t30 - t265 * t29) * t256;
t5 = t38 * t258 + (t262 * t28 - t265 * t27) * t256;
t4 = t237 * t29 + t239 * t30 - t314 * t39;
t3 = t237 * t27 + t239 * t28 - t314 * t38;
t1 = t206 * t27 + t208 * t28 + t224 * t38;
t83 = [Icges(2,3) + (-t168 - t169 - t176) * t314 + m(7) * (t55 ^ 2 + t56 ^ 2) + m(6) * (t72 ^ 2 + t73 ^ 2) + m(5) * (t94 ^ 2 + t95 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2) + m(3) * (t164 ^ 2 + t165 ^ 2) + m(2) * (t246 ^ 2 + t247 ^ 2) + t43 + t301 + t302 + t328; t42 + t77 + t76 + m(7) * (t23 * t56 + t24 * t55) + m(6) * (t35 * t73 + t36 * t72) + m(5) * (t45 * t95 + t46 * t94) + m(4) * (t101 * t82 + t102 * t81) + m(3) * (t151 * t165 + t152 * t164) + ((-t96 / 0.2e1 - t65 / 0.2e1 - t74 / 0.2e1 - t104 / 0.2e1 - t268) * t265 + (t97 / 0.2e1 + t66 / 0.2e1 + t75 / 0.2e1 + t105 / 0.2e1 + t267) * t262) * t256 + t326; m(7) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t26 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t37 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(4) * (t67 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(3) * (t106 ^ 2 + t151 ^ 2 + t152 ^ 2) + (t6 + t17 + t15 + ((t136 * t239 + t214 * t138 + t215 * t140) * t262 - (t239 * t135 + t214 * t137 + t215 * t139) * t265) * t256 + (t183 * t315 - t185 * t239 + t240 * t187) * t315) * t315 + (-t5 - t16 - t14 - ((t237 * t136 + t212 * t138 + t213 * t140) * t262 - (t135 * t237 + t137 * t212 + t139 * t213) * t265) * t256 + (-t182 * t313 - t184 * t237 + t186 * t238) * t313 + (-t182 * t315 + t183 * t313 + t239 * t184 + t237 * t185 - t240 * t186 - t238 * t187) * t315) * t313 + (t9 + t21 + t22 + (t105 + t75) * t315 + (-t104 - t74) * t313 + ((-t65 - t96) * t265 + (t66 + t97) * t262) * t256 + t326) * t258; m(7) * (t237 * t56 + t239 * t55) + m(6) * (t237 * t73 + t239 * t72) + m(5) * (t237 * t95 + t239 * t94) + m(4) * (t101 * t239 + t102 * t237); m(7) * (-t20 * t314 + t23 * t237 + t239 * t24) + m(6) * (t237 * t35 + t239 * t36 - t26 * t314) + m(5) * (t237 * t45 + t239 * t46 - t314 * t37) + m(4) * (t237 * t81 + t239 * t82 - t314 * t67); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + t295) * (t256 ^ 2 * t264 ^ 2 + t237 ^ 2 + t239 ^ 2); (-t43 - t79 - t80) * t314 + m(7) * (t33 * t55 + t34 * t56) + m(6) * (t61 * t72 + t62 * t73) + m(5) * (t84 * t94 + t85 * t95) + t267 * t239 + t268 * t237; (t8 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1) * t258 + (t6 / 0.2e1 + t17 / 0.2e1 + t15 / 0.2e1) * t239 + (t5 / 0.2e1 + t14 / 0.2e1 + t16 / 0.2e1) * t237 + m(7) * (t20 * t25 + t23 * t34 + t24 * t33) + m(6) * (t26 * t41 + t35 * t62 + t36 * t61) + m(5) * (t37 * t78 + t45 * t85 + t46 * t84) + ((-t3 / 0.2e1 - t12 / 0.2e1 - t10 / 0.2e1) * t265 + (-t9 / 0.2e1 - t21 / 0.2e1 - t22 / 0.2e1) * t264 + (t4 / 0.2e1 + t13 / 0.2e1 + t11 / 0.2e1) * t262) * t256; m(5) * (t237 * t85 + t239 * t84 - t314 * t78) + m(6) * (t237 * t62 + t239 * t61 - t314 * t41) + m(7) * (t237 * t34 + t239 * t33 - t25 * t314); (-t18 - t19 - t8) * t314 + (t4 + t13 + t11) * t239 + (t3 + t10 + t12) * t237 + m(7) * (t25 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t41 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t78 ^ 2 + t84 ^ 2 + t85 ^ 2); m(7) * (t205 * t56 + t207 * t55) + m(6) * (t205 * t73 + t207 * t72); m(7) * (t20 * t223 + t205 * t23 + t207 * t24) + m(6) * (t205 * t35 + t207 * t36 + t223 * t26); (t205 * t237 + t207 * t239 - t223 * t314) * t327; m(7) * (t205 * t34 + t207 * t33 + t223 * t25) + m(6) * (t205 * t62 + t207 * t61 + t223 * t41); (t205 ^ 2 + t207 ^ 2 + t223 ^ 2) * t327; m(7) * (t55 * t63 + t56 * t64) + t40 + t293 * t208 + t294 * t206; m(7) * (t20 * t44 + t23 * t64 + t24 * t63) + t5 * t324 + t9 * t322 + t6 * t323 + t258 * t7 / 0.2e1 + (t262 * t325 - t265 * t1 / 0.2e1) * t256; m(7) * (t237 * t64 + t239 * t63 - t314 * t44); -t7 * t314 / 0.2e1 + m(7) * (t25 * t44 + t33 * t63 + t34 * t64) + t8 * t322 + t3 * t324 + t4 * t323 + t239 * t325 + t237 * t1 / 0.2e1; m(7) * (t205 * t64 + t207 * t63 + t223 * t44); t208 * t2 + t206 * t1 + t224 * t7 + m(7) * (t44 ^ 2 + t63 ^ 2 + t64 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t83(1) t83(2) t83(4) t83(7) t83(11) t83(16); t83(2) t83(3) t83(5) t83(8) t83(12) t83(17); t83(4) t83(5) t83(6) t83(9) t83(13) t83(18); t83(7) t83(8) t83(9) t83(10) t83(14) t83(19); t83(11) t83(12) t83(13) t83(14) t83(15) t83(20); t83(16) t83(17) t83(18) t83(19) t83(20) t83(21);];
Mq  = res;
