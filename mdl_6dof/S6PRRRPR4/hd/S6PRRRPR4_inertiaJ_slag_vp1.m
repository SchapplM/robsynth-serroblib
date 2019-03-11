% Calculate joint inertia matrix for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:42
% EndTime: 2019-03-08 23:15:54
% DurationCPUTime: 5.66s
% Computational Cost: add. (30599->620), mult. (58081->875), div. (0->0), fcn. (74215->14), ass. (0->290)
t331 = m(6) + m(7);
t260 = sin(pkin(11));
t262 = cos(pkin(11));
t269 = cos(qJ(2));
t263 = cos(pkin(6));
t267 = sin(qJ(2));
t316 = t263 * t267;
t240 = t260 * t269 + t262 * t316;
t266 = sin(qJ(3));
t261 = sin(pkin(6));
t326 = cos(qJ(3));
t282 = t261 * t326;
t229 = t240 * t266 + t262 * t282;
t242 = -t260 * t316 + t262 * t269;
t231 = t242 * t266 - t260 * t282;
t318 = t261 * t266;
t243 = -t263 * t326 + t267 * t318;
t230 = t240 * t326 - t262 * t318;
t315 = t263 * t269;
t239 = t260 * t267 - t262 * t315;
t259 = qJ(4) + pkin(12);
t256 = qJ(6) + t259;
t251 = sin(t256);
t252 = cos(t256);
t194 = -t230 * t251 + t239 * t252;
t195 = t230 * t252 + t239 * t251;
t129 = Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t229;
t131 = Icges(7,4) * t195 + Icges(7,2) * t194 + Icges(7,6) * t229;
t133 = Icges(7,1) * t195 + Icges(7,4) * t194 + Icges(7,5) * t229;
t232 = t242 * t326 + t260 * t318;
t241 = t260 * t315 + t262 * t267;
t196 = -t232 * t251 + t241 * t252;
t197 = t232 * t252 + t241 * t251;
t67 = t129 * t231 + t131 * t196 + t133 * t197;
t130 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t231;
t132 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t231;
t134 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t231;
t68 = t130 * t231 + t132 * t196 + t134 * t197;
t244 = t263 * t266 + t267 * t282;
t317 = t261 * t269;
t224 = -t244 * t251 - t252 * t317;
t225 = t244 * t252 - t251 * t317;
t163 = Icges(7,5) * t225 + Icges(7,6) * t224 + Icges(7,3) * t243;
t164 = Icges(7,4) * t225 + Icges(7,2) * t224 + Icges(7,6) * t243;
t165 = Icges(7,1) * t225 + Icges(7,4) * t224 + Icges(7,5) * t243;
t88 = t163 * t231 + t164 * t196 + t165 * t197;
t8 = t229 * t67 + t231 * t68 + t243 * t88;
t330 = t8 / 0.2e1;
t329 = t229 / 0.2e1;
t328 = t231 / 0.2e1;
t327 = t243 / 0.2e1;
t268 = cos(qJ(4));
t324 = t268 * pkin(4);
t265 = sin(qJ(4));
t322 = t239 * t265;
t321 = t241 * t265;
t320 = t260 * t261;
t319 = t261 * t262;
t254 = sin(t259);
t281 = pkin(5) * t254;
t255 = cos(t259);
t302 = pkin(5) * t255;
t118 = pkin(10) * t229 + t230 * t302 + t239 * t281;
t136 = rSges(7,1) * t195 + rSges(7,2) * t194 + rSges(7,3) * t229;
t314 = t118 + t136;
t119 = pkin(10) * t231 + t232 * t302 + t241 * t281;
t137 = rSges(7,1) * t197 + rSges(7,2) * t196 + rSges(7,3) * t231;
t313 = t119 + t137;
t138 = pkin(4) * t322 + qJ(5) * t229 + t230 * t324;
t192 = t230 * pkin(3) + t229 * pkin(9);
t183 = t241 * t192;
t312 = t241 * t138 + t183;
t139 = pkin(4) * t321 + qJ(5) * t231 + t232 * t324;
t200 = -t232 * t254 + t241 * t255;
t201 = t232 * t255 + t241 * t254;
t147 = rSges(6,1) * t201 + rSges(6,2) * t200 + rSges(6,3) * t231;
t311 = -t139 - t147;
t204 = -t232 * t265 + t241 * t268;
t205 = t232 * t268 + t321;
t157 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t231;
t193 = t232 * pkin(3) + t231 * pkin(9);
t310 = -t157 - t193;
t159 = pkin(10) * t243 + t244 * t302 - t281 * t317;
t166 = rSges(7,1) * t225 + rSges(7,2) * t224 + rSges(7,3) * t243;
t309 = t159 + t166;
t227 = -t244 * t254 - t255 * t317;
t228 = t244 * t255 - t254 * t317;
t170 = rSges(6,1) * t228 + rSges(6,2) * t227 + rSges(6,3) * t243;
t289 = t265 * t317;
t171 = -pkin(4) * t289 + qJ(5) * t243 + t244 * t324;
t308 = -t170 - t171;
t233 = -t244 * t265 - t268 * t317;
t234 = t244 * t268 - t289;
t184 = rSges(5,1) * t234 + rSges(5,2) * t233 + rSges(5,3) * t243;
t226 = t244 * pkin(3) + t243 * pkin(9);
t307 = -t184 - t226;
t306 = t192 * t317 + t239 * t226;
t223 = pkin(2) * t242 + pkin(8) * t241;
t221 = t263 * t223;
t305 = t263 * t193 + t221;
t222 = pkin(2) * t240 + pkin(8) * t239;
t304 = -t192 - t222;
t303 = t222 * t320 + t223 * t319;
t79 = t129 * t243 + t131 * t224 + t133 * t225;
t80 = t130 * t243 + t132 * t224 + t134 * t225;
t96 = t163 * t243 + t164 * t224 + t165 * t225;
t28 = t229 * t79 + t231 * t80 + t243 * t96;
t65 = t129 * t229 + t131 * t194 + t133 * t195;
t66 = t130 * t229 + t132 * t194 + t134 * t195;
t87 = t163 * t229 + t164 * t194 + t165 * t195;
t7 = t229 * t65 + t231 * t66 + t243 * t87;
t300 = t229 * t7 + t231 * t8 + t243 * t28;
t198 = -t230 * t254 + t239 * t255;
t199 = t230 * t255 + t239 * t254;
t140 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t229;
t142 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t229;
t144 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t229;
t69 = t140 * t229 + t142 * t198 + t144 * t199;
t141 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t231;
t143 = Icges(6,4) * t201 + Icges(6,2) * t200 + Icges(6,6) * t231;
t145 = Icges(6,1) * t201 + Icges(6,4) * t200 + Icges(6,5) * t231;
t70 = t141 * t229 + t143 * t198 + t145 * t199;
t167 = Icges(6,5) * t228 + Icges(6,6) * t227 + Icges(6,3) * t243;
t168 = Icges(6,4) * t228 + Icges(6,2) * t227 + Icges(6,6) * t243;
t169 = Icges(6,1) * t228 + Icges(6,4) * t227 + Icges(6,5) * t243;
t89 = t167 * t229 + t168 * t198 + t169 * t199;
t19 = t69 * t239 + t70 * t241 - t317 * t89;
t202 = -t230 * t265 + t239 * t268;
t203 = t230 * t268 + t322;
t150 = Icges(5,5) * t203 + Icges(5,6) * t202 + Icges(5,3) * t229;
t152 = Icges(5,4) * t203 + Icges(5,2) * t202 + Icges(5,6) * t229;
t154 = Icges(5,1) * t203 + Icges(5,4) * t202 + Icges(5,5) * t229;
t75 = t150 * t229 + t152 * t202 + t154 * t203;
t151 = Icges(5,5) * t205 + Icges(5,6) * t204 + Icges(5,3) * t231;
t153 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t231;
t155 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t231;
t76 = t151 * t229 + t153 * t202 + t155 * t203;
t180 = Icges(5,5) * t234 + Icges(5,6) * t233 + Icges(5,3) * t243;
t181 = Icges(5,4) * t234 + Icges(5,2) * t233 + Icges(5,6) * t243;
t182 = Icges(5,1) * t234 + Icges(5,4) * t233 + Icges(5,5) * t243;
t93 = t180 * t229 + t181 * t202 + t182 * t203;
t29 = t75 * t239 + t76 * t241 - t317 * t93;
t298 = t19 / 0.2e1 + t29 / 0.2e1;
t71 = t140 * t231 + t142 * t200 + t144 * t201;
t72 = t141 * t231 + t143 * t200 + t145 * t201;
t90 = t167 * t231 + t168 * t200 + t169 * t201;
t20 = t71 * t239 + t72 * t241 - t317 * t90;
t77 = t150 * t231 + t152 * t204 + t154 * t205;
t78 = t151 * t231 + t153 * t204 + t155 * t205;
t94 = t180 * t231 + t181 * t204 + t182 * t205;
t30 = t77 * t239 + t78 * t241 - t317 * t94;
t297 = t20 / 0.2e1 + t30 / 0.2e1;
t21 = t89 * t263 + (t260 * t70 - t262 * t69) * t261;
t33 = t93 * t263 + (t260 * t76 - t262 * t75) * t261;
t296 = t21 / 0.2e1 + t33 / 0.2e1;
t22 = t90 * t263 + (t260 * t72 - t262 * t71) * t261;
t34 = t94 * t263 + (t260 * t78 - t262 * t77) * t261;
t295 = t22 / 0.2e1 + t34 / 0.2e1;
t15 = t229 * t69 + t231 * t70 + t243 * t89;
t23 = t229 * t75 + t231 * t76 + t243 * t93;
t294 = t23 / 0.2e1 + t15 / 0.2e1;
t16 = t229 * t71 + t231 * t72 + t243 * t90;
t24 = t229 * t77 + t231 * t78 + t243 * t94;
t293 = t24 / 0.2e1 + t16 / 0.2e1;
t81 = t140 * t243 + t142 * t227 + t144 * t228;
t82 = t141 * t243 + t143 * t227 + t145 * t228;
t99 = t167 * t243 + t168 * t227 + t169 * t228;
t38 = t81 * t239 + t82 * t241 - t317 * t99;
t105 = t180 * t243 + t181 * t233 + t182 * t234;
t84 = t150 * t243 + t152 * t233 + t154 * t234;
t85 = t151 * t243 + t153 * t233 + t155 * t234;
t41 = -t105 * t317 + t84 * t239 + t85 * t241;
t292 = t38 / 0.2e1 + t41 / 0.2e1;
t36 = t229 * t81 + t231 * t82 + t243 * t99;
t40 = t105 * t243 + t229 * t84 + t231 * t85;
t291 = -t40 / 0.2e1 - t36 / 0.2e1;
t39 = t99 * t263 + (t260 * t82 - t262 * t81) * t261;
t42 = t105 * t263 + (t260 * t85 - t262 * t84) * t261;
t290 = t42 / 0.2e1 + t39 / 0.2e1;
t288 = -t139 - t313;
t287 = t263 * t139 + t305;
t286 = -t138 + t304;
t285 = -t193 + t311;
t284 = -t171 - t309;
t283 = -t226 + t308;
t218 = t244 * rSges(4,1) - t243 * rSges(4,2) - rSges(4,3) * t317;
t245 = (pkin(2) * t267 - pkin(8) * t269) * t261;
t280 = (-t218 - t245) * t261;
t279 = -t193 + t288;
t278 = t138 * t317 + t239 * t171 + t306;
t277 = -t226 + t284;
t276 = t192 * t320 + t193 * t319 + t303;
t275 = (-t245 + t307) * t261;
t13 = t65 * t239 + t66 * t241 - t317 * t87;
t14 = t67 * t239 + t68 * t241 - t317 * t88;
t32 = t79 * t239 + t80 * t241 - t317 * t96;
t274 = t14 * t328 - t28 * t317 / 0.2e1 + t239 * t7 / 0.2e1 + t32 * t327 + t241 * t330 + t13 * t329;
t17 = t87 * t263 + (t260 * t66 - t262 * t65) * t261;
t18 = t88 * t263 + (t260 * t68 - t262 * t67) * t261;
t37 = t96 * t263 + (t260 * t80 - t262 * t79) * t261;
t273 = t320 * t330 + t17 * t329 + t18 * t328 + t263 * t28 / 0.2e1 + t37 * t327 - t7 * t319 / 0.2e1;
t272 = (-t245 + t283) * t261;
t271 = t138 * t320 + t139 * t319 + t276;
t270 = (-t245 + t277) * t261;
t238 = t263 * rSges(3,3) + (rSges(3,1) * t267 + rSges(3,2) * t269) * t261;
t237 = Icges(3,5) * t263 + (Icges(3,1) * t267 + Icges(3,4) * t269) * t261;
t236 = Icges(3,6) * t263 + (Icges(3,4) * t267 + Icges(3,2) * t269) * t261;
t235 = Icges(3,3) * t263 + (Icges(3,5) * t267 + Icges(3,6) * t269) * t261;
t217 = Icges(4,1) * t244 - Icges(4,4) * t243 - Icges(4,5) * t317;
t216 = Icges(4,4) * t244 - Icges(4,2) * t243 - Icges(4,6) * t317;
t215 = Icges(4,5) * t244 - Icges(4,6) * t243 - Icges(4,3) * t317;
t214 = rSges(3,1) * t242 - rSges(3,2) * t241 + rSges(3,3) * t320;
t213 = rSges(3,1) * t240 - rSges(3,2) * t239 - rSges(3,3) * t319;
t212 = Icges(3,1) * t242 - Icges(3,4) * t241 + Icges(3,5) * t320;
t211 = Icges(3,1) * t240 - Icges(3,4) * t239 - Icges(3,5) * t319;
t210 = Icges(3,4) * t242 - Icges(3,2) * t241 + Icges(3,6) * t320;
t209 = Icges(3,4) * t240 - Icges(3,2) * t239 - Icges(3,6) * t319;
t208 = Icges(3,5) * t242 - Icges(3,6) * t241 + Icges(3,3) * t320;
t207 = Icges(3,5) * t240 - Icges(3,6) * t239 - Icges(3,3) * t319;
t187 = -t213 * t263 - t238 * t319;
t186 = t214 * t263 - t238 * t320;
t179 = rSges(4,1) * t232 - rSges(4,2) * t231 + rSges(4,3) * t241;
t178 = rSges(4,1) * t230 - rSges(4,2) * t229 + rSges(4,3) * t239;
t177 = Icges(4,1) * t232 - Icges(4,4) * t231 + Icges(4,5) * t241;
t176 = Icges(4,1) * t230 - Icges(4,4) * t229 + Icges(4,5) * t239;
t175 = Icges(4,4) * t232 - Icges(4,2) * t231 + Icges(4,6) * t241;
t174 = Icges(4,4) * t230 - Icges(4,2) * t229 + Icges(4,6) * t239;
t173 = Icges(4,5) * t232 - Icges(4,6) * t231 + Icges(4,3) * t241;
t172 = Icges(4,5) * t230 - Icges(4,6) * t229 + Icges(4,3) * t239;
t162 = (t213 * t260 + t214 * t262) * t261;
t160 = t229 * t171;
t158 = t229 * t166;
t156 = rSges(5,1) * t203 + rSges(5,2) * t202 + rSges(5,3) * t229;
t149 = -t179 * t317 - t241 * t218;
t148 = t178 * t317 + t239 * t218;
t146 = rSges(6,1) * t199 + rSges(6,2) * t198 + rSges(6,3) * t229;
t125 = t243 * t139;
t123 = -t215 * t317 - t243 * t216 + t244 * t217;
t122 = t243 * t137;
t121 = t231 * t138;
t120 = t231 * t136;
t117 = t178 * t241 - t179 * t239;
t116 = (-t178 - t222) * t263 + t262 * t280;
t115 = t263 * t179 + t260 * t280 + t221;
t114 = t215 * t241 - t216 * t231 + t217 * t232;
t113 = t215 * t239 - t216 * t229 + t217 * t230;
t112 = (t178 * t260 + t179 * t262) * t261 + t303;
t111 = t157 * t243 - t184 * t231;
t110 = -t156 * t243 + t184 * t229;
t109 = -t173 * t317 - t243 * t175 + t244 * t177;
t108 = -t172 * t317 - t243 * t174 + t244 * t176;
t107 = -t166 * t231 + t122;
t106 = -t136 * t243 + t158;
t104 = t173 * t241 - t175 * t231 + t177 * t232;
t103 = t172 * t241 - t174 * t231 + t176 * t232;
t102 = t173 * t239 - t175 * t229 + t177 * t230;
t101 = t172 * t239 - t174 * t229 + t176 * t230;
t100 = t156 * t231 - t157 * t229;
t98 = t241 * t307 + t310 * t317;
t97 = t156 * t317 + t239 * t184 + t306;
t95 = -t137 * t229 + t120;
t92 = (-t156 + t304) * t263 + t262 * t275;
t91 = t263 * t157 + t260 * t275 + t305;
t86 = t241 * t156 + t239 * t310 + t183;
t83 = (t156 * t260 + t157 * t262) * t261 + t276;
t74 = t243 * t147 + t231 * t308 + t125;
t73 = t229 * t170 + t160 + (-t138 - t146) * t243;
t64 = t241 * t283 + t285 * t317;
t63 = t146 * t317 + t239 * t170 + t278;
t62 = (-t146 + t286) * t263 + t262 * t272;
t61 = t263 * t147 + t260 * t272 + t287;
t60 = t231 * t146 + t229 * t311 + t121;
t59 = t123 * t263 + (-t108 * t262 + t109 * t260) * t261;
t58 = t108 * t239 + t109 * t241 - t123 * t317;
t57 = t241 * t146 + t239 * t285 + t312;
t56 = (t146 * t260 + t147 * t262) * t261 + t271;
t55 = t114 * t263 + (-t103 * t262 + t104 * t260) * t261;
t54 = t113 * t263 + (-t101 * t262 + t102 * t260) * t261;
t53 = t103 * t239 + t104 * t241 - t114 * t317;
t52 = t101 * t239 + t102 * t241 - t113 * t317;
t51 = t243 * t119 + t231 * t284 + t122 + t125;
t50 = t229 * t159 + t158 + t160 + (-t138 - t314) * t243;
t49 = t241 * t277 + t279 * t317;
t48 = t239 * t309 + t314 * t317 + t278;
t47 = (t286 - t314) * t263 + t262 * t270;
t46 = t260 * t270 + t263 * t313 + t287;
t45 = t231 * t118 + t229 * t288 + t120 + t121;
t44 = t239 * t279 + t241 * t314 + t312;
t43 = (t260 * t314 + t262 * t313) * t261 + t271;
t1 = [m(4) + m(5) + m(2) + m(3) + t331; m(3) * t162 + m(4) * t112 + m(5) * t83 + m(6) * t56 + m(7) * t43; m(7) * (t43 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t56 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t83 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t112 ^ 2 + t115 ^ 2 + t116 ^ 2) + m(3) * (t162 ^ 2 + t186 ^ 2 + t187 ^ 2) + (t18 + t34 + t22 + t55 + (t208 * t320 - t210 * t241 + t212 * t242) * t320) * t320 + (-t17 - t33 - t21 - t54 + (-t207 * t319 - t209 * t239 + t211 * t240) * t319 + (-t207 * t320 + t208 * t319 + t209 * t241 + t210 * t239 - t211 * t242 - t212 * t240) * t320) * t319 + ((t235 * t320 - t241 * t236 + t242 * t237) * t320 - (-t235 * t319 - t239 * t236 + t240 * t237) * t319 + t37 + t39 + t42 + t59 + ((t210 * t269 + t212 * t267) * t260 - (t209 * t269 + t211 * t267) * t262) * t261 ^ 2 + ((-t207 * t262 + t208 * t260 + t236 * t269 + t237 * t267) * t261 + t263 * t235) * t263) * t263; m(4) * t117 + m(5) * t86 + m(6) * t57 + m(7) * t44; (t32 / 0.2e1 + t58 / 0.2e1 + t292) * t263 + (t18 / 0.2e1 + t55 / 0.2e1 + t295) * t241 + (t17 / 0.2e1 + t54 / 0.2e1 + t296) * t239 + m(7) * (t43 * t44 + t46 * t49 + t47 * t48) + m(6) * (t56 * t57 + t61 * t64 + t62 * t63) + m(5) * (t83 * t86 + t91 * t98 + t92 * t97) + m(4) * (t112 * t117 + t115 * t149 + t116 * t148) + ((-t37 / 0.2e1 - t59 / 0.2e1 - t290) * t269 + (-t13 / 0.2e1 - t52 / 0.2e1 - t298) * t262 + (t14 / 0.2e1 + t53 / 0.2e1 + t297) * t260) * t261; (-t32 - t38 - t41 - t58) * t317 + (t14 + t30 + t20 + t53) * t241 + (t13 + t19 + t29 + t52) * t239 + m(7) * (t44 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t57 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t86 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(4) * (t117 ^ 2 + t148 ^ 2 + t149 ^ 2); m(5) * t100 + m(6) * t60 + m(7) * t45; -t291 * t263 + t290 * t243 + t295 * t231 + t296 * t229 + m(7) * (t43 * t45 + t46 * t51 + t47 * t50) + m(6) * (t56 * t60 + t61 * t74 + t62 * t73) + m(5) * (t100 * t83 + t110 * t92 + t111 * t91) + (t260 * t293 - t262 * t294) * t261 + t273; t291 * t317 + t292 * t243 + t293 * t241 + t294 * t239 + t297 * t231 + t298 * t229 + m(7) * (t44 * t45 + t48 * t50 + t49 * t51) + m(6) * (t57 * t60 + t63 * t73 + t64 * t74) + m(5) * (t100 * t86 + t110 * t97 + t111 * t98) + t274; (t40 + t36) * t243 + (t16 + t24) * t231 + (t15 + t23) * t229 + m(7) * (t45 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t60 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t100 ^ 2 + t110 ^ 2 + t111 ^ 2) + t300; t243 * t331; m(7) * (t229 * t46 + t231 * t47 + t243 * t43) + m(6) * (t229 * t61 + t231 * t62 + t243 * t56); m(7) * (t229 * t49 + t231 * t48 + t243 * t44) + m(6) * (t229 * t64 + t231 * t63 + t243 * t57); m(7) * (t229 * t51 + t231 * t50 + t243 * t45) + m(6) * (t229 * t74 + t231 * t73 + t243 * t60); (t229 ^ 2 + t231 ^ 2 + t243 ^ 2) * t331; m(7) * t95; m(7) * (t106 * t47 + t107 * t46 + t43 * t95) + t273; m(7) * (t106 * t48 + t107 * t49 + t44 * t95) + t274; m(7) * (t106 * t50 + t107 * t51 + t45 * t95) + t300; m(7) * (t106 * t231 + t107 * t229 + t243 * t95); m(7) * (t106 ^ 2 + t107 ^ 2 + t95 ^ 2) + t300;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
