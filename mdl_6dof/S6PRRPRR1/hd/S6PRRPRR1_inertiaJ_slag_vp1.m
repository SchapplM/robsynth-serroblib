% Calculate joint inertia matrix for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:51:06
% EndTime: 2019-03-08 21:51:19
% DurationCPUTime: 5.79s
% Computational Cost: add. (29652->561), mult. (40123->825), div. (0->0), fcn. (50014->14), ass. (0->270)
t258 = qJ(3) + pkin(12);
t285 = qJ(5) + t258;
t252 = sin(t285);
t262 = cos(pkin(6));
t278 = cos(t285);
t260 = sin(pkin(6));
t266 = sin(qJ(2));
t317 = t260 * t266;
t232 = t252 * t317 - t262 * t278;
t274 = t260 * t278;
t233 = t262 * t252 + t266 * t274;
t269 = cos(qJ(2));
t315 = t260 * t269;
t190 = Icges(6,5) * t233 - Icges(6,6) * t232 - Icges(6,3) * t315;
t191 = Icges(6,4) * t233 - Icges(6,2) * t232 - Icges(6,6) * t315;
t192 = Icges(6,1) * t233 - Icges(6,4) * t232 - Icges(6,5) * t315;
t259 = sin(pkin(11));
t261 = cos(pkin(11));
t313 = t262 * t266;
t241 = t259 * t269 + t261 * t313;
t216 = t241 * t252 + t261 * t274;
t319 = t260 * t261;
t217 = t241 * t278 - t252 * t319;
t312 = t262 * t269;
t240 = t259 * t266 - t261 * t312;
t101 = t190 * t240 - t191 * t216 + t192 * t217;
t242 = t259 * t312 + t261 * t266;
t264 = sin(qJ(6));
t267 = cos(qJ(6));
t186 = -t217 * t264 + t240 * t267;
t187 = t217 * t267 + t240 * t264;
t119 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t216;
t121 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t216;
t123 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t216;
t58 = t119 * t216 + t121 * t186 + t123 * t187;
t243 = -t259 * t313 + t261 * t269;
t320 = t259 * t260;
t219 = t243 * t278 + t252 * t320;
t188 = -t219 * t264 + t242 * t267;
t189 = t219 * t267 + t242 * t264;
t218 = t243 * t252 - t259 * t274;
t120 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t218;
t122 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t218;
t124 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t218;
t59 = t120 * t216 + t122 * t186 + t124 * t187;
t222 = -t233 * t264 - t267 * t315;
t223 = t233 * t267 - t264 * t315;
t141 = Icges(7,5) * t223 + Icges(7,6) * t222 + Icges(7,3) * t232;
t142 = Icges(7,4) * t223 + Icges(7,2) * t222 + Icges(7,6) * t232;
t143 = Icges(7,1) * t223 + Icges(7,4) * t222 + Icges(7,5) * t232;
t70 = t141 * t216 + t142 * t186 + t143 * t187;
t11 = t58 * t240 + t59 * t242 - t315 * t70;
t145 = Icges(6,5) * t217 - Icges(6,6) * t216 + Icges(6,3) * t240;
t147 = Icges(6,4) * t217 - Icges(6,2) * t216 + Icges(6,6) * t240;
t149 = Icges(6,1) * t217 - Icges(6,4) * t216 + Icges(6,5) * t240;
t77 = t145 * t240 - t147 * t216 + t149 * t217;
t146 = Icges(6,5) * t219 - Icges(6,6) * t218 + Icges(6,3) * t242;
t148 = Icges(6,4) * t219 - Icges(6,2) * t218 + Icges(6,6) * t242;
t150 = Icges(6,1) * t219 - Icges(6,4) * t218 + Icges(6,5) * t242;
t78 = t146 * t240 - t148 * t216 + t150 * t217;
t338 = -t101 * t315 + t77 * t240 + t78 * t242 + t11;
t102 = t190 * t242 - t191 * t218 + t192 * t219;
t60 = t119 * t218 + t121 * t188 + t123 * t189;
t61 = t120 * t218 + t122 * t188 + t124 * t189;
t71 = t141 * t218 + t142 * t188 + t143 * t189;
t12 = t60 * t240 + t61 * t242 - t315 * t71;
t79 = t145 * t242 - t147 * t218 + t149 * t219;
t80 = t146 * t242 - t148 * t218 + t150 * t219;
t337 = -t102 * t315 + t79 * t240 + t80 * t242 + t12;
t15 = t70 * t262 + (t259 * t59 - t261 * t58) * t260;
t336 = t15 + t101 * t262 + (t259 * t78 - t261 * t77) * t260;
t16 = t71 * t262 + (t259 * t61 - t261 * t60) * t260;
t335 = t16 + t102 * t262 + (t259 * t80 - t261 * t79) * t260;
t106 = -t190 * t315 - t232 * t191 + t233 * t192;
t62 = t119 * t232 + t121 * t222 + t123 * t223;
t63 = t120 * t232 + t122 * t222 + t124 * t223;
t76 = t141 * t232 + t142 * t222 + t143 * t223;
t21 = t62 * t240 + t63 * t242 - t315 * t76;
t87 = -t145 * t315 - t232 * t147 + t233 * t149;
t88 = -t146 * t315 - t232 * t148 + t233 * t150;
t334 = -t106 * t315 + t87 * t240 + t88 * t242 + t21;
t23 = t76 * t262 + (t259 * t63 - t261 * t62) * t260;
t333 = t23 + t106 * t262 + (t259 * t88 - t261 * t87) * t260;
t125 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t216;
t310 = pkin(5) * t217 + pkin(10) * t216 + t125;
t144 = rSges(7,1) * t223 + rSges(7,2) * t222 + rSges(7,3) * t232;
t332 = pkin(5) * t233 + pkin(10) * t232 + t144;
t331 = t260 ^ 2;
t330 = t216 / 0.2e1;
t329 = t218 / 0.2e1;
t328 = t232 / 0.2e1;
t327 = t240 / 0.2e1;
t326 = t242 / 0.2e1;
t325 = t262 / 0.2e1;
t268 = cos(qJ(3));
t323 = t268 * pkin(3);
t265 = sin(qJ(3));
t318 = t260 * t265;
t316 = t260 * t268;
t314 = t262 * t265;
t311 = t310 * t242;
t126 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t218;
t309 = pkin(5) * t219 + pkin(10) * t218 + t126;
t254 = sin(t258);
t284 = pkin(4) * t254;
t255 = cos(t258);
t298 = pkin(4) * t255;
t134 = pkin(9) * t240 + t241 * t298 - t284 * t319;
t294 = t261 * t318;
t167 = -pkin(3) * t294 + qJ(4) * t240 + t241 * t323;
t138 = t242 * t167;
t308 = t242 * t134 + t138;
t135 = pkin(9) * t242 + t243 * t298 + t284 * t320;
t295 = t259 * t318;
t168 = pkin(3) * t295 + qJ(4) * t242 + t243 * t323;
t307 = -t135 - t168;
t210 = pkin(3) * t314 + (-qJ(4) * t269 + t266 * t323) * t260;
t305 = t167 * t315 + t240 * t210;
t221 = t243 * pkin(2) + t242 * pkin(8);
t215 = t262 * t221;
t304 = t262 * t168 + t215;
t226 = -t243 * t254 + t255 * t320;
t227 = t243 * t255 + t254 * t320;
t165 = rSges(5,1) * t227 + rSges(5,2) * t226 + rSges(5,3) * t242;
t303 = -t165 - t168;
t220 = t241 * pkin(2) + t240 * pkin(8);
t302 = -t167 - t220;
t153 = rSges(6,1) * t217 - rSges(6,2) * t216 + rSges(6,3) * t240;
t193 = t233 * rSges(6,1) - t232 * rSges(6,2) - rSges(6,3) * t315;
t115 = t153 * t315 + t240 * t193;
t183 = t284 * t262 + (-pkin(9) * t269 + t266 * t298) * t260;
t301 = -t183 - t210;
t237 = -t254 * t317 + t255 * t262;
t238 = t254 * t262 + t255 * t317;
t197 = t238 * rSges(5,1) + t237 * rSges(5,2) - rSges(5,3) * t315;
t300 = -t197 - t210;
t299 = t220 * t320 + t221 * t319;
t296 = -m(5) - m(6) - m(7);
t293 = t262 * t135 + t304;
t292 = -t134 + t302;
t154 = rSges(6,1) * t219 - rSges(6,2) * t218 + rSges(6,3) * t242;
t291 = -t154 + t307;
t290 = -t193 + t301;
t287 = -t315 / 0.2e1;
t286 = t240 * t338 + t337 * t242;
t244 = t262 * t268 - t265 * t317;
t245 = t266 * t316 + t314;
t211 = t245 * rSges(4,1) + t244 * rSges(4,2) - rSges(4,3) * t315;
t246 = (pkin(2) * t266 - pkin(8) * t269) * t260;
t283 = (-t211 - t246) * t260;
t282 = t307 - t309;
t281 = t134 * t315 + t240 * t183 + t305;
t280 = t301 - t332;
t279 = t167 * t320 + t168 * t319 + t299;
t72 = t240 * t332 + t310 * t315;
t277 = (-t246 + t300) * t260;
t18 = t216 * t62 + t218 * t63 + t232 * t76;
t3 = t216 * t58 + t218 * t59 + t232 * t70;
t4 = t216 * t60 + t218 * t61 + t232 * t71;
t276 = t11 * t330 + t12 * t329 + t18 * t287 + t21 * t328 + t3 * t327 + t4 * t326;
t275 = (-t246 + t290) * t260;
t273 = t134 * t320 + t135 * t319 + t279;
t272 = (-t246 + t280) * t260;
t271 = -t315 * t334 + t286;
t270 = t336 * t327 + t335 * t326 + t334 * t325 + t337 * t320 / 0.2e1 - t338 * t319 / 0.2e1 + t333 * t287;
t239 = t262 * rSges(3,3) + (rSges(3,1) * t266 + rSges(3,2) * t269) * t260;
t236 = Icges(3,5) * t262 + (Icges(3,1) * t266 + Icges(3,4) * t269) * t260;
t235 = Icges(3,6) * t262 + (Icges(3,4) * t266 + Icges(3,2) * t269) * t260;
t234 = Icges(3,3) * t262 + (Icges(3,5) * t266 + Icges(3,6) * t269) * t260;
t231 = t243 * t268 + t295;
t230 = -t243 * t265 + t259 * t316;
t229 = t241 * t268 - t294;
t228 = -t241 * t265 - t261 * t316;
t225 = t241 * t255 - t254 * t319;
t224 = -t241 * t254 - t255 * t319;
t209 = Icges(4,1) * t245 + Icges(4,4) * t244 - Icges(4,5) * t315;
t208 = Icges(4,4) * t245 + Icges(4,2) * t244 - Icges(4,6) * t315;
t207 = Icges(4,5) * t245 + Icges(4,6) * t244 - Icges(4,3) * t315;
t206 = rSges(3,1) * t243 - rSges(3,2) * t242 + rSges(3,3) * t320;
t205 = rSges(3,1) * t241 - rSges(3,2) * t240 - rSges(3,3) * t319;
t204 = Icges(3,1) * t243 - Icges(3,4) * t242 + Icges(3,5) * t320;
t203 = Icges(3,1) * t241 - Icges(3,4) * t240 - Icges(3,5) * t319;
t202 = Icges(3,4) * t243 - Icges(3,2) * t242 + Icges(3,6) * t320;
t201 = Icges(3,4) * t241 - Icges(3,2) * t240 - Icges(3,6) * t319;
t200 = Icges(3,5) * t243 - Icges(3,6) * t242 + Icges(3,3) * t320;
t199 = Icges(3,5) * t241 - Icges(3,6) * t240 - Icges(3,3) * t319;
t196 = Icges(5,1) * t238 + Icges(5,4) * t237 - Icges(5,5) * t315;
t195 = Icges(5,4) * t238 + Icges(5,2) * t237 - Icges(5,6) * t315;
t194 = Icges(5,5) * t238 + Icges(5,6) * t237 - Icges(5,3) * t315;
t181 = -t205 * t262 - t239 * t319;
t180 = t206 * t262 - t239 * t320;
t177 = rSges(4,1) * t231 + rSges(4,2) * t230 + rSges(4,3) * t242;
t176 = rSges(4,1) * t229 + rSges(4,2) * t228 + rSges(4,3) * t240;
t175 = Icges(4,1) * t231 + Icges(4,4) * t230 + Icges(4,5) * t242;
t174 = Icges(4,1) * t229 + Icges(4,4) * t228 + Icges(4,5) * t240;
t173 = Icges(4,4) * t231 + Icges(4,2) * t230 + Icges(4,6) * t242;
t172 = Icges(4,4) * t229 + Icges(4,2) * t228 + Icges(4,6) * t240;
t171 = Icges(4,5) * t231 + Icges(4,6) * t230 + Icges(4,3) * t242;
t170 = Icges(4,5) * t229 + Icges(4,6) * t228 + Icges(4,3) * t240;
t164 = rSges(5,1) * t225 + rSges(5,2) * t224 + rSges(5,3) * t240;
t163 = Icges(5,1) * t227 + Icges(5,4) * t226 + Icges(5,5) * t242;
t162 = Icges(5,1) * t225 + Icges(5,4) * t224 + Icges(5,5) * t240;
t161 = Icges(5,4) * t227 + Icges(5,2) * t226 + Icges(5,6) * t242;
t160 = Icges(5,4) * t225 + Icges(5,2) * t224 + Icges(5,6) * t240;
t159 = Icges(5,5) * t227 + Icges(5,6) * t226 + Icges(5,3) * t242;
t158 = Icges(5,5) * t225 + Icges(5,6) * t224 + Icges(5,3) * t240;
t139 = (t205 * t259 + t206 * t261) * t260;
t137 = t242 * t153;
t128 = -t177 * t315 - t242 * t211;
t127 = t176 * t315 + t240 * t211;
t117 = -t207 * t315 + t244 * t208 + t245 * t209;
t116 = -t154 * t315 - t242 * t193;
t113 = t176 * t242 - t177 * t240;
t112 = (-t176 - t220) * t262 + t261 * t283;
t111 = t262 * t177 + t259 * t283 + t215;
t110 = -t194 * t315 + t237 * t195 + t238 * t196;
t109 = t207 * t242 + t208 * t230 + t209 * t231;
t108 = t207 * t240 + t208 * t228 + t209 * t229;
t107 = -t154 * t240 + t137;
t105 = t194 * t242 + t195 * t226 + t196 * t227;
t104 = t194 * t240 + t195 * t224 + t196 * t225;
t103 = (t176 * t259 + t177 * t261) * t260 + t299;
t100 = -t171 * t315 + t244 * t173 + t245 * t175;
t99 = -t170 * t315 + t244 * t172 + t245 * t174;
t98 = -t159 * t315 + t237 * t161 + t238 * t163;
t97 = -t158 * t315 + t237 * t160 + t238 * t162;
t96 = t242 * t300 + t303 * t315;
t95 = t164 * t315 + t240 * t197 + t305;
t94 = t171 * t242 + t173 * t230 + t175 * t231;
t93 = t170 * t242 + t172 * t230 + t174 * t231;
t92 = t171 * t240 + t173 * t228 + t175 * t229;
t91 = t170 * t240 + t172 * t228 + t174 * t229;
t90 = t126 * t232 - t144 * t218;
t89 = -t125 * t232 + t144 * t216;
t86 = (-t164 + t302) * t262 + t261 * t277;
t85 = t262 * t165 + t259 * t277 + t304;
t84 = t159 * t242 + t161 * t226 + t163 * t227;
t83 = t158 * t242 + t160 * t226 + t162 * t227;
t82 = t159 * t240 + t161 * t224 + t163 * t225;
t81 = t158 * t240 + t160 * t224 + t162 * t225;
t75 = t125 * t218 - t126 * t216;
t74 = t242 * t164 + t240 * t303 + t138;
t73 = -t242 * t332 - t309 * t315;
t69 = (t164 * t259 + t165 * t261) * t260 + t279;
t68 = -t240 * t309 + t311;
t67 = t242 * t290 + t291 * t315;
t66 = t281 + t115;
t65 = (-t153 + t292) * t262 + t261 * t275;
t64 = t262 * t154 + t259 * t275 + t293;
t57 = t240 * t291 + t137 + t308;
t56 = (t153 * t259 + t154 * t261) * t260 + t273;
t55 = t117 * t262 + (t100 * t259 - t261 * t99) * t260;
t54 = t100 * t242 - t117 * t315 + t99 * t240;
t53 = t242 * t280 + t282 * t315;
t52 = t72 + t281;
t51 = (t292 - t310) * t262 + t261 * t272;
t50 = t259 * t272 + t262 * t309 + t293;
t49 = t110 * t262 + (t259 * t98 - t261 * t97) * t260;
t48 = t109 * t262 + (t259 * t94 - t261 * t93) * t260;
t47 = t108 * t262 + (t259 * t92 - t261 * t91) * t260;
t46 = -t110 * t315 + t97 * t240 + t98 * t242;
t45 = -t109 * t315 + t93 * t240 + t94 * t242;
t44 = -t108 * t315 + t91 * t240 + t92 * t242;
t39 = t105 * t262 + (t259 * t84 - t261 * t83) * t260;
t38 = t104 * t262 + (t259 * t82 - t261 * t81) * t260;
t37 = -t105 * t315 + t83 * t240 + t84 * t242;
t36 = -t104 * t315 + t81 * t240 + t82 * t242;
t35 = t240 * t282 + t308 + t311;
t34 = (t259 * t310 + t261 * t309) * t260 + t273;
t1 = [m(2) + m(3) + m(4) - t296; m(3) * t139 + m(4) * t103 + m(5) * t69 + m(6) * t56 + m(7) * t34; m(7) * (t34 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t56 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (t69 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t103 ^ 2 + t111 ^ 2 + t112 ^ 2) + m(3) * (t139 ^ 2 + t180 ^ 2 + t181 ^ 2) + (t48 + t39 + (t200 * t320 - t202 * t242 + t204 * t243) * t320 + t335) * t320 + (-t47 - t38 + (-t199 * t319 - t201 * t240 + t203 * t241) * t319 + (-t199 * t320 + t200 * t319 + t201 * t242 + t202 * t240 - t203 * t243 - t204 * t241) * t320 - t336) * t319 + (-(-t234 * t319 - t240 * t235 + t241 * t236) * t319 + (t234 * t320 - t235 * t242 + t236 * t243) * t320 + t49 + t55 + ((t202 * t269 + t204 * t266) * t259 - (t201 * t269 + t203 * t266) * t261) * t331 + ((-t199 * t261 + t200 * t259 + t235 * t269 + t236 * t266) * t260 + t262 * t234) * t262 + t333) * t262; m(4) * t113 + m(5) * t74 + m(6) * t57 + m(7) * t35; (t39 / 0.2e1 + t48 / 0.2e1) * t242 + (t38 / 0.2e1 + t47 / 0.2e1) * t240 + ((-t49 / 0.2e1 - t55 / 0.2e1) * t269 + (-t36 / 0.2e1 - t44 / 0.2e1) * t261 + (t37 / 0.2e1 + t45 / 0.2e1) * t259) * t260 + m(4) * (t103 * t113 + t111 * t128 + t112 * t127) + m(6) * (t56 * t57 + t64 * t67 + t65 * t66) + m(7) * (t34 * t35 + t50 * t53 + t51 * t52) + m(5) * (t69 * t74 + t85 * t96 + t86 * t95) + (t46 / 0.2e1 + t54 / 0.2e1) * t262 + t270; (t45 + t37) * t242 + (t36 + t44) * t240 + (-t46 - t54 - t334) * t315 + m(7) * (t35 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t57 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t74 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(4) * (t113 ^ 2 + t127 ^ 2 + t128 ^ 2) + t286; t296 * t315; m(7) * (t240 * t50 + t242 * t51 - t315 * t34) + m(6) * (t240 * t64 + t242 * t65 - t315 * t56) + m(5) * (t240 * t85 + t242 * t86 - t315 * t69); m(7) * (t240 * t53 + t242 * t52 - t315 * t35) + m(6) * (t240 * t67 + t242 * t66 - t315 * t57) + m(5) * (t240 * t96 + t242 * t95 - t315 * t74); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t269 ^ 2 * t331 + t240 ^ 2 + t242 ^ 2); m(6) * t107 + m(7) * t68; m(7) * (t34 * t68 + t50 * t73 + t51 * t72) + m(6) * (t107 * t56 + t115 * t65 + t116 * t64) + t270; m(7) * (t35 * t68 + t52 * t72 + t53 * t73) + m(6) * (t107 * t57 + t115 * t66 + t116 * t67) + t271; m(6) * (-t107 * t315 + t115 * t242 + t116 * t240) + m(7) * (t73 * t240 + t72 * t242 - t315 * t68); m(7) * (t68 ^ 2 + t72 ^ 2 + t73 ^ 2) + m(6) * (t107 ^ 2 + t115 ^ 2 + t116 ^ 2) + t271; m(7) * t75; m(7) * (t34 * t75 + t50 * t90 + t51 * t89) + t18 * t325 + t23 * t328 + t16 * t329 + t15 * t330 + (t259 * t4 / 0.2e1 - t261 * t3 / 0.2e1) * t260; m(7) * (t35 * t75 + t52 * t89 + t53 * t90) + t276; m(7) * (t90 * t240 + t89 * t242 - t315 * t75); m(7) * (t68 * t75 + t72 * t89 + t73 * t90) + t276; t218 * t4 + t216 * t3 + t232 * t18 + m(7) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
