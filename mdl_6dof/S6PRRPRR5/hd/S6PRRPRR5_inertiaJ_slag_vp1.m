% Calculate joint inertia matrix for
% S6PRRPRR5
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:25
% EndTime: 2019-03-08 22:16:36
% DurationCPUTime: 5.25s
% Computational Cost: add. (27452->590), mult. (50301->849), div. (0->0), fcn. (64175->14), ass. (0->273)
t314 = m(5) + m(6) + m(7);
t252 = sin(pkin(11));
t255 = cos(pkin(11));
t260 = cos(qJ(2));
t256 = cos(pkin(6));
t259 = sin(qJ(2));
t298 = t256 * t259;
t231 = t252 * t260 + t255 * t298;
t258 = sin(qJ(3));
t253 = sin(pkin(6));
t307 = cos(qJ(3));
t275 = t253 * t307;
t222 = t231 * t258 + t255 * t275;
t313 = t222 / 0.2e1;
t233 = -t252 * t298 + t255 * t260;
t224 = t233 * t258 - t252 * t275;
t312 = t224 / 0.2e1;
t297 = t256 * t260;
t230 = t252 * t259 - t255 * t297;
t311 = t230 / 0.2e1;
t232 = t252 * t297 + t255 * t259;
t310 = t232 / 0.2e1;
t300 = t253 * t258;
t234 = -t256 * t307 + t259 * t300;
t309 = t234 / 0.2e1;
t308 = t256 / 0.2e1;
t254 = cos(pkin(12));
t305 = t254 * pkin(4);
t251 = sin(pkin(12));
t304 = t230 * t251;
t303 = t232 * t251;
t302 = t252 * t253;
t301 = t253 * t255;
t299 = t253 * t260;
t223 = t231 * t307 - t255 * t300;
t250 = pkin(12) + qJ(5);
t245 = sin(t250);
t273 = pkin(5) * t245;
t246 = cos(t250);
t283 = pkin(5) * t246;
t112 = pkin(10) * t222 + t283 * t223 + t273 * t230;
t247 = qJ(6) + t250;
t242 = sin(t247);
t243 = cos(t247);
t185 = -t223 * t242 + t230 * t243;
t186 = t223 * t243 + t230 * t242;
t128 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t222;
t295 = t112 + t128;
t225 = t233 * t307 + t252 * t300;
t113 = pkin(10) * t224 + t283 * t225 + t273 * t232;
t187 = -t225 * t242 + t232 * t243;
t188 = t225 * t243 + t232 * t242;
t129 = rSges(7,1) * t188 + rSges(7,2) * t187 + rSges(7,3) * t224;
t294 = t113 + t129;
t130 = pkin(4) * t304 + pkin(9) * t222 + t305 * t223;
t183 = t223 * pkin(3) + t222 * qJ(4);
t175 = t232 * t183;
t293 = t232 * t130 + t175;
t131 = pkin(4) * t303 + pkin(9) * t224 + t305 * t225;
t184 = t225 * pkin(3) + t224 * qJ(4);
t292 = -t131 - t184;
t195 = -t225 * t251 + t232 * t254;
t196 = t225 * t254 + t303;
t149 = rSges(5,1) * t196 + rSges(5,2) * t195 + rSges(5,3) * t224;
t291 = -t149 - t184;
t235 = t256 * t258 + t259 * t275;
t151 = pkin(10) * t234 + t283 * t235 - t273 * t299;
t216 = -t235 * t242 - t243 * t299;
t217 = t235 * t243 - t242 * t299;
t157 = rSges(7,1) * t217 + rSges(7,2) * t216 + rSges(7,3) * t234;
t290 = t151 + t157;
t280 = t251 * t299;
t162 = -pkin(4) * t280 + pkin(9) * t234 + t305 * t235;
t215 = t235 * pkin(3) + t234 * qJ(4);
t289 = -t162 - t215;
t220 = -t235 * t251 - t254 * t299;
t221 = t235 * t254 - t280;
t174 = rSges(5,1) * t221 + rSges(5,2) * t220 + rSges(5,3) * t234;
t288 = -t174 - t215;
t287 = t183 * t299 + t230 * t215;
t214 = pkin(2) * t233 + pkin(8) * t232;
t212 = t256 * t214;
t286 = t256 * t184 + t212;
t213 = pkin(2) * t231 + pkin(8) * t230;
t285 = -t183 - t213;
t284 = t213 * t302 + t214 * t301;
t122 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t222;
t124 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t222;
t126 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t222;
t73 = t122 * t234 + t124 * t216 + t126 * t217;
t123 = Icges(7,5) * t188 + Icges(7,6) * t187 + Icges(7,3) * t224;
t125 = Icges(7,4) * t188 + Icges(7,2) * t187 + Icges(7,6) * t224;
t127 = Icges(7,1) * t188 + Icges(7,4) * t187 + Icges(7,5) * t224;
t74 = t123 * t234 + t125 * t216 + t127 * t217;
t154 = Icges(7,5) * t217 + Icges(7,6) * t216 + Icges(7,3) * t234;
t155 = Icges(7,4) * t217 + Icges(7,2) * t216 + Icges(7,6) * t234;
t156 = Icges(7,1) * t217 + Icges(7,4) * t216 + Icges(7,5) * t234;
t90 = t154 * t234 + t155 * t216 + t156 * t217;
t26 = t222 * t73 + t224 * t74 + t234 * t90;
t61 = t122 * t222 + t124 * t185 + t126 * t186;
t62 = t123 * t222 + t125 * t185 + t127 * t186;
t81 = t154 * t222 + t155 * t185 + t156 * t186;
t7 = t222 * t61 + t224 * t62 + t234 * t81;
t63 = t122 * t224 + t124 * t187 + t126 * t188;
t64 = t123 * t224 + t125 * t187 + t127 * t188;
t82 = t154 * t224 + t155 * t187 + t156 * t188;
t8 = t222 * t63 + t224 * t64 + t234 * t82;
t281 = t222 * t7 + t224 * t8 + t234 * t26;
t279 = t256 * t131 + t286;
t278 = -t130 + t285;
t191 = -t225 * t245 + t232 * t246;
t192 = t225 * t246 + t232 * t245;
t139 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t224;
t277 = -t139 + t292;
t218 = -t235 * t245 - t246 * t299;
t219 = t235 * t246 - t245 * t299;
t161 = rSges(6,1) * t219 + rSges(6,2) * t218 + rSges(6,3) * t234;
t276 = -t161 + t289;
t274 = -t299 / 0.2e1;
t209 = t235 * rSges(4,1) - t234 * rSges(4,2) - rSges(4,3) * t299;
t236 = (pkin(2) * t259 - pkin(8) * t260) * t253;
t272 = (-t209 - t236) * t253;
t270 = t292 - t294;
t269 = t130 * t299 + t230 * t162 + t287;
t268 = t289 - t290;
t267 = t183 * t302 + t184 * t301 + t284;
t266 = (-t236 + t288) * t253;
t13 = t61 * t230 + t62 * t232 - t81 * t299;
t14 = t63 * t230 + t64 * t232 - t82 * t299;
t32 = t73 * t230 + t74 * t232 - t90 * t299;
t265 = t13 * t313 + t14 * t312 + t26 * t274 + t32 * t309 + t8 * t310 + t7 * t311;
t17 = t81 * t256 + (t252 * t62 - t255 * t61) * t253;
t18 = t82 * t256 + (t252 * t64 - t255 * t63) * t253;
t35 = t90 * t256 + (t252 * t74 - t255 * t73) * t253;
t264 = t17 * t313 + t18 * t312 + t26 * t308 + t35 * t309 + t8 * t302 / 0.2e1 - t7 * t301 / 0.2e1;
t263 = (-t236 + t276) * t253;
t262 = t130 * t302 + t131 * t301 + t267;
t261 = (-t236 + t268) * t253;
t229 = t256 * rSges(3,3) + (rSges(3,1) * t259 + rSges(3,2) * t260) * t253;
t228 = Icges(3,5) * t256 + (Icges(3,1) * t259 + Icges(3,4) * t260) * t253;
t227 = Icges(3,6) * t256 + (Icges(3,4) * t259 + Icges(3,2) * t260) * t253;
t226 = Icges(3,3) * t256 + (Icges(3,5) * t259 + Icges(3,6) * t260) * t253;
t208 = Icges(4,1) * t235 - Icges(4,4) * t234 - Icges(4,5) * t299;
t207 = Icges(4,4) * t235 - Icges(4,2) * t234 - Icges(4,6) * t299;
t206 = Icges(4,5) * t235 - Icges(4,6) * t234 - Icges(4,3) * t299;
t205 = rSges(3,1) * t233 - rSges(3,2) * t232 + rSges(3,3) * t302;
t204 = rSges(3,1) * t231 - rSges(3,2) * t230 - rSges(3,3) * t301;
t203 = Icges(3,1) * t233 - Icges(3,4) * t232 + Icges(3,5) * t302;
t202 = Icges(3,1) * t231 - Icges(3,4) * t230 - Icges(3,5) * t301;
t201 = Icges(3,4) * t233 - Icges(3,2) * t232 + Icges(3,6) * t302;
t200 = Icges(3,4) * t231 - Icges(3,2) * t230 - Icges(3,6) * t301;
t199 = Icges(3,5) * t233 - Icges(3,6) * t232 + Icges(3,3) * t302;
t198 = Icges(3,5) * t231 - Icges(3,6) * t230 - Icges(3,3) * t301;
t194 = t223 * t254 + t304;
t193 = -t223 * t251 + t230 * t254;
t190 = t223 * t246 + t230 * t245;
t189 = -t223 * t245 + t230 * t246;
t178 = -t204 * t256 - t229 * t301;
t177 = t205 * t256 - t229 * t302;
t173 = rSges(4,1) * t225 - rSges(4,2) * t224 + rSges(4,3) * t232;
t172 = rSges(4,1) * t223 - rSges(4,2) * t222 + rSges(4,3) * t230;
t171 = Icges(5,1) * t221 + Icges(5,4) * t220 + Icges(5,5) * t234;
t170 = Icges(5,4) * t221 + Icges(5,2) * t220 + Icges(5,6) * t234;
t169 = Icges(5,5) * t221 + Icges(5,6) * t220 + Icges(5,3) * t234;
t168 = Icges(4,1) * t225 - Icges(4,4) * t224 + Icges(4,5) * t232;
t167 = Icges(4,1) * t223 - Icges(4,4) * t222 + Icges(4,5) * t230;
t166 = Icges(4,4) * t225 - Icges(4,2) * t224 + Icges(4,6) * t232;
t165 = Icges(4,4) * t223 - Icges(4,2) * t222 + Icges(4,6) * t230;
t164 = Icges(4,5) * t225 - Icges(4,6) * t224 + Icges(4,3) * t232;
t163 = Icges(4,5) * t223 - Icges(4,6) * t222 + Icges(4,3) * t230;
t160 = Icges(6,1) * t219 + Icges(6,4) * t218 + Icges(6,5) * t234;
t159 = Icges(6,4) * t219 + Icges(6,2) * t218 + Icges(6,6) * t234;
t158 = Icges(6,5) * t219 + Icges(6,6) * t218 + Icges(6,3) * t234;
t153 = (t204 * t252 + t205 * t255) * t253;
t150 = t222 * t157;
t148 = rSges(5,1) * t194 + rSges(5,2) * t193 + rSges(5,3) * t222;
t147 = Icges(5,1) * t196 + Icges(5,4) * t195 + Icges(5,5) * t224;
t146 = Icges(5,1) * t194 + Icges(5,4) * t193 + Icges(5,5) * t222;
t145 = Icges(5,4) * t196 + Icges(5,2) * t195 + Icges(5,6) * t224;
t144 = Icges(5,4) * t194 + Icges(5,2) * t193 + Icges(5,6) * t222;
t143 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t224;
t142 = Icges(5,5) * t194 + Icges(5,6) * t193 + Icges(5,3) * t222;
t141 = -t173 * t299 - t232 * t209;
t140 = t172 * t299 + t230 * t209;
t138 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t222;
t137 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t224;
t136 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t222;
t135 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t224;
t134 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t222;
t133 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t224;
t132 = Icges(6,5) * t190 + Icges(6,6) * t189 + Icges(6,3) * t222;
t116 = -t206 * t299 - t234 * t207 + t235 * t208;
t115 = t234 * t129;
t114 = t224 * t128;
t111 = t172 * t232 - t173 * t230;
t110 = (-t172 - t213) * t256 + t255 * t272;
t109 = t256 * t173 + t252 * t272 + t212;
t108 = t206 * t232 - t207 * t224 + t208 * t225;
t107 = t206 * t230 - t207 * t222 + t208 * t223;
t106 = (t172 * t252 + t173 * t255) * t253 + t284;
t105 = -t164 * t299 - t234 * t166 + t235 * t168;
t104 = -t163 * t299 - t234 * t165 + t235 * t167;
t103 = t139 * t234 - t161 * t224;
t102 = -t138 * t234 + t161 * t222;
t101 = -t157 * t224 + t115;
t100 = -t128 * t234 + t150;
t99 = t169 * t234 + t170 * t220 + t171 * t221;
t98 = t164 * t232 - t166 * t224 + t168 * t225;
t97 = t163 * t232 - t165 * t224 + t167 * t225;
t96 = t164 * t230 - t166 * t222 + t168 * t223;
t95 = t163 * t230 - t165 * t222 + t167 * t223;
t94 = t158 * t234 + t159 * t218 + t160 * t219;
t93 = t138 * t224 - t139 * t222;
t92 = t288 * t232 + t291 * t299;
t91 = t148 * t299 + t230 * t174 + t287;
t89 = -t129 * t222 + t114;
t88 = t169 * t224 + t170 * t195 + t171 * t196;
t87 = t169 * t222 + t170 * t193 + t171 * t194;
t86 = (-t148 + t285) * t256 + t255 * t266;
t85 = t256 * t149 + t252 * t266 + t286;
t84 = t158 * t224 + t159 * t191 + t160 * t192;
t83 = t158 * t222 + t159 * t189 + t160 * t190;
t80 = t232 * t148 + t291 * t230 + t175;
t79 = t143 * t234 + t145 * t220 + t147 * t221;
t78 = t142 * t234 + t144 * t220 + t146 * t221;
t77 = (t148 * t252 + t149 * t255) * t253 + t267;
t76 = t133 * t234 + t135 * t218 + t137 * t219;
t75 = t132 * t234 + t134 * t218 + t136 * t219;
t72 = t143 * t224 + t145 * t195 + t147 * t196;
t71 = t142 * t224 + t144 * t195 + t146 * t196;
t70 = t143 * t222 + t145 * t193 + t147 * t194;
t69 = t142 * t222 + t144 * t193 + t146 * t194;
t68 = t133 * t224 + t135 * t191 + t137 * t192;
t67 = t132 * t224 + t134 * t191 + t136 * t192;
t66 = t133 * t222 + t135 * t189 + t137 * t190;
t65 = t132 * t222 + t134 * t189 + t136 * t190;
t60 = t276 * t232 + t277 * t299;
t59 = t138 * t299 + t230 * t161 + t269;
t58 = t234 * t113 - t290 * t224 + t115;
t57 = t222 * t151 - t295 * t234 + t150;
t56 = (-t138 + t278) * t256 + t255 * t263;
t55 = t256 * t139 + t252 * t263 + t279;
t54 = t116 * t256 + (-t104 * t255 + t105 * t252) * t253;
t53 = t104 * t230 + t105 * t232 - t116 * t299;
t52 = t224 * t112 - t294 * t222 + t114;
t51 = t232 * t138 + t277 * t230 + t293;
t50 = (t138 * t252 + t139 * t255) * t253 + t262;
t49 = t108 * t256 + (t252 * t98 - t255 * t97) * t253;
t48 = t107 * t256 + (t252 * t96 - t255 * t95) * t253;
t47 = -t108 * t299 + t97 * t230 + t98 * t232;
t46 = -t107 * t299 + t95 * t230 + t96 * t232;
t45 = t268 * t232 + t270 * t299;
t44 = t290 * t230 + t295 * t299 + t269;
t43 = (t278 - t295) * t256 + t255 * t261;
t42 = t252 * t261 + t294 * t256 + t279;
t41 = t270 * t230 + t295 * t232 + t293;
t40 = (t295 * t252 + t294 * t255) * t253 + t262;
t39 = t99 * t256 + (t252 * t79 - t255 * t78) * t253;
t38 = t78 * t230 + t79 * t232 - t99 * t299;
t37 = t94 * t256 + (t252 * t76 - t255 * t75) * t253;
t36 = t75 * t230 + t76 * t232 - t94 * t299;
t34 = t222 * t75 + t224 * t76 + t234 * t94;
t31 = t88 * t256 + (t252 * t72 - t255 * t71) * t253;
t30 = t87 * t256 + (t252 * t70 - t255 * t69) * t253;
t28 = t71 * t230 + t72 * t232 - t88 * t299;
t27 = t69 * t230 + t70 * t232 - t87 * t299;
t22 = t84 * t256 + (t252 * t68 - t255 * t67) * t253;
t21 = t83 * t256 + (t252 * t66 - t255 * t65) * t253;
t20 = t67 * t230 + t68 * t232 - t84 * t299;
t19 = t65 * t230 + t66 * t232 - t83 * t299;
t16 = t222 * t67 + t224 * t68 + t234 * t84;
t15 = t222 * t65 + t224 * t66 + t234 * t83;
t1 = [m(2) + m(3) + m(4) + t314; m(3) * t153 + m(4) * t106 + m(5) * t77 + m(6) * t50 + m(7) * t40; m(7) * (t40 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t50 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t77 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t106 ^ 2 + t109 ^ 2 + t110 ^ 2) + m(3) * (t153 ^ 2 + t177 ^ 2 + t178 ^ 2) + (t18 + t22 + t49 + t31 + (t199 * t302 - t201 * t232 + t203 * t233) * t302) * t302 + (-t17 - t21 - t48 - t30 + (-t198 * t301 - t200 * t230 + t202 * t231) * t301 + (-t198 * t302 + t199 * t301 + t200 * t232 + t201 * t230 - t202 * t233 - t203 * t231) * t302) * t301 + ((t226 * t302 - t232 * t227 + t233 * t228) * t302 - (-t226 * t301 - t230 * t227 + t231 * t228) * t301 + t35 + t37 + t54 + t39 + ((t201 * t260 + t203 * t259) * t252 - (t200 * t260 + t202 * t259) * t255) * t253 ^ 2 + ((-t198 * t255 + t199 * t252 + t227 * t260 + t228 * t259) * t253 + t256 * t226) * t256) * t256; m(4) * t111 + m(5) * t80 + m(6) * t51 + m(7) * t41; (t32 / 0.2e1 + t36 / 0.2e1 + t38 / 0.2e1 + t53 / 0.2e1) * t256 + (t18 / 0.2e1 + t22 / 0.2e1 + t31 / 0.2e1 + t49 / 0.2e1) * t232 + (t17 / 0.2e1 + t21 / 0.2e1 + t30 / 0.2e1 + t48 / 0.2e1) * t230 + m(7) * (t40 * t41 + t42 * t45 + t43 * t44) + m(6) * (t50 * t51 + t55 * t60 + t56 * t59) + m(5) * (t77 * t80 + t85 * t92 + t86 * t91) + m(4) * (t106 * t111 + t109 * t141 + t110 * t140) + ((-t35 / 0.2e1 - t37 / 0.2e1 - t39 / 0.2e1 - t54 / 0.2e1) * t260 + (-t13 / 0.2e1 - t19 / 0.2e1 - t27 / 0.2e1 - t46 / 0.2e1) * t255 + (t14 / 0.2e1 + t20 / 0.2e1 + t28 / 0.2e1 + t47 / 0.2e1) * t252) * t253; (-t32 - t36 - t38 - t53) * t299 + (t14 + t20 + t47 + t28) * t232 + (t13 + t19 + t27 + t46) * t230 + m(7) * (t41 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t51 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t80 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t111 ^ 2 + t140 ^ 2 + t141 ^ 2); t234 * t314; m(7) * (t222 * t42 + t224 * t43 + t234 * t40) + m(6) * (t222 * t55 + t224 * t56 + t234 * t50) + m(5) * (t222 * t85 + t224 * t86 + t234 * t77); m(7) * (t222 * t45 + t224 * t44 + t234 * t41) + m(6) * (t222 * t60 + t224 * t59 + t234 * t51) + m(5) * (t222 * t92 + t224 * t91 + t234 * t80); (t222 ^ 2 + t224 ^ 2 + t234 ^ 2) * t314; m(6) * t93 + m(7) * t52; t34 * t308 + t37 * t309 + t22 * t312 + t21 * t313 + (t252 * t16 / 0.2e1 - t255 * t15 / 0.2e1) * t253 + m(7) * (t40 * t52 + t42 * t58 + t43 * t57) + m(6) * (t102 * t56 + t103 * t55 + t50 * t93) + t264; t34 * t274 + t15 * t311 + t19 * t313 + t20 * t312 + t16 * t310 + t36 * t309 + m(7) * (t41 * t52 + t44 * t57 + t45 * t58) + m(6) * (t102 * t59 + t103 * t60 + t51 * t93) + t265; m(6) * (t102 * t224 + t103 * t222 + t234 * t93) + m(7) * (t222 * t58 + t224 * t57 + t234 * t52); t222 * t15 + t224 * t16 + t234 * t34 + m(7) * (t52 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t102 ^ 2 + t103 ^ 2 + t93 ^ 2) + t281; m(7) * t89; m(7) * (t100 * t43 + t101 * t42 + t40 * t89) + t264; m(7) * (t100 * t44 + t101 * t45 + t41 * t89) + t265; m(7) * (t100 * t224 + t101 * t222 + t234 * t89); m(7) * (t100 * t57 + t101 * t58 + t52 * t89) + t281; m(7) * (t100 ^ 2 + t101 ^ 2 + t89 ^ 2) + t281;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
