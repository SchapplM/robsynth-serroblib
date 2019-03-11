% Calculate joint inertia matrix for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:47
% EndTime: 2019-03-08 20:31:57
% DurationCPUTime: 5.44s
% Computational Cost: add. (28046->531), mult. (36492->780), div. (0->0), fcn. (45622->14), ass. (0->247)
t243 = cos(pkin(6));
t242 = cos(pkin(11));
t248 = cos(qJ(2));
t285 = t248 * t242;
t239 = sin(pkin(11));
t246 = sin(qJ(2));
t287 = t246 * t239;
t221 = -t243 * t285 + t287;
t286 = t248 * t239;
t288 = t242 * t246;
t223 = t243 * t286 + t288;
t240 = sin(pkin(6));
t289 = t240 * t248;
t222 = t243 * t288 + t286;
t237 = pkin(12) + qJ(4);
t261 = qJ(5) + t237;
t231 = sin(t261);
t257 = cos(t261);
t291 = t240 * t242;
t198 = t222 * t257 - t231 * t291;
t245 = sin(qJ(6));
t247 = cos(qJ(6));
t165 = -t198 * t245 + t221 * t247;
t166 = t198 * t247 + t221 * t245;
t253 = t240 * t257;
t197 = t222 * t231 + t242 * t253;
t103 = Icges(7,5) * t166 + Icges(7,6) * t165 + Icges(7,3) * t197;
t105 = Icges(7,4) * t166 + Icges(7,2) * t165 + Icges(7,6) * t197;
t107 = Icges(7,1) * t166 + Icges(7,4) * t165 + Icges(7,5) * t197;
t51 = t103 * t197 + t105 * t165 + t107 * t166;
t224 = -t243 * t287 + t285;
t292 = t239 * t240;
t200 = t224 * t257 + t231 * t292;
t167 = -t200 * t245 + t223 * t247;
t168 = t200 * t247 + t223 * t245;
t199 = t224 * t231 - t239 * t253;
t104 = Icges(7,5) * t168 + Icges(7,6) * t167 + Icges(7,3) * t199;
t106 = Icges(7,4) * t168 + Icges(7,2) * t167 + Icges(7,6) * t199;
t108 = Icges(7,1) * t168 + Icges(7,4) * t167 + Icges(7,5) * t199;
t52 = t104 * t197 + t106 * t165 + t108 * t166;
t212 = t231 * t243 + t246 * t253;
t201 = -t212 * t245 - t247 * t289;
t202 = t212 * t247 - t245 * t289;
t290 = t240 * t246;
t211 = t231 * t290 - t243 * t257;
t122 = Icges(7,5) * t202 + Icges(7,6) * t201 + Icges(7,3) * t211;
t123 = Icges(7,4) * t202 + Icges(7,2) * t201 + Icges(7,6) * t211;
t124 = Icges(7,1) * t202 + Icges(7,4) * t201 + Icges(7,5) * t211;
t62 = t122 * t197 + t123 * t165 + t124 * t166;
t11 = t221 * t51 + t223 * t52 - t289 * t62;
t128 = Icges(6,5) * t198 - Icges(6,6) * t197 + Icges(6,3) * t221;
t130 = Icges(6,4) * t198 - Icges(6,2) * t197 + Icges(6,6) * t221;
t132 = Icges(6,1) * t198 - Icges(6,4) * t197 + Icges(6,5) * t221;
t70 = t128 * t221 - t130 * t197 + t132 * t198;
t129 = Icges(6,5) * t200 - Icges(6,6) * t199 + Icges(6,3) * t223;
t131 = Icges(6,4) * t200 - Icges(6,2) * t199 + Icges(6,6) * t223;
t133 = Icges(6,1) * t200 - Icges(6,4) * t199 + Icges(6,5) * t223;
t71 = t129 * t221 - t131 * t197 + t133 * t198;
t169 = Icges(6,5) * t212 - Icges(6,6) * t211 - Icges(6,3) * t289;
t170 = Icges(6,4) * t212 - Icges(6,2) * t211 - Icges(6,6) * t289;
t171 = Icges(6,1) * t212 - Icges(6,4) * t211 - Icges(6,5) * t289;
t87 = t169 * t221 - t170 * t197 + t171 * t198;
t314 = t221 * t70 + t223 * t71 - t289 * t87 + t11;
t53 = t103 * t199 + t105 * t167 + t107 * t168;
t54 = t104 * t199 + t106 * t167 + t108 * t168;
t63 = t122 * t199 + t123 * t167 + t124 * t168;
t12 = t221 * t53 + t223 * t54 - t289 * t63;
t72 = t128 * t223 - t130 * t199 + t132 * t200;
t73 = t129 * t223 - t131 * t199 + t133 * t200;
t88 = t169 * t223 - t170 * t199 + t171 * t200;
t313 = t221 * t72 + t223 * t73 - t289 * t88 + t12;
t15 = t243 * t62 + (t239 * t52 - t242 * t51) * t240;
t312 = t15 + t243 * t87 + (t239 * t71 - t242 * t70) * t240;
t16 = t243 * t63 + (t239 * t54 - t242 * t53) * t240;
t311 = t16 + t243 * t88 + (t239 * t73 - t242 * t72) * t240;
t55 = t103 * t211 + t105 * t201 + t107 * t202;
t56 = t104 * t211 + t106 * t201 + t108 * t202;
t67 = t122 * t211 + t123 * t201 + t124 * t202;
t21 = t221 * t55 + t223 * t56 - t289 * t67;
t80 = -t128 * t289 - t130 * t211 + t132 * t212;
t81 = -t129 * t289 - t131 * t211 + t133 * t212;
t91 = -t169 * t289 - t170 * t211 + t171 * t212;
t310 = t221 * t80 + t223 * t81 - t289 * t91 + t21;
t23 = t243 * t67 + (t239 * t56 - t242 * t55) * t240;
t309 = t23 + t243 * t91 + (t239 * t81 - t242 * t80) * t240;
t109 = rSges(7,1) * t166 + rSges(7,2) * t165 + rSges(7,3) * t197;
t283 = pkin(5) * t198 + pkin(10) * t197 + t109;
t125 = rSges(7,1) * t202 + rSges(7,2) * t201 + rSges(7,3) * t211;
t308 = pkin(5) * t212 + pkin(10) * t211 + t125;
t307 = t240 ^ 2;
t238 = sin(pkin(12));
t241 = cos(pkin(12));
t219 = -t238 * t290 + t241 * t243;
t293 = t238 * t243;
t220 = t241 * t290 + t293;
t178 = Icges(4,5) * t220 + Icges(4,6) * t219 - Icges(4,3) * t289;
t214 = Icges(3,6) * t243 + (Icges(3,4) * t246 + Icges(3,2) * t248) * t240;
t306 = t178 - t214;
t305 = t197 / 0.2e1;
t304 = t199 / 0.2e1;
t303 = t211 / 0.2e1;
t302 = t221 / 0.2e1;
t301 = t223 / 0.2e1;
t300 = t239 / 0.2e1;
t299 = -t242 / 0.2e1;
t298 = t243 / 0.2e1;
t296 = t241 * pkin(3);
t294 = t283 * t223;
t110 = rSges(7,1) * t168 + rSges(7,2) * t167 + rSges(7,3) * t199;
t282 = pkin(5) * t200 + pkin(10) * t199 + t110;
t233 = sin(t237);
t260 = pkin(4) * t233;
t234 = cos(t237);
t274 = pkin(4) * t234;
t116 = pkin(9) * t221 + t222 * t274 - t260 * t291;
t163 = t260 * t243 + (-pkin(9) * t248 + t246 * t274) * t240;
t281 = t116 * t289 + t163 * t221;
t117 = pkin(9) * t223 + t224 * t274 + t260 * t292;
t135 = rSges(6,1) * t200 - rSges(6,2) * t199 + rSges(6,3) * t223;
t280 = -t117 - t135;
t271 = t238 * t292;
t148 = pkin(3) * t271 + pkin(8) * t223 + t224 * t296;
t196 = pkin(2) * t224 + qJ(3) * t223;
t194 = t243 * t196;
t278 = t148 * t243 + t194;
t270 = t238 * t291;
t147 = -pkin(3) * t270 + pkin(8) * t221 + t222 * t296;
t195 = pkin(2) * t222 + qJ(3) * t221;
t277 = -t147 - t195;
t134 = rSges(6,1) * t198 - rSges(6,2) * t197 + rSges(6,3) * t221;
t172 = rSges(6,1) * t212 - rSges(6,2) * t211 - rSges(6,3) * t289;
t98 = t134 * t289 + t172 * t221;
t225 = (pkin(2) * t246 - qJ(3) * t248) * t240;
t276 = -pkin(3) * t293 - (-pkin(8) * t248 + t246 * t296) * t240 - t225;
t275 = t195 * t292 + t196 * t291;
t272 = -m(4) - m(5) - m(6) - m(7);
t269 = -t117 - t282;
t268 = t117 * t243 + t278;
t267 = -t116 + t277;
t266 = -t163 + t276;
t263 = -t289 / 0.2e1;
t262 = t221 * t314 + t223 * t313;
t259 = (-rSges(4,1) * t220 - rSges(4,2) * t219 + rSges(4,3) * t289 - t225) * t240;
t258 = t147 * t292 + t148 * t291 + t275;
t64 = t221 * t308 + t283 * t289;
t216 = -t233 * t290 + t234 * t243;
t217 = t233 * t243 + t234 * t290;
t176 = rSges(5,1) * t217 + rSges(5,2) * t216 - rSges(5,3) * t289;
t256 = (-t176 + t276) * t240;
t18 = t197 * t55 + t199 * t56 + t211 * t67;
t3 = t197 * t51 + t199 * t52 + t211 * t62;
t4 = t197 * t53 + t199 * t54 + t211 * t63;
t255 = t11 * t305 + t12 * t304 + t18 * t263 + t21 * t303 + t3 * t302 + t301 * t4;
t254 = (-t172 + t266) * t240;
t252 = t116 * t292 + t117 * t291 + t258;
t251 = (t266 - t308) * t240;
t250 = -t289 * t310 + t262;
t249 = t312 * t302 + t311 * t301 + t310 * t298 + t313 * t292 / 0.2e1 - t314 * t291 / 0.2e1 + t309 * t263;
t218 = t243 * rSges(3,3) + (rSges(3,1) * t246 + rSges(3,2) * t248) * t240;
t215 = Icges(3,5) * t243 + (Icges(3,1) * t246 + Icges(3,4) * t248) * t240;
t213 = Icges(3,3) * t243 + (Icges(3,5) * t246 + Icges(3,6) * t248) * t240;
t210 = t224 * t241 + t271;
t209 = -t224 * t238 + t241 * t292;
t208 = t222 * t241 - t270;
t207 = -t222 * t238 - t241 * t291;
t206 = t224 * t234 + t233 * t292;
t205 = -t224 * t233 + t234 * t292;
t204 = t222 * t234 - t233 * t291;
t203 = -t222 * t233 - t234 * t291;
t190 = rSges(3,1) * t224 - rSges(3,2) * t223 + rSges(3,3) * t292;
t189 = rSges(3,1) * t222 - rSges(3,2) * t221 - rSges(3,3) * t291;
t186 = Icges(3,1) * t224 - Icges(3,4) * t223 + Icges(3,5) * t292;
t185 = Icges(3,1) * t222 - Icges(3,4) * t221 - Icges(3,5) * t291;
t184 = Icges(3,4) * t224 - Icges(3,2) * t223 + Icges(3,6) * t292;
t183 = Icges(3,4) * t222 - Icges(3,2) * t221 - Icges(3,6) * t291;
t182 = Icges(3,5) * t224 - Icges(3,6) * t223 + Icges(3,3) * t292;
t181 = Icges(3,5) * t222 - Icges(3,6) * t221 - Icges(3,3) * t291;
t180 = Icges(4,1) * t220 + Icges(4,4) * t219 - Icges(4,5) * t289;
t179 = Icges(4,4) * t220 + Icges(4,2) * t219 - Icges(4,6) * t289;
t175 = Icges(5,1) * t217 + Icges(5,4) * t216 - Icges(5,5) * t289;
t174 = Icges(5,4) * t217 + Icges(5,2) * t216 - Icges(5,6) * t289;
t173 = Icges(5,5) * t217 + Icges(5,6) * t216 - Icges(5,3) * t289;
t161 = -t189 * t243 - t218 * t291;
t160 = t190 * t243 - t218 * t292;
t157 = rSges(4,1) * t210 + rSges(4,2) * t209 + rSges(4,3) * t223;
t156 = rSges(4,1) * t208 + rSges(4,2) * t207 + rSges(4,3) * t221;
t154 = Icges(4,1) * t210 + Icges(4,4) * t209 + Icges(4,5) * t223;
t153 = Icges(4,1) * t208 + Icges(4,4) * t207 + Icges(4,5) * t221;
t152 = Icges(4,4) * t210 + Icges(4,2) * t209 + Icges(4,6) * t223;
t151 = Icges(4,4) * t208 + Icges(4,2) * t207 + Icges(4,6) * t221;
t150 = Icges(4,5) * t210 + Icges(4,6) * t209 + Icges(4,3) * t223;
t149 = Icges(4,5) * t208 + Icges(4,6) * t207 + Icges(4,3) * t221;
t146 = rSges(5,1) * t206 + rSges(5,2) * t205 + rSges(5,3) * t223;
t145 = rSges(5,1) * t204 + rSges(5,2) * t203 + rSges(5,3) * t221;
t143 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t223;
t142 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t221;
t141 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t223;
t140 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t221;
t139 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t223;
t138 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t221;
t120 = (t189 * t239 + t190 * t242) * t240;
t119 = t223 * t134;
t111 = t223 * t116;
t102 = -t146 * t289 - t176 * t223;
t101 = t145 * t289 + t176 * t221;
t99 = -t135 * t289 - t223 * t172;
t96 = (-t156 - t195) * t243 + t242 * t259;
t95 = t157 * t243 + t239 * t259 + t194;
t94 = -t173 * t289 + t174 * t216 + t175 * t217;
t93 = t145 * t223 - t146 * t221;
t92 = -t135 * t221 + t119;
t90 = t173 * t223 + t174 * t205 + t175 * t206;
t89 = t173 * t221 + t174 * t203 + t175 * t204;
t86 = (t156 * t239 + t157 * t242) * t240 + t275;
t85 = -t139 * t289 + t141 * t216 + t143 * t217;
t84 = -t138 * t289 + t140 * t216 + t142 * t217;
t83 = t110 * t211 - t125 * t199;
t82 = -t109 * t211 + t125 * t197;
t79 = t139 * t223 + t141 * t205 + t143 * t206;
t78 = t138 * t223 + t140 * t205 + t142 * t206;
t77 = t139 * t221 + t141 * t203 + t143 * t204;
t76 = t138 * t221 + t140 * t203 + t142 * t204;
t75 = (-t145 + t277) * t243 + t242 * t256;
t74 = t146 * t243 + t239 * t256 + t278;
t69 = t280 * t289 + (-t163 - t172) * t223;
t68 = t98 + t281;
t66 = t109 * t199 - t110 * t197;
t65 = -t223 * t308 - t282 * t289;
t61 = (t145 * t239 + t146 * t242) * t240 + t258;
t60 = t221 * t280 + t111 + t119;
t59 = -t221 * t282 + t294;
t58 = (-t134 + t267) * t243 + t242 * t254;
t57 = t135 * t243 + t239 * t254 + t268;
t50 = t269 * t289 + (-t163 - t308) * t223;
t49 = t64 + t281;
t48 = (t134 * t239 + t135 * t242) * t240 + t252;
t47 = (t267 - t283) * t243 + t242 * t251;
t46 = t239 * t251 + t243 * t282 + t268;
t45 = t221 * t269 + t111 + t294;
t44 = t243 * t94 + (t239 * t85 - t242 * t84) * t240;
t43 = t221 * t84 + t223 * t85 - t289 * t94;
t38 = t243 * t90 + (t239 * t79 - t242 * t78) * t240;
t37 = t243 * t89 + (t239 * t77 - t242 * t76) * t240;
t36 = t221 * t78 + t223 * t79 - t289 * t90;
t35 = t221 * t76 + t223 * t77 - t289 * t89;
t32 = (t239 * t283 + t242 * t282) * t240 + t252;
t1 = [m(2) + m(3) - t272; m(3) * t120 + m(4) * t86 + m(5) * t61 + m(6) * t48 + m(7) * t32; m(7) * (t32 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t48 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(5) * (t61 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t86 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(3) * (t120 ^ 2 + t160 ^ 2 + t161 ^ 2) + (t38 + ((t150 * t223 + t152 * t209 + t154 * t210) * t239 - (t149 * t223 + t151 * t209 + t153 * t210) * t242) * t240 + (t182 * t292 - t184 * t223 + t186 * t224) * t292 + t311) * t292 + (-t37 - ((t150 * t221 + t152 * t207 + t154 * t208) * t239 - (t149 * t221 + t151 * t207 + t153 * t208) * t242) * t240 + (-t181 * t291 - t183 * t221 + t185 * t222) * t291 + (-t181 * t292 + t182 * t291 + t183 * t223 + t184 * t221 - t185 * t224 - t186 * t222) * t292 - t312) * t291 + (t44 + ((t184 * t248 + t186 * t246) * t239 - (t183 * t248 + t185 * t246) * t242) * t307 + (-t178 * t289 + t219 * t179 + t220 * t180 + (-t181 * t242 + t182 * t239 + t214 * t248 + t215 * t246) * t240 + t243 * t213) * t243 + (-t150 * t289 + t219 * t152 + t220 * t154 + t179 * t209 + t180 * t210 + t213 * t292 + t215 * t224 + t223 * t306) * t292 + (t149 * t289 - t219 * t151 - t220 * t153 - t179 * t207 - t180 * t208 + t213 * t291 - t215 * t222 - t221 * t306) * t291 + t309) * t243; t272 * t289; m(7) * (t221 * t46 + t223 * t47 - t289 * t32) + m(6) * (t221 * t57 + t223 * t58 - t289 * t48) + m(5) * (t221 * t74 + t223 * t75 - t289 * t61) + m(4) * (t221 * t95 + t223 * t96 - t289 * t86); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t248 ^ 2 * t307 + t221 ^ 2 + t223 ^ 2); m(5) * t93 + m(6) * t60 + m(7) * t45; m(7) * (t32 * t45 + t46 * t50 + t47 * t49) + m(6) * (t48 * t60 + t57 * t69 + t58 * t68) + m(5) * (t101 * t75 + t102 * t74 + t61 * t93) + (-t248 * t44 / 0.2e1 + t36 * t300 + t35 * t299) * t240 + t249 + t43 * t298 + t37 * t302 + t38 * t301; m(5) * (t101 * t223 + t102 * t221 - t289 * t93) + m(6) * (t221 * t69 + t223 * t68 - t289 * t60) + m(7) * (t221 * t50 + t223 * t49 - t289 * t45); t221 * t35 + t223 * t36 + (-t43 - t310) * t289 + m(7) * (t45 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t60 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(5) * (t101 ^ 2 + t102 ^ 2 + t93 ^ 2) + t262; m(6) * t92 + m(7) * t59; m(7) * (t32 * t59 + t46 * t65 + t47 * t64) + m(6) * (t48 * t92 + t57 * t99 + t58 * t98) + t249; m(6) * (t221 * t99 + t223 * t98 - t289 * t92) + m(7) * (t221 * t65 + t223 * t64 - t289 * t59); m(7) * (t45 * t59 + t49 * t64 + t50 * t65) + m(6) * (t60 * t92 + t68 * t98 + t69 * t99) + t250; m(7) * (t59 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(6) * (t92 ^ 2 + t98 ^ 2 + t99 ^ 2) + t250; m(7) * t66; t23 * t303 + t18 * t298 + t15 * t305 + t16 * t304 + m(7) * (t32 * t66 + t46 * t83 + t47 * t82) + (t299 * t3 + t300 * t4) * t240; m(7) * (t221 * t83 + t223 * t82 - t289 * t66); m(7) * (t45 * t66 + t49 * t82 + t50 * t83) + t255; m(7) * (t59 * t66 + t64 * t82 + t65 * t83) + t255; t199 * t4 + t197 * t3 + t211 * t18 + m(7) * (t66 ^ 2 + t82 ^ 2 + t83 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
