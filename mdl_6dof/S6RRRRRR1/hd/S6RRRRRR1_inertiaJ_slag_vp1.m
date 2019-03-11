% Calculate joint inertia matrix for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:27:52
% EndTime: 2019-03-10 03:28:01
% DurationCPUTime: 3.52s
% Computational Cost: add. (15366->443), mult. (11139->643), div. (0->0), fcn. (11561->12), ass. (0->229)
t206 = sin(qJ(1));
t209 = cos(qJ(1));
t203 = qJ(2) + qJ(3);
t190 = sin(t203);
t191 = cos(t203);
t239 = Icges(4,5) * t191 - Icges(4,6) * t190;
t130 = -Icges(4,3) * t209 + t239 * t206;
t131 = Icges(4,3) * t206 + t239 * t209;
t202 = t209 ^ 2;
t281 = Icges(4,4) * t191;
t243 = -Icges(4,2) * t190 + t281;
t133 = Icges(4,6) * t206 + t243 * t209;
t282 = Icges(4,4) * t190;
t247 = Icges(4,1) * t191 - t282;
t135 = Icges(4,5) * t206 + t247 * t209;
t231 = -t133 * t190 + t135 * t191;
t132 = -Icges(4,6) * t209 + t243 * t206;
t134 = -Icges(4,5) * t209 + t247 * t206;
t232 = t132 * t190 - t134 * t191;
t192 = qJ(4) + t203;
t186 = sin(t192);
t187 = cos(t192);
t238 = Icges(5,5) * t187 - Icges(5,6) * t186;
t122 = -Icges(5,3) * t209 + t238 * t206;
t123 = Icges(5,3) * t206 + t238 * t209;
t279 = Icges(5,4) * t187;
t242 = -Icges(5,2) * t186 + t279;
t125 = Icges(5,6) * t206 + t242 * t209;
t280 = Icges(5,4) * t186;
t246 = Icges(5,1) * t187 - t280;
t127 = Icges(5,5) * t206 + t246 * t209;
t233 = -t125 * t186 + t127 * t187;
t124 = -Icges(5,6) * t209 + t242 * t206;
t126 = -Icges(5,5) * t209 + t246 * t206;
t234 = t124 * t186 - t126 * t187;
t189 = qJ(5) + t192;
t182 = sin(t189);
t183 = cos(t189);
t237 = Icges(6,5) * t183 - Icges(6,6) * t182;
t108 = -Icges(6,3) * t209 + t237 * t206;
t109 = Icges(6,3) * t206 + t237 * t209;
t277 = Icges(6,4) * t183;
t241 = -Icges(6,2) * t182 + t277;
t111 = Icges(6,6) * t206 + t241 * t209;
t278 = Icges(6,4) * t182;
t245 = Icges(6,1) * t183 - t278;
t113 = Icges(6,5) * t206 + t245 * t209;
t235 = -t111 * t182 + t113 * t183;
t110 = -Icges(6,6) * t209 + t241 * t206;
t112 = -Icges(6,5) * t209 + t245 * t206;
t236 = t110 * t182 - t112 * t183;
t207 = cos(qJ(6));
t270 = t207 * t209;
t204 = sin(qJ(6));
t273 = t204 * t206;
t145 = -t183 * t273 - t270;
t271 = t206 * t207;
t272 = t204 * t209;
t146 = t183 * t271 - t272;
t276 = t182 * t206;
t73 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t276;
t75 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t276;
t77 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t276;
t18 = t145 * t75 + t146 * t77 + t73 * t276;
t147 = -t183 * t272 + t271;
t148 = t183 * t270 + t273;
t275 = t182 * t209;
t74 = Icges(7,5) * t148 + Icges(7,6) * t147 + Icges(7,3) * t275;
t76 = Icges(7,4) * t148 + Icges(7,2) * t147 + Icges(7,6) * t275;
t78 = Icges(7,1) * t148 + Icges(7,4) * t147 + Icges(7,5) * t275;
t19 = t145 * t76 + t146 * t78 + t74 * t276;
t8 = -t18 * t209 + t19 * t206;
t302 = -t202 * t108 - (t235 * t206 + (-t109 + t236) * t209) * t206 - t8;
t262 = -t202 * t122 - (t233 * t206 + (-t123 + t234) * t209) * t206 + t302;
t308 = -t202 * t130 - (t231 * t206 + (-t131 + t232) * t209) * t206 + t262;
t201 = t206 ^ 2;
t307 = t206 * pkin(7);
t274 = t183 * t209;
t80 = t148 * rSges(7,1) + t147 * rSges(7,2) + rSges(7,3) * t275;
t306 = pkin(5) * t274 + pkin(11) * t275 + t80;
t251 = rSges(5,1) * t187 - rSges(5,2) * t186;
t252 = rSges(4,1) * t191 - rSges(4,2) * t190;
t210 = -pkin(8) - pkin(7);
t305 = t206 / 0.2e1;
t304 = -t209 / 0.2e1;
t20 = t147 * t75 + t148 * t77 + t73 * t275;
t21 = t147 * t76 + t148 * t78 + t74 * t275;
t9 = -t20 * t209 + t206 * t21;
t303 = (t201 * t109 + t9 + (t236 * t209 + (-t108 + t235) * t206) * t209) * t206;
t205 = sin(qJ(2));
t301 = pkin(2) * t205;
t300 = pkin(3) * t190;
t299 = pkin(4) * t186;
t298 = pkin(5) * t183;
t208 = cos(qJ(2));
t188 = t208 * pkin(2) + pkin(1);
t167 = pkin(3) * t191 + t188;
t149 = pkin(4) * t187 + t167;
t136 = t209 * t149;
t161 = t209 * t167;
t297 = t209 * (t136 - t161) + (t149 - t167) * t201;
t180 = t209 * t188;
t296 = t209 * (t161 - t180) + (t167 - t188) * t201;
t295 = rSges(3,1) * t208;
t292 = rSges(3,2) * t205;
t97 = -Icges(7,6) * t183 + (Icges(7,4) * t207 - Icges(7,2) * t204) * t182;
t289 = t204 * t97;
t288 = t209 * rSges(3,3);
t26 = -t183 * t73 + (-t204 * t75 + t207 * t77) * t182;
t287 = t26 * t209;
t27 = -t183 * t74 + (-t204 * t76 + t207 * t78) * t182;
t286 = t27 * t206;
t99 = -rSges(7,3) * t183 + (rSges(7,1) * t207 - rSges(7,2) * t204) * t182;
t285 = -pkin(5) * t182 + pkin(11) * t183 - t99;
t284 = Icges(3,4) * t205;
t283 = Icges(3,4) * t208;
t222 = rSges(6,1) * t274 - rSges(6,2) * t275 + t206 * rSges(6,3);
t250 = rSges(6,1) * t183 - rSges(6,2) * t182;
t61 = t206 * (-rSges(6,3) * t209 + t250 * t206) + t209 * t222;
t223 = t206 * rSges(5,3) + t251 * t209;
t64 = t206 * (-rSges(5,3) * t209 + t251 * t206) + t209 * t223;
t199 = t209 * pkin(7);
t269 = t206 * (t199 + (-pkin(1) + t188) * t206) + t209 * (-t209 * pkin(1) + t180 - t307);
t224 = t206 * rSges(4,3) + t252 * t209;
t81 = t206 * (-t209 * rSges(4,3) + t252 * t206) + t209 * t224;
t267 = t206 * rSges(3,3) + t209 * t295;
t264 = t201 + t202;
t200 = -pkin(9) + t210;
t263 = t206 * (t201 * t123 + (t234 * t209 + (-t122 + t233) * t206) * t209) + t303;
t261 = t206 * (t201 * t131 + (t232 * t209 + (-t130 + t231) * t206) * t209) + t263;
t166 = rSges(4,1) * t190 + rSges(4,2) * t191;
t260 = -t166 - t301;
t160 = rSges(5,1) * t186 + rSges(5,2) * t187;
t259 = -t160 - t300;
t154 = rSges(6,1) * t182 + rSges(6,2) * t183;
t258 = -t154 - t299;
t197 = -pkin(10) + t200;
t257 = -t197 * t206 + t136;
t249 = -rSges(7,1) * t146 - rSges(7,2) * t145;
t79 = rSges(7,3) * t276 - t249;
t32 = t206 * t79 + t201 * (pkin(11) * t182 + t298) + t306 * t209;
t31 = t61 + t297;
t37 = t64 + t296;
t96 = -Icges(7,3) * t183 + (Icges(7,5) * t207 - Icges(7,6) * t204) * t182;
t98 = -Icges(7,5) * t183 + (Icges(7,1) * t207 - Icges(7,4) * t204) * t182;
t35 = t145 * t97 + t146 * t98 + t96 * t276;
t3 = -t183 * t35 + (t18 * t206 + t19 * t209) * t182;
t36 = t147 * t97 + t148 * t98 + t96 * t275;
t4 = -t183 * t36 + (t20 * t206 + t209 * t21) * t182;
t256 = t3 * t304 + t4 * t305 - t183 * (t286 - t287) / 0.2e1 + t8 * t276 / 0.2e1 + t9 * t275 / 0.2e1;
t255 = t285 - t299;
t254 = -t299 - t300;
t253 = -t292 + t295;
t248 = Icges(3,1) * t208 - t284;
t244 = -Icges(3,2) * t205 + t283;
t240 = Icges(3,5) * t208 - Icges(3,6) * t205;
t152 = Icges(6,2) * t183 + t278;
t153 = Icges(6,1) * t182 + t277;
t228 = -t152 * t182 + t153 * t183;
t158 = Icges(5,2) * t187 + t280;
t159 = Icges(5,1) * t186 + t279;
t227 = -t158 * t186 + t159 * t187;
t164 = Icges(4,2) * t191 + t282;
t165 = Icges(4,1) * t190 + t281;
t226 = -t164 * t190 + t165 * t191;
t225 = t302 * t209 + t303;
t14 = t32 + t297;
t15 = t31 + t296;
t221 = t259 - t301;
t220 = -t154 + t254;
t219 = t254 + t285;
t151 = Icges(6,5) * t182 + Icges(6,6) * t183;
t218 = t286 / 0.2e1 - t287 / 0.2e1 + (t111 * t183 + t113 * t182 + t151 * t206 + t209 * t228 + t36) * t305 + (t110 * t183 + t112 * t182 - t151 * t209 + t206 * t228 + t35) * t304;
t217 = t254 - t301;
t216 = t262 * t209 + t263;
t12 = t14 + t296;
t215 = -t154 + t217;
t214 = t209 * t308 + t261;
t213 = t217 + t285;
t157 = Icges(5,5) * t186 + Icges(5,6) * t187;
t212 = t218 + (t125 * t187 + t127 * t186 + t157 * t206 + t209 * t227) * t305 + (t124 * t187 + t126 * t186 - t157 * t209 + t206 * t227) * t304;
t163 = Icges(4,5) * t190 + Icges(4,6) * t191;
t211 = t212 + (t133 * t191 + t135 * t190 + t163 * t206 + t209 * t226) * t305 + (t132 * t191 + t134 * t190 - t163 * t209 + t206 * t226) * t304;
t179 = rSges(2,1) * t209 - rSges(2,2) * t206;
t178 = -rSges(2,1) * t206 - rSges(2,2) * t209;
t177 = rSges(3,1) * t205 + rSges(3,2) * t208;
t140 = Icges(3,3) * t206 + t240 * t209;
t139 = -Icges(3,3) * t209 + t240 * t206;
t129 = t260 * t209;
t128 = t260 * t206;
t115 = t307 + (pkin(1) - t292) * t209 + t267;
t114 = t288 + t199 + (-pkin(1) - t253) * t206;
t107 = t259 * t209;
t106 = t259 * t206;
t101 = t258 * t209;
t100 = t258 * t206;
t95 = t221 * t209;
t94 = t221 * t206;
t93 = -t206 * t210 + t180 + t224;
t92 = (rSges(4,3) - t210) * t209 + (-t188 - t252) * t206;
t89 = t220 * t209;
t88 = t220 * t206;
t87 = t182 * t207 * t98;
t86 = t209 * (-t209 * t292 + t267) + (t253 * t206 - t288) * t206;
t85 = -t206 * t200 + t161 + t223;
t84 = (rSges(5,3) - t200) * t209 + (-t167 - t251) * t206;
t83 = t215 * t209;
t82 = t215 * t206;
t70 = t285 * t209;
t69 = t285 * t206;
t68 = t222 + t257;
t67 = (rSges(6,3) - t197) * t209 + (-t149 - t250) * t206;
t58 = t255 * t209;
t57 = t255 * t206;
t52 = t219 * t209;
t51 = t219 * t206;
t48 = t213 * t209;
t47 = t213 * t206;
t44 = -t183 * t80 - t99 * t275;
t43 = t183 * t79 + t99 * t276;
t42 = t257 + t306;
t41 = -t197 * t209 + (-t298 - t149 + (-rSges(7,3) - pkin(11)) * t182) * t206 + t249;
t40 = t81 + t269;
t39 = -t182 * t289 - t183 * t96 + t87;
t38 = (-t206 * t80 + t209 * t79) * t182;
t28 = t37 + t269;
t13 = t15 + t269;
t10 = t12 + t269;
t1 = [t187 * t158 + t186 * t159 + t191 * t164 + t190 * t165 + t208 * (Icges(3,2) * t208 + t284) + t205 * (Icges(3,1) * t205 + t283) + Icges(2,3) + t87 + (t152 - t96) * t183 + (t153 - t289) * t182 + m(7) * (t41 ^ 2 + t42 ^ 2) + m(6) * (t67 ^ 2 + t68 ^ 2) + m(5) * (t84 ^ 2 + t85 ^ 2) + m(4) * (t92 ^ 2 + t93 ^ 2) + m(3) * (t114 ^ 2 + t115 ^ 2) + m(2) * (t178 ^ 2 + t179 ^ 2); m(3) * (-t114 * t209 - t115 * t206) * t177 + ((-Icges(3,6) * t209 + t244 * t206) * t208 + (-Icges(3,5) * t209 + t248 * t206) * t205) * t304 + ((Icges(3,6) * t206 + t244 * t209) * t208 + (Icges(3,5) * t206 + t248 * t209) * t205) * t305 + t211 + (t202 / 0.2e1 + t201 / 0.2e1) * (Icges(3,5) * t205 + Icges(3,6) * t208) + m(7) * (t41 * t48 + t42 * t47) + m(6) * (t67 * t83 + t68 * t82) + m(5) * (t84 * t95 + t85 * t94) + m(4) * (t128 * t93 + t129 * t92); m(7) * (t10 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t13 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(5) * (t28 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(4) * (t128 ^ 2 + t129 ^ 2 + t40 ^ 2) + m(3) * (t264 * t177 ^ 2 + t86 ^ 2) + t206 * t201 * t140 + t261 + (-t202 * t139 + (-t206 * t139 + t209 * t140) * t206 + t308) * t209; m(7) * (t41 * t52 + t42 * t51) + m(6) * (t67 * t89 + t68 * t88) + m(5) * (t106 * t85 + t107 * t84) + m(4) * (-t206 * t93 - t209 * t92) * t166 + t211; m(7) * (t10 * t12 + t47 * t51 + t48 * t52) + m(6) * (t15 * t13 + t82 * t88 + t83 * t89) + m(5) * (t106 * t94 + t107 * t95 + t37 * t28) + m(4) * (t40 * t81 + (-t128 * t206 - t129 * t209) * t166) + t214; m(7) * (t12 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(6) * (t15 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(5) * (t106 ^ 2 + t107 ^ 2 + t37 ^ 2) + m(4) * (t264 * t166 ^ 2 + t81 ^ 2) + t214; m(7) * (t41 * t58 + t42 * t57) + m(6) * (t100 * t68 + t101 * t67) + t212 + m(5) * (-t206 * t85 - t209 * t84) * t160; m(7) * (t10 * t14 + t47 * t57 + t48 * t58) + m(6) * (t100 * t82 + t101 * t83 + t31 * t13) + m(5) * (t64 * t28 + (-t206 * t94 - t209 * t95) * t160) + t216; m(7) * (t12 * t14 + t51 * t57 + t52 * t58) + m(6) * (t100 * t88 + t101 * t89 + t31 * t15) + m(5) * (t64 * t37 + (-t106 * t206 - t107 * t209) * t160) + t216; m(7) * (t14 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(6) * (t100 ^ 2 + t101 ^ 2 + t31 ^ 2) + m(5) * (t264 * t160 ^ 2 + t64 ^ 2) + t216; m(7) * (t41 * t70 + t42 * t69) + m(6) * (-t206 * t68 - t209 * t67) * t154 + t218; m(7) * (t10 * t32 + t47 * t69 + t48 * t70) + m(6) * (t61 * t13 + (-t206 * t82 - t209 * t83) * t154) + t225; m(7) * (t12 * t32 + t51 * t69 + t52 * t70) + m(6) * (t61 * t15 + (-t206 * t88 - t209 * t89) * t154) + t225; m(7) * (t14 * t32 + t57 * t69 + t58 * t70) + m(6) * (t61 * t31 + (-t100 * t206 - t101 * t209) * t154) + t225; m(6) * (t264 * t154 ^ 2 + t61 ^ 2) + m(7) * (t32 ^ 2 + t69 ^ 2 + t70 ^ 2) + t225; m(7) * (t41 * t43 + t42 * t44) - t39 * t183 + ((t27 / 0.2e1 + t36 / 0.2e1) * t209 + (t35 / 0.2e1 + t26 / 0.2e1) * t206) * t182; m(7) * (t10 * t38 + t43 * t48 + t44 * t47) + t256; m(7) * (t12 * t38 + t43 * t52 + t44 * t51) + t256; m(7) * (t14 * t38 + t43 * t58 + t44 * t57) + t256; m(7) * (t32 * t38 + t43 * t70 + t44 * t69) + t256; t183 ^ 2 * t39 + m(7) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + (t209 * t4 + t206 * t3 - t183 * (t26 * t206 + t27 * t209)) * t182;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
