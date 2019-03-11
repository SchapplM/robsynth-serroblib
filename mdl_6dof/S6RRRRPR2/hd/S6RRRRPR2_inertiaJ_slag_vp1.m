% Calculate joint inertia matrix for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:25
% EndTime: 2019-03-09 21:57:31
% DurationCPUTime: 3.18s
% Computational Cost: add. (13035->461), mult. (10303->668), div. (0->0), fcn. (10831->12), ass. (0->227)
t210 = sin(qJ(1));
t203 = t210 ^ 2;
t306 = t210 * pkin(7);
t207 = cos(pkin(11));
t187 = pkin(5) * t207 + pkin(4);
t208 = -pkin(10) - qJ(5);
t206 = sin(pkin(11));
t271 = t206 * t210;
t205 = qJ(2) + qJ(3);
t195 = qJ(4) + t205;
t189 = cos(t195);
t212 = cos(qJ(1));
t277 = t189 * t212;
t188 = sin(t195);
t278 = t188 * t212;
t201 = pkin(11) + qJ(6);
t192 = cos(t201);
t273 = t192 * t210;
t191 = sin(t201);
t274 = t191 * t212;
t142 = -t189 * t274 + t273;
t272 = t192 * t212;
t275 = t191 * t210;
t143 = t189 * t272 + t275;
t77 = t143 * rSges(7,1) + t142 * rSges(7,2) + rSges(7,3) * t278;
t305 = pkin(5) * t271 + t187 * t277 - t208 * t278 + t77;
t193 = sin(t205);
t194 = cos(t205);
t246 = rSges(4,1) * t194 - rSges(4,2) * t193;
t235 = Icges(4,5) * t194 - Icges(4,6) * t193;
t134 = -Icges(4,3) * t212 + t235 * t210;
t135 = Icges(4,3) * t210 + t235 * t212;
t204 = t212 ^ 2;
t282 = Icges(4,4) * t194;
t238 = -Icges(4,2) * t193 + t282;
t137 = Icges(4,6) * t210 + t238 * t212;
t283 = Icges(4,4) * t193;
t241 = Icges(4,1) * t194 - t283;
t139 = Icges(4,5) * t210 + t241 * t212;
t230 = -t137 * t193 + t139 * t194;
t136 = -Icges(4,6) * t212 + t238 * t210;
t138 = -Icges(4,5) * t212 + t241 * t210;
t231 = t136 * t193 - t138 * t194;
t234 = Icges(5,5) * t189 - Icges(5,6) * t188;
t121 = -Icges(5,3) * t212 + t234 * t210;
t122 = Icges(5,3) * t210 + t234 * t212;
t268 = t207 * t212;
t152 = -t189 * t271 - t268;
t269 = t207 * t210;
t270 = t206 * t212;
t153 = t189 * t269 - t270;
t280 = Icges(5,4) * t189;
t237 = -Icges(5,2) * t188 + t280;
t124 = Icges(5,6) * t210 + t237 * t212;
t281 = Icges(5,4) * t188;
t240 = Icges(5,1) * t189 - t281;
t126 = Icges(5,5) * t210 + t240 * t212;
t232 = -t124 * t188 + t126 * t189;
t123 = -Icges(5,6) * t212 + t237 * t210;
t125 = -Icges(5,5) * t212 + t240 * t210;
t233 = t123 * t188 - t125 * t189;
t279 = t188 * t210;
t140 = -t189 * t275 - t272;
t141 = t189 * t273 - t274;
t70 = Icges(7,5) * t141 + Icges(7,6) * t140 + Icges(7,3) * t279;
t72 = Icges(7,4) * t141 + Icges(7,2) * t140 + Icges(7,6) * t279;
t74 = Icges(7,1) * t141 + Icges(7,4) * t140 + Icges(7,5) * t279;
t18 = t140 * t72 + t141 * t74 + t70 * t279;
t71 = Icges(7,5) * t143 + Icges(7,6) * t142 + Icges(7,3) * t278;
t73 = Icges(7,4) * t143 + Icges(7,2) * t142 + Icges(7,6) * t278;
t75 = Icges(7,1) * t143 + Icges(7,4) * t142 + Icges(7,5) * t278;
t19 = t140 * t73 + t141 * t75 + t71 * t279;
t8 = -t18 * t212 + t19 * t210;
t84 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t279;
t154 = -t189 * t270 + t269;
t155 = t189 * t268 + t271;
t85 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t278;
t86 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t279;
t87 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t278;
t88 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t279;
t89 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t278;
t258 = -t8 - t204 * t121 + (t152 * t86 + t153 * t88 + t84 * t279) * t212 + (-(-t122 + t233) * t212 - t152 * t87 - t153 * t89 - t85 * t279 - t232 * t210) * t210;
t304 = -t204 * t134 - (t230 * t210 + (-t135 + t231) * t212) * t210 + t258;
t106 = -rSges(7,3) * t189 + (rSges(7,1) * t192 - rSges(7,2) * t191) * t188;
t267 = qJ(5) + t208;
t294 = -pkin(4) + t187;
t303 = -t294 * t188 - t267 * t189 - t106;
t302 = t189 ^ 2;
t301 = m(6) / 0.2e1;
t300 = m(7) / 0.2e1;
t213 = -pkin(8) - pkin(7);
t299 = t210 / 0.2e1;
t298 = -t212 / 0.2e1;
t209 = sin(qJ(2));
t297 = pkin(2) * t209;
t296 = pkin(3) * t193;
t295 = pkin(4) * t189;
t211 = cos(qJ(2));
t190 = t211 * pkin(2) + pkin(1);
t169 = pkin(3) * t194 + t190;
t163 = t212 * t169;
t182 = t212 * t190;
t293 = t212 * (t163 - t182) + (t169 - t190) * t203;
t292 = rSges(3,1) * t211;
t290 = rSges(3,2) * t209;
t288 = t212 * rSges(3,3);
t24 = -t189 * t70 + (-t191 * t72 + t192 * t74) * t188;
t287 = t24 * t212;
t25 = -t189 * t71 + (-t191 * t73 + t192 * t75) * t188;
t286 = t25 * t210;
t285 = Icges(3,4) * t209;
t284 = Icges(3,4) * t211;
t104 = -Icges(7,6) * t189 + (Icges(7,4) * t192 - Icges(7,2) * t191) * t188;
t276 = t191 * t104;
t110 = -rSges(6,3) * t189 + (rSges(6,1) * t207 - rSges(6,2) * t206) * t188;
t160 = pkin(4) * t188 - qJ(5) * t189;
t266 = -t110 - t160;
t223 = rSges(5,1) * t277 - rSges(5,2) * t278 + t210 * rSges(5,3);
t245 = rSges(5,1) * t189 - rSges(5,2) * t188;
t80 = t210 * (-rSges(5,3) * t212 + t245 * t210) + t212 * t223;
t200 = t212 * pkin(7);
t265 = t210 * (t200 + (-pkin(1) + t190) * t210) + t212 * (-t212 * pkin(1) + t182 - t306);
t224 = t210 * rSges(4,3) + t246 * t212;
t83 = t210 * (-t212 * rSges(4,3) + t246 * t210) + t212 * t224;
t263 = pkin(4) * t277 + qJ(5) * t278;
t264 = t203 * (qJ(5) * t188 + t295) + t212 * t263;
t262 = t210 * rSges(3,3) + t212 * t292;
t260 = t203 + t204;
t20 = t142 * t72 + t143 * t74 + t70 * t278;
t21 = t142 * t73 + t143 * t75 + t71 * t278;
t9 = -t20 * t212 + t21 * t210;
t259 = (t9 + t203 * t122 + (t154 * t87 + t155 * t89 + t85 * t278) * t210 + (-t154 * t86 - t155 * t88 - t84 * t278 + (-t121 + t232) * t210 + t233 * t212) * t212) * t210;
t257 = -t160 + t303;
t256 = t155 * rSges(6,1) + t154 * rSges(6,2) + rSges(6,3) * t278;
t255 = t210 * (t203 * t135 + (t231 * t212 + (-t134 + t230) * t210) * t212) + t259;
t168 = rSges(4,1) * t193 + rSges(4,2) * t194;
t254 = -t168 - t297;
t253 = -t160 - t296;
t161 = rSges(5,1) * t188 + rSges(5,2) * t189;
t252 = -t161 - t296;
t202 = -pkin(9) + t213;
t251 = -t210 * t202 + t163;
t244 = -rSges(6,1) * t153 - rSges(6,2) * t152;
t35 = t210 * (rSges(6,3) * t279 - t244) + t212 * t256 + t264;
t40 = t80 + t293;
t103 = -Icges(7,3) * t189 + (Icges(7,5) * t192 - Icges(7,6) * t191) * t188;
t105 = -Icges(7,5) * t189 + (Icges(7,1) * t192 - Icges(7,4) * t191) * t188;
t36 = t103 * t279 + t104 * t140 + t105 * t141;
t3 = -t189 * t36 + (t18 * t210 + t19 * t212) * t188;
t37 = t103 * t278 + t104 * t142 + t105 * t143;
t4 = -t189 * t37 + (t20 * t210 + t21 * t212) * t188;
t250 = t3 * t298 + t4 * t299 - t189 * (t286 - t287) / 0.2e1 + t8 * t279 / 0.2e1 + t9 * t278 / 0.2e1;
t249 = -t110 + t253;
t248 = -t296 - t297;
t247 = -t290 + t292;
t243 = -rSges(7,1) * t141 - rSges(7,2) * t140;
t242 = Icges(3,1) * t211 - t285;
t239 = -Icges(3,2) * t209 + t284;
t236 = Icges(3,5) * t211 - Icges(3,6) * t209;
t158 = Icges(5,2) * t189 + t281;
t159 = Icges(5,1) * t188 + t280;
t227 = -t158 * t188 + t159 * t189;
t166 = Icges(4,2) * t194 + t283;
t167 = Icges(4,1) * t193 + t282;
t226 = -t166 * t193 + t167 * t194;
t225 = t253 + t303;
t76 = rSges(7,3) * t279 - t243;
t16 = t264 + (-t263 + t305) * t212 + (-pkin(5) * t270 + t76 + (-t267 * t188 + t294 * t189) * t210) * t210;
t17 = t35 + t293;
t222 = -t160 + t248;
t221 = -t161 + t248;
t219 = -t110 + t222;
t218 = t258 * t212 + t259;
t14 = t16 + t293;
t217 = t222 + t303;
t216 = t212 * t304 + t255;
t107 = -Icges(6,3) * t189 + (Icges(6,5) * t207 - Icges(6,6) * t206) * t188;
t108 = -Icges(6,6) * t189 + (Icges(6,4) * t207 - Icges(6,2) * t206) * t188;
t109 = -Icges(6,5) * t189 + (Icges(6,1) * t207 - Icges(6,4) * t206) * t188;
t157 = Icges(5,5) * t188 + Icges(5,6) * t189;
t215 = t286 / 0.2e1 - t287 / 0.2e1 + (t107 * t278 + t108 * t154 + t109 * t155 + t157 * t210 + t227 * t212 + t37 + (t124 - t85) * t189 + (-t206 * t87 + t207 * t89 + t126) * t188) * t299 + (t107 * t279 + t108 * t152 + t109 * t153 - t157 * t212 + t227 * t210 + t36 + (-t84 + t123) * t189 + (-t206 * t86 + t207 * t88 + t125) * t188) * t298;
t165 = Icges(4,5) * t193 + Icges(4,6) * t194;
t214 = t215 + (t137 * t194 + t139 * t193 + t165 * t210 + t226 * t212) * t299 + (t136 * t194 + t138 * t193 - t165 * t212 + t226 * t210) * t298;
t181 = rSges(2,1) * t212 - rSges(2,2) * t210;
t180 = -rSges(2,1) * t210 - rSges(2,2) * t212;
t179 = rSges(3,1) * t209 + rSges(3,2) * t211;
t147 = Icges(3,3) * t210 + t236 * t212;
t146 = -Icges(3,3) * t212 + t236 * t210;
t133 = t254 * t212;
t132 = t254 * t210;
t116 = t306 + (pkin(1) - t290) * t212 + t262;
t115 = t288 + t200 + (-pkin(1) - t247) * t210;
t114 = t252 * t212;
t113 = t252 * t210;
t102 = t221 * t212;
t101 = t221 * t210;
t99 = -t210 * t213 + t182 + t224;
t98 = (rSges(4,3) - t213) * t212 + (-t190 - t246) * t210;
t95 = t188 * t192 * t105;
t94 = t212 * (-t212 * t290 + t262) + (t247 * t210 - t288) * t210;
t93 = t223 + t251;
t92 = (rSges(5,3) - t202) * t212 + (-t169 - t245) * t210;
t91 = t266 * t212;
t90 = t266 * t210;
t69 = t249 * t212;
t68 = t249 * t210;
t61 = t219 * t212;
t60 = t219 * t210;
t55 = t251 + t256 + t263;
t54 = -t212 * t202 + (-t295 - t169 + (-rSges(6,3) - qJ(5)) * t188) * t210 + t244;
t53 = t257 * t212;
t52 = t257 * t210;
t51 = t225 * t212;
t50 = t225 * t210;
t49 = t83 + t265;
t48 = -t106 * t278 - t189 * t77;
t47 = t106 * t279 + t189 * t76;
t46 = t217 * t212;
t45 = t217 * t210;
t44 = t251 + t305;
t43 = (pkin(5) * t206 - t202) * t212 + (-t187 * t189 - t169 + (-rSges(7,3) + t208) * t188) * t210 + t243;
t42 = -t189 * t103 - t188 * t276 + t95;
t41 = (-t210 * t77 + t212 * t76) * t188;
t28 = t40 + t265;
t15 = t17 + t265;
t11 = t14 + t265;
t1 = [t194 * t166 + t193 * t167 + t211 * (Icges(3,2) * t211 + t285) + t209 * (Icges(3,1) * t209 + t284) + Icges(2,3) + t95 + (-t103 + t158 - t107) * t189 + (-t206 * t108 + t207 * t109 + t159 - t276) * t188 + m(7) * (t43 ^ 2 + t44 ^ 2) + m(6) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t92 ^ 2 + t93 ^ 2) + m(4) * (t98 ^ 2 + t99 ^ 2) + m(3) * (t115 ^ 2 + t116 ^ 2) + m(2) * (t180 ^ 2 + t181 ^ 2); t214 + (t204 / 0.2e1 + t203 / 0.2e1) * (Icges(3,5) * t209 + Icges(3,6) * t211) + m(7) * (t43 * t46 + t44 * t45) + m(6) * (t54 * t61 + t55 * t60) + m(5) * (t101 * t93 + t102 * t92) + m(4) * (t132 * t99 + t133 * t98) + m(3) * (-t115 * t212 - t116 * t210) * t179 + ((Icges(3,6) * t210 + t239 * t212) * t211 + (Icges(3,5) * t210 + t242 * t212) * t209) * t299 + ((-Icges(3,6) * t212 + t239 * t210) * t211 + (-Icges(3,5) * t212 + t242 * t210) * t209) * t298; m(7) * (t11 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t101 ^ 2 + t102 ^ 2 + t28 ^ 2) + m(4) * (t132 ^ 2 + t133 ^ 2 + t49 ^ 2) + t210 * t203 * t147 + m(3) * (t260 * t179 ^ 2 + t94 ^ 2) + t255 + (-t204 * t146 + (-t210 * t146 + t212 * t147) * t210 + t304) * t212; t214 + m(7) * (t43 * t51 + t44 * t50) + m(6) * (t54 * t69 + t55 * t68) + m(5) * (t113 * t93 + t114 * t92) + m(4) * (-t210 * t99 - t212 * t98) * t168; m(7) * (t11 * t14 + t45 * t50 + t46 * t51) + m(6) * (t15 * t17 + t60 * t68 + t61 * t69) + m(5) * (t101 * t113 + t102 * t114 + t40 * t28) + m(4) * (t49 * t83 + (-t132 * t210 - t133 * t212) * t168) + t216; m(7) * (t14 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(6) * (t17 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(5) * (t113 ^ 2 + t114 ^ 2 + t40 ^ 2) + m(4) * (t260 * t168 ^ 2 + t83 ^ 2) + t216; t215 + m(7) * (t43 * t53 + t44 * t52) + m(6) * (t54 * t91 + t55 * t90) + m(5) * (-t210 * t93 - t212 * t92) * t161; m(7) * (t11 * t16 + t45 * t52 + t46 * t53) + m(6) * (t35 * t15 + t60 * t90 + t61 * t91) + m(5) * (t28 * t80 + (-t101 * t210 - t102 * t212) * t161) + t218; m(7) * (t14 * t16 + t50 * t52 + t51 * t53) + m(6) * (t35 * t17 + t68 * t90 + t69 * t91) + m(5) * (t40 * t80 + (-t113 * t210 - t114 * t212) * t161) + t218; m(7) * (t16 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t35 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(5) * (t260 * t161 ^ 2 + t80 ^ 2) + t218; 0.2e1 * ((t210 * t44 + t212 * t43) * t300 + (t210 * t55 + t212 * t54) * t301) * t188; m(7) * (-t11 * t189 + (t210 * t45 + t212 * t46) * t188) + m(6) * (-t15 * t189 + (t210 * t60 + t212 * t61) * t188); m(7) * (-t14 * t189 + (t210 * t50 + t212 * t51) * t188) + m(6) * (-t17 * t189 + (t210 * t68 + t212 * t69) * t188); m(7) * (-t16 * t189 + (t210 * t52 + t212 * t53) * t188) + m(6) * (-t189 * t35 + (t210 * t90 + t212 * t91) * t188); 0.2e1 * (t301 + t300) * (t260 * t188 ^ 2 + t302); m(7) * (t43 * t47 + t44 * t48) - t42 * t189 + ((t37 / 0.2e1 + t25 / 0.2e1) * t212 + (t36 / 0.2e1 + t24 / 0.2e1) * t210) * t188; m(7) * (t11 * t41 + t45 * t48 + t46 * t47) + t250; m(7) * (t14 * t41 + t47 * t51 + t48 * t50) + t250; m(7) * (t16 * t41 + t47 * t53 + t48 * t52) + t250; m(7) * (-t189 * t41 + (t210 * t48 + t212 * t47) * t188); t302 * t42 + m(7) * (t41 ^ 2 + t47 ^ 2 + t48 ^ 2) + (t212 * t4 + t210 * t3 - t189 * (t210 * t24 + t212 * t25)) * t188;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
