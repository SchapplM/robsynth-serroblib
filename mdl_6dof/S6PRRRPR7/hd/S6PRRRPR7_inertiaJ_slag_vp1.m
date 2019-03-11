% Calculate joint inertia matrix for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:39:05
% EndTime: 2019-03-08 23:39:23
% DurationCPUTime: 8.04s
% Computational Cost: add. (53649->639), mult. (140458->908), div. (0->0), fcn. (185237->16), ass. (0->293)
t319 = m(6) + m(7);
t247 = sin(pkin(12));
t250 = cos(pkin(12));
t255 = sin(qJ(2));
t251 = cos(pkin(6));
t256 = cos(qJ(2));
t300 = t251 * t256;
t235 = -t247 * t255 + t250 * t300;
t301 = t251 * t255;
t236 = t247 * t256 + t250 * t301;
t254 = sin(qJ(3));
t248 = sin(pkin(6));
t308 = sin(pkin(7));
t277 = t248 * t308;
t309 = cos(pkin(7));
t312 = cos(qJ(3));
t204 = t236 * t312 + (t235 * t309 - t250 * t277) * t254;
t278 = t248 * t309;
t224 = -t235 * t308 - t250 * t278;
t253 = sin(qJ(4));
t311 = cos(qJ(4));
t190 = t204 * t253 - t224 * t311;
t237 = -t247 * t300 - t250 * t255;
t238 = -t247 * t301 + t250 * t256;
t206 = t238 * t312 + (t237 * t309 + t247 * t277) * t254;
t225 = -t237 * t308 + t247 * t278;
t192 = t206 * t253 - t225 * t311;
t223 = t251 * t308 * t254 + (t254 * t256 * t309 + t255 * t312) * t248;
t234 = t251 * t309 - t256 * t277;
t207 = t223 * t253 - t234 * t311;
t191 = t204 * t311 + t224 * t253;
t261 = t312 * t308;
t259 = t248 * t261;
t262 = t309 * t312;
t203 = -t235 * t262 + t236 * t254 + t250 * t259;
t245 = pkin(13) + qJ(6);
t243 = sin(t245);
t244 = cos(t245);
t155 = -t191 * t243 + t203 * t244;
t156 = t191 * t244 + t203 * t243;
t103 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t190;
t105 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t190;
t107 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t190;
t40 = t103 * t190 + t105 * t155 + t107 * t156;
t193 = t206 * t311 + t225 * t253;
t205 = -t237 * t262 + t238 * t254 - t247 * t259;
t157 = -t193 * t243 + t205 * t244;
t158 = t193 * t244 + t205 * t243;
t104 = Icges(7,5) * t158 + Icges(7,6) * t157 + Icges(7,3) * t192;
t106 = Icges(7,4) * t158 + Icges(7,2) * t157 + Icges(7,6) * t192;
t108 = Icges(7,1) * t158 + Icges(7,4) * t157 + Icges(7,5) * t192;
t41 = t104 * t190 + t106 * t155 + t108 * t156;
t208 = t223 * t311 + t234 * t253;
t302 = t248 * t255;
t222 = -t248 * t256 * t262 - t251 * t261 + t254 * t302;
t185 = -t208 * t243 + t222 * t244;
t186 = t208 * t244 + t222 * t243;
t128 = Icges(7,5) * t186 + Icges(7,6) * t185 + Icges(7,3) * t207;
t129 = Icges(7,4) * t186 + Icges(7,2) * t185 + Icges(7,6) * t207;
t130 = Icges(7,1) * t186 + Icges(7,4) * t185 + Icges(7,5) * t207;
t65 = t128 * t190 + t129 * t155 + t130 * t156;
t1 = t190 * t40 + t192 * t41 + t207 * t65;
t318 = t1 / 0.2e1;
t42 = t103 * t192 + t105 * t157 + t107 * t158;
t43 = t104 * t192 + t106 * t157 + t108 * t158;
t66 = t128 * t192 + t129 * t157 + t130 * t158;
t2 = t190 * t42 + t192 * t43 + t207 * t66;
t317 = t2 / 0.2e1;
t54 = t103 * t207 + t105 * t185 + t107 * t186;
t55 = t104 * t207 + t106 * t185 + t108 * t186;
t72 = t128 * t207 + t129 * t185 + t130 * t186;
t15 = t190 * t54 + t192 * t55 + t207 * t72;
t316 = t15 / 0.2e1;
t315 = t190 / 0.2e1;
t314 = t192 / 0.2e1;
t313 = t207 / 0.2e1;
t249 = cos(pkin(13));
t310 = pkin(5) * t249;
t246 = sin(pkin(13));
t307 = t203 * t246;
t306 = t205 * t246;
t305 = t222 * t246;
t304 = t247 * t248;
t303 = t248 * t250;
t109 = rSges(7,1) * t156 + rSges(7,2) * t155 + rSges(7,3) * t190;
t298 = pkin(5) * t307 + pkin(11) * t190 + t191 * t310 + t109;
t110 = rSges(7,1) * t158 + rSges(7,2) * t157 + rSges(7,3) * t192;
t297 = pkin(5) * t306 + pkin(11) * t192 + t193 * t310 + t110;
t159 = -t191 * t246 + t203 * t249;
t160 = t191 * t249 + t307;
t117 = rSges(6,1) * t160 + rSges(6,2) * t159 + rSges(6,3) * t190;
t152 = t191 * pkin(4) + t190 * qJ(5);
t296 = -t117 - t152;
t161 = -t193 * t246 + t205 * t249;
t162 = t193 * t249 + t306;
t118 = rSges(6,1) * t162 + rSges(6,2) * t161 + rSges(6,3) * t192;
t153 = t193 * pkin(4) + t192 * qJ(5);
t295 = -t118 - t153;
t131 = rSges(7,1) * t186 + rSges(7,2) * t185 + rSges(7,3) * t207;
t294 = pkin(5) * t305 + pkin(11) * t207 + t208 * t310 + t131;
t188 = -t208 * t246 + t222 * t249;
t189 = t208 * t249 + t305;
t135 = rSges(6,1) * t189 + rSges(6,2) * t188 + rSges(6,3) * t207;
t183 = t208 * pkin(4) + t207 * qJ(5);
t293 = -t135 - t183;
t142 = rSges(5,1) * t191 - rSges(5,2) * t190 + rSges(5,3) * t203;
t181 = pkin(3) * t204 + pkin(10) * t203;
t292 = -t142 - t181;
t175 = t225 * t181;
t291 = t225 * t152 + t175;
t182 = pkin(3) * t206 + pkin(10) * t205;
t177 = t234 * t182;
t290 = t234 * t153 + t177;
t174 = rSges(5,1) * t208 - rSges(5,2) * t207 + rSges(5,3) * t222;
t200 = pkin(3) * t223 + pkin(10) * t222;
t289 = -t174 - t200;
t187 = t224 * t200;
t288 = t224 * t183 + t187;
t211 = t238 * pkin(2) + pkin(9) * t225;
t209 = t251 * t211;
t287 = t251 * t182 + t209;
t210 = t236 * pkin(2) + pkin(9) * t224;
t286 = t210 * t304 + t211 * t303;
t284 = -t152 - t298;
t283 = -t153 - t297;
t282 = -t181 + t296;
t281 = -t183 - t294;
t280 = -t200 + t293;
t279 = t251 * t153 + t287;
t197 = rSges(4,1) * t223 - rSges(4,2) * t222 + rSges(4,3) * t234;
t226 = pkin(2) * t302 + pkin(9) * t234;
t276 = (-t197 - t226) * t248;
t136 = Icges(5,5) * t191 - Icges(5,6) * t190 + Icges(5,3) * t203;
t138 = Icges(5,4) * t191 - Icges(5,2) * t190 + Icges(5,6) * t203;
t140 = Icges(5,1) * t191 - Icges(5,4) * t190 + Icges(5,5) * t203;
t74 = t136 * t203 - t138 * t190 + t140 * t191;
t137 = Icges(5,5) * t193 - Icges(5,6) * t192 + Icges(5,3) * t205;
t139 = Icges(5,4) * t193 - Icges(5,2) * t192 + Icges(5,6) * t205;
t141 = Icges(5,1) * t193 - Icges(5,4) * t192 + Icges(5,5) * t205;
t75 = t137 * t203 - t139 * t190 + t141 * t191;
t169 = Icges(5,5) * t208 - Icges(5,6) * t207 + Icges(5,3) * t222;
t170 = Icges(5,4) * t208 - Icges(5,2) * t207 + Icges(5,6) * t222;
t171 = Icges(5,1) * t208 - Icges(5,4) * t207 + Icges(5,5) * t222;
t84 = t169 * t203 - t170 * t190 + t171 * t191;
t22 = t203 * t74 + t205 * t75 + t222 * t84;
t3 = t203 * t40 + t205 * t41 + t222 * t65;
t111 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t190;
t113 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t190;
t115 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t190;
t46 = t111 * t190 + t113 * t159 + t115 * t160;
t112 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t192;
t114 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t192;
t116 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t192;
t47 = t112 * t190 + t114 * t159 + t116 * t160;
t132 = Icges(6,5) * t189 + Icges(6,6) * t188 + Icges(6,3) * t207;
t133 = Icges(6,4) * t189 + Icges(6,2) * t188 + Icges(6,6) * t207;
t134 = Icges(6,1) * t189 + Icges(6,4) * t188 + Icges(6,5) * t207;
t67 = t132 * t190 + t133 * t159 + t134 * t160;
t7 = t203 * t46 + t205 * t47 + t222 * t67;
t275 = t3 / 0.2e1 + t7 / 0.2e1 + t22 / 0.2e1;
t76 = t136 * t205 - t138 * t192 + t140 * t193;
t77 = t137 * t205 - t139 * t192 + t141 * t193;
t85 = t169 * t205 - t170 * t192 + t171 * t193;
t23 = t203 * t76 + t205 * t77 + t222 * t85;
t4 = t203 * t42 + t205 * t43 + t222 * t66;
t48 = t111 * t192 + t113 * t161 + t115 * t162;
t49 = t112 * t192 + t114 * t161 + t116 * t162;
t68 = t132 * t192 + t133 * t161 + t134 * t162;
t8 = t203 * t48 + t205 * t49 + t222 * t68;
t274 = t4 / 0.2e1 + t23 / 0.2e1 + t8 / 0.2e1;
t273 = -t181 + t284;
t272 = -t200 + t281;
t271 = t181 * t304 + t182 * t303 + t286;
t11 = t224 * t46 + t225 * t47 + t234 * t67;
t24 = t224 * t74 + t225 * t75 + t234 * t84;
t5 = t224 * t40 + t225 * t41 + t234 * t65;
t270 = t5 / 0.2e1 + t11 / 0.2e1 + t24 / 0.2e1;
t12 = t224 * t48 + t225 * t49 + t234 * t68;
t25 = t224 * t76 + t225 * t77 + t234 * t85;
t6 = t224 * t42 + t225 * t43 + t234 * t66;
t269 = t6 / 0.2e1 + t12 / 0.2e1 + t25 / 0.2e1;
t13 = t67 * t251 + (t247 * t47 - t250 * t46) * t248;
t26 = t84 * t251 + (t247 * t75 - t250 * t74) * t248;
t9 = t65 * t251 + (t247 * t41 - t250 * t40) * t248;
t268 = t9 / 0.2e1 + t26 / 0.2e1 + t13 / 0.2e1;
t10 = t66 * t251 + (t247 * t43 - t250 * t42) * t248;
t14 = t68 * t251 + (t247 * t49 - t250 * t48) * t248;
t27 = t85 * t251 + (t247 * t77 - t250 * t76) * t248;
t267 = t10 / 0.2e1 + t27 / 0.2e1 + t14 / 0.2e1;
t16 = t203 * t54 + t205 * t55 + t222 * t72;
t56 = t111 * t207 + t113 * t188 + t115 * t189;
t57 = t112 * t207 + t114 * t188 + t116 * t189;
t73 = t132 * t207 + t133 * t188 + t134 * t189;
t18 = t203 * t56 + t205 * t57 + t222 * t73;
t78 = t136 * t222 - t138 * t207 + t140 * t208;
t79 = t137 * t222 - t139 * t207 + t141 * t208;
t95 = t169 * t222 - t170 * t207 + t171 * t208;
t30 = t203 * t78 + t205 * t79 + t222 * t95;
t266 = t16 / 0.2e1 + t30 / 0.2e1 + t18 / 0.2e1;
t17 = t224 * t54 + t225 * t55 + t234 * t72;
t20 = t224 * t56 + t225 * t57 + t234 * t73;
t32 = t224 * t78 + t225 * t79 + t234 * t95;
t265 = t17 / 0.2e1 + t20 / 0.2e1 + t32 / 0.2e1;
t19 = t72 * t251 + (t247 * t55 - t250 * t54) * t248;
t21 = t73 * t251 + (t247 * t57 - t250 * t56) * t248;
t33 = t95 * t251 + (t247 * t79 - t250 * t78) * t248;
t264 = t19 / 0.2e1 + t33 / 0.2e1 + t21 / 0.2e1;
t263 = (-t226 + t289) * t248;
t260 = (-t226 + t280) * t248;
t258 = t152 * t304 + t153 * t303 + t271;
t257 = (-t226 + t272) * t248;
t232 = t251 * rSges(3,3) + (rSges(3,1) * t255 + rSges(3,2) * t256) * t248;
t231 = Icges(3,5) * t251 + (Icges(3,1) * t255 + Icges(3,4) * t256) * t248;
t230 = Icges(3,6) * t251 + (Icges(3,4) * t255 + Icges(3,2) * t256) * t248;
t229 = Icges(3,3) * t251 + (Icges(3,5) * t255 + Icges(3,6) * t256) * t248;
t219 = rSges(3,1) * t238 + rSges(3,2) * t237 + rSges(3,3) * t304;
t218 = rSges(3,1) * t236 + rSges(3,2) * t235 - rSges(3,3) * t303;
t217 = Icges(3,1) * t238 + Icges(3,4) * t237 + Icges(3,5) * t304;
t216 = Icges(3,1) * t236 + Icges(3,4) * t235 - Icges(3,5) * t303;
t215 = Icges(3,4) * t238 + Icges(3,2) * t237 + Icges(3,6) * t304;
t214 = Icges(3,4) * t236 + Icges(3,2) * t235 - Icges(3,6) * t303;
t213 = Icges(3,5) * t238 + Icges(3,6) * t237 + Icges(3,3) * t304;
t212 = Icges(3,5) * t236 + Icges(3,6) * t235 - Icges(3,3) * t303;
t199 = -t218 * t251 - t232 * t303;
t198 = t219 * t251 - t232 * t304;
t196 = Icges(4,1) * t223 - Icges(4,4) * t222 + Icges(4,5) * t234;
t195 = Icges(4,4) * t223 - Icges(4,2) * t222 + Icges(4,6) * t234;
t194 = Icges(4,5) * t223 - Icges(4,6) * t222 + Icges(4,3) * t234;
t184 = (t218 * t247 + t219 * t250) * t248;
t173 = rSges(4,1) * t206 - rSges(4,2) * t205 + rSges(4,3) * t225;
t172 = rSges(4,1) * t204 - rSges(4,2) * t203 + rSges(4,3) * t224;
t168 = Icges(4,1) * t206 - Icges(4,4) * t205 + Icges(4,5) * t225;
t167 = Icges(4,1) * t204 - Icges(4,4) * t203 + Icges(4,5) * t224;
t166 = Icges(4,4) * t206 - Icges(4,2) * t205 + Icges(4,6) * t225;
t165 = Icges(4,4) * t204 - Icges(4,2) * t203 + Icges(4,6) * t224;
t164 = Icges(4,5) * t206 - Icges(4,6) * t205 + Icges(4,3) * t225;
t163 = Icges(4,5) * t204 - Icges(4,6) * t203 + Icges(4,3) * t224;
t154 = t203 * t183;
t146 = t222 * t153;
t144 = t205 * t152;
t143 = rSges(5,1) * t193 - rSges(5,2) * t192 + rSges(5,3) * t205;
t126 = t173 * t234 - t197 * t225;
t125 = -t172 * t234 + t197 * t224;
t124 = (-t172 - t210) * t251 + t250 * t276;
t123 = t251 * t173 + t247 * t276 + t209;
t122 = t194 * t234 - t195 * t222 + t196 * t223;
t121 = t172 * t225 - t173 * t224;
t120 = t194 * t225 - t195 * t205 + t196 * t206;
t119 = t194 * t224 - t195 * t203 + t196 * t204;
t102 = (t172 * t247 + t173 * t250) * t248 + t286;
t99 = t143 * t222 - t174 * t205;
t98 = -t142 * t222 + t174 * t203;
t97 = t164 * t234 - t166 * t222 + t168 * t223;
t96 = t163 * t234 - t165 * t222 + t167 * t223;
t94 = t164 * t225 - t166 * t205 + t168 * t206;
t93 = t163 * t225 - t165 * t205 + t167 * t206;
t92 = t164 * t224 - t166 * t203 + t168 * t204;
t91 = t163 * t224 - t165 * t203 + t167 * t204;
t90 = t142 * t205 - t143 * t203;
t89 = t143 * t234 + t225 * t289 + t177;
t88 = t224 * t174 + t234 * t292 + t187;
t87 = (-t210 + t292) * t251 + t250 * t263;
t86 = t251 * t143 + t247 * t263 + t287;
t83 = t110 * t207 - t131 * t192;
t82 = -t109 * t207 + t131 * t190;
t81 = t142 * t225 + t175 + (-t143 - t182) * t224;
t80 = (t142 * t247 + t143 * t250) * t248 + t271;
t71 = t109 * t192 - t110 * t190;
t70 = t118 * t222 + t205 * t293 + t146;
t69 = t135 * t203 + t222 * t296 + t154;
t64 = (-t210 + t282) * t251 + t250 * t260;
t63 = t251 * t118 + t247 * t260 + t279;
t62 = t118 * t234 + t225 * t280 + t290;
t61 = t224 * t135 + t234 * t282 + t288;
t60 = t122 * t251 + (t247 * t97 - t250 * t96) * t248;
t59 = t117 * t205 + t203 * t295 + t144;
t58 = t122 * t234 + t224 * t96 + t225 * t97;
t53 = (t117 * t247 + t118 * t250) * t248 + t258;
t52 = t120 * t251 + (t247 * t94 - t250 * t93) * t248;
t51 = t119 * t251 + (t247 * t92 - t250 * t91) * t248;
t50 = t117 * t225 + (-t182 + t295) * t224 + t291;
t45 = t120 * t234 + t224 * t93 + t225 * t94;
t44 = t119 * t234 + t224 * t91 + t225 * t92;
t39 = t205 * t281 + t222 * t297 + t146;
t38 = t203 * t294 + t222 * t284 + t154;
t37 = (-t210 + t273) * t251 + t250 * t257;
t36 = t247 * t257 + t251 * t297 + t279;
t35 = t225 * t272 + t234 * t297 + t290;
t34 = t224 * t294 + t234 * t273 + t288;
t31 = t203 * t283 + t205 * t298 + t144;
t29 = (t247 * t298 + t250 * t297) * t248 + t258;
t28 = t298 * t225 + (-t182 + t283) * t224 + t291;
t100 = [m(2) + m(3) + m(4) + m(5) + t319; m(3) * t184 + m(4) * t102 + m(5) * t80 + m(6) * t53 + m(7) * t29; m(7) * (t29 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t53 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t80 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(4) * (t102 ^ 2 + t123 ^ 2 + t124 ^ 2) + m(3) * (t184 ^ 2 + t198 ^ 2 + t199 ^ 2) + (t10 + t27 + t14 + t52 + (t213 * t304 + t215 * t237 + t217 * t238) * t304) * t304 + (-t9 - t26 - t13 - t51 + (-t212 * t303 + t214 * t235 + t216 * t236) * t303 + (-t212 * t304 + t213 * t303 - t214 * t237 - t215 * t235 - t216 * t238 - t217 * t236) * t304) * t303 + ((t229 * t304 + t230 * t237 + t231 * t238) * t304 - (-t229 * t303 + t235 * t230 + t236 * t231) * t303 + t19 + t33 + t21 + t60 + ((t215 * t256 + t217 * t255) * t247 - (t214 * t256 + t216 * t255) * t250) * t248 ^ 2 + ((-t212 * t250 + t213 * t247 + t230 * t256 + t231 * t255) * t248 + t251 * t229) * t251) * t251; m(4) * t121 + m(5) * t81 + m(6) * t50 + m(7) * t28; (t58 / 0.2e1 + t265) * t251 + (t60 / 0.2e1 + t264) * t234 + (t52 / 0.2e1 + t267) * t225 + (t51 / 0.2e1 + t268) * t224 + m(7) * (t28 * t29 + t34 * t37 + t35 * t36) + m(6) * (t50 * t53 + t61 * t64 + t62 * t63) + m(5) * (t80 * t81 + t86 * t89 + t87 * t88) + m(4) * (t102 * t121 + t123 * t126 + t124 * t125) + ((-t44 / 0.2e1 - t270) * t250 + (t45 / 0.2e1 + t269) * t247) * t248; (t17 + t20 + t32 + t58) * t234 + (t6 + t12 + t25 + t45) * t225 + (t5 + t11 + t24 + t44) * t224 + m(7) * (t28 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t50 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t81 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t121 ^ 2 + t125 ^ 2 + t126 ^ 2); m(5) * t90 + m(6) * t59 + m(7) * t31; t266 * t251 + t264 * t222 + t267 * t205 + t268 * t203 + m(7) * (t29 * t31 + t36 * t39 + t37 * t38) + m(6) * (t53 * t59 + t63 * t70 + t64 * t69) + m(5) * (t80 * t90 + t86 * t99 + t87 * t98) + (t247 * t274 - t250 * t275) * t248; t266 * t234 + t274 * t225 + t275 * t224 + t265 * t222 + t269 * t205 + t270 * t203 + m(7) * (t28 * t31 + t34 * t38 + t35 * t39) + m(6) * (t50 * t59 + t61 * t69 + t62 * t70) + m(5) * (t81 * t90 + t88 * t98 + t89 * t99); (t16 + t30 + t18) * t222 + (t4 + t23 + t8) * t205 + (t3 + t7 + t22) * t203 + m(7) * (t31 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(6) * (t59 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t90 ^ 2 + t98 ^ 2 + t99 ^ 2); t207 * t319; m(7) * (t190 * t36 + t192 * t37 + t207 * t29) + m(6) * (t190 * t63 + t192 * t64 + t207 * t53); m(7) * (t190 * t35 + t192 * t34 + t207 * t28) + m(6) * (t190 * t62 + t192 * t61 + t207 * t50); m(7) * (t190 * t39 + t192 * t38 + t207 * t31) + m(6) * (t190 * t70 + t192 * t69 + t207 * t59); (t190 ^ 2 + t192 ^ 2 + t207 ^ 2) * t319; m(7) * t71; m(7) * (t29 * t71 + t36 * t83 + t37 * t82) + t251 * t316 + t9 * t315 + t19 * t313 + t10 * t314 + (t247 * t317 - t250 * t1 / 0.2e1) * t248; m(7) * (t28 * t71 + t34 * t82 + t35 * t83) + t234 * t316 + t225 * t317 + t17 * t313 + t5 * t315 + t6 * t314 + t224 * t318; m(7) * (t31 * t71 + t38 * t82 + t39 * t83) + t4 * t314 + t16 * t313 + t3 * t315 + t203 * t318 + t222 * t316 + t205 * t317; m(7) * (t190 * t83 + t192 * t82 + t207 * t71); t192 * t2 + t190 * t1 + t207 * t15 + m(7) * (t71 ^ 2 + t82 ^ 2 + t83 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t100(1) t100(2) t100(4) t100(7) t100(11) t100(16); t100(2) t100(3) t100(5) t100(8) t100(12) t100(17); t100(4) t100(5) t100(6) t100(9) t100(13) t100(18); t100(7) t100(8) t100(9) t100(10) t100(14) t100(19); t100(11) t100(12) t100(13) t100(14) t100(15) t100(20); t100(16) t100(17) t100(18) t100(19) t100(20) t100(21);];
Mq  = res;
