% Calculate joint inertia matrix for
% S6PRRRPR5
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:23:05
% EndTime: 2019-03-08 23:23:22
% DurationCPUTime: 7.55s
% Computational Cost: add. (52003->639), mult. (125591->908), div. (0->0), fcn. (164772->16), ass. (0->293)
t319 = m(6) + m(7);
t244 = sin(pkin(12));
t246 = cos(pkin(12));
t252 = sin(qJ(2));
t247 = cos(pkin(6));
t255 = cos(qJ(2));
t300 = t247 * t255;
t235 = -t244 * t252 + t246 * t300;
t301 = t247 * t252;
t236 = t244 * t255 + t246 * t301;
t251 = sin(qJ(3));
t245 = sin(pkin(6));
t308 = sin(pkin(7));
t277 = t245 * t308;
t309 = cos(pkin(7));
t312 = cos(qJ(3));
t204 = t236 * t312 + (t235 * t309 - t246 * t277) * t251;
t278 = t245 * t309;
t224 = -t235 * t308 - t246 * t278;
t286 = qJ(4) + pkin(13);
t243 = sin(t286);
t276 = cos(t286);
t183 = t204 * t243 - t224 * t276;
t237 = -t244 * t300 - t246 * t252;
t238 = -t244 * t301 + t246 * t255;
t206 = t238 * t312 + (t237 * t309 + t244 * t277) * t251;
t225 = -t237 * t308 + t244 * t278;
t185 = t206 * t243 - t225 * t276;
t223 = t247 * t308 * t251 + (t251 * t255 * t309 + t252 * t312) * t245;
t234 = t247 * t309 - t255 * t277;
t199 = t223 * t243 - t234 * t276;
t184 = t204 * t276 + t224 * t243;
t260 = t312 * t308;
t258 = t245 * t260;
t261 = t309 * t312;
t203 = -t235 * t261 + t236 * t251 + t246 * t258;
t249 = sin(qJ(6));
t253 = cos(qJ(6));
t150 = -t184 * t249 + t203 * t253;
t151 = t184 * t253 + t203 * t249;
t100 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t183;
t102 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t183;
t104 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t183;
t40 = t100 * t183 + t102 * t150 + t104 * t151;
t186 = t206 * t276 + t225 * t243;
t205 = -t237 * t261 + t238 * t251 - t244 * t258;
t152 = -t186 * t249 + t205 * t253;
t153 = t186 * t253 + t205 * t249;
t101 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t185;
t103 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t185;
t105 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t185;
t41 = t101 * t183 + t103 * t150 + t105 * t151;
t200 = t223 * t276 + t234 * t243;
t302 = t245 * t252;
t222 = -t245 * t255 * t261 - t247 * t260 + t251 * t302;
t181 = -t200 * t249 + t222 * t253;
t182 = t200 * t253 + t222 * t249;
t126 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t199;
t127 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t199;
t128 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t199;
t55 = t126 * t183 + t127 * t150 + t128 * t151;
t1 = t183 * t40 + t185 * t41 + t199 * t55;
t318 = t1 / 0.2e1;
t42 = t100 * t185 + t102 * t152 + t104 * t153;
t43 = t101 * t185 + t103 * t152 + t105 * t153;
t56 = t126 * t185 + t127 * t152 + t128 * t153;
t2 = t183 * t42 + t185 * t43 + t199 * t56;
t317 = t2 / 0.2e1;
t48 = t100 * t199 + t102 * t181 + t104 * t182;
t49 = t101 * t199 + t103 * t181 + t105 * t182;
t62 = t126 * t199 + t127 * t181 + t128 * t182;
t9 = t183 * t48 + t185 * t49 + t199 * t62;
t316 = t9 / 0.2e1;
t315 = t183 / 0.2e1;
t314 = t185 / 0.2e1;
t313 = t199 / 0.2e1;
t254 = cos(qJ(4));
t311 = pkin(4) * t254;
t250 = sin(qJ(4));
t307 = t224 * t250;
t306 = t225 * t250;
t305 = t234 * t250;
t304 = t244 * t245;
t303 = t245 * t246;
t106 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t183;
t299 = pkin(5) * t184 + pkin(11) * t183 + t106;
t107 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t185;
t298 = pkin(5) * t186 + pkin(11) * t185 + t107;
t124 = pkin(4) * t307 + qJ(5) * t203 + t204 * t311;
t178 = t204 * pkin(3) + t203 * pkin(10);
t171 = t225 * t178;
t297 = t225 * t124 + t171;
t125 = pkin(4) * t306 + qJ(5) * t205 + t206 * t311;
t179 = t206 * pkin(3) + t205 * pkin(10);
t173 = t234 * t179;
t296 = t234 * t125 + t173;
t136 = rSges(6,1) * t184 - rSges(6,2) * t183 + rSges(6,3) * t203;
t295 = -t124 - t136;
t137 = rSges(6,1) * t186 - rSges(6,2) * t185 + rSges(6,3) * t205;
t294 = -t125 - t137;
t129 = rSges(7,1) * t182 + rSges(7,2) * t181 + rSges(7,3) * t199;
t293 = pkin(5) * t200 + pkin(11) * t199 + t129;
t188 = -t204 * t250 + t224 * t254;
t189 = t204 * t254 + t307;
t144 = rSges(5,1) * t189 + rSges(5,2) * t188 + rSges(5,3) * t203;
t292 = -t144 - t178;
t154 = pkin(4) * t305 + qJ(5) * t222 + t223 * t311;
t198 = t223 * pkin(3) + t222 * pkin(10);
t187 = t224 * t198;
t291 = t224 * t154 + t187;
t158 = rSges(6,1) * t200 - rSges(6,2) * t199 + rSges(6,3) * t222;
t290 = -t154 - t158;
t207 = -t223 * t250 + t234 * t254;
t208 = t223 * t254 + t305;
t170 = rSges(5,1) * t208 + rSges(5,2) * t207 + rSges(5,3) * t222;
t289 = -t170 - t198;
t211 = t238 * pkin(2) + pkin(9) * t225;
t209 = t247 * t211;
t288 = t247 * t179 + t209;
t210 = t236 * pkin(2) + pkin(9) * t224;
t287 = t210 * t304 + t211 * t303;
t284 = -t124 - t299;
t283 = -t125 - t298;
t282 = t247 * t125 + t288;
t281 = -t178 + t295;
t280 = -t154 - t293;
t279 = -t198 + t290;
t195 = rSges(4,1) * t223 - rSges(4,2) * t222 + rSges(4,3) * t234;
t226 = pkin(2) * t302 + pkin(9) * t234;
t275 = (-t195 - t226) * t245;
t274 = -t178 + t284;
t273 = -t198 + t280;
t272 = t178 * t304 + t179 * t303 + t287;
t130 = Icges(6,5) * t184 - Icges(6,6) * t183 + Icges(6,3) * t203;
t132 = Icges(6,4) * t184 - Icges(6,2) * t183 + Icges(6,6) * t203;
t134 = Icges(6,1) * t184 - Icges(6,4) * t183 + Icges(6,5) * t203;
t63 = t130 * t203 - t132 * t183 + t134 * t184;
t131 = Icges(6,5) * t186 - Icges(6,6) * t185 + Icges(6,3) * t205;
t133 = Icges(6,4) * t186 - Icges(6,2) * t185 + Icges(6,6) * t205;
t135 = Icges(6,1) * t186 - Icges(6,4) * t185 + Icges(6,5) * t205;
t64 = t131 * t203 - t133 * t183 + t135 * t184;
t155 = Icges(6,5) * t200 - Icges(6,6) * t199 + Icges(6,3) * t222;
t156 = Icges(6,4) * t200 - Icges(6,2) * t199 + Icges(6,6) * t222;
t157 = Icges(6,1) * t200 - Icges(6,4) * t199 + Icges(6,5) * t222;
t81 = t155 * t203 - t156 * t183 + t157 * t184;
t13 = t203 * t63 + t205 * t64 + t222 * t81;
t138 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t203;
t140 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t203;
t142 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t203;
t69 = t138 * t203 + t140 * t188 + t142 * t189;
t190 = -t206 * t250 + t225 * t254;
t191 = t206 * t254 + t306;
t139 = Icges(5,5) * t191 + Icges(5,6) * t190 + Icges(5,3) * t205;
t141 = Icges(5,4) * t191 + Icges(5,2) * t190 + Icges(5,6) * t205;
t143 = Icges(5,1) * t191 + Icges(5,4) * t190 + Icges(5,5) * t205;
t70 = t139 * t203 + t141 * t188 + t143 * t189;
t165 = Icges(5,5) * t208 + Icges(5,6) * t207 + Icges(5,3) * t222;
t166 = Icges(5,4) * t208 + Icges(5,2) * t207 + Icges(5,6) * t222;
t167 = Icges(5,1) * t208 + Icges(5,4) * t207 + Icges(5,5) * t222;
t83 = t165 * t203 + t166 * t188 + t167 * t189;
t17 = t203 * t69 + t205 * t70 + t222 * t83;
t3 = t203 * t40 + t205 * t41 + t222 * t55;
t271 = t3 / 0.2e1 + t13 / 0.2e1 + t17 / 0.2e1;
t65 = t130 * t205 - t132 * t185 + t134 * t186;
t66 = t131 * t205 - t133 * t185 + t135 * t186;
t82 = t155 * t205 - t156 * t185 + t157 * t186;
t14 = t203 * t65 + t205 * t66 + t222 * t82;
t71 = t138 * t205 + t140 * t190 + t142 * t191;
t72 = t139 * t205 + t141 * t190 + t143 * t191;
t84 = t165 * t205 + t166 * t190 + t167 * t191;
t18 = t203 * t71 + t205 * t72 + t222 * t84;
t4 = t203 * t42 + t205 * t43 + t222 * t56;
t270 = t4 / 0.2e1 + t14 / 0.2e1 + t18 / 0.2e1;
t15 = t224 * t63 + t225 * t64 + t234 * t81;
t21 = t224 * t69 + t225 * t70 + t234 * t83;
t5 = t224 * t40 + t225 * t41 + t234 * t55;
t269 = t5 / 0.2e1 + t15 / 0.2e1 + t21 / 0.2e1;
t16 = t224 * t65 + t225 * t66 + t234 * t82;
t22 = t224 * t71 + t225 * t72 + t234 * t84;
t6 = t224 * t42 + t225 * t43 + t234 * t56;
t268 = t6 / 0.2e1 + t16 / 0.2e1 + t22 / 0.2e1;
t19 = t81 * t247 + (t244 * t64 - t246 * t63) * t245;
t23 = t83 * t247 + (t244 * t70 - t246 * t69) * t245;
t7 = t55 * t247 + (t244 * t41 - t246 * t40) * t245;
t267 = t7 / 0.2e1 + t23 / 0.2e1 + t19 / 0.2e1;
t20 = t82 * t247 + (t244 * t66 - t246 * t65) * t245;
t24 = t84 * t247 + (t244 * t72 - t246 * t71) * t245;
t8 = t56 * t247 + (t244 * t43 - t246 * t42) * t245;
t266 = t8 / 0.2e1 + t24 / 0.2e1 + t20 / 0.2e1;
t10 = t203 * t48 + t205 * t49 + t222 * t62;
t73 = t130 * t222 - t132 * t199 + t134 * t200;
t74 = t131 * t222 - t133 * t199 + t135 * t200;
t89 = t155 * t222 - t156 * t199 + t157 * t200;
t25 = t203 * t73 + t205 * t74 + t222 * t89;
t75 = t138 * t222 + t140 * t207 + t142 * t208;
t76 = t139 * t222 + t141 * t207 + t143 * t208;
t95 = t165 * t222 + t166 * t207 + t167 * t208;
t27 = t203 * t75 + t205 * t76 + t222 * t95;
t265 = t10 / 0.2e1 + t25 / 0.2e1 + t27 / 0.2e1;
t11 = t224 * t48 + t225 * t49 + t234 * t62;
t26 = t224 * t73 + t225 * t74 + t234 * t89;
t29 = t224 * t75 + t225 * t76 + t234 * t95;
t264 = t11 / 0.2e1 + t29 / 0.2e1 + t26 / 0.2e1;
t12 = t62 * t247 + (t244 * t49 - t246 * t48) * t245;
t28 = t89 * t247 + (t244 * t74 - t246 * t73) * t245;
t30 = t95 * t247 + (t244 * t76 - t246 * t75) * t245;
t263 = t12 / 0.2e1 + t28 / 0.2e1 + t30 / 0.2e1;
t262 = (-t226 + t289) * t245;
t259 = (-t226 + t279) * t245;
t257 = t124 * t304 + t125 * t303 + t272;
t256 = (-t226 + t273) * t245;
t233 = t247 * rSges(3,3) + (rSges(3,1) * t252 + rSges(3,2) * t255) * t245;
t232 = Icges(3,5) * t247 + (Icges(3,1) * t252 + Icges(3,4) * t255) * t245;
t231 = Icges(3,6) * t247 + (Icges(3,4) * t252 + Icges(3,2) * t255) * t245;
t230 = Icges(3,3) * t247 + (Icges(3,5) * t252 + Icges(3,6) * t255) * t245;
t219 = rSges(3,1) * t238 + rSges(3,2) * t237 + rSges(3,3) * t304;
t218 = rSges(3,1) * t236 + rSges(3,2) * t235 - rSges(3,3) * t303;
t217 = Icges(3,1) * t238 + Icges(3,4) * t237 + Icges(3,5) * t304;
t216 = Icges(3,1) * t236 + Icges(3,4) * t235 - Icges(3,5) * t303;
t215 = Icges(3,4) * t238 + Icges(3,2) * t237 + Icges(3,6) * t304;
t214 = Icges(3,4) * t236 + Icges(3,2) * t235 - Icges(3,6) * t303;
t213 = Icges(3,5) * t238 + Icges(3,6) * t237 + Icges(3,3) * t304;
t212 = Icges(3,5) * t236 + Icges(3,6) * t235 - Icges(3,3) * t303;
t197 = -t218 * t247 - t233 * t303;
t196 = t219 * t247 - t233 * t304;
t194 = Icges(4,1) * t223 - Icges(4,4) * t222 + Icges(4,5) * t234;
t193 = Icges(4,4) * t223 - Icges(4,2) * t222 + Icges(4,6) * t234;
t192 = Icges(4,5) * t223 - Icges(4,6) * t222 + Icges(4,3) * t234;
t180 = (t218 * t244 + t219 * t246) * t245;
t169 = rSges(4,1) * t206 - rSges(4,2) * t205 + rSges(4,3) * t225;
t168 = rSges(4,1) * t204 - rSges(4,2) * t203 + rSges(4,3) * t224;
t164 = Icges(4,1) * t206 - Icges(4,4) * t205 + Icges(4,5) * t225;
t163 = Icges(4,1) * t204 - Icges(4,4) * t203 + Icges(4,5) * t224;
t162 = Icges(4,4) * t206 - Icges(4,2) * t205 + Icges(4,6) * t225;
t161 = Icges(4,4) * t204 - Icges(4,2) * t203 + Icges(4,6) * t224;
t160 = Icges(4,5) * t206 - Icges(4,6) * t205 + Icges(4,3) * t225;
t159 = Icges(4,5) * t204 - Icges(4,6) * t203 + Icges(4,3) * t224;
t146 = t203 * t154;
t145 = rSges(5,1) * t191 + rSges(5,2) * t190 + rSges(5,3) * t205;
t120 = t169 * t234 - t195 * t225;
t119 = -t168 * t234 + t195 * t224;
t116 = t222 * t125;
t115 = t205 * t124;
t114 = (-t168 - t210) * t247 + t246 * t275;
t113 = t247 * t169 + t244 * t275 + t209;
t112 = t192 * t234 - t193 * t222 + t194 * t223;
t111 = t168 * t225 - t169 * t224;
t110 = t192 * t225 - t193 * t205 + t194 * t206;
t109 = t192 * t224 - t193 * t203 + t194 * t204;
t108 = (t168 * t244 + t169 * t246) * t245 + t287;
t99 = t145 * t222 - t170 * t205;
t98 = -t144 * t222 + t170 * t203;
t97 = t160 * t234 - t162 * t222 + t164 * t223;
t96 = t159 * t234 - t161 * t222 + t163 * t223;
t94 = t160 * t225 - t162 * t205 + t164 * t206;
t93 = t159 * t225 - t161 * t205 + t163 * t206;
t92 = t160 * t224 - t162 * t203 + t164 * t204;
t91 = t159 * t224 - t161 * t203 + t163 * t204;
t90 = t144 * t205 - t145 * t203;
t88 = t234 * t145 + t225 * t289 + t173;
t87 = t224 * t170 + t234 * t292 + t187;
t86 = (-t210 + t292) * t247 + t246 * t262;
t85 = t247 * t145 + t244 * t262 + t288;
t80 = t107 * t199 - t129 * t185;
t79 = -t106 * t199 + t129 * t183;
t78 = t225 * t144 + t171 + (-t145 - t179) * t224;
t77 = (t144 * t244 + t145 * t246) * t245 + t272;
t68 = t222 * t137 + t205 * t290 + t116;
t67 = t203 * t158 + t222 * t295 + t146;
t61 = t106 * t185 - t107 * t183;
t60 = (-t210 + t281) * t247 + t246 * t259;
t59 = t247 * t137 + t244 * t259 + t282;
t58 = t234 * t137 + t225 * t279 + t296;
t57 = t224 * t158 + t234 * t281 + t291;
t54 = t205 * t136 + t203 * t294 + t115;
t53 = t112 * t247 + (t244 * t97 - t246 * t96) * t245;
t52 = (t136 * t244 + t137 * t246) * t245 + t257;
t51 = t112 * t234 + t224 * t96 + t225 * t97;
t50 = t225 * t136 + (-t179 + t294) * t224 + t297;
t47 = t110 * t247 + (t244 * t94 - t246 * t93) * t245;
t46 = t109 * t247 + (t244 * t92 - t246 * t91) * t245;
t45 = t110 * t234 + t224 * t93 + t225 * t94;
t44 = t109 * t234 + t224 * t91 + t225 * t92;
t39 = t205 * t280 + t222 * t298 + t116;
t38 = t203 * t293 + t222 * t284 + t146;
t37 = (-t210 + t274) * t247 + t246 * t256;
t36 = t244 * t256 + t247 * t298 + t282;
t35 = t225 * t273 + t234 * t298 + t296;
t34 = t224 * t293 + t234 * t274 + t291;
t33 = t203 * t283 + t205 * t299 + t115;
t32 = (t244 * t299 + t246 * t298) * t245 + t257;
t31 = t299 * t225 + (-t179 + t283) * t224 + t297;
t117 = [m(3) + m(4) + m(5) + m(2) + t319; m(3) * t180 + m(4) * t108 + m(5) * t77 + m(6) * t52 + m(7) * t32; m(7) * (t32 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t52 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t77 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t108 ^ 2 + t113 ^ 2 + t114 ^ 2) + m(3) * (t180 ^ 2 + t196 ^ 2 + t197 ^ 2) + (t8 + t24 + t20 + t47 + (t213 * t304 + t215 * t237 + t217 * t238) * t304) * t304 + (-t7 - t23 - t19 - t46 + (-t212 * t303 + t214 * t235 + t216 * t236) * t303 + (-t212 * t304 + t213 * t303 - t214 * t237 - t215 * t235 - t216 * t238 - t217 * t236) * t304) * t303 + ((t230 * t304 + t237 * t231 + t238 * t232) * t304 - (-t230 * t303 + t235 * t231 + t236 * t232) * t303 + t12 + t28 + t30 + t53 + ((t215 * t255 + t217 * t252) * t244 - (t214 * t255 + t216 * t252) * t246) * t245 ^ 2 + ((-t212 * t246 + t213 * t244 + t231 * t255 + t232 * t252) * t245 + t247 * t230) * t247) * t247; m(4) * t111 + m(5) * t78 + m(6) * t50 + m(7) * t31; (t51 / 0.2e1 + t264) * t247 + (t53 / 0.2e1 + t263) * t234 + (t47 / 0.2e1 + t266) * t225 + (t46 / 0.2e1 + t267) * t224 + m(7) * (t31 * t32 + t34 * t37 + t35 * t36) + m(6) * (t50 * t52 + t57 * t60 + t58 * t59) + m(5) * (t77 * t78 + t85 * t88 + t86 * t87) + m(4) * (t108 * t111 + t113 * t120 + t114 * t119) + ((-t44 / 0.2e1 - t269) * t246 + (t45 / 0.2e1 + t268) * t244) * t245; (t11 + t29 + t26 + t51) * t234 + (t6 + t16 + t22 + t45) * t225 + (t5 + t15 + t21 + t44) * t224 + m(7) * (t31 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t50 ^ 2 + t57 ^ 2 + t58 ^ 2) + m(5) * (t78 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t111 ^ 2 + t119 ^ 2 + t120 ^ 2); m(5) * t90 + m(6) * t54 + m(7) * t33; t265 * t247 + t263 * t222 + t266 * t205 + t267 * t203 + m(7) * (t32 * t33 + t36 * t39 + t37 * t38) + m(6) * (t52 * t54 + t59 * t68 + t60 * t67) + m(5) * (t77 * t90 + t85 * t99 + t86 * t98) + (t244 * t270 - t246 * t271) * t245; t265 * t234 + t270 * t225 + t271 * t224 + t264 * t222 + t268 * t205 + t269 * t203 + m(7) * (t31 * t33 + t34 * t38 + t35 * t39) + m(6) * (t50 * t54 + t57 * t67 + t58 * t68) + m(5) * (t78 * t90 + t87 * t98 + t88 * t99); (t10 + t27 + t25) * t222 + (t4 + t18 + t14) * t205 + (t3 + t17 + t13) * t203 + m(7) * (t33 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(6) * (t54 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(5) * (t90 ^ 2 + t98 ^ 2 + t99 ^ 2); t222 * t319; m(7) * (t203 * t36 + t205 * t37 + t222 * t32) + m(6) * (t203 * t59 + t205 * t60 + t222 * t52); m(7) * (t203 * t35 + t205 * t34 + t222 * t31) + m(6) * (t203 * t58 + t205 * t57 + t222 * t50); m(7) * (t203 * t39 + t205 * t38 + t222 * t33) + m(6) * (t203 * t68 + t205 * t67 + t222 * t54); (t203 ^ 2 + t205 ^ 2 + t222 ^ 2) * t319; m(7) * t61; m(7) * (t32 * t61 + t36 * t80 + t37 * t79) + t7 * t315 + t8 * t314 + t12 * t313 + t247 * t316 + (t244 * t317 - t246 * t1 / 0.2e1) * t245; t225 * t317 + t234 * t316 + t11 * t313 + t5 * t315 + t224 * t318 + t6 * t314 + m(7) * (t31 * t61 + t34 * t79 + t35 * t80); m(7) * (t33 * t61 + t38 * t79 + t39 * t80) + t222 * t316 + t10 * t313 + t3 * t315 + t4 * t314 + t205 * t317 + t203 * t318; m(7) * (t203 * t80 + t205 * t79 + t222 * t61); t183 * t1 + t199 * t9 + t185 * t2 + m(7) * (t61 ^ 2 + t79 ^ 2 + t80 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t117(1) t117(2) t117(4) t117(7) t117(11) t117(16); t117(2) t117(3) t117(5) t117(8) t117(12) t117(17); t117(4) t117(5) t117(6) t117(9) t117(13) t117(18); t117(7) t117(8) t117(9) t117(10) t117(14) t117(19); t117(11) t117(12) t117(13) t117(14) t117(15) t117(20); t117(16) t117(17) t117(18) t117(19) t117(20) t117(21);];
Mq  = res;
