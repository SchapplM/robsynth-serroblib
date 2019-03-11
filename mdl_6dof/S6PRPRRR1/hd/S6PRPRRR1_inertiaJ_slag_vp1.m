% Calculate joint inertia matrix for
% S6PRPRRR1
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:27
% EndTime: 2019-03-08 20:22:40
% DurationCPUTime: 6.14s
% Computational Cost: add. (30628->529), mult. (62562->781), div. (0->0), fcn. (82020->14), ass. (0->249)
t240 = sin(pkin(11));
t242 = cos(pkin(11));
t290 = sin(pkin(12));
t291 = cos(pkin(12));
t296 = sin(qJ(2));
t297 = cos(qJ(2));
t231 = -t296 * t290 + t297 * t291;
t243 = cos(pkin(6));
t249 = t243 * t231;
t250 = t290 * t297 + t291 * t296;
t208 = -t240 * t250 + t242 * t249;
t210 = -t240 * t249 - t242 * t250;
t241 = sin(pkin(6));
t221 = t231 * t241;
t223 = t250 * t243;
t209 = t223 * t242 + t231 * t240;
t284 = qJ(4) + qJ(5);
t238 = sin(t284);
t264 = cos(t284);
t288 = t241 * t242;
t183 = t209 * t264 - t238 * t288;
t244 = sin(qJ(6));
t246 = cos(qJ(6));
t150 = -t183 * t244 - t208 * t246;
t151 = t183 * t246 - t208 * t244;
t256 = t241 * t264;
t182 = t209 * t238 + t242 * t256;
t103 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t182;
t105 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t182;
t107 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t182;
t50 = t103 * t182 + t105 * t150 + t107 * t151;
t211 = -t223 * t240 + t231 * t242;
t289 = t240 * t241;
t185 = t211 * t264 + t238 * t289;
t152 = -t185 * t244 - t210 * t246;
t153 = t185 * t246 - t210 * t244;
t184 = t211 * t238 - t240 * t256;
t104 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t184;
t106 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t184;
t108 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t184;
t51 = t104 * t182 + t106 * t150 + t108 * t151;
t222 = t250 * t241;
t213 = t222 * t264 + t243 * t238;
t180 = -t213 * t244 - t221 * t246;
t181 = t213 * t246 - t221 * t244;
t212 = t222 * t238 - t243 * t264;
t132 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t212;
t133 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t212;
t134 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t212;
t63 = t132 * t182 + t133 * t150 + t134 * t151;
t11 = -t208 * t50 - t210 * t51 - t221 * t63;
t123 = Icges(6,5) * t183 - Icges(6,6) * t182 - Icges(6,3) * t208;
t125 = Icges(6,4) * t183 - Icges(6,2) * t182 - Icges(6,6) * t208;
t127 = Icges(6,1) * t183 - Icges(6,4) * t182 - Icges(6,5) * t208;
t70 = -t123 * t208 - t125 * t182 + t127 * t183;
t124 = Icges(6,5) * t185 - Icges(6,6) * t184 - Icges(6,3) * t210;
t126 = Icges(6,4) * t185 - Icges(6,2) * t184 - Icges(6,6) * t210;
t128 = Icges(6,1) * t185 - Icges(6,4) * t184 - Icges(6,5) * t210;
t71 = -t124 * t208 - t126 * t182 + t128 * t183;
t164 = Icges(6,5) * t213 - Icges(6,6) * t212 - Icges(6,3) * t221;
t165 = Icges(6,4) * t213 - Icges(6,2) * t212 - Icges(6,6) * t221;
t166 = Icges(6,1) * t213 - Icges(6,4) * t212 - Icges(6,5) * t221;
t88 = -t164 * t208 - t165 * t182 + t166 * t183;
t315 = -t208 * t70 - t210 * t71 - t221 * t88 + t11;
t52 = t103 * t184 + t105 * t152 + t107 * t153;
t53 = t104 * t184 + t106 * t152 + t108 * t153;
t64 = t132 * t184 + t133 * t152 + t134 * t153;
t12 = -t208 * t52 - t210 * t53 - t221 * t64;
t72 = -t123 * t210 - t125 * t184 + t127 * t185;
t73 = -t124 * t210 - t126 * t184 + t128 * t185;
t89 = -t164 * t210 - t165 * t184 + t166 * t185;
t314 = -t208 * t72 - t210 * t73 - t221 * t89 + t12;
t15 = t63 * t243 + (t240 * t51 - t242 * t50) * t241;
t313 = t15 + t88 * t243 + (t240 * t71 - t242 * t70) * t241;
t16 = t64 * t243 + (t240 * t53 - t242 * t52) * t241;
t312 = t16 + t89 * t243 + (t240 * t73 - t242 * t72) * t241;
t57 = t103 * t212 + t105 * t180 + t107 * t181;
t58 = t104 * t212 + t106 * t180 + t108 * t181;
t74 = t132 * t212 + t133 * t180 + t134 * t181;
t22 = -t208 * t57 - t210 * t58 - t221 * t74;
t80 = -t123 * t221 - t125 * t212 + t127 * t213;
t81 = -t124 * t221 - t126 * t212 + t128 * t213;
t94 = -t164 * t221 - t165 * t212 + t166 * t213;
t311 = -t208 * t80 - t210 * t81 - t221 * t94 + t22;
t24 = t74 * t243 + (t240 * t58 - t242 * t57) * t241;
t310 = t24 + t94 * t243 + (t240 * t81 - t242 * t80) * t241;
t109 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t182;
t282 = pkin(5) * t183 + pkin(10) * t182 + t109;
t110 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t184;
t281 = pkin(5) * t185 + pkin(10) * t184 + t110;
t135 = rSges(7,1) * t181 + rSges(7,2) * t180 + rSges(7,3) * t212;
t279 = pkin(5) * t213 + pkin(10) * t212 + t135;
t156 = Icges(4,5) * t209 + Icges(4,6) * t208 - Icges(4,3) * t288;
t268 = t243 * t297;
t226 = -t240 * t296 + t242 * t268;
t267 = t243 * t296;
t227 = t240 * t297 + t242 * t267;
t197 = Icges(3,5) * t227 + Icges(3,6) * t226 - Icges(3,3) * t288;
t309 = -t197 - t156;
t157 = Icges(4,5) * t211 + Icges(4,6) * t210 + Icges(4,3) * t289;
t228 = -t240 * t268 - t242 * t296;
t229 = -t240 * t267 + t242 * t297;
t198 = Icges(3,5) * t229 + Icges(3,6) * t228 + Icges(3,3) * t289;
t308 = t198 + t157;
t307 = Icges(4,5) * t222 + Icges(4,6) * t221 + (Icges(3,5) * t296 + Icges(3,6) * t297) * t241 + (Icges(4,3) + Icges(3,3)) * t243;
t306 = t182 / 0.2e1;
t305 = t184 / 0.2e1;
t304 = -t208 / 0.2e1;
t303 = -t210 / 0.2e1;
t302 = t212 / 0.2e1;
t301 = -t221 / 0.2e1;
t300 = t240 / 0.2e1;
t299 = -t242 / 0.2e1;
t298 = t243 / 0.2e1;
t295 = t297 * pkin(2);
t247 = cos(qJ(4));
t294 = pkin(4) * t247;
t292 = t282 * t210;
t245 = sin(qJ(4));
t287 = t241 * t245;
t286 = t241 * t247;
t285 = t243 * t245;
t283 = t281 * t221;
t280 = t279 * t208;
t172 = t211 * pkin(3) - t210 * pkin(8);
t263 = pkin(2) * t267 - qJ(3) * t241;
t206 = -t240 * t263 + t242 * t295;
t196 = t243 * t206;
t278 = t243 * t172 + t196;
t171 = t209 * pkin(3) - t208 * pkin(8);
t205 = t240 * t295 + t242 * t263;
t277 = -t171 - t205;
t276 = t205 * t289 + t206 * t288;
t232 = pkin(2) * t241 * t296 + t243 * qJ(3);
t275 = -t222 * pkin(3) + t221 * pkin(8) - t232;
t274 = m(4) + m(5) + m(6) + m(7);
t273 = t240 * t287;
t272 = t242 * t287;
t121 = pkin(4) * t273 - pkin(9) * t210 + t211 * t294;
t271 = t243 * t121 + t278;
t120 = -pkin(4) * t272 - pkin(9) * t208 + t209 * t294;
t270 = -t120 + t277;
t155 = pkin(4) * t285 - pkin(9) * t221 + t222 * t294;
t269 = -t155 + t275;
t262 = (-rSges(4,1) * t222 - rSges(4,2) * t221 - rSges(4,3) * t243 - t232) * t241;
t261 = t171 * t289 + t172 * t288 + t276;
t214 = -t222 * t245 + t243 * t247;
t215 = t222 * t247 + t285;
t176 = rSges(5,1) * t215 + rSges(5,2) * t214 - rSges(5,3) * t221;
t260 = (-t176 + t275) * t241;
t18 = t182 * t57 + t184 * t58 + t212 * t74;
t5 = t182 * t50 + t184 * t51 + t212 * t63;
t6 = t182 * t52 + t184 * t53 + t212 * t64;
t257 = t11 * t306 + t12 * t305 + t18 * t301 + t22 * t302 + t6 * t303 + t5 * t304;
t255 = -t208 * t315 - t314 * t210 - t311 * t221;
t167 = rSges(6,1) * t213 - rSges(6,2) * t212 - rSges(6,3) * t221;
t254 = (-t167 + t269) * t241;
t253 = t120 * t289 + t121 * t288 + t261;
t252 = (t269 - t279) * t241;
t251 = t313 * t304 + t312 * t303 + t310 * t301 + t311 * t298 + t314 * t289 / 0.2e1 - t315 * t288 / 0.2e1;
t220 = t243 * rSges(3,3) + (rSges(3,1) * t296 + rSges(3,2) * t297) * t241;
t219 = Icges(3,5) * t243 + (Icges(3,1) * t296 + Icges(3,4) * t297) * t241;
t218 = Icges(3,6) * t243 + (Icges(3,4) * t296 + Icges(3,2) * t297) * t241;
t204 = rSges(3,1) * t229 + rSges(3,2) * t228 + rSges(3,3) * t289;
t203 = rSges(3,1) * t227 + rSges(3,2) * t226 - rSges(3,3) * t288;
t202 = Icges(3,1) * t229 + Icges(3,4) * t228 + Icges(3,5) * t289;
t201 = Icges(3,1) * t227 + Icges(3,4) * t226 - Icges(3,5) * t288;
t200 = Icges(3,4) * t229 + Icges(3,2) * t228 + Icges(3,6) * t289;
t199 = Icges(3,4) * t227 + Icges(3,2) * t226 - Icges(3,6) * t288;
t194 = Icges(4,1) * t222 + Icges(4,4) * t221 + Icges(4,5) * t243;
t193 = Icges(4,4) * t222 + Icges(4,2) * t221 + Icges(4,6) * t243;
t189 = t211 * t247 + t273;
t188 = -t211 * t245 + t240 * t286;
t187 = t209 * t247 - t272;
t186 = -t209 * t245 - t242 * t286;
t178 = -t203 * t243 - t220 * t288;
t177 = t204 * t243 - t220 * t289;
t175 = Icges(5,1) * t215 + Icges(5,4) * t214 - Icges(5,5) * t221;
t174 = Icges(5,4) * t215 + Icges(5,2) * t214 - Icges(5,6) * t221;
t173 = Icges(5,5) * t215 + Icges(5,6) * t214 - Icges(5,3) * t221;
t163 = rSges(4,1) * t211 + rSges(4,2) * t210 + rSges(4,3) * t289;
t162 = rSges(4,1) * t209 + rSges(4,2) * t208 - rSges(4,3) * t288;
t161 = Icges(4,1) * t211 + Icges(4,4) * t210 + Icges(4,5) * t289;
t160 = Icges(4,1) * t209 + Icges(4,4) * t208 - Icges(4,5) * t288;
t159 = Icges(4,4) * t211 + Icges(4,2) * t210 + Icges(4,6) * t289;
t158 = Icges(4,4) * t209 + Icges(4,2) * t208 - Icges(4,6) * t288;
t154 = (t203 * t240 + t204 * t242) * t241;
t146 = t208 * t167;
t145 = t208 * t155;
t143 = rSges(5,1) * t189 + rSges(5,2) * t188 - rSges(5,3) * t210;
t142 = rSges(5,1) * t187 + rSges(5,2) * t186 - rSges(5,3) * t208;
t141 = Icges(5,1) * t189 + Icges(5,4) * t188 - Icges(5,5) * t210;
t140 = Icges(5,1) * t187 + Icges(5,4) * t186 - Icges(5,5) * t208;
t139 = Icges(5,4) * t189 + Icges(5,2) * t188 - Icges(5,6) * t210;
t138 = Icges(5,4) * t187 + Icges(5,2) * t186 - Icges(5,6) * t208;
t137 = Icges(5,5) * t189 + Icges(5,6) * t188 - Icges(5,3) * t210;
t136 = Icges(5,5) * t187 + Icges(5,6) * t186 - Icges(5,3) * t208;
t130 = rSges(6,1) * t185 - rSges(6,2) * t184 - rSges(6,3) * t210;
t129 = rSges(6,1) * t183 - rSges(6,2) * t182 - rSges(6,3) * t208;
t122 = t221 * t130;
t115 = t210 * t129;
t114 = t221 * t121;
t113 = (-t162 - t205) * t243 + t242 * t262;
t112 = t243 * t163 + t240 * t262 + t196;
t111 = t210 * t120;
t101 = (t162 * t240 + t163 * t242) * t241 + t276;
t100 = -t143 * t221 + t176 * t210;
t99 = t142 * t221 - t176 * t208;
t98 = t167 * t210 - t122;
t97 = t129 * t221 - t146;
t95 = -t173 * t221 + t174 * t214 + t175 * t215;
t93 = -t142 * t210 + t143 * t208;
t92 = -t173 * t210 + t174 * t188 + t175 * t189;
t91 = -t173 * t208 + t174 * t186 + t175 * t187;
t90 = t130 * t208 - t115;
t87 = (-t142 + t277) * t243 + t242 * t260;
t86 = t243 * t143 + t240 * t260 + t278;
t85 = t110 * t212 - t135 * t184;
t84 = -t109 * t212 + t135 * t182;
t83 = -t137 * t221 + t139 * t214 + t141 * t215;
t82 = -t136 * t221 + t138 * t214 + t140 * t215;
t79 = (t142 * t240 + t143 * t242) * t241 + t261;
t78 = -t137 * t210 + t139 * t188 + t141 * t189;
t77 = -t136 * t210 + t138 * t188 + t140 * t189;
t76 = -t137 * t208 + t139 * t186 + t141 * t187;
t75 = -t136 * t208 + t138 * t186 + t140 * t187;
t69 = -t114 - t122 + (t155 + t167) * t210;
t68 = -t145 - t146 - (-t120 - t129) * t221;
t67 = t109 * t184 - t110 * t182;
t66 = t210 * t279 - t283;
t65 = t221 * t282 - t280;
t62 = (-t129 + t270) * t243 + t242 * t254;
t61 = t243 * t130 + t240 * t254 + t271;
t60 = -t111 - t115 + (t121 + t130) * t208;
t59 = t208 * t281 - t292;
t56 = (t129 * t240 + t130 * t242) * t241 + t253;
t55 = -t114 + (t155 + t279) * t210 - t283;
t54 = -t145 - (-t120 - t282) * t221 - t280;
t49 = (t270 - t282) * t243 + t242 * t252;
t48 = t240 * t252 + t243 * t281 + t271;
t47 = -t111 + (t121 + t281) * t208 - t292;
t46 = t95 * t243 + (t240 * t83 - t242 * t82) * t241;
t45 = (t240 * t282 + t242 * t281) * t241 + t253;
t44 = -t208 * t82 - t210 * t83 - t221 * t95;
t38 = t92 * t243 + (t240 * t78 - t242 * t77) * t241;
t37 = t91 * t243 + (t240 * t76 - t242 * t75) * t241;
t36 = -t208 * t77 - t210 * t78 - t221 * t92;
t35 = -t208 * t75 - t210 * t76 - t221 * t91;
t1 = [m(2) + m(3) + t274; m(3) * t154 + m(4) * t101 + m(5) * t79 + m(6) * t56 + m(7) * t45; m(7) * (t45 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t56 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t79 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(4) * (t101 ^ 2 + t112 ^ 2 + t113 ^ 2) + m(3) * (t154 ^ 2 + t177 ^ 2 + t178 ^ 2) + (t46 + (t221 * t193 + t222 * t194 + t307 * t243) * t243 + ((t221 * t159 + t222 * t161) * t240 - (t221 * t158 + t222 * t160) * t242 + (-t156 * t242 + t157 * t240 + t218 * t297 + t219 * t296) * t243) * t241 + t310) * t243 + (t38 + (t210 * t159 + t211 * t161 + t228 * t200 + t229 * t202 + t308 * t289) * t289 + ((t200 * t297 + t202 * t296) * t241 + t228 * t218 + t229 * t219 + t210 * t193 + t211 * t194 + t307 * t289 + t243 * t198) * t243 + t312) * t289 + (-t37 + (t208 * t158 + t209 * t160 + t226 * t199 + t227 * t201 + t309 * t288) * t288 + (-(t199 * t297 + t201 * t296) * t241 - t226 * t218 - t227 * t219 - t208 * t193 - t209 * t194 + t307 * t288 - t243 * t197) * t243 + (-t210 * t158 - t208 * t159 - t211 * t160 - t209 * t161 - t228 * t199 - t226 * t200 - t229 * t201 - t227 * t202 + t308 * t288 + t309 * t289) * t289 - t313) * t288; t274 * t243; m(7) * (t243 * t45 + (t240 * t49 - t242 * t48) * t241) + m(6) * (t243 * t56 + (t240 * t62 - t242 * t61) * t241) + m(5) * (t243 * t79 + (t240 * t87 - t242 * t86) * t241) + m(4) * (t243 * t101 + (-t112 * t242 + t113 * t240) * t241); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1 + m(4) / 0.2e1) * (t243 ^ 2 + (t240 ^ 2 + t242 ^ 2) * t241 ^ 2); m(5) * t93 + m(6) * t60 + m(7) * t47; (t299 * t35 + t300 * t36) * t241 + m(7) * (t45 * t47 + t48 * t55 + t49 * t54) + m(6) * (t56 * t60 + t61 * t69 + t62 * t68) + m(5) * (t100 * t86 + t79 * t93 + t87 * t99) + t251 + t37 * t304 + t38 * t303 + t46 * t301 + t44 * t298; m(5) * (t93 * t243 + (-t100 * t242 + t240 * t99) * t241) + m(6) * (t60 * t243 + (t240 * t68 - t242 * t69) * t241) + m(7) * (t47 * t243 + (t240 * t54 - t242 * t55) * t241); -t208 * t35 - t210 * t36 - t221 * t44 + m(7) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t60 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(5) * (t100 ^ 2 + t93 ^ 2 + t99 ^ 2) + t255; m(6) * t90 + m(7) * t59; t251 + m(7) * (t45 * t59 + t48 * t66 + t49 * t65) + m(6) * (t56 * t90 + t61 * t98 + t62 * t97); m(6) * (t90 * t243 + (t240 * t97 - t242 * t98) * t241) + m(7) * (t59 * t243 + (t240 * t65 - t242 * t66) * t241); m(7) * (t47 * t59 + t54 * t65 + t55 * t66) + m(6) * (t60 * t90 + t68 * t97 + t69 * t98) + t255; m(7) * (t59 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t90 ^ 2 + t97 ^ 2 + t98 ^ 2) + t255; m(7) * t67; m(7) * (t45 * t67 + t48 * t85 + t49 * t84) + t15 * t306 + t16 * t305 + t18 * t298 + t24 * t302 + (t299 * t5 + t300 * t6) * t241; m(7) * (t67 * t243 + (t240 * t84 - t242 * t85) * t241); m(7) * (t47 * t67 + t54 * t84 + t55 * t85) + t257; m(7) * (t59 * t67 + t65 * t84 + t66 * t85) + t257; t182 * t5 + t212 * t18 + m(7) * (t67 ^ 2 + t84 ^ 2 + t85 ^ 2) + t184 * t6;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
