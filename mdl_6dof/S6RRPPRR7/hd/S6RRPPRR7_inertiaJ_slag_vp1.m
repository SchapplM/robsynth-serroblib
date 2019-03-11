% Calculate joint inertia matrix for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:12
% EndTime: 2019-03-09 09:17:22
% DurationCPUTime: 4.78s
% Computational Cost: add. (9654->547), mult. (24617->785), div. (0->0), fcn. (30558->10), ass. (0->247)
t229 = sin(pkin(6));
t230 = cos(pkin(6));
t233 = sin(qJ(2));
t236 = cos(qJ(2));
t183 = -Icges(5,3) * t230 + (-Icges(5,5) * t236 - Icges(5,6) * t233) * t229;
t185 = Icges(3,3) * t230 + (Icges(3,5) * t233 + Icges(3,6) * t236) * t229;
t187 = Icges(4,2) * t230 + (Icges(4,4) * t233 - Icges(4,6) * t236) * t229;
t188 = Icges(3,6) * t230 + (Icges(3,4) * t233 + Icges(3,2) * t236) * t229;
t190 = Icges(4,4) * t230 + (Icges(4,1) * t233 - Icges(4,5) * t236) * t229;
t191 = Icges(3,5) * t230 + (Icges(3,1) * t233 + Icges(3,4) * t236) * t229;
t283 = t229 * t236;
t285 = t229 * t233;
t306 = t188 * t283 + (t190 + t191) * t285 + (-t183 + t185 + t187) * t230;
t237 = cos(qJ(1));
t276 = t236 * t237;
t234 = sin(qJ(1));
t279 = t233 * t234;
t210 = -t230 * t276 + t279;
t277 = t234 * t236;
t278 = t233 * t237;
t211 = t230 * t278 + t277;
t282 = t229 * t237;
t129 = Icges(4,5) * t211 - Icges(4,6) * t282 + Icges(4,3) * t210;
t137 = Icges(3,4) * t211 - Icges(3,2) * t210 - Icges(3,6) * t282;
t139 = Icges(5,1) * t210 - Icges(5,4) * t211 + Icges(5,5) * t282;
t305 = t129 - t137 + t139;
t133 = Icges(5,4) * t210 - Icges(5,2) * t211 + Icges(5,6) * t282;
t141 = Icges(4,1) * t211 - Icges(4,4) * t282 + Icges(4,5) * t210;
t143 = Icges(3,1) * t211 - Icges(3,4) * t210 - Icges(3,5) * t282;
t304 = t133 - t141 - t143;
t212 = t230 * t277 + t278;
t213 = -t230 * t279 + t276;
t284 = t229 * t234;
t134 = Icges(5,4) * t212 - Icges(5,2) * t213 - Icges(5,6) * t284;
t142 = Icges(4,1) * t213 + Icges(4,4) * t284 + Icges(4,5) * t212;
t144 = Icges(3,1) * t213 - Icges(3,4) * t212 + Icges(3,5) * t284;
t303 = t134 - t142 - t144;
t130 = Icges(4,5) * t213 + Icges(4,6) * t284 + Icges(4,3) * t212;
t138 = Icges(3,4) * t213 - Icges(3,2) * t212 + Icges(3,6) * t284;
t140 = Icges(5,1) * t212 - Icges(5,4) * t213 - Icges(5,5) * t284;
t302 = t140 - t138 + t130;
t228 = t229 ^ 2;
t301 = 0.2e1 * t229;
t300 = m(5) / 0.2e1;
t299 = m(6) / 0.2e1;
t298 = m(7) / 0.2e1;
t232 = sin(qJ(5));
t291 = cos(qJ(5));
t257 = t229 * t291;
t172 = t210 * t232 - t237 * t257;
t174 = t212 * t232 + t234 * t257;
t175 = t212 * t291 - t232 * t284;
t231 = sin(qJ(6));
t235 = cos(qJ(6));
t120 = -t175 * t231 + t213 * t235;
t121 = t175 * t235 + t213 * t231;
t173 = t210 * t291 + t232 * t282;
t118 = -t173 * t231 + t211 * t235;
t119 = t173 * t235 + t211 * t231;
t68 = Icges(7,5) * t119 + Icges(7,6) * t118 + Icges(7,3) * t172;
t70 = Icges(7,4) * t119 + Icges(7,2) * t118 + Icges(7,6) * t172;
t72 = Icges(7,1) * t119 + Icges(7,4) * t118 + Icges(7,5) * t172;
t19 = t120 * t70 + t121 * t72 + t174 * t68;
t69 = Icges(7,5) * t121 + Icges(7,6) * t120 + Icges(7,3) * t174;
t71 = Icges(7,4) * t121 + Icges(7,2) * t120 + Icges(7,6) * t174;
t73 = Icges(7,1) * t121 + Icges(7,4) * t120 + Icges(7,5) * t174;
t20 = t120 * t71 + t121 * t73 + t174 * t69;
t208 = -t230 * t291 + t232 * t283;
t209 = -t230 * t232 - t236 * t257;
t170 = -t209 * t231 + t235 * t285;
t171 = t209 * t235 + t231 * t285;
t95 = Icges(7,5) * t171 + Icges(7,6) * t170 - Icges(7,3) * t208;
t96 = Icges(7,4) * t171 + Icges(7,2) * t170 - Icges(7,6) * t208;
t97 = Icges(7,1) * t171 + Icges(7,4) * t170 - Icges(7,5) * t208;
t28 = t120 * t96 + t121 * t97 + t174 * t95;
t2 = t172 * t19 + t174 * t20 - t208 * t28;
t297 = t2 / 0.2e1;
t23 = t170 * t70 + t171 * t72 - t208 * t68;
t24 = t170 * t71 + t171 * t73 - t208 * t69;
t35 = t170 * t96 + t171 * t97 - t208 * t95;
t32 = t35 * t208;
t7 = t23 * t172 + t24 * t174 - t32;
t296 = t7 / 0.2e1;
t295 = -pkin(2) - pkin(3);
t294 = t172 / 0.2e1;
t293 = t174 / 0.2e1;
t292 = -t208 / 0.2e1;
t290 = t173 * pkin(5);
t289 = t210 * rSges(5,1);
t248 = -t119 * rSges(7,1) - t118 * rSges(7,2);
t74 = t172 * rSges(7,3) - t248;
t288 = t172 * pkin(10) + t290 + t74;
t75 = t121 * rSges(7,1) + t120 * rSges(7,2) + t174 * rSges(7,3);
t287 = t175 * pkin(5) + pkin(10) * t174 + t75;
t98 = rSges(7,1) * t171 + rSges(7,2) * t170 - rSges(7,3) * t208;
t286 = pkin(5) * t209 - pkin(10) * t208 + t98;
t186 = -Icges(5,6) * t230 + (-Icges(5,4) * t236 - Icges(5,2) * t233) * t229;
t280 = t233 * t186;
t195 = t210 * qJ(3);
t159 = t211 * pkin(2) + t195;
t160 = t213 * pkin(2) + qJ(3) * t212;
t275 = t159 * t284 + t160 * t282;
t156 = t230 * t160;
t205 = t213 * pkin(3);
t265 = qJ(4) * t284;
t180 = t205 - t265;
t274 = t230 * t180 + t156;
t264 = qJ(4) * t282;
t179 = t211 * pkin(3) + t264;
t273 = -t159 - t179;
t271 = t212 * rSges(5,1) - t213 * rSges(5,2);
t214 = (pkin(2) * t233 - qJ(3) * t236) * t229;
t270 = -pkin(3) * t285 + qJ(4) * t230 - t214;
t269 = t237 * pkin(1) + pkin(8) * t284;
t184 = Icges(4,6) * t230 + (Icges(4,5) * t233 - Icges(4,3) * t236) * t229;
t189 = -Icges(5,5) * t230 + (-Icges(5,1) * t236 - Icges(5,4) * t233) * t229;
t268 = (-t184 * t283 + (-t189 * t236 - t280) * t229 + t306) * t230;
t27 = t118 * t96 + t119 * t97 + t172 * t95;
t267 = t23 / 0.2e1 + t27 / 0.2e1;
t266 = t28 / 0.2e1 + t24 / 0.2e1;
t123 = Icges(6,5) * t209 + Icges(6,6) * t208 + Icges(6,3) * t285;
t124 = Icges(6,4) * t209 + Icges(6,2) * t208 + Icges(6,6) * t285;
t125 = Icges(6,1) * t209 + Icges(6,4) * t208 + Icges(6,5) * t285;
t59 = t123 * t285 + t208 * t124 + t209 * t125;
t128 = Icges(5,5) * t212 - Icges(5,6) * t213 - Icges(5,3) * t284;
t132 = Icges(3,5) * t213 - Icges(3,6) * t212 + Icges(3,3) * t284;
t136 = Icges(4,4) * t213 + Icges(4,2) * t284 + Icges(4,6) * t212;
t263 = -t128 + t132 + t136;
t127 = Icges(5,5) * t210 - Icges(5,6) * t211 + Icges(5,3) * t282;
t131 = Icges(3,5) * t211 - Icges(3,6) * t210 - Icges(3,3) * t282;
t135 = Icges(4,4) * t211 - Icges(4,2) * t282 + Icges(4,6) * t210;
t262 = -t135 + t127 - t131;
t162 = t212 * pkin(4) + pkin(9) * t213;
t261 = t230 * t162 + t274;
t161 = t210 * pkin(4) + t211 * pkin(9);
t260 = -t161 + t273;
t106 = t175 * rSges(6,1) - t174 * rSges(6,2) + t213 * rSges(6,3);
t149 = t213 * rSges(4,1) + rSges(4,2) * t284 + t212 * rSges(4,3);
t150 = t213 * rSges(3,1) - t212 * rSges(3,2) + rSges(3,3) * t284;
t258 = -(-pkin(4) * t236 + pkin(9) * t233) * t229 + t270;
t256 = -t234 * pkin(1) + pkin(8) * t282;
t255 = t229 * (-rSges(5,3) - qJ(4));
t254 = t229 * (-rSges(4,2) * t230 - (rSges(4,1) * t233 - rSges(4,3) * t236) * t229 - t214);
t253 = t300 + t299 + t298;
t252 = t179 * t284 + t180 * t282 + t275;
t251 = -t195 + t256;
t250 = t229 * (rSges(5,3) * t230 - (-rSges(5,1) * t236 - rSges(5,2) * t233) * t229 + t270);
t249 = -t173 * rSges(6,1) + t172 * rSges(6,2);
t126 = rSges(6,1) * t209 + rSges(6,2) * t208 + rSges(6,3) * t285;
t247 = t229 * (-t126 + t258);
t246 = t160 + t269;
t101 = Icges(6,4) * t173 - Icges(6,2) * t172 + Icges(6,6) * t211;
t103 = Icges(6,1) * t173 - Icges(6,4) * t172 + Icges(6,5) * t211;
t99 = Icges(6,5) * t173 - Icges(6,6) * t172 + Icges(6,3) * t211;
t44 = t101 * t208 + t103 * t209 + t99 * t285;
t49 = t123 * t211 - t124 * t172 + t125 * t173;
t245 = t49 / 0.2e1 + t44 / 0.2e1 + t267;
t100 = Icges(6,5) * t175 - Icges(6,6) * t174 + Icges(6,3) * t213;
t102 = Icges(6,4) * t175 - Icges(6,2) * t174 + Icges(6,6) * t213;
t104 = Icges(6,1) * t175 - Icges(6,4) * t174 + Icges(6,5) * t213;
t45 = t100 * t285 + t102 * t208 + t104 * t209;
t50 = t123 * t213 - t124 * t174 + t125 * t175;
t244 = t50 / 0.2e1 + t45 / 0.2e1 + t266;
t243 = rSges(4,2) * t282 - t210 * rSges(4,3);
t242 = t161 * t284 + t162 * t282 + t252;
t241 = t229 * (t258 - t286);
t240 = t205 + t246;
t147 = t211 * rSges(3,1) - t210 * rSges(3,2) - rSges(3,3) * t282;
t239 = -t161 + t251 - t264;
t238 = t162 + t240 - t265;
t218 = rSges(2,1) * t237 - t234 * rSges(2,2);
t217 = -t234 * rSges(2,1) - rSges(2,2) * t237;
t193 = rSges(3,3) * t230 + (rSges(3,1) * t233 + rSges(3,2) * t236) * t229;
t148 = -rSges(5,3) * t284 + t271;
t146 = t211 * rSges(4,1) - t243;
t145 = -t211 * rSges(5,2) + rSges(5,3) * t282 + t289;
t114 = t150 + t269;
t113 = -t147 + t256;
t108 = -t230 * t147 - t193 * t282;
t107 = t150 * t230 - t193 * t284;
t105 = t211 * rSges(6,3) - t249;
t91 = t246 + t149;
t90 = (-rSges(4,1) - pkin(2)) * t211 + t243 + t251;
t88 = (t147 * t234 + t150 * t237) * t229;
t87 = t185 * t284 - t188 * t212 + t191 * t213;
t86 = t184 * t212 + t187 * t284 + t190 * t213;
t85 = -t183 * t284 - t186 * t213 + t189 * t212;
t84 = -t185 * t282 - t210 * t188 + t211 * t191;
t83 = t210 * t184 - t187 * t282 + t211 * t190;
t82 = t183 * t282 - t211 * t186 + t210 * t189;
t79 = t234 * t255 + t240 + t271;
t78 = -t289 + t237 * t255 + (rSges(5,2) + t295) * t211 + t251;
t77 = (-t146 - t159) * t230 + t237 * t254;
t76 = t149 * t230 + t234 * t254 + t156;
t67 = t106 * t285 - t126 * t213;
t66 = -t105 * t285 + t126 * t211;
t65 = -t128 * t230 + (-t134 * t233 - t140 * t236) * t229;
t64 = -t127 * t230 + (-t133 * t233 - t139 * t236) * t229;
t63 = t132 * t230 + (t138 * t236 + t144 * t233) * t229;
t62 = t131 * t230 + (t137 * t236 + t143 * t233) * t229;
t61 = t136 * t230 + (-t130 * t236 + t142 * t233) * t229;
t60 = t135 * t230 + (-t129 * t236 + t141 * t233) * t229;
t58 = t59 * t230;
t57 = t59 * t285;
t56 = (t146 * t234 + t149 * t237) * t229 + t275;
t55 = (-t145 + t273) * t230 + t237 * t250;
t54 = t148 * t230 + t234 * t250 + t274;
t53 = t238 + t106;
t52 = (-rSges(6,3) + t295) * t211 + t239 + t249;
t51 = t105 * t213 - t106 * t211;
t48 = (t145 * t234 + t148 * t237) * t229 + t252;
t47 = -t174 * t98 - t208 * t75;
t46 = t172 * t98 + t208 * t74;
t43 = (-t105 + t260) * t230 + t237 * t247;
t42 = t106 * t230 + t234 * t247 + t261;
t41 = t238 + t287;
t40 = -t290 + t295 * t211 + (-rSges(7,3) - pkin(10)) * t172 + t239 + t248;
t39 = t100 * t213 - t102 * t174 + t104 * t175;
t38 = -t101 * t174 + t103 * t175 + t213 * t99;
t37 = t100 * t211 - t102 * t172 + t104 * t173;
t36 = -t101 * t172 + t103 * t173 + t211 * t99;
t34 = t35 * t230;
t33 = t35 * t285;
t31 = -t172 * t75 + t174 * t74;
t30 = -t286 * t213 + t287 * t285;
t29 = t286 * t211 - t288 * t285;
t26 = (t105 * t234 + t106 * t237) * t229 + t242;
t25 = -t287 * t211 + t288 * t213;
t22 = (t260 - t288) * t230 + t237 * t241;
t21 = t287 * t230 + t234 * t241 + t261;
t18 = t118 * t71 + t119 * t73 + t172 * t69;
t17 = t118 * t70 + t119 * t72 + t172 * t68;
t16 = (t288 * t234 + t287 * t237) * t229 + t242;
t15 = t58 + (t45 * t234 - t44 * t237) * t229;
t14 = t44 * t211 + t45 * t213 + t57;
t13 = t50 * t230 + (t234 * t39 - t237 * t38) * t229;
t12 = t49 * t230 + (t234 * t37 - t237 * t36) * t229;
t11 = t211 * t38 + t213 * t39 + t50 * t285;
t10 = t211 * t36 + t213 * t37 + t49 * t285;
t9 = t34 + (-t23 * t237 + t24 * t234) * t229;
t8 = t23 * t211 + t24 * t213 + t33;
t6 = t28 * t230 + (-t19 * t237 + t20 * t234) * t229;
t5 = t27 * t230 + (-t17 * t237 + t18 * t234) * t229;
t4 = t19 * t211 + t20 * t213 + t28 * t285;
t3 = t17 * t211 + t18 * t213 + t27 * t285;
t1 = t17 * t172 + t174 * t18 - t208 * t27;
t80 = [Icges(2,3) + (-t280 + (-t184 - t189) * t236) * t229 + m(7) * (t40 ^ 2 + t41 ^ 2) + m(6) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t90 ^ 2 + t91 ^ 2) + m(5) * (t78 ^ 2 + t79 ^ 2) + m(3) * (t113 ^ 2 + t114 ^ 2) + m(2) * (t217 ^ 2 + t218 ^ 2) + t59 + t35 + t306; t58 + t34 + m(7) * (t21 * t41 + t22 * t40) + m(6) * (t42 * t53 + t43 * t52) + m(5) * (t54 * t79 + t55 * t78) + m(4) * (t76 * t91 + t77 * t90) + m(3) * (t107 * t114 + t108 * t113) + ((-t60 / 0.2e1 - t64 / 0.2e1 - t62 / 0.2e1 - t82 / 0.2e1 - t83 / 0.2e1 - t84 / 0.2e1 - t245) * t237 + (t61 / 0.2e1 + t65 / 0.2e1 + t63 / 0.2e1 + t85 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1 + t244) * t234) * t229 + t268; (t9 + t15 + t268) * t230 + m(7) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(6) * (t26 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t48 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t56 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2 + t88 ^ 2) + ((-t12 - t5 + (t262 * t228 * t237 + (t210 * t305 - t211 * t304) * t229) * t237 + (-t60 - t62 - t64 - t82 - t83 - t84) * t230) * t237 + (t6 + t13 + (t263 * t228 * t234 + (t212 * t302 - t213 * t303) * t229) * t234 + (t85 + t87 + t86 + t61 + t65 + t63) * t230 + ((t262 * t234 + t263 * t237) * t229 + t304 * t213 - t305 * t212 + t303 * t211 - t302 * t210) * t282) * t234) * t229; m(7) * (t210 * t41 + t212 * t40) + m(6) * (t210 * t53 + t212 * t52) + m(4) * (t210 * t91 + t212 * t90) + m(5) * (t210 * t79 + t212 * t78); m(7) * (-t16 * t283 + t21 * t210 + t212 * t22) + m(6) * (t210 * t42 + t212 * t43 - t26 * t283) + m(5) * (t210 * t54 + t212 * t55 - t48 * t283) + m(4) * (t210 * t76 + t212 * t77 - t56 * t283); 0.2e1 * (m(4) / 0.2e1 + t253) * (t228 * t236 ^ 2 + t210 ^ 2 + t212 ^ 2); ((-t234 * t40 + t237 * t41) * t298 + (-t234 * t52 + t237 * t53) * t299 + (-t234 * t78 + t237 * t79) * t300) * t301; m(7) * (-t230 * t16 + (t21 * t237 - t22 * t234) * t229) + m(6) * (-t230 * t26 + (-t234 * t43 + t237 * t42) * t229) + m(5) * (-t230 * t48 + (-t234 * t55 + t237 * t54) * t229); t253 * (t210 * t237 - t212 * t234 + t230 * t236) * t301; 0.2e1 * t253 * (t230 ^ 2 + (t234 ^ 2 + t237 ^ 2) * t228); t33 + t57 + m(7) * (t29 * t40 + t30 * t41) + m(6) * (t52 * t66 + t53 * t67) + t244 * t213 + t245 * t211; (t8 / 0.2e1 + t14 / 0.2e1) * t230 + (t6 / 0.2e1 + t13 / 0.2e1) * t213 + (t5 / 0.2e1 + t12 / 0.2e1) * t211 + m(7) * (t16 * t25 + t21 * t30 + t22 * t29) + m(6) * (t26 * t51 + t42 * t67 + t43 * t66) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t237 + (t4 / 0.2e1 + t11 / 0.2e1) * t234 + (t9 / 0.2e1 + t15 / 0.2e1) * t233) * t229; m(6) * (t210 * t67 + t212 * t66 - t51 * t283) + m(7) * (t210 * t30 + t212 * t29 - t25 * t283); m(6) * (-t51 * t230 + (-t234 * t66 + t237 * t67) * t229) + m(7) * (-t25 * t230 + (-t234 * t29 + t237 * t30) * t229); (t14 + t8) * t285 + (t4 + t11) * t213 + (t3 + t10) * t211 + m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t51 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * (t40 * t46 + t41 * t47) - t32 + t266 * t174 + t267 * t172; m(7) * (t16 * t31 + t21 * t47 + t22 * t46) + t5 * t294 + t6 * t293 + t230 * t296 + t9 * t292 + (-t237 * t1 / 0.2e1 + t234 * t297) * t229; m(7) * (t210 * t47 + t212 * t46 - t283 * t31); m(7) * (-t31 * t230 + (-t234 * t46 + t237 * t47) * t229); m(7) * (t25 * t31 + t29 * t46 + t30 * t47) + t8 * t292 + t211 * t1 / 0.2e1 + t4 * t293 + t3 * t294 + t213 * t297 + t285 * t296; t174 * t2 + t172 * t1 - t208 * t7 + m(7) * (t31 ^ 2 + t46 ^ 2 + t47 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t80(1) t80(2) t80(4) t80(7) t80(11) t80(16); t80(2) t80(3) t80(5) t80(8) t80(12) t80(17); t80(4) t80(5) t80(6) t80(9) t80(13) t80(18); t80(7) t80(8) t80(9) t80(10) t80(14) t80(19); t80(11) t80(12) t80(13) t80(14) t80(15) t80(20); t80(16) t80(17) t80(18) t80(19) t80(20) t80(21);];
Mq  = res;
