% Calculate joint inertia matrix for
% S6RRPPRR9
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:26
% EndTime: 2019-03-09 09:28:38
% DurationCPUTime: 4.86s
% Computational Cost: add. (9622->543), mult. (24540->772), div. (0->0), fcn. (30449->10), ass. (0->241)
t234 = sin(pkin(6));
t235 = cos(pkin(6));
t238 = sin(qJ(2));
t241 = cos(qJ(2));
t185 = Icges(3,3) * t235 + (Icges(3,5) * t238 + Icges(3,6) * t241) * t234;
t186 = Icges(3,6) * t235 + (Icges(3,4) * t238 + Icges(3,2) * t241) * t234;
t187 = Icges(3,5) * t235 + (Icges(3,1) * t238 + Icges(3,4) * t241) * t234;
t188 = Icges(5,5) * t235 + (-Icges(5,6) * t241 + Icges(5,3) * t238) * t234;
t192 = Icges(5,1) * t235 + (-Icges(5,4) * t241 + Icges(5,5) * t238) * t234;
t193 = Icges(4,1) * t235 + (-Icges(4,4) * t238 - Icges(4,5) * t241) * t234;
t283 = t234 * t241;
t285 = t234 * t238;
t303 = t186 * t283 + (t187 + t188) * t285 + (t185 + t192 + t193) * t235;
t258 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t302 = 0.2e1 * t258;
t242 = cos(qJ(1));
t277 = t241 * t242;
t239 = sin(qJ(1));
t280 = t238 * t239;
t211 = -t235 * t277 + t280;
t278 = t239 * t241;
t279 = t238 * t242;
t212 = t235 * t279 + t278;
t282 = t234 * t242;
t128 = -Icges(5,5) * t282 + Icges(5,6) * t211 + Icges(5,3) * t212;
t134 = -Icges(4,4) * t282 - Icges(4,2) * t212 + Icges(4,6) * t211;
t143 = Icges(3,1) * t212 - Icges(3,4) * t211 - Icges(3,5) * t282;
t301 = t128 - t134 + t143;
t130 = -Icges(4,5) * t282 - Icges(4,6) * t212 + Icges(4,3) * t211;
t132 = -Icges(5,4) * t282 + Icges(5,2) * t211 + Icges(5,6) * t212;
t141 = Icges(3,4) * t212 - Icges(3,2) * t211 - Icges(3,6) * t282;
t300 = t130 + t132 - t141;
t213 = t235 * t278 + t279;
t214 = -t235 * t280 + t277;
t284 = t234 * t239;
t129 = Icges(4,5) * t284 - Icges(4,6) * t214 + Icges(4,3) * t213;
t131 = Icges(5,4) * t284 + Icges(5,2) * t213 + Icges(5,6) * t214;
t142 = Icges(3,4) * t214 - Icges(3,2) * t213 + Icges(3,6) * t284;
t299 = -t142 + t129 + t131;
t127 = Icges(5,5) * t284 + Icges(5,6) * t213 + Icges(5,3) * t214;
t133 = Icges(4,4) * t284 - Icges(4,2) * t214 + Icges(4,6) * t213;
t144 = Icges(3,1) * t214 - Icges(3,4) * t213 + Icges(3,5) * t284;
t298 = t144 - t133 + t127;
t233 = t234 ^ 2;
t237 = sin(qJ(5));
t292 = cos(qJ(5));
t171 = -t214 * t292 + t237 * t284;
t173 = t212 * t292 + t237 * t282;
t261 = t234 * t292;
t174 = t212 * t237 - t242 * t261;
t236 = sin(qJ(6));
t240 = cos(qJ(6));
t120 = -t174 * t236 - t211 * t240;
t121 = t174 * t240 - t211 * t236;
t172 = t214 * t237 + t239 * t261;
t118 = -t172 * t236 - t213 * t240;
t119 = t172 * t240 - t213 * t236;
t68 = Icges(7,5) * t119 + Icges(7,6) * t118 + Icges(7,3) * t171;
t70 = Icges(7,4) * t119 + Icges(7,2) * t118 + Icges(7,6) * t171;
t72 = Icges(7,1) * t119 + Icges(7,4) * t118 + Icges(7,5) * t171;
t19 = t120 * t70 + t121 * t72 - t173 * t68;
t69 = Icges(7,5) * t121 + Icges(7,6) * t120 - Icges(7,3) * t173;
t71 = Icges(7,4) * t121 + Icges(7,2) * t120 - Icges(7,6) * t173;
t73 = Icges(7,1) * t121 + Icges(7,4) * t120 - Icges(7,5) * t173;
t20 = t120 * t71 + t121 * t73 - t173 * t69;
t209 = t235 * t237 - t238 * t261;
t210 = t235 * t292 + t237 * t285;
t168 = -t210 * t236 + t240 * t283;
t169 = t210 * t240 + t236 * t283;
t95 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t209;
t96 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t209;
t97 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t209;
t28 = t120 * t96 + t121 * t97 - t173 * t95;
t2 = t171 * t19 - t173 * t20 + t209 * t28;
t297 = -t2 / 0.2e1;
t23 = t168 * t70 + t169 * t72 + t209 * t68;
t24 = t168 * t71 + t169 * t73 + t209 * t69;
t35 = t168 * t96 + t169 * t97 + t209 * t95;
t32 = t35 * t209;
t7 = t23 * t171 - t24 * t173 + t32;
t296 = t7 / 0.2e1;
t295 = t171 / 0.2e1;
t294 = -t173 / 0.2e1;
t293 = t209 / 0.2e1;
t291 = pkin(5) * t174;
t290 = t212 * pkin(2);
t74 = t119 * rSges(7,1) + t118 * rSges(7,2) + t171 * rSges(7,3);
t289 = t172 * pkin(5) + pkin(10) * t171 + t74;
t254 = -rSges(7,1) * t121 - rSges(7,2) * t120;
t75 = -rSges(7,3) * t173 - t254;
t288 = -pkin(10) * t173 + t291 + t75;
t98 = rSges(7,1) * t169 + rSges(7,2) * t168 + rSges(7,3) * t209;
t287 = pkin(5) * t210 + pkin(10) * t209 + t98;
t286 = qJ(3) * t213;
t191 = Icges(4,4) * t235 + (-Icges(4,2) * t238 - Icges(4,6) * t241) * t234;
t281 = t238 * t191;
t198 = t211 * qJ(3);
t157 = t198 + t290;
t208 = t214 * pkin(2);
t158 = t208 + t286;
t276 = t157 * t284 + t158 * t282;
t155 = t235 * t158;
t178 = pkin(3) * t284 + qJ(4) * t214;
t275 = t235 * t178 + t155;
t179 = -pkin(3) * t282 + t212 * qJ(4);
t274 = -t157 - t179;
t215 = (pkin(2) * t238 - qJ(3) * t241) * t234;
t272 = -pkin(3) * t235 - qJ(4) * t285 - t215;
t271 = t242 * pkin(1) + pkin(8) * t284;
t181 = -pkin(4) * t282 - t211 * pkin(9);
t189 = Icges(4,5) * t235 + (-Icges(4,6) * t238 - Icges(4,3) * t241) * t234;
t190 = Icges(5,4) * t235 + (-Icges(5,2) * t241 + Icges(5,6) * t238) * t234;
t270 = (-t190 * t283 + (-t189 * t241 - t281) * t234 + t303) * t235;
t27 = t118 * t96 + t119 * t97 + t171 * t95;
t269 = t23 / 0.2e1 + t27 / 0.2e1;
t268 = -t28 / 0.2e1 - t24 / 0.2e1;
t123 = Icges(6,5) * t210 - Icges(6,6) * t209 + Icges(6,3) * t283;
t124 = Icges(6,4) * t210 - Icges(6,2) * t209 + Icges(6,6) * t283;
t125 = Icges(6,1) * t210 - Icges(6,4) * t209 + Icges(6,5) * t283;
t59 = t123 * t283 - t209 * t124 + t210 * t125;
t136 = -Icges(5,1) * t282 + Icges(5,4) * t211 + Icges(5,5) * t212;
t138 = -Icges(4,1) * t282 - Icges(4,4) * t212 + Icges(4,5) * t211;
t139 = Icges(3,5) * t212 - Icges(3,6) * t211 - Icges(3,3) * t282;
t267 = -t139 - t138 - t136;
t135 = Icges(5,1) * t284 + Icges(5,4) * t213 + Icges(5,5) * t214;
t137 = Icges(4,1) * t284 - Icges(4,4) * t214 + Icges(4,5) * t213;
t140 = Icges(3,5) * t214 - Icges(3,6) * t213 + Icges(3,3) * t284;
t266 = t140 + t137 + t135;
t227 = pkin(4) * t284;
t180 = -pkin(9) * t213 + t227;
t265 = t235 * t180 + t275;
t264 = -t181 + t274;
t105 = t172 * rSges(6,1) - t171 * rSges(6,2) - t213 * rSges(6,3);
t150 = t214 * rSges(3,1) - t213 * rSges(3,2) + rSges(3,3) * t284;
t262 = -pkin(4) * t235 - pkin(9) * t283 + t272;
t145 = rSges(5,1) * t284 + t213 * rSges(5,2) + t214 * rSges(5,3);
t146 = rSges(4,1) * t284 - t214 * rSges(4,2) + t213 * rSges(4,3);
t260 = -t239 * pkin(1) + pkin(8) * t282;
t259 = t234 * (-t235 * rSges(4,1) - (-rSges(4,2) * t238 - rSges(4,3) * t241) * t234 - t215);
t257 = t178 * t282 + t179 * t284 + t276;
t256 = -t198 + t260;
t255 = t234 * (-t235 * rSges(5,1) - (-rSges(5,2) * t241 + rSges(5,3) * t238) * t234 + t272);
t126 = rSges(6,1) * t210 - rSges(6,2) * t209 + rSges(6,3) * t283;
t253 = t234 * (-t126 + t262);
t101 = Icges(6,4) * t172 - Icges(6,2) * t171 - Icges(6,6) * t213;
t103 = Icges(6,1) * t172 - Icges(6,4) * t171 - Icges(6,5) * t213;
t99 = Icges(6,5) * t172 - Icges(6,6) * t171 - Icges(6,3) * t213;
t44 = -t101 * t209 + t103 * t210 + t99 * t283;
t49 = -t123 * t213 - t124 * t171 + t125 * t172;
t252 = -t44 / 0.2e1 - t49 / 0.2e1 - t269;
t100 = Icges(6,5) * t174 + Icges(6,6) * t173 - Icges(6,3) * t211;
t102 = Icges(6,4) * t174 + Icges(6,2) * t173 - Icges(6,6) * t211;
t104 = Icges(6,1) * t174 + Icges(6,4) * t173 - Icges(6,5) * t211;
t45 = t100 * t283 - t102 * t209 + t104 * t210;
t50 = -t123 * t211 + t124 * t173 + t125 * t174;
t251 = -t50 / 0.2e1 - t45 / 0.2e1 + t268;
t250 = rSges(4,1) * t282 - t211 * rSges(4,3);
t249 = rSges(5,1) * t282 - t211 * rSges(5,2);
t248 = t180 * t282 + t181 * t284 + t257;
t247 = t256 - t179;
t246 = t234 * (t262 - t287);
t245 = t178 + t208 + t271;
t106 = rSges(6,1) * t174 + rSges(6,2) * t173 - rSges(6,3) * t211;
t149 = t212 * rSges(3,1) - t211 * rSges(3,2) - rSges(3,3) * t282;
t244 = -t181 + t247 - t290;
t243 = t227 + (-pkin(9) + qJ(3)) * t213 + t245;
t219 = rSges(2,1) * t242 - t239 * rSges(2,2);
t218 = -t239 * rSges(2,1) - rSges(2,2) * t242;
t194 = t235 * rSges(3,3) + (rSges(3,1) * t238 + rSges(3,2) * t241) * t234;
t148 = -t212 * rSges(4,2) - t250;
t147 = t212 * rSges(5,3) - t249;
t115 = t150 + t271;
t114 = -t149 + t260;
t108 = -t235 * t149 - t194 * t282;
t107 = t150 * t235 - t194 * t284;
t91 = t158 + t146 + t271;
t90 = (rSges(4,2) - pkin(2)) * t212 + t250 + t256;
t88 = (t149 * t239 + t150 * t242) * t234;
t87 = t185 * t284 - t186 * t213 + t187 * t214;
t86 = -t185 * t282 - t211 * t186 + t212 * t187;
t85 = t211 * t189 - t212 * t191 - t193 * t282;
t84 = t212 * t188 + t211 * t190 - t192 * t282;
t83 = t189 * t213 - t191 * t214 + t193 * t284;
t82 = t188 * t214 + t190 * t213 + t192 * t284;
t79 = t145 + t245 + t286;
t78 = (-rSges(5,3) - pkin(2)) * t212 + t247 + t249;
t77 = (-t148 - t157) * t235 + t242 * t259;
t76 = t146 * t235 + t239 * t259 + t155;
t67 = t105 * t283 + t126 * t213;
t66 = -t106 * t283 - t126 * t211;
t65 = t138 * t235 + (-t130 * t241 - t134 * t238) * t234;
t64 = t137 * t235 + (-t129 * t241 - t133 * t238) * t234;
t63 = t136 * t235 + (t128 * t238 - t132 * t241) * t234;
t62 = t135 * t235 + (t127 * t238 - t131 * t241) * t234;
t61 = t140 * t235 + (t142 * t241 + t144 * t238) * t234;
t60 = t139 * t235 + (t141 * t241 + t143 * t238) * t234;
t58 = t243 + t105;
t57 = -t106 + t244;
t56 = t59 * t235;
t55 = t59 * t283;
t54 = (t146 * t242 + t148 * t239) * t234 + t276;
t53 = (-t147 + t274) * t235 + t242 * t255;
t52 = t145 * t235 + t239 * t255 + t275;
t51 = t105 * t211 - t106 * t213;
t48 = (t145 * t242 + t147 * t239) * t234 + t257;
t47 = -t173 * t98 - t209 * t75;
t46 = -t171 * t98 + t209 * t74;
t43 = (-t106 + t264) * t235 + t242 * t253;
t42 = t105 * t235 + t239 * t253 + t265;
t41 = t243 + t289;
t40 = -t291 + (rSges(7,3) + pkin(10)) * t173 + t244 + t254;
t39 = -t100 * t211 + t102 * t173 + t104 * t174;
t38 = t101 * t173 + t103 * t174 - t211 * t99;
t37 = -t100 * t213 - t102 * t171 + t104 * t172;
t36 = -t101 * t171 + t103 * t172 - t213 * t99;
t34 = t35 * t235;
t33 = t35 * t283;
t31 = t171 * t75 + t173 * t74;
t30 = t287 * t213 + t289 * t283;
t29 = -t287 * t211 - t288 * t283;
t26 = (t105 * t242 + t106 * t239) * t234 + t248;
t25 = t289 * t211 - t288 * t213;
t22 = (t264 - t288) * t235 + t242 * t246;
t21 = t289 * t235 + t239 * t246 + t265;
t18 = t118 * t71 + t119 * t73 + t171 * t69;
t17 = t118 * t70 + t119 * t72 + t171 * t68;
t16 = (t288 * t239 + t289 * t242) * t234 + t248;
t15 = t56 + (t44 * t239 - t45 * t242) * t234;
t14 = -t45 * t211 - t44 * t213 + t55;
t13 = t50 * t235 + (t239 * t38 - t242 * t39) * t234;
t12 = t49 * t235 + (t239 * t36 - t242 * t37) * t234;
t11 = -t211 * t39 - t213 * t38 + t50 * t283;
t10 = -t211 * t37 - t213 * t36 + t49 * t283;
t9 = t34 + (t23 * t239 - t24 * t242) * t234;
t8 = -t24 * t211 - t23 * t213 + t33;
t6 = t28 * t235 + (t19 * t239 - t20 * t242) * t234;
t5 = t27 * t235 + (t17 * t239 - t18 * t242) * t234;
t4 = -t19 * t213 - t20 * t211 + t28 * t283;
t3 = -t17 * t213 - t18 * t211 + t27 * t283;
t1 = t17 * t171 - t173 * t18 + t209 * t27;
t80 = [Icges(2,3) + (-t281 + (-t189 - t190) * t241) * t234 + m(7) * (t40 ^ 2 + t41 ^ 2) + m(6) * (t57 ^ 2 + t58 ^ 2) + m(5) * (t78 ^ 2 + t79 ^ 2) + m(4) * (t90 ^ 2 + t91 ^ 2) + m(3) * (t114 ^ 2 + t115 ^ 2) + m(2) * (t218 ^ 2 + t219 ^ 2) + t59 + t35 + t303; t56 + t34 + m(7) * (t21 * t41 + t22 * t40) + m(6) * (t42 * t58 + t43 * t57) + m(5) * (t52 * t79 + t53 * t78) + m(4) * (t76 * t91 + t77 * t90) + m(3) * (t107 * t115 + t108 * t114) + ((-t60 / 0.2e1 - t65 / 0.2e1 - t63 / 0.2e1 - t84 / 0.2e1 - t85 / 0.2e1 - t86 / 0.2e1 + t251) * t242 + (t61 / 0.2e1 + t64 / 0.2e1 + t62 / 0.2e1 + t82 / 0.2e1 + t83 / 0.2e1 + t87 / 0.2e1 - t252) * t239) * t234 + t270; (t9 + t15 + t270) * t235 + m(7) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(6) * (t26 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t48 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(4) * (t54 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2 + t88 ^ 2) + ((-t13 - t6 + (t267 * t233 * t242 + (t211 * t300 + t212 * t301) * t234) * t242 + (-t60 - t63 - t65 - t84 - t85 - t86) * t235) * t242 + (t5 + t12 + (t266 * t233 * t239 + (t213 * t299 + t214 * t298) * t234) * t239 + (t87 + t83 + t82 + t61 + t64 + t62) * t235 + ((t239 * t267 + t242 * t266) * t234 - t301 * t214 - t300 * t213 - t298 * t212 - t299 * t211) * t282) * t239) * t234; m(7) * (t211 * t41 + t213 * t40) + m(6) * (t211 * t58 + t213 * t57) + m(5) * (t211 * t79 + t213 * t78) + m(4) * (t211 * t91 + t213 * t90); m(7) * (-t16 * t283 + t21 * t211 + t213 * t22) + m(6) * (t211 * t42 + t213 * t43 - t26 * t283) + m(5) * (t211 * t52 + t213 * t53 - t283 * t48) + m(4) * (t211 * t76 + t213 * t77 - t283 * t54); 0.2e1 * (m(4) / 0.2e1 + t258) * (t233 * t241 ^ 2 + t211 ^ 2 + t213 ^ 2); m(7) * (t212 * t41 + t214 * t40) + m(6) * (t212 * t58 + t214 * t57) + m(5) * (t212 * t79 + t214 * t78); m(7) * (t16 * t285 + t21 * t212 + t214 * t22) + m(6) * (t212 * t42 + t214 * t43 + t26 * t285) + m(5) * (t212 * t52 + t214 * t53 + t285 * t48); (-t233 * t238 * t241 + t211 * t212 + t213 * t214) * t302; (t233 * t238 ^ 2 + t212 ^ 2 + t214 ^ 2) * t302; t33 + t55 + m(7) * (t29 * t40 + t30 * t41) + m(6) * (t57 * t66 + t58 * t67) + t252 * t213 + t251 * t211; (t8 / 0.2e1 + t14 / 0.2e1) * t235 + (-t5 / 0.2e1 - t12 / 0.2e1) * t213 + (-t6 / 0.2e1 - t13 / 0.2e1) * t211 + m(7) * (t16 * t25 + t21 * t30 + t22 * t29) + m(6) * (t26 * t51 + t42 * t67 + t43 * t66) + ((-t4 / 0.2e1 - t11 / 0.2e1) * t242 + (t9 / 0.2e1 + t15 / 0.2e1) * t241 + (t3 / 0.2e1 + t10 / 0.2e1) * t239) * t234; m(6) * (t211 * t67 + t213 * t66 - t283 * t51) + m(7) * (t211 * t30 + t213 * t29 - t25 * t283); m(6) * (t212 * t67 + t214 * t66 + t285 * t51) + m(7) * (t212 * t30 + t214 * t29 + t25 * t285); (t14 + t8) * t283 + (-t3 - t10) * t213 + (-t4 - t11) * t211 + m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t51 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * (t40 * t47 + t41 * t46) + t32 + t268 * t173 + t269 * t171; t6 * t294 + t9 * t293 + m(7) * (t16 * t31 + t21 * t46 + t22 * t47) + t235 * t296 + t5 * t295 + (t239 * t1 / 0.2e1 + t242 * t297) * t234; m(7) * (t211 * t46 + t213 * t47 - t283 * t31); m(7) * (t212 * t46 + t214 * t47 + t285 * t31); m(7) * (t25 * t31 + t29 * t47 + t30 * t46) + t211 * t297 + t3 * t295 - t213 * t1 / 0.2e1 + t4 * t294 + t8 * t293 + t283 * t296; t171 * t1 - t173 * t2 + t209 * t7 + m(7) * (t31 ^ 2 + t46 ^ 2 + t47 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t80(1) t80(2) t80(4) t80(7) t80(11) t80(16); t80(2) t80(3) t80(5) t80(8) t80(12) t80(17); t80(4) t80(5) t80(6) t80(9) t80(13) t80(18); t80(7) t80(8) t80(9) t80(10) t80(14) t80(19); t80(11) t80(12) t80(13) t80(14) t80(15) t80(20); t80(16) t80(17) t80(18) t80(19) t80(20) t80(21);];
Mq  = res;
