% Calculate joint inertia matrix for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:19
% EndTime: 2019-03-09 09:07:29
% DurationCPUTime: 4.89s
% Computational Cost: add. (9654->548), mult. (24617->779), div. (0->0), fcn. (30558->10), ass. (0->246)
t228 = sin(pkin(6));
t229 = cos(pkin(6));
t232 = sin(qJ(2));
t235 = cos(qJ(2));
t184 = -Icges(5,3) * t229 + (Icges(5,5) * t232 - Icges(5,6) * t235) * t228;
t185 = Icges(4,6) * t229 + (Icges(4,5) * t232 - Icges(4,3) * t235) * t228;
t186 = Icges(3,3) * t229 + (Icges(3,5) * t232 + Icges(3,6) * t235) * t228;
t187 = -Icges(5,6) * t229 + (Icges(5,4) * t232 - Icges(5,2) * t235) * t228;
t188 = Icges(4,2) * t229 + (Icges(4,4) * t232 - Icges(4,6) * t235) * t228;
t189 = Icges(3,6) * t229 + (Icges(3,4) * t232 + Icges(3,2) * t235) * t228;
t190 = -Icges(5,5) * t229 + (Icges(5,1) * t232 - Icges(5,4) * t235) * t228;
t191 = Icges(4,4) * t229 + (Icges(4,1) * t232 - Icges(4,5) * t235) * t228;
t192 = Icges(3,5) * t229 + (Icges(3,1) * t232 + Icges(3,4) * t235) * t228;
t278 = t228 * t235;
t280 = t228 * t232;
t301 = (t189 - t185 - t187) * t278 + (t190 + t192 + t191) * t280 + (-t184 + t186 + t188) * t229;
t236 = cos(qJ(1));
t273 = t235 * t236;
t233 = sin(qJ(1));
t276 = t232 * t233;
t211 = -t229 * t273 + t276;
t274 = t233 * t235;
t275 = t232 * t236;
t212 = t229 * t275 + t274;
t277 = t228 * t236;
t129 = Icges(4,5) * t212 - Icges(4,6) * t277 + Icges(4,3) * t211;
t133 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t277;
t137 = Icges(3,4) * t212 - Icges(3,2) * t211 - Icges(3,6) * t277;
t300 = t129 + t133 - t137;
t213 = t229 * t274 + t275;
t214 = -t229 * t276 + t273;
t279 = t228 * t233;
t130 = Icges(4,5) * t214 + Icges(4,6) * t279 + Icges(4,3) * t213;
t134 = Icges(5,4) * t214 + Icges(5,2) * t213 - Icges(5,6) * t279;
t138 = Icges(3,4) * t214 - Icges(3,2) * t213 + Icges(3,6) * t279;
t299 = t130 + t134 - t138;
t139 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t277;
t141 = Icges(4,1) * t212 - Icges(4,4) * t277 + Icges(4,5) * t211;
t143 = Icges(3,1) * t212 - Icges(3,4) * t211 - Icges(3,5) * t277;
t298 = t139 + t141 + t143;
t140 = Icges(5,1) * t214 + Icges(5,4) * t213 - Icges(5,5) * t279;
t142 = Icges(4,1) * t214 + Icges(4,4) * t279 + Icges(4,5) * t213;
t144 = Icges(3,1) * t214 - Icges(3,4) * t213 + Icges(3,5) * t279;
t297 = t142 + t140 + t144;
t227 = t228 ^ 2;
t296 = 0.2e1 * t228;
t295 = m(5) / 0.2e1;
t294 = m(6) / 0.2e1;
t293 = m(7) / 0.2e1;
t231 = sin(qJ(5));
t286 = cos(qJ(5));
t173 = t212 * t286 + t231 * t277;
t230 = sin(qJ(6));
t234 = cos(qJ(6));
t118 = -t173 * t230 - t211 * t234;
t119 = t173 * t234 - t211 * t230;
t255 = t228 * t286;
t172 = t212 * t231 - t236 * t255;
t68 = Icges(7,5) * t119 + Icges(7,6) * t118 + Icges(7,3) * t172;
t70 = Icges(7,4) * t119 + Icges(7,2) * t118 + Icges(7,6) * t172;
t72 = Icges(7,1) * t119 + Icges(7,4) * t118 + Icges(7,5) * t172;
t17 = t118 * t70 + t119 * t72 + t172 * t68;
t174 = t214 * t231 + t233 * t255;
t175 = t214 * t286 - t231 * t279;
t120 = -t175 * t230 - t213 * t234;
t121 = t175 * t234 - t213 * t230;
t69 = Icges(7,5) * t121 + Icges(7,6) * t120 + Icges(7,3) * t174;
t71 = Icges(7,4) * t121 + Icges(7,2) * t120 + Icges(7,6) * t174;
t73 = Icges(7,1) * t121 + Icges(7,4) * t120 + Icges(7,5) * t174;
t18 = t118 * t71 + t119 * t73 + t172 * t69;
t209 = t229 * t286 + t231 * t280;
t210 = -t229 * t231 + t232 * t255;
t170 = -t210 * t230 + t234 * t278;
t171 = t210 * t234 + t230 * t278;
t95 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t209;
t96 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t209;
t97 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t209;
t27 = t118 * t96 + t119 * t97 + t172 * t95;
t1 = t17 * t172 + t174 * t18 + t209 * t27;
t292 = -t1 / 0.2e1;
t23 = t170 * t70 + t171 * t72 + t209 * t68;
t24 = t170 * t71 + t171 * t73 + t209 * t69;
t35 = t170 * t96 + t171 * t97 + t209 * t95;
t32 = t35 * t209;
t7 = t23 * t172 + t24 * t174 + t32;
t291 = t7 / 0.2e1;
t290 = -pkin(2) - pkin(3);
t289 = t172 / 0.2e1;
t288 = t174 / 0.2e1;
t287 = t209 / 0.2e1;
t285 = t173 * pkin(5);
t284 = t211 * rSges(5,2);
t246 = -t119 * rSges(7,1) - t118 * rSges(7,2);
t74 = t172 * rSges(7,3) - t246;
t283 = t172 * pkin(10) + t285 + t74;
t75 = t121 * rSges(7,1) + t120 * rSges(7,2) + t174 * rSges(7,3);
t282 = t175 * pkin(5) + pkin(10) * t174 + t75;
t98 = rSges(7,1) * t171 + rSges(7,2) * t170 + rSges(7,3) * t209;
t281 = pkin(5) * t210 + pkin(10) * t209 + t98;
t196 = t211 * qJ(3);
t159 = t212 * pkin(2) + t196;
t207 = t214 * pkin(2);
t160 = qJ(3) * t213 + t207;
t272 = t159 * t279 + t160 * t277;
t156 = t229 * t160;
t206 = t214 * pkin(3);
t181 = -qJ(4) * t279 + t206;
t271 = t229 * t181 + t156;
t262 = qJ(4) * t277;
t180 = t212 * pkin(3) + t262;
t270 = -t159 - t180;
t268 = t214 * rSges(5,1) + t213 * rSges(5,2);
t215 = (pkin(2) * t232 - qJ(3) * t235) * t228;
t267 = -pkin(3) * t280 + qJ(4) * t229 - t215;
t266 = t236 * pkin(1) + pkin(8) * t279;
t265 = t301 * t229;
t264 = t23 / 0.2e1 + t27 / 0.2e1;
t28 = t120 * t96 + t121 * t97 + t174 * t95;
t263 = t28 / 0.2e1 + t24 / 0.2e1;
t123 = Icges(6,5) * t210 - Icges(6,6) * t209 + Icges(6,3) * t278;
t124 = Icges(6,4) * t210 - Icges(6,2) * t209 + Icges(6,6) * t278;
t125 = Icges(6,1) * t210 - Icges(6,4) * t209 + Icges(6,5) * t278;
t59 = t123 * t278 - t209 * t124 + t210 * t125;
t127 = Icges(5,5) * t212 + Icges(5,6) * t211 + Icges(5,3) * t277;
t131 = Icges(3,5) * t212 - Icges(3,6) * t211 - Icges(3,3) * t277;
t135 = Icges(4,4) * t212 - Icges(4,2) * t277 + Icges(4,6) * t211;
t261 = -t135 + t127 - t131;
t128 = Icges(5,5) * t214 + Icges(5,6) * t213 - Icges(5,3) * t279;
t132 = Icges(3,5) * t214 - Icges(3,6) * t213 + Icges(3,3) * t279;
t136 = Icges(4,4) * t214 + Icges(4,2) * t279 + Icges(4,6) * t213;
t260 = t136 - t128 + t132;
t205 = t214 * pkin(4);
t162 = -pkin(9) * t213 + t205;
t259 = t229 * t162 + t271;
t204 = t211 * pkin(9);
t161 = t212 * pkin(4) - t204;
t258 = -t161 + t270;
t106 = t175 * rSges(6,1) - t174 * rSges(6,2) - t213 * rSges(6,3);
t149 = t214 * rSges(4,1) + rSges(4,2) * t279 + t213 * rSges(4,3);
t150 = t214 * rSges(3,1) - t213 * rSges(3,2) + rSges(3,3) * t279;
t256 = -(pkin(4) * t232 + pkin(9) * t235) * t228 + t267;
t254 = -t233 * pkin(1) + pkin(8) * t277;
t253 = t228 * (-rSges(5,3) - qJ(4));
t252 = t228 * (-rSges(4,2) * t229 - (rSges(4,1) * t232 - rSges(4,3) * t235) * t228 - t215);
t250 = t295 + t294 + t293;
t249 = t180 * t279 + t181 * t277 + t272;
t248 = -t196 + t254;
t247 = t228 * (rSges(5,3) * t229 - (rSges(5,1) * t232 - rSges(5,2) * t235) * t228 + t267);
t126 = rSges(6,1) * t210 - rSges(6,2) * t209 + rSges(6,3) * t278;
t245 = t228 * (-t126 + t256);
t244 = t160 + t266;
t101 = Icges(6,4) * t173 - Icges(6,2) * t172 - Icges(6,6) * t211;
t103 = Icges(6,1) * t173 - Icges(6,4) * t172 - Icges(6,5) * t211;
t99 = Icges(6,5) * t173 - Icges(6,6) * t172 - Icges(6,3) * t211;
t44 = -t101 * t209 + t103 * t210 + t278 * t99;
t49 = -t123 * t211 - t124 * t172 + t125 * t173;
t243 = -t49 / 0.2e1 - t44 / 0.2e1 - t264;
t100 = Icges(6,5) * t175 - Icges(6,6) * t174 - Icges(6,3) * t213;
t102 = Icges(6,4) * t175 - Icges(6,2) * t174 - Icges(6,6) * t213;
t104 = Icges(6,1) * t175 - Icges(6,4) * t174 - Icges(6,5) * t213;
t45 = t100 * t278 - t102 * t209 + t104 * t210;
t50 = -t123 * t213 - t124 * t174 + t125 * t175;
t242 = -t50 / 0.2e1 - t45 / 0.2e1 - t263;
t241 = rSges(4,2) * t277 - t211 * rSges(4,3);
t240 = t161 * t279 + t162 * t277 + t249;
t239 = t228 * (t256 - t281);
t105 = t173 * rSges(6,1) - t172 * rSges(6,2) - t211 * rSges(6,3);
t147 = t212 * rSges(3,1) - t211 * rSges(3,2) - rSges(3,3) * t277;
t238 = t205 + t207 + (-pkin(9) + qJ(3)) * t213 + t181 + t266;
t237 = t204 + (-pkin(4) + t290) * t212 - t262 + t248;
t219 = rSges(2,1) * t236 - t233 * rSges(2,2);
t218 = -t233 * rSges(2,1) - rSges(2,2) * t236;
t195 = rSges(3,3) * t229 + (rSges(3,1) * t232 + rSges(3,2) * t235) * t228;
t148 = -rSges(5,3) * t279 + t268;
t146 = t212 * rSges(4,1) - t241;
t145 = t212 * rSges(5,1) + rSges(5,3) * t277 + t284;
t114 = t150 + t266;
t113 = -t147 + t254;
t108 = -t229 * t147 - t195 * t277;
t107 = t150 * t229 - t195 * t279;
t91 = t244 + t149;
t90 = (-rSges(4,1) - pkin(2)) * t212 + t241 + t248;
t88 = (t147 * t233 + t150 * t236) * t228;
t87 = t186 * t279 - t189 * t213 + t192 * t214;
t86 = t185 * t213 + t188 * t279 + t191 * t214;
t85 = -t184 * t279 + t187 * t213 + t190 * t214;
t84 = -t186 * t277 - t211 * t189 + t212 * t192;
t83 = t211 * t185 - t188 * t277 + t212 * t191;
t82 = t184 * t277 + t211 * t187 + t212 * t190;
t79 = t233 * t253 + t206 + t244 + t268;
t78 = -t284 + t236 * t253 + (-rSges(5,1) + t290) * t212 + t248;
t77 = (-t146 - t159) * t229 + t236 * t252;
t76 = t149 * t229 + t233 * t252 + t156;
t67 = t106 * t278 + t126 * t213;
t66 = -t105 * t278 - t126 * t211;
t65 = t132 * t229 + (t138 * t235 + t144 * t232) * t228;
t64 = t131 * t229 + (t137 * t235 + t143 * t232) * t228;
t63 = t136 * t229 + (-t130 * t235 + t142 * t232) * t228;
t62 = t135 * t229 + (-t129 * t235 + t141 * t232) * t228;
t61 = -t128 * t229 + (-t134 * t235 + t140 * t232) * t228;
t60 = -t127 * t229 + (-t133 * t235 + t139 * t232) * t228;
t58 = t59 * t229;
t57 = t59 * t278;
t56 = (t146 * t233 + t149 * t236) * t228 + t272;
t55 = (-t145 + t270) * t229 + t236 * t247;
t54 = t148 * t229 + t233 * t247 + t271;
t53 = t238 + t106;
t52 = -t105 + t237;
t51 = -t105 * t213 + t106 * t211;
t48 = (t145 * t233 + t148 * t236) * t228 + t249;
t47 = -t174 * t98 + t209 * t75;
t46 = t172 * t98 - t209 * t74;
t43 = (-t105 + t258) * t229 + t236 * t245;
t42 = t106 * t229 + t233 * t245 + t259;
t41 = t238 + t282;
t40 = -t285 + (-rSges(7,3) - pkin(10)) * t172 + t237 + t246;
t39 = -t100 * t213 - t102 * t174 + t104 * t175;
t38 = -t101 * t174 + t103 * t175 - t213 * t99;
t37 = -t100 * t211 - t102 * t172 + t104 * t173;
t36 = -t101 * t172 + t103 * t173 - t211 * t99;
t34 = t35 * t229;
t33 = t35 * t278;
t31 = -t172 * t75 + t174 * t74;
t30 = t213 * t281 + t278 * t282;
t29 = -t211 * t281 - t278 * t283;
t26 = (t105 * t233 + t106 * t236) * t228 + t240;
t25 = t211 * t282 - t213 * t283;
t22 = (t258 - t283) * t229 + t236 * t239;
t21 = t229 * t282 + t233 * t239 + t259;
t20 = t120 * t71 + t121 * t73 + t174 * t69;
t19 = t120 * t70 + t121 * t72 + t174 * t68;
t16 = (t233 * t283 + t236 * t282) * t228 + t240;
t15 = t58 + (t45 * t233 - t44 * t236) * t228;
t14 = -t44 * t211 - t45 * t213 + t57;
t13 = t50 * t229 + (t233 * t39 - t236 * t38) * t228;
t12 = t49 * t229 + (t233 * t37 - t236 * t36) * t228;
t11 = -t211 * t38 - t213 * t39 + t278 * t50;
t10 = -t211 * t36 - t213 * t37 + t278 * t49;
t9 = t34 + (-t23 * t236 + t24 * t233) * t228;
t8 = -t23 * t211 - t24 * t213 + t33;
t6 = t28 * t229 + (-t19 * t236 + t20 * t233) * t228;
t5 = t27 * t229 + (-t17 * t236 + t18 * t233) * t228;
t4 = -t19 * t211 - t20 * t213 + t278 * t28;
t3 = -t17 * t211 - t18 * t213 + t27 * t278;
t2 = t172 * t19 + t174 * t20 + t209 * t28;
t80 = [Icges(2,3) + m(7) * (t40 ^ 2 + t41 ^ 2) + m(6) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t90 ^ 2 + t91 ^ 2) + m(5) * (t78 ^ 2 + t79 ^ 2) + m(3) * (t113 ^ 2 + t114 ^ 2) + m(2) * (t218 ^ 2 + t219 ^ 2) + t59 + t35 + t301; t58 + t34 + m(7) * (t21 * t41 + t22 * t40) + m(6) * (t42 * t53 + t43 * t52) + m(5) * (t54 * t79 + t55 * t78) + m(4) * (t76 * t91 + t77 * t90) + m(3) * (t107 * t114 + t108 * t113) + ((-t62 / 0.2e1 - t60 / 0.2e1 - t82 / 0.2e1 - t83 / 0.2e1 - t84 / 0.2e1 - t64 / 0.2e1 + t243) * t236 + (t61 / 0.2e1 + t85 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1 + t65 / 0.2e1 + t63 / 0.2e1 - t242) * t233) * t228 + t265; (t9 + t15 + t265) * t229 + m(7) * (t16 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(6) * (t26 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t48 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t56 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2 + t88 ^ 2) + ((-t12 - t5 + (t261 * t227 * t236 + (t211 * t300 + t212 * t298) * t228) * t236 + (-t60 - t62 - t64 - t82 - t83 - t84) * t229) * t236 + (t6 + t13 + (t260 * t227 * t233 + (t213 * t299 + t214 * t297) * t228) * t233 + (t86 + t85 + t87 + t65 + t63 + t61) * t229 + ((t233 * t261 + t236 * t260) * t228 - t298 * t214 - t300 * t213 - t297 * t212 - t299 * t211) * t277) * t233) * t228; m(7) * (t211 * t41 + t213 * t40) + m(6) * (t211 * t53 + t213 * t52) + m(4) * (t211 * t91 + t213 * t90) + m(5) * (t211 * t79 + t213 * t78); m(7) * (-t16 * t278 + t21 * t211 + t213 * t22) + m(6) * (t211 * t42 + t213 * t43 - t26 * t278) + m(5) * (t211 * t54 + t213 * t55 - t278 * t48) + m(4) * (t211 * t76 + t213 * t77 - t278 * t56); 0.2e1 * (m(4) / 0.2e1 + t250) * (t227 * t235 ^ 2 + t211 ^ 2 + t213 ^ 2); ((-t233 * t40 + t236 * t41) * t293 + (-t233 * t52 + t236 * t53) * t294 + (-t233 * t78 + t236 * t79) * t295) * t296; m(7) * (-t229 * t16 + (t21 * t236 - t22 * t233) * t228) + m(6) * (-t229 * t26 + (-t233 * t43 + t236 * t42) * t228) + m(5) * (-t229 * t48 + (-t233 * t55 + t236 * t54) * t228); t250 * (t211 * t236 - t213 * t233 + t229 * t235) * t296; 0.2e1 * t250 * (t229 ^ 2 + (t233 ^ 2 + t236 ^ 2) * t227); t33 + t57 + m(7) * (t29 * t40 + t30 * t41) + m(6) * (t52 * t66 + t53 * t67) + t242 * t213 + t243 * t211; (t8 / 0.2e1 + t14 / 0.2e1) * t229 + (-t6 / 0.2e1 - t13 / 0.2e1) * t213 + (-t5 / 0.2e1 - t12 / 0.2e1) * t211 + m(7) * (t16 * t25 + t21 * t30 + t22 * t29) + m(6) * (t26 * t51 + t42 * t67 + t43 * t66) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t236 + (t9 / 0.2e1 + t15 / 0.2e1) * t235 + (t4 / 0.2e1 + t11 / 0.2e1) * t233) * t228; m(6) * (t211 * t67 + t213 * t66 - t278 * t51) + m(7) * (t211 * t30 + t213 * t29 - t25 * t278); m(6) * (-t51 * t229 + (-t233 * t66 + t236 * t67) * t228) + m(7) * (-t25 * t229 + (-t233 * t29 + t236 * t30) * t228); (t14 + t8) * t278 + (-t4 - t11) * t213 + (-t3 - t10) * t211 + m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t51 ^ 2 + t66 ^ 2 + t67 ^ 2); m(7) * (t40 * t46 + t41 * t47) + t32 + t263 * t174 + t264 * t172; m(7) * (t16 * t31 + t21 * t47 + t22 * t46) + t229 * t291 + t9 * t287 + t5 * t289 + t6 * t288 + (t233 * t2 / 0.2e1 + t236 * t292) * t228; m(7) * (t211 * t47 + t213 * t46 - t278 * t31); m(7) * (-t31 * t229 + (-t233 * t46 + t236 * t47) * t228); m(7) * (t25 * t31 + t29 * t46 + t30 * t47) + t4 * t288 - t213 * t2 / 0.2e1 + t8 * t287 + t3 * t289 + t211 * t292 + t278 * t291; t174 * t2 + t172 * t1 + t209 * t7 + m(7) * (t31 ^ 2 + t46 ^ 2 + t47 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t80(1) t80(2) t80(4) t80(7) t80(11) t80(16); t80(2) t80(3) t80(5) t80(8) t80(12) t80(17); t80(4) t80(5) t80(6) t80(9) t80(13) t80(18); t80(7) t80(8) t80(9) t80(10) t80(14) t80(19); t80(11) t80(12) t80(13) t80(14) t80(15) t80(20); t80(16) t80(17) t80(18) t80(19) t80(20) t80(21);];
Mq  = res;
