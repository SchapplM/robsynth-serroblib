% Calculate joint inertia matrix for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:53
% EndTime: 2019-03-09 09:59:01
% DurationCPUTime: 3.69s
% Computational Cost: add. (5336->452), mult. (8499->640), div. (0->0), fcn. (9016->8), ass. (0->214)
t299 = Icges(4,1) + Icges(3,3);
t200 = sin(qJ(2));
t203 = cos(qJ(2));
t298 = (-Icges(4,4) + Icges(3,5)) * t203 + (Icges(4,5) - Icges(3,6)) * t200;
t193 = qJ(4) + pkin(9);
t184 = sin(t193);
t204 = cos(qJ(1));
t185 = cos(t193);
t201 = sin(qJ(1));
t259 = t201 * t185;
t125 = t184 * t204 + t200 * t259;
t263 = t184 * t201;
t127 = -t185 * t204 + t200 * t263;
t291 = rSges(7,3) + qJ(6);
t292 = rSges(7,1) + pkin(5);
t297 = t291 * t125 - t292 * t127;
t296 = -t201 / 0.2e1;
t295 = t201 / 0.2e1;
t294 = -t204 / 0.2e1;
t293 = t204 / 0.2e1;
t107 = Icges(7,6) * t200 + (-Icges(7,5) * t184 + Icges(7,3) * t185) * t203;
t108 = Icges(6,3) * t200 + (-Icges(6,5) * t184 - Icges(6,6) * t185) * t203;
t109 = Icges(7,2) * t200 + (-Icges(7,4) * t184 + Icges(7,6) * t185) * t203;
t199 = sin(qJ(4));
t202 = cos(qJ(4));
t128 = Icges(5,3) * t200 + (-Icges(5,5) * t199 - Icges(5,6) * t202) * t203;
t261 = t185 * t203;
t290 = t107 * t261 + (t108 + t109 + t128) * t200;
t110 = Icges(6,6) * t200 + (-Icges(6,4) * t184 - Icges(6,2) * t185) * t203;
t111 = Icges(7,4) * t200 + (-Icges(7,1) * t184 + Icges(7,5) * t185) * t203;
t112 = Icges(6,5) * t200 + (-Icges(6,1) * t184 - Icges(6,4) * t185) * t203;
t131 = Icges(5,6) * t200 + (-Icges(5,4) * t199 - Icges(5,2) * t202) * t203;
t134 = Icges(5,5) * t200 + (-Icges(5,1) * t199 - Icges(5,4) * t202) * t203;
t289 = -t202 * t131 - t199 * t134 - t185 * t110 + (-t111 - t112) * t184;
t288 = t200 / 0.2e1;
t287 = -t298 * t201 + t299 * t204;
t286 = t299 * t201 + t298 * t204;
t260 = t200 * t204;
t123 = -t185 * t260 + t263;
t124 = t184 * t260 + t259;
t254 = t203 * t204;
t63 = Icges(7,5) * t124 + Icges(7,6) * t254 + Icges(7,3) * t123;
t67 = Icges(7,4) * t124 + Icges(7,2) * t254 + Icges(7,6) * t123;
t71 = Icges(7,1) * t124 + Icges(7,4) * t254 + Icges(7,5) * t123;
t19 = t123 * t63 + t124 * t71 + t254 * t67;
t256 = t201 * t203;
t64 = Icges(7,5) * t127 + Icges(7,6) * t256 - Icges(7,3) * t125;
t68 = Icges(7,4) * t127 + Icges(7,2) * t256 - Icges(7,6) * t125;
t72 = Icges(7,1) * t127 + Icges(7,4) * t256 - Icges(7,5) * t125;
t20 = t123 * t64 + t124 * t72 + t254 * t68;
t65 = Icges(6,5) * t124 - Icges(6,6) * t123 + Icges(6,3) * t254;
t69 = Icges(6,4) * t124 - Icges(6,2) * t123 + Icges(6,6) * t254;
t73 = Icges(6,1) * t124 - Icges(6,4) * t123 + Icges(6,5) * t254;
t21 = -t123 * t69 + t124 * t73 + t254 * t65;
t66 = Icges(6,5) * t127 + Icges(6,6) * t125 + Icges(6,3) * t256;
t70 = Icges(6,4) * t127 + Icges(6,2) * t125 + Icges(6,6) * t256;
t74 = Icges(6,1) * t127 + Icges(6,4) * t125 + Icges(6,5) * t256;
t22 = -t123 * t70 + t124 * t74 + t254 * t66;
t255 = t202 * t204;
t258 = t201 * t199;
t152 = t200 * t255 - t258;
t242 = t199 * t260;
t257 = t201 * t202;
t153 = t242 + t257;
t88 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t254;
t90 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t254;
t92 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t254;
t31 = t152 * t90 + t153 * t92 + t254 * t88;
t154 = t199 * t204 + t200 * t257;
t155 = t200 * t258 - t255;
t89 = Icges(5,5) * t155 + Icges(5,6) * t154 + Icges(5,3) * t256;
t91 = Icges(5,4) * t155 + Icges(5,2) * t154 + Icges(5,6) * t256;
t93 = Icges(5,1) * t155 + Icges(5,4) * t154 + Icges(5,5) * t256;
t32 = t152 * t91 + t153 * t93 + t254 * t89;
t42 = t123 * t107 + t109 * t254 + t124 * t111;
t43 = t108 * t254 - t123 * t110 + t124 * t112;
t48 = t128 * t254 + t152 * t131 + t153 * t134;
t285 = ((t19 + t21 + t31) * t204 + (t20 + t22 + t32) * t201) * t203 + (t42 + t43 + t48) * t200;
t23 = -t125 * t63 + t127 * t71 + t256 * t67;
t24 = -t125 * t64 + t127 * t72 + t256 * t68;
t25 = t125 * t69 + t127 * t73 + t256 * t65;
t26 = t125 * t70 + t127 * t74 + t256 * t66;
t33 = t154 * t90 + t155 * t92 + t256 * t88;
t34 = t154 * t91 + t155 * t93 + t256 * t89;
t44 = -t107 * t125 + t109 * t256 + t111 * t127;
t45 = t108 * t256 + t110 * t125 + t112 * t127;
t49 = t128 * t256 + t131 * t154 + t134 * t155;
t284 = ((t23 + t25 + t33) * t204 + (t24 + t26 + t34) * t201) * t203 + (t44 + t45 + t49) * t200;
t27 = t200 * t67 + (-t184 * t71 + t185 * t63) * t203;
t29 = t200 * t65 + (-t184 * t73 - t185 * t69) * t203;
t40 = t200 * t88 + (-t199 * t92 - t202 * t90) * t203;
t283 = t27 + t29 + t40;
t28 = t200 * t68 + (-t184 * t72 + t185 * t64) * t203;
t30 = t200 * t66 + (-t184 * t74 - t185 * t70) * t203;
t41 = t200 * t89 + (-t199 * t93 - t202 * t91) * t203;
t282 = t28 + t30 + t41;
t195 = t201 ^ 2;
t197 = t204 ^ 2;
t281 = 0.2e1 * t203;
t280 = m(4) / 0.2e1;
t279 = m(5) / 0.2e1;
t278 = m(6) / 0.2e1;
t277 = m(7) / 0.2e1;
t273 = pkin(4) * t199;
t272 = rSges(7,2) * t254 + t123 * t291 + t124 * t292;
t271 = rSges(7,2) * t256 - t297;
t270 = t204 * rSges(4,1);
t269 = t204 * rSges(3,3);
t268 = Icges(3,4) * t200;
t267 = Icges(3,4) * t203;
t266 = Icges(4,6) * t200;
t265 = Icges(4,6) * t203;
t264 = qJ(3) * t200;
t252 = t200 * rSges(7,2) + (-t184 * t292 + t185 * t291) * t203;
t248 = pkin(2) * t254 + qJ(3) * t260;
t251 = t195 * (pkin(2) * t203 + t264) + t204 * t248;
t166 = pkin(2) * t200 - qJ(3) * t203;
t250 = rSges(4,2) * t200 + rSges(4,3) * t203 - t166;
t183 = pkin(4) * t202 + pkin(3);
t198 = -qJ(5) - pkin(8);
t249 = -t204 * t183 - t198 * t256;
t247 = t204 * pkin(1) + t201 * pkin(7);
t246 = t201 * pkin(3) + pkin(8) * t254;
t245 = t195 + t197;
t244 = t278 + t277;
t243 = (t203 * t289 + t290) * t200;
t76 = t124 * rSges(6,1) - t123 * rSges(6,2) + rSges(6,3) * t254;
t94 = t153 * rSges(5,1) + t152 * rSges(5,2) + rSges(5,3) * t254;
t190 = t204 * pkin(7);
t241 = t190 - t249;
t240 = -Icges(4,4) * t200 / 0.2e1 + Icges(3,5) * t288 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t203;
t239 = -pkin(8) * t200 - t166;
t191 = t204 * pkin(3);
t238 = t201 * (pkin(8) * t256 - t191) + t204 * t246 + t251;
t237 = t247 + t248;
t144 = t200 * rSges(5,3) + (-rSges(5,1) * t199 - rSges(5,2) * t202) * t203;
t236 = -t144 + t239;
t147 = -t203 * t273 + (-pkin(8) - t198) * t200;
t235 = -t147 + t239;
t234 = rSges(3,1) * t203 - rSges(3,2) * t200;
t233 = -t155 * rSges(5,1) - t154 * rSges(5,2);
t232 = -t127 * rSges(6,1) - t125 * rSges(6,2);
t114 = t147 * t256;
t99 = t191 + (-pkin(8) * t203 + t200 * t273) * t201 + t249;
t16 = t114 + t252 * t256 + (-t99 - t271) * t200;
t212 = pkin(4) * t242 + t201 * t183 - t198 * t254;
t98 = t212 - t246;
t85 = t200 * t98;
t17 = t85 + t272 * t200 + (-t147 - t252) * t254;
t231 = t16 * t204 + t17 * t201;
t116 = t200 * rSges(6,3) + (-rSges(6,1) * t184 - rSges(6,2) * t185) * t203;
t78 = rSges(6,3) * t256 - t232;
t36 = t116 * t256 + t114 + (-t78 - t99) * t200;
t37 = t200 * t76 + t85 + (-t116 - t147) * t254;
t230 = t201 * t37 + t204 * t36;
t209 = t235 - t252;
t53 = t209 * t201;
t54 = t209 * t204;
t229 = t201 * t53 + t204 * t54;
t215 = -t116 + t235;
t61 = t215 * t201;
t62 = t215 * t204;
t228 = t201 * t61 + t204 * t62;
t227 = Icges(3,1) * t203 - t268;
t226 = -Icges(3,2) * t200 + t267;
t223 = -Icges(4,2) * t203 + t266;
t222 = Icges(4,3) * t200 - t265;
t221 = t123 * t204 - t125 * t201;
t214 = rSges(3,1) * t254 - rSges(3,2) * t260 + t201 * rSges(3,3);
t213 = t201 * rSges(4,1) - rSges(4,2) * t254 + rSges(4,3) * t260;
t211 = t201 * t99 + t204 * t98 + t238;
t210 = (-qJ(3) - t273) * t200 - pkin(1);
t208 = t40 / 0.2e1 + t29 / 0.2e1 + t27 / 0.2e1 + t48 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1;
t207 = t49 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 + t41 / 0.2e1 + t30 / 0.2e1 + t28 / 0.2e1;
t206 = t212 + t237;
t38 = ((-rSges(7,2) - pkin(2)) * t203 + t210) * t201 + t241 + t297;
t39 = t206 + t272;
t46 = ((-rSges(6,3) - pkin(2)) * t203 + t210) * t201 + t232 + t241;
t47 = t206 + t76;
t205 = (t201 * t39 + t204 * t38) * t277 + (t201 * t47 + t204 * t46) * t278;
t196 = t203 ^ 2;
t194 = t200 ^ 2;
t170 = rSges(2,1) * t204 - t201 * rSges(2,2);
t169 = -t201 * rSges(2,1) - rSges(2,2) * t204;
t168 = rSges(3,1) * t200 + rSges(3,2) * t203;
t106 = t250 * t204;
t105 = t250 * t201;
t102 = t214 + t247;
t101 = t269 + t190 + (-pkin(1) - t234) * t201;
t97 = t213 + t237;
t96 = t270 + t190 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t203 + (-rSges(4,3) - qJ(3)) * t200) * t201;
t95 = rSges(5,3) * t256 - t233;
t84 = t236 * t204;
t83 = t236 * t201;
t80 = t99 * t254;
t79 = t204 * t214 + (t201 * t234 - t269) * t201;
t60 = -t144 * t254 + t200 * t94;
t59 = t144 * t256 - t200 * t95;
t58 = t237 + t94 + t246;
t57 = t190 + t191 + (-t264 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(8)) * t203) * t201 + t233;
t55 = t204 * t213 + (-t270 + (-rSges(4,2) * t203 + rSges(4,3) * t200) * t201) * t201 + t251;
t52 = (-t201 * t94 + t204 * t95) * t203;
t35 = t201 * t95 + t204 * t94 + t238;
t18 = t80 + (t204 * t78 + (-t76 - t98) * t201) * t203;
t15 = t201 * t78 + t204 * t76 + t211;
t14 = t80 + (t271 * t204 + (-t98 - t272) * t201) * t203;
t13 = t33 * t201 - t204 * t34;
t12 = t31 * t201 - t204 * t32;
t11 = t201 * t271 + t204 * t272 + t211;
t10 = t25 * t201 - t204 * t26;
t9 = t23 * t201 - t204 * t24;
t8 = t21 * t201 - t204 * t22;
t7 = t19 * t201 - t20 * t204;
t1 = [Icges(2,3) + m(7) * (t38 ^ 2 + t39 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t96 ^ 2 + t97 ^ 2) + m(3) * (t101 ^ 2 + t102 ^ 2) + m(2) * (t169 ^ 2 + t170 ^ 2) + (t266 + t268 + (Icges(4,3) + Icges(3,2)) * t203 + t289) * t203 + (t265 + t267 + (Icges(3,1) + Icges(4,2)) * t200) * t200 + t290; (t240 * t204 + (Icges(4,5) * t294 + Icges(3,6) * t293 + t222 * t295 + t226 * t296) * t203 + (Icges(4,4) * t294 + Icges(3,5) * t293 + t223 * t295 + t227 * t296) * t200 - t207) * t204 + (t240 * t201 + (Icges(4,5) * t296 + Icges(3,6) * t295 + t222 * t294 + t226 * t293) * t203 + (Icges(4,4) * t296 + Icges(3,5) * t295 + t223 * t294 + t227 * t293) * t200 + t208) * t201 + m(7) * (t38 * t54 + t39 * t53) + m(6) * (t46 * t62 + t47 * t61) + m(5) * (t57 * t84 + t58 * t83) + m(4) * (t105 * t97 + t106 * t96) + m(3) * (-t101 * t204 - t102 * t201) * t168; m(7) * (t11 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t15 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t35 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t105 ^ 2 + t106 ^ 2 + t55 ^ 2) + m(3) * (t168 ^ 2 * t245 + t79 ^ 2) + (t287 * t197 - t10 - t13 - t9) * t204 + (t12 + t7 + t8 + t286 * t195 + (t287 * t201 + t286 * t204) * t204) * t201; 0.2e1 * ((t201 * t58 + t204 * t57) * t279 + (t201 * t97 + t204 * t96) * t280 + t205) * t200; m(7) * (-t203 * t11 + t200 * t229) + m(6) * (-t203 * t15 + t200 * t228) + m(5) * (-t203 * t35 + (t201 * t83 + t204 * t84) * t200) + m(4) * (-t203 * t55 + (t105 * t201 + t106 * t204) * t200); 0.2e1 * (t280 + t279 + t244) * (t194 * t245 + t196); m(7) * (t16 * t38 + t17 * t39) + m(6) * (t36 * t46 + t37 * t47) + m(5) * (t57 * t59 + t58 * t60) + (t201 * t207 + t204 * t208) * t203 + t243; m(7) * (t11 * t14 + t16 * t54 + t17 * t53) + m(6) * (t15 * t18 + t36 * t62 + t37 * t61) + m(5) * (t35 * t52 + t59 * t84 + t60 * t83) + ((t12 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1) * t204 + (t13 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t201) * t203 + (t283 * t201 - t282 * t204) * t288 + t285 * t295 + t284 * t294; m(5) * (-t52 * t203 + (t201 * t60 + t204 * t59) * t200) + m(6) * (-t18 * t203 + t200 * t230) + m(7) * (-t14 * t203 + t200 * t231); t243 * t200 + m(7) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(6) * (t18 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t52 ^ 2 + t59 ^ 2 + t60 ^ 2) + (t285 * t204 + t284 * t201 + (t282 * t201 + t283 * t204) * t200) * t203; t205 * t281; m(7) * (t200 * t11 + t203 * t229) + m(6) * (t200 * t15 + t203 * t228); t244 * (-0.1e1 + t245) * t200 * t281; m(7) * (t200 * t14 + t203 * t231) + m(6) * (t200 * t18 + t203 * t230); 0.2e1 * t244 * (t196 * t245 + t194); m(7) * (t123 * t38 - t125 * t39); m(7) * (t11 * t261 + t123 * t54 - t125 * t53); m(7) * (-t196 * t185 + t200 * t221); m(7) * (t123 * t16 - t125 * t17 + t14 * t261); m(7) * (t185 * t200 + t221) * t203; m(7) * (t185 ^ 2 * t196 + t123 ^ 2 + t125 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
