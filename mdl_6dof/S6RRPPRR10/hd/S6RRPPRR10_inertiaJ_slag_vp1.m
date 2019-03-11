% Calculate joint inertia matrix for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:11
% EndTime: 2019-03-09 09:34:20
% DurationCPUTime: 3.79s
% Computational Cost: add. (7274->461), mult. (8573->662), div. (0->0), fcn. (9046->10), ass. (0->224)
t299 = Icges(4,1) + Icges(3,3);
t204 = sin(qJ(2));
t206 = cos(qJ(2));
t298 = (-Icges(4,4) + Icges(3,5)) * t206 + (Icges(4,5) - Icges(3,6)) * t204;
t205 = sin(qJ(1));
t297 = -t205 / 0.2e1;
t287 = t205 / 0.2e1;
t207 = cos(qJ(1));
t286 = -t207 / 0.2e1;
t296 = t207 / 0.2e1;
t288 = t204 / 0.2e1;
t295 = t205 * t299 + t298 * t207;
t294 = -t298 * t205 + t207 * t299;
t198 = t205 ^ 2;
t200 = t207 ^ 2;
t293 = 0.2e1 * t206;
t292 = m(4) / 0.2e1;
t291 = m(5) / 0.2e1;
t290 = m(6) / 0.2e1;
t289 = m(7) / 0.2e1;
t167 = rSges(3,1) * t204 + rSges(3,2) * t206;
t285 = m(3) * t167;
t201 = sin(pkin(10));
t284 = pkin(4) * t201;
t203 = -pkin(8) - qJ(4);
t202 = cos(pkin(10));
t183 = t202 * pkin(4) + pkin(3);
t195 = -pkin(9) + t203;
t255 = -t195 + t203;
t272 = t204 * t207;
t250 = t201 * t272;
t261 = -pkin(4) * t250 - t205 * t183;
t196 = pkin(10) + qJ(5);
t185 = cos(t196);
t155 = pkin(5) * t185 + t183;
t184 = sin(t196);
t156 = pkin(5) * t184 + t284;
t262 = t205 * t155 + t156 * t272;
t264 = t206 * t207;
t60 = t255 * t264 + t261 + t262;
t186 = qJ(6) + t196;
t182 = cos(t186);
t181 = sin(t186);
t271 = t205 * t181;
t114 = t182 * t272 - t271;
t270 = t205 * t182;
t115 = t181 * t272 + t270;
t72 = t115 * rSges(7,1) + t114 * rSges(7,2) + rSges(7,3) * t264;
t283 = t60 + t72;
t245 = -t156 + t284;
t265 = t205 * t206;
t259 = t207 * t183 + t203 * t265;
t273 = t155 * t207;
t61 = -t273 + (-t195 * t206 - t204 * t245) * t205 + t259;
t116 = t181 * t207 + t204 * t270;
t117 = -t182 * t207 + t204 * t271;
t233 = -t117 * rSges(7,1) - t116 * rSges(7,2);
t73 = rSges(7,3) * t265 - t233;
t282 = t61 + t73;
t281 = t207 * rSges(4,1);
t280 = t207 * rSges(3,3);
t106 = rSges(7,3) * t204 + (-rSges(7,1) * t181 - rSges(7,2) * t182) * t206;
t95 = t204 * t255 + t206 * t245;
t279 = -t106 - t95;
t278 = Icges(3,4) * t204;
t277 = Icges(3,4) * t206;
t276 = Icges(4,6) * t204;
t275 = Icges(4,6) * t206;
t274 = qJ(3) * t204;
t269 = t205 * t184;
t268 = t205 * t185;
t267 = t205 * t201;
t266 = t205 * t202;
t258 = pkin(2) * t264 + qJ(3) * t272;
t263 = t198 * (pkin(2) * t206 + t274) + t207 * t258;
t165 = pkin(2) * t204 - qJ(3) * t206;
t260 = rSges(4,2) * t204 + rSges(4,3) * t206 - t165;
t257 = t207 * pkin(1) + t205 * pkin(7);
t256 = t205 * pkin(3) + qJ(4) * t264;
t254 = t198 + t200;
t66 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t264;
t68 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t264;
t70 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t264;
t32 = t204 * t66 + (-t181 * t70 - t182 * t68) * t206;
t67 = Icges(7,5) * t117 + Icges(7,6) * t116 + Icges(7,3) * t265;
t69 = Icges(7,4) * t117 + Icges(7,2) * t116 + Icges(7,6) * t265;
t71 = Icges(7,1) * t117 + Icges(7,4) * t116 + Icges(7,5) * t265;
t33 = t204 * t67 + (-t181 * t71 - t182 * t69) * t206;
t102 = Icges(7,6) * t204 + (-Icges(7,4) * t181 - Icges(7,2) * t182) * t206;
t103 = Icges(7,5) * t204 + (-Icges(7,1) * t181 - Icges(7,4) * t182) * t206;
t220 = -t102 * t182 - t103 * t181;
t101 = Icges(7,3) * t204 + (-Icges(7,5) * t181 - Icges(7,6) * t182) * t206;
t97 = t204 * t101;
t49 = (t206 * t220 + t97) * t204;
t20 = t114 * t68 + t115 * t70 + t264 * t66;
t21 = t114 * t69 + t115 * t71 + t264 * t67;
t39 = t101 * t264 + t114 * t102 + t115 * t103;
t5 = t39 * t204 + (t20 * t207 + t205 * t21) * t206;
t22 = t116 * t68 + t117 * t70 + t265 * t66;
t23 = t116 * t69 + t117 * t71 + t265 * t67;
t40 = t101 * t265 + t102 * t116 + t103 * t117;
t6 = t40 * t204 + (t205 * t23 + t207 * t22) * t206;
t253 = t5 * t264 + t6 * t265 + t204 * (t49 + (t205 * t33 + t207 * t32) * t206);
t124 = t185 * t272 - t269;
t125 = t184 * t272 + t268;
t74 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t264;
t76 = Icges(6,4) * t125 + Icges(6,2) * t124 + Icges(6,6) * t264;
t78 = Icges(6,1) * t125 + Icges(6,4) * t124 + Icges(6,5) * t264;
t34 = t204 * t74 + (-t184 * t78 - t185 * t76) * t206;
t107 = Icges(6,3) * t204 + (-Icges(6,5) * t184 - Icges(6,6) * t185) * t206;
t108 = Icges(6,6) * t204 + (-Icges(6,4) * t184 - Icges(6,2) * t185) * t206;
t109 = Icges(6,5) * t204 + (-Icges(6,1) * t184 - Icges(6,4) * t185) * t206;
t44 = t107 * t264 + t124 * t108 + t125 * t109;
t252 = t44 / 0.2e1 + t34 / 0.2e1;
t126 = t184 * t207 + t204 * t268;
t127 = -t185 * t207 + t204 * t269;
t75 = Icges(6,5) * t127 + Icges(6,6) * t126 + Icges(6,3) * t265;
t77 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t265;
t79 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t265;
t35 = t204 * t75 + (-t184 * t79 - t185 * t77) * t206;
t45 = t107 * t265 + t108 * t126 + t109 * t127;
t251 = t45 / 0.2e1 + t35 / 0.2e1;
t80 = t125 * rSges(6,1) + t124 * rSges(6,2) + rSges(6,3) * t264;
t148 = t202 * t272 - t267;
t149 = t250 + t266;
t249 = t149 * rSges(5,1) + t148 * rSges(5,2) + rSges(5,3) * t264;
t248 = t265 / 0.2e1;
t247 = t264 / 0.2e1;
t246 = -Icges(4,4) * t204 / 0.2e1 + Icges(3,5) * t288 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t206;
t244 = -qJ(4) * t204 - t165;
t243 = t291 + t290 + t289;
t193 = t207 * pkin(3);
t242 = t205 * (qJ(4) * t265 - t193) + t207 * t256 + t263;
t241 = t257 + t258;
t12 = t20 * t205 - t207 * t21;
t13 = t22 * t205 - t207 * t23;
t240 = t12 * t247 + t13 * t248 + t6 * t286 + t5 * t287 + (t32 * t205 - t33 * t207) * t288;
t239 = t49 + (t33 + t40) * t248 + (t32 + t39) * t247;
t238 = -rSges(5,3) * t204 - (-rSges(5,1) * t201 - rSges(5,2) * t202) * t206 + t244;
t237 = t206 * t284 - (-qJ(4) - t203) * t204 + t244;
t236 = rSges(3,1) * t206 - rSges(3,2) * t204;
t150 = t201 * t207 + t204 * t266;
t151 = -t202 * t207 + t204 * t267;
t235 = -rSges(5,1) * t151 - rSges(5,2) * t150;
t234 = -rSges(6,1) * t127 - rSges(6,2) * t126;
t96 = t106 * t265;
t24 = -t204 * t282 + t265 * t95 + t96;
t65 = t204 * t72;
t25 = t204 * t60 + t264 * t279 + t65;
t232 = t205 * t25 + t207 * t24;
t209 = t237 + t279;
t47 = t209 * t205;
t48 = t209 * t207;
t231 = t205 * t47 + t207 * t48;
t53 = -t204 * t73 + t96;
t54 = -t106 * t264 + t65;
t230 = t205 * t54 + t207 * t53;
t113 = rSges(6,3) * t204 + (-rSges(6,1) * t184 - rSges(6,2) * t185) * t206;
t81 = rSges(6,3) * t265 - t234;
t55 = t113 * t265 - t204 * t81;
t56 = -t113 * t264 + t204 * t80;
t229 = t205 * t56 + t207 * t55;
t211 = -t113 + t237;
t62 = t211 * t205;
t63 = t211 * t207;
t228 = t205 * t62 + t207 * t63;
t83 = t238 * t205;
t84 = t238 * t207;
t227 = t205 * t83 + t207 * t84;
t226 = Icges(3,1) * t206 - t278;
t225 = -Icges(3,2) * t204 + t277;
t222 = -Icges(4,2) * t206 + t276;
t221 = Icges(4,3) * t204 - t275;
t219 = -t108 * t185 - t109 * t184;
t214 = rSges(3,1) * t264 - rSges(3,2) * t272 + t205 * rSges(3,3);
t213 = t205 * rSges(4,1) - rSges(4,2) * t264 + rSges(4,3) * t272;
t212 = -t203 * t264 - t261;
t210 = t205 * (t193 + (-qJ(4) * t206 + t204 * t284) * t205 - t259) + t207 * (t212 - t256) + t242;
t192 = t207 * pkin(7);
t42 = t273 + t192 + (-pkin(1) + (-qJ(3) - t156) * t204 + (-rSges(7,3) - pkin(2) + t195) * t206) * t205 + t233;
t43 = -t195 * t264 + t241 + t262 + t72;
t50 = t192 + (-pkin(1) + (-rSges(6,3) - pkin(2)) * t206 + (-qJ(3) - t284) * t204) * t205 + t234 + t259;
t51 = t212 + t241 + t80;
t58 = t192 + t193 + (-t274 - pkin(1) + (-rSges(5,3) - pkin(2) - qJ(4)) * t206) * t205 + t235;
t59 = t241 + t249 + t256;
t208 = (t205 * t43 + t207 * t42) * t289 + (t205 * t51 + t207 * t50) * t290 + (t205 * t59 + t207 * t58) * t291;
t199 = t206 ^ 2;
t197 = t204 ^ 2;
t169 = rSges(2,1) * t207 - t205 * rSges(2,2);
t168 = -t205 * rSges(2,1) - rSges(2,2) * t207;
t122 = Icges(5,5) * t204 + (-Icges(5,1) * t201 - Icges(5,4) * t202) * t206;
t121 = Icges(5,6) * t204 + (-Icges(5,4) * t201 - Icges(5,2) * t202) * t206;
t105 = t260 * t207;
t104 = t260 * t205;
t100 = t204 * t107;
t99 = t214 + t257;
t98 = t280 + t192 + (-pkin(1) - t236) * t205;
t94 = t213 + t241;
t93 = t281 + t192 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t206 + (-rSges(4,3) - qJ(3)) * t204) * t205;
t90 = Icges(5,1) * t151 + Icges(5,4) * t150 + Icges(5,5) * t265;
t89 = Icges(5,1) * t149 + Icges(5,4) * t148 + Icges(5,5) * t264;
t88 = Icges(5,4) * t151 + Icges(5,2) * t150 + Icges(5,6) * t265;
t87 = Icges(5,4) * t149 + Icges(5,2) * t148 + Icges(5,6) * t264;
t86 = Icges(5,5) * t151 + Icges(5,6) * t150 + Icges(5,3) * t265;
t85 = Icges(5,5) * t149 + Icges(5,6) * t148 + Icges(5,3) * t264;
t82 = t207 * t214 + (t205 * t236 - t280) * t205;
t64 = t73 * t264;
t57 = t207 * t213 + (-t281 + (-rSges(4,2) * t206 + rSges(4,3) * t204) * t205) * t205 + t263;
t52 = (t206 * t219 + t100) * t204;
t46 = (-t205 * t80 + t207 * t81) * t206;
t41 = -t265 * t72 + t64;
t36 = t205 * (rSges(5,3) * t265 - t235) + t207 * t249 + t242;
t29 = t126 * t77 + t127 * t79 + t265 * t75;
t28 = t126 * t76 + t127 * t78 + t265 * t74;
t27 = t124 * t77 + t125 * t79 + t264 * t75;
t26 = t124 * t76 + t125 * t78 + t264 * t74;
t19 = t205 * t81 + t207 * t80 + t210;
t18 = t64 + (-t205 * t283 + t207 * t61) * t206;
t16 = t28 * t205 - t207 * t29;
t15 = t26 * t205 - t207 * t27;
t14 = t205 * t282 + t207 * t283 + t210;
t8 = t45 * t204 + (t205 * t29 + t207 * t28) * t206;
t7 = t44 * t204 + (t205 * t27 + t207 * t26) * t206;
t1 = [Icges(2,3) + t100 + t97 + m(7) * (t42 ^ 2 + t43 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t93 ^ 2 + t94 ^ 2) + m(3) * (t98 ^ 2 + t99 ^ 2) + m(2) * (t168 ^ 2 + t169 ^ 2) + (-t121 * t202 - t122 * t201 + t219 + t220 + t276 + t278 + (Icges(3,2) + Icges(4,3)) * t206) * t206 + (t277 + t275 + (-Icges(5,5) * t201 - Icges(5,6) * t202) * t206 + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t204) * t204; m(4) * (t104 * t94 + t105 * t93) + m(5) * (t58 * t84 + t59 * t83) + m(6) * (t50 * t63 + t51 * t62) + m(7) * (t42 * t48 + t43 * t47) + (-t150 * t121 / 0.2e1 - t151 * t122 / 0.2e1 - t98 * t285 - t33 / 0.2e1 - t40 / 0.2e1 + t246 * t207 + (Icges(3,5) * t296 + t226 * t297 + Icges(4,4) * t286 + t222 * t287 - t86 / 0.2e1) * t204 - t251) * t207 + (t148 * t121 / 0.2e1 + t149 * t122 / 0.2e1 - t99 * t285 + t32 / 0.2e1 + t39 / 0.2e1 + t246 * t205 + (t85 / 0.2e1 + Icges(3,5) * t287 + t226 * t296 + Icges(4,4) * t297 + t222 * t286) * t204 + t252) * t205 + ((Icges(3,6) * t296 + t225 * t297 + Icges(4,5) * t286 + t221 * t287 + t201 * t90 / 0.2e1 + t202 * t88 / 0.2e1) * t207 + (Icges(3,6) * t287 + t225 * t296 + Icges(4,5) * t297 + t221 * t286 - t201 * t89 / 0.2e1 - t202 * t87 / 0.2e1) * t205) * t206; m(7) * (t14 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t19 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(5) * (t36 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t104 ^ 2 + t105 ^ 2 + t57 ^ 2) + m(3) * (t167 ^ 2 * t254 + t82 ^ 2) + (-t13 - t16 + (t150 * t88 + t151 * t90 + t86 * t265) * t207 + t294 * t200) * t207 + (t12 + t15 + (t148 * t87 + t149 * t89 + t85 * t264) * t205 + t295 * t198 + (-t148 * t88 - t149 * t90 - t150 * t87 - t151 * t89 + t294 * t205 + t295 * t207 - t86 * t264 - t85 * t265) * t207) * t205; 0.2e1 * ((t205 * t94 + t207 * t93) * t292 + t208) * t204; m(7) * (-t206 * t14 + t204 * t231) + m(6) * (-t206 * t19 + t204 * t228) + m(5) * (t204 * t227 - t206 * t36) + m(4) * (-t206 * t57 + (t104 * t205 + t105 * t207) * t204); 0.2e1 * (t292 + t243) * (t197 * t254 + t199); t208 * t293; m(7) * (t204 * t14 + t206 * t231) + m(6) * (t204 * t19 + t206 * t228) + m(5) * (t204 * t36 + t206 * t227); t243 * (-0.1e1 + t254) * t204 * t293; 0.2e1 * t243 * (t199 * t254 + t197); t52 + m(7) * (t24 * t42 + t25 * t43) + m(6) * (t50 * t55 + t51 * t56) + (t205 * t251 + t207 * t252) * t206 + t239; t8 * t286 + t7 * t287 + (t34 * t205 - t35 * t207) * t288 + (t15 * t296 + t16 * t287) * t206 + m(7) * (t14 * t18 + t24 * t48 + t25 * t47) + m(6) * (t19 * t46 + t55 * t63 + t56 * t62) + t240; m(6) * (t204 * t229 - t46 * t206) + m(7) * (-t18 * t206 + t204 * t232); m(6) * (t46 * t204 + t206 * t229) + m(7) * (t18 * t204 + t206 * t232); t204 * t52 + m(7) * (t18 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t46 ^ 2 + t55 ^ 2 + t56 ^ 2) + (t207 * t7 + t205 * t8 + t204 * (t205 * t35 + t207 * t34)) * t206 + t253; m(7) * (t42 * t53 + t43 * t54) + t239; m(7) * (t14 * t41 + t47 * t54 + t48 * t53) + t240; m(7) * (t204 * t230 - t41 * t206); m(7) * (t41 * t204 + t206 * t230); m(7) * (t18 * t41 + t24 * t53 + t25 * t54) + t253; m(7) * (t41 ^ 2 + t53 ^ 2 + t54 ^ 2) + t253;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
