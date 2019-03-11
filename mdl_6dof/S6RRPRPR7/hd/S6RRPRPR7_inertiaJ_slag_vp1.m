% Calculate joint inertia matrix for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:21
% EndTime: 2019-03-09 10:45:31
% DurationCPUTime: 4.15s
% Computational Cost: add. (7618->415), mult. (11383->607), div. (0->0), fcn. (13192->10), ass. (0->209)
t307 = Icges(3,1) + Icges(4,1);
t304 = Icges(3,5) + Icges(4,4);
t302 = Icges(6,3) + Icges(5,3);
t183 = sin(qJ(2));
t306 = (Icges(3,4) - Icges(4,5)) * t183;
t239 = qJ(4) + pkin(10);
t171 = sin(t239);
t188 = cos(qJ(1));
t227 = cos(t239);
t217 = t183 * t227;
t187 = cos(qJ(2));
t246 = t187 * t188;
t115 = t171 * t246 - t188 * t217;
t138 = t183 * t171 + t187 * t227;
t116 = t138 * t188;
t186 = cos(qJ(4));
t182 = sin(qJ(4));
t249 = t182 * t187;
t144 = t183 * t186 - t249;
t136 = t144 * t188;
t250 = t182 * t183;
t204 = t186 * t187 + t250;
t137 = t204 * t188;
t184 = sin(qJ(1));
t305 = -Icges(5,5) * t137 - Icges(6,5) * t116 - Icges(5,6) * t136 + Icges(6,6) * t115 + t302 * t184;
t303 = Icges(4,2) + Icges(3,3);
t301 = t304 * t187 + (-Icges(3,6) + Icges(4,6)) * t183;
t300 = t307 * t187 - t306;
t299 = -t184 / 0.2e1;
t292 = t184 / 0.2e1;
t298 = -t188 / 0.2e1;
t297 = t188 / 0.2e1;
t139 = -t187 * t171 + t217;
t296 = Icges(5,5) * t144 + Icges(6,5) * t139 - Icges(5,6) * t204 - Icges(6,6) * t138;
t295 = t115 / 0.2e1;
t294 = t138 / 0.2e1;
t290 = -t301 * t184 + t303 * t188;
t289 = t303 * t184 + t301 * t188;
t113 = t139 * t184;
t114 = t138 * t184;
t134 = t144 * t184;
t135 = t204 * t184;
t179 = t188 ^ 2;
t254 = Icges(6,6) * t188;
t255 = Icges(5,6) * t188;
t256 = Icges(6,2) * t113;
t257 = Icges(5,2) * t134;
t258 = Icges(6,5) * t188;
t259 = Icges(5,5) * t188;
t262 = Icges(6,4) * t114;
t263 = Icges(5,4) * t135;
t266 = Icges(6,1) * t114;
t267 = Icges(5,1) * t135;
t181 = sin(qJ(6));
t185 = cos(qJ(6));
t96 = -t116 * t181 - t184 * t185;
t97 = t116 * t185 - t181 * t184;
t48 = Icges(7,5) * t97 + Icges(7,6) * t96 + Icges(7,3) * t115;
t49 = Icges(7,4) * t97 + Icges(7,2) * t96 + Icges(7,6) * t115;
t50 = Icges(7,1) * t97 + Icges(7,4) * t96 + Icges(7,5) * t115;
t94 = -t114 * t181 + t185 * t188;
t95 = t114 * t185 + t181 * t188;
t12 = -t113 * t48 + t49 * t94 + t50 * t95;
t252 = Icges(7,3) * t113;
t273 = Icges(7,5) * t95;
t216 = -t252 + 0.2e1 * t273;
t272 = Icges(7,2) * t94;
t274 = Icges(7,4) * t95;
t219 = t272 + 0.2e1 * t274;
t271 = Icges(7,6) * t94;
t275 = Icges(7,1) * t95 ^ 2;
t3 = t12 * t184 - (t275 + t219 * t94 - (t216 + 0.2e1 * t271) * t113) * t188;
t67 = Icges(6,4) * t116 - Icges(6,2) * t115 - Icges(6,6) * t184;
t68 = Icges(6,1) * t116 - Icges(6,4) * t115 - Icges(6,5) * t184;
t77 = Icges(5,4) * t137 + Icges(5,2) * t136 - Icges(5,6) * t184;
t78 = Icges(5,1) * t137 + Icges(5,4) * t136 - Icges(5,5) * t184;
t238 = -t3 + ((0.2e1 * t258 + t266) * t114 - (-0.2e1 * t254 - t256 - 0.2e1 * t262) * t113 + (0.2e1 * t259 + t267) * t135 + (0.2e1 * t255 + t257 + 0.2e1 * t263) * t134 + t302 * t179) * t188 + (-t113 * t67 - t114 * t68 - t134 * t77 - t135 * t78 + t305 * t188) * t184;
t194 = t254 + t256 + t262;
t195 = t255 + t257 + t263;
t196 = Icges(6,4) * t113 + t258 + t266;
t197 = Icges(5,4) * t134 + t259 + t267;
t13 = t115 * t48 + t49 * t96 + t50 * t97;
t198 = -t252 + t271 + t273;
t253 = Icges(7,6) * t113;
t199 = -t253 + t272 + t274;
t200 = Icges(7,1) * t95 + Icges(7,4) * t94 - Icges(7,5) * t113;
t189 = t115 * t198 + t199 * t96 + t200 * t97;
t4 = t13 * t184 - t188 * t189;
t237 = t4 + (t115 * t194 - t116 * t196 - t136 * t195 - t137 * t197) * t188 + (-t115 * t67 + t116 * t68 + t136 * t77 + t137 * t78 + t305 * t184 + (Icges(5,5) * t135 + Icges(6,5) * t114 + Icges(5,6) * t134 + Icges(6,6) * t113 + t302 * t188) * t188) * t184;
t288 = -t237 * t184 - t238 * t188;
t180 = -qJ(5) - pkin(8);
t169 = t188 * t180;
t166 = pkin(3) * t246;
t228 = -pkin(8) * t184 + t166;
t170 = pkin(4) * t186 + pkin(3);
t235 = pkin(4) * t250;
t232 = t170 * t246 + t184 * t180 + t188 * t235;
t281 = -pkin(3) + t170;
t282 = pkin(8) * t188;
t276 = t184 * (-t282 - t169 + (t187 * t281 + t235) * t184) + t188 * (-t228 + t232);
t52 = t97 * rSges(7,1) + t96 * rSges(7,2) + t115 * rSges(7,3);
t279 = t116 * pkin(5) + pkin(9) * t115 + t52;
t283 = pkin(5) * t114;
t222 = -rSges(7,1) * t95 - rSges(7,2) * t94;
t51 = -rSges(7,3) * t113 - t222;
t9 = -t184 * (-pkin(9) * t113 + t283 + t51) - t188 * t279 - t276;
t178 = t184 ^ 2;
t287 = m(4) / 0.2e1;
t286 = m(6) / 0.2e1;
t285 = m(7) / 0.2e1;
t284 = -rSges(5,3) - pkin(8);
t62 = Icges(7,3) * t138 + (Icges(7,5) * t185 - Icges(7,6) * t181) * t139;
t64 = Icges(7,5) * t138 + (Icges(7,1) * t185 - Icges(7,4) * t181) * t139;
t278 = t139 * t185 * t64 + t138 * t62;
t65 = rSges(7,3) * t138 + (rSges(7,1) * t185 - rSges(7,2) * t181) * t139;
t277 = pkin(5) * t139 + pkin(9) * t138 + t65;
t63 = Icges(7,6) * t138 + (Icges(7,4) * t185 - Icges(7,2) * t181) * t139;
t270 = t181 * t63;
t269 = t188 * rSges(4,2);
t268 = t188 * rSges(3,3);
t264 = Icges(3,4) * t187;
t260 = Icges(4,5) * t187;
t251 = qJ(3) * t183;
t248 = t183 * t188;
t245 = t137 * rSges(5,1) + t136 * rSges(5,2);
t242 = pkin(2) * t246 + qJ(3) * t248;
t244 = t178 * (pkin(2) * t187 + t251) + t188 * t242;
t155 = pkin(2) * t183 - qJ(3) * t187;
t243 = -rSges(4,1) * t183 + rSges(4,3) * t187 - t155;
t241 = t188 * pkin(1) + t184 * pkin(7);
t240 = t178 + t179;
t236 = t286 + t285;
t16 = t138 * t198 + (-t181 * t199 + t185 * t200) * t139;
t18 = -t113 * t62 + t63 * t94 + t64 * t95;
t234 = t18 / 0.2e1 + t16 / 0.2e1;
t17 = t138 * t48 + (-t181 * t49 + t185 * t50) * t139;
t19 = t115 * t62 + t63 * t96 + t64 * t97;
t233 = t19 / 0.2e1 + t17 / 0.2e1;
t231 = rSges(4,1) * t246 + t184 * rSges(4,2) + rSges(4,3) * t248;
t230 = t304 * t183 / 0.2e1 + (Icges(3,6) / 0.2e1 - Icges(4,6) / 0.2e1) * t187;
t229 = -pkin(3) * t183 - t155;
t226 = t184 * (t184 * t187 * pkin(3) + t282) + t188 * t228 + t244;
t225 = t241 + t242;
t90 = rSges(6,1) * t139 - rSges(6,2) * t138;
t224 = t229 - t90;
t102 = rSges(5,1) * t144 - rSges(5,2) * t204;
t223 = -t102 + t229;
t70 = t116 * rSges(6,1) - t115 * rSges(6,2) - rSges(6,3) * t184;
t221 = rSges(3,1) * t187 - rSges(3,2) * t183;
t220 = -t135 * rSges(5,1) - t134 * rSges(5,2);
t71 = t223 * t184;
t72 = t223 * t188;
t218 = t184 * t71 + t188 * t72;
t47 = t184 * t220 - t188 * t245;
t215 = t229 - t277;
t212 = -Icges(3,2) * t183 + t264;
t209 = Icges(4,3) * t183 + t260;
t203 = rSges(3,1) * t246 - rSges(3,2) * t248 + t184 * rSges(3,3);
t175 = t188 * pkin(7);
t58 = t175 + t284 * t188 + (-t251 - pkin(1) + (-pkin(2) - pkin(3)) * t187) * t184 + t220;
t59 = t184 * t284 + t166 + t225 + t245;
t202 = m(5) * (t184 * t59 + t188 * t58);
t69 = t114 * rSges(6,1) + t113 * rSges(6,2) + rSges(6,3) * t188;
t25 = -t184 * t69 - t188 * t70 - t276;
t201 = t225 + t232;
t100 = Icges(5,4) * t144 - Icges(5,2) * t204;
t101 = Icges(5,1) * t144 - Icges(5,4) * t204;
t88 = Icges(6,4) * t139 - Icges(6,2) * t138;
t89 = Icges(6,1) * t139 - Icges(6,4) * t138;
t193 = -t138 * t194 / 0.2e1 + t139 * t196 / 0.2e1 + t134 * t100 / 0.2e1 + t135 * t101 / 0.2e1 + t113 * t88 / 0.2e1 + t114 * t89 / 0.2e1 - t204 * t195 / 0.2e1 + t144 * t197 / 0.2e1 + t234 + t296 * t297;
t192 = -t100 * t136 / 0.2e1 - t101 * t137 / 0.2e1 + t88 * t295 - t116 * t89 / 0.2e1 + t204 * t77 / 0.2e1 - t144 * t78 / 0.2e1 + t67 * t294 - t139 * t68 / 0.2e1 - t233 + t296 * t292;
t191 = t169 + t175 + (-pkin(1) + (-pkin(2) - t170) * t187 + (-pkin(4) * t182 - qJ(3)) * t183) * t184;
t1 = t12 * t115 + t18 * t138 - (t275 - t216 * t113 + (t219 - 0.2e1 * t253) * t94) * t113;
t2 = -t113 * t189 + t13 * t115 + t19 * t138;
t190 = -t113 * t3 / 0.2e1 + t4 * t295 + (-t16 * t188 + t17 * t184) * t294 + t2 * t292 + t1 * t298;
t159 = rSges(2,1) * t188 - t184 * rSges(2,2);
t158 = -t184 * rSges(2,1) - rSges(2,2) * t188;
t157 = rSges(3,1) * t183 + rSges(3,2) * t187;
t117 = -pkin(4) * t249 + t183 * t281;
t108 = t188 * t117;
t107 = t184 * t117;
t106 = t243 * t188;
t105 = t243 * t184;
t104 = t203 + t241;
t103 = t268 + t175 + (-pkin(1) - t221) * t184;
t84 = t225 + t231;
t83 = t269 + t175 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t187 + (-rSges(4,3) - qJ(3)) * t183) * t184;
t75 = t188 * t203 + (t184 * t221 - t268) * t184;
t61 = t188 * t90 + t108;
t60 = t184 * t90 + t107;
t56 = t188 * t231 + (-t269 + (rSges(4,1) * t187 + rSges(4,3) * t183) * t184) * t184 + t244;
t54 = t188 * t224 - t108;
t53 = t184 * t224 - t107;
t46 = t201 + t70;
t45 = t191 - t69;
t40 = t188 * t277 + t108;
t39 = t184 * t277 + t107;
t36 = t188 * t215 - t108;
t35 = t184 * t215 - t107;
t34 = -t47 + t226;
t29 = t201 + t279;
t28 = -t283 - (-rSges(7,3) - pkin(9)) * t113 + t191 + t222;
t27 = -t115 * t65 + t138 * t52;
t26 = -t113 * t65 - t138 * t51;
t22 = t113 * t52 + t115 * t51;
t21 = -t25 + t226;
t20 = (-t139 * t270 + t278) * t138;
t6 = t226 - t9;
t5 = [-t204 * t100 + t144 * t101 - t138 * t88 + Icges(2,3) + (t89 - t270) * t139 + m(7) * (t28 ^ 2 + t29 ^ 2) + m(6) * (t45 ^ 2 + t46 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t83 ^ 2 + t84 ^ 2) + m(3) * (t103 ^ 2 + t104 ^ 2) + m(2) * (t158 ^ 2 + t159 ^ 2) + t278 + ((Icges(3,2) + Icges(4,3)) * t187 + t306) * t187 + (t307 * t183 - t260 + t264) * t183; (t230 * t188 + (Icges(3,6) * t297 + Icges(4,6) * t298 + t209 * t292 + t212 * t299) * t187 + (t304 * t297 + t300 * t299) * t183 - t193) * t188 + (t230 * t184 + (Icges(3,6) * t292 + Icges(4,6) * t299 + t209 * t298 + t212 * t297) * t187 + (t304 * t292 + t300 * t297) * t183 - t192) * t184 + m(7) * (t28 * t36 + t29 * t35) + m(6) * (t45 * t54 + t46 * t53) + m(5) * (t58 * t72 + t59 * t71) + m(4) * (t105 * t84 + t106 * t83) + m(3) * (-t103 * t188 - t104 * t184) * t157; m(7) * (t35 ^ 2 + t36 ^ 2 + t6 ^ 2) + m(6) * (t21 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t34 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t105 ^ 2 + t106 ^ 2 + t56 ^ 2) + m(3) * (t157 ^ 2 * t240 + t75 ^ 2) + (t290 * t179 + t238) * t188 + (t289 * t178 + (t290 * t184 + t289 * t188) * t188 + t237) * t184; 0.2e1 * ((t184 * t29 + t188 * t28) * t285 + (t184 * t46 + t188 * t45) * t286 + t202 / 0.2e1 + (t184 * t84 + t188 * t83) * t287) * t183; m(7) * (-t187 * t6 + (t184 * t35 + t188 * t36) * t183) + m(6) * (-t187 * t21 + (t184 * t53 + t188 * t54) * t183) + m(5) * (t183 * t218 - t187 * t34) + m(4) * (-t187 * t56 + (t105 * t184 + t106 * t188) * t183); 0.2e1 * (t287 + m(5) / 0.2e1 + t236) * (t183 ^ 2 * t240 + t187 ^ 2); t193 * t188 + t192 * t184 + m(7) * (t28 * t40 + t29 * t39) + m(6) * (t45 * t61 + t46 * t60) + t102 * t202; m(7) * (t35 * t39 + t36 * t40 + t6 * t9) + m(6) * (t21 * t25 + t53 * t60 + t54 * t61) + m(5) * (t102 * t218 + t47 * t34) + t288; m(5) * (t102 * t183 * t240 - t187 * t47) + m(6) * (-t25 * t187 + (t184 * t60 + t188 * t61) * t183) + m(7) * (-t9 * t187 + (t184 * t39 + t188 * t40) * t183); m(7) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(6) * (t25 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t102 ^ 2 * t240 + t47 ^ 2) - t288; m(7) * (-t184 * t28 + t188 * t29) + m(6) * (-t184 * t45 + t188 * t46); m(7) * (-t184 * t36 + t188 * t35) + m(6) * (-t184 * t54 + t188 * t53); 0; m(7) * (-t184 * t40 + t188 * t39) + m(6) * (-t184 * t61 + t188 * t60); 0.2e1 * t236 * t240; t20 + m(7) * (t26 * t28 + t27 * t29) + t233 * t115 - t234 * t113; m(7) * (t22 * t6 + t26 * t36 + t27 * t35) + t190; m(7) * (-t22 * t187 + (t184 * t27 + t188 * t26) * t183); m(7) * (t22 * t9 + t26 * t40 + t27 * t39) - t190; m(7) * (-t26 * t184 + t188 * t27); t115 * t2 - t113 * t1 + t138 * (-t16 * t113 + t17 * t115 + t20) + m(7) * (t22 ^ 2 + t26 ^ 2 + t27 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
