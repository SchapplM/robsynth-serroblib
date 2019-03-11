% Calculate joint inertia matrix for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:24
% EndTime: 2019-03-09 15:21:33
% DurationCPUTime: 4.08s
% Computational Cost: add. (8937->426), mult. (7135->615), div. (0->0), fcn. (7378->12), ass. (0->209)
t302 = Icges(4,3) + Icges(5,3);
t198 = qJ(2) + qJ(3);
t186 = pkin(10) + t198;
t180 = sin(t186);
t181 = cos(t186);
t187 = sin(t198);
t188 = cos(t198);
t301 = Icges(4,5) * t188 + Icges(5,5) * t181 - Icges(4,6) * t187 - Icges(5,6) * t180;
t203 = sin(qJ(1));
t205 = cos(qJ(1));
t300 = t203 * t302 + t301 * t205;
t299 = -t301 * t203 + t205 * t302;
t267 = Icges(5,4) * t181;
t227 = -Icges(5,2) * t180 + t267;
t115 = Icges(5,6) * t203 + t227 * t205;
t268 = Icges(5,4) * t180;
t230 = Icges(5,1) * t181 - t268;
t117 = Icges(5,5) * t203 + t230 * t205;
t269 = Icges(4,4) * t188;
t228 = -Icges(4,2) * t187 + t269;
t130 = Icges(4,6) * t203 + t228 * t205;
t270 = Icges(4,4) * t187;
t231 = Icges(4,1) * t188 - t270;
t132 = Icges(4,5) * t203 + t231 * t205;
t298 = t115 * t180 - t117 * t181 + t130 * t187 - t132 * t188;
t114 = -Icges(5,6) * t205 + t227 * t203;
t116 = -Icges(5,5) * t205 + t230 * t203;
t129 = -Icges(4,6) * t205 + t228 * t203;
t131 = -Icges(4,5) * t205 + t231 * t203;
t297 = t114 * t180 - t116 * t181 + t129 * t187 - t131 * t188;
t196 = t203 ^ 2;
t296 = t203 * pkin(7);
t200 = cos(pkin(11));
t182 = pkin(5) * t200 + pkin(4);
t201 = -pkin(9) - qJ(5);
t199 = sin(pkin(11));
t259 = t199 * t203;
t264 = t181 * t205;
t265 = t180 * t205;
t195 = pkin(11) + qJ(6);
t185 = cos(t195);
t261 = t185 * t203;
t184 = sin(t195);
t262 = t184 * t205;
t135 = -t181 * t262 + t261;
t260 = t185 * t205;
t263 = t184 * t203;
t136 = t181 * t260 + t263;
t73 = t136 * rSges(7,1) + t135 * rSges(7,2) + rSges(7,3) * t265;
t295 = pkin(5) * t259 + t182 * t264 - t201 * t265 + t73;
t294 = Icges(4,5) * t187 + Icges(5,5) * t180 + Icges(4,6) * t188 + Icges(5,6) * t181;
t150 = Icges(5,2) * t181 + t268;
t151 = Icges(5,1) * t180 + t267;
t158 = Icges(4,2) * t188 + t270;
t159 = Icges(4,1) * t187 + t269;
t293 = -t150 * t180 + t151 * t181 - t158 * t187 + t159 * t188;
t237 = rSges(4,1) * t188 - rSges(4,2) * t187;
t255 = qJ(5) + t201;
t282 = -pkin(4) + t182;
t99 = -t181 * rSges(7,3) + (rSges(7,1) * t185 - rSges(7,2) * t184) * t180;
t292 = -t282 * t180 - t255 * t181 - t99;
t256 = t200 * t205;
t145 = -t181 * t259 - t256;
t257 = t200 * t203;
t258 = t199 * t205;
t146 = t181 * t257 - t258;
t197 = t205 ^ 2;
t266 = t180 * t203;
t78 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t266;
t147 = -t181 * t258 + t257;
t148 = t181 * t256 + t259;
t79 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t265;
t133 = -t181 * t263 - t260;
t134 = t181 * t261 - t262;
t66 = Icges(7,5) * t134 + Icges(7,6) * t133 + Icges(7,3) * t266;
t68 = Icges(7,4) * t134 + Icges(7,2) * t133 + Icges(7,6) * t266;
t70 = Icges(7,1) * t134 + Icges(7,4) * t133 + Icges(7,5) * t266;
t17 = t133 * t68 + t134 * t70 + t66 * t266;
t67 = Icges(7,5) * t136 + Icges(7,6) * t135 + Icges(7,3) * t265;
t69 = Icges(7,4) * t136 + Icges(7,2) * t135 + Icges(7,6) * t265;
t71 = Icges(7,1) * t136 + Icges(7,4) * t135 + Icges(7,5) * t265;
t18 = t133 * t69 + t134 * t71 + t67 * t266;
t8 = -t17 * t205 + t18 * t203;
t80 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t266;
t81 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t265;
t82 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t266;
t83 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t265;
t291 = -t8 + (t145 * t80 + t146 * t82 + t78 * t266) * t205 + t299 * t197 + (-t145 * t81 - t146 * t83 - t79 * t266 + t298 * t203 + (-t297 + t300) * t205) * t203;
t290 = t181 ^ 2;
t289 = m(6) / 0.2e1;
t288 = m(7) / 0.2e1;
t206 = -pkin(8) - pkin(7);
t287 = t203 / 0.2e1;
t286 = -t205 / 0.2e1;
t202 = sin(qJ(2));
t285 = pkin(2) * t202;
t284 = pkin(3) * t187;
t283 = pkin(4) * t181;
t204 = cos(qJ(2));
t183 = t204 * pkin(2) + pkin(1);
t161 = pkin(3) * t188 + t183;
t155 = t205 * t161;
t174 = t205 * t183;
t281 = t205 * (t155 - t174) + (t161 - t183) * t196;
t280 = rSges(3,1) * t204;
t278 = rSges(3,2) * t202;
t97 = -Icges(7,6) * t181 + (Icges(7,4) * t185 - Icges(7,2) * t184) * t180;
t276 = t184 * t97;
t275 = t205 * rSges(3,3);
t23 = -t181 * t66 + (-t184 * t68 + t185 * t70) * t180;
t274 = t23 * t205;
t24 = -t181 * t67 + (-t184 * t69 + t185 * t71) * t180;
t273 = t24 * t203;
t272 = Icges(3,4) * t202;
t271 = Icges(3,4) * t204;
t193 = t205 * pkin(7);
t254 = t203 * (t193 + (-pkin(1) + t183) * t203) + t205 * (-t205 * pkin(1) + t174 - t296);
t215 = t203 * rSges(4,3) + t205 * t237;
t84 = t203 * (-t205 * rSges(4,3) + t237 * t203) + t205 * t215;
t253 = pkin(4) * t264 + qJ(5) * t265;
t252 = t203 * rSges(3,3) + t205 * t280;
t250 = t196 + t197;
t249 = t289 + t288;
t248 = t148 * rSges(6,1) + t147 * rSges(6,2) + rSges(6,3) * t265;
t19 = t135 * t68 + t136 * t70 + t66 * t265;
t20 = t135 * t69 + t136 * t71 + t67 * t265;
t9 = -t19 * t205 + t20 * t203;
t247 = (t9 + (t147 * t81 + t148 * t83 + t79 * t265) * t203 + t300 * t196 + (-t147 * t80 - t148 * t82 - t78 * t265 + (-t298 + t299) * t203 + t297 * t205) * t205) * t203;
t160 = rSges(4,1) * t187 + rSges(4,2) * t188;
t246 = -t160 - t285;
t152 = t180 * pkin(4) - t181 * qJ(5);
t245 = -t152 - t284;
t153 = rSges(5,1) * t180 + rSges(5,2) * t181;
t244 = -t153 - t284;
t194 = -qJ(4) + t206;
t243 = -t203 * t194 + t155;
t242 = t196 * (qJ(5) * t180 + t283) + t205 * t253 + t281;
t214 = rSges(5,1) * t264 - rSges(5,2) * t265 + t203 * rSges(5,3);
t236 = rSges(5,1) * t181 - rSges(5,2) * t180;
t39 = t203 * (-t205 * rSges(5,3) + t236 * t203) + t205 * t214 + t281;
t96 = -Icges(7,3) * t181 + (Icges(7,5) * t185 - Icges(7,6) * t184) * t180;
t98 = -Icges(7,5) * t181 + (Icges(7,1) * t185 - Icges(7,4) * t184) * t180;
t34 = t133 * t97 + t134 * t98 + t96 * t266;
t3 = -t34 * t181 + (t17 * t203 + t18 * t205) * t180;
t35 = t135 * t97 + t136 * t98 + t96 * t265;
t4 = -t35 * t181 + (t19 * t203 + t20 * t205) * t180;
t241 = t3 * t286 + t4 * t287 - t181 * (t273 - t274) / 0.2e1 + t8 * t266 / 0.2e1 + t9 * t265 / 0.2e1;
t103 = -t181 * rSges(6,3) + (rSges(6,1) * t200 - rSges(6,2) * t199) * t180;
t240 = -t103 + t245;
t239 = -t284 - t285;
t238 = -t278 + t280;
t235 = -t146 * rSges(6,1) - t145 * rSges(6,2);
t234 = -t134 * rSges(7,1) - t133 * rSges(7,2);
t233 = t245 + t292;
t232 = Icges(3,1) * t204 - t272;
t229 = -Icges(3,2) * t202 + t271;
t226 = Icges(3,5) * t204 - Icges(3,6) * t202;
t16 = t203 * (rSges(6,3) * t266 - t235) + t205 * t248 + t242;
t213 = -t152 + t239;
t212 = -t153 + t239;
t210 = -t103 + t213;
t72 = rSges(7,3) * t266 - t234;
t14 = t242 + (-t253 + t295) * t205 + (-pkin(5) * t258 + t72 + (-t255 * t180 + t282 * t181) * t203) * t203;
t209 = t213 + t292;
t208 = t205 * t291 + t247;
t100 = -Icges(6,3) * t181 + (Icges(6,5) * t200 - Icges(6,6) * t199) * t180;
t101 = -Icges(6,6) * t181 + (Icges(6,4) * t200 - Icges(6,2) * t199) * t180;
t102 = -Icges(6,5) * t181 + (Icges(6,1) * t200 - Icges(6,4) * t199) * t180;
t207 = -t274 / 0.2e1 + t273 / 0.2e1 + (t100 * t265 + t101 * t147 + t102 * t148 + t130 * t188 + t132 * t187 + t35 + t293 * t205 + t294 * t203 + (-t79 + t115) * t181 + (-t199 * t81 + t200 * t83 + t117) * t180) * t287 + (t100 * t266 + t101 * t145 + t102 * t146 + t129 * t188 + t131 * t187 + t34 - t294 * t205 + t293 * t203 + (-t78 + t114) * t181 + (-t199 * t80 + t200 * t82 + t116) * t180) * t286;
t173 = rSges(2,1) * t205 - rSges(2,2) * t203;
t172 = -rSges(2,1) * t203 - rSges(2,2) * t205;
t171 = rSges(3,1) * t202 + rSges(3,2) * t204;
t140 = Icges(3,3) * t203 + t226 * t205;
t139 = -Icges(3,3) * t205 + t226 * t203;
t126 = t246 * t205;
t125 = t246 * t203;
t109 = t296 + (pkin(1) - t278) * t205 + t252;
t108 = t275 + t193 + (-pkin(1) - t238) * t203;
t107 = t244 * t205;
t106 = t244 * t203;
t95 = t212 * t205;
t94 = t212 * t203;
t93 = -t203 * t206 + t174 + t215;
t92 = (rSges(4,3) - t206) * t205 + (-t183 - t237) * t203;
t88 = t205 * (-t205 * t278 + t252) + (t238 * t203 - t275) * t203;
t87 = t180 * t185 * t98;
t86 = t214 + t243;
t85 = (rSges(5,3) - t194) * t205 + (-t161 - t236) * t203;
t63 = t240 * t205;
t62 = t240 * t203;
t57 = t210 * t205;
t56 = t210 * t203;
t51 = t243 + t248 + t253;
t50 = -t205 * t194 + (-t283 - t161 + (-rSges(6,3) - qJ(5)) * t180) * t203 + t235;
t49 = t84 + t254;
t48 = t233 * t205;
t47 = t233 * t203;
t46 = -t181 * t73 - t99 * t265;
t45 = t181 * t72 + t99 * t266;
t44 = t209 * t205;
t43 = t209 * t203;
t42 = t243 + t295;
t41 = (pkin(5) * t199 - t194) * t205 + (-t181 * t182 - t161 + (-rSges(7,3) + t201) * t180) * t203 + t234;
t40 = -t180 * t276 - t181 * t96 + t87;
t38 = (-t203 * t73 + t205 * t72) * t180;
t27 = t39 + t254;
t15 = t16 + t254;
t11 = t14 + t254;
t1 = [t188 * t158 + t187 * t159 + t204 * (Icges(3,2) * t204 + t272) + t202 * (Icges(3,1) * t202 + t271) + Icges(2,3) + t87 + (-t96 + t150 - t100) * t181 + (-t101 * t199 + t102 * t200 + t151 - t276) * t180 + m(7) * (t41 ^ 2 + t42 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2) + m(5) * (t85 ^ 2 + t86 ^ 2) + m(4) * (t92 ^ 2 + t93 ^ 2) + m(3) * (t108 ^ 2 + t109 ^ 2) + m(2) * (t172 ^ 2 + t173 ^ 2); m(3) * (-t108 * t205 - t109 * t203) * t171 + m(7) * (t41 * t44 + t42 * t43) + m(6) * (t50 * t57 + t51 * t56) + m(5) * (t85 * t95 + t86 * t94) + m(4) * (t125 * t93 + t126 * t92) + (t196 / 0.2e1 + t197 / 0.2e1) * (Icges(3,5) * t202 + Icges(3,6) * t204) + t207 + (t204 * (Icges(3,6) * t203 + t229 * t205) + t202 * (Icges(3,5) * t203 + t232 * t205)) * t287 + (t204 * (-Icges(3,6) * t205 + t229 * t203) + t202 * (-Icges(3,5) * t205 + t232 * t203)) * t286; m(7) * (t11 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t15 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t27 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(4) * (t125 ^ 2 + t126 ^ 2 + t49 ^ 2) + m(3) * (t250 * t171 ^ 2 + t88 ^ 2) + t203 * t196 * t140 + t247 + (-t197 * t139 + (-t203 * t139 + t205 * t140) * t203 + t291) * t205; m(7) * (t41 * t48 + t42 * t47) + m(6) * (t50 * t63 + t51 * t62) + m(5) * (t106 * t86 + t107 * t85) + m(4) * (-t203 * t93 - t205 * t92) * t160 + t207; m(7) * (t11 * t14 + t43 * t47 + t44 * t48) + m(6) * (t15 * t16 + t56 * t62 + t57 * t63) + m(5) * (t106 * t94 + t107 * t95 + t27 * t39) + m(4) * (t84 * t49 + (-t125 * t203 - t126 * t205) * t160) + t208; m(7) * (t14 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(5) * (t106 ^ 2 + t107 ^ 2 + t39 ^ 2) + m(6) * (t16 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(4) * (t250 * t160 ^ 2 + t84 ^ 2) + t208; m(7) * (t203 * t41 - t205 * t42) + m(6) * (t203 * t50 - t205 * t51) + m(5) * (t203 * t85 - t205 * t86); m(7) * (t203 * t44 - t205 * t43) + m(6) * (t203 * t57 - t205 * t56) + m(5) * (t203 * t95 - t205 * t94); m(7) * (t203 * t48 - t205 * t47) + m(5) * (-t106 * t205 + t107 * t203) + m(6) * (t203 * t63 - t205 * t62); 0.2e1 * (m(5) / 0.2e1 + t249) * t250; 0.2e1 * ((t203 * t42 + t205 * t41) * t288 + (t203 * t51 + t205 * t50) * t289) * t180; m(7) * (-t181 * t11 + (t203 * t43 + t205 * t44) * t180) + m(6) * (-t181 * t15 + (t203 * t56 + t205 * t57) * t180); m(7) * (-t181 * t14 + (t203 * t47 + t205 * t48) * t180) + m(6) * (-t181 * t16 + (t203 * t62 + t205 * t63) * t180); 0; 0.2e1 * t249 * (t250 * t180 ^ 2 + t290); m(7) * (t41 * t45 + t42 * t46) - t40 * t181 + ((t35 / 0.2e1 + t24 / 0.2e1) * t205 + (t34 / 0.2e1 + t23 / 0.2e1) * t203) * t180; m(7) * (t11 * t38 + t43 * t46 + t44 * t45) + t241; m(7) * (t14 * t38 + t45 * t48 + t46 * t47) + t241; m(7) * (t203 * t45 - t205 * t46); m(7) * (-t38 * t181 + (t203 * t46 + t205 * t45) * t180); t290 * t40 + m(7) * (t38 ^ 2 + t45 ^ 2 + t46 ^ 2) + (t205 * t4 + t203 * t3 - t181 * (t203 * t23 + t205 * t24)) * t180;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
