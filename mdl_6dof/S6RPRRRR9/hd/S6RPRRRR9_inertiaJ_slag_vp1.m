% Calculate joint inertia matrix for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:04
% EndTime: 2019-03-09 07:23:12
% DurationCPUTime: 4.14s
% Computational Cost: add. (11494->477), mult. (13324->669), div. (0->0), fcn. (14418->10), ass. (0->238)
t232 = cos(qJ(3));
t315 = Icges(4,5) * t232;
t229 = sin(qJ(3));
t314 = Icges(4,6) * t229;
t313 = t315 / 0.2e1 - t314 / 0.2e1;
t233 = cos(qJ(1));
t312 = (rSges(4,1) * t229 + rSges(4,2) * t232) * t233;
t230 = sin(qJ(1));
t225 = t230 ^ 2;
t226 = t233 ^ 2;
t311 = -pkin(1) - pkin(7);
t234 = -pkin(9) - pkin(8);
t227 = qJ(4) + qJ(5);
t219 = qJ(6) + t227;
t213 = cos(t219);
t281 = t233 * t213;
t212 = sin(t219);
t291 = t230 * t212;
t158 = -t229 * t291 + t281;
t282 = t233 * t212;
t290 = t230 * t213;
t159 = t229 * t290 + t282;
t285 = t230 * t232;
t101 = Icges(7,1) * t159 + Icges(7,4) * t158 - Icges(7,5) * t285;
t97 = Icges(7,5) * t159 + Icges(7,6) * t158 - Icges(7,3) * t285;
t99 = Icges(7,4) * t159 + Icges(7,2) * t158 - Icges(7,6) * t285;
t33 = t159 * t101 + t158 * t99 - t285 * t97;
t160 = t229 * t282 + t290;
t161 = -t229 * t281 + t291;
t284 = t232 * t233;
t100 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t284;
t102 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t284;
t98 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t284;
t34 = t158 * t100 + t159 * t102 - t285 * t98;
t145 = Icges(7,3) * t229 + (Icges(7,5) * t213 - Icges(7,6) * t212) * t232;
t146 = Icges(7,6) * t229 + (Icges(7,4) * t213 - Icges(7,2) * t212) * t232;
t147 = Icges(7,5) * t229 + (Icges(7,1) * t213 - Icges(7,4) * t212) * t232;
t62 = -t145 * t285 + t158 * t146 + t159 * t147;
t4 = t62 * t229 + (-t230 * t33 + t233 * t34) * t232;
t216 = cos(t227);
t279 = t233 * t216;
t215 = sin(t227);
t289 = t230 * t215;
t174 = -t229 * t289 + t279;
t280 = t233 * t215;
t288 = t230 * t216;
t175 = t229 * t288 + t280;
t106 = Icges(6,5) * t175 + Icges(6,6) * t174 - Icges(6,3) * t285;
t108 = Icges(6,4) * t175 + Icges(6,2) * t174 - Icges(6,6) * t285;
t110 = Icges(6,1) * t175 + Icges(6,4) * t174 - Icges(6,5) * t285;
t40 = -t106 * t285 + t174 * t108 + t175 * t110;
t176 = t229 * t280 + t288;
t177 = -t229 * t279 + t289;
t107 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t284;
t109 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t284;
t111 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t284;
t41 = -t107 * t285 + t174 * t109 + t175 * t111;
t151 = Icges(6,3) * t229 + (Icges(6,5) * t216 - Icges(6,6) * t215) * t232;
t152 = Icges(6,6) * t229 + (Icges(6,4) * t216 - Icges(6,2) * t215) * t232;
t153 = Icges(6,5) * t229 + (Icges(6,1) * t216 - Icges(6,4) * t215) * t232;
t68 = -t151 * t285 + t174 * t152 + t175 * t153;
t9 = t68 * t229 + (-t230 * t40 + t233 * t41) * t232;
t310 = -t4 - t9;
t308 = t229 / 0.2e1;
t307 = t230 / 0.2e1;
t305 = t233 / 0.2e1;
t46 = t229 * t97 + (t101 * t213 - t212 * t99) * t232;
t47 = t229 * t98 + (-t100 * t212 + t102 * t213) * t232;
t35 = t161 * t101 + t160 * t99 + t284 * t97;
t36 = t160 * t100 + t161 * t102 + t284 * t98;
t63 = t145 * t284 + t160 * t146 + t161 * t147;
t5 = t63 * t229 + (-t230 * t35 + t233 * t36) * t232;
t270 = t232 * t213 * t147 + t229 * t145;
t296 = t212 * t146;
t75 = (-t232 * t296 + t270) * t229;
t304 = t5 * t284 + t229 * (t75 + (-t230 * t46 + t233 * t47) * t232);
t200 = t232 * rSges(4,1) - t229 * rSges(4,2);
t303 = m(4) * t200;
t228 = sin(qJ(4));
t302 = t228 * pkin(4);
t231 = cos(qJ(4));
t214 = t231 * pkin(4) + pkin(3);
t103 = t159 * rSges(7,1) + t158 * rSges(7,2) - rSges(7,3) * t285;
t278 = t233 * t228;
t283 = t232 * t234;
t293 = t229 * t230;
t256 = pkin(4) * t278 + t214 * t293 + t230 * t283;
t191 = pkin(5) * t216 + t214;
t194 = pkin(5) * t215 + t302;
t224 = -pkin(10) + t234;
t257 = t191 * t293 + t233 * t194 + t224 * t285;
t90 = -t256 + t257;
t301 = -t103 - t90;
t244 = -t161 * rSges(7,1) - t160 * rSges(7,2);
t104 = rSges(7,3) * t284 - t244;
t292 = t229 * t233;
t264 = t214 * t292 + t233 * t283;
t297 = t191 * t229;
t300 = t104 + (-t224 * t232 - t297) * t233 + (t194 - t302) * t230 + t264;
t149 = t229 * rSges(7,3) + (rSges(7,1) * t213 - rSges(7,2) * t212) * t232;
t80 = t229 * t103 + t149 * t285;
t295 = t215 * t152;
t167 = Icges(5,6) * t229 + (Icges(5,4) * t231 - Icges(5,2) * t228) * t232;
t294 = t228 * t167;
t287 = t230 * t228;
t286 = t230 * t231;
t277 = t233 * t231;
t112 = t175 * rSges(6,1) + t174 * rSges(6,2) - rSges(6,3) * t285;
t209 = pkin(3) * t293;
t188 = -pkin(8) * t285 + t209;
t127 = -t188 + t256;
t276 = -t112 - t127;
t113 = t177 * rSges(6,1) + t176 * rSges(6,2) + rSges(6,3) * t284;
t211 = pkin(3) * t292;
t247 = pkin(8) * t284 - t211;
t128 = pkin(4) * t287 - t247 - t264;
t275 = -t113 - t128;
t150 = (-pkin(3) + t214) * t232 + (-pkin(8) - t234) * t229;
t274 = t229 * t127 + t150 * t285;
t181 = t233 * t247;
t273 = t233 * t128 + t181;
t129 = (t191 - t214) * t232 + (-t224 + t234) * t229;
t137 = t149 * t284;
t272 = t129 * t284 + t137;
t271 = t129 + t149;
t269 = t232 * t216 * t153 + t229 * t151;
t156 = t229 * rSges(6,3) + (rSges(6,1) * t216 - rSges(6,2) * t215) * t232;
t84 = t229 * t112 + t156 * t285;
t203 = t232 * pkin(3) + t229 * pkin(8);
t190 = t230 * t203;
t268 = t230 * t150 + t190;
t164 = Icges(5,3) * t229 + (Icges(5,5) * t231 - Icges(5,6) * t228) * t232;
t170 = Icges(5,5) * t229 + (Icges(5,1) * t231 - Icges(5,4) * t228) * t232;
t267 = t232 * t231 * t170 + t229 * t164;
t266 = -t150 - t203;
t184 = -t229 * t287 + t277;
t185 = t229 * t286 + t278;
t265 = t185 * rSges(5,1) + t184 * rSges(5,2);
t263 = t233 * pkin(1) + t230 * qJ(2);
t262 = t225 + t226;
t119 = Icges(5,5) * t185 + Icges(5,6) * t184 - Icges(5,3) * t285;
t121 = Icges(5,4) * t185 + Icges(5,2) * t184 - Icges(5,6) * t285;
t123 = Icges(5,1) * t185 + Icges(5,4) * t184 - Icges(5,5) * t285;
t58 = t229 * t119 + (-t121 * t228 + t123 * t231) * t232;
t76 = -t164 * t285 + t184 * t167 + t185 * t170;
t261 = -t58 / 0.2e1 - t76 / 0.2e1;
t186 = t229 * t278 + t286;
t187 = -t229 * t277 + t287;
t120 = Icges(5,5) * t187 + Icges(5,6) * t186 + Icges(5,3) * t284;
t122 = Icges(5,4) * t187 + Icges(5,2) * t186 + Icges(5,6) * t284;
t124 = Icges(5,1) * t187 + Icges(5,4) * t186 + Icges(5,5) * t284;
t59 = t229 * t120 + (-t122 * t228 + t124 * t231) * t232;
t77 = t164 * t284 + t186 * t167 + t187 * t170;
t260 = t59 / 0.2e1 + t77 / 0.2e1;
t259 = -t127 + t301;
t258 = -t128 - t300;
t255 = rSges(4,1) * t293 + rSges(4,2) * t285 + t233 * rSges(4,3);
t254 = t233 * pkin(7) + t263;
t253 = (-rSges(5,3) - pkin(8)) * t232;
t252 = -t285 / 0.2e1;
t251 = t284 / 0.2e1;
t42 = t106 * t284 + t176 * t108 + t177 * t110;
t43 = t107 * t284 + t176 * t109 + t177 * t111;
t69 = t151 * t284 + t176 * t152 + t177 * t153;
t10 = t69 * t229 + (-t230 * t42 + t233 * t43) * t232;
t50 = t229 * t106 + (-t108 * t215 + t110 * t216) * t232;
t51 = t229 * t107 + (-t109 * t215 + t111 * t216) * t232;
t79 = (-t232 * t295 + t269) * t229;
t250 = t10 * t284 + t229 * (t79 + (-t230 * t50 + t233 * t51) * t232) + t304;
t38 = t129 * t285 + t229 * t90 + t80;
t17 = t34 * t230 + t33 * t233;
t18 = t36 * t230 + t35 * t233;
t249 = t17 * t252 + t18 * t251 + t4 * t305 + t5 * t307 + (t47 * t230 + t46 * t233) * t308;
t248 = t75 + (t46 + t62) * t252 + (t47 + t63) * t251;
t245 = -t187 * rSges(5,1) - t186 * rSges(5,2);
t243 = -t285 * t4 + t304;
t240 = Icges(4,5) * t229 + Icges(4,6) * t232;
t237 = t285 * t310 + t250;
t21 = t41 * t230 + t40 * t233;
t22 = t43 * t230 + t42 * t233;
t236 = t10 * t307 + t21 * t252 + t22 * t251 + t9 * t305 + t249 + (t51 * t230 + t50 * t233) * t308;
t235 = t248 + t79 + (t50 + t68) * t252 + (t51 + t69) * t251;
t218 = t233 * qJ(2);
t201 = t233 * rSges(2,1) - t230 * rSges(2,2);
t199 = -t230 * rSges(2,1) - t233 * rSges(2,2);
t196 = -t314 + t315;
t179 = -t233 * rSges(3,2) + t230 * rSges(3,3) + t263;
t178 = t233 * rSges(3,3) + t218 + (rSges(3,2) - pkin(1)) * t230;
t173 = t229 * rSges(5,3) + (rSges(5,1) * t231 - rSges(5,2) * t228) * t232;
t166 = Icges(4,3) * t230 - t233 * t240;
t165 = Icges(4,3) * t233 + t230 * t240;
t142 = t156 * t284;
t139 = t150 * t284;
t134 = t254 + t255;
t133 = t218 + t312 + (-rSges(4,3) + t311) * t230;
t131 = (-t173 - t203) * t233;
t130 = t230 * t173 + t190;
t126 = rSges(5,3) * t284 - t245;
t125 = -rSges(5,3) * t285 + t265;
t114 = -t230 * t255 + (t230 * rSges(4,3) - t312) * t233;
t95 = (-t156 + t266) * t233;
t94 = t230 * t156 + t268;
t93 = t230 * t253 + t209 + t254 + t265;
t92 = t230 * t311 + t233 * t253 + t211 + t218 + t245;
t89 = -t229 * t126 + t173 * t284;
t88 = t229 * t125 + t173 * t285;
t86 = (-t232 * t294 + t267) * t229;
t85 = -t229 * t113 + t142;
t83 = t112 + t254 + t256;
t82 = t218 + (-t302 + t311) * t230 - t113 + t264;
t81 = -t229 * t104 + t137;
t78 = (-t125 * t233 - t126 * t230) * t232;
t74 = t103 + t254 + t257;
t73 = t218 + (t297 + (-rSges(7,3) + t224) * t232) * t233 + (-t194 + t311) * t230 + t244;
t72 = (t266 - t271) * t233;
t71 = t230 * t271 + t268;
t70 = (-t112 * t233 - t113 * t230) * t232;
t67 = (-t103 * t233 - t104 * t230) * t232;
t64 = t233 * t126 + t181 + (-t125 - t188) * t230;
t57 = t229 * t275 + t139 + t142;
t56 = t84 + t274;
t55 = t120 * t284 + t186 * t122 + t187 * t124;
t54 = t119 * t284 + t186 * t121 + t187 * t123;
t53 = -t120 * t285 + t184 * t122 + t185 * t124;
t52 = -t119 * t285 + t184 * t121 + t185 * t123;
t39 = -t229 * t300 + t272;
t37 = (t230 * t275 + t233 * t276) * t232;
t32 = t233 * t113 + (-t188 + t276) * t230 + t273;
t31 = (-t230 * t300 + t233 * t301) * t232;
t30 = t229 * t258 + t139 + t272;
t29 = t38 + t274;
t28 = t55 * t230 + t54 * t233;
t27 = t53 * t230 + t52 * t233;
t25 = (t230 * t258 + t233 * t259) * t232;
t23 = t300 * t233 + (-t188 + t259) * t230 + t273;
t16 = t77 * t229 + (-t230 * t54 + t233 * t55) * t232;
t15 = t76 * t229 + (-t230 * t52 + t233 * t53) * t232;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t232 - t294 - t295 - t296) * t232 + m(7) * (t73 ^ 2 + t74 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t92 ^ 2 + t93 ^ 2) + m(4) * (t133 ^ 2 + t134 ^ 2) + m(3) * (t178 ^ 2 + t179 ^ 2) + m(2) * (t199 ^ 2 + t201 ^ 2) + t267 + t269 + t270 + (-0.2e1 * Icges(4,4) * t232 + Icges(4,2) * t229) * t229; m(7) * (t230 * t73 - t233 * t74) + m(6) * (t230 * t82 - t233 * t83) + m(5) * (t230 * t92 - t233 * t93) + m(4) * (t230 * t133 - t233 * t134) + m(3) * (t230 * t178 - t233 * t179); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t262; m(7) * (t71 * t73 + t72 * t74) + m(6) * (t82 * t94 + t83 * t95) + m(5) * (t130 * t92 + t131 * t93) + (t62 / 0.2e1 + t46 / 0.2e1 + t68 / 0.2e1 + t50 / 0.2e1 - t134 * t303 + t196 * t305 - t261 + t313 * t233) * t233 + (t47 / 0.2e1 + t63 / 0.2e1 + t69 / 0.2e1 + t51 / 0.2e1 + t133 * t303 + t196 * t307 + t260 + t313 * t230) * t230; m(5) * (t130 * t230 - t131 * t233) + m(6) * (t94 * t230 - t95 * t233) + m(7) * (t71 * t230 - t72 * t233) + t262 * t303; m(7) * (t23 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t32 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(5) * (t130 ^ 2 + t131 ^ 2 + t64 ^ 2) + m(4) * (t200 ^ 2 * t262 + t114 ^ 2) + (t226 * t165 + t17 + t21 + t27) * t233 + (t225 * t166 + t18 + t22 + t28 + (t230 * t165 + t233 * t166) * t233) * t230; t86 + (t230 * t261 + t233 * t260) * t232 + m(7) * (t29 * t74 + t30 * t73) + m(6) * (t56 * t83 + t57 * t82) + m(5) * (t88 * t93 + t89 * t92) + t235; m(5) * (t89 * t230 - t88 * t233) + m(6) * (t57 * t230 - t56 * t233) + m(7) * (t30 * t230 - t29 * t233); t236 + m(7) * (t23 * t25 + t29 * t72 + t30 * t71) + m(6) * (t32 * t37 + t56 * t95 + t57 * t94) + m(5) * (t130 * t89 + t131 * t88 + t64 * t78) + (t28 * t305 - t230 * t27 / 0.2e1) * t232 + t16 * t307 + (t59 * t230 + t58 * t233) * t308 + t15 * t305; t229 * t86 + m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t37 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t78 ^ 2 + t88 ^ 2 + t89 ^ 2) + ((t229 * t59 + t16) * t233 + (-t229 * t58 - t15 + t310) * t230) * t232 + t250; m(7) * (t38 * t74 + t39 * t73) + m(6) * (t82 * t85 + t83 * t84) + t235; m(6) * (t85 * t230 - t84 * t233) + m(7) * (t39 * t230 - t38 * t233); m(7) * (t23 * t31 + t38 * t72 + t39 * t71) + m(6) * (t32 * t70 + t84 * t95 + t85 * t94) + t236; m(7) * (t25 * t31 + t29 * t38 + t30 * t39) + m(6) * (t37 * t70 + t56 * t84 + t57 * t85) + t237; m(7) * (t31 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(6) * (t70 ^ 2 + t84 ^ 2 + t85 ^ 2) + t237; m(7) * (t73 * t81 + t74 * t80) + t248; m(7) * (t81 * t230 - t80 * t233); m(7) * (t23 * t67 + t71 * t81 + t72 * t80) + t249; m(7) * (t25 * t67 + t29 * t80 + t30 * t81) + t243; m(7) * (t31 * t67 + t38 * t80 + t39 * t81) + t243; m(7) * (t67 ^ 2 + t80 ^ 2 + t81 ^ 2) + t243;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
