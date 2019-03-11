% Calculate joint inertia matrix for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:21
% EndTime: 2019-03-09 07:12:31
% DurationCPUTime: 4.03s
% Computational Cost: add. (16351->482), mult. (13330->682), div. (0->0), fcn. (14430->12), ass. (0->243)
t222 = pkin(11) + qJ(3);
t215 = sin(t222);
t313 = Icges(4,5) * t215;
t312 = t313 / 0.2e1;
t226 = qJ(4) + qJ(5);
t219 = qJ(6) + t226;
t212 = sin(t219);
t213 = cos(t219);
t231 = sin(qJ(1));
t216 = cos(t222);
t233 = cos(qJ(1));
t286 = t216 * t233;
t171 = -t212 * t286 + t213 * t231;
t172 = t212 * t231 + t213 * t286;
t289 = t215 * t233;
t110 = t172 * rSges(7,1) + t171 * rSges(7,2) + rSges(7,3) * t289;
t232 = cos(qJ(4));
t214 = t232 * pkin(4) + pkin(3);
t218 = cos(t226);
t196 = pkin(5) * t218 + t214;
t217 = sin(t226);
t230 = sin(qJ(4));
t197 = pkin(4) * t230 + pkin(5) * t217;
t311 = t196 * t286 + t231 * t197 + t110;
t224 = t231 ^ 2;
t225 = t233 ^ 2;
t234 = -pkin(9) - pkin(8);
t290 = t215 * t231;
t287 = t216 * t231;
t169 = -t212 * t287 - t213 * t233;
t170 = -t212 * t233 + t213 * t287;
t103 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t290;
t105 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t290;
t107 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t290;
t38 = t103 * t290 + t105 * t169 + t107 * t170;
t104 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t289;
t106 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t289;
t108 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t289;
t39 = t104 * t290 + t106 * t169 + t108 * t170;
t147 = -Icges(7,3) * t216 + (Icges(7,5) * t213 - Icges(7,6) * t212) * t215;
t148 = -Icges(7,6) * t216 + (Icges(7,4) * t213 - Icges(7,2) * t212) * t215;
t149 = -Icges(7,5) * t216 + (Icges(7,1) * t213 - Icges(7,4) * t212) * t215;
t64 = t147 * t290 + t148 * t169 + t149 * t170;
t5 = -t64 * t216 + (t231 * t38 + t233 * t39) * t215;
t40 = t103 * t289 + t105 * t171 + t107 * t172;
t41 = t104 * t289 + t106 * t171 + t108 * t172;
t65 = t147 * t289 + t148 * t171 + t149 * t172;
t6 = -t65 * t216 + (t231 * t40 + t233 * t41) * t215;
t310 = t6 * t289 + t5 * t290;
t309 = -t216 / 0.2e1;
t308 = t231 / 0.2e1;
t307 = -t233 / 0.2e1;
t306 = t233 / 0.2e1;
t193 = rSges(4,1) * t215 + rSges(4,2) * t216;
t305 = m(4) * t193;
t304 = pkin(3) * t216;
t303 = pkin(8) * t215;
t302 = -pkin(3) + t214;
t138 = t215 * t213 * t149;
t293 = t148 * t212;
t75 = -t216 * t147 - t215 * t293 + t138;
t298 = t75 * t216;
t48 = -t216 * t103 + (-t105 * t212 + t107 * t213) * t215;
t49 = -t216 * t104 + (-t106 * t212 + t108 * t213) * t215;
t13 = -t298 + (t231 * t48 + t233 * t49) * t215;
t282 = t218 * t233;
t285 = t217 * t231;
t177 = -t216 * t285 - t282;
t283 = t218 * t231;
t284 = t217 * t233;
t178 = t216 * t283 - t284;
t112 = Icges(6,5) * t178 + Icges(6,6) * t177 + Icges(6,3) * t290;
t114 = Icges(6,4) * t178 + Icges(6,2) * t177 + Icges(6,6) * t290;
t116 = Icges(6,1) * t178 + Icges(6,4) * t177 + Icges(6,5) * t290;
t54 = -t216 * t112 + (-t114 * t217 + t116 * t218) * t215;
t179 = -t216 * t284 + t283;
t180 = t216 * t282 + t285;
t113 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t289;
t115 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t289;
t117 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t289;
t55 = -t216 * t113 + (-t115 * t217 + t117 * t218) * t215;
t153 = -Icges(6,5) * t216 + (Icges(6,1) * t218 - Icges(6,4) * t217) * t215;
t140 = t215 * t218 * t153;
t151 = -Icges(6,3) * t216 + (Icges(6,5) * t218 - Icges(6,6) * t217) * t215;
t152 = -Icges(6,6) * t216 + (Icges(6,4) * t218 - Icges(6,2) * t217) * t215;
t292 = t152 * t217;
t80 = -t216 * t151 - t215 * t292 + t140;
t301 = -t13 + t80 * t216 - (t231 * t54 + t233 * t55) * t215;
t300 = -t75 - t80;
t223 = -pkin(10) + t234;
t280 = t230 * t233;
t288 = t215 * t234;
t267 = pkin(4) * t280 + t231 * t288;
t268 = t196 - t214;
t97 = -t233 * t197 + (-t215 * t223 + t216 * t268) * t231 + t267;
t246 = -t170 * rSges(7,1) - t169 * rSges(7,2);
t109 = rSges(7,3) * t290 - t246;
t99 = t109 * t289;
t299 = t97 * t289 + t99;
t297 = rSges(3,3) + qJ(2);
t265 = t223 - t234;
t281 = t230 * t231;
t269 = -pkin(4) * t281 - t214 * t286;
t296 = -t265 * t289 + t269 + t311;
t294 = Icges(4,4) * t216;
t158 = -Icges(5,6) * t216 + (Icges(5,4) * t232 - Icges(5,2) * t230) * t215;
t291 = t158 * t230;
t279 = t231 * t232;
t278 = t232 * t233;
t229 = -pkin(7) - qJ(2);
t277 = t233 * t229;
t119 = t180 * rSges(6,1) + t179 * rSges(6,2) + rSges(6,3) * t289;
t238 = -t233 * t288 - t269;
t266 = pkin(3) * t286 + pkin(8) * t289;
t126 = t238 - t266;
t276 = -t119 - t126;
t125 = (t216 * t302 - t303) * t231 - t267;
t146 = (pkin(8) + t234) * t216 + t302 * t215;
t275 = t216 * t125 + t146 * t290;
t137 = t215 * t268 + t216 * t265;
t150 = -t216 * rSges(7,3) + (rSges(7,1) * t213 - rSges(7,2) * t212) * t215;
t274 = -t137 - t150;
t82 = t216 * t109 + t150 * t290;
t247 = -t178 * rSges(6,1) - t177 * rSges(6,2);
t118 = rSges(6,3) * t290 - t247;
t154 = -t216 * rSges(6,3) + (rSges(6,1) * t218 - rSges(6,2) * t217) * t215;
t87 = t216 * t118 + t154 * t290;
t273 = -t146 - t154;
t160 = -t216 * rSges(5,3) + (rSges(5,1) * t232 - rSges(5,2) * t230) * t215;
t195 = t215 * pkin(3) - t216 * pkin(8);
t272 = -t160 - t195;
t271 = t224 * (t303 + t304) + t233 * t266;
t264 = t224 + t225;
t186 = -t216 * t280 + t279;
t187 = t216 * t278 + t281;
t128 = Icges(5,5) * t187 + Icges(5,6) * t186 + Icges(5,3) * t289;
t130 = Icges(5,4) * t187 + Icges(5,2) * t186 + Icges(5,6) * t289;
t132 = Icges(5,1) * t187 + Icges(5,4) * t186 + Icges(5,5) * t289;
t61 = -t216 * t128 + (-t130 * t230 + t132 * t232) * t215;
t157 = -Icges(5,3) * t216 + (Icges(5,5) * t232 - Icges(5,6) * t230) * t215;
t159 = -Icges(5,5) * t216 + (Icges(5,1) * t232 - Icges(5,4) * t230) * t215;
t77 = t157 * t289 + t158 * t186 + t159 * t187;
t263 = t61 / 0.2e1 + t77 / 0.2e1;
t184 = -t216 * t281 - t278;
t185 = t216 * t279 - t280;
t127 = Icges(5,5) * t185 + Icges(5,6) * t184 + Icges(5,3) * t290;
t129 = Icges(5,4) * t185 + Icges(5,2) * t184 + Icges(5,6) * t290;
t131 = Icges(5,1) * t185 + Icges(5,4) * t184 + Icges(5,5) * t290;
t60 = -t216 * t127 + (-t129 * t230 + t131 * t232) * t215;
t76 = t157 * t290 + t158 * t184 + t159 * t185;
t262 = t76 / 0.2e1 + t60 / 0.2e1;
t261 = -t126 - t296;
t44 = t112 * t290 + t114 * t177 + t116 * t178;
t45 = t113 * t290 + t115 * t177 + t117 * t178;
t69 = t151 * t290 + t152 * t177 + t153 * t178;
t11 = -t69 * t216 + (t231 * t44 + t233 * t45) * t215;
t46 = t112 * t289 + t114 * t179 + t116 * t180;
t47 = t113 * t289 + t115 * t179 + t117 * t180;
t70 = t151 * t289 + t152 * t179 + t153 * t180;
t12 = -t70 * t216 + (t231 * t46 + t233 * t47) * t215;
t260 = t11 * t290 + t12 * t289 + t310;
t259 = -t146 + t274;
t258 = -t195 + t273;
t134 = t187 * rSges(5,1) + t186 * rSges(5,2) + rSges(5,3) * t289;
t257 = t290 / 0.2e1;
t256 = t289 / 0.2e1;
t255 = (t48 + t64) * t257 + (t49 + t65) * t256;
t228 = cos(pkin(11));
t211 = pkin(2) * t228 + pkin(1);
t254 = t233 * t211 - t231 * t229;
t36 = t137 * t290 + t216 * t97 + t82;
t253 = -t216 * t13 + t310;
t252 = t231 * t125 + t233 * t126 + t271;
t251 = -t195 + t259;
t19 = t231 * t39 - t233 * t38;
t20 = t231 * t41 - t233 * t40;
t250 = t19 * t257 + t20 * t256 + t5 * t307 + t6 * t308 + (t49 * t231 - t48 * t233) * t309;
t249 = rSges(4,1) * t216 - rSges(4,2) * t215;
t248 = -t185 * rSges(5,1) - t184 * rSges(5,2);
t244 = -Icges(4,2) * t215 + t294;
t243 = Icges(4,5) * t216 - Icges(4,6) * t215;
t240 = rSges(4,1) * t286 - rSges(4,2) * t289 + t231 * rSges(4,3);
t227 = sin(pkin(11));
t239 = rSges(3,1) * t228 - rSges(3,2) * t227 + pkin(1);
t237 = t216 * t301 + t260;
t236 = t255 + (t54 + t69) * t257 + (t55 + t70) * t256;
t23 = t231 * t45 - t233 * t44;
t24 = t231 * t47 - t233 * t46;
t235 = t11 * t307 + t12 * t308 + t23 * t257 + t24 * t256 + t250 + (t55 * t231 - t54 * t233) * t309;
t200 = rSges(2,1) * t233 - rSges(2,2) * t231;
t199 = -rSges(2,1) * t231 - rSges(2,2) * t233;
t189 = Icges(4,6) * t216 + t313;
t164 = Icges(4,3) * t231 + t233 * t243;
t163 = -Icges(4,3) * t233 + t231 * t243;
t156 = t231 * t297 + t233 * t239;
t155 = -t231 * t239 + t233 * t297;
t145 = t215 * t232 * t159;
t143 = t240 + t254;
t142 = (rSges(4,3) - t229) * t233 + (-t211 - t249) * t231;
t136 = t272 * t233;
t135 = t272 * t231;
t133 = rSges(5,3) * t290 - t248;
t122 = t233 * t240 + (-t233 * rSges(4,3) + t231 * t249) * t231;
t111 = t125 * t289;
t101 = t118 * t289;
t96 = t254 + t134 + t266;
t95 = -t277 + (-t304 - t211 + (-rSges(5,3) - pkin(8)) * t215) * t231 + t248;
t92 = t258 * t233;
t91 = t258 * t231;
t90 = -t134 * t216 - t160 * t289;
t89 = t133 * t216 + t160 * t290;
t88 = -t216 * t119 - t154 * t289;
t86 = -t216 * t157 - t215 * t291 + t145;
t85 = t238 + t254 + t119;
t84 = -t277 + (-rSges(6,3) * t215 - t214 * t216 - t211) * t231 + t247 + t267;
t83 = -t216 * t110 - t150 * t289;
t81 = (t133 * t233 - t134 * t231) * t215;
t79 = -t223 * t289 + t254 + t311;
t78 = (t197 - t229) * t233 + (-t196 * t216 - t211 + (-rSges(7,3) + t223) * t215) * t231 + t246;
t74 = -t119 * t290 + t101;
t73 = -t110 * t290 + t99;
t72 = t251 * t233;
t71 = t251 * t231;
t68 = t133 * t231 + t134 * t233 + t271;
t59 = t128 * t289 + t130 * t186 + t132 * t187;
t58 = t127 * t289 + t129 * t186 + t131 * t187;
t57 = t128 * t290 + t130 * t184 + t132 * t185;
t56 = t127 * t290 + t129 * t184 + t131 * t185;
t53 = t216 * t276 + t273 * t289;
t52 = t87 + t275;
t37 = -t216 * t296 + t274 * t289;
t35 = t276 * t290 + t101 + t111;
t34 = t118 * t231 + t119 * t233 + t252;
t33 = -t290 * t296 + t299;
t32 = t216 * t261 + t259 * t289;
t31 = t36 + t275;
t30 = t231 * t59 - t233 * t58;
t29 = t231 * t57 - t233 * t56;
t27 = t261 * t290 + t111 + t299;
t25 = t296 * t233 + (t109 + t97) * t231 + t252;
t16 = -t77 * t216 + (t231 * t58 + t233 * t59) * t215;
t15 = -t76 * t216 + (t231 * t56 + t233 * t57) * t215;
t1 = [Icges(3,2) * t228 ^ 2 + Icges(2,3) + t138 + t140 + t145 + (Icges(3,1) * t227 + 0.2e1 * Icges(3,4) * t228) * t227 + (Icges(4,4) * t215 + Icges(4,2) * t216 - t147 - t151 - t157) * t216 + (Icges(4,1) * t215 - t291 - t292 - t293 + t294) * t215 + m(7) * (t78 ^ 2 + t79 ^ 2) + m(6) * (t84 ^ 2 + t85 ^ 2) + m(5) * (t95 ^ 2 + t96 ^ 2) + m(4) * (t142 ^ 2 + t143 ^ 2) + m(3) * (t155 ^ 2 + t156 ^ 2) + m(2) * (t199 ^ 2 + t200 ^ 2); m(7) * (t231 * t78 - t233 * t79) + m(6) * (t231 * t84 - t233 * t85) + m(5) * (t231 * t95 - t233 * t96) + m(4) * (t142 * t231 - t143 * t233) + m(3) * (t155 * t231 - t156 * t233); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t264; m(7) * (t71 * t79 + t72 * t78) + m(6) * (t84 * t92 + t85 * t91) + m(5) * (t135 * t96 + t136 * t95) + (-t48 / 0.2e1 - t64 / 0.2e1 - t54 / 0.2e1 - t69 / 0.2e1 - t142 * t305 + (-Icges(4,6) * t233 + t231 * t244) * t309 + t233 * t312 + t189 * t306 - t262) * t233 + (t65 / 0.2e1 + t49 / 0.2e1 + t55 / 0.2e1 + t70 / 0.2e1 - t143 * t305 + t216 * (Icges(4,6) * t231 + t233 * t244) / 0.2e1 + t231 * t312 + t189 * t308 + t263) * t231; m(5) * (-t135 * t233 + t136 * t231) + m(6) * (t231 * t92 - t233 * t91) + m(7) * (t231 * t72 - t233 * t71); m(7) * (t25 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t34 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t135 ^ 2 + t136 ^ 2 + t68 ^ 2) + m(4) * (t193 ^ 2 * t264 + t122 ^ 2) + (-t225 * t163 - t19 - t23 - t29) * t233 + (t224 * t164 + t20 + t24 + t30 + (-t231 * t163 + t233 * t164) * t233) * t231; (-t86 + t300) * t216 + m(7) * (t31 * t78 + t32 * t79) + m(6) * (t52 * t84 + t53 * t85) + m(5) * (t89 * t95 + t90 * t96) + (t231 * t262 + t233 * t263) * t215 + t236; m(5) * (t231 * t89 - t233 * t90) + m(6) * (t231 * t52 - t233 * t53) + m(7) * (t231 * t31 - t233 * t32); m(7) * (t25 * t27 + t31 * t72 + t32 * t71) + m(6) * (t34 * t35 + t52 * t92 + t53 * t91) + m(5) * (t135 * t90 + t136 * t89 + t68 * t81) + (t29 * t308 + t30 * t306) * t215 + t16 * t308 + (t61 * t231 - t60 * t233) * t309 + t15 * t307 + t235; (t86 * t216 + t301) * t216 + (t233 * t16 + t231 * t15 - t216 * (t231 * t60 + t233 * t61)) * t215 + m(7) * (t27 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t35 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(5) * (t81 ^ 2 + t89 ^ 2 + t90 ^ 2) + t260; t300 * t216 + m(7) * (t36 * t78 + t37 * t79) + m(6) * (t84 * t87 + t85 * t88) + t236; m(6) * (t231 * t87 - t233 * t88) + m(7) * (t231 * t36 - t233 * t37); m(7) * (t25 * t33 + t36 * t72 + t37 * t71) + m(6) * (t34 * t74 + t87 * t92 + t88 * t91) + t235; m(7) * (t27 * t33 + t31 * t36 + t32 * t37) + m(6) * (t35 * t74 + t52 * t87 + t53 * t88) + t237; m(7) * (t33 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t74 ^ 2 + t87 ^ 2 + t88 ^ 2) + t237; m(7) * (t78 * t82 + t79 * t83) - t298 + t255; m(7) * (t231 * t82 - t233 * t83); m(7) * (t25 * t73 + t71 * t83 + t72 * t82) + t250; m(7) * (t27 * t73 + t31 * t82 + t32 * t83) + t253; m(7) * (t33 * t73 + t36 * t82 + t37 * t83) + t253; m(7) * (t73 ^ 2 + t82 ^ 2 + t83 ^ 2) + t253;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
