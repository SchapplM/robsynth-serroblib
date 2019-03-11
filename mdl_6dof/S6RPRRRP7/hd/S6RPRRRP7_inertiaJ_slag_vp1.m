% Calculate joint inertia matrix for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:40
% EndTime: 2019-03-09 06:18:51
% DurationCPUTime: 4.69s
% Computational Cost: add. (11601->425), mult. (11058->612), div. (0->0), fcn. (12034->10), ass. (0->220)
t210 = pkin(10) + qJ(3);
t206 = cos(t210);
t213 = qJ(4) + qJ(5);
t208 = cos(t213);
t220 = cos(qJ(1));
t268 = t208 * t220;
t207 = sin(t213);
t218 = sin(qJ(1));
t270 = t207 * t218;
t172 = t206 * t270 + t268;
t263 = t218 * t208;
t269 = t207 * t220;
t173 = t206 * t263 - t269;
t205 = sin(t210);
t275 = t205 * t218;
t101 = Icges(7,4) * t173 + Icges(7,2) * t275 + Icges(7,6) * t172;
t99 = Icges(6,5) * t173 - Icges(6,6) * t172 + Icges(6,3) * t275;
t320 = t101 + t99;
t103 = Icges(6,4) * t173 - Icges(6,2) * t172 + Icges(6,6) * t275;
t97 = Icges(7,5) * t173 + Icges(7,6) * t275 + Icges(7,3) * t172;
t319 = -t103 + t97;
t174 = t206 * t269 - t263;
t175 = t206 * t268 + t270;
t274 = t205 * t220;
t104 = Icges(6,4) * t175 - Icges(6,2) * t174 + Icges(6,6) * t274;
t98 = Icges(7,5) * t175 + Icges(7,6) * t274 + Icges(7,3) * t174;
t318 = -t104 + t98;
t100 = Icges(6,5) * t175 - Icges(6,6) * t174 + Icges(6,3) * t274;
t102 = Icges(7,4) * t175 + Icges(7,2) * t274 + Icges(7,6) * t174;
t317 = t100 + t102;
t105 = Icges(7,1) * t173 + Icges(7,4) * t275 + Icges(7,5) * t172;
t107 = Icges(6,1) * t173 - Icges(6,4) * t172 + Icges(6,5) * t275;
t316 = t105 + t107;
t106 = Icges(7,1) * t175 + Icges(7,4) * t274 + Icges(7,5) * t174;
t108 = Icges(6,1) * t175 - Icges(6,4) * t174 + Icges(6,5) * t274;
t315 = t106 + t108;
t294 = rSges(7,3) + qJ(6);
t302 = rSges(7,1) + pkin(5);
t314 = -t294 * t172 - t302 * t173;
t313 = t319 * t172 + t316 * t173 + t320 * t275;
t312 = t318 * t172 + t315 * t173 + t317 * t275;
t311 = t319 * t174 + t316 * t175 + t320 * t274;
t310 = t318 * t174 + t315 * t175 + t317 * t274;
t144 = -Icges(7,6) * t206 + (Icges(7,5) * t208 + Icges(7,3) * t207) * t205;
t146 = -Icges(7,2) * t206 + (Icges(7,4) * t208 + Icges(7,6) * t207) * t205;
t148 = -Icges(7,4) * t206 + (Icges(7,1) * t208 + Icges(7,5) * t207) * t205;
t69 = t144 * t172 + t146 * t275 + t148 * t173;
t145 = -Icges(6,3) * t206 + (Icges(6,5) * t208 - Icges(6,6) * t207) * t205;
t147 = -Icges(6,6) * t206 + (Icges(6,4) * t208 - Icges(6,2) * t207) * t205;
t149 = -Icges(6,5) * t206 + (Icges(6,1) * t208 - Icges(6,4) * t207) * t205;
t70 = t145 * t275 - t147 * t172 + t149 * t173;
t309 = -t70 - t69;
t71 = t144 * t174 + t146 * t274 + t148 * t175;
t72 = t145 * t274 - t147 * t174 + t149 * t175;
t308 = -t71 - t72;
t271 = t207 * t147;
t277 = t205 * t207;
t304 = t144 * t277 + (t148 + t149) * t205 * t208;
t306 = -t145 - t146;
t295 = -t205 * t271 + t306 * t206 + t304;
t307 = t295 * t206;
t305 = Icges(4,5) * t205;
t303 = t305 / 0.2e1;
t301 = t309 * t206 + (t313 * t218 + t312 * t220) * t205;
t300 = t308 * t206 + (t311 * t218 + t310 * t220) * t205;
t299 = t312 * t218 - t313 * t220;
t298 = t310 * t218 - t311 * t220;
t50 = -t101 * t206 + (t105 * t208 + t207 * t97) * t205;
t52 = -t206 * t99 + (-t103 * t207 + t107 * t208) * t205;
t297 = -t50 - t52;
t51 = -t102 * t206 + (t106 * t208 + t207 * t98) * t205;
t53 = -t100 * t206 + (-t104 * t207 + t108 * t208) * t205;
t296 = t51 + t53;
t293 = rSges(7,2) * t275 - t314;
t292 = -rSges(7,2) * t206 + (t294 * t207 + t302 * t208) * t205;
t211 = t218 ^ 2;
t212 = t220 ^ 2;
t291 = -t206 / 0.2e1;
t290 = t218 / 0.2e1;
t289 = -t220 / 0.2e1;
t288 = t220 / 0.2e1;
t186 = rSges(4,1) * t205 + rSges(4,2) * t206;
t287 = m(4) * t186;
t286 = pkin(3) * t206;
t285 = pkin(8) * t205;
t219 = cos(qJ(4));
t204 = pkin(4) * t219 + pkin(3);
t284 = -pkin(3) + t204;
t283 = t307 + (t218 * t297 - t220 * t296) * t205;
t281 = rSges(3,3) + qJ(2);
t280 = t293 * t274;
t234 = -rSges(6,1) * t173 + rSges(6,2) * t172;
t110 = rSges(6,3) * t275 - t234;
t151 = -rSges(6,3) * t206 + (rSges(6,1) * t208 - rSges(6,2) * t207) * t205;
t84 = t206 * t110 + t151 * t275;
t278 = Icges(4,4) * t206;
t221 = -pkin(9) - pkin(8);
t273 = t205 * t221;
t272 = t206 * t220;
t216 = -pkin(7) - qJ(2);
t267 = t216 * t220;
t217 = sin(qJ(4));
t155 = -Icges(5,6) * t206 + (Icges(5,4) * t219 - Icges(5,2) * t217) * t205;
t266 = t217 * t155;
t265 = t217 * t218;
t264 = t217 * t220;
t262 = t218 * t219;
t261 = t219 * t220;
t260 = rSges(7,2) * t274 + t294 * t174 + t175 * t302;
t112 = t175 * rSges(6,1) - t174 * rSges(6,2) + rSges(6,3) * t274;
t226 = pkin(4) * t265 + t204 * t272 - t220 * t273;
t251 = pkin(3) * t272 + pkin(8) * t274;
t120 = t226 - t251;
t259 = -t112 - t120;
t252 = -pkin(4) * t264 - t218 * t273;
t119 = (t284 * t206 - t285) * t218 + t252;
t142 = (pkin(8) + t221) * t206 + t284 * t205;
t258 = t206 * t119 + t142 * t275;
t256 = -t142 - t151;
t157 = -rSges(5,3) * t206 + (rSges(5,1) * t219 - rSges(5,2) * t217) * t205;
t188 = t205 * pkin(3) - t206 * pkin(8);
t254 = -t157 - t188;
t253 = t211 * (t285 + t286) + t220 * t251;
t250 = t211 + t212;
t178 = -t206 * t265 - t261;
t179 = t206 * t262 - t264;
t121 = Icges(5,5) * t179 + Icges(5,6) * t178 + Icges(5,3) * t275;
t123 = Icges(5,4) * t179 + Icges(5,2) * t178 + Icges(5,6) * t275;
t125 = Icges(5,1) * t179 + Icges(5,4) * t178 + Icges(5,5) * t275;
t60 = -t121 * t206 + (-t123 * t217 + t125 * t219) * t205;
t154 = -Icges(5,3) * t206 + (Icges(5,5) * t219 - Icges(5,6) * t217) * t205;
t156 = -Icges(5,5) * t206 + (Icges(5,1) * t219 - Icges(5,4) * t217) * t205;
t74 = t154 * t275 + t155 * t178 + t156 * t179;
t249 = t74 / 0.2e1 + t60 / 0.2e1;
t180 = -t206 * t264 + t262;
t181 = t206 * t261 + t265;
t122 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t274;
t124 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t274;
t126 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t274;
t61 = -t122 * t206 + (-t124 * t217 + t126 * t219) * t205;
t75 = t154 * t274 + t155 * t180 + t156 * t181;
t248 = t75 / 0.2e1 + t61 / 0.2e1;
t247 = t274 * t300 + t275 * t301;
t246 = -t120 - t260;
t245 = -t142 - t292;
t244 = -t188 + t256;
t128 = t181 * rSges(5,1) + t180 * rSges(5,2) + rSges(5,3) * t274;
t243 = t275 / 0.2e1;
t242 = t274 / 0.2e1;
t215 = cos(pkin(10));
t203 = pkin(2) * t215 + pkin(1);
t241 = t220 * t203 - t218 * t216;
t240 = -t204 * t206 - t203;
t58 = t206 * t293 + t275 * t292;
t239 = t218 * t119 + t220 * t120 + t253;
t238 = -t188 + t245;
t237 = -t252 - t267;
t236 = rSges(4,1) * t206 - rSges(4,2) * t205;
t235 = -rSges(5,1) * t179 - rSges(5,2) * t178;
t232 = -Icges(4,2) * t205 + t278;
t231 = Icges(4,5) * t206 - Icges(4,6) * t205;
t228 = rSges(4,1) * t272 - rSges(4,2) * t274 + t218 * rSges(4,3);
t214 = sin(pkin(10));
t227 = rSges(3,1) * t215 - rSges(3,2) * t214 + pkin(1);
t225 = t283 * t206 + t247;
t224 = (-t297 - t309) * t243 + (t296 - t308) * t242;
t223 = (t218 * t296 + t220 * t297) * t291 + t300 * t290 + t301 * t289 + t299 * t243 + t298 * t242;
t222 = t226 + t241;
t191 = rSges(2,1) * t220 - rSges(2,2) * t218;
t190 = -rSges(2,1) * t218 - rSges(2,2) * t220;
t183 = Icges(4,6) * t206 + t305;
t159 = Icges(4,3) * t218 + t231 * t220;
t158 = -Icges(4,3) * t220 + t231 * t218;
t153 = t281 * t218 + t227 * t220;
t152 = -t227 * t218 + t281 * t220;
t141 = t205 * t219 * t156;
t138 = t228 + t241;
t137 = (rSges(4,3) - t216) * t220 + (-t203 - t236) * t218;
t130 = t254 * t220;
t129 = t254 * t218;
t127 = rSges(5,3) * t275 - t235;
t115 = t220 * t228 + (-t220 * rSges(4,3) + t236 * t218) * t218;
t96 = t119 * t274;
t93 = t110 * t274;
t91 = t241 + t128 + t251;
t90 = -t267 + (-t286 - t203 + (-rSges(5,3) - pkin(8)) * t205) * t218 + t235;
t89 = t244 * t220;
t88 = t244 * t218;
t87 = -t128 * t206 - t157 * t274;
t86 = t127 * t206 + t157 * t275;
t85 = -t112 * t206 - t151 * t274;
t83 = -t206 * t154 - t205 * t266 + t141;
t82 = t222 + t112;
t81 = (-rSges(6,3) * t205 + t240) * t218 + t234 + t237;
t80 = t238 * t220;
t79 = t238 * t218;
t78 = (t127 * t220 - t128 * t218) * t205;
t73 = -t112 * t275 + t93;
t68 = t127 * t218 + t128 * t220 + t253;
t67 = t222 + t260;
t66 = (-rSges(7,2) * t205 + t240) * t218 + t237 + t314;
t59 = -t260 * t206 - t274 * t292;
t57 = t122 * t274 + t124 * t180 + t126 * t181;
t56 = t121 * t274 + t123 * t180 + t125 * t181;
t55 = t122 * t275 + t124 * t178 + t126 * t179;
t54 = t121 * t275 + t123 * t178 + t125 * t179;
t49 = t259 * t206 + t256 * t274;
t48 = t258 + t84;
t35 = -t260 * t275 + t280;
t34 = t259 * t275 + t93 + t96;
t33 = t110 * t218 + t112 * t220 + t239;
t32 = t246 * t206 + t245 * t274;
t31 = t58 + t258;
t30 = t246 * t275 + t280 + t96;
t29 = t218 * t293 + t260 * t220 + t239;
t28 = t218 * t57 - t220 * t56;
t27 = t218 * t55 - t220 * t54;
t16 = -t206 * t75 + (t218 * t56 + t220 * t57) * t205;
t15 = -t206 * t74 + (t218 * t54 + t220 * t55) * t205;
t1 = [Icges(3,2) * t215 ^ 2 + Icges(2,3) + t141 + (Icges(3,1) * t214 + 0.2e1 * Icges(3,4) * t215) * t214 + (Icges(4,1) * t205 - t266 - t271 + t278) * t205 + (Icges(4,4) * t205 + Icges(4,2) * t206 - t154 + t306) * t206 + m(7) * (t66 ^ 2 + t67 ^ 2) + m(6) * (t81 ^ 2 + t82 ^ 2) + m(5) * (t90 ^ 2 + t91 ^ 2) + m(4) * (t137 ^ 2 + t138 ^ 2) + m(3) * (t152 ^ 2 + t153 ^ 2) + m(2) * (t190 ^ 2 + t191 ^ 2) + t304; m(7) * (t218 * t66 - t220 * t67) + m(6) * (t218 * t81 - t220 * t82) + m(5) * (t218 * t90 - t220 * t91) + m(4) * (t137 * t218 - t138 * t220) + m(3) * (t152 * t218 - t153 * t220); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t250; m(7) * (t66 * t80 + t67 * t79) + m(6) * (t81 * t89 + t82 * t88) + m(5) * (t129 * t91 + t130 * t90) + (-t70 / 0.2e1 - t69 / 0.2e1 - t52 / 0.2e1 - t50 / 0.2e1 + (-Icges(4,6) * t220 + t232 * t218) * t291 + t220 * t303 - t137 * t287 + t183 * t288 - t249) * t220 + (t53 / 0.2e1 + t51 / 0.2e1 + t72 / 0.2e1 + t71 / 0.2e1 + (Icges(4,6) * t218 + t232 * t220) * t206 / 0.2e1 + t218 * t303 - t138 * t287 + t183 * t290 + t248) * t218; m(5) * (-t129 * t220 + t130 * t218) + m(6) * (t218 * t89 - t220 * t88) + m(7) * (t218 * t80 - t220 * t79); m(7) * (t29 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t33 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(5) * (t129 ^ 2 + t130 ^ 2 + t68 ^ 2) + m(4) * (t250 * t186 ^ 2 + t115 ^ 2) + (-t212 * t158 - t27 - t299) * t220 + (t211 * t159 + t28 + (-t218 * t158 + t220 * t159) * t220 + t298) * t218; (-t83 - t295) * t206 + m(7) * (t31 * t66 + t32 * t67) + m(6) * (t48 * t81 + t49 * t82) + m(5) * (t86 * t90 + t87 * t91) + (t249 * t218 + t248 * t220) * t205 + t224; m(5) * (t218 * t86 - t220 * t87) + m(6) * (t218 * t48 - t220 * t49) + m(7) * (t218 * t31 - t220 * t32); (t61 * t218 - t60 * t220) * t291 + m(7) * (t30 * t29 + t31 * t80 + t32 * t79) + m(6) * (t34 * t33 + t48 * t89 + t49 * t88) + m(5) * (t129 * t87 + t130 * t86 + t68 * t78) + (t27 * t290 + t28 * t288) * t205 + t223 + t15 * t289 + t16 * t290; (t83 * t206 + t283) * t206 + (t220 * t16 + t218 * t15 - t206 * (t218 * t60 + t220 * t61)) * t205 + m(7) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t34 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t78 ^ 2 + t86 ^ 2 + t87 ^ 2) + t247; -t307 + m(7) * (t58 * t66 + t59 * t67) + m(6) * (t81 * t84 + t82 * t85) + t224; m(6) * (t218 * t84 - t220 * t85) + m(7) * (t218 * t58 - t220 * t59); m(7) * (t35 * t29 + t58 * t80 + t59 * t79) + m(6) * (t33 * t73 + t84 * t89 + t85 * t88) + t223; m(7) * (t30 * t35 + t31 * t58 + t32 * t59) + m(6) * (t34 * t73 + t48 * t84 + t49 * t85) + t225; m(7) * (t35 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(6) * (t73 ^ 2 + t84 ^ 2 + t85 ^ 2) + t225; m(7) * (t172 * t67 + t174 * t66); m(7) * (-t172 * t220 + t174 * t218); m(7) * (t172 * t79 + t174 * t80 + t29 * t277); m(7) * (t172 * t32 + t174 * t31 + t30 * t277); m(7) * (t172 * t59 + t174 * t58 + t35 * t277); m(7) * (t205 ^ 2 * t207 ^ 2 + t172 ^ 2 + t174 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
