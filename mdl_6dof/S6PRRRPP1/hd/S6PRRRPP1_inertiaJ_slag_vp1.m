% Calculate joint inertia matrix for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:54
% EndTime: 2019-03-08 22:43:05
% DurationCPUTime: 5.11s
% Computational Cost: add. (23245->587), mult. (49503->824), div. (0->0), fcn. (63373->12), ass. (0->264)
t290 = m(6) + m(7);
t292 = rSges(7,1) + pkin(5);
t291 = rSges(7,3) + qJ(6);
t289 = cos(qJ(3));
t234 = cos(qJ(4));
t288 = pkin(4) * t234;
t229 = cos(pkin(6));
t228 = cos(pkin(10));
t235 = cos(qJ(2));
t277 = t235 * t228;
t226 = sin(pkin(10));
t233 = sin(qJ(2));
t279 = t233 * t226;
t214 = -t229 * t277 + t279;
t231 = sin(qJ(4));
t286 = t214 * t231;
t278 = t235 * t226;
t280 = t228 * t233;
t216 = t229 * t278 + t280;
t285 = t216 * t231;
t227 = sin(pkin(6));
t284 = t226 * t227;
t283 = t227 * t228;
t232 = sin(qJ(3));
t282 = t227 * t232;
t281 = t227 * t235;
t215 = t229 * t280 + t278;
t255 = t227 * t289;
t202 = t215 * t232 + t228 * t255;
t203 = t215 * t289 - t228 * t282;
t108 = pkin(4) * t286 + qJ(5) * t202 + t203 * t288;
t171 = pkin(3) * t203 + pkin(9) * t202;
t161 = t216 * t171;
t276 = t216 * t108 + t161;
t217 = -t229 * t279 + t277;
t204 = t217 * t232 - t226 * t255;
t205 = t217 * t289 + t226 * t282;
t109 = pkin(4) * t285 + qJ(5) * t204 + t205 * t288;
t264 = qJ(4) + pkin(11);
t225 = sin(t264);
t254 = cos(t264);
t175 = t205 * t225 - t216 * t254;
t176 = t205 * t254 + t216 * t225;
t125 = rSges(6,1) * t176 - rSges(6,2) * t175 + rSges(6,3) * t204;
t275 = -t109 - t125;
t173 = t203 * t225 - t214 * t254;
t174 = t203 * t254 + t214 * t225;
t274 = rSges(7,2) * t202 + t173 * t291 + t174 * t292;
t273 = rSges(7,2) * t204 + t175 * t291 + t176 * t292;
t179 = -t205 * t231 + t216 * t234;
t180 = t205 * t234 + t285;
t135 = rSges(5,1) * t180 + rSges(5,2) * t179 + rSges(5,3) * t204;
t172 = pkin(3) * t205 + pkin(9) * t204;
t272 = -t135 - t172;
t219 = t229 * t232 + t233 * t255;
t200 = t219 * t225 + t254 * t281;
t201 = t219 * t254 - t225 * t281;
t218 = -t229 * t289 + t233 * t282;
t271 = rSges(7,2) * t218 + t200 * t291 + t201 * t292;
t148 = rSges(6,1) * t201 - rSges(6,2) * t200 + rSges(6,3) * t218;
t262 = t231 * t281;
t149 = -pkin(4) * t262 + qJ(5) * t218 + t219 * t288;
t270 = -t148 - t149;
t206 = -t219 * t231 - t234 * t281;
t207 = t219 * t234 - t262;
t162 = rSges(5,1) * t207 + rSges(5,2) * t206 + rSges(5,3) * t218;
t199 = t219 * pkin(3) + t218 * pkin(9);
t269 = -t162 - t199;
t268 = t171 * t281 + t214 * t199;
t198 = pkin(2) * t217 + pkin(8) * t216;
t196 = t229 * t198;
t267 = t229 * t172 + t196;
t197 = pkin(2) * t215 + pkin(8) * t214;
t266 = -t171 - t197;
t265 = t197 * t284 + t198 * t283;
t261 = t229 * t109 + t267;
t260 = -t108 + t266;
t259 = -t109 - t273;
t258 = -t172 + t275;
t257 = -t149 - t271;
t256 = -t199 + t270;
t193 = t219 * rSges(4,1) - t218 * rSges(4,2) - rSges(4,3) * t281;
t220 = (pkin(2) * t233 - pkin(8) * t235) * t227;
t253 = (-t193 - t220) * t227;
t110 = Icges(7,5) * t174 + Icges(7,6) * t202 + Icges(7,3) * t173;
t114 = Icges(7,4) * t174 + Icges(7,2) * t202 + Icges(7,6) * t173;
t118 = Icges(7,1) * t174 + Icges(7,4) * t202 + Icges(7,5) * t173;
t50 = t110 * t173 + t114 * t202 + t118 * t174;
t111 = Icges(7,5) * t176 + Icges(7,6) * t204 + Icges(7,3) * t175;
t115 = Icges(7,4) * t176 + Icges(7,2) * t204 + Icges(7,6) * t175;
t119 = Icges(7,1) * t176 + Icges(7,4) * t204 + Icges(7,5) * t175;
t51 = t111 * t173 + t115 * t202 + t119 * t174;
t141 = Icges(7,5) * t201 + Icges(7,6) * t218 + Icges(7,3) * t200;
t143 = Icges(7,4) * t201 + Icges(7,2) * t218 + Icges(7,6) * t200;
t145 = Icges(7,1) * t201 + Icges(7,4) * t218 + Icges(7,5) * t200;
t72 = t141 * t173 + t143 * t202 + t145 * t174;
t1 = t202 * t50 + t204 * t51 + t218 * t72;
t177 = -t203 * t231 + t214 * t234;
t178 = t203 * t234 + t286;
t128 = Icges(5,5) * t178 + Icges(5,6) * t177 + Icges(5,3) * t202;
t130 = Icges(5,4) * t178 + Icges(5,2) * t177 + Icges(5,6) * t202;
t132 = Icges(5,1) * t178 + Icges(5,4) * t177 + Icges(5,5) * t202;
t60 = t128 * t202 + t130 * t177 + t132 * t178;
t129 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t204;
t131 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t204;
t133 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t204;
t61 = t129 * t202 + t131 * t177 + t133 * t178;
t158 = Icges(5,5) * t207 + Icges(5,6) * t206 + Icges(5,3) * t218;
t159 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t218;
t160 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t218;
t78 = t158 * t202 + t159 * t177 + t160 * t178;
t13 = t202 * t60 + t204 * t61 + t218 * t78;
t112 = Icges(6,5) * t174 - Icges(6,6) * t173 + Icges(6,3) * t202;
t116 = Icges(6,4) * t174 - Icges(6,2) * t173 + Icges(6,6) * t202;
t120 = Icges(6,1) * t174 - Icges(6,4) * t173 + Icges(6,5) * t202;
t52 = t112 * t202 - t116 * t173 + t120 * t174;
t113 = Icges(6,5) * t176 - Icges(6,6) * t175 + Icges(6,3) * t204;
t117 = Icges(6,4) * t176 - Icges(6,2) * t175 + Icges(6,6) * t204;
t121 = Icges(6,1) * t176 - Icges(6,4) * t175 + Icges(6,5) * t204;
t53 = t113 * t202 - t117 * t173 + t121 * t174;
t142 = Icges(6,5) * t201 - Icges(6,6) * t200 + Icges(6,3) * t218;
t144 = Icges(6,4) * t201 - Icges(6,2) * t200 + Icges(6,6) * t218;
t146 = Icges(6,1) * t201 - Icges(6,4) * t200 + Icges(6,5) * t218;
t73 = t142 * t202 - t144 * t173 + t146 * t174;
t2 = t202 * t52 + t204 * t53 + t218 * t73;
t252 = -t13 / 0.2e1 - t2 / 0.2e1 - t1 / 0.2e1;
t62 = t128 * t204 + t130 * t179 + t132 * t180;
t63 = t129 * t204 + t131 * t179 + t133 * t180;
t79 = t158 * t204 + t159 * t179 + t160 * t180;
t14 = t202 * t62 + t204 * t63 + t218 * t79;
t54 = t110 * t175 + t114 * t204 + t118 * t176;
t55 = t111 * t175 + t115 * t204 + t119 * t176;
t74 = t141 * t175 + t143 * t204 + t145 * t176;
t3 = t202 * t54 + t204 * t55 + t218 * t74;
t56 = t112 * t204 - t116 * t175 + t120 * t176;
t57 = t113 * t204 - t117 * t175 + t121 * t176;
t75 = t142 * t204 - t144 * t175 + t146 * t176;
t4 = t202 * t56 + t204 * t57 + t218 * t75;
t251 = t14 / 0.2e1 + t4 / 0.2e1 + t3 / 0.2e1;
t15 = t60 * t214 + t61 * t216 - t281 * t78;
t5 = t50 * t214 + t51 * t216 - t281 * t72;
t6 = t52 * t214 + t53 * t216 - t281 * t73;
t250 = t15 / 0.2e1 + t6 / 0.2e1 + t5 / 0.2e1;
t16 = t62 * t214 + t63 * t216 - t281 * t79;
t7 = t54 * t214 + t55 * t216 - t281 * t74;
t8 = t56 * t214 + t57 * t216 - t281 * t75;
t249 = t16 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1;
t248 = t108 * t281 + t214 * t149 + t268;
t247 = -t172 + t259;
t246 = -t199 + t257;
t245 = t171 * t284 + t172 * t283 + t265;
t10 = t229 * t73 + (t226 * t53 - t228 * t52) * t227;
t17 = t229 * t78 + (t226 * t61 - t228 * t60) * t227;
t9 = t229 * t72 + (t226 * t51 - t228 * t50) * t227;
t244 = t10 / 0.2e1 + t9 / 0.2e1 + t17 / 0.2e1;
t11 = t229 * t74 + (t226 * t55 - t228 * t54) * t227;
t12 = t229 * t75 + (t226 * t57 - t228 * t56) * t227;
t18 = t229 * t79 + (t226 * t63 - t228 * t62) * t227;
t243 = t12 / 0.2e1 + t11 / 0.2e1 + t18 / 0.2e1;
t64 = t110 * t200 + t114 * t218 + t118 * t201;
t65 = t111 * t200 + t115 * t218 + t119 * t201;
t82 = t141 * t200 + t143 * t218 + t145 * t201;
t19 = t202 * t64 + t204 * t65 + t218 * t82;
t66 = t112 * t218 - t116 * t200 + t120 * t201;
t67 = t113 * t218 - t117 * t200 + t121 * t201;
t83 = t142 * t218 - t144 * t200 + t146 * t201;
t20 = t202 * t66 + t204 * t67 + t218 * t83;
t69 = t128 * t218 + t130 * t206 + t132 * t207;
t70 = t129 * t218 + t131 * t206 + t133 * t207;
t89 = t158 * t218 + t159 * t206 + t160 * t207;
t25 = t202 * t69 + t204 * t70 + t218 * t89;
t242 = -t20 / 0.2e1 - t19 / 0.2e1 - t25 / 0.2e1;
t21 = t64 * t214 + t65 * t216 - t281 * t82;
t22 = t66 * t214 + t67 * t216 - t281 * t83;
t26 = t69 * t214 + t70 * t216 - t281 * t89;
t241 = t21 / 0.2e1 + t26 / 0.2e1 + t22 / 0.2e1;
t23 = t229 * t82 + (t226 * t65 - t228 * t64) * t227;
t24 = t229 * t83 + (t226 * t67 - t228 * t66) * t227;
t27 = t229 * t89 + (t226 * t70 - t228 * t69) * t227;
t240 = t27 / 0.2e1 + t24 / 0.2e1 + t23 / 0.2e1;
t239 = (-t220 + t269) * t227;
t238 = (-t220 + t256) * t227;
t237 = t108 * t284 + t109 * t283 + t245;
t236 = (-t220 + t246) * t227;
t213 = t229 * rSges(3,3) + (rSges(3,1) * t233 + rSges(3,2) * t235) * t227;
t212 = Icges(3,5) * t229 + (Icges(3,1) * t233 + Icges(3,4) * t235) * t227;
t211 = Icges(3,6) * t229 + (Icges(3,4) * t233 + Icges(3,2) * t235) * t227;
t210 = Icges(3,3) * t229 + (Icges(3,5) * t233 + Icges(3,6) * t235) * t227;
t192 = Icges(4,1) * t219 - Icges(4,4) * t218 - Icges(4,5) * t281;
t191 = Icges(4,4) * t219 - Icges(4,2) * t218 - Icges(4,6) * t281;
t190 = Icges(4,5) * t219 - Icges(4,6) * t218 - Icges(4,3) * t281;
t189 = rSges(3,1) * t217 - rSges(3,2) * t216 + rSges(3,3) * t284;
t188 = rSges(3,1) * t215 - rSges(3,2) * t214 - rSges(3,3) * t283;
t187 = Icges(3,1) * t217 - Icges(3,4) * t216 + Icges(3,5) * t284;
t186 = Icges(3,1) * t215 - Icges(3,4) * t214 - Icges(3,5) * t283;
t185 = Icges(3,4) * t217 - Icges(3,2) * t216 + Icges(3,6) * t284;
t184 = Icges(3,4) * t215 - Icges(3,2) * t214 - Icges(3,6) * t283;
t183 = Icges(3,5) * t217 - Icges(3,6) * t216 + Icges(3,3) * t284;
t182 = Icges(3,5) * t215 - Icges(3,6) * t214 - Icges(3,3) * t283;
t165 = -t188 * t229 - t213 * t283;
t164 = t189 * t229 - t213 * t284;
t157 = rSges(4,1) * t205 - rSges(4,2) * t204 + rSges(4,3) * t216;
t156 = rSges(4,1) * t203 - rSges(4,2) * t202 + rSges(4,3) * t214;
t155 = Icges(4,1) * t205 - Icges(4,4) * t204 + Icges(4,5) * t216;
t154 = Icges(4,1) * t203 - Icges(4,4) * t202 + Icges(4,5) * t214;
t153 = Icges(4,4) * t205 - Icges(4,2) * t204 + Icges(4,6) * t216;
t152 = Icges(4,4) * t203 - Icges(4,2) * t202 + Icges(4,6) * t214;
t151 = Icges(4,5) * t205 - Icges(4,6) * t204 + Icges(4,3) * t216;
t150 = Icges(4,5) * t203 - Icges(4,6) * t202 + Icges(4,3) * t214;
t140 = (t188 * t226 + t189 * t228) * t227;
t138 = t202 * t149;
t134 = rSges(5,1) * t178 + rSges(5,2) * t177 + rSges(5,3) * t202;
t127 = -t157 * t281 - t216 * t193;
t126 = t156 * t281 + t214 * t193;
t123 = rSges(6,1) * t174 - rSges(6,2) * t173 + rSges(6,3) * t202;
t103 = t218 * t109;
t101 = -t190 * t281 - t218 * t191 + t219 * t192;
t100 = t204 * t108;
t99 = t156 * t216 - t157 * t214;
t98 = (-t156 - t197) * t229 + t228 * t253;
t97 = t157 * t229 + t226 * t253 + t196;
t96 = t190 * t216 - t191 * t204 + t192 * t205;
t95 = t190 * t214 - t191 * t202 + t192 * t203;
t94 = (t156 * t226 + t157 * t228) * t227 + t265;
t93 = t135 * t218 - t162 * t204;
t92 = -t134 * t218 + t162 * t202;
t91 = -t151 * t281 - t218 * t153 + t219 * t155;
t90 = -t150 * t281 - t218 * t152 + t219 * t154;
t88 = t151 * t216 - t153 * t204 + t155 * t205;
t87 = t150 * t216 - t152 * t204 + t154 * t205;
t86 = t151 * t214 - t153 * t202 + t155 * t203;
t85 = t150 * t214 - t152 * t202 + t154 * t203;
t84 = t134 * t204 - t135 * t202;
t81 = t216 * t269 + t272 * t281;
t80 = t134 * t281 + t214 * t162 + t268;
t77 = (-t134 + t266) * t229 + t228 * t239;
t76 = t135 * t229 + t226 * t239 + t267;
t71 = t134 * t216 + t214 * t272 + t161;
t68 = (t134 * t226 + t135 * t228) * t227 + t245;
t59 = t125 * t218 + t204 * t270 + t103;
t58 = t148 * t202 + t138 + (-t108 - t123) * t218;
t49 = t216 * t256 + t258 * t281;
t48 = t123 * t281 + t214 * t148 + t248;
t47 = (-t123 + t260) * t229 + t228 * t238;
t46 = t125 * t229 + t226 * t238 + t261;
t45 = t123 * t204 + t202 * t275 + t100;
t44 = t101 * t229 + (t226 * t91 - t228 * t90) * t227;
t43 = -t101 * t281 + t90 * t214 + t91 * t216;
t42 = t204 * t257 + t218 * t273 + t103;
t41 = t138 + t271 * t202 + (-t108 - t274) * t218;
t40 = t123 * t216 + t214 * t258 + t276;
t39 = (t123 * t226 + t125 * t228) * t227 + t237;
t38 = t216 * t246 + t247 * t281;
t37 = t214 * t271 + t274 * t281 + t248;
t36 = (t260 - t274) * t229 + t228 * t236;
t35 = t226 * t236 + t229 * t273 + t261;
t34 = t229 * t96 + (t226 * t88 - t228 * t87) * t227;
t33 = t229 * t95 + (t226 * t86 - t228 * t85) * t227;
t32 = t87 * t214 + t88 * t216 - t281 * t96;
t31 = t85 * t214 + t86 * t216 - t281 * t95;
t30 = t202 * t259 + t204 * t274 + t100;
t29 = t214 * t247 + t216 * t274 + t276;
t28 = (t226 * t274 + t228 * t273) * t227 + t237;
t102 = [m(2) + m(3) + m(4) + m(5) + t290; m(3) * t140 + m(4) * t94 + m(5) * t68 + m(6) * t39 + m(7) * t28; m(7) * (t28 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(6) * (t39 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(5) * (t68 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(4) * (t94 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(3) * (t140 ^ 2 + t164 ^ 2 + t165 ^ 2) + (t18 + t12 + t11 + t34 + (t183 * t284 - t216 * t185 + t217 * t187) * t284) * t284 + (-t17 - t10 - t9 - t33 + (-t182 * t283 - t184 * t214 + t186 * t215) * t283 + (-t182 * t284 + t183 * t283 + t184 * t216 + t185 * t214 - t186 * t217 - t187 * t215) * t284) * t283 + ((t210 * t284 - t211 * t216 + t212 * t217) * t284 - (-t210 * t283 - t211 * t214 + t212 * t215) * t283 + t23 + t27 + t24 + t44 + ((t185 * t235 + t187 * t233) * t226 - (t184 * t235 + t186 * t233) * t228) * t227 ^ 2 + ((-t182 * t228 + t183 * t226 + t211 * t235 + t212 * t233) * t227 + t229 * t210) * t229) * t229; m(4) * t99 + m(5) * t71 + m(6) * t40 + m(7) * t29; (t43 / 0.2e1 + t241) * t229 + (t34 / 0.2e1 + t243) * t216 + (t33 / 0.2e1 + t244) * t214 + m(7) * (t28 * t29 + t35 * t38 + t36 * t37) + m(6) * (t39 * t40 + t46 * t49 + t47 * t48) + m(5) * (t68 * t71 + t76 * t81 + t77 * t80) + m(4) * (t126 * t98 + t127 * t97 + t94 * t99) + ((-t44 / 0.2e1 - t240) * t235 + (-t31 / 0.2e1 - t250) * t228 + (t32 / 0.2e1 + t249) * t226) * t227; (-t21 - t22 - t26 - t43) * t281 + (t16 + t8 + t7 + t32) * t216 + (t5 + t6 + t15 + t31) * t214 + m(7) * (t29 ^ 2 + t37 ^ 2 + t38 ^ 2) + m(6) * (t40 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t71 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2 + t99 ^ 2); m(5) * t84 + m(6) * t45 + m(7) * t30; -t242 * t229 + t240 * t218 + t243 * t204 + t244 * t202 + m(7) * (t28 * t30 + t35 * t42 + t36 * t41) + m(6) * (t39 * t45 + t46 * t59 + t47 * t58) + m(5) * (t68 * t84 + t76 * t93 + t77 * t92) + (t226 * t251 + t228 * t252) * t227; t242 * t281 + t241 * t218 + t251 * t216 - t252 * t214 + t249 * t204 + t250 * t202 + m(7) * (t29 * t30 + t37 * t41 + t38 * t42) + m(6) * (t40 * t45 + t48 * t58 + t49 * t59) + m(5) * (t71 * t84 + t80 * t92 + t81 * t93); (t25 + t20 + t19) * t218 + (t3 + t14 + t4) * t204 + (t2 + t13 + t1) * t202 + m(7) * (t30 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t45 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(5) * (t84 ^ 2 + t92 ^ 2 + t93 ^ 2); t218 * t290; m(7) * (t202 * t35 + t204 * t36 + t218 * t28) + m(6) * (t202 * t46 + t204 * t47 + t218 * t39); m(7) * (t202 * t38 + t204 * t37 + t218 * t29) + m(6) * (t202 * t49 + t204 * t48 + t218 * t40); m(7) * (t202 * t42 + t204 * t41 + t218 * t30) + m(6) * (t202 * t59 + t204 * t58 + t218 * t45); (t202 ^ 2 + t204 ^ 2 + t218 ^ 2) * t290; m(7) * t200; m(7) * (t173 * t35 + t175 * t36 + t200 * t28); m(7) * (t173 * t38 + t175 * t37 + t200 * t29); m(7) * (t173 * t42 + t175 * t41 + t200 * t30); m(7) * (t173 * t202 + t175 * t204 + t200 * t218); m(7) * (t173 ^ 2 + t175 ^ 2 + t200 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t102(1) t102(2) t102(4) t102(7) t102(11) t102(16); t102(2) t102(3) t102(5) t102(8) t102(12) t102(17); t102(4) t102(5) t102(6) t102(9) t102(13) t102(18); t102(7) t102(8) t102(9) t102(10) t102(14) t102(19); t102(11) t102(12) t102(13) t102(14) t102(15) t102(20); t102(16) t102(17) t102(18) t102(19) t102(20) t102(21);];
Mq  = res;
