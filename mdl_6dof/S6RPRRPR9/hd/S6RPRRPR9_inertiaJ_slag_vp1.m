% Calculate joint inertia matrix for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:27:00
% EndTime: 2019-03-09 05:27:14
% DurationCPUTime: 5.60s
% Computational Cost: add. (35699->571), mult. (85654->790), div. (0->0), fcn. (112336->16), ass. (0->270)
t249 = sin(qJ(1));
t244 = cos(pkin(6));
t252 = cos(qJ(1));
t307 = sin(pkin(12));
t273 = t252 * t307;
t309 = cos(pkin(12));
t276 = t249 * t309;
t260 = t244 * t276 + t273;
t243 = sin(pkin(6));
t310 = cos(pkin(7));
t279 = t243 * t310;
t308 = sin(pkin(7));
t220 = t249 * t279 + t260 * t308;
t274 = t252 * t309;
t275 = t249 * t307;
t261 = -t244 * t274 + t275;
t219 = -t252 * t279 + t261 * t308;
t322 = m(7) / 0.2e1;
t323 = m(6) / 0.2e1;
t289 = t323 + t322;
t327 = 0.2e1 * t289;
t326 = m(3) / 0.2e1;
t325 = m(4) / 0.2e1;
t324 = m(5) / 0.2e1;
t228 = t244 * t273 + t276;
t248 = sin(qJ(3));
t259 = t261 * t310;
t278 = t243 * t308;
t315 = cos(qJ(3));
t210 = t228 * t315 + (-t252 * t278 - t259) * t248;
t290 = qJ(4) + pkin(13);
t240 = sin(t290);
t272 = cos(t290);
t180 = t210 * t240 - t219 * t272;
t229 = -t244 * t275 + t274;
t257 = t260 * t310;
t212 = t229 * t315 + (t249 * t278 - t257) * t248;
t182 = t212 * t240 - t220 * t272;
t268 = t310 * t309;
t218 = t244 * t308 * t248 + (t248 * t268 + t307 * t315) * t243;
t227 = t244 * t310 - t278 * t309;
t196 = t218 * t240 - t227 * t272;
t181 = t210 * t272 + t219 * t240;
t269 = t315 * t308;
t265 = t243 * t269;
t209 = t228 * t248 + t252 * t265 + t259 * t315;
t246 = sin(qJ(6));
t250 = cos(qJ(6));
t146 = -t181 * t246 + t209 * t250;
t147 = t181 * t250 + t209 * t246;
t83 = Icges(7,5) * t147 + Icges(7,6) * t146 + Icges(7,3) * t180;
t85 = Icges(7,4) * t147 + Icges(7,2) * t146 + Icges(7,6) * t180;
t87 = Icges(7,1) * t147 + Icges(7,4) * t146 + Icges(7,5) * t180;
t28 = t146 * t85 + t147 * t87 + t180 * t83;
t183 = t212 * t272 + t220 * t240;
t211 = t229 * t248 - t249 * t265 + t257 * t315;
t148 = -t183 * t246 + t211 * t250;
t149 = t183 * t250 + t211 * t246;
t84 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t182;
t86 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t182;
t88 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t182;
t29 = t146 * t86 + t147 * t88 + t180 * t84;
t197 = t218 * t272 + t227 * t240;
t277 = t243 * t307;
t217 = -t243 * t268 * t315 - t244 * t269 + t248 * t277;
t173 = -t197 * t246 + t217 * t250;
t174 = t197 * t250 + t217 * t246;
t108 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t196;
t109 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t196;
t110 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t196;
t36 = t108 * t180 + t109 * t146 + t110 * t147;
t1 = t180 * t28 + t182 * t29 + t196 * t36;
t321 = t1 / 0.2e1;
t30 = t148 * t85 + t149 * t87 + t182 * t83;
t31 = t148 * t86 + t149 * t88 + t182 * t84;
t37 = t108 * t182 + t109 * t148 + t110 * t149;
t2 = t180 * t30 + t182 * t31 + t196 * t37;
t320 = t2 / 0.2e1;
t32 = t173 * t85 + t174 * t87 + t196 * t83;
t33 = t173 * t86 + t174 * t88 + t196 * t84;
t43 = t196 * t108 + t173 * t109 + t174 * t110;
t40 = t43 * t196;
t7 = t32 * t180 + t33 * t182 + t40;
t319 = t7 / 0.2e1;
t318 = t180 / 0.2e1;
t317 = t182 / 0.2e1;
t316 = t196 / 0.2e1;
t314 = t181 * pkin(5);
t251 = cos(qJ(4));
t239 = pkin(4) * t251 + pkin(3);
t313 = -pkin(3) + t239;
t266 = -t147 * rSges(7,1) - t146 * rSges(7,2);
t89 = rSges(7,3) * t180 - t266;
t312 = pkin(11) * t180 + t314 + t89;
t90 = t149 * rSges(7,1) + t148 * rSges(7,2) + t182 * rSges(7,3);
t311 = t183 * pkin(5) + t182 * pkin(11) + t90;
t245 = -qJ(5) - pkin(10);
t306 = t209 * t245;
t247 = sin(qJ(4));
t305 = t219 * t247;
t304 = t220 * t247;
t303 = t227 * t247;
t302 = t243 * t249;
t301 = t243 * t252;
t203 = t209 * pkin(10);
t288 = pkin(4) * t305;
t106 = t210 * t313 - t203 + t288 - t306;
t170 = pkin(3) * t210 + t203;
t164 = t220 * t170;
t300 = t220 * t106 + t164;
t171 = t212 * pkin(3) + t211 * pkin(10);
t281 = pkin(4) * t304 - t211 * t245 + t212 * t239;
t107 = -t171 + t281;
t166 = t227 * t171;
t299 = t227 * t107 + t166;
t267 = -t181 * rSges(6,1) + t180 * rSges(6,2);
t118 = rSges(6,3) * t209 - t267;
t298 = -t106 - t118;
t119 = t183 * rSges(6,1) - t182 * rSges(6,2) + t211 * rSges(6,3);
t297 = -t107 - t119;
t111 = rSges(7,1) * t174 + rSges(7,2) * t173 + rSges(7,3) * t196;
t296 = pkin(5) * t197 + pkin(11) * t196 + t111;
t186 = -t210 * t247 + t219 * t251;
t187 = t210 * t251 + t305;
t126 = rSges(5,1) * t187 + rSges(5,2) * t186 + rSges(5,3) * t209;
t295 = -t126 - t170;
t188 = -t212 * t247 + t220 * t251;
t189 = t212 * t251 + t304;
t127 = t189 * rSges(5,1) + t188 * rSges(5,2) + t211 * rSges(5,3);
t294 = -t127 - t171;
t140 = pkin(4) * t303 + t313 * t218 + (-pkin(10) - t245) * t217;
t194 = t218 * pkin(3) + t217 * pkin(10);
t179 = t219 * t194;
t293 = t219 * t140 + t179;
t150 = rSges(6,1) * t197 - rSges(6,2) * t196 + rSges(6,3) * t217;
t292 = -t140 - t150;
t291 = t252 * pkin(1) + qJ(2) * t302;
t287 = t36 / 0.2e1 + t32 / 0.2e1;
t286 = t37 / 0.2e1 + t33 / 0.2e1;
t285 = -t106 - t312;
t284 = -t107 - t311;
t283 = -t140 - t296;
t143 = Icges(6,5) * t197 - Icges(6,6) * t196 + Icges(6,3) * t217;
t144 = Icges(6,4) * t197 - Icges(6,2) * t196 + Icges(6,6) * t217;
t145 = Icges(6,1) * t197 - Icges(6,4) * t196 + Icges(6,5) * t217;
t72 = t217 * t143 - t196 * t144 + t197 * t145;
t207 = -t218 * t247 + t227 * t251;
t208 = t218 * t251 + t303;
t152 = Icges(5,5) * t208 + Icges(5,6) * t207 + Icges(5,3) * t217;
t153 = Icges(5,4) * t208 + Icges(5,2) * t207 + Icges(5,6) * t217;
t154 = Icges(5,1) * t208 + Icges(5,4) * t207 + Icges(5,5) * t217;
t76 = t217 * t152 + t207 * t153 + t208 * t154;
t190 = Icges(4,5) * t218 - Icges(4,6) * t217 + Icges(4,3) * t227;
t191 = Icges(4,4) * t218 - Icges(4,2) * t217 + Icges(4,6) * t227;
t192 = Icges(4,1) * t218 - Icges(4,4) * t217 + Icges(4,5) * t227;
t282 = t227 * t190 - t217 * t191 + t218 * t192;
t163 = t212 * rSges(4,1) - t211 * rSges(4,2) + t220 * rSges(4,3);
t280 = -t249 * pkin(1) + qJ(2) * t301;
t162 = rSges(4,1) * t210 - rSges(4,2) * t209 + rSges(4,3) * t219;
t264 = -t228 * pkin(2) - t219 * pkin(9) + t280;
t112 = Icges(6,5) * t181 - Icges(6,6) * t180 + Icges(6,3) * t209;
t114 = Icges(6,4) * t181 - Icges(6,2) * t180 + Icges(6,6) * t209;
t116 = Icges(6,1) * t181 - Icges(6,4) * t180 + Icges(6,5) * t209;
t55 = t112 * t217 - t114 * t196 + t116 * t197;
t120 = Icges(5,5) * t187 + Icges(5,6) * t186 + Icges(5,3) * t209;
t122 = Icges(5,4) * t187 + Icges(5,2) * t186 + Icges(5,6) * t209;
t124 = Icges(5,1) * t187 + Icges(5,4) * t186 + Icges(5,5) * t209;
t57 = t120 * t217 + t122 * t207 + t124 * t208;
t64 = t143 * t209 - t144 * t180 + t145 * t181;
t66 = t152 * t209 + t153 * t186 + t154 * t187;
t263 = t66 / 0.2e1 + t64 / 0.2e1 + t57 / 0.2e1 + t55 / 0.2e1 + t287;
t113 = Icges(6,5) * t183 - Icges(6,6) * t182 + Icges(6,3) * t211;
t115 = Icges(6,4) * t183 - Icges(6,2) * t182 + Icges(6,6) * t211;
t117 = Icges(6,1) * t183 - Icges(6,4) * t182 + Icges(6,5) * t211;
t56 = t113 * t217 - t115 * t196 + t117 * t197;
t121 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t211;
t123 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t211;
t125 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t211;
t58 = t121 * t217 + t123 * t207 + t125 * t208;
t65 = t143 * t211 - t144 * t182 + t145 * t183;
t67 = t152 * t211 + t153 * t188 + t154 * t189;
t262 = t58 / 0.2e1 + t56 / 0.2e1 + t67 / 0.2e1 + t65 / 0.2e1 + t286;
t255 = -t210 * t239 + t264 - t288;
t254 = t229 * pkin(2) + t220 * pkin(9) + t291;
t253 = t254 + t281;
t235 = rSges(2,1) * t252 - t249 * rSges(2,2);
t234 = -t249 * rSges(2,1) - rSges(2,2) * t252;
t206 = t229 * rSges(3,1) - rSges(3,2) * t260 + rSges(3,3) * t302 + t291;
t205 = -t228 * rSges(3,1) + rSges(3,2) * t261 + rSges(3,3) * t301 + t280;
t193 = rSges(4,1) * t218 - rSges(4,2) * t217 + rSges(4,3) * t227;
t161 = Icges(4,1) * t212 - Icges(4,4) * t211 + Icges(4,5) * t220;
t160 = Icges(4,1) * t210 - Icges(4,4) * t209 + Icges(4,5) * t219;
t159 = Icges(4,4) * t212 - Icges(4,2) * t211 + Icges(4,6) * t220;
t158 = Icges(4,4) * t210 - Icges(4,2) * t209 + Icges(4,6) * t219;
t157 = Icges(4,5) * t212 - Icges(4,6) * t211 + Icges(4,3) * t220;
t156 = Icges(4,5) * t210 - Icges(4,6) * t209 + Icges(4,3) * t219;
t155 = rSges(5,1) * t208 + rSges(5,2) * t207 + rSges(5,3) * t217;
t134 = t254 + t163;
t133 = -t162 + t264;
t130 = t209 * t140;
t105 = t163 * t227 - t193 * t220;
t104 = -t162 * t227 + t193 * t219;
t101 = t217 * t107;
t100 = t211 * t106;
t98 = t162 * t220 - t163 * t219;
t97 = t282 * t227;
t94 = t190 * t220 - t191 * t211 + t192 * t212;
t93 = t190 * t219 - t191 * t209 + t192 * t210;
t92 = t254 - t294;
t91 = t264 + t295;
t82 = t253 + t119;
t81 = (-rSges(6,3) + t245) * t209 + t255 + t267;
t80 = t127 * t217 - t155 * t211;
t79 = -t126 * t217 + t155 * t209;
t78 = t157 * t227 - t159 * t217 + t161 * t218;
t77 = t156 * t227 - t158 * t217 + t160 * t218;
t75 = t126 * t211 - t127 * t209;
t74 = t76 * t227;
t73 = t76 * t217;
t71 = t72 * t227;
t70 = t72 * t217;
t69 = t227 * t127 + t166 + (-t155 - t194) * t220;
t68 = t219 * t155 + t227 * t295 + t179;
t63 = t253 + t311;
t62 = -t314 + t306 + (-rSges(7,3) - pkin(11)) * t180 + t255 + t266;
t61 = t220 * t126 + t219 * t294 + t164;
t60 = -t111 * t182 + t196 * t90;
t59 = t111 * t180 - t196 * t89;
t54 = t121 * t211 + t123 * t188 + t125 * t189;
t53 = t120 * t211 + t122 * t188 + t124 * t189;
t52 = t121 * t209 + t123 * t186 + t125 * t187;
t51 = t120 * t209 + t122 * t186 + t124 * t187;
t50 = t113 * t211 - t115 * t182 + t117 * t183;
t49 = t112 * t211 - t114 * t182 + t116 * t183;
t48 = t113 * t209 - t115 * t180 + t117 * t181;
t47 = t112 * t209 - t114 * t180 + t116 * t181;
t46 = t217 * t119 + t211 * t292 + t101;
t45 = t209 * t150 + t217 * t298 + t130;
t44 = -t180 * t90 + t182 * t89;
t42 = t43 * t227;
t41 = t43 * t217;
t39 = t227 * t119 + (-t194 + t292) * t220 + t299;
t38 = t219 * t150 + (-t170 + t298) * t227 + t293;
t35 = t211 * t118 + t209 * t297 + t100;
t34 = t220 * t118 + (-t171 + t297) * t219 + t300;
t27 = t283 * t211 + t217 * t311 + t101;
t26 = t209 * t296 + t217 * t285 + t130;
t25 = t311 * t227 + (-t194 + t283) * t220 + t299;
t24 = t296 * t219 + (-t170 + t285) * t227 + t293;
t23 = t284 * t209 + t211 * t312 + t100;
t22 = t312 * t220 + (-t171 + t284) * t219 + t300;
t21 = t57 * t219 + t58 * t220 + t74;
t20 = t57 * t209 + t58 * t211 + t73;
t19 = t55 * t219 + t56 * t220 + t71;
t18 = t55 * t209 + t56 * t211 + t70;
t17 = t219 * t53 + t220 * t54 + t227 * t67;
t16 = t219 * t51 + t220 * t52 + t227 * t66;
t15 = t209 * t53 + t211 * t54 + t217 * t67;
t14 = t209 * t51 + t211 * t52 + t217 * t66;
t13 = t219 * t49 + t220 * t50 + t227 * t65;
t12 = t219 * t47 + t220 * t48 + t227 * t64;
t11 = t209 * t49 + t211 * t50 + t217 * t65;
t10 = t209 * t47 + t211 * t48 + t217 * t64;
t9 = t32 * t219 + t33 * t220 + t42;
t8 = t32 * t209 + t33 * t211 + t41;
t6 = t219 * t30 + t220 * t31 + t227 * t37;
t5 = t219 * t28 + t220 * t29 + t227 * t36;
t4 = t209 * t30 + t211 * t31 + t217 * t37;
t3 = t209 * t28 + t211 * t29 + t217 * t36;
t95 = [t282 + (Icges(3,5) * t244 + (Icges(3,1) * t307 + Icges(3,4) * t309) * t243) * t277 + m(7) * (t62 ^ 2 + t63 ^ 2) + m(6) * (t81 ^ 2 + t82 ^ 2) + m(5) * (t91 ^ 2 + t92 ^ 2) + m(4) * (t133 ^ 2 + t134 ^ 2) + m(3) * (t205 ^ 2 + t206 ^ 2) + m(2) * (t234 ^ 2 + t235 ^ 2) + t244 * (Icges(3,3) * t244 + (Icges(3,5) * t307 + Icges(3,6) * t309) * t243) + t243 * t309 * (Icges(3,6) * t244 + (Icges(3,4) * t307 + Icges(3,2) * t309) * t243) + Icges(2,3) + t76 + t72 + t43; 0.2e1 * ((t249 * t62 - t252 * t63) * t322 + (t249 * t81 - t252 * t82) * t323 + (t249 * t91 - t252 * t92) * t324 + (t133 * t249 - t134 * t252) * t325 + (t205 * t249 - t206 * t252) * t326) * t243; 0.2e1 * (t326 + t325 + t324 + t289) * (t244 ^ 2 + (t249 ^ 2 + t252 ^ 2) * t243 ^ 2); t42 + t71 + t74 + t97 + m(7) * (t24 * t62 + t25 * t63) + m(6) * (t38 * t81 + t39 * t82) + m(5) * (t68 * t91 + t69 * t92) + m(4) * (t104 * t133 + t105 * t134) + (t94 / 0.2e1 + t78 / 0.2e1 + t262) * t220 + (t93 / 0.2e1 + t77 / 0.2e1 + t263) * t219; m(4) * (t98 * t244 + (t104 * t249 - t105 * t252) * t243) + m(5) * (t61 * t244 + (t249 * t68 - t252 * t69) * t243) + m(6) * (t34 * t244 + (t249 * t38 - t252 * t39) * t243) + m(7) * (t22 * t244 + (t24 * t249 - t25 * t252) * t243); (t19 + t21 + t9 + t97) * t227 + (t13 + t17 + t6 + (t157 * t220 - t159 * t211 + t161 * t212) * t220 + (t78 + t94) * t227) * t220 + (t5 + t16 + t12 + (t156 * t219 - t158 * t209 + t160 * t210) * t219 + (t93 + t77) * t227 + (t156 * t220 + t157 * t219 - t158 * t211 - t159 * t209 + t160 * t212 + t161 * t210) * t220) * t219 + m(7) * (t22 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t34 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t61 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(4) * (t104 ^ 2 + t105 ^ 2 + t98 ^ 2); t41 + t70 + t73 + m(7) * (t26 * t62 + t27 * t63) + m(6) * (t45 * t81 + t46 * t82) + m(5) * (t79 * t91 + t80 * t92) + t262 * t211 + t263 * t209; m(5) * (t75 * t244 + (t249 * t79 - t252 * t80) * t243) + m(6) * (t35 * t244 + (t249 * t45 - t252 * t46) * t243) + m(7) * (t23 * t244 + (t249 * t26 - t252 * t27) * t243); (t8 / 0.2e1 + t18 / 0.2e1 + t20 / 0.2e1) * t227 + (t4 / 0.2e1 + t11 / 0.2e1 + t15 / 0.2e1) * t220 + (t3 / 0.2e1 + t14 / 0.2e1 + t10 / 0.2e1) * t219 + (t9 / 0.2e1 + t19 / 0.2e1 + t21 / 0.2e1) * t217 + (t6 / 0.2e1 + t13 / 0.2e1 + t17 / 0.2e1) * t211 + (t5 / 0.2e1 + t16 / 0.2e1 + t12 / 0.2e1) * t209 + m(7) * (t22 * t23 + t24 * t26 + t25 * t27) + m(6) * (t34 * t35 + t38 * t45 + t39 * t46) + m(5) * (t61 * t75 + t68 * t79 + t69 * t80); (t8 + t18 + t20) * t217 + (t4 + t15 + t11) * t211 + (t3 + t14 + t10) * t209 + m(7) * (t23 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t35 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t75 ^ 2 + t79 ^ 2 + t80 ^ 2); m(7) * (t209 * t63 + t211 * t62) + m(6) * (t209 * t82 + t211 * t81); (t217 * t244 + (-t209 * t252 + t211 * t249) * t243) * t327; m(7) * (t209 * t25 + t211 * t24 + t217 * t22) + m(6) * (t209 * t39 + t211 * t38 + t217 * t34); m(7) * (t209 * t27 + t211 * t26 + t217 * t23) + m(6) * (t209 * t46 + t211 * t45 + t217 * t35); (t209 ^ 2 + t211 ^ 2 + t217 ^ 2) * t327; t40 + m(7) * (t59 * t62 + t60 * t63) + t286 * t182 + t287 * t180; m(7) * (t44 * t244 + (t249 * t59 - t252 * t60) * t243); t9 * t316 + m(7) * (t22 * t44 + t24 * t59 + t25 * t60) + t6 * t317 + t5 * t318 + t227 * t319 + t219 * t321 + t220 * t320; m(7) * (t23 * t44 + t26 * t59 + t27 * t60) + t217 * t319 + t4 * t317 + t3 * t318 + t209 * t321 + t8 * t316 + t211 * t320; m(7) * (t209 * t60 + t211 * t59 + t217 * t44); t182 * t2 + t180 * t1 + t196 * t7 + m(7) * (t44 ^ 2 + t59 ^ 2 + t60 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t95(1) t95(2) t95(4) t95(7) t95(11) t95(16); t95(2) t95(3) t95(5) t95(8) t95(12) t95(17); t95(4) t95(5) t95(6) t95(9) t95(13) t95(18); t95(7) t95(8) t95(9) t95(10) t95(14) t95(19); t95(11) t95(12) t95(13) t95(14) t95(15) t95(20); t95(16) t95(17) t95(18) t95(19) t95(20) t95(21);];
Mq  = res;
