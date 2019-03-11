% Calculate joint inertia matrix for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:24
% EndTime: 2019-03-08 22:02:44
% DurationCPUTime: 9.51s
% Computational Cost: add. (55621->627), mult. (151133->909), div. (0->0), fcn. (200373->16), ass. (0->293)
t316 = m(5) + m(6) + m(7);
t299 = sin(pkin(13));
t300 = sin(pkin(7));
t257 = t300 * t299;
t301 = cos(pkin(13));
t258 = t301 * t300;
t305 = sin(qJ(3));
t307 = cos(qJ(3));
t224 = t257 * t307 + t258 * t305;
t240 = cos(pkin(7));
t249 = t299 * t307 + t301 * t305;
t226 = t249 * t240;
t237 = sin(pkin(12));
t241 = cos(pkin(6));
t244 = sin(qJ(2));
t239 = cos(pkin(12));
t246 = cos(qJ(2));
t293 = t246 * t239;
t230 = -t237 * t244 + t241 * t293;
t294 = t246 * t237;
t295 = t241 * t244;
t231 = t239 * t295 + t294;
t235 = -t299 * t305 + t301 * t307;
t238 = sin(pkin(6));
t297 = t238 * t239;
t179 = -t224 * t297 + t226 * t230 + t231 * t235;
t243 = sin(qJ(5));
t296 = t238 * t240;
t250 = t230 * t300 + t239 * t296;
t306 = cos(qJ(5));
t172 = t179 * t243 + t250 * t306;
t232 = -t239 * t244 - t241 * t294;
t233 = -t237 * t295 + t293;
t298 = t237 * t238;
t181 = t224 * t298 + t226 * t232 + t233 * t235;
t251 = t232 * t300 - t237 * t296;
t174 = t181 * t243 + t251 * t306;
t193 = t241 * t224 + (t226 * t246 + t235 * t244) * t238;
t229 = -t238 * t246 * t300 + t241 * t240;
t182 = t193 * t243 - t229 * t306;
t173 = t179 * t306 - t243 * t250;
t225 = t235 * t240;
t248 = -t257 * t305 + t258 * t307;
t247 = t238 * t248;
t178 = t225 * t230 - t231 * t249 - t239 * t247;
t242 = sin(qJ(6));
t245 = cos(qJ(6));
t129 = -t173 * t242 - t178 * t245;
t130 = t173 * t245 - t178 * t242;
t88 = Icges(7,5) * t130 + Icges(7,6) * t129 + Icges(7,3) * t172;
t90 = Icges(7,4) * t130 + Icges(7,2) * t129 + Icges(7,6) * t172;
t92 = Icges(7,1) * t130 + Icges(7,4) * t129 + Icges(7,5) * t172;
t28 = t129 * t90 + t130 * t92 + t172 * t88;
t175 = t181 * t306 - t243 * t251;
t180 = t232 * t225 - t233 * t249 + t237 * t247;
t131 = -t175 * t242 - t180 * t245;
t132 = t175 * t245 - t180 * t242;
t89 = Icges(7,5) * t132 + Icges(7,6) * t131 + Icges(7,3) * t174;
t91 = Icges(7,4) * t132 + Icges(7,2) * t131 + Icges(7,6) * t174;
t93 = Icges(7,1) * t132 + Icges(7,4) * t131 + Icges(7,5) * t174;
t29 = t129 * t91 + t130 * t93 + t172 * t89;
t183 = t193 * t306 + t229 * t243;
t192 = t241 * t248 + (t225 * t246 - t244 * t249) * t238;
t150 = -t183 * t242 - t192 * t245;
t151 = t183 * t245 - t192 * t242;
t102 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t182;
t103 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t182;
t104 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t182;
t41 = t102 * t172 + t103 * t129 + t104 * t130;
t1 = t172 * t28 + t174 * t29 + t182 * t41;
t313 = -t1 / 0.2e1;
t315 = pkin(9) * t250;
t314 = pkin(9) * t251;
t30 = t131 * t90 + t132 * t92 + t174 * t88;
t31 = t131 * t91 + t132 * t93 + t174 * t89;
t42 = t102 * t174 + t103 * t131 + t104 * t132;
t2 = t172 * t30 + t174 * t31 + t182 * t42;
t312 = t2 / 0.2e1;
t34 = t150 * t90 + t151 * t92 + t182 * t88;
t35 = t150 * t91 + t151 * t93 + t182 * t89;
t45 = t102 * t182 + t103 * t150 + t104 * t151;
t9 = t172 * t34 + t174 * t35 + t182 * t45;
t311 = t9 / 0.2e1;
t310 = t172 / 0.2e1;
t309 = t174 / 0.2e1;
t308 = t182 / 0.2e1;
t304 = pkin(3) * t307;
t94 = rSges(7,1) * t130 + rSges(7,2) * t129 + rSges(7,3) * t172;
t303 = pkin(5) * t173 + pkin(11) * t172 + t94;
t95 = rSges(7,1) * t132 + rSges(7,2) * t131 + rSges(7,3) * t174;
t302 = pkin(5) * t175 + pkin(11) * t174 + t95;
t105 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t182;
t292 = pkin(5) * t183 + pkin(11) * t182 + t105;
t139 = rSges(5,1) * t179 + rSges(5,2) * t178 - rSges(5,3) * t250;
t259 = t300 * t305;
t227 = pkin(3) * t259 + (pkin(9) + qJ(4)) * t240;
t267 = t300 * pkin(9);
t268 = t240 * t305;
t228 = pkin(3) * t268 - qJ(4) * t300 - t267;
t170 = -t227 * t297 + t230 * t228 + t231 * t304 + t315;
t291 = -t139 - t170;
t146 = pkin(4) * t179 - pkin(10) * t178;
t152 = t251 * t170;
t290 = -t146 * t251 - t152;
t147 = pkin(4) * t181 - pkin(10) * t180;
t171 = t227 * t298 + t232 * t228 + t233 * t304 + t314;
t157 = t229 * t171;
t289 = t229 * t147 + t157;
t288 = -t146 - t170;
t287 = -t147 - t171;
t166 = pkin(4) * t193 - pkin(10) * t192;
t188 = (-t240 * pkin(9) + t227) * t241 + ((t267 + t228) * t246 + t304 * t244) * t238;
t176 = t250 * t188;
t286 = -t166 * t250 - t176;
t156 = rSges(5,1) * t193 + rSges(5,2) * t192 + rSges(5,3) * t229;
t285 = -t156 - t188;
t284 = -t166 - t188;
t202 = t233 * pkin(2) - t314;
t200 = t241 * t202;
t283 = t241 * t171 + t200;
t201 = t231 * pkin(2) - t315;
t282 = t201 * t298 + t202 * t297;
t109 = Icges(6,5) * t173 - Icges(6,6) * t172 - Icges(6,3) * t178;
t111 = Icges(6,4) * t173 - Icges(6,2) * t172 - Icges(6,6) * t178;
t113 = Icges(6,1) * t173 - Icges(6,4) * t172 - Icges(6,5) * t178;
t54 = -t109 * t178 - t111 * t172 + t113 * t173;
t110 = Icges(6,5) * t175 - Icges(6,6) * t174 - Icges(6,3) * t180;
t112 = Icges(6,4) * t175 - Icges(6,2) * t174 - Icges(6,6) * t180;
t114 = Icges(6,1) * t175 - Icges(6,4) * t174 - Icges(6,5) * t180;
t55 = -t110 * t178 - t112 * t172 + t114 * t173;
t125 = Icges(6,5) * t183 - Icges(6,6) * t182 - Icges(6,3) * t192;
t126 = Icges(6,4) * t183 - Icges(6,2) * t182 - Icges(6,6) * t192;
t127 = Icges(6,1) * t183 - Icges(6,4) * t182 - Icges(6,5) * t192;
t67 = -t125 * t178 - t126 * t172 + t127 * t173;
t13 = -t178 * t54 - t180 * t55 - t192 * t67;
t3 = -t178 * t28 - t180 * t29 - t192 * t41;
t281 = t3 / 0.2e1 + t13 / 0.2e1;
t56 = -t109 * t180 - t111 * t174 + t113 * t175;
t57 = -t110 * t180 - t112 * t174 + t114 * t175;
t68 = -t125 * t180 - t126 * t174 + t127 * t175;
t14 = -t178 * t56 - t180 * t57 - t192 * t68;
t4 = -t178 * t30 - t180 * t31 - t192 * t42;
t280 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t229 * t67 - t250 * t54 - t251 * t55;
t5 = t229 * t41 - t250 * t28 - t251 * t29;
t279 = -t5 / 0.2e1 - t15 / 0.2e1;
t16 = t229 * t68 - t250 * t56 - t251 * t57;
t6 = t229 * t42 - t250 * t30 - t251 * t31;
t278 = -t6 / 0.2e1 - t16 / 0.2e1;
t17 = t241 * t67 + (t237 * t55 - t239 * t54) * t238;
t7 = t241 * t41 + (t237 * t29 - t239 * t28) * t238;
t277 = -t7 / 0.2e1 - t17 / 0.2e1;
t18 = t241 * t68 + (t237 * t57 - t239 * t56) * t238;
t8 = t241 * t42 + (t237 * t31 - t239 * t30) * t238;
t276 = -t8 / 0.2e1 - t18 / 0.2e1;
t10 = -t178 * t34 - t180 * t35 - t192 * t45;
t59 = -t109 * t192 - t111 * t182 + t113 * t183;
t60 = -t110 * t192 - t112 * t182 + t114 * t183;
t69 = -t125 * t192 - t126 * t182 + t127 * t183;
t19 = -t178 * t59 - t180 * t60 - t192 * t69;
t275 = t10 / 0.2e1 + t19 / 0.2e1;
t11 = t229 * t45 - t250 * t34 - t251 * t35;
t20 = t229 * t69 - t250 * t59 - t251 * t60;
t274 = -t11 / 0.2e1 - t20 / 0.2e1;
t12 = t241 * t45 + (t237 * t35 - t239 * t34) * t238;
t21 = t241 * t69 + (t237 * t60 - t239 * t59) * t238;
t273 = -t12 / 0.2e1 - t21 / 0.2e1;
t115 = rSges(6,1) * t173 - rSges(6,2) * t172 - rSges(6,3) * t178;
t272 = -t115 + t288;
t128 = rSges(6,1) * t183 - rSges(6,2) * t182 - rSges(6,3) * t192;
t271 = -t128 + t284;
t270 = t241 * t147 + t283;
t269 = t240 * t307;
t260 = t307 * t300;
t213 = t241 * t260 + (-t244 * t305 + t246 * t269) * t238;
t214 = t241 * t259 + (t244 * t307 + t246 * t268) * t238;
t187 = rSges(4,1) * t214 + rSges(4,2) * t213 + rSges(4,3) * t229;
t218 = t238 * t244 * pkin(2) + pkin(9) * t229;
t266 = (-t187 - t218) * t238;
t264 = t288 - t303;
t263 = t284 - t292;
t262 = t170 * t298 + t171 * t297 + t282;
t261 = (-t218 + t285) * t238;
t256 = (-t218 + t271) * t238;
t255 = t238 * t260;
t254 = t238 * t259;
t253 = t146 * t298 + t147 * t297 + t262;
t252 = (-t218 + t263) * t238;
t222 = t241 * rSges(3,3) + (rSges(3,1) * t244 + rSges(3,2) * t246) * t238;
t221 = Icges(3,5) * t241 + (Icges(3,1) * t244 + Icges(3,4) * t246) * t238;
t220 = Icges(3,6) * t241 + (Icges(3,4) * t244 + Icges(3,2) * t246) * t238;
t219 = Icges(3,3) * t241 + (Icges(3,5) * t244 + Icges(3,6) * t246) * t238;
t210 = rSges(3,1) * t233 + rSges(3,2) * t232 + rSges(3,3) * t298;
t209 = rSges(3,1) * t231 + rSges(3,2) * t230 - rSges(3,3) * t297;
t208 = Icges(3,1) * t233 + Icges(3,4) * t232 + Icges(3,5) * t298;
t207 = Icges(3,1) * t231 + Icges(3,4) * t230 - Icges(3,5) * t297;
t206 = Icges(3,4) * t233 + Icges(3,2) * t232 + Icges(3,6) * t298;
t205 = Icges(3,4) * t231 + Icges(3,2) * t230 - Icges(3,6) * t297;
t204 = Icges(3,5) * t233 + Icges(3,6) * t232 + Icges(3,3) * t298;
t203 = Icges(3,5) * t231 + Icges(3,6) * t230 - Icges(3,3) * t297;
t199 = t232 * t268 + t233 * t307 + t237 * t254;
t198 = t232 * t269 - t233 * t305 + t237 * t255;
t197 = t230 * t268 + t231 * t307 - t239 * t254;
t196 = t230 * t269 - t231 * t305 - t239 * t255;
t191 = -t209 * t241 - t222 * t297;
t190 = t210 * t241 - t222 * t298;
t186 = Icges(4,1) * t214 + Icges(4,4) * t213 + Icges(4,5) * t229;
t185 = Icges(4,4) * t214 + Icges(4,2) * t213 + Icges(4,6) * t229;
t184 = Icges(4,5) * t214 + Icges(4,6) * t213 + Icges(4,3) * t229;
t177 = (t209 * t237 + t210 * t239) * t238;
t165 = rSges(4,1) * t199 + rSges(4,2) * t198 - rSges(4,3) * t251;
t164 = rSges(4,1) * t197 + rSges(4,2) * t196 - rSges(4,3) * t250;
t163 = Icges(4,1) * t199 + Icges(4,4) * t198 - Icges(4,5) * t251;
t162 = Icges(4,1) * t197 + Icges(4,4) * t196 - Icges(4,5) * t250;
t161 = Icges(4,4) * t199 + Icges(4,2) * t198 - Icges(4,6) * t251;
t160 = Icges(4,4) * t197 + Icges(4,2) * t196 - Icges(4,6) * t250;
t159 = Icges(4,5) * t199 + Icges(4,6) * t198 - Icges(4,3) * t251;
t158 = Icges(4,5) * t197 + Icges(4,6) * t196 - Icges(4,3) * t250;
t155 = Icges(5,1) * t193 + Icges(5,4) * t192 + Icges(5,5) * t229;
t154 = Icges(5,4) * t193 + Icges(5,2) * t192 + Icges(5,6) * t229;
t153 = Icges(5,5) * t193 + Icges(5,6) * t192 + Icges(5,3) * t229;
t140 = rSges(5,1) * t181 + rSges(5,2) * t180 - rSges(5,3) * t251;
t138 = Icges(5,1) * t181 + Icges(5,4) * t180 - Icges(5,5) * t251;
t137 = Icges(5,1) * t179 + Icges(5,4) * t178 - Icges(5,5) * t250;
t136 = Icges(5,4) * t181 + Icges(5,2) * t180 - Icges(5,6) * t251;
t135 = Icges(5,4) * t179 + Icges(5,2) * t178 - Icges(5,6) * t250;
t134 = Icges(5,5) * t181 + Icges(5,6) * t180 - Icges(5,3) * t251;
t133 = Icges(5,5) * t179 + Icges(5,6) * t178 - Icges(5,3) * t250;
t122 = t165 * t229 + t187 * t251;
t121 = -t164 * t229 - t187 * t250;
t120 = (-t164 - t201) * t241 + t239 * t266;
t119 = t165 * t241 + t237 * t266 + t200;
t118 = t184 * t229 + t185 * t213 + t186 * t214;
t117 = -t164 * t251 + t165 * t250;
t116 = rSges(6,1) * t175 - rSges(6,2) * t174 - rSges(6,3) * t180;
t108 = -t184 * t251 + t185 * t198 + t186 * t199;
t107 = -t184 * t250 + t185 * t196 + t186 * t197;
t106 = (t164 * t237 + t165 * t239) * t238 + t282;
t101 = t159 * t229 + t161 * t213 + t163 * t214;
t100 = t158 * t229 + t160 * t213 + t162 * t214;
t99 = -t159 * t251 + t161 * t198 + t163 * t199;
t98 = -t158 * t251 + t160 * t198 + t162 * t199;
t97 = -t159 * t250 + t161 * t196 + t163 * t197;
t96 = -t158 * t250 + t160 * t196 + t162 * t197;
t87 = t153 * t229 + t154 * t192 + t155 * t193;
t86 = t140 * t229 - t251 * t285 + t157;
t85 = -t156 * t250 + t229 * t291 - t176;
t84 = -t153 * t251 + t154 * t180 + t155 * t181;
t83 = -t153 * t250 + t154 * t178 + t155 * t179;
t82 = (-t201 + t291) * t241 + t239 * t261;
t81 = t140 * t241 + t237 * t261 + t283;
t80 = -t116 * t192 + t128 * t180;
t79 = t115 * t192 - t128 * t178;
t78 = t134 * t229 + t136 * t192 + t138 * t193;
t77 = t133 * t229 + t135 * t192 + t137 * t193;
t76 = -t139 * t251 - t152 - (-t140 - t171) * t250;
t75 = (t139 * t237 + t140 * t239) * t238 + t262;
t74 = -t134 * t251 + t136 * t180 + t138 * t181;
t73 = -t133 * t251 + t135 * t180 + t137 * t181;
t72 = -t134 * t250 + t136 * t178 + t138 * t179;
t71 = -t133 * t250 + t135 * t178 + t137 * t179;
t70 = -t115 * t180 + t116 * t178;
t66 = -t105 * t174 + t182 * t95;
t65 = t105 * t172 - t182 * t94;
t64 = (-t201 + t272) * t241 + t239 * t256;
t63 = t116 * t241 + t237 * t256 + t270;
t62 = t116 * t229 - t251 * t271 + t289;
t61 = -t128 * t250 + t229 * t272 + t286;
t58 = t118 * t241 + (-t100 * t239 + t101 * t237) * t238;
t53 = -t100 * t250 - t101 * t251 + t118 * t229;
t52 = -t172 * t95 + t174 * t94;
t51 = (t115 * t237 + t116 * t239) * t238 + t253;
t50 = -t115 * t251 - (-t116 + t287) * t250 + t290;
t49 = t108 * t241 + (t237 * t99 - t239 * t98) * t238;
t48 = t107 * t241 + (t237 * t97 - t239 * t96) * t238;
t47 = t108 * t229 - t250 * t98 - t251 * t99;
t46 = t107 * t229 - t250 * t96 - t251 * t97;
t44 = t180 * t292 - t192 * t302;
t43 = -t178 * t292 + t192 * t303;
t40 = t178 * t302 - t180 * t303;
t39 = (-t201 + t264) * t241 + t239 * t252;
t38 = t237 * t252 + t241 * t302 + t270;
t37 = t229 * t302 - t251 * t263 + t289;
t36 = t229 * t264 - t250 * t292 + t286;
t33 = t241 * t87 + (t237 * t78 - t239 * t77) * t238;
t32 = t229 * t87 - t250 * t77 - t251 * t78;
t27 = t241 * t84 + (t237 * t74 - t239 * t73) * t238;
t26 = t241 * t83 + (t237 * t72 - t239 * t71) * t238;
t25 = (t237 * t303 + t239 * t302) * t238 + t253;
t24 = t229 * t84 - t250 * t73 - t251 * t74;
t23 = t229 * t83 - t250 * t71 - t251 * t72;
t22 = -t303 * t251 - (t287 - t302) * t250 + t290;
t123 = [m(2) + m(3) + m(4) + t316; m(3) * t177 + m(4) * t106 + m(5) * t75 + m(6) * t51 + m(7) * t25; m(7) * (t25 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(6) * (t51 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t75 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t106 ^ 2 + t119 ^ 2 + t120 ^ 2) + m(3) * (t177 ^ 2 + t190 ^ 2 + t191 ^ 2) + (t8 + t18 + t49 + t27 + (t204 * t298 + t206 * t232 + t208 * t233) * t298) * t298 + (-t7 - t17 - t48 - t26 + (-t203 * t297 + t205 * t230 + t207 * t231) * t297 + (-t203 * t298 + t204 * t297 - t205 * t232 - t206 * t230 - t207 * t233 - t208 * t231) * t298) * t297 + ((t219 * t298 + t220 * t232 + t221 * t233) * t298 - (-t219 * t297 + t220 * t230 + t221 * t231) * t297 + t12 + t21 + t33 + t58 + ((t206 * t246 + t208 * t244) * t237 - (t205 * t246 + t207 * t244) * t239) * t238 ^ 2 + ((-t203 * t239 + t204 * t237 + t220 * t246 + t221 * t244) * t238 + t241 * t219) * t241) * t241; m(4) * t117 + m(5) * t76 + m(6) * t50 + m(7) * t22; (t32 / 0.2e1 + t53 / 0.2e1 - t274) * t241 + (t33 / 0.2e1 + t58 / 0.2e1 - t273) * t229 - (t27 / 0.2e1 + t49 / 0.2e1 - t276) * t251 - (t26 / 0.2e1 + t48 / 0.2e1 - t277) * t250 + m(7) * (t22 * t25 + t36 * t39 + t37 * t38) + m(6) * (t50 * t51 + t61 * t64 + t62 * t63) + m(5) * (t75 * t76 + t81 * t86 + t82 * t85) + m(4) * (t106 * t117 + t119 * t122 + t120 * t121) + ((-t23 / 0.2e1 - t46 / 0.2e1 + t279) * t239 + (t24 / 0.2e1 + t47 / 0.2e1 - t278) * t237) * t238; (t11 + t20 + t53 + t32) * t229 - (t6 + t16 + t47 + t24) * t251 - (t5 + t15 + t46 + t23) * t250 + m(7) * (t22 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t50 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t76 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(4) * (t117 ^ 2 + t121 ^ 2 + t122 ^ 2); t229 * t316; m(7) * (t229 * t25 - t250 * t38 - t251 * t39) + m(6) * (t229 * t51 - t250 * t63 - t251 * t64) + m(5) * (t229 * t75 - t250 * t81 - t251 * t82); m(7) * (t22 * t229 - t250 * t37 - t251 * t36) + m(6) * (t229 * t50 - t250 * t62 - t251 * t61) + m(5) * (t229 * t76 - t250 * t86 - t251 * t85); (t229 ^ 2 + t250 ^ 2 + t251 ^ 2) * t316; m(6) * t70 + m(7) * t40; t275 * t241 + t273 * t192 + t276 * t180 + t277 * t178 + m(7) * (t25 * t40 + t38 * t44 + t39 * t43) + m(6) * (t51 * t70 + t63 * t80 + t64 * t79) + (t237 * t280 - t239 * t281) * t238; t275 * t229 - t280 * t251 - t281 * t250 + t274 * t192 + t278 * t180 + t279 * t178 + m(7) * (t22 * t40 + t36 * t43 + t37 * t44) + m(6) * (t50 * t70 + t61 * t79 + t62 * t80); m(6) * (t229 * t70 - t250 * t80 - t251 * t79) + m(7) * (t229 * t40 - t250 * t44 - t251 * t43); (-t10 - t19) * t192 + (-t4 - t14) * t180 + (-t3 - t13) * t178 + m(7) * (t40 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t70 ^ 2 + t79 ^ 2 + t80 ^ 2); m(7) * t52; t7 * t310 + t12 * t308 + m(7) * (t25 * t52 + t38 * t66 + t39 * t65) + t241 * t311 + t8 * t309 + (t237 * t312 + t239 * t313) * t238; t6 * t309 + t250 * t313 + t229 * t311 + t5 * t310 - t251 * t312 + t11 * t308 + m(7) * (t22 * t52 + t36 * t65 + t37 * t66); m(7) * (t229 * t52 - t250 * t66 - t251 * t65); m(7) * (t40 * t52 + t43 * t65 + t44 * t66) + t4 * t309 - t180 * t2 / 0.2e1 + t3 * t310 - t192 * t9 / 0.2e1 + t10 * t308 + t178 * t313; t174 * t2 + t172 * t1 + t182 * t9 + m(7) * (t52 ^ 2 + t65 ^ 2 + t66 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t123(1) t123(2) t123(4) t123(7) t123(11) t123(16); t123(2) t123(3) t123(5) t123(8) t123(12) t123(17); t123(4) t123(5) t123(6) t123(9) t123(13) t123(18); t123(7) t123(8) t123(9) t123(10) t123(14) t123(19); t123(11) t123(12) t123(13) t123(14) t123(15) t123(20); t123(16) t123(17) t123(18) t123(19) t123(20) t123(21);];
Mq  = res;
