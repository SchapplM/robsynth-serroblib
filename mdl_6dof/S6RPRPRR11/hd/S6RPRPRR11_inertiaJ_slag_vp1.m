% Calculate joint inertia matrix for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_inertiaJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:26
% EndTime: 2019-03-09 04:11:39
% DurationCPUTime: 5.48s
% Computational Cost: add. (31396->532), mult. (74065->743), div. (0->0), fcn. (97242->16), ass. (0->250)
t232 = sin(qJ(1));
t228 = cos(pkin(6));
t234 = cos(qJ(1));
t282 = sin(pkin(12));
t256 = t234 * t282;
t284 = cos(pkin(12));
t259 = t232 * t284;
t242 = t228 * t259 + t256;
t226 = sin(pkin(6));
t285 = cos(pkin(7));
t262 = t226 * t285;
t283 = sin(pkin(7));
t202 = t232 * t262 + t242 * t283;
t257 = t234 * t284;
t258 = t232 * t282;
t243 = -t228 * t257 + t258;
t201 = -t234 * t262 + t243 * t283;
t231 = sin(qJ(3));
t250 = t285 * t284;
t296 = cos(qJ(3));
t200 = t228 * t283 * t231 + (t231 * t250 + t282 * t296) * t226;
t261 = t226 * t283;
t209 = t228 * t285 - t284 * t261;
t225 = sin(pkin(13));
t227 = cos(pkin(13));
t182 = -t200 * t225 + t209 * t227;
t278 = t209 * t225;
t183 = t200 * t227 + t278;
t251 = t296 * t283;
t260 = t226 * t282;
t199 = -t226 * t296 * t250 - t228 * t251 + t231 * t260;
t134 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t199;
t135 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t199;
t136 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t199;
t172 = Icges(4,5) * t200 - Icges(4,6) * t199 + Icges(4,3) * t209;
t173 = Icges(4,4) * t200 - Icges(4,2) * t199 + Icges(4,6) * t209;
t174 = Icges(4,1) * t200 - Icges(4,4) * t199 + Icges(4,5) * t209;
t311 = t182 * t135 + t183 * t136 + t209 * t172 + t200 * t174 + (t134 - t173) * t199;
t303 = m(7) / 0.2e1;
t304 = m(6) / 0.2e1;
t305 = m(5) / 0.2e1;
t254 = t305 + t304 + t303;
t310 = 0.2e1 * t254;
t210 = t228 * t256 + t259;
t241 = t243 * t285;
t192 = t210 * t296 + (-t234 * t261 - t241) * t231;
t168 = -t192 * t225 + t201 * t227;
t280 = t201 * t225;
t169 = t192 * t227 + t280;
t247 = t226 * t251;
t191 = t210 * t231 + t234 * t247 + t296 * t241;
t103 = Icges(5,5) * t169 + Icges(5,6) * t168 + Icges(5,3) * t191;
t140 = Icges(4,4) * t192 - Icges(4,2) * t191 + Icges(4,6) * t201;
t309 = t103 - t140;
t211 = -t228 * t258 + t257;
t239 = t242 * t285;
t194 = t211 * t296 + (t232 * t261 - t239) * t231;
t170 = -t194 * t225 + t202 * t227;
t279 = t202 * t225;
t171 = t194 * t227 + t279;
t193 = t211 * t231 - t232 * t247 + t296 * t239;
t104 = Icges(5,5) * t171 + Icges(5,6) * t170 + Icges(5,3) * t193;
t141 = Icges(4,4) * t194 - Icges(4,2) * t193 + Icges(4,6) * t202;
t308 = t104 - t141;
t307 = m(3) / 0.2e1;
t306 = m(4) / 0.2e1;
t270 = pkin(13) + qJ(5);
t222 = sin(t270);
t255 = cos(t270);
t162 = t192 * t222 - t201 * t255;
t164 = t194 * t222 - t202 * t255;
t178 = t200 * t222 - t209 * t255;
t163 = t192 * t255 + t201 * t222;
t230 = sin(qJ(6));
t233 = cos(qJ(6));
t128 = -t163 * t230 + t191 * t233;
t129 = t163 * t233 + t191 * t230;
t68 = Icges(7,5) * t129 + Icges(7,6) * t128 + Icges(7,3) * t162;
t70 = Icges(7,4) * t129 + Icges(7,2) * t128 + Icges(7,6) * t162;
t72 = Icges(7,1) * t129 + Icges(7,4) * t128 + Icges(7,5) * t162;
t19 = t128 * t70 + t129 * t72 + t162 * t68;
t165 = t194 * t255 + t202 * t222;
t130 = -t165 * t230 + t193 * t233;
t131 = t165 * t233 + t193 * t230;
t69 = Icges(7,5) * t131 + Icges(7,6) * t130 + Icges(7,3) * t164;
t71 = Icges(7,4) * t131 + Icges(7,2) * t130 + Icges(7,6) * t164;
t73 = Icges(7,1) * t131 + Icges(7,4) * t130 + Icges(7,5) * t164;
t20 = t128 * t71 + t129 * t73 + t162 * t69;
t179 = t200 * t255 + t209 * t222;
t155 = -t179 * t230 + t199 * t233;
t156 = t179 * t233 + t199 * t230;
t91 = Icges(7,5) * t156 + Icges(7,6) * t155 + Icges(7,3) * t178;
t92 = Icges(7,4) * t156 + Icges(7,2) * t155 + Icges(7,6) * t178;
t93 = Icges(7,1) * t156 + Icges(7,4) * t155 + Icges(7,5) * t178;
t27 = t128 * t92 + t129 * t93 + t162 * t91;
t1 = t162 * t19 + t164 * t20 + t178 * t27;
t302 = t1 / 0.2e1;
t21 = t130 * t70 + t131 * t72 + t164 * t68;
t22 = t130 * t71 + t131 * t73 + t164 * t69;
t28 = t130 * t92 + t131 * t93 + t164 * t91;
t2 = t162 * t21 + t164 * t22 + t178 * t28;
t301 = t2 / 0.2e1;
t23 = t155 * t70 + t156 * t72 + t178 * t68;
t24 = t155 * t71 + t156 * t73 + t178 * t69;
t36 = t155 * t92 + t156 * t93 + t178 * t91;
t33 = t36 * t178;
t7 = t23 * t162 + t24 * t164 + t33;
t300 = t7 / 0.2e1;
t299 = t162 / 0.2e1;
t298 = t164 / 0.2e1;
t297 = t178 / 0.2e1;
t295 = pkin(5) * t163;
t221 = pkin(4) * t227 + pkin(3);
t294 = -pkin(3) + t221;
t293 = t311 * t209;
t248 = -rSges(7,1) * t129 - rSges(7,2) * t128;
t74 = rSges(7,3) * t162 - t248;
t292 = pkin(11) * t162 + t295 + t74;
t75 = t131 * rSges(7,1) + t130 * rSges(7,2) + t164 * rSges(7,3);
t291 = t165 * pkin(5) + t164 * pkin(11) + t75;
t181 = t191 * qJ(4);
t152 = pkin(3) * t192 + t181;
t146 = t202 * t152;
t269 = pkin(4) * t280;
t229 = -pkin(10) - qJ(4);
t281 = t191 * t229;
t89 = t294 * t192 - t181 + t269 - t281;
t290 = t202 * t89 + t146;
t153 = t194 * pkin(3) + t193 * qJ(4);
t148 = t209 * t153;
t264 = pkin(4) * t279 - t193 * t229 + t194 * t221;
t90 = -t153 + t264;
t289 = t209 * t90 + t148;
t94 = rSges(7,1) * t156 + rSges(7,2) * t155 + rSges(7,3) * t178;
t288 = pkin(5) * t179 + pkin(11) * t178 + t94;
t287 = -t152 - t89;
t286 = -t153 - t90;
t277 = t226 * t232;
t276 = t226 * t234;
t109 = rSges(5,1) * t169 + rSges(5,2) * t168 + rSges(5,3) * t191;
t275 = -t109 - t152;
t110 = t171 * rSges(5,1) + t170 * rSges(5,2) + t193 * rSges(5,3);
t274 = -t110 - t153;
t122 = pkin(4) * t278 + t294 * t200 + (-qJ(4) - t229) * t199;
t176 = pkin(3) * t200 + qJ(4) * t199;
t161 = t201 * t176;
t273 = t201 * t122 + t161;
t272 = -t122 - t176;
t271 = t234 * pkin(1) + qJ(2) * t277;
t268 = t27 / 0.2e1 + t23 / 0.2e1;
t267 = t28 / 0.2e1 + t24 / 0.2e1;
t125 = Icges(6,5) * t179 - Icges(6,6) * t178 + Icges(6,3) * t199;
t126 = Icges(6,4) * t179 - Icges(6,2) * t178 + Icges(6,6) * t199;
t127 = Icges(6,1) * t179 - Icges(6,4) * t178 + Icges(6,5) * t199;
t59 = t199 * t125 - t178 * t126 + t179 * t127;
t102 = t165 * rSges(6,1) - t164 * rSges(6,2) + t193 * rSges(6,3);
t145 = t194 * rSges(4,1) - t193 * rSges(4,2) + t202 * rSges(4,3);
t263 = -t232 * pkin(1) + qJ(2) * t276;
t249 = -rSges(6,1) * t163 + rSges(6,2) * t162;
t95 = Icges(6,5) * t163 - Icges(6,6) * t162 + Icges(6,3) * t191;
t97 = Icges(6,4) * t163 - Icges(6,2) * t162 + Icges(6,6) * t191;
t99 = Icges(6,1) * t163 - Icges(6,4) * t162 + Icges(6,5) * t191;
t42 = -t178 * t97 + t179 * t99 + t199 * t95;
t51 = t125 * t191 - t126 * t162 + t127 * t163;
t246 = t51 / 0.2e1 + t42 / 0.2e1 + t268;
t100 = Icges(6,1) * t165 - Icges(6,4) * t164 + Icges(6,5) * t193;
t96 = Icges(6,5) * t165 - Icges(6,6) * t164 + Icges(6,3) * t193;
t98 = Icges(6,4) * t165 - Icges(6,2) * t164 + Icges(6,6) * t193;
t43 = t100 * t179 - t178 * t98 + t199 * t96;
t52 = t125 * t193 - t126 * t164 + t127 * t165;
t245 = t52 / 0.2e1 + t43 / 0.2e1 + t267;
t144 = rSges(4,1) * t192 - rSges(4,2) * t191 + rSges(4,3) * t201;
t244 = -pkin(2) * t210 - t201 * pkin(9) + t263;
t237 = -t192 * t221 + t244 - t269;
t236 = t211 * pkin(2) + t202 * pkin(9) + t271;
t235 = t236 + t264;
t217 = rSges(2,1) * t234 - t232 * rSges(2,2);
t216 = -t232 * rSges(2,1) - rSges(2,2) * t234;
t190 = t211 * rSges(3,1) - t242 * rSges(3,2) + rSges(3,3) * t277 + t271;
t189 = -t210 * rSges(3,1) + t243 * rSges(3,2) + rSges(3,3) * t276 + t263;
t175 = rSges(4,1) * t200 - rSges(4,2) * t199 + rSges(4,3) * t209;
t143 = Icges(4,1) * t194 - Icges(4,4) * t193 + Icges(4,5) * t202;
t142 = Icges(4,1) * t192 - Icges(4,4) * t191 + Icges(4,5) * t201;
t139 = Icges(4,5) * t194 - Icges(4,6) * t193 + Icges(4,3) * t202;
t138 = Icges(4,5) * t192 - Icges(4,6) * t191 + Icges(4,3) * t201;
t137 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t199;
t132 = rSges(6,1) * t179 - rSges(6,2) * t178 + rSges(6,3) * t199;
t116 = t236 + t145;
t115 = -t144 + t244;
t108 = Icges(5,1) * t171 + Icges(5,4) * t170 + Icges(5,5) * t193;
t107 = Icges(5,1) * t169 + Icges(5,4) * t168 + Icges(5,5) * t191;
t106 = Icges(5,4) * t171 + Icges(5,2) * t170 + Icges(5,6) * t193;
t105 = Icges(5,4) * t169 + Icges(5,2) * t168 + Icges(5,6) * t191;
t101 = rSges(6,3) * t191 - t249;
t88 = t145 * t209 - t175 * t202;
t87 = -t144 * t209 + t175 * t201;
t83 = t144 * t202 - t145 * t201;
t79 = t172 * t202 - t173 * t193 + t174 * t194;
t78 = t172 * t201 - t173 * t191 + t174 * t192;
t77 = t236 - t274;
t76 = t244 + t275;
t67 = t235 + t102;
t66 = (-rSges(6,3) + t229) * t191 + t237 + t249;
t65 = t102 * t199 - t132 * t193;
t64 = -t101 * t199 + t132 * t191;
t63 = t139 * t209 - t141 * t199 + t143 * t200;
t62 = t138 * t209 - t140 * t199 + t142 * t200;
t60 = t101 * t193 - t102 * t191;
t58 = t59 * t209;
t57 = t59 * t199;
t56 = t110 * t209 + t148 + (-t137 - t176) * t202;
t55 = t137 * t201 + t275 * t209 + t161;
t54 = t134 * t193 + t135 * t170 + t136 * t171;
t53 = t134 * t191 + t135 * t168 + t136 * t169;
t50 = t235 + t291;
t49 = -t295 + t281 + (-rSges(7,3) - pkin(11)) * t162 + t237 + t248;
t48 = t109 * t202 + t274 * t201 + t146;
t47 = -t164 * t94 + t178 * t75;
t46 = t162 * t94 - t178 * t74;
t45 = t104 * t199 + t106 * t182 + t108 * t183;
t44 = t103 * t199 + t105 * t182 + t107 * t183;
t41 = t100 * t165 - t164 * t98 + t193 * t96;
t40 = -t164 * t97 + t165 * t99 + t193 * t95;
t39 = t100 * t163 - t162 * t98 + t191 * t96;
t38 = -t162 * t97 + t163 * t99 + t191 * t95;
t37 = -t162 * t75 + t164 * t74;
t35 = t36 * t209;
t34 = t36 * t199;
t32 = t102 * t209 + (-t132 + t272) * t202 + t289;
t31 = t132 * t201 + (-t101 + t287) * t209 + t273;
t30 = -t288 * t193 + t291 * t199;
t29 = t288 * t191 - t292 * t199;
t26 = t101 * t202 + (-t102 + t286) * t201 + t290;
t25 = -t291 * t191 + t292 * t193;
t18 = t291 * t209 + (t272 - t288) * t202 + t289;
t17 = t288 * t201 + (t287 - t292) * t209 + t273;
t16 = t292 * t202 + (t286 - t291) * t201 + t290;
t15 = t42 * t201 + t43 * t202 + t58;
t14 = t42 * t191 + t43 * t193 + t57;
t13 = t201 * t40 + t202 * t41 + t209 * t52;
t12 = t201 * t38 + t202 * t39 + t209 * t51;
t11 = t191 * t40 + t193 * t41 + t199 * t52;
t10 = t191 * t38 + t193 * t39 + t199 * t51;
t9 = t23 * t201 + t24 * t202 + t35;
t8 = t23 * t191 + t24 * t193 + t34;
t6 = t201 * t21 + t202 * t22 + t209 * t28;
t5 = t19 * t201 + t20 * t202 + t209 * t27;
t4 = t191 * t21 + t193 * t22 + t199 * t28;
t3 = t19 * t191 + t193 * t20 + t199 * t27;
t61 = [t226 * t284 * (Icges(3,6) * t228 + (t282 * Icges(3,4) + t284 * Icges(3,2)) * t226) + (Icges(3,5) * t228 + (t282 * Icges(3,1) + t284 * Icges(3,4)) * t226) * t260 + Icges(2,3) + m(7) * (t49 ^ 2 + t50 ^ 2) + m(6) * (t66 ^ 2 + t67 ^ 2) + m(5) * (t76 ^ 2 + t77 ^ 2) + m(4) * (t115 ^ 2 + t116 ^ 2) + m(3) * (t189 ^ 2 + t190 ^ 2) + m(2) * (t216 ^ 2 + t217 ^ 2) + t228 * (Icges(3,3) * t228 + (t282 * Icges(3,5) + t284 * Icges(3,6)) * t226) + t59 + t36 + t311; 0.2e1 * ((t232 * t49 - t234 * t50) * t303 + (t232 * t66 - t234 * t67) * t304 + (t232 * t76 - t234 * t77) * t305 + (t115 * t232 - t116 * t234) * t306 + (t189 * t232 - t190 * t234) * t307) * t226; 0.2e1 * (t306 + t307 + t254) * (t228 ^ 2 + (t232 ^ 2 + t234 ^ 2) * t226 ^ 2); t35 + t58 + m(7) * (t17 * t49 + t18 * t50) + m(6) * (t31 * t66 + t32 * t67) + m(5) * (t55 * t76 + t56 * t77) + m(4) * (t115 * t87 + t116 * t88) + (t79 / 0.2e1 + t54 / 0.2e1 + t63 / 0.2e1 + t45 / 0.2e1 + t245) * t202 + (t62 / 0.2e1 + t44 / 0.2e1 + t78 / 0.2e1 + t53 / 0.2e1 + t246) * t201 + t293; m(4) * (t83 * t228 + (t232 * t87 - t234 * t88) * t226) + m(5) * (t48 * t228 + (t232 * t55 - t234 * t56) * t226) + m(6) * (t26 * t228 + (t232 * t31 - t234 * t32) * t226) + m(7) * (t16 * t228 + (t17 * t232 - t18 * t234) * t226); (t15 + t9 + t293) * t209 + (t13 + t6 + (t106 * t170 + t108 * t171 + t139 * t202 + t143 * t194 + t308 * t193) * t202 + (t45 + t54 + t63 + t79) * t209) * t202 + (t5 + t12 + (t105 * t168 + t107 * t169 + t138 * t201 + t142 * t192 + t309 * t191) * t201 + (t78 + t62 + t44 + t53) * t209 + (t105 * t170 + t106 * t168 + t107 * t171 + t108 * t169 + t138 * t202 + t139 * t201 + t142 * t194 + t143 * t192 + t308 * t191 + t309 * t193) * t202) * t201 + m(7) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) + m(6) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t48 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(4) * (t83 ^ 2 + t87 ^ 2 + t88 ^ 2); m(7) * (t191 * t50 + t193 * t49) + m(6) * (t191 * t67 + t193 * t66) + m(5) * (t191 * t77 + t193 * t76); (t199 * t228 + (-t191 * t234 + t193 * t232) * t226) * t310; m(7) * (t16 * t199 + t17 * t193 + t18 * t191) + m(6) * (t191 * t32 + t193 * t31 + t199 * t26) + m(5) * (t191 * t56 + t193 * t55 + t199 * t48); (t191 ^ 2 + t193 ^ 2 + t199 ^ 2) * t310; t34 + t57 + m(7) * (t29 * t49 + t30 * t50) + m(6) * (t64 * t66 + t65 * t67) + t245 * t193 + t246 * t191; m(6) * (t60 * t228 + (t232 * t64 - t234 * t65) * t226) + m(7) * (t25 * t228 + (t232 * t29 - t234 * t30) * t226); (t8 / 0.2e1 + t14 / 0.2e1) * t209 + (t4 / 0.2e1 + t11 / 0.2e1) * t202 + (t3 / 0.2e1 + t10 / 0.2e1) * t201 + (t9 / 0.2e1 + t15 / 0.2e1) * t199 + (t6 / 0.2e1 + t13 / 0.2e1) * t193 + (t5 / 0.2e1 + t12 / 0.2e1) * t191 + m(7) * (t16 * t25 + t17 * t29 + t18 * t30) + m(6) * (t26 * t60 + t31 * t64 + t32 * t65); m(6) * (t191 * t65 + t193 * t64 + t199 * t60) + m(7) * (t191 * t30 + t193 * t29 + t199 * t25); (t8 + t14) * t199 + (t4 + t11) * t193 + (t3 + t10) * t191 + m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t60 ^ 2 + t64 ^ 2 + t65 ^ 2); t33 + m(7) * (t46 * t49 + t47 * t50) + t267 * t164 + t268 * t162; m(7) * (t37 * t228 + (t232 * t46 - t234 * t47) * t226); m(7) * (t16 * t37 + t17 * t46 + t18 * t47) + t6 * t298 + t5 * t299 + t202 * t301 + t201 * t302 + t209 * t300 + t9 * t297; m(7) * (t191 * t47 + t193 * t46 + t199 * t37); t199 * t300 + m(7) * (t25 * t37 + t29 * t46 + t30 * t47) + t193 * t301 + t8 * t297 + t191 * t302 + t3 * t299 + t4 * t298; t164 * t2 + t162 * t1 + t178 * t7 + m(7) * (t37 ^ 2 + t46 ^ 2 + t47 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t61(1) t61(2) t61(4) t61(7) t61(11) t61(16); t61(2) t61(3) t61(5) t61(8) t61(12) t61(17); t61(4) t61(5) t61(6) t61(9) t61(13) t61(18); t61(7) t61(8) t61(9) t61(10) t61(14) t61(19); t61(11) t61(12) t61(13) t61(14) t61(15) t61(20); t61(16) t61(17) t61(18) t61(19) t61(20) t61(21);];
Mq  = res;
