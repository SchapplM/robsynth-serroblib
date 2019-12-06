% Calculate time derivative of joint inertia matrix for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:45
% EndTime: 2019-12-05 15:59:58
% DurationCPUTime: 7.35s
% Computational Cost: add. (12107->512), mult. (20170->808), div. (0->0), fcn. (19189->8), ass. (0->263)
t212 = sin(qJ(2));
t214 = cos(qJ(2));
t350 = ((Icges(4,5) - Icges(3,6)) * t214 + (Icges(4,4) - Icges(3,5)) * t212) * qJD(2);
t209 = sin(pkin(8));
t205 = t209 ^ 2;
t210 = cos(pkin(8));
t206 = t210 ^ 2;
t341 = t205 + t206;
t349 = -t341 + 0.1e1;
t340 = t341 * qJD(2);
t348 = m(3) * (rSges(3,1) * t212 + rSges(3,2) * t214) * t340;
t347 = t350 * t209;
t346 = t350 * t210;
t211 = sin(qJ(4));
t189 = -t214 * t211 * pkin(4) + pkin(7) * t212;
t337 = 2 * m(5);
t336 = 2 * m(6);
t335 = 0.2e1 * qJD(2);
t334 = m(4) / 0.2e1;
t333 = m(5) / 0.2e1;
t332 = m(6) / 0.2e1;
t331 = t209 / 0.2e1;
t330 = -t210 / 0.2e1;
t329 = t212 / 0.2e1;
t213 = cos(qJ(4));
t328 = pkin(4) * t213;
t208 = qJ(4) + qJ(5);
t203 = sin(t208);
t204 = cos(t208);
t305 = t209 * t212;
t169 = t203 * t210 + t204 * t305;
t170 = t203 * t305 - t204 * t210;
t304 = t209 * t214;
t111 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t304;
t113 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t304;
t115 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t304;
t303 = t210 * t212;
t167 = -t203 * t209 + t204 * t303;
t168 = t203 * t303 + t204 * t209;
t302 = t210 * t214;
t47 = t111 * t302 + t113 * t167 + t115 * t168;
t326 = t209 * t47;
t300 = t212 * t213;
t195 = t209 * t300 + t210 * t211;
t301 = t211 * t212;
t196 = t209 * t301 - t210 * t213;
t122 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t304;
t124 = Icges(5,4) * t196 + Icges(5,2) * t195 + Icges(5,6) * t304;
t126 = Icges(5,1) * t196 + Icges(5,4) * t195 + Icges(5,5) * t304;
t193 = -t209 * t211 + t210 * t300;
t194 = t209 * t213 + t210 * t301;
t55 = t122 * t302 + t124 * t193 + t126 * t194;
t325 = t209 * t55;
t110 = Icges(6,5) * t168 + Icges(6,6) * t167 + Icges(6,3) * t302;
t112 = Icges(6,4) * t168 + Icges(6,2) * t167 + Icges(6,6) * t302;
t114 = Icges(6,1) * t168 + Icges(6,4) * t167 + Icges(6,5) * t302;
t48 = t110 * t304 + t112 * t169 + t114 * t170;
t324 = t210 * t48;
t121 = Icges(5,5) * t194 + Icges(5,6) * t193 + Icges(5,3) * t302;
t123 = Icges(5,4) * t194 + Icges(5,2) * t193 + Icges(5,6) * t302;
t125 = Icges(5,1) * t194 + Icges(5,4) * t193 + Icges(5,5) * t302;
t56 = t121 * t304 + t123 * t195 + t125 * t196;
t323 = t210 * t56;
t116 = rSges(6,1) * t168 + rSges(6,2) * t167 + rSges(6,3) * t302;
t290 = qJD(2) * t212;
t278 = t209 * t290;
t207 = qJD(4) + qJD(5);
t289 = qJD(2) * t214;
t231 = t207 * t210 + t209 * t289;
t281 = t207 * t305;
t139 = t203 * t231 + t204 * t281;
t140 = -t203 * t281 + t204 * t231;
t81 = rSges(6,1) * t139 + rSges(6,2) * t140 - rSges(6,3) * t278;
t322 = t116 * t278 + t81 * t302;
t287 = qJD(4) * t213;
t283 = pkin(4) * t287;
t216 = -t189 * qJD(2) + t212 * t283;
t288 = qJD(4) * t211;
t284 = pkin(4) * t288;
t118 = t209 * t216 + t210 * t284;
t321 = t118 + t81;
t119 = -t209 * t284 + t210 * t216;
t232 = -t207 * t209 + t210 * t289;
t280 = t207 * t303;
t137 = t203 * t232 + t204 * t280;
t138 = -t203 * t280 + t204 * t232;
t277 = t210 * t290;
t80 = rSges(6,1) * t137 + rSges(6,2) * t138 - rSges(6,3) * t277;
t320 = t119 + t80;
t317 = Icges(5,4) * t211;
t316 = Icges(5,4) * t213;
t315 = Icges(6,4) * t203;
t314 = Icges(6,4) * t204;
t254 = Icges(6,5) * t203 + Icges(6,6) * t204;
t306 = t207 * t214;
t311 = ((-Icges(6,5) * t204 + Icges(6,6) * t203) * t306 + (Icges(6,3) * t214 + t212 * t254) * qJD(2)) * t212;
t255 = Icges(5,5) * t211 + Icges(5,6) * t213;
t286 = qJD(4) * t214;
t310 = ((-Icges(5,5) * t213 + Icges(5,6) * t211) * t286 + (Icges(5,3) * t214 + t212 * t255) * qJD(2)) * t212;
t155 = Icges(6,3) * t212 - t214 * t254;
t309 = t155 * t212;
t171 = Icges(5,3) * t212 - t214 * t255;
t308 = t171 * t212;
t267 = rSges(6,1) * t203 + rSges(6,2) * t204;
t109 = (-rSges(6,1) * t204 + rSges(6,2) * t203) * t306 + (rSges(6,3) * t214 + t212 * t267) * qJD(2);
t220 = pkin(4) * t301 + pkin(7) * t214;
t144 = t220 * qJD(2) - t214 * t283;
t299 = -t109 - t144;
t134 = t328 * t209 + t220 * t210;
t298 = t116 + t134;
t117 = rSges(6,1) * t170 + rSges(6,2) * t169 + rSges(6,3) * t304;
t135 = t220 * t209 - t328 * t210;
t297 = t117 + t135;
t199 = pkin(2) * t212 - qJ(3) * t214;
t296 = t341 * (-t199 * qJD(2) + qJD(3) * t212);
t158 = rSges(6,3) * t212 - t214 * t267;
t295 = -t158 - t189;
t264 = pkin(2) * t214 + qJ(3) * t212;
t294 = t341 * t264;
t190 = t264 * qJD(2) - qJD(3) * t214;
t266 = -rSges(4,2) * t214 + rSges(4,3) * t212;
t293 = -t266 * qJD(2) - t190;
t265 = rSges(4,2) * t212 + rSges(4,3) * t214;
t292 = -t199 + t265;
t282 = t116 * t289 + t158 * t277 + t212 * t80;
t256 = Icges(6,2) * t204 + t315;
t107 = (Icges(6,2) * t203 - t314) * t306 + (Icges(6,6) * t214 + t212 * t256) * qJD(2);
t259 = Icges(6,1) * t203 + t314;
t108 = (-Icges(6,1) * t204 + t315) * t306 + (Icges(6,5) * t214 + t212 * t259) * qJD(2);
t156 = Icges(6,6) * t212 - t214 * t256;
t157 = Icges(6,5) * t212 - t214 * t259;
t251 = t112 * t204 + t114 * t203;
t74 = Icges(6,5) * t137 + Icges(6,6) * t138 - Icges(6,3) * t277;
t76 = Icges(6,4) * t137 + Icges(6,2) * t138 - Icges(6,6) * t277;
t78 = Icges(6,1) * t137 + Icges(6,4) * t138 - Icges(6,5) * t277;
t17 = (t251 * qJD(2) + t74) * t212 + (qJD(2) * t110 + (-t114 * t207 - t76) * t204 + (t112 * t207 - t78) * t203) * t214;
t250 = t113 * t204 + t115 * t203;
t75 = Icges(6,5) * t139 + Icges(6,6) * t140 - Icges(6,3) * t278;
t77 = Icges(6,4) * t139 + Icges(6,2) * t140 - Icges(6,6) * t278;
t79 = Icges(6,1) * t139 + Icges(6,4) * t140 - Icges(6,5) * t278;
t18 = (t250 * qJD(2) + t75) * t212 + (qJD(2) * t111 + (-t115 * t207 - t77) * t204 + (t113 * t207 - t79) * t203) * t214;
t246 = t156 * t204 + t157 * t203;
t52 = t110 * t212 - t214 * t251;
t53 = t111 * t212 - t214 * t250;
t263 = t53 * t209 + t52 * t210;
t236 = -t110 * t290 + t214 * t74;
t19 = t112 * t138 + t114 * t137 + t167 * t76 + t168 * t78 + t210 * t236;
t235 = -t111 * t290 + t214 * t75;
t20 = t113 * t138 + t115 * t137 + t167 * t77 + t168 * t79 + t210 * t235;
t46 = t110 * t302 + t112 * t167 + t114 * t168;
t65 = t155 * t302 + t156 * t167 + t157 * t168;
t5 = (t107 * t167 + t108 * t168 + t137 * t157 + t138 * t156) * t212 + (t20 * t209 + (t19 + t311) * t210) * t214 + (t65 * t214 + (-t326 + (-t46 - t309) * t210) * t212) * qJD(2);
t21 = t112 * t140 + t114 * t139 + t169 * t76 + t170 * t78 + t209 * t236;
t22 = t113 * t140 + t115 * t139 + t169 * t77 + t170 * t79 + t209 * t235;
t49 = t111 * t304 + t113 * t169 + t115 * t170;
t66 = t155 * t304 + t156 * t169 + t157 * t170;
t6 = (t107 * t169 + t108 * t170 + t139 * t157 + t140 * t156) * t212 + (t21 * t210 + (t22 + t311) * t209) * t214 + (t66 * t214 + (-t324 + (-t49 - t309) * t209) * t212) * qJD(2);
t82 = -t214 * t246 + t309;
t279 = t5 * t302 + t6 * t304 + (t212 * t82 + t214 * t263) * t289 + t212 * ((t311 + (t212 * t246 - t263) * qJD(2)) * t212 + (t17 * t210 + t18 * t209 + (qJD(2) * t155 + (-t157 * t207 - t107) * t204 + (t156 * t207 - t108) * t203) * t212 + t82 * qJD(2)) * t214);
t276 = t211 * t289;
t275 = t213 * t289;
t273 = -t212 * pkin(6) - t199;
t272 = t294 + (t209 * t304 + t210 * t302) * pkin(6);
t268 = rSges(5,1) * t211 + rSges(5,2) * t213;
t174 = rSges(5,3) * t212 - t214 * t268;
t271 = -t174 + t273;
t270 = -pkin(6) * t289 - t190;
t249 = t123 * t213 + t125 * t211;
t61 = t121 * t212 - t214 * t249;
t248 = t124 * t213 + t126 * t211;
t62 = t122 * t212 - t214 * t248;
t262 = t62 * t209 + t61 * t210;
t260 = Icges(5,1) * t211 + t316;
t257 = Icges(5,2) * t213 + t317;
t129 = rSges(5,1) * t194 + rSges(5,2) * t193 + rSges(5,3) * t302;
t130 = rSges(5,1) * t196 + rSges(5,2) * t195 + rSges(5,3) * t304;
t247 = t129 * t209 - t130 * t210;
t172 = Icges(5,6) * t212 - t214 * t257;
t173 = Icges(5,5) * t212 - t214 * t260;
t243 = t172 * t213 + t173 * t211;
t238 = t273 + t295;
t136 = (-rSges(5,1) * t213 + rSges(5,2) * t211) * t286 + (rSges(5,3) * t214 + t212 * t268) * qJD(2);
t237 = -t136 + t270;
t148 = qJD(4) * t193 + t210 * t276;
t149 = -qJD(4) * t194 + t210 * t275;
t91 = Icges(5,5) * t148 + Icges(5,6) * t149 - Icges(5,3) * t277;
t234 = -t121 * t290 + t214 * t91;
t150 = qJD(4) * t195 + t209 * t276;
t151 = -qJD(4) * t196 + t209 * t275;
t92 = Icges(5,5) * t150 + Icges(5,6) * t151 - Icges(5,3) * t278;
t233 = -t122 * t290 + t214 * t92;
t230 = t270 + t299;
t13 = t19 * t209 - t20 * t210;
t14 = t209 * t21 - t210 * t22;
t229 = t6 * t330 + t5 * t331 + (t17 * t209 - t18 * t210) * t329 + t14 * t304 / 0.2e1 + t13 * t302 / 0.2e1 + (t209 * t52 - t210 * t53) * t289 / 0.2e1 - (t209 * (t209 * t48 - t210 * t49) + t210 * (t209 * t46 - t210 * t47)) * t290 / 0.2e1;
t219 = -pkin(6) * t290 * t341 + t296;
t23 = t212 * t65 + (t210 * t46 + t326) * t214;
t24 = t212 * t66 + (t209 * t49 + t324) * t214;
t217 = (-t209 * t24 - t210 * t23) * t290 + t279;
t154 = t292 * t210;
t153 = t292 * t209;
t147 = t158 * t304;
t142 = t293 * t210;
t141 = t293 * t209;
t133 = (-Icges(5,1) * t213 + t317) * t286 + (Icges(5,5) * t214 + t212 * t260) * qJD(2);
t132 = (Icges(5,2) * t211 - t316) * t286 + (Icges(5,6) * t214 + t212 * t257) * qJD(2);
t128 = t271 * t210;
t127 = t271 * t209;
t105 = t212 * t116;
t103 = t117 * t302;
t102 = t109 * t304;
t100 = t238 * t210;
t99 = t238 * t209;
t98 = rSges(5,1) * t150 + rSges(5,2) * t151 - rSges(5,3) * t278;
t97 = rSges(5,1) * t148 + rSges(5,2) * t149 - rSges(5,3) * t277;
t96 = Icges(5,1) * t150 + Icges(5,4) * t151 - Icges(5,5) * t278;
t95 = Icges(5,1) * t148 + Icges(5,4) * t149 - Icges(5,5) * t277;
t94 = Icges(5,4) * t150 + Icges(5,2) * t151 - Icges(5,6) * t278;
t93 = Icges(5,4) * t148 + Icges(5,2) * t149 - Icges(5,6) * t277;
t90 = t129 * t212 - t174 * t302;
t89 = -t130 * t212 + t174 * t304;
t88 = t237 * t210;
t87 = t237 * t209;
t86 = -t214 * t243 + t308;
t85 = -t158 * t302 + t105;
t84 = -t117 * t212 + t147;
t83 = t341 * t266 + t294;
t71 = t171 * t304 + t172 * t195 + t173 * t196;
t70 = t171 * t302 + t172 * t193 + t173 * t194;
t69 = t247 * t214;
t68 = t265 * t340 + t296;
t67 = -t116 * t304 + t103;
t64 = t230 * t210;
t63 = t230 * t209;
t60 = t134 * t212 + t295 * t302 + t105;
t59 = t189 * t304 - t297 * t212 + t147;
t58 = t129 * t210 + t130 * t209 + t272;
t57 = t122 * t304 + t124 * t195 + t126 * t196;
t54 = t121 * t302 + t123 * t193 + t125 * t194;
t51 = -t136 * t302 + t212 * t97 + (t129 * t214 + t174 * t303) * qJD(2);
t50 = t136 * t304 - t212 * t98 + (-t130 * t214 - t174 * t305) * qJD(2);
t45 = t103 + (t135 * t210 - t298 * t209) * t214;
t44 = t209 * t98 + t210 * t97 + t219;
t43 = -t109 * t302 + t282;
t42 = -t212 * t81 + t102 + (-t117 * t214 - t158 * t305) * qJD(2);
t41 = (-t209 * t97 + t210 * t98) * t214 + t247 * t290;
t40 = t297 * t209 + t298 * t210 + t272;
t39 = -t117 * t277 - t80 * t304 + t322;
t37 = t321 * t209 + t320 * t210 + t219;
t36 = t134 * t289 + t119 * t212 + (t189 * t290 + t299 * t214) * t210 + t282;
t35 = t144 * t304 + t102 - t321 * t212 + (-t297 * t214 + t295 * t305) * qJD(2);
t32 = t124 * t151 + t126 * t150 + t195 * t94 + t196 * t96 + t209 * t233;
t31 = t123 * t151 + t125 * t150 + t195 * t93 + t196 * t95 + t209 * t234;
t30 = t124 * t149 + t126 * t148 + t193 * t94 + t194 * t96 + t210 * t233;
t29 = t123 * t149 + t125 * t148 + t193 * t93 + t194 * t95 + t210 * t234;
t27 = (t248 * qJD(2) + t92) * t212 + (qJD(2) * t122 - t211 * t96 - t213 * t94 + (t124 * t211 - t126 * t213) * qJD(4)) * t214;
t26 = (t249 * qJD(2) + t91) * t212 + (qJD(2) * t121 - t211 * t95 - t213 * t93 + (t123 * t211 - t125 * t213) * qJD(4)) * t214;
t25 = (t118 * t210 - t320 * t209) * t214 + (t134 * t209 - t297 * t210) * t290 + t322;
t16 = t209 * t31 - t210 * t32;
t15 = t209 * t29 - t210 * t30;
t9 = (t132 * t195 + t133 * t196 + t150 * t173 + t151 * t172) * t212 + (t31 * t210 + (t32 + t310) * t209) * t214 + (t71 * t214 + (-t323 + (-t57 - t308) * t209) * t212) * qJD(2);
t8 = (t132 * t193 + t133 * t194 + t148 * t173 + t149 * t172) * t212 + (t30 * t209 + (t29 + t310) * t210) * t214 + (t70 * t214 + (-t325 + (-t54 - t308) * t210) * t212) * qJD(2);
t1 = [0; m(4) * t68 + m(5) * t44 + m(6) * t37 - t348; (t100 * t64 + t37 * t40 + t63 * t99) * t336 + (t127 * t87 + t128 * t88 + t44 * t58) * t337 + 0.2e1 * m(4) * (t141 * t153 + t142 * t154 + t68 * t83) + (-t347 * t206 - t14 - t16) * t210 + (t13 + t15 + t346 * t205 + (-t347 * t209 + t346 * t210) * t210) * t209 + 0.2e1 * t349 * (rSges(3,1) * t214 - rSges(3,2) * t212) * t348; (m(4) + m(5) + m(6)) * t290; m(6) * (-t214 * t37 + t64 * t303 + t63 * t305) + m(5) * (-t214 * t44 + t88 * t303 + t87 * t305) + m(4) * (t141 * t305 + t142 * t303 - t214 * t68) + ((t100 * t302 + t212 * t40 + t99 * t304) * t332 + (t127 * t304 + t128 * t302 + t212 * t58) * t333 + (t153 * t304 + t154 * t302 + t212 * t83) * t334) * t335; -0.4e1 * (t334 + t333 + t332) * t349 * t212 * t289; m(5) * t41 + m(6) * t25; t8 * t331 + t9 * t330 + (t209 * t26 - t210 * t27) * t329 + (t16 * t331 + t210 * t15 / 0.2e1) * t214 + m(6) * (t100 * t35 + t25 * t40 + t36 * t99 + t37 * t45 + t59 * t64 + t60 * t63) + m(5) * (t127 * t51 + t128 * t50 + t41 * t58 - t44 * t69 + t87 * t90 + t88 * t89) + (t214 * (t209 * t61 - t210 * t62) / 0.2e1 + (-t209 * (t209 * t56 - t210 * t57) / 0.2e1 + (t209 * t54 - t210 * t55) * t330) * t212) * qJD(2) + t229; m(5) * (-t214 * t41 + t50 * t303 + t51 * t305) + m(6) * (-t214 * t25 + t35 * t303 + t36 * t305) + ((-t212 * t69 + t89 * t302 + t90 * t304) * t333 + (t212 * t45 + t59 * t302 + t60 * t304) * t332) * t335; t9 * t304 + t8 * t302 + (t25 * t45 + t35 * t59 + t36 * t60) * t336 + (t212 * t86 + t214 * t262) * t289 + t212 * ((t310 + (t212 * t243 - t262) * qJD(2)) * t212 + (t26 * t210 + t27 * t209 + (qJD(2) * t171 - t132 * t213 - t133 * t211 + t172 * t288 - t173 * t287) * t212 + t86 * qJD(2)) * t214) + (-t41 * t69 + t50 * t89 + t51 * t90) * t337 + t279 + (-t24 - t212 * t71 - (t209 * t57 + t323) * t214) * t278 + (-t23 - t212 * t70 - (t210 * t54 + t325) * t214) * t277; m(6) * t39; m(6) * (t100 * t42 + t37 * t67 + t39 * t40 + t43 * t99 + t63 * t85 + t64 * t84) + t229; m(6) * (-t214 * t39 + (t209 * t43 + t210 * t42) * t212 + (t212 * t67 + (t209 * t85 + t210 * t84) * t214) * qJD(2)); m(6) * (t25 * t67 + t35 * t84 + t36 * t85 + t39 * t45 + t42 * t59 + t43 * t60) + t217; (t39 * t67 + t42 * t84 + t43 * t85) * t336 + t217;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
