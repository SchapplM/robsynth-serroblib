% Calculate time derivative of joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:33
% DurationCPUTime: 5.59s
% Computational Cost: add. (4812->408), mult. (5516->573), div. (0->0), fcn. (4274->8), ass. (0->240)
t166 = sin(qJ(1));
t168 = cos(qJ(1));
t161 = qJ(2) + pkin(7);
t151 = qJ(4) + t161;
t143 = sin(t151);
t284 = rSges(5,2) * t143;
t144 = cos(t151);
t287 = rSges(5,1) * t144;
t216 = -t284 + t287;
t84 = t166 * rSges(5,3) + t168 * t216;
t193 = Icges(5,5) * t144 - Icges(5,6) * t143;
t77 = -Icges(5,3) * t168 + t166 * t193;
t313 = qJD(1) * t77;
t283 = rSges(5,2) * t144;
t117 = rSges(5,1) * t143 + t283;
t160 = qJD(2) + qJD(4);
t183 = t117 * t160;
t165 = sin(qJ(2));
t167 = cos(qJ(2));
t275 = Icges(3,4) * t167;
t200 = -Icges(3,2) * t165 + t275;
t104 = Icges(3,6) * t166 + t168 * t200;
t276 = Icges(3,4) * t165;
t205 = Icges(3,1) * t167 - t276;
t106 = Icges(3,5) * t166 + t168 * t205;
t191 = t104 * t165 - t106 * t167;
t312 = t166 * t191;
t149 = sin(t161);
t150 = cos(t161);
t273 = Icges(4,4) * t150;
t198 = -Icges(4,2) * t149 + t273;
t92 = Icges(4,6) * t166 + t168 * t198;
t274 = Icges(4,4) * t149;
t203 = Icges(4,1) * t150 - t274;
t94 = Icges(4,5) * t166 + t168 * t203;
t210 = t149 * t92 - t150 * t94;
t311 = t166 * t210;
t271 = Icges(5,4) * t144;
t196 = -Icges(5,2) * t143 + t271;
t80 = Icges(5,6) * t166 + t168 * t196;
t272 = Icges(5,4) * t143;
t201 = Icges(5,1) * t144 - t272;
t82 = Icges(5,5) * t166 + t168 * t201;
t214 = t143 * t80 - t144 * t82;
t310 = t166 * t214;
t157 = t167 * pkin(2);
t145 = t157 + pkin(1);
t292 = pkin(1) - t145;
t309 = t166 * t292;
t103 = -Icges(3,6) * t168 + t166 * t200;
t105 = -Icges(3,5) * t168 + t166 * t205;
t192 = t103 * t165 - t105 * t167;
t308 = t168 * t192;
t91 = -Icges(4,6) * t168 + t166 * t198;
t93 = -Icges(4,5) * t168 + t166 * t203;
t211 = t149 * t91 - t150 * t93;
t307 = t168 * t211;
t79 = -Icges(5,6) * t168 + t166 * t196;
t81 = -Icges(5,5) * t168 + t166 * t201;
t215 = t143 * t79 - t144 * t81;
t306 = t168 * t215;
t164 = -qJ(3) - pkin(5);
t159 = -pkin(6) + t164;
t245 = t159 - t164;
t305 = t245 * t168;
t155 = t166 * rSges(4,3);
t285 = rSges(4,2) * t149;
t304 = -t168 * t285 + t155;
t115 = Icges(5,2) * t144 + t272;
t116 = Icges(5,1) * t143 + t271;
t190 = t115 * t143 - t116 * t144;
t303 = qJD(1) * t190 + t193 * t160;
t194 = Icges(4,5) * t150 - Icges(4,6) * t149;
t89 = -Icges(4,3) * t168 + t166 * t194;
t195 = Icges(3,5) * t167 - Icges(3,6) * t165;
t101 = -Icges(3,3) * t168 + t166 * t195;
t302 = 2 * m(3);
t301 = 2 * m(4);
t300 = 2 * m(5);
t162 = t166 ^ 2;
t163 = t168 ^ 2;
t299 = t166 / 0.2e1;
t298 = -t168 / 0.2e1;
t135 = rSges(3,1) * t165 + rSges(3,2) * t167;
t297 = m(3) * t135;
t296 = pkin(2) * t165;
t295 = pkin(3) * t149;
t294 = pkin(3) * t150;
t293 = t166 * pkin(5);
t158 = t168 * pkin(5);
t291 = -pkin(5) - t164;
t83 = -rSges(5,3) * t168 + t166 * t216;
t36 = t166 * t83 + t168 * t84;
t87 = t164 * t168 + t158 - t309;
t137 = t168 * t145;
t88 = -pkin(1) * t168 + t166 * t291 + t137;
t290 = t166 * t87 + t168 * t88;
t289 = rSges(3,1) * t167;
t288 = rSges(4,1) * t150;
t286 = rSges(3,2) * t165;
t282 = rSges(3,3) * t168;
t281 = t149 * t93;
t280 = t149 * t94;
t279 = t150 * t91;
t278 = t150 * t92;
t156 = t166 * rSges(3,3);
t277 = rSges(5,3) - t159;
t78 = Icges(5,3) * t166 + t168 * t193;
t261 = qJD(1) * t78;
t90 = Icges(4,3) * t166 + t168 * t194;
t260 = qJD(1) * t90;
t259 = t103 * t167;
t258 = t104 * t167;
t257 = t105 * t165;
t256 = t106 * t165;
t255 = t115 * t160;
t254 = t116 * t160;
t253 = t143 * t160;
t252 = t144 * t160;
t251 = t160 * t168;
t250 = t166 * t159;
t241 = qJD(1) * t168;
t242 = qJD(1) * t166;
t249 = rSges(5,3) * t241 + t242 * t284;
t248 = rSges(4,3) * t241 + t242 * t285;
t238 = qJD(2) * t165;
t231 = pkin(2) * t238;
t247 = t164 * t242 + t166 * t231;
t246 = t168 * t289 + t156;
t244 = t162 + t163;
t102 = Icges(3,3) * t166 + t168 * t195;
t243 = qJD(1) * t102;
t240 = qJD(2) * t149;
t239 = qJD(2) * t150;
t237 = qJD(2) * t166;
t236 = qJD(2) * t167;
t235 = t166 * (t84 * qJD(1) - t166 * t183) + t168 * (-t251 * t283 + (-t143 * t251 - t144 * t242) * rSges(5,1) + t249) + t83 * t241;
t152 = qJD(3) * t166;
t221 = t168 * t231;
t153 = qJD(3) * t168;
t230 = -t153 - t247;
t234 = t166 * ((-t168 * t292 - t293) * qJD(1) + t230) + t168 * (-t221 + t152 + (t168 * t291 + t309) * qJD(1)) + t87 * t241;
t233 = t168 * t286;
t229 = t165 * t242;
t123 = rSges(4,1) * t149 + rSges(4,2) * t150;
t226 = -t123 - t296;
t39 = -qJD(1) * t79 - t168 * t255;
t225 = t160 * t82 + t39;
t40 = qJD(1) * t80 - t166 * t255;
t224 = t160 * t81 + t40;
t41 = -qJD(1) * t81 - t168 * t254;
t223 = -t160 * t80 + t41;
t42 = qJD(1) * t82 - t166 * t254;
t222 = t160 * t79 - t42;
t14 = -t166 * t215 - t168 * t77;
t15 = -t168 * t78 - t310;
t16 = t166 * t77 - t306;
t17 = t166 * t78 - t168 * t214;
t114 = Icges(5,5) * t143 + Icges(5,6) * t144;
t178 = t160 * t114;
t37 = -t168 * t178 - t313;
t38 = -t166 * t178 + t261;
t220 = -t168 * ((t168 * t38 + (t15 + t306) * qJD(1)) * t168 + (t14 * qJD(1) + (-t143 * t39 + t144 * t41 - t252 * t80 - t253 * t82 + t261) * t166 + (-t37 + t222 * t144 + t224 * t143 + (-t214 - t77) * qJD(1)) * t168) * t166) + t166 * ((t166 * t37 + (t16 + t310) * qJD(1)) * t166 + (t17 * qJD(1) + (t143 * t40 - t144 * t42 + t252 * t79 + t253 * t81 - t313) * t168 + (-t38 + t223 * t144 - t225 * t143 + (-t215 + t78) * qJD(1)) * t166) * t168) + (-t14 * t168 + t15 * t166) * t242 + (-t16 * t168 + t17 * t166) * t241;
t219 = -t295 - t296;
t218 = -t286 + t289;
t217 = -t285 + t288;
t127 = t145 + t294;
t186 = -t127 - t216;
t53 = t166 * t186 + t168 * t277;
t118 = t168 * t127;
t54 = t118 + t84 - t250;
t207 = t166 * t54 + t168 * t53;
t188 = -t117 + t219;
t67 = t188 * t166;
t68 = t188 * t168;
t206 = t166 * t67 + t168 * t68;
t204 = Icges(3,1) * t165 + t275;
t202 = Icges(4,1) * t149 + t273;
t199 = Icges(3,2) * t167 + t276;
t197 = Icges(4,2) * t150 + t274;
t99 = t168 * t288 + t304;
t189 = -pkin(1) - t218;
t187 = -t145 - t217;
t96 = t196 * t160;
t97 = t201 * t160;
t169 = qJD(1) * t114 + (t97 - t255) * t144 + (-t96 - t254) * t143;
t184 = (t143 * t223 + t144 * t225 + t166 * t303 + t169 * t168) * t299 + (-t143 * t222 + t144 * t224 + t169 * t166 - t168 * t303) * t298 + (-t114 * t168 + t143 * t81 + t144 * t79 - t166 * t190) * t242 / 0.2e1 + (t166 * t114 + t143 * t82 + t144 * t80 - t168 * t190) * t241 / 0.2e1;
t182 = qJD(2) * t135;
t181 = qJD(2) * t123;
t177 = qJD(2) * t204;
t176 = qJD(2) * t202;
t175 = qJD(2) * t199;
t174 = qJD(2) * t197;
t173 = qJD(2) * (-Icges(3,5) * t165 - Icges(3,6) * t167);
t172 = qJD(2) * (-Icges(4,5) * t149 - Icges(4,6) * t150);
t100 = t216 * t160;
t171 = -t100 + (-t157 - t294) * qJD(2);
t170 = rSges(3,2) * t229 + rSges(3,3) * t241 - t168 * t182;
t140 = pkin(2) * t229;
t128 = t218 * qJD(2);
t122 = t219 * qJD(2);
t113 = t217 * qJD(2);
t109 = t168 * t122;
t108 = -t233 + t246;
t107 = t166 * t218 - t282;
t98 = -rSges(4,3) * t168 + t166 * t217;
t86 = t226 * t168;
t85 = t226 * t166;
t76 = t293 + (pkin(1) - t286) * t168 + t246;
t75 = t166 * t189 + t158 + t282;
t66 = -t166 * t245 + t118 - t137;
t65 = t305 + (t127 - t145) * t166;
t64 = -t166 * t164 + t137 + t99;
t63 = (rSges(4,3) - t164) * t168 + t187 * t166;
t58 = t166 * t173 + t243;
t57 = -qJD(1) * t101 + t168 * t173;
t48 = t166 * t172 + t260;
t47 = -qJD(1) * t89 + t168 * t172;
t46 = t135 * t237 + ((-rSges(3,3) - pkin(5)) * t166 + t189 * t168) * qJD(1);
t45 = (t158 + (-pkin(1) - t289) * t166) * qJD(1) + t170;
t44 = -t123 * t241 - t166 * t113 + (-t165 * t241 - t166 * t236) * pkin(2);
t43 = t123 * t242 + t140 + (-pkin(2) * t236 - t113) * t168;
t31 = t166 * t102 - t191 * t168;
t30 = t166 * t101 - t308;
t29 = -t102 * t168 - t312;
t28 = -t101 * t168 - t166 * t192;
t25 = t123 * t237 + (t168 * t187 - t155) * qJD(1) - t230;
t24 = t152 + (-t145 - t288) * t242 + (-qJD(1) * t164 + qJD(2) * t226) * t168 + t248;
t23 = qJD(1) * t68 + t166 * t171;
t22 = t140 + (t117 + t295) * t242 + t171 * t168;
t21 = t166 * t90 - t210 * t168;
t20 = t166 * t89 - t307;
t19 = -t168 * t90 - t311;
t18 = -t166 * t211 - t168 * t89;
t13 = t153 + (-t122 + t183) * t166 + (-t166 * t277 + t168 * t186) * qJD(1);
t12 = t109 + t152 - t168 * t183 + (-t159 * t168 + (-t127 - t287) * t166) * qJD(1) + t249;
t11 = -t242 * t84 + t235;
t10 = t166 * t65 + t168 * t66 + t290 + t36;
t3 = t166 * (t166 * t122 + t247) + t168 * (t109 + t221) + ((t65 - t305) * t168 + (-t66 - t84 - t88 - t250) * t166) * qJD(1) + t234 + t235;
t1 = [(t45 * t76 + t46 * t75) * t302 + (t24 * t64 + t25 * t63) * t301 + t116 * t252 + t143 * t97 - t115 * t253 + t144 * t96 + (t12 * t54 + t13 * t53) * t300 + (t203 - t197) * t240 + (t202 + t198) * t239 + (t205 - t199) * t238 + (t204 + t200) * t236; t184 + m(4) * (t24 * t85 + t25 * t86 + t43 * t63 + t44 * t64) + m(5) * (t12 * t67 + t13 * t68 + t22 * t53 + t23 * t54) + ((t280 / 0.2e1 + t278 / 0.2e1 + t258 / 0.2e1 + t256 / 0.2e1 - t76 * t297) * t168 + (t75 * t297 + t259 / 0.2e1 + t257 / 0.2e1 + t281 / 0.2e1 + t279 / 0.2e1) * t166) * qJD(1) + m(3) * ((-t166 * t45 - t168 * t46) * t135 + (-t166 * t76 - t168 * t75) * t128) + (t149 * (-qJD(1) * t93 - t168 * t176) + t150 * (-qJD(1) * t91 - t168 * t174) + t165 * (-qJD(1) * t105 - t168 * t177) + t167 * (-qJD(1) * t103 - t168 * t175) + (-t191 - t210) * qJD(2)) * t299 + (t149 * (qJD(1) * t94 - t166 * t176) + t150 * (qJD(1) * t92 - t166 * t174) + t165 * (qJD(1) * t106 - t166 * t177) + t167 * (qJD(1) * t104 - t166 * t175) + (-t192 - t211) * qJD(2)) * t298 + (t195 + t194) * qJD(2) * (t162 / 0.2e1 + t163 / 0.2e1); (t10 * t3 + t22 * t68 + t23 * t67) * t300 + (t86 * t43 + t85 * t44 + (t166 * t98 + t168 * t99 + t290) * ((qJD(1) * t98 - t168 * t181 + t248) * t168 + (-t166 * t181 + (-t88 - t99 + t304) * qJD(1)) * t166 + t234)) * t301 + t166 * ((t166 * t57 + (t30 + t312) * qJD(1)) * t166 + (t31 * qJD(1) + (t103 * t236 + t105 * t238) * t168 + (-t58 + (-t256 - t258) * qJD(2) + (t102 - t192) * qJD(1)) * t166) * t168) - t168 * ((t168 * t48 + (t19 + t307) * qJD(1)) * t168 + (t18 * qJD(1) + (-t239 * t92 - t240 * t94 + t260) * t166 + (-t47 + (t279 + t281) * qJD(2) - t210 * qJD(1)) * t168) * t166) - t168 * ((t168 * t58 + (t29 + t308) * qJD(1)) * t168 + (t28 * qJD(1) + (-t104 * t236 - t106 * t238 + t243) * t166 + (-t57 + (t257 + t259) * qJD(2) - t191 * qJD(1)) * t168) * t166) + t166 * ((t166 * t47 + (t20 + t311) * qJD(1)) * t166 + (t21 * qJD(1) + (t239 * t91 + t240 * t93) * t168 + (-t48 + (-t278 - t280) * qJD(2) + (-t211 + t90) * qJD(1)) * t166) * t168) + ((t166 * t107 + t108 * t168) * ((qJD(1) * t107 + t170) * t168 + (-t166 * t182 + (-t108 - t233 + t156) * qJD(1)) * t166) + t244 * t135 * t128) * t302 + t220 + ((-t18 - t28) * t168 + (t19 + t29) * t166) * t242 + ((-t20 - t30) * t168 + (t21 + t31) * t166) * t241; m(4) * (t166 * t25 - t168 * t24 + (t166 * t64 + t168 * t63) * qJD(1)) + m(5) * (qJD(1) * t207 - t12 * t168 + t166 * t13); m(5) * (qJD(1) * t206 + t166 * t22 - t168 * t23) + m(4) * (t166 * t43 - t168 * t44 + (t166 * t85 + t168 * t86) * qJD(1)); 0; m(5) * (-t207 * t100 + (-t12 * t166 - t13 * t168 + (t166 * t53 - t168 * t54) * qJD(1)) * t117) + t184; m(5) * (t11 * t10 + t36 * t3 - t206 * t100 + (-t166 * t23 - t168 * t22 + (t166 * t68 - t168 * t67) * qJD(1)) * t117) + t220; 0; (t100 * t117 * t244 + t11 * t36) * t300 + t220;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
