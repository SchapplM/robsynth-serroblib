% Calculate time derivative of joint inertia matrix for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:17
% DurationCPUTime: 6.11s
% Computational Cost: add. (7580->386), mult. (5870->531), div. (0->0), fcn. (4504->8), ass. (0->232)
t335 = Icges(4,3) + Icges(5,3);
t168 = qJ(3) + pkin(9);
t160 = sin(t168);
t162 = cos(t168);
t170 = sin(qJ(3));
t171 = cos(qJ(3));
t334 = Icges(4,5) * t171 + Icges(5,5) * t162 - Icges(4,6) * t170 - Icges(5,6) * t160;
t166 = pkin(8) + qJ(2);
t159 = sin(t166);
t161 = cos(t166);
t332 = t335 * t159 + t334 * t161;
t333 = t334 * t159 - t335 * t161;
t280 = Icges(4,4) * t171;
t203 = -Icges(4,2) * t170 + t280;
t102 = Icges(4,6) * t159 + t161 * t203;
t281 = Icges(4,4) * t170;
t208 = Icges(4,1) * t171 - t281;
t104 = Icges(4,5) * t159 + t161 * t208;
t278 = Icges(5,4) * t162;
t201 = -Icges(5,2) * t160 + t278;
t92 = Icges(5,6) * t159 + t161 * t201;
t279 = Icges(5,4) * t160;
t206 = Icges(5,1) * t162 - t279;
t94 = Icges(5,5) * t159 + t161 * t206;
t310 = t102 * t170 - t104 * t171 + t160 * t92 - t162 * t94;
t331 = t159 * t310;
t330 = t310 * t161;
t163 = qJ(5) + t168;
t154 = sin(t163);
t289 = rSges(6,2) * t154;
t155 = cos(t163);
t292 = rSges(6,1) * t155;
t219 = -t289 + t292;
t88 = t159 * rSges(6,3) + t161 * t219;
t329 = (-Icges(4,5) * t170 - Icges(5,5) * t160 - Icges(4,6) * t171 - Icges(5,6) * t162) * qJD(3);
t101 = -Icges(4,6) * t161 + t159 * t203;
t103 = -Icges(4,5) * t161 + t159 * t208;
t91 = -Icges(5,6) * t161 + t159 * t201;
t93 = -Icges(5,5) * t161 + t159 * t206;
t328 = t101 * t170 - t103 * t171 + t160 * t91 - t162 * t93;
t326 = t332 * qJD(2);
t325 = t333 * t159 - t331 + (-t328 - t332) * t161;
t157 = t159 ^ 2;
t158 = t161 ^ 2;
t306 = 2 * m(5);
t126 = rSges(5,1) * t160 + rSges(5,2) * t162;
t185 = qJD(3) * t126;
t148 = qJD(4) * t159;
t240 = qJD(3) * t170;
t234 = pkin(3) * t240;
t224 = t161 * t234;
t149 = qJD(4) * t161;
t169 = -qJ(4) - pkin(6);
t245 = qJD(2) * t169;
t252 = (t234 + t245) * t159;
t233 = -t149 - t252;
t246 = qJD(2) * t161;
t296 = -pkin(6) - t169;
t164 = t171 * pkin(3);
t156 = t164 + pkin(2);
t297 = pkin(2) - t156;
t298 = pkin(6) * t159;
t314 = t159 * t297;
t153 = t161 * pkin(6);
t79 = t161 * t169 + t153 - t314;
t237 = t159 * ((-t161 * t297 - t298) * qJD(2) + t233) + t161 * (-t224 + t148 + (t161 * t296 + t314) * qJD(2)) + t79 * t246;
t247 = qJD(2) * t159;
t290 = rSges(5,2) * t160;
t253 = rSges(5,3) * t246 + t247 * t290;
t151 = t159 * rSges(5,3);
t311 = -t161 * t290 + t151;
t135 = t161 * t156;
t80 = -pkin(2) * t161 + t159 * t296 + t135;
t293 = rSges(5,1) * t162;
t220 = -t290 + t293;
t95 = -rSges(5,3) * t161 + t159 * t220;
t96 = t161 * t293 + t311;
t4 = (qJD(2) * t95 - t161 * t185 + t253) * t161 + (-t159 * t185 + (-t80 - t96 + t311) * qJD(2)) * t159 + t237;
t324 = t306 * t4;
t323 = t328 * t159 + t333 * t161;
t320 = t332 * t159 - t330;
t319 = -t333 * qJD(2) + t329 * t161;
t318 = -t329 * t159 - t326;
t291 = rSges(4,2) * t170;
t294 = rSges(4,1) * t171;
t221 = -t291 + t294;
t285 = t161 * rSges(4,3);
t105 = t159 * t221 - t285;
t236 = t161 * t291;
t152 = t159 * rSges(4,3);
t251 = t161 * t294 + t152;
t106 = -t236 + t251;
t143 = t170 * rSges(4,1) + rSges(4,2) * t171;
t186 = qJD(3) * t143;
t244 = qJD(2) * t170;
t232 = t159 * t244;
t173 = rSges(4,2) * t232 + rSges(4,3) * t246 - t161 * t186;
t13 = (qJD(2) * t105 + t173) * t161 + (-t159 * t186 + (-t106 - t236 + t152) * qJD(2)) * t159;
t307 = 2 * m(4);
t317 = t13 * t307;
t196 = Icges(6,5) * t155 - Icges(6,6) * t154;
t81 = -Icges(6,3) * t161 + t159 * t196;
t316 = qJD(2) * t81;
t288 = rSges(6,2) * t155;
t120 = rSges(6,1) * t154 + t288;
t167 = qJD(3) + qJD(5);
t187 = t120 * t167;
t276 = Icges(6,4) * t155;
t199 = -Icges(6,2) * t154 + t276;
t84 = Icges(6,6) * t159 + t161 * t199;
t277 = Icges(6,4) * t154;
t204 = Icges(6,1) * t155 - t277;
t86 = Icges(6,5) * t159 + t161 * t204;
t217 = t154 * t84 - t155 * t86;
t315 = t159 * t217;
t83 = -Icges(6,6) * t161 + t159 * t199;
t85 = -Icges(6,5) * t161 + t159 * t204;
t218 = t154 * t83 - t155 * t85;
t313 = t161 * t218;
t165 = -pkin(7) + t169;
t249 = t165 - t169;
t312 = t249 * t161;
t118 = Icges(6,2) * t155 + t277;
t119 = Icges(6,1) * t154 + t276;
t193 = t118 * t154 - t119 * t155;
t308 = qJD(2) * t193 + t196 * t167;
t305 = 2 * m(6);
t304 = t159 / 0.2e1;
t303 = -t161 / 0.2e1;
t302 = m(4) * t143;
t301 = pkin(3) * t170;
t300 = pkin(4) * t160;
t299 = pkin(4) * t162;
t295 = t159 * t79 + t161 * t80;
t87 = -rSges(6,3) * t161 + t159 * t219;
t36 = t159 * t87 + t161 * t88;
t287 = t160 * t93;
t286 = t160 * t94;
t284 = t162 * t91;
t283 = t162 * t92;
t282 = rSges(6,3) - t165;
t82 = Icges(6,3) * t159 + t161 * t196;
t266 = qJD(2) * t82;
t264 = t101 * t171;
t263 = t102 * t171;
t262 = t118 * t167;
t261 = t119 * t167;
t260 = t154 * t167;
t259 = t155 * t167;
t258 = t159 * t165;
t257 = t161 * t167;
t256 = t170 * t103;
t255 = t170 * t104;
t254 = rSges(6,3) * t246 + t247 * t289;
t250 = t157 + t158;
t243 = qJD(3) * t159;
t242 = qJD(3) * t160;
t241 = qJD(3) * t162;
t239 = qJD(3) * t171;
t238 = t159 * (t88 * qJD(2) - t159 * t187) + t161 * (-t257 * t288 + (-t154 * t257 - t155 * t247) * rSges(6,1) + t254) + t87 * t246;
t229 = -t126 - t301;
t41 = -qJD(2) * t83 - t161 * t262;
t228 = t167 * t86 + t41;
t42 = qJD(2) * t84 - t159 * t262;
t227 = t167 * t85 + t42;
t43 = -qJD(2) * t85 - t161 * t261;
t226 = -t167 * t84 + t43;
t44 = qJD(2) * t86 - t159 * t261;
t225 = t167 * t83 - t44;
t14 = -t159 * t218 - t161 * t81;
t15 = -t161 * t82 - t315;
t16 = t159 * t81 - t313;
t17 = t159 * t82 - t161 * t217;
t117 = Icges(6,5) * t154 + Icges(6,6) * t155;
t181 = t167 * t117;
t39 = -t161 * t181 - t316;
t40 = -t159 * t181 + t266;
t223 = -t161 * ((t161 * t40 + (t15 + t313) * qJD(2)) * t161 + (t14 * qJD(2) + (-t154 * t41 + t155 * t43 - t259 * t84 - t260 * t86 + t266) * t159 + (-t39 + t225 * t155 + t227 * t154 + (-t217 - t81) * qJD(2)) * t161) * t159) + t159 * ((t159 * t39 + (t16 + t315) * qJD(2)) * t159 + (t17 * qJD(2) + (t154 * t42 - t155 * t44 + t259 * t83 + t260 * t85 - t316) * t161 + (-t40 + t226 * t155 - t228 * t154 + (-t218 + t82) * qJD(2)) * t159) * t161) + (-t14 * t161 + t15 * t159) * t247 + (t159 * t17 - t16 * t161) * t246;
t222 = -t300 - t301;
t132 = t156 + t299;
t189 = -t132 - t219;
t57 = t159 * t189 + t161 * t282;
t112 = t161 * t132;
t58 = t112 + t88 - t258;
t216 = t159 * t58 + t161 * t57;
t191 = -t120 + t222;
t69 = t191 * t159;
t70 = t191 * t161;
t215 = t159 * t69 + t161 * t70;
t207 = Icges(4,1) * t170 + t280;
t205 = Icges(5,1) * t160 + t278;
t202 = Icges(4,2) * t171 + t281;
t200 = Icges(5,2) * t162 + t279;
t192 = -pkin(2) - t221;
t190 = -t156 - t220;
t108 = t199 * t167;
t109 = t204 * t167;
t172 = qJD(2) * t117 + (t109 - t262) * t155 + (-t108 - t261) * t154;
t184 = (t154 * t226 + t155 * t228 + t159 * t308 + t172 * t161) * t304 + (-t154 * t225 + t155 * t227 + t172 * t159 - t161 * t308) * t303 + (-t117 * t161 + t154 * t85 + t155 * t83 - t159 * t193) * t247 / 0.2e1 + (t117 * t159 + t154 * t86 + t155 * t84 - t161 * t193) * t246 / 0.2e1;
t180 = qJD(3) * t207;
t179 = qJD(3) * t205;
t178 = qJD(3) * t202;
t177 = qJD(3) * t200;
t110 = t219 * t167;
t174 = -t110 + (-t164 - t299) * qJD(3);
t138 = pkin(3) * t232;
t133 = t221 * qJD(3);
t124 = t222 * qJD(3);
t116 = t220 * qJD(3);
t111 = t161 * t124;
t98 = t229 * t161;
t97 = t229 * t159;
t76 = t298 + (pkin(2) - t291) * t161 + t251;
t75 = t159 * t192 + t153 + t285;
t68 = -t159 * t249 + t112 - t135;
t67 = t312 + t159 * (t132 - t156);
t66 = -t159 * t169 + t135 + t96;
t65 = (rSges(5,3) - t169) * t161 + t190 * t159;
t48 = -t126 * t246 - t159 * t116 + (-t159 * t239 - t161 * t244) * pkin(3);
t47 = t126 * t247 + t138 + (-pkin(3) * t239 - t116) * t161;
t46 = t143 * t243 + ((-rSges(4,3) - pkin(6)) * t159 + t192 * t161) * qJD(2);
t45 = (t153 + (-pkin(2) - t294) * t159) * qJD(2) + t173;
t31 = qJD(2) * t70 + t159 * t174;
t30 = t138 + (t120 + t300) * t247 + t174 * t161;
t29 = t126 * t243 + (t161 * t190 - t151) * qJD(2) - t233;
t28 = t148 + (-t156 - t293) * t247 + (qJD(3) * t229 - t245) * t161 + t253;
t19 = t149 + (-t124 + t187) * t159 + (-t159 * t282 + t161 * t189) * qJD(2);
t18 = t111 + t148 - t161 * t187 + (-t161 * t165 + (-t132 - t292) * t159) * qJD(2) + t254;
t12 = -t247 * t88 + t238;
t7 = t159 * t67 + t161 * t68 + t295 + t36;
t3 = t159 * (t124 * t159 + t252) + t161 * (t111 + t224) + ((t67 - t312) * t161 + (-t68 - t80 - t88 - t258) * t159) * qJD(2) + t237 + t238;
t1 = [0; 0; t119 * t259 + t154 * t109 - t118 * t260 + t155 * t108 + (t18 * t58 + t19 * t57) * t305 + (t28 * t66 + t29 * t65) * t306 + (t45 * t76 + t46 * t75) * t307 + (t206 - t200) * t242 + (t205 + t201) * t241 + (t208 - t202) * t240 + (t207 + t203) * t239; m(4) * t13 + m(5) * t4 + m(6) * t3; m(4) * ((-t159 * t45 - t161 * t46) * t143 + (-t159 * t76 - t161 * t75) * t133) + ((-t76 * t302 + t286 / 0.2e1 + t283 / 0.2e1 + t263 / 0.2e1 + t255 / 0.2e1) * t161 + (t75 * t302 + t287 / 0.2e1 + t284 / 0.2e1 + t264 / 0.2e1 + t256 / 0.2e1) * t159) * qJD(2) + m(5) * (t28 * t97 + t29 * t98 + t47 * t65 + t48 * t66) + m(6) * (t18 * t69 + t19 * t70 + t30 * t57 + t31 * t58) + t184 + (-qJD(3) * t310 + t160 * (-qJD(2) * t93 - t161 * t179) + t162 * (-qJD(2) * t91 - t161 * t177) + t170 * (-qJD(2) * t103 - t161 * t180) + t171 * (-qJD(2) * t101 - t161 * t178)) * t304 + (-qJD(3) * t328 + t160 * (qJD(2) * t94 - t159 * t179) + t162 * (qJD(2) * t92 - t159 * t177) + t170 * (qJD(2) * t104 - t159 * t180) + t171 * (qJD(2) * t102 - t159 * t178)) * t303 + t334 * qJD(3) * (t158 / 0.2e1 + t157 / 0.2e1); (t7 * t3 + t30 * t70 + t31 * t69) * t305 + (t295 * t4 + t98 * t47 + t97 * t48) * t306 + t250 * t143 * t133 * t307 + t223 + (t106 * t317 + t318 * t158 + t323 * t247 + t96 * t324 + (-t328 * t161 - t325) * t246) * t161 + (t95 * t324 + t105 * t317 + t319 * t157 + t320 * t246 + ((t102 * t239 + t104 * t240 + t241 * t92 + t94 * t242 + t318 - t326) * t159 + (t101 * t239 + t103 * t240 + t241 * t91 + t242 * t93 + t319) * t161 + ((-t283 - t286 - t255 - t263) * t159 + (-t284 - t287 - t256 - t264) * t161) * qJD(3) + ((-t328 + t332) * t159 + t330 + t320 + t323) * qJD(2)) * t161 + (t325 + t331) * t247) * t159; 0; m(6) * (qJD(2) * t216 + t159 * t19 - t161 * t18) + m(5) * (t159 * t29 - t161 * t28 + (t159 * t66 + t161 * t65) * qJD(2)); m(6) * (qJD(2) * t215 + t159 * t30 - t161 * t31) + m(5) * (t159 * t47 - t161 * t48 + (t159 * t97 + t161 * t98) * qJD(2)); 0; m(6) * t12; m(6) * (-t216 * t110 + (-t159 * t18 - t161 * t19 + (t159 * t57 - t161 * t58) * qJD(2)) * t120) + t184; m(6) * (t12 * t7 + t36 * t3 - t215 * t110 + (-t159 * t31 - t161 * t30 + (t159 * t70 - t161 * t69) * qJD(2)) * t120) + t223; 0; (t110 * t120 * t250 + t36 * t12) * t305 + t223;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
