% Calculate time derivative of joint inertia matrix for
% S4RRPR7
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR7_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR7_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:58
% EndTime: 2019-12-31 17:06:10
% DurationCPUTime: 6.41s
% Computational Cost: add. (8801->529), mult. (13533->773), div. (0->0), fcn. (12719->8), ass. (0->279)
t185 = sin(qJ(1));
t188 = cos(qJ(1));
t184 = sin(qJ(2));
t187 = cos(qJ(2));
t308 = Icges(3,4) * t187;
t218 = -Icges(3,2) * t184 + t308;
t127 = Icges(3,6) * t185 + t188 * t218;
t309 = Icges(3,4) * t184;
t223 = Icges(3,1) * t187 - t309;
t129 = Icges(3,5) * t185 + t188 * t223;
t207 = t127 * t184 - t129 * t187;
t343 = t185 * t207;
t181 = qJ(2) + pkin(7);
t174 = sin(t181);
t175 = cos(t181);
t306 = Icges(4,4) * t175;
t216 = -Icges(4,2) * t174 + t306;
t119 = Icges(4,6) * t185 + t188 * t216;
t307 = Icges(4,4) * t174;
t221 = Icges(4,1) * t175 - t307;
t121 = Icges(4,5) * t185 + t188 * t221;
t209 = t119 * t174 - t121 * t175;
t342 = t185 * t209;
t171 = pkin(2) * t187 + pkin(1);
t318 = pkin(1) - t171;
t341 = t185 * t318;
t126 = -Icges(3,6) * t188 + t185 * t218;
t128 = -Icges(3,5) * t188 + t185 * t223;
t208 = t126 * t184 - t128 * t187;
t340 = t188 * t208;
t118 = -Icges(4,6) * t188 + t185 * t216;
t120 = -Icges(4,5) * t188 + t185 * t221;
t210 = t118 * t174 - t120 * t175;
t339 = t188 * t210;
t178 = t185 * rSges(4,3);
t285 = t174 * t188;
t338 = -rSges(4,2) * t285 + t178;
t183 = sin(qJ(4));
t186 = cos(qJ(4));
t304 = Icges(5,4) * t186;
t214 = -Icges(5,2) * t183 + t304;
t109 = -Icges(5,6) * t175 + t174 * t214;
t305 = Icges(5,4) * t183;
t219 = Icges(5,1) * t186 - t305;
t110 = -Icges(5,5) * t175 + t174 * t219;
t337 = -t109 * t186 - t110 * t183;
t212 = Icges(4,5) * t175 - Icges(4,6) * t174;
t116 = -Icges(4,3) * t188 + t185 * t212;
t213 = Icges(3,5) * t187 - Icges(3,6) * t184;
t124 = -Icges(3,3) * t188 + t185 * t213;
t242 = qJD(1) * t175 - qJD(4);
t264 = qJD(2) * t188;
t249 = t174 * t264;
t336 = t185 * t242 + t249;
t211 = Icges(5,5) * t186 - Icges(5,6) * t183;
t108 = -Icges(5,3) * t175 + t174 * t211;
t266 = qJD(2) * t186;
t270 = qJD(2) * t174;
t297 = t109 * t183;
t263 = qJD(4) * t174;
t67 = (-Icges(5,5) * t183 - Icges(5,6) * t186) * t263 + (Icges(5,3) * t174 + t175 * t211) * qJD(2);
t69 = (-Icges(5,1) * t183 - t304) * t263 + (Icges(5,5) * t174 + t175 * t219) * qJD(2);
t335 = t174 * t186 * t69 + t108 * t270 + (-qJD(2) * t297 + t110 * t266 - t67) * t175;
t334 = 2 * m(3);
t333 = 2 * m(4);
t332 = 2 * m(5);
t331 = t185 ^ 2;
t330 = t188 ^ 2;
t329 = -t175 / 0.2e1;
t328 = t185 / 0.2e1;
t326 = t188 / 0.2e1;
t325 = -rSges(5,3) - pkin(6);
t160 = rSges(3,1) * t184 + rSges(3,2) * t187;
t324 = m(3) * t160;
t323 = pkin(2) * t184;
t322 = pkin(3) * t174;
t321 = pkin(3) * t175;
t320 = t185 * pkin(5);
t180 = t188 * pkin(5);
t319 = qJD(1) / 0.2e1;
t182 = -qJ(3) - pkin(5);
t317 = -pkin(5) - t182;
t316 = rSges(3,1) * t187;
t315 = rSges(4,1) * t175;
t314 = rSges(3,2) * t184;
t313 = rSges(3,3) * t188;
t68 = (-Icges(5,2) * t186 - t305) * t263 + (Icges(5,6) * t174 + t175 * t214) * qJD(2);
t312 = t183 * t68;
t179 = t185 * rSges(3,3);
t239 = pkin(6) * t174 + t321;
t279 = t186 * t188;
t281 = t185 * t183;
t136 = -t175 * t281 - t279;
t280 = t185 * t186;
t282 = t183 * t188;
t137 = t175 * t280 - t282;
t236 = -t137 * rSges(5,1) - t136 * rSges(5,2);
t286 = t174 * t185;
t89 = rSges(5,3) * t286 - t236;
t311 = t239 * t185 + t89;
t284 = t175 * t188;
t138 = -t175 * t282 + t280;
t139 = t175 * t279 + t281;
t90 = t139 * rSges(5,1) + t138 * rSges(5,2) + rSges(5,3) * t285;
t310 = pkin(3) * t284 + pkin(6) * t285 + t90;
t294 = t118 * t175;
t293 = t119 * t175;
t292 = t120 * t174;
t291 = t121 * t174;
t290 = t126 * t187;
t289 = t127 * t187;
t288 = t128 * t184;
t287 = t129 * t184;
t283 = t182 * t188;
t114 = t180 + t283 - t341;
t165 = t188 * t171;
t115 = -pkin(1) * t188 + t185 * t317 + t165;
t278 = t185 * t114 + t188 * t115;
t235 = rSges(5,1) * t186 - rSges(5,2) * t183;
t111 = -rSges(5,3) * t175 + t174 * t235;
t149 = -pkin(6) * t175 + t322;
t277 = t111 + t149;
t272 = qJD(1) * t185;
t252 = t174 * t272;
t271 = qJD(1) * t188;
t276 = rSges(4,2) * t252 + rSges(4,3) * t271;
t275 = t188 * t316 + t179;
t117 = Icges(4,3) * t185 + t188 * t212;
t274 = qJD(1) * t117;
t125 = Icges(3,3) * t185 + t188 * t213;
t273 = qJD(1) * t125;
t269 = qJD(2) * t175;
t268 = qJD(2) * t184;
t267 = qJD(2) * t185;
t265 = qJD(2) * t187;
t176 = qJD(3) * t185;
t258 = pkin(2) * t268;
t253 = qJD(3) * t188 + t182 * t272 + t185 * t258;
t262 = t114 * t271 + t185 * ((-t188 * t318 - t320) * qJD(1) - t253) + t188 * (-t188 * t258 + t176 + (t188 * t317 + t341) * qJD(1));
t247 = t175 * t264;
t243 = -qJD(4) * t175 + qJD(1);
t205 = t243 * t188;
t70 = t336 * t183 + t186 * t205;
t71 = t183 * t205 - t336 * t186;
t261 = t71 * rSges(5,1) + t70 * rSges(5,2) + rSges(5,3) * t247;
t260 = t188 * t314;
t257 = pkin(2) * t265;
t85 = Icges(5,4) * t137 + Icges(5,2) * t136 + Icges(5,6) * t286;
t87 = Icges(5,1) * t137 + Icges(5,4) * t136 + Icges(5,5) * t286;
t232 = -t183 * t85 + t186 * t87;
t83 = Icges(5,5) * t137 + Icges(5,6) * t136 + Icges(5,3) * t286;
t29 = t174 * t232 - t175 * t83;
t41 = t108 * t286 + t109 * t136 + t110 * t137;
t256 = t29 / 0.2e1 + t41 / 0.2e1;
t86 = Icges(5,4) * t139 + Icges(5,2) * t138 + Icges(5,6) * t285;
t88 = Icges(5,1) * t139 + Icges(5,4) * t138 + Icges(5,5) * t285;
t231 = -t183 * t86 + t186 * t88;
t84 = Icges(5,5) * t139 + Icges(5,6) * t138 + Icges(5,3) * t285;
t30 = t174 * t231 - t175 * t84;
t42 = t108 * t285 + t138 * t109 + t139 * t110;
t255 = t30 / 0.2e1 + t42 / 0.2e1;
t251 = t184 * t272;
t248 = t175 * t267;
t246 = t269 / 0.2e1;
t148 = rSges(4,1) * t174 + rSges(4,2) * t175;
t245 = -t148 - t323;
t244 = -t185 * t182 + t165;
t241 = -t277 - t323;
t72 = t243 * t280 + (t174 * t267 - t188 * t242) * t183;
t73 = t242 * t279 + (-t174 * t266 + t183 * t243) * t185;
t240 = t73 * rSges(5,1) + t72 * rSges(5,2);
t238 = -t314 + t316;
t237 = -rSges(4,2) * t174 + t315;
t25 = t136 * t85 + t137 * t87 + t286 * t83;
t26 = t136 * t86 + t137 * t88 + t286 * t84;
t17 = t26 * t185 - t188 * t25;
t228 = t185 * t25 + t188 * t26;
t27 = t138 * t85 + t139 * t87 + t285 * t83;
t28 = t138 * t86 + t139 * t88 + t285 * t84;
t18 = t28 * t185 - t188 * t27;
t227 = t185 * t27 + t188 * t28;
t226 = t30 * t185 - t188 * t29;
t225 = t185 * t29 + t188 * t30;
t224 = -t185 * t90 + t188 * t89;
t222 = Icges(3,1) * t184 + t308;
t220 = Icges(4,1) * t174 + t306;
t217 = Icges(3,2) * t187 + t309;
t215 = Icges(4,2) * t175 + t307;
t74 = (-rSges(5,1) * t183 - rSges(5,2) * t186) * t263 + (rSges(5,3) * t174 + t175 * t235) * qJD(2);
t206 = -t239 * qJD(2) - t257 - t74;
t123 = rSges(4,1) * t284 + t338;
t204 = -pkin(1) - t238;
t66 = t241 * t188;
t203 = -t171 - t237;
t201 = qJD(2) * t160;
t200 = qJD(2) * t148;
t199 = qJD(2) * t222;
t198 = qJD(2) * t220;
t197 = qJD(2) * t217;
t196 = qJD(2) * t215;
t195 = qJD(2) * (-Icges(3,5) * t184 - Icges(3,6) * t187);
t194 = qJD(2) * (-Icges(4,5) * t174 - Icges(4,6) * t175);
t193 = t174 * t325 - t171 - t321;
t192 = t174 * t271 + t248;
t191 = t247 - t252;
t190 = rSges(3,2) * t251 + rSges(3,3) * t271 - t188 * t201;
t189 = t185 * t193 - t283;
t168 = pkin(2) * t251;
t156 = pkin(6) * t247;
t153 = t238 * qJD(2);
t143 = t237 * qJD(2);
t133 = -t260 + t275;
t132 = t185 * t238 - t313;
t122 = -rSges(4,3) * t188 + t185 * t237;
t113 = t245 * t188;
t112 = t245 * t185;
t107 = t320 + (pkin(1) - t314) * t188 + t275;
t106 = t185 * t204 + t180 + t313;
t100 = t123 + t244;
t99 = (rSges(4,3) - t182) * t188 + t203 * t185;
t94 = t185 * t195 + t273;
t93 = -qJD(1) * t124 + t188 * t195;
t78 = t185 * t194 + t274;
t77 = -qJD(1) * t116 + t188 * t194;
t76 = t160 * t267 + ((-rSges(3,3) - pkin(5)) * t185 + t204 * t188) * qJD(1);
t75 = (t180 + (-pkin(1) - t316) * t185) * qJD(1) + t190;
t65 = t241 * t185;
t62 = -t148 * t271 - t185 * t143 + (-t184 * t271 - t185 * t265) * pkin(2);
t61 = t148 * t272 + t168 + (-t143 - t257) * t188;
t59 = t185 * t125 - t207 * t188;
t58 = t185 * t124 - t340;
t57 = -t125 * t188 - t343;
t56 = -t124 * t188 - t185 * t208;
t55 = t148 * t267 + (t188 * t203 - t178) * qJD(1) + t253;
t54 = t176 + (-t171 - t315) * t272 + (-qJD(1) * t182 + qJD(2) * t245) * t188 + t276;
t53 = t244 + t310;
t52 = t189 + t236;
t51 = -t111 * t285 - t175 * t90;
t50 = t111 * t286 + t175 * t89;
t49 = t185 * t117 - t209 * t188;
t48 = t185 * t116 - t339;
t47 = -t117 * t188 - t342;
t46 = -t116 * t188 - t185 * t210;
t45 = -t108 * t175 + (t110 * t186 - t297) * t174;
t44 = t45 * t270;
t43 = t224 * t174;
t40 = rSges(5,3) * t192 + t240;
t39 = -rSges(5,3) * t252 + t261;
t38 = Icges(5,1) * t73 + Icges(5,4) * t72 + Icges(5,5) * t192;
t37 = Icges(5,1) * t71 + Icges(5,4) * t70 + Icges(5,5) * t191;
t36 = Icges(5,4) * t73 + Icges(5,2) * t72 + Icges(5,6) * t192;
t35 = Icges(5,4) * t71 + Icges(5,2) * t70 + Icges(5,6) * t191;
t34 = Icges(5,5) * t73 + Icges(5,6) * t72 + Icges(5,3) * t192;
t33 = Icges(5,5) * t71 + Icges(5,6) * t70 + Icges(5,3) * t191;
t32 = qJD(1) * t66 + t185 * t206;
t31 = t188 * t206 + t272 * t277 + t168;
t24 = (t175 * t325 + t322) * t267 + t193 * t271 - t240 + t253;
t23 = t156 + t176 + (-t322 - t323) * t264 + t189 * qJD(1) + t261;
t22 = t185 * t311 + t188 * t310 + t278;
t21 = (t111 * t267 + t40) * t175 + (-qJD(2) * t89 + t111 * t271 + t185 * t74) * t174;
t20 = (-t111 * t264 - t39) * t175 + (qJD(2) * t90 + t111 * t272 - t188 * t74) * t174;
t19 = (t337 * qJD(4) - t312) * t174 + t335;
t16 = t108 * t192 + t72 * t109 + t73 * t110 + t136 * t68 + t137 * t69 + t286 * t67;
t15 = t108 * t191 + t70 * t109 + t71 * t110 + t138 * t68 + t139 * t69 + t285 * t67;
t14 = t224 * t269 + (-t185 * t39 + t188 * t40 + (-t185 * t89 - t188 * t90) * qJD(1)) * t174;
t13 = t174 * t227 - t42 * t175;
t12 = t174 * t228 - t41 * t175;
t11 = (-pkin(3) * t249 + t156 + t39) * t188 + (-t149 * t267 + t40) * t185 + (t311 * t188 + (-t115 - t310) * t185) * qJD(1) + t262;
t10 = (qJD(2) * t231 - t33) * t175 + (qJD(2) * t84 - t183 * t35 + t186 * t37 + (-t183 * t88 - t186 * t86) * qJD(4)) * t174;
t9 = (qJD(2) * t232 - t34) * t175 + (qJD(2) * t83 - t183 * t36 + t186 * t38 + (-t183 * t87 - t186 * t85) * qJD(4)) * t174;
t8 = t84 * t248 + t136 * t35 + t137 * t37 + t72 * t86 + t73 * t88 + (t185 * t33 + t271 * t84) * t174;
t7 = t83 * t248 + t136 * t36 + t137 * t38 + t72 * t85 + t73 * t87 + (t185 * t34 + t271 * t83) * t174;
t6 = t84 * t247 + t138 * t35 + t139 * t37 + t70 * t86 + t71 * t88 + (t188 * t33 - t272 * t84) * t174;
t5 = t83 * t247 + t138 * t36 + t139 * t38 + t70 * t85 + t71 * t87 + (t188 * t34 - t272 * t83) * t174;
t4 = qJD(1) * t228 + t8 * t185 - t188 * t7;
t3 = qJD(1) * t227 + t6 * t185 - t188 * t5;
t2 = (qJD(2) * t228 - t16) * t175 + (-qJD(1) * t17 + qJD(2) * t41 + t185 * t7 + t188 * t8) * t174;
t1 = (qJD(2) * t227 - t15) * t175 + (-qJD(1) * t18 + qJD(2) * t42 + t185 * t5 + t188 * t6) * t174;
t60 = [(t106 * t76 + t107 * t75) * t334 + (t100 * t54 + t55 * t99) * t333 - t174 * t312 + (t23 * t53 + t24 * t52) * t332 + (t221 - t215) * t270 + (t220 + t216) * t269 + (t223 - t217) * t268 + (t222 + t218) * t265 + t337 * t263 + t335; m(5) * (t23 * t65 + t24 * t66 + t31 * t52 + t32 * t53) + m(4) * (t100 * t62 + t112 * t54 + t113 * t55 + t61 * t99) + m(3) * ((-t185 * t75 - t188 * t76) * t160 + (-t106 * t188 - t107 * t185) * t153) + ((t293 / 0.2e1 + t291 / 0.2e1 - t107 * t324 + t289 / 0.2e1 + t287 / 0.2e1 + t255) * t188 + (t294 / 0.2e1 + t292 / 0.2e1 + t290 / 0.2e1 + t288 / 0.2e1 + t106 * t324 + t256) * t185) * qJD(1) + (t213 + t212) * qJD(2) * (t331 / 0.2e1 + t330 / 0.2e1) + (t174 * (-t120 * qJD(1) - t188 * t198) + t175 * (-t118 * qJD(1) - t188 * t196) + t184 * (-t128 * qJD(1) - t188 * t199) + t187 * (-t126 * qJD(1) - t188 * t197) + t10 + t15 + (-t207 - t209) * qJD(2)) * t328 - (t174 * (qJD(1) * t121 - t185 * t198) + t175 * (qJD(1) * t119 - t185 * t196) + t184 * (qJD(1) * t129 - t185 * t199) + t187 * (qJD(1) * t127 - t185 * t197) + t16 + t9 + (-t208 - t210) * qJD(2)) * t188 / 0.2e1; (t11 * t22 + t31 * t66 + t32 * t65) * t332 - t188 * t4 + t185 * t3 + (t113 * t61 + t112 * t62 + (t185 * t122 + t123 * t188 + t278) * ((qJD(1) * t122 - t188 * t200 + t276) * t188 + (-t185 * t200 + (-t115 - t123 + t338) * qJD(1)) * t185 + t262)) * t333 + t185 * ((t185 * t77 + (t48 + t342) * qJD(1)) * t185 + (t49 * qJD(1) + (t118 * t269 + t120 * t270) * t188 + (-t78 + (-t291 - t293) * qJD(2) + (t117 - t210) * qJD(1)) * t185) * t188) + t185 * ((t185 * t93 + (t58 + t343) * qJD(1)) * t185 + (t59 * qJD(1) + (t126 * t265 + t128 * t268) * t188 + (-t94 + (-t287 - t289) * qJD(2) + (t125 - t208) * qJD(1)) * t185) * t188) - t188 * ((t188 * t78 + (t47 + t339) * qJD(1)) * t188 + (t46 * qJD(1) + (-t119 * t269 - t121 * t270 + t274) * t185 + (-t77 + (t292 + t294) * qJD(2) - t209 * qJD(1)) * t188) * t185) + ((t185 * t132 + t133 * t188) * ((qJD(1) * t132 + t190) * t188 + (-t185 * t201 + (-t133 - t260 + t179) * qJD(1)) * t185) + (t330 + t331) * t160 * t153) * t334 - t188 * ((t188 * t94 + (t57 + t340) * qJD(1)) * t188 + (t56 * qJD(1) + (-t127 * t265 - t129 * t268 + t273) * t185 + (-t93 + (t288 + t290) * qJD(2) - t207 * qJD(1)) * t188) * t185) + (t17 + (-t46 - t56) * t188 + (t47 + t57) * t185) * t272 + (t18 + (-t48 - t58) * t188 + (t49 + t59) * t185) * t271; m(4) * (t185 * t55 - t188 * t54 + (t100 * t185 + t188 * t99) * qJD(1)) + m(5) * (t185 * t24 - t188 * t23 + (t185 * t53 + t188 * t52) * qJD(1)); m(5) * (t185 * t31 - t188 * t32 + (t185 * t65 + t188 * t66) * qJD(1)) + m(4) * (t185 * t61 - t188 * t62 + (t112 * t185 + t113 * t188) * qJD(1)); 0; t44 + m(5) * (t20 * t53 + t21 * t52 + t23 * t51 + t24 * t50) + (-t19 + (t185 * t256 + t188 * t255) * qJD(2)) * t175 + ((t10 / 0.2e1 + t15 / 0.2e1) * t188 + (t9 / 0.2e1 + t16 / 0.2e1) * t185 + (-t185 * t255 + t188 * t256) * qJD(1)) * t174; m(5) * (t11 * t43 + t14 * t22 + t20 * t65 + t21 * t66 + t31 * t50 + t32 * t51) + (-t2 / 0.2e1 + t18 * t246 + (qJD(1) * t30 - t9) * t329 + t13 * t319) * t188 + (t12 * t319 + (qJD(1) * t29 + t10) * t329 + t1 / 0.2e1 + t17 * t246) * t185 + (t3 * t326 + qJD(2) * t226 / 0.2e1 + t4 * t328 + (-t185 * t18 / 0.2e1 + t17 * t326) * qJD(1)) * t174; m(5) * (t21 * t185 - t188 * t20 + (t185 * t51 + t188 * t50) * qJD(1)); (t14 * t43 + t20 * t51 + t21 * t50) * t332 + (t19 * t175 - t44 + (t185 * t12 + t188 * t13 - t175 * t225) * qJD(2)) * t175 + (t188 * t1 + t185 * t2 - t175 * (t10 * t188 + t185 * t9) + (t174 * t225 - t45 * t175) * qJD(2) + (t188 * t12 - t185 * t13 + t175 * t226) * qJD(1)) * t174;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t60(1), t60(2), t60(4), t60(7); t60(2), t60(3), t60(5), t60(8); t60(4), t60(5), t60(6), t60(9); t60(7), t60(8), t60(9), t60(10);];
Mq = res;
