% Calculate time derivative of joint inertia matrix for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR7_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:55
% EndTime: 2019-12-31 18:19:09
% DurationCPUTime: 8.03s
% Computational Cost: add. (14315->514), mult. (14079->732), div. (0->0), fcn. (13129->10), ass. (0->273)
t370 = Icges(4,3) + Icges(5,3);
t186 = qJ(3) + pkin(9);
t181 = sin(t186);
t183 = cos(t186);
t190 = sin(qJ(3));
t193 = cos(qJ(3));
t369 = Icges(4,5) * t193 + Icges(5,5) * t183 - Icges(4,6) * t190 - Icges(5,6) * t181;
t187 = qJ(1) + pkin(8);
t182 = sin(t187);
t184 = cos(t187);
t365 = t182 * t370 + t369 * t184;
t316 = Icges(5,4) * t183;
t224 = -Icges(5,2) * t181 + t316;
t115 = Icges(5,6) * t182 + t184 * t224;
t317 = Icges(5,4) * t181;
t229 = Icges(5,1) * t183 - t317;
t117 = Icges(5,5) * t182 + t184 * t229;
t318 = Icges(4,4) * t193;
t226 = -Icges(4,2) * t190 + t318;
t127 = Icges(4,6) * t182 + t184 * t226;
t319 = Icges(4,4) * t190;
t231 = Icges(4,1) * t193 - t319;
t130 = Icges(4,5) * t182 + t184 * t231;
t347 = t115 * t181 - t117 * t183 + t127 * t190 - t130 * t193;
t368 = t347 * t182;
t367 = t347 * t184;
t366 = t369 * t182 - t184 * t370;
t364 = (-Icges(4,5) * t190 - Icges(5,5) * t181 - Icges(4,6) * t193 - Icges(5,6) * t183) * qJD(3);
t114 = -Icges(5,6) * t184 + t182 * t224;
t116 = -Icges(5,5) * t184 + t182 * t229;
t126 = -Icges(4,6) * t184 + t182 * t226;
t129 = -Icges(4,5) * t184 + t182 * t231;
t363 = t114 * t181 - t116 * t183 + t126 * t190 - t129 * t193;
t362 = t365 * qJD(1);
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t314 = Icges(6,4) * t192;
t222 = -Icges(6,2) * t189 + t314;
t125 = -Icges(6,6) * t183 + t181 * t222;
t315 = Icges(6,4) * t189;
t227 = Icges(6,1) * t192 - t315;
t128 = -Icges(6,5) * t183 + t181 * t227;
t361 = -t125 * t192 - t128 * t189;
t360 = t366 * t182 - t368 + (-t363 - t365) * t184;
t219 = Icges(6,5) * t192 - Icges(6,6) * t189;
t122 = -Icges(6,3) * t183 + t181 * t219;
t274 = qJD(3) * t192;
t276 = qJD(3) * t189;
t280 = qJD(3) * t181;
t272 = qJD(5) * t181;
t93 = (-Icges(6,5) * t189 - Icges(6,6) * t192) * t272 + (Icges(6,3) * t181 + t183 * t219) * qJD(3);
t99 = (-Icges(6,1) * t189 - t314) * t272 + (Icges(6,5) * t181 + t183 * t227) * qJD(3);
t359 = t181 * t192 * t99 + t122 * t280 + (-t125 * t276 + t128 * t274 - t93) * t183;
t342 = t182 ^ 2;
t341 = t184 ^ 2;
t358 = t363 * t182 + t366 * t184;
t355 = t365 * t182 - t367;
t354 = -t366 * qJD(1) + t364 * t184;
t353 = -t364 * t182 - t362;
t180 = pkin(3) * t193 + pkin(2);
t163 = t184 * t180;
t188 = -qJ(4) - pkin(6);
t327 = -pkin(6) - t188;
t110 = -pkin(2) * t184 + t182 * t327 + t163;
t325 = rSges(5,1) * t183;
t246 = -rSges(5,2) * t181 + t325;
t118 = -rSges(5,3) * t184 + t182 * t246;
t293 = t183 * t184;
t177 = t182 * rSges(5,3);
t296 = t181 * t184;
t349 = -rSges(5,2) * t296 + t177;
t119 = rSges(5,1) * t293 + t349;
t150 = rSges(5,1) * t181 + rSges(5,2) * t183;
t209 = qJD(3) * t150;
t179 = t184 * pkin(6);
t292 = t184 * t188;
t328 = pkin(2) - t180;
t350 = t182 * t328;
t109 = t179 + t292 - t350;
t175 = qJD(4) * t182;
t275 = qJD(3) * t190;
t267 = pkin(3) * t275;
t283 = qJD(1) * t182;
t262 = qJD(4) * t184 + t182 * t267 + t188 * t283;
t282 = qJD(1) * t184;
t330 = pkin(6) * t182;
t271 = t109 * t282 + t182 * ((-t184 * t328 - t330) * qJD(1) - t262) + t184 * (-t184 * t267 + t175 + (t184 * t327 + t350) * qJD(1));
t261 = t181 * t283;
t287 = rSges(5,2) * t261 + rSges(5,3) * t282;
t19 = (qJD(1) * t118 - t184 * t209 + t287) * t184 + (-t182 * t209 + (-t110 - t119 + t349) * qJD(1)) * t182 + t271;
t344 = 2 * m(5);
t352 = t19 * t344;
t324 = rSges(4,2) * t190;
t326 = rSges(4,1) * t193;
t247 = -t324 + t326;
t323 = rSges(4,3) * t184;
t132 = t182 * t247 - t323;
t269 = t184 * t324;
t178 = t182 * rSges(4,3);
t286 = t184 * t326 + t178;
t133 = -t269 + t286;
t171 = rSges(4,1) * t190 + rSges(4,2) * t193;
t210 = qJD(3) * t171;
t281 = qJD(1) * t190;
t260 = t182 * t281;
t196 = rSges(4,2) * t260 + rSges(4,3) * t282 - t184 * t210;
t32 = (qJD(1) * t132 + t196) * t184 + (-t182 * t210 + (-t133 - t269 + t178) * qJD(1)) * t182;
t345 = 2 * m(4);
t351 = t32 * t345;
t334 = sin(qJ(1)) * pkin(1);
t348 = t179 - t334;
t343 = 2 * m(6);
t340 = t182 / 0.2e1;
t339 = -t183 / 0.2e1;
t337 = t184 / 0.2e1;
t336 = -rSges(6,3) - pkin(7);
t335 = m(4) * t171;
t333 = pkin(3) * t190;
t332 = pkin(4) * t181;
t331 = pkin(4) * t183;
t185 = cos(qJ(1)) * pkin(1);
t329 = qJD(1) / 0.2e1;
t96 = (-Icges(6,2) * t192 - t315) * t272 + (Icges(6,6) * t181 + t183 * t222) * qJD(3);
t322 = t189 * t96;
t250 = pkin(7) * t181 + t331;
t290 = t184 * t192;
t295 = t182 * t189;
t138 = -t183 * t295 - t290;
t291 = t184 * t189;
t294 = t182 * t192;
t139 = t183 * t294 - t291;
t245 = -rSges(6,1) * t139 - rSges(6,2) * t138;
t297 = t181 * t182;
t89 = rSges(6,3) * t297 - t245;
t321 = t250 * t182 + t89;
t140 = -t183 * t291 + t294;
t141 = t183 * t290 + t295;
t90 = t141 * rSges(6,1) + t140 * rSges(6,2) + rSges(6,3) * t296;
t320 = pkin(4) * t293 + pkin(7) * t296 + t90;
t307 = t114 * t183;
t306 = t115 * t183;
t305 = t116 * t181;
t304 = t117 * t181;
t302 = t126 * t193;
t301 = t127 * t193;
t299 = t129 * t190;
t298 = t130 * t190;
t289 = t182 * t109 + t184 * t110;
t244 = rSges(6,1) * t192 - rSges(6,2) * t189;
t131 = -rSges(6,3) * t183 + t181 * t244;
t151 = -pkin(7) * t183 + t332;
t288 = t131 + t151;
t279 = qJD(3) * t182;
t278 = qJD(3) * t183;
t277 = qJD(3) * t184;
t273 = qJD(3) * t193;
t257 = t183 * t277;
t254 = -qJD(5) * t183 + qJD(1);
t198 = t181 * t276 + t192 * t254;
t253 = qJD(1) * t183 - qJD(5);
t68 = t184 * t198 + t253 * t295;
t197 = -t181 * t274 + t189 * t254;
t69 = t184 * t197 - t253 * t294;
t270 = t69 * rSges(6,1) + t68 * rSges(6,2) + rSges(6,3) * t257;
t266 = pkin(3) * t273;
t85 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t297;
t87 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t297;
t235 = -t189 * t85 + t192 * t87;
t83 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t297;
t30 = t181 * t235 - t183 * t83;
t44 = t122 * t297 + t125 * t138 + t128 * t139;
t265 = t30 / 0.2e1 + t44 / 0.2e1;
t86 = Icges(6,4) * t141 + Icges(6,2) * t140 + Icges(6,6) * t296;
t88 = Icges(6,1) * t141 + Icges(6,4) * t140 + Icges(6,5) * t296;
t234 = -t189 * t86 + t192 * t88;
t84 = Icges(6,5) * t141 + Icges(6,6) * t140 + Icges(6,3) * t296;
t31 = t181 * t234 - t183 * t84;
t45 = t122 * t296 + t125 * t140 + t128 * t141;
t264 = t31 / 0.2e1 + t45 / 0.2e1;
t258 = t182 * t278;
t256 = t278 / 0.2e1;
t255 = -t150 - t333;
t252 = -t288 - t333;
t70 = t182 * t198 - t253 * t291;
t71 = t182 * t197 + t253 * t290;
t251 = t71 * rSges(6,1) + t70 * rSges(6,2);
t249 = -t182 * t188 + t163 + t185;
t121 = t255 * t184;
t243 = -t292 - t334;
t26 = t138 * t85 + t139 * t87 + t297 * t83;
t27 = t138 * t86 + t139 * t88 + t297 * t84;
t15 = t182 * t27 - t184 * t26;
t240 = t182 * t26 + t184 * t27;
t28 = t140 * t85 + t141 * t87 + t296 * t83;
t29 = t140 * t86 + t141 * t88 + t296 * t84;
t16 = t182 * t29 - t184 * t28;
t239 = t182 * t28 + t184 * t29;
t238 = t31 * t182 - t30 * t184;
t237 = t182 * t30 + t184 * t31;
t236 = -t182 * t90 + t184 * t89;
t230 = Icges(4,1) * t190 + t318;
t228 = Icges(5,1) * t181 + t316;
t225 = Icges(4,2) * t193 + t319;
t223 = Icges(5,2) * t183 + t317;
t102 = (-rSges(6,1) * t189 - rSges(6,2) * t192) * t272 + (rSges(6,3) * t181 + t183 * t244) * qJD(3);
t214 = -t250 * qJD(3) - t102 - t266;
t213 = -pkin(2) - t247;
t81 = t252 * t184;
t212 = -t180 - t246;
t208 = qJD(3) * t230;
t207 = qJD(3) * t228;
t206 = qJD(3) * t225;
t205 = qJD(3) * t223;
t202 = t181 * t336 - t180 - t331;
t200 = t181 * t282 + t258;
t199 = t257 - t261;
t195 = t182 * t202 + t243;
t166 = pkin(3) * t260;
t158 = t247 * qJD(3);
t154 = pkin(7) * t257;
t145 = t246 * qJD(3);
t120 = t255 * t182;
t104 = t330 + t185 + (pkin(2) - t324) * t184 + t286;
t103 = t182 * t213 + t323 + t348;
t92 = t119 + t249;
t91 = -t334 + (rSges(5,3) - t188) * t184 + t212 * t182;
t80 = t252 * t182;
t67 = -t150 * t282 - t145 * t182 + (-t182 * t273 - t184 * t281) * pkin(3);
t66 = t150 * t283 + t166 + (-t145 - t266) * t184;
t63 = t171 * t279 + (-t185 + (-rSges(4,3) - pkin(6)) * t182 + t213 * t184) * qJD(1);
t62 = ((-pkin(2) - t326) * t182 + t348) * qJD(1) + t196;
t57 = -t122 * t183 + (-t125 * t189 + t128 * t192) * t181;
t56 = -t131 * t296 - t183 * t90;
t55 = t131 * t297 + t183 * t89;
t54 = t150 * t279 + (t184 * t212 - t177 - t185) * qJD(1) + t262;
t53 = t175 + qJD(3) * t121 + ((-t180 - t325) * t182 + t243) * qJD(1) + t287;
t52 = t249 + t320;
t51 = t195 + t245;
t50 = t57 * t280;
t43 = t236 * t181;
t42 = qJD(1) * t81 + t182 * t214;
t41 = t184 * t214 + t283 * t288 + t166;
t40 = rSges(6,3) * t200 + t251;
t39 = -rSges(6,3) * t261 + t270;
t38 = Icges(6,1) * t71 + Icges(6,4) * t70 + Icges(6,5) * t200;
t37 = Icges(6,1) * t69 + Icges(6,4) * t68 + Icges(6,5) * t199;
t36 = Icges(6,4) * t71 + Icges(6,2) * t70 + Icges(6,6) * t200;
t35 = Icges(6,4) * t69 + Icges(6,2) * t68 + Icges(6,6) * t199;
t34 = Icges(6,5) * t71 + Icges(6,6) * t70 + Icges(6,3) * t200;
t33 = Icges(6,5) * t69 + Icges(6,6) * t68 + Icges(6,3) * t199;
t25 = (t183 * t336 + t332) * t279 + (t184 * t202 - t185) * qJD(1) - t251 + t262;
t24 = t154 + t175 + (-t332 - t333) * t277 + t195 * qJD(1) + t270;
t23 = t182 * t321 + t184 * t320 + t289;
t22 = (t361 * qJD(5) - t322) * t181 + t359;
t21 = (t131 * t279 + t40) * t183 + (-qJD(3) * t89 + t102 * t182 + t131 * t282) * t181;
t20 = (-t131 * t277 - t39) * t183 + (qJD(3) * t90 - t102 * t184 + t131 * t283) * t181;
t18 = t122 * t200 + t70 * t125 + t128 * t71 + t138 * t96 + t139 * t99 + t297 * t93;
t17 = t122 * t199 + t68 * t125 + t128 * t69 + t140 * t96 + t141 * t99 + t296 * t93;
t14 = t236 * t278 + (-t182 * t39 + t184 * t40 + (-t182 * t89 - t184 * t90) * qJD(1)) * t181;
t13 = t181 * t239 - t183 * t45;
t12 = t181 * t240 - t183 * t44;
t11 = (-t277 * t332 + t154 + t39) * t184 + (-t151 * t279 + t40) * t182 + (t321 * t184 + (-t110 - t320) * t182) * qJD(1) + t271;
t10 = (qJD(3) * t234 - t33) * t183 + (qJD(3) * t84 - t189 * t35 + t192 * t37 + (-t189 * t88 - t192 * t86) * qJD(5)) * t181;
t9 = (qJD(3) * t235 - t34) * t183 + (qJD(3) * t83 - t189 * t36 + t192 * t38 + (-t189 * t87 - t192 * t85) * qJD(5)) * t181;
t8 = t84 * t258 + t138 * t35 + t139 * t37 + t70 * t86 + t71 * t88 + (t182 * t33 + t282 * t84) * t181;
t7 = t83 * t258 + t138 * t36 + t139 * t38 + t70 * t85 + t71 * t87 + (t182 * t34 + t282 * t83) * t181;
t6 = t84 * t257 + t140 * t35 + t141 * t37 + t68 * t86 + t69 * t88 + (t184 * t33 - t283 * t84) * t181;
t5 = t83 * t257 + t140 * t36 + t141 * t38 + t68 * t85 + t69 * t87 + (t184 * t34 - t283 * t83) * t181;
t4 = qJD(1) * t240 + t182 * t8 - t184 * t7;
t3 = qJD(1) * t239 + t182 * t6 - t184 * t5;
t2 = (qJD(3) * t240 - t18) * t183 + (-qJD(1) * t15 + qJD(3) * t44 + t182 * t7 + t184 * t8) * t181;
t1 = (qJD(3) * t239 - t17) * t183 + (-qJD(1) * t16 + qJD(3) * t45 + t182 * t5 + t184 * t6) * t181;
t46 = [-t181 * t322 + (t24 * t52 + t25 * t51) * t343 + (t53 * t92 + t54 * t91) * t344 + (t103 * t63 + t104 * t62) * t345 + (t229 - t223) * t280 + (t228 + t224) * t278 + (t231 - t225) * t275 + (t230 + t226) * t273 + t361 * t272 + t359; 0; 0; m(4) * ((-t182 * t62 - t184 * t63) * t171 + (-t103 * t184 - t104 * t182) * t158) + m(5) * (t120 * t53 + t121 * t54 + t66 * t91 + t67 * t92) + m(6) * (t24 * t80 + t25 * t81 + t41 * t51 + t42 * t52) + ((t301 / 0.2e1 + t298 / 0.2e1 + t306 / 0.2e1 + t304 / 0.2e1 - t104 * t335 + t264) * t184 + (t302 / 0.2e1 + t299 / 0.2e1 + t103 * t335 + t307 / 0.2e1 + t305 / 0.2e1 + t265) * t182) * qJD(1) + t369 * qJD(3) * (t341 / 0.2e1 + t342 / 0.2e1) + (-t347 * qJD(3) + (-t129 * qJD(1) - t184 * t208) * t190 + t181 * (-t116 * qJD(1) - t184 * t207) + t183 * (-t114 * qJD(1) - t184 * t205) + t193 * (-t126 * qJD(1) - t184 * t206) + t10 + t17) * t340 - (-t363 * qJD(3) + (qJD(1) * t130 - t182 * t208) * t190 + t181 * (qJD(1) * t117 - t182 * t207) + t183 * (qJD(1) * t115 - t182 * t205) + t193 * (qJD(1) * t127 - t182 * t206) + t18 + t9) * t184 / 0.2e1; m(4) * t32 + m(5) * t19 + m(6) * t11; t15 * t283 + (t11 * t23 + t41 * t81 + t42 * t80) * t343 + t16 * t282 + (t120 * t67 + t121 * t66 + t289 * t19) * t344 + (t341 + t342) * t171 * t158 * t345 + (t119 * t352 + t133 * t351 + t358 * t283 + t353 * t341 - t4 + (-t184 * t363 - t360) * t282) * t184 + (t3 + t118 * t352 + t132 * t351 + t354 * t342 + t355 * t282 + ((t115 * t278 + t117 * t280 + t127 * t273 + t130 * t275 + t353 - t362) * t182 + (t114 * t278 + t116 * t280 + t126 * t273 + t129 * t275 + t354) * t184 + ((-t304 - t306 - t298 - t301) * t182 + (-t305 - t307 - t299 - t302) * t184) * qJD(3) + ((-t363 + t365) * t182 + t367 + t355 + t358) * qJD(1)) * t184 + (t360 + t368) * t283) * t182; m(6) * (t182 * t25 - t184 * t24 + (t182 * t52 + t184 * t51) * qJD(1)) + m(5) * (t182 * t54 - t184 * t53 + (t182 * t92 + t184 * t91) * qJD(1)); 0; m(6) * (t182 * t41 - t184 * t42 + (t182 * t80 + t184 * t81) * qJD(1)) + m(5) * (t182 * t66 - t184 * t67 + (t120 * t182 + t121 * t184) * qJD(1)); 0; t50 + m(6) * (t20 * t52 + t21 * t51 + t24 * t56 + t25 * t55) + (-t22 + (t182 * t265 + t184 * t264) * qJD(3)) * t183 + ((t10 / 0.2e1 + t17 / 0.2e1) * t184 + (t9 / 0.2e1 + t18 / 0.2e1) * t182 + (-t182 * t264 + t184 * t265) * qJD(1)) * t181; m(6) * t14; m(6) * (t11 * t43 + t14 * t23 + t20 * t80 + t21 * t81 + t41 * t55 + t42 * t56) + (-t2 / 0.2e1 + t13 * t329 + (qJD(1) * t31 - t9) * t339 + t16 * t256) * t184 + (t12 * t329 + t1 / 0.2e1 + t15 * t256 + (qJD(1) * t30 + t10) * t339) * t182 + (t4 * t340 + qJD(3) * t238 / 0.2e1 + t3 * t337 + (t15 * t337 - t182 * t16 / 0.2e1) * qJD(1)) * t181; m(6) * (t182 * t21 - t184 * t20 + (t182 * t56 + t184 * t55) * qJD(1)); (t14 * t43 + t20 * t56 + t21 * t55) * t343 + (t22 * t183 - t50 + (t182 * t12 + t184 * t13 - t183 * t237) * qJD(3)) * t183 + (t184 * t1 + t182 * t2 - t183 * (t10 * t184 + t182 * t9) + (t181 * t237 - t57 * t183) * qJD(3) + (t184 * t12 - t182 * t13 + t183 * t238) * qJD(1)) * t181;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t46(1), t46(2), t46(4), t46(7), t46(11); t46(2), t46(3), t46(5), t46(8), t46(12); t46(4), t46(5), t46(6), t46(9), t46(13); t46(7), t46(8), t46(9), t46(10), t46(14); t46(11), t46(12), t46(13), t46(14), t46(15);];
Mq = res;
