% Calculate time derivative of joint inertia matrix for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:27
% EndTime: 2019-12-05 15:53:45
% DurationCPUTime: 7.15s
% Computational Cost: add. (18997->534), mult. (21402->850), div. (0->0), fcn. (20313->10), ass. (0->275)
t227 = sin(pkin(8));
t222 = t227 ^ 2;
t229 = cos(pkin(8));
t223 = t229 ^ 2;
t361 = t222 + t223;
t365 = 0.1e1 - t361;
t231 = sin(qJ(2));
t232 = cos(qJ(2));
t359 = qJD(2) * t361;
t364 = m(3) * (t231 * rSges(3,1) + rSges(3,2) * t232) * t359;
t226 = sin(pkin(9));
t228 = cos(pkin(9));
t281 = rSges(4,1) * t228 - rSges(4,2) * t226;
t360 = t232 * rSges(4,3) - t231 * t281;
t224 = pkin(9) + qJ(4);
t218 = cos(t224);
t307 = pkin(4) * t218;
t154 = -pkin(7) * t232 + t231 * t307;
t344 = t228 * pkin(3);
t358 = pkin(6) * t232 - t231 * t344;
t355 = 2 * m(5);
t354 = 2 * m(6);
t353 = 0.2e1 * qJD(2);
t352 = m(4) / 0.2e1;
t351 = m(5) / 0.2e1;
t350 = m(6) / 0.2e1;
t349 = t227 / 0.2e1;
t348 = -t229 / 0.2e1;
t347 = t229 / 0.2e1;
t346 = -t232 / 0.2e1;
t219 = qJ(5) + t224;
t214 = sin(t219);
t215 = cos(t219);
t324 = t227 * t232;
t183 = -t214 * t324 - t229 * t215;
t184 = -t229 * t214 + t215 * t324;
t325 = t227 * t231;
t116 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t325;
t118 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t325;
t120 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t325;
t321 = t229 * t232;
t185 = -t214 * t321 + t227 * t215;
t186 = t227 * t214 + t215 * t321;
t322 = t229 * t231;
t51 = t116 * t322 + t118 * t185 + t120 * t186;
t343 = t227 * t51;
t217 = sin(t224);
t191 = -t217 * t324 - t229 * t218;
t192 = -t229 * t217 + t218 * t324;
t125 = Icges(5,5) * t192 + Icges(5,6) * t191 + Icges(5,3) * t325;
t127 = Icges(5,4) * t192 + Icges(5,2) * t191 + Icges(5,6) * t325;
t129 = Icges(5,1) * t192 + Icges(5,4) * t191 + Icges(5,5) * t325;
t193 = -t217 * t321 + t227 * t218;
t194 = t227 * t217 + t218 * t321;
t57 = t125 * t322 + t127 * t193 + t129 * t194;
t342 = t227 * t57;
t117 = Icges(6,5) * t186 + Icges(6,6) * t185 + Icges(6,3) * t322;
t119 = Icges(6,4) * t186 + Icges(6,2) * t185 + Icges(6,6) * t322;
t121 = Icges(6,1) * t186 + Icges(6,4) * t185 + Icges(6,5) * t322;
t50 = t117 * t325 + t119 * t183 + t121 * t184;
t341 = t229 * t50;
t126 = Icges(5,5) * t194 + Icges(5,6) * t193 + Icges(5,3) * t322;
t128 = Icges(5,4) * t194 + Icges(5,2) * t193 + Icges(5,6) * t322;
t130 = Icges(5,1) * t194 + Icges(5,4) * t193 + Icges(5,5) * t322;
t56 = t126 * t325 + t128 * t191 + t130 * t192;
t340 = t229 * t56;
t303 = qJD(4) * t217;
t300 = pkin(4) * t303;
t233 = -qJD(2) * t154 - t232 * t300;
t302 = qJD(4) * t218;
t299 = pkin(4) * t302;
t225 = qJD(4) + qJD(5);
t305 = qJD(2) * t231;
t291 = t229 * t305;
t252 = -t225 * t227 + t291;
t296 = t225 * t321;
t152 = t214 * t252 - t215 * t296;
t153 = -t214 * t296 - t215 * t252;
t304 = qJD(2) * t232;
t290 = t229 * t304;
t87 = t153 * rSges(6,1) + t152 * rSges(6,2) + rSges(6,3) * t290;
t338 = t227 * t299 + t229 * t233 + t87;
t122 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t325;
t293 = t227 * t305;
t251 = t225 * t229 + t293;
t297 = t225 * t324;
t150 = t214 * t251 - t215 * t297;
t151 = -t214 * t297 - t215 * t251;
t292 = t227 * t304;
t86 = t151 * rSges(6,1) + t150 * rSges(6,2) + rSges(6,3) * t292;
t337 = t122 * t290 + t86 * t322;
t334 = Icges(5,4) * t217;
t333 = Icges(5,4) * t218;
t332 = Icges(6,4) * t214;
t331 = Icges(6,4) * t215;
t268 = Icges(6,5) * t215 - Icges(6,6) * t214;
t173 = -Icges(6,3) * t232 + t231 * t268;
t330 = t173 * t232;
t269 = Icges(5,5) * t218 - Icges(5,6) * t217;
t178 = -Icges(5,3) * t232 + t231 * t269;
t329 = t178 * t232;
t327 = t225 * t231;
t326 = t227 * t226;
t323 = t229 * t226;
t320 = t232 * ((-Icges(6,5) * t214 - Icges(6,6) * t215) * t327 + (Icges(6,3) * t231 + t232 * t268) * qJD(2));
t301 = qJD(4) * t231;
t319 = t232 * ((-Icges(5,5) * t217 - Icges(5,6) * t218) * t301 + (Icges(5,3) * t231 + t232 * t269) * qJD(2));
t240 = pkin(7) * t231 + t232 * t307;
t288 = pkin(4) * t217;
t103 = t227 * t240 - t229 * t288;
t317 = t103 + t122;
t104 = t227 * t288 + t229 * t240;
t123 = rSges(6,1) * t186 + rSges(6,2) * t185 + rSges(6,3) * t322;
t316 = t104 + t123;
t279 = rSges(6,1) * t215 - rSges(6,2) * t214;
t115 = (-rSges(6,1) * t214 - rSges(6,2) * t215) * t327 + (rSges(6,3) * t231 + t232 * t279) * qJD(2);
t124 = qJD(2) * t240 - t231 * t300;
t315 = -t115 - t124;
t176 = -rSges(6,3) * t232 + t231 * t279;
t314 = -t154 - t176;
t88 = t232 * t122 + t176 * t325;
t278 = pkin(2) * t232 + qJ(3) * t231;
t204 = qJD(2) * t278 - qJD(3) * t232;
t241 = pkin(6) * t231 + t232 * t344;
t313 = -t241 * qJD(2) - t204;
t212 = t231 * pkin(2) - qJ(3) * t232;
t312 = t361 * (-qJD(2) * t212 + qJD(3) * t231);
t311 = -t212 + t358;
t310 = -(rSges(4,3) * t231 + t232 * t281) * qJD(2) - t204;
t309 = -t212 + t360;
t308 = t361 * t278;
t298 = t115 * t325 + t176 * t292 + t232 * t86;
t280 = rSges(5,1) * t218 - rSges(5,2) * t217;
t136 = (-rSges(5,1) * t217 - rSges(5,2) * t218) * t301 + (rSges(5,3) * t231 + t232 * t280) * qJD(2);
t295 = -t136 + t313;
t182 = -rSges(5,3) * t232 + t231 * t280;
t294 = -t182 + t311;
t287 = t316 * t232;
t286 = t313 + t315;
t285 = t227 * (-pkin(3) * t323 + t227 * t241) + t229 * (pkin(3) * t326 + t229 * t241) + t308;
t284 = t358 * t359 + t312;
t283 = t311 + t314;
t267 = -t118 * t214 + t120 * t215;
t59 = -t116 * t232 + t231 * t267;
t266 = -t119 * t214 + t121 * t215;
t60 = -t117 * t232 + t231 * t266;
t277 = t59 * t227 + t60 * t229;
t265 = -t127 * t217 + t129 * t218;
t61 = -t125 * t232 + t231 * t265;
t264 = -t128 * t217 + t130 * t218;
t62 = -t126 * t232 + t231 * t264;
t276 = t61 * t227 + t62 * t229;
t274 = Icges(5,1) * t218 - t334;
t273 = Icges(6,1) * t215 - t332;
t271 = -Icges(5,2) * t217 + t333;
t270 = -Icges(6,2) * t214 + t331;
t131 = rSges(5,1) * t192 + rSges(5,2) * t191 + rSges(5,3) * t325;
t132 = rSges(5,1) * t194 + rSges(5,2) * t193 + rSges(5,3) * t322;
t263 = t131 * t229 - t132 * t227;
t174 = -Icges(6,6) * t232 + t231 * t270;
t175 = -Icges(6,5) * t232 + t231 * t273;
t261 = t174 * t214 - t175 * t215;
t179 = -Icges(5,6) * t232 + t231 * t271;
t180 = -Icges(5,5) * t232 + t231 * t274;
t260 = t179 * t217 - t180 * t218;
t80 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t292;
t256 = t116 * t304 + t231 * t80;
t81 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t290;
t255 = t117 * t304 + t231 * t81;
t157 = -qJD(4) * t192 + t217 * t293;
t158 = qJD(4) * t191 - t218 * t293;
t93 = Icges(5,5) * t158 + Icges(5,6) * t157 + Icges(5,3) * t292;
t254 = t125 * t304 + t231 * t93;
t159 = -qJD(4) * t194 + t217 * t291;
t160 = qJD(4) * t193 - t218 * t291;
t94 = Icges(5,5) * t160 + Icges(5,6) * t159 + Icges(5,3) * t290;
t253 = t126 * t304 + t231 * t94;
t113 = (-Icges(6,2) * t215 - t332) * t327 + (Icges(6,6) * t231 + t232 * t270) * qJD(2);
t114 = (-Icges(6,1) * t214 - t331) * t327 + (Icges(6,5) * t231 + t232 * t273) * qJD(2);
t82 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t292;
t84 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t292;
t18 = (qJD(2) * t267 - t80) * t232 + (qJD(2) * t116 + (-t118 * t225 + t84) * t215 + (-t120 * t225 - t82) * t214) * t231;
t83 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t290;
t85 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t290;
t19 = (qJD(2) * t266 - t81) * t232 + (qJD(2) * t117 + (-t119 * t225 + t85) * t215 + (-t121 * t225 - t83) * t214) * t231;
t49 = t116 * t325 + t118 * t183 + t120 * t184;
t20 = t150 * t118 + t151 * t120 + t183 * t82 + t184 * t84 + t227 * t256;
t21 = t150 * t119 + t151 * t121 + t183 * t83 + t184 * t85 + t227 * t255;
t66 = t173 * t325 + t174 * t183 + t175 * t184;
t5 = -(t183 * t113 + t184 * t114 + t150 * t174 + t151 * t175) * t232 + (t21 * t229 + (t20 - t320) * t227) * t231 + (t66 * t231 + (t341 + (t49 - t330) * t227) * t232) * qJD(2);
t52 = t117 * t322 + t119 * t185 + t121 * t186;
t22 = t152 * t118 + t153 * t120 + t185 * t82 + t186 * t84 + t229 * t256;
t23 = t152 * t119 + t153 * t121 + t185 * t83 + t186 * t85 + t229 * t255;
t67 = t173 * t322 + t174 * t185 + t175 * t186;
t6 = -(t185 * t113 + t186 * t114 + t152 * t174 + t153 * t175) * t232 + (t22 * t227 + (t23 - t320) * t229) * t231 + (t67 * t231 + (t343 + (t52 - t330) * t229) * t232) * qJD(2);
t79 = -t231 * t261 - t330;
t250 = -t232 * ((t320 + (t232 * t261 + t277) * qJD(2)) * t232 + (t19 * t229 + t18 * t227 - (qJD(2) * t173 + (-t174 * t225 + t114) * t215 + (-t175 * t225 - t113) * t214) * t232 + t79 * qJD(2)) * t231) + t6 * t322 + t5 * t325 + (-t232 * t66 + (t227 * t49 + t341) * t231) * t292 + (-t232 * t67 + (t229 * t52 + t343) * t231) * t290 + (t231 * t277 - t232 * t79) * t305;
t13 = -t20 * t229 + t21 * t227;
t14 = -t22 * t229 + t227 * t23;
t249 = t5 * t348 + t6 * t349 + (-t18 * t229 + t19 * t227) * t346 + t13 * t325 / 0.2e1 + t14 * t322 / 0.2e1 + (t227 * t60 - t229 * t59) * t305 / 0.2e1 + (t227 * (t227 * t50 - t229 * t49) + t229 * (t227 * t52 - t229 * t51)) * t304 / 0.2e1;
t242 = qJD(2) * (-Icges(3,5) * t231 - Icges(3,6) * t232);
t236 = qJD(2) * (Icges(4,5) * t232 + (-Icges(4,1) * t228 + Icges(4,4) * t226) * t231);
t235 = qJD(2) * (Icges(4,6) * t232 + (-Icges(4,4) * t228 + Icges(4,2) * t226) * t231);
t208 = t228 * t321 + t326;
t207 = -t226 * t321 + t227 * t228;
t206 = t228 * t324 - t323;
t205 = -t226 * t324 - t229 * t228;
t199 = t229 * t242;
t198 = t227 * t242;
t167 = t229 * t236;
t166 = t227 * t236;
t165 = t229 * t235;
t164 = t227 * t235;
t156 = t309 * t229;
t155 = t309 * t227;
t138 = t310 * t229;
t137 = t310 * t227;
t135 = (-Icges(5,1) * t217 - t333) * t301 + (Icges(5,5) * t231 + t232 * t274) * qJD(2);
t134 = (-Icges(5,2) * t218 - t334) * t301 + (Icges(5,6) * t231 + t232 * t271) * qJD(2);
t110 = t123 * t305;
t109 = t122 * t322;
t106 = t294 * t229;
t105 = t294 * t227;
t101 = t227 * t233 - t229 * t299;
t100 = t160 * rSges(5,1) + t159 * rSges(5,2) + rSges(5,3) * t290;
t99 = t158 * rSges(5,1) + t157 * rSges(5,2) + rSges(5,3) * t292;
t98 = Icges(5,1) * t160 + Icges(5,4) * t159 + Icges(5,5) * t290;
t97 = Icges(5,1) * t158 + Icges(5,4) * t157 + Icges(5,5) * t292;
t96 = Icges(5,4) * t160 + Icges(5,2) * t159 + Icges(5,6) * t290;
t95 = Icges(5,4) * t158 + Icges(5,2) * t157 + Icges(5,6) * t292;
t92 = -t132 * t232 - t182 * t322;
t91 = t131 * t232 + t182 * t325;
t90 = -t231 * t260 - t329;
t89 = -t123 * t232 - t176 * t322;
t76 = t295 * t229;
t75 = t295 * t227;
t74 = t283 * t229;
t73 = t283 * t227;
t72 = t263 * t231;
t71 = t178 * t322 + t179 * t193 + t180 * t194;
t70 = t178 * t325 + t179 * t191 + t180 * t192;
t69 = -t123 * t325 + t109;
t68 = t359 * t360 + t312;
t65 = t227 * (rSges(4,1) * t206 + rSges(4,2) * t205 + rSges(4,3) * t325) + t229 * (rSges(4,1) * t208 + rSges(4,2) * t207 + rSges(4,3) * t322) + t308;
t64 = t286 * t229;
t63 = t286 * t227;
t58 = t126 * t322 + t128 * t193 + t130 * t194;
t55 = t125 * t325 + t127 * t191 + t129 * t192;
t54 = t314 * t322 - t287;
t53 = t103 * t232 + t154 * t325 + t88;
t48 = -t136 * t322 - t100 * t232 + (t132 * t231 - t182 * t321) * qJD(2);
t47 = t136 * t325 + t232 * t99 + (-t131 * t231 + t182 * t324) * qJD(2);
t46 = t131 * t227 + t132 * t229 + t285;
t45 = -t232 * t87 + t110 + (-t115 * t231 - t176 * t304) * t229;
t44 = -t122 * t305 + t298;
t43 = t109 + (t103 * t229 - t227 * t316) * t231;
t42 = (-t100 * t227 + t229 * t99) * t231 + t263 * t304;
t41 = t100 * t229 + t227 * t99 + t284;
t40 = (-t123 * t304 - t231 * t87) * t227 + t337;
t38 = t227 * t317 + t229 * t316 + t285;
t34 = t338 * t229 + (t101 + t86) * t227 + t284;
t33 = t104 * t305 + t110 - t338 * t232 + (t231 * t315 + t304 * t314) * t229;
t32 = t124 * t325 + t101 * t232 + (t154 * t324 - t231 * t317) * qJD(2) + t298;
t31 = t159 * t128 + t160 * t130 + t193 * t96 + t194 * t98 + t229 * t253;
t30 = t159 * t127 + t160 * t129 + t193 * t95 + t194 * t97 + t229 * t254;
t29 = t157 * t128 + t158 * t130 + t191 * t96 + t192 * t98 + t227 * t253;
t28 = t157 * t127 + t158 * t129 + t191 * t95 + t192 * t97 + t227 * t254;
t27 = (qJD(2) * t264 - t94) * t232 + (qJD(2) * t126 - t217 * t96 + t218 * t98 + (-t128 * t218 - t130 * t217) * qJD(4)) * t231;
t26 = (qJD(2) * t265 - t93) * t232 + (qJD(2) * t125 - t217 * t95 + t218 * t97 + (-t127 * t218 - t129 * t217) * qJD(4)) * t231;
t17 = (t101 * t231 + t103 * t304) * t229 + (-qJD(2) * t287 - t231 * t338) * t227 + t337;
t16 = t227 * t31 - t229 * t30;
t15 = t227 * t29 - t229 * t28;
t9 = -(t193 * t134 + t194 * t135 + t159 * t179 + t160 * t180) * t232 + (t30 * t227 + (t31 - t319) * t229) * t231 + (t71 * t231 + (t342 + (t58 - t329) * t229) * t232) * qJD(2);
t8 = -(t191 * t134 + t192 * t135 + t157 * t179 + t158 * t180) * t232 + (t29 * t229 + (t28 - t319) * t227) * t231 + (t70 * t231 + (t340 + (t55 - t329) * t227) * t232) * qJD(2);
t1 = [0; m(4) * t68 + m(5) * t41 + m(6) * t34 - t364; (t34 * t38 + t63 * t73 + t64 * t74) * t354 + (t105 * t75 + t106 * t76 + t41 * t46) * t355 + 0.2e1 * m(4) * (t137 * t155 + t138 * t156 + t65 * t68) + (-t223 * t198 - t13 - t15 + (t205 * t164 + t206 * t166) * t229) * t229 + (t14 + t16 + t222 * t199 + (t207 * t165 + t208 * t167) * t227 + (-t207 * t164 - t205 * t165 - t208 * t166 - t206 * t167 - t227 * t198 + t229 * t199) * t229) * t227 + 0.2e1 * t365 * (rSges(3,1) * t232 - rSges(3,2) * t231) * t364; (m(4) + m(5) + m(6)) * t305; m(6) * (-t232 * t34 + t322 * t64 + t325 * t63) + m(5) * (-t232 * t41 + t322 * t76 + t325 * t75) + m(4) * (t137 * t325 + t138 * t322 - t232 * t68) + ((t231 * t38 + t321 * t74 + t324 * t73) * t350 + (t105 * t324 + t106 * t321 + t231 * t46) * t351 + (t155 * t324 + t156 * t321 + t231 * t65) * t352) * t353; -0.4e1 * (t352 + t351 + t350) * t365 * t231 * t304; m(5) * t42 + m(6) * t17; (t227 * t27 - t229 * t26) * t346 + t9 * t349 + t8 * t348 + (t15 * t349 + t16 * t347) * t231 + m(6) * (t17 * t38 + t32 * t74 + t33 * t73 + t34 * t43 + t53 * t64 + t54 * t63) + m(5) * (t105 * t48 + t106 * t47 + t41 * t72 + t42 * t46 + t75 * t92 + t76 * t91) + (t231 * (t227 * t62 - t229 * t61) / 0.2e1 + ((t227 * t56 - t229 * t55) * t349 + (t227 * t58 - t229 * t57) * t347) * t232) * qJD(2) + t249; m(5) * (-t232 * t42 + t322 * t47 + t325 * t48) + m(6) * (-t17 * t232 + t32 * t322 + t325 * t33) + ((t231 * t72 + t321 * t91 + t324 * t92) * t351 + (t231 * t43 + t321 * t53 + t324 * t54) * t350) * t353; (-t232 * t70 + (t227 * t55 + t340) * t231) * t292 + t8 * t325 + (-t232 * t71 + (t229 * t58 + t342) * t231) * t290 + t9 * t322 + (t17 * t43 + t32 * t53 + t33 * t54) * t354 + (t231 * t276 - t232 * t90) * t305 - t232 * ((t319 + (t232 * t260 + t276) * qJD(2)) * t232 + (t27 * t229 + t26 * t227 - (qJD(2) * t178 - t134 * t217 + t135 * t218 - t179 * t302 - t180 * t303) * t232 + t90 * qJD(2)) * t231) + (t42 * t72 + t47 * t91 + t48 * t92) * t355 + t250; m(6) * t40; m(6) * (t34 * t69 + t38 * t40 + t44 * t74 + t45 * t73 + t63 * t89 + t64 * t88) + t249; m(6) * (-t232 * t40 + (t227 * t45 + t229 * t44) * t231 + (t231 * t69 + (t227 * t89 + t229 * t88) * t232) * qJD(2)); m(6) * (t17 * t69 + t32 * t88 + t33 * t89 + t40 * t43 + t44 * t53 + t45 * t54) + t250; (t40 * t69 + t44 * t88 + t45 * t89) * t354 + t250;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
