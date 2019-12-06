% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:58
% EndTime: 2019-12-05 17:31:22
% DurationCPUTime: 12.64s
% Computational Cost: add. (9940->533), mult. (26586->687), div. (0->0), fcn. (30448->10), ass. (0->237)
t237 = sin(pkin(9));
t240 = cos(pkin(9));
t241 = cos(pkin(7));
t239 = sin(pkin(7));
t340 = cos(pkin(8));
t296 = t239 * t340;
t198 = t237 * t296 + t241 * t240;
t164 = qJD(5) * t198 + qJD(1);
t243 = sin(qJ(1));
t295 = t243 * t340;
t238 = sin(pkin(8));
t245 = cos(qJ(1));
t329 = t245 * t238;
t201 = -t241 * t295 + t329;
t332 = t239 * t243;
t154 = t201 * t240 - t237 * t332;
t294 = t245 * t340;
t330 = t243 * t238;
t200 = t241 * t330 + t294;
t242 = sin(qJ(5));
t244 = cos(qJ(5));
t115 = -t154 * t242 - t200 * t244;
t116 = t154 * t244 - t200 * t242;
t153 = t201 * t237 + t240 * t332;
t56 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t153;
t339 = Icges(6,4) * t116;
t59 = Icges(6,2) * t115 + Icges(6,6) * t153 + t339;
t109 = Icges(6,4) * t115;
t62 = Icges(6,1) * t116 + Icges(6,5) * t153 + t109;
t14 = t115 * t59 + t116 * t62 + t153 * t56;
t203 = t241 * t294 + t330;
t331 = t239 * t245;
t156 = t203 * t240 + t237 * t331;
t202 = t241 * t329 - t295;
t118 = t156 * t244 + t202 * t242;
t119 = t156 * t242 - t202 * t244;
t155 = t203 * t237 - t240 * t331;
t57 = Icges(6,5) * t118 - Icges(6,6) * t119 + Icges(6,3) * t155;
t338 = Icges(6,4) * t118;
t61 = Icges(6,2) * t119 - Icges(6,6) * t155 - t338;
t110 = Icges(6,4) * t119;
t64 = -Icges(6,1) * t118 - Icges(6,5) * t155 + t110;
t15 = -t115 * t61 - t116 * t64 + t153 * t57;
t276 = t14 * t153 + t15 * t155;
t199 = -t241 * t237 + t240 * t296;
t333 = t238 * t239;
t151 = -t199 * t242 + t244 * t333;
t152 = t199 * t244 + t242 * t333;
t88 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t198;
t337 = Icges(6,4) * t152;
t89 = Icges(6,2) * t151 + Icges(6,6) * t198 + t337;
t142 = Icges(6,4) * t151;
t90 = Icges(6,1) * t152 + Icges(6,5) * t198 + t142;
t28 = t115 * t89 + t116 * t90 + t153 * t88;
t5 = qJD(5) * t276 + t28 * t164;
t313 = qJD(5) * t155;
t298 = t313 / 0.2e1;
t16 = t118 * t62 - t119 * t59 + t155 * t56;
t17 = -t118 * t64 + t119 * t61 + t57 * t155;
t275 = t153 * t16 + t155 * t17;
t21 = -t151 * t61 - t152 * t64 + t198 * t57;
t91 = rSges(6,1) * t152 + rSges(6,2) * t151 + rSges(6,3) * t198;
t383 = t155 * t91;
t66 = rSges(6,1) * t118 - t119 * rSges(6,2) + rSges(6,3) * t155;
t382 = t164 * t66;
t29 = t118 * t90 - t119 * t89 + t155 * t88;
t371 = t153 * t91;
t293 = -t203 * pkin(3) - t202 * qJ(4);
t336 = qJ(2) * t243;
t370 = t293 - t336;
t130 = t201 * pkin(3) - qJ(4) * t200;
t351 = pkin(2) * t241;
t278 = qJ(3) * t239 + t351;
t204 = t278 * t243;
t316 = qJD(3) * t245;
t223 = t239 * t316;
t234 = t245 * qJ(2);
t212 = -pkin(1) * t243 + t234;
t232 = qJD(2) * t243;
t324 = qJD(1) * t212 + t232;
t290 = -qJD(1) * t204 + t223 + t324;
t315 = qJD(4) * t202;
t362 = qJD(1) * t130 + t290 + t315;
t358 = -t203 * rSges(4,1) + t202 * rSges(4,2);
t361 = -rSges(4,3) * t331 + t358;
t106 = pkin(4) * t156 + pkin(6) * t155;
t93 = rSges(5,1) * t156 - rSges(5,2) * t155 + t202 * rSges(5,3);
t178 = t203 * qJD(1);
t318 = qJD(1) * t245;
t304 = t239 * t318;
t138 = -t178 * t237 + t240 * t304;
t139 = -t178 * t240 - t237 * t304;
t359 = t139 * pkin(4) + pkin(6) * t138;
t105 = t154 * pkin(4) + pkin(6) * t153;
t65 = t116 * rSges(6,1) + t115 * rSges(6,2) + t153 * rSges(6,3);
t357 = -qJD(1) * t105 - t164 * t65 - t362;
t92 = t154 * rSges(5,1) - t153 * rSges(5,2) - t200 * rSges(5,3);
t51 = qJD(1) * t92 + t362;
t177 = t202 * qJD(1);
t356 = t139 * rSges(5,1) - t138 * rSges(5,2) - t177 * rSges(5,3);
t355 = t241 ^ 2;
t98 = qJD(5) * t138 + qJDD(5) * t153;
t354 = t98 / 0.2e1;
t176 = t201 * qJD(1);
t319 = qJD(1) * t243;
t305 = t239 * t319;
t136 = t176 * t237 + t240 * t305;
t99 = qJD(5) * t136 + qJDD(5) * t155;
t353 = t99 / 0.2e1;
t352 = -m(5) - m(6);
t349 = g(2) * t245;
t20 = t151 * t59 + t152 * t62 + t198 * t56;
t348 = t20 * t98;
t347 = t21 * t99;
t163 = qJDD(5) * t198 + qJDD(1);
t140 = t151 * qJD(5);
t141 = t152 * qJD(5);
t94 = Icges(6,5) * t140 - Icges(6,6) * t141;
t95 = Icges(6,4) * t140 - Icges(6,2) * t141;
t96 = Icges(6,1) * t140 - Icges(6,4) * t141;
t19 = t140 * t90 - t141 * t89 + t151 * t95 + t152 * t96 + t198 * t94;
t32 = t151 * t89 + t152 * t90 + t198 * t88;
t345 = t32 * t163 + t19 * t164;
t344 = rSges(3,1) * t241;
t343 = -rSges(3,3) - qJ(2);
t342 = -Icges(6,2) * t152 + t142 + t90;
t341 = Icges(6,1) * t151 - t337 - t89;
t328 = -t178 * rSges(4,1) + t177 * rSges(4,2);
t327 = t201 * rSges(4,1) + t200 * rSges(4,2);
t205 = t278 * t245;
t214 = pkin(1) * t245 + t336;
t325 = -t205 - t214;
t226 = rSges(3,2) * t331;
t267 = -rSges(3,3) * t243 - t245 * t344;
t173 = -t226 - t267;
t323 = -t214 - t173;
t322 = rSges(3,2) * t332 + t245 * rSges(3,3);
t321 = -pkin(1) * t319 + t232;
t211 = qJD(1) * t214;
t233 = qJD(2) * t245;
t320 = t233 - t211;
t317 = qJD(3) * t243;
t314 = qJD(5) * t153;
t183 = t200 * qJD(4);
t312 = -m(4) + t352;
t311 = qJDD(3) * t239;
t73 = -qJD(5) * t116 - t139 * t242 - t177 * t244;
t74 = qJD(5) * t115 + t139 * t244 - t177 * t242;
t40 = t74 * rSges(6,1) + t73 * rSges(6,2) + t138 * rSges(6,3);
t308 = t293 + t325;
t306 = t325 + t361;
t162 = -qJ(3) * t305 - t319 * t351 + t223;
t303 = t239 * t317;
t302 = -pkin(1) - t344;
t300 = t314 / 0.2e1;
t299 = -t313 / 0.2e1;
t292 = -t93 + t308;
t206 = -qJDD(3) * t241 + qJDD(4) * t333;
t291 = qJDD(1) * t212 + qJDD(2) * t243 + (t233 + t320) * qJD(1);
t289 = -t106 + t308;
t175 = t200 * qJD(1);
t108 = pkin(3) * t176 - t175 * qJ(4) + t315;
t137 = t176 * t240 - t237 * t305;
t288 = -pkin(4) * t137 - pkin(6) * t136;
t71 = -qJD(5) * t118 - t137 * t242 - t175 * t244;
t72 = -qJD(5) * t119 + t137 * t244 - t175 * t242;
t33 = Icges(6,5) * t72 + Icges(6,6) * t71 + Icges(6,3) * t136;
t34 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t138;
t35 = Icges(6,4) * t72 + Icges(6,2) * t71 + Icges(6,6) * t136;
t36 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t138;
t37 = Icges(6,1) * t72 + Icges(6,4) * t71 + Icges(6,5) * t136;
t38 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t138;
t287 = (t118 * t38 - t119 * t36 + t136 * t56 + t155 * t34 + t59 * t71 + t62 * t72) * t153 + t155 * (t118 * t37 - t119 * t35 + t136 * t57 + t155 * t33 - t61 * t71 - t64 * t72);
t286 = t153 * (t115 * t36 + t116 * t38 + t138 * t56 + t153 * t34 + t59 * t73 + t62 * t74) + t155 * (t115 * t35 + t116 * t37 + t138 * t57 + t153 * t33 - t61 * t73 - t64 * t74);
t7 = t140 * t62 - t141 * t59 + t151 * t36 + t152 * t38 + t198 * t34;
t8 = -t140 * t64 + t141 * t61 + t151 * t35 + t152 * t37 + t198 * t33;
t285 = t153 * t7 + t155 * t8;
t284 = t233 - t303;
t197 = qJ(2) * t318 + t321;
t283 = -t162 - t197 - t232;
t215 = rSges(2,1) * t245 - t243 * rSges(2,2);
t282 = rSges(2,1) * t243 + rSges(2,2) * t245;
t281 = -rSges(3,2) * t239 + t344;
t280 = -t176 * rSges(4,1) - t175 * rSges(4,2);
t231 = qJDD(2) * t245;
t277 = -t243 * t311 + t231;
t39 = rSges(6,1) * t72 + rSges(6,2) * t71 + rSges(6,3) * t136;
t274 = t153 * t39 - t155 * t40;
t273 = t153 * t66 - t155 * t65;
t272 = t153 * (Icges(6,5) * t115 - Icges(6,6) * t116) + t155 * (-Icges(6,5) * t119 - Icges(6,6) * t118);
t270 = -t162 - t321;
t269 = -pkin(1) - t278;
t268 = -t183 + t284;
t264 = -t178 * pkin(3) - qJ(4) * t177 - t183;
t263 = -qJD(1) * t205 - t211 + t284;
t262 = -qJD(4) * t177 - qJDD(4) * t200 + t277;
t261 = -t351 - pkin(1) + (-rSges(4,3) - qJ(3)) * t239;
t260 = qJD(1) * (-t278 * t318 - t303) - qJDD(1) * t204 + t245 * t311 + t291;
t259 = -rSges(5,1) * t137 + rSges(5,2) * t136 + rSges(5,3) * t175;
t257 = (Icges(6,1) * t115 - t339 - t59) * t153 + (-Icges(6,1) * t119 - t338 + t61) * t155;
t256 = (-Icges(6,2) * t116 + t109 + t62) * t153 + (-Icges(6,2) * t118 - t110 - t64) * t155;
t255 = t269 * t349;
t254 = -qJ(2) * t319 + t284;
t253 = qJD(1) * t293 - t183 + t263;
t252 = -t108 + t270;
t251 = -t108 + t283 - t223;
t249 = qJD(1) * t264 - qJD(4) * t175 + qJDD(1) * t130 + qJDD(4) * t202 + t260;
t248 = t243 * t269 + t130 + t234;
t247 = t254 + t264;
t219 = rSges(3,2) * t304;
t172 = -t243 * t344 + t322;
t127 = qJD(1) * t323 + t233;
t126 = qJD(1) * t172 + t324;
t121 = -rSges(4,3) * t332 + t327;
t104 = rSges(6,1) * t151 - rSges(6,2) * t152;
t101 = Icges(6,5) * t151 - Icges(6,6) * t152;
t97 = rSges(6,1) * t140 - rSges(6,2) * t141;
t87 = qJD(1) * t306 + t284;
t86 = qJD(1) * t121 + t290;
t85 = t231 + t323 * qJDD(1) + (-rSges(3,3) * t318 - t197 + (qJD(1) * t281 - qJD(2)) * t243) * qJD(1);
t84 = qJDD(1) * t172 + (qJD(1) * t267 + t219) * qJD(1) + t291;
t82 = -rSges(6,1) * t119 - rSges(6,2) * t118;
t81 = rSges(6,1) * t115 - rSges(6,2) * t116;
t52 = qJD(1) * t292 + t268;
t43 = t306 * qJDD(1) + ((rSges(4,3) * t319 - t316) * t239 + t280 + t283) * qJD(1) + t277;
t42 = qJDD(1) * t121 + ((-rSges(4,3) * t318 - t317) * t239 + t328) * qJD(1) + t260;
t30 = -qJD(3) * t241 + qJD(4) * t333 + qJD(5) * t273;
t26 = t292 * qJDD(1) + (t251 + t259) * qJD(1) + t262;
t25 = (-t303 + t356) * qJD(1) + t249 + qJDD(1) * t92;
t24 = qJD(1) * t289 + t313 * t91 + t268 - t382;
t23 = -t314 * t91 - t357;
t13 = t115 * t95 + t116 * t96 + t138 * t88 + t153 * t94 + t73 * t89 + t74 * t90;
t12 = t118 * t96 - t119 * t95 + t136 * t88 + t155 * t94 + t71 * t89 + t72 * t90;
t11 = qJD(5) * t274 - t65 * t99 + t66 * t98 + t206;
t10 = t97 * t313 - t163 * t66 - t164 * t39 + t99 * t91 + t289 * qJDD(1) + (t251 + t288) * qJD(1) + t262;
t9 = t249 + (-t303 + t359) * qJD(1) - t97 * t314 - t98 * t91 + qJDD(1) * t105 + t163 * t65 + t164 * t40;
t1 = [t29 * t353 + t28 * t354 + t348 / 0.2e1 + t347 / 0.2e1 + t345 - m(2) * (-g(2) * t215 - g(3) * t282) + t5 * t299 + (t13 + t7) * t300 + (-t255 + t24 * (t252 + t288 - t39) + t23 * (t247 + t40 + t359) + (t10 * t269 + (-t24 * qJ(2) + t23 * t269) * qJD(1)) * t245 - t24 * t357 - t23 * (-qJD(1) * t106 + t253 - t382) - (t23 * t383 + t24 * t371) * qJD(5) + (-g(3) + t9) * (t105 + t248 + t65) + (-g(2) + t10) * (-t106 + t370 - t66)) * m(6) + (-(-qJD(1) * t93 + t253 - t52) * t51 - t255 + t52 * (t252 + t259) + t51 * (t247 + t356) + (t26 * t269 + (-t52 * qJ(2) + t269 * t51) * qJD(1)) * t245 + (-g(3) + t25) * (t248 + t92) + (-g(2) + t26) * (t370 - t93)) * m(5) + (-(qJD(1) * t361 + t263 - t87) * t86 + t87 * (rSges(4,3) * t305 + t270 + t280) + t86 * (t254 + t328) + (t43 * t261 + (-t87 * qJ(2) + t86 * (-rSges(4,3) * t239 + t269)) * qJD(1)) * t245 - t261 * t349 + (t43 - g(2)) * (t358 - t336) + (t42 - g(3)) * (t243 * t261 + t234 + t327)) * m(4) + (-t127 * t321 + t126 * (t219 + t233) + ((t126 * t302 + t127 * t343) * t245 + (t126 * t343 + t127 * t281) * t243) * qJD(1) - (-qJD(1) * t173 - t127 + t320) * t126 + (t85 - g(2)) * (t343 * t243 + t302 * t245 + t226) + (t84 - g(3)) * (t243 * t302 + t234 + t322)) * m(3) + (m(2) * (t215 ^ 2 + t282 ^ 2) + (Icges(5,1) * t199 + 0.2e1 * Icges(5,5) * t333) * t199 + (-0.2e1 * Icges(5,4) * t199 + Icges(5,2) * t198 - 0.2e1 * Icges(5,6) * t333) * t198 + Icges(3,2) * t355 + Icges(2,3) + (-Icges(4,5) * t296 + Icges(4,6) * t333 + Icges(4,3) * t241) * t241 + ((Icges(4,1) * t340 - Icges(4,4) * t238) * t296 - (Icges(4,4) * t340 - Icges(4,2) * t238) * t333 + (Icges(5,3) * t238 ^ 2 + Icges(3,1)) * t239 + (-Icges(4,5) * t340 + Icges(4,6) * t238 + (2 * Icges(3,4))) * t241) * t239) * qJDD(1) + (t5 + t12 + t8) * t298; (-m(3) + t312) * (g(3) * t243 + t349) + m(3) * (t243 * t84 + t245 * t85) + m(4) * (t243 * t42 + t245 * t43) + m(5) * (t243 * t25 + t245 * t26) + m(6) * (t10 * t245 + t243 * t9); t312 * (-g(1) * t241 + (-g(2) * t243 + g(3) * t245) * t239) + m(4) * (qJDD(3) * t355 + t331 * t42 - t332 * t43) + m(5) * (-t206 * t241 + t25 * t331 - t26 * t332) + m(6) * (-t10 * t332 - t11 * t241 + t331 * t9); t352 * (g(1) * t333 - g(2) * t200 + g(3) * t202) + m(5) * (-t175 * t51 - t177 * t52 - t200 * t26 + t202 * t25 + t206 * t333) + m(6) * (-t10 * t200 + t11 * t333 - t175 * t23 - t177 * t24 + t202 * t9) + 0.2e1 * (-m(5) * (-t200 * t51 - t202 * t52) / 0.2e1 - m(6) * (-t200 * t23 - t202 * t24) / 0.2e1) * qJD(1); t198 * (t285 * qJD(5) + t345 + t347 + t348) / 0.2e1 + t163 * (t153 * t20 + t155 * t21 + t198 * t32) / 0.2e1 + t164 * (t136 * t21 + t138 * t20 + t19 * t198 + t285) / 0.2e1 + t138 * t5 / 0.2e1 + t153 * (t286 * qJD(5) + t13 * t164 + t14 * t98 + t15 * t99 + t163 * t28) / 0.2e1 + (t198 * t28 + t276) * t354 + (t13 * t198 + t136 * t15 + t138 * t14 + t286) * t300 + t136 * (qJD(5) * t275 + t164 * t29) / 0.2e1 + t155 * (t287 * qJD(5) + t12 * t164 + t16 * t98 + t163 * t29 + t17 * t99) / 0.2e1 + (t198 * t29 + t275) * t353 + (t12 * t198 + t136 * t17 + t138 * t16 + t287) * t298 - t164 * ((t198 * t101 + t342 * t151 + t341 * t152) * t164 + (t151 * t256 + t152 * t257 + t198 * t272) * qJD(5)) / 0.2e1 - ((t153 * t101 + t115 * t342 + t116 * t341) * t164 + (t115 * t256 + t116 * t257 + t153 * t272) * qJD(5)) * t314 / 0.2e1 + ((t155 * t101 + t118 * t341 - t119 * t342) * t164 + (t118 * t257 - t119 * t256 + t155 * t272) * qJD(5)) * t299 + (t11 * t273 + t30 * (-t136 * t65 + t138 * t66 + t274) + t10 * (-t198 * t66 + t383) + t24 * (t136 * t91 + t155 * t97 - t198 * t39) + t9 * (t198 * t65 - t371) + t23 * (-t138 * t91 - t153 * t97 + t198 * t40) - (t23 * t81 - t24 * t82) * t164 - (t30 * (t153 * t82 - t155 * t81) + (-t153 * t23 + t155 * t24) * t104) * qJD(5) - g(1) * t104 - g(2) * t81 - g(3) * t82) * m(6);];
tau = t1;
