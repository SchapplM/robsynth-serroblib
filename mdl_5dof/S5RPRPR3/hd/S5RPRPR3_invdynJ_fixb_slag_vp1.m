% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:56
% EndTime: 2020-01-03 11:36:14
% DurationCPUTime: 6.84s
% Computational Cost: add. (11923->439), mult. (9980->577), div. (0->0), fcn. (9392->10), ass. (0->213)
t215 = qJ(1) + pkin(8);
t209 = qJ(3) + t215;
t200 = sin(t209);
t183 = qJD(4) * t200;
t201 = cos(t209);
t214 = qJD(1) + qJD(3);
t296 = t201 * t214;
t283 = qJ(4) * t296 + t183;
t299 = t200 * t214;
t106 = pkin(3) * t299 - t283;
t216 = sin(pkin(9));
t275 = qJD(5) * t214;
t113 = (-qJDD(5) * t201 + t200 * t275) * t216;
t217 = cos(pkin(9));
t218 = sin(qJ(5));
t220 = cos(qJ(5));
t137 = -rSges(6,3) * t217 + (rSges(6,1) * t220 - rSges(6,2) * t218) * t216;
t213 = qJDD(1) + qJDD(3);
t181 = -qJDD(5) * t217 + t213;
t182 = -qJD(5) * t217 + t214;
t221 = cos(qJ(1));
t301 = pkin(1) * qJDD(1);
t203 = t221 * t301;
t208 = cos(t215);
t222 = qJD(1) ^ 2;
t300 = pkin(2) * qJDD(1);
t207 = sin(t215);
t219 = sin(qJ(1));
t211 = t219 * pkin(1);
t327 = pkin(2) * t207 + t211;
t232 = t208 * t300 - t222 * t327 + t203;
t229 = t214 * t183 + t232;
t152 = (-rSges(6,1) * t218 - rSges(6,2) * t220) * t216;
t143 = qJD(5) * t152;
t274 = qJD(5) * t216;
t247 = -t143 * t274 - qJDD(4);
t294 = t201 * t217;
t295 = t201 * t216;
t131 = pkin(4) * t294 + pkin(7) * t295;
t146 = t201 * pkin(3) + t200 * qJ(4);
t286 = t146 + t131;
t317 = pkin(4) * t217;
t292 = t217 * t218;
t124 = -t200 * t292 - t201 * t220;
t291 = t217 * t220;
t127 = t200 * t218 + t201 * t291;
t87 = qJD(5) * t127 + t124 * t214;
t125 = t200 * t291 - t201 * t218;
t126 = -t200 * t220 + t201 * t292;
t88 = qJD(5) * t126 + t125 * t214;
t259 = rSges(6,1) * t88 + rSges(6,2) * t87;
t293 = t214 * t216;
t272 = t200 * t293;
t49 = rSges(6,3) * t272 + t259;
t81 = t127 * rSges(6,1) - t126 * rSges(6,2) + rSges(6,3) * t295;
t17 = t113 * t137 + t181 * t81 - t182 * t49 + t286 * t213 + t247 * t201 + (-t106 + (-pkin(7) * t216 - t317) * t299) * t214 + t229;
t348 = t17 - g(2);
t159 = rSges(3,1) * t208 - t207 * rSges(3,2);
t212 = t221 * pkin(1);
t347 = t159 + t212;
t72 = Icges(6,5) * t127 - Icges(6,6) * t126 + Icges(6,3) * t295;
t117 = Icges(6,4) * t127;
t76 = Icges(6,2) * t126 - Icges(6,6) * t295 - t117;
t116 = Icges(6,4) * t126;
t78 = Icges(6,1) * t127 + Icges(6,5) * t295 - t116;
t31 = -(t218 * t76 + t220 * t78) * t216 + t217 * t72;
t252 = t126 * t76 + t127 * t78;
t27 = t295 * t72 + t252;
t345 = t201 * t27;
t344 = t81 + t286;
t328 = t200 * rSges(4,1) + t201 * rSges(4,2);
t122 = t328 * t214;
t323 = t327 * qJD(1);
t102 = t323 + t122;
t269 = t137 * t274;
t343 = t182 * t81 + (-qJD(4) - t269) * t201;
t342 = -t124 * t76 + t125 * t78;
t133 = -Icges(6,3) * t217 + (Icges(6,5) * t220 - Icges(6,6) * t218) * t216;
t302 = Icges(6,4) * t220;
t134 = -Icges(6,6) * t217 + (-Icges(6,2) * t218 + t302) * t216;
t303 = Icges(6,4) * t218;
t135 = -Icges(6,5) * t217 + (Icges(6,1) * t220 - t303) * t216;
t231 = -t126 * t134 + t127 * t135 + t133 * t295;
t338 = t231 * t182;
t112 = (qJDD(5) * t200 + t201 * t275) * t216;
t297 = t200 * t217;
t298 = t200 * t216;
t130 = pkin(4) * t297 + pkin(7) * t298;
t193 = t200 * pkin(3);
t144 = -qJ(4) * t201 + t193;
t276 = qJD(4) * t201;
t284 = pkin(3) * t296 + qJ(4) * t299;
t256 = -t276 + t284;
t199 = pkin(2) * t208;
t279 = t222 * t212 + t219 * t301;
t265 = t222 * t199 + t207 * t300 + t279;
t246 = t213 * t144 + t214 * t256 + t265;
t270 = t214 * t294;
t271 = t201 * t293;
t329 = pkin(4) * t270 + pkin(7) * t271;
t89 = -qJD(5) * t125 - t126 * t214;
t90 = qJD(5) * t124 + t127 * t214;
t50 = t90 * rSges(6,1) + t89 * rSges(6,2) + rSges(6,3) * t271;
t80 = t125 * rSges(6,1) + t124 * rSges(6,2) + rSges(6,3) * t298;
t16 = -t112 * t137 + t213 * t130 + t181 * t80 + t182 * t50 + (-t276 + t329) * t214 + t247 * t200 + t246;
t336 = t16 - g(3);
t170 = rSges(5,1) * t297;
t260 = -rSges(5,2) * t298 + t170;
t108 = -rSges(5,3) * t201 + t260;
t233 = rSges(5,1) * t270 + rSges(5,3) * t299 + (-rSges(5,2) * t293 - qJD(4)) * t201;
t32 = -qJDD(4) * t200 + t213 * t108 + t214 * t233 + t246;
t335 = t32 - g(3);
t109 = rSges(5,1) * t294 - rSges(5,2) * t295 + rSges(5,3) * t200;
t100 = t146 + t109;
t285 = rSges(5,2) * t272 + rSges(5,3) * t296;
t33 = -qJDD(4) * t201 + (-t170 * t214 - t106 + t285) * t214 + t100 * t213 + t229;
t334 = t33 - g(2);
t123 = rSges(4,1) * t296 - rSges(4,2) * t299;
t333 = t123 * t214 + t213 * t328 - g(3) + t265;
t147 = rSges(4,1) * t201 - t200 * rSges(4,2);
t332 = -t122 * t214 + t147 * t213 - g(2) + t232;
t277 = qJD(1) * t208;
t310 = pkin(1) * qJD(1);
t282 = pkin(2) * t277 + t221 * t310;
t36 = t214 * t286 + t282 + t343;
t331 = t36 * (-t317 - pkin(3) + (-rSges(6,3) - pkin(7)) * t216) * t299;
t330 = t214 * (t108 + t144);
t158 = t207 * rSges(3,1) + t208 * rSges(3,2);
t136 = t214 * t146;
t325 = -t214 * t131 - t136 + t256 + t329 - t343 + t50;
t322 = -t214 * t109 - t136 + t233 + t284;
t289 = t130 + t144;
t321 = t182 * t80 - t200 * t269 + t214 * t289;
t115 = Icges(6,4) * t124;
t77 = Icges(6,1) * t125 + Icges(6,5) * t298 + t115;
t235 = t200 * (-Icges(6,2) * t125 + t115 + t77) - t201 * (Icges(6,2) * t127 + t116 - t78);
t304 = Icges(6,4) * t125;
t74 = Icges(6,2) * t124 + Icges(6,6) * t298 + t304;
t320 = t200 * (-Icges(6,1) * t124 + t304 + t74) - t201 * (-Icges(6,1) * t126 - t117 + t76);
t319 = t112 / 0.2e1;
t318 = t113 / 0.2e1;
t149 = (-Icges(6,5) * t218 - Icges(6,6) * t220) * t216;
t138 = qJD(5) * t149;
t150 = (-Icges(6,2) * t220 - t303) * t216;
t139 = qJD(5) * t150;
t151 = (-Icges(6,1) * t218 - t302) * t216;
t140 = qJD(5) * t151;
t41 = -t138 * t217 + (-t139 * t218 + t140 * t220 + (-t134 * t220 - t135 * t218) * qJD(5)) * t216;
t63 = -t133 * t217 + (-t134 * t218 + t135 * t220) * t216;
t316 = t63 * t181 + t41 * t182;
t315 = -t126 * t74 + t127 * t77;
t308 = t200 * t72;
t71 = Icges(6,5) * t125 + Icges(6,6) * t124 + Icges(6,3) * t298;
t307 = t201 * t71;
t30 = -t217 * t71 + (-t218 * t74 + t220 * t77) * t216;
t306 = t30 * t112;
t305 = t31 * t113;
t288 = t134 - t151;
t287 = t135 + t150;
t280 = t199 + t212;
t278 = qJD(1) * t207;
t273 = m(3) + m(4) + m(5);
t24 = t124 * t74 + t125 * t77 + t71 * t298;
t25 = -t298 * t72 - t342;
t267 = -t274 / 0.2e1;
t266 = t274 / 0.2e1;
t264 = t200 * t267;
t263 = t200 * t266;
t262 = t201 * t267;
t261 = t201 * t266;
t103 = t147 * t214 + t282;
t255 = -t276 + t282;
t178 = rSges(2,1) * t221 - t219 * rSges(2,2);
t177 = rSges(2,1) * t219 + rSges(2,2) * t221;
t234 = t323 - t183;
t35 = t234 + t321;
t251 = -t200 * t35 - t201 * t36;
t250 = t200 * t49 + t201 * t50;
t249 = -t200 * t81 + t201 * t80;
t248 = t200 * (Icges(6,5) * t124 - Icges(6,6) * t125) - t201 * (Icges(6,5) * t126 + Icges(6,6) * t127);
t240 = (t200 * t24 - t201 * t25) * t216;
t26 = -t295 * t71 - t315;
t239 = (t200 * t26 - t345) * t216;
t238 = -t259 + t283;
t54 = t80 + t289;
t99 = t193 + (-rSges(5,3) - qJ(4)) * t201 + t260;
t225 = (-rSges(5,1) * t217 - pkin(3)) * t299 + t283 + t285;
t10 = qJD(5) * t239 - t338;
t44 = Icges(6,5) * t90 + Icges(6,6) * t89 + Icges(6,3) * t271;
t46 = Icges(6,4) * t90 + Icges(6,2) * t89 + Icges(6,6) * t271;
t48 = Icges(6,1) * t90 + Icges(6,4) * t89 + Icges(6,5) * t271;
t13 = -t217 * t44 + (-t218 * t46 + t220 * t48 + (-t218 * t77 - t220 * t74) * qJD(5)) * t216;
t43 = Icges(6,5) * t88 + Icges(6,6) * t87 + Icges(6,3) * t272;
t45 = Icges(6,4) * t88 + Icges(6,2) * t87 + Icges(6,6) * t272;
t47 = Icges(6,1) * t88 + Icges(6,4) * t87 + Icges(6,5) * t272;
t14 = -t217 * t43 + (-t218 * t45 + t220 * t47 + (t218 * t78 - t220 * t76) * qJD(5)) * t216;
t21 = t126 * t139 - t127 * t140 + t134 * t87 + t135 * t88 + (t133 * t299 - t138 * t201) * t216;
t22 = t124 * t139 + t125 * t140 + t134 * t89 + t135 * t90 + (t133 * t296 + t138 * t200) * t216;
t51 = t124 * t134 + t125 * t135 + t133 * t298;
t42 = t51 * t182;
t9 = qJD(5) * t240 + t42;
t223 = (t42 + ((t24 - t252 + t27) * t200 + (t26 + (t307 - t308) * t216 - t25 + t315) * t201) * t274) * t261 + t306 / 0.2e1 + t305 / 0.2e1 + t51 * t319 - t231 * t318 + t316 + (t13 + t22) * t263 + (Icges(5,2) * t217 ^ 2 + (Icges(5,1) * t216 + 0.2e1 * Icges(5,4) * t217) * t216 + Icges(4,3)) * t213 + (t338 + (t345 + (t315 + t25 + (t307 + t308) * t216 + t342) * t200) * t274 + t10) * t264 + (t14 + t21 + t9) * t262;
t98 = rSges(6,1) * t126 + rSges(6,2) * t127;
t97 = rSges(6,1) * t124 - rSges(6,2) * t125;
t65 = t100 * t214 + t255;
t64 = t234 + t330;
t39 = t249 * t274 + qJD(2);
t15 = -t112 * t81 - t113 * t80 + t250 * t274 + qJDD(2);
t6 = t124 * t45 + t125 * t47 + t76 * t89 - t78 * t90 + (t200 * t43 - t296 * t72) * t216;
t5 = t124 * t46 + t125 * t48 + t74 * t89 + t77 * t90 + (t200 * t44 + t296 * t71) * t216;
t4 = t126 * t45 - t127 * t47 + t76 * t87 - t78 * t88 + (-t201 * t43 - t299 * t72) * t216;
t3 = t126 * t46 - t127 * t48 + t74 * t87 + t77 * t88 + (-t201 * t44 + t299 * t71) * t216;
t1 = [t223 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t332 * (t147 + t280) + t333 * (t327 + t328) + (-t103 + t123 + t282) * t102) * m(4) + ((qJDD(1) * t159 - g(2) + t203) * t347 + (-t222 * t347 + qJDD(1) * t158 - g(3) + t279 + (0.2e1 * rSges(3,1) * t277 - 0.2e1 * rSges(3,2) * t278 - qJD(1) * t159) * qJD(1)) * (t211 + t158)) * m(3) + ((t177 ^ 2 + t178 ^ 2) * qJDD(1) - g(2) * t178 - g(3) * t177) * m(2) + (t36 * (-pkin(2) * t278 - t219 * t310 + t238) + t336 * (t54 + t327) + (t36 + t325) * t35 + t331 + t348 * (t344 + t280)) * m(6) + (t65 * (-t323 + t225) + t334 * (t100 + t280) + t335 * (t99 + t327) + (t282 - t255 + t65 + t322) * t64) * m(5); m(6) * t15 + t273 * qJDD(2) + (-m(6) - t273) * g(1); t223 + (t336 * t54 + (-t183 + t238 + t321) * t36 + t325 * t35 + t331 + t348 * t344) * m(6) + (t335 * t99 + (-t183 + t225 + t330) * t65 + (t276 + t322) * t64 + t334 * t100) * m(5) + (t102 * t123 - t103 * t122 + (-t102 * t214 + t332) * t147 + (t103 * t214 + t333) * t328) * m(4); (-m(5) - m(6)) * (-g(2) * t201 - g(3) * t200) + m(5) * (-t200 * t32 - t201 * t33) + m(6) * (-t16 * t200 - t17 * t201); -t217 * (t306 + t305 + (t13 * t200 - t14 * t201) * t274 + t316) / 0.2e1 + t181 * (-t217 * t63 + (t200 * t30 - t201 * t31) * t216) / 0.2e1 + t182 * (-t217 * t41 + ((t214 * t30 - t14) * t201 + (t214 * t31 + t13) * t200) * t216) / 0.2e1 + (t112 * t24 + t113 * t25 + t181 * t51 + t182 * t22 + (t200 * t5 - t201 * t6) * t274) * t298 / 0.2e1 + (-t51 * t217 + t240) * t319 + (-t217 * t22 + ((t214 * t24 - t6) * t201 + (t214 * t25 + t5) * t200) * t216) * t263 - (t112 * t26 + t113 * t27 - t181 * t231 + t182 * t21 + (t200 * t3 - t201 * t4) * t274) * t295 / 0.2e1 + (t217 * t231 + t239) * t318 + (-t21 * t217 + ((t214 * t26 - t4) * t201 + (t214 * t27 + t3) * t200) * t216) * t262 - t182 * (-t217 * t149 * t182 + ((-t218 * t287 - t220 * t288) * t182 + ((-t218 * t235 - t220 * t320) * t216 - t248 * t217) * qJD(5)) * t216) / 0.2e1 + ((t124 * t287 - t125 * t288 + t149 * t298) * t182 + (t124 * t235 - t125 * t320 + t248 * t298) * t274) * t264 + ((t126 * t287 + t127 * t288 - t149 * t295) * t182 + (t235 * t126 + t127 * t320 - t248 * t295) * t274) * t261 + (t200 * t10 + t201 * t9) * t293 / 0.2e1 + ((-t16 * t80 - t17 * t81 - t35 * t50 + t36 * t49) * t217 + (t15 * t249 + t39 * (-t296 * t81 - t299 * t80 + t250) + t251 * t143 + ((-t214 * t35 - t17) * t201 + (t214 * t36 - t16) * t200) * t137) * t216 - (t35 * t97 - t36 * t98) * t182 - (t39 * (t200 * t98 + t201 * t97) + t251 * t152) * t274 - g(1) * t152 - g(2) * t97 - g(3) * t98) * m(6);];
tau = t1;
