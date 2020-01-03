% Calculate vector of inverse dynamics joint torques for
% S5RPPPR1
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:05
% EndTime: 2020-01-03 11:20:24
% DurationCPUTime: 8.83s
% Computational Cost: add. (9352->489), mult. (9298->636), div. (0->0), fcn. (8638->10), ass. (0->235)
t208 = qJ(1) + pkin(7);
t203 = cos(t208);
t201 = sin(t208);
t313 = qJ(3) * t201;
t152 = pkin(2) * t203 + t313;
t215 = cos(qJ(1));
t206 = t215 * pkin(1);
t199 = qJD(1) * t206;
t287 = qJD(3) * t203;
t289 = qJD(1) * t203;
t290 = qJD(1) * t201;
t243 = pkin(2) * t289 + qJ(3) * t290 - t287;
t365 = -qJD(1) * t152 + t199 + t243;
t210 = sin(pkin(8));
t212 = cos(pkin(8));
t336 = pkin(3) * t212;
t240 = qJ(4) * t210 + t336;
t142 = t240 * t203;
t364 = -qJD(1) * t142 + t365;
t207 = pkin(9) + qJ(5);
t200 = sin(t207);
t202 = cos(t207);
t307 = t203 * t212;
t115 = t200 * t307 - t201 * t202;
t116 = t200 * t201 + t202 * t307;
t308 = t203 * t210;
t62 = Icges(6,5) * t116 - Icges(6,6) * t115 + Icges(6,3) * t308;
t101 = Icges(6,4) * t116;
t66 = Icges(6,2) * t115 - Icges(6,6) * t308 - t101;
t100 = Icges(6,4) * t115;
t68 = Icges(6,1) * t116 + Icges(6,5) * t308 - t100;
t26 = -(t200 * t66 + t202 * t68) * t210 + t212 * t62;
t286 = qJD(4) * t210;
t168 = t201 * t286;
t189 = -qJD(5) * t212 + qJD(1);
t112 = -rSges(6,3) * t212 + (rSges(6,1) * t202 - rSges(6,2) * t200) * t210;
t285 = qJD(5) * t210;
t268 = t112 * t285;
t71 = t116 * rSges(6,1) - t115 * rSges(6,2) + rSges(6,3) * t308;
t363 = t189 * t71 + t168 + t199 + (-qJD(3) - t268) * t203;
t310 = t201 * t212;
t113 = -t200 * t310 - t202 * t203;
t114 = -t203 * t200 + t202 * t310;
t311 = t201 * t210;
t61 = Icges(6,5) * t114 + Icges(6,6) * t113 + Icges(6,3) * t311;
t317 = Icges(6,4) * t114;
t64 = Icges(6,2) * t113 + Icges(6,6) * t311 + t317;
t99 = Icges(6,4) * t113;
t67 = Icges(6,1) * t114 + Icges(6,5) * t311 + t99;
t16 = t113 * t64 + t114 * t67 + t61 * t311;
t239 = t115 * t66 + t116 * t68;
t362 = -t239 + t16;
t151 = t201 * rSges(3,1) + t203 * rSges(3,2);
t214 = sin(qJ(1));
t205 = t214 * pkin(1);
t143 = t205 + t151;
t260 = t206 + t313;
t211 = cos(pkin(9));
t304 = t211 * t212;
t209 = sin(pkin(9));
t306 = t209 * t212;
t312 = t201 * t209;
t355 = -(-t201 * t211 + t203 * t306) * rSges(5,2) + (t203 * t304 + t312) * rSges(5,1);
t274 = rSges(5,3) * t308 + t355;
t361 = t260 + t274;
t359 = -t113 * t66 + t114 * t68;
t109 = -Icges(6,3) * t212 + (Icges(6,5) * t202 - Icges(6,6) * t200) * t210;
t315 = Icges(6,4) * t202;
t110 = -Icges(6,6) * t212 + (-Icges(6,2) * t200 + t315) * t210;
t316 = Icges(6,4) * t200;
t111 = -Icges(6,5) * t212 + (Icges(6,1) * t202 - t316) * t210;
t219 = t109 * t308 - t110 * t115 + t111 * t116;
t356 = t219 * t189;
t354 = -rSges(4,1) * t310 + rSges(4,2) * t311;
t309 = t203 * t209;
t353 = -(t201 * t304 - t309) * rSges(5,1) - (-t201 * t306 - t203 * t211) * rSges(5,2);
t140 = pkin(3) * t310 + qJ(4) * t311;
t194 = t201 * pkin(2);
t150 = -qJ(3) * t203 + t194;
t262 = t150 + t205;
t252 = t140 + t262;
t79 = rSges(5,3) * t311 - t353;
t51 = t252 + t79;
t266 = qJD(1) * t286;
t281 = qJDD(4) * t210;
t347 = t201 * t281 + t203 * t266;
t327 = rSges(4,3) * t203;
t346 = t327 + t354;
t221 = t201 * (-Icges(6,2) * t114 + t67 + t99) - t203 * (Icges(6,2) * t116 + t100 - t68);
t344 = t201 * (-Icges(6,1) * t113 + t317 + t64) - t203 * (-Icges(6,1) * t115 - t101 + t66);
t216 = qJD(1) ^ 2;
t342 = -m(5) - m(6);
t283 = qJD(1) * qJD(5);
t131 = (qJDD(5) * t201 + t203 * t283) * t210;
t341 = t131 / 0.2e1;
t132 = (-qJDD(5) * t203 + t201 * t283) * t210;
t340 = t132 / 0.2e1;
t339 = t201 / 0.2e1;
t338 = -t203 / 0.2e1;
t337 = -t212 / 0.2e1;
t335 = pkin(4) * t209;
t334 = g(2) * t203;
t188 = -qJDD(5) * t212 + qJDD(1);
t137 = (-Icges(6,5) * t200 - Icges(6,6) * t202) * t210;
t121 = qJD(5) * t137;
t138 = (-Icges(6,2) * t202 - t316) * t210;
t122 = qJD(5) * t138;
t139 = (-Icges(6,1) * t200 - t315) * t210;
t123 = qJD(5) * t139;
t28 = -t121 * t212 + (-t122 * t200 + t123 * t202 + (-t110 * t202 - t111 * t200) * qJD(5)) * t210;
t50 = -t109 * t212 + (-t110 * t200 + t111 * t202) * t210;
t333 = t50 * t188 + t28 * t189;
t332 = -t115 * t64 + t116 * t67;
t326 = t201 * t62;
t169 = t203 * t286;
t190 = qJD(3) * t201;
t295 = -t169 - t190;
t70 = t114 * rSges(6,1) + t113 * rSges(6,2) + rSges(6,3) * t311;
t195 = pkin(4) * t211 + pkin(3);
t213 = -pkin(6) - qJ(4);
t305 = t210 * t213;
t248 = t195 * t310 - t201 * t305;
t280 = pkin(4) * t309;
t77 = t248 - t280 - t140;
t23 = -t201 * t268 + t189 * t70 + (t252 + t77) * qJD(1) + t295;
t325 = t203 * t23;
t324 = t203 * t61;
t301 = t142 + t152;
t261 = (qJ(4) + t213) * t210;
t299 = -pkin(4) * t312 - t195 * t307;
t78 = (t261 + t336) * t203 + t299;
t277 = -t78 + t301;
t24 = t277 * qJD(1) + t363;
t323 = t24 * t201;
t25 = -t212 * t61 + (-t200 * t64 + t202 * t67) * t210;
t322 = t25 * t131;
t321 = t26 * t132;
t320 = rSges(4,3) + qJ(3);
t293 = qJ(3) * t289 + t190;
t118 = pkin(2) * t290 - t293;
t319 = -t240 * t290 - t118 + t169;
t176 = rSges(4,2) * t308;
t98 = rSges(4,1) * t307 + rSges(4,3) * t201 - t176;
t318 = t152 + t98;
t314 = pkin(1) * qJDD(1);
t303 = t110 - t139;
t302 = t111 + t138;
t269 = t212 * t289;
t278 = qJD(1) * t335;
t300 = t195 * t269 + t201 * t278;
t288 = qJD(1) * t210;
t270 = t203 * t288;
t298 = pkin(3) * t269 + qJ(4) * t270;
t297 = rSges(4,1) * t269 + rSges(4,3) * t290;
t292 = t194 + t205;
t291 = t216 * t206 + t214 * t314;
t284 = -m(4) + t342;
t282 = qJDD(3) * t201;
t17 = -t311 * t62 - t359;
t83 = -qJD(1) * t115 - qJD(5) * t114;
t84 = qJD(1) * t116 + qJD(5) * t113;
t45 = t84 * rSges(6,1) + t83 * rSges(6,2) + rSges(6,3) * t270;
t276 = t274 + t301;
t275 = rSges(5,3) * t270 + qJD(1) * t355;
t273 = t168 + t298;
t272 = t169 + t293;
t271 = t201 * t288;
t267 = rSges(4,1) * t212 + pkin(2);
t264 = -t285 / 0.2e1;
t263 = t285 / 0.2e1;
t153 = rSges(3,1) * t203 - t201 * rSges(3,2);
t259 = t199 - t287;
t258 = qJD(1) * t243 + qJDD(1) * t150 + t291;
t187 = -qJDD(4) * t212 + qJDD(2);
t257 = t201 * t264;
t256 = t201 * t263;
t255 = t203 * t264;
t254 = t203 * t263;
t253 = -t205 * t216 + t215 * t314;
t81 = qJD(1) * t113 + qJD(5) * t116;
t82 = qJD(1) * t114 + qJD(5) * t115;
t251 = rSges(6,1) * t82 + rSges(6,2) * t81;
t38 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t271;
t39 = Icges(6,5) * t84 + Icges(6,6) * t83 + Icges(6,3) * t270;
t40 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t271;
t41 = Icges(6,4) * t84 + Icges(6,2) * t83 + Icges(6,6) * t270;
t42 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t271;
t43 = Icges(6,1) * t84 + Icges(6,4) * t83 + Icges(6,5) * t270;
t249 = (t115 * t41 - t116 * t43 + t64 * t81 + t67 * t82 + (-t203 * t39 + t290 * t61) * t210) * t201 - (t115 * t40 - t116 * t42 + t66 * t81 - t68 * t82 + (-t203 * t38 - t290 * t62) * t210) * t203;
t247 = t201 * (t113 * t41 + t114 * t43 + t64 * t83 + t67 * t84 + (t201 * t39 + t289 * t61) * t210) - t203 * (t113 * t40 + t114 * t42 + t66 * t83 - t68 * t84 + (t201 * t38 - t289 * t62) * t210);
t218 = qJD(1) * t273 + qJDD(1) * t140 + t201 * t266 - t203 * t281 + t258;
t141 = (-rSges(6,1) * t200 - rSges(6,2) * t202) * t210;
t124 = qJD(5) * t141;
t233 = -t124 * t285 - qJDD(3);
t7 = qJDD(1) * t77 - t131 * t112 + t188 * t70 + t189 * t45 + t233 * t201 + ((-t213 * t288 - qJD(3)) * t203 - t298 + t300) * qJD(1) + t218;
t232 = qJD(1) * t190 + t253;
t44 = rSges(6,3) * t271 + t251;
t8 = t132 * t112 + t188 * t71 - t189 * t44 + t233 * t203 + t277 * qJDD(1) + (t203 * t278 + ((pkin(3) - t195) * t212 + t261) * t290 + t319) * qJD(1) + t232 + t347;
t246 = -t7 * t201 - t8 * t203;
t10 = -t212 * t38 + (-t200 * t40 + t202 * t42 + (t200 * t68 - t202 * t66) * qJD(5)) * t210;
t9 = -t212 * t39 + (-t200 * t41 + t202 * t43 + (-t200 * t67 - t202 * t64) * qJD(5)) * t210;
t245 = -t10 * t203 + t201 * t9;
t244 = t168 + t259;
t178 = rSges(2,1) * t215 - t214 * rSges(2,2);
t177 = rSges(2,1) * t214 + rSges(2,2) * t215;
t241 = t353 * qJD(1);
t238 = t16 * t203 + t17 * t201;
t237 = -t201 * t23 - t203 * t24;
t236 = t201 * t44 + t203 * t45;
t235 = -t201 * t71 + t203 * t70;
t234 = t201 * (Icges(6,5) * t113 - Icges(6,6) * t114) - t203 * (Icges(6,5) * t115 + Icges(6,6) * t116);
t231 = pkin(2) + t240;
t224 = (t16 * t201 - t17 * t203) * t210;
t18 = -t308 * t61 - t332;
t19 = t308 * t62 + t239;
t223 = (t18 * t201 - t19 * t203) * t210;
t220 = -qJDD(3) * t203 + t232;
t144 = t153 + t206;
t92 = rSges(6,1) * t115 + rSges(6,2) * t116;
t91 = rSges(6,1) * t113 - rSges(6,2) * t114;
t74 = qJD(1) * t318 + t259;
t47 = qJD(1) * t276 + t244;
t34 = t318 * qJDD(1) + (qJD(1) * t346 - t118) * qJD(1) + t220;
t33 = -qJDD(1) * t346 - t282 + ((-rSges(4,2) * t288 - qJD(3)) * t203 + t297) * qJD(1) + t258;
t31 = t109 * t311 + t110 * t113 + t111 * t114;
t30 = t31 * t189;
t29 = -qJD(4) * t212 + t235 * t285 + qJD(2);
t15 = t276 * qJDD(1) + (-rSges(5,3) * t271 + t241 + t319) * qJD(1) + t220 + t347;
t14 = qJDD(1) * t79 - t282 + (t275 - t287) * qJD(1) + t218;
t13 = t110 * t83 + t111 * t84 + t113 * t122 + t114 * t123 + (t109 * t289 + t121 * t201) * t210;
t12 = t110 * t81 + t111 * t82 + t115 * t122 - t116 * t123 + (t109 * t290 - t121 * t203) * t210;
t11 = -t131 * t71 - t132 * t70 + t236 * t285 + t187;
t6 = qJD(5) * t223 - t356;
t5 = qJD(5) * t224 + t30;
t1 = [t333 + t321 / 0.2e1 - t219 * t340 + (t30 + ((t19 + t362) * t201 + (t18 + (t324 - t326) * t210 - t17 + t332) * t203) * t285) * t254 + t31 * t341 + t322 / 0.2e1 - m(2) * (g(2) * t178 + g(3) * t177) + (t13 + t9) * t256 + ((-t151 * t216 - g(2) + t253) * t144 + (t291 - g(3) + (0.2e1 * rSges(3,1) * t289 - 0.2e1 * rSges(3,2) * t290 - qJD(1) * t153) * qJD(1)) * t143) * m(3) + (t356 + ((t332 + t359) * t201 - t362 * t203 + ((t324 + t326) * t201 + t203 ^ 2 * t62) * t210 + t238) * t285 + t6) * t257 + (t10 + t12 + t5) * t255 + (t24 * (-t251 + t272) + (t24 * (t280 - t205) - t305 * t325 + (-t195 * t212 - pkin(2) + (-rSges(6,3) + t213) * t210) * t323) * qJD(1) + (t7 - g(3)) * ((-qJ(3) - t335) * t203 + t248 + t70 + t292) + (t8 - g(2)) * ((pkin(2) - t305) * t203 + t260 - t299 + t71) + (qJD(1) * t78 + t168 + t24 + t300 - t363 + t364 + t45) * t23) * m(6) + (-g(2) * t361 - t231 * t334 + (t241 + t272 + (-t205 + (-t336 - pkin(2) + (-rSges(5,3) - qJ(4)) * t210) * t201) * qJD(1)) * t47 + (-g(3) + t14) * t51 + (t231 * t203 + t361) * t15 + (-qJD(1) * t274 - t244 + t273 + t275 + t364 + t47) * (qJD(1) * t51 + t295)) * m(5) + ((t34 - g(2)) * (t201 * t320 + t203 * t267 - t176 + t206) + (t33 - g(3)) * (-t203 * t320 + t292 - t354) + (t293 + (t327 - t205 + (rSges(4,2) * t210 - t267) * t201) * qJD(1)) * t74 + (-t259 + t297 + t74 + (-t176 - t98) * qJD(1) + t365) * (-t190 + (t262 - t346) * qJD(1))) * m(4) + (m(3) * (t143 * t151 + t144 * t153) + m(2) * (t177 ^ 2 + t178 ^ 2) + Icges(2,3) + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t212 ^ 2 + ((Icges(4,1) + Icges(5,1) * t211 ^ 2 + (-0.2e1 * Icges(5,4) * t211 + Icges(5,2) * t209) * t209) * t210 + 0.2e1 * (-Icges(5,5) * t211 + Icges(5,6) * t209 + Icges(4,4)) * t212) * t210) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t187 + m(6) * t11 + (-m(3) + t284) * g(1); t284 * (-g(3) * t201 - t334) + m(4) * (-t201 * t33 - t203 * t34) + m(5) * (-t14 * t201 - t15 * t203) + m(6) * t246; t342 * (-g(1) * t212 + (g(2) * t201 - g(3) * t203) * t210) + 0.2e1 * (t11 * t337 + (t338 * t7 + t339 * t8) * t210) * m(6) + 0.2e1 * (t187 * t337 + (t14 * t338 + t15 * t339) * t210) * m(5); (t245 * t285 + t321 + t322 + t333) * t337 + t188 * (-t212 * t50 + (t201 * t25 - t203 * t26) * t210) / 0.2e1 + t189 * (-t212 * t28 + ((t26 * t201 + t203 * t25) * qJD(1) + t245) * t210) / 0.2e1 + (t13 * t189 + t131 * t16 + t132 * t17 + t188 * t31 + t247 * t285) * t311 / 0.2e1 + (-t212 * t31 + t224) * t341 + (-t13 * t212 + (qJD(1) * t238 + t247) * t210) * t256 - (t12 * t189 + t131 * t18 + t132 * t19 - t188 * t219 + t249 * t285) * t308 / 0.2e1 + (t212 * t219 + t223) * t340 + (-t12 * t212 + ((t18 * t203 + t19 * t201) * qJD(1) + t249) * t210) * t255 - t189 * (-t212 * t137 * t189 + ((-t200 * t302 - t202 * t303) * t189 + ((-t200 * t221 - t202 * t344) * t210 - t234 * t212) * qJD(5)) * t210) / 0.2e1 + ((t113 * t302 - t114 * t303 + t137 * t311) * t189 + (t113 * t221 - t114 * t344 + t234 * t311) * t285) * t257 + ((t115 * t302 + t116 * t303 - t137 * t308) * t189 + (t221 * t115 + t116 * t344 - t234 * t308) * t285) * t254 + (t201 * t6 + t203 * t5) * t288 / 0.2e1 + ((-t23 * t45 + t24 * t44 - t7 * t70 - t71 * t8) * t212 + (t11 * t235 + t29 * (-t289 * t71 - t290 * t70 + t236) + t237 * t124 + ((t323 - t325) * qJD(1) + t246) * t112) * t210 - (t23 * t91 - t24 * t92) * t189 - (t29 * (t201 * t92 + t203 * t91) + t237 * t141) * t285 - g(1) * t141 - g(2) * t91 - g(3) * t92) * m(6);];
tau = t1;
