% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:03
% EndTime: 2019-12-05 18:20:16
% DurationCPUTime: 7.73s
% Computational Cost: add. (12170->472), mult. (10102->595), div. (0->0), fcn. (9456->10), ass. (0->244)
t202 = sin(pkin(9));
t203 = cos(pkin(9));
t204 = sin(qJ(5));
t206 = cos(qJ(5));
t145 = -rSges(6,3) * t203 + (rSges(6,1) * t206 - rSges(6,2) * t204) * t202;
t282 = qJD(5) * t202;
t350 = t145 * t282;
t201 = qJ(1) + qJ(2);
t194 = pkin(8) + t201;
t189 = sin(t194);
t190 = cos(t194);
t200 = qJD(1) + qJD(2);
t283 = qJD(5) * t200;
t116 = (qJDD(5) * t190 - t189 * t283) * t202;
t199 = qJDD(1) + qJDD(2);
t181 = -qJDD(5) * t203 + t199;
t182 = -qJD(5) * t203 + t200;
t198 = t200 ^ 2;
t208 = qJD(1) ^ 2;
t207 = cos(qJ(1));
t326 = pkin(1) * t207;
t205 = sin(qJ(1));
t327 = pkin(1) * t205;
t253 = -qJDD(1) * t326 + t208 * t327;
t195 = sin(t201);
t325 = pkin(2) * t195;
t235 = t198 * t325 + t253;
t229 = qJDD(4) * t190 + t235;
t323 = pkin(4) * t203;
t252 = pkin(7) * t202 + t323;
t134 = t252 * t190;
t301 = qJ(4) * t189;
t153 = pkin(3) * t190 + t301;
t196 = cos(t201);
t324 = pkin(2) * t196;
t265 = -t153 - t324;
t255 = -t134 + t265;
t185 = qJD(4) * t189;
t299 = t189 * t200;
t173 = pkin(3) * t299;
t284 = t185 - t173;
t297 = t190 * t200;
t263 = -qJ(4) * t297 - t185 - t284;
t159 = (-rSges(6,1) * t204 - rSges(6,2) * t206) * t202;
t149 = qJD(5) * t159;
t270 = t149 * t282;
t274 = t203 * t299;
t293 = t200 * t202;
t275 = t189 * t293;
t286 = pkin(4) * t274 + pkin(7) * t275;
t292 = t203 * t204;
t130 = t189 * t292 + t190 * t206;
t291 = t203 * t206;
t133 = t189 * t204 + t190 * t291;
t89 = -qJD(5) * t133 + t130 * t200;
t131 = t189 * t291 - t190 * t204;
t132 = -t189 * t206 + t190 * t292;
t90 = -qJD(5) * t132 - t131 * t200;
t49 = rSges(6,1) * t90 + rSges(6,2) * t89 - rSges(6,3) * t275;
t289 = -t133 * rSges(6,1) + t132 * rSges(6,2);
t296 = t190 * t202;
t81 = rSges(6,3) * t296 - t289;
t17 = t190 * t270 + t116 * t145 - t181 * t81 - t182 * t49 + (t263 + t286) * t200 + t255 * t199 + t229;
t349 = -g(2) + t17;
t348 = -t182 * t81 + t190 * t350;
t249 = rSges(3,1) * t195 + rSges(3,2) * t196;
t143 = t249 * t200;
t310 = pkin(1) * qJD(1);
t277 = t205 * t310;
t126 = t143 + t277;
t245 = -t301 - t324;
t346 = t245 + t289;
t121 = Icges(6,4) * t133;
t75 = -Icges(6,2) * t132 + Icges(6,6) * t296 + t121;
t120 = Icges(6,4) * t132;
t79 = -Icges(6,1) * t133 - Icges(6,5) * t296 + t120;
t320 = t130 * t75 + t131 * t79;
t243 = -t132 * t75 - t133 * t79;
t72 = Icges(6,5) * t133 - Icges(6,6) * t132 + Icges(6,3) * t296;
t31 = -t203 * t72 + (-t204 * t75 - t206 * t79) * t202;
t115 = (-qJDD(5) * t189 - t190 * t283) * t202;
t186 = qJD(4) * t190;
t228 = (-qJDD(1) * t205 - t207 * t208) * pkin(1);
t213 = t228 + (-t195 * t199 - t196 * t198) * pkin(2);
t300 = qJ(4) * t190;
t244 = -pkin(3) * t189 + t300;
t141 = t200 * t153;
t333 = t186 - t141;
t211 = qJDD(4) * t189 + t199 * t244 + t213 + (t186 + t333) * t200;
t273 = t190 * t293;
t91 = qJD(5) * t131 + t132 * t200;
t92 = qJD(5) * t130 - t133 * t200;
t315 = t92 * rSges(6,1) + t91 * rSges(6,2);
t50 = -rSges(6,3) * t273 + t315;
t290 = -t131 * rSges(6,1) + t130 * rSges(6,2);
t298 = t189 * t202;
t80 = -rSges(6,3) * t298 + t290;
t16 = -t115 * t145 + t181 * t80 + t182 * t50 - t198 * t134 + (-t199 * t252 + t270) * t189 + t211;
t345 = t16 - g(3);
t162 = rSges(5,2) * t273;
t313 = rSges(5,1) * t203;
t247 = rSges(5,2) * t202 - t313;
t311 = rSges(5,3) * t190;
t217 = t189 * t247 + t311;
t312 = rSges(5,3) * t189;
t234 = -t190 * t313 - t312;
t32 = t199 * t217 + (t200 * t234 + t162) * t200 + t211;
t344 = t32 - g(3);
t161 = rSges(5,1) * t274;
t174 = rSges(5,2) * t296;
t112 = -t174 - t234;
t256 = -t112 + t265;
t279 = rSges(5,2) * t298;
t33 = t256 * t199 + (t161 + (-t279 - t311) * t200 + t263) * t200 + t229;
t343 = t33 - g(2);
t170 = rSges(4,2) * t299;
t248 = -rSges(4,1) * t189 - rSges(4,2) * t190;
t342 = t199 * t248 + t200 * (-rSges(4,1) * t297 + t170) + t213 - g(3);
t187 = t189 * rSges(4,2);
t314 = rSges(4,1) * t190;
t154 = -t187 + t314;
t264 = -t154 - t324;
t285 = -rSges(4,1) * t299 - rSges(4,2) * t297;
t341 = t199 * t264 - t200 * t285 - g(2) + t235;
t294 = t196 * t200;
t295 = t195 * t200;
t144 = -rSges(3,1) * t294 + rSges(3,2) * t295;
t340 = t144 * t200 - t199 * t249 - g(3) + t228;
t166 = rSges(3,1) * t196 - t195 * rSges(3,2);
t339 = t143 * t200 - t166 * t199 - g(2) + t253;
t337 = -t200 * t112 - t162;
t251 = -t314 - t324;
t336 = -t170 + (-t154 - t251) * t200;
t280 = pkin(2) * t294;
t257 = t186 - t280;
t335 = qJ(4) * t299 - t200 * t134 - t257 - t315 + t348;
t227 = -t323 - pkin(3) + (-rSges(6,3) - pkin(7)) * t202;
t322 = g(2) * t190;
t139 = t200 * t244;
t281 = pkin(2) * t295;
t226 = t277 + t281;
t220 = -t139 - t185 + t226;
t332 = -t182 * t80 - t189 * t350 + t252 * t299;
t35 = t220 + t332;
t276 = t207 * t310;
t254 = t186 - t276;
t36 = t200 * t255 + t254 + t348;
t334 = (t17 * t227 + (-t36 * qJ(4) - t35 * (-rSges(6,3) * t202 - pkin(3) - t252)) * t200) * t190 - t227 * t322;
t71 = -Icges(6,5) * t131 + Icges(6,6) * t130 - Icges(6,3) * t298;
t304 = Icges(6,4) * t131;
t74 = Icges(6,2) * t130 - Icges(6,6) * t298 - t304;
t119 = Icges(6,4) * t130;
t77 = -Icges(6,1) * t131 - Icges(6,5) * t298 + t119;
t26 = -t132 * t74 + t133 * t77 + t71 * t296;
t221 = t189 * (-Icges(6,1) * t130 - t304 + t74) + t190 * (-Icges(6,1) * t132 - t121 - t75);
t223 = t189 * (Icges(6,2) * t131 + t119 + t77) - t190 * (-Icges(6,2) * t133 - t120 - t79);
t330 = -m(5) - m(6);
t329 = t115 / 0.2e1;
t328 = t116 / 0.2e1;
t302 = Icges(6,4) * t206;
t137 = -Icges(6,6) * t203 + (-Icges(6,2) * t204 + t302) * t202;
t303 = Icges(6,4) * t204;
t138 = -Icges(6,5) * t203 + (Icges(6,1) * t206 - t303) * t202;
t156 = (-Icges(6,5) * t204 - Icges(6,6) * t206) * t202;
t146 = qJD(5) * t156;
t157 = (-Icges(6,2) * t206 - t303) * t202;
t147 = qJD(5) * t157;
t158 = (-Icges(6,1) * t204 - t302) * t202;
t148 = qJD(5) * t158;
t41 = -t146 * t203 + (-t147 * t204 + t148 * t206 + (-t137 * t206 - t138 * t204) * qJD(5)) * t202;
t136 = -Icges(6,3) * t203 + (Icges(6,5) * t206 - Icges(6,6) * t204) * t202;
t67 = -t136 * t203 + (-t137 * t204 + t138 * t206) * t202;
t321 = t67 * t181 + t41 * t182;
t319 = -t130 * t74 + t131 * t77;
t52 = -t132 * t137 + t133 * t138 + t136 * t296;
t309 = t182 * t52;
t30 = -t203 * t71 + (-t204 * t74 + t206 * t77) * t202;
t308 = t30 * t115;
t307 = t31 * t116;
t306 = rSges(5,3) + qJ(4);
t288 = t137 - t158;
t287 = t138 + t157;
t272 = t161 - t284;
t268 = -pkin(3) - t313;
t267 = -t282 / 0.2e1;
t266 = t282 / 0.2e1;
t262 = t189 * t267;
t261 = t189 * t266;
t260 = t190 * t267;
t259 = t190 * t266;
t258 = -t185 + t281;
t178 = rSges(2,1) * t207 - t205 * rSges(2,2);
t250 = rSges(2,1) * t205 + rSges(2,2) * t207;
t242 = -t189 * t35 + t190 * t36;
t241 = -t189 * t49 - t190 * t50;
t240 = -t189 * t81 - t190 * t80;
t239 = t189 * (Icges(6,5) * t130 + Icges(6,6) * t131) - t190 * (-Icges(6,5) * t132 - Icges(6,6) * t133);
t238 = t298 * t71 + t319;
t237 = -t139 + t258;
t236 = -t141 + t257;
t129 = t187 + t251;
t127 = -t166 * t200 - t276;
t25 = -t298 * t72 + t320;
t231 = (t189 * t238 + t190 * t25) * t202;
t27 = t296 * t72 + t243;
t230 = (-t189 * t26 + t190 * t27) * t202;
t225 = t276 + t280;
t128 = t248 - t325;
t219 = t225 - t333;
t94 = -t189 * t306 + t190 * t268 + t174 - t324;
t214 = t173 + t258 - t49 + t286;
t93 = -t325 + t306 * t190 + (-pkin(3) + t247) * t189;
t10 = qJD(5) * t230 + t309;
t44 = Icges(6,5) * t92 + Icges(6,6) * t91 - Icges(6,3) * t273;
t46 = Icges(6,4) * t92 + Icges(6,2) * t91 - Icges(6,6) * t273;
t48 = Icges(6,1) * t92 + Icges(6,4) * t91 - Icges(6,5) * t273;
t13 = -t203 * t44 + (-t204 * t46 + t206 * t48 + (-t204 * t77 - t206 * t74) * qJD(5)) * t202;
t43 = Icges(6,5) * t90 + Icges(6,6) * t89 - Icges(6,3) * t275;
t45 = Icges(6,4) * t90 + Icges(6,2) * t89 - Icges(6,6) * t275;
t47 = Icges(6,1) * t90 + Icges(6,4) * t89 - Icges(6,5) * t275;
t14 = -t203 * t43 + (-t204 * t45 + t206 * t47 + (t204 * t79 - t206 * t75) * qJD(5)) * t202;
t21 = -t132 * t147 + t133 * t148 + t137 * t89 + t138 * t90 + (-t136 * t299 + t146 * t190) * t202;
t22 = t130 * t147 - t131 * t148 + t137 * t91 + t138 * t92 + (-t136 * t297 - t146 * t189) * t202;
t51 = t130 * t137 - t131 * t138 - t136 * t298;
t42 = t51 * t182;
t9 = qJD(5) * t231 + t42;
t212 = (t42 + (t320 * t190 + (t238 + t243 - t27) * t189) * t282) * t260 + t308 / 0.2e1 + t307 / 0.2e1 + t51 * t329 + t52 * t328 + t321 + (-t309 + (-(-t26 - t320) * t189 + t238 * t190 + (-t243 - t319) * t190 - t25 * t189 + (-t72 * t189 ^ 2 + (-t189 * t71 - t190 * t72) * t190) * t202) * t282 + t10) * t261 + (t13 + t22) * t262 + (t14 + t21 + t9) * t259 + (Icges(5,2) * t203 ^ 2 + (Icges(5,1) * t202 + 0.2e1 * Icges(5,4) * t203) * t202 + Icges(4,3) + Icges(3,3)) * t199;
t54 = t189 * t227 + t290 + t300 - t325;
t106 = t200 * t217;
t65 = -t106 + t220;
t66 = t200 * t256 + t254;
t209 = (t66 * (-t279 + t325) - t65 * (t245 - t312) + (-t268 * t65 - t306 * t66) * t190) * t200;
t140 = t200 * t248;
t104 = t200 * t264 - t276;
t103 = -t140 + t226;
t102 = -rSges(6,1) * t132 - rSges(6,2) * t133;
t101 = rSges(6,1) * t130 + rSges(6,2) * t131;
t39 = t240 * t282 + qJD(3);
t15 = t115 * t81 - t116 * t80 + t241 * t282 + qJDD(3);
t6 = t130 * t45 - t131 * t47 + t75 * t91 - t79 * t92 + (-t189 * t43 - t297 * t72) * t202;
t5 = t130 * t46 - t131 * t48 + t74 * t91 + t77 * t92 + (-t189 * t44 - t297 * t71) * t202;
t4 = -t132 * t45 + t133 * t47 + t75 * t89 - t79 * t90 + (t190 * t43 - t299 * t72) * t202;
t3 = -t132 * t46 + t133 * t48 + t74 * t89 + t77 * t90 + (t190 * t44 - t299 * t71) * t202;
t1 = [Icges(2,3) * qJDD(1) + t212 + (t339 * (-t166 - t326) + t340 * (-t249 - t327) + (t127 - t144 + t276) * t126) * m(3) + ((qJDD(1) * t250 + g(3)) * t250 + (qJDD(1) * t178 + g(2)) * t178) * m(2) + (t36 * (t214 + t277) + t345 * (t54 - t327) + (t276 - t219 - t36 + t335) * t35 + t334 + t349 * (t346 - t326)) * m(6) + (t66 * (t272 + t277) + t209 + t343 * (t94 - t326) + t344 * (t93 - t327) + (-t219 - t66 - t254 + t337) * t65) * m(5) + (t104 * (t226 - t285) + t341 * (t129 - t326) + t342 * (t128 - t327) + (-t104 - t225 + t276 + t336) * t103) * m(4); t212 + (t345 * t54 + (t214 - t237 - t332) * t36 + (t236 + t335) * t35 + t334 + t349 * t346) * m(6) + (t209 + t343 * t94 + t344 * t93 + (t106 - t237 + t272) * t66 + (t236 - t186 + t337) * t65) * m(5) + (t341 * t129 + t342 * t128 + (-t280 + t336) * t103 + (t140 - t285) * t104) * m(4) + (-t126 * t144 + t127 * t143 + (-t126 * t200 - t339) * t166 - (t127 * t200 + t340) * t249) * m(3); m(6) * t15 + (m(4) + m(5)) * qJDD(3) + (-m(4) + t330) * g(1); t330 * (g(3) * t189 + t322) + m(5) * (t189 * t32 + t190 * t33) + m(6) * (t16 * t189 + t17 * t190); -t203 * (t308 + t307 + (-t13 * t189 + t14 * t190) * t282 + t321) / 0.2e1 + t181 * (-t203 * t67 + (-t189 * t30 + t190 * t31) * t202) / 0.2e1 + t182 * (-t203 * t41 + ((-t200 * t30 + t14) * t190 + (-t200 * t31 - t13) * t189) * t202) / 0.2e1 - (-t115 * t238 + t116 * t25 + t181 * t51 + t182 * t22 + (-t189 * t5 + t190 * t6) * t282) * t298 / 0.2e1 + (-t203 * t51 + t231) * t329 + (-t203 * t22 + ((t200 * t238 + t6) * t190 + (-t200 * t25 - t5) * t189) * t202) * t262 + (t115 * t26 + t116 * t27 + t181 * t52 + t182 * t21 + (-t189 * t3 + t190 * t4) * t282) * t296 / 0.2e1 + (-t203 * t52 + t230) * t328 + (-t203 * t21 + ((-t200 * t26 + t4) * t190 + (-t200 * t27 - t3) * t189) * t202) * t259 - t182 * (-t203 * t156 * t182 + ((-t204 * t287 - t206 * t288) * t182 + ((t204 * t223 + t206 * t221) * t202 + t239 * t203) * qJD(5)) * t202) / 0.2e1 + ((t130 * t287 + t131 * t288 - t156 * t298) * t182 + (-t130 * t223 - t131 * t221 + t239 * t298) * t282) * t261 + ((-t132 * t287 - t133 * t288 + t156 * t296) * t182 + (t132 * t223 + t133 * t221 - t239 * t296) * t282) * t260 - (t10 * t189 + t190 * t9) * t293 / 0.2e1 + ((-t16 * t80 + t17 * t81 + t35 * t50 + t36 * t49) * t203 + (t15 * t240 + t39 * (-t297 * t81 + t299 * t80 + t241) + t242 * t149 + ((-t200 * t35 + t17) * t190 + (-t200 * t36 + t16) * t189) * t145) * t202 - (-t101 * t35 - t102 * t36) * t182 - (t39 * (-t101 * t190 - t102 * t189) + t242 * t159) * t282 - g(1) * t159 - g(2) * t101 - g(3) * t102) * m(6);];
tau = t1;
