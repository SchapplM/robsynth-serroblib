% Calculate vector of inverse dynamics joint torques for
% S5RRPPR1
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:43
% EndTime: 2020-01-03 11:55:50
% DurationCPUTime: 4.21s
% Computational Cost: add. (9636->417), mult. (5845->495), div. (0->0), fcn. (4304->10), ass. (0->244)
t199 = qJ(1) + qJ(2);
t189 = pkin(8) + t199;
t176 = sin(t189);
t177 = cos(t189);
t120 = t176 * rSges(4,1) + t177 * rSges(4,2);
t190 = sin(t199);
t179 = pkin(2) * t190;
t331 = t179 + t120;
t198 = qJD(1) + qJD(2);
t191 = cos(t199);
t332 = t190 * rSges(3,1) + t191 * rSges(3,2);
t111 = t332 * t198;
t203 = sin(qJ(1));
t307 = pkin(1) * qJD(1);
t259 = t203 * t307;
t102 = t259 + t111;
t197 = pkin(9) + qJ(5);
t187 = sin(t197);
t188 = cos(t197);
t172 = Icges(6,4) * t188;
t232 = -Icges(6,2) * t187 + t172;
t333 = Icges(6,1) * t187 + t172;
t278 = t333 + t232;
t298 = Icges(6,4) * t187;
t129 = Icges(6,2) * t188 + t298;
t132 = Icges(6,1) * t188 - t298;
t279 = t129 - t132;
t349 = (t187 * t278 + t188 * t279) * t198;
t135 = rSges(6,1) * t187 + rSges(6,2) * t188;
t100 = t135 * t177;
t171 = t176 * pkin(3);
t119 = -qJ(4) * t177 + t171;
t254 = t119 + t179;
t201 = cos(pkin(9));
t178 = pkin(4) * t201 + pkin(3);
t202 = -pkin(7) - qJ(4);
t275 = t176 * t178 + t177 * t202;
t290 = t176 * t188;
t144 = rSges(6,1) * t290;
t291 = t176 * t187;
t260 = rSges(6,2) * t291;
t85 = -rSges(6,3) * t177 + t144 - t260;
t338 = t275 + t85;
t314 = -t119 + t338;
t234 = t254 + t314;
t167 = qJD(4) * t176;
t268 = qJD(5) * t177;
t245 = -t135 * t268 + t167;
t27 = t198 * t234 - t245 + t259;
t204 = cos(qJ(1));
t186 = t204 * t307;
t282 = t191 * t198;
t163 = pkin(2) * t282;
t270 = qJD(4) * t177;
t252 = t163 - t270;
t269 = qJD(5) * t176;
t221 = -t135 * t269 + t252;
t316 = pkin(3) * t177;
t121 = qJ(4) * t176 + t316;
t143 = t177 * t178;
t73 = t316 - t143 + (qJ(4) + t202) * t176;
t287 = t177 * t187;
t145 = rSges(6,2) * t287;
t286 = t177 * t188;
t262 = rSges(6,1) * t286;
t248 = -t145 + t262;
t86 = rSges(6,3) * t176 + t248;
t264 = t121 - t73 + t86;
t28 = t198 * t264 + t186 + t221;
t99 = t135 * t176;
t348 = -t100 * t28 - t27 * t99;
t265 = qJD(5) * t198;
t107 = -qJDD(5) * t176 - t177 * t265;
t136 = rSges(6,1) * t188 - rSges(6,2) * t187;
t118 = t136 * qJD(5);
t285 = t177 * t198;
t150 = qJ(4) * t285;
t196 = qJDD(1) + qJDD(2);
t193 = t203 * pkin(1);
t205 = qJD(1) ^ 2;
t293 = pkin(1) * qJDD(1);
t249 = -t193 * t205 + t204 * t293;
t317 = pkin(2) * t196;
t195 = t198 ^ 2;
t318 = pkin(2) * t195;
t213 = -t190 * t318 + t191 * t317 + t249;
t210 = -qJDD(4) * t177 + t198 * t167 + t213;
t266 = qJD(5) * t188;
t267 = qJD(5) * t187;
t281 = rSges(6,3) * t285 + t198 * t260;
t289 = t176 * t198;
t55 = rSges(6,2) * t177 * t266 + (t177 * t267 + t188 * t289) * rSges(6,1) - t281;
t273 = t150 + t167;
t88 = pkin(3) * t289 - t273;
t8 = -t118 * t269 + t107 * t135 + t264 * t196 + (-t88 - t150 - t55 + (-t275 + t171) * t198) * t198 + t210;
t347 = -g(2) + t8;
t108 = -qJDD(5) * t177 + t176 * t265;
t123 = t178 * t285;
t274 = pkin(3) * t285 + qJ(4) * t289;
t244 = -t270 + t274;
t194 = t204 * pkin(1);
t271 = t205 * t194 + t203 * t293;
t251 = t190 * t317 + t191 * t318 + t271;
t212 = -qJDD(4) * t176 + t196 * t119 + t198 * t244 + t251;
t288 = t176 * t202;
t280 = rSges(6,3) * t289 + t198 * t262;
t56 = -rSges(6,1) * t176 * t267 + (-t176 * t266 - t187 * t285) * rSges(6,2) + t280;
t7 = t118 * t268 - t108 * t135 + t314 * t196 + (-t198 * t288 + t123 - t270 - t274 + t56) * t198 + t212;
t346 = -g(3) + t7;
t309 = rSges(5,1) * t201;
t157 = t176 * t309;
t200 = sin(pkin(9));
t308 = rSges(5,2) * t200;
t261 = t198 * t308;
t277 = rSges(5,3) * t285 + t176 * t261;
t158 = t177 * t308;
t263 = t177 * t309;
t90 = rSges(5,3) * t176 - t158 + t263;
t300 = t121 + t90;
t23 = (-t157 * t198 + t277 - t88) * t198 + t300 * t196 + t210;
t345 = -g(2) + t23;
t168 = t176 * rSges(4,2);
t122 = rSges(4,1) * t177 - t168;
t344 = -t120 * t195 + t196 * t122 - g(2) + t213;
t173 = t190 * rSges(3,2);
t138 = rSges(3,1) * t191 - t173;
t343 = -t111 * t198 + t138 * t196 - g(2) + t249;
t276 = rSges(5,3) * t289 + t198 * t263;
t247 = -t176 * t308 + t157;
t89 = -rSges(5,3) * t177 + t247;
t22 = t196 * t89 + ((-qJD(4) - t261) * t177 + t276) * t198 + t212;
t342 = -g(3) + t22;
t155 = rSges(4,1) * t285;
t341 = -g(3) + t196 * t120 + t198 * (-rSges(4,2) * t289 + t155) + t251;
t112 = rSges(3,1) * t282 - t173 * t198;
t340 = t112 * t198 + t196 * t332 - g(3) + t271;
t337 = t198 * (t254 + t89);
t336 = t198 * t331;
t110 = t198 * t122;
t335 = t155 - t110;
t334 = -t171 - t179;
t109 = t198 * t121;
t74 = t198 * t86;
t330 = t198 * t73 - t109 + t123 - t221 + t252 + t280 - t74;
t329 = -t198 * t90 - t109 + t163 + t244 - t252 + t276;
t230 = -t187 * t129 + t188 * t333;
t80 = -Icges(6,6) * t177 + t176 * t232;
t82 = -Icges(6,5) * t177 + t132 * t176;
t40 = t187 * t82 + t188 * t80;
t226 = t232 * t198;
t326 = -Icges(6,6) * t198 + qJD(5) * t129;
t52 = -t176 * t326 + t177 * t226;
t227 = t132 * t198;
t323 = -Icges(6,5) * t198 + qJD(5) * t333;
t54 = -t176 * t323 + t177 * t227;
t128 = Icges(6,5) * t188 - Icges(6,6) * t187;
t78 = -Icges(6,3) * t177 + t128 * t176;
t328 = qJD(5) * t40 + t187 * t52 - t188 * t54 - t198 * t78;
t127 = Icges(6,5) * t187 + Icges(6,6) * t188;
t327 = -Icges(6,3) * t198 + qJD(5) * t127;
t116 = t232 * qJD(5);
t117 = t132 * qJD(5);
t231 = t129 * t188 + t187 * t333;
t325 = qJD(5) * t231 + t116 * t187 - t117 * t188 - t127 * t198;
t81 = Icges(6,4) * t286 - Icges(6,2) * t287 + Icges(6,6) * t176;
t141 = Icges(6,4) * t287;
t83 = Icges(6,1) * t286 + Icges(6,5) * t176 - t141;
t236 = t187 * t83 + t188 * t81;
t51 = t176 * t226 + t177 * t326;
t53 = t176 * t227 + t177 * t323;
t79 = Icges(6,5) * t286 - Icges(6,6) * t287 + Icges(6,3) * t176;
t324 = qJD(5) * t236 - t187 * t51 + t188 * t53 - t198 * t79;
t321 = -m(5) - m(6);
t320 = t107 / 0.2e1;
t319 = t108 / 0.2e1;
t313 = t176 * t333 + t80;
t312 = t177 * t333 + t81;
t311 = -t129 * t176 + t82;
t310 = -Icges(6,2) * t286 - t141 + t83;
t306 = t187 * t80;
t305 = t187 * t81;
t304 = t188 * t82;
t303 = t198 * t28;
t292 = t127 * t176;
t48 = t177 * t230 + t292;
t302 = t48 * t198;
t301 = rSges(5,3) + qJ(4);
t94 = t127 * t177;
t225 = t128 * t198;
t258 = -t269 / 0.2e1;
t257 = t269 / 0.2e1;
t256 = -t268 / 0.2e1;
t255 = t268 / 0.2e1;
t103 = t138 * t198 + t186;
t246 = -t167 + t259;
t180 = pkin(2) * t191;
t105 = t122 + t180;
t166 = rSges(2,1) * t204 - t203 * rSges(2,2);
t165 = rSges(2,1) * t203 + rSges(2,2) * t204;
t241 = -t176 * t28 + t177 * t27;
t67 = t82 * t290;
t29 = -t177 * t78 - t291 * t80 + t67;
t68 = t83 * t290;
t30 = t177 * t79 + t291 * t81 - t68;
t240 = -t30 * t176 - t29 * t177;
t69 = t80 * t287;
t31 = -t176 * t78 - t286 * t82 + t69;
t235 = -t188 * t83 + t305;
t32 = t176 * t79 - t235 * t177;
t239 = -t32 * t176 - t31 * t177;
t238 = t176 * t85 + t177 * t86;
t237 = t304 - t306;
t220 = t187 * t311 + t188 * t313;
t219 = t187 * t310 + t188 * t312;
t216 = -t176 * t225 - t177 * t327 + t198 * t235;
t215 = t176 * t327 - t177 * t225 + t198 * t237;
t214 = -t128 * qJD(5) + t198 * t230;
t59 = t179 + t338;
t60 = t143 + t180 + (rSges(6,3) - t202) * t176 + t248;
t65 = -t177 * t301 + t247 - t334;
t66 = -t158 + t180 + (pkin(3) + t309) * t177 + t301 * t176;
t75 = t259 + t336;
t76 = t186 + t163 + t110;
t209 = (-t75 * t168 - t76 * t331) * t198;
t47 = t176 * t230 - t94;
t46 = t47 * t198;
t11 = qJD(5) * t240 + t46;
t12 = qJD(5) * t239 - t302;
t16 = qJD(5) * t237 + t187 * t54 + t188 * t52;
t17 = qJD(5) * t235 + t187 * t53 + t188 * t51;
t20 = t214 * t176 + t177 * t325;
t21 = -t176 * t325 + t214 * t177;
t208 = (t46 + ((t31 + t68 - t69 + (t78 - t305) * t176) * t176 + (-t67 - t32 + (-t235 + t78) * t177 + (t304 + t306) * t176) * t177) * qJD(5)) * t257 - t236 * t320 - t107 * t48 / 0.2e1 + (qJD(5) * t230 + t116 * t188 + t117 * t187) * t198 + (t40 + t47) * t319 + (t12 + t302 + ((-t30 + t69 + (t79 - t304) * t177) * t177 + (t29 - t67 + (t79 + t306) * t176) * t176) * qJD(5)) * t255 + (t17 + t20 + t11) * t258 + (t16 + t21) * t256 + (Icges(4,3) + Icges(3,3) + t231 + Icges(5,2) * t201 ^ 2 + (Icges(5,1) * t200 + 0.2e1 * Icges(5,4) * t201) * t200) * t196;
t44 = t246 + t337;
t45 = t198 * t300 + t186 + t252;
t207 = (t45 * (-t157 + t334) - t44 * t158) * t198;
t206 = (t28 * (-t144 - t275 - t179) + t27 * (-t145 - t288)) * t198 + t348 * qJD(5);
t39 = qJD(5) * t238 + qJD(3);
t13 = -t107 * t85 - t108 * t86 + qJDD(3) + (t176 * t56 - t177 * t55) * qJD(5);
t6 = t176 * t324 + t216 * t177;
t5 = -t176 * t328 + t215 * t177;
t4 = t216 * t176 - t177 * t324;
t3 = t215 * t176 + t177 * t328;
t1 = [Icges(2,3) * qJDD(1) + t208 + (t343 * (t138 + t194) + t340 * (t193 + t332) + (-t103 + t112 + t186) * t102) * m(3) + ((t165 ^ 2 + t166 ^ 2) * qJDD(1) - g(2) * t166 - g(3) * t165) * m(2) + (t28 * (-t246 + t281) + t206 + t347 * (t194 + t60) + t346 * (t193 + t59) + (t28 + t330) * t27) * m(6) + (t45 * (t150 - t246 + t277) + t207 + t345 * (t194 + t66) + t342 * (t193 + t65) + (t45 + t329) * t44) * m(5) + (-t76 * t259 + t209 + t344 * (t105 + t194) + t341 * (t193 + t331) + (t76 + t335) * t75) * m(4); t208 + (t234 * t303 + t206 + t347 * t60 + t346 * t59 + (-t245 + t167 + t281) * t28 + t330 * t27) * m(6) + (t207 + t345 * t66 + t342 * t65 + t329 * t44 + (-t167 + t273 + t277 + t337) * t45) * m(5) + (t344 * t105 + t331 * t341 + t335 * t75 + t336 * t76 + t209) * m(4) + (t102 * t112 - t103 * t111 + (-t102 * t198 + t343) * t138 + (t103 * t198 + t340) * t332) * m(3); m(6) * t13 + (m(4) + m(5)) * qJDD(3) + (-m(4) + t321) * g(1); t321 * (-g(2) * t177 - g(3) * t176) + m(5) * (-t176 * t22 - t177 * t23) + m(6) * (-t176 * t7 - t177 * t8); t196 * (t176 * t236 - t177 * t40) / 0.2e1 + t198 * ((t198 * t236 - t16) * t177 + (t198 * t40 - t17) * t176) / 0.2e1 + t11 * t289 / 0.2e1 - t177 * (t107 * t30 + t29 * t108 + t47 * t196 + t21 * t198 + (-t176 * t6 - t177 * t5) * qJD(5)) / 0.2e1 + t240 * t319 + ((-t198 * t30 - t5) * t177 + (t198 * t29 - t6) * t176) * t256 - t12 * t285 / 0.2e1 - t176 * (t32 * t107 + t31 * t108 - t48 * t196 + t20 * t198 + (-t176 * t4 - t177 * t3) * qJD(5)) / 0.2e1 + t239 * t320 + ((-t198 * t32 - t3) * t177 + (t198 * t31 - t4) * t176) * t258 - t198 * ((-t187 * t279 + t188 * t278) * t198 + ((t176 * t310 - t177 * t311) * t188 + (-t176 * t312 + t177 * t313) * t187) * qJD(5)) / 0.2e1 + ((-t268 * t292 - t225) * t177 + (-t349 + (-t219 * t176 + (t94 + t220) * t177) * qJD(5)) * t176) * t255 + ((t94 * t269 - t225) * t176 + (t349 + (-t220 * t177 + (-t292 + t219) * t176) * qJD(5)) * t177) * t257 + (t13 * t238 + t39 * ((t198 * t85 - t55) * t177 + (t56 - t74) * t176) + t241 * t118 + ((t7 - t303) * t177 + (-t198 * t27 - t8) * t176) * t135 - t348 * t198 - (t39 * (-t100 * t177 - t176 * t99) + t241 * t136) * qJD(5) - g(1) * t136 + g(2) * t99 - g(3) * t100) * m(6);];
tau = t1;
