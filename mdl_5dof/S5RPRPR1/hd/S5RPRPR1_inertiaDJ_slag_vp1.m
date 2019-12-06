% Calculate time derivative of joint inertia matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:19
% DurationCPUTime: 6.19s
% Computational Cost: add. (5119->451), mult. (6088->621), div. (0->0), fcn. (4648->8), ass. (0->251)
t184 = sin(qJ(1));
t186 = cos(qJ(1));
t179 = qJ(3) + pkin(8);
t167 = qJ(5) + t179;
t161 = sin(t167);
t162 = cos(t167);
t232 = rSges(6,1) * t161 + rSges(6,2) * t162;
t80 = t186 * rSges(6,3) + t184 * t232;
t183 = sin(qJ(3));
t185 = cos(qJ(3));
t264 = qJD(3) * t185;
t269 = qJD(1) * t186;
t353 = t183 * t269 + t184 * t264;
t165 = sin(t179);
t166 = cos(t179);
t233 = rSges(5,1) * t165 + rSges(5,2) * t166;
t328 = pkin(3) * t183;
t132 = pkin(4) * t165 + t328;
t352 = t184 * t132 + t80;
t305 = Icges(4,4) * t185;
t219 = Icges(4,1) * t183 + t305;
t101 = Icges(4,5) * t186 + t184 * t219;
t306 = Icges(4,4) * t183;
t214 = Icges(4,2) * t185 + t306;
t99 = Icges(4,6) * t186 + t184 * t214;
t221 = t101 * t183 + t185 * t99;
t351 = t186 * t221;
t304 = Icges(5,4) * t165;
t212 = Icges(5,2) * t166 + t304;
t86 = Icges(5,6) * t186 + t184 * t212;
t303 = Icges(5,4) * t166;
t217 = Icges(5,1) * t165 + t303;
t88 = Icges(5,5) * t186 + t184 * t217;
t227 = t165 * t88 + t166 * t86;
t350 = t186 * t227;
t302 = Icges(6,4) * t161;
t211 = Icges(6,2) * t162 + t302;
t76 = Icges(6,6) * t186 + t184 * t211;
t301 = Icges(6,4) * t162;
t216 = Icges(6,1) * t161 + t301;
t78 = Icges(6,5) * t186 + t184 * t216;
t231 = t161 * t78 + t162 * t76;
t349 = t186 * t231;
t197 = t186 * t233;
t234 = rSges(4,1) * t183 + rSges(4,2) * t185;
t198 = t186 * t234;
t251 = t185 * t269;
t253 = t353 * rSges(4,1) + rSges(4,2) * t251;
t261 = -rSges(4,3) - pkin(1) - pkin(6);
t266 = qJD(3) * t183;
t274 = qJ(2) * t269 + qJD(2) * t184;
t34 = (-rSges(4,2) * t266 + qJD(1) * t261) * t184 + t253 + t274;
t320 = rSges(4,2) * t183;
t141 = rSges(4,1) * t185 - t320;
t170 = qJD(2) * t186;
t263 = qJD(3) * t186;
t35 = t170 + t141 * t263 + (t261 * t186 + (-qJ(2) - t234) * t184) * qJD(1);
t348 = t184 * t35 - t186 * t34;
t206 = Icges(6,5) * t161 + Icges(6,6) * t162;
t347 = -Icges(6,3) * t184 + t186 * t206;
t207 = Icges(5,5) * t165 + Icges(5,6) * t166;
t346 = -Icges(5,3) * t184 + t186 * t207;
t209 = Icges(4,5) * t183 + Icges(4,6) * t185;
t345 = -Icges(4,3) * t184 + t186 * t209;
t344 = -Icges(6,6) * t184 + t186 * t211;
t343 = -Icges(5,6) * t184 + t186 * t212;
t342 = -Icges(4,6) * t184 + t186 * t214;
t341 = -Icges(6,5) * t184 + t186 * t216;
t340 = -Icges(5,5) * t184 + t186 * t217;
t339 = -Icges(4,5) * t184 + t186 * t219;
t114 = Icges(6,5) * t162 - Icges(6,6) * t161;
t115 = -Icges(6,2) * t161 + t301;
t116 = Icges(6,1) * t162 - t302;
t178 = qJD(3) + qJD(5);
t91 = t211 * t178;
t92 = t216 * t178;
t338 = t161 * (t115 * t178 + t92) - t162 * (t116 * t178 - t91) + qJD(1) * t114;
t337 = 2 * m(4);
t336 = 2 * m(5);
t335 = 2 * m(6);
t180 = t184 ^ 2;
t181 = t186 ^ 2;
t334 = t184 / 0.2e1;
t333 = t186 / 0.2e1;
t332 = rSges(3,2) - pkin(1);
t331 = -rSges(5,3) - pkin(1);
t330 = -rSges(6,3) - pkin(1);
t329 = m(4) * t141;
t327 = pkin(3) * t185;
t326 = pkin(4) * t166;
t325 = pkin(6) * t186;
t182 = -qJ(4) - pkin(6);
t323 = rSges(5,1) * t166;
t321 = rSges(6,1) * t162;
t319 = rSges(5,2) * t165;
t316 = t165 * t86;
t315 = t165 * t343;
t314 = t166 * t88;
t313 = t166 * t340;
t312 = t183 * t99;
t311 = t184 * rSges(4,3);
t310 = t184 * rSges(5,3);
t175 = t186 * rSges(4,3);
t174 = t186 * rSges(5,3);
t280 = t183 * t184;
t158 = pkin(3) * t280;
t108 = t158 + (-pkin(6) - t182) * t186;
t93 = t233 * t184 + t174;
t307 = -t108 - t93;
t74 = Icges(6,3) * t186 + t184 * t206;
t291 = qJD(1) * t74;
t84 = Icges(5,3) * t186 + t184 * t207;
t290 = qJD(1) * t84;
t97 = Icges(4,3) * t186 + t184 * t209;
t289 = qJD(1) * t97;
t288 = t183 * t342;
t287 = t101 * t185;
t286 = t185 * t339;
t285 = t161 * t178;
t284 = t162 * t178;
t177 = -pkin(7) + t182;
t283 = t177 * t186;
t282 = t178 * t184;
t281 = t178 * t186;
t279 = t184 * t114;
t278 = t184 * t185;
t159 = pkin(3) * t278;
t257 = pkin(3) * t263;
t277 = qJD(1) * t159 + t183 * t257;
t276 = t182 * t269 + t185 * t257;
t275 = t184 * t182 + t186 * t328;
t273 = t186 * pkin(1) + t184 * qJ(2);
t272 = t180 + t181;
t270 = qJD(1) * t184;
t268 = qJD(3) * t165;
t267 = qJD(3) * t166;
t265 = qJD(3) * t184;
t262 = qJD(4) * t184;
t260 = t177 + t330;
t259 = -t108 + t158 - (-t177 + t182) * t186 - t352;
t258 = rSges(6,2) * t285;
t236 = t326 + t327;
t123 = t236 * qJD(3);
t256 = -t184 * t123 - t132 * t269 - t177 * t270;
t255 = t232 * t269 + t282 * t321;
t254 = -t233 * t269 - t265 * t323;
t252 = t353 * pkin(3) + t182 * t270;
t103 = rSges(4,1) * t280 + rSges(4,2) * t278 + t175;
t117 = -rSges(6,2) * t161 + t321;
t246 = t117 + t326;
t40 = qJD(1) * t76 - t115 * t281;
t245 = t178 * t341 - t40;
t41 = qJD(1) * t344 + t115 * t282;
t244 = t178 * t78 + t41;
t42 = qJD(1) * t78 - t116 * t281;
t243 = -t178 * t344 - t42;
t43 = qJD(1) * t341 + t116 * t282;
t242 = -t178 * t76 + t43;
t95 = t232 * t178;
t241 = t272 * t95;
t129 = t234 * qJD(3);
t238 = t272 * t129;
t168 = qJD(4) * t186;
t237 = t168 + t252;
t235 = t262 - t276;
t124 = -t319 + t323;
t230 = -t161 * t341 - t162 * t344;
t226 = -t165 * t340 - t166 * t343;
t172 = t186 * qJ(2);
t202 = t132 + t232;
t45 = t184 * t260 + t186 * t202 + t172;
t46 = t273 - t283 + t352;
t223 = t184 * t45 - t186 * t46;
t64 = t184 * t246 + t159;
t65 = (-t117 - t236) * t186;
t222 = t184 * t64 - t186 * t65;
t220 = Icges(4,1) * t185 - t306;
t218 = Icges(5,1) * t166 - t304;
t215 = -Icges(4,2) * t183 + t305;
t213 = -Icges(5,2) * t165 + t303;
t210 = Icges(4,5) * t185 - Icges(4,6) * t183;
t208 = Icges(5,5) * t166 - Icges(5,6) * t165;
t205 = -t183 * t339 - t185 * t342;
t204 = t115 * t162 + t116 * t161;
t14 = t184 * t231 + t186 * t74;
t196 = t230 * t184;
t15 = -t186 * t347 + t196;
t16 = t184 * t74 - t349;
t17 = -t184 * t347 - t186 * t230;
t38 = -t114 * t281 + t291;
t39 = qJD(1) * t347 + t178 * t279;
t203 = -t270 * (t14 * t186 + t15 * t184) + t184 * ((t184 * t38 + (-t16 + t196) * qJD(1)) * t184 + (t17 * qJD(1) + (-t161 * t43 - t162 * t41 - t284 * t78 + t285 * t76 + t291) * t186 + (t231 * qJD(1) + t243 * t161 + t245 * t162 + t39) * t184) * t186) + t186 * ((t186 * t39 + (t15 + t349) * qJD(1)) * t186 + (-t14 * qJD(1) + (t161 * t42 + t162 * t40 - t284 * t341 + t285 * t344) * t184 + (t38 + t244 * t162 + t242 * t161 + (t230 - t74) * qJD(1)) * t186) * t184) + (t16 * t186 + t17 * t184) * t269;
t200 = rSges(3,3) * t186 + t184 * t332;
t189 = qJD(1) * t204 - t206 * t178;
t199 = (t161 * t245 - t162 * t243 + t189 * t184 + t338 * t186) * t334 + (-t161 * t244 + t162 * t242 - t184 * t338 + t189 * t186) * t333 - (t114 * t186 - t161 * t76 + t162 * t78 + t184 * t204) * t270 / 0.2e1 + (t161 * t344 - t162 * t341 - t186 * t204 + t279) * t269 / 0.2e1;
t195 = t226 * t184;
t194 = t205 * t184;
t193 = qJD(3) * t220;
t192 = qJD(3) * t218;
t191 = qJD(3) * t215;
t190 = qJD(3) * t213;
t12 = t168 + (qJD(1) * t330 - t258) * t184 + t255 - t256 + t274;
t13 = -t262 + t170 + (t117 * t178 + t123) * t186 + (t260 * t186 + (-qJ(2) - t202) * t184) * qJD(1);
t188 = -t12 * t186 + t184 * t13 + (t184 * t46 + t186 * t45) * qJD(1);
t24 = t117 * t270 + t186 * t95 + (t165 * t263 + t166 * t270) * pkin(4) + t277;
t153 = pkin(3) * t251;
t25 = t153 + t246 * t269 + (-qJD(3) * t132 - t95) * t184;
t187 = t25 * t184 - t186 * t24 + (t184 * t65 + t186 * t64) * qJD(1);
t113 = t233 * qJD(3);
t107 = -t184 * pkin(6) - t275;
t106 = -rSges(3,2) * t186 + t184 * rSges(3,3) + t273;
t105 = t172 + t200;
t104 = t311 - t198;
t96 = t186 * t107;
t94 = t310 - t197;
t83 = (-t124 - t327) * t186;
t82 = t124 * t184 + t159;
t81 = t184 * rSges(6,3) - t186 * t232;
t73 = t170 + (t332 * t186 + (-rSges(3,3) - qJ(2)) * t184) * qJD(1);
t72 = qJD(1) * t200 + t274;
t71 = t186 * t81;
t70 = t103 + t273 + t325;
t69 = t184 * t261 + t172 + t198;
t68 = pkin(6) * t270 + t237;
t66 = -t132 * t186 - t184 * t177 + t275;
t63 = t186 * ((t158 - t325) * qJD(1) + t235);
t58 = qJD(1) * t345 + t210 * t265;
t57 = -t210 * t263 + t289;
t56 = -t182 * t186 + t158 + t273 + t93;
t55 = t184 * t331 + t172 + t197 + t275;
t50 = qJD(1) * t346 + t208 * t265;
t49 = -t208 * t263 + t290;
t48 = t124 * t269 + t153 + (-pkin(3) * t266 - t113) * t184;
t47 = t113 * t186 + t124 * t270 + t277;
t44 = (-rSges(6,3) * qJD(1) - t258) * t184 + t255;
t37 = -t184 * t80 + t71;
t36 = t186 * (qJD(1) * t80 - t117 * t281);
t31 = -t184 * t345 - t186 * t205;
t30 = t184 * t97 - t351;
t29 = -t186 * t345 + t194;
t28 = t221 * t184 + t186 * t97;
t23 = -t184 * t346 - t186 * t226;
t22 = t184 * t84 - t350;
t21 = -t186 * t346 + t195;
t20 = t227 * t184 + t186 * t84;
t19 = t170 + t124 * t263 + (t331 * t186 + (-qJ(2) - t233 - t328) * t184) * qJD(1) - t235;
t18 = (-rSges(5,2) * t268 + qJD(1) * t331) * t184 + t237 - t254 + t274;
t11 = t184 * t259 + t186 * t66 + t71 + t96;
t10 = -t184 * t44 + t36 + (-t184 * t81 - t186 * t80) * qJD(1);
t3 = t63 + t186 * (-t123 * t186 + t276) + t36 + (-t68 - t44 + t252 + t256) * t184 + ((t259 - t283) * t186 + (-t107 - t66 + t186 * (t132 - t328) - t81) * t184) * qJD(1);
t1 = [-t116 * t285 - t162 * t92 - t115 * t284 + t161 * t91 + (t18 * t56 + t19 * t55) * t336 + (t12 * t46 + t13 * t45) * t335 - t183 * t193 - t219 * t264 - t185 * t191 + t214 * t266 - t165 * t192 - t217 * t267 - t166 * t190 + t212 * t268 + (t34 * t70 + t35 * t69) * t337 + 0.2e1 * m(3) * (t105 * t73 + t106 * t72); m(5) * (-t18 * t186 + t184 * t19 + (t184 * t56 + t186 * t55) * qJD(1)) + m(6) * t188 + m(4) * ((t184 * t70 + t186 * t69) * qJD(1) + t348) + m(3) * (t184 * t73 - t186 * t72 + (t105 * t186 + t106 * t184) * qJD(1)); 0; ((t316 / 0.2e1 - t314 / 0.2e1 - t287 / 0.2e1 + t312 / 0.2e1 + t70 * t329) * t184 + (t315 / 0.2e1 - t313 / 0.2e1 + t69 * t329 + t288 / 0.2e1 - t286 / 0.2e1) * t186) * qJD(1) + t199 + m(4) * (t348 * t141 - (t184 * t69 - t186 * t70) * t129) + m(5) * (t18 * t83 + t19 * t82 + t47 * t56 + t48 * t55) + m(6) * (t12 * t65 + t13 * t64 + t24 * t46 + t25 * t45) + (-t165 * (qJD(1) * t86 - t213 * t263) + t166 * (qJD(1) * t88 - t218 * t263) - t183 * (qJD(1) * t99 - t215 * t263) + t185 * (qJD(1) * t101 - t220 * t263) + (-t205 - t226) * qJD(3)) * t334 + (-t165 * (qJD(1) * t343 + t184 * t190) + t166 * (qJD(1) * t340 + t184 * t192) - t183 * (qJD(1) * t342 + t184 * t191) + t185 * (qJD(1) * t339 + t184 * t193) + (-t221 - t227) * qJD(3)) * t333 + (-t209 - t207) * qJD(3) * (t181 / 0.2e1 + t180 / 0.2e1); m(5) * (t48 * t184 - t186 * t47 + (t184 * t83 + t186 * t82) * qJD(1)) + m(6) * t187 - m(4) * t238; (t11 * t3 + t24 * t65 + t25 * t64) * t335 + (t82 * t48 + t83 * t47 + (t184 * t307 + t186 * t94 + t96) * (t63 + (-t68 + t254) * t184 + (-t124 * t181 + t180 * t319) * qJD(3) + ((t307 + t174) * t186 + (-t107 + t197 - t94 + t310) * t184) * qJD(1))) * t336 + t184 * ((t184 * t49 + (-t22 + t195) * qJD(1)) * t184 + (t23 * qJD(1) + (-t267 * t88 + t268 * t86 + t290) * t186 + (t50 + (t313 - t315) * qJD(3) + t227 * qJD(1)) * t184) * t186) + t186 * ((t186 * t58 + (t29 + t351) * qJD(1)) * t186 + (-t28 * qJD(1) + (-t264 * t339 + t266 * t342) * t184 + (t57 + (t287 - t312) * qJD(3) + (t205 - t97) * qJD(1)) * t186) * t184) + t184 * ((t184 * t57 + (-t30 + t194) * qJD(1)) * t184 + (t31 * qJD(1) + (-t101 * t264 + t266 * t99 + t289) * t186 + (t58 + (t286 - t288) * qJD(3) + t221 * qJD(1)) * t184) * t186) + t186 * ((t186 * t50 + (t21 + t350) * qJD(1)) * t186 + (-t20 * qJD(1) + (-t267 * t340 + t268 * t343) * t184 + (t49 + (t314 - t316) * qJD(3) + (t226 - t84) * qJD(1)) * t186) * t184) + ((-t184 * t103 + t104 * t186) * (-t184 * t253 + (-t141 * t181 + t180 * t320) * qJD(3) + ((-t103 + t175) * t186 + (-t104 + t198 + t311) * t184) * qJD(1)) - t141 * t238) * t337 + t203 + ((-t20 - t28) * t186 + (-t21 - t29) * t184) * t270 + ((t22 + t30) * t186 + (t23 + t31) * t184) * t269; m(5) * (t184 * t18 + t186 * t19 + (-t184 * t55 + t186 * t56) * qJD(1)) + m(6) * (-qJD(1) * t223 + t184 * t12 + t13 * t186); 0; m(6) * (-qJD(1) * t222 + t184 * t24 + t186 * t25) + m(5) * (t184 * t47 + t186 * t48 + (-t184 * t82 + t186 * t83) * qJD(1)); 0; m(6) * (t117 * t188 - t223 * t95) + t199; -m(6) * t241; m(6) * (t10 * t11 + t117 * t187 - t222 * t95 + t37 * t3) + t203; 0; (t10 * t37 - t117 * t241) * t335 + t203;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
